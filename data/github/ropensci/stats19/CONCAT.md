---
title: 'stats19: A package for working with open road crash data'
tags:
  - stats19
  - crashes
  - dft
  - road safety
authors:
 - name: R Lovelace
   orcid: 0000-0001-5679-6536
   affiliation: 1
 - name: M Morgan
   affiliation: 1
 - name: L Hama
   orcid: 0000-0003-1912-4890
   affiliations: 1
 - name: M Padgham
   orcid: 0000-0003-2172-5265
   affiliations: 2
affiliations:
 - name: Institute for Transport Studies (ITS) and Leeds Institute for Data Analytics (LIDA), University of Leeds
   index: 1
 - name: ATFutures GmbH.
   index: 2
date: 19 December 2018
bibliography: paper.bib
---

# Summary

**stats19** provides functions for downloading and formatting road crash data.
Specifically, it enables access to the UK's official road traffic casualty database, [STATS19](https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data) (the name comes from the form used by the police to record car crashes and other incidents resulting in casualties on the roads).
Finding, reading-in and formatting the data for research can be a time consuming process subject to human error, leading to previous (incomplete) attempts to facilitate the processes with open source software [@lovelace_stplanr_Inpress].
**stats19** speeds-up these data access and cleaning stages by streamlining the work into 3 stages:

1. **Download** the data, by `year`, `type` and/or `filename`. An interactive menu of options is provided if there are multiple matches for a particular year.
2. **Read** the data in and with appropriate formatting of columns.
3. **Format** the data so that labels are added to the raw integer values for each column.

Functions for each stage are named `dl_stats19()`, `read_*()` and `format_*()`, with `*` representing the type of data to be read-in: STATS19 data consists of `accidents`, `casualties` and `vehicles` tables, which correspond to incident records, people injured or killed, and vehicles involved, respectively.

The package is needed because currently downloading and formatting STATS19 data is a time-consuming and error-prone process.
By abstracting the process to its fundamental steps (download, read, format), the package makes it easy to get the data into appropriate formats (of classes `tbl`, `data.frame` and `sf`), ready for for further processing and analysis steps.
We developed the package for road safety research, building on a clear need for reproducibility in the field [@lovelace_who_2016] and the importance of the geo-location in STATS19 data for assessing the effectiveness of interventions aimed to make roads safer and save lives [@sarkar_street_2018].
A useful feature of the package is that it enables creation of geographic representations of the data, geo-referenced to the correct coordinate reference system, in a single function call (`format_sf()`).
The package will be of use and interest to road safety data analysts working at local authority and national levels in the UK.
The datasets generated will also be of interest to academics and educators as an open, reproducible basis for analysing large point pattern data on an underlying route network, and for teaching on geography, transport and road safety courses.

# References 


<!-- badges: start -->
<!-- [![Travis build status](https://travis-ci.org/ropensci/stats19.svg?branch=master)](https://travis-ci.org/ropensci/stats19) -->

[![](http://www.r-pkg.org/badges/version/stats19)](https://www.r-pkg.org/pkg/stats19)
[![R-CMD-check](https://github.com/ropensci/stats19/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/stats19/actions)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/stats19)](https://www.r-pkg.org/pkg/stats19)
[![Life
cycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![](https://badges.ropensci.org/266_status.svg)](https://github.com/ropensci/software-review/issues/266)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.01181/status.svg)](https://doi.org/10.21105/joss.01181)
![codecov](https://codecov.io/gh/ropensci/stats19/branch/master/graph/badge.svg)
<!-- badges: end -->

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2540781.svg)](https://doi.org/10.5281/zenodo.2540781) -->
<!-- [![Gitter chat](https://badges.gitter.im/ITSLeeds/stats19.png)](https://gitter.im/stats19/Lobby?source=orgpage) -->
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stats19 <a href='https://docs.ropensci.org/stats19/'><img src='https://raw.githubusercontent.com/ropensci/stats19/master/man/figures/logo.png' align="right" height=215/></a>

**stats19** provides functions for downloading and formatting road crash
data. Specifically, it enables access to the UK’s official road traffic
casualty database,
[STATS19](https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data).
(The name comes from the form used by the police to record car crashes
and other incidents resulting in casualties on the roads.)

A full overview of STATS19 variables be found in a
[document](https://data.dft.gov.uk/road-accidents-safety-data/Brief-guide-to%20road-accidents-and-safety-data.doc)
provided by the UK’s Department for Transport (DfT).

The raw data is provided as a series of `.csv` files that contain
integers and which are stored in dozens of `.zip` files. Finding,
reading-in and formatting the data for research can be a time consuming
process subject to human error. **stats19** speeds up these vital but
boring and error-prone stages of the research process with a single
function: `get_stats19()`. By allowing public access to properly
labelled road crash data, **stats19** aims to make road safety research
more reproducible and accessible.

For transparency and modularity, each stage can be undertaken
separately, as documented in the [stats19
vignette](https://itsleeds.github.io/stats19/articles/stats19.html).

## Installation

Install and load the latest version with:

``` r
remotes::install_github("ropensci/stats19")
```

``` r
library(stats19)
#> Data provided under OGL v3.0. Cite the source and link to:
#> www.nationalarchives.gov.uk/doc/open-government-licence/version/3/
```

You can install the released version of stats19 from
[CRAN](https://cran.r-project.org/package=stats19) with:

``` r
install.packages("stats19")
```

## get_stats19()

`get_stats19()` requires `year` and `type` parameters, mirroring the
provision of STATS19 data files, which are categorised by year (from
1979 onward) and type (with separate tables for crashes, casualties and
vehicles, as outlined below). The following command, for example, gets
crash data from 2017 (**note**: we follow the “crash not accident”
campaign of
[RoadPeace](https://www.roadpeace.org/take-action/crash-not-accident/)
in naming crashes, although the DfT refers to the relevant tables as
‘accidents’ data):

``` r
crashes = get_stats19(year = 2017, type = "accident")
#> Files identified: dft-road-casualty-statistics-accident-2017.csv
#>    https://data.dft.gov.uk/road-accidents-safety-data/dft-road-casualty-statistics-accident-2017.csv
#> Data already exists in data_dir, not downloading
#> Data saved at ~/stats19-data/dft-road-casualty-statistics-accident-2017.csv
#> Reading in:
#> ~/stats19-data/dft-road-casualty-statistics-accident-2017.csv
#> Rows: 129982 Columns: 36
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr   (8): accident_index, accident_reference, longitude, latitude, date, lo...
#> dbl  (27): accident_year, location_easting_osgr, location_northing_osgr, pol...
#> time  (1): time
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> date and time columns present, creating formatted datetime column
```

What just happened? For the `year` 2017 we read-in crash-level
(`type = "accident"`) data on all road crashes recorded by the police
across Great Britain. The dataset contains 37 columns (variables) for
129,982 crashes. We were not asked to download the file (by default you
are asked to confirm the file that will be downloaded). The contents of
this dataset, and other datasets provided by **stats19**, are outlined
below and described in more detail in the [stats19
vignette](https://itsleeds.github.io/stats19/articles/stats19.html).

We will see below how the function also works to get the corresponding
casualty and vehicle datasets for 2017. The package also allows STATS19
files to be downloaded and read-in separately, allowing more control
over what you download, and subsequently read-in, with
`read_accidents()`, `read_casualties()` and `read_vehicles()`, as
described in the vignette.

## Data download

Data files can be downloaded without reading them in using the function
`dl_stats19()`. If there are multiple matches, you will be asked to
choose from a range of options. Providing just the year, for example,
will result in the following options:

``` r
dl_stats19(year = 2020, data_dir = tempdir())
#> Files identified: dft-road-casualty-statistics-accident-2020.csv
#>    https://data.dft.gov.uk/road-accidents-safety-data/dft-road-casualty-statistics-accident-2020.csv
#> Attempt downloading from: https://data.dft.gov.uk/road-accidents-safety-data/dft-road-casualty-statistics-accident-2020.csv
#> checking...
#> https://data.dft.gov.uk/road-accidents-safety-data/dft-road-casualty-statistics-accident-2020.csv not available currently
#> NULL
```

    Multiple matches. Which do you want to download?

    1: dft-road-casualty-statistics-casualty-2020.csv
    2: dft-road-casualty-statistics-vehicle-2020.csv
    3: dft-road-casualty-statistics-accident-2020.csv
    4: dft-road-casualty-statistics-vehicle-e-scooter-2020.csv
    5: road-safety-data-mid-year-provisional-estimates-2020.csv

    Selection: 
    Enter an item from the menu, or 0 to exit

## Using the data

STATS19 data consists of 3 main tables:

-   Accidents, the main table which contains information on the crash
    time, location and other variables (32 columns in total)
-   Casualties, containing data on people hurt or killed in each crash
    (16 columns in total)
-   Vehicles, containing data on vehicles involved in or causing each
    crash (23 columns in total)

The contents of each is outlined below.

### Crash data

Crash data was downloaded and read-in using the function
`get_stats19()`, as described above.

``` r
nrow(crashes)
#> [1] 129982
ncol(crashes)
#> [1] 37
```

Some of the key variables in this dataset include:

``` r
key_column_names = grepl(pattern = "severity|speed|pedestrian|light_conditions", x = names(crashes))
crashes[key_column_names]
#> # A tibble: 129,982 × 5
#>    accident_severity speed_limit pedestrian_crossing… pedestrian_crossing_physi…
#>    <chr>                   <dbl> <chr>                <chr>                     
#>  1 Fatal                      30 None within 50 metr… No physical crossing faci…
#>  2 Slight                     30 None within 50 metr… No physical crossing faci…
#>  3 Slight                     30 None within 50 metr… No physical crossing faci…
#>  4 Slight                     30 None within 50 metr… Pelican, puffin, toucan o…
#>  5 Serious                    20 None within 50 metr… Pedestrian phase at traff…
#>  6 Slight                     30 None within 50 metr… No physical crossing faci…
#>  7 Slight                     40 None within 50 metr… Pelican, puffin, toucan o…
#>  8 Slight                     30 Control by other au… Pedestrian phase at traff…
#>  9 Serious                    50 None within 50 metr… Pelican, puffin, toucan o…
#> 10 Serious                    30 None within 50 metr… No physical crossing faci…
#> # … with 129,972 more rows, and 1 more variable: light_conditions <chr>
```

For the full list of columns, run `names(crashes)` or see the
[vignette](https://github.com/ropensci/stats19/blob/master/vignettes/stats19.Rmd).

<!-- This means `crashes` is much more usable than `crashes_raw`, as shown below, which shows three records and some key variables in the messy and clean datasets: -->

### Casualties data

As with `crashes`, casualty data for 2017 can be downloaded, read-in and
formatted as follows:

``` r
casualties = get_stats19(year = 2017, type = "casualty", ask = FALSE)
#> Files identified: dft-road-casualty-statistics-casualty-2017.csv
#>    https://data.dft.gov.uk/road-accidents-safety-data/dft-road-casualty-statistics-casualty-2017.csv
#> Data already exists in data_dir, not downloading
#> Data saved at ~/stats19-data/dft-road-casualty-statistics-casualty-2017.csv
#> Rows: 170993 Columns: 18
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr  (2): accident_index, accident_reference
#> dbl (16): accident_year, vehicle_reference, casualty_reference, casualty_cla...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
nrow(casualties)
#> [1] 170993
ncol(casualties)
#> [1] 18
```

The results show that there were 170,993 casualties reported by the
police in the STATS19 dataset in 2017, and 18 columns (variables).
Values for a sample of these columns are shown below:

``` r
casualties[c(4, 5, 6, 14)]
#> # A tibble: 170,993 × 4
#>    vehicle_reference casualty_reference casualty_class  bus_or_coach_passenger  
#>                <dbl>              <dbl> <chr>           <chr>                   
#>  1                 1                  1 Passenger       Not a bus or coach pass…
#>  2                 2                  2 Driver or rider Not a bus or coach pass…
#>  3                 2                  3 Passenger       Not a bus or coach pass…
#>  4                 1                  1 Passenger       Not a bus or coach pass…
#>  5                 3                  1 Driver or rider Not a bus or coach pass…
#>  6                 1                  1 Passenger       Not a bus or coach pass…
#>  7                 1                  1 Pedestrian      Not a bus or coach pass…
#>  8                 2                  1 Driver or rider Not a bus or coach pass…
#>  9                 1                  1 Driver or rider Not a bus or coach pass…
#> 10                 2                  2 Driver or rider Not a bus or coach pass…
#> # … with 170,983 more rows
```

The full list of column names in the `casualties` dataset is:

``` r
names(casualties)
#>  [1] "accident_index"                     "accident_year"                     
#>  [3] "accident_reference"                 "vehicle_reference"                 
#>  [5] "casualty_reference"                 "casualty_class"                    
#>  [7] "sex_of_casualty"                    "age_of_casualty"                   
#>  [9] "age_band_of_casualty"               "casualty_severity"                 
#> [11] "pedestrian_location"                "pedestrian_movement"               
#> [13] "car_passenger"                      "bus_or_coach_passenger"            
#> [15] "pedestrian_road_maintenance_worker" "casualty_type"                     
#> [17] "casualty_home_area_type"            "casualty_imd_decile"
```

### Vehicles data

Data for vehicles involved in crashes in 2017 can be downloaded, read-in
and formatted as follows:

``` r
vehicles = get_stats19(year = 2017, type = "vehicle", ask = FALSE)
#> Files identified: dft-road-casualty-statistics-vehicle-2017.csv
#>    https://data.dft.gov.uk/road-accidents-safety-data/dft-road-casualty-statistics-vehicle-2017.csv
#> Data already exists in data_dir, not downloading
#> Data saved at ~/stats19-data/dft-road-casualty-statistics-vehicle-2017.csv
#> Rows: 238926 Columns: 27
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr  (2): accident_index, accident_reference
#> dbl (25): accident_year, vehicle_reference, vehicle_type, towing_and_articul...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
nrow(vehicles)
#> [1] 238926
ncol(vehicles)
#> [1] 27
```

The results show that there were 238,926 vehicles involved in crashes
reported by the police in the STATS19 dataset in 2017, with 27 columns
(variables). Values for a sample of these columns are shown below:

``` r
vehicles[c(3, 14:16)]
#> # A tibble: 238,926 × 4
#>    accident_reference vehicle_leaving_car… hit_object_off_car… first_point_of_i…
#>    <chr>              <chr>                <chr>               <chr>            
#>  1 010001708          Did not leave carri… None                Front            
#>  2 010001708          Did not leave carri… None                Back             
#>  3 010009342          Did not leave carri… None                Back             
#>  4 010009342          Did not leave carri… None                Front            
#>  5 010009344          Did not leave carri… None                Front            
#>  6 010009344          Did not leave carri… None                Front            
#>  7 010009344          Did not leave carri… None                Front            
#>  8 010009348          Did not leave carri… None                Front            
#>  9 010009348          Did not leave carri… None                Offside          
#> 10 010009350          Did not leave carri… None                Offside          
#> # … with 238,916 more rows
```

The full list of column names in the `vehicles` dataset is:

``` r
names(vehicles)
#>  [1] "accident_index"                   "accident_year"                   
#>  [3] "accident_reference"               "vehicle_reference"               
#>  [5] "vehicle_type"                     "towing_and_articulation"         
#>  [7] "vehicle_manoeuvre"                "vehicle_direction_from"          
#>  [9] "vehicle_direction_to"             "vehicle_location_restricted_lane"
#> [11] "junction_location"                "skidding_and_overturning"        
#> [13] "hit_object_in_carriageway"        "vehicle_leaving_carriageway"     
#> [15] "hit_object_off_carriageway"       "first_point_of_impact"           
#> [17] "vehicle_left_hand_drive"          "journey_purpose_of_driver"       
#> [19] "sex_of_driver"                    "age_of_driver"                   
#> [21] "age_band_of_driver"               "engine_capacity_cc"              
#> [23] "propulsion_code"                  "age_of_vehicle"                  
#> [25] "generic_make_model"               "driver_imd_decile"               
#> [27] "driver_home_area_type"
```

## Creating geographic crash data

An important feature of STATS19 data is that the “accidents” table
contains geographic coordinates. These are provided at \~10m resolution
in the UK’s official coordinate reference system (the Ordnance Survey
National Grid, EPSG code 27700). **stats19** converts the non-geographic
tables created by `format_accidents()` into the geographic data form of
the [`sf` package](https://cran.r-project.org/package=sf) with the
function `format_sf()` as follows:

``` r
crashes_sf = format_sf(crashes)
#> Warning: One or more parsing issues, see `problems()` for details
#> 19 rows removed with no coordinates
```

The note arises because `NA` values are not permitted in `sf`
coordinates, and so rows containing no coordinates are automatically
removed. Having the data in a standard geographic form allows various
geographic operations to be performed on it. The following code chunk,
for example, returns all crashes within the boundary of West Yorkshire
(which is contained in the object
[`police_boundaries`](https://itsleeds.github.io/stats19/reference/police_boundaries.html),
an `sf` data frame containing all police jurisdictions in England and
Wales).

``` r
library(sf)
#> Linking to GEOS 3.9.1, GDAL 3.3.2, PROJ 7.2.1
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
wy = filter(police_boundaries, pfa16nm == "West Yorkshire")
#> old-style crs object detected; please recreate object with a recent sf::st_crs()
crashes_wy = crashes_sf[wy, ]
nrow(crashes_sf)
#> [1] 129963
nrow(crashes_wy)
#> [1] 4371
```

This subsetting has selected the 4,371 crashes which occurred within
West Yorkshire in 2017.

## Joining tables

The three main tables we have just read-in can be joined by shared key
variables. This is demonstrated in the code chunk below, which subsets
all casualties that took place in Leeds, and counts the number of
casualties by severity for each crash:

``` r
sel = casualties$accident_index %in% crashes_wy$accident_index
casualties_wy = casualties[sel, ]
names(casualties_wy)
#>  [1] "accident_index"                     "accident_year"                     
#>  [3] "accident_reference"                 "vehicle_reference"                 
#>  [5] "casualty_reference"                 "casualty_class"                    
#>  [7] "sex_of_casualty"                    "age_of_casualty"                   
#>  [9] "age_band_of_casualty"               "casualty_severity"                 
#> [11] "pedestrian_location"                "pedestrian_movement"               
#> [13] "car_passenger"                      "bus_or_coach_passenger"            
#> [15] "pedestrian_road_maintenance_worker" "casualty_type"                     
#> [17] "casualty_home_area_type"            "casualty_imd_decile"
cas_types = casualties_wy %>%
  select(accident_index, casualty_type) %>%
  mutate(n = 1) %>%
  group_by(accident_index, casualty_type) %>%
  summarise(n = sum(n)) %>%
  tidyr::spread(casualty_type, n, fill = 0)
cas_types$Total = rowSums(cas_types[-1])
cj = left_join(crashes_wy, cas_types, by = "accident_index")
```

What just happened? We found the subset of casualties that took place in
West Yorkshire with reference to the `accident_index` variable. Then we
used functions from the **tidyverse** package **dplyr** (and `spread()`
from **tidyr**) to create a dataset with a column for each casualty
type. We then joined the updated casualty data onto the `crashes_wy`
dataset. The result is a spatial (`sf`) data frame of crashes in Leeds,
with columns counting how many road users of different types were hurt.
The original and joined data look like this:

``` r
crashes_wy %>%
  select(accident_index, accident_severity) %>% 
  st_drop_geometry()
#> # A tibble: 4,371 × 2
#>    accident_index accident_severity
#>  * <chr>          <chr>            
#>  1 2017120009776  Slight           
#>  2 2017120010412  Slight           
#>  3 2017120111341  Serious          
#>  4 2017120135780  Slight           
#>  5 2017120223550  Slight           
#>  6 2017130086428  Slight           
#>  7 20171333D0295  Slight           
#>  8 2017133AP0313  Serious          
#>  9 2017133BE0850  Slight           
#> 10 2017134110858  Slight           
#> # … with 4,361 more rows
cas_types[1:2, c("accident_index", "Cyclist")]
#> # A tibble: 2 × 2
#> # Groups:   accident_index [2]
#>   accident_index Cyclist
#>   <chr>            <dbl>
#> 1 2017120009776        0
#> 2 2017120010412        1
cj[1:2, c(1, 5, 34)] %>% st_drop_geometry()
#> # A tibble: 2 × 3
#>   accident_index latitude  lsoa_of_accident_location
#> * <chr>          <chr>     <chr>                    
#> 1 2017120009776  53.644355 E01027923                
#> 2 2017120010412  53.929575 E01027735
```

## Mapping crashes

The join operation added a geometry column to the casualty data,
enabling it to be mapped (for more advanced maps, see the
[vignette](https://itsleeds.github.io/stats19/articles/stats19.html)):

``` r
cex = cj$Total / 3
plot(cj["speed_limit"], cex = cex)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

The spatial distribution of crashes in West Yorkshire clearly relates to
the region’s geography. Crashes tend to happen on busy Motorway roads
(with a high speed limit, of 70 miles per hour, as shown in the map
above) and city centres, of Leeds and Bradford in particular. The
severity and number of people hurt (proportional to circle width in the
map above) in crashes is related to the speed limit.

STATS19 data can be used as the basis of road safety research. The map
below, for example, shows the results of an academic paper on the
social, spatial and temporal distribution of bike crashes in West
Yorkshire, which estimated the number of crashes per billion km cycled
based on commuter cycling as a proxy for cycling levels overall (more
sophisticated measures of cycling levels are now possible thanks to new
data sources) (Lovelace, Roberts, and Kellar 2016):

<img src="https://ars.els-cdn.com/content/image/1-s2.0-S136984781500039X-gr9.jpg" width="100%" />

## Time series analysis

We can also explore seasonal trends in crashes by aggregating crashes by
day of the year:

``` r
library(ggplot2)
crashes_dates = cj %>% 
  st_set_geometry(NULL) %>% 
  group_by(date) %>% 
  summarise(
    walking = sum(Pedestrian),
    cycling = sum(Cyclist),
    passenger = sum(`Car occupant`)
    ) %>% 
  tidyr::gather(mode, casualties, -date)
ggplot(crashes_dates, aes(date, casualties)) +
  geom_smooth(aes(colour = mode), method = "loess") +
  ylab("Casualties per day")
#> `geom_smooth()` using formula 'y ~ x'
```

<img src="man/figures/README-crash-date-plot-1.png" width="100%" />

Different types of crashes also tend to happen at different times of
day. This is illustrated in the plot below, which shows the times of day
when people who were travelling by different modes were most commonly
injured.

``` r
library(stringr)

crash_times = cj %>% 
  st_set_geometry(NULL) %>% 
  group_by(hour = as.numeric(str_sub(time, 1, 2))) %>% 
  summarise(
    walking = sum(Pedestrian),
    cycling = sum(Cyclist),
    passenger = sum(`Car occupant`)
    ) %>% 
  tidyr::gather(mode, casualties, -hour)

ggplot(crash_times, aes(hour, casualties)) +
  geom_line(aes(colour = mode))
```

<img src="man/figures/README-crash-time-plot-1.png" width="100%" />

Note that cycling manifests distinct morning and afternoon peaks (see
Lovelace, Roberts, and Kellar 2016 for more on this).

## Usage in research and policy contexts

The package has now been peer reviewed and is stable, and has been
published in the Journal of Open Source Software (Lovelace et al. 2019).
Please tell people about the package, link to it and cite it if you use
it in your work.

Examples of how the package can been used for policy making include:

-   Use of the package in a web app created by the library service of
    the UK Parliament. See
    [commonslibrary.parliament.uk](https://commonslibrary.parliament.uk/constituency-data-traffic-accidents/),
    screenshots of which from December 2019 are shown below, for
    details.

![](https://user-images.githubusercontent.com/1825120/70164249-bf730080-16b8-11ea-96d8-ec92c0b5cc69.png)

-   Use of methods taught in the
    [stats19-training](https://docs.ropensci.org/stats19/articles/stats19-training.html)
    vignette by road safety analysts at Essex Highways and the Safer
    Essex Roads Partnership ([SERP](https://saferessexroads.org/)) to
    inform the deployment of proactive front-line police enforcement in
    the region (credit: Will Cubbin).

-   Mention of road crash data analysis based on the package in an
    [article](https://www.theguardian.com/cities/2019/oct/07/a-deadly-problem-should-we-ban-suvs-from-our-cities)
    on urban SUVs. The question of how vehicle size and type relates to
    road safety is an important area of future research. A starting
    point for researching this topic can be found in the
    [`stats19-vehicles`](https://docs.ropensci.org/stats19/articles/stats19-vehicles.html)
    vignette, representing a possible next step in terms of how the data
    can be used.

## Next steps

There is much important research that needs to be done to help make the
transport systems in many cities safer. Even if you’re not working with
UK data, we hope that the data provided by **stats19** data can help
safety researchers develop new methods to better understand the reasons
why people are needlessly hurt and killed on the roads.

The next step is to gain a deeper understanding of **stats19** and the
data it provides. Then it’s time to pose interesting research questions,
some of which could provide an evidence-base in support policies that
save lives (e.g. Sarkar, Webster, and Kumari 2018). For more on these
next steps, see the package’s introductory
[vignette](https://itsleeds.github.io/stats19/articles/stats19.html).

## Further information

The **stats19** package builds on previous work, including:

-   code in the [bikeR](https://github.com/Robinlovelace/bikeR) repo
    underlying an academic paper on collisions involving cyclists
-   functions in [**stplanr**](https://docs.ropensci.org/stplanr/) for
    downloading Stats19 data
-   updated functions related to the
    [CyIPT](https://github.com/cyipt/stats19) project

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-lovelace_stats19_2019" class="csl-entry">

Lovelace, Robin, Malcolm Morgan, Layik Hama, Mark Padgham, and M
Padgham. 2019. “Stats19 A Package for Working with Open Road Crash
Data.” *Journal of Open Source Software* 4 (33): 1181.
<https://doi.org/10.21105/joss.01181>.

</div>

<div id="ref-lovelace_who_2016" class="csl-entry">

Lovelace, Robin, Hannah Roberts, and Ian Kellar. 2016. “Who, Where,
When: The Demographic and Geographic Distribution of Bicycle Crashes in
West Yorkshire.” *Transportation Research Part F: Traffic Psychology and
Behaviour*, Bicycling and bicycle safety, 41, Part B.
<https://doi.org/10.1016/j.trf.2015.02.010>.

</div>

<div id="ref-sarkar_street_2018" class="csl-entry">

Sarkar, Chinmoy, Chris Webster, and Sarika Kumari. 2018. “Street
Morphology and Severity of Road Casualties: A 5-Year Study of Greater
London.” *International Journal of Sustainable Transportation* 12 (7):
510–25. <https://doi.org/10.1080/15568318.2017.1402972>.

</div>

</div>
# stats19 2.0.0 2020-10

* Major changes to the datasets provided by the DfT have led to major changes to the package. See (#212) for details.
* To reduce code complexity the package no longer supports reading in multiple years
* This puts the onus on the user of the package to understand the input data, rather than relying on clever coding to join everything together. Note: you can easily join different years, e.g. with the command `purrr::map_dfr()`.


# stats19 1.5.0 2021-10

* Support new https download links (#208)
* Package tests now pass when wifi is turned off
* URLs have been fixed

# stats19 1.4.3 2021-07-21

* Use 1st edition of `readr` on Windows to prevent errors on reading data (#205)

# stats19 1.4.2 2021-07

* Fix CRAN checks associated with access to online resources (#204)
* [Fix](https://github.com/ropensci/stats19/commit/826a1d0ed3b9fbcf80675b64fd5731ae8b7b0498) issues associated with `get_ULEZ()` and `get_MOT()` functions


# stats19 1.4.1

* New function `get_ULEZ()` to get data on vehicles from a number plate (thanks to Ivo Wengraf)
* Added a test to prevent rare failures in `get_stats19()` when `data_dir` points to the working directory


# stats19 1.4.0

* Add `get_stats19_adjustments()` function
* Use GH Actions for CI (#177)
* Fixed a problem with `get_stats19()` and multiple years that could be linked with the same data file (#168)
* Fix issues with vignettes for CRAN (#190)

# stats19 1.3.0

* Support for 2019 data (#171)

# stats19 1.2.0

* Tests now pass on the development version of R (4.0.0)
* The package now has a hex sticker! See https://github.com/ropensci/stats19/issues/132 for discussion
* The output of formatted crash datasets gains a new column, `datetime` that is a properly formatted date-time (`POSIXct`) object in the correct timezone (`Europe/London`) (#146)
* Enables the download of multiple years as per https://github.com/ropensci/stats19/issues/99, thanks to Layik Hama
* Users can now set the default data download directory with STATS19_DOWNLOAD_DIRECTORY=/path/to/data in your .Renviron file: https://github.com/ropensci/stats19/issues/141
* `get_stats19()` gains a new argument `output_format()` that enables results to be returned as an `sf` object or a `ppp` object for use the the `spatstat` package thanks to work by Andrea Gilardi https://github.com/ropensci/stats19/pull/136

# stats19 1.1.0

* Now enables the download of 2018 data
* Various bug fixes, see https://github.com/ropensci/stats19/issues
* Update website link: https://docs.ropensci.org/stats19/
* New work-in-progress vignette on vehicles data: https://docs.ropensci.org/stats19/articles/stats19-vehicles.html

# stats19 1.0.0

* Major change to `dl_stats19()`: it is now much easier to download STATS19 data. By default `ask = FALSE` in `get_stats19()` and `dl_stats19()`.

# stats19 0.2.1

* Fixed issue with column labels not being there - see [#82](https://github.com/ropensci/stats19/issues/92)

# stats19 0.2.0

* `get_stats19()` gains an `ask` argument (`TRUE` by default, set as `FALSE` to make road crash data access even more automated!)
* The `date` column now is of the correct class after formatting `POSIXct`. See [#86](https://github.com/ropensci/stats19/issues/86)
* Added a `NEWS.md` file to track changes to the package.
This document was adapted from @peterdesmet's [template](https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c)

# Contributing to stats19

First of all, thanks for considering contributing to stats19! 👍 It's people like you that make it rewarding for us - the project maintainers - to work on stats19. 😊

stats19 is an open source project, maintained by people who care. We are not directly funded to do so.

[repo]: https://github.com/ITSLeeds/stats19
[issues]: https://github.com/ITSLeeds/stats19/issues
[new_issue]: https://github.com/ITSLeeds/stats19/issues/new
[website]: https://ITSLeeds.github.io/stats19
[citation]: https://ITSLeeds.github.io/stats19/authors.html
[email]: maintainer_email

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Share the love ❤️

Think stats19 is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using stats19 for a paper you are writing? Consider [citing it][citation].

### Ask a question ⁉️

Using stats19 and got stuck? Browse the [documentation][website] to see if you can find a solution. Still stuck? Post your question as an [issue on GitHub][new_issue]. While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [email][email].

### Propose an idea 💡

Have an idea for a new stats19 feature? Take a look at the [documentation][website] and [issue list][issues] to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub][new_issue]. While we can't promise to implement your idea, it helps to:

* Explain in detail how it would work.
* Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a bug 🐛

Using stats19 and discovered a bug? That's annoying! Don't let others have the same experience and report it as an [issue on GitHub][new_issue] so we can fix it. A good bug report makes it easier for us to do so, so please include:

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

Care to fix bugs or implement new functionality for stats19? Awesome! 👏 Have a look at the [issue list][issues] and leave a comment on the things you want to work on. See also the development guidelines below.

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
I just realised that the errant URL was still there in the previous submission.

Apologies, this one should fix it! Many thanks.


## Test environments
* local R installation, R 4.1.1
* ubuntu 20.04 (GitHub Actions), R 4.1.1
* win-builder (devel)

Your package stats19_2.0.0.tar.gz has been built (if working) and checked for Windows.
Please check the log files and (if working) the binary package at:
https://win-builder.r-project.org/GFSUVmV38fB2
The files will be removed after roughly 72 hours.
Installation time in seconds: 6
Check time in seconds: 241
Status: 1 NOTE
R Under development (unstable) (2021-10-18 r81071)

## R CMD check results

0 errors | 0 warnings | 0 notes
This folder contains the scripts to create the `stats19_schema` and
`stats19_variables` data, which contain metadata on stats19 data.
`stats19_schema` is a look-up table matching codes provided in the raw stats19
dataset with character strings.
---
output: github_document
bibliography: vignettes/references.bib
---

<!-- badges: start -->
<!-- [![Travis build status](https://travis-ci.org/ropensci/stats19.svg?branch=master)](https://travis-ci.org/ropensci/stats19) -->
[![](http://www.r-pkg.org/badges/version/stats19)](https://www.r-pkg.org/pkg/stats19)
[![R-CMD-check](https://github.com/ropensci/stats19/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/stats19/actions)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/stats19)](https://www.r-pkg.org/pkg/stats19)
[![Life cycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![](https://badges.ropensci.org/266_status.svg)](https://github.com/ropensci/software-review/issues/266)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.01181/status.svg)](https://doi.org/10.21105/joss.01181)
![codecov](https://codecov.io/gh/ropensci/stats19/branch/master/graph/badge.svg)
<!-- badges: end -->

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2540781.svg)](https://doi.org/10.5281/zenodo.2540781) -->
<!-- [![Gitter chat](https://badges.gitter.im/ITSLeeds/stats19.png)](https://gitter.im/stats19/Lobby?source=orgpage) -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# stats19 <a href='https://docs.ropensci.org/stats19/'><img src='https://raw.githubusercontent.com/ropensci/stats19/master/man/figures/logo.png' align="right" height=215/></a>

**stats19** provides functions for downloading and formatting road crash data.
Specifically, it enables access to the UK's official road traffic casualty database, [STATS19](https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data). (The name comes from the form used by the police to record car crashes and other incidents resulting in casualties on the roads.)

A full overview of STATS19 variables be found in a [document](https://data.dft.gov.uk/road-accidents-safety-data/Brief-guide-to%20road-accidents-and-safety-data.doc) provided by the UK's Department for Transport (DfT).

The raw data is provided as a series of `.csv` files that contain integers and which are stored in dozens of `.zip` files.
Finding, reading-in and formatting the data for research can be a time consuming process subject to human error.
**stats19** speeds up these vital but boring and error-prone stages of the research process with a single function: `get_stats19()`.
By allowing public access to properly labelled road crash data, **stats19** aims to make road safety research more reproducible and accessible.

For transparency and modularity, each stage can be undertaken separately, as documented in the [stats19 vignette](https://itsleeds.github.io/stats19/articles/stats19.html).

## Installation

Install and load the latest version with:

```{r eval=FALSE}
remotes::install_github("ropensci/stats19")
```

```{r attach}
library(stats19)
```

You can install the released version of stats19 from [CRAN](https://cran.r-project.org/package=stats19) with:

```r
install.packages("stats19")
```

## get_stats19()

`get_stats19()` requires `year` and `type` parameters, mirroring the provision of STATS19 data files, which are categorised by year (from 1979 onward) and type (with separate tables for crashes, casualties and vehicles, as outlined below).
The following command, for example, gets crash data from 2017 (**note**: we follow the "crash not accident" campaign of [RoadPeace](https://www.roadpeace.org/take-action/crash-not-accident/) in naming crashes, although the DfT refers to the relevant tables as 'accidents' data):

```{r}
crashes = get_stats19(year = 2017, type = "accident")
```

What just happened?
For the `year` 2017 we read-in crash-level (`type = "accident"`) data on all road crashes recorded by the police across Great Britain.
The dataset contains `r ncol(crashes)` columns (variables) for `r format(nrow(crashes), big.mark = ",")` crashes.
We were not asked to download the file (by default you are asked to confirm the file that will be downloaded).
The contents of this dataset, and other datasets provided by **stats19**, are outlined below and described in more detail in the [stats19 vignette](https://itsleeds.github.io/stats19/articles/stats19.html).

We will see below how the function also works to get the corresponding casualty and vehicle datasets for 2017.
The package also allows STATS19 files to be downloaded and read-in separately, allowing more control over what you download, and subsequently read-in, with `read_accidents()`, `read_casualties()` and `read_vehicles()`, as described in the vignette.


## Data download

Data files can be downloaded without reading them in using the function `dl_stats19()`.
If there are multiple matches, you will be asked to choose from a range of options.
Providing just the year, for example, will result in the following options:

```{r dl2017-all}
dl_stats19(year = 2020, data_dir = tempdir())
```

```
Multiple matches. Which do you want to download?

1: dft-road-casualty-statistics-casualty-2020.csv
2: dft-road-casualty-statistics-vehicle-2020.csv
3: dft-road-casualty-statistics-accident-2020.csv
4: dft-road-casualty-statistics-vehicle-e-scooter-2020.csv
5: road-safety-data-mid-year-provisional-estimates-2020.csv

Selection: 
Enter an item from the menu, or 0 to exit
```

## Using the data

STATS19 data consists of 3 main tables:

- Accidents, the main table which contains information on the crash time, location and other variables (`r ncol(accidents_sample)` columns in total)
- Casualties, containing data on people hurt or killed in each crash (`r ncol(casualties_sample)` columns in total)
- Vehicles, containing data on vehicles involved in or causing each crash (`r ncol(vehicles_sample)` columns in total)

The contents of each is outlined below.

### Crash data

Crash data was downloaded and read-in using the function `get_stats19()`, as described above.

```{r read2017-raw-format}
nrow(crashes)
ncol(crashes)
```

Some of the key variables in this dataset include:

```{r crashes2017-columns}
key_column_names = grepl(pattern = "severity|speed|pedestrian|light_conditions", x = names(crashes))
crashes[key_column_names]
```

For the full list of columns, run `names(crashes)` or see the [vignette](https://github.com/ropensci/stats19/blob/master/vignettes/stats19.Rmd).

<!-- This means `crashes` is much more usable than `crashes_raw`, as shown below, which shows three records and some key variables in the messy and clean datasets: -->

### Casualties data

As with `crashes`, casualty data for 2017 can be downloaded, read-in and formatted as follows:

```{r 2017-cas}
casualties = get_stats19(year = 2017, type = "casualty", ask = FALSE)
nrow(casualties)
ncol(casualties)
```

The results show that there were `r format(nrow(casualties), big.mark=",")` casualties reported by the police in the STATS19 dataset in 2017, and `r ncol(casualties)` columns (variables).
Values for a sample of these columns are shown below:

```{r 2017-cas-columns}
casualties[c(4, 5, 6, 14)]
```

The full list of column names in the `casualties` dataset is:

```{r 2017-cas-columns-all}
names(casualties)
```

### Vehicles data

Data for vehicles involved in crashes in 2017 can be downloaded, read-in and formatted as follows:

```{r dl2017-vehicles}
vehicles = get_stats19(year = 2017, type = "vehicle", ask = FALSE)
nrow(vehicles)
ncol(vehicles)
```

The results show that there were `r format(nrow(vehicles), big.mark=",")` vehicles involved in crashes reported by the police in the STATS19 dataset in 2017, with `r ncol(vehicles)` columns (variables).
Values for a sample of these columns are shown below:

```{r 2017-veh-columns}
vehicles[c(3, 14:16)]
```

The full list of column names in the `vehicles` dataset is:

```{r 2017-veh-columns-all}
names(vehicles)
```

## Creating geographic crash data

An important feature of STATS19 data is that the "accidents" table contains geographic coordinates.
These are provided at ~10m resolution in the UK's official coordinate reference system (the Ordnance Survey National Grid, EPSG code 27700).
**stats19** converts the non-geographic tables created by `format_accidents()` into the geographic data form of the [`sf` package](https://cran.r-project.org/package=sf) with the function `format_sf()` as follows:

```{r format-crashes-sf}
crashes_sf = format_sf(crashes)
```

The note arises because `NA` values are not permitted in `sf` coordinates, and so rows containing no coordinates are automatically removed.
Having the data in a standard geographic form allows various geographic operations to be performed on it.
The following code chunk, for example, returns all crashes within the boundary of West Yorkshire (which is contained in the object [`police_boundaries`](https://itsleeds.github.io/stats19/reference/police_boundaries.html), an `sf` data frame containing all police jurisdictions in England and Wales). 

```{r crashes-leeds}
library(sf)
library(dplyr)
wy = filter(police_boundaries, pfa16nm == "West Yorkshire")
crashes_wy = crashes_sf[wy, ]
nrow(crashes_sf)
nrow(crashes_wy)
```

This subsetting has selected the `r format(nrow(crashes_wy), big.mark = ",")`
crashes which occurred within West Yorkshire in 2017.


## Joining tables

The three main tables we have just read-in can be joined by shared key variables.
This is demonstrated in the code chunk below, which subsets all casualties that took place in Leeds, and counts the number of casualties by severity for each crash:

```{r table-join, message = FALSE}
sel = casualties$accident_index %in% crashes_wy$accident_index
casualties_wy = casualties[sel, ]
names(casualties_wy)
cas_types = casualties_wy %>%
  select(accident_index, casualty_type) %>%
  mutate(n = 1) %>%
  group_by(accident_index, casualty_type) %>%
  summarise(n = sum(n)) %>%
  tidyr::spread(casualty_type, n, fill = 0)
cas_types$Total = rowSums(cas_types[-1])
cj = left_join(crashes_wy, cas_types, by = "accident_index")
```

What just happened? We found the subset of casualties that took place in West Yorkshire with reference to the `accident_index` variable.
Then we used functions from the **tidyverse** package **dplyr** (and `spread()` from **tidyr**) to create a dataset with a column for each casualty type.
We then joined the updated casualty data onto the `crashes_wy` dataset.
The result is a spatial (`sf`) data frame of crashes in Leeds, with columns counting how many road users of different types were hurt.
The original and joined data look like this:

```{r table-join-examples}
crashes_wy %>%
  select(accident_index, accident_severity) %>% 
  st_drop_geometry()
cas_types[1:2, c("accident_index", "Cyclist")]
cj[1:2, c(1, 5, 34)] %>% st_drop_geometry()
```

## Mapping crashes

The join operation added a geometry column to the casualty data, enabling it to be mapped (for more advanced maps, see the [vignette](https://itsleeds.github.io/stats19/articles/stats19.html)):

```{r}
cex = cj$Total / 3
plot(cj["speed_limit"], cex = cex)
```

The spatial distribution of crashes in West Yorkshire clearly relates to the region's geography.
Crashes tend to happen on busy Motorway roads (with a high speed limit, of 70 miles per hour, as shown in the map above) and city centres, of Leeds and Bradford in particular.
The severity and number of people hurt (proportional to circle width in the map above) in crashes is related to the speed limit.

STATS19 data can be used as the basis of road safety research.
The map below, for example, shows the results of an academic paper on the social, spatial and temporal distribution of bike crashes in West Yorkshire, which estimated the number of crashes per billion km cycled based on commuter cycling as a proxy for cycling levels overall (more sophisticated measures of cycling levels are now possible thanks to new data sources) [@lovelace_who_2016]:

```{r, echo=FALSE}
knitr::include_graphics("https://ars.els-cdn.com/content/image/1-s2.0-S136984781500039X-gr9.jpg")
```

## Time series analysis

We can also explore seasonal trends in crashes by aggregating crashes by day of the year:

```{r crash-date-plot}
library(ggplot2)
crashes_dates = cj %>% 
  st_set_geometry(NULL) %>% 
  group_by(date) %>% 
  summarise(
    walking = sum(Pedestrian),
    cycling = sum(Cyclist),
    passenger = sum(`Car occupant`)
    ) %>% 
  tidyr::gather(mode, casualties, -date)
ggplot(crashes_dates, aes(date, casualties)) +
  geom_smooth(aes(colour = mode), method = "loess") +
  ylab("Casualties per day")
```


Different types of crashes also tend to happen at different times of day.
This is illustrated in the plot below, which shows the times of day when people who were travelling by different modes were most commonly injured.

```{r crash-time-plot}
library(stringr)

crash_times = cj %>% 
  st_set_geometry(NULL) %>% 
  group_by(hour = as.numeric(str_sub(time, 1, 2))) %>% 
  summarise(
    walking = sum(Pedestrian),
    cycling = sum(Cyclist),
    passenger = sum(`Car occupant`)
    ) %>% 
  tidyr::gather(mode, casualties, -hour)

ggplot(crash_times, aes(hour, casualties)) +
  geom_line(aes(colour = mode))
```

Note that cycling manifests distinct morning and afternoon peaks [see @lovelace_who_2016 for more on this].

## Usage in research and policy contexts

The package has now been peer reviewed and is stable, and has been published in the Journal of Open Source Software [@lovelace_stats19_2019].
Please tell people about the package, link to it and cite it if you use it in your work.

Examples of how the package can been used for policy making include:

- Use of the package in a web app created by the library service of the UK Parliament. See [commonslibrary.parliament.uk](https://commonslibrary.parliament.uk/constituency-data-traffic-accidents/), screenshots of which from December 2019 are shown below, for details.

![](https://user-images.githubusercontent.com/1825120/70164249-bf730080-16b8-11ea-96d8-ec92c0b5cc69.png)

- Use of methods taught in the [stats19-training](https://docs.ropensci.org/stats19/articles/stats19-training.html) vignette by road safety analysts at Essex Highways and the Safer Essex Roads Partnership ([SERP](https://saferessexroads.org/)) to inform the deployment of proactive front-line police enforcement in the region (credit: Will Cubbin).

- Mention of road crash data analysis based on the package in an [article](https://www.theguardian.com/cities/2019/oct/07/a-deadly-problem-should-we-ban-suvs-from-our-cities) on urban SUVs. 
The question of how vehicle size and type relates to road safety is an important area of future research.
A starting point for researching this topic can be found in the [`stats19-vehicles`](https://docs.ropensci.org/stats19/articles/stats19-vehicles.html) vignette, representing a possible next step in terms of how the data can be used.


## Next steps

There is much important research that needs to be done to help make the transport systems in many cities safer.
Even if you're not working with UK data, we hope that the data provided by **stats19** data can help safety researchers develop new methods to better understand the reasons why people are needlessly hurt and killed on the roads.

The next step is to gain a deeper understanding of **stats19** and the data it provides.
Then it's time to pose interesting research questions, some of which could provide an evidence-base in support policies that save lives [e.g. @sarkar_street_2018].
For more on these next steps, see the package's introductory [vignette](https://itsleeds.github.io/stats19/articles/stats19.html).

## Further information

The **stats19** package builds on previous work, including:

- code in the [bikeR](https://github.com/Robinlovelace/bikeR) repo underlying an academic paper on collisions involving cyclists
- functions in [**stplanr**](https://docs.ropensci.org/stplanr/) for downloading Stats19 data
- updated functions related to the [CyIPT](https://github.com/cyipt/stats19) project

 [![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

## References

---
output: github_document
---

## Responses to review 1 of stats19

Thanks for the review.
We've had a chance, after making some changes and fixes to the package, to take-in and act on each of the comments.
The code-base has evolved substantially since the review, but the fundamental design of the package, with its 3 stage API mirroring workflows that happened before the package was developed, remains unchanged.
That is:

- `dl_stats19()` downloads files from the DfT. Good news: we have spoken to the relevant people at the Department for Transport and they assured us that the endpoints are stable. The function now uses `menu()` to provide a menu of download options for any year/type combinations and now finds files outside those explicitly mentioned in the file names.
E.g.:

```{r, eval=FALSE}
dl_stats19(year = 2017)
# Multiple matches. Which do you want to download?
# 
# 1: dftRoadSafetyData_Vehicles_2017.zip
# 2: dftRoadSafetyData_Casualties_2017.zip
# 3: dftRoadSafetyData_Accidents_2017.zip
dl_stats19(year = 2017, type = "ac")
# Files identified: dftRoadSafetyData_Accidents_2017.zip
# 
# Wanna do it (y = enter, n = esc)? 
dl_stats19(year = 1985)
# Year not in range, changing to match 1979:2004 data
# This file is over 240 MB in size.
# Once unzipped it is over 1.8 GB.
# Files identified: Stats19-Data1979-2004.zip
# 
# Wanna do it (y = enter, n = esc)?
```

- `read_*()` these functions remain unchanged, except the order of arguments has changed.
Like `dl_stats19()`, `year` is now the first argument, which is more intuitive.

- `format_*()` functions have been refactored. Each now uses `format_stats19()` behind the scenes reducing duplication.
The results are now better: more variables are now labelled.

We'll focus on areas flagged in the review for the rest of this response:

> I would tease a bit more of what's in these data sets. I wasn't entirely sure until I downloaded and opened the supporting documentation. If I were searching for this kind of data, and I didn't know what STATS19 was, I'd like to know I'm in the right place after scanning the README. Maybe a map?

We have added a map (well technically 9 maps!) and a couple of time series plots showing the scale of the data.
Also show a sample of the additional casualty and vehicle tables has been added to show more clearly the richness of data provided.

> I couldn't load the vignette from the console:

We also could not see the vignette when installing using `devtools::install_github(build_vignettes = TRUE`. But we can see the vignette if we install locally.

This was the code we ran:

```{r}
devtools::install(build_vignettes = TRUE)
vignette(package = "stats19")
```

> Several of the examples failed:

These have now been fixed - thanks for testing and reporting.

> I couldn't find any explicit contributing guidelines in the README, and there is no CONTRIBUTING document.

A CONTRIBUTING is added now. Thank you.

> The package has an obvious research application according to JOSS's definition

> There is no paper.md.

One is added with:

- A short summary describing the high-level functionality of the software
- Authors: A list of authors with their affiliations
- A statement of need clearly stating problems the software is designed to solve and its target audience.
- References: with DOIs for all those that have one (e.g. papers, datasets, software).

Review Comments

> A superb and essential package--we need this data and we need it in these formats. The download-format-read-explore workflow is intuitive and relatively frictionless. I have only some brief comments:

Thank you.

> I wonder you could possibly merge the formatting and reading step with a raw = TRUE or format = TRUE argument in the read_* functions. But perhaps that's my tendency towards abstraction. Something like ac = read_accidents(year = 2017, format = TRUE)

Done, appreciate your input.

> My personal preference would be to have the schema from dl_schema lazily loaded with the package.

DESCRIPTION: has the line LazyData which means stats19_schema is lazy loaded.

> According to the vignette, the dl_* functions are interactive, although the interactivity is commented out in the code. Will the interactivity be returning? Or does the vignette need to be updated?

Back in, as stated above.

> Out of curiosity, what's happening with https://github.com/cyipt/stats19? It was updated recently.

@mem48 answered this: cyipt/stats19 is not actually a proper R package. It is a repo containing scripts for CyIPT project, it has different sources (UK DS), and usage so there is no current need to adapt the use to this package. Malcolm is one of the contributors to this package.

> I confess I wish the package name was more expressive--stats19 sounds like an introductory statistics class.

This a reasonable point that we have thought of and discussed.
We are open minded about changing the name but, as with so many things, there are +s and -s (outlined for some options below):

- **stats19data**
  - + clarifies that it's about data
  - - longer, suffers from some of the same issues that **stats19** suffers from, the package is more about data formatting than data provision

- **roadcrashesUK**
  - + explicit, makes region of data access transparent
  - - there are other types of road crash data, also the data currently provided is technically for Great Britain, but **roadcrashesGB** doesn't work so well, and we may want to add data access options for Northern Ireland at some point also

- **roadSafetyData**
  - + Matches DfT's [webpage](https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data) title on the topic
  - - longer and, again, is less specific.
  
The main benefit we can see of changing the name would be making the package easier to find.
We think good documentation and clear description and some write-ups of the package and what it does could address these issues.
We've explored **stat19** name and it links directly to (and is almost synonymous with) road crash data.
See https://en.wikipedia.org/wiki/STATS19 for an excellent example (we plan to add this link to the README)

so the name is OK for we think, but we're open minded to alternative names mentioned above and perhaps names we've not thought of.


> This data will be used to make many maps. I personally would love a nudge in that direction in either the README or the vignette.

Definitely. Thank you very much for your input.
---
output: github_document
---

## Responses to review 2 of stats19

Thanks for the detailed review.
We think all the suggestions make sense, and link to previous discussions about combining the 3 stage (`dl`, `read`, `format`) process into a single function call, which we were thinking of calling `get_stats19()`.
There reason for splitting the process up is to ensure maximum transparency and to give the user control over what the package is doing.
However, as long as it is properly documents, we think the benefits of a `get_stats19()` function will outweigh any possible negatives we can think of so we plan to go ahead and create this function.

We would very much welcome a pull request that address some of the other issues you mention.

Responses to other issues/questions/comments are provided below.
We suspect that some of the comments refer to an older version of the package, which is completely understandable (the package has evolved since the initial submission!) and explains why some of the responses are short / questions.

- [x]Clearly much effort has gone into this package. I greatly support the sentiment behind making these data available in R having done the same for other data sets, myself. This package should be a great benefit to researchers in this area. I really appreciate the slim dependencies. It will make this package much easier to maintain into the future.

Thanks for the comments, we have indeed tried to keep dependencies to a minimum but consider `readr` and `tibble` worthwhile.
`readxl` and `curl` have been demoted to `Suggests`, as detailed in another comment.

- [x]I found the package to be well documented, the vignette is helpful in illustrating some use cases for the data along with how to access it and the code is clear.

Thanks. If you think of other was we can communicate the value of the data, do let us know (I think the second mapping figure could be improved...).

- [ ] Some of the functionality I find to be mildly confusing like downloading and then importing the files into the R session and then formatting. As a user I'd prefer it all in one step, but there are likely use cases I'm not aware of that mean that this is useful so some examples of this would be useful I think.

We have long been planning to add a `get_stats19()` function as per https://github.com/ITSLeeds/stats19/issues/11 
The review comment, combined with further discussion, has triggered us to re-prioritise it.
It's been beneficial to polish each of the component functions first, however, and good to document each stage for maximum transparency, however, so we plan to keep the `dl`, `read` and `format` functions exported.

- [ ] My general comments on the code follow and sections for the DESCRTIPTION and Vignette as well. I've commented quite a bit on grammar and spelling as I think that the polish of a package is important as it lends to the perception of the quality.

Agreed.

--------------------------------------------------------------------------------

- [ ] Per rOpenSci policy, avoid start-up messages, rely on the documentation for
citation information:
https://ropensci.github.io/dev_guide/building.html#recommended-scaffolding.

The guidance is to 'Only use package startup messages when necessary'.
A case can be made that this is necessary.
As with `osmdata`, the package provides access to data that has a license that requires it to be cited.
The `osmdata` load message is as follows:

```{r}
library(osmdata)
```

We fully agree with the reasoning behind remove package startup messages however.
As a compromise, we've shortened the startup from 4 lines to 2:

```{r}
# before:
# Data provided under the conditions of the Open Government License.
# If you use data from this package, mention the source
# (UK Department for Transport), cite the package, and link to:
# www.nationalarchives.gov.uk/doc/open-government-licence/version/3/.
```

```{r}
# after:
library(stats19)
```

- [ ] Avoid long lines >80 chars

running `goodpractice::gp()` found the following lines with > 80 lines:

```
    R/format.R:62:1
    R/format.R:67:1
    R/read.R:141:1
    R/utils.R:167:1
```

All these have been fixed.

- [ ] Inconsistent use of white spaces in code, see `find_file_name()` `if` statements for examples.

- [ ] The package does not pass `R CMD check`. *curl*, *readxl* and *tibble* are all listed as an Imports in DESCRIPTION but not imported from. With *curl* being used in tests, this means it should be in Suggests, I think. The others should be removed.

`curl` is used in the tests and `readxl` is used in the examples.
These have been demoted to `Suggests`.
`tibble` has been removed from the DESCRIPTION file.

- [ ] I don't think it's good form to have an example that won't work on Windows in the help file for `stats19_schema`, from data.R - line 17? Most of what I see there would be better served in a `data_raw` folder showing how the data were created with the documentation actually documenting what the variables are not how they were created, see <http://r-pkgs.had.co.nz/data.html> and for an example, <https://github.com/ropensci/GSODR/tree/master/data-raw>.

- [ ] I would suggest to use proper formatting in help files, when naming packages, e.g. \pkg{stats19} and when referring to documented functions or data, e.g. \code{\link{stats19_schema}}, or with single quotes around abbreviations, e.g. 'DfT'. @ColinFay has an excellent page that outlines the formatting options and when/how to use them, <https://colinfay.me/writing-r-extensions/writing-r-documentation-files.html>. This will greatly enhance the users' experience when using the help files by making them more readable.

- [ ] I also would suggest making use of `@seealso` in documentation. For example, the `dl_stats19()` example works great in the help files, but from there I have the data but it's not in R. Using the `@seealso` you can let the user know about the `read_*()` functions.

- [ ] I downloaded files using `dl_stats19()`, selecting "Casualties", and then ran `read_accidents()` and got 
Is it possible to be more descriptive and say that I've used the wrong `read_*()` based on the file/data found and offer to import it?

- [ ] Missing "." after "e.g." in dl.R on lines 8 and 9, there may be others that I didn't spy.

- [ ] Capitalisation in help files is inconsistent, e.g. lines 123-125 of read.R, parameter descriptions are mixed upper and lower case for first word after parameter itself. This applies to other functions where the descriptions are given in all lower case for other functions or upper case.

- [ ] Testing the functionality, I get this, when I expect it to tell me that `deaths` is not a valid input. But then when I hit escape, I expect it simply exit, not provide a warning message on the way out as well.
```r
dl_stats19(year = 1979, type = "deaths")
No files of that type found for that year.
This will download 240 MB+ (1.8 GB unzipped).
Files identified: Stats19-Data1979-2004.zip

Download now (y = enter, n = esc)? 

Warning message:
In find_file_name(years = year, type = type) :
  Coordinates unreliable in this data.
```

- [ ] I got caught out when using the interactive features. I read "y = enter" but hit "y" thinking that would work as well as hitting "enter", but R cancelled the operation anyway just as if I'd hit "esc"

- [ ] Per a recent conversation with CRAN, you should use `donttest()` rather than `dontrun()` for examples you don't want to be run on CRAN. Then set .travis.yml to run them by using `r_check_args: --as-cran --run-donttest`. **This may not be appropriate in all cases, e.g. interactive functions.**


- [ ] When validating user inputs and using `stop()` it's nice to use `call. = FALSE` to simplify the error message that the user receives.

- [ ] Consider using [`hoardr`](https://ropensci.github.io/hoardr/) for managing user-saved files on disk that aren't in `tempdir()`?

- [ ] When using `utils::download.file()`, you should use `mode = "wb"` or Windows users may end up with corrupted downloads in my experience. `curl::curl_download()` does the same thing but uses more updated ways of doing it and defaults to using a binary mode (wb).

- [ ] I don't think that there is much need for the `Attempting download from` or `Reading in: ` message. If it takes that long, I would suggest to use a progress bar to show progress. But this is just a personal observation.

- [ ] Consider setting up a `pkgdown` site? It's easy to do and you can automate deployment with your Travis-CI so it's less to remember.

#### Tests

- [ ] I'm unclear how the interactive portion of the package functions is handled in testing? There are ways to handle this, but I don't see any implemented and when I run `devtools::test()` I'm asked to provide my own input.

- [ ] Suggest using `skip_on_cran()` since some of the tests can take some time to execute due to download times.

#### DESCRIPTION File

- [ ] In the DESCRIPTION file, Mark's author entry is missing his ORCID.

- [ ] More information in the DESCRIPTION's Description field would be desirable, a link to the data's website or other information to give more background perhaps.

- [ ] STATS19 should be in "'" in DESCRIPTION for CRAN, i.e., 'STATS19', I think.

- [ ] Check spelling in DESCRIPTION file, see: "analysie"

- [ ] The Description should include a link to the DfT website.

- [ ] Language field should be set, `Language: en-GB`

#### README File(s)

- [ ] Use `remotes::install_github()` in place of `devtools::install_github()` in README.

- [ ] The code style is inconsistent in the README.Rmd file in the code chunks, e.g. line 85 is missing space around `=`.

- [ ] The example in the README showing two steps seems necessarily confusing to new users. If there is a good reason for having the raw data in R, document in a vignette why this is useful and show the two-step process, but if normal users won't do this, I wouldn't show it in the quick-start.

- [ ] Line 43 of README uses inconsistent "(" around the phrases with the other two `read_*` function description.

#### Vignette

- [ ] Run spell-check on it.

- [ ] The term "attach"" has a specific meaning in R. Suggest rewording the portion about installation and loading the package to omit the use of "attach", since you're not using `attach()` in the R sense (and really shouldn't use it anyway).

- [ ] I would describe why a user might want or need to install the Development version from GitHub in the vignette. Presumably if they are reading the vignette, they've already installed the package from CRAN (in the future).

- [ ] Try to consistently use `function()` to identify functions in the vignette text. This also means that if/when you use pkgdown to build a site, the functions are linked to the help file.

- [ ] In the introduction, the description of why there are `read_*()` and `format_*()` functions is confusing. To me, it reads as if `format` is only a parameter for `read_*()` in the introduction. I was left wondering why it's documented there or why the `format_*()`s even exist until I reached the end of the vignette.

- [ ] There is a comma out of place in Vignette,

- [ ] Format: Each of the read_*() functions has a format parameter which, when TRUE, adds

should be 

- [ ] Format: Each of the read_*() functions has a format parameter, which, when TRUE, adds 

- [x] I'm unsure about including a package that's not on CRAN in the vignette (`ukboundaries`), something like this should be listed in Suggests, but it's not on CRAN, @sckott do you have any thoughts?

This is a good point.
Fixed, by adding a much more useful dataset, representing the juristictions of polic forces across England and Wales. 

- [ ] The first figures in the `sf` section after the join aren't immediately clear to me. The axis lack labels, I'm not really sure what I'm looking at.

#### Meta

- [ ] The contributing guidelines mention a `pkgdown` website, this does not exist



This function generates the data object `stats19_schema`.

The function also generates `stats19_variables` (see the function's source code
for details).

```{r}
library(tidyverse)
```


# Load stats19 schema and save variable names

```{r}
schema_url = "https://data.dft.gov.uk/road-accidents-safety-data/Road-Safety-Open-Dataset-Data-Guide.xlsx"
schema_f = basename(schema_url)
schema_saved = file.path(stats19::get_data_directory(), schema_f)
download.file(schema_url, destfile = schema_saved)
schema_dft = readxl::read_excel(schema_saved)
schema_dft
stats19_variables
readr::write_csv(stats19_variables, "data-raw/stats19_variables.csv")
stats19_variables_dft = schema_dft %>% 
  rename(variable = `field name`) %>% 
  group_by(table, variable) %>% 
  summarise(
    note = first(note)
    ) 
stats19_variables_dft
summary(stats19_variables$table %in% stats19_variables_dft$table)
summary(in_original <- stats19_variables$column_name %in% stats19_variables_dft$variable)
stats19_variables$column_name[!in_original]
# [1] "latitude"                    "was_vehicle_left_hand_drive"
summary(in_new <- stats19_variables_dft$variable %in% stats19_variables$column_name)
stats19_variables_dft$variable[!in_new]
#  [1] "accident_reference"           "accident_year"                "Latitude"                    
#  [4] "local_authority_ons_district" "trunk_road_flag"              "accident_reference"          
#  [7] "accident_year"                "accident_reference"           "accident_year"               
# [10] "vehicle_text"                 "accident_reference"           "accident_year"               
# [13] "generic_make_model"           "vehicle_direction_from"       "vehicle_direction_to"        
# [16] "vehicle_left_hand_drive"   
stats19_variables_dft$column_name = snakecase::to_snake_case(stats19_variables_dft$variable)
stats19_variables_minimal = stats19_variables %>% 
  select(column_name, type)
stats19_variables_joined = left_join(stats19_variables_dft, stats19_variables_minimal)
table(stats19_variables_joined$type)
# character      date  location   numeric     other      time 
#        40         1         3         9        28         1 
stats19_variables_joined %>% 
  filter(is.na(type))
stats19_variables_joined$type = "character"
stats19_variables_joined$type[
  grepl(pattern = "year", x = stats19_variables_joined$column_name) 
] = "numeric"

```

# Save the schema

```{r}
stats19_schema
readr::write_csv(stats19_schema, "data-raw/stats19_schema.csv")
table(schema_dft$`code/format`)
stats19_schema_dft = schema_dft %>% 
  rename(code = `code/format`, variable = `field name`)
stats19_schema_joined = left_join(stats19_schema_dft, stats19_variables_joined) 
stats19_schema_joined = stats19_schema_joined %>% 
  rename(variable_formatted = column_name) %>% 
  filter(!is.na(as.numeric(code)))
```

# Tests

```{r}
s = stats19_schema
s
s %>% 
  filter(variable == "vehicle_type")
s_na = s %>% 
  filter(is.na(variable_formatted)) %>% 
  select(variable_formatted, variable)
View(s_na)
stats19_schema$variable_formatted[
  is.na(stats19_schema$variable_formatted)
] = stats19_schema$variable[
  is.na(stats19_schema$variable_formatted)
] 
```

```{r}
stats19_schema = stats19_schema %>% 
  filter(variable != "speed_limit")
```



# Update the schemas

```{r}
stats19_variables_old = stats19_variables
stats19_variables = stats19_variables_joined
readr::write_csv(stats19_variables, "data-raw/stats19_variables.csv")

stats19_schema_old = stats19_schema
stats19_schema = stats19_schema_joined
readr::write_csv(stats19_schema, "data-raw/stats19_schema.csv")

usethis::use_data(stats19_variables, overwrite = TRUE)
usethis::use_data(stats19_schema, overwrite = TRUE)
```


To regenerate the file names use the following script:

```{r}
library(rvest)
library(stringr)
u = paste0(
  "https://data.gov.uk/dataset/",
  "cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data")
u = "https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data"
page = read_html(u)
all_links = page %>%
  html_nodes("a") %>%       # find all links
  html_attr("href")

zips = all_links %>% str_subset("\\.zip")
csvs = all_links %>% str_subset("\\.csv")

r = c(zips, csvs)
dr = c()
for(i in 1:length(r)) {
  dr[i] = sub("https://data.dft.gov.uk/road-accidents-safety-data/",
              "", URLdecode(r[i]))
  dr[i] = sub("https://data.dft.gov.uk/road-accidents-safety-data/",
              "", dr[i])
}

file_names = setNames(as.list(dr), dr)
file_name_df = tibble::tibble(
  file_name = unlist(file_names),
  url = r
  )
usethis::use_data(file_names, overwrite = TRUE)

# get original stats19 data filenames
download.file("https://github.com/ropensci/stats19/raw/master/data/file_names.rda", "file_names_old.rda")
load("file_names_old.rda")
file_names_old = file_names
usethis::use_data(file_names_old, overwrite = TRUE)
file_names = setNames(as.list(dr), dr)
# compare the objects
waldo::compare(file_names_old, file_names)
class(file_names_old)
file_names_char = unlist(file_names)
writeLines(file_names_char, "data-raw/file_names.txt")
readr::write_csv(file_name_df, "data-raw/file_name_df.csv")
file.edit("data-raw/file_names.txt")
file.remove("file_names_old.rda")

file_names$`accident-and-casualty-adjustment-2004-to-2019.zip` 
file_names$`accident-and-casualty-adjustment-2004-to-2019.zip` = NULL
file_names$`accident-and-casualty-adjustment-2004-to-2019.zip` 
usethis::use_data(file_names, overwrite = TRUE)
```

The `accidents_sample_raw` can be (re)generated using:

```{r}
# Obtained with:
dl_stats19(year = 2017, type = "Accide")
accidents_2017_raw = read_accidents(year = 2017)
set.seed(350)
sel = sample(nrow(accidents_2017_raw), 3)
accidents_sample_raw = accidents_2017_raw[sel, ]
accidents_sample = format_accidents(accidents_sample_raw)
```

Similarly for casualites, use:

```{r}
# Obtained with:
dl_stats19(year = 2017, type = "cas")
casualties_2017_raw = read_casualties(year = 2017)
set.seed(350)
sel = sample(nrow(casualties_2017_raw), 3)
casualties_sample_raw = casualties_2017_raw[sel, ]
casualties_sample = format_casualties(casualties_sample_raw)
```

and for vehicles, use:
```{r}
# Obtained with:
dl_stats19(year = 2017, type = "veh")
vehicles_2017_raw = read_vehicles(year = 2017)
set.seed(350)
sel = sample(nrow(vehicles_2017_raw), 3)
vehicles_sample_raw = vehicles_2017_raw[sel, ]
vehicles_sample = format_vehicles(vehicles_sample_raw)
```

Aslso, to (re)generate the `police_boundaries` data, use:
```{r}
# Obtained with:
library(sf)
u = "https://opendata.arcgis.com/datasets/3e5a096a8c7c456fb6d3164a3f44b005_3.geojson"
police_boundaries_wgs = sf::st_read(u)
names(police_boundaries_wgs)
police_boundaries = st_transform(police_boundaries_wgs, 27700)
names(police_boundaries)
police_boundaries = police_boundaries[c("pfa16cd", "pfa16nm")]
```

WARNING: this does not work on Windows.

This function generates the data object `stats19_schema` in a reproducible way
using DfT's schema definition (see function [dl_schema()]).

The function also generates `stats19_variables` (see the function's source code
for details).

# Load stats19 schema

```{r}
read_schema = function(
  data_dir = tempdir(),
  filename = "Road-Accident-Safety-Data-Guide.xls",
  sheet = NULL
) {

  file_path = file.path(data_dir, filename)
  if(!file.exists(file_path)) {
    dl_schema()
  }
  if(is.null(sheet)) {
    export_variables = readxl::read_xls(path = file_path,
                                        sheet = 2,
                                        skip = 2)
    export_variables_accidents = tibble::tibble(
      table = "accidents",
      variable = export_variables$`Accident Circumstances`
    )
    export_variables_vehicles =  tibble::tibble(
      table = "vehicles",
      variable = export_variables$Vehicle
    )
    export_variables_casualties =  tibble::tibble(
      table = "casualties",
      variable = export_variables$Casualty
    )
    export_variables_long = rbind(
      export_variables_accidents,
      export_variables_casualties,
      export_variables_vehicles
    )
    stats19_variables = stats::na.omit(export_variables_long)
    stats19_variables$type = stats19_vtype(stats19_variables$variable)

    # add variable linking to names in formatted data
    names_acc = names(accidents_sample)
    names_veh = names(vehicles_sample)
    names_cas = names(casualties_sample)
    names_all = c(names_acc, names_veh, names_cas)

    variables_lower = schema_to_variable(stats19_variables$variable)
    # test result
    # variables_lower[!variables_lower %in% names_all]
    # names_all[!names_all %in% variables_lower]
    stats19_variables$column_name = variables_lower
    # head(stats19_variables)

    #' export result: usethis::use_data(stats19_variables, overwrite = TRUE)

    sel_character = stats19_variables$type == "character"
    character_vars = stats19_variables$variable[sel_character]
    character_vars = stats19_vname_switch(character_vars)

    schema_list = lapply(
      X = seq_along(character_vars),
      FUN = function(i) {
        x = readxl::read_xls(path = file_path, sheet = character_vars[i])
        names(x) = c("code", "label")
        x
      }
    )

    stats19_schema = do.call(what = rbind, args = schema_list)
    n_categories = vapply(schema_list, nrow, FUN.VALUE = integer(1))
    stats19_schema$variable = rep(character_vars, n_categories)

    character_cols = stats19_variables$column_name[sel_character]
    stats19_schema$variable_formatted = rep(character_cols, n_categories)

  } else {
    stats19_schema = readxl::read_xls(path = file_path, sheet = sheet)
  }
  stats19_schema
}
```

# Download schema

```{r}
dl_schema = function(data_dir = tempdir()) {
u = paste0(
  "https://data.dft.gov.uk/road-accidents-safety-data/",
  "Road-Accident-Safety-Data-Guide.xls"
)
destfile = file.path(data_dir, "Road-Accident-Safety-Data-Guide.xls")
utils::download.file(u, destfile = destfile)
# download and unzip the data if it's not present
}
```

# Function to covert schema names

```{r}
schema_to_variable = function(x) {
  x = format_column_names(x)
  x = gsub(pattern = " ", replacement = "_", x = x)
  x = gsub(pattern = "_null_if_not_known", replacement = "", x)
  x = gsub(pattern = "_dd/mm/yyyy|_hh:mm", replacement = "", x)
  x = gsub(pattern = "highway_authority___ons_code", replacement = "highway", x)
  x = gsub(pattern = "lower_super_ouput_area", replacement = "lsoa", x)
  x = gsub(pattern = "_england_&_wales_only", replacement = "", x)
  x = gsub(pattern = "_cc", replacement = "", x)
  x = gsub(pattern = "vehicle_propulsion_code", replacement = "propulsion_code", x)
  x = gsub(pattern = "pedestrian_road_maintenance_worker_from_2011",
           replacement = "pedestrian_road_maintenance_worker", x)
  x = gsub(pattern = "engine_capacity", replacement = "engine_capacity_cc", x)
  x = gsub(pattern = "age_of_vehicle_manufacture", replacement = "age_of_vehicle", x)
  x
}
```

# Return type of variable of stats19 data

```{r}
stats19_vtype = function(x) {
  variable_types = rep("character", length(x))
  sel_numeric = grepl(pattern = "Number|Speed|Age*.of|Capacity", x = x)
  variable_types[sel_numeric] = "numeric"
  sel_date = grepl(pattern = "^Date", x = x)
  variable_types[sel_date] = "date"
  sel_time = grepl(pattern = "^Time", x = x)
  variable_types[sel_time] = "time"
  sel_location = grepl(pattern = "^Location|Longi|Lati", x = x)
  variable_types[sel_location] = "location"
  #' remove other variables with no lookup: no weather ?!
  sel_other = grepl(
    pattern = paste0(
      "Did|Lower|Accident*.Ind|Reference|Restricted|",
      "Leaving|Hit|Age*.Band*.of*.D|Driver*.H"
    ),
    x = x
  )
  variable_types[sel_other] = "other"
  variable_types
}

stats19_vname_switch = function(x) {
  x = gsub(pattern = " Authority - ONS code", "", x = x)
  x = gsub(
    pattern = "Pedestrian Crossing-Human Control",
    "Ped Cross - Human",
    x = x
    )
  x = gsub(
    pattern = "Pedestrian Crossing-Physical Facilities",
    "Ped Cross - Physical",
    x = x
    )
  x = gsub(pattern = "r Conditions", "r", x = x)
  x = gsub(pattern = "e Conditions", "e", x = x)
  x = gsub(pattern = " Area|or ", "", x = x)
  x = gsub(pattern = "Age Band of Casualty", "Age Band", x = x)
  x = gsub(pattern = "Pedestrian", "Ped", x = x)
  x = gsub(pattern = "Bus Coach Passenger", "Bus Passenger", x = x)
  x = gsub(pattern = " \\(From 2011\\)", "", x = x)
  x = gsub(pattern = "Casualty Home Type", "Home Area Type", x = x)
  x = gsub(pattern = "Casualty IMD Decile", "IMD Decile", x = x)
  x = gsub(pattern = "Driver IMD Decile", "IMD Decile", x = x)
  x = gsub(pattern = "Journey Purpose of Driver", "Journey Purpose", x = x)
  x
}
```

Informal test:

```{r}
variable_types = stats19_vtype(stats19_variables$variable)
names(variable_types) = stats19_variables$variable
variable_types
x = names(get_stats19(year = 2017, type = "accidents", ask = FALSE))
n = stats19_vtype(x)
names(n) = x
n
```

The data bundled with the package can then be (re-)generated with

```{r}
stats19_schema = read_schema()
View(stats19_schema)
```
---
title: "SUVs and pedestrian/cycle safety"
subtitle: "`r emojifont::emoji(c('walking', 'bike', 'car', 'truck'))`<br/>empowering planners with data"
author: "Robin Lovelace"
institute: "University of Leeds"
date: "2019-12-05 (updated: `r Sys.Date()`)"
output:
  xaringan::moon_reader:
    css: [default, robot, robot-fonts]
    lib_dir: libs
    nature: 
      beforeInit: "https://platform.twitter.com/widgets.js"
      highlightStyle: github
      highlightLines: true
      countIncrementalSlides: false
---

background-image: url(https://www.jato.com/wp-content/uploads/2019/11/Chart-4-1024x401.jpg)

--

```{r, eval=FALSE, echo=FALSE}
# to reproduce these slides do:
pkgs = c("rgdal", "sf", "geojsonsf")
install.packages(pkgs)
```

--


```{r setup, include=FALSE, echo=FALSE}
options(htmltools.dir.version = FALSE)

```

```{r, eval=FALSE, echo=FALSE}
# get ods data
if(!file.exists("veh0105.ods")) {
  download.file("https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/830778/veh0101.ods", "veh0101.ods")
  download.file("https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/794455/veh0103.ods", "veh0103.ods")
  download.file("https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/830779/veh0104.ods", "veh0104.ods")
  download.file("https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/794433/veh0105.ods", "veh0105.ods")
  download.file("https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/830786/veh0128.ods", "veh0128.ods")
}
# veh0101 = readODS::read.ods("veh0101.ods")
# veh0103 = readODS::read.ods("veh0103.ods")
# veh0104 = readODS::read.ods("veh0104.ods")
# veh0105 = readODS::read.ods("veh0105.ods")
# veh0128 = readODS::read.ods("veh0128.ods")
library(tidyverse)
veh0128 = readr::read_csv("https://github.com/ropensci/stats19/releases/download/1.1.0/veh0128.csv")
veh_names = "QASHQAI|fiesta")
veh0128_sub = veh0128 %>% 
  filter(str_match(string = `Generic model 1`, veh_names))
```


--

## Source [Jato: SUVs make up 40% of total car registrations](https://www.jato.com/suvs-make-up-40-of-total-car-registrations-as-the-european-market-records-its-best-october-result-since-2009/)

--


???

---



background-image: url(https://cdn1.carbuyer.co.uk/sites/carbuyer_d7/files/styles/insert_main_wide_image/public/car_images/peugeot_3008_cutout_best_mid_suv_2019bestbuys.jpg?itok=GMiNNkDy)

--

# Largely unnoticed

In tax band E, D or even A - C for electrics

![](https://user-images.githubusercontent.com/1825120/69948517-b4b54180-14e7-11ea-8dca-ed2172d5789e.png)

---

# Defining SUVs

.pull-left[

No standard [definition](https://en.wikipedia.org/wiki/Sport_utility_vehicle#British_English)

- In USA: "rugged automotive vehicle similar to a station wagon but built on a light-truck chassis"
- In UK: "powerful vehicle with four-wheel drive that can be driven over rough ground. The abbreviation SUV is often used."
]

.pull-right[

Data on weight and height is hard to come by

https://www.epa.gov/compliance-and-fuel-economy-data/data-cars-used-testing-fuel-economy

MOT data excludes weight: https://www.gov.uk/transport/car-motorcycle-and-van-mot-tests

But we do have good data on tax band

]

## Suggestion: define on weight (energy, size) + height (colision height)

---



---

background-image: url(https://informedforlife.org/demos/FCKeditor/UserFiles/Image/weight%20graph%20with%20ave.%20passenger%20vehicle.png)

Source `https://www.bogleheads.org/forum/viewtopic.php?t=281280`

???

Previously my cars have always been midsize sedans. Chevrolets, Mazdas, etc. My intent was and may be still to purchase a more economical commuter. I test drove a Civic, Accord, Elantra, Sonata, Camry, Mazda3, Mazda6, Corolla, and Camry. I'm nothing if not thorough about figuring this out. I've noticed a huge uptick in the average size of vehicles around me during my 2 hours of round trip commuting each day (this is a newer job).

I really noticed this from the drivers seat of a Honda Civic. It seems like only about 10% of the vehicles in my area are actual cars anymore. Pickups, SUVs, Crossovers, Minivans. I never used to feel like safety was a major concern in a smaller vehicle, but it just seems like in the event of a crash I would be much more likely to be on the... squished... end of the accident. Is this unreasonable?
---

Huge cobenefits of downsizing

```{r, echo=FALSE, out.width="70%"}
knitr::include_graphics("https://royalsocietypublishing.org/cms/attachment/6cf8790d-1cc3-4e01-b140-e913718b3f2e/rsta20160364f03.jpg")
```

Serrenho et al (2017) The impact of reducing car weight on global emissions: the future fleet in Great Britain  https://doi.org/10.1098/rsta.2016.0364


???

2008 book by Scott Burkun aimed at software architects

Too techy for me, about managing 100s of people!

---

### A *Sport* Utility Vehicle for the 21st Century

<iframe width="560" height="315" src="https://www.youtube.com/embed/gZ7iGl3j10w?start=10" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

Source: [Quicab](https://www.crowdfunder.co.uk/quicab-e-assisted-cycle-taxis-for-london)

- Slower
- Lighter
- Physical activity
- Better air quality

--

### Saves lives


???

2015 book by Daniel Levitin

Too techy for me, about managing 100s of people!

---

.pull-left[

```{r, echo=FALSE}
knitr::include_graphics("https://ec.europa.eu/easme/sites/easme-site/files/getreal_growinggap_e_year2017.jpg")
```

]

--

.pull-right[

## Real world emissions have diverged

## EU tests show size of gap

See the [Get Real report](https://www.get-real.org/wp-content/uploads/2019/03/Get-Real-CO2-report.pdf)

]

---

## Road crash data

### From open STATS19 data (stats19 R package)

![](https://docs.ropensci.org/stats19/articles/stats19-vehicles_files/figure-html/unnamed-chunk-7-1.png)


---

### Next steps

![](https://raw.githubusercontent.com/Robinlovelace/stats19-gisruk/master/README_files/figure-gfm/unnamed-chunk-25-1.png)
.pull-left[

Crash data vs size (engine + tax data)

]

--

.pull-right[

Estimate potential lives saved

]

---

class: center, middle

# Thanks!

Contact me at r. lovelace at leeds ac dot uk (email), `@robinlovelace` (twitter + github)

--

Check-out my repos at https://github.com/robinlovelace/

--

For more information on stats19 data, see [*stats19 on GitHub*](https://docs.ropensci.org/stats19/)

--

Thanks to all the R developers who made this possible, including (for this presentation):

[remark.js](https://remarkjs.com), [**knitr**](http://yihui.name/knitr), and [R Markdown](https://rmarkdown.rstudio.com).

Slides created via the R package [**xaringan**](https://github.com/yihui/xaringan).

--

Thanks to everyone for building a open and collaborative communities!
---
title: "Introduction to R for Road Safety"
author: "Robin Lovelace"
institute: "Institute for Transport Studies, University of Leeds"
date: "`r Sys.Date()`"
output:
  xaringan::moon_reader:
    lib_dir: libs
    nature:
      highlightStyle: github
      highlightLines: true
      countIncrementalSlides: false
---

background-image: url(https://user-images.githubusercontent.com/1825120/61066243-12dc6d80-a3fd-11e9-8805-826a47c553f6.png)

# About the course team and location

- Robin Lovelace
- Haruko
- Pangiotis
- Martyna 

---

# About the course

- What it is...
  - An opportunity to learn R with support
  - Based on recent, cutting edge software
  - By people experienced teachers/researchers/developers
  
--
  
- What it is not
  - A course on statistical modelling


---

# About the package

See https://docs.ropensci.org/stats19/

![](https://docs.ropensci.org/stats19/reference/figures/README-unnamed-chunk-2-1.png)

---

# R and RStudio demo

Actions speak louder than words...

![](https://raw.githubusercontent.com/ITSLeeds/TDS/master/courses/2day/images/rstudio-ui.png)

---

# Over to you

Work through Section 2 of the exercises

- Projects and scripts
- Writing and running code
- Vewing Objects
- Autocompletion
- Getting help

---

# Time data

Section 3.4 shows how to change classes in R

Time is represented by special classes:

```{r}
library(lubridate)
x = today()
class(x)
day(x)
month(x)
year(x)
weekdays(x)
```

---

# Time representations

```{r}
as.Date("2019-10-17") # works
# as.Date("2019 10 17") # fails
ymd("2019 10 17") # works
dmy("17/10/2019") # works
```

---

# Subsetting time objects

```{r}
c_sample = stats19::accidents_sample
c_sample$date
c_sample$date_formatted = dmy(c_sample$date)
c_sample$date_formatted > ymd("2017-08-01")
```

---
# Next up

- Try running the code and answer the questions in Section 5
- Continue working through Sections 2 and 3

## After lunch

Spatial data (maps!)
---

# Spatial data and maps

```{r, message=FALSE}
library(sf)
crashes = readr::read_csv("https://github.com/ropensci/stats19/releases/download/1.0.0/crashes.csv")
crashes_sf = crashes # create copy of crashes dataset
crashes_sf$longitude = c(-1.3, -1.2, -1.1)
crashes_sf$latitude = c(50.7, 50.7, 50.68)
crashes_sf = st_as_sf(crashes_sf, coords = c("longitude", "latitude"), crs = 4326)
plot(crashes_sf)
```

---

# Doing spatial data with R

- Context: spatial ecosystem (see section [1.4 of Geocomputation with R - package ecosystem](https://geocompr.robinlovelace.net/intro.html#rs-spatial-ecosystem))
- [Exercises](https://geocompr.robinlovelace.net/attr.html#exercises-1): Section 6 of the handout

- Further reading: [Section 3.2 to 3.2.2](https://geocompr.robinlovelace.net/attr.html#vector-attribute-manipulation) of handouts

---

# Practical

- Work through Section 6.1 and try exercises 1 to 3
- Take a short read of [Section 1.4](https://geocompr.robinlovelace.net/intro.html#rs-spatial-ecosystem) of Geocomputation with R (first challenge, find the resource online!)
- Return to the practical

---

# Afternoon session

- Talk on Road Safety 1

- Applying the methods to stats19 data - live demo

    How to access data with stats19
    Key stats19 functions
    Excercises: analysing road crash data on the Isle of Wight

- Continue working through the handouts
- Talk on Road Safety 2
- Homework - ensure that you have at least read-over the handout Sections 1:5

---

# Day 2

- Prioirty: consolidate knowledge from day 1
- Continuing to work through the handouts starting on Section 6
- Demonstrations of counting number of crashes on roads
- Talk on network analysis 
- Talk on R in a professional setting
- Working on your own data

--

- Missed homework! Ensure that you have at least read-over the handout Sections 1:5


---

# The R learning curve

See video: https://www.youtube.com/watch?v=7oyiPBjLAWY&feature=youtu.be&t=357

![](https://i.imgur.com/CBdpkug.jpg)

---

# Agenda

09:30-11:00 Point pattern analysis

- Visualising data with tmap
- Spatial and temporal subsetting
- Aggregation

11:15-12:30 Road network data

- Desire lines: using origin-destination data
- Downloading road network data from OSM
- Buffers on road networks

**Lunch**

13:30-15:00 Analysing crash data on road network

**Break**
  
15:15-15:30: Talk on Road Safety 3

15:30-16:30 Applying the methods to your own data

---

# Bonus extras

![](final-figure.png)

--

Merging, forecasting, network analysis

---

# Visualising spatial data practical


- Foundations of sf: try Sections 6.2 and 6.3 (building on new subsetting knowledge)
- Making maps: try out the exercises in Section 7

- From break until lunch, either:
  - Continue working through the practicals (Sections 2 to 5), or
  - Dive into processing large stats19 dataset (Section 8)

- Sections 1:5 R vital
- Section 3 on subsetting (vital)
- Section 4 on packages (important)
- Section 5 on times (can do later)

--

- Live demo + questions

---

# Further information
---
title: "Introduction to R for road safety: an introduction to R and practical exercises"
subtitle: "Practical exercises based on UK data from the stats19 package, developed by the Institute for Transport Studies, University of Leeds"
# output: pagedown::html_paged
# output: bookdown::html_document2
   # github_document
output:
  pdf_document: 
    highlight: tango
    number_sections: yes
# output: rmarkdown::html_vignette
# vignette: >
#   %\VignetteIndexEntry{R for road safety: exercises}
#   %\VignetteEngine{knitr::rmarkdown}
author: Robin Lovelace, Malcolm Morgan & Andrea Gilardi
bibliography: ../vignettes/references.bib
---

```{r cache, include=FALSE}
knitr::opts_chunk$set(cache = TRUE)
```

```{r upload, eval=FALSE, echo=FALSE}
# Upload exercises
file.copy("inst/stats-19-exercises.pdf", ".", overwrite = TRUE)
piggyback::pb_upload("stats-19-exercises.pdf")
piggyback::pb_download_url("stats-19-exercises.pdf")
# to run this for the first time:
download.file("https://github.com/ITSLeeds/TDS/archive/master.zip", "tds-master.zip")
unzip("tds-master.zip")
file.rename("TDS-master/courses/2day/images/", "inst/images")
```


# Introduction

This document provides information, code and, vitally, exercises to test and improve your R skills.
It starts with introductory R skills that will be of use in any domain but the focus is on R for Road Safety, in support of a [2 day course](https://www.racfoundation.org/introduction-to-r-for-road-safety).
The course is based on open road crash records from the **stats19** package [@lovelace_stats19_2019].
Code and data supporting the content can be found in the package's GitHub repo at [github.com/ropensci/stats19](https://github.com/ropensci/stats19/).
The '[issue tracker](https://github.com/ropensci/stats19/issues)' associated with that repo is a good place to ask questions about the course. 

Course pre-requisites are outlined in the [stats19-training-setup](https://docs.ropensci.org/stats19/articles/stats19-training-setup.html) hosted at [docs.ropensci.org/stats19](https://docs.ropensci.org/stats19).
It makes use of a number of packages which can be installed with `install.packages()` and loaded as follows:

```{r, message=FALSE, warning=FALSE, eval=FALSE}
library(pct)      # access travel data from DfT-funded PCT project 
library(sf)       # spatial vector data classes
library(stats19)  # get stats19 data
library(stplanr)  # transport planning tools
library(tidyverse)# packages for 'data science'
library(tmap)     # interactive maps
```

You should type, run and ensure you understand each line of code in this document.

# R and RStudio

The learning outcomes of this first session are to learn: 
RStudio main features and scripts,
R objects and functions,
subsetting,
basic plotting, and
getting help.

The first exercise is to open up RStudio and take a look around and identify the main components, shown in the figure below.
**Explore each of the main components of RStudio.**
Try changing the Global Settings (in the Tools menu) and see RStudio's short cuts by pressing `Alt-Shift-K` (or `Option+Shift+K` on Mac).

```{r rstudioui, echo=FALSE, out.width="70%"}
knitr::include_graphics("images/rstudio-ui.png")
```

## Projects and scripts

Projects are a way to organise related work together. Each project has its own folder and Rproj file. **Advice: always working from projects will make your life easier!** Start a new project with:

> File > New Project
You can choose to create a new directory (folder) or associate a project with an existing directory. Make a new project called stats1-course and save it in a sensible place on your computer. Notice that stats1-course now appears in the top right of RStudio.

Scripts are the files where R code is stored.
**Keeping your code in sensibly named, well organised and reproducible scripts will make your life easier:**
you could simply type all our code into the console, but that require retyping commands each time you run it.
Instead, code that you want to keep and share should be saved script files, plain text files that have the `.R` extension.

Make a new script with Flie > New File > Rscript or Ctrl+Shift+N

Save the script and give it a sensible name like `stats19-lesson-1.R` with File > Save, the save button on the toolbar, or Ctrl+S.

**Pro tip:** You can also create new R scripts by typing and running this command in the R console:

```{r edit, eval=FALSE}
file.edit("stats19-lesson-1.R")
```

Keeping scripts and other files associated with a project in a single folder per project (in an RStudio project) will help you find things you need and develop an efficient workflow.

## Writing and running code

Let's start with some basic R operations.
Write this code into your new `stats19-lesson-1.R` R script and execute the result line-by-line by pressing Ctrl+Enter

```{r, eval=FALSE}
x = 1:5
y = c(0, 1, 3, 9, 18)
plot(x, y)
```

This code creates two objects, both are vectors of 5 elements, and then plots them (bonus: check their length using the `length()` function).
Save the script by pressing Ctrl+S.

There are several ways to run code within a script and it is worth becoming familiar with each.
Try running the the code you saved in the previous section using each of these methods:

1. Place the cursor in different places on each line of code and press `Ctrl+Enter` to run that line of code.
1. Highlight a block of code or part of a line of code and press `Ctrl+Enter` to run the highlighted code.
1. Press `Ctrl+Shift+Enter` to run all the code in a script.
1. Press the Run button on the toolbar to run all the code in a script.
1. Use the function `source()` to run all the code in a script e.g. `source("stats19-lesson-1.R")`
<!-- (but don't create an infinite loop!) -->

**Pro tip:** Try jumping between the console and the source editor by pressing Ctl+1 and Ctl+2.

## Viewing Objects

Create new objects by typing and running the following code chunk in a new script, e.g. called `objects.R`.

```{r}
vehicle_type = c("car", "bus", "tank")
casualty_type = c("pedestrian", "cyclist", "cat")
casualty_age = seq(from = 20, to = 60, by = 20)
set.seed(1)
dark = sample(x = c(TRUE, FALSE), size = 3, replace = TRUE)
small_matrix = matrix(1:24, nrow = 12)
crashes = data.frame(vehicle_type, casualty_type, casualty_age, dark)
```

We can view the objects in a range of ways:

1. Type the name of the object into the console e.g. `crashes` and `small_matrix`. Scroll up to see the numbers that didn't fit on the screen.
1. Use the `head()` function to view just the first 6 rows e.g. `head(small_matrix)`
1. Bonus: use the `n` argument in the previous function call to show only the first 2 rows of `small_matrix`
1. Click on the `crashes` object in the environment tab to View it in a spreadsheet.
1. Run the command `View(vehicle_type)`. What just happened?

We can also get an overview of an object using a range of functions, including 
`summary()`,
`class()`,
`typeof()`,
`dim()`, and
`length()`.

You can, for example, view a summary of the `casualty_age` variable by running the following line of code:

```{r summary}
summary(casualty_age)
```

**Exercise** try these functions on each of the objects, what results do they give?

```{r summary-answers, echo=FALSE, eval=FALSE}
summary(vehicle_type)
class(vehicle_type)
typeof(vehicle_type)
dim(vehicle_type)
length(vehicle_type)
```

**Bonus**: Find out the class of the column `vehicle_type` in the data frame `crashes` with the command `class(crashes$vehicle_type)`.
Why has it changed? 
Create a new object called `crashes_char` that keeps the class of the character vectors intact by using the function `tibble::tibble()` (see [tibble.tidyverse.org](https://tibble.tidyverse.org/) and Section 4 for details).

```{r tibble1, echo=FALSE, eval=FALSE}
tibble::tibble(
  vehicle_type,
  casualty_type,
  casualty_age,
  dark
)
```

## Autocompletion

RStudio can help you write code by autocompleting it. RStudio will look for similar objects and functions after typing the first three letters of a name.

```{r autocomp, echo=FALSE}
knitr::include_graphics("images/autocomplete.jpg")
```

When there is more than one option you can select from the list using the mouse or arrow keys.
Within a function, you can get a list of arguments by pressing Tab.

```{r help, echo=FALSE}
knitr::include_graphics("images/fucntionhelp.jpg")
```

## Getting help

Every function in R has a help page. You can view the help using `?` for example `?sum`. Many packages also contain vignettes, these are long form help documents containing examples and guides. `vignette()` will show a list of all the vignettes available, or you can show a specific vignette for example `vignette(topic = "sf1", package = "sf")`.

## Commenting Code

It is good practice to use comments in your code to explain what it does. You can comment code using `#`

For example:

```{r}
# Create vector objects (a whole line comment)
x = 1:5 # a seqence of consecutive integers (inline comment)
y = c(0, 1, 3, 9, 18.1) 
```

You can comment/uncomment a whole block of text by selecting it and using `Ctrl+Shift+C`.
<!-- not sure about the next statement so commenting out (RL) -->
<!-- and you can reformat a block of code using `Ctrl+Shift  + /`.  -->

**Pro tip:** You can add a comment section using Ctrl + Shift + R


## The global environment

The Environment tab shows all the objects in your environment, this includes datasets, parameters, and any functions you have created.
By default, new objects appear in the Global Environment but you can see other environments with the drop-down menu.
For example, each package has its own environment.

Sometimes you wish to remove things from your environment, perhaps because you no longer need them or things are getting cluttered.

You can remove an object with the `rm()` function e.g. `rm(x)` or `rm(x, y)` or you can clear your whole environment with the broom button on the Environment Tab.

1. Remove the object `x` that was created in a previous section.
1. What happens when you try to print the `x` by entering it into the console?
1. Try running the following commands in order: `save.image(); rm(list = ls()); load(".RData")`. What happened?
1. How big (how many bytes) is the `.RData` file in your project's folder?
1. Tidy up by removing the `.Rdata` file with `file.remove(".Rdata")`.

## Debugging Code

All the code shown so far is reproducible.
To test RStudio's debugging features, let's write some code that fails, as illustrated in the figure below.

```{r debug, echo=FALSE, out.width="60%"}
knitr::include_graphics("rstudio-autocomplete.png")
```

1. What is the problem with the code shown in the figure?
1. Create other types of error in the code you have run (e.g. no symetrical brackets and other typos)
1. Does RStudio pick up on the errors? And what happens when you try to run buggy code?

**Always address debugging prompts to ensure your code is reproducible**

## Saving R objects

We have already seen that you can save R scripts.
You can also save individual R objects in the RDS format.

```{r}
saveRDS(crashes, "crashes.Rds")
```

We can also read back in our data.

```{r}
crashes2 = readRDS("crashes.Rds")
identical(crashes, crashes2)
```

R also supports many other formats, including CSV files, which can be created and imported with the functions `readr::read_csv()` and `readr::write_csv()` (see also the [readr](https://readr.tidyverse.org/) package).

```{r readr-write, eval=FALSE}
readr::write_csv(crashes, "crashes.csv")
crashes3 = readr::read_csv("crashes.csv")
identical(crashes3, crashes) 
```

Notice that `crashes3` and `crashes` are not identical, what has changed? Hint: read the help page associated with `?readr::write_csv`.

# Manipulating R objects

## Subsetting by index or name

Subsetting returns part of an R object. 
It can be done by providing numbers representing the positions of the elements we want (e.g. the 2^nd^ element) or with a logical vector, with values associated with `TRUE` returned. 
Two dimension object such as matrices and data frames can be subset by rows and columns.
Subsetting in base R is done with square brackets `[]` after the name of an object. **Run the following commands to practice subsetting.**

```{r, eval=FALSE}
casualty_age[2:3] # second and third casualty_age
crashes[c(1, 2), ] # first and second row of crashes
crashes$vehicle_type # returns just one column
crashes[, c("casualty_type", "casualty_age")] # first and third columns
```

```{r, eval=FALSE, echo=FALSE}
crashes[, c(1, 3)] # first and third column of crashes by positional numbers
crashes[c(2), c(3)]
crashes[c(2), c(2, 3)]
class(crashes[, c(1, 3)])
class(crashes[c(2), c(3)])
```

1. Use the `$` operator to print the `dark` column of `crashes`.
1. Subset the crashes with the `[,]` syntax so that only the first and third columns of `crashes` are returned.
1. Return the 2^nd^ row and the 3^rd^ column of the `crashes` dataset. 
1. Return the 2^nd^ row and the columns 2:3 of the `crashes` dataset. 
1. **Bonus**: what class resulted from each of the previous exercises?

## Subsetting by values

It is also possible to subset objects by the values of their elements.
This works because the `[` operator accepts logical vectors returned by queries such as 'is it less than 3?' (`x < 3` in R) and 'was it light?' (`crashes$dark == FALSE`), as demonstrated below:

```{r, eval=FALSE}
x[c(TRUE, FALSE, TRUE, FALSE, TRUE)] # 1st, 3rd, and 5th element in x
x[x == 5] # only when x == 5 (notice the use of double equals)
x[x < 3] # less than 3
x[x < 3] = 0 # assign specific elements
casualty_age[casualty_age %% 6 == 0] # just the ages that are a multiple of 6
crashes[crashes$dark == FALSE, ]
```

1. Subset the `casualty_age` object using the inequality (`<`) so that only elements less than 50 are returned.
1. Subset the `crashes` data frame so that only tanks are returned using the `==` operator.
1. **Bonus**: assign the age of all all tanks to 61.

```{r, eval=FALSE, echo=FALSE}
casualty_age[casualty_age < 50] # the  casualty_age less than 50
crashes[crashes$vehicle_type == "tank", ] # rows where the name is tank
crashes$casualty_age[crashes$vehicle_type == "tank"] = 61
```

## Dealing with NAs and recoding

R objects can have a value of NA. This is how R represents missing data.

```{r, eval=FALSE}
z = c(4, 5, NA, 7)
```

NA values are common in real-world data but can cause trouble, for example

```{r, eval=FALSE}
sum(z) # result is NA
```

Some functions can be told to ignore NA values.

```{r, eval=FALSE}
sum(z, na.rm = TRUE) # result is equal to 4 + 5 + 7
```

You can find NAs using the `is.na()` function, and then remove them

```{r, eval=FALSE}
is.na(z)
z_nona = z[!is.na(z)] # note the use of the not operator !
sum(z)
```

If you remove records with NAs be warned: the average of a value excluding NAs may not be representative.

## Changing class

Sometimes you may want to change the class of an object.
This is called class coercion, and can be done with functions such as `as.logical()`, `as.numeric()` and `as.matrix()`.

1. Coerce the `vehicle_type` column of `crashes` to the class `character`.
1. Coerce the `crashes` object into a matrix. What happened to the values?
1. **Bonus:** What is the difference between the output of `summary()` on `character` and `factor` variables?

```{r, echo=FALSE, eval=FALSE}
crashes$vehicle_type = as.character(crashes$vehicle_type)
as.matrix(crashes)
```

## Recoding values

Often it is useful to 'recode' values.
In the raw STATS19 files, for example, -1 means NA.
There are many ways to recode values in R, the simplest and most mature of which is the use of factors, as shown below:

```{r}
z = c(1, 2, -1, 1, 3)
l = c(NA, "a", "b", "c") # labels in ascending order
z_factor = factor(z, labels = l)
z_charcter = as.character(z_factor)
z_charcter
```

1. Recode `z` to Slight, Serious and Fatal for 1:3 respectively.
1. Bonus: read the help file at `?dplyr::case_when` and try to recode the values using this function.

## Now you are ready to use R

**Bonus: reproduce the following plot**

```{r smile, out.width="30%", fig.align="center"}
eyes = c(2.3, 4, 3.7, 4)
eyes = matrix(eyes, ncol = 2, byrow = T)
mouth = c(2, 2, 2.5, 1.3, 3, 1, 3.5, 1.3, 4, 2)
mouth = matrix(mouth, ncol = 2, byrow = T)
plot(eyes, type = "p", main = "RRR!", cex = 2, xlim = c(1, 5), ylim = c(0, 5))
lines(mouth, type = "l", col = "red")
```

\newpage

# R Packages

## What are packages?

R has over 15,000 packages (effectively plugins for base R), extending it in almost every direction of statistics and computing.
Packages provide additional functions, data and documentation. They are very often written by subject-matter experts and therefore tend to fit well with the workflow of the analyst in that particular specialism.
There are two main stages to using a package: installing it and loading it.
A third stage is updating it, this is also important.

<!-- installing it... -->
Install new packages from [The Comprehensive R Archive Network](https://cran.r-project.org/) with the command `install.packages()` (or `remotes::install_github()` to install from GitHub).
Update packages with the command `update.package()` or in Tools > Check for Package Updates in RStudio.
You only need to install a package once.
<!-- **Note: avoid `install.packages()` within a script** -->
<!-- Packages only need to be installed once. -->
<!-- You can use `remotes::install_cran()` or `remotes::install_github()` to only install a package if it is not yet installed and up-to-date (note: you only need to use one of these): -->

```{r, eval=FALSE}
install.packages("sf")
# remotes::install_github("r-spatial/sf")
```

<!-- now talk about loading packages -->
Installed packages are loaded with the command `library()`.
Usually, the package will load silently.
In some cases the package will provide a message, as illustrated below.

```{r}
library(sf)
```

To use a function in a package without first loading the package, use double colons, as shown below (this calls the `tibble()` function from the `tibble` package).

```{r tibble2, eval=FALSE}
crashes_tibble = tibble::tibble(
  vehicle_type,
  casualty_type,
  casualty_age,
  dark
)
```

1. Take a look in the Packages tab in the Files pane in RStudio (bottom right by default).
1. What version of the `stats19` package is installed on your computer?
1. Run the command `update.packages()`. What happens? Why?

## ggplot2

Let's take a look at a particular package.
`ggplot2` is a generic plotting package that is part of the ['tidyverse'](https://www.tidyverse.org/) meta-package, which is an "opinionated collection of R packages designed for data science". 
All packages in the tidyverse "share an underlying design philosophy, grammar, and data structures". 
`ggplot2` is flexible, popular, and has dozens of add-on packages which build on it, such as `gganimate`.
To plot non-spatial data, it works as follows (see figure below, left for result):

```{r, message=FALSE, out.width="40%", eval=FALSE}
library(ggplot2)
ggplot(crashes) + geom_point(aes(x = casualty_type, y = casualty_age))
```

Note that the `+` operator adds layers onto one another.

1. Install a package that build on `ggplot2` that begins with with `gg`. Hint: enter `install.packages(gg)` and hit Tab when your cursor is between the `g` and the `)`.
1. Open a help page in the newly installed package with the `?package_name::function()` syntax.
1. Attach the package.
1. **Bonus:** try using functionality from the new 'gg' package building on the example above (hint: the right plot below uses the economist theme from the `ggthemes` package, try other themes).

```{r gg-extend, echo=FALSE, message=FALSE, eval=FALSE}
library(ggplot2)
g1 = ggplot(crashes) + geom_point(aes(x = casualty_type, y = casualty_age)) ] th
# install.packages("ggthemes")
g2 = ggplot(crashes) + geom_point(aes(x = casualty_type, y = casualty_age)) +
  ggthemes::theme_economist()
g3 = cowplot::plot_grid(g1, g2)
ggsave(filename = "inst/ggtheme-plot.png", width = 8, height = 2, dpi = 80)
```

```{r gg2, echo=FALSE, out.width="80%", fig.align="center"}
library(ggplot2)
knitr::include_graphics("ggtheme-plot.png")
```

## dplyr and pipes

Another useful package in the tidyverse is `dplyr`.
It provides functions for manipulating data frames and using the pipe operator ` %>% `. 
The pipe puts the output of one command into the first argument of the next, as shown below (note the results are the same):

```{r}
library(dplyr)
class(crashes)       
crashes %>% class()
```

Useful `dplyr` functions are demonstrated below.

```{r, eval=FALSE}
crashes %>%
  filter(casualty_age > 50) # filter rows
crashes %>%
  select(casualty_type) # select just one column
crashes %>%
  group_by(dark) %>% 
  summarise(mean_age = mean(casualty_age))
```

1. Use `dplyr` to filter row in which `casualty_age` is less than 18, and then 28.
1. Use the `arrange` function to sort the `crashes` object in descending order of age (hint: see the `?arrange` help page).
1. Read the help page of `dplyr::mutate()`. What does the function do?
1. Use the mutate function to create a new variable, `birth_year`, the current year minus their age.
1. **Bonus:** Use the ` %>% ` operator to filter the output from the previous exercise so that only observations with `birth_year` after 1969 are returned.

```{r dplyr, eval=FALSE, echo=FALSE}
# answers
crashes %>% 
  arrange(desc(casualty_age))
crashes %>% filter(casualty_age > 21)
crashes %>% 
  mutate(birth_year = 2019 - casualty_age) %>% 
  filter(birth_year > 1969)
```

# Temporal data

For the analysis and manipulation of temporal data we will first load the R package `lubridate`:

```{r, message=FALSE}
library(lubridate)
```

The simplest example of a Date object that we can analyze is just the current date, i.e.

```{r}
today()
```

We can manipulate this object using several `lubridate` functions to extract the current day, month, year, weekday and so on...

```{r, eval=FALSE}
x = today()
day(x)
month(x)
year(x)
weekdays(x)
```

Exercises: 

1. Look at the help page of the function `month` to see how it is possible to extract the current month as character vector 
1. Look at other functions in lubridate to extract the current weekday as a number, the week of year and the day of the year

Date variables are often stored simply as a character vectors.
This is a problem, since R is not always smart enough to distinguish between character vectors representing Dates.
`lubridate` provides functions that can translate a wide range of date encodings such as `ymd()`, which extracts the Year Month and Day from a character string, as demonstrated below.

```{r, eval=FALSE}
as.Date("2019-10-17") # works
as.Date("2019 10 17") # fails
ymd("2019 10 17") # works
dmy("17/10/2019") # works
```

Import function such as `read_csv` try to recognize the Date variables.
Sometimes this fails.
You can manually create Date objects, as shown below.

```{r}
x = c("2009-01-01", "2009-02-02", "2009-03-03")
x_date = ymd(x)
x_date
```

Exercises: 

1. Extract the day, the year-day, the month and the weekday (as a non-abbreviated character vector) of each element of `x_date`. 
1. Convert `"09/09/93"` into a date object and extract its weekday. 
1. **Bonus:** Read the help page of `as.Date` and `strptime` for further details on the format argument in base R. 

```{r, echo=FALSE, eval=FALSE}
# 1. Extract the day, the year-day, the month and the weekday (as a non-abbreviated character vector) of each element of `x_date`. 
day(x_date)
yday(x_date)
month(x_date)
weekdays(x_date, abbreviate = FALSE)
# 1. Modify the previous example to parse the following character string: `"09/09/1993"` and extract its weekday. 
weekdays(dmy("09/09/93"))
wday(dmy("09/09/93"))
```

We can use Dates also for subsetting events in a dataframe. For example, if we define `x_date` as before and add it to the `crash` dataset, i.e.

```{r}
crashes$casualty_day = x_date
```

then we can subset events using Dates. For example

```{r}
filter(crashes, day(casualty_day) < 7) # the events that ocurred in the first week of the month
filter(crashes, weekdays(casualty_day) == "Monday") # the events occurred on monday
```

Exercises: 

1. Select only the events (rows in `crashes`) that occurred in January
1. Select only the events that ocurred in an odd year-day 
1. Select only the events that ocurred in a leap-year (HINT: check the function `leap_year`)
1. Select only the events that ocurred during the weekend or in June
1. Select only the events that ocurred during the weekend and in June
1. Count how many events that ocurred during each day of the week. 

Now we'll take a look at the time components of a Date. Using the function `hms` (acronym for Hour Minutes Seconds) and its subfunctions such as `hm` or `ms`, we can parse a character vector representing several times as an Hour object (which is tecnically called a Period object). 

```{r}
x = c("18:23:35", "00:00:01", "12:34:56")
x_hour = hms(x)
x_hour
```

We can manipulate these objects using several `lubridate` functions to extract the hour component, the minutes and so on:

```{r}
hour(x_hour)
minute(x_hour)
second(x_hour)
```

If the Hour data do not specify the seconds, then we just have to use a subfunction of `hms`, namely `hm`, and everything works as before. 

```{r}
x = c("18:23", "00:00", "12:34")
(x_hour = hm(x))
```

We can use Hour data also for subsetting events, like we did for Dates. Let's add a new column to crashes data, 

```{r}
crashes$casualty_hms = hms(c("18:23:35", "00:00:01", "12:34:56"))
crashes$casualty_hour = hour(crashes$casualty_hms)
```

Exercises: 

1. Filter only the events that ocurred after midday (i.e. the PM events). Hint: your answer may include `>= 12`.
1. Filter only the events that ocurred between 15:00 and 19:00
<!-- 1. Round all hours to the next hour. Hint: Look at the help page of the `round_date` function.  -->
1. **Bonus:** (difficult): run the following code, which downloades data for car crashes occurred during 2017.

```{r, eval=FALSE}
library(stats19)
crashes_2017 = stats19::get_stats19(year = 2017, type = "ac")
crashes_2017
```

Extract the weekday from the variable called `date`.
How many crashes happened on Monday?

**Advanced challenge:** calculate how many crashes occurred for each day of the week. Then plot it with ggplot2. Repeat the same exercises extracting the hour of the car accident from the variable called time. How would you combine the two informations in a single plot? 

```{r, eval=FALSE, echo=FALSE}
# solutions
crashes %>% filter(casualty_hour >= 12)
crashes %>% filter(casualty_hour > 15 & casualty_hour < 19)

crashes_2017 %>% 
  mutate(my_weekdays = weekdays(date)) %>%
  filter(my_weekdays == "Monday") %>% 
  nrow()
crashes_2017 %>% 
  mutate(my_weekdays = weekdays(date)) %>%
  filter(my_weekdays == "Friday") %>% 
  nrow()

crashes_2017 %>% 
  mutate(my_weekdays = weekdays(date)) %>% 
  group_by(my_weekdays) %>% 
  summarize(n = n()) %>% 
  ggplot() + 
  geom_col(aes(x = my_weekdays, y = n))

crashes_2017 %>% 
  mutate(my_hours = hour(hm(time))) %>% 
  group_by(my_hours) %>% 
  summarize(n = n()) %>% 
  ggplot() + 
  geom_col(aes(x = my_hours, y = n))

crashes_2017 %>% 
  mutate(my_weekdays = weekdays(date), my_hours = hour(hm(time))) %>% 
  group_by(my_weekdays, my_hours) %>% 
  summarise(n = n()) %>% 
  ggplot() + 
  geom_line(aes(x = my_hours, y = n, col = my_weekdays), size = 1.05)
# the legend needs some reordering
```

# Spatial data

## sf objects

All road crashes happen somewhere and, in the UK at least, all collisions recorded by the police are given geographic coordinates, something that can help prioritise interventions to save lives by intervening in and around 'crash hotspots'.
R has strong geographic data capabilities, with the `sf` package provides a generic class for spatial vector data: points, lines and polygons, are represented in `sf` objects as a special 'geometry column', typically called 'geom' or 'geometry', extending the data frame class we've already seen in `crashes`.

Create an `sf` data frame called `crashes_sf` as follows:

```{r crashes-sf, fig.height=2, fig.width=3}
library(sf) # load the sf package for working with spatial data
crashes_sf = crashes # create copy of crashes dataset
crashes_sf$longitude = c(-1.3, -1.2, -1.1)
crashes_sf$latitude = c(50.7, 50.7, 50.68)
crashes_sf = st_as_sf(crashes_sf, coords = c("longitude", "latitude"), crs = 4326)
# plot(crashes_sf[1:4]) # basic plot
# mapview::mapview(crashes_sf) # for interactive map
```

1. Plot only the geometry of `crashes_sf` (hint: the solution may contain `$geometry`). If the result is like the figure below, congratulations, it worked!).
1. Plot `crashes_sf`, only showing the age variable.
1. Plot the 2^nd^ and 3^rd^ crashes, showing which happened in the dark.
1. **Bonus**: How far are the points apart (hint: `sf` functions begin with `st_`)?
1. **Bonus**: Near which settlement did the tank runover the cat?

```{r crashes-sf-ex, echo=FALSE, out.width="30%", fig.show='hold'}
plot(crashes_sf$geometry)
plot(crashes_sf["casualty_age"])
plot(crashes_sf[2:3, "dark"])
# st_distance(crashes_sf)
# Bembridge

# # updload geographic crash data
# write_sf(crashes_sf, "crashes_sf.geojson")
# piggyback::pb_upload("crashes_sf.geojson")
```

## Reading and writing spatial data

You can read and write spatial data with `read_sf()` and `write_sf()`, as shown below (see `?read_sf`).

```{r, eval=FALSE}
write_sf(zones, "zones.geojson") # save geojson file
write_sf(zones, "zmapinfo", driver = "MapInfo file")
read_sf("zmapinfo") # read in mapinfo file
```

See [Chapter 6](https://geocompr.robinlovelace.net/read-write.html) of Geocomputation with R for further information.

## sf polygons

`sf` objects can also represent administrative zones.
This is illustrated below with reference to `zones`, a spatial object representing the Isle of Wight, that we will download using the `pct` package (note: the `[1:9]` appended to the function selects only the first 9 columns).

```{r}
zones = pct::get_pct_zones("isle-of-wight")[1:9]
```

1. What is the class of the `zones` object?
1. What are its column names?
1. Print its first 2 rows and columns 6:8 (the result is below).

```{r, echo=FALSE}
# class(zones)
# names(zones)
zones[1:2, c(1, 5, 6, 7, 8)]
```

## Spatial subsetting and sf plotting

Like index and value subsetting, spatial subsetting can be done with the `[` notation.
Subset the `zones` that contain features in `crashes_sf` as follows:

```{r, message=FALSE}
zones_containing_crashes = zones[crashes_sf, ]
```

To plot a new layer on top of an existing `sf` plot, use the `add = TRUE` argument.
Remember to plot only the `geometry` column of objects to avoid multiple maps.
Colours can be set with the `col` argument.

1. Plot the geometry of the zones, with the zones containing crashes overlaid on top in red.
1. Plot the zone containing the 2^nd^ crash in blue.
1. **Bonus:** plot all zones that intersect with a zone containing crashes, with the actual crash points plotted in black.

```{r sp-ex, echo=FALSE, out.width="33%", fig.show='hold', message=FALSE, warning=FALSE}
plot(zones$geometry)
plot(zones_containing_crashes$geometry, col = "red", add = TRUE)
plot(zones$geometry)
plot(zones[crashes_sf[2, ], ], col = "blue", add = TRUE)
plot(zones$geometry)
plot(zones[zones_containing_crashes, ], col = "yellow", add = TRUE)
plot(crashes_sf$geometry, pch = 20, add = TRUE)
```

## Geographic joins

Geographic joins involve assigning values from one object to a new column in another, based on the geographic relationship between them.
With `sf` objects it works as follows:

```{r, message=FALSE}
zones_joined = st_join(zones[NULL], crashes_sf)
```

1. Plot the `casualty_age` variable of the new `zones_joined` object (see the figure below to verify the result).
1. How many zones are returned in the previous command? 
1. Select only the `geo_code` column from the `zones` and the `dark` column from `crashes_sf` and use the `left = FALSE` argument to return only zones in which crashes occured. Plot the result.

See [Chapter 4](https://geocompr.robinlovelace.net/spatial-operations.html#spatial-joining) of Geocomputation with R [@lovelace_geocomputation_2019] for further information on geographic joins.

```{r joinf, echo=FALSE, out.width="40%", fig.show='hold', message=FALSE}
plot(zones_joined["casualty_age"])
zjd = st_join(zones[1], crashes_sf["dark"], left = FALSE)
plot(zjd)
```


## CRSs

Get and set Coordinate Reference Systems (CRSs) with the command `st_crs()`.
Transform CRSs with the command `st_transform()`, as demonstrated in the code chunk below, which converts the 'lon/lat' geographic CRS of `crashes_sf` into the projected CRS of the British National Grid:

```{r crs1}
crashes_osgb = st_transform(crashes_sf, 27700)
```

1. Try to subset the zones with the `crashes_osgb`. What does the error message say?
1. Create `zones_osgb` by transforming the `zones` object.
1. **Bonus:** use `st_crs()` to find out the units measurement of the British National Grid?

For more information on CRSs see [Chapter 6](https://geocompr.robinlovelace.net/reproj-geo-data.html) of Geocompuation with R.

## Buffers

Buffers are polygons surrounding geometries of a (usually) fixed distance.
Currently buffer operations in R only work on objects with projected CRSs.

1. Find out and read the help page of `sf`'s buffer function.
1. Create an object called `crashes_1km_buffer` representing the area within 1 km of the crashes.
1. **Bonus:** try creating buffers on the geographic version of the `crashes_sf` object. What happens?

## Attribute operations on sf objects

Because `sf` objects are `data.frame`s, we can do non-spatial operations on them.
Try the following attribute operations on the `zones` data.

```{r, eval=TRUE}
# load example dataset if it doesn't already exist
zones = pct::get_pct_zones("isle-of-wight")
sel = zones$all > 3000  # create a subsetting object
zones_large = zones[sel, ] # subset areas with a popualtion over 100,000
zones_2 = zones[zones$geo_name == "Isle of Wight 002",] # subset based on 'equality' query
zones_first_and_third_column = zones[c(1, 3)]
zones_just_all = zones["all"]
```


1. Practice subsetting techniques you have learned on the `sf data.frame` object `zones`:
     1. Create an object called `zones_small` which contains only regions with less than 3000 people in the `all` column
     1. Create a selection object called `sel_high_car` which is `TRUE` for regions with above median numbers of people who travel by car and `FALSE` otherwise
     1. Create an object called `zones_foot` which contains only the foot attribute from `zones`
     1. Bonus: plot `zones_foot` using the function `plot` to show where walking is a popular mode of travel to work
     1. Bonus: bulding on your answers to previous questions, use `filter()` from the `dplyr` package to subset small regions where car use is high. 
1. Bonus: What is the population density of each region (hint: you may need to use the functions `st_area()`, `as.numeric()` and use the 'all' column)?
1. Bonus: Which zone has the highest percentage of people who cycle?

```{r, echo=FALSE, eval=FALSE}
# 1. Practice subsetting techniques you have learned on the `sf data.frame` object `zones`:
#      1. Create an object called `zones_small` which contains only regions with less than 3000 people in the `all` column
# in base R
zones_small = zones[zones$all < 3000, ]
# with dplyr
zones_small = zones %>% 
  filter(all < 3000)
#      1. Create a selection object called `sel_high_car` which is `TRUE` for regions with above median numbers of people who travel by car and `FALSE` otherwise
median_car = median(zones$car_driver)
sel_high_car = zones$car_driver > median_car 
#      1. How many regions have the number '1' in the column 'geo_name'? What percentage of the regions in the Isle of Wight is this?
sel_region_name_contains_1 = grepl("1", x = zones$geo_name)
sum(sel_region_name_contains_1) / nrow(zones)
#      1. Create an object called `zones_foot` which contains only the foot attribute from `zones`
# using base R
zones_foot = zones["foot"]
# dplyr
zones_foot = zones %>% 
  select(foot)
#      1. Bonus: plot the result to show where walking is a popular mode of travel to work
plot(zones_foot)
#      1. Bonus: bulding on your answers to previous questions, use `filter()` from the `dplyr` package to subset small regions where high car use is high
zones_small_car_high = zones %>% 
  filter(all < 3000, car_driver > median_car)
# 1. Bonus: What is the population density of each region (hint: you may need to use the functions `st_area()`, `as.numeric()` and use the 'all' column)?
zones$area_km2 = as.numeric(st_area(zones)) /1000000
zones$population_density = zones$all / zones$area_km2
plot(zones["population_density"])
# in dplyr
zones_density = zones %>% 
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1000000) %>% 
  mutate(population_density = all / area_km2)
plot(zones_density %>% select(population_density))
# 1. Bonus: Which zone has the highest percentage who cycle?
zones %>% 
  mutate(pcycle = bicycle / all) %>% 
  top_n(n = 1, wt = pcycle)
# 1. Bonus: Find the proportion of people who drive to work (`car_driver`) in areas in which more than 500 people walk to work
zones %>% 
  group_by(foot > 500) %>% 
  summarise(mean_car = sum(car_driver) / sum(all) )
```

## Matching roads to crashes



# Visualising spatial datasets

So far we have used the `plot()` function to make maps.
That's fine for basic visualisation, but for publication-quality maps, we recommend using `tmap` (see Chapter 8 of Geocomputation with R for reasons and alternatives).
Load the package as follows:

```{r}
library(tmap)
tmap_mode("plot")
```

1. Create the following plots using `plot()` and `tm_shape() + tm_polygons()` functions (note: the third figure relies on setting `tmap_mode("view")`.
1. Add an additional layer to the interactive map showing the location of crashes, using marker and dot symbols.
1. Bonus: Change the default basemap (hint: you may need to search in the package documentation or online for the solution).

```{r plot3, fig.show='hold', out.width="33%", echo=FALSE}
plot(zones[c("all", "bicycle")])
tm_shape(zones) + 
  tm_polygons(c("all", "bicycle"))
tmap_mode("view")
m = tm_shape(zones_joined) + 
  tm_polygons(c("casualty_type")) +
  tm_scale_bar()
knitr::include_graphics("tmap-zones-interactive.png")
```

# Analysing point data from stats19

Based on the saying "don't run before you can walk", we've learned the vital foundations of R before tackling a real dataset.
Temporal and spatial attributes are key to road crash data, hence the emphasis on `lubridate` and `sf`.
Visualisation is key to understanding and policy influence, which is where `tmap` comes in.
With these solid foundations, plus knowledge of how to ask for help (ask R's internal help, colleagues, online forums/GitHub, generally in that order of priority), you are ready to test the methods on some real data.

Before doing so, take a read of the `stats19` vignette, which can be launched as follows:

```{r, eval=FALSE}
vignette(package = "stats19") # view all vignettes available on stats19
vignette("stats19") # view the introductory vignette
```

This should now be sufficient to tackle the following exercises:

1. Download and plot all crashes reported in Great Britain in 2018 (hint: see [the stats19 vignette](https://cran.r-project.org/web/packages/stats19/vignettes/stats19.html))
1. Find the function in the `stats19` package that converts a `data.frame` object into an `sf` data frame. Use this function to convert the road crashes into an `sf` object, called `crashes_sf`, for example.
1. Filter crashes that happened in the Isle of Wight based on attribute data (hint: the relevant column contains the word `local`)
1. Filter crashes happened in the Isle of Wight using geographic subsetting (hint: remember `st_crs()`?)
1. Bonus: Which type of subsetting yielded more results and why? 
1. Bonus: how many crashes happened in each zone?
1. Create a new column called `month` in the crash data using the function `lubridate::month()` and the `date` column.
1. Create an object called `a_zones_may` representing all the crashes that happened in the Isle of Wight in the month of May
1. Bonus: Calculate the average (`mean`) speed limit associated with each crash that happened in May across the zones of the Isle of Wight (the result is shown in the map)


```{r, echo=FALSE, results='hide', message=FALSE, eval=FALSE}
library(stats19)
library(dplyr)
library(sf)
a = get_stats19(2018, "ac")
asf = format_sf(a)
a_zones = asf %>% 
  filter(local_authority_district == "Isle of Wight")
nrow(a_zones)
zones = pct::get_pct_zones(region = "isle-of-wight")
zones_osbg = st_transform(zones, 27700)
a_zones_sf = a_zones[zones_osbg, ]
nrow(a_zones_sf)
# mapview::mapview(zones) +
#   mapview::mapview(a_zones)
class(a$date)
class(a$time)
a_zones$month = lubridate::month(a_zones$date)
a_zones_may = a_zones %>% 
  filter(month == 5)
a_agg = aggregate(a_zones_may["speed_limit"], zones_osbg, mean)
plot(a_agg)
class(a$date)
```

# Analysing crash data on road networks

Road network data can be accessed from a range of sources, including OpenStreetMap (OSM) and Ordnance Survey.
We will use some OSM data from the Ilse of Wight, which can be loaded as follows:

```{r}
u = "https://github.com/ropensci/stats19/releases/download/1.1.0/roads_key.Rds"
roads_wgs = readRDS(url(u))
roads = roads_wgs %>% st_transform(crs = 27700)
```

You should already have road crashes for the Isle of Wight from the previous stage.
If not, load crash data as follows:

```{r}
u = "https://github.com/ropensci/stats19/releases/download/1.1.0/car_accidents_2017_iow.Rds"
crashes_iow = readRDS(url(u))
```

1. Plot the roads with the crashes overlaid.
2. Create a buffer around the roads with a distance of 200 m.
3. How many crashes fall outside the buffered roads?
3. Bonus: Use the `aggregate()` function to identify how many crashes happened per segment and plot the result (hint: see `?aggregate.sf` and take a read of Section [4.2.5](https://geocompr.robinlovelace.net/spatial-operations.html#spatial-aggr) of Geocomputation with R) with `tmap` and plot the crashes that happened outside the road buffers on top.

```{r, echo=FALSE, out.width="49%", fig.show='hold', message=FALSE}
plot(roads$geometry)
plot(crashes_iow["accident_severity"], add = TRUE)
roads_buffer = st_buffer(roads, 200, endCapStyle = "FLAT")
crashes_outside_roads = crashes_iow[roads_buffer, , op = sf::st_disjoint]
roads_agg = aggregate(crashes_iow[1], by = roads_buffer, FUN = length)
# plot(roads_agg, border = NA, main = "")
names(roads_agg)[1] = "N. Crashes"
tmap_mode("plot")
tm_shape(roads_agg) + tm_fill("N. Crashes") +
  tm_shape(crashes_outside_roads) + tm_dots(col = "blue")
```

\newpage

# Bonus exercises

Identify a region and zonal units of interest from http://geoportal.statistics.gov.uk/ or from the object `police_boundaries` in the `stats19` package.

1. Read them into R as an `sf` object
1. Create a map showing the number of crashes in each zone
1. Identify the average speed limit associated with crashes in each zone
1. Identify an interesting question you can ask of the data and use exploratory data analysis to find answers
1. Take a look at the code in [the file iow_example.R in the inst directory in the stats19 repo](https://github.com/ropensci/stats19/blob/master/inst/iow_example.R). Run it to reproduce the figure below.

```{r final-plot, echo=FALSE, out.width="100%"}
# knitr::include_graphics("final-figure.png")
```


# References
---
title: "Introducing stats19"
author: 
  - "R Lovelace, M Morgan, L Hama and M Padgham"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{stats19-intro}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: references.bib
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%",
  eval = curl::has_internet()
)
```

## Introduction

**stats19** enables access to and processing of Great Britain's official road traffic casualty database, [STATS19](https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data). 
A description of variables in the database can be found in a [document](https://data.dft.gov.uk/road-accidents-safety-data/Brief-guide-to%20road-accidents-and-safety-data.doc) provided by the UK's Department for Transport (DfT).
The datasets are collectively called STATS19 after the form used to report them, which can be found [here](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/995422/stats19.pdf).
This vignette focuses on how to use the **stats19** package to work with STATS19 data.

**Note**: The Department for Transport refers to "accidents", but "crashes" is a more appropriate term, as emphasised in the "crash not accident" arguments of road safety advocacy groups such as [RoadPeace](https://www.roadpeace.org/take-action/crash-not-accident/).
We use the term "accidents" only in reference to nomenclature within the data as provided.

The development version is hosted on [GitHub](https://github.com/ITSLeeds/stats19) and can be installed and loaded as follows:

```{r, eval=FALSE}
# from CRAN
install.packages("stats19")
# you can install the latest development (discoraged) using:
remotes::install_github("ITSLeeds/stats19")
```

```{r}
library(stats19)
```

## Functions

The easiest way to get STATS19 data is with `get_stats19()`.
This function takes 2 main arguments, `year` and `type`.
The year can be any year between 1979 and 202x where x is the current year minus one or two due to the delay in publishing STATS19 statistics.
The type can be one of `accidents`, `casualties` and `vehicles`, described below.
`get_stats19()` performs 3 jobs, corresponding to three main types of functions:

- **Download**: A `dl_stats19()` function accepts `year`, `type` and `filename` arguments to make it easy to find the right file to download only.

- **Read**: STATS19 data is provided in a particular format that benefits from being read-in with pre-specified column types. This is taken care of with `read_*()` functions providing access to the 3 main tables in STATS19 data:

    - `read_accidents()` reads-in the crash data (which has one row per incident)
    - `read_casualties()` reads-in the casualty data (which has one row per person injured or killed)
    - `read_vehicles()` reads-in the vehicles table, which contains information on the vehicles involved in the crashes (and has one row per vehicle)

- **Format**: There are corresponding `format_*()` functions for each of the `read_*()` functions. These have been exported for convenience, as the two sets of functions are closely related, there is also a `format` parameter for the `read_*()` functions, which by default is `TRUE`, adds labels to the tables. 
The raw data provided by the DfT contains only integers. Running `read_*(..., format = TRUE)` converts these integer values to the corresponding character variables for each of the three tables.
For example, `read_accidents(format = TRUE)` converts values in the `accident_severity` column from `1`, `2` and `3` to `Slight`, `Serious` and `Fatal` using `fromat_accidents()` function.
To read-in raw data without formatting, set `format = FALSE`.

Multiple functions (`read_*` and `format_*`) are needed for each step because of the structure of STATS19 data, which are divided into 3 tables:

1. "accident circumstances, with details about location, severity, weather, etc.; 
2. casualties, referencing knowledge about the victims; and
3. vehicles, which contains more information about the vehicle type and manoeuvres, as well the some information about the driver."

Data files containing multiple years worth of data can be downloaded.
Datasets since 1979 are broadly consistent, meaning that STATS19 data represents a rich historic geographic record of road casualties at a national level, as stated in the DfT's road casualties report in [2017](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/744077/reported-road-casualties-annual-report-2017.pdf):

> The current set of definitions and detail of information 
goes back to 1979, providing a long period for comparison.

## Download STATS19 data

**stats19** enables download of raw STATS19 data with `dl_*` functions.
The following code chunk, for example, downloads and unzips a .zip file containing STATS19 data from 2017:

```{r dl2017-accidents}
dl_stats19(year = 2017, type = "accident", ask = FALSE)
```

Note that in the previous command, `ask = FALSE`, meaning you will not be asked.
By default you are asked to confirm, before downloading large files.
Currently, these files are downloaded to a default location of `tempdir` which is a platform independent "safe" but temporary location to download the data in. Once downloaded, they are unzipped under original DfT file names.
The `dl_stats19()` function prints out the location and final file name(s) of unzipped files(s) as shown above.

`dl_stats19()` takes three parameters.
Supplying a `file_name` is interpreted to mean that the user is aware of what to download and the other two parameters will be ignored.
You can also use `year` and `type` to "search" through the file names, which are stored in a lazy-loaded dataset called `stats19::file_names`.

<!-- In versions above `1.3.0`, function `dl_stats19()` behaves more like `get_stats19()`, that is parameter `type` is by default `accidents` and it can also download multiple years more efficiently. -->

You can find out the names of files that can be downloaded with `names(stats19::file_names)`, an example of which is shown below:

```{r}
stats19::file_names$DigitalBreathTestData2013.zip
```

To see how `file_names` was created, see `?file_names`.
Data files from other years can be selected interactively.
Just providing a year, for example, presents the user with multiple options (from `file_names`), illustrated below:

```{r dl2017-all, eval=FALSE}
dl_stats19(year = 2017)
```

```
Multiple matches. Which do you want to download?

1: dft-road-casualty-statistics-casualty-2017.csv
2: dft-road-casualty-statistics-vehicle-2017.csv
3: dft-road-casualty-statistics-accident-2017.csv

Selection: 
Enter an item from the menu, or 0 to exit
```

When R is running interactively, you can select which of the 3 matching files to download:
those relating to vehicles, casualties or accidents in 2017.

## Read STATS19 data

In a similar approach to the download section before, we can read files downloaded using a `data_dir` location of the file and the `filename` to read.
The code below will download the `dftRoadSafetyData_Accidents_2017.zip` file from the DfT servers and read its content.
Files are saved by default in `tempdir()`, but this can be overridden to ensure permanent storage in a user-defined location.

```{r dl2017-read}
crashes_2017_raw = get_stats19(year = 2017, type = "acc", format = FALSE)
```

```{r, echo=FALSE}
# skip vignettes if resource unavailable
if(object.size(crashes_2017_raw) < 1000) {
  knitr::opts_chunk$set(eval = FALSE)
}
```


**stats19** imports data with `readr::read_csv()` which results in a 'tibble' object: a data frame with more user-friendly printing and a few other features.

```{r crashes2017-class}
class(crashes_2017_raw)
dim(crashes_2017_raw)
```

There are three `read_*()` functions, corresponding to the three different classes of data provided by the DfT:
1. `read_accidents()`
2. `read_casualties()`
3. `read_vehicles()`

In all cases, a default parameter `read_*(format = TRUE)` returns the data in formatted form, as described above. Data can also be imported in the form directly provided by the DfT by passing `format = FALSE`, and then subsequently formatted with additional `format_*()` functions, as described in a final section of this vignette. Each of these `read_*()` functions is now described in more detail.


### Crash data

After raw data files have been downloaded as described in the previous section, they can then be read-in as follows:

```{r read2017-raw-format}
crashes_2017_raw = read_accidents(year = 2017, format = FALSE)
crashes_2017 = format_accidents(crashes_2017_raw)
nrow(crashes_2017_raw)
ncol(crashes_2017_raw)
nrow(crashes_2017)
ncol(crashes_2017)
```

What just happened?
We read-in data on all road crashes recorded by the police in 2017 across Great Britain.
The dataset contains
`r # ncol(crashes_2017_raw)`
32
columns (variables) for 
`r # format(nrow(crashes_2017_raw), big.mark = ",")`
129,982 crashes.

This work was done by `read_accidents(format = FALSE)`, which imported the "raw" STATS19 data without cleaning messy column names or re-categorising the outputs.
`format_accidents()` function automates the process of matching column names with variable names and labels in a [`.xls` file](https://data.dft.gov.uk/road-accidents-safety-data/variable%20lookup.xls) provided by the DfT.
This means `crashes_2017` is much more usable than `crashes_2017_raw`, as shown below, which shows some key variables in the messy and clean datasets:

```{r crashes2017-columns}
crashes_2017_raw[c(7, 18, 23, 25)]
crashes_2017[c(7, 18, 23, 25)]
```

By default, `format = TRUE`, meaning that the two stages of `read_accidents(format = FALSE)` and `format_accidents()` yield the same result as `read_accidents(format = TRUE)`.
For the full list of columns, run `names(crashes_2017)`.

<!-- This means `crashes_2017` is much more usable than `crashes_2017_raw`, as shown below, which shows three records and some key variables in the messy and clean datasets: -->

```{r, echo=FALSE, eval=FALSE}
# commented out as confusing...
key_patt = "severity|speed|light|human"
key_vars = grep(key_patt, x = names(crashes_2017_raw), ignore.case = TRUE)
random_n = sample(x = nrow(crashes_2017_raw), size = 3)
crashes_2017_raw[random_n, key_vars]
crashes_2017[random_n, key_vars]
```

**Note**: As indicated above, the term "accidents" is only used as directly provided by the DfT; "crashes" is a more appropriate term, hence we call our resultant datasets `crashes_*`.

## Format STATS19 data

It is also possible to import the "raw" data as provided by the DfT.
A [`.xls` file](https://data.dft.gov.uk/road-accidents-safety-data/variable%20lookup.xls) provided by the DfT defines the column names for the datasets provided.
The packaged datasets `stats19_variables` and `stats19_schema` provide summary information about the contents of this data guide.
These contain the full variable names in the guide (`stats19_variables`) and a complete look up table relating integer values to the `.csv` files provided by the DfT and their labels (`stats19_schema`).
The first rows of each dataset are shown below:

```{r variables-and-schema}
stats19_variables
stats19_schema
```

The code that generated these small datasets can be found in their help pages (accessed with `?stats19_variables` and `?stats19_schema` respectively).
`stats19_schema` is used internally to automate the process of formatting the downloaded `.csv` files.
Column names are formatted by the function `format_column_names()`, as illustrated below:

```{r format-col-names}
format_column_names(stats19_variables$variable[1:3])
```

Previous approaches to data formatting `STATS19` data involved hard-coding results. This more automated approach to data cleaning is more consistent and fail-safe.
The three functions: `format_accidents()`, `format_vehicles()` and
`format_casualties()` do the data formatting on the respective data frames, as illustrated below:

```{r format-main}
crashes_2017 = format_accidents(crashes_2017_raw)

# vehicle data for 2017
dl_stats19(year = 2017, type = "vehicle", ask = FALSE)
vehicles_2017_raw = read_vehicles(year = 2017)
vehicles_2017 = format_vehicles(vehicles_2017_raw)

# casualties data for 2017
dl_stats19(year = 2017, type = "casualty", ask = FALSE)
casualties_2017 = read_casualties(year = 2017)
```

The package automates this two-step `read_*` and `format_*`
process by defaulting in all cases to `data_year = read_*(year, format = TRUE)`.
`read_*` functions return, by default, formatted data.
The two-step process may nevertheless be important for reference to the official nomenclature and values as provided by the DfT.

A summary of the outputs for each of the three tables is shown below.

```{r summarise-stats19}
summarise_stats19 = function(x) {
  data.frame(row.names = 1:length(x),
    name = substr(names(x), 1, 19),
    class = sapply(x, function(v) class(v)[1]),
    n_unique = sapply(x, function(v) length(unique(v))),
    first_label = sapply(x, function(v) substr(unique(v)[1], 1, 16)),
    most_common_value = sapply(x, function(v) 
      substr(names(sort(table(v), decreasing = TRUE)[1]), 1, 16)[1])
  )
}
```

```{r summarise-crashes}
knitr::kable(summarise_stats19(crashes_2017), 
             caption = "Summary of formatted crash data.")
```

```{r summarise-vehicles}
knitr::kable(summarise_stats19(vehicles_2017), 
             caption = "Summary of formatted vehicles data.")
```

```{r summarise-casualties}
knitr::kable(summarise_stats19(casualties_2017), 
             caption = "Summary of formatted casualty data.")
```

For testing and other purposes, a sample from the accidents table is provided in the package.
A few columns from the two-row sample is shown below:

```{r, echo=FALSE, results='asis'}
key_patt = "severity|speed|light|human"
key_vars = grep(key_patt, x = names(stats19::accidents_sample_raw), ignore.case = TRUE)
knitr::kable(stats19::accidents_sample_raw[, key_vars])
```

## Casualties data

As with `crashes_2017`, casualty data for 2017 can be downloaded, read-in and formatted as follows:

```{r 2017-cas}
dl_stats19(year = 2017, type = "casualty", ask = FALSE)
casualties_2017 = read_casualties(year = 2017)
nrow(casualties_2017)
ncol(casualties_2017)
```

The results show that there were 
`r # format(nrow(casualties_2017), big.mark=",")`
170,993
casualties reported by the police in the STATS19 dataset in 2017, and 
`r # ncol(casualties_2017)`
16
columns (variables).
Values for a sample of these columns are shown below:

```{r 2017-cas-columns}
casualties_2017[c(4, 5, 6, 14)]
```

The full list of column names in the `casualties` dataset is:

```{r 2017-cas-columns-all}
names(casualties_2017)
```

## Vehicles data

Data for vehicles involved in crashes in 2017 can be downloaded, read-in and formatted as follows:

```{r dl2017-vehicles}
dl_stats19(year = 2017, type = "vehicle", ask = FALSE)
vehicles_2017 = read_vehicles(year = 2017)
nrow(vehicles_2017)
ncol(vehicles_2017)
```

The results show that there were 
`r # format(nrow(vehicles_2017), big.mark=",")`
238,926
vehicles involved in crashes reported by the police in the STATS19 dataset in 2017, with 
`r # ncol(vehicles_2017)`
23
columns (variables).
Values for a sample of these columns are shown below:

```{r 2017-veh-columns}
vehicles_2017[c(3, 14:16)]
```

The full list of column names in the `vehicles` dataset is:

```{r 2017-veh-columns-all}
names(vehicles_2017)
```

<!-- More data can be read-in as follows: -->

```{r, eval=FALSE, echo=FALSE}
# old code to be up-dated
d14 = "Stats19_Data_2005-2014"
crashes_2005_2014 = read_accidents(data_dir = d14)
crashes_2005_2014_f = format_stats19_2005_2014_ac(crashes_2005_2014)
d15 = "RoadSafetyData_2015"
crashes_2015 = read_accidents(data_dir = d15, filename = "Accidents_2015.csv")
crashes_2015_f = format_stats19_2015_ac(crashes_2015)
d16 = "dftRoadSafety_Accidents_2016"
crashes_2016 = read_accidents(data_dir = d16, filename = "dftRoadSafety_Accidents_2016.csv")
crashes_2016_f = format_stats19_2016_ac(crashes_2016)
all_crashes = rbind(crashes_2015_f, crashes_2016_f, crashes_2017_f)
table(ac$Accident_Severity)
```

## Creating geographic crash data

An important feature of STATS19 data is that the "accidents" table contains geographic coordinates.
These are provided at ~10m resolution in the UK's official coordinate reference system (the Ordnance Survey National Grid, EPSG code 27700).
**stats19** converts the non-geographic tables created by `format_accidents()` into the geographic data form of the [`sf` package](https://cran.r-project.org/package=sf) with the function `format_sf()` as follows:

```{r format-crashes-sf}
crashes_sf = format_sf(crashes_2017)
```

The note arises because `NA` values are not permitted in `sf` coordinates, and so rows containing no coordinates are automatically removed.
Having the data in a standard geographic form allows various geographic operations to be performed on it.
Spatial operations, such as spatial subsetting and spatial aggregation, can be performed, to show the relationship between STATS19 data and other geographic objects, such as roads, schools and administrative zones.

An example of an administrative zone dataset of relevance to STATS19 data is the boundaries of police forces in England, which is provided in the packaged dataset `police_boundaries`.
The following code chunk demonstrates the kind of spatial operations that can be performed on geographic STATS19 data, by counting and plotting the number of fatalities per police force:

```{r nfatalities}
library(sf)
library(dplyr)
crashes_sf %>% 
  filter(accident_severity == "Fatal") %>% 
  select(n_fatalities = accident_index) %>% 
  aggregate(by = police_boundaries, FUN = length) %>% 
  plot()
```

Of course, one should not draw conclusions from such analyses without care.
In this case, denominators are needed to infer anything about road safety in any of the police regions.
After suitable denominators have been included, performance metrics such as 'health risk' (fatalities per 100,000 people), 'traffic risk' (fatalities per billion km, f/bkm) and 'exposure risk' (fatalities per million hours, f/mh) can be calculated [@feleke_comparative_2018; @elvik_handbook_2009].

The following code chunk, for example, returns all crashes within the jurisdiction of [West Yorkshire Police](https://en.wikipedia.org/wiki/West_Yorkshire_Police):

```{r ukboundaries}
west_yorkshire =
  police_boundaries[police_boundaries$pfa16nm == "West Yorkshire", ]
```


```{r crashes-west_yorkshire}
crashes_wy = crashes_sf[west_yorkshire, ]
nrow(crashes_sf)
nrow(crashes_wy)
```

This subsetting has selected the 
`r # format(nrow(crashes_wy), big.mark = ",")`
4,371
crashes which occurred in West Yorkshire.


## Joining tables

The three main tables we have just read-in can be joined by shared key variables.
This is demonstrated in the code chunk below, which subsets all casualties that took place in West Yorkshire, and counts the number of casualties by severity for each crash:

```{r table-join, message = FALSE}
library(tidyr)
library(dplyr)
sel = casualties_2017$accident_index %in% crashes_wy$accident_index
casualties_wy = casualties_2017[sel, ]
cas_types = casualties_wy %>% 
  select(accident_index, casualty_type) %>% 
  group_by(accident_index) %>% 
  summarise(
    Total = n(),
    walking = sum(casualty_type == "Pedestrian"),
    cycling = sum(casualty_type == "Cyclist"),
    passenger = sum(casualty_type == "Car occupant")
    ) 
cj = left_join(crashes_wy, cas_types)
```

What just happened? We found the subset of casualties that took place in West Yorkshire with reference to the `accident_index` variable.
Then we used the **dplyr** function `summarise()`, to find the number of people who were in a car, cycling, and walking when they were injured.
This new casualty dataset is joined onto the `crashes_wy` dataset.
The result is a spatial (`sf`) data frame of crashes in West Yorkshire, with columns counting how many road users of different types were hurt.
The joined data has additional variables:

```{r table-join-examples}
base::setdiff(names(cj), names(crashes_wy))
```

As a simple spatial plot, we can map all the crashes that have happened in West Yorkshire in 2017, with the colour related to the total number of people hurt in each crash.
Placing this plot next to a map of West Yorkshire provides context:

```{r, out.width="90%", fig.show='hold'}
plot(
  cj[cj$cycling > 0, "speed_limit", ],
  cex = cj$Total[cj$cycling > 0] / 3,
  main = "Speed limit (cycling)"
  )
plot(
  cj[cj$passenger > 0, "speed_limit", ],
  cex = cj$Total[cj$passenger > 0] / 3,
  main = "Speed limit (passenger)"
  )
```

The spatial distribution of crashes in West Yorkshire clearly relates to the region's geography.
Car crashes tend to happen on fast roads, including busy Motorway roads, displayed in yellow above.
Cycling is as an urban activity, and the most bike crashes can be found in near Leeds city centre, which has a comparatively high level of cycling (compared with the low baseline of 3%).
This can be seen by comparing the previous map with an overview of the area, from an academic paper on the social, spatial and temporal distribution of bike crashes [@lovelace_who_2016]:

```{r, echo=FALSE}
knitr::include_graphics("wy-overview.jpg")
```

In addition to the `Total` number of people hurt/killed, `cj` contains a column for each type of casualty (cyclist, car occupant, etc.), and a number corresponding to the number of each type hurt in each crash.
It also contains the `geometry` column from `crashes_sf`.
In other words, joins allow the casualties and vehicles tables to be geo-referenced.
We can then explore the spatial distribution of different casualty types.
The following figure, for example, shows the spatial distribution of pedestrians and car passengers hurt in car crashes across West Yorkshire in 2017:

```{r sfplot, fig.show='hold', out.width="100%", fig.cap="Spatial distribution of serious and fatal crashes in West Yorkshire, for cycling, walking, being a car passenger and other modes of travel. Colour is related to the speed limit where the crash happened (red is faster) and size is proportional to the total number of people hurt in each crash (legend not shown).", fig.width=9, fig.height=7}
library(ggplot2)
crashes_types = cj %>% 
  filter(accident_severity != "Slight") %>% 
  mutate(type = case_when(
    walking > 0 ~ "Walking",
    cycling > 0 ~ "Cycling",
    passenger > 0 ~ "Passenger",
    TRUE ~ "Other"
  ))
table(crashes_types$speed_limit)
ggplot(crashes_types, aes(size = Total, colour = speed_limit)) +
  geom_sf(show.legend = "point", alpha = 0.3) +
  facet_grid(vars(type), vars(accident_severity)) +
  scale_size(
    breaks = c(1:3, 12),
    labels = c(1:2, "3+", 12)
    ) +
  scale_color_gradientn(colours = c("blue", "yellow", "red")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
```

It is clear that different types of road users tend to get hurt in different places.
Car occupant casualties (labelled 'passengers' in the map above), for example, are comparatively common on the outskirts of cities such as Leeds, where speed limits tend to be higher and where there are comparatively higher volumes of motor traffic.
Casualties to people on foot tend to happen in the city centres.
That is not to say that cities centres are more dangerous per unit distance (typically casualties per billion kilometres, bkm, is the unit used) walked:
there is more walking in city centres (you need a denominator to estimate risk).

To drill down further, we can find the spatial distribution of all pedestrian casualties, broken-down by seriousness of casualty, and light conditions.
This can be done with **tidyvers** functions follows:

```{r ggplot-ped-severity, fig.height=5, fig.width=6}
table(cj$light_conditions)
cj %>% 
  filter(walking > 0) %>% 
  mutate(light = case_when(
    light_conditions == "Daylight" ~ "Daylight",
    light_conditions == "Darkness - lights lit" ~ "Lit",
    TRUE ~ "Other/Unlit"
  )) %>% 
  ggplot(aes(colour = speed_limit)) +
  geom_sf() +
  facet_grid(vars(light), vars(accident_severity)) +
  scale_color_continuous(low = "blue", high = "red") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
```

<!-- These figures show the issue of pedestrian casualties in Leeds and Bradford city centres, which both have busy roads going through or around popular shopping and other destinations. -->

## Time series analysis

We can also explore seasonal and daily trends in crashes by aggregating crashes by day of the year:

```{r crash-date-plot, fig.width=5, fig.height=5}
crashes_dates = cj %>% 
  st_set_geometry(NULL) %>% 
  group_by(date) %>% 
  summarise(
    walking = sum(walking),
    cycling = sum(cycling),
    passenger = sum(passenger)
    ) %>% 
  gather(mode, casualties, -date)
ggplot(crashes_dates, aes(date, casualties)) +
  geom_smooth(aes(colour = mode), method = "loess") +
  ylab("Casualties per day")
```


Different types of crashes also tend to happen at different times of day.
This is illustrated in the plot below, which shows the times of day when people who were travelling by different modes were most commonly injured.

```{r crash-time-plot, fig.width=5, fig.height=5}
library(stringr)

crash_times = cj %>% 
  st_set_geometry(NULL) %>% 
  group_by(hour = as.numeric(str_sub(time, 1, 2))) %>% 
  summarise(
    walking = sum(walking),
    cycling = sum(cycling),
    passenger = sum(passenger)
    ) %>% 
  gather(mode, casualties, -hour)

ggplot(crash_times, aes(hour, casualties)) +
  geom_line(aes(colour = mode))
```

Note that bike crashes tend to have distinct morning and afternoon peaks, in-line with previous research [@lovelace_who_2016].
A disproportionate number of car crashes appear to happen in the afternoon.

## Further work

There is much potential to extend the package beyond downloading, reading and formatting STATS19 data.
The greatest potential is to provide functions that will help with analysis of STATS19 data, to help with road safety research.
Much academic research has been done using the data, a few examples of which are highlighted below to demonstrate the wide potential for further work.

- Research exploring the effectiveness of road safety policies such as speed limits. An example in this area is this [paper on 20 mph zones](https://pubmed.ncbi.nlm.nih.gov/20007666/) who found that areas with 20mph speed limits were safer. This raises the question: can the same result be repeated using reproducible methods? Does the finding hold for more recent 20 mph zones? Is the recent finding of the Department for Transport's ([2018](https://www.gov.uk/government/publications/20-mph-speed-limits-on-roads)) research, that 20 mph zones alone do not reduce crash rates, supported by reproducible analysis? What are the factors that make speed limits more or less effective [see @sarkar_street_2018 for example]?
- Research into weather as a contributing factor to road traffic casualties [e.g. @edwards_relationship_1998]. This raises the question: could matching crash data from the STATS19 data with historic weather data from other R packages help advance knowledge in this area?
- Assessment of crash rates normalised by estimated exposure rates (risk). An example of this type of research by an author of the package found substantial spatial variation in the number of cyclist casualties across West Yorkshire [@lovelace_who_2016]. This raises the questions: are similar spatial differences found in other regions? What are the factors leading to relatively high and low rates of different types of crash? 

The broader point is that the **stats19** package could help road safety research, by making open access data on road crashes more accessible to researchers worldwide.
By easing the data download and cleaning stages of research, it could also encourage reproducible analysis in the field.

There is great potential to add value to and gain insight from the data by joining the datasets with open data, for example from the Consumer Data Research Centre ([CDRC](https://www.cdrc.ac.uk/), which funded this research), OpenStreetMap and the UK's Ordnance Survey.
If you have any suggestions on priorities for these future directions of (hopefully safe) travel, please get in touch on at [github.com/ITSLeeds/stats19/issues](https://github.com/ITSLeeds/stats19/issues).

## References
---
title: "An introduction to road safety analysis with R: setup notes"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{stats19-training-setup}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "50%"
)
```

If you are not experienced with R, it is strongly advised that you read-up on and more importantly **test out** R and RStudio before attempting analyse road crash data with R.

To read up on R, we recommend reading Chapter 1 Getting Started with Data in R of the online book Statistical Inference via Data Science, which can be found here: https://moderndive.netlify.app/1-getting-started.html

Reading sections 1.1 to 1.3 of that book and trying a few of the examples are considered **essential prerequisites**, unless you are already experienced with R.

Optionally, if you want a more interactive learning environment, you can try getting started with online resources, such as those found at [education.rstudio.com/learn](https://education.rstudio.com/learn/beginner/).

And for more information on how R can be used for transport research, the Transportation chapter of Geocomputation with R is a good place to start: https://geocompr.robinlovelace.net/transport.html

**Your computer should also have the necessary software installed.**

To ensure your computer is ready for the course, you should have a recent (3.6.0 or later) version of R or RStudio installed. You should have installed packages stats19, tidyverse and a few others shown below.
To check you have the necessary packages installed, try running the following line of code:

```{r, eval=FALSE}
source("https://git.io/JeaZH")
```

That does some basic checks.
For more comprehensive checkes, and to get used to typing in R code, you can also test your setup by typing and executing the following lines in the RStudio console (this will install the packages you need if they are not already installed):

```{r, eval=FALSE}
install.packages("remotes")
pkgs = c(
  "pct",         # package for getting travel data in the UK
  "sf",          # spatial data package
  "stats19",     # downloads and formats open stats19 crash data
  "stplanr",     # for working with origin-destination and route data
  "tidyverse",   # a package for user friendly data science
  "tmap"         # for making maps
)
remotes::install_cran(pkgs)
# remotes::install_github("ITSLeeds/pct")
```

To test your computer is ready to work with road crash data in R, try running the following commands from RStudio (which should result in the map below):

 <!-- method for helping people set up their computers. Type this single line into the console and follow the instructions.  -->

```{r message=FALSE, eval=FALSE}
library(stats19)
library(tidyverse)
library(tmap) # installed alongside mapview
crashes = get_stats19(year = 2017, type = "ac")
crashes_iow = crashes %>% 
  filter(local_authority_district == "Isle of Wight") %>% 
  format_sf()
  
# basic plot
plot(crashes_iow)
```

You should see results like those shown in the map here: https://github.com/ropensci/stats19/issues/105

If you cannot create that map by running the code above before the course, get in touch with us, e.g. by writing a comment under that github issue page (Note: You will need a github account). 

# Time

Perhaps the most important pre-requisite is time.
You'll need to find time to work-through these materials, either in one go (see suggested agenda below) or in chunks of perhaps 1 hour per week over a 2 month period.
I think ~8 hours is a good amount of time to spend on this course but it can be done in small pieces, e.g.:

- If you're totally new to R a 3 hour session with a one hour break to do sections 1 to 4
- If you're an intermediant user, you could skim through sections 1:3 and focus on sections 4:6 to gain understanding of the temporal and spatial aspects of the data in a few hourse
- If you're an advanced user, feel free to skip ahead and work through the sections you find most interesting.

# 2 day course agenda

For the more structured 2 day course for R beginners, a preliminary agenda is as follows:

## Day 1: An introduction to R and RStudio for spatial and temporal data

09:00-09:30 Arrival and set-up

09:30-11:00 Introduction to the course and software

- Introduction to R + coding [video](https://youtu.be/7oyiPBjLAWY?t=357)/[links](https://github.com/jennybc/code-smells-and-feels)
- R installation questions/debugging
- How to use RStudio (practical in groups of 2)
- R classes and working with data frames (CC)

**Break**

11:15-12:30 Working with temporal data

- Time classes
- Filtering by time of crash
- Aggregating over time
- Forecasting crashes over time

**Lunch**

13:30-15:00 Working with spatial data

- Spatial data in R
- Context: spatial ecosystem (see section [1.4 of Geocomputation with R - package ecosystem](https://geocompr.robinlovelace.net/intro.html#rs-spatial-ecosystem))
- [Exercises](https://geocompr.robinlovelace.net/attr.html#exercises-1): Section 6 of the handout

- Further reading: [Section 3.2 to 3.2.2](https://geocompr.robinlovelace.net/attr.html#vector-attribute-manipulation) of handouts

**Break**
  
15:15-15:30 Talk on Road Safety 1

15:30-16:15 Practical - Applying the methods to stats19 data - a taster

- How to access data with **stats19**
- Key **stats19** functions
- Excercises: analysing road crash data on the Isle of Wight

16:15-16:30 Talk on Road Safety 2

## Day 2 road safety analysis with R

09:30-11:00 Point pattern analysis

- Visualising data with tmap
- Spatial and temporal subsetting
- Aggregation

11:15-12:30 Road network data

- Desire lines: using origin-destination data
- Downloading road network data from OSM
- Buffers on road networks

**Lunch**

13:30-15:00 Analysing crash data on road network

**Break**
  
15:15-15:30: Talk on Road Safety 3

15:30-16:30 Applying the methods to your own data

---
title: "stats19: a package for road safety research"
author: 
  - "Layik Hama and Robin Lovelace"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{stats19: a package for road safety research}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: references.bib
always_allow_html: yes
---

#### Introduction

`stats19` is a new R package that enables access to and processing of
Great Britain’s official road traffic casualty database,
[STATS19](https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data).

We started the package in late 2018 following three main motivations:

1.  The release of the 2017 road crash statistics, which showed
    worsening road safety in some areas, increasing the importance of
    making the data more accessible.
2.  The realisation that many researchers were writing add-hoc code to
    clean the data, with a huge amount of duplicated (wasted) effort and
    potential for mistakes to lead to mistakes in the labelling of the
    data (more on that below).
3.  An understanding of the concept of ‘modularity’ in software design,
    following the [Unix
    philosophy](https://en.wikipedia.org/wiki/Unix_philosophy) that
    programs should ‘do one thing and do it well’. This realisation has
    led to code inside the rOpenSci-hosted package
    [`stplanr`](https://github.com/ropensci/stplanr) being split-out
    into 2 separate packages already:
    [`cyclestreets`](https://github.com/Robinlovelace/cyclestreets) and
    [`stats19`](https://github.com/ropensci/stats19).

We have a wider motivation: we want the roads to be safer. By making
data on the nature of road crashes in the UK more publicly accessible
(in a usable format), we hope this package saves lives.

It has been peer reviewed thanks to rOpenSci and is now published in the
Journal of Open Source Software
([JOSS](https://joss.theoj.org/papers/10.21105/joss.01181)) (Lovelace et
al. 2019). For installation and the code, see its home in rOpenSci:
<https://github.com/ropensci/stats19>, and the package documentation at
<https://itsleeds.github.io/stats19/>.

In this post, we’ll provide a bit of context, show how the package works, and provide ideads for future work building on the experience. Now is a good time to report on the package: version [`0.2.0`](https://cran.r-project.org/package=stats19) has just been release on CRAN, which contains a few improvements, some of which are used in this blog post (`ask = FALSE` in `get_stats19()`, for example, which makes it even quicker to get data with this package).

#### Short history and affiliations

The main authors are based at the [Institute for Transport Studies
(ITS)](https://environment.leeds.ac.uk/transport) and the [Leeds
Institute for Data Analytics (LIDA)](https://lida.leeds.ac.uk/),
institutions focussed on transport and data-intensive research. We have
prior experience writing code to work with road crash data: Robin wrote
code for an academic paper on cycle safety in Yorkshire based on STATS19
data (Lovelace, Roberts, and Kellar 2016), and put it in the
[robinlovelace/bikeR](https://github.com/Robinlovelace/bikeR) repo for
posterity/reproducibility. Package contributor, Dr Malcolm Morgan, wrote
code processing different STATS19 data for the Cycling Infrastructure
Prioritisation Toolkit ([CyIPT](https://www.cyipt.bike)) and put it in
[cyipt/stats19](https://github.com/cyipt/stats19). 

The large and complex STATS19 data from the UK's Department for Transport, which is open access but difficult-to-use, represented a perfect opportunity for us to get stuck into a chunky data processing challenge.

#### What is STATS19 anyway?

One of the reviewer comments from the [rOpenSci review process (which
has 77 comments - lots of knowledge
shared)](https://github.com/ropensci/software-review/issues/266) alluded
to the package’s esoteric name:

> I confess I wish the package name was more expressive–stats19 sounds
> like an introductory statistics class.

We agree\! However, the priority with the package is to remain faithful
to the data, and alternative name options, such as stats19data,
roadcrashesUK and roadSafetyData were not popular. Furthermore, the term
‘stats19’ is strongly associated with road crash data online. The URL
<https://en.wikipedia.org/wiki/STATS19> resolves to
<https://en.wikipedia.org/wiki/Reported_Road_Casualties_Great_Britain>,
for example (this page is provides and excellent introduction to
STATS19).

The name comes from a UK police form called
[STATS19](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/775149/Operational_Metrics_Manual.pdf)
(note the capital letters). There is another document called [STATS20](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/995423/stats20-2011.pdf), which is the guidance to STATS19 officers filling in STATS19 form. 

An important point is that the dataset omits crashes in which
nobody was hurt.
The Department for Transport (DfT) also names the dataset [STATS19 on the main web page that links to open access road crash data](https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data).

The importance of road safety and informed decision making based on crash data cannot be overstated. Deliberately avoiding the matter of life and death of road safety, two numbers from a strategy [document](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/8146/strategicframework.pdf) by the UK government (2011) are worth mentioning to show the scale of the numbers: 

> The economic welfare costs [of road collisions] are estimated at around £16 billion a year while insurance payouts for motoring claims alone are now over £12 billion a year.

Even more shocking are the [global statistics](https://www.who.int/data/gho/data/themes/road-safety), as summarised by an open access and reproducible academic [paper that uses data from the package to explore car-pedestrian crashes](https://github.com/Robinlovelace/stats19-gisruk):

> While so many people die on the roads each year in the UK (1,793 people in 2017, 3 deaths per 100,000) and worldwide (1,250,000 people in 2015, 17 deaths per 100,000) and ‘vision zero’ remains a Swedish dream (Johansson 2009), we urge people researching STATS19 and other road safety datasets to focus on a more urgent question: how to stop this carnage?

#### The road crash data in stats19

There are three main different types of CSV files released by the DfT: `accidents`, `vehicles` and `casualties` tables.
There is a schema covering these tables but a good amount of work is needed to understand it, let alone be able to the data contained within the files and convert the integers they contain into meaningful data.

The annual statistics have not been released in a consistent way either, making it hard for people to download, or even find, the relevant files.
For example, there are separate files for each of the above tables for certain years (e.g. 2016, 2017) but not for all of 1979 - 2017 or 2018 now.
The largest chunk is the 1979 - 2004 data, which is made available in a huge ZIP file ([link](https://data.dft.gov.uk/road-accidents-safety-data/Stats19-Data1979-2004.zip)).
Unzipped this contains the following 3 files, which occupy almost 2 GB on your hard drive: 

```sh
721M Apr  3  2013 Accidents7904.csv
344M Apr  3  2013 Casualty7904.csv
688M Apr  3  2013 Vehicles7904.csv
# total 1.753 GB data
```

#### Note

As of February 2021 we do not run any of the code in this vignette due to CRAN check issues.
See https://ropensci.org/blog/2019/02/26/stats19/ for the final version of this vignette.

```{r}
knitr::opts_chunk$set(eval = FALSE)
```


#### How stats19 works

With those introductions out of the way, lets see what the package can do and how it can empower people (including target audiences of academics, policy makers and road safety campaigners and R programmers/analysts wanting to test their skills on a large spatio-temporal dataset) to access STATS19 data, back to 1979.
First install the package in the usual way:

```{r, eval=FALSE, message=FALSE}
# release version - currently 0.2.0
install.packages("stats19") 
# dev version
# remotes::install_github("ropensci/stats19") 
```

Attach the package as follows:

```{r}
library(stats19)
```

The easiest way to get STATS19 data is with `get_stats19()`.
This function takes 2 main arguments, `year` and `type`.
The year can be any year between 1979 and 2017.

```{r dl2017-accidents, message=FALSE}
crashes_2017 = get_stats19(year = 2017, type = "Accidents", ask = FALSE)
nrow(crashes_2017)
```

What just happened?
We just downloaded, cleaned and read-in data on all road crashes recorded by the police in 2017 across Great Britain. We can explore the `crashes_2017` object a little more:

```{r crashes_2017-explore}
column_names = names(crashes_2017)
length(column_names)
head(column_names)
class(crashes_2017)
kableExtra::kable(head(crashes_2017[, c(1, 4, 5, 7, 10)]))
```

The package contains the names of all "zip" files released by the DfT and hosted on Amazon servers to download. These file names have been included in the package and can be found under `file_names` variable name. for example:

```{r file_names}
stats19::file_names$dftRoadSafetyData_Vehicles_2017.zip
```

You can also get the raw data (if you really want!) to see how much more useful the data is after it has been cleaned and labelled by the `stats19` package, compared with the data provided by government:

```{r}
crashes_2017_raw = get_stats19(year = 2017, type = "Accidents", ask = FALSE, format = FALSE)
```
The first two columns are raw read, the next two are formatted by `stats19` package:

```{r crashes_2017-raw}
kableExtra::kable(cbind(head(crashes_2017_raw[1:2, c(7, 10)]), head(crashes_2017[1:2, c(7, 10)])))
class(crashes_2017_raw$Date)
class(crashes_2017$date)
```

Note: the severity type is not labelled (this problem affects dozens of columns), the column names are inconsistent, and the dates have note been cleaned and converted into a user-friendly date (`POSIXct`) class:


```{r}
class(crashes_2017$date)
class(crashes_2017_raw$Date)
```

#### Creating geographic crash data

An important feature of STATS19 data is that the "accidents" table contains geographic coordinates.
These are provided at ~10m resolution in the UK's official coordinate reference system (the Ordnance Survey National Grid, EPSG code 27700).
**stats19** converts the non-geographic tables created by `format_accidents()` into the geographic data form of the [`sf` package](https://cran.r-project.org/package=sf) with the function `format_sf()` as follows:

```{r format-crashes-sf}
crashes_sf = format_sf(crashes_2017)
# crashes_sf = format_sf(crashes_2017, lonlat = TRUE) # provides the data in lon/lat format
```

An example of an administrative zone dataset of relevance to STATS19 data is the boundaries of police forces in England, which is provided in the packaged dataset `police_boundaries`.
The following code chunk demonstrates the kind of spatial operations that can be performed on geographic STATS19 data, by counting and plotting the number of fatalities per police force:

```{r nfatalities, message=FALSE}
library(sf)
library(dplyr)
crashes_sf %>% 
  filter(accident_severity == "Fatal") %>% 
  select(n_fatalities = accident_index) %>% 
  aggregate(by = police_boundaries, FUN = length) %>% 
  plot()
```

```{r ukboundaries}
west_yorkshire =
  police_boundaries[police_boundaries$pfa16nm == "West Yorkshire", ]
```


```{r crashes-west_yorkshire}
crashes_wy = crashes_sf[west_yorkshire, ]
nrow(crashes_wy) # which is 3.36%
```

#### The big picture: road safety

We can combine the three sets of tables to analyse the data further. Lets read the datasets first:

```{r dl2017-vehcas, message=FALSE}
#crashes_2017 = get_stats19(year = 2017, type = "Accidents", ask = FALSE)
casualties_2017 = get_stats19(year = 2017, type = "casualty", ask = FALSE)
nrow(casualties_2017)
vehicles_2017 = get_stats19(year = 2017, type = "vehicle", ask = FALSE)
nrow(vehicles_2017)
```

Lets now read in casualties that took place in West Yorkshire (using `crashes_wy` object above), and count the number of casualties by severity for each crash:

```{r table-join, message = FALSE}
library(tidyr)
library(dplyr)
sel = casualties_2017$accident_index %in% crashes_wy$accident_index
casualties_wy = casualties_2017[sel, ]
cas_types = casualties_wy %>% 
  select(accident_index, casualty_type) %>% 
  group_by(accident_index) %>% 
  summarise(
    Total = n(),
    walking = sum(casualty_type == "Pedestrian"),
    cycling = sum(casualty_type == "Cyclist"),
    passenger = sum(casualty_type == "Car occupant")
    ) 
cj = left_join(crashes_wy, cas_types)
```

What just happened? 

We found the subset of casualties that took place in West Yorkshire with reference to the `accident_index` variable in the `accidents` table.
Then we used the **dplyr** function `summarise()`, to find the number of people who were in a car, cycling, and walking when they were injured.
This new casualty dataset is joined onto the `crashes_wy` dataset.
The result is a spatial (`sf`) data frame of crashes in West Yorkshire, with columns counting how many road users of different types were hurt.
The joined data has additional variables:

```{r table-join-examples}
base::setdiff(names(cj), names(crashes_wy))
```

In addition to the `Total` number of people hurt/killed, `cj` contains a column for each type of casualty (cyclist, car occupant, etc.), and a number corresponding to the number of each type hurt in each crash.
It also contains the `geometry` column from `crashes_sf`.
In other words, joins allow the casualties and vehicles tables to be geo-referenced.
We can then explore the spatial distribution of different casualty types.
The following figure, for example, shows the spatial distribution of pedestrians and car passengers hurt in car crashes across West Yorkshire in 2017:

```{r sfplot, fig.show='hold', out.width="100%", fig.cap="Spatial distribution of serious and fatal collisions in which people who were walking on the road network ('pedestrians') were hit by a car or other vehicle.", fig.width=9, fig.height=7}
library(ggplot2)
crashes_types = cj %>% 
  filter(accident_severity != "Slight") %>% 
  mutate(type = case_when(
    walking > 0 ~ "Walking",
    cycling > 0 ~ "Cycling",
    passenger > 0 ~ "Passenger",
    TRUE ~ "Other"
  ))
ggplot(crashes_types, aes(size = Total, colour = speed_limit)) +
  geom_sf(show.legend = "point", alpha = 0.3) +
  facet_grid(vars(type), vars(accident_severity)) +
  scale_size(
    breaks = c(1:3, 12),
    labels = c(1:2, "3+", 12)
    ) +
  scale_color_gradientn(colours = c("blue", "yellow", "red")) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
```

To show what is possible when the data are in this form, and allude to next steps, let's create an interactive map.
We will plot crashes in which pedestrians were hurt, from the `crashes_types`, using `leaflet` package:

```{r crashes-map, fig.show='hold', out.width="100%", fig.cap="Spatial distribution of all collisions in which people who were walking on the road network ('pedestrians') were hit by a car or other vehicle in 2017 within West Yorkshire boundary.", fig.width=9, fig.height=7}
library(leaflet)
crashes_pedestrians = crashes_types %>% 
  filter(walking > 0)
# convert to lon lat CRS
crashes_pedestrians_lonlat = st_transform(crashes_pedestrians, crs = 4326)
pal = colorFactor(palette = "Reds", domain = crashes_pedestrians_lonlat$accident_severity, reverse = TRUE)
map = leaflet(data = crashes_pedestrians_lonlat, height = "280px") %>%
  addProviderTiles(provider = providers$OpenStreetMap.BlackAndWhite) %>%
  addCircleMarkers(radius = 0.5, color = ~pal(accident_severity)) %>% 
  addLegend(pal = pal, values = ~accident_severity) %>% 
  leaflet::addMiniMap(toggleDisplay = TRUE)
# map # if you like to see the leaflet version
```

```{r custom-leaflet}
library(geojsonsf)
library(htmltools)
geojson = sf_geojson(
  crashes_pedestrians_lonlat[,c("accident_severity")])
template = paste0('
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.4.0/dist/leaflet.css" integrity="sha512-puBpdR0798OZvTTbP4A8Ix/l+A4dHDD0DGqYW6RQ+9jxkRFclaxxQb/SJAWZfWAkuyeQUytO7+7N4QKrDh+drA==" crossorigin=""/>
<script src="https://unpkg.com/leaflet@1.4.0/dist/leaflet.js" integrity="sha512-QVftwZFqvtRNi0ZyCtsznlKSWOStnDORoefr1enyq5mVL4tmKB3S/EnC3rRJcxCPavG10IcrVGSmPh6Qw5lwrg==" crossorigin=""></script>
<div id="mapid" style="width: 100%; height: 400px;">
<script>
	var map = L.map("mapid");
	L.tileLayer("https://api.tiles.mapbox.com/v4/{id}/{z}/{x}/{y}.png?access_token=pk.eyJ1IjoibWFwYm94IiwiYSI6ImNpejY4NXVycTA2emYycXBndHRqcmZ3N3gifQ.rJcFIG214AriISLbB6B5aw", {
		maxZoom: 18,
		attribution: \'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors, <a href="https://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>, Imagery © <a href="https://www.mapbox.com/">Mapbox</a>\',
		id: "mapbox.streets"
  }).addTo(map);   
  var json = ', geojson, ';', 
  '
  var geojsonMarkerOptions = {
    radius: 6,
    color: "#000",
    weight: 1,
    opacity: 1,
    fillOpacity: 0.5
  };
  var layer = L.geoJSON(json, {
    pointToLayer: function (feature, latlng) {
        return L.circleMarker(latlng, geojsonMarkerOptions);
    },
    style: function(feature) {
      switch (feature.properties.accident_severity) {
          case "Serious": return {color: "#FEB24C"};
          case "Fatal":   return {color: "#BD0026"};
       }
    }
  }).addTo(map);
  map.fitBounds(layer.getBounds());
  var legend = L.control({position: "bottomright"});
	legend.onAdd = function (map) {
		var div = L.DomUtil.create("div", "info legend"), labels = [];
    labels.push("<i style=\'background:#FEB24C\'></i>Serious");
    labels.push("<i style=\'background:#BD0026\'></i>Fatal");
		div.innerHTML = labels.join("<br>");
		return div;
	};
  legend.addTo(map);
  // control that shows state info on hover
	var info = L.control();
	info.onAdd = function (map) {
		this._div = L.DomUtil.create("div", "info");
		this.update();
		return this._div;
	};
	info.update = function (props) {
		this._div.innerHTML = "<h6>Crashes in West Yorkshire (2017)</h6>";
	};
	info.addTo(map);
</script>
<style>
.info { padding: 6px 8px; font: 14px/16px Arial, Helvetica, sans-serif; background: white; background: rgba(255,255,255,0.8); box-shadow: 0 0 15px rgba(0,0,0,0.2); border-radius: 5px; } .info h4 { margin: 0 0 5px; color: #777; }
.legend { text-align: left; line-height: 18px; color: #555; } .legend i { width: 18px; height: 18px; float: left; margin-right: 8px; opacity: 0.7; }</style>
</div>')
path = file.path(tempdir(), "temp.html")
write(template, path)
includeHTML(path)
```

#### Conclusion

[`stats19`](https://github.com/ropensci/stats19) provides access to a reliable and official road safety dataset.
As covered in this post, it helps with data discovery, download, cleaning and formatting, allowing you to focus on the real work of analysis (see the package's introductory [vignette](https://itsleeds.github.io/stats19/articles/stats19.html) for more on this).
In our experience, 80% of time spent using STATS19 data was spent on data cleaning.
Hopefully, now the data is freely available in a more useful form, 100% of the time can be spent on analysis!
We think it could help many people, especially, including campaigners, academics and local authority planners aiming to make the roads of the UK and the rest of the world a safe place for all of us.

There are many possible next steps, including:

- Comparing these datasets with interventions such as [20 mph zones](https://pubmed.ncbi.nlm.nih.gov/20007666/) and links with street morphology [@sarkar_street_2018].
- The creation of more general software for accessing and working with road crash data worldwide.
- Making the data even more available by provide the data as part of an interactive web application, a technique successfully used in the Propensity to Cycle Tool (PCT) project hosted at [www.pct.bike/](https://www.pct.bike/) (this would likely take further times/resources beyond what we can provide in our spare time!).

For now, however, we want to take the opportunity to celebrate the release of `stats19` 🎉, thank rOpenSci for a great review process 🙏 and let you know: the package and data are now out there, and are ready to be used 🚀.

#### References
---
title: "Introduction to R for road safety: an introduction to R and practical exercises"
subtitle: "Practical exercises based on UK data from the stats19 package, developed by the Institute for Transport Studies, University of Leeds"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introducing stats19}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
author: Robin Lovelace, Malcolm Morgan & Andrea Gilardi
bibliography:
  - references.bib
  - packages.bib
---

# Note

This introductory practical vignette has largely been replaced by the workbook "Reproducible Road Safety Research with R", which can be found at https://itsleeds.github.io/rrsrr/ and as a PDF at [racfoundation.org](https://www.racfoundation.org/wp-content/uploads/Reproducible_road_safety_research_with_R_Lovelace_December_2020.pdf)  [@lovelace_reproducible_2020].
That workbook is more comprehensive than the content in this tutorial and is the recommended place to learn about using R for reproducible road safety research.


# Introduction

This document provides information, code and exercises to test and improve your R skills with an emphasis on road safety research.
It was initially developed to support a [2 day course](https://www.racfoundation.org/introduction-to-r-for-road-safety).
The course is based on open road crash records from the **stats19** package [@lovelace_stats19_2019].
However, the content should be of use for anyone working with road crash data that has (at a minimum):

- A timestamp
- A location (or address that can be geocoded)
- Attribute data, such as severity of crash, mode of vehicles involved etc.

You should type, run and ensure you understand each line of code in this document.

Code and data supporting the content can be found in the package's GitHub repo at [github.com/ropensci/stats19](https://github.com/ropensci/stats19/).
The '[issue tracker](https://github.com/ropensci/stats19/issues)' associated with that repo is a good place to ask questions about the course. 

## Prerequisites

If you are not experienced with R, it is strongly advised that you read-up on and more importantly **test out** R and RStudio before attempting analyse road crash data with R.
See the `stats19-training-setup` vignette at https://docs.ropensci.org/stats19/articles/stats19-training-setup.html for guidance on getting started with R, RStudio and installing R packages.

The completing the course requires that the following packages, which can be installed with `install.packages()`, can be loaded as follows:

```{r, message=FALSE, warning=FALSE, eval=FALSE}
library(pct)      # access travel data from DfT-funded PCT project 
library(sf)       # spatial vector data classes
library(stats19)  # get stats19 data
library(stplanr)  # transport planning tools
library(tidyverse)# packages for 'data science'
library(tmap)     # interactive maps
```

You should type, run and ensure you understand each line of code in this document.


```{r, include = FALSE}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "50%",
  eval = curl::has_internet()
)
```

```{r pkgs, warning=FALSE, echo=FALSE}
pkgs = c(
  "sf",          # spatial data package
  "stats19",     # downloads and formats open stats19 crash data
  "dplyr",       # a package data manipulation, part of the tidyverse
  "tmap"         # for making maps
)
```

```{r cite, echo=FALSE}
knitr::write_bib(x = pkgs, "packages.bib")
```

```{r, eval=FALSE, echo=FALSE}
remotes::install_cran(pkgs)
# remotes::install_github("ITSLeeds/pct")
```

# R and RStudio

The learning outcomes of this first session are to learn: 
RStudio main features and scripts,
R objects and functions,
subsetting,
basic plotting, and
getting help.

The first exercise is to open up RStudio and take a look around and identify the main components, shown in the figure below.
**Explore each of the main components of RStudio.**
Try changing the Global Settings (in the Tools menu) and see RStudio's short cuts by pressing `Alt-Shift-K` (or `Option+Shift+K` on Mac).

```{r rstudioui, echo=FALSE, out.width="70%"}
knitr::include_graphics("https://raw.githubusercontent.com/ITSLeeds/TDS/master/courses/2day/images/rstudio-ui.png")
```

## Projects and scripts

Projects are a way to organise related work together. Each project has its own folder and Rproj file. **Advice: always working from projects will make your life easier!** Start a new project with:

> File > New Project
You can choose to create a new directory (folder) or associate a project with an existing directory. Make a new project called stats1-course and save it in a sensible place on your computer. Notice that stats1-course now appears in the top right of RStudio.

Scripts are the files where R code is stored.
**Keeping your code in sensibly named, well organised and reproducible scripts will make your life easier:**
you could simply type all our code into the console, but that require retyping commands each time you run it.
Instead, code that you want to keep and share should be saved script files, plain text files that have the `.R` extension.

Make a new script with Flie > New File > Rscript or Ctrl+Shift+N

Save the script and give it a sensible name like `stats19-lesson-1.R` with File > Save, the save button on the toolbar, or Ctrl+S.

**Pro tip:** You can also create new R scripts by typing and running this command in the R console:

```{r edit, eval=FALSE}
file.edit("stats19-lesson-1.R")
```

Keeping scripts and other files associated with a project in a single folder per project (in an RStudio project) will help you find things you need and develop an efficient workflow.

## Writing and running code

Let's start with some basic R operations.
Write this code into your new `stats19-lesson-1.R` R script and execute the result line-by-line by pressing Ctrl+Enter

```{r, eval=FALSE}
x = 1:5
y = c(0, 1, 3, 9, 18)
plot(x, y)
```

This code creates two objects, both are vectors of 5 elements, and then plots them (bonus: check their length using the `length()` function).
Save the script by pressing Ctrl+S.

There are several ways to run code within a script and it is worth becoming familiar with each.
Try running the code you saved in the previous section using each of these methods:

1. Place the cursor in different places on each line of code and press `Ctrl+Enter` to run that line of code.
1. Highlight a block of code or part of a line of code and press `Ctrl+Enter` to run the highlighted code.
1. Press `Ctrl+Shift+Enter` to run all the code in a script.
1. Press the Run button on the toolbar to run all the code in a script.
1. Use the function `source()` to run all the code in a script e.g. `source("stats19-lesson-1.R")`
<!-- (but don't create an infinite loop!) -->

**Pro tip:** Try jumping between the console and the source editor by pressing Ctl+1 and Ctl+2.

## Viewing Objects

Create new objects by typing and running the following code chunk in a new script, e.g. called `objects.R`.

```{r}
vehicle_type = c("car", "bus", "tank")
casualty_type = c("pedestrian", "cyclist", "cat")
casualty_age = seq(from = 20, to = 60, by = 20)
set.seed(1)
dark = sample(x = c(TRUE, FALSE), size = 3, replace = TRUE)
small_matrix = matrix(1:24, nrow = 12)
crashes = data.frame(vehicle_type, casualty_type, casualty_age, dark)
```

We can view the objects in a range of ways:

1. Type the name of the object into the console, e.g. `crashes` and `small_matrix`, and run that code. Scroll up to see the numbers that didn't fit on the screen.
1. Use the `head()` function to view just the first 6 rows e.g. `head(small_matrix)`
1. Bonus: use the `n` argument in the previous function call to show only the first 2 rows of `small_matrix`
1. Click on the `crashes` object in the environment tab to View it in a spreadsheet.
1. Run the command `View(vehicle_type)`. What just happened?

We can also get an overview of an object using a range of functions, including 
`summary()`,
`class()`,
`typeof()`,
`dim()`, and
`length()`.

You can, for example, view a summary of the `casualty_age` variable by running the following line of code:

```{r summary}
summary(casualty_age)
```

**Exercise** try these functions on each of the objects, what results do they give?

```{r summary-answers, echo=FALSE, eval=FALSE}
summary(vehicle_type)
class(vehicle_type)
typeof(vehicle_type)
dim(vehicle_type)
length(vehicle_type)
```

**Bonus**: Find out the class of the column `vehicle_type` in the data frame `crashes` with the command `class(crashes$vehicle_type)`.
Why has it changed? 
Create a new object called `crashes_char` that keeps the class of the character vectors intact by using the function `tibble::tibble()` (see [tibble.tidyverse.org](https://tibble.tidyverse.org/) and Section 4 for details).

```{r tibble1, echo=FALSE, eval=FALSE}
tibble::tibble(
  vehicle_type,
  casualty_type,
  casualty_age,
  dark
)
```

## Autocompletion

RStudio can help you write code by autocompleting it. RStudio will look for similar objects and functions after typing the first three letters of a name.

```{r autocomp, echo=FALSE}
knitr::include_graphics("https://raw.githubusercontent.com/ITSLeeds/TDS/master/courses/2day/images/autocomplete.jpg")
```

When there is more than one option you can select from the list using the mouse or arrow keys.
Within a function, you can get a list of arguments by pressing Tab.

```{r help, echo=FALSE}
knitr::include_graphics("https://raw.githubusercontent.com/ITSLeeds/TDS/master/courses/2day/images/fucntionhelp.jpg")
```

## Getting help

Every function in R has a help page. You can view the help using `?` for example `?sum`. Many packages also contain vignettes, these are long form help documents containing examples and guides. `vignette()` will show a list of all the vignettes available, or you can show a specific vignette for example `vignette(topic = "sf1", package = "sf")`.

## Commenting Code

It is good practice to use comments in your code to explain what it does. You can comment code using `#`

For example:

```{r}
# Create vector objects (a whole line comment)
x = 1:5 # a seqence of consecutive integers (inline comment)
y = c(0, 1, 3, 9, 18.1) 
```

You can comment/uncomment a whole block of text by selecting it and using `Ctrl+Shift+C`.
<!-- not sure about the next statement so commenting out (RL) -->
<!-- and you can reformat a block of code using `Ctrl+Shift  + /`.  -->

**Pro tip:** You can add a comment section using Ctrl + Shift + R


## The global environment

The Environment tab shows all the objects in your environment, this includes datasets, parameters, and any functions you have created.
By default, new objects appear in the Global Environment but you can see other environments with the drop-down menu.
For example, each package has its own environment.

Sometimes you wish to remove things from your environment, perhaps because you no longer need them or things are getting cluttered.

You can remove an object with the `rm()` function e.g. `rm(x)` or `rm(x, y)` or you can clear your whole environment with the broom button on the Environment Tab.

1. Remove the object `x` that was created in a previous section.
1. What happens when you try to print the `x` by entering it into the console?
1. Try running the following commands in order: `save.image(); rm(list = ls()); load(".RData")`. What happened?
1. How big (how many bytes) is the `.RData` file in your project's folder?
1. Tidy up by removing the `.Rdata` file with `file.remove(".Rdata")`.

## Debugging Code

All the code shown so far is reproducible.
To test RStudio's debugging features, let's write some code that fails, as illustrated in the figure below.

```{r debug, echo=FALSE, out.width="60%"}
knitr::include_graphics("https://raw.githubusercontent.com/ropensci/stats19/master/inst/rstudio-autocomplete.png")
```

1. What is the problem with the code shown in the figure?
1. Create other types of error in the code you have run (e.g. no symetrical brackets and other typos)
1. Does RStudio pick up on the errors? And what happens when you try to run buggy code?

**Always address debugging prompts to ensure your code is reproducible**

## Saving R objects

We have already seen that you can save R scripts.
You can also save individual R objects in the RDS format.

```{r}
saveRDS(crashes, "crashes.Rds")
```

We can also read back in our data.

```{r}
crashes2 = readRDS("crashes.Rds")
identical(crashes, crashes2)
```

R also supports many other formats, including CSV files, which can be created and imported with the functions `readr::read_csv()` and `readr::write_csv()` (see also the [readr](https://readr.tidyverse.org/) package).

```{r readr-write, eval=FALSE}
readr::write_csv(crashes, "crashes.csv")
crashes3 = readr::read_csv("crashes.csv")
identical(crashes3, crashes) 
```

Notice that `crashes3` and `crashes` are not identical, what has changed? Hint: read the help page associated with `?readr::write_csv`.

# Manipulating R objects

## Subsetting by index or name

Subsetting returns part of an R object. 
It can be done by providing numbers representing the positions of the elements we want (e.g. the 2^nd^ element) or with a logical vector, with values associated with `TRUE` returned. 
Two dimension object such as matrices and data frames can be subset by rows and columns.
Subsetting in base R is done with square brackets `[]` after the name of an object. **Run the following commands to practice subsetting.**

```{r, eval=FALSE}
casualty_age[2:3] # second and third casualty_age
crashes[c(1, 2), ] # first and second row of crashes
crashes$vehicle_type # returns just one column
crashes[, c("casualty_type", "casualty_age")] # first and third columns
```

```{r, eval=FALSE, echo=FALSE}
crashes[, c(1, 3)] # first and third column of crashes by positional numbers
crashes[c(2), c(3)]
crashes[c(2), c(2, 3)]
class(crashes[, c(1, 3)])
class(crashes[c(2), c(3)])
```

1. Use the `$` operator to print the `dark` column of `crashes`.
1. Subset the crashes with the `[,]` syntax so that only the first and third columns of `crashes` are returned.
1. Return the 2^nd^ row and the 3^rd^ column of the `crashes` dataset. 
1. Return the 2^nd^ row and the columns 2:3 of the `crashes` dataset. 
1. **Bonus**: what is the `class()` of the objects created by each of the previous exercises? 

## Subsetting by values

It is also possible to subset objects by the values of their elements.
This works because the `[` operator accepts logical vectors returned by queries such as 'is it less than 3?' (`x < 3` in R) and 'was it light?' (`crashes$dark == FALSE`), as demonstrated below:

```{r, eval=FALSE}
x[c(TRUE, FALSE, TRUE, FALSE, TRUE)] # 1st, 3rd, and 5th element in x
x[x == 5] # only when x == 5 (notice the use of double equals)
x[x < 3] # less than 3
x[x < 3] = 0 # assign specific elements
casualty_age[casualty_age %% 6 == 0] # just the ages that are a multiple of 6
crashes[crashes$dark == FALSE, ]
```

1. Subset the `casualty_age` object using the inequality (`<`) so that only elements less than 50 are returned.
1. Subset the `crashes` data frame so that only tanks are returned using the `==` operator.
1. **Bonus**: assign the age of all tanks to 61.

```{r, eval=FALSE, echo=FALSE}
casualty_age[casualty_age < 50] # the  casualty_age less than 50
crashes[crashes$vehicle_type == "tank", ] # rows where the name is tank
crashes$casualty_age[crashes$vehicle_type == "tank"] = 61
```

## Dealing with NAs and recoding

R objects can have a value of NA. This is how R represents missing data.

```{r, eval=FALSE}
z = c(4, 5, NA, 7)
```

NA values are common in real-world data but can cause trouble, for example

```{r, eval=FALSE}
sum(z) # result is NA
```

Some functions can be told to ignore NA values.

```{r, eval=FALSE}
sum(z, na.rm = TRUE) # result is equal to 4 + 5 + 7
```

You can find NAs using the `is.na()` function, and then remove them

```{r, eval=FALSE}
is.na(z)
z_nona = z[!is.na(z)] # note the use of the not operator !
sum(z)
```

If you remove records with NAs be warned: the average of a value excluding NAs may not be representative.

## Changing class

Sometimes you may want to change the class of an object.
This is called class coercion, and can be done with functions such as `as.logical()`, `as.numeric()` and `as.matrix()`.

1. Coerce the `vehicle_type` column of `crashes` to the class `character`.
1. Coerce the `crashes` object into a matrix. What happened to the values?
1. **Bonus:** What is the difference between the output of `summary()` on `character` and `factor` variables?

```{r, echo=FALSE, eval=FALSE}
crashes$vehicle_type = as.character(crashes$vehicle_type)
as.matrix(crashes)
```

## Recoding values

Often it is useful to 'recode' values.
In the raw STATS19 files, for example, -1 means NA.
There are many ways to recode values in R, the simplest and most mature of which is the use of factors, as shown below:

```{r}
z = c(1, 2, -1, 1, 3)
l = c(NA, "a", "b", "c") # labels in ascending order
z_factor = factor(z, labels = l)
z_charcter = as.character(z_factor)
z_charcter
```

1. Recode `z` to Slight, Serious and Fatal for 1:3 respectively.
1. Bonus: read the help file at `?dplyr::case_when` and try to recode the values using this function.

## Now you are ready to use R

**Bonus: reproduce the following plot**

```{r smile, out.width="30%", fig.align="center"}
eyes = c(2.3, 4, 3.7, 4)
eyes = matrix(eyes, ncol = 2, byrow = T)
mouth = c(2, 2, 2.5, 1.3, 3, 1, 3.5, 1.3, 4, 2)
mouth = matrix(mouth, ncol = 2, byrow = T)
plot(eyes, type = "p", main = "RRR!", cex = 2, xlim = c(1, 5), ylim = c(0, 5))
lines(mouth, type = "l", col = "red")
```

\newpage

# R Packages

## What are packages?

R has over 15,000 packages (effectively plugins for base R), extending it in almost every direction of statistics and computing.
Packages provide additional functions, data and documentation. They are very often written by subject-matter experts and therefore tend to fit well with the workflow of the analyst in that particular specialism.
There are two main stages to using a package: installing it and loading it.
A third stage is updating it, this is also important.

<!-- installing it... -->
Install new packages from [The Comprehensive R Archive Network](https://cran.r-project.org/) with the command `install.packages()` (or `remotes::install_github()` to install from GitHub).
Update packages with the command `update.package()` or in Tools > Check for Package Updates in RStudio.
You only need to install a package once.
<!-- **Note: avoid `install.packages()` within a script** -->
<!-- Packages only need to be installed once. -->
<!-- You can use `remotes::install_cran()` or `remotes::install_github()` to only install a package if it is not yet installed and up-to-date (note: you only need to use one of these): -->

```{r, eval=FALSE}
install.packages("sf")
# remotes::install_github("r-spatial/sf")
```

<!-- now talk about loading packages -->
Installed packages are loaded with the command `library()`.
Usually, the package will load silently.
In some cases the package will provide a message, as illustrated below.

```{r}
library(sf)
```

To use a function in a package without first loading the package, use double colons, as shown below (this calls the `tibble()` function from the `tibble` package).

```{r tibble2, eval=FALSE}
crashes_tibble = tibble::tibble(
  vehicle_type,
  casualty_type,
  casualty_age,
  dark
)
```

1. Take a look in the Packages tab in the Files pane in RStudio (bottom right by default).
1. What version of the `stats19` package is installed on your computer?
1. Run the command `update.packages()`. What happens? Why?

## ggplot2

Let's take a look at a particular package.
`ggplot2` is a generic plotting package that is part of the ['tidyverse'](https://www.tidyverse.org/) meta-package, which is an "opinionated collection of R packages designed for data science". 
All packages in the tidyverse "share an underlying design philosophy, grammar, and data structures". 
`ggplot2` is flexible, popular, and has dozens of add-on packages which build on it, such as `gganimate`.
To plot non-spatial data, it works as follows (see figure below, left for result):

```{r, message=FALSE, out.width="40%", eval=FALSE}
library(ggplot2)
ggplot(crashes) + geom_point(aes(x = casualty_type, y = casualty_age))
```

Note that the `+` operator adds layers onto one another.

1. Install a package that build on `ggplot2` that begins with with `gg`. Hint: enter `install.packages(gg)` and hit Tab when your cursor is between the `g` and the `)`.
1. Open a help page in the newly installed package with the `?package_name::function()` syntax.
1. Attach the package.
1. **Bonus:** try using functionality from the new 'gg' package building on the example above to create plots like those shown below (hint: the right plot below uses the economist theme from the `ggthemes` package, try other themes).

```{r gg-extend, echo=FALSE, message=FALSE, eval=FALSE}
library(ggplot2)
# install.packages("ggthemes")
g1 = ggplot(crashes) + geom_point(aes(x = casualty_type, y = casualty_age)) 
g2 = ggplot(crashes) + geom_point(aes(x = casualty_type, y = casualty_age)) +
  ggthemes::theme_economist()
g3 = cowplot::plot_grid(g1, g2)
ggsave(filename = "inst/ggtheme-plot.png", width = 8, height = 2, dpi = 80)
```

```{r gg2, echo=FALSE, out.width="80%", fig.align="center"}
library(ggplot2)
knitr::include_graphics("https://raw.githubusercontent.com/ropensci/stats19/b4c40ad4c134853007493a9eac116b00acd4ec5a/inst/ggtheme-plot.png")
```

## dplyr and pipes

Another useful package in the tidyverse is `dplyr`.
It provides functions for manipulating data frames and using the pipe operator ` %>% `. 
The pipe puts the output of one command into the first argument of the next, as shown below (note the results are the same):

```{r}
library(dplyr)
class(crashes)       
crashes %>% class()
```

Useful `dplyr` functions are demonstrated below.

```{r, eval=FALSE}
crashes %>%
  filter(casualty_age > 50) # filter rows
crashes %>%
  select(casualty_type) # select just one column
crashes %>%
  group_by(dark) %>% 
  summarise(mean_age = mean(casualty_age))
```

1. Use `dplyr` to filter row in which `casualty_age` is less than 18, and then 28.
1. Use the `arrange` function to sort the `crashes` object in descending order of age (hint: see the `?arrange` help page).
1. Read the help page of `dplyr::mutate()`. What does the function do?
1. Use the mutate function to create a new variable, `birth_year`, in the `crashes` data.frame which is defined as the current year minus their age.
1. **Bonus:** Use the ` %>% ` operator to filter the output from the previous exercise so that only observations with `birth_year` after 1969 are returned.

```{r dplyr, eval=FALSE, echo=FALSE}
# answers
crashes %>% 
  arrange(desc(casualty_age))
crashes %>% filter(casualty_age > 21)
crashes %>% 
  mutate(birth_year = 2019 - casualty_age) %>% 
  filter(birth_year > 1969)
```

# Temporal data

For the analysis and manipulation of temporal data we will first load the R package `lubridate`:

```{r, message=FALSE}
library(lubridate)
```

The simplest example of a Date object that we can analyze is just the current date, i.e.

```{r}
today()
```

We can manipulate this object using several `lubridate` functions to extract the current day, month, year, weekday and so on...

```{r, eval=FALSE}
x = today()
day(x)
month(x)
year(x)
weekdays(x)
```

Exercises: 

1. Look at the help page of the function `month` to see how it is possible to extract the current month as character vector 
1. Look at other functions in lubridate to extract the current weekday as a number, the week of year and the day of the year

Date variables are often stored simply as a character vectors.
This is a problem, since R is not always smart enough to distinguish between character vectors representing Dates.
`lubridate` provides functions that can translate a wide range of date encodings such as `ymd()`, which extracts the Year Month and Day from a character string, as demonstrated below.

```{r, eval=FALSE}
as.Date("2019-10-17") # works
as.Date("2019 10 17") # fails
ymd("2019 10 17") # works
dmy("17/10/2019") # works
```

Import function such as `read_csv` try to recognize the Date variables.
Sometimes this fails.
You can manually create Date objects, as shown below.

```{r}
x = c("2009-01-01", "2009-02-02", "2009-03-03")
x_date = ymd(x)
x_date
```

Exercises: 

1. Extract the day, the year-day, the month and the weekday (as a non-abbreviated character vector) of each element of `x_date`. 
1. Convert `"09/09/93"` into a date object and extract its weekday. 
1. **Bonus:** Read the help page of `as.Date` and `strptime` for further details on base R functions for dates. 
1. **Bonus:** Read the Chapter 16 of [R for Data Science book](https://r4ds.had.co.nz/dates-and-times.html) for further details on `lubridate` package. 

```{r, echo=FALSE, eval=FALSE}
# 1. Extract the day, the year-day, the month and the weekday (as a non-abbreviated character vector) of each element of `x_date`. 
day(x_date)
yday(x_date)
month(x_date)
weekdays(x_date, abbreviate = FALSE)
# 1. Modify the previous example to parse the following character string: `"09/09/1993"` and extract its weekday. 
weekdays(dmy("09/09/93"))
wday(dmy("09/09/93"))
```

We can use Dates also for subsetting events in a dataframe. For example, if we define `x_date` as before and add it to the `crash` dataset, i.e.

```{r}
crashes$casualty_day = x_date
```

then we can subset events using Dates. For example

```{r}
filter(crashes, day(casualty_day) < 7) # the events that ocurred in the first week of the month
filter(crashes, weekdays(casualty_day) == "Monday") # the events occurred on monday
```

Exercises: 

1. Select only the events (rows in `crashes`) that occurred in January
1. Select only the events that ocurred in an odd year-day 
1. Select only the events that ocurred in a leap-year (HINT: check the function `leap_year`)
1. Select only the events that ocurred during the weekend or in June
1. Select only the events that ocurred during the weekend and in June
1. Count how many events ocurred during each day of the week. 

Now we'll take a look at the time components of a Date. Using the function `hms` (acronym for Hour Minutes Seconds) and its subfunctions such as `hm` or `ms`, we can parse a character vector representing several times as an Hour object (which is tecnically called a Period object). 

```{r}
x = c("18:23:35", "00:00:01", "12:34:56")
x_hour = hms(x)
x_hour
```

We can manipulate these objects using several `lubridate` functions to extract the hour component, the minutes and so on:

```{r}
hour(x_hour)
minute(x_hour)
second(x_hour)
```

If the Hour data do not specify the seconds, then we just have to use a subfunction of `hms`, namely `hm`, and everything works as before. 

```{r}
x = c("18:23", "00:00", "12:34")
(x_hour = hm(x))
```

We can use Hour data also for subsetting events, like we did for Dates. Let's add a new column to crashes data, 

```{r}
crashes$casualty_hms = hms(c("18:23:35", "00:00:01", "12:34:56"))
crashes$casualty_hour = hour(crashes$casualty_hms)
```

Exercises: 

1. Filter only the events that ocurred after midday (i.e. the PM events). Hint: your answer may include `>= 12`.
1. Filter only the events that ocurred between 15:00 and 19:00
<!-- 1. Round all hours to the next hour. Hint: Look at the help page of the `round_date` function.  -->
1. **Bonus (difficult):** run the following code, which downloades data for car crashes occurred during 2017.

```{r, eval=FALSE}
library(stats19)
crashes_2017 = stats19::get_stats19(year = 2017, type = "ac")
crashes_2017
```

Extract the weekday from the variable called `date`.
How many crashes happened on Monday?

**Advanced challenge:** calculate how many crashes occurred for each day of the week. Then plot it with ggplot2. Repeat the same exercises extracting the hour of the car accident from the variable called time. How would you combine the two informations in a single plot? 

```{r, eval=FALSE, echo=FALSE}
# solutions
crashes %>% filter(casualty_hour >= 12)
crashes %>% filter(casualty_hour > 15 & casualty_hour < 19)

crashes_2017 %>% 
  mutate(my_weekdays = weekdays(date)) %>%
  filter(my_weekdays == "Monday") %>% 
  nrow()
crashes_2017 %>% 
  mutate(my_weekdays = weekdays(date)) %>%
  filter(my_weekdays == "Friday") %>% 
  nrow()

crashes_2017 %>% 
  mutate(my_weekdays = weekdays(date)) %>% 
  group_by(my_weekdays) %>% 
  summarize(n = n()) %>% 
  ggplot() + 
  geom_col(aes(x = my_weekdays, y = n))

crashes_2017 %>% 
  mutate(my_hours = hour(hm(time))) %>% 
  group_by(my_hours) %>% 
  summarize(n = n()) %>% 
  ggplot() + 
  geom_col(aes(x = my_hours, y = n))

crashes_2017 %>% 
  mutate(my_weekdays = weekdays(date), my_hours = hour(hm(time))) %>% 
  group_by(my_weekdays, my_hours) %>% 
  summarise(n = n()) %>% 
  ggplot() + 
  geom_line(aes(x = my_hours, y = n, col = my_weekdays), size = 1.05)
# the legend needs some reordering
```

# Spatial data

## sf objects

All road crashes happen somewhere and, in the UK at least, all collisions recorded by the police are given geographic coordinates, something that can help prioritise interventions to save lives by intervening in and around 'crash hotspots'.
R has strong geographic data capabilities, with the `sf` package provides a generic class for spatial vector data: points, lines and polygons, are represented in `sf` objects as a special 'geometry column', typically called 'geom' or 'geometry', extending the data frame class we've already seen in `crashes`.

Create an `sf` data frame called `crashes_sf` as follows:

```{r crashes-sf, fig.height=2, fig.width=3}
library(sf) # load the sf package for working with spatial data
crashes_sf = crashes # create copy of crashes dataset
crashes_sf$longitude = c(-1.3, -1.2, -1.1)
crashes_sf$latitude = c(50.7, 50.7, 50.68)
crashes_sf = st_as_sf(crashes_sf, coords = c("longitude", "latitude"), crs = 4326)
# plot(crashes_sf[1:4]) # basic plot
# mapview::mapview(crashes_sf) # for interactive map
```

1. Plot only the geometry column of `crashes_sf` (hint: the solution may contain `$geometry`). If the result is like the figure below, congratulations, it worked!).
1. Plot `crashes_sf`, only showing the age variable.
1. Plot the 2^nd^ and 3^rd^ crashes, showing which happened in the dark.
1. **Bonus**: How far are the points apart (hint: `sf` functions begin with `st_`)?
1. **Bonus**: Near which settlement did the tank runover the cat?

```{r crashes-sf-ex, echo=FALSE, out.width="30%", fig.show='hold'}
plot(crashes_sf$geometry)
plot(crashes_sf["casualty_age"])
plot(crashes_sf[2:3, "dark"])
# st_distance(crashes_sf)
# Bembridge

# # updload geographic crash data
# write_sf(crashes_sf, "crashes_sf.geojson")
# piggyback::pb_upload("crashes_sf.geojson")
```

## Reading and writing spatial data

You can read and write spatial data with `read_sf()` and `write_sf()`, as shown below (see `?read_sf`).

```{r, eval=FALSE}
write_sf(zones, "zones.geojson") # save geojson file
write_sf(zones, "zmapinfo", driver = "MapInfo file")
read_sf("zmapinfo") # read in mapinfo file
```

See [Chapter 6](https://geocompr.robinlovelace.net/read-write.html) of Geocomputation with R for further information.

## sf polygons

Note: the code beyond this point is not evaluated in the vignette:

```{r}
knitr::opts_chunk$set(eval = FALSE)
```


`sf` objects can also represent administrative zones.
This is illustrated below with reference to `zones`, a spatial object representing the Isle of Wight, that we will download using the `pct` package (note: the `[1:9]` appended to the function selects only the first 9 columns).

```{r}
zones = pct::get_pct_zones("isle-of-wight")[1:9]
```

1. What is the class of the `zones` object?
1. What are its column names?
1. Print its first 2 rows and columns 6:8 (the result is below).

```{r, echo=FALSE}
# class(zones)
# names(zones)
zones[1:2, c(1, 5, 6, 7, 8)]
```

## Spatial subsetting and sf plotting

Like index and value subsetting, spatial subsetting can be done with the `[` notation.
Subset the `zones` that contain features in `crashes_sf` as follows:

```{r, message=FALSE}
zones_containing_crashes = zones[crashes_sf, ]
```

To plot a new layer on top of an existing `sf` plot, use the `add = TRUE` argument.
Remember to plot only the `geometry` column of objects to avoid multiple maps.
Colours can be set with the `col` argument.

1. Plot the geometry of the zones, with the zones containing crashes overlaid on top in red.
1. Plot the zone containing the 2^nd^ crash in blue.
1. **Bonus:** plot all zones that intersect with a zone containing crashes, with the actual crash points plotted in black.

```{r sp-ex, echo=FALSE, out.width="33%", fig.show='hold', message=FALSE, warning=FALSE}
plot(zones$geometry)
plot(zones_containing_crashes$geometry, col = "red", add = TRUE)
plot(zones$geometry)
plot(zones[crashes_sf[2, ], ], col = "blue", add = TRUE)
plot(zones$geometry)
plot(zones[zones_containing_crashes, ], col = "yellow", add = TRUE)
plot(crashes_sf$geometry, pch = 20, add = TRUE)
```

## Geographic joins

Geographic joins involve assigning values from one object to a new column in another, based on the geographic relationship between them.
With `sf` objects it works as follows:

```{r, message=FALSE}
zones_joined = st_join(zones[1], crashes_sf)
```

1. Plot the `casualty_age` variable of the new `zones_joined` object (see the figure below to verify the result).
1. How many zones are returned in the previous command? 
1. Select only the `geo_code` column from the `zones` and the `dark` column from `crashes_sf` and use the `left = FALSE` argument to return only zones in which crashes occured. Plot the result.

See [Chapter 4](https://geocompr.robinlovelace.net/spatial-operations.html#spatial-joining) of Geocomputation with R [@lovelace_geocomputation_2019] for further information on geographic joins.

```{r joinf, echo=FALSE, out.width="40%", fig.show='hold', message=FALSE}
plot(zones_joined["casualty_age"])
zjd = st_join(zones[1], crashes_sf["dark"], left = FALSE)
plot(zjd)
```


## CRSs

Get and set Coordinate Reference Systems (CRSs) with the command `st_crs()`.
Transform CRSs with the command `st_transform()`, as demonstrated in the code chunk below, which converts the 'lon/lat' geographic CRS of `crashes_sf` into the projected CRS of the British National Grid:

```{r crs1}
crashes_osgb = st_transform(crashes_sf, 27700)
```

1. Try to subset the zones with the `crashes_osgb`. What does the error message say?
1. Create `zones_osgb` by transforming the `zones` object.
1. **Bonus:** use `st_crs()` to find out the units measurement of the British National Grid?

For more information on CRSs see [Chapter 6](https://geocompr.robinlovelace.net/reproj-geo-data.html) of Geocompuation with R.

## Buffers

Buffers are polygons surrounding geometries of a (usually) fixed distance.
Currently buffer operations in R only work on objects with projected CRSs.

1. Find out and read the help page of `sf`'s buffer function.
1. Create an object called `crashes_1km_buffer` representing the area within 1 km of the crashes.
1. **Bonus:** try creating buffers on the geographic version of the `crashes_sf` object. What happens?

## Attribute operations on sf objects

Because `sf` objects are `data.frame`s, we can do non-spatial operations on them.
Try the following attribute operations on the `zones` data.

```{r}
# load example dataset if it doesn't already exist
zones = pct::get_pct_zones("isle-of-wight")
sel = zones$all > 3000  # create a subsetting object
zones_large = zones[sel, ] # subset areas with a popualtion over 100,000
zones_2 = zones[zones$geo_name == "Isle of Wight 002",] # subset based on 'equality' query
zones_first_and_third_column = zones[c(1, 3)]
zones_just_all = zones["all"]
```


1. Practice subsetting techniques you have learned on the `sf data.frame` object `zones`:
     1. Create an object called `zones_small` which contains only regions with less than 3000 people in the `all` column
     1. Create a selection object called `sel_high_car` which is `TRUE` for regions with above median numbers of people who travel by car and `FALSE` otherwise
     1. Create an object called `zones_foot` which contains only the foot attribute from `zones`
     1. Bonus: plot `zones_foot` using the function `plot` to show where walking is a popular mode of travel to work
     1. Bonus: bulding on your answers to previous questions, use `filter()` from the `dplyr` package to subset small regions where car use is high. 
1. Bonus: What is the population density of each region (hint: you may need to use the functions `st_area()`, `as.numeric()` and use the 'all' column)?
1. Bonus: Which zone has the highest percentage of people who cycle?

```{r, echo=FALSE, eval=FALSE}
# 1. Practice subsetting techniques you have learned on the `sf data.frame` object `zones`:
#      1. Create an object called `zones_small` which contains only regions with less than 3000 people in the `all` column
# in base R
zones_small = zones[zones$all < 3000, ]
# with dplyr
zones_small = zones %>% 
  filter(all < 3000)
#      1. Create a selection object called `sel_high_car` which is `TRUE` for regions with above median numbers of people who travel by car and `FALSE` otherwise
median_car = median(zones$car_driver)
sel_high_car = zones$car_driver > median_car 
#      1. How many regions have the number '1' in the column 'geo_name'? What percentage of the regions in the Isle of Wight is this?
sel_region_name_contains_1 = grepl("1", x = zones$geo_name)
sum(sel_region_name_contains_1) / nrow(zones)
#      1. Create an object called `zones_foot` which contains only the foot attribute from `zones`
# using base R
zones_foot = zones["foot"]
# dplyr
zones_foot = zones %>% 
  select(foot)
#      1. Bonus: plot the result to show where walking is a popular mode of travel to work
plot(zones_foot)
#      1. Bonus: bulding on your answers to previous questions, use `filter()` from the `dplyr` package to subset small regions where high car use is high
zones_small_car_high = zones %>% 
  filter(all < 3000, car_driver > median_car)
# 1. Bonus: What is the population density of each region (hint: you may need to use the functions `st_area()`, `as.numeric()` and use the 'all' column)?
zones$area_km2 = as.numeric(st_area(zones)) /1000000
zones$population_density = zones$all / zones$area_km2
plot(zones["population_density"])
# in dplyr
zones_density = zones %>% 
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1000000) %>% 
  mutate(population_density = all / area_km2)
plot(zones_density %>% select(population_density))
# 1. Bonus: Which zone has the highest percentage who cycle?
zones %>% 
  mutate(pcycle = bicycle / all) %>% 
  top_n(n = 1, wt = pcycle)
# 1. Bonus: Find the proportion of people who drive to work (`car_driver`) in areas in which more than 500 people walk to work
zones %>% 
  group_by(foot > 500) %>% 
  summarise(mean_car = sum(car_driver) / sum(all) )
```

## Matching roads to crashes

I think you forgot something here. For example we could introduce `st_nearest_feature`? Or counting using `st_within` and `st_buffer`. 

# Visualising spatial datasets

So far we have used the `plot()` function to make maps.
That's fine for basic visualisation, but for publication-quality maps, we recommend using `tmap` (see Chapter 8 of Geocomputation with R for reasons and alternatives).
Load the package as follows:

```{r}
library(tmap)
tmap_mode("plot")
```

1. Create the following plots using `plot()` and `tm_shape() + tm_polygons()` functions (note: the third figure relies on setting `tmap_mode("view")`.
1. Add an additional layer to the interactive map showing the location of crashes, using marker and dot symbols.
1. Bonus: Change the default basemap (hint: you may need to search in the package documentation or online for the solution).

```{r plot3, fig.show='hold', out.width="33%", echo=FALSE}
plot(zones[c("all", "bicycle")])
# tm_shape(zones) + 
#   tm_polygons(c("all", "bicycle"))
# tmap_mode("view")
# m = tm_shape(zones_joined) + 
#   tm_polygons(c("casualty_type")) +
#   tm_scale_bar()
# m
# knitr::include_graphics("tmap-zones-interactive.png")
# piggyback::pb_upload("zones_joined.Rds")
# create bug report:
# See https://github.com/mtennekes/tmap/issues/551
# piggyback::pb_download_url("zones_joined.Rds")
# "https://github.com/ropensci/stats19/releases/download/1.3.0/zones_joined.Rds"
# library(tmap)
# u = "https://github.com/ropensci/stats19/releases/download/1.3.0/zones_joined.Rds"
# zones_joined = readRDS(url(u))
# qtm(zones_joined)
```

# Analysing point data from stats19

Based on the saying "don't run before you can walk", we've learned the vital foundations of R before tackling a real dataset.
Temporal and spatial attributes are key to road crash data, hence the emphasis on `lubridate` and `sf`.
Visualisation is key to understanding and policy influence, which is where `tmap` comes in.
With these solid foundations, plus knowledge of how to ask for help (read R's internal help functions, ask colleagues, create new comments on online forums/GitHub, generally in that order of priority), you are ready to test the methods on some real data.

Before doing so, take a read of the `stats19` vignette, which can be launched as follows:

```{r, eval=FALSE}
vignette(package = "stats19") # view all vignettes available on stats19
vignette("stats19") # view the introductory vignette
```

This should now be sufficient to tackle the following exercises:

1. Download and plot all crashes reported in Great Britain in 2018 (hint: see [the stats19 vignette](https://docs.ropensci.org/stats19/articles/stats19.html))
1. Find the function in the `stats19` package that converts a `data.frame` object into an `sf` data frame. Use this function to convert the road crashes into an `sf` object, called `crashes_sf`, for example.
1. Filter crashes that happened in the Isle of Wight based on attribute data (hint: the relevant column contains the word `local`)
1. Filter crashes happened in the Isle of Wight using geographic subsetting (hint: remember `st_crs()`?)
1. **Bonus:** Which type of subsetting yielded more results and why? 
1. **Bonus:** how many crashes happened in each zone?
1. Create a new column called `month` in the crash data using the function `lubridate::month()` and the `date` column.
1. Create an object called `a_zones_may` representing all the crashes that happened in the Isle of Wight in the month of May
1. Bonus: Calculate the average (`mean`) speed limit associated with each crash that happened in May across the zones of the Isle of Wight (the result is shown in the map)


```{r, echo=FALSE, results='hide', message=FALSE, eval=FALSE}
library(stats19)
library(dplyr)
library(sf)
a = get_stats19(2018, "ac")
asf = format_sf(a)
a_zones = asf %>% 
  filter(local_authority_district == "Isle of Wight")
nrow(a_zones)
zones = pct::get_pct_zones(region = "isle-of-wight")
zones_osbg = st_transform(zones, 27700)
a_zones_sf = a_zones[zones_osbg, ]
nrow(a_zones_sf)
# mapview::mapview(zones) +
#   mapview::mapview(a_zones)
class(a$date)
class(a$time)
a_zones$month = lubridate::month(a_zones$date)
a_zones_may = a_zones %>% 
  filter(month == 5)
a_agg = aggregate(a_zones_may["speed_limit"], zones_osbg, mean)
plot(a_agg)
class(a$date)
```

# Analysing crash data on road networks

Road network data can be accessed from a range of sources, including OpenStreetMap (OSM) and Ordnance Survey.
We will use some OSM data from the Ilse of Wight, which can be loaded as follows:

```{r}
u = "https://github.com/ropensci/stats19/releases/download/1.1.0/roads_key.Rds"
roads_wgs = readRDS(url(u))
roads = roads_wgs %>% st_transform(crs = 27700)
```

You should already have road crashes for the Isle of Wight from the previous stage.
If not, load crash data as follows:

```{r}
u = "https://github.com/ropensci/stats19/releases/download/1.1.0/car_accidents_2017_iow.Rds"
crashes_iow = readRDS(url(u))
```

1. Plot the roads with the crashes overlaid.
2. Create a buffer around the roads with a distance of 200 m.
3. How many crashes fall outside the buffered roads?
3. Bonus: Use the `aggregate()` function to identify how many crashes happened per segment and plot the result (hint: see `?aggregate.sf` and take a read of Section [4.2.5](https://geocompr.robinlovelace.net/spatial-operations.html#spatial-aggr) of Geocomputation with R) with `tmap` and plot the crashes that happened outside the road buffers on top.

```{r, echo=FALSE, out.width="49%", fig.show='hold', message=FALSE}
plot(roads$geometry)
plot(crashes_iow["accident_severity"], add = TRUE)
roads_buffer = st_buffer(roads, 200, endCapStyle = "FLAT")
crashes_outside_roads = crashes_iow[roads_buffer, , op = sf::st_disjoint]
roads_agg = aggregate(crashes_iow[1], by = roads_buffer, FUN = length)
# plot(roads_agg, border = NA, main = "")
names(roads_agg)[1] = "N. Crashes"
# tmap_mode("plot")
# tm_shape(roads_agg) + tm_fill("N. Crashes") +
#   tm_shape(crashes_outside_roads) + tm_dots(col = "blue")
```

\newpage

# Bonus exercises

Identify a region and zonal units of interest from http://geoportal.statistics.gov.uk/ or from the object `police_boundaries` in the `stats19` package.

1. Read them into R as an `sf` object
1. Create a map showing the number of crashes in each zone
1. Identify the average speed limit associated with crashes in each zone
1. Identify an interesting question you can ask to the data and use exploratory data analysis to find answers
1. Check another [related project](https://github.com/agila5/leeds_seminar) for further information on smoothing techniques of counts on a linear network. 
<!-- 1. Take a look at the code in [the file iow_example.R in the inst directory in the stats19 repo](https://github.com/ropensci/stats19/blob/master/inst/iow_example.R). Run it to create smoothed estimates of crash frequency on the road network (see [code in the  GitHub repo](https://github.com/agila5/leeds_seminar/blob/master/examples.R) for further information on these preliminary methods). -->

```{r final-plot, echo=FALSE, out.width="100%"}
# knitr::include_graphics("final-figure.png")
```


# References
---
title: "Researching vehicles involved in collisions with STATS19 data"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Researching vehicles involved in collisions with STATS19 data}
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

```{r setup, message=FALSE}
library(stats19)
library(dplyr)
```

# Vehicle level variables in the STATS19 datasets

Of the three dataset types in STATS19, the vehicle tables are perhaps the most revealing yet under-explored.
They look like this:

```{r}
v = get_stats19(year = 2018, type = "vehicle")
names(v)
v
```

We will categorise the vehicle types to simplify subsequent results:

```{r}
v = v %>% mutate(vehicle_type2 = case_when(
  grepl(pattern = "motorcycle", vehicle_type, ignore.case = TRUE) ~ "Motorbike",
  grepl(pattern = "Car", vehicle_type, ignore.case = TRUE) ~ "Car",
  grepl(pattern = "Bus", vehicle_type, ignore.case = TRUE) ~ "Bus",
  grepl(pattern = "cycle", vehicle_type, ignore.case = TRUE) ~ "Cycle",
  # grepl(pattern = "Van", vehicle_type, ignore.case = TRUE) ~ "Van",
  grepl(pattern = "Goods", vehicle_type, ignore.case = TRUE) ~ "Goods",
  
  TRUE ~ "Other"
))
# barplot(table(v$vehicle_type2))
```


All of these variables are of potential interest to road safety researchers.
Let's take a look at summaries of a few of them:

```{r}
table(v$vehicle_type2)
summary(v$age_of_driver)
summary(v$engine_capacity_cc)
table(v$propulsion_code)
summary(v$age_of_vehicle)
```

The output shows vehicle type (a wide range of vehicles are represented), age of driver (with young and elderly drivers often seen as more risky), engine capacity and populsion (related to vehicle type and size) and age of vehicle.
In addition to these factors appearing in prior road safety research and debate, they are also things that policy makers can influence, e.g by:

- Encouraging modal shift away more dangerous modes and towards safer modes
- Incentivising people in particular risk categories to use safer modes
- Encouraging use of certain (safer) kinds of vehicle, e.g. with tax policies

# Relationships between vehicle type and crash severity

To explore the relationship between vehicles and crash severity, we must first join on the 'accidents' table:

```{r}
a = get_stats19(year = 2018, type = "accidents")
va = dplyr::inner_join(v, a)
```

Now we have additional variables available to us:

```{r}
dim(v)
dim(va)
names(va)
```

Let's see how crash severity relates to the variables of interest mentioned above:

```{r, out.width="100%"}
xtabs(~vehicle_type2 + accident_severity, data = va) %>% prop.table()
xtabs(~vehicle_type2 + accident_severity, data = va) %>% prop.table() %>% plot()
```

As expected, crashes involving large vehicles such as buses and good vehicles tend to be more serious (involve proportionally more deaths) than crashes involving smaller vehicles.

To focus only on cars, we can filter the `va` table as follows:

```{r}
vac = va %>% filter(vehicle_type2 == "Car")
```

The best proxy we have for car type in the open STATS19 data (there are non-open versions of the data with additional columns) is engine capacity, measured in cubic centimetres (cc).
The distribution of engine cc's in the cars dataset created above is shown below.

```{r}
summary(vac$engine_capacity_cc)
```

The output shows that there are some impossible values in the data, likely due to recording error. 
Very few cars have an engine capacity above 5 litres (5000 cc) and we can be confident that none have an engine capacity below 300 cc.
We'll identify these records and remove them as follows:

```{r, echo=FALSE, eval=FALSE}
library(tidyverse)

max_engine_size = 5000
min_engine_size = 300

sel_too_big = vac$engine_capacity_cc > max_engine_size
sel_too_small = vac$engine_capacity_cc < min_engine_size
sum(sel_too_big) / nrow(vac)
sum(sel_too_small) / nrow(vac)
vac$engine_capacity_cc[sel_too_big | sel_too_small] = NA
```

We have set the anomolous vehicle size data to NA meaning it will not be used in the subsequent analysis.

```{r, eval=FALSE, echo=FALSE}

vac = vac %>% filter(val_size)
vac %>% 
  mutate(age = formatC(age_band_of_driver, digits = 2, flag = "0")) %>% 
  ggplot() +
  geom_violin(aes(age, engine_capacity_cc)) 

vac$sev_factor = factor(vac$accident_severity, labels = 3:1)
vac$sev_numeric = vac$sev_factor %>% as.character() %>%  as.numeric()
summary(vac$sev_factor)
summary(vac$sev_numeric)

m = lm(sev_numeric ~ engine_capacity_cc + age_of_driver + speed_limit, data = vac)
summary(m)
```



```{r, echo=FALSE, eval=FALSE}
table_vehicle_type = xtabs(cbind(accident_severity, vehicle_type) ~ accident_severity, data = va)
group_totals = va %>% 
  group_by(accident_severity) %>% 
  summarise(n = n())

# fails
fit = glm(data = va, factor(accident_severity) ~
            engine_capacity_cc +
            age_of_driver +
            engine_capacity_cc +
            factor(propulsion_code) +
            age_of_vehicle
          )
# works but result does not have probabilities
?nnet::multinom
mod = nnet::multinom(formula = accident_severity ~
            engine_capacity_cc +
            age_of_driver +
            engine_capacity_cc +
            propulsion_code +
            age_of_vehicle, data = va)
mod
summary(mod)
mod$nunits
class(mod$fitted.values)
dim(mod$fitted.values)
colnames(mod$fitted.values)
summary(mod$fitted.values) # result!

probs = as.data.frame(mod$fitted.values)
head(probs)
head(rowSums(probs))
colSums(probs) / group_totals$n
nrow(probs)
nrow(va)

# install.packages("mlogit")
install.packages("AER")
vignette(package = "mlogit")
vignette("c2.formula.data")
library(mlogit)
data("TravelMode", package = "AER")
head(TravelMode)
?TravelMode
summary(TravelMode$choice)
summary(TravelMode$mode)
TM = mlogit.data(TravelMode, choice = "choice", shape = "long",
                 alt.levels = c("air", "train", "bus", "car"))

vamld = mlogit.data(va, choice = "accident_severity", alt.levels = c("Slight", "Serious", "Fatal"), shape = "wide")

mlogit(accident_severity ~ speed_limit | 0, vamld[1:999, ])

vamld = mlogit::mlogit.data(va, choice = "accident_severity", shape = "wide")
head(vamld)


# fails
# vm = mlogit(accident_severity ~ engine_capacity_cc + speed_limit | 0, vamld[1:999, ])
# apply(fitted(vm, outcome = FALSE), 2, mean)
# vm
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{check_input_file}
\alias{check_input_file}
\title{Local helper to be reused.}
\usage{
check_input_file(filename = NULL, type = NULL, data_dir = NULL, year = NULL)
}
\arguments{
\item{filename}{Character string of the filename of the .csv to read, if this is given, type and
years determine whether there is a target to read, otherwise disk scan would be needed.}

\item{type}{The type of file to be downloaded (e.g. 'Accidents', 'Casualties' or
'Vehicles'). Not case sensitive and searches using regular expressions ('acc' will work).}

\item{data_dir}{Where sets of downloaded data would be found.}

\item{year}{Single year for which data are to be read}
}
\description{
Local helper to be reused.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{set_data_directory}
\alias{set_data_directory}
\title{Set data download dir}
\usage{
set_data_directory(data_path)
}
\arguments{
\item{data_path}{valid existing path to save downloaded files in.}
}
\description{
Handy function to manage \code{stats19} package underlying environment
variable. If run interactively it makes sure user does not change
directory by mistatke.
}
\examples{
# set_data_directory("MY_PATH")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_url}
\alias{get_url}
\title{Convert file names to urls}
\usage{
get_url(
  file_name = "",
  domain = "https://data.dft.gov.uk",
  directory = "road-accidents-safety-data"
)
}
\arguments{
\item{file_name}{Optional file name to add to the url returned (empty by default)}

\item{domain}{The domain from where the data will be downloaded}

\item{directory}{The subdirectory of the url}
}
\description{
Convert file names to urls
}
\details{
This function returns urls that allow data to be downloaded from the pages:

https://data.dft.gov.uk/road-accidents-safety-data/RoadSafetyData_2015.zip

Last updated: October 2020.
Files available from the s3 url in the default \code{domain} argument.
}
\examples{
# get_url(find_file_name(1985))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_data_directory}
\alias{get_data_directory}
\title{Get data download dir}
\usage{
get_data_directory()
}
\description{
Get data download dir
}
\examples{
# get_data_directory()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vehicles_sample}
\alias{vehicles_sample}
\alias{vehicles_sample_raw}
\title{Sample of stats19 data (2017 vehicles)}
\format{
A data frame
}
\description{
Sample of stats19 data (2017 vehicles)
}
\note{
These were generated using the script in the
\code{data-raw} directory (\code{misc.Rmd} file).
}
\examples{
\donttest{
nrow(vehicles_sample_raw)
vehicles_sample_raw
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{find_file_name}
\alias{find_file_name}
\title{Find file names within stats19::file_names.}
\usage{
find_file_name(years = NULL, type = NULL)
}
\arguments{
\item{years}{Year for which data are to be found}

\item{type}{One of 'Accidents', 'Casualties', 'Vehicles'; defaults to 'Accidents', ignores case.}
}
\description{
Currently, there are 52 file names to download/read data from.
}
\examples{
find_file_name(2016)
find_file_name(2016, type = "accident")
find_file_name(1985, type = "accident")
find_file_name(type = "cas")
find_file_name(type = "accid")
find_file_name(2016:2017) # warning when multiple years requested
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format.R
\name{format_accidents}
\alias{format_accidents}
\title{Format STATS19 'accidents' data}
\usage{
format_accidents(x)
}
\arguments{
\item{x}{Data frame created with \code{read_accidents()}}
}
\description{
Format STATS19 'accidents' data
}
\section{Details}{

This is a helper function to format raw STATS19 data
}

\examples{
\donttest{
if(curl::has_internet()) {
dl_stats19(year = 2017, type = "accident")
x = read_accidents(year = 2017, format = FALSE)
if(nrow(x) > 0) {
x[1:3, 1:12]
crashes = format_accidents(x)
crashes[1:3, 1:12]
summary(crashes$datetime)
}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ulez.R
\name{get_ULEZ}
\alias{get_ULEZ}
\title{Download DVLA-based vehicle data from the TfL API using VRM.}
\usage{
get_ULEZ(vrm)
}
\arguments{
\item{vrm}{A list of VRMs as character strings.}
}
\description{
Download DVLA-based vehicle data from the TfL API using VRM.
}
\section{Details}{

This function takes a character vector of vehicle registrations (VRMs) and returns DVLA-based vehicle data from TfL's API, included ULEZ eligibility.
It returns a data frame of those VRMs which were successfully used with the TfL API.  Vehicles are either compliant, non-compliant or exempt.  ULEZ-exempt vehicles will not have all vehicle details returned - they will simply be marked "exempt".

Be aware that the API has usage limits.  The function will therefore limit API calls to below 50 per minute - this is the maximum rate before an API key is required.
}

\examples{
\donttest{
if(curl::has_internet()) {
vrm = c("1RAC","P1RAC")
get_ULEZ(vrm = vrm)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{schema_original}
\alias{schema_original}
\title{Schema for stats19 data (UKDS)}
\format{
A data frame
}
\description{
Schema for stats19 data (UKDS)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{casualties_sample}
\alias{casualties_sample}
\alias{casualties_sample_raw}
\title{Sample of stats19 data (2017 casualties)}
\format{
A data frame
}
\description{
Sample of stats19 data (2017 casualties)
}
\note{
These were generated using the script in the
\code{data-raw} directory (\code{misc.Rmd} file).
}
\examples{
\donttest{
nrow(casualties_sample_raw)
casualties_sample_raw
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{file_names}
\alias{file_names}
\alias{file_names_old}
\title{stats19 file names for easy access}
\format{
A named list
}
\description{
URL decoded file names. Currently there are 52 file names
released by the DfT (Department for Transport) and the details include
how these were obtained and would be kept up to date.
}
\note{
These were generated using the script in the
\code{data-raw} directory (\code{misc.Rmd} file).
}
\examples{
\dontrun{
 length(file_names)
 file_names$dftRoadSafetyData_Vehicles_2017.zip
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.R
\name{get_stats19}
\alias{get_stats19}
\title{Download, read and format STATS19 data in one function.}
\usage{
get_stats19(
  year = NULL,
  type = "accident",
  data_dir = get_data_directory(),
  file_name = NULL,
  format = TRUE,
  ask = FALSE,
  silent = FALSE,
  output_format = "tibble",
  ...
)
}
\arguments{
\item{year}{A year matching file names on the STATS19
\href{https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data}{data release page}
e.g. \code{2020}}

\item{type}{One of 'Accident', 'Casualty', 'Vehicle'; defaults to 'Accident'.
Or any variation of to search the file names with such as "acc" or "accid".}

\item{data_dir}{Parent directory for all downloaded files. Defaults to \code{tempdir()}.}

\item{file_name}{The file name (DfT named) to download.}

\item{format}{Switch to return raw read from file, default is \code{TRUE}.}

\item{ask}{Should you be asked whether or not to download the files? \code{TRUE} by default.}

\item{silent}{Boolean. If \code{FALSE} (default value), display useful progress
messages on the screen.}

\item{output_format}{A string that specifies the desired output format. The
default value is \code{"tibble"}. Other possible values are \code{"data.frame"}, \code{"sf"}
and \code{"ppp"}, that, respectively, returns objects of class \code{\link{data.frame}},
\code{\link[sf:sf]{sf::sf}} and \code{\link[spatstat.geom:ppp]{spatstat.geom::ppp}}. Any other string is ignored and a tibble
output is returned. See details and examples.}

\item{...}{Other arguments be passed to \code{\link[=format_sf]{format_sf()}} or
\code{\link[=format_ppp]{format_ppp()}} functions. Read and run the examples.}
}
\description{
Download, read and format STATS19 data in one function.
}
\section{Details}{

This function uses gets STATS19 data. Behind the scenes it uses
\code{dl_stats19()} and \verb{read_*} functions, returning a
\code{tibble} (default), \code{data.frame}, \code{sf} or \code{ppp} object, depending on the
\code{output_format} parameter.
The function returns data for a specific year (e.g. \code{year = 2017})

Note: for years before 2016 the function may return data from more years than are
requested due to the nature of the files hosted at
\href{https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data}{data.gov.uk}.

As this function uses \code{dl_stats19} function, it can download many MB of data,
so ensure you have a sufficient disk space.

If \code{output_format = "data.frame"} or \code{output_format = "sf"} or \code{output_format = "ppp"} then the output data is transformed into a data.frame, sf or ppp
object using the \code{\link[=as.data.frame]{as.data.frame()}} or \code{\link[=format_sf]{format_sf()}} or \code{\link[=format_ppp]{format_ppp()}}
functions, as shown in the examples.
}

\examples{
\donttest{
if(curl::has_internet()) {
# default tibble output
x = get_stats19(2019)
class(x)
x = get_stats19(2017, silent = TRUE)

# data.frame output
x = get_stats19(2017, silent = TRUE, output_format = "data.frame")
class(x)

# Run tests only if endpoint is alive:
if(nrow(x) > 0) {

# sf output
x_sf = get_stats19(2017, silent = TRUE, output_format = "sf")

# sf output with lonlat coordinates
x_sf = get_stats19(2017, silent = TRUE, output_format = "sf", lonlat = TRUE)
sf::st_crs(x_sf)

if (requireNamespace("spatstat.core", quietly = TRUE)) {
# ppp output
x_ppp = get_stats19(2017, silent = TRUE, output_format = "ppp")

# We can use the window parameter of format_ppp function to filter only the
# events occurred in a specific area. For example we can create a new bbox
# of 5km around the city center of Leeds

leeds_window = spatstat.geom::owin(
xrange = c(425046.1, 435046.1),
yrange = c(428577.2, 438577.2)
)

leeds_ppp = get_stats19(2017, silent = TRUE, output_format = "ppp", window = leeds_window)
spatstat.geom::plot.ppp(leeds_ppp, use.marks = FALSE, clipwin = leeds_window)

# or even more fancy examples where we subset all the events occurred in a
# pre-defined polygon area

# The following example requires osmdata package
# greater_london_sf_polygon = osmdata::getbb(
# "Greater London, UK",
# format_out = "sf_polygon"
# )
# spatstat works only with planar coordinates
# greater_london_sf_polygon = sf::st_transform(greater_london_sf_polygon, 27700)
# then we extract the coordinates and create the window object.
# greater_london_polygon = sf::st_coordinates(greater_london_sf_polygon)[, c(1, 2)]
# greater_london_window = spatstat.geom::owin(poly = greater_london_polygon)

# greater_london_ppp = get_stats19(2017, output_format = "ppp", window = greater_london_window)
# spatstat.geom::plot.ppp(greater_london_ppp, use.marks = FALSE, clipwin = greater_london_window)
}
}
}
}
}
\seealso{
\code{\link[=dl_stats19]{dl_stats19()}}

\code{\link[=read_accidents]{read_accidents()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{select_file}
\alias{select_file}
\title{Interactively select from options}
\usage{
select_file(fnames)
}
\arguments{
\item{fnames}{File names to select from}
}
\description{
Interactively select from options
}
\examples{
# fnames = c("f1", "f2")
# stats19:::select_file(fnames)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_vehicles}
\alias{read_vehicles}
\title{Read in stats19 road safety data from .csv files downloaded.}
\usage{
read_vehicles(
  year = NULL,
  filename = "",
  data_dir = get_data_directory(),
  format = TRUE
)
}
\arguments{
\item{year}{Single year for which data are to be read}

\item{filename}{Character string of the filename of the .csv to read, if this is given, type and
years determine whether there is a target to read, otherwise disk scan would be needed.}

\item{data_dir}{Where sets of downloaded data would be found.}

\item{format}{Switch to return raw read from file, default is \code{TRUE}.}
}
\description{
Read in stats19 road safety data from .csv files downloaded.
}
\section{Details}{

The function returns a data frame, in which each record is a reported vehicle in the
STATS19 dataset for the data_dir and filename provided.
}

\examples{
\donttest{
if(curl::has_internet()) {
dl_stats19(year = 2019, type = "vehicle")
ve = read_vehicles(year = 2019)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{accidents_sample}
\alias{accidents_sample}
\alias{accidents_sample_raw}
\title{Sample of stats19 data (2017 accidents)}
\format{
A data frame
}
\description{
Sample of stats19 data (2017 accidents)
}
\note{
These were generated using the script in the
\code{data-raw} directory (\code{misc.Rmd} file).
}
\examples{
\donttest{
nrow(accidents_sample_raw)
accidents_sample_raw
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_casualties}
\alias{read_casualties}
\title{Read in STATS19 road safety data from .csv files downloaded.}
\usage{
read_casualties(
  year = NULL,
  filename = "",
  data_dir = get_data_directory(),
  format = TRUE
)
}
\arguments{
\item{year}{Single year for which data are to be read}

\item{filename}{Character string of the filename of the .csv to read, if this is given, type and
years determine whether there is a target to read, otherwise disk scan would be needed.}

\item{data_dir}{Where sets of downloaded data would be found.}

\item{format}{Switch to return raw read from file, default is \code{TRUE}.}
}
\description{
Read in STATS19 road safety data from .csv files downloaded.
}
\section{Details}{

The function returns a data frame, in which each record is a reported casualty
in the STATS19 dataset.
}

\examples{
\donttest{
if(curl::has_internet()) {
dl_stats19(year = 2017, type = "casualty")
casualties = read_casualties(year = 2017)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{phrase}
\alias{phrase}
\title{Generate a phrase for data download purposes}
\usage{
phrase()
}
\description{
Generate a phrase for data download purposes
}
\examples{
stats19:::phrase()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{locate_one_file}
\alias{locate_one_file}
\title{Pin down a file on disk from four parameters.}
\usage{
locate_one_file(
  filename = NULL,
  data_dir = get_data_directory(),
  year = NULL,
  type = NULL
)
}
\arguments{
\item{filename}{Character string of the filename of the .csv to read, if this
is given, type and years determine whether there is a target to read,
otherwise disk scan would be needed.}

\item{data_dir}{Where sets of downloaded data would be found.}

\item{year}{Single year for which file is to be found.}

\item{type}{One of: 'Accidents', 'Casualties', 'Vehicles'; ignores case.}
}
\value{
One of: path for one file, a message \verb{More than one file found} or error if none found.
}
\description{
Pin down a file on disk from four parameters.
}
\examples{
\donttest{
locate_one_file()
locate_one_file(filename = "Cas.csv")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format.R
\name{format_sf}
\alias{format_sf}
\title{Format convert STATS19 data into spatial (sf) object}
\usage{
format_sf(x, lonlat = FALSE)
}
\arguments{
\item{x}{Data frame created with \code{read_accidents()}}

\item{lonlat}{Should the results be returned in longitude/latitude?
By default \code{FALSE}, meaning the British National Grid (EPSG code: 27700)
is used.}
}
\description{
Format convert STATS19 data into spatial (sf) object
}
\examples{
x_sf = format_sf(accidents_sample)
sf:::plot.sf(x_sf)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format.R
\name{format_casualties}
\alias{format_casualties}
\title{Format STATS19 casualties}
\usage{
format_casualties(x)
}
\arguments{
\item{x}{Data frame created with \code{read_casualties()}}
}
\description{
Format STATS19 casualties
}
\section{Details}{

This function formats raw STATS19 data
}

\examples{
\donttest{
if(curl::has_internet()) {
dl_stats19(year = 2017, type = "casualty")
x = read_casualties(year = 2017)
casualties = format_casualties(x)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{police_boundaries}
\alias{police_boundaries}
\title{Police force boundaries in England (2016)}
\format{
An sf data frame
}
\description{
This dataset represents the 43 police forces in England and Wales.
These are described on the
\href{https://en.wikipedia.org/wiki/List_of_police_forces_of_the_United_Kingdom}{Wikipedia page}.
on UK police forces.
}
\details{
The geographic boundary data were taken from the UK government's
official geographic data portal.
See http://geoportal.statistics.gov.uk/
}
\note{
These were generated using the script in the
\code{data-raw} directory (\code{misc.Rmd} file) in the package's GitHub repo:
\href{https://github.com/ITSLeeds/stats19}{github.com/ITSLeeds/stats19}.
}
\examples{
nrow(police_boundaries)
police_boundaries[police_boundaries$pfa16nm == "West Yorkshire", ]
sf:::plot.sf(police_boundaries)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format.R
\name{format_ppp}
\alias{format_ppp}
\title{Convert STATS19 data into ppp (spatstat) format.}
\usage{
format_ppp(data, window = NULL, ...)
}
\arguments{
\item{data}{A STATS19 dataframe to be converted into ppp format.}

\item{window}{A windows of observation, an object of class \code{owin()}. If
\code{window = NULL} (i.e. the default) then the function creates an approximate
bounding box covering the whole UK. It can also be used to filter only the
events occurring in a specific region of UK (see the examples of
\code{\link{get_stats19}}).}

\item{...}{Additional parameters that should be passed to
\code{\link[spatstat.geom:ppp]{spatstat.geom::ppp()}} function. Read the help page of that function
for a detailed description of the available parameters.}
}
\value{
A ppp object.
}
\description{
This function is a wrapper around the \code{\link[spatstat.geom:ppp]{spatstat.geom::ppp()}} function and
it is used to transform STATS19 data into a ppp format.
}
\examples{
if (requireNamespace("spatstat.core", quietly = TRUE)) {
  x_ppp = format_ppp(accidents_sample)
  x_ppp
}

}
\seealso{
\code{\link{format_sf}} for an analogous function used to convert
data into sf format and \code{\link[spatstat.geom:ppp]{spatstat.geom::ppp()}} for the original
spatstat.core function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format.R
\name{format_vehicles}
\alias{format_vehicles}
\title{Format STATS19 vehicles data}
\usage{
format_vehicles(x)
}
\arguments{
\item{x}{Data frame created with \code{read_vehicles()}}
}
\description{
Format STATS19 vehicles data
}
\section{Details}{

This function formats raw STATS19 data
}

\examples{
\donttest{
if(curl::has_internet()) {
dl_stats19(year = 2017, type = "vehicle", ask = FALSE)
x = read_vehicles(year = 2017, format = FALSE)
vehicles = format_vehicles(x)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{stats19_schema}
\alias{stats19_schema}
\alias{stats19_variables}
\title{Stats19 schema and variables}
\description{
\code{stats19_schema} and \code{stats19_variables} contain
metadata on \pkg{stats19} data.
\code{stats19_schema} is a look-up table matching
codes provided in the raw stats19 dataset with
character strings.
}
\note{
The schema data can be (re-)generated using the script in the
\code{data-raw} directory.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjustments.R
\name{get_stats19_adjustments}
\alias{get_stats19_adjustments}
\title{Download and read-in severity adjustment factors}
\usage{
get_stats19_adjustments(
  data_dir = get_data_directory(),
  u = paste0("https://data.dft.gov.uk/road-accidents-safety-data/",
    "accident-and-casualty-adjustment-2004-to-2019.zip"),
  filename = "cas_adjustment_lookup_2019.csv",
  adj_folder = "adjustment-data"
)
}
\arguments{
\item{data_dir}{Where sets of downloaded data would be found.}

\item{u}{The URL of the zip file with adjustments to download}

\item{filename}{The file name of the .csv file in the unzipped folder to read in}

\item{adj_folder}{The folder name where R will look for the unzipped adjustment files}
}
\description{
See the DfT's documentation on adjustment factors
\href{https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/833813/annex-update-severity-adjustments-methodology.pdf}{Annex: Update to severity adjustments methodology}.
}
\details{
See \href{https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/820588/severity-reporting-methodology-final-report.odt}{Estimating and adjusting for changes in the method of severity reporting for road accidents and casualty data: final report}
for details.
}
\examples{
\donttest{
if(curl::has_internet()) {
adjustment = get_stats19_adjustments()
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dl.R
\name{dl_stats19}
\alias{dl_stats19}
\title{Download STATS19 data for a year}
\usage{
dl_stats19(
  year = NULL,
  type = NULL,
  data_dir = get_data_directory(),
  file_name = NULL,
  ask = FALSE,
  silent = FALSE
)
}
\arguments{
\item{year}{A year matching file names on the STATS19
\href{https://data.gov.uk/dataset/cb7ae6f0-4be6-4935-9277-47e5ce24a11f/road-safety-data}{data release page}
e.g. \code{2020}}

\item{type}{One of 'Accident', 'Casualty', 'Vehicle'; defaults to 'Accident'.
Or any variation of to search the file names with such as "acc" or "accid".}

\item{data_dir}{Parent directory for all downloaded files. Defaults to \code{tempdir()}.}

\item{file_name}{The file name (DfT named) to download.}

\item{ask}{Should you be asked whether or not to download the files? \code{TRUE} by default.}

\item{silent}{Boolean. If \code{FALSE} (default value), display useful progress
messages on the screen.}
}
\description{
Download STATS19 data for a year
}
\section{Details}{

This function downloads and unzips UK road crash data.
It results in unzipped .csv files that are put
in the temporary directory specified by \code{get_data_directory()} or provided \code{data_dir}.

The file downloaded would be for a specific year (e.g. 2017).
It could also be a file containing data for a range of two (e.g. 2005-2014).

The \verb{dl_*} functions can download many MB of data so ensure you
have a sufficient internet access and hard disk space.
}

\examples{
\donttest{
if(curl::has_internet()) {
# type by default is accidents table
dl_stats19(year = 2017)
# try another year
dl_stats19(year = 2018)
}
}
}
\seealso{
\code{\link[=get_stats19]{get_stats19()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mot.R
\name{get_MOT}
\alias{get_MOT}
\title{Download vehicle data from the DVSA MOT API using VRM.}
\usage{
get_MOT(vrm, apikey)
}
\arguments{
\item{vrm}{A list of VRMs as character strings.}

\item{apikey}{Your API key as a character string.}
}
\description{
Download vehicle data from the DVSA MOT API using VRM.
}
\section{Details}{

This function takes a a character vector of vehicle registrations (VRMs) and returns vehicle data from MOT records.
It returns a data frame of those VRMs which were successfully used with the DVSA MOT API.

Information on the DVSA MOT API is available here:
https://dvsa.github.io/mot-history-api-documentation/

The DVSA MOT API requires a registration.  The function therefore requires the API key provided by the DVSA.
Be aware that the API has usage limits.  The function will therefore limit lists with more than 150,000 VRMs.
}

\examples{
\donttest{
vrm = c("1RAC","P1RAC")
apikey = Sys.getenv("MOTKEY")
if(nchar(apikey) > 0) {
  get_MOT(vrm = vrm, apikey = apikey)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{locate_files}
\alias{locate_files}
\title{Locate a file on disk}
\usage{
locate_files(
  data_dir = get_data_directory(),
  type = NULL,
  years = NULL,
  quiet = FALSE
)
}
\arguments{
\item{data_dir}{Super directory where dataset(s) were first downloaded to.}

\item{type}{One of 'Accidents', 'Casualties', 'Vehicles'; defaults to 'Accidents', ignores case.}

\item{years}{Years for which data are to be found}

\item{quiet}{Print out messages (files found)}
}
\value{
Character string representing the full path of a single file found,
list of directories where data from the Department for Transport
(stats19::filenames) have been downloaded, or NULL if no files were found.
}
\description{
Helper function to locate files. Given below params, the function
returns 0 or more files found at location/names given.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_accidents}
\alias{read_accidents}
\title{Read in STATS19 road safety data from .csv files downloaded.}
\usage{
read_accidents(
  year = NULL,
  filename = "",
  data_dir = get_data_directory(),
  format = TRUE,
  silent = FALSE
)
}
\arguments{
\item{year}{Single year for which data are to be read}

\item{filename}{Character string of the filename of the .csv to read, if this is given, type and
years determine whether there is a target to read, otherwise disk scan would be needed.}

\item{data_dir}{Where sets of downloaded data would be found.}

\item{format}{Switch to return raw read from file, default is \code{TRUE}.}

\item{silent}{Boolean. If \code{FALSE} (default value), display useful progress
messages on the screen.}
}
\description{
Read in STATS19 road safety data from .csv files downloaded.
}
\section{Details}{

This is a wrapper function to access and load stats 19 data in a user-friendly way.
The function returns a data frame, in which each record is a reported incident in the
STATS19 data.
}

\examples{
\donttest{
if(curl::has_internet()) {
dl_stats19(year = 2019, type = "accident")
ac = read_accidents(year = 2019)

dl_stats19(year = 2019, type = "accident")
ac_2019 = read_accidents(year = 2019)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format.R
\name{format_column_names}
\alias{format_column_names}
\title{Format column names of raw STATS19 data}
\usage{
format_column_names(column_names)
}
\arguments{
\item{column_names}{Column names to be cleaned}
}
\value{
Column names cleaned.
}
\description{
This function takes messy column names and returns clean ones that work well with
R by default. Names that are all lower case with no R-unfriendly characters
such as spaces and \code{-} are returned.
}
\examples{
\donttest{
if(curl::has_internet()) {
crashes_raw = read_accidents(year = 2017)
column_names = names(crashes_raw)
column_names
format_column_names(column_names = column_names)
}
}
}
