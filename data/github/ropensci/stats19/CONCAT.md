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
