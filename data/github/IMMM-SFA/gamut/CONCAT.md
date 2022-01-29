
<!-- README.md is generated from README.Rmd. Please edit that file -->

![build](https://github.com/IMMM-SFA/gamut/workflows/build/badge.svg)
[![codecov](https://codecov.io/gh/IMMM-SFA/gamut/branch/main/graph/badge.svg?token=uF3EvxwvCO)](https://codecov.io/gh/IMMM-SFA/gamut)
[![DOI](https://zenodo.org/badge/203447802.svg)](https://zenodo.org/badge/latestdoi/203447802)

# `gamut`

## **G**eospatial **A**nalytics for **M**ultisector **U**rban **T**eleconnections

## Description

`gamut` is a tool for exploring *teleconnections* between cities of the
United States and human activities that occur in their associated water
supply catchments. A *teleconnection* is a causal connection or
correlation between human and environmental phenomena that occur a long
distance apart.

<p align="center">
<img src="https://github.com/IMMM-SFA/gamut/blob/main/inst/extdata/gamut_plot.png?raw=true" width="800" height="550">
</p>

## Get Started with `gamut`

`gamut` can be installed remotely from the repository using the R
`devtools` package. From an R prompt, run the command:

``` r
install.packages("devtools")
library(devtools)
devtools::install_github('IMMM-SFA/gamut')
library(gamut)
```

If you run into problems with the remote installation, you may also try
these other options to install gamut:

1.  Save the package file by clicking
    [here](https://api.github.com/repos/IMMM-SFA/gamut/tarball/HEAD),
    then run `install_local()` as shown below:

``` r
install_local('path/to/package')
```

2.  Clone the repo to your computer using
    `git clone "https://github.com/IMMM-SFA/gamut"`. You can then load
    this project into your RStudio and install it.

NOTE: Depending on your version of R, you may need to install Rtools to
retrieve the package. If you have trouble installing it with
`install.packages("Rtools")`, you can find the install file
[here](https://cran.r-project.org/bin/windows/Rtools/). Depending on
your version of `sf`, also may need to install the package `Rcpp` in
order for `gamut` to build correctly.

## Data Files

To download all of the `gamut` input datasets, visit the [Zenodo data
repository](https://zenodo.org/record/5554939#.YV9dnNrMJPY) and download
the zipped data files to your preferred directory. Make sure to combine
the energy, land, misc, and water folders into a single directory. The
table below shows all the files used within the `gamut` software
package.

| gamut Name       | Sub Folder Location | Full File Path                                                                  | Data Source                                                                                              |
|:-----------------|:--------------------|:--------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------|
| watersheds       | water               | water/CWM\_v2\_2/World\_Watershed8.shp                                          | <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF1J67DWR>                                           |
| withdrawal       | water               | water/CWM\_v2\_2/Snapped\_Withdrawal\_Points.shp                                | <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF1J67DWR>                                           |
| citypoint        | water               | water/CWM\_v2\_2/City\_Centroid.shp                                             | <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF1J67DWR>                                           |
| powerplants      | water               | water/UCS-EW3-Energy-Water-Database.xlsx                                        | <https://www.ucsusa.org/resources/ucs-ew3-energy-water-database>                                         |
| crop             | land                | land/2016\_90m\_cdls/cdl\_lowres\_usa.img                                       | <https://www.nass.usda.gov/Research_and_Science/Cropland/Release/>                                       |
| crop\_attributes | land                | land/2016\_90m\_cdls/cdl\_lowres\_usa.img.vat.dbf                               | <https://www.nass.usda.gov/Research_and_Science/Cropland/Release/>                                       |
| irrigation       | land                | land/Version2\_USA\_Demeter.csv                                                 | GCAM Demeter Data                                                                                        |
| nlud             | land                | land/usa\_nlud\_LR.tif                                                          | <https://drive.google.com/file/d/1vmNfwjcaLf0sZTYJ1wsB3liG37sN8gyC/view>                                 |
| hydro            | energy              | energy/EHA\_Public\_PlantFY2019\_GIS\_6/ORNL\_EHAHydroPlant\_FY2020revised.xlsx | <https://hydrosource.ornl.gov/node/250>                                                                  |
| climate          | land                | land/kop\_climate\_classes.tif                                                  | <http://koeppen-geiger.vu-wien.ac.at/present.htm>                                                        |
| HUC4             | water               | water/USA\_HUC4/huc4\_to\_huc2.shp                                              | <http://prd-tnm.s3-website-us-west-2.amazonaws.com/?prefix=StagedProducts/Hydrography/WBD/National/GDB/> |
| population       | land                | land/pden2010\_block/pden2010\_60m.tif                                          | <https://www.sciencebase.gov/catalog/item/57753ebee4b07dd077c70868>                                      |
| runoff           | water               | water/UWSCatCH/Historical\_Mean\_Runoff/USA\_Mean\_Runoff.tif                   | <https://zenodo.org/record/4315195>                                                                      |
| nhd\_flow        | water               | water/UWSCatCH/Watershed\_Flow\_Contributions/UWB\_Intake\_Flows.shp            | <https://zenodo.org/record/4315195>                                                                      |
| contributions    | water               | water/UWSCatCH/Watershed\_Flow\_Contributions/Watershed\_Contributions.csv      | <https://zenodo.org/record/4315195>                                                                      |

## Usage

Once all the necessary data is organized in the correct format, the
package is ready to be used. The main function is
`count_watershed_teleconnections`. Within this function, set up your
data directory and select cities. The city names need to be in the
format of `City | State`, and for multiple cities, place them inside
`c()`. **Available cities can be found on the [Available Cities wiki
page](https://github.com/IMMM-SFA/gamut/wiki/Available-Cities).**

Here is an example of what you would type into your console:

``` r
count_watershed_teleconnections(data_dir = "your/gamut/data_dir", cities = c("Portland | OR", "Knoxville | TN", "New York | NY", "Indianapolis | IN", "Seattle | WA"))
```

The package will cycle through each city and their respective
watersheds, and produce a table with several columns of information. To
learn what each of these variables mean, scroll to the bottom of this
section and see the table of variables. The result of this function will
look something like this:

| city               | city\_population | n\_watersheds | n\_other\_cities | dependent\_city\_pop | watershed\_area\_sqkm | storage\_BCM | yield\_BCM | irr\_cons\_BCM | n\_climate\_zones | n\_hydro\_plants | n\_thermal\_plants | n\_fac\_agcrop | n\_fac\_aglivestock | n\_fac\_cnsmnf | n\_fac\_mining | n\_fac\_oilgas | n\_fac\_total | hydro\_gen\_MWh | thermal\_gen\_MWh | thermal\_cons\_BCM | thermal\_with\_BCM | n\_utilities | n\_ba | n\_crop\_classes | cropland\_fraction | developed\_fraction | ag\_runoff\_max | ag\_runoff\_av\_exgw | ag\_runoff\_av | dev\_runoff\_max | dev\_runoff\_av\_exgw | dev\_runoff\_av | np\_runoff\_max | np\_runoff\_av\_exgw | np\_runoff\_av\_exgw\_unweighted | np\_runoff\_av | n\_economic\_sectors | max\_withdr\_dist\_km | avg\_withdr\_dis\_km | n\_treatment\_plants | watershed\_pop | pop\_cons\_m3sec | av\_fl\_sur\_conc\_pct | av\_fl\_sur\_conc\_pct\_unweighted | av\_ro\_sur\_conc\_pct | av\_fl\_all\_conc\_pct | av\_ro\_all\_conc\_pct | av\_fl\_max\_conc\_pct | av\_ro\_max\_conc\_pct | surface\_contribution\_pct | importance\_of\_worst\_watershed\_pct |
|:-------------------|-----------------:|--------------:|-----------------:|---------------------:|----------------------:|-------------:|-----------:|---------------:|------------------:|-----------------:|-------------------:|---------------:|--------------------:|---------------:|---------------:|---------------:|--------------:|----------------:|------------------:|-------------------:|-------------------:|-------------:|------:|-----------------:|-------------------:|--------------------:|----------------:|---------------------:|---------------:|-----------------:|----------------------:|----------------:|----------------:|---------------------:|---------------------------------:|---------------:|---------------------:|----------------------:|---------------------:|---------------------:|---------------:|-----------------:|-----------------------:|-----------------------------------:|-----------------------:|-----------------------:|-----------------------:|-----------------------:|-----------------------:|---------------------------:|--------------------------------------:|
| Portland \| OR     |           653115 |             1 |                1 |               653115 |              280.7526 |    0.0970009 |  0.0930801 |      0.0000000 |                 1 |                1 |                  0 |              0 |                   0 |              0 |              0 |              0 |             0 |        55263.24 |               0.0 |          0.0000000 |          0.0000000 |            0 |     0 |                0 |          0.0000000 |           0.0044529 |       0.0000000 |            0.0000000 |      0.0000000 |        0.0000172 |             0.0000172 |       0.0000172 |       0.0000172 |            0.0000172 |                        0.0000172 |      0.0000172 |                    4 |             40.923236 |            40.923236 |                    0 |   1.963243e+01 |        0.0000454 |              0.0000000 |                          0.0000000 |              0.0000000 |              0.0000000 |              0.0000000 |               0.000000 |               0.000000 |                        100 |                                   100 |
| Knoxville \| TN    |           187500 |             1 |                2 |               159230 |            23196.0276 |    5.1667173 |  0.2792402 |      0.0132724 |                 2 |               13 |                  4 |              0 |                   7 |            149 |             20 |              0 |           569 |      1633153.23 |         8595013.6 |          0.0105083 |          0.9841422 |            4 |     4 |                7 |          0.0392887 |           0.1217025 |       0.0017350 |            0.0017350 |      0.0017350 |        0.0187317 |             0.0187317 |       0.0187317 |       0.0204667 |            0.0204667 |                        0.0204667 |      0.0204667 |                   13 |              5.942854 |             5.942854 |                    0 |   1.559519e+06 |        3.6100047 |              1.4421179 |                          1.4421179 |              1.7276066 |              1.4421179 |              1.7276066 |               1.442118 |               1.727607 |                        100 |                                   100 |
| New York \| NY     |          8398748 |             8 |                1 |              8398748 |             5203.8799 |    2.3480124 |  0.7446405 |      0.0046190 |                 3 |                5 |                  0 |              0 |                   3 |             45 |              7 |              0 |           377 |       155873.63 |               0.0 |          0.0000000 |          0.0000000 |            0 |     0 |                7 |          0.0473016 |           0.0822482 |       0.0207107 |            0.0034037 |      0.0034037 |        0.3966929 |             0.0263228 |       0.0263228 |       0.3966933 |            0.0297265 |                        0.0667699 |      0.0297265 |                   13 |            191.509637 |           130.711473 |                    0 |   2.566443e+05 |        0.5940853 |              0.1757086 |                          0.3212001 |              0.1946581 |              0.1757086 |              0.1946581 |               1.363024 |               1.150314 |                        100 |                                     5 |
| Indianapolis \| IN |           867125 |             2 |                1 |               867125 |             4564.3384 |    0.1120111 |  0.1861601 |      0.0085113 |                 2 |                0 |                  2 |              0 |                   3 |            226 |             17 |              0 |          2590 |            0.00 |          184897.6 |          0.0001779 |          0.0002461 |            2 |     1 |                7 |          0.6037272 |           0.2275309 |       0.8135638 |            0.7875826 |      0.6536936 |        0.1440485 |             0.1319367 |       0.1095075 |       0.9432512 |            0.9195193 |                        0.8674918 |      0.7632010 |                   11 |             13.022947 |             7.316844 |                    0 |   9.299683e+05 |        2.1527092 |              5.8410899 |                          5.9615696 |              5.1381545 |              4.8481046 |              4.2646682 |               6.137005 |               5.636823 |                         83 |                                    13 |
| Seattle \| WA      |           744955 |             2 |                1 |               744955 |              400.0307 |    0.2083348 |  0.1861601 |      0.0000000 |                 2 |                1 |                  0 |              0 |                   0 |              0 |              0 |              0 |             0 |        74449.12 |               0.0 |          0.0000000 |          0.0000000 |            0 |     0 |                1 |          0.0000223 |           0.0760191 |       0.0000000 |            0.0000000 |      0.0000000 |        0.0066085 |             0.0059192 |       0.0059192 |       0.0066085 |            0.0059192 |                        0.0054596 |      0.0059192 |                    7 |             48.461821 |            42.872284 |                    0 |   1.245807e+03 |        0.0028838 |              0.0000000 |                          0.0000000 |              0.0000000 |              0.0000000 |              0.0000000 |               0.000000 |               0.000000 |                        100 |                                    70 |

This table can be used to compare different variables between multiple
cities. Below is a graph comparing how much developed land are in
cities’ watersheds.

<p align="center">
<img src="https://github.com/IMMM-SFA/gamut/blob/main/inst/extdata/gamut_graph.png?raw=true" width = "650" height="500" >
</p>

The table below shows explanations for each of these variables that are
created through this function:

| Variable Name                         | Description                                                                            | Units                 |
|:--------------------------------------|:---------------------------------------------------------------------------------------|:----------------------|
| city\_population                      | The population of the city being analyzed                                              | people                |
| n\_watersheds                         | Number of watersheds that city uses to source drinking water                           | watersheds            |
| n\_other\_cities                      | Number of other cities pulling off the same watersheds                                 | cities                |
| dependent\_city\_pop                  | Total population of people dependent on that city’s watersheds                         | people                |
| watershed\_area\_sqkm                 | Combined area of all the source watersheds of a city                                   | square kilometers     |
| storage\_BCM                          | Combined storage capacity of all the city catchments                                   | billion cubic meters  |
| yield\_BCM                            | Combined yield capacity of all the city catchments                                     | billion cubic meters  |
| irr\_cons\_BCM                        | Combined water consumption that is used for irrigation with the watersheds             | billion cubic meters  |
| n\_climate\_zones                     | Number of climate zones that the source watersheds cover                               | zones                 |
| n\_hydro\_plants                      | Number of hydroelectric power plants operating within the source watersheds            | plants                |
| n\_thermal\_plants                    | Number of thermal power plants operating within the source watersheds                  | plants                |
| n\_fac\_agcrop                        | Number of agricultural crop facilities within the source watersheds                    | facilities            |
| n\_fac\_aglivestock                   | Number of agicultural livestock facilities within the source watersheds                | facilities            |
| n\_fac\_cnsmnf                        | Number of construction and manufacturing facilities within the source watersheds       | facilities            |
| n\_fac\_mining                        | Number of mining facilities within the source watersheds                               | facilities            |
| n\_fac\_oilgas                        | Number of oil and gas facilities within the source watersheds                          | facilities            |
| n\_fac\_total                         | Total number of facilities operating within the source watersheds                      | facilities            |
| hydro\_gen\_MWh                       | Combined hydroelectric generation from all the facilities within the source watersheds | megawatt-hours        |
| thermal\_gen\_MWh                     | Combined thermal generation from all the facilities within the source watersheds       | megawatt-hours        |
| thermal\_cons\_BCM                    | Combined water consumption that is used for thermal generation                         | billion cubic meters  |
| thermal\_with\_BCM                    | Combined water withdrawal for thermal generation                                       | billion cubic meters  |
| n\_utilities                          | Number of electric utilities within the source watersheds                              | utilities             |
| n\_ba                                 | Number of balancing authorities within the source watersheds                           | balancing authorities |
| n\_crop\_classes                      | Total number of different types of crops within the source watersheds                  | crops                 |
| cropland\_fraction                    | Fraction of land that is used for crops within the source watersheds                   | fraction              |
| developed\_fraction                   | Fraction of land that is developed within the source watersheds                        | fraction              |
| ag\_runoff\_max                       | Max amount of agricultural runoff within the source watersheds                         | fraction              |
| ag\_runoff\_av\_exgw                  | Average agricultural runoff (excluding ground water)                                   | fraction              |
| ag\_runoff\_av                        | Average runoff from agricultural lands                                                 | fraction              |
| dev\_runof\_max                       | Max amount of agricultural runoff within the source watersheds                         | fraction              |
| dev\_runof\_av\_exgw                  | Average developed runoff (excluding ground water)                                      | fraction              |
| dev\_runof\_av                        | Average runoff from developed lands                                                    | fraction              |
| np\_runoff\_max                       | Max amount of non-point source runoff within the source watersheds                     | fraction              |
| np\_runoff\_av\_exgw                  | Average non-point runoff (excluding ground water)                                      | fraction              |
| np\_runoff\_av\_exgw\_unweighted      | Average non-point runoff unweighted (excluding ground water)                           | fraction              |
| np\_runoff\_av                        | Average non-point source runoff.                                                       | fraction              |
| n\_economic\_sectors                  | Total number of different economic sectors within the source watersheds                | sectors               |
| max\_withdr\_dis\_km                  | Maximum distance between a city’s intake points                                        | kilometers            |
| avg\_withdr\_dis\_km                  | Average distance between a city’s intake points                                        | kilometers            |
| n\_treatment\_plants                  | Total number of waste water treatment plants operating within the source watersheds    | plants                |
| watershed\_pop                        | Total number of people living within the source watershed boundaries                   | people                |
| pop\_cons\_m3sec                      | Combined water consumption from the source watersheds that is used for people          | m3/sec                |
| av\_fl\_sur\_conc\_pct                | Average surface flow concentration                                                     | %                     |
| av\_fl\_sur\_conc\_pct\_unweighted    | Average surface flow concentration unweighted                                          | %                     |
| av\_ro\_sur\_conc\_pct                | Average surface runoff concentration                                                   | %                     |
| av\_fl\_all\_conc\_pct                | Average flow concentration                                                             | %                     |
| av\_ro\_all\_conc\_pct                | Average runoff concentration                                                           | %                     |
| av\_fl\_max\_conc\_pct                | Max average flow concentration                                                         | %                     |
| av\_ro\_max\_conc\_pct                | Max average runoff concentration                                                       | %                     |
| surface\_contribution\_pct            | Surface contribution                                                                   | %                     |
| importance\_of\_worst\_watershed\_pct | Measures the importance of the watershed with the worst contamination                  | %                     |

### Dependencies

`gamut` relies on functionality from the following R packages:
clisymbols, crayon, dplyr, dams, exactextractr, foreign, geosphere,
ggplot2, lwgeom, magrittr, purrr, raster, readxl, reservoir, rgdal,
rgeos, sf, sp, spex, stringr, tibble, tidyr, vroom, testthat, knitr,
rmarkdown, knitr.

## Support

For any questions about the package, please contact any of the
contributors below:

Kristian Nelson: <kristian.nelson@pnnl.gov>

Sean Turner: <sean.turner@pnnl.gov>

Chris Vernon <chris.vernon@pnnl.gov>

## Authors and Acknowledgement

Authors: Kristian Nelson, Sean Turner, Chris Vernon, Jennie Rice, Casey
Burleyson, Ryan McManamay, Kerim Dickson

This research was supported by the US Department of Energy, Office of
Science, as part of research in the MultiSector Dynamics, Earth and
Environmental System Modeling Program.
# `gamut` Software Contributions

[repository]: https://github.com/IMMM-SFA/gamut
[issues]: https://github.com/IMMM-SFA/gamut/issues
[new_issue]: https://github.com/IMMM-SFA/gamut/issues/new
[readme]: https://github.com/IMMM-SFA/gamut#readme
[email]: kristian.nelson@pnnl.gov


### Software Questions

If you have a question while using the `gamut` package, first look through the [documentation][readme] to see if the answer is there. If the documentation doesn't answer your question, you can open up an [issue on Github][new_issue]. Explain your question, and the package maintainer will get back to you as soon as possible. If you would like to contact the package maintainer directly, you can do so by [email][email].

### Software Contributions and Guidlines

If you would like to suggest new ideas or functionality for the `gamut` package, go to the [issue page][issues] and create issues for new ideas you might have. If you would like to fix bugs or add in new functions yourself, you may do so by following the development guidelines below.

It is best to follow the standard Github workflow for development.

1. Clone [the repo][repository] to your computer using `git clone "https://github.com/IMMM-SFA/gamut"`. 
2. Open the RStudio `gamut` project file (`.Rproj`), and create a new branch off of the `dev` branch.
3. Once you are on a new branch, make your changes:
    * Write your code.
    * Test your code by running this function: `count_watershed_teleconnections(data_dir = "your/data/dir", cities = c("New York | NY", "Portland | OR"))`. This will make sure that cities with single and multiple watersheds work properly. 
    * Document your code (Ctrl+Shift+D).
    * Check your code with `devtools::check()` and aim for 0 errors. `gamut` currently has warnings from outside packages, so warnings can be ignored for now. 
5. Commit and push your changes.
6. Submit a [pull request](https://github.com/IMMM-SFA/gamut/pulls). A package maintainer will review your pull request and either merge or close it. 

### Report a Bug

If you discover a bug while using the `gamut` package, please create an [issue on GitHub][new_issue] so we can fix it. If possible please include the following in your issue:

* Your operating system name and version.
* Steps of code that reproduce the bug.
---
title: 'gamut: A Geospatial R Package to Analyze Multisectoral Urban Teleconnections'
tags:
- R
- Multisector Dynamics
- Water
- Energy
- Land
- Urban
- Geospatial
date: "30 September 2021"
output:
  html_document:
    df_print: paged
authors:
- name: Kristian D. Nelson
  orcid: 0000-0002-6745-167X
  affiliation: 1
- name: Sean W. Turner
  orcid: 0000-0003-4400-9800
  affiliation: 1
- name: Chris R. Vernon
  orcid: 0000-0002-3406-6214
  affiliation: 1
- name: Jennie S. Rice
  orcid: 0000-0002-7833-9456
  affiliation: 1
bibliography: paper.bib
affiliations:
- name: Pacific Northwest National Laboratory, Richland, Washington, USA
  index: 1
---

### Summary

Most cities in the United States withdraw surface water to meet public water supply needs. The lands on which this water is generated are often developed for human activities&mdash;such as agriculture, mining, and industry&mdash;that may compete for water resources or contaminate water supplies. Cities are thereby connected to other sectors through their water supply catchments. These connections are an example of a multisectoral urban teleconnection, or an interdependency to a geographically disparate region from a source region where events in one (e.g., land use changes) often impact the other [@Seto:2012]. This term was brought about to bring greater understanding of the connections between urbanization and land use changes [@Seto:2012]. The Geospatial Analytics for Multisectoral Urban Teleconnections (`gamut`) package provides national-scale information on these urban teleconnections for 235 cities by combining land use data with hydrological analysis to characterize and quantify urban source watershed human interactions across the conterminous United States (Figure 1).
 
![The `gamut` package analyzes urban cities and their watersheds across the conterminous U.S. As shown in the figure, it can look at characteristics like land use inside watershed boundaries.](gamut_figure.png){ width=95% }

The `gamut` package computes dozens of city-level metrics that inform on the geographical nature of surface water supply catchments and the presence, intensity, and impact of human activities in those catchments. The package cycles through a 3-step process. First, it connects cities to their drinking water resources by using the Urban Water Blueprint dataset [@McDonald:2014]. This dataset is combined with an enhancement dataset of source contribution estimates, river flow statistics, and high resolution runoff data[@Nelson:2021]. After linking the cities to the watersheds, the second step is using the watershed boundary as a mask for the input geospatial layers. Additionally, `gamut` relies heavily on the use of the `st_intersection` function from the [sf](https://cran.r-project.org/web/packages/sf/index.html) package to intersect features and watersheds, which joins the information from the geospatial layers. These layers encompass a wide variety of watershed characteristics including land use, land cover, power generation, infrastructure, hydrological data (irrigation, stream flow, runoff), and population data. Similar to other multisector packages, `gamut` produces new data from a combination of multiple datasets, which are then used to create new statistics. Non-geospatial layers in the form of data tables are joined through city names and municipal IDs. The input layers used in this package have been combined into an open-source dataset and can be accessed [here](https://zenodo.org/record/5554939#.YV9dnNrMJPY) [@Nelson2:2021]. The last step in the process is the creation of output metrics, which is done through the calculation of statistics based on the masked input data for each city.

Creating the links between cities and watershed characteristics enables the `gamut` package to calculate numerous metrics that may be used for multiple types of city-level multisector dynamics research. Metrics reported by `gamut` fall into four main categories: geographical characteristics of watersheds (e.g., climate zones, land area, distance from city, hydrology), potential water contamination concentrations (nonpoint and point), withdrawal/consumption of water from other sectors, and presence/intensity of multisectoral land uses. Table 1 shows the metrics that are created and includes descriptions and units. An R vignette is provided to help users to get started with `gamut` and may be accessed [here](https://immm-sfa.github.io/gamut/).

Table 1: Metrics reported in `gamut`

| Metric Name                            | Description                                                                             | Units                 |
| :------------------------------------ | :-------------------------------------------------------------------------------------- | :-------------------- |
| city\_population                      | The population of the city being analyzed                                               | people                |
| n\_watersheds                         | Number of watersheds that city uses to source drinking water                            | watersheds            |
| n\_other\_cities                      | Number of other cities pulling off the same watersheds                                  | cities                |
| dependent\_city\_pop                  | Total population of people dependent on that city’s watersheds                          | people                |
| watershed\_area\_sqkm                 | Combined area of all the source watersheds of a city                                    | square kilometers     |
| storage\_BCM                          | Combined storage capacity of all the city catchments                                    | billion cubic meters  |
| yield\_BCM                            | Combined yield capacity of all the city catchments                                      | billion cubic meters  |
| irr\_cons\_BCM                        | Combined water consumption that is used for irrigation with the watersheds              | billion cubic meters  |
| n\_climate\_zones                     | Number of climate zones that the source watersheds cover                                | zones                 |
| n\_hydro\_plants                      | Number of hydroelectric power plants operating within the source watersheds            | plants                |
| n\_thermal\_plants                    | Number of thermal power plants operating within the source watersheds                   | plants                |
| n\_fac\_agcrop                        | Number of agricultural crop facilities within the source watersheds                     | facilities            |
| n\_fac\_aglivestock                   | Number of agicultural livestock facilities within the source watersheds                 | facilities            |
| n\_fac\_cnsmnf                        | Number of construction and manufacturing facilities within the source watersheds        | facilities            |
| n\_fac\_mining                        | Number of mining facilities within the source watersheds                                | facilities            |
| n\_fac\_oilgas                        | Number of oil and gas facilities within the source watersheds                           | facilities            |
| n\_fac\_total                         | Total number of facilities operating within the source watersheds                       | facilities            |
| hydro\_gen\_MWh                       | Combined hydroelectric generation from all the facilities within the source watersheds | megawatt-hours         |
| thermal\_gen\_MWh                     | Combined thermal generation from all the facilities within the source watersheds        | megawatt-hours         |
| thermal\_cons\_BCM                    | Combined water consumption that is used for thermal generation                          | billion cubic meters  |
| thermal\_with\_BCM                    | Combined water withdrawal for thermal generation                                        | billion cubic meters  |
| n\_utilities                          | Number of electric utilities within the source watersheds                               | utilities             |
| n\_ba                                 | Number of balancing authorities within the source watersheds                            | balancing authorities |
| n\_crop\_classes                      | Total number of different types of crops within the source watersheds                   | crops                 |
| cropland\_fraction                      | Fraction of land that is used for crops                | fraction                     |
| developed\_fraction                       | Fraction of land that is developed                | fraction                     |
| ag\_runoff\_max                       | Agricultural runoff as proportion of total runoff (worst-case watershed)                | fraction                     |
| ag\_runoff\_av\_exgw                  | Agricultural runoff as proportion of total runoff in supply (exc. groundwater)          | fraction                     |                                
| ag\_runoff\_av                        | Agricultural runoff as proportion of total runoff in supply (inc. groundwater)          | fraction                     |
| dev\_runof\_max                       | Urban runoff as proportion of total runoff (worst-case watershed)                       | fraction                     |
| dev\_runof\_av\_exgw                  | Urban runoff as proportion of total runoff in supply (exc. groundwater)                 | fraction                     |   
| dev\_runof\_av                        | Urban runoff as proportion of total runoff in supply (inc. groundwater)                 | fraction                     |
| np\_runoff\_max                       | Max amount of non-point source runoff within the source watersheds                      | fraction                     |
| np\_runoff\_av\_exgw                  | Nonpoint Proportion of Potentially Contaminated Supply (PPCS) (exc. groundwater)        | fraction                     |
| np\_runoff\_av\_ exgw\_unweighted      | Nonpoint supply contamination averaged across watersheds                                | fraction                     |
| np\_runoff\_av                        | Nonpoint Proportion of Potentially Contaminated Supply (PPCS)                           | fraction                     |
| n\_economic\_sectors                  | Total number of different economic sectors within the source watersheds                 | sectors               |
| max\_withdr\_dis\_km                  | Maximum distance between a city’s intake points                                         | kilometers            |
| avg\_withdr\_dis\_km                  | Average distance between a city’s intake points                                         | kilometers            |
| n\_treatment\_plants                  | Total number of waste water treatment plants operating within the source watersheds     | plants                |
| watershed\_pop                        | Total number of people living within the source watershed boundaries                    | people                |
| pop\_cons\_m3sec                      | Combined water consumption from the source watersheds that is used for people           | m3/sec                |
| av\_fl\_sur\_conc\_pct                | Point PPCS (surface water only, based on flow)                                          | %                     |
| av\_fl\_sur\_ conc\_pct\_unweighted    | Point PPCS (surface water only, based on flow, not weighted by source importance)       | %                     |
| av\_ro\_sur\_conc\_pct                | Point PPCS (surface water only, based on runoff)                                        | %                     |
| av\_fl\_all\_conc\_pct                | Point PPCS (based on flow)                                                              | %                     |
| av\_ro\_all\_conc\_pct                | Point PPCS (based on runoff)                                                            | %                     |
| av\_fl\_max\_conc\_pct                | Point PPCS (based on flow, worst-case catchment only)                                   | %                     |
| av\_ro\_max\_conc\_pct                | Point PPCS (based on runoff, worst-case catchment only)                                 | %                     |
| surface\_contribution\_pct            | Proportion of total average supply made up from surface water                           | %                     |
| importance\_of\_worst\_ watershed\_pct | Proportion of total average supply made up from most heavily contamined watershed       | % |
          |

### Statement of Need

Multisector Dynamics (MSD) research is the study of the co-evolution of human and natural systems. This research requires infrastructure expansion and land use scenarios, resource demand projections, and multisectoral modeling to capture the impacts of trends and shocks on human systems. The `gamut` package offers new data that meet a number of MSD needs. The package may be used to infer possible water resources expansion strategies for major cities in the United States. For example, cities found to be heavily exposed to potential contamination may be more likely to seek alternative means of supply (e.g., water transfers) or invest in water reuse facilities. In a study by @Rice:2013 which looked at de facto wastewater reuse across the US, it was found that there had been an increase in wastewater concentrations in drinking water treatment plants from a 1980 EPA report, especially at low flow conditions. The `gamut` package has the ability to look at wastewater discharge and average flow to find these concentrations at a much larger scale, showing that this package could be useful in studies like this in the future.

In addition to water contamination analysis, the `gamut` package has the ability to reveal which source watersheds are heavily protected by receiving cites. This information can inform land use and energy expansion scenarios applied in MSD research, for example by preventing significant expansion of human developments in protected source watersheds. `gamut` may also be used in large-scale hydrological modeling to correctly assign urban water demands to specific intakes. Whether research is being done on water scarcity, water pollution, or urbanization effects, the `gamut` package provides useful data that can brings greater understanding of anthropogenic impacts on urban source watersheds.

The `gamut` package is open source and may be downloaded using the [devtools](https://devtools.r-lib.org/) package with the code below [@Wickham:2020]. Further instructions on package download can be found in the [documentation](https://github.com/IMMM-SFA/gamut#readme).

```r
install.packages("devtools")
library(devtools)
devtools::install_github('IMMM-SFA/gamut')
library(gamut)
```

### Dependencies

`gamut` relies on functionality from the following R packages:
    clisymbols [@Csardi2:2017],
    crayon [@Csardi:2017],
    dplyr [@Henry2:2020],
    dams [@Goteti:2020],
    exactextractr [@Baston:2020],
    foreign [@R-Core-Team:2020],
    geosphere [@Hijmans:2019],
    ggplot2 [@Wickham:2016],
    lwgeom [@Pebesma:2020],
    magrittr [@Bache:2014],
    purrr [@Henry:2020],
    raster [@Hijmans:2020],
    readxl [@Wickham:2019],
    reservoir [@Turner:2016],
    rgdal [@Bivand2:2020],
    rgeos [@Bivand:2020],
    sf [@Pebesma:2018],
    sp [@Bivand:2013],
    spex [@Sumner:2020],
    stringr [@Wickham2:2019],
    tibble [@Muller:2020],
    tidyr [@Wickham2:2020],
    vroom [@Hester:2021],
    testthat [@Wickham:2011],
    rmarkdown [@Xie:2018],
    knitr [@Xie:2014].


### Acknowledgements

This research was supported by the U.S. Department of Energy, Office of Science, as part of research in [MultisectorDynamics, Earth and Environmental System Modeling Program](https://climatemodeling.science.energy.gov/program/multisector-dynamics).

### References
