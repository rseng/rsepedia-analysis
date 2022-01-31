
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- used devtools::build_readme() to update the md -->

[![CRAN Status
Badge](http://www.r-pkg.org/badges/version/rnaturalearth)](https://cran.r-project.org/package=rnaturalearth)
[![CRAN Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/rnaturalearth)](https://cran.r-project.org/package=rnaturalearth)
[![CRAN Monthly
Downloads](http://cranlogs.r-pkg.org/badges/rnaturalearth)](https://cran.r-project.org/package=rnaturalearth)
[![](https://badges.ropensci.org/22_status.svg)](https://github.com/ropensci/onboarding/issues/22)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![Travis-CI Build
Status](https://travis-ci.org/ropensci/rnaturalearth.svg?branch=master)](https://travis-ci.org/ropensci/rnaturalearth)
[![Build
status](https://ci.appveyor.com/api/projects/status/yp26qgeb1iligrpp?svg=true)](https://ci.appveyor.com/project/AndySouth/rnaturalearth)

# rnaturalearth

An R package to hold and facilitate interaction with [Natural
Earth](http://www.naturalearthdata.com/) map data.

### Provides :

1.  access to a pre-downloaded subset of Natural Earth v4.1.0 (March
    2018) vector data commonly used in world mapping
2.  easy subsetting by countries and regions
3.  functions to download other Natural Earth vector and raster data
4.  a simple, reproducible and sustainable workflow from Natural Earth
    data to rnaturalearth enabling updating as new versions become
    available
5.  clarification of differences in world maps classified by countries,
    sovereign states and map units
6.  consistency with Natural Earth naming conventions so that
    rnaturalearth users can use Natural Earth documentation
7.  data in ‘sf’ or ‘sp’ formats

The [Natural Earth](http://www.naturalearthdata.com/) website structures
vector data by scale, category and type. These determine the filenames
of downloads. rnaturalearth uses this structure to facilitate download
(like an API).

### Install rnaturalearth

Install from CRAN :

    install.packages("rnaturalearth")

or install the development version from GitHub using
[devtools](https://github.com/hadley/devtools).

    devtools::install_github("ropensci/rnaturalearth")

Data to support much of the package functionality are stored in two data
packages that you will be prompted to install when required if you do
not do so here.

    devtools::install_github("ropensci/rnaturalearthdata")
    devtools::install_github("ropensci/rnaturalearthhires")

### First Usage

Here using `sp::plot` as a simple, quick way to plot maps. Maps could
also be made with `ggplot2`, `tmap` or other options. All retrieval
functions accept an argument `returnclass='sf'` to return package `sf`
(Simple Features) objects.

``` r
library(rnaturalearth)
library(sp)
Warning: package 'sp' was built under R version 4.0.5

#world countries
sp::plot(ne_countries())
Warning in wkt(obj): CRS object has no comment
```

![](tools/README-unnamed-chunk-2-1.png)<!-- -->

``` r
#uk
sp::plot(ne_countries(country = 'united kingdom'))
Warning in wkt(obj): CRS object has no comment
```

![](tools/README-unnamed-chunk-2-2.png)<!-- -->

``` r
#states, admin level1 boundaries
sp::plot(ne_states(country = 'spain')) 
Warning in wkt(obj): CRS object has no comment
```

![](tools/README-unnamed-chunk-2-3.png)<!-- -->

### Introductory vignette

``` r
vignette('rnaturalearth', package='rnaturalearth')
```

### To download Natural Earth data not already in the package

There are a wealth of other data available at the [Natural
Earth](http://www.naturalearthdata.com/) website. `rnaturalearth` has
functions to help with download of these data.

The data available are outlined in the two tables below and online
[here](http://www.naturalearthdata.com/downloads/50m-physical-vectors/).

``` 

category   cultural 
                                type scale110 scale50 scale10
1                          countries     TRUE    TRUE    TRUE
2                          map_units     TRUE    TRUE    TRUE
3                       map_subunits    FALSE    TRUE    TRUE
4                        sovereignty     TRUE    TRUE    TRUE
5                     tiny_countries     TRUE    TRUE    TRUE
6                             states    FALSE    TRUE    TRUE
7                   populated_places     TRUE    TRUE    TRUE
8                boundary_lines_land     TRUE    TRUE    TRUE
9                  pacific_groupings     TRUE    TRUE    TRUE
10          breakaway_disputed_areas    FALSE    TRUE    TRUE
11     boundary_lines_disputed_areas    FALSE    TRUE    TRUE
12 boundary_lines_maritime_indicator    FALSE    TRUE    TRUE
13                          airports    FALSE    TRUE    TRUE
14                             ports    FALSE    TRUE    TRUE
15                       urban_areas    FALSE    TRUE    TRUE
16                             roads    FALSE   FALSE    TRUE
17                         railroads    FALSE   FALSE    TRUE

category   physical 
                                 type scale110 scale50 scale10
1                           coastline     TRUE    TRUE    TRUE
2                                land     TRUE    TRUE    TRUE
3                               ocean     TRUE    TRUE    TRUE
4             rivers_lake_centerlines     TRUE    TRUE    TRUE
5  rivers_lake_centerlines_scale_rank    FALSE    TRUE    TRUE
6                       rivers_europe    FALSE   FALSE    TRUE
7                rivers_north_america    FALSE   FALSE    TRUE
8                               lakes     TRUE    TRUE    TRUE
9                     glaciated_areas     TRUE    TRUE    TRUE
10        antarctic_ice_shelves_polys     TRUE    TRUE    TRUE
11                   geographic_lines     TRUE    TRUE    TRUE
12                       graticules_1     TRUE    TRUE    TRUE
13                       graticules_5     TRUE    TRUE    TRUE
14                      graticules_10     TRUE    TRUE    TRUE
15                      graticules_15     TRUE    TRUE    TRUE
16                      graticules_20     TRUE    TRUE    TRUE
17                      graticules_30     TRUE    TRUE    TRUE
18                 wgs84_bounding_box     TRUE    TRUE    TRUE
19                             playas    FALSE    TRUE    TRUE
20                      minor_islands    FALSE   FALSE    TRUE
21                              reefs    FALSE   FALSE    TRUE
```

Specify the `scale`, `category` and `type` of the vector you want as in
the examples below.

``` r
#lakes
lakes110 <- ne_download(scale = 110, type = 'lakes', category = 'physical')
sp::plot(lakes110)

#rivers
rivers50 <- ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
sp::plot(rivers50)
```

### Details of different country definitions and scales

``` r
vignette('what-is-a-country', package='rnaturalearth')
```

## Reproducible download of Natural Earth data into the package

[Script](https://github.com/ropensci/rnaturalearthdata/blob/master/data-raw/data_download_script.r)
used to get data into the accompanying data packages.

## Acknowledgements

Thanks to [Lincoln Mullen](https://github.com/lmullen) for code
structure inspiration from
[USAboundaries](https://github.com/ropensci/USAboundaries), [Hadley
Wickham](https://github.com/hadley) for comments and prompting, [Bob
Rudis](https://github.com/hrbrmstr) for answers to stackoverflow
questions about downloading Natural Earth data into R. The [Natural
Earth team](http://www.naturalearthdata.com/about/contributors/) and
[Nathan Kelso](https://github.com/nvkelso) for providing such a great
resource.

## Potential future work

### potential additional data

1.  Country synonyms lookup
      - dataframe with ISO3 and country synonyms
      - similar to
        <https://github.com/AndySouth/rworldmap/blob/master/data/countrySynonyms.rda>
2.  Country larger regions lookup
      - dataframe with ISO3 and membership of different regional
        groupings, e.g. continent, least developed countries etc.
      - similar to
        <https://github.com/AndySouth/rworldmap/blob/master/data/countryRegions.rda>

### potential additional functions

1.  facilitate joining of user data to country boundaries
      - similar to
        <https://github.com/AndySouth/rworldmap/blob/master/R/joinCountryData2Map.R>
      - … but with a better name
      - similar allowing of join by ISO codes or names, with attempted
        synonym matching
      - similar reporting of country joining success and failure
2.  facilitate subsetting by country groupings
      - e.g. least developed countries etc.

[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
rnaturalearth 0.3.0 2021-10-11
===================

* fix rnaturalearthhires installation #47 thankyou Ian Taylor for #43


rnaturalearth 0.2.0
===================

* add to river options in ne_download() by adding to data_list_physical.csv fixing [#23](https://github.com/ropensci/rnaturalearth/issues/23)
* update data to new version [Natural Earth v4.1](https://www.naturalearthdata.com/blog/miscellaneous/natural-earth-v4-1-0-release-notes/) released May 2018.


rnaturalearth 0.1.0  CRAN
=========================

* Initial release
* sf supportFirst CRAN release 

---
 
  
## Test environments
* local windows R3.3.3
* travis-ci R Under development (unstable) (2017-03-17 r72361) Platform: x86_64-pc-linux-gnu (64-bit)
* win-builder (devel) using R Under development (unstable) (2017-03-17 r72361) platform: x86_64-w64-mingw32 (64-bit)

## R CMD check results

0 errors | 0 warnings

2 notes on Winbuilder
1 New submission
2 * checking package dependencies ... NOTE
Package suggested but not available for checking: 'rnaturalearthhires'

Availability using Additional_repositories specification:
  rnaturalearthhires   yes   http://packages.ropensci.org


## Downstream dependencies
none
First CRAN release after acceptance at rOpenSci.

---
 
   
## Test environments
* local Windows install, R 3.3.2
* win-builder (devel and release)
* ubuntu 12.04 (on travis-ci), R 3.3.0

## R CMD check results

### local
0 errors | 0 warnings | 1 notes

* checking DESCRIPTION meta-information ... NOTE
Authors@R field gives persons with no valid roles:
Lincoln Mullen [rev]
Robin Lovelace [rev]

My understanding is that [rev] is a valid code for reviewer according to the MARC bibliographic standards.

### win-builder



## Downstream dependencies
none



---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- used devtools::build_readme() to update the md -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "",
  fig.path = "tools/README-"
)
```

[![CRAN Status Badge](http://www.r-pkg.org/badges/version/rnaturalearth)](https://cran.r-project.org/package=rnaturalearth)
[![CRAN Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/rnaturalearth)](https://cran.r-project.org/package=rnaturalearth)
[![CRAN Monthly Downloads](http://cranlogs.r-pkg.org/badges/rnaturalearth)](https://cran.r-project.org/package=rnaturalearth)
[![](https://badges.ropensci.org/22_status.svg)](https://github.com/ropensci/onboarding/issues/22)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![Travis-CI Build Status](https://travis-ci.org/ropensci/rnaturalearth.svg?branch=master)](https://travis-ci.org/ropensci/rnaturalearth)
[![Build status](https://ci.appveyor.com/api/projects/status/yp26qgeb1iligrpp?svg=true)](https://ci.appveyor.com/project/AndySouth/rnaturalearth)


# rnaturalearth

An R package to hold and facilitate interaction with [Natural Earth](http://www.naturalearthdata.com/) map data.

### Provides :
1. access to a pre-downloaded subset of Natural Earth v4.1.0 (March 2018) vector data commonly used in world mapping
1. easy subsetting by countries and regions
1. functions to download other Natural Earth vector and raster data
1. a simple, reproducible and sustainable workflow from Natural Earth data to rnaturalearth enabling updating as new versions become available
1. clarification of differences in world maps classified by countries, sovereign states and map units
1. consistency with Natural Earth naming conventions so that rnaturalearth users can use Natural Earth documentation
1. data in 'sf' or 'sp' formats

The [Natural Earth](http://www.naturalearthdata.com/) website structures vector data by scale, category and type. These determine the filenames of downloads. rnaturalearth uses this structure to facilitate download (like an API). 



### Install rnaturalearth

Install from CRAN :

```
install.packages("rnaturalearth")
```

or install the development version from GitHub using [devtools](https://github.com/hadley/devtools). 

```
devtools::install_github("ropensci/rnaturalearth")

```

Data to support much of the package functionality are stored in two data packages that you will be prompted to install when required if you do not do so here.

```
devtools::install_github("ropensci/rnaturalearthdata")
devtools::install_github("ropensci/rnaturalearthhires")
```

### First Usage
Here using `sp::plot` as a simple, quick way to plot maps. Maps could also be made with `ggplot2`, `tmap` or other options. All retrieval functions accept an argument `returnclass='sf'` to return package `sf` (Simple Features) objects. 
```{r, eval=TRUE}
library(rnaturalearth)
library(sp)

#world countries
sp::plot(ne_countries())
#uk
sp::plot(ne_countries(country = 'united kingdom'))
#states, admin level1 boundaries
sp::plot(ne_states(country = 'spain')) 

```

### Introductory vignette
```{r, eval=FALSE}
vignette('rnaturalearth', package='rnaturalearth')
```

### To download Natural Earth data not already in the package
There are a wealth of other data available at the [Natural Earth](http://www.naturalearthdata.com/) website. `rnaturalearth` has functions to help with download of these data.

The data available are outlined in the two tables below and online [here](http://www.naturalearthdata.com/downloads/50m-physical-vectors/).

```{r, eval=TRUE, echo=FALSE}
library(knitr)
for(category in c('cultural','physical'))
{
  df_data <- read.csv( system.file("extdata", paste0("data_list_", category, ".csv"), package = "rnaturalearth") )
  cat("\ncategory  ",category,"\n")

  print(df_data)
  #kable(df_data)
}

```

Specify the `scale`, `category` and `type` of the vector you want as in the examples below.

```{r, eval=FALSE}
#lakes
lakes110 <- ne_download(scale = 110, type = 'lakes', category = 'physical')
sp::plot(lakes110)

#rivers
rivers50 <- ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
sp::plot(rivers50)
```

### Details of different country definitions and scales
```{r, eval=FALSE}
vignette('what-is-a-country', package='rnaturalearth')
```

## Reproducible download of Natural Earth data into the package
[Script](https://github.com/ropensci/rnaturalearthdata/blob/master/data-raw/data_download_script.r) used to get data into the accompanying data packages.

## Acknowledgements
Thanks to [Lincoln Mullen](https://github.com/lmullen) for code structure inspiration from [USAboundaries](https://github.com/ropensci/USAboundaries), [Hadley Wickham](https://github.com/hadley) for comments and prompting, [Bob Rudis](https://github.com/hrbrmstr) for answers to stackoverflow questions about downloading Natural Earth data into R. The [Natural Earth team](http://www.naturalearthdata.com/about/contributors/) and [Nathan Kelso](https://github.com/nvkelso) for providing such a great resource.


## Potential future work

### potential additional data

1. Country synonyms lookup
    + dataframe with ISO3 and country synonyms
    + similar to https://github.com/AndySouth/rworldmap/blob/master/data/countrySynonyms.rda
    
1. Country larger regions lookup
    + dataframe with ISO3 and membership of different regional groupings, e.g. continent, least developed countries etc.
    + similar to https://github.com/AndySouth/rworldmap/blob/master/data/countryRegions.rda


### potential additional functions

1. facilitate joining of user data to country boundaries
    + similar to https://github.com/AndySouth/rworldmap/blob/master/R/joinCountryData2Map.R
    + ... but with a better name
    + similar allowing of join by ISO codes or names, with attempted synonym matching
    + similar reporting of country joining success and failure

1. facilitate subsetting by country groupings
    + e.g. least developed countries etc.

    
[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)    
---
title: "Language support in rnaturalearth"
author: "Andy South"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
#to produce a pdf
#output: rmarkdown::pdf_document
vignette: >
  %\VignetteIndexEntry{Language Support}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

(using languages but also about using with tmap and ggplot)

This vignette shows how [rnaturalearth](https://github.com/ropensci/rnaturalearth) makes it easier to make maps with labels in languages other than Engish, and to make thematic maps if you have data that is referenced in languages other than Engish.

[rnaturalearth](https://github.com/ropensci/rnaturalearth) is an R package to hold and facilitate interaction with natural earth vector map data.

[Natural Earth](http://www.naturalearthdata.com/) is a public domain map dataset including vector country boundaries. 


#### load required packages
```{r, eval=TRUE, echo=TRUE, message=FALSE}
library(rnaturalearth)
library(sp)
library(sf)
library(ggplot2)
library(ggrepel)
#library(tmap)
library(knitr)
```


### testing showing a table of available data
```{r echo = FALSE, results = 'asis'}

#category <- 'physical'

#df_data <- read.csv( system.file("extdata", paste0("data_list_", category, ".csv"), package = "rnaturalearth") )

#convert true,false to 0,1, but loses the name row
#I should just store as 0/1
#df_data <- as.data.frame(lapply(df_data[2:4],as.numeric))

knitr::kable(df_layers_physical, caption = "physical vector data available via ne_download()")

knitr::kable(df_layers_cultural, caption = "cultural vector data available via ne_download()")
```


### Country maps with labels in other languages

```{r, eval=FALSE, echo=TRUE, message=FALSE}
#eval FALSE while testing

# Africa
sp::plot(ne_countries(continent = 'africa'))

sfaf <- ne_countries(continent = 'africa', returnclass = 'sf')
sfafc <- st_centroid(sfaf)

sfaf <- cbind(sfaf, st_coordinates(st_centroid(sfaf$geometry)))

#G = cbind(G, st_coordinates(st_centroid(G$geometry)))

#this adds centroids in the middle of countries
ggplot(sfaf) +
  geom_sf() +
  geom_sf(data=sfafc)

#trying labels in the middle of countries, doesn't quite work needs x,y,label
#but once x & y added on with st_coordinates ...
#seems getting the coords might not be necessary for much longer https://github.com/slowkow/ggrepel/issues/111
#cool nearly there ...
ggplot(sfaf) +
  geom_sf() +
  geom_text_repel(aes(x=X, y=Y, label=name))
  #geom_text_repel(data=sfafc, aes(x=X, y=Y, label=name))
  #geom_label_repel(data=sfafc, aes(x=X, y=Y, label=name))
  #geom_text(data=sfafc, label='name_es', x='X', y='Y')

#getting there, labels still overlap a bit 
#maybe make map bigger to allow space for labels
ggplot(sfaf) +
     geom_sf() +
     geom_text_repel(aes(x=X, y=Y, label=name_es))

# point.padding=NA allows labels to overlap the centroid
ggplot(sfaf) +
     geom_sf() +
     geom_text_repel(aes(x=X, y=Y, label=name_es), point.padding = NA)

#Africa labels just down left & right sides
#works pretty well I think
ggplot(sfaf) +
  geom_sf() +
  xlim(-28,61) +
  geom_text_repel(aes(x=X, y=Y, label=name_es),
    data          = subset(sfaf, X > 21),
    nudge_x       = 60 - subset(sfaf, X > 21)$X,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "y",
    hjust         = 0
  ) +
  geom_text_repel(aes(x=X, y=Y, label=name_es),
    data          = subset(sfaf, X < 21),
    nudge_x       = -19 - subset(sfaf, X < 21)$X,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "y",
    hjust         = 1
  )

#french labels
ggplot(sfaf) +
  geom_sf() +
  xlim(-28,61) +
  geom_text_repel(aes(x=X, y=Y, label=name_fr),
    data          = subset(sfaf, X > 21),
    nudge_x       = 60 - subset(sfaf, X > 21)$X,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "y",
    hjust         = 0
  ) +
  geom_text_repel(aes(x=X, y=Y, label=name_fr),
    data          = subset(sfaf, X < 21),
    nudge_x       = -19 - subset(sfaf, X < 21)$X,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "y",
    hjust         = 1
  )



#tmap good but labels currently overlap

#english labels
tm_shape(sfaf) +
  tm_borders() +
  tm_text("name")

#spanish labels
tm_shape(sfaf) +
     tm_borders() +
     tm_text("name_es")

#other languages de, fr, nl, 
tm_shape(sfaf) +
     tm_borders() +
     tm_text("name_de", size=0.5)

#chinese labels don't currently work, because attribute data is filled with NA
#I suspect I may need to do something different on reading the data in ?
# tm_shape(sfaf) +
#      tm_borders() +
#      tm_text("name_zh")

# The full list of
#    languages is: name_ar, name_bn, name_de, name_en, name_es, name_fr, name_el,
#    name_hi, name_hu, name_id, name_it, name_ja, name_ko, name_nl, name_pl,
#    name_pt, name_ru, name_sv, name_tr, name_vi, and name_zh.
# A 2-character language code decoder ring is here:https://en.wikipedia.org/wiki/List_of_ISO_639-2_codes. 



```



---
title: "Introduction to rnaturalearth."
author: "Andy South"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
#to produce a pdf
#output: rmarkdown::pdf_document
vignette: >
  %\VignetteIndexEntry{Introduction to rnaturalearth.}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

This vignette is an introduction to [rnaturalearth](https://github.com/ropensci/rnaturalearth), an R package to hold and facilitate interaction with natural earth vector map data. `rnaturalearth` is a data package designed to provide map data that can be visualised using other R packages.

[Natural Earth](http://www.naturalearthdata.com/) is a public domain map dataset including vector country and other administrative boundaries. 

[rnaturalearth](https://github.com/ropensci/rnaturalearth) does two main things.

1. Contains pre-downloaded vector maps for : 
    + countries `ne_countries()`
    + states `ne_states()`
    + coastline `ne_coastline()`
    
1. Has `ne_download()` function to facilitate download of other vector and raster maps.


This vignette uses `sp::plot` as a simple, quick way to show how different data can be accessed.`rnaturalearth` is designed to provide data allowing creation of more elaborate maps in other visualisation packages (e.g. `ggplot2`, `tmap` and `choroplethr`).



#### load required packages
```{r, eval=TRUE, echo=TRUE, message=FALSE}
library(rnaturalearth)
library(sp)
```


## 1. Maps in the package.
Pre-downloaded maps can be accessed with : 

1. `ne_countries()` for country (admin-0) boundaries 
1. `ne_states()` for boundaries within countries (admin-1)
1. `ne_coastline()` for world coastline


```{r, eval=TRUE, echo=TRUE, message=FALSE}

# world at small scale (low resolution)
sp::plot(ne_countries(type = 'countries', scale = 'small'))

# countries, UK undivided
sp::plot(ne_countries(country = 'united kingdom', type='countries'))
# map_units, UK divided into England, Scotland, Wales and Northern Ireland
sp::plot(ne_countries(country = 'united kingdom', type='map_units'))

     
# countries, small scale
sp::plot(ne_countries(country = 'united kingdom', scale = 'small'))   

# countries, medium scale
sp::plot(ne_countries(country = 'united kingdom', scale = 'medium'))

```


```{r, eval=FALSE, echo=TRUE, message=FALSE}
# not evaluated because rely on rnaturalearthhires data which are on rOpenSci so CRAN check likely to fail

# countries, large scale
sp::plot(ne_countries(country = 'united kingdom', scale = 'large'))

# states country='united kingdom'
sp::plot(ne_states(country = 'united kingdom'))  
# states geounit='england'
sp::plot(ne_states(geounit = 'england')) 

# states country='france'
sp::plot(ne_states(country = 'france'))

```

```{r, eval=TRUE, echo=TRUE, message=FALSE}

# coastline of the world
# subsetting of coastline is not possible because the Natural Earth data are not attributed in that way
sp::plot(ne_coastline())

```




## 2. Downloading other Natural Earth vectors with ne_download().  

Each [Natural Earth](http://www.naturalearthdata.com/) dataset is characterised on the website according to `scale`, `type` and `category`. [rnaturalearth](https://github.com/ropensci/rnaturalearth) allows you to specify `scale`, `type` and `category` and will construct the url and download the corresponding file.

```{r, eval=FALSE, echo=TRUE, message=FALSE}

# lakes
lakes110 <- ne_download(scale = 110, type = 'lakes', category = 'physical')
sp::plot(lakes110, col = 'blue')

# rivers
rivers110 <- ne_download(scale = 110, type = 'rivers_lake_centerlines', category = 'physical')
sp::plot(rivers110, col = 'blue')

```

### Tables of vector layers available via `ne_download(type=[layer_name], scale=)`
1=available, 0=not

```{r echo = FALSE, results = 'asis'}

knitr::kable(df_layers_physical, caption = "category='physical' vector data available via ne_download()")
```

```{r echo = FALSE, results = 'asis'}

knitr::kable(df_layers_cultural, caption = "category='cultural' vector data available via ne_download()")
```
---
title: "What is a country ?"
author: "Andy South"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
#to produce a pdf
#output: rmarkdown::pdf_document
vignette: >
  %\VignetteIndexEntry{What is a country ?}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

This vignette shows how [rnaturalearth](https://github.com/ropensci/rnaturalearth) allows mapping countries using different definitions of what a country is. What a country is can be more complicated than you might expect.
    
For example, from my own parochial perspective, it allows mapping the UK as a whole or separating out England, Scotland, Wales and Northern Ireland. It also allows you to exclude far away places like the Falkland Islands, or not. Mapping France it allows the inclusion or exclusion of French Guiana and islands in the South Pacific.

[rnaturalearth](https://github.com/ropensci/rnaturalearth) is an R package to hold and facilitate interaction with natural earth vector map data.

[Natural Earth](http://www.naturalearthdata.com/) is a public domain map dataset including vector country boundaries. 

This vignette uses `sp::plot` as a simple, quick way to plot the data obtained using `rnaturalearth`. `rnaturalearth` data can also be used to make more elaborate maps with `ggplot2`, `tmap` and other options.


#### load required packages
```{r, eval=TRUE, echo=TRUE, message=FALSE}
library(rnaturalearth)
library(sp)
```


### Country types : countries, map_units and sovereignty.
Natural Earth data are classified by `countries`, `map_units` and `sovereignty`. Below you will see that specifying `united kingdom` for 

1. `countries` gives the UK undivided   
1. `map_units` gives England, Scotland, Wales and Northern Ireland   
1. `sovereignty` includes the Falkland Islands   

Filtering by `geounit` can give finer control, e.g. to plot Scotland alone, or France without French Guiana.

```{r, eval=TRUE, echo=TRUE, message=FALSE}

# countries, UK undivided
sp::plot(ne_countries(country = 'united kingdom', type = 'countries'))
# map_units, UK divided into England, Scotland, Wales and Northern Ireland
sp::plot(ne_countries(country = 'united kingdom', type = 'map_units'))
# map_units, select by geounit to plot Scotland alone
sp::plot(ne_countries(geounit = 'scotland', type = 'map_units'))
# sovereignty, Falkland Islands included in UK
sp::plot(ne_countries(country = 'united kingdom', type = 'sovereignty'), col = 'red')
sp::plot(ne_coastline(scale = 110), col = 'lightgrey', lty = 3, add = TRUE)

# France, country includes French Guiana
sp::plot(ne_countries(country = 'france'))
# France map_units includes French Guiana too
sp::plot(ne_countries(country = 'france', type = 'map_units'))
# France filter map_units by geounit to exclude French Guiana
sp::plot(ne_countries(geounit = 'france', type = 'map_units'))
# France sovereignty includes South Pacicic islands
sp::plot(ne_countries(country = 'france', type = 'sovereignty'), col = 'red')
sp::plot(ne_coastline(scale = 110), col = 'lightgrey', lty = 3, add = TRUE)

```

### Country scales : small, medium and large.
The different definitions of a country outlined above are available at different scales.
```{r, eval=FALSE, echo=TRUE, message=FALSE}

# countries, large scale
sp::plot(ne_countries(country = 'united kingdom', scale='large'))

# countries, medium scale
# temporarily commented out to avoid utf conversion travis error with Curacao
# sp::plot(ne_countries(country = 'united kingdom', scale = 'medium'))
     
# countries, small scale
sp::plot(ne_countries(country = 'united kingdom', scale = 'small'))     
     
     
```


### States, admin level 1, select by country or geounit. 
```{r, eval=FALSE, echo=TRUE, message=FALSE}

# states country='united kingdom'
sp::plot(ne_states(country = 'united kingdom'))  
# states geounit='england'
sp::plot(ne_states(geounit = 'england')) 

# states country='france'
sp::plot(ne_states(country = 'france'))  
# states geounit='france'
sp::plot(ne_states(geounit = 'france')) 

```


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install-rnaturalearthhires.r
\name{check_rnaturalearthhires}
\alias{check_rnaturalearthhires}
\title{Check whether to install rnaturalearthhires and install if necessary}
\usage{
check_rnaturalearthhires()
}
\description{
If the rnaturalearthhires package is not installed, install it from GitHub using
devtools. If it is not up to date, reinstall it.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install-rnaturalearthdata.r
\name{check_rnaturalearthdata}
\alias{check_rnaturalearthdata}
\title{Check whether to install rnaturalearthdata and install if necessary}
\usage{
check_rnaturalearthdata()
}
\description{
If the rnaturalearthdata package is not installed, install it from GitHub using
devtools. If it is not up to date, reinstall it.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_data.R
\name{ne_git_layer_names}
\alias{ne_git_layer_names}
\title{Create list of layer names and metadata links}
\usage{
ne_git_layer_names(x, scale, getmeta)
}
\arguments{
\item{x}{object returned by ne_git_contents}

\item{scale}{one of \code{110}, \code{50}, \code{10}}

\item{getmeta}{whether to get url of the metadata for each layer}
}
\value{
list of lists with layer names and metadata links.
}
\description{
Parses Natural Earth Github folder content for layer names and metadata
links.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ne_states.r
\name{ne_states}
\alias{ne_states}
\alias{ne_admin1}
\title{Get natural earth world state (admin level 1) polygons}
\usage{
ne_states(
  country = NULL,
  geounit = NULL,
  iso_a2 = NULL,
  spdf = NULL,
  returnclass = c("sp", "sf")
)
}
\arguments{
\item{country}{a character vector of country names.}

\item{geounit}{a character vector of geounit names.}

\item{iso_a2}{a character vector of iso_a2 country codes}

\item{spdf}{an optional alternative states map}

\item{returnclass}{'sp' default or 'sf' for Simple Features}
}
\value{
\code{SpatialPolygonsDataFrame} or \code{sf}
}
\description{
returns state polygons (administrative level 1) for specified countries
}
\examples{

# comparing using country and geounit to filter
if (requireNamespace("rnaturalearthhires")) {
  spdf_france_country <- ne_states(country = 'france')
  spdf_france_geounit <- ne_states(geounit = 'france')
  if (require(sp)) {
     plot(spdf_france_country)
     plot(spdf_france_geounit) 
  
     plot(ne_states(country = 'united kingdom'))  
     plot(ne_states(geounit = 'england'))  
  }
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_data_exist.r
\name{check_data_exist}
\alias{check_data_exist}
\title{check whether the requested data exist on Natural Earth}
\usage{
check_data_exist(
  scale = 110,
  type,
  category = c("cultural", "physical", "raster")
)
}
\arguments{
\item{scale}{scale of map to return, one of \code{110}, \code{50}, \code{10} or \code{'small'}, \code{'medium'}, \code{'large'}}

\item{type}{type of natural earth file to download one of 'countries', 'map_units', 'map_subunits', 'sovereignty', 'states'
OR the portion of any natural earth vector url after the scale and before the . 
e.g. for 'ne_50m_urban_areas.zip' this would be 'urban_areas'
OR the raster filename e.g. for 'MSR_50M.zip' this would be 'MSR_50M'}

\item{category}{one of natural earth categories : 'cultural', 'physical', 'raster'}
}
\value{
TRUE or FALSE
}
\description{
checks from a list dependent on type, category and scale. If it returns FALSE the data may still exist on the website.
Doesn't yet do checking on raster names because I found the naming convention too tricky.
}
\examples{
check_data_exist( scale = 110, category = 'cultural', type = 'countries' )
# type not in list for this category
check_data_exist( scale = 110, category = 'physical', type = 'airports' )
# type in list but scale shows FALSE
check_data_exist( scale = 110, category = 'cultural', type = 'airports' )

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ne_coastline.r
\name{ne_coastline}
\alias{ne_coastline}
\title{Get natural earth world coastline}
\usage{
ne_coastline(scale = 110, returnclass = c("sp", "sf"))
}
\arguments{
\item{scale}{scale of map to return, one of \code{110}, \code{50}, \code{10} or \code{'small'}, \code{'medium'}, \code{'large'}}

\item{returnclass}{'sp' default or 'sf' for Simple Features}
}
\value{
\code{SpatialLinesDataFrame} or \code{sf}
}
\description{
returns world coastline at specified scale
}
\examples{

if (requireNamespace("rnaturalearthdata")) {
   sldf_coast <- ne_coastline()

   if (require(sp)) {
     plot(sldf_coast)
   }
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.r
\docType{data}
\name{df_layers_physical}
\alias{df_layers_physical}
\title{list of physical layers available from Natural Earth}
\format{
A \code{DataFrame}

An object of class \code{data.frame} with 29 rows and 4 columns.
}
\usage{
df_layers_physical
}
\description{
list of physical layers available from Natural Earth
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ne_file_name.r
\name{ne_file_name}
\alias{ne_file_name}
\title{return a natural earth filename based on arguments}
\usage{
ne_file_name(
  scale = 110,
  type = "countries",
  category = c("cultural", "physical", "raster"),
  full_url = FALSE
)
}
\arguments{
\item{scale}{scale of map to return, one of \code{110}, \code{50}, \code{10} or \code{'small'}, \code{'medium'}, \code{'large'}}

\item{type}{type of natural earth file to download one of 'countries', 'map_units', 'map_subunits', 'sovereignty', 'states'
OR the portion of any natural earth vector url after the scale and before the . 
e.g. for 'ne_50m_urban_areas.zip' this would be 'urban_areas'
OR the raster filename e.g. for 'MSR_50M.zip' this would be 'MSR_50M'}

\item{category}{one of natural earth categories : 'cultural', 'physical', 'raster'}

\item{full_url}{whether to return just the filename [default] or the full URL needed for download}
}
\value{
string
}
\description{
returns a string that can then be used to download the file.
}
\examples{
ne_name <- ne_file_name( scale = 110, type = 'countries' )
ne_url  <- ne_file_name( scale = 110, type = 'countries', full_url = TRUE )

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.r
\docType{data}
\name{countries}
\alias{countries}
\alias{countries110}
\title{world country polygons from Natural Earth}
\format{
A \code{SpatialPolygonsDataFrame}

An object of class \code{SpatialPolygonsDataFrame} with 177 rows and 94 columns.
}
\source{
\url{http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries.zip}
}
\usage{
countries110
}
\description{
at 1:110m scale (small). Other data and resolutions are in the packages rnaturalearthdata and rnaturalearthhires.
}
\section{Slots}{

\describe{
\item{\code{data}}{A data frame with country attributes}
}}

\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ne_download.r
\name{ne_download}
\alias{ne_download}
\title{download data from Natural Earth and (optionally) read into R}
\usage{
ne_download(
  scale = 110,
  type = "countries",
  category = c("cultural", "physical", "raster"),
  destdir = tempdir(),
  load = TRUE,
  returnclass = c("sp", "sf")
)
}
\arguments{
\item{scale}{scale of map to return, one of \code{110}, \code{50}, \code{10} or \code{'small'}, \code{'medium'}, \code{'large'}}

\item{type}{type of natural earth file to download one of 'countries', 'map_units', 'map_subunits', 'sovereignty', 'states'
OR the portion of any natural earth vector url after the scale and before the . 
e.g. for 'ne_50m_urban_areas.zip' this would be 'urban_areas'. See Details.
OR the raster filename e.g. for 'MSR_50M.zip' this would be 'MSR_50M'}

\item{category}{one of natural earth categories : 'cultural', 'physical', 'raster'}

\item{destdir}{where to save files, defaults to \code{tempdir()}, \code{getwd()} is also possible.}

\item{load}{TRUE/FALSE whether to load file into R and return}

\item{returnclass}{'sp' default or 'sf' for Simple Features}
}
\value{
A \code{Spatial} object depending on the data (points, lines, polygons or raster), 
   unless load=FALSE in which case it returns the name of the downloaded shapefile (without extension).
}
\description{
returns downloaded data as a spatial object or the filename if \code{load=FALSE}.  
if \code{destdir} is specified the data can be reloaded in a later R session using \code{\link{ne_load}}
with the same arguments.
}
\details{
A non-exhaustive list of datasets available according to \code{scale} specified by the \code{type} param 
  \tabular{lccc}{
         	                   \tab scale = 'small'	\tab scale = 'medium'	\tab scale = 'large' \cr
  category = 'physical', type = '[below]' \cr
  coastline	                 \tab y	        \tab y      	\tab y        \cr
  land     	                 \tab y	        \tab y      	\tab y        \cr
  ocean     	                 \tab y	        \tab y      	\tab y        \cr
  rivers_lake_centerlines     \tab y	        \tab y      	\tab y        \cr
  lakes     	                 \tab y	        \tab y      	\tab y        \cr 
  glaciated_areas     	       \tab y	        \tab y      	\tab y        \cr
  antarctic_ice_shelves_polys \tab -	        \tab y      	\tab y        \cr
  geographic_lines            \tab y	        \tab y      	\tab y        \cr
  graticules_1     	         \tab y	        \tab y      	\tab y        \cr
  graticules_30     	         \tab y	        \tab y      	\tab y        \cr
  wgs84_bounding_box     	   \tab y	        \tab y      	\tab y        \cr
  playas     	               \tab -	        \tab y      	\tab y        \cr
  minor_islands      	       \tab -	        \tab -      	\tab y        \cr
  reefs              	       \tab -	        \tab -      	\tab y        \cr 
  category = 'cultural', type = '[below]'             \cr    
  populated_places        	   \tab y	        \tab y      	\tab y        \cr
  boundary_lines_land \tab y	        \tab y      	\tab y        \cr
  breakaway_disputed_areas \tab -	        \tab y      	\tab y        \cr
  airports              	     \tab -	        \tab y      	\tab y        \cr
  ports              	       \tab -	        \tab y      	\tab y        \cr
  urban_areas              	 \tab -	        \tab y      	\tab y        \cr
  roads              	       \tab -	        \tab -      	\tab y        \cr   
  railroads              	   \tab -	        \tab -      	\tab y        \cr       
  }
}
\examples{
\dontrun{
spdf_world <- ne_download( scale = 110, type = 'countries' )

if (require(sp)) {
  plot(spdf_world)
  plot(ne_download(type = 'populated_places'))
}

# reloading from the saved file in the same session with same arguments
spdf_world2 <-    ne_load( scale = 110, type = 'countries' )

# download followed by load from specified directory will work across sessions
spdf_world <- ne_download( scale = 110, type = 'countries', destdir = getwd() )
spdf_world2 <-    ne_load( scale = 110, type = 'countries', destdir = getwd() )

# for raster, here an example with Manual Shaded Relief (MSR)
# download & load
rst <- ne_download(scale = 50, type = 'MSR_50M', category = 'raster', destdir = getwd())

# load after having downloaded
rst <- ne_load(scale = 50, type = 'MSR_50M', category = 'raster', destdir = getwd())

# plot
library(raster)
raster::plot(rst)
} # end dontrun
}
\seealso{
\code{\link{ne_load}}, pre-downloaded data are available using \code{\link{ne_countries}}, \code{\link{ne_states}}.
Other geographic data are available in the raster package : \code{\link[raster]{getData}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaturalearth-package.r
\docType{package}
\name{rnaturalearth}
\alias{rnaturalearth}
\title{rnaturalearth : world map data from Natural Earth}
\description{
Facilitates world mapping by making \href{http://www.naturalearthdata.com/}{Natural Earth} map data more easily available to R users.
}
\seealso{
\code{\link{ne_countries}} \code{\link{ne_states}} \code{\link{ne_download}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ne_countries.r
\name{ne_countries}
\alias{ne_countries}
\alias{ne_admin0}
\title{Get natural earth world country polygons}
\usage{
ne_countries(
  scale = 110,
  type = "countries",
  continent = NULL,
  country = NULL,
  geounit = NULL,
  sovereignty = NULL,
  returnclass = c("sp", "sf")
)
}
\arguments{
\item{scale}{scale of map to return, one of \code{110}, \code{50}, \code{10} or \code{'small'}, \code{'medium'}, \code{'large'}}

\item{type}{country type, one of 'countries', 'map_units', 'sovereignty', 'tiny_countries'}

\item{continent}{a character vector of continent names to get countries from.}

\item{country}{a character vector of country names.}

\item{geounit}{a character vector of geounit names.}

\item{sovereignty}{a character vector of sovereignty names.}

\item{returnclass}{'sp' default or 'sf' for Simple Features}
}
\value{
\code{SpatialPolygonsDataFrame},\code{SpatialPointsDataFrame} or \code{sf}
}
\description{
returns world country polygons at a specified scale, or points of tiny_countries
}
\examples{
spdf_world <- ne_countries()
spdf_africa <- ne_countries(continent = 'africa')
spdf_france <- ne_countries(country = 'france')

if (require(sp)) {
  plot(spdf_world)
  plot(spdf_africa)
  plot(spdf_france)
}

# get as sf
if (require(sf)) { 
  sf_world <- ne_countries(returnclass='sf')
  plot(sf_world)
}

if (require(rnaturalearthdata) & require(sp)) {
  spdf_tiny_countries <- ne_countries(type = 'tiny_countries', scale = 50) 
  plot(spdf_tiny_countries)  
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.r
\name{get_data}
\alias{get_data}
\title{Get data from within the package}
\usage{
get_data(
  scale = 110,
  type = c("countries", "map_units", "sovereignty", "tiny_countries")
)
}
\arguments{
\item{scale}{scale of map to return, one of \code{110}, \code{50}, \code{10}, \code{'small'}, \code{'medium'}, \code{'large'}}

\item{type}{country type, one of 'countries', 'map_units', 'sovereignty', 'tiny_countries'}
}
\value{
A \code{SpatialPolygonsDataFrame} object.
}
\description{
returns world country polygons at a specified scale, used by ne_countries()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.r
\docType{data}
\name{df_layers_cultural}
\alias{df_layers_cultural}
\title{list of cultural layers available from Natural Earth}
\format{
A \code{DataFrame}

An object of class \code{data.frame} with 43 rows and 4 columns.
}
\usage{
df_layers_cultural
}
\description{
list of cultural layers available from Natural Earth
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ne_as_sf.r
\name{ne_as_sf}
\alias{ne_as_sf}
\title{coerce return object to sf if option set}
\usage{
ne_as_sf(x, returnclass = c("sp", "sf"))
}
\arguments{
\item{x}{scale of map to return, one of \code{110}, \code{50}, \code{10} or \code{'small'}, \code{'medium'}, \code{'large'}}

\item{returnclass}{'sp' default or 'sf' for Simple Features

#not exported}
}
\value{
an sf or sp object
}
\description{
coerce return object to sf if option set
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_data.R
\name{ne_git_contents}
\alias{ne_git_contents}
\title{Return contents of Natural Earth Github directory}
\usage{
ne_git_contents(path)
}
\arguments{
\item{path}{string, one of: \code{'110m_physical'}, \code{'110m_cultural'},
\code{'50m_physical'}, \code{'50m_cultural'}, \code{'10m_physical'},
\code{'10m_cultural'}}
}
\value{
list. Includes parsed json content, http path, and response
  code.
}
\description{
Uses the Github API to return contents of Natural Earth Github directories.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_scale.r
\name{check_scale}
\alias{check_scale}
\title{check that this scale is present in Natural Earth}
\usage{
check_scale(x)
}
\arguments{
\item{x}{scale of map to return, one of \code{110}, \code{50}, \code{10} or \code{'small'}, \code{'medium'}, \code{'large'}}
}
\value{
integer scale of map
}
\description{
check name or numeric scale representations, return numeric one
}
\examples{
# commented out because not exported
# check_scale(110)
# check_scale("small")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install-rnaturalearthdata.r
\name{install_rnaturalearthdata}
\alias{install_rnaturalearthdata}
\title{Install the naturalearthdata package after checking with the user}
\usage{
install_rnaturalearthdata()
}
\description{
Install the naturalearthdata package after checking with the user
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_data.R
\name{ne_find_vector_data}
\alias{ne_find_vector_data}
\title{Return a dataframe of available vector layers on Natural Earth}
\usage{
ne_find_vector_data(
  scale = 110,
  category = c("cultural", "physical"),
  getmeta = FALSE
)
}
\arguments{
\item{scale}{scale of map to return, one of \code{110}, \code{50}, \code{10} or \code{'small'}, \code{'medium'}, \code{'large'}}

\item{category}{one of natural earth categories : 'cultural', 'physical'}

\item{getmeta}{whether to get url of the metadata for each layer}
}
\value{
dataframe with two variables: layer and metadata
}
\description{
Checks the Natural Earth Github repository for current vector layers and
provides the file name required in the type argument of ne_download.
}
\examples{
\dontrun{
ne_find_vector_data(scale = 10, category = "physical")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ne_load.r
\name{ne_load}
\alias{ne_load}
\title{load a Natural Earth vector that has already been downloaded to R using \code{\link{ne_download}}}
\usage{
ne_load(
  scale = 110,
  type = "countries",
  category = c("cultural", "physical", "raster"),
  destdir = tempdir(),
  file_name = NULL,
  returnclass = c("sp", "sf")
)
}
\arguments{
\item{scale}{scale of map to return, one of \code{110}, \code{50}, \code{10} or \code{'small'}, \code{'medium'}, \code{'large'}}

\item{type}{type of natural earth file one of 'countries', 'map_units', 'map_subunits', 'sovereignty', 'states'
OR the portion of any natural earth vector url after the scale and before the . 
e.g. for 'ne_50m_urban_areas.zip' this would be 'urban_areas'
OR the raster filename e.g. for 'MSR_50M.zip' this would be 'MSR_50M'}

\item{category}{one of natural earth categories : 'cultural', 'physical', 'raster'}

\item{destdir}{folder to load files from, default=tempdir()}

\item{file_name}{OPTIONAL name of file (excluding path) instead of natural earth attributes}

\item{returnclass}{'sp' default or 'sf' for Simple Features}
}
\value{
A \code{Spatial} object depending on the data (points, lines, polygons or raster).
}
\description{
returns loaded data as a spatial object.
}
\examples{
\dontrun{
# download followed by load from tempdir() works in same R session
spdf_world <- ne_download( scale = 110, type = 'countries' )
spdf_world2 <-    ne_load( scale = 110, type = 'countries' )

# download followed by load from specified directory works between R sessions
spdf_world <- ne_download( scale = 110, type = 'countries', destdir = getwd() )
spdf_world2 <-    ne_load( scale = 110, type = 'countries', destdir = getwd() )

# for raster
# download & load
rst <- ne_download(scale = 50, type = 'OB_50M', category = 'raster', destdir = getwd())

# load after having downloaded
rst <- ne_load(scale = 50, type = 'OB_50M', category = 'raster', destdir = getwd())

# plot
library(raster)
plot(rst)
} # end dontrun

}
\seealso{
\code{\link{ne_download}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install-rnaturalearthhires.r
\name{install_rnaturalearthhires}
\alias{install_rnaturalearthhires}
\title{Install the naturalearthhires package after checking with the user}
\usage{
install_rnaturalearthhires()
}
\description{
Install the naturalearthhires package after checking with the user
}
