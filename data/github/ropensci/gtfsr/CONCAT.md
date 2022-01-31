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
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/

[![Build
Status](http://travis-ci.org/ropensci/gtfsr.svg?branch=master)](http://travis-ci.org/ropensci/gtfsr)
[![codecov.io](http://codecov.io/github/ropensci/gtfsr/coverage.svg?branch=master)](http://codecov.io/github/ropensci/gtfsr?branch=master)
[![](http://badges.ropensci.org/55_status.svg)](https://github.com/ropensci/onboarding/issues/55)

## Description

`gtfsr` is an R package for easily importing, validating, and mapping
transit data that follows the [General Transit Feed Specification
(GTFS)](https://developers.google.com/transit/gtfs/) format.

The `gtfsr` package provides functions for converting files following
the GTFS format into a single `gtfs` data objects. A `gtfs` object can
then be validated for proper data formatting (i.e. if the source data is
properly structured and formatted as a GTFS feed) or have any spatial
data for stops and routes mapped using `leaflet`. The `gtfsr` package
also provides API wrappers for the popular public GTFS feed sharing site
[TransitFeeds](https://transitfeeds.com/), allowing users quick, easy
access to hundreds of GTFS feeds from within R.

## Installation

You can install this package from GitHub using the devtools package:

    if (!require(devtools)) {
        install.packages('devtools')
    }
    devtools::install_github('ropensci/gtfsr')

If you have already installed `gtfsr`, you can get the latest version by
running

    remove.packages('gtfsr')
    devtools::install_github('ropensci/gtfsr')

If you’d like to build the accompanying vignette, then run

    devtools::install_github('ropensci/gtfsr', build_vignettes = TRUE)

## Example Usage

``` r
library(gtfsr)
library(magrittr)
library(dplyr)

# set the API key
# set_api_key() # uncomment to set api key

# get the feedlist dataframe and filter out NYC subway
feedlist_df <- get_feedlist() %>%
  filter(grepl('NYC Subway GTFS', t, ignore.case= TRUE))

# import NYC gtfs feed by sending the url to `import_gtfs`
NYC <- import_gtfs(feedlist_df$url_d)
#> [1] "agency.txt"         "calendar_dates.txt" "calendar.txt"      
#> [4] "routes.txt"         "shapes.txt"         "stop_times.txt"    
#> [7] "stops.txt"          "transfers.txt"      "trips.txt"

# get line (routes) A and B
routes <- NYC[['routes_df']] %>%
  slice(which(grepl('a|b', route_id, ignore.case=TRUE))) %>%
  '$'('route_id')

# take the NYC `gtfs` object and map routes. includes stops by default.
NYC %>% map_gtfs(route_ids = routes)
```

<img src="README/README-readme-body-1.png" width="100%" />

``` r

# gtfs will plot ALL shapes for a given route_ids. These can be reduced using the `service_ids` option.
ids <- NYC$trips_df %>%
  select(route_id, service_id, shape_id) %>%
  distinct() %>%
  filter(route_id %in% routes)
ids %>% head(5) # see all unique combos of ids
#> # A tibble: 5 x 3
#>   route_id service_id   shape_id
#>   <chr>    <chr>        <chr>   
#> 1 A        B20171105WKD A..N43R 
#> 2 A        B20171105WKD A..S43R 
#> 3 A        B20171105WKD A..N85R 
#> 4 A        B20171105WKD A..N54R 
#> 5 A        B20171105WKD A..N65R

# lets map just the the first row
route_ids <- ids$route_id[1]
service_ids <- ids$service_id[1]
shape_ids <- ids$shape_id[1]

# lets map the specific data with some other options enabled.
NYC %>%
  map_gtfs(route_ids = route_ids,
    service_ids = service_ids,
    shape_ids = shape_ids,
    route_colors = 'green', # set the route color
    stop_details = TRUE, # get more stop details on click
    route_opacity = .5) # change the route opacity
```

<img src="README/README-readme-body-2.png" width="100%" />

[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# gtfsr 1.0.3 (2016-07-09)

## Major Changes

- Added conversion function to translate `gtfs_obj` routes to `sf`. New function is called `convert_gtfs_routes_to_sf` (#24)
- All plots now have an overlay feature. The user can enable/disable route/stop layers. (#33)

I did a very poor job logging changes for all the prior releases. My bad.

# gtfsr 0.1.0 (2016-06-23)

- Release of Minimum Viable Product (MVP)
	+ Importing from API or any link
	+ Validating GTFS feed file structure
	+ Basic plotting of stops/routesLeaflet-providers
=================
An extension to [Leaflet](http://leafletjs.com/) that contains configurations for various free<sup>[1](#what-is-free)</sup> tile providers.

# Usage
Leaflet-providers [providers](#providers) are refered to with a `provider[.<variant>]`-string. Let's say you want to add the nice [Watercolor](http://maps.stamen.com/#watercolor/) style from Stamen to your map, you pass `Stamen.Watercolor` to the `L.tileLayer.provider`-constructor, which will return a [L.TileLayer](http://leafletjs.com/reference.html#tilelayer) instance for Stamens Watercolor tile layer.

```Javascript
// add Stamen Watercolor to map.
L.tileLayer.provider('Stamen.Watercolor').addTo(map);
```

## Protocol relativity (`https://`-urls)

Leaflet-providers tries to use `https://` if the page uses `https://` and the provider supports it.
You can force the use of `http://` by passing `force_http: true` in the options argument.

## Retina tiles

Some providers have retina tiles for which the URL only needs to be slightly adjusted, e.g. `-----@2x.png`. For this, add the retina option in the URL, e.g. `-----{retina}.png`, and set a retina value in the options, e.g. `retina: '@2x'`. If Leaflet detects a retina screen (`L.Browser.retina`), the retina option passed to the tileLayer is set to the value supplied, otherwise it's replaced by an empty string.

# Providers

Leaflet-providers provides tile layers from different providers, including *OpenStreetMap*, *Stamen*, *Esri* and *OpenWeatherMap*. The full listing of free to use layers can be [previewed](http://leaflet-extras.github.io/leaflet-providers/preview/index.html). The page will show you the name to use with `leaflet-providers.js` and the code to use it without dependencies.

## Providers requiring registration

In addition to the providers you are free<b id="what-is-free">1</b> to use, we support some layers which require registration.

### HERE (formerly Nokia).

In order to use HERE layers, you must [register](http://developer.here.com/). Once registered, you can create an `app_id` and `app_code` which you have to pass to `L.tileLayer.provider` in the options:

```Javascript
L.tileLayer.provider('HERE.terrainDay', {
    app_id: '<insert ID here>',
    app_code: '<insert ID here>'
}).addTo(map);
```

[Available HERE layers](http://leaflet-extras.github.io/leaflet-providers/preview/#filter=HERE)

### Mapbox

In order to use Mapbox maps, you must [register](https://tiles.mapbox.com/signup). You can get map ID and ACCESS_TOKEN from [Mapbox projects](https://www.mapbox.com/projects):
```JavaScript
L.tileLayer.provider('MapBox', {id: 'ID', accessToken: 'ACCESS_TOKEN'}).addTo(map);
```

### Esri/ArcGIS

In order to use ArcGIS maps, you must [register](https://developers.arcgis.com/en/sign-up/) and abide by the [terms of service](https://developers.arcgis.com/en/terms/). No special syntax is required.

[Available Esri layers](http://leaflet-extras.github.io/leaflet-providers/preview/#filter=Esri)

# Attribution

This work was inspired from <https://gist.github.com/1804938>, and originally created by [Stefan Seelmann](https://github.com/seelmann).

### What do we mean by *free*?
<b id="what-is-free">1</b>
We try to maintain leaflet-providers in such a way that you'll be able to use the layers we include without paying money.
This doesn't mean no limits apply, you should always check before using these layers for anything serious.
So you want to add a layer?
=======

Yay! go add it to the leaflet-providers.js as long as it follows the following 
rules:

- Don't violate a providers TOS (if it exists, include a link to it)
- Don't pre-populate api keys with working keys.
- It should be a basic tile source, no exteral libraries etc.
- The owner hasn't asked us to remove it (hasn't happened yet)Leaflet-providers
=================
An extension to [Leaflet](http://leafletjs.com/) that contains configurations for various free<sup>[1](#what-is-free)</sup> tile providers.

# Usage
Leaflet-providers [providers](#providers) are refered to with a `provider[.<variant>]`-string. Let's say you want to add the nice [Watercolor](http://maps.stamen.com/#watercolor/) style from Stamen to your map, you pass `Stamen.Watercolor` to the `L.tileLayer.provider`-constructor, which will return a [L.TileLayer](http://leafletjs.com/reference.html#tilelayer) instance for Stamens Watercolor tile layer.

```Javascript
// add Stamen Watercolor to map.
L.tileLayer.provider('Stamen.Watercolor').addTo(map);
```

## Protocol relativity (`https://`-urls)

Leaflet-providers tries to use `https://` if the page uses `https://` and the provider supports it.
You can force the use of `http://` by passing `force_http: true` in the options argument.

## Retina tiles

Some providers have retina tiles for which the URL only needs to be slightly adjusted, e.g. `-----@2x.png`. For this, add the retina option in the URL, e.g. `-----{retina}.png`, and set a retina value in the options, e.g. `retina: '@2x'`. If Leaflet detects a retina screen (`L.Browser.retina`), the retina option passed to the tileLayer is set to the value supplied, otherwise it's replaced by an empty string.

# Providers

Leaflet-providers provides tile layers from different providers, including *OpenStreetMap*, *Stamen*, *Esri* and *OpenWeatherMap*. The full listing of free to use layers can be [previewed](http://leaflet-extras.github.io/leaflet-providers/preview/index.html). The page will show you the name to use with `leaflet-providers.js` and the code to use it without dependencies.

## Providers requiring registration

In addition to the providers you are free<b id="what-is-free">1</b> to use, we support some layers which require registration.

### HERE (formerly Nokia).

In order to use HERE layers, you must [register](http://developer.here.com/). Once registered, you can create an `app_id` and `app_code` which you have to pass to `L.tileLayer.provider` in the options:

```Javascript
L.tileLayer.provider('HERE.terrainDay', {
    app_id: '<insert ID here>',
    app_code: '<insert ID here>'
}).addTo(map);
```

[Available HERE layers](http://leaflet-extras.github.io/leaflet-providers/preview/#filter=HERE)

### Mapbox

In order to use Mapbox maps, you must [register](https://tiles.mapbox.com/signup). You can get map ID and ACCESS_TOKEN from [Mapbox projects](https://www.mapbox.com/projects):
```JavaScript
L.tileLayer.provider('MapBox', {id: 'ID', accessToken: 'ACCESS_TOKEN'}).addTo(map);
```

### Esri/ArcGIS

In order to use ArcGIS maps, you must [register](https://developers.arcgis.com/en/sign-up/) and abide by the [terms of service](https://developers.arcgis.com/en/terms/). No special syntax is required.

[Available Esri layers](http://leaflet-extras.github.io/leaflet-providers/preview/#filter=Esri)

# Attribution

This work was inspired from <https://gist.github.com/1804938>, and originally created by [Stefan Seelmann](https://github.com/seelmann).

### What do we mean by *free*?
<b id="what-is-free">1</b>
We try to maintain leaflet-providers in such a way that you'll be able to use the layers we include without paying money.
This doesn't mean no limits apply, you should always check before using these layers for anything serious.
So you want to add a layer?
=======

Yay! go add it to the leaflet-providers.js as long as it follows the following 
rules:

- Don't violate a providers TOS (if it exists, include a link to it)
- Don't pre-populate api keys with working keys.
- It should be a basic tile source, no exteral libraries etc.
- The owner hasn't asked us to remove it (hasn't happened yet)---
output: github_document
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  cache = FALSE,
  comment = "#>",
  message = FALSE,
  error = FALSE,
  warning = FALSE,
  fig.path = "README/README-",
  fig.width=7.3,
  fig.height=5,
  out.width = '100%'
)
```

[![Build Status](http://travis-ci.org/ropensci/gtfsr.svg?branch=master)](http://travis-ci.org/ropensci/gtfsr)
[![codecov.io](http://codecov.io/github/ropensci/gtfsr/coverage.svg?branch=master)](http://codecov.io/github/ropensci/gtfsr?branch=master)
[![](http://badges.ropensci.org/55_status.svg)](https://github.com/ropensci/onboarding/issues/55)


Description
----------

`gtfsr` is an R package for easily importing, validating, and mapping transit data that follows the [General Transit Feed Specification (GTFS)](https://developers.google.com/transit/gtfs/) format.

The `gtfsr` package provides functions for converting files following the GTFS format into a single `gtfs` data objects. A `gtfs` object can then be validated for proper data formatting (i.e. if the source data is properly structured and formatted as a GTFS feed) or have any spatial data for stops and routes mapped using `leaflet`. The `gtfsr` package also provides API wrappers for the popular public GTFS feed sharing site [TransitFeeds](https://transitfeeds.com/), allowing users quick, easy access to hundreds of GTFS feeds from within R.

Installation
-----------------

You can install this package from GitHub using the devtools package:

```
if (!require(devtools)) {
    install.packages('devtools')
}
devtools::install_github('ropensci/gtfsr')
```

If you have already installed `gtfsr`, you can get the latest version by running

```
remove.packages('gtfsr')
devtools::install_github('ropensci/gtfsr')
```

If you'd like to build the accompanying vignette, then run

```
devtools::install_github('ropensci/gtfsr', build_vignettes = TRUE)
```

Example Usage
------------------

```{r readme-body}
library(gtfsr)
library(magrittr)
library(dplyr)

# set the API key
# set_api_key() # uncomment to set api key

# get the feedlist dataframe and filter out NYC subway
feedlist_df <- get_feedlist() %>%
  filter(grepl('NYC Subway GTFS', t, ignore.case= TRUE))

# import NYC gtfs feed by sending the url to `import_gtfs`
NYC <- import_gtfs(feedlist_df$url_d)

# get line (routes) A and B
routes <- NYC[['routes_df']] %>%
  slice(which(grepl('a|b', route_id, ignore.case=TRUE))) %>%
  '$'('route_id')

# take the NYC `gtfs` object and map routes. includes stops by default.
NYC %>% map_gtfs(route_ids = routes)

# gtfs will plot ALL shapes for a given route_ids. These can be reduced using the `service_ids` option.
ids <- NYC$trips_df %>%
  select(route_id, service_id, shape_id) %>%
  distinct() %>%
  filter(route_id %in% routes)
ids %>% head(5) # see all unique combos of ids

# lets map just the the first row
route_ids <- ids$route_id[1]
service_ids <- ids$service_id[1]
shape_ids <- ids$shape_id[1]

# lets map the specific data with some other options enabled.
NYC %>%
  map_gtfs(route_ids = route_ids,
    service_ids = service_ids,
    shape_ids = shape_ids,
    route_colors = 'green', # set the route color
    stop_details = TRUE, # get more stop details on click
    route_opacity = .5) # change the route opacity
```


[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "Getting GTFS Data and Mapping with gtfsr"
author: "Danton Noriega <danton.noriega@gmail.com>"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using gtfsr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r init, warning=FALSE, message=FALSE, echo=FALSE, eval=TRUE}
knitr::opts_chunk$set(collapse = TRUE, cache=FALSE, comment = "#>", fig.width=7.3, fig.height=5)
```

# `gtfsr` (v1.0.3)

`gtfsr` is an R package for easily importing, validating, and mapping transit data that follows the [General Transit Feed Specification (GTFS)](https://developers.google.com/transit/gtfs/) format.

The `gtfsr` package provides functions for converting files following the GTFS format into a single `gtfs` data objects. A `gtfs` object can then be validated for proper data formatting (i.e. if the source data is properly structured and formatted as a GTFS feed) or have any spatial data for stops and routes mapped using `leaflet`. The `gtfsr` package also provides API wrappers for the popular public GTFS feed sharing site [TransitFeeds](https://transitfeeds.com/), allowing users quick, easy access to hundreds of GTFS feeds from within R.

## 1. Get an GTFS API key

This package can get data from a user-specified URL and is also able to get GTFS data from the [TransitFeeds API](http://transitfeeds.com/api/). This vignette will focus on the case where GTFS data is extracted from the TransitFeed API. Below are the steps needed to get a API key (note: requires a GitHub account), including a YouTube (click the GIF to see the YouTube video) that visually guides you through the steps.

1. *Go to [http://transitfeeds.com/](http://transitfeeds.com/)*
2. *Click* "Sign in with GitHub" *in the top-right corner.*
    - If it is your first time visiting the site, it will ask you to sign in (and likely every time if you do not have cookies enabled).
3. *Once signed in, click your profile icon in the top-right and select* "API Keys" *from the drop-down menu.*
    - Your GitHub profile icon and username replaces "Sign in with GitHub".
4. *Fill in* "Enter a description" *and then click the* "Create Key" *button*.
5. *Copy your new API Key to your clipboard.*

[![vid-gif](https://j.gifs.com/kRNVY5.gif)](https://youtu.be/ufM67FoIMho)



## 2. Use `gtfsr` package to download feed list

First things first, load the `gtfsr` package and set your key to access the TransitFeeds API. This example also using the `dplyr` package to manage data frames and `magrittr` for piping.

```{r setup, warning=TRUE, message=FALSE, echo=TRUE, eval=TRUE}
library(gtfsr)
library(dplyr)
options(dplyr.width = Inf) # I like to see all the columns
library(magrittr)

# set_api_key() # input your API key here

```

### Getting full list of available GTFS feeds

With a valid API key loaded, you can easily get the full list of GTFS feeds using the `get_feedlist` function. What we care most about are the feed GTFS data urls contained in column `url_d` of the feed list. Since we are interested in acquiring the GTFS data (not just the feedlist), we can use the `filter_feedlist()` function to return a data frame containing only valid feed urls.

_By default, `filter_feedlist()` only checks to make sure each links starts with `http[s]://`. To check the link is actually working, use option `test_url = TRUE`. But beware, this can take a while!_

```{r get-objs, eval = TRUE, echo = FALSE}
feedlist_df <- readRDS(here::here("data-raw/feedlist_df"))
gtfs_obj    <- readRDS(here::here("data-raw/gtfs_obj"))
```


```{r feedlist, warning=TRUE, message=TRUE, echo=TRUE, eval=FALSE}
feedlist_df <- get_feedlist() # create a data frame of all feeds
```
```{r feedlist-load, echo = TRUE, eval = TRUE}
feedlist_df <- feedlist_df %>% filter_feedlist() # filter the feedlist

feedlist_df %>% select(url_d) %>% head(5) # show first 5 feed urls
```

Here is a map of all available locations.

```{r transitfeeds_map, warning=TRUE, message=TRUE, echo=TRUE, eval=FALSE}
leaflet::leaflet() %>% leaflet::addTiles() %>%
    leaflet::addCircleMarkers(data = feedlist_df, lat = ~loc_lat, lng = ~loc_lng, popup = ~paste(sep = "<br/>", t, loc_t))
```

### Subsetting the GTFS feedlist

If we want only the data for a specific location (or locations), we can get then search the feedlist for feeds of interest.

Assume we are interested in getting all the GTFS data from *Australian* feeds (i.e. we search for location names for the word 'australia'). We can match Australian agencies by name (filter on `loc_t`) and extract the corresponding url feeds (select `url_d`).

```{r aussie, warning=TRUE, message=TRUE, echo=TRUE, eval=TRUE}
## get australian feeds
aussie_df <- feedlist_df %>%
    filter(grepl('australia', loc_t, ignore.case = TRUE)) # filter out locations with "australia" in name

aussie_df %>% select(loc_t) %>% head(5) # look at location names

aussie_urls <- aussie_df %>% select(url_d) # get aussie urls
```

Once we have the urls for the feeds of interest, we can download and extract all the GTFS data into a list of `gtfs` objects using the `import_gtfs` function.

```{r import_gtfs, warning=FALSE, message=FALSE, echo=TRUE, eval=FALSE}
gtfs_objs <- aussie_urls %>% slice(c(6,9)) %>% import_gtfs()
```

### Inspecting Parsing Errors/Warnings

During the import of the any feed url, you will see the following message:

```
NOTE: Parsing errors and warnings while importing data can be extracted from any given data frame with `attr(df, "problems")`.
```

This output was suppressed in the last section to save space given how verbose it is. But the highlighted `NOTE` explains that *if one observes an error or warning during the import process*, one can extract a data frame of problems, which is stored as an attribute for any data frame contained within any `gtfs` object that had a warning output.

As an example, let's extract the gtfs data and problems data for a url with parsing errors/warnings.  You can use `import gtfs` without going through transitfeeds.com if you choose too.

```{r problems, warning=FALSE, message=FALSE, echo=TRUE, eval=TRUE, results='hide'}
url <- 'http://www.co.fairbanks.ak.us/transportation/MACSDocuments/GTFS.zip'

gtfs_obj <- url %>% import_gtfs()
```

If you look at the console output when creating the `gtfs_obj` object, you could see this kind of warning.

```
...
Reading calendar.txt
Warning: 2 parsing failures.
row col   expected    actual
  3  -- 10 columns 1 columns
  4  -- 10 columns 1 columns
...
```

To understand the problem, let's extract the data frame `calendar_df`. Recall that `import_gtfs` returns either a single `gtfs` list object (if one url is provided) or a list of `gtfs` objects.

```{r calendar, echo=TRUE, eval=TRUE}
# extract `calendar_df` from gtfs_obj
df <- gtfs_obj$calendar_df

df

attr(df, 'problems')
```

From inspecting the output from `attr(df, 'problems')` and comparing it to `df`, it appears the problems for this particular `calendar_df` stem from the empty rows added to the end of the original text file. Not a big deal and easily cleaned to fit the standard but we leave such specific fixes to the user to correct.


## 3. Mapping networks, routes, and stops using `gtfsr`

The `gtfsr` has mapping functions designed to help users quickly map spatial data that is found within most GTFS feeds. These functions input `gtfs` objects and then map the desired datum or data (stop, route, route networks).

There are two mapping functions:

1. `map_gtfs` is flexible function used for mapping route shapes and stops. Once can specify the agency (there can be more than one per feed) and/or specific routes by route ID.
2. `map_gtfs_stop` is a simple function used for mapping *a single stop*.


### Example: Duke University

Let's investigate Duke University's transit system.

First, we convert its GTFS transit feed into a `gtfs` object.

```{r duke-extract, include=TRUE}
duke_gtfs_obj <- feedlist_df %>%
    filter(grepl('duke', t, ignore.case=TRUE) & # note, we search `t` (agency name)
           grepl('NC, USA', loc_t, ignore.case=TRUE)) %>%  # get NC agencies
    select(url_d) %>%   # get duke university feed url
    import_gtfs(quiet=TRUE)     # suppress import messages and prints
```

### Mapping an agency route network

We can get visualize all of the routes that make up Duke University's Transit system using `map_gtfs` and just passing the `gtfs` objected `duke_gtfs_obj`. This is because the Duke University Transit system is made of only one agency (`duke_agency_name = "Duke Transit"`) and, when you pass a single `gtfs` object, the default behavior of `map_gtfs` is to take the *first* observed agency name and plot all it's routes.

```{r duke-map2, warning=FALSE, message=FALSE, echo=TRUE, eval=TRUE, cache=FALSE}
map_gtfs(gtfs_obj = duke_gtfs_obj) # map all routes of agency with stops
```


```{r duke-map2a, warning=FALSE, message=FALSE, echo=TRUE, eval=FALSE, cache=FALSE}
# below is equivalent because duke only has a single agency.
duke_agency_name <- duke_gtfs_obj[['agency_df']]$agency_name[1]
map_gtfs(gtfs_obj = duke_gtfs_obj, agency_name = duke_agency_name)
```

If desired, we can also omit stops for every route in the network by using option `include_stops = FALSE` (this option is `include_stops = TRUE` by default).

```{r duke-map2b, warning=FALSE, message=FALSE, echo=TRUE, eval=TRUE, cache=FALSE}
duke_agency_name <- duke_gtfs_obj[['agency_df']]$agency_name[1]
map_gtfs(gtfs_obj = duke_gtfs_obj, agency_name = duke_agency_name, include_stops = FALSE) # map all routes of agency, with no stops
```

### Mapping routes and route stops

Let's get more specific and map out all stops and the shape of the popular *C1 East-West Loop* bus route. We need only find the `route_id` before mapping all the stops using `map_gtfs(..., only_stops = TRUE)` and the shape using `map_gtfs(..., only_stops = FALSE)`.

```{r duke-map1a, warning=FALSE, message=FALSE, echo=TRUE, eval=TRUE, cache=FALSE}
C1_route_id <- duke_gtfs_obj[['routes_df']] %>%
    slice(which(grepl('C1', route_short_name, ignore.case=TRUE))) %>% # search for "C1"
    extract2('route_id') # extract just the datum in route_id

map_gtfs(gtfs_obj = duke_gtfs_obj, route_ids = C1_route_id) # map route shape with stops, the default
map_gtfs(gtfs_obj = duke_gtfs_obj, route_ids = C1_route_id, include_stops = FALSE) # map just the route shape, no stops
map_gtfs(gtfs_obj = duke_gtfs_obj, route_ids = C1_route_id, only_stops = TRUE) # map all stops along route using `only_stops = TRUE`
```

We can also map more than one route *shape* at a time by passing 2 or more route IDs. Let's add the Central Campus Express `CCX`. (Note this feature does not exists for route stops but it's coming soon.)

```{r duke-map1b, warning=FALSE, message=FALSE, echo=TRUE, eval=TRUE, cache=FALSE}
C1_CCX_route_ids <- duke_gtfs_obj[['routes_df']] %>%
    slice(which(grepl('C1|CCX', route_short_name, ignore.case=TRUE))) %>% # search for "C1"
    extract2('route_id') # extract just the datum in route_id

map_gtfs(gtfs_obj = duke_gtfs_obj, route_ids = C1_CCX_route_ids) # pass multiple route IDS and map route shapes with stops (the default)
```

### Mapping a single stop

Sometimes, one wants to see a single stop. For example, the *C1* idles at one of the busiest stops at Duke---the "West Campus Chapel" stop. (This bus stop is located in front of Duke University's iconic gothic Chapel, Duke's most famous landmark.) Let's isolate this stop and map it.


We can search the required field `stop_name` for something that matches "West Campus Chapel" with a combination of `dplyr::slice` plus `which` and `grepl`.

```{r duke-match, include=TRUE}
# look for west chapel stop
west_chapel_stop_id <- duke_gtfs_obj[['stops_df']] %>%
    slice(which(grepl('west campus chapel', stop_name, ignore.case=TRUE))) %>%
    extract2('stop_id') # extract just the stop_id

west_chapel_stop_id
```

Now, we can map the stop using the function `map_gtfs_stop()`.

```{r duke-map, warning=FALSE, message=FALSE, echo=TRUE, eval=FALSE, cache=FALSE}
map_gtfs_stop(gtfs_obj = duke_gtfs_obj, stop_id = west_chapel_stop_id, stop_color = 'blue')
```

## 4. Validating the file and fields structure of a GTFS feed

GTFS feeds contain *required* and *optional* files. And within each of these files, there are also *required* and *optional* fields (For more detailed information, please see Google's [GTFS Feed Specification Reference](https://developers.google.com/transit/gtfs/reference). Information on non-standard GTFS files---specifically `timetables-new.txt` and `timetable_stop_order-new.txt`---can be found at the [GTFS-to-HTML repo](https://github.com/brendannee/gtfs-to-html).

After one has successfully downloaded and unpacked a transit feed, there is no guarantee that it satisfies the requirements of a valid GTFS feed. For example, an unpacked directory may contain all the properly named text files (e.g. `agency.txt`, `stops.txt`, etc), but it could be that within each text file there is no data or that some of the required fields (or variables) (e.g. `stop_id`) are missing.

The `gtfsr` package can quickly check the file and field structure of a GTFS feed and inform you if all required files and fields have been found. Additional information about optional files and fields is also provided. The function is called `validate_gtfs_structure()`. It inputs an object of class `gtfs` (the output of functions `import_gtfs()` or `read_gtfs()`) and by default, attaches the `validate` attribute (i.e. `attr(gtfs_obj, 'validate')`) to the `gtfs` object. The `validate` attribute is just a list of validation information. Set the option `return_gtfs_obj = FALSE` if you only want this validation list.

Let's take a look at an example, using transit feed data from agencies in Durham, NC, USA.


```{r validate, warning=FALSE, message=FALSE, echo=TRUE, eval=TRUE}
nc <- feedlist_df %>%
    filter(grepl('NC, USA', loc_t, ignore.case=TRUE)) # get NC agencies

durham_urls <- nc %>%
    filter(grepl('durham', loc_t, ignore.case=TRUE)) %>%
    select(url_d) # get durham urls

gtfs_objs <- durham_urls %>% import_gtfs(quiet=TRUE) # quietly import

sapply(gtfs_objs, class) # verify that each object of is a `gtfs` object

# validate file and field structures ----------
# attach `validate` data as attribute
gtfs_objs_w_validate <- lapply(gtfs_objs, validate_gtfs_structure)

# extract `validate` attribute data
validate_list_attr <- lapply(gtfs_objs_w_validate, attr, which = 'validate')

# extract validation data directly
validate_list_direct <- lapply(gtfs_objs, validate_gtfs_structure, return_gtfs_obj = FALSE)

# both methods work. option `return_gtfs_obj = FALSE` is more direct
identical(validate_list_attr, validate_list_direct)
```

The `validate` attribute (or list) will always contain 4 elements:

- `all_req_files` a logical value which checks if all *required* files have been found
- `all_req_fields_in_req_files` a logical value which checks if all *required* fields *within required files* have been found
- `all_req_fields_in_opt_files` a logical value which checks if all *required* fields *within any __optional__ files* have been found (i.e. `FALSE` if an optional file is provided but is missing a *required* field)
- `validate_df` a data frame containing all files and fields found plus their status

There can also be 3 other elements:

- `problem_req_files` a data frame which highlights problematic *required* files (required files that are either missing or have missing required fields)
- `problem_opt_files` a data frame which highlights problematic *optional* files (optional files that are missing *required fields*)
- `extra_files` a data frame of any extra files found (i.e. non-standard GTFS feed files not listed as optional or required)

Taking a closer look, we can see that *not* all Durham agencies provide all required files. The second object, `gtfs_objs[[2]]`, is `NULL` given that the link doesn't connect to a valid feed. (The link connects you to [Go Transit NC's Developer Resources page](https://gotransitnc.org/developer-resources/gtfs) but not directly to any feeds.)

The two valid gtfs objects, `gtfs_objs[[1]]` and `gtfs_objs[[3]]`, contain all required fields. However, these agencies provided *optional files* that are missing *required fields*.

```{r durham-validate1, include=TRUE}
validate_list_attr %>% sapply(. %>% extract2('all_req_files'))
validate_list_attr %>% sapply(. %>% extract2('all_req_fields_in_req_files'))
validate_list_attr %>% sapply(. %>% extract2('all_req_fields_in_opt_files'))

# OR, without piping
# sapply(validate_list_attr, '[[', 'all_req_files')
# sapply(validate_list_attr, '[[', 'all_req_fields_in_req_files')
# sapply(validate_list_attr, '[[', 'all_req_fields_in_opt_files')
```

We can get more detail about the problematic optional files by extracting the element `problem_opt_fields`.

```{r durham-opt-files, include=TRUE}
# extract the `problem_opt_files` from the validation list
validate_list_attr[[3]]$problem_opt_files
```

We can see that the optional `frequencies.txt` file was provided but all of the *required fields* were empty.

It is important to recall that GTFS feed files and fields can contain **optional** fields. Therefore, while  it is useful to know any potential problems with optional files provided by a given feed, we can still proceed with interesting analyses as long as we have all the required files and fields.





% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map-gtfs.R
\name{map_gtfs}
\alias{map_gtfs}
\title{General mapping function. Specify a map type and/or a route id.}
\usage{
map_gtfs(gtfs_obj, route_ids = NULL, service_ids = NULL, shape_ids = NULL,
  agency_name = NULL, include_stops = TRUE, only_stops = FALSE,
  stop_details = FALSE, stop_opacity = 0.5, route_opacity = 0.75,
  route_colors = NULL)
}
\arguments{
\item{gtfs_obj}{A GTFS list object with components agency_df, etc.}

\item{route_ids}{Vector (Character). IDs for routes of interest.}

\item{service_ids}{Vector (Character). Service IDs. NULL by Default.}

\item{shape_ids}{Vector (Character). Shape IDs. NULL by Default.}

\item{agency_name}{Character. Provide the name of the agency whose routes are being mapped.}

\item{include_stops}{Boolean. Whether to layer on stops to the route shape. Default is TRUE.}

\item{only_stops}{Boolean. Whether to map only stops, no routes. Overrides \code{include_stops}. Default is FALSE.}

\item{stop_details}{Boolean. Whether to generate detail stop information. Default is FALSE.}

\item{stop_opacity}{Numeric. Value must be between 0 and 1. Defaults is 0.5.}

\item{route_opacity}{Numeric. Value must be between 0 and 1. Default is .75.}

\item{route_colors}{Character. Names of colors (e.g. "blue") or hex values (e.g. '#000000'). Default is NULL.}
}
\value{
Leaflet map object with all stop lat/long values plotted for a route.
}
\description{
General mapping function. Specify a map type and/or a route id.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate-gtfs-structure.R
\name{validate_gtfs_structure}
\alias{validate_gtfs_structure}
\title{Create validation list for a gtfs_obj. It provides an overview of the structure of all files that were imported.}
\usage{
validate_gtfs_structure(gtfs_obj, return_gtfs_obj = TRUE, quiet = FALSE)
}
\arguments{
\item{gtfs_obj}{A GTFS list object with components agency_df, etc.}

\item{return_gtfs_obj}{Boolean. If TRUE, returns gtfs_obj list with a 'validate' attribute appended (TRUE by default),
if FALSE, returns validate list only}

\item{quiet}{Boolean. Option to suppress any messages, prints, etc}
}
\value{
A gtfs_obj list object with attribute 'validate' or just a list containing validation data
}
\description{
Create validation list for a gtfs_obj. It provides an overview of the structure of all files that were imported.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\docType{data}
\name{gtfs_obj}
\alias{gtfs_obj}
\title{Example GTFS data}
\description{
Data obtained from
\url{http://data.trilliumtransit.com/gtfs/duke-nc-us/duke-nc-us.zip}.
}
\seealso{
convert_gtfs_routes_to_sf
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-gtfs-feeds.R
\name{get_feed}
\alias{get_feed}
\title{Download a zipped GTFS feed file from a url}
\usage{
get_feed(url, path = NULL, quiet = FALSE)
}
\arguments{
\item{url}{Character URL of GTFS feed.}

\item{path}{Character. Folder into which to put zipped file. If NULL, then save a tempfile}

\item{quiet}{Boolean. Whether to see file download progress. FALSE by default.}
}
\value{
File path
}
\description{
Download a zipped GTFS feed file from a url
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read-gtfs-files.R
\name{unzip_gtfs_files}
\alias{unzip_gtfs_files}
\title{Unzip GTFS file and delete zip}
\usage{
unzip_gtfs_files(zipfile, delete_zip = FALSE, move_path = NULL,
  quiet = FALSE)
}
\arguments{
\item{zipfile}{path to zipped file}

\item{delete_zip}{Boolean. whether to delete the zipped file after extraction.  Deletes by default.}

\item{move_path}{Character. full file path to desire new location}

\item{quiet}{Boolean. Whether to output files found in folder.}
}
\value{
file path to directory with gtfs .txt files
}
\description{
Unzip GTFS file and delete zip
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert-gtfs.R
\name{routes_df_as_sf}
\alias{routes_df_as_sf}
\title{Get a \code{sf} dataframe for gtfs routes}
\usage{
routes_df_as_sf(gtfs_obj)
}
\arguments{
\item{gtfs_obj}{gtfsr object}
}
\value{
an sf dataframe for gtfs routes with a multilinestring column
}
\description{
Get a \code{sf} dataframe for gtfs routes
}
\examples{
library(gtsf)
routes_sf <- routes_df_as_sf(gtfs_obj)
plot(routes_sf[1,])
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert-gtfs.R
\name{shapes_df_as_sfg}
\alias{shapes_df_as_sfg}
\title{return an sf multilinestring with lat and long from gtfs for a route}
\usage{
shapes_df_as_sfg(df)
}
\arguments{
\item{df}{the shapes_df dataframe from a gtfsr object}
}
\value{
a multilinestring simple feature geometry (sfg) for the routes
}
\description{
return an sf multilinestring with lat and long from gtfs for a route
}
\examples{
shapes_sfg <- shapes_df_as_sfg(gtfs_obj$shapes_df)
plot(shapes_sfg[[1]])
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import-gtfs.R
\name{import_gtfs}
\alias{import_gtfs}
\title{Get a Dataframes of GTFS data.}
\usage{
import_gtfs(paths, local = FALSE, quiet = FALSE)
}
\arguments{
\item{paths}{Character. url links to zip files OR paths to local zip files. if to local path, then option \code{local} must be set to TRUE.}

\item{local}{Boolean. If the paths are searching locally or not. Default is FALSE (that is, urls).}

\item{quiet}{Boolean. Whether to see file download progress and files extract. FALSE by default.}
}
\value{
Dataframes of GTFS data.
}
\description{
Get a Dataframes of GTFS data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\docType{package}
\name{gtfsr-package}
\alias{gtfsr-package}
\alias{gtfsr}
\title{A package for plotting GTFS data.}
\description{
This package allows one to quickly map GTFS data. It also provides a quick and simple way to validate the structure of GTFS data.
}
\examples{

 library(dplyr)
 url <- "http://data.trilliumtransit.com/gtfs/duke-nc-us/duke-nc-us.zip"
 gtfs_obj <- url \%>\% import_gtfs(quiet=TRUE)
 route <- "1693"
 map_gtfs(gtfs_obj, route)
}
\author{
Danton Noriega-Goodwin \email{danton.noriega@gmail.com},
Elaine McVey \email{elaine@transloc.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert-gtfs.R
\name{shape_route_service}
\alias{shape_route_service}
\title{Join the shapes, trips and routes tables together - also checks on some potential errors in the data and warns accordingly}
\usage{
shape_route_service(gtfs_obj, route_ids = NULL, service_ids = NULL)
}
\arguments{
\item{gtfs_obj}{a gtfs object}

\item{route_ids}{the routes for which to join the tables together - required, but not sure why this can't just be any/all routes in routes_df}

\item{service_ids}{\itemize{
\item an optional filter for a certain service-default NULL
}}
}
\value{
shapes_routes_service_df - a dataframe in which routes, services, and shape_ids are all joined
}
\description{
Join the shapes, trips and routes tables together - also checks on some potential errors in the data and warns accordingly
}
\examples{
df <- shape_route_service(gtfs_obj)
#get a summary of the number of shapes and services for a route
library(magrittr)
library(dplyr)
routes_shapes_services <- df \%>\% 
          group_by(route_id) \%>\% 
          summarize(shapes = length(unique(shape_id)), 
          services= length(unique(service_id)))
summary(routes_shapes_services)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-gtfs-feeds.R
\name{clear_api_key}
\alias{clear_api_key}
\title{Clear the API key.}
\usage{
clear_api_key()
}
\description{
Clear the API key.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-gtfs-feeds.R
\name{filter_feedlist}
\alias{filter_feedlist}
\title{Filter a feedlist to include only valid urls (ending in .zip and connect)}
\usage{
filter_feedlist(feedlist_df, test_url = FALSE)
}
\arguments{
\item{feedlist_df}{A dataframe of feed metadata such as output from get_feedlist}

\item{test_url}{Boolean. Whether to test if the url connects or not. FALSE by default (can take a while).}
}
\value{
A dataframe of feed metadata for all feeds in input that are downloadable
}
\description{
Filter a feedlist to include only valid urls (ending in .zip and connect)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map-gtfs-stop.R
\name{map_gtfs_stop}
\alias{map_gtfs_stop}
\title{map a single stop}
\usage{
map_gtfs_stop(gtfs_obj, stop_id, stop_color = NULL)
}
\arguments{
\item{gtfs_obj}{A GTFS list object with components agency_df, etc.}

\item{stop_id}{Character. A single ID for a stop of interest.}

\item{stop_color}{Character. An R color or hex value. Expects single value or NULL. Default is NULL. If length(stop_color) > 1, then it reverts to \code{stop_color <- stop_color[1]}}
}
\value{
Leaflet map object with point plotted at stop lat/long value
}
\description{
map a single stop
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map-gtfs-helpers.R
\name{get_routes_sldf}
\alias{get_routes_sldf}
\title{Get shapes spatial data for given route ids}
\usage{
get_routes_sldf(gtfs_obj, route_ids, service_ids, shape_ids, route_opacity,
  route_colors)
}
\arguments{
\item{gtfs_obj}{A GTFS list object with components agency_df, etc.}

\item{route_ids}{Vector (Character). IDs for routes of interest.}

\item{service_ids}{Vector (Character). Service IDs. NULL by Default.}

\item{shape_ids}{Vector (Character). Shape IDs. NULL by Default.}

\item{route_opacity}{Numeric. Value must be between 0 and 1. Default is NULL.}

\item{route_colors}{Character. Names of colors (e.g. "blue") or hex values (e.g. '#000000').}
}
\value{
Environment containing spatial data, labels, colorings used for plotting
}
\description{
Get shapes spatial data for given route ids
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-gtfs-feeds.R
\name{set_api_key}
\alias{set_api_key}
\title{Set API key for recall}
\usage{
set_api_key()
}
\description{
Set API key for recall
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert-gtfs.R
\name{stops_df_as_sf}
\alias{stops_df_as_sf}
\title{Get a \code{sf} dataframe for gtfs stops}
\usage{
stops_df_as_sf(stops_df)
}
\arguments{
\item{stops_df}{a gtfsr$stops_df dataframe}
}
\value{
an sf dataframe for gtfs routes with a point column
}
\description{
Get a \code{sf} dataframe for gtfs stops
}
\examples{
library(gtsf)
some_stops <- gtfs_obj$stops_df[sample(nrow(gtfs_obj$stops_df), 40),]
some_stops_sf <- stops_df_as_sf(some_stops)
plot(some_stops_sf)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-gtfs-feeds.R
\name{get_api_key}
\alias{get_api_key}
\title{Get API key}
\usage{
get_api_key()
}
\description{
Get API key
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-gtfs-feeds.R
\name{has_api_key}
\alias{has_api_key}
\title{Make sure API key string is not empty}
\usage{
has_api_key()
}
\value{
logical TRUE if key is not empty
}
\description{
Make sure API key string is not empty
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-gtfs-feeds.R
\name{get_feedlist}
\alias{get_feedlist}
\title{Get list of all available feeds from transitfeeds API}
\usage{
get_feedlist()
}
\value{
a data frame with the result of httr::GET
}
\description{
Get list of all available feeds from transitfeeds API
}
