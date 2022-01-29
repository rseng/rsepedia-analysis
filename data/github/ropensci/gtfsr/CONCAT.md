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
- The owner hasn't asked us to remove it (hasn't happened yet)