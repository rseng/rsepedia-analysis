# lingtypology

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![R build status](https://github.com/ropensci/lingtypology/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/lingtypology/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/lingtypology/master.svg)](https://codecov.io/github/ropensci/lingtypology?branch=master)

[![CRAN version](http://www.r-pkg.org/badges/version/lingtypology)](https://cran.r-project.org/package=lingtypology)
[![](http://cranlogs.r-pkg.org/badges/grand-total/lingtypology)](https://CRAN.R-project.org/package=lingtypology)
[![](https://badges.ropensci.org/95_status.svg)](https://github.com/ropensci/onboarding/issues/95)
[![Research software impact](http://depsy.org/api/package/cran/lingtypology/badge.svg)](http://depsy.org/package/r/lingtypology)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.815028.svg)](https://doi.org/10.5281/zenodo.815028)


`lingtypology` package connects R with the [Glottolog database (v. 4.5)](https://glottolog.org/) and provides additional functionality for linguistic mapping. The Glottolog database contains the catalogue of the world's languages. This package helps researchers to make linguistic maps, using philosophy of [the Cross-Linguistic Linked Data project](https://clld.org/), which uniform access to the data across publications. This package is based on [`leaflet` package](https://rstudio.github.io/leaflet/), so `lingtypology` package is a package for linguistic interactive mapping. You also might be interested in looking into some alternatives to lingtypology:

* [lingtypology](https://pypi.org/project/lingtypology/) in Python by Michael Voronov;
* [glottospace](https://github.com/SietzeN/glottospace) -- R package for the geospatial analysis based on Glottolog by Sietze Norder et al;
* [`lingtypr`](https://gitlab.com/laurabecker/lingtypr) -- R package which partially intersects with `lingtypology` functionality by Laura Becker;
* [`glottoTrees`](https://github.com/erichround/glottoTrees) -- R package for visualising and modifing glottolog trees by Erich Round

## Installation

Get the stable version from CRAN:
```R
install.packages("lingtypology")
```
… or get the development version from GitHub:
```R
install.packages("devtools")
devtools::install_github("ropensci/lingtypology")
```

Sometimes installation failed because of the absence of the package `crosstalk`. Just install it using command `install.packages("crosstalk")`. 

Load a library:
```R
library(lingtypology)
```

For a detailed tutorial see [GitHub pages](https://docs.ropensci.org/lingtypology/).

You can contribute to `lingtypology`, but look through [contribution info](https://github.com/ropensci/lingtypology/blob/master/CONTRIBUTING.md) before.

---

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# Contributing to `lingtypology`

## Issues

When filing an issue, the most important thing is to include a minimal reproducible example so that I can quickly verify the problem. So please include:

* required packages
* package versions
* data
* code

There are some additional information that I put in the issue template.

## Pull requests

To contribute a change to `lingtypology`, you follow these steps:

1. Create a branch in git and make your changes. If you add a new function, do not forget to put some tests and provide your data as contrubuter in the `DESCRIPTION` file.
2. Push branch to github and issue pull request.
3. Discuss the pull request.
# Location
* from: github.com/schloerke/leaflet-providers@urlProtocol

* Inspiration taken from https://github.com/leaflet-extras/leaflet-providers/commit/dea786a3219f9cc824b8e96903a17f46ca9a5afc to use the 'old' relative url protocols and to 'upgrade' them at js runtime.



# Notes...

* Copy/paste provider information into `providers.json`
```js
var providers = L.TileLayer.Provider.providers;
JSON.stringify(providers, null, "  ");
```
  * `./data-raw/providerNames.R` was re-ran to update to the latest providers

* Some providers had their protocols turned into '//'.
  * This allows browsers to pick the protocol
  * To stop files from the protocols staying as files, a ducktape patch was applied to `L.TileLayer.prototype.initialize` and `L.TileLayer.WMS.prototype.initialize`
# Location
* from: github.com/schloerke/leaflet-providers@urlProtocol

* Inspiration taken from https://github.com/leaflet-extras/leaflet-providers/commit/dea786a3219f9cc824b8e96903a17f46ca9a5afc to use the 'old' relative url protocols and to 'upgrade' them at js runtime.



# Notes...

* Copy/paste provider information into `providers.json`
```js
var providers = L.TileLayer.Provider.providers;
JSON.stringify(providers, null, "  ");
```
  * `./data-raw/providerNames.R` was re-ran to update to the latest providers

* Some providers had their protocols turned into '//'.
  * This allows browsers to pick the protocol
  * To stop files from the protocols staying as files, a ducktape patch was applied to `L.TileLayer.prototype.initialize` and `L.TileLayer.WMS.prototype.initialize`
Leaflet-providers
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