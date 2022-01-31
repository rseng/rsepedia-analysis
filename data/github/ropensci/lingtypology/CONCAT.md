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
- The owner hasn't asked us to remove it (hasn't happened yet)---
title: "`lingtypology`: Typological databases API"
author: "George Moroz"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{`lingtypology`: Typological databases API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```{r, include=FALSE}
library(lingtypology)
knitr::opts_chunk$set(message= FALSE, eval=FALSE)
```

`lingtypology` provides an ability to download data from these typological databases

* [World Atlas of Language Structures](https://wals.info/)
* [AUTOTYP](https://github.com/autotyp/autotyp-data#the-autotyp-database)
* [PHOIBLE](https://phoible.org/)
* [Affix Borrowing database](https://afbo.info)
* [South American Indigenous Language Structures](https://sails.clld.org/)
* [Austronesian Basic Vocabulary Database](https://abvd.shh.mpg.de/austronesian/)

All database function names have identical structure: **database_name.feature**. All functions have as first argument `feature`. All functions create dataframe with column `language` that can be used in `map.feature()` function. It should be noted that all functions cut out the data that can't be maped, so if you want to prevent functions from this behaviour set argument `na.rm` to `FALSE`.

### 1. WALS
The names of the WALS features can be typed in a lower case. This function preserves coordinates from WALS, so you can map coordinates from the WALS or use coordinates from `lingtypology`.
```{r}
df <- wals.feature(c("1a", "20a"))
head(df)
map.feature(df$language,
            features = df$`1a`,
            latitude = df$latitude,
            longitude = df$longitude,
            label = df$language,
            title = "Consonant Inventories")
```

### 2. AUTOTYP
The AUTOTYP features are listed on [the GitHub page](https://github.com/autotyp/autotyp-data#the-autotyp-database). You can use more human way with spaces.
```{r}
df <- autotyp.feature(c('Gender', 'Numeral classifiers'))
head(df)
map.feature(df$language,
            features = df$NumClass.Presence,
            label = df$language,
            title = "Presence of Numeral Classifiers")
```

### 3. PHOIBLE
I used only four features from PHOIBLE: the number of phonemes, the number of consonants, the number of tones and the number of vowels. If you need only a set of them, just specify it in the `features` argument. Since there is a lot of doubling information in the PHOIBLE database, there is an argument `source`.
```{r}
df <- phoible.feature(source = "UPSID")
head(df)
```

### 4. AfBo
The AfBo database has a lot of features that distinguish affix functions, but again you can use a bare function without any arguments to download the whole database. There will be no difference in time, since this function downloads the whole database to your PC. The main destinction is that this database provides recipient and donor languages, so other column names should be used.

```{r}
df <- afbo.feature(c("adjectivizer", "adverbializer"))
head(df)
map.feature(df$Recipient.name,
            features = df$adjectivizer,
            label = df$Recipient.name,
            title = "Borrowed adjectivizer affixes")
```

### 5. SAILS
The SAILS database provide a lot of [features](https://sails.clld.org/parameters), so the function work with their ids:
```{r}
df <- sails.feature(features = "ics10")
head(df)
map.feature(df$language,
            features = df$ics10_description,
            longitude = df$longitude,
            latitude = df$latitude,
            label = df$language,
            title = "Are there numeral classifiers?")
```

### 6. ABVD
The ABVD database is a lexical database, so it is different from clld databases. First of all, ABVD has its own language classification ids. The information about the same language from different sources can be received from these database different ids. So I select several languages and map them coloring by word with the meaning 'hand'.
```{r}
df <- abvd.feature(c(292, 7))
head(df)
new_df <- df[df$word == "hand",]
map.feature(new_df$language,
            features = new_df$item,
            label = new_df$language)
```

### 7. UraLex

`uralex.feature` downloads data from UraLex basic vocabulary dataset. Original language names are stored in the `language` variable. Converted language names for `map.feature` are stored in the `language2` variable.

```{r}
df <- uralex.feature()
df <- df[df$uralex_mng == "crush",]

map.feature(df$language2,
            label = df$item,
            title = "crush")
```

---
title: "`lingtypology`: creating maps"
author: "George Moroz"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEncoding{UTF-8}
  %\VignetteIndexEntry{`lingtypology`: creating maps}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---
```{r, include=FALSE}
library(lingtypology)
knitr::opts_chunk$set(eval = FALSE)
```

### 1. Base map
The most important part of the `lingtypology` package is the function `map.feature`. This function allows you to produce maps similar to known projects within [the Cross-Linguistic Linked Data philosophy](https://clld.org/), such as [WALS](https://wals.info/) and [Glottolog](https://glottolog.org/):
```{r}
map.feature(c("Adyghe", "Kabardian", "Polish", "Russian", "Bulgarian"))
```

As shown in the picture above, this function generates an interactive Leaflet map. All specific points on the map have a pop-up box that appears when markers are clicked (see section 3.3 for more information about editing pop-up boxes). By default, they contain language names linked to the glottolog site.

If for some reasons you are not using RStudio or you want to automatically create and save a lot of maps, you can save a map to a variable and use the `htmlwidgets` package for saving created maps to an .html file. I would like to thank Timo Roettger for mentioning this problem.

```{r, eval = FALSE}
m <- map.feature(c("Adyghe", "Korean"))
# install.packages("htmlwidgets")
library(htmlwidgets)
saveWidget(m, file="TYPE_FILE_PATH/m.html")
```

There is an export button in RStudio, but for some reason it is not so easy to save a map as a .png or.jpg file using code. [Here](https://stackoverflow.com/a/34672309/6056442/) is a possible solution.

### 2. Set features
The goal of this package is to allow typologists (or any other linguists) to map language features. A list of languages and correspondent features can be stored in a `data.frame` as follows:
```{r}
df <- data.frame(language = c("Adyghe", "Kabardian", "Polish", "Russian", "Bulgarian"),
               features = c("polysynthetic", "polysynthetic", "fusional", "fusional", "fusional"))
df
```

Now we can draw a map:
```{r}
map.feature(languages = df$language,
            features = df$features)
```

If you have a lot of features and they appear in the legend in a senseless order(by default it is ordered alphabetically), you can reorderthem using factors (a vector with ordered levels, for more information see `?factor`). for example, I want the feature polysynthetic to be listed first, followed by fusional:

```{r}
df$features <- factor(df$features, levels = c("polysynthetic", "fusional"))
map.feature(languages = df$language, features = df$features)
```

Like in most functions, it is not necessary to name all arguments, so the same result can be obtained by:
```{r}
map.feature(df$language, df$features)
```

As shown in the picture above, all points are grouped by feature, colored and counted. As before, a pop-up box appears when markers are clicked. A control feature allows users to toggle the visibility of points grouped by feature.

There are several types of variables in R and `map.feature` works differently depending on the variable type. I will use a build in data set `phonological_profiles` that contains 19 languages from [UPSyD database](https://lapsyd.huma-num.fr/lapsyd/). This dataset have three variables: the categorical variable `ejectives` indicates whether some language has any ejective sound, numeric variables `consonants` and `vowels` that contains information about the number of consonants and vowels (based on UPSyD database). We can create two maps with categorical variable and with numeric variable:

```{r}
map.feature(phonological_profiles$language,
            phonological_profiles$ejectives) # categorical
map.feature(phonological_profiles$language,
            phonological_profiles$consonants) # numeric
```

The main point is that for creating a correct map, you should correctly define the type of the variable.

This dataset also can be used to show one other parameter of the `map.feature` function. There are two possible ways to show the World map: with the Atlantic sea or with the Pacific sea in the middle. If you don't need default Pacific view use the `map.orientation` parameter(thanks @languageSpaceLabs and @tzakharko for that idea):
```{r}
map.feature(phonological_profiles$language,
            phonological_profiles$consonants,
            map.orientation = "Atlantic")
```

### 3. Set pop-up boxes
Sometimes it is a good idea to add some additional information (e.g. language affiliation, references or even examples) to pop-up boxes that appear when points are clicked. In order to do so, first of all we need to create an extra vector of strings in our dataframe:
```{r}
df$popup <- aff.lang(df$language)
```

The function `aff.lang()` creates a vector of genealogical affiliations that can be easily mapped:
```{r}
map.feature(languages = df$language, features = df$features, popup = df$popup)
```

Pop-up strings can contain HTML tags, so it is easy to insert a link, a couple of lines, a table or even a video and sound. Here is how pop-up boxes can demonstrate language examples:
```{r}
# change a df$popup vector
df$popup <- c("sɐ s-ɐ-k'ʷɐ<br> 1sg 1sg.abs-dyn-go<br>'I go'",
              "sɐ s-o-k'ʷɐ<br> 1sg 1sg.abs-dyn-go<br>'I go'",
              "id-ę<br> go-1sg.npst<br> 'I go'",
              "ya id-u<br> 1sg go-1sg.npst <br> 'I go'",
              "id-a<br> go-1sg.prs<br> 'I go'")
# create a map

map.feature(df$language,
            features = df$features,
            popup = df$popup)
```

How to say _moon_ in Sign Languages? Here is an example:
```{r}
# Create a dataframe with links to video
sign_df <- data.frame(languages = c("American Sign Language", "Russian-Tajik Sign Language", "French Sign Language"),
                      popup = c("https://media.spreadthesign.com/video/mp4/13/48600.mp4", "https://media.spreadthesign.com/video/mp4/12/17639.mp4", "https://media.spreadthesign.com/video/mp4/10/17638.mp4"))

# Change popup to an HTML code
sign_df$popup <- paste("<video width='200' height='150' controls> <source src='",
                       as.character(sign_df$popup),
                       "' type='video/mp4'></video>", sep = "")
# create a map
map.feature(languages = sign_df$languages, popup = sign_df$popup)
```

### 4. Set labels
An alternative way to add some short text to a map is to use the `label` option.
```{r}
map.feature(df$language, df$features,
            label = df$language)
```

There are some additional arguments for customization: `label.fsize` for setting font size, `label.position` for controlling the label position, and `label.hide` to control the appearance of the label: if `TRUE`, the labels are displayed on mouse over(as on the previous map), if `FALSE`, the labels are always displayed (as on the next map).

```{r}
map.feature(df$language, df$features,
            label = df$language,
            label.fsize = 20,
            label.position = "left",
            label.hide = FALSE)
```

There is an additional tool for emphasis of some points on the map. The argument `label.emphasize` allows to emphasize selected points with the color specified by a user.

```{r}
map.feature(df$language, df$features,
            label = df$language,
            label.fsize = 20,
            label.position = "left",
            label.hide = FALSE,
            label.emphasize = list(2:4, "red"))
```

In this example the first vector of the list in the `label.emphasize` argument is vector `2:4` that produce elements `2`, `3` and `4`. You can create youro wn selected rows. e. g. `c(1, 3, 4)`. The second vector of the list is the string with a color.

### 5. Set coordinates
You can set your own coordinates using the arguments `latitude` and `longitude`. It is important to note, that `lingtypology` works only with decimal degrees (something like this: 0.1), not with degrees, minutes and seconds (something like this: 0° 06′ 0″). I will illustrate this with the dataset `circassian` built into the `lingtypology` package. This dataset comes from fieldwork collected during several expeditions in the period 2011-2016 and contains a list of Circassian villages:
```{r}
head(circassian)
```

In this dataframe you can find variables `latitude` and `longitude` that could be used:
```{r}
map.feature(languages = circassian$language,
            features = circassian$dialect,
            popup = circassian$village,
            latitude = circassian$latitude,
            longitude = circassian$longitude)
```

It is possible to collapse multiple dots into clusters:

```{r}
map.feature(languages = circassian$language,
            features = circassian$dialect,
            popup = circassian$village,
            latitude = circassian$latitude,
            longitude = circassian$longitude, 
            point.cluster = TRUE)
```

### 6. Set colors
You can set your own colors using the argument `color`:
```{r}
df <- data.frame(language = c("Adyghe", "Kabardian", "Polish", "Russian", "Bulgarian"),
                 features = c("polysynthetic", "polysynthetic", "fusional", "fusional", "fusional"))
map.feature(languages = df$language,
            features = df$features,
            color= c("yellowgreen", "navy"))
```

Arguments from [RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer) or [viridis](https://CRAN.R-project.org/package=viridis) also can be used as a color argument:

```{r}
map.feature(phonological_profiles$language,
            phonological_profiles$consonants,
            color= "magma")
```

### 7. Set shapes
For some scientific papers it is not possible to use colors for destinguishing features. In that cases it is posible to use `shape` argument:

```{r}
map.feature(languages = circassian$language,
            features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            shape = TRUE)
```

The argument `shape = TRUE` works fine only with 6 or less levels in `features` variable. If there are more levels in `fetures` argument, user need to provide a vector with corresponding shapes:

```{r}
map.feature(languages = circassian$language,
            features = circassian$dialect,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            shape = 1:10,
            shape.size = 14)
```

Arguments `shape.size` and `shape.color` help to change corresponding features of markers.

### 8. Set control box
The package can generate a control box that allows users to toggle the visibility of some points. To enable it, there is an argument `control` in the `map.feature` function:

```{r}
map.feature(languages = df$language,
            features = df$features,
            control = c("a", "b", "b", "b", "a"))
```

As you can see the `control` and `features` arguments are independent of each other.

### 9. Set an additional set of features using strokes
The `map.feature` function has an additional argument `stroke.features`. Using this argument it becomes possible to show two independent sets of features on one map. By default strokes are colored in grey (so for two levels it will be black and white, for three --- black, grey, white and so on), but you can set your own colors using the argument `stroke.color`:
```{r}
map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude)
```

It is important to note that `stroke.features` can work with `NA` values. The function won't plot anything if there is an `NA` value. Let's set a language value to `NA` in all Baksan villages from the `circassian` dataset.

```{r, message= F}
# create newfeature variable
newfeature <- circassian[,c(5,6)]
# set language feature of the Baksan villages to NA and reduce newfeature from dataframe to vector
newfeature <- replace(newfeature$language, newfeature$language == "Baksan", NA)
# create a map

map.feature(circassian$language,
            features = circassian$dialect,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            stroke.features = newfeature)
```

### 10. Set width and an opacity feature
All markers have their own width and opacity, so you can set it. Just use the arguments `width`, `stroke.radius`, `opacity` and `stroke.opacity`:
```{r}
map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            width = 7, stroke.radius = 13)

map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            opacity = 0.7, stroke.opacity = 0.6)
```

### 11. Customizing legends
By default the legend appears in the top right corner. If there are stroke features, two legends are generated. There are additional arguments that control the appearence and the title of the legends.
```{r}
map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            legend = FALSE, stroke.legend = TRUE)

map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            title = "Circassian dialects", stroke.title = "Languages")
```
The arguments `legend.position` and `stroke.legend.position` allow you to change a legend's position using "topright", "bottomright", "bottomleft" or"topleft" strings.

### 12. Set scale bar
A scale bar is automatically added to a map, but you can control its appearance (set `scale.bar` argument to `TRUE` or`FALSE`) and its position (use `scale.bar.position` argument values "topright", "bottomright", "bottomleft" or"topleft").
```{r}
map.feature(c("Adyghe", "Polish", "Kabardian", "Russian"),
            scale.bar= TRUE,
            scale.bar.position = "topright")
```

### 13. Set layouts
It is possible to use different tiles on the same map using the `tile` argument. For more tiles see [here](https://leaflet-extras.github.io/leaflet-providers/preview/index.html).
```{r}
df <- data.frame(lang = c("Adyghe", "Kabardian", "Polish", "Russian", "Bulgarian"),
                 feature = c("polysynthetic", "polysynthetic", "fusion", "fusion", "fusion"),
                 popup = c("Adyghe", "Adyghe", "Slavic", "Slavic", "Slavic"))

map.feature(df$lang, df$feature, df$popup,
            tile = "Stamen.TonerLite")
```

It is possible to use different map tiles on the same map. Just add a vector with tiles.
```{r}
map.feature(df$lang, df$feature, df$popup,
            tile = c("OpenStreetMap", "Stamen.TonerLite"))
```

It is possible to name tiles using the `tile.name` argument.
```{r}
map.feature(df$lang, df$feature, df$popup,
            tile = c("OpenStreetMap", "Stamen.TonerLite"),
            tile.name = c("colored", "b & w"))
```

It is possible to combine the tiles' control box with the features' control box.
```{r}
map.feature(df$lang, df$feature, df$popup,
            tile = c("OpenStreetMap", "Stamen.TonerLite"),
            control = TRUE)
```

### 14. Add a minimap to a map
It is possible to add a minimap to a map.
```{r}
map.feature(c("Adyghe", "Polish", "Kabardian", "Russian"),
            minimap = TRUE)
```

You can control its appearance (by setting the `minimap` argument to `TRUE` or `FALSE`), its position (by using the values "topright", "bottomright", "bottomleft" or"topleft" of the `minimap.position` argument) and its height and width (with the arguments `minimap.height` and `minimap.width`).
```{r}
map.feature(c("Adyghe", "Polish", "Kabardian", "Russian"),
            minimap = TRUE,
            minimap.position = "topright",
            minimap.height = 100,
            minimap.width = 100)
```

### 15. Add minicharts instead of points
This part is created using the beutifull [`leaflet.minicharts` library](https://CRAN.R-project.org/package=leaflet.minicharts). The argument `minichart` allows you to add piecharts or barplots instead of standard point markers. In this part I will use a build in data set `phonological_profiles` that contains 19 languages from [UPSyD database](https://lapsyd.huma-num.fr/lapsyd/). Here is an example of barplot:

```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("vowels", "consonants")])
```

Here is an example of piechart:
```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("vowels", "consonants")],
            minichart = "pie")
```

Colors and opacity could be changed, legend moved:
```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("vowels", "consonants")],
            color= c("yellowgreen", "navy"),
            opacity = 0.7,
            label = phonological_profiles$language,
            legend.position = "topleft")
```

It is possible to add values using argument `minichart.labels`:
```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("vowels", "consonants")],
            minichart = "pie",
            minichart.labels = TRUE)
```

It is also possible to use pie chart in non-convenient way: just indicating with `TRUE` or `FALSE` of pressence of some feature (thanks to Diana Forker for the task!):

```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("tone", "long_vowels", "stress", "ejectives")],
            minichart = "pie", 
            width = 3)
```

Unfortunately this kind of visualisation doesn't work, when you have some lines in your dataset that contain only `FALSE` values. This is **non-convenient** way of category visualisation, so visualisation experts could have a negative opinion about it. This kind of visualisation is also bad for huge number of variables.

### 16. Add a rectangle to a map
It is possible to highlight some part of your map with a rectangle. You need to provide a latitude and longitude of the diagonal (`rectangle.lat` and `rectangel.lng`) and color of the rectangle (`rectangle.color`):

```{r}
map.feature(circassian$language,
            circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            rectangle.lng = c(42.7, 45),
            rectangle.lat = c(42.7, 44.4),
            rectangle.color= "green")
```

### 17. Add a density contourplot to a map
Sometimes it is easier to look at a density contourplot. It can be created using `density.estimation` argument. There are two possibility for creation a density contourplot in `lingtypology`:

* `density.method = "fixed distance"`. First algorithm creates circle polygons with fixed radius around each point and then merge all polygons that are overlapped. It has only one parameter that should be estimated: radius of the circle (`density.width`).
* `density.method = "kernal density estimation"`. Second algorithm uses a kernal density estimation and has two parameters that should be estimated: latitude and longitude bandwidths (`density.width[1]` and (`density.width[2]`))


```{r}
map.feature(circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = circassian$language,
            density.width = 0.15)
```

Density estimation plot can be separated by `features` variable:
```{r}
map.feature(circassian$language,
            features = circassian$dialect,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = circassian$language,
            density.width = 0.15)
```

It is possible to remove points and display only the kernal density estimation plot, using the `density.points` argument:

```{r}
map.feature(circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = circassian$language,
            density.width = 0.15,
            density.points = FALSE)
```

It is possible to change kernal density estimation plot opacity using the`density.estimation.opacity` argument:

```{r}
map.feature(circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = circassian$language,
            density.width = 0.15,
            density.estimation.opacity = 0.2)
```

If you want to use kernal density estimation, you need to change method type and provide a vector of parameters that increase/decrease area:

```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = "Circassian",
            density.method = "kernal density estimation",
            density.width = c(0.3, 0.3), 
            color= c("darkgreen", "blue"))
```
```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = "Circassian",
            density.method = "kernal density estimation",
            density.width = c(0.7, 0.7), 
            color= c("darkgreen", "blue"))
```
```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = "Circassian",
            density.method = "kernal density estimation",
            density.width = c(1.3, 0.9), 
            color= c("darkgreen", "blue"))
```

It is important to note, that this type of visualization have some shortcomings. The kernel density estimation is calculated without any adjustment, so longitude and latitude values used as a values in Cartesian coordinate system. To reduce consequences of that solution it is better to use a different coordinate projection. That allows not to treat Earth as a flat object.

### 18. Add isoglosses
It is possible to try to catch isoglosses, using the kernel density estimation algorithm. The `map.feature` argument `isogloss` recieves a dataframe with set of features:

```{r}
map.feature(languages = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            features = circassian$dialect,
            label = circassian$dialect,
            legend = TRUE,
            isogloss = as.data.frame(circassian[,"dialect"]),
            isogloss.width = 0.15)
```

It is possible to create true isoglosses by hand, see tools for it [here](https://ropensci.github.io/lingtypology/lingtypology_dplyr.html#2_leaflet).

### 19. Add lines
It is possible to show some lines on the map using coordinates (`line.lng` and `line.lat` arguments).
```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            line.lng = c(39, 43),
            line.lat = c(44.5, 43))
```

If there are more then two coordinates, multiple lines will appear. It is also possible to change the color of the line using the `line.color` argument.

```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            line.lng = c(43, 39, 38.5),
            line.lat = c(43, 44.5, 45),
            line.color= "green")
```

If there are two levels in the `features` variable, it is possible to draw a boundary line between point clusters (the logistic regression is used for calculation).
```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            line.type = "logit")
```

It is possible to add a graticule to a map.
```{r}
map.feature(c("Russian", "Adyghe"),
            graticule = 5)
```

### 20. `ggplot`
Some journals and book publishers are not happy with the resolution of `lingtypology` maps. In order to obtain maps with high resolution in `lingtypology` I need to implement multiple things, and I only started this work. For now only this type of maps are available:

```{r}
ggmap.feature(phonological_profiles$language, phonological_profiles$ejectives)
ggmap.feature(phonological_profiles$language, phonological_profiles$consonants)
```

There will be more functionality in the future.
---
title: "`lingtypology`: introduction and installation"
author: "George Moroz"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{`lingtypology`: introduction and installation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

### What is `lingtypology`?
The `lingtypology` package connects R with the [Glottolog database (v. 4.4)](https://glottolog.org/) and provides an additional functionality for linguistic typology. The Glottolog database contains a catalogue of the world's languages. This package helps researchers to make linguistic maps, using the philosophy of [the Cross-Linguistic Linked Data project](https://clld.org/), which is creating a uniform access to linguistic data across publications. This package is based on the [leaflet package](https://rstudio.github.io/leaflet/), so `lingtypology` is a package for interactive linguistic mapping. In addition, the package provides an ability to download data from typological databases such as WALS, AUTOTYP and others (see section 4). The functionality of this package intersects with the package [`lingtypr`](https://gitlab.com/laurabecker/lingtypr) by Laura Becker. I would like to thank Natalya Tyshkevich, Samira Verhees and Eugenya  Klyagina for reading and correcting some versions of this vignette.

### 1. Installation
Since `lingtypology` is an R package, you should install [R (version >= 3.1.0)](https://www.r-project.org/) on your PC if you haven't already done so. To install the `lingtypology` package, run the following command at your R IDE, so you get the stable version from CRAN:
```{r, eval=FALSE}
install.packages("lingtypology", dependencies = TRUE)
```

You can also get the development version from GitHub:
```{r, eval= FALSE}
install.packages("devtools")
devtools::install_github("ropensci/lingtypology")
```

Load the package:
```{r}
library(lingtypology)
```

### 2. Citing `lingtyplogy`
It is important to cite R and R packages when you use them. For this purpose use the `citation` function:
```{r}
citation("lingtypology")
```
---
title: "`lingtypology` and other packages"
author: "George Moroz"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{`lingtypology` and other packages}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---
```{r, include=FALSE}
library(lingtypology)
knitr::opts_chunk$set(message= FALSE, eval=FALSE)
```

### 1. `dplyr` and `magritr` %>%

It is possible to use [`dplyr`](https://github.com/tidyverse/dplyr) functions and pipes with `lingtypology`. It is widely used, so I will give some examples, how to use it with the`lingtypology` package. Using query "list of languages csv" I found Vincent Garnier's [languages-list repository](https://github.com/forxer/languages-list). Let’s download and map all the languages from that set. First download the data:
 
```{r}
new_data <- read.csv("https://goo.gl/GgscBE")
tail(new_data)
```

As we see, some values of the `Language.name` variable contain more than one language name. Some of the names probably have different names in our database. Imagine that we want to map all languages from Africa. So that the following examples work correctly, use `library(dplyr)`.

```{r, message= FALSE}
library(dplyr)
new_data %>%
  mutate(Language.name = gsub(pattern = " ", replacement = "", Language.name)) %>% 
  filter(is.glottolog(Language.name) == TRUE) %>% 
  filter(area.lang(Language.name) == "Africa") %>% 
  select(Language.name) %>% 
  map.feature()
```

We start with a dataframe, here a `new_data`. First we remove spaces at the end of each string. Then we check, whether the language names are in the glottolog database. Then we select only rows that contain languages of Africa. Then we select the `Language.name` variable. And the last line maps all selected languages.

By default, the values that came from the pipe are treated as the first argument of a function. But when there are some additional arguments, point sign specify what exact position should be piped to. Let’s produce the same map with a minimap.

```{r, message= FALSE}
new_data %>%
  mutate(Language.name = gsub(pattern = " ", replacement = "", Language.name)) %>% 
  filter(is.glottolog(Language.name) == TRUE) %>% 
  filter(area.lang(Language.name) == "Africa") %>% 
  select(Language.name) %>% 
  map.feature(., minimap = TRUE)
```

### 2. `leaflet`, `leaflet.extras`, `mapview`, `mapedit`

There is also a possibility to use `lingtypology` with other [`leaflet`](https://github.com/rstudio/leaflet) functions (thanks to [Niko Partanen](https://github.com/nikopartanen) for the idea):

```{r}
library(leaflet)
map.feature(c("French", "Occitan")) %>% 
  fitBounds(0, 40, 10, 50) %>% 
  addPopups(2, 48, "Great day!")
```

If you add `leaflet` arguments befor `map.feature` function, you need to use argument `pipe.data = .`:

```{r}
leaflet() %>% 
  fitBounds(0, 40, 10, 50) %>% 
  addPopups(2, 48, "Great day!") %>% 
  map.feature(c("French", "Occitan"), pipe.data = .)
```

The other usage of this `pipe.data` argument is to put there a variable with a `leaflet` object:

```{r}
m <- leaflet() %>% 
  fitBounds(0, 40, 10, 50) %>% 
  addPopups(2, 48, "Great day!")

map.feature(c("French", "Occitan"), pipe.data = m)
```

If you want to define tiles in `leaflet` part, you need to change tile argument in `map.feature` function, because the default value for the `tile` argument is "OpenStreetMap.Mapnik".

```{r}
leaflet()  %>% 
  addProviderTiles("Stamen.TonerLite") %>% 
  fitBounds(0, 40, 10, 50) %>% 
  addPopups(2, 48, "Great day!") %>% 
  map.feature(c("French", "Occitan"), pipe.data = ., tile = "none")
```

It is also possible to use some tools provided by [`leaflet.extras` package](https://github.com/RCura/leaflet.extras):

```{r}
map.feature(c("French", "Occitan")) %>% 
  leaflet.extras::addDrawToolbar()  %>%
  leaflet.extras::addStyleEditor()
map.feature(c("French", "Occitan")) %>% 
  leaflet.extras::addFullscreenControl()
```

Also there is a nice package `mapedit` that provide a possibility of creating and editing of leaflet objects by hand:
```{r, eval = FALSE}
map.feature(c("Adyghe", "Russian")) %>% 
  mapedit::editMap() ->
  my_polygone

map.feature(c("Adyghe", "Russian")) %>% 
  leaflet::addPolygons(data = my_polygone$finished)
```
![](use_mapedit.gif)

### 3. Combining maps in a grid and facetisation with `mapview`

The [`leafsync` package](https://github.com/r-spatial/leafsync) provides a possibility to create a multiple maps in a grid and even synchronise them. There are two functions for that: `latticeview()` and `sync()`. Facetisation is a really powerfull tool (look for `facet_grid()` and `facet_wrap()` functions from `ggplot2`). `lingtypology` doesn't provide a facetisation itself, but the `facet` argument of the `map.feature()` function create a list of maps based on this variable. The result of the work of this function then is changed: instead of creating a map in Viewer pane it will return a list that could be used in `latticeview()` and `sync()` functions from the `leafsync` package.

```{r}
faceted <- map.feature(circassian$language,
                       latitude = circassian$latitude,
                       longitude = circassian$longitude,
                       features = circassian$dialect,
                       facet = circassian$language)
library(leafsync)
sync(faceted, no.initial.sync = FALSE)
```

As you can see we provided a `circassian$language` to the `facet` argument, so it returned a list of two maps that stored in `faceted` variable.

It is also possible to combine any maps that were created, just store them in a variable, and combine them in `latticeview()` and `sync()` functions

```{r, eval=FALSE}
m1 <- map.feature(lang.aff("Tsezic"), label = lang.aff("Tsezic"))
m2 <- map.feature(lang.aff("Avar-Andi"), label = lang.aff("Avar-Andi"))
sync(m1, m2)
```

### 4. Get data from OpenStreetMap with `overpass`
This section is inspired by talk with [Niko Partanen](https://github.com/nikopartanen) and his [gist](https://gist.github.com/nikopartanen/f5b4a325808ea8993bfb14b9f81cdfc1). [Overpass](https://github.com/hrbrmstr/overpass) is a packge with tools to work with the OpenStreetMap (OSM) [Overpass API](https://wiki.openstreetmap.org/wiki/Overpass_API). Explore simple Overpass queries with [overpass turbo](https://overpass-turbo.eu/). Imagine that we need to get all settlements from Ingushetia, Daghestan and Chechnya. So, first, load a library:

```{r, eval=FALSE}
library(overpass)
```

Create a query:

```{r, eval=FALSE}
settlements <- 'area[name~"Дагестан|Ингушетия|Чечня"];
(node["place"~"city|village|town|hamlet"](area););
out;'
```

Pass the query to `overpass_query()` function and change the input result to dataframe:
```{r, eval=FALSE}
query_result <- overpass_query(settlements)
settlement_data <- as.data.frame(query_result[, c("id", "lon", "lat", "name")])
```

Some values could be `NA`, so I profer clean it with `complete.cases()` function:
```{r, eval=FALSE}
settlement_data <- settlement_data[complete.cases(settlement_data),]
```

On the last step, I will use a "fake"  language argument to avoid the creation of some Glottolog links:

```{r, eval=FALSE}
map.feature(language = "fake",
            latitude = settlement_data$lat,
            longitude = settlement_data$lon,
            label = settlement_data$name)
```

Results are not ideal: there are some villages Дагестанская and Красный Дагестан in Adygeya and Krasnodarskiy district, but the most points are correct. It is also possible to get all data from some polygone created with `mapedit` (see previous section).

### 5. Create your own atlas with `rmarkdown`
This section is inspired by talk with [Niko Partanen](https://github.com/nikopartanen). It is possible to create an atlas website using `lingtypology` and [`rmarkdown`](https://github.com/rstudio/rmarkdown) packages. The function `atlas.database()` creates a folder in the working directory that contains an `rmarkdown` template for a web-site.

First, lets create a `dataframe` with some data.
```{r}
df <- wals.feature(c("1a", "20a"))
```

Second we can create a website using `atlas.database()` function:

* `languages` argument is a language list
* `features` argument is a data.frame with corresponding features
* `latitude` and `longitude` arguments are optional

```{r}
atlas.database(languages = df$language,
               features = df[,c(4:5)],
               latitude = df$latitude,
               longitude = df$longitude,
               atlas.name = "Some WALS features",
               author = "Author Name")
```

We can see that this function creates a subfolder with following files:
```{r}
list.files("./atlas_Some_WALS_features/")
```

The last step is to run a command:
```{r, eval=FALSE}
rmarkdown::render_site("./atlas_Some_WALS_features/")
```

Then the atlas website will be created (here is [a result](https://agricolamz.github.io/lingtypology_atlas_example/index.html)). If you want to change something in the website, just change some files:

* write information about atlas in index.Rmd file
* list citation information
* change any `.Rmd` file
* ...
* and on the end rerun the `rmarkdown::render_site("./atlas_Some_WALS_features/")` command.

```{r, include=FALSE}
unlink("./atlas_Some_WALS_features/", recursive = TRUE)
```

### 6. Create .kml file using `sp` and `rgdal`
.kml file is a common file type for geospatial data. This kind of files are used in [Google Earth](https://en.wikipedia.org/wiki/Google_Earth), [Gabmap](https://gabmap.nl/) (a web application that visualizes dialect variations) and others. In order to produce a .kml file you need to have a dataset with coordinates such as `circassian`:
```{r, eval= FALSE}
sp::coordinates(circassian) <- ~longitude+latitude
sp::proj4string(circassian) <- sp::CRS("+proj=longlat +datum=WGS84")
rgdal::writeOGR(circassian["village"],
                "circassian.kml", 
                layer="village",
                driver="KML")
```
---
title: "`lingtypology`: Glottolog functions"
author: "George Moroz"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{`lingtypology`: Glottolog functions}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r, include=FALSE}
library(lingtypology)
```

This package is based on the [Glottolog database (v. 4.4)](https://glottolog.org/), so `lingtypology` has several functions for accessing data from that database.

### 1. Command name's syntax
Most of the functions in `lingtypology` have the same syntax: **what you need.what you have**. Most of them are based on _language name_.

* **aff.lang()** — get affiliation by language
* **area.lang()** — get macro area by language
* **country.lang()** — get country by language
* **iso.lang()** — get ISO 639-3 code by language
* **gltc.lang()** — get glottocode (identifier for a language in the Glottolog database) code by language
* **lat.lang()** — get latitude by language
* **long.lang()** — get longitude by language
* **level.lang()** — get level by language
* **subc.lang()** — get subclassification in the Newick tree format by language

Some of them help to define a vector of languages.

* **lang.aff()** — get language by affiliation
* **lang.iso()** — get language by ISO 639-3 code
* **lang.gltc()** — get language by glottocode

Additionally there are some functions to convert glottocodes to ISO 639-3 codes and vice versa:

* **gltc.iso()** — get glottocode by ISO 639-3 code
* **iso.gltc()** — get ISO 639-3 code by glottocode

The most important functionality of `lingtypology` is the ability to create interactive maps based on features and sets of languages (see the third section):

* **map.feature()**

[Glottolog database (v. 4.1)](https://glottolog.org/) provides `lingtypology` with language names, ISO codes, glottocodes, affiliation, macro area, coordinates, and much information. This set of functions doesn't have a goal to cover all possible combinations of functions. Check out additional information that is preserved in the version of the Glottolog database used in `lingtypology`:

```{r}
names(glottolog)
```

Using R functions for data manipulation you can create your own database for your purpose.

### 2. Using base functions
All functions introduced in the previous section are regular functions, so they can take the following objects as input:

* a regular string
```{r}
iso.lang("Adyghe")
lang.iso("ady")
lang.aff("West Caucasian")
```

I would like to point out that you can create strings in R using single or double quotes. Since inserting single quotes in a string created with single quotes causes an error in R, I use double quotes in my tutorial. You can use single quotes, but be careful and remember that `'Ma'ya'` is an incorrect string in R.

* a vector of strings
```{r}
area.lang(c("Adyghe", "Aduge"))
lang <- c("Adyghe", "Russian")
aff.lang(lang)
```

*  other functions. For example, let's try to get a vector of ISO codes for the Circassian languages
```{r}
iso.lang(lang.aff("Circassian"))
```

If you are new to R, it is important to mention that you can create a table with languages, features and other parametres with any spreadsheet software you used to work. Then you can import the created file to R using standard tools.

### 3. Spell Checker: look carefully at warnings!

All functions which take a vector of languages are enriched with a kind of a spell checker. If a language from a query is absent in the database, functions return a warning message containing a set of candidates with the minimal Levenshtein distance to the language from the query.

```{r}
aff.lang("Adyge")
```

### 4. `subc.lang()` function

The `subc.lang()` function returns language subclassification in the Newick tree format.

```{r}
subc.lang("Lechitic")
```

This format is hard to interpret by itself, but there are some tools in R that make it possible to visualise those subclassifications:

```{r}
library(ape)
plot(read.tree(text = subc.lang("Lechitic")))
```

It is possible to specify colors of tips in case you want to emphasize some nodes:

```{r}
plot(read.tree(text = subc.lang("Lechitic")),
     tip.color = c("red", "black", "black", "black"))
```

As you can see nodes are counted from bottom to top.

For more sophisticated tree visualization you can look into [`ggtree` package](https://bioconductor.org/packages/release/bioc/html/ggtree.html). There are several linguistic packages that provide some functionality for creating glottolog trees:

* [`glottoTrees`](https://github.com/erichround/glottoTrees) package by Erich Round
* [`lingtypr`](https://gitlab.com/laurabecker/lingtypr) package by Laura Becker
---
title: "`lingtypology`: Typological databases API"
author: "[George Moroz](mailto:agricolamz@gmail.com)"
editor_options: 
  chunk_output_type: console
---

```{r, include=FALSE}
library(lingtypology)
knitr::opts_chunk$set(fig.width=7.5, fig.height=3.4, comment = "")
```

`lingtypology` provides an ability to download data from these typological databases

* [World Atlas of Language Structures](http://wals.info/)
* [AUTOTYP](https://github.com/autotyp/autotyp-data#the-autotyp-database)
* [PHOIBLE](http://phoible.org/)
* [Affix Borrowing database](http://afbo.info)
* [South American Indigenous Language Structures](http://sails.clld.org/)
* [Austronesian Basic Vocabulary Database](https://abvd.shh.mpg.de/austronesian/)
* [UraLex](https://github.com/lexibank/uralex)

All database function names have identical structure: **database_name.feature**. All functions have as first argument `feature`. All functions create dataframe with column `language` that can be used in `map.feature()` function. It should be noted that all functions cut out the data that can't be maped, so if you want to prevent functions from this behaviour set argument `na.rm` to `FALSE`.

### 1. WALS
The names of the WALS features can be typed in a lower case. This function preserves coordinates from WALS, so you can map coordinates from the WALS or use coordinates from `lingtypology`.
```{r}
df <- wals.feature(c("1a", "20a"))
head(df)
map.feature(df$language,
            features = df$`1a`,
            latitude = df$latitude,
            longitude = df$longitude,
            label = df$language,
            title = "Consonant Inventories")
```

### 2. AUTOTYP
The AUTOTYP features are listed on [the GitHub page](https://github.com/autotyp/autotyp-data#the-autotyp-database). You can use more human way with spaces.
```{r}
df <- autotyp.feature(c('Gender', 'Numeral classifiers'))
head(df)
map.feature(df$language,
            features = df$NumClass.Presence,
            label = df$language,
            title = "Presence of Numeral Classifiers")
```

### 3. PHOIBLE
I used only four features from PHOIBLE: the number of phonemes, the number of consonants, the number of tones and the number of vowels. If you need only a set of them, just specify it in the `features` argument. Since there is a lot of doubling information in the PHOIBLE database, there is an argument `source`.
```{r}
df <- phoible.feature()
head(df)
```

### 4. AfBo
The AfBo database has a lot of features that distinguish affix functions, but again you can use a bare function without any arguments to download the whole database. There will be no difference in time, since this function downloads the whole database to your PC. The main destinction is that this database provides recipient and donor languages, so other column names should be used.

```{r}
df <- afbo.feature(c("adjectivizer", "adverbializer"))
head(df)
map.feature(df$Recipient.name,
            features = df$adjectivizer,
            label = df$Recipient.name,
            title = "Borrowed adjectivizer affixes")
```

### 5. SAILS
The SAILS database provide a lot of [features](http://sails.clld.org/parameters), so the function work with their ids:
```{r}
df <- sails.feature(features = "ics10")
head(df)
map.feature(df$language,
            features = df$ics10_description,
            longitude = df$longitude,
            latitude = df$latitude,
            label = df$language,
            title = "Are there numeral classifiers?")
```

### 6. ABVD
The ABVD database is a lexical database, so it is different from clld databases. First of all, ABVD has its own language classification ids. The information about the same language from different sources can be received from these database different ids. So I select several languages and map them coloring by word with the meaning 'hand'.
```{r}
df <- abvd.feature(c(292, 7))
head(df)
new_df <- df[df$word == "hand",]
map.feature(new_df$language,
            features = new_df$item,
            label = new_df$language)
```

### 7. UraLex

`uralex.feature` downloads data from UraLex basic vocabulary dataset. Original language names are stored in the `language` variable. Converted language names for `map.feature` are stored in the `language2` variable.

```{r}
df <- uralex.feature()
df <- df[df$uralex_mng == "crush",]

map.feature(df$language2,
            label = df$item,
            title = "crush")
```

---
title: "`lingtypology`: creating maps"
author: "[George Moroz](mailto:agricolamz@gmail.com)"
editor_options: 
  chunk_output_type: console
  output: html_document
---

```{r, include=FALSE}
library(lingtypology)
knitr::opts_chunk$set(fig.width=7.5, fig.height=3.4)
```


### 1. Base map
The most important part of the `lingtypology` package is the function `map.feature`. This function allows you to produce maps similar to known projects within [the Cross-Linguistic Linked Data philosophy](http://clld.org/), such as [WALS](http://wals.info/) and [Glottolog](http://glottolog.org/):
```{r}
map.feature(c("Adyghe", "Kabardian", "Polish", "Russian", "Bulgarian"))
```

As shown in the picture above, this function generates an interactive Leaflet map. All specific points on the map have a pop-up box that appears when markers are clicked (see section 3.3 for more information about editing pop-up boxes). By default, they contain language names linked to the glottolog site.

If for some reasons you are not using RStudio or you want to automatically create and save a lot of maps, you can save a map to a variable and use the `htmlwidgets` package for saving created maps to an .html file. I would like to thank Timo Roettger for mentioning this problem.

```{r, eval = FALSE}
m <- map.feature(c("Adyghe", "Korean"))
# install.packages("htmlwidgets")
library(htmlwidgets)
saveWidget(m, file="TYPE_FILE_PATH/m.html")
```

There is an export button in RStudio, but for some reason it is not so easy to save a map as a .png or.jpg file using code. [Here](http://stackoverflow.com/a/34672309/6056442) is a possible solution.

### 2. Set features
The goal of this package is to allow typologists (or any other linguists) to map language features. A list of languages and correspondent features can be stored in a `data.frame` as follows:
```{r}
df <- data.frame(language = c("Adyghe", "Kabardian", "Polish", "Russian", "Bulgarian"),
               features = c("polysynthetic", "polysynthetic", "fusional", "fusional", "fusional"))
df
```

Now we can draw a map:
```{r}
map.feature(languages = df$language,
            features = df$features)
```

If you have a lot of features and they appear in the legend in a senseless order(by default it is ordered alphabetically), you can reorderthem using factors (a vector with ordered levels, for more information see `?factor`). for example, I want the feature polysynthetic to be listed first, followed by fusional:

```{r}
df$features <- factor(df$features, levels = c("polysynthetic", "fusional"))
map.feature(languages = df$language, features = df$features)
```

Like in most functions, it is not necessary to name all arguments, so the same result can be obtained by:
```{r}
map.feature(df$language, df$features)
```

As shown in the picture above, all points are grouped by feature, colored and counted. As before, a pop-up box appears when markers are clicked. A control feature allows users to toggle the visibility of points grouped by feature.

There are several types of variables in R and `map.feature` works differently depending on the variable type. I will use a build in data set `phonological_profiles` that contains 19 languages from [UPSyD database](https://lapsyd.huma-num.fr/lapsyd/). This dataset have three variables: the categorical variable `ejectives` indicates whether some language has any ejective sound, numeric variables `consonants` and `vowels` that contains information about the number of consonants and vowels (based on UPSyD database). We can create two maps with categorical variable and with numeric variable:

```{r}
map.feature(phonological_profiles$language,
            phonological_profiles$ejectives) # categorical
map.feature(phonological_profiles$language,
            phonological_profiles$consonants) # numeric
```

The main point is that for creating a correct map, you should correctly define the type of the variable.

This dataset also can be used to show one other parameter of the `map.feature` function. There are two possible ways to show the World map: with the Atlantic sea or with the Pacific sea in the middle. If you don't need default Pacific view use the `map.orientation` parameter(thanks @languageSpaceLabs and @tzakharko for that idea):
```{r}
map.feature(phonological_profiles$language,
            phonological_profiles$consonants,
            map.orientation = "Atlantic")
```

Sometimes the result image legend are too small. It is possible to increase size of legend using some trick described here. In order to do this assign your map to a variable:

```{r}
m <- map.feature(languages = df$language,
                 features = df$features)
```

Map will not appear, but now you can get it calling the variable `m`. After this in order to increase the legend font one can use code like this:

```{r}
library(htmltools)
browsable(
  tagList(
    list(
      tags$head(
        tags$style(
          ".leaflet .legend {
                 line-height: 20px;
                 font-size: 20px;
                 }",
          ".leaflet .legend i{
                width: 20px;
                height: 20px;
                 }"
        )
      ),
      m)))
```


### 3. Set pop-up boxes
Sometimes it is a good idea to add some additional information (e.g. language affiliation, references or even examples) to pop-up boxes that appear when points are clicked. In order to do so, first of all we need to create an extra vector of strings in our dataframe:
```{r}
df$popup <- aff.lang(df$language)
```

The function `aff.lang()` creates a vector of genealogical affiliations that can be easily mapped:
```{r}
map.feature(languages = df$language, features = df$features, popup = df$popup)
```

Pop-up strings can contain HTML tags, so it is easy to insert a link, a couple of lines, a table or even a video and sound. Here is how pop-up boxes can demonstrate language examples:
```{r}
# change a df$popup vector
df$popup <- c("sɐ s-ɐ-k'ʷɐ<br> 1sg 1sg.abs-dyn-go<br>'I go'",
              "sɐ s-o-k'ʷɐ<br> 1sg 1sg.abs-dyn-go<br>'I go'",
              "id-ę<br> go-1sg.npst<br> 'I go'",
              "ya id-u<br> 1sg go-1sg.npst <br> 'I go'",
              "id-a<br> go-1sg.prs<br> 'I go'")
# create a map

map.feature(df$language,
            features = df$features,
            popup = df$popup)
```

How to say _moon_ in Sign Languages? Here is an example:
```{r}
# Create a dataframe with links to video
sign_df <- data.frame(languages = c("American Sign Language", "Russian-Tajik Sign Language", "French Sign Language"),
                      popup = c("https://media.spreadthesign.com/video/mp4/13/48600.mp4", "https://media.spreadthesign.com/video/mp4/12/17639.mp4", "https://media.spreadthesign.com/video/mp4/10/17638.mp4"))

# Change popup to an HTML code
sign_df$popup <- paste("<video width='200' height='150' controls> <source src='",
                       as.character(sign_df$popup),
                       "' type='video/mp4'></video>", sep = "")
# create a map
map.feature(languages = sign_df$languages, popup = sign_df$popup)
```

### 4. Set labels
An alternative way to add some short text to a map is to use the `label` option.
```{r}
map.feature(df$language, df$features,
            label = df$language)
```

There are some additional arguments for customization: `label.fsize` for setting font size, `label.position` for controlling the label position, and `label.hide` to control the appearance of the label: if `TRUE`, the labels are displayed on mouse over(as on the previous map), if `FALSE`, the labels are always displayed (as on the next map).

```{r}
map.feature(df$language, df$features,
            label = df$language,
            label.fsize = 20,
            label.position = "left",
            label.hide = FALSE)
```

There is an additional tool for emphasis of some points on the map. The argument `label.emphasize` allows to emphasize selected points with the color specified by a user.

```{r}
map.feature(df$language, df$features,
            label = df$language,
            label.fsize = 20,
            label.position = "left",
            label.hide = FALSE,
            label.emphasize = list(2:4, "red"))
```

In this example the first vector of the list in the `label.emphasize` argument is vector `2:4` that produce elements `2`, `3` and `4`. You can create youro wn selected rows. e. g. `c(1, 3, 4)`. The second vector of the list is the string with a color.

### 5. Set coordinates
You can set your own coordinates using the arguments `latitude` and `longitude`. It is important to note, that `lingtypology` works only with decimal degrees (something like this: 0.1), not with degrees, minutes and seconds (something like this: 0° 06′ 0″). I will illustrate this with the dataset `circassian` built into the `lingtypology` package. This dataset comes from fieldwork collected during several expeditions in the period 2011-2016 and contains a list of Circassian villages:
```{r}
head(circassian)
```

In this dataframe you can find variables `latitude` and `longitude` that could be used:
```{r}
map.feature(languages = circassian$language,
            features = circassian$dialect,
            popup = circassian$village,
            latitude = circassian$latitude,
            longitude = circassian$longitude)
```

It is possible to collapse multiple dots into clusters:

```{r}
map.feature(languages = circassian$language,
            features = circassian$dialect,
            popup = circassian$village,
            latitude = circassian$latitude,
            longitude = circassian$longitude, 
            point.cluster = TRUE)
```

### 6. Set colors
You can set your own colors using the argument `color`:
```{r}
df <- data.frame(language = c("Adyghe", "Kabardian", "Polish", "Russian", "Bulgarian"),
                 features = c("polysynthetic", "polysynthetic", "fusional", "fusional", "fusional"))
map.feature(languages = df$language,
            features = df$features,
            color= c("yellowgreen", "navy"))
```

Arguments from [RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer) or [viridis](https://CRAN.R-project.org/package=viridis) also can be used as a color argument:

```{r}
map.feature(phonological_profiles$language,
            phonological_profiles$consonants,
            color= "magma")
```

### 7. Set shapes
For some scientific papers it is not possible to use colors for destinguishing features. In that cases it is posible to use `shape` argument:

```{r}
map.feature(languages = circassian$language,
            features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            shape = TRUE)
```

The argument `shape = TRUE` works fine only with 6 or less levels in `features` variable. If there are more levels in `fetures` argument, user need to provide a vector with corresponding shapes:

```{r}
map.feature(languages = circassian$language,
            features = circassian$dialect,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            shape = 1:10,
            shape.size = 14)
```

Arguments `shape.size` and `shape.color` help to change corresponding features of markers.

### 8. Set control box
The package can generate a control box that allows users to toggle the visibility of some points. To enable it, there is an argument `control` in the `map.feature` function:

```{r}
map.feature(languages = df$language,
            features = df$features,
            control = c("a", "b", "b", "b", "a"))
```

As you can see the `control` and `features` arguments are independent of each other.

### 9. Set an additional set of features using strokes
The `map.feature` function has an additional argument `stroke.features`. Using this argument it becomes possible to show two independent sets of features on one map. By default strokes are colored in grey (so for two levels it will be black and white, for three --- black, grey, white and so on), but you can set your own colors using the argument `stroke.color`:
```{r}
map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude)
```

It is important to note that `stroke.features` can work with `NA` values. The function won't plot anything if there is an `NA` value. Let's set a language value to `NA` in all Baksan villages from the `circassian` dataset.

```{r, message= F}
# create newfeature variable
newfeature <- circassian[,c(5,6)]
# set language feature of the Baksan villages to NA and reduce newfeature from dataframe to vector
newfeature <- replace(newfeature$language, newfeature$language == "Baksan", NA)
# create a map

map.feature(circassian$language,
            features = circassian$dialect,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            stroke.features = newfeature)
```

### 10. Set width and an opacity feature
All markers have their own width and opacity, so you can set it. Just use the arguments `width`, `stroke.radius`, `opacity` and `stroke.opacity`:
```{r}
map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            width = 7, stroke.radius = 13)

map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            opacity = 0.7, stroke.opacity = 0.6)
```

### 11. Customizing legends
By default the legend appears in the top right corner. If there are stroke features, two legends are generated. There are additional arguments that control the appearence and the title of the legends.
```{r}
map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            legend = FALSE, stroke.legend = TRUE)

map.feature(circassian$language,
            features = circassian$dialect,
            stroke.features = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            title = "Circassian dialects", stroke.title = "Languages")
```
The arguments `legend.position` and `stroke.legend.position` allow you to change a legend's position using "topright", "bottomright", "bottomleft" or"topleft" strings.

### 12. Set scale bar
A scale bar is automatically added to a map, but you can control its appearance (set `scale.bar` argument to `TRUE` or`FALSE`) and its position (use `scale.bar.position` argument values "topright", "bottomright", "bottomleft" or"topleft").
```{r}
map.feature(c("Adyghe", "Polish", "Kabardian", "Russian"),
            scale.bar= TRUE,
            scale.bar.position = "topright")
```

### 13. Set layouts
It is possible to use different tiles on the same map using the `tile` argument. For more tiles see [here](https://leaflet-extras.github.io/leaflet-providers/preview/index.html).
```{r}
df <- data.frame(lang = c("Adyghe", "Kabardian", "Polish", "Russian", "Bulgarian"),
                 feature = c("polysynthetic", "polysynthetic", "fusion", "fusion", "fusion"),
                 popup = c("Adyghe", "Adyghe", "Slavic", "Slavic", "Slavic"))

map.feature(df$lang, df$feature, df$popup,
            tile = "Stamen.TonerLite")
```

It is possible to use different map tiles on the same map. Just add a vector with tiles.
```{r}
map.feature(df$lang, df$feature, df$popup,
            tile = c("OpenStreetMap", "Stamen.TonerLite"))
```

It is possible to name tiles using the `tile.name` argument.
```{r}
map.feature(df$lang, df$feature, df$popup,
            tile = c("OpenStreetMap", "Stamen.TonerLite"),
            tile.name = c("colored", "b & w"))
```

It is possible to combine the tiles' control box with the features' control box.
```{r}
map.feature(df$lang, df$feature, df$popup,
            tile = c("OpenStreetMap", "Stamen.TonerLite"),
            control = TRUE)
```

### 14. Add a minimap to a map
It is possible to add a minimap to a map.
```{r}
map.feature(c("Adyghe", "Polish", "Kabardian", "Russian"),
            minimap = TRUE)
```

You can control its appearance (by setting the `minimap` argument to `TRUE` or `FALSE`), its position (by using the values "topright", "bottomright", "bottomleft" or"topleft" of the `minimap.position` argument) and its height and width (with the arguments `minimap.height` and `minimap.width`).
```{r}
map.feature(c("Adyghe", "Polish", "Kabardian", "Russian"),
            minimap = TRUE,
            minimap.position = "topright",
            minimap.height = 100,
            minimap.width = 100)
```

### 15. Add minicharts instead of points
This part is created using the beutifull [`leaflet.minicharts` library](https://CRAN.R-project.org/package=leaflet.minicharts). The argument `minichart` allows you to add piecharts or barplots instead of standard point markers. In this part I will use a build in data set `phonological_profiles` that contains 19 languages from [UPSyD database](https://lapsyd.huma-num.fr/lapsyd/). Here is an example of barplot:

```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("vowels", "consonants")])
```

Here is an example of piechart:
```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("vowels", "consonants")],
            minichart = "pie")
```

Colors and opacity could be changed, legend moved:
```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("vowels", "consonants")],
            color= c("yellowgreen", "navy"),
            opacity = 0.7,
            label = phonological_profiles$language,
            legend.position = "topleft")
```

It is possible to add values using argument `minichart.labels`:
```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("vowels", "consonants")],
            minichart = "pie",
            minichart.labels = TRUE)
```

It is also possible to use pie chart in non-convenient way: just indicating with `TRUE` or `FALSE` of pressence of some feature (thanks to Diana Forker for the task!):

```{r}
map.feature(languages = phonological_profiles$language,
            minichart.data = phonological_profiles[, c("tone", "long_vowels", "stress", "ejectives")],
            minichart = "pie", 
            width = 3)
```

Unfortunately this kind of visualisation doesn't work, when you have some lines in your dataset that contain only `FALSE` values. This is **non-convenient** way of category visualisation, so visualisation experts could have a negative opinion about it. This kind of visualisation is also bad for huge number of variables.

### 16. Add a rectangle to a map
It is possible to highlight some part of your map with a rectangle. You need to provide a latitude and longitude of the diagonal (`rectangle.lat` and `rectangel.lng`) and color of the rectangle (`rectangle.color`):

```{r}
map.feature(circassian$language,
            circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            rectangle.lng = c(42.7, 45),
            rectangle.lat = c(42.7, 44.4),
            rectangle.color= "green")
```

### 17. Add a density contourplot to a map
Sometimes it is easier to look at a density contourplot. It can be created using `density.estimation` argument. There are two possibility for creation a density contourplot in `lingtypology`:

* `density.method = "fixed distance"`. First algorithm creates circle polygons with fixed radius around each point and then merge all polygons that are overlapped. It has only one parameter that should be estimated: radius of the circle (`density.width`).
* `density.method = "kernal density estimation"`. Second algorithm uses a kernal density estimation and has two parameters that should be estimated: latitude and longitude bandwidths (`density.width[1]` and (`density.width[2]`))


```{r}
map.feature(circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = circassian$language,
            density.width = 0.15)
```

Density estimation plot can be separated by `features` variable:
```{r}
map.feature(circassian$language,
            features = circassian$dialect,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = circassian$language,
            density.width = 0.15)
```

It is possible to remove points and display only the kernal density estimation plot, using the `density.points` argument:

```{r}
map.feature(circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = circassian$language,
            density.width = 0.15,
            density.points = FALSE)
```

It is possible to change kernal density estimation plot opacity using the`density.estimation.opacity` argument:

```{r}
map.feature(circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = circassian$language,
            density.width = 0.15,
            density.estimation.opacity = 0.2)
```

If you want to use kernal density estimation, you need to change method type and provide a vector of parameters that increase/decrease area:

```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = "Circassian",
            density.method = "kernal density estimation",
            density.width = c(0.3, 0.3), 
            color= c("darkgreen", "blue"))
```
```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = "Circassian",
            density.method = "kernal density estimation",
            density.width = c(0.7, 0.7), 
            color= c("darkgreen", "blue"))
```
```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            density.estimation = "Circassian",
            density.method = "kernal density estimation",
            density.width = c(1.3, 0.9), 
            color= c("darkgreen", "blue"))
```

It is important to note, that this type of visualization have some shortcomings. The kernel density estimation is calculated without any adjustment, so longitude and latitude values used as a values in Cartesian coordinate system. To reduce consequences of that solution it is better to use a different coordinate projection. That allows not to treat Earth as a flat object.

### 18. Add isoglosses
It is possible to try to catch isoglosses, using the kernel density estimation algorithm. The `map.feature` argument `isogloss` recieves a dataframe with set of features:

```{r}
map.feature(languages = circassian$language,
            latitude = circassian$latitude,
            longitude = circassian$longitude,
            features = circassian$dialect,
            label = circassian$dialect,
            legend = TRUE,
            isogloss = as.data.frame(circassian[,"dialect"]),
            isogloss.width = 0.15)
```

It is possible to create true isoglosses by hand, see tools for it [here](https://ropensci.github.io/lingtypology/lingtypology_dplyr.html#2_leaflet).

### 19. Add lines
It is possible to show some lines on the map using coordinates (`line.lng` and `line.lat` arguments).
```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            line.lng = c(39, 43),
            line.lat = c(44.5, 43))
```

If there are more then two coordinates, multiple lines will appear. It is also possible to change the color of the line using the `line.color` argument.

```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            line.lng = c(43, 39, 38.5),
            line.lat = c(43, 44.5, 45),
            line.color= "green")
```

If there are two levels in the `features` variable, it is possible to draw a boundary line between point clusters (the logistic regression is used for calculation).
```{r}
map.feature(circassian$language,
            features = circassian$language,
            longitude = circassian$longitude,
            latitude = circassian$latitude,
            line.type = "logit")
```

It is possible to add a graticule to a map.
```{r}
map.feature(c("Russian", "Adyghe"),
            graticule = 5)
```


### 20. `ggplot`
Some journals and book publishers are not happy with the resolution of `lingtypology` maps. In order to obtain maps with high resolution in `lingtypology` I need to implement multiple things, and I only started this work. For now only this type of maps are available:

```{r}
ggmap.feature(phonological_profiles$language, phonological_profiles$ejectives)
ggmap.feature(phonological_profiles$language, phonological_profiles$consonants)
```

There will be more functionality in the future.
---
title: "`lingtypology` and other packages"
author: "[George Moroz](mailto:agricolamz@gmail.com)"
editor_options: 
  chunk_output_type: console
---

```{r, include=FALSE}
library(lingtypology)
knitr::opts_chunk$set(fig.width=7.5, fig.height=3.4)
```

### 1. `dplyr` and `magritr` %>%

It is possible to use [`dplyr`](https://github.com/tidyverse/dplyr) functions and pipes with `lingtypology`. It is widely used, so I will give some examples, how to use it with the`lingtypology` package. Using query "list of languages csv" I found Vincent Garnier's [languages-list repository](https://github.com/forxer/languages-list). Let’s download and map all the languages from that set. First download the data:
 
```{r}
new_data <- read.csv("https://goo.gl/GgscBE")
tail(new_data)
```

As we see, some values of the `Language.name` variable contain more than one language name. Some of the names probably have different names in our database. Imagine that we want to map all languages from Africa. So that the following examples work correctly, use `library(dplyr)`.

```{r, message= FALSE}
library(dplyr)
new_data %>%
  mutate(Language.name = gsub(pattern = " ", replacement = "", Language.name)) %>% 
  filter(is.glottolog(Language.name) == TRUE) %>% 
  filter(area.lang(Language.name) == "Africa") %>% 
  select(Language.name) %>% 
  map.feature()
```

We start with a dataframe, here a `new_data`. First we remove spaces at the end of each string. Then we check, whether the language names are in the glottolog database. Then we select only rows that contain languages of Africa. Then we select the `Language.name` variable. And the last line maps all selected languages.

By default, the values that came from the pipe are treated as the first argument of a function. But when there are some additional arguments, point sign specify what exact position should be piped to. Let’s produce the same map with a minimap.

```{r, message= FALSE}
new_data %>%
  mutate(Language.name = gsub(pattern = " ", replacement = "", Language.name)) %>% 
  filter(is.glottolog(Language.name) == TRUE) %>% 
  filter(area.lang(Language.name) == "Africa") %>% 
  select(Language.name) %>% 
  map.feature(., minimap = TRUE)
```

### 2. `leaflet`, `leaflet.extras`, `mapview`, `mapedit`

There is also a possibility to use `lingtypology` with other [`leaflet`](https://github.com/rstudio/leaflet) functions (thanks to [Niko Partanen](https://github.com/nikopartanen) for the idea):

```{r}
library(leaflet)
map.feature(c("French", "Occitan")) %>% 
  fitBounds(0, 40, 10, 50) %>% 
  addPopups(2, 48, "Great day!")
```

If you add `leaflet` arguments befor `map.feature` function, you need to use argument `pipe.data = .`:

```{r}
leaflet() %>% 
  fitBounds(0, 40, 10, 50) %>% 
  addPopups(2, 48, "Great day!") %>% 
  map.feature(c("French", "Occitan"), pipe.data = .)
```

The other usage of this `pipe.data` argument is to put there a variable with a `leaflet` object:

```{r}
m <- leaflet() %>% 
  fitBounds(0, 40, 10, 50) %>% 
  addPopups(2, 48, "Great day!")

map.feature(c("French", "Occitan"), pipe.data = m)
```

If you want to define tiles in `leaflet` part, you need to change tile argument in `map.feature` function, because the default value for the `tile` argument is "OpenStreetMap.Mapnik".

```{r}
leaflet()  %>% 
  addProviderTiles("Stamen.TonerLite") %>% 
  fitBounds(0, 40, 10, 50) %>% 
  addPopups(2, 48, "Great day!") %>% 
  map.feature(c("French", "Occitan"), pipe.data = ., tile = "none")
```

It is also possible to use some tools provided by [`leaflet.extras` package](https://github.com/RCura/leaflet.extras):

```{r}
map.feature(c("French", "Occitan")) %>% 
  leaflet.extras::addDrawToolbar()  %>%
  leaflet.extras::addStyleEditor()
map.feature(c("French", "Occitan")) %>% 
  leaflet.extras::addFullscreenControl()
```

Also there is a nice package `mapedit` that provide a possibility of creating and editing of leaflet objects by hand:
```{r, eval = FALSE}
map.feature(c("Adyghe", "Russian")) %>% 
  mapedit::editMap() ->
  my_polygone

map.feature(c("Adyghe", "Russian")) %>% 
  leaflet::addPolygons(data = my_polygone$finished)
```
![](https://raw.githubusercontent.com/ropensci/lingtypology/master/docs/use_mapedit.gif)

### 3. Combining maps in a grid and facetisation with `mapview`

The [`leafsync` package](https://github.com/r-spatial/leafsync) provides a possibility to create a multiple maps in a grid and even synchronise them. There are two functions for that: `latticeview()` and `sync()`. Facetisation is a really powerfull tool (look for `facet_grid()` and `facet_wrap()` functions from `ggplot2`). `lingtypology` doesn't provide a facetisation itself, but the `facet` argument of the `map.feature()` function create a list of maps based on this variable. The result of the work of this function then is changed: instead of creating a map in Viewer pane it will return a list that could be used in `latticeview()` and `sync()` functions from the `leafsync` package.

```{r}
faceted <- map.feature(circassian$language,
                       latitude = circassian$latitude,
                       longitude = circassian$longitude,
                       features = circassian$dialect,
                       facet = circassian$language)
library(leafsync)
sync(faceted, no.initial.sync = FALSE)
```

As you can see we provided a `circassian$language` to the `facet` argument, so it returned a list of two maps that stored in `faceted` variable.

It is also possible to combine any maps that were created, just store them in a variable, and combine them in `latticeview()` and `sync()` functions

```{r, eval=FALSE}
m1 <- map.feature(lang.aff("Tsezic"), label = lang.aff("Tsezic"))
m2 <- map.feature(lang.aff("Avar-Andi"), label = lang.aff("Avar-Andi"))
sync(m1, m2)
```

### 4. Get data from OpenStreetMap with `overpass`
This section is inspired by talk with [Niko Partanen](https://github.com/nikopartanen) and his [gist](https://gist.github.com/nikopartanen/f5b4a325808ea8993bfb14b9f81cdfc1). [Overpass](https://github.com/hrbrmstr/overpass) is a packge with tools to work with the OpenStreetMap (OSM) [Overpass API](http://wiki.openstreetmap.org/wiki/Overpass_API). Explore simple Overpass queries with [overpass turbo](http://overpass-turbo.eu/). Imagine that we need to get all settlements from Ingushetia, Daghestan and Chechnya. So, first, load a library:

```{r, eval=FALSE}
library(overpass)
```

Create a query:

```{r, eval=FALSE}
settlements <- 'area[name~"Дагестан|Ингушетия|Чечня"];
(node["place"~"city|village|town|hamlet"](area););
out;'
```

Pass the query to `overpass_query()` function and change the input result to dataframe:
```{r, eval=FALSE}
query_result <- overpass_query(settlements)
settlement_data <- as.data.frame(query_result[, c("id", "lon", "lat", "name")])
```

Some values could be `NA`, so I profer clean it with `complete.cases()` function:
```{r, eval=FALSE}
settlement_data <- settlement_data[complete.cases(settlement_data),]
```

On the last step, I will use a "fake"  language argument to avoid the creation of some Glottolog links:

```{r, eval=FALSE}
map.feature(language = "fake",
            latitude = settlement_data$lat,
            longitude = settlement_data$lon,
            label = settlement_data$name)
```

Results are not ideal: there are some villages Дагестанская and Красный Дагестан in Adygeya and Krasnodarskiy district, but the most points are correct. It is also possible to get all data from some polygone created with `mapedit` (see previous section).

### 5. Create your own atlas with `rmarkdown`
This section is inspired by talk with [Niko Partanen](https://github.com/nikopartanen). It is possible to create an atlas website using `lingtypology` and [`rmarkdown`](https://github.com/rstudio/rmarkdown) packages. The function `atlas.database()` creates a folder in the working directory that contains an `rmarkdown` template for a web-site.

First, lets create a `dataframe` with some data.
```{r}
df <- wals.feature(c("1a", "20a"))
```

Second we can create a website using `atlas.database()` function:

* `languages` argument is a language list
* `features` argument is a data.frame with corresponding features
* `latitude` and `longitude` arguments are optional

```{r}
atlas.database(languages = df$language,
               features = df[,c(4:5)],
               latitude = df$latitude,
               longitude = df$longitude,
               atlas.name = "Some WALS features",
               author = "Author Name")
```

We can see that this function creates a subfolder with following files:
```{r}
list.files("./atlas_Some_WALS_features/")
```

The last step is to run a command:
```{r, eval=FALSE}
rmarkdown::render_site("./atlas_Some_WALS_features/")
```

Then the atlas website will be created (here is [a result](https://agricolamz.github.io/lingtypology_atlas_example/index.html)). If you want to change something in the website, just change some files:

* write information about atlas in index.Rmd file
* list citation information
* change any `.Rmd` file
* ...
* and on the end rerun the `rmarkdown::render_site("./atlas_Some_WALS_features/")` command.

```{r, include=FALSE}
unlink("./atlas_Some_WALS_features/", recursive = TRUE)
```

### 6. Create .kml file using `sp` and `rgdal`
.kml file is a common file type for geospatial data. This kind of files are used in [Google Earth](https://en.wikipedia.org/wiki/Google_Earth), [Gabmap](http://www.gabmap.nl/) (a web application that visualizes dialect variations) and others. In order to produce a .kml file you need to have a dataset with coordinates such as `circassian`:
```{r, eval= FALSE}
sp::coordinates(circassian) <- ~longitude+latitude
sp::proj4string(circassian) <- sp::CRS("+proj=longlat +datum=WGS84")
rgdal::writeOGR(circassian["village"],
                "circassian.kml", 
                layer="village",
                driver="KML")
```
---
title: "`lingtypology`: Glottolog functions"
author: "[George Moroz](mailto:agricolamz@gmail.com)"
---

```{r, include=FALSE}
library(lingtypology)
```

This package is based on the [Glottolog database (v. 4.4)](https://glottolog.org/), so `lingtypology` has several functions for accessing data from that database.

### 1. Command name's syntax
Most of the functions in `lingtypology` have the same syntax: **what you need.what you have**. Most of them are based on _language name_.

* **aff.lang()** — get affiliation by language
* **area.lang()** — get macro area by language
* **country.lang()** — get country by language
* **iso.lang()** — get ISO 639-3 code by language
* **gltc.lang()** — get glottocode (identifier for a language in the Glottolog database) code by language
* **lat.lang()** — get latitude by language
* **long.lang()** — get longitude by language
* **level.lang()** — get level by language
* **subc.lang()** — get subclassification in the Newick tree format by language

Some of them help to define a vector of languages.

* **lang.aff()** — get language by affiliation
* **lang.iso()** — get language by ISO 639-3 code
* **lang.gltc()** — get language by glottocode

Additionally there are some functions to convert glottocodes to ISO 639-3 codes and vice versa:

* **gltc.iso()** — get glottocode by ISO 639-3 code
* **iso.gltc()** — get ISO 639-3 code by glottocode

The most important functionality of `lingtypology` is the ability to create interactive maps based on features and sets of languages (see the third section):

* **map.feature()**

[Glottolog database (v. 4.1)](https://glottolog.org/) provides `lingtypology` with language names, ISO codes, glottocodes, affiliation, macro area, coordinates, and much information. This set of functions doesn't have a goal to cover all possible combinations of functions. Check out additional information that is preserved in the version of the Glottolog database used in `lingtypology`:

```{r}
names(glottolog)
```

Using R functions for data manipulation you can create your own database for your purpose.

### 2. Using base functions
All functions introduced in the previous section are regular functions, so they can take the following objects as input:

* a regular string
```{r}
iso.lang("Adyghe")
lang.iso("ady")
lang.aff("West Caucasian")
```

I would like to point out that you can create strings in R using single or double quotes. Since inserting single quotes in a string created with single quotes causes an error in R, I use double quotes in my tutorial. You can use single quotes, but be careful and remember that `'Ma'ya'` is an incorrect string in R.

* a vector of strings
```{r}
area.lang(c("Adyghe", "Aduge"))
lang <- c("Adyghe", "Russian")
aff.lang(lang)
```

*  other functions. For example, let's try to get a vector of ISO codes for the Circassian languages
```{r}
iso.lang(lang.aff("Circassian"))
```

If you are new to R, it is important to mention that you can create a table with languages, features and other parametres with any spreadsheet software you used to work. Then you can import the created file to R using standard tools.

### 3. Spell Checker: look carefully at warnings!

All functions which take a vector of languages are enriched with a kind of a spell checker. If a language from a query is absent in the database, functions return a warning message containing a set of candidates with the minimal Levenshtein distance to the language from the query.

```{r}
aff.lang("Adyge")
```

### 4. `subc.lang()` function

The `subc.lang()` function returns language subclassification in the Newick tree format.

```{r}
subc.lang("Lechitic")
```

This format is hard to interpret by itself, but there are some tools in R that make it possible to visualise those subclassifications:

```{r}
library(ape)
plot(read.tree(text = subc.lang("Lechitic")))
```

It is possible to specify colors of tips in case you want to emphasize some nodes:

```{r}
plot(read.tree(text = subc.lang("Lechitic")),
     tip.color = c("red", "black", "black", "black"))
```

As you can see nodes are counted from bottom to top.

For more sophisticated tree visualization you can look into [`ggtree` package](https://bioconductor.org/packages/release/bioc/html/ggtree.html). There are several linguistic packages that provide some functionality for creating glottolog trees:

* [`glottoTrees`](https://github.com/erichround/glottoTrees) package by Erich Round
* [`lingtypr`](https://gitlab.com/laurabecker/lingtypr) package by Laura Becker
---
title: "`lingtypology`: introduction and installation"
author: "[George Moroz](mailto:agricolamz@gmail.com), [NRU HSE Linguistic Convergence Laboratory](https://ilcl.hse.ru/en/)"
---

### What is `lingtypology`?
The `lingtypology` package connects R with the [Glottolog database (v. 4.4)](http://glottolog.org/) and provides an additional functionality for linguistic typology. The Glottolog database contains a catalogue of the world's languages. This package helps researchers to make linguistic maps, using the philosophy of [the Cross-Linguistic Linked Data project](http://clld.org/), which is creating a uniform access to linguistic data across publications. This package is based on the [leaflet package](https://rstudio.github.io/leaflet/), so `lingtypology` is a package for interactive linguistic mapping. In addition, the package provides an ability to download data from typological databases such as WALS, AUTOTYP and others (see section 4). The functionality of this package intersects with the package [`lingtypr`](https://gitlab.com/laurabecker/lingtypr) by Laura Becker. I would like to thank Natalya Tyshkevich, Samira Verhees and Eugenya  Klyagina for reading and correcting some versions of this vignette.

### 1. Installation
Since `lingtypology` is an R package, you should install [R (version >= 3.1.0)](https://www.r-project.org/) on your PC if you haven't already done so. To install the `lingtypology` package, run the following command at your R IDE, so you get the stable version from CRAN:
```{r, eval=FALSE}
install.packages("lingtypology", dependencies = TRUE)
```

You can also get the development version from GitHub:
```{r, eval= F}
install.packages("devtools")
devtools::install_github("ropensci/lingtypology")
```

Load the package:
```{r}
library(lingtypology)
```

### 2. Citing `lingtyplogy`
It is important to cite R and R packages when you use them. For this purpose use the `citation` function:
```{r}
citation("lingtypology")
```
---
title: "Introduction to leaflet.minicharts"
author: "Francois Guillem"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

For a few years now, it has become very to create interactive maps with R thanks to the package `leaflet` by the Rstudio team. Nevertheless, it only provides only a few functions to create basic shapes on a map, so the information that can be represented on a single map is limited: if you have some data associated to some points, you can only represent at most two variables by drawing circles and changing their radius and color according to data.

`leaflet.minicharts` is an R package that provides two functions to add and update small charts on an interactive maps created with the package `leaflet`. These charts can be used to represent as many variables as desired associated to geographical points. Currently, three types of chart are supported: barcharts (the default), pie charts and polar area charts.

let's have a look to a concrete example.

## Data

The package provides a table that contains the electric production, consumption and exchanges of France from january 2010 to february 2017 and of 12 french regions from january 2013 to february 2017. 

In addition to the total production, the table contains one column for each type of production. The table also contains the latitude and longitude of the center of the regions.

```{r}
library(leaflet.minicharts)
data("eco2mix")
head(eco2mix)
```

## Renewable productions in 2016

Nowadays, France has an objective of 23% of renewable energies in the consumption of the country by 2020. Are the country close to its objective. Is the share of renewable energies similar in all regions? 

To answer this question let us focus on the year 2016 We first prepare the required data with package `dplyr`:

```{r message=FALSE}
library(dplyr)

prod2016 <- eco2mix %>%
  mutate(
    renewable = bioenergy + solar + wind + hydraulic,
    non_renewable = total - bioenergy - solar - wind - hydraulic
  ) %>%
  filter(grepl("2016", month) & area != "France") %>%
  select(-month) %>%
  group_by(area, lat, lng) %>%
  summarise_all(sum) %>%
  ungroup()

head(prod2016)
```

We also create a base map that will be used in all the following examples

```{r message=FALSE, results='hide'}
library(leaflet)

tilesURL <- "http://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}"

basemap <- leaflet(width = "100%", height = "400px") %>%
  addTiles(tilesURL)
```

We now add to the base map a pie chart for each region that represents the share of renewable energies. We also change the width of the pie charts so their area is proportional to the total production of the corresponding region.

```{r}
colors <- c("#4fc13c", "#cccccc")

basemap %>%
  addMinicharts(
    prod2016$lng, prod2016$lat,
    type = "pie",
    chartdata = prod2016[, c("renewable", "non_renewable")], 
    colorPalette = colors, 
    width = 60 * sqrt(prod2016$total) / sqrt(max(prod2016$total)), transitionTime = 0
  )
```

We can see that the three south east regions exceed the target of 23%, but most regions are far from this objective. Globally, renewable energies represented only 19% percent of the production of 2016. 

Now let's represent the different types of renewable production using bar charts.

```{r}
renewable2016 <- prod2016 %>% select(hydraulic, solar, wind)
colors <- c("#3093e5", "#fcba50", "#a0d9e8")
basemap %>%
  addMinicharts(
    prod2016$lng, prod2016$lat,
    chartdata = renewable2016,
    colorPalette = colors,
    width = 45, height = 45
  )
```

Hydraulic production is far more important than solar and wind. Without surprise, solar production is more important in south while wind production is more important in the north.

## Representing a single variable

`leaflet.minicharts` has been designed to represent multiple variables at once, but you still may want to use it to represent a single variable. In the next example, we represent the total load of each french region in 2016. When data passed to `addMinicharts` contains a single column, it automatically represents it with circle which area is proportional to the corresponding value. In the example we also use the parameter `showLabels` to display rounded values of the variable inside the circles. 

```{r}
basemap %>%
  addMinicharts(
    prod2016$lng, prod2016$lat,
    chartdata = prod2016$load,
    showLabels = TRUE,
    width = 45
  )
```

This is nice, isn't it?

## Animated maps

Until now, we have only represented aggregated data but it would be nice to create a map that represents the evolution over time of some variables. It is actually easy with `leaflet.minicharts`. The first step is to construct a table containing longitude, latitude, a time column and the variables we want to represent. The table `eco2mix` already has all these columns. We only need to filter the rows containing data for the entire country.

```{r}
prodRegions <- eco2mix %>% filter(area != "France")
```

Now we can create our animated map by using the argument "time":

```{r}
basemap %>% 
  addMinicharts(
    prodRegions$lng, prodRegions$lat, 
    chartdata = prodRegions[, c("hydraulic", "solar", "wind")],
    time = prodRegions$month,
    colorPalette = colors,
    width = 45, height = 45
  )
```

## Represent flows

Since version 0.2, `leaflet.minicharts` has also functions to represent flows between points and their evolution. To illustrate this, let's represent the evolution of electricity exchanges between France and Neighboring countries. 

To do that, we use function `addFlows`. It requires coordinates of two points for each flow and the value of the flow. Other arguments are similar to `addMinicharts`.

```{r}
data("eco2mixBalance")
bal <- eco2mixBalance
basemap %>%
  addFlows(
    bal$lng0, bal$lat0, bal$lng1, bal$lat1,
    flow = bal$balance,
    time = bal$month
  )
```

Of course, you can represent flows and minicharts on the same map!

## Use in shiny web applications

In shiny applications, you can create nice transition effects by using functions `leafletproxy` and `updateMinicharts`/`updateFlows`. In the server function you first need to initialize the map and the minicharts. The important thing here is to use parameter `layerId` so that `updateMinicharts` can know which chart to update with which values. 

```{r eval = FALSE}
server <- function(input, output, session) {
  # Initialize map
  output$mymap <- renderLeaflet(
    leaflet() %>% addTiles() %>%
      addMinicharts(lon, lat, layerId = uniqueChartIds)
  )
}
```

Then use `leafletProxy()` and `updateMinicharts` in your reactive code:

```{r eval = FALSE} 
server <- function(input, output, session) {
  # Initialize map
  ...
  
  # Update map
  observe({
    newdata <- getData(input$myinput)
    
    leafletProxy("mymap") %>% 
      updateMinicharts(uniqueChartIds, chartdata = newdata, ...)
  })
}
```

You can find a [live example here](https://francoisguillem.shinyapps.io/shiny-demo/).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/soundcomparisons.feature.R
\name{soundcomparisons.feature}
\alias{soundcomparisons.feature}
\title{Download SOUNDCOMPARISONS data}
\usage{
soundcomparisons.feature(word)
}
\arguments{
\item{word}{A character vector that define with a feature ids from SOUNDCOMPARISONS (e. g. "one", "sharp_fem", "near_neut", "on_the_left", "    I_will_give", "write_ipv_sg", "your_pl_pl").}
}
\description{
This function downloads data from SOUNDCOMPARISONS (\url{https://soundcomparisons.com/}) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# soundcomparisons.feature(c("sun", "house"))
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}},  \code{\link{eurasianphonology.feature}}, \code{\link{eurasianphonology.feature}}

\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}
}
\author{
Anna Smirnova
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phoible.R
\docType{data}
\name{phoible}
\alias{phoible}
\title{Phoible glottolog - language correspondencies}
\format{
A data frame with 2185 rows and 2 variables:
\describe{
  \item{language}{language}
  \item{Glottocode}{Glottocode}
}
}
\usage{
phoible
}
\description{
Language correspondencies for Phoible (\url{https://phoible.org/}). This dataset is created for \code{\link{phoible.feature}} function.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oto_mangueanIC.feature.R
\name{oto_mangueanIC.feature}
\alias{oto_mangueanIC.feature}
\title{Download Oto-Manguean Inflectional Class Database data}
\usage{
oto_mangueanIC.feature()
}
\description{
This function downloads data from Oto-Manguean Inflectional Class Database (\url{https://oto-manguean.surrey.ac.uk/}) and creates a language column with the names from lingtypology database. You need the internet connection.
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{wals.feature}}

\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}
# oto_mangueanIC.feature()
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eurasianphonology.feature.R
\name{eurasianphonology.feature}
\alias{eurasianphonology.feature}
\title{Opens data from the database of Eurasian phonological inventories}
\usage{
eurasianphonology.feature()
}
\description{
This function opens downloaded data from the database of Eurasian phonological inventories (\url{https://eurphon.info}).
}
\examples{

eurasianphonology.feature()

}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}
}
\author{
Kirill Koncha <majortomblog@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autotyp.R
\docType{data}
\name{autotyp}
\alias{autotyp}
\title{AUTOTYP's Language identifiers}
\format{
An object of class \code{data.frame} with 2853 rows and 2 columns.
}
\usage{
autotyp
}
\description{
Language identifiers from AUTOTYP v. 0.1.4 (\url{https://github.com/autotyp/autotyp-data/}). This dataset is created for \code{\link{autotyp.feature}} function.
}
\details{
#' @format A data frame with 2853 rows and 2 variables:
\describe{
  \item{LID}{language identifier}
  \item{Glottocode}{Glottocode}
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sails.feature.R
\name{sails.feature}
\alias{sails.feature}
\title{Download SAILS data}
\usage{
sails.feature(features, na.rm = TRUE)
}
\arguments{
\item{features}{A character vector that define with a feature ids from SAILS (e. g. "and1", "argex4-1-3").}

\item{na.rm}{Logical. If TRUE function removes all languages not available in lingtypology database. By default is TRUE.}
}
\description{
This function downloads data from SAILS (\url{https://sails.clld.org/}) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# sails.feature(c("and1", "and11"))
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}

\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{wals.feature}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glottolog.R
\docType{data}
\name{glottolog}
\alias{glottolog}
\title{Catalogue of languages of the world}
\format{
A data frame with 26101 rows and 10 variables:
\describe{
  \item{glottocode}{languoid code from Glottolog 4.5}
  \item{language}{name of the language}
  \item{iso}{code based on ISO 639--3 \url{https://iso639-3.sil.org/}}
  \item{level}{languoid type: dialect or language (possible values are dialect, language, family, bookkeeping, pseudo family, sign language, unclassifiable, pidgin, unattested, artificial language, speech register, mixed language)}
  \item{area}{have six values Africa, Australia, Eurasia, North America, Papunesia, South America}
  \item{latitude}{latitude}
  \item{longitude}{longitude}
  \item{countries}{list of countries, where the language is spoken}
  \item{affiliation}{genealogical affiliation}
  \item{subclassification}{subclassification in a Newick format}
}
}
\source{
\url{https://glottolog.org/}
}
\usage{
glottolog
}
\description{
A dataset containes the original catalogue of languages of the world
involving genealogical affiliation, macro-area, country, iso code,
and coordinates.
}
\details{
Hammarström, Harald & Forkel, Robert & Haspelmath, Martin & Bank, Sebastian. 2021.
Glottolog 4.5.
Leipzig: Max Planck Institute for Evolutionary Anthropology.
https://doi.org/10.5281/zenodo.5772642
(Available online at http://glottolog.org, Accessed on 2021-12-10.)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phoible.feature.R
\name{phoible.feature}
\alias{phoible.feature}
\title{Download PHOIBLE data}
\usage{
phoible.feature(source = "all", na.rm = TRUE)
}
\arguments{
\item{source}{A character vector that define with a source names from PHOIBLE (possible values: "all", "aa", "gm", "ph", "ra", "saphon", "spa", "upsid").}

\item{na.rm}{Logical. If TRUE function removes all languages not available in lingtypology database. By default is TRUE.}
}
\description{
This function downloads data from PHOIBLE (\url{https://phoible.org/}) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# phoible.feature()
# phoible.feature(c('consonants', 'vowels'), source = "UPSID")
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}

\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{sails.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{wals.feature}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/atlas.database.R
\name{atlas.database}
\alias{atlas.database}
\title{Create an atlas}
\usage{
atlas.database(
  languages,
  latitude,
  longitude,
  features,
  atlas.name = "",
  author = ""
)
}
\arguments{
\item{languages}{character vector of languages (can be written in lower case)}

\item{latitude}{numeric vector of latitudes (optional)}

\item{longitude}{numeric vector of longitudes (optional)}

\item{features}{dataframe where each column is a feature set}

\item{atlas.name}{string with an atlas name}

\item{author}{string with the authors list}
}
\description{
This function creates an rmarkdown based atlas from data provided by users. This function creates the template, after it should be rendered by rmarkdown package. The DT package is required during the rendering.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bivaltyp.feature.R
\name{bivaltyp.feature}
\alias{bivaltyp.feature}
\title{Download BivalTyp data}
\usage{
bivaltyp.feature()
}
\description{
This function downloads data from BivalTyp (\url{https://www.bivaltyp.info/}) and changes language names to the names from lingtypology database. You need the internet connection.
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{valpal.feature}}, \code{\link{wals.feature}}

\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}
# bivaltyp.feature()
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iso.gltc.R
\name{iso.gltc}
\alias{iso.gltc}
\title{Get ISO 639--3 code by Glottocode}
\usage{
iso.gltc(x)
}
\arguments{
\item{x}{A character vector of Glottocodes.}
}
\description{
Takes any vector of Glotocodes and returns ISO code.
}
\examples{
iso.gltc('adyg1241')
iso.gltc(c('adyg1241', 'udii1243'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{area.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subc.lang.R
\name{subc.lang}
\alias{subc.lang}
\title{Get subclassification by language}
\usage{
subc.lang(x)
}
\arguments{
\item{x}{A character vector of the languoids (can be written in lower case)}
}
\description{
Takes any vector of languoids and returns subclassification in the Newick tree format.
}
\examples{
subc.lang('Korean')
subc.lang(c('Korean', 'Lechitic'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{area.lang}}, \code{\link{country.lang}}, \code{\link{gltc.lang}}, \code{\link{iso.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.glottolog.R
\name{is.glottolog}
\alias{is.glottolog}
\title{Are these languages in glottolog?}
\usage{
is.glottolog(x, response = FALSE)
}
\arguments{
\item{x}{A character vector of languages (can be written in lower case)or ISO codes}

\item{response}{logical. If TRUE, when language is absent, return warnings with a possible candidates.}
}
\description{
Takes any vector of languages or ISO codes and returns a logical vector.
}
\examples{
is.glottolog(c('Adyghe', 'Russian'))
is.glottolog('Buyaka')

# Add warning message with sugestions
is.glottolog(c('Adygey', 'Russian'), response = TRUE)
# > FALSE TRUE
# Warning message:
# In is.glottolog(c('Adyge', 'Russian'), response = TRUE) :
# Language Adyge is absent in our version of the Glottolog database. Did you mean Aduge, Adyghe?

}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abvd.R
\docType{data}
\name{abvd}
\alias{abvd}
\title{ABVD's Language identifiers}
\format{
A data frame with 1468 rows and 2 variables:
\describe{
  \item{id}{language identifier}
  \item{glottocode}{Glottocode}
}
}
\usage{
abvd
}
\description{
Language identifiers from ABVD (\url{https://abvd.shh.mpg.de/austronesian/}). This dataset is created for \code{\link{abvd.feature}} function.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lang.aff.R
\name{lang.aff}
\alias{lang.aff}
\title{Get languages by affiliation}
\usage{
lang.aff(x, include.dialects = FALSE, list = FALSE)
}
\arguments{
\item{x}{A character vector of the affiliations (can be written in lower case)}

\item{include.dialects}{logical. If TRUE, it returns all langauges and dialects, if FALSE it returns only languages.}

\item{list}{logical. If TRUE, it returns a list of languages, if FALSE it returns a named vector.}
}
\description{
Takes any vector of affiliations and returns languages.
}
\examples{
lang.aff('Slavic')
lang.aff(c('Slavic', 'Celtic'))
lang.aff(c('Slavic', 'Celtic'), list = TRUE)
}
\seealso{
\code{\link{lang.iso}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygon.points_kde.R
\name{polygon.points_kde}
\alias{polygon.points_kde}
\title{Get kernel density estimation poligon from coordinates}
\usage{
polygon.points_kde(latitude, longitude, latitude.width, longitude.width)
}
\arguments{
\item{latitude}{numeric vector of latitudes}

\item{longitude}{numeric vector of longitudes}

\item{latitude.width}{bandwidths for latitude values. Defaults to normal reference bandwidth (see \link{bandwidth.nrd}).}

\item{longitude.width}{bandwidths for longitude values. Defaults to normal reference bandwidth (see \link{bandwidth.nrd}).}
}
\description{
This function is based on this answer: https://gis.stackexchange.com/a/203623/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggmap.feature.R
\name{ggmap.feature}
\alias{ggmap.feature}
\title{Create a map with ggplot2}
\usage{
ggmap.feature(
  languages,
  features = "",
  latitude = NA,
  longitude = NA,
  color = NULL,
  title = NULL,
  legend = TRUE,
  width = 2,
  opacity = 1,
  map.orientation = "Atlantic"
)
}
\arguments{
\item{languages}{character vector of languages (can be written in lower case).}

\item{features}{character vector of features.}

\item{latitude}{numeric vector of latitudes.}

\item{longitude}{numeric vector of longitudes.}

\item{color}{vector of colors or palette. The color argument can be (1) a character vector of RGM or named colors; (2) the name of an RColorBrewer palette; (3) the full name of a viridis palette; (4) a function that receives a single value between 0 and 1 and returns a color. For more examples see \code{\link{colorNumeric}}.}

\item{title}{title of a legend.}

\item{legend}{logical. If TRUE, function show legend. By default is TRUE.}

\item{width}{a numeric vector of radius for circles or width for barcharts in minicharts.}

\item{opacity}{a numeric vector of marker opacity.}

\item{map.orientation}{a character verctor with values "Pacific" and "Atlantic". It distinguishes Pacific-centered and Atlantic-centered maps. By default is "Atlantic".}
}
\description{
Map a set of languages and color them by feature.
}
\examples{
ggmap.feature(c("Adyghe", "Russian"))

}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lang.iso.R
\name{lang.iso}
\alias{lang.iso}
\title{Get language by ISO 639--3 code}
\usage{
lang.iso(x)
}
\arguments{
\item{x}{A character vector of the ISO codes.}
}
\description{
Takes any vector of ISO codes and returns languages.
}
\examples{
lang.iso('ady')
lang.iso(c('ady', 'rus'))
}
\seealso{
\code{\link{lang.aff}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circassian.R
\docType{data}
\name{circassian}
\alias{circassian}
\title{Circassian villages in Russia}
\format{
A data frame with 158 rows and 6 variables:
\describe{
  \item{longitude}{longitude}
  \item{latitude}{latitude}
  \item{village}{name of the village}
  \item{district}{names of the subjects of the Russian Federation: kbr --- Kabardino--Balkar Republic, kch --- Karachay--Cherkess Republic, kk --- Krasnodar Krai, ra --- Republic of Adygea, stv --- Stavropol Krai}
  \item{dialect}{names of the Circassian dialects}
  \item{language}{according standard Circassian devision there are Adyghe and Kabardian languages}
}
}
\usage{
circassian
}
\description{
A dataset containes the list of the Circassian villages in Russia
with genealogical affiliation, coordinates and district names. Most
data collected during the fieldworks (2011--2018).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/url.lang.R
\name{url.lang}
\alias{url.lang}
\title{Make a url-link to glottolog page for a language}
\usage{
url.lang(x, popup = "")
}
\arguments{
\item{x}{A character vector of languages (can be written in lower case)}

\item{popup}{character vector of strings that will appear in pop-up window of the function map.feature}
}
\description{
Takes any vector of languages and returns links to glottolog pages.
}
\examples{
url.lang('Korean')
url.lang(c('Gangou', 'Hachijo', 'Adyghe', 'Ganai'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{area.lang}}, \code{\link{country.lang}}, \code{\link{gltc.lang}}, \code{\link{iso.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}, \code{\link{subc.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wals.R
\docType{data}
\name{wals}
\alias{wals}
\title{WALS's Language identifiers}
\format{
A data frame with 2950 rows and 2 variables:
\describe{
  \item{wals.code}{WALS language identifier}
  \item{glottocode}{Glottocode}
}
}
\usage{
wals
}
\description{
Language identifiers from WALS (\url{https://wals.info/}). This dataset is created for \code{\link{wals.feature}} function.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oto_mangueanIC.R
\docType{data}
\name{oto_mangueanIC}
\alias{oto_mangueanIC}
\title{Oto-Manguean Inflectional Class Database Language identifiers}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 20 rows and 2 columns.
}
\usage{
oto_mangueanIC
}
\description{
Language identifiers from Oto-Manguean Inflectional Class Database (\url{https://oto-manguean.surrey.ac.uk/}). This dataset is created for \code{\link{oto_mangueanIC.feature}} function.
}
\details{
#' @format A data frame with 20 rows and 2 variables:
\describe{
  \item{Language.name}{Languaeg names from Oto-Manguean Inflectional Class Database}
  \item{language}{Language names from Glottolog database}
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/soundcomparisons.R
\docType{data}
\name{soundcomparisons}
\alias{soundcomparisons}
\title{SOUNDCOMPARISONS's Language identifiers}
\format{
An object of class \code{data.frame} with 556 rows and 3 columns.
}
\usage{
soundcomparisons
}
\description{
Language identifiers from SOUNDCOMPARISONS (\url{https://soundcomparisons.com/}). This dataset is created for \code{\link{soundcomparisons.feature}} function.
}
\details{
#' @format A data frame with 556 rows and 2 variables:
\describe{
  \item{LanguageName}{SOUNDCOMPARISONS language identifier}
  \item{LanduageId}{Language Id}
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abvd.feature.R
\name{abvd.feature}
\alias{abvd.feature}
\title{Download ABVD data}
\usage{
abvd.feature(feature)
}
\arguments{
\item{feature}{A character vector that define a language id from ABVD (e. g. "1", "292").}
}
\description{
This function downloads data from ABVD (\url{https://abvd.shh.mpg.de/austronesian/}) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# abvd.feature(c(292, 7))
}
\seealso{
\code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gltc.iso.R
\name{gltc.iso}
\alias{gltc.iso}
\title{Get Glottocode by ISO 639--3 code}
\usage{
gltc.iso(x)
}
\arguments{
\item{x}{A character vector of the Glottocodes.}
}
\description{
Takes any vector of ISO 639--3 codes and returns Glottocodes.
}
\examples{
gltc.iso('ady')
gltc.iso(c('ady', 'rus'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{area.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gltc.lang.R
\name{gltc.lang}
\alias{gltc.lang}
\title{Get Glottocode by language}
\usage{
gltc.lang(x)
}
\arguments{
\item{x}{A character vector of the languages (can be written in lower case)}
}
\description{
Takes any vector of languages and returns Glottocode.
}
\examples{
gltc.lang('Adyghe')
gltc.lang(c('Adyghe', 'Udi'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{area.lang}}, \code{\link{country.lang}}, \code{\link{iso.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}, \code{\link{subc.lang}}, \code{\link{url.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/providers.R
\docType{data}
\name{providers}
\alias{providers}
\title{Providers}
\format{
A list of characters
}
\source{
\url{https://github.com/leaflet-extras/leaflet-providers/blob/master/leaflet-providers.js}
}
\usage{
providers
}
\description{
List of all providers with their variations taken from leaflet package
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uralex.feature.R
\name{uralex.feature}
\alias{uralex.feature}
\title{Download UraLex data}
\usage{
uralex.feature(na.rm = TRUE)
}
\arguments{
\item{na.rm}{Logical. If TRUE function removes all languages not available in lingtypology database. By default is TRUE.}
}
\description{
This function downloads data from UraLex (https://github.com/lexibank/uralex/) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# uralex.feature()
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lang.country.R
\name{lang.country}
\alias{lang.country}
\title{Get language by country}
\usage{
lang.country(x, list = FALSE)
}
\arguments{
\item{x}{character vector of the countries (in alpha-2 ISO codes)}

\item{list}{logical. If TRUE, it returns a list of languages, if FALSE it returns a named vector.}
}
\description{
Takes any vector of countries and returns languages.
}
\examples{
lang.country('AD')
lang.country(c('AD', 'AE'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{country.lang}}, \code{\link{gltc.lang}}, \code{\link{iso.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}, \code{\link{subc.lang}}, \code{\link{url.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/long.lang.R
\name{long.lang}
\alias{long.lang}
\title{Get longitude by language}
\usage{
long.lang(x, map.orientation = "Pacific")
}
\arguments{
\item{x}{A character vector of the languages (can be written in lower case)}

\item{map.orientation}{A character verctor with values "Pacific" and "Atlantic". It distinguishes Pacific-centered and Atlantic-centered maps. By default is "Pacific".}
}
\description{
Takes any vector of languages and returns longitude.
}
\examples{
lat.lang('Adyghe')
long.lang('Adyghe')
lat.lang(c('Adyghe', 'Russian'))
long.lang(c('Adyghe', 'Russian'))
long.lang(c('Adyghe', 'Aleut'), map.orientation = "Pacific")
}
\seealso{
\code{\link{aff.lang}}, \code{\link{area.lang}}, \code{\link{country.lang}}, \code{\link{gltc.lang}}, \code{\link{iso.lang}}, \code{\link{lat.lang}},  \code{\link{subc.lang}}, \code{\link{url.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/level.lang.R
\name{level.lang}
\alias{level.lang}
\title{Get a level of language by language}
\usage{
level.lang(x)
}
\arguments{
\item{x}{character vector of the languages (can be written in lower case)}
}
\description{
Takes any vector of languages and returns a level of language.
}
\examples{
level.lang('Russian Sign Language')
level.lang(c('Archi', 'Chechen'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{country.lang}}, \code{\link{gltc.lang}}, \code{\link{iso.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}, \code{\link{subc.lang}}, \code{\link{url.lang}}
}
\author{
Sasha Shakhnova
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phonological_profiles.R
\docType{data}
\name{phonological_profiles}
\alias{phonological_profiles}
\title{Number of consonants and presence of ejectives}
\format{
A data frame with 19 rows and 4 variables:
\describe{
  \item{language}{language name}
  \item{consonants}{number of consonants. Based on UPSID database.}
  \item{vowels}{number of vowels. Based on UPSID database.}
  \item{ejectives}{presence of ejective sounds.}
  \item{tone}{presence of tone.}
  \item{stress}{presence of stress.}
  \item{long_vowels}{presence of long vowels.}
}
}
\usage{
phonological_profiles
}
\description{
Number of consonants and presence of ejectives
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autotyp.feature.R
\name{autotyp.feature}
\alias{autotyp.feature}
\title{Download AUTOTYP data}
\usage{
autotyp.feature(features, na.rm = TRUE)
}
\arguments{
\item{features}{A character vector that define with a feature names from AUTOTYP.}

\item{na.rm}{Logical. If TRUE function removes all languages not available in lingtypology database. By default is TRUE.}
}
\description{
This function downloads data from AUTOTYP (https://github.com/autotyp/autotyp-data#the-autotyp-database) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# autotyp.feature(c('Gender', 'Numeral classifiers'))
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}

\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{wals.feature}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/afbo.feature.R
\name{afbo.feature}
\alias{afbo.feature}
\title{Download AfBo data}
\usage{
afbo.feature(features = "all", na.rm = TRUE)
}
\arguments{
\item{features}{A character vector that define with an affix functions from AfBo (e. g. "all", "adjectivizer", "focus").}

\item{na.rm}{Logical. If TRUE function removes all languages not available in lingtypology database. By default is TRUE.}
}
\description{
This function downloads data from AfBo (\url{https://afbo.info/}) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# afbo.feature()
# afbo.feature(c("adjectivizer", "adverbializer"))
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}

\code{\link{abvd.feature}}, \code{\link{autotyp.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{wals.feature}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polygon.points_fd.R
\name{polygon.points_fd}
\alias{polygon.points_fd}
\title{Get poligons from fixed distance circles around coordinates}
\usage{
polygon.points_fd(latitude, longitude, width)
}
\arguments{
\item{latitude}{numeric vector of latitudes}

\item{longitude}{numeric vector of longitudes}

\item{width}{radius for creating poligons around points}
}
\description{
This function is based on this answer: https://www.r-bloggers.com/merging-spatial-buffers-in-r/
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eurasianphonology.R
\docType{data}
\name{eurasianphonology}
\alias{eurasianphonology}
\title{Eurasianphonology data}
\format{
A data frame with 19825 rows and 19 variables:
\describe{
  \item{id}{Language id}
  \item{iso}{ISO code}
  \item{name}{Another language name}
  \item{type}{Language or dialect}
  \item{language}{Language name}
  \item{latitude}{latitude}
  \item{longitude}{longitude}
  \item{gen1}{Language Family}
  \item{gen2}{Language Family}
  \item{tones}{Inventory of tones}
  \item{syllab}{Syllab structure}
  \item{cluster}{Cluster}
  \item{finals}{Finals}
  \item{source}{Source}
  \item{comment}{Comment}
  \item{contr}{Contributor}
  \item{segment_type}{Vowels or consonants}
  \item{segments}{Segments}
  \item{glottocode}{Glottocode}
}
}
\usage{
eurasianphonology
}
\description{
Data from The database of Eurasian phonological inventories (\url{https://eurphon.info}). This dataset is created for \code{\link{eurasianphonology.feature}} function.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vanuatu.feature.R
\name{vanuatu.feature}
\alias{vanuatu.feature}
\title{Download Vanuatu Voices data}
\usage{
vanuatu.feature(features, na.rm = TRUE)
}
\arguments{
\item{features}{A vector with parameters from Concepts (\url{https://vanuatuvoices.clld.org/parameters}))}

\item{na.rm}{Logical. If TRUE function removes all languages not available in lingtypology database. By default is TRUE.}
}
\description{
This function downloads data from Vanuatu Voices (\url{https://vanuatuvoices.clld.org/}). You need the internet connection.
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{wals.feature}}
}
\author{
Mikhail Leonov
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iso.lang.R
\name{iso.lang}
\alias{iso.lang}
\title{Get ISO 639--3 code by language}
\usage{
iso.lang(x)
}
\arguments{
\item{x}{A character vector of the languages (can be written in lower case)}
}
\description{
Takes any vector of languages and returns ISO code.
}
\examples{
iso.lang('Adyghe')
iso.lang(c('Adyghe', 'Udi'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{area.lang}}, \code{\link{country.lang}}, \code{\link{gltc.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}, \code{\link{subc.lang}}, \code{\link{url.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wals.feature.R
\name{wals.feature}
\alias{wals.feature}
\title{Download WALS data}
\usage{
wals.feature(features, na.rm = TRUE)
}
\arguments{
\item{features}{A character vector that define with a feature ids from WALS (e. g. "1a", "21b").}

\item{na.rm}{Logical. If TRUE function removes all languages not available in lingtypology database. By default is TRUE.}
}
\description{
This function downloads data from WALS (\url{https://wals.info/}) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# wals.feature(c("1a", "20a"))
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}, \code{\link{vanuatu.feature}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amap.R
\docType{data}
\name{amap}
\alias{amap}
\title{Atlantic center template for ggmap.feature() function}
\format{
A list with 9 variables.
}
\usage{
amap
}
\description{
.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bantu.R
\docType{data}
\name{bantu}
\alias{bantu}
\title{BANTU's Language identifiers}
\format{
A data frame with 430 rows and 2 variables:
\describe{
  \item{id}{BANTU word id}
  \item{word}{word}
}
}
\usage{
bantu
}
\description{
Language identifiers from BANTU (\url{https://abvd.shh.mpg.de/bantu/index.php}). This dataset is created for \code{\link{bantu.feature}} function.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/valpal.feature.R
\name{valpal.feature}
\alias{valpal.feature}
\title{Download ValPaL data}
\usage{
valpal.feature(na.rm = FALSE)
}
\arguments{
\item{na.rm}{Logical. If TRUE function removes all languages not available in lingtypology database. By default is FALSE.}
}
\description{
This function downloads data from ValPal (www.valpal.info/) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# valpal.feature()
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{bivaltyp.feature}}, \code{\link{eurasianphonology.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{soundcomparisons.feature}}, \code{\link{uralex.feature}}, \code{\link{vanuatu.feature}}, \code{\link{wals.feature}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uralex.R
\docType{data}
\name{uralex}
\alias{uralex}
\title{UraLex's Language identifiers}
\format{
A data frame with 27 rows and 3 variables:
\describe{
  \item{language}{language name from database}
  \item{Glottocode}{Glottocodes}
  \item{language2}{language from lingtypology}
}
}
\usage{
uralex
}
\description{
Language identifiers from UraLex (\url{https://github.com/lexibank/uralex/}). This dataset is created for \code{\link{uralex.feature}} function.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/country.lang.R
\name{country.lang}
\alias{country.lang}
\title{Get country by language}
\usage{
country.lang(x)
}
\arguments{
\item{x}{A character vector of the languages (can be written in lower case)}
}
\description{
Takes any vector of languages and returns countries where those languages are used as ISO 3166-1 alpha-2 codes.
}
\examples{
country.lang('Korean')
country.lang(c('Korean', 'Polish'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{area.lang}}, \code{\link{gltc.lang}}, \code{\link{iso.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}, \code{\link{subc.lang}}, \code{\link{url.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lang.gltc.R
\name{lang.gltc}
\alias{lang.gltc}
\title{Get language by Glottocode}
\usage{
lang.gltc(x)
}
\arguments{
\item{x}{A character vector of the Glottocodes.}
}
\description{
Takes any vector of Glottocodes and returns languages.
}
\examples{
lang.gltc('adyg1241')
lang.gltc(c('adyg1241', 'udii1243'))
}
\seealso{
\code{\link{lang.aff}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bantu.feature.R
\name{bantu.feature}
\alias{bantu.feature}
\title{Download BANTU data}
\usage{
bantu.feature(features)
}
\arguments{
\item{features}{A character vector that define with a feature ids from BANTU ('house', 'cat').}
}
\description{
This function downloads data from Bantu Basic Vocabulary Database (\url{https://abvd.shh.mpg.de/bantu/index.php}) and changes language names to the names from lingtypology database. You need the internet connection.
}
\examples{
# bantu.feature(c('house', 'cat'))
}
\seealso{
\code{\link{abvd.feature}}, \code{\link{afbo.feature}}, \code{\link{autotyp.feature}}, \code{\link{oto_mangueanIC.feature}}, \code{\link{phoible.feature}}, \code{\link{sails.feature}}, \code{\link{uralex.feature}}, \code{\link{valpal.feature}}
}
\author{
Anna Smirnova <annedadaa@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aff.lang.R
\name{aff.lang}
\alias{aff.lang}
\title{Get affiliation by language}
\usage{
aff.lang(x)
}
\arguments{
\item{x}{A character vector of the languages (can be written in lower case)}
}
\description{
Takes any vector of languages and returns affiliation.
}
\examples{
aff.lang('Korean')
aff.lang(c('Korean', 'Polish'))
}
\seealso{
\code{\link{area.lang}}, \code{\link{country.lang}}, \code{\link{gltc.lang}}, \code{\link{iso.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}, \code{\link{subc.lang}}, \code{\link{url.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lat.lang.R
\name{lat.lang}
\alias{lat.lang}
\title{Get latitude by language}
\usage{
lat.lang(x)
}
\arguments{
\item{x}{A character vector of the languages (can be written in lower case)}
}
\description{
Takes any vector of languages and returns latitude.
}
\examples{
lat.lang('Adyghe')
long.lang('Adyghe')
lat.lang(c('Adyghe', 'Russian'))
long.lang(c('Adyghe', 'Russian'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{area.lang}}, \code{\link{country.lang}}, \code{\link{gltc.lang}}, \code{\link{iso.lang}}, \code{\link{long.lang}}, \code{\link{subc.lang}}, \code{\link{url.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/area.lang.R
\name{area.lang}
\alias{area.lang}
\title{Get macro area by language}
\usage{
area.lang(x)
}
\arguments{
\item{x}{character vector of the languages (can be written in lower case)}
}
\description{
Takes any vector of languages and returns macro area.
}
\examples{
area.lang('Adyghe')
area.lang(c('Adyghe', 'Aduge'))
}
\seealso{
\code{\link{aff.lang}}, \code{\link{country.lang}}, \code{\link{gltc.lang}}, \code{\link{iso.lang}}, \code{\link{lat.lang}}, \code{\link{long.lang}}, \code{\link{subc.lang}}, \code{\link{url.lang}}
}
\author{
George Moroz <agricolamz@gmail.com>
}
\name{imports}
\alias{\%>\%}
\docType{import}
\title{Objects imported from other packages}
\description{
  These objects are imported from other packages. Follow the links to their documentation.
  \describe{
    \item{magrittr}{\code{\link[magrittr:\%>\%]{\%>\%}}}
  }}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map.feature.R
\name{map.feature}
\alias{map.feature}
\title{Create a map}
\usage{
map.feature(
  languages,
  features = "",
  label = "",
  popup = "",
  latitude = NA,
  longitude = NA,
  label.hide = TRUE,
  label.fsize = 15,
  label.font = "sans-serif",
  label.position = "right",
  label.emphasize = list(NULL, "black"),
  shape = NULL,
  shape.size = 20,
  pipe.data = NULL,
  shape.color = "black",
  stroke.features = NULL,
  point.cluster = FALSE,
  density.estimation = NULL,
  density.method = "fixed distance",
  density.estimation.color = NULL,
  density.estimation.opacity = 0.6,
  density.points = TRUE,
  density.width = NULL,
  density.legend = TRUE,
  density.legend.opacity = 1,
  density.legend.position = "bottomleft",
  density.title = "",
  density.control = FALSE,
  isogloss = NULL,
  isogloss.color = "black",
  isogloss.opacity = 0.2,
  isogloss.line.width = 3,
  isogloss.width = NULL,
  color = NULL,
  stroke.color = NULL,
  image.url = NULL,
  image.width = 100,
  image.height = 100,
  image.X.shift = 0,
  image.Y.shift = 0,
  title = NULL,
  stroke.title = NULL,
  control = "",
  legend = TRUE,
  legend.opacity = 1,
  legend.position = "topright",
  stroke.legend = TRUE,
  stroke.legend.opacity = 1,
  stroke.legend.position = "bottomleft",
  width = 5,
  stroke.radius = 9.5,
  opacity = 1,
  stroke.opacity = 1,
  scale.bar = TRUE,
  scale.bar.position = "bottomleft",
  minimap = FALSE,
  minimap.position = "bottomright",
  minimap.width = 150,
  minimap.height = 150,
  facet = NULL,
  tile = "OpenStreetMap.Mapnik",
  tile.name = NULL,
  tile.opacity = 1,
  zoom.control = FALSE,
  zoom.level = NULL,
  rectangle.lng = NULL,
  rectangle.lat = NULL,
  rectangle.color = "black",
  line.lng = NULL,
  line.lat = NULL,
  line.type = "standard",
  line.color = "black",
  line.opacity = 0.8,
  line.label = NULL,
  line.width = 3,
  graticule = NULL,
  minichart = "bar",
  minichart.data = NULL,
  minichart.time = NULL,
  minichart.labels = FALSE,
  map.orientation = "Pacific",
  radius = NULL
)
}
\arguments{
\item{languages}{character vector of languages (can be written in lower case)}

\item{features}{character vector of features}

\item{label}{character vector of strings that will appear near points}

\item{popup}{character vector of strings that will appear in pop-up window}

\item{latitude}{numeric vector of latitudes}

\item{longitude}{numeric vector of longitudes}

\item{label.hide}{logical. If FALSE, labels are displayed allways. If TRUE, labels are displayed on mouse over. By default is TRUE.}

\item{label.fsize}{numeric value of the label font size. By default is 14.}

\item{label.font}{string with values of generic family: "serif", "sans-serif", "monospace", or font name e. g. "Times New Roman"}

\item{label.position}{the position of labels: "left", "right", "top", "bottom"}

\item{label.emphasize}{is the list. First argument is a vector of points in datframe that should be emphasized. Second argument is a string with a color for emphasis.}

\item{shape}{\enumerate{ \item if TRUE, creates icons (up to five categories) for values in the \code{features} variable; \item it also could be a vector of any strings that represents the levels of the  \code{features} variable; \item it also could be a string vector that represents the number of observations in dataset.}}

\item{shape.size}{size of the \code{shape} icons}

\item{pipe.data}{this variable is important, when you use map.feature with dplyr pipes. Expected usage: pipe.data = .}

\item{shape.color}{color of the \code{shape} icons}

\item{stroke.features}{additional independent stroke features}

\item{point.cluster}{logical. If TRUE, points will be united into clusters.}

\item{density.estimation}{additional independent features, used for density estimation}

\item{density.method}{string with one of the two methods: "kernal density estimation" or "fixed distance" (default)}

\item{density.estimation.color}{vector of density polygons' colors}

\item{density.estimation.opacity}{a numeric vector of density polygons opacity.}

\item{density.points}{logical. If FALSE, it doesn't show points in polygones.}

\item{density.width}{for density.method = "fixed distance" it is a numeric measure (1 is 1km). For density.method = "kernal density estimation" it is a vector with two meausures (first is latitude, secong is longitude). Defaults are normal reference bandwidth (see \link{bandwidth.nrd}).}

\item{density.legend}{logical. If TRUE, function show legend for density features. By default is FALSE.}

\item{density.legend.opacity}{a numeric vector of density-legend opacity.}

\item{density.legend.position}{the position of the legend: "topright", "bottomright", "bottomleft","topleft"}

\item{density.title}{title of a density-feature legend}

\item{density.control}{logical. If TRUE, function show layer control buttons for density plot. By default is FALSE}

\item{isogloss}{dataframe with corresponding features}

\item{isogloss.color}{vector of isoglosses' colors}

\item{isogloss.opacity}{a numeric vector of density polygons opacity.}

\item{isogloss.line.width}{a numeric value for line width}

\item{isogloss.width}{for density.method = "fixed distance" it is a numeric measure (1 is 1km). For density.method = "kernal density estimation" it is a vector with two meausures (first is latitude, secong is longitude). Defaults are normal reference bandwidth (see \link{bandwidth.nrd}).}

\item{color}{vector of colors or palette. The color argument can be (1) a character vector of RGM or named colors; (2) the name of an RColorBrewer palette; (3) the full name of a viridis palette; (4) a function that receives a single value between 0 and 1 and returns a color. For more examples see \code{\link{colorNumeric}}}

\item{stroke.color}{vector of stroke colors}

\item{image.url}{character vector of URLs with an images}

\item{image.width}{numeric vector of image widths}

\item{image.height}{numeric vector of image heights}

\item{image.X.shift}{numeric vector of image's X axis shift relative to the latitude-longitude point}

\item{image.Y.shift}{numeric vector of image's Y axis shift relative to the latitude-longitude point}

\item{title}{title of a legend.}

\item{stroke.title}{title of a stroke-feature legend.}

\item{control}{vector of grouping values that make it possible to create control panel that can turn off/on some points on the map.}

\item{legend}{logical. If TRUE, function show legend. By default is TRUE.}

\item{legend.opacity}{a numeric vector of legend opacity.}

\item{legend.position}{the position of the legend: "topright", "bottomright", "bottomleft","topleft"}

\item{stroke.legend}{logical. If TRUE, function show stroke.legend. By default is FALSE.}

\item{stroke.legend.opacity}{a numeric vector of stroke.legend opacity.}

\item{stroke.legend.position}{the position of the stroke.legend: "topright", "bottomright", "bottomleft","topleft"}

\item{width}{a numeric vector of radius for circles or width for barcharts in minicharts.}

\item{stroke.radius}{a numeric vector of stroke radii for the circles.}

\item{opacity}{a numeric vector of marker opacity.}

\item{stroke.opacity}{a numeric vector of stroke opacity.}

\item{scale.bar}{logical. If TRUE, function shows scale-bar. By default is TRUE.}

\item{scale.bar.position}{the position of the scale-bar: "topright", "bottomright", "bottomleft","topleft"}

\item{minimap}{logical. If TRUE, function shows mini map. By default is FALSE.}

\item{minimap.position}{the position of the minimap: "topright", "bottomright", "bottomleft","topleft"}

\item{minimap.width}{The width of the minimap in pixels.}

\item{minimap.height}{The height of the minimap in pixels.}

\item{facet}{character vector that provide a grouping variable. If it is no \code{NULL}, then as a result a list of leaflets for \code{sync} or \code{latticeView} functions from \code{mapview} package is returned.}

\item{tile}{a character verctor with a map tiles, popularized by Google Maps. See \href{https://leaflet-extras.github.io/leaflet-providers/preview/index.html}{here} for the complete set.}

\item{tile.name}{a character verctor with a user's map tiles' names.}

\item{tile.opacity}{numeric value from 0 to 1 denoting opacity of the tile.}

\item{zoom.control}{logical. If TRUE, function shows zoom controls. By default is FALSE.}

\item{zoom.level}{a numeric value of the zoom level.}

\item{rectangle.lng}{vector of two longitude values for rectangle.}

\item{rectangle.lat}{vector of two latitude values for rectangle.}

\item{rectangle.color}{vector of rectangle border color.}

\item{line.lng}{vector of two (or more) longitude values for line.}

\item{line.lat}{vector of two (or more) latitude values for line.}

\item{line.type}{a character string indicating which type of line is to be computed. One of "standard" (default), or "logit". The first one should be combined with the arguments line.lat and line.lng and provide simple lines. Other variant "logit" is the decision boundary of the logistic regression made using longitude and latitude coordinates (works only if feature argument have two levels).}

\item{line.color}{vector of line color.}

\item{line.opacity}{a numeric vector of line opacity.}

\item{line.label}{character vector that will appear near the line.}

\item{line.width}{a numeric vector of line width.}

\item{graticule}{a numeric vector for graticule spacing in map units between horizontal and vertical lines.}

\item{minichart}{citation from leaflet.minicharts package: "Possible values are "bar" for bar charts, "pie" for pie charts, "polar-area" and "polar-radius"."}

\item{minichart.data}{citation from leaflet.minicharts package: "A numeric matrix with number of rows equal to the number of elements in lng or lat and number of column equal to the number of variables to represent. If parameter time is set, the number of rows must be equal to the length of lng times the number of unique time steps in the data."}

\item{minichart.time}{citation from leaflet.minicharts package: "A vector with length equal to the number of rows in chartdata and containing either numbers representing time indices or dates or datetimes. Each unique value must appear as many times as the others. This parameter can be used when one wants to represent the evolution of some variables on a map."}

\item{minichart.labels}{citation from leaflet.minicharts package: "Should values be displayed above chart elements."}

\item{map.orientation}{a character verctor with values "Pacific" and "Atlantic". It distinguishes Pacific-centered and Atlantic-centered maps. By default is "Pacific".}

\item{radius}{deprecated argument}
}
\description{
Map a set of languages and color them by feature or two sets of features.
}
\examples{
map.feature(c("Adyghe", "Russian"))

}
\author{
George Moroz <agricolamz@gmail.com>
}
