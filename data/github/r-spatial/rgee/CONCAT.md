<h1 align="center">
  <br>
  <a href="https://github.com/r-spatial/rgee/"><img src="https://user-images.githubusercontent.com/16768318/118376965-5f7dca80-b5cb-11eb-9a82-47876680a3e6.png" alt="Markdownify" width="200"></a>
  <a href="https://github.com/r-earthengine/rgeeExtra/"><img src="https://user-images.githubusercontent.com/16768318/118376968-63a9e800-b5cb-11eb-83e7-3f36299e17cb.png" alt="Markdownify" width="200"></a>
  <a href="https://github.com/r-earthengine/rgeebook/"><img src="https://user-images.githubusercontent.com/16768318/118376966-60aef780-b5cb-11eb-8df2-ca70dcfe04c5.png" alt="Markdownify" width="200"></a>
  <br>
  rgee: Google Earth Engine for R
  <br>
</h1>


<h4 align="center">rgee is an R binding package for calling <a href="https://developers.google.com/earth-engine/" target="_blank">Google Earth Engine API</a> from within R. <a href="https://r-spatial.github.io/rgee/reference/index.html" target="_blank">Various functions</a> are implemented to simplify the connection with the R spatial ecosystem.</h4>

<p align="center">
<a href="https://colab.research.google.com/github/r-spatial/rgee/blob/examples/rgee_colab.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open in Colab" title="Open and Execute in Google Colaboratory"></a>
<a href="https://github.com/r-spatial/rgee/actions"><img src="https://github.com/r-spatial/rgee/workflows/R-CMD-check/badge.svg" alt="R build status"></a>
<a href="https://www.repostatus.org/#active"><img src="https://www.repostatus.org/badges/latest/active.svg" alt="Project Status: Active – The project has reached a stable, usable
state and is being actively
developed."></a>
<a href="https://codecov.io/gh/r-spatial/rgee"><img src="https://codecov.io/gh/r-spatial/rgee/branch/master/graph/badge.svg" alt="codecov"></a>
<a href="https://opensource.org/licenses/Apache-2.0"><img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg" alt="License"></a>
<a href="https://www.tidyverse.org/lifecycle/#maturing"><img src="https://img.shields.io/badge/lifecycle-maturing-blue.svg" alt="lifecycle"></a>
<a href="https://joss.theoj.org/papers/aea42ddddd79df480a858bc1e51857fc"><img src="https://joss.theoj.org/papers/aea42ddddd79df480a858bc1e51857fc/status.svg" alt="status"></a>
<a href="https://cran.r-project.org/package=rgee"><img src="https://www.r-pkg.org/badges/version/rgee" alt="CRAN
status"></a>
<a href="https://doi.org/10.5281/zenodo.3945409"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3945409.svg" alt="DOI"></a>
<br>
<a href="https://www.buymeacoffee.com/csay" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>
</p>

<p align="center">
  •
  <a href="#installation">Installation</a> &nbsp;•
  <a href="#hello-world">Hello World</a> &nbsp;•
  <a href="#how-does-rgee-work">How does rgee work?</a> &nbsp;•
  <a href="#quick-start-users-guide-for-rgee">Guides</a> &nbsp;•
  <a href="#contributing-guide">Contributing</a> &nbsp;•
  <a href="#share-the-love">Citation</a> &nbsp;•
  <a href="#credits">Credits</a>
</p>

## What is Google Earth Engine?

[Google Earth Engine](https://earthengine.google.com/) is a cloud-based platform that lets users access a petabyte-scale archive of remote sensing data and run geospatial analysis on Google's infrastructure. Currently, Google offers support only for Python and JavaScript. `rgee` fills the gap **by providing support for R!**. Below you will find the comparison between the syntax of `rgee` and the two other Google-supported client libraries.

<table>
<tr>
<th> JS (Code Editor) </th>
<th> Python </th>
<th> R </th>
</tr>
<tr>
<td>

``` javascript
var db = 'CGIAR/SRTM90_V4'
var image = ee.Image(db)
print(image.bandNames())
#> 'elevation'
```

</td>
<td>

``` python
import ee
ee.Initialize()
db = 'CGIAR/SRTM90_V4'
image = ee.Image(db)
image.bandNames().getInfo()
#> [u'elevation']
```

</td>
<td>

``` r
library(rgee)
ee_Initialize()
db <- 'CGIAR/SRTM90_V4'
image <- ee$Image(db)
image$bandNames()$getInfo()
#> [1] "elevation"
```
</td>
</tr>
</table>


**Quite similar, isn't it?**. However, additional more minor changes should be considered when using Google Earth Engine with R. Please check the [consideration section](https://r-spatial.github.io/rgee/articles/rgee02.html) before you start coding!

## Installation

Install from CRAN with:

``` r
install.packages("rgee")
```

Install the development versions from github with

``` r
library(remotes)
install_github("r-spatial/rgee")
```

Additionally, `rgee` depends on the [Python packages](https://rstudio.github.io/reticulate/articles/package.html): [numpy](https://pypi.org/project/numpy/) and [ee](https://pypi.org/project/earthengine-api/). To install them, users can follow any of these three methods:

  1.  Use [**ee_install**](https://r-spatial.github.io/rgee/reference/ee_install.html) (Highly recommended for users with no experience with Python environments)

``` r
rgee::ee_install()
```

2.  Use [**ee_install_set_pyenv**](https://r-spatial.github.io/rgee/reference/ee_install_set_pyenv.html) (Recommended for users with experience with Python environments)

``` r
rgee::ee_install_set_pyenv(
  py_path = "/home/csaybar/.virtualenvs/rgee/bin/python", # Change it for your own Python PATH
  py_env = "rgee" # Change it for your own Python ENV
)
```

Take into account that the Python PATH you set must have installations of the Earth Engine Python API and numpy. The use of **miniconda/anaconda is mandatory for Windows users,** Linux and MacOS users could also use virtualenv. See [reticulate](https://rstudio.github.io/reticulate/articles/python_packages.html) documentation for more details.

Another option, only possible for MacOS and Linux, is just set the Python PATH:

``` r
rgee::ee_install_set_pyenv(
  py_path = "/usr/bin/python3",
  py_env = NULL
)
```

However, [**rgee::ee_install_upgrade**](https://r-spatial.github.io/rgee/reference/ee_install_upgrade.html) and [**reticulate::py_install**](https://rstudio.github.io/reticulate/reference/py_install.html) will not work until you set a Python ENV.

3.  Use the Python PATH setting support that offer [Rstudio v.1.4 \>](https://blog.rstudio.com/2020/10/07/rstudio-v1-4-preview-python-support/). See this [tutorial](https://github.com/r-spatial/rgee/tree/help/rstudio/).

After install `Python dependencies` (and Restart R!!), you might use the function below for checking the status of rgee.

``` r
ee_check() # Check non-R dependencies
```

## Sync rgee with other Python packages

Integrate [rgee](https://r-spatial.github.io/rgee/) with [geemap](https://geemap.org/).

``` r
library(reticulate)
library(rgee)

# 1. Initialize the Python Environment
ee_Initialize()

# 2. Install geemap in the same Python ENV that use rgee
py_install("geemap")
gm <- import("geemap")
```

Upgrade the [earthengine-api](https://pypi.org/project/earthengine-api/)

``` r
library(rgee)
ee_Initialize()
ee_install_upgrade()
```

## Package Conventions

-   All `rgee` functions have the prefix ee\_. Auto-completion is your friend :).
-   Full access to the Earth Engine API with the prefix [**ee\$...**](https://developers.google.com/earth-engine/).
-   Authenticate and Initialize the Earth Engine R API with [**ee_Initialize**](https://r-spatial.github.io/rgee/reference/ee_Initialize.html). It is necessary once per session!.
-   `rgee` is "pipe-friendly"; we re-export %\>% but do not require its use.

## Hello World

### 1. Compute the trend of night-time lights ([JS version](https://github.com/google/earthengine-api/))

Authenticate and Initialize the Earth Engine R API.

``` r
library(rgee)
ee_Initialize()
```

Adds a band containing image date as years since 1991.

``` r
createTimeBand <-function(img) {
  year <- ee$Date(img$get('system:time_start'))$get('year')$subtract(1991L)
  ee$Image(year)$byte()$addBands(img)
}
```

Map the time band creation helper over the [night-time lights collection](https://developers.google.com/earth-engine/datasets/catalog/NOAA_DMSP-OLS_NIGHTTIME_LIGHTS/).

``` r
collection <- ee$
  ImageCollection('NOAA/DMSP-OLS/NIGHTTIME_LIGHTS')$
  select('stable_lights')$
  map(createTimeBand)
```

Compute a linear fit over the series of values at each pixel, visualizing the y-intercept in green, and positive/negative slopes as red/blue.

``` r
col_reduce <- collection$reduce(ee$Reducer$linearFit())
col_reduce <- col_reduce$addBands(
  col_reduce$select('scale'))
ee_print(col_reduce)
```

Create an interactive visualization!

``` r
Map$setCenter(9.08203, 47.39835, 3)
Map$addLayer(
  eeObject = col_reduce,
  visParams = list(
    bands = c("scale", "offset", "scale"),
    min = 0,
    max = c(0.18, 20, -0.18)
  ),
  name = "stable lights trend"
)
```

![rgee_01](https://user-images.githubusercontent.com/16768318/71565699-51e4a500-2aa9-11ea-83c3-9e1d32c82ba6.png)

### 2. Extract precipitation values

Install and load `tidyverse` and `sf` R packages, and initialize the Earth Engine R API.

``` r
library(tidyverse)
library(rgee)
library(sf)

ee_Initialize()
```

Read the `nc` shapefile.

``` r
nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
```

Map each image from 2001 to extract the monthly precipitation (Pr) from the [Terraclimate dataset](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE/)

``` r
terraclimate <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") %>%
  ee$ImageCollection$filterDate("2001-01-01", "2002-01-01") %>%
  ee$ImageCollection$map(function(x) x$select("pr")) %>% # Select only precipitation bands
  ee$ImageCollection$toBands() %>% # from imagecollection to image
  ee$Image$rename(sprintf("PP_%02d",1:12)) # rename the bands of an image
```

Extract monthly precipitation values from the Terraclimate ImageCollection through `ee_extract`. `ee_extract` works similar to `raster::extract`, you just need to define: the ImageCollection object (x), the geometry (y), and a function to summarize the values (fun).

``` r
ee_nc_rain <- ee_extract(x = terraclimate, y = nc["NAME"], sf = FALSE)
```

Use ggplot2 to generate a beautiful static plot!

``` r
ee_nc_rain %>%
  pivot_longer(-NAME, names_to = "month", values_to = "pr") %>%
  mutate(month, month=gsub("PP_", "", month)) %>%
  ggplot(aes(x = month, y = pr, group = NAME, color = pr)) +
  geom_line(alpha = 0.4) +
  xlab("Month") +
  ylab("Precipitation (mm)") +
  theme_minimal()
```

<p align="center">

  <img src="https://user-images.githubusercontent.com/16768318/81945044-2cbd8280-95c3-11ea-9ef5-fd9f6fd5fe89.png" width="80%"/>

  </p>

  ### 3. Create an NDVI-animation ([JS version](https://developers.google.com/earth-engine/tutorials/community/modis-ndvi-time-series-animation/))

  Install and load `sf`, after that, initialize the Earth Engine R API.

``` r
library(magick)
library(rgee)
library(sf)

ee_Initialize()
```

Define the regional bounds of animation frames and a mask to clip the NDVI data by.

``` r
mask <- system.file("shp/arequipa.shp", package = "rgee") %>%
  st_read(quiet = TRUE) %>%
  sf_as_ee()
region <- mask$geometry()$bounds()
```

Retrieve the MODIS Terra Vegetation Indices 16-Day Global 1km dataset as an `ee.ImageCollection` and select the NDVI band.

``` r
col <- ee$ImageCollection('MODIS/006/MOD13A2')$select('NDVI')
```

Group images by composite date

``` r
col <- col$map(function(img) {
  doy <- ee$Date(img$get('system:time_start'))$getRelative('day', 'year')
  img$set('doy', doy)
})
distinctDOY <- col$filterDate('2013-01-01', '2014-01-01')
```

Define a filter that identifies which images from the complete collection match the DOY from the distinct DOY collection.

``` r
filter <- ee$Filter$equals(leftField = 'doy', rightField = 'doy')
```

Define a join; convert the resulting FeatureCollection to an ImageCollection.

``` r
join <- ee$Join$saveAll('doy_matches')
joinCol <- ee$ImageCollection(join$apply(distinctDOY, col, filter))
```

Apply median reduction among matching DOY collections.

``` r
comp <- joinCol$map(function(img) {
  doyCol = ee$ImageCollection$fromImages(
    img$get('doy_matches')
  )
  doyCol$reduce(ee$Reducer$median())
})
```

Define RGB visualization parameters.

``` r
visParams = list(
  min = 0.0,
  max = 9000.0,
  bands = "NDVI_median",
  palette = c(
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  )
)
```

Create RGB visualization images for use as animation frames.

``` r
rgbVis <- comp$map(function(img) {
  do.call(img$visualize, visParams) %>%
    ee$Image$clip(mask)
})
```

Define GIF visualization parameters.

``` r
gifParams <- list(
  region = region,
  dimensions = 600,
  crs = 'EPSG:3857',
  framesPerSecond = 10
)
```

Get month names

``` r
dates_modis_mabbr <- distinctDOY %>%
  ee_get_date_ic %>% # Get Image Collection dates
  '[['("time_start") %>% # Select time_start column
  lubridate::month() %>% # Get the month component of the datetime
  '['(month.abb, .) # subset around month abbreviations
```

Use ee_utils_gif\_\* functions to render the GIF animation and add some texts.

``` r
animation <- ee_utils_gif_creator(rgbVis, gifParams, mode = "wb")
animation %>%
  ee_utils_gif_annotate(
    text = "NDVI: MODIS/006/MOD13A2",
    size = 15, color = "white",
    location = "+10+10"
  ) %>%
  ee_utils_gif_annotate(
    text = dates_modis_mabbr,
    size = 30,
    location = "+290+350",
    color = "white",
    font = "arial",
    boxcolor = "#000000"
  ) # -> animation_wtxt

# ee_utils_gif_save(animation_wtxt, path = "raster_as_ee.gif")
```

<p align="center">

  <img src="https://user-images.githubusercontent.com/16768318/77121867-203e0300-6a34-11ea-97ba-6bed74ef4300.gif"/>

  </p>

  ## How does rgee work?


  `rgee` is **not** a native Earth Engine API like the Javascript or Python client. Developing an Earth Engine API from scratch would create too much maintenance burden, especially considering that the API is in [active development](https://github.com/google/earthengine-api). So, how is it possible to run Earth Engine using R? the answer is [reticulate](https://rstudio.github.io/reticulate/). `reticulate` is an R package designed to allow seamless interoperability between R and Python. When an Earth Engine **request** is created in R, `reticulate` will translate this request into Python and pass it to the `Earth Engine Python API`, which  converts the request to a `JSON` format. Finally, the request is received by the GEE Platform through a Web REST API. The **response** will follow the same path in reverse.

![workflow](https://user-images.githubusercontent.com/16768318/71569603-3341d680-2ac8-11ea-8787-4dd1fbba326f.png)

## Quick Start User's Guide for rgee

<a href="https://bit.ly/35W0pwa"><img src="https://user-images.githubusercontent.com/16768318/88080933-5b1c8880-cb45-11ea-9546-f0640da13997.png" height="99"/></a> <a href="https://bit.ly/3iyxvYi"><img src="https://user-images.githubusercontent.com/16768318/86457619-8ef45300-bce9-11ea-9f08-7c1ee14071fb.png" height="100"/></a> <a href="https://ambarja.github.io/Handbook_rgee/pdf/vol01.pdf"><img src="https://user-images.githubusercontent.com/16768318/86457622-90be1680-bce9-11ea-92f0-78cfb915c4bc.png" height="101"/></a>

  **Created by:** - EN and POR: Andres Luiz Lima Costa <https://bit.ly/3p1DFm9> - SPA: Antony Barja Ingaruca <https://ambarja.github.io/>

  ## Code of Conduct

  Please note that the `rgee` project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

## Contributing Guide

👍 Thanks for taking the time to contribute! 🎉👍 Please review our [Contributing Guide](CONTRIBUTING.md).

## Share the love

Think **rgee** is useful? Let others discover it, by telling them in person via Twitter or a blog post.

Using **rgee** for a paper you are writing? Consider citing it

``` r
citation("rgee")
To cite rgee in publications use:

  C Aybar, Q Wu, L Bautista, R Yali and A Barja (2020) rgee: An R
  package for interacting with Google Earth Engine Journal of Open
  Source Software URL https://github.com/r-spatial/rgee/.

A BibTeX entry for LaTeX users is

@Article{,
  title = {rgee: An R package for interacting with Google Earth Engine},
  author = {Cesar Aybar and Quisheng Wu and Lesly Bautista and Roy Yali and Antony Barja},
  journal = {Journal of Open Source Software},
  year = {2020},
}
```

## Credits

We want to offer a **special thanks** :raised_hands: :clap: to [**Justin Braaten**](https://github.com/jdbcode) for his wise and helpful comments in the whole development of **rgee**. As well, we would like to mention the following third-party R/Python packages for contributing indirectly to the improvement of rgee:

-   [**gee_asset_manager - Lukasz Tracewski**](https://github.com/tracek/gee_asset_manager/)
-   [**geeup - Samapriya Roy**](https://github.com/samapriya/geeup/)
-   [**geeadd - Samapriya Roy**](https://github.com/samapriya/gee_asset_manager_addon/)
-   [**cartoee - Kel Markert**](https://github.com/KMarkert/cartoee/)
-   [**geetools - Rodrigo E. Principe**](https://github.com/gee-community/gee_tools/)
-   [**landsat-extract-gee - Loïc Dutrieux**](https://github.com/loicdtx/landsat-extract-gee/)
-   [**earthEngineGrabR - JesJehle**](https://github.com/JesJehle/earthEngineGrabR/)
-   [**sf - Edzer Pebesma**](https://github.com/r-spatial/sf/)
-   [**stars - Edzer Pebesma**](https://github.com/r-spatial/stars/)
-   [**gdalcubes - Marius Appel**](https://github.com/appelmar/gdalcubes/)
---
title: "NEWS"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{NEWS}
  %\VignetteEncoding{UTF-8}
---
# rgee 1.1.2.9000

- Better `Map$addLayer` support to COG.
- predefinedAcl='bucketLevel' is set as deault in stars_as_ee, sf_as_ee, raster_as_ee, and local_to_gcs.
- ee_utils_sak_validate and ee_utils_sak_copy added to ee_utils.R
- vignette to describe Shiny & rgee sync is added.
- vignette that describe how to set up a SaK.
- vignette that describe how to integrate rgee and markdown.


# rgee 1.1.2

- Fix an error in 'ee_check' warning message.
- Now 'ee_install' install Python 3.8 by default for Windows.
- Now users can control the access of buckets and objects. See predefinedAcl argument in local_to_gcs.

# rgee 1.1.1

- Deprecated.R file deleted.
- ee_help Rstudio addin critical bug solved.
- ee_clean argument name changed to 'user' rather than 'email'.
- New test inside ee_Initialize checks if the user token has enough permissions to read or modify files from GD.
- DESCRIPTION file modified: rgee always must use googledrive 2.0.0>=
- rstudioapi package moved from Suggests to Imports.

# rgee 1.1.0

- re-coded the Map and R6Map modules to simplify the maintenance. Many bugs were solved.
- ee_utils_get_crs now call to the web.archive.org if spatialreference is shut down.
- Math module and subsetting module were migrated to [rgeeExtra](https://github.com/r-earthengine/rgeeExtra).
- In ee_Initialize, the "email" parameter was renamed to "user".
- ee_get is now an internal function of rgee. A new/faster version of ee_get is available in rgeeExtra.
- Support to display COG resources. See Map or R6Map examples.
- rgee now supports googledrive version 2.0.0.
- Obtain COG metadata (`ee_utils_cog_metadata`).
- New test unit tests.
- Solve a bug in `ee_help`.
- `Map$addLegend` supports categorical legends in leaflet interactive maps.
- New logos :) 


# rgee 1.0.9

- Accessing the Earth Engine Data Catalog via '$' thanks to the [Earth-Engine-Datasets-List](https://github.com/samapriya/Earth-Engine-Datasets-List) created and supported by [@samapriya](https://github.com/samapriya).
- Math  functions (abs, sign, sqrt, ceiling, cummax, cummin, cumprod, cumsum, 
log, log10, log1p, acos, floor, asin, atan, exp, expm1, cos, cosh, sin, sinh, 
tan, and tanh.) to `ee$Image`.
- Summary functions (max, mean, min, range, sum, product) to `ee$Image`.
- Comparison operators (==, !=, >, <, <=, >=) to `ee$Image`.
- Logic operators (!, &, |) to `ee$Image`.
- Arithmetic operators (+, -, *, /, ^, %%, %/%) to `ee$Image`.
- Subsetting operators ('[[<-', '[[') to `ee$Image` and `ee$ImageCollection`.
- GH Action to automatically updated the Earth Engine Python API.
- `ee_as_sf(..., via = "getInfo")` does not write in temp folder.
- `ee_as_sf` now returns by default a GeoJSON instead of a ESRI shapefile.
- When EarthEngineMaps have the same name, a random hex string is added to the second map.
- Fix a bug in `sf_as_ee` that add `id` column to the results.
- `ee_extract` now supports lazy evaluation and containers `drive` and `gcs`.
- R6, class to display Earth Engine (EE) spatial objects, added.
- Map is now a `R6` object instead of a environment.
- Vignettes documentation upgrade.
- `ee_install_set_pyenv` support local .Renviron. 
- Fix a bug for new tokens in {googledrive} #139.
- Fix a bug in ee_print that sometimes make see the warning: *ee_utils_py_to_r(.) : restarting interrupted promise evaluation*.
- Fix a bug in `ee_install_set_pyenv` when py_env=NULL (#118, thanks @MatthieuStigler).
- GH Action test-coverage removed.
- Fix a bug in `ee_extract` that changes column names when starts with numbers (#119, thanks @joshualerickson).


# rgee 1.0.8

- Unit testing enhanced.
- Fix a bug in Map\$addLayer(..., legend=TRUE) when `eeobject` is an constant image  (i.e. ee\$Image(0)).
- Stop message in `ee_Initialize` to does not allow the use of rgee when the folder ".../rgee/python/" does not exist.
- Info messages when rgee make changes to.Renviron.
- Earth Engine Python API updated to 0.1.247.

# rgee 1.0.7

- Unit testing enhanced.
- More documentation related to credentials.
- Smoother connection with Python (reticulate).
- Now Map$... functions only depend on {leaflet}.
- Public argument added to `ee_as_sf`, `ee_as_raster`, `ee_as_stars`, `ee_imagecollection_to_local`, `ee_drive_to_local` and `ee_gcs_to_local` which permit to create a public link to the resources generated.
- A new argument "**metadata**" is added to `ee_as_sf`, `ee_as_raster`, `ee_as_stars`, `ee_drive_to_local`, `ee_imagecollection_to_local`, and `ee_gcs_to_local`. If TRUE, the metadata related to the export of the images will be added to raster/stars objects.
- Fix a bug in Rstudio `ee_help` addins.
- Fix a bug in `ee_extract` which adds the `system:index` to the colnames when the `x` argument is an `ee$ImageCollection`. 
- Fix a bug that does not permit to `ee_as_raster` and `ee_as_stars` change the fileNamePrefix (#112).
- a stop added in `sf_as_ee` since {geojsonio} does not support POSIXt objects (#113).
- Lazy evaluation support to `ee_imagecollection_to_local`, `ee_as_sf`, `ee_as_raster` and `ee_as_stars`.
- Export images via 'getInfo' was removed from `ee_as_raster` and `ee_as_stars` to avoid problems related to geometric offset.
- Now `ee_monitoring` can also be invoked with the ID of a EE task started.
- `ee_search` module deprecated, it will be removed of rgee in version 1.0.8.
- New functions: `ee_utils_search_display` that display the website related to the Earth Engine dataset, and `ee_utils_future_value` that helps to run a {future} container.
- Earth Engine Python API test updated to 0.1.246.

# rgee 1.0.6
- Class method chaining (i.e. `x$size()$getInfo()`) were changed by pipes (i.e. ee_x %>% `ee$FeatureCollection$size() %>% ee$Number()`) in all the `rgee` functions. This solve the problem "OverflowError: python int too large to convert to C long" on Window systems.
- rgee functions has a cleaner method to run system processes, {**processx**} 
instead of **base::system**. 
- `rgee` I/O functions now check argument before to start to upload/download data.
- Map operators (**+** and **|**) now support EarthEnginemap objects with the 
same name.
- Now `Map$addLayers` only display the legend of the first image.
- Fix a bug in `rgee:::ee_image_local` which makes do not work when all bands have not the same crs and crsTransform.
- "getInfo" method in download raster functions was deprecated and will be removed in v.1.0.8.
- Fix a bug in `sf_as_ee` and `ee_as_sf` now both support SR-ORG CRS codes.
- `ee_users` returns a data.frame.
- `ee_monitoring` counts the processing time.
- Fix a bug in `ee_utils_gif_creator` which makes don't work in windows.
- Several changes in `ee_extract`, now is faster and code is cleaner.
- Fix a bug in name creator in `ee_imagecollection_to_local`.
- A new message more detailed when the Python path does not have the earth-engine Python API.
- Earth Engine Python API updated to 0.1.235.

# rgee 1.0.5
- Important changes in the low level API to upload raster and vector with GCS. However, high upload API (`sf_as_ee`, `stars_as_ee`, and `raster_as_ee`) continue working in the same way.
- Add the functions: `ee_utils_create_manifest_image` and `ee_utils_create_manifest_table`
to help users to create a JSON file with all the upload parameters ("manifest", see https://developers.google.com/earth-engine/image_manifest/).
- Add the functions: `ee_utils_gif_creator`, `ee_utils_gif_annotate` and
`ee_utils_gif_save` to help users to read, add text, and save gif files.
- New | operator inspired in mapview 2.9.0! try: m4 | m5
- Fix several typos.
- Earth Engine Python API updated to 0.1.232.

# rgee 1.0.4
- Add `ee_help` a new Rstudio addins that mimics the help Rstudio interface (F1).
- Fix a bug that makes that `ee_as_sf` only supports `GeoJSON` format.
- If `dsn` is not specified in `ee_as_sf`, it will create a temporary shapefile (in \tmp dir).
- Fix a bug in `ee_imagecollection_to_local` (#87 Thanks @cedlfc44)
- Fix a bug in `ee_image_local` (#88 Thanks @cedlfc44)
- Fix a bug in `ee_create_credentials_drive` (#90 #78 Thanks @cedlfc44)

# rgee 1.0.3
- getPass library removed from `ee_Initialize`.
- New argument `display` in `ee_Initialize` to return the authentication URI. Useful for `rgee` colab users.
- Changes in some diagnostic messages to make possible to use `rgee` in colab.
- `ee_help` returns a HTML file rather than TRUE. It also now supports characters (e.g. `ee_help("ee$Image")`).
- Fix a strange bug when `ee_Initialize` tries to connect to reticulate the first time.
- Fix small bugs in `ee_user_info` and `ee_users`

# rgee 1.0.2
- Earth Engine Python API updated to 0.1.229.
- Fix a bug in `ee_Initialize`, that does not permit users to use `ee_createAssetHome` to define their *Earth Engine Assets home root folder*

# rgee 1.0.1
- Earth Engine Python API updated to 0.1.228.
- `ee_Initialize` now set the global env "RETICULATE_PYTHON" rather than .onLoad

# rgee 1.0.0
- We implement `ee_createAssetHome` to help users to define their *Earth Engine Assets home root folder* without leaving ee_Initialize(). (#70 Thanks @jhollist)
- Fix a bug in `ee_Initialize(drive = TRUE, gcs = TRUE)` which do not permit users save credentials. (#72 Thanks @appelmar).
- Removed `check_ring_dir` argument from `sf_as_ee`. Now `rgee` expect that users fix potential geometry problems before upload geometries to their Earth Engine assets.
- Changes in "welcome to rgee"  message (located in `ee_Initialize`). We delete `stop("Initialization aborted")` and implement `ee_search_init_message`. It will permit to `rgee` to know if the user accepted to create the global var "EARTHENGINE_INIT_MESSAGE" without need to restart the R session.
- Fix minor bugs in `sf_as_ee` and `gcs_to_ee_table`. The argument `command_line_tool_path` was added to give users the option to set the path of the Earth Engine command linetool. This new argument is only relevant to upload files using Google Cloud Storage. New diagnostic message were added.
- Now `sf_as_ee` returns an `ee$Geometry$...` when `x` is a `sfg` or a single sfc object. If `x` is a sf or a `sfc` object with multiple geometries will return an `ee$FeatureCollection`. New unit test for `sf_as_ee`.
- Changes in the documentation of `ee_as_stars` and `ee_as_raster` (#72 Thanks @appelmar).
- Fix minor bugs in `raster_as_ee`, `stars_as_ee` and `gcs_to_ee_image`. The argument `command_line_tool_path` was added to give users the option to set the path of the Earth Engine command linetool. New diagnostic message were added. New unit test added.
- Fixed a bug in `stars_as_ee` that didn't allow to read single-band images.
- `ee_manage_asset_access` has a better approach to determine the user owner of the asset.
- Add a new logical argument called 'strict' to `ee_manage_assetlist`, `ee_manage_copy`,
`ee_manage_move`, `ee_manage_set_properties` and `ee_manage_delete_properties`. If TRUE, 
the existence of the asset will be evaluate before to perform the task. By default TRUE.
- If the `path_asset` is not specified in `ee_manage_assetlist`, rgee will assume that the
path_asset is `ee_get_assethome`.
- `raster_as_ee.R` was created in the /R folder to maintain an order between functions
to upload and download images.
- Fix a bug in the documentation of `print.ee.computedobject.ComputedObject()` (#72 Thanks @appelmar).
- Fix a bug in `ee_install_set_pyenv` now users can set `py_path` without set 
`py_env` and  vice versa.
- `ee_extract` was adapted to work well with changes in `sf_as_ee`.
- R CMD check is more friendly with users, the use of `--run-dontrun` is also 
available (#72 Thanks @appelmar).
- Fix a bug in `ee_get_date_ic`.
- `Map$addLayer` now could display a legend when `eeObject` is an one-band `ee$Image` (By default legend = TRUE).
- New function `ee_get`: Return the element at the specified position in a Earth Engine Collection.
- New function `Map$addLayers`: Create interactive visualizations of ImageCollections.
- Fix a bug: "OverflowError: Python int too large to convert to C long" in `ee_monitoring` (#79 Thanks @andreatitolo).
- Earth Engine Python API updated to 0.1.227.

# rgee 0.6.2
- Earth Engine Python API updated to 0.1.226.
- Fix some typos.
- Fix a minor bug in ee_monitoring.
- Users can mix mapview and EarthEnginemap objects in the same pipeline (see Examples in `Map$addLayer`).
- Add `ee_as_mapview`, a function to convert `EarthEnginemap` objects to  `mapview` objects.
- add a new logical argument called 'strict' to **ee_manage_delete**. If TRUE, the existence of the asset will be evaluate before to perform the task.
- Fix a bug in ee_Initialize, now users without an Earth Engine Assets home root will see a message.
- Fix a minor bug when ee_Initialize change of user, now before to change of user the GCS and GD credentials will be deleted.
- ee_check completely renovated. 
  - New diagnostic messages.
  - Fix a minor bug when testing GCS credentials.
  - The file ee_check.py was deleted.
- Roy Samapriya added as a contributor.  

# rgee 0.6.1
- Fix some typos.
- rgee website update.
- Add citation package option .
- Additional export arguments add to ee_as_stars,ee_as_raster,  ee_imagecollection_to_local and ee_as_sf.

# rgee 0.6.0 
- Earth Engine Python API updated to 0.1.225.
- Fix some typos.
- DESCRIPTION: Moving leaflet, mapview, geojsonio, sf and stars from Import to Suggest. Now users with installation problems can equally use the Earth Engine API although with less operability.
- The 'EarthEngineMap' S4 class was created to avoid incompatibilities with mapview.
- Fix a critical bug in **ee_install** due to the lack of breaks in repeat bucles.
- New function **ee_install_upgrade**.
- New global environment EARTHENGINE_PYTHON was created to help **ee_install_upgrade** to identify the Python environment used by rgee.

# rgee 0.5.4
- Earth Engine Python API updated to 0.1.224.
- Fix a Map typo.
- Fix a bug in ee_as_thumbnail, now the vizparams are checked before to pass to ee$Image$getThumbURL(...).
- ee is now an internal rgee environment.
- ee_reattach was deleted.
- ee_print now display the ee$Image properties: system:id, system:time_start and system:time_end.
- rgee now asks users if they would like to save EARTHENGINE_PYTHON in the .Renviron.

# rgee 0.5.3
- Fix a bug in ee_check_python_packages.
- \donttest changed by \dontrun in the documentation.
- gdal and v8 system dependencies added to GH actions.
- Fix a bug in the paper (`ee_extract` example).
- Fix a minor bug in `ee_extract`, now the argument ... works.
- `ee_table_to_drive`: changed the initial value of the argument folder from NULL to "rgee_backup".
- Fix a minor bug in `ee_monitoring`.
- Fix a minor bug in `ee_manage_cancel_all_running_task`.
- Fix a minor bug in `ee_manage_cancel_all_running_task`.
- Improvement in the documentation of `ee_install`.
- Changes in vignettes, `Best Practices` vignette added.

# rgee 0.5.2
- DESCRIPTION: single quotes in title and description.
- DESCRIPTION: A more compressible description of what rgee does.
- DESCRIPTION: Added web reference to the Earth Engine API.
- \dontrun changed by \donttest in all our examples.
- Added "#' @return" to several functions.
- Added "#' @family" to all the functions.
- Added 'quiet' argument to all the functions that needed.
- Added new contributors to rgee (Kevin Ushey, Tim Appelhans, JJ Allaire, Yuan Tang).
- New environmental variable for rgee "EARTHENGINE_INIT_MESSAGE". It will be used to display a message to new users.
- Earth Engine Python API updated to 0.1.223.
- Documentation updated for ee_print and ee_manage_*.
- Fix a bug in ee_install_set_pyenv that did not permit to create properly
the .Renviron_backup file.

# rgee 0.5.1
- ee_install_* functions were deprecated and replaced by ee_install. ee_install create an isolated Python virtual environment with all rgee dependencies.
- ee_install_set_pyenv can be used to set the EARTHENGINE_PYTHON variable.

# rgee 0.5.0
- **First attempt to submit to CRAN**
- Several typos fixed.
- rgee paper added.
- GitHub actions for automated testing and build the website.
- Due the changes in latest `reticulate` version (1.1.5), the functions `ee_install_earthengine_upgrade` and `ee_install_python_packages` were deprecated and both will remove in rgee 0.5.3.
- Config/reticulate added to DESCRIPTION file.
- .onLoad IMPORTANT CHANGES: Now `rgee` set EARTHENGINE_PYTHON instead of RETICULATE_PYTHON directly and RETICULATE_PYTHON_ENV is no longer required. This change will permit users to avoid problems with other R packages that use Python in the backend (such as tensorflow or keras).
- `ee_search_display` function added.
- Several typos fixed in all the documentation.
- Minor changes in `ee_as_sf` to support ee$FeatureCollections without elements.
- `data.colec.fbf` eliminated from all the examples.
- `rgee` now pass all `goodpractice` checks.
- `ee_get_img_date` and `ee_get_ic_date` are now `ee_get_date_img` and `ee_get_date_ic`.
- A new group of functions was created at `ee_utils.R`. `ee_pyfunc` is now `ee_utils_pyfunc`.
- Added new examples in `README.R`.
# Contributor Covenant Code of Conduct

<!-- This CODE_OF_CONDUCT.md is adapted from 
https://github.com/ropensci/dbparser/commit/474e50a10bb1e3e23ce6e5ec703ca7ba8ebd4adf#diff-a1ee87dafebc22cbd96979f1b2b7e837/-->

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education, socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

  * Using welcoming and inclusive language
  * Being respectful of differing viewpoints and experiences
  * Gracefully accepting constructive criticism
  * Focusing on what is best for the community
  * Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

  * The use of sexualized language or imagery and unwelcome sexual attention or
advances
  * Trolling, insulting/derogatory comments, and personal or political attacks
  * Public or private harassment
  * Publishing others' private information, such as a physical or electronic address, without explicit permission
  * Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at csaybar@gmail.com. All complaints will be reviewed and investigated and will result in a response that is deemed necessary and appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html 

[homepage]: https://www.contributor-covenant.org/

For answers to common questions about this code of conduct, see https://www.contributor-covenant.org/faq/
# Contributing to rgee

<!-- This CONTRIBUTING.md is adapted from https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c/ -->

First of all, thanks for considering contributing to rgee! 👍 It's people like you that make it rewarding for us - the project maintainers - to work on rgee. 😊

rgee is an open source project, maintained by people who care. We are not directly funded to do so.

[repo]: https://github.com/r-spatial/rgee/
[issues]: https://github.com/r-spatial/rgee/issues/
[new_issue]: https://github.com/r-spatial/rgee/issues/new/
[website]: https://r-spatial.github.io/rgee/
[citation]: https://r-spatial.github.io/rgee/authors.html
[email]: csaybar@gmail.com

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Share the love ❤️

Think **rgee** is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using **rgee** for a paper you are writing? Consider [citing it][citation].

### Ask a question ⁉️

Using our_package and got stuck? Browse the [documentation][website] to see if you can find a solution. Still stuck? Post your question as an [issue on GitHub][new_issue]. While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [email][email].

### Propose an idea 💡

Have an idea for a new rgee feature? Take a look at the [documentation][website] and [issue list][issues] to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub][new_issue]. While we can't promise to implement your idea, it helps to:

  * Explain in detail how it would work.
  * Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a bug 🐛

Using rgee and discovered a bug? That's annoying! Don't let others have the same experience and report it as an [issue on GitHub][new_issue] so we can fix it. A good bug report makes it easier for us to do so, so please include:

  * Your operating system name and version (e.g. Mac OS 10.13.6).
  * Any details about your local setup that might be helpful in troubleshooting.
  * Detailed steps to reproduce the bug.

### Improve the documentation 📖

Noticed a typo on the website? Think a function could use a better example? Good documentation makes all the difference, so your help to improve it is very welcome!

#### The website

[This website][website] is generated with [`pkgdown`](https://pkgdown.r-lib.org/). That means we don't have to write any html: content is pulled together from documentation in the code, vignettes, [Markdown](https://guides.github.com/features/mastering-markdown/) files, the package `DESCRIPTION` and `_pkgdown.yml` settings. If you know your way around `pkgdown`, you can [propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to improve documentation. If not, [report an issue][new_issue] and we can point you in the right direction.

#### Function documentation

Functions are described as comments near their code and translated to documentation using [`roxygen2`](https://klutometis.github.io/roxygen/). If you want to improve a function description:

1. Go to `R/` directory in the [code repository][repo].
2. Look for the file with the name of the function.
3. [Propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to update the function documentation in the roxygen comments (starting with `#'`).

### Contribute code 📝

Care to fix bugs or implement new functionality for our_package? Awesome! 👏 Have a look at the [issue list][issues] and leave a comment on the things you want to work on. See also the development guidelines below.

## Development guidelines

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo][repo] and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream master`.
3. Open the RStudio project file (`.Rproj`).
5. Make your changes:
  * Write your code.
* Test your code (bonus points for adding unit tests).
* Document your code (see function documentation above).
* Check your code with `devtools::check()` and aim for 0 errors and warnings.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request/).
---
title: 'rgee: An R package for interacting with Google Earth Engine'
bibliography: paper.bib
date: \today
output: pdf_document
tags:
  - R
  - Earth Engine
  - Earth Observation
  - spatial analysis
authors:
  - name: Cesar Aybar
    orcid: 0000-0003-2745-9535
    affiliation: 1
  - name: Qiusheng Wu
    affiliation: 2
    orcid: 0000-0001-5437-4073
  - name: Lesly Bautista
    affiliation: 3
    orcid: 0000-0003-3523-8687
  - name: Roy Yali
    affiliation: 3
    orcid: 0000-0003-4542-3755
  - name: Antony Barja
    affiliation: 3
    orcid: 0000-0001-5921-2858    
affiliations:
  - name: Department of Geoinformatics – Z_GIS, University of Salzburg, Austria
    index: 1
  - name: Department of Geography, University of Tennessee, Knoxville, TN 37996, USA
    index: 2
  - name: Universidad Nacional Mayor de San Marcos, Lima, Lima 15081, Peru
    index: 3
---

# Summary
Google Earth Engine [@gorelick2017google] is a cloud computing platform designed for planetary-scale environmental data analysis. Its multi-petabyte data catalog and computation services are accessed via an Internet-accessible API. The API is exposed through JavaScript and Python client libraries. Google provides a browser-based IDE for the JavaScript API, and while convenient and useful for rapid data exploration and script development, it does not allow third-party package integration, relying solely on Google Maps and Google Charts for data visualization, and proprietary systems for metadata viewing and asset management. In contrast, the Python and Node.js distributions offer much flexibility for developers to integrate with third-party libraries. However, without the structure of a dedicated IDE, casual users can be left directionless and daunted. A significant gap exists between these two offerings (Google-supported JavaScript IDE and base client libraries) where convenience and flexibility meet. We propose to fill this gap with an R package that wraps the Earth Engine Python API to provide R users with a familiar interface, rapid development features, and flexibility to analyze data using open-source, third-party packages.

rgee is an Earth Engine (EE) client library for R that allows users to leverage the strengths of the R spatial ecosystem and Google Earth Engine in the same workflow. All of the Earth Engine Python API classes, modules, and functions are made available through the reticulate package [@reticulate], which embeds a Python session within an R session, enabling seamless interoperability. Additionally, rgee adds several new features such as (i) new I/O design, (ii) interactive map display,  (iii) easy extraction of time series, (iv) asset management interface, and (v) metadata display. In addition, rgee also makes it possible to execute Earth Engine Python code from within R, making the translation of large Python projects unnecessary.

# Features

## Enhanced I/O 

rgee implements several functions to support download/upload of spatial objects (Table \ref{table:1} and Table \ref{table:2}). For instance, to download vector (image) files one can use `ee_as_sf` (`ee_as_raster` or `ee_as_stars`). In rgee, all the functions from server to local side have the option to fetch data using an intermediate container (Google Drive or Google Cloud Storage) or through a REST call ("\$getInfo"). Although the latter option performs a quick download, there is a request limit of 262144 pixels for `ee$Image` and 5000 elements for `ee$FeatureCollection` which makes it unsuitable for large objects. Other download functions, from server-side to others (see Table \ref{table:1}), are implemented to enable more customized download workflows. For example, using `ee_image_to_drive` and `ee_drive_to_local` users could create scripts which save results in the `.TFRecord` rather than the `.GeoTIFF` format. The upload process follows the same logic (Table \ref{table:2}). rgee includes `raster_as_ee` and `stars_as_ee` for uploading images and `sf_as_ee` for uploading vector data. Large uploads are only possible with an active Google Cloud Storage account.


|         	|                   	|      FROM      	|       TO      	|       RETURN       	|
|---------	|-------------------	|:--------------:	|:-------------:	|:------------------:	|
| Image   	| ee_image_to_drive 	|    EE server   	|     Drive     	|   Unstarted task   	|
|         	| ee_image_to_gcs   	|    EE server   	| Cloud Storage 	|   Unstarted task   	|
|         	| ee_image_to_asset 	|    EE server   	|    EE asset   	|   Unstarted task   	|
|         	| ee_as_raster      	|    EE server   	|     Local     	| RasterStack object 	|
|         	| ee_as_stars       	|    EE server   	|     Local     	| Proxy-stars object 	|
| Table   	| ee_table_to_drive 	|    EE server   	|     Drive     	|   Unstarted task   	|
|         	| ee_table_to_gcs   	|    EE server   	| Cloud Storage 	|   Unstarted task   	|
|         	| ee_table_to_asset 	|    EE server   	|    EE asset   	|   Unstarted task   	|
|         	| ee_as_sf          	|    EE server   	|     Local     	|      sf object     	|
| Generic 	| ee_drive_to_local 	|      Drive     	|     Local     	|   object filename  	|
|         	| ee_gcs_to_local   	|  Cloud Storage 	|     Local     	|     GCS filename  	|

: Download functions provided by the rgee package. \label{table:1}


|         	|                 	|      FROM     	|       TO      	|            RETURN           	|
|---------	|-----------------	|:-------------:	|:-------------:	|:---------------------------:	|
| Image   	| gcs_to_ee_image 	| Cloud Storage 	|    EE asset   	|          EE Asset ID       	  |
|         	| raster_as_ee    	|     Local     	|    EE asset   	|          EE Asset ID       	  |
|         	| stars_as_ee     	|     Local     	|    EE asset   	|          EE Asset ID       	  |
| Table   	| gcs_to_ee_table 	| Cloud Storage 	|    EE asset   	|          EE Asset ID       	  |
|         	| sf_as_ee        	|     Local     	|    EE asset   	|          EE Asset ID       	  |
| Generic 	| local_to_gcs    	|     Local     	| Cloud Storage 	|         GCS filename        	|

: Upload functions provided by the rgee package. \label{table:2}

The following example illustrates how to integrate the rgee I/O module and ggplot2 [@wickham2011ggplot2] to download and visualize metadata for the [BLM AIM TerrestrialAIM](https://developers.google.com/earth-engine/datasets/catalog/BLM_AIM_v1_TerrADat_TerrestrialAIM#description/) dataset.

```r
library(tidyverse)
library(rgee)
library(sf)

ee_Initialize()

# Define a Region of interest
roi <- ee$Geometry$Point(-120.06227, 40.64189)$buffer(25000)

# Load TerrADat TerrestrialAIM Dataset
blocks <- ee$FeatureCollection("BLM/AIM/v1/TerrADat/TerrestrialAIM")
subset <- blocks$filterBounds(roi)

# Move an Earth Engine FeatureCollection to their local env
sf_subset <- ee_as_sf(x = subset)

# Create a boxplot with ggplot2
gapPct <- c("_25_50" = "GapPct_25_50","_51_100"="GapPct_51_100",
            "101_200" = "GapPct_101_200","200_>" = "GapPct_200_plus")

sf_subset[gapPct] %>% 
  st_set_geometry(NULL) %>% 
  as_tibble() %>% 
  rename(!!gapPct) %>% 
  pivot_longer(seq_along(gapPct), names_to = "Range") %>% 
  ggplot(aes(x = Range, y = value, fill = Range)) +
  geom_boxplot() +
  xlab("") + ylab("% of the plot's soil surface") +
  theme_minimal()
```
![Gaps percentage between plant canopies of different sizes in a place near to Carson City, Nevada, USA. \label{fig:AIM}](rgee_paper_00.png){ width=65% }

## Interactive Map Display
rgee offers interactive map display through  `Map$addLayer`, an R function mimicking the mapping module of the Earth Engine JavaScript Code Editor. `Map$addLayer` takes advantage of the `getMapId` EE method to fetch and return an ID dictionary being used to create layers in a mapview [@mapview] object. Users can specify visualization parameters to `Map$addLayer` by using the visParams argument, as demostrated below:

```r
library(rgee)
ee_Initialize()

# Load an ee$Image
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")

# Centers the map view
Map$centerObject(image)

# Display the ee$Image
Map$addLayer(
  eeObject = image, 
  visParams = list(bands = c("B4", "B3", "B2"), max = 10000), 
  name = "SF"
)
```

![Landsat 8 false color composite of San Francisco bay area, California, USA.](rgee_paper_mapview.png){ width=70% }

## Extraction of time series

rgee can extract values from `ee$Image` and `ee$ImageCollection` objects at a certain location based on `ee$Geometry`, `ee$Feature`, `ee$FeatureCollection` and `sf` objects. If the geometry is a polygon, users can summarize the values using built-in Earth Engine reducer functions. The code below explains how to extract the average areal rainfall from North Carolina counties using the [TerraClimate](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE/) dataset.

```r
library(ggplot2)
library(tidyr)
library(dplyr)
library(rgee)
library(sf)

ee_Initialize()

# Filter the terraclimate dataset by dates, reproject
# and select only the band "pr".
terraclimate <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")$
  filterDate("2001-01-01", "2002-01-01")$
  map(function(x) x$reproject("EPSG:4326")$select("pr"))

# Define a geometry
nc <- st_read(system.file("shape/nc.shp", package = "sf"))

# Extract the average areal rainfall
ee_nc_rain <- ee_extract(terraclimate, nc, sf = FALSE)
colnames(ee_nc_rain) <- sprintf("%02d", 1:12)
ee_nc_rain$name <- nc$NAME

# Create a data frame in a tidy format and display rainfall values
ee_nc_rain %>%
  pivot_longer(-name, names_to = "month", values_to = "pr") %>%
  ggplot(aes(x = month, y = pr, group = name, color = pr)) +
  geom_line(alpha = 0.4) +
  xlab("Month") +
  ylab("Precipitation (mm)") +
  theme_minimal()
```

![Average areal rainfall in counties of North Carolina for the year 2001 according to the TerraClimate dataset.](rgee_paper_01.png){ width=80% }


## Asset Management Interface

rgee implements an interface to batch actions on assets extending capabilities of the existing EE data module (`ee$data$*`). The interface allows users to create and eliminate folders, move and copy assets, set and delete properties, handle access control lists, and manage or cancel tasks. For example, users can copy a Landsat 8 image to their personal EE assets as follows:

```r
library(rgee)
ee_Initialize()

server_path <- "LANDSAT/LC08/C01/T1/"
user_asset_path <- ee_get_assethome()

ee_manage_copy(
  path_asset = paste0(server_path,"/LC08_044034_20140318"),
  final_path = paste0(user_asset_path,"/LC08_044034_20140318")
)
```

## Metadata display
The `ee_print` function can save and display all metadata related to EE spatial objects. With `ee_print`, users can retrieve information about the number of images or features, number of bands or geometries, number of pixels, geotransform, datatype, properties and approximate object size. `ee_print` can be used inside debugging pipelines (e.g. linking with `ee$Image$aside`).

```r
library(rgee)

ee_Initialize()
l8 <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")
ee_print(l8)
```
![Metadata for a Landsat 8 Image.](rgee_paper_02.png){ width=90% }

# Availability

rgee is an open-source software package made available under the Apache 2.0 license. It can be installed through GitHub repository using the remotes package: remotes::install_github("r-spatial/rgee"). A series of examples for using rgee are available at [https://r-spatial.github.io/rgee/](https://r-spatial.github.io/rgee/).

# Acknowledgments

The authors would like to thank Justin Braaten for his reviewing and helpful comments during the preparation of this manuscript and development of rgee.

# References
The paper can be compiled locally using the following command:

pandoc paper.md --bibliography paper.bib -o rgee_paper.pdf* rgee version:
* R version:
* Operating System:

#### At submit an issue, please attached the following information of your `rgee` session:

- [ ] You have the Python API installed (from terminal):
```
earthengine -h
```

- [ ] You can find the credentials file on your system: 
```r
library(rgee)
ee_path <- path.expand("~/.config/earthengine/credentials")
file.exists(ee_path)
```
- [ ] You can run a simple EE command from R: 

```r
library(rgee)

# Initialize the Earth Engine module.
ee_Initialize()

# Print metadata for a DEM dataset.
print(ee$Image('USGS/SRTMGL1_003')$getInfo())
``` 

Attach your Python (reticulate) configuration:

```r
library(reticulate)
py_config()
```

### Description

Describe what you were trying to get done.
Tell us what happened, what went wrong, and what you expected to happen.

### What I Did

```
Paste the command(s) you ran and the output.
If there was a crash, please include the traceback here.
```
