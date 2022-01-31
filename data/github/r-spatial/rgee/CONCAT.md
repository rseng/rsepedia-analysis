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
---
title: "**Introduction to rgee**"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{1. Introduction}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}    
---

<link rel="preconnect" href="https://fonts.gstatic.com">
<link href="https://fonts.googleapis.com/css2?family=Roboto&display=swap" rel="stylesheet">
<style type="text/css">
  body{
    font-size: 15pt;
    font-family: 'Roboto', sans-serif;
  }
  pre code{
    font-size: 12pt;
  }
  .list-group-item:last-child{
    font-size: 11pt;
    font-weight: bold;
  }
</style>

## **1. What is Google Earth Engine ?**

**Google Earth Engine** is a computing platform that allows users to run geospatial analysis on Google's infrastructure. There are several ways to interact with the platform:

 * Explorer
 * Code Editor
 * Javascript client library
 * Python client library
 * **R client library**   

This website is focused on the last one, you can use the R client library to send/receive messages to/from the Earth Engine server and develop **[web applications](https://github.com/MVanDenburg92/RGEE_Shiny/)**. 


```{r, setup, include=FALSE}
knitr::opts_chunk$set(
  comment = '', fig.width = 6, fig.height = 6
)
```


## **2. The purpose of Earth Engine is to:**

 * Perform highly interactive algorithm development at global scale

 * Push the envelope for big data in remote sensing

 * Enable high-impact, data-driven science

 * Make substantive progress on global challenges that involve large geospatial datasets

## **3. Components:**

The main components of Earth Engine are:

**Datasets**: A petabyte-scale archive of publicly available remotely sensed imagery and other data. [Explore the data catalog](https://developers.google.com/earth-engine/datasets).

**Compute power**: Google’s computational infrastructure optimized for parallel processing of geospatial data.

**WEB REST API/Client libraries**: For making requests to the Earth Engine servers.

## 4. Meet Earth Engine 

The Earth Engine API and advanced Earth Engine functionality **are experimental and subject to change**. Access is limited and requires requesting access via the **[form](https://earthengine.google.com/signup/)**. See **[Earth Engine official website](https://earthengine.google.com/)** to obtain more information.

<br>

<center>
  <iframe name="Stack" src="https://www.youtube.com/embed/gKGOeTFHnKY/" style='height: 450px; width: 80%;
  'frameborder="0" scrolling="no" id="iframe">
  </iframe>
</center>

## 5. Why rgee instead of code editor (Javascript)?

A short comparison based on **[Tyler Erickson presentation](https://docs.google.com/presentation/d/1MVVeyCdm-FrMVRPop6wB3iyd85TAlwB-F9ygTQZ8S1w/pub?slide=id.g1e419debf0_1_205/)**.

<center>
<table style="border-collapse:collapse;border-color:#ccc;border-spacing:0;border:none" class="tg"><thead><tr><th style="background-color:#f0f0f0;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;font-weight:bold;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal">Code Editor</th><th style="background-color:#f0f0f0;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;font-weight:bold;overflow:hidden;padding:10px 5px;text-align:center;vertical-align:top;word-break:normal">R</th></tr></thead><tbody><tr><td style="background-color:#fff;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">Easy to get started</td><td style="background-color:#f9f9f9;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">Easy to share code between scripts.</td></tr><tr><td style="background-color:#fff;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">Trivial to share scripts</td><td style="background-color:#f9f9f9;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">Easier transition to a web application (<span style="font-weight:bold">Shiny</span>).</td></tr><tr><td style="background-color:#fff;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">Built in authentication</td><td style="background-color:#f9f9f9;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">An I/O API more friendly with R users.</td></tr><tr><td style="background-color:#fff;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal"><span style="font-weight:bold;color:#FE0000">Limited input/output functionality</span></td><td style="background-color:#f9f9f9;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal">Many, many plotting options</td></tr><tr><td style="background-color:#fff;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal"><span style="font-weight:bold;color:#FE0000">Integration with other JS libraries is not possible</span></td><td style="background-color:#f9f9f9;border-color:inherit;border-style:solid;border-width:0px;color:#333;font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:10px 5px;text-align:left;vertical-align:top;word-break:normal"><span style="font-weight:bold;color:#FE0000">Some assembly (&amp; maintenance) required!.</span></td></tr></tbody></table>
</center>

## 6. Installation

**rgee** depends on **[reticulate](https://cran.r-project.org/package=reticulate)**, **[R6](https://cran.r-project.org/package=R6)** and **[processx](https://cran.r-project.org/package=processx)**. To install **rgee** run:

Stable version:

```{r eval = FALSE}
install.packages("rgee")
```

Dev version:

```{r eval = FALSE}
remotes::install_github("r-spatial/rgee")
```

**rgee** has two types of dependencies: <span
style="color:#b52b09">**strict dependencies**</span> that must be
present before **rgee** initialization (i.e. **[ee_Initialize()](https://r-spatial.github.io/rgee/reference/ee_Initialize.html)**) and the <span style="color:#857e04"><b>credentials dependencies</b></span> that 
unlock all **rgee** I/0 functionality with Google Drive (GD) and Google Cloud Storage (GCS).

If the strict dependencies are not installed, **rgee just will not work**. These dependencies are:

  * <span style="color:#b52b09"><b> Google account with Earth Engine
    activated </b></span>
    
  * <span style="color:#b52b09"><b> Python >= v3.5 </b></span>

  * <span style="color:#b52b09"><b> EarthEngine Python API (Python package) </b></span>

The activation of an **Earth Engine account** depends on each user, check
the official website of [Google Earth Engine](https://earthengine.google.com/) for more details. If you do not have a Python environment or a version of the EarthEngine Python API, we strongly recommend you run:

```{r eval = FALSE}
library(rgee)
ee_install(py_env = "rgee") # It is just necessary once!
```

This function will perform the following six tasks:

  1. If you do not have a Python environment, it will display an
  interactive  menu to install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
  (a free minimal installer for conda).
  
  2. Remove the previous Python environment defined with the same name if it exists.
  
  3. Create a new Python environment.

  4. Set the environmental variables EARTHENGINE_PYTHON and EARTHENGINE_ENV. These variables
  will be used to define the reticulate environmental variable [RETICULATE_PYTHON](https://rstudio.github.io/reticulate/articles/versions.html#providing-hints) when 
  rgee is loaded. 
  
  5. Install rgee Python dependencies: [Earth Engine Python API](https://pypi.org/project/earthengine-api/) and
  [numpy](https://pypi.org/project/numpy/).
  
  6. Ask to restart the R session in order to see changes.

**However, the use of rgee::ee_install() is not mandatory; you can instead use your own custom installation.** If you are an Rstudio v.1.4 > user, this [**tutorial**](https://github.com/r-spatial/rgee/tree/help/rstudio/) will help you to properly set a Python Environment with your R session without **rgee::ee_install()**.  Take into account that the Python Environment you set must have installed the  [Earth Engine Python API](https://pypi.org/project/earthengine-api/) and [Numpy](https://pypi.org/project/numpy/).

On the other hand, the <span style="color:#857e04"><b>credentials dependencies</b></span>
are only needed to move data from Google Drive and Google Cloud Storage to your local environment.
These dependencies are not mandatory. However, they will help you to create a seamless connection between
R and Earth Engine.  These dependencies are:

  -   <span style="color:#857e04">**Google Cloud Storage
      credential**</span>
      
  -   <span style="color:#857e04">**Google Drive credential**</span>

See the next section to learn how to correctly set both credentials.

## 7. Authentication

As we have seen previously, **rgee** deals with three different Google API's:

  -   Google Earth Engine
  
  -   Google Drive
  
  -   Google Cloud Storage

To authenticate/initialize either Google Drive or Google Cloud Storage, you just need to run:

```{r eval = FALSE}
library(rgee)
#ee_reattach() # reattach ee as a reserve word
# Initialize just Earth Engine
ee_Initialize() 
ee_Initialize(user = 'csaybar@gmail.com') # Use the argument email is not mandatory, but it's helpful to change of EE user.
# Initialize Earth Engine and GD
ee_Initialize(user = 'csaybar@gmail.com', drive = TRUE)
# Initialize Earth Engine and GCS
ee_Initialize(user = 'csaybar@gmail.com', gcs = TRUE)
# Initialize Earth Engine, GD and GCS
ee_Initialize(user = 'csaybar@gmail.com', drive = TRUE, gcs = TRUE)
```

If the Google account is verified and the permission is granted, you
will be directed to an authentication token. Copy this token and paste
it in your R console. Unlike Earth Engine and Google Drive, Google Cloud 
Storage needs to set up its credential manually ([link1](https://code.markedmondson.me/googleCloudStorageR/articles/googleCloudStorageR.html) and [link2](https://github.com/r-spatial/rgee/tree/help/gcs/)). In all cases, the user's credentials will be stored in: 

``` {r eval = FALSE}
ee_get_earthengine_path()
```

Remember you only have to authorize once, for subsequent sessions it will not be necessary.

## 8. Hello World

we know that installing **rgee** can be frustrating sometimes :( , so, congratulations if you've gotten this far :D :D. In this small example will show you how to display SRTM elevation values worldwide!

```{r eval = FALSE}
library(rgee)
ee_Initialize()
srtm <- ee$Image("USGS/SRTMGL1_003")
```

**Define visualization parameters**

```{r}
viz <- list(
  max = 4000,
  min = 0,
  palette = c("#000000","#5AAD5A","#A9AD84","#FFFFFF")
)
```

**Use Map$addLayer to visualize the map interactively**

```{r eval=FALSE}
Map$addLayer(
  eeObject = srtm,
  visParams =  viz,
  name = 'SRTM',
  legend = TRUE
)
```

<center>
<img src="https://user-images.githubusercontent.com/16768318/103684232-8bd56a80-4f8b-11eb-8bbb-fac3e1c2aada.png" width="90%">
</center>

## 9. Checking

The `ee_check()` function will help you for checking the sanity of your 
`rgee` installation and dependencies. 

-   `ee_check_python()` - Python version
-   `ee_check_credentials()` - Google Drive and GCS credentials
-   `ee_check_python_packages()` - Python packages

``` r
library(rgee)
ee_check()

ee_check_python()
ee_check_credentials()
ee_check_python_packages()
```
---
title: "Integrate Google Cloud Storage and rgee"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{5. Integrate Google Cloud Storage and rgee.}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---

<link rel="preconnect" href="https://fonts.gstatic.com"> <link href="https://fonts.googleapis.com/css2?family=Roboto&display=swap" rel="stylesheet">

```{=html}
<style type="text/css">
  body{
    font-size: 15pt;
    font-family: 'Roboto', sans-serif;
  }
  pre code{
    font-size: 12pt;
  }
  .list-group-item:last-child{
    font-size: 11pt;
    font-weight: bold;
  }
</style>
```
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)
```

*Google Cloud Platform may slightly change its website. If this occurs, help us by informing at  [rgee issues](https://github.com/r-spatial/rgee/issues/new). Based on [Interacting with Google Storage in R](https://bookdown.org/hegghammer/interacting_with_google_storage_in_r/interacting_with_google_storage.html) by [Thomas Hegghammer](https://github.com/Hegghammer).*
 

This tutorial explains how to integrate rgee and Google Cloud Storage (GCS) step by step. In rgee, GCS is used as an **intermediary container** for massive downloading/uploading of files which is more flexible than Google Drive. At today's date (December 2021), GCS is free for [uses of less than 0.5 GB](https://cloud.google.com/artifact-registry/pricing).

<img src=https://user-images.githubusercontent.com/16768318/146637221-f0c7d9c1-d779-47e1-8de2-14fd902df320.png>

## **Create a Service Account Key**

The bulk of GCS & rgee sync issues is related to the creation of a [service accounts key (SaK)](https://cloud.google.com/iam/docs/creating-managing-service-account-keys) with not enough privileges for
writing/reading in GCS buckets. In order to create, configure and locally store your SaK, perform as follow:


#### **1. Create a Google Cloud Platform (GCP) account**

Go to https://cloud.google.com/, and create an account. You will have to add your address and a credit card.

<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146672237-ff02d553-5d0c-4297-99ba-583522b6056f.png>
</center>
<br>

#### **2. Create a new project**

Go to console.

<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146672301-f117b742-f2ee-4057-8927-ecb6e99ef5c9.png>
</center>
<br>

If it is your first time using GCP, you are assigned a project named **_My first project_**. If you want to create a new project (**OPTIONAL**). First, click on 'My first project' in the top blue bar, just to the right of 'Google Cloud Platform', and then '**NEW PROJECT**'. 

<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146672770-4ace34e3-cae8-4157-8127-221f7248e4d2.png>
</center>
<br>

Create a project with the desired name.

<center>
<img src=https://user-images.githubusercontent.com/16768318/146675593-7646d367-e0dc-4eac-b84e-0362f6c57c5f.png>
</center>

#### **3: Activate GCS API**

By default, the GCS API is activated. Make sure of it by typing "Cloud Storage API" in the search browser.

<br>
<center>
<img src=https://user-images.githubusercontent.com/24624959/146641268-da151477-e993-4614-ad5c-e09129f673dd.png>
</center>
<br>

Click on 'Enable APIs and Services' (if the API is turned off). A green checkmark as the image below means
that the API is activated!

<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146673009-ec91a12b-ecba-4e4a-80a9-6925fa939903.png>
</center>
<br>


#### **4: Set up a service account**

Now we are going to create a [**service account**](https://cloud.google.com/iam/docs/service-accounts). A service account is used to 'sign in' applications that require one or many GCP services. In our case, **rgee** (the app) needs an account with **GCS admin privileges**. To create a service account, search for ‘Cloud Storage’ in the browser, and click on ‘Cloud Storage’ product.


<br>
<center>
<img src=https://user-images.githubusercontent.com/24624959/146641879-71809544-e171-4ea2-8509-0b58a602b119.png>
</center>
<br>


Click on 'settings'.

<br>
<center>
<img src=https://user-images.githubusercontent.com/24624959/146642118-51c496c9-f03e-4784-aaf9-6124da880123.png>
</center>
<br>


Click on 'INTEROPERABILITY', and then click on 'CREATE A KEY FOR A SERVICE ACCOUNT'.

<br>
<center>
<img src=https://user-images.githubusercontent.com/24624959/146642177-bb8329c5-4a08-4da3-af60-14f5fab23a4c.png>
</center>
<br>

Click on 'CREATE NEW ACCOUNT'.

<br>
<center>
<img src=https://user-images.githubusercontent.com/24624959/146642293-886c4e18-2f00-4e16-ac7c-9f607582d1e2.png>
</center>
<br>


Set a name (1), and define **Storage Admin** (**DON'T FORGET THIS STEP!!**) in role (3). 

<br>
<center>
<img src=https://user-images.githubusercontent.com/24624959/146642396-6b0b4e5f-4147-4229-b3f3-98830785d64d.png>
</center>
<br>

#### **5: Create and download a SaK as a json file.**
 
Once create the service account, we have to download the **S**ervice **a**ccount **K**ey to use GCS outside the Google Cloud Platform console. A SaK is just a JSON file with your public/private RSA keys. First, click the small edit icon on the bottom right (three horizontal lines). Then, go to API & Services and click on **credentials**.

<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146674751-b9733151-00e7-411c-8422-e500f8deaf19.png>
</center>

<br>

On the next page, click on the service account name.

<br>
<center>
<img src=https://user-images.githubusercontent.com/24624959/146643392-308cc7e8-127d-4566-86d2-a5554729edba.png>
</center>

click on 'KEYS', 'ADD KEY', and 'Create new key'.
<center>
<img src=https://user-images.githubusercontent.com/24624959/146643435-18ac2316-2b61-40a2-9f1a-c470df2ce071.png>
</center>

<br>
Then, select JSON format, and click 'create'.
<br>

<center>
<img src=https://user-images.githubusercontent.com/24624959/146643466-5bf5b7a7-6953-4520-9c23-ae453931a59d.png>
</center>
<br>

This should prompt a save file window. Save the file to your hard drive. You can change the name to something more memorable if you like (**but keep the “.json” extension**). Also, please take note of where you stored it. Now we are done in the Google Cloud Console and can finally start working in RStudio.

## **Copy the SaK in your system**

From rgee v.1.2.9000 we added [ee_utils_sak_copy](https://r-spatial.github.io/rgee/reference/ee_utils_sak_copy.html) and [ee_utils_sak_validate](https://r-spatial.github.io/rgee/reference/ee_utils_sak_validate.html) to help you to validate and store your SaK. Please run as follow to properly set your SaK in your system.

```{r}
# remotes::install_github("r-spatial/rgee") Install rgee v.1.3
library(rgee)

ee_Initialize("csaybar")

SaK_file <- "/home/csaybar/Downloads/SaK_rgee.json" # PUT HERE THE FULLNAME OF YOUR SAK.

# Assign the SaK to a EE user.
ee_utils_sak_copy(
  sakfile =  SaK_file,
  users = c("csaybar", "ryali93") # Unlike GD, we can use the same SaK for multiple users.
)

# Validate your SaK
ee_utils_sak_validate()
```

`ee_utils_sak_validate` evaluate if **rgee** with your SaK can: (1) create buckets, (2) write objects, (3) read objects, and (4) connect GEE and GCS. If it does not retrieve an error, `ee_Initialize(..., gcs = TRUE)` will work like a charm!. The next step is create your own **GCS bucket**. Consider that the bucket name you set must be **globally unique**. In other words, two buckets can not exist with the same name in Google Cloud Storage. 

```{r}
library(rgee)
library(jsonlite)
library(googleCloudStorageR)

ee_Initialize("csaybar", gcs = TRUE)

# Create your own container
project_id <- ee_get_earthengine_path() %>% 
  list.files(., "\\.json$", full.names = TRUE) %>% 
  jsonlite::read_json() %>% 
  '$'(project_id) # Get the Project ID

googleCloudStorageR::gcs_create_bucket("CHOOSE_A_BUCKET_NAME", projectId = project_id)

```

## **ERROR: Cannot insert legacy ACL for an object when uniform bucket-level access is enabled**

This is a [common issue](https://github.com/cloudyr/googleCloudStorageR/issues/111) related to the control access of buckets. GCS offers two systems for granting users permission to access your buckets and objects: [IAM](https://cloud.google.com/storage/docs/access-control/iam) (recommended, used in all Google Cloud services) and [Access Control Lists (ACL)](https://cloud.google.com/storage/docs/access-control/lists) (legacy access control system, only available in GCS).

If you use the Google Cloud Platform console to create a bucket, it will use **Uniform access** (**IAM is used alone to manage permissions**) by default.

<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146841722-6da1a14a-8203-4ad8-83bb-6903730a3efe.png>
</center>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146684785-a82ae364-bf77-48e8-8b36-974e73792453.png>
</center>
<br>


On the contrary if you use [`googleCloudStorageR::gcs_create_bucket`](https://code.markedmondson.me/googleCloudStorageR/reference/gcs_create_bucket.html), it will use fine-grained access (**IAM and ACLs to manage permissions**).


<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146684836-e38a3414-d6b6-441d-becc-b2697df8d4d1.png>
</center>
<br>


Why is this important?. It is important, because if you create a bucket using the first option an error message will arise: **"Cannot insert legacy ACL for an object when uniform bucket-level access is enabled. Read more at https://cloud.google.com/storage/docs/uniform-bucket-level-access"**.

```{r}
demo_data <- data.frame(a = 1:10, b = 1:10)

# Bad --------------------------------------------------
googleCloudStorageR::gcs_upload(
  file = demo_data,
  name = "demo_data.csv",
  bucket = "demo_0002" # Bucket with uniform control access
)
#  Error: Insert legacy ACL for an object when uniform bucket-level access
#  is enabled. Read more at https://cloud.google.com/storage/docs/uniform-bucket-level-access


# Good -------------------------------------------------
googleCloudStorageR::gcs_upload(
  file = demo_data,
  name = "demo_data.csv",
  bucket = "demo_0002", # Bucket with uniform control access
  predefinedAcl = "bucketLevel"
)
```

It happens due that `googleCloudStorageaR` by default expects buckets created with fine-grained access (ACL support, see [cloudyr/googleCloudStorageR#111](https://github.com/cloudyr/googleCloudStorageR/issues/111)). 
To avoid this issue, from rgee  v.1.2.9000 we opt to change the default [`predefinedAcl`](https://cran.r-project.org/web/packages/googleCloudStorageR/googleCloudStorageR.pdf) argument from 'private' to 'bucketLevel'. This simple change should avoid users dealing with access control issues. However, if for some reason a user needs to change the access control policy (maybe to reduce the data exposure), from rgee v.1.2.0 all the rgee GCS functions ([sf_as_ee](https://r-spatial.github.io/rgee/reference/sf_as_ee.html), [local_to_gcs](https://r-spatial.github.io/rgee/reference/local_to_gcs.html), [raster_as_ee](https://r-spatial.github.io/rgee/reference/raster_as_ee.html), and [stars_as_ee](https://r-spatial.github.io/rgee/reference/stars_as_ee.html)) support the `predefinedAcl` argument too (Thanks to @[jsocolar](https://github.com/jsocolar)).

## **ERROR in Earth Engine servers: Unable to write to bucket demo_0001 (permission denied).**

This error arises when GEE tries to send an exported task results but your **EE USER** does not have enough privileges to write/read the bucket. **Why does this occur if I have successfully configured my SaK in my local system?**. Well, the SaK ensures a smooth connection between your local environment and GCS, but **not between GEE and GCS**.

<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146687641-dc28511e-2ca3-4c45-a1d8-50ba01bbb64d.png>
</center>
<br>

For instance, imagine that you have access to 2 Google user accounts, one personal and one for your work (in this example, we will call them David and Cesar). Both accounts have access to GEE. David creates a SaK and sends it to Cesar. As a result of this action, both Cesar and David can work together in the same bucket, downloading and creating GCS objects (Local <-> GCS). If David tries to use `ee_as_raster(..., via='gcs')` (from GEE -> GCS -> Local), GEE will recognize that the owner of the bucket is David so it will allow the procedure (GEE -> GCS) and thanks to the SaK there will be no problem to send the information to his local environment (GCS -> Local). However, if Cesar tries to do the same, it will get an error when passing the information from GEE -> GCS because **GEE does not know that Cesar has David SaK in their local system**.

```{r}
library(rgee)

ee_Initialize(gcs = TRUE)


# Define an image.
img <- ee$Image("LANDSAT/LC08/C01/T1_SR/LC08_038029_20180810")$
  select(c("B4", "B3", "B2"))$
  divide(10000)

# Define an area of interest.
geometry <- ee$Geometry$Rectangle(
  coords = c(-110.8, 44.6, -110.6, 44.7),
  proj = "EPSG:4326",
  geodesic = FALSE
)

img_03 <- ee_as_raster(
  image = img,
  region = geometry,
  container = "demo_0001",
  via = "gcs",
  scale = 1000
)

# ERROR in Earth Engine servers: Unable to write to bucket demo_0001 (permission denied).
```

This error is quite easy to fix. Just go to the bucket in your Google Cloud Platform console.



<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146689541-4c3bf4af-ae0f-417b-b1e8-49950f617d36.png>
</center>

Click on 'PERMISSIONS', and click on 'ADD'.

<center>
<img src=https://user-images.githubusercontent.com/16768318/146689562-6b115a25-a6f5-4f13-840d-10f88df3424e.png>
</center>

Finally, add the Google user account of the **EE USER**. **Do not forget to add 'STORAGE ADMIN' in the role!**

<center>
<img src=https://user-images.githubusercontent.com/16768318/146689576-03b4af8f-e8a8-4506-9c77-b8e06b338ae3.png>
</center>
<br>

## **Conclusion**

Setting up a SaK for GCS can be quite frustrating but it is definitely worth it!. If you are still having problems setting up your SaK, feel free to clearly detail your problem in [rgee isues](https://github.com/r-spatial/rgee/issues).
---
title: "How to create and deploy rgee Shiny apps - I"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{4. rgee Apps}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---

<link rel="preconnect" href="https://fonts.gstatic.com"> <link href="https://fonts.googleapis.com/css2?family=Roboto&display=swap" rel="stylesheet">

```{=html}
<style type="text/css">
  body{
    font-size: 15pt;
    font-family: 'Roboto', sans-serif;
  }
  pre code{
    font-size: 12pt;
  }
  .list-group-item:last-child{
    font-size: 11pt;
    font-weight: bold;
  }
</style>
```
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)
```

Google Earth Engine (GEE) allows users to create apps by three different approaches: **user tokens**, [service accounts](https://developers.google.com/earth-engine/guides/service_account) and [client-side authentication](https://github.com/google/earthengine-api/tree/master/demos/client-auth). In this tutorial, we will get focus on the first option.

## **What is an EE user token (EEtk)?**

An EEtk is a 100-character text string (OAuth2 credential) stored on your local system that is used to identify and authenticate users (See Figure below). In other words, it will permit to connect the EE Web REST API with their local system.

<br>

<center>
<img src="vig_4.png" width="80%"/>
</center>

<br>

In **rgee** the authentication procedure is triggered internally by [**ee_Initialize**](https://r-spatial.github.io/rgee/reference/ee_Initialize.html). This [function](https://github.com/google/earthengine-api/blob/b4c068924b1c7574ee717761cb9fe0499a3b932b/python/ee/data.py#L210%60) will search for the '**_credentials_**' file  (it stores the EEtk) on the path: **\~/.config/earthengine/**. 

<br>

```{r}
library(rgee)
sprintf("%s/credentials", dirname(rgee::ee_get_earthengine_path()))
```

<br>

If the file exists, then an [Oauth2 Credential](**https://google-auth.readthedocs.io/en/stable/reference/google.oauth2.credentials.html**) object is created using a [refresh token](https://developers.google.com/identity/protocols/oauth2) grant. This *refresh token* must come from a Google account registered in GEE if not a *Bad Request error* will be invoked. Once the Oauth2 Credential is successfully loaded, it is dynamically passed to all EE methods (See [Initialize](https://github.com/google/earthengine-api/blob/b4c068924b1c7574ee717761cb9fe0499a3b932b/python/ee/__init__.py#L125)) in order to realize independent calls to the Web REST API (https://earthengine.googleapis.com). As you realize, the **_credentials_** file is crucial to interact with the EE API and if it does not exist on your system, it will simply not be possible to use rgee.


<center>
<br>
**Please never show or share your token with anyone, it will give full access to all your Earth Engine resources.**
</center>

## **Deploying a simple rgee shiny app on shinyapps.io**

Deploying a rgee application can be a bit tricky, as you must perform the following task:

1. Install R packages.
2. Install Python and third-party packages.
3. **Set the _credentials_ file in the path _~/.config/earthengine/_**.

The first step is automatically accomplished for [**shinyapps.io**](https://www.shinyapps.io/). On the other hand, the second and third steps need to configure manually the virtual machine. To make the process straightforward we create [shiny_rgee](https://github.com/csaybar/shiny_rgee_template) template.

<br>
<center>
<img src="https://user-images.githubusercontent.com/16768318/146476146-6196f937-ad22-4687-bff0-1ca85f86806d.png" width="90%"/>
</center>
<br>

To use shiny_rgee template, first download it by running on terminal:

```
git clone https://github.com/csaybar/shiny_rgee_template.git
```

Load the rgeeApp.Rproj and modify the .Renviron file according to their personal token user information. It is available in your shinyapps profile https://www.shinyapps.io/admin/#/tokens.

<br>
<center>
<img src="https://user-images.githubusercontent.com/16768318/146477123-6bec91ed-d5bb-42e6-aa6d-5f6b08d08d7f.png" width="90%"/>
</center>

<br>

[**.Renviron**](https://github.com/csaybar/shiny_rgee_template/blob/main/.Renviron)

```
    SHINY_ACC_NAME="your_account_name"
    TOKEN="a_token_you_got_from_shinyapps.io"
    SECRET="a_secret_you_recieved_fromshinyapps.io"
    MASTERNAME="name_of_the_shiny_app"
```
<br>

Finally run the [**deploy.R**](https://github.com/csaybar/shiny_rgee_template/blob/main/deploy.R) file.

```{r}
library(reticulate)
library(rsconnect)
library(rgee)

# 1. Create the credentials file
ee_Initialize()

# 2. Copy credentials file to the project folder
file_credentials <- sprintf("%s/credentials", dirname(rgee::ee_get_earthengine_path()))
file.copy(file_credentials, to = ".")

# 3. Set ShinyApps account info
# FIRST MODIFY LOCAL .Renviron!!
error_on_missing_name <- function(name){
  var <- Sys.getenv(name, unset=NA)
  if(is.na(var)){
    stop(paste0("cannot find ",name),call. = FALSE)
  }
  gsub("\"", '',var)
}

setAccountInfo(name   = error_on_missing_name("SHINY_ACC_NAME"),
               token  = error_on_missing_name("TOKEN"),
               secret = error_on_missing_name("SECRET"))


# 4. Run the application
deployApp(
  appFiles = c("app.R", "utils.R", "credentials"),
  appTitle = "rgee_app_demo",
  lint = FALSE
)

# 5. Delete EE credentials file
file.remove("credentials")
```

<br>

After a couple of minutes, the app will be available on shinyapps.io. 
See our live demo in **https://cesar-aybar.shinyapps.io/rgee_app_demo/**.
---
title: "**Best Practices**"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{3. Best Practices}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}    
---

<link rel="preconnect" href="https://fonts.gstatic.com">
<link href="https://fonts.googleapis.com/css2?family=Roboto&display=swap" rel="stylesheet">
<style type="text/css">
  body{
    font-size: 15pt;
    font-family: 'Roboto', sans-serif;
  }
  pre code{
    font-size: 12pt;
  }
  .list-group-item:last-child{
    font-size: 11pt;
    font-weight: bold;
  }
</style>


  ```{r setup, include=FALSE}
  knitr::opts_chunk$set(eval = FALSE)
  ```


**Adapted from [Google Earth Engine Documentation](https://developers.google.com/earth-engine/guides/best_practices).**

This doc describes coding practices that are intended to maximize the chance of success for complex or expensive Earth Engine computations. 

## **Avoid mixing client functions and objects with server functions and objects**

Earth Engine server objects are objects with constructors that start with `ee` (e.g. ee\$Image, ee\$Reducer) and any methods on such objects are server functions. Any object not constructed in this manner is a client object. Client objects may come from the R Earth Engine client (e.g. Map) or the R language (e.g. date, data.frame, c(), list()).

To avoid unintended behavior, do not mix client and server functions in your script as discussed [here](https://developers.google.com/earth-engine/guides/debugging). See [this page](https://developers.google.com/earth-engine/guides/client_server) for in-depth explanation of client vs. server in Earth Engine. The following example illustrates the dangers of mixing client and server functionality:


<img src="thumb_down.png" width=25px>
&nbsp; **Error** — This code doesn't work!



```{r}
# Won't work.
for (i in seq_len(table$size())) {
  print('No!') 
}
```


Can you spot the error? Note that `table$size()` is a server method on a server
object and can not be used with client-side functionality such as the `seq_len` function.

A situation in which you may want to use for-loops is with to display results 
with `Map`, since the Map object and methods are client-side.

<img src="thumb_up.png" width=25px>
&nbsp; **Good** — Use client functions for display Earth Engine spatial objects.


```{r}
l8_ts <- sprintf(
  "LANDSAT/LC08/C01/T1/LC08_044034_%s",
  c("20140318", "20140403","20140419","20140505")
)

display_l8ts <- list()
for (l8 in l8_ts) {
  ee_l8 <- ee$Image(l8)
  display_l8ts[[l8]] <- Map$addLayer(ee_l8)
}

Map$centerObject(ee_l8)
Reduce('+', display_l8ts)
```

Conversely, `map()` is a server function and client functionality won't work 
inside the function passed to map(). For example:

<img src="thumb_down.png" width=25px>
&nbsp; **Error** — This code doesn't work!

```{r}
table <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')

# Error:
foobar <- table$map(function(f) {
  print(f); # Can't use a client function here.
  # Can't Export, either.
})
```

<img src="thumb_up.png" width=25px>
&nbsp; **Good** — Use `map()` `set()`.

```{r}
table <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')

# Do something to every element of a collection.
withMoreProperties = table$map(function(f) {
  # Set a property.
  f$set("area_sq_meters", f$area())
})
print(withMoreProperties$first()$get("area_sq_meters")$getInfo())
```

You can also `filter()` the collection based on computed or existing properties 
and `print()` the result. Note that you can not print a collection with more 
5000 elements. If you get the "Collection query aborted after accumulating over
5000 elements" error, `filter()` or `limit()` the collection before printing.


## **Avoid converting to list unnecessarily**

Collections in Earth Engine are processed using optimizations that are broken by converting the collection to a `List` or `Array` type. Unless you need random access to collection elements (i.e. you need to get the i'th element of a collection), use filters on the collection to access individual collection elements. The following example illustrates the difference between type conversion (not recommended) and filtering (recommended) to access an element in a collection:

<img src="thumb_down.png" width=25px>
&nbsp; **Bad** — Don't convert to list unnecessarily!

```{r}
table <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017');

# Do NOT do this!!
list <- table$toList(table$size())
print(list$get(13)$getInfo()) # User memory limit exceeded.
```

Note that you can easily trigger errors by converting a collection to a list unnecessarily. The safer way is to use `filter()`:

<img src="thumb_up.png" width=25px>
&nbsp; **Good** — Use `filter()`.

```{r}
print(table$filter(ee$Filter$eq('country_na', 'Niger'))$first()$getInfo())
```

Note that you should [use filters as early as possible in your analysis](https://developers.google.com/earth-engine/guides/best_practices).

## **Avoid ee.Algorithms.If()**

Do not use `ee.Algorithms.If()` to implement branching logic, especially in a mapped function. As the following example illustrates, `ee.Algorithms.If()` can be memory intensive and is not recommended:

<img src="thumb_down.png" width=25px>
&nbsp; **Bad** — Don't use `If()`:

```{r}
table <- ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')

# Do NOT do this!
veryBad = table$map(function(f) {
  ee$Algorithms$If(
    condition = ee$String(f$get('country_na'))$compareTo('Chad')$gt(0),
    trueCase = f,      # Do something.
    falseCase = NULL   # Do something else.
  )
}, TRUE)
print(veryBad$getInfo()) # User memory limit exceeded.

# If() may evaluate both the true and false cases.
```


Note that the second argument to `map()` is `TRUE`. This means that the mapped 
function may return nulls and they will be dropped in the resultant collection. 
That can be useful (without `If()`), but here the easiest solution is to use a
filter:

<img src="thumb_up.png" width=25px>
&nbsp; **Good** — Use `filter()`.

```{r}
print(table$filter(ee$Filter$eq('country_na', 'Chad')))
```

As shown in [this tutorial](https://developers.google.com/earth-engine/tutorials/tutorial_js_03), a functional programming approach using filters is the correct way to apply one logic to some elements of a collection and another logic to the other elements of the collection.

## **Avoid reproject()**

Don't use `reproject` unless absolutely necessary. One reason you might want to use reproject() is to force `Map` display computations to happen at a specific scale so you can examine the results at your desired scale of analysis. In the next example, patches of hot pixels are computed and the count of pixels in each patch is computed. Run the example and click on one of the patches. Note that the count of pixels differs between the reprojected data the data that has not been reprojected.

```{r}
l8sr <- ee$ImageCollection("LANDSAT/LC08/C01/T1_SR")
sf <- ee$Geometry$Point(c(-122.405, 37.786))
Map$centerObject(sf, 13)

# A reason to reproject - counting pixels and exploring interactively.
image <- l8sr$filterBounds(sf)$
  filterDate("2019-06-01", "2019-12-31")$
  first()
Map$addLayer(image, list(bands = "B10", min = 2800, max = 3100), "image")

hotspots <- image$select("B10")$
  gt(3100)$
  selfMask()$
  rename("hotspots")

objectSize <- hotspots$connectedPixelCount(256)
# Beware of reproject!  Don't zoom out on reprojected data.
reprojected <- objectSize$reproject(hotspots$projection())
Map$addLayer(objectSize, list(min = 1, max = 256), "Size No Reproject", FALSE) +
Map$addLayer(reprojected, list(min = 1, max = 256), "Size Reproject", FALSE)
```

The reason for the discrepancy is because the [scale of analysis](https://developers.google.com/earth-engine/guides/scale) is set by the Code Editor zoom level. By calling `reproject()` you set the scale of the computation instead of the Map display. Use `reproject()` with extreme caution for reasons described in [this doc](https://developers.google.com/earth-engine/guides/projections).


## **Filter and select() first**

In general, filter input collections by time, location and/or metadata prior to doing anything else with the collection. Apply more selective filters before less selective filters. Spatial and/or temporal filters are often more selective. For example, note that `select()` and `filter()` are applied before `map()`:

```{r}
images <- ee$ImageCollection("COPERNICUS/S2_SR")
sf <- ee$Geometry$Point(c(-122.463, 37.768))

# Expensive function to reduce the neighborhood of an image.
reduceFunction <- function(image) {
  image$reduceNeighborhood(
    reducer = ee$Reducer$mean(),
    kernel = ee$Kernel$square(4)
  )
}

bands <- list("B4", "B3", "B2")
# Select and filter first!
reasonableComputation <- images$select(bands)$
  filterBounds(sf)$
  filterDate("2018-01-01", "2019-02-01")$
  filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", 1))$
  aside(ee_print)$ # Useful for debugging.
  map(reduceFunction)$
  reduce('mean')$
  rename(bands)

viz <- list(bands = bands, min = 0, max = 10000)
Map$addLayer(reasonableComputation, viz, "resonableComputation")
```

## **Use updateMask() instead of mask()**

The difference between `updateMask()` and `mask()` is that the former does a logical `and()` of the argument (the new mask) and the existing image mask whereas `mask()` simply replaces the image mask with the argument. The danger of the latter is that you can unmask pixels unintentionally. In this example, the goal is to mask pixels less than or equal to 300 meters elevation. As you can see (zoom out), using `mask()` causes a lot of pixels to become unmasked, pixels that don't belong in the image of interest:

```{r}
l8sr <- ee$ImageCollection("LANDSAT/LC08/C01/T1_SR")
sf <- ee$Geometry$Point(c(-122.40554461769182, 37.786807309873716))
aw3d30 <- ee$Image("JAXA/ALOS/AW3D30_V1_1")

Map$centerObject(sf, 7)

image <- l8sr$filterBounds(sf)$filterDate("2019-06-01", "2019-12-31")$first()
vis <- list(bands = c("B4", "B3", "B2"), min = 0, max = 3000)

Map$addLayer(image, vis, "image", FALSE)

mask <- aw3d30$select("AVE")$gt(300)
Map$addLayer(mask, {}, 'mask', FALSE)

# NO!  Don't do this!
badMask <- image$mask(mask)
Map$addLayer(badMask, vis, "badMask")

goodMask <- image.updateMask(mask)
Map$addLayer(goodMask, vis, "goodMask", FALSE)
```

## **Combine reducers**

If you need multiple statistics (e.g. mean and standard deviation) from a single input (e.g. an image region), it is more efficient to combine reducers. For example, to get image statistics, combine reducers as follows:

```{r}
image <- ee$Image('COPERNICUS/S2/20150821T111616_20160314T094808_T30UWU')

# Get mean and SD in every band by combining reducers.
stats <- image$reduceRegion(
  reducer = ee$Reducer$mean()$combine(
    reducer2 = ee$Reducer$stdDev(),
    sharedInputs = TRUE
  ),
  geometry = ee$Geometry$Rectangle(c(-2.15, 48.55, -1.83, 48.72)),
  scale = 10,
  bestEffort = TRUE # Use maxPixels if you care about scale.
)

print(stats$getInfo())

# Extract means and SDs to images.
meansImage <- stats$toImage()$select('.*_mean')
sdsImage <- stats$toImage()$select('.*_stdDev')
```

In this example, note that the mean reducer is combined with the standard deviation reducer and `sharedInputs` is true to enable a single pass through the input pixels. In the output dictionary, the name of the reducer is appended to the band name. To get mean and SD images (for example to normalize the input image), you can turn the values into an image and use regexes to extract means and SDs individually as demonstrated in the example.

## **Use Export**

For computations that result in "User memory limit exceeded" or "Computation timed out" errors in the Code Editor, the same computations may be able to succeed by using `Export`. This is because the timeouts are longer and the allowable memory footprint is larger when running in the batch system (where exports run). (There are other approaches you may want to try first as detailed in the debugging doc). Continuing the previous example, suppose that dictionary returned an error. You could obtain the results by doing something like:

```{r}
link <- '86836482971a35a5e735a17e93c23272'
task <- ee$batch$Export$table$toDrive(
  collection = ee$FeatureCollection(ee$Feature(NULL, stats)),
  description = paste0("exported_stats_demo_", link),
  fileFormat = "CSV"
)

# Using rgee I/O
task <- ee_table_to_drive(
  collection = ee$FeatureCollection(ee$Feature(NULL, stats)),
  description = paste0("exported_stats_demo_", link),
  fileFormat = "CSV"
)
task$start()
ee_monitoring(task)

exported_stats <- ee_drive_to_local(task = task,dsn = "exported_stats.csv")
read.csv(exported_stats)
```


Note that the link is embedded into the asset name, for reproducibility. Also note that if you want to export `toAsset`, you will need to supply a geometry, which can be anything, for example the image centroid, which is small and cheap to compute. (i.e. don't use a complex geometry if you don't need it).

See the debugging page for examples of using `Export` to resolve [Computation](https://developers.google.com/earth-engine/guides/debugging) timed out and [Too many concurrent aggregations](https://developers.google.com/earth-engine/guides/debugging). See [this doc](https://developers.google.com/earth-engine/guides/exporting) for details on exporting in general.

## **If you don't need to clip, don't use clip() **

Using `clip()` unnecessarily will increase computation time. Avoid `clip()` unless it's necessary to your analysis. If you're not sure, don't clip. An example of a bad use of clip:

<img src="thumb_down.png" width=25px>
&nbsp; **Bad** — Don't clip inputs unnecessarily!

```{r}
table <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')
l8sr <- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')

chad <- table$filter(ee$Filter$eq('country_na', 'Chad'))$first()

# Do NOT clip unless you need to.
unnecessaryClip <- l8sr$
  select('B4')$                           # Good.
  filterBounds(chad$geometry())$          # Good.
  filterDate('2019-01-01', '2019-12-31')$ # Good.
  map(function(image) {
    image$clip(chad$geometry())   # NO! Bad! Not necessary.
  })$
  median()$
  reduceRegion(
    reducer = ee$Reducer$mean(),
    geometry = chad$geometry(),
    scale = 30,
    maxPixels = 1e10
  )
print(unnecessaryClip$getInfo())
```

Clipping the input images can be skipped entirely, because the region is specified in the `reduceRegion()` call:

<img src="thumb_up.png" width=25px>
&nbsp; **Good** — Specify the region on the output.

```{r}
noClipNeeded <- l8sr$
  select('B4')$                          # Good.
  filterBounds(chad$geometry())$          # Good.
  filterDate('2019-01-01', '2019-12-31')$ # Good.
  median()$
  reduceRegion(
    reducer = ee$Reducer$mean(),
    geometry = chad$geometry(), # Geometry is specified here.
    scale = 30,
    maxPixels = 1e10
  )
print(noClipNeeded$getInfo())
```

If this computation times out, `Export` it as in [this example](https://developers.google.com/earth-engine/guides/best_practices).

## **If you need to clip with a complex collection, use clipToCollection()**

If you really need to clip something, and the geometries you want to use for clipping are in a collection, use `clipToCollection()`:

```{r}
ecoregions <- ee$FeatureCollection('RESOLVE/ECOREGIONS/2017')
image <- ee$Image('JAXA/ALOS/AW3D30_V1_1')

complexCollection <- ecoregions$
  filter(
    ee$Filter$eq(
      'BIOME_NAME',
      'Tropical & Subtropical Moist Broadleaf Forests'
    )
  )

Map$addLayer(complexCollection, {}, 'complexCollection')

clippedTheRightWay <- image$select('AVE')$
  clipToCollection(complexCollection)

Map$addLayer(clippedTheRightWay, {}, 'clippedTheRightWay', FALSE)
```

Do NOT use `featureCollection.geometry()` or `featureCollection.union()` on 
large and/or complex collections, which can be more memory intensive.

## **Don't use a complex collection as the region for a reducer**

If you need to do a spatial reduction such that the reducer pools inputs from multiple regions in a `FeatureCollection`, don't supply `featureCollection.geometry()` as the `geometry` input to the reducer. Instead, use `clipToCollection()` and a region large enough to encompass the bounds of the collection. For example:


```{r}
ecoregions <- ee$FeatureCollection('RESOLVE/ECOREGIONS/2017')
image <- ee$Image('JAXA/ALOS/AW3D30_V1_1')
complexCollection <- ecoregions$filter(
  ee$Filter$eq('BIOME_NAME', 'Tropical & Subtropical Moist Broadleaf Forests')
)

clippedTheRightWay <- image$select('AVE')$clipToCollection(complexCollection)
Map$addLayer(clippedTheRightWay, {}, 'clippedTheRightWay')

reduction <- clippedTheRightWay$reduceRegion(
  reducer = ee$Reducer$mean(),
  geometry = ee$Geometry$Rectangle(
    coords = c(-179.9, -50, 179.9, 50),  # Almost global.
    geodesic = FALSE
  ),
  scale = 30,
  maxPixels = 1e12
)

print(reduction$getInfo()) # If this times out, export it.
```

## **Use a non-zero errorMargin**

For possibly expensive geometry operations, use the largest error margin possible given the required precision of the computation. The error margin specifies the maximum allowable error (in meters) permitted during operations on geometries (e.g. during reprojection). Specifying a small error margin can result in the need to densify geometries (with coordinates), which can be memory intensive. It's good practice to specify as large an error margin as possible for your computation:

```{r}
ecoregions <- ee$FeatureCollection("RESOLVE/ECOREGIONS/2017")
complexCollection <- ecoregions$limit(10)

Map$centerObject(complexCollection)
Map$addLayer(complexCollection)

expensiveOps <- complexCollection$map(function(f) {
  f$buffer(10000, 200)$bounds(200)
})

Map$addLayer(expensiveOps, {}, 'expensiveOps')
```


## **Don't use a ridiculously small scale with reduceToVectors()**

If you want to convert a raster to a vector, use an appropriate scale. Specifying a very small scale can substantially increase computation cost. Set scale as high as possible give the required precision. For example, to get polygons representing global land masses:

```{r}
etopo <- ee$Image('NOAA/NGDC/ETOPO1')

# Approximate land boundary.
bounds <- etopo$select(0)$gt(-100)

# Non-geodesic polygon.
almostGlobal <- ee$Geometry$Polygon(
  coords = list(
    c(-180, -80),
    c(180, -80),
    c(180, 80),
    c(-180, 80),
    c(-180, -80)
  ),
  proj = "EPSG:4326",
  geodesic = FALSE
)

Map$addLayer(almostGlobal, {}, "almostGlobal")

vectors <- bounds$selfMask()$reduceToVectors(
  reducer = ee$Reducer$countEvery(),
  geometry = almostGlobal,
  # Set the scale to the maximum possible given
  # the required precision of the computation.
  scale = 50000
)

Map$addLayer(vectors, {}, "vectors")
```
In the previous example, note the use of a non-geodesic polygon for use in global reductions.

## **Don't use reduceToVectors() with reduceRegions()**

Don't use a `FeatureCollection` returned by `reduceToVectors()` as an input to `reduceRegions()`. Instead, add the bands you want to reduce before calling `reduceToVectors()`:

```{r}
etopo <- ee$Image('NOAA/NGDC/ETOPO1')
mod11a1 <- ee$ImageCollection('MODIS/006/MOD11A1')

# Approximate land boundary.
bounds <- etopo$select(0)$gt(-100)

# Non-geodesic polygon.
almostGlobal <- ee$Geometry$Polygon(
  coords = list(c(-180, -80), c(180, -80), c(180, 80), c(-180, 80), c(-180, -80)),
  proj = "EPSG:4326",
  geodesic = FALSE
)

lst <- mod11a1$first()$select(0)
means <- bounds$selfMask()$addBands(lst)$reduceToVectors(
  reducer = ee$Reducer$mean(),
  geometry = almostGlobal,
  scale = 1000,
  maxPixels = 1e10
)
print(means$limit(10)$getInfo())
```

Note that other ways of reducing pixels of one image within zones of another include [reduceConnectedCommponents()](https://developers.google.com/earth-engine/api_docs#ee.image.reduceconnectedcomponents/) and/or [grouping reducers](https://developers.google.com/earth-engine/api_docs#ee.image.reduceconnectedcomponents/).

## **Use fastDistanceTransform() for neighborhood operations**

For some convolution operations, `fastDistanceTransform()` may be more efficient than `reduceNeighborhood()` or `convolve()`. For example, to do erosion and/or dilation of binary inputs:

```{r}
aw3d30 <- ee$Image("JAXA/ALOS/AW3D30_V1_1")

# Make a simple binary layer from a threshold on elevation.
mask <- aw3d30$select("AVE")$gt(300)
Map$setCenter(-122.0703, 37.3872, 11)
Map$addLayer(mask, {}, "mask")

# Distance in pixel units.
distance <- mask$fastDistanceTransform()$sqrt()
# Threshold on distance (three pixels) for a dilation.
dilation <- distance$lt(3)
Map$addLayer(dilation, {}, "dilation")

# Do the reverse for an erosion.
notDistance <- mask$Not()$fastDistanceTransform()$sqrt()
erosion <- notDistance$gt(3)
Map$addLayer(erosion, {}, 'erosion')
```

## **Use the optimizations in reduceNeighborhood()**

If you need to perform a convolution and can't use `fastDistanceTransform()`, use the optimizations in `reduceNeighborhood()`.

```{r}
l8raw <- ee$ImageCollection('LANDSAT/LC08/C01/T1_RT')
composite <- ee$Algorithms$Landsat$simpleComposite(l8raw)

bands <- c('B4', 'B3', 'B2')

optimizedConvolution <- composite$select(bands)$reduceNeighborhood(
  reducer = ee$Reducer$mean(),
  kernel = ee$Kernel$square(3),
  optimization = "boxcar" # Suitable optimization for mean.
)$rename(bands)

viz <- list(bands = bands, min = 0, max = 72)
Map$setCenter(-122.0703, 37.3872, 11)
Map$addLayer(composite, viz, "composite") +
Map$addLayer(optimizedConvolution, viz, "optimizedConvolution")
```

## **Don't sample more data than you need**

Resist the urge to increase your training dataset size unnecessarily. Although increasing the amount of training data is an effective machine learning strategy in some circumstances, it can also increase computational cost with no corresponding increase in accuracy. (For an understanding of when to increase training dataset size, see [this reference](https://www.deeplearning.ai/programs/)). The following example demonstrates how requesting too much training data can result in the dreaded "Computed value is too large" error:

<img src="thumb_down.png" width=25px>
&nbsp; **Bad** — Don't sample too much data!

```{r}
l8raw <- ee$ImageCollection('LANDSAT/LC08/C01/T1_RT')
composite <- ee$Algorithms$Landsat$simpleComposite(l8raw)
labels <- ee$FeatureCollection('projects/google/demo_landcover_labels')

# No!  Not necessary.  Don't do this:
labels <- labels$map(function(f){f$buffer(100000, 1000)})

bands <- c('B2', 'B3', 'B4', 'B5', 'B6', 'B7')

training <- composite$select(bands)$sampleRegions(
  collection = labels,
  properties = list("landcover"),
  scale = 30
)

classifier <- ee$Classifier$smileCart()$train(
  features = training,
  classProperty = "landcover",
  inputProperties = bands
)

print(classifier$explain()) # Computed value is too large
```

The better approach is to start with a moderate amount of data and tune the hyperparameters of the classifier to determine if you can achieve your desired accuracy:

<img src="thumb_up.png" width=25px>
&nbsp; **Good** — Tune hyperparameters.

```{r}
l8raw <- ee$ImageCollection("LANDSAT/LC08/C01/T1_RT")
composite <- ee$Algorithms$Landsat$simpleComposite(l8raw)
labels <- ee$FeatureCollection("projects/google/demo_landcover_labels")

# Increase the data a little bit, possibly introducing noise.
labels <- labels$map(function(f) {f$buffer(100, 10)})

bands <- c('B2', 'B3', 'B4', 'B5', 'B6', 'B7')

data <- composite$select(bands)$sampleRegions(
  collection = labels,
  properties = list("landcover"),
  scale = 30
)

# Add a column of uniform random numbers called 'random'.
data <- data$randomColumn()

# Partition into training and testing.
training <- data$filter(ee$Filter$lt("random", 0.5))
testing <- data$filter(ee$Filter$gte("random", 0.5))

# Tune the minLeafPopulation parameter.
minLeafPops <- ee$List$sequence(1, 10)

accuracies <- minLeafPops$map(
  ee_utils_pyfunc(
    function(p) {
      classifier <- ee$Classifier$smileCart(minLeafPopulation = p)$
        train(
          features = training,
          classProperty = "landcover",
          inputProperties = bands
        )
      
      testing$
        classify(classifier)$
        errorMatrix("landcover", "classification")$
        accuracy()
    }
  )
)

minLeafPopulation_array <- accuracies$getInfo()
plot(
  x = minLeafPopulation_array,
  type = "b", 
  col = "blue",
  lwd = 2,
  ylab = "accuracy",
  xlim = c(0,10),
  xlab = "value",
  main = "Hyperparameter tunning (minLeafPopulation)"
)

```

In this example, the classifier is already very accurate, so there's not much tuning to do. You might want to choose the smallest tree possible (i.e. largest `minLeafPopulation`) that still has the required accuracy.

## **Export intermediate results**

Suppose your objective is to take samples from a relatively complex computed image. It is often more efficient to `Export` the image `toAsset()`, load the exported image, then sample. For example:

```{r}
image <- ee$Image('UMD/hansen/global_forest_change_2018_v1_6')
geometry <- ee$Geometry$Polygon(
  coords = list(
    c(-76.64069800085349, 5.511777325802095),
    c(-76.64069800085349, -20.483938229362376),
    c(-35.15632300085349, -20.483938229362376),
    c(-35.15632300085349, 5.511777325802095)
  ),
  proj =  "EPSG:4326",
  geodesic =  FALSE
)

testRegion <- ee$Geometry$Polygon(
  coords = list(
    c(-48.86726050085349, -3.0475996402515717),
    c(-48.86726050085349, -3.9248707849303295),
    c(-47.46101050085349, -3.9248707849303295),
    c(-47.46101050085349, -3.0475996402515717)
  ),
  proj = "EPSG:4326",
  geodesic = FALSE
)

# Forest loss in 2016, to stratify a sample.
loss <- image$select("lossyear")
loss16 <- loss$eq(16)$rename("loss16")

# Cloud masking function.
maskL8sr <- function(image) {
  cloudShadowBitMask <- bitwShiftL(1, 3)
  cloudsBitMask <- bitwShiftL(1, 5)
  qa <- image$select('pixel_qa')
  mask <- qa$bitwiseAnd(cloudShadowBitMask)$eq(0)$
    And(qa$bitwiseAnd(cloudsBitMask)$eq(0))
  
  image$updateMask(mask)$
    divide(10000)$
    select("B[0-9]*")$
    copyProperties(image, list("system:time_start"))
}

collection <- ee$ImageCollection("LANDSAT/LC08/C01/T1_SR")$map(maskL8sr)

# Create two annual cloud-free composites.
composite1 <- collection$filterDate('2015-01-01', '2015-12-31')$median()
composite2 <- collection$filterDate('2017-01-01', '2017-12-31')$median()

# We want a strtatified sample of this stack.
stack <- composite1$addBands(composite2)$float() # Export the smallest size possible.

# Export the image.  This block is commented because the export is complete.
# link <- "0b8023b0af6c1b0ac7b5be649b54db06"
# desc <- paste0(ee_get_assethome(), "/Logistic_regression_stack_", link)
# 
# #ee_image_info(stack)
# task <- ee_image_to_asset(
#   image = stack,
#   description = link,
#   assetId = desc,
#   region = geometry,
#   scale = 100,
#   maxPixels = 1e10
# )

  
# Load the exported image.
exportedStack <- ee$Image(
  "projects/google/Logistic_regression_stack_0b8023b0af6c1b0ac7b5be649b54db06"
)

# Take a very small sample first, to debug.
testSample <- exportedStack$addBands(loss16)$stratifiedSample(
  numPoints = 1,
  classBand = "loss16",
  region = testRegion,
  scale = 30,
  geometries = TRUE
)

print(testSample$getInfo()) # Check this in the console.

# Take a large sample.
sample <- exportedStack$addBands(loss16)$stratifiedSample(
  numPoints = 10000,
  classBand = "loss16",
  region = geometry,
  scale = 30
)

# Export the large sample...
```

In this example, note that the imagery is exported as float. Don't export at double precision unless absolutely necessary.

Once the export is completed, reload the asset and proceed with sampling from it. Note that a very small sample over a very small test area is run first, for debugging. When that is shown to succeed, take a larger sample and export it. Such large samples typically need to be exported. Do not expect such samples to be available interactively (for example through `print()`) or useable (for example as input to a classifier) without exporting them first.

## **Join vs. map-filter**

Suppose you want to join collections based on time, location or some metadata property. Generally, this is most efficiently accomplished with a join. The following example does a spatio-temporal join between the Landsat 8 and Sentinel-2 collections:

```{r}
s2 <- ee$ImageCollection("COPERNICUS/S2")$
  filterBounds(ee$Geometry$Point(c(-2.0205, 48.647)))

l8 <- ee$ImageCollection("LANDSAT/LC08/C01/T1_SR")

joined <- ee$Join$saveAll("landsat")$apply(
  primary = s2,
  secondary = l8,
  condition = ee$Filter$And(
    ee$Filter$maxDifference(
      difference = 1000 * 60 * 60 * 24, # One day in milliseconds
      leftField = "system:time_start",
      rightField = "system:time_start"
    ),
    ee$Filter$intersects(
      leftField = ".geo",
      rightField = ".geo"
    )
  )
)

print(joined$first()$getInfo())
```

Although you should try a join first (`Export` if needed), occasionally a `filter()` within a `map()` can also be effective, particularly for very large collections.


```{r}
s2 <- ee$ImageCollection("COPERNICUS/S2")$
  filterBounds(ee$Geometry$Point(c(-2.0205, 48.647)))

l8 <- ee$ImageCollection("LANDSAT/LC08/C01/T1_SR")

mappedFilter <- s2$map(function(image) {
  date <- image$date()
  landsat <- l8$
    filterBounds(image$geometry())$
    filterDate(date$advance(-1, "day"), date$advance(1, "day"))
    # Return the input image with matching scenes in a property.
  image$set(
    list(
      landsat = landsat,
      size = landsat$size()
    )
  )
})$filter(ee$Filter$gt("size", 0))

print(mappedFilter$first()$getInfo())
```


## **reduceRegion() vs. reduceRegions() vs. for-loop**

Calling `reduceRegions()` with a very large or complex
`FeatureCollection` as input may result in the dreaded “Computed value
is too large” error. One potential solution is to map `reduceRegion()`
over the `FeatureCollection` instead. Another potential solution is to
use a (gasp) for-loop. Although this is strongly discouraged in Earth
Engine as described [here](https://developers.google.com/earth-engine/guides/getstarted), [here](https://developers.google.com/earth-engine/tutorials/tutorial_js_03) and [here](https://developers.google.com/earth-engine/guides/client_server), `reduceRegion()` can be implemented in a for-loop to perform large reductions.

Suppose your objective is to obtain the mean of pixels (or any statistic) in each feature in a `FeatureCollection` for each image in an `ImageCollection`. The following example compares the three approaches previously described:

```{r}
# Table of countries.
countriesTable <- ee$FeatureCollection("USDOS/LSIB_SIMPLE/2017")

# Time series of images.
mod13a1 <- ee$ImageCollection("MODIS/006/MOD13A1")

# MODIS vegetation indices (always use the most recent version).
band <- "NDVI"
imagery <- mod13a1$select(band)

# Option 1: reduceRegions()
testTable <- countriesTable$limit(1) # Do this outside map()s and loops.

data <- imagery$map(function(image) {
  image$reduceRegions(
    collection = testTable,
    reducer = ee$Reducer$mean(),
    scale = 500
  )$map(function(f) {
    f$set(
      list(
        time = image$date()$millis(),
        date = image$date()$format()
      )
    )
  })
})$flatten()

print(data$first()$getInfo())

# Option 2: mapped reduceRegion()
data <- countriesTable$map(function(feature) {
  imagery$map(
    function(image) {
      ee$Feature(
        feature$geometry()$centroid(100),
        image$reduceRegion(
          reducer = ee$Reducer$mean(),
          geometry = feature$geometry(),
          scale = 500
        )
      )$set(
        list(
          time = image$date()$millis(),
          date = image$date()$format()
        )
      )$copyProperties(feature)
    }
  )
})$flatten()

print(data$first()$getInfo())

# Option 3: for-loop (WATCH OUT!)
size <- countriesTable$size()
print(size$getInfo()) # 312

countriesList <- countriesTable$toList(1) # Adjust size.
data <- ee$FeatureCollection(list()) # Empty table.

for (j in (seq_len(countriesList$length()$getInfo()) - 1)) {
  feature <- ee$Feature(countriesList$get(j))
  # Convert ImageCollection > FeatureCollection
  fc <- ee$FeatureCollection(
    imagery$map(
      function(image) {
        ee$Feature(
          feature$geometry()$centroid(100),
          image$reduceRegion(
            reducer = ee$Reducer$mean(),
            geometry = feature$geometry(),
            scale = 500
          )
        )$set(
          list(
            time = image$date()$millis(),
            date = image$date()$format()
          )
        )$copyProperties(feature)
      }
    )
  )
  data <- data$merge(fc)
}
print(data$first()$getInfo())
```

Note that the `first()` thing from each collection is printed, for debugging purposes. You should not expect that the complete result will be available interactively: you'll need to `Export`. Also note that for-loops should be used with extreme caution and only as a last resort. Finally, the for-loop requires manually obtaining the size of the input collection and hardcoding that in the appropriate locations. If any of that sounds unclear to you, don't use a for-loop.

## **Use forward differencing for neighbors in time**

Suppose you have a temporally sorted `ImageCollection` (i.e. a time series) and you want to compare each image to the previous (or next) image. Rather than use `iterate()` for this purpose, it may be more efficient to use an array-based forward differencing. The following example uses this method to de-duplicate the Sentinel-2 collection, where duplicates are defined as images with the same day of year:

```{r}
sentinel2 <- ee$ImageCollection("COPERNICUS/S2")
sf <- ee$Geometry$Point(c(-122.47555371521855, 37.76884708376152))
s2 <- sentinel2$
  filterBounds(sf)$
  filterDate("2018-01-01", "2019-12-31")

withDoys <- s2$map(function(image) {
  ndvi <- image$normalizedDifference(c("B4", "B8"))$rename("ndvi")
  date <- image$date()
  doy <- date$getRelative("day", "year")
  time <- image$metadata("system:time_start")
  doyImage <- ee$Image(doy)$
    rename("doy")$
    int()
  
  ndvi$
    addBands(doyImage)$
    addBands(time)$
    clip(image$geometry()) # Appropriate use of clip.
})

array <- withDoys$toArray()
timeAxis <- 0
bandAxis <- 1

dedup <- function(array) {
  time <- array$arraySlice(bandAxis, -1)
  sorted <- array$arraySort(time)
  doy <- sorted$arraySlice(bandAxis, -2, -1)
  left <- doy$arraySlice(timeAxis, 1)
  right <- doy$arraySlice(timeAxis, 0, -1)
  mask <- ee$Image(ee$Array(list(list(1))))$
    arrayCat(left$neq(right), timeAxis)
  array$arrayMask(mask)
}

deduped <- dedup(array)

# Inspect these outputs to confirm that duplicates have been removed.
print(array$reduceRegion("first", sf, 10)$getInfo())
print(deduped$reduceRegion("first", sf, 10)$getInfo())
```

---
title: "How to integrate Rmarkdown and rgee"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{5. Interactive maps display in Rmarkdown.}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---

<link rel="preconnect" href="https://fonts.gstatic.com"> <link href="https://fonts.googleapis.com/css2?family=Roboto&display=swap" rel="stylesheet">

```{=html}
<style type="text/css">
  body{
    font-size: 15pt;
    font-family: 'Roboto', sans-serif;
  }
  pre code{
    font-size: 12pt;
  }
  .list-group-item:last-child{
    font-size: 11pt;
    font-weight: bold;
  }
</style>
```
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)
```

## **1. The problem**

GEE offers on-the-fly computation for rendering EE spatial objects:

```{r}
library(rgee)
library(rgeeExtra)

ee_Initialize()

img <- ee$Image$Dataset$CGIAR_SRTM90_V4
Map$addLayer(log1p(img), list(min = 0, max = 7))
```

<br>
<center>
<img src=https://user-images.githubusercontent.com/16768318/146937907-17ba9a1c-c37d-4de1-9c61-5c7ee2d18aab.png>
</center>
<br>

However, this interactive map service **is temporary**, disappearing after a short period of time (~ 4 hours). This makes `Map$addLayer` unusable for report generation. In this vignette, we will learn to create a **permanent interactive map**.

## **2. A tentative workaround**

Instead of using GEE API for creating interactive maps, we will use [**titiler**](https://github.com/developmentseed/titiler). titiler creates web map tiles dynamically based on COG (STAC) resources. Since an exported EE task to retrieve images can return a COG, we just have to move these results to a **storage web service with [HTTP GET range requests](https://www.cogeo.org/)**.

<br>
<center>
<img src=https://user-images.githubusercontent.com/10407788/98711823-7f4de080-2353-11eb-9c8a-8a46550651ae.png>
</center>
<br>


Fortunately, [GCS counts with this feature](https://cloud.google.com/storage/docs/json_api/v1/objects/get), so if we manage to move our results to GCS, the work would be already done :)

```
GET /OBJECT_NAME HTTP/1.1
Host: BUCKET_NAME.storage.googleapis.com
Content-Length: 0
Authorization: AUTHENTICATION_STRING
Range: bytes=BYTE_RANGE
If-Match: ENTITY_TAG
If-Modified-Since: DATE
If-None-Match: ENTITY_TAG
If-Unmodified-Since: DATE
```

## 3. Show me the code!

First, load `rgee` and `googleCloudStorageR` and initialize the EE API. You must have correctly configured a service account key, if not check our tutorial "[**how to integrate Google Cloud Storage and rgee**](https://r-spatial.github.io/rgee/articles/rgee05.html)".

```{r}
library(rgee)
library(googleCloudStorageR)

# Init the EE API
ee_Initialize("csaybar", gcs = TRUE)

# Validate your SaK
# ee_utils_sak_validate(bucket = "rgee_examples")
```


Define your study area.

```{r}
# Define an study area
EE_geom <- ee$Geometry$Point(c(-70.06240, -6.52077))$buffer(5000)
```


Select an `ee$Image`, for instance, a Landsat-8 image.

```{r}
l8img <- ee$ImageCollection$Dataset$LANDSAT_LC08_C02_T2_L2 %>% 
  ee$ImageCollection$filterDate('2021-06-01', '2021-12-01') %>% 
  ee$ImageCollection$filterBounds(EE_geom) %>% 
  ee$ImageCollection$first()
```

Move `l8img` from EE to GCS.

```{r}
gcs_l8_name  <- "l8demo2" # name of the image in GCS.
BUCKET_NAME <- "rgee_examples" # set here your bucket name
task <- ee_image_to_gcs(
  image = l8img$select(sprintf("SR_B%s",1:5)),
  region = EE_geom,
  fileNamePrefix = gcs_l8_name,
  timePrefix = FALSE,
  bucket = BUCKET_NAME,
  scale = 10,
  formatOptions = list(cloudOptimized = TRUE) # Return a COG rather than a TIFF file.
)
task$start()
ee_monitoring()
```

Titiler needs resources downloadable for anyone. Therefore, **we recommend you to work with GCS buckets with fine-grained access**. In this way, you can decide individually which objects to make public. On the other hand, if you decide to work with buckets with uniform access, you will have to expose the entire bucket!. The code below makes a specific object in your bucket **public to internet**.


```{r}
# Make PUBLIC the GCS object 
googleCloudStorageR::gcs_update_object_acl(
  object_name = paste0(gcs_l8_name, ".tif"),
  bucket = BUCKET_NAME,
  entity_type = "allUsers"
)
```

Finally, use `Map$addLayer` to display the COG resource. By default, `Map$addLayer` use the open endpoint [`https://api.cogeo.xyz/docs`](`https://api.cogeo.xyz/docs`). 


```{r eval=TRUE, echo=FALSE}
library(rgee)
gcs_l8_name  <- "l8demo2" # name of the image in GCS.
BUCKET_NAME <- "rgee_examples" # set here your bucket name
```


```{r eval=TRUE}
img_id <- sprintf("https://storage.googleapis.com/%s/%s.tif", BUCKET_NAME, gcs_l8_name)
visParams <- list(bands=c("SR_B4","SR_B3","SR_B2"), min = 8000, max = 20000, nodata = 0)
Map$centerObject(img_id)
Map$addLayer(
  eeObject = img_id, 
  visParams = visParams,
  name = "My_first_COG",
  titiler_server = "https://api.cogeo.xyz/"
)
```

If you prefer to use [titiler syntax](https://api.cogeo.xyz/docs), set the parameter 
`titiler_viz_convert` as FALSE.


```{r eval=TRUE}
visParams <- list(expression = "B4,B3,B2", rescale = "8000, 20000", resampling_method = "cubic")
Map$addLayer(
  eeObject = img_id, 
  visParams = visParams,
  name = "My_first_COG",
  titiler_server = "https://api.cogeo.xyz/",
  titiler_viz_convert = FALSE
)
```

---
title: "**Considerations**"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 2
vignette: >
  %\VignetteIndexEntry{2. Considerations}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}    
---

<link rel="preconnect" href="https://fonts.gstatic.com">
<link href="https://fonts.googleapis.com/css2?family=Roboto&display=swap" rel="stylesheet">
<style type="text/css">
  body{
    font-size: 15pt;
    font-family: 'Roboto', sans-serif;
  }
  pre code{
    font-size: 12pt;
  }
  .list-group-item:last-child{
    font-size: 11pt;
    font-weight: bold;
  }
</style>

Thanks to **reticulate** is it possible to embed a Python session within an R session. The **[Earth Engine Python API](https://pypi.org/project/earthengine-api/)** and **rgee** share the **same modules, classes, functions, and methods**. In other words, the logic of the syntax is the same and just as fast (just change **.** by a **$**). Notwithstanding, differences in the language design of R and Python might cause some problems in specific scenarios. We identify **four** bug-potential cases. Each of them is explained in-depth below.


```{r, setup, include=FALSE}
knitr::opts_chunk$set(
  comment = '', fig.width = 6, fig.height = 6
)
```

### **1. The map message error:**

This issue happens when the **map** method is used while: (1) running a reticulate version
lower than &lt; 1.14 (please update it!); or (2) leading
with **ee$List** objects. For instance:

``` r
library(rgee)
ee$Initialize()
mylist = ee$List$sequence(10)
mylist$map(function(x) ee$Number(x)$add(1))
#> Error in py_call_impl(callable, dots$args, dots$keywords): RuntimeError: Evaluation error: argument "x" is missing, with no default.
#> 
#> Detailed traceback: 
#>   File "/home/aybarpc01/.virtualenvs/r-reticulate/lib/python3.7/site-packages/ee/apifunction.py", line 205, in <lambda>
#>     return lambda *args, **kwargs: func.call(*args, **kwargs)  # pylint: disable=unnecessary-lambda
#>   File "/home/aybarpc01/.virtualenvs/r-reticulate/lib/python3.7/site-packages/ee/function.py", line 67, in call
#>     return self.apply(self.nameArgs(args, kwargs))
#>   File "/home/aybarpc01/.virtualenvs/r-reticulate/lib/python3.7/site-packages/ee/function.py", line 80, in apply
#>     result = computedobject.ComputedObject(self, self.promoteArgs(named_args))
#>   File "/home/aybarpc01/.virtualenvs/r-reticulate/lib/python3.7/site-packages/ee/function.py", line 107, in promoteArgs
#>     promoted_args[name] = Function._promoter(args[name], spec['type'])
#>   File "/home/aybarpc01/.virtualenvs/r-reticulate/lib/python3.7/site-packages/ee/__init__.py", line 242, in _Promote
#>     return CustomFunction.create(arg, 'Object', ['Object'] * args_count)
#>   File "/home/aybarpc01/.virtualenvs/r-reticulate/lib/python3.7/site-packages/ee/customfunction.py", line 121, in create
#>     return CustomFunction(signature, func)
#>   File "/home/aybarpc01/.virtualenvs/r-reticulate/lib/python3.7/site-packages/ee/customfunction.py", line 47, in __init__
#>     self._body = body(*variables)
#>   File "/home/aybarpc01/R/x86_64-pc-linux-gnu-library/3.6/reticulate/python/rpytools/call.py", line 21, in python_function
#>     raise RuntimeError(res[kErrorKey])
```

The code before is perfectly valid but `rgee` will produce an error.
This problem should be easily solved by adding the function **ee_utils_pyfunc**.
It will permit to wrap R functions before to send it to `reticulate`. Let’s see:

``` r
library(rgee)
ee$Initialize()
mylist = ee$List$sequence(0,10)
mynewlist = mylist$map(
  ee_utils_pyfunc(
    function(x) ee$Number(x)$add(1)   
  )
)
mynewlist$getInfo()
#>  [1]  1  2  3  4  5  6  7  8  9 10 11
```

### **2. Do not forget the L**

By default, when you define a number in R it will produce a **double
precision** value. This does not happen in Python because, by default it
will create an **int** value.

**Python**

``` python
type(1)
#> <class 'int'>
```

**R**

``` r
class(1)
#> [1] "numeric"
```

But why does this matter? Let's explain with an example:

**Python**

``` python
ee.Initialize()
and_bitwise = ee.Number(32).bitwiseAnd(100)
and_bitwise.getInfo()
#> 32
```

**R**

``` r
and_bitwise = ee$Number(32)$bitwiseAnd(100) #caution: silent error
and_bitwise$getInfo()
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/aybarpc01/.local/lib/python3.7/site-packages/ee/computedobject.py", line 95, in getInfo
    return data.computeValue(self)
  File "/home/aybarpc01/.local/lib/python3.7/site-packages/ee/data.py", line 490, in computeValue
    return send_('/value', ({'json': obj.serialize(), 'json_format': 'v2'}))
  File "/home/aybarpc01/.local/lib/python3.7/site-packages/ee/data.py", line 1186, in send_
    raise ee_exception.EEException(json_content['error']['message'])
ee.ee_exception.EEException: Number.bitwiseAnd: Bitwise operands must be integer only.
```

Users need to take into consideration that most of the arguments of the
Earth Engine methods are strict to admit only **integer values**. The
creation of integers in R is quite simple; you just need to add the
letter **L** to the end of the specific number or employ the
function `as.integer`. The **correct code** in R would be:

``` r
and_bitwise = ee$Number(32L)$bitwiseAnd(100L)
and_bitwise$getInfo()
#> [1] 32
```

### **3. Be careful with ee$Date**

This problem also appears due to differences between the design of R and
Python as programming languages. Currently, R only supports integer data 
type of 32 bits. Such integers can only count up to about 2 billion. Unfortunately, 
this range is insufficient to deal with [Google Earth
Engine timestamp](https://developers.google.com/earth-engine/glossary/)
which is saved in milliseconds since the [UNIX epoch](https://en.wikipedia.org/wiki/Unix_time).

**Python**

``` python
my_date = ee.Date('1990-01-01')
my_date.getInfo()
#> {'type': 'Date', 'value': 631152000000} # greater than 2 billion
```

**R**

``` r
my_date <- ee$Date('1990-01-01')
my_date$getInfo()
#> $type
#> [1] "Date"
#> 
#> $value
#> [1] -208192512
```

The problems with `ee$Date` just appear in the last mile (Python to R or
vice-versa, `reticulate`), and they should not be too severe if treated
with care. `rgee` implements two functions to deal with Earth Engine
dates: `eedate_to_rdate` and `rdate_to_eedate`.

``` r
# Era5 dataset
era_img <- ee$ImageCollection("ECMWF/ERA5/DAILY")$
  filterDate("2019-01-01", "2019-12-31")$
  first()
# Extracting init date
ee_date <- era_img$get('system:time_start')
ee_date$getInfo() # Silent error
#> [1] 112573440
eedate_to_rdate(ee_date = ee_date, timestamp = TRUE)
#> [1] 1.546301e+12
```

### **4. Take into consideration reserved words in R**

A reserved word is a word that cannot be used as an identifier, such as the name
of a variable or a function. According with `?reserved`, the reserved words in R's parser
are: `if`, `else`, **`repeat`**, `while`, `function`, `for`, `in`, `next`, `break`, `TRUE`, `FALSE`, `NULL`,
`Inf`, `NaN`, `NA`, `NA_integer_`, `NA_real_`, `NA_complex_`, `NA_character_`. Of these words,
the only one that is part of the Earth Engine API is **repeat**.

We can find **repeat** as a
method for an Earth Engine List object. See **[`ee$List$repeat(value, count)`](https://developers.google.com/earth-engine/apidocs/ee-list-repeat)**:

``` r
library(rgee)
ee_Initialize()
ee_list <- ee$List(1:10)
ee_list$repeat(10,2)$getInfo()
#> Error: unexpected 'repeat' in "ee_list$repeat"
```

To avoid this error use backticks/quotation marks:

``` r
library(rgee)
ee_Initialize()
ee_list <- ee$List(1:10)
ee_list$'repeat'(10,2)$getInfo()
#> 10 10
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_Date.R
\name{eedate_to_rdate}
\alias{eedate_to_rdate}
\title{Pass an Earth Engine date object to R}
\usage{
eedate_to_rdate(ee_date, timestamp = FALSE)
}
\arguments{
\item{ee_date}{ee$date object (ee$Date)}

\item{timestamp}{Logical. If TRUE, return the date in milliseconds
from the Unix Epoch (1970-01-01 00:00:00 UTC). Otherwise, return the
date as a POSIXct object. By default FALSE.}
}
\value{
\code{eedate_to_rdate} will return either a numeric timestamp or
a POSIXct object depending on the \code{timestamp} argument.
}
\description{
Pass an Earth Engine date object to R
}
\details{
\code{eedate_to_rdate} is essential to avoid potential errors that
might appear when users call to retrieve dates. Currently,
R integer only supports 32 bit signed (such integers can only
count up to about 2 billion). This range is notably insufficient for dealing
with GEE date objects represented by timestamps in milliseconds since the
UNIX epoch. \code{eedate_to_rdate} uses Python in the backend to obtain the
date and convert it in float before exporting to R.
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()

eeDate <- ee$Date$fromYMD(2010,1,1)
eedate_to_rdate(eeDate,timestamp = TRUE) # good
eeDate$getInfo()$value # bad
}
}
\seealso{
Other date functions: 
\code{\link{ee_get_date_ic}()},
\code{\link{ee_get_date_img}()},
\code{\link{rdate_to_eedate}()}
}
\concept{date functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-download.R
\name{ee_monitoring}
\alias{ee_monitoring}
\title{Monitoring Earth Engine task progress}
\usage{
ee_monitoring(task, task_time = 5, eeTaskList = FALSE, quiet = FALSE)
}
\arguments{
\item{task}{List generated after a task is started (i.e., after run
\code{ee$batch$Task$start()}) or a character that represents the ID of a EE
task started.}

\item{task_time}{Numeric. How often (in seconds) should a task be polled?}

\item{eeTaskList}{Logical. If \code{TRUE}, all Earth Engine tasks will be
listed.}

\item{quiet}{Logical. Suppress info message}
}
\value{
An \code{ee$batch$Task} object with a state "COMPLETED" or "FAILED"
according to the Earth Engine server's response.
}
\description{
Monitoring Earth Engine task progress
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()
ee_monitoring(eeTaskList = TRUE)
}
}
\seealso{
Other helper functions: 
\code{\link{ee_help}()},
\code{\link{ee_print}()}
}
\concept{helper functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_utils_sak_copy}
\alias{ee_utils_sak_copy}
\title{Stores a Service account key (SaK) inside the EE folder}
\usage{
ee_utils_sak_copy(sakfile, users = NULL, delete = FALSE, quiet = FALSE)
}
\arguments{
\item{sakfile}{Character. SaK filename. If missing, the SaK of the first user is used.}

\item{users}{Character. The user related to the SaK file. A SaK
file can be related to multiple users.}

\item{delete}{Logical. If TRUE, the SaK filename is deleted after copy.}

\item{quiet}{Logical. Suppress info message}
}
\description{
Copy SaK in the ~/.config/earthengine/$USER.
}
\examples{
\dontrun{
library(rgee)

ee_Initialize()

sakfile <- "/home/rgee_dev/sak_file.json"
# Copy sakfile to the users 'csaybar' and 'ndef'
ee_utils_sak_copy(sakfile = sakfile, users = c("csaybar", "ndef"))

# Copy the sakfile of the user1 to the user2 and user3.
ee_utils_sak_copy(users = c("csaybar", "ndef", "ryali93"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6Map.R
\docType{data}
\name{Map}
\alias{Map}
\title{R6 object (Map) to display Earth Engine (EE) spatial objects}
\format{
An object of class environment with the
following functions:
\itemize{
\item  \strong{addLayer(eeObject, visParams, name = NULL, shown = TRUE,
opacity = 1, titiler_viz_convert = TRUE,
titiler_server = "https://api.cogeo.xyz/")}: Adds a given EE object to the
map as a layer. \cr
\itemize{
\item \strong{eeObject:} The object to add to the interactive map.\cr
\item \strong{visParams:} List of parameters for visualization.
See details.\cr
\item \strong{name:} The name of the layer.\cr
\item \strong{shown:} A flag indicating whether the
layer should be on by default. \cr
\item \strong{opacity:} The layer's opacity is represented as a number
between 0 and 1. Defaults to 1. \cr
\item \strong{titiler_viz_convert:} Logical. If it is TRUE, Map$addLayer
will transform the visParams to titiler style. Ignored if eeObject is
not a COG file. \cr
\item \strong{titiler_server:} TiTiler endpoint. Defaults to
"https://api.cogeo.xyz/".
}
\item  \strong{addLayers(eeObject, visParams, name = NULL, shown = TRUE,
opacity = 1)}: Adds a given ee$ImageCollection to the map
as multiple layers. \cr
\itemize{
\item \strong{eeObject:} The ee$ImageCollection to add to the interactive map.\cr
\item \strong{visParams:} List of parameters for visualization.
See details.\cr
\item \strong{name:} The name of layers.\cr
\item \strong{shown:} A flag indicating whether
layers should be on by default. \cr
\item \strong{opacity:} The layer's opacity is represented as a number
between 0 and 1. Defaults to 1. \cr
\item \strong{nmax:} Numeric. The maximum number of images to display.
By default 5.
}

\item  \strong{addLegend(visParams, name = "Legend", position = c("bottomright",
"topright", "bottomleft", "topleft"), color_mapping= "numeric", opacity = 1, ...)}:
Adds a given ee$ImageCollection to the map as multiple layers. \cr
\itemize{
\item \strong{visParams:} List of parameters for visualization.\cr
\item \strong{name:} The title of the legend.\cr
\item \strong{position:} Character. The position of the legend. By default bottomright. \cr
\item \strong{color_mapping:} Map data values (numeric or factor/character) to
colors according to a given palette. Use "numeric" ("discrete") for continuous
(categorical) data. For display characters use "character" and add to visParams
the element "values" containing the desired character names. \cr
\item \strong{opacity:} The legend's opacity is represented as a number between 0
and 1. Defaults to 1. \cr
\item \strong{...:} Extra legend creator arguments. See \link[leaflet]{addLegend}. \cr
}

\item \strong{setCenter(lon = 0, lat = 0, zoom = NULL)}: Centers the map
view at the given coordinates with the given zoom level. If no zoom level
is provided, it uses 1 by default.
\itemize{
\item \strong{lon:} The longitude of the center, in degrees.\cr
\item \strong{lat:} The latitude of the center, in degrees.\cr
\item \strong{zoom:} The zoom level, from 1 to 24.
}
\item \strong{setZoom(zoom = NULL)}: Sets the zoom level of the map.
\itemize{
\item \strong{zoom:} The zoom level, from 1 to 24.
}
\item \strong{centerObject(eeObject, zoom = NULL,
maxError = ee$ErrorMargin(1))}: Centers the
map view on a given object. If no zoom level is provided, it will
be predicted according to the bounds of the Earth Engine object specified.
\itemize{
\item \strong{eeObject:} EE object.\cr
\item \strong{zoom:} The zoom level, from 1 to 24.
\item \strong{maxError:} 	Max error when input
image must be reprojected to an explicitly
requested result projection or geodesic state.
}
}
}
\usage{
Map
}
\value{
Object of class leaflet, with the following extra parameters: tokens, name,
opacity, shown, min, max, palette, and legend. Use the $ method to retrieve
the data (e.g. m$rgee$min).
}
\description{
Create interactive visualizations of spatial EE objects
(ee$FeatureCollection, ee$ImageCollection, ee$Geometry, ee$Feature, and
ee$Image.) using \code{leaflet} in the backend.
}
\details{
\code{Map} use the Earth Engine method
\href{https://developers.google.com/earth-engine/api_docs#ee.data.getmapid/}{
getMapId} to fetch and return an ID dictionary being used to create
layers in a \code{leaflet} object. Users can specify visualization
parameters to Map$addLayer by using the visParams argument. Each Earth
Engine spatial object has a specific format. For
\code{ee$Image}, the
\href{https://developers.google.com/earth-engine/guides/image_visualization}{
parameters} available are:

\tabular{lll}{
\strong{Parameter}\tab \strong{Description}  \tab \strong{Type}\cr
\strong{bands}    \tab  Comma-delimited list of three band (RGB) \tab  list \cr
\strong{min}      \tab  Value(s) to map to 0 \tab  number or list of three
numbers, one for each band \cr
\strong{max}      \tab  Value(s) to map to 1 \tab  number or list of three
numbers, one for each band \cr
\strong{gain}     \tab  Value(s) by which to multiply each pixel value \tab
number or list of three numbers, one for each band \cr
\strong{bias}     \tab  Value(s) to add to each Digital Number
value \tab number or list of three numbers, one for each band \cr
\strong{gamma}    \tab  Gamma correction factor(s) \tab  number or list of
three numbers, one for each band \cr
\strong{palette}  \tab  List of CSS-style color strings
(single-band only) \tab  comma-separated list of hex strings \cr
\strong{opacity}   \tab  The opacity of the layer (from 0 to 1)  \tab  number \cr
}

If you add an \code{ee$Image} to Map$addLayer without any additional
parameters, by default it assigns the first three bands to red,
green, and blue bands, respectively. The default stretch is based on the
min-max range. On the other hand, the available parameters for
\code{ee$Geometry}, \code{ee$Feature}, and \code{ee$FeatureCollection}
are:

\itemize{
\item \strong{color}: A hex string in the format RRGGBB specifying the
color to use for drawing the features. By default #000000.
\item \strong{pointRadius}: The radius of the point markers. By default 3.
\item \strong{strokeWidth}: The width of lines and polygon borders. By
default 3.
}
}
\examples{
\dontrun{
library(rgeeExtra)
library(rgee)
library(sf)

ee_Initialize()

# Case 1: Geometry*
geom1 <- ee$Geometry$Point(list(-73.53, -15.75))
Map$centerObject(geom1, zoom = 8)
m1 <- Map$addLayer(
  eeObject = geom1,
  visParams = list(
    pointRadius = 10,
    color = "FF0000"
  ),
  name = "Geometry-Arequipa"
)

# Case 2: Feature
feature_arq <- ee$Feature(ee$Geometry$Point(list(-72.53, -15.75)))
m2 <- Map$addLayer(
  eeObject = feature_arq,
  name = "Feature-Arequipa"
)
m2 + m1

# Case 4: Image
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")
Map$centerObject(image)
m4 <- Map$addLayer(
  eeObject = image,
  visParams = list(
    bands = c("B4", "B3", "B2"),
    max = 10000
  ),
  name = "SF"
)

# Case 5: ImageCollection
nc <- st_read(system.file("shape/nc.shp", package = "sf")) \%>\%
  st_transform(4326) \%>\%
  sf_as_ee()

ee_s2 <- ee$ImageCollection("COPERNICUS/S2")$
  filterDate("2016-01-01", "2016-01-31")$
  filterBounds(nc) \%>\%
  ee_get(0:4)
Map$centerObject(nc$geometry())
m5 <- Map$addLayers(ee_s2)
m5

# Case 6: Map comparison
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")
Map$centerObject(image)
m_ndvi <- Map$addLayer(
  eeObject = image$normalizedDifference(list("B5", "B4")),
  visParams = list(min = 0, max = 0.7),
  name = "SF_NDVI"
) + Map$addLegend(list(min = 0, max = 0.7), name = "NDVI", position = "bottomright", bins = 4)
m6 <- m4 | m_ndvi
m6

# Case 7: digging up the metadata
m6$rgee$tokens
m5$rgee$tokens

# Case 8: COG support
# See parameters here: https://api.cogeo.xyz/docs

server <- "https://storage.googleapis.com/pdd-stac/disasters/"
file <- "hurricane-harvey/0831/20170831_172754_101c_3B_AnalyticMS.tif"
resource <- paste0(server, file)
visParams <- list(bands = c("B3", "B2", "B1"), min = 3000, max = 13500, nodata = 0)
Map$centerObject(resource)
Map$addLayer(resource, visParams = visParams, shown = TRUE)
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_utils_shp_to_zip}
\alias{ee_utils_shp_to_zip}
\title{Create a zip file from an sf object}
\usage{
ee_utils_shp_to_zip(
  x,
  filename,
  SHP_EXTENSIONS = c("dbf", "prj", "shp", "shx")
)
}
\arguments{
\item{x}{sf object}

\item{filename}{data source name}

\item{SHP_EXTENSIONS}{file extension of the files to save
into the zip file. By default: "dbf", "prj", "shp", "shx".}
}
\value{
Character. The full path of the created zip file.
}
\description{
Create a zip file from an sf object
}
\examples{
\dontrun{
library(rgee)
library(sf)
ee_Initialize(gcs = TRUE)

# Create sf object
nc <- st_read(system.file("shape/nc.shp", package="sf"))
zipfile <- ee_utils_shp_to_zip(nc)
}
}
\seealso{
Other ee_utils functions: 
\code{\link{ee_utils_py_to_r}()},
\code{\link{ee_utils_pyfunc}()}
}
\concept{ee_utils functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_download.R
\name{ee_table_to_asset}
\alias{ee_table_to_asset}
\title{Creates a task to export a FeatureCollection to an EE table asset.}
\usage{
ee_table_to_asset(
  collection,
  description = "myExportTableTask",
  assetId = NULL,
  overwrite = FALSE
)
}
\arguments{
\item{collection}{The feature collection to be exported.}

\item{description}{Human-readable name of the task.}

\item{assetId}{The destination asset ID. **kwargs: Holds other
keyword arguments that may have been deprecated.}

\item{overwrite}{Logical. If TRUE, the assetId will be overwritten if
it exists.}
}
\value{
An unstarted Task that exports the table to Earth Engine Asset.
}
\description{
Creates a task to export a FeatureCollection to an EE table asset.
This function is a wrapper around \code{ee$batch$Export$table$toAsset(...)}.
}
\examples{
\dontrun{
library(rgee)
library(stars)
library(sf)

ee_users()
ee_Initialize()

# Define study area (local -> earth engine)
# Communal Reserve Amarakaeri - Peru
rlist <- list(xmin = -71.13, xmax = -70.95,ymin = -12.89, ymax = -12.73)
ROI <- c(rlist$xmin, rlist$ymin,
         rlist$xmax, rlist$ymin,
         rlist$xmax, rlist$ymax,
         rlist$xmin, rlist$ymax,
         rlist$xmin, rlist$ymin)
ee_ROI <- matrix(ROI, ncol = 2, byrow = TRUE) \%>\%
  list() \%>\%
  st_polygon() \%>\%
  st_sfc() \%>\%
  st_set_crs(4326) \%>\%
  sf_as_ee()

amk_fc <- ee$FeatureCollection(
  list(ee$Feature(ee_ROI, list(name = "Amarakaeri")))
)

assetid <- paste0(ee_get_assethome(), '/geom_Amarakaeri')
task_vector <- ee_table_to_asset(
  collection = amk_fc,
  overwrite = TRUE,
  assetId = assetid
)
task_vector$start()
ee_monitoring(task_vector) # optional

ee_fc <- ee$FeatureCollection(assetid)
Map$centerObject(ee_fc)
Map$addLayer(ee_fc)
}
}
\seealso{
Other vector export task creator: 
\code{\link{ee_table_to_drive}()},
\code{\link{ee_table_to_gcs}()}
}
\concept{vector export task creator}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_get.R
\name{ee_get_assethome}
\alias{ee_get_assethome}
\title{Get the Asset home name}
\usage{
ee_get_assethome()
}
\value{
Character. The name of the Earth Engine Asset home
(e.g. users/datacolecfbf)
}
\description{
Get the Asset home name
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()
ee_get_assethome()
}
}
\seealso{
Other path utils: 
\code{\link{ee_get_earthengine_path}()}
}
\concept{path utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_help.R
\name{ee_help}
\alias{ee_help}
\title{Documentation for Earth Engine Objects}
\usage{
ee_help(eeobject, browser = FALSE)
}
\arguments{
\item{eeobject}{Earth Engine Object to print documentation.}

\item{browser}{Logical. Display documentation in the browser.}
}
\value{
No return value, called for displaying Earth Engine documentation.
}
\description{
Documentation for Earth Engine Objects
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()

ee$Image()$geometry()$centroid \%>\% ee_help()
ee$Image()$geometry() \%>\% ee_help()
ee$Geometry$Rectangle(c(-110.8, 44.6, -110.6, 44.7)) \%>\% ee_help()
ee$Image \%>\% ee_help()
ee$Image \%>\% ee_help(browser = TRUE)
}
}
\seealso{
Other helper functions: 
\code{\link{ee_monitoring}()},
\code{\link{ee_print}()}
}
\concept{helper functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sf_as_ee.R
\name{sf_as_ee}
\alias{sf_as_ee}
\title{Convert an sf object to an EE object}
\usage{
sf_as_ee(
  x,
  via = "getInfo",
  assetId = NULL,
  bucket = NULL,
  predefinedAcl = "bucketLevel",
  command_line_tool_path = NULL,
  overwrite = TRUE,
  monitoring = TRUE,
  proj = "EPSG:4326",
  evenOdd = TRUE,
  geodesic = NULL,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{x}{object of class sf, sfc or sfg.}

\item{via}{Character. Upload method for sf objects. Three methods are
implemented: 'getInfo', 'getInfo_to_asset' and 'gcs_to_asset'. See details.}

\item{assetId}{Character. Destination asset ID for the uploaded file. Ignore
if \code{via} argument is "getInfo".}

\item{bucket}{Character. Name of the bucket (GCS) to save intermediate files
(ignore if \code{via} is not defined as "gcs_to_asset").}

\item{predefinedAcl}{Specify user access to object. Passed to
\code{googleCloudStorageR::gcs_upload}.}

\item{command_line_tool_path}{Character. Path to the Earth Engine command line
tool (CLT). If NULL, rgee assumes that CLT is set in the system PATH.
(ignore if \code{via} is not defined as "gcs_to_asset").}

\item{overwrite}{A boolean argument that indicates indicating
whether "filename" should be overwritten. Ignore if \code{via} argument
is "getInfo". By default TRUE.}

\item{monitoring}{Logical. Ignore if via is not set as
\code{getInfo_to_asset} or \code{gcs_to_asset}. If TRUE the exportation task
will be monitored.}

\item{proj}{Integer or character. Coordinate Reference System (CRS) for the
EE object, defaults to "EPSG:4326" (x=longitude, y=latitude).}

\item{evenOdd}{Logical. Ignored if \code{x} is not a Polygon. If TRUE,
polygon interiors will be determined by the even/odd rule, where a point
is inside if it crosses an odd number of edges to reach a point at infinity.
Otherwise polygons use the left-inside rule, where interiors are on the
left side of the shell's edges when walking the vertices in the given order.
If unspecified, defaults to TRUE.}

\item{geodesic}{Logical. Ignored if \code{x} is not a Polygon or LineString.
Whether line segments should be interpreted as spherical geodesics. If
FALSE, indicates that line segments should be interpreted as planar lines
in the specified CRS. If absent, defaults to TRUE if the CRS is geographic
(including the default EPSG:4326), or to FALSE if the CRS is projected.}

\item{quiet}{Logical. Suppress info message.}

\item{...}{\code{\link{ee_utils_create_manifest_table}} arguments might be included.}
}
\value{
When \code{via} is "getInfo" and \code{x} is either an sf or sfc object
with multiple geometries will return an \code{ee$FeatureCollection}. For
single sfc and sfg objects will return an \code{ee$Geometry$...}.

If \code{via} is either "getInfo_to_asset" or "gcs_to_asset" always will
return an \code{ee$FeatureCollection}.
}
\description{
Load an sf object to Earth Engine.
}
\details{
\code{sf_as_ee} supports the upload of \code{sf} objects by three different
options: "getInfo" (default), "getInfo_to_asset", and "gcs_to_asset". \code{getInfo}
transforms sf objects (sfg, sfc, or sf) to GeoJSON (using \code{geojsonio::geojson_json})
and then encrusted them in an HTTP request using the server-side objects that are
implemented in the Earth Engine API (i.e. ee$Geometry$...). If the sf object is too
large (~ >1Mb) is likely to cause bottlenecks since it is a temporary
file that is not saved in your EE Assets (server-side). The second option implemented
is 'getInfo_to_asset'. It is similar to the previous one, with the difference
that after create the server-side object will save it in your Earth Engine
Assets. For dealing with very large spatial objects is preferable to use the
third option 'gcs_to_asset'. This option firstly saves the sf object as a *.shp
file in the /temp directory. Secondly, using the function \code{local_to_gcs}
will move the shapefile from local to Google Cloud Storage. Finally, using
the function \code{gcs_to_ee_table} the ESRI shapefile will be loaded
to their EE Assets. See \href{https://developers.google.com/earth-engine/guides/table_upload/}{Importing
table data} documentation for more details.
}
\examples{
\dontrun{
library(rgee)
library(sf)
ee_Initialize()

# 1. Handling geometry parameters
# Simple
ee_x <- st_read(system.file("shape/nc.shp", package = "sf")) \%>\%
  sf_as_ee()

Map$centerObject(eeObject = ee_x)
Map$addLayer(ee_x)

# Create a right-inside polygon.
toy_poly <- matrix(data = c(-35,-10,-35,10,35,10,35,-10,-35,-10),
                   ncol = 2,
                   byrow = TRUE) \%>\%
  list() \%>\%
  st_polygon()

holePoly <- sf_as_ee(x = toy_poly, evenOdd = FALSE)

# Create an even-odd version of the polygon.
evenOddPoly <- sf_as_ee(toy_poly, evenOdd = TRUE)

# Create a point to test the insideness of the polygon.
pt <- ee$Geometry$Point(c(1.5, 1.5))

# Check insideness with a contains operator.
print(holePoly$contains(pt)$getInfo() \%>\% ee_utils_py_to_r())
print(evenOddPoly$contains(pt)$getInfo() \%>\% ee_utils_py_to_r())

# 2. Upload small geometries to EE asset
assetId <- sprintf("\%s/\%s", ee_get_assethome(), 'toy_poly')
eex <- sf_as_ee(
 x = toy_poly,
 overwrite = TRUE,
 assetId = assetId,
 via = "getInfo_to_asset")
# 3. Upload large geometries to EE asset
ee_Initialize(gcs = TRUE)
assetId <- sprintf("\%s/\%s", ee_get_assethome(), 'toy_poly_gcs')
eex <- sf_as_ee(
  x = toy_poly,
  overwrite = TRUE,
  assetId = assetId,
  bucket = 'rgee_dev',
  monitoring = FALSE,
  via = 'gcs_to_asset'
)
ee_monitoring()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_Initialize.R
\name{ee_users}
\alias{ee_users}
\title{Display the credentials of all users as a table}
\usage{
ee_users(quiet = FALSE)
}
\arguments{
\item{quiet}{Logical. Suppress info messages.}
}
\value{
A data.frame with credential information of all users.
}
\description{
Display Earth Engine, Google Drive, and Google Cloud Storage Credentials as
a table.
}
\examples{
\dontrun{
library(rgee)
ee_users()
}
}
\seealso{
Other session management functions: 
\code{\link{ee_Initialize}()},
\code{\link{ee_user_info}()},
\code{\link{ee_version}()}
}
\concept{session management functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_get.R
\name{ee_get_date_ic}
\alias{ee_get_date_ic}
\title{Get the date of a EE ImageCollection}
\usage{
ee_get_date_ic(x, time_end = FALSE)
}
\arguments{
\item{x}{ee$ImageCollection object}

\item{time_end}{Logical. If TRUE, the \code{system:time_end} property is
also returned. See details.}
}
\value{
A data.frame with the columns: \code{id} (ID of the image),
\code{time_start}, and \code{time_end} (only if the argument \code{time_end} is
set as TRUE). The number of rows depends on the number of images
(\code{ee$ImageCollection$size}).
}
\description{
Get the date of a EE ImageCollection
}
\details{
\code{system:time_start} Sets the start period of data acquisition while
\code{system:time_end} does the same for the end period. See the
\href{https://developers.google.com/earth-engine/glossary/}{Earth Engine glossary}
for getting more information.
}
\examples{
\dontrun{
library(rgee)
library(sf)
ee_Initialize()

nc <- st_read(system.file("shape/nc.shp", package = "sf")) \%>\%
  st_transform(4326) \%>\%
  sf_as_ee()

ee_s2 <- ee$ImageCollection("COPERNICUS/S2")$
  filterDate("2016-01-01", "2016-01-31")$
  filterBounds(nc)

ee_get_date_ic(ee_s2)
}
}
\seealso{
Other date functions: 
\code{\link{ee_get_date_img}()},
\code{\link{eedate_to_rdate}()},
\code{\link{rdate_to_eedate}()}
}
\concept{date functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_download.R
\name{ee_gcs_to_local}
\alias{ee_gcs_to_local}
\title{Move results from Google Cloud Storage to a local directory}
\usage{
ee_gcs_to_local(
  task,
  dsn,
  public = FALSE,
  metadata = FALSE,
  overwrite = TRUE,
  quiet = FALSE
)
}
\arguments{
\item{task}{List generated after finished an EE task correctly. See details.}

\item{dsn}{Character. Output filename. If missing, a temporary
file (i.e. \code{tempfile()}) is assigned.}

\item{public}{Logical. If TRUE, a public link to Google Cloud Storage
resource is created.}

\item{metadata}{Logical. If TRUE, export the metadata related to the Google
Cloud Storage resource. See details.}

\item{overwrite}{A boolean argument that indicates indicating
whether "filename" should be overwritten. By default TRUE.}

\item{quiet}{Logical. Suppress info message}
}
\value{
If \code{metadata} is FALSE, will return the filename of the Google
Cloud Storage resource on their system. Otherwise, a list with two elements
(\code{dns} and \code{metadata}) is returned.
}
\description{
Move results of an EE task saved in Google Cloud Storage to a local
directory.
}
\details{
The task argument needs "COMPLETED" task state to work due to that the parameters
necessaries to locate the file into Google Cloud Storage are obtained from \cr
\code{ee$batch$Export$*$toCloudStorage(...)$start()$status()}.

If the argument \code{metadata} is TRUE, a list with the
following elements is exported join with the output filename (dsn):

\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{gcs_name: }}{Name of the Table in Google Cloud Storage.}
\item{\bold{gcs_bucket: }}{Name of the bucket.}
\item{\bold{gcs_fileFormat: }}{Format of the table.}
\item{\bold{gcs_public_link: }}{Download link to the table.}
\item{\bold{gcs_URI: }}{gs:// link to the table.}
}
}
\examples{
\dontrun{
library(rgee)
library(stars)
library(sf)

ee_users()
ee_Initialize(gcs = TRUE)

# Define study area (local -> earth engine)
# Communal Reserve Amarakaeri - Peru
rlist <- list(xmin = -71.13, xmax = -70.95,ymin = -12.89, ymax = -12.73)
ROI <- c(rlist$xmin, rlist$ymin,
         rlist$xmax, rlist$ymin,
         rlist$xmax, rlist$ymax,
         rlist$xmin, rlist$ymax,
         rlist$xmin, rlist$ymin)
ee_ROI <- matrix(ROI, ncol = 2, byrow = TRUE) \%>\%
  list() \%>\%
  st_polygon() \%>\%
  st_sfc() \%>\%
  st_set_crs(4326) \%>\%
  sf_as_ee()


# Get the mean annual NDVI for 2011
cloudMaskL457 <- function(image) {
  qa <- image$select("pixel_qa")
  cloud <- qa$bitwiseAnd(32L)$
    And(qa$bitwiseAnd(128L))$
    Or(qa$bitwiseAnd(8L))
  mask2 <- image$mask()$reduce(ee$Reducer$min())
  image <- image$updateMask(cloud$Not())$updateMask(mask2)
  image$normalizedDifference(list("B4", "B3"))
}

ic_l5 <- ee$ImageCollection("LANDSAT/LT05/C01/T1_SR")$
  filterBounds(ee$FeatureCollection(ee_ROI))$
  filterDate("2011-01-01", "2011-12-31")$
  map(cloudMaskL457)

# Create simple composite
mean_l5 <- ic_l5$mean()$rename("NDVI")
mean_l5 <- mean_l5$reproject(crs = "EPSG:4326", scale = 500)
mean_l5_Amarakaeri <- mean_l5$clip(ee_ROI)

# Move results from Earth Engine to Drive
task_img <- ee_image_to_gcs(
   image = mean_l5_Amarakaeri,
   bucket = "rgee_dev",
   fileFormat = "GEO_TIFF",
   region = ee_ROI,
   fileNamePrefix = "my_image_demo"
)

task_img$start()
ee_monitoring(task_img)

# Move results from Drive to local
img <- ee_gcs_to_local(task = task_img)
}
}
\seealso{
Other generic download functions: 
\code{\link{ee_drive_to_local}()}
}
\concept{generic download functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_install.R
\name{ee_install_set_pyenv}
\alias{ee_install_set_pyenv}
\title{Configure which version of Python to use with rgee}
\usage{
ee_install_set_pyenv(
  py_path = NULL,
  py_env = NULL,
  Renviron = "global",
  quiet = FALSE
)
}
\arguments{
\item{py_path}{The path to a Python interpreter}

\item{py_env}{The name of the environment}

\item{Renviron}{Character. If it is "global" the environment variables are set in
the .Renviron located in the Sys.getenv("HOME") folder. On the other hand,  if
it is "local" the environment variables are set in the .Renviron on the
working directory (getwd()). Finally, users can also enter a specific path
(see examples).}

\item{quiet}{Logical. Suppress info message}
}
\value{
no return value, called for setting EARTHENGINE_PYTHON in .Renviron
}
\description{
Configure which version of Python to use with rgee. This function creates two
environment variables: 'EARTHENGINE_PYTHON' and 'EARTHENGINE_ENV' both will be
saved into the file .Renviron.
}
\examples{
\dontrun{
library(rgee)
# ee_install_set_pyenv(py_path = "/usr/bin/python3", Renviron = "local")
# ee_install_set_pyenv(py_path = "/usr/bin/python3", Renviron = "/home/zgis/")
}
}
\seealso{
Other ee_install functions: 
\code{\link{ee_install_upgrade}()},
\code{\link{ee_install}()}
}
\concept{ee_install functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_version.R
\name{ee_version}
\alias{ee_version}
\title{Earth Engine API version}
\usage{
ee_version()
}
\value{
Character. Earth Engine Python API version used to build rgee.
}
\description{
Earth Engine API version
}
\seealso{
Other session management functions: 
\code{\link{ee_Initialize}()},
\code{\link{ee_user_info}()},
\code{\link{ee_users}()}
}
\concept{session management functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_Date.R
\name{rdate_to_eedate}
\alias{rdate_to_eedate}
\title{Pass an R date object to Earth Engine}
\usage{
rdate_to_eedate(date, timestamp = FALSE)
}
\arguments{
\item{date}{R date object}

\item{timestamp}{Logical. If TRUE, return the date in milliseconds
from the Unix Epoch (1970-01-01 00:00:00 UTC). Otherwise return a
EE date object. By default, FALSE.}
}
\value{
\code{rdate_to_eedate} will return either a numeric timestamp or
an ee$Date depending on the \code{timestamp} argument.
}
\description{
Pass an R date object ("Date", "Numeric", "character", "POSIXt",
and "POSIXct") to Google Earth Engine (ee$Date).
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()
rdate_to_eedate('2000-01-01')
rdate_to_eedate(315532800000) # float number
}
}
\seealso{
Other date functions: 
\code{\link{ee_get_date_ic}()},
\code{\link{ee_get_date_img}()},
\code{\link{eedate_to_rdate}()}
}
\concept{date functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_download.R
\name{ee_image_to_drive}
\alias{ee_image_to_drive}
\title{Creates a task to export an EE Image to Drive.}
\usage{
ee_image_to_drive(
  image,
  description = "myExportImageTask",
  folder = "rgee_backup",
  fileNamePrefix = NULL,
  timePrefix = TRUE,
  dimensions = NULL,
  region = NULL,
  scale = NULL,
  crs = NULL,
  crsTransform = NULL,
  maxPixels = NULL,
  shardSize = NULL,
  fileDimensions = NULL,
  skipEmptyTiles = NULL,
  fileFormat = NULL,
  formatOptions = NULL
)
}
\arguments{
\item{image}{The image to be exported.}

\item{description}{Human-readable name of the task.}

\item{folder}{The name of a folder in their Drive account to be
exported into. By default "rgee-backup".}

\item{fileNamePrefix}{The Google Drive filename for the export. Defaults to
the name of the task.}

\item{timePrefix}{Add current date and time as a prefix to files to export.}

\item{dimensions}{The dimensions of the exported image. It takes either a
single positive integer as the maximum dimension or "WIDTHxHEIGHT" where
WIDTH and HEIGHT are each positive integers.}

\item{region}{The lon,lat coordinates for a LinearRing or Polygon specifying
the region to export. It can be specified as nested lists of numbers or a
serialized string. Defaults to the image's region.}

\item{scale}{The resolution in meters per pixel. Defaults to the native
resolution of the image asset unless a crsTransform is specified.}

\item{crs}{The coordinate reference system of the exported image's
projection. Defaults to the image's default projection.}

\item{crsTransform}{A comma-separated string of 6 numbers describing
the affine transform of the coordinate reference system of the exported
image's projection, in the order: xScale, xShearing, xTranslation,
yShearing, yScale, and yTranslation. Defaults to the image's native
CRS transform.}

\item{maxPixels}{The maximum allowed number of pixels in the
exported image. The task will fail if the exported region covers
more pixels in the specified projection. Defaults to 100,000,000.}

\item{shardSize}{Size in pixels of the shards in which this image
will be computed. Defaults to 256.}

\item{fileDimensions}{The dimensions in pixels of each image file,
if the image is too large to fit in a single file. May specify a
single number to indicate a square shape, or a list of two dimensions
to indicate (width, height). Note that the image will still be clipped
to the overall image dimensions. Must be a multiple of shardSize.}

\item{skipEmptyTiles}{If TRUE, skip writing empty (i.e., fully-masked)
image tiles. Defaults to FALSE.}

\item{fileFormat}{The string file format to which the image is exported.
Currently only 'GeoTIFF' and 'TFRecord' are supported, defaults to 'GeoTIFF'.}

\item{formatOptions}{A dictionary of string keys to format-specific
options. **kwargs: Holds other keyword arguments that may have been
deprecated, such as 'crs_transform', 'driveFolder', and 'driveFileNamePrefix'.}
}
\value{
An unstarted Task that exports the image to Drive.
}
\description{
Creates a task to export an EE Image to Drive.
This function is a wrapper around \cr
\code{ee$batch$Export$image$toDrive(...)}.
}
\examples{
\dontrun{
library(rgee)
library(stars)
library(sf)

ee_users()
ee_Initialize(drive = TRUE)

# Define study area (local -> earth engine)
# Communal Reserve Amarakaeri - Peru
rlist <- list(xmin = -71.13, xmax = -70.95,ymin = -12.89, ymax = -12.73)
ROI <- c(rlist$xmin, rlist$ymin,
         rlist$xmax, rlist$ymin,
         rlist$xmax, rlist$ymax,
         rlist$xmin, rlist$ymax,
         rlist$xmin, rlist$ymin)

ee_ROI <- matrix(ROI, ncol = 2, byrow = TRUE) \%>\%
  list() \%>\%
  st_polygon() \%>\%
  st_sfc() \%>\%
  st_set_crs(4326) \%>\%
  sf_as_ee()


# Get the mean annual NDVI for 2011
cloudMaskL457 <- function(image) {
  qa <- image$select("pixel_qa")
  cloud <- qa$bitwiseAnd(32L)$
    And(qa$bitwiseAnd(128L))$
    Or(qa$bitwiseAnd(8L))
  mask2 <- image$mask()$reduce(ee$Reducer$min())
  image <- image$updateMask(cloud$Not())$updateMask(mask2)
  image$normalizedDifference(list("B4", "B3"))
}

ic_l5 <- ee$ImageCollection("LANDSAT/LT05/C01/T1_SR")$
  filterBounds(ee$FeatureCollection(ee_ROI))$
  filterDate("2011-01-01", "2011-12-31")$
  map(cloudMaskL457)

# Create simple composite
mean_l5 <- ic_l5$mean()$rename("NDVI")
mean_l5 <- mean_l5$reproject(crs = "EPSG:4326", scale = 500)
mean_l5_Amarakaeri <- mean_l5$clip(ee_ROI)

# Move results from Earth Engine to Drive
task_img <- ee_image_to_drive(
  image = mean_l5_Amarakaeri,
  fileFormat = "GEO_TIFF",
  region = ee_ROI,
  fileNamePrefix = "my_image_demo"
)

task_img$start()
ee_monitoring(task_img)

# Move results from Drive to local
ee_drive_to_local(task = task_img)
}
}
\seealso{
Other image export task creator: 
\code{\link{ee_image_to_asset}()},
\code{\link{ee_image_to_gcs}()}
}
\concept{image export task creator}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-upload.R
\name{local_to_gcs}
\alias{local_to_gcs}
\title{Upload local files to Google Cloud Storage}
\usage{
local_to_gcs(x, bucket = NULL, predefinedAcl = "bucketLevel", quiet = FALSE)
}
\arguments{
\item{x}{Character. filename.}

\item{bucket}{bucket name you are uploading to}

\item{predefinedAcl}{Specify user access to object. Passed to
\code{googleCloudStorageR::gcs_upload}.}

\item{quiet}{Logical. Suppress info message.}
}
\value{
Character that represents the full path of the object in the GCS
bucket specified.
}
\description{
Upload images or tables to Google Cloud Storage
}
\examples{
\dontrun{
library(rgee)
library(stars)

# Initialize a specific Earth Engine account and
# Google Cloud Storage credentials
ee_Initialize(gcs = TRUE)

# # Define an image.
tif <- system.file("tif/L7_ETMs.tif", package = "stars")
local_to_gcs(x = tif, bucket = 'rgee_dev')
}
}
\seealso{
Other generic upload functions: 
\code{\link{ee_utils_create_manifest_image}()},
\code{\link{ee_utils_create_manifest_table}()}
}
\concept{generic upload functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-upload.R
\name{ee_utils_create_manifest_image}
\alias{ee_utils_create_manifest_image}
\title{Create a manifest to upload an image}
\usage{
ee_utils_create_manifest_image(
  gs_uri,
  assetId,
  properties = NULL,
  start_time = "1970-01-01",
  end_time = "1970-01-01",
  pyramiding_policy = "MEAN",
  returnList = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{gs_uri}{Character. GCS full path of the image to upload to Earth Engine assets,
e.g. gs://rgee_dev/l8.tif}

\item{assetId}{Character. How to call the file once uploaded
to the Earth Engine Asset. e.g. users/datacolecfbf/l8.}

\item{properties}{List. Set of parameters to be set up as properties
of the EE object.}

\item{start_time}{Character. Sets the start time property (system:time_start).
It could be a number (timestamp) or a date.}

\item{end_time}{Character. Sets the end time property (system:time_end).
It could be a number (timestamp) or a date.}

\item{pyramiding_policy}{Character. The pyramid reduction policy to use.}

\item{returnList}{Logical. If TRUE will return the "manifest" as a list. Otherwise,
will return a JSON file.}

\item{quiet}{Logical. Suppress info message.}
}
\value{
If \code{returnList} is TRUE, a list otherwise a JSON file.
}
\description{
Create a manifest to upload a GeoTIFF to Earth Engine asset folder. The
"manifest" is simply a JSON file that describe all the upload parameters. See
\url{https://developers.google.com/earth-engine/guides/image_manifest} to get more
details.
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()

tif <- system.file("tif/L7_ETMs.tif", package = "stars")

# Return a JSON file
ee_utils_create_manifest_image(
  gs_uri = "gs://rgee_dev/l8.tif",
  assetId = "users/datacolecfbf/l8"
)

# Return a list
ee_utils_create_manifest_image(
  gs_uri = "gs://rgee_dev/l8.tif",
  assetId = "users/datacolecfbf/l8",
  returnList = TRUE
)
}
}
\seealso{
Other generic upload functions: 
\code{\link{ee_utils_create_manifest_table}()},
\code{\link{local_to_gcs}()}
}
\concept{generic upload functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-upload.R
\name{gcs_to_ee_table}
\alias{gcs_to_ee_table}
\title{Move a zipped shapefile from GCS to their EE Assets}
\usage{
gcs_to_ee_table(
  manifest,
  command_line_tool_path = NULL,
  overwrite = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{manifest}{Character. manifest upload file. See \code{\link{ee_utils_create_manifest_table}}.}

\item{command_line_tool_path}{Character. Path to the Earth Engine command line
tool (CLT). If NULL, rgee assumes that CLT is set in the system PATH.
(ignore if \code{via} is not defined as "gcs_to_asset").}

\item{overwrite}{Logical. If TRUE, the assetId will be overwritten if
it exists.}

\item{quiet}{Logical. Suppress info message.}
}
\value{
Character. The Earth Engine asset ID.
}
\description{
Move a zipped shapefile from GCS to their EE Assets
}
\examples{
\dontrun{
library(rgee)
library(sf)
ee_Initialize(gcs = TRUE)

# 1. Read dataset and create a output filename
x <- st_read(system.file("shape/nc.shp", package = "sf"))
assetId <- sprintf("\%s/\%s", ee_get_assethome(), 'toy_poly_gcs')

# 2. From sf to .shp
shp_dir <- sprintf("\%s.shp", tempfile())
geozip_dir <- ee_utils_shp_to_zip(x, shp_dir)

# 3. From local to gcs
gcs_filename <- local_to_gcs(
 x = geozip_dir,
 bucket = "rgee_dev" # Insert your own bucket here!
)

# 4. Create Table Manifest
manifest <- ee_utils_create_manifest_table(
 gs_uri = gcs_filename,
 assetId = assetId
)

# 5. From GCS to Earth Engine
ee_nc <- gcs_to_ee_table(manifest, overwrite = TRUE)
ee_monitoring()
Map$addLayer(ee$FeatureCollection(ee_nc))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_as_thumbnail.R
\name{ee_as_thumbnail}
\alias{ee_as_thumbnail}
\title{Create an R spatial gridded object from an EE thumbnail image}
\usage{
ee_as_thumbnail(
  image,
  region,
  dimensions,
  vizparams = NULL,
  raster = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{image}{EE Image object to be converted into a stars object.}

\item{region}{EE Geometry Rectangle (\code{ee$Geometry$Rectangle}) specifying
the region to export.The CRS needs to be the same as the \code{x} argument.
Otherwise, it will be forced.}

\item{dimensions}{Numeric vector of length 2. Thumbnail dimensions in pixel
units. If a single integer is provided, it defines the size of the
image's larger aspect dimension and scales the smaller dimension
proportionally. Defaults to 512 pixels for the larger image aspect dimension.}

\item{vizparams}{A list that contains the visualization parameters.
See details.}

\item{raster}{Logical. Should the thumbnail image be saved as a
RasterStack object?}

\item{quiet}{logical; suppress info messages.}
}
\value{
An stars or Raster object depending on the \code{raster} argument.
}
\description{
Wrapper function around \code{ee$Image$getThumbURL} to create a stars or
RasterLayer R object from a
\href{ https://developers.google.com/earth-engine/guides/image_visualization}{EE thumbnail image}.
}
\details{
\code{vizparams} set up the details of the thumbnail image. With
\code{ee_as_thumbnail} only is possible to export one-band (G) or three-band
(RGB) images. Several parameters can be passed on to control color,
intensity, the maximum and minimum values, etc. The table below provides
all the parameters that admit \code{ee_as_thumbnail}.

\tabular{lll}{
\strong{Parameter}\tab \strong{Description}  \tab \strong{Type}\cr
\strong{bands}    \tab  Comma-delimited list of three band (RGB) \tab  list \cr
\strong{min}      \tab  Value(s) to map to 0 \tab  number or list of three
numbers, one for each band \cr
\strong{max}      \tab  Value(s) to map to 1 \tab  number or list of three
numbers, one for each band \cr
\strong{gain}     \tab  Value(s) by which to multiply each pixel value \tab
number or list of three numbers, one for each band \cr
\strong{bias}     \tab  Value(s) to add to each Digital Number
value \tab number or list of three numbers, one for each band \cr
\strong{gamma}    \tab  Gamma correction factor(s) \tab  number or list of
three numbers, one for each band \cr
\strong{palette}  \tab  List of CSS-style color strings
(single-band only) \tab  comma-separated list of hex strings \cr
\strong{opacity}   \tab  The opacity of the layer (from 0 to 1)  \tab  number \cr
}
}
\examples{
\dontrun{
library(raster)
library(stars)
library(rgee)

ee_Initialize()

nc <- st_read(system.file("shp/arequipa.shp", package = "rgee"))
dem_palette <- c(
  "#008435", "#1CAC17", "#48D00C", "#B3E34B", "#F4E467",
  "#F4C84E", "#D59F3C", "#A36D2D", "#C6A889", "#FFFFFF"
)

## DEM data -SRTM v4.0
image <- ee$Image("CGIAR/SRTM90_V4")
world_region <- ee$Geometry$Rectangle(
  coords = c(-180,-60,180,60),
  proj = "EPSG:4326",
  geodesic = FALSE
)

## world - elevation
world_dem <- ee_as_thumbnail(
  image = image,
  region = world_region,
  dimensions = 1024,
  vizparams = list(min = 0, max = 5000)
)

world_dem[world_dem <= 0] <- NA
world_dem <- world_dem * 5000

plot(
  x = world_dem, col = dem_palette, breaks = "equal",
  reset = FALSE, main = "SRTM - World"
)

## Arequipa-Peru
arequipa_region <- nc \%>\%
  st_bbox() \%>\%
  st_as_sfc() \%>\%
  sf_as_ee()

arequipa_dem <- ee_as_thumbnail(
  image = image,
  region = arequipa_region$buffer(1000)$bounds(),
  dimensions = 512,
  vizparams = list(min = 0, max = 5000)
)

arequipa_dem <- arequipa_dem * 5000
st_crs(arequipa_dem) <- 4326
plot(
  x = arequipa_dem[nc], col = dem_palette, breaks = "equal",
  reset = FALSE, main = "SRTM - Arequipa"
)

suppressWarnings(plot(
  x = nc, col = NA, border = "black", add = TRUE,
  lwd = 1.5
))
dev.off()

## LANDSAT 8
img <- ee$Image("LANDSAT/LC08/C01/T1_SR/LC08_038029_20180810")$
  select(c("B4", "B3", "B2"))
Map$centerObject(img)
Map$addLayer(img, list(min = 0, max = 5000, gamma = 1.5))

## Teton Wilderness
l8_img <- ee_as_thumbnail(
  image = img,
  region = img$geometry()$bounds(),
  dimensions = 1024,
  vizparams = list(min = 0, max = 5000, gamma = 1.5),
  raster = TRUE
)
crs(l8_img) <- "+proj=longlat +datum=WGS84 +no_defs"
plotRGB(l8_img, stretch = "lin")
}
}
\seealso{
Other image download functions: 
\code{\link{ee_as_raster}()},
\code{\link{ee_as_stars}()},
\code{\link{ee_imagecollection_to_local}()}
}
\concept{image download functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-upload.R
\name{ee_utils_create_manifest_table}
\alias{ee_utils_create_manifest_table}
\title{Create a manifest to upload a table}
\usage{
ee_utils_create_manifest_table(
  gs_uri,
  assetId,
  start_time = "1970-01-01",
  end_time = "1970-01-01",
  properties = NULL,
  returnList = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{gs_uri}{Character. GCS full path of the table to upload to Earth Engine assets
e.g. gs://rgee_dev/nc.zip}

\item{assetId}{Character. How to call the file once uploaded
to the Earth Engine Asset. e.g. users/datacolecfbf/nc.}

\item{start_time}{Character. Sets the start time property (system:time_start).
It could be a number (timestamp) or a date.}

\item{end_time}{Character. Sets the end time property (system:time_end).
It could be a number (timestamp) or a date.}

\item{properties}{List. Set of parameters to be set up as properties
of the EE object.}

\item{returnList}{Logical. If TRUE will return the "manifest" as a list otherwise
will return a JSON file.}

\item{quiet}{Logical. Suppress info message.}
}
\value{
If \code{returnList} is TRUE, a list otherwise a JSON file.
}
\description{
Create a manifest to upload a zipped shapefile to Earth Engine assets folder. The
"manifest" is simply a JSON file that describe all the upload parameters. See
\url{https://developers.google.com/earth-engine/guides/image_manifest} to get more
details.
}
\examples{
\dontrun{
library(rgee)
library(sf)
ee_Initialize(gcs = TRUE)

x <- st_read(system.file("shape/nc.shp", package = "sf"))
shp_dir <- sprintf("\%s.shp", tempfile())
geozip_dir <- ee_utils_shp_to_zip(x, shp_dir)

# Return a JSON file
manifest <- ee_utils_create_manifest_table(
  gs_uri = "gs://rgee_dev/nc.zip",
  assetId = "users/datacolecfbf/nc"
)

# Return a list
ee_utils_create_manifest_table(
  gs_uri = "gs://rgee_dev/nc.zip",
  assetId = "users/datacolecfbf/nc",
  returnList = TRUE
)
}
}
\seealso{
Other generic upload functions: 
\code{\link{ee_utils_create_manifest_image}()},
\code{\link{local_to_gcs}()}
}
\concept{generic upload functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_as_sf.R
\name{ee_as_sf}
\alias{ee_as_sf}
\title{Convert an Earth Engine table in an sf object}
\usage{
ee_as_sf(
  x,
  dsn,
  overwrite = TRUE,
  via = "getInfo",
  container = "rgee_backup",
  crs = NULL,
  maxFeatures = 5000,
  selectors = NULL,
  lazy = FALSE,
  public = TRUE,
  add_metadata = TRUE,
  timePrefix = TRUE,
  quiet = FALSE
)
}
\arguments{
\item{x}{Earth Engine table (ee$FeatureCollection) to be converted in an sf
object.}

\item{dsn}{Character. Output filename. In case \code{dsn} is missing,
a shapefile is created in the \code{tmp()} directory.}

\item{overwrite}{Logical. Delete data source \code{dsn} before attempting
to write?.}

\item{via}{Character. Method to export the image. Three method are
implemented: "getInfo", "drive", "gcs". See details.}

\item{container}{Character. Name of the folder ('drive') or bucket ('gcs')
to be exported into (ignore if \code{via} is not defined as "drive" or
"gcs").}

\item{crs}{Integer or Character. Coordinate Reference System (CRS)
for the EE table. If it is NULL, \code{ee_as_sf} will take the CRS of
the first element.}

\item{maxFeatures}{Numeric. The maximum allowed number of features to
export (ignore if \code{via} is not set as "getInfo"). The task will fail
if the exported region covers more features than the specified in
\code{maxFeatures}. Defaults to 5000.}

\item{selectors}{The list of properties to include in the output, as a
list/vector of strings or a comma-separated string. By default, all properties are
included.}

\item{lazy}{Logical. If TRUE, a \code{\link[future:sequential]{
future::sequential}} object is created to evaluate the task in the future.
Ignore if \code{via} is set as "getInfo". See details.}

\item{public}{Logical. If TRUE, a public link to the file is created.
See details.}

\item{add_metadata}{Add metadata to the sf object. See details.}

\item{timePrefix}{Logical. Add current date and time (\code{Sys.time()}) as
a prefix to export files. This parameter helps to avoid exported files
with the same name. By default TRUE.}

\item{quiet}{logical. Suppress info message.}
}
\value{
An sf object.
}
\description{
Convert an Earth Engine table in an sf object
}
\details{
\code{ee_as_sf} supports the download of \code{ee$Geometry}, \code{ee$Feature},
and \code{ee$FeatureCollection} by three different options:
"getInfo" (which make an REST call to retrieve the data), "drive"
(which use \href{https://CRAN.R-project.org/package=googledrive}{Google Drive})
and "gcs" (which use \href{https://CRAN.R-project.org/package=googleCloudStorageR}{
Google Cloud Storage}). The advantage of use "getInfo" is a
direct and faster download. However, there is a limitation of 5000 features by
request, making it not recommendable for large FeatureCollection. Instead of
"getInfo", the options: "drive" and "gcs" are suitable for large FeatureCollections
due to the use of an intermediate container. When via is set as "drive" or "gcs"
\code{ee_as_sf} perform the following steps:
\itemize{
\item{1. }{A task is started (i.e., \code{ee$batch$Task$start()}) to
move the EE Table from Earth Engine to the file storage system (Google Drive
or Google Cloud Storage) specified in the argument \code{via}.}
\item{2. }{If the argument \code{lazy} is TRUE, the task will not be
monitored. This is useful to lunch several tasks simultaneously and
calls them later using \code{\link{ee_utils_future_value}} or
\code{\link[future:value]{future::value}}. At the end of this step,
the EE Table is stored on the path specified in the argument
\code{dsn}.}
\item{3. }{Finally, if the argument \code{add_metadata} is TRUE, a list
with the following elements is added to the sf object.
\itemize{
\item{\bold{if via is "drive":}}
\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{drive_name: }}{Name of the Table in Google Drive.}
\item{\bold{drive_id: }}{Id of the Table in Google Drive.}
\item{\bold{drive_download_link: }}{Download link to the table.}
}
}
\itemize{
\item{\bold{if via is "gcs":}}
\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{gcs_name: }}{Name of the Table in Google Cloud Storage.}
\item{\bold{gcs_bucket: }}{Name of the bucket.}
\item{\bold{gcs_fileFormat: }}{Format of the table.}
\item{\bold{gcs_public_link: }}{Download link to the table.}
\item{\bold{gcs_URI: }}{gs:// link to the table.}
}
}
Run \code{attr(sf, "metadata")} to get the list.
}
}

For getting more information about exporting data from Earth Engine, take
a look at the
\href{https://developers.google.com/earth-engine/guides/exporting}{Google
Earth Engine Guide - Export data}.
}
\examples{
\dontrun{
library(rgee)

ee_Initialize(drive = TRUE, gcs = TRUE)

# Region of interest
roi <- ee$Geometry$Polygon(list(
  c(-122.275, 37.891),
  c(-122.275, 37.868),
  c(-122.240, 37.868),
  c(-122.240, 37.891)
))

# TIGER: US Census Blocks Dataset
blocks <- ee$FeatureCollection("TIGER/2010/Blocks")
subset <- blocks$filterBounds(roi)
sf_subset <- ee_as_sf(x = subset)
plot(sf_subset)

# Create Random points in Earth Engine
region <- ee$Geometry$Rectangle(-119.224, 34.669, -99.536, 50.064)
ee_help(ee$FeatureCollection$randomPoints)
ee_randomPoints <- ee$FeatureCollection$randomPoints(region, 100)

# Download via GetInfo
sf_randomPoints <- ee_as_sf(ee_randomPoints)
plot(sf_randomPoints)

# Download via drive
sf_randomPoints_drive <- ee_as_sf(
  x = ee_randomPoints,
  via = 'drive'
)

# Download via GCS
sf_randomPoints_gcs <- ee_as_sf(
 x = subset,
 via = 'gcs',
 container = 'rgee_dev' #GCS bucket name
)
}
}
\concept{vector download functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_imagecollection.R
\name{ee_imagecollection_to_local}
\alias{ee_imagecollection_to_local}
\title{Save an EE ImageCollection in their local system}
\usage{
ee_imagecollection_to_local(
  ic,
  region,
  dsn = NULL,
  via = "drive",
  container = "rgee_backup",
  scale = NULL,
  maxPixels = 1e+09,
  lazy = FALSE,
  public = TRUE,
  add_metadata = TRUE,
  timePrefix = TRUE,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{ic}{ee$ImageCollection to be saved in the system.}

\item{region}{EE Geometry (ee$Geometry$Polygon). The
CRS needs to be the same that the \code{ic} argument. Otherwise, it will be
forced.}

\item{dsn}{Character. Output filename. If missing, a temporary file will
be created for each image.}

\item{via}{Character. Method to export the image. Two methods are implemented:
"drive", "gcs". See details.}

\item{container}{Character. Name of the folder ('drive') or bucket ('gcs')
to be exported into (ignored if \code{via} is not defined as "drive" or
"gcs").}

\item{scale}{Numeric. The resolution in meters per pixel. Defaults
to the native resolution of the image.}

\item{maxPixels}{Numeric. The maximum allowed number of pixels in the
exported image. The task will fail if the exported region covers
more pixels in the specified projection. Defaults to 100,000,000.}

\item{lazy}{Logical. If TRUE, a \code{\link[future:sequential]{
future::sequential}} object is created to evaluate the task in the future.
See details.}

\item{public}{Logical. If TRUE, a public link to the image is created.}

\item{add_metadata}{Add metadata to the stars_proxy object. See details.}

\item{timePrefix}{Logical. Add current date and time (\code{Sys.time()}) as
a prefix to export files. This parameter helps to avoid exported files
with the same name. By default TRUE.}

\item{quiet}{Logical. Suppress info message}

\item{...}{Extra exporting argument. See \link{ee_image_to_drive} and}
}
\value{
If add_metadata is FALSE, \code{ee_imagecollection_to_local} will
return a character vector containing the filename of the images downloaded.
Otherwise, if add_metadata is TRUE, will return a list with extra information
related to the exportation (see details).
}
\description{
Save an EE ImageCollection in their local system
}
\details{
\code{ee_imagecollection_to_local} supports the download of \code{ee$Images}
by two different options: "drive"
(\href{https://CRAN.R-project.org/package=googledrive}{Google Drive}) and "gcs"
(\href{https://CRAN.R-project.org/package=googleCloudStorageR}{
Google Cloud Storage}). In both cases, \code{ee_imagecollection_to_local}
works as follow:
\itemize{
\item{1. }{A task is started (i.e., \code{ee$batch$Task$start()}) to
move the \code{ee$Image} from Earth Engine to the intermediate container
specified in the argument \code{via}.}
\item{2. }{If the argument \code{lazy} is TRUE, the task will not be
monitored. This is useful to lunch several tasks simultaneously and
calls them later using \code{\link{ee_utils_future_value}} or
\code{\link[future:value]{future::value}}. At the end of this step,
the \code{ee$Images} are stored on the path specified in the argument
\code{dsn}.}
\item{3. }{Finally, if the argument \code{add_metadata} is TRUE, a list
with the following elements will be added to the argument \code{dsn}.
\itemize{
\item{\bold{if via is "drive":}}
\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{drive_name: }}{Name of the Image in Google Drive.}
\item{\bold{drive_id: }}{Id of the Image in Google Drive.}
\item{\bold{drive_download_link: }}{Download link to the image.}
}
}
\itemize{
\item{\bold{if via is "gcs":}}
\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{gcs_name: }}{Name of the Image in Google Cloud Storage.}
\item{\bold{gcs_bucket: }}{Name of the bucket.}
\item{\bold{gcs_fileFormat: }}{Format of the image.}
\item{\bold{gcs_public_link: }}{Download link to the image.}
\item{\bold{gcs_URI: }}{gs:// link to the image.}
}
}
}
}

For getting more information about exporting data from Earth Engine, take
a look at the
\href{https://developers.google.com/earth-engine/guides/exporting}{Google
Earth Engine Guide - Export data}.
}
\examples{
\dontrun{
library(rgee)
library(raster)
ee_Initialize(drive = TRUE, gcs = TRUE)

# USDA example
loc <- ee$Geometry$Point(-99.2222, 46.7816)
collection <- ee$ImageCollection('USDA/NAIP/DOQQ')$
  filterBounds(loc)$
  filterDate('2008-01-01', '2020-01-01')$
  filter(ee$Filter$listContains("system:band_names", "N"))

# From ImageCollection to local directory
ee_crs <- collection$first()$projection()$getInfo()$crs
geometry <- collection$first()$geometry(proj = ee_crs)$bounds()
tmp <- tempdir()

## Using drive
# one by once
ic_drive_files_1 <- ee_imagecollection_to_local(
  ic = collection,
  region = geometry,
  scale = 250,
  dsn = file.path(tmp, "drive_")
)

# all at once
ic_drive_files_2 <- ee_imagecollection_to_local(
  ic = collection,
  region = geometry,
  scale = 250,
  lazy = TRUE,
  dsn = file.path(tmp, "drive_")
)

# From Google Drive to client-side
doqq_dsn <- ic_drive_files_2 \%>\% ee_utils_future_value()
sapply(doqq_dsn, '[[', 1)
}
}
\seealso{
Other image download functions: 
\code{\link{ee_as_raster}()},
\code{\link{ee_as_stars}()},
\code{\link{ee_as_thumbnail}()}
}
\concept{image download functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_utils_sak_validate}
\alias{ee_utils_sak_validate}
\title{Validate a Service account key (SaK)}
\usage{
ee_utils_sak_validate(sakfile, bucket = NULL, quiet = FALSE)
}
\arguments{
\item{sakfile}{Character. SaK filename.}

\item{bucket}{Character. Name of the GCS bucket. If bucket is not set,
rgee will tries to create a bucket using \code{googleCloudStorageR::gcs_create_bucket}.}

\item{quiet}{Logical. Suppress info message}
}
\description{
Validate a Service account key (SaK). local_to_gcs, raster_as_ee,
stars_as_ee, and sf_as_ee(via = "gcs_to_asset", ...) need that the SaK
have privileges to write/read objects in a GCS bucket.
}
\examples{
\dontrun{
library(rgee)

ee_Initialize(gcs = TRUE)

# Check a specific SaK
sakfile <- "/home/rgee_dev/sak_file.json"
ee_utils_sak_validate(sakfile, bucket = "rgee_dev")

# Check the SaK for the current user
ee_utils_sak_validate()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_extract.R
\name{ee_extract}
\alias{ee_extract}
\title{Extract values from EE Images or ImageCollections objects}
\usage{
ee_extract(
  x,
  y,
  fun = ee$Reducer$mean(),
  scale = NULL,
  sf = FALSE,
  via = "getInfo",
  container = "rgee_backup",
  lazy = FALSE,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{x}{ee$Image.}

\item{y}{ee$Geometry$*, ee$Feature, ee$FeatureCollection, sfc or sf objects.}

\item{fun}{ee$Reducer object. Function to summarize the values. The function
must take a single numeric value as an argument and return a single value.
See details.}

\item{scale}{A nominal scale in meters of the Image projection to work in.
By default 1000.}

\item{sf}{Logical. Should return an sf object?}

\item{via}{Character. Method to export the image. Three method are
implemented: "getInfo", "drive", "gcs".}

\item{container}{Character. Name of the folder ('drive') or bucket ('gcs')
to be exported into (ignore if \code{via} is not defined as "drive" or
"gcs").}

\item{lazy}{Logical. If TRUE, a \code{\link[future:sequential]{
future::sequential}} object is created to evaluate the task in the future.
Ignore if \code{via} is set as "getInfo". See details.}

\item{quiet}{Logical. Suppress info message.}

\item{...}{ee$Image$reduceRegions additional parameters. See
\code{ee_help(ee$Image$reduceRegions)} for more details.}
}
\value{
A data.frame or an sf object depending on the sf argument.
Column names are extracted from band names. Use \code{ee$Image$rename} to
rename the bands of an \code{ee$Image}. See \code{ee_help(ee$Image$rename)}.
}
\description{
Extract values from an \code{ee$Image} at the
locations of a geometry object. Users can use \code{ee$Geometry$*},
\code{ee$Feature}, \code{ee$FeatureCollection}, sf or sfc object to filter
spatially. This function mimicking how \link[raster]{extract} currently works.
}
\details{
The reducer functions that return one value are:
\itemize{
\item  \strong{allNonZero}: Returns a Reducer that returns 1 if all of its
inputs are non-zero, 0 otherwise. \cr
\item \strong{anyNonZero}: Returns a Reducer that returns 1 if any of its
inputs are non-zero, 0 otherwise. \cr
\item \strong{bitwiseAnd}: Returns a Reducer that computes the bitwise-and
summation of its inputs.
\item \strong{bitwiseOr}: Returns a Reducer that computes the bitwise-or
summation of its inputs.
\item \strong{count}: Returns a Reducer that computes the number of
non-null inputs.
\item \strong{first}: Returns a Reducer that returns the first of its inputs.
\item \strong{firstNonNull}: Returns a Reducer that returns the first of
its non-null inputs.
\item \strong{kurtosis}: Returns a Reducer that Computes the kurtosis of
its inputs.
\item \strong{last}: Returns a Reducer that returns the last of its inputs.
\item \strong{lastNonNull}: Returns a Reducer that returns the last of its
non-null inputs.
\item \strong{max}: Creates a reducer that outputs the maximum value of
its (first) input. If numInputs is greater than one, also outputs the
corresponding values of the additional inputs.
\item \strong{mean}: Returns a Reducer that computes the (weighted)
arithmetic mean of its inputs.
\item \strong{median}: Create a reducer that will compute the median of
the inputs. For small numbers of inputs (up to maxRaw) the median will be
computed directly; for larger numbers of inputs the median will be derived
from a histogram.
\item \strong{min}: Creates a reducer that outputs the minimum value
of its (first) input.  If numInputs is greater than one, also outputs
additional inputs.
\item \strong{mode}: Create a reducer that will compute the mode of the
inputs.  For small numbers of inputs (up to maxRaw) the mode will be
computed directly; for larger numbers of inputs the mode will be derived
from a histogram.
\item \strong{product}: Returns a Reducer that computes the product of
its inputs.
\item \strong{sampleStdDev}: Returns a Reducer that computes the sample
standard deviation of its inputs.
\item \strong{sampleVariance}: Returns a Reducer that computes the sample
variance of its inputs.
\item \strong{stdDev}: Returns a Reducer that computes the standard
deviation of its inputs.
\item \strong{sum}: Returns a Reducer that computes the (weighted) sum
of its inputs.
\item \strong{variance}: Returns a Reducer that computes the variance
of its inputs.
}
}
\examples{
\dontrun{
library(rgee)
library(sf)

ee_Initialize(gcs = TRUE, drive = TRUE)

# Define a Image or ImageCollection: Terraclimate
terraclimate <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") \%>\%
 ee$ImageCollection$filterDate("2001-01-01", "2002-01-01") \%>\%
ee$ImageCollection$map(
   function(x) {
     date <- ee$Date(x$get("system:time_start"))$format('YYYY_MM_dd')
     name <- ee$String$cat("Terraclimate_pp_", date)
     x$select("pr")$rename(name)
   }
 )

# Define a geometry
nc <- st_read(
 dsn = system.file("shape/nc.shp", package = "sf"),
 stringsAsFactors = FALSE,
 quiet = TRUE
)


#Extract values - getInfo
ee_nc_rain <- ee_extract(
 x = terraclimate,
 y = nc["NAME"],
 scale = 250,
 fun = ee$Reducer$mean(),
 sf = TRUE
)

# Extract values - drive (lazy = TRUE)
ee_nc_rain <- ee_extract(
 x = terraclimate,
 y = nc["NAME"],
 scale = 250,
 fun = ee$Reducer$mean(),
 via = "drive",
 lazy = TRUE,
 sf = TRUE
)
ee_nc_rain <- ee_nc_rain \%>\% ee_utils_future_value()

# Extract values - gcs (lazy = FALSE)
ee_nc_rain <- ee_extract(
 x = terraclimate,
 y = nc["NAME"],
 scale = 250,
 fun = ee$Reducer$mean(),
 via = "gcs",
 container = "rgee_dev",
 sf = TRUE
)

# Spatial plot
plot(
 ee_nc_rain["X200101_Terraclimate_pp_2001_01_01"],
 main = "2001 Jan Precipitation - Terraclimate",
 reset = FALSE
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_install.R
\name{ee_install}
\alias{ee_install}
\title{Create an isolated Python virtual environment with all rgee dependencies.}
\usage{
ee_install(
  py_env = "rgee",
  earthengine_version = ee_version(),
  python_version = "3.8",
  confirm = interactive()
)
}
\arguments{
\item{py_env}{Character. The name, or full path, of the Python environment
to be used by rgee.}

\item{earthengine_version}{Character. The Earth Engine Python API version
to install. By default \code{rgee::ee_version()}.}

\item{python_version}{Only windows users. The version of Python to be used in
this conda environment. The associated Python package from conda will be requested
as python={python_version}. When NULL, the default python package will be
used instead. For example, use python_version = "3.6" to request that the
conda environment be created with a copy of Python 3.6. This argument will be
ignored if python is specified as part of the packages argument, for backwards
compatibility.}

\item{confirm}{Logical. Confirm before restarting R?.}
}
\value{
No return value, called for installing non-R dependencies.
}
\description{
Create an isolated Python virtual environment with all rgee dependencies.
\code{ee_install} realize the following six (6) tasks:
\itemize{
\item{1. }{If you do not count with a Python environment, it will display
an interactive menu to install \href{https://docs.conda.io/en/latest/miniconda.html}{Miniconda}
(a free minimal installer for conda).}
\item{2. }{Remove the previous Python environment defined in \code{py_env} if
it exist.}
\item{3. }{Create a new Python environment (See \code{py_env}).}
\item{4. }{ Set the environment variable EARTHENGINE_PYTHON. It is used to
define RETICULATE_PYTHON when the library is loaded. See this
\href{https://rstudio.github.io/reticulate/articles/versions.html}{article}
for further details.
}
\item{5. }{Install rgee Python dependencies. Using
\code{reticulate::py_install}.}
\item{6. }{Interactive menu to confirm if restart the R session to see
changes.}
}
}
\seealso{
Other ee_install functions: 
\code{\link{ee_install_set_pyenv}()},
\code{\link{ee_install_upgrade}()}
}
\concept{ee_install functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_module.R
\docType{data}
\name{ee}
\alias{ee}
\title{Main Earth Engine module}
\format{
Earth Engine module
}
\usage{
ee
}
\description{
Interface to main Earth Engine module. Provides access to the top level
classes and functions as well as sub-modules (e.g. \code{ee$Image},
\code{ee$FeatureCollection$first}, etc.).
}
\examples{
\dontrun{
library(rgee)

ee_Initialize()

ee_img <- ee$Image(0)
ee_ic <- ee$ImageCollection(ee_img)

print(ee_img$getInfo())
print(ee_ic$getInfo())
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-upload.R
\name{gcs_to_ee_image}
\alias{gcs_to_ee_image}
\title{Move a GeoTIFF image from GCS to their EE assets}
\usage{
gcs_to_ee_image(
  manifest,
  overwrite = FALSE,
  command_line_tool_path = NULL,
  quiet = FALSE
)
}
\arguments{
\item{manifest}{Character. Manifest upload file. See \code{\link{ee_utils_create_manifest_image}}.}

\item{overwrite}{Logical. If TRUE, the assetId will be overwritten if
it exists.}

\item{command_line_tool_path}{Character. Path to the Earth Engine command line
tool (CLT). If NULL, rgee assumes that CLT is set in the system PATH.
(ignore if \code{via} is not defined as "gcs_to_asset").}

\item{quiet}{Logical. Suppress info message.}
}
\value{
Character. The Earth Engine asset ID.
}
\description{
Move a GeoTIFF image from GCS to their EE assets
}
\examples{
\dontrun{
library(rgee)
library(stars)
ee_Initialize("csaybar", gcs = TRUE)

# 1. Read GeoTIFF file and create a output filename
tif <- system.file("tif/L7_ETMs.tif", package = "stars")
x <- read_stars(tif)
assetId <- sprintf("\%s/\%s",ee_get_assethome(),'stars_l7')

# 2. From local to gcs
gs_uri <- local_to_gcs(
  x = tif,
  bucket = 'rgee_dev' # Insert your own bucket here!
)

# 3. Create an Image Manifest
manifest <- ee_utils_create_manifest_image(gs_uri, assetId)

# 4. From GCS to Earth Engine
gcs_to_ee_image(
  manifest = manifest,
  overwrite = TRUE
)

# OPTIONAL: Monitoring progress
ee_monitoring()

# OPTIONAL: Display results
ee_stars_01 <- ee$Image(assetId)
ee_stars_01$bandNames()$getInfo()

Map$centerObject(ee_stars_01)
Map$addLayer(ee_stars_01, list(min = 0, max = 255, bands = c("b3", "b2", "b1")))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_check.R
\name{ee_check-tools}
\alias{ee_check-tools}
\alias{ee_check}
\alias{ee_check_python}
\alias{ee_check_python_packages}
\alias{ee_check_credentials}
\title{Interface to check Python and non-R dependencies}
\usage{
ee_check(user = NULL, quiet = FALSE)

ee_check_python(quiet = FALSE)

ee_check_python_packages(quiet = FALSE)

ee_check_credentials(quiet = FALSE)
}
\arguments{
\item{user}{Character. User to check credentials. If it is not defined,
ee_check will skip the check of credentials.}

\item{quiet}{Logical. Suppress info message}
}
\value{
No return value, called for checking non-R rgee dependencies.
}
\description{
R functions for checking Google credentials (Google Earth Engine,
Google Drive and Google Cloud Storage), Python environment and
Third-Party Python Packages used by rgee.
}
\examples{
\dontrun{
library(rgee)

ee_check_python()
ee_check_python_packages()
ee_check_credentials()
ee_check() # put them all together
}
}
\concept{ee_check functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6Map.R
\name{R6Map}
\alias{R6Map}
\title{R6 class to display Earth Engine (EE) spatial objects}
\value{
Object of class \code{leaflet} and \code{EarthEngineMap}, with the
following extra parameters: tokens, name, opacity, shown, min, max, palette,
position, and legend. Use the $ method to retrieve the data (e.g., m$rgee$min).
}
\description{
Create interactive visualizations of spatial EE objects
(ee$Geometry, ee$Image, ee$Feature, and ee$FeatureCollection)
using \code{leaflet}.
}
\details{
\code{R6Map} uses the Earth Engine method
\href{https://developers.google.com/earth-engine/api_docs#ee.data.getmapid/}{
getMapId} to fetch and return an ID dictionary used to create
layers in a \code{leaflet} object. Users can specify visualization
parameters to Map$addLayer by using the visParams argument. Each Earth
Engine spatial object has a specific format. For
\code{ee$Image}, the
\href{https://developers.google.com/earth-engine/guides/image_visualization}{
parameters} available are:

\tabular{lll}{
\strong{Parameter}\tab \strong{Description}  \tab \strong{Type}\cr
\strong{bands}    \tab  Comma-delimited list of three band (RGB) \tab  list \cr
\strong{min}      \tab  Value(s) to map to 0 \tab  number or list of three
numbers, one for each band \cr
\strong{max}      \tab  Value(s) to map to 1 \tab  number or list of three
numbers, one for each band \cr
\strong{gain}     \tab  Value(s) by which to multiply each pixel value \tab
number or list of three numbers, one for each band \cr
\strong{bias}     \tab  Value(s) to add to each Digital Number
value \tab number or list of three numbers, one for each band \cr
\strong{gamma}    \tab  Gamma correction factor(s) \tab  number or list of
three numbers, one for each band \cr
\strong{palette}  \tab  List of CSS-style color strings
(single-band only) \tab  comma-separated list of hex strings \cr
\strong{opacity}   \tab  The opacity of the layer (from 0 to 1)  \tab  number \cr
}

If you add an \code{ee$Image} to Map$addLayer without any additional
parameters. By default it assigns the first three bands to red,
green, and blue bands, respectively. The default stretch is based on the
min-max range. On the other hand, the available parameters for
\code{ee$Geometry}, \code{ee$Feature}, and \code{ee$FeatureCollection}
are:

\itemize{
\item \strong{color}: A hex string in the format RRGGBB specifying the
color to use for drawing the features. By default #000000.
\item \strong{pointRadius}: The radius of the point markers. By default 3.
\item \strong{strokeWidth}: The width of lines and polygon borders. By
default 3.
}
}
\examples{

## ------------------------------------------------
## Method `R6Map$reset`
## ------------------------------------------------

\dontrun{
library(rgee)
ee_Initialize()

# Load an Image
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")

# Create
Map <- R6Map$new()
Map$centerObject(image)

# Simple display: Map just will
Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max = 10000, bands = c("B4", "B3", "B2")),
  name = "l8_01"
)
Map # display map

Map$reset() # Reset arguments
Map
}

## ------------------------------------------------
## Method `R6Map$setCenter`
## ------------------------------------------------

\dontrun{
library(rgee)

ee_Initialize()

Map <- R6Map$new()
Map$setCenter(lon = -76, lat = 0, zoom = 5)
Map

# Map$lat
# Map$lon
# Map$zoom
}

## ------------------------------------------------
## Method `R6Map$setZoom`
## ------------------------------------------------

\dontrun{
library(rgee)

ee_Initialize()

Map <- R6Map$new()
Map$setZoom(zoom = 4)
Map

# Map$lat
# Map$lon
# Map$zoom
}

## ------------------------------------------------
## Method `R6Map$centerObject`
## ------------------------------------------------

\dontrun{
library(rgee)

ee_Initialize()

Map <- R6Map$new()
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")
Map$centerObject(image)
Map
}

## ------------------------------------------------
## Method `R6Map$addLayer`
## ------------------------------------------------

\dontrun{
library(rgee)
ee_Initialize()

# Load an Image
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")

# Create
Map <- R6Map$new()
Map$centerObject(image)

# Simple display: Map just will
Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max = 10000, bands = c("B4", "B3", "B2")),
  name = "l8_01"
)

Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max = 20000, bands = c("B4", "B3", "B2")),
  name = "l8_02"
)

# Simple display: Map just will (if the position is not specified it will
# be saved on the right side)
Map$reset() # Reset Map to the initial arguments.
Map$centerObject(image)
Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max=10000, bands = c("B4", "B3", "B2")),
  name = "l8_left",
  position = "left"
)

Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max=20000, bands = c("B4", "B3", "B2")),
  name = "l8_right"
)

Map$reset()
}

## ------------------------------------------------
## Method `R6Map$addLayers`
## ------------------------------------------------

\dontrun{
library(sf)
library(rgee)
library(rgeeExtra)

ee_Initialize()

Map <- R6Map$new()

nc <- st_read(system.file("shape/nc.shp", package = "sf")) \%>\%
  st_transform(4326) \%>\%
  sf_as_ee()

ee_s2 <- ee$ImageCollection("COPERNICUS/S2")$
  filterDate("2016-01-01", "2016-01-31")$
  filterBounds(nc) \%>\%
  ee_get(0:2)

Map$centerObject(nc$geometry())
Map$addLayers(eeObject = ee_s2,position = "right")

# digging up the metadata
Map$previous_map_right$rgee$tokens

Map$reset()
}

## ------------------------------------------------
## Method `R6Map$addLegend`
## ------------------------------------------------

\dontrun{
library(leaflet)
library(rgee)
ee_Initialize()

Map$reset()

# Load MODIS ImageCollection
imgcol <- ee$ImageCollection$Dataset$MODIS_006_MOD13Q1

# Parameters for visualization
labels <- c("good", "marginal", "snow", "cloud")
cols   <- c("#999999", "#00BFC4", "#F8766D", "#C77CFF")
vis_qc <- list(min = 0, max = 3, palette = cols, bands = "SummaryQA", values = labels)

# Create interactive map
m_qc <- Map$addLayer(imgcol$median(), vis_qc, "QC")

# continous palette
Map$addLegend(vis_qc)

# categorical palette
Map$addLegend(vis_qc, name = "Legend1", color_mapping = "discrete")

# character palette
Map$addLegend(vis_qc, name = "Legend2", color_mapping = "character")
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{lon}}{The longitude of the center, in degrees.}

\item{\code{lat}}{The latitude of the center, in degrees.}

\item{\code{zoom}}{The zoom level, from 1 to 24.}

\item{\code{save_maps}}{Should \code{R6Map} save the previous maps?. If TRUE, Map
will work in an OOP style. Otherwise it will be a functional programming
style.}

\item{\code{previous_map_left}}{Container on maps in the left side.}

\item{\code{previous_map_right}}{Container on maps in the right side.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{R6Map$new()}}
\item \href{#method-reset}{\code{R6Map$reset()}}
\item \href{#method-print}{\code{R6Map$print()}}
\item \href{#method-setCenter}{\code{R6Map$setCenter()}}
\item \href{#method-setZoom}{\code{R6Map$setZoom()}}
\item \href{#method-centerObject}{\code{R6Map$centerObject()}}
\item \href{#method-addLayer}{\code{R6Map$addLayer()}}
\item \href{#method-addLayers}{\code{R6Map$addLayers()}}
\item \href{#method-addLegend}{\code{R6Map$addLegend()}}
\item \href{#method-clone}{\code{R6Map$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Constructor of R6Map.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$new(lon = 0, lat = 0, zoom = 1, save_maps = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lon}}{The longitude of the center, in degrees. By default -76.942478.}

\item{\code{lat}}{The latitude of the center, in degrees. By default -12.172116.}

\item{\code{zoom}}{The zoom level, from 1 to 24. By default 18.}

\item{\code{save_maps}}{Should \code{R6Map} save previous maps?.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{EarthEngineMap} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-reset"></a>}}
\if{latex}{\out{\hypertarget{method-reset}{}}}
\subsection{Method \code{reset()}}{
Reset to initial arguments.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$reset(lon = 0, lat = 0, zoom = 1, save_maps = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lon}}{The longitude of the center, in degrees. By default -76.942478.}

\item{\code{lat}}{The latitude of the center, in degrees. By default -12.172116.}

\item{\code{zoom}}{The zoom level, from 1 to 24. By default 18.}

\item{\code{save_maps}}{Should \code{R6Map} save previous maps?.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{EarthEngineMap} object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
library(rgee)
ee_Initialize()

# Load an Image
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")

# Create
Map <- R6Map$new()
Map$centerObject(image)

# Simple display: Map just will
Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max = 10000, bands = c("B4", "B3", "B2")),
  name = "l8_01"
)
Map # display map

Map$reset() # Reset arguments
Map
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
Display a \code{EarthEngineMap} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$print()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
An \code{EarthEngineMap} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-setCenter"></a>}}
\if{latex}{\out{\hypertarget{method-setCenter}{}}}
\subsection{Method \code{setCenter()}}{
Centers the map view at the given coordinates with the given zoom level. If
no zoom level is provided, it uses 10 by default.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$setCenter(lon = 0, lat = 0, zoom = 10)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lon}}{The longitude of the center, in degrees. By default -76.942478.}

\item{\code{lat}}{The latitude of the center, in degrees. By default -12.172116.}

\item{\code{zoom}}{The zoom level, from 1 to 24. By default 18.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
No return value, called to set initial coordinates and zoom.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
library(rgee)

ee_Initialize()

Map <- R6Map$new()
Map$setCenter(lon = -76, lat = 0, zoom = 5)
Map

# Map$lat
# Map$lon
# Map$zoom
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-setZoom"></a>}}
\if{latex}{\out{\hypertarget{method-setZoom}{}}}
\subsection{Method \code{setZoom()}}{
Sets the zoom level of the map.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$setZoom(zoom = 10)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{zoom}}{The zoom level, from 1 to 24. By default 10.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
No return value, called to set zoom.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
library(rgee)

ee_Initialize()

Map <- R6Map$new()
Map$setZoom(zoom = 4)
Map

# Map$lat
# Map$lon
# Map$zoom
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-centerObject"></a>}}
\if{latex}{\out{\hypertarget{method-centerObject}{}}}
\subsection{Method \code{centerObject()}}{
Centers the map view on a given object. If no zoom level is provided, it
will be predicted according to the bounds of the Earth Engine object
specified.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$centerObject(
  eeObject,
  zoom = NULL,
  maxError = ee$ErrorMargin(1),
  titiler_server = "https://api.cogeo.xyz/"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{eeObject}}{Earth Engine spatial object.}

\item{\code{zoom}}{The zoom level, from 1 to 24. By default NULL.}

\item{\code{maxError}}{Max error when input image must be reprojected to an
explicitly requested result projection or geodesic state.}

\item{\code{titiler_server}}{TiTiler endpoint. Defaults to "https://api.cogeo.xyz/".}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
No return value, called to set zoom.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
library(rgee)

ee_Initialize()

Map <- R6Map$new()
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")
Map$centerObject(image)
Map
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-addLayer"></a>}}
\if{latex}{\out{\hypertarget{method-addLayer}{}}}
\subsection{Method \code{addLayer()}}{
Adds a given Earth Engine spatial object to the map as a layer
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$addLayer(
  eeObject,
  visParams = NULL,
  name = NULL,
  shown = TRUE,
  opacity = 1,
  position = NULL,
  titiler_viz_convert = TRUE,
  titiler_server = "https://api.cogeo.xyz/"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{eeObject}}{The Earth Engine spatial object to display in the interactive map.}

\item{\code{visParams}}{List of parameters for visualization. See details.}

\item{\code{name}}{The name of layers.}

\item{\code{shown}}{A flag indicating whether layers should be on by default.}

\item{\code{opacity}}{The layer's opacity is represented as a number between 0 and 1. Defaults to 1.}

\item{\code{position}}{Character. Activate panel creation. If "left" the map will be displayed in
the left panel. Otherwise, if it is "right" the map will be displayed in the right panel.
By default NULL (No panel will be created).}

\item{\code{titiler_viz_convert}}{Logical. If it is TRUE, Map$addLayer will transform the
visParams to titiler style. Ignored if eeObject is not a COG file.}

\item{\code{titiler_server}}{TiTiler endpoint. Defaults to "https://api.cogeo.xyz/".}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An \code{EarthEngineMap} object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
library(rgee)
ee_Initialize()

# Load an Image
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")

# Create
Map <- R6Map$new()
Map$centerObject(image)

# Simple display: Map just will
Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max = 10000, bands = c("B4", "B3", "B2")),
  name = "l8_01"
)

Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max = 20000, bands = c("B4", "B3", "B2")),
  name = "l8_02"
)

# Simple display: Map just will (if the position is not specified it will
# be saved on the right side)
Map$reset() # Reset Map to the initial arguments.
Map$centerObject(image)
Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max=10000, bands = c("B4", "B3", "B2")),
  name = "l8_left",
  position = "left"
)

Map$addLayer(
  eeObject = image,
  visParams = list(min=0, max=20000, bands = c("B4", "B3", "B2")),
  name = "l8_right"
)

Map$reset()
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-addLayers"></a>}}
\if{latex}{\out{\hypertarget{method-addLayers}{}}}
\subsection{Method \code{addLayers()}}{
Adds a given ee$ImageCollection to the map as multiple layers.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$addLayers(
  eeObject,
  visParams = NULL,
  nmax = 5,
  name = NULL,
  shown = TRUE,
  position = NULL,
  opacity = 1
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{eeObject}}{ee$ImageCollection to display in the interactive map.}

\item{\code{visParams}}{List of parameters for visualization. See details.}

\item{\code{nmax}}{Numeric. The maximum number of images to display. By default 5.}

\item{\code{name}}{The name of layers.}

\item{\code{shown}}{A flag indicating whether layers should be on by default.}

\item{\code{position}}{Character. Activate panel creation. If "left" the map will be displayed in
the left panel. Otherwise, if it is "right" the map will be displayed in the right panel.
By default NULL (No panel will be created).}

\item{\code{opacity}}{The layer's opacity is represented as a number between 0 and 1. Defaults to 1.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{EarthEngineMap} object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
library(sf)
library(rgee)
library(rgeeExtra)

ee_Initialize()

Map <- R6Map$new()

nc <- st_read(system.file("shape/nc.shp", package = "sf")) \%>\%
  st_transform(4326) \%>\%
  sf_as_ee()

ee_s2 <- ee$ImageCollection("COPERNICUS/S2")$
  filterDate("2016-01-01", "2016-01-31")$
  filterBounds(nc) \%>\%
  ee_get(0:2)

Map$centerObject(nc$geometry())
Map$addLayers(eeObject = ee_s2,position = "right")

# digging up the metadata
Map$previous_map_right$rgee$tokens

Map$reset()
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-addLegend"></a>}}
\if{latex}{\out{\hypertarget{method-addLegend}{}}}
\subsection{Method \code{addLegend()}}{
Adds a color legend to an EarthEngineMap.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$addLegend(
  visParams,
  name = "Legend",
  position = c("bottomright", "topright", "bottomleft", "topleft"),
  color_mapping = "numeric",
  opacity = 1,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{visParams}}{List of parameters for visualization.}

\item{\code{name}}{The title of the legend.}

\item{\code{position}}{Character. The position of the legend. By default bottomright.}

\item{\code{color_mapping}}{Map data values (numeric or factor/character) to
colors according to a given palette. Use "numeric" ("discrete") for continuous
(categorical) data. For display characters use "character" and add to visParams
the element "values" containing the desired character names.}

\item{\code{opacity}}{The legend's opacity is represented as a number between 0
and 1. Defaults to 1.}

\item{\code{...}}{Extra legend creator arguments. See \link[leaflet]{addLegend}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{EarthEngineMap} object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
library(leaflet)
library(rgee)
ee_Initialize()

Map$reset()

# Load MODIS ImageCollection
imgcol <- ee$ImageCollection$Dataset$MODIS_006_MOD13Q1

# Parameters for visualization
labels <- c("good", "marginal", "snow", "cloud")
cols   <- c("#999999", "#00BFC4", "#F8766D", "#C77CFF")
vis_qc <- list(min = 0, max = 3, palette = cols, bands = "SummaryQA", values = labels)

# Create interactive map
m_qc <- Map$addLayer(imgcol$median(), vis_qc, "QC")

# continous palette
Map$addLegend(vis_qc)

# categorical palette
Map$addLegend(vis_qc, name = "Legend1", color_mapping = "discrete")

# character palette
Map$addLegend(vis_qc, name = "Legend2", color_mapping = "character")
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{R6Map$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-upload.R
\name{ee_utils_create_json}
\alias{ee_utils_create_json}
\title{Convert an R list into a JSON file in the temp() file}
\usage{
ee_utils_create_json(x)
}
\arguments{
\item{x}{List to convert into a JSON file.}
}
\value{
A JSON file saved in a /tmp dir.
}
\description{
Convert an R list into a JSON file in the temp() file
}
\examples{
\dontrun{
library(rgee)
ee_utils_create_json(list(a=10,b=10))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_manage.R
\name{ee_manage-tools}
\alias{ee_manage-tools}
\alias{ee_manage_create}
\alias{ee_manage_delete}
\alias{ee_manage_assetlist}
\alias{ee_manage_quota}
\alias{ee_manage_copy}
\alias{ee_manage_move}
\alias{ee_manage_set_properties}
\alias{ee_manage_delete_properties}
\alias{ee_manage_asset_access}
\alias{ee_manage_task}
\alias{ee_manage_cancel_all_running_task}
\alias{ee_manage_asset_size}
\title{Interface to manage the Earth Engine Asset}
\usage{
ee_manage_create(path_asset, asset_type = "Folder", quiet = FALSE)

ee_manage_delete(path_asset, quiet = FALSE, strict = TRUE)

ee_manage_assetlist(path_asset, quiet = FALSE, strict = TRUE)

ee_manage_quota(quiet = FALSE)

ee_manage_copy(path_asset, final_path, strict = TRUE, quiet = FALSE)

ee_manage_move(path_asset, final_path, strict = TRUE, quiet = FALSE)

ee_manage_set_properties(path_asset, add_properties, strict = TRUE)

ee_manage_delete_properties(path_asset, del_properties = "ALL", strict = TRUE)

ee_manage_asset_access(
  path_asset,
  owner = NULL,
  editor = NULL,
  viewer = NULL,
  all_users_can_read = TRUE,
  quiet = FALSE
)

ee_manage_task(cache = FALSE)

ee_manage_cancel_all_running_task()

ee_manage_asset_size(path_asset, quiet = FALSE)
}
\arguments{
\item{path_asset}{Character. Name of the EE asset (Table, Image, Folder or
ImageCollection).}

\item{asset_type}{Character. The asset type to create ('Folder' or
'ImageCollection').}

\item{quiet}{Logical. Suppress info message.}

\item{strict}{Character vector. If TRUE, the existence of the asset will be
evaluated before performing the task.}

\item{final_path}{Character. Output filename
(e.g users/datacolecfbf/ic_moved)}

\item{add_properties}{List. Set of parameters to established as a property
of an EE object. See details.}

\item{del_properties}{Character. Names of properties to be deleted. See
details.}

\item{owner}{Character vector. Define owner user in the IAM Policy.}

\item{editor}{Character vector. Define editor users in the IAM Policy.}

\item{viewer}{Character vector. Define viewer users in the IAM Policy.}

\item{all_users_can_read}{Logical. All users can see the asset element.}

\item{cache}{Logical. If TRUE, the task report will be saved
in the /temp directory and used when the function.}
}
\description{
R functions to manage the Earth Engine Asset. The interface allows users
to create and eliminate folders, move and copy assets, set and delete
properties, handle access control lists, and manage and/or cancel tasks.
}
\details{
If the argument \code{del_properties} is 'ALL',
\link[=rgee]{ee_manage_delete_properties} will delete all
the properties.
}
\examples{
\dontrun{
library(rgee)

ee_Initialize()
ee_user_info()

# Change datacolecfbf by your EE user to be able to reproduce
user <- ee_get_assethome()
addm <- function(x) sprintf("\%s/\%s",user, x)
# 1. Create a folder or Image Collection
# Change path asset according to your specific user
ee_manage_create(addm("rgee"))

# 1. List all the elements inside a folder or a ImageCollection
ee_manage_assetlist(path_asset = addm("rgee"))

# 2. Create a Folder or a ImageCollection
ee_manage_create(
  path_asset = addm("rgee/rgee_folder"),
  asset_type = "Folder"
)

ee_manage_create(
  path_asset = addm("rgee/rgee_ic"),
  asset_type = "ImageCollection"
)

ee_manage_assetlist(path_asset = addm("rgee"))

# 3. Shows Earth Engine quota
ee_manage_quota()

# 4. Move an EE object to another folder
ee_manage_move(
  path_asset = addm("rgee/rgee_ic"),
  final_path = addm("rgee/rgee_folder/rgee_ic_moved")
)

ee_manage_assetlist(path_asset = addm("rgee/rgee_folder"))

# 5. Set properties to an EE object.
ee_manage_set_properties(
  path_asset = addm("rgee/rgee_folder/rgee_ic_moved"),
  add_properties = list(message = "hello-world", language = "R")
)

ic_id <- addm("rgee/rgee_folder/rgee_ic_moved")
test_ic <- ee$ImageCollection(ic_id)
test_ic$getInfo()

# 6. Delete properties
ee_manage_delete_properties(
  path_asset = addm("rgee/rgee_folder/rgee_ic_moved"),
  del_properties = c("message", "language")
)
test_ic$getInfo()

# 7. Create a report based on all the tasks
# that are running or have already been completed.
ee_manage_task()

# 8. Cancel all the running task
ee_manage_cancel_all_running_task()

# 9. Delete EE objects or folders
ee_manage_delete(addm("rgee/"))
}
}
\author{
Samapriya Roy, adapted to R and improved by csaybar.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_image.R
\name{ee_image_info}
\alias{ee_image_info}
\title{Approximate size of an EE Image object}
\usage{
ee_image_info(image, getsize = TRUE, compression_ratio = 20, quiet = FALSE)
}
\arguments{
\item{image}{Single-band EE Image object.}

\item{getsize}{Logical. If TRUE, the size of the object
is estimated.}

\item{compression_ratio}{Numeric. Measurement of the relative data size
reduction produced by a data compression algorithm (ignored if
\code{getsize} is FALSE). By default is 20.}

\item{quiet}{Logical. Suppress info message}
}
\value{
A list containing information about the number of rows (nrow),
number of columns (ncol), total number of pixels (total_pixel), and image
size (image_size).
}
\description{
Get the approximate number of rows, cols, and size of a single-band
Earth Engine Image.
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()

# World SRTM
srtm <- ee$Image("CGIAR/SRTM90_V4")
ee_image_info(srtm)

# Landast8
l8 <- ee$Image("LANDSAT/LC08/C01/T1_SR/LC08_038029_20180810")$select("B4")
ee_image_info(l8)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_clean.R
\name{ee_clean_pyenv}
\alias{ee_clean_pyenv}
\title{Remove rgee system variables from .Renviron}
\usage{
ee_clean_pyenv(Renviron = "global")
}
\arguments{
\item{Renviron}{Character. If it is "global" the environment variables in
the .Renviron located in the Sys.getenv("HOME") folder will be deleted. On the
other hand, if it is "local" the environment variables in the .Renviron on the
working directory (getwd()) will be deleted. Finally, users can also enter a
specific absolute path (see examples).}
}
\value{
No return value, called for cleaning environmental variables in
their system.
}
\description{
Remove rgee system variables from .Renviron
}
\seealso{
Other ee_clean functions: 
\code{\link{ee_clean_container}()},
\code{\link{ee_clean_credentials}()}
}
\concept{ee_clean functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.ee.computedobject.ComputedObject}
\alias{print.ee.computedobject.ComputedObject}
\alias{print}
\title{print Earth Engine object}
\usage{
\method{print}{ee.computedobject.ComputedObject}(x, ..., type = getOption("rgee.print.option"))
}
\arguments{
\item{x}{Earth Engine spatial object.}

\item{...}{ignored}

\item{type}{Character. What to show about the x object?. Three options are
supported: "json", "simply", "ee_print". By default "simply".}
}
\value{
No return value, called for displaying Earth Engine objects.
}
\description{
print Earth Engine object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_utils_py_to_r}
\alias{ee_utils_py_to_r}
\title{Convert between Python and R objects}
\usage{
ee_utils_py_to_r(x)
}
\arguments{
\item{x}{A python object}
}
\value{
An R object
}
\description{
Convert between Python and R objects
}
\seealso{
Other ee_utils functions: 
\code{\link{ee_utils_pyfunc}()},
\code{\link{ee_utils_shp_to_zip}()}
}
\concept{ee_utils functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_download.R
\name{ee_table_to_drive}
\alias{ee_table_to_drive}
\title{Creates a task to export a FeatureCollection to Google Drive.}
\usage{
ee_table_to_drive(
  collection,
  description = "myExportTableTask",
  folder = "rgee_backup",
  fileNamePrefix = NULL,
  timePrefix = TRUE,
  fileFormat = NULL,
  selectors = NULL
)
}
\arguments{
\item{collection}{The feature collection to be exported.}

\item{description}{Human-readable name of the task.}

\item{folder}{The name of a unique folder in your Drive
account to export into. Defaults to the root of the drive.}

\item{fileNamePrefix}{The Google Drive filename for the
export. Defaults to the name of the task.}

\item{timePrefix}{Add current date and time as a prefix to files to export.}

\item{fileFormat}{The output format: "CSV" (default), "GeoJSON",
"KML", "KMZ", "SHP", or "TFRecord".}

\item{selectors}{The list of properties to include in the output,
as a list of strings or a comma-separated string. By default, all
properties are included. **kwargs: Holds other keyword arguments
that may have been deprecated such as 'driveFolder' and
'driveFileNamePrefix'.}
}
\value{
An unstarted Task that exports the table to Google Drive.
}
\description{
Creates a task to export a FeatureCollection to Google Drive.
This function is a wrapper around \code{ee$batch$Export$table$toDrive(...)}.
}
\examples{
\dontrun{
library(rgee)
library(stars)
library(sf)

ee_users()
ee_Initialize(drive = TRUE)


# Define study area (local -> earth engine)
# Communal Reserve Amarakaeri - Peru
rlist <- list(xmin = -71.13, xmax = -70.95,ymin = -12.89, ymax = -12.73)
ROI <- c(rlist$xmin, rlist$ymin,
         rlist$xmax, rlist$ymin,
         rlist$xmax, rlist$ymax,
         rlist$xmin, rlist$ymax,
         rlist$xmin, rlist$ymin)
ee_ROI <- matrix(ROI, ncol = 2, byrow = TRUE) \%>\%
  list() \%>\%
  st_polygon() \%>\%
  st_sfc() \%>\%
  st_set_crs(4326) \%>\%
  sf_as_ee()

amk_fc <- ee$FeatureCollection(
  list(ee$Feature(ee_ROI, list(name = "Amarakaeri")))
)

task_vector <- ee_table_to_drive(
  collection = amk_fc,
  fileFormat = "GEO_JSON",
  fileNamePrefix = "geom_Amarakaeri"
)
task_vector$start()
ee_monitoring(task_vector) # optional
ee_drive_to_local(task = task_vector)
}
}
\seealso{
Other vector export task creator: 
\code{\link{ee_table_to_asset}()},
\code{\link{ee_table_to_gcs}()}
}
\concept{vector export task creator}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_utils_pyfunc}
\alias{ee_utils_pyfunc}
\title{Wrap an R function in a Python function with the same signature.}
\usage{
ee_utils_pyfunc(f)
}
\arguments{
\item{f}{An R function}
}
\value{
A Python function that calls the R function \code{f} with the same
signature.
}
\description{
This function could wrap an R function in a Python
function with the same signature. Note that the signature of the
R function must not contain esoteric Python-incompatible constructs.
}
\note{
\code{\link[reticulate]{py_func}} has been renamed to ee_utils_pyfunc
just to maintain the rgee functions name's style. All recognition
for this function must always be given to \pkg{reticulate}.
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()

# Earth Engine List
ee_SimpleList <- ee$List$sequence(0, 12)
ee_NewList <- ee_SimpleList$map(
  ee_utils_pyfunc(
    function(x) {
      ee$Number(x)$add(x)
    }
  )
)

ee_NewList$getInfo()

# Earth Engine ImageCollection
constant1 <- ee$Image(1)
constant2 <- ee$Image(2)
ee_ic <- ee$ImageCollection(c(constant2, constant1))
ee_newic <- ee_ic$map(
  ee_utils_pyfunc(
    function(x) ee$Image(x)$add(x)
  )
)
ee_newic$mean()$getInfo()$type
}
}
\seealso{
Other ee_utils functions: 
\code{\link{ee_utils_py_to_r}()},
\code{\link{ee_utils_shp_to_zip}()}
}
\author{
Yuan Tang and J.J. Allaire
}
\concept{ee_utils functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_download.R
\name{ee_drive_to_local}
\alias{ee_drive_to_local}
\title{Move results from Google Drive to a local directory}
\usage{
ee_drive_to_local(
  task,
  dsn,
  overwrite = TRUE,
  consider = TRUE,
  public = FALSE,
  metadata = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{task}{List generated after finished an EE task correctly. See details.}

\item{dsn}{Character. Output filename. If missing, a temporary
file will be assigned.}

\item{overwrite}{A boolean argument that indicates
whether "filename" should be overwritten. By default TRUE.}

\item{consider}{Interactive. See details.}

\item{public}{Logical. If TRUE, a public link to the Google Drive resource is
created.}

\item{metadata}{Logical. If TRUE, export the metadata related to the Google
Drive resource. See details.}

\item{quiet}{logical. Suppress info message.}
}
\value{
If \code{metadata} is FALSE, will return the filename of the Google
Drive resource on their system. Otherwise, a list with two elements
(\code{dns} and \code{metadata}) is returned.
}
\description{
Move results of an EE task saved in Google Drive to a local directory.
}
\details{
The task argument needs a status as task "COMPLETED" to work, due that the
parameters necessary to identify EE objects into Google Drive are obtained
from \cr \code{ee$batch$Export$*$toDrive(...)$start()$status()}.

\code{consider} argument is necessary due that Google Drive permits users to
create files with the same name. \code{consider} uses an interactive R
session by default to help users identify just the files that they want to
download. Additionally, the options "last" and "all" are implemented. "last"
will download just the last file saved in Google Drive while with "all" all
files will be downloaded.

Finally, if the argument \code{metadata} is TRUE, a list with the
following elements are exported join with the output filename (dsn):

\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{drive_name: }}{Name of the Table in Google Drive.}
\item{\bold{drive_id: }}{Id of the Table in Google Drive.}
\item{\bold{drive_download_link: }}{Download link to the table.}
}
}
\examples{
\dontrun{
library(rgee)
library(stars)
library(sf)

ee_users()
ee_Initialize(drive = TRUE)

# Define study area (local -> earth engine)
# Communal Reserve Amarakaeri - Peru
rlist <- list(xmin = -71.13, xmax = -70.95,ymin = -12.89, ymax = -12.73)
ROI <- c(rlist$xmin, rlist$ymin,
         rlist$xmax, rlist$ymin,
         rlist$xmax, rlist$ymax,
         rlist$xmin, rlist$ymax,
         rlist$xmin, rlist$ymin)

ee_ROI <- matrix(ROI, ncol = 2, byrow = TRUE) \%>\%
  list() \%>\%
  st_polygon() \%>\%
  st_sfc() \%>\%
  st_set_crs(4326) \%>\%
  sf_as_ee()


# Get the mean annual NDVI for 2011
cloudMaskL457 <- function(image) {
  qa <- image$select("pixel_qa")
  cloud <- qa$bitwiseAnd(32L)$
    And(qa$bitwiseAnd(128L))$
    Or(qa$bitwiseAnd(8L))
  mask2 <- image$mask()$reduce(ee$Reducer$min())
  image <- image$updateMask(cloud$Not())$updateMask(mask2)
  image$normalizedDifference(list("B4", "B3"))
}

ic_l5 <- ee$ImageCollection("LANDSAT/LT05/C01/T1_SR")$
  filterBounds(ee$FeatureCollection(ee_ROI))$
  filterDate("2011-01-01", "2011-12-31")$
  map(cloudMaskL457)

# Create simple composite
mean_l5 <- ic_l5$mean()$rename("NDVI")
mean_l5 <- mean_l5$reproject(crs = "EPSG:4326", scale = 500)
mean_l5_Amarakaeri <- mean_l5$clip(ee_ROI)

# Move results from Earth Engine to Drive
task_img <- ee_image_to_drive(
  image = mean_l5_Amarakaeri,
  folder = "Amarakaeri",
  fileFormat = "GEO_TIFF",
  region = ee_ROI,
  fileNamePrefix = "my_image_demo"
)

task_img$start()
ee_monitoring(task_img)

# Move results from Drive to local
img <- ee_drive_to_local(task = task_img)
}
}
\seealso{
Other generic download functions: 
\code{\link{ee_gcs_to_local}()}
}
\concept{generic download functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_Initialize.R
\name{ee_user_info}
\alias{ee_user_info}
\title{Display the credentials and general info of the initialized user}
\usage{
ee_user_info(quiet = FALSE)
}
\arguments{
\item{quiet}{Logical. Suppress info messages.}
}
\value{
A list with information about the Earth Engine user.
}
\description{
Display the credentials and general info of the initialized user
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()
ee_user_info()
}
}
\seealso{
Other session management functions: 
\code{\link{ee_Initialize}()},
\code{\link{ee_users}()},
\code{\link{ee_version}()}
}
\concept{session management functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_utils_future_value}
\alias{ee_utils_future_value}
\title{The value of a future or the values of all elements in a container}
\usage{
ee_utils_future_value(future, stdout = TRUE, signal = TRUE, ...)
}
\arguments{
\item{future, }{x A Future, an environment, a list, or a list environment.}

\item{stdout}{If TRUE, standard output captured while resolving futures
is relayed, otherwise not.}

\item{signal}{If TRUE, \link[base]{conditions} captured while resolving
futures are relayed, otherwise not.}

\item{\dots}{All arguments used by the S3 methods.}
}
\value{
\code{value()} of a Future object returns the value of the future, which can
be any type of \R object.

\code{value()} of a list, an environment, or a list environment returns an
object with the same number of elements and of the same class.
Names and dimension attributes are preserved, if available.
All future elements are replaced by their corresponding \code{value()} values.
For all other elements, the existing object is kept as-is.

If \code{signal} is TRUE and one of the futures produces an error, then
that error is produced.
}
\description{
Gets the value of a future or the values of all elements (including futures)
in a container such as a list, an environment, or a list environment.
If one or more futures is unresolved, then this function blocks until all
queried futures are resolved.
}
\author{
Henrik Bengtsson \url{https://github.com/HenrikBengtsson/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raster_as_ee.R
\name{raster_as_ee}
\alias{raster_as_ee}
\title{Convert a Raster* object into an EE Image object}
\usage{
raster_as_ee(
  x,
  assetId,
  bucket = NULL,
  predefinedAcl = "bucketLevel",
  command_line_tool_path = NULL,
  overwrite = FALSE,
  monitoring = TRUE,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{x}{RasterLayer, RasterStack or RasterBrick object to be converted into
an ee$Image.}

\item{assetId}{Character. Destination asset ID for the uploaded file.}

\item{bucket}{Character. Name of the GCS bucket.}

\item{predefinedAcl}{Specify user access to object. Passed to
\code{googleCloudStorageR::gcs_upload}.}

\item{command_line_tool_path}{Character. Path to the Earth Engine command line
tool (CLT). If NULL, rgee assumes that CLT is set in the system PATH.
(ignore if \code{via} is not defined as "gcs_to_asset").}

\item{overwrite}{Logical. If TRUE, the assetId will be overwritten.}

\item{monitoring}{Logical. If TRUE the exportation task will be monitored.}

\item{quiet}{Logical. Suppress info message.}

\item{...}{parameter(s) passed on to \code{\link{ee_utils_create_manifest_image}}}
}
\value{
An ee$Image object
}
\description{
Convert a Raster* object into an EE Image object
}
\examples{
\dontrun{
library(raster)
library(stars)
library(rgee)

ee_Initialize(gcs = TRUE)

# Get the filename of a image
tif <- system.file("tif/L7_ETMs.tif", package = "stars")
x <- stack(tif)
assetId <- sprintf("\%s/\%s",ee_get_assethome(),'raster_l7')

# Method 1
# 1. Move from local to gcs
gs_uri <- local_to_gcs(x = tif, bucket = 'rgee_dev')

# 2. Create a manifest
manifest <- ee_utils_create_manifest_image(gs_uri, assetId)

# 3. Pass from gcs to asset
gcs_to_ee_image(
 manifest = manifest,
 overwrite = TRUE
)

# OPTIONAL: Monitoring progress
ee_monitoring()

# OPTIONAL: Display results
ee_stars_01 <- ee$Image(assetId)
Map$centerObject(ee_stars_01)
Map$addLayer(ee_stars_01, list(min = 0, max = 255))

# Method 2
ee_stars_02 <- raster_as_ee(
 x = x,
 overwrite = TRUE,
 assetId = assetId,
 bucket = "rgee_dev"
)
Map$centerObject(ee_stars_02)
Map$addLayer(ee_stars_02, list(min = 0, max = 255))
}
}
\seealso{
Other image upload functions: 
\code{\link{stars_as_ee}()}
}
\concept{image upload functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_image.R
\name{ee_as_stars}
\alias{ee_as_stars}
\title{Convert an Earth Engine (EE) image in a stars object}
\usage{
ee_as_stars(
  image,
  region = NULL,
  dsn = NULL,
  via = "drive",
  container = "rgee_backup",
  scale = NULL,
  maxPixels = 1e+09,
  lazy = FALSE,
  public = TRUE,
  add_metadata = TRUE,
  timePrefix = TRUE,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{image}{ee$Image to be converted into a stars object.}

\item{region}{EE Geometry (ee$Geometry$Polygon) which specifies the region
to export. CRS needs to be the same that the argument \code{image}.
Otherwise, it will be forced. If not specified, image bounds are taken.}

\item{dsn}{Character. Output filename. If missing, a temporary file is
created.}

\item{via}{Character. Method to export the image. Two methods are
implemented: "drive", "gcs". See details.}

\item{container}{Character. Name of the folder ('drive') or bucket ('gcs')
to be exported.}

\item{scale}{Numeric. The resolution in meters per pixel. Defaults
to the native resolution of the image.}

\item{maxPixels}{Numeric. The maximum allowed number of pixels in the
exported image. The task will fail if the exported region covers
more pixels in the specified projection. Defaults to 100,000,000.}

\item{lazy}{Logical. If TRUE, a \code{\link[future:sequential]{
future::sequential}} object is created to evaluate the task in the future.
See details.}

\item{public}{Logical. If TRUE, a public link to the image is created.}

\item{add_metadata}{Add metadata to the stars_proxy object. See details.}

\item{timePrefix}{Logical. Add current date and time (\code{Sys.time()}) as
a prefix to export files. This parameter helps to avoid exported files
with the same name. By default TRUE.}

\item{quiet}{Logical. Suppress info message}

\item{...}{Extra exporting argument. See \link{ee_image_to_drive} and
\link{ee_image_to_gcs}.}
}
\value{
A stars-proxy object
}
\description{
Convert an ee$Image in a stars object.
}
\details{
\code{ee_as_stars} supports the download of \code{ee$Images}
by two different options: "drive"
(\href{https://CRAN.R-project.org/package=googledrive}{Google Drive}) and "gcs"
(\href{https://CRAN.R-project.org/package=googleCloudStorageR}{
Google Cloud Storage}). In both cases, \code{ee_as_stars} works as follow:
\itemize{
\item{1. }{A task is started (i.e. \code{ee$batch$Task$start()}) to
move the \code{ee$Image} from Earth Engine to the intermediate container
specified in the argument \code{via}.}
\item{2. }{If the argument \code{lazy} is TRUE, the task will not be
monitored. This is useful to lunch several tasks simultaneously and
calls them later using \code{\link{ee_utils_future_value}} or
\code{\link[future:value]{future::value}}. At the end of this step,
the \code{ee$Image} is stored on the path specified in the argument
\code{dsn}.}
\item{3. }{Finally, if the argument \code{add_metadata} is TRUE, a list
with the following elements is added to the stars-proxy object.
\itemize{
\item{\bold{if via is "drive":}}
\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{drive_name: }}{Name of the Image in Google Drive.}
\item{\bold{drive_id: }}{Id of the Image in Google Drive.}
\item{\bold{drive_download_link: }}{Download link to the image.}
}
}
\itemize{
\item{\bold{if via is "gcs":}}
\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{gcs_name: }}{Name of the Image in Google Cloud Storage.}
\item{\bold{gcs_bucket: }}{Name of the bucket.}
\item{\bold{gcs_fileFormat: }}{Format of the image.}
\item{\bold{gcs_public_link: }}{Download link to the image.}
\item{\bold{gcs_URI: }}{gs:// link to the image.}
}
}
Run \code{attr(stars, "metadata")} to get the list.
}
}

For getting more information about exporting data from Earth Engine, take
a look at the
\href{https://developers.google.com/earth-engine/guides/exporting}{Google
Earth Engine Guide - Export data}.
}
\examples{
\dontrun{
library(rgee)

ee_Initialize(drive = TRUE, gcs = TRUE)
ee_user_info()

# Define an image.
img <- ee$Image("LANDSAT/LC08/C01/T1_SR/LC08_038029_20180810")$
  select(c("B4", "B3", "B2"))$
  divide(10000)

# OPTIONAL display it using Map
Map$centerObject(eeObject = img)
Map$addLayer(eeObject = img, visParams = list(max = 0.4,gamma=0.1))

# Define an area of interest.
geometry <- ee$Geometry$Rectangle(
  coords = c(-110.8, 44.6, -110.6, 44.7),
  proj = "EPSG:4326",
  geodesic = FALSE
)

## drive - Method 01
# Simple
img_02 <- ee_as_stars(
  image = img,
  region = geometry,
  via = "drive"
)

# Lazy
img_02 <- ee_as_stars(
  image = img,
  region = geometry,
  via = "drive",
  lazy = TRUE
)

img_02_result <- img_02 \%>\% ee_utils_future_value()
attr(img_02_result, "metadata") # metadata

## gcs - Method 02
# Simple
img_03 <- ee_as_stars(
  image = img,
  region = geometry,
  container = "rgee_dev",
  via = "gcs"
)

# Lazy
img_03 <- ee_as_stars(
  image = img,
  region = geometry,
  container = "rgee_dev",
  lazy = TRUE,
  via = "gcs"
)

img_03_result <- img_03 \%>\% ee_utils_future_value()
attr(img_03_result, "metadata") # metadata

# OPTIONAL: clean containers
# ee_clean_container(name = "rgee_backup", type = "drive")
# ee_clean_container(name = "rgee_dev", type = "gcs")
}
}
\seealso{
Other image download functions: 
\code{\link{ee_as_raster}()},
\code{\link{ee_as_thumbnail}()},
\code{\link{ee_imagecollection_to_local}()}
}
\concept{image download functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_download.R
\name{ee_image_to_asset}
\alias{ee_image_to_asset}
\title{Creates a task to export an EE Image to their EE Assets.}
\usage{
ee_image_to_asset(
  image,
  description = "myExportImageTask",
  assetId = NULL,
  overwrite = FALSE,
  pyramidingPolicy = NULL,
  dimensions = NULL,
  region = NULL,
  scale = NULL,
  crs = NULL,
  crsTransform = NULL,
  maxPixels = NULL
)
}
\arguments{
\item{image}{The image to be exported.}

\item{description}{Human-readable name of the task.}

\item{assetId}{The destination asset ID.}

\item{overwrite}{Logical. If TRUE, the assetId will be overwritten if
it exists.}

\item{pyramidingPolicy}{The pyramiding policy to apply to each band
in the image, a dictionary keyed by band name. Values must be one
of: "mean", "sample", "min", "max", or "mode". Defaults to "mean".
A special key, ".default", may be used to change the default for all bands.}

\item{dimensions}{The dimensions of the exported image. It takes either a
single positive integer as the maximum dimension or "WIDTHxHEIGHT" where
WIDTH and HEIGHT are each positive integers.}

\item{region}{The lon,lat coordinates for a LinearRing or Polygon
specifying the region to export. It can be specified as nested lists
of numbers or a serialized string. Defaults to the image's region.}

\item{scale}{The resolution in meters per pixel. Defaults to the native
resolution of the image asset unless a crsTransform is specified.}

\item{crs}{The coordinate reference system of the exported image's
projection. Defaults to the image's default projection.}

\item{crsTransform}{A comma-separated string of 6 numbers describing
the affine transform of the coordinate reference system of the exported
image's projection, in the order:
xScale, xShearing, xTranslation, yShearing, yScale, and yTranslation.
Defaults to the image's native CRS transform.}

\item{maxPixels}{The maximum allowed number of pixels in the exported
image. The task will fail if the exported region covers more pixels
in the specified projection. Defaults to 100,000,000. **kwargs: Holds
other keyword arguments that may have been deprecated, such
as 'crs_transform'.}
}
\value{
An unstarted task
}
\description{
Creates a task to export an EE Image to their EE Assets.
This function is a wrapper around \code{ee$batch$Export$image$toAsset(...)}.
}
\examples{
\dontrun{
library(rgee)
library(stars)
library(sf)

ee_users()
ee_Initialize()

# Define study area (local -> earth engine)
# Communal Reserve Amarakaeri - Peru
rlist <- list(xmin = -71.13, xmax = -70.95,ymin = -12.89, ymax = -12.73)
ROI <- c(rlist$xmin, rlist$ymin,
         rlist$xmax, rlist$ymin,
         rlist$xmax, rlist$ymax,
         rlist$xmin, rlist$ymax,
         rlist$xmin, rlist$ymin)
ee_ROI <- matrix(ROI, ncol = 2, byrow = TRUE) \%>\%
  list() \%>\%
  st_polygon() \%>\%
  st_sfc() \%>\%
  st_set_crs(4326) \%>\%
  sf_as_ee()


# Get the mean annual NDVI for 2011
cloudMaskL457 <- function(image) {
  qa <- image$select("pixel_qa")
  cloud <- qa$bitwiseAnd(32L)$
    And(qa$bitwiseAnd(128L))$
    Or(qa$bitwiseAnd(8L))
  mask2 <- image$mask()$reduce(ee$Reducer$min())
  image <- image$updateMask(cloud$Not())$updateMask(mask2)
  image$normalizedDifference(list("B4", "B3"))
}

ic_l5 <- ee$ImageCollection("LANDSAT/LT05/C01/T1_SR")$
  filterBounds(ee$FeatureCollection(ee_ROI))$
  filterDate("2011-01-01", "2011-12-31")$
  map(cloudMaskL457)

# Create simple composite
mean_l5 <- ic_l5$mean()$rename("NDVI")
mean_l5 <- mean_l5$reproject(crs = "EPSG:4326", scale = 500)
mean_l5_Amarakaeri <- mean_l5$clip(ee_ROI)

# Move results from Earth Engine to Drive
assetid <- paste0(ee_get_assethome(), '/l5_Amarakaeri')
task_img <- ee_image_to_asset(
  image = mean_l5_Amarakaeri,
  assetId = assetid,
  overwrite = TRUE,
  scale = 500,
  region = ee_ROI
)

task_img$start()
ee_monitoring(task_img)

ee_l5 <- ee$Image(assetid)
Map$centerObject(ee_l5)
Map$addLayer(ee_l5)
}
}
\seealso{
Other image export task creator: 
\code{\link{ee_image_to_drive}()},
\code{\link{ee_image_to_gcs}()}
}
\concept{image export task creator}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Map_operators.R
\name{map-operator}
\alias{map-operator}
\alias{+.EarthEngineMap}
\alias{|.EarthEngineMap}
\alias{|,}
\alias{EarthEngineMap,}
\alias{EarthEngineMap-method}
\title{EarthEngineMap + EarthEngineMap; adds data from the second map to the first}
\usage{
\method{+}{EarthEngineMap}(e1, e2)

\method{|}{EarthEngineMap}(e1, e2)
}
\arguments{
\item{e1}{an EarthEngineMap object.}

\item{e2}{an EarthEngineMap object.}
}
\description{
EarthEngineMap + EarthEngineMap; adds data from the second map to the first

EarthEngineMap | EarthEngineMap provides a slider in the middle to compare two maps.
}
\author{
tim-salabim. Adapted from mapview code.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-upload.R
\name{ee_utils_get_crs}
\alias{ee_utils_get_crs}
\title{Convert EPSG, ESRI or SR-ORG code into a OGC WKT}
\usage{
ee_utils_get_crs(code)
}
\arguments{
\item{code}{The projection code.}
}
\value{
A character which represents the same projection in WKT2 string.
}
\description{
Convert EPSG, ESRI or SR-ORG code into a OGC WKT
}
\examples{
\dontrun{
library(rgee)

ee_utils_get_crs("SR-ORG:6864")
ee_utils_get_crs("EPSG:4326")
ee_utils_get_crs("ESRI:37002")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_print.R
\name{ee_print}
\alias{ee_print}
\alias{ee_print.ee.geometry.Geometry}
\alias{ee_print.ee.feature.Feature}
\alias{ee_print.ee.featurecollection.FeatureCollection}
\alias{ee_print.ee.image.Image}
\alias{ee_print.ee.imagecollection.ImageCollection}
\title{Print and return metadata about Spatial Earth Engine Objects}
\usage{
ee_print(eeobject, ...)

\method{ee_print}{ee.geometry.Geometry}(eeobject, ..., clean = FALSE, quiet = FALSE)

\method{ee_print}{ee.feature.Feature}(eeobject, ..., clean = FALSE, quiet = FALSE)

\method{ee_print}{ee.featurecollection.FeatureCollection}(eeobject, ..., f_index = 0, clean = FALSE, quiet = FALSE)

\method{ee_print}{ee.image.Image}(
  eeobject,
  ...,
  img_band,
  time_end = TRUE,
  compression_ratio = 20,
  clean = FALSE,
  quiet = FALSE
)

\method{ee_print}{ee.imagecollection.ImageCollection}(
  eeobject,
  ...,
  time_end = TRUE,
  img_index = 0,
  img_band,
  compression_ratio = 20,
  clean = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{eeobject}{Earth Engine Object. Available for: Geometry, Feature,
FeatureCollection, Image or ImageCollection.}

\item{...}{ignored}

\item{clean}{Logical. If TRUE, the cache will be cleaned.}

\item{quiet}{Logical. Suppress info message}

\item{f_index}{Numeric. Index of the \code{ee$FeatureCollection} to fetch.
Relevant just for \code{ee$FeatureCollection} objects.}

\item{img_band}{Character. Band name of the \code{ee$Image} to fetch.
Relevant just for \code{ee$ImageCollection} and \code{ee$Image} objects.}

\item{time_end}{Logical. If TRUE, the system:time_end property in ee$Image
is also returned. See \code{rgee::ee_get_date_img} for details.}

\item{compression_ratio}{Numeric. Measurement of the relative data size
reduction produced by a data compression algorithm (ignored if \code{eeobject}
is not an \code{ee$Image} or \code{ee$ImageCollection}). By default is 20.}

\item{img_index}{Numeric. Index of the \code{ee$ImageCollection} to fetch.
Relevant just for \code{ee$ImageCollection} objects.}
}
\value{
A list with the metadata of the Earth Engine object.
}
\description{
Print and return metadata about Spatial Earth Engine Objects.
\code{ee_print} can retrieve information about the number of images
or features, number of bands or geometries, number of pixels, geotransform,
data type, properties, and object size.
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()

# Geometry
geom <- ee$Geometry$Rectangle(-10,-10,10,10)
Map$addLayer(geom)
ee_print(geom)

# Feature
feature <- ee$Feature(geom, list(rgee = "ee_print", data = TRUE))
ee_print(feature)

# FeatureCollection
featurecollection <- ee$FeatureCollection(feature)
ee_print(featurecollection)

# Image
srtm <- ee$Image("CGIAR/SRTM90_V4")
ee_print(srtm)

srtm_clip <- ee$Image("CGIAR/SRTM90_V4")$clip(geom)
srtm_metadata <- ee_print(srtm_clip)
srtm_metadata$img_bands_names

# ImageCollection
object <- ee$ImageCollection("LANDSAT/LC08/C01/T1_TOA")$
  filter(ee$Filter()$eq("WRS_PATH", 44))$
  filter(ee$Filter()$eq("WRS_ROW", 34))$
  filterDate("2014-03-01", "2014-08-01")$
  aside(ee_print)
}
}
\seealso{
Other helper functions: 
\code{\link{ee_help}()},
\code{\link{ee_monitoring}()}
}
\concept{helper functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_clean.R
\name{ee_clean_container}
\alias{ee_clean_container}
\title{Delete files from either a Folder (Google Drive) or a Bucket (GCS)}
\usage{
ee_clean_container(name = "rgee_backup", type = "drive", quiet = FALSE)
}
\arguments{
\item{name}{Character. Name of the folder (Google Drive) or bucket (GCS)
to delete all files.}

\item{type}{Character. Name of the file storage web service. 'drive'
and 'gcs' are supported.}

\item{quiet}{logical. Suppress info message}
}
\value{
No return value, called for cleaning Google Drive or Google
Cloud Storage container.
}
\description{
Delete all files from a folder (Google Drive) or a bucket
(Google Cloud Storage). Caution: This will permanently delete
their backup files generated by using \code{ee_as_stars} and \code{ee_as_sf}.
}
\seealso{
Other ee_clean functions: 
\code{\link{ee_clean_credentials}()},
\code{\link{ee_clean_pyenv}()}
}
\concept{ee_clean functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_Initialize.R
\name{ee_Initialize}
\alias{ee_Initialize}
\title{Authenticate and Initialize Earth Engine}
\usage{
ee_Initialize(
  user = NULL,
  drive = FALSE,
  gcs = FALSE,
  display = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{user}{Character (optional, e.g. \code{data.colec.fbf}). The user
argument is used to create a folder inside the path
\code{~/.config/earthengine/} that save all the credentials for a specific
Google identity.}

\item{drive}{Logical (optional). If TRUE, the drive credential
is cached in the path \code{~/.config/earthengine/}.}

\item{gcs}{Logical (optional). If TRUE, the Google Cloud Storage
credential is cached in the path \code{~/.config/earthengine/}.}

\item{display}{Logical. If TRUE, display the earthengine authentication URL.}

\item{quiet}{Logical. Suppress info messages.}
}
\value{
No return value, called for initializing the earthengine-api.
}
\description{
Authorize rgee to manage Earth Engine resources, Google
Drive, and Google Cloud Storage. The \code{ee_initialize()} via
web-browser will ask users to sign into your Google account and
allows you to grant permission to manage resources. This function is
a wrapper around \code{rgee::ee$Initialize()}.
}
\details{
\code{ee_Initialize(...)} can manage Google Drive, and Google
Cloud Storage resources using the R packages googledrive and
googlecloudStorageR, respectively. By default, rgee does not require
them. These are only necessary to enable rgee I/O functionality.
All user credentials are saved in the directory
\code{~/.config/earthengine/}. If a user does not specify the "user"
argument, all user credentials are saved in the the subdirectory
\code{~/.config/earthengine/ndef}.
}
\examples{
\dontrun{
library(rgee)

# Simple init - Load just the Earth Engine credential
ee_Initialize()
ee_user_info()
}
}
\seealso{
Other session management functions: 
\code{\link{ee_user_info}()},
\code{\link{ee_users}()},
\code{\link{ee_version}()}
}
\concept{session management functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_download.R
\name{ee_table_to_gcs}
\alias{ee_table_to_gcs}
\title{Creates a task to export a FeatureCollection to Google Cloud Storage.}
\usage{
ee_table_to_gcs(
  collection,
  description = "myExportTableTask",
  bucket = NULL,
  fileNamePrefix = NULL,
  timePrefix = TRUE,
  fileFormat = NULL,
  selectors = NULL
)
}
\arguments{
\item{collection}{The feature collection to be exported.}

\item{description}{Human-readable name of the task.}

\item{bucket}{The name of a Cloud Storage bucket for the export.}

\item{fileNamePrefix}{Cloud Storage object name prefix
for the export. Defaults to the name of the task.}

\item{timePrefix}{Add current date and time as a prefix to files to export.}

\item{fileFormat}{The output format: "CSV" (default),
"GeoJSON", "KML", "KMZ", "SHP", or "TFRecord".}

\item{selectors}{The list of properties to include in the output,
as a list of strings or a comma-separated string. By default, all
properties are included. **kwargs: Holds other keyword arguments
that may have been deprecated such as 'outputBucket'.}
}
\value{
An unstarted Task that exports the table to Google Cloud Storage.
}
\description{
Creates a task to export a FeatureCollection to Google Cloud Storage.
This function is a wrapper around
\code{ee$batch$Export$table$toCloudStorage(...)}.
}
\examples{
\dontrun{
library(rgee)
library(stars)
library(sf)

ee_users()
ee_Initialize(gcs = TRUE)

# Define study area (local -> earth engine)
# Communal Reserve Amarakaeri - Peru
rlist <- list(xmin = -71.13, xmax = -70.95,ymin = -12.89, ymax = -12.73)
ROI <- c(rlist$xmin, rlist$ymin,
         rlist$xmax, rlist$ymin,
         rlist$xmax, rlist$ymax,
         rlist$xmin, rlist$ymax,
         rlist$xmin, rlist$ymin)
ee_ROI <- matrix(ROI, ncol = 2, byrow = TRUE) \%>\%
  list() \%>\%
  st_polygon() \%>\%
  st_sfc() \%>\%
  st_set_crs(4326) \%>\%
  sf_as_ee()

amk_fc <- ee$FeatureCollection(
  list(ee$Feature(ee_ROI, list(name = "Amarakaeri")))
)

task_vector <- ee_table_to_gcs(
  collection = amk_fc,
  bucket = "rgee_dev",
  fileFormat = "SHP",
  fileNamePrefix = "geom_Amarakaeri"
)
task_vector$start()
ee_monitoring(task_vector) # optional
amk_geom <- ee_gcs_to_local(task = task_vector)
plot(sf::read_sf(amk_geom[3]), border = "red", lwd = 10)
}
}
\seealso{
Other vector export task creator: 
\code{\link{ee_table_to_asset}()},
\code{\link{ee_table_to_drive}()}
}
\concept{vector export task creator}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgee-package.R
\docType{package}
\name{rgee-package}
\alias{rgee}
\alias{rgee-package}
\title{rgee: An R package for interacting with Google Earth Engine}
\description{
Google Earth Engine (Gorelick et al., 2017) is a cloud computing platform
designed for planetary-scale environmental data analysis that only can be
accessed via the Earth Engine code editor, third-party web apps, and the
JavaScript and Python client libraries. \code{rgee} is a non-official
client library for R that uses \code{reticulate} to wrap the Earth Engine
Python API and provide R users with a familiar interface, rapid development
features, and flexibility to analyze data using open-source, R third-party
packages.
}
\details{
The package implements and supports:

\itemize{
\item Earth Engine Module
\item Install or set all rgee dependencies
\item Check non-R dependencies
\item Clean non-R dependencies
\item Session management
\item Transform an R Date to an EE Date or vice versa
\item Create Interactive visualization Maps
\item Image download
\item Vector download
\item Generic download
\item Assets management
\item Upload raster
\item Upload vector
\item Upload generic
\item Extract values
\item Helper functions
\item Utils functions
}
}
\section{I. Earth Engine Module}{


Interface to main Earth Engine module. Provides access to top level classes
and functions as well as sub-modules (e.g. ee$Image, ee$FeatureCollection$first, etc.).

\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee}}\tab Main Earth Engine module. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{II. Install or set non-R rgee dependencies}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_install}}\tab Create an isolated Python virtual environment with all rgee dependencies. \cr
\code{\link{ee_install_set_pyenv}}\tab Configure which version of Python to use with rgee. \cr
\code{\link{ee_install_upgrade}}\tab Upgrade the Earth Engine Python API. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{III. Check non-R dependencies}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_check}}\tab Check all non-R dependencies. \cr
\code{\link{ee_check_python}}\tab Check Python environment. \cr
\code{\link{ee_check_credentials}}\tab Check Google credentials. \cr
\code{\link{ee_check_python_packages}}\tab Check Python packages: earthengine-api and numpy. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{IV. Clean container, credentials, or rgee system variables}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_clean_container}}\tab Delete files from either a Folder or a Bucket. \cr
\code{\link{ee_clean_credentials}}\tab Delete Credentials. \cr
\code{\link{ee_clean_pyenv}}\tab Remove rgee system variables from .Renviron. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{V. Session management}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_Initialize}}\tab Authenticate and Initialize Earth Engine. \cr
\code{\link{ee_version}}\tab Earth Engine API version. \cr
\code{\link{ee_user_info}}\tab Display the credentials and general info of the initialized user. \cr
\code{\link{ee_users}}\tab Display the credentials of all users as a table. \cr
\code{\link{ee_get_assethome}}\tab Get the Asset home name. \cr
\code{\link{ee_get_earthengine_path}}\tab Get the path where the credentials are stored. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{VII. Transform an R Date to an EE Date or vice versa}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{eedate_to_rdate}}\tab Pass an Earth Engine date object to R. \cr
\code{\link{rdate_to_eedate}}\tab Pass an R date object to Earth Engine. \cr
\code{\link{ee_get_date_img}}\tab Get the date of a EE Image. \cr
\code{\link{ee_get_date_ic}}\tab Get the date of a EE ImageCollection. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{VIII. Visualization Map}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{Map}}\tab R6 object (Map) to display Earth Engine (EE) spatial objects. \cr
\code{\link{R6Map}}\tab R6 class to display Earth Engine (EE) spatial objects. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{IX. Image download}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_as_raster}}\tab Convert an Earth Engine (EE) image in a raster object. \cr
\code{\link{ee_as_stars}}\tab Convert an Earth Engine (EE) image in a stars object. \cr
\code{\link{ee_as_thumbnail}}\tab Create an R spatial gridded object from an EE thumbnail image. \cr
\code{\link{ee_image_to_asset}}\tab Creates a task to export an EE Image to their EE Assets. \cr
\code{\link{ee_image_to_drive}}\tab Creates a task to export an EE Image to Drive. \cr
\code{\link{ee_image_to_gcs}}\tab Creates a task to export an EE Image to Google Cloud Storage. \cr
\code{\link{ee_image_info}}\tab Approximate size of an EE Image object. \cr
\code{\link{ee_imagecollection_to_local}}\tab Save an EE ImageCollection in their local system. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{X. Vector download}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_as_sf}}\tab Convert an Earth Engine table in an sf object. \cr
\code{\link{ee_table_to_asset}}\tab Creates a task to export a FeatureCollection to an EE table asset. \cr
\code{\link{ee_table_to_drive}}\tab Creates a task to export a FeatureCollection to Google Drive. \cr
\code{\link{ee_table_to_gcs}}\tab Creates a task to export a FeatureCollection to Google Cloud Storage. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{XI. Generic download}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_drive_to_local}}\tab Move results from Google Drive to a local directory. \cr
\code{\link{ee_gcs_to_local}}\tab Move results from Google Cloud Storage to a local directory. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{XII. Assets management}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_manage-tools}}\tab Interface to manage the Earth Engine Asset. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{XIII. Upload raster}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{stars_as_ee}}\tab Convert a stars or stars-proxy object into an EE Image object. \cr
\code{\link{raster_as_ee}}\tab Convert a Raster* object into an EE Image object. \cr
\code{\link{gcs_to_ee_image}}\tab Move a GeoTIFF image from GCS to their EE assets. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{XIV. Upload vector}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{gcs_to_ee_table}}\tab Move a zipped shapefile from GCS to their EE Assets. \cr
\code{\link{sf_as_ee}}\tab Convert an sf object to an EE object. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{XV. Upload generic}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{local_to_gcs}}\tab Upload local files to Google Cloud Storage. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{XVI. Extract values}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_extract}}\tab Extract values from EE Images or ImageCollections objects. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{XVII. Helper functions}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_help}}\tab Documentation for Earth Engine Objects. \cr
\code{\link{ee_print}}\tab Print and return metadata about Spatial Earth Engine Objects. \cr
\code{\link{ee_monitoring}}\tab Monitoring Earth Engine task progress. \cr
\code{\link{print}}\tab Print Earth Engine objects. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{XVIII. Utils functions}{


\tabular{ll}{
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
\code{\link{ee_utils_py_to_r}}\tab Convert between Python and R objects. \cr
\code{\link{ee_utils_pyfunc}}\tab Wrap an R function in a Python function with the same signature. \cr
\code{\link{ee_utils_shp_to_zip}}\tab Create a zip file from an sf object. \cr
\code{\link{ee_utils_create_json}}\tab Convert a R list into a JSON file. \cr
\code{\link{ee_utils_create_manifest_image}}\tab Create a manifest to upload an image. \cr
\code{\link{ee_utils_create_manifest_table}}\tab Create a manifest to upload a table. \cr
\code{\link{ee_utils_get_crs}}\tab Convert EPSG, ESRI or SR-ORG code into a OGC WKT. \cr
\code{\link{ee_utils_future_value}}\tab The value of a future or the values of all elements in a container. \cr
\code{\link{ee_utils_dataset_display}}\tab Search into the Earth Engine Data Catalog. \cr
---------------------------\tab --------------------------------------------------------------------------------------------------- \cr
}
}

\section{Acknowledgments}{


We want to offer a special thanks to Justin Braaten for his wise and
helpful comments in the whole development of rgee. As well, we would like
to mention the following third-party R/Python packages for contributing
indirectly to the improvement of rgee:

\itemize{
\item gee_asset_manager - Lukasz Tracewski
\item geeup - Samapriya Roy
\item geeadd - Samapriya Roy
\item eemont - David Montero Loaiza
\item cartoee - Kel Markert
\item geetools - Rodrigo E. Principe
\item landsat-extract-gee - Loïc Dutrieux
\item earthEngineGrabR - JesJehle
\item sf - Edzer Pebesma
\item stars - Edzer Pebesma
\item gdalcubes - Marius Appel
}
}

\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/r-spatial/rgee/}
  \item \url{https://r-spatial.github.io/rgee/}
  \item \url{https://github.com/google/earthengine-api/}
  \item Report bugs at \url{https://github.com/r-spatial/rgee/issues/}
}

}
\author{
\strong{Maintainer}: Cesar Aybar \email{csaybar@gmail.com} (\href{https://orcid.org/0000-0003-2745-9535}{ORCID})

Other contributors:
\itemize{
  \item Wu Qiusheng \email{giswqs@gmail.com} (\href{https://orcid.org/0000-0001-5437-4073}{ORCID}) [contributor]
  \item Lesly Bautista \email{leslyarcelly.213@gmail.com} (\href{https://orcid.org/0000-0003-3523-8687}{ORCID}) [contributor]
  \item Roy Yali \email{ryali93@gmail.com} (\href{https://orcid.org/0000-0003-4542-3755}{ORCID}) [contributor]
  \item Antony Barja \email{antony.barja8@gmail.com} (\href{https://orcid.org/0000-0001-5921-2858}{ORCID}) [contributor]
  \item Kevin Ushey \email{kevin@rstudio.com} [contributor]
  \item Jeroen Ooms \email{jeroen@berkeley.edu} (\href{https://orcid.org/0000-0002-4035-0289}{ORCID}) [contributor]
  \item Tim Appelhans \email{tim.appelhans@gmail.com} [contributor]
  \item JJ Allaire \email{jj@rstudio.com} [contributor]
  \item Yuan Tang \email{terrytangyuan@gmail.com} [contributor]
  \item Samapriya Roy \email{samapriya.roy@gmail.com} [contributor]
  \item MariaElena Adauto \email{2a.mariaelena@gmail.com} (\href{https://orcid.org/0000-0002-2154-2429}{ORCID}) [contributor]
  \item Gabriel Carrasco \email{gabriel.carrasco@upch.pe} (\href{https://orcid.org/0000-0002-6945-0419}{ORCID}) [contributor]
  \item Henrik Bengtsson \email{henrikb@braju.com} [contributor]
  \item Jeffrey Hollister \email{hollister.jeff@epa.gov} (Hollister reviewed the package for JOSS, see 
                        https://github.com/openjournals/joss-reviews/issues/2272/) [reviewer]
  \item Gennadii Donchyts (Gena reviewed the package for JOSS, see 
                        https://github.com/openjournals/joss-reviews/issues/2272/) [reviewer]
  \item Marius Appel \email{marius.appel@uni-muenster.de} (Appel reviewed the package for JOSS, see 
                        https://github.com/openjournals/joss-reviews/issues/2272/) [reviewer]
}

}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_download.R
\name{ee_image_to_gcs}
\alias{ee_image_to_gcs}
\title{Creates a task to export an EE Image to Google Cloud Storage.}
\usage{
ee_image_to_gcs(
  image,
  description = "myExportImageTask",
  bucket = NULL,
  fileNamePrefix = NULL,
  timePrefix = TRUE,
  dimensions = NULL,
  region = NULL,
  scale = NULL,
  crs = NULL,
  crsTransform = NULL,
  maxPixels = NULL,
  shardSize = NULL,
  fileDimensions = NULL,
  skipEmptyTiles = NULL,
  fileFormat = NULL,
  formatOptions = NULL
)
}
\arguments{
\item{image}{The image to be exported.}

\item{description}{Human-readable name of the task.}

\item{bucket}{The name of a Cloud Storage bucket for the export.}

\item{fileNamePrefix}{Cloud Storage object name prefix for the export.
Defaults to the name of the task.}

\item{timePrefix}{Add current date and time as a prefix to files to export.}

\item{dimensions}{The dimensions of the exported image. Takes either a
single positive integer as the maximum dimension or "WIDTHxHEIGHT"
where WIDTH and HEIGHT are each positive integers.}

\item{region}{The lon,lat coordinates for a LinearRing or Polygon
specifying the region to export. It can be specified as nested lists
of numbers or a serialized string. Defaults to the image's region.}

\item{scale}{The resolution in meters per pixel. Defaults to the native
resolution of the image assset unless a crsTransform is specified.}

\item{crs}{The coordinate reference system of the exported image's
projection. Defaults to the image's default projection.}

\item{crsTransform}{A comma-separated string of 6 numbers describing
the affine transform of the coordinate reference system of the exported
image's projection, in the order:
xScale, xShearing, xTranslation, yShearing, yScale, and yTranslation.
Defaults to the image's native CRS transform.}

\item{maxPixels}{The maximum allowed number of pixels in the
exported image. The task will fail if the exported region covers more
pixels in the specified projection. Defaults to 100,000,000.}

\item{shardSize}{Size in pixels of the shards in which this image
will be computed. Defaults to 256.}

\item{fileDimensions}{The dimensions in pixels of each image file, if
the image is too large to fit in a single file. May specify a single
number to indicate a square shape, or a list of two dimensions to
indicate (width, height). Note that the image will still be clipped to
the overall image dimensions. Must be a multiple of shardSize.}

\item{skipEmptyTiles}{If TRUE, skip writing empty (i.e., fully-masked)
image tiles. Defaults to FALSE.}

\item{fileFormat}{The string file format to which the image is exported.
Currently only 'GeoTIFF' and 'TFRecord' are supported, defaults
to 'GeoTIFF'.}

\item{formatOptions}{A dictionary of string keys to format-specific
options. **kwargs: Holds other keyword arguments that may have been
deprecated, such as 'crs_transform'.}
}
\value{
An unstarted Task that exports the image to Google Cloud Storage.
}
\description{
Creates a task to export an EE Image to Google Cloud Storage.
This function is a wrapper around
\code{ee$batch$Export$image$toCloudStorage(...)}.
}
\examples{
\dontrun{
library(rgee)
library(stars)
library(sf)

ee_users()
ee_Initialize(gcs = TRUE)

# Define study area (local -> earth engine)
# Communal Reserve Amarakaeri - Peru
rlist <- list(xmin = -71.13, xmax = -70.95,ymin = -12.89, ymax = -12.73)
ROI <- c(rlist$xmin, rlist$ymin,
         rlist$xmax, rlist$ymin,
         rlist$xmax, rlist$ymax,
         rlist$xmin, rlist$ymax,
         rlist$xmin, rlist$ymin)
ee_ROI <- matrix(ROI, ncol = 2, byrow = TRUE) \%>\%
  list() \%>\%
  st_polygon() \%>\%
  st_sfc() \%>\%
  st_set_crs(4326) \%>\%
  sf_as_ee()


# Get the mean annual NDVI for 2011
cloudMaskL457 <- function(image) {
  qa <- image$select("pixel_qa")
  cloud <- qa$bitwiseAnd(32L)$
    And(qa$bitwiseAnd(128L))$
    Or(qa$bitwiseAnd(8L))
  mask2 <- image$mask()$reduce(ee$Reducer$min())
  image <- image$updateMask(cloud$Not())$updateMask(mask2)
  image$normalizedDifference(list("B4", "B3"))
}

ic_l5 <- ee$ImageCollection("LANDSAT/LT05/C01/T1_SR")$
  filterBounds(ee$FeatureCollection(ee_ROI))$
  filterDate("2011-01-01", "2011-12-31")$
  map(cloudMaskL457)

# Create simple composite
mean_l5 <- ic_l5$mean()$rename("NDVI")
mean_l5 <- mean_l5$reproject(crs = "EPSG:4326", scale = 500)
mean_l5_Amarakaeri <- mean_l5$clip(ee_ROI)

# Move results from Earth Engine to GCS
task_img <- ee_image_to_gcs(
 image = mean_l5_Amarakaeri,
 bucket = "rgee_dev",
 fileFormat = "GEO_TIFF",
 region = ee_ROI,
 fileNamePrefix = "my_image_demo"
)

task_img$start()
ee_monitoring(task_img)

# Move results from GCS to local
ee_gcs_to_local(task = task_img)

}
}
\seealso{
Other image export task creator: 
\code{\link{ee_image_to_asset}()},
\code{\link{ee_image_to_drive}()}
}
\concept{image export task creator}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_get.R
\name{ee_get_earthengine_path}
\alias{ee_get_earthengine_path}
\title{Get the path where the credentials are stored}
\usage{
ee_get_earthengine_path()
}
\value{
A character that represents the path credential of a specific
user
}
\description{
Get the path where the credentials are stored
}
\seealso{
Other path utils: 
\code{\link{ee_get_assethome}()}
}
\concept{path utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_utils_cog_metadata}
\alias{ee_utils_cog_metadata}
\title{Return metadata of a COG tile server}
\usage{
ee_utils_cog_metadata(
  resource,
  visParams,
  titiler_server = "https://api.cogeo.xyz/"
)
}
\arguments{
\item{resource}{Character that represents a COG tile server file.}

\item{visParams}{Visualization parameters see "https://api.cogeo.xyz/docs".}

\item{titiler_server}{TiTiler endpoint. Defaults to "https://api.cogeo.xyz/".}
}
\value{
A metadata list for a COG file.
}
\description{
Return metadata of a COG tile server
}
\examples{
\dontrun{
 library(rgee)

server <- "https://s3-us-west-2.amazonaws.com/planet-disaster-data/hurricane-harvey/"
file <- "SkySat_Freeport_s03_20170831T162740Z3.tif"
resource <- paste0(server, file)
visParams <- list(nodata = 0, expression = "B3, B2, B1", rescale = "3000, 13500")
ee_utils_cog_metadata(resource, visParams)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_install.R
\name{ee_install_upgrade}
\alias{ee_install_upgrade}
\title{Upgrade the Earth Engine Python API}
\usage{
ee_install_upgrade(
  version = NULL,
  earthengine_env = Sys.getenv("EARTHENGINE_ENV")
)
}
\arguments{
\item{version}{Character. The Earth Engine Python API version to upgrade.
By default \code{rgee::ee_version()}.}

\item{earthengine_env}{Character. The name, or full path, of the
environment in which the earthengine-api packages are to be installed.}
}
\value{
no return value, called to upgrade the earthengine-api Python package
}
\description{
Upgrade the Earth Engine Python API
}
\seealso{
Other ee_install functions: 
\code{\link{ee_install_set_pyenv}()},
\code{\link{ee_install}()}
}
\concept{ee_install functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_image.R
\name{ee_as_raster}
\alias{ee_as_raster}
\title{Convert an Earth Engine (EE) image in a raster object}
\usage{
ee_as_raster(
  image,
  region = NULL,
  dsn = NULL,
  via = "drive",
  container = "rgee_backup",
  scale = NULL,
  maxPixels = 1e+09,
  lazy = FALSE,
  public = TRUE,
  add_metadata = TRUE,
  timePrefix = TRUE,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{image}{ee$Image to be converted into a raster object.}

\item{region}{EE Geometry (ee$Geometry$Polygon) which specifies the region
to export. CRS needs to be the same that the argument \code{image}.
Otherwise, it will be forced. If not specified, image bounds are taken.}

\item{dsn}{Character. Output filename. If missing, a temporary file is
created.}

\item{via}{Character. Method to export the image. Two methods are
implemented: "drive", "gcs". See details.}

\item{container}{Character. Name of the folder ('drive') or bucket ('gcs')
to be exported.}

\item{scale}{Numeric. The resolution in meters per pixel. Defaults
to the native resolution of the image.}

\item{maxPixels}{Numeric. The maximum allowed number of pixels in the
exported image. The task will fail if the exported region covers
more pixels in the specified projection. Defaults to 100,000,000.}

\item{lazy}{Logical. If TRUE, a \code{\link[future:sequential]{
future::sequential}} object is created to evaluate the task in the future.
See details.}

\item{public}{Logical. If TRUE, a public link to the image is created.}

\item{add_metadata}{Add metadata to the stars_proxy object. See details.}

\item{timePrefix}{Logical. Add current date and time (\code{Sys.time()}) as
a prefix to files to export. This parameter helps to avoid exported files
with the same name. By default TRUE.}

\item{quiet}{Logical. Suppress info message}

\item{...}{Extra exporting argument. See \link{ee_image_to_drive} and
\link{ee_image_to_gcs}.}
}
\value{
A RasterStack object
}
\description{
Convert an ee$Image in a raster object
}
\details{
\code{ee_as_raster} supports the download of \code{ee$Images}
by two different options: "drive"
(\href{https://CRAN.R-project.org/package=googledrive}{Google Drive}) and "gcs"
(\href{https://CRAN.R-project.org/package=googleCloudStorageR}{
Google Cloud Storage}). In both cases, \code{ee_as_stars} works as follow:
\itemize{
\item{1. }{A task is started (i.e., \code{ee$batch$Task$start()}) to
move the \code{ee$Image} from Earth Engine to the intermediate container
specified in the argument \code{via}.}
\item{2. }{If the argument \code{lazy} is TRUE, the task is not be
monitored. This is useful to lunch several tasks simultaneously and
calls them later using \code{\link{ee_utils_future_value}} or
\code{\link[future:value]{future::value}}. At the end of this step,
the \code{ee$Image} is stored on the path specified in the argument
\code{dsn}.}
\item{3. }{Finally, if the argument \code{add_metadata} is TRUE, a list
with the following elements are added to the stars-proxy object.
\itemize{
\item{\bold{if via is "drive":}}
\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{drive_name: }}{Name of the Image in Google Drive.}
\item{\bold{drive_id: }}{Id of the Image in Google Drive.}
\item{\bold{drive_download_link: }}{Download link to the image.}
}
}
\itemize{
\item{\bold{if via is "gcs":}}
\itemize{
\item{\bold{ee_id: }}{Name of the Earth Engine task.}
\item{\bold{gcs_name: }}{Name of the Image in Google Cloud Storage.}
\item{\bold{gcs_bucket: }}{Name of the bucket.}
\item{\bold{gcs_fileFormat: }}{Format of the image.}
\item{\bold{gcs_public_link: }}{Download link to the image.}
\item{\bold{gcs_URI: }}{gs:// link to the image.}
}
}
Run \code{raster@history@metadata} to get the list.
}
}

For getting more information about exporting data from Earth Engine, take
a look at the
\href{https://developers.google.com/earth-engine/guides/exporting}{Google
Earth Engine Guide - Export data}.
}
\examples{
\dontrun{
library(rgee)

ee_Initialize(drive = TRUE, gcs = TRUE)
ee_user_info()

# Define an image.
img <- ee$Image("LANDSAT/LC08/C01/T1_SR/LC08_038029_20180810")$
  select(c("B4", "B3", "B2"))$
  divide(10000)

# OPTIONAL display it using Map
Map$centerObject(eeObject = img)
Map$addLayer(eeObject = img, visParams = list(max = 0.4,gamma=0.1))

# Define an area of interest.
geometry <- ee$Geometry$Rectangle(
  coords = c(-110.8, 44.6, -110.6, 44.7),
  proj = "EPSG:4326",
  geodesic = FALSE
)

## drive - Method 01
# Simple
img_02 <- ee_as_raster(
  image = img,
  region = geometry,
  via = "drive"
)

# Lazy
img_02 <- ee_as_raster(
  image = img,
  region = geometry,
  via = "drive",
  lazy = TRUE
)

img_02_result <- img_02 \%>\% ee_utils_future_value()
img_02_result@history$metadata # metadata

## gcs - Method 02
# Simple
img_03 <- ee_as_raster(
 image = img,
 region = geometry,
 container = "rgee_dev",
 via = "gcs"
)

# Lazy
img_03 <- ee_as_raster(
 image = img,
 region = geometry,
 container = "rgee_dev",
 lazy = TRUE,
 via = "gcs"
)

img_03_result <- img_03 \%>\% ee_utils_future_value()
img_03_result@history$metadata # metadata

# OPTIONAL: clean containers
# ee_clean_container(name = "rgee_backup", type = "drive")
# ee_clean_container(name = "rgee_dev", type = "gcs")
}
}
\seealso{
Other image download functions: 
\code{\link{ee_as_stars}()},
\code{\link{ee_as_thumbnail}()},
\code{\link{ee_imagecollection_to_local}()}
}
\concept{image download functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raster_as_ee.R
\name{stars_as_ee}
\alias{stars_as_ee}
\title{Convert a stars or stars-proxy object into an EE Image object}
\usage{
stars_as_ee(
  x,
  assetId,
  bucket = NULL,
  predefinedAcl = "bucketLevel",
  command_line_tool_path = NULL,
  overwrite = FALSE,
  monitoring = TRUE,
  quiet = FALSE,
  ...
)
}
\arguments{
\item{x}{stars or stars-proxy object to be converted into an ee$Image.}

\item{assetId}{Character. Destination asset ID for the uploaded file.}

\item{bucket}{Character. Name of the GCS bucket.}

\item{predefinedAcl}{Specify user access to object. Passed to
\code{googleCloudStorageR::gcs_upload}.}

\item{command_line_tool_path}{Character. Path to the Earth Engine command line
tool (CLT). If NULL, rgee assumes that CLT is set in the system PATH.
(ignore if \code{via} is not defined as "gcs_to_asset").}

\item{overwrite}{Logical. If TRUE, the assetId will be overwritten.}

\item{monitoring}{Logical. If TRUE the exportation task will be monitored.}

\item{quiet}{Logical. Suppress info message.}

\item{...}{parameter(s) passed on to \code{\link{ee_utils_create_manifest_image}}}
}
\value{
An ee$Image object
}
\description{
Convert a stars or stars-proxy object into an EE Image object
}
\examples{
\dontrun{
library(rgee)
library(stars)
ee_Initialize(gcs = TRUE)

# Get the filename of a image
tif <- system.file("tif/L7_ETMs.tif", package = "stars")
x <- read_stars(tif)
assetId <- sprintf("\%s/\%s",ee_get_assethome(),'stars_l7')

# # Method 1
# 1. Move from local to gcs
gs_uri <- local_to_gcs(x = tif, bucket = 'rgee_dev')

# 2. Create a manifest
manifest <- ee_utils_create_manifest_image(gs_uri, assetId)

# 3. Pass from gcs to asset
gcs_to_ee_image(
  manifest = manifest,
  overwrite = TRUE
)

# OPTIONAL: Monitoring progress
ee_monitoring()

# OPTIONAL: Display results
ee_stars_01 <- ee$Image(assetId)
Map$centerObject(ee_stars_01)
Map$addLayer(ee_stars_01, list(min = 0, max = 255))

# Method 2
ee_stars_02 <- stars_as_ee(
 x = x,
 overwrite = TRUE,
 assetId = assetId,
 bucket = "rgee_dev"
)
Map$centerObject(ee_stars_02)
Map$addLayer(ee_stars_02, list(min = 0, max = 255))
}
}
\seealso{
Other image upload functions: 
\code{\link{raster_as_ee}()}
}
\concept{image upload functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_utils_dataset_display}
\alias{ee_utils_dataset_display}
\title{Search into the Earth Engine Data Catalog}
\usage{
ee_utils_dataset_display(ee_search_dataset)
}
\arguments{
\item{ee_search_dataset}{Character that represents the EE dataset ID.}
}
\value{
No return value, called for displaying the Earth Engine dataset in the browser.
}
\description{
Search into the Earth Engine Data Catalog
}
\examples{
\dontrun{
 library(rgee)

 ee_datasets <- c("WWF/HydroSHEDS/15DIR", "WWF/HydroSHEDS/03DIR")
 ee_utils_dataset_display(ee_datasets)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_get.R
\name{ee_get_date_img}
\alias{ee_get_date_img}
\title{Get the date of a EE Image}
\usage{
ee_get_date_img(x, time_end = FALSE)
}
\arguments{
\item{x}{ee$Image or ee$ImageCollection object}

\item{time_end}{Logical. If TRUE, the \code{system:time_end} property is
also returned. See details.}
}
\value{
An List object with the elements: id, time_start and time_end
(only if the \code{time_end} argument is TRUE).
}
\description{
Get the date of a EE Image
}
\details{
\code{system:time_start} sets the start period of data acquisition while
\code{system:time_end} does the same for the end period. See the
\href{https://developers.google.com/earth-engine/glossary/}{Earth Engine glossary}
for getting more information.
}
\examples{
\dontrun{
library(rgee)
ee_Initialize()

l8 <- ee$Image('LANDSAT/LC08/C01/T1_TOA/LC08_044034_20140318')
ee_get_date_img(l8)
srtm <- ee$Image('CGIAR/SRTM90_V4')
ee_get_date_img(srtm, time_end = TRUE)
}
}
\seealso{
Other date functions: 
\code{\link{ee_get_date_ic}()},
\code{\link{eedate_to_rdate}()},
\code{\link{rdate_to_eedate}()}
}
\concept{date functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_clean.R
\name{ee_clean_credentials}
\alias{ee_clean_credentials}
\title{Delete Credentials}
\usage{
ee_clean_credentials(user = "not_defined", quiet = FALSE)
}
\arguments{
\item{user}{Character. Earth Engine user (e.g. \code{data.colec.fbf}).}

\item{quiet}{Logical. Suppress info messages.}
}
\value{
No return value, called for cleaning Google Drive, Google Cloud Storage,
and/or Earth Engine credentials.
}
\description{
Delete all the credentials according to a specific user. The credentials
(Google Earth Engine, Google Drive and Google Cloud Storage) are created
after running \code{ee_Initialize(...)} successfully. They are saved in
the path \code{rgee::ee_get_earthengine_path()}.
}
\seealso{
Other ee_clean functions: 
\code{\link{ee_clean_container}()},
\code{\link{ee_clean_pyenv}()}
}
\concept{ee_clean functions}
