
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `prism`

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/prism)](https://cran.r-project.org/package=prism)
[![R build
status](https://github.com/ropensci/prism/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/prism/actions)
[![codecov.io](https://codecov.io/github/ropensci/prism/coverage.svg?branch=master)](https://codecov.io/github/ropensci/prism?branch=master)

This package allows users to access and visualize data from the [Oregon
State PRISM project](https://prism.nacse.org). Data are all in the form
of gridded rasters for the continental US at 4 different temporal
scales: daily, monthly, annual, and 30 year normals. Please see their
webpage for a full description of the data products, or [see their
overview](https://www.prism.oregonstate.edu/documents/PRISM_datasets_aug2013.pdf).

## Installation

prism is available on CRAN:

``` r
install.packages("prism")
```

Or the development version can be installed from GitHub with devtools:

``` r
# install.packages("devtools")
library(devtools)
install_github("ropensci/prism")
```

## Quickstart

The overall work flow in the prism package is (links go to details on
this page):

1.  [Set the download directory](#downloading-data), i.e., the folder on
    your computer that prism data will be saved to:
    `prism_set_dl_dir()`. This is now referred to as the “prism
    archive”.
2.  [Download prism data to the archive:](#download-30-year-normal-data)
    `get_prism_*()`. Each folder, or variable, timestep, day/month/year
    is stored in a single folder in the archive and referred to as prism
    data (`pd`).
3.  [Interact with the prism
    archive:](#interact-with-the-archive-and-prism-data)
    `prism_archive_*()`. Or interact with the prism data: `pd_*()`.

The remainder of this README provides examples following this work flow.

## prism data and parameters

Data are available in 4 different temporal scales as mentioned above. At
each temporal scale, there are 7 different parameters/variables
available. Keep in mind these are modeled parameters, not measured.
Please see the [full
description](https://www.prism.oregonstate.edu/documents/Daly2008_PhysiographicMapping_IntJnlClim.pdf)
for how they are calculated.

| Parameter name | Description                          |
| :------------- | :----------------------------------- |
| *tmean*        | Mean temperature                     |
| *tmax*         | Maximum temperature                  |
| *tmin*         | Minimum temperature                  |
| *tdmean*       | Mean dew point temperature           |
| *ppt*          | Total precipitation (rain and snow)  |
| *vpdmin*       | Daily minimum vapor pressure deficit |
| *vpdmax*       | Daily maximum vapor pressure deficit |

## Downloading data

Before downloading any data, set the directory that the prism data will
be saved to:

``` r
library(prism)
#> Be sure to set the download folder using `prism_set_dl_dir()`.
prism_set_dl_dir("~/prismtmp")
```

This is now referred to as the “prism archive”. The `prism_archive_*()`
functions allow the user to search through the archive. The prism
archive contains “prism data”. The prism data are referred to by their
folder names, even though the “real” data are the .bil, .txt, and other
files that exist in the folder. The prism data (`pd`) can be accessed
using the `pd_*()` functions.

### Download 30-year normal data

Normals are based on the latest 30-year period; currently 1981 - 2010.
Normals can be downloaded in two resolutions, 4km and 800m, and a
resolution must be specified. They can be downloaded for a given month,
vector of months, or annual averages for all 30 years.

``` r
# Download the January - June 30-year averages at 4km resolution
get_prism_normals(type="tmean", resolution = "4km", mon = 1:6, keepZip = FALSE)

# Download the 30-year annual average precip and annual average temperature
get_prism_normals("ppt", "4km", annual = TRUE, keepZip = FALSE)
get_prism_normals("tmean", "4km", annual = TRUE, keepZip = FALSE)
```

If the archive has not already been set, calling any of the
`get_prism_*()` functions will prompt the user to specify the directory.
prism data are downloaded as zip files and then unzipped. If the
`keepZip` argument is `TRUE` the zip file will remain on your machine,
otherwise it will be automatically deleted.

### Download daily, monthly, and annual data

Let us download daily average temperatures from June 1 to June 14, 2013.
We can also download January average temperature data from 1982 to 2014.
Finally, we will download annual average precipitation for 2000 to 2015.

``` r
get_prism_dailys(
  type = "tmean", 
  minDate = "2013-06-01", 
  maxDate = "2013-06-14", 
  keepZip = FALSE
)
get_prism_monthlys(type = "tmean", year = 1982:2014, mon = 1, keepZip = FALSE)
get_prism_annual("ppt", years = 2000:2015, keepZip = FALSE)
```

Note that for daily data you need to give a well formed date string in
the form of “YYYY-MM-DD”.

## Interact with the archive and prism data

You can view all the prism data you have downloaded with a simple
command: `prism_archive_ls()`. This function gives a list of folder
names, i.e., prism data (`pd`). All the functions in the prism package
work off of one or more of these folder names (`pd`).

``` r
## Truncated to keep file list short
prism_archive_ls()
#>  [1] "PRISM_ppt_30yr_normal_4kmM2_annual_bil"  
#>  [2] "PRISM_ppt_30yr_normal_800mM2_02_bil"     
#>  [3] "PRISM_ppt_stable_4kmD2_19810101_bil"     
#>  [4] "PRISM_ppt_stable_4kmD2_19820101_bil"     
#>  [5] "PRISM_ppt_stable_4kmD2_19830101_bil"     
#>  [6] "PRISM_ppt_stable_4kmD2_20120101_bil"     
#>  [7] "PRISM_ppt_stable_4kmM3_2000_bil"         
#>  [8] "PRISM_ppt_stable_4kmM3_2001_bil"         
#>  [9] "PRISM_ppt_stable_4kmM3_2002_bil"         
#> [10] "PRISM_ppt_stable_4kmM3_2003_bil"         
....
```

While prism functions use this folder format, other files may need an
absolute path (e.g. the `raster` package). The `pd_to_file()` function
conveniently returns the absolute path. Alternatively, you may want to
see what the normal name for the product is (not the file name), and we
can get that with the `pd_get_name()` function.

``` r
## Truncated to keep file list short
pd_to_file(prism_archive_ls())
#>  [1] "C:\\Users\\RAButler\\Documents\\prismtmp\\PRISM_ppt_30yr_normal_4kmM2_annual_bil\\PRISM_ppt_30yr_normal_4kmM2_annual_bil.bil"    
#>  [2] "C:\\Users\\RAButler\\Documents\\prismtmp\\PRISM_ppt_30yr_normal_800mM2_02_bil\\PRISM_ppt_30yr_normal_800mM2_02_bil.bil"          
#>  [3] "C:\\Users\\RAButler\\Documents\\prismtmp\\PRISM_ppt_stable_4kmD2_19810101_bil\\PRISM_ppt_stable_4kmD2_19810101_bil.bil"          
#>  [4] "C:\\Users\\RAButler\\Documents\\prismtmp\\PRISM_ppt_stable_4kmD2_19820101_bil\\PRISM_ppt_stable_4kmD2_19820101_bil.bil"          
#>  [5] "C:\\Users\\RAButler\\Documents\\prismtmp\\PRISM_ppt_stable_4kmD2_19830101_bil\\PRISM_ppt_stable_4kmD2_19830101_bil.bil"          
....

pd_get_name(prism_archive_ls())
#>  [1] "Annual 30-year normals - 4km resolution - Precipitation"      
#>  [2] "Feb 30-year normals - 800m resolution - Precipitation"        
#>  [3] "Jan 01 1981 - 4km resolution - Precipitation"                 
#>  [4] "Jan 01 1982 - 4km resolution - Precipitation"                 
#>  [5] "Jan 01 1983 - 4km resolution - Precipitation"                 
....
```

Finally, `prism_archive_subset()` is a convenient way to search for
specific parameters, time steps, days, months, years, or ranges of days,
months, years.

``` r
# we know we have downloaded June 2013 daily data, so lets search for those 
prism_archive_subset("tmean", "daily", mon = 6)
#>  [1] "PRISM_tmean_stable_4kmD2_20130601_bil"
#>  [2] "PRISM_tmean_stable_4kmD2_20130602_bil"
#>  [3] "PRISM_tmean_stable_4kmD2_20130603_bil"
#>  [4] "PRISM_tmean_stable_4kmD2_20130604_bil"
#>  [5] "PRISM_tmean_stable_4kmD2_20130605_bil"
#>  [6] "PRISM_tmean_stable_4kmD2_20130606_bil"
#>  [7] "PRISM_tmean_stable_4kmD2_20130607_bil"
#>  [8] "PRISM_tmean_stable_4kmD2_20130608_bil"
#>  [9] "PRISM_tmean_stable_4kmD2_20130609_bil"
#> [10] "PRISM_tmean_stable_4kmD2_20130610_bil"
#> [11] "PRISM_tmean_stable_4kmD2_20130611_bil"
#> [12] "PRISM_tmean_stable_4kmD2_20130612_bil"
#> [13] "PRISM_tmean_stable_4kmD2_20130613_bil"
#> [14] "PRISM_tmean_stable_4kmD2_20130614_bil"

# or we can look for days between June 7 and June 10
prism_archive_subset(
  "tmean", "daily", minDate = "2013-06-07", maxDate = "2013-06-10"
)
#> [1] "PRISM_tmean_stable_4kmD2_20130607_bil"
#> [2] "PRISM_tmean_stable_4kmD2_20130608_bil"
#> [3] "PRISM_tmean_stable_4kmD2_20130609_bil"
#> [4] "PRISM_tmean_stable_4kmD2_20130610_bil"
```

### Raster plots

You can easily make a quick plot of your data using the output of
`prism_archive_ls()` or `prism_archive_subset()` with `pd_image()`.

``` r
# Plot the January 30-year average temperatures
jmean <- prism_archive_subset(
  "tmean", "monthly normals", mon = 1, resolution = "4km"
)
pd_image(jmean)
```

![](man/figures/README-quick_plot-1.png)<!-- -->

It is easy to load the prism data with the raster package. This time we
will look at January temperature anomalies. To do this we will examine
the difference between January 2013 and the January 30-year normals.
Conveniently, we already downloaded these data. We just need to grab
them out of our archive.

``` r
library(raster)
#> Loading required package: sp
# knowing the name of the files you are after allows you to find them in the 
# list of all files that exist
# jnorm_name <- "PRISM_tmean_30yr_normal_4kmM2_01_bil"
# j2013_name <- "PRISM_tmean_stable_4kmM3_201301_bil"
# but we will use prism_archive_subset() to find the files we need

jnorm <- prism_archive_subset(
  "tmean", "monthly normals", mon = 1, resolution = "4km"
)
j2013 <- prism_archive_subset("tmean", "monthly", years = 2013, mon = 1)

# raster needs a full path, not the "short" prism data name
jnorm <- pd_to_file(jnorm)
j2013 <- pd_to_file(j2013)

## Now we'll load the rasters.
jnorm_rast <- raster(jnorm)
j2013_rast <- raster(j2013)

# Now we can do simple subtraction to get the anomaly by subtracting 2014 
# from the 30 year normal map
anomCalc <- function(x, y) {
  return(x - y)
}

anom_rast <- raster::overlay(j2013_rast,jnorm_rast,fun = anomCalc)

plot(anom_rast)
```

![](man/figures/README-raster_math-1.png)<!-- -->

The plot shows that January 2013 was warmer than the average over the
last 30 years. It also shows how easy it is to use the raster library to
work with prism data. The package provides a simple framework to work
with a large number of rasters that you can easily download and
visualize or use with other data sets.

### Single grid cell plot

You can also visualize a single point across multiple prism data files
(slice) using `pd_plot_slice()`. This procedure will take a set of
rasters, create a “`raster::stack`”, extract data at a point, and then
create a ggplot2 object.

Let’s now make a plot of January temperatures in Boulder between 1982
and 2014. First we’ll grab all the data from the US (downloaded in the
previous step), and then give our function a point to get data from. The
point must be a vector in the form of longitude, latitude. Because
`pd_plot_slice()` returns a gg object, it can be combined with other
ggplot functions.

``` r
library(ggplot2)
# data already exist in the prism dl dir
boulder <- c(-105.2797, 40.0176)

# prism_archive_subset() will return prism data that matches the specified 
# variable, time step, years, months, days, etc.
to_slice <- prism_archive_subset("tmean", "monthly", mon = 1)
p <- pd_plot_slice(to_slice, boulder)

# add a linear average and title
p + 
  stat_smooth(method="lm", se = FALSE) + 
  theme_bw() + 
  ggtitle("Average January temperature in Boulder, CO 1982-2014")
#> `geom_smooth()` using formula 'y ~ x'
```

![](man/figures/README-plot_Boulder-1.png)<!-- -->

### leaflet map

Finally, the prism data are in a form that can be used with leaflet maps
(with the help of the raster package). The [leaflet
package](https://CRAN.R-project.org/package=leaflet) allows you to
easily make JavaScript maps using the [leaflet](https://leafletjs.com/)
mapping framework using prism data. These can easily be hosted on
websites like [Rpubs](https://rpubs.com/) or your own site. Here is a
simple example of plotting the [30-year normal for annual
temperature](https://rpubs.com/DistribEcology/122453). If you run this
code you will have an interactive map, instead of just the screen shot
shown here.

``` r
library(leaflet)
library(raster)
library(prism)

# 30-year normal average temperature have already been downloaded for 
norm <- prism_archive_subset(
  "tmean", "annual normals", resolution = "4km"
)
rast <- raster(pd_to_file(norm))

# Create color palette and plot
pal <- colorNumeric(
  c("#0000FF", "#FFFF00", "#FF0000"), 
  values(rast),
  na.color = "transparent"
)

leaflet() %>% 
  addTiles(
    urlTemplate = 'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'
  ) %>% 
  addRasterImage(rast, colors = pal, opacity=.65) %>% 
  addLegend(pal = pal, values = values(rast), title = "Deg C")
```

\[![leaflet example figure](vignettes/leaflet_example.png)
# prism 0.2.0

**Published December 5, 2020**

Thanks to @jsta for updating the README and helping with several other under the hood fixes. 

## Breaking Changes

* `prism_webservice()` is no longer exported as it is wrapped by `get_prism_*()` functions, and requires a correctly specified url. It can still be called with `prism:::prism_webservice()` if users really need it. (#83)
* `pr_parse()` is no longer exported. Use `pd_get_name()` or `pd_get_date()` instead. 


## Major Updates

There are two overall major updates with this release. (1) All functions should work with all temporal periods and all variables; previously some functions only worked with daily data and not all variables were able to be downloaded from the prism website. (2) A new API was implemented that results in many functions being deprecated in favor of the updated naming convention. This change was intended to provide consistent names for functions that apply to different steps in the work flow implemented in this package. The details of these changes are:

* Users can now download vpdmin, vpdmax, and tdmean variables for 30-year normals, daily, monthly, and annual data. (#68)
* There are now several functions (`prism_*_dl_dir()`) that set and check the prism download directory. These are hopefully easier to remember than using the base R `options()` and `getOption()` functions and the prism option variable name "prism.path". 
  * `prism_check_dl_dir()` replaces `check_path()`, which is deprecated and will be removed in the next release.
* The prism data are downloaded to this directory and then referred to as the "prism archive". These are the `prism_archive_*()` functions.
  * New function: `prism_archive_subset()`. This makes it much easier to get the data for a specific type/temporal period from the prism archive. (#69)
  * `prism_archive_ls()` replaces `ls_prism_data()`, which will be removed in a future release. 
    * The return type and options changed from `ls_prism_data()` to `prism_archive_ls()`, which now always returns only folder names as a vector, instead of a data.frame that could have between 1 and 3 columns.
    * The previous behavior can be achieved by applying new functions (`pd_get_name()` and `pd_to_file()`) to the vector returned by `prism_archive_ls()`. 
  * `prism_archive_clean()` replaces `del_early_prov()`, which will be removed in a future release. It also now works with all time steps and prompts user to select which folders will be removed before removing them (when R is in interactive mode). (#89)
  * `prism_archive_verify()` replaces `check_corrupt()`, which will be removed in a future release. It also now works with time steps other than daily and it gains a `download_corrupt` argument that controls whether corrupt files are automatically re-downloaded.
* `prism_archive_ls()` and `prism_archive_pd()` both return vectors of prism data folder names, i.e., prism data, i.e., `pd`. There are a number of functions that act on the prism data. These are the `pd_*()` functions. 
  * `pd_image()` replaces `prism_image()`.
  * `pd_plot_slice()` replaces `prism_slice()`.
  * `pd_stack()` replaces `prism_stack()`.
  * `pd_get_station_md()` replaces `get_prism_station_md()`. (#88)
  * `prism_md()` will be removed in a future release. It is replaced by:
    * `pd_get_name()`, which is equivalent to `prism_md(f, FALSE)`.
    * `pd_get_date()`, which is equivalent to `prism_md(f, TRUE)`.
  * Two new functions were added that convert the prism data to a full absolute path (`pd_to_file()`) and get the type (parameter) of the prism data (`pd_get_type()`).
  

## Other Changes

* The way pre 1981 data is handled has been updated. 
  * `get_prism_annual()` and `get_prism_monthlys()` gain a `keep_pre81_months` parameter. This lets the user determine if all of the monthly and annual data are kept, since the download includes all 12 months + the annual data for years before 1981. If this is `TRUE` then all monthly data are kept, instead of only those that were specified in the current call to `get_prism_*()`. (#82)
  * Because pre-1981 data might already have been downloaded based on `keep_pre81_months` parameter in previous downloads, the download functions now check that pre-1981 data does not exist before downloading it. To do this, `prism_webservice()` and `prism_check()` gain a `pre81_months` parameter, which allows the functions to know which months were requested for downloading. (#81)
* The prism web service only allows a user to download the same data twice in a 24-hour period. The download functions now report when the user has exceeded the allowable number of attempts to download the same file (in one day). If a user tries to download the same file more than two times in one day, the message is posted as a warning and the returned text file is saved in the prism archive. This is not posted as an error so that if a query of multiple files runs into this issue, it does not abort the full query. Another warning posts if the unzipped folder is empty. (#80)
* `prism_check()` is deprecated and will be no longer be exported in the next release.
* `get_prism_normals()` will now error if neither monthly nor annual data are specified to be downloaded and will download monthly and annual data simultaneously if asked to do so. (#77)
* all `get_prism_*()` functions are documented in same help page. (#79)
* Help pages for non-exported functions have been removed.
* More tests were added. Up to 50% coverage now. Some tests are only run locally to ensure prism download limits are not exceeded. 
* `pd_get_station_md()` (formerly `get_prism_station_md()`) now reports a warning if not all requested dates exist in the metadata data frame. (#87 related). It also now works for monthly and normals; not solely daily prism data. 
* `pd_get_md()` was added to parse .info.txt metadata, by converting an existing internal function. (#88)
* `prism_archive_clean()` (formerly `del_early_prov()`) now invisibly returns the folders that it removes.
* `pd_image()` (formerly `prism_image()`) invisibly returns the `gg` object it creates. It also shows the units for the prism variable in the fill legend. (#99)


# prism 0.1.0

## Minor changes

* New functions
    - `del_early_prov()` searches the download folder for duplicated PRISM data and keeps only the newest version.
    - `get_prism_station_md()` extracts metadata from daily PRISM data.
* `get_prism_dailys()` gains a `check` parameter that allows the user to specify how prism files are checked.


## Bug fixes

* `get_prism_monthlys()` can now download 1981 data. (@sdtaylor #59, #63)
* `get_prism_annual()` can now download pre 1981 data by itself. (@rabutler #64)
* `get_prism_dailys()` now correctly sets the progress bar.
* fixed bug in `gen_dates()` so that `get_prism_dailys()` works with only the `dates` parameter specified. (@rabutler #66)

## Under the hood

* added internal `gen_dates()` function for determining the specified dates (either from `minDate` and `maxDate`. or `dates`) used by the `get_prism_*()` functions.
* added tests for `gen_dates()` .

# prism 0.0.7

### Changes

* Changed acquisition method to use prism webservice.  Prior to this change the FTP acquisition method often caused the server to time out when download request volume became too high.  

* Changed method of metadata extraction.  Originally we parsed XML metadata to get information about raster files.  However the XML for historical data is particularly sparse. The new method relies on parsing file names instead. While more universal, it may be less stable

### Bug fixes

* Fixed FTP time out error by switching to webservice

* Fixed `prism_stack` by adjusting to new metadata extraction methodRe-submission of submission from 2020-11-10. The tarball is now < 5 MB.

## Test environment

* Local: Windows 10 Enterprise: R 3.6.3
* Ubuntu 20.04 (using GitHub Actions): R 3.6.3, 4.0.3, and R-devel
* MacOS-latest (using GitHub Actions): R 3.6.3, 4.0.3, and R-devel
* Windows (using GitHub Actions): R 3.6.3, R 4.0.3, and R-devel; Win-builder R 4.0.3 and R-devel
* R-Hub: Windows, R-devel; Ubuntu Linux 16.04, R-release; Fedora Linux, R-devel

## R CMD check results

There were no ERRORs or WARNINGs or NOTEs.

## Downstream dependencies

There are no downstream dependencies.
ppt - 2012-01-01 is corrupt/bad data in the bil file. (It was intentionally corrupted to test one of the functions.)prism/

- files in this folder are downloaded from the prism archive

small_prism/

- files in this folder are created using create_test_data.R by cropping the original files to only include California based lat/long extents. These files are then used in the tests. To get around the RCMD checks on vignette which would fail b/c the PRISM data
only exists locally, I added a local PRISM_AUTHOR environment variable. 

This environment variable must be set to "true" locally, if you wish to build 
the vignette. 

Otherwise, the figures and output will not be correctly rendered. ---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-"
)

# knitr hook to truncate long output. See:
# https://stackoverflow.com/q/23114654/3362993
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- unlist(strsplit(x, "\n"))
    if (length(x) > n) {
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})
```

# `prism`

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/prism)](https://cran.r-project.org/package=prism)
[![R build status](https://github.com/ropensci/prism/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/prism/actions)
[![codecov.io](https://codecov.io/github/ropensci/prism/coverage.svg?branch=master)](https://codecov.io/github/ropensci/prism?branch=master)


This package allows users to access and visualize data from the [Oregon State PRISM project](https://prism.nacse.org).  Data are all in the form of gridded rasters for the continental US at 4 different temporal scales: daily, monthly, annual, and 30 year normals.  Please see their webpage for a full description of the data products, or [see their overview](https://www.prism.oregonstate.edu/documents/PRISM_datasets_aug2013.pdf).

## Installation 

prism is available on CRAN:

```{r eval=FALSE}
install.packages("prism")
```

Or the development version can be installed from GitHub with devtools:
```{r start, eval=FALSE}
# install.packages("devtools")
library(devtools)
install_github("ropensci/prism")
```

## Quickstart

The overall work flow in the prism package is (links go to details on this page):

1. [Set the download directory](#downloading-data), i.e., the folder on your computer that prism data will be saved to: `prism_set_dl_dir()`. This is now referred to as the "prism archive". 
2. [Download prism data to the archive:](#download-30-year-normal-data) `get_prism_*()`. Each folder, or variable, timestep, day/month/year is stored in a single folder in the archive and referred to as prism data (`pd`). 
3. [Interact with the prism archive:](#interact-with-the-archive-and-prism-data) `prism_archive_*()`. Or interact with the prism data: `pd_*()`. 

The remainder of this README provides examples following this work flow.

## prism data and parameters

Data are available in 4 different temporal scales as mentioned above. At each temporal scale, there are 7 different parameters/variables available. Keep in mind these are modeled parameters, not measured.  Please see the [full description](https://www.prism.oregonstate.edu/documents/Daly2008_PhysiographicMapping_IntJnlClim.pdf) for how they are calculated.

| Parameter name| Description           |
|:---------------|:-------------|
| *tmean*      | Mean temperature |
| *tmax*      | Maximum temperature      |
| *tmin* | Minimum temperature      |
| *tdmean* | Mean dew point temperature |
| *ppt*  | Total precipitation (rain and snow)|
| *vpdmin* | Daily minimum vapor pressure deficit |
| *vpdmax* |Daily maximum vapor pressure deficit |

## Downloading data

Before downloading any data, set the directory that the prism data will be saved to:

```{r prism setup}
library(prism)
prism_set_dl_dir("~/prismtmp")
```

This is now referred to as the "prism archive". The `prism_archive_*()` functions allow the user to search through the archive. The prism archive contains "prism data". The prism data are referred to by their folder names, even though the "real" data are the .bil, .txt, and other files that exist in the folder. The prism data (`pd`) can be accessed using the `pd_*()` functions. 

### Download 30-year normal data

Normals are based on the latest 30-year period; currently 1981 - 2010. Normals can be downloaded in two resolutions, 4km and 800m, and a resolution must be specified.  They can be downloaded for a given month, vector of months, or annual averages for all 30 years.

```{r get normals,results=FALSE, eval=FALSE}
# Download the January - June 30-year averages at 4km resolution
get_prism_normals(type="tmean", resolution = "4km", mon = 1:6, keepZip = FALSE)

# Download the 30-year annual average precip and annual average temperature
get_prism_normals("ppt", "4km", annual = TRUE, keepZip = FALSE)
get_prism_normals("tmean", "4km", annual = TRUE, keepZip = FALSE)
```

If the archive has not already been set, calling any of the `get_prism_*()` functions will prompt the user to specify the directory. prism data are downloaded as zip files and then unzipped. If the `keepZip` argument is `TRUE` the zip file will remain on your machine, otherwise it will be automatically deleted.

### Download daily, monthly, and annual data

Let us download daily average temperatures from June 1 to June 14, 2013. We can also download January average temperature data from 1982 to 2014. Finally, we will download annual average precipitation for 2000 to 2015. 

```{r get daily monthly, message=FALSE, results=FALSE, eval=FALSE}
get_prism_dailys(
  type = "tmean", 
  minDate = "2013-06-01", 
  maxDate = "2013-06-14", 
  keepZip = FALSE
)
get_prism_monthlys(type = "tmean", year = 1982:2014, mon = 1, keepZip = FALSE)
get_prism_annual("ppt", years = 2000:2015, keepZip = FALSE)
```

Note that for daily data you need to give a well formed date string in the form of "YYYY-MM-DD".

## Interact with the archive and prism data

You can view all the prism data you have downloaded with a simple command: `prism_archive_ls()`.  This function gives a list of folder names, i.e., prism data (`pd`).  All the functions in the prism package work off of one or more of these folder names (`pd`).

```{r listingFiles, out.lines=10}
## Truncated to keep file list short
prism_archive_ls()
```

While prism functions use this folder format, other files may need an absolute path (e.g. the `raster` package). The `pd_to_file()` function conveniently returns the absolute path.  Alternatively, you may want to see what the normal name for the product is (not the file name), and we can get that with the `pd_get_name()` function.

```{r moreListing, out.lines=5}
## Truncated to keep file list short
pd_to_file(prism_archive_ls())

pd_get_name(prism_archive_ls())
```

Finally, `prism_archive_subset()` is a convenient way to search for specific parameters, time steps, days, months, years, or ranges of days, months, years. 

```{r}
# we know we have downloaded June 2013 daily data, so lets search for those 
prism_archive_subset("tmean", "daily", mon = 6)

# or we can look for days between June 7 and June 10
prism_archive_subset(
  "tmean", "daily", minDate = "2013-06-07", maxDate = "2013-06-10"
)
```

### Raster plots

You can easily make a quick plot of your data using the output of `prism_archive_ls()` or `prism_archive_subset()` with `pd_image()`. 

```{r quick_plot,fig.height=5,fig.width=7}
# Plot the January 30-year average temperatures
jmean <- prism_archive_subset(
  "tmean", "monthly normals", mon = 1, resolution = "4km"
)
pd_image(jmean)
```

It is easy to load the prism data with the raster package. This time we will look at January temperature anomalies.  To do this we will examine the difference between January 2013 and the January 30-year normals. Conveniently, we already downloaded these data. We just need to grab them out of our archive.

```{r raster_math,fig.height=5,fig.width=7}
library(raster)
# knowing the name of the files you are after allows you to find them in the 
# list of all files that exist
# jnorm_name <- "PRISM_tmean_30yr_normal_4kmM2_01_bil"
# j2013_name <- "PRISM_tmean_stable_4kmM3_201301_bil"
# but we will use prism_archive_subset() to find the files we need

jnorm <- prism_archive_subset(
  "tmean", "monthly normals", mon = 1, resolution = "4km"
)
j2013 <- prism_archive_subset("tmean", "monthly", years = 2013, mon = 1)

# raster needs a full path, not the "short" prism data name
jnorm <- pd_to_file(jnorm)
j2013 <- pd_to_file(j2013)

## Now we'll load the rasters.
jnorm_rast <- raster(jnorm)
j2013_rast <- raster(j2013)

# Now we can do simple subtraction to get the anomaly by subtracting 2014 
# from the 30 year normal map
anomCalc <- function(x, y) {
  return(x - y)
}

anom_rast <- raster::overlay(j2013_rast,jnorm_rast,fun = anomCalc)

plot(anom_rast)
```

The plot shows that January 2013 was warmer than the average over the last 30 years. It also shows how easy it is to use the raster library to work with prism data. The package provides a simple framework to work with a large number of rasters that you can easily download and visualize or use with other data sets.

### Single grid cell plot

You can also visualize a single point across multiple prism data files (slice) using `pd_plot_slice()`. This procedure will take a set of rasters, create a "`raster::stack`", extract data at a point, and then create a ggplot2 object.

Let's now make a plot of January temperatures in Boulder between 1982 and 2014. First we'll grab all the data from the US (downloaded in the previous step), and then give our function a point to get data from. The point must be a vector in the form of longitude, latitude. Because `pd_plot_slice()` returns a gg object, it can be combined with other ggplot functions.

```{r plot_Boulder,fig.height=5,fig.width=7, results=FALSE}
library(ggplot2)
# data already exist in the prism dl dir
boulder <- c(-105.2797, 40.0176)

# prism_archive_subset() will return prism data that matches the specified 
# variable, time step, years, months, days, etc.
to_slice <- prism_archive_subset("tmean", "monthly", mon = 1)
p <- pd_plot_slice(to_slice, boulder)

# add a linear average and title
p + 
  stat_smooth(method="lm", se = FALSE) + 
  theme_bw() + 
  ggtitle("Average January temperature in Boulder, CO 1982-2014")
```

### leaflet map

Finally, the prism data are in a form that can be used with leaflet maps (with the help of the raster package). The [leaflet package](https://CRAN.R-project.org/package=leaflet) allows you to easily make JavaScript maps using the [leaflet](https://leafletjs.com/) mapping framework using prism data.  These can easily be hosted on websites like [Rpubs](https://rpubs.com/) or your own site.  Here is a simple example of plotting the [30-year normal for annual temperature](https://rpubs.com/DistribEcology/122453). If you run this code you will have an interactive map, instead of just the screen shot shown here. 

```{r leaflet,eval=F}
library(leaflet)
library(raster)
library(prism)

# 30-year normal average temperature have already been downloaded for 
norm <- prism_archive_subset(
  "tmean", "annual normals", resolution = "4km"
)
rast <- raster(pd_to_file(norm))

# Create color palette and plot
pal <- colorNumeric(
  c("#0000FF", "#FFFF00", "#FF0000"), 
  values(rast),
  na.color = "transparent"
)

leaflet() %>% 
  addTiles(
    urlTemplate = 'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'
  ) %>% 
  addRasterImage(rast, colors = pal, opacity=.65) %>% 
  addLegend(pal = pal, values = values(rast), title = "Deg C")
```

[![leaflet example figure](vignettes/leaflet_example.png)
---
title: "Download and plot PRISM data"
author: "Alan Butler"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Download and plot PRISM data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = Sys.getenv("PRISM_AUTHOR") == "true"
)

# knitr hook to truncate long output. See:
# https://stackoverflow.com/q/23114654/3362993
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- unlist(strsplit(x, "\n"))
    if (length(x) > n) {
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})
```

**Previous Versions:**

- Edmund Hart - November 13, 2015


The prism package allows users to access and visualize data from the [Oregon State PRISM project](https://prism.nacse.org).  Data are all in the form of gridded rasters for the continental US at 4 different temporal scales: daily, monthly, annual, and 30 year normals.  Please see their webpage for a full description of the data products, or [see their overview](https://www.prism.oregonstate.edu/documents/PRISM_datasets_aug2013.pdf).

## Installation 

prism is available on CRAN:

```{r eval=FALSE}
install.packages("prism")
```

Or the development version can be installed from GitHub with devtools:
```{r start, eval=FALSE}
# install.packages("devtools")
library(devtools)
install_github("ropensci/prism")
```

## Quickstart

The overall work flow in the prism package is (links go to details on this page):

1. [Set the download directory](#downloading-data), i.e., the folder on your computer that prism data will be saved to: `prism_set_dl_dir()`. This is now referred to as the "prism archive". 
2. [Download prism data to the archive:](#download-30-year-normal-data) `get_prism_*()`. Each folder, or variable, timestep, day/month/year is stored in a single folder in the archive and referred to as prism data (`pd`). 
3. [Interact with the prism archive:](#interact-with-the-archive-and-prism-data) `prism_archive_*()`. Or interact with the prism data: `pd_*()`. 

The remainder of this README provides examples following this work flow.

## prism data and parameters

Data are available in 4 different temporal scales as mentioned above. At each temporal scale, there are 7 different parameters/variables available. Keep in mind these are modeled parameters, not measured.  Please see the [full description](https://www.prism.oregonstate.edu/documents/Daly2008_PhysiographicMapping_IntJnlClim.pdf) for how they are calculated.

| Parameter name| Description           |
|:---------------|:-------------|
| *tmean*      | Mean temperature |
| *tmax*      | Maximum temperature      |
| *tmin* | Minimum temperature      |
| *tdmean* | Mean dew point temperature |
| *ppt*  | Total precipitation (rain and snow)|
| *vpdmin* | Daily minimum vapor pressure deficit |
| *vpdmax* |Daily maximum vapor pressure deficit |

## Downloading data

Before downloading any data, set the directory that the prism data will be saved to:

```{r prism setup}
library(prism)
prism_set_dl_dir("~/prismtmp")
```

This is now referred to as the "prism archive". The `prism_archive_*()` functions allow the user to search through the archive. The prism archive contains "prism data". The prism data are referred to by their folder names, even though the "real" data are the .bil, .txt, and other files that exist in the folder. The prism data (`pd`) can be accessed using the `pd_*()` functions. 

### Download 30-year normal data

Normals are based on the latest 30-year period; currently 1981 - 2010. Normals can be downloaded in two resolutions, 4km and 800m, and a resolution must be specified.  They can be downloaded for a given month, vector of months, or annual averages for all 30 years.

```{r get normals,results=FALSE, eval=FALSE}
# Download the January - June 30-year averages at 4km resolution
get_prism_normals(type="tmean", resolution = "4km", mon = 1:6, keepZip = FALSE)

# Download the 30-year annual average precip and annual average temperature
get_prism_normals("ppt", "4km", annual = TRUE, keepZip = FALSE)
get_prism_normals("tmean", "4km", annual = TRUE, keepZip = FALSE)
```

If the archive has not already been set, calling any of the `get_prism_*()` functions will prompt the user to specify the directory. prism data are downloaded as zip files and then unzipped. If the `keepZip` argument is `TRUE` the zip file will remain on your machine, otherwise it will be automatically deleted.

### Download daily, monthly, and annual data

Let us download daily average temperatures from June 1 to June 14, 2013. We can also download January average temperature data from 1982 to 2014. Finally, we will download annual average precipitation for 2000 to 2015. 

```{r get daily monthly, message=FALSE, results=FALSE, eval=FALSE}
get_prism_dailys(
  type = "tmean", 
  minDate = "2013-06-01", 
  maxDate = "2013-06-14", 
  keepZip = FALSE
)
get_prism_monthlys(type = "tmean", year = 1982:2014, mon = 1, keepZip = FALSE)
get_prism_annual("ppt", years = 2000:2015, keepZip = FALSE)
```

Note that for daily data you need to give a well formed date string in the form of "YYYY-MM-DD".

## Interact with the archive and prism data

You can view all the prism data you have downloaded with a simple command: `prism_archive_ls()`.  This function gives a list of folder names, i.e., prism data (`pd`).  All the functions in the prism package work off of one or more of these folder names (`pd`).

```{r listingFiles, out.lines=10}
## Truncated to keep file list short
prism_archive_ls()
```

While prism functions use this folder format, other files may need an absolute path (e.g. the `raster` package). The `pd_to_file()` function conveniently returns the absolute path.  Alternatively, you may want to see what the normal name for the product is (not the file name), and we can get that with the `pd_get_name()` function.

```{r moreListing, out.lines=5}
## Truncated to keep file list short
pd_to_file(prism_archive_ls())

pd_get_name(prism_archive_ls())
```

Finally, `prism_archive_subset()` is a convenient way to search for specific parameters, time steps, days, months, years, or ranges of days, months, years. 

```{r}
# we know we have downloaded June 2013 daily data, so lets search for those 
prism_archive_subset("tmean", "daily", mon = 6)

# or we can look for days between June 7 and June 10
prism_archive_subset(
  "tmean", "daily", minDate = "2013-06-07", maxDate = "2013-06-10"
)
```

### Raster plots

You can easily make a quick plot of your data using the output of `prism_archive_ls()` or `prism_archive_subset()` with `pd_image()`. 

```{r quick_plot,fig.height=5,fig.width=7}
# Plot the January 30-year average temperatures
# grab only the first value, just in case multiple values are returned
jmean <- prism_archive_subset(
  "tmean", "monthly normals", mon = 1, resolution = "4km"
)
pd_image(jmean)
```

It is easy to load the prism data with the raster package. This time we will look at January temperature anomalies.  To do this we will examine the difference between January 2013 and the January 30-year normals. Conveniently, we already downloaded these data. We just need to grab them out of our archive.

```{r raster_math,fig.height=5,fig.width=7}
library(raster)
# knowing the name of the files you are after allows you to find them in the 
# list of all files that exist
# jnorm_name <- "PRISM_tmean_30yr_normal_4kmM2_01_bil"
# j2013_name <- "PRISM_tmean_stable_4kmM3_201301_bil"
# but we will use prism_archive_subset() to find the files we need

jnorm <- prism_archive_subset(
  "tmean", "monthly normals", mon = 1, resolution = "4km"
)
j2013 <- prism_archive_subset("tmean", "monthly", years = 2013, mon = 1)

# raster needs a full path, not the "short" prism data name
jnorm <- pd_to_file(jnorm)
j2013 <- pd_to_file(j2013)

## Now we'll load the rasters.
jnorm_rast <- raster(jnorm)
j2013_rast <- raster(j2013)

# Now we can do simple subtraction to get the anomaly by subtracting 2014 
# from the 30 year normal map
anomCalc <- function(x, y) {
  return(x - y)
}

anom_rast <- raster::overlay(j2013_rast,jnorm_rast,fun = anomCalc)

plot(anom_rast)
```

The plot shows that January 2013 was warmer than the average over the last 30 years. It also shows how easy it is to use the raster library to work with prism data. The package provides a simple framework to work with a large number of rasters that you can easily download and visualize or use with other data sets.

### Single grid cell plot

You can also visualize a single point across multiple prism data files (slice) using `pd_plot_slice()`. This procedure will take a set of rasters, create a "`raster::stack`", extract data at a point, and then create a ggplot2 object.

Let's now make a plot of January temperatures in Boulder between 1982 and 2014. First we'll grab all the data from the US (downloaded in the previous step), and then give our function a point to get data from. The point must be a vector in the form of longitude, latitude. Because `pd_plot_slice()` returns a gg object, it can be combined with other ggplot functions.

```{r plot_Boulder,fig.height=5,fig.width=7, results=FALSE}
library(ggplot2)
# data already exist in the prism dl dir
boulder <- c(-105.2797, 40.0176)

# prism_archive_subset() will return prism data that matches the specified 
# variable, time step, years, months, days, etc.
to_slice <- prism_archive_subset("tmean", "monthly", mon = 1)
p <- pd_plot_slice(to_slice, boulder)

# add a linear average and title
p + 
  stat_smooth(method="lm", se = FALSE) + 
  theme_bw() + 
  ggtitle("Average January temperature in Boulder, CO 1982-2014")
```

### leaflet map

Finally, the prism data are in a form that can be used with leaflet maps (with the help of the raster package). The [leaflet package](https://CRAN.R-project.org/package=leaflet) allows you to easily make JavaScript maps using the [leaflet](https://leafletjs.com/) mapping framework using prism data.  These can easily be hosted on websites like [Rpubs](https://rpubs.com/) or your own site.  Here is a simple example of plotting the [30-year normal for annual temperature](https://rpubs.com/DistribEcology/122453). If you run this code you will have an interactive map, instead of just the screen shot shown here. 

```{r leaflet,eval=F}
library(leaflet)
library(raster)
library(prism)

# 30-year normal average temperature have already been downloaded for 
norm <- prism_archive_subset(
  "tmean", "annual normals", resolution = "4km"
)
rast <- raster(pd_to_file(norm))

# Create color palette and plot
pal <- colorNumeric(
  c("#0000FF", "#FFFF00", "#FF0000"), 
  values(rast),
  na.color = "transparent"
)

leaflet() %>% 
  addTiles(
    urlTemplate = 'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'
  ) %>% 
  addRasterImage(rast, colors = pal, opacity=.65) %>% 
  addLegend(pal = pal, values = values(rast), title = "Deg C")
```

[![leaflet example figure](leaflet_example.png)

---
title: "Example prism Analysis"
author: "Alan Butler"
date: "4/23/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

This document provides an example of how to use the prism package, combined with the raster package, to compute an example analyis. In this example, we will download monthly prism data, clip it to the Upper Colorado River Basin area, spatially average it, and then create water year totals. 

```{r, message = FALSE, warning=FALSE}
library(prism)
library(raster)
library(rgdal)
library(sp)
library(dplyr)
library(stringr)

prism_set_dl_dir("./prism_example")
```

## Download the Data

We will download all average temperature data for 2017-2020 so that we can have full water year (October - September) data for water years 2018-2020. 

```{r, eval=FALSE}
get_prism_monthlys("tmean", years = 2017:2020, mon = 1:12)
```
Then, create a raster stack of the data, using only October 2017 - September 2020. Need two partial years of data, so have to call `prism_archive_subset()` multiple times. 

```{r}
ond2017 <- prism_archive_subset("tmean", "monthly", years = 2017, mon = 10:12)
full_yrs <- prism_archive_subset("tmean", "monthly", years = 2018:2019, mon = 1:12)
js2019 <- prism_archive_subset("tmean", "monthly", years = 2020, mon = 1:9)

ps <- pd_stack(c(ond2017, full_yrs, js2019))
```

## Crop to CRB

Currently, each month of data is for CONUS. We want all of these months to only be for the Colorado River Basin. First, read in a shapefile of the CRB:

```{r}
crb <- readOGR("data/LC_UC_Basin 2011.shp")[2,] %>% # UB only but in UTM
  # transform to the same projection that the raster data is using
  spTransform("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
```
Then, we need a function to mask and crop a single layer (single month of prism data). Then, we can apply that function to every layer.

```{r}
mask_and_crop <- function(i, rasterStack, basin, progress_bar) {
  zz <- mask(crop(rasterStack[[i]], extent(basin)), basin)
  setTxtProgressBar(progress_bar, i)
  zz
}

# initialize the progress bar
pb <- txtProgressBar(min = 0, max = nlayers(ps) , style = 3)


ps_crop <- lapply(
    seq_len(nlayers(ps)), 
    mask_and_crop, 
    rasterStack = ps, 
    basin = crb,
    progress_bar = pb
  )

# this returns a list, and we want it to be a raster stack again
# convert the list to a stack
ps_crop <- raster::stack(ps_crop)
```

Here's a plot of the first month of data before and after masking and cropping:
```{r, echo = FALSE}
par(mfrow = c(2, 1))
plot(ps[[1]])
title("October 2017 tmean - CONUS")

plot(ps_crop[[1]])
title("October 2017 - tmean CRB")
```

## Spatial Stats

Finally, we want to compute a spatial average to get a CRB basin-wide average monthly temperature:

```{r}
crb_avg <- cellStats(ps_crop, "mean")
# this results in a named vector. We know that the PRISM naming convention uses
# PRISM_var_stable_4kmM3_yyyymm_bil, where var is tmean in this example, 
# so we can parse this to create a data frame with year, month, and value

df <- data.frame(ts = names(crb_avg), value = unname(crb_avg)) %>%
    # we want just the yyyymm_bil part, and then remove _bil
    mutate(
      ts = str_extract(ts, "\\d{6}_bil"),
      ts = str_remove(ts, "_bil"),
      # then get yyyy as year and mm as month,
      yrs = substr(ts, 1, 4),
      mm = substr(ts, 5, 6)
    )
```

And now you have a data frame of monthly average temperatures in the CRB:

```{r}
head(df)
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd_get.R
\name{pd_get_name}
\alias{pd_get_name}
\alias{pd_get_date}
\alias{pd_get_type}
\alias{prism_md}
\alias{pd_to_file}
\title{Perform action on "prism data"}
\usage{
pd_get_name(pd)

pd_get_date(pd)

pd_get_type(pd)

prism_md(f, returnDate = FALSE)

pd_to_file(pd)
}
\arguments{
\item{pd}{prism data character vector.}

\item{f}{1 or more prism directories name or .bil files.}

\item{returnDate}{TRUE or FALSE. If TRUE, an ISO date is returned.  By
default years will come back with YYYY-01-01 and months as YYYY-MM-01}
}
\value{
\code{pd_get_name()} and \code{pd_get_date()} return a character vector of
names/dates.
}
\description{
"prism data", i.e., \code{pd} are the folder names returned by
\code{\link[=prism_archive_ls]{prism_archive_ls()}} or \code{\link[=prism_archive_subset]{prism_archive_subset()}}. These functions get the
name or date from these data, or convert these data to a file name.

\code{pd_get_date()} extracts the date from the prism data.
Date is returned in yyyy-mm-dd format. For monthly data, dd is 01 and
for annual data mm is also 01. For normals, an empty character is returned.

\code{pd_get_type()} parses the variable from the prism data.

\code{prism_md()} is a deprecated function that has been replaced with
\code{pd_get_name()} and \code{pd_get_date()}

\code{pd_to_file()} converts prism data  to a fully specified .bil file, i.e., the
full path to the file in the prism archive. A warning is posted if the
file does not exist in the local prism archive.
}
\details{
\code{pd_get_name()} extracts a long, human readable name from the prism
data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd_stack.R
\name{pd_stack}
\alias{pd_stack}
\alias{prism_stack}
\title{Stack prism data}
\usage{
pd_stack(pd)

prism_stack(prismfile)
}
\arguments{
\item{pd, prismfile}{A vector of prism data returned by \code{\link[=prism_archive_ls]{prism_archive_ls()}}
or \code{\link[=prism_archive_subset]{prism_archive_subset()}}.}
}
\description{
\code{pd_stack()} creates a raster stack from prism data. It is up to the user to
ensure that \code{pd} is of the expected variable and temporal period, i.e., the
function does no checking and will stack data with different variables or
temporal periods.

\code{prism_stack()} is the deprecated version of \code{pd_stack()}.
}
\examples{
\dontrun{
get_prism_dailys(
  type="tmean", 
  minDate = "2013-06-01", 
  maxDate = "2013-06-14", 
  keepZip = FALSE
)
# get a raster stack of June 1-14 daily tmean
mystack <- prism_stack(prism_archive_subset(
  "tmean", 
  minDate = "2013-06-01", 
  maxDate = "2013-06-14"
))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prism_archive_clean.R
\name{prism_archive_clean}
\alias{prism_archive_clean}
\alias{del_early_prov}
\title{Clean the prism data by removing early and provisional data}
\usage{
prism_archive_clean(
  type,
  temp_period,
  years = NULL,
  mon = NULL,
  minDate = NULL,
  maxDate = NULL,
  dates = NULL
)

del_early_prov(type, minDate = NULL, maxDate = NULL, dates = NULL)
}
\arguments{
\item{type}{The type of data you want to subset. Must be "ppt", "tmean",
"tmin", "tmax", "tdmean", "vpdmin", or "vpdmax".}

\item{temp_period}{The temporal period to subset. Must be "annual",
"monthly", "daily", "monthly normals", or "annual normals".}

\item{years}{Valid numeric year, or vector of years.}

\item{mon}{Valid numeric month, or vector of months.}

\item{minDate}{Date to start subsetting daily data. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{maxDate}{Date to end subsetting daily data.  Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{dates}{A vector of daily dates to subset. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}
}
\value{
Invisibly returns vector of all deleted folders.
}
\description{
\code{prism_archive_clean()} 'cleans' the prism download data by removing early
and/or provisional data if newer (provisional or stable) data also exist
for the same variable and temporal period. Stable data are newer than
provisional data that are newer than early data; only the newest data are
kept when the "clean" is performed.

\code{del_early_prov()} is a deprecated version of \code{prism_archive_clean()} that
only works for daily data, and does not prompt the user to confirm which
folders should be removed.
}
\details{
\code{prism_archive_clean()} prompts the user to verify the folders that should be
removed when R is running in interactive mode. Otherwise, all data that are
identified to be older than the newest available data are removed.

Daily data are considered "early" for the current month. The previous six
months are provisional data. After six months data are considered stable.
Thus early data only exist for daily data, while there can be monthly (and
presumably yearly) provisional data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd_get_md.R
\name{pd_get_md}
\alias{pd_get_md}
\title{Get prism metadata}
\usage{
pd_get_md(pd)
}
\arguments{
\item{pd}{prism data character vector.}
}
\value{
data.frame containing metadata for all specified prism data.
}
\description{
Retrieves prism metadata from the specified prism data. "prism data", i.e.,
\code{pd} are the folder names returned by \code{\link[=prism_archive_ls]{prism_archive_ls()}} or
\code{\link[=prism_archive_subset]{prism_archive_subset()}}. These functions get the name or date from these
data, or convert these data to a file name. A warning is provided if the
specified prism data do not exist in the archive.
}
\details{
The metadata includes the following variables from the .info.txt file for
daily, monthly, and annual data:
\itemize{
\item PRISM_DATASET_FILENAME
\item PRISM_DATASET_CREATE_DATE
\item PRISM_DATASET_TYPE
\item PRISM_DATASET_VERSION
\item PRISM_CODE_VERSION
\item PRISM_DATASET_REMARKS
}

Additionally, two local variables are added identifying where the file is
located on the local system:
\itemize{
\item file_path
\item folder_path
}

The annual and monthly normals data includes different keys in
the .info.txt, so they are renamed to be the same as those found in the
other temporal data. The keys/variables are renamed as follows:
\itemize{
\item PRISM_FILENAME --> PRISM_DATASET_FILENAME
\item PRISM_CREATE_DATE --> PRISM_DATASET_CREATE_DATE
\item PRISM_DATASET --> PRISM_DATASET_TYPE
\item PRISM_VERSION --> PRISM_CODE_VERSION
\item PRISM_REMARKS --> PRISM_DATASET_REMARKS
}

Additionally, the normals does not include PRISM_DATASET_VERSION, so that
variable is added with \code{NA} values.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd_get_station_md.R
\name{pd_get_station_md}
\alias{pd_get_station_md}
\alias{get_prism_station_md}
\title{Extract prism station metadata}
\usage{
pd_get_station_md(pd)

get_prism_station_md(type, minDate = NULL, maxDate = NULL, dates = NULL)
}
\arguments{
\item{pd}{prism data character vector.}

\item{type}{The type of data you want to subset. Must be "ppt", "tmean",
"tmin", "tmax", "tdmean", "vpdmin", or "vpdmax".}

\item{minDate}{Date to start subsetting daily data. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{maxDate}{Date to end subsetting daily data.  Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{dates}{A vector of daily dates to subset. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}
}
\value{
A \code{tbl_df} containing metadata on the stations used for the specified
day and variable. The data frame contains the following columns:
"date", "prism_data", "type", "station", "name", "longitude",
"latitude", "elevation", "network", "stnid"

The "date" column is a character representation of the data. Monthly and
annual data are given first day of month, and first month of year for
reporting here. Monthly and annual normals are empty strings.
}
\description{
\code{pd_get_station_md()} extracts prism metadata on the stations used to
generate the prism data. \strong{The data must already be downloaded
and available in the prism download folder.} "prism data", i.e., \code{pd} are
the folder names returned by \code{\link[=prism_archive_ls]{prism_archive_ls()}} or
\code{\link[=prism_archive_subset]{prism_archive_subset()}}.

\code{get_prism_station_md()} is a deprecated version of
\code{pd_get_station_md()} that only works with daily prism data.
}
\details{
Note that station metadata does not exist for "tmean" type or for any
"annual" temporal periods.

See \code{\link[=prism_archive_subset]{prism_archive_subset()}} for further details
on specifying ranges of dates for different temporal periods.
}
\examples{
\dontrun{
# download and then get meta data for January 1, 2010 precipitation
get_prism_dailys("ppt", dates = "2010-01-01")
pd <- prism_archive_subset("ppt", "daily", dates = "2010-01-01")

# will warn that 2010-01-02 is not found:
pd_get_station_md(pd)
}
  
}
\seealso{
\code{\link[=prism_archive_subset]{prism_archive_subset()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd_image.R
\name{pd_image}
\alias{pd_image}
\alias{prism_image}
\title{Quick spatial image of prism data}
\usage{
pd_image(pd, col = "heat")

prism_image(prismfile, col = "heat")
}
\arguments{
\item{pd, prismfile}{the name of a single file to be plotted, this is most
easily found through \code{\link[=prism_archive_ls]{prism_archive_ls()}} or \code{\link[=prism_archive_subset]{prism_archive_subset()}}.}

\item{col}{the color pattern to use.  The default is heat, the other valid
option is "redblue".}
}
\value{
Invisibly returns \code{gg} object of the image.
}
\description{
\code{pd_image()} makes a spatial image plot of the specified prism
data (single variable and time step.). It is meant for rapid visualization,
but more detailed plots will require other methods.

\code{prism_image()} is the deprecated version of \code{pd_image()}.
}
\examples{
\dontrun{
get_prism_dailys(
  type = "tmean",
  minDate = "2013-06-01",
  maxDate = "2013-06-14",
  keepZip = FALSE
)

# get June 5th
pd <- prism_archive_subset("tmean", "daily", dates = "2013-06-05")

# and plot it
pd_image(pd)
}

}
\seealso{
\code{\link[=prism_archive_ls]{prism_archive_ls()}}, \code{\link[=prism_archive_subset]{prism_archive_subset()}},
\code{\link[ggplot2:geom_tile]{ggplot2::geom_raster()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{prism_check}
\alias{prism_check}
\title{Check if prism files exist}
\usage{
prism_check(prismfiles, lgl = FALSE, pre81_months = NULL)
}
\arguments{
\item{prismfiles}{a list of full prism file names ending in ".zip".}

\item{lgl}{\code{TRUE} returns a logical vector indicating those
not yet downloaded; \code{FALSE} returns the file names that are not yet
downloaded.}

\item{pre81_months}{Numeric vector of months that will be downloaded, if
downloading data before 1981. This is so that the existence of the data can
be correctly checked, as the file includes all monthly data for a given
year.}
}
\value{
a character vector of file names that are not yet downloaded
or a logical vector indication those not yet downloaded.
}
\description{
Helper function to check if files already exist in the prism download
directory. Determines if files have \strong{not} been downloaded yet, i.e.,
returns \code{TRUE} if they do not exist.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prism.R
\docType{package}
\name{prism-package}
\alias{prism}
\alias{prism-package}
\title{prism: Access Data from the Oregon State Prism Climate Project}
\description{
Allows users to access the Oregon State Prism climate data
    (<https://www.prism.oregonstate.edu/>). Using the web service API data
    can easily downloaded in bulk and loaded into R for spatial analysis.
    Some user friendly visualizations are also provided.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/prism/}
  \item \url{https://github.com/ropensci/prism}
  \item Report bugs at \url{https://github.com/ropensci/prism/issues}
}

}
\author{
\strong{Maintainer}: Alan Butler \email{rabutler@usbr.gov} [contributor]

Authors:
\itemize{
  \item Hart Edmund \email{Edmund.m.hart@gmail.com} [conceptor]
  \item Kendon Bell \email{kmb56@berkeley.edu}
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd_plot_slice.R
\name{pd_plot_slice}
\alias{pd_plot_slice}
\alias{prism_slice}
\title{Plot a slice of a raster stack}
\usage{
pd_plot_slice(pd, location)

prism_slice(location, prismfile)
}
\arguments{
\item{pd, prismfile}{a vector of output from \code{\link[=prism_archive_ls]{prism_archive_ls()}} or
\code{\link[=prism_archive_subset]{prism_archive_subset()}} giving a list of prism files to extract data from
and plot. The latter is preferred as it will help ensure the prism data
are from the same variable and temporal period.}

\item{location}{a vector of a single location in the form of long,lat}
}
\value{
A \code{gg} object of the plot for the requested \code{location}.
}
\description{
\code{pd_plot_slice()} plots a slice of data at a single point location from the
specified prism data.

\code{prism_slice()} is the deprecated version of \code{pd_plot_slice()}.
}
\details{
The user should ensure the prism data comes from a continuous data
set and is made up of the same temporal period. Otherwise the plot will look
erratic and incorrect.
}
\examples{
\dontrun{
### Assumes you have a clean prism directory
get_prism_dailys(
  type="tmean", 
  minDate = "2013-06-01", 
  maxDate = "2013-06-14",
  keepZip = FALSE
)
p <- pd_plot_slice(
  prism_archive_subset("tmean", "daily", year = 2020), 
  c(-73.2119,44.4758)
)
print(p)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_prism_annual.R, R/get_prism_dailys.R,
%   R/get_prism_monthlys.R, R/get_prism_normals.R
\name{get_prism_annual}
\alias{get_prism_annual}
\alias{get_prism_dailys}
\alias{get_prism_monthlys}
\alias{get_prism_normals}
\title{Download prism data}
\usage{
get_prism_annual(type, years = NULL, keepZip = TRUE, keep_pre81_months = FALSE)

get_prism_dailys(
  type,
  minDate = NULL,
  maxDate = NULL,
  dates = NULL,
  keepZip = TRUE,
  check = "httr"
)

get_prism_monthlys(
  type,
  years = NULL,
  mon = NULL,
  keepZip = TRUE,
  keep_pre81_months = TRUE
)

get_prism_normals(type, resolution, mon = NULL, annual = FALSE, keepZip = TRUE)
}
\arguments{
\item{type}{The type of data to download. Must be "ppt", "tmean", "tmin",
"tmax", "tdmean", "vpdmin", or "vpdmax". Note that \code{tmean ==  mean(tmin, tmax)}.}

\item{years}{a valid numeric year, or vector of years, to download data for.
If no month is specified, year averages for that year will be downloaded.}

\item{keepZip}{if \code{TRUE}, leave the downloaded zip files in your
'prism.path', if \code{FALSE}, they will be deleted.}

\item{keep_pre81_months}{The pre-1981 data includes all monthly data and the
annual data for the specified year. If you need annual and monthly data it
is advantageous to keep all the monthly data when downloading the annual
data so you don't have to download the zip file again. When downloading
annual data, this defaults to \code{FALSE}. When downloading monthly data, this
defaults to \code{TRUE}.}

\item{minDate}{Date to start downloading daily data. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{maxDate}{Date to end downloading daily data.  Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{dates}{A vector of dates to download daily data. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{check}{One of "httr" or "internal". See details.}

\item{mon}{a valid numeric month, or vector of months.}

\item{resolution}{The spatial resolution of the data, must be either "4km"
or "800m".}

\item{annual}{if \code{TRUE} download annual normals.}
}
\description{
Download grid cell data from the
\href{https://prism.nacse.org/}{prism project}. Temperature (min, max,
and mean), mean dewpoint temperature, precipitation, and vapor pressure
deficit (min and max) can be downloaded for annual (\code{get_prism_annual()}),
monthly (\code{get_prism_monthlys()}), daily (\code{get_prism_dailys()}), and 30-year
averages (\code{get_prism_normals()}). Data are at 4km resolution, except for the
normals which can also be downloaded at 800m resolution.

Download data from the prism project for 30 year normals at 4km
or 800m grid cell resolution for precipitation, mean, min and max
temperature
}
\details{
A valid download directory must exist before downloading any prism data. This
can be set using \code{\link[=prism_set_dl_dir]{prism_set_dl_dir()}} and can be verified using
\code{\link[=prism_check_dl_dir]{prism_check_dl_dir()}}.

For the \code{check} parameter, "httr", the default, checks the file name using
the web service, and downloads if that file name is not in the file system.
"internal" (much faster) only attempts to download layers that are not
already in the file system as stable. "internal" should be used with caution
as it is not robust to changes in version or file names.
}
\section{Annual and Monthly}{


Annual and monthly prism data are available from 1891 to present. For
1891-1980 data, monthly and annual data are grouped together in one download
file; \code{keep_pre81_months} determines if the other months/yearly data are kept
after the download.  Data will be downloaded for all specified months (\code{mon})
in all the \code{years} in the supplied vectors.
}

\section{Daily}{


Daily prism data are available beginning on January 1, 1981. To download the
daily data, dates must be in the proper format or downloading will not work
properly. Dates can be specified using either a date range via \code{minDate} and
\code{maxDate}, or a vector of \code{dates}, but not both.
}

\section{Normals}{


30-year normals are currently computed using 1981-2010 and are available at
4km and 800m resolution. See
\url{https://prism.nacse.org/normals/}.
If \code{mon} is specified and \code{annual} is \code{TRUE}, then monthly and annual normal
data will be downloaded.
}

\examples{
\dontrun{
# Get all annual average temperature data from 1990 to 2000
get_prism_annual(type = "tmean", year = 1990:2000, keepZip = FALSE)
}

\dontrun{
# get daily average temperature data for June 1 - 14, 2013
get_prism_dailys(
  type = "tmean", 
  minDate = "2013-06-01", 
  maxDate = "2013-06-14", 
  keepZip=FALSE
)

# get precipitation datat for June 1, 2013
get_prism_dailys(type = "ppt", dates = "2013/06/01", keepZip = FALSE)

# get average temperature for three specific days
get_prism_dailys(
  type="tmean", 
  dates = as.Date("2013-06-01", "2013-06-14", "2014-06-30"), 
  keepZip=FALSE
)

# will fail:
get_prism_dailys(
  type = "ppt", 
  minDate = "2013-06-01", 
  dates = "2013-06-14", 
  keepZip = FALSE
)

get_prism_dailys(
  type = "ppt", 
  minDate = "2013-06-01", 
  keepZip=FALSE
)
}

\dontrun{
# Get all the precipitation data for January from 1990 to 2000
get_prism_monthlys(type = "ppt", years = 1990:2000, mon = 1, keepZip = FALSE)

# Get January-December 2005 monthly precipitation
get_prism_monthlys(type = "ppt", years = 2005, mon = 1:12, keepZip = FALSE)
}

\dontrun{
# Get 30 year normal values for January rainfall
get_prism_normals(type = "ppt", resolution = "4km", mon = 1, keepZip = FALSE)

# Get monthly (every month) and annual 30-year normals for mean temperature
get_prism_normals(
  type = "tmean", 
  resolution = "800m", 
  mon = 1:12, 
  annual = TRUE,
  keepZip = FALSE
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prism_archive_subset.R
\name{prism_archive_subset}
\alias{prism_archive_subset}
\title{Subsets PRISM folders on the disk}
\usage{
prism_archive_subset(
  type,
  temp_period,
  years = NULL,
  mon = NULL,
  minDate = NULL,
  maxDate = NULL,
  dates = NULL,
  resolution = NULL
)
}
\arguments{
\item{type}{The type of data you want to subset. Must be "ppt", "tmean",
"tmin", "tmax", "tdmean", "vpdmin", or "vpdmax".}

\item{temp_period}{The temporal period to subset. Must be "annual",
"monthly", "daily", "monthly normals", or "annual normals".}

\item{years}{Valid numeric year, or vector of years.}

\item{mon}{Valid numeric month, or vector of months.}

\item{minDate}{Date to start subsetting daily data. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{maxDate}{Date to end subsetting daily data.  Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{dates}{A vector of daily dates to subset. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{resolution}{The spatial resolution of the data, must be either "4km" or
"800m". Should only be specified for \code{temp_period} of "normals".}
}
\value{
A character vector of the folders that meet the type and temporal
period specified. \code{character(0)} is returned if no folders are found that
meet the specifications.
}
\description{
\code{prism_archive_subset()} subsets the PRISM folders stored on disk by type,
temporal period, and date. It looks through all of the PRISM data that have
been downloaded in the prism archive (\code{\link[=prism_get_dl_dir]{prism_get_dl_dir()}}) and returns the
subset based on the specified \code{type}, \code{temp_period}, and dates.
}
\details{
\code{temp_period} must be specified so the function can distinguish between
wanting annual data or wanting monthly data for a specified year. For example
\code{prism_archive_subset("tmean", "annual", years = 2012)} would provide only
one folder: the annual average temperature for 2012. However,
\code{prism_archive_subset("tmean", "monthly", years = 2012)} would provide 12
folders: each monthly tmean folder for 2012.

\code{temp_period}, \code{years}, and \code{mon} can be combined in various different ways
to obtain different groupings of data. \code{years}, \code{mon}, and the daily
specifiers (\code{minDate}/\code{maxDate} or \code{dates}) are optional. Not specifying any
of those would result in getting all annual, monthly, or daily data.

\code{minDate}/\code{maxDate} or \code{dates} should only be specified for a \code{temp_period}
of "daily". Additionally, only \code{dates}, or \code{minDate} and \code{maxDate}, should be
specified, but all three should not be specified. Nor should the daily
arguments be combined with \code{years} and/or \code{mon}. For example, if daily
folders are desired, then specify \code{years} and/or \code{mon} to get all days for
those years and months \strong{or} specify the specific dates using
\code{minDate}/\code{maxDate} or \code{dates}
}
\examples{
\dontrun{
# get all annual tmin
prism_archive_subset("tmin", "annual")
# get only 2000-2015 annual tmin
prism_subset_folder("tmin", "annual", years = 2000-2015)

# get monthly precipitation for 2000-2010
prism_archive_subset("ppt", "monthly", years = 2000-2010)
# get only June-August monthly precip data for 2000-2010
prism_archive_subset("ppt", "monthly", years = 2000-2010, mon = 6:8)

# get all daily tmax for July-August in 2010
prism_archive_subset("tmax", "daily", years = 2010, mon = 7:8)
# same as:
prism_archive_subset(
  "tmax", 
  "daily", 
  minDate = "2010-07-01", 
  maxDate = "2010-08-31"
)

# get the 4km 30-year average precip for January and February
prism_archive_subset("ppt", "monthly normals", mon = 1:2, resolution = "4km")
}

}
\seealso{
\code{\link[=prism_archive_ls]{prism_archive_ls()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prism_set_dl_dir.R
\name{prism_set_dl_dir}
\alias{prism_set_dl_dir}
\alias{prism_get_dl_dir}
\alias{prism_check_dl_dir}
\alias{path_check}
\title{Set, check, and get prism download directory}
\usage{
prism_set_dl_dir(path, create = TRUE)

prism_get_dl_dir()

prism_check_dl_dir()

path_check()
}
\arguments{
\item{path}{The path that prism data will be unzipped into.}

\item{create}{Boolean that determines if the \code{path} will be created if it
does not already exist.}
}
\description{
\code{prism_set_dl_dir()} sets the directory that downloaded prism data will be
saved to. The prism download directory is saved in the "prism.path" option.

\code{prism_get_dl_dir()} gets the folder that prism data will be saved to. It is
a wrapper around \code{getOption("prism.path")} so the user does not have to
remember the option name.

\code{prism_check_dl_dir()} checks that prism download folder has been set. If it
has not been set, and in interactive mode, then prompt user to specify the
download location. If not in interactive mode, and it has not been set, then
set to "~/prismtmp".

\code{path_check()} is a deprecated version of \code{prism_check_dl_dir()}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prism_archive_ls.R
\name{prism_archive_ls}
\alias{prism_archive_ls}
\alias{ls_prism_data}
\title{List available prism data}
\usage{
prism_archive_ls()

ls_prism_data(absPath = FALSE, name = FALSE)
}
\arguments{
\item{absPath}{TRUE if you want to return the absolute path.}

\item{name}{TRUE if you want file names and titles of data products.}
}
\value{
\code{prism_archive_ls()} returns a character vector.

\code{ls_prism_data()} returns a data frame. It can have 1-3 columns, but
always has the \code{files} column. \code{abs_path} and \code{product_name} columns are
added if \code{absPath} and \code{name} are \code{TRUE}, respectively.
}
\description{
\code{prism_archive_ls()} lists all available prism data (all variables and all
temporal periods) that are available in the local archive, i.e., they
have already been downloaded and are available in \code{\link[=prism_get_dl_dir]{prism_get_dl_dir()}}.
\code{\link[=prism_archive_subset]{prism_archive_subset()}} can be used to subset the archive based on specified
variables and temporal periods.

\code{ls_prism_data()} is a deprecated version of \code{prism_data_ls()}.
}
\details{
\code{prism_archive_ls()} only returns the values found in the \code{files} column as
returned by \code{ls_prism_data()}. To replicate
the behavior of \code{ls_prism_data()}, use \code{\link[=pd_get_name]{pd_get_name()}} and
\code{\link[=pd_to_file]{pd_to_file()}} with the output of \code{prism_archive_ls()}
}
\examples{
\dontrun{
# Get prism data names, used in many other prism* functions 
get_prism_dailys(
  type="tmean", 
  minDate = "2013-06-01", 
  maxDate = "2013-06-14", 
  keepZip = FALSE
)
prism_archive_ls()
}

}
\seealso{
\code{\link[=prism_archive_subset]{prism_archive_subset()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prism_archive_verify.R
\name{prism_archive_verify}
\alias{prism_archive_verify}
\alias{check_corrupt}
\title{Check the integrity of downloaded PRISM data}
\usage{
prism_archive_verify(
  type,
  temp_period,
  years = NULL,
  mon = NULL,
  minDate = NULL,
  maxDate = NULL,
  dates = NULL,
  download_corrupt = TRUE,
  keepZip = TRUE
)

check_corrupt(type, minDate = NULL, maxDate = NULL, dates = NULL)
}
\arguments{
\item{type}{The type of data you want to subset. Must be "ppt", "tmean",
"tmin", "tmax", "tdmean", "vpdmin", or "vpdmax".}

\item{temp_period}{The temporal period to subset. Must be "annual",
"monthly", "daily", "monthly normals", or "annual normals".}

\item{years}{Valid numeric year, or vector of years.}

\item{mon}{Valid numeric month, or vector of months.}

\item{minDate}{Date to start subsetting daily data. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{maxDate}{Date to end subsetting daily data.  Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{dates}{A vector of daily dates to subset. Must be specified in
a valid iso-8601 (e.g. YYYY-MM-DD) format. May be provided as either a
character or \link[base:Dates]{base::Date} class.}

\item{download_corrupt}{If \code{TRUE}, then any unreadable prism data are
automatically re-downloaded.}

\item{keepZip}{If \code{TRUE}, leave the downloaded zip files in your
'prism.path', if \code{FALSE}, they will be deleted.}
}
\value{
\code{prism_archive_verify()} returns \code{TRUE} if all data are readable.
Any prism data that are not readable are returned (folder names), whether
they are re-downloaded or not.

\code{check_corrupt()} returns \code{logical} indicating whether the process
succeeded.
}
\description{
\code{prism_archive_verify()} checks the data in the prism archive to ensure it
is valid, or at least can be read into R, i.e., it is not corrupt. The
prism variable type, time period, etc. is specified the same as for
\code{\link[=prism_archive_subset]{prism_archive_subset()}}. Any files that are not readable can automatically
be re-downloaded.

\code{check_corrupt()} is the deprecated version of
\code{prism_archive_verify()}
}
\details{
Under the hood, it uses \code{raster::stack()} and then \code{raster::rasterToPoints()}
to determine if the bil files are readable. If both those files are able
to successfully read the files, they are assumed to be valid/readable.
}
