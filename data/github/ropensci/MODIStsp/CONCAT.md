
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![](https://www.r-pkg.org/badges/version-ago/MODIStsp)](https://cran.rstudio.com/web/packages/MODIStsp/index.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.290683.svg)](https://doi.org/10.5281/zenodo.290683)
[![Downloads](https://cranlogs.r-pkg.org/badges/MODIStsp?color=orange)](https://cran.rstudio.com/web/packages/MODIStsp/index.html)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Travis-CI Build
Status](https://travis-ci.org/ropensci/MODIStsp.svg?branch=master)](https://travis-ci.org/ropensci/MODIStsp)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/MODIStsp/master.svg)](https://codecov.io/github/ropensci/MODIStsp?branch=master)

# <i class="fa fa-globe" aria-hidden="true"></i> MODIStsp <img src="man/figures/logo.png" width="100" height="100" align="right"/>

<!-- # MODIStsp <img src='man/figures/logo.png' align="right" height="139" /> -->

**`{MODIStsp}`** is a `R` package devoted to automatizing the creation
of time series of raster images derived from MODIS Land Products data.

**`{MODIStsp}`** allows performing several preprocessing steps (e.g.,
download, mosaicing, reprojection, resize, data extraction) on MODIS
data available within a given time period. Users have the ability to
select which specific layers of the original MODIS HDF files they want
to process. They also can select which additional **Quality Indicators**
should be extracted from the aggregated MODIS Quality Assurance layers
and, in the case of Surface Reflectance products, which **Spectral
Indexes** should be computed from the original reflectance bands.

All processing parameters can be easily selected with a **powerful and
user-friendly GUI**, although non-interactive execution exploiting a
previously created Options File is possible. Stand-alone execution
outside an `R` environment is also possible, allowing using scheduled
execution of MODIStsp to automatically update time series related to a
MODIS product and extent whenever a new image is available.

For each output layer, outputs are saved as **single-band raster** files
corresponding to each available acquisition date. **Virtual files**,
allowing accessing to the entire time series as a single file, can be
also created.

<a href="http://www.irea.cnr.it/en/">
<img src="man/figures/irea_logo.png" height="40" align="left" /></a>

<span style="font-style:italic;font-weight:bold;">`{MODIStsp}` was
developed by Lorenzo Busetto and Luigi Ranghetti, [Institute of Remote
Sensing of Environment](http://www.irea.cnr.it/en/) - National Research
Council - Italy (CNR-IREA). [It is dedicated to the memory of
Lorenzo](https://docs.ropensci.org/MODIStsp/articles/lorenzo.html).</span>

------------------------------------------------------------------------

## <i class="fa fa-newspaper-o" aria-hidden="true"></i> What’s New

-   29/10/2021 - `{MODIStsp}` (GitHub version 2.0.6.9000) supports
    products with version 061. Version 006 will remain the default
    product version until its decommission will be announced. Version
    061 can be specified through the argument `prod_version` of function
    `MODIStsp()` or by selecting it in the GUI.

-   10/12/2020 - `{MODIStsp}` was resubmitted to CRAN after the
    maintainer’s death. Now `{MODIStsp}` is dedicated to Lorenzo Busetto
    (<https://docs.ropensci.org/MODIStsp/articles/lorenzo>).

-   01/09/2020 - `{MODIStsp}` 2.0.0 is out. Replaces the old gWidgets
    GUI with a new one based on Shiny, enhances support for CLI usage
    and enhances support/provides bug fixing for datasets with multiple
    NoData values when applying scale/offset.

-   09/05/2020 - `{MODIStsp}` 1.4.0 is out. Switches to use of
    GDAL3/PROJ6 WKTs for projection representation and usage of `{sf}`
    for all internal work on vector data. Adds support for products
    MCD19A1 and MCD19A2 products.

-   07/06/2019 - `{MODIStsp}` 1.3.9 is out. Fixes a bug causing crashes
    on MOD14A1 product, adds support for product MCD12Q2 and removes
    support for no longer available version 5 of some products.

-   05/03/2019 - `{MODIStsp}` 1.3.8 is out. Fixes an issue causing
    incorrect application of scale/offset values on GDAL versions > 2.3
    (<https://github.com/ropensci/MODIStsp/issues/163>) and adds support
    for products `MOD21A1D.006 MOD21A1N.006 MOD21A2.006`.

-   29/11/2018 - We recently discovered a nasty bug in the computation
    of some custom spectral indices (those including additions /
    subtractions on reflectance values, such as in
    $\\frac{(b1\_{NIR}+0.1)}{b2\_{Red}}$ (with *ρ* being a reflectance).
    See
    [here](https://docs.ropensci.org/MODIStsp/articles/discovered_bug.html)
    for further details. The bug is now fixed on the GitHub version. A
    patched release will be made available on CRAN as soon as possible.

-   07/08/2018 - We are glad to report that `{MODIStsp}` is now included
    in the [rOpenSci](https://ropensci.org/about/) packages ecosystem.
    We thank reviewers Leah Wasser and Jeffrey Hanson for their valuable
    reviews, which helped us  
    further improving the package.

-   10/07/2018 - `{MODIStsp}` v. 1.3.6 is out (check out the [Release
    Notes](https://github.com/ropensci/MODIStsp/releases/tag/1.3.6) for
    further details).

-   20/06/2018 - `{MODIStsp}` v. 1.3.5 is out (check out the [Release
    Notes](https://github.com/ropensci/MODIStsp/releases/tag/v1.3.5) for
    further details).

-   11/04/2018 - Due to new NASA Policies the MODIS FTP servers were
    shut down starting, April 2, 2018. **FTP download is therefore no
    longer working** and will be removed in the next MODIStsp version.

-   11/04/2018 - [**Decommissioning of MODIS Version 5 Land Data
    Products**](https://lpdaac.usgs.gov/news/decommissioning-modis-version-51-land-cover-type-data-products-january-7-2019/).
    As per NASA notice above, MODIS v005 products are going to be
    decommissioned, and will soon be no longer available for download.
    Support for those products will be removed in the next MODIStsp
    version.

-   11/08/2017 - `{MODIStsp}` 1.3.3 was released today. It provides
    improvements in processing speed, as well as the usual bug fixes.
    See our [\<i class=“fa fa-newspaper-o aria-hidden=”true”></i>
    news](news/index.html) page for a detailed changelog.

------------------------------------------------------------------------

## <i class="fa fa-cog" aria-hidden="true"></i> Getting Started

-   **To install `{MODIStsp}`**, please follow instructions reported
    [here](articles/installation.html), both for
    [<i class="fa fa-windows" aria-hidden="true"></i>](articles/installation.html#installing-on-windows)
    ,
    [<i class="fa fa-linux" aria-hidden="true"></i>](articles/installation.html#installing-on-linux-systems)
    and
    [<i class="fa fa-apple" aria-hidden="true"></i>](articles/installation.html#installing-on-mac).

-   **`{MODIStsp}`** can be used either in [interactive
    mode](articles/interactive_execution.html) exploiting its
    user-friendly GUI, or in [non-interactive
    mode](articles/noninteractive_execution.html) from within `R`
    scripts.

-   The list of **currently supported MODIS products and versions** can
    be found [here](articles/products_list.html).

-   [Scheduled
    Processing](articles/noninteractive_execution.html#scheduled-processing)
    allows automatic updating of MODIS time series through scheduled
    jobs, both on
    [<i class="fa fa-windows" aria-hidden="true"></i>](articles/standalone_execution.html#on-windows)
    and
    [<i class="fa fa-linux" aria-hidden="true"></i>](articles/standalone_execution.html#on-linux).

-   Solutions to common **installation, data download and processing
    problems** can be found in our
    [<i class="fa fa-question-circle-o" aria-hidden="true"></i>
    faq](https://docs.ropensci.org/MODIStsp/articles/faq.html).

-   Please **report any issues** you may encounter in our [issues page
    on github
    <i class="fa fa-github-square" aria-hidden="true"></i>](https://github.com/ropensci/MODIStsp/issues).

------------------------------------------------------------------------

## <i class="fa fa-pencil" aria-hidden="true"></i>Citation

To cite MODIStsp in publications, please use:

L. Busetto, L. Ranghetti (2016) MODIStsp: An R package for automatic
preprocessing of MODIS Land Products time series, Computers &
Geosciences, Volume 97, Pages 40-48, ISSN 0098-3004,
<https://doi.org/10.1016/j.cageo.2016.08.020>, URL:
<https://github.com/ropensci/MODIStsp>.

A BibTeX entry for LaTeX users is:

      @Article{,
        author  = {Lorenzo Busetto and Luigi Ranghetti},
        title   = {MODIStsp: an R package for preprocessing of MODIS Land Products time series},
        journal = {Computers & Geosciences},
        year    = {2016},
        volume  = {97},
        pages   = {40-48},
        issn    = {0098-3004},
        doi     = {10.1016/j.cageo.2016.08.020},
        url     = {https://github.com/ropensci/MODIStsp},
      }



[![](https://www.r-pkg.org/badges/version-ago/MODIStsp)](https://cran.r-project.org/package=MODIStsp)
[![](https://cranlogs.r-pkg.org/badges/MODIStsp?color=orange)](https://cran.r-project.org/package=MODIStsp)
[![Travis-CI Build
Status](https://travis-ci.org/ropensci/MODIStsp.svg?branch=master)](https://travis-ci.org/ropensci/MODIStsp)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1972039.svg)](https://doi.org/10.5281/zenodo.1972039)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/MODIStsp/master.svg)](https://codecov.io/github/ropensci/MODIStsp?branch=master)
[![](https://badges.ropensci.org/184_status.svg)](https://github.com/ropensci/software-review/issues/184)

# MODIStsp <img src='man/figures/logo.png' align="right" height="139" />

[`{MODIStsp}`](https://docs.ropensci.org/MODIStsp/) is a `R` package
devoted to automatizing the creation of time series of rasters derived
from MODIS Land Products data. `{MODIStsp}` allows performing several
preprocessing steps (e.g., download, mosaicing, reprojection and resize)
on MODIS data available within a given time period. Users have the
ability to select which specific layers of the original MODIS HDF files
they want to process. They also can select which additional Quality
Indicators should be extracted from the aggregated MODIS Quality
Assurance layers and, in the case of Surface Reflectance products, which
Spectral Indexes should be computed from the original reflectance bands.
For each output layer, outputs are saved as single-band raster files
corresponding to each available acquisition date. Virtual files allowing
access to the entire time series as a single file can be also created.
All processing parameters can be easily selected with a user-friendly
GUI, although non-interactive execution exploiting a previously created
Options File is possible. Stand-alone execution outside an “R”
environment is also possible, allowing to use scheduled execution of
MODIStsp to automatically update time series related to a MODIS product
and extent whenever a new image is available.

<a href="http://www.irea.cnr.it/en/">
<img src="man/figures/irea_logo.png" height="40" align="left" /></a>

<span style="font-style:italic;font-weight:bold;">`{MODIStsp}` was
developed by Lorenzo Busetto and Luigi Ranghetti, [Institute of Remote
Sensing of Environment](http://www.irea.cnr.it/en/) - National Research
Council - Italy (CNR-IREA). [It is dedicated to the memory of
Lorenzo](https://docs.ropensci.org/MODIStsp/articles/lorenzo.html).</span>

## Citation

To cite `{MODIStsp}` please use:

L. Busetto, L. Ranghetti (2016) MODIStsp: An R package for automatic
preprocessing of MODIS Land Products time series, Computers &
Geosciences, Volume 97, Pages 40-48, ISSN 0098-3004,
<https://doi.org/10.1016/j.cageo.2016.08.020>, URL:
<https://github.com/ropensci/MODIStsp>.

## Website

For more information, documentation and examples of use, **see also the
`{MODIStsp}` website at
[docs.ropensci.org/MODIStsp](https://docs.ropensci.org/MODIStsp/)**.

## Important News

-   29/10/2021 - `{MODIStsp}` (GitHub version 2.0.6.9000) supports
    products with version 061. Version 006 will remain the default
    product version until its decommission will be announced. Version
    061 can be specified through the argument `prod_version` of function
    `MODIStsp()` or by selecting it in the GUI.

-   10/12/2020 - `{MODIStsp}` was resubmitted to CRAN after the
    maintainer’s death. Now `{MODIStsp}` is dedicated to Lorenzo Busetto
    (<https://docs.ropensci.org/MODIStsp/articles/lorenzo>).

-   01/09/2020 - `{MODIStsp}` 2.0.0 is out. Provides a new GUI interface
    based on Shiny, getting rid of the archived dependencies on
    gWidgets/gWidgetsRGtk2. Also provides much easier usage from the
    CLI, by allowing to set all processing arguments also from the CLI.
    **Note:** due to the introduced changes, options files created with
    previous versions of `{MODIStsp}` will no longer work. Also,
    processing scripts using `{MODIStsp}` may need to be slightly
    adapted.

-   09/05/2020 - `{MODIStsp}` 1.4.0 is out. Switches to use of
    GDAL3/PROJ6 WKTs for projection representation and usage of `{sf}`
    for all internal work on vector data. Adds support for products
    MCD19A1 and MCD19A2 products.

-   07/06/2019 - `{MODIStsp}` 1.3.9 is out. Fixes a bug causing crashes
    on MOD14A1 product, adds support for product MCD12Q2 and removes
    support for no longer available version 5 of some products.

-   05/03/2019 - `{MODIStsp}` 1.3.8 is out. Fixes an issue causing
    incorrect application of scale/offset values on GDAL versions > 2.3
    (<https://github.com/ropensci/MODIStsp/issues/163>) and adds support
    for products MOD21A1D.006, MOD21A1N.006 and MOD21A2.006.

-   29/11/2018 - We recently discovered a nasty bug in the computation
    of some custom spectral indices (those including additions /
    subtractions on reflectance values, such as in (b1_NIR+0.1) /
    b2_Red. See
    [here](https://docs.ropensci.org/MODIStsp/articles/discovered_bug.html)
    for further details. The bug is fixed as of version 1.3.7.

-   07/08/2018 - We are glad to report that `{MODIStsp}` is now included
    in the [rOpenSci](https://ropensci.org/about/) packages ecosystem.
    We thank reviewers Leah Wasser and Jeffrey Hanson for their valuable
    reviews, which helped us to further improve the package.

-   10/07/2018 - `{MODIStsp}` v. 1.3.6 is out (check out the [Release
    Notes](https://github.com/ropensci/MODIStsp/releases/tag/1.3.6) for
    further details).

-   20/06/2018 - `{MODIStsp}` v. 1.3.5 is out (check out the [Release
    Notes](https://github.com/ropensci/MODIStsp/releases/tag/v1.3.5) for
    further details).

-   11/04/2018 - Due to new NASA Policies the MODIS FTP servers were
    shut down starting, April 2, 2018. **FTP download is therefore no
    longer working** and will be removed in the next MODIStsp version.

-   11/04/2018 - [**Decommissioning of MODIS Version 5 Land Data
    Products**](https://lpdaac.usgs.gov/news/decommissioning-modis-version-51-land-cover-type-data-products-january-7-2019/).
    As per NASA notice above, MODIS v005 products are going to be
    decommissioned, and will soon be no longer available for download.
    Support for those products will be removed in the next MODIStsp
    version.

-   11/08/2017 - `{MODIStsp}` 1.3.3 was released today. It provides
    improvements in processing speed, as well as the usual bug fixes
    (thanks to all the users that signaled problems). Check the [Release
    Notes](https://github.com/ropensci/MODIStsp/releases/tag/v1.3.3) for
    further details.

-   25/07/2017 - As of today, **most of the content related to
    `{MODIStsp}` has been moved to our new website at
    [docs.ropensci.org/MODIStsp](https://docs.ropensci.org/MODIStsp/)
    **, which provides a much better user interface and ease of access
    to MODIStsp-related information. From now on, please **consult the
    new website for detailed and updated information on the package**.

-   Also our previous FAQ page on GitHub containing info for solving
    common installation, downloading and processing problems and issues
    was discontinued and **migrated at
    [docs.ropensci.org/MODIStsp/articles/faq.html](https://docs.ropensci.org/MODIStsp/articles/faq.html)**.

## Problems and Issues

-   Please **report any issues** you may encounter in our [issues page
    on github
    <i class="fa fa-github-square" aria-hidden="true"></i>](https://github.com/ropensci/MODIStsp/issues).

## <i class="fa fa-desktop" aria-hidden="true"></i> System Requirements

`{MODIStsp}` requires [`R`](https://cran.r-project.org) v \>= 3.6.3.

------------------------------------------------------------------------

# Installation Instructions

## <i class="fa fa-windows" aria-hidden="true"></i> Installing on Windows

You can install the stable version of `{MODIStsp}` from CRAN:

`install.packages("MODIStsp")`

, or the development version (containing the latest improvements and bug
fixes) from GitHub:

``` r
install.packages("remotes")
library(remotes)
install_github("ropensci/MODIStsp")
```

## <i class="fa fa-linux" aria-hidden="true"></i> Installing on Linux Systems

To install `{MODIStsp}` on Linux, you need to be able to install the
`{sf}` package, which requires several dependencies. See
[here](https://github.com/r-spatial/sf#installing) if you have trouble
installing `{sf}`.

Then, you can install the stable version of MODIStsp from CRAN:

``` r
install.packages("MODIStsp")
```

, or the development version (containing the latest improvements and bug
fixes) from GitHub:

``` r
library(devtools)
install_github("ropensci/MODIStsp")
```

## <i class="fa fa-apple" aria-hidden="true"></i> Installing on Mac

To install `{MODIStsp}` on Mac, you need to be able to install the
`{sf}` package, which requires several dependencies. See
[here](https://github.com/r-spatial/sf#installing) if you have trouble
installing `{sf}`.

Then, you can install the stable version of `{MODIStsp}` from CRAN:

``` r
install.packages("MODIStsp")
```

, or the development version (containing the latest improvements and bug
fixes) from GitHub:

``` r
library(devtools)
install_github("ropensci/MODIStsp")
```

# Usage

The easiest way to use `{MODIStsp}` is to use its powerful GUI
(Graphical User Interface) for selection of processing options, and then
run the processing.

To open the GUI, load the package and launch the MODIStsp function, with
no parameters:

``` r
library(MODIStsp)
MODIStsp()
```

This **opens a Shiny GUI** from which processing options can be
specified (and eventually saved or loaded). After specifying all
required parameters, clicking on “Start” will start the processing (see
[here](https://docs.ropensci.org/MODIStsp/articles/interactive_execution.html)
for more detailed instructions).

`{MODIStsp}` can also be launched in non-interactive mode within an `R`
session or script by setting the optional `GUI` parameter to FALSE, and
the `opts_file` parameter to the path of a previously saved JSON Options
file. This allows to exploit `{MODIStsp}` functionalities within generic
“R” processing scripts.

``` r
library(MODIStsp) 
# --> Specify the path to a valid options file saved in advance from MODIStsp GUI 
opts_file <- "X:/yourpath/youroptions.json" 
  
# --> Launch the processing
MODIStsp(gui = FALSE, opts_file = opts_file)
```

Finally, `{MODIStsp}` can be run by manually specifying all processing
arguments, or by overwriting some of the arguments contained in a saved
json file in the call to the package, such as in:

``` r
library(MODIStsp) 
# --> Specify the path to a valid options file saved in advance from MODIStsp GUI 
opts_file <- "X:/yourpath/youroptions.json" 
  
# --> Launch the processing
MODIStsp(gui        = FALSE, 
         opts_file  = opts_file, 
         start_date = "2020.05.01", 
         end_date   = "2020.08.01", 
         spatmeth   = "file", 
         spafile    = "X:/path_to/spatial_extent_file.gpkg")
```

, where we are overwriting the options related to spatial and temporal
extent contained in the options file with new values. This allows easily
running processing based on the same main options (e.g., product layers,
output format, etc.) but changing on the fly the desired ones.

See
[here](https://docs.ropensci.org/MODIStsp/articles/noninteractive_execution.html)
for more detailed instructions and examples.

# Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/MODIStsp/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.


# MODIStsp 2.0.6

## Minor changes
- Replace `M*D17A3H` with `M*D17A3HGF` and add `M*D17A2HGF` (#237)
- Avoid errors in case of missing internet connection

## Bug fixes 
- Patch for bbox loaded from json in case of drawn extent (#228)
- Fix #226
- Fix #232 
- Fix #234
- Fix CRAN notes


# MODIStsp 2.0.5

## Main changes

- Edit documentation related to the change of maintainer
    (see https://docs.ropensci.org/MODIStsp/articles/lorenzo).
    
- Add the argument `parallel` to function `MODIStsp()` and `MODSIStsp_process()`
    to allow running the processing in single core modality.
    
## Minor changes

- Fix Travis tests

- Bug fix (#222)


# MODIStsp 2.0.3

## Main changes

- This submission should fix errors on Debian CRAN builds, due to improper 
trigger of an internal function leading to writing in the user's 
lib folder. 

- Fixes a bug leading to crash when using scale_val = TRUE and change_no_data = FALSE

- Fixes a bug leading to the GUI crashing rather than giving info messages in 
  case not all input parameters are specified

- Implements redirection to MODIS products web pages when pressing 
 the corresponding button
 
- Modifies slightly the Shiny GUI


# MODIStsp 2.0.0

## Main changes

- Replace the old gWidgets-based GUI with a new one based on Shiny;

- Enhances support for CLI usage. Now all parameters can be passed to 
the `MODIStsp` function. If also a opts_file is passed, values specified
explicitly in the call override those in the options file;

- Fixes problems in retrieval of corners for MODIS products in 4008 projection (fixes #204);

- Fixes problems/improves support for datasets with multiple NoData values.
Now, all NoData values are kept to original values if NoData change is set
to FALSE. Also, Scale/Offset are no longer wrongly applied also to NoData values when scaleval = TRUE;


# MODIStsp 1.4.0

## Main changes

- switch to use of GDAL3/PROJ6 WKTs for projection representation, using sf::gdal_utils to perform gdalwarp/gdaltranslate instead of gdalUtils on external GDAL.

- switch to sf for all internal work on vector data.

- Remove sp, rgdal, rgeos, pacman, gdalUtils dependencies

- Adds support for products MCD19A1 and MCD19A2 products

# MODIStsp 1.3.9

## Main changes

- Fixes a bug causing crashes on MOD14A1 product

- Adds support for product MCD12Q2 and removes support for no longer available version 5 of some products.

- Updates web links to MODIS products description pages

# MODIStsp 1.3.8

## Main changes

- Fixed an issue causing incorrect application of scale/offset values on GDAL 
versions > 2.3 (introduced by change of behaviour of gdal_buildvrt - https://trac.osgeo.org/gdal/ticket/3221#comment:5) - see https://github.com/ropensci/MODIStsp/issues/163

- Added support for the following products: MOD21A1D.006 MOD21A1N.006 MOD21A2.006/

## Bug fixing

- Fixed an issue preventing correct processing of some products in offline mode (
https://github.com/ropensci/MODIStsp/issues/142)


# MODIStsp 1.3.7

## Main changes

- Fixed a bug leading to incorrect computation of custom spectral indices containing "additive"
parameters (e.g., (NIR+0.1)/(Red-0.2)) when scale_val == FALSE

## Bug fixing

- Fixed a bug leading to not properly reassigning NoData values when processing 
  full tiles and using change_nodata = TRUE

- Fixed inconsistencies in definition of characteristics of products
  MOD/MYD13C2/C1 and MOD/MYD13A3 (erroneous layers in xml)

- Fixed a bug leading to help messages in the select layers portion of the GUI 
to not render

- Updated MOD44B specifications to allow download of v006 version

# MODIStsp 1.3.6

## Main changes

Maintenance release to solve CRAN build errors on debian, due to the test_resetindexes
test. The test is now skipped on CRAN. Additionally, the MODIStsp_addindex 
function was modified to require explicit permission by the user to write on 
the MODIStsp_previous.json file

## Bug fixing

- Fixed a bug leading to errors if only "Aqua" platform was selected for download [#133](https://github.com/ropensci/MODIStsp/issues/133)

# MODIStsp 1.3.5

## Main changes

Maintenance release to solve CRAN build errors on debian, due to the test_addindex
test. The test is now skipped on CRAN. Additionally, the MODIStsp_addindex 
function was modified to require explicit permission by the user to write on 
the MODIStsp_previous.json file

## Bug fixing

- Fixed bug leading to errors in processing extent when switching products with different Native projection (4008 vs sinusoidal), the projection string was not properly updated. [77f5693e9](https://github.com/ropensci/MODIStsp/commit/77f5693e9e1e180f05efaa04fa031567e782ba89)

- Fixed warnings on check for uniqueness in http addresses

# MODIStsp 1.3.4

## Main changes

### Breaking changes

- Due to improvements and changes in the GUI (see below), `MODIStsp` .json options 
  files saved with older versions are no longer supported. Users will be informed of 
  this if trying to use an obsolete json file.
  
- **Removed support for FTP download** due to switch-off of NASA servers.

### Updates in supported products

- **Removed all v005 and earlier products**, due to discontinuation of their 
distribution by NASA servers

- **Added support** for the following products:

  MCD64A1; MCD12C1; MCD18A1; MCD18A2; MCD12Q1; MOD44B; MOD44W; MCD12C1; MCD12Q1;
  MOD12A2; MOD12A3

### Improvements in download functions

- **Improvements in GUI**. It is now possible to set the processing extent interactively
  using the "Select on Map" button. This opens a browser window allowing to select
  and extent. 

- Use of `httr::RETRY` to improve behavior while navigating the servers to
  retrieve available files and while downloading hdf file (when use_aria == FALSE), 
  thus removing dependency to RCurl. 
  
### Improvements in processing functions

 - **Improved functionality for dealing with NoData** for products with multiple 
   fill-values. If "Change NoData" is set to "Yes", then in case a layer 
   has multiple Nodata values all those values are set to NA in the output 
   (see github.com/ropensci/MODIStsp#113)

### Extensive code refactoring for submission to ropensci. 

- Long functions were simplified/split into smaller functions to allow for 
  easier maintenance
- GUI event handlers were moved into dedicated "R" files
- Extensive code linting to abide to ropensci standards
- Switch to jsonlite/xml2 to abide to ropensci standards
- Removal of some less-used dependencies (e.g., hash)

### Improvements in documentation and website

- More detailed documentation for the main functions
- Improvements in pkgdown articles to include some "worked examples" (e.g., 
  MODIStsp_Extract)
- New article in pkgdown showing a list of supported products

### Improvements in test coverage

- Several new tests added, bringing coverage above 90%

### New functions

- Added `MODIStsp_resetindexes` to remove all custom indexes from a MODIStsp 
  json options file and `MODIStsp_reset_options` to reset MODIStsp options to 
  default.
  
### Bug fixing

- Fixed bug affecting extent selection when working with non-tiled (MCD) products
https://github.com/ropensci/MODIStsp/issues/122

- Fixed bugs affecting the "Seasonal" time series download 
  
________________________________________________________________________________

  
# MODIStsp 1.3.3.1 Release Notes

This was mostly a maintenance release - fixing a bug in 1.3.3 submission related
to a missing import in NAMESPACE

## Minor Changes

Improved organization of Virtual Raster files and RData files in the "Time_Series"
output subfolder. Now virtual files and RData files are organized by sensor and 
layer to facilitate access. 

____________________________________________________________________________________

# MODIStsp 1.3.3 Release Notes

v1.3.3 was released on 10/08/2017

## Major Changes

-  Improved speed in computation of spectral indexes, quality indicators and in 
   computation of scaled variables by using `raster::calc()` and `raster::overlay`
   (commits [0f5d76d](https://github.com/ropensci/MODIStsp/commit/0f5d76de1661958cd5cbaa79f8115035cb9c348e),     [0f5d76d](https://github.com/ropensci/MODIStsp/commit/0f5d76de1661958cd5cbaa79f8115035cb9c348e), [e462721](https://github.com/ropensci/MODIStsp/commit/e462721a06a079185ec5a84270ea0c8bd8edf54d))
   
-  Added functionality for unit testing using `testthat` and codecov integration.
   (commit [0c00fc6](https://github.com/ropensci/MODIStsp/commit/0c00fc6bf07aed046b2b198e0278ab3264e5298a)
   and others)
   
-  Added "testing mode" to allow users to test proper functioning. Now, running 
   `MODIStsp(test = X)` (with X in (0,6)) runs the processing using default processing
   parameters  (commit [0c00fc6](https://github.com/ropensci/MODIStsp/commit/0c00fc6bf07aed046b2b198e0278ab3264e5298a) and others)

## Minor Changes

-  Suppression of verbose messages and (useless) warning messages while parsing the NASA
servers and downloading data using "ftp" ( [3775d60](https://github.com/ropensci/MODIStsp/commit/3775d6099bc359925d3dcbd96c2ffe8455502648));

## Bug fixing

-   Fixed a bug preventing the "last" choice (or that present in the json file) from 
    correctly showing in the GUI upon launch/restore of a saved json file (commit
    [633c2dd](https://github.com/ropensci/MODIStsp/commit/633c2dddd29d45c618e4ca121112000ceefe91e3))

-   Fixed a bug affecting MODIS layers coded as Unsigned Integer - 32 bit (Thanks
    to Rob Critchlow for signaling this). The bug was due to improper handling of 
    UInt32 data in `gdalbuildvrt`, causing sometimes an incorrect translation from 
    HDF to output formats ([#72](https://github.com/ropensci/MODIStsp/issues/72)).

     **M\*D09A1** - 500m Reflectance Band Quality (V005 and V006); 
     **M\*DO9CMG** - Coarse Resolution QA (V005 and V006);
     **M\*D09CMG** - Coarse Resolution Number Mapping (V006); 
     **M\*D09GA** - 500m Reflectance Band Quality (V005 and V006);
     **M\*DOCGA** - Band quality for MODIS bands 8-15 (V006);
     **M\*D11C3** - Days with clear-sky conditions and validated LSTs;
                    Nights with clear-sky conditions and validated LSTs (V005 and V006); 
     **MCD43A2** - BRDF\_Albedo\_Band\_Quality (V005 and V006).

- Fixed a bug affecting creation of time series files (RData and virtual rasters) 
  on all MCD products ([#77](https://github.com/ropensci/MODIStsp/issues/77))

- Fixed a bug a error on creation of "burn_date" layers for MCD45A1 product 
  ([#77](https://github.com/ropensci/MODIStsp/issues/77))

- Fixed bugs on specifying spatial extent files on non-interactive execution
  ([#75](https://github.com/ropensci/MODIStsp/issues/75))

____________________________________________________________________________________

## 17/04/2017 - MODIStsp is now on CRAN !

MODIStsp was recently accepted on CRAN. From now on, you can install it simply using

`install.packages("MODIStsp")`

You'll however still be able to install the `development` version from github,
containing the last improvements and bug fixing using:

`install_github("ropensci/MODIStsp", ref = "master")`

____________________________________________________________________________________

# MODIStsp 1.3.2 Release Notes

v1.3.2 was released on 22/03/2017

## Major Changes:

- Added functionality to apply scale and offset coefficients on MODIS original values
  according with the specifications of single MODIS products.

### Details:

- MODIS hdf datasets are always stored as integer values, with scales factor and/or
  offsets to apply in order to convert to the indicated measure units reported in
  the products' documentation.

- Starting from v1.3.2: 

    - Leaving the "Scale output values" option to "No", output files are left as
      provided by NASA, and additional indices are produced as integer values with
      a 10000 factor scale; 
    - Setting the "Scale output values" option to "Yes", scale factor and offsets 
      are applied if existing (for example, in this case Land Surface Temperature
      values in the output raster will be in °K), and  spectral indices are floating
      point values (for example, NDVI is between -1 and 1, etc.).

## Minor Changes:

- Some product names and output layer names were modified to reduce the length of
  output file names, homogenize the names of the outputs and correct some errors.
__For compatibility with already created output files__ (versions up to 1.3.1),
  the old "XML" file specifying output files format is still available in
  `inst/ExtData/MODIStsp_ProdOpts_old_v1.3.1.xml`. To use the old file naming
  conventions, users have to:

    1. delete `inst/ExtData/MODIStsp_ProdOpts.xml` and rename `MODIStsp_ProdOpts_old_v1.3.1.xml`
       to `MODIStsp_ProdOpts.xml`.
    2. delete `MODIStsp_ProdOpts.RData` from the `Previous` folder within 
       `your_R-library_path/MODIStsp/Previous`
    3. Restart `MODIStsp`  
    
    <br>
    
- Timeouts of httr/ftp requests were increased to prevent problems on download on
  slow connections

## Bug fixing:

- Fixed bug on FTP download speed (Issue [#65](https://github.com/ropensci/MODIStsp/issues/65)) 
- Fixed bug on download of tile 0, preventing download of images with DOY 001 and
  of all "yearly based" products (e.g., MOD17)(Issue [#64](https://github.com/ropensci/MODIStsp/issues/64)) 
- Fixed other bugs affecting FTP download (https://github.com/ropensci/MODIStsp/commit/efbf1b469e7518ffc8a7ec6d9922242d6a5c228f, https://github.com/ropensci/MODIStsp/commit/1dc53a5ff5b355965acec86678a3104bd2d27fd9, https://github.com/ropensci/MODIStsp/commit/fa6c7b42eadce516a2f781604c9db28418120f36)

____________________________________________________________________________________

# MODIStsp 1.3.1 Release Notes

v1.3.1 was released on 13/02/2017

## Major Changes

- Added functionality for processing of Snow Cover datasets: MOD10A1, MOD10A2, 
  MOD10C1, MOD10C2, MOD10CM (Issue
[#55](https://github.com/ropensci/MODIStsp/issues/55)) on devel

- Added functionality for downloading "partial" years (Issue
  [#54](https://github.com/ropensci/MODIStsp/issues/54)) on devel

- Added functionality for computing vegetation indexes on MCD43A4 (v5-v6),
  MCD43B4 (v5), MCD43C4 (v5-v6) (Issue [#59](https://github.com/ropensci/MODIStsp/issues/59))
  on master/devel

- Added functionality for accelerating download using aria2c (Issue 
  [#55](https://github.com/ropensci/MODIStsp/issues/55)) on devel

## Bug fixing

- Fixed bug on download with aria, throwing an error on partial download on http
  download with aria ([6fbc875](https://github.com/ropensci/MODIStsp/commit/6fbc87547b6214b500afc0291c02166c0b855c78))

- Fixed bug on M*D15A2 processing (Issue [#60](https://github.com/ropensci/MODIStsp/issues/60))
  on devel/master

- Fixed bug on MCD12Q1 processing (Issue [#58](https://github.com/ropensci/MODIStsp/issues/58))
  on devel/master

- Fixed bug on MOD13C2 processing (Issue [#52](https://github.com/ropensci/MODIStsp/issues/52))
  on devel/master

- Fixed bug on insertion of custom projection (Issue [#57](https://github.com/ropensci/MODIStsp/issues/57))
  on devel/master

- Fixed bug on selection of custom index (Issue [#53](https://github.com/ropensci/MODIStsp/issues/53))
  on devel/master

____________________________________________________________________________________

# MODIStsp 1.3.0  Release Notes

v1.3.0 was released on 11/05/2016

## Major Changes

- Added functionality for downloading and preprocessing MODIS collection 006 datasets.
  For products with both 005 and 006 collections, the user can select the version
  using a new droplist in the GUI.

- Added functionality for off-line processing. This allows both to _i)_ reprocessing
  already downloaded data (for example, to create time series for an additional
  layer) without the need to connect to NASA servers, and _ii)_ process HDF files 
  downloaded outside _MODIStsp_ (e.g., directly from NASA ftp) and stored on the
  user's PC, without the need of an active internet connection. 

- Improved the way in which options are saved. Much more readable. JSON files are
  now used instead than .RData. User options are no longer saved alongside products
  characteristics. This will allow to re-use an "old" options file even if changes 
  are made on the XML file describing the products.

- Improved the GUI interface for specifying additional Spectral Indexes. Hints are 
  now showed to the user, and multiple indexes can be added in the same session.

## Minor Changes

- General improvements in the GUI interface. Products are now grouped by categories,
  to allow easier identification and selection.

- Improvements in the README file and vignettes, providing more instructions on package
  use.

- Improved functionality for checking for "complete" download, by comparing the 
  size of the downloaded files with that of files on the server. 

- Added "configure" file for Linux installation.

- Temporary files necessary for processing (e.g., vrt files) are now created (and
  destroyed) within the "R" temporary folder.

- Miscellaneous bug-fixing


____________________________________________________________________________________

# MODIStsp 1.2.1 Release Notes

v1.2.1 was released on 11/05/2016
 
## Major Changes

1. Modified format of "R" output time series from _rts_ objects to _RasterStack_ 
   objects with temporal information added in the "z" attribute via setZ()

2. Major changes/improvements in _MODIStsp\_extract_ function:
    * Use of plain `RasterStack` with "z" attribute instead than `rasterstackts`
    * Use of gdal\_rasterize (_gdalUtils_) instead of rasterize (_rgdal_) to improve 
      speed. Temporary shapes and rasters necessary are saved in "R" temporary folder
      and removed automatically
    * Fixed bugs on functionality for point/lines shapefiles, according to what 
      specified by the "small" and "small_method" parameters
    * Added functionality for retrieving data for small polygons
    * Added out_format selection - xts or plain data.frame
    * Added possibility to use a shp filename as input directly
    * Added conformity checks on inputs
    * Added functionality to run without specifying start and end dates
    * Added id_field parameter for choosing which column of the input SP object
      should be used for "naming" the columns of the output
  
3. Removed possibility to use "complex" resampling methods when reprojecting (e.g.,
   bilinear, cubic, etc.) to avoid incorrect resampling on categorical variables
   and "contamination" of good pixels data. 

## Minor Changes

* Changed the input method for starting and ending dates selection in the GUI.
  Now a text field is used 
* Added functionality for writing data ignore value on ENVI files
* Removed automatic deletion of XML files created by writeRaster to keep metadata
  information
* Changed names of products in the GUI for products with both TERRA and AQUA
  dataset to M\*D09A1, M\*D13Q1, etc...
* Modified code syntax to satisfy R code styling guidelines
* Modified roxygen parameters so that only required functions are imported from 
  dependent packages
* Updated and corrected the list of dependencies
* Updated required "R" version to 3.2, and minimum versions for dependent packages 
  to current versions.
* Added Welcome message
* Updated links to LPDAAC product description pages
* Changed all "print" and "cat" calls to show messages/warnings to "message" or
  "warning" to allow easy disabling MODIStsp verbose messages
* Using "R" tempfile/tempdir to save vrt files

## Bug Fixes

* Corrected a bug that threw an error in case incorrect bounding box specified

________________________________________________________________________________


# MODIStsp 1.2.0 Release Notes

v1.2.0 was released on 29/07/2015

## Major changes

First stable release of advanced implementation of MODIStsp !
We know it should be 1.0.0, but that's it !
MODIStsp 2.0.6
================

* Windows 10 on local install, R 4.1.0
* ArchLinux on local install, R 4.1.0
* win-builder (R-devel, R-release, R-oldrelease)

## R CMD check results

There were no ERRORs, WARNINGs nor NOTEs.


MODIStsp 2.0.5
================

* Windows 10 on local install, R 4.0.3
* Ubuntu 18.04 on local install, R 4.0.3
* ArchLinux on local install, R 4.0.3
* win-builder (R-devel, R-release, R-oldrelease)

This submission should fix errors on Debian builds, due to improper 
trigger of an helper function leading to writing in the user's 
lib folder. 

Also fixes a couple of bugs.

## R CMD check results

There were no ERRORs nor WARNINGs.

NOTE:
```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Luigi Ranghetti <luigi@ranghetti.info>'

New submission

Package was archived on CRAN
```
Unfortunately Lorenzo Busetto, who maintained the package, suddently passed away
(https://docs.ropensci.org/MODIStsp/articles/lorenzo).
After that, package was automatically archived for policy violation, since 
communications sent to him were not read.
I now started maintaining MODIStsp (I already authored parts of code).

```
CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2020-10-31 for policy violation.

  Attempts writing to the user library.
```
In this submission I fixed a possible tentative of writing in the user library
(commit ebdfa2a6ebbadf1ab42bc4679b48d2224ddf333f).

```
Possibly mis-spelled words in DESCRIPTION:
  Busetto (40:5)
  HDF (29:24)
  MODIS (2:35, 25:10, 29:18, 30:60, 38:18)
  Ranghetti (40:17)
  mosaicking (26:43)
  rasters (24:63)
  reflectance (31:50, 32:65)
  reprojecting (26:55)
  resizing (26:72)
```
All these words are correctly spelled.

## CRAN resubmission review

Please find below the answers to the CRAN reviewer.

> Please reduce the length of the title to less than 65 characters.

The title was changed from "A Tool for Automating Download and Preprocessing of 
MODIS Land Products Data" to "Find, Download and Process MODIS Land Products".


> Please only capitalize sentence beginnings and names in the description text.
> e.g.  Quality Indicators --> quality indicators
>          Spectral Indexes --> spectral indexes
>          Surface Reflectance --> surface reflectance

All capitalized letters were changed (unless sentence beginnings).


> Please change:
> products ,
> -->
> products,

Done.


> Please provide a link to the used webservices to the description field
> of your DESCRIPTION file in the form
> <http:...> or <https:...>
> with angle brackets for auto-linking and no space after 'http:' and
> 'https:'.

No webservices were used (the MODIStsp Shiny GUI is launched from local PC).
I added a reference to the paper which describes MODIStsp:
Busetto and Ranghetti (2016) <doi:10.1016/j.cageo.2016.08.020>


> You have examples for unexported functions.
> Please either omit these examples or export these functions.
> Used ::: in documentation:
>       man/split_nodata_values.Rd:
>          MODIStsp:::split_nodata_values(c("255", "250,254:255"))
>       man/split_nodata_values.Rd:
>          MODIStsp:::split_nodata_values(c("255", "250,254:255"),
> take_all = FALSE)
>       man/split_nodata_values.Rd:
>          MODIStsp:::create_nodata_rcl(c("255", "250,254:255"), c("255",
> "255"))

Examples of unexported functions were removed.


> Some code lines in examples are commented out in MODIStsp.Rd.
> 
> Please never do that. Ideally find toy examples that can be regularly
> executed and checked. Lengthy examples (> 5 sec), can be wrapped in
> \donttest.

Examples were changed following these indications.
All lenghy examples were wrapped in \donttest{}.
\dontrun{} is no more used.


> Please ensure that you do not use more than 2 cores in your examples,
> vignettes, etc.

Argument `parallel` was added to functions `MODIStsp()`, MODIStsp_process()`
and `MODIStsp_process_bands()` (the only functions exploiting multicore 
computation): if `parallel = FALSE`, single core is used.
All examples and vignettes were edited setting `parallel = FALSE` everywhere.


> Please always add all authors, contributors and copyright holders in the
> Authors@R field with the appropriate roles.
> e.g.: Babak Naimi

Babak Naimi's contribution was added to the DESCRIPTION.


MODIStsp 2.0.5
================

* Windows 10 on local install, R 4.0
* Ubuntu 18.04 on local install, R 3.6.3
* win-builder (R-devel, R-release)

## R CMD check results

There were no ERRORs, WARNINGs and NOTES

MODIStsp 2.0.4
================

* Windows 10 on local install, R 4.0
* Ubuntu 18.04 on local install, R 3.6.3
* win-builder (R-devel, R-release)

## R CMD check results

There were no ERRORs, WARNINGs

There is one note related to nuw submission (package was archived before I had time to fix it)

This submission should fix errors due to trying accessing temporarily non available
urls when building vignettes, which led to CRAN archival. Sorry for that. 

Also limits number of retries on HTTR GET to 5 (previously 20) as suggested by
DR Riply in his mail

MODIStsp 2.0.3
================

* Windows 10 on local install, R 4.0
* Ubuntu 18.04 on local install, R 3.6.3
* win-builder (R-devel, R-release)

## R CMD check results

There were no ERRORs, WARNINGs and NOTES

This submission should fix errors on Debian builds, due to improper 
trigger of an helper function leading to writing in the user's 
lib folder. 

Also fixes a couple of bugs.


MODIStsp 2.0.1
================
* Windows 10 on local install, R 3.6.3
* Windows 10 on local install, R 4.0
* Ubuntu 18.04 on local install, R 4.0
* win-builder (R-devel, R-release)

## R CMD check results

There were no ERRORs, WARNINGs and NOTES

This is a resubmission, fixing notes related to redirected URLS. Sorry for that.

MODIStsp 2.0.0
================
* Windows 10 on local install, R 3.6.3
* Windows 10 on local install, R 4.0
* Ubuntu 18.04 on local install, R 4.0
* win-builder (R-devel, R-release)

## R CMD check results

There were no ERRORs, WARNINGs

There is a NOTE related to inability to retreive future timestamps.
There is a NOTE related to possible invalid urls on some builds, but the url is valid.

MODIStsp 1.4.0
================

## Test environments

* Windows 10 on local install, R 3.6.3
* Windows 10 on local install, R 4.0
* Ubuntu 16.04 on travis-ci (devel and release)
* win-builder (R-devel, R-release)

## R CMD check results

There were no ERRORs, WARNINGs

There is a NOTE related to suggesting a orphaned package (gWidgets)
I received a mail from Prof. Ripley about that recently. As suggested, I temporarily moved gWidgets to suggests and use it conditionally for the time being, with the intention of removing the dependency in the next MODIstsp release by switching to a Shiny-based GUI. 

## Note on installation on macos 

Build on macos currently fails (also for MODIStsp 1.3.8) with error

"ERROR Package required but not available: ‘gWidgetsRGtk2’"

due to erroring in build of package "RGtk2". 

MODIStsp can be however installed by installing RGtk2 beforehand, following
instructions reported here: 

http://lbusett.github.io/MODIStsp/articles/installation.html#installing-on-mac


MODIStsp 1.3.9
================

## Test environments

* ubuntu 18.04 on local install, R 3.6.0
* ubuntu 14.04 on travis-ci (devel and release)
* win-builder (R-devel, R-release)
* windows 10 on local install, R 3.6.0

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs

## Note on installation on macos 

Build on macos currently fails (also for MODIStsp 1.3.8) with error

"ERROR Package required but not available: ‘gWidgetsRGtk2’"

due to erroring in build of package "RGtk2". 

MODIStsp can be however installed by installing RGtk2 beforehand, following
instructions reported here: 

http://lbusett.github.io/MODIStsp/articles/installation.html#installing-on-mac



## Test environments

* ubuntu 18.04 on local install, R 3.5.2
* ubuntu 14.04 on travis-ci
* win-builder (R-devel, R-release)
* windows 10 on local install, R 3.5.2

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs

## Note on installation on macos 

Build on macos currently fails (also for MODIStsp 1.3.3.1) with error

"ERROR Package required but not available: ‘gWidgetsRGtk2’"

due to erroring in build of package "RGtk2". 

MODIStsp can be however installed by installing RGtk2 beforehand, following
instructions reported here: 

http://lbusett.github.io/MODIStsp/articles/installation.html#installing-on-mac


MODIStsp 1.3.7
================

## Test environments

* ubuntu 18.04 on local install, R 3.5.0
* ubuntu 14.04 on travis-ci (https://travis-ci.org/ropensci/MODIStsp/builds/462820855)
* win-builder (R-devel, R-release)
* windows 10 on local install, R 3.5.1

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs

## Note on installation on macos 

Build on macos currently fails (also for MODIStsp 1.3.3.1) with error

"ERROR Package required but not available: ‘gWidgetsRGtk2’"

due to erroring in build of package "RGtk2". 

MODIStsp can be however installed by installing RGtk2 beforehand, following
instructions reported here: 

http://lbusett.github.io/MODIStsp/articles/installation.html#installing-on-mac


MODIStsp 1.3.6
================

Maintenance release to solve CRAN build errors on debian, due to the test_resetindexes
test (I'm sorry I missed this in the previous resubmission). 
The test is now skipped on CRAN. Additionally, the MODIStsp_resetindexes
function was modified to require explicit permission by the user to write on 
the MODIStsp_previous.json file

## Test environments

* ubuntu 18.04 on local install, R 3.4.4
* ubuntu 14.04 on travis-ci (https://travis-ci.org/lbusett/MODIStsp/jobs/402298232)
  (R-devel, R-release)
* win-builder (R-devel, R-release) (r-release and devel)
* windows 10 on local install, R 3.5

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs

MODIStsp 1.3.5
================

Maintenance release to solve CRAN build errors on debian, due to the test_addindex
test. The test is now skipped on CRAN. Additionally, the MODIStsp_addindex 
function was modified to require explicit permission by the user to write on 
the MODIStsp_previous.json file

## Test environments

* ubuntu 18.04 on local install, R 3.4.4
* ubuntu 14.04 on travis-ci (https://travis-ci.org/lbusett/MODIStsp/builds/394612258)
  (R-devel, R-release)
* win-builder (R-devel, R-release)
* windows 10 on local install, R 3.5

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs


MODIStsp 1.3.4
================

## Resubmission 
 
This is a resubmission. In this version I have: 
 
* Re-added package "testthat" in Suggests and removed an unnecessary call to  
library(testthat) to avoid the following WARNING on Debian: 
 
  "* checking for unstated dependencies in ‘tests’ ... WARNING 
  '::' or ':::' import not declared from: ‘testthat’ 
  'library' or 'require' call not declared from: ‘testthat’" 
  
## Test environments

* ubuntu 18.04 on local install, R 3.4.4
* ubuntu 14.04 on travis-ci (https://travis-ci.org/lbusett/MODIStsp/jobs/383285994#L6)
* win-builder (R-devel, R-release)
* windows 10 on local install, R 3.5

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs

## Note on installation on macos 

Build on macos currently fails (also for MODIStsp 1.3.3.1) with error

"ERROR Package required but not available: ‘gWidgetsRGtk2’"

due to erroring in build of package "RGtk2". 

Note that MODIStsp can be however installed by installing RGtk2 beforehand, following
instructions reported here: 

http://lbusett.github.io/MODIStsp/articles/installation.html#installing-on-mac


MODIStsp 1.3.3.1
================

re-submission of v 1.3.3 patching an issue related to missing import of 
gWidgetsRGtk2, leading to NOTES on Solaris, Fedora and OSX builds (see this 
thread on r-package-devel for details: 
https://stat.ethz.ch/pipermail/r-package-devel/2017q3/001790.html)

There were no ERRORs or WARNING. There was 1 NOTE:

   * checking CRAN incoming feasibility ... NOTE
   Maintainer: 'Lorenzo Busetto <lbusett@gmail.com>'

   Days since last update: 4

due to close resubmission. Sorry for this.

For comments on build environments etc. see the notes below concerning v 1.3.3. 

MODIStsp 1.3.3
==============

## Test environments
* ubuntu 17.10 on local install, R 3.4.1
* ubuntu 16.04 on travis-ci (https://travis-ci.org/lbusett/MODIStsp/builds/263768847)
* win-builder (R-devel)
* windows 10 on local install, R 3.4.1 (As for the previous reelease, 
  R CMD check passes if GTK+ library is properly installed and on Windows PATH,
  otherwise the check causes an endless GTK+ installation loop. This seems a 
  common behaviour for packages relying on RGtk2/gWidgetsRGtk2)
* MacOs on rhub:
  
  Build on Maverick works (https://builder.r-hub.io/status/MODIStsp_1.3.3.tar.gz-3ec4e4680884dc9ca4ebca061af09da3)
  
  Build on El Capitain or Sierra is not possible at the moment, due to the lack
  of updated binaries for RGtk2 (https://builder.r-hub.io/status/MODIStsp_1.3.3.tar.gz-c0e71d5c8f0a42ef85157b346255b4bb)
  A workaround solution using oldrel (Maverick) binaries is provided in MODIStsp README and website


## R CMD check results

There were no ERRORs, WARNINGs or NOTEs


## Downstream dependencies

This package has no downstream dependencies.

--------------------------------------------------------------------------------
MODIStsp 1.3.2
==============

## Test environments
* ubuntu 16.10 on local install, R 3.3.3
* ubuntu 14.04 on travis-ci, R 3.3.3 (https://travis-ci.org/lbusett/MODIStsp/builds/222347443)
* win-builder (R-devel)
* windows 10 on local install, R 3.3.3 (R CMD check passes if GTK+ library is
  properly installed and on Windows PATH, otherwise the check causes an endless
  GTK+ installation loop. This seems a common behaviour for packages relying on 
  gWidgetsRGtk2)
* local OS X install, R 3.3.3

## R CMD check results

There were no ERRORs, WARNINGs 

There was 1 NOTE, related to the fact that this is a first submission.

There was a warning in win_builder about the following (possibly) invalid URL:

https://notehub.org/fctdn

I checked it, and it's working.

## Downstream dependencies

This package has no downstream dependencies.
---
name: Help needed / question
about: Use this in case of problems or doubts using {MODIStsp}
title: ''
labels: assistance
assignees: ''

---

<!--

IMPORTANT NOTE: the following template must not be deleted, but read and filled with the required information. Issues written without following the template will be marked as "malformed" and ignored.

Use this template if you need assistance running {MODIStsp} on your code (e.g., in case of errors which are not a bug / you are not sure if they are a bug), or if in doubt about which template you should use. Please use this method instead than sending private email to the authors.

Before opening a new issue please check if the problem was already been mentioned (if so but the found issue is closed, open a new issue citing the old task instead than reopening it).

Ensure that your {MODIStsp} version is update with the last CRAN version:
install.packages("MODIStsp")

Please take particular care with code reproducibility (follow indications provided in the template).

Due to the limited time available to address all the issues, the developer will preferentially address requests finalised to publish scientific works; in this case, please provide the required information in the template below.

IMPORTANT NOTES ABOUT NETIQUETTE
1. Please remember that {MODIStsp} is not a commercial tool, so the developer is not obliged to provide assistance: please be polite, be patient if developers will not answer you instantly and respect the Code of Conduct (https://ropensci.org/code-of-conduct/).
2. Your are required to answer when details (generally outputs of R commands) are required, and to provide a feedback after opening an issue, even after solving your problem or if you are not yet interested in solving it. In the case of missing feedback, the developer reserve the right to ignore your future requests.
3. Tasks can be closed after 10 days of inactivity (you can reopen it if you need further help).
-->

**Issue description**
<!-- Add here a clear and concise description of what the problem is about. -->

**Reproducible example**
<!-- Please provide here a reproducible example in the chunk below. -->

```r
## PLEASE DELETE AND WRITE YOUR OWN
library(MODIStsp)
MODIStsp(gui = FALSE, opts_file = "/path/of/the/parameter_file.json")
# file parameter_file.json _must_ be attached or copied as text,
# as well as referenced files (e.g. spafile)
```

**Expected and actual behavior**
<!-- Provide here the full output of the provided example and describe what is going wrong- -->

```
## PASTE HERE THE OUTPUT OF YOUR EXAMPLE CODE
```

**System information**
<!-- Provide here the output of the following R commands:
sessionInfo()
packageVersion("MODIStsp")
 -->

```
# PASTE HERE YOUR OUTPUT
```

**Additional context**
<!-- Add here any other context about the problem here (for example, the content of the output folder in case the error appears during a subsequent code execution. -->

**Scientific publication**
<!-- If your work is finalised to publish a scientific paper, please provide here the following details:
1. which is the aim of your work;
2. available details - even if subjected to be modified - about publication (title, authors, candidate journal).
Please remember to cite {MODIStsp} in your work (see https://docs.ropensci.org/MODIStsp/authors.html ). -->
---
name: Missing / unclear documentation
about: Question about general {MODIStsp} functionalities
title: ''
labels: question
assignees: ''

---

<!--
Use this template if you need information about {sen2r} funcionalities which are not [well] documented. You can use it also for GENERAL questions (e.g. "Is it possible to [...] with sen2r?"); instead, for questions related to use cases (or if in doubt about which template to use) please use the "Help needed" template.

Before opening a new issue:
1. please read the online documentation at https://sen2r.ranghetti.info/ (https://sen2r.ranghetti.info/reference/ in case of a question about a specific function);
2. if your question is related to the {sen2r} GUI, read the embedded documentation ("?" marks in the GUI);
3. check if the question was already been mentioned as a GitHub issue.

If your question is not general but related with your specific use case, please use the "Help needed" template.

IMPORTANT NOTES
1. Please remember that {sen2r} is not a commercial tool, so the developer is not obliged to provide assistance: please be polite, be patient if noone will answer you instantly and respect the Code of Conduct (https://sen2r.ranghetti.info/CODE-OF-CONDUCT.html)
2. Your are required to answer when details (generally outputs of R commands) are required, and to provide a feedback after opening an issue, even after solving your problem or if you are not yet interested in solving it. In the case of missing feedback, the developer reserve the right to ignore your future requests.
3. Tasks can be closed after 10 days of inactivity (you can reopen it if you need further help).
-->

**Functionality**
<!-- Add here a clear and concise description of the topic of your request (what you want to do and you do not know how to do it). -->

**Documentation reference**
<!-- Describe here which documentation you checked (e.g. specific pages from https://docs.ropensci.org/MODIStsp/, blog pages, paper at https://doi.org/10.1016/j.cageo.2016.08.020, StackOverflow questions) without finding any answer. -->
---
name: Bug report
about: Report a code bug
title: ''
labels: bug
assignees: ''

---

<!--
Use this template to report a bug. Please use this method instead than sending private email to the authors. In case you are not sure if your error is due to a code bug, please open the "Help needed" template.

Before opening a new issue please check if the problem was already been mentioned (if so but the found issue is closed, open a new issue citing the old task instead than reopening it).

Ensure that your {MODIStsp} version is update with the newest GitHub master branch:
install.packages("remotes")
remotes::install_github("ropensci/MODIStsp")

Please take particular care with code reproducibility (follow indications provided in the template).
-->

**Bug description**
<!-- Add here a clear and concise description of what the bug is. -->

**Reproducible example**
<!-- Please provide here a reproducible example in the chunk below. -->

```r
## PLEASE DELETE AND WRITE YOUR OWN
library(MODIStsp)
MODIStsp(gui = FALSE, opts_file = "/path/of/the/parameter_file.json")
# file parameter_file.json _must_ be attached or copied as text,
# as well as referenced files (e.g. spafile)
```

**Expected and actual behavior**
<!-- Provide here the full output of the provided example and describe what is going wrong- -->

```
## PASTE HERE THE OUTPUT OF YOUR EXAMPLE CODE
```

**System information**
<!-- Provide here the output of the following R commands:
sessionInfo()
packageVersion("MODIStsp")
 -->

```
# PASTE HERE YOUR OUTPUT
```

**Additional context**
<!-- Add here any other context about the problem here (for example, the content of the output folder in case the error appears during a subsequent code execution). -->
