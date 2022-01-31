
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
---
output:
  github_document:
    toc: no
    toc_depth: 2
---


```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```
[![](https://www.r-pkg.org/badges/version-ago/MODIStsp)](https://cran.r-project.org/package=MODIStsp)
[![](https://cranlogs.r-pkg.org/badges/MODIStsp?color=orange)](https://cran.r-project.org/package=MODIStsp)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/MODIStsp.svg?branch=master)](https://travis-ci.org/ropensci/MODIStsp)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1972039.svg)](https://doi.org/10.5281/zenodo.1972039)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/MODIStsp/master.svg)](https://codecov.io/github/ropensci/MODIStsp?branch=master)
[![](https://badges.ropensci.org/184_status.svg)](https://github.com/ropensci/software-review/issues/184)

# MODIStsp <img src='man/figures/logo.png' align="right" height="139" />

[`{MODIStsp}`](https://docs.ropensci.org/MODIStsp/) is a `R` package devoted to
automatizing the creation of time series of rasters derived from MODIS Land Products
data.
`{MODIStsp}` allows performing several preprocessing steps (e.g., download, mosaicing,
reprojection and resize) on MODIS data available within a given time period. 
Users have the ability to select which specific layers of the original MODIS HDF
files they want to process. They also can select which additional Quality Indicators 
should be extracted from the aggregated MODIS Quality Assurance layers and, in the 
case of Surface Reflectance products, which Spectral Indexes should be computed 
from the original reflectance bands. For each output layer, outputs are saved as
single-band raster files corresponding to each available acquisition date. 
Virtual files allowing access to the entire time series as a single file can be 
also created. All processing parameters can be easily selected with a user-friendly
GUI, although non-interactive execution exploiting a previously created Options 
File is possible. Stand-alone execution outside an "R" environment is also possible,
allowing to use scheduled execution of MODIStsp to automatically update time series 
related to a MODIS product and extent whenever a new image is available. 

<a href="http://www.irea.cnr.it/en/"> <img src="man/figures/irea_logo.png" height="40" align="left" /></a> 

<span style='font-style:italic;font-weight:bold;'>`{MODIStsp}` was developed by Lorenzo Busetto and Luigi Ranghetti, 
[Institute of Remote Sensing of Environment](http://www.irea.cnr.it/en/) - National Research Council - Italy (CNR-IREA).
[It is dedicated to the memory of Lorenzo](https://docs.ropensci.org/MODIStsp/articles/lorenzo.html).</span>

## Citation
  
To cite `{MODIStsp}` please use:


L. Busetto, L. Ranghetti (2016) MODIStsp: An R package for automatic preprocessing of MODIS
  Land Products time series, Computers & Geosciences, Volume 97, Pages
  40-48, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2016.08.020, URL: https://github.com/ropensci/MODIStsp. 
  

## Website

For more information, documentation and examples of use, __see also the `{MODIStsp}` website at [docs.ropensci.org/MODIStsp](https://docs.ropensci.org/MODIStsp/)__. 


## Important News

- 29/10/2021 - `{MODIStsp}` (GitHub version 2.0.6.9000) supports products with version 061.
Version 006 will remain the default product version until its decommission
will be announced.
Version 061 can be specified through the argument `prod_version` of function
`MODIStsp()` or by selecting it in the GUI.

- 10/12/2020 - `{MODIStsp}` was resubmitted to CRAN after the maintainer's death.
Now `{MODIStsp}` is dedicated to Lorenzo Busetto (https://docs.ropensci.org/MODIStsp/articles/lorenzo).

- 01/09/2020 - `{MODIStsp}` 2.0.0 is out. Provides a new GUI interface based on Shiny, getting rid
of the archived dependencies on gWidgets/gWidgetsRGtk2. Also provides much easier usage from 
the CLI, by allowing to set all processing arguments also from the CLI. __Note:__ due to the 
introduced changes, options files created with previous versions of `{MODIStsp}` will no 
longer work. Also, processing scripts using `{MODIStsp}` may need to be slightly adapted.

- 09/05/2020 - `{MODIStsp}` 1.4.0 is out. Switches to use of GDAL3/PROJ6 WKTs for projection representation and usage of `{sf}` for all internal work on vector data. Adds support for products MCD19A1 and MCD19A2 products.

- 07/06/2019 - `{MODIStsp}` 1.3.9 is out. Fixes a bug causing crashes on MOD14A1 product, adds support for product MCD12Q2 and removes support for no longer available version 5 of some products.

- 05/03/2019 - `{MODIStsp}` 1.3.8 is out. Fixes an issue causing incorrect application of scale/offset values on GDAL versions > 2.3 (https://github.com/ropensci/MODIStsp/issues/163) and adds support for products MOD21A1D.006, MOD21A1N.006 and MOD21A2.006.

- 29/11/2018 - We recently discovered a nasty bug in the computation of some custom spectral indices (those including additions / subtractions on reflectance values, such as in (b1_NIR+0.1) / b2_Red. See [here](https://docs.ropensci.org/MODIStsp/articles/discovered_bug.html) for further details. The bug is fixed as of version 1.3.7. 

- 07/08/2018 - We are glad to report that `{MODIStsp}` is now included in the 
[rOpenSci](https://ropensci.org/about/) packages ecosystem. We thank reviewers
Leah Wasser and Jeffrey Hanson for their valuable reviews, which helped us to 
further improve the package.

- 10/07/2018 - `{MODIStsp}` v. 1.3.6 is out (check out the [Release Notes](https://github.com/ropensci/MODIStsp/releases/tag/1.3.6) for further details).

- 20/06/2018 - `{MODIStsp}` v. 1.3.5 is out (check out the [Release Notes](https://github.com/ropensci/MODIStsp/releases/tag/v1.3.5) for further details).

- 11/04/2018 - Due to new NASA Policies the MODIS FTP servers were shut down 
starting, April 2, 2018. **FTP download is therefore no longer working** and will
be removed in the next MODIStsp version.

- 11/04/2018 - [**Decommissioning of MODIS Version 5 Land Data Products**]( https://lpdaac.usgs.gov/news/decommissioning-modis-version-51-land-cover-type-data-products-january-7-2019/). As per NASA notice above, MODIS v005 products are going to be 
decommissioned, and will soon be no longer available for download. Support for those
products will be removed in the next MODIStsp version.

- 11/08/2017 - `{MODIStsp}` 1.3.3 was released today. It provides improvements in
processing speed, as well as the usual bug fixes
(thanks to all the users that signaled problems). Check the [Release Notes](https://github.com/ropensci/MODIStsp/releases/tag/v1.3.3) for further details.

-   25/07/2017 - As of today, **most of the content related to `{MODIStsp}` has been
moved to our new website at [docs.ropensci.org/MODIStsp](https://docs.ropensci.org/MODIStsp/)
**, which provides a much better user interface and ease of access to MODIStsp-related 
information. From now on, please **consult the new website for detailed and updated
information on the package**. 

-   Also our previous FAQ page on GitHub containing info for solving common 
installation, downloading and processing problems and issues was discontinued
and **migrated at [docs.ropensci.org/MODIStsp/articles/faq.html](https://docs.ropensci.org/MODIStsp/articles/faq.html)**.  

## Problems and Issues

- Please **report any issues** you may encounter in our [issues page on github <i class="fa fa-github-square" aria-hidden="true"></i>](https://github.com/ropensci/MODIStsp/issues).

## <i class="fa fa-desktop" aria-hidden="true"></i> System Requirements 

`{MODIStsp}` requires [`R`](https://cran.r-project.org) v >= 3.6.3.

____________________________________________________________________________________

# Installation Instructions

## <i class="fa fa-windows" aria-hidden="true"></i> Installing on Windows

You can install the stable version of `{MODIStsp}` from CRAN: 

`install.packages("MODIStsp")`

, or the development version (containing the latest improvements and bug fixes) 
from GitHub:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
install.packages("remotes")
library(remotes)
install_github("ropensci/MODIStsp")
```

## <i class="fa fa-linux" aria-hidden="true"></i> Installing on Linux Systems

To install `{MODIStsp}` on Linux, you need to be able to install the `{sf}` package, 
which requires several dependencies. See [here](https://github.com/r-spatial/sf#installing)
if you have trouble installing `{sf}`. 

Then, you can install the stable version of MODIStsp from CRAN:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
install.packages("MODIStsp")
```
, or the development version (containing the latest improvements and bug fixes) 
from GitHub:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
library(devtools)
install_github("ropensci/MODIStsp")
```

## <i class="fa fa-apple" aria-hidden="true"></i> Installing on Mac


To install `{MODIStsp}` on Mac, you need to be able to install the `{sf}` package, 
which requires several dependencies. See [here](https://github.com/r-spatial/sf#installing)
if you have trouble installing `{sf}`. 

Then, you can install the stable version of `{MODIStsp}` from CRAN:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
install.packages("MODIStsp")
```
, or the development version (containing the latest improvements and bug fixes) 
from GitHub:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
library(devtools)
install_github("ropensci/MODIStsp")
```

# Usage

The easiest way to use `{MODIStsp}` is to use its powerful GUI (Graphical User Interface)
for selection of processing options, and then run the processing. 

To open the GUI, load the package and launch the MODIStsp function, with no parameters:
```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, caption=FALSE}
library(MODIStsp)
MODIStsp()
```
This **opens a Shiny GUI** from which processing options can be specified (and eventually 
saved or loaded). After specifying all required parameters, clicking on "Start" 
will start
the processing (see [here](https://docs.ropensci.org/MODIStsp/articles/interactive_execution.html) 
for more detailed instructions).

`{MODIStsp}` can also be launched in non-interactive mode within an `R` session or
script by setting the optional `GUI` parameter to FALSE, and the `opts_file`
parameter to the path of a previously saved JSON Options file. This allows to 
exploit `{MODIStsp}` functionalities within generic "R" processing scripts.

```{r, eval=FALSE}
library(MODIStsp) 
# --> Specify the path to a valid options file saved in advance from MODIStsp GUI 
opts_file <- "X:/yourpath/youroptions.json" 
  
# --> Launch the processing
MODIStsp(gui = FALSE, opts_file = opts_file)
```

Finally, `{MODIStsp}` can be run by manually specifying all processing arguments, 
or by overwriting some of the arguments contained in a saved json file in the call 
to the package, such as in: 

```{r, eval=FALSE}
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
, where we are overwriting the options related to spatial and temporal extent
contained in the options file with new values. This allows easily running
processing based on the same main options (e.g., product layers, output format, etc.)
but changing on the fly the desired ones. 

See [here](https://docs.ropensci.org/MODIStsp/articles/noninteractive_execution.html) 
for more detailed instructions and examples.


# Code of Conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.


---
output: github_document
---
<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

[![](https://www.r-pkg.org/badges/version-ago/MODIStsp)](https://cran.rstudio.com/web/packages/MODIStsp/index.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.290683.svg)](https://doi.org/10.5281/zenodo.290683)
[![Downloads](https://cranlogs.r-pkg.org/badges/MODIStsp?color=orange)](https://cran.rstudio.com/web/packages/MODIStsp/index.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/MODIStsp.svg?branch=master)](https://travis-ci.org/ropensci/MODIStsp)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/MODIStsp/master.svg)](https://codecov.io/github/ropensci/MODIStsp?branch=master)


# <i class="fa fa-globe" aria-hidden="true"></i> MODIStsp <img src="man/figures/logo.png" width="100" height="100" align="right"/>

<!-- # MODIStsp <img src='man/figures/logo.png' align="right" height="139" /> -->

**`{MODIStsp}`** is a `R` package devoted to automatizing the creation of time 
series of raster images derived from MODIS Land Products data.

**`{MODIStsp}`** allows performing several preprocessing steps (e.g., download, mosaicing, 
reprojection, resize, data extraction) on MODIS data available within a given time period.
Users have the ability to select which specific layers of the original MODIS HDF files they 
want to process. They also can select which additional **Quality Indicators** should be 
extracted from the aggregated MODIS Quality Assurance layers and, in the case of Surface 
Reflectance products, which **Spectral Indexes** should be computed from the original reflectance 
bands. 

All processing parameters can be easily selected with a **powerful and user-friendly GUI**, 
although non-interactive execution exploiting a previously created Options File is possible. 
Stand-alone execution outside an `R` environment is also possible, allowing using scheduled 
execution of MODIStsp to automatically update time series related to a MODIS product and extent 
whenever a new image is available. 

For each output layer, outputs are saved as **single-band raster** files corresponding to 
each available acquisition date. **Virtual files**, allowing accessing to the entire time series 
as a single file, can be also created. 


<a href="http://www.irea.cnr.it/en/"> <img src="man/figures/irea_logo.png" height="40" align="left" /></a>


<span style='font-style:italic;font-weight:bold;'>`{MODIStsp}` was developed by Lorenzo Busetto and Luigi Ranghetti, 
[Institute of Remote Sensing of Environment](http://www.irea.cnr.it/en/) - National Research Council - Italy (CNR-IREA).
[It is dedicated to the memory of Lorenzo](https://docs.ropensci.org/MODIStsp/articles/lorenzo.html).</span>

____________________________________________________________________________________

## <i class="fa fa-newspaper-o" aria-hidden="true"></i> What's New 

- 29/10/2021 - `{MODIStsp}` (GitHub version 2.0.6.9000) supports products with version 061.
Version 006 will remain the default product version until its decommission
will be announced.
Version 061 can be specified through the argument `prod_version` of function
`MODIStsp()` or by selecting it in the GUI.

- 10/12/2020 - `{MODIStsp}` was resubmitted to CRAN after the maintainer's death.
Now `{MODIStsp}` is dedicated to Lorenzo Busetto (https://docs.ropensci.org/MODIStsp/articles/lorenzo).

- 01/09/2020 - `{MODIStsp}` 2.0.0 is out. Replaces the old gWidgets GUI with a new one 
based on Shiny, enhances support for CLI usage and enhances support/provides bug fixing for datasets with 
multiple NoData values when applying scale/offset.

- 09/05/2020 - `{MODIStsp}` 1.4.0 is out. Switches to use of GDAL3/PROJ6 WKTs for projection representation and usage of `{sf}` for all internal work on vector data. Adds support for products MCD19A1 and MCD19A2 products.

- 07/06/2019 - `{MODIStsp}` 1.3.9 is out. Fixes a bug causing crashes on MOD14A1 product, adds support for product MCD12Q2 and removes support for no longer available version 5 of some products.

- 05/03/2019 - `{MODIStsp}` 1.3.8 is out. Fixes an issue causing incorrect application of scale/offset values on GDAL versions > 2.3 (https://github.com/ropensci/MODIStsp/issues/163) and adds support for products `MOD21A1D.006 MOD21A1N.006 MOD21A2.006`.

- 29/11/2018 - We recently discovered a nasty bug in the computation of some custom spectral indices (those including additions / subtractions on reflectance values, such as in $\frac{(b1_{NIR}+0.1)}{b2_{Red}}$ (with $\rho$ being a reflectance). See [here](https://docs.ropensci.org/MODIStsp/articles/discovered_bug.html) for further details. The bug is now fixed on the 
GitHub version. A patched release will be made available on CRAN as soon as possible. 

- 07/08/2018 - We are glad to report that `{MODIStsp}` is now included in the 
[rOpenSci](https://ropensci.org/about/) packages ecosystem. We thank reviewers
Leah Wasser and Jeffrey Hanson for their valuable reviews, which helped us  
further improving the package.

- 10/07/2018 - `{MODIStsp}` v. 1.3.6 is out (check out the [Release Notes](https://github.com/ropensci/MODIStsp/releases/tag/1.3.6) for further details).

- 20/06/2018 - `{MODIStsp}` v. 1.3.5 is out (check out the [Release Notes](https://github.com/ropensci/MODIStsp/releases/tag/v1.3.5) for further details).

- 11/04/2018 - Due to new NASA Policies the MODIS FTP servers were shut down 
starting, April 2, 2018. **FTP download is therefore no longer working** and will
be removed in the next MODIStsp version.

- 11/04/2018 - [**Decommissioning of MODIS Version 5 Land Data Products**]( https://lpdaac.usgs.gov/news/decommissioning-modis-version-51-land-cover-type-data-products-january-7-2019/). As per NASA notice above, MODIS v005 products are going to be 
decommissioned, and will soon be no longer available for download. Support for those
products will be removed in the next MODIStsp version.

- 11/08/2017 - `{MODIStsp}` 1.3.3 was released today. It provides improvements in processing speed, as well as the usual bug fixes.  See our [<i class="fa fa-newspaper-o aria-hidden="true"></i> news](news/index.html) page for a detailed changelog.

____________________________________________________________________________________

## <i class="fa fa-cog" aria-hidden="true"></i> Getting Started

- __To install `{MODIStsp}`__, please follow instructions reported [here](articles/installation.html),
both for [<i class="fa fa-windows" aria-hidden="true"></i>](articles/installation.html#installing-on-windows) , [<i class="fa fa-linux" aria-hidden="true"></i>](articles/installation.html#installing-on-linux-systems) and [<i class="fa fa-apple" aria-hidden="true"></i>](articles/installation.html#installing-on-mac).

- **`{MODIStsp}`** can be used either in [interactive mode](articles/interactive_execution.html) 
exploiting its user-friendly GUI, or in [non-interactive mode](articles/noninteractive_execution.html) 
from within `R` scripts.

- The list of __currently supported MODIS products and versions__ can be found [here](articles/products_list.html).

- [Scheduled Processing](articles/noninteractive_execution.html#scheduled-processing) allows 
automatic updating of MODIS time series through scheduled jobs, both on [<i class="fa fa-windows" aria-hidden="true"></i>](articles/standalone_execution.html#on-windows) and [<i class="fa fa-linux" aria-hidden="true"></i>](articles/standalone_execution.html#on-linux).

- Solutions to common **installation, data download and processing problems** can be found 
in our [<i class="fa fa-question-circle-o" aria-hidden="true"></i> faq](https://docs.ropensci.org/MODIStsp/articles/faq.html).

- Please **report any issues** you may encounter in our [issues page on github <i class="fa fa-github-square" aria-hidden="true"></i>](https://github.com/ropensci/MODIStsp/issues).

____________________________________________________________________________________

## <i class="fa fa-pencil" aria-hidden="true"></i>Citation


To cite MODIStsp in publications, please use:

L. Busetto, L. Ranghetti (2016) MODIStsp: An R package for automatic preprocessing of MODIS
  Land Products time series, Computers & Geosciences, Volume 97, Pages
  40-48, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2016.08.020, URL: https://github.com/ropensci/MODIStsp. 
  
A BibTeX entry for LaTeX users is:

```
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
```

---
title: "Accessing and Analyzing Processed Data from R"
output: 
  github_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(rmarkdown)
library(knitr)
library(raster)
library(base)
library(sf)
library(MODIStsp)
library(dplyr)
test_folder <-  system.file("testdata/VI_16Days_500m_v6/NDVI",
                            package = "MODIStsp")
dir.create(file.path(tempdir(), "MODIStsp/VI_16Days_500m_v6/NDVI"), 
           showWarnings = FALSE, recursive = TRUE)
file.copy(list.files(test_folder, full.names = TRUE),
          file.path(tempdir(), "MODIStsp/VI_16Days_500m_v6/NDVI"), 
          recursive = T)
Sys.setlocale("LC_TIME", "en_US.UTF-8")
```

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Accessing the raster time series__

Preprocessed MODIS data can be retrieved within `R` either by accessing the 
single-date raster files, or by loading the saved `RasterStack` objects (see 
[here](output.html) for a description of `{MODIStsp}` output folder
structure and naming conventions).

To test this functionality, you can run the following example: 

```{r eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE}
library(MODIStsp)
opts_file <- system.file("testdata/test_extract.json", package = "MODIStsp")
MODIStsp(opts_file = opts_file, gui = FALSE, verbose = FALSE)
```

This will download a yearly time series of MODIS NDVI data and subset it over the
region of the Como Lake in Italy.

After the download and processing finishes (it will take a while, depending on your
network speed), the MODIS time series will be placed in subfolder `MODIStsp/VI_16Days_500m_v6` 
of `R` `tempdir()`. In this case, I have them in "`file.path(tempdir(), "VI_16Days_500m_v6")`". 
If you want to save them elsewhere, run `MODIStsp(opts_file = opts_file, gui = TRUE)`
and select a different output folder.

Single-date processed rasters can be accessed by simply opening them with a
`raster` command: 

```{R echo=TRUE, message=FALSE, warning=FALSE, highlight=TRUE, tidy=TRUE}
modistsp_file <- file.path(tempdir(),"MODIStsp/VI_16Days_500m_v6/NDVI",
                           "MOD13A1_NDVI_2016_177.tif")
modistsp_file

my_raster     <- raster(modistsp_file)
plot(my_raster)
```

`RasterStack` (or GDAL vrt) time series containing all the processed data for a
given parameter are saved in the `Time Series/RData` subfolder, and can be accessed
by: 

```{r echo=TRUE, message=FALSE, warning=FALSE}
# Load the NDVI time series
ts_folder <- file.path(tempdir(), "MODIStsp/VI_16Days_500m_v6/Time_Series")
in_virtual_file <- file.path(ts_folder, "RData/Terra/NDVI/MOD13A1_NDVI_1_2016_353_2016_RData.RData")
in_virtual_file

ts_data   <- get(load(in_virtual_file))
ts_data
```

, which gives us a 23-band `RasterStack` in the `ts_data` variable.

This `RasterStack` can be analyzed using the functionalities for raster/raster 
time series analysis, extraction and plotting provided for example by the 
``{raster}`` or   ``{rasterVis}`` packages: 

```{r echo=TRUE, message=FALSE, warning=FALSE}
# plot some dates
plot(ts_data[[c(1,5,10,15,20)]], axes = FALSE, horizontal = T)

# Extract one date from the stack

mydate   <- as.Date("2016-01-01")
substack <- subset(ts_data, which(getZ(ts_data) == mydate)) %>% setZ(mydate)
substack  

# Extract multiple dates from the stack

mindate   <- as.Date("2016-01-01")
maxdate   <- as.Date("2016-04-01")
substack  <- subset(ts_data, 
                    which(getZ(ts_data) >= mindate & getZ(ts_data) <= maxdate))
substack  

# Compute monthly averages

month_avg <- stackApply(ts_data, months(getZ(ts_data)), fun = mean)
month_avg

plot(month_avg, axes = FALSE, horizontal = T)
```

________________________________________________________________________________

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Extracting Time Series Data on Areas of Interest__

`{MODIStsp}` provides an efficient function (`MODIStsp_extract()`) for extracting 
time series data at specific locations. The function takes as input a `RasterStack` 
virtual object created by `MODIStsp()` (see above), the starting and ending dates 
for the extraction and a standard `{sp}` object (or an ESRI shapefile name) 
specifying the locations (points, lines or polygons) of interest, and provides 
as output a `xts` object or `data.frame` containing time series data for those
locations. 

If the input is of class `SpatialPoints`, the output object contains one column 
for each point specified, and one row for each date. If it is of class
`SpatialPolygons` (or `SpatialLines`), it contains one column for each polygon 
(or each line), with values obtained applying the function specified as the "FUN" 
argument (e.g., mean, standard deviation, etc.) on pixels belonging to the polygon
(or touched by the line), and one row for each date. 

To test `MODIStsp_extract()` we can start by loading the NDVI `RasterStack`
time series created in the previous example:

```{r echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE}

ts_folder <- file.path(tempdir(), "MODIStsp/VI_16Days_500m_v6/Time_Series")
in_virtual_file <- file.path(ts_folder, "RData/Terra/NDVI/MOD13A1_NDVI_1_2016_353_2016_RData.RData")
ts_data   <- get(load(in_virtual_file))
ts_data
```

, which again gives us the 23-band `RasterStack` in the `ts_data` variable. 

To extract the NDVI data over selected areas, we can then use the `MODIStsp_extract()`
function as follows: 

```{r echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE}

# load a shapefile containing polygons --> Here we use a test shapefile in 
# the testdata folder of MODIStsp.

polygons <- sf::st_read(system.file("testdata/extract_polys.shp",
                                       package = "MODIStsp"),
                           quiet = TRUE)
polygons

# Now extract the average values for each polygon and date and plot the results

out_dataavg <- MODIStsp_extract(ts_data, polygons, id_field = "lc_type", 
                                small = FALSE)

head(out_dataavg)

# Other summarization functions can be used, by specifying the "FUN" argument. 
# Compute the Standard Deviation over each polygon: 
out_datasd <- MODIStsp_extract(ts_data, polygons, id_field = "lc_type",
                               FUN = "sd", small = FALSE)

```

The output is a `xts` object, with one column for each polygon of the input 
shapefile, and one row for each date. We can easily plot the computed averages 
using: 

```{r echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE}
library(xts)
plot(out_dataavg, legend.loc = "topleft")
```

We can also transform the output dataset to a long format and use `{ggplot2}` for
plotting and other tidyverse-friendly functions for analysis:

```{r echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE}
library(ggplot2)
out_dataavg_df <- data.frame(date = index(out_dataavg), coredata(out_dataavg)) %>% 
  tidyr::gather(Polygon, Value, -date)

ggplot2::ggplot(out_dataavg_df, aes(x = date, y = Value, color = Polygon)) + 
  geom_line() + theme_light()

```

```{r echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE}

require(tibbletime)
require(dplyr)

# Compute monthly averages using tibbletime: 

out_dataavg_tt <- tibbletime::as_tbl_time(out_dataavg_df, index = date) 

month_avg <- out_dataavg_tt %>% 
  dplyr::group_by(Polygon) %>% 
  tibbletime::as_period(period = "monthly")

ggplot2::ggplot(month_avg, aes(x = date, y = Value, color = Polygon)) + 
  geom_line() + theme_light()

```

**Important note:** `MODIStsp_extract()` is usually much faster than the 
standard `raster::extract()` function, but does not deal with overlapping polygons.
If your polygons are overlapping, please use `raster::extract()` instead. 
---
title: "Running the tool in Interactive Mode: the MODIStsp GUI"
bibliography: MODIStsp.bib
output: 
  github_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(DT)
```

The easiest way to use `{MODIStsp}` is to use its powerful GUI (Graphical User Interface) 
for selection of processing options, and then run the processing. 

To open the GUI, load the package and launch the `{MODIStsp}` function, with no parameters:
```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, caption=FALSE}
library(MODIStsp)
MODIStsp()
```
This **opens a GUI** from which processing options can be specified and eventually 
saved (or loaded from a previously saved file). 

The  GUI allows selecting all processing options required for the creation of the 
desired MODIS time series. The GUI uses a dashboard structure, divided in the 
following tabs. The available processing options configurable in each tab are described 
in the following.

____________________________________________________________________________________


# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Selecting Processing Parameters__

<br>

## _Product and Layers Tab_


```{r GUIfig, echo=FALSE, fig.align="center", fig.width=10, message=FALSE, warning=FALSE}
  library(png)
  library(grid)
  library(knitr)
  img <- readPNG("GUI_1.PNG")
  grid.raster(img)
```

The top-most tab allow to specify details of the desired MODIS Product and Layers
to be processed:

1. **"Category"** and **"Product"**: selects the MODIS product of interest;
2. **MODIS platform(s)**: selects if only TERRA, only AQUA or Both MODIS platforms 
should be considered for download and creation of the time series;
3. **MODIS layers to be processed**: the user **must** 
select which MODIS original layers and/or derived Quality Indexes (QI) and Spectral
Indexes (SI) layers should be processed: 
    - the left-hand selector allows to select which _original MODIS layers_ 
        should be processed;
    - the central selector allows to select which _Quality Indicators should be extracted_ 
        from the original MODIS Quality Assurance layers;
    - for MODIS products containing surface reflectance data, the right-hand selector 
        allows selecting which additional _Spectral Indexes should be computed_. 

The following commonly used Spectral Indexes are available for computation by default: 

<br>

```{r xtable, echo=FALSE, paged.print=TRUE, results="asis"}
# library(xtable)
tab <- tibble::tribble(
 ~"Acronym"     ,~"Index Name and reference", ~"Index Formula",                     
 "NDVI"         , "Normalized Difference Vegetation Index (Rouse, 1973)"          , "(NIR - RED)/(NIR + RED)",        
 "EVI"          , "Enhanced Vegetation Index (Huete, 2002)"                       , "2.5 * (NIR - RED)/(NIR + 6 * RED - 7.5 * BLUE + 1",
 "SR"           , "Simple Ratio[@Tucker1979]"                                    , "NIR / RED",
 "NDFI"         , "Normalized Difference Flood Index (Boschetti, 2014)"           , "(NIR - SWIR1) / (NIR + SWIR1)", 
 "NDII6 (NDWI6)" , "Normalized Difference Infrared Index - Band 6 (Hunt, 1989)" , "(NIR - SWIR1) / (NIR + SWIR1)",
 "NDII7 (NDWI7)" , "Normalized Difference Infrared Index - Band 7 (Hunt, 1989)" , "(NIR - SWIR2) / (NIR + SWIR2)",
 "SAVI"         , "Soil Adjusted Vegetation Index  (Huete, 1988)"                 , "((NIR - RED) / (NIR + RED + 0.5)) * (1 + 0.5)",
 "NDSI"         , "Normalized Difference Snow Index (Hall, 2002)"                 , "(GREEN - SWIR1) / GREEN + SWIR1)",
 "GNDVI"        , "Green Normalized Difference Vegetation Index (Gitelson, 1998)" ,  "(NIR - GREEN)/(NIR + GREEN)",       
 "RGRI"         , "Red Green Ratio Index (Gamon, 1999)"                          , "RED / GREEN",
 "GRVI"         , "Green-Red ratio Vegetation Index  (Tucker, 1979)"              , "(RED - GREEN) / (RED + GREEN)"       
)

DT::datatable(tab, rownames = FALSE, style = "bootstrap", 
              options = list(dom = 'tip', pageLength = 11))
```

You can however **specify other SIs to be computed without modifying MODIStsp source code** 
by clicking on the _**"Add New Spectral Index"**_ button, which allow providing info related 
to the new desired SI using a simple GUI interface. 

```{r indexfig, echo=FALSE, message=FALSE, fig.width=6, warning=FALSE, fig.align="center"}
  library(png)
  library(grid)
  img <- readPNG('GUI_newind.PNG')
  grid.raster(img)
```

Provided information (e.g., correct band-names, computable formula, etc...) is 
automatically checked upon clicking "Set New Index". On success, the new index 
is added in the list of available ones for all products allowing its computation. 
Clicking "Done!" returns to the main.

__Note:__ all custom defined indexes can be removed by using the `MODIStsp_resetindexes()` 
function.

## _Spatial/Temporal Options Tab_

The middle tab allow specifying details about the temporal and spatial extent of 
the analysis. 

```{r , echo=FALSE, message=FALSE, warning=FALSE, fig.width=10, fig.align="center"}
  library(png)
  library(grid)
  img <- readPNG('GUI_2_tiles.PNG')
  grid.raster(img)
```

### _Temporal Extent_

Specify the starting and ending dates to be considered for the creation of the 
time in the series corresponding fields. 

The **Date Range Type** drop-down menu allows to choose between two options:

1.  **full**: all available images between the starting and ending dates are 
downloaded and processed;

2.  **seasonal**: data is downloaded only for one part of the year, but for 
multiple years. For example, if the starting date is 2005-03-01 and the ending is
2010-06-01, only the images of March, April and May for the years between 2005
and 2010 will be downloaded. This allows to easily process data concerning a 
particular season of interest.

### _Output Projection_

Specify the options to be used for reprojecting and resizing the MODIS images. 

- **"Output Projection"**: select either the Native MODIS projection (Default) 
or specify a user-defined one. To specify a user selected projection, select 
"Change" and then insert a valid "EPSG" code or WKT string in the pop-up window. 
Validity of the new projection string is automatically checked, and error messages issued 
if the check fails.

- **"Output Resolution"**, **"Pixel Size"** and **"Resampling Method"**: specify 
whether output images should inherit their spatial resolution from the original
MODIS files, or be resampled to a user-defined resolution. In the latter case, 
output spatial resolution must be specified in the measure units of the selected
output projection. Resampling method can  be chosen among the ones available for
the "gdalwarp" routine. 
__Note:__ resampling methods different than Nearest Neighbour" and "Mode" (Useful for down-sampling purposes) should used carefully. Other resampling methods (e.g., 
bilinear, cubic) i) cannot be used for resampling of categorical variables such as the QA and QI layers, and ii) using them on continuous variable (e.g., reflectance, VI values) without performing an 
a-priori data cleaning would risk to contaminate the values of high-quality
observations with those of low-quality ones.

### _Spatial Extent_

Allows defining the area of interest for the processing. Four main options are
possible, and can be selected using a dropdown menu. 

1.  **Select Tiles**: specify which MODIS tiles need to be processed either 
by: 

    a. Using the "Start" and "End" horizontal and vertical sliders in the 
    _Required MODIS Tiles_ frame.  
    b. Selecting the __"From Map"__ option in the dropdown and clicking on "change selection". 
    A map will open, allowing interactive selection of the required tiles
    
Note that during processing, data from the different tiles is mosaicked, and a single file 
covering the total area is produced for each acquisition date. For this reason, 
selected tiles must cover a rectangular area. 

2.  **Select Bounding Box**: manually insert the coordinates of the Upper Left and Lower Right corners
    of the area of interest in the __Bounding Box__ frame. _Coordinates of the corners 
    must be provided in the coordinate system of the selected output projection_.

2.  **Load From File**: click the __"Browse"__ button and select a raster or vector 
    spatial file that will be used to compute the required bounding box in output
    projection coordinates.

4.  **Draw on Map**: click the __"Draw Extent"__ button a map will open in a  
    window, allowing interactive selection of the spatial extent using the tools
    on the left.

## _Output Format, Options and Folders Tab_

The last tab allows specifying some options concerning processing, and the output folders. 


```{r , echo=FALSE, message=FALSE, warning=FALSE, fig.width=10, fig.align="center"}
  library(png)
  library(grid)
  img <- readPNG('GUI_3.PNG')
  grid.raster(img)
```

### _Download Method_

Select the method to be used for download. Available choices are: 

1.  **http**: download through http from NASA lpdaac http archive (http://e4ftl01.cr.usgs.gov). 
This requires providing a user name and password, which can be obtained by registering 
an account at the address [https://urs.earthdata.nasa.gov/profile](https://urs.earthdata.nasa.gov/profile);

2.  **offline**: this option allows to process/reprocess HDF files already available 
on the user's PC without downloading from NASA -- useful if the user already has an
archive of HDF images, or to reprocess data already downloaded via `MODIStsp()` to create
time series for an additional layer (_it is fundamental that the HDFs are those 
directly downloaded from NASA servers_; see [here](https://docs.ropensci.org/MODIStsp/articles/faq.html#working-with-already-downloaded-hdf-files) 
for additional details). 

A second dropdown menu allows selecting if using standard http download to access
NASA servers, or using the **aria2c** downloader, which may speed-up the download. This requires however that that the "aria2c" software is installed in your system. To download and install it, see [aria2.github.io](https://aria2.github.io/).

### _Output Options_

Several processing options can be set using check-boxes/dropdowns:

- **Output Files Format**: two of the most commonly formats used in remote
sensing applications are available at the moment: ENVI binary and GeoTiff. If
GeoTiff is selected, the type of file compression can be also specified among
"None", "PACKBITS", "LZW" and "DEFLATE".

- **Save Time Series as**: specify if virtual multitemporal files should be 
created. These virtual files allow access to the entire time series of images as 
a single file without the need of creating large multitemporal raster images.
Available virtual files formats are "R" rasterStacks, ENVI meta-files and GDAL 
"vrt" files. In particular, `R` RasterStacks may be useful in order to easily 
access the preprocessed MODIS data within `R` scripts 
(see also [here](output.html)).

- **Apply Scale/Offset**: specify if scale and offset values of the different 
MODIS layers should be applied. If selected, outputs are appropriately rescaled
on the fly, and saved in the true "measure units" of the selected parameter (e.g., 
spectral indexes are saved as floating point values; Land Surface Temperature is 
saved in degrees Kelvin, etc.). 

- **Modify No Data Values**: specify if NoData values of MODIS layers should be kept
at their original values, or changed to those specified within the `MODIStsp_Products_Opts`
XML file. By selecting "Yes" in the "Change Original NoData values" check-box, 
NoData of outputs are set to the largest integer value possible for the data type 
of the processed layer (e.g., for 8-bit unsigned integer layers, NoData is set 
always to 255, for 16-bit signed  integer layers to 32767, and  for 16-bit unsigned
integer layers to 65535). Information about the new NoData values is stored both 
in the output rasters, and in the XML files associated with them. __Note:__ some 
MODIS layers have multiple NoData (a.k.a. _fill_) values. if _Modify No Data_ is 
set to "Yes", `MODIStsp()` will convert all _fill_ values to a common output NoData
value.

### _Output Folders_

#### _Main MODIStsp Output Folder_

Select the main folder where the pre-processed time series data will be stored. 
All `MODIStsp()` outputs **will be placed in specific sub-folders of this main folder** 
(see [here](output.html) for details on `MODIStsp` naming conventions).

The **"Reprocess"** selector allows to decide if images already 
available should be reprocessed if a new run of `MODIStsp()` is launched with the same 
output folder. If set to "No", `MODIStsp()` skips dates for which output files following
the `{MODIStsp}` naming conventions are already present in the output folder. This
allows to incrementally extend MODIS time series without reprocessing already available
dates. 

#### _Output Folder for storage of original MODIS HDF_

Select the folder where downloaded **original MODIS HDF files** downloaded from 
NASA servers will be stored. 

The **"delete HDF"** selector allows also to decide if the
downloaded images should be deleted from the file system at the end of the 
processing. To avoid accidental file deletion, this is always set to "No" by default, 
and a warning is issued before execution whenever the selection is changed to "Yes".

____________________________________________________________________________________

<br>

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Saving and Loading Processing Options__


```{r , echo=FALSE, message=FALSE, warning=FALSE, fig.width=6, fig.align="center"}
  library(png)
  library(grid)
  img <- readPNG('GUI_bar.PNG')
  grid.raster(img)
```

Specified processing parameters can be saved to a JSON file for later use by clicking
on the _**Save Options**_ button in the sidebar.

Previously saved options can be restored clicking on the _**Load Options**_ button
and navigating to the previously saved JSON file.

____________________________________________________________________________________

<br>

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Starting the processing__

Once you are happy with your choices, click on **Run MODIStsp**. `MODIStsp()` 
will start accessing NASA servers to download and process the MODIS data corresponding 
to your choices.

For each date of the specified time period, `MODIStp()` downloads and preprocesses
all HDF images required to cover the desired spatial extent. Informative messages
concerning the status of the processing are provided on the console, as well as on
a self-updating progress window. 

The processed time series are saved in specific subfolders of the main selected
output folder, as explained in detail [here](output.html).

________________________________________________________________________________
---
title: "Stand-alone execution and scheduled processing"
output: 
  github_document: default
---

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Stand-alone execution__ 
`{MODIStsp}` can be also executed as a "stand-alone" application (i.e., without having 
to open R/RStudio). 

To be able to do that, from R launch the function `install_MODIStsp_launcher()`.

- In a Linux operating system this function creates a desktop entry (accessible from 
the menu in the sections "Science" and "Geography") and a symbolic link in a known path
(default: `/usr/bin/MODIStsp`).
- In Windows, a link in the Start Menu and optionally a desktop shortcut are created.
See `?install_MODIStsp_launcher` for details and path customization.

Double-clicking these files or launching them from a shell without parameters will
launch `MODIStsp()` in interactive mode, by opening the GUI. Non-interactive mode 
can be triggered by adding the "-g"  argument to the call, and specifying the path
to a valid Options File as `-s` argument:

  - Linux: `MODIStsp -g -s "/yourpath/youroptions.json"`
  (see `MODIStsp -h` for details).
  
  - Windows:`your_r_library\MODIStsp\ExtData\Launcher\MODIStsp.bat -g -s "yourpath/youroptions.json"`
  (see `C:\Users\you\Desktop\MODIStsp -h` for details).

If you do not want to install any link, launchers can be found in the subdirectory
`MODIStsp/ExtData/Launcher` of your library path.

____________________________________________________________________________________

<br>

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Scheduled Processing__

Stand-alone non-interactive execution can be exploited to periodically and automatically
update the time series of a selected product over a given study area. To do that, you
should simply follow these steps.

1.	Open the `MODIStsp()` GUI, define the parameters of the processing specifying a 
date in the future as the "Ending Date" and save the processing options; then quit
the program.
 
2. schedule non-interactive execution of the `{MODIStsp}` launcher script installed 
as seen before (or located in the subdirectory `MODIStsp/ExtData/Launcher` of your
`R` library path) as a windows scheduled task or linux "cron" job according to a 
specified time schedule, specifying the path of a previously saved Options file
as additional argument.

## _On Linux_

3. Edit your crontab by opening a terminal and typing:

    ```bash
    crontab -e
    ```
 
4. Add an entry for execution of MODIStsp launcher. For example, if you have 
installed it in /usr/bin and you want to run the tool every day at 23.00, add the
following row:
        
    ```bash
      0 23 * * * /bin/bash /usr/bin/MODIStsp -g -s "/yourpath/youroptions.json"
    ```
      
## _On Windows_

3. Create a Task following <a href="https://docs.microsoft.com/en-us/previous-versions/windows/it-pro/windows-server-2008-R2-and-2008/cc748993(v=ws.11)?redirectedfrom=MSDN" target="_blank">this instructions</a>.

4. Add the path of the MODIStsp.bat launcher as Action (point 6), 
specifying  `-g -s "X:/yourpath/youroptions.json"` as argument.
---
title: "Outputs Format and Naming Conventions"
output: 
  github_document: default
---

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Single-band outputs__

Output raster files are saved in specific subfolders of the main output folder. 
In particular, **a separate subfolder** is created for each processed original 
MODIS layer, Quality Indicator or Spectral Index. Each subfolder contains one image
for each processed date, created according to the following naming conventions: 

```
myoutfolder/"ProdName"/"Layer"/"ProdCode"_"Layer"_"YYYY"_"DOY"."ext"
```

, where: 

  - **_ProdName_** is a short name describing the MODIS product from which the 
    datasets were derived (e.g., VI_16Days_1Km_v6);
  - **_Layer_** is a short name describing the dataset (e.g., b1_Red, NDII, UI);
  - **_ProdCode_** is the code name of the MODIS product from which the image was
  derived (e.g., MOD13Q1 - Note that if you choose to process both "Terra" and 
    "Aqua" data, data coming from either platform will be placed in the same folder
    and you will have for example both "MOD13Q1..." and "MYD13Q1..." files in the 
    output folder);
  - **_YYYY_** and **_DOY_** correspond to the year and DOY (Day of the Year) of 
  acquisition of the original MODIS image;
  - **_ext_** is the file extension (.tif for GTiff outputs, or .dat for ENVI outputs). 
  
So, for example, if you process layers "NDVI" and "EVI" for MODIS product MOD13A2 
you will find the resulting GTiff of ENVI output single-date rasters in:

```
/your_out_folder/VI_16Days_1Km_v6/NDVI/
```

and

```
/your_out_folder/VI_16Days_1Km_v6/EVI/
```

____________________________________________________________________________________

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Virtual multi-band outputs__

ENVI and/or GDAL virtual time series files and _RasterStack_ RData objects are 
instead stored **in the "Time\_Series" subfolder** if required.

Naming conventions for these files is as follow:

```
<path_of_out_folder>/Time_Series/<vrt_type>/<Sensor>/<Layer>/<ProdCode>_<Layer>_<StartDOY>_<StartYear>_<EndDOY>_<EndYear>_<suffix>.<ext> 
```

, where: 

  - `<vrt_type>` indicates the type of virtual file (`"RData"`, `"GDAL"` or `"ENVI_META"`);
  - `<Sensor>` indicates to which MODIS sensor the time series belongs (`"Terra"`,
  `"Aqua"`, `"Mixed"` or `"Combined"` (for MCD* products));
  - `<Layer>` is a short name describing the dataset (e.g., `"b1_Red"`, `"NDII"`, `"UI"`);
  - `<ProdCode>` is the code name of the MODIS product from which the image was 
  derived (e.g., `"MOD13Q1"`);
  - `<StartDOY>`, `<StartYear>`, `<EndDOY>` and `<EndYear>` indicate the 
  temporal extent of the time serie created; 
  - `<suffix>` indicates the type of virtual file (`"ENVI"`, `"GDAL"` or `"RData"`);
  - `<ext>` is the file extension (`".vrt"` for GDAL virtual files, `"META"` for 
  ENVI meta files or `"RData"` for `R` raster stacks). 

So, for example, if you process layer "NDVI" for MODIS product MOD13A2 and platform 
"Terra", and ask to create "R rasterStacks" and "GDAL vrt" time series you will 
find the resulting GTiff virtual files in:

```
/your_out_folder/VI_16Days_1Km_v6/Time_Series/RData/Terra/NDVI
```

and 

```
/your_out_folder/VI_16Days_1Km_v6/Time_Series/GDAL/Terra/NDVI
```
---
title: Installation
output: 
  github_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## <i class="fa fa-windows" aria-hidden="true"></i> Installing on Windows

You can install the stable version of `{MODIStsp}` from CRAN: 

`install.packages("MODIStsp")`

, or the development version (containing the latest improvements and bug fixes) 
from GitHub:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
install.packages("remotes")
library(remotes)
install_github("ropensci/MODIStsp")
```

## <i class="fa fa-linux" aria-hidden="true"></i> Installing on Linux Systems

To install `{MODIStsp}` on Linux, you need to be able to install the `{sf}` package, 
which requires several dependencies. See [HERE](https://github.com/r-spatial/sf#installing)
if you have trouble installing `{sf}`.

In addition, you need to install dependencies
required by the `{protolite}` package, required by `{geojson}`. See [HERE](https://github.com/jeroen/protolite/) for instructions
on installing them. 

Then, you can install the stable version of `{MODIStsp}` from CRAN:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
install.packages("MODIStsp")
```
, or the development version (containing the latest improvements and bug fixes) 
from GitHub:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
library(devtools)
install_github("ropensci/MODIStsp")
```

## <i class="fa fa-apple" aria-hidden="true"></i> Installing on Mac

To install `{MODIStsp}` on MacOS, you need to be able to install the `{sf}` package, 
which requires gdal to be installed. See [HERE](https://github.com/r-spatial/sf#installing)
if you have trouble installing `{sf}`. 

Then, you can install the stable version of `{MODIStsp}` from CRAN:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
install.packages("MODIStsp")
```
, or the development version (containing the latest improvements and bug fixes) 
from GitHub:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
library(devtools)
install_github("ropensci/MODIStsp")
```
---
title: "Frequently Asked Questions"
output: 
  github_document:
    toc_depth: 1
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Here you can find possible **solutions for common `MODIStsp` problems**. We will add to this page as soon as additional problems are reported. If these solutions does not work for you, or you have a different problem, signal it at: https://github.com/ropensci/MODIStsp/issues

<a name="Installation Problems"> </a>

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> Installation Problems 
    
-   [On Linux, installation fails due to some packages missing](#Missing-packages)
-   [On Windows, the program hangs at first interactive execution, while attempting to install `gWisgets2RGtk2`](#gWidgets-error)
-   [`MODIStsp` hangs at first execution, while looking for `gdal` installation](#gdal-error)
-   [Installing from github fails due to missing dependencies](#github-error)

<a name="gWidgets-error"></a> 

#### - On Windows, `MODIStsp` fails to install, signaling problems in installation of either `gwidgetsRGtk2` or `RGtk2`   

 Please try to install package `gwidgetsRGtk2` beforehand using: 
  
```
install.packages("gwidgetsRGtk2")
```

You will probably get some errors/problems during installation: please see instructions on the next point !

#### - On Windows, `MODIStsp` hangs at first interactive execution, while attempting to install `gWisgetsRGtk2`

 At first interactive execution (i.e., with "gui = TRUE) of `MODIStsp`, an error window may appear, signaling that  _libatk-1.0-0.dll_ is missing from your system. This is due to the fact that library "GTK+"" is not yet installed on your system and needs to be installed. To do so, press "OK". A new window dialog window will appear, asking if you want to install "GTK+". Select *"Install GTK+"* and then *"OK"*. Windows will download and install the GTK+ library. When it finishes, the R Session should be restarted and next `MODIStsp` should go well ! In case RStudio does not automatically restart after installing GTK+, simply kill it form "Task Manager" and reload RStudio.
 
 If you still have problems, try to install the GTK+ library independently, following instructions reported here:

https://www.stat.auckland.ac.nz/%7Ekimihia/rgtk2

, starting from "Download the GTK+ all-in-one bundle. I downloaded version 2.22 for Windows 64-bit.". Then try to install `gWisgetsRGtk2` before `MODIStsp` as explained in the point above.
 
<a name="Missing-packages"></a> 

#### - On Linux, installation fails due to some packages missing

Please Install the following required dependencies: 
    * Cairo >= 1.0.0, ATK >= 1.10.0, Pango >= 1.10.0, GTK+ >= 2.8.0, GLib >= 2.8.0 (required by package ```RGtk2```)
    * Curl (required by package ```curl```)
    * GDAL >= 1.6.3, PROJ.4 >= 4.4.9 (required by package ```rgdal```)
    
On **Debian and Ubuntu-based systems**, to install packages open a terminal and type:  
```
  sudo apt-get install r-cran-cairodevice r-cran-rgtk2 libcairo2-dev libatk1.0-dev libpango1.0-dev 
  libgtk2.0-dev libglib2.0-dev libcurl4-openssl-dev libgdal-dev libproj-dev
```

On **rpm-based systems**, to install packages open a terminal and type;  
```
   sudo yum install libcairo2-devel libatk1.0-devel libpango1.0-devel gtk2 gtk2-devel 
   glib2-devel libcurl4-devel gdal-devel proj-devel
```
On **others distros**: we did not test installation on other distros yet - if you experience problems please contact us.


<a name="gdal-error"></a>

#### - MODIStsp hangs at first execution, while looking for `gdal` installation

At the first execution `MODIStsp` searches for a valid `GDAL` installation. If nothing happens for a long time (e.g., several minutes), `MODIStsp` (and in particular the `gdalUtils` package on which it relies) is not finding a valid `GDAL` installation in the more common locations. To solve the problem:

1. Ensure that `GDAL` is properly installed in your system. See the main `MODIStsp` github page for simple instructions
2. (On Windows) If it is installed, verify that `GDAL` is in your system PATH, and that the _GDAL\_DATA_ environment variable is correctly set (You can find simple instructions [HERE](http://gisforthought.com/setting-up-your-gdal-and-ogr-environmental-variables/)) (If gdal is not correctly installed or the path not set, then opening a windows shell ("cmd") and issuing the "gdalinfo" command will result in an error !)

If nothing works, please report the issue here: https://github.com/ropensci/MODIStsp/issues


<a name="github-error"></a>

#### - Installing from github fails due to missing dependencies

There are currently some problems in installing `MODIStsp` **via `install_github`** on R >= 3.3.1 due to not correct installation of dependencies (related to a bug in CRAN version of `install_github`). Installing the development version of `devtools` should solve the issue. To do so, on a **clean** R/RStudio session do:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
  install.packages("devtools")
  devtools::install_github("hadley/devtools")
  library(devtools)
```

, then continue with standard `MODIStsp` installation. 

If you have problems in installing the "devel" version of `devtools`, manually installing all the dependencies should also solve the issue. To do so, please try doing: 
  
  ```
  install.packages(c("bitops", "data.table" , "gdalUtilities", "gWidgets", "gWidgetsRGtk2",
    "httr" , "jsonlite", "parallel", "raster", "sf", "stringr", "xts"))
  ```

 , then continue with standard `MODIStsp` installation.
 
If you still don't succeed, please contact us ! 

____________________________________________________________________________________

<a name="download-errors"></a>

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Data download Problems__

-   [On 'http' download, files are not getting downloaded, OR you get 'timeout' problems](#timeout-error)


<a name="timeout-error"></a>

#### -  On 'http' download, files are not getting downloaded, OR you get 'timeout' problems

1. Visit your [earthdata "profile" page](https://urs.earthdata.nasa.gov/profile), click on "My Applications" and *ensure that "LP DAAC Data Pool" is authorized*. If not, click on "Approve More Applications", search for it in the list and approve it.
2. Verify that the **username and password** you provided are correct (Those can be obtained by **registering an account** at: https://urs.earthdata.nasa.gov/profile.);
3. _Is it Wednesday ???_ If so, NASA http server may be down for maintenance. **Try switching to "ftp" download**;
3. In some cases, access to the http server seems to be not allowed (we don't know why - maybe firewalling). **Try switching to ftp download**.

#### - Files are not getting downloaded on either "http" or "ftp", and you get continuous 'timeout' problems

1. If you're connecting to the internet **via a proxy**, download will fail. To solve the problem, identify the IP address and port of you proxy, and before running MODIStsp, run the following instructions: 

```
    library(httr)
    set_config(use_proxy(url="XXX.XXX.XXX.XXX", port=YYYY))
```

(substitute XXX.XXX.XXX.XXX and YYYY with the IP address and port number of the proxy, respectively)

2. If your connection is very slow, you could get frequent timeout problems while `MODIStsp` fetches the list of files available on NASA archives. This problem is more pronounced for ftp download, so switching to http may be a good idea !

#### - No luck ! Nothing worked 

If nothing of the above solves your problem, please report the issue here: https://github.com/ropensci/MODIStsp/issues

____________________________________________________________________________________


<a name="processing-errors"></a> 

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Processing Problems__ 

-   [`MODIStsp` fails while processing product **MxDyyy**](#product-error)
-   [I can't find the processed files](#whereare-error)
-   [Working with already downloaded hdf files](#old-hdfs)


<a name="product-error"></a> 

#### - MODIStsp fails while processing product **MxDyyy**

Although we tried to test functionality for all products, some bugs may be still present for specific products, due to the complexity and variability of MODIS hdfs structure. In that case, please report the issue here: https://github.com/ropensci/MODIStsp/issues 


<a name="whereare-error"></a>

#### - Where/How can I access the processed data ? 

Please see [here](output.html)

<a name="old-hdfs"></a> 

#### - Working with already downloaded hdf files

If you wish to use `MODIStsp` to process MODIS hdf images that you already downloaded from NASA servers, you should proceed like this: 

1. Place all the hdf files in a folder of your choice (e.g., "D:\\myfolder\\mydir\\hdf\_modis"). All the images must reside in the root of that folder (i.e., no subfolders)

2. Open MODIStsp GUI and set the processing parameters for your analysis as you would do if you had still to download the data. In particular, be sure to:

    - select the product _**corresponding to your hdf images**_ on the top of the GUI (e.g., MOD13Q1 v006), and set the processing layers you wish to analyze; 
    - set the _**processing period**_ to that of your available images; 
    - Set the _**spatial extent**_ to that of your available images. A minimum requirement is to set the Horizontal and Vertical "Required MODIS tiles" so to correspond to your imagery. If you also specify a bounding box for the output, that will be considered while creating the outputs. 
    
    
3. Set the "download method" to either: 

   - _**http/ftp**_: if you want `MODIStsp` to check if all the images for the selected product and time period are already available on your PC. `MODIStsp` will  take care of downloading any missing image (i.e., if you forgot to download something, or the processing period you provided is larger than that "covered" by the images you already have). 

   - _**offline**_: if you just want to process the images that you have. In that case, `MODIStsp` will not even connect to NASA servers, and just process all images that it finds in your "hdf folder" **_and which satisfy the processing parameters you specified_** (to be clear: if you have a folder with MOD13Q1 images, but the GUI is set as to process MOD09A1 data, nothing will happen. Analogously, if you have imagery from 2003 to 2009, but you set the processing period from 2010 to 2015 nothing will happen)
 
3. Set the "Output Folder for Original HDF files download" to the folder containing the hdf images, and the "Main Output Folder for Time Series storage to the folder where you want to store the results of the processing. 
4. Start the processing. 

____________________________________________________________________________________


<a name="other-errors"></a> 

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Other Problems__

Please report any other issues at https://github.com/ropensci/MODIStsp/issues
---
title: "Non-Interactive Execution from within R"
output: 
  github_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

`MODIStsp()` can be launched in non-interactive mode within an `R` session by setting the optional `GUI` parameter to FALSE, and either providing the desired processing argument in the call to the function, or providing a previously saved `opts_file` specifying the path to a JSON Options file previously saved through the GUI. 
This allows to exploit `{MODIStsp}` functionalities within generic `R` processing scripts.

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Specifying the processing parameters in the function call__

All processing parameters can be set in the call to `MODIStsp()`. Mandatory 
parameters are `selprod` (specifying the MODIS product), (one of) `bandsel`, `quality_bandsel` or `indexes_bandsel` (that specify the desired output layers), `out_folder`, `start_date` and `end_date`. 
`user` and `password` are also needed if `download_server` is not equal to `"offline"`.

The new function `MODIStsp_get_prodlayers()` allows easily retrieving the 
names of products and available layers based on product code, such as in: 

```{r echo=TRUE, message=FALSE, warning=FALSE}
library(MODIStsp) 
MODIStsp_get_prodlayers("M*D13Q1")
```

The other parameters are set automatically to default values (see [`MODIStsp()` documentation](../reference/MODIStsp.html) for details on the different available function arguments). 

For example, the following code processes layers __NDVI__ and __EVI__ and quality indicator __usefulness__ of product __M*D13Q1__, considering both Terra and Aqua platforms, for dates comprised between 2020-06-01 and 2020-06-15 and saves output to `R` `tempdir()`: 

```{r echo=TRUE, message=FALSE, warning=FALSE}
library(MODIStsp) 

# **NOTE** Output files of examples are saved to file.path by setting out_folder to "$tempdir".

# --> See name and available layers for product M*D13Q1
# MODIStsp_get_prodlayers("M*D13A2")

# --> Launch the processing
MODIStsp(gui             = FALSE, 
         out_folder      = "$tempdir", 
         selprod         = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
         bandsel         = c("EVI", "NDVI"), 
         quality_bandsel = "QA_usef", 
         indexes_bandsel = "SR", 
         user            = "mstp_test" ,
         password        = "MSTP_test_01",
         start_date      = "2020.06.01", 
         end_date        = "2020.06.15", 
         verbose         = FALSE)

# Outputs are in this case in subfolder "MODIStsp/VI_16Days_1Km_v6" of 
# `base::tempdir()`: 

out_fold <- file.path(tempdir(), "MODIStsp/VI_16Days_1Km_v6/") 
list.files(out_fold)
list.files(file.path(out_fold ,"EVI"))
list.files(file.path(out_fold ,"QA_usef"))
```

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Launching MODIStsp using a saved "Options file"__

Alternatively, you can run `MODIStsp()` without opening the GUI by specifying a previously saved options file: 

```{r echo=TRUE, message=FALSE, warning=FALSE}
library(MODIStsp) 

# **NOTE** Output files of examples are saved to file.path(tempdir(), "MODIStsp").

# --> Specify the path to a valid options file saved in advance from MODIStsp GUI 
# Here we use a test json file saved in MODIStsp installation folder which
# downloads and processed 3 MOD13A2 images over the Como Lake (Lombardy, Italy)
# and retrieves NDVI and EVI data, plus the Usefulness Index Quality Indicator.

opts_file <- system.file("testdata/test_MOD13A2.json", package = "MODIStsp")

# --> Launch the processing
MODIStsp(gui = FALSE, opts_file = opts_file, verbose = FALSE)

# Outputs are in this case in subfolder "MODIStsp/VI_16Days_1Km_v6" of 
# `base::tempdir()`: 

out_fold <- file.path(tempdir(), "MODIStsp/VI_16Days_1Km_v6") 
list.files(out_fold)
list.files(file.path(out_fold ,"EVI"))
```

## _Looping over different Options files_

If you need to process different MODIS products, you can prepare beforehand different
`MODIStsp()` options files by using the GUI, and then loop over them like this:

```{r echo=TRUE, message=FALSE, warning=FALSE}

opts_files <- c(system.file("testdata/test_MOD13A2.json", package = "MODIStsp"), 
                system.file("testdata/test_MOD10A2.json", package = "MODIStsp"))

for (opts_file in opts_files) {
  MODIStsp(gui       = FALSE, 
           opts_file = opts_file, 
           verbose   = FALSE)
}

# MOD13A2 ouputs
out_fold <- file.path(tempdir(), "MODIStsp/VI_16Days_1Km_v6") 
list.files(out_fold)
list.files(file.path(out_fold ,"NDVI"))

# MOD10A2 ouputs
out_fold <- file.path(tempdir(), "MODIStsp/Surf_Temp_8Days_1Km_v6") 
list.files(out_fold)
list.files(file.path(out_fold ,"LST_Night_1km"))
```

________________________________________________________________________________

<br>

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Specifying the processing parameters using a previously saved options file and overwriting some parameters__

Finally, it is possible to both specify a previously saved options file AND setting some parameters in the call to the function. This allows __easily performing similar processings, by only updating the required arguments__, as in the examples below. 

## __Looping over different spatial extents__

Specifying the `spafile` parameter while setting the `spatmeth` parameter to "file" overrides for example the output extent of the selected Options File. This allows performing the same preprocessing on different extents using a single Options File. For example:

```{r echo=TRUE, message=FALSE, warning=FALSE}
# Run the tool using the settings previously saved in a specific option file
# and specifying the extent from a spatial file allows to re-use the same
# processing settings to perform download and reprocessing on a different area
library(MODIStsp) 
opts_file    <- system.file("testdata/test_MOD13A2.json", package = "MODIStsp")
spatial_file <- system.file("testdata/lakeshapes/garda_lake.shp", package = "MODIStsp")
MODIStsp(gui       = FALSE, 
         opts_file = opts_file,
         spatmeth  = "file", 
         spafile   = spatial_file, 
         verbose   = FALSE)

# --> Create a character array containing a list of shapefiles (or other spatial files)
extent_list <- list.files(system.file("testdata/lakeshapes/", package = "MODIStsp"), full.names = TRUE, "\\.shp$")  

extent_list

# --> Loop on the list of spatial files and run MODIStsp using each of them to 
# automatically define the output extent (A separate output folder is created for
# each input spatial file).

for (single_shape in extent_list) {
  MODIStsp(gui       = FALSE, 
           opts_file = opts_file, 
           spatmeth  = "file", 
           spafile   = single_shape, 
           verbose   = FALSE)
}

# output files are placed in separate folders: 

outfiles_garda <- list.files(
  file.path(tempdir(), "MODIStsp/garda_lake/VI_16Days_1Km_v6/EVI"),
  full.names = TRUE)
outfiles_garda

library(raster)       
plot(raster(outfiles_garda[1]))

outfiles_iseo <- list.files(
  file.path(tempdir(), "MODIStsp/iseo_lake/VI_16Days_1Km_v6/EVI"),
  full.names = TRUE)
outfiles_iseo

plot(raster(outfiles_iseo[1]))

```
---
title: "List of supported MODIS products"
output: 
  github_document: default
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
webshot::install_phantomjs()
library(data.table)
library(dplyr)
library(knitr)
prodopts_file <- system.file("ExtData/", "MODIStsp_ProdOpts.RData", 
                             package = "MODIStsp")
load(prodopts_file)
out_list = list()
ind <- 1 
for (prod in seq_along(prod_opt_list)) {
  prod_data <- prod_opt_list[[prod]]
  for (ver in seq_along(prod_data)) {
    prod_verdata <- prod_data[[ver]]
    out_list[[ind]] <- data.frame(cat01 = prod_verdata[["cat01"]], 
                                  cat02 = prod_verdata[["cat02"]],
                                  Version   = prod_verdata[["v_number"]],
                                  Name = prod_verdata[["prod_fullname"]], 
                                  Resolution = prod_verdata[["native_res"]],
                                  Info  = prod_verdata[["www"]])
    ind <- ind + 1
  }
}
prod_df <- data.table::rbindlist(out_list)
```

# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Radiation Budget Variables__

```{r cat1, message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE}
cat_products <- subset(prod_df, cat01 == "Radiation Budget Variables")
```

<br>

## _Land Surface Reflectance_


```{r tab1, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Land Surface Reflectance") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution = as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>"))
library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________

<br>

##  _Snow Cover_

```{r tab2, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Snow Cover") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution =  as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>")) 

library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________

<br>

##  _Land Surface Temperature/Emissivity_

```{r tab3, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Land Surface Temperature/Emissivity") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>% 
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution =  as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>")) 

library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________

<br>

## _BRDF and Albedo_

```{r tab4, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "BRDF and Albedo") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>% 
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution =  as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>%
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>")) 

# dplyr::select(Name, Version, Resolution, Info)
# knitr::kable(cat2_products)
# tab <- knitr::kable(tab, digits = 1, align = c("l","l","c","c","l"),
#                     col.names = c("Product Code", "Product Name", "Version",
#                                   "Native Resolution", "Info"))
# print(tab, type = "html")
library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

## _Radiation_

```{r tabrad, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Radiation") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>% 
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution =  as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>%
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>")) 

# dplyr::select(Name, Version, Resolution, Info)
# knitr::kable(cat2_products)
# tab <- knitr::kable(tab, digits = 1, align = c("l","l","c","c","l"),
#                     col.names = c("Product Code", "Product Name", "Version",
#                                   "Native Resolution", "Info"))
# print(tab, type = "html")
library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________


# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Ecosystem Variables__

```{r cat2, message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE}
cat_products <- subset(prod_df, cat01 == "Ecosystem Variables")
```

<br>

## _Vegetation Indices_

```{r tab5, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Vegetation Indices") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution = as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>"))
library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________

<br>

## _LAI/FPAR_

```{r tab6, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "LAI/FPAR") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution = as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>"))

library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________

<br>

## _Evapotranspiration_

```{r tab7, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Evapotranspiration") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution = as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>")) 
library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________

<br>

## _Gross Primary Productivity_ 

```{r tab8, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Gross Primary Productivity") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution = as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>")) 

library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________

<br>

## _Net Primary Productivity_ 

```{r tab9, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Net Primary Productivity") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution = as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>")) 

library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________

<br>

## _Vegetation Continuous Cover/Fields_ 

```{r tab10, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Vegetation Continuous Cover/Fields") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution = as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>")) 

library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________


# <i class="fa fa-arrow-circle-o-right" aria-hidden="true"></i> __Land Cover Characteristics__

```{r cat3, message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE}
cat_products <- subset(prod_df, cat01 == "Land Cover Characteristics")
```

<br>

## _Thermal Anomalies and Fire_

```{r tab11, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Thermal Anomalies and Fire") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution = as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>"))
library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________

<br>

## _Land Cover_

```{r tab12, echo=FALSE, message=FALSE, warning=FALSE, paged.print=TRUE, results="asis"}
library(xtable)
tab <- cat_products %>% 
  dplyr::filter(cat02 == "Land Cover") %>% 
  dplyr::mutate(Code = as.factor(stringr::str_split_fixed(Name, ":", 2)[,1])) %>% 
  dplyr::mutate(Name = as.factor(stringr::str_split_fixed(Name, ":", 3)[,2])) %>%
  dplyr::select(Code, Name, Version, Resolution, Info) %>% 
  dplyr::mutate(Resolution = as.numeric(as.character(Resolution))) %>%  
  dplyr::mutate(Resolution = as.factor(format(Resolution, digits = 2, 
                                              nsmall = 1))) %>% 
  dplyr::mutate(Info = paste0("<a href='", Info,"'>", "Link","</a>"))
library(DT)
DT::datatable(tab, rownames = FALSE, filter = "top", style = "bootstrap",
              escape = -5, options = list(pageLength = 10, searching = TRUE))
```

________________________________________________________________________________
---
title: In memoriam
output: 
  github_document:
    toc_depth: 1
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

<div style='display: block; margin: auto;'>
<img src='https://s.gravatar.com/avatar/6f1e3990eec2eb3b4f8a1233b8098c4c?s=500'>
</div>

This package is dedicated to the memory of
<a href="https://lbusett.netlify.app/" target="_blank">Lorenzo Busetto</a>.

Lorenzo firmly believed in the open source model: he designed and built
this package in the framework of the activities conducted at the 
Italian National Research Council (CNR-IREA), but releasing it as a free software 
so that anyone can use it for free.
As an expert in using R for spatial analysis and remote sensing applications, 
he actively contributed to the geospatial community.

Lorenzo maintained {MODIStsp} until 21st October 2020, when he suddenly
passed away.
He left an enormous void among his collaborators and friends.

---
title: 'MODIStsp: A Tool for Automatic Preprocessing of MODIS Time Series - v2.0.5'
author: "Lorenzo Busetto, 
  Luigi Ranghetti ([ranghetti.l@irea.cnr.it](mailto:ranghetti.l@irea.cnr.it))"
bibliography: MODIStsp.bib
date: '`r Sys.Date()`'
output:
  rmarkdown::html_vignette:
    fig_caption: yes
    number_section: yes
    toc: no
    toc_depth: 2
  rmarkdown::pdf_document:
    fig_caption: yes
    number_sections: yes
    toc: yes
    toc_depth: 2
  rmarkdown::html_document:
    fig_caption: yes
    number_section: yes
    toc: no
    toc_depth: 1
urlcolor: blue
linkcolor: blue
vignette: >
  \usepackage[utf8]{inputenc}
  %\VignetteIndexEntry{MODIStsp: A Tool for Automatic Preprocessing of MODIS Time Series - v2.0.5}
  %\VignetteEngine{knitr::rmarkdown}
---

# Introduction

`{MODIStsp}` is an `R` package allowing to automatize the creation of time series
of rasters derived from MODIS Land Products data. It allows performing several 
preprocessing steps on MODIS data available within a given time period. 

Development of `{MODIStsp}` started from modifications of the `ModisDownload` `R` script
by Thomas Hengl [-@Hengl2010], and successive adaptations by Babak Naimi [-@Naimi2014]. 
The basic functionalities for download and preprocessing of MODIS datasets provided 
by these scripts were gradually incremented with the aim of: 

* developing a stand-alone application allowing to perform several preprocessing
steps (e.g., download, mosaicing, reprojection and resize) on all available MODIS 
land products by exploiting  a powerful and user-friendly GUI front-end;
* allowing the creation of time series of both MODIS original layers and additional 
Quality Indicators (e.g., data acquisition quality, cloud/snow presence, algorithm 
used for data production, etc. ) extracted from the aggregated bit-field QA layers;
* allowing the automatic calculation and creation of time series of several additional
Spectral Indexes starting form MODIS surface reflectance products.

All processing parameters can be easily set with a user-friendly GUI, although
non-interactive execution exploiting a previously created Options File is possible. 
Stand-alone execution outside an `R` environment is also possible, allowing to use
scheduled execution of MODIStsp to automatically update time series related to a 
MODIS product and extent whenever a new image is available. 

Required MODIS HDF files are automatically downloaded from NASA servers and resized,
reprojected, resampled and processed according to user's choices. For each desired 
output layer, outputs are saved as single-band rasters corresponding to each acquisition 
date available for the selected MODIS product within the specified time period. 
"R" _RasterStack_ objects with temporal information as well as Virtual raster files
(GDAL vrt and ENVI META files) facilitating access to the entire time series can
be also created.

# Installation

`{MODIStsp}` requires [R](https://cran.r-project.org) v >= 3.6.3.

## On Windows

You can install the stable version of `{MODIStsp}` from CRAN: 

`install.packages("MODIStsp")`

, or the development version (containing the latest improvements and bug fixes) 
from GitHub:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
install.packages("remotes")
library(remotes)
install_github("ropensci/MODIStsp")
```

    
## On Linux systems

To install `{MODIStsp}` on Linux, you need to be able to install the `{sf}` package, 
which requires several dependencies. See [here](https://github.com/r-spatial/sf#installing)
if you have trouble installing `{sf}`.

In addition, you need to install dependencies
required by the `{protolite}` package, required by `{geojson}`. See [here](https://github.com/jeroen/protolite/) for instructions
on installing them. 

Then, you can install the stable version of `{MODIStsp}` from CRAN:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
install.packages("MODIStsp")
```
, or the development version (containing the latest improvements and bug fixes) 
from GitHub:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
library(devtools)
install_github("ropensci/MODIStsp")
```

## On Mac OS

To install `{MODIStsp}` on MacOS, you need to be able to install the `{sf}` package, which requires several dependencies. See [HERE](https://github.com/r-spatial/sf#installing)
if you have trouble installing `{sf}`. 

Then, you can install the stable version of `{MODIStsp}` from CRAN:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
install.packages("MODIStsp")
```
, or the development version (containing the latest improvements and bug fixes) 
from GitHub:

```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE}
library(devtools)
install_github("ropensci/MODIStsp")
```

# Running the tool in Interactive Mode: the MODIStsp GUI

The easiest way to use `{MODIStsp}` is to use its powerful GUI (Graphical User 
Interface) for selection of processing options, and then run the processing. 

To open the GUI, load the package and launch the MODIStsp function, with no parameters:
```{r, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, caption=FALSE}
library(MODIStsp)
MODIStsp()
```
This **opens a GUI** from which processing options can be specified and eventually
saved (or loaded).

The  GUI allows selecting all processing options required for the creation of the
desired MODIS time series. The main available processing options are described in
detail in the following.

```{r GUIfig, echo=FALSE, fig.align="center", fig.width=10, message=FALSE, warning=FALSE}
  library(png)
  library(grid)
  library(knitr)
  img <- readPNG("GUI_1.PNG")
  grid.raster(img)
```

See (https://docs.ropensci.org/MODIStsp/articles/interactive_execution.html) for further
info and instructions regarding the usage of the GUI. 

# Non-Interactive Execution from within R

`MODIStsp()` can be launched in non-interactive mode within an `R` session by setting the optional `GUI` parameter to FALSE, and either providing the desired processing argument in the call to the function, or providing a previously saved `opts_file` specifying the path to a JSON Options file previously saved through the GUI. 
This allows exploiting `{MODIStsp}` functionalities within generic `R` processing scripts.

## __Specifying the processing parameters in the function call__

All processing parameters can be set in the call to `MODIStsp()`. Mandatory 
parameters are `selprod` (specifying the MODIS product), (one of) `bandsel`, `quality_bandsel` or `indexes_bandsel`, (that specify the desired output layers), `out_folder` and `start_date` and `end_date`. 
`user` and `password` are also needed if `download_server` is not equal to `"offline"`.

The function `MODIStsp_get_prodlayers()` allows to easily retrieve the 
names of products and available layers based on product code, such as in: 

```{r eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE}
library(MODIStsp) 
MODIStsp_get_prodlayers("M*D13Q1")
```

The other parameters are set automatically to default values (see [the function reference](https://docs.ropensci.org/MODIStsp/reference/MODIStsp_get_prodlayers.html) for details on the different available function arguments). 

For example, the following code processes layers __NDVI__ and __EVI__ and quality indicator __usefulness__ of product __M*D13Q1__, considering both Terra and Aqua platforms, for dates comprised between 2020-06-01 and 2020-06-15 and saves output to `R` `tempdir()`: 

```{r eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE}
library(MODIStsp) 

# **NOTE** Output files of examples are saved to file.path by setting out_folder to $tempdir.

# --> See name and available layers for product M*D13Q1
# MODIStsp_get_prodlayers("M*D13A2")

# --> Launch the processing
MODIStsp(
  gui             = FALSE, 
  out_folder      = "$tempdir", 
  selprod         = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
  bandsel         = c("EVI", "NDVI"), 
  quality_bandsel = "QA_usef", 
  indexes_bandsel = "SR", 
  user            = "mstp_test" ,
  password        = "MSTP_test_01",
  start_date      = "2020.06.01", 
  end_date        = "2020.06.15", 
  verbose         = FALSE,
  parallel        = FALSE
)

# Outputs are in this case in subfolder "MODIStsp/VI_16Days_1Km_v6" of 
# `base::tempdir()`: 

out_fold <- file.path(tempdir(), "MODIStsp/VI_16Days_1Km_v6/") 
list.files(out_fold)
list.files(file.path(out_fold ,"EVI"))
list.files(file.path(out_fold ,"QA_usef"))
```

_Note that this example, as well as the following ones, is run in single
core to follow CRAN policies, by setting `parallel = FALSE`;
users can exploit multicore functionalities skipping to set this argument._

## __Launching MODIStsp using a saved "Options file"__

Alternatively, you can run `MODIStsp()` without opening the GUI by specifying a previously saved options file: 

```{r eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE}
library(MODIStsp) 

# **NOTE** Output files of examples are saved to file.path(tempdir(), "MODIStsp").

# --> Specify the path to a valid options file saved in advance from MODIStsp GUI 
# Here we use a test json file saved in MODIStsp installation folder which
# downloads and processed 3 MOD13A2 images over the Como Lake (Lombardy, Italy)
# and retrieves NDVI and EVI data, plus the Usefulness Index Quality Indicator.

opts_file <- system.file("testdata/test_MOD13A2.json", package = "MODIStsp")

# --> Launch the processing
MODIStsp(gui = FALSE, opts_file = opts_file, verbose = FALSE, parallel = FALSE)

# Outputs are in this case in subfolder "MODIStsp/VI_16Days_1Km_v6" of 
# tempdir(): 

out_fold <- file.path(tempdir(), "MODIStsp/VI_16Days_1Km_v6") 
list.files(out_fold)
list.files(file.path(out_fold ,"EVI"))
```

### _Looping over different Options files_

If you need to process different MODIS products, you can prepare beforehand different
MODIStsp options files by using the GUI, and then loop over them like this:

```{r eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE}

opts_files <- c(system.file("testdata/test_MOD13A2.json", package = "MODIStsp"), 
                system.file("testdata/test_MOD10A2.json", package = "MODIStsp"))

for (opts_file in opts_files) {
  MODIStsp(gui = FALSE, opts_file = opts_file, verbose = FALSE, parallel = FALSE)
}

# MOD13A2 ouptuts
out_fold <- file.path(tempdir(), "MODIStsp/VI_16Days_1Km_v6") 
list.files(out_fold)
list.files(file.path(out_fold ,"NDVI"))

# MOD10A2 ouptuts
out_fold <- file.path(tempdir(), "MODIStsp/Surf_Temp_8Days_1Km_v6") 
list.files(out_fold)
list.files(file.path(out_fold ,"LST_Night_1km"))
```

## __Specifying the processing parameters using a previously saved options file and overwriting some parameters__

Finally, it is possible to both specify a previously saved options file AND setting some parameters in the call to the function. This allows __easily performing similar processings, by only updating the required arguments__, as in the examples below. 

## __Looping over different spatial extents__

Specifying the `spafile` parameter while setting the `spatmeth` parameter to `"file"` overrides for example the output extent of the selected Options File. This allows performing the same preprocessing on different extents using a single Options File. For example:

```{r eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE}
# Run the tool using the settings previously saved in a specific option file
# and specifying the extent from a spatial file allows to re-use the same
# processing settings to perform download and reprocessing on a different area
library(MODIStsp) 
opts_file    <- system.file("testdata/test_MOD13A2.json", package = "MODIStsp")
spatial_file <- system.file("testdata/lakeshapes/garda_lake.shp", package = "MODIStsp")
MODIStsp(
  gui = FALSE, 
  opts_file = opts_file,
  spatmeth = "file", 
  spafile = spatial_file, 
  verbose = FALSE, 
  parallel = FALSE
)

# --> Create a character array containing a list of shapefiles (or other spatial files)
extent_list <- list.files(system.file("testdata/lakeshapes/", package = "MODIStsp"), full.names = TRUE, "\\.shp$")  

extent_list

# --> Loop on the list of spatial files and run MODIStsp using each of them to 
# automatically define the output extent (A separate output folder is created for
# each input spatial file).

for (single_shape in extent_list) {
  MODIStsp(
    gui = FALSE, 
    opts_file = opts_file, 
    spatmeth = "file", 
    spafile = single_shape, 
    verbose = FALSE, 
    parallel = FALSE
  )
}

# output files are placed in separate folders: 

outfiles_garda <- list.files(
  file.path(tempdir(), "MODIStsp/garda_lake/VI_16Days_1Km_v6/EVI"),
  full.names = TRUE)
outfiles_garda

library(raster)       
plot(raster(outfiles_garda[1]))

outfiles_iseo <- list.files(
  file.path(tempdir(), "MODIStsp/iseo_lake/VI_16Days_1Km_v6/EVI"),
  full.names = TRUE)
outfiles_iseo

plot(raster(outfiles_iseo[1]))

```

# Scheduled Processing

Stand-alone non-interactive execution can be used to periodically and automatically 
update the time series of a selected product over a given study area. To do that, 
you should simply:

1.	Open the `{MODIStsp}` GUI, define the parameters of the processing specifying 
a date in the future as the "Ending Date" and save the processing options. Then
quit the program.
 
2. Schedule non-interactive execution of the launcher installed as seen before
(or located in the subdirectory `"MODIStsp/ExtData/Launcher"` of your library path)
as windows scheduled task (or linux "cron" job) according to a specified time 
schedule, specifying the path of a previously saved Options file as additional
argument.

#### On Linux

3. Edit your crontab by opening a terminal and type:

```bash
  crontab -e
```
 
4. Add an entry for the launcher. For example, if you have installed it in 
`/usr/bin` and you want to run the tool every day at 23.00, add the following row:
        
```bash
  0 23 * * * /bin/bash /usr/bin/MODIStsp -g -s "/yourpath/youroptions.json"
```
      
#### On Windows

3. Create a Task following <a href="https://docs.microsoft.com/en-us/previous-versions/windows/it-pro/windows-server-2008-R2-and-2008/cc748993(v=ws.11)?redirectedfrom=MSDN" target="_blank">these instructions</a>; add the path of the `MODIStsp.bat` launcher as Action (point 6), and specifying  `-g -s "X:/yourpath/youroptions.json"` as argument.

# Outputs Format and Naming Conventions

## Single-band outputs

Output raster files are saved in specific subfolders of the main output folder. 
In particular, **a separate subfolder** is created for each processed original 
MODIS layer, Quality Indicator or Spectral Index. Each subfolder contains one image
for each processed date, created according to the following naming conventions: 

```
myoutfolder/"Layer"/"ProdCode"_"Layer"_"YYYY"_"DOY"."ext"
```

<font size="2"> _(e.g., myoutfolder/NDVI/MOD13Q1\_NDVI\_2000\_065.dat)_ </font size="2">

, where: 

  - **_Layer_** is a short name describing the dataset (e.g., b1_Red, NDII, UI);
  - **_ProdCode_** is the code name of the MODIS product from which the image was
  derived (e.g., MOD13Q1);
  - **_YYYY_** and **_DOY_** correspond to the year and DOY (Day of the Year) of 
  acquisition of the original MODIS image;
  - **_ext_** is the file extension (.tif for GTiff outputs, or .dat for ENVI outputs). 

____________________________________________________________________________________

## Virtual multi-band outputs

ENVI and/or GDAL virtual time series files and _RasterStack_ RData objects are 
instead stored **in the "Time\_Series" subfolder** if required.

Naming convention for these files is as follow:

```
<path_of_out_folder>/Time_Series/<vrt_type>/<Sensor>/<Layer>/<ProdCode>_<Layer>_<StartDOY>_<StartYear>_<EndDOY>_<EndYear>_<suffix>.<ext> 
```
<font size="2"> _(e.g., `/my/out/folder/Time_Series/RData/Terra/NDVI/MOD13Q1_MYD13Q1_NDVI_49_2000_353_2015_RData.RData`)_ </font size="2"> 
             
, where: 

  - `<vrt_type>` indicates the type of virtual file (`"RData"`, `"GDAL"` or `"ENVI_META"`);
  - `<Sensor>` indicates to which MODIS sensor the time series belongs (`"Terra"`,
  `"Aqua"`, `"Mixed"` or `"Combined"` (for MCD* products));
  - `<Layer>` is a short name describing the dataset (e.g., `"b1_Red"`, `"NDII"`, `"UI"`);
  - `<ProdCode>` is the code name of the MODIS product from which the image was 
  derived (e.g., `"MOD13Q1"`);
  - `<StartDOY>`, `<StartYear>`, `<EndDOY>` and `<EndYear>` indicate the 
  temporal extent of the time serie created; 
  - `<suffix>` indicates the type of virtual file (`"ENVI"`, `"GDAL"` or `"RData"`);
  - `<ext>` is the file extension (`".vrt"` for GDAL virtual files, `"META"` for 
  ENVI meta files or `"RData"` for `R` raster stacks). 


# Accessing the processed time series from R

Preprocessed MODIS data can be retrieved within `R` either by accessing the single-date raster files, or by loading the saved _RasterStack_ objects. 

Any single-date image can be accessed by simply opening it with a `raster` command: 

```{R eval=FALSE, tidy=TRUE, highlight=TRUE}
library(raster)
modistsp_file <- "/my_outfolder/EVI/MOD13Q1_2005_137_EVI.tif"
my_raster     <- raster(modistsp_file)
```

`rasterStack` time series containing all the processed data for a given parameter (saved in the `"Time Series/RData"` subfolder - see [here](https://docs.ropensci.org/MODIStsp/articles/output.html) for details) can be opened by: 

```{r eval=FALSE}
in_virtual_file <- "/my_outfolder/Time_Series/RData/Terra/EVI/MOD13Q1_MYD13Q1_NDVI_49_2000_353_2015_RData.RData" 
indata          <- get(load(in_virtual_file))
```

This second option allows accessing the complete data stack and analyzing it using the functionalities for raster/raster time series analysis, extraction and plotting provided for example by the `{raster}` or   `{rasterVis}` packages.


## Extracting Time Series Data on Areas of Interest

`{MODIStsp}` provides an efficient function (`MODIStsp_extract()`) for extracting time series data at specific locations. The function takes as input a _RasterStack_ virtual object created by `MODIStsp()` (see above), the starting and ending dates for the extraction and a standard _Sp*_ object (or an ESRI shapefile name) specifying the locations (points, lines or polygons) of interest, and provides as output a `xts` object or `data.frame` containing time series data for those locations. 

If the input is of class _SpatialPoints_, the output object contains one column for each point specified, and one row for each date. If it is of class _SpatialPolygons_ (or _SpatialLines_), it contains one column for each polygon (or each line), with values obtained applying the function specified as the `FUN` argument (e.g., mean, standard deviation, etc.) on pixels belonging to the polygon (or touched by the line), and one row for each date. 

As an example the following code:

```{r, eval=FALSE}
  #Set the input paths to raster and shape file
  infile    <- 'myoutfolder/Time_Series/RData/Mixed/MOD13Q1_MYD13Q1_NDVI_49_2000_353_2015_RData.RData'  
  shpname   <- 'path_to_file/rois.shp'  
  #Set the start/end dates for extraction
  startdate <- as.Date("2010-01-01")  
  enddate   <- as.Date("2014-12-31")    
  #Load the RasterStack
  inrts     <- get(load(infile)) 
  # Compute average and St.dev
  dataavg   <- MODIStsp_extract(inrts, shpname, startdate, enddate, FUN = 'mean', na.rm = T)
  datasd    <- MODIStsp_extract (inrts, shpname, startdate, enddate, FUN = 'sd', na.rm = T)
  # Plot average time series for the polygons
  plot.xts(dataavg) 
```

  loads a _RasterStack_ object containing 8-days 250 m resolution time series for the 2000-2015 period and extracts time series of average and standard deviation values over the different polygons of a user's selected shapefile on the 2010-2014 period.
  
# Problems and Issues

Solutions to some common **installation and processing problems** can be found in `{MODIStsp}` FAQ: 
https://docs.ropensci.org/MODIStsp/articles/faq.html

Please **report any issues** you may encounter in our issues page on GitHub:
https://github.com/ropensci/MODIStsp/issues

  
# Citation
  
To cite MODIStsp please use:

L. Busetto, L. Ranghetti (2016) MODIStsp: An R package for automatic preprocessing of MODIS
  Land Products time series, Computers & Geosciences, Volume 97, Pages
  40-48, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2016.08.020, URL: https://github.com/ropensci/MODIStsp. 
  


# References
---
title: "Bug discovered in spectral indices computation"
output: 
  github_document: default
---

We are sorry to report that we recently discovered a nasty bug (or rather, a stupid mistake...) in the [MODIStsp](https://github.com/ropensci/MODIStsp) package.

The bug led to improper computation of custom spectral indices in the case that their formula included addition or subtraction operations on reflectance values (e.g., something like $\frac{(\rho_{NIR}+0.1)}{\rho_{Red}}$, with $\rho$ indicating a reflectance).

## What is affected

* Values of the following *Additional Spectral Indices* selectable using the MODIStsp GUI:

  - EVI
  - SAVI
  
  , *in the case that the **Apply Scale/Offset** option was set to "No"*
  
* Values of any *custom spectral indexes* added by the user, in case they included additive or subtractive coefficients. 

  , *in the case that the **Apply Scale/Offset** option was set to "No"*

## What is NOT affected  

* Values of spectral indexes available in MODIS HDF images as original sds layers (e.g., EVI in MOD13Q1)
  
* Values of any additional / custom spectral indexes in case they did not include additive or 
subtractive coefficients, or the **Apply Scale/Offset** option was set to "Yes"


## What to do if you are affected 

The bug is now fixed on the GitHub version. A patched release will be made available on CRAN as soon as possible. 
Unfortunately, if you have time series processed with the old version falling in the "What is affected" category, there's nothing you can do, save for reprocessing them. 

**We are truly sorry for the problem**, which somehow managed to slip under the radar until now. 
We hope it will not bring you too much trouble!

## What exactly was the problem? 

This is **so basic that can easily go unnoticed.** So it's better to document it...

MODIS reflectances are stored in HDF layers as integers with a 10000 scale factor (e.g., a 0.1 reflectance is stored as 1000). If you need to "manually" compute and index such as SAVI: 

$SAVI = \frac{(\rho_{NIR} - \rho_{Red})}{(\rho_{NIR} + \rho_{Red} + 0.5)} * (1 + 0.5)$).

starting from MODIS reflectances, you must take care of multiplying the MODIS data by 10E-4 beforehand. Your formula then becomes: 

$SAVI = \frac{(0.0001 * b2_{NIR} - 0.0001 * b1_{Red})}{0.0001 * b2_{NIR}  + 0.0001 * b1_{Red} + 0.5} * (1 + 0.5)$).

, otherwise the additive constants (in this case, the $+ 0.5$ in the denominator) would be made practically irrelevant.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_GUI.R
\name{MODIStsp_GUI}
\alias{MODIStsp_GUI}
\title{Build and manage the MODIStsp GUI}
\usage{
MODIStsp_GUI()
}
\value{
the function is called for its side effects - opening
the GUI and allowing to set, save, load options and eventually
launch the processing.
}
\description{
Function used to generate and handle the GUI used to allow selection of
MODIStsp processing parameters.
}
\note{
License: GPL 3.0
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_nodata_values.R
\name{split_nodata_values}
\alias{split_nodata_values}
\alias{create_nodata_rcl}
\title{Split NODATA values or create matrix for reclassification}
\usage{
split_nodata_values(nodata_in, take_all = TRUE)

create_nodata_rcl(nodata_in, nodata_out)
}
\arguments{
\item{nodata_in}{Character vector corresponding to input NoData values
as saved in the xml product file (one or more values per band).}

\item{take_all}{Logical: if TRUE (default), all the NoData values are considered;
if FALSE, only the last one is taken.
See "details" for the meaning of this parameter.}

\item{nodata_out}{Character vector corresponding to output NoData values
as saved in the xml product file (one single value per band).}
}
\value{
\link{split_nodata_values}  returns a list with the same length
of \code{nodata_in} vector, in which each element
is a vector with all the NoData values.

\link{create_nodata_rcl} returns a list of matrices in the format
specified for parameter \code{rcl} in \link[raster:reclassify]{raster::reclassify}.
The parameter \code{right} is intended to be used as \code{right = NA}.
}
\description{
Internal functions:
\link{split_nodata_values} splits the ranges of NODATA saved
in the xml product file to a readable vector of NoData values;
\link{create_nodata_rcl} creates the matrix for the reclassification of NODATA
values to be used with \link[raster:reclassify]{raster::reclassify} function.
}
\details{
MODIS products can have more than one NoData values (sometimes
with different meanings, e.g. 255 = "fill" and 254 = "detector saturated"
in \href{https://nsidc.org/data/mod10a1}{MOD09A1} product). By setting
"Change NoData values" to "Yes" in the GUI, all the NoData values are
coerced to one single new NoData value; conversely, setting it to "No"
only one value is assumed to be NoData.
The parameter \code{take_all} is assumed to be used in this way, by using this
function with \code{take_all = TRUE} with "Change NoData values" = "Yes" and
\code{take_all = FALSE} with "Change NoData values" = "No".

In the xml product file, NoData ranges are set as:
\itemize{
\item \code{x} for products with single NoData values;
\item \verb{x,y,z} for products with a vector of NoData values;
\item \code{x:y} for products with a range of NoData values;
\item \verb{x:y,z} for a combination of NoData ranges and/or values.
}

In \link{split_nodata_values} \emph{NoData values are assumed to be integer}:
this means that intervals are splitted in integer values
(e.g. "250:255" becomes "250 251 252 253 254 255").
Conversely, function \link{create_nodata_rcl} creates intervals, so it can also
manage float values (in practice, this should not make difference within
MODIS products, since NoData values are always integer values).

This function interprets these strings and convert them in vectors with
single values.
Notice that the first NoData value is the only one which is considered if
'Change NoData values' was set to 'No'.
}
\examples{
MODIStsp:::create_nodata_rcl(c("255","250,254:255"), c("255","255"))
}
\author{
Luigi Ranghetti, phD (2018) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_resetindexes.R
\name{MODIStsp_resetindexes}
\alias{MODIStsp_resetindexes}
\title{Remove custom spectral indexes}
\usage{
MODIStsp_resetindexes()
}
\value{
The function is called for its side effects. On success, the
MODIStsp_indexes.json file is modified so to
remove all previously custom-specified Spectral Indexes.
}
\description{
Function used to remove all user-defined Spectral Indexes from
MODIStsp, thus resetting the list of available indexes
to the default ones.
}
\note{
License: GPL 3.0
}
\examples{
\dontrun{
# Remove all custom-defined spectral indexes from an options file

# Add a custom index for testing purposes
library(jsonlite)
opts_jsfile = system.file("testdata/test_addindex.json",
                             package = "MODIStsp")
 MODIStsp_addindex(
   opts_jsfile = opts_jsfile,
   gui = FALSE,
   new_indexbandname = paste0("Index_", as.character(sample(10000, 1))),
   new_indexformula = "b1_Red - b2_NIR",
   new_indexfullname = paste0("Index_", as.character(sample(10000, 1)))
   )

 opts <- jsonlite::fromJSON(indexes_file)
 opts$custom_indexes[1]

 # Now remove all custom indexes
 MODIStsp_resetindexes()
 opts <- jsonlite::fromJSON(opts_jsfile)
 opts$custom_indexes[1]
 }

}
\seealso{
\link{MODIStsp_addindex}
}
\author{
Lorenzo Busetto, phD (2014-2017) \email{busetto.l@irea.cnr.it}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_files_existance.R
\name{check_files_existence}
\alias{check_files_existence}
\title{Check if all files required for a given date already exist}
\usage{
check_files_existence(
  out_prod_folder,
  file_prefix,
  yy,
  DOY,
  bandnames,
  bandsel_orig_choice,
  indexes_bandnames,
  indexes_bandsel,
  quality_bandnames,
  quality_bandsel,
  out_format
)
}
\arguments{
\item{out_prod_folder}{\code{character} MODIStsp output folder}

\item{file_prefix}{\code{character} File prefix of the processed product
(e.g., MOD13Q1)}

\item{yy}{\code{character} year}

\item{DOY}{\code{character} doy}

\item{bandnames}{\verb{character array} Bandnames of the MODIS product}

\item{bandsel_orig_choice}{\verb{numeric 0/1 array} Indicates which original MODIS
layers were selected for processing (does not contain names of bands needed
to compute SIs but not selected by the user!)}

\item{indexes_bandnames}{\verb{character array} Names of available spectral
indexes (standard + custom) available for the currently processed product}

\item{indexes_bandsel}{\verb{numeric 0/1 array} Indicates which spectral indexes
were selected for processing}

\item{quality_bandnames}{\verb{character array} Name of available Quality
Indicators for the currently processed product}

\item{quality_bandsel}{\verb{numeric 0/1 array} Indicates which Quality Indicators
were selected}

\item{out_format}{\code{character} GTiff or ENVI}
}
\value{
check - logical = 1 if all expected output files are already existing
}
\description{
Accessory function used to see if all expected out files for the
selected date are already present in the output folder. If all expected out
files are already present, check_files is set to TRUE, and the date is
skipped in MODIStsp_process.
}
\note{
License: GPL 3.0
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_proc_opts.R
\name{check_proc_opts}
\alias{check_proc_opts}
\title{check_proc_opts}
\usage{
check_proc_opts(proc_opts)
}
\arguments{
\item{proc_opts}{data frame of parameters passed by \code{MODIStsp}}
}
\value{
NULL - processing interrupted if any condition is not met
}
\description{
helper function used to check consistency of processing
parameters
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_yeardates.R
\name{get_yeardates}
\alias{get_yeardates}
\title{identify dates to be processed for a year}
\usage{
get_yeardates(download_range, yy, start_year, end_year, start_date, end_date)
}
\arguments{
\item{download_range}{\code{character ["Full" | "Seasonal"]} If "full", all the
available images between the starting and the ending dates are downloaded;
If "seasonal", only the images included in the season are downloaded
(e.g: if the starting date is 2005-12-01 and the ending is 2010-02-31, only
the images of December, January and February from 2005 to 2010 - excluding
2005-01, 2005-02 and 2010-12 - are downloaded), Default: Full}

\item{yy}{\code{numeric} year for which the processing dates need to be identified}

\item{start_year}{\code{numeric} start year of current \code{MODIStsp_process} run}

\item{end_year}{\code{numeric} end year of current \code{MODIStsp_process} run.}

\item{start_date}{\verb{character`` Start date for images download and preprocessing (yyyy.mm.dd) of current }MODIStsp_process` run.}

\item{end_date}{\code{character} Start date for images download and preprocessing
(yyyy.mm.dd) of current \code{MODIStsp_process} run.}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
helper function needed to identify the ranges of dates to
be processed for a given year as a function of \code{download_range} selection
and starting/ending dates and years
}
\author{
Lorenzo Busetto, phD (2017)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_process_indexes.R
\name{MODIStsp_process_indexes}
\alias{MODIStsp_process_indexes}
\title{MODIStsp helper for computing spectral indexes}
\usage{
MODIStsp_process_indexes(
  out_filename,
  out_prod_folder,
  formula,
  bandnames,
  nodata_out,
  indexes_nodata_out,
  file_prefix,
  compress,
  yy,
  out_format,
  DOY,
  scale_val
)
}
\arguments{
\item{out_filename}{\code{character} basename of the file in to which save results}

\item{out_prod_folder}{\code{character} output folder for the product used to retrieve filenames
of rasters of original bands to be used in computations}

\item{formula}{\code{character} Index formula, as derived from XML file and stored in prod_opts
within previous_file}

\item{bandnames}{\code{character} array of names of original HDF layer. Used to identify the
bands required for index computation}

\item{nodata_out}{\code{character} array of NoData values of reflectance bands}

\item{indexes_nodata_out}{\code{character} NoData value for resulting raster}

\item{file_prefix}{\code{character} used to retrieve filenames of rasters of original bands
to be used in computations}

\item{compress}{\code{character} compression option for GTiff files}

\item{yy}{\code{character} year string used to retrieve filenames of rasters of original bands
to be used in computations}

\item{out_format}{\code{character} string used to retrieve filenames of rasters of original bands
to be used in computations}

\item{DOY}{\code{character} doy string used to retrieve filenames of rasters of original bands to be
used in computations}

\item{scale_val}{\code{character} (Yes/No) if Yes, output values in are computed as float -1 - 1,
otherwise integer -10000 - 10000}
}
\value{
NULL - new raster file saved in out_filename
}
\description{
function used to compute spectral indexes, given the index formula
}
\details{
the function parses the index formula to identify the required bands. On the basis
of identified bands, it retrieves the reflectance bands required, gets the data into R raster
objects, performs the computation and stores results in a GeoTiff or ENVI raster file
}
\note{
License: GPL 3.0
}
\author{
Lorenzo Busetto, phD (2017)

Luigi Ranghetti, phD (2017) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mod_dirs.R
\name{get_mod_dirs}
\alias{get_mod_dirs}
\title{Get list of MODIS data folders from http server}
\usage{
get_mod_dirs(
  http,
  download_server,
  user,
  password,
  yy,
  n_retries,
  gui,
  out_folder_mod
)
}
\arguments{
\item{http}{\code{character} http site on lpdaac corresponding to the selected MODIS
product}

\item{download_server}{\code{character ["http" | "offline"]} download service
to be used; if NA, the script tries to download with http.}

\item{user}{\code{character} username for earthdata http server}

\item{password}{\code{character} password for earthdata http server}

\item{yy}{\code{character} Year for which the folder containing HDF images are to
be identified}

\item{n_retries}{\code{numeric} number of times the access to the http server
should be retried in case of error before quitting, Default: 20}

\item{gui}{`logical`` indicates if processing was called from the GUI
environment or not. If not, processing messages are sent to a log file
instead than to the console/GTK progress windows.}

\item{out_folder_mod}{\code{character} output folder for MODIS HDF storage}
}
\value{
\verb{character array} listing all available folders (a.k.a. dates) for
the requested MODIS product on lpdaac http archive, for the years
included in the time range selected for processing.
}
\description{
Accessory function to get the full list of directories on the
lpdaac http site containing data included in the time range selected for
processing (modified after Barry Rowlingson function):
}
\note{
License: GPL 3.0
}
\author{
Original code by Babak Naimi (\code{.getModisList}, in
\href{http://r-gis.net/?q=ModisDownload}{ModisDownload.R})
modified to adapt it to MODIStsp scheme and to http archive (instead than old
FTP) by:

Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2016-2017) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_install_launcher.R
\name{install_MODIStsp_launcher}
\alias{install_MODIStsp_launcher}
\title{Install a launcher for MODIStsp}
\usage{
install_MODIStsp_launcher(
  bin_dir = NA,
  rscript_dir = NA,
  desktop_dir = NA,
  desktop_shortcut = TRUE,
  sudo = FALSE
)
}
\arguments{
\item{bin_dir}{\itemize{
\item on Linux, directory in which the link to the bash script should be
placed, Default: "/usr/bin" - use of a path included in the PATH environment variable is
suggested;
\item on Windows, directory where to place the menu entry in the Start Menu,
Default: Start Menu -> Programs -> MODIStsp.
}}

\item{rscript_dir}{\code{character} in Windows only, the path of the directory in which
Rscript is installed (usually is \"C:/Progra~1/R/R-\code{version}/bin/\code{x64}").
Edit this parameter if R is installed in a custom directory.}

\item{desktop_dir}{\code{character}
\itemize{
\item on Linux, directory in which the desktop entry should be placed, Default: /usr/share/applications;
\item on Windows, directory where to place the desktop entry, Default: "Desktop"
(Ignored if desktop_shortcut = FALSE).
}}

\item{desktop_shortcut}{\code{logical} indicates if the desktop entry or the
desktop shortcut should be created, Default: TRUE.}

\item{sudo}{(Linux only) \code{logical}  indicates if administrator rights have to
be used to write within bin_dir and desktop_dir, If FALSE the root password is requested
when launching the function. Note that using default values of bin_dir and desktop_dir
requires to set this option to TRUE (or to launch the script in a root session of R),
Default: FALSE}
}
\value{
The function is called for its side effects.
}
\description{
Function which allows to use MODIStsp in batch mode by creating links
}
\details{
MODIStsp can be used also as a stand-alone tool (i.e., without opening RStudio
or R-GUI) by launching a bash/batch script, which is stored in the installation
folder (/ExtData/Launcher)
To allow to easily find it, this function creates a desktop entry and a symbolic link
to the bash script (on Linux) or a link in the Start Menu to the batch script and a
shortcut on the desktop (on Windows).
\strong{Note that}, if the packages MODIStsp is installed in a version-dependent directory
(as the default one is), this function should be re-executed after an R upgrade,
otherwise the links would continue to point to the old package version!
}
\note{
License: GPL 3.0
}
\examples{
# Linux: common installation (script in /usr/bin,
# desktop entry in /usr/share/applications)
# (requires administrator permissions)
\dontrun{
# the administrator password is asked interactively
install_MODIStsp_launcher(sudo = TRUE)
}

# Linux: installation in a directory which does not require administrator
# permissions
\dontrun{
install_MODIStsp_launcher(bin_dir = "~/bin", desktop_dir = "~/Desktop")
}

# Windows: common installation
# (script in the Start Menu and shortcut on the desktop)
\dontrun{
install_MODIStsp_launcher()
}
}
\author{
Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_reqbands.R
\name{get_reqbands}
\alias{get_reqbands}
\title{Identify the MODIS original bands needed for a given processing run}
\usage{
get_reqbands(
  bands_indexes_matrix,
  indexes_bandsel,
  indexes_bandnames,
  quality_bandsel,
  quality_bandnames,
  out_prod_folder,
  file_prefix,
  yy,
  DOY,
  out_format,
  reprocess
)
}
\arguments{
\item{bands_indexes_matrix}{\code{matrix} built by \code{set_bandind_matrix}}

\item{indexes_bandsel}{\verb{character array}Spectral Indexes to be computed starting from reflectance bands.
You can get a list of available quality layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$indexes_bandnames),
Default: NULL}

\item{indexes_bandnames}{names of all indexes available for the product being processed}

\item{quality_bandsel}{\verb{character array} Quality Indicators to be computed starting from
bit fields of original MODIS layers. You can get a list of available quality layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$quality_bandnames),
Default: NULL}

\item{quality_bandnames}{names of all quality indicators available for the product being processed}

\item{out_prod_folder}{\code{character} Main folder where the MODIStsp processed
raster will be stored. Used to check if a given processed image already exists.}

\item{file_prefix}{File prefix corresponding to the MODIS product being
processed. Used to check if a given processed image already exists.}

\item{yy}{Year corresponding to the image being processed. Used to check if
a given processed image already exists.}

\item{DOY}{DOY corresponding to the image being processed. Used to check if
a given processed image already exists.. Used to check if a given processed image already exists.}

\item{out_format}{\code{character ["ENVI" | "GTiff"]} Desired output format.}

\item{reprocess}{\code{logical} If TRUE, reprocess data for already existing dates.}
}
\value{
req_bands_indexes
}
\description{
Helper function used in MODIStsp_process to identify which
MODIS hdf layers are required for the current process. The required layers
include all MODIS original layers selected by the user, plus all those
required to compute the Spectral Indexes and Quality Indicators selected
by the user
}
\author{
Lorenzo Busetto, phD (2017)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_bandind_matrix.R
\name{set_bandind_matrix}
\alias{set_bandind_matrix}
\title{Helper function to determine the bands needed to compute SIs and QIs}
\usage{
set_bandind_matrix(
  bandnames,
  bandsel,
  indexes_bandnames,
  indexes_bandsel,
  indexes_formula,
  quality_bandnames,
  quality_bandsel,
  quality_source
)
}
\arguments{
\item{bandnames}{names of all layers available for the product being processed}

\item{bandsel}{\verb{character array} Original MODIS layers to be processed.
You can get a list of available layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$bandnames),
Default: NULL}

\item{indexes_bandnames}{names of all indexes available for the product being processed}

\item{indexes_bandsel}{\verb{character array}Spectral Indexes to be computed starting from reflectance bands.
You can get a list of available quality layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$indexes_bandnames),
Default: NULL}

\item{indexes_formula}{formulas of all indexes available for the product being processed}

\item{quality_bandnames}{names of all quality indicators available for the product being processed}

\item{quality_bandsel}{\verb{character array} Quality Indicators to be computed starting from
bit fields of original MODIS layers. You can get a list of available quality layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$quality_bandnames),
Default: NULL}

\item{quality_source}{sources of data (original layers) of all quality indicators
available for the product being processed}
}
\value{
\code{matrix} containing info on which bands are needed for computing
each available QI or SI
}
\description{
FUNCTION_DESCRIPTION
}
\author{
Lorenzo Busetto, phD (2017)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reproj_bbox.R
\name{reproj_bbox}
\alias{reproj_bbox}
\title{Reproject a bounding box}
\usage{
reproj_bbox(bbox, in_proj, out_proj, enlarge = TRUE)
}
\arguments{
\item{bbox}{The input bounding box (it can be a matrix obtained from \code{sp::bbox()},
or a numeric vector in the format (xmin, ymin, xmax, ymax)).}

\item{in_proj}{(\code{crs} | \code{character}) crs of the input projection,
or string coercible to it using \code{sf::st_crs()} (e.g., WKT or numeric
EPSG code)}

\item{out_proj}{\code{crs} \code{crs} of the output projection, or string coercible to
it using \code{sf::st_crs()} (e.g., WKT or numeric EPSG code)}

\item{enlarge}{`logical`` if TRUE, the reprojected bounding box is the
one which completely include the original one; if FALSE, it is simply the
one obtained by reprojecting the upper-left and the lower-right corners.}
}
\description{
Helper function used to reproject bounding boxes; setting the parameter
'enlarge' allows to choose if the new one would be the one which completely
includes the original extent in the output projection, or if is simply the
one obtained by reprojecting the upper-left and the lower-right corners.
}
\note{
License: GPL 3.0
}
\author{
Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_download.R
\name{MODIStsp_download}
\alias{MODIStsp_download}
\title{MODIStsp download function}
\usage{
MODIStsp_download(
  modislist,
  out_folder_mod,
  download_server,
  http,
  n_retries,
  use_aria,
  date_dir,
  year,
  DOY,
  user,
  password,
  sens_sel,
  date_name,
  gui,
  verbose
)
}
\arguments{
\item{modislist}{\verb{character array} List of MODIS images to be downloaded for
the selected date (as returned from \code{get_mod_filenames}). Can be a single
image, or a list of images in case different tiles are needed!}

\item{out_folder_mod}{\code{character} Folder where the hdfs are to be stored}

\item{download_server}{\code{character ["http"]} Server to be used.}

\item{http}{\code{character} Address of the http server for the selected product.}

\item{n_retries}{\code{numeric} Max number of retry attempts on download. If
download fails more that n_retries times consecutively, abort}

\item{use_aria}{\code{logical} if TRUE, download using aria2c}

\item{date_dir}{\verb{character array} Sub-folder where the different images
can be found (element of the list returned from \code{get_mod_dirs}, used in case
of http download to generate the download addresses).}

\item{year}{\code{character} Acquisition year of the images to be downloaded}

\item{DOY}{\verb{character array} Acquisition doys of the images to be downloaded}

\item{user}{\code{character} Username for http download}

\item{password}{\code{character} Password for http download}

\item{sens_sel}{\code{character ["terra" | "aqua"]} Selected sensor.}

\item{date_name}{\code{character} Date of acquisition of the images to be downloaded.}

\item{gui}{\code{logical} Indicates if on an interactive or non-interactive execution
(only influences where the log messages are sent).}

\item{verbose}{\code{logical} If FALSE, suppress processing messages, Default: TRUE}
}
\value{
The function is called for its side effects
}
\description{
Internal function dealing with download of MODIS hdfs from
http remote server for a given date.
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_process_bands.R
\name{MODIStsp_process_bands}
\alias{MODIStsp_process_bands}
\title{MODIStsp helper for processing original HDF layers}
\usage{
MODIStsp_process_bands(
  out_folder_mod,
  modislist,
  outproj_str,
  mod_proj_str,
  sens_sel,
  band,
  bandname,
  date_name,
  datatype,
  nodata_in,
  nodata_out,
  full_ext,
  bbox,
  scale_val,
  scale_factor,
  offset,
  out_format,
  outrep_file,
  compress,
  out_res_sel,
  out_res,
  resampling,
  nodata_change,
  gui,
  verbose,
  parallel
)
}
\arguments{
\item{out_folder_mod}{\code{character} Output folder for original HDF storage.
If \code{"$tempdir"} (default), a temporary directory is used.}

\item{modislist}{\verb{character array} List of MODIS images to be downloaded for
the selected date (as returned from \code{get_mod_filenames}). Can be a single
image, or a list of images in case different tiles are needed!}

\item{outproj_str}{\code{character} EPSG or WKT of output projection.}

\item{mod_proj_str}{\code{character} EPSG or WKT of MODIS projection.}

\item{sens_sel}{\code{character ["terra" | "aqua"]} Selected sensor.}

\item{band}{\code{numeric} band number corresponding to the HDF layer to be
processed}

\item{bandname}{\code{character} Name of the HDF layer to be processed.}

\item{date_name}{\code{character} Date of acquisition of the images to be
downloaded.}

\item{datatype}{\code{character} Datatype to the HDF layer to be processed.}

\item{nodata_in}{\code{numeric} Original nodata value to the HDF layer to be
processed.}

\item{nodata_out}{\code{numeric} Output nodata value to the HDF layer to be
processed.}

\item{full_ext}{\code{logical} If TRUE, process full tiles, if FALSE, process
bbox}

\item{bbox}{\code{numeric(4)} Output bounding box (xmin, ymin, xmax, ymax) in
out_proj coordinate system. Ignored if spatmeth == "tiles", Default: NULL}

\item{scale_val}{\code{logical} If TRUE,  scale and offset are applied to
original MODIS layers, and Spectral Indexes are saved as floating point. If
FALSE, no rescaling is done and Spectral Indexes are saved as integer, with a
10000 scaling factor.}

\item{scale_factor}{\code{numeric} Scale factor to be applied to the HDF layer
to be processed (Ignored if \code{scale_val} == FALSE).}

\item{offset}{\code{numeric} Offset to be applied to the HDF layer
to be processed (Ignored if \code{scale_val} == FALSE).}

\item{out_format}{\code{character ["ENVI" | "GTiff"]} Desired output format.}

\item{outrep_file}{\code{character} Full path of the file where results of the
processing are to be stored (created in \code{MODIStsp_process})}

\item{compress}{\code{character ["None" | "PACKBITS" | "LZW" | "DEFLATE"]}
Compression method for GTiff outputs (Ignored if \code{out_format == ENVI})}

\item{out_res_sel}{\verb{character ["Native", "User Defined}]. If "Native", the
outputs keep the original resolution of MODIS HDF images. Otherwise, the value
set in "out_res" is used.}

\item{out_res}{\code{float} Output resolution (in output projection measurement
unit). Ignored if out_res_sel == "Native".}

\item{resampling}{\verb{character ["near" | "bilinear" | "cubic" | "cubicspline", |lanczos"|, "average"|, "mode", |"max"|, |"min"|, |"q1"|, |"q3"|, |"sum"|]}
Resampling method to be used by \code{gdalwarp}.}

\item{nodata_change}{\code{logical} if TRUE, NoData values are set to the max value
of the datatype of the layer on the MODIStsp output rasters. NOTE: If multiple
nodata values are reported for a layer, all are reset to the new value.}

\item{gui}{\code{logical} if TRUE: the GUI is opened before processing. If FALSE:
processing parameters are retrieved from the provided \code{opts_file}
argument), Default: TRUE}

\item{verbose}{\code{logical} If FALSE, suppress processing messages, Default: TRUE}

\item{parallel}{\code{logical} If TRUE, the function is run using parallel
processing, to speed-up the computation for large rasters (with a maximum
of 8 cores).
The number of cores is automatically determined; specifying it is also
possible (e.g. \code{parallel = 4}). In this case, more than 8 cores can be
specified. If FALSE (default), single core processing is used.}
}
\value{
The function is called for its side effects
}
\description{
Internal function used to perform the required spatial
processing on MODIS original hdf layers (reprojection, resizing, resampling,
mosaicing, computation of scaling factors). The function is based on the
use of \code{gdal} routines.
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbox_from_file.R
\name{bbox_from_file}
\alias{bbox_from_file}
\title{Retrieve bbox from a spatial file}
\usage{
bbox_from_file(file_path, crs_out)
}
\arguments{
\item{file_path}{\code{character} path of a spatial file.}

\item{crs_out}{(\code{crs} | \code{character}) crs of the desired output projection,
or string coercible to it using \code{sf::st_crs()} (e.g., WKT or numeric
EPSG code)}
}
\description{
Helper function used to retrieve the bounding box of a specified spatial file
recognized by  \code{sf} or \code{raster}: the function reads the extent using \code{sf::st_bbox()}
}
\note{
License: GPL 3.0
}
\author{
Lorenzo Busetto, phD (2017)

Luigi Ranghetti, phD (2017) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_extract.R
\name{MODIStsp_extract}
\alias{MODIStsp_extract}
\title{Extract data from MODIStsp time series}
\usage{
MODIStsp_extract(
  in_rts,
  sf_object,
  start_date = NULL,
  end_date = NULL,
  id_field = NULL,
  FUN = "mean",
  out_format = "xts",
  small = TRUE,
  small_method = "centroids",
  na.rm = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{in_rts}{A \code{RasterStack} bject created by MODIStsp
(it MUST contain acquisition dates in the "Z" attribute)}

\item{sf_object}{"sf" object OR name of an GDAL-readable vector file specifying the
"area" from which data has to be extracted.
\itemize{
\item If \code{sf_object} represents lines, the output object contains one column for
each line, containing values obtained applying the function specified
as the FUN argument over all pixels touched by the line, and one line for
each date.
\item If \code{sf_object} represents points, the output object contains one column
for each point, containing values of the cells corresponding to the point,
and one line for each date.
\item If \code{sf_object} represents polygons, the output object contains one column
for each polygon, containing values obtained applying the function
specified as the FUN argument over all pixels belonging to the polygon,
and one line for each date
}}

\item{start_date}{object of class \code{Date}, \code{POSIXct} or \code{POSIXlt} OR \code{character}
coercible to Date class (format = "yyyy-mm-dd")Starting date of the period
to be considered for data extraction . If not provided, the first date of
the RasterStack is used.}

\item{end_date}{object of class \code{Date}, \code{POSIXct} or \code{POSIXlt} OR \code{character}
coercible to Date class (format = "yyyy-mm-dd"). Ending date of the period
to be considered for data extraction . If not provided, the last date of
the RasterStack is used.}

\item{id_field}{\code{character} name of the column of the input sp object or
shapefile to be used in the data extraction. Values contained in the column
MUST be unique. The names of the columns of the output are taken from this
column. If not provided, or an invalid value is provided, then the names
of the columns of the output reflect the number of the feature in
\code{sf_object}.}

\item{FUN}{function to summarize the values (e.g. mean) on polygon data frames.
The function should take a single numeric vector as argument and return a
single value (e.g. mean, min or max), and accept a na.rm argument. Thus,
standard R functions not including an na.rm argument must  be wrapped as in
this example: fun=function(x,...)length(x). Defaults to "mean"}

\item{out_format}{\code{character ["xts" | "dframe"]} If dframe, the output is a
data frame with dates in the first column and extracted data in the others,
otherwise it is a \code{xts} object, Default: "xts"}

\item{small}{\code{logical} If TRUE, and input is polygons, then values are
returned also for polygons not covering at least one raster cell. "Included"
cells in this case depend on the values of the "small_method" parameter.}

\item{small_method}{\code{character ["centroids" | "full"]} If small == TRUE and
input is polygons, controls which cells are "extracted" for small polygons.
If set to "centroids" (default), then only the cells corresponding to polygon
centroid are considered (faster, may have problems on strangely shaped
polygons). If set to "full", then all cells intersected by the small polygon
are extracted and used in calculations, Default: "centroids"}

\item{na.rm}{\code{logical} If TRUE, and sf_object is a polygon, then na.rm = TRUE
is used when applying FUN to the different pixels of the polygon, Default = TRUE.}

\item{verbose}{\code{logical} If TRUE, messages on processing status are sent
to the console. Default = TRUE.}
}
\value{
data.frame or xts object. Each column of data corresponds to one
point or one polygon, each row to a date.
}
\description{
function used to extract time series data from rts files created by MODIStsp
on spatial locations provided in the form of "R" spatial objects (SpatialPoints,
SpatialPolygons, etc.)
}
\details{
The function takes as input a RasterStack object containing time information
in the "z" attribute (set by \code{raster::setZ}), a starting and ending date
and a standard "R" spatial object, and returns the time series for the spatial locations
specified in the spatial object in the form of a "R" xts object OR a plain data.frame
with a "date" column in first position.
If the input spatial object is a "point" or "line" one, the  output object
contains one column for each specified point, or for each cell intersecting
the line, and one line for each date. If the input spatial object is a "polygon"
one, the output object contains one column for each polygon, containing values
obtained applying the function specified as the FUN argument over all pixels
belonging to the polygon, and one line for each date.
}
\note{
License: GPL 3.0
}
\examples{
\dontrun{
# Extract average and standard deviation values from a rts object created by
# MODIStsp for each polygon of a shapefile, for each date in the period
# between 2001-01-01 and 2014-12-31

# The example uses tif files in testdata/VI_16Days_500m_v6 to build
# a MODIStsp rasterStack corresponding to the 2016 time series of the NDVI index
# over the Como Lake (Italy). It then extracts data on polygons corresponding
# to different land cover classes saved in testdata/extract_polys.shp

# First, prepare the test dataset.
# __NOTE__ To avoid redownloading, here we copy some test data from MODIStsp
# installation folder to tempdir and use it to create a test time series.

test_folder <-  system.file("testdata/VI_16Days_500m_v6/NDVI",
                            package = "MODIStsp")
dir.create(file.path(tempdir(), "MODIStsp/VI_16Days_500m_v6/NDVI/"),
           showWarnings = FALSE, recursive = TRUE)
file.copy(list.files(test_folder, full.names = TRUE),
          file.path(tempdir(), "MODIStsp/VI_16Days_500m_v6/NDVI/"))

opts_file <- system.file("testdata/test_extract.json", package = "MODIStsp")
MODIStsp(opts_file = opts_file, gui = FALSE, verbose = FALSE)

# Now load the MODIStsp stack: This is a MODIS NDVI time series ranging between
# 2016-01-01 and 2016-12-18
# __NOTE__: MODIStsp rasterStack files are always saved in the "Time_Series\/RData"
# subfolder of your main output folder - see
# "https://docs.ropensci.org/MODIStsp/articles/output.html")

# Specify the filename of the RData RasterStack of interest
stack_file  <- file.path(tempdir(),
 "MODIStsp/VI_16Days_500m_v6/Time_Series/RData/Terra/NDVI",
  "MOD13A1_NDVI_1_2016_353_2016_RData.RData")
basename(stack_file)

ts_data <- get(load(stack_file))
ts_data

# Now load a shapefile containing polygons from which we want to extract data

polygons <- sf::st_read(system.file("testdata/extract_polys.shp",
                           package = "MODIStsp"), quiet = TRUE)
polygons

# Finally, extract the average values for each polygon and date and plot the
# results

out_dataavg <- suppressMessages(MODIStsp_extract(ts_data, polygons, id_field = "lc_type",
                             small = FALSE))
head(out_dataavg)

plot(out_dataavg, legend.loc = "topleft")

# use a different summarization function

out_datasd <- MODIStsp_extract(ts_data, polygons, id_field = "lc_type",
                              FUN = "sd", small = FALSE)
head(out_datasd)

# (See also https://docs.ropensci.org/MODIStsp/articles/Analyze.html for a
# worked-out example)
}
}
\author{
Lorenzo Busetto, phD (2015 - 2017)
email: busetto.l@irea.cnr.it
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp-package.R
\docType{package}
\name{MODIStsp-package}
\alias{MODIStsp-package}
\title{MODIStsp: a package to automatize the creation of time series of raster
images derived from MODIS Land Products}
\description{
MODIStsp allows automating the creation of time series of rasters derived
from MODIS Satellite Land Products data. It performs several typical
preprocessing steps such as download, mosaicking, reprojection and resize
of data acquired on a specified time period. All processing parameters
can be set using a user-friendly GUI. Users can select which layers of
the original MODIS HDF files they want to process, which additional
Quality Indicators should be extracted from aggregated MODIS Quality
Assurance layers and, in the case of Surface Reflectance products
, which Spectral Indexes should be computed from the original reflectance
bands. For each output layer, outputs are saved as single-band raster
files corresponding to each available acquisition date. Virtual files
allowing access to the entire time series as a single file are also created.
Command-line execution exploiting a previously saved processing options
file is also possible, allowing to automatically update time series
related to a MODIS product whenever a new image is available.
}
\seealso{
\url{https://docs.ropensci.org/MODIStsp/}

\url{https://github.com/ropensci/MODIStsp}
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015-2017) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mod_dates.R
\name{get_mod_dates}
\alias{get_mod_dates}
\title{Find MODIS dates included in selected processing period}
\usage{
get_mod_dates(dates, date_dirs)
}
\arguments{
\item{dates}{2- element string array specifying start/end dates (yyyy.mm.dd)
for which the http addresses of folders in lpdaac should be retrieved
(e.g., c("2015.1.1", "2015.12.31)}

\item{date_dirs}{data frame full list of folders in lpdaac archive for product of interest}
}
\value{
array of folder names containing data for the MODIS product acquired in
the period specified by "dates"
}
\description{
Accessory function to find the folders corresponding to the requested
dates period within the full list retrieved by get_moddirs
}
\note{
License: GPL 3.0
}
\author{
Luigi Ranghetti, phD (2016) \email{luigi@ranghetti.info}

Lorenzo Busetto, phD (2017)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_vrt_create.R
\name{MODIStsp_vrt_create}
\alias{MODIStsp_vrt_create}
\title{Create MODIStsp virtual files}
\usage{
MODIStsp_vrt_create(
  sensor,
  out_prod_folder,
  bandnames,
  bandsel,
  nodata_out,
  indexes_bandnames,
  indexes_bandsel,
  indexes_nodata_out,
  quality_bandnames,
  quality_bandsel,
  quality_nodata_out,
  file_prefixes,
  ts_format,
  out_format,
  verbose
)
}
\arguments{
\item{sensor}{\code{character ["Terra"| "Aqua" | "Both"]} MODIS platform to be considered.
(Ignored for MCD* products). Default: "Both"}

\item{out_prod_folder}{\code{character} Main output folder.}

\item{bandnames}{names of all layers available for the product being processed}

\item{bandsel}{\verb{character array} Original MODIS layers to be processed.
You can get a list of available layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$bandnames),
Default: NULL}

\item{nodata_out}{\code{numeric} Output nodata value to be used in vrt files}

\item{indexes_bandnames}{names of all indexes available for the product being processed}

\item{indexes_bandsel}{\verb{character array}Spectral Indexes to be computed starting from reflectance bands.
You can get a list of available quality layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$indexes_bandnames),
Default: NULL}

\item{indexes_nodata_out}{nodata value for indexes vrts}

\item{quality_bandnames}{names of all quality indicators available for the product being processed}

\item{quality_bandsel}{\verb{character array} Quality Indicators to be computed starting from
bit fields of original MODIS layers. You can get a list of available quality layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$quality_bandnames),
Default: NULL}

\item{quality_nodata_out}{nodata value for quality vrts}

\item{file_prefixes}{\verb{character array (2)} file_prefixes for TERRA and AQUA -
used to identify the files corresponding to each sensor}

\item{ts_format}{\code{character ["ENVI" | "GDAL" | "Both"]} Required output format
for virtual file.}

\item{out_format}{\code{character ["ENVI" | "GTiff"]} Format of images used as
"input" for the vrt and contained in out_prod_folder/band folders.}

\item{verbose}{\code{logical} If FALSE, suppress processing messages, Default: TRUE}
}
\value{
NULL - the function is called for its side effects
}
\description{
Function used to create virtual files from time series of single-band
files corresponding to different acquisition dates. The function takes as input
the folder in which the single-band files are stored, and creates a ENVI Meta
file and/or a GDAL vrt file that allows access to the full time series as if
it was a single physical file.
Created virtual files are stored in the "Time Series" subfolder of `out_prod_folder``
}
\note{
License: GPL 3.0
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mod_filenames.R
\name{get_mod_filenames}
\alias{get_mod_filenames}
\title{Find the names of MODIS images corresponding to the selected dates}
\usage{
get_mod_filenames(
  http,
  used_server,
  user,
  password,
  n_retries,
  date_dir,
  v,
  h,
  tiled,
  out_folder_mod,
  gui
)
}
\arguments{
\item{http}{\code{character} url of http site on lpdaac corresponding to a given MODIS
product.}

\item{used_server}{\code{character} can assume values "http"; it cannot be NA.}

\item{user}{\code{character} username for earthdata server.}

\item{password}{\code{character} password for earthdata server.}

\item{n_retries}{\code{numeric} number of times the access to the http server
should be retried in case of error before quitting, Default: 20.}

\item{date_dir}{\verb{character array} array of folder names corresponding to acquisition
containing dates where MODIS files to be downloaded are to be identified
(return array from \code{get_mod_dates}).}

\item{v}{\verb{integer array} containing a sequence of the vertical tiles of interest
(e.g., c(18,19)).}

\item{h}{\verb{integer array} containing a sequence of the horizontal tiles of interest
(e.g., c(3,4)).}

\item{tiled}{\code{numeric [0/1]} indicates if the product to be downloaded is
tiled or not tiled. 1 = tiled product; 0 = non-tiled product (resolution 0.05 deg).}

\item{out_folder_mod}{\code{character} folder where hdf files have to be stored.}

\item{gui}{\code{logical} indicates if processing was called within the GUI environment
or not. If not, processing messages are redirected direct to the log file.}
}
\value{
\verb{character array} containing names of HDF images corresponding to the
requested tiles available for the product in the selected date
}
\description{
Accessory function to find the names of HDF images corresponding
to a given date and interval of spatial tiles within the lpdaac archive.
}
\note{
License: GPL 3.0
}
\author{
Original code by Babak Naimi (\code{.getModisList}, in
\href{http://r-gis.net/?q=ModisDownload}{ModisDownload.R})
modified to adapt it to MODIStsp scheme and to http archive (instead than old
FTP) by:

Lorenzo Busetto, phD (2014-2016)

Luigi Ranghetti, phD (2016) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp.R
\name{MODIStsp}
\alias{MODIStsp}
\title{MODIStsp main function}
\usage{
MODIStsp(
  ...,
  gui = TRUE,
  out_folder = NULL,
  out_folder_mod = NULL,
  opts_file = NULL,
  selprod = NULL,
  prod_version = NULL,
  bandsel = NULL,
  quality_bandsel = NULL,
  indexes_bandsel = NULL,
  sensor = NULL,
  download_server = NULL,
  downloader = NULL,
  user = NULL,
  password = NULL,
  download_range = NULL,
  start_date = NULL,
  end_date = NULL,
  spatmeth = NULL,
  start_x = NULL,
  end_x = NULL,
  start_y = NULL,
  end_y = NULL,
  bbox = NULL,
  spafile = NULL,
  out_projsel = NULL,
  output_proj = NULL,
  out_res_sel = NULL,
  out_res = NULL,
  resampling = NULL,
  reprocess = NULL,
  delete_hdf = NULL,
  nodata_change = NULL,
  scale_val = NULL,
  ts_format = NULL,
  out_format = NULL,
  compress = NULL,
  test = NULL,
  n_retries = 5,
  verbose = TRUE,
  parallel = TRUE
)
}
\arguments{
\item{...}{not used for values, forces later arguments to bind by name}

\item{gui}{\code{logical} if TRUE: the GUI is opened before processing. If FALSE:
processing parameters are retrieved from the provided \code{opts_file}
argument), Default: TRUE}

\item{out_folder}{\code{character} Main output folder, default: NULL.}

\item{out_folder_mod}{\code{character} Output folder for original HDF storage.
If \code{"$tempdir"} (default), a temporary directory is used.}

\item{opts_file}{\code{character} full path to a JSON file
containing MODIStsp processing options saved from the GUI, Default: NULL}

\item{selprod}{\code{character} Name of selected MODIS product (e.g.,
Vegetation Indexes_16Days_250m (M*D13Q1)). You can get
a list of available product names using function \code{MODIStsp_get_prodnames},
Default: NULL}

\item{prod_version}{Version of the selected MODIS product.
Currently versions \code{"006"} and/or \code{"061"} can be chosen.
Default value is \code{"006"} until decommission of this version will be
announced by USGS.
Products with version \verb{"061} are experimental: in case users would encounter
an error in the encoding of bands or quality flags they are encouraged
to report it by opening a new issue on GitHub at
\url{https://github.com/ropensci/MODIStsp/issues}.}

\item{bandsel}{\verb{character array} Original MODIS layers to be processed.
You can get a list of available layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$bandnames),
Default: NULL}

\item{quality_bandsel}{\verb{character array} Quality Indicators to be computed starting from
bit fields of original MODIS layers. You can get a list of available quality layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$quality_bandnames),
Default: NULL}

\item{indexes_bandsel}{\verb{character array}Spectral Indexes to be computed starting from reflectance bands.
You can get a list of available quality layers for a given product
using function \code{MODIStsp_get_prodlayers} (e.g., MODIStsp_get_prodlayers("M*D13Q1")$indexes_bandnames),
Default: NULL}

\item{sensor}{\code{character ["Terra"| "Aqua" | "Both"]} MODIS platform to be considered.
(Ignored for MCD* products). Default: "Both"}

\item{download_server}{\code{character ["http" | "offline"]} service to be used for
download. Default: "http"}

\item{downloader}{download_server \code{character ["http" | "aria2"]} downloader to be used,
Default: "http"}

\item{user}{\code{character} Username for NASA http server.
(\href{https://urs.earthdata.nasa.gov/home}{urs.earthdata.nasa.gov/home}).}

\item{password}{\code{character} Password for NASA http server
(\href{https://urs.earthdata.nasa.gov/home}{urs.earthdata.nasa.gov/home}).}

\item{download_range}{\code{character ["Full" | "Seasonal"]} If "full", all the
available images between the starting and the ending dates are downloaded;
If "seasonal", only the images included in the season are downloaded
(e.g: if the starting date is 2005-12-01 and the ending is 2010-02-31, only
the images of December, January and February from 2005 to 2010 - excluding
2005-01, 2005-02 and 2010-12 - are downloaded), Default: Full}

\item{start_date}{\code{character} Start date for images download and preprocessing
(yyyy.mm.dd), Default: NULL}

\item{end_date}{\code{character} End date for images download and preprocessing
(yyyy.mm.dd), Default: NULL}

\item{spatmeth}{\code{character ["tiles" | "bbox" | "file"]}, indicates how the processing
extent is retrieves. if "tiles", use the specified tiles (start_x....).
If "file", retrieve extent from spatial file specifies in \code{spafile}. If
"bbox", use the specified bounding box, Default: "tiles"}

\item{start_x}{\code{integer [0-35]} Start MODIS horizontal tile defining spatial extent.
Ignored if spatmeth != "tiles", Default: 18}

\item{end_x}{\code{integer [0-35]} End MODIS horizontal tile defining spatial extent.
Ignored if spatmeth != "tiles", Default: 18}

\item{start_y}{\code{integer [0-17]} Start MODIS vertical tile defining spatial extent.
Ignored if spatmeth != "tiles", Default: 4}

\item{end_y}{\code{integer [0-17]} End MODIS vertical tile defining spatial extent.
Ignored if spatmeth != "tiles", Default: 4}

\item{bbox}{\code{numeric(4)} Output bounding box (xmin, ymin, xmax, ymax) in
out_proj coordinate system. Ignored if spatmeth == "tiles", Default: NULL}

\item{spafile}{\code{character} (optional) full path of a spatial file
to use to derive the processing extent. If not NULL, the processing options
which define the extent, the selected tiles and the "Full Tile / Custom"
in the JSON options file are overwritten and new files are created on the
extent of the provided spatial file. Ignored if spatmeth != "file", Default: NULL}

\item{out_projsel}{\verb{character ["Native", "User Defined}] If "Native", the
outputs keep the original resolution of MODIS HDF images. Otherwise, the value
set in "out_res" is used, Default:Native}

\item{output_proj}{\code{character} either equal to "MODIS Sinusoidal",
or to the code of a valid EPSG or to a WKT projection string.
Ignored if outproj_sel == "Native", Default: NULL}

\item{out_res_sel}{\verb{character ["Native", "User Defined}]. If "Native", the
outputs keep the original resolution of MODIS HDF images. Otherwise, the value
set in "out_res" is used.}

\item{out_res}{\code{float} Output resolution (in output projection measurement
unit). Ignored if out_res_sel == "Native".}

\item{resampling}{\verb{character ["near" | "bilinear" | "cubic" | "cubicspline", |lanczos"|, "average"|, "mode", |"max"|, |"min"|, |"q1"|, |"q3"|, |"sum"|]}
Resampling method to be used by \code{gdalwarp}.}

\item{reprocess}{\code{logical} If TRUE, reprocess data for already existing dates.}

\item{delete_hdf}{\code{logical} If TRUE, delete downloaded HDF files after completion.}

\item{nodata_change}{\code{logical} if TRUE, NoData values are set to the max value
of the datatype of the layer on the MODIStsp output rasters. NOTE: If multiple
nodata values are reported for a layer, all are reset to the new value.}

\item{scale_val}{\code{logical} If TRUE,  scale and offset are applied to
original MODIS layers, and Spectral Indexes are saved as floating point. If
FALSE, no rescaling is done and Spectral Indexes are saved as integer, with a
10000 scaling factor.}

\item{ts_format}{\verb{character array including ["R RasterStack" | "ENVI Meta Files" | "GDAL VRT" | "ENVI and GDAL"]} Selected virtual time series format.}

\item{out_format}{\code{character ["ENVI" | "GTiff"]} Desired output format.}

\item{compress}{\code{character ["None" | "PACKBITS" | "LZW" | "DEFLATE"]}
Compression method for GTiff outputs (Ignored if \code{out_format == ENVI})}

\item{test}{\code{integer | character  (e.g., "01a")} if set, MODIStsp is executed in
"test mode", using a preset Options File instead than opening the GUI or accepting the
\code{opts_file} parameter. This allows both to check correct installation on
user's machines, and to implement unit testing.}

\item{n_retries}{\code{numeric} maximum number of retries on download functions.
In case any download function fails more than \code{n_retries} times consecutively,
MODIStsp_process will abort, Default: 20}

\item{verbose}{\code{logical} If FALSE, suppress processing messages,
Default: TRUE}

\item{parallel}{\code{logical} If TRUE (default), the function is run using parallel
processing, to speed-up the computation for large rasters (with a maximum
of 8 cores).
The number of cores is automatically determined; specifying it is also
possible (e.g. \code{parallel = 4}). In this case, more than 8 cores can be
specified. If FALSE (default), single core processing is used.}
}
\description{
Main function for the MODIS Time Series Processing Tool
(MODIStsp)
}
\details{
The function is used to:
\itemize{
\item initialize the processing (folder names, packages, etc.);
\item launch the GUI (\code{\link[=MODIStsp_GUI]{MODIStsp_GUI()}}) on interactive
execution, or load an options file to set processing arguments and/or
retrieve CLI inputs and run processing on non-interactive execution;
\item launch the routines for downloading and processing the requested datasets.
(\code{\link[=MODIStsp_process]{MODIStsp_process()}})
\item launching the function with GUI = FALSE and without specifying a opts_file
initializes arguments with default values. This allows making a test run.
}
}
\note{
License: GPL 3.0
}
\examples{
\donttest{

#' # - Running the tool using the GUI
# Running the tool without any option will start the GUI with the default or
# last used settings, in interactive mode (i.e., with gui = TRUE).
if (interactive()) {
  MODIStsp()
}


#' # - Running the tool specifying processing arguments in the call

# **NOTE** Output files of examples are saved to file.path(tempdir(), "MODIStsp").

# Here we process layers __NDVI__ and __EVI__ and quality indicator __usefulness__
# of product __M*D13Q1__, considering both Terra and Aqua platforms, for dates
# comprised between 2020-06-01 and 2020-06-15 and saves output to R tempdir
# --> See name and available layers for product M*D13Q1.
# Note that this example (as well as the following ones) is run in single
# core to follow CRAN policies, by setting parallel = FALSE.
# Users can exploit multicore functionalities skipping to set this argument.

MODIStsp_get_prodlayers("M*D13A2")
MODIStsp(
  gui = FALSE,
  out_folder = "$tempdir",
  selprod = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
  bandsel = c("EVI", "NDVI"),
  quality_bandsel = "QA_usef",
  indexes_bandsel = "SR",
  user = "mstp_test" ,
  password = "MSTP_test_01",
  start_date = "2020.06.01",
  end_date = "2020.06.15",
  verbose = FALSE,
  parallel = FALSE
)


#' # - Running the tool using the settings previously saved in a specific options file

# **NOTE** Output files of examples are saved to file.path(tempdir(), "MODIStsp").
# You can run the examples with `gui = TRUE` to set a different output folder!

# Here we use a test json file saved in MODIStsp installation folder which
# downloads and processed 3 MOD13A2 images over the Como Lake (Lombardy, Italy)
# and retrieves NDVI and EVI data, plus the Usefulness Index Quality Indicator.

opts_file <- system.file("testdata/test_MOD13A2.json", package = "MODIStsp")

MODIStsp(gui = FALSE, opts_file = opts_file, verbose = TRUE, parallel = FALSE)


# Running the tool using the settings previously saved in a specific option file
# and specifying the extent from a spatial file allows to re-use the same
# processing settings to perform download and reprocessing on a different area

opts_file <- system.file("testdata/test_MOD13A2.json", package = "MODIStsp")
spatial_file <- system.file("testdata/lakeshapes/garda_lake.shp", package = "MODIStsp")
MODIStsp(
  gui = FALSE, 
  opts_file = opts_file,
  spatmeth = "file",
  spafile = spatial_file, 
  verbose = TRUE,
  parallel = FALSE
)


# Running the tool using the settings previously saved in a
# specific options file and specifying each time the extent from a different
# spatial file (e.g., to perform the same processing on several extents)
# Note that you can also put all your extent files in a specific folder and
# create the extent list using for example.

extent_list = list.files(
  system.file("testdata/lakeshapes/", package = "MODIStsp"),
  "\\\\.shp$", 
  full.names = TRUE
)
extent_list
opts_file <- system.file("testdata/test_MOD13A2.json", package = "MODIStsp")

for (single_shape in extent_list) {
  MODIStsp(
    gui = FALSE, 
    opts_file = opts_file,
    spatmeth = "file",
    spafile = single_shape, 
    verbose = TRUE,
    parallel = FALSE
  )
}

# output files are placed in separate folders:
outfiles_garda <- list.files(
  file.path(tempdir(), "MODIStsp/garda_lake/VI_16Days_1Km_v6/NDVI"),
  full.names = TRUE
)
outfiles_garda
require(raster)
if (length(outfiles_garda) > 0) {
  plot(raster(outfiles_garda[1] ))
}

outfiles_iseo <- list.files(
  file.path(tempdir(), "MODIStsp/iseo_lake/VI_16Days_1Km_v6/NDVI"),
  full.names = TRUE
)
outfiles_iseo
if (length(outfiles_garda) > 0) {
  plot(raster(outfiles_iseo[1]))
}

# See also https://docs.ropensci.org/MODIStsp/articles/noninteractive_execution.html
}
}
\seealso{
\code{\link[=MODIStsp_GUI]{MODIStsp_GUI()}}, \code{\link[=MODIStsp_process]{MODIStsp_process()}}
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015-2017) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_prodopts.R
\name{load_prodopts}
\alias{load_prodopts}
\title{Load characteristics of the different MODIS products}
\usage{
load_prodopts()
}
\value{
OUTPUT_DESCRIPTION
}
\description{
FUNCTION_DESCRIPTION
}
\details{
Load characteristics of the different MODIS products from \code{prodopts_file}
}
\author{
Lorenzo Busetto, phD (2017)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_projection.R
\name{check_projection}
\alias{check_projection}
\alias{check_projection.default}
\alias{check_projection.numeric}
\alias{check_projection.character}
\alias{check_projection.crs}
\title{Check the validity of the input projection}
\usage{
check_projection(projection, abort = FALSE, verbose = TRUE)

\method{check_projection}{default}(projection, abort = FALSE, verbose = TRUE)

\method{check_projection}{numeric}(projection, abort = FALSE, verbose = TRUE)

\method{check_projection}{character}(projection, abort = FALSE, verbose = TRUE)

\method{check_projection}{crs}(projection, abort = FALSE, verbose = TRUE)
}
\arguments{
\item{projection}{\code{character} or \code{integer} corresponding to the
an EPSG code, a UTM zone (e.g. "32N") or a WKT representation of  a projection;}

\item{abort}{\code{logical} if TRUE, the function aborts in case an invalid invalid
projection is passed. Otherwise, the function returns "NA", Default: TRUE}

\item{verbose}{\code{logical} if TRUE, return messages}
}
\value{
\code{character} proj4string of the object or file
}
\description{
helper function used to check that the input projection
(passed as UTM zone, EPSG code, WKT string) is a valid projection for MODIStsp.
}
\note{
This function was forked from package \code{sprawl}, version 0.3.0.
}
\examples{

\dontrun{
check_projection("32632")

check_projection("32631")

check_projection(32633)

check_projection(30, abort = FALSE)

check_projection("example of invalid string", abort = FALSE)

proj_wkt <- sf::st_as_text(sf::st_crs(32632))
check_projection(proj_wkt)
}
}
\author{
Lorenzo Busetto, phD (2017)

Luigi Ranghetti, phD (2017) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_addindex.R
\name{MODIStsp_addindex}
\alias{MODIStsp_addindex}
\title{Add custom spectral indexes}
\usage{
MODIStsp_addindex(
  new_indexbandname = "",
  new_indexfullname = "",
  new_indexformula = "",
  new_indexnodata_out = "32767"
)
}
\arguments{
\item{new_indexbandname}{\code{character} short name (acronym) of the new
spectral index (Ignored if gui == TRUE), Default: NULL}

\item{new_indexfullname}{\code{character} extended name (acronym) of the new
spectral index (Ignored if gui == TRUE), Default: NULL}

\item{new_indexformula}{\code{character} string containing the formula of
the new spectral indexes (Ignored if gui == TRUE). Variables allowed in
the formula are the names of the bands:
b1_Red, b2_NIR, b3_Blue, b4_Green, b5_SWIR, b6_SWIR and b7_SWIR.
Default: NULL}

\item{new_indexnodata_out}{\code{character} nodata value to use for rasters
containing the new index}
}
\value{
The function is called for its side effects. On success, the
MODIStsp_indexes.json
is modified so to allow computation of the additional indexes.
}
\description{
Function used to add a user-defined Spectral Index to the
default list of computable spectral indexes. Execution without the GUI
(i.e., to add a new index from a script) is also possible (see examples).
}
\details{
\itemize{
\item The function asks the user to provide the info related to the new desired
Spectral Index, checks for correctness of provided
information (e.g., correct bandnames, computable formula, etc...).
If the index is legit, it modifies the MODIStsp_addindex.json file
so to allow computation of the additional index within MODIStsp for all
products containing the required reflectance bands.
\item To remove all custom-added spectral indexes, run MODIStsp_resetindexes()
}
}
\note{
License: GPL 3.0
}
\examples{
# Run the GUI to interactively define a new index
 \dontrun{
 MODIStsp_addindex()}

# Define the new index in non-interactive execution

\dontrun{
MODIStsp_addindex(new_indexbandname = "SSI",
  new_indexfullname = "Simple Useless Index",
  new_indexformula = "b2_NIR+b1_Red")
}
}
\seealso{
\link{MODIStsp_resetindexes}
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_get_prodlayers.R
\name{MODIStsp_get_prodnames}
\alias{MODIStsp_get_prodnames}
\title{Retrieve the names of all available product}
\usage{
MODIStsp_get_prodnames()
}
\value{
character array of product names
}
\description{
Function used to retrieve the names of available MODIS products
}
\note{
License: GPL 3.0
}
\examples{
# Get MODIStsp product names
MODIStsp_get_prodnames()
}
\author{
Lorenzo Busetto, phD (2014-2020)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_read_xml.R
\name{MODIStsp_read_xml}
\alias{MODIStsp_read_xml}
\title{Read MODIS products characteristics from XML}
\usage{
MODIStsp_read_xml(prodopts_file, xml_file)
}
\arguments{
\item{prodopts_file}{string filename of the RData in which to store the data
parsed from the XML file}

\item{xml_file}{string filename of the XML file containing the MODIS products
characteristics}
}
\value{
NULL - retrieved data are stored in the specified RData file
}
\description{
function used to parse the XML file used to store the characteristics of
MODIS Land Products and store them in the "prod_opts" data frame
}
\details{
The function parses the XML file product by product, stores data in a data frame
and saves the data frame within the "MODIStsp_previous" RData file as a list of lists
}
\note{
License: GPL 3.0
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_process.R
\name{MODIStsp_process}
\alias{MODIStsp_process}
\title{MODIStsp main processing function}
\usage{
MODIStsp_process(proc_opts, n_retries, verbose = TRUE, parallel = TRUE)
}
\arguments{
\item{proc_opts}{\code{data.frame} containing all processing parameters, as
passed from the MODIStsp GUI, or created in \code{MODIStsp} by joining
explicitly passed arguments with a (not mandatory) options file.}

\item{n_retries}{\code{numeric} maximum number of retries on download functions.
In case any download function fails more than \code{n_retries} times consecutively,
MODIStsp_process will abort, Default: 20}

\item{verbose}{\code{logical} If FALSE, suppress processing messages, Default: TRUE}

\item{parallel}{\code{logical} If TRUE (default), the function is run using parallel
processing, to speed-up the computation for large rasters (with a maximum
of 8 cores).
The number of cores is automatically determined; specifying it is also
possible (e.g. \code{parallel = 4}). In this case, more than 8 cores can be
specified. If FALSE (default), single core processing is used.}
}
\value{
The function is called for its side effects.
}
\description{
Main processing function of MODIStsp. Takes as input processing
parameters specified by the user and performs all required processing.
}
\details{
After retrieving the input processing options, the function
\enumerate{
\item Accesses NASA http archive to determine the list of files to be
downloaded/processed (or, in case of offline processing, get the list
of already available hdf files present in \code{out_mod_folder});
\item Performs all required processing steps on each date (download,
reprojection, resize, mosaicing, Spectral Indexes and Quality indicators
computation);
\item Creates virtual files of the processed time series.
}

Reprojection and resize is dealt with by accessing gdal routines through the
\href{https://CRAN.R-project.org/package=gdalUtilities}{\code{gdalUtilities}}
package.
Extraction of bitfields from Quality layers is done though bitwise computation
Checks are done in order to not re-download already existing HDF images, and not
reprocess already processed dates (if the user did not specify that)
}
\note{
Thanks Tomislav Hengl and Babak Naimi, whose scripts made the starting point for
development of this function (\href{http://r-gis.net/?q=ModisDownload}{ModisDownload};
\href{https://en.wikipedia.org/wiki/Regression-kriging?title=Download_and_resampling_of_MODIS_images}{Download_and_resampling_of_MODIS_images})

License: GPL 3.0
}
\author{
Lorenzo Busetto, phD (2014-2017)

Luigi Ranghetti, phD (2015) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_get_prodlayers.R
\name{MODIStsp_get_prodlayers}
\alias{MODIStsp_get_prodlayers}
\title{Retrieve the names of MODIS layers for a product}
\usage{
MODIStsp_get_prodlayers(prodname, prodver = "006")
}
\arguments{
\item{prodname}{character containing the code of the desired MODIS product.
NOTE: for products available separately for Terra and Aqua (e.g., MOD13Q1/MYD13Q1),
use the code M\emph{D_code_ (e.g., M}D13Q1)}

\item{prodver}{character containing the version of the desired MODIS product.
Currently versions \code{"006"} (default) and \code{"061"} can be chosen.}
}
\value{
list, containing the slots: \code{prodname}, \code{bandnames}, \code{quality_bandnames} and
\code{indexes_bandnames}, \code{band_fullnames}, \code{quality_fullnames}, \code{indexes_fullnames}
}
\description{
Function used to retrieve the names of original MODIS layers, quality layers
and eventually available spectral indexes given a MODIS product code.
It is useful to identify the names of the layers to be processed when
launching MODIStsp from the CLI.
}
\note{
License: GPL 3.0
}
\examples{

# Get layers of product M*13Q1 based on code
MODIStsp_get_prodlayers("M*13Q1")

# Get layers of product M*13Q1 based on full name
MODIStsp_get_prodlayers("Vegetation Indexes_16Days_250m (M*D13Q1)")

# Get indexes names of product M*13Q1 based on full name
MODIStsp_get_prodlayers("MCD43C4")$indexes_bandnames


}
\author{
Lorenzo Busetto, phD (2014-2020)

Luigi Ranghetti, phD (2021)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIStsp_process_QA_bits.R
\name{MODIStsp_process_QA_bits}
\alias{MODIStsp_process_QA_bits}
\title{MODIStsp helper function to compute Quality Indicators from HDF bit-field layers}
\usage{
MODIStsp_process_QA_bits(
  out_filename,
  in_source_file,
  bitN,
  out_format,
  nodata_source,
  nodata_qa_in,
  nodata_qa_out,
  compress
)
}
\arguments{
\item{out_filename}{\code{character} file name of the output raster files
containing QI values}

\item{in_source_file}{\code{character} name of the file created by MODIStsp
containing the data required to compute the quality indicator}

\item{bitN}{\code{character} position of the bits corresponding to the quality
indicator of interest (e.g., 0-1 = first two bits; 2-5: bits from 2 to 5,
etc.)}

\item{out_format}{output format (ENVI or GTiff)}

\item{nodata_source}{\code{character} NoData values of the MODIS band containing
data from which the bit field corresponding to the quality indicator must
be extracted}

\item{nodata_qa_in}{\code{character} in NoData for quality bands ("255")}

\item{nodata_qa_out}{\code{character} out NoData for quality bands ("255")}

\item{compress}{\code{character} compression option for GTiff files}
}
\description{
function used to extract quality indicator from MODIS aggregated
quality layers
}
\details{
On the basis of the name of the image containing the aggregated quality information
(\verb{in_source_file``) and of the position of the bit fields corresponding to the QI of interest in the bitfield representation (}bitN``), the function extracts the correct information exploiting
bitwise operators, and save the result in a new raster image
}
\note{
License: GPL 3.0
Based on the \code{modis.qc.R} script by Yann Chemin (2008) (\url{https://goo.gl/7Fhreo})
license GPL 3.0
}
\author{
Lorenzo Busetto, phD (2017)

Luigi Ranghetti, phD (2017) \email{luigi@ranghetti.info}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_message.R
\name{process_message}
\alias{process_message}
\title{Spawn processing update messages}
\usage{
process_message(mess_text, verbose = TRUE)
}
\arguments{
\item{mess_text}{\code{character} text to be shown in the processing windows
and/or the console}

\item{verbose}{\code{logical} If FALSE, suppress processing messages, Default: TRUE}
}
\value{
The function is called for its side effects
}
\description{
helper function to provide processing messages
}
\author{
Lorenzo Busetto, phD (2017)
}
