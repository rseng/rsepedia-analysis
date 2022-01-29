rinat: Access iNaturalist data with R
================
Edmund Hart, Stéphane Guillou

[![Build
Status](https://api.travis-ci.org/ropensci/rinat.png)](https://travis-ci.org/ropensci/rinat)
[![Build
status](https://ci.appveyor.com/api/projects/status/gv7s9um107bep4na/branch/master)](https://ci.appveyor.com/project/sckott/rinat/branch/master)
[![codecov.io](https://codecov.io/github/ropensci/rinat/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rinat?branch=master)
[![](https://cranlogs.r-pkg.org/badges/rinat)](https://CRAN.R-project.org/package=rinat)

R wrapper for iNaturalist APIs for accessing the observations. The
detailed documentation of the API is available on the [iNaturalist
website](https://www.inaturalist.org/pages/api+reference) and is part of
our larger species occurrence searching packages
[SPOCC](https://github.com/ropensci/spocc).

## Installation

You can install the latest version available on CRAN with:

``` r
install.packages("rinat")
```

Alternatively, you can install the development version from Github with:

``` r
remotes::install_github("ropensci/rinat")
```

## Usage

### Get observations

#### Text search

You can search for observations by either common or scientific name. It
will search the entire iNaturalist database, so the search below will
return all entries that *mention* Monarch butterflies, not just Monarch
observations.

``` r
library(rinat)
monarchs <- get_inat_obs(query = "Monarch Butterfly")
unique(monarchs$scientific_name)
```

    ## [1] "Danaus plexippus" "Danaina"

> Note that `get_inat_obs()` will return 100 observations by default.
> This can be controlled with the `maxresults` argument.

Another use for a fuzzy search is searching for a habitat,
e.g. searching for all observations that might happen in a vernal pool.
We can then see all the taxon names found.

``` r
vp_obs <- get_inat_obs(query = "vernal pool")
# see the first few taxa
head(vp_obs$scientific_name)
```

    ## [1] "Micrathetis triplex"     "Sphaeropthalma unicolor"
    ## [3] "Lepidoptera"             "Lepidoptera"            
    ## [5] "Synchlora faseolaria"    "Lepidoptera"

#### Taxon search

To return only records of a specific species or taxonomic group, use the
`taxon_name` argument. For example, to return observations of anything
from the Nymphalidae family, and restricting the search to the year
2015:

``` r
nymphalidae <- get_inat_obs(taxon_name  = "Nymphalidae", year = 2015)
# how many unique taxa?
length(unique(nymphalidae$scientific_name))
```

    ## [1] 79

And to return only the Monarch butterfly observations that also mention
the term “chrysalis”:

``` r
monarch_chrysalis <- get_inat_obs(taxon_name = "Danaus plexippus", query = "chrysalis")
```

#### Bounding box search

You can also search within a bounding box by giving a simple set of
coordinates.

``` r
## Search by area
bounds <- c(38.44047, -125, 40.86652, -121.837)
deer <- get_inat_obs(query = "Mule Deer", bounds = bounds)
plot(deer$longitude, deer$latitude)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Other functions

More functions are available, notably to access:

  - observations in a project with `get_inat_obs_project()`
  - details of a single observation with `get_inat_obs_id()`
  - observations from a single user with `get_inat_obs_user()`
  - taxa statistics with `get_inat_taxon_stats()`
  - user statistics with `get_inat_user_stats()`

More detailed examples are included in the vignette:

``` r
vignette("rinat-intro", package = "rinat")
```

#### Mapping

Basic maps can be created as well to quickly visualize search results.
Maps can either be plotted automatically with `plot = TRUE` (the
default), or simply return a ggplot2 object with `plot = FALSE`. This
works well with single species data, but more complicated plots are best
made from scratch.

``` r
library(ggplot2)

## Map 100 spotted salamanders
a_mac <- get_inat_obs(taxon_name = "Ambystoma maculatum")
salamander_map <- inat_map(a_mac, plot = FALSE)

### Now we can modify the returned map
salamander_map + borders("state") + theme_bw()
```

<img src="README_files/figure-gfm/unnamed-chunk-9-1.png" width="672" />

`inat_map()` is useful for quickly mapping data obtained with rinat.
Here is an example of customised map that does not make use of it. (Not
the use of `quality = "research"` to restrict the search to the more
reliable observations.)

``` r
## A more elaborate map of Colibri sp.
colibri <- get_inat_obs(taxon_name = "Colibri",
                        quality = "research",
                        maxresults = 500)
ggplot(data = colibri, aes(x = longitude,
                         y = latitude,
                         colour = scientific_name)) +
  geom_polygon(data = map_data("world"),
                   aes(x = long, y = lat, group = group),
                   fill = "grey95",
                   color = "gray40",
                   size = 0.1) +
  geom_point(size = 0.7, alpha = 0.5) +
  coord_fixed(xlim = range(colibri$longitude, na.rm = TRUE),
              ylim = range(colibri$latitude, na.rm = TRUE)) +
  theme_bw()
```

<img src="README_files/figure-gfm/unnamed-chunk-10-1.png" width="672" />

-----

[![](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
# rinat 0.1.8

* Properly cater for mismatch in project observation numbers to stay on CRAN (attempt in previous version did not cover all cases)
* Improve console messages in `get_inat_obs_project()`
* Skip tests that use the iNaturalist API on CRAN
* Don't error when iNaturalist or Internet is not available, as per CRAN Repository Policies (#49)

# rinat 0.1.7

* Cater for a corner case in which a mismatch between a project's reported number of observations and the actual number of observations returned produced an error when merging data frames in `get_inat_obs_project()`

# rinat 0.1.6

* New maintainer: Stéphane Guillou (@stragu)
* Improve documentation
* Clean up code according to CRAN feedback
* Note that 0.1.x versions of rinat will only fix existing issues, before a move to the new iNaturalist API with version 0.2 (which is likely to introduce breaking changes).

## New features

* `get_inat_obs()` can now use objects of class `Spatial` or `sf` as a bounding box (#38, @Martin-Jung and #40, @LDalby; fixes #30)
* Allow the use of an iNat `place_id` when searching for observations (commit 1d7f14f, @LDalby)

## Bug fixes

* Lower result number limit when bounding box present (#20)
* Fix result pagination in `get_obs_project()` for cases when the number of results is a multiple of 100 (#1, @vijaybarve)
* Stop `get_inat_obs()` if no search parameter(s) provided (#32, @LDalby)
* Avoid API error by limiting search results to 10000 (#32, @LDalby)
* Fix argument name in `inat_handle()` (#32, @LDalby)
* Use JSON endpoint instead of CSV to avoid `get_inat_obs_project()` failing on some projects (#35, @LDalby and #38, @Martin-Jung; fixes #28 and #37)
* Fix code according to `devtools::check()` feedback in preparation for CRAN release

# rinat 0.1.5

## Bug fixes

* Fixed bug where an error occurred when >20K records were requested and code now throws a warning to not hammer the API
* Fixed warning thrown when building vignettes on Linux systems
* Fixed bug where example code parameter names were different than actual parameter names

## New features

* Added NEWS file.
* Added a full suite of tests
* Added new vignette that builds with markdown and not hacky prebuilt PDF## Test environments

* local Ubuntu 18.04, R 4.0.4
* win-builder R-devel with `devtools::check_win_devel()`
* win-builder R-release with `devtools::check_win_release()`

## R CMD check results

There were no ERRORs, no WARNINGs

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'St�phane Guillou <stephane.guillou@member.fsf.org>'
  New submission
  Package was archived on CRAN
  Possibly mis-spelled words in DESCRIPTION:
  APIs (3:42)
  CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2021-03-04 for policy violation.
  
Misspelling is irrelevant, and reason for archival is addressed in this release.