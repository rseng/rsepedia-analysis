
[![Travis build
status](https://travis-ci.org/ropensci/landscapetools.svg?branch=master)](https://travis-ci.org/ropensci/landscapetools)
[![Build
status](https://ci.appveyor.com/api/projects/status/aehfkxfb5r4vjlm9?svg=true)](https://ci.appveyor.com/project/ropensci/landscapetools)
[![codecov](https://codecov.io/gh/ropensci/landscapetools/branch/develop/graph/badge.svg)](https://codecov.io/gh/ropensci/landscapetools)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
status](https://www.r-pkg.org/badges/version/landscapetools)](https://cran.r-project.org/package=landscapetools)
[![](http://cranlogs.r-pkg.org/badges/grand-total/landscapetools)](http://cran.rstudio.com/web/packages/landscapetools/index.html)
[![](https://badges.ropensci.org/188_status.svg)](https://github.com/ropensci/onboarding/issues/188)
[![DOI:10.1111/2041-210X.13076](https://zenodo.org/badge/DOI/10.1111/2041-210X.13076.svg)](https://doi.org/10.1111/2041-210X.13076)

# landscapetools

`landscapetools` provides utility functions for some of the
less-glamorous tasks involved in landscape analysis:

#### Utilities:

  - `util_binarize`: Binarize continuous raster values, if \> 1 breaks
    are given, return a RasterBrick.
  - `util_classify`: Classify a raster into proportions based upon a
    vector of class weightings.
  - `util_merge`: Merge a primary raster with other rasters weighted by
    scaling factors.
  - `util_raster2tibble`, `util_tibble2raster`: Coerce raster\* objects
    to tibbles and vice versa.
  - `util_rescale`: Linearly rescale element values in a raster to a
    range between 0 and 1.
  - `util_writeESRI`: Export raster objects as ESRI asciis (with Windows
    linebreaks).

#### Visualization

  - `show_landscape`: Plot a Raster\* object with the landscapetools
    default theme (as ggplot) or multiple raster (RasterStack, -brick or
    list of raster) side by side as facets.
  - `show_shareplot`: Plot the landscape share in subsequential buffers
    around a/multiple point(s) of interest

#### Themes:

  - `theme_nlm`, `theme_nlm_grey`: Opinionated ggplot2 theme to
    visualize raster (continuous data).
  - `theme_nlm_discrete`, `theme_nlm_grey_discrete`: Opinionated ggplot2
    theme to visualize raster (discrete data).
  - `theme_faceplot`: Opinionated ggplot2 theme to visualize raster in a
    facet wrap.

## Installation

You can install the released version from CRAN with:

``` r
install.packages("landscapetools")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/landscapetools")
```

## Utilities

### Classify

``` r
# Classify the landscape into land uses
classified_landscape <- util_classify(fractal_landscape,
                                      n = 3,
                                      level_names = c("Land Use 1", 
                                                      "Land Use 2",
                                                      "Land Use 3"))

show_landscape(classified_landscape, discrete = TRUE)
```

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

### Merge

``` r
# Merge all landscapes into one
merged_landscape <- util_merge(fractal_landscape,
                               c(gradient_landscape, random_landscape),
                               scalingfactor = 1)

# Plot an overview
merge_vis <- list(
    "1) Primary" = fractal_landscape,
    "2) Secondary 1" = gradient_landscape,
    "3) Secondary 2" = random_landscape,
    "4) Result" = merged_landscape
)

show_landscape(merge_vis)
#> Warning: Removed 1196 rows containing missing values (geom_raster).
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

## See also

In the examples above we make heavy use of the `NLMR` package. Both
packages were developed together until we split them into pure landscape
functionality and utility tools. If you are interested in generating
neutral landscapes via a multitude of available algorithms take a closer
look at the [NLMR](https://github.com/ropensci/NLMR/) package.

## Meta

  - Please [report any issues or
    bugs](https://github.com/ropensci/landscapetools/issues/new/).
  - License: GPL3
  - Get citation information for `landscapetools` in R doing
    `citation(package = 'landscapetools')`
  - We are very open to contributions - if you are interested check
    [Contributing](CONTRIBUTING.md).
      - Please note that this project is released with a [Contributor
        Code of Conduct](CODE_OF_CONDUCT.md). By participating in this
        project you agree to abide by its
terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# landscapetools 0.6.2
- Bugfix in `util_classify`

# landscapetools 0.6.0
- `util_raster2tibble` can now return a wide tibble
- New function `show_shareplot`
- `util_as_integer` now returns integer values from 1:n instead of rounding numeric values

# landscapetools 0.5.0
- new interface for `util_classify`
    - now takes argument n to specify number of classes
    - n argument implemented in C++
- Removed Roboto font and `util_import_roboto`
- Removed `util_plot_grey`
- Renamed:
    - `util_plot` to `show_landscape`
- new function `util_writeESRI` that produces a replica of esris ascii file format

# landscapetools 0.4.0

* minor bug fixes
* util_facetplot now better handles lists of raster
* improved theme_facetplot
* util_classify can now reclassify based on real landscapes, the classification then overwrites the weightings with the proportions from this landscape
* util_classify now has an mask argument, that allows for the classification only outside this mask
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# CONTRIBUTING #

### Please contribute!

We love collaboration.

### Bugs?

* Submit an issue on the Issues page [here](https://github.com/marcosci/landscapetools/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/landscapetools.git`
* Make sure to track progress upstream (i.e., on our version of `landscapetools` at `marcosci/landscapetools`) by doing `git remote add upstream https://github.com/marcosci/landscapetools.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new branch)
* If you alter package functionality at all (e.g., the code itself, not just documentation)
please do write some tests to cover the new functionality.
* Push up to your account
* Submit a pull request to home base at `marcosci/landscapetools`

### Questions? Get in touch: [sciaini.marco@gmail.com](mailto:sciaini.marco@gmail.com)

### Thanks for contributing!
## Submission

New version includes massive dependency trimming and new functionality.

## Test environments

* local Ubuntu Linux 16.04 LTS install, R 3.4.1
* Ubuntu 14.04 (on travis-ci), R 3.4.1
* Windows Server 2012 R2 x64 (build 9600) (on appveyor), R 3.4.2
* Rhub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Fedora Linux, R-devel, clang, gfortran
* macOS 10.11 El Capitan, R-release
* macOS 10.9 Mavericks, R-oldrel
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

## Reverse dependencies

There are currently no reverse dependencies.
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/{owner}/{repo}/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/{repo}.git`
* Make sure to track progress upstream (i.e., on our version of `{repo}` at `{owner}/{repo}`) by doing `git remote add upstream https://github.com/{owner}/{repo}.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `{owner}/{repo}`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email?

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
