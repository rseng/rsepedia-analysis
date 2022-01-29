
<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![peer-review](https://badges.ropensci.org/237_status.svg)](https://github.com/ropensci/software-review/issues/237)
[![CRAN](https://www.r-pkg.org/badges/version/spatsoc)](https://cran.r-project.org/package=spatsoc)
[![](https://img.shields.io/badge/devel%20version-0.1.16-blue.svg)](https://github.com/robitalec/spatsoc)
[![cran
checks](https://cranchecks.info/badges/summary/spatsoc)](https://cran.r-project.org/web/checks/check_results_spatsoc.html)
[![codecov](https://codecov.io/gh/ropensci/spatsoc/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/spatsoc)
[![R-CMD-check](https://github.com/ropensci/spatsoc/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/spatsoc/actions)
<!-- badges: end -->

# spatsoc

### [News](#news) - [Installation](#installation) - [Usage](#usage) - [Contributing](#contributing)

`spatsoc` is an R package for detecting spatial and temporal groups in
GPS relocations. It can be used to convert GPS relocations to
gambit-of-the-group format to build proximity-based social networks with
grouping and edge-list generating functions. In addition, the
`randomizations` function provides data-stream randomization methods
suitable for GPS data and the `get_gbi` function generates group by
individual matrices useful for building networks with
`asnipe::get_network`.

See below for [installation](#installation) and basic [usage](#usage).

For more details, see the [blog
post](https://ropensci.org/blog/2018/12/04/spatsoc/) and vignettes:

-   [Introduction to
    spatsoc](https://docs.ropensci.org/spatsoc/articles/intro-spatsoc.html)
-   [Frequently asked
    questions](https://docs.ropensci.org/spatsoc/articles/faq.html)
-   [Using spatsoc in social network
    analysis](https://docs.ropensci.org/spatsoc/articles/using-in-sna.html)
-   [Using edge list and dyad id
    functions](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html)

## News

We wrote a [`targets`](https://github.com/ropensci/targets) workflow,
available at
[github.com/robitalec/targets-spatsoc-networks](https://github.com/ropensci/targets).
`targets` is an incredible package for designing workflows in R and,
with it, we can reproducibly run all steps from raw telemetry data to
output networks and metrics. Check it out and let us know how it works
for you!

Edge-list generating functions added:

-   `edge_nn`
-   `edge_dist`

and dyad id function:

-   `dyad_id`

(feedback welcome as always!)

Both documented further in a vignette: [Using edge list and dyad id
functions](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html).

Also, our article describing `spatsoc` is published at Methods in
Ecology and Evolution. [Link
here](https://doi.org/10.1111/2041-210X.13215). Thanks to reviewers and
editors at
[rOpenSci](https://github.com/ropensci/software-review/issues/237) and
at [MEE](https://besjournals.onlinelibrary.wiley.com/journal/2041210x).

More detailed news
[here](https://docs.ropensci.org/spatsoc/news/index.html).

## Installation

``` r
# Stable release
install.packages('spatsoc')

# Development version
remotes::install_github('ropensci/spatsoc')
```

`spatsoc` depends on `rgeos` and requires
[GEOS](https://trac.osgeo.org/geos/) installed on the system.

-   Debian/Ubuntu: `apt-get install libgeos-dev`
-   Arch: `pacman -S geos`
-   Fedora: `dnf install geos geos-devel`
-   Mac: `brew install geos`
-   Windows: see [here](https://trac.osgeo.org/osgeo4w/)

## Usage

### Load package, import data

`spatsoc` expects a `data.table` for all of its functions. If you have a
`data.frame`, you can use `data.table::setDT()` to convert it by
reference. If your data is a text file (e.g.: CSV), you can use
`data.table::fread()` to import it as a `data.table`.

``` r
library(spatsoc)
library(data.table)
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]
```

### Temporal grouping

`group_times` groups rows temporally using a threshold defined in units
of minutes (B), hours (C) or days (D).

<img src="man/figures/fig1.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

### Spatial grouping

`group_pts` groups points spatially using a distance matrix (B) and a
spatial threshold defined by the user (50m in this case). Combined with
`group_times`, the returned ‘group’ column represents spatiotemporal,
point based groups (D).

<img src="man/figures/fig2.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

`group_lines` groups sequences of points (forming a line) spatially by
buffering each line (A) by the user defined spatial threshold. Combined
with `group_times`, the returned ‘group’ column represents
spatiotemporal, line overlap based groups (B).

<img src="man/figures/fig3.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

`group_polys` groups home ranges by spatial and proportional overlap.
Combined with `group_times`, the returned ‘group’ column represents
spatiotemporal, polygon overlap based groups.

<img src="man/figures/fig4.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

### Edge-list generating functions

`edge_dist` and `edge_nn` generate edge-lists. `edge_dist` measures the
spatial distance between individuals (A) and returns all pairs within
the user specified distance threshold (B). `edge_nn` measures the
distance between individuals (C) and returns the nearest neighbour to
each individual (D).

<img src="man/figures/fig5.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

### Social network analysis functions

`randomizations` for data-stream randomization and `get_gbi` for
generating group by individual matrices.

# Contributing

Please note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree
to abide by its terms.

Development of `spatsoc` welcomes contribution of feature requests, bug
reports and suggested improvements through the [issue
board](https://github.com/ropensci/spatsoc/issues).

See details in [CONTRIBUTING.md](CONTRIBUTING.md).

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# v 0.1.17 (unreleased)

* added a link to our `spatsoc` + `targets` workflow example
* changed the error and underlying check for `group_polys` from alphanumeric to
spaces in input DT's id column

# v 0.1.16 (2021-03-23)
* added an option for `edge_dist` to handle threshold = NULL. If NULL, `edge_dist` will return all neighbours observed (eg. useful if one wanted to calculated mean nearest neighbour distance at each timegroup). 
* updated EPSG argument according to newest recommendations in tests, man and vignettes ([PR 38](https://github.com/ropensci/spatsoc/pull/38)
* removed expect_silent tests ([PR 37](https://github.com/ropensci/spatsoc/pull/37))
* switched CI for tests and code coverage to GitHub Actions ([PR 36](https://github.com/ropensci/spatsoc/pull/36))


# v 0.1.15 (2020-10-21)
* fix TZ=UTC data.table tests ([Issue 32](https://github.com/ropensci/spatsoc/issues/32))

# v 0.1.14 (2020-07-03)
* updated tests, man and vignettes following new handling of projections in sp ([PR 31](https://github.com/ropensci/spatsoc/pull/31), [R spatial information](https://www.r-spatial.org/r/2020/03/17/wkt.html))
* clarified explicit drop of NAs in dyadID in edge list vignette

# v 0.1.13 (2020-03-25)
* added `dyad_id` function for generating dyad IDs with edge functions ([PR 27](https://github.com/ropensci/spatsoc/pull/25))
* added a vignette describing `edge_dist`, `edge_nn` and `dyad_id` functions [here](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html) ([PR 14](https://github.com/ropensci/spatsoc/pull/14))

# v 0.1.12 (2020-03-02)
* fixed `data.table` error in `edge_dist` and `edge_nn` ([PR 25](https://github.com/ropensci/spatsoc/pull/25))


# v 0.1.11 (2020-02-20)
* removed default NULL from 'timegroup' arguments in `group_pts`, `edge_dist` and `edge_nn` ([PR 24](https://github.com/ropensci/spatsoc/pull/24))


# v 0.1.10 (2019-06-06)
* added optional return of distance between individuals with `edge_dist` ([PR 19](https://github.com/ropensci/spatsoc/pull/19)) and `edge_nn` ([PR 21](https://github.com/ropensci/spatsoc/pull/21))


# v 0.1.9 (2019-05-14)
* fixed bug for randomizations type 'step' and 'daily' ([PR 13](https://github.com/ropensci/spatsoc/pull/13)). 
* clarified `SIMPLIFY=FALSE` in SNA vignette. 


# v 0.1.8 (2019-04-05)
* update [FAQ](https://docs.ropensci.org/spatsoc/articles/faq.html) and [Introduction to spatsoc](https://docs.ropensci.org/spatsoc/articles/intro-spatsoc.html) vignettes adding entries for edge list generating functions. 
* added edge list generating function `edge_nn` ([PR 11](https://github.com/ropensci/spatsoc/pull/12))
* added edge list generating function `edge_dist` ([PR 11](https://github.com/ropensci/spatsoc/pull/11))


# v 0.1.7 (2019-03-26)
* fix inconsistent blocks across years ([PR 10](https://github.com/ropensci/spatsoc/pull/10))
* update FAQ: remove old randomizations notes, clarify group_times block


# v 0.1.6 (2019-01-10)
* fix bug 'group_times misses nearest hour with mins threshold' ([#5](https://github.com/ropensci/spatsoc/issues/5) and [PR 6](https://github.com/ropensci/spatsoc/pull/6))

# v 0.1.5 (2018-12-04)
* update issue labels and contributing
* change over issue board location from GitLab to rOpenSci repository on GitHub
* added preprint CITATION
* added "https://" to `pkgdown` URL ([PR 1](https://github.com/ropensci/spatsoc/pull/1))

# v 0.1.4 (2018-10-26)
* fin [rOpenSci onboarding process](https://github.com/ropensci/software-review/issues/237)
* fixed bug couldn't provide percent to kernel type `build_polys` or `group_polys`([!3](https://gitlab.com/robit.a/spatsoc/-/merge_requests/3))


# v 0.1.3 
* added `get_gbi` to generate group by individual matrices for better integrating `spatsoc` in social network analysis workflows ([!2](https://gitlab.com/robit.a/spatsoc/-/merge_requests/2))


# v 0.1.2

* **major change to randomizations**: when `iterations = 1`, `randomizations` no longer returns the DT with appended columns. Regardless of the value of iterations, `randomizations` always returns observed rows followed by randomized rows in a long `data.table` ([!1](https://gitlab.com/robit.a/spatsoc/-/merge_requests/1)). 

# v 0.1.1 (2018-09-17)

* improvements to package, function documentation
* [FAQ](https://docs.ropensci.org/spatsoc/articles/faq.html) vignette added
* fixed `build_lines` ordering bug to ensure rows are ordered by date time when building lines
* added CODE_OF_CONDUCT.md and CONTRIBUTING.md
* [Using spatsoc in social network analysis](https://docs.ropensci.org/spatsoc/articles/using-in-sna.html) vignette added

# v 0.1.0 (2018-07-20)

## Initial release

* temporal grouping function: `group_times`
* spatial grouping functions: `group_pts`, `group_lines`, `group_polys`
* data-stream randomization function: `randomizations`
* spatial build functions: `build_lines`, `build_polys`
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
# Contributing

Development of `spatsoc` welcomes contribution of feature requests, bug reports and suggested improvements through the [issue board](https://github.com/ropensci/spatsoc/issues). 

### Issues

Use the labels available on the Issue Board:

#### priority
 
| priority  | description | colour                                          	| 
|-----------|-------------|---------------------------------------------------|
| critical 	|             | <span style="background:#e50000">#e50000</span> 	|
| high 	    |             | <span style="background:#fc7436">#fc7436</span> 	|
| medium 	  |             | <span style="background:#ffbb5c">#ffbb5c</span> 	|
| low 	    |             | <span style="background:#dfdfdf">#dfdfdf</span> 	|


#### status

| status      	| description                                	| colour                                            |
|-------------	|--------------------------------------------	|--------------------------------------------------	|
| completed   	|                                            	| <span style="background:#34495E">#34495E</span> 	|
| confirmed   	| bug is reproducible                        	| <span style="background:#e50000">#e50000</span> 	|
| in progress 	| [workinonit](https://youtu.be/5nO7IA1DeeI) 	| <span style="background:#aefd3d">#aefd3d</span> 	|
| duplicate   	|                                            	| <span style="background:#dfdfdf">#dfdfdf</span> 	|

#### type

|  type         	| description                        | colour                                           |
|---------------	|----------------------------------- |-------------------------------------------------	|
| documentation 	|             	                     | <span style="background:#639f3b">#639f3b</span> 	|
| support       	| usage or implementation troubles   | <span style="background:#e48ec0">#e48ec0</span>  |
| bug           	| not working as described, intended | <span style="background:#e50000">#e50000</span> 	|
| discussion    	|             	                     | <span style="background:#47b6a3">#47b6a3</span> 	|
| enhancement   	|             	                     | <span style="background:#8935b8">#8935b8</span> 	|
| typo          	|             	                     | <span style="background:#dfdfdf">#dfdfdf</span> 	|


### Prefer to Email? 

Email Alec Robitaille (see DESCRIPTION). 

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Code of Conduct

Please note that the `spatsoc` project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md
## v 0.1.16
See NEWS.md for all updates. 

* added an option for `edge_dist` to handle threshold = NULL. If NULL, `edge_dist` will return all neighbours observed (eg. useful if one wanted to calculated mean nearest neighbour distance at each timegroup). 
* updated EPSG argument according to newest recommendations in tests, man and vignettes ([PR 38](https://github.com/ropensci/spatsoc/pull/38)
* removed expect_silent tests ([PR 37](https://github.com/ropensci/spatsoc/pull/37))
* switched CI for tests and code coverage to GitHub Actions ([PR 36](https://github.com/ropensci/spatsoc/pull/36))

## Test environments
* windows-latest (release) 
* macOS-latest (release)
* ubuntu-20.04 (release)
* ubuntu-20.04 (devel)
* ubuntu 16.04 (3.5)

## R CMD check results

There were no ERRORs, WARNING or NOTES
