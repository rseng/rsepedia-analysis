
<!-- README.md is generated from README.Rmd. Please edit that file -->

# osfr <a href="https://docs.ropensci.org/osfr"><img src="man/figures/logo.png" align="right" height="139" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/osfr)](https://CRAN.R-project.org/package=osfr)
[![Build
Status](https://travis-ci.com/ropensci/osfr.svg)](https://travis-ci.com/ropensci/osfr)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ropensci/osfr?branch=master&svg=true)](https://ci.appveyor.com/project/aaronwolen/osfr)
[![Coverage
status](https://codecov.io/gh/ropensci/osfr/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/osfr?branch=master)
[![](https://badges.ropensci.org/279_status.svg)](https://github.com/ropensci/software-review/issues/279)
[![](https://joss.theoj.org/papers/d5398fc36ea92794a20914143d3fcdc4/status.svg)](https://joss.theoj.org/papers/d5398fc36ea92794a20914143d3fcdc4)
[![DOI](https://zenodo.org/badge/42329785.svg)](https://zenodo.org/badge/latestdoi/42329785)

## Overview

osfr provides a suite of functions for interacting with the Open Science
Framework ([OSF](https://osf.io "Open Science Framework")).

**What is OSF?**

*OSF is a free and [open
source](https://github.com/CenterForOpenScience/osf.io "OSF's GitHub Repository")
project management repository designed to support researchers across
their entire project lifecycle. The service includes unlimited cloud
storage and file version history, providing a centralized location for
all your research materials that can be kept private, shared with select
collaborators, or made publicly available with citable DOIs.*

## Installation

You can install the current release of osfr from CRAN (*recommended*):

``` r
install.packages("osfr")
```

Or the development version from GitHub with the *remotes* package:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/osfr")
```

## Usage Examples

*Note: You need to [setup an OSF personal access token
(PAT)](https://docs.ropensci.org/osfr/articles/auth) to use osfr to
manage projects or upload files.*

### Accessing Open Research Materials

Many researchers use OSF to archive and share their work. You can use
osfr to explore publicly accessible projects and download the associated
files—all you need to get started is the project’s URL or GUID (global
unique identifier).

Every user, project, component, and file on OSF is assigned a GUID that
is embedded in the corresponding entity’s URL. For example, you can
access the main OSF project for the *Cancer Reproducibility Project* at
<https://osf.io/e81xl/>. The GUID for this project is `e81xl`.

We can then use osfr to *retrieve* this project and load it into R by
providing the GUID:

``` r
library(osfr)

cr_project <- osf_retrieve_node("e81xl")
cr_project
#> # A tibble: 1 x 3
#>   name                                    id    meta            
#>   <chr>                                   <chr> <list>          
#> 1 Reproducibility Project: Cancer Biology e81xl <named list [3]>
```

This returns an `osf_tbl` object with a single row representing the
retrieved project. Let’s list the files that have been uploaded to this
project.

``` r
osf_ls_files(cr_project)
#> # A tibble: 4 x 3
#>   name                                     id                     meta          
#>   <chr>                                    <chr>                  <list>        
#> 1 Adjustment of 50 studies to 37 studies.… 565602398c5e4a3877d72… <named list […
#> 2 papers_and_keywords.xlsx                 553e671b8c5e4a219919e… <named list […
#> 3 Full_dataset_of_papers_formatted.xls     553e671b8c5e4a219919e… <named list […
#> 4 METHOD_to_select_papers.txt              553e671b8c5e4a219919e… <named list […
```

This returns another `osf_tbl` with 1 row for each of the files and
directories in the project. We can examine any of these files directly
on OSF with `osf_open()`, which opens the corresponding file’s view in
your default browser.

This project contains 2 ***components***: *Replication Studies* and
*Data collection and publishing guidelines*. We can list these
components with osfr using `osf_ls_nodes()`.

``` r
osf_ls_nodes(cr_project)
#> # A tibble: 2 x 3
#>   name                                      id    meta            
#>   <chr>                                     <chr> <list>          
#> 1 Replication Studies                       p7ayb <named list [3]>
#> 2 Data collection and publishing guidelines a5imq <named list [3]>
```

osfr is compatible with the [pipe
operator](https://magrittr.tidyverse.org) and
[dplyr](https://dplyr.tidyverse.org), providing a powerful set of tools
for working with `osf_tbl`s. Here, we’re listing the sub-components
nested within the *Replication Studies* component, filtering for a
specific study ([*Study 19*](https://osf.io/7zqxp/)) and then listing
the files uploaded to that study’s component.

``` r
library(dplyr)

cr_project %>%
  osf_ls_nodes() %>%
  filter(name == "Replication Studies") %>%
  osf_ls_nodes(pattern = "Study 19") %>%
  osf_ls_files()
#> # A tibble: 6 x 3
#>   name                                    id                     meta           
#>   <chr>                                   <chr>                  <list>         
#> 1 Replication_Study_19.docx               57c9e8ed594d9001e7a24… <named list [3…
#> 2 Replication_Study_19.Rmd                578e2b23594d9001f4816… <named list [3…
#> 3 Replication_Study_19_track_changes.docx 581a27b76c613b0223322… <named list [3…
#> 4 Replication_Study_19_track_changes_2.d… 58714d46594d9001f801f… <named list [3…
#> 5 Response_letter_Replication_Study_19.d… 58755747b83f6901ff066… <named list [3…
#> 6 Study_19_Correction_Letter.docx         5a56569125719b000ff28… <named list [3…
```

We could continue this pattern of exploration and even download local
copies of project files using `osf_download()`. Or, if you come across a
publication that directly references a file’s OSF URL, you could quickly
download it to your project directory by providing the URL or simply the
GUID:

``` r
osf_retrieve_file("https://osf.io/btgx3/") %>%
  osf_download()
#> # A tibble: 1 x 4
#>   name                id                    local_path            meta          
#>   <chr>               <chr>                 <chr>                 <list>        
#> 1 Study_19_Figure_1.… 5751d71d9ad5a1020793… ./Study_19_Figure_1.… <named list […
```

### Managing Projects

You can use osfr to [create
projects](https://docs.ropensci.org/osfr/reference/osf_create), [add
sub-components](https://docs.ropensci.org/osfr/reference/osf_create) or
[directories](https://docs.ropensci.org/osfr/reference/osf_mkdir), and
[upload files](https://docs.ropensci.org/osfr/reference/osf_upload). See
[Getting
Started](https://docs.ropensci.org/osfr/articles/getting_started) to
learn more about building projects with osfr, but here is a quick
example in which we:

1.  Create a new project called *Motor Trend Car Road Tests*
2.  Create a sub-component called *Car Data*
3.  Create a directory named *rawdata*
4.  Upload a file (`mtcars.csv`) to the new directory
5.  Open the uploaded file on OSF

<!-- end list -->

``` r
# create an external data file
write.csv(mtcars, "mtcars.csv")

osf_create_project(title = "Motor Trend Car Road Tests") %>%
  osf_create_component("Car Data") %>%
  osf_mkdir("rawdata") %>%
  osf_upload("mtcars.csv") %>%
  osf_open()
```

![Screenshot of the uploaded file on OSF](man/figures/screen-shot.png)

## Details on `osf_tbls`

There are 3 main types of OSF entities that osfr can work with:

1.  **nodes:** both
    [projects](https://help.osf.io/hc/en-us/articles/360019737594-Create-a-Project "OSF: Create a Project")
    and
    [components](https://help.osf.io/hc/en-us/articles/360019737614-Create-Components "OSF: Create a Component")
    (i.e., sub-projects) are referred to as nodes
2.  **files:** this includes both files *and* folders stored on OSF
3.  **users:** individuals with OSF accounts

osfr represents these entities within `osf_tbl`s—specialized data frames
built on the tibble class that provide useful information about the
entities like their `name` and unique `id` for users, and API data in
the `meta` column that’s necessary for osfr’s internal functions.
Otherwise, they’re just `data.frames` and can be manipulated using
standard functions from base R or dplyr.

## Acknowledgments

OSF is developed by the [Center for Open
Science](https://cos.io "Center for Open Science") in Charlottesville,
VA.

The original version of osfr was developed by [Chris
Chartgerink](https://github.com/chartgerink) and further developed by
[Brian Richards](https://github.com/bgrich) and [Ryan
Hafen](https://github.com/hafen). The current version was developed by
[Aaron Wolen](https://github.com/aaronwolen) and is *heavily* inspired
by [Jennifer Bryan](https://github.com/jennybc) and [Lucy D’Agostino
McGowan](https://github.com/lucymcgowan)’s excellent
[googledrive](https://googledrive.tidyverse.org) package. Seriously, we
borrowed a lot of great ideas from them. Other important resources
include [http testing](https://books.ropensci.org/http-testing/) by
Scott Chamberlain and [R Packages](http://r-pkgs.had.co.nz) by Hadley
Wickham. Development was also greatly facilitated by OSF’s excellent
[API documentation](https://developer.osf.io "OSF API Documentation").

Big thanks to Rusty Speidel for designing our logo and [Tim
Errington](https://github.com/timerrington) for his feedback during
development.

## Contributing

Check out the [Contributing
Guidelines](https://github.com/ropensci/osfr/blob/master/.github/CONTRIBUTING.md)
to get started with osfr development and note that by contributing to
this project, you agree to abide by the terms outlined in the
[Contributor Code of
Conduct](https://github.com/ropensci/osfr/blob/master/.github/CODE_OF_CONDUCT.md).

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

<!-- links -->
# osfr (development version)

## Minor changes

* tibble v3.0.0 is now the minimum required version

## Fixes

* Fixed bug preventing uploads directly to OSF directories that contained conflicting files (#121, #129)

# osfr 0.2.8

* Initial CRAN release
* Publication of accompanying paper in the [Journal of Open Source Software](http://joss.theoj.org/) that can be cited in papers using osfr, see `citation("osfr")` for details

## Minor changes

* Add rOpenSci reviewers to DESCRIPTION
* Remove deleted URLs from vignette
* Add badges for zenodo and JOSS
* Add `Makefile` for common developer tasks

# osfr 0.2.7

## Important changes

osfr is now part of rOpenSci and the documentation website has moved to a new URL:  <https://docs.ropensci.org/osfr>.

## New features

* New `osf_cp()` function for copying files to new locations (@tpyork, #114)

## Other changes

* The *getting started* vignette was overhauled to better leverage multi-file transfers and is now precomputed
* Encoded HTML symbols in node titles are now handled properly (#117)
* `osf_rm()` argument `recursive` been renamed to `recurse` in order to be consistent with other functions
* Internal links now point to the ropensci repository and new documentation URL

# osfr 0.2.6

## Improved uploading

* New approach uses a file manifest to compare local and remote files
* We now search for conflicting remote files within each directory queued for upload, this avoids the issue reported in (#108, thanks @tpyork)

## osfr 0.2.6.1

* Fix pkgdown deployment

## osfr 0.2.6.2

* Fix project creation on windows (#110)
* Fix tests when PAT is undefined

## osfr 0.2.6.3

* Add JOSS paper

# osfr 0.2.5

## Multi-file transfers!

`osf_download()` and `osf_upload()` are now vectorized, making the process of adding files to or retrieving files from OSF much more convenient. This functionality required significant refactoring and brings with it several notable breaking changes (see below).

## Other new features

* `osf_download()` and `osf_upload()` gain the option to display progress bars.
* New `osf_refresh()` to update an existing `osf_tbl`.
* Devs can now enable logging API requests and responses by defining`OSF_LOG` (see Contributing for more information).

## Breaking changes

* `osf_download()` and `osf_upload()`'s `overwrite` argument has been replaced with `conflicts`, which can be set to `"error"` (the default), `"skip"`, or `"overwrite"`.
* `osf_upload()`'s `name` argument has been removed, so it is no longer possible to upload a file *and* change it's OSF name.
* `osf_download()`'s `path` argument must point to an existing directory where all downloaded files will be saved.
* `osf_download()`'s `decompress` argument has been removed. The zip file downloaded from OSF is always decompressed in a temp directory where the enclosed files are selectively copied to the specified `path`.

## Minor changes

* Better error message when user attempts to upload directly to a file
(#102, @tiernanmartin).
* crul v0.7.4 is now the minimum required version.
* The waterbutler client will now re-attempt failed requests 3 times.
* Consolidated internal client constructors.
* Increased wait time on travis to avoid time outs during testing.

# osfr 0.2.4

## Minor fixes

* Listing files within a specified `path` would fail if sibling directories
shared a common substring in their names (#95)
* Setting `verbose=TRUE` now works properly for `osf_upload()`
* A startup message is printed when `OSF_SERVER` is defined
* Improved documentation for `n_max`, GUIDs and the mysterious `meta` column

# osfr 0.2.3

## New features

* Failed OSF API requests are now re-attempted 3 times (requires crul v0.7.0)

## Minor fixes

* Fix incorrect column name in empty `osf_tbl`s (#88, @machow)
* No longer importing `modify_at()`
* Add rOpenSci badge (#89, @maelle)
* Don't build vignettes on travis

# osfr 0.2.2

## New features

* Added `osf_mv()` to move files and directories to a new project, component, or
subdirectory
* `osf_rm()` can now delete files and directories

## Minor improvements and fixes

* Restructured tests to better handle environments in which `OSF_PAT` and/or `OSF_SERVER` are not defined

# osfr 0.2.1

* Minor tweaks to the website
* `osf_retrieve_file()` will no longer retrieve files on 3rd-party storage
providers, since other osfr functions currently only support OSF storage

# osfr 0.2.0

**NOTE:** This version of osfr is a rewrite of the original codebase. It is
effectively an entirely different package and provides no backwards
compatibility with functions in versions < 0.2.0. The last version of the
previous package can be installed with the *remotes* package:

```r
remotes::install_github("ropensci/osfr@v0.1.1")
```

See <https://docs.ropensci.org/osfr> for details about the new
package.
## Test environments
* local OS X install, R 3.6
* ubuntu 14.04 (on travis-ci),  R 3.5.1, devel and release
* win-builder (devel and release)
* Windows Server 2008 R2 SP1 (on R-hub), R-devel
* Ubuntu Linux 16.04 LTS + GCD (on R-hub), R-release
* Fedora Linux + clang + gfortran (on R-hub), R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission. In this version I have:

- replaced github links in README.md with fully specified URLs
- updated the TITLE and DESCRIPTION to place 'osfr', 'OSF', and 'Open Science Framework' in quotes
- removed the full license from builds
- fixed invalid README URLs, which resulted from migrating to rOpenSci's docs server
- added a CITATION for the recently published JOSS paper
- switched to the MIT license

## Reverse dependencies

There are no reverse dependencies.
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
# Contributing to osfr

This outlines how to propose a change to osfr.

### Development environment

To get started with osfr development you'll need to generate a personal access token (PAT) on OSF's *testing* server. The following steps will get you setup:

1. Create an account on <https://test.osf.io/>.
2. Generate a PAT for the new account with all permission scopes enabled (this is necessary to run osfr's unit tests). (Under "Settings"->Personal access token)
3. Fork the osfr repository and clone a local copy.
4. Create a `.Renviron` file in the root of your project directory that defines the `OSF_PAT` and `OSF_SERVER` environment variables. You can easily create or edit an existing `.Renviron` file by running `usethis::edit_r_environ(scope = "project")`. The end result should look like this:

   ```
   OSF_PAT=<YOUR PAT GOES HERE>
   OSF_SERVER=test
   ```

5. Restart R to load your `.Renviron` file or use `readRneviron(".Renviron")`.  Then load your local copy of osfr with `devtools::load_all()` and verify that `osf_open(osf_retrieve_user("me"))` opens your user profile on the `test.osf.io` domain.

Once this is setup correctly, you should be able to run osfr's tests without error (`devtools::test()`).

You can also enable logging by defining `OSF_LOG` to point to a logfile. For example:

```
OSF_PAT=osfr.log
```

This will log all API requests to `osfr.log` for inspection.

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.
*  New code should follow the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with
your PR.
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html),
for documentation.
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Using the included `Makefile`

Run `make docs` to:
* rebuild `README.md` from `README.Rmd`
* regenerate the precomputed vignette `getting_started.Rmd` from `getting_started.Rmd.orig`
* rebuild the vignettes
* rebuild package documentation

Run `make` to
* perform all of the documentation steps noted above
* build the package
* check the package as CRAN but without running tests (this is temporary until mock tests are implemented)

**Helpers:**

* `make clean` to remove build/check files
* `make test` to run unit tests
* `make tag` to tag the last git commit with the current version

---
title: 'osfr: An R Interface to the Open Science Framework'
tags:
  - R
  - open science
  - reproducible research
  - project management
authors:
  - name: Aaron R. Wolen
    orcid: 0000-0003-2542-2202
    affiliation: 1
  - name: Chris H.J. Hartgerink
    orcid: 0000-0003-1050-6809
    affiliation: 2
  - name: Ryan Hafen
    affiliation: 3
  - name: Brian G. Richards
    affiliation: 4
  - name: Courtney K. Soderberg
    orcid: 0000-0003-1227-7042
    affiliation: 5
  - name: Timothy P. York
    orcid: 0000-0003-4068-4286
    affiliation: 6
affiliations:
  - name: Transplant Research Institute, Department of Surgery, University of Tennessee Health Science Center
    index: 1
  - name: Liberate Science GmbH
    index: 2
  - name: Department of Statistics, Purdue University
    index: 3
  - name: Merkle Group Inc.
    index: 4
  - name: Center for Open Science
    index: 5
  - name: Data Science Lab, Department of Human and Molecular Genetics, Virginia Commonwealth University
    index: 6
date: 27 November 2019
bibliography: paper.bib
---

# Background

Reproducible research requires effective project management workflows that promote consistency and facilitate transparency. Hallmarks of effective workflows include strategies for tracking the provenance of results, recording intermediate changes, conveniently documenting procedures, and working with collaborators without duplicating effort [@Sandve:2013]. For technically skilled researchers, the combination of version control software (VCS) such as git and cloud-based project repositories (e.g., GitHub and GitLab) enable highly effective workflows [@Ram:2013de] that facilitate automation and computational reproducibility. However, these tools have a steep learning curve, especially for researchers whose training is far removed from software development. Alternatively, the Open Science Framework (OSF) offers much of the same functionality through an intuitive point-and-click web-based interface, significantly lowering the barrier to adopting best practices for researchers of all skill levels, or research groups composed of individuals with different levels of computational expertise [@Sullivan:2019]. Yet the increase in accessibility comes at the cost of limiting opportunities for automation. `osfr` fills this gap for R users by allowing them to programmatically interact with OSF through a suite of functions for managing their projects and files.

# Functional Overview

On OSF, individual repositories are referred to as *projects* and serve as the top-level unit of content organization. New projects can be created with osfr using `osf_create_project()`, which allows you to specify the project's title, provide a description, and indicate whether it should be private (the default) or publicly accessible. Every OSF project includes a cloud-based storage bucket where files can be stored and organized into directories. You can use `osf_mkdir()` to add directories and `osf_upload()` to populate the project with files. osfr supports recursively uploading nested directories, making it possible to easily mirror the contents of a local project on OSF. For existing projects and project files, osfr provides functions for most common file operations such as copying (`osf_cp()`), moving (`osf_mv()`), deleting (`osf_rm()`), and downloading (`osf_download()`).

A key organizational feature of OSF is the ability to augment a project's structure with sub-projects, which are referred to as *components* on OSF, and can be added with the `osf_create_component()` function. Like top-level projects, every component is assigned a unique URL upon creation and contains its own cloud-based storage bucket, activity log, wiki, and user permissions. This allows users to create arbitrarily nested projects that can easily scale to meet the needs of even large, multi-institutional collaborations (see the [Cancer Biology Reproducibility Project][cbrp] for a great example). Whatever the scale of your work, adopting a consistent structure across projects creates predictable expectations, facilitates understanding for you and your collaborators [@Wilson:2017], and makes it easier to stay organized as a project inevitably grows in complexity over time. Maintaining a consistent structure can be a challenge, especially if implemented in an *ad hoc* process, but osfr enables you to codify your preferred organizational structure of components and directories in a simple script that can be run at the outset of every new project.

## Implementation and Design

osfr is built on the OSF public REST API, available at https://developer.osf.io, and uses rOpenSci's HTTP client, crul [@crul], for API Communication. In order to provide an interface that feels natural to R users, items retrieved from the OSF are represented as `data.frame`-like objects called `osf_tbl`s. The `osf_tbl` class is built on top of the [tibble package][tibble] [@tibble] and, like [googledrive][]'s dribble class [@googledrive], uses a list-column to encapsulate JSON responses from the API. These deeply nested structures are rarely of interest to the end user but are essential for the package's internal methods. The vast majority of osfr functions return `osf_tbl`s as output and expect them as input, so that method chaining is possible using [magrittr][]'s pipe operator [@magrittr]. 

Exported osfr functions all start with the prefix, `osf_`, following the `<prefix>_<verb>` naming convention used in packages like [stringr][] [@stringr], which facilitates auto-completion in supported IDEs (like RStudio) and avoids namespace clashes with other packages that perform similar file-based operations. Where possible, we adopt the names of common Unix utilities that perform analogous tasks (e.g., `osf_cp()`, `osf_mkdir()`).

# Summary

`osfr` provides an idiomatic R interface to OSF (Open Science Framework, https://www.osf.io), a free and open source web application that is part open-access repository and part collaborative project management tool.

# References

<!-- link -->
[cbrp]: https://osf.io/e81xl/ "Reproducibility Project: Cancer Biology"
[googledrive]: https://googledrive.tidyverse.org
[magrittr]: https://magrittr.tidyverse.org
[stringr]: https://stringr.tidyverse.org
[tibble]: https://tibble.tidyverse.org
This file is here so that R CMD CHECK doesn't complain about an empty "inst" directory, but this directory must be there because we need the .lintr file so that codecov works.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)

Sys.unsetenv(c("OSF_SERVER", "OSF_PAT"))
```

# osfr <a href="https://docs.ropensci.org/osfr"><img src="man/figures/logo.png" align="right" height="139" /></a>

[![CRAN status](https://www.r-pkg.org/badges/version/osfr)](https://CRAN.R-project.org/package=osfr)
[![Build Status](https://travis-ci.com/ropensci/osfr.svg)](https://travis-ci.com/ropensci/osfr)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/osfr?branch=master&svg=true)](https://ci.appveyor.com/project/aaronwolen/osfr)
[![Coverage status](https://codecov.io/gh/ropensci/osfr/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/osfr?branch=master)
[![](https://badges.ropensci.org/279_status.svg)](https://github.com/ropensci/software-review/issues/279)
[![](https://joss.theoj.org/papers/d5398fc36ea92794a20914143d3fcdc4/status.svg)](https://joss.theoj.org/papers/d5398fc36ea92794a20914143d3fcdc4)
[![DOI](https://zenodo.org/badge/42329785.svg)](https://zenodo.org/badge/latestdoi/42329785)

## Overview

osfr provides a suite of functions for interacting with the Open Science Framework ([OSF][osf]).

**What is OSF?**

*OSF is a free and [open source][osf-gh] project management repository designed to support researchers across their entire project lifecycle. The service includes unlimited cloud storage and file version history, providing a centralized location for all your research materials that can be kept private, shared with select collaborators, or made publicly available with citable DOIs.*

## Installation

You can install the current release of osfr from CRAN (*recommended*):

```r
install.packages("osfr")
```

Or the development version from GitHub with the *remotes* package:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/osfr")
```

## Usage Examples

*Note: You need to [setup an OSF personal access token (PAT)][auth] to use osfr to manage projects or upload files.*

### Accessing Open Research Materials

Many researchers use OSF to archive and share their work. You can use osfr to explore publicly accessible projects and download the associated files---all you need to get started is the project's URL or GUID (global unique identifier).

Every user, project, component, and file on OSF is assigned a GUID that is embedded in the corresponding entity's URL. For example, you can access the main OSF project for the *Cancer Reproducibility Project* at <https://osf.io/e81xl/>. The GUID for this project is `e81xl`.

We can then use osfr to *retrieve* this project and load it into R by providing the GUID:

```{r message=FALSE}
library(osfr)

cr_project <- osf_retrieve_node("e81xl")
cr_project
```

This returns an `osf_tbl` object with a single row representing the retrieved project. Let's list the files that have been uploaded to this project.

```{r}
osf_ls_files(cr_project)
```

This returns another `osf_tbl` with 1 row for each of the files and directories in the project. We can examine any of these files directly on OSF with `osf_open()`, which opens the corresponding file's view in your default browser.

This project contains 2 ***components***: *Replication Studies* and *Data collection and publishing guidelines*. We can list these components with osfr using `osf_ls_nodes()`.

```{r}
osf_ls_nodes(cr_project)
```

osfr is compatible with the [pipe operator][magrittr] and [dplyr][], providing a powerful set of tools for working with `osf_tbl`s. Here, we're listing the sub-components nested within the *Replication Studies* component, filtering for a specific study ([*Study 19*](https://osf.io/7zqxp/)) and then listing the files uploaded to that study's component.

```{r message=FALSE}
library(dplyr)

cr_project %>%
  osf_ls_nodes() %>%
  filter(name == "Replication Studies") %>%
  osf_ls_nodes(pattern = "Study 19") %>%
  osf_ls_files()
```

We could continue this pattern of exploration and even download local copies of project files using `osf_download()`. Or, if you come across a publication that  directly references a file's OSF URL, you could quickly download it to your project directory by providing the URL or simply the GUID:

```{r}
osf_retrieve_file("https://osf.io/btgx3/") %>%
  osf_download()
```


### Managing Projects

You can use osfr to [create projects][osf-create], [add sub-components][osf-create] or [directories][osf-mkdir], and [upload files][osf-upload]. See [Getting Started][getting-started] to learn more about building projects with osfr, but here is a quick example in which we:

1. Create a new project called *Motor Trend Car Road Tests*
2. Create a sub-component called *Car Data*
3. Create a directory named *rawdata*
4. Upload a file (`mtcars.csv`) to the new directory
5. Open the uploaded file on OSF

```{r eval=FALSE}
# create an external data file
write.csv(mtcars, "mtcars.csv")

osf_create_project(title = "Motor Trend Car Road Tests") %>%
  osf_create_component("Car Data") %>%
  osf_mkdir("rawdata") %>%
  osf_upload("mtcars.csv") %>%
  osf_open()
```

![Screenshot of the uploaded file on OSF](man/figures/screen-shot.png)

## Details on `osf_tbls`

There are 3 main types of OSF entities that osfr can work with:

1. **nodes:** both [projects][help-proj] and [components][help-comp] (i.e., sub-projects) are referred to as nodes
2. **files:** this includes both files *and* folders stored on OSF
3. **users:** individuals with OSF accounts

osfr represents these entities within `osf_tbl`s---specialized data frames built on the tibble class that provide useful information about the entities like their `name` and unique `id` for users, and API data in the `meta` column that's necessary for osfr's internal functions. Otherwise, they're just `data.frames` and can be manipulated using standard functions from base R or dplyr.

## Acknowledgments

OSF is developed by the [Center for Open Science][cos] in Charlottesville, VA.

The original version of osfr was developed by [Chris Chartgerink][chris] and further developed by [Brian Richards][brian] and [Ryan Hafen][ryan]. The current version was developed by [Aaron Wolen][aaron] and is *heavily* inspired by [Jennifer Bryan][jenny] and [Lucy D'Agostino McGowan][lucy]'s excellent [googledrive][] package. Seriously, we borrowed a lot of great ideas from them. Other important resources include [http testing](https://books.ropensci.org/http-testing/) by Scott Chamberlain and [R Packages](http://r-pkgs.had.co.nz) by Hadley Wickham. Development was also greatly facilitated by OSF's excellent [API documentation][osf-api].

Big thanks to Rusty Speidel for designing our logo and [Tim Errington][tim] for his feedback during development.

## Contributing

Check out the [Contributing Guidelines][contrib] to get started with osfr development and note that by contributing to this project, you agree to abide by the terms outlined in the [Contributor Code of Conduct][coc].

```{r cleanup, include=FALSE}
unlink("Study_19_Figure_1.pdf")
```

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

<!-- links -->
[osf]: https://osf.io "Open Science Framework"
[cos]: https://cos.io "Center for Open Science"
[osf-gh]: https://github.com/CenterForOpenScience/osf.io "OSF's GitHub Repository"
[osf-api]: https://developer.osf.io "OSF API Documentation"
[help]: http://help.osf.io "OSF Support"
[help-proj]: https://help.osf.io/hc/en-us/articles/360019737594-Create-a-Project "OSF: Create a Project"
[help-comp]: https://help.osf.io/hc/en-us/articles/360019737614-Create-Components "OSF: Create a Component"
[magrittr]: https://magrittr.tidyverse.org
[dplyr]: https://dplyr.tidyverse.org
[googledrive]: https://googledrive.tidyverse.org
[tibble]: https://tibble.tidyverse.org

[chris]: https://github.com/chartgerink
[brian]: https://github.com/bgrich
[ryan]: https://github.com/hafen
[aaron]: https://github.com/aaronwolen
[jenny]: https://github.com/jennybc
[lucy]: https://github.com/lucymcgowan
[tim]: https://github.com/timerrington

[getting-started]: https://docs.ropensci.org/osfr/articles/getting_started
[auth]: https://docs.ropensci.org/osfr/articles/auth

[osf-create]: https://docs.ropensci.org/osfr/reference/osf_create
[osf-mkdir]: https://docs.ropensci.org/osfr/reference/osf_mkdir
[osf-upload]: https://docs.ropensci.org/osfr/reference/osf_upload

[contrib]: https://github.com/ropensci/osfr/blob/master/.github/CONTRIBUTING.md
[coc]: https://github.com/ropensci/osfr/blob/master/.github/CODE_OF_CONDUCT.md
---
title: "Authenticating osfr"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Authenticating osfr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval = FALSE, comment = "#>", collapse = TRUE)
```

## osfr

This vignette covers how to obtain and use an OSF *personal access token* (PAT) for use with `osfr`.

PATs are required to upload files, create projects/components, access information about your private projects, or download files in your private projects. PATs are not required for accessing information about public projects or downloading public files, but authentication with a PAT will increase the rate limit on the API.

Unauthenticated requests are limited to 100 per hour; Authenticated requests are limited to 10,000 per day. So, if you are planning to make a large number of calls to the API, we suggest authenticating with a PAT even if one is not required for the functions you are running.

## Creating an OSF PAT

To create an OSF PAT, log into OSF through your browser, navigate to the [OSF settings page](https://osf.io/settings/tokens/), and click the **create token** button. You will need to create a PAT with at least `osf.read_full` scope; if you want to be able to edit information through the package, you will also need the PAT to have `osf.write_full` permissions. Once the PAT is created, save it in a safe space.

## Using your PAT

You can authenticate in two different ways:

1. Call the `osf_auth()` function in the console at the start of each new session and paste in your PAT. 

   ```r
   osf_auth("ThIsIsNoTaReAlPaTbUtYoUgEtIt")
   ```

   Your PAT functions like a password, so it should **not** be hard coded into any scripts, ever.

2. To authenticate automatically, store the PAT as an environment variable named `OSF_PAT`, which osfr will detect upon being loaded. One way to do this is to create a `.Renviron` file in your home or project working directory that defines `OSF_PAT`. 

    If you'd like to learn more about `.Renviron` files, the [R Startup chapter](https://rstats.wtf/r-startup.html) of *What They Forgot to Teach You About R* is highly recommended.

## Removing a PAT

If your PAT has accidentally been publicly released in the world, you should deactivate that PAT. To do this, navigate to the [OSF settings page](https://osf.io/settings/tokens/) and click on the :x: to the right of the PAT you want to deactivate. You can then create a new PAT and reauthenticate.
---
title: "Getting Started with osfr"
date: "2020-02-06"
output:
  rmarkdown::html_vignette:
    toc: TRUE
vignette: >
  %\VignetteIndexEntry{Getting Started with osfr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



This vignette provides a quick tour of the *osfr* package.


```r
library(osfr)
```



## What is OSF?

[OSF][osf] is a free and open source web application that provides a space for researchers to collaboratively store, manage, and share their research materials (e.g. data, code, protocols).

Most work on OSF is organized around ***projects***, which include a cloud-based storage bucket where files can be stored and organized into directories. Note there is no storage limit on the size of projects but individual files must be < 5Gb. Projects can be kept private, shared with a specific group of collaborators, or made publicly available with citable DOIs so you can get credit for their work.

If you'd like to learn more about OSF the Center for Open Science has published an excellent series of [guides](http://help.osf.io/) to help get you started. We'll provide links to specific guides throughout this vignette.
Here are a few relevant topics:

* [Creating and Managing Projects][osf-projects]
* [Adding sub-projects (i.e., components) to a project][osf-components]
* [Understanding project permissions][osf-permissions]

## Accessing OSF projects

Let's check out an example project containing materials for an analysis of the 2012 American National Election Survey (ANES). You can access the OSF project in your browser by navigating to its URL: <https://osf.io/jgyxm/>.

Let's load this project into R with `osfr::osf_retrieve_node()`:


```r
anes_project <- osf_retrieve_node("https://osf.io/jgyxm")
anes_project
#> # A tibble: 1 x 3
#>   name                                id    meta            
#>   <chr>                               <chr> <list>          
#> 1 Political identification and gender jgyxm <named list [3]>
```

This returns an `osf_tbl` object, which is the `data.frame`-like class *osfr* uses to represent items retrieved from OSF. You can now use `anes_project` to perform a variety of project related tasks by passing it to different osfr functions.

### Downloading files

Let's list all of the files that have been uploaded to the project:


```r
anes_files <- osf_ls_files(anes_project)
anes_files
#> # A tibble: 5 x 3
#>   name                 id                       meta            
#>   <chr>                <chr>                    <list>          
#> 1 cleaning.R           5e20d22bedceab002d82e0f1 <named list [3]>
#> 2 Questionnaire.docx   5e20d22bedceab002b82dc3f <named list [3]>
#> 3 raw_data.csv         5e20d22c675e0e00096b4de8 <named list [3]>
#> 4 Data Dictionary.docx 5e20d22c675e0e000e6b4b18 <named list [3]>
#> 5 analyses.R           5e20d22c675e0e000a6b4bd3 <named list [3]>
```

This returns another `osf_tbl` but this one contains 5 rows;  one for each of the project *files* stored on OSF. A nice feature of OSF is it provides rendered views for a wide variety of file formats, so it's not necessary to actually download and open a file if you just want to quickly examine it. Let's open the Word Document containing the project's data dictionary by extracting the relevant row from `anes_tbl` and passing it to `osf_open()`:


```r
osf_open(anes_files[4, ])
```

Because `osf_tbl`s are just specialized `data.frame`s, we could also `subset()` or `dplyr::filter()` to achieve the same result.

*__Note:__ If an `osf_tbl` with multiple entities is passed to an non-vectorized osfr function like `osf_open()`, the default behavior is to use the entity in the first row and warn that all other entities are ignored.*

We can also download local copies of these files by passing `anes_files` to `osf_download()`.


```r
osf_download(anes_files)
#> # A tibble: 5 x 4
#>   name                id                     local_path           meta          
#>   <chr>               <chr>                  <chr>                <list>        
#> 1 cleaning.R          5e20d22bedceab002d82e… ./cleaning.R         <named list […
#> 2 Questionnaire.docx  5e20d22bedceab002b82d… ./Questionnaire.docx <named list […
#> 3 raw_data.csv        5e20d22c675e0e00096b4… ./raw_data.csv       <named list […
#> 4 Data Dictionary.do… 5e20d22c675e0e000e6b4… ./Data Dictionary.d… <named list […
#> 5 analyses.R          5e20d22c675e0e000a6b4… ./analyses.R         <named list […
```

We'll use these files in the next section for creating a new project.

### Pipes

As you've likely noticed, `osf_tbl` objects are central to osfr's functionality. Indeed, nearly all of its functions both expect an `osf_tbl` as input and return an `osf_tbl` as output. As such, osfr functions can be chained together using the [pipe operator][magrittr] (`%>%`), allowing for the creation of pipelines to automate OSF-based tasks.

Here is a short example that consolidates all of the steps we've performed so far:


```r
osf_retrieve_node("jgyxm") %>%
  osf_ls_files() %>%
  osf_download()
```

## Project management

Now let's see how to use osfr to create and manage your own projects. The goal for this section is to create your own version of the *Political Identification and Gender* project but with a better organizational structure. To follow along with this section you'll need to authenticate osfr using a personal access token (PAT). See the `?osf_auth()` function documentation or the `auth` vignette for more information.

### Creating a project

First you will need to create a new private project on OSF to store all the files related to the project. Here, we're giving the new project a title (required) and description (optional).


```r
my_project <- osf_create_project(
  title = "Political Identification and Gender: Re-examined",
  description = "A re-analysis of the original study's results."
)
my_project
#> # A tibble: 1 x 3
#>   name                                             id    meta            
#>   <chr>                                            <chr> <list>          
#> 1 Political Identification and Gender: Re-examined kqy2p <named list [3]>
```

The GUID for this new project is `kqy2p`, but yours will be something different. You can check out the project on OSF by opening it's URL (`https://www.osf.io/<GUID>`), or, more conveniently: `osf_open(my_project)`.

### Adding structure with components

A key organizational feature of OSF is the ability to augment a project's structure with sub-projects, which are referred to as *components* on OSF. Like top-level projects, every component is assigned a unique URL and contains its own cloud-based storage bucket. They can also have different privacy settings from the parent project.

We are going to create two nested *components*, one for the raw data and one for the analysis scripts.


```r
data_comp <- osf_create_component(my_project, title = "Raw Data")
script_comp <- osf_create_component(my_project, title = "Analysis Scripts")

# verify the components were created
# osf_open(my_project)
```

If you refresh the OSF project in your browser the *Components* widget should now contain two entries for each of our newly created components.

### Uploading files

Now that our project components are in place we can start to populate them with files. Let's start with the csv file containing our raw data.


```r
data_file  <- osf_upload(my_project, path = "raw_data.csv")
data_file
#> # A tibble: 1 x 3
#>   name         id                       meta            
#>   <chr>        <chr>                    <list>          
#> 1 raw_data.csv 5e3c4af3032a4d00dee71f02 <named list [3]>
```

Oh no! Instead of uploading `raw_data.csv` to the *Raw Data* component, we uploaded it to the parent project instead.

Fear not. We can easily fix this contrived mistake by simply moving the file to its intended location.


```r
data_file <- osf_mv(data_file, to = data_comp)
```

Crisis averted. Now if you open *Raw Data* on OSF (`osf_open(data_comp)`), it should contain the csv file.

Our next step is to upload the R scripts into the *Analysis Scripts* component. Rather than upload each file individually, we'll take advantage of `osf_upload()`'s ability to handle multiple files/directories and use `list.files()` to identify all `.R` files in the working directory:


```r
r_files <- osf_upload(script_comp, path = list.files(pattern = ".R$"))
r_files
#> # A tibble: 3 x 3
#>   name         id                       meta            
#>   <chr>        <chr>                    <list>          
#> 1 analyses.R   5e3c4af6f1369e01158aaf3b <named list [3]>
#> 2 cleaning.R   5e3c4af9032a4d00e8e70915 <named list [3]>
#> 3 precompile.R 5e3c4afaf1369e010e8acaf7 <named list [3]>
```

### Putting it all together

Finally, let's repeat the process for the 2 `.docx` file containing the survey and accompanying data dictionary. This time we'll use a more succinct approach that leverages pipes to create and populate the component in one block:


```r
my_project %>%
  osf_create_component("Research Materials") %>%
  osf_upload(path = list.files(pattern = "\\.docx$"))
#> # A tibble: 2 x 3
#>   name                 id                       meta            
#>   <chr>                <chr>                    <list>          
#> 1 Data Dictionary.docx 5e3c4afef1369e01108ad708 <named list [3]>
#> 2 Questionnaire.docx   5e3c4b00032a4d00ece6faab <named list [3]>
```

We can verify the project is now structured the way we wanted by listing the components we have under the main project.


```r
osf_ls_nodes(my_project)
#> # A tibble: 3 x 3
#>   name               id    meta            
#>   <chr>              <chr> <list>          
#> 1 Research Materials hvmpr <named list [3]>
#> 2 Analysis Scripts   85tck <named list [3]>
#> 3 Raw Data           mzq49 <named list [3]>
```

which gives us an `osf_tbl` with one row for each of the project's components.

### Updating files

OSF provides automatic and unlimited file versioning. Let's see how this works with osfr. Make a small edit to your local copy of `cleaning.R` and save. Now, if we attempt to upload this new version to the *Analysis Scripts* component, osfr will throw a conflict error:


```r
osf_upload(script_comp, path = "cleaning.R")
```

```
Error: Can't upload file 'cleaning.R'.
  * A file with the same name already exists at the destination.
  * Use the `conflicts` argument to avoid this error in the future.
```

 As the error indicates, we need to use the `conflicts` argument to instruct `osf_upload()` how to handle the conflict. In this case, we want to overwrite the original copy with our new version:

 
 ```r
 osf_upload(script_comp, path = "cleaning.R", conflicts = "overwrite")
 ```

Learn more about file versioning on OSF [here][osf-versioning].

### Sharing

Remember, new OSF projects are *always* private by default. You can change this by opening the project's settings page on OSF and making it public. See the following guides more information about OSF permissions and how to optionally generate a DOI so other can cite your project.

* [Privacy Settings][osf-privacy]
* [Sharing Projects][osf-sharing]
* [Generate DOIs][osf-doi]

## A few details about files on OSF

On OSF, files can exist in projects, components, and/or directories. Files can be stored on *OSF's Storage* or in another service that is connected to an OSF project (e.g. GitHub, Dropbox, or Google Drive). However, `osfr` currently only supports interacting with files on OSF Storage.

We can download files from any public or private node that we have access to and can identify files to download in two different ways:

1. If we know where the file is located, but don't remember its GUID, you can use the `osf_ls_files` function to filter by filename within a specified node and then pipe the results to `osf_download()`.

    
    ```r
    anes_project %>%
      osf_ls_files(pattern = ) %>%
      osf_download(conflicts = "overwrite")
    ```

2. For a public file that was referenced in a published article, you may already have the GUID, and so can retrieve the file directly before downloading it. For example, let's download Daniel Laken's helpful spreadsheet for calculating effect sizes (available from <https://osf.io/vbdah/>).

    
    ```r
    osf_retrieve_file("vbdah") %>%
      osf_download(excel_file)
    ```



## Additional resources

For more information on OSF and `osfr` check out:

* [OSF][osf]
* [OSF API Documentation][osf-api]
* [OSF Support](https://osf.io/support/)
* [osfr GitHub Repository](https://github.com/ropensci/osfr)

<!-- links -->
[osf]: https://osf.io
[cos]: https://cos.io
[osf-api]: https://developer.osf.io
[magrittr]: https://magrittr.tidyverse.org
[tibble]: https://tibble.tidyverse.org

[osf-projects]: https://help.osf.io/hc/en-us/categories/360001495973-Creating-and-Managing-Projects
[osf-components]: https://help.osf.io/hc/en-us/articles/360019737614-Create-Components
[osf-privacy]: https://help.osf.io/hc/en-us/articles/360018981414-Control-Your-Privacy-Settings
[osf-permissions]: https://help.osf.io/hc/en-us/articles/360019737774-Edit-Contributor-Permissions
[osf-doi]: https://help.osf.io/hc/en-us/articles/360019931013-Create-DOIs
[osf-versioning]: https://help.osf.io/hc/en-us/articles/360019738694-File-Revisions-and-Version-Control
[osf-sharing]: https://help.osf.io/hc/en-us/categories/360001530614-Sharing-Projects
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_tbl.R
\name{osf_tbl}
\alias{osf_tbl}
\alias{osf_tbl_node}
\alias{osf_tbl_file}
\alias{osf_tbl_user}
\title{OSF Tibbles}
\description{
Items retrieved from OSF are represented as \code{osf_tbl} objects, specialized
data frames based on the \link[tibble:tibble-package]{tibble} class. See below
for additional details.
}
\details{
Each row of an \code{osf_tbl} represents a single OSF entity. This could be a
user, project, component, directory, or file. An \code{osf_tbl} must include
the following 3 columns:
\enumerate{
\item \code{name}: indicates the name or title of the entity.
\item \code{id}: the unique identifier assigned by OSF.
\item \code{meta}: a list-column that stores the processed response returned by OSF's
API. See the \emph{Meta column} section below for more information.
}
}
\section{Subclasses}{


\code{osf_tbl} is the parent class of 3 subclasses that are used to represent
each of OSF's main entities:
\enumerate{
\item \code{osf_tbl_user} for users.
\item \code{osf_tbl_file} for files and directories.
\item \code{osf_tbl_node} for projects and components.
}
}

\section{OSF nodes}{


Projects and components are both implemented as \emph{nodes} on OSF. The only
distinction between the two is that a project is a top-level node, and a
component must have a parent node (i.e., must be a sub-component of another
project or component). Because projects and components are functionally
identical, osfr uses the same \code{\link{osf_tbl_node}} class to represent both.
}

\section{Meta column}{


The \code{meta} column contains all of the information returned from OSF's API for
a single entity, structured as a named list with 3 elements:
\enumerate{
\item \code{attributes} contains metadata about the entity (e.g., names,
descriptions, tags, etc).
\item \code{links} contains urls to API endpoints with alternative representations of
the entity or actions that may be performed on the entity.
\item \code{relationships} contains URLs to other entities with relationships to the
entity (e.g., collaborators attached to a project).
}

This information is critical for \code{osfr}'s internal functions and should not
be altered by users. For even more information about these elements, see
\href{https://developer.osf.io/#tag/Entities-and-Entity-Collections}{OSF's API documentation}.
}

\section{Acknowledgments}{


Our implementation of the \code{osf_tbl} class is based on \code{dribble} objects from
the \href{https://googledrive.tidyverse.org}{googledrive} package.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_ls_files.R
\name{osf_ls_files}
\alias{osf_ls_files}
\title{List files and directories on OSF}
\usage{
osf_ls_files(
  x,
  path = NULL,
  type = "any",
  pattern = NULL,
  n_max = 10,
  verbose = FALSE
)
}
\arguments{
\item{x}{One of the following:
\itemize{
\item An \code{\link{osf_tbl_node}} with a single project or component.
\item An \code{\link{osf_tbl_file}} with a single directory.
}}

\item{path}{List files within the specified subdirectory path.}

\item{type}{Filter query by type. Set to \code{"file"} to list only files, or
\code{"folder"}to list only folders}

\item{pattern}{Character string used to filter for results that contain the
substring \code{"pattern"} in their name. \emph{Note:} this is a fixed, case-insensitive
search.}

\item{n_max}{Maximum number of results to return from OSF (default is 10).
Set to \code{Inf} to return \emph{all} results.}

\item{verbose}{Logical, indicating whether to print informative messages
about interactions with the OSF API (default \code{FALSE}).}
}
\value{
An \code{\link{osf_tbl_file}} with one row for each file or directory, ordered
by modification time.
}
\description{
List the files and directories in the top-level of an OSF project, component,
or directory. Specify a \code{path} to list the contents of a particular
subdirectory.
}
\examples{
\dontrun{
# Retrieve the Psychology Reproducibility Project from OSF
psych_rp <- osf_retrieve_node("ezum7")

# List all files and directories
osf_ls_files(psych_rp)

# ...only the directories
osf_ls_files(psych_rp, type = "folder")

# ...only PDF files
osf_ls_files(psych_rp, type = "file", pattern = "pdf")

# List the contents of the first directory
osf_ls_files(psych_rp, path = "RPP_SI_Figures")
}
}
\seealso{
\code{\link[=osf_ls_nodes]{osf_ls_nodes()}} to generate a list of projects and components.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_create.R
\name{osf_create}
\alias{osf_create}
\alias{osf_create_project}
\alias{osf_create_component}
\title{Create a new project or component on OSF}
\usage{
osf_create_project(
  title,
  description = NULL,
  public = FALSE,
  category = "project"
)

osf_create_component(
  x,
  title,
  description = NULL,
  public = FALSE,
  category = NULL
)
}
\arguments{
\item{title, description}{Set a title (required) and, optionally, a
description.}

\item{public}{Logical, should it be publicly available (\code{TRUE}) or private
(\code{FALSE}, the default)?}

\item{category}{Character string, specify a category to change the icon
displayed on OSF. The defaults are \code{"project"} for projects and
\code{"uncategorized"} for components. The specified category can be easily
changed later on OSF. Valid category options include:
\itemize{
\item analysis
\item communication
\item data
\item hypothesis
\item instrumentation
\item methods and measures
\item procedure
\item project
\item software
\item other
}}

\item{x}{An \code{\link{osf_tbl_node}} with a single OSF project or component that will
serve as the new sub-component's parent node.}
}
\value{
An \code{\link{osf_tbl_node}} containing the new project or component.
}
\description{
Use \code{osf_create_project()} to create a new top-level project on OSF. A nested
component can be created by providing an \code{\link{osf_tbl_node}} containing an
existing project or component to \code{osf_create_component()}'s \code{x} argument.
}
\section{OSF nodes}{


Projects and components are both implemented as \emph{nodes} on OSF. The only
distinction between the two is that a project is a top-level node, and a
component must have a parent node (i.e., must be a sub-component of another
project or component). Because projects and components are functionally
identical, osfr uses the same \code{\link{osf_tbl_node}} class to represent both.
}

\examples{
\dontrun{
# create a new public project
project <- osf_create_project(title = "Private OSF Project", public = TRUE)

# add a private component to the new project
component <- osf_create_component(project, title = "Project Data")
}
}
\references{
\enumerate{
\item OSF Guides: Create a Project.
\url{https://help.osf.io/hc/en-us/articles/360019737594-Create-a-Project}.
\item OSF Guides: Create a Component.
\url{https://help.osf.io/hc/en-us/articles/360019737614-Create-Components}.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_retrieve.R
\name{osf_retrieve}
\alias{osf_retrieve}
\alias{osf_retrieve_user}
\alias{osf_retrieve_node}
\alias{osf_retrieve_file}
\title{Retrieve an entity from OSF}
\usage{
osf_retrieve_user(id)

osf_retrieve_node(id)

osf_retrieve_file(id)
}
\arguments{
\item{id}{An OSF identifier corresponding to an OSF user, project, component,
or file. Set \code{id = "me"} to retrieve your own OSF profile.}
}
\value{
An \code{\link{osf_tbl_user}}, \code{\link{osf_tbl_node}}, or \code{\link{osf_tbl_file}} containing
the corresponding OSF entity.
}
\description{
Create an \code{\link{osf_tbl}} representation of an existing OSF project, component,
file, or user based on the associated unique identifier. Usually this is a
5-character global unique identifier (GUID) but for files or directories, it
could also be an 11-character Waterbutler ID. See below for details.
}
\section{OSF identifiers}{
 A 5-character GUID is assigned to every user,
project, component, and file on OSF and forms the basis for the service's
URL scheme. For example the GUID for a project accessible at
\url{https://osf.io/ezum7} is simply \code{ezum7}. You can learn more about GUIDs
\href{https://help.osf.io/hc/en-us/articles/360019737894-FAQs}{here}.

An important detail is that files and directories are handled internally on
OSF by another serviced called \href{http://www.waterbutler.io/}{Waterbutler},
which uses 11-character identifiers. Although Waterbutler IDs are largely
hidden from users on \url{https://osf.io}, they represent the primary method for
identifying files/directories by the OSF API. In fact, files do not receive a
GUID until it is viewed directly on \url{https://osf.io} and directories never
receive a GUID. Therefore, osfr relies on Waterbutler IDs for files and
directories, and always includes them (rather than GUIDs) in \code{\link{osf_tbl_file}}
objects.
}

\section{Retrieving OSF objects}{

To begin using osfr to interact with resources on OSF you must use one of the
following \emph{retrieve} functions to create an \code{\link{osf_tbl}} that contains
the entity of interest. Note the functions are entity-type specific, use:
\itemize{
\item \code{osf_retrieve_node()} to retrieve a project or component
\item \code{osf_retrieve_file()} to retrieve a file or directory
\item \code{osf_retrieve_user()} to retrieve a user
}
}

\section{A note on 3rd-party storage providers}{

While OSF supports integration with a variety of 3rd-party cloud storage
providers, osfr can currently only access files stored on the default OSF
storage service. Support for additional storage providers is planned for a
future release.
}

\examples{
\dontrun{
 # retrieve your own OSF user profile (must be authenticated, ?osf_auth)
 osf_retrieve_user("me")

# retrieve the Psychology Reproducibility Project (P:RP, osf.io/ezum7)
osf_retrieve_node("ezum7")

# get the first figure from the P:RP
osf_retrieve_file("https://osf.io/7js8c")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_rm.R
\name{osf_rm}
\alias{osf_rm}
\title{Delete an entity from OSF}
\usage{
osf_rm(x, recurse = FALSE, verbose = FALSE, check = TRUE)
}
\arguments{
\item{x}{One of the following:
\itemize{
\item An \code{\link{osf_tbl_node}} with a single OSF project or component.
\item An \code{\link{osf_tbl_file}} containing a single directory or file.
}}

\item{recurse}{Remove all sub-components before deleting the top-level
entity. This only applies when deleting projects or components.}

\item{verbose}{Logical, indicating whether to print informative messages
about interactions with the OSF API (default \code{FALSE}).}

\item{check}{If \code{FALSE} deletion will proceed without opening the item or
requesting verification---this effectively removes your safety net.}
}
\value{
Invisibly returns \code{TRUE} if deletion was successful.
}
\description{
Use \code{osf_rm()} to \strong{permanently} delete a project, component, file or
directory from OSF, including any uploaded files, wiki content, or comments
contained therein. Because this process is \strong{irreversible}, osfr will first
open the item in your web browser so you can verify what is about to be
deleted before proceeding.

If the project or component targeted for deletion contains sub-components,
those must be deleted first. Setting \code{recurse = TRUE} will attempt to
remove the hierarchy of sub-components before deleting the top-level entity.

\emph{Note: This functionality is limited to contributors with admin-level
permissions.}
}
\examples{
\dontrun{
project <- osf_create_project("My Short-Lived Project")
osf_rm(project)
}

}
\seealso{
Other OSF file operations: 
\code{\link{osf_cp}()},
\code{\link{osf_mkdir}()},
\code{\link{osf_mv}()}
}
\concept{OSF file operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_open.R
\name{osf_open}
\alias{osf_open}
\title{Open on OSF}
\usage{
osf_open(x)
}
\arguments{
\item{x}{one of the following:
\itemize{
\item an OSF URL, or a generic string containing a GUID or Waterbutler ID.
\item an \code{\link{osf_tbl_node}} with a single project or component.
\item an \code{\link{osf_tbl_file}} with a single file or directory.
\item an \code{\link{osf_tbl_user}} with a single OSF user.
}}
}
\description{
View a project, component, file, or user profile on OSF with your default
web browser.
}
\examples{
\dontrun{
# Navigate to a project based on its GUID
osf_open("e81xl")

# You can also provide an osf_tbl subclass
crp_file <- osf_retrieve_file("ucpye")
osf_open(crp_file)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osfr-package.R
\docType{package}
\name{osfr-package}
\alias{osfr}
\alias{osfr-package}
\title{osfr: R interface to OSF}
\description{
osfr provides a suite of functions for interacting with the Open Science
Framework (OSF; \url{https://www.osf.io}).
}
\section{What is OSF?}{


OSF is a free and open source project management repository designed to
support researchers across their entire project lifecycle. The service
includes free cloud storage and file version history, providing a
centralized location for all your research materials that can be kept
private, shared with select collaborators, or made publicly available with
citable DOIs.

Most work on OSF is organized around \emph{\strong{projects}}. Projects can contain
\emph{files}, groups of files in \emph{directories}, and/or files in sub-projects
called \emph{\strong{components}}. Note there is no storage limit on the size of
projects but individual files must be < 5Gb.
}

\section{Resources}{

\itemize{
\item To learn more about OSF check out the helpful series of guides published by
the Center for Open Science: \url{https://help.osf.io}
\item See the vignette for an overview of osfr's features:
\code{vignette("getting_started", package = "osfr")}
}
}

\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/osfr}
  \item \url{https://github.com/ropensci/osfr}
  \item Report bugs at \url{https://github.com/ropensci/osfr/issues}
}

}
\author{
\strong{Maintainer}: Aaron Wolen \email{aaron@wolen.com} (\href{https://orcid.org/0000-0003-2542-2202}{ORCID})

Authors:
\itemize{
  \item Chris Hartgerink \email{chjh@protonmail.com} (\href{https://orcid.org/0000-0003-1050-6809}{ORCID})
}

Other contributors:
\itemize{
  \item Timothy York \email{timothypyork@gmail.com} (\href{https://orcid.org/0000-0003-4068-4286}{ORCID}) [contributor]
  \item Ryan Hafen \email{rhafen@purdue.edu} [contributor]
  \item Brian Richards \email{brian.g.richards@gmail.com} [contributor]
  \item Courtney Soderberg \email{courtney@cos.io} (\href{https://orcid.org/0000-0003-1227-7042}{ORCID}) [contributor]
  \item Carl Boettiger \email{cboettig@gmail.com} (\href{https://orcid.org/0000-0002-1642-628X}{ORCID}) [reviewer]
  \item Heidi Seibold \email{heidi.seibold@stat.uni-muenchen.de} [reviewer]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_auth.R
\name{osf_auth}
\alias{osf_auth}
\title{Authenticate osfr with a personal access token}
\usage{
osf_auth(token = NULL)
}
\arguments{
\item{token}{OSF personal access token.}
}
\value{
Invisibly returns your OSF PAT along with a message indicating it was
registered.
}
\description{
Authorize osfr to interact with your OSF data and OSF account by passing a
personal access token (PAT) to \code{osf_auth()}. If no \code{token} is provided,
\code{osf_auth()} will attempt to obtain a PAT from the \code{OSF_PAT} environment
variable. However, since osfr checks for the presence of \code{OSF_PAT} on
start-up, this is only necessary if the variable was created or redefined in
the middle of a session. See below for additional details and instructions
for generating and utilizing your PAT.
}
\details{
Out of the box osfr can only access publicly available projects, components,
and files on OSF. In order for osfr to view and manage your private resources
you must provide a Personal Access Token (PAT). The following instructions
will walk you through the process of generating a PAT and using it to
authenticate osfr.
}
\section{Creating an OSF personal access token}{

\itemize{
\item Navigate to \url{https://osf.io/settings/tokens/}
\item Click the \emph{New token} button and provide a descriptive name
\item Select the scopes (i.e., permissions) you'd like to grant osfr
\item Click the \emph{Create} button to generate your PAT
\item If successful, your 70 character token will be displayed along with several
important warnings you should definitely read over carefully
\item You read those warnings, right?
\item Copy your token and keep it in a safe place
}
}

\section{Using your PAT with osfr}{


There are two possible approaches for authenticating osfr with your PAT.

The simpler approach is to call the \code{osf_auth()} function at the start of
every new R session and manually paste in your token. Note that your PAT
should be treated like a password and, as such, should not be hardcoded into
your script.

A better approach is to store your PAT as an environment variable called
\code{OSF_PAT}. Doing so will allow osfr to detect and utilize the token
automatically without any need to manually call \code{osf_auth()}. One way to
accomplish this is by creating an \code{.Renviron} file in your home or working
directory that defines the \code{OSF_PAT} variable. For example:\preformatted{OSF_PAT=bdEEFMCuBtaBoSK11YzyjOdjUjKtWIj2FWxHl6kTBRax7uaeyBghALumTO1kT8RA
}

For new users we suggest adding the \code{.Renviron} to your home directory so it
is automatically applied to all your projects. To verify this was done
correctly, restart R and run \code{Sys.getenv("OSF_PAT")}, which should return
your PAT.
}

\examples{
\dontrun{
# manually authenticate with a PAT
osf_auth("bdEEFMCuBtaBoSK11YzyjOdjUjKtWIj2FWxHl6kTBRax7uaeyBghALumTO1kT8RA")
}

}
\references{
\enumerate{
\item Colin Gillespie and Robin Lovelace (2017). \emph{Efficient R programming}.
O'Reilly Press. \url{https://csgillespie.github.io/efficientR}.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_cp.R
\name{osf_cp}
\alias{osf_cp}
\title{Copy a file or directory}
\usage{
osf_cp(x, to, overwrite = FALSE, verbose = FALSE)
}
\arguments{
\item{x}{An \code{\link{osf_tbl_file}} containing a single file or directory.}

\item{to}{Destination where the file or directory will be copied.
This can be one of the following:
\itemize{
\item An \code{\link{osf_tbl_node}} with a single project or component.
\item An \code{\link{osf_tbl_file}} with a single directory.
}}

\item{overwrite}{Logical, if a file or directory with the same name already
exists at the destination should it be replaced with \code{x}?}

\item{verbose}{Logical, indicating whether to print informative messages
about interactions with the OSF API (default \code{FALSE}).}
}
\value{
An \code{\link{osf_tbl_file}} containing the updated OSF file.
}
\description{
Use \code{osf_cp()} to make a copy of a file or directory in a new location.
}
\details{
Note that a file (or directory) cannot be moved or copied onto itself, even
if \code{overwrite = TRUE}.
}
\examples{
\dontrun{
project <- osf_create_project("Flower Data")

write.csv(iris, file = "iris.csv")
data_file <- osf_upload(project,"iris.csv")

# Create a new directory to copy our file to
data_dir <- osf_mkdir(project, "data")

# Copy the file to our data directory
data_file <- osf_cp(data_file, to = data_dir)

# Copy directory to new component
data_comp <- osf_create_component(project, title = "data", category = "data")
data_dir \%>\%
  osf_cp(to = data_comp) \%>\%
  osf_open()
}

}
\seealso{
Other OSF file operations: 
\code{\link{osf_mkdir}()},
\code{\link{osf_mv}()},
\code{\link{osf_rm}()}
}
\concept{OSF file operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
Pipe operator
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_mv.R
\name{osf_mv}
\alias{osf_mv}
\title{Move a file or directory}
\usage{
osf_mv(x, to, overwrite = FALSE, verbose = FALSE)
}
\arguments{
\item{x}{An \code{\link{osf_tbl_file}} containing a single file or directory.}

\item{to}{Destination where the file or directory will be moved. This
can be one of the following:
\itemize{
\item An \code{\link{osf_tbl_node}} with a single project or component.
\item An \code{\link{osf_tbl_file}} with a single directory.
}}

\item{overwrite}{Logical, if a file or directory with the same name already
exists at the destination should it be replaced with \code{x}?}

\item{verbose}{Logical, indicating whether to print informative messages
about interactions with the OSF API (default \code{FALSE}).}
}
\value{
An \code{\link{osf_tbl_file}} containing the updated OSF file.
}
\description{
Use \code{osf_mv()} to move a file or directory to a new project, component, or
subdirectory.
}
\details{
Note that a file (or directory) cannot be moved or copied onto itself, even
if \code{overwrite = TRUE}.
}
\examples{
\dontrun{
# Create an example file to upload to our example project
project <- osf_create_project("Flower Data")

write.csv(iris, file = "iris.csv")
data_file <- osf_upload(project,"iris.csv")

# Create a new directory to move our file to
data_dir <- osf_mkdir(project, "data")

# Move the file to our data directory
data_file <- osf_mv(data_file, to = data_dir)

# Move our data directory to a new component
data_comp <- osf_create_component(project, title = "data", category = "data")
data_dir \%>\%
  osf_mv(to = data_comp) \%>\%
  osf_open()
}

}
\seealso{
Other OSF file operations: 
\code{\link{osf_cp}()},
\code{\link{osf_mkdir}()},
\code{\link{osf_rm}()}
}
\concept{OSF file operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_refresh.R
\name{osf_refresh}
\alias{osf_refresh}
\title{Refresh an OSF entity}
\usage{
osf_refresh(x)
}
\arguments{
\item{x}{an \code{\link{osf_tbl}}.}
}
\description{
Use \code{osf_refresh()} to update one or more entities in an \code{\link[=osf_tbl]{osf_tbl()}} with
the latest information from OSF.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_ls_nodes.R
\name{osf_ls_nodes}
\alias{osf_ls_nodes}
\title{List projects or components on OSF}
\usage{
osf_ls_nodes(x, pattern = NULL, n_max = 10, verbose = FALSE)
}
\arguments{
\item{x}{one of the following:
\itemize{
\item An \code{\link{osf_tbl_node}} with a single project or component.
\item An \code{\link{osf_tbl_user}} with a single OSF user.
}}

\item{pattern}{Character string used to filter for results that contain the
substring \code{"pattern"} in their name. \emph{Note:} this is a fixed, case-insensitive
search.}

\item{n_max}{Maximum number of results to return from OSF (default is 10).
Set to \code{Inf} to return \emph{all} results.}

\item{verbose}{Logical, indicating whether to print informative messages
about interactions with the OSF API (default \code{FALSE}).}
}
\value{
An \code{\link{osf_tbl_node}} with one row for each OSF project or component,
ordered by modification time.
}
\description{
List the projects or components associated with a user or contained in the
top-level of another OSF project or component.
}
\examples{
\dontrun{
# List your recent projects and components
osf_retrieve_user("me") \%>\%
  osf_ls_nodes()

# List the first 10 components in the #ScanAllFish project
fish_ctscans <- osf_retrieve_node("ecmz4")
osf_ls_nodes(fish_ctscans)

# Now just the components with scans of species from the Sphyrna genus
osf_ls_nodes(fish_ctscans, pattern = "Sphyrna")
}
}
\seealso{
\code{\link[=osf_ls_files]{osf_ls_files()}} to generate a list of files and files.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_id.R
\name{as_id}
\alias{as_id}
\title{Extract OSF identifiers}
\usage{
as_id(x)
}
\arguments{
\item{x}{An \code{osf_tbl}, OSF URL, or a generic string containing a GUID or
Waterbutler ID.}
}
\value{
A character vector with class \code{osf_id}.
}
\description{
Extract OSF GUIDs and Waterbutler IDs from various types of inputs. Valid
looking IDs are returned as \code{osf_id} objects.
}
\section{Identifier types}{

There are 2 types of identifiers you'll encounter on OSF. The first is the
globally unique identifier, or GUID, that OSF assigns to every entity. A
valid OSF GUID consists of 5 alphanumeric characters. The second type of
identifier is specific to files stored on OSF. All file operations on OSF are
handled via Waterbutler. A valid Waterbutler ID consists of 24 alphanumeric
characters.
}

\examples{
\dontrun{
# extract a GUID from an OSF URL
proj_id <- as_id("https://osf.io/7zqxp/")

# extract waterbutler IDs from an `osf_tbl_file`` with multiple files
osf_retrieve_node(proj_id) \%>\%
  osf_ls_files() \%>\%
  as_id()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_mkdir.R
\name{osf_mkdir}
\alias{osf_mkdir}
\title{Create directories on OSF}
\usage{
osf_mkdir(x, path, verbose = FALSE)
}
\arguments{
\item{x}{One of the following:
\itemize{
\item An \code{\link{osf_tbl_node}} with a single OSF project or component.
\item An \code{\link{osf_tbl_file}} containing a single directory.
}}

\item{path}{Name of the new directory or a path ending with the new directory.}

\item{verbose}{Logical, indicating whether to print informative messages
about interactions with the OSF API (default \code{FALSE}).}
}
\value{
An \code{\link{osf_tbl_file}} with one row containing the leaf directory
specified in \code{path}.
}
\description{
Use \code{osf_mkdir()} to add new directories to projects, components, or nested
within existing OSF directories. If \code{path} contains multiple directory
levels (e.g., \code{"data/rawdata"}) the intermediate-level directories are
created automatically. If the directory you're attempting to create already
exists on OSF it will be silently ignored and included in the output.
}
\examples{
\dontrun{
proj <- osf_create_project("Directory Example")

# add directory to the top-level of the Directory Example project
data_dir <- osf_mkdir(proj, path = "data")

# add a subdirectory nested within data/
osf_mkdir(data_dir, path = "rawdata")

# recursively create multiple directory levels within data/
osf_mkdir(data_dir, path = "samples/pcr/qc")
}
}
\seealso{
Other OSF file operations: 
\code{\link{osf_cp}()},
\code{\link{osf_mv}()},
\code{\link{osf_rm}()}
}
\concept{OSF file operations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_download.R
\name{osf_download}
\alias{osf_download}
\title{Download files and directories from OSF}
\usage{
osf_download(
  x,
  path = NULL,
  recurse = FALSE,
  conflicts = "error",
  verbose = FALSE,
  progress = FALSE
)
}
\arguments{
\item{x}{An \code{\link{osf_tbl_file}} containing a single file or directory.}

\item{path}{Path pointing to a local directory where the downloaded files
will be saved. Default is to use the current working directory.}

\item{recurse}{Applies only to OSF directories. If \code{TRUE}, a directory is
fully recursed and all nested files and subdirectories are downloaded.
Alternatively, a positive number will determine the number of levels
to recurse.}

\item{conflicts}{This determines what happens when a file with the same name
exists at the specified destination. Can be one of the following:
\itemize{
\item \code{"error"} (the default): throw an error and abort the file transfer operation.
\item \code{"skip"}: skip the conflicting file(s) and continue transferring the
remaining files.
\item \code{"overwrite"}: replace the existing file with the transferred copy.
}}

\item{verbose}{Logical, indicating whether to print informative messages
about interactions with the OSF API (default \code{FALSE}).}

\item{progress}{Logical, if \code{TRUE} progress bars are displayed for each file
transfer. Mainly useful for transferring large files. For tracking lots of
small files, setting \code{verbose = TRUE} is more informative.}
}
\value{
The same \code{\link{osf_tbl_file}} passed to \code{x} with a new column,
\code{"local_path"}, containing paths to the local files.
}
\description{
Files stored on OSF can be downloaded locally by passing an \code{\link{osf_tbl_file}}
that contains the files and folders of interest. Use \code{path} to specify
\emph{where} the files should be downloaded, otherwise they are downloaded to
your working directory by default.
}
\section{Implementation details}{

Directories are always downloaded from OSF as zip files that contain its
entire contents. The logic for handling conflicts and recursion is
implemented locally, acting on these files in a temporary location and
copying them to \code{path} as needed. This creates a \emph{gotcha} if you're
downloading directories with large files and assuming that setting \code{conflicts = "skip"} and/or limiting recursion will reduce the number of files you're
downloading. In such a case, a better strategy would be to use
\code{osf_ls_files()} to list the contents of the directory and pass that output
to \code{osf_download()}.
}

\section{A note about synchronization}{

While \code{osf_download()} and \code{osf_upload()} can be used to conveniently shuttle
files back and forth between OSF and your local machine, it's important to
note that \strong{they are not file synchronization functions}. In contrast to
something like \href{https://rsync.samba.org}{\code{rsync}},
\code{osf_download()}/\code{osf_upload()} do not take into account a file's contents or
modification time. Whether you're uploading or downloading, if \code{overwrite = TRUE}, osfr will overwrite an existing file regardless of whether the
existing file is the more recent copy. You have been warned.
}

\examples{
\dontrun{
# download a single file
analysis_plan <- osf_retrieve_file("2ryha") \%>\%
  osf_download()

# verify the file was downloaded locally
file.exists(analysis_plan$local_path)
}
}
\seealso{
\itemize{
\item \code{\link[=osf_upload]{osf_upload()}} for uploading files to OSF.
\item \code{\link[=osf_ls_files]{osf_ls_files()}} for listing files and directories on OSF.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/osf_upload.R
\name{osf_upload}
\alias{osf_upload}
\title{Upload files to OSF}
\usage{
osf_upload(
  x,
  path,
  recurse = FALSE,
  conflicts = "error",
  progress = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{x}{The upload destination on OSF. Can be one of the following:
\itemize{
\item An \code{\link{osf_tbl_node}} with a single project or component.
\item An \code{\link{osf_tbl_file}} with a single directory.
}}

\item{path}{A character vector of paths pointing to existing
local files and/directories.}

\item{recurse}{If \code{TRUE}, fully recurse directories included in \code{path}. You
can also control the number of levels to recurse by specifying a positive
number.}

\item{conflicts}{This determines what happens when a file with the same name
exists at the specified destination. Can be one of the following:
\itemize{
\item \code{"error"} (the default): throw an error and abort the file transfer operation.
\item \code{"skip"}: skip the conflicting file(s) and continue transferring the
remaining files.
\item \code{"overwrite"}: replace the existing file with the transferred copy.
}}

\item{progress}{Logical, if \code{TRUE} progress bars are displayed for each file
transfer. Mainly useful for transferring large files. For tracking lots of
small files, setting \code{verbose = TRUE} is more informative.}

\item{verbose}{Logical, indicating whether to print informative messages
about interactions with the OSF API (default \code{FALSE}).}
}
\value{
An \code{\link{osf_tbl_file}} containing the uploaded files and directories
that were directly specified in \code{path}.
}
\description{
Upload local files to a project, component, or directory on OSF.
}
\section{File and directory paths}{

The \code{x} argument indicates \emph{where} on OSF the files will be uploaded (\emph{i.e.},
the destination). The \code{path} argument indicates \emph{what} will be uploaded,
which can include a combination of files \emph{and} directories.

When \code{path} points to a local file, the file is uploaded to the \emph{root} of the
specified OSF destination, regardless of where it's on your local machine
(\emph{i.e.}, the intermediate paths are not preserved). For example, the
following would would upload both \code{a.txt} and \code{b.txt} to the root of
\code{my_proj}:\preformatted{osf_upload(my_proj, c("a.txt", "subdir/b.txt"))`
}

When \code{path} points to a local directory, a corresponding directory will be
created at the root of the OSF destination, \code{x}, and any files within the
local directory are uploaded to the new OSF directory. Therefore, we could
maintain the directory structure in the above example by passing \code{b.txt}'s
parent directory to \code{path} instead of the file itself:\preformatted{osf_upload(my_proj, c("a.txt", "subdir2"))
}

Likewise, \code{osf_upload(my_proj, path = ".")} will upload your entire current
working directory to the specified OSF destination.
}

\section{Uploading to subdirectories}{

In order to upload directly to an existing OSF directory you would first need
to retrieve the directory as an \code{\link{osf_tbl_file}}. This can be accomplished by
passing the directory's unique identifier to \code{\link[=osf_retrieve_file]{osf_retrieve_file()}}, or, if
you don't have the ID handy, you can use \code{\link[=osf_ls_files]{osf_ls_files()}} to retrieve the
directory by name.\preformatted{# search for the 'rawdata' subdirectory within top-level 'data' directory
target_dir <- osf_ls_files(my_proj, path = "data", pattern = "rawdata")
# upload 'a.txt' to data/rawdata/ on OSF
osf_upload(target_dir, path = "a.txt")
}
}

\section{A note about synchronization}{

While \code{osf_download()} and \code{osf_upload()} can be used to conveniently shuttle
files back and forth between OSF and your local machine, it's important to
note that \strong{they are not file synchronization functions}. In contrast to
something like \href{https://rsync.samba.org}{\code{rsync}},
\code{osf_download()}/\code{osf_upload()} do not take into account a file's contents or
modification time. Whether you're uploading or downloading, if \code{overwrite = TRUE}, osfr will overwrite an existing file regardless of whether the
existing file is the more recent copy. You have been warned.
}

\examples{
\dontrun{
# Create an example file to upload to our example project
write.csv(iris, file = "iris.csv")
project <- osf_create_project("Flower Data")

# Upload the first version
osf_upload(project,"iris.csv")

# Modify the data file, upload version 2, and view it on OSF
write.csv(subset(iris, Species != "setosa"), file = "iris.csv")
project \%>\%
  osf_upload("iris.csv", overwrite = TRUE) \%>\%
  osf_open()
}

}
\seealso{
\itemize{
\item \code{\link[=osf_download]{osf_download()}} for downloading files and directories from OSF.
\item \code{\link[=osf_ls_files]{osf_ls_files()}} for listing files and directories on OSF.
}
}
