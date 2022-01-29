
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
