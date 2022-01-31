README
================

# DataPackageR

DataPackageR is used to reproducibly process raw data into packaged,
analysis-ready data sets.

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/DataPackageR)](https://CRAN.R-project.org/package=DataPackageR)
[![R-CMD-check](https://github.com/ropensci/DataPackageR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/DataPackageR/actions)
[![Coverage
status](https://codecov.io/gh/ropensci/DataPackageR/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/DataPackageR?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://badges.ropensci.org/230_status.svg)](https://github.com/ropensci/software-review/issues/230)
[![DOI](https://zenodo.org/badge/29267435.svg)](https://doi.org/10.5281/zenodo.1292095)
<!-- badges: end -->

-   [yaml configuration
    guide](https://github.com/ropensci/DataPackageR/blob/main/vignettes/YAML_CONFIG.md)
-   [a more detailed technical
    vignette](https://github.com/ropensci/DataPackageR/blob/main/vignettes/usingDataPackageR.md)

> **Important Note**: [datapack](https://github.com/ropensci/datapack)
> is a *different package* that is used to “create, send and load data
> from common repositories such as DataONE into the R environment.”

> **This package** is for processing raw data into tidy data sets and
> bundling them into R packages.

## What problems does DataPackageR tackle?

You have diverse raw data sets that you need to preprocess and tidy in
order to:

-   Perform data analysis
-   Write a report
-   Publish a paper
-   Share data with colleagues and collaborators
-   Save time in the future when you return to this project but have
    forgotten all about what you did.

### Why package data sets?

**Definition:** A *data package* is a formal R package whose sole
purpose is to contain, access, and / or document data sets.

-   **Reproducibility.**

    As described [elsewhere](https://github.com/ropensci/rrrpkg),
    packaging your data promotes reproducibility. R’s packaging
    infrastructure promotes unit testing, documentation, a reproducible
    build system, and has many other benefits. Coopting it for packaging
    data sets is a natural fit.

-   **Collaboration.**

    A data set packaged in R is easy to distribute and share amongst
    collaborators, and is easy to install and use. All the hard work
    you’ve put into documenting and standardizing the tidy data set
    comes right along with the data package.

-   **Documentation.**

    R’s package system allows us to document data objects. What’s more,
    the `roxygen2` package makes this very easy to do with [markup
    tags](https://r-pkgs.org/data.html). That documentation is the
    equivalent of a data dictionary and can be extremely valuable when
    returning to a project after a period of time.

-   **Convenience.**

    Data pre-processing can be time consuming, depending on the data
    type and raw data sets may be too large to share conveniently in a
    packaged format. Packaging and sharing the small, tidied data saves
    the users computing time and time spent waiting for downloads.

## Challenges.

-   **Package size limits.**

    R packages have a 5MB size limit, at least on CRAN. BioConductor has
    explicit [data
    package](https://www.bioconductor.org/developers/package-guidelines/#package-types)
    types that can be larger and use git LFS for very large files.

    Sharing large volumes of raw data in an R package format is still
    not ideal, and there are public biological data repositories better
    suited for raw data: e.g., [GEO](https://www.ncbi.nlm.nih.gov/geo/),
    [SRA](https://www.ncbi.nlm.nih.gov/sra),
    [ImmPort](https://www.immport.org:443/shared/immport-open/public/home/home),
    [ImmuneSpace](https://immunespace.org/),
    [FlowRepository](https://flowrepository.org/).

    Tools like [datastorr](https://github.com/ropenscilabs/datastorr)
    can help with this and we hope to integrate the into DataPackageR in
    the future.

-   **Manual effort**

    There is still a substantial manual effort to set up the correct
    directory structures for an R data package. This can dissuade many
    individuals, particularly new users who have never built an R
    package, from going this route.

-   **Scale**

    Setting up and building R data packages by hand is a workable
    solution for a small project or a small number of projects, but when
    dealing with many projects each involving many data sets, tools are
    needed to help automate the process.

## DataPackageR

DataPackageR provides a number of benefits when packaging your data.

-   It aims to automate away much of the tedium of packaging data sets
    without getting too much in the way, and keeps your processing
    workflow reproducible.

-   It sets up the necessary package structure and files for a data
    package.

-   It allows you to keep the large, raw data and only ship the packaged
    tidy data, saving space and time consumers of your data set need to
    spend downloading and re-processing it.

-   It maintains a reproducible record (vignettes) of the data
    processing along with the package. Consumers of the data package can
    verify how the processing was done, increasing confidence in your
    data.

-   It automates construction of the documentation and maintains a data
    set version and an md5 fingerprint of each data object in the
    package. If the data changes and the package is rebuilt, the data
    version is automatically updated.

## Similar work

There are a number of tools out there that address similar and
complementary problems:

-   **datastorr** [github
    repo](https://github.com/ropenscilabs/datastorr)

    Simple data retrieval and versioning using GitHub to store data.

    -   Caches downloads and uses github releases to version data.
    -   Deal consistently with translating the file stored online into a
        loaded data object
    -   Access multiple versions of the data at once

    `datastorrr` could be used with DataPackageR to store / access
    remote raw data sets, remotely store / access tidied data that are
    too large to fit in the package itself.

-   **fst** [github repo](https://github.com/fstpackage/fst)

    `fst` provides lightning fast serialization of data frames.

-   **The modern data package**
    [pdf](https://github.com/noamross/2018-04-18-rstats-nyc/blob/master/Noam_Ross_ModernDataPkg_rstatsnyc_2018-04-20.pdf)

    A presentation from @noamross touching on modern tools for open
    science and reproducibility. Discusses `datastorr` and `fst` as well
    as standardized metadata and documentation.

-   **rrrpkg** [github repo](https://github.com/ropensci/rrrpkg)

    A document from ropensci describing using an R package as a research
    compendium. Based on ideas originally introduced by Robert Gentleman
    and Duncan Temple Lang (Gentleman and Lang
    (2004)<!--@Gentleman2004-oj-->)

-   **template** [github repo](https://github.com/ropensci/rrrpkg)

    An R package template for data packages.

See the [publication](#publication) for further discussion.

## Installation

You can install the latest version of DataPackageR from
[github](https://github.com/ropensci/DataPackageR) with:

``` r
library(devtools)
devtools::install_github("ropensci/DataPackageR")
```

## Blog Post - building packages interactively.

See this [rOpenSci blog
post](https://ropensci.org/blog/2018/09/18/datapackager/) on how to
build data packages interactively using DataPackageR. This uses several
new interfaces: `use_data_object()`, `use_processing_script()` and
`use_raw_dataset()` to build up a data package, rather than assuming the
user has all the code and data ready to go for `datapackage_skeleton()`.

## Example (assuming all code and data are available)

``` r
library(DataPackageR)

# Let's reproducibly package up
# the cars in the mtcars dataset
# with speed > 20.
# Our dataset will be called cars_over_20.
# There are three steps:

# 1. Get the code file that turns the raw data
# into our packaged and processed analysis-ready dataset.
# This is in a file called subsetCars.Rmd located in exdata/tests of the DataPackageR package.
# For your own projects you would write your own Rmd processing file.
processing_code <- system.file(
  "extdata", "tests", "subsetCars.Rmd", package = "DataPackageR"
)

# 2. Create the package framework.
# We pass in the Rmd file in the `processing_code` variable and the names of the data objects it creates (called "cars_over_20")
# The new package is called "mtcars20"
datapackage_skeleton(
  "mtcars20", force = TRUE, 
  code_files = processing_code, 
  r_object_names = "cars_over_20", 
  path = tempdir()) 

# 3. Run the preprocessing code to build the cars_over_20 data set 
# and reproducibly enclose it in the mtcars20 package.
# packageName is the full path to the package source directory created at step 2.
# You'll be prompted for a text description (one line) of the changes you're making.
# These will be added to the NEWS.md file along with the DataVersion in the package source directory.
# If the build is run in non-interactive mode, the description will read
# "Package built in non-interactive mode". You may update it later.
dir.create(file.path(tempdir(),"lib"))
package_build(packageName = file.path(tempdir(),"mtcars20"), install = TRUE, lib = file.path(tempdir(),"lib"))
#> Warning: package 'mtcars20' is in use and will not be installed

# Update the autogenerated roxygen documentation in data-raw/documentation.R. 
# edit(file.path(tempdir(),"mtcars20","R","mtcars20.R"))

# 4. Rebuild the documentation.
document(file.path(tempdir(),"mtcars20"), install = TRUE, lib = file.path(tempdir(),"lib"))
#> Warning: package 'mtcars20' is in use and will not be installed

# Let's use the package we just created.
install.packages(file.path(tempdir(),"mtcars20_1.0.tar.gz"), type = "source", repos = NULL)
#> Warning: package 'mtcars20' is in use and will not be installed
library(mtcars20)
data("cars_over_20") # load the data
cars_over_20  # Now we can use it.
?cars_over_20 # See the documentation you wrote in data-raw/documentation.R.

# We have our dataset!
# Since we preprocessed it,
# it is clean and under the 5 MB limit for data in packages.
cars_over_20

# We can easily check the version of the data
data_version("mtcars20")

# You can use an assert to check the data version in  reports and
# analyses that use the packaged data.
assert_data_version(data_package_name = "mtcars20",
                    version_string = "0.1.0",
                    acceptable = "equal")
```

### Reading external data from within R / Rmd processing scripts.

When creating a data package, your processing scripts will need to read
your raw data sets in order to process them. These data sets can be
stored in `inst/extdata` of the data package source tree, or elsewhere
outside the package source tree. In order to have portable and
reproducible code, you should not use absolute paths to the raw data.
Instead, `DataPackageR` provides several APIs to access the data package
project root directory, the `inst/extdata` subdirectory, and the `data`
subdirectory.

``` r
# This returns the datapackage source 
# root directory. 
# In an R or Rmd processing script this can be used to build a path to a directory that is exteral to the package, for 
# example if we are dealing with very large data sets where data cannot be packaged.
DataPackageR::project_path()

# This returns the   
# inst/extdata directory. 
# Raw data sets that are included in the package should be placed there.
# They can be read from that location, which is returned by: 
DataPackageR::project_extdata_path()

# This returns the path to the datapackage  
# data directory. This can be used to access 
# stored data objects already created and saved in `data` from 
# other processing scripts.
DataPackageR::project_data_path()
```

## Preprint and publication. <a id = "publication"></a>

The publication describing the package, (Finak *et* *al.,*
2018)<!--@10.12688/gatesopenres.12832.2-->, is now available at [Gates
Open Research](https://gatesopenresearch.org/articles/2-31/v1) .

The preprint is on [biorxiv](https://doi.org/10.1101/342907).

## Code of conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/DataPackageR/blob/main/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.

### References

1.  Gentleman, Robert, and Duncan Temple Lang. 2004. “Statistical
    Analyses and Reproducible Research.” Bioconductor Project Working
    Papers, Bioconductor project working papers,. bepress.

2.  Finak G, Mayer B, Fulp W et al. DataPackageR: Reproducible data
    preprocessing, standardization and sharing using R/Bioconductor for
    collaborative data analysis \[version 1; referees: 1 approved with
    reservations\]. Gates Open Res 2018, 2:31 (doi:
    10.12688/gatesopenres.12832.1)

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# DataPackageR 0.15.8
* Fix to datapackager_object_read that was causing a test to break. `get` needs to have `inherits=FALSE`. 
* Other fixes for `usethis` 1.6.0
* Fixes to tests that were failing on CRAN
* In `package_build`, remove `devtools::reload` and put `devtools::unload` and in front of install.packages
* in `document`, remove `devtools::reload` and put in `devtools::unload` and install.packages

# DataPackageR 0.15.7
* Fix test and vignette bugs related to upcoming version of usethis (1.5)


# DataPackageR 0.15.6
* Fix bug in vignette and code that writes to user space during CRAN checks.

# DataPackageR 0.15.4.900
* Fix a bug in update_news.
* Create news files if it doesn't exist.


# DataPackageR 0.15.4
* New CRAN Release

# DataPackageR 0.15.3.9000

## Features and enhancements
* Reduce the console output from logging. (ropensci/DataPackageR/issues/50)
* Create a new logger that logs at different thresholds to console and to file (ropensci/DataPackageR/issues/50)
* Default on build is not to install.
* Hide console output from Rmd render.
* Nicer messages describing data sets that are created (ropensci/DataPackageR/issues/51)
* Write deleted, changed, and added data objects to the NEWS file automatically.
* Add option to overwrite (or not) via use_processing_script. Provide warning.
* Add use_ignore() to ignore files and data sets in .Rbuildignore and .gitignore and added ignore argument to use_raw_dataset().

## Bug fixes
* code argument no longer required for construct_yml_config
* Fix the documentation for datapackager_object_read() and "Migrating old packages".
* Copy over vignettes generated as pdfs into the package inst/doc
* Data objects are incrementally stored during the build process, into the render_root directory specified in the datapackager.yml config file.

# DataPackageR 0.15.3
* conditional tests when pandoc is missing (ropensci/DataPackager/issues/46)
* add use_data_object and use_processing_script (ropensci/DataPackager/issues/44)
* allow datapacakge_skeleton to be called without files or data objects for interactive construction. (ropensci/DataPackager/issues/44)

# DataPackageR 0.15.2
* Add  pandoc to SystemRequirements (ropensci/DataPackager/issues/46)
* Add use_raw_dataset() method (and tests) to add data sets to inst/extdata. interactively. (ropensci/DataPackager/issues/44)

# DataPackageR 0.15.1.9000
* Development version

# DataPackageR 0.15.1
- Fix CRAN notes.

# DataPackageR 0.15.0
- Prepare for CRAN submission.


# DataPackageR 0.14.9

- Moving towards ropensci compliance
- NEWS.md updated with description of changes to data sets when version is bumped (or new package is created).
- Output of "next steps" for user when pakcage is built
- New `document()` function to rebuild docs from `documentation.R` in `data-raw` without rebuilding the whole package.
- Improved package test.
- R scripts processed properly into vignettes.
- Packages installed and loaded after build to make vignettes and data sets accessible in same R session.
- 

# DataPackageR 0.13.6

- Added a NEWS file.
- Cleaned up the examples.
- Snake case for all exported functions.

# DataPackageR 0.13.3

- Added the `render_root` property to the YAML configuration. Specifies where `render()` processing is done, instead of the `data-raw` directory.
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
## Submission 0.15.8
* Fix tests that were failing. This should resolve the failures identified by CRAN
* Change maintainer to Ellis Hughes from Greg Finak
* update filepaths/URLS based on request from CRAN


## Test environments
* local R installation, R 4.0.3
* github actions/ windows-latest - R Release
* github actions/ macOS-latest - R Release
* github actions/ ubuntu-20.04 - R Release
* github actions/ ubuntu-20.04 - R devel
* github actions/ ubuntu-20.04 - R 3.6
* github actions/ ubuntu-20.04 - R 3.5
* rhub/windows-x86_64- R devel
* rhub/ubuntu-gcc - R release
* rhub/fedora-clang - R devel
* win-builder (devel)


## R CMD check results

There were no ERRORs or WARNINGs

One NOTE:

Maintainer: 'Ellis Hughes <ehhughes@scharp.org>'

New submission

Package was archived on CRAN

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2020-04-25 as check problems were not
    corrected in time.

This package was archived earlier due to not resolving check errors. With this
upload, this should resolve any errors.


## Downstream dependencies

The package has no reverse dependencies.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
# CONTRIBUTING #

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
*  We recommend the tidyverse [style guide](https://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the DataPackageR project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
title: "The DataPackageR YAML configuration file."
author: "Greg Finak <gfinak@fredhutch.org>"
date: "2018-10-24"
output: 
  rmarkdown::html_vignette:
    keep_md: TRUE
    toc: yes
  bibliography: bibliography.bib
vignette: >
  %\VignetteIndexEntry{DataPackageR YAML configuration.}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
  %\usepackage{graphicx}
editor_options: 
  chunk_output_type: inline
---

# Configuring and controlling DataPackageR builds.

Data package builds are controlled using the `datapackager.yml` file. 

This file is created in the package source tree when the user creates a package using `datapackage_skeleton()`. 

It is automatically populated with the names of the `code_files` and `data_objects` the passed in to datapackage_skeleton.

## The `datapackager.yml` file.

The structure of a correctly formatted `datapackager.yml` file is shown below:




```
configuration:
  files:
    subsetCars.Rmd:
      enabled: yes
  objects: cars_over_20
  render_root:
    tmp: '450393'
```

## YAML config file properties.

The main section of the file is the `configuration:` section.

It has three properties:

- `files:` 
  
  The files (`R` or `Rmd`) to be processed by DataPackageR. They are processed in the order shown. Users running multi-script workflows with dependencies between the scripts need to ensure the files are processed in the correct order. 
  
  Here `subsetCars.Rmd` is the only file to process. The name is transformed to an absolute path within the package.
  
  Each file itself has just one property:
  
  <!-- - `name:`  -->
     <!-- The name of the file. This is transformed to an absolute path within the package. -->
    
  - `enabled:` 
     A logical `yes`, `no` flag indicating whether the file should be rendered during the build, or whether it should be skipped.
    This is useful for 'turning off' long running processing tasks if they have not changed. Disabling processing of a file will not overwrite existing documentation or data objecs created during previous builds.
  
- `objects:` 
  
  The names of the data objects created by the processing files, to be stored in the package. These names are compared against the objects created in the render environment by each file. They names must match.
- `render_root:` 
  
  The directory where the `Rmd` or `R` files will be rendered. Defaults to a randomly named subdirectory of `tempdir()`. Allows workflows that use multiple scripts and create file system artifacts to function correctly by simply writing to and reading from the working directory.

## Editing the YAML config file.

The structure of the YAML is simple enough to understand but complex enough that it can be a pain to edit by hand.

DataPackageR provides a number of API calls to construct, read, modify, and write the yaml config file.

### API calls

#### `construct_yml_config` 

  Make an r object representing a YAML config file.
  
##### Example
  The YAML config shown above was created by: 

```r
# Note this is done by the datapackage_skeleton. 
# The user doesn't usually need to call 
# construct_yml_config()
yml <- DataPackageR::construct_yml_config(
  code = "subsetCars.Rmd",
  data = "cars_over_20"
  )
```
  
  
#### `yml_find` 

  Read a yaml config file from a package path into an r object.
  
##### Example
  Read the YAML config file from the `mtcars20` example.
  

```r
# returns an r object representation of
# the config file.
mtcars20_config <- yml_find(
  file.path(tempdir(),"mtcars20")
  )
```

#### `yml_list_objects` 

  List the `objects` in a config read by `yml_find`.
  

##### Example


```r
  yml_list_objects(yml)
```

```

cars_over_20
```
  
#### `yml_list_files` 

  List the `files` in a config read by `yml_find`.

##### Example


```r
  yml_list_files(yml)
```

```

subsetCars.Rmd
```
  
#### `yml_disable_compile` 

  Disable compilation of named files in a config read by `yml_find`.

##### Example


```r
yml_disabled <- yml_disable_compile(
    yml,
    filenames = "subsetCars.Rmd")
```

```
configuration:
  files:
    subsetCars.Rmd:
      enabled: no
  objects: cars_over_20
  render_root:
    tmp: '912178'
```

#### `yml_enable_compile` 

  Enable compilation of named files in a config read by `yml_find`.

##### Example


```r
yml_enabled <- yml_enable_compile(
    yml,
    filenames = "subsetCars.Rmd")
```

```
configuration:
  files:
    subsetCars.Rmd:
      enabled: yes
  objects: cars_over_20
  render_root:
    tmp: '912178'
```

#### `yml_add_files` 
  
  Add named files to a config read by `yml_find`.
  
##### Example


```r
yml_twofiles <- yml_add_files(
    yml,
    filenames = "anotherFile.Rmd")
```

```
configuration:
  files:
    subsetCars.Rmd:
      enabled: yes
    anotherFile.Rmd:
      enabled: yes
  objects: cars_over_20
  render_root:
    tmp: '912178'
```

#### `yml_add_objects` 

  Add named objects to a config read by `yml_find`.

##### Example


```r
yml_twoobj <- yml_add_objects(
    yml_twofiles,
    objects = "another_object")
```

```
configuration:
  files:
    subsetCars.Rmd:
      enabled: yes
    anotherFile.Rmd:
      enabled: yes
  objects:
  - cars_over_20
  - another_object
  render_root:
    tmp: '912178'
```

#### `yml_remove_files` 

  Remove named files from a config read by `yml_find`.

##### Example


```r
yml_twoobj <- yml_remove_files(
    yml_twoobj,
    filenames = "anotherFile.Rmd")
```

```
configuration:
  files:
    subsetCars.Rmd:
      enabled: yes
  objects:
  - cars_over_20
  - another_object
  render_root:
    tmp: '912178'
```

#### `yml_remove_objects` 

  Remove named objects from a config read by `yml_find`.

##### Example


```r
yml_oneobj <- yml_remove_objects(
    yml_twoobj,
    objects = "another_object")
```

```
configuration:
  files:
    subsetCars.Rmd:
      enabled: yes
  objects: cars_over_20
  render_root:
    tmp: '912178'
```

#### `yml_write` 

  Write a modified config to its package path.

##### Example


```r
yml_write(yml_oneobj, path = "path_to_package")
```

The `yml_oneobj` read by `yml_find()` carries an attribute
that is the path to the package. The user doesn't need to pass a `path` to `yml_write` if the config has been read by `yml_find`.
---
title: "Using DataPackageR"
author: "Greg Finak <gfinak@fredhutch.org>"
date: "2019-03-11"
output: 
  rmarkdown::html_vignette:
    keep_md: TRUE
    toc: yes
vignette: >
  %\VignetteIndexEntry{A Guide to using DataPackageR}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
  %\usepackage{graphicx}
---



## Purpose

This vignette demonstrates how to use DataPackageR to build a data package. 

DataPackageR aims to simplify data package construction.

It provides mechanisms for reproducibly preprocessing and tidying raw data into into documented, versioned, and packaged analysis-ready data sets. 

Long-running or computationally intensive data processing can be decoupled from the usual `R CMD build` process while maintinaing [data lineage](https://en.wikipedia.org/wiki/Data_lineage).

In this vignette we will subset and package the `mtcars` data set.

## Set up a new data package.

We'll set up a new data package based on `mtcars` example in the [README](https://github.com/RGLab/DataPackageR/blob/master/README.md).
The `datapackage_skeleton()` API is used to set up a new package. 
The user needs to provide:

- R or Rmd code files that do data processing.
- A list of R object names created by those code files.
- Optionally a path to a directory of raw data (will be copied into the package).
- Optionally a list of additional code files that may be dependencies of your R scripts. 



```r
library(DataPackageR)

# Let's reproducibly package up
# the cars in the mtcars dataset
# with speed > 20.
# Our dataset will be called cars_over_20.

# Get the code file that turns the raw data
# to our packaged and processed analysis-ready dataset.
processing_code <-
  system.file("extdata", 
              "tests",
              "subsetCars.Rmd",
              package = "DataPackageR")

# Create the package framework.
DataPackageR::datapackage_skeleton(name = "mtcars20",
  force = TRUE,
  code_files = processing_code,
  r_object_names = "cars_over_20",
  path = tempdir() 
  #dependencies argument is empty
  #raw_data_dir argument is empty.
  ) 
```

### What's in the package skeleton structure?

This has created a datapackage source tree named "mtcars20" (in a temporary directory). 
For a real use case you would pick a `path` on your filesystem where you could then initialize a new github repository for the package.

The contents of `mtcars20` are:


```
Registered S3 methods overwritten by 'ggplot2':
  method         from 
  [.quosures     rlang
  c.quosures     rlang
  print.quosures rlang
                levelName
1  mtcars20              
2   ¦--DESCRIPTION       
3   ¦--R                 
4   ¦--Read-and-delete-me
5   ¦--data              
6   ¦--data-raw          
7   ¦   °--subsetCars.Rmd
8   ¦--datapackager.yml  
9   ¦--inst              
10  ¦   °--extdata       
11  °--man               
```

You should fill out the `DESCRIPTION` file to describe your data package. 
It contains a new `DataVersion` string that will be automatically incremented when the data package is built *if the packaged data has changed*. 

The user-provided code files reside in `data-raw`. They are executed during the data package build process.

### A few words about the YAML config file

A `datapackager.yml` file is used to configure and control the build process.

The contents are:


```
configuration:
  files:
    subsetCars.Rmd:
      enabled: yes
  objects: cars_over_20
  render_root:
    tmp: '69649'
```

The two main pieces of information in the configuration are a list of the files to be processed and the data sets the package will store.

This example packages an R data set named `cars_over_20` (the name was passed in to `datapackage_skeleton()`).
It is created by the `subsetCars.Rmd` file. 


The objects must be listed in the yaml configuration file. `datapackage_skeleton()`  ensures this is done for you automatically. 

DataPackageR provides an API for modifying this file, so it does not need to be done by hand. 

Further information on the contents of the YAML configuration file, and the API are in the [YAML Configuration Details](YAML_CONFIG.md)

### Where do I put my raw datasets?

Raw data (provided the size is not prohibitive) can be placed in `inst/extdata`.

The `datapackage_skeleton()` API has the `raw_data_dir` argument, which will copy the contents of `raw_data_dir`  (and its subdirectories) into `inst/extdata` automatically. 

In this example we are reading the `mtcars` data set that is already in memory, rather than from the file system.

### An API to read raw data sets from within an R or Rmd procesing script.

As stated in the README, in order for your processing scripts to be portable, you should not use absolute paths to files. 
DataPackageR provides an API to point to the data package root directory and the `inst/extdata` and `data` subdirectories.
These are useful for constructing portable paths in your code to read files from these locations.

For example: to construct a path to a file named "mydata.csv" located in `inst/extdata` in your data package source tree:

- use `DataPackageR::project_extdata_path("mydata.csv")` in your `R` or `Rmd` file. This would return: e.g., /var/folders/jh/x0h3v3pd4dd497g3gtzsm8500000gn/T//Rtmp765KFQ/mtcars20/inst/extdata/mydata.csv

Similarly: 

- `DataPackageR::project_path()`  constructs a path to the data package root directory. (e.g., /var/folders/jh/x0h3v3pd4dd497g3gtzsm8500000gn/T//Rtmp765KFQ/mtcars20)
- `DataPackageR::project_data_path()` constructs a path to the data package `data` subdirectory. (e.g., /var/folders/jh/x0h3v3pd4dd497g3gtzsm8500000gn/T//Rtmp765KFQ/mtcars20/data)

Raw data sets that are stored externally (outside the data package source tree) can be constructed relative to the `project_path()`.

### YAML header metadata for R files and Rmd files.

If your processing scripts are Rmd files, the usual yaml header for rmarkdown documents should be present.

If you have Rmd files, you can still include a yaml header, but it should be commented with `#'` and it should be at the top of your R file. For example, a test R file in the DataPackageR package looks as follows:

```
#'---
#'title: Sample report  from R script
#'author: Greg Finak
#'date: August 1, 2018
#'---
data <- runif(100)
```

This will be converted to an Rmd file with a proper yaml header, which will then be turned into a vignette and indexed in the built package.


## Build the data package.

Once the skeleton framework is set up, 


```r
# Run the preprocessing code to build cars_over_20
# and reproducibly enclose it in a package.
DataPackageR:::package_build(file.path(tempdir(),"mtcars20"))

✔ Setting active project to '/private/var/folders/jh/x0h3v3pd4dd497g3gtzsm8500000gn/T/Rtmp765KFQ/mtcars20'
✔ 1 data set(s) created by subsetCars.Rmd
• cars_over_20
☘ Built  all datasets!
Non-interactive NEWS.md file update.

✔ Creating 'vignettes/'
✔ Creating 'inst/doc/'
First time using roxygen2. Upgrading automatically...
Updating roxygen version in /private/var/folders/jh/x0h3v3pd4dd497g3gtzsm8500000gn/T/Rtmp765KFQ/mtcars20/DESCRIPTION
Writing NAMESPACE
Loading mtcars20
Writing mtcars20.Rd
Writing cars_over_20.Rd
  
   checking for file ‘/private/var/folders/jh/x0h3v3pd4dd497g3gtzsm8500000gn/T/Rtmp765KFQ/mtcars20/DESCRIPTION’ ...
  
✔  checking for file ‘/private/var/folders/jh/x0h3v3pd4dd497g3gtzsm8500000gn/T/Rtmp765KFQ/mtcars20/DESCRIPTION’

  
─  preparing ‘mtcars20’:

  
   checking DESCRIPTION meta-information ...
  
✔  checking DESCRIPTION meta-information

  
─  checking for LF line-endings in source and make files and shell scripts

  
─  checking for empty or unneeded directories

  
─  looking to see if a ‘data/datalist’ file should be added

  
     NB: this package now depends on R (>= 3.5.0)

  
     WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects: 'mtcars20/data/cars_over_20.rda'

  
─  building 'mtcars20_1.0.tar.gz'

  
   

Next Steps 
1. Update your package documentation. 
   - Edit the documentation.R file in the package source data-raw subdirectory and update the roxygen markup. 
   - Rebuild the package documentation with  document() . 
2. Add your package to source control. 
   - Call  git init .  in the package source root directory. 
   -  git add  the package files. 
   -  git commit  your new package. 
   - Set up a github repository for your pacakge. 
   - Add the github repository as a remote of your local package repository. 
   -  git push  your local repository to gitub. 
[1] "/private/var/folders/jh/x0h3v3pd4dd497g3gtzsm8500000gn/T/Rtmp765KFQ/mtcars20_1.0.tar.gz"
```

### Documenting your data set changes in NEWS.md

When you build a package in interactive mode, you will be
prompted to input text describing the changes to your data package (one line). 

These will appear in the NEWS.md file in the following format:

```
DataVersion: xx.yy.zz
========
A description of your changes to the package

[The rest of the file]
```


### Why not just use R CMD build?

If the processing script is time consuming or the data set is particularly large, then `R CMD build` would run the code each time the package is installed. In such cases, raw data may not be available, or the environment to do the data processing may not be set up for each user of the data. DataPackageR decouples data processing from package building/installation for data consumers.

### A log of the build process

DataPackageR uses the `futile.logger` package to log progress. 

If there are errors in the processing, the script will notify you via logging to console and to  `/private/tmp/Test/inst/extdata/Logfiles/processing.log`. Errors should be corrected and the build repeated.

If everything goes smoothly, you will have a new package built in the parent directory. 

In this case we have a new package 
`mtcars20_1.0.tar.gz`. 


### A note about the package source directory after building.

The pacakge source directory changes after the first build.


```
                         levelName
1  mtcars20                       
2   ¦--DATADIGEST                 
3   ¦--DESCRIPTION                
4   ¦--NAMESPACE                  
5   ¦--NEWS.md                    
6   ¦--R                          
7   ¦   °--mtcars20.R             
8   ¦--Read-and-delete-me         
9   ¦--data                       
10  ¦   °--cars_over_20.rda       
11  ¦--data-raw                   
12  ¦   ¦--documentation.R        
13  ¦   ¦--subsetCars.R           
14  ¦   °--subsetCars.Rmd         
15  ¦--datapackager.yml           
16  ¦--inst                       
17  ¦   ¦--doc                    
18  ¦   ¦   ¦--subsetCars.Rmd     
19  ¦   ¦   °--subsetCars.html    
20  ¦   °--extdata                
21  ¦       °--Logfiles           
22  ¦           ¦--processing.log 
23  ¦           °--subsetCars.html
24  ¦--man                        
25  ¦   ¦--cars_over_20.Rd        
26  ¦   °--mtcars20.Rd            
27  °--vignettes                  
28      °--subsetCars.Rmd         
```

### Update the autogenerated documentation. 

After the first build, the `R` directory contains `mtcars.R` that has autogenerated `roxygen2` markup documentation for the data package and for the packaged data `cars_over20`. 

The processed `Rd` files can be found in `man`. 

The autogenerated documentation source is in the `documentation.R` file in `data-raw`. 

You should update this file to properly document your objects. Then rebuild the documentation:


```r
document(file.path(tempdir(),"mtcars20"))

✔ Setting active project to '/private/var/folders/jh/x0h3v3pd4dd497g3gtzsm8500000gn/T/Rtmp765KFQ/mtcars20'
Updating mtcars20 documentation
Loading mtcars20
[1] TRUE
```

This is done without reprocessing the data.

#### Dont' forget to rebuild the package.

You should update the documentation in `R/mtcars.R`, then call `package_build()` again.


## Installing and using the new data package


### Accessing vignettes, data sets, and data set documentation. 

The package source also contains files in the `vignettes` and `inst/doc` directories that provide a log of the data processing. 

When the package is installed, these will be accessible via the `vignette()` API. 

The vignette will detail the processing performed by the `subsetCars.Rmd` processing script. 

The data set documentation will be accessible via `?cars_over_20`, and the data sets via `data()`. 


```r
# create a temporary library to install into.
dir.create(file.path(tempdir(),"lib"))
# Let's use the package we just created.
install.packages(file.path(tempdir(),"mtcars20_1.0.tar.gz"), type = "source", repos = NULL, lib = file.path(tempdir(),"lib"))
if(!"package:mtcars20"%in%search())
  attachNamespace('mtcars20') #use library() in your code
data("cars_over_20") # load the data

cars_over_20 # now we can use it.
   speed dist
44    22   66
45    23   54
46    24   70
47    24   92
48    24   93
49    24  120
50    25   85
?cars_over_20 # See the documentation you wrote in data-raw/documentation.R.
  
vignettes <- vignette(package = "mtcars20")
vignettes$results
      Package   
Topic "mtcars20"
      LibPath                                                         
Topic "/Library/Frameworks/R.framework/Versions/3.6/Resources/library"
      Item         Title                                            
Topic "subsetCars" "A Test Document for DataPackageR (source, html)"
```


### Using the DataVersion

Your downstream data analysis can depend on a specific version of the data in your data package by testing the DataVersion string in the DESCRIPTION file. 

We provide an API for this:


```r
# We can easily check the version of the data
DataPackageR::data_version("mtcars20")
[1] '0.1.0'

# You can use an assert to check the data version in  reports and
# analyses that use the packaged data.
assert_data_version(data_package_name = "mtcars20",
                    version_string = "0.1.0",
                    acceptable = "equal")  #If this fails, execution stops
                                           #and provides an informative error.
```


# Migrating old data packages.

Version 1.12.0 has moved away from controlling the build process using `datasets.R` and an additional `masterfile` argument. 

The build process is now controlled via a `datapackager.yml` configuration file located in the package root directory.  (see [YAML Configuration Details](YAML_CONFIG.md))

### Create a datapackager.yml file

You can migrate an old package by constructing such a config file using the `construct_yml_config()` API.


```r
# assume I have file1.Rmd and file2.R located in /data-raw, 
# and these create 'object1' and 'object2' respectively.

config <- construct_yml_config(code = c("file1.Rmd", "file2.R"),
                              data = c("object1", "object2"))
cat(yaml::as.yaml(config))
configuration:
  files:
    file1.Rmd:
      enabled: yes
    file2.R:
      enabled: yes
  objects:
  - object1
  - object2
  render_root:
    tmp: '251712'
```

`config` is a newly constructed yaml configuration object. It can be written to the package directory:


```r
path_to_package <- tempdir() #e.g., if tempdir() was the root of our package.
yml_write(config, path = path_to_package)
```

Now the package at `path_to_package` will build with version 1.12.0 or greater.

### Reading data sets from Rmd files

In versions prior to 1.12.1 we would read data sets from `inst/extdata` in an `Rmd` script using paths relative to
`data-raw` in the data package source tree. 

For example:

#### The old way

```r
# read 'myfile.csv' from inst/extdata relative to data-raw where the Rmd is rendered.
read.csv(file.path("../inst/extdata","myfile.csv"))
```

Now `Rmd` and `R` scripts are processed in `render_root` defined in the yaml config.

To read a raw data set we can get the path to the package source directory using an API call:


#### The new way

```r
# DataPackageR::project_extdata_path() returns the path to the data package inst/extdata subdirectory directory.
# DataPackageR::project_path() returns the path to the data package root directory.
# DataPackageR::project_data_path() returns the path to the data package data subdirectory directory.
read.csv(
    DataPackageR::project_extdata_path("myfile.csv")
    )
```

# Partial builds

We can also perform partial builds of a subset of files in a package by toggling the `enabled` key in the config file.

This can be done with the following API:


```r
config <- yml_disable_compile(config,filenames = "file2.R")
yml_write(config, path = path_to_package) # write modified yml to the package.
configuration:
  files:
    file1.Rmd:
      enabled: yes
    file2.R:
      enabled: no
  objects:
  - object1
  - object2
  render_root:
    tmp: '251712'
```

Note that the modified configuration needs to be written back to the package source directory in order for the 
changes to take effect. 

The consequence of toggling a file to `enable: no` is that it will be skipped when the package is rebuilt, 
but the data will still be retained in the package, and the documentation will not be altered. 

This is useful in situations where we have multiple data sets, and want to re-run one script to update a specific data set, but
not the other scripts because they may be too time consuming, for example.

# Multi-script pipelines.

We may have situations where we have mutli-script pipelines. There are two ways to share data among scripts. 

1. filesystem artifacts
2. data objects passed to subsequent scripts.

### File system artifacts

The yaml configuration property `render_root` specifies the working directory where scripts will be rendered.

If a script writes files to the working directory, that is where files will appear. These can be read by subsequent scripts.

### Passing data objects to subsequent scripts.

A script (e.g., `script2.Rmd`) running after `script1.Rmd` can access a stored data object named `script1_dataset` created by `script1.Rmd` by calling

`script1_dataset <- DataPackageR::datapackager_object_read("script1_dataset")`. 

Passing of data objects amongst scripts can be turned off via:

`package_build(deps = FALSE)`

# Next steps 

We recommend the following once your package is created.

## Place your package under source control

You now have a data package source tree. 

- **Place your package under version control**
    1. Call `git init` in the package source root to initialize a new git repository.
    2. [Create a new repository for your data package on github](https://help.github.com/articles/create-a-repo/).
    3. Push your local package repository to `github`. [see step 7](https://help.github.com/articles/adding-an-existing-project-to-github-using-the-command-line/)


This will let you version control your data processing code, and provide a mechanism for sharing your package with others.


For more details on using git and github with R, there is an excellent guide provided by Jenny Bryan: [Happy Git and GitHub for the useR](http://happygitwithr.com/) and Hadley Wickham's [book on R packages](http://r-pkgs.had.co.nz/).

# Additional Details

We provide some additional details for the interested.

### Fingerprints of stored data objects

DataPackageR calculates an md5 checksum of each data object it stores, and keeps track of them in a file
called `DATADIGEST`.

- Each time the package is rebuilt, the md5 sums of the new data objects are compared against the DATADIGEST.
- If they don't match, the build process checks that the `DataVersion` string has been incremented in the `DESCRIPTION` file.
- If it has not the build process will exit and produce an error message.

#### DATADIGEST


The `DATADIGEST` file contains the following:


```
DataVersion: 0.1.0
cars_over_20: 3ccb5b0aaa74fe7cfc0d3ca6ab0b5cf3
```


#### DESCRIPTION

The description file has the new `DataVersion` string.


```
Package: mtcars20
Title: What the Package Does (One Line, Title Case)
Version: 1.0
Authors@R: 
    person(given = "First",
           family = "Last",
           role = c("aut", "cre"),
           email = "first.last@example.com")
Description: What the package does (one paragraph).
License: What license it uses
Encoding: UTF-8
LazyData: true
DataVersion: 0.1.0
Roxygen: list(markdown = TRUE)
Date: 2019-03-11
Suggests: 
    knitr,
    rmarkdown
VignetteBuilder: knitr
RoxygenNote: 6.1.1
```




---
title: README
output: github_document
bibliography: bibliography.bib
editor_options: 
  chunk_output_type: inline
---


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# DataPackageR

DataPackageR is used to reproducibly process raw data into packaged, analysis-ready data sets.

<!-- badges: start -->
[![CRAN](https://www.r-pkg.org/badges/version/DataPackageR)]( https://CRAN.R-project.org/package=DataPackageR)
[![R-CMD-check](https://github.com/ropensci/DataPackageR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/DataPackageR/actions)
[![Coverage status](https://codecov.io/gh/ropensci/DataPackageR/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/DataPackageR?branch=master)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://badges.ropensci.org/230_status.svg)](https://github.com/ropensci/software-review/issues/230)
[![DOI](https://zenodo.org/badge/29267435.svg)](https://doi.org/10.5281/zenodo.1292095)
<!-- badges: end -->

- [yaml configuration guide](https://github.com/ropensci/DataPackageR/blob/main/vignettes/YAML_CONFIG.md)
- [a more detailed technical vignette](https://github.com/ropensci/DataPackageR/blob/main/vignettes/usingDataPackageR.md)

> **Important Note**: [datapack](https://github.com/ropensci/datapack) is a *different package* that is used to "create, send and load data from common repositories such as DataONE into the R environment".

> **This package** is for processing raw data into tidy data sets and bundling them into R packages.

## What problems does DataPackageR tackle?

You have diverse raw data sets that you need to preprocess and tidy in order to:

- Perform data analysis
- Write a report
- Publish a paper 
- Share data with colleagues and collaborators
- Save time in the future when you return to this project but have forgotten all about what you did.

### Why package data sets?

**Definition:** A *data package* is a formal R package whose sole purpose is to contain, access, and / or document data sets.

- **Reproducibility.**

  As described [elsewhere](https://github.com/ropensci/rrrpkg), packaging your data promotes reproducibility.
  R's packaging infrastructure promotes unit testing, documentation, a reproducible build system, and has many other benefits.
  Coopting it for packaging data sets is a natural fit.
  
- **Collaboration.**
  
  A data set packaged in R is easy to distribute and share amongst collaborators, and is easy to install and use.
  All the hard work you've put into documenting and standardizing the tidy data set comes right along with the data package.

- **Documentation.**

  R's package system allows us to document data objects. What's more, the `roxygen2` package makes this very easy to do with [markup tags](https://r-pkgs.org/data.html). 
  That documentation is the equivalent of a data dictionary and can be extremely valuable when returning to a project after a period of time.
  
- **Convenience.**

  Data pre-processing can be time consuming, depending on the data type and raw data sets may be too large to share conveniently in a packaged format. 
  Packaging and sharing the small, tidied data saves the users computing time and time spent waiting for downloads.
     
## Challenges.

- **Package size limits.**

  R packages have a 5MB size limit, at least on CRAN. BioConductor has explicit [data package](https://www.bioconductor.org/developers/package-guidelines/#package-types) types that can be larger and use git LFS for very large files. 
  
  Sharing large volumes of raw data in an R package format is still not ideal, and there are public biological data repositories better suited for raw data: e.g.,  [GEO](https://www.ncbi.nlm.nih.gov/geo/), [SRA](https://www.ncbi.nlm.nih.gov/sra), [ImmPort](https://www.immport.org:443/shared/immport-open/public/home/home), [ImmuneSpace](https://immunespace.org/), [FlowRepository](https://flowrepository.org/).
  
  Tools like [datastorr](https://github.com/ropenscilabs/datastorr) can help with this and we hope to integrate the into DataPackageR in the future.

- **Manual effort**

  There is still a substantial manual effort to set up the correct directory structures for an R data package. This can dissuade many individuals, particularly new users who have never built an R package, from going this route.
  
- **Scale** 

  Setting up and building R data packages by hand is a workable solution for a small project or a small number of projects, but when dealing with many projects each involving many data sets, tools are needed to help automate the process.
  
## DataPackageR

DataPackageR provides a number of benefits when packaging your data.

- It aims to automate away much of the tedium of packaging data sets without getting too much in the way, and keeps your processing workflow reproducible.

- It sets up the necessary package structure and files for a data package. 

- It allows you to keep the large, raw data and only ship the packaged tidy data, saving space and time consumers of your data set need to spend downloading and re-processing it.

- It maintains a reproducible record (vignettes) of the data processing along with the package. Consumers of the data package can verify how the processing was done, increasing confidence in your data. 

- It automates construction of the documentation and maintains a data set version and an md5 fingerprint of each data object in the package. If the data changes and the package is rebuilt, the data version is automatically updated.

## Similar work

There are a number of tools out there that address similar and complementary problems: 

- **datastorr**
  [github repo](https://github.com/ropenscilabs/datastorr)
  
  Simple data retrieval and versioning using GitHub to store data.
  
    - Caches downloads and uses github releases to version data.
    - Deal consistently with translating the file stored online into a loaded data object
    - Access multiple versions of the data at once 
    
  `datastorrr` could be used with DataPackageR to store / access remote raw data sets, remotely store / access tidied data that are too large to fit in the package itself.

- **fst**
  [github repo](https://github.com/fstpackage/fst)

  `fst` provides lightning fast serialization of data frames.

- **The modern data package**
  [pdf](https://github.com/noamross/2018-04-18-rstats-nyc/blob/master/Noam_Ross_ModernDataPkg_rstatsnyc_2018-04-20.pdf)

  A presentation from \@noamross touching on modern tools for open science and reproducibility. Discusses `datastorr` and `fst` as well as standardized metadata and documentation.
 
- **rrrpkg**
  [github repo](https://github.com/ropensci/rrrpkg)

  A document from ropensci describing using an R package as a research compendium. Based on ideas originally introduced by Robert Gentleman and Duncan Temple Lang (Gentleman and Lang (2004)<!--@Gentleman2004-oj-->)

- **template**
  [github repo](https://github.com/ropensci/rrrpkg)
  
  An R package template for data packages.
 
 See the [publication](#publication) for further discussion.
 
## Installation

You can install the latest version of DataPackageR from [github](https://github.com/ropensci/DataPackageR) with:

```{r, eval=FALSE}
library(devtools)
devtools::install_github("ropensci/DataPackageR")
```

## Blog Post - building packages interactively.

See this [rOpenSci blog post](https://ropensci.org/blog/2018/09/18/datapackager/) on how to build data packages interactively using DataPackageR. 
This uses several new interfaces: `use_data_object()`, `use_processing_script()` and `use_raw_dataset()` to build up a data package, rather than assuming
the user has all the code and data ready to go for `datapackage_skeleton()`. 

## Example (assuming all code and data are available)

```{r minimal_example, results='hide', message=FALSE}
library(DataPackageR)

# Let's reproducibly package up
# the cars in the mtcars dataset
# with speed > 20.
# Our dataset will be called cars_over_20.
# There are three steps:

# 1. Get the code file that turns the raw data
# into our packaged and processed analysis-ready dataset.
# This is in a file called subsetCars.Rmd located in exdata/tests of the DataPackageR package.
# For your own projects you would write your own Rmd processing file.
processing_code <- system.file(
  "extdata", "tests", "subsetCars.Rmd", package = "DataPackageR"
)

# 2. Create the package framework.
# We pass in the Rmd file in the `processing_code` variable and the names of the data objects it creates (called "cars_over_20")
# The new package is called "mtcars20"
datapackage_skeleton(
  "mtcars20", force = TRUE, 
  code_files = processing_code, 
  r_object_names = "cars_over_20", 
  path = tempdir()) 

# 3. Run the preprocessing code to build the cars_over_20 data set 
# and reproducibly enclose it in the mtcars20 package.
# packageName is the full path to the package source directory created at step 2.
# You'll be prompted for a text description (one line) of the changes you're making.
# These will be added to the NEWS.md file along with the DataVersion in the package source directory.
# If the build is run in non-interactive mode, the description will read
# "Package built in non-interactive mode". You may update it later.
dir.create(file.path(tempdir(),"lib"))
package_build(packageName = file.path(tempdir(),"mtcars20"), install = TRUE, lib = file.path(tempdir(),"lib"))

# Update the autogenerated roxygen documentation in data-raw/documentation.R. 
# edit(file.path(tempdir(),"mtcars20","R","mtcars20.R"))

# 4. Rebuild the documentation.
document(file.path(tempdir(),"mtcars20"), install = TRUE, lib = file.path(tempdir(),"lib"))

# Let's use the package we just created.
install.packages(file.path(tempdir(),"mtcars20_1.0.tar.gz"), type = "source", repos = NULL)
library(mtcars20)
data("cars_over_20") # load the data
cars_over_20  # Now we can use it.
?cars_over_20 # See the documentation you wrote in data-raw/documentation.R.

# We have our dataset!
# Since we preprocessed it,
# it is clean and under the 5 MB limit for data in packages.
cars_over_20

# We can easily check the version of the data
data_version("mtcars20")

# You can use an assert to check the data version in  reports and
# analyses that use the packaged data.
assert_data_version(data_package_name = "mtcars20",
                    version_string = "0.1.0",
                    acceptable = "equal")
```

### Reading external data from within R / Rmd processing scripts.

When creating a data package, your processing scripts will need to read your raw data sets in order to process them.
These data sets can be stored in `inst/extdata` of the data package source tree, or elsewhere outside the package source tree.
In order to have portable and reproducible code, you should not use absolute paths to the raw data.
Instead, `DataPackageR` provides several APIs to access the data package project root directory, the `inst/extdata` subdirectory, and the `data` subdirectory.

```{r, eval = FALSE}
# This returns the datapackage source 
# root directory. 
# In an R or Rmd processing script this can be used to build a path to a directory that is exteral to the package, for 
# example if we are dealing with very large data sets where data cannot be packaged.
DataPackageR::project_path()

# This returns the   
# inst/extdata directory. 
# Raw data sets that are included in the package should be placed there.
# They can be read from that location, which is returned by: 
DataPackageR::project_extdata_path()

# This returns the path to the datapackage  
# data directory. This can be used to access 
# stored data objects already created and saved in `data` from 
# other processing scripts.
DataPackageR::project_data_path()
```


## Preprint and publication. <a id = "publication"></a>

The publication describing the package, (Finak *et* *al.,* 2018)<!--@10.12688/gatesopenres.12832.2-->, is now available at   [Gates Open Research](https://gatesopenresearch.org/articles/2-31/v1) .


The preprint is on [biorxiv](https://doi.org/10.1101/342907).


## Code of conduct 

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/DataPackageR/blob/main/CODE_OF_CONDUCT.md).
  By participating in this project you agree to abide by its terms.

### References

1. Gentleman, Robert, and Duncan Temple Lang. 2004. “Statistical Analyses and Reproducible Research.” Bioconductor Project Working Papers, Bioconductor project working papers,. bepress.

2. Finak G, Mayer B, Fulp W et al. DataPackageR: Reproducible data preprocessing, standardization and sharing using R/Bioconductor for collaborative data analysis [version 1; referees: 1 approved with reservations]. Gates Open Res 2018, 2:31
(doi: 10.12688/gatesopenres.12832.1)



  [![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)


---
title: "A Test Document for DataPackageR"
author: "Greg Finak"
output: html_document
---

```{r include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

This is a simple Rmd file that demonstrates how DataPackageR processes Rmarkdown files and creates data sets
that are then stored in an R data package. 

In the `config.yml` for this example, this file is listed first, and therefore processed first.

This particular document simply subsets the `cars` data set:

```{r cars}
summary(cars)
dim(cars)
```

`cars` consists of a data frame of 50 rows and two columns. The `?cars` documentation specifies that it consists of speed and stopping distances of cars.

Let's say, for some reason, we are only interested in the stopping distances of cars traveling greater than 20 miles per hour.

```{r}
cars_over_20 = subset(cars, speed > 20)
```

The data frame `cars_over_20` now holds this information. 

# Storing data set objects and making making accessible to other processing scripts. 

When DataPackageR processes this file, it creates this `cars_over_20` object. After processing the file it does several things:

1. It compares the objects in the rmarkdown render environment of `subsetCars.Rmd` against the objects listed in the `config.yml` file `objects` property. 
2. It finds `cars_over_20` is listed there, so it stores it in a new environment.
3. That environment is passed to subsequent R and Rmd files. Specifically when the `extra.rmd` file is processed, it has access to an environment object that holds all the `objects` (defined in the yaml config) that have already been created and processed. This environment is passed into subsequent scripts at the `render()` call. 

All of the above is done automatically. The user only needs to list the objects to be stored and passed to other scripts in the `config.yml` file. 

The `datapackager_object_read()` API can be used to retrieve these objects from the environment. 

### Storing objects in the data package

In addition to passing around an environment to subsequent scripts, the `cars_over_20` object is stored in the data package `/data` directory as an `rda` file. 

Note that this is all done automatically. The user does not need to explicitly save anything, they only need to list the objects to be store in the `config.yml`. 

This object is then accessible in the resulting package via the `data()` API, and its documentation is accessible via `?cars_over_20`. 

### Data object documentation

The documentation for the `cars_over_20` object is created in a `subsetCars.R` file in the `/R` directory of the data package. 

While the data object document stub is created automatically, it must be edited by the user to provide additional details about the data object.


---
title: "The DataPackageR YAML configuration file."
author: "Greg Finak <gfinak@fredhutch.org>"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    keep_md: TRUE
    toc: yes
  bibliography: bibliography.bib
vignette: >
  %\VignetteIndexEntry{The DataPackageR YAML configuration file.}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
  %\usepackage{graphicx}
editor_options: 
  chunk_output_type: inline
---

# Configuring and controlling DataPackageR builds.

Data package builds are controlled using the `datapackager.yml` file. 

This file is created in the package source tree when the user creates a package using `datapackage_skeleton()`. 

It is automatically populated with the names of the `code_files` and `data_objects` the passed in to datapackage_skeleton.

## The `datapackager.yml` file.

The structure of a correctly formatted `datapackager.yml` file is shown below:

```{r, echo = FALSE, results = 'hide', eval = rmarkdown::pandoc_available()}
library(DataPackageR)
library(yaml)
yml <- DataPackageR::construct_yml_config(code = "subsetCars.Rmd", data = "cars_over_20")
```

```{r, echo = FALSE, comment="", eval = rmarkdown::pandoc_available()}
cat(yaml::as.yaml(yml))
```

## YAML config file properties.

The main section of the file is the `configuration:` section.

It has three properties:

- `files:` 
  
  The files (`R` or `Rmd`) to be processed by DataPackageR. They are processed in the order shown. Users running multi-script workflows with dependencies between the scripts need to ensure the files are processed in the correct order. 
  
  Here `subsetCars.Rmd` is the only file to process. The name is transformed to an absolute path within the package.
  
  Each file itself has just one property:
  
  <!-- - `name:`  -->
     <!-- The name of the file. This is transformed to an absolute path within the package. -->
    
  - `enabled:` 
     A logical `yes`, `no` flag indicating whether the file should be rendered during the build, or whether it should be skipped.
    This is useful for 'turning off' long running processing tasks if they have not changed. Disabling processing of a file will not overwrite existing documentation or data objecs created during previous builds.
  
- `objects:` 
  
  The names of the data objects created by the processing files, to be stored in the package. These names are compared against the objects created in the render environment by each file. They names must match.
- `render_root:` 
  
  The directory where the `Rmd` or `R` files will be rendered. Defaults to a randomly named subdirectory of `tempdir()`. Allows workflows that use multiple scripts and create file system artifacts to function correctly by simply writing to and reading from the working directory.

## Editing the YAML config file.

The structure of the YAML is simple enough to understand but complex enough that it can be a pain to edit by hand.

DataPackageR provides a number of API calls to construct, read, modify, and write the yaml config file.

### API calls

#### `construct_yml_config` 

  Make an r object representing a YAML config file.
  
##### Example
  The YAML config shown above was created by: 
```{r, eval = rmarkdown::pandoc_available()}
# Note this is done by the datapackage_skeleton. 
# The user doesn't usually need to call 
# construct_yml_config()
yml <- DataPackageR::construct_yml_config(
  code = "subsetCars.Rmd",
  data = "cars_over_20"
  )
```
  
  
#### `yml_find` 

  Read a yaml config file from a package path into an r object.
  
##### Example
  Read the YAML config file from the `mtcars20` example.
  
```{r eval=FALSE}
# returns an r object representation of
# the config file.
mtcars20_config <- yml_find(
  file.path(tempdir(),"mtcars20")
  )
```

#### `yml_list_objects` 

  List the `objects` in a config read by `yml_find`.
  

##### Example

```{r, comment="", eval = rmarkdown::pandoc_available()}
  yml_list_objects(yml)
```
  
#### `yml_list_files` 

  List the `files` in a config read by `yml_find`.

##### Example

```{r, comment="", eval = rmarkdown::pandoc_available()}
  yml_list_files(yml)
```
  
#### `yml_disable_compile` 

  Disable compilation of named files in a config read by `yml_find`.

##### Example

```{r, comment="", echo = 1, eval = rmarkdown::pandoc_available()}
yml_disabled <- yml_disable_compile(
    yml,
    filenames = "subsetCars.Rmd")
cat(as.yaml(yml_disabled))
```

#### `yml_enable_compile` 

  Enable compilation of named files in a config read by `yml_find`.

##### Example

```{r, comment="", echo = 1, eval = rmarkdown::pandoc_available()}
yml_enabled <- yml_enable_compile(
    yml,
    filenames = "subsetCars.Rmd")
cat(as.yaml(yml_enabled))
```

#### `yml_add_files` 
  
  Add named files to a config read by `yml_find`.
  
##### Example

```{r, comment="", echo = 1, eval = rmarkdown::pandoc_available()}
yml_twofiles <- yml_add_files(
    yml,
    filenames = "anotherFile.Rmd")
# cat(as.yaml(yml_twofiles))
```

#### `yml_add_objects` 

  Add named objects to a config read by `yml_find`.

##### Example

```{r, comment="", echo = 1, eval = rmarkdown::pandoc_available()}
yml_twoobj <- yml_add_objects(
    yml_twofiles,
    objects = "another_object")
# cat(as.yaml(yml_twoobj))
```

#### `yml_remove_files` 

  Remove named files from a config read by `yml_find`.

##### Example

```{r, comment="", echo = 1, eval = rmarkdown::pandoc_available()}
yml_twoobj <- yml_remove_files(
    yml_twoobj,
    filenames = "anotherFile.Rmd")
# cat(as.yaml(yml_twoobj))
```

#### `yml_remove_objects` 

  Remove named objects from a config read by `yml_find`.

##### Example

```{r, comment="", echo = 1, eval = rmarkdown::pandoc_available()}
yml_oneobj <- yml_remove_objects(
    yml_twoobj,
    objects = "another_object")
# cat(as.yaml(yml_oneobj))
```

#### `yml_write` 

  Write a modified config to its package path.

##### Example

```{r, eval = FALSE}
yml_write(yml_oneobj, path = "path_to_package")
```

The `yml_oneobj` read by `yml_find()` carries an attribute
that is the path to the package. The user doesn't need to pass a `path` to `yml_write` if the config has been read by `yml_find`.
---
title: "Using DataPackageR"
author: "Greg Finak <gfinak@fredhutch.org>"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    keep_md: TRUE
    toc: yes
vignette: >
  %\VignetteIndexEntry{A Guide to using DataPackageR}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
  %\usepackage{graphicx}
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "",
  eval = TRUE
)
```

## Purpose

This vignette demonstrates how to use DataPackageR to build a data package. 

DataPackageR aims to simplify data package construction.

It provides mechanisms for reproducibly preprocessing and tidying raw data into into documented, versioned, and packaged analysis-ready data sets. 

Long-running or computationally intensive data processing can be decoupled from the usual `R CMD build` process while maintinaing [data lineage](https://en.wikipedia.org/wiki/Data_lineage).

In this vignette we will subset and package the `mtcars` data set.

## Set up a new data package.

We'll set up a new data package based on `mtcars` example in the [README](https://github.com/ropensci/DataPackageR/blob/master/README.md).
The `datapackage_skeleton()` API is used to set up a new package. 
The user needs to provide:

- R or Rmd code files that do data processing.
- A list of R object names created by those code files.
- Optionally a path to a directory of raw data (will be copied into the package).
- Optionally a list of additional code files that may be dependencies of your R scripts. 


```{r minimal_example, results='hide', eval = rmarkdown::pandoc_available()}
library(DataPackageR)

# Let's reproducibly package up
# the cars in the mtcars dataset
# with speed > 20.
# Our dataset will be called cars_over_20.

# Get the code file that turns the raw data
# to our packaged and processed analysis-ready dataset.
processing_code <-
  system.file("extdata", 
              "tests",
              "subsetCars.Rmd",
              package = "DataPackageR")

# Create the package framework.
DataPackageR::datapackage_skeleton(name = "mtcars20",
  force = TRUE,
  code_files = processing_code,
  r_object_names = "cars_over_20",
  path = tempdir() 
  #dependencies argument is empty
  #raw_data_dir argument is empty.
  )
```

### What's in the package skeleton structure?

This has created a datapackage source tree named "mtcars20" (in a temporary directory). 
For a real use case you would pick a `path` on your filesystem where you could then initialize a new github repository for the package.

The contents of `mtcars20` are:

```{r dirstructure,echo=FALSE, eval = rmarkdown::pandoc_available()}
library(data.tree)
df <- data.frame(pathString = file.path(
  "mtcars20",
  list.files(
  file.path(tempdir(), "mtcars20"),
  include.dirs = TRUE,
  recursive = TRUE
  )
  ))
as.Node(df)
```

You should fill out the `DESCRIPTION` file to describe your data package. 
It contains a new `DataVersion` string that will be automatically incremented when the data package is built *if the packaged data has changed*. 

The user-provided code files reside in `data-raw`. They are executed during the data package build process.

### A few words about the YAML config file

A `datapackager.yml` file is used to configure and control the build process.

The contents are:

```{r, echo=FALSE, eval = rmarkdown::pandoc_available()}
cat(yaml::as.yaml(yaml::yaml.load_file(file.path(tempdir(),"mtcars20","datapackager.yml"))))
```

The two main pieces of information in the configuration are a list of the files to be processed and the data sets the package will store.

This example packages an R data set named `cars_over_20` (the name was passed in to `datapackage_skeleton()`).
It is created by the `subsetCars.Rmd` file. 


The objects must be listed in the yaml configuration file. `datapackage_skeleton()`  ensures this is done for you automatically. 

DataPackageR provides an API for modifying this file, so it does not need to be done by hand. 

Further information on the contents of the YAML configuration file, and the API are in the [YAML Configuration Details](YAML_CONFIG.html)

### Where do I put my raw datasets?

Raw data (provided the size is not prohibitive) can be placed in `inst/extdata`.

The `datapackage_skeleton()` API has the `raw_data_dir` argument, which will copy the contents of `raw_data_dir`  (and its subdirectories) into `inst/extdata` automatically. 

In this example we are reading the `mtcars` data set that is already in memory, rather than from the file system.

### An API to read raw data sets from within an R or Rmd procesing script.

As stated in the README, in order for your processing scripts to be portable, you should not use absolute paths to files. 
DataPackageR provides an API to point to the data package root directory and the `inst/extdata` and `data` subdirectories.
These are useful for constructing portable paths in your code to read files from these locations.

For example: to construct a path to a file named "mydata.csv" located in `inst/extdata` in your data package source tree:

- use `DataPackageR::project_extdata_path("mydata.csv")` in your `R` or `Rmd` file. This would return: e.g., `r file.path(tempdir(),"mtcars20","inst","extdata","mydata.csv")`

Similarly: 

- `DataPackageR::project_path()`  constructs a path to the data package root directory. (e.g., `r file.path(tempdir(),"mtcars20")`)
- `DataPackageR::project_data_path()` constructs a path to the data package `data` subdirectory. (e.g., `r file.path(tempdir(),"mtcars20","data")`)

Raw data sets that are stored externally (outside the data package source tree) can be constructed relative to the `project_path()`.

### YAML header metadata for R files and Rmd files.

If your processing scripts are Rmd files, the usual yaml header for rmarkdown documents should be present.

If you have Rmd files, you can still include a yaml header, but it should be commented with `#'` and it should be at the top of your R file. For example, a test R file in the DataPackageR package looks as follows:

```
#'---
#'title: Sample report  from R script
#'author: Greg Finak
#'date: August 1, 2018
#'---
data <- runif(100)
```

This will be converted to an Rmd file with a proper yaml header, which will then be turned into a vignette and indexed in the built package.


## Build the data package.

Once the skeleton framework is set up, 

```{r , eval = rmarkdown::pandoc_available()}
# Run the preprocessing code to build cars_over_20
# and reproducibly enclose it in a package.
dir.create(file.path(tempdir(),"lib"))
DataPackageR:::package_build(file.path(tempdir(),"mtcars20"), install = TRUE,  lib = file.path(tempdir(),"lib"))
```

### Documenting your data set changes in NEWS.md

When you build a package in interactive mode, you will be
prompted to input text describing the changes to your data package (one line). 

These will appear in the NEWS.md file in the following format:

```
DataVersion: xx.yy.zz
========
A description of your changes to the package

[The rest of the file]
```


### Why not just use R CMD build?

If the processing script is time consuming or the data set is particularly large, then `R CMD build` would run the code each time the package is installed. In such cases, raw data may not be available, or the environment to do the data processing may not be set up for each user of the data. DataPackageR decouples data processing from package building/installation for data consumers.

### A log of the build process

DataPackageR uses the `futile.logger` package to log progress. 

If there are errors in the processing, the script will notify you via logging to console and to  `/private/tmp/Test/inst/extdata/Logfiles/processing.log`. Errors should be corrected and the build repeated.

If everything goes smoothly, you will have a new package built in the parent directory. 

In this case we have a new package 
`mtcars20_1.0.tar.gz`. 


### A note about the package source directory after building.

The pacakge source directory changes after the first build.

```{r, echo=FALSE, eval = rmarkdown::pandoc_available()}
df <- data.frame(pathString = file.path(
  "mtcars20",
  list.files(
  file.path(tempdir(), "mtcars20"),
  include.dirs = TRUE,
  recursive = TRUE
  )
  ))
  as.Node(df)
```

### Update the autogenerated documentation. 

After the first build, the `R` directory contains `mtcars.R` that has autogenerated `roxygen2` markup documentation for the data package and for the packaged data `cars_over20`. 

The processed `Rd` files can be found in `man`. 

The autogenerated documentation source is in the `documentation.R` file in `data-raw`. 

You should update this file to properly document your objects. Then rebuild the documentation:

```{r rebuild_docs, eval = rmarkdown::pandoc_available()}
dir.create(file.path(tempdir(),"lib")) # a temporary library directory
document(file.path(tempdir(),"mtcars20"), lib = file.path(tempdir(),"lib"))
```

This is done without reprocessing the data.

#### Dont' forget to rebuild the package.

You should update the documentation in `R/mtcars.R`, then call `package_build()` again.


## Installing and using the new data package


### Accessing vignettes, data sets, and data set documentation. 

The package source also contains files in the `vignettes` and `inst/doc` directories that provide a log of the data processing. 

When the package is installed, these will be accessible via the `vignette()` API. 

The vignette will detail the processing performed by the `subsetCars.Rmd` processing script. 

The data set documentation will be accessible via `?cars_over_20`, and the data sets via `data()`. 

```{r, eval = rmarkdown::pandoc_available()}
# create a temporary library to install into.
dir.create(file.path(tempdir(),"lib"))
# Let's use the package we just created.
install.packages(file.path(tempdir(),"mtcars20_1.0.tar.gz"), type = "source", repos = NULL, lib = file.path(tempdir(),"lib"))
lns <- loadNamespace    
if (!"package:mtcars20"%in%search())
  attachNamespace(lns('mtcars20',lib.loc = file.path(tempdir(),"lib"))) #use library() in your code
data("cars_over_20") # load the data

cars_over_20 # now we can use it.
?cars_over_20 # See the documentation you wrote in data-raw/documentation.R.
  
vignettes <- vignette(package = "mtcars20", lib.loc = file.path(tempdir(),"lib"))
vignettes$results
```


### Using the DataVersion

Your downstream data analysis can depend on a specific version of the data in your data package by testing the DataVersion string in the DESCRIPTION file. 

We provide an API for this:

```{r, eval = rmarkdown::pandoc_available()}
# We can easily check the version of the data
DataPackageR::data_version("mtcars20", lib.loc = file.path(tempdir(),"lib"))

# You can use an assert to check the data version in  reports and
# analyses that use the packaged data.
assert_data_version(data_package_name = "mtcars20",
                    version_string = "0.1.0",
                    acceptable = "equal",
                    lib.loc = file.path(tempdir(),"lib"))  #If this fails, execution stops
                                           #and provides an informative error.
```


# Migrating old data packages.

Version 1.12.0 has moved away from controlling the build process using `datasets.R` and an additional `masterfile` argument. 

The build process is now controlled via a `datapackager.yml` configuration file located in the package root directory.  (see [YAML Configuration Details](YAML_CONFIG.html))

### Create a datapackager.yml file

You can migrate an old package by constructing such a config file using the `construct_yml_config()` API.

```{r  construct_config, eval = rmarkdown::pandoc_available()}
# assume I have file1.Rmd and file2.R located in /data-raw, 
# and these create 'object1' and 'object2' respectively.

config <- construct_yml_config(code = c("file1.Rmd", "file2.R"),
                              data = c("object1", "object2"))
cat(yaml::as.yaml(config))
```

`config` is a newly constructed yaml configuration object. It can be written to the package directory:

```{r, eval = rmarkdown::pandoc_available()}
path_to_package <- tempdir() #e.g., if tempdir() was the root of our package.
yml_write(config, path = path_to_package)
```

Now the package at `path_to_package` will build with version 1.12.0 or greater.

### Reading data sets from Rmd files

In versions prior to 1.12.1 we would read data sets from `inst/extdata` in an `Rmd` script using paths relative to
`data-raw` in the data package source tree. 

For example:

#### The old way

```r
# read 'myfile.csv' from inst/extdata relative to data-raw where the Rmd is rendered.
read.csv(file.path("../inst/extdata","myfile.csv"))
```

Now `Rmd` and `R` scripts are processed in `render_root` defined in the yaml config.

To read a raw data set we can get the path to the package source directory using an API call:


#### The new way

```r
# DataPackageR::project_extdata_path() returns the path to the data package inst/extdata subdirectory directory.
# DataPackageR::project_path() returns the path to the data package root directory.
# DataPackageR::project_data_path() returns the path to the data package data subdirectory directory.
read.csv(
    DataPackageR::project_extdata_path("myfile.csv")
    )
```

# Partial builds

We can also perform partial builds of a subset of files in a package by toggling the `enabled` key in the config file.

This can be done with the following API:

```{r echo=1:2, eval = rmarkdown::pandoc_available()}
config <- yml_disable_compile(config,filenames = "file2.R")
yml_write(config, path = path_to_package) # write modified yml to the package.
cat(yaml::as.yaml(config))
```

Note that the modified configuration needs to be written back to the package source directory in order for the 
changes to take effect. 

The consequence of toggling a file to `enable: no` is that it will be skipped when the package is rebuilt, 
but the data will still be retained in the package, and the documentation will not be altered. 

This is useful in situations where we have multiple data sets, and want to re-run one script to update a specific data set, but
not the other scripts because they may be too time consuming, for example.

# Multi-script pipelines.

We may have situations where we have mutli-script pipelines. There are two ways to share data among scripts. 

1. filesystem artifacts
2. data objects passed to subsequent scripts.

### File system artifacts

The yaml configuration property `render_root` specifies the working directory where scripts will be rendered.

If a script writes files to the working directory, that is where files will appear. These can be read by subsequent scripts.

### Passing data objects to subsequent scripts.

A script (e.g., `script2.Rmd`) running after `script1.Rmd` can access a stored data object named `script1_dataset` created by `script1.Rmd` by calling

`script1_dataset <- DataPackageR::datapackager_object_read("script1_dataset")`. 

Passing of data objects amongst scripts can be turned off via:

`package_build(deps = FALSE)`

# Next steps 

We recommend the following once your package is created.

## Place your package under source control

You now have a data package source tree. 

- **Place your package under version control**
    1. Call `git init` in the package source root to initialize a new git repository.
    2. [Create a new repository for your data package on github](https://help.github.com/articles/create-a-repo/).
    3. Push your local package repository to `github`. [see step 7](https://help.github.com/articles/adding-an-existing-project-to-github-using-the-command-line/)


This will let you version control your data processing code, and provide a mechanism for sharing your package with others.


For more details on using git and github with R, there is an excellent guide provided by Jenny Bryan: [Happy Git and GitHub for the useR](https://happygitwithr.com/) and Hadley Wickham's [book on R packages]( https://r-pkgs.org/).

# Additional Details

We provide some additional details for the interested.

### Fingerprints of stored data objects

DataPackageR calculates an md5 checksum of each data object it stores, and keeps track of them in a file
called `DATADIGEST`.

- Each time the package is rebuilt, the md5 sums of the new data objects are compared against the DATADIGEST.
- If they don't match, the build process checks that the `DataVersion` string has been incremented in the `DESCRIPTION` file.
- If it has not the build process will exit and produce an error message.

#### DATADIGEST


The `DATADIGEST` file contains the following:

```{r, echo=FALSE, eval = rmarkdown::pandoc_available()}
cat(readLines(file.path(tempdir(),"mtcars20","DATADIGEST")),sep="\n")
```


#### DESCRIPTION

The description file has the new `DataVersion` string.

```{r echo=FALSE, eval = rmarkdown::pandoc_available()}
cat(readLines(file.path(tempdir(),"mtcars20","DESCRIPTION")),sep="\n")
```




% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use.R
\name{use_processing_script}
\alias{use_processing_script}
\title{Add a processing script to a data package.}
\usage{
use_processing_script(
  file = NULL,
  title = NULL,
  author = NULL,
  overwrite = FALSE
)
}
\arguments{
\item{file}{\code{character} path to an existing file or name of a new R or Rmd file to create.}

\item{title}{\code{character} title of the processing script for the yaml header. Used only if file is being created.}

\item{author}{\code{character} author name for the yaml header. Used only if the file is being created.}

\item{overwrite}{\code{logical} default FALSE. Overwrite existing file of the same name.}
}
\value{
invisibly returns TRUE for success. Stops on failure.
}
\description{
The Rmd or R file or directory specified by \code{file} will be moved into
the data-raw directory. It will also be added to the yml configuration file.
Any existing file by that name will be overwritten when overwrite is set to TRUE
}
\examples{
if(rmarkdown::pandoc_available()){
myfile <- tempfile()
file <- system.file("extdata", "tests", "extra.rmd",
                     package = "DataPackageR")
datapackage_skeleton(
  name = "datatest",
  path = tempdir(),
  code_files = file,
  force = TRUE,
  r_object_names = "data")
use_processing_script(file = "newScript.Rmd",
    title = "Processing a new dataset",
    author = "Y.N. Here.")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataversion.r
\name{data_version}
\alias{data_version}
\alias{dataVersion}
\title{Get the DataVersion for a package}
\usage{
data_version(pkg, lib.loc = NULL)

dataVersion(pkg, lib.loc = NULL)
}
\arguments{
\item{pkg}{\code{character} the package name}

\item{lib.loc}{\code{character} path to library location.}
}
\description{
Retrieves the DataVersion of a package if available
}
\note{
\code{dataVersion()} has been renamed to \code{data_version()}
}
\examples{
if(rmarkdown::pandoc_available()){
f <- tempdir()
f <- file.path(f,"foo.Rmd")
con <- file(f)
writeLines("```{r}\n tbl = table(sample(1:10,1000,replace=TRUE)) \n```\n",con=con)
close(con)
pname <- basename(tempfile())
datapackage_skeleton(name = pname,
   path=tempdir(),
   force = TRUE,
   r_object_names = "tbl",
   code_files = f)

   package_build(file.path(tempdir(),pname), install = FALSE)

   devtools::load_all(file.path(tempdir(),pname))
   data_version(pname)
}
}
\seealso{
\code{\link[utils]{packageVersion}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use.R
\name{use_raw_dataset}
\alias{use_raw_dataset}
\title{Add a raw data set to inst/extdata}
\usage{
use_raw_dataset(path = NULL, ignore = FALSE)
}
\arguments{
\item{path}{\code{character} path to file or directory.}

\item{ignore}{\code{logical} whether to ignore the path or file in git and R build.}
}
\value{
invisibly returns TRUE for success. Stops on failure.
}
\description{
The file or directory specified by \code{path} will be moved into
the inst/extdata directory.
}
\examples{
if(rmarkdown::pandoc_available()){
myfile <- tempfile()
file <- system.file("extdata", "tests", "extra.rmd",
                     package = "DataPackageR")
raw_data <- system.file("extdata", "tests", "raw_data",
                        package = "DataPackageR")
datapackage_skeleton(
  name = "datatest",
  path = tempdir(),
  code_files = file,
  force = TRUE,
  r_object_names = "data")
use_raw_dataset(raw_data)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{keepDataObjects-defunct}
\alias{keepDataObjects-defunct}
\alias{keepDataObjects}
\title{These functions are no longer available.}
\usage{
keepDataObjects(...)
}
\arguments{
\item{...}{arguments}
}
\description{
These functions are no longer available.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processData.R
\docType{package}
\name{DataPackageR-package}
\alias{DataPackageR-package}
\title{DataPackageR}
\description{
A framework to automate the processing, tidying and packaging of raw data into analysis-ready
data sets as R packages.
}
\details{
DataPackageR will automate running of data processing code,
storing tidied data sets in an R package, producing
data documentation stubs, tracking data object finger prints (md5 hash)
and tracking and incrementing a "DataVersion" string
in the DESCRIPTION file of the package when raw data or data
objects change.
Code to perform the data processing is passed to DataPackageR by the user.
The user also specifies the names of the tidy data objects to be stored,
documented and tracked in the final package. Raw data should be read from
"inst/extdata" but large raw data files can be read from sources external
to the package source tree.

Configuration is controlled via the config.yml file created at the package root.
Its properties include a list of R and Rmd files that are to be rendered / sourced and
which read data and do the actual processing.
It also includes a list of r object names created by those files. These objects
are stored in the final package and accessible via the \code{data()} API.
The documentation for these objects is accessible via "?object-name", and md5
fingerprints of these objects are created and tracked.

The Rmd and R files used to process the objects are transformed into vignettes
accessible in the final package so that the processing is fully documented.

A DATADIGEST file in the package source keeps track of the data object fingerprints.
A DataVersion string is added to the package DESCRIPTION file and updated when these
objects are updated or changed on subsequent builds.

Once the package is built and installed, the data objects created in the package are accessible via
the \code{data()} API, and
Calling \code{datapackage_skeleton()} and passing in R / Rmd file names, and r object names
constructs a skeleton data package source tree and an associated \code{config.yml} file.

Calling \code{build_package()} sets the build process in motion.
}
\examples{
# A simple Rmd file that creates one data object
# named "tbl".
if(rmarkdown::pandoc_available()){
f <- tempdir()
f <- file.path(f,"foo.Rmd")
con <- file(f)
writeLines("```{r}\n tbl = table(sample(1:10,1000,replace=TRUE)) \n```\n",con=con)
close(con)

# construct a data package skeleton named "MyDataPackage" and pass
# in the Rmd file name with full path, and the name of the object(s) it
# creates.

pname <- basename(tempfile())
datapackage_skeleton(name=pname,
   path=tempdir(),
   force = TRUE,
   r_object_names = "tbl",
   code_files = f)

# call package_build to run the "foo.Rmd" processing and
# build a data package.
package_build(file.path(tempdir(), pname), install = FALSE)

# "install" the data package
devtools::load_all(file.path(tempdir(), pname))

# read the data version
data_version(pname)

# list the data sets in the package.
data(package = pname)

# The data objects are in the package source under "/data"
list.files(pattern="rda", path = file.path(tempdir(),pname,"data"), full = TRUE)

# The documentation that needs to be edited is in "/R"
list.files(pattern="R", path = file.path(tempdir(), pname,"R"), full = TRUE)
readLines(list.files(pattern="R", path = file.path(tempdir(),pname,"R"), full = TRUE))
# view the documentation with
?tbl
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skeleton.R
\name{datapackage_skeleton}
\alias{datapackage_skeleton}
\alias{datapackage.skeleton}
\title{Create a Data Package skeleton for use with DataPackageR.}
\usage{
datapackage_skeleton(
  name = NULL,
  path = ".",
  force = FALSE,
  code_files = character(),
  r_object_names = character(),
  raw_data_dir = character(),
  dependencies = character()
)

datapackage.skeleton(
  name = NULL,
  list = character(),
  environment = .GlobalEnv,
  path = ".",
  force = FALSE,
  code_files = character(),
  r_object_names = character()
)
}
\arguments{
\item{name}{\code{character} name of the package to create.}

\item{path}{A \code{character} path where the package is located. See \code{\link[utils]{package.skeleton}}}

\item{force}{\code{logical} Force the package skeleton to be recreated even if it exists. see \code{\link[utils]{package.skeleton}}}

\item{code_files}{Optional \code{character} vector of paths to Rmd files that process raw data
into R objects.}

\item{r_object_names}{\code{vector} of quoted r object names , tables, etc. created when the files in \code{code_files} are run.}

\item{raw_data_dir}{\code{character} pointing to a raw data directory. Will be moved with all its subdirectories to "inst/extdata"}

\item{dependencies}{\code{vector} of \code{character}, paths to R files that will be moved to "data-raw" but not included in the yaml config file. e.g., dependency scripts.}

\item{list}{Not used.}

\item{environment}{Not used.}
}
\description{
Creates a package skeleton directory structure for use with DataPackageR.
Adds the DataVersion string to DESCRIPTION, creates the DATADIGEST file, and the data-raw directory.
Updates the Read-and-delete-me file to reflect the additional necessary steps.
}
\note{
renamed \code{datapackage.skeleton()} to \code{datapackage_skeleton()}.
}
\examples{
if(rmarkdown::pandoc_available()){
f <- tempdir()
f <- file.path(f,"foo.Rmd")
con <- file(f)
writeLines("```{r}\n tbl = table(sample(1:10,1000,replace=TRUE)) \n```\n",con=con)
close(con)
pname <- basename(tempfile())
datapackage_skeleton(name = pname,
   path = tempdir(),
   force = TRUE,
   r_object_names = "tbl",
   code_files = f)
   }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processData.R
\name{DataPackageR}
\alias{DataPackageR}
\title{Process data generation code in 'data-raw'}
\usage{
DataPackageR(arg = NULL, deps = TRUE)
}
\arguments{
\item{arg}{\code{character} name of the package to build.}

\item{deps}{\code{logical} should scripts pass data objects to each other (default=TRUE)}
}
\value{
logical TRUE if successful, FALSE, if not.
}
\description{
Assumes .R files in 'data-raw' generate rda files to be stored in 'data'.
Sources datasets.R which can source other R files.
R files sourced by datasets.R must invoke \code{sys.source('myRfile.R',env=topenv())}.
Meant to be called before R CMD build.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use.R
\name{use_data_object}
\alias{use_data_object}
\title{Add a data object to a data package.}
\usage{
use_data_object(object_name = NULL)
}
\arguments{
\item{object_name}{Name of the data object. Should be created by a processing script in data-raw. \code{character} vector of length 1.}
}
\value{
invisibly returns TRUE for success.
}
\description{
The data object will be added to the yml configuration file.
}
\examples{
if(rmarkdown::pandoc_available()){
myfile <- tempfile()
file <- system.file("extdata", "tests", "extra.rmd",
                     package = "DataPackageR")
datapackage_skeleton(
  name = "datatest",
  path = tempdir(),
  code_files = file,
  force = TRUE,
  r_object_names = "data")
use_data_object(object_name = "newobject")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{package_build}
\alias{package_build}
\title{Pre-process, document and build a data package}
\usage{
package_build(
  packageName = NULL,
  vignettes = FALSE,
  log = INFO,
  deps = TRUE,
  install = FALSE,
  ...
)
}
\arguments{
\item{packageName}{\code{character} path to package source directory. Defaults to the current path when NULL.}

\item{vignettes}{\code{logical} specify whether to build vignettes. Default FALSE.}

\item{log}{log level \code{INFO,WARN,DEBUG,FATAL}}

\item{deps}{\code{logical} should we pass data objects into subsequent scripts? Default TRUE}

\item{install}{\code{logical} automatically install and load the package after building. (default TRUE)}

\item{...}{additional arguments passed to \code{install.packages} when \code{install=TRUE}.}
}
\description{
Combines the preprocessing, documentation, and build steps into one.
}
\examples{
if(rmarkdown::pandoc_available()){
f <- tempdir()
f <- file.path(f,"foo.Rmd")
con <- file(f)
writeLines("```{r}\n tbl = table(sample(1:10,1000,replace=TRUE)) \n```\n",con=con)
close(con)
pname <- basename(tempfile())
datapackage_skeleton(name=pname,
   path=tempdir(),
   force = TRUE,
   r_object_names = "tbl",
   code_files = f)

package_build(file.path(tempdir(),pname), install = FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processData.R
\name{project_data_path}
\alias{project_data_path}
\title{Get DataPackageR data path}
\usage{
project_data_path(file = NULL)
}
\arguments{
\item{file}{\code{character} or \code{NULL} (default).}
}
\value{
\code{character}
}
\description{
Get DataPackageR data path
}
\details{
Returns the path to the data package data subdirectory, or
constructs a path to a file in the data subdirectory from the
file argument.
}
\examples{
if(rmarkdown::pandoc_available()){
project_data_path( file = "data.rda" )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/environments.R
\name{datapackager_object_read}
\alias{datapackager_object_read}
\title{Read an object created in a previously run processing script.}
\usage{
datapackager_object_read(name)
}
\arguments{
\item{name}{\code{character} the name of the object. Must be a
name available in the configuration objects. Other objects are not saved.}
}
\value{
An R object.
}
\description{
Read an object created in a previously run processing script.
}
\details{
This function is only accessible within an R or Rmd file processed by DataPackageR.
It searches for an environment named \code{ENVS} within the current environment,
that holds the object with the given \code{name}. Such an environment is constructed and populated
with objects specified in the yaml \code{objects} property and passed along
to subsequent R and Rmd files as DataPackageR processes them in order.
}
\examples{
if(rmarkdown::pandoc_available()){
ENVS <- new.env() # ENVS would be in the environment
                 # where the data processing is run. It is
                 # handled automatically by the package.
assign("find_me", 100, ENVS) #This is done automatically by DataPackageR

find_me <- datapackager_object_read("find_me") # This would appear in an Rmd processed by
                                    # DataPackageR to access the object named "find_me" created
                                    # by a previous script. "find_me" would also need to
                                    # appear in the objects property of config.yml
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yamlR.R
\name{construct_yml_config}
\alias{construct_yml_config}
\title{Construct a datapackager.yml configuration}
\usage{
construct_yml_config(code = NULL, data = NULL, render_root = NULL)
}
\arguments{
\item{code}{A vector of filenames}

\item{data}{A vector of quoted object names}

\item{render_root}{The root directory where the package data processing code will be rendered.
Defaults to is set to a randomly generated named subdirectory of \code{tempdir()}.}
}
\value{
a datapackager.yml configuration represented as an R object
}
\description{
Constructs a datapackager.yml configuration object from a vector of file names and a vector of object names (all quoted).
Can be written to disk via \code{yml_write}.
\code{render_root} is set to a randomly generated named subdirectory of \code{tempdir()}.
}
\examples{
conf <- construct_yml_config(code = c('file1.rmd','file2.rmd'), data=c('object1','object2'))
tmp <- normalizePath(tempdir(), winslash = "/")
yml_write(conf,path=tmp)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ignore.R
\name{use_ignore}
\alias{use_ignore}
\title{Ignore specific files by git and R build.}
\usage{
use_ignore(file = NULL, path = NULL)
}
\arguments{
\item{file}{\code{character} File to ignore.}

\item{path}{\code{character} Path to the file.}
}
\value{
invisibly returns 0.
}
\description{
Ignore specific files by git and R build.
}
\examples{
datapackage_skeleton(name="test",path = tempdir())
use_ignore("foo", ".")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yamlR.R
\name{yml_find}
\alias{yml_find}
\alias{yml_add_files}
\alias{yml_disable_compile}
\alias{yml_enable_compile}
\alias{yml_add_objects}
\alias{yml_list_objects}
\alias{yml_list_files}
\alias{yml_remove_objects}
\alias{yml_remove_files}
\alias{yml_write}
\title{Edit DataPackageR yaml configuration}
\usage{
yml_find(path)

yml_add_files(config, filenames)

yml_disable_compile(config, filenames)

yml_enable_compile(config, filenames)

yml_add_objects(config, objects)

yml_list_objects(config)

yml_list_files(config)

yml_remove_objects(config, objects)

yml_remove_files(config, filenames)

yml_write(config, path = NULL)
}
\arguments{
\item{path}{Path to the data package source or path to write config file (for \code{yml_write})}

\item{config}{an R representation of the datapackager.yml config, returned by yml_find, or a path to the package root.}

\item{filenames}{A vector of filenames.}

\item{objects}{A vector of R object names.}
}
\value{
A yaml configuration structured as an R nested list.
}
\description{
Edit a yaml configuration file via an API.
}
\details{
Add, remove files and objects, enable or disable parsing of specific files,  list objects or files in a yaml config, or write a config back to a package.
}
\examples{
if(rmarkdown::pandoc_available()){
f <- tempdir()
f <- file.path(f,"foo.Rmd")
con <- file(f)
writeLines("```{r}\n tbl = table(sample(1:10,1000,replace=TRUE)) \n```\n",con=con)
close(con)
pname <- basename(tempfile())
datapackage_skeleton(name=pname,
   path = tempdir(),
   force = TRUE,
   r_object_names = "tbl",
   code_files = f)
yml <- yml_find(file.path(tempdir(),pname))
cat(yaml::as.yaml(yml))
yml <- yml_add_files(yml,"foo.Rmd")
yml_list_files(yml)
yml <- yml_disable_compile(yml,"foo.Rmd")
cat(yaml::as.yaml(yml))
yml <- yml_enable_compile(yml,"foo.Rmd")
cat(yaml::as.yaml(yml))
yml <- yml_add_objects(yml,"data1")
yml_list_objects(yml)
yml <- yml_remove_objects(yml,"data1")
yml <- yml_remove_files(yml,"foo.Rmd")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processData.R
\name{project_path}
\alias{project_path}
\title{Get DataPackageR Project Root Path}
\usage{
project_path(file = NULL)
}
\arguments{
\item{file}{\code{character} or \code{NULL} (default).}
}
\value{
\code{character}
}
\description{
Get DataPackageR Project Root Path
}
\details{
Returns the path to the data package project root, or
constructs a path to a file in the project root from the
file argument.
}
\examples{
if(rmarkdown::pandoc_available()){
project_path( file = "DESCRIPTION" )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataversion.r
\name{assert_data_version}
\alias{assert_data_version}
\title{Assert that a data version in a data package matches an expectation.}
\usage{
assert_data_version(
  data_package_name = NULL,
  version_string = NULL,
  acceptable = "equal",
  ...
)
}
\arguments{
\item{data_package_name}{\code{character} Name of the package.}

\item{version_string}{\code{character} Version string in "x.y.z" format.}

\item{acceptable}{\code{character} one of "equal", "equal_or_greater", describing what version match is acceptable.}

\item{...}{additional arguments passed to data_version (such as lib.loc)}
}
\value{
invisible \code{logical} TRUE if success, otherwise stop on mismatch.
}
\description{
Assert that a data version in a data package matches an expectation.
}
\details{
Tests the DataVersion string in \code{data_package_name} against \code{version_string} testing the major, minor and revision portion.

Tests "data_package_name version equal version_string" or "data_package_name version equal_or_greater version_string".
}
\examples{
if(rmarkdown::pandoc_available()){
f <- tempdir()
f <- file.path(f, "foo.Rmd")
con <- file(f)
writeLines("```{r}\n tbl = table(sample(1:10,1000,replace=TRUE)) \n```\n",con = con)
close(con)
pname <- basename(tempfile())
datapackage_skeleton(name = pname,
   path=tempdir(),
   force = TRUE,
   r_object_names = "tbl",
   code_files = f)
package_build(file.path(tempdir(),pname), install = FALSE)

devtools::load_all(file.path(tempdir(),pname))

assert_data_version(data_package_name = pname,version_string = "0.1.0",acceptable = "equal")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processData.R
\name{document}
\alias{document}
\title{Build documentation for a data package using DataPackageR.}
\usage{
document(path = ".", install = TRUE, ...)
}
\arguments{
\item{path}{\code{character} the path to the data package source root.}

\item{install}{\code{logical} install and reload the package. (default TRUE)}

\item{...}{additional arguments to \code{install}}
}
\description{
Build documentation for a data package using DataPackageR.
}
\examples{
# A simple Rmd file that creates one data object
# named "tbl".
if(rmarkdown::pandoc_available()){
f <- tempdir()
f <- file.path(f,"foo.Rmd")
con <- file(f)
writeLines("```{r}\n tbl = table(sample(1:10,100,replace=TRUE)) \n```\n",con=con)
close(con)

# construct a data package skeleton named "MyDataPackage" and pass
# in the Rmd file name with full path, and the name of the object(s) it
# creates.

pname <- basename(tempfile())
datapackage_skeleton(name=pname,
   path=tempdir(),
   force = TRUE,
   r_object_names = "tbl",
   code_files = f)

# call package_build to run the "foo.Rmd" processing and
# build a data package.
package_build(file.path(tempdir(), pname), install = FALSE)
document(path = file.path(tempdir(), pname), install=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processData.R
\name{project_extdata_path}
\alias{project_extdata_path}
\title{Get DataPackageR extdata path}
\usage{
project_extdata_path(file = NULL)
}
\arguments{
\item{file}{\code{character} or \code{NULL} (default).}
}
\value{
\code{character}
}
\description{
Get DataPackageR extdata path
}
\details{
Returns the path to the data package extdata subdirectory, or
constructs a path to a file in the extdata subdirectory from the
file argument.
}
\examples{
if(rmarkdown::pandoc_available()){
project_extdata_path(file = "mydata.csv")
}
}
