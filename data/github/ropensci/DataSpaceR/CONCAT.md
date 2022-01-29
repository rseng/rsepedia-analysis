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
(http:contributor-covenant.org), version 1.0.0, available at
http://contributor-covenant.org/version/1/0/0/

<!-- README.md is generated from README.Rmd. Please edit that file -->

# DataSpaceR <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/DataSpaceR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/DataSpaceR/actions)
[![codecov](https://codecov.io/gh/ropensci/DataSpaceR/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/DataSpaceR/branch/main)
[![CRAN
Status](https://www.r-pkg.org/badges/version/DataSpaceR)](https://cran.r-project.org/package=DataSpaceR)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/DataSpaceR)](https://www.rpackages.io/package/DataSpaceR)
[![monthly](https://cranlogs.r-pkg.org/badges/DataSpaceR)](https://www.rpackages.io/package/DataSpaceR)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![](https://badges.ropensci.org/261_status.svg)](https://github.com/ropensci/software-review/issues/261)
<!-- badges: end -->

`DataSpaceR` is an R interface to [the CAVD
DataSpace](https://dataspace.cavd.org), a data sharing and discovery
tool that facilitates exploration of HIV immunological data from
pre-clinical and clinical HIV vaccine studies.

This package is intended for use by immunologists, bioinformaticians, and
statisticians in HIV vaccine research, or anyone interested in the
analysis of HIV immunological data across assays, studies, and time.

This package simplifies access to the database by taking advantage of
the standardization of the database to hide all the
[Rlabkey](https://cran.r-project.org/package=Rlabkey) specific code away
from the user, and it allows the users to access the study-specific
datasets via [an object-oriented
paradigm](https://cran.r-project.org/package=R6/readme/README.html).

## Examples & Documentation

For more detailed examples and detailed documentation, see [the
introductory
vignette](https://docs.ropensci.org/DataSpaceR/articles/DataSpaceR.html)
and [the pkgdown site](https://docs.ropensci.org/DataSpaceR/).

For a quick guide of how to use the API, see our [cheat sheet](https://dataspace.cavd.org/_webdav/static/%40files/documents/dataspacer_cheat_sheet.pdf) .

## Installation

Install from CRAN:

``` r
install.packages("DataSpaceR")
```

You can install the latest development version from
[GitHub](https://github.com/ropensci/DataSpaceR) with
[devtools](https://cran.r-project.org/package=devtools):

``` r
# install.packages("devtools")
devtools::install_github("ropensci/DataSpaceR")
```

## Register and set DataSpace credential

The database is accessed with the user’s credentials. A netrc file
storing login and password information is ***required***.

1.  [Create an account](https://dataspace.cavd.org/) and read the terms
    of use
2.  On your R console, create a netrc file using a function from
    `DataSpaceR`:

``` r
library(DataSpaceR)
writeNetrc(
  login = "yourEmail@address.com", 
  password = "yourSecretPassword",
  netrcFile = "/your/home/directory/.netrc" # use getNetrcPath() to get the default path 
)
```

This will create a netrc file in your home directory.

***Alternatively***, you can manually create a netrc file in the
computer running R.

-   On Windows, this file should be named `_netrc`
-   On UNIX, it should be named `.netrc`
-   The file should be located in the user’s home directory, and the
    permissions on the file should be unreadable for everybody except
    the owner
-   To determine home directory, run `Sys.getenv("HOME")` in R

The following three lines must be included in the `.netrc` or `_netrc`
file either separated by white space (spaces, tabs, or newlines) or
commas. Multiple such blocks can exist in one file.

    machine dataspace.cavd.org
    login myuser@domain.com
    password supersecretpassword

See
[here](https://www.labkey.org/wiki/home/Documentation/page.view?name=netrc)
for more information about `netrc`.

## Usage

The general idea is that the user:

1.  creates an instance of `DataSpaceConnection` class via `connectDS`
2.  browses available studies and groups in the instance via
    `availableStudies` and `availableGroups`
3.  creates a connection to a specific study via `getStudy` or a group
    via `getGroup`
4.  retrieves datasets by name via `getDataset`

### for example:

``` r
library(DataSpaceR)
#> By exporting data from the CAVD DataSpace, you agree to be bound by the Terms of Use available on the CAVD DataSpace sign-in page at https://dataspace.cavd.org

con <- connectDS()
con
#> <DataSpaceConnection>
#>   URL: https://dataspace.cavd.org
#>   User: jkim2345@scharp.org
#>   Available studies: 273
#>     - 77 studies with data
#>     - 5049 subjects
#>     - 423195 data points
#>   Available groups: 6
#>   Available publications: 1530
#>     - 12 publications with data
```

`connectDS()` will create a connection to DataSpace.

### available studies can be listed by `availableStudies` field

``` r
knitr::kable(head(con$availableStudies))
```

| study\_name | short\_name                    | title                                                                                                                                                                                                                                                              | type               | status   | stage            | species            | start\_date | strategy                             | network | data\_availability | ni\_data\_availability |
|:------------|:-------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------|:---------|:-----------------|:-------------------|:------------|:-------------------------------------|:--------|:-------------------|:-----------------------|
| cor01       | NA                             | The correlate of risk targeted intervention study (CORTIS): A randomized, partially-blinded, clinical trial of isoniazid and rifapentine (3HP) therapy to prevent pulmonary tuberculosis in high-risk individuals identified by a transcriptomic correlate of risk | Phase III          | Inactive | Assays Completed | Human              | NA          | NA                                   | GH-VAP  | NA                 | NA                     |
| cvd232      | Parks\_RV\_232                 | ​Limiting Dose Vaginal SIVmac239 Challenge of RhCMV-SIV vaccinated Indian rhesus macaques.                                                                                                                                                                         | Pre-Clinical NHP   | Inactive | Assays Completed | Rhesus macaque     | 2009-11-24  | Vector vaccines (viral or bacterial) | CAVD    | NA                 | NA                     |
| cvd234      | Zolla-Pazner\_Mab\_test1 Study | Zolla-Pazner\_Mab\_Test1                                                                                                                                                                                                                                           | Antibody Screening | Inactive | Assays Completed | Non-Organism Study | 2009-02-03  | Prophylactic neutralizing Ab         | CAVD    | NA                 | NA                     |
| cvd235      | mAbs potency                   | Weiss mAbs potency                                                                                                                                                                                                                                                 | Antibody Screening | Inactive | Assays Completed | Non-Organism Study | 2008-08-21  | Prophylactic neutralizing Ab         | CAVD    | NA                 | NA                     |
| cvd236      | neutralization assays          | neutralization assays                                                                                                                                                                                                                                              | Antibody Screening | Active   | In Progress      | Non-Organism Study | 2009-02-03  | Prophylactic neutralizing Ab         | CAVD    | NA                 | NA                     |
| cvd238      | Gallo\_PA\_238                 | HIV-1 neutralization responses in chronically infected individuals                                                                                                                                                                                                 | Antibody Screening | Inactive | Assays Completed | Non-Organism Study | 2009-01-08  | Prophylactic neutralizing Ab         | CAVD    | NA                 | NA                     |

### available groups can be listed by `availableGroups` field

``` r
knitr::kable(con$availableGroups)
```

| group\_id | label                              | original\_label                    | description                                                                                                               | created\_by | shared |   n | studies                        |
|----------:|:-----------------------------------|:-----------------------------------|:--------------------------------------------------------------------------------------------------------------------------|:------------|:-------|----:|:-------------------------------|
|       216 | mice                               | mice                               | NA                                                                                                                        | readjk      | FALSE  |  75 | cvd468, cvd483, cvd316, cvd331 |
|       217 | CAVD 242                           | CAVD 242                           | This is a fake group for CAVD 242                                                                                         | readjk      | FALSE  |  30 | cvd242                         |
|       220 | NYVAC durability comparison        | NYVAC\_durability                  | Compare durability in 4 NHP studies using NYVAC-C (vP2010) and NYVAC-KC-gp140 (ZM96) products.                            | ehenrich    | TRUE   |  78 | cvd281, cvd434, cvd259, cvd277 |
|       224 | cvd338                             | cvd338                             | NA                                                                                                                        | readjk      | FALSE  |  36 | cvd338                         |
|       228 | HVTN 505 case control subjects     | HVTN 505 case control subjects     | Participants from HVTN 505 included in the case-control analysis                                                          | drienna     | TRUE   | 189 | vtn505                         |
|       230 | HVTN 505 polyfunctionality vs BAMA | HVTN 505 polyfunctionality vs BAMA | Compares ICS polyfunctionality (CD8+, Any Env) to BAMA mfi-delta (single Env antigen) in the HVTN 505 case control cohort | drienna     | TRUE   | 170 | vtn505                         |

***Note***: A group is a curated collection of participants from
filtering of treatments, products, studies, or species, and it is
created in [the DataSpace
App](https://dataspace.cavd.org/cds/CAVD/app.view).

Check out [the reference
page](https://docs.ropensci.org/DataSpaceR/reference/DataSpaceConnection.html)
of `DataSpaceConnection` for all available fields and methods.

### create an instance of `cvd408`

``` r
cvd408 <- con$getStudy("cvd408")
cvd408
#> <DataSpaceStudy>
#>   Study: cvd408
#>   URL: https://dataspace.cavd.org/CAVD/cvd408
#>   Available datasets:
#>     - Binding Ab multiplex assay
#>     - Demographics
#>     - Intracellular Cytokine Staining
#>     - Neutralizing antibody
#>   Available non-integrated datasets:
class(cvd408)
#> [1] "DataSpaceStudy" "R6"
```

### available datasets can be listed by `availableDatasets` field

``` r
knitr::kable(cvd408$availableDatasets)
```

| name         | label                           |    n | integrated |
|:-------------|:--------------------------------|-----:|:-----------|
| BAMA         | Binding Ab multiplex assay      | 1080 | TRUE       |
| Demographics | Demographics                    |   20 | TRUE       |
| ICS          | Intracellular Cytokine Staining | 3720 | TRUE       |
| NAb          | Neutralizing antibody           |  540 | TRUE       |

which will print names of available datasets.

### Neutralizing Antibody dataset (`NAb`) can be retrieved by:

``` r
NAb <- cvd408$getDataset("NAb")
dim(NAb)
#> [1] 540  33
colnames(NAb)
#>  [1] "participant_id"      "participant_visit"   "visit_day"          
#>  [4] "assay_identifier"    "summary_level"       "specimen_type"      
#>  [7] "antigen"             "antigen_type"        "virus"              
#> [10] "virus_type"          "virus_insert_name"   "clade"              
#> [13] "neutralization_tier" "tier_clade_virus"    "target_cell"        
#> [16] "initial_dilution"    "titer_ic50"          "titer_ic80"         
#> [19] "response_call"       "nab_lab_source_key"  "lab_code"           
#> [22] "exp_assayid"         "titer_id50"          "titer_id80"         
#> [25] "nab_response_id50"   "nab_response_id80"   "slope"              
#> [28] "vaccine_matched"     "study_prot"          "virus_full_name"    
#> [31] "virus_species"       "virus_host_cell"     "virus_backbone"
```

Check out [the reference
page](https://docs.ropensci.org/DataSpaceR/reference/DataSpaceStudy.html)
of `DataSpaceStudy` for all available fields and methods.

***Note***: The package uses a
[R6](https://cran.r-project.org/package=R6) class to represent the
connection to a study and get around some of R’s copy-on-change
behavior.

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/DataSpaceR/issues).
-   License: GPL-3
-   Get citation information for `DataSpaceR` in R doing
    `citation(package = 'DataSpaceR')`
-   Please note that this project is released with a [Contributor Code
    of
    Conduct](https://github.com/ropensci/DataSpaceR/blob/main/CONDUCT.md).
    By participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# DataSpaceR 0.7.5

* Fixed bug concerning missing cookie session name for app reports which will update this package to work when run on CDS as in a report.
* Show NI data in available datasets active binding which will allow users to look for NI (non-integrated) data in DSR.
* Set active binding names to snake case which is in effort to standardize fields names used in the API.
* Add all authors to publication table which will allow users to find publication by authors who were not listed as the primary author.

# DataSpaceR 0.7.4

* Added `virusMetadata` field to `DataSpaceConnection` which shows virus metadata for viruses used in NAb assays (thanks @helenmiller16 #26)
* Added `availablePublications` field to `DataSpaceConnection` which summarizes all available publications  (thanks @helenmiller16 #27)
* Added `downloadPublicationData` method to `DataSpaceConnection` which will download available publication data for a specified publication (thanks @helenmiller16 #27)
* Added a vignette: `Accessing Publication Data` (thanks @helenmiller16 #27)
* Made `PKMAb` dataset available to retrieve in `DataSpaceStudy` (thanks @helenmiller16 and @jmtaylor-fhcrc #27)
* Relaxed some assumptions when pulling non-integrated data, which allows users to pull non-integrated mab data, like in cvd812. (thanks @helenmiller16 #28)
* Fixed a test for `mabGridSummary$geometric_mean_curve_ic50` calculation to reflect updated data. (thanks @helenmiller16 #28)
* Updated documentation using the R6 documentation syntax in roxygen

# DataSpaceR 0.7.3

* Remove `tolower` in functions that check study names. (thanks @jmtaylor-fhcrc #23)
* Update DataSpaceStudy methods for non-integrated datasets. (thanks @helenmiller16 #24)

# DataSpaceR 0.7.2

* Fix broken and invalid URLs for CRAN submission

# DataSpaceR 0.7.1

* Fix broken and invalid URLs for CRAN submission

# DataSpaceR 0.7.0

* Added `DataSpaceMab` class and several methods and fields in `DataSpaceConnection` to allow access to monoclonal antibody data (thanks @jmtaylor-fhcrc #14, #19, #21)

# DataSpaceR 0.6.3

* Prepare for CRAN submission

# DataSpaceR 0.6.2

* Write the netrc in temporary directory as default in `writeNetrc`. (CRAN requirement)

# DataSpaceR 0.6.1

* Remove LICENSE file for CRAN submission

# DataSpaceR 0.6.0

* rOpenSci submission: https://github.com/ropensci/software-review/issues/261

# DataSpaceR 0.5.2

* Modified `getDataset` method in `DataSpaceStudy` to take `mergeExtra` as an argument in order to merge extra information (demographics and treatment arm). #5

# DataSpaceR 0.5.1

* Adjusted the package to use the latest version of Rlabkey (v2.2) and httr packages.

# DataSpaceR 0.5.0

* Modified `DataSpaceConnection` and `DataSpaceStudy` classes to convert data.frame objects to data.table.
* Added a package startup message on the terms of use.
* Created `getGroup` method in `DataSpaceConnection` class and deprecated `groupId` in `getStudy` method.

# DataSpaceR 0.4.2

* Added `refresh` method for both connection and study classes.
* Added `studyInfo` field in study class.

# DataSpaceR 0.4.1

* Included additional columns (`short_name`, `type`, `status`, `stage`, `species`, `start_date`, `strategy`) in `availableStudies`.
* Updated the introductory vignette.

# DataSpaceR 0.4.0

* Added `availableGroups` field to `DataSpaceConnection` class.
* Modified `getStudy` method in `DataSpaceConnection` to take `groupId` as an argument in order to create a `DataSpaceStudy` object for a particular group.

# DataSpaceR 0.3.0

* Updated `DataSpaceConnection` class and `connectDS` constructor to be not study-specific.

# DataSpaceR 0.2.1

* Modified `getUserEmail()` to get email from netrc file.
* Updated license to GPL-3.

# DataSpaceR 0.2.0

* Added `getVariableInfo` method to `DataSpaceConnection` class.
* Changed variable names in `getAvailableDatasets` to lowercase.
* Added `treatmentArm` field to `DataSpaceConnection` class.
* Removed `getAvailableDatasets` method (now a private method).
* Changed the default connection from staging (`dataspace-staging.cavd.org`) to production (`dataspace.cavd.org`).
* Added an option to connect to the staging server.
* Renamed `write_netrc` and `check_netrc` to `writeNetrc` and `checkNetrc`.

# DataSpaceR 0.1.0

* Initialized the package with `connectDS()` and other basic functions.
* Added a test framework on `tests` to test the package.
* Added `.travis.yml` to build and check the package automatically. (not used yet)
* Added `codecov.yml` to track code coverage. (not used yet)
* Added `README.Rmd` to introduce the package.
* Added `NEWS.md` to track changes to the package.
* Added `_pkgdown.yml` and `/docs` to build a package website.
* Added a introductory vignette called `Intro_to_DataSpaceR.Rmd`.
## Release summary

* This is a new release.
* Jason Taylor (jmtaylor@fredhutch.org) is the new maintainer.

## Test environments
* local R installation, R 4.1.1
* ubuntu 18.04 and macOS (on GitHub Actions), R 4.1.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Jason Taylor <jmtaylor@fredhutch.org>'

  Jason Taylor <jmtaylor@fredhutch.org>
New maintainer:
Old maintainer(s):
  Ju Yeong Kim <jkim2345@fredhutch.org>
```
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
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
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

Please note that the DataSpaceR project is released with a
[Contributor Code of Conduct](CONDUCT.md). By contributing to this
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
