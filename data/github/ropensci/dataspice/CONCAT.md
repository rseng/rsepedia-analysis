# dataspice

![CRAN Version](https://www.r-pkg.org/badges/version/dataspice)
![CI](https://github.com/ropensci/dataspice/workflows/R-CMD-check/badge.svg)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/dataspice/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/dataspice?branch=main)
[![](https://badges.ropensci.org/426_status.svg)](https://github.com/ropensci/software-review/issues/426)

The goal of `dataspice` is to make it easier for researchers to create
basic, lightweight, and concise metadata files for their datasets by
editing the kind of files they’re probably most familiar with: CSVs. To
spice up their data with a dash of metadata. These metadata files can
then be used to:

-   Make useful information available during analysis.
-   Create a helpful dataset README webpage for your data similar to how
    [pkgdown](https://pkgdown.r-lib.org/) creates websites for R
    packages.
-   Produce more complex metadata formats for richer description of your
    datasets and to aid dataset discovery.

Metadata fields are based on
[Schema.org/Dataset](https://schema.org/Dataset) and other [metadata
standards](#resources) and represent a lowest common denominator which
means converting between formats should be relatively straightforward.

## Example

An basic example repository for demonstrating what using `dataspice`
might look like can be found at
[https://github.com/amoeba/dataspice-example](https://github.com/amoeba/dataspice-example/).
From there, you can also check out a preview of the HTML `dataspice`
generates at
[https://amoeba.github.io/dataspice-example](https://amoeba.github.io/dataspice-example/)
and how Google sees it at
<https://search.google.com/test/rich-results?url=https%3A%2F%2Famoeba.github.io%2Fdataspice-example%2F>.

A much more detailed example has been created by [Anna
Krystalli](https://annakrystalli.me) at
<https://annakrystalli.me/dataspice-tutorial/> ([GitHub
repo](https://github.com/annakrystalli/dataspice-tutorial)).

## Installation

You can install the latest version from
[CRAN](https://cran.r-project.org):

``` r
install.packages("dataspice")
```

## Workflow

``` r
create_spice()
# Then fill in template CSV files, more on this below
write_spice()
build_site() # Optional
```

![diagram showing a workflow for using
dataspice](man/figures/dataspice_workflow.png)

### Create spice

`create_spice()` creates template metadata spreadsheets in a folder (by
default created in the `data` folder in the current working directory).

The template files are:

-   **biblio.csv** - for title, abstract, spatial and temporal coverage,
    etc.
-   **creators.csv** - for data authors
-   **attributes.csv** - explains each of the variables in the dataset
-   **access.csv** - for files, file types, and download URLs (if
    appropriate)

### Fill in templates

The user needs to fill in the details of the four template files. These
csv files can be directly modified, or they can be edited using either
the associated helper function and/or
[Shiny](https://shiny.rstudio.com/) app.

#### Helper functions

-   `prep_attributes()` populates the **`fileName`** and
    **`variableName`** columns of the `attributes.csv` file using the
    header row of the data files.

-   `prep_access()` populates the **`fileName`**, **`name`** and
    **`encodingFormat`** columns of the `access.csv` file from the files
    in the folder containing the data.

To see an example of how `prep_attributes()` works, load the data files
that ship with the package:

``` r
data_files <- list.files(system.file("example-dataset/", package = "dataspice"),
  pattern = ".csv",
  full.names = TRUE
)
```

This function assumes that the metadata templates are in a folder called
`metadata` within a `data` folder.

``` r
attributes_path <- file.path("data", "metadata", "attributes.csv")
```

Using `purrr::map()`, this function can be applied over multiple files
to populate the header names

``` r
data_files %>%
  purrr::map(~ prep_attributes(.x, attributes_path),
    attributes_path = attributes_path
  )
```

The output of `prep_attributes()` has the first two columns filled out:

<table>
<thead>
<tr>
<th style="text-align:left;">
fileName
</th>
<th style="text-align:left;">
variableName
</th>
<th style="text-align:left;">
description
</th>
<th style="text-align:left;">
unitText
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BroodTables.csv
</td>
<td style="text-align:left;">
Stock.ID
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
BroodTables.csv
</td>
<td style="text-align:left;">
Species
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
BroodTables.csv
</td>
<td style="text-align:left;">
Stock
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
BroodTables.csv
</td>
<td style="text-align:left;">
Ocean.Region
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
BroodTables.csv
</td>
<td style="text-align:left;">
Region
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
BroodTables.csv
</td>
<td style="text-align:left;">
Sub.Region
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
</tbody>
</table>

#### Shiny helper apps

Each of the metadata templates can be edited interactively using a
[Shiny](https://shiny.rstudio.com/) app:

-   `edit_attributes()` opens a Shiny app that can be used to edit
    `attributes.csv`. The Shiny app displays the current `attributes`
    table and lets the user fill in an informative description and units
    (e.g. meters, hectares, etc.) for each variable.
-   `edit_access()` opens an editable version of `access.csv`
-   `edit_creators()` opens an editable version of `creators.csv`
-   `edit_biblio()` opens an editable version of `biblio.csv`

![edit\_attributes Shiny app](man/figures/edit_attributes.png)

Remember to click on **Save** when finished editing.

#### Completed metadata files

The first few rows of the completed metadata tables in this example will
look like this:

`access.csv` has one row for each file

| fileName        | name            | contentUrl | encodingFormat |
|:----------------|:----------------|:-----------|:---------------|
| StockInfo.csv   | StockInfo.csv   | NA         | CSV            |
| BroodTables.csv | BroodTables.csv | NA         | CSV            |
| SourceInfo.csv  | SourceInfo.csv  | NA         | CSV            |

`attributes.csv` has one row for each variable in each file

| fileName        | variableName | description                                      | unitText |
|:----------------|:-------------|:-------------------------------------------------|:---------|
| BroodTables.csv | Stock.ID     | Unique stock identifier                          | NA       |
| BroodTables.csv | Species      | species of stock                                 | NA       |
| BroodTables.csv | Stock        | Stock name, generally river where stock is found | NA       |
| BroodTables.csv | Ocean.Region | Ocean region                                     | NA       |
| BroodTables.csv | Region       | Region of stock                                  | NA       |
| BroodTables.csv | Sub.Region   | Sub.Region of stock                              | NA       |

`biblio.csv` is one row containing descriptors including spatial and
temporal coverage

| title                                                                 | description                                                                                                                                                                                            | datePublished       | citation | keywords                   | license | funder | geographicDescription | northBoundCoord | eastBoundCoord | southBoundCoord | westBoundCoord | wktString | startDate           | endDate             |
|:----------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------|:---------|:---------------------------|:--------|:-------|:----------------------|----------------:|---------------:|----------------:|---------------:|:----------|:--------------------|:--------------------|
| Compiled annual statewide Alaskan salmon escapement counts, 1921-2017 | The number of mature salmon migrating from the marine environment to freshwater streams is defined as escapement. Escapement data are the enumeration of these migrating fish as they pass upstream, … | 2018-02-12 08:00:00 | NA       | salmon, alaska, escapement | NA      | NA     | NA                    |              78 |           -131 |              47 |           -171 | NA        | 1921-01-01 08:00:00 | 2017-01-01 08:00:00 |

`creators.csv` has one row for each of the dataset authors

| id  | name           | affiliation                                           | email                      |
|:----|:---------------|:------------------------------------------------------|:---------------------------|
| NA  | Jeanette Clark | National Center for Ecological Analysis and Synthesis | <jclark@nceas.ucsb.edu>    |
| NA  | Rich,Brenner   | Alaska Department of Fish and Game                    | richard.brenner.alaska.gov |

### Save JSON-LD file

`write_spice()` generates a json-ld file (“linked data”) to aid in
[dataset
discovery](https://developers.google.com/search/docs/data-types/dataset),
creation of more extensive metadata
(e.g. [EML](https://eml.ecoinformatics.org)), and creating a website.

Here’s a view of the `dataspice.json` file of the example data:

![listviewer pack output showing an example dataspice JSON
file](man/figures/listviewer.png)

### Build website

-   `build_site()` creates a bare-bones `index.html` file in the
    repository `docs` folder with a simple view of the dataset with the
    metadata and an interactive map. For example, this
    [repository](https://github.com/amoeba/dataspice-example/) results
    in this [website](https://amoeba.github.io/dataspice-example/)

![dataspice-website](man/figures/website_example.png)

### Convert to EML

The metadata fields `dataspice` uses are based largely on their
compatibility with terms from [Schema.org](https://schema.org). However,
`dataspice` metadata can be converted to Ecological Metadata Language
(EML), a much richer schema. The conversion isn’t perfect but
`dataspice` will do its best to convert your `dataspice` metadata to
EML:

``` r
library(dataspice)

# Load an example dataspice JSON that comes installed with the package
spice <- system.file(
  "examples", "annual-escapement.json",
  package = "dataspice"
)

# Convert it to EML
eml_doc <- spice_to_eml(spice)
#> Warning: variableMeasured not crosswalked to EML because we don't have enough
#> information. Use `crosswalk_variables` to create the start of an EML attributes
#> table. See ?crosswalk_variables for help.
#> You might want to run EML::eml_validate on the result at this point and fix what validations errors are produced. You will commonly need to set `packageId`, `system`, and provide `attributeList` elements for each `dataTable`.
```

You may receive warnings depending on which `dataspice` fields you
filled in and this process will very likely produce an invalid EML
record which is totally fine:

``` r
library(EML)
#> 
#> Attaching package: 'EML'
#> The following object is masked from 'package:magrittr':
#> 
#>     set_attributes

eml_validate(eml_doc)
#> [1] FALSE
#> attr(,"errors")
#> [1] "Element '{https://eml.ecoinformatics.org/eml-2.2.0}eml': The attribute 'packageId' is required but missing."                                  
#> [2] "Element '{https://eml.ecoinformatics.org/eml-2.2.0}eml': The attribute 'system' is required but missing."                                     
#> [3] "Element 'dataTable': Missing child element(s). Expected is one of ( physical, coverage, methods, additionalInfo, annotation, attributeList )."
#> [4] "Element 'dataTable': Missing child element(s). Expected is one of ( physical, coverage, methods, additionalInfo, annotation, attributeList )."
#> [5] "Element 'dataTable': Missing child element(s). Expected is one of ( physical, coverage, methods, additionalInfo, annotation, attributeList )."
```

This is because some fields in `dataspice` store information in
different structures and because EML requires many fields that
`dataspice` doesn’t have fields for. At this point, you should look over
the validation errors produced by `EML::eml_validate` and fix those.
Note that this will likely require familiarity with the [EML
Schema](https://eml.ecoinformatics.org/) and the [EML
package](https://github.com/ropensci/eml).

Once you’re done, you can write out an EML XML file:

``` r
out_path <- tempfile()
write_eml(eml_doc, out_path)
#> NULL
```

### Convert from EML

Like converting `dataspice` to EML, we can convert an existing EML
record to a set of `dataspice` metadata tables which we can then work
from within `dataspice`:

``` r
library(EML)

eml_path <- system.file("example-dataset/broodTable_metadata.xml", package = "dataspice")
eml <- read_eml(eml_path)
```

``` r
# Creates four CSVs files in the `data/metadata` directory
my_spice <- eml_to_spice(eml, "data/metadata")
```

## Resources

A few existing tools & data standards to help users in specific domains:

-   [Darwin Core](http://rs.tdwg.org/dwc/)
-   [Ecological Metadata
    Language](https://knb.ecoinformatics.org/#external//emlparser/docs/index.html)
    (EML) (& [`EML`](https://github.com/ropensci/EML))
-   [ISO 19115](https://www.iso.org/standard/53798.html) - Geographic
    Information Metadata
-   [ISO 19139](https://www.iso.org/standard/32557.html) - Geographic
    Info Metadata XML schema
-   [Minimum Information for Biological and Biomedical
    Investigations](https://fairsharing.org/collection/MIBBI) (MIBBI)

…And others indexed in [Fairsharing.org](https://fairsharing.org) & the
[RDA metadata
directory](http://rd-alliance.github.io/metadata-directory/standards/).

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## Contributors

This package was developed at rOpenSci’s 2018 unconf by (in alphabetical
order):

-   [Carl Boettiger](https://github.com/cboettig)
-   [Scott Chamberlain](https://github.com/sckott)
-   [Auriel Fournier](https://github.com/aurielfournier)
-   [Kelly Hondula](https://github.com/khondula)
-   [Anna Krystalli](https://github.com/annakrystalli)
-   [Bryce Mecum](https://github.com/amoeba)
-   [Maëlle Salmon](https://github.com/maelle)
-   [Kate Webbink](https://github.com/magpiedin)
-   [Kara Woo](https://github.com/karawoo)
-   [Irene Steves](https://github.com/isteves)
# NEWS

## dataspice 1.1.0

This release is for the version of the package after going through [rOpenSci Onboarding](https://github.com/ropensci/software-review/issues/426).
No breaking user-facing changes were made.

NEW

- `dataspice` helpfully now has a help page: `?dataspice`

CHANGED

- Test coverage dramatically increased
- Various documentation (README, vignettes, help pages) now much more consistent

## dataspice 1.0.0

Initial CRAN release
# cran-comments

This is a minor release submission and is largely made of changes to documentation and tests.

## Test environments

On GitHub Actions:

- windows (release)
- macOS (release)
- ubuntu-20.04 (release)
- ubuntu-20.04 (devel)

## R CMD check

0 errors, 0 warnings, 0 notes
---
output:
  md_document:
    variant: gfm
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
library(knitr)
library(kableExtra)
library(magrittr)
library(stringr)
```
# dataspice

![CRAN Version](https://www.r-pkg.org/badges/version/dataspice)
![CI](https://github.com/ropensci/dataspice/workflows/R-CMD-check/badge.svg)
[![Codecov test coverage](https://codecov.io/gh/ropensci/dataspice/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/dataspice?branch=main)
[![](https://badges.ropensci.org/426_status.svg)](https://github.com/ropensci/software-review/issues/426)

The goal of `dataspice` is to make it easier for researchers to create basic, lightweight, and concise metadata files for their datasets by editing the kind of files they're probably most familiar with: CSVs. To spice up their data with a dash of metadata. These metadata files can then be used to:

- Make useful information available during analysis.
- Create a helpful dataset README webpage for your data similar to how [pkgdown](https://pkgdown.r-lib.org/) creates websites for R packages.
- Produce more complex metadata formats for richer description of your datasets and to aid dataset discovery.

Metadata fields are based on [Schema.org/Dataset](https://schema.org/Dataset) and other [metadata standards](#resources) and represent a lowest common denominator which means converting between formats should be relatively straightforward.

## Example

An basic example repository for demonstrating what using `dataspice` might look like can be found at [https://github.com/amoeba/dataspice-example](https://github.com/amoeba/dataspice-example/). From there, you can also check out a preview of the HTML `dataspice` generates at [https://amoeba.github.io/dataspice-example](https://amoeba.github.io/dataspice-example/) and how Google sees it at [https://search.google.com/test/rich-results?url=https%3A%2F%2Famoeba.github.io%2Fdataspice-example%2F
](https://search.google.com/test/rich-results?url=https%3A%2F%2Famoeba.github.io%2Fdataspice-example%2F).

A much more detailed example has been created by [Anna Krystalli](https://annakrystalli.me) at [https://annakrystalli.me/dataspice-tutorial/](https://annakrystalli.me/dataspice-tutorial/) ([GitHub repo](https://github.com/annakrystalli/dataspice-tutorial)).

## Installation

You can install the latest version from [CRAN](https://cran.r-project.org):

```r
install.packages("dataspice")
```

## Workflow

```{r example, eval=FALSE}
create_spice()
# Then fill in template CSV files, more on this below
write_spice()
build_site() # Optional
```

![diagram showing a workflow for using dataspice](man/figures/dataspice_workflow.png)

### Create spice

`create_spice()` creates template metadata spreadsheets in a folder (by default created in the `data` folder in the current working directory).

The template files are:

* **biblio.csv** - for title, abstract, spatial and temporal coverage, etc.
* **creators.csv** - for data authors
* **attributes.csv** - explains each of the variables in the dataset
* **access.csv** - for files, file types, and download URLs (if appropriate)

### Fill in templates

The user needs to fill in the details of the four template files. These csv files can be directly modified, or they can be edited using either the associated helper function and/or [Shiny](https://shiny.rstudio.com/) app.

#### Helper functions

* `prep_attributes()` populates the **`fileName`** and **`variableName`** columns of the `attributes.csv` file using the header row of the data files.

* `prep_access()` populates the **`fileName`**, **`name`** and **`encodingFormat`** columns of the `access.csv` file from the files in the folder containing the data.

To see an example of how `prep_attributes()` works, load the data files that ship with the package:

```{r, eval=FALSE}
data_files <- list.files(system.file("example-dataset/", package = "dataspice"),
  pattern = ".csv",
  full.names = TRUE
)
```

This function assumes that the metadata templates are in a folder called `metadata` within a `data` folder.

```{r, eval = FALSE}
attributes_path <- file.path("data", "metadata", "attributes.csv")
```

Using `purrr::map()`, this function can be applied over multiple files to populate the header names

```{r, eval=FALSE, message=FALSE}
data_files %>%
  purrr::map(~ prep_attributes(.x, attributes_path),
    attributes_path = attributes_path
  )
```

The output of `prep_attributes()` has the first two columns filled out:

```{r, echo=FALSE, message=FALSE}
readr::read_csv(
  system.file(
    "metadata-tables/attributes_prepped.csv",
    package = "dataspice"
  )
) %>%
  head() %>%
  kable()
```

#### Shiny helper apps

Each of the metadata templates can be edited interactively using a [Shiny](https://shiny.rstudio.com/) app:

* `edit_attributes()` opens a Shiny app that can be used to edit `attributes.csv`. The Shiny app displays the current `attributes` table and lets the user fill in an informative description and units (e.g. meters, hectares, etc.) for each variable.
* `edit_access()` opens an editable version of `access.csv`
* `edit_creators()` opens an editable version of `creators.csv`
* `edit_biblio()` opens an editable version of `biblio.csv`

![edit_attributes Shiny app](man/figures/edit_attributes.png)

Remember to click on **Save** when finished editing.

#### Completed metadata files

The first few rows of the completed metadata tables in this example will look like this:

`access.csv` has one row for each file

```{r, echo=FALSE, message=FALSE}
readr::read_csv(system.file("metadata-tables/access.csv", package = "dataspice")) %>%
  head() %>%
  kable(format = "markdown")
```

`attributes.csv` has one row for each variable in each file

```{r, echo=FALSE, message=FALSE}
readr::read_csv(system.file("metadata-tables/attributes.csv", package = "dataspice")) %>%
  head() %>%
  kable(format = "markdown")
```

`biblio.csv` is one row containing descriptors including spatial and temporal coverage

```{r, echo=FALSE, message=FALSE, warning=FALSE}
readr::read_csv(system.file("metadata-tables/biblio.csv", package = "dataspice")) %>%
  dplyr::mutate(description = str_trunc(description, 200, side = "right")) %>%
  kable(format = "markdown")
```

`creators.csv` has one row for each of the dataset authors

```{r, echo=FALSE, message=FALSE}
readr::read_csv(system.file("metadata-tables/creators.csv", package = "dataspice")) %>%
  kable(format = "markdown")
```

### Save JSON-LD file

`write_spice()` generates a json-ld file ("linked data") to aid in [dataset discovery](https://developers.google.com/search/docs/data-types/dataset), creation of more extensive metadata (e.g. [EML](https://eml.ecoinformatics.org)), and creating a website.

Here's a view of the `dataspice.json` file of the example data:

![listviewer pack output showing an example dataspice JSON file](man/figures/listviewer.png)

### Build website

* `build_site()` creates a bare-bones `index.html` file in the repository `docs` folder with a simple view of the dataset with the metadata and an interactive map. For example, this [repository](https://github.com/amoeba/dataspice-example/) results in this [website](https://amoeba.github.io/dataspice-example/)

![dataspice-website](man/figures/website_example.png)

### Convert to EML

The metadata fields `dataspice` uses are based largely on their compatibility with terms from [Schema.org](https://schema.org).
However, `dataspice` metadata can be converted to Ecological Metadata Language (EML), a much richer schema.
The conversion isn't perfect but `dataspice` will do its best to convert your `dataspice` metadata to EML:

```{r, as_eml}
library(dataspice)

# Load an example dataspice JSON that comes installed with the package
spice <- system.file(
  "examples", "annual-escapement.json",
  package = "dataspice"
)

# Convert it to EML
eml_doc <- spice_to_eml(spice)
```

You may receive warnings depending on which `dataspice` fields you filled in and this process will very likely produce an invalid EML record which is totally fine:

```{r, eml_validate}
library(EML)

eml_validate(eml_doc)
```

This is because some fields in `dataspice` store information in different structures and because EML requires many fields that `dataspice` doesn't have fields for.
At this point, you should look over the validation errors produced by `EML::eml_validate` and fix those.
Note that this will likely require familiarity with the [EML Schema](https://eml.ecoinformatics.org/) and the [EML package](https://github.com/ropensci/eml).

Once you're done, you can write out an EML XML file:

```{r, write_eml}
out_path <- tempfile()
write_eml(eml_doc, out_path)
```

### Convert from EML

Like converting `dataspice` to EML, we can convert an existing EML record to a set of `dataspice` metadata tables which we can then work from within `dataspice`:

```{r}
library(EML)

eml_path <- system.file("example-dataset/broodTable_metadata.xml", package = "dataspice")
eml <- read_eml(eml_path)
```

```{r, eval = FALSE}
# Creates four CSVs files in the `data/metadata` directory
my_spice <- eml_to_spice(eml, "data/metadata")
```

## Resources

A few existing tools & data standards to help users in specific domains:

* [Darwin Core](http://rs.tdwg.org/dwc/)
* [Ecological Metadata Language](https://knb.ecoinformatics.org/#external//emlparser/docs/index.html) (EML) (& [`EML`](https://github.com/ropensci/EML))
* [ISO 19115](https://www.iso.org/standard/53798.html) - Geographic Information Metadata
* [ISO 19139](https://www.iso.org/standard/32557.html) - Geographic Info Metadata XML schema
* [Minimum Information for Biological and Biomedical Investigations](https://fairsharing.org/collection/MIBBI) (MIBBI)

...And others indexed in [Fairsharing.org](https://fairsharing.org) & the [RDA metadata directory](http://rd-alliance.github.io/metadata-directory/standards/).

## Code of Conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/).
By contributing to this project, you agree to abide by its terms.

## Contributors

This package was developed at rOpenSci's 2018 unconf by (in alphabetical order):

* [Carl Boettiger](https://github.com/cboettig)
* [Scott Chamberlain](https://github.com/sckott)
* [Auriel Fournier](https://github.com/aurielfournier)
* [Kelly Hondula](https://github.com/khondula)
* [Anna Krystalli](https://github.com/annakrystalli)
* [Bryce Mecum](https://github.com/amoeba)
* [Maëlle Salmon](https://github.com/maelle)
* [Kate Webbink](https://github.com/magpiedin)
* [Kara Woo](https://github.com/karawoo)
* [Irene Steves](https://github.com/isteves)
---
title: "dataspice Overview"
author: "dataspice team"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{dataspice Overview}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

The goal of the `dataspice` package is to make it easier for researchers to create basic, lightweight, and concise metadata files for their datasets. These basic files can then be used to:

 - Make useful information available during analysis.
 - Create a helpful dataset README webpage.
 - Convert to more complex metadata formats to aid dataset discovery.

The `dataspice` metadata fields are based on [Schema.org](https://schema.org/) and and other, richer metadata standards such as [Ecological Metadata Language](https://eml.ecoinformatics.org/).

## Step 1 - Start with one or more data files

To start the user will have one or more datafiles in a common directory. We currently support rectangular data with headers in a .csv file or spatial data with attributes.

## Step 2 - Fill in Templates

`create_spice()` reads in the files from that directory and creates a set of template metadata files for the user to populate. These metadata templates will include some data extracted from the user's datafiles, including things like file names and measured variable names to aid the user in populating the metadata files.

Once these are created, you need to fill in each template file with the remaining metadata.
Templates can be filled either:

1. With the integrated [Shiny](https://shiny.rstudio.com/) applications (`edit_attributes()`, `edit_access()`, `edit_creators()`, `edit_biblio()`)
2. Manually (i.e., with a text or CSV/spreadsheet editor)

For new and even advanced users, option (1) is probably easier and more friendly.

### Metadata Files

- `creators.csv`: One row for each creator, and gives their affiliation, contact email, [ORCID](https://orcid.org/), etc.
- `attributes.csv`: This is where most of the user data entry will take place. For each variable, its name, units, and a written description are filled in.
- `biblio.csv`: Citation information about the project, as much or as little data as possible can be included, but if things like bounding box coordinates are not included, then when the website is generated there will not be a bounding box map generated.
- `access.csv`: Includes a row for each file that was read in, and documents the name of each file and its format.

<img style="width: 100%" src="https://user-images.githubusercontent.com/6106733/40335756-444642e0-5d1a-11e8-8c33-0d3053441a79.png" />

## Step 3 - Save metadata in JSON

Once the metadata templates are filled in, you can create a single JSON-LD metadata record for your metadata with `write_spice()`.

## Step 4 - Create a dataspice website

A dataset README website is an interactive representation of the JSON information about the data. Assuming sufficient information is provided in the JSON it will include a map of the points and a bounding box of the area of study.

The output from `write_spice()` is fed into `build_site()` to create the website.

<img style="width: 100%" src="../man/figures/website_example.png" />
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edit_file.R
\name{edit_creators}
\alias{edit_creators}
\title{Shiny apps for editing dataspice metadata tables}
\usage{
edit_creators(metadata_dir = file.path("data", "metadata"))
}
\arguments{
\item{metadata_dir}{the directory containing the \code{dataspice} metadata \code{.csv}
files. Defaults to \verb{data/metadata/} directory in \strong{current project root}.}
}
\description{
Launch Shiny app for editing individual \code{dataspice} metadata tables. Use
\verb{edit_*()} where \code{*} is one of the four \code{dataspice} metadata tables
\strong{\code{attributes}, \code{biblio}, \code{access} or \code{creators}}.
}
\examples{
\dontrun{
edit_attributes()
edit_biblio()
edit_access()
edit_creators()

# Specifying a different dataspice metadata directory
edit_attributes(metadata_dir = "analysis/data/metadata/"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_to_spice.R
\name{es_creators}
\alias{es_creators}
\title{Get creators from EML}
\usage{
es_creators(eml, path = NULL)
}
\arguments{
\item{eml}{(emld) an EML object}

\item{path}{(character) folder path for saving the table to disk}
}
\description{
Return EML creators in the dataspice creators.csv format.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_spice.R
\name{write_spice}
\alias{write_spice}
\title{Write spice}
\usage{
write_spice(path = file.path("data", "metadata"), ...)
}
\arguments{
\item{path}{location of metadata files}

\item{...}{additional arguments to \code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}}}
}
\value{
A JSON-LD file at the path specified
}
\description{
Write out your metadata as a dataspice JSON-LD document
}
\examples{
\dontrun{
# First create your metadata templates
create_spice()

# Then fill in the template files however you like

# Then write out your dataspice file
write_spice()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spice_to_eml.R
\name{spice_to_eml}
\alias{spice_to_eml}
\title{Convert \code{dataspice} metadata to EML}
\usage{
spice_to_eml(spice = file.path("data", "metadata", "dataspice.json"))
}
\arguments{
\item{spice}{(list) Your \code{dataspice} metadata. Uses
\code{data/metadata/dataspice.json} by default.}
}
\value{
(emld) The crosswalked \code{emld} object
}
\description{
Performs an (imperfect) conversion of \code{dataspice} metadata to EML. It's
very likely you will get validation errors and need to fix them afterwards
but \code{spice_to_eml} is a good way to a richer metadata schema (EML) when
you're already using \code{dataspice} but need a richer metadata schema.
}
\examples{
# Load an example dataspice JSON that comes installed with the package
spice <- system.file(
  "examples", "annual-escapement.json",
  package = "dataspice"
)

# And crosswalk it to EML
spice_to_eml(spice)

# We can also create dataspice metadata from scratch and crosswalk it to EML
myspice <- list(
  name = "My example spice",
  creator = "Me",
  contact = "Me"
)
spice_to_eml(myspice)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_crosswalk.R
\name{crosswalk_distribution}
\alias{crosswalk_distribution}
\title{Crosswalk a Schema.org/distribution}
\usage{
crosswalk_distribution(distribution)
}
\arguments{
\item{distribution}{(list) A distribution}
}
\description{
Crosswalk a Schema.org/distribution
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jsonld_to_mustache.R
\name{parse_GeoShape_box}
\alias{parse_GeoShape_box}
\title{Parse spatialCoverage$geo$box section for use in a Leaflet map}
\usage{
parse_GeoShape_box(box)
}
\arguments{
\item{box}{(list) spatialCoverage$geo$box section of the JSONLD}
}
\value{
(list) Template-specific variables for Leaflet
}
\description{
Parse spatialCoverage$geo$box section for use in a Leaflet map
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_crosswalk.R
\name{crosswalk_variables}
\alias{crosswalk_variables}
\title{Crosswalk \code{dataspice} variables to EML}
\usage{
crosswalk_variables(spice)
}
\arguments{
\item{spice}{(list) Your \code{dataspice} metadata}
}
\value{
(data.frame) A partial EML attributes table
}
\description{
See \code{\link[EML]{set_attributes}} for more information on what must be
filled out after this is run in order to get a valid EML \code{attributeList}.
}
\examples{
\dontrun{
# Load an example dataspice JSON that comes installed with the package
spice <- system.file(
  "examples", "annual-escapement.json",
  package = "dataspice")

# Convert it to EML (notice the warning)
eml_doc <- suppressWarnings({spice_to_eml(spice)})
attributes <- crosswalk_variables(spice)

# Now fill in the attributes data.frame. See `EML::set_attributes`.

# And last, set the attributes on our EML document
eml_doc$dataset$dataTable[[1]]$attributeList <-
  EML::set_attributes(attributes)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_metadata.R
\name{validate_creators}
\alias{validate_creators}
\title{Validate creators.csv}
\usage{
validate_creators(creators)
}
\arguments{
\item{creators}{(data.frame) A \code{data.frame} read in from
\code{creators.csv}}
}
\value{
Nothing. Side-effect: Can \code{stop} execution if validation fails.
}
\description{
Validate creators.csv
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_crosswalk.R
\name{crosswalk_datetime}
\alias{crosswalk_datetime}
\title{Convert a date(time) of unknown format into EML}
\usage{
crosswalk_datetime(input)
}
\arguments{
\item{input}{(character) Some unknown date(time) input}
}
\value{
(list) A \code{list} with members \code{calendarDate} and \code{time}. \code{time} will
be \code{NULL} if parsing fails or if the time string inside \code{input} isn't
ISO8601
}
\description{
A quick and dirty crosswalk of an unknown date(time) input to EML that really
only works for ISO8601 input. All other formats will fail and be returned
as-is as a \code{calendarDate}. From there the user will need to do a conversion
themselves.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_spice.R
\name{create_spice}
\alias{create_spice}
\title{Put metadata templates within a \code{metadata} subdirectory}
\usage{
create_spice(dir = "data")
}
\arguments{
\item{dir}{Directory containing data, within which a \code{metadata} subdirectory
will be created. Defaults to \code{data}.}
}
\description{
Put metadata templates within a \code{metadata} subdirectory
}
\examples{
\dontrun{
create_spice()

# Create templates from the data in a folder other than `data`
create_spice("my_data")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_attributes.R
\name{prep_attributes}
\alias{prep_attributes}
\title{Prepare attributes}
\usage{
prep_attributes(
  data_path = "data",
  attributes_path = file.path("data", "metadata", "attributes.csv"),
  ...
)
}
\arguments{
\item{data_path}{character vector of either:
\enumerate{
\item path(s) to the data file(s).
\item single path to directory containing data file(s).
Currently only tabular \code{.csv} and \code{.tsv} files are supported. Alternatively
attributes returned using \code{names()} can be extracted from r object, stored as
\code{.rds} files.
}}

\item{attributes_path}{path to the \verb{attributes.csv`` file. Defaults to }data/metadata/attributes.csv`.}

\item{...}{parameters passed to \code{list.files()}. For example, use
\code{recursive = TRUE} to list files in a folder recursively or use \code{pattern} to
filter files for patterns.}
}
\value{
\code{prep_attributes()} updates the \code{attributes.csv} and writes to
\code{attributes_path}.
}
\description{
Extract \code{variableNames} from data file(s) and add them to \code{attributes.csv}.
The helper \code{\link{validate_file_paths}} can be used to create vectors of
valid file paths that can be checked and then passed as \code{data_path} argument
to \code{\link{prep_attributes}}.
}
\examples{
\dontrun{
create_spice()
# extract attributes from all `csv`, `tsv`, `rds` files in the data folder
# (non recursive)
prep_attributes()
# recursive
prep_attributes(recursive = TRUE)
# extract attributes from a single file using file path
data_path <- system.file("example-dataset","BroodTables.csv",
                         package = "dataspice")
prep_attributes(data_path)
# extract attributes from a single file by file path pattern matching
data_path <- system.file("example-dataset", package = "dataspice")
prep_attributes(data_path, pattern = "StockInfo")
# extract from a folder using folder path
data_path <- system.file("example-dataset", package = "dataspice")
prep_attributes(data_path)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_site.R
\name{build_site}
\alias{build_site}
\title{Build a dataspice site}
\usage{
build_site(
  path = file.path("data", "metadata", "dataspice.json"),
  template_path = system.file("template.html5", package = "dataspice"),
  out_path = file.path("docs", "index.html")
)
}
\arguments{
\item{path}{(character) Path to a JSON+LD file with dataspice metadata}

\item{template_path}{(character) Optional. Path to a template for
\code{\link[whisker]{whisker.render}}}

\item{out_path}{(character) Optional. Path to write the site's \code{index.html}
to. Defaults to \code{docs/index.html}.}
}
\value{
Nothing. Creates/overwrites \code{docs/index.html}
}
\description{
Build a dataspice site
}
\examples{
\dontrun{
# Create JSON+LD from a set of metadata templates
json <- write_json(biblio, access, attributes, creators)
build_site(json)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_jsonld.R
\name{as_jsonld}
\alias{as_jsonld}
\title{Convert a list object to JSON-LD}
\usage{
as_jsonld(
  x,
  context = "http://schema.org",
  pretty = TRUE,
  auto_unbox = TRUE,
  ...
)
}
\arguments{
\item{x}{the object to be encoded.}

\item{context}{JSON-LD context; "http://schema.org".}

\item{pretty}{Whether or not to prettify output. See
\code{\link[jsonlite]{toJSON}}.}

\item{auto_unbox}{Whether or not to automatically unbox output. See
\code{\link[jsonlite]{toJSON}}.}

\item{...}{Other arguments to be passed to
\code{\link[jsonlite]{toJSON}}.}
}
\description{
Convert a list object to JSON-LD
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_to_spice.R
\name{eml_to_spice}
\alias{eml_to_spice}
\title{Create dataspice metadata tables from EML}
\usage{
eml_to_spice(eml, path = NULL)
}
\arguments{
\item{eml}{(emld) An EML object}

\item{path}{(character) Folder path for saving the tables to disk}
}
\value{
A list with names \code{attributes}, \code{access}, \code{biblio}, and \code{creators}.
Optionally, if \code{path} is specified, saves the four tables as \code{CSV} files.
}
\description{
Create dataspice metadata tables from EML
}
\examples{
\dontrun{
# First, load up an example EML record
library(EML)

eml_path <- system.file(
 file.path("example-dataset", "broodTable_metadata.xml"),
 package = "dataspice")
eml <- read_eml(eml_path)


# Generate the four dataspice tables
my_spice <- eml_to_spice(eml)

# Or save them as a file
# Generate the four dataspice tables
eml_to_spice(eml, ".")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_jsonld.R
\name{write_jsonld}
\alias{write_jsonld}
\title{Write a list out as object to JSON-LD}
\usage{
write_jsonld(
  x,
  path,
  context = "http://schema.org",
  pretty = TRUE,
  auto_unbox = TRUE,
  ...
)
}
\arguments{
\item{x}{an object to be serialized to JSON}

\item{path}{file on disk}

\item{context}{JSON-LD context; "http://schema.org"}

\item{pretty}{adds indentation whitespace to JSON output. Can be TRUE/FALSE or a number specifying the number of spaces to indent. See \code{\link[jsonlite]{prettify}}}

\item{auto_unbox}{automatically \code{\link[jsonlite]{unbox}} all atomic vectors of length 1. It is usually safer to avoid this and instead use the \code{\link[jsonlite]{unbox}} function to unbox individual elements.
An exception is that objects of class \code{AsIs} (i.e. wrapped in \code{I()}) are not automatically unboxed. This is a way to mark single values as length-1 arrays.}

\item{...}{additional conversion arguments, see also \link[jsonlite]{toJSON} or \link[jsonlite]{fromJSON}}
}
\description{
Write a list out as object to JSON-LD
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_metadata.R
\name{validate_access}
\alias{validate_access}
\title{Validate access.csv}
\usage{
validate_access(access)
}
\arguments{
\item{access}{(data.frame) A \code{data.frame} read in from \code{access.csv}}
}
\value{
Nothing. Side-effect: Can \code{stop} execution if validation fails.
}
\description{
Validate access.csv
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_to_spice.R
\name{es_attributes}
\alias{es_attributes}
\title{Get attributes from EML}
\usage{
es_attributes(eml, path = NULL)
}
\arguments{
\item{eml}{(emld) an EML object}

\item{path}{(character) folder path for saving the table to disk}
}
\description{
Return EML attributes in the dataspice attributes.csv format.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_crosswalk.R
\name{crosswalk_creator}
\alias{crosswalk_creator}
\title{Crosswalk a Schema.org/creator}
\usage{
crosswalk_creator(creator)
}
\arguments{
\item{creator}{(list) A creator}
}
\description{
Crosswalk a Schema.org/creator
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jsonld_to_mustache.R
\name{jsonld_to_mustache}
\alias{jsonld_to_mustache}
\title{Convert JSONLD to a list suitable for Mustache templating}
\usage{
jsonld_to_mustache(path)
}
\arguments{
\item{path}{(character) Path to file on disk to convert}
}
\value{
(list) Mustache-appropriate list
}
\description{
Convert JSONLD to a list suitable for Mustache templating
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataspice-package.R
\docType{package}
\name{dataspice-package}
\alias{dataspice}
\alias{dataspice-package}
\title{dataspice: Create Lightweight Schema.org Descriptions of Data}
\description{
The goal of 'dataspice' is to make it easier for researchers to
  create basic, lightweight, and concise metadata files for their datasets.
  These basic files can then be used to make useful information available during
  analysis, create a helpful dataset "README" webpage, and produce more complex
  metadata formats to aid dataset discovery. Metadata fields are based on
  the 'Schema.org' and 'Ecological Metadata Language' standards.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/dataspice}
  \item Report bugs at \url{https://github.com/ropensci/dataspice/issues}
}

}
\author{
\strong{Maintainer}: Bryce Mecum \email{brycemecum@gmail.com} (https://github.com/amoeba)

Authors:
\itemize{
  \item Carl Boettiger (https://github.com/cboettig)
  \item Scott Chamberlain (https://github.com/sckott)
  \item Auriel Fournier (https://github.com/aurielfournier)
  \item Kelly Hondula (https://github.com/khondula)
  \item Anna Krystalli (https://github.com/annakrystalli)
  \item Maëlle Salmon (https://github.com/maelle)
  \item Kate Webbink (https://github.com/magpiedin)
  \item Kara Woo (https://github.com/karawoo)
}

Other contributors:
\itemize{
  \item Irene Steves (https://github.com/isteves) [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serve_site.R
\name{serve_site}
\alias{serve_site}
\title{Serve site}
\usage{
serve_site(path = "docs")
}
\arguments{
\item{path}{(character) Optional. Directory to serve. Defaults to
\code{docs}.}
}
\value{
Nothing.
}
\description{
Serve site
}
\examples{
\dontrun{
# Build your site
json <- write_json(biblio, access, attributes, creators)
build_site(json)

# Serve it
serve_site()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_crosswalk.R
\name{crosswalk_Person}
\alias{crosswalk_Person}
\title{Crosswalk functions for \code{as_eml}
Crosswalk a Schema.org/Person}
\usage{
crosswalk_Person(creator)
}
\arguments{
\item{creator}{(list) A creator}
}
\description{
Crosswalk functions for \code{as_eml}
Crosswalk a Schema.org/Person
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_attributes.R
\name{validate_file_paths}
\alias{validate_file_paths}
\title{Validate file paths}
\usage{
validate_file_paths(data_path = "data", ...)
}
\arguments{
\item{data_path}{character vector of either:
\enumerate{
\item path(s) to the data file(s).
\item single path to directory containing data file(s).
Currently only tabular \code{.csv} and \code{.tsv} files are supported. Alternatively
attributes returned using \code{names()} can be extracted from r object, stored as
\code{.rds} files.
}}

\item{...}{parameters passed to \code{list.files()}. For example, use
\code{recursive = TRUE} to list files in a folder recursively or use \code{pattern} to
filter files for patterns.}
}
\value{
One or more data file paths
}
\description{
Helper function to return a set of file paths for use in other functions
}
\examples{
\dontrun{
# Assuming some data files in "./data"
my_files <- validate_file_paths()

# If your data files are in `another_folder`
my_files <- validate_file_paths("another_folder")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edit_file.R
\name{edit_attributes}
\alias{edit_attributes}
\title{Shiny apps for editing dataspice metadata tables}
\usage{
edit_attributes(metadata_dir = file.path("data", "metadata"))
}
\arguments{
\item{metadata_dir}{the directory containing the \code{dataspice} metadata \code{.csv}
files. Defaults to \verb{data/metadata/} directory in \strong{current project root}.}
}
\description{
Launch Shiny app for editing individual \code{dataspice} metadata tables. Use
\verb{edit_*()} where \code{*} is one of the four \code{dataspice} metadata tables
\strong{\code{attributes}, \code{biblio}, \code{access} or \code{creators}}.
}
\examples{
\dontrun{
edit_attributes()
edit_biblio()
edit_access()
edit_creators()

# Specifying a different dataspice metadata directory
edit_attributes(metadata_dir = "analysis/data/metadata/"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_crosswalk.R
\name{crosswalk}
\alias{crosswalk}
\title{Crosswalk a term}
\usage{
crosswalk(doc, term)
}
\arguments{
\item{doc}{(list) A \code{dataspice} document as a \code{list}}

\item{term}{(character) The term to crosswalk.}
}
\value{
(list) The result of the crosswalk. May be an empty \code{list} on
failure.
}
\description{
Crosswalk a term
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_to_spice.R
\name{es_access}
\alias{es_access}
\title{Get access from EML}
\usage{
es_access(eml, path = NULL)
}
\arguments{
\item{eml}{(emld) an EML object}

\item{path}{(character) folder path for saving the table to disk}
}
\description{
Return EML access in the dataspice access.csv format.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edit_file.R
\name{edit_biblio}
\alias{edit_biblio}
\title{Shiny apps for editing dataspice metadata tables}
\usage{
edit_biblio(metadata_dir = file.path("data", "metadata"))
}
\arguments{
\item{metadata_dir}{the directory containing the \code{dataspice} metadata \code{.csv}
files. Defaults to \verb{data/metadata/} directory in \strong{current project root}.}
}
\description{
Launch Shiny app for editing individual \code{dataspice} metadata tables. Use
\verb{edit_*()} where \code{*} is one of the four \code{dataspice} metadata tables
\strong{\code{attributes}, \code{biblio}, \code{access} or \code{creators}}.
}
\examples{
\dontrun{
edit_attributes()
edit_biblio()
edit_access()
edit_creators()

# Specifying a different dataspice metadata directory
edit_attributes(metadata_dir = "analysis/data/metadata/"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edit_file.R
\name{edit_access}
\alias{edit_access}
\title{Shiny apps for editing dataspice metadata tables}
\usage{
edit_access(metadata_dir = file.path("data", "metadata"))
}
\arguments{
\item{metadata_dir}{the directory containing the \code{dataspice} metadata \code{.csv}
files. Defaults to \verb{data/metadata/} directory in \strong{current project root}.}
}
\description{
Launch Shiny app for editing individual \code{dataspice} metadata tables. Use
\verb{edit_*()} where \code{*} is one of the four \code{dataspice} metadata tables
\strong{\code{attributes}, \code{biblio}, \code{access} or \code{creators}}.
}
\examples{
\dontrun{
edit_attributes()
edit_biblio()
edit_access()
edit_creators()

# Specifying a different dataspice metadata directory
edit_attributes(metadata_dir = "analysis/data/metadata/"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_to_spice.R
\name{es_biblio}
\alias{es_biblio}
\title{Get biblio from EML}
\usage{
es_biblio(eml, path = NULL)
}
\arguments{
\item{eml}{(emld) an EML object}

\item{path}{(character) folder path for saving the table to disk}
}
\description{
Return EML biblio in the dataspice biblio.csv format.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jsonld_to_mustache.R
\name{parse_spatialCoverage}
\alias{parse_spatialCoverage}
\title{Parse spatialCoverage section for use in a Leaflet map}
\usage{
parse_spatialCoverage(spatialCoverage)
}
\arguments{
\item{spatialCoverage}{(list) spatialCoverage section of the JSONLD}
}
\value{
(list) Template-specific variables for Leaflet
}
\description{
Parse spatialCoverage section for use in a Leaflet map
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_metadata.R
\name{validate_biblio}
\alias{validate_biblio}
\title{Validate biblio.csv}
\usage{
validate_biblio(biblio)
}
\arguments{
\item{biblio}{(data.frame) A \code{data.frame} read in from \code{biblio.csv}}
}
\value{
Nothing. Side-effect: Can \code{stop} execution if validation fails.
}
\description{
Validate biblio.csv
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_access.R
\name{prep_access}
\alias{prep_access}
\title{Prepare access}
\usage{
prep_access(
  data_path = "data",
  access_path = file.path("data", "metadata", "access.csv"),
  ...
)
}
\arguments{
\item{data_path}{character vector of either:
\enumerate{
\item path(s) to the data file(s).
\item single path to directory containing data file(s).
Currently only tabular \code{.csv} and \code{.tsv} or \code{.rds} files are supported.
}}

\item{access_path}{path to the \code{access.csv} file. Defaults to
\code{data/metadata/access.csv}.}

\item{...}{parameters passed to \code{list.files()}. For example, use
\code{recursive = TRUE} to list files in a folder recursively or use \code{pattern} to
filter files for patterns.}
}
\value{
Updates \code{access.csv} and writes to \code{access_path}.
}
\description{
Extract \code{fileNames} from data file(s) and add them to \code{access.csv}. The
helper \code{\link{validate_file_paths}} can be used to create vectors of
valid file paths that can be checked and then passed as \code{data_path} argument
to \code{\link{prep_access}}.
}
\examples{
\dontrun{
# First create the metadata tempaltes
create_spice()

# Then begin filling them in from your data files
prep_access()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_crosswalk.R
\name{crosswalk_Organization}
\alias{crosswalk_Organization}
\title{Crosswalk a Schema.org/Organization}
\usage{
crosswalk_Organization(creator)
}
\arguments{
\item{creator}{(list) A creator}
}
\description{
Crosswalk a Schema.org/Organization
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_metadata.R
\name{validate_attributes}
\alias{validate_attributes}
\title{Validate attributes.csv}
\usage{
validate_attributes(attributes)
}
\arguments{
\item{attributes}{(data.frame) A \code{data.frame} read in from
\code{attributes.csv}}
}
\value{
Nothing. Side-effect: Can \code{stop} execution if validation fails.
}
\description{
Validate attributes.csv
}
