
<!-- README.md is generated from README.Rmd. Please edit that file -->

# citesdb

Authors: *Noam Ross, Evan A. Eskew and Mauricio Vargas*

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![rOpensci\_Badge](https://badges.ropensci.org/292_status.svg)](https://github.com/ropensci/software-review/issues/292)
[![Published in the Journal of Open Source
Software](http://joss.theoj.org/papers/10.21105/joss.01483/status.svg)](https://doi.org/10.21105/joss.01483)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2630836.svg)](https://doi.org/10.5281/zenodo.2630836)
[![CircleCI](https://circleci.com/gh/ropensci/citesdb/tree/master.svg?style=shield)](https://circleci.com/gh/ropensci/citesdb)
[![codecov](https://codecov.io/gh/ropensci/citesdb/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/citesdb)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

<!-- badges: end -->

**citesdb** is an R package to conveniently analyze the full CITES
shipment-level wildlife trade database, available at
<https://trade.cites.org/>. This data consists of over 40 years and 20
million records of reported shipments of wildlife and wildlife products
subject to oversight under the [Convention on International Trade in
Endangered Species of Wild Fauna and Flora](https://www.cites.org). The
source data are maintained by the [UN Environment World Conservation
Monitoring Centre](https://www.unep-wcmc.org/).

## Installation

Install the **citesdb** package with this command:

``` r
devtools::install_github("ropensci/citesdb")
```

Note that since **citesdb** installs a source dependency from GitHub,
you will need [package build
tools](http://stat545.com/packages01_system-prep.html).

## Usage

### Getting the data

When you first load the package, you will see a message like this:

    library(citesdb)
    #> Local CITES database empty or corrupt. Download with cites_db_download()

Not to worry, just do as it says and run `cites_db_download()`. This
will fetch the most recent database from online, an approximately 158 MB
download. It will expand to over 1 GB in the local database. During the
download and database building, up to 3.5 GB of disk space may be used
temporarily.

### Using the database

Once you fetch the data, you can connect to the database with the
`cites_db()` command. The `cites_shipments()` command loads a remote
`tibble` that is backed by the database but is not loaded into R. You
can use this command to analyze CITES data without ever loading it into
memory, gathering your results with the `dplyr` function `collect()`.
For example:

``` r
library(citesdb)
library(dplyr)

start <- Sys.time()

cites_shipments() %>%
  group_by(Year) %>%
  summarize(n_records = n()) %>%
  arrange(desc(Year)) %>%
  collect()
#> # A tibble: 45 x 2
#>     Year n_records
#>    <int>     <dbl>
#>  1  2019     12610
#>  2  2018   1143044
#>  3  2017   1246684
#>  4  2016   1293178
#>  5  2015   1299183
#>  6  2014   1109877
#>  7  2013   1127377
#>  8  2012   1096664
#>  9  2011    950148
#> 10  2010    894115
#> # … with 35 more rows

stop <- Sys.time()
```

(*Note that running `collect()` on all of `cites_shipments()` will load
a \>3 GB data frame into memory\!*)

The back-end database, [duckdb](https://duckdb.org/), is very fast and
powerful, making analyses on such large data quite snappy using normal
desktops and laptops. Here’s the timing of the above query, which
processes over 20 million records:

``` r
stop - start
#> Time difference of 0.4658868 secs
```

If you are using a recent version of RStudio interactively, loading the
CITES package also brings up a browsable pane in the “Connections” tab
that lets you explore and preview the database, as well as interact with
it directly via SQL commands.

If you don’t need any of the bells and whistles of this package, you can
download the raw data as a single compressed TSV file from the [releases
page](https://github.com/ropensci/citesdb/releases), or as a `.zip` file
of many CSV files from the original source at
<https://trade.cites.org/>.

### Metadata

The package database also contains tables of field metadata, codes used,
and CITES countries. This information comes from [“A guide to using the
CITES Trade
Database”](https://trade.cites.org/cites_trade_guidelines/en-CITES_Trade_Database_Guide.pdf),
on the CITES website. Convenience functions `cites_metadata()`,
`cites_codes()`, and `cites_parties()` access this information:

``` r
head(cites_metadata())
#> # A tibble: 6 x 2
#>   variable description                                 
#>   <chr>    <chr>                                       
#> 1 Year     year in which trade occurred                
#> 2 Appendix CITES Appendix of taxon concerned           
#> 3 Taxon    scientific name of animal or plant concerned
#> 4 Class    scientific name of animal or plant concerned
#> 5 Order    scientific name of animal or plant concerned
#> 6 Family   scientific name of animal or plant concerned

head(cites_codes())
#> # A tibble: 6 x 3
#>   field   code  description                                    
#>   <chr>   <chr> <chr>                                          
#> 1 Purpose B     Breeding in captivity or artificial propagation
#> 2 Purpose E     Educational                                    
#> 3 Purpose G     Botanical garden                               
#> 4 Purpose H     Hunting trophy                                 
#> 5 Purpose L     Law enforcement / judicial / forensic          
#> 6 Purpose M     Medical (including biomedical research)

head(cites_parties())
#> # A tibble: 6 x 6
#>   country        code  former_code non_ISO_code date       data_source                                                  
#>   <chr>          <chr> <lgl>       <lgl>        <chr>      <chr>                                                        
#> 1 Afghanistan    AF    FALSE       FALSE        1986-01-28 'A guide to using the CITES Trade Database', Version 8, Anne…
#> 2 Africa         XF    FALSE       TRUE         <NA>       'A guide to using the CITES Trade Database', Version 8, Anne…
#> 3 Åland Islands  AX    FALSE       FALSE        <NA>       'A guide to using the CITES Trade Database', Version 8, Anne…
#> 4 Albania        AL    FALSE       FALSE        2003-09-25 'A guide to using the CITES Trade Database', Version 8, Anne…
#> 5 Algeria        DZ    FALSE       FALSE        1984-02-21 'A guide to using the CITES Trade Database', Version 8, Anne…
#> 6 American Samoa AS    FALSE       FALSE        <NA>       'A guide to using the CITES Trade Database', Version 8, Anne…
```

More information on the release of shipment-level CITES data can be
found in the `?guidance` help file.

## Related work

The [**rcites**](https://github.com/ropensci/rcites) package provides
access to the Speciesplus/CITES Checklist API, which includes metadata
about species and their protected status through time.

## Citation

If you use **citesdb** in a publication, please cite both the package
and source data:

Ross, Noam, Evan A. Eskew, and Nicolas Ray. 2019. citesdb: An R package
to support analysis of CITES Trade Database shipment-level data. Journal
of Open Source Software, 4(37), 1483,
<https://doi.org/10.21105/joss.01483>

UNEP-WCMC (Comps.) 2019. Full CITES Trade Database Download. Version
2019.2. CITES Secretariat, Geneva, Switzerland. Compiled by UNEP-WCMC,
Cambridge, UK. Available at: <https://trade.cites.org>.

## Contributing

Have feedback or want to contribute? Great\! Please take a look at the
[contributing
guidelines](https://github.com/ropensci/citesdb/blob/master/.github/CONTRIBUTING.md)
before filing an issue or pull request.

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/citesdb/blob/master/.github/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.

[![Created by EcoHealth
Alliance](https://raw.githubusercontent.com/ropensci/citesdb/master/vignettes/figures/eha-footer.png)](https://www.ecohealthalliance.org/)
# citesdb 0.1.0

* First release with public data
# Contributor Covenant Code of Conduct

[Bosnian](http://contributor-covenant.org/version/1/4/bs/) 
| [Deutsch](http://contributor-covenant.org/version/1/4/de/) 
| [ελληνικά](http://contributor-covenant.org/version/1/4/el/) 
| [English](http://contributor-covenant.org/version/1/4/) 
| [Español](http://contributor-covenant.org/version/1/4/es/) 
| [Français](http://contributor-covenant.org/version/1/4/fr/) 
| [Italiano](http://contributor-covenant.org/version/1/3/0/it/) 
| [日本語](http://contributor-covenant.org/version/1/3/0/ja/) 
| [Magyar](http://contributor-covenant.org/version/1/3/0/hu/) 
| [Nederlands](http://contributor-covenant.org/version/1/4/nl/) 
| [Polski](http://contributor-covenant.org/version/1/4/pl/) 
| [Português](http://contributor-covenant.org/version/1/4/pt/) 
| [Português do Brasil](http://contributor-covenant.org/version/1/4/pt_br/) 
| [Pусский](http://contributor-covenant.org/version/1/3/0/ru/) 
| [Română](http://contributor-covenant.org/version/1/4/ro/) 
| [Svenska](http://contributor-covenant.org/version/1/4/sv/) 
| [Slovenščina](http://contributor-covenant.org/version/1/4/sl/) 
| [Türkçe](http://contributor-covenant.org/version/1/4/tr/) 
| [Українська](http://contributor-covenant.org/version/1/4/uk/) 
| [한국어](http://contributor-covenant.org/version/1/4/ko/)


## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at noam.ross@gmail.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Contributing to the citesdb R package

You want to contribute to **citesdb**? Great! 

Please submit questions, bug reports, and requests in the [issues tracker](https://github.com/ropensci/citesdb/issues). Please submit bug
reports with a minimal  [reprex](https://www.tidyverse.org/help/#reprex).

If you plan to contribute code, go ahead and fork the repo and submit a pull request. A few notes:

-   This package is released with a [Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.  Why? We want contribution to be enjoyable and rewarding for everyone!
-   If you have large change, please open an issue first to discuss.
-   I'll generally include contributors as authors in the DESCRIPTION file (with
their permission) for most contributions that go beyond small typos in code or documentation.
-   This package generally uses the [rOpenSci packaging guidelines](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md) for style and structure.
-   Documentation is generated by **roxygen2**. Please write documentation in code files and let it auto-generate documentation files.  We use a recent version so documentation my be [written in markdown](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html) 
-   We aim for testing that has high coverage and is robust.  Include tests with
   any major contribution to code. Test your changes the package with [**goodpractice**](https://cran.r-project.org/web/packages/goodpractice/index.html) before
submitting your change, and run spelling::spell_check_package() and lintr::lint_package(), both of which are tested.


---
title: 'citesdb: An R package to support analysis of CITES Trade Database shipment-level data'
tags:
  - R
  - database
  - wildlife
  - trade
  - conservation
  - sustainability
authors:
 - name: Noam Ross
   orcid: 0000-0002-2136-0000
   affiliation: 1
 - name: Evan A. Eskew
   orcid: 0000-0002-1153-5356
   affiliation: 1
 - name: Nicolas Ray
   affiliation: "2, 3"
affiliations:
 - name: EcoHealth Alliance, New York, New York, USA
   index: 1
 - name: GeoHealth Group, Institute of Global Health, Faculty of Medicine, University of Geneva, Geneva, Switzerland
   index: 2
 - name: Institute for Environmental Sciences, University of Geneva, Geneva, Switzerland
   index: 3
date: 21 May 2019
bibliography: paper.bib
---

# Summary

International trade is a significant threat to wildlife globally [@Bennett_2002; @Lenzen_2012; @Bush_2014; @Tingley_2017]. Consequently, high-quality, widely accessible data on the wildlife trade are urgently needed to generate effective conservation strategies and action [@Joppa_2016]. The [Convention on International Trade in Endangered Species of Wild Fauna and Flora](https://www.cites.org) (CITES) provides a key wildlife trade dataset for conservationists, the CITES Trade Database, which is maintained by the [UN Environment World Conservation Monitoring Centre](https://www.unep-wcmc.org/). Broadly, CITES is a trade oversight mechanism which aims to limit the negative effects of overharvesting, and the CITES Trade Database represents compiled data from CITES Parties regarding the trade of wildlife or wildlife products listed under the CITES Appendices. Despite data complexities that can complicate interpretation [@Harrington_2015; @Lopes_2017; @Berec_2018; @Robinson_2018; @Eskew_2019], the CITES Trade Database remains a critically important resource for evaluating the extent and impact of the legal, international wildlife trade [@Harfoot_2018].

`citesdb` is an R package designed to support analysis of the recently released shipment-level CITES Trade Database [@tradedb]. Currently, the database contains over 40 years and 20 million records of shipments of wildlife and wildlife products subject to reporting under CITES, including individual shipment permit IDs that have been anonymized by hashing, and accompanying metadata. @Harfoot_2018 provide a recent overview of broad temporal and spatial trends in this data. To facilitate further analysis of this large dataset, the `citesdb` package imports the CITES Trade Database into a local, on-disk embedded database [@duckdb]. This avoids the need for users to pre-process the data or load the multi-gigabyte dataset into memory. The DuckDB back-end allows high-performance querying and is accessible via a `DBI`- and `dplyr`-compatible interface familiar to most R users [@DBI; @dplyr]. For users of the RStudio integrated development environment [@rstudio], the package also provides an interactive pane for exploring the database and previewing data. `citesdb` has undergone [code review at rOpenSci](https://github.com/ropensci/software-review/issues/292).

# Acknowledgements

Authors N. Ross and E. A. Eskew were funded by the generous support of the American people through the United States Agency for International Development (USAID) Emerging Pandemic Threats PREDICT project.

# References
---
slug: "citesdb"
title: "
package_version: 0.1.0
authors:
  - Noam Ross
  - Evan A. Eskew
date: 2019-06-01
categories: blog
topicid:
tags:
- Software Peer Review
- R
- community
- software
- packages
- citesdb
- wildlife
- databases
---

# topics

At EHA, we work with wildlife trade data to help understand risks to wildlife
and their potential as vectors for infectious disease.  Critical in controlling,
monitoring, and understanding this trade is the Convention on Trade in
Endangered Species, a treaty and organization that regulates and monitors the
international wildlife trade.

Recently, CITES, together with WCMC, began publishing their data at a new level
of detail - individual shipment-level records of trade going back 40 years.


Importance of context and metadata for analyzing data.

I've been an editor for rOpenSci for five years, but this is the first time I
put a package of my own through the process!

citesdb was also an an opportunity to put into practice some  [ideas
for data package design](https://www.youtube.com/watch?v=zsEsh5QpN0U).
The CITES trade database falls into a common category of "biggish" data - larger
than RAM for many users, but still small enough to analyze on a single machine.

Thanks to reviewers and editors.
# Response to reviewers

Thanks to both reviewers for their comments and issues that have 
helped us to make a more user-friendly package.

## Changes implemented

- We have changed the installation method in the `README` to use **devtools**.
  
- We have added information to the `cites_status` table indicating the location 
  of the database files on disk.

- We now clear out and disconnect the database tables after updating, which
  should trigger disk cleanup and avoid doubling database size.
  
- Thanks for the reviewers for sticking with us through a long and merry chase
  of a bug that was preventing package-building in 
  https://github.com/ecohealthalliance/citesdb/issues/1, which had its origins
  in a missing token for use with the **rcites** package as well as a low-level
  database lock issue. We now cache this the **rcites** information to avoid
  having to make remote calls in vignette-building, and also have resolved the
  DB locking conflict.

- We have modified our linter tests to avoid the false positives shown in
  https://github.com/ecohealthalliance/citesdb/issues/2.
  
- We have renamed the internal function `check_status()` to `cites_check_status()`
  to be less generic.
  
- We have made more elaborate help and examples for the `cites_db()`, 
  `cites_shipments()`, and metadata functions to illustrate their use and
  distinguish between **dplyr**- and **DBI**-based workflows.

## Changes not implemented / justifications

- We have opted not to use the **glue** package and instead stick with the base
  R `paste()` functions in the interest of limiting dependencies. We believe
  this is a minor trade-off.
  
- The low test coverage shown in https://github.com/ecohealthalliance/citesdb/issues/4
  is due to tests skipped on CRAN. Setting the environment variable to 
  `NOT_CRAN=true` shows that our test coverage is 65%. This is lower than
  typical, but as we note in the issue above, this is largely due to the 
  verbose code for interacting with the RStudio connection pane, which cannot be
  tested except in an interactive session. Other code coverage is 
  [greater than 90%](https://codecov.io/gh/ecohealthalliance/citesdb/tree/master/R).
  We believe that all core functionality is tested, including important
  edge-cases and conditions not reflected in the coverage statistic, such as
  error handling in multiple sessions and changing up upstream data sources.
  
A final note: we have learned from the maintainers of **MonetDBLite** that it will 
not be returning to CRAN (its current iteration fails on R-devel), but they are
working on a successor embedded database package that will replace it and go
to CRAN later this year. So, for now, we will host this package on GitHub and
replace the database back-end and send to CRAN when the successor package is ready.
We've added a note to the README showing that users need package build tools for
the current version.
