# rcites <img src="man/figures/rcites_logo.png" width="130" height="150" align="right"/>

[![R CMD Check](https://github.com/ropensci/rcites/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ropensci/rcites/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/ropensci/rcites/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/rcites)
[![status](https://tinyverse.netlify.com/badge/rcites)](https://CRAN.R-project.org/package=rcites)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![ROpenSci status](https://badges.ropensci.org/244_status.svg)](https://github.com/ropensci/software-review/issues/244)
[![CRAN status](https://www.r-pkg.org/badges/version/rcites)](https://CRAN.R-project.org/package=rcites)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/rcites)](https://cran.r-project.org/package=rcites)
[![Zenodo DOI](https://zenodo.org/badge/113842199.svg)](https://zenodo.org/badge/latestdoi/113842199)
[![JOSS DOI](http://joss.theoj.org/papers/10.21105/joss.01091/status.svg)](https://doi.org/10.21105/joss.01091)


An R package to access information from the [Speciesplus](https://speciesplus.net/) database via the [Speciesplus/CITES Checklist API](https://api.speciesplus.net/documentation/v1.html). The package is available for download from [CRAN](https://cran.r-project.org/package=rcites) (stable version) and [Github](https://github.com/ropensci/rcites) (development version).

Please see the [release paper](https://doi.org/10.21105/joss.01091) for background information about the Convention on International Trade in Endangered Species of Wild Fauna and Flora ([CITES](https://cites.org/eng)), the Speciesplus database and basic information about the aim of the package.


### Installation

The package can be installed from CRAN:

```R
install.packages("rcites")
```

The development version can be installed via the [`remotes`](https://CRAN.R-project.org/package=remotes) :package:

```R
remotes::install_github("ropensci/rcites")
```


### Setup requirements and use

To set up a connection to the CITES Speciesplus database, a personal authentication token is required. Please see the vignette for details how to get a token and how to set the token for package use: [Get started with rcites](https://docs.ropensci.org/rcites/articles/a_get_started.html)

Additional information about specific use examples are provided for the
[African bush elephant (*Loxodonta africana*)](https://docs.ropensci.org/rcites/articles/b_elephant.html).
The package usage for querying multiple species is described in another
vignette entitled ['Bulk analysis with rcites'](https://docs.ropensci.org/rcites/articles/c_bulk_analysis.html).


### Key features

Once the token is set, the package has five key features:

- `spp_taxonconcept()`: [access the Speciesplus taxon concept](https://api.speciesplus.net/documentation/v1/taxon_concepts/index.html) and retrieve a taxon id
- `spp_cites_legislation()`: [access CITES legislation data](https://api.speciesplus.net/documentation/v1/cites_legislation/index.html)
- `spp_eu_legislation()`: [access EU legislation data](https://api.speciesplus.net/documentation/v1/eu_legislation/index.html)
- `spp_distributions()`: [access a taxon distribution data](https://api.speciesplus.net/documentation/v1/distributions/index.html)
- `spp_references()`: [access a listing reference data](https://api.speciesplus.net/documentation/v1/references/index.html)


### Prefix information

The package functions have three different prefixes:

- `set_` for `set_token()` to initially set the API token
- `spp_` for the key features
- `rcites_` for helper functions that are called within the key features


### Citation information

When citing, please refer to both the [package citation](https://docs.ropensci.org/rcites/authors.html) and the [release paper](https://doi.org/10.21105/joss.01091).


## Contributors

- [Main contributors](https://github.com/ropensci/rcites/graphs/contributors)

- Reviewers of the package:
  - [Noam Ross](https://github.com/noamross)
  - [Margaret Siple](https://github.com/mcsiple)

- Editor: [Scott Chamberlain](https://github.com/sckott)

- Reporting issue(s):
  - @FVFaleiro
  - @eveskew
  - @fleurhierink
  - @wajra



## Resources

Another package dealing with data from and about CITES, providing access to its
wildlife trade database: [cites](https://github.com/ecohealthalliance/cites/)

While creating this package, we greatly benefited from:

1. [taxize](https://github.com/ropensci/taxize) that inspired the structuring of this repository/package;

2. the `httr` vignette: [Managing secrets](https://CRAN.R-project.org/package=httr/vignettes/secrets.html), which is extremely helpful for packages dealing with API.



## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://docs.ropensci.org/rcites/CONDUCT.html).
By participating in this project you agree to abide by its terms.

Also, please read the Terms and Conditions of Use of Speciesplus Data:
https://speciesplus.net/terms-of-use.


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# rcites 1.2.0

* Default branch is now set to `main`.
* `synonyms` are properly formatted (see #65).
* `print()` methods are tested (see #57).
* `curl::curl_escape(x)` is used to encode some URL parts (see #63).
* Consistently uses `message()` for console (see #60).
* Vignettes are not precomputed (see #58).
* Tests now uses `vcr` (now listed in `Suggests`, see #56).
* Classes are now tested properly (see #54).
* `rcites_res()` gains arguments `verbose` and `raw` (see #43 and #62).
* Request status now are reported by `warn_for_status()` rather than by `stop_for_status()`, this prevents fast-failing in batch mode (see #62).


# rcites 1.1.0

* Internal function `rcites_simplify_distributions()` has been re-written to fix a bug that made `spp_distributions()` throw an error for `taxon_id` with only one distribution entry (see #53).
* `spp_*()` functions gain an argument `pause` (see #50 and #51 following the issue reported by @fleurhierink in #49).
* Minor text editions through the documentation.
* Return an empty data frame when there is no listing available for a given species (fix #47 reported by @eveskew).


# rcites 1.0.1

### New features :art:

* rcites now imports [cli](https://github.com/r-lib/cli) :package: to:
  1. clarify notes reported when downloading material;
  2. color titles in the default print methods.


### Bugs fixed :bug:

* Fix a bug that prevented `spp_taxonconcept()` from downloading all the taxon
concepts, see [#42](https://github.com/ropensci/rcites/issues/42).

* Fix a bug that generated infinte recursion error when using non-interactively without token, see [#44](https://github.com/ropensci/rcites/issues/44).


# rcites 1.0.0

### New features :art:

  - `spp_taxonconcept()` now includes an auto-pagination that allows to retrieve
  all entries for queries that have more than 500 results;
  - `spp_taxonconcept()`, `spp_eu_legislation()`, `spp_cites_legislation()` and  `spp_references()` now supports vectors as `taxon_id` argument which allows
  bulk analysis.
  - Functions `spp_*` now returns S3 objects:
    - `spp_taxonconcept()` returns an object of class `spp_taxon`;
    - `spp_cites_legislation()` returns an object of class `spp_cites_leg` or `spp_cites_leg_multi`;
    - `spp_eu_legislation()` returns an object of class `spp_eu_leg` or `spp_eu_leg_multi`;
    - `spp_distributions()` returns an object of class `spp_distr` or `spp_distr_multi`;
    - `spp_references()` returns an object of class `spp_refs` or `spp_refs_multi`;
  - Moreover:
    - `spp_raw` class is defined for all functions including the argument `raw`,
    it returns the parsed output from the API as a list object.   


### New function arguments

  - `taxonomy`, `with_descendants`, `language`, `updated_since`, `per_page`,
  `seq_page`, `raw`, `verbose` in `spp_taxonconcept()`;
  - `scope`, `language` and `raw` to `spp_cites_legislation()`;
  - `scope`, `language` and `raw` to `spp_eu_legislation()`;
  - `language` and `raw` to `spp_distributions()`.

### Function renamed

  - `set_token()` instead of `sppplus_login()`;
  - `spp_taxonconcept()` instead of `sppplus_taxonconcept()`;
  - `spp_cites_legislation()` instead of `taxon_cites_legislation()`;
  - `spp_eu_legislation()` instead of `taxon_eu_legislation()`;
  - `spp_distributions()` instead of `taxon_distribution()`;
  - `spp_references()` instead of `taxon_references()`;
  - `rcites_simplify()` instead of `sppplus_simplify()`;
  - `rcites_` instead of `sppplus_` for helper functions.

### Using [goodpractice](https://github.com/MangoTheCat/goodpractice)

  - use '<-' for assignment instead of '=',
  - omit "Date" in DESCRIPTION,
  - avoid `1:length(...)`, `1:nrow(...)`, `1:ncol(...)`.



# rcites 0.1.0

### NB: this was the first release :one:


### Features :art:

- spp_taxonconcept( ): [access the Speciesplus taxon concept](https://api.speciesplus.net/documentation/v1/taxon_concepts/index.html)
- spp_cites_legislation( ): [access CITES legislation data](https://api.speciesplus.net/documentation/v1/cites_legislation/index.html)
- spp_eu_legislation( ): [access EU legislation data](https://api.speciesplus.net/documentation/v1/eu_legislation/index.html)
- spp_distributions( ): [access a taxon distribution data](https://api.speciesplus.net/documentation/v1/distributions/index.html)
- spp_references( ): [access a listing reference data](https://api.speciesplus.net/documentation/v1/references/index.html)
This is a minor release that mainly improves tests (tests have been reviewed, http requests have been recorded using vcr) and the testing environment (see below). Minor bugs have been squashed along the way. 

## Test environments

* GitHub Actions, Ubuntu 20.04: R-oldrel,
* GitHub Actions, Ubuntu 20.04: R-release,
* GitHub Actions, Ubuntu 20.04: R-devel,
* GitHub Actions, macOS 11.6: R-release,
* GitHub Actions, macOS 11.6: R-devel,
* GitHub Actions, Windows Server 2019 (10.0.17763): R-release,
* GitHub Actions, Windows Server 2019 (10.0.17763): R-devel,
* GitHub Actions, Windows Server 2019 (10.0.17763): R-release,
* GitHub Actions, win-builder (R-release and R-devel),
* local Debian 11 (Kernel: 5.14.0-2-amd64 x86_64), R-4.1.1.


## R CMD check results.

0 ERRORs | 0 WARNINGs | 0 NOTES.


## Downstream dependencies

There are currently no downstream dependencies for this package.
# test token

    Code
      set_token("hackme")
    Message <simpleMessage>
      i Authentication token stored for the session.

# spp_taxonconcept() defaults work

    Code
      print(res)
    Output
      
      -- General info - CITES ($general): --------------------------------------------
      # A tibble: 1 x 8
        id    full_name    author_year    rank  name_status updated_at          active
      * <chr> <chr>        <chr>          <chr> <chr>       <dttm>              <lgl> 
      1 4521  Loxodonta a~ (Blumenbach, ~ SPEC~ A           2021-10-13 13:12:58 TRUE  
      # ... with 1 more variable: cites_listing <chr>
      
      -- Classification ($higher_taxa): ----------------------------------------------
      # A tibble: 1 x 6
        id    kingdom  phylum   class    order       family      
        <chr> <chr>    <chr>    <chr>    <chr>       <chr>       
      1 4521  Animalia Chordata Mammalia Proboscidea Elephantidae
      
      -- Synonyms ($synonyms): -------------------------------------------------------
      # A tibble: 1 x 4
           id full_name          author_year      rank   
        <int> <chr>              <chr>            <chr>  
      1 37069 Loxodonta cyclotis (Matschie, 1900) SPECIES
      
      -- Common names ($common_names): -----------------------------------------------
      # A tibble: 37 x 3
            id name               language
         <int> <chr>              <chr>   
       1  4521 Ndovo              SW      
       2  4521 Tembo              SW      
       3  4521 Haathi             UR      
       4  4521 Elefante           PT      
       5  4521 Slon               RU      
       6  4521 Elefant            NO      
       7  4521 Olifant            NL      
       8  4521 Afrikaanse olifant NL      
       9  4521 Elefante africano  ES      
      10  4521 afrikansk elefant  SV      
      # ... with 27 more rows
      
      Information available: $all_id, $general, $higher_taxa, $accepted_names, $common_names, $synonyms, $cites_listings 

# spp_references() defaults works

    Code
      print(res)
    Output
      
      -- References ($references): ---------------------------------------------------
      # A tibble: 15 x 3
         id    citation                                             is_standard
         <chr> <chr>                                                <chr>      
       1 10265 Anon. 1978. Red data book: Mammalia. IUC [truncated] FALSE      
       2 6344  Barnes, R. F., Agnagna, M., Alers, M. P. [truncated] FALSE      
       3 17013 Blanc, J.J., Thouless, C.R., Hart, J.A., [truncated] FALSE      
       4 6371  Burton, M. P. 1994. Alternative projecti [truncated] FALSE      
       5 6532  Douglas-Hamilton, I. 1987. African Eleph [truncated] FALSE      
       6 6534  Douglas-Hamilton, I. 1987. African Eleph [truncated] FALSE      
       7 6825  Jackson, P. 1982. Elephants and rhinos i [truncated] FALSE      
       8 7224  Meester, J. and Setzer, H. W (eds.) 1974 [truncated] FALSE      
       9 7609  Parker, I. and Amin, M. 1983. Ivory cris [truncated] FALSE      
      10 19397 Parker, I.S.C. and Martin, E.B. 1982. Ho [truncated] FALSE      
      11 19572 Riddle, H.S., Schulte, B.A., Desai, A.A. [truncated] FALSE      
      12 15783 Roca, A. L., Georgiadis, N., Pecon-Slatt [truncated] FALSE      
      13 7852  Said, M. and Change, R. 1994. African el [truncated] FALSE      
      14 43172 Wilson, D.E. and Reeder, D.M. (Eds.) 199 [truncated] TRUE       
      15 17655 Wilson, D.E. and Reeder, D.M. (Eds.). 20 [truncated] FALSE      

# spp_cites_legislation() defaults work

    Code
      print(res)
    Output
      
      -- Cites listings ($cites_listings): -------------------------------------------
      # A tibble: 10 x 6
         id    taxon_concept_id is_current appendix change_type effective_at
         <chr> <chr>            <lgl>      <chr>    <chr>       <chr>       
       1 30344 4521             TRUE       I        +           2017-01-02  
       2 30115 4521             TRUE       II       +           2019-11-26  
       3 32160 4521             TRUE       II       R+          2019-11-26  
       4 32161 4521             TRUE       II       R+          2019-11-26  
       5 32156 4521             TRUE       II       R+          2019-11-26  
       6 32158 4521             TRUE       II       R+          2019-11-26  
       7 32154 4521             TRUE       II       R+          2019-11-26  
       8 32159 4521             TRUE       II       R+          2019-11-26  
       9 32157 4521             TRUE       II       R+          2019-11-26  
      10 32155 4521             TRUE       II       R+          2019-11-26  
      
      -- Cites quotas ($cites_quotas): -----------------------------------------------
      # A tibble: 38 x 10
         id    taxon_concept_id quota publication_date public_display is_current unit 
         <chr> <chr>            <chr> <chr>            <lgl>          <lgl>      <chr>
       1 25337 4521             0     2021-02-03       TRUE           TRUE       <NA> 
       2 25348 4521             0     2021-02-03       TRUE           TRUE       <NA> 
       3 25355 4521             0     2021-02-03       TRUE           TRUE       <NA> 
       4 25358 4521             0     2021-02-03       TRUE           TRUE       <NA> 
       5 25375 4521             0     2021-02-03       TRUE           TRUE       <NA> 
       6 25390 4521             0     2021-02-03       TRUE           TRUE       <NA> 
       7 25414 4521             0     2021-02-03       TRUE           TRUE       <NA> 
       8 25431 4521             0     2021-02-03       TRUE           TRUE       <NA> 
       9 25554 4521             100   2021-02-03       TRUE           TRUE       <NA> 
      10 25555 4521             300   2021-02-03       TRUE           TRUE       <NA> 
      # ... with 28 more rows, and 3 more variables: geo_entity.iso_code2 <chr>,
      #   geo_entity.name <chr>, geo_entity.type <chr>
      Field(s) not printed:  notes, url 
      
      -- Cites suspensions ($cites_suspensions): -------------------------------------
      # A tibble: 13 x 8
         id    taxon_concept_id start_date is_current applies_to_import
         <chr> <chr>            <chr>      <lgl>      <lgl>            
       1 17621 4521             2014-08-11 TRUE       TRUE             
       2 17620 4521             2014-08-11 TRUE       TRUE             
       3 17686 4521             2014-10-10 TRUE       TRUE             
       4 18709 4521             2010-08-16 TRUE       TRUE             
       5 15983 <NA>             2011-01-19 TRUE       FALSE            
       6 22079 <NA>             2018-01-30 TRUE       FALSE            
       7 22076 <NA>             2018-01-22 TRUE       FALSE            
       8 22132 4521             2018-03-19 TRUE       FALSE            
       9 23168 <NA>             2019-07-04 TRUE       FALSE            
      10 24947 4521             2020-05-26 TRUE       TRUE             
      11 25328 4521             2020-12-11 TRUE       FALSE            
      12 26123 <NA>             2021-05-06 TRUE       FALSE            
      13 26275 4521             2018-06-01 TRUE       TRUE             
      # ... with 3 more variables: geo_entity.iso_code2 <chr>, geo_entity.name <chr>,
      #   geo_entity.type <chr>
      Field(s) not printed:  notes, start_notification.name, start_notification.date, start_notification.url 

# spp_distributions() defaults work

    Code
      print(res)
    Output
      
      -- Distributions ($distributions): ---------------------------------------------
      # A tibble: 42 x 5
            id iso_code2 name              type    tags 
         <int> <chr>     <chr>             <chr>   <chr>
       1  1778 ML        Mali              COUNTRY ""   
       2  1923 GQ        Equatorial Guinea COUNTRY ""   
       3  4429 RW        Rwanda            COUNTRY ""   
       4  4491 GH        Ghana             COUNTRY ""   
       5  5628 SD        Sudan             COUNTRY ""   
       6  6724 ET        Ethiopia          COUNTRY ""   
       7  8995 GA        Gabon             COUNTRY ""   
       8 12983 AO        Angola            COUNTRY ""   
       9 15554 CM        Cameroon          COUNTRY ""   
      10 17060 BJ        Benin             COUNTRY ""   
      # ... with 32 more rows
      
      -- References ($references): ---------------------------------------------------
      # A tibble: 146 x 2
            id reference                       
         <int> <chr>                           
       1  1778 Kingdon, J., Happold [truncated]
       2  1923 Basilio, A. 1962. La [truncated]
       3  4429 Jackson, P. 1982. El [truncated]
       4  4429 Monfort, A. 1992. Pr [truncated]
       5  4491 Grubb, P., Jones, T. [truncated]
       6  4491 Jackson, P. 1982. El [truncated]
       7  5628 Hameed, S.M.A. and E [truncated]
       8  6724 Bolton, M. 1973. Not [truncated]
       9  6724 Largen, M. J. and Ya [truncated]
      10  6724 Meester, J. and Setz [truncated]
      # ... with 136 more rows

# spp_distributions() works when no info available

    Code
      print(res)
    Output
      
      -- Distributions ($distributions): ---------------------------------------------
      No records available.
      
      -- References ($references): ---------------------------------------------------
      No records available.

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
# CONTRIBUTING

## We love collaboration!

## Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/rcites/issues)
- be sure to include R session information and a reproducible example (repex).


## Code contributions

### Broad overview of contributing workflow

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rcites.git`
* Make sure to track progress upstream (i.e., on our version of `rcites` at `ropensci/rcites`) by doing `git remote add upstream https://github.com/ropensci/rcites.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs (see Tests below)
* Push up to your account
* Submit a pull request to home base at `ropensci/rcites`

### Tests

To add tests, go to the folder `tests/testthat/`. Tests are generally organized as individual files for each exported function from the package (that is, listed as an export in the `NAMESPACE` file). If you are adding a new exported function, add a new test file. If you are changing an existing function, work in the tests file for that function, unless it doesn't have tests, in which case make a new test file.

The book R packages book provides [a chapter on testing in general](http://r-pkgs.had.co.nz/tests.html). Do consult that first if you aren't familiar with testing in R.

The easiest set up to run tests is from within an R session:

```r
library(devtools)
library(testthat)
# loads the package
load_all()
```

To test an individual test file

```r
test_file("tests/testthat/test-foobar.R")
```

To run all tests

```r
devtools::test()
```

If you are running tests that have `skip_on_cran()` in them, set `Sys.setenv(NOT_CRAN = "true")` prior to running tests.


### Making changes

In addition to changing the code, do make sure to update the documentation if applicable. The R packages book book has a [chapter on documentation](http://r-pkgs.had.co.nz/man.html) you should read if you aren't familiar.

After code and documentation has been changed, update documentation by running either `devtools::document()` or `roxygen2::roxygenise()`.

Make sure if you change what packages or even functions within packages are
imported, most likely add the package to Imports in the DESCRIPTION file.

Be conservative about adding new dependencies.


### Style

* Make sure code, documentation, and comments are no more than 80 characters in width;
* Use `<-` instead of `=` for assignment;
* Always use `snake_case` (and all lowercase) instead of `camelCase`.



## Thanks for contributing!
<!--
If this issue relates to usage of the package, whether a question, bug or similar,
along with your query, please paste your devtools::session_info() or sessionInfo()
into the code block below. If not, delete all this and proceed :)
-->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
title: 'rcites: An R package to access the CITES Speciesplus database'
authors:
- affiliation: 1
  name: Jonas Geschke
  orcid: 0000-0002-5654-9313
- affiliation: 2
  name: Kevin Cazelles
  orcid: 0000-0001-6619-9874
- affiliation: 3
  name: Ignasi Bartomeus
  orcid: 0000-0001-7893-4389
date: "12 August 2018"
bibliography: paper.bib
tags:
- CITES Speciesplus
- taxonomy
- endangered species
- illegal wildlife trade
- species legislation
- species distribution
affiliations:
- index: 1
  name: Museum für Naturkunde Berlin - Leibniz Institute for Research on Evolution and Biodiversity, Berlin, Germany
- index: 2
  name: Department of Integrative Biology, University Of Guelph, Guelph, Ontario, Canada
- index: 3
  name: Estación Biológica de Doñana (EBD-CSIC), Sevilla, Spain
---


# Introduction

The conservation of biodiversity is a complex problem strongly tight to political actions. CITES, the Convention on International Trade in Endangered Species of Wild Fauna and Flora, is a multilateral environmental agreement that was established in 1975 [@CITES_about] and aims to monitor and regulate the trade of endangered species so that their trade does not threaten the survival of the species in the wild [@CITES_about]. CITES is one of the eight main international agreements relevant to biodiversity [@CBD_biodiv-conv] and constitutes a key tool for conservationists, scientists and policy makers.

# The Speciesplus database

In 2013, the UNEP World Conservation Monitoring Centre (UNEP-WCMC) and the CITES Secretariat initiated a partnership funded by UNEP, the European Commission and the CITES Secretariat. Together, they created Speciesplus (or Species+), a comprehensive database of not only CITES listed species and their regulation status within CITES but also the species' status within the EU legislation and the species' status within the Convention on the Conservation of Migratory Species of Wild Animals (CMS) [@Speciesplus_about]. Speciesplus is publicly available at https://speciesplus.net [@UNEP].

# The ``rcites`` R package

With ``rcites`` we provide an R [@R] client to the Speciesplus/CITES Checklist API, giving access to the Speciesplus database. The ability to query the Speciesplus database directly from one of the most used programming languages for data analyses will improve the efficiency and reproducibility of biodiversity conservation analysis workflows.

We provide functions to:

1. access the Speciesplus taxon concept, and thereafter
2. get a species' legislation status, both from CITES and from the European Union,
3. get a species' country-wise distribution range, as listed in Speciesplus, and
4. get the references on which a Speciesplus listing is based.

Overall information about `rcites` and the package vignette can be found at https://docs.ropensci.org/rcites/. The package is available for download from CRAN (stable version; https://cran.r-project.org/package=rcites) and Github (development version; https://github.com/ropensci/rcites).

``rcites`` will support researchers and national authorities in more efficiently dealing with information of endangered species and their legislation status.

Recent publications with data extraction from the Speciesplus database illustrate in what kind of research the package may be of help [@Hinsley:2018; @Robinson:2018]. In the spirit of CITES, research in regard to illegal wildlife trade may be pointed out especially [@Zhou:2016; @Ingram:2015; @Zhou:2015].


# Acknowledgments

Many thanks to the British Ecological Society for bringing us together during the "Ecology Hackathon: Developing R Packages for Accessing, Synthesizing and Analyzing Ecological Data", part of the BES, GfÖ, NecoV and EEF Joint Annual Meeting 2017 "Ecology Across Borders". Also, thanks to our hackathon group and package reviewers for the fruitful discussions and contributions to the package.


# References
---
title: "Study case: the African bush elephant (*Loxodonta africana*)"
author: "rcites team"
date: 10-08-2018
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Study case: the African bush elephant (*Loxodonta africana*)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




# Introduction and setup

In the vignettes "Get started with rcites", we explained how to get a token and
set it up for general access to the CITES Species+ database. Also, we very
briefly introduced to how to code the key features of `rcites`. With this
article, we aim to further introduce to the functionality and workflows of
`rcites`. For this, we use the African bush elephant (*Loxodonta africana*,
hereafter "elephant") as a case study. The elephant not only is a highly
endangered species that is illegally traded globally but also a flagship species
of nature conservation and the logo species of CITES.

We start with a basic set up: we load the package and set the token:


```r
library(rcites)
set_token("8QW6Qgh57sBG2k0gtt")
```

# Retrieve the taxon id

In order to access information about the elephant, we first need to retrieve its
Species+ taxon identifier. For this, we use the `spp_taxonconcept()` function
and the elephant's scientific name, *Loxodonta africana*, as `query_taxon` argument.


```r
elephant_taxonconcept <- spp_taxonconcept(query_taxon = "Loxodonta africana")
```

```
#>  ℹ Retrieving info from page 1 ........................ ✔
```

```r
elephant_taxonconcept
```

```
#>  
#>  ── General info - CITES ($general): ──────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 8
#>    id    full_name          author_year        rank    name_status updated_at          active cites_listing
#>    <chr> <chr>              <chr>              <chr>   <chr>       <dttm>              <lgl>  <chr>        
#>  1 4521  Loxodonta africana (Blumenbach, 1797) SPECIES A           2021-10-13 13:12:58 TRUE   I/II         
#>  
#>  ── Classification ($higher_taxa): ────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 6
#>    id    kingdom  phylum   class    order       family      
#>    <chr> <chr>    <chr>    <chr>    <chr>       <chr>       
#>  1 4521  Animalia Chordata Mammalia Proboscidea Elephantidae
#>  
#>  ── Synonyms ($synonyms): ─────────────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 4
#>       id full_name          author_year      rank   
#>    <int> <chr>              <chr>            <chr>  
#>  1 37069 Loxodonta cyclotis (Matschie, 1900) SPECIES
#>  
#>  ── Common names ($common_names): ─────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 3
#>        id name               language
#>     <int> <chr>              <chr>   
#>   1  4521 Ndovo              SW      
#>   2  4521 Tembo              SW      
#>   3  4521 Haathi             UR      
#>   4  4521 Elefante           PT      
#>   5  4521 Slon               RU      
#>   6  4521 Elefant            NO      
#>   7  4521 Olifant            NL      
#>   8  4521 Afrikaanse olifant NL      
#>   9  4521 Elefante africano  ES      
#>  10  4521 afrikansk elefant  SV      
#>  -------truncated-------
#>  
#>  Information available: $all_id, $general, $higher_taxa, $accepted_names, $common_names, $synonyms, $cites_listings
```


As the first column of the output shows, the taxon identifier of the elephant is 4521. This `taxon_id` will be used for all next function coding.

Beyond the taxon identifier, the output also provides information about the taxon classification and other names, both synonyms and common names if any, in different languages.


## Map the elephant's distribution

Before giving more insights into the legislation status of the elephant, we have a look at where the elephant actually occurs naturally. For this, we can access the elephant's distribution information with the `spp_distributions()` function. Thereafter, we can map the distribution with the help of the `rworldmap` package.


```r
library(rworldmap)

par(las = 1)
elephant_distr <- spp_distributions(taxon_id = "4521",
                                    verbose = FALSE)$distributions

map2 <- joinCountryData2Map(elephant_distr,
                            joinCode="ISO2",
                            nameJoinColumn = "iso_code2",
                            nameCountryColumn = "name")
```

```
#>  42 codes from your data successfully matched countries in the map
#>  0 codes from your data failed to match with a country code in the map
#>  201 codes from the map weren't represented in your data
```

```r
plot(c(-23, 62), c(45, -40),
     type = "n",
     main = "Loxodonta africana",
     xlab = "Longitude",
     ylab = "Latitude")
plot(map2, add = TRUE)
plot(map2[!is.na(map2$iso_code2),], col = "#cba74d", add = TRUE)
```

<img src="img/map-1.png" title="plot of chunk map" alt="plot of chunk map" style="display: block; margin: auto;" />


## Access the legislation status

The functions `spp_cites_legislation()` and `spp_eu_legislation()` provide access to
the legislation status information of the elephant.

First, we have a look at the CITES legislation status:


```r
elephant_cites <- spp_cites_legislation(taxon_id = "4521")
```

```
#>  ℹ Now processing taxon_id '4521'...................... ✔
```

```r
elephant_cites
```

```
#>  
#>  ── Cites listings ($cites_listings): ─────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 6
#>     id    taxon_concept_id is_current appendix change_type effective_at
#>     <chr> <chr>            <lgl>      <chr>    <chr>       <chr>       
#>   1 30344 4521             TRUE       I        +           2017-01-02  
#>   2 30115 4521             TRUE       II       +           2019-11-26  
#>   3 32160 4521             TRUE       II       R+          2019-11-26  
#>   4 32161 4521             TRUE       II       R+          2019-11-26  
#>   5 32156 4521             TRUE       II       R+          2019-11-26  
#>   6 32158 4521             TRUE       II       R+          2019-11-26  
#>   7 32154 4521             TRUE       II       R+          2019-11-26  
#>   8 32159 4521             TRUE       II       R+          2019-11-26  
#>   9 32157 4521             TRUE       II       R+          2019-11-26  
#>  10 32155 4521             TRUE       II       R+          2019-11-26  
#>  
#>  ── Cites quotas ($cites_quotas): ─────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 10
#>     id    taxon_concept_id quota publication_date public_display is_current unit  geo_entity.iso_code2 geo_entity.name 
#>     <chr> <chr>            <chr> <chr>            <lgl>          <lgl>      <chr> <chr>                <chr>           
#>   1 25337 4521             0     2021-02-03       TRUE           TRUE       <NA>  KE                   Kenya           
#>   2 25348 4521             0     2021-02-03       TRUE           TRUE       <NA>  LR                   Liberia         
#>   3 25355 4521             0     2021-02-03       TRUE           TRUE       <NA>  MW                   Malawi          
#>   4 25358 4521             0     2021-02-03       TRUE           TRUE       <NA>  ML                   Mali            
#>   5 25375 4521             0     2021-02-03       TRUE           TRUE       <NA>  MZ                   Mozambique      
#>   6 25390 4521             0     2021-02-03       TRUE           TRUE       <NA>  AO                   Angola          
#>   7 25414 4521             0     2021-02-03       TRUE           TRUE       <NA>  NE                   Niger           
#>   8 25431 4521             0     2021-02-03       TRUE           TRUE       <NA>  NG                   Nigeria         
#>   9 25554 4521             100   2021-02-03       TRUE           TRUE       <NA>  TZ                   United Republic…
#>  10 25555 4521             300   2021-02-03       TRUE           TRUE       <NA>  ZA                   South Africa    
#>  # … with 1 more variable: geo_entity.type <chr>
#>  -------truncated-------
#>  Field(s) not printed:  notes, url 
#>  
#>  ── Cites suspensions ($cites_suspensions): ───────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 8
#>     id    taxon_concept_id start_date is_current applies_to_import geo_entity.iso_code2 geo_entity.name geo_entity.type
#>     <chr> <chr>            <chr>      <lgl>      <lgl>             <chr>                <chr>           <chr>          
#>   1 17621 4521             2014-08-11 TRUE       TRUE              US                   United States … COUNTRY        
#>   2 17620 4521             2014-08-11 TRUE       TRUE              US                   United States … COUNTRY        
#>   3 17686 4521             2014-10-10 TRUE       TRUE              US                   United States … COUNTRY        
#>   4 18709 4521             2010-08-16 TRUE       TRUE              ZW                   Zimbabwe        COUNTRY        
#>   5 15983 <NA>             2011-01-19 TRUE       FALSE             DJ                   Djibouti        COUNTRY        
#>   6 22079 <NA>             2018-01-30 TRUE       FALSE             DJ                   Djibouti        COUNTRY        
#>   7 22076 <NA>             2018-01-22 TRUE       FALSE             LR                   Liberia         COUNTRY        
#>   8 22132 4521             2018-03-19 TRUE       FALSE             AU                   Australia       COUNTRY        
#>   9 23168 <NA>             2019-07-04 TRUE       FALSE             SO                   Somalia         COUNTRY        
#>  10 24947 4521             2020-05-26 TRUE       TRUE              CN                   China           COUNTRY        
#>  -------truncated-------
#>  Field(s) not printed:  notes, start_notification.name, start_notification.date, start_notification.url
```


We can do the same for the elephant's legislation status in the European Union:


```r
elephant_eu <- spp_eu_legislation(taxon_id = "4521")
```

```
#>  ℹ Now processing taxon_id '4521'...................... ✔
```

```r
elephant_eu
```

```
#>  
#>  ── EU listings ($eu_listings): ───────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 2 × 6
#>    id    taxon_concept_id is_current annex change_type effective_at
#>    <chr> <chr>            <lgl>      <chr> <chr>       <chr>       
#>  1 31788 4521             TRUE       A     +           2019-12-14  
#>  2 31876 4521             TRUE       B     +           2019-12-14  
#>  
#>  ── EU decisions ($eu_decisions): ─────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 15
#>     id    taxon_concept_id start_date is_current eu_decision_type… eu_decision_type… geo_entity.iso_c… geo_entity.name 
#>     <chr> <chr>            <chr>      <lgl>      <chr>             <chr>             <chr>             <chr>           
#>   1 26285 4521             2015-04-09 TRUE       Positive          POSITIVE_OPINION  ZW                Zimbabwe        
#>   2 25508 4521             2014-09-03 TRUE       Positive          POSITIVE_OPINION  BW                Botswana        
#>   3 11682 4521             2012-02-23 TRUE       Positive          POSITIVE_OPINION  NA                Namibia         
#>   4 24825 4521             2014-05-28 TRUE       Positive          POSITIVE_OPINION  ZA                South Africa    
#>   5 27017 4521             2015-09-15 TRUE       Positive          POSITIVE_OPINION  ZM                Zambia          
#>   6 27360 4521             2016-06-27 TRUE       Negative          NEGATIVE_OPINION  MZ                Mozambique      
#>   7 35567 4521             2020-03-03 TRUE       <NA>              <NA>              MZ                Mozambique      
#>   8 30377 4521             2017-06-21 TRUE       Negative          NEGATIVE_OPINION  TZ                United Republic…
#>   9 30553 4521             2017-06-21 TRUE       Positive          POSITIVE_OPINION  TZ                United Republic…
#>  10 32143 4521             2019-10-17 TRUE       Suspension (a)    SUSPENSION        CM                Cameroon        
#>  # … with 7 more variables: geo_entity.type <chr>, start_event.name <chr>, start_event.date <chr>, source.code <chr>,
#>  #   source.name <chr>, term.code <chr>, term.name <chr>
#>  Field(s) not printed:  notes, eu_decision_type.description, start_event.url
```


With the combination of `map2` and the legislation data, one might be able to illustrate the elephant's trade directions. This and other use examples of the `rcites` data output will be added bit by bit.


## Access the elephant's Species+ reference data

Last but not least, it is important to identify which references the Species+ information about the elephant is based on. For this, we can access the Species+ reference data with the `spp_references()` function.


```r
elephant_refs <- spp_references(taxon_id = "4521", verbose = FALSE)
elephant_refs
```

```
#>  
#>  ── References ($references): ─────────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 3
#>     id    citation                                             is_standard
#>     <chr> <chr>                                                <chr>      
#>   1 10265 Anon. 1978. Red data book: Mammalia. IUC [truncated] FALSE      
#>   2 6344  Barnes, R. F., Agnagna, M., Alers, M. P. [truncated] FALSE      
#>   3 17013 Blanc, J.J., Thouless, C.R., Hart, J.A., [truncated] FALSE      
#>   4 6371  Burton, M. P. 1994. Alternative projecti [truncated] FALSE      
#>   5 6532  Douglas-Hamilton, I. 1987. African Eleph [truncated] FALSE      
#>   6 6534  Douglas-Hamilton, I. 1987. African Eleph [truncated] FALSE      
#>   7 6825  Jackson, P. 1982. Elephants and rhinos i [truncated] FALSE      
#>   8 7224  Meester, J. and Setzer, H. W (eds.) 1974 [truncated] FALSE      
#>   9 7609  Parker, I. and Amin, M. 1983. Ivory cris [truncated] FALSE      
#>  10 19397 Parker, I.S.C. and Martin, E.B. 1982. Ho [truncated] FALSE      
#>  -------truncated-------
```

```r
dim(elephant_refs$references)
```

```
#>  [1] 15  3
```
---
title: "Get started with rcites"
author: "rcites team"
date: 10-08-2018
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Get started with rcites}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---





# Set up a connection to the Species+/CITES Checklist API

## Get your personal token

To set up a connection to the Species+/CITES Checklist API, an authentication
token is required. Each user should obtain his or her own personal token to run
the code below (see <https://api.speciesplus.net/documentation> for more
details). To obtain a token, sign up on the [Species+ API
website](http://api.speciesplus.net/users/sign_up).

## Set the token

Now, we assume that you already have a token. For illustrative purposes,
we will use the generic token value `8QW6Qgh57sBG2k0gtt` from the API
documentation.

A token is mandatory and needs to be passed to the header of all URL
requests. They are three different ways to set the token in `rcites`:

1. set an environment variable `SPECIESPLUS_TOKEN` in your [`.Renviron`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Startup.html)
file (preferred for frequent users);

2. use `set_token()` to set up the token as a character string (with quotes). If
not entered in the function, the token can be passed without quotes (not as a
character string) after the prompt. This way, the token `SPECIESPLUS_TOKEN` is
interactively set up only for the current R session, meaning a login will be
required for the future R sessions (preferred);


```r
library(rcites)
set_token("8QW6Qgh57sBG2k0gtt")
```

3. use the `token` argument inside the functions, *i.e.* the token is passed manually to each function call.

Note that if you set a wrong token and you wish to set it again interactively,
you must first forget the previous token with `forget_token()`.

# Use of key features

## Retrieve a taxon identifiers with spp_taxonconcept()

### Basic calls to spp_taxonconcept()

In order to efficiently query information from the CITES database, you
first need to retrieve the unique taxon identifier `taxon_id` from the [Species+ taxon concept](https://api.speciesplus.net/documentation/v1/taxon_concepts/index.html).
To do so, you should first call `spp_taxonconcept()` and provide the scientific
name of the taxon you are looking for. Let us start by requesting the identifier of the African bush elephant, *i.e.* *Loxodonta africana*.


```r
res1 <- spp_taxonconcept(query_taxon = "Loxodonta africana")
```

Note that if you have decide to set your token using the third option, then the code should look like the one below:


```r
res1 <- spp_taxonconcept(query_taxon = "Loxodonta africana", token = "8QW6Qgh57sBG2k0gtt")
```

### Object `spp_taxonconcept`

`res1` is an S3 object of class `spp_taxon`:


```r
attributes(res1)
```

```
#>  $names
#>  [1] "all_id"         "general"        "higher_taxa"    "accepted_names" "common_names"   "synonyms"      
#>  [7] "cites_listings"
#>  
#>  $class
#>  [1] "spp_taxon"
#>  
#>  $taxonomy
#>  [1] "CITES"
```

that contains information sorted into several data frames (see `?spp_taxonconcept`
for further details):



```r
res1
```

```
#>  
#>  ── General info - CITES ($general): ──────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 8
#>    id    full_name          author_year        rank    name_status updated_at          active cites_listing
#>    <chr> <chr>              <chr>              <chr>   <chr>       <dttm>              <lgl>  <chr>        
#>  1 4521  Loxodonta africana (Blumenbach, 1797) SPECIES A           2021-10-13 13:12:58 TRUE   I/II         
#>  
#>  ── Classification ($higher_taxa): ────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 6
#>    id    kingdom  phylum   class    order       family      
#>    <chr> <chr>    <chr>    <chr>    <chr>       <chr>       
#>  1 4521  Animalia Chordata Mammalia Proboscidea Elephantidae
#>  
#>  ── Synonyms ($synonyms): ─────────────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 4
#>       id full_name          author_year      rank   
#>    <int> <chr>              <chr>            <chr>  
#>  1 37069 Loxodonta cyclotis (Matschie, 1900) SPECIES
#>  
#>  ── Common names ($common_names): ─────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 3
#>        id name               language
#>     <int> <chr>              <chr>   
#>   1  4521 Ndovo              SW      
#>   2  4521 Tembo              SW      
#>   3  4521 Haathi             UR      
#>   4  4521 Elefante           PT      
#>   5  4521 Slon               RU      
#>   6  4521 Elefant            NO      
#>   7  4521 Olifant            NL      
#>   8  4521 Afrikaanse olifant NL      
#>   9  4521 Elefante africano  ES      
#>  10  4521 afrikansk elefant  SV      
#>  -------truncated-------
#>  
#>  Information available: $all_id, $general, $higher_taxa, $accepted_names, $common_names, $synonyms, $cites_listings
```

For some taxa, there are more than one taxon identifier available. In `general` only
*active* identifiers are listed, but the full list of identifiers are
available in `all_id`:


```r
res3 <- spp_taxonconcept(query = "Amazilia versicolor")
```

```
#>  ℹ Retrieving info from page 1 ........................ ✔
```

```r
res3$general
```

```
#>  # A tibble: 1 × 8
#>    id    full_name           author_year      rank    name_status updated_at          active cites_listing
#>  * <chr> <chr>               <chr>            <chr>   <chr>       <dttm>              <lgl>  <chr>        
#>  1 3210  Amazilia versicolor (Vieillot, 1818) SPECIES A           2015-05-07 15:10:59 TRUE   II
```

```r
res3$all_id
```

```
#>  # A tibble: 2 × 7
#>    id    full_name           author_year      rank    name_status updated_at          active
#>    <chr> <chr>               <chr>            <chr>   <chr>       <dttm>              <lgl> 
#>  1 3210  Amazilia versicolor (Vieillot, 1818) SPECIES A           2015-05-07 15:10:59 TRUE  
#>  2 65789 Amazilia versicolor (Vieillot, 1818) SPECIES S           2016-09-23 15:30:46 FALSE
```

Also, if the taxon is not listed, a warning message should come up.


```r
res4 <- spp_taxonconcept(query = "Homo Sapiens")
```

```
#>  ℹ Retrieving info from page 1 ........................ ✔
```

```
#>  Warning in spp_taxonconcept(query = "Homo Sapiens"): Taxon not listed.
```

### Custom calls to spp_taxonconcept()

`spp_taxonconcept()` includes several arguments to retrieve a specific subset
of information (see `?spp_taxonconcept` for more details):


```r
args('spp_taxonconcept')
```

```
#>  function (query_taxon, taxonomy = "CITES", with_descendants = FALSE, 
#>      language = NULL, updated_since = NULL, per_page = 500, pages = NULL, 
#>      raw = FALSE, token = NULL, verbose = TRUE, pause = 1, ...) 
#>  NULL
```

Most importantly, the argument `taxonomy` allows a selection between the two
databases (CITES or CMS):


```r
spp_taxonconcept(query = "Amazilia versicolor", taxonomy = "CMS")
```

```
#>  ℹ Retrieving info from page 1 ........................ ✔
```

```
#>  Warning in spp_taxonconcept(query = "Amazilia versicolor", taxonomy = "CMS"): Taxon not listed.
```

```
#>  NULL
```

```r
spp_taxonconcept(query = "Loxodonta africana", taxonomy = "CMS")
```

```
#>  ℹ Retrieving info from page 1 ........................ ✔
```

```
#>  
#>  ── General info - CMS ($general): ────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 7
#>    id    full_name          author_year        rank    name_status updated_at          active
#>    <chr> <chr>              <chr>              <chr>   <chr>       <dttm>              <lgl> 
#>  1 11691 Loxodonta africana (Blumenbach, 1797) SPECIES A           2021-05-06 12:44:13 TRUE  
#>  
#>  ── Classification ($higher_taxa): ────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 6
#>    id    kingdom  phylum   class    order       family      
#>    <chr> <chr>    <chr>    <chr>    <chr>       <chr>       
#>  1 11691 Animalia Chordata Mammalia Proboscidea Elephantidae
#>  
#>  ── Synonyms ($synonyms): ─────────────────────────────────────────────────────────────────────────────────────────────
#>  No records available.
#>  
#>  ── Common names ($common_names): ─────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 3
#>        id name               language
#>     <int> <chr>              <chr>   
#>   1 11691 Tembo              SW      
#>   2 11691 Ndovo              SW      
#>   3 11691 Haathi             UR      
#>   4 11691 Elefante           PT      
#>   5 11691 Slon               RU      
#>   6 11691 Elefant            NO      
#>   7 11691 Olifant            NL      
#>   8 11691 Afrikaanse olifant NL      
#>   9 11691 Elefante africano  ES      
#>  10 11691 afrikansk elefant  SV      
#>  -------truncated-------
#>  
#>  Information available: $all_id, $general, $higher_taxa, $accepted_names, $common_names, $synonyms
```

`language` and `updated_since` are convenient filters for
the written language of common names (must be a
two-letters code, see [ISO 3166-1 alpha-2](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2))
and the last update of the entries, respectively:


```r
spp_taxonconcept(query_taxon = "Loxodonta africana", language = 'EN',
  updated_since = "2016-01-01")
```

```
#>  ℹ Retrieving info from page 1 ........................ ✔
```

```
#>  
#>  ── General info - CITES ($general): ──────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 8
#>    id    full_name          author_year        rank    name_status updated_at          active cites_listing
#>    <chr> <chr>              <chr>              <chr>   <chr>       <dttm>              <lgl>  <chr>        
#>  1 4521  Loxodonta africana (Blumenbach, 1797) SPECIES A           2021-10-13 13:12:58 TRUE   I/II         
#>  
#>  ── Classification ($higher_taxa): ────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 6
#>    id    kingdom  phylum   class    order       family      
#>    <chr> <chr>    <chr>    <chr>    <chr>       <chr>       
#>  1 4521  Animalia Chordata Mammalia Proboscidea Elephantidae
#>  
#>  ── Synonyms ($synonyms): ─────────────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 1 × 4
#>       id full_name          author_year      rank   
#>    <int> <chr>              <chr>            <chr>  
#>  1 37069 Loxodonta cyclotis (Matschie, 1900) SPECIES
#>  
#>  ── Common names ($common_names): ─────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 2 × 3
#>       id name                      language
#>    <int> <chr>                     <chr>   
#>  1  4521 African Elephant          EN      
#>  2  4521 African Savannah Elephant EN      
#>  
#>  Information available: $all_id, $general, $higher_taxa, $accepted_names, $common_names, $synonyms, $cites_listings
```

## Retrieve a taxon identifiers available data

In order to use the four `spp_*` functions, one needs to use the
*active taxon identifier* of a given species. For instance, for the two
species we used as examples above we use the value indicated in the table below:

|name                | taxon identifier|
|:-------------------|:----------------|
|Loxodonta africana  |       4521      |
|Amazilia versicolor |       3210      |



### CITES legislation data

First, we can retrieve current CITES appendix listings and reservations, CITES quotas, and CITES suspensions for a given taxon concept.


```r
spp_cites_legislation(taxon_id = "4521", verbose = FALSE)
```

```
#>  
#>  ── Cites listings ($cites_listings): ─────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 6
#>     id    taxon_concept_id is_current appendix change_type effective_at
#>     <chr> <chr>            <lgl>      <chr>    <chr>       <chr>       
#>   1 30344 4521             TRUE       I        +           2017-01-02  
#>   2 30115 4521             TRUE       II       +           2019-11-26  
#>   3 32160 4521             TRUE       II       R+          2019-11-26  
#>   4 32161 4521             TRUE       II       R+          2019-11-26  
#>   5 32156 4521             TRUE       II       R+          2019-11-26  
#>   6 32158 4521             TRUE       II       R+          2019-11-26  
#>   7 32154 4521             TRUE       II       R+          2019-11-26  
#>   8 32159 4521             TRUE       II       R+          2019-11-26  
#>   9 32157 4521             TRUE       II       R+          2019-11-26  
#>  10 32155 4521             TRUE       II       R+          2019-11-26  
#>  
#>  ── Cites quotas ($cites_quotas): ─────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 10
#>     id    taxon_concept_id quota publication_date public_display is_current unit  geo_entity.iso_code2 geo_entity.name 
#>     <chr> <chr>            <chr> <chr>            <lgl>          <lgl>      <chr> <chr>                <chr>           
#>   1 25337 4521             0     2021-02-03       TRUE           TRUE       <NA>  KE                   Kenya           
#>   2 25348 4521             0     2021-02-03       TRUE           TRUE       <NA>  LR                   Liberia         
#>   3 25355 4521             0     2021-02-03       TRUE           TRUE       <NA>  MW                   Malawi          
#>   4 25358 4521             0     2021-02-03       TRUE           TRUE       <NA>  ML                   Mali            
#>   5 25375 4521             0     2021-02-03       TRUE           TRUE       <NA>  MZ                   Mozambique      
#>   6 25390 4521             0     2021-02-03       TRUE           TRUE       <NA>  AO                   Angola          
#>   7 25414 4521             0     2021-02-03       TRUE           TRUE       <NA>  NE                   Niger           
#>   8 25431 4521             0     2021-02-03       TRUE           TRUE       <NA>  NG                   Nigeria         
#>   9 25554 4521             100   2021-02-03       TRUE           TRUE       <NA>  TZ                   United Republic…
#>  10 25555 4521             300   2021-02-03       TRUE           TRUE       <NA>  ZA                   South Africa    
#>  # … with 1 more variable: geo_entity.type <chr>
#>  -------truncated-------
#>  Field(s) not printed:  notes, url 
#>  
#>  ── Cites suspensions ($cites_suspensions): ───────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 8
#>     id    taxon_concept_id start_date is_current applies_to_import geo_entity.iso_code2 geo_entity.name geo_entity.type
#>     <chr> <chr>            <chr>      <lgl>      <lgl>             <chr>                <chr>           <chr>          
#>   1 17621 4521             2014-08-11 TRUE       TRUE              US                   United States … COUNTRY        
#>   2 17620 4521             2014-08-11 TRUE       TRUE              US                   United States … COUNTRY        
#>   3 17686 4521             2014-10-10 TRUE       TRUE              US                   United States … COUNTRY        
#>   4 18709 4521             2010-08-16 TRUE       TRUE              ZW                   Zimbabwe        COUNTRY        
#>   5 15983 <NA>             2011-01-19 TRUE       FALSE             DJ                   Djibouti        COUNTRY        
#>   6 22079 <NA>             2018-01-30 TRUE       FALSE             DJ                   Djibouti        COUNTRY        
#>   7 22076 <NA>             2018-01-22 TRUE       FALSE             LR                   Liberia         COUNTRY        
#>   8 22132 4521             2018-03-19 TRUE       FALSE             AU                   Australia       COUNTRY        
#>   9 23168 <NA>             2019-07-04 TRUE       FALSE             SO                   Somalia         COUNTRY        
#>  10 24947 4521             2020-05-26 TRUE       TRUE              CN                   China           COUNTRY        
#>  -------truncated-------
#>  Field(s) not printed:  notes, start_notification.name, start_notification.date, start_notification.url
```

### EU legislation data

Similarly, we can also retrieve current EU annex listings, SRG opinions, and EU
suspensions with `spp_eu_legislation`. Both legislation functions have a `scope`
argument that sets the time scope of legislation and take one value among
`current`, `historic` and `all` (default is set to `current`). For instance, one
can get all information pertaining to EU annex listing for *Amazilia versicolor*
with the following command line:


```r
spp_eu_legislation(taxon_id = "3210", scope = "all", verbose = FALSE)
```

```
#>  
#>  ── EU listings ($eu_listings): ───────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 6
#>     id    taxon_concept_id is_current annex change_type effective_at
#>     <chr> <chr>            <lgl>      <chr> <chr>       <chr>       
#>   1 17611 3210             FALSE      B     +           1997-06-01  
#>   2 17610 3210             FALSE      B     +           2000-12-18  
#>   3 17609 3210             FALSE      B     +           2003-08-30  
#>   4 17608 3210             FALSE      B     +           2005-08-22  
#>   5 17607 3210             FALSE      B     +           2008-04-11  
#>   6 17606 3210             FALSE      B     +           2009-05-22  
#>   7 17605 3210             FALSE      B     +           2010-08-15  
#>   8 17604 3210             FALSE      B     +           2012-02-14  
#>   9 17603 3210             FALSE      B     +           2012-12-15  
#>  10 17602 3210             FALSE      B     +           2013-08-10  
#>  -------truncated-------
#>  
#>  ── EU decisions ($eu_decisions): ─────────────────────────────────────────────────────────────────────────────────────
#>  No records available.
```


### Distribution data

Distribution data at the country level is also available for a given taxon concept:


```r
spp_distributions("4521", verbose = FALSE)
```

```
#>  
#>  ── Distributions ($distributions): ───────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 5
#>        id iso_code2 name              type    tags 
#>     <int> <chr>     <chr>             <chr>   <chr>
#>   1  1778 ML        Mali              COUNTRY ""   
#>   2  1923 GQ        Equatorial Guinea COUNTRY ""   
#>   3  4429 RW        Rwanda            COUNTRY ""   
#>   4  4491 GH        Ghana             COUNTRY ""   
#>   5  5628 SD        Sudan             COUNTRY ""   
#>   6  6724 ET        Ethiopia          COUNTRY ""   
#>   7  8995 GA        Gabon             COUNTRY ""   
#>   8 12983 AO        Angola            COUNTRY ""   
#>   9 15554 CM        Cameroon          COUNTRY ""   
#>  10 17060 BJ        Benin             COUNTRY ""   
#>  -------truncated-------
#>  
#>  ── References ($references): ─────────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 2
#>        id reference                       
#>     <int> <chr>                           
#>   1  1778 Kingdon, J., Happold [truncated]
#>   2  1923 Basilio, A. 1962. La [truncated]
#>   3  4429 Jackson, P. 1982. El [truncated]
#>   4  4429 Monfort, A. 1992. Pr [truncated]
#>   5  4491 Grubb, P., Jones, T. [truncated]
#>   6  4491 Jackson, P. 1982. El [truncated]
#>   7  5628 Hameed, S.M.A. and E [truncated]
#>   8  6724 Bolton, M. 1973. Not [truncated]
#>   9  6724 Largen, M. J. and Ya [truncated]
#>  10  6724 Meester, J. and Setz [truncated]
#>  -------truncated-------
```

### References

Finally, we can retrieve all available references for a given taxa.


```r
spp_references("4521", verbose = FALSE)
```

```
#>  
#>  ── References ($references): ─────────────────────────────────────────────────────────────────────────────────────────
#>  # A tibble: 10 × 3
#>     id    citation                                             is_standard
#>     <chr> <chr>                                                <chr>      
#>   1 10265 Anon. 1978. Red data book: Mammalia. IUC [truncated] FALSE      
#>   2 6344  Barnes, R. F., Agnagna, M., Alers, M. P. [truncated] FALSE      
#>   3 17013 Blanc, J.J., Thouless, C.R., Hart, J.A., [truncated] FALSE      
#>   4 6371  Burton, M. P. 1994. Alternative projecti [truncated] FALSE      
#>   5 6532  Douglas-Hamilton, I. 1987. African Eleph [truncated] FALSE      
#>   6 6534  Douglas-Hamilton, I. 1987. African Eleph [truncated] FALSE      
#>   7 6825  Jackson, P. 1982. Elephants and rhinos i [truncated] FALSE      
#>   8 7224  Meester, J. and Setzer, H. W (eds.) 1974 [truncated] FALSE      
#>   9 7609  Parker, I. and Amin, M. 1983. Ivory cris [truncated] FALSE      
#>  10 19397 Parker, I.S.C. and Martin, E.B. 1982. Ho [truncated] FALSE      
#>  -------truncated-------
```
---
title: "Bulk analysis with rcites"
author: "rcites team"
date: 15-08-2018
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Bulk analysis with rcites}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
editor_options:
  chunk_output_type: console
---




## Broad taxon concept queries

If you want to query all taxa, you can use `spp_taxonconcept()` with
`query_taxon = ""` (assuming your token is already set up):


```r
res_cms <- spp_taxonconcept("", taxonomy = "CMS") #slow
```

```
#>  ℹ Retrieving info from page 1 ........................ ✔
#>  ℹ 10 pages available, retrieving info from 9 more
#>  ℹ Retrieving info from page 2 ........................ ✔
#>  ℹ Retrieving info from page 3 ........................ ✔
#>  ℹ Retrieving info from page 4 ........................ ✔
#>  ℹ Retrieving info from page 5 ........................ ✔
#>  ℹ Retrieving info from page 6 ........................ ✔
#>  ℹ Retrieving info from page 7 ........................ ✔
#>  ℹ Retrieving info from page 8 ........................ ✔
#>  ℹ Retrieving info from page 9 ........................ ✔
#>  ℹ Retrieving info from page 10 ....................... ✔
```

```r
dim(res_cms$general)
```

```
#>  [1] 2541    7
```

Alternatively, you can retrieve, for example, the first three pages of results
returned by the API.


```r
res_cites <- spp_taxonconcept("", page = 1:2)
```

```
#>  ℹ Retrieving info from page 1 ........................ ✔
#>  ℹ 167 pages available, retrieving info from 1 more
#>  ℹ Retrieving info from page 2 ........................ ✔
```

```r
dim(res_cites$general)
```

```
#>  [1] 1000    8
```

## Retrieving information for a set of taxon_concept ID

All `spp_` functions (i.e. `spp_distributions()`, `spp_eu_legislation()`,
`spp_cites_legislation()` and `spp_references()`) can handle a vector of
taxon_id which allows bulk analysis.
Below we exemplify this feature for the four functions.

### spp_distributions()


```r
vc_txn <- c('4521', '3210', '10255')
res1 <- spp_distributions(taxon_id = vc_txn, verbose = FALSE)
## Number of countries concerned per taxon ID
table(res1$distributions$taxon_id)
```

```
#>  
#>  10255  3210  4521 
#>     15     8    42
```



### spp_cites_legislation()


```r
res2 <- spp_cites_legislation(taxon_id = vc_txn, verbose = FALSE)
res2$cites_listings
```

```
#>  # A tibble: 12 × 7
#>     taxon_id id    taxon_concept_id is_current appendix change_type effective_at
#>     <chr>    <chr> <chr>            <lgl>      <chr>    <chr>       <chr>       
#>   1 4521     30344 4521             TRUE       I        +           2017-01-02  
#>   2 4521     30115 4521             TRUE       II       +           2019-11-26  
#>   3 4521     32160 4521             TRUE       II       R+          2019-11-26  
#>   4 4521     32161 4521             TRUE       II       R+          2019-11-26  
#>   5 4521     32156 4521             TRUE       II       R+          2019-11-26  
#>   6 4521     32158 4521             TRUE       II       R+          2019-11-26  
#>   7 4521     32154 4521             TRUE       II       R+          2019-11-26  
#>   8 4521     32159 4521             TRUE       II       R+          2019-11-26  
#>   9 4521     32157 4521             TRUE       II       R+          2019-11-26  
#>  10 4521     32155 4521             TRUE       II       R+          2019-11-26  
#>  11 3210     4661  3210             TRUE       II       +           1987-10-22  
#>  12 10255    4645  10255            TRUE       I        +           2005-01-12
```

### spp_eu_legislation()




```r
res3 <- spp_eu_legislation(taxon_id = vc_txn, verbose = FALSE)
res3$eu_listings
```

```
#>  # A tibble: 4 × 7
#>    taxon_id id    taxon_concept_id is_current annex change_type effective_at
#>    <chr>    <chr> <chr>            <lgl>      <chr> <chr>       <chr>       
#>  1 4521     31788 4521             TRUE       A     +           2019-12-14  
#>  2 4521     31876 4521             TRUE       B     +           2019-12-14  
#>  3 3210     30578 3210             TRUE       B     +           2019-12-14  
#>  4 10255    30516 10255            TRUE       A     +           2019-12-14
```


### spp_references()




```r
res4 <- spp_references(taxon_id = vc_txn, verbose = FALSE)
## Number of references per taxon ID
table(res4$references$taxon_id)
```

```
#>  
#>  10255  3210  4521 
#>      1     3    15
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spp_taxonconcept.R
\name{spp_taxonconcept}
\alias{spp_taxonconcept}
\title{Get taxon concepts for a search term}
\usage{
spp_taxonconcept(
  query_taxon,
  taxonomy = "CITES",
  with_descendants = FALSE,
  language = NULL,
  updated_since = NULL,
  per_page = 500,
  pages = NULL,
  raw = FALSE,
  token = NULL,
  verbose = TRUE,
  pause = 1,
  ...
)
}
\arguments{
\item{query_taxon}{a character string containing the query (e.g. species).
Scientific taxa only (max 255 characters).}

\item{taxonomy}{filter taxon concepts by taxonomy, accepts either 'CITES' or
'CMS' as its value. Default sets to 'CITES'.}

\item{with_descendants}{a logical. Should the search by name be broadened to
include higher taxa?}

\item{language}{filter languages returned for common names. Value should be a
vector of character strings including one or more country codes (two-letters
country code
\href{https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}{ISO 3166-1 alpha-2}).
Default is set to \code{NULL}, showing all available languages.}

\item{updated_since}{a timestamp. Only entries updated after (and including)
this timestamp will be pulled.}

\item{per_page}{an integer that indicates how many objects are returned per
page for paginated responses. Default set to 500 which is the maximum.}

\item{pages}{a vector of integer that contains page numbers. Default is
set to \code{NULL}, i.e. all pages are accessed.}

\item{raw}{a logical. Should raw data be returned?}

\item{token}{a character string containing the authentification token, see
\url{https://api.speciesplus.net/documentation}. Default is set to
\code{NULL} and requires the environment variable \code{SPECIESPLUS_TOKEN} to be set
directly in \code{Renviron}. Alternatively, \code{set_token()} can be used to set
\code{SPECIESPLUS_TOKEN} for the current session.}

\item{verbose}{a logical. Should extra information be reported on progress?}

\item{pause}{a duration (in second) to suspend execution for (see
\code{\link[=Sys.sleep]{Sys.sleep()}}). This was added cause the web API returns a 404 error too many
requests in a short time interval.}

\item{...}{Further named parameters, see \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
If \code{raw = TRUE}, then a object of class \code{spp_raw} is returned, which is
a list of lists. If \code{raw = FALSE}, then an object of class \code{spp_taxon} is
returned, it is a collection of seven data frames:
\enumerate{
\item \code{all_id}: general information for all entries, including non-active taxon
concepts,
\item \code{general}: includes general information for active taxon concepts,
\item \code{higher_taxa}: includes taxonomy information,
\item \code{accepted_names}: list of accepted names (only for synonyms),
\item \code{common_names}: list of common names (only for accepted names),
\item \code{synonyms}: list of synonyms (only for accepted names),
\item \code{cites_listing}: list of current CITES listings with annotations
(missing if \code{taxonomy == 'CMS'}).
}
}
\description{
Retrieve the taxon concept of a specific taxon (scientific name).
}
\examples{
\donttest{
res1 <- spp_taxonconcept(query_taxon = 'Loxodonta africana')
res2 <- spp_taxonconcept(query_taxon = 'Amazilia versicolor', raw = TRUE)
res3 <- spp_taxonconcept(query_taxon = '', taxonomy = 'CMS', pages = c(1, 3),
 language = 'EN', config = httr::progress())
res4 <- spp_taxonconcept(query_taxon = '', per_page = 20, pages = 44)
}
}
\references{
\url{https://api.speciesplus.net/documentation/v1/taxon_concepts/index.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spp_cites_legislation.R
\name{spp_cites_legislation}
\alias{spp_cites_legislation}
\title{Get current CITES appendix listings and reservations}
\usage{
spp_cites_legislation(
  taxon_id,
  scope = "current",
  language = "en",
  raw = FALSE,
  token = NULL,
  verbose = TRUE,
  pause = 1,
  ...
)
}
\arguments{
\item{taxon_id}{a vector of character strings containing species' taxon
concept identifiers (see \code{\link[=spp_taxonconcept]{spp_taxonconcept()}}).}

\item{scope}{vector of character strings indicating the time scope of
legislation, values are taken among \code{current}, \code{historic} and \code{all}. Default
is \code{current}.}

\item{language}{vector of character strings indicating the language for the
text of legislation notes, values are taken among \code{en} (English),
\code{fr} (French) and \code{es} (Spanish). Default is \code{en}.}

\item{raw}{a logical. Should raw data be returned?}

\item{token}{a character string containing the authentification token, see
\url{https://api.speciesplus.net/documentation}. Default is set to
\code{NULL} and requires the environment variable \code{SPECIESPLUS_TOKEN} to be
set directly in \code{Renviron}. Alternatively, \code{\link[=set_token]{set_token()}} can
be used to set \code{SPECIESPLUS_TOKEN} for the current session.}

\item{verbose}{a logical. Should extra information be reported on progress?}

\item{pause}{a duration (in second) to suspend execution for (see
\code{\link[=Sys.sleep]{Sys.sleep()}}). This was added cause the web API returns a 404 error too many
requests in a short time interval.}

\item{...}{Further named parameters, see \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
If \code{raw} is set to \code{TRUE} then an object of class \code{spp_raw} (or
\code{spp_raw_multi} if \code{length(taxon_id) > 1}) is returned which is essentially
a list of lists (see option \code{as = 'parsed'} in \code{\link[httr:content]{httr::content()}}).
Otherwise, an object of class \code{spp_cites_leg} (or \code{spp_cites_leg_multi} if
\code{length(taxon_id)>1}) is returned which is a list of three data frames:
\enumerate{
\item \code{cites_listings}: lists CITES annex listings EU suspensions,
\item \code{cites_quotas}: lists CITES quotas,
\item \code{cites_suspensions}: lists CITES suspensions.
}
}
\description{
Retrieve current CITES appendix listings and reservations, CITES quotas, and
CITES suspensions for a given taxon concept.
}
\examples{
\donttest{
res1 <- spp_cites_legislation(taxon_id = 4521)
res2 <- spp_cites_legislation(taxon_id = c('4521', '3210', '10255'))
res3 <- spp_cites_legislation(taxon_id = 4521, scope = 'all',
verbose = FALSE, config=httr::verbose())
res4 <- spp_cites_legislation(taxon_id = 4521, language = 'fr')
}
}
\references{
\url{https://api.speciesplus.net/documentation/v1/cites_legislation/index.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spp_references.R
\name{spp_references}
\alias{spp_references}
\title{Get references for a given taxon concept}
\usage{
spp_references(
  taxon_id,
  raw = FALSE,
  token = NULL,
  verbose = TRUE,
  pause = 1,
  ...
)
}
\arguments{
\item{taxon_id}{a vector of character strings containing species' taxon
concept identifiers (see \code{\link[=spp_taxonconcept]{spp_taxonconcept()}}).}

\item{raw}{a logical. Should raw data be returned?}

\item{token}{a character string containing the authentification token, see
\url{https://api.speciesplus.net/documentation}. Default is set to
\code{NULL} and requires the environment variable \code{SPECIESPLUS_TOKEN} to be
set directly in \code{Renviron}. Alternatively, \code{\link[=set_token]{set_token()}} can be used to set
\code{SPECIESPLUS_TOKEN} for the current session.}

\item{verbose}{a logical. Should extra information be reported on progress?}

\item{pause}{a duration (in second) to suspend execution for (see
\code{\link[=Sys.sleep]{Sys.sleep()}}). This was added cause the web API returns a 404 error too many
requests in a short time interval.}

\item{...}{Further named parameters, see \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
If \code{raw} is set to \code{TRUE} then an object of class \code{spp_raw} (or
\code{spp_raw_multi} if \code{length(taxon_id) > 1}) is returned which is essentially
a list of lists (see option \code{as = 'parsed'} in \code{\link[httr:content]{httr::content()}}).
Otherwise, an object of class \code{spp_refs} (or \code{spp_refs_multi} if
\code{length(taxon_id) > 1}) is returned which is a list of one
data frame:
\itemize{
\item \code{references} that includes the identifier of the reference and the
corresponding citation.
}
}
\description{
Retrieve available references for a given taxon concept.
}
\examples{
\donttest{
res1 <- spp_references(taxon_id = '4521')
res2 <- spp_references(c('4521', '3210', '10255'))
res3 <- spp_references(taxon_id = '4521', raw = TRUE, verbose = FALSE,
 config = httr::progress())
}
}
\references{
\url{https://api.speciesplus.net/documentation/v1/references/index.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.spp.R
\name{print.spp}
\alias{print.spp}
\alias{print.spp_raw}
\alias{print.spp_raw_multi}
\alias{print.spp_cites_leg}
\alias{print.spp_cites_leg_multi}
\alias{print.spp_distr}
\alias{print.spp_distr_multi}
\alias{print.spp_eu_leg}
\alias{print.spp_eu_leg_multi}
\alias{print.spp_refs}
\alias{print.spp_refs_multi}
\alias{print.spp_taxon}
\title{Print methods for objects of class \verb{spp_raw*}.}
\usage{
\method{print}{spp_raw}(x, ...)

\method{print}{spp_raw_multi}(x, ...)

\method{print}{spp_cites_leg}(x, ...)

\method{print}{spp_cites_leg_multi}(x, ...)

\method{print}{spp_distr}(x, ...)

\method{print}{spp_distr_multi}(x, ...)

\method{print}{spp_eu_leg}(x, ...)

\method{print}{spp_eu_leg_multi}(x, ...)

\method{print}{spp_refs}(x, ...)

\method{print}{spp_refs_multi}(x, ...)

\method{print}{spp_taxon}(x, ...)
}
\arguments{
\item{x}{an object of class \verb{spp_raw*}.}

\item{...}{ignored.}
}
\value{
The JSON result.
}
\description{
Print the outputs of a Species+ API call.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{package}
\name{rcites}
\alias{rcites}
\alias{rcites-package}
\title{rcites}
\description{
A programmatic interface to the Species+ \url{https://speciesplus.net/} database
via the Species+/CITES Checklist API \url{https://api.speciesplus.net/}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/rcites/}
  \item \url{https://github.com/ropensci/rcites}
  \item Report bugs at \url{https://github.com/ropensci/rcites/issues}
}

}
\author{
\strong{Maintainer}: Kevin Cazelles \email{kevin.cazelles@gmail.com} (\href{https://orcid.org/0000-0001-6619-9874}{ORCID})

Authors:
\itemize{
  \item Jonas Geschke \email{jonas.e.geschke@gmail.com} (\href{https://orcid.org/0000-0002-5654-9313}{ORCID})
  \item Ignasi Bartomeus (\href{https://orcid.org/0000-0001-7893-4389}{ORCID})
}

Other contributors:
\itemize{
  \item Jonathan Goldenberg [contributor]
  \item Marie-Bé Leduc [contributor]
  \item Yasmine Verzelen [contributor]
  \item Noam Ross [reviewer]
  \item Margaret Siple [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spp_distributions.R
\name{spp_distributions}
\alias{spp_distributions}
\title{Get distributions data available for a given taxon concept}
\usage{
spp_distributions(
  taxon_id,
  language = "en",
  raw = FALSE,
  token = NULL,
  verbose = TRUE,
  pause = 1,
  ...
)
}
\arguments{
\item{taxon_id}{a vector of character strings containing species'
taxon concept identifiers (see \code{\link[=spp_taxonconcept]{spp_taxonconcept()}}).}

\item{language}{vector of character strings indicating the language for the
names of distributions, values are taken among \code{en} (English),
\code{fr} (French) and \code{es} (Spanish). Default is \code{en}.}

\item{raw}{a logical. Should raw data be returned?}

\item{token}{a character string containing the authentification token, see
\url{https://api.speciesplus.net/documentation}. Default is set to
\code{NULL} and requires the environment variable \code{SPECIESPLUS_TOKEN} to be
set directly in \code{Renviron}. Alternatively, \code{set_token()} can
be used to set \code{SPECIESPLUS_TOKEN} for the current session.}

\item{verbose}{a logical. Should extra information be reported on progress?}

\item{pause}{a duration (in second) to suspend execution for (see
\code{\link[=Sys.sleep]{Sys.sleep()}}). This was added cause the web API returns a 404 error too many
requests in a short time interval.}

\item{...}{Further named parameters, see \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
If \code{raw} is set to \code{TRUE} then an object of class \code{spp_raw} (or
\code{spp_raw_multi} if \code{length(taxon_id)>1}) is returned which is essentially
a list of lists (see option \code{as = 'parsed'} in \code{\link[httr:content]{httr::content()}}).
Otherwise, an object of class \code{spp_distr} (or \code{spp_distr_multi} if
\code{length(taxon_id) > 1}) is returned which is a list of two data frames:
\enumerate{
\item \code{distributions}: lists distributions for a given taxon concept,
\item \code{references}: lists the corresponding references.
In case \code{taxon_id} includes several elements
}
}
\description{
Retrieve distributions data available for a given taxon concept for which the
the taxon identifier is known.
}
\examples{
\donttest{
 res1 <- spp_distributions(taxon_id = '4521')
 res2 <- spp_distributions(taxon_id = c('4521', '3210', '10255'))
 res3 <- spp_distributions(taxon_id = '4521', raw = TRUE)
 res4 <- spp_distributions(taxon_id = '4521', language = 'fr',
 verbose = FALSE, config = httr::progress())
}
}
\references{
\url{https://api.speciesplus.net/documentation/v1/distributions/index.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_token.R
\name{set_token}
\alias{set_token}
\alias{forget_token}
\title{Login helper function}
\usage{
set_token(token = NULL)

forget_token()
}
\arguments{
\item{token}{a character string (with quotes) containing your token. If
\code{NULL}, then the token can be passed without quotes (not as character
string) after a prompt.}
}
\description{
Set and forget the authentification token for the current session.
}
\section{Functions}{
\itemize{
\item \code{set_token}: set the environment variable \code{SPECIESPLUS_TOKEN}.

\item \code{forget_token}: forget the environment variable \code{SPECIESPLUS_TOKEN}.
}}

\examples{
\dontrun{
 # NB: the token below is not working
 set_token('8QW6Qgh57sBG2k0gtt')
 # interactively
 set_token()
}
}
\references{
\url{https://api.speciesplus.net/documentation}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spp_eu_legislation.R
\name{spp_eu_legislation}
\alias{spp_eu_legislation}
\title{Get current EU annex listings, SRG opinions, and EU suspensions}
\usage{
spp_eu_legislation(
  taxon_id,
  scope = "current",
  language = "en",
  raw = FALSE,
  token = NULL,
  verbose = TRUE,
  pause = 1,
  ...
)
}
\arguments{
\item{taxon_id}{a vector of character strings containing species' taxon
concept identifiers (see \code{\link[=spp_taxonconcept]{spp_taxonconcept()}}).}

\item{scope}{vector of character strings indicating the time scope of
legislation, values are taken among \code{current}, \code{historic} and \code{all}.
Default is set to \code{current}.}

\item{language}{vector of character strings indicating the language for the
text of legislation notes, values are taken among \code{en} (English),
\code{fr} (French) and \code{es} (Spanish). Default is \code{en}.}

\item{raw}{a logical. Should raw data be returned?}

\item{token}{a character string containing the authentification token, see
\url{https://api.speciesplus.net/documentation}. Default is set to
\code{NULL} and requires the environment variable \code{SPECIESPLUS_TOKEN} to be
set directly in \code{Renviron}. Alternatively, \code{\link[=set_token]{set_token()}} can
be used to set \code{SPECIESPLUS_TOKEN} for the current session.}

\item{verbose}{a logical. Should extra information be reported on progress?}

\item{pause}{a duration (in second) to suspend execution for (see
\code{\link[=Sys.sleep]{Sys.sleep()}}). This was added cause the web API returns a 404 error too many
requests in a short time interval.}

\item{...}{Further named parameters, see \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
If \code{raw} is set to \code{TRUE} then an object of class \code{spp_raw} (or
\code{spp_raw_multi} if \code{length(taxon_id)>1}) is returned which is essentially
a list of lists (see option \code{as = 'parsed'} in \code{\link[httr:content]{httr::content()}}).
Otherwise, an object of class \code{spp_eu_leg} (or \code{spp_eu_leg_multi} if
\code{length(taxon_id)>1}) is returned which is a list of two data frames:
\enumerate{
\item \code{eu_listings}: lists EU annex listings EU suspensions,
\item \code{eu_decisions}: lists EU decisions
}
}
\description{
Retrieve current EU annex listings, SRG opinions, and EU suspensions for a
given taxon concept (identifier must be known).
}
\examples{
\donttest{
res1 <- spp_eu_legislation(taxon_id = '4521')
res2 <- spp_eu_legislation(taxon_id = c('4521', '3210', '10255'))
res3 <- spp_eu_legislation(taxon_id = '4521', scope = 'historic')
res4 <- spp_eu_legislation(taxon_id = '4521', scope = 'all', language='fr',
 verbose = FALSE, config=httr::verbose())
}
}
\references{
\url{https://api.speciesplus.net/documentation/v1/eu_legislation/index.html}
}
