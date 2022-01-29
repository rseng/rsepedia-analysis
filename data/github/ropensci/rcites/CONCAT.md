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
