
# rfishbase <img src="man/figures/logo.svg" align="right" alt="" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/ropensci/rfishbase/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rfishbase/actions)
[![cran
checks](https://cranchecks.info/badges/worst/rfishbase)](https://cranchecks.info/pkgs/rfishbase)
[![Coverage
status](https://codecov.io/gh/ropensci/rfishbase/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/rfishbase?branch=master)
[![Onboarding](https://badges.ropensci.org/137_status.svg)](https://github.com/ropensci/software-review/issues/137)
[![CRAN
status](https://www.r-pkg.org/badges/version/rfishbase)](https://cran.r-project.org/package=rfishbase)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rfishbase)](https://github.com/r-hub/cranlogs.app)
<!-- badges: end -->

Welcome to `rfishbase 4`. This is the fourth rewrite of the original
`rfishbase` package described in [Boettiger et
al. (2012)](https://doi.org/10.1111/j.1095-8649.2012.03464.x).

-   `rfishbase 1.0` relied on parsing of XML pages served directly from
    Fishbase.org.  
-   `rfishbase 2.0` relied on calls to a ruby-based API, `fishbaseapi`,
    that provided access to SQL snapshots of about 20 of the more
    popular tables in FishBase or SeaLifeBase.
-   `rfishbase 3.0` side-stepped the API by making queries which
    directly downloaded compressed csv tables from a static web host.
    This substantially improved performance a reliability, particularly
    for large queries. The release largely remained backwards compatible
    with 2.0, and added more tables.
-   `rfishbase 4.0` extends the static model and interface. Static
    tables are distributed in parquet and accessed through a
    provenance-based identifier. While old functions are retained, a new
    interface is introduced to provide easy access to all fishbase
    tables.

We welcome any feedback, issues or questions that users may encounter
through our issues tracker on GitHub:
<https://github.com/ropensci/rfishbase/issues>

## Installation

``` r
remotes::install_github("ropensci/rfishbase")
```

``` r
library("rfishbase")
library("dplyr") # convenient but not required
```

## Getting started

## Generic table interface

All fishbase tables can be accessed by name using the `fb_tbl()`
function:

``` r
fb_tbl("ecosystem")
```

    # A tibble: 155,792 × 18
       autoctr E_CODE EcosystemRefno Speccode Stockcode Status CurrentPresence
         <int>  <int>          <int>    <int>     <int> <chr>  <chr>          
     1       1      1          50628      549       565 native Present        
     2       2      1            189      552       568 native Present        
     3       3      1            189      554       570 native Present        
     4       4      1          79732      873       889 native Present        
     5       5      1           5217      948       964 native Present        
     6       7      1          39852      956       972 native Present        
     7       8      1          39852      957       973 native Present        
     8       9      1          39852      958       974 native Present        
     9      10      1            188     1526      1719 native Present        
    10      11      1            188     1626      1819 native Present        
    # … with 155,782 more rows, and 11 more variables: Abundance <chr>,
    #   LifeStage <chr>, Remarks <chr>, Entered <int>, Dateentered <dttm>,
    #   Modified <int>, Datemodified <dttm>, Expert <int>, Datechecked <dttm>,
    #   WebURL <chr>, TS <dttm>

You can see all the tables using `fb_tables()` to see a list of all the
table names (specify `sealifebase` if desired). Careful, there are a lot
of them! The fishbase databases have grown a lot in the decades, and
were not intended to be used directly by most end-users, so you may have
considerable work to determine what’s what. Keep in mind that many
variables can be estimated in different ways (e.g. trophic level), and
thus may report different values in different tables. Also note that
species is name (or SpecCode) is not always the primary key for a table
– many tables are specific to stocks or even individual samples, and
some tables are reference lists that are not species focused at all, but
meant to be joined to other tables (`faoareas`, etc). Compare tables
against what you see on fishbase.org, or ask on our issues forum for
advice!

``` r
fish <- c("Oreochromis niloticus", "Salmo trutta")

fb_tbl("species") %>% 
  mutate(sci_name = paste(Genus, Species)) %>%
  filter(sci_name %in% fish) %>% 
  select(sci_name, FBname, Length)
```

    # A tibble: 2 × 3
      sci_name              FBname       Length
      <chr>                 <chr>         <dbl>
    1 Oreochromis niloticus Nile tilapia     60
    2 Salmo trutta          Sea trout       140

## SeaLifeBase

SeaLifeBase.org is maintained by the same organization and largely
parallels the database structure of Fishbase. As such, almost all
`rfishbase` functions can instead be instructed to address the

``` r
fb_tbl("species", "sealifebase")
```

    # A tibble: 97,220 × 109
       SpecCode Genus  Species Author  SpeciesRefNo FBname FamCode Subfamily GenCode
          <int> <chr>  <chr>   <chr>          <int> <chr>    <int> <chr>       <int>
     1    32307 Aapto… americ… (Pilsb…           19 <NA>       815 <NA>        27838
     2    32306 Aapto… brinto… Newman…        81749 <NA>       815 <NA>        27838
     3    32308 Aapto… callis… (Pilsb…           19 <NA>       815 <NA>        27838
     4    32304 Aapto… leptod… Newman…           19 <NA>       815 <NA>        27838
     5    32305 Aapto… trider… Newman…           19 <NA>       815 <NA>        27838
     6    51720 Aaptos aaptos  (Schmi…           19 <NA>      2630 <NA>         9253
     7   165941 Aaptos bergma… de Lau…       108813 <NA>      2630 <NA>         9253
     8   105687 Aaptos ciliata (Wilso…         3477 <NA>      2630 <NA>         9253
     9   139407 Aaptos duchas… (Topse…        85482 <NA>      2630 <NA>         9253
    10   130070 Aaptos laxosu… (Solla…        81108 <NA>      2630 <NA>         9253
    # … with 97,210 more rows, and 100 more variables: TaxIssue <int>,
    #   Remark <chr>, PicPreferredName <chr>, PicPreferredNameM <chr>,
    #   PicPreferredNameF <chr>, PicPreferredNameJ <chr>, Source <chr>,
    #   AuthorRef <int>, SubGenCode <int>, Fresh <int>, Brack <int>,
    #   Saltwater <int>, Land <int>, BodyShapeI <chr>, DemersPelag <chr>,
    #   AnaCat <chr>, MigratRef <int>, DepthRangeShallow <int>,
    #   DepthRangeDeep <int>, DepthRangeRef <int>, DepthRangeComShallow <int>, …

## Versions and importing all tables

By default, tables are downloaded the first time they are used.
`rfishbase` defaults to download the latest available snapshot; be aware
that the most recent snapshot may be months behind the latest data on
fishbase.org. Check available releases:

``` r
available_releases()
```

    [1] "21.06" "19.04"

Users can trigger a one-time download of all fishbase tables (or a list
of desired tables) using `fb_import()`. This will ensure later use of
any function can operate smoothly even when no internet connection is
available. Any table already downloaded will not be re-downloaded.
(Note: `fb_import()` also returns a remote duckdb database connection to
the tables, for users who prefer to work with the remote data objects.)

``` r
conn <- fb_import()
```

## Low-memory environments

If you have very limited RAM (e.g. &lt;= 2 GB available) it may be
helpful to use `fishbase` tables in remote form by setting
`collect = FALSE`. This allows the tables to remain on disk, while the
user is still able to use almost all `dplyr` functions (see the `dbplyr`
vignette). Once the table is appropriately subset, the user will need to
call `dplyr::collect()` to use generic non-dplyr functions, such as
plotting commmands.

``` r
fb_tbl("occurrence", collect = FALSE)
```

    # Source:   table<occurrence> [?? x 106]
    # Database: duckdb_connection
       catnum2 OccurrenceRefNo SpecCode Syncode Stockcode GenusCol       SpeciesCol 
         <int>           <int>    <int>   <int>     <int> <chr>          <chr>      
     1   34424           36653      227   22902       241 "Megalops"     "cyprinoid…
     2   95154           45880       NA      NA        NA ""             ""         
     3   97606           45880       NA      NA        NA ""             ""         
     4  100025           45880     5520   25676      5809 "Johnius"      "belangeri…
     5   98993           45880     5676   16650      5969 "Chromis"      "retrofasc…
     6   99316           45880      454   23112       468 "Drepane"      "punctata" 
     7   99676           45880     5388  145485      5647 "Gymnothorax"  "boschi"   
     8   99843           45880    16813  119925     15264 "Hemiramphus"  "balinensi…
     9  100607           45880     8288   59635      8601 "Ostracion"    "rhinorhyn…
    10  101529           45880       NA      NA        NA "Scomberoides" "toloo-par…
    # … with more rows, and 99 more variables: ColName <chr>, PicName <chr>,
    #   CatNum <chr>, URL <chr>, Station <chr>, Cruise <chr>, Gazetteer <chr>,
    #   LocalityType <chr>, WaterDepthMin <dbl>, WaterDepthMax <dbl>,
    #   AltitudeMin <int>, AltitudeMax <int>, LatitudeDeg <int>, LatitudeMin <dbl>,
    #   NorthSouth <chr>, LatitudeDec <dbl>, LongitudeDeg <int>,
    #   LongitudeMIn <dbl>, EastWest <chr>, LongitudeDec <dbl>, Accuracy <chr>,
    #   Salinity <chr>, LatitudeTo <dbl>, LongitudeTo <dbl>, LatitudeDegTo <int>, …

## Interactive RStudio pane

RStudio users can also browse all fishbase tables interactively in the
RStudio connection browser by using the function `fisbase_pane()`. Note
that this function will first download a complete set of the fishbase
tables.

## Backwards compatibility

`rfishbase` 4.0 tries to maintain as much backwards compatibility as
possible with rfishbase 3.0. Because parquet preserves native data
types, some encoded types may differ from earlier versions. As before,
these are not always the native type – e.g. fishbase encodes some
boolean (logical TRUE/FALSE) values as integer (-1, 0) or character
types. Use `as.logical()` to coerce into the appropriate type in that
case.

Toggling between fishbase and sealifebase servers using an environmental
variable, `FISHBASE_API`, is now deprecated.

Note that fishbase will store downloaded files by hash in the app
directory, given by `db_dir()`. The default location can be set by
configuring the desired path in the environmental variable,
`FISHBASE_HOME`.

------------------------------------------------------------------------

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
NEWS
====

For more fine-grained list of changes or to report a bug, consult 

* [The issues log](https://github.com/ropensci/rfishbase/issues)
* [The commit log](https://github.com/ropensci/rfishbase/commits/master)

Versioning
----------

Releases will be numbered with the following semantic versioning format:

<major>.<minor>.<patch>

And constructed with the following guidelines:

* Breaking backward compatibility bumps the major (and resets the minor 
  and patch)
* New additions without breaking backward compatibility bumps the minor 
  (and resets the patch)
* Bug fixes and misc changes bumps the patch

For more information on SemVer, please visit http://semver.org/.

v 4.0.0
--------

* Major upgrade that introduces content-identifier based downloads and parquet-backed database interface.
  Provides improved access to all tables and improved performance.
* See README or details

v 3.1.9
-------

* avoid GitHub API calls to determine versions.

v 3.1.8
-------

use `collect()` on taxa_tbl()

v 3.1.7
-------

- Avoid needless warning about `arrange()`

v 3.1.6
-------

- ensure any test needing internet connection is "fails gracefully" with
  no warnings or errors whenever tests are run without internet or the 
  CRAN test environment where internet connectivity may be unreliable.

v 3.1.5
-------

- replace rappdirs with base R tools

v 3.1.4
-------

- Uses `arkdb` with `duckdb` as database backend
- resolve compatibility issues

v 3.0.5
--------

- `validate_names()` has been rewritten [#170]

v 3.0.4
--------

- use latest version by default [#164]

v 3.0.3
------

- fix bug in sealifebase name resolution [#154]

v 3.0.2
--------

- fix missing function endpoint for `diet_items()` [#158]


v 3.0.1
--------

- patch for upcoming R release with staged installs
- patch common_names to allow omitting species_list for full table (#156)

v 3.0.0
------

v 3.0.0 accesses a new static API for `fishbase` with in-memory
memoization that significantly improves performance, particularly
for large queries.  

- Functions no longer have default limits on returns, so pagination
  is never involved -- all functions now return full set of available
  results.  
- Almost all functions can be called without arguments (e.g. without
  a species list) to return the complete record of the requested table.
- Various minor issues in some functions have been resolved, see 
  <https://github.com/ropensci/rfishbase/issues/> for details.


v2.2.0 
-------

(not released to CRAN, rolled into 3.0 release)

* bugfix for `validate_names()` ([#121](https://github.com/ropensci/rfishbase/issues/121))
* bugfix for `faoareas()` ([#123](https://github.com/ropensci/rfishbase/issues/123))
* add `genetics()` endpoint ([#122](https://github.com/ropensci/rfishbase/issues/122))
* add `taxonomy()` endpoint ([#126](https://github.com/ropensci/rfishbase/issues/126))

v2.1.2   (2017-04-19)
------

* bugfix avoid spurious warning when using http instead of https API
* bugfix to for taxa table as used on sealifebase


v2.1.1
-------


* bugfix for endpoints with inconsistent spelling of SpecCode column 
(e.g. maturity, diet, ecosystem, morphology, morphometrics, popchar, poplf).
* Now properly queries by input species list
* Minor bug fixes for issues #77 #82 #83 #89 #93 #94 #95 #97 #99 #100 (see https://github.com/ropensci/rfishbase/issues)

v2.1.0
------

* improve sealifebase operations
* Updates for optimized version of fishbaseapi
* More stable tests

v2.0.3
------

* Use `dontrun` instead of `\donttest` in the documentation.  Apparently CRAN has decided to run donttest blocks in the examples, which can fail on their test servers when transferring data files over a network connection (which is why they were marked `donttest` in the first place.)

v2.0.2
------

* bugfix for package unit tests on some platforms

v2.0.1
------

First official release of the new rfishbase package, a non-backwards
compatible replacement to all earlier versions of FishBase. See the
package vignette for a more detailed overview and introduction.
Dear CRAN Maintainers,

This release provides the changes detailed in NEWS.md.
I have checked reverse dependencies.
I have also checked on M1 Mac and Solaris platforms using `rhub`.



Sincerely,

Carl Boettiger<!--- Provide a general summary of your changes in the Title above -->

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

### Please contribute!

We love collaboration.

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/rfishbase/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rfishbase.git`
* Make sure to track progress upstream (i.e., on our version of `rfishbase` at `ropensci/rfishbase`) by doing `git remote add upstream https://github.com/ropensci/rfishbase.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/rfishbase`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Thanks for contributing!
<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- Please ensure you are using the most recent version of `rfishbase` before reporting a bug.-->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
output: github_document
---

# rfishbase <img src="man/figures/logo.svg" align="right" alt="" width="120" />


<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/rfishbase/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rfishbase/actions)
[![cran checks](https://cranchecks.info/badges/worst/rfishbase)](https://cranchecks.info/pkgs/rfishbase)
[![Coverage status](https://codecov.io/gh/ropensci/rfishbase/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/rfishbase?branch=master)
[![Onboarding](https://badges.ropensci.org/137_status.svg)](https://github.com/ropensci/software-review/issues/137)
[![CRAN status](https://www.r-pkg.org/badges/version/rfishbase)](https://cran.r-project.org/package=rfishbase)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rfishbase)](https://github.com/r-hub/cranlogs.app)
<!-- badges: end -->



Welcome to `rfishbase 4`. This is the fourth rewrite of the original `rfishbase` package described in [Boettiger et al. (2012)](https://doi.org/10.1111/j.1095-8649.2012.03464.x).   

- `rfishbase 1.0` relied on parsing of XML pages served directly from Fishbase.org.  
- `rfishbase 2.0` relied on calls to a ruby-based API, `fishbaseapi`, that provided access to SQL snapshots of about 20 of the more popular tables in FishBase or SeaLifeBase.
- `rfishbase 3.0` side-stepped the API by making queries which directly downloaded compressed csv tables from a static web host. This substantially improved performance a reliability, particularly for large queries. The release largely remained backwards compatible with 2.0, and added more tables.
- `rfishbase 4.0` extends the static model and interface. Static tables are distributed in parquet and accessed through a provenance-based identifier. While old functions are retained, a new interface is introduced to provide easy access to all fishbase tables.

We welcome any feedback, issues or questions that users may encounter through our issues tracker on GitHub: <https://github.com/ropensci/rfishbase/issues>




```{r include=FALSE}
knitr::opts_chunk$set(warning=FALSE, comment=NA)
```


## Installation



```{r message=FALSE, warning=FALSE, results="hide", eval=FALSE}
remotes::install_github("ropensci/rfishbase")
```


```{r message=FALSE, warning=FALSE, results="hide"}
library("rfishbase")
library("dplyr") # convenient but not required
```

## Getting started

## Generic table interface

All fishbase tables can be accessed by name using the `fb_tbl()` function:

```{r}
fb_tbl("ecosystem")
```


You can see all the tables using `fb_tables()` to see a list of all the table names (specify `sealifebase` if desired). Careful, there are a lot of them! The fishbase databases have grown a lot in the decades, and were not intended to be used directly by most end-users, so you may have considerable work to determine what's what. Keep in mind that many variables can be estimated in different ways (e.g. trophic level), and thus may report different values in different tables.  Also note that species is name (or SpecCode) is not always the primary key for a table -- many tables are specific to stocks or even individual samples, and some tables are reference lists that are not species focused at all, but meant to be joined to other tables (`faoareas`, etc).  Compare tables against what you see on fishbase.org, or ask on our issues forum for advice!


```{r}
fish <- c("Oreochromis niloticus", "Salmo trutta")

fb_tbl("species") %>% 
  mutate(sci_name = paste(Genus, Species)) %>%
  filter(sci_name %in% fish) %>% 
  select(sci_name, FBname, Length)

```


## SeaLifeBase

SeaLifeBase.org is maintained by the same organization and largely parallels the database structure of Fishbase. As such, almost all `rfishbase` functions can instead be instructed to address the 


```{r}
fb_tbl("species", "sealifebase")
```

## Versions and importing all tables

By default, tables are downloaded the first time they are used.  `rfishbase` defaults to download the latest available snapshot; be aware that the most recent snapshot may be months behind the latest data on fishbase.org. Check available releases:

```{r}
available_releases()
```

Users can trigger a one-time download of all fishbase tables (or a list of desired tables) using `fb_import()`. This will ensure later use of any function can operate smoothly even when no internet connection is available. Any table already downloaded will not be re-downloaded. (Note: `fb_import()` also returns a remote duckdb database connection to the tables, for users who prefer to work with the remote data objects.) 

```{r}
conn <- fb_import()
```
 
## Low-memory environments

If you have very limited RAM (e.g. <= 2 GB available) it may be helpful to use `fishbase` tables in remote form by setting `collect = FALSE`.  This allows the tables to remain on disk, while the user is still able to use almost all `dplyr` functions (see the `dbplyr` vignette).  Once the table is appropriately subset, the user will need to call `dplyr::collect()` to use generic non-dplyr functions, such as plotting commmands.  

```{r}
fb_tbl("occurrence", collect = FALSE)
```


## Interactive RStudio pane

RStudio users can also browse all fishbase tables interactively in the RStudio connection browser by using the function `fisbase_pane()`.  Note that this function will first download a complete set of the fishbase tables.  

## Backwards compatibility


`rfishbase` 4.0 tries to maintain as much backwards compatibility as possible with rfishbase 3.0. Because parquet preserves native data types, some encoded types may differ from earlier versions. As before, these are not always the native type -- e.g. fishbase encodes some boolean (logical TRUE/FALSE) values as integer (-1, 0) or character types. Use `as.logical()` to coerce into the appropriate type in that case. 

Toggling between fishbase and sealifebase servers using an environmental variable, `FISHBASE_API`, is now deprecated.  

Note that fishbase will store downloaded files by hash in the app directory, given by `db_dir()`.  The default location can be set by configuring the desired path in the environmental variable, `FISHBASE_HOME`. 




-----------

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.


[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)


---

---


```{r}
library("rfishbase")
```

For several species simultaneously, yes this is possible.  Just give a string with multiple species scientific names to the function, e.g.:

```{r}
     df <- length_weight(c("Oreochromis niloticus", "Salmo trutta"))
```

Very long lists of species can take quite a while to run.  rfishbase can help you build a long list of species names; for instance, if you want every species in the genus Labroides: 

```{r}
     fish <- species_list(Genus = "Labroides") 
     df <- length_weight(fish)
```


Note that the scientific name must be the one recognized by fishbase; see the help page ?synonyms for a function to help get the right scientific name.  

If you are interested in only certain columns of the table, just give those columns in the `fields` argument; for example:

```{r}
    lw <- length_weight(c("Oreochromis niloticus", "Salmo trutta"), 
                           fields = c("SpecCode", "a", "b"))
```

(Almost all functions work the same way, taking a list of species and optionally a list of the desired columns).  Note that the species identity is in SpecCode, and we have many different values of a and b taken from different specimens of the same species. We can use the function `speciesnames` to replace the SpecCode with the actual species name, e.g. 

```{r}
    lw$SpecCode = speciesnames(lw$SpecCode)  
```

Trophic position data is available from the ecology table I believe, so you could do:

```{r}
    troph <- ecology(c("Oreochromis niloticus", "Salmo trutta"),   
                       fields = c("SpecCode", "FoodTroph", 
                                  "FoodSeTroph", "DietTroph", "DietSeTroph"))
```

You may want to use `merge` to combine these tables.  Note that the ecology table has only one row per species (SpecCode), while the length_weight table has many rows.  

So, what about the Bayesian method?  Right now, you'd have to that manually.  I can look into adding a feature based on the R code they provide, but there are many challenges to doing this responsibly.  I'm not an expert in this area, but there are many ways to go from the sample of "a" and "b" values to a single "best estimate" for the species.  Right now rfishbase is aiming at providing the raw data to better enable researchers to use their own analysis methods and compare different analysis methods, rather than provide button click defaults. (sorry, didn't mean to lecture).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_taxa.R
\docType{data}
\name{sealifebase}
\alias{sealifebase}
\title{A table of all the the species found in SeaLifeBase, including taxonomic
classification and the Species Code (SpecCode) by which the species is
identified in SeaLifeBase}
\description{
A table of all the the species found in SeaLifeBase, including taxonomic
classification and the Species Code (SpecCode) by which the species is
identified in SeaLifeBase
}
\author{
Carl Boettiger \email{carl@ropensci.org}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/species_names.R
\name{species_names}
\alias{species_names}
\title{species names}
\usage{
species_names(
  codes,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db()
)
}
\arguments{
\item{codes}{a vector of speccodes (e.g. column from a table)}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}
}
\value{
A character vector of species names for the SpecCodes
}
\description{
returns species names given FishBase's SpecCodes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fb_tbl.R
\name{fb_tbl}
\alias{fb_tbl}
\title{Access a fishbase or sealifebase table}
\usage{
fb_tbl(
  tbl,
  server = c("fishbase", "sealifebase"),
  version = "latest",
  db = fb_conn(server, version),
  collect = TRUE
)
}
\arguments{
\item{tbl}{table name, as it appears in the database. See [fb_tables()]
for a list.}

\item{server}{fishbase or sealifebase}

\item{version}{release version}

\item{db}{A cachable duckdb database connection}

\item{collect}{should we return an in-memory table? Generally best to leave
as TRUE unless RAM is too limited.  A remote table can be used with most
dplyr functions (filter, select, joins, etc) to further refine.}
}
\description{
Access a fishbase or sealifebase table
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
fb_tbl("species")
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{ecosystem}
\alias{ecosystem}
\title{ecosystem}
\usage{
ecosystem(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species ecosystems data
}
\description{
ecosystem
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
ecosystem("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fishbase.R
\docType{package}
\name{rfishbase-package}
\alias{rfishbase-package}
\alias{rfishbase}
\title{The new R interface to Fishbase, v2.0}
\description{
A programmatic interface to FishBase, re-written
based on an accompanying 'RESTful' API. Access tables describing over 30,000
species of fish, their biology, ecology, morphology, and more. This package also
supports experimental access to SeaLifeBase data, which contains
nearly 200,000 species records for all types of aquatic life not covered by
FishBase.'
}
\author{
Carl Boettiger \email{carl@ropensci.org}

Scott Chamberlain \email{scott@ropensci.org}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trophic_ecology.R
\name{fooditems}
\alias{fooditems}
\title{fooditems}
\usage{
fooditems(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species fooditems
}
\description{
fooditems
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
fooditems("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/fishbasethe_food_items_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synonyms.R
\name{synonyms}
\alias{synonyms}
\title{synonyms}
\usage{
synonyms(
  species_list = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
A table with information about the synonym. Will generally be only a single
row if a species name is given.  If a FishBase SpecCode is given, all synonyms matching
that SpecCode are shown, and the table indicates which one is Valid for FishBase. This may
or may not match the valid name for Catalog of Life (Col), also shown in the table. See examples for details.
}
\description{
Check for alternate versions of a scientific name
}
\details{
For further information on fields returned, see:
http://www.fishbase.org/manual/english/fishbasethe_synonyms_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population_dynamics.R
\name{popchar}
\alias{popchar}
\title{popchar}
\usage{
popchar(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\description{
Table of maximum length (Lmax), weight (Wmax) and age (tmax)
}
\details{
See references for official documentation.  From FishBase.org:
This table presents information on maximum length (Lmax), 
weight (Wmax) and age (tmax) from various localities where a species
occurs. The largest values from this table are also entered in the
SPECIES table. The POPCHAR table also indicates whether the Lmax,
Wmax and tmax values or various combinations thereof refer to the
same individual fish.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
popchar("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/fishbasethe_popchar_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trophic_ecology.R
\name{predators}
\alias{predators}
\title{predators}
\usage{
predators(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of predators
}
\description{
predators
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
predators("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/fishbasethe_predators_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/species_list.R
\name{species_list}
\alias{species_list}
\title{species_list}
\usage{
species_list(
  Class = NULL,
  Order = NULL,
  Family = NULL,
  Subfamily = NULL,
  Genus = NULL,
  Species = NULL,
  SpecCode = NULL,
  SuperClass = NULL,
  server = getOption("FISHBASE_API", FISHBASE_API)
)
}
\arguments{
\item{Class}{Request all species in this taxonomic Class}

\item{Order}{Request all species in this taxonomic Order}

\item{Family}{Request all species in this taxonomic Family}

\item{Subfamily}{Request all species in this taxonomic SubFamily}

\item{Genus}{Request all species in this taxonomic Genus}

\item{Species}{Request all species in this taxonomic Species}

\item{SpecCode}{Request species name of species matching this SpecCode}

\item{SuperClass}{Request all species of this Superclass}

\item{server}{fishbase or sealifebase}
}
\description{
Return the a species list given a taxonomic group
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
## All species in the Family 
  species_list(Family = 'Scaridae')
## All species in the Genus 
  species_list(Genus = 'Labroides')
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_db.R
\name{db_dir}
\alias{db_dir}
\title{show fishbase directory}
\usage{
db_dir()
}
\description{
show fishbase directory
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trophic_ecology.R
\name{ecology}
\alias{ecology}
\title{ecology}
\usage{
ecology(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species ecology data
}
\description{
ecology
}
\details{
By default, will only return one entry (row) per species.  Increase limit to
get multiple returns for different stocks of the same species, though often data is either
identical to the first or simply missing in the additional stocks.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

\dontrun{
ecology("Oreochromis niloticus")

## trophic levels and standard errors for a list of species
ecology(c("Oreochromis niloticus", "Salmo trutta"),
        fields=c("SpecCode", "FoodTroph", "FoodSeTroph", "DietTroph", "DietSeTroph"))
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/fishbasethe_ecology_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refs.R
\name{references}
\alias{references}
\title{references}
\usage{
references(
  codes = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  ...
)
}
\arguments{
\item{codes}{One or more Fishbase reference numbers, matching the RefNo 
field}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a tibble (data.frame) of reference data
}
\description{
references
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
references(codes = 1)
references(codes = 1:6)
references(codes = 1:6, fields = c('Author', 'Year', 'Title'))
references() # all references
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio_pane.R
\name{fishbase_pane}
\alias{fishbase_pane}
\title{Open database connection pane in RStudio}
\usage{
fishbase_pane()
}
\description{
This function launches the RStudio "Connection" pane to interactively
explore the database.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

if (!is.null(getOption("connectionObserver"))) fishbase_pane()
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{occurrence}
\alias{occurrence}
\title{occurrence}
\usage{
occurrence()
}
\description{
occurrence
}
\details{
THE OCCURRENCE TABLE HAS BEEN DROPPED BY FISHBASE - THIS
FUNCTION NOW RETURNS A STOP MESSAGE.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{countrysubref}
\alias{countrysubref}
\title{countrysubref}
\usage{
countrysubref(
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  ...
)
}
\arguments{
\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\description{
return a table of countrysubref
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
countrysubref()
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population_dynamics.R
\name{length_length}
\alias{length_length}
\alias{popll}
\title{length_length}
\usage{
length_length(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of lengths
}
\description{
return a table of lengths
}
\details{
This table contains relationships for the conversion of one length type to another for over 8,000
species of fish, derived from different publications, e.g. Moutopoulos and Stergiou (2002) and
Gaygusuz et al (2006), or from fish pictures, e.g. Collette and Nauen (1983), Compagno (1984)
and Randall (1997). The relationships, which always refer to centimeters, may consist either of a
regression linking two length types, of the form:
 Length type (2) = a + b x Length type (1)
Length type (2) = b' x Length type (1)
The available length types are, as elsewhere in FishBase,
TL = total length;
FL = fork length;
SL = standard length;
WD = width (in rays);
OT = other type (to be specified in the Comment field).
When a version of equation (1) is presented, the length range, the number of fish used in the regression,
the sex and the correlation coefficient are presented, if available.
When a version of equation (2) is presented, the range and the correlation coefficient are omitted,
as the ratio in (2) will usually be estimated from a single specimen, or a few fish covering a narrow
range of lengths.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
length_length("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/PDF/FB_Book_CBinohlan_Length-Length_RF_JG.pdf
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parquet_db.R
\name{fb_tables}
\alias{fb_tables}
\title{fb_tables
list tables}
\usage{
fb_tables(server = c("fishbase", "sealifebase"), version = "latest")
}
\arguments{
\item{server}{fishbase or sealifebase}

\item{version}{release version}
}
\description{
fb_tables
list tables
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trophic_ecology.R
\name{diet_items}
\alias{diet_items}
\title{diet_items}
\usage{
diet_items(...)
}
\arguments{
\item{...}{additional arguments (not used)}
}
\value{
a table of diet_items
}
\description{
diet_items
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
diet_items()
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trophic_ecology.R
\name{diet}
\alias{diet}
\title{diet}
\usage{
diet(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species diet
}
\description{
diet
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
diet()
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/fishbasethe_diet_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parquet_db.R
\name{fb_import}
\alias{fb_import}
\title{Import tables to local store}
\usage{
fb_import(
  server = c("fishbase", "sealifebase"),
  version = "latest",
  db = fb_conn(server, version),
  tables = NULL
)
}
\arguments{
\item{server}{fishbase or sealifebase}

\item{version}{release version}

\item{db}{A cachable duckdb database connection}

\item{tables}{list of tables to import. Default `NULL` will
import all tables.}
}
\description{
Import tables to local store
}
\details{
Downloads and stores tables from the requested version of 
fishbase or sealifebase.  If the table is already downloaded, it will
not be re-downloaded.  Imported tables are added to the active duckdb
connection. Note that there is no need to call this
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
conn <- fb_import()
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{countrysub}
\alias{countrysub}
\title{countrysub}
\usage{
countrysub(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\description{
return a table of countrysub for the requested species
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
countrysub(species_list(Genus='Labroides'))
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_taxa.R
\docType{data}
\name{fishbase}
\alias{fishbase}
\title{A table of all the the species found in FishBase, including taxonomic
classification and the Species Code (SpecCode) by which the species is
identified in FishBase.}
\description{
A table of all the the species found in FishBase, including taxonomic
classification and the Species Code (SpecCode) by which the species is
identified in FishBase.
}
\author{
Carl Boettiger \email{carl@ropensci.org}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trophic_ecology.R
\name{ration}
\alias{ration}
\title{ration}
\usage{
ration(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species ration
}
\description{
ration
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
ration("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/fishbasethe_ration_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{introductions}
\alias{introductions}
\title{introductions}
\usage{
introductions(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species introductions data
}
\description{
introductions
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
introductions("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population_dynamics.R
\name{popgrowth}
\alias{popgrowth}
\title{popgrowth}
\usage{
popgrowth(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of population growth information by species; see details
}
\description{
This table contains information on growth, natural mortality and length
at first maturity, which serve as inputs to many fish stock assessment
models. The data can also be used to generate empirical relationships
between growth parameters or natural mortality estimates, and their
correlates (e.g., body shape, temperature, etc.), a line of research
that is useful both for stock assessment and for increasing understanding
of the evolution of life-history strategies
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
popgrowth("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/fishbasethe_popgrowth_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population_dynamics.R
\name{length_weight}
\alias{length_weight}
\alias{poplw}
\title{length_weight}
\usage{
length_weight(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of length_weight information by species; see details
}
\description{
The LENGTH-WEIGHT table presents the a and b values of over 5,000
length-weight relationships of the form W = a x Lb, pertaining to about over 2,000 fish species.
}
\details{
See references for official documentation.  From FishBase.org:
Length-weight relationships are important in fisheries science, 
notably to raise length-frequency samples to total catch, or to
estimate biomass from underwater length observations. 
The units of length and weight in FishBase are centimeter and gram, respectively. 
Thus when length-weight relationships are not in cm-g, the intercept 'a' 
is transformed as follows:

a'(cm, g) = a (mm, g)*10^b
a'(cm, g) = a (cm, kg)*1000
a'(cm, g) = a (mm, mg)*10^b/1000
a'(cm, g) = a (mm, kg)*10^b*1000

However, published length-weight relationships are sometimes difficult to use,
as they may be based on a length measurement type (e.g., fork length) different
from ones length measurements (expressed e.g., as total length).
Therefore, to facilitate conversion between length types, an additional
LENGTH-LENGTH table, #' presented below, was devised which presents linear
regressions or ratios linking length types (e.g., FL vs. TL). 
We included a calculated field with the weight of a 10 cm fish (which
should be in the order of 10 g for normal, fusiform shaped fish),
to allow identification of gross errors, given knowledge of the body
form of a species.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
length_weight("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/fishbasethe_length_weight_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population_dynamics.R
\name{length_freq}
\alias{length_freq}
\alias{poplf}
\title{length_freq}
\usage{
length_freq(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of length_freq information by species; see details
}
\description{
return a table of species fooditems
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
length_freq("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/lengthfrequency.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs.R
\name{docs}
\alias{docs}
\title{docs}
\usage{
docs(table = NULL, server = NULL, ...)
}
\arguments{
\item{table}{the table for which the documentation should be displayed.  If no table is given,
documentation summarizing all available tables is shown.}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{...}{unused; for backwards compatibility only}
}
\value{
A data.frame which lists the name of each table (if no table argument is given), along with a description
of the table and a URL linking to further information about the table.  If a specific table is named in the
table argument, then the function will return a data.frame listing all the fields (columns) found in that table, 
a description of what the field label means, and the units in which the field is measured.  These descriptions of the
columns are not made available by FishBase and must be manually generated and curated by FishBase users. 
At this time, many fields are still missing.  Please take a moment to fill in any fields you use in the source
table here: https://github.com/ropensci/fishbaseapi/tree/master/docs/docs-sources
}
\description{
documentation of tables and fields
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
tables <- docs()
# Describe the fecundity table
dplyr::filter(tables, table == "fecundity")$description
## See fields in fecundity table
docs("fecundity")
## Note: only 
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reproduction.R
\name{spawning}
\alias{spawning}
\title{spawning}
\usage{
spawning(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species spawning
}
\description{
spawning
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
spawning("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common_names.R
\name{common_names}
\alias{common_names}
\alias{sci_to_common}
\title{common names}
\usage{
common_names(
  species_list = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  Language = "English",
  fields = NULL
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{Language}{a string specifying the language for the common name, e.g. "English"}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}
}
\value{
a data.frame of common names by species queried. If multiple species are queried,
The resulting data.frames are concatenated.
}
\description{
Return a table of common names
}
\details{
Note that there are many common names for a given sci name, so sci_to_common doesn't make sense
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{faoareas}
\alias{faoareas}
\title{faoareas}
\usage{
faoareas(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a tibble, empty tibble if no results found
}
\description{
return a table of species locations as reported in FishBASE.org FAO location data
}
\details{
currently this is ~ FAO areas table (minus "note" field)
e.g. http://www.fishbase.us/Country/FaoAreaList.php?ID=5537
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
  faoareas()
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/available_releases.R
\name{available_releases}
\alias{available_releases}
\title{List available releases}
\usage{
available_releases(server = c("fishbase", "sealifebase"))
}
\arguments{
\item{server}{fishbase or sealifebase}
}
\description{
List available releases
}
\details{
Lists all available releases (year.month format).  
To use a specific release, set the desired release using
`options(FISHBASE_VERSION=)`, as shown in the examples. 
Otherwise, rfishbase will use the latest available version if this
option is unset.  NOTE: it will be necessary 
to clear the cache with `clear_cache()` or by restarting the R session
with a fresh environment.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
available_releases()
options(FISHBASE_VERSION="19.04")
## unset
options(FISHBASE_VERSION=NULL)
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morpho_physio.R
\name{morphometrics}
\alias{morphometrics}
\title{morphometrics}
\usage{
morphometrics(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species morphometrics data
}
\description{
morphometrics
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
morphometrics("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common_names.R
\name{common_to_sci}
\alias{common_to_sci}
\title{common_to_sci}
\usage{
common_to_sci(
  x,
  Language = "English",
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db()
)
}
\arguments{
\item{x}{a common name or list of common names}

\item{Language}{a string specifying the language for the common name, e.g. "English"}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}
}
\value{
a character vector of scientific names
}
\description{
Return a list of scientific names corresponding to given the common name(s).
}
\details{
If more than one scientific name matches the common name (e.g. "trout"), the function
will simply return a list of all matching scientific names.  If given more than one common name,
the resulting strings of matching scientific names are simply concatenated.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\donttest{
common_to_sci(c("Bicolor cleaner wrasse", "humphead parrotfish"), Language="English")
common_to_sci(c("Coho Salmon", "trout"))
}
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link{species_list}}, \code{\link{synonyms}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_taxa.R
\name{load_taxa}
\alias{load_taxa}
\title{load_taxa}
\usage{
load_taxa(
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  collect = TRUE,
  ...
)
}
\arguments{
\item{server}{Either "fishbase" (the default) or "sealifebase"}

\item{version}{the version of the database you want. Will default to the
latest avialable; see [available_releases()].}

\item{db}{A remote database connection. Will default to the best available
system, see [default_db()].}

\item{collect}{return a data.frame if TRUE, otherwise, a DBI connection to
the table in the database}

\item{...}{for compatibility with previous versions}
}
\value{
the taxa list
}
\description{
load_taxa
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reproduction.R
\name{fecundity}
\alias{fecundity}
\title{fecundity}
\usage{
fecundity(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species fecundity
}
\description{
fecundity
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
fecundity("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morpho_physio.R
\name{morphology}
\alias{morphology}
\title{morphology}
\usage{
morphology(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species morphology data
}
\description{
morphology
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
 \dontrun{
morphology("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fb_tbl.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synonyms.R
\name{validate_names}
\alias{validate_names}
\title{validate_names}
\usage{
validate_names(
  species_list,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a string of the validated names
}
\description{
Check for alternate versions of a scientific name and return
the scientific names FishBase recognizes as valid
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

  \donttest{
validate_names("Abramites ternetzi")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/species.R
\name{species}
\alias{species}
\alias{species_info}
\title{species}
\usage{
species(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a data.frame with rows for species and columns for the fields returned by the query (FishBase 'species' table)
}
\description{
Provide wrapper to work with species lists.
}
\details{
The Species table is the heart of FishBase. This function provides a convenient way to 
query, tidy, and assemble data from that table given an entire list of species.
For details, see: http://www.fishbase.org/manual/english/fishbasethe_species_table.htm

Species scientific names are defined according to fishbase taxonomy and nomenclature.
}
\examples{
\dontrun{

species(c("Labroides bicolor", "Bolbometopon muricatum")) 
species(c("Labroides bicolor", "Bolbometopon muricatum"), fields = species_fields$habitat) 

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{species_by_ecosystem}
\alias{species_by_ecosystem}
\title{Species list by ecosystem}
\usage{
species_by_ecosystem(
  ecosystem,
  species_list = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  ...
)
}
\arguments{
\item{ecosystem}{(character) an ecosystem name}

\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species ecosystems data
}
\description{
Species list by ecosystem
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
species_by_ecosystem(ecosystem = "Arctic", server = "sealifebase")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trophic_ecology.R
\name{estimate}
\alias{estimate}
\title{estimate}
\usage{
estimate(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of estimates from some models on trophic levels
}
\description{
estimate
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
estimate("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.us/manual/English/FishbaseThe_FOOD_ITEMS_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reproduction.R
\name{larvae}
\alias{larvae}
\title{larvae}
\usage{
larvae(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of larval data
}
\description{
larvae
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
larvae("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reproduction.R
\name{reproduction}
\alias{reproduction}
\title{reproduction}
\usage{
reproduction(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species reproduction
}
\description{
reproduction
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
reproduction("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{c_code}
\alias{c_code}
\title{c_code}
\usage{
c_code(
  c_code = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  ...
)
}
\arguments{
\item{c_code}{a C_Code or list of C_Codes (FishBase country code)}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\description{
return a table of country information for the requested c_code, as reported in FishBASE.org
}
\details{
e.g. http://www.fishbase.us/Country
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
c_code(440)
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morpho_physio.R
\name{swimming}
\alias{swimming}
\title{swimming}
\usage{
swimming(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species swimming data
}
\description{
swimming
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
swimming("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/species.R
\docType{data}
\name{species_fields}
\alias{species_fields}
\title{A list of the species_fields available}
\description{
A list of the species_fields available
}
\author{
Carl Boettiger \email{carl@ropensci.org}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_db.R
\name{fb_conn}
\alias{fb_conn}
\title{Cacheable database connection}
\usage{
fb_conn(server = c("fishbase", "sealifebase"), version = "latest")
}
\arguments{
\item{server}{fishbase or sealifebase}

\item{version}{release version}
}
\description{
Cacheable database connection
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{distribution}
\alias{distribution}
\title{distribution}
\usage{
distribution(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\description{
return a table of species locations as reported in FishBASE.org FAO location data
}
\details{
currently this is ~ FAO areas table (minus "note" field)
e.g. http://www.fishbase.us/Country/FaoAreaList.php?ID=5537
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
distribution(species_list(Genus='Labroides'))
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{stocks}
\alias{stocks}
\title{stocks}
\usage{
stocks(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species stocks data
}
\description{
stocks
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
stocks("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morpho_physio.R
\name{speed}
\alias{speed}
\title{speed}
\usage{
speed(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species speed data
}
\description{
speed
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
speed("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morpho_physio.R
\name{oxygen}
\alias{oxygen}
\title{oxygen}
\usage{
oxygen(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species oxygen data
}
\description{
oxygen
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
 \dontrun{
oxygen("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morpho_physio.R
\name{genetics}
\alias{genetics}
\title{genetics}
\usage{
genetics(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species genetics data
}
\description{
genetics
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
genetics("Oreochromis niloticus")
genetics("Labroides dimidiatus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brains.R
\name{brains}
\alias{brains}
\title{brains}
\usage{
brains(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species brains
}
\description{
brains
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
 \dontrun{
brains("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reproduction.R
\name{maturity}
\alias{maturity}
\title{maturity}
\usage{
maturity(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species maturity
}
\description{
maturity
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
maturity("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trophic_ecology.R
\name{popqb}
\alias{popqb}
\title{popqb}
\usage{
popqb(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\value{
a table of species popqb
}
\description{
popqb
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
popqb("Oreochromis niloticus")
}
\dontshow{\}) # examplesIf}
}
\references{
http://www.fishbase.org/manual/english/fishbasethe_popqb_table.htm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution.R
\name{country}
\alias{country}
\title{country}
\usage{
country(
  species_list = NULL,
  fields = NULL,
  server = getOption("FISHBASE_API", "fishbase"),
  version = get_latest_release(),
  db = default_db(server, version),
  ...
)
}
\arguments{
\item{species_list}{A vector of scientific names (each element as "genus species"). If empty, a table for all fish will be returned.}

\item{fields}{a character vector specifying which fields (columns) should be returned. By default,
all available columns recognized by the parser are returned.  Mostly for backwards compatibility as users can subset by column later}

\item{server}{can be set to either "fishbase" or "sealifebase" to switch between databases. NOTE: it is usually
easier to leave this as NULL and set the source instead using the environmental variable `FISHBASE_API`, e.g.
`Sys.setenv(FISHBASE_API="sealifebase")`.}

\item{version}{a version string for the database, will default to the latest release. see [get_releases()] for details.}

\item{db}{the}

\item{...}{unused; for backwards compatibility only}
}
\description{
return a table of country for the requested species, as reported in FishBASE.org
}
\details{
e.g. http://www.fishbase.us/Country
}
\examples{
\dontshow{if (interactive() ) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
\dontrun{
country(species_list(Genus='Labroides'))
}
\dontshow{\}) # examplesIf}
}
