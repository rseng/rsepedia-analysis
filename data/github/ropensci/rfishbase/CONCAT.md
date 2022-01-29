
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
