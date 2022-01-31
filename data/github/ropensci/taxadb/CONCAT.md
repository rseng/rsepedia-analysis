
<!-- README.md is generated from README.Rmd. Please edit that file -->

# taxadb <img src="man/figures/logo.svg" align="right" alt="" width="120" />

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/taxadb/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/taxadb/actions)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Coverage
status](https://codecov.io/gh/ropensci/taxadb/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/taxadb?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/taxadb)](https://cran.r-project.org/package=taxadb)
[![DOI](https://zenodo.org/badge/130153207.svg)](https://zenodo.org/badge/latestdoi/130153207)

<!-- [![peer-review](https://badges.ropensci.org/344_status.svg)](https://github.com/ropensci/software-review/issues/344) -->

<!-- badges: end -->

The goal of `taxadb` is to provide *fast*, *consistent* access to
taxonomic data, supporting common tasks such as resolving taxonomic
names to identifiers, looking up higher classification ranks of given
species, or returning a list of all species below a given rank. These
tasks are particularly common when synthesizing data across large
species assemblies, such as combining occurrence records with trait
records.

Existing approaches to these problems typically rely on web APIs, which
can make them impractical for work with large numbers of species or in
more complex pipelines. Queries and returned formats also differ across
the different taxonomic authorities, making tasks that query multiple
authorities particularly complex. `taxadb` creates a *local* database of
most readily available taxonomic authorities, each of which is
transformed into consistent, standard, and researcher-friendly tabular
formats.

## Install and initial setup

To get started, install from CRAN

``` r
install.pacakges("taxadb")
```

or install the development version directly from GitHub:

``` r
devtools::install_github("ropensci/taxadb")
```

``` r
library(taxadb)
library(dplyr) # Used to illustrate how a typical workflow combines nicely with `dplyr`
```

Create a local copy of the Catalogue of Life (2018) database:

``` r
td_create("col", overwrite=FALSE)
```

Read in the species list used by the Breeding Bird Survey:

``` r
bbs_species_list <- system.file("extdata/bbs.tsv", package="taxadb")
bbs <- read.delim(bbs_species_list)
```

## Getting names and ids

Two core functions are `get_ids()` and `get_names()`. These functions
take a vector of names or ids (respectively), and return a vector of ids
or names (respectively). For instance, we can use this to attempt to
resolve all the bird names in the Breeding Bird Survey against the
Catalogue of Life:

``` r
birds <- bbs %>% 
  select(species) %>% 
  mutate(id = get_ids(species, "col"))

head(birds, 10)
#>                          species          id
#> 1         Dendrocygna autumnalis COL:3177882
#> 2            Dendrocygna bicolor COL:3177881
#> 3                Anser canagicus COL:3178026
#> 4             Anser caerulescens COL:3178024
#> 5  Chen caerulescens (blue form)        <NA>
#> 6                   Anser rossii COL:3178025
#> 7                Anser albifrons COL:3178017
#> 8                Branta bernicla COL:3178037
#> 9      Branta bernicla nigricans COL:3185200
#> 10             Branta hutchinsii COL:3178039
```

Note that some names cannot be resolved to an identifier. This can occur
because of miss-spellings, non-standard formatting, or the use of a
synonym not recognized by the naming provider. Names that cannot be
uniquely resolved because they are known synonyms of multiple different
species will also return `NA`. The `filter_name` filtering functions can
help us resolve this last case (see below).

`get_ids()` returns the IDs of accepted names, that is
`dwc:AcceptedNameUsageID`s. We can resolve the IDs into accepted names:

``` r
birds %>% 
  mutate(accepted_name = get_names(id, "col")) %>% 
  head()
#>                         species          id          accepted_name
#> 1        Dendrocygna autumnalis COL:3177882 Dendrocygna autumnalis
#> 2           Dendrocygna bicolor COL:3177881    Dendrocygna bicolor
#> 3               Anser canagicus COL:3178026          Chen canagica
#> 4            Anser caerulescens COL:3178024      Chen caerulescens
#> 5 Chen caerulescens (blue form)        <NA>                   <NA>
#> 6                  Anser rossii COL:3178025            Chen rossii
```

This illustrates that some of our names, e.g. *Dendrocygna bicolor* are
accepted in the Catalogue of Life, while others, *Anser canagicus* are
**known synonyms** of a different accepted name: **Chen canagica**.
Resolving synonyms and accepted names to identifiers helps us avoid the
possible miss-matches we could have when the same species is known by
two different names.

## Taxonomic Data Tables

Local access to taxonomic data tables lets us do much more than look up
names and ids. A family of `filter_*` functions in `taxadb` help us work
directly with subsets of the taxonomic data. As we noted above, this can
be useful in resolving certain ambiguous names.

For instance, *Trochalopteron henrici gucenense* does not resolve to an
identifier in ITIS:

``` r
get_ids("Trochalopteron henrici gucenense") 
#> Warning:   Found 2 possible identifiers for Trochalopteron henrici gucenense.
#>   Returning NA. Try filter_id('Trochalopteron henrici gucenense', 'itis') to resolve manually.
#> [1] NA
```

Using `filter_name()`, we find this is because the name resolves not to
zero matches, but to more than one match:

``` r
filter_name("Trochalopteron henrici gucenense") 
#> # A tibble: 2 x 17
#>    sort taxonID   scientificName      taxonRank acceptedNameUsa… taxonomicStatus
#>   <int> <chr>     <chr>               <chr>     <chr>            <chr>          
#> 1     1 ITIS:924… Trochalopteron hen… subspeci… ITIS:916116      synonym        
#> 2     1 ITIS:924… Trochalopteron hen… subspeci… ITIS:916117      synonym        
#> # … with 11 more variables: update_date <chr>, kingdom <chr>, phylum <chr>,
#> #   class <chr>, order <chr>, family <chr>, genus <chr>, specificEpithet <chr>,
#> #   infraspecificEpithet <chr>, vernacularName <chr>, input <chr>
```

``` r
filter_name("Trochalopteron henrici gucenense")  %>%
  mutate(acceptedNameUsage = get_names(acceptedNameUsageID)) %>% 
  select(scientificName, taxonomicStatus, acceptedNameUsage, acceptedNameUsageID)
#> # A tibble: 2 x 4
#>   scientificName          taxonomicStatus acceptedNameUsage    acceptedNameUsag…
#>   <chr>                   <chr>           <chr>                <chr>            
#> 1 Trochalopteron henrici… synonym         Trochalopteron elli… ITIS:916116      
#> 2 Trochalopteron henrici… synonym         Trochalopteron henr… ITIS:916117
```

Similar functions `filter_id`, `filter_rank`, and `filter_common` take
IDs, scientific ranks, or common names, respectively. Here, we can get
taxonomic data on all bird names in the Catalogue of Life:

``` r
filter_rank(name = "Aves", rank = "class", provider = "col")
#> # A tibble: 36,336 x 21
#>     sort taxonID   scientificName    acceptedNameUsag… taxonomicStatus taxonRank
#>    <int> <chr>     <chr>             <chr>             <chr>           <chr>    
#>  1     1 COL:3148… Nisaetus nanus    COL:3148416       accepted        species  
#>  2     1 COL:3148… Circaetus beaudo… COL:3148666       accepted        species  
#>  3     1 COL:3148… Cariama cristata  COL:3148731       accepted        species  
#>  4     1 COL:3148… Chunga burmeiste… COL:3148732       accepted        species  
#>  5     1 COL:3148… Eurypyga helias   COL:3148733       accepted        species  
#>  6     1 COL:3148… Rhynochetos juba… COL:3148734       accepted        species  
#>  7     1 COL:3148… Leptosomus disco… COL:3148735       accepted        species  
#>  8     1 COL:3148… Neotis heuglinii  COL:3148736       accepted        species  
#>  9     1 COL:3148… Neotis ludwigii   COL:3148737       accepted        species  
#> 10     1 COL:3148… Neotis denhami    COL:3148738       accepted        species  
#> # … with 36,326 more rows, and 15 more variables: kingdom <chr>, phylum <chr>,
#> #   class <chr>, order <chr>, family <chr>, genus <chr>, specificEpithet <chr>,
#> #   infraspecificEpithet <chr>, taxonConceptID <chr>, isExtinct <chr>,
#> #   nameAccordingTo <chr>, namePublishedIn <chr>,
#> #   scientificNameAuthorship <chr>, vernacularName <chr>, input <chr>
```

Combining these with `dplyr` functions can make it easy to explore this
data: for instance, which families have the most species?

``` r
filter_rank(name = "Aves", rank = "class", provider = "col") %>%
  filter(taxonomicStatus == "accepted", taxonRank=="species") %>% 
  group_by(family) %>%
  count(sort = TRUE) %>% 
  head()
#> # A tibble: 6 x 2
#> # Groups:   family [6]
#>   family           n
#>   <chr>        <int>
#> 1 Tyrannidae     401
#> 2 Thraupidae     374
#> 3 Psittacidae    370
#> 4 Columbidae     344
#> 5 Trochilidae    338
#> 6 Muscicapidae   314
```

## Using the database connection directly

`filter_*` functions by default return in-memory data frames. Because
they are filtering functions, they return a subset of the full data
which matches a given query (names, ids, ranks, etc), so the returned
data.frames are smaller than the full record of a naming provider.
Working directly with the SQL connection to the MonetDBLite database
gives us access to all the data. The `taxa_tbl()` function provides this
connection:

``` r
taxa_tbl("col")
#> # Source:   table<2020_dwc_col> [?? x 19]
#> # Database: duckdb_connection
#>    taxonID  scientificName    acceptedNameUsa… taxonomicStatus taxonRank kingdom
#>    <chr>    <chr>             <chr>            <chr>           <chr>     <chr>  
#>  1 COL:3738 Lobesia triacant… COL:3738         accepted        species   Animal…
#>  2 COL:4116 Syncollesis tril… COL:4116         accepted        species   Animal…
#>  3 COL:4118 Anisodes anablem… COL:4118         accepted        species   Animal…
#>  4 COL:4122 Cyclophora carol… COL:4122         accepted        species   Animal…
#>  5 COL:4127 Morchella magnis… COL:4127         accepted        species   Fungi  
#>  6 COL:4128 Streptothrix eff… COL:4128         accepted        species   Fungi  
#>  7 COL:4344 Aplosporella fau… COL:4344         accepted        species   Fungi  
#>  8 COL:9466 Synalus angustus  COL:9466         accepted        species   Animal…
#>  9 COL:9467 Synalus terrosus  COL:9467         accepted        species   Animal…
#> 10 COL:9468 Synema spinosum   COL:9468         accepted        species   Animal…
#> # … with more rows, and 13 more variables: phylum <chr>, class <chr>,
#> #   order <chr>, family <chr>, genus <chr>, specificEpithet <chr>,
#> #   infraspecificEpithet <chr>, taxonConceptID <chr>, isExtinct <chr>,
#> #   nameAccordingTo <chr>, namePublishedIn <chr>,
#> #   scientificNameAuthorship <chr>, vernacularName <chr>
```

We can still use most familiar `dplyr` verbs to perform common tasks.
For instance: which species has the most known synonyms?

``` r
taxa_tbl("col") %>% 
  count(acceptedNameUsageID, sort=TRUE)
#> # Source:     lazy query [?? x 2]
#> # Database:   duckdb_connection
#> # Ordered by: desc(n)
#>    acceptedNameUsageID     n
#>    <chr>               <dbl>
#>  1 COL:274062            456
#>  2 COL:353741            373
#>  3 COL:3778950           329
#>  4 COL:2535424           328
#>  5 COL:2921616           322
#>  6 COL:2532677           307
#>  7 COL:3779182           302
#>  8 COL:353740            296
#>  9 COL:2531203           253
#> 10 COL:1585420           252
#> # … with more rows
```

However, unlike the `filter_*` functions which return convenient
in-memory tables, this is still a remote connection. This means that
direct access using the `taxa_tbl()` function (or directly accessing the
database connection using `td_connect()`) is more low-level and requires
greater care. For instance, we cannot just add a `%>%
mutate(acceptedNameUsage = get_names(acceptedNameUsageID))` to the
above, because `get_names` does not work on a remote collection.
Instead, we would first need to use a `collect()` to pull the summary
table into memory. Users familiar with remote databases in `dplyr` will
find using `taxa_tbl()` directly to be convenient and fast, while other
users may find the `filter_*` approach to be more intuitive.

## Learn more

  - See richer examples the package
    [Tutorial](https://docs.ropensci.org/taxadb/articles/articles/intro.html).

  - Learn about the underlying data sources and formats in [Data
    Sources](https://docs.ropensci.org/taxadb/articles/data-sources.html)

  - Get better performance by selecting an alternative [database
    backend](https://docs.ropensci.org/taxadb/articles/backends.html)
    engines.

-----

Please note that this project is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By participating in
this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# taxadb 0.1.4

* bugfix in `get_ids()` when multiple English common names are accepted for the species.
* export `taxadb_dir()`, making it easier to purge the DB after `duckdb` upgrades
* All imports must be used
* Improve testing in `db=NULL` case.
* Require R.utils, to ensure compressed files can be expanded

# taxadb 0.1.3

* more robust testing

# taxadb 0.1.2

* avoid erroneous messages when installing providers that lack common names.

# taxadb 0.1.1

* introduce `tl_import` to import taxonomic databases [#79]
* make `duckdb` the default backend
* bugfix to possible ordering problem in `get_names` [#78]
* Added a `NEWS.md` file to track changes to the package.
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
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
Dear CRAN Maintainers,

Changes in this release are described in NEWS.md

Cheers,

Carl
# data-raw

The scripts used to update the cache are now being maintained separately at
<https://github.com/boettiger-lab/taxadb-cache>

Please consult that repository for the latest scripts, functions, and metadata regarding the backend cache for `taxadb`.

---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = TRUE,
  comment = "#>"
)
taxadb:::td_disconnect()
```

# taxadb  <img src="man/figures/logo.svg" align="right" alt="" width="120" />

<!-- badges: start -->
[![R build status](https://github.com/ropensci/taxadb/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/taxadb/actions)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Coverage status](https://codecov.io/gh/ropensci/taxadb/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/taxadb?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/taxadb)](https://cran.r-project.org/package=taxadb)
[![DOI](https://zenodo.org/badge/130153207.svg)](https://zenodo.org/badge/latestdoi/130153207)

<!-- [![peer-review](https://badges.ropensci.org/344_status.svg)](https://github.com/ropensci/software-review/issues/344) -->
<!-- badges: end -->


The goal of `taxadb` is to provide *fast*, *consistent* access to taxonomic data, supporting common tasks such as resolving taxonomic names to identifiers, looking up higher classification ranks of given species, or returning a list of all species below a given rank. These tasks are particularly common when synthesizing data across large species assemblies, such as combining occurrence records with trait records. 

Existing approaches to these problems typically rely on web APIs, which can make them impractical for work with large numbers of species or in more complex pipelines.  Queries and returned formats also differ across the different taxonomic authorities, making tasks that query multiple authorities particularly complex. `taxadb` creates a *local* database of most readily available taxonomic authorities, each of which is transformed into consistent, standard, and researcher-friendly tabular formats.  


## Install and initial setup

To get started, install from CRAN

``` r
install.pacakges("taxadb")
```

or install the development version directly from GitHub:

``` r
devtools::install_github("ropensci/taxadb")
```


```{r message = FALSE}
library(taxadb)
library(dplyr) # Used to illustrate how a typical workflow combines nicely with `dplyr`
```

Create a local copy of the Catalogue of Life (2018) database: 

```{r }
td_create("col", overwrite=FALSE)
```


Read in the species list used by the Breeding Bird Survey:

```{r, message = FALSE}
bbs_species_list <- system.file("extdata/bbs.tsv", package="taxadb")
bbs <- read.delim(bbs_species_list)
```

## Getting names and ids

Two core functions are `get_ids()` and `get_names()`.  These functions take a vector of names or ids (respectively), and return a vector of ids or names (respectively).  For instance, we can use this to attempt to resolve all the bird names in the Breeding Bird Survey against the Catalogue of Life:


```{r}
birds <- bbs %>% 
  select(species) %>% 
  mutate(id = get_ids(species, "col"))

head(birds, 10)
```

Note that some names cannot be resolved to an identifier.  This can occur because of miss-spellings, non-standard formatting, or the use of a synonym not recognized by the naming provider.  Names that cannot be uniquely resolved because they are known synonyms of multiple different species will also return `NA`.  The `filter_name` filtering functions can help us resolve this last case (see below).

`get_ids()` returns the IDs of accepted names, that is `dwc:AcceptedNameUsageID`s.  We can resolve the IDs into accepted names:


```{r}
birds %>% 
  mutate(accepted_name = get_names(id, "col")) %>% 
  head()
```

This illustrates that some of our names, e.g. *Dendrocygna bicolor* are accepted in the Catalogue of Life, while others, *Anser canagicus* are **known synonyms** of a different accepted name: **Chen canagica**.  Resolving synonyms and accepted names to identifiers helps us avoid the possible miss-matches we could have when the same species is known by two different names.


## Taxonomic Data Tables

Local access to taxonomic data tables lets us do much more than look up names and ids.  A family of `filter_*` functions in `taxadb` help us work directly with subsets of the taxonomic data.  As we noted above, this can be useful in resolving certain ambiguous names.  

For instance, *Trochalopteron henrici gucenense* does not resolve to an identifier in ITIS:

```{r}
get_ids("Trochalopteron henrici gucenense") 
```

Using `filter_name()`, we find this is because the name resolves not to zero matches, but to more than one match:

```{r}
filter_name("Trochalopteron henrici gucenense") 
```


```{r}
filter_name("Trochalopteron henrici gucenense")  %>%
  mutate(acceptedNameUsage = get_names(acceptedNameUsageID)) %>% 
  select(scientificName, taxonomicStatus, acceptedNameUsage, acceptedNameUsageID)
```


Similar functions `filter_id`, `filter_rank`, and `filter_common` take IDs, scientific ranks, or common names, respectively.  Here, we can get taxonomic data on all bird names in the Catalogue of Life:


```{r}
filter_rank(name = "Aves", rank = "class", provider = "col")
```

Combining these with `dplyr` functions can make it easy to explore this data: for instance, which families have the most species?


```{r}
filter_rank(name = "Aves", rank = "class", provider = "col") %>%
  filter(taxonomicStatus == "accepted", taxonRank=="species") %>% 
  group_by(family) %>%
  count(sort = TRUE) %>% 
  head()
```

## Using the database connection directly

`filter_*` functions by default return in-memory data frames.  Because they are filtering functions, they return a subset of the full data which matches a given query (names, ids, ranks, etc), so the returned data.frames are smaller than the full record of a naming provider.  Working directly with the SQL connection to the MonetDBLite database gives us access to all the data. The `taxa_tbl()` function provides this connection:

```{r}
taxa_tbl("col")
```

We can still use most familiar `dplyr` verbs to perform common tasks.  For instance: which species has the most known synonyms?

```{r}
taxa_tbl("col") %>% 
  count(acceptedNameUsageID, sort=TRUE)
```

However, unlike the `filter_*` functions which return convenient in-memory tables, this is still a remote connection.  This means that direct access using the `taxa_tbl()` function (or directly accessing the database connection using `td_connect()`) is more low-level and requires greater care.  For instance, we cannot just add a `%>% mutate(acceptedNameUsage = get_names(acceptedNameUsageID))` to the above, because `get_names` does not work on a remote collection.  Instead, we would first need to use a `collect()` to pull the summary table into memory.  Users familiar with remote databases in `dplyr` will find using `taxa_tbl()` directly to be convenient and fast, while other users may find the `filter_*` approach to be more intuitive.



## Learn more

- See richer examples the package [Tutorial](https://docs.ropensci.org/taxadb/articles/articles/intro.html).

- Learn about the underlying data sources and formats in [Data Sources](https://docs.ropensci.org/taxadb/articles/data-sources.html)

- Get better performance by selecting an alternative [database backend](https://docs.ropensci.org/taxadb/articles/backends.html) engines.



```{r include=FALSE}
taxadb:::td_disconnect()

if(require(codemetar)) codemetar::write_codemeta()
```

----

Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/).
By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "taxadb: A High-Performance Local Taxonomic Database Interface"
author:
  - name: "Kari E A Norman"
    affiliation: ucb
  - name: "Scott Chamberlain"
    affiliation: ropensci
  - name: "Carl Boettiger"
    affiliation: ucb, ropensci
address:
  - code: ucb
    address: "Dept of Environmental Science, Policy, and Management, University of California Berkeley, Berkeley CA 94720-3114, USA"
  - code: ropensci
    address: "The rOpenSci Project, University of California Berkeley, Berkeley CA 94720-3114, USA" 
abstract: |
  1)  A familiar and growing challenge in ecological and evolutionary research is that of establishing consistent taxonomy when combining data from separate sources. While this problem is already well understood and numerous naming authorities have been created to address the issue, most researchers lack a fast, consistent, and intuitive way to retrieve taxonomic names. 
  2) We present `taxadb` R package which creates a local database, managed automatically from within R, to provide fast operations on millions of taxonomic names. 
  3) `taxadb` provides access to established naming authorities to resolve synonyms, taxonomic identifiers, and hierarchical classification in a consistent and intuitive data format.
  4) `taxadb` makes operation on millions of taxonomic names fast and manageable.
  
journal: "Methods in Ecology & Evolution"
date: "`r Sys.Date()`"
bibliography: refs.bib
layout: 3p
header-includes:
   - \usepackage{lineno}
   - \linenumbers
output: 
  rticles::elsevier_article:
    includes:
      in_header: preamble.tex
---

```{r message=FALSE, include = FALSE}
library(kableExtra)
library(magrittr)
library(tidyverse)
library(taxadb)
library(printr)
library(rticles)
library(knitr)

## see https://blog.rstudio.com/2017/06/13/dplyr-0-7-0/, 
## https://community.rstudio.com/t/why-does-na-match-na-when-joining-two-dataframes/28785/3
pkgconfig::set_config("dplyr::na_matches" = "never")

## Display NAs as "-" because "NA" looks like a category
options(knitr.kable.NA = '-') 

printtable <- function(df, ...){
  df %>%
    kableExtra::kable("latex", booktabs=T, ...) %>%
    kableExtra::kable_styling(full_width = F, latex_options = "hold_position")
}
knitr::opts_chunk$set(cache=FALSE, message=FALSE, warning=FALSE)

# no sci notation integers pleases
options(scipen=999)

taxadb::td_disconnect()
```



As ecologists and evolutionary biologists synthesize datasets across larger and larger assemblies of species,
we face a continual challenge of maintaining consistent taxonomy. How many species are in the combined data?
Do the studies use the same names for the same species, or do they use different synonyms for the same species?
Failing to correct for such differences can lead to significant inflation of species counts and miss-aligned datasets.
These challenges have become particularly acute as it becomes increasingly common for researchers to work
across a larger number and diversity of species in any given analysis, 
which may preclude the resources or substantive taxonomic expertise for all clades 
needed to resolve scientific names [@Patterson2010].  

While these issues have long been recognized in the literature [@boyle2013; @dayrat2005; @bortolus2008; @maldonado2015; @remsen2016],
and a growing number of databases and tools have emerged over the past few decades
[e.g. @itis; @ncbi; @col; @rees2014; @alvarez2018; @wagner2016; @foster2018; @gries2014],
it remains difficult to resolve taxonomic names to a common authority in a transparent, efficient, and automatable manner.
Here, we present an R package, `taxadb`, which seeks to address this gap.

Databases of taxonomic names such as the Integrated Taxonomic Information System [ITIS; @itis], 
the National Center for Biological Information's (NCBI) Taxonomy database [@ncbi], 
the Catalogue of Life [COL; @col], and over one hundred other providers have sought to address these problems
by providing expert-curated lists of accepted taxonomic names, synonyms, associated taxonomic rank,
hierarchical classifications, and scientific authority (e.g. author and date) establishing a scientific name.
The R language [@R] is widely used in ecology and evolution [@Lai2019] and the `taxize` package [@Chamberlain2013]
has become a popular way for R users to interact with naming providers and name resolution services. 
`taxize` implements bindings to the web APIs (Application Programming Interface)
hosted by many popular taxonomic name providers. 
Nevertheless, this means that functions in the `taxize` are impacted by several major drawbacks
that are inherent in the implementation of these central API servers, such as:

- Queries require internet access at all times.
- Queries are slow and inefficient to implement and perform; frequently requiring separate API calls for each taxonomic name.
- The type of query is highly limited by the API design. For instance, it is usually impossible to make queries across the entire corpus of names, such as "which accepted name has the most known synonyms?".
- Both query formats and responses differ substantially across different naming providers, making it difficult to apply a script designed for one provider to different provider.
- Most queries are not reproducible, as the results depend on the state of the central server (and potentially the quality of the internet connection)[@rees2017].  Many names providers update the server data either continuously or at regular intervals, including both revising existing names (for spelling or changes in accepted name designation) and adding new names.

Instead of binding existing web APIs, `taxadb` is built around a set of compressed text files which are 
automatically downloaded, imported, and stored on a local database by `taxadb`. 
The largest of the taxonomic naming providers today contain under 6 million name records with uncompressed file sizes
over a GB, which can be compressed to around 50 MB and downloaded in under a minute on a 1 MB/s connection. 
By using a local database as the backend, `taxadb` allows R users to interact with large data files without large memory (RAM) requirements.  A query for a single name over the web API requires a remote server to respond, execute the query, 
and serialize the response, which can take several seconds. Thus it does not take many taxa before transferring the 
entire data set to query locally is more efficient.  Moreover, this local copy can be cached on the user's machine, 
requiring only the one-time setup, and enabling offline use and reproducible queries.  Rather than returning data
in whatever format is given by the provider, `taxadb` provides a data structure following a consistent,
standardized layout or schema following Darwin Core, which provides standard terms for biodiversity data [@Wieczorek2012].  Table 1 summarizes the list of all naming providers currently accessed by `taxadb`.  More details are provided in the Data Sources Vignette, <https://docs.ropensci.org/taxadb/articles/data-sources.html>.



```{r echo = FALSE, results = "asis"}
#provider descriptions
desc <- c(itis = "originally formed to standardize taxonomic name usage across many agencies in the United States federal government",
          ncbi = "nomenclature for sequences in the International Nucleotide Sequence Database Collaboration database",
          col = "comprehensive taxonomic effort, includes some other providers (e.g. itis)",
          gbif = "taxonomic backbone of the GBIF database, assembled from other sources including COL",
          fb = "nomenclature for global database of fishes",
          ott = "comprehensive tree of life based on phylogenetic trees and taxonomic data",
          iucn = "taxonomy for classification of species status")

tibble(provider = c("Integrated Taxonomic Information System (ITIS 2019)",
                                 "National Center for Biological Information's Taxonomy database (Biotechnology Information 2019)",
                                 "Catalogue of Life (Roskov Y. 2018)",
                                 "Global Biodiversity Information Facility Taxonomic Backbone (GBIF 2019)",
                                 "FishBase (Froese and Pauly 2019)",
                                 "Open Tree Taxonomy (J. A. Rees and Cranston 2017)",
                                 "International Union for Conservation of Nature and Natural Resources (IUCN 2019)"),
                    abbreviation = c("itis", "ncbi", "col", "gbif", "fb", "ott", "iucn"),
                    total_identifiers = map(abbreviation,
                                            ~taxa_tbl(.x) %>%
                                              select(acceptedNameUsageID) %>% pull() %>% n_distinct()),
                    description = desc
) %>%
  knitr::kable(format = "latex", col.names = c("Provider", "Abbreviation", "Number of \nIdentifiers", "Description"), escape = FALSE, 
               caption = "Descriptions of the providers supported by taxadb with their reference abbreviation and the total number of identifiers contained by each provider.") %>%
  kable_styling(latex_options= c("scale_down", "hold_position")) %>%
  kableExtra::column_spec(1, width = "5cm") %>%
  kableExtra::column_spec(4, width = "10cm")


```

# Package Overview

```{r message = FALSE, warning=FALSE}
library(tidyverse)
library(taxadb)
```

After loading our package and the tidyverse package for ease in manipulating function output, we look up the taxonomic identifier for Atlantic Cod, *Gadus morhua*, and the compliment:


```{r}
get_ids("Gadus morhua")
get_names("ITIS:164712")
```


Our first call to any `taxadb` functions will automatically set up a local, persistent database if one has not yet been created. This one-time setup will download, extract, and import the compressed data into persistent database storage (using the appropriate location specified by the operating system [see @rappdirs], or configured using the environmental variable `TAXADB_HOME`).  The example above searches for names in ITIS, the default provider, which can be configured using the `provider` argument. Any future function calls to this function or any other function using data from the same provider will be able to access this data rapidly without the need for processing or an internet connection.  

Users can also explicitly trigger this one-time setup using `td_create()` and specifying the provider abbreviation (see Table 1), or simply using `all` to install all available providers:


```{r eval=FALSE, message=FALSE}
td_create("all")
```


`taxadb` functions like `get_ids()` and `td_create()` take an optional argument, `db`, to an external database connection.  `taxadb` will work with most DBI-compliant databases such as MySQL or Postgres, but will be much faster when using a column-oriented database engine such as `duckdb` or `MonetDBLite`.  These latter options are also much easier for most users, since each can be installed directly as an R package. `taxadb` will default to the fastest available option.  `taxadb` can also run without a database backend by setting `db=NULL`, though some functions will require a lot (2-20 GB) of free RAM for this to work with many of the larger providers.  

`taxadb` uses the widely known SQLite database by default, but users are encouraged to install the optional, suggested database backends by passing the option `dependencies = TRUE` to the install command.  This installs a MonetDBLite database instance [@monetdblite], a columnar-oriented relational database requiring no additional installation while also providing persistent disk-based storage.  This also installs `duckdb`, another local columnar database which is rapidly emerging as an alternative to MonetDB and SQLite. `taxadb` will automatically detect and use these database engines if available, and automatically handles opening, caching, and closing the database connection. For large queries, MonetDBLite or `duckdb` deliver impressive improvements.  Our benchmark on resolving the 750 species names in the Breeding Bird Survey against over 3 million names known in the 2019 Catalogue of Life takes 8 minutes in SQLite but less than a second in MonetDBLite.  
  
Functions in `taxadb` are organized into several families: 

- queries that return vectors: `get_ids()` and it's complement, `get_names()`,
- queries that filter the underlying taxonomic data frames: `filter_name()`, `filter_rank()`, `filter_id()`, and `filter_common()`,
- database functions  `td_create()`, `td_connect()` and `taxa_tbl()`,
- and helper utilities, such as `clean_names()`. 


## Taxonomic Identifiers

Taxonomic identifiers provide a fundamental abstraction which lies at the heart of managing taxonomic names. For instance, by resolving scientific names to identifiers, we can identify which names are synonyms -- different scientific names used to describe the same species -- and which names are not recognized. Each naming authority provides its own identifiers for the names it recognizes. For example, the name Homo sapiens has the identifier 9606 in NCBI and 180092 in ITIS.  To avoid possible confusion, taxadb always prefixes the naming provider, e.g. NCBI:9606. Some taxonomic naming providers include separate identifiers for synonyms, see Box 1. Unmatched names may indicate an error in data entry or otherwise warrant further investigation. Taxon identifiers are also easily resolved to the original authority (scientific publication) establishing the name. The common practice of appending an author and year to a scientific name, e.g. *Poa annua annua* (Smith 1912), serves a valuable role in disambiguating different uses of the same name but can be notoriously harder to resolve to the appropriate reference, while variation in this convention creates many distinct versions of the same name [@Patterson2010].  

These issues are best illustrated using a minimal example.  We'll consider the task of combining data on bird extinction risk as assessed by the IUCN [@iucn] with data on average adult biomass, as estimated in the Elton Traits v1.0 database [@elton-traits].  To keep the example concise enough for for visual presentation we will focus on a subset involving just 10 species (Table 2, 3).

```{r message=FALSE}
trait_data <- read_tsv(system.file("extdata", "trait_data.tsv", package="taxadb"))
status_data <- read_tsv(system.file("extdata", "status_data.tsv", package="taxadb"))
```

```{r iucn_table, echo=FALSE, cache = FALSE}
status_data %>% printtable(caption = "The subset of the IUCN status data used for subsequent taxonomic identifier examples.") %>%
  column_spec(1, italic = TRUE)
```


```{r trait_table, echo = FALSE, cache = FALSE }
trait_data %>% printtable(caption = "The subset of the Elton trait data used for subsequent taxonomic identifier examples.") %>%
  column_spec(1, italic = TRUE)
```


If we attempted to join these data directly on the species names provided by each table, we would find very little overlap, with only one species name having both a body mass and an IUCN threat status resolved (Table 4). 

```{r}
joined <- full_join(trait_data, status_data, by = c("elton_name" = "iucn_name")) 
```


```{r echo = FALSE, cache = FALSE}
joined %>%
  printtable(caption = "Example IUCN and trait data joined directly on scientific name showing only one match. While common, joining on scientific name does not account for nomenclatural and taxonomic inconsistencies between databases and therefore results in seemingly very little overlap in species representation between the two.") %>%
  column_spec(1, italic = TRUE)
```

If we first resolve names used in each data set into shared identifiers, (for instance, using the Catalogue of Life), we discover that there is far more overlap in the species coverage than we might have initially realized. First, we just add an ID column to each table by looking up the Catalog of Life identifier for the name provided:

```{r}

traits <- trait_data %>% mutate(id = get_ids(elton_name, "col"))
status <- status_data %>% mutate(id = get_ids(iucn_name, "col"))
```

We can now join on the `id` column instead of names directly:

```{r}
joined <- full_join(traits, status, by = "id") 
```

```{r cache = FALSE, echo = FALSE}
## Just for pretty-printing
joined %>%  
  tidyr::replace_na(list(category = "-", elton_name = "-", iucn_name = "-")) %>%
  select(elton_name, iucn_name, mass, category, id) %>%
  printtable(caption = "Example IUCN and trait data joined on taxonomic ID. Multiple species have a different scientific name in the Elton and IUCN Redlist databases but can be match based on their COL taxonomic ID.") %>%
  column_spec(c(1,2), italic = TRUE)
```


This results in many more matches (Table 5), as different scientific names are recognized by the naming provider (Catalog of Life 2018 in this case), as *synonyms* for the same species, and thus resolve to the same taxonomic identifier.  While we have focused on a small example for visual clarity here, the `get_ids()` function in `taxadb` can quickly resolve hundreds of thousands of species names to unique identifiers, thanks to the performance of fast joins in a local MonetDBLite database.



\dummy{\Begin{tcolorbox}[title= Box 1: Taxonomic Identifiers and Synonyms, lower separated=false]}

`get_ids()` returns the `acceptedNameUsageID`, the identifier associated with the *accepted* name.  Some naming providers, such as ITIS and NCBI, provide taxonomic identifiers to both synonyms and accepted names.  Other providers, such as COL and GBIF, only provide identifiers for accepted names.  Common practice in Darwin Core archives is to provide an `acceptedNameUsageID` only for names which are synonyms, and otherwise to provide a `taxonID`.  For accepted names, the `acceptedNameUsageID` is then given as missing (`NA`), while for synonyms, the `taxonID` may be missing (`NA`).  In contrast, `taxadb` lists the `acceptedNameUsageID` for accepted names (where it matches the `taxonID`), as well as known synonyms.  This is semantically identical, but also more convenient for database interfaces, since it allows a name to mapped to its accepted identifier (or an identifier to map to it's accepted name usage) without the additional logic.  For consistency, we will use the term "identifier" to mean the `acceptedNameUsageID` rather than the more ambiguous `taxonID` (which is undefined for synonyms listed by many providers), unless explicitly stated otherwise.

\End{tcolorbox}

## Unresolved names

`get_ids` offers a first pass at matching scientific names to id, but names may remain unresolved for a number of reasons. First, a name may match to multiple accepted names, as in the case of a species that has been split. By design, these cases are left to be resolved by the researcher using the `filter_` functions to filter underlying taxonomic tables for additional information. A name may also be unresolved due to typos or improper formatting. `clean_names` addresses common formatting issues such as the inclusion of missing species epithets (e.g. `Accipiter sp.`) that prevent matches to the Genus, or intraspecific epithets such as `Colaptes auratus cafer` that prevent matches to the binomial name. These modifications are not appropriate in all settings and should be used with care. Spell check of input names is outside the scope of `taxadb`, however existing tools such as those developed by the Global Names Architecture (http://globalnames.org/apps/) could be incorporated into a `taxadb` workflow.

Names may also have an ambiguous resolution wherein a name may be resolved by a different provider than the one specified, either as an accepted name or a synonym. Mapping between providers represent a meaningful scientific statement requiring an understanding of the underlying taxonomic concepts of each provider [@franz2009; @Franz2018; @lepage2014]. The spirit of taxadb is not to automate steps that require expert knowledge, but provide access to multiple potential "taxonomic theories".


## `filter_` functions for access to underlying tables

Underlying data tables can be accessed through the family of `filter_` functions, which filter by certain attributes such as scientific name, id, common name, and rank. These functions allow us to ask general questions such as, how many bird species are there? 

```{r}
filter_rank("Aves", rank="class", provider = "col") %>%
  filter(taxonomicStatus == "accepted", taxonRank == "species") %>%
  pull(taxonID) %>%
  n_distinct()
```

We can also use this to gain a detailed look at specific species or ids.
For example, we can explore why `get_ids` fails to resolve a seemingly common species:

```{r}
multi_match <- filter_name("Abies menziesii", provider = "col")
```

```{r echo = FALSE, cache = FALSE}
multi_match %>%
  select(1:5, genus, specificEpithet) %>%
  unite(acceptedScientificName, genus, specificEpithet, sep = " ") %>%
  printtable(caption = "Some names may not resolve to an identifier using get\\_ids() because they match to more than one accepted ID. In such cases filter\\_ functions give further detail, as in the example of *Abies menziesii* below which has two accepted ID matches.", escape = FALSE) %>%
  column_spec(3, italic = TRUE)
```

We see that *Abies menziesii* is a synonym for three accepted names which the user will have to choose between (Table 6).
This is an example of how `taxadb` seeks to provide users with information from existing authorities and names providers,
rather than make a potentially arbitrary decision.  Because they return `data.frame`s,  `filter_` functions provide both potential matches.  Note that the simpler `get_` functions (`get_ids()`) consider multiple name matches as `NA` for the `id`, making them suitable for automated pipelines where manual resolution of duplicates is not an option. 

## Direct database access

The full taxonomic record in the database can also be directly accessed by `taxa_tbl()`, allowing for whole-database queries that are not possible through the API or web interface of many providers.  For example, we can easily check the coverage of accepted species names in each of the classes of vertebrates within the Catalogue of Life (Table 7):

```{r cache = FALSE}
verts <- taxa_tbl("col") %>%
  filter(taxonomicStatus == "accepted", phylum == "Chordata", taxonRank == "species") %>% 
  count(class, sort = TRUE)
```  


```{r echo = FALSE}
 verts %>%
  printtable(caption = "taxadb also provides direct access to the database, allowing dplyr or SQL queries which can compute across the entire dataset, such as counting accepted species in all vertebrate classes shown here.  This kind of query is effectively impossible in most REST API-based interfaces.")
```




\dummy{\Begin{tcolorbox}[title= Box 2: Common Names, lower separated=false]}

`taxadb` can also resolve common names to their identifier by mapping common name to the accepted scientific name. Common names have many of the same issues as scientific names but even more frequent (e.g. matching to more than one accepted name, non-standardized formatting). Common names are accessed via `filter_common` which takes a vector of common names. The user can then resolve discrepancies. 

\End{tcolorbox}


# Discussion


Some taxonomic name providers (e.g. OTT, COL, NCBI) offer periodic releases of a static names list, while many other providers (e.g. ITIS, FB, IUCN) offer name data on a rolling basis (i.e. the data returned by a given download URL is updated continuously or at arbitrary intervals without any additional indication if and how that data has changed.)  `taxadb`'s `td_create()` function downloads and stores cached snapshots from each provider, which follow an annual release model to support reproducible analyses.   All taxadb functions that download or access data include an optional argument `version` to indicate which version of the provider data should be used.  By default, `taxadb` will determine the latest version available (at the time of writing this is version `2019`). Appropriate metadata is stored with each snapshot, including scripts used to access and reformat the data files, as described in the "Data Sources" vignette, <https://docs.ropensci.org/taxadb/articles/data-sources.html>. 

Taxonomic identifiers are an essential first step for maintaining taxonomic consistency, a key task for a wide variety of applications. Despite multiple taxonomic standardization efforts, resolving names to taxonomic identifiers is often not a standard step in the research work flow due to difficulty in accessing providers and the time consuming API queries necessary for resolving even moderately sized data sets. `taxadb` fills an important gap between existing tools and typical research patterns by providing a fast, reproducible approach for matching names to taxonomic identifiers. It could also be used to verify that conclusions were robust to the choice of naming provider. `taxadb` is not intended as an improvement or replacement for any existing approaches to taxonomic name resolution. In particular, `taxadb` is not a replacement for the APIs or databases provided, but merely an interface to taxonomic naming information contained within that data. 

Lastly, we note that local database design used in `taxadb` is not unique to taxonomic names.  Despite the rapid expansion of REST API-based interfaces to ecological data [@ropensci], in our experience, much of the data relevant to ecologists and evolutionary biologists today would be also be amenable to the local database design.  The local database approach is much easier for data providers (who can leverage static scientific database repositories instead of maintaining REST servers) and often much faster for data consumers. 

# Acknowledgments
We thank the many researchers who contributed to the data and infrastructure of the various taxonomic providers we access through our package. Support for the development of this package was provided by United States Department of Energy through the Computational Sciences Graduate Fellowship (DOE CSGF) under grant number DE-FG02-97ER25308 awarded to K.E.A.N..

# Data Availability
Code for the R package can be found on GitHub at <https://github.com/ropensci/taxadb> and is archived on Zenodo at DOI:10.5281/zenodo.3903858 [@taxadb]. The taxonomic database is also stored on Github at <https://github.com/boettiger-lab/taxadb-cache>. The original taxonomic data are stored by the individual provider, see "Catalogue of Life", <http://www.catalogueoflife.org/> [@col], "ITIS", <https://www.itis.gov> [@itis], "NCBI", <https://www.ncbi.nlm.nih.gov/taxonomy> [@ncbi], "GBIF", <https://gbif.org> [@gbif], "Fishbase", <https://fishbase.se> [@fishbase], "Open Tree Taxonomy", <https://tree.opentreeoflife.org> [@Rees2017], "IUCN", <https://www.iucnredlist.org/resources/tax-sources> [@iucn]. 

# Authors' Contributions
K.E.A.N., S.C., and C.B. contributed to conceptual development of the package. K.E.A.N. and C.B. developed the package and contributed to the manuscript. 

\pagebreak


```{r include = FALSE}
td_disconnect()
```

# References
---
author: Kari Norman
return-address: 
  - Dept of Environmental Science Policy and Management 
  - \linebreak University of California Berkeley
  - \linebreak 130 Hilgard Way
  - Berkeley, CA 94709

opening: Dear Editors,
closing: Thank you for your consideration,
signature: Kari Norman, Scott Chamberlain, Carl Boettiger
return-email: kari.norman@berkeley.edu

output: komaletter::komaletter
---

We are submitting our manuscript titled “taxadb: A High-Performance Local Taxonomic Database Interface” to be considered for publication as an Application in Methods in Ecology and Evolution. 

Tracking and reconciling changes in taxonomy across time and between datasets is a core but often overlooked challenge in Ecology and Evolutionary biology. As combining larger and larger datasets becomes increasingly common, the potential for serious species mismatches has important implications for the scientific conclusions drawn. While multiple efforts to standardize taxonomic data exist and and some tools to access those providers have been developed, it remains difficult to resolve names to a taxonomic authority in quick, reproducible way. Here we present a new R package for accessing  taxonomic data using infrastructure that enables easy incorporation of taxonomic data to existing workflows and further exploration of the taxonomic providers. With our software it is possible to resolve hundreds of names in under a second making taxonomic resolution realistic for large datasets for the first time. This package was reviewed and accepted by rOpenSci and can be found at www.github.com/ropensci/taxadb. 

```{r}
library(tidyverse)
library(taxadb)
```

Populate the `taxadb` database (only required right after installing `taxadb`, otherwise this will just note the previously installed tables).

```{r}
td_create("all")
```


```{r}
algae <- read_csv("~/projects/algae-names/algae_uncleanednames_NAdropped.csv")
```


```{r}
## Input table with clean names
algae <- algae %>% mutate(input = clean_names(species), sort = 1:length(input))

## Let's get some matches
taxa <- taxa_tbl("ott") %>% 
        mutate_db(clean_names, "scientificName", "input") %>%
        right_join(algae, copy=TRUE, by="input") %>% 
        arrange(sort)  %>% 
        collect()

## lots of duplicate matches, pick the first one for now:
matched <- taxa %>% select(acceptedNameUsageID, sort) %>% distinct() %>% 
  group_by(sort) %>% top_n(1, acceptedNameUsageID)


# 46,045 / 57,700 have been matched!

## Who is unmatched? (sort id is in algae table but not in the matched table)
unmatched <- anti_join(algae, matched, by="sort")

# 11,655 are still unmatched.  Many appear to be known synonyms to Algaebase...
unmatched %>% count(source) %>% arrange(desc(n))

```




-------




I'm using names given for Red, Green, and Brown Algae from [Guiry (2012)](https://doi.org/10.1111/j.1529-8817.2012.01222.x):

Open Tree Taxonomy (OTT) has 32,347 names recognized as belonging to one the three phyla:

```{r}
ott_phyla <- bind_rows(
  descendants(name = "Cyanobacteria", rank = "phylum", authority = "ott"),
  descendants(name = "Rhodophyta", rank = "phylum", authority = "ott"),
  descendants(name = "Phaeophyceae", rank = "phylum", authority = "ott")
)

ott_phyla
```


GBIF has 7,412 recognized names 

```{r}
gbif_phyla <- bind_rows(
  descendants(name = "Cyanobacteria", rank = "phylum", authority = "gbif"),
  descendants(name = "Rhodophyta", rank = "phylum", authority = "gbif"),
  descendants(name = "Phaeophyceae", rank = "phylum", authority = "gbif")
)

gbif_phyla
```


How many GBIF names also match the names given in OTT? Looks like only 2,923 exact matches.

```{r}
gbif_in_ott <- gbif_phyla %>% 
  select(gbif_id = taxonID, scientificName, taxonRank) %>%
  inner_join(ott_phyla)
```
---
title: "fuzzy-matching.Rmd"
author: "Carl Boettiger"
date: "1/3/2019"
output: html_document
---


```{r}
library(tidyverse)
library(taxadb)
```

```{r}
td_create("itis")
```


```{r}
bbs <- read_tsv(system.file("extdata/bbs.tsv", package="taxadb"))

```



```{r}
name <- bbs$species
authority <- "itis"
db <- td_connect()
match <- "contains"
```


# Strategy 1: SQL-based fuzzy match

```{r}
name_pattern <- switch(match,
                         starts_with = paste0(name, "%"),
                         contains =  paste0("%", name, "%")
  )

  system.time({
  out <- purrr::map_dfr(name_pattern,
          function(pattern)
            taxa_tbl(authority, "taxonid", db) %>%
            filter(name %like% pattern) %>% collect()
         )
  })
out
```


# Strategy 2: Extract a smaller table

```{r}
  ## Strategy: extract all potential matches by Genus alone.  assumes first name is a genus name!
  only_genus <- function(name)  stringi::stri_extract_first_words(name)
  id_tbl <- ids(only_genus(name), authority = authority, db = db, collect = FALSE) %>%
    select(name) %>%
    inner_join(select(taxa_tbl(authority, "hierarchy"), id, genus), by = c(name = "genus")) %>%
    select(id) %>% inner_join(taxa_tbl(authority, "taxonid"), by = "id") %>%
    distinct() %>%
    collect()

  name_regex <- switch(match,
                         starts_with = paste0(name, ".*"),
                         contains =  paste0(".*", name, ".*")
  )
  
  id_tbl <- collect(taxa_tbl(authority, "taxonid"))
  
  ## Using the genus subset -- a much smaller list of matches -- is this good or bad?
  system.time({
    out2 <- purrr::map_dfr(name_regex, function(pattern)
      filter(id_tbl, grepl(pattern, name))
    )
  })
  
  

```



```{r}   
  ## In memory, even slower!! 
  system.time({
    id_tbl <-  collect(taxa_tbl(authority, "taxonid"))
    out2 <- purrr::map_dfr(name_regex, function(pattern)
      filter(id_tbl, grepl(pattern, name))
    )
  })

```

Load libraries and data to get started

```{r}
library(tidyverse)
library(taxadb)
```

```{r}
taxadb::td_create("all") # only needed once.
```

```{r}
algae <- read_csv("~/projects/algae-names/algae_uncleanednames_NAdropped.csv")
```



## Summarizing the data

Before we begin, let's take a look around the data to get a sense of what we have:

- `r dim(algae)[[1]]` rows, `r names(algae)` columns.

```{r}
dim(algae)[[1]] # 57,700 rows. 
names(algae) # cols are: source, family, species, is_source_herbaria
```

```{r}
algae %>% count(family) %>% arrange(desc(n))  # 1,126 families (uncleaned)
```

```{r}
# 25 sources, including MACROALGAE, iDigBio, GBIF, OBIS, ...
algae %>% count(source) %>% arrange(desc(n))
```

lots of duplicate names across & sources: only 37,505 unique species.

```{r}
algae %>% count(species) %>% arrange(desc(n))
```

However, note that some of these duplicated names have different families though:

```{r}
unique_sp_family <- algae %>% count(species, family) %>% arrange(desc(n))
unique_sp_family  # 42,242 rows
```

We can see this is sometimes due to missing data or differences in capitalization, but also see cases where the same species name is assigned to different families:  

```{r}
species_multiple_families <- 
  unique_sp_family %>% 
  filter(species != "sp.", species != "Indet. sp.") %>%
  select(species) %>% 
  count(species) %>% 
  arrange(desc(n)) %>%
  left_join(unique_sp_family %>% select(-n))
species_multiple_families %>%  arrange(desc(n))
```

We will assume these differences correspond to changes in taxonomic group assignment at the family level, and not to two distinct species with the same scientific name (genus+specific epithet) belonging to different families.  

----------

## Resolving names to IDs

So we'll begin by focusing on resolving the unique species names:

```{r}
names <- unique(algae$species) # 37,505
```

Okay, without futher ado, let's start matching names.  We'll resolve against the Open Tree Taxonomy (OTT) first, because it is a assembly that already includes names from WORMS, GBIF, and others.  

```{r}
ott_ids <- names %>% ids("ott")
```

The resulting data frame has at least one row per name provided, though possibily all NA.  

This database includes synonyms, which do not have their own OTT identifiers in `taxonID`, but can be resolved to accepted names with accepted identifiers.  Thus, we want to look at `acceptedNameUsageID` column to see who matched and who didn't.  (This nomenclature for column names probably feels clumsy, but it comes from the Darwin Core standard and thus lets us be explicitly consistent with how these terms are used elsewhere.)

The `ids` function also automatically normalizes all strings to lowercase before matching against the (lowercased) scientificNames in the database.  `ids` returns these strings in the `input` column.  Standardizing between upper and lowercase reveals that many of the 37,505 names are still duplicates:

```{r}
length(ott_ids$input)
length(unique(ott_ids$input))
```

The `sort` column indicates the position of the original input data. Without this column, input names that are identical after correcting for capitalization result in identical rows.  We can drop these by filtering for only the distinct columns:

```{r}
clean_ott_ids <- ott_ids %>% select(-sort) %>% distinct()
clean_ott_ids # 29,878 rows
```


We can now check how many unique names were matched:


```{r}
matched <- clean_ott_ids %>% filter(!is.na(acceptedNameUsageID))
length(unique(matched$input)) # 18,568
length(unique(matched$acceptedNameUsageID)) #17,766

matched %>% filter(taxonomicStatus != "accepted") # 1,766
```

We matched 18,568 unique names (out of `r length(unique(clean_ott_ids$input))`), which resolved to 17,766 unique IDs in OTT, since 1,766 names were recognized synonyms to OTT.  Meanwhile, in the unmatched names set, this leaves us with:

```{r}
unmatched <- clean_ott_ids %>% filter(is.na(acceptedNameUsageID)) %>% pull(input)
length(unique(unmatched))
```

So 10,911 unmatched names still to resolve (a little more than one third of the unique names).  

```{r}
wordcount <- 
  data.frame(input = unmatched) %>% 
  mutate(n = stringi::stri_count_regex(input, "\\s")) %>% 
  arrange(desc(n))

wordcount
```

Wow, so some of our entries have a whole lot more words than a species name: the top 6 entries have dozens of words.  Many more (1,551) have 4 words (three spaces), indicating subspecies and varities:

```{r}
wordcount %>% filter(n==3)
```

If we restrict ourselves to matching on the first two-word names, we have a reasonable chance of resolving names to the genus or species level.  `clean_names` does up to three transformations: we drop missing specific epithet indication `sp`, leaving only the genus name, we resolve the first two names, and we standardize delimiters between names.  


```{r}
clean_unmatched <- unmatched %>% clean_names() %>% unique()
length(clean_unmatched)
```
This gives us `r length(clean_unmatched)` unmatched names to resolve.  

```{r}
ott_ids2 <- clean_unmatched %>% ids("ott")
unmatched2 <- ott_ids2 %>% filter(is.na(acceptedNameUsageID)) %>% pull(input)
```

We can bind the new matches to our existing matched names table:

```{r}
matched2 <- ott_ids2 %>% filter(!is.na(acceptedNameUsageID)) %>% bind_rows(matched)
```

With these clean names, another `r length(clean_unmatched) - length(umatched2)` can be matched to OTT ids, leaving `r length(umatched2)` unmatched:  


```{r}
length(unmatched2)
head(unmatched2, 10)
```

Working from the names directly:

```{r}

df <- db_mutate(r_fn = clean_names,
          tbl = "ott",
          db = td_connect(),
          col = "scientificName",
          new_column = "input") 


```


```{r}
sp <- tibble(input = unmatched2)

match_binomial <- right_join(df, sp, copy=TRUE, by = "input")  %>% collect()
match_binomial %>% filter(is.na(acceptedNameUsageID)) %>% distinct()

dim(match_binomial)
```


We can find at least some matches for these remaining names in most of the other databases:

```{r}
unmatched2 %>% ids("slb") %>% summarise(matched = sum(!is.na(acceptedNameUsageID)))
unmatched2 %>% ids("col") %>% summarise(matched = sum(!is.na(acceptedNameUsageID)))
unmatched2 %>% ids("gbif") %>% summarise(matched = sum(!is.na(acceptedNameUsageID)))
unmatched2 %>% ids("itis") %>% summarise(matched = sum(!is.na(acceptedNameUsageID)))
unmatched2 %>% ids("ncbi") %>% summarise(matched = sum(!is.na(acceptedNameUsageID)))
unmatched2 %>% ids("wd") %>% summarise(matched = sum(!is.na(acceptedNameUsageID)))
unmatched2 %>% ids("iucn") %>% summarise(matched = sum(!is.na(acceptedNameUsageID)))
unmatched2 %>% ids("tpl") %>% summarise(matched = sum(!is.na(acceptedNameUsageID)))

```
GBIF does the best among these, resolving another 1614 names exactly.  

```{r}
gbif_ids <- unmatched2 %>% ids("gbif")
unmatched3 <- gbif_ids %>% filter(is.na(acceptedNameUsageID)) %>% pull(input)
length(names(umatched3))
```


FIXME also try matching against `clean_names` versions of scientificName column on the database side; in particular, on the synonyms.  
---
title: "Crosswalking The Authorities"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{crosswalk}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

taxadb:::td_disconnect()
MonetDBLite::monetdblite_shutdown()

```



```{r message = FALSE}
library(taxadb)
library(dplyr)
library(tidyr)
```

Work in progress -- the mismatches of `itis`, `col` and `wd` here illustrate the complexity of name matching.  Some names resolve multiple times.  Some authorities (ITIS) recognize certain synonyms that they still haven't mapped to accepted id.  

```{r}
fish <- taxa_tbl("fb") %>%  
  select(fb = id, species) %>% 
  collect() 
species <- pull(fish, species)

itis <- ids(species, "itis")
wd <- ids(species, "wd")
col <- ids(species, "col")

## Note each of these authorities all return more rows than we had in the input table!
sapply(list(fb = fish, itis = itis, wd = wd, col = col), function(x) length(x[[1]]))

```


```{r}
my_taxa <- fish %>%

  mutate(#itis = get_ids(species, "itis"),
         ncbi = get_ids(species, "ncbi"),
         #col = get_ids(species, "col"),
         gbif = get_ids(species, "gbif"),
         #wd = get_ids(species, "wd"),
         tpl = get_ids(species, "tpl"),
        # fb = get_ids(species, "fb"),
         slb = get_ids(species, "slb")) 
```


```{r}
my_taxa %>% 
  select(-species) %>% 
  purrr::map_dbl(function(x) sum(!is.na(x)))
```

Looks like three plants have matching scientific names to some of our fish:

```{r}
dup <- fish %>% pull(species) %>% ids(authority = "tpl") %>% filter(!is.na(id))
dup
```

This also probably explains why `col` and `wd` are returning the wrong-length matches:

```{r}
species <- pull(fish, species)

col_fb <- ids(species, "col") 
dim(col_fb)[1] - length(species)

wd_fb <- ids(species, "wd") 
dim(wd_fb)[1] - length(species)

```

```{r}

#col_hierarchy <- classification(id = col_fb$id) %>% filter(kingdom == "Animalia") 

taxa_tbl("col", "hierarchy") %>%
  filter(kingdom == "Animalia") %>%
  semi_join(select(col_fb, id), copy = TRUE) %>%
  select(id, species) %>% 
  collect()
```



```{r}
has_match <- my_taxa %>% 
  select(-species, -fb) %>%
  purrr::map_dfc(function(x) !is.na(x)) %>% 
  rowSums() > 0

my_taxa %>% filter(!has_match)
```



```{r include=FALSE}
taxadb:::td_disconnect()
```

---
title: "taxadb overview"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{taxadb}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

taxadb:::td_disconnect()
MonetDBLite::monetdblite_shutdown()

```



```{r message = FALSE}
library(taxadb)
library(dplyr)
library(tidyr)
```



# Introduction: the problem of name matching

We frequently want to combine different data sources using species names.  Perhaps we have one table which gives us occurrance information for a given list of species, and another table that contains trait information for species, and maybe a third source that provides a phylogenetic tree.  If our system of scientific species names was perfect, we could simply do a "table join" using species name as the joining `key`, and all would be well. 

Unfortunately, as anyone who has attempted this kind of exercise over more than a handful of species has discovered, this often works for most but very rarely for all of the species names involved.  There are several reasons this process runs into problems.  

**Synonyms**.  One of the most common problems is the existence of synonyms: different names or different spellings of a given scientific name.  While this could include definite miss-spellings, for there are many species for which authorities can differ in the preferred name -- essentially, two or more different names correspond to the same species, and should be treated as such in our species join.  


For instance, we consider the 233 primate species names used in the `geiger` package `primates` data.  To avoid installing that dependency-heavy package, a copy of these names has been cached in `taxadb` and can be loaded as follows: 

```{r}
ex <- system.file("extdata", "primates.tsv.bz2", package = "taxadb", mustWork = TRUE)
primates <- readr::read_tsv(ex)
primates
```

Our goal will be to associate each name with a definitive taxonomic ID of an existing naming authority.  This will allow us to merge on IDs directly, rather than names.  By mapping both recognized names and synonyms to the same corresponding taxonomic ID, we can be sure that we can join the relevant data correctly.


## Setup

To get started, we'll install all the data sources available to `taxadb` into a local database. This may take a while, particularly over a slow internet connection, but it needs to be done only once.  The downloaded size of all data is around 3 GB.  Once this task completes, our subsequent operations should all be quite fast.  

```{r message=FALSE}
td_create(authorities = "all")
```


## Resolving names


Match a list of 233 species names against a naming authority:  


```{r}
my_taxa <- primates %>%  # 233 taxa
  mutate(id = get_ids(species, format = "prefix")) 
```


```{r}
my_taxa %>% filter(is.na(id))
```

Only 3 species are missing ids, so we have managed to resolve 230 of the 233 species.  Not bad!  Actually,two of these do appear in ITIS records, but map to no accepted name according to ITIS:

```{r}
unmatched <- my_taxa %>% filter(is.na(id)) %>% pull(species) 
synonyms(unmatched)

```

Under the hood, `get_ids` checks the names against synonyms known to the authority as well. If we only looke for direct matches, we would have faired far worse.  

Mimicking `taxize` function of the same name `get_ids` is returning just a vector of `ids`, which is useful for operations such as `mutate` above.  However, it is often convenient to return the full table of matched names.  Here we can see that many of the names we have matched against are actually synonyms, and that `get_ids` has been able to resolve them into accepted ids.  

```{r}
ids(primates$species)
```

Note that `name` is the name we provided (species names in this case, but it they do not have to be).  Because `synonyms-check` is `TRUE` by default, we get a list of accepted names recognized (by ITIS in this case, our default authority).  

We aren't stuck with ITIS though, we can work with other naming authorities as well:

```{r}
my_taxa <- primates %>%  
  mutate(itis = get_ids(species, "itis"),
         ncbi = get_ids(species, "ncbi"),
         col = get_ids(species, "col"),
         gbif = get_ids(species, "gbif"),
         wd = get_ids(species, "wd"),
         tpl = get_ids(species, "tpl"),
         fb = get_ids(species, "fb"),
         slb = get_ids(species, "slb")) 
my_taxa

```

Can any single authority resolve all species names?

```{r}
my_taxa %>% 
  select(-species) %>% 
  purrr::map_dbl(function(x) sum(!is.na(x)))
```

Looks like `itis` has the most matches with 230.  (Of course taxon-specific authorities like The Plant List, FishBase, and SeaLifeBase contain no primates).  

We can also ask: do any species have no match in any authority?

```{r}
has_match <- my_taxa %>% 
  select(-species) %>%
  purrr::map_dfc(function(x) !is.na(x)) %>% 
  rowSums() > 0

my_taxa %>% filter(!has_match)
```

One species still fails to resolve under any authority. A web search does resolve this as a misspelling, to `Leontopithecus chrysomelas`, but not one that is recongized automatically by the ITIS known synonyms. 

## Hierarchy

Once we have resolved our ids, we can also get full classification information. This example uses the default authority, ITIS:

```{r}
primate_ids <- get_ids(primates$species, format = "prefix")
classification(id = primate_ids)
```

the `classification` function can work on species names directly, but this will not resolve syonyms to identifiers first.  As a result, only the 180 already recognized species names can be resolved in this case:

```{r}
classification(species = primates$species)
```



```{r include=FALSE}
taxadb:::td_disconnect()
MonetDBLite::monetdblite_shutdown()
```

---
title: "Backends for taxadb"
author: "Carl Boettiger, Kari Norman"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{intro}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

```

`taxadb` is designed to work with a variety of different "backends" -- software that works under the hood to store and retrieve the requested data.  `taxadb` has an intelligent default method selector which will attempt to use the best method available on your system, which means you can use `taxadb` without having to worry about these details.  However, to improve performance of `taxadb`, becoming familiar with these backends can yield significant improvements in performance.

# RSQLite

`RSQLite` is the default database backend if no suggested backend is detected.  `RSQLite` has no external software dependencies and will be automatically installed with `taxadb` (it is a hard dependency as an imported rather than suggested package). The term `Lite` indicates that SQLite does not require the separate "server" and "client" software model found on traditional databases such as MySQL, and SQLite is widely used in consumer software everywhere.  RSQLite packages SQLite for R.  It enables persistent local storage for R applications but will be slower than the alternatives.  For certain operations it can be significantly slower.

# MonetDBLite & duckdb

`MonetDBLite` is a modern alternative to `RSQLite`.  `MonetDBLite` is both more powerful than SQLite (in supporting a greater array of operations), and can run much faster.  Filtering joins in particular can be much faster even than the in-memory operations of `dplyr`. Because filtering joins lie at the heart of many `taxadb` functions this can yield substantial improvements in performance.  Unfortunately, the R interface, `MonetDBLite` was removed from CRAN in April 2019. The package can still be installed from GitHub by running `devtools::install_github("hannesmuehleisen/MonetDBLite-R")`, though this requires the appropriate compilers.  The developer plans to replace MonetDBLite with `duckdb`, (see <https://github.com/duckdb/duckdb>), but this is not yet feature complete and thus not yet fully compatible for `taxadb` use.  Because installation is more difficult, `MonetDBLite` is not a required dependency, but will be used by default if `taxadb` detects an existing installation.  `duckdb` support will be switched on as the first priority in the method waterfall.   

# in-memory

`taxadb` can also be set to use in-memory only, without a backend.  (Note that this is distinct from using `RSQlite` or `MonetDBLite` with over `in-memory` mode, because it uses only native R `data.frame`s to store data).  This will tend to be faster that `RSQLite` but slower than `MonetDBLite` or `duckdb`.  In this mode, data will persist over a single session but not between sessions (since memory is cleared when the user quits out of R).  Note that many taxonomic tables are quite large when uncompressed, and users with less than 8-16GB of free RAM may find their machine becomes slow or unresponsive in this mode.  

# Manual control of the backend engine

Users can override the automatic preferences of `taxadb` by setting the environmental variable `TAXADB_DRIVER`.  For example, running `Sys.setenv(TAXADB_DRIVER="RSQLite")` will make `RSQLite` the default driver, even if `MonetDBLite` is installed.  

# Local storage

The first time `taxadb` accesses a data source, it will download and store the full dataset from that provider.  Users can trigger a download ahead of time by running `td_create()`, e.g. `td_create("fb")` will create a local copy of the FishBase taxonomy.  If a user does not call `td_create()` first, `taxadb` simply downloads the data the first time that provider is queried -- e.g. `filter_name("Homo sapiens", "gibf")` will first download and install GBIF if that has not been done already.  These download and install operations may be slow depending on your internet connection, but need be performed only once.  Downloaded data is stored on your local harddisk and will persist between R sessions. The default location depends on the default set by your operating system (see the `rappdirs` package).  Users can configure this location by setting the environmental variable `TAXADB_HOME`.  For example, all unit tests in the package use temporary storage by setting  `Sys.setenv(TAXADB_HOME=tempdir())`, which is cleared out after the R session ends. 

A user can install all available name providers up front with `td_create("all")`.  An overview of the available scientific name providers is found in the providers vignette.  


# Other backends

`taxadb` will work just as well with any `DBI`-compatible database backend (Postgres, MariaDB, etc).  All `taxadb` functions take an argument `taxadb_db`, which is just a `DBI` connection used by `dplyr`.  For example, we can create an in-memory RSQLite connection and use that to store data for a single session:

```r
con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
taxadb::get_ids("Homo sapiens", taxadb_db = con)
```

Users can also call the `td_connect()` function to connect to `taxadb`'s default databases.  Running `td_connect()` with no arguments will return the current default connection.  This is a convenient way to confirm that your system is using the database engine you intended it to use.  You can also use that connection to interact directly with the `taxadb` databases (e.g. using `dplyr` or `DBI` functions).  
---
title: "Name Providers and schema used in taxadb"
author: "Carl Boettiger, Kari Norman"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{schema}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```



`taxadb` relies on a set of pre-assembled tables following a set of standardized schema layouts using Darwin Core vocabulary, as outlined below.  The database dumps provided by providers supported in `taxadb` at this time are:

`taxadb` abbreviation | name 
----------------------|-------------------
`itis`   | The Integrated Taxonomic Information System, `https://www.itis.gov/`
`col`  | [The Catalogue of Life](http://www.catalogueoflife.org/)
`ncbi` |  [The National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/)
`gbif` | [The Global Biodiversity Information Facility](https://www.gbif.org/)
`tpl`  | [The Plant List](http://www.theplantlist.org/)
`fb`   | [FishBase](https://www.fishbase.de/)
`slb`  | [SeaLifeBase](https://www.sealifebase.ca/)
`wd`   | WikiData,  (wikidata.org)
`iucn` | The IUCN Red List of endangered species status, `https://www.iucnredlist.org`
`ott`  | [Open Tree of Life taxonomy](https://tree.opentreeoflife.org/about/taxonomy-version/ott3.1).  


***Please Note***: `taxadb` advises against uncritically combining data from multiple providers.  The same name is frequently used by different providers to mean different things -- some providers consider two names synonyms that other providers consider distinct species.  *It is crucial to recognize that taxonomic name providers represent independent taxonomic theories*, and not merely additional observations of the same immutable reality ([Franz & Sterner (2018)](https://doi.org/10.1093/database/bax100 "Nico M Franz, Beckett W Sterner, To increase trust, change the social design behind aggregated biodiversity data, Database, Volume 2018, 2018, bax100")). You cannot just merge two databases of taxonomic names like you can two databases of, say, plant traits to get a bigger and more complete sample, because the former can contain meaningful contradictions.  


At the same time, it is also important to note that `col`, `gbif`, `ott`, are explicitly synthesis projects integrating the databases of names from a range of (many) other providers, while `itis`, `iucn`, `ncbi`, `tpl`, `fb`, and `slb` are independent name providers.  The synthetic or integrated name lists are not simple merges, but the product of considerably expert opinion, and occasional nonsense automation. As such, they too represent novel (justified or otherwise) assertions of taxonomy, and are in no way a complete substitute for the databases they integrate, owing to both differences in how up-to-date the relative records are as well as to either expert disagreements or algorithmic miss-matches.  `taxadb` makes no attempt to provide an opinion or reconciliation mechanism to any of these issues, but only to provide convenient access to data and functions for manipulating these records in a fast and consistent manner.  (In fact, it is easy to use `taxadb` to verify that many of the names recognized in, say, ITIS, are not in fact included at all in Catalogue of Life or other databases that claim to derive from ITIS).

These providers also distribute taxonomic data in a wide range of database formats using a wide range of data layouts (schemas), not all of which are particularly easy to use or interpret (e.g. hierarchies are often but not always specified in `taxon_id,parent_id` pairs.)  To make it faster and easier to work across these providers, `taxadb` defines a common set of table schemas outlined below that are particularly suited for efficient computation of common tasks.  The `taxadb` format follows a strict interpretation of a subset of [Darwin Core](http://rs.tdwg.org/dwc).  `taxadb` pre-processes and publicly archives compressed, flat tables corresponding to each of these schema for each of these providers. The providers vary widely in the frequency at which they update their records, as well as whether they provide immutable versioned releases (e.g. `col`, `ott`), direct access to a database that is updated on a dynamic/continual basis without any log of the changes (`itis`, `ncbi`, others), or is simply unknown.  The `taxadb` maintainers take semi-annual snapshots and distribute versioned releases of the underlying data.  

Most common operations can be expressed in terms of standard database operations, such as simple filtering joins in SQL.  To implement these, `taxadb` imports the compressed flat files into a local, column-oriented database, which can be installed entirely as an R package with no additional server setup required.  This provides a persistent store, and ensures that operations can be performed on disk since the taxonomic tables considered here are frequently too large to store in active memory.  The columnar structure enables blazingly fast joins.  Once the database is created, `taxadb` simply wraps a set of user-friendly R functions around common `SQL` queries, implemented in the popular `dplyr` syntax.  By default, `taxadb` will always collect the results of these queries to return familiar, in-memory objects to the R user.  Optional arguments allow more direct access the database queries.  

## Data Schema

`taxadb` relies on the Simple Darwin Core Namespace for Taxon objects, <http://rs.tdwg.org/dwc/terms/> [@dwc].  This is the mostly widely recognized format for exchange of taxonomic information.

- `taxonID`: a unique id for the name (including provider prefix).  Note that some providers do not assign IDs to synonyms, but only to accepted names.  In this case, the `taxonID` should be `NA`, and the ID to the accepted name should be specified in `acceptedNameUsageID`.  
- `scientificName`: a Latin name, either accepted or known synonym, at the lowest resolved level for the taxon.  While DWC encourages the use of authorship citations, these are intentionally omitted in most tables as inconsistency in abbreviations and formatting make names with authors much harder to resolve.  When available, this information is provided in the additional optional columns using the corresponding Darwin Core terms.  ***Please note***: `scientificName` includes names at all taxonomic rank levels, it does not mean just "genus + specific epithet".  For example, "Animalia" is also a scientific name.  The `taxonRank` column indicates the associated taxonomic rank.  
- `taxonRank`: the rank (as given by the provider) of this taxon. **Please note**: While DarwinCore specifies seven ranks as separate columns (see below), many providers recognize many more of possible `taxonRank` values, such as "superclass", "superorder."  For example, NCBI (`ncbi`) and OpenTree Taxonomy (`ott`) recognize over 40 different ranks, many of which are unnamed, while Catalogue of Life (`col`), GBIF an others recognize only the seven principle ranks. Conflicting claims between naming providers about what rank a given name belongs to or what species are included in which rank are common.  
- `acceptedNameUsageID` the accepted identifier.  For synonyms, the scientificName of the row with the corresponding `taxonID` gives the accepted name, according to the data provider in question.  For accepted names, this is identical to the `taxonID` for the name. If not provided, it is assumed this is the same as the `taxonID`.  
- `taxonomicStatus` Either "accepted", for an accepted scientific name, or a term indicating if the name is a known synonym, common misspelling, etc.

Some providers may report additional optional columns, see below.  


## Hierarchy Terms

Darwin Core defines several commonly recognized ranks as possible Taxon properties as well: `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `specificEpithet`, and `intraspecificEpithet`.  Additionally, the taxonomic rank of any scientific name can be specified under `taxonRank`, whether or not it is one of these names.  

Semantically (specifically in the RDF sense), treating ranks as properties seems somewhat crude.  Database providers (and thus different experts) disagree both about what rank levels they recognize and what names belong in what ranks.  NCBI recognizes over 40 named ranks and numerous unnamed ranks.  OTT, in true cladistic fashion, identifies all mammals as being not only in the class "Mammalia", but also in the "class" of lobe-finned-fish, Sarcopterygii.  To distinguish between these different treatments, it would be semantically most consistent to associate a (or multiple) `taxonRankID` with each taxonomic entry, rather than a a taxonRank. This ID could be specific to the data provider, and indicate the rank name that provider associates with that rank.  Few (wikidata, with its strong RDF roots, is an exception) providers associate IDs with rank levels though. 

In practice, treating ranks as properties (i.e. as column headings) is far more consistent with typical scientific usage and convenient for common applications, such as generating a list of all birds or all frogs by a simple filter on names in a column. 


## Synonyms

The `taxonomicStatus` value indicates if the name provided is a synonym, misspelling or an accepted name.  `taxadb` does not enforce any controlled vocabulary on the use of these terms beyond using the term `accepted` to indicate that the `scientificName` is an accepted name (i.e. the `dwc:acceptedNameUsage`) for the taxon.  Including both accepted names and synonyms in the `scientificName` column greatly facilitates taxonomic name resolution: a user can just perform an SQL filtering join from a given list of names and the taxadb table in order to resolve names to identifiers (`acceptedNameUsageID`s).  


## Common names

Common names are available from several providers, but tidy tables for `taxadb` have not yet been implemented.  Common names tables are expected to follow the following schema:

- `id` The taxonomic identifier for the species (or possibly other rank)
- `name` The common name / vernacular name
- `language` The language in which the common name is given, if known. (all lowercase)
- `language_code` the two-letter language code.

## Linked Data formats

`taxadb` tables can easily be interpreted as semantic data and will be made available as RDF triples.  This permits the richer SPARQL-based queries of taxonomic information, in addition to the SQL-based queries.  This data format will be the focus of a separate R package interface `taxald`.  

## Conventions

- Identifiers use the integer identifier defined by the provider, prefixed by the provider abbreviation in all capital letters: `ITIS:`, `GBIF:`, etc.
- Rank names are always lower case without hyphens or spaces. Rank names should be mapped
  to a table of standard accepted rank names (i.e. those recognized by ITIS, NCBI, Wikidata),
  and rank names should have 
- Encoding is UTF-8

# Data Processing

A set of R scripts for pre-processing data from each of the names providers is included in the source code of `taxadb`, in the `data-raw/` sub-directory.  These scripts automate the process from download to generation of the cached copy accessed by the package.  While specific processing steps vary across providers, the most of the scripts focus on extracting a variety of formats (mostly SQLite and various text formats) and combining tables into a consistent implementation of Darwin Core following the schema and conventions outlined above, and writing this data out as compressed (bz2) tab-separated value files -- a cross-platform standard format that requires little specialized software to work with. Metadata regarding the provenance of each data file are also provided, including sha256 hashes of raw data (uncompressed data) are generated for cryptographic verification of data integrity.  

# Data Versioning

The above scripts are intended to be rerun annually to generate updated snapshots of the each of the data providers.  Each snapshot is then distributed as described above, as a separate cache release. All `taxadb` functions interacting with the provided taxonomic names data can specify which version (year) snapshot should be used, which facilitates reproducibility and easy comparisons across versions.  The scripts required to generate the data may be adjusted as needed if any of the naming providers change there own format over time.  Additional names providers may be added as opportunity presents. 
---
title: "Tutorial for taxadb"
author: "Carl Boettiger, Kari Norman"
date: "2020-02-06"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{intro}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




The goal of `taxadb` is to provide *fast*, *consistent* access to taxonomic data, supporting common tasks such as resolving taxonomic names to identifiers, looking up higher classification ranks of given species, or returning a list of all species below a given rank. These tasks are particularly common when synthesizing data across large species assemblies, such as combining occurrence records with trait records. 

Existing approaches to these problems typically rely on web APIs, which can make them impractical for work with large numbers of species or in more complex pipelines.  Queries and returned formats also differ across the different taxonomic authorities, making tasks that query multiple authorities particularly complex. `taxadb` creates a *local* database of most readily available taxonomic authorities, each of which is transformed into consistent, standard, and researcher-friendly tabular formats.  


## Install and initial setup

To get started, install the development version directly from GitHub:


```r
devtools::install_github("ropensci/taxadb")
```



```r
library(taxadb)
library(dplyr) # Used to illustrate how a typical workflow combines nicely with `dplyr`
```

Create a local copy of the Catalogue of Life (2018) database: 


```r
td_create("col")
#> Importing 2019_common_col.tsv.bz2 in 100000 line chunks:
#> 
[-] chunk 2
[\] chunk 3
[|] chunk 4
[/] chunk 5	...Done! (in 9.042107 secs)
```


Read in the species list used by the Breeding Bird Survey:


```r
bbs_species_list <- system.file("extdata/bbs.tsv", package="taxadb")
bbs <- read.delim(bbs_species_list)
```

## Getting names and ids

Two core functions are `get_ids()` and `get_names()`.  These functions take a vector of names or ids (respectively), and return a vector of ids or names (respectively).  For instance, we can use this to attempt to resolve all the bird names in the Breeding Bird Survey against the Catalogue of Life:



```r
birds <- bbs %>% 
  select(species) %>% 
  mutate(id = get_ids(species, "col"))

head(birds, 10)
#>                          species           id
#> 1         Dendrocygna autumnalis COL:35517330
#> 2            Dendrocygna bicolor COL:35517332
#> 3                Anser canagicus COL:35517329
#> 4             Anser caerulescens COL:35517325
#> 5  Chen caerulescens (blue form)         <NA>
#> 6                   Anser rossii COL:35517328
#> 7                Anser albifrons COL:35517308
#> 8                Branta bernicla COL:35517301
#> 9      Branta bernicla nigricans COL:35537100
#> 10             Branta hutchinsii COL:35536445
```

Note that some names cannot be resolved to an identifier.  This can occur because of miss-spellings, non-standard formatting, or the use of a synonym not recognized by the naming provider.  Names that cannot be uniquely resolved because they are known synonyms of multiple different species will also return `NA`.  The `filter_name` filtering functions can help us resolve this last case (see below).

`get_ids()` returns the IDs of accepted names, that is `dwc:AcceptedNameUsageID`s.  We can resolve the IDs into accepted names:



```r
birds %>% 
  mutate(accepted_name = get_names(id, "col")) %>% 
  head()
#>                         species           id        accepted_name
#> 1        Dendrocygna autumnalis COL:35517330      Tringa flavipes
#> 2           Dendrocygna bicolor COL:35517332    Picoides dorsalis
#> 3               Anser canagicus COL:35517329   Setophaga castanea
#> 4            Anser caerulescens COL:35517325  Bombycilla cedrorum
#> 5 Chen caerulescens (blue form)         <NA>       Icteria virens
#> 6                  Anser rossii COL:35517328 Somateria mollissima
```

This illustrates that some of our names, e.g. *Dendrocygna bicolor* are accepted in the Catalogue of Life, while others, *Anser canagicus* are **known synonyms** of a different accepted name: **Chen canagica**.  Resolving synonyms and accepted names to identifiers helps us avoid the possible miss-matches we could have when the same species is known by two different names.


## Taxonomic Data Tables

Local access to taxonomic data tables lets us do much more than look up names and ids.  A family of `filter_*` functions in `taxadb` help us work directly with subsets of the taxonomic data.  As we noted above, this can be useful in resolving certain ambiguous names.  

For instance, *Trochalopteron henrici gucenense* does not resolve to an identifier in ITIS:


```r
get_ids("Trochalopteron henrici gucenense") 
#> [1] NA
```

Using `filter_name()`, we find this is because the name resolves not to zero matches, but to more than one match:


```r
filter_name("Trochalopteron henrici gucenense") 
#> # A tibble: 2 x 17
#>    sort taxonID scientificName taxonRank acceptedNameUsa… taxonomicStatus update_date kingdom phylum class order family genus
#>   <int> <chr>   <chr>          <chr>     <chr>            <chr>           <chr>       <chr>   <chr>  <chr> <chr> <chr>  <chr>
#> 1     1 ITIS:9… Trochaloptero… subspeci… ITIS:916117      synonym         <NA>        Animal… Chord… Aves  Pass… Leiot… Troc…
#> 2     1 ITIS:9… Trochaloptero… subspeci… ITIS:916116      synonym         <NA>        Animal… Chord… Aves  Pass… Leiot… Troc…
#> # … with 4 more variables: specificEpithet <chr>, vernacularName <chr>, infraspecificEpithet <chr>, input <chr>
```



```r
filter_name("Trochalopteron henrici gucenense")  %>%
  mutate(acceptedNameUsage = get_names(acceptedNameUsageID)) %>% 
  select(scientificName, taxonomicStatus, acceptedNameUsage, acceptedNameUsageID)
#> # A tibble: 2 x 4
#>   scientificName                   taxonomicStatus acceptedNameUsage       acceptedNameUsageID
#>   <chr>                            <chr>           <chr>                   <chr>              
#> 1 Trochalopteron henrici gucenense synonym         Trochalopteron henrici  ITIS:916117        
#> 2 Trochalopteron henrici gucenense synonym         Trochalopteron elliotii ITIS:916116
```


Similar functions `filter_id`, `filter_rank`, and `filter_common` take IDs, scientific ranks, or common names, respectively.  Here, we can get taxonomic data on all bird names in the Catalogue of Life:



```r
filter_rank(name = "Aves", rank = "class", provider = "col")
#> # A tibble: 35,398 x 21
#>     sort taxonID scientificName acceptedNameUsa… taxonomicStatus taxonRank kingdom phylum class order family genus
#>    <int> <chr>   <chr>          <chr>            <chr>           <chr>     <chr>   <chr>  <chr> <chr> <chr>  <chr>
#>  1     1 COL:35… Sturnella mag… COL:35520416     accepted        species   Animal… Chord… Aves  Pass… Icter… Stur…
#>  2     1 COL:35… Tauraco porph… COL:35530219     accepted        infraspe… Animal… Chord… Aves  Muso… Musop… Taur…
#>  3     1 COL:35… Pyroderus scu… COL:35534370     accepted        infraspe… Animal… Chord… Aves  Pass… Cotin… Pyro…
#>  4     1 COL:35… Dromaius minor COL:35552206     synonym         infraspe… Animal… Chord… Aves  Casu… Droma… Drom…
#>  5     1 COL:35… Lepidocolapte… COL:35525495     accepted        species   Animal… Chord… Aves  Pass… Furna… Lepi…
#>  6     1 COL:35… Casuarius pap… COL:35552204     synonym         infraspe… Animal… Chord… Aves  Casu… Casua… Casu…
#>  7     1 COL:35… Forpus modest… COL:35536431     accepted        species   Animal… Chord… Aves  Psit… Psitt… Forp…
#>  8     1 COL:35… Pterocnemia p… COL:35552203     synonym         infraspe… Animal… Chord… Aves  Rhei… Rheid… Rhea 
#>  9     1 COL:35… Ceyx lepidus … COL:35532279     accepted        infraspe… Animal… Chord… Aves  Cora… Alced… Ceyx 
#> 10     1 COL:35… Rhea tarapace… COL:35552202     synonym         infraspe… Animal… Chord… Aves  Rhei… Rheid… Rhea 
#> # … with 35,388 more rows, and 9 more variables: specificEpithet <chr>, infraspecificEpithet <chr>, taxonConceptID <chr>,
#> #   isExtinct <chr>, nameAccordingTo <chr>, namePublishedIn <chr>, scientificNameAuthorship <chr>, vernacularName <chr>,
#> #   input <chr>
```

Combining these with `dplyr` functions can make it easy to explore this data: for instance, which families have the most species?



```r
filter_rank(name = "Aves", rank = "class", provider = "col") %>%
  filter(taxonomicStatus == "accepted", taxonRank=="species") %>% 
  group_by(family) %>%
  count(sort = TRUE) %>% 
  head()
#> # A tibble: 6 x 2
#> # Groups:   family [6]
#>   family           n
#>   <chr>        <int>
#> 1 Tyrannidae     401
#> 2 Thraupidae     374
#> 3 Psittacidae    370
#> 4 Trochilidae    338
#> 5 Muscicapidae   314
#> 6 Columbidae     312
```

## Using the database connection directly

`filter_*` functions by default return in-memory data frames.  Because they are filtering functions, they return a subset of the full data which matches a given query (names, ids, ranks, etc), so the returned data.frames are smaller than the full record of a naming provider.  Working directly with the SQL connection to the MonetDBLite database gives us access to all the data. The `taxa_tbl()` function provides this connection:


```r
taxa_tbl("col")
#> # Source:   table<2019_dwc_col> [?? x 19]
#> # Database: duckdb_connection
#>    taxonID scientificName acceptedNameUsa… taxonomicStatus taxonRank kingdom phylum class order family genus specificEpithet
#>    <chr>   <chr>          <chr>            <chr>           <chr>     <chr>   <chr>  <chr> <chr> <chr>  <chr> <chr>          
#>  1 COL:31… Limacoccus br… COL:316423       accepted        species   Animal… Arthr… Inse… Hemi… Beeso… Lima… brasiliensis   
#>  2 COL:31… Coccus bromel… COL:316424       accepted        species   Animal… Arthr… Inse… Hemi… Cocci… Cocc… bromeliae      
#>  3 COL:31… Apiomorpha po… COL:316425       accepted        species   Animal… Arthr… Inse… Hemi… Erioc… Apio… pomaphora      
#>  4 COL:31… Eriococcus ch… COL:316441       accepted        species   Animal… Arthr… Inse… Hemi… Erioc… Erio… chaoticus      
#>  5 COL:31… Eriococcus ch… COL:316442       accepted        species   Animal… Arthr… Inse… Hemi… Erioc… Erio… chathamensis   
#>  6 COL:31… Eriococcus ch… COL:316443       accepted        species   Animal… Arthr… Inse… Hemi… Erioc… Erio… chilensis      
#>  7 COL:31… Eriococcus ci… COL:316444       accepted        species   Animal… Arthr… Inse… Hemi… Erioc… Erio… cingulatus     
#>  8 COL:31… Eriococcus ci… COL:316445       accepted        species   Animal… Arthr… Inse… Hemi… Erioc… Erio… cistacearum    
#>  9 COL:31… Eriococcus bu… COL:316447       accepted        species   Animal… Arthr… Inse… Hemi… Erioc… Erio… busariae       
#> 10 COL:31… Eriococcus au… COL:316450       accepted        species   Animal… Arthr… Inse… Hemi… Erioc… Erio… australis      
#> # … with more rows, and 7 more variables: infraspecificEpithet <chr>, taxonConceptID <chr>, isExtinct <chr>,
#> #   nameAccordingTo <chr>, namePublishedIn <chr>, scientificNameAuthorship <chr>, vernacularName <chr>
```

We can still use most familiar `dplyr` verbs to perform common tasks.  For instance: which species has the most known synonyms?


```r
taxa_tbl("col") %>% 
  count(acceptedNameUsageID, sort=TRUE)
#> # Source:     lazy query [?? x 2]
#> # Database:   duckdb_connection
#> # Ordered by: desc(n)
#>    acceptedNameUsageID     n
#>    <chr>               <dbl>
#>  1 COL:43082445          456
#>  2 COL:43081989          373
#>  3 COL:43124375          329
#>  4 COL:43353659          328
#>  5 COL:43223150          322
#>  6 COL:43337824          307
#>  7 COL:43124158          302
#>  8 COL:43081973          296
#>  9 COL:43333057          253
#> 10 COL:23162697          252
#> # … with more rows
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_rank.R
\name{filter_rank}
\alias{filter_rank}
\title{Get all members (descendants) of a given rank level}
\usage{
filter_rank(
  name,
  rank,
  provider = getOption("taxadb_default_provider", "itis"),
  version = latest_version(),
  collect = TRUE,
  ignore_case = TRUE,
  db = td_connect()
)
}
\arguments{
\item{name}{taxonomic scientific name (e.g. "Aves")}

\item{rank}{taxonomic rank name. (e.g. "class")}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{collect}{logical, default \code{TRUE}. Should we return an in-memory
data.frame (default, usually the most convenient), or a reference to
lazy-eval table on disk (useful for very large tables on which we may
first perform subsequent filtering operations.)}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}

\item{db}{a connection to the taxadb database. See details.}
}
\value{
a data.frame in the Darwin Core tabular format containing the
matching taxonomic entities.
}
\description{
Get all members (descendants) of a given rank level
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }

filter_rank("Aves", "class")

}

}
\seealso{
Other filter_by: 
\code{\link{filter_by}()},
\code{\link{filter_common}()},
\code{\link{filter_id}()},
\code{\link{filter_name}()}
}
\concept{filter_by}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/td_create.R
\name{td_create}
\alias{td_create}
\title{create a local taxonomic database}
\usage{
td_create(
  provider = getOption("taxadb_default_provider", "itis"),
  schema = c("dwc", "common"),
  version = latest_version(),
  overwrite = TRUE,
  lines = 1e+05,
  dbdir = taxadb_dir(),
  db = td_connect(dbdir)
)
}
\arguments{
\item{provider}{a list (character vector) of provider to be included in the
database. By default, will install \code{itis}.  See details for a list of
recognized provider. Use \code{provider="all"} to install all
available provider automatically.}

\item{schema}{One of "dwc" (for Darwin Core data) or "common"
(for the Common names table.)}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{overwrite}{Should we overwrite existing tables? Default is \code{TRUE}.
Change to "ask" for interactive interface, or \code{TRUE} to force overwrite
(i.e. updating a local database upon new release.)}

\item{lines}{number of lines that can be safely read in to memory at once.
Leave at default or increase for faster importing if you have
plenty of spare RAM.}

\item{dbdir}{a location on your computer where the database
should be installed. Defaults to user data directory given by
\verb{[rappdirs::user_data_dir]}.}

\item{db}{connection to a database.  By default, taxadb will set up its own
fast database connection.}
}
\value{
path where database has been installed (invisibly)
}
\description{
create a local taxonomic database
}
\details{
Authorities currently recognized by taxadb are:
\itemize{
\item \code{itis}: Integrated Taxonomic Information System, \verb{https://www.itis.gov}
\item \code{ncbi}:  National Center for Biotechnology Information,
\url{https://www.ncbi.nlm.nih.gov/taxonomy}
\item \code{col}: Catalogue of Life, \url{http://www.catalogueoflife.org/}
\item \code{tpl}: The Plant List, \url{http://www.theplantlist.org/}
\item \code{gbif}: Global Biodiversity Information Facility, \url{https://www.gbif.org/}
\item \code{fb}: FishBase, \url{https://www.fishbase.de/}
\item \code{slb}: SeaLifeBase, \url{http://sealifebase.org}
\item \code{wd}: Wikidata: https://www.wikidata.org
\item \code{ott}: OpenTree Taxonomy:
\url{https://github.com/OpenTreeOfLife/reference-taxonomy}
\item \code{iucn}: IUCN Red List, https://iucnredlist.org
\item \code{itis_test}: a small subset of ITIS, cached locally with the package for testing purposes only
}
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")

  }
  ## Install the ITIS database
  td_create()

  ## force re-install:
  td_create( overwrite = TRUE)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuzzy_filter.R
\name{name_starts_with}
\alias{name_starts_with}
\title{scientific name starts with}
\usage{
name_starts_with(
  name,
  provider = getOption("taxadb_default_provider", "itis"),
  version = latest_version(),
  db = td_connect(),
  ignore_case = TRUE
)
}
\arguments{
\item{name}{vector of names (scientific or common, see \code{by})
to be matched against.}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{db}{a connection to the taxadb database. See details.}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}
}
\description{
scientific name starts with
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }
name_starts_with("Chera")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tl_import.R
\name{tl_import}
\alias{tl_import}
\title{Import taxonomic database tables}
\usage{
tl_import(
  provider = getOption("tl_default_provider", "itis"),
  schema = c("dwc", "common"),
  version = latest_version(),
  prov = paste0("https://raw.githubusercontent.com/",
    "boettiger-lab/taxadb-cache/master/prov.json")
)
}
\arguments{
\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{schema}{One of "dwc" (for Darwin Core data) or "common"
(for the Common names table.)}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{prov}{Address (URL) to provenance record}
}
\value{
path(s) to the downloaded files in the cache
}
\description{
Downloads the requested taxonomic data tables and return a local path
to the data in \code{tsv.gz} format.  Downloads are cached and identified by
content hash so that \code{tl_import} will not attempt to download the
same file multiple times.
}
\details{
\code{tl_import} parses a DCAT2/PROV-O record to determine the correct version
to download. If offline, \code{tl_import} will attempt to resolve against
it's own provenance cache. Users can also examine / parse the prov
JSON-LD file directly to determine the provenance of the data products
used.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/td_connect.R
\name{td_disconnect}
\alias{td_disconnect}
\title{Disconnect from the taxadb database.}
\usage{
td_disconnect(db = td_connect())
}
\arguments{
\item{db}{database connection}
}
\description{
Disconnect from the taxadb database.
}
\details{
This function manually closes a connection to the \code{taxadb} database.
}
\examples{
\donttest{

## Disconnect from the database:
td_disconnect()

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutate_db.R
\name{mutate_db}
\alias{mutate_db}
\title{Add new variables to a database}
\usage{
mutate_db(.data, r_fn, col, new_column, n = 5000L, ...)
}
\arguments{
\item{.data}{A \link[dplyr:tbl]{dplyr::tbl} that uses a database connection, \code{tbl_dbi} class.}

\item{r_fn}{any R function that can be called on a vector (column)
of the table}

\item{col}{the name of the column to which the R function is applied.
(Note, \code{\link[dplyr:mutate]{dplyr::mutate()}} can operate on an arbitrary list of columns,
this function only operates on a single column at this time...)}

\item{new_column}{column name for the new column.}

\item{n}{the number of rows included in each chunk, see \code{\link[DBI:dbFetch]{DBI::dbFetch()}}}

\item{...}{named arguments to be passed to \code{r_fn}}
}
\value{
a dplyr tbl connection to the temporary table in the database
}
\description{
\code{\link[dplyr:mutate]{dplyr::mutate()}} cannot pass arbitrary R functions over a
database connection. This function provides a way to work
around this, by querying the data in chunks
and applying the function to each chunk, which is then
appended back out to a temporary table.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_id.R
\name{filter_id}
\alias{filter_id}
\title{Return a taxonomic table matching the requested ids}
\usage{
filter_id(
  id,
  provider = getOption("taxadb_default_provider", "itis"),
  type = c("taxonID", "acceptedNameUsageID"),
  version = latest_version(),
  collect = TRUE,
  db = td_connect()
)
}
\arguments{
\item{id}{taxonomic id, in prefix format}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{type}{id type.  Can be \code{taxonID} or \code{acceptedNameUsageID},
see details.}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{collect}{logical, default \code{TRUE}. Should we return an in-memory
data.frame (default, usually the most convenient), or a reference to
lazy-eval table on disk (useful for very large tables on which we may
first perform subsequent filtering operations.)}

\item{db}{a connection to the taxadb database. See details.}
}
\value{
a data.frame with id and name of all matching species
}
\description{
Return a taxonomic table matching the requested ids
}
\details{
Use \code{type="acceptedNameUsageID"} to return all rows
for which this ID is the accepted ID, including both synonyms and
and accepted names (since both all synonyms of a name share the
same \code{acceptedNameUsageID}.) Use \code{taxonID} (default) to only return
those rows for which the Scientific name corresponds to the \code{taxonID.}

Some providers (e.g. ITIS) assign taxonIDs to synonyms, most others
only assign IDs to accepted names.  In the latter case, this means
requesting \code{taxonID} will only match accepted names, while requesting
matches to the \code{acceptedNameUsageID} will also return any known synonyms.
See examples.
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }

filter_id(c("ITIS:1077358", "ITIS:175089"))
filter_id("ITIS:1077358", type="acceptedNameUsageID")

}
}
\seealso{
Other filter_by: 
\code{\link{filter_by}()},
\code{\link{filter_common}()},
\code{\link{filter_name}()},
\code{\link{filter_rank}()}
}
\concept{filter_by}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuzzy_filter.R
\name{name_contains}
\alias{name_contains}
\title{return all taxa in which scientific name contains the text provided}
\usage{
name_contains(
  name,
  provider = getOption("taxadb_default_provider", "itis"),
  version = latest_version(),
  db = td_connect(),
  ignore_case = TRUE
)
}
\arguments{
\item{name}{vector of names (scientific or common, see \code{by})
to be matched against.}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{db}{a connection to the taxadb database. See details.}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}
}
\description{
return all taxa in which scientific name contains the text provided
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }
name_contains("Chera")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_names.R
\name{clean_names}
\alias{clean_names}
\title{Clean taxonomic names}
\usage{
clean_names(
  names,
  fix_delim = TRUE,
  binomial_only = TRUE,
  remove_sp = TRUE,
  ascii_only = TRUE,
  lowercase = TRUE,
  remove_punc = FALSE
)
}
\arguments{
\item{names}{a character vector of taxonomic names (usually species names)}

\item{fix_delim}{Should we replace separators \code{.}, \verb{_}, \code{-}
with spaces? e.g. 'Homo.sapiens' becomes 'Homo sapiens'.
logical, default TRUE.}

\item{binomial_only}{Attempt to prune name to a binomial name, e.g.
Genus and species (specific epithet), e.g. \verb{Homo sapiens sapiens}
becomes \verb{Homo sapiens}. logical, default \link{TRUE}.}

\item{remove_sp}{Should we drop unspecified species epithet designations?
e.g. \verb{Homo sp.} becomes \code{Homo} (thus only matching against genus level ids).
logical, default \link{TRUE}.}

\item{ascii_only}{should we coerce strings to ascii characters?
(see \code{\link[stringi:stri_trans_general]{stringi::stri_trans_general()}})}

\item{lowercase}{should names be coerced to lower-case to provide
case-insensitive matching?}

\item{remove_punc}{replace all punctuation but apostrophes with a space,
remove apostrophes}
}
\description{
A utility to sanitize taxonomic names to increase
probability of resolving names.
}
\details{
Current implementation is limited to handling a few
common cases. Additional extensions may be added later.
A goal of the \code{clean_names} function is that any
modification rule of the name strings be precise, atomic, and
toggle-able, rather than relying on clever but more opaque rules and
arbitrary scores. This utility should always be used with care, as
indiscriminate modification of names may result in successful but inaccurate
name matching. A good pattern is to only apply this function to the subset
of names that cannot be directly matched.
}
\examples{
clean_names(c("Homo sapiens sapiens", "Homo.sapiens", "Homo sp."))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_by.R
\name{filter_by}
\alias{filter_by}
\title{Creates a data frame with column name given by \code{by}, and values given
by the vector \code{x}, and then uses this table to do a filtering join,
joining on the \code{by} column to return all rows matching the \code{x} values
(scientificNames, taxonIDs, etc).}
\usage{
filter_by(
  x,
  by,
  provider = getOption("taxadb_default_provider", "itis"),
  schema = c("dwc", "common"),
  version = latest_version(),
  collect = TRUE,
  db = td_connect(),
  ignore_case = TRUE
)
}
\arguments{
\item{x}{a vector of values to filter on}

\item{by}{a column name in the taxa_tbl (following Darwin Core Schema terms).
The filtering join is executed with this column as the joining variable.}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{schema}{One of "dwc" (for Darwin Core data) or "common"
(for the Common names table.)}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{collect}{logical, default \code{TRUE}. Should we return an in-memory
data.frame (default, usually the most convenient), or a reference to
lazy-eval table on disk (useful for very large tables on which we may
first perform subsequent filtering operations.)}

\item{db}{a connection to the taxadb database. See details.}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}
}
\value{
a data.frame in the Darwin Core tabular format containing the
matching taxonomic entities.
}
\description{
Creates a data frame with column name given by \code{by}, and values given
by the vector \code{x}, and then uses this table to do a filtering join,
joining on the \code{by} column to return all rows matching the \code{x} values
(scientificNames, taxonIDs, etc).
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }

sp <- c("Trochalopteron henrici gucenense",
        "Trochalopteron elliotii")
filter_by(sp, "scientificName")

filter_by(c("ITIS:1077358", "ITIS:175089"), "taxonID")

filter_by("Aves", "class")

}

}
\seealso{
Other filter_by: 
\code{\link{filter_common}()},
\code{\link{filter_id}()},
\code{\link{filter_name}()},
\code{\link{filter_rank}()}
}
\concept{filter_by}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa_tbl.R
\name{taxa_tbl}
\alias{taxa_tbl}
\title{Return a reference to a given table in the taxadb database}
\usage{
taxa_tbl(
  provider = getOption("taxadb_default_provider", "itis"),
  schema = c("dwc", "common"),
  version = latest_version(),
  db = td_connect()
)
}
\arguments{
\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{schema}{One of "dwc" (for Darwin Core data) or "common"
(for the Common names table.)}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{db}{a connection to the taxadb database. See details.}
}
\description{
Return a reference to a given table in the taxadb database
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }

  ## default schema is the dwc table
  taxa_tbl()

  ## common names table
  taxa_tbl(schema = "common")



}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_common.R
\name{filter_common}
\alias{filter_common}
\title{Look up taxonomic information by common name}
\usage{
filter_common(
  name,
  provider = getOption("taxadb_default_provider", "itis"),
  version = latest_version(),
  collect = TRUE,
  ignore_case = TRUE,
  db = td_connect()
)
}
\arguments{
\item{name}{a character vector of common (vernacular English) names,
e.g. "Humans"}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{collect}{logical, default \code{TRUE}. Should we return an in-memory
data.frame (default, usually the most convenient), or a reference to
lazy-eval table on disk (useful for very large tables on which we may
first perform subsequent filtering operations.)}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}

\item{db}{a connection to the taxadb database. See details.}
}
\value{
a data.frame in the Darwin Core tabular format containing the
matching taxonomic entities.
}
\description{
Look up taxonomic information by common name
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   options("taxadb_default_provider"="itis_test")
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
  }

filter_common("Pied Tamarin")

}

  \dontshow{
   ## All examples use a temporary directory
   Sys.unsetenv("TAXADB_HOME")
   options("taxadb_default_provider" = NULL)
  }
}
\seealso{
Other filter_by: 
\code{\link{filter_by}()},
\code{\link{filter_id}()},
\code{\link{filter_name}()},
\code{\link{filter_rank}()}
}
\concept{filter_by}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxadb_dir.R
\name{taxadb_dir}
\alias{taxadb_dir}
\title{Show the taxadb directory}
\usage{
taxadb_dir()
}
\description{
Show the taxadb directory
}
\details{
NOTE: after upgrading \code{duckdb}, a user may need to delete any
existing databases created with the previous version. An efficient
way to do so is \code{unlink(taxadb::taxadb_dir(), TRUE)}.
}
\examples{
## show the directory
taxadb_dir()
## Purge the local db
unlink(taxadb::taxadb_dir(), TRUE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuzzy_filter.R
\name{common_contains}
\alias{common_contains}
\title{common name starts with}
\usage{
common_contains(
  name,
  provider = getOption("taxadb_default_provider", "itis"),
  version = latest_version(),
  db = td_connect(),
  ignore_case = TRUE
)
}
\arguments{
\item{name}{vector of names (scientific or common, see \code{by})
to be matched against.}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{db}{a connection to the taxadb database. See details.}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}
}
\description{
common name starts with
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }
common_contains("monkey")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuzzy_filter.R
\name{common_starts_with}
\alias{common_starts_with}
\title{common name starts with}
\usage{
common_starts_with(
  name,
  provider = getOption("taxadb_default_provider", "itis"),
  version = latest_version(),
  db = td_connect(),
  ignore_case = TRUE
)
}
\arguments{
\item{name}{vector of names (scientific or common, see \code{by})
to be matched against.}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{db}{a connection to the taxadb database. See details.}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}
}
\description{
common name starts with
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }
common_starts_with("monkey")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuzzy_filter.R
\name{fuzzy_filter}
\alias{fuzzy_filter}
\title{Match names that start or contain a specified text string}
\usage{
fuzzy_filter(
  name,
  by = c("scientificName", "vernacularName"),
  provider = getOption("taxadb_default_provider", "itis"),
  match = c("contains", "starts_with"),
  version = latest_version(),
  db = td_connect(),
  ignore_case = TRUE,
  collect = TRUE
)
}
\arguments{
\item{name}{vector of names (scientific or common, see \code{by})
to be matched against.}

\item{by}{a column name in the taxa_tbl (following Darwin Core Schema terms).
The filtering join is executed with this column as the joining variable.}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{match}{should we match by names starting with the term or containing
the term anywhere in the name?}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{db}{a connection to the taxadb database. See details.}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}

\item{collect}{logical, default \code{TRUE}. Should we return an in-memory
data.frame (default, usually the most convenient), or a reference to
lazy-eval table on disk (useful for very large tables on which we may
first perform subsequent filtering operations.)}
}
\description{
Match names that start or contain a specified text string
}
\details{
Note that fuzzy filter will be fast with an single or small number
of names, but will be slower if given a very large vector of
names to match, as unlike other \code{filter_} commands,
fuzzy matching requires separate SQL calls for each name.
As fuzzy matches should all be confirmed manually in any event, e.g.
not every common name containing "monkey" belongs to a primate species.

This method utilizes the database operation \verb{\%like\%} to filter tables without
loading into memory.  Note that this does not support the use of regular
expressions at this time.
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }

## match any common name containing:
name <- c("woodpecker", "monkey")
fuzzy_filter(name, "vernacularName")

## match scientific name
fuzzy_filter("Chera", "scientificName",
             match = "starts_with")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ids.R
\name{get_ids}
\alias{get_ids}
\title{get_ids}
\usage{
get_ids(
  names,
  db = getOption("taxadb_default_provider", "itis"),
  format = c("prefix", "bare", "uri"),
  version = latest_version(),
  taxadb_db = td_connect(),
  ignore_case = TRUE,
  warn = TRUE,
  ...
)
}
\arguments{
\item{names}{a list of scientific names (which may
include higher-order ranks in most authorities).}

\item{db}{abbreviation code for the provider.  See details.}

\item{format}{Format for the returned identifier, one of
\itemize{
\item \code{prefix} (e.g. \code{NCBI:9606}, the default), or
\item \code{bare} (e.g. \code{9606}, used in \code{taxize::get_ids()}),
\item \code{uri} (e.g.
\verb{http://ncbi.nlm.nih.gov/taxonomy/9606}).
}}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  see \verb{[avialable_releases()]} for details.}

\item{taxadb_db}{Connection to from \verb{[td_connect()]}.}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}

\item{warn}{should we display warnings on NAs resulting from multiply-resolved matches?
(Unlike unmatched names, these NAs can usually be resolved manually via \code{\link[=filter_id]{filter_id()}})}

\item{...}{additional arguments (currently ignored)}
}
\value{
a vector of IDs, of the same length as the input names Any
unmatched names or multiply-matched names will return as \link{NA}s.
To resolve multi-matched names, use \verb{[filter_name()]} instead to return
a table with a separate row for each separate match of the input name.
}
\description{
A drop-in replacement for \verb{[taxize::get_ids()]}
}
\details{
Note that some taxize authorities: \code{nbn}, \code{tropicos}, and \code{eol},
are not recognized by taxadb and will throw an error here. Meanwhile,
taxadb recognizes several authorities not known to \verb{[taxize::get_ids()]}.
Both include \code{itis}, \code{ncbi}, \code{col}, and \code{gbif}.

Like all taxadb functions, this function will run
fastest if a local copy of the provider is installed in advance
using \verb{[td_create()]}.
}
\examples{
\donttest{

  \dontshow{
   ## All examples use a temporary directory
   options("taxadb_default_provider"="itis_test")
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
  }

get_ids("Midas bicolor")
get_ids(c("Midas bicolor", "Homo sapiens"), format = "prefix")
get_ids("Midas bicolor", format = "uri")

}

  \dontshow{
   Sys.unsetenv("TAXADB_HOME")
  }

}
\seealso{
filter_name

Other get: 
\code{\link{get_names}()}
}
\concept{get}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_names.R
\name{get_names}
\alias{get_names}
\title{get_names}
\usage{
get_names(
  id,
  db = getOption("taxadb_default_provider", "itis"),
  version = latest_version(),
  format = c("guess", "prefix", "bare", "uri"),
  taxadb_db = td_connect()
)
}
\arguments{
\item{id}{a list of taxonomic identifiers.}

\item{db}{abbreviation code for the provider.  See details.}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  see \verb{[avialable_releases()]} for details.}

\item{format}{Format for the returned identifier, one of
\itemize{
\item \code{prefix} (e.g. \code{NCBI:9606}, the default), or
\item \code{bare} (e.g. \code{9606}, used in \code{taxize::get_ids()}),
\item \code{uri} (e.g.
\verb{http://ncbi.nlm.nih.gov/taxonomy/9606}).
}}

\item{taxadb_db}{Connection to from \verb{[td_connect()]}.}
}
\value{
a vector of names, of the same length as the input ids. Any
unmatched IDs will return as \link{NA}s.
}
\description{
Translate identifiers into scientific names
}
\details{
Like all taxadb functions, this function will run
fastest if a local copy of the provider is installed in advance
using \verb{[td_create()]}.
}
\examples{
\donttest{

\dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME = file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }

get_names(c("ITIS:1025094", "ITIS:1025103"), format = "prefix")

}

}
\seealso{
Other get: 
\code{\link{get_ids}()}
}
\concept{get}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/td_connect.R
\name{td_connect}
\alias{td_connect}
\title{Connect to the taxadb database}
\usage{
td_connect(
  dbdir = taxadb_dir(),
  driver = Sys.getenv("TAXADB_DRIVER", "duckdb"),
  read_only = FALSE
)
}
\arguments{
\item{dbdir}{Path to the database.}

\item{driver}{Default driver, one of "duckdb", "MonetDBLite", "RSQLite".
\code{taxadb} will select the first one of those it finds available if a
driver is not set. This fallback can be overwritten either by explicit
argument or by setting the environmental variable \code{TAXADB_DRIVER}.}

\item{read_only}{logical, should the database be opened read_only? Prevents
importing but will allow concurrent access from multiple sessions.}
}
\value{
Returns a DBI \code{connection} to the default duckdb database
}
\description{
Connect to the taxadb database
}
\details{
This function provides a default database connection for
\code{taxadb}. Note that you can use \code{taxadb} with any DBI-compatible database
connection  by passing the connection object directly to \code{taxadb}
functions using the \code{db} argument. \code{td_connect()} exists only to provide
reasonable automatic defaults based on what is available on your system.

\code{duckdb} or \code{MonetDBLite} will give the best performance, and regular users
\code{taxadb} will work with the built-in \code{RSQlite}, and with other database
connections such as Postgres or MariaDB, but queries (filtering joins)
will be much slower on these non-columnar databases.

For performance reasons, this function will also cache and restore the
existing database connection, making repeated calls to \code{td_connect()} much
faster and more failsafe than repeated calls to \link[DBI:dbConnect]{DBI::dbConnect}
}
\examples{
\donttest{
## OPTIONAL: you can first set an alternative home location,
## such as a temporary directory:
Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))

## Connect to the database:
db <- td_connect()

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_name.R
\name{filter_name}
\alias{filter_name}
\title{Look up taxonomic information by scientific name}
\usage{
filter_name(
  name,
  provider = getOption("taxadb_default_provider", "itis"),
  version = latest_version(),
  collect = TRUE,
  ignore_case = TRUE,
  db = td_connect()
)
}
\arguments{
\item{name}{a character vector of scientific names, e.g. "Homo sapiens"}

\item{provider}{from which provider should the hierarchy be returned?
Default is 'itis', which can also be configured using \verb{options(default_taxadb_provider=...")}.
See \verb{[td_create]} for a list of recognized providers.}

\item{version}{Which version of the taxadb provider database should we use?
defaults to latest.  See \link{tl_import} for details.}

\item{collect}{logical, default \code{TRUE}. Should we return an in-memory
data.frame (default, usually the most convenient), or a reference to
lazy-eval table on disk (useful for very large tables on which we may
first perform subsequent filtering operations.)}

\item{ignore_case}{should we ignore case (capitalization) in matching names?
default is \code{TRUE}.}

\item{db}{a connection to the taxadb database. See details.}
}
\value{
a data.frame in the Darwin Core tabular format containing the
matching taxonomic entities.
}
\description{
Look up taxonomic information by scientific name
}
\details{
Most but not all authorities can match against both species level and
higher-level (or lower, e.g. subspecies or variety) taxonomic names.
The rank level is indicated by \code{taxonRank} column.

Most authorities include both known synonyms and accepted names in the
\code{scientificName} column, (with the status indicated by \code{taxonomicStatus}).
This is convenient, as users will typically not know if the names they
have are synonyms or accepted names, but will want to get the match to the
accepted name and accepted ID in either case.
}
\examples{
\donttest{
  \dontshow{
   ## All examples use a temporary directory
   Sys.setenv(TAXADB_HOME=file.path(tempdir(), "taxadb"))
   options("taxadb_default_provider"="itis_test")
  }

sp <- c("Trochalopteron henrici gucenense",
        "Trochalopteron elliotii")
filter_name(sp)

}

}
\seealso{
Other filter_by: 
\code{\link{filter_by}()},
\code{\link{filter_common}()},
\code{\link{filter_id}()},
\code{\link{filter_rank}()}
}
\concept{filter_by}
