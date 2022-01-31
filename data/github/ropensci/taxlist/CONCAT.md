
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- Use snippet 'render_markdown' for it -->

# taxlist <img src="man/figures/taxlist_logo.png" height="150" align="right" />

<!-- Budges -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/taxlist)](https://cran.r-project.org/package=taxlist)
[![](https://badges.ropensci.org/233_status.svg)](https://github.com/ropensci/software-review/issues/233)
[![Rdoc](http://www.rdocumentation.org/badges/version/taxlist)](http://www.rdocumentation.org/packages/taxlist)
[![DOI](https://zenodo.org/badge/54913161.svg)](https://zenodo.org/badge/latestdoi/54913161)
<br>
[![R-CMD-check](https://github.com/ropensci/taxlist/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/taxlist/actions)
[![codecov](https://codecov.io/gh/ropensci/taxlist/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/taxlist)
<br>
[![CRAN\_downloads](http://cranlogs.r-pkg.org/badges/taxlist)](https://cran.r-project.org/package=taxlist)
[![total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/taxlist)](https://cran.r-project.org/package=taxlist)

<!-- [![DOI](https://zenodo.org/badge/54913161.svg)](https://zenodo.org/badge/latestdoi/54913161) -->

## Introduction

`taxlist` is a package designed to handle and assess taxonomic lists in
**R**, providing an object class and functions in `S4` language. The
homonymous object class `taxlist` was originally designed as a module
for taxa recorded in vegetation-plot observations (see
[`vegtable`](https://github.com/ropensci/vegtable)), but became as an
independent object with the ability of contain not only lists of species
but also synonymy, hierarchical taxonomy, and functional traits
(attributes of taxa).

The main aim of this package is to keep consistence in taxonomic lists
(a set of rules are checked by the function `validObject()`), to enable
the re-arrangement of such data, and to statistically assess functional
traits and other attributes, for instance taxonomy itself (function
`tax2traits()` set taxonomic information as trait).

While this package only includes a function for the import of taxonomic
lists from [Turboveg](https://www.synbiosys.alterra.nl/turboveg/),
almost any data source can be structured as `taxlist` object, so far the
information is imported into data frames in an R session and the
consistency rules are respected (validity).

The use of `taxlist` is recommended for people cleaning raw data before
importing it to relational databases, either in the context of taxonomic
work or biodiversity assessments. The other way around, people having
relational databases or clean and structured taxonomic lists may use
`taxlist` as recipient of this information in R sessions in order to
carry out further statistical assessments. Finally, the function
`print_name()` makes `taxlist` suitable for its implementation in
interactive documents using `rmarkdonw` and `knitr` (e.g. reports,
manuscripts and check-lists).

The structure of `taxlist` objects is inspired on the structure of data
handled by [Turboveg](https://www.synbiosys.alterra.nl/turboveg/) and
relational databases.

![](man/figures/taxlist_model.png)<br/> **Figure:** Relational model for
`taxlist` objects (see [Alvarez & Luebert
2018](https://doi.org/10.3897/BDJ.6.e23635)).

## Installing taxlist

This package is available from the Comprehensive R Archive Network
(**CRAN**) and can be directly installed within an R-session:

``` r
install.packages("taxlist", dependencies=TRUE)
```

Alternatively, the current development version is available from
[GitHub](https://github.com/ropensci/taxlist) and can be installed using
the package `devtools`:

``` r
library(devtools)
install_github("ropensci/taxlist", build_vignette=TRUE)
```

A vignette is installed with this package introducing to the work with
`taxlist` and can be accessed by following command in your R-session:

``` r
vignette("taxlist-intro")
```

## Building taxlist Objects

Objects can be built step-by-step as in the following example. For it,
we will use as reference the “Ferns of Chile” (original in Spanish:
“Helechos de Chile”) by **Gunkel (1984)**. We will create an empty
`taxlist` object using the function `new()`:

``` r
library(taxlist)

Fern <- new("taxlist")
Fern
#> object size: 5.1 Kb 
#> validation of 'taxlist' object: TRUE 
#> 
#> number of taxon usage names: 0 
#> number of taxon concepts: 0 
#> trait entries: 0 
#> number of trait variables: 0 
#> taxon views: 0
```

Then we have to set the respective taxonomic ranks. In such case, the
levels have to be provided from the lowest to highest hierarchical
level:

``` r
levels(Fern) <- c("variety","species","genus")
```

For convenience, we start inserting taxa with their respective names in
a top-down direction. We will use the function `add_concept()` to add a
new taxon. Note that the arguments `TaxonName`, `AuthorName`, and
`Level` are used to provide the name of the taxon, the authority of the
name and the taxonomic rank, respectively.

``` r
Fern <- add_concept(Fern, TaxonName="Asplenium", AuthorName="L.", Level="genus")
summary(Fern, "all")
#> ------------------------------ 
#> concept ID: 1 
#> view ID: none 
#> level: genus 
#> parent: none 
#> 
#> # accepted name: 
#> 1 Asplenium L. 
#> ------------------------------
```

As you see, the inserted genus got the concept ID **1** (see
`TaxonConceptID` in the previous figure). To insert a species of this
genus, we use again the function `add_concept()`, but this time we will
also provide the ID of the parent taxon with the argument `Parent`.

``` r
Fern <- add_concept(Fern, TaxonName="Asplenium obliquum", AuthorName="Forster",
    Level="species", Parent=1)
summary(Fern, "Asplenium obliquum")
#> ------------------------------ 
#> concept ID: 2 
#> view ID: none 
#> level: species 
#> parent: 1 Asplenium L. 
#> 
#> # accepted name: 
#> 2 Asplenium obliquum Forster 
#> ------------------------------
```

In the same way, we can add now two varieties of the inserted species:

``` r
Fern <- add_concept(Fern,
    TaxonName=c("Asplenium obliquum var. sphenoides",
        "Asplenium obliquum var. chondrophyllum"),
    AuthorName=c("(Kunze) Espinosa",
        "(Bertero apud Colla) C. Christense & C. Skottsberg"),
    Level="variety", Parent=c(2,2))
```

You may have realized that the function `summary()` is applied to
provide on the one side a display of meta-information for the whole
`taxlist` object, and on the other side to show a detail of the taxa
included in the object. In the later case adding the keyword `"all"` as
second argument, the summary will show a detailed information for every
taxon included in the object.

``` r
Fern
#> object size: 6.2 Kb 
#> validation of 'taxlist' object: TRUE 
#> 
#> number of taxon usage names: 4 
#> number of taxon concepts: 4 
#> trait entries: 0 
#> number of trait variables: 0 
#> taxon views: 0 
#> 
#> concepts with parents: 3 
#> concepts with children: 2 
#> 
#> hierarchical levels: variety < species < genus 
#> number of concepts in level variety: 2
#> number of concepts in level species: 1
#> number of concepts in level genus: 1

summary(Fern, "all")
#> ------------------------------ 
#> concept ID: 1 
#> view ID: none 
#> level: genus 
#> parent: none 
#> 
#> # accepted name: 
#> 1 Asplenium L. 
#> ------------------------------ 
#> concept ID: 2 
#> view ID: none 
#> level: species 
#> parent: 1 Asplenium L. 
#> 
#> # accepted name: 
#> 2 Asplenium obliquum Forster 
#> ------------------------------ 
#> concept ID: 3 
#> view ID: none 
#> level: variety 
#> parent: 2 Asplenium obliquum Forster 
#> 
#> # accepted name: 
#> 3 Asplenium obliquum var. sphenoides (Kunze) Espinosa 
#> ------------------------------ 
#> concept ID: 4 
#> view ID: none 
#> level: variety 
#> parent: 2 Asplenium obliquum Forster 
#> 
#> # accepted name: 
#> 4 Asplenium obliquum var. chondrophyllum (Bertero apud Colla) C. Christense & C. Skottsberg 
#> ------------------------------
```

## Indented lists

A feature implemented in version 0.2.1 is the function
`indented_list()`, which provides a better display on the hierarchical
strucutre of `taxlist` objects.

``` r
indented_list(Fern)
#> Asplenium L.
#>  Asplenium obliquum Forster
#>   Asplenium obliquum var. sphenoides (Kunze) Espinosa
#>   Asplenium obliquum var. chondrophyllum (Bertero apud Colla) C. Christense & C. Skottsberg
```

## From data frame to taxlist

A more convenient way is to create an object from a data frame including
both, the taxon concepts with their accepted names and the taxonomic
ranks with parent-child relationships. In the case of the last example,
the required data frame looks like this one:

``` r
Fern_df <- data.frame(
        TaxonConceptID=1:4,
        TaxonUsageID=1:4,
        TaxonName=c("Asplenium", "Asplenium obliquum",
                "Asplenium obliquum var. sphenoides",
                "Asplenium obliquum var. chondrophyllum"),
        AuthorName=c("L.", "Forster", "(Kunze) Espinosa",
                "(Bertero apud Colla) C. Christense & C. Skottsberg"),
        Level=c("genus", "species", "variety", "variety"),
        Parent=c(NA, 1, 2, 2),
        stringsAsFactors=FALSE)
Fern_df
#>   TaxonConceptID TaxonUsageID                              TaxonName
#> 1              1            1                              Asplenium
#> 2              2            2                     Asplenium obliquum
#> 3              3            3     Asplenium obliquum var. sphenoides
#> 4              4            4 Asplenium obliquum var. chondrophyllum
#>                                           AuthorName   Level Parent
#> 1                                                 L.   genus     NA
#> 2                                            Forster species      1
#> 3                                   (Kunze) Espinosa variety      2
#> 4 (Bertero apud Colla) C. Christense & C. Skottsberg variety      2
```

This kind of tables can be written in a spreadsheet application and
imported to your R session. The two first columns correspond to the IDs
of the taxon concept and the respective accepted name. They can be
custom IDs but are restricted to integers in `taxlist`. For the use of
the function `df2taxlist()`, the two first columns are mandatory. Also
note that the column **Parent** is pointing to the concept IDs of the
respective parent taxon. To get the object, we just use the
`df2taxlist()` indicating the sequence of taxonomic ranks in the
argument `levels`.

``` r
Fern2 <- df2taxlist(Fern_df, levels=c("variety", "species", "genus"))
Fern2
#> object size: 6.2 Kb 
#> validation of 'taxlist' object: TRUE 
#> 
#> number of taxon usage names: 4 
#> number of taxon concepts: 4 
#> trait entries: 0 
#> number of trait variables: 0 
#> taxon views: 0 
#> 
#> concepts with parents: 3 
#> concepts with children: 2 
#> 
#> hierarchical levels: variety < species < genus 
#> number of concepts in level variety: 2
#> number of concepts in level species: 1
#> number of concepts in level genus: 1
```

## Similar Packages

The package `taxlist` shares similar objectives with the package
[`taxa`](https://github.com/ropensci/taxa), but uses different
approaches for object oriented programming in **R**, namely `taxlist`
applies **S4** while `taxa` uses **R6**. Additionally, `taxa` is rather
developer-oriented, while `taxlist` is rather a user-oriented package.

In following cases you may prefer to use `taxlist`:

  - When you need an automatic check on the consistency of information
    regarding taxonomic ranks and parent-child relationships (parents
    have to be of a higher rank then children), as well as
    non-duplicated combinations of names and authors. Such checks are
    done by the function `validObject()`.
  - When you foresee statistical assessments on taxonomy diversity or
    taxon properties (chorology, conservation status, functional traits,
    etc.).
  - When you seek to produce documents using **rmarkdown**, for instance
    guide books or check-lists. Also in article manuscripts taxonomic
    names referring to a taxon concept can easily get formatted by the
    function `print_name()`.
  - When importing taxonomic lists from databases stored in
    [**Turboveg 2**](http://www.synbiosys.alterra.nl/turboveg/).
  - When you seek to implement the package
    [`vegtable`](https://CRAN.R-project.org/package=vegtable) for
    handling and assessing biodiversity records, especially
    vegetation-plot data. In that case, taxonomic lists will be
    formatted by `taxlist` as a slot within a `vegtable` object.

## Rmarkdown Integration

As mentioned before, `taxlist` objects can be also used for writing
rmarkdown documents (see [this
poster](https://dx.doi.org/10.13140/RG.2.2.36713.90728)). For instance
you can insert your objects at the beginning of the document with a
hidden chunk:

```` markdown
```{r  echo=FALSE, message=FALSE, warning=FALSE}
library(taxlist)
data(Easplist)
```
````

To mention a taxon, you can write in-line codes, such as <code>\`r
print\_name(Easplist, 206)\`</code>, which will insert *Cyperus papyrus*
L. in your document (note that the number is the ID of the taxon concept
in `Easplist`). Fort a second mention of the same species, you can then
use <code>\`r print\_name(Easplist, 206, second\_mention=TRUE)\`</code>,
which will insert *C. papyrus* in your text.

## Descriptive Statistics

Information located in the slot **taxonTraits** are suitable for
statistical assessments. For instance, in the installed object
`Easplist` a column called **lf\_behn\_2018** includes a classification
of macrophytes into different life forms. To know the frequency of these
life forms in the `Easplist`, we can use the function `count_taxa()`:

``` r
# how man taxa in 'Easplist'
count_taxa(Easplist)
#> [1] 3887

# frequency of life forms
count_taxa(~ lf_behn_2018, Easplist)
#>          lf_behn_2018 taxa_count
#> 1    acropleustophyte          8
#> 2         chamaephyte         25
#> 3      climbing_plant         25
#> 4  facultative_annual         20
#> 5     obligate_annual        114
#> 6        phanerophyte         26
#> 7    pleustohelophyte          8
#> 8          reed_plant         14
#> 9       reptant_plant         19
#> 10      tussock_plant         52
```

Furthermore, taxonomic information can be also transferred to this slot
using the function `tax2traits()`. By this way we will make taxonomic
ranks suitable for frequency calculations.

``` r
Easplist <- tax2traits(Easplist, get_names=TRUE)
head(Easplist@taxonTraits)
#>   TaxonConceptID       lf_behn_2018 form variety subspecies
#> 1              7       phanerophyte <NA>    <NA>       <NA>
#> 2              9       phanerophyte <NA>    <NA>       <NA>
#> 3             18 facultative_annual <NA>    <NA>       <NA>
#> 4             20 facultative_annual <NA>    <NA>       <NA>
#> 5             21    obligate_annual <NA>    <NA>       <NA>
#> 6             22        chamaephyte <NA>    <NA>       <NA>
#>                  species complex        genus        family
#> 1        Acacia mearnsii    <NA>       Acacia   Leguminosae
#> 2     Acacia polyacantha    <NA>       Acacia   Leguminosae
#> 3     Achyranthes aspera    <NA>  Achyranthes Amaranthaceae
#> 4     Acmella caulirhiza    <NA>      Acmella    Compositae
#> 5      Acmella uliginosa    <NA>      Acmella    Compositae
#> 6 Aeschynomene schimperi    <NA> Aeschynomene   Leguminosae
```

Note that the respective parental ranks are inserted in the table
**taxonTraits**, which contains the attributes of the taxa. In the two
next command lines, we will produce a subset with only members of the
family Cyperaceae and then calculate the frequency of species per
genera.

``` r
Cype <- subset(Easplist, family == "Cyperaceae", slot="taxonTraits")
Cype_stat <- count_taxa(species ~ genus, Cype)
```

Now, we can sort them to produce a nice bar plot.

``` r
Cype_stat <- Cype_stat[order(Cype_stat$species_count, decreasing=TRUE), ]

par(las=2, mar=c(10,5,1,1))
with(Cype_stat, barplot(species_count, names.arg=genus,
                ylab="Number of Species"))
```

![](man/figures/genera_bar-1.png)<!-- -->

## Acknowledgements

The author thanks **Stephan Hennekens**, developer of
[Turboveg](http://www.synbiosys.alterra.nl/turboveg/), for his patience
and great support finding a common language between **R** and
**Turboveg**, as well as for his advices on formatting our taxonomic
list **EA-Splist**.

Also thanks to **Federico Luebert** for the fruitful discussions
regarding the terminology used in this project.
taxlist 0.2.3
=============

### New Features

* Function `match_names()` allows to sort output data frame in the
  `'character,taxlist-method'`

### Improvements

* Slot **taxonViews** allowing class `lib_df` from package `biblio`
* New style of scripts using the package `styler`

taxlist 0.2.2
=============

### Bug Fixes

* Functions `taxlist2taxmap()` and `taxmap2taxlist()` temporarily
  deprecated due to conflicts with release of [taxa v. 0.4.0](https://github.com/ropensci/taxa/releases/tag/v0.4.0)

taxlist 0.2.1
=============

### New Features

* New function `indented_list()` to print taxonomic ranks in indented lists.

### Improvements

* New argument `repaste` in function `dissect_name()` for re-pasting
  dissected names.
* Function `replace_idx()` setting by default `idx1 = x`.
* Functions `replace_idx()` and `replace_na()` setting by default `idx2 = idx1`.
* Special characters corrected in data set *Cyperus*.
* Validation allowing taxa without rank but parents.

taxlist 0.2.0
=============

### Improvements

* Several improvements to meet **ROpenSci** requirements documented [here](https://github.com/ropensci/software-review/issues/233#issuecomment-652846890).


taxlist 0.1.9
=============

### Bug Fixes

* Problems with encoding of data set `Easplist`


taxlist 0.1.8
=============

### New Features

* Function `taxlist2taxmap()` for the interaction between packages `taxlist` and `taxmap`.
* Function `taxmap2taxlist()` for the conversion of `Taxmap` objects into `taxlist` ones.

### Improvements

* Roxygenized version.
* Method `formula` for function `count_taxa()`.
* New argument `fext` in function `backup_object()` setting the extension of the backup file.

taxlist 0.1.7
=============

### New Features

* Method for character values in function `match_names()`.
* Set of functions for data manipulation, namely `replace_x()`, `replace_idx()`, `replace_na()`, and `insert_rows()`.
* Function `clean()` with new argument **times** for repeat cleaning of `taxlist` objects.

### Improvements

* Warning in function `tax2traits()` for objects without taxonomic ranks.
* Second argument in function `[` applies only to slot **taxonTraits**.
* Replacement method for functions `[` and `$` deprecated.
* Method for function `$` matches all taxon concepts when retrieving information from slot **taxonTraits**.
* Missing argument **idx2** will be set as **idx1** in functions `replace_idx()` and `replace_na()`.
* Function `replace_view()` deprecated.
* Example data set cleaned (specifically author names)

### Bug Fixes

* Function `match_names()` was not properly working for the option `accepted_only=TRUE`.
* Function `merge_taxa()` caused orphaned children of replaced taxon concepts.
* Function `clean()` not working for deleted names.


taxlist 0.1.6
=============

### New Features

* New function `count_taxa()`

### Improvements
* A new option `style="knitr"` for function `print_name()` (See [this issue](https://stackoverflow.com/questions/51092103/formatted-scientific-names-from-r-to-latex-using-sweave-or-knitr) at **Stack Overflow**).
* In function `backup_object()`, the message will be done after successful saving and not before.
* New argument `accepted_only` in function `match_names()`, for comparing strings only with accepted names.
* Error message for NA's in argument `x` at function `match_names()`

### Bugs Fixes
* Function `add_synonym()` was not properly working for incomplete entries (missing variables in the replacement values.)
* Function `load_last()` was not properly working for values of `file` without mention of subfolder.
* Function `accepted_name()` with option `show_traits=TRUE` was not displaying taxa with no entries for taxon traits.
* Prototype for object `taxlist` wrongly included a slot **hierarchy.**

taxlist 0.1.5
=============

### New Features

* A **CITATION** file is included in the installation.
* New method `replace_view`.
* New method `print_name` for formatting taxon names to italic style.
* New method `update_name`, for updating information in slot `taxonNames`.
* New method `synonyms` retrieving synonyms for indicated concepts.
* New method `delete_name` for deleting synonyms in `taxlist` objects.
* New method `basionym` for handling basionyms.

### Improvements

* Function `accepted_name` retrieves also information on `Level` (taxonomic rank) and traits (optional in argument `show_traits`).
* Function `summary` for single taxon is displaying the name of the parent taxon (accepted name) and optional a string for the taxon view.
* Function `backup_object` prints a message in the console.
* Related functions will join documentation files.
* Data set `Easplist` adapted to new state of database **SWEA-Dataveg**.
* Function `match_names` counts multiple best matchings and includes a new argument `show_concepts` for displaying the respective accepted names and taxon concept ID.

### Bugs Fixes
* Function `load_last` was not working for single files with suffix, neither for absolute path or paths with underscores.
* Function `summary` for single taxa was not displaying names that are homonyms to the accepted name.
* Re-organized documentation.

taxlist 0.1.4
=============

### New Features

* New function `load_last` to load last backup in an R-session.
* File **inst/ChangeLog** replaced by **NEWS.md**.
* New function `dissect_name` for splitting names into their parts.
* New function `match_names` matching character vectors with names of a `taxlist` object.

### Improvements

* Function `backup_object` is also working with relative paths.

### Bugs Fixes

* Function `add_view` was not adding new columns in the respective slot.
* Function `tv2taxlist` does not modify slot `taxonViews` in prototype.
* Function `load_last` was not working with values of `filename` having underscores.

taxlist 0.1.3
=============

### New Features

* New function: `add_trait`.
* New function: `tax2traits`.

### Improvements

* Argument `level` inserted in function `merge_taxa`.
* Function `clean` also set keys to class `integer`.
* Validation checks for the existence of accepted names in names list.

### Bugs Fixes

* Bug in `add_concept`: wrong assignment of `AcceptedName`.

taxlist 0.1.2
=============

### New Features

* new function `merge_taxa`.

### Improvements

* Argument `ConceptID` in `summary` (`taxlis-method`) can be a character vector matching `TaxonName`.

taxlist 0.1.1
=============

### New Features

* New vignette `taxlist-intro`.

### Improvements

* Package `vegdata` moved from Depends to Imports.
* Function `df2taxlist` adapted to species lists with duplicated names.
* Arguments `keep_parents` and `keep_children` implemented in function `subset`.

taxlist 0.1.0
=============

### New Features

* Released to **CRAN** (https://cran.r-project.org/package=taxlist).

# Contributor Covenant Code of Conduct

## Our Pledge

We as contributors and maintainers of this project, pledge to make participation
in our community a harassment-free experience for everyone, regardless of age,
body size, visible or invisible disability, ethnicity, sex characteristics,
gender identity and expression, level of experience, education,
socio-economic status, nationality, personal appearance, race, religion, or
sexual identity and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org/),
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.
## Test environments
* ubuntu (on travis-ci)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 
### TODO

- [ ] Change the default name (**TaxonConceptID**) in output of `count_taxa()`.
- [ ] Improve documentation in `print_name()` and test different options.
- [ ] Solve problems with `as.phylo()` (communication with developer of package `ape`).
- [ ] Check for function `tnrs()` and improve it.
- [ ] Make function `check_names()` more efficient for big species lists.
- [x] An argument in `match_names()` to compare only among accepted names.
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

Before you make a substantial pull request, you should always file an issue at
the packages [BugReports](https://github.com/ropensci/taxlist/issues)
and make sure someone from the team agrees that it’s a problem. If you’ve found
a bug, create an associated issue and illustrate the bug with a minimal 
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

Please note that the taxlist project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, an use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer form discussion.

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
output:
  github_document:
    html_preview: false
---

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- Use snippet 'render_markdown' for it -->

```{r,echo=FALSE}
knitr::opts_chunk$set(
  collapse=TRUE,
  comment="#>",
  fig.path="man/figures/"
)
```

# taxlist <img src="man/figures/taxlist_logo.png" height="150" align="right" />

<!-- Budges -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/taxlist)](https://cran.r-project.org/package=taxlist)
[![](https://badges.ropensci.org/233_status.svg)](https://github.com/ropensci/software-review/issues/233)
[![DOI](https://zenodo.org/badge/54913161.svg)](https://zenodo.org/badge/latestdoi/54913161)
<br>
[![R-CMD-check](https://github.com/ropensci/taxlist/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/taxlist/actions)
[![codecov](https://codecov.io/gh/ropensci/taxlist/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/taxlist)
<br>
[![CRAN_downloads](http://cranlogs.r-pkg.org/badges/taxlist)](https://cran.r-project.org/package=taxlist)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/taxlist)](https://cran.r-project.org/package=taxlist)

<!-- [![DOI](https://zenodo.org/badge/54913161.svg)](https://zenodo.org/badge/latestdoi/54913161) -->

## Introduction

`taxlist` is a package designed to handle and assess taxonomic lists in **R**,
providing an object class and functions in `S4` language.
The homonymous object class `taxlist` was originally designed as a module for
taxa recorded in vegetation-plot observations (see
[`vegtable`](https://github.com/ropensci/vegtable)), but became as an independent
object with the ability of contain not only lists of species but also synonymy,
hierarchical taxonomy, and functional traits (attributes of taxa).

The main aim of this package is to keep consistence in taxonomic lists (a set of
rules are checked by the function `validObject()`), to enable the
re-arrangement of such data, and to statistically assess functional traits and
other attributes, for instance taxonomy itself (function `tax2traits()` set
taxonomic information as trait).

While this package only includes a function for the import of taxonomic lists
from [Turboveg](https://www.synbiosys.alterra.nl/turboveg/), almost any data
source can be structured as `taxlist` object, so far the information is
imported into data frames in an R session and the consistency rules are
respected (validity).

The use of `taxlist` is recommended for people cleaning raw data before
importing it to relational databases, either in the context of taxonomic work
or biodiversity assessments.
The other way around, people having relational databases or clean and
structured taxonomic lists may use `taxlist` as recipient of this information
in R sessions in order to carry out further statistical assessments.
Finally, the function `print_name()` makes `taxlist` suitable for its
implementation in interactive documents using `rmarkdonw` and `knitr` (e.g.
reports, manuscripts and check-lists).

The structure of `taxlist` objects is inspired on the structure of data handled
by [Turboveg](https://www.synbiosys.alterra.nl/turboveg/) and relational
databases.

![](man/figures/taxlist_model.png)<br/>
**Figure:** Relational model for `taxlist` objects (see [Alvarez & Luebert
2018](https://doi.org/10.3897/BDJ.6.e23635)).

## Installing taxlist

This package is available from the Comprehensive R Archive Network (**CRAN**)
and can be directly installed within an R-session:

```{r eval=FALSE}
install.packages("taxlist", dependencies=TRUE)
```

Alternatively, the current development version is available from
[GitHub](https://github.com/ropensci/taxlist) and can be installed using the
package `devtools`:

```{r, eval=FALSE}
library(devtools)
install_github("ropensci/taxlist", build_vignette=TRUE)
```

A vignette is installed with this package introducing to the work with
`taxlist` and can be accessed by following command in your R-session:

```{r eval=FALSE}
vignette("taxlist-intro")
```

## Building taxlist Objects

Objects can be built step-by-step as in the following example. For it, we will
use as reference the "Ferns of Chile" (original in Spanish: "Helechos de Chile")
by **Gunkel (1984)**.
We will create an empty `taxlist` object using the function `new()`:

```{r}
library(taxlist)

Fern <- new("taxlist")
Fern
```

Then we have to set the respective taxonomic ranks. In such case, the levels
have to be provided from the lowest to highest hierarchical level:

```{r}
levels(Fern) <- c("variety","species","genus")
```

For convenience, we start inserting taxa with their respective names in a
top-down direction. We will use the function `add_concept()` to add a new
taxon. Note that the arguments `TaxonName`, `AuthorName`, and `Level` are
used to provide the name of the taxon, the authority of the name and the
taxonomic rank, respectively.

```{r}
Fern <- add_concept(Fern, TaxonName="Asplenium", AuthorName="L.", Level="genus")
summary(Fern, "all")
```

As you see, the inserted genus got the concept ID **1** (see `TaxonConceptID`
in the previous figure).
To insert a species of this genus, we use again the function `add_concept()`,
but this time we will also provide the ID of the parent taxon with the argument
`Parent`.

```{r}
Fern <- add_concept(Fern, TaxonName="Asplenium obliquum", AuthorName="Forster",
	Level="species", Parent=1)
summary(Fern, "Asplenium obliquum")
```

In the same way, we can add now two varieties of the inserted species:

```{r}
Fern <- add_concept(Fern,
	TaxonName=c("Asplenium obliquum var. sphenoides",
		"Asplenium obliquum var. chondrophyllum"),
	AuthorName=c("(Kunze) Espinosa",
		"(Bertero apud Colla) C. Christense & C. Skottsberg"),
	Level="variety", Parent=c(2,2))
```

You may have realized that the function `summary()` is applied to provide on
the one side a display of meta-information for the whole `taxlist` object, and
on the other side to show a detail of the taxa included in the object. In the
later case adding the keyword `"all"` as second argument, the summary will show
a detailed information for every taxon included in the object.

```{r}
Fern

summary(Fern, "all")
```

## Indented lists

A feature implemented in version 0.2.1 is the function `indented_list()`, which
provides a better display on the hierarchical strucutre of `taxlist` objects.

```{r}
indented_list(Fern)
```

## From data frame to taxlist

A more convenient way is to create an object from a data frame including both,
the taxon concepts with their accepted names and the taxonomic ranks with
parent-child relationships.
In the case of the last example, the required data frame looks like this one:

```{r}
Fern_df <- data.frame(
		TaxonConceptID=1:4,
		TaxonUsageID=1:4,
		TaxonName=c("Asplenium", "Asplenium obliquum",
				"Asplenium obliquum var. sphenoides",
				"Asplenium obliquum var. chondrophyllum"),
		AuthorName=c("L.", "Forster", "(Kunze) Espinosa",
				"(Bertero apud Colla) C. Christense & C. Skottsberg"),
		Level=c("genus", "species", "variety", "variety"),
		Parent=c(NA, 1, 2, 2),
		stringsAsFactors=FALSE)
Fern_df
```

This kind of tables can be written in a spreadsheet application and imported to
your R session.
The two first columns correspond to the IDs of the taxon concept and the
respective accepted name. They can be custom IDs but are restricted to integers
in `taxlist`.
For the use of the function `df2taxlist()`, the two first columns are
mandatory.
Also note that the column **Parent** is pointing to the concept IDs of the
respective parent taxon.
To get the object, we just use the `df2taxlist()` indicating the sequence of
taxonomic ranks in the argument `levels`.

```{r}
Fern2 <- df2taxlist(Fern_df, levels=c("variety", "species", "genus"))
Fern2
```

## Similar Packages

The package `taxlist` shares similar objectives with the package
[`taxa`](https://github.com/ropensci/taxa), but uses different approaches for
object oriented programming in **R**, namely `taxlist` applies **S4** while
`taxa` uses **R6**. Additionally, `taxa` is rather developer-oriented, while
`taxlist` is rather a user-oriented package.

In following cases you may prefer to use `taxlist`:

- When you need an automatic check on the consistency of information regarding
  taxonomic ranks and parent-child relationships (parents have to be of a
  higher rank then children), as well as non-duplicated combinations of names
  and authors. Such checks are done by the function `validObject()`.
- When you foresee statistical assessments on taxonomy diversity or taxon
  properties (chorology, conservation status, functional traits, etc.).
- When you seek to produce documents using **rmarkdown**, for instance guide
  books or check-lists. Also in article manuscripts taxonomic names referring
  to a taxon concept can easily get formatted by the function `print_name()`.
- When importing taxonomic lists from databases stored in
  [**Turboveg 2**](http://www.synbiosys.alterra.nl/turboveg/).
- When you seek to implement the package
  [`vegtable`](https://CRAN.R-project.org/package=vegtable) for handling and
  assessing biodiversity records, especially vegetation-plot data. In that
  case, taxonomic lists will be formatted by `taxlist` as a slot within a
  `vegtable` object.


## Rmarkdown Integration

As mentioned before, `taxlist` objects can be also used for writing rmarkdown
documents (see [this poster](https://dx.doi.org/10.13140/RG.2.2.36713.90728)).
For instance you can insert your objects at the beginning of the document with a
hidden chunk:

````markdown
`r ''````{r  echo=FALSE, message=FALSE, warning=FALSE}
library(taxlist)
data(Easplist)
```
````

To mention a taxon, you can write in-line codes, such as
<code>&grave;r print_name(Easplist, 206)&grave;</code>, which will insert
`r print_name(Easplist, 206)` in your document (note that the number is the ID
of the taxon concept in `Easplist`).
Fort a second mention of the same species, you can then use
<code>&grave;r print_name(Easplist, 206, second_mention=TRUE)&grave;</code>,
which will insert `r print_name(Easplist, 206, second_mention=TRUE)` in your
text.


## Descriptive Statistics

Information located in the slot **taxonTraits** are suitable for statistical
assessments.
For instance, in the installed object `Easplist` a column called
**lf_behn_2018** includes a classification of macrophytes into different life
forms.
To know the frequency of these life forms in the `Easplist`, we can use the
function `count_taxa()`:

```{r}
# how man taxa in 'Easplist'
count_taxa(Easplist)

# frequency of life forms
count_taxa(~ lf_behn_2018, Easplist)
```

Furthermore, taxonomic information can be also transferred to this slot using
the function `tax2traits()`.
By this way we will make taxonomic ranks suitable for frequency calculations.

```{r}
Easplist <- tax2traits(Easplist, get_names=TRUE)
head(Easplist@taxonTraits)
```

Note that the respective parental ranks are inserted in the table
**taxonTraits**, which contains the attributes of the taxa.
In the two next command lines, we will produce a subset with only members of
the family Cyperaceae and then calculate the frequency of species per genera.

```{r}
Cype <- subset(Easplist, family == "Cyperaceae", slot="taxonTraits")
Cype_stat <- count_taxa(species ~ genus, Cype)
```

Now, we can sort them to produce a nice bar plot.

```{r genera_bar}
Cype_stat <- Cype_stat[order(Cype_stat$species_count, decreasing=TRUE), ]

par(las=2, mar=c(10,5,1,1))
with(Cype_stat, barplot(species_count, names.arg=genus,
				ylab="Number of Species"))
```


## Acknowledgements

The author thanks **Stephan Hennekens**, developer of
[Turboveg](http://www.synbiosys.alterra.nl/turboveg/), for his patience and
great support finding a common language between **R** and **Turboveg**, as well
as for his advices on formatting our taxonomic list **EA-Splist**.

Also thanks to **Federico Luebert** for the fruitful discussions regarding the
terminology used in this project.
---
title: Structure of taxlist objects
author: Miguel Alvarez
output: pdf_document
header-includes:
  - \usepackage{fontawesome}
  - \usepackage{tikz}
  - \usetikzlibrary{shapes}
  - \tikzset{font={\fontsize{10pt}{12}\selectfont}}
  #- \usepackage[usenames,dvipsnames]{xcolor} % colors
  - \definecolor{Grey1}{RGB}{191,191,191}
  - \definecolor{Grey2}{RGB}{222,222,222}
---

\tikzstyle{table}=[rectangle,draw=black,rounded corners,anchor=north west,
text width=3.5cm,rectangle split,rectangle split parts=2,
rectangle split part fill={Grey1,Grey2}]

\begin{tikzpicture}
    % Rectangles
    \node[table](taxonNames){
        \textbf{taxonNames}
        \nodepart{second}TaxonConceptID
        \newline TaxonUsageID \faKey
        \newline TaxonName
        \newline AuthorName
        \newline \ldots   
    };
    \node[table](taxonRelations) at ([xshift=1.5cm]taxonNames.north east){
        \textbf{taxonRelations}
        \nodepart{second}TaxonConceptID \faKey
        \newline AcceptedName
        \newline Basionym
        \newline Parent
        \newline Level
        \newline ViewID
        \newline \ldots
    };
    \node[table](taxonTraits) at ([xshift=1.5cm]taxonRelations.north east){
        \textbf{taxonTraits}
        \nodepart{second}TaxonConceptID \faKey
        \newline \ldots
    };
    \node[table](taxonViews) at ([yshift=-1cm]taxonTraits.south west){
        \textbf{taxonViews}
        \nodepart{second}ViewID \faKey
        \newline \ldots
    };
    % Arrows
    \draw([yshift=-0.7cm]taxonNames.north east)--
    ([yshift=-0.7cm]taxonRelations.north west);
    \draw([yshift=-1.2cm]taxonNames.north east)--
    ([yshift=-1.2cm]taxonRelations.north west);
    \draw([yshift=-0.7cm]taxonRelations.north east)--
    ([yshift=-0.7cm]taxonTraits.north west);
    % Auxiliar nodes and break line
    \node(ctrl5) at ([xshift=0.12cm,yshift=-1.4cm]taxonNames.north east);
    \node(ctrl6) at ([xshift=0.9cm,yshift=-1.4cm]taxonNames.north east);
    \node(ctrl7) at ([xshift=0.9cm,yshift=-1.72cm]taxonNames.north east);
    \node(ctrl8) at ([xshift=1.38cm,yshift=-1.72cm]taxonNames.north east);
        \draw(ctrl5.north west)--(ctrl6.north west)--(ctrl7.north west)--
        (ctrl8.north east);
    \node(ctrl9) at ([xshift=0.12cm,yshift=-0.9cm]taxonRelations.north east);
    \node(ctrl10) at ([xshift=0.5cm,yshift=-0.9cm]taxonRelations.north east);
    \node(ctrl11) at ([xshift=0.5cm,yshift=-2.2cm]taxonRelations.north east);
    \node(ctrl12) at ([xshift=0.12cm,yshift=-2.2cm]taxonRelations.north east);
        \draw(ctrl9.north west)--(ctrl10.north west)--(ctrl11.north west)--
        (ctrl12.north west);  
    \node(ctrl1) at ([xshift=0.12cm,yshift=-3.04cm]taxonRelations.north east);
    \node(ctrl2) at ([xshift=0.7cm,yshift=-3.04cm]taxonRelations.north east);
    \node(ctrl3) at ([xshift=0.7cm,yshift=-3.32cm]taxonRelations.north east);
    \node(ctrl4) at ([xshift=1.38cm,yshift=-3.32cm]taxonRelations.north east);
        \draw(ctrl1.north west)--(ctrl2.north west)--(ctrl3.north west)--
        (ctrl4.north east);
    % Names of arrows
    \node at([xshift=0.75cm,yshift=-0.5cm]taxonNames.north east){n:1};
    \node at([xshift=0.75cm,yshift=-0.95cm]taxonNames.north east){1:1};
    \node at([xshift=0.75cm,yshift=-1.9cm]taxonNames.north east){1:1};
    \node at([xshift=0.75cm,yshift=-0.5cm]taxonRelations.north east){1:1};
    \node at([xshift=0.75cm,yshift=-1.4cm]taxonRelations.north east){1:n};
    \node at([xshift=0.9cm,yshift=-2.9cm]taxonRelations.north east){n:1};
\end{tikzpicture}
---
title: "Applying taxlist to species lists on diversity records"
author: "Miguel Alvarez"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Applying taxlist to species lists on diversity records}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# 1. Getting started

The package `taxlist` aims to implement an object class and functions (methods)
for handling taxonomic data in **R**.
The homonymous object class `taxlist` can be further linked to biodiversity
records (e.g. for observations in vegetation plots).

The `taxlist` package is developed on the repository **GitHub**
([https://github.com/ropensci/taxlist](https://github.com/ropensci/taxlist)) and can
be installed in your R-session using the package `devtools`:

```{r install_github, eval=FALSE}
library(devtools)
install_github("ropensci/taxlist", build_vignettes=TRUE)
```

Since this package is already available in the Comprehensive R Archive Network
(CRAN), it is also possible to install it using the function
`install.packages`:

```{r install_cran, eval=FALSE}
install.packages("taxlist", dependencies=TRUE)
```

Of course, you have to load `taxlist` into your R-session.

```{r load_taxlist, message=FALSE}
library(taxlist)
```

For accessing to this vignette, use following command:

```{r call_vignette, eval=FALSE}
vignette("taxlist-intro")
```

# 2. Extracting a species list from a vegetation table

## 2.1 Example data

One of the main tasks of `taxlist` is to structure taxonomic information for a
further linkage to biodiversity records.
This structure have to be on the one side
consistent with taxonomic issues (e.g. synonyms, hierarchies, etc.), on the
other side have to be flexible for containing different depth of information
availability (from plain species lists to hierarchical structures).

In this guide, we will work with a species list from phytosociological
relev&eacute;s collected at the borderline between the
**Democratic Republic of the Congo** and **Rwanda** (Mullenders 1953
*Vegetatio* 4(2): 73--83).

![](Siegesbeckia_orientalis.png)

The digitized data can be loaded by following command:

```{r load_example_table}
load(file.path(path.package("taxlist"), "Cross.rda"))
```

The data is formatted as `data.frame` in **R**, including the names of the
species in the first column:

```{r head_example}
head(Cross[ ,1:8])
```

## 2.2 From plain list to taxlist

As already mentioned, the first column in the cross table contains the names
of the species occurring in the observed plots.
Thus, we can use this character vector to construct a `taxlist` object.
This can be achieved through the function `df2taxlist`.

```{r character2taxlist}
Splist <- Cross[ ,"TaxonName"]
Splist <- df2taxlist(Splist)
summary(Splist)
```

Note that the function `summary` provides a quick overview in the content of
the resulting object.
This function can be also applied to a specific taxon:

```{r summary_character}
summary(Splist, "Erigeron floribundus")
```

# 3. Built-in data set

## 3.1 Easplist

The installation of `taxlist` includes the data `Easplist`, which is formatted
as a `taxlist` object.
This data is a subset of the species list used by the database **SWEA-Dataveg**
([GIVD ID AF-006](http://www.givd.info/ID/AF-00-006 "SWEA-Dataveg")):


```{r load_easplist}
data(Easplist)
Easplist
```


## 3.2 Access to slots

The common ways to access to the content of slots in `S4` objects are either
using the function `slot(object, name)` or the symbol `@` (i.e. `object@name`).
Additional functions, which are specific for `taxlist` objects are
`taxon_names`, `taxon_relations`, `taxon_traits` and `taxon_views` (see the help
documentation).

Additionally, it is possible to use the methods `$` and `[` , the first for
access to information in the slot `taxonTraits`, while the second can be also
used for other slots in the object.

```{r summary_life_forms}
summary(as.factor(Easplist$lf_behn_2018))
```


## 3.3 Subsets

Methods for the function `subset` are also implemented in this package.
Such subsets usually apply pattern matching (for character vectors) or logical
operations and are analogous to query building in relational databases.
The `subset` method can be apply to any slot by setting the value of the
argument `slot`.


```{r papyrus_otp1, results="hide"}
Papyrus <- subset(x=Easplist, subset=grepl("papyrus", TaxonName), slot="names")
summary(Papyrus, "all")
```

Or the very same results:

```{r papyrus_opt2, results="hide"}
Papyrus <- subset(x=Easplist, subset=TaxonConceptID == 206, slot="relations")
summary(Papyrus, "all")
```

Similarly, you can look for a specific name.

```{r phragmites, results="hide"}
Phraaus <- subset(x=Easplist,
		subset=charmatch("Phragmites australis", TaxonName), slot="names")
summary(Phraaus, "all")
```


## 3.4 Hierarchical structure

Objects belonging to the class `taxlist` can optionally content parent-child
relationships and taxonomic levels.
Such information is also included in the data `Easplist`, as shown in the
summary output.

```{r summary_again}
Easplist
```

Note that such information can get lost once `subset()` has been applied, since
the respective parents or children from the original data set are not anymore in
the subset.
May you like to recover parents and children, you can use the functions
`get_parents` or `get_children`, respectively.

```{r recover_parents}
summary(Papyrus, "all")
Papyrus <- get_parents(Easplist, Papyrus)
summary(Papyrus, "all")
```

# 4. Applying taxlist to syntaxonomic schemes

## 4.1 Example of a phytosociological classification

To illustrate the flexibility of the `taxlist` objects, the next example will
handle a syntaxonomical scheme.
As example it will be used a scheme proposed by the author for
aquatic and semi-aquatic vegetation in Tanzania (Alvarez 2017 *Phytocoenologia*
in review).
The scheme includes 10 associations classified into 4 classes:

![](wetlands_syntax.png)

## 4.2 Building the taxlist object

The content for the taxonomic list is included in a data frame and can be
downloaded by following command:

```{r load_syntax}
load(file.path(path.package("taxlist"), "wetlands_syntax.rda"))
```

The data frame `Concepts` contains the list of syntaxon names that are
considered as accepted in the previous scheme.
This list will be used to insert the new concepts in the `taxlist` object.

```{r prototype}
head(Concepts)

Syntax <- new("taxlist")

levels(Syntax) <- c("association","alliance","order","class")

taxon_views(Syntax) <- data.frame(ViewID=1, Secundum="Alvarez (2017)",
		Author="Alvarez M", Year=2017,
        Title="Classification of aquatic and semi-aquatic vegetation in East Africa",
        stringsAsFactors=FALSE)

Syntax <- add_concept(Syntax, TaxonName=Concepts$TaxonName,
		AuthorName=Concepts$AuthorName, Parent=Concepts$Parent,
		Level=Concepts$Level, ViewID=rep(1, nrow(Concepts)))

Syntax
```

Note that the function `new` created an empty object (prototype), while
`levels` insert the custom levels (syntaxonomical hierarchies).
For the later function, the levels have to be inserted from the lower to the
higher ranks.
Furthermore the reference defining the concepts included in the syntaxonomic
scheme was inserted in the object using the function `taxon_views` and finally
the concepts were inserted by the function `add_concept`.

The next step will be inserting those names that are considered as synonyms for
the respective syntaxa.
Synonyms are included in the data frame `Synonyms`.

```{r adding_synonyms}
head(Synonyms)
Syntax <- add_synonym(Syntax, ConceptID=Synonyms$TaxonConceptID,
		TaxonName=Synonyms$TaxonName, AuthorName=Synonyms$AuthorName)
```

Finally, the codes provided for the associations will be inserted as traits
properties) of them in the slot `taxonTraits`.

```{r adding_traits}
head(Codes)
taxon_traits(Syntax) <- Codes
Syntax
```

For instance, you may like to get the parental chain from an association (e.g.
for *Nymphaeetum loti*).

```{r get_nymplot}
Nymplot <- subset(Syntax, charmatch("Nymphaeetum", TaxonName), slot="names")
summary(Nymplot, "all")
```

Note that there is the logical arguments `keep_parents` and `keep_children` to
preserve hierarchical information in the subset:

```{r get_nymplot_2}
Nymplot <- subset(Syntax, charmatch("Nymphaeetum", TaxonName), slot="names",
	keep_parents=TRUE)
summary(Nymplot, "all")
```

By using the function `subset` we just created a new object containing only the
association *Nymphaeetum loti* and its parental chain.
This subset was then used to extract the parental chain from `Syntax`.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary}
\alias{summary}
\alias{summary,taxlist-method}
\alias{show,taxlist-method}
\alias{print,taxlist-method}
\alias{print}
\title{Print overviews for taxlist Objects and their content}
\usage{
\S4method{summary}{taxlist}(
  object,
  ConceptID,
  units = "Kb",
  check_validity = TRUE,
  display = "both",
  maxsum = 5,
  secundum = NULL,
  ...
)

\S4method{show}{taxlist}(object)

\S4method{print}{taxlist}(x, ...)
}
\arguments{
\item{object, x}{A \linkS4class{taxlist} object.}

\item{ConceptID}{IDs of concepts to be displayed in the summary.}

\item{units}{Character value indicating the units shown in the object's
allocated space.}

\item{check_validity}{Logical value indicating whether the validity of
\code{object} should be checked or not.}

\item{display}{Character value indicating the field to be displayed (see
details).}

\item{maxsum}{Integer indicating the maximum number of displayed taxa.}

\item{secundum}{A character value indicating the column from slot\code{taxonViews}
to be displayed in the summary.}

\item{...}{Further arguments passed to or from another methods.}
}
\description{
A method to display either an overview of the content of
\linkS4class{taxlist} objects or an overview of selected taxa.
}
\details{
A general overview indicating number of names, concepts and taxon views
included in \linkS4class{taxlist} objects.
If argument \code{ConceptID} is a vector with concept IDs or names to be matched
by \code{\link[=grepl]{grepl()}}, then a display of all names included in each concept will be
produced.
Alternative you can use \code{taxon="all"} in order to get the listing of names
for all concepts included in the object (truncated to the input number of
\code{maxsum}).

For summaries applied to concepts, there are three alternative displays of
names using the argument \code{display}.
Use \code{display="name"} to show the value \code{TaxonName}, \code{display="author"} to
show the value \code{AuthorName} or \code{display="both"} to show both values.
Such values are taken from slot \code{taxonNames}.

For big objects it will be recommended to set \code{units="Mb"} (see also
\code{\link[=object.size]{object.size()}} for further alternatives).
}
\examples{
## summary of the object
summary(Easplist, units = "Mb")

## the same output
summary(Easplist)
show(Easplist)
print(Easplist)
Easplist

## summary for two taxa
summary(Easplist, c(51128, 51140))

## summary for a name
summary(Easplist, "Acmella")

## summary for the first 10 taxa
summary(object = Easplist, ConceptID = "all", maxsum = 10)
}
\seealso{
\linkS4class{taxlist}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_children.R
\name{get_children}
\alias{get_children}
\alias{get_children,taxlist,numeric-method}
\alias{get_children,taxlist,taxlist-method}
\alias{get_parents}
\alias{get_parents,taxlist,numeric-method}
\alias{get_parents,taxlist,taxlist-method}
\title{Retrieve children or parents of taxon concepts}
\usage{
get_children(taxlist, ConceptID, ...)

\S4method{get_children}{taxlist,numeric}(taxlist, ConceptID, ...)

\S4method{get_children}{taxlist,taxlist}(taxlist, ConceptID, ...)

get_parents(taxlist, ConceptID, ...)

\S4method{get_parents}{taxlist,numeric}(taxlist, ConceptID, ...)

\S4method{get_parents}{taxlist,taxlist}(taxlist, ConceptID, ...)
}
\arguments{
\item{taxlist}{A \linkS4class{taxlist} object.}

\item{ConceptID}{Concept IDs for selecting parents or children or a subset of
\code{taxlist}.}

\item{...}{Further arguments passed among methods.}
}
\value{
A \linkS4class{taxlist} object with a subset including
requested concepts with children or parents.
}
\description{
Retrieve all children or all parents of a queried taxon concept.
}
\details{
This function produces subsets of \linkS4class{taxlist} objects
including all children or parents of queried taxon concepts.
Multiple concepts can be queried in these function.
The argument \code{ConceptID} can be a vector of concept IDs or a subset of
the input \code{taxlist} object.
}
\examples{
## Subset with family Ebenaceae and children
Ebenaceae <- subset(Easplist, charmatch("Ebenaceae", TaxonName))
Ebenaceae <- get_children(Easplist, Ebenaceae)

summary(Ebenaceae)
summary(object = Ebenaceae, ConceptID = "all", maxsum = 100)

## Get parents of Diospyros tricolor
Diostri <- subset(
  x = Easplist, subset = TaxonConceptID == 52403,
  slot = "relations"
)
Diostri <- get_parents(Easplist, Diostri)

summary(Diostri)
summary(Diostri, "all")
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_names.R
\name{match_names}
\alias{match_names}
\alias{match_names,character,character-method}
\alias{match_names,character,taxlist-method}
\title{Search matchings between character and taxlist objects}
\usage{
match_names(x, object, ...)

\S4method{match_names}{character,character}(x, object, best = 5, clean = TRUE, decreasing = TRUE, ...)

\S4method{match_names}{character,taxlist}(
  x,
  object,
  clean = TRUE,
  output = "data.frame",
  best = 5,
  show_concepts = FALSE,
  accepted_only = FALSE,
  method = "lcs",
  sort_by,
  order_args = list(),
  ...
)
}
\arguments{
\item{x}{A character vector with names to be compared.}

\item{object}{An object of class \linkS4class{taxlist} to be compared
with.}

\item{best}{Integer value indicating how many from the best matches have to
be displayed (only working for \code{output="list"}).}

\item{clean}{Logical value, whether leading, tailing and double blanks should
be deleted from \code{x}.}

\item{decreasing}{Logical value indicating whether retrieved names should be
sorted by decreasing or increasing similarity value. In the character
method, the sorting corresponds to similarities between the queried value
and the reference vector (argument \code{object}). In the taxlist method using
\verb{'output = "data.frame'"}, the order corresponds to the similarity of the
best match (by default, no sorting is done). This argument is passed to
\code{\link[=order]{order()}}.}

\item{output}{Character value indicating the type of output. Alternative
values are "list" (taxon concepts ID's sorted by similarity for each
queried name) or "data.frame" (a table including the best match for every
queried name).}

\item{show_concepts}{Logical value, whether respective concepts should be
displayed in output or not.}

\item{accepted_only}{Logical value, whether only accepted names should be
matched or all usage names (including synonyms).}

\item{method, ...}{Further arguments passed to \code{\link[=stringsim]{stringsim()}}.}

\item{sort_by}{A character vector including the output columns used for
sorting the output table. Used only in 'character,taxlist-method'. The
function checks the presence of these values as columns in the output
data frame.}

\item{order_args}{A named list including arguments passed to \code{\link[=order]{order()}} in the
'character,taxlist-method'.}
}
\description{
Names provided in a character vector will be compared with names stored in
slot \code{taxonNames} within an object of class \linkS4class{taxlist} by
using the function \code{\link[=stringsim]{stringsim()}}.
}
\examples{
## Names to be compared
species <- c("Cperus papyrus", "Typha australis", "Luke skywalker")

## Comparing character vectors
match_names("Cyperus paper", species)

## Retrieve taxon usage names
match_names(species, Easplist)

## Display accepted names in output
match_names(x = species, object = Easplist, show_concepts = TRUE)
}
\seealso{
\code{\link[=stringsim]{stringsim()}}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dissect_name.R
\name{dissect_name}
\alias{dissect_name}
\title{Dissect Scientific Names into their Elements}
\usage{
dissect_name(x, split = " ", fixed = TRUE, repaste, ...)
}
\arguments{
\item{x}{A character vector containing taxon names.}

\item{split, fixed, ...}{Arguments passed to \code{\link[=strsplit]{strsplit()}}.}

\item{repaste}{An integer vector indicating the elements of the name selected
for the output.}
}
\value{
A character matrix with as many rows as names in the input vector.
If \code{repaste} is indicated, then the output will be a character vector.
}
\description{
Depending the degree of resolution and specific roles of nomenclature,
strings containing taxon usage names (scientific names) are constructed with
different parts.
A string with names can be consequently split into those elements, meanwhile
the number of elements may suggest the taxonomic ranks.

This function is a wrapper of \code{\link[=strsplit]{strsplit()}}, while name element can be
re-pasted if indicated in argument \code{repaste}.
}
\examples{
Easplist <- subset(x = Easplist, subset = Level == "variety", slot = "relations")
Easplist <- accepted_name(Easplist)[c(1:10), "TaxonName"]

# split name
dissect_name(Easplist)

# re-paste the two first words
dissect_name(Easplist, repaste = c(1:2))
}
\seealso{
\code{\link[=strsplit]{strsplit()}}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_name.R
\name{print_name}
\alias{print_name}
\alias{print_name,taxlist,numeric-method}
\title{Format usage names for publications}
\usage{
print_name(object, id, ...)

\S4method{print_name}{taxlist,numeric}(
  object,
  id,
  concept = TRUE,
  second_mention = FALSE,
  include_author = TRUE,
  secundum,
  style = "markdown",
  ...
)
}
\arguments{
\item{object}{An object of class \linkS4class{taxlist}.}

\item{id}{Integer containing either a concept or a name ID.}

\item{...}{Further arguments passed among methods.}

\item{concept}{Logical value, whether \code{id} corresponds to a concept ID
or a taxon usage name ID.}

\item{second_mention}{Logical value, whether the genus name should be
abbreviated or not.}

\item{include_author}{Logical value, whether authors of the name should be
mentioned or not.}

\item{secundum}{Character value indicating the column in slot \code{taxonViews}
that will be mentioned as \emph{secundum} (according to).}

\item{style}{Character value indicating the alternative format for italics
(at the moment only markdown and html implemented).}
}
\value{
A character value including format to italic font.
}
\description{
When writing on bio-diversity, usage names could be automatically inserted in
documents including the typical italic format for different elements of a
scientific name.
The function \code{print_name} can be applied either in markdown documents or
for graphics.
}
\details{
In \strong{Rmarkdown} documents use \code{*Cyperus papyrus* L.} for
inserting a formatted a species name.
}
\examples{
summary(Easplist, 363, secundum = "secundum")

## Empty plot
plot(
  x = NA, xlim = c(0, 5), ylim = c(7, 1), bty = "n", xaxt = "n", xlab = "",
  ylab = "options"
)

## Accepted name with author
text(x = 0, y = 1, labels = print_name(Easplist, 363, style = "expression"), pos = 4)

## Including taxon view
text(x = 0, y = 2, labels = print_name(Easplist, 363,
  style = "expression",
  secundum = "secundum"
), pos = 4)

## Second mention in text
text(x = 0, y = 3, labels = print_name(Easplist, 363,
  style = "expression",
  second_mention = TRUE
), pos = 4)

## Using synonym
text(x = 0, y = 4, labels = print_name(Easplist, 50037,
  style = "expression",
  concept = FALSE
), pos = 4)

## Markdown style
text(0, 5, labels = print_name(Easplist, 363, style = "markdown"), pos = 4)

## HTML style
text(0, 6, labels = print_name(Easplist, 363, style = "html"), pos = 4)

## LaTeX style for knitr
text(x = 0, y = 7, labels = print_name(Easplist, 363, style = "knitr"), pos = 4)
}
\seealso{
\code{\link[ape:mixedFontLabel]{ape::mixedFontLabel()}}.
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/df2taxlist.R
\name{df2taxlist}
\alias{df2taxlist}
\alias{df2taxlist,data.frame,logical-method}
\alias{df2taxlist,data.frame,missing-method}
\alias{df2taxlist,character,missing-method}
\title{Convert data frames into taxlist objects}
\usage{
df2taxlist(x, AcceptedName, ...)

\S4method{df2taxlist}{data.frame,logical}(x, AcceptedName, levels, ...)

\S4method{df2taxlist}{data.frame,missing}(x, AcceptedName, ...)

\S4method{df2taxlist}{character,missing}(x, AcceptedName, ...)
}
\arguments{
\item{x}{A data frame or a character vector with taxon names.}

\item{AcceptedName}{A logical vector indicating accepted names with value
\code{TRUE}.}

\item{...}{Additional vectors to be added as columns in slot\code{taxonNames}.}

\item{levels}{A vector with the names of the taxonomic ranks. This argument
is passed to \code{\link[=levels]{levels()}}.}
}
\value{
A \linkS4class{taxlist} object.
}
\description{
Taxon lists may be provided in data frame format, which will be converted to
a \linkS4class{taxlist} object.
}
\details{
In the method \code{data.frame}, the input data frame must have following columns:
\describe{
\item{TaxonUsageID}{Numeric code for the name.}
\item{TaxonConceptID}{Numeric code for the concept.}
\item{TaxonName}{Full name (usage), excluding author name.}
\item{AuthorName}{Author of the combination (taxon name).}
}

If the argument \code{AcceptedName} is missing, all names will be assumed as
accepted names.
In the alternative \code{character} method, author names have to be added as
additional vectors.

Be aware that the resulting object misses any information on taxon views,
basionyms, parent concepts, hierarchical levels and taxon traits.
All those elements can be added \emph{a posteriori} by further functions
provided in this package.
}
\examples{
## Read the table with names of Cyperus species
Cyperus <- read.csv(file = file.path(
  path.package("taxlist"), "cyperus",
  "names.csv"
), stringsAsFactors = FALSE)
head(Cyperus)

## Convert to 'taxlist' object
Cyperus <- df2taxlist(Cyperus, AcceptedName = !Cyperus$SYNONYM)
summary(Cyperus)

## Create a 'taxlist' object from character vectors
Plants <- df2taxlist(c("Triticum aestivum", "Zea mays"), AuthorName = "L.")
summary(Plants, "all")
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_traits.R
\name{taxon_traits}
\alias{taxon_traits}
\alias{taxon_traits,taxlist-method}
\alias{taxon_traits<-}
\alias{taxon_traits<-,taxlist,data.frame-method}
\alias{update_trait}
\alias{update_trait,taxlist,numeric-method}
\title{Manipulation of taxon traits in taxlist objects.}
\usage{
taxon_traits(taxlist, ...)

\S4method{taxon_traits}{taxlist}(taxlist, ...)

taxon_traits(taxlist) <- value

\S4method{taxon_traits}{taxlist,data.frame}(taxlist) <- value

update_trait(taxlist, ConceptID, ...)

\S4method{update_trait}{taxlist,numeric}(taxlist, ConceptID, ...)
}
\arguments{
\item{taxlist}{A \linkS4class{taxlist} object.}

\item{...}{Further arguments to be passed among methods.}

\item{value}{Data frame to be set as slot \code{taxonTraits}.}

\item{ConceptID}{A numeric vector with the respective taxon concept IDs.}
}
\description{
The slot \code{taxonTraits} in \linkS4class{taxlist} objects contains
attributes of taxon concepts (e.g. functional traits).
These functions are suitable for replacing, retrieving and appending trait
information in taxonomic lists.
}
\details{
Taxon traits are contained in a data frame at the slot \code{taxonTraits} in
\linkS4class{taxlist} objects.
To optimise space, this data frame contain only entries for those concepts
with information, while taxa with no information are skipped from this table.
Thus appending new variables may also have to include new rows in this slot,
which is automatically carried out by this function.

The replacement method \verb{taxon_traits<-} should be only used when
constructing \linkS4class{taxlist} objects from an empty one.
}
\examples{
head(taxon_traits(Easplist))

## Updating traits for Launaea cornuta
summary(Easplist, "Launaea cornuta")
accepted_name(taxlist = Easplist, ConceptID = 355, show_traits = TRUE)

# Update
Easplist <- update_trait(
  taxlist = Easplist, ConceptID = 355,
  lf_behn_2018 = "annual"
)
accepted_name(taxlist = Easplist, ConceptID = 355, show_traits = TRUE)
}
\seealso{
\linkS4class{taxlist}.
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{add_nacolumn}
\alias{add_nacolumn}
\title{Filling missed columns with NAs}
\usage{
add_nacolumn(x, y)
}
\arguments{
\item{x}{(\code{data.frame}) The data frame to be compared.}

\item{y}{(\code{data.frame}) The data frame used as reference.}
}
\value{
A \code{data.frame}.
}
\description{
If columns of \code{y} are missed in \code{x}, the later gets these columns filled with
\code{NA} values.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{overview_taxlist}
\alias{overview_taxlist}
\title{Function producing the overview of whole object.}
\usage{
overview_taxlist(object, units, check_validity)
}
\description{
Function producing the overview of whole object.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated-functions.R
\name{Deprecated-functions}
\alias{Deprecated-functions}
\alias{add_parent}
\alias{add_trait}
\alias{add_level}
\alias{replace_view}
\alias{taxlist2taxmap}
\alias{taxmap2taxlist}
\title{Deprecated functions}
\usage{
add_parent()

add_trait()

add_level()

replace_view()

taxlist2taxmap()

taxmap2taxlist()
}
\description{
Most of those functions have been replaced by alternative 'update' ones.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_relations.R
\name{taxon_relations}
\alias{taxon_relations}
\alias{taxon_relations,taxlist-method}
\alias{taxon_relations<-}
\alias{taxon_relations<-,taxlist,data.frame-method}
\alias{add_concept}
\alias{add_concept,taxlist,character-method}
\alias{add_concept,taxlist,taxlist-method}
\alias{update_concept}
\alias{update_concept,taxlist,numeric-method}
\title{Retrieve or replace slot taxonRelations in taxlist objects}
\usage{
taxon_relations(taxlist, ...)

\S4method{taxon_relations}{taxlist}(taxlist, ...)

taxon_relations(taxlist) <- value

\S4method{taxon_relations}{taxlist,data.frame}(taxlist) <- value

add_concept(taxlist, TaxonName, ...)

\S4method{add_concept}{taxlist,character}(taxlist, TaxonName, Level, ...)

\S4method{add_concept}{taxlist,taxlist}(taxlist, TaxonName, insert_view, ...)

update_concept(taxlist, ConceptID, ...)

\S4method{update_concept}{taxlist,numeric}(taxlist, ConceptID, ...)
}
\arguments{
\item{taxlist}{A \linkS4class{taxlist} object.}

\item{...}{Further arguments passed among methods.}

\item{value}{A \code{data.frame} object to be set as slot \code{taxonRelations}.}

\item{TaxonName}{Character vector with the accepted name for the new taxon
concepts.}

\item{Level}{Character vector indicating the level of the concept in the
list.}

\item{insert_view}{A numeric (integer) vector, indicating the views to be
inserted in \code{taxlist} or the value \code{TRUE} (see details).}

\item{ConceptID}{Concept IDs to be updated.}
}
\value{
An object of class \linkS4class{taxlist} with added names and
concepts.
}
\description{
Retrieve the content of slot \code{taxonRelations} from a
\linkS4class{taxlist} object or replace it by a new data frame.
}
\details{
The replacement method \verb{taxon_relations<-} should be only used when
constructing \linkS4class{taxlist} objects from an empty one
(prototype).

New concepts should be first added to a \linkS4class{taxlist} object
using their respective accepted names.
Synonyms can be further provided using the function \code{\link[=add_synonym]{add_synonym()}}.

Additional named vectors can be provided to be included in slot \code{taxonNames},
in the cases where those variables already exist, otherwise they will be
ignored.

It is recommended also to provide a concept view as \code{ViewID} (see
\code{\link[=taxon_views]{taxon_views()}}).
For adding a new view, use \code{\link[=add_view]{add_view()}}.
}
\examples{
## Subset for the genus Euclea and display of slot 'taxonNames'
Euclea <- subset(
  x = Easplist, subset = charmatch("Euclea", TaxonName),
  slot = "names"
)
Euclea <- get_children(Easplist, Euclea)

Euclea
taxon_relations(Euclea)

## Subset with family Ebenaceae and children
Ebenaceae <- subset(Easplist, charmatch("Ebenaceae", TaxonName))
Ebenaceae <- get_children(Easplist, Ebenaceae)

Ebenaceae
summary(object = Ebenaceae, ConceptID = "all", maxsum = 100)

## Adding a new concept
Ebenaceae <- add_concept(
  taxlist = Ebenaceae, TaxonName = "Euclea acutifolia",
  AuthorName = "E. Mey. ex A. DC.", Level = "species", Parent = 55707, ViewID = 1
)

## A summary again
Ebenaceae
summary(Ebenaceae, "all", maxsum = 100)

## Display two Typha species
summary(Easplist, c("Typha domingensis", "Typha latifolia"))

## Update a concept
summary(Easplist, "Corchorus olitorius")
Easplist <- update_concept(
  taxlist = Easplist, ConceptID = 155,
  Level = "subspecies"
)
summary(Easplist, "Corchorus olitorius")
}
\seealso{
\linkS4class{taxlist}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tnrs.R
\name{tnrs}
\alias{tnrs}
\alias{tnrs,character-method}
\alias{tnrs,taxlist-method}
\title{Taxonomic Name Resolution Service}
\usage{
tnrs(query, ...)

\S4method{tnrs}{character}(query, ...)

\S4method{tnrs}{taxlist}(query, min_score = 0.8, source = "iPlant_TNRS", ...)
}
\arguments{
\item{query}{Either a character vector or a taxlist object with names to
search.}

\item{...}{Further arguments passed to \code{\link[taxize:tnrs-defunct]{taxize::tnrs()}}.}

\item{min_score}{Minimum value of score for considering accepted names as
suggested by the output.}

\item{source}{Source database.}
}
\value{
A data frame or an object of class \linkS4class{taxlist}.
}
\description{
Methods of \code{\link[taxize:tnrs-defunct]{taxize::tnrs()}} for \linkS4class{taxlist} objects.
}
\details{
This function checks for matching of taxon names in \linkS4class{taxlist}
objects with the Taxonomic Name Resolution Service (TNRS).
Misspelled names as well as author names will be replaced in the the new
object and new accepted names will be inserted.

A method for character vectors is defined for the original function.
}
\seealso{
\code{\link[taxize:tnrs-defunct]{taxize::tnrs()}}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/levels.R
\name{levels}
\alias{levels}
\alias{levels,taxlist-method}
\alias{levels<-,taxlist-method}
\alias{levels<-}
\title{Set and retrieves hierarchical levels}
\usage{
\S4method{levels}{taxlist}(x)

\S4method{levels}{taxlist}(x) <- value
}
\arguments{
\item{x}{A \linkS4class{taxlist} object.}

\item{value}{A character vector with replacement values for levels o \code{x}.}
}
\value{
A \code{character} vector or a \linkS4class{taxlist} object with
added or modified taxonomic levels.
}
\description{
Taxonomic hierarchies can be set as levels in \linkS4class{taxlist}
objects, ordered from lower to higher levels.

Add taxonomic levels for specific taxon concepts in a
\linkS4class{taxlist} object.
Also changes in concept circumscription may implicate changes in its
taxonomic hierarchy.
}
\details{
Taxonomic levels will be handled as factors in the
\linkS4class{taxlist} objects.
Those levels are useful for creating subsets of related groups (e.g. by
functions \code{\link[=get_children]{get_children()}} or \code{\link[=get_parents]{get_parents()}}).

Levels in combination to parent-child relationships will be further used for
checking consistency of taxonomic lists.

A replacement method of the form \code{levels(x) <- value} it is also implemented.
}
\examples{
## Get levels of species list
taxlist::levels(Easplist)

## Add aggregate as new taxonomic level
levels(Easplist) <- c(
  "form", "variety", "subspecies", "species",
  "complex", "aggregate", "genus", "family"
)
summary(Easplist)
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{sort_backups}
\alias{sort_backups}
\title{Sorting files with time stamp and suffix}
\usage{
sort_backups(file, f_timestamp = "\%Y-\%m-\%d", fext = ".rda")
}
\arguments{
\item{file}{A character value indicating the root of the backup's name. It
may include also the path relative to the working directory when the
backup is stored in a sub-folder.}

\item{fext}{A character value indicating the file extension (including the
dot symbol).}
}
\value{
A data frame including the sorted names of backup files from the
oldest to the newest.
}
\description{
When backups have been saved with a time stamp in the file's name and a
suffix, in the case of more than one backup in one day.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/backup_object.R
\name{backup_object}
\alias{backup_object}
\alias{load_last}
\title{Make and load backups of R objects}
\usage{
backup_object(
  ...,
  objects = character(),
  file,
  stamp = TRUE,
  overwrite = FALSE
)

load_last(file, fext = ".rda")
}
\arguments{
\item{...}{Names of the objects to be saved (either symbols or character
strings).}

\item{objects}{A character vector indicating the names of objects to be
included in the backup file.}

\item{file}{A character value indicating the name of the backup file, without
the extension.}

\item{stamp}{A logical value indicating whether time should be stamped in the
backup name or not.}

\item{overwrite}{A logical value indicating whether existing files must be
overwritten or not.}

\item{fext}{A character value indicating the file extension (including the
dot symbol).}
}
\value{
An R image with extension \bold{*.rda}.
}
\description{
When work with data becomes risky, the best practice is to produce backup
files.
The function of \code{backup_object} is a wrapper of \code{\link[=save]{save()}}, adding a
time stamp and a suffix to the name of the resulting file (an R image file
with extension \bold{*.rda}).
The function \code{load_last} is adapted to this style, loading the newest
version to the session.
}
\details{
In both functions the argument \code{file} may include either the path
relative to the working directory or the absolute path to the file, excluding
stamps and extension.
For \code{overwrite=FALSE} (the default), a numeric suffix will be added to
the backup's name, if another backup was produced at the same day.
For \code{overwrite=TRUE} no suffix will be included in the file and existing
files will be overwritten.

The function \code{load_last()} will load the newest version among backups
stored in the same folder, provided that the backup name includes a time
stamp.
}
\examples{
\dontrun{
## A subset with Pseudognaphalium and relatives
Pseudognaphalium <- subset(x = Easplist, subset = grepl(
  "Pseudognaphalium",
  TaxonName
), slot = "names")
Pseudognaphalium <- get_parents(Easplist, Pseudognaphalium)

## Create a backup with date stamp
backup_object(Pseudognaphalium, file = "Pseudonaphalium")

## The same
backup_object(objects = "Pseudognaphalium", file = "Pseudonaphalium")

## To load the last backup into a session
load_last("Pseudognaphalium")
}

## Load pre-installed backup
load_last(file.path(path.package("taxlist"), "extdata", "Podocarpus"))
}
\seealso{
\code{\link{save}} \code{\link{load}}.
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.list.R
\name{as.list}
\alias{as.list}
\alias{S4_to_list}
\alias{as.list,taxlist-method}
\title{Coerce an S4 object to a list.}
\usage{
S4_to_list(x)

\S4method{as.list}{taxlist}(x, ...)
}
\arguments{
\item{x}{An object of class \linkS4class{taxlist} or any S4 class.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An object of class \link{list}.
}
\description{
Coercion of S4 objects to lists can be applied to explore their content,
avoiding errors caused by their validation.
}
\details{
The function \code{S4_to_list} transforms any S4 object to a list setting
slots to elements of the list and it is running internally in the method
\code{as.list} for \linkS4class{taxlist} objects.
}
\examples{
Easplist <- as.list(Easplist)
class(Easplist)
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indented_list.R
\name{indented_list}
\alias{indented_list}
\alias{indented_list,taxlist-method}
\title{Print hierarchical structure in indented lists}
\usage{
indented_list(object, ...)

\S4method{indented_list}{taxlist}(
  object,
  filter,
  keep_children = TRUE,
  keep_parents = TRUE,
  rankless_as,
  indent = " ",
  lead_br = "",
  print = TRUE,
  author = TRUE,
  level = FALSE,
  synonyms = FALSE,
  syn_encl = c("= ", ""),
  secundum,
  alphabetical = FALSE,
  ...
)
}
\arguments{
\item{object}{A \linkS4class{taxlist} object containing taxonomic concepts.}

\item{...}{Further arguments (not used yet).}

\item{filter}{A character value (optional) that will be matched with the
taxon usage names to produce a subset of 'object'. Note that this filter
will be also applied to synonyms, independent of the argument applied in
parameter 'synonyms'.}

\item{keep_children}{A logical value indicating whether children of matched
concept should be included in the result.}

\item{keep_parents}{A logical value indicating whether parents of matched
concept should be included in the result.}

\item{rankless_as}{A character vector indicating a level (taxonomic rank) to
which rankless taxa may be set before doing the list.}

\item{indent}{Symbol used for indentation. This symbol will be multiplied by
the depth of the taxonomic rank. The default is a blank space.
This can be also provided as a named vector, with a different indentation
symbol for the respective taxonomic ranks.}

\item{lead_br}{Optional line break symbol leading before the indentation.
It may be required for r-markdown documents.}

\item{print}{A logical value indicating whether the indented list should be
printed in the console or not (default = TRUE).}

\item{author}{A logical value indicating whether the author should be printed
with the name (default = TRUE).}

\item{level}{A logical value indicating whether the name of the level
(taxonomic rank) should be included before the name or not
(default = FALSE).}

\item{synonyms}{A logical value indicating whether the synonyms should be
included after accepted names or not (default = FALSE).}

\item{syn_encl}{A character vector of length 2 including the symbols used to
enclose synonyms. First value will be set before the synonyms and second
value, after the synonyms.}

\item{secundum}{A character value matching a name in slot 'taxonViews', which
will be printed as secundum (taxon view). It is not printed by default.}

\item{alphabetical}{A logical value indicating whether taxa may be sorted by
names or by IDs. The default is FALSE, thus taxa are sorted by IDs.
Note that argument TRUE may not work properly if the object contains
homonymous taxa.}
}
\value{
If 'print = TRUE', the indented list is printed in the console. The result,
which is a data frame with the elements used to format the names, can be also
assigned to an object.
}
\description{
Print taxonomic hierarchies (ranks and parent-child relationships) from
\linkS4class{taxlist} objects in an indented list.
}
\examples{
## Show taxonomy of papyrus
indented_list(Easplist, "papyrus")

## Include synonyms and taxon views
indented_list(Easplist, "papyrus",
  level = TRUE, synonyms = TRUE,
  secundum = "secundum"
)
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_taxa.R
\name{merge_taxa}
\alias{merge_taxa}
\alias{merge_taxa,taxlist,numeric,missing-method}
\alias{merge_taxa,taxlist,missing,character-method}
\alias{change_concept<-}
\alias{change_concept<-,taxlist-method}
\title{Merge concepts or move names}
\usage{
merge_taxa(object, concepts, level, ...)

\S4method{merge_taxa}{taxlist,numeric,missing}(object, concepts, print_output = FALSE, ...)

\S4method{merge_taxa}{taxlist,missing,character}(object, concepts, level, ...)

change_concept(taxlist, UsageID) <- value

\S4method{change_concept}{taxlist}(taxlist, UsageID) <- value
}
\arguments{
\item{object, taxlist}{Object of class \linkS4class{taxlist}.}

\item{concepts}{Numeric (integer) vector including taxon concepts to be
merged.}

\item{level}{Character vector indicating the lowest level for merging.}

\item{...}{Further arguments to be passed to or from other methods.}

\item{print_output}{Logical value indicating whether the merged concept
should be displayed in the console.}

\item{UsageID}{Numeric vector with taxon usage IDs to be changed from
concept.}

\item{value}{Numeric vector with taxon concept IDs to be assigned to the
names.}
}
\value{
An object of class \linkS4class{taxlist}.
}
\description{
Merge taxon concepts form a \linkS4class{taxlist} object into single ones.
}
\details{
Taxon concepts indicated in argument \code{concepts} will be merged into a
single concept.
The new concept inherits the ID and respective attributes from slots
\code{taxonRelations} and \code{taxonTraits} from the first taxon concept
indicated in argument \code{concepts}.

For convenience the resulting concept can be displayed by setting
\code{print_output=TRUE} but only when using argument \code{concepts}.

An alternative application of this function is implemented through the
argument \code{level}, where all lower rank taxa will be merged to the indicated
level or higher (if parent of merged taxa are at a higher rank).
}
\examples{
## Merge Cyperus papyrus and Cyperus dives
summary(Easplist, c(206, 197))

Easplist <- merge_taxa(
  object = Easplist, concepts = c(206, 197),
  print_output = TRUE
)

## Move the name Typha aethiopica to concept 573 (T. latifolia)
change_concept(Easplist, 53130) <- 573
summary(Easplist, c(50105, 573))
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{overview_taxon}
\alias{overview_taxon}
\title{Function producing the overview of single taxon concepts.}
\usage{
overview_taxon(object, ConceptID, display, maxsum, secundum = NULL)
}
\description{
Function producing the overview of single taxon concepts.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Extract.R
\name{Extract}
\alias{Extract}
\alias{[}
\alias{[,taxlist-method}
\alias{$,taxlist-method}
\alias{$}
\title{Extract or Replace Parts of taxlist Objects}
\usage{
\S4method{[}{taxlist}(x, i, j, drop = FALSE)

\S4method{$}{taxlist}(x, name)
}
\arguments{
\item{x}{Object of class \linkS4class{taxlist}.}

\item{i}{Integer or logical vector used as index for access to taxon
concepts, referring to the rows in slot 'taxonRelations'. These indices can
be used to produce a object with a subset of taxon concepts. It is not
recommended to use character values for this index.}

\item{j}{Integer, logical or character vector used as index for access to
variables in slot 'taxonTraits'. These indices can be used to reduce the
number of variables in the mentioned slot.}

\item{drop}{A logical value passed to \code{\link[base]{Extract}}.}

\item{name}{A symbol or character value for the method \code{$}, corresponding to
a variable either at slot 'taxonTraits' or slot 'taxonRelations'.}
}
\value{
The method \code{$} retrieves a vector, while \code{[} retrieves a subset
of the input \linkS4class{taxlist} object.
}
\description{
Quick access to slots \code{taxonTraits} and \code{taxonRelations} within
\linkS4class{taxlist} objects.
}
\examples{
## Statistics on life forms
summary(as.factor(Easplist$lf_behn_2018))

## First ten concepts in this list
summary(Easplist[1:10, ], "all")
}
\seealso{
\linkS4class{taxlist} \code{\link[taxlist]{subset}}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{add_suffix}
\alias{add_suffix}
\title{Add a suffix when strings are identical}
\usage{
add_suffix(x, y, sep = "_")
}
\arguments{
\item{x}{(\code{character} of length 1) The name to be compared.}

\item{y}{(\code{character}) Existing names, with or without suffixes.}

\item{sep}{(\code{character} of length 1) Symbol to be placed between the original
name and the suffix.}
}
\value{
A \code{character} value, either \code{x} or \code{x} with suffix if already in \code{y}.
}
\description{
An integer will be added as suffix in \code{x} if an identical value is in \code{y}.
If a value with suffix is already in \code{y}, the function will search for the
highest suffix to avoid duplicated values.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.R
\name{clean_once_taxlist}
\alias{clean_once_taxlist}
\title{One run clean function}
\usage{
clean_once_taxlist(object)
}
\description{
One run clean function
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_strings.R
\name{clean_strings}
\alias{clean_strings}
\alias{clean_strings,character-method}
\alias{clean_strings,factor-method}
\alias{clean_strings,data.frame-method}
\title{Cleaning character strings.}
\usage{
clean_strings(x, ...)

\S4method{clean_strings}{character}(x, from = "utf8", to = "utf8", ...)

\S4method{clean_strings}{factor}(x, from = "utf8", to = "utf8", ...)

\S4method{clean_strings}{data.frame}(x, from = "utf8", to = "utf8", ...)
}
\arguments{
\item{x}{Object to be cleaned.}

\item{...}{Further arguments passed among methods (not yet in use).}

\item{from, to}{Arguments passed to \code{\link[=iconv]{iconv()}}.}
}
\value{
The same as input \code{x}.
}
\description{
Multiple, leading and trailing white spaces as well as wrong encodings may
cause serious problems in information dealing with taxonomic names.
The function \code{clean_strings} get rid of them.
}
\details{
This function automatically deletes leading, trailing and multiple white
spaces, either in strings (method \code{character}), levels (method
\code{factor} or in single columns (method \code{data.frame}).
}
\examples{
library(taxlist)
clean_strings(" Cyperus    papyrus L.     ")
}
\author{
Miguel Alvarez.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxlist-package.R
\docType{package}
\name{taxlist-package}
\alias{taxlist-package}
\title{taxlist: Handling taxonomic lists.}
\description{
The class \linkS4class{taxlist} is defined in this package using the
S4 language.
The main task of \code{taxlist} objects is to keep the taxonomic coherence of
information included in taxonomic lists and to implement functions (methods)
for a proper data handling.
Objects of class \linkS4class{taxlist} can be included in further
objects, for instance in biodiversity records as done in the package
\href{https://cran.r-project.org/package=vegtable}{vegtable}.
}
\details{
The class \linkS4class{taxlist} is defined in this package using the
S4 language.
The main task of \code{taxlist} objects is to keep the taxonomic coherence of
information included in taxonomic lists and to implement functions (methods)
for a proper data handling.
Objects of class \linkS4class{taxlist} can be included in further
objects, for instance in biodiversity records as done in the package
\href{https://cran.r-project.org/package=vegtable}{vegtable}.

For a more detailed description of this package, see Alvarez & Luebert
(2018).
}
\references{
\bold{Alvarez M, Luebert F (2018).} The taxlist package: managing plant
taxonomic lists in R. \emph{Biodiversity Data Journal} 6: e23635.
\doi{10.3897/bdj.6.e23635}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_views.R
\name{taxon_views}
\alias{taxon_views}
\alias{taxon_views,taxlist-method}
\alias{taxon_views<-}
\alias{taxon_views<-,taxlist,data.frame-method}
\alias{add_view}
\alias{add_view,taxlist-method}
\title{Management of concept views in taxonomic lists.}
\usage{
taxon_views(taxlist, ...)

\S4method{taxon_views}{taxlist}(taxlist, ...)

taxon_views(taxlist) <- value

\S4method{taxon_views}{taxlist,data.frame}(taxlist) <- value

add_view(taxlist, ...)

\S4method{add_view}{taxlist}(taxlist, ...)
}
\arguments{
\item{taxlist}{A \linkS4class{taxlist} object.}

\item{...}{Further arguments to be passed among methods.}

\item{value}{An object of class \link{data.frame} containing the references
used to define the circumscription of taxon concepts included in
\code{taxlist}.}
}
\value{
An object of class \linkS4class{taxlist} with added views.
}
\description{
Retrieve or replace slot \code{taxonViews} in an object of class \linkS4class{taxlist}
}
\details{
Taxon views indicate in \linkS4class{taxlist} objects the references
determining the circumscription of the respective taxon concepts.
When adding a new concept (see \code{\link[=add_concept]{add_concept()}}), the respective
reference may not yet occur in the input \linkS4class{taxlist} object.

The term taxon view was introduced by \strong{Zhong et al. (1996)} and
corresponds to the reference used for the definition of a concept.

This function retrieves the slot \code{taxonViews} from objects of the class
\linkS4class{taxlist}.

The replacement method \verb{taxon_views<-} replaces the whole content of slot
\code{taxonViews} and it is only recommended to use when constructing a new
\linkS4class{taxlist} object from an empty prototype.
}
\examples{
## See existing views
taxon_views(Easplist)

## Add a new view
Easplist <- add_view(
  taxlist = Easplist, secundum = "Beentje et al. (1952)",
  Title = "Flora of Tropical East Africa",
  URL = "http://www.kew.org/science/directory/projects/FloraTropEAfrica.html"
)

taxon_views(Easplist)
}
\references{
\bold{Zhong Y, Jung S, Pramanik S, Beaman JH (1996).} Data model and
comparison and query methods for interacting classifications in a taxonomic
database.
\emph{Taxon} 45: 223--241.
\doi{10.1093/bioinformatics/15.2.149}
}
\seealso{
\linkS4class{taxlist}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset.R
\name{subset}
\alias{subset}
\alias{subset,taxlist-method}
\title{Subset method for taxlist objects}
\usage{
\S4method{subset}{taxlist}(
  x,
  subset,
  slot = "names",
  keep_children = FALSE,
  keep_parents = FALSE,
  ...
)
}
\arguments{
\item{x}{Object of class \linkS4class{taxlist}.}

\item{subset}{Logical vector or logical operation to apply as subset.}

\item{slot}{Character value indicating the slot to be used for the subset.}

\item{keep_children}{Logical value applied to hierarchical structures.}

\item{keep_parents}{Logical value applied to hierarchical structures.}

\item{...}{Further arguments to be passed to or from other methods.}
}
\value{
An object of class \linkS4class{taxlist}.
}
\description{
Subset of \linkS4class{taxlist} objects will be done applying either
logical operations or pattern matchings.
Subsets can be referred to information contained either in the slot
\code{taxonNames}, \code{taxonRelations} or \code{taxonTraits}.
}
\details{
The argument \code{subset} will be applied to the slot specified in argument
\code{slot}.
This argument also allows partial matchings.

Arguments \code{keep_children} and \code{keep_parents} are applied to objects
including parent-child relationships.
When those arguments are set as \code{FALSE} (the default), children or parents
of selected taxon concepts will not be included in the subset.

Be aware that \code{\link[=subset]{subset()}} won't work properly inside of function definitions.
}
\examples{
Easplist <- subset(
  x = Easplist, subset = lf_behn_2018 == "reed_plant",
  slot = "traits"
)
summary(Easplist)

summary(as.factor(Easplist$lf_behn_2018))
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accepted_name.R
\name{accepted_name}
\alias{accepted_name}
\alias{accepted_name,taxlist,numeric-method}
\alias{accepted_name,taxlist,missing-method}
\alias{accepted_name<-}
\alias{accepted_name<-,taxlist,numeric,numeric-method}
\alias{synonyms}
\alias{synonyms,taxlist,numeric-method}
\alias{synonyms,taxlist,missing-method}
\alias{basionym}
\alias{basionym,taxlist,numeric-method}
\alias{basionym,taxlist,missing-method}
\alias{basionym<-}
\alias{basionym<-,taxlist,numeric,numeric-method}
\title{Manage accepted names, synonyms and basionyms}
\usage{
accepted_name(taxlist, ConceptID, ...)

\S4method{accepted_name}{taxlist,numeric}(taxlist, ConceptID, show_traits = FALSE, ...)

\S4method{accepted_name}{taxlist,missing}(taxlist, ConceptID, ...)

accepted_name(taxlist, ConceptID) <- value

\S4method{accepted_name}{taxlist,numeric,numeric}(taxlist, ConceptID) <- value

synonyms(taxlist, ConceptID, ...)

\S4method{synonyms}{taxlist,numeric}(taxlist, ConceptID, ...)

\S4method{synonyms}{taxlist,missing}(taxlist, ConceptID, ...)

basionym(taxlist, ConceptID, ...)

\S4method{basionym}{taxlist,numeric}(taxlist, ConceptID, ...)

\S4method{basionym}{taxlist,missing}(taxlist, ConceptID, ...)

basionym(taxlist, ConceptID) <- value

\S4method{basionym}{taxlist,numeric,numeric}(taxlist, ConceptID) <- value
}
\arguments{
\item{taxlist}{An object of class \linkS4class{taxlist}.}

\item{ConceptID}{Integer containing concept IDs where to request or set names
for one category.}

\item{...}{Further arguments passed among methods.}

\item{show_traits}{Logical value, whether traits should be included in the
output of \code{accepted_name} or not.}

\item{value}{Integer containing usage IDs to be set to the respective
category in the respective taxon concept.}
}
\value{
Most of the methods return information in data frames, while
replacement methods do it as \linkS4class{taxlist} objects.
}
\description{
Taxon usage names for a taxon concept can be divided into three categories:
accepted names, basionyms and synonyms.
Each single taxon concept may at least have an accepted name, while basionym
and synonyms are optional.
The functions \code{accepted_name}, \code{basionym} and \code{synonyms}  can be used either
to display the respective usage names or to set usage names in one of those
categories.
}
\details{
The function \code{accepted_name} retrieves the accepted names for the indicated
taxon concepts or for the whole \linkS4class{taxlist} object.
By using \code{show_traits=TRUE}, the respective taxon traits will be
displayed as well, providing an overview of taxa included in the object.
The replacement method for this function will set the respective usage name
IDs as accepted names for the respective taxon concept, provided that these
names are already set as synonyms in the respective concepts.

The function \code{synonyms} is working in a similar way as \code{accepted_name}, but
this function does not include taxon traits in the output and there is no
replacing method for \code{synonyms}.
Alternatives for inserting new synonyms into a taxon concept are either
moving synonyms from other taxa by using \link{change_concept<-} or
inserting new names in the object by using \code{\link[=add_synonym]{add_synonym()}}.

The function \code{basionym} is retrieving and setting basionyms in the
respective taxon concepts similarly to \code{accepted_name}, but this function
does not retrieve any information on taxon traits, either.
}
\examples{
## Set a different accepted name for Cyclosorus interruptus
summary(Easplist, "Cyclosorus interruptus")
accepted_name(Easplist, 50074) <- 53097
summary(Easplist, 50074)

## Inserting a new name first
summary(Easplist, "Basella alba")
Easplist <- add_synonym(
  taxlist = Easplist, ConceptID = 68,
  TaxonName = "Basella cordifolia", AuthorName = "Lam."
)
summary(Easplist, 68)
accepted_name(Easplist, 68) <- 56139
summary(Easplist, 68)
}
\seealso{
\code{\link[=add_synonym]{add_synonym()}} \link{change_concept<-}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxlist-class.R
\docType{class}
\name{taxlist-class}
\alias{taxlist-class}
\title{An S4 class to represent taxonomic lists.}
\description{
Class for taxonomic lists including synonyms, hierarchical ranks,
parent-child relationships, taxon views and taxon traits.

Note that each taxon becomes an identifier, represented by the column
\strong{TaxonConceptID} in the slot \strong{taxonRelations}, analogous to a primary key
in a relational database.
This identifier is restricted to an integer in \code{taxlist} and is specific for
the object.

In the same way, each taxon usage name has an identifier in the column
\strong{TaxonUsageID}, slot \strong{taxonNames}.
The column \strong{ViewID} in slot \strong{taxonViews} is the identifier of the taxon
view.
}
\section{Slots}{

\describe{
\item{\code{taxonNames}}{(\code{data.frame}) Table of taxon usage names (accepted names
and synonyms).}

\item{\code{taxonRelations}}{(\code{data.frame}) Relations between concepts, accepted
names, basionyms, parents and hierarchical level.}

\item{\code{taxonTraits}}{Table of taxon traits.}

\item{\code{taxonViews}}{References used to determine the respective concept
circumscription.}
}}

\examples{
library(taxlist)

showClass("taxlist")

## Create an empty object
Splist <- new("taxlist")
}
\references{
\bold{Alvarez M, Luebert F (2018).} The taxlist package: managing plant
taxonomic lists in R. \emph{Biodiversity Data Journal} 6: e23635.
\doi{10.3897/bdj.6.e23635}
}
\author{
Miguel Alvarez
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.R
\name{clean}
\alias{clean}
\alias{clean,taxlist-method}
\title{Delete orphaned records}
\usage{
clean(object, ...)

\S4method{clean}{taxlist}(object, times = 2, ...)
}
\arguments{
\item{object}{A \linkS4class{taxlist} object.}

\item{...}{Further arguments passed from or to other methods.}

\item{times}{An integer indicating how many times the cleaning should be
repeated.}
}
\value{
A clean \linkS4class{taxlist} object.
}
\description{
Manipulation of slots may generate orphaned entries in
\linkS4class{taxlist} objects.
The function \code{clean} deletes such entries and restores the consistency
of the objects.
}
\details{
Cleaning of objects will follow the deletion of orphaned names, orphaned
taxon trait entries, and orphaned parent entries.
}
\examples{
## Direct manipulation of slot taxonRelations generates an invalid object
Easplist@taxonRelations <- Easplist@taxonRelations[1:5, ]
\dontrun{
summary(Easplist)
}

## Now apply cleaning
Easplist <- clean(Easplist)
summary(Easplist)
}
\author{
Miguel Alvarez.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax2traits.R
\name{tax2traits}
\alias{tax2traits}
\alias{tax2traits,taxlist-method}
\title{Set taxonomic information as taxon traits}
\usage{
tax2traits(object, ...)

\S4method{tax2traits}{taxlist}(object, get_names = FALSE, ...)
}
\arguments{
\item{object}{An object of class \linkS4class{taxlist}.}

\item{...}{Further arguments to be passed among methods.}

\item{get_names}{Logical value indicating whether taxon names should be
retrieved instead of taxon IDs.}
}
\value{
An object of class \linkS4class{taxlist} with taxonomy added
as traits.
}
\description{
Taxonomic classification can be included in \linkS4class{taxlist}
objects within the information provided at slot \code{taxonRelations}.
Nevertheless, for statistical analyses it may be more convenient to insert
such information in the slot \code{taxonTraits}.
}
\details{
This function can only be applied to objects containing parent-child
relationships and information on taxonomic levels.
}
\examples{
## Family Acanthaceae with children
Acanthaceae <- subset(
  x = Easplist, subset = TaxonName == "Acanthaceae",
  slot = "names", keep_children = TRUE
)
summary(Acanthaceae)

## Insert taxonomy to taxon traits
Acanthaceae <- tax2traits(Acanthaceae, get_names = TRUE)
head(taxon_traits(Acanthaceae))
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxon_names.R
\name{taxon_names}
\alias{taxon_names}
\alias{taxon_names,taxlist-method}
\alias{taxon_names<-}
\alias{taxon_names<-,taxlist,data.frame-method}
\alias{add_synonym}
\alias{add_synonym,taxlist-method}
\alias{update_name}
\alias{update_name,taxlist,numeric-method}
\alias{delete_name}
\alias{delete_name,taxlist,numeric-method}
\title{Handle information on taxon usage names.}
\usage{
taxon_names(taxlist, ...)

\S4method{taxon_names}{taxlist}(taxlist, ...)

taxon_names(taxlist) <- value

\S4method{taxon_names}{taxlist,data.frame}(taxlist) <- value

add_synonym(taxlist, ConceptID, ...)

\S4method{add_synonym}{taxlist}(taxlist, ConceptID, TaxonName, AuthorName, ...)

update_name(taxlist, UsageID, ...)

\S4method{update_name}{taxlist,numeric}(taxlist, UsageID, ...)

delete_name(taxlist, UsageID, ...)

\S4method{delete_name}{taxlist,numeric}(taxlist, UsageID, ...)
}
\arguments{
\item{taxlist}{A \linkS4class{taxlist} object to be modified.}

\item{...}{Further arguments passed among methods. In \code{update_name} are
vectors including the variables to be updated for the respective taxon
usage ID.}

\item{value}{A data frame used as new slot \code{taxonNames} in \code{taxlist}.}

\item{ConceptID}{Numeric vector indicating the concept ID to which the
synonyms will be added.}

\item{TaxonName, AuthorName}{Character values used for the new names
(synonyms).}

\item{UsageID}{Numeric vector indicating the taxon usage IDs to be updated.}
}
\value{
A data frame or, in the case of the replacement method, a
\linkS4class{taxlist} object with modified slot \code{taxonNames}.
}
\description{
The slot \code{taxonNames} in \linkS4class{taxlist} objects contains
taxon usage names for the respective taxon.
These functions assist on the access and modification of entries for names.
}
\details{
The replacement method \verb{taxon_names<-} is a quick alternative to include
names in empty \linkS4class{taxlist} objects.

The function \code{add_synonym()} works only for adding names to existing
taxon concepts.
For adding new taxon concepts as well you should use \code{\link[=add_concept]{add_concept()}}.
}
\examples{
## Display of slot 'taxonNames'
Euclea <- subset(
  x = Easplist, subset = charmatch("Euclea", TaxonName),
  slot = "names", keep_children = TRUE
)
Euclea
taxon_names(Euclea)

## Insert a synonym to Diospyros scabra
summary(Easplist, "Diospyros scabra")
Easplist <- add_synonym(
  taxlist = Easplist, ConceptID = 51793,
  TaxonName = "Maba scabra", AuthorName = "Chiov."
)
summary(Easplist, "Diospyros scabra")

## Delete a synonym of Launaea cornuta
summary(Easplist, "Launaea cornuta")
Easplist <- delete_name(Easplist, 53821)
summary(Easplist, "Launaea cornuta")

## Hypothetical correction in author name in Launaea cornuta
Easplist <- update_name(taxlist = Easplist, UsageID = 355, AuthorName = "L.")
summary(Easplist, "Launaea cornuta")
}
\seealso{
\linkS4class{taxlist}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tv2taxlist.R
\name{tv2taxlist}
\alias{tv2taxlist}
\title{Import species lists from Turboveg databases}
\usage{
tv2taxlist(taxlist, tv_home = tv.home())
}
\arguments{
\item{taxlist}{The name of a species list in Turboveg as character value.}

\item{tv_home}{Character value indicating the path to the main Turboveg
folder.}
}
\value{
An object of class \linkS4class{taxlist}.
}
\description{
Importing species lists from Turboveg
\url{https://www.synbiosys.alterra.nl/turboveg/} databases into an object of
class \linkS4class{taxlist}.
}
\details{
This function imports species lists using the function \code{\link[=read.dbf]{read.dbf()}}.
When available, also taxon traits will be imported into the output object
(usually the file \strong{ecodbase.dbf}).
During import of taxon traits, duplicated entries for a same concept will
be discarded as well as entries for non-existing concepts.

By default \code{tv_home} will be set by the function \code{\link[=tv.home]{tv.home()}} from the
package \link{vegdata-package}.

By default, the name of the database will be set as concept view for all
concepts included in the species list.
If this is not correct, consider setting it manually by using the functions
\code{\link[=taxon_views]{taxon_views()}} and \code{\link[=add_view]{add_view()}}.
}
\examples{
## Cyperus data set installed as Turboveg species list
Cyperus <- tv2taxlist(
  taxlist = "cyperus",
  tv_home = file.path(path.package("taxlist"), "tv_data")
)

summary(Cyperus)
}
\seealso{
\linkS4class{taxlist}
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/replace_x.R
\name{replace_x}
\alias{replace_x}
\alias{replace_idx}
\alias{replace_na}
\alias{insert_rows}
\title{Data manipulation.}
\usage{
replace_x(x, old, new)

replace_idx(x, idx1 = x, idx2 = idx1, new)

replace_na(x, idx1, idx2 = idx1, new)

insert_rows(x, y)
}
\arguments{
\item{x}{A vector to be modified. In the case of \code{insert_rows()}, \code{x} is a
data frame.}

\item{old}{A vector with values to be replaced by \code{replace_x()} in a vector.}

\item{new}{A vector containing values to be inserted, either comparing values
or using indices.}

\item{idx1, idx2}{Indices applied for value replacements to match \code{x} with
\code{new}, respectively. If \code{idx2} is not provided, it will be assumed as
equivalent to \code{idx1}.}

\item{y}{A data frame including rows (and columns) to be inserted in \code{x}.}
}
\value{
A vector or data frame with the modified values.
}
\description{
This is a series of functions designed for a fast coding of replacements
both, as internal functions and in workflows dealing with information stored
in vectors and data frames.
Such functions are especially useful when handling with functional traits
stored in \linkS4class{taxlist} objects.

\code{replace_x()} is used to exchange values in vectors.
\code{replace_idx()} changes values in vectors by matching indices or conditions.
The function \code{replace_na()} works in the same way as \code{replace_idx()} but will
only insert values in empty elements (NAs).

The function \code{insert_rows()} will add rows and columns at the same time.
This function will be used when a new table is appended to another but
sharing only part of the columns.
}
\examples{
## Replace values in vector
replace_x(x = letters, old = c("b", "p", "f"), new = c("bee", "pork", "fungus"))

## Replace values using indices
replace_idx(
  x = letters, idx1 = 1:length(letters), idx2 = c(2, 7, 17),
  new = c("second", "seventh", "seventeenth")
)

## Replace values if they are NAs
letters[2] <- NA
replace_na(
  x = letters, idx1 = 1:length(letters), idx2 = c(1:3),
  new = c("alpha", "beta", "zeta")
)

## The same applications but this time for functional traits
summary(as.factor(Easplist$lf_behn_2018))

# Merge annuals
Easplist@taxonTraits$lifeform <- replace_x(
  x = Easplist@taxonTraits$lf_behn_2018,
  old = c("obligate_annual", "facultative_annual"),
  new = c("annual", "annual")
)
summary(as.factor(Easplist$lifeform))

# The same effect
Easplist@taxonTraits$lifeform <- replace_idx(
  x = Easplist@taxonTraits$lf_behn_2018,
  idx1 = grepl("annual", Easplist@taxonTraits$lf_behn_2018),
  idx2 = TRUE,
  new = "annual"
)
summary(as.factor(Easplist$lifeform))

## Merge data frames including new columns
data(iris)
iris$Species <- paste(iris$Species)
new_iris <- data.frame(
  Species = rep("humilis", 2), Height = c(15, 20),
  stringsAsFactors = FALSE
)
insert_rows(iris, new_iris)
}
\author{
Miguel Alvarez.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_taxa.R
\name{count_taxa}
\alias{count_taxa}
\alias{count_taxa,character,missing-method}
\alias{count_taxa,factor,missing-method}
\alias{count_taxa,taxlist,missing-method}
\alias{count_taxa,formula,taxlist-method}
\title{Count taxa within a taxlist object.}
\usage{
count_taxa(object, data, ...)

\S4method{count_taxa}{character,missing}(object, na.rm = TRUE, ...)

\S4method{count_taxa}{factor,missing}(object, na.rm = TRUE, ...)

\S4method{count_taxa}{taxlist,missing}(object, level, ...)

\S4method{count_taxa}{formula,taxlist}(object, data, include_na = FALSE, suffix = "_count", ...)
}
\arguments{
\item{object}{An object containing a taxonomic list or a formula.}

\item{data}{An object of class \linkS4class{taxlist} in the \code{formula} method.}

\item{...}{further arguments passed among methods.}

\item{na.rm}{Logical value, whether NAs have to be removed from the input
vector or not.}

\item{level}{Character value indicating the taxonomic rank of counted taxa.}

\item{include_na}{Logical value indicating whether \code{NA} values in a taxon
trait should be considered for counting taxa or just ignored (only
used in \code{formula} method).}

\item{suffix}{Character value used as suffix for the counted rank in the
output data frame (only used in \code{formula} method).}
}
\value{
An integer with the number of taxa.
}
\description{
Counting number of taxa within \linkS4class{taxlist} objects or
character vectors containing taxon names.
}
\details{
This function is written by convenience in order to reduce code for counting
taxa within \linkS4class{taxlist} objects and it is just a wrapper of \code{\link[=length]{length()}}.
}
\examples{
## factor method
count_taxa(iris$Species)

## taxlist method
count_taxa(Easplist)
count_taxa(Easplist, level = "species")

## using a formula
count_taxa(~lf_behn_2018, Easplist)
}
\author{
Miguel Alvarez \email{kamapu78@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{two2one_df}
\alias{two2one_df}
\title{Inserting new rows and columns by merging two data frames}
\usage{
two2one_df(x, y)
}
\arguments{
\item{x, y}{(\code{data.frame}) The data frames to be merged.}
}
\description{
Two data frames sharing partially columns will be merged including the sum
of all variables.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Easplist-data.R
\docType{data}
\name{Easplist-data}
\alias{Easplist-data}
\alias{Easplist}
\title{List of vascular plants from East Africa}
\format{
An object of class \linkS4class{taxlist}.
}
\source{
\href{http://www.ville-ge.ch/musinfo/bd/cjb/africa/recherche.php}{African
Plant Database},
\href{http://www.givd.info/ID/AF-00-006}{SWEA-Dataveg}.
}
\usage{
Easplist
}
\description{
Example of an incomplete taxonomic list including taxa recorded in East
Africa.
}
\details{
This list is a subset of the taxonomic list implemented in the database
\href{http://www.givd.info/ID/AF-00-006}{SWEA-Dataveg}.
Since this list is being complemented regarding stored vegetation plots, it
is an incomplete list.
}
\examples{
summary(Easplist)
}
\keyword{datasets}
