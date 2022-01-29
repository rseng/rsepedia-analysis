
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
