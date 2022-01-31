
# EML <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Build
Status](https://travis-ci.org/ropensci/EML.svg?branch=master)](https://travis-ci.org/ropensci/EML)
[![Windows build
status](https://ci.appveyor.com/api/projects/status/u2gw24yfkvxgny96?svg=true)](https://ci.appveyor.com/project/cboettig/eml)
[![codecov](https://codecov.io/gh/ropensci/EML/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/EML)
[![CRAN
status](https://www.r-pkg.org/badges/version/EML)](https://cran.r-project.org/package=EML)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/RNeXML)
[![DOI](https://zenodo.org/badge/10894022.svg)](https://zenodo.org/badge/latestdoi/10894022)

<!-- README.md is generated from README.Rmd. Please edit that file -->

EML is a widely used metadata standard in the ecological and
environmental sciences. We strongly recommend that interested users
visit the [EML Homepage](https://eml.ecoinformatics.org/) for an
introduction and thorough documentation of the standard. Additionally,
the scientific article *[The New Bioinformatics: Integrating Ecological
Data from the Gene to the Biosphere (Jones et
al 2006)](https://doi.org/10.1146/annurev.ecolsys.37.091305.110031)*
provides an excellent introduction into the role EML plays in building
metadata-driven data repositories to address the needs of highly
heterogeneous data that cannot be easily reduced to a traditional
vertically integrated database. At this time, the `EML` R package
provides support for the serializing and parsing of all low-level EML
concepts, but still assumes some familiarity with the EML standard,
particularly for users seeking to create their own EML files. We hope to
add more higher-level functions which will make such familiarity less
essential in future development.

## Notes on the EML v2.0 Release

`EML` v2.0 is a complete re-write which aims to provide both a drop-in
replacement for the higher-level functions of the existing EML package
while also providing additional functionality. This new `EML` version
uses only simple and familiar list structures (S3 classes) instead of
the more cumbersome use of S4 found in the original `EML`. While the
higher-level functions are identical, this makes it easier to for most
users and developers to work with `eml` objects and also to write their
own functions for creating and manipulating EML objects. Under the hood,
`EML` relies on the [emld](https://github.com/ropensci/emld/) package,
which uses a Linked Data representation for EML. It is this approach
which lets us combine the simplicity of lists with the specificity
required by the XML schema.

This revision also supports the **[recently released EML 2.2.0
specification](https://eml.ecoinformatics.org/whats-new-in-eml-2-2-0.html)**.

# Creating EML

``` r
library(EML)
```

## A minimal valid EML document:

``` r
me <- list(individualName = list(givenName = "Carl", surName = "Boettiger"))
my_eml <- list(dataset = list(
              title = "A Minimal Valid EML Dataset",
              creator = me,
              contact = me)
            )


write_eml(my_eml, "ex.xml")
#> NULL
eml_validate("ex.xml")
#> [1] TRUE
#> attr(,"errors")
#> character(0)
```

## A Richer Example

Here we show the creation of a relatively complete EML document using
`EML`. This closely parallels the function calls shown in the original
EML [R-package
vignette](https://docs.ropensci.org/EML/articles/creating-EML.html).

## `set_*` methods

The original EML R package defines a set of higher-level `set_*` methods
to facilitate the creation of complex metadata structures. `EML`
provides these same methods, taking the same arguments for
`set_coverage`, `set_attributes`, `set_physical`, `set_methods` and
`set_textType`, as illustrated here:

### Coverage metadata

``` r
geographicDescription <- "Harvard Forest Greenhouse, Tom Swamp Tract (Harvard Forest)"
coverage <- 
  set_coverage(begin = '2012-06-01', end = '2013-12-31',
               sci_names = "Sarracenia purpurea",
               geographicDescription = geographicDescription,
               west = -122.44, east = -117.15, 
               north = 37.38, south = 30.00,
               altitudeMin = 160, altitudeMaximum = 330,
               altitudeUnits = "meter")
```

### Reading in text from Word and Markdown

We read in detailed methods written in a Word doc. This uses EML’s
docbook-style markup to preserve formatting of paragraphs, lists,
titles, and so forth. (This is a drop-in replacement for EML
`set_method()`)

``` r
methods_file <- system.file("examples/hf205-methods.docx", package = "EML")
methods <- set_methods(methods_file)
```

We can also read in text that uses Markdown for markup elements:

``` r
abstract_file <-  system.file("examples/hf205-abstract.md", package = "EML")
abstract <- set_TextType(abstract_file)
```

### Attribute Metadata from Tables

Attribute metadata can be verbose, and is often defined in separate
tables (e.g. separate Excel sheets or `.csv` files). Here we use
attribute metadata and factor definitions as given from `.csv` files.

``` r
attributes <- read.table(system.file("extdata/hf205_attributes.csv", package = "EML"))
factors <- read.table(system.file("extdata/hf205_factors.csv", package = "EML"))
attributeList <- 
  set_attributes(attributes, 
                 factors, 
                 col_classes = c("character", 
                                 "Date",
                                 "Date",
                                 "Date",
                                 "factor",
                                 "factor",
                                 "factor",
                                 "numeric"))
```

### Data file format

Though the `physical` metadata specifying the file format is extremely
flexible, the `set_physical` function provides defaults appropriate for
`.csv` files. DEVELOPER NOTE: ideally the `set_physical` method should
guess the appropriate metadata structure based on the file extension.

``` r
physical <- set_physical("hf205-01-TPexp1.csv")
```

## Generic construction

In the `EML` R package, objects for which there is no `set_` method are
constructed using the `new()` S4 constructor. This provided an easy way
to see the list of available slots. In `eml2`, all objects are just
lists, and so there is no need for special methods. We can create any
object directly by nesting lists with names corresponding to the EML
elements. Here we create a `keywordSet` from scratch:

``` r
keywordSet <- list(
    list(
        keywordThesaurus = "LTER controlled vocabulary",
        keyword = list("bacteria",
                    "carnivorous plants",
                    "genetics",
                    "thresholds")
        ),
    list(
        keywordThesaurus = "LTER core area",
        keyword =  list("populations", "inorganic nutrients", "disturbance")
        ),
    list(
        keywordThesaurus = "HFR default",
        keyword = list("Harvard Forest", "HFR", "LTER", "USA")
        ))
```

Of course, this assumes that we have some knowledge of what the possible
terms permitted in an EML keywordSet are\! Not so useful for novices. We
can get a preview of the elements that any object can take using the
`emld::template()` option, but this involves a two-part workflow.
Instead, `eml2` provides generic `construct` methods for all objects.

## Constructor methods

For instance, the function `eml$creator()` has function arguments
corresponding to each possible slot for a creator. This means we can
rely on *tab completion* (and/or autocomplete previews in RStudio) to
see what the possible options are. `eml$` functions exist for all
complex types. If `eml$` does not exist for an argument (e.g. there is
no `eml$givenName`), then the field takes a simple string argument.

### Creating parties (creator, contact, publisher)

``` r
aaron <- eml$creator(
  individualName = eml$individualName(
    givenName = "Aaron", 
    surName = "Ellison"),
  electronicMailAddress = "fakeaddress@email.com")
```

``` r
HF_address <- eml$address(
                  deliveryPoint = "324 North Main Street",
                  city = "Petersham",
                  administrativeArea = "MA",
                  postalCode = "01366",
                  country = "USA")
```

``` r
publisher <- eml$publisher(
                 organizationName = "Harvard Forest",
                 address = HF_address)
```

``` r
contact <- 
  list(
    individualName = aaron$individualName,
    electronicMailAddress = aaron$electronicMailAddress,
    address = HF_address,
    organizationName = "Harvard Forest",
    phone = "000-000-0000")
```

### Putting it all together

``` r
my_eml <- eml$eml(
           packageId = uuid::UUIDgenerate(),  
           system = "uuid",
           dataset = eml$dataset(
               title = "Thresholds and Tipping Points in a Sarracenia",
               creator = aaron,
               pubDate = "2012",
               intellectualRights = "http://www.lternet.edu/data/netpolicy.html.",
               abstract = abstract,
               keywordSet = keywordSet,
               coverage = coverage,
               contact = contact,
               methods = methods,
               dataTable = eml$dataTable(
                 entityName = "hf205-01-TPexp1.csv",
                 entityDescription = "tipping point experiment 1",
                 physical = physical,
                 attributeList = attributeList)
               ))
```

## Serialize and validate

We can also validate first and then serialize:

``` r
eml_validate(my_eml)
#> [1] TRUE
#> attr(,"errors")
#> character(0)
write_eml(my_eml, "eml.xml")
#> NULL
```

## Setting the version

EML will use the latest EML specification by default. To switch to a
different version, use `emld::eml_version()`

``` r
emld::eml_version("eml-2.1.1")
#> [1] "eml-2.1.1"
```

Switch back to the 2.2.0 release:

``` r
emld::eml_version("eml-2.2.0")
#> [1] "eml-2.2.0"
```
# EML 2.0.5

* migrate upstream namespace changes units::as.units -> units::as_units

# EML 2.0.4

* bugfix for CRAN testing

# EML 2.0.3

* Note recent improvements to validation have been inherited through the release of `emld` 0.5.0,
  package dependency now requires upgrading `emld` as well.
* Fixed a bug in `set_attributes` causing an error when specifying an interval `measurementScale`. (#293)
* Updated test suite to account for the switch from taxize to taxadb
* Updated test suite to match recent changes in `emld` 0.5.0 regarding unit definitions (See https://github.com/ropensci/emld/issues/56)

# EML 2.0.2

* minor bugfix to documentation
* Moves to taxadb in place of taxize for optional species classification
* Note recent improvements to validation have been inherited through the release of `emld` 0.4.0,
  package dependency now requires upgrading `emld` as well.

# EML 2.0.1

* Improve error message for get_attributes (#286)
* Avoid `ifelse()` for portability (#283)
* Avoid edge case that can create invalid EML in `set_taxonomicNames()` (#280)
* Add documentation regarding the use of dimensionless units (#276)
* Avoid test errors on systems for which pandoc cannot be installed (#290)

# EML 2.0.0

* EML 2.0.0 is a ground-up rewrite of EML 1.x package.  The primary difference
  is that EML 2.0.0 is built on S3 (list) objects instead of S4 object system.
  This makes the package interface easier to use and extend.  Under the hood, this
  approach relies on the `emld` package, which uses a JSON-LD representation of EML
  which provides a natural translation into the list-based format.  
  
  While most high level functions for creating EML have been preserved, the change to
  S3 means that this package will not be backwards-compatible with  many scripts
  which relied on the S4 system. 

* Added a `NEWS.md` file to track changes to the package.
This release fixes an issue created by the upcoming 0.5.0 release of `emld`.  Checks
will pass once that version (emld 0.5.0) is accepted to CRAN (currently in pending
with no errors, notes, or warnings, but does trigger the revdep check for this package
to fail).  Thanks.


## R CMD check results

0 errors | 0 warnings | 0 notes

The primary goal of this project is to determine
experimentally the amount of lead time required to prevent a state
change. To achieve this goal, we will (1) experimentally induce state
changes in a natural aquatic ecosystem - the Sarracenia microecosystem;
(2) use proteomic analysis to identify potential indicators of states
and state changes; and (3) test whether we can forestall state changes
by experimentally intervening in the system. This work uses state-of-the
art molecular tools to identify early warning indicators in the field
of aerobic to anaerobic state changes driven by nutrient enrichment
in an aquatic ecosystem. The study tests two general hypotheses: (1)
proteomic biomarkers can function as reliable indicators of impending
state changes and may give early warning before increasing variances
and statistical flickering of monitored variables; and (2) well-timed
intervention based on proteomic biomarkers can avert future state changes
in ecological systems.
## General Protocols

Field methods. All experiments will be carried out in the greenhouse at Harvard Forest. We have developed an instrumentation system that allows us to collect continuous dissolved [O2] measurements: dedicated micro-probes (DO-166MT; Lazar Research Laboratories: http://www.shelfscientific.com/) connected to multiplexers and data loggers (AM16/32B multiplexer, CR-1000 datalogger and control system [Campbell Scientific: http://www.cambellsci.com]). The initial ecosystem composition in all experimental plants will be standardized by seeding each pitcher with a 10-ml inoculum of liquid collected from pitchers growing at Tom Swamp.  In all experiments, prey will be supplied to pitchers as standardized aliquots of dried and finely ground bald-faced hornets (Dolichovespula maculata; Hymenoptera: Vespidae), which we collect in quantity throughout New England. Both hornets and ants (the latter are the dominant prey of S. purpurea) are hymenoptera, and have nearly identical C:N ratios (hornets: 3.97; common bog-dwelling ants [Tapinoma sessile and Myrmica lobifrons]: 3.37), but on average hornets have greater than 100 times the dry mass of these ants, and are easier to collect and process as a standardized food source. Additions of prey, either as large "pulses" or chronic "presses" are analogous to the enrichment and eutrophication that occur in aquatic "green" food webs in which phytoplankton abundance is boosted through addition of limiting nutrients. In "brown" food webs such as the Sarracenia microecosystem, detritus - not primary production - is at the base of the web, and our treatments boost this material as would happen through increases in arthropod prey capture78 or through nitrogen-enriched precipitation.

Proteomic analysis. Proteomic profiles of microbial communities are determined after separating the microbial fraction from the pitcher fluid, prey, and other detritus. The microbial "pellet" is subjected to SDS-PAGE (sodium dodecyl sulfate polyacrylamide gel) electrophoresis; bands are cut out and digested in-gel with trypsin. Tryptic peptides are subjected to LC-MS/MS (liquid chromatography tandem mass spectrometry) for peptide and protein identification. Absolute abundance of peptides and proteins are quantified using AQUA (Absolute QUAntification) analysis109.
      
## Specific Experiments

Experiment #1. Effects of nutrient enrichment on state changes and [O2] profiles. This experiment alters nutrient enrichment rates to characterize the [O2] profile and the transition to the anaerobic state. The experimental design is a one-way layout with 5 treatment groups: one control (no enrichment) and 4 enrichment levels (0.125, 0.25, 0.5, 1.0 mg prey added ml-1 d-1). One plant is assigned to each treatment group, and the entire set is replicated 6 times over successive weeks. [O2] is monitored continuously for 4 days to characterize state changes and tipping points under different enrichment rates. This experiment tracks a continuous [O2] profile but does not include proteomic analysis. The purpose of Experiment #1 is to identify an enrichment rate E that generates a long pre-tipping period before transition time T to the anaerobic state. This enrichment rate will be used in Experiments #2 - #4.

Experiment #2. Identification of early intervention time and characterization of aerobic and anaerobic proteomes. This experiment will use the single enrichment rate E determined from Experiment #1 and impose different intervention times I at which nutrient enrichment will be terminated. Thus, this experiment will identify the latest time I* at which it is possible to intervene and stop the transition to the anaerobic state by halting enrichment. The [O2] profile will again be monitored continuously over 10 days to measure the state of the system. From Experiment #1, the transition time T to the anaerobic state with no intervention will be known. We will use one control group (no prey addition) and ten levels of intervention time (all with the same enrichment rate E) as a proportion of T (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0). Six plants will be assigned randomly to each of the 11 treatments in a randomized one-way layout and [O2] profiles will be monitored continuously. In addition to the [O2] profiles, we will also characterize the protein profiles of aerobic and anaerobic pitchers in all 11 treatment groups at the end of the experiment.

After the plants are harvested, we will create proteomic profiles of the predominantly bacterial portion (centrifuged pellet) of the pitcher fluid from each plant, as described in General Protocols. Thus, 66 separate pellet-fraction samples will be analyzed by SDS-PAGE. After examining the SDS-PAGE profiles, approximately ten proteins that show dynamic patterns consistent with state change and five that do not change will be cut from the gel, subjected to in-gel tryptic digestion and a portion of the tryptic peptides will be analyzed by LC-MS/MS. Using these data, we will choose three identified peptides from each protein for peptide synthesis such that each synthesized peptide contains stable isotope labels (AQUA peptides) to distinguish them by mass from the native peptides. We will then quantify all 45 of the native peptides from the original samples using a known amount of each AQUA peptide spiked into the tryptic digest. The AQUA analysis of proteins that do not show changes will be used for normalization between samples. These data will be used to independently identify the current state of the system and forecast the time-to-transition.

We will use Sequest searches for initial identification of peptides; relevant scores including Xcorr and Delta Cn values will be given for each peptide. Other peptides will be identified by de novo sequencing using PepNovo; all PepNovo scores will likewise be given including any N- or C-terminal gaps. Mass error in ppm will be reported for each precursor ion. We will use standard multivariate analysis to search for distinctive proteomic profiles114 that characterize aerobic and anaerobic ecosystems, as well as ecosystems that developed with and without inputs of photosynthetic O2 and plant metabolites.
        
Experiment #3. Identification of diagnostic proteins. Using Experiments #1 and 2, we will have identified an enrichment rate E with a long pre-tipping period and an intervention time I* before which mitigation and termination of enrichment will prevent eutrophication. Experiment #3 will characterize the mean and variance of the protein profile before and after I*. We are especially interested in identifying proteins that increase rapidly in abundance (or variance) well before the onset of flickering in [O2] and before the transition time T from the aerobic to the anaerobic state.

A cohort of 100 plants all will be fed at rate E (determined from Experiment #1), with intervention time I* determined from Experiment #2, although no intervention will be used in this "press" experiment so that we can contrast proteins before and after the state change. At seven times before I* and three times after I*, we will harvest 10 randomly chosen plants. At each prescribed harvest time, we will measure [O2] and collect samples from each plant for proteomic screening using both SDS-PAGE and AQUA analysis. This experiment will identify proteins that rise quickly in abundance during the pre- I* period and can be used as early indicators of a future tipping point. Because different plants will be harvested at each time period, this is a one-way ANOVA design, with pre-and post- I* a priori contrasts. A randomization test will be used to determine whether variances in protein expression differ through time. During these analyses we will use the data from the AQUA peptides and from known amounts of protein standards, such as bovine serum albumin, to approximate the amount of protein in a given coomassie-stained SDS-PAGE gel band. The reason for doing this is to provide a fast "real-time" assay based just on expression in the SDS-PAGE. This rapid assay will be used in Experiment #4.

Experiment # 4. Proof-of-application. This experiment will provide a benchmark test of our methods and their ability to correctly identify tipping points. A cohort of 100 plants will each be fitted with [O2] probes and started on the enrichment regime. Two times per day, we will collect 3 plants each, pool their contents, and conduct a rapid screen in the lab with SDS-PAGE for the diagnostic proteins that were identified in Experiment #3. We will use the protein expression in the gel to delineate an "early" and a "late" mitigation strategy. As soon as diagnostic proteins measured in the SDS-gels are at abundances that signal we are at 0.5×I* - approximately one-half of the way to the latest intervention time - we will randomly select one third of the remaining plants for mitigation and termination of enrichment (the "early" mitigation strategy). We will continue to harvest plants from the remainder of the cohort and monitor proteins. As soon as diagnostic proteins signal we are at 0.75 times I*, we will randomly select one half of the remaining plants for mitigation and termination of enrichment (the "late" mitigation strategy). The remaining plants (approximately one sixth to one third of the original cohort) will continue to be enriched. We will monitor [O2] in all 3 groups (no-mitigation control, early mitigation, late mitigation) until all plants reach a new [O2] equilibrium. If the protein markers are successful, the proportion of food webs that remain aerobic will be significantly higher in the two mitigation treatments than in the no-mitigation control.
---
output: github_document
---

# EML <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Build Status](https://travis-ci.org/ropensci/EML.svg?branch=master)](https://travis-ci.org/ropensci/EML)
[![Windows build status](https://ci.appveyor.com/api/projects/status/u2gw24yfkvxgny96?svg=true)](https://ci.appveyor.com/project/cboettig/eml)
[![codecov](https://codecov.io/gh/ropensci/EML/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/EML)
[![CRAN status](https://www.r-pkg.org/badges/version/EML)](https://cran.r-project.org/package=EML)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/RNeXML)
[![DOI](https://zenodo.org/badge/10894022.svg)](https://zenodo.org/badge/latestdoi/10894022)

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```


```{r include=FALSE}
has_pandoc <-  rmarkdown::pandoc_available()
```




EML is a widely used metadata standard in the ecological and environmental sciences. We strongly recommend that interested users visit the [EML Homepage](https://eml.ecoinformatics.org/) for an introduction and thorough documentation of the standard. Additionally, the scientific article *[The New Bioinformatics: Integrating Ecological Data from the Gene to the Biosphere (Jones et al 2006)](https://doi.org/10.1146/annurev.ecolsys.37.091305.110031)* provides an excellent introduction into the role EML plays in building metadata-driven data repositories to address the needs of highly heterogeneous data that cannot be easily reduced to a traditional vertically integrated database. At this time, the `EML` R package provides support for the serializing and parsing of all low-level EML concepts, but still assumes some familiarity with the EML standard, particularly for users seeking to create their own EML files. We hope to add more higher-level functions which will make such familiarity less essential in future development.


## Notes on the EML v2.0 Release

`EML` v2.0 is a complete re-write which aims to provide both a drop-in replacement for the higher-level functions of the existing EML package while also providing additional functionality.  This new `EML` version uses only simple and familiar list structures (S3 classes) instead of the more cumbersome use of S4 found in the original `EML`.  While the higher-level functions are identical, this makes it easier to for most users and developers to work with `eml` objects and also to write their own functions for creating and manipulating EML objects.  Under the hood, `EML` relies on the [emld](https://github.com/ropensci/emld/) package, which uses a Linked Data representation for EML.  It is this approach which lets us combine the simplicity of lists with the specificity required by the XML schema.    

This revision also supports the **[recently released EML 2.2.0 specification](https://eml.ecoinformatics.org/whats-new-in-eml-2-2-0.html)**.  

# Creating EML


```{r message=FALSE, warning=FALSE}
library(EML)
```



## A minimal valid EML document:

```{r}
me <- list(individualName = list(givenName = "Carl", surName = "Boettiger"))
my_eml <- list(dataset = list(
              title = "A Minimal Valid EML Dataset",
              creator = me,
              contact = me)
            )


write_eml(my_eml, "ex.xml")
eml_validate("ex.xml")
```



## A Richer Example

Here we show the creation of a relatively complete EML document using `EML`.  This closely parallels the function calls shown in the original EML [R-package vignette](https://docs.ropensci.org/EML/articles/creating-EML.html).  


## `set_*` methods

The original EML R package defines a set of higher-level `set_*` methods to facilitate the creation of complex metadata structures.  `EML` provides these same methods, taking the same arguments for `set_coverage`, `set_attributes`, `set_physical`, `set_methods` and `set_textType`, as illustrated here:

### Coverage metadata

```{r}
geographicDescription <- "Harvard Forest Greenhouse, Tom Swamp Tract (Harvard Forest)"
coverage <- 
  set_coverage(begin = '2012-06-01', end = '2013-12-31',
               sci_names = "Sarracenia purpurea",
               geographicDescription = geographicDescription,
               west = -122.44, east = -117.15, 
               north = 37.38, south = 30.00,
               altitudeMin = 160, altitudeMaximum = 330,
               altitudeUnits = "meter")
```


### Reading in text from Word and Markdown

We read in detailed methods written in a Word doc.  This uses EML's docbook-style markup to preserve formatting of paragraphs, lists, titles, and so forth. (This is a drop-in replacement for EML `set_method()`)

```{r eval=has_pandoc}
methods_file <- system.file("examples/hf205-methods.docx", package = "EML")
methods <- set_methods(methods_file)
```

We can also read in text that uses Markdown for markup elements:

```{r eval=has_pandoc}
abstract_file <-  system.file("examples/hf205-abstract.md", package = "EML")
abstract <- set_TextType(abstract_file)
```



### Attribute Metadata from Tables

Attribute metadata can be verbose, and is often defined in separate tables (e.g. separate Excel sheets or `.csv` files).
Here we use attribute metadata and factor definitions as given from `.csv` files.

```{r}
attributes <- read.table(system.file("extdata/hf205_attributes.csv", package = "EML"))
factors <- read.table(system.file("extdata/hf205_factors.csv", package = "EML"))
attributeList <- 
  set_attributes(attributes, 
                 factors, 
                 col_classes = c("character", 
                                 "Date",
                                 "Date",
                                 "Date",
                                 "factor",
                                 "factor",
                                 "factor",
                                 "numeric"))
```


### Data file format 

Though the `physical` metadata specifying the file format is extremely flexible, the `set_physical` function provides defaults appropriate for `.csv` files. 
DEVELOPER NOTE: ideally the `set_physical` method should guess the appropriate metadata structure based on the file extension.  

```{r}
physical <- set_physical("hf205-01-TPexp1.csv")
```

## Generic construction

In the `EML` R package, objects for which there is no `set_` method are constructed using the `new()` S4 constructor.  This provided an easy way to see the list of available slots.  In `eml2`, all objects are just lists, and so there is no need for special methods.  We can create any object directly by nesting lists with names corresponding to the EML elements.  Here we create a `keywordSet` from scratch:

```{r}
keywordSet <- list(
    list(
        keywordThesaurus = "LTER controlled vocabulary",
        keyword = list("bacteria",
                    "carnivorous plants",
                    "genetics",
                    "thresholds")
        ),
    list(
        keywordThesaurus = "LTER core area",
        keyword =  list("populations", "inorganic nutrients", "disturbance")
        ),
    list(
        keywordThesaurus = "HFR default",
        keyword = list("Harvard Forest", "HFR", "LTER", "USA")
        ))
```

Of course, this assumes that we have some knowledge of what the possible terms permitted in an EML keywordSet are!   Not so useful for novices.  We can get a preview of the elements that any object can take using the `emld::template()` option, but this involves a two-part workflow.  Instead, `eml2` provides generic `construct` methods for all objects. 

## Constructor methods


For instance, the function `eml$creator()` has function arguments corresponding to each possible slot for a creator.  This means we can rely on *tab completion* (and/or autocomplete previews in RStudio) to see what the possible options are.  `eml$` functions exist for all complex types.  If `eml$` does not exist for an argument (e.g. there is no `eml$givenName`), then the field takes a simple string argument. 

### Creating parties (creator, contact, publisher)

```{r}
aaron <- eml$creator(
  individualName = eml$individualName(
    givenName = "Aaron", 
    surName = "Ellison"),
  electronicMailAddress = "fakeaddress@email.com")
```


```{r}
HF_address <- eml$address(
                  deliveryPoint = "324 North Main Street",
                  city = "Petersham",
                  administrativeArea = "MA",
                  postalCode = "01366",
                  country = "USA")
```

```{r}
publisher <- eml$publisher(
                 organizationName = "Harvard Forest",
                 address = HF_address)
```


```{r}              
contact <- 
  list(
    individualName = aaron$individualName,
    electronicMailAddress = aaron$electronicMailAddress,
    address = HF_address,
    organizationName = "Harvard Forest",
    phone = "000-000-0000")

```

### Putting it all together


```{r}
my_eml <- eml$eml(
           packageId = uuid::UUIDgenerate(),  
           system = "uuid",
           dataset = eml$dataset(
               title = "Thresholds and Tipping Points in a Sarracenia",
               creator = aaron,
               pubDate = "2012",
               intellectualRights = "http://www.lternet.edu/data/netpolicy.html.",
               abstract = abstract,
               keywordSet = keywordSet,
               coverage = coverage,
               contact = contact,
               methods = methods,
               dataTable = eml$dataTable(
                 entityName = "hf205-01-TPexp1.csv",
                 entityDescription = "tipping point experiment 1",
                 physical = physical,
                 attributeList = attributeList)
               ))

```

## Serialize and validate

We can also validate first and then serialize:

```{r}
eml_validate(my_eml)
write_eml(my_eml, "eml.xml")
```



## Setting the version

EML will use the latest EML specification by default.  To switch to a different version, use `emld::eml_version()` 

```{r}
emld::eml_version("eml-2.1.1")
```

Switch back to the 2.2.0 release: 

```{r}
emld::eml_version("eml-2.2.0")
```


```{r include = FALSE}
unlink("eml.xml")
unlink("ex.xml")
codemetar::write_codemeta()
```

---
title: Working with Ecological Metadata Language Files from R
author:
  - name: Carl Boettiger
    email: cboettig@berkeley.edu
    affiliation: ucb
    footnote: Corresponding Author
address:
  - code: ucb
    address: Department of Environmental Science, Policy, and Management, University of California, Berkeley, Berkeley, CA
abstract: |
  This is the abstract.

  It consists of two paragraphs.

journal: "Methods in Ecology and Evolution"
date: "`r Sys.Date()`"
bibliography: references.bib
output: rticles::elsevier_article
---

---
title: "Project Description"
subtitle: "Realizing the Promise of Machine-Readable Metadata: Ecoinformatics for the Rest of Us"
# author: "Carl Boettiger"
bibliography: "references.bib"
output: pdf_document
header-includes:
   - \usepackage{multicol}
   - \usepackage{bold-extra}

---



Despite a sea change in the availability of public data collected in ecological research, broad synthesis of ecological data remains a rare and largely manual endeavor. Realizing the potential of these unprecedented new data sources requires thorough and accurate machine-readable data and metadata files, and the tools to make use of them.  Here, I describe user friendly R tools that will enable more ecologists to both take advantage of and contribute to growth of well-annotated public ecological data. 

Many factors drive the rapid expansion of data availability in ecological research: from drones to satellite images, new technology for gene sequencing, large-scale simulations of climate and other processes [@Overpeck2011], to cultural changes such as the adoption of public data archiving [@Fairbairn2010; @Whitlock2010] and the significant investment in a national network of ecological [NEON, @Keller2008] and ocean [OOI, @Witze2013] observatories.  These factors also reflect the range and heterogeneity in ecological data, from the microscopic to the planetary.  How can researchers combine laboratory measurements of species environmental tolerances, aerial estimates of species occurrence, and simulated projections of changing climate?  How do we even find if the relevant data exist in these rapidly expanding archives? 

Rich, machine-readable metadata greatly facilitates discovery, integration and synthesis of the kind of heterogeneous data inherent to ecological research. Metadata provides all of the context required to understand raw data files -- without it data may be little more than meaningless columns of numbers.  Metadata tells us what column labels mean, what units were used. Metadata tells us which species, geographic areas, and time periods are covered by the study, who performed the research, and how.  Data search and discovery requires good metadata: we are looking for measurements of body mass on species found in this region.  Data integration and synthesis require good metadata: we want to combine nine different data sets on different species found in the same region into a single data table, so we must know the species involved, the units of measurement, and so forth.  

This metadata is most effective in a format that can be parsed and interpreted by a computer program.  While conventions have primarily relied on scientific papers or the occasional README file to provide this information, machine-readable metadata formats can help ensure these records are complete and accurate and can automate the often tedious process of recording metadata. A machine-readable format can be indexed by specialized search engines, allowing a user to identify all data sets that include measurements of a specific variable on a specific species or from a certain geographic region.  A machine-readable format can even help automate the process of data integration and synthesis, such as standardizing units when combining data columns from different studies.  Such automation is helpful even in small studies, it is crucial if ecological research is to keep up with our rapidly expanding data.

### What is EML?

The Ecological Metadata Language [EML, @Jones2006], originally designed around previous work from the Ecological Society of America [@Michener1997], seeks to tackle this challenge.  EML is an is an XML-based schema which provides a rich, extensible, machine-readable format to describe a wide range of data types, including tabular data, spatial raster and vector data, methods, protocols, software or citations.  EML files can be *validated* against the schema, ensuring this information conforms to the predictable structure and can thus be read by any computer software implementing the EML specification.  Recognizing the importance of thorough, machine-readable metadata, key providers of ecological data such as the NSF-supported Long Term Ecological Research (LTER) sites [@LTER_EML], and the NSF's recently launched National Ecological Observatory Network [NEON, @Keller2008; @NEON] now provide rich EML metadata accompanying their data products.   

More than XX individual datasets have been descibed in EML over the past two decades.  EML metadata are required for all data produced by LTER projects (representing about about $12 million NSF support annually), by all NEON sites (about $469 million NSF support) and all projects supported by the NSF Arctic programs (about $40 million annually) (https://www.nsf.gov/pubs/2016/nsf16055/nsf16055.jsp).  

Unfortunately, this XML-based format remains ignored or inaccessible to many ecologists, who lack tools and training to capitalize on these advances in computation and data management.  Few are familiar with the skills required to create EML descriptions of their own data, or to fully benefit from the rich information provided in EML [@Hernandez2012].  NEON relies on dedicated professionals to generate EML metadata, while ecologists working with LTER sites benefit from the support of informatics staff in creating EML files describing the data they contribute to LTER data repositories.  Consequently, LTER sites contribute the majority of EML files found in the DataONE Network (which indexes most major repositories of ecological and environmental data), Fig 1.  This highlights the need for tools for ecologists to create EML without assistance from experts or prior expertise in working with XML documents and schema. 


![Number of EML metadata files by member nodes for the DataONE Network of the major ecological and environmental data repositories.  EML is used across a wide range of data repositories shown here, all though repositories in which ecoinformatics staff assist in preparation of metadata, such as LTER, provide significantly more metadata than those which primarly host direct community contributions, such as KNB.](img/repo_counts.pdf)


### Relation to existing software tools 

<!--
3. How the proposed software compares to alternative or existing elements (including other commercial and research solutions) and what are the limitations of these existing elements.
-->

No extensible tool or utility for reading and writing EML files exists that integrates into the workflow of most ecologists.  Researchers can use the Java application `Morpho`, which provides a graphical user interface (GUI) for creating EML files describing their own data -- this provides an excellent introduction to creating EML, but can be tedious and the interface prohibits automation and templating possible with a programmatic tool. Moreover, `Morpho` can only write new EML files, but not help users parse and analyze existing ones.  The open-source software platform `MetaCat` is available for data repositories which wish to serve and index data files described in EML [@metacat].  This allows data hosts to present users with a nice web-based search interface for discovering data (`MetaCat-UI`) which is already in use by repositories such as KNB and DataONE (see Fig 2).  However, this web-based platform is ultimately software for data repositories to use and deploy, not software for end users.  The environment cannot be easily deployed or extended by research ecologists looking to to work with a generic collection of EML files, cannot perform arbitrary queries or assist in data integration and synthesis tasks.  @Weitz2016 identify these factors (a focus on servers, limited extensibility, and a software stack that is foreign to most researchers in the field of application) as the leading reasons why cyberinfrastructure ends up under-utilized by communities they seek to serve.  Most importantly, `Morpho` and `MetaCat` do little to empower the researcher to take advantage of machine-readable metadata in their own work: instead, these tools work in concert to allow other researchers to take advantage of the user's data -- a noble objective that is often seen as self-sacrificing and a disincentive to data sharing [@VanNoorden2013].  Researchers need an extensible tool that can easily be incorporated into their existing workflow and facilitate their own research.





![A powerful EML metadata-driven search interface from DataONE, based on `MetaCat` platform.  The web interface identifies data files containing measurements of mass of *Cyprinidae* fishes in the geographical regions indicated by the grid. Such a search requires metadata that is thorough, granular, and machine readable. An R package for EML would both make it easier for more researchers to create this EML metadata, perform search queries like that shown here and extend both search and integration tasks through a programmable interface](img/dataone_search.png)




## An R package for EML

This proposal seeks to bridge this gap between the skills and workflows of most research ecologists and the powerful expressiveness of EML-annotated data being generated by NSF's major investments in informatics-supported data generation such as NEON and LTER through the development of an R package for reading, writing and manipulating EML files.  The R language is widely used throughout the ecological research community and has proven an effective platform to engage both users and developers [@Mair2015; @ropensci].  This proposal would support the development and maintenance of the EML R package.  A note on semantics: following the convention for naming R packages dedicated to parsing and serializing a given format, the proposed package name and the data format are both referred to by the acronym EML.  Throughout, I will refer to the data format as simply "EML", and the software package specifically as "the EML R package." 


EML parsing in R would enable users to analyze metadata in the same programmatic environment many already use for data analysis.  This kind of programmable environment is essential for users to truly leverage the power of machine-readable metadata in large scale analyses. For instance, this would make it possible for a user to create scripts which can combine data frames coming from different studies which may use different column names or different units for the same variable.  Currently such data integration tasks require manual inspection of the metadata for each data file to verify things like column names and units, a process which does not easily scale to arbitrarily large numbers of separate data files. Or a user might search for a particular phrase being used in the protocol section of a corpus of EML files, or automatically extract citations from metadata files into format recognized by a reference manager.  Once metadata is available in a predictably structured R objects, users can take advantage of the familiar, programmable environment provided by R to automate these tasks.  

The ability to create EML from within the R environment is equally powerful.  Again, this is a natural setting to generate metadata since many users are already working with their data in R.  Metadata creation and curation is infamously tedious, making a scripted approach that can automate repetitive tasks a natural choice. Reuse of common metadata elements can avoid tedious data entry (including individual and institutional information, citations, and common data formats), automatically detect metadata from existing data files (such as species coverage or geographic range).  Helper routines within the package can also improve the quality of metadata created. For instance, EML includes a description of the taxonomic coverage of the data being described.  This allows search engines such as MetaCat or the EML R package search functions to identify all data on a given species.  This metadata is most useful when it includes a complete taxonomic description, e.g. taxonomic ranks involved, such that a search for *Aves* returns all bird species.  The EML R package can leverage existing functions in the R package `taxize` to automatically generate the rank classification metadata [@Chamberlain2013].   Package users have the greatest incentive to generate careful and thorough metadata when they can benefit directly from these annotations.  By generating EML files for their own data, users can later leverage the EML parsing functions to help search and synthesize old data sets.  The EML package also integrates with `dataone` R package to automate uploading data and metadata files to a public or private repository archive in the DataONE network, such as the KNB.  This gives users a convenient way to back-up, archive, and share their data files using configurable permissions options and authentication provided by DataONE, and to release data publicly under a permanent DOI, facilitating discovery and citation by other researchers.  Once users can realize these benefits for themselves, opting in to public data sharing can be only a click away, with data already annotated. In this approach, the public benefits of data annotation and data sharing come as a by-product of a tool that puts the user's own needs first. 
These use cases are summarized in Box 1.  This proposal outlines the project plan to implement this functionality into a robust, user-friendly, sustainable software tool.  


# Software Architecture and Engineering Process #

The software design follows the format of an R package, the definitive community standard for packaging software to R users. This provides a widely recognized layout for the project, cross platform support, and well defined norms for documentation, testing, distribution, release and evaluation described here.  

**Documentation**: The project uses the widely adopted `roxygen2` framework for composing and maintaining all documentation in the software files, consistent with literate programming practice. Keeping code and documentation together helps ensure documentation is up-to-date and consistent with the software.  The automated testing suite also checks the documentation for consistency with the code.  R generates both HTML and PDF forms of documentation from these files.  A package README, Code of Conduct, and Contributing Guide provide further documentation relevant to other developers.  The project also uses the `pkgdown` package to generate a GitHub-hosted website from the package documentation.  

**Testing and validation**: The package uses `testthat` [@testthat], a widely used unit-testing framework for R packages, and `covr` [@covr] to measure test coverage.  The test suite is fully automated both locally and through the Travis Continuous Integration (CI) service.


**Security and Trustworthiness**: Open source development practices ensure users and developers can inspect all aspects of the code and tests.  EML files conform to open standards that permit validation and encryption of data. 

**Provenance**: Semantic versioning, a NEWS log summarizing changes, and the automatic generation of archival DOIs for each version help researchers both access and cite specific version of the EML package that they use.  Additionally, the EML software itself helps researchers track and document provenance of their data and methods in the machine-readable EML format.  

**Reproducibilty**: Our development process is committed to making sure analyses using the EML package remain as reproducible as possible in the future. All effort is made to maintain a backwards-compatible function interface by avoiding changes to the function signature that could break existing code whenever possible.  An increase of the major version number indicates when a new version may break backwards compatibility, and older versions of the software are archived by CRAN and in the Zenodo data repository at the DOI indicated for the version. 



\pagebreak

# References
---
title: "Reading and parsing EML"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Parsing EML}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


# listviewer

If the `listviewer` package is installed it can provide a convenient way to view and edit EML:

```{r, eval = FALSE}
f <- system.file("tests/eml.xml", package = "emld")
eml <- read_eml(f)
listviewer::jsonedit(eml)
```

# Parsing an EML file

```{r}
f <- system.file("xsd/test/eml-i18n.xml", package = "EML")
eml <- read_eml(f)
```

Here we request all `temporalCoverage` elements occurring in the anywhere in the `eml` document:

```{r}
temporalCoverage <- eml_get(eml, "temporalCoverage")
temporalCoverage
```


Any EML element can be extracted in this way. Let's try an example metadata file for a `dataset` that documents 11 seperate `dataTables`: 

```{r}
hf001 <- system.file("examples/hf001.xml", package="EML") 
eml_HARV <- read_eml(hf001)
```

How many `dataTable` entities are there in this dataset?

```{r}
dt <- eml_get(eml_HARV, "dataTable")
length(dt)
```

We can iterate over our list of `dataTable` elements to extract relevant metadata, such as the `entityName` or the download `url`:

```{r}
entities <- eml_get(eml_HARV, "dataTable.entityName")
urls <- sapply(dt, eml_get, "url")
```


Note that the latter example is the same as providing the more verbose arbument that specificies exactly where the `url` of interest is located:

```{r}
urls <- sapply(dt, function(x) x@physical[[1]]@distribution[[1]]@online@url)
```

this verbose syntax can be useful if there are multiple `url` elements in each `dataTable` metadata, and we are trying to get only certain ones and not others. Specifying the exact path in this way can also improve the speed of the command.  For these reasons, programmatic use should consider this format, while the much simpler `eml_get` example shown above is practical for most interactive applications.

Although the default return type for `eml_get` is just the S4 object (whose `print` method displays the corresponding XML structure used to represent that metadata), for a few commonly accessed complex elements, `eml_get` returns a more convenient `data.frame`.  For instance, the `attributeList` describing the metadata for every column in an EML document is returned as a pair of `data.frame`s, one for all the attributes, and an second optional `data.frame` defnining the levels for the factors, if any are used.  Let's take a look: 

Here we get the `attributeList` for each `dataTable` in the dataset. We check the length to confirm we get one `attributeList` for each `dataTable`

```{r}
attrs <- eml_get(dt, "attributeList") 
length(attrs)
attrs[[1]]
```

(Note, we could have passed this argument the original `eml_HARV` instead of `dt` here, since we know all `attributeList` elements are inside `dataTable` elements, but this is a bit more explicit and a bit faster.)

This returned `data.frame` object containing the attribute metadata for the first table (hence the `[[1]]`, though `attrs` contains this metadata for all 11 tables now.)  This is the same result we would have gotten using the more explicit call to the helper function `get_attributes()`:

```{r}
get_attributes(eml_HARV@dataset@dataTable[[1]]@attributeList)
```

---
title: "Creating EML"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Creating EML}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r include=FALSE}
has_pandoc <-  rmarkdown::pandoc_available()
```


```{r}
library(EML)
library(emld)
```

Here we construct a common EML file, including: 

- Constructing more complete lists of authors, publishers and contact.
- Summarizing the geographic, temporal, and taxonomic coverage of the dataset
- Reading in pages of methods descriptions from a Word document
- Adding arbitrary additional metadata
- Indicating the canonical citation to the paper that should be acknowledged when the data is re-used.  
- Conversion between EML and other metadata formats, such as NCBII and ISO standards.

In so doing we will take a modular approach that will allow us to build up our metadata from reusable components, while also providing a more fine-grained control over the resulting output fields and files.  

### Overview of the EML hierarchy

A basic knowledge of the components of an EML metadata file is essential to being able to take full advantage of the language. While more complete information can be found in the [official schema documentation](https://knb.ecoinformatics.org/#external//emlparser/docs/eml-2.1.1/index.html), here we provide a general overview of commonly used metadata elements most relevant to describing data tables.  

This schematic shows each of the metadata elements we will generate.  Most these elements have sub-components (e.g. a 'publisher' may have a name, address, and so forth) which are not shown for simplicity.  Other optional fields we will not be generating in this example are also not shown.   

```yaml
- eml
  - dataset
    - creator
    - title
    - publisher
    - pubDate
    - keywords
    - abstract 
    - intellectualRights
    - contact
    - methods
    - coverage
      - geographicCoverage
      - temporalCoverage
      - taxonomicCoverage
    - dataTable
      - entityName
      - entityDescription
      - physical
      - attributeList

```



### Our example: (re)-creating hf205.xml

In this example, we will use R to re-generate
the EML metadata originally published by [Ellison _et al_
(2012), HF205](https://harvardforest.fas.harvard.edu/harvard-forest-data-archive)
through the Harvard Forest Long Term Ecological Research
Center, accompanying the PNAS paper [Sirota _et al_
(2013)](https://doi.org/10.1073/pnas.1221037110). We have made only
a few modifications to simplify the presentation of this tutorial, 
so our resulting EML will not be perfectly identical to the original.


### Our strategy

We will build this EML file from the bottom up, starting with the two main components
of a `dataTable` indicated above: the `attributeList` and the `physical` file type.
We will then slip these two pieces into place inside a `dataTable` element, and slip
that into our `eml` element along with the rest of the generic metadata, much like
building a puzzle or nesting a set of Russian dolls.  

The original metadata file was created in association with the publication in PNAS
based on a Microsoft Word document template that Harvard Forest provides
to the academic researchers.  Metadata from this template is then read off
by hand and an EML file is generated using a combination of a commercial
XML editing platform (Oxygen) for commonly used higher-level elements,
and the Java platform `Morpho` provided by the EML development team for
lower level attribute metadata.


## Attribute Metadata

A fundamental part of EML metadata is a description of the attributes
(usually columns) of a text file (usually a csv file) containing the
data being described.  This is the heart of many EML files.

```{r}
attributes <-
tibble::tribble(
~attributeName, ~attributeDefinition,                                                 ~formatString, ~definition,        ~unit,   ~numberType,
  "run.num",    "which run number (=block). Range: 1 - 6. (integer)",                 NA,            "which run number", NA,       NA,
  "year",       "year, 2012",                                                         "YYYY",        NA,                 NA,       NA,
  "day",        "Julian day. Range: 170 - 209.",                                      "DDD",         NA,                 NA,       NA,
  "hour.min",   "hour and minute of observation. Range 1 - 2400 (integer)",           "hhmm",        NA,                 NA,       NA,
  "i.flag",     "is variable Real, Interpolated or Bad (character/factor)",           NA,            NA,                 NA,       NA,
  "variable",   "what variable being measured in what treatment (character/factor).", NA,            NA,                 NA,       NA,
  "value.i",    "value of measured variable for run.num on year/day/hour.min.",       NA,            NA,                 NA,       NA,
  "length",    "length of the species in meters (dummy example of numeric data)",     NA,            NA,                 "meter",  "real")

```

Every column (attribute) in the dataset needs an `attributeName` (column name, as it appears in the CSV file) and `attributeDefinition`, a longer description of what the column contains. Additional information required depends on the data type:

**Strings** (character vectors) data just needs a "definition" value, often the same as the `attributeDefinition` in this case.

**Numeric** data needs a `numberType` (e.g. "real", "integer"), and a unit. 

**Dates** need a date format.

**Factors** (enumerated domains) need to specify definitions for each of the code terms appearing in the data columns.  This does not fit so nicely in the above table, where each attribute is a single row, so if data uses factors (instead of non-enumerated strings), these definitions must be provided in a separate table.  The format expected of this table has three columns: `attributeName` (as before), `code`, and `definition`.  Note that `attributeName` is simply repeated for all codes belonging to a common attribute.  

In this case we have three attributes that are factors,  To make the code below more readable (aligning code and definitions side by side), we define these first as named character vectors, and convert that to a `data.frame`. (The `dplyr::frame_data` function also permits this more readable way to define data.frames inline).

```{r}
i.flag <- c(R = "real",
            I = "interpolated",
            B = "bad")
variable <- c(
  control  = "no prey added",
  low      = "0.125 mg prey added ml-1 d-1",
  med.low  = "0,25 mg prey added ml-1 d-1",
  med.high = "0.5 mg prey added ml-1 d-1",
  high     = "1.0 mg prey added ml-1 d-1",
  air.temp = "air temperature measured just above all plants (1 thermocouple)",
  water.temp = "water temperature measured within each pitcher",
  par       = "photosynthetic active radiation (PAR) measured just above all plants (1 sensor)"
)

value.i <- c(
  control  = "% dissolved oxygen",
  low      = "% dissolved oxygen",
  med.low  = "% dissolved oxygen",
  med.high = "% dissolved oxygen",
  high     = "% dissolved oxygen",
  air.temp = "degrees C",
  water.temp = "degrees C",
  par      = "micromoles m-1 s-1"
)

## Write these into the data.frame format
factors <- rbind(
data.frame(
  attributeName = "i.flag",
  code = names(i.flag),
  definition = unname(i.flag)
),
data.frame(
  attributeName = "variable",
  code = names(variable),
  definition = unname(variable)
),
data.frame(
  attributeName = "value.i",
  code = names(value.i),
  definition = unname(value.i)
)
)
```

With these two data frames in place, we are ready to create our `attributeList` element:

```{r}
attributeList <- set_attributes(attributes, factors, col_classes = c("character", "Date", "Date", "Date", "factor", "factor", "factor", "numeric"))
```


## Data file format 

The documentation of a `dataTable` also requires a description of the file format itself.  From where can the data file be downloaded?  Is it in CSV format, or TSV (tab-separated), or some other format? Are there header lines that should be skipped? This information documents the physical file itself, and is provided using the `physical` child element to the `dataTable`.  To assist in documenting common file types such as CSV files, the `EML` R package provides the function `set_physical`, which takes as arguments many of these common options.  By default these options are already set to document a standard `csv` formatted object, so we do not need to specify delimiters and so forth if our file conforms to that.  We simply provide the name of the file, which is used as the `objectName`.  (See examples for `set_physical()` for reading other common variations, analogous to the options covered in R's `read.table()` function.)

```{r}
physical <- set_physical("hf205-01-TPexp1.csv")
```

## Assembling the `dataTable`

Once we have defined the `attributeList` and `physical` file, we can now assemble the `dataTable` element itself. Unlike the old `EML` R package, in `EML` version 2.0 there is no need to call `new()` to create elements.  Everything is just a list.  Template lists for a given class can be viewed with the `emld::template()` function. 

```{r}
dataTable <- list(
                 entityName = "hf205-01-TPexp1.csv",
                 entityDescription = "tipping point experiment 1",
                 physical = physical,
                 attributeList = attributeList)
```


## Coverage metadata

One of the most common and useful types of metadata is coverage
information, specifying the temporal, taxonomic, and geographic coverage
of the data.  This kind of metadata is frequently indexed by data
repositories, allowing users to search for all data about a specific
region, time, or species.  In EML, these descriptions can take many forms,
allowing for detailed descriptions as well as more general terms when such
precision is not possible (such as geological epoch instead of date range,
or higher taxonomic rank information in place of species definitions.)

Most common specifications can
be made using the more convenient but less flexible `set_coverage()`
function in EML.  This function takes a date range or list of specific
dates, a list of scientific names, a geographic description and bounding
boxes, as shown here:


```{r}
geographicDescription <- "Harvard Forest Greenhouse, Tom Swamp Tract (Harvard Forest)"


coverage <- 
  set_coverage(begin = '2012-06-01', end = '2013-12-31',
               sci_names = "Sarracenia purpurea",
               geographicDescription = geographicDescription,
               west = -122.44, east = -117.15, 
               north = 37.38, south = 30.00,
               altitudeMin = 160, altitudeMaximum = 330,
               altitudeUnits = "meter")
```



## Creating methods

Careful documentation of the methods involved in the experimental design, measurement and collection of data are a key part
of metadata.  Though frequently documented in scientific papers, such method sections may be too brief or incomplete, and may become more readily disconnected from the data file itself.  Such documentation is usually written using word-processing software such as MS Word, LaTeX or markdown.  Users with `pandoc` installed (which ships as part of RStudio) can install the `rmarkdown` package to take advantage of its automatic conversion into the DocBook XML format used by EML.  Here we open a MS Word file with the methods and read this into our methods element using the helper function `set_methods()`.  While not used in this example, note that the `set_methods()` function also includes many optional arguments for documenting additional information about sampling, or relevant citations.

```{r eval=has_pandoc}
methods_file <- system.file("examples/hf205-methods.docx", package = "EML")
methods <- set_methods(methods_file)
```

```{r include=FALSE, eval=!has_pandoc}
## placeholder if pandoc is not installed
methods <- NULL
```


## Creating parties

Individuals and organizations appear in many capacities in an EML document.  Meanwhile, R already has a native object class, `person` for describing individuals, which it uses in citations and package descriptions, among other things. We can use native R function `person()` to create an R `person` object.  Often it is more convenient to use R's coercion function, `as.person()`, to turn a string with standardized notation into a `person` class (Though this is not always reliable, for instance, in surnames containing whitespace).  However it is constructed, a `person` class can then be coerced into the appropriate EML object like so:

```{r}
R_person <- person("Aaron", "Ellison", ,"fakeaddress@email.com", "cre", 
                  c(ORCID = "0000-0003-4151-6081"))
aaron <- as_emld(R_person)
```

Likewise this method can be applied to a list of `person` objects:   

```{r}
others <- c(as.person("Benjamin Baiser"), as.person("Jennifer Sirota"))
associatedParty <- as_emld(others)
associatedParty[[1]]$role <- "Researcher"
associatedParty[[2]]$role <- "Researcher"
```

Note that R only permits certain codes such as `ctb` be be given in square brackets or as the `role` slot in a `person` object.  

We can instead always use the list approach to create any of these elements, instead of the 
shorthand coercion methods shown above.  This permits a bit more flexibility, particularly for constructing elements
where we want to include more metadata than R's `person` object knows about.  Here we define an `address` element
first, since we can then re-use that element in defining both the contact person and publisher of the dataset:

```{r}
HF_address <- list(
                  deliveryPoint = "324 North Main Street",
                  city = "Petersham",
                  administrativeArea = "MA",
                  postalCode = "01366",
                  country = "USA")
```


```{r}
publisher <- list(
                 organizationName = "Harvard Forest",
                 address = HF_address)
```


```{r}              
contact <- 
  list(
    individualName = aaron$individualName,
    electronicMailAddress = aaron$electronicMailAddress,
    address = HF_address,
    organizationName = "Harvard Forest",
    phone = "000-000-0000")

```


## Creating a `keywordSet`

Constructing the `keywordSet` is just a list of lists.  Note that everything is a list.

```{r}
keywordSet <- list(
    list(
        keywordThesaurus = "LTER controlled vocabulary",
        keyword = list("bacteria",
                    "carnivorous plants",
                    "genetics",
                    "thresholds")
        ),
    list(
        keywordThesaurus = "LTER core area",
        keyword =  list("populations", "inorganic nutrients", "disturbance")
        ),
    list(
        keywordThesaurus = "HFR default",
        keyword = list("Harvard Forest", "HFR", "LTER", "USA")
        ))
```


Lastly, some of the elements needed for `eml` object can simply be given as text strings.

```{r}
pubDate <- "2012" 

title <- "Thresholds and Tipping Points in a Sarracenia 
Microecosystem at Harvard Forest since 2012"

abstract <- "The primary goal of this project is to determine
  experimentally the amount of lead time required to prevent a state
change. To achieve this goal, we will (1) experimentally induce state
changes in a natural aquatic ecosystem - the Sarracenia microecosystem;
(2) use proteomic analysis to identify potential indicators of states
and state changes; and (3) test whether we can forestall state changes
by experimentally intervening in the system. This work uses state-of-the
art molecular tools to identify early warning indicators in the field
of aerobic to anaerobic state changes driven by nutrient enrichment
in an aquatic ecosystem. The study tests two general hypotheses: (1)
proteomic biomarkers can function as reliable indicators of impending
state changes and may give early warning before increasing variances
and statistical flickering of monitored variables; and (2) well-timed
intervention based on proteomic biomarkers can avert future state changes
in ecological systems."  

intellectualRights <- "This dataset is released to the public and may be freely
  downloaded. Please keep the designated Contact person informed of any
plans to use the dataset. Consultation or collaboration with the original
investigators is strongly encouraged. Publications and data products
that make use of the dataset must include proper acknowledgement. For
more information on LTER Network data access and use policies, please
see: http://www.lternet.edu/data/netpolicy.html."
```


Many of these text fields can instead be read in from an external file that has richer formatting, such as we did with
the `set_methods()` step.  Any text field containing a slot named `section` can import text data from a MS Word `.docx` file, markdown file, or other file format recognized by [Pandoc](https://pandoc.org) into that element.  For instance, here we import the same paragraph of text shown above for `abstract` from an external file (this time, a markdown-formatted file) instead:

```{r eval=has_pandoc}
abstract_file <-  system.file("examples/hf205-abstract.md", package = "EML")
abstract <- set_TextType(abstract_file)
```


We are now ready to add each of these elements we have created so far into our `dataset` element, like so:


```{r}
dataset <- list(
               title = title,
               creator = aaron,
               pubDate = pubDate,
               intellectualRights = intellectualRights,
               abstract = abstract,
               associatedParty = associatedParty,
               keywordSet = keywordSet,
               coverage = coverage,
               contact = contact,
               methods = methods,
               dataTable = dataTable)
```





With the `dataset` in place, we are ready to declare our root `eml` element.  In addition to our `dataset` element we have already built, all we need is a packageId code and the system on which it is based.  Here we have generated a unique id using the standard `uuid` algorithm, which is available in the R package `uuid`.  

```{r}
eml <- list(
           packageId = uuid::UUIDgenerate(),
           system = "uuid", # type of identifier
           dataset = dataset)

```



With our `eml` object fully constructed in R, we can now check that it is valid, conforming to all criteria set forth in the EML Schema.  This will ensure that other researchers and other software can readily parse and understand the contents of our metadata file:


```{r}
write_eml(eml, "eml.xml")
```


```{r}
eml_validate("eml.xml")
```


The validator returns a status `0` to indicate success.  Otherwise, the first error message encountered will be displayed. The most common reason for an error is probably the omission of a required metadata field. 


To take the greatest advantage of EML, we should consider depositing our file in a [Metacat](https://knb.ecoinformatics.org/knb/docs/intro.html)-enabled repository, which we discuss in the next vignette on using EML with data repositories.

```{r include=FALSE}
unlink("eml.xml")
```
---
title: "Working with Units"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Working with Units}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Overview


One essential role of EML metadata is in precisely defining the units in which data is measured. To make sure these units can be understood by (and thus potentially converted to other units by) a computer, it's necessary to be rather precise about our choice of units.  EML knows about a lot of commonly used units, referred to as "standardUnits," already.  If a unit is in EML's `standardUnit` dictionary, we can refer to it without further explanation as long as we're careful to use the precise `id` for that unit, as we will see below.

Sometimes data involves a unit that is not in `standardUnit` dictionary.  In this case, the metadata must provide additional information about the unit, including how to convert the unit into the SI system.  EML uses an existing standard, [stmml](http://www.ch.ic.ac.uk/rzepa/codata2/), to represent this information, which must be given in the `additionalMetadata` section of an EML file.  The `stmml` standard is also used to specify EML's own `standardUnit` definitions.


## Add a custom unit to EML


```{r}
library("EML")
```




```{r}
custom_units <- 
  data.frame(id = "speciesPerSquareMeter", 
             unitType = "arealDensity", 
             parentSI = "numberPerSquareMeter", 
             multiplierToSI = 1, 
             description = "number of species per square meter")


unitList <- set_unitList(custom_units)
```             


Start with a minimal EML document
```{r}
me <- list(individualName = list(givenName = "Carl", surName = "Boettiger"))
my_eml <- list(dataset = list(
              title = "A Minimal Valid EML Dataset",
              creator = me,
              contact = me),
              additionalMetadata = list(metadata = list(
                unitList = unitList
              ))
            )

```


```{r}
write_eml(my_eml, "eml-with-units.xml")
eml_validate("eml-with-units.xml")
```

Note: Custom units are widely misunderstood and misused in EML.  See examples from [custom-units](https://gist.github.com/amoeba/67a23818dfca49904c7a54b0632d76bc#file-all-arcticdata-customUnits)





## Reading EML: parsing unit information, including custom unit types


Let us start by examining the numeric attributes in an example EML file.  First we read in the file:


```{r}
f <- system.file("tests", emld::eml_version(), "eml-datasetWithUnits.xml", package = "emld")
eml <- read_eml(f)
```

We extract the `attributeList`, and examine the numeric attributes (e.g. those which have units): 



```{r include = FALSE}
# clean up
unlink("eml-with-units.xml")
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_eml.R
\name{read_eml}
\alias{read_eml}
\title{read_eml}
\usage{
read_eml(x, from = "xml")
}
\arguments{
\item{x}{path to an EML file}

\item{from}{explicit type for the input format. Possible values:
"xml", "json", "list", or "guess" with "xml" as the default.}
}
\value{
an emld object (list / S3 object)
}
\description{
Read an EML file into R as an emld object.
}
\examples{
f <- system.file("extdata", "example.xml", package = "emld")
eml <- read_eml(f)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_validate.R
\name{eml_validate}
\alias{eml_validate}
\title{eml_validate}
\usage{
eml_validate(eml, encoding = "UTF-8", schema = NULL)
}
\arguments{
\item{eml}{file path, xml_document,}

\item{encoding}{optional encoding for files, default UTF-8.}

\item{schema}{path to schema}
}
\value{
Whether the document is valid (logical)
}
\description{
eml_validate processes an EML document using the XSD schema for the
appropriate version of EML and determines if the document is schema-valid
as defined by the XSD specification
}
\note{
this function is simply an alias to `eml_validate` in `emld` package
}
\examples{
\donttest{

f <- system.file("extdata", "example.xml", package = "emld")

## validate file directly from disk:
eml_validate(f)

## validate an eml object:
eml <- read_eml(f)
eml_validate(eml)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_unitList.R
\name{set_unitList}
\alias{set_unitList}
\title{set_unitList}
\usage{
set_unitList(units, unitTypes = NULL, as_metadata = FALSE)
}
\arguments{
\item{units}{a data.frame describing the custom units, see details.}

\item{unitTypes}{optional, a data.frame defining any additional unitTypes not already defined}

\item{as_metadata}{logical, default FALSE. If true, returns an `additionalMetadata` element, see below.}
}
\value{
unitList list object
}
\description{
Define custom units, including new unitTypes.  Note that it is not necessary to define
most common units.
}
\details{
The units data.frame must have the following columns:
 - id: the referenced name of unit (singular). e.g. 'meter', 'second'
 - unitType: the base type of unit, e.g. 'length'.  If not from a standard type, a new unitType must be provided
 - multiplierToSI: the multiplicative constant to convert to the SI unit.
 - parentSI: the name of the parent SI unit, e.g. second.
 - description: a text string describing the unit of measure.
 The following columns are optional:
 - name: usually the same as the id of the unit, e.g. second
 - abbreviation: common abbreviation, e.g. s
 - constantToSI: an additive constant to convert to the equivalent SI unit. If not given, default is "0"

In practice, researchers may save these tables of custom units they frequently use in an external .csv
or other format and read them in to R for ready re-use.

The unitType table must have the following columns:
 - id: the name by which the unitType is referred to.
 - name: optional, default is same as the id
 - dimension: name of a base dimension of the unit
 - power: the power to which the dimension is raised (NA implies power of 1)
}
\examples{
## create the "unitType" table for custom unit
id <- c("speed", "speed", "acceleration", "acceleration", "frequency")
dimension <- c("length", "time", "length", "time", "time")
power <- c(NA, "-1", NA, "-2", "-1")
unitTypes <- data.frame(
  id = id, dimension = dimension,
  power = power, stringsAsFactors = FALSE
)

## Create the units table
id <- c("minute", "centimeter")
unitType <- c("time", "length")
parentSI <- c("second", "meter")
multiplierToSI <- c("0.0166", "1")
description <- c("one minute is 60 seconds", "centimeter is a 100th of a meter")
units <- data.frame(
  id = id, unitType = unitType, parentSI = parentSI,
  multiplierToSI = multiplierToSI, description = description,
  stringsAsFactors = FALSE
)

unitList <- set_unitList(units, unitTypes)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eml_get.R
\name{eml_get}
\alias{eml_get}
\title{eml_get}
\usage{
eml_get(x, element, from = "list", ...)
}
\arguments{
\item{x}{an EML object or child/descendant object}

\item{element}{name of the element to be extracted.
If multiple occurrences are found, will extract all}

\item{from}{explicit type for the input format. Possible values:
"xml", "json", "list", or "guess" with "list" as the default.}

\item{...}{additional arguments}
}
\description{
eml_get
}
\examples{
\donttest{
f <- system.file("tests", emld::eml_version(), "eml-datasetWithUnits.xml", package = "emld")
eml <- read_eml(f)
eml_get(eml, "physical")
eml_get(eml, "attributeList")

## The first argument need not be an "eml" class, it could be a child element; e.g.
eml_get(eml$dataset$dataTable, "physical")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_attributes.R
\name{table_to_r}
\alias{table_to_r}
\title{handsontable to r}
\usage{
table_to_r(table)
}
\arguments{
\item{table}{input table}
}
\description{
Takes a handsontable and converts to r data.frame for shiny app
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_physical.R
\name{detect_delim}
\alias{detect_delim}
\title{Automatically detect line delimiters in a text file}
\usage{
detect_delim(path, nchar = 1000)
}
\arguments{
\item{path}{(character) File to search for a delimiter}

\item{nchar}{(numeric) Maximum number of characters to read from disk when
searching}
}
\value{
(character) If found, the delimiter, it not, \\r\\n
}
\description{
This helper function was written expressly for \code{\link{set_physical}} to
be able to automate its \code{recordDelimiter} argument.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_TextType.R
\name{set_TextType}
\alias{set_TextType}
\title{set_TextType}
\usage{
set_TextType(file = NULL, text = NULL)
}
\arguments{
\item{file}{path to a file providing formatted input text, see details.}

\item{text}{a plain text character string which will be used directly as the content
of the node if no file is given}
}
\value{
a TextType object that can be coerced into any element inheriting from TextType, see examples
}
\description{
For any EML element of class TextType, this function can be used to generate
 the appropriate EML from a markdown-formatted file.
}
\details{
If the `rmarkdown` package is installed, then the input file can
be a Microsoft Word (.docx) file, a markdown file, or other file
recognized by Pandoc (see https://pandoc.org), which will automate the conversion
to a docbook. Otherwise, the input file should already be in docbook format (with
.xml or .dbk extension).  Note that pandoc comes pre-installed in RStudio and is
required for the rmarkdown package.
}
\examples{
\donttest{
## using a simple character string
a <- set_TextType(text = "This is the abstract")

## Using an external markdown file
f <- system.file("examples/hf205-abstract.md", package = "EML")
a <- set_TextType(f)

## Can also import from methods written in a .docx MS Word file.
f <- system.file("examples/hf205-abstract.docx", package = "EML")
a <- set_TextType(f)

## Documents with title headings use `section` instead of `para` notation
f <- system.file("examples/hf205-methods.docx", package = "EML")
d <- set_TextType(f)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_unitList.R
\name{get_unitList}
\alias{get_unitList}
\title{get_unitList}
\usage{
get_unitList(x = NULL)
}
\arguments{
\item{x}{an emld object}
}
\value{
a list with two data.frames: "units", a table defining unit names, types, and conversions to SI,
and "unitTypes", defining the type of unit. For instance, the unit table could define "Hertz" as a unit
of unitType frequency, and the unitType define frequency as a type whose dimension is 1/time.
}
\description{
get_unitList
}
\details{
If no unitList is provided, the function reads in the eml-unitDictionary defining all standard
units and unitTypes.  This provides a convenient way to look up standard units and their EML-recognized names
when defining metadata, e.g. in the table passed to `set_attributes()`.
}
\examples{

# Read in additional units defined in a EML file
\donttest{
f <- system.file("tests", emld::eml_version(),
  "eml-datasetWithUnits.xml",
  package = "emld"
)
eml <- read_eml(f)
unitList <- get_unitList(eml)

## Read in the definitions of standard units:
get_unitList()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_attributes.R
\name{set_attributes}
\alias{set_attributes}
\title{set_attributes}
\usage{
set_attributes(
  attributes,
  factors = NULL,
  col_classes = NULL,
  missingValues = NULL
)
}
\arguments{
\item{attributes}{a joined table of all attribute metadata}

\item{factors}{a table with factor code-definition pairs; see details}

\item{col_classes}{optional, list of R column classes ('ordered', 'numeric', 'factor', 'Date', or 'character', case sensitive)
will let the function infer missing 'domain' and 'measurementScale' values for attributes column.
Should be in same order as attributeNames in the attributes table, or be a named list with names corresponding to attributeNames
in the attributes table.}

\item{missingValues}{optional, a table with missing value code-deinition pairs; see details}
}
\value{
an eml "attributeList" object
}
\description{
set_attributes
}
\details{
The attributes data frame must use only the recognized column
headers shown here.  The attributes data frame must contain columns for required metadata.
These are:

\strong{For all data:}

\itemize{
  \item attributeName (required, free text field)

  \item attributeDefinition (required, free text field)

  \item measurementScale (required, "nominal", "ordinal", "ratio", "interval", or "dateTime",
 case sensitive) but it can be inferred from col_classes.
 
  \item domain (required, "numericDomain", "textDomain", "enumeratedDomain", or "dateTimeDomain",
 case sensitive) but it can be inferred from col_classes.
}

\strong{For numeric (ratio or interval) data:}
\itemize{
  \item unit (required). Unitless values should use "dimensionless" as the unit.
}

\strong{For character (textDomain) data:}
\itemize{
  \item definition (required)
}

\strong{For dateTime data:}
\itemize{
  \item formatString (required)
}

Other optional allowed columns in the attributes table are:
source, pattern, precision, numberType, missingValueCode, missingValueCodeExplanation,
attributeLabel, storageType, minimum, maximum

The \strong{factors} data frame, required for attributes in an enumerated domain, must use only the
 following recognized column headers:
\itemize{
  \item attributeName (required)
  \item code (required)
  \item definition (required)
}

The \strong{missingValues} data frame, optional, can be used in the case that multiple missing value codes
need to be set for the same attribute. This table must contain the following recognized column
headers.
\itemize{
  \item attributeName (required)
  \item code (required)
  \item definition (required)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_attributes.R
\name{get_numberType}
\alias{get_numberType}
\title{Get EML numberType}
\usage{
get_numberType(values)
}
\arguments{
\item{values}{(numeric/character) a vector of values, if vector is non-numeric will return NA}
}
\value{
the numberType of \code{values} (either 'real', 'integer', 'whole', or 'natural').
}
\description{
returns the EML numberType (either 'real', 'integer', 'whole', or 'natural') of input values
}
\examples{
\dontrun{
# To get numberType for each column in a data.frame:

unlist(lapply(df, function(x) get_numberType(x)))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_coverage.R
\name{set_taxonomicCoverage}
\alias{set_taxonomicCoverage}
\title{set_taxonomicCoverage}
\usage{
set_taxonomicCoverage(sci_names, expand = FALSE, db = "itis")
}
\arguments{
\item{sci_names}{string (space separated) or list or data frame of scientific names for species covered.}

\item{expand}{Set to TRUE to use `[taxadb]` to expand sci_names into full taxonomic classifications}

\item{db}{The taxonomic database to query (when expand is set to \code{TRUE}). See `[taxadb::filter_name]` for valid options. Defaults to 'itis'.}
}
\value{
a taxonomicCoverage object for EML
}
\description{
set_taxonomicCoverage
}
\details{
Turn a data.frame or a list of scientific names into a taxonomicCoverage block
sci_names can be a space-separated character string or a data frame with column names as rank name
or a list of user-defined taxonomicClassification
}
\note{
If "sci_names" is a data frame, column names of the data frame are rank names.
For user-defined "sci_names", users must make sure that the order of rank names
they specify is from high to low.
Ex. "Kingdom","Phylum","Class","Order","Family","Genus","Species","Common"
EML permits any rank names provided they go in descending order.
}
\examples{

taxon_coverage <- set_taxonomicCoverage("Macrocystis pyrifera")

sci_names <- data.frame(
  Kingdom = "Plantae",
  Phylum = "Phaeophyta",
  Class = "Phaeophyceae",
  Order = "Laminariales",
  Family = "Lessoniaceae",
  Genus = "Macrocystis",
  specificEpithet = "pyrifera"
)
taxon_coverage <- set_taxonomicCoverage(sci_names)

\donttest{ # Examples that may take > 5s

## use a list of lists for multiple species
sci_names <- list(list(
  Kingdom = "Plantae",
  Phylum = "Phaeophyta",
  Class = "Phaeophyceae",
  Order = "Laminariales",
  Family = "Lessoniaceae",
  Genus = "Macrocystis",
  specificEpithet = "pyrifera"
))
set_taxonomicCoverage(sci_names)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_responsibleParty.R
\name{set_responsibleParty}
\alias{set_responsibleParty}
\title{set_responsibleParty}
\usage{
set_responsibleParty(
  givenName = NULL,
  surName = NULL,
  organizationName = NULL,
  positionName = NULL,
  address = NULL,
  phone = NULL,
  electronicMailAddress = NULL,
  onlineUrl = NULL,
  userId = NULL,
  id = NULL,
  email = NULL
)
}
\arguments{
\item{givenName}{individual's given names (list or vector for multiple names).  OR a person object.}

\item{surName}{individual name}

\item{organizationName}{if party is an organization instead of an individual, name for the org}

\item{positionName}{individual's position, i.e. "Researcher", "Graduate Student", "Professor"}

\item{address}{address object, see `eml$address` to build an address object}

\item{phone}{individual or organization phone number}

\item{electronicMailAddress}{email address (alternatively, can use 'email' argument)}

\item{onlineUrl}{a URL to the homepage of the individual or organization}

\item{userId}{the user's ID, usually within a particular system (KNB, DataONE)}

\item{id}{Identifier for this block, ideally an ORCID id (optional)}

\item{email}{alias for electronicMailAddress}
}
\value{
A emld object for any responsibleParty (e.g. creator, contact, etc)
}
\description{
set_responsibleParty
}
\examples{
carl <- set_responsibleParty(as.person("Carl Boettiger <cboettig@ropensci.org>"))
matt <- set_responsibleParty("Matthew", "Jones", email = "mbjones@nceas.ucsb.edu")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_methods.R
\name{set_methods}
\alias{set_methods}
\title{set_methods}
\usage{
set_methods(
  methods_file,
  instrumentation = character(),
  software = NULL,
  sampling_file = NULL,
  sampling_coverage = NULL,
  sampling_citation = NULL,
  qualityControl_file = NULL
)
}
\arguments{
\item{methods_file}{Path to a file (markdown or .docx) containing a description of the methods used}

\item{instrumentation}{optional, text describing instrumentation used in methods}

\item{software}{optional, an EML software node describing software used in methods}

\item{sampling_file}{optional, Path to a file (.md or .docx) describing sampling method}

\item{sampling_coverage}{optional, coverage node for methods, e.g. set_coverage()}

\item{sampling_citation}{optional, a citation element describing the sampling protocol}

\item{qualityControl_file}{optional, path to a file (.md or .docx) describing quality control methods}
}
\value{
A methods object
}
\description{
set_methods
}
\examples{
\donttest{
f <- system.file("examples/hf205-methods.md", package = "EML")
set_methods(methods_file = f)

## Can also import from methods written in a .docx MS Word file.
f <- system.file("examples/hf205-methods.docx", package = "EML")
set_methods(methods_file = f)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_attributes.R
\name{htmlwidgets_attributes}
\alias{htmlwidgets_attributes}
\title{Launch attributes htmlwidget}
\usage{
htmlwidgets_attributes(df, type = NULL)
}
\arguments{
\item{df}{(data.frame) the data.frame of data that needs an attribute table}

\item{type}{(character) either "attributes", "units", or "factors"}
}
\description{
Used to call handsontable html widget to build attributes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_attributes.R
\name{is_standardUnit}
\alias{is_standardUnit}
\title{is_standardUnit}
\usage{
is_standardUnit(x)
}
\arguments{
\item{x}{name of unit to check}
}
\value{
TRUE if unit is exact match to the id of
 a unit in the Standard Units Table, FALSE otherwise.
}
\description{
is_standardUnit
}
\examples{
is_standardUnit("amperePerMeter") # TRUE
is_standardUnit("speciesPerSquareMeter") # FALSE
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_attributes.R
\name{build_factors}
\alias{build_factors}
\title{build factor table}
\usage{
build_factors(att_table, data)
}
\arguments{
\item{att_table}{(data.frame) input attributes table}

\item{data}{(data.frame) input data}
}
\description{
builds factor table for shiny app
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_attributes.R
\name{build_units_table}
\alias{build_units_table}
\title{build units table}
\usage{
build_units_table(in_units, eml_units)
}
\arguments{
\item{in_units}{input units}

\item{eml_units}{eml units}
}
\description{
builds unit table for shiny app
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_attributes.R
\name{shiny_attributes}
\alias{shiny_attributes}
\title{Create/Edit EML attributes}
\usage{
shiny_attributes(data = NULL, attributes = NULL)
}
\arguments{
\item{data}{(data.frame) the data.frame of data that needs an attribute table}

\item{attributes}{(data.frame) an existing attributes table}
}
\description{
Create/edit EML attributes, custom units, and factors in a shiny environment.
}
\details{
Attributes can be created from scratch using \code{shiny_attributes()}.
Or an existing attribute table can be edited using \code{shiny_attributes(NULL, attributes)}.
Or new attributes can be created from a data table using \code{shiny_attributes(data, NULL)}.
If attributes are created from a data table, fields such as `attributeName` and `numberType` will be automatically
completed based on the attributes within the data table.
If both existing attributes and data table are entered (i.e. \code{shiny_attributes(data, attributes)}),
any automatically generated fields based attributes within the data table **will not** override any non-empty fields in the
entered attributes
}
\examples{
\dontrun{
# from scratch
out <- shiny_attributes(NULL, NULL)

# from data
data <- iris
out <- shiny_attributes(data, NULL)

# from exisiting attributes
file <- system.file("tests", emld::eml_version(),
  "eml-datasetWithAttributelevelMethods.xml",
  package = "emld"
)
eml <- read_eml(file)
x <- eml$dataset$dataTable$attributeList
df <- get_attributes(x, eml)
out <- shiny_attributes(NULL, df$attributes)

# from attributes and data
out <- shiny_attributes(data, df$attributes)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_person.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{as_emld}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{emld}{\code{\link[emld]{as_emld}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_attributes.R
\name{get_attributes}
\alias{get_attributes}
\title{get_attributes}
\usage{
get_attributes(x, eml = NULL)
}
\arguments{
\item{x}{an "attributeList" element from an emld object}

\item{eml}{The full eml document, needed only if <references> outside of attributes must be resolved.}
}
\value{
a data frame whose rows are the attributes (names of each column in the data file)
and whose columns describe metadata about those attributes.  By default separate tables
are given for each type
}
\description{
get_attributes
}
\details{
EML metadata can use "references" elements which allow one attribute to use metadata
declared elsewhere in the document.  This function will automatically resolve these references
and thus infer the correct metadata.
}
\examples{
f <- system.file("tests", emld::eml_version(), 
  "eml-datasetWithAttributelevelMethods.xml", package = "emld")
eml <- read_eml(f)
get_attributes(eml$dataset$dataTable$attributeList)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_coverage.R
\name{set_coverage}
\alias{set_coverage}
\title{set_coverage}
\usage{
set_coverage(
  beginDate = character(),
  endDate = character(),
  date = character(),
  sci_names = character(),
  geographicDescription = character(),
  westBoundingCoordinate = numeric(),
  eastBoundingCoordinate = numeric(),
  northBoundingCoordinate = numeric(),
  southBoundingCoordinate = numeric(),
  altitudeMinimum = numeric(),
  altitudeMaximum = numeric(),
  altitudeUnits = character()
)
}
\arguments{
\item{beginDate}{Starting date for temporal coverage range.}

\item{endDate}{End date for temporal coverage range}

\item{date}{give a single date, or vector of single dates covered (instead of beginDate and endDate)}

\item{sci_names}{string (space separated) or list or data frame of scientific names for species covered.  See details}

\item{geographicDescription}{text string describing the geographic location}

\item{westBoundingCoordinate}{Decimal longitude for west edge bounding box}

\item{eastBoundingCoordinate}{Decimal longitude for east edge bounding box}

\item{northBoundingCoordinate}{Decimal latitude value for north of bounding box}

\item{southBoundingCoordinate}{Decimal latitude value for south edge of bounding box}

\item{altitudeMinimum}{minimum altitude covered by the data (optional)}

\item{altitudeMaximum}{maximum altitude covered by the data (optional)}

\item{altitudeUnits}{name of the units used to measure altitude, if given}
}
\value{
a coverage object for EML
}
\description{
set_coverage
}
\details{
set_coverage provides a simple and concise way to specify most common temporal,
taxonomic, and geographic coverage metadata. For certain studies this will not be
well suited, and users will need the more flexible but more verbose construction using
"new()" methods; for instance, to specify temporal coverage in geological epoch instead
of calendar dates, or to specify taxonomic coverage in terms of other ranks or identifiers.
}
\note{
If "sci_names" is a data frame, column names of the data frame are rank names.
For user-defined "sci_names", users must make sure that the order of rank names
they specify is from high to low.
Ex. "Kingdom","Phylum","Class","Order","Family","Genus","Species","Common"
}
\examples{
coverage <-
  set_coverage(
    begin = "2012-06-01", end = "2013-12-31",
    sci_names = "Sarracenia purpurea",
    geographicDescription = "California coast, down through Baja, Mexico",
    west = -122.44, east = -117.15,
    north = 37.38, south = 30.00
  )
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_physical.R
\name{set_physical}
\alias{set_physical}
\title{set_physical}
\usage{
set_physical(
  objectName,
  id = character(),
  numHeaderLines = character(),
  numFooterLines = character(),
  recordDelimiter = detect_delim(objectName),
  fieldDelimiter = ",",
  collapseDelimiters = logical(),
  literalCharacter = character(),
  quoteCharacter = character(),
  attributeOrientation = "column",
  size = NULL,
  sizeUnit = "bytes",
  authentication = NULL,
  authMethod = NULL,
  characterEncoding = character(),
  encodingMethod = character(),
  compressionMethod = character(),
  url = character()
)
}
\arguments{
\item{objectName}{name for the object, usually a filename like "hf205-1.csv"}

\item{id}{optional, an id value for the <physical> element in EML, for use in referencing}

\item{numHeaderLines}{Number of header lines preceding data. Lines are determined by the physicalLineDelimiter, or if it is absent, by the recordDelimiter. This value indicated the number of header lines that should be skipped before starting to parse the data.}

\item{numFooterLines}{Number of footer lines following data. Lines are determined by the physicalLineDelimiter, or if it is absent, by the recordDelimiter. This value indicated the number of footer lines that should be skipped after parsing the data. If this value is omitted, parsers should assume the data continues to the end of the data stream.}

\item{recordDelimiter}{This element specifies the record delimiter character when the format is text. The record delimiter is usually a linefeed (\\n) on UNIX, a carriage return (\\r) on MacOS, or both (\\r\\n) on Windows/DOS. Multiline records are usually delimited with two line ending characters, for example on UNIX it would be two linefeed characters (\\n\\n). As record delimiters are often non-printing characters, one can use either the special value "\\n" to represent a linefeed (ASCII 0x0a) and "\\r" to represent a carriage return (ASCII 0x0d). Alternatively, one can use the hex value to represent character values (e.g., 0x0a).}

\item{fieldDelimiter}{"," character by default (for csv files). This element specifies a character to be used in the object for indicating the ending column for an attribute. The delimiter character itself is not part of the attribute value, but rather is present in the column following the last character of the value. Typical delimiter characters include commas, tabs, spaces, and semicolons. The only time the fieldDelimiter character is not interpreted as a delimiter is if it is contained in a quoted string (see quoteCharacter) or is immediately preceded by a literalCharacter. Non-printable quote characters can be provided as their hex values, and for tab characters by its ASCII string "\\t". Processors should assume that the field starts in the column following the previous field if the previous field was fixed, or in the column following the delimiter from the previous field if the previous field was delimited.}

\item{collapseDelimiters}{The collapseDelimiters element specifies whether sequential delimiters should be treated as a single delimiter or multiple delimiters. An example is when a space delimiter is used; often there may be several repeated spaces that should be treated as a single delimiter, but not always. The valid values are yes or no. If it is set to yes, then consecutive delimiters will be collapsed to one. If set to no or absent, then consecutive delimiters will be treated as separate delimiters. Default behavior is no; hence, consecutive delimiters will be treated as separate delimiters, by default.}

\item{literalCharacter}{This element specifies a character to be used for escaping special character values so that they are treated as literal values. This allows "escaping" for special characters like quotes, commas, and spaces when they are intended to be used in an attribute value rather than being intended as a delimiter. The literalCharacter is typically a \\.}

\item{quoteCharacter}{This element specifies a character to be used in the object for quoting values so that field delimiters can be used within the value. This basically allows delimiter "escaping". The quoteChacter is typically a " or '. When a processor encounters a quote character, it should not interpret any following characters as a delimiter until a matching quote character has been encountered (i.e., quotes come in pairs). It is an error to not provide a closing quote before the record ends. Non-printable quote characters can be provided as their hex values.}

\item{attributeOrientation}{Specifies whether the attributes described in the physical stream are found in columns or rows. The valid values are column or row. If set to 'column', then the attributes are in columns. If set to 'row', then the attributes are in rows. Row orientation is rare.}

\item{size}{This element contains information of the physical size of the entity, by default represented in bytes unless the sizeUnit attribute is provided to change the units.}

\item{sizeUnit}{the unit in which size is measured; default is 'bytes'}

\item{authentication}{This element describes authentication procedures or techniques, typically by giving a checksum value for the object. The method used to compute the authentication value (e.g., MD5) is listed in the method attribute.}

\item{authMethod}{the method for authentication checksum, e.g. MD5}

\item{characterEncoding}{This element contains the name of the character encoding. This is typically ASCII or UTF-8, or one of the other common encodings.}

\item{encodingMethod}{This element lists a encoding method used to encode the object, such as base64, BinHex.}

\item{compressionMethod}{This element lists a compression method used to compress the object, such as zip, compress, etc. Compression and encoding methods must be listed in the order in which they were applied, so that decompression and decoding should occur in the reverse order of the listing. For example, if a file is compressed using zip and then encoded using MIME base64, the compression method would be listed first and the encoding method second.}

\item{url}{optional. The complete url from which the data file can be downloaded, if possible.}
}
\value{
an EML physical object, such as used in a dataTable element to define the format of the data file.
}
\description{
Will calculate the file size, checksum, and checksum authentication method
algorithm automatically if the argument \code{objectName} is a file that exists.
}
\examples{
set_physical("hf205-01-TPexp1.csv")
# FIXME set recordDelimiter based on user's system?
# FIXME richer distribution options? use set_distribution at top level?
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_software.R
\name{set_software}
\alias{set_software}
\title{set_software}
\usage{
set_software(codemeta)
}
\arguments{
\item{codemeta}{codemeta object, see examples}
}
\value{
an eml software element
}
\description{
set_software
}
\examples{
cm <- jsonlite::read_json(system.file("extdata/codemeta.json", package = "EML"))
software <- set_software(cm)
my_eml <- eml$eml(packageId = "eml-1.2", system = "knb", software = software)

# write_eml(my_eml, "test.xml")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_eml.R
\name{write_eml}
\alias{write_eml}
\title{write_eml}
\usage{
write_eml(eml, file, namespaces = NULL, ns = "eml", ...)
}
\arguments{
\item{eml}{an emld class object}

\item{file}{file name to write XML.}

\item{namespaces}{named character vector of additional XML namespaces to use.}

\item{ns}{root namespace abbreviation}

\item{...}{additional arguments to \code{\link{write_xml}}}
}
\value{
If file is not specified, the result is a character string containing
   the resulting XML content. Otherwise return silently.
}
\description{
write_eml
}
\examples{
f <- system.file("extdata", "example.xml", package = "emld")
eml <- read_eml(f)
write_eml(eml, "test.xml")
eml_validate("test.xml")
unlink("test.xml") # clean up
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/template_constructor.R
\docType{data}
\name{eml}
\alias{eml}
\title{eml}
\format{
A list with constructor functions
}
\description{
eml
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_unit_id.R
\name{get_unit_id}
\alias{get_unit_id}
\title{get_unit_id}
\usage{
get_unit_id(input_units, eml_version = emld::eml_version())
}
\arguments{
\item{input_units}{(character|vector) input units that needs
valid EML unit ids}

\item{eml_version}{(character) the eml schema version desired
(there is a change in the way eml units are named from eml-2.1.1
 to eml-2.2.0)}
}
\value{
(character) A valid EML unit id. If no valid EML unit id
 can be found, the function will output a warning, along with a
 preformatted custom unit id.
}
\description{
A function to assist with getting valid EML unit ids
 (see examples). Warning: ensure returned unit is equivalent to
 input unit (for example "pH" will return "picohenry" which may
 or may not be equivalent to the input unit "pH").
}
\examples{
\dontrun{
# The following all return the same id
get_unit_id("kilometersPerSquareSecond")
get_unit_id("kilometerPerSecondSquared")
get_unit_id("Kilometers per seconds squared")
get_unit_id("km/s^2")
get_unit_id("km s-2")
get_unit_id("s-2 /     kilometers-1") # this works but is not advised
}
}
