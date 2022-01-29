
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
