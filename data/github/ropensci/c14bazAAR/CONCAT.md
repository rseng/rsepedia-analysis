---
title: 'c14bazAAR: An R package for downloading and preparing C14 dates from different source databases'
tags:
  - R
  - archaeology
  - radiocarbon dates
  - data retrieval
authors:
  - name: Clemens Schmid
    orcid: 0000-0003-3448-5715
    affiliation: 1
  - name: Dirk Seidensticker
    orcid: 0000-0002-8155-7702
    affiliation: 2
  - name: Martin Hinz
    orcid: 0000-0002-9904-6548
    affiliation: 3
affiliations:
  - name: Max Planck Institute for the Science of Human History
    index: 1
  - name: Faculty of Arts and Philosophy, Department of Languages and Cultures, Ghent University
    index: 2
  - name: Institute of Archaeological Sciences, University of Bern
    index: 3
date: 17 November 2019
bibliography: paper.bib
---

# Background

Radiocarbon dating is one of the most important methods for absolute and relative chronological reconstruction of cultural development in prehistoric archaeology [@Taylor:2016]. This is true for the intrasite level, for regional comparisons and also in cases, where processes that have a large spatial and temporal reach are to be investigated. Prominent examples for the latter include the 'neolithization' of Europe or the 'Bantu expansion' in sub-Saharan Africa. They have been extensively analysed with radiocarbon data (e.g. @Ammerman:1971, @Garcin:2018, @Jerardino:2014, @Lemmen:2011, @Oslisly:2013, @Pinhasi:2005, @Russell:2014, @Weninger:2009).

Data selection for such models is complex and requires a thorough understanding of the archaeological questions. Most of the time it is not sufficient to include every date that vaguely fits into the context. Some dates have to be deliberately omitted. In order to make this important process of data filtering as transparent and reproducible as possible, the criteria for selection and especially the original data set must be generally accessible and well contextualized. Otherwise, peers can not evaluate the results in a meaningful way.

Fortunately, many archaeological institutions and individual authors are sharing their radiocarbon collections online (e.g. @Hinz:2012, @Kneisel:2013, @Mustaphi:2016, @Seidensticker:2016) -- some of them are archives with a long tradition of quality control and maintenance. Also, the boom of the Open Data movement in recent years has led to an increase of publications with raw data in archaeology (e.g. @Palmisano:2017). These collections are an important archive for future research questions.

However, the entire data basis is currently highly decentralized and lacks basic standardisation. That results in an effective loss of the possible added value that could result from the intersection of data sets in terms of searchability, error checking and further analysis. The creation of a world-wide and centralised database of radiocarbon dates which could solve these issues is not to be expected for the near future.

# Code Summary

``c14bazAAR`` is an R package that attempts to tackle the problem at hand by providing an independent interface to access radiocarbon data and make it available for a reproducible research process: from modelling to publication to scientific discourse. It queries openly available ^14^C data archives, but not those behind pay- or login-walls.

The package includes download functions (accessible with the main interface `c14bazAAR::get_c14data()`) that -- first of all -- acquire databases from different sources online. They then reduce the tables to a set of common variables and store them in a dedicated R S3 class: `c14_date_list`. The `c14_date_list` is based on `tibble::tibble` to integrate well into the R [tidyverse](https://www.tidyverse.org/) ecosystem. It also establishes standardised data types for the most important variables usually defined to describe radiocarbon data.

Beyond the download functions, ``c14bazAAR`` contains a multitude of useful helpers that can be applied to objects of class `c14_date_list`. These include methods for bulk calibration of radiocarbon dates with the Bchron R package [@Haslett:2008], removal of duplicates, estimation of coordinate precision, or conversion to other useful R data types (e.g. `sf::sf` [@Pebesma:2018]). For the classification of sample material ``c14bazAAR`` provides a manually curated reference list that maps the inconsistent attributions in the source databases to a standardized set of material classes. Such a reference list exists as well to fix the country attribution value of dates -- which is especially important in case of missing coordinate information. Methods to determine the source country based on coordinates fail on such dates.

``c14bazAAR`` was already used for data acquisition and preparation in at least one research paper: @Schmid:2019.

# Acknowledgements

The package got valuable code input from several members of the [ISAAKiel group](https://isaakiel.github.io) (Initiative for Statistical Analysis in Archaeology Kiel), most notably: Daniel Knitter, David Matzig, Wolfgang Hamer, Kay Schmütz and Nils Müller-Scheeßel.

# References
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
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/<p align="center">
  <img src="man/figures/logo.png" width = 200>
</p>

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/c14bazAAR/actions/workflows/check-release.yaml/badge.svg)](https://github.com/ropensci/c14bazAAR/actions/workflows/check-release.yaml)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/c14bazAAR/master.svg)](https://codecov.io/github/ropensci/c14bazAAR?branch=master)
[![license](https://img.shields.io/badge/license-GPL%202-B50B82.svg)](https://www.r-project.org/Licenses/GPL-2)
[![DOI](https://img.shields.io/badge/DOI-10.17605%2FOSF.IO%2F3DS6A-blue)](https://doi.org/10.17605/OSF.IO/3DS6A)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.01914/status.svg)](https://doi.org/10.21105/joss.01914)

c14bazAAR is an R package to query different openly accessible radiocarbon date databases. It allows basic data cleaning, calibration and merging. If you're not familiar with R other tools (such as [GoGet](https://www.ibercrono.org/goget/index.php)) to search for radiocarbon dates might be better suited for your needs.

- [**Installation**](#installation)
- [**How to use**](#how-to-use) ([Download](#download), [Calibration](#calibration), [Country attribution](#country-attribution), [Duplicates](#duplicates), [Conversion](#conversion), [Technical functions](#technical-functions), [Plotting and visualization](#plotting-radiocarbon-data), [Interaction with other radiocarbon data packages](#other-radiocarbon-packages))
- [**Databases**](#databases)
- [**Contributing**](#contributing) ([Adding database getter functions](#adding-database-getter-functions), [Pre-submision testing](#pre-submision-testing), [Versioning](#versioning))
- [**Citation**](#citation)
- [**License**](#license)

If you want to use data downloaded with c14bazAAR for your research, you have to cite the respective source databases. Most databases have a preferred way of citation that also may change over time with new versions and publications. Please check the [relevant homepages](#databases) to find out more. The output of c14bazAAR does not contain the full citations of the individual dates, but only a short reference tag. For further information you have to consult the source databases.

### Installation

We recommend to install the stable version from the [R-universe](https://r-universe.dev/) repository of [rOpenSci](https://ropensci.org/) with the following command (in your R console):

```
install.packages("c14bazAAR", repos = c(ropensci = "https://ropensci.r-universe.dev"))
```

The development version can be installed from github with the following command (in your R console):

```
if(!require('remotes')) install.packages('remotes')
remotes::install_github("ropensci/c14bazAAR")
```

Both versions are up-to-date and include all databases and features. Installing the development version on Windows requires the toolchain bundle [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

The package needs a lot of other packages -- many of them only necessary for specific tasks. Functions that require certain packages you don't have installed yet will stop and ask you to enable them. Please do so with [`install.packages()`](https://www.r-bloggers.com/installing-r-packages/) to download and install the respective packages from CRAN.

### How to use

The package contains a set of getter functions (see below) to query the databases. Thereby not every available variable from every archive is downloaded. Instead c14bazAAR focuses on a [selection](https://github.com/ropensci/c14bazAAR/blob/master/data-raw/variable_definition.csv) of the most important and most common variables to achieve a certain degree of standardization. The downloaded dates are stored in the custom S3 class `c14_date_list` which acts as a wrapper around the [tibble](https://tibble.tidyverse.org/) class and provides specific class methods.

A workflow to download and prepare all dates could look like this:

```
library(c14bazAAR)
library(magrittr)

get_c14data("all") %>%
  remove_duplicates() %>%
  calibrate() %>%
  determine_country_by_coordinate()
```

It takes quite some time to run all of this and it's probably not necessary for your use case. Here's a list of the main tasks c14bazAAR can handle. That allows you to pick what you need:

#### Download

c14bazAAR contains a growing selection of getter functions to download radiocarbon date databases. [Here's](#databases) a list of all available databases. You can download all dates at once with [`get_c14data("all")`](https://github.com/ropensci/c14bazAAR/blob/master/R/get_c14data.R). The getters download the data, adjust the variable selection according to a defined [variable key](https://github.com/ropensci/c14bazAAR/blob/master/data-raw/variable_reference.csv) and transform the resulting list into a `c14_date_list`. 

See `?get_c14data` for more information.

```
x <- get_c14data("all")
```

#### Calibration

The [`calibrate()`](https://github.com/ropensci/c14bazAAR/blob/master/R/c14_date_list_calibrate.R) function calibrates all valid dates in a `c14_date_list` individually with [`Bchron::BchronCalibrate()`](https://github.com/andrewcparnell/Bchron/blob/master/R/BchronCalibrate.R). It provides two different types of output: calprobdistr and calrange.

See `?calibrate` for more information.

```
x %>% calibrate()
```

#### Country attribution

Filtering 14C dates by country is useful for spatial filtering. Most databases provide the variable country, but they don't rely on a unified naming convention and therefore use various terms to represent the same entity. The function [`determine_country_by_coordinate()`](https://github.com/ropensci/c14bazAAR/blob/master/R/c14_date_list_spatial_determine_country_by_coordinate.R) determines the country a date is coming from by intersecting its spatial coordinates with polygons from [`rworldxtra::countriesHigh`](https://github.com/AndySouth/rworldxtra).

See `?country_attribution` for more information.

```
x %>% determine_country_by_coordinate()
```

#### Duplicates

Some of the source databases already contain duplicated dates and for sure you'll have some if you combine different databases. As a result of the long history of these archives, which includes even mutual absorption, duplicates make up a significant proportion of combined datasets. It's not trivial to find and deal with theses duplicates, because they are not exactly identical between databases: Sometimes they are linked to conflicting and mutually exclusive context information.

For an automatic search and removal based on identical lab numbers we wrote [`remove_duplicates()`](https://github.com/ropensci/c14bazAAR/blob/master/R/c14_date_list_duplicates_remove.R). This functions offers several options on how exactly duplicates should be treated.

If you call `remove_duplicates()` with the option `mark_only = TRUE` then no data is removed, but you can inspect the duplicate groups identified.

See `?duplicates` for more information.

```
x %>% remove_duplicates()
```

#### Conversion

A c14_date_list can be directly converted to other R data structures. So far only [`as.sf()`](https://github.com/ropensci/c14bazAAR/blob/master/R/c14_date_list_convert.R) is implemented. The sf package provides great tools to manipulate and plot spatial vector data. This simplifies certain spatial operations with the date point cloud.

See `?as.sf` for more information.

```
x %>% as.sf()
```

#### Technical functions

c14_date_lists are constructed with [`as.c14_date_list`](https://github.com/ropensci/c14bazAAR/blob/master/R/c14_date_list_basic.R). This function takes data.frames or tibbles and adds the c14_date_list class tag. It also calls [`order_variables()`](https://github.com/ropensci/c14bazAAR/blob/master/R/c14_date_list_order_variables.R) to establish a certain variable order and [`enforce_types()`](https://github.com/ropensci/c14bazAAR/blob/master/R/c14_date_list_enforce_types.R) which converts all variables to the correct data type. There are custom `print()`, `format()` and `plot()` methods for c14_date_lists.

The [`fuse()`](https://github.com/ropensci/c14bazAAR/blob/master/R/c14_date_list_fuse.R) function allows to rowbind multiple c14_date_lists.

See `?as.c14_date_list` and `?fuse`.

```
x1 <- data.frame(
  c14age = 2000,
  c14std = 30
) %>% as.c14_date_list()

x2 <- fuse(x1, x1)
```

#### Plotting radiocarbon data

c14bazAAR only provides a very basic `plot` function for `c14_date_list`s. The [simple plotting vignette](https://github.com/ropensci/c14bazAAR/blob/master/vignettes/simple_plotting.Rmd) introduces some techniques to help you get started with more sophisticated visualization.

#### Other radiocarbon packages

There are several R packages that provide functions to calibrate, analyze or model radiocarbon dates: e.g. [oxcAAR](https://github.com/ISAAKiel/oxcAAR), [rcarbon](https://github.com/ahb108/rcarbon), [Bchron](https://github.com/andrewcparnell/Bchron)

They usually have a simple, vector based interface and you can use `c14_date_list` columns as input.

```
rcarbon::calibrate(x = x$c14age, error = x$c14std)
```

### Databases

<p align="center">
  <img src="man/figures/README_map_figure.jpeg" width = 100%>
</p>

To suggest other archives to be queried you can join the discussion [here](https://github.com/ropensci/c14bazAAR/issues/2).

* [`get_c14data("14cpalaeolithic")`](R/get_14cpalaeolithic.R) [**14cpalaeolithic**](https://ees.kuleuven.be/geography/projects/14c-palaeolithic): Radiocarbon Palaeolithic Europe Database V28, June 2021.
* [`get_c14data("14sea")`](R/get_14sea.R) [**14sea**](http://www.14sea.org/) 14C database for Southeast Europe and Anatolia (10,000–3000 calBC).
* [`get_c14data("adrac")`](R/get_adrac.R) [**adrac**](https://github.com/dirkseidensticker/aDRAC): Archives des datations radiocarbone d'Afrique centrale by Dirk Seidensticker.
* [`get_c14data("agrichange")`](R/get_agrichange.R) [**agrichange**](https://zenodo.org/record/4541470): Radiocarbon dates associated to Neolithic contexts (ca. 5900 – 2000 cal BC) from the northwestern Mediterranean Arch to the High Rhine area by [Héctor Martínez-Grau, Berta Morell-Rovira & Ferran Antolín](https://openarchaeologydata.metajnl.com/articles/10.5334/joad.72/) (2021).
* [`get_c14data("aida")`](R/get_aida.R) [**AIDA**](https://github.com/apalmisano82/AIDA): Archive of Italian radiocarbon DAtes Alessio Palmisano, Andrew Bevan, Alex Kabelindde, Neil Roberts & Stephen Shennan (2021).
* [`get_c14data("austarch")`](R/get_austarch.R) [**austarch**](https://archaeologydataservice.ac.uk/archives/view/austarch_na_2014/): A Database of 14C and Luminescence Ages from Archaeological Sites in Australia by [Alan N. Williams, Sean Ulm, Mike Smith, Jill Reid](https://intarch.ac.uk/journal/issue36/6/williams.html).
* [`get_c14data("bda")`](R/get_bda.R) [**BDA**](https://nakala.fr/10.34847/nkl.dde9fnm8): Base de Données Archéologiques by Thomas Perrin (1994).
* [`get_c14data("calpal")`](R/get_calpal.R) [**calpal**](https://uni-koeln.academia.edu/BernhardWeninger/CalPal): Radiocarbon Database of the CalPal software package by Bernhard Weninger. See [nevrome/CalPal-Database](https://github.com/nevrome/CalPal-Database) for an interface.
* [`get_c14data("caribbean")`](R/get_caribbean.R) [**caribbean**](https://github.com/philriris/caribbean-14C/): A compilation of 2147 anthropogenic radiocarbon (14C) dates for the Caribbean region from 504 sites across 57 islands by Phil Riris (2021).
* [`get_c14data("context")`](R/get_context.R) [**context**](http://context-database.uni-koeln.de/): Collection of radiocarbon dates from sites in the Near East and neighboring regions (20.000 - 5.000 calBC) by Utz Böhner and Daniel Schyle.
* [`get_c14data("eubar")`](R/get_eubar.R) [**eubar**](https://telearchaeology.org/eubar-c14-database/): A database of 14C measurements for the European Bronze Age by [Gacomo Capuzzo](https://telearchaeology.org/EUBAR/).
* [`get_c14data("euroevol")`](R/get_euroevol.R) [**euroevol**](https://discovery.ucl.ac.uk/1469811/): Cultural Evolution of Neolithic Europe Dataset by [Katie Manning, Sue Colledge, Enrico Crema, Stephen Shennan and Adrian Timpson](https://openarchaeologydata.metajnl.com/articles/10.5334/joad.40/).
* [`get_c14data("irdd")`](R/get_irdd.R) [**irdd**](https://sites.google.com/site/chapplearchaeology/irish-radiocarbon-dendrochronological-dates): [Robert M Chapple](https://doi.org/10.5281/zenodo.3367518)'s Catalogue of Radiocarbon Determinations & Dendrochronology Dates is a free-to-download resource for Irish archaeology.
* [`get_c14data("jomon")`](R/get_jomon.R) [**jomon**](https://github.com/ercrema/jomonPhasesAndPopulation): A multi-proxy inference of Jōmon population dynamics using bayesian phase models, residential data, and summed probability distribution of 14C dates [Enrico R. Crema and Ken'ichi Kobayashi](https://www.sciencedirect.com/science/article/pii/S0305440320300583) (2020).
* [`get_c14data("katsianis")`](R/get_katsianis.R) [**katsianis**](https://rdr.ucl.ac.uk/articles/Dataset_for_An_Aegean_history_and_archaeology_written_through_radiocarbon_dates/12489137/1): An Aegean History and Archaeology Written through Radiocarbon Dates [Markos Katsianis, Andrew Bevan, Giorgos Styliaras & Yannis Maniatis](https://openarchaeologydata.metajnl.com/articles/10.5334/joad.65/) (2020).
* [`get_c14data("kiteeastafrica")`](R/get_kiteeastafrica.R) [**kiteeastafrica**](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/NJLNRJ): Radiocarbon dates from eastern Africa in the CARD2.0 format by [Colin Courtney Mustaphi, Rob Marchant](https://www.openquaternary.com/articles/10.5334/oq.22/).
* [`get_c14data("medafricarbon")`](R/get_medafricarbon.R) [**MedAfriCarbon**](https://zenodo.org/record/3689716#.XnSHp4hKiUk): The MedAfriCarbon Radiocarbon Database and [Web Application](https://theia.arch.cam.ac.uk/MedAfriCarbon/). Archaeological Dynamics in Mediterranean Africa, ca. 9600–700 BC by [Giulio Lucarini, Toby Wilkinson, Enrico R. Crema, Augusto Palombini, Andrew Bevan and Cyprian Broodbank](https://openarchaeologydata.metajnl.com/articles/10.5334/joad.60/) (2020).
* [`get_c14data("nerd")`](R/get_nerd.R) [**NERD**](https://github.com/apalmisano82/NERD): Near East Radiocarbon Dates Alessio Palmisano, Andrew Bevan, Dan Lawrence & Stephen Shennan (2021).
* [`get_c14data("pacea")`](R/get_pacea.R) [**pacea**](http://www.paleoanthro.org/media/journal/content/PA20110001_S01.zip): PACEA Geo-Referenced Radiocarbon Database for the late Middle Paleolithic, Upper Paleolithic, and initial Holocene in Europe by [Francesco D'Errico, William E. Banks, Marian Vanhaeren, Véronique Laroulandie and Mathieu Langlais](http://www.paleoanthro.org/media/journal/content/PA20110001.pdf) (2011). 
* [`get_c14data("palmisano")`](R/get_palmisano.R) [**palmisano**](https://dx.doi.org/10.14324/000.ds.1575442): Regional Demographic Trends and Settlement Patterns in Central Italy: Archaeological Sites and Radiocarbon Dates by [Alessio Palmisano, Andrew Bevan and Stephen Shennan](https://openarchaeologydata.metajnl.com/articles/10.5334/joad.43/) (2018).
* [`get_c14data("radon")`](R/get_radon.R) [**radon**](https://radon.ufg.uni-kiel.de/): Central European and Scandinavian database of 14C dates for the Neolithic and Early Bronze Age by [Dirk Raetzel-Fabian, Martin Furholt, Martin Hinz, Johannes Müller, Christoph Rinne, Karl-Göran Sjögren und Hans-Peter Wotzka](https://www.jna.uni-kiel.de/index.php/jna/article/view/65).
* [`get_c14data("radonb")`](R/get_radonb.R) [**radonb**](https://radon-b.ufg.uni-kiel.de/): Database for European 14C dates for the Bronze and Early Iron Age by Jutta Kneisel, Martin Hinz, Christoph Rinne.
* [`get_c14data("rxpand")`](R/get_rxpand.R) [**rxpand**](https://github.com/jgregoriods/rxpand): Radiocarbon dates for the spread of farming and ceramics in tropical South America by Jonas Gregorio de Souza.
* [`get_c14data("sard")`](R/get_sard.R) [**sard**](https://github.com/emmaloftus/Southern-African-Radiocarbon-Database): Southern-African-Radiocarbon-Database by [Emma Loftus, Peter J. Mitchell & Christopher Bronk Ramsey](https://www.cambridge.org/core/journals/antiquity/article/abs/an-archaeological-radiocarbon-database-for-southern-africa/26FE99E995C4507015704D552CB0C196).

### Contributing

If you would like to contribute to this project, please start by reading our [Guide to Contributing](https://github.com/ropensci/c14bazAAR/blob/master/CONTRIBUTING.md). Please note that this project is released with a Contributor [Code of Conduct](https://github.com/ropensci/c14bazAAR/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

#### Adding database getter functions

If you want to add another radiocarbon database to c14bazAAR (maybe from the list [here](https://github.com/ropensci/c14bazAAR/issues/2)) you can follow this checklist to apply all the necessary changes to the package:

1. Add your database to the [variable_reference table](https://github.com/ropensci/c14bazAAR/blob/master/data-raw/variable_reference.csv) and map the database variables to the variables of c14bazAAR and other databases.
2. Write the getter function `get_[The Database Name]` in an own script file: **get_[the database name].R**. For the script file names we used a lowercase version of the database name. The getter functions have a standardized layout and always yield an object of the class `c14_date_list`. Please look at some of the available functions to get an idea how it is supposed to look like and which checks it has to include. Make sure not to store data outside of `tempdir()`. Some databases include non-radiocarbon dates: Make sure to filter them out -- c14bazAAR so far only works with radiocarbon dates.
3. Add the following roxygen2 tags above the function definition to include it in the package documentation.

```
#' @rdname db_getter_backend
#' @export
```

4. Update the package documentation with roxygen2.
5. Add the database url(s) to the [url_reference table](https://github.com/ropensci/c14bazAAR/blob/master/data-raw/db_info_table.csv) to make `get_db_url("[the database name]")` work.
6. Run the data-raw/data_prep.R script to update the data objects in the package. Only this enables the changes made in step 5. You should test your changes now by running the respective getter function.
7. Add the getter function your wrote in 2 to the functions vector in [`get_all_parser_functions()`](https://github.com/ropensci/c14bazAAR/blob/master/R/get_c14data.R#L128).
8. Document the addition of the new function in the NEWS.md file.
9. Add the new database to the list of *Currently available databases* in the DESCRIPTION file.
10. Add your function to the database list in the README file [here](https://github.com/ropensci/c14bazAAR#databases).
11. Update the README map figure by running the script [README_map_figure.R](https://github.com/ropensci/c14bazAAR/blob/master/figures/README_map_figure.R).

#### Pre-submision testing

Before submitting patches or new getter functions via a pull request, we ask you to check the following items:

1. The package works and all functions are usable
2. The package documentation is up-to-date and represents the functions correctly
3. The test coverage of the package functions is sufficient
4. `DESCRIPTION` is up-to-date with the latest version number and database list
5. `README.md` is up-to-date
6. `NEWS.md` is up-to-date and includes the latest changes
7. **Package checks ran and did not yield any ERRORS, WARNINGS or NOTES (or at least the NOTES are addressed in the cran-comments.md)**
	- **locally (`devtools::check()`)**
	- rhub (`devtools::check_rhub(email = ...)`)
	- winbuilder (`devtools::check_win_release(email = ....)` + `devtools::check_win_devel(email = ....)`)
8. Spellcheck with `devtools::spell_check()` ran and did yield not only false-positives
9. codemeta.json is up-to-date (can be updated with `codemetar::write_codemeta()`)
10. `inst/CITATION` is up-to-date
11. The package does not make external changes without explicit user permission. It does not write to the file system, change options, install packages, quit R, send information over the internet, open external software, etc.
12. No reverse dependencies break because of the new package version (`devtools::revdep_check()`)

Please make sure to run the tests listed above and pay special attention to the highlighted items.

#### Versioning

Version numbers (releases) follow the [semantic versioning schema](https://semver.org/) and consist of mayor and minor releases as well as patches.

* **x**.y.z: a **mayor** release will be made once an existing function is radically changed or removed and thus the package API is changed.
* x.**y**.z: a **minor** release contains new parsers and auxiliary functions.
* x.y.**z**: a **patch** updates existing parsers and functions.

### Citation

Schmid et al., (2019). c14bazAAR: An R package for downloading and preparing C14 dates from different source databases. Journal of Open Source Software, 4(43), 1914, https://doi.org/10.21105/joss.01914

```
@Article{Schmid2019,
  title = {{c14bazAAR}: An {R} package for downloading and preparing {C14} dates from different source databases},
  author = {Clemens Schmid and Dirk Seidensticker and Martin Hinz},
  journal = {Journal of Open Source Software},
  volume = {4},
  number = {43},
  pages = {1914},
  month = {nov},
  year = {2019},
  doi = {10.21105/joss.01914},
  url = {https://doi.org/10.21105/joss.01914},
}
```

### License

For the code in this project apply the terms and conditions of GNU GENERAL PUBLIC LICENSE Version 2. The source databases are published under different licenses.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# 3.2.1

- simplified many getter functions by replacing the explicit column classes definition with a general `colClasses = "character"`, which should do exactly the same (#153)

# 3.2.0

- added getter function for sard database: `get_sard`

# 3.1.0

- added getter function for aida database: `get_aida`
- adjusted getter function for the nerd database (`get_nerd`)

# 3.0.1

- updated `get_14cpalaeolithic` to be compatible with the latest respective database version v.28

# 3.0.0

- deprecated `classify_material`
- deprecated `fix_database_country_name`

# 2.4.2

- fixed an encoding issue in `get_pacea`

# 2.4.1

- some tiny fixes
- added an informative figure to the README

# 2.4.0

- added getter function for agrichange database: `get_agrichange`
- added getter function for caribbean database: `get_caribbean`

# 2.3.0

- deprecated `get_emedyd`, because the emedyd database was superseded by the nerd database

# 2.2.0

- added getter function for rxpand database: `get_rxpand`

# 2.1.0

- added getter function for NERD database: `get_nerd`
- added getter function for BDA database: `get_bda`

# 2.0.0

- the default lookup tables are now stored in the package and not downloaded from Github any more (#128)
- some changes in the internal functioning of the lookup methods (#128)
- introduced `inspect_lookup_country` and `inspect_lookup_material` which replace the automatic printing in `fix_database_country_name` and `classify_material` (#128)

# 1.3.3

- modified variable reference table layout (#126) 

# 1.3.2

- added submission ToDo list to README (#124)
- defined version update schema in README (#124)

# 1.3.1

- filtered out TL dates from aDRAC in parser function (#123)

# 1.3.0

## general changes
- added a CITATION file (see citation("c14bazAAR"))
- deprecated `mark_duplicates` to get rid of this extra step. You can get the same result now with `remove_duplicates(mark_only = TRUE)` (#100)
- deprecated `coordinate_precision`. The functionality was not essential and the calculated precision values probably frequently misleading. Beyond that the name was confusing (#96) (#106)
- deprecated `finalize_country_name`. This wrapper function was rather confusing and the functionality can be very easily be reimplemented if necessary (#96) (#106)
- renamed `standardize_country_name` to `fix_database_country_name` to make more clear what it does (#96) (#106)
- switched from openxlsx to readxl & writexl for handling xlsx files (#105) (#111)
- allowed for different calibration curves in calibrate (#118)

## new getter functions
- added getter function for MesoRAD database: `get_mesorad` (#112)
- added getter function for Katsianis et al. database: `get_katsianis` (#103)
- added getter function for PACEA database: `get_pacea` (#90)
- added getter function for 14C-Palaeolithic database: `get_14cpalaeolithic` (#90)
- added getter function for MedAfriCarbon database: `get_medafricarbon` (#95)
- added getter function for Jōmon population dynamics database: `get_jomon` (#95)
- added getter function for emedyd database: `get_emedyd` (#102)

## database updates
- updated the CalPal database from version 2017_04_27 to 2020_08_20 (#108)

## bugfixes
- `lwgeom::st_make_valid` was replaced by `sf::st_make_valid` (#99)
- enabled UTF-8 characters in country thesaurus (#96) (#104)
- `plot.c14_date_list` does not choke any more on c14_date_lists without coordinate columns (#112)
- fixed some entries in the country thesaurus

# 1.2.0

## general changes
- unified database names in all functions, tables, variables and documentation (#86)
- new logo and some layout changes in the README (#81)

## new features
- added a basic plot function for c14_date_lists (#82)
- added a basic write function for c14_date_lists: `write_c14()` (#84)
- added a version column that documents from which database version a certain date is pulled (#85)

# 1.1.0

## general changes
- [ROpenSci review](https://github.com/ropensci/software-review/issues/333)
- moved main development repository to github/ropensci/c14bazAAR (e0e6827f0381be04c50380eec277c01cad44ac7d)
- created [c14bazAAR project](https://doi.org/10.17605/OSF.IO/3DS6A) on the OSF platform with a DOI
- more work on an article for the Journal of Open Source Software (paper.md + paper.bib)
- changed citation to JOSS article after its publication

## new features
- new download interface as suggested by Enrico Crema in the ROpenSci review: `get_c14data()` (#76)
- replaced hard coded URLs with arguments to get helper functions (caabcb7b)

## new getter functions
- added getter function for irdd database: `get_irdd` (#79)

# 1.0.3

## general changes
- reformatted authors in DESCRIPTION and added ORCIDs (#72)
- added a [citation](https://github.com/ropensci/c14bazAAR#citation) section to the README
- added a [checklist](https://github.com/ropensci/c14bazAAR#adding-database-getter-functions) to the README on how to add new getter functions
- work on an article for the Journal of Open Source Software (paper.md + paper.bib)
- added a vignette with some plotting workflows
- created a completely artificial example dataset that replaces the sampled version

## new getter functions
- added getter function for Palmisano et al. database: `get_palmisano` (#59)
- added getter function for eubar database: `get_eubar` (#64)

## new features
- added new options for the deduplication function (see `?duplicates`) (#63)
- added an internal function `clean_labnr` to the `as.c14_date_list` workflow to fix certain syntactically wrong representations of lab numbers in several input databases as part of the downloading process (#61)
- better implementation of the `c14_date_list` as a subclass of tibble for a better integration of the subclass into the tidyverse (#67)

## bugfixes
- small file path construction fix in context getter function
- replaced some deprecated functions by other packages (dplyr::funs & tibble::as.tibble)
- replaced `RCurl::url.exists` with `httr::http_error` in `check_connection_to_url` (#68)
- fixed `as.sf` error that occurred when date lists with dates without coordinates were transformed

## removed objects
- data objects `c14bazAAR::country_thesaurus`, `c14bazAAR::material_thesaurus`, `c14bazAAR::variable-reference` have been removed from the package -- they are queried from [here](https://github.com/ropensci/c14bazAAR/tree/master/data-raw) anyway and it's not necessary to put them into the package
- some helper functions have been made internal
- .Rd files for unexported, internal objects have been removed (@noRd)

# 1.0.2

Release version
# Contributing

We love pull requests from everyone. By participating in this project, you
agree to abide by our [code of conduct](CONDUCT.md).

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free). If you are not familar with git and GitHub, take a look at <http://happygitwithr.com/> to get started.
* [Submit a post for your issue](https://github.com/ropensci/c14bazAAR/issues/), assuming one does not already exist.
  * Clearly describe your issue, including steps to reproduce when it is a bug, or some justification for a proposed improvement.
* [Fork](https://github.com/ropensci/c14bazAAR/#fork-destination-box) the repository on GitHub to make a copy of the repository on your account. Or use this line in your shell terminal:

    `git clone git@github.com:your-username/c14bazAAR.git`
    
## Making changes

* Edit the files, save often, and make commits of logical units, where each commit indicates one concept
* Follow a good [style guide](http://adv-r.had.co.nz/Style.html).
* Make sure you write [good commit messages](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html).
* Make sure you have added the necessary tests for your code changes.
* Run _all_ the tests using `devtools::check()` to assure nothing else was accidentally broken.
* If you need help or unsure about anything, post an update to [your issue](https://github.com/ropensci/c14bazAAR/issues/).

## Submitting your changes

Push to your fork and [submit a pull request](https://github.com/ropensci/c14bazAAR/compare/).

At this point you're waiting on us. We like to at least comment on pull requests
within a few days. We may suggest some changes or improvements or alternatives.
## Submission

This is the first submission of an updated version of the package.

## Test environments

* Ubuntu 18.04.3 LTS (Bionic Beaver), R 3.6.2
* Ubuntu Trusty 16.04.6 (on travis-ci), R 3.6.2
* local Windows 10 install, R 4.0.2
* win-builder (devel and release)

## R CMD check results

There were no NOTES, ERRORs or WARNINGs.
# Changes (19.06.2021)

- added additional terms for the agrichange database
- added additional terms for the caribbean database

# Changes (25.04.2020)

- added additional terms for the rxpand database

# Changes (30.08.2020)

- added additional terms from the new calpal version
- combined the categories NA, other and unknown to other/unknown - the distinction was not very clear
- changed some individual attributions

# Changes (21.8.2020)

- added additional terms from the Katsianis database
 
# Changes (18.08.2020)

- added values from emedyd

# Changes (13.3.2018)

- had to add all the values from 14SEA.
- decided to further reduce the amount of categories. 

## defined classes (13.3.2018)

* antler
* bone
* ceramics
* charcoal
* fabric
* food residues
* other
* plant remains
* shell
* soil
* teeth
* unknown
* wood

# Changes (21.2.2018)

* applied some of Martin's suggestions (see github 21.2.2018):
    * aragonite, ash, glass and quartz are now __"other"__. 
* Martin suggested to differentiate between human and animal bone. This would lead to three categories, since some bones are not described in more detail: e.g. "human bone", "animal bone", "bone"?
    * I don't know how important it is to differentiate between them. But in my opinion it would make the whole thing less easy to understand and work with? This is similar to the differentiation between "seawater shell" and "freshwater shell", which you, @MartinHinz, opted for to discard.
    *   for consistency reasons I suggest to either summarize all bones under the category "bone" and all shells under the category "shell", or to break them both down into smaller, more precise categories.
        * __categorized all kinds of shell as "shell"__

* "Baked clay from heart", "mud-brick wall", "Clayball from hearth"
    * as __"burnt clay"__
* "calcium carbonate" (eggs) as __"other"__
* charred fruit, reed, seeds, (vekoold), etc. are now __"charcoal"__

## defined classes (as of 21.2.2018)

* antler
* bark
* beeswax
* bone
* burnt clay
* ceramics
* charcoal
* fabric
* faeces
* food
* hair
* humic acid
* micro organism
* other
* plant remains
* pollen
* shell
* snail
* soil
* teeth
* tissue
* unknown
* wood


# Changes (19.2.2018)

* applied Clemens' suggestions (see github: 15.2.18: _"...Right now we have two different values for missing/bad/unclear information: NA and undefined. I suggest to replace undefined with other and NA with unknown. other: We have an information, but it's too strange/unique/specific to add it to the other established categories. unkown: We don't have any infomation."_)
* new categories: beeswax, quartz, glass, pollen, humic acid, 
* open questions: 
    * what happens with "burned" or "charred" samples? 
    * Is a marine snail (mangrove gastropod) different to a terrestrial snail (e.g. reservoir-effect)? 
    * can "leather" be categorized as "tissue"?
    * where to put: 
        * coal, graphite
        * clayball from hearth, mud-brick wall (ceramics?)
        * marine coral 


## defined classes (as of 19.2.2018)

* antler
* aragonite
* ash
* bark
* beeswax
* bone
* calcium carbonate
* ceramics
* charcoal
* fabric
* faeces
* food
* freshwater shell
* glass
* hair
* humic acid
* micro organism
* other
* plant remains
* pollen
* quartz
* shell
* snail
* soil
* teeth
* tissue
* unknown
* wood

# Changes (10.2.2018)

## defined classes (as of 10.2.2018)

* ash
* bark
* brain -> now tissue
* bone
* ceramics
* charcoal
* fabric
* faeces
* food
* hair
* horn
* plant remains
* shell
* soil
* teeth
* tissue
* NA
* wood
* aragonite

# Changes (3.5.2019)

* added additional terms for the Palmisano database

# Changes (25.06.2019)

* added additional terms for the EUBAR database 

# Changes (27.11.2019)

* added additional terms for the IRDD (Ireland) database

# Changes (23.03.2020)

* added additional terms from the MedAfriCarbon & Jomon database

## not yet defined, uncertain

* "Ash and charcoal" as ?
* "Ash and shell" as ? 
* "ashy soil" as "ash"?
* "soot" as "ash"?

* "Austomylitus rostratus shell" non-existant?
* "baked clay from hearth" and "mud-brick wall" as "ceramics"?
* "beeswax"
* can "bulk..." be classified? or is it NA because of uncertainties? (e.g. "bulk sediment")
* what about "burnt" or "charred" plant remains or soil? is it important that it's burnt?
* "calcium carbonte" as own category. Eggshell is classified as "calcium carbonate"
* "carbon" can be charcoal in french, or just carbon? how to classify?
* differentiation between marine and freshwater shell?
* can "glass" even be dated?
* "human" has been classified as "bone" but that's very imprecise, isn't it? Changed it to "NA"
* "pig" has been classified as "bone" -> incorrect?
* what about "human blood protein"? how can this even be 14C-dated?
* "humid acid" was defined as "soil". changed it to "undefined"
* "humic..." as new soil categorie?
* "leather" was defined as "fabric". V. imprecise. changed it to "undefined"
* "pollen" as own category?
* "stone" has been defined as "soil" -> incorrect?
* "turtle carapace" as "bone"




### db_info_table.csv

A database version and url lookup table.

- **db**: database name
- **version**: database version
- **url_num**: url number (some databases are spread over multiple files)
- **url**: file url where the database can be downloaded

### material_thesaurus.csv

A thesaurus for material classes. More info about the version history in material_thesaurus_comments.md.

- **cor**: fixed name
- **var**: variations

### variable_definition.csv

Definition of the variables in a c14_date_list

- **c14bazAAR_variable**: name of variable in c14bazAAR
- **type**: data type of variable in R
- **definition**: meaning of variable
- **source**: is the variable imported (databases) or generated (c14bazAAR)

### variable_reference.csv

Which variables in a c14_date_list equal the ones in the source databases?

- **c14bazAAR_variable**: name of variable in c14bazAAR
- **database**: source database
- **database_variable**: name of variable in the respective source database

