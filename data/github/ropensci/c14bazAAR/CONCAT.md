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

---
title: "Simple plotting options for radiocarbon dates in c14_date_lists"
author: "Clemens Schmid"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    fig_width: 7
    fig_height: 5
vignette: >
  %\VignetteIndexEntry{Simple plotting options for radiocarbon dates in c14_date_lists}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

This vignette shows some basic workflows to plot radiocarbon dates in `c14_date_list`s. This is only a short compilation to get you started.

So let's begin by loading the main packages for this endeavour: c14bazAAR and [ggplot2](https://CRAN.R-project.org/package=ggplot2). And of course [magrittr](https://CRAN.R-project.org/package=magrittr) to enable the pipe (`%>%`) operator. We will use functions from other packages as well, but address them individually with the `::` operator. 

```{r}
library(c14bazAAR)
library(ggplot2)
library(magrittr)
```

The basis for this example code is the [adrac data collection](https://github.com/dirkseidensticker/aDRAC) by Dirk Seidensticker. So let's download the data with the relevant c14bazAAR getter function: 

```{r, include = FALSE}
adrac <- get_c14data("adrac")
```

```{r, eval = FALSE}
adrac <- get_c14data("adrac")
```

## Temporal plotting of radiocarbon ages

Radiocarbon dating is a method to determine the absolute age of samples. Therefore one of the main aims for plotting naturally is to display temporal information. Let's select the dates of one individual site -- *Batalimo* in Central Africa -- as a subset to reproduce two of the most common types of radiocarbon date plots.

```{r}
Batalimo <- adrac %>%
  dplyr::filter(site == "Batalimo")
```

If age modelling and date plotting on a local or regional scale is the major aim of your analysis, you might want to take a look at the [oxcAAR](https://CRAN.R-project.org/package=oxcAAR) package. It serves as an R interface to [OxCal](https://c14.arch.ox.ac.uk/oxcal.html) and provides powerful default plotting methods

### Ridgeplots of density distributions

One way to plot radiocarbon ages is to show the probability density distribution of individual calibrated dates as ridgeplots. To produce a plot like that, we first of all need the age-probability information for each date. We can calculate that with the function `c14bazAAR::calibrate()`. 

```{r, include = FALSE}
Batalimo_calibrated <- Batalimo %>%
  calibrate(choices = "calprobdistr")
```

```{r, eval = FALSE}
Batalimo_calibrated <- Batalimo %>%
  calibrate(choices = "calprobdistr")
```

This adds a list column `calprobdistr` to the input `c14_date_list`. The list column contains a nested data.frame for each date with its probability distribution. 

```{r, echo=FALSE}
Batalimo_calibrated
```

With `tidyr::unnest()` the list column can be dissolved ("unnested") and integrated into the initial `c14_date_list`. Of course the latter looses its original structure and meaning with this step. Each row now represents the probability for one date and year.

```{r}
Batalimo_cal_dens <- Batalimo_calibrated %>% tidyr::unnest(cols = c("calprobdistr"))
```

A table like that can be used for plotting a ridgeplot. 

```{r, warning=FALSE}
Batalimo_cal_dens %>%
  ggplot() +
  # a special geom for ridgeplots is provided by the ggridges package
  ggridges::geom_ridgeline(
    # the relevant variables that have to be mapped for this geom are 
    # x (the time -- here the calibrated age transformed to calBC), 
    # y (the individual lab number of the dates) and
    # height (the probability for each year and date) 
    aes(x = -calage + 1950, y = labnr, height = density),
    # ridgeplots lack a scientifically clear y axis for each 
    # distribution plot and we can adjust the scaling to our needs
    scale = 300
  ) +
  xlab("age calBC/calAD") +
  ylab("dates")
```

### Calcurve plot

Another way to plot radiocarbon dates is to project them onto a calibration curve. The [Bchron](https://CRAN.R-project.org/package=Bchron) R package contains a data.frame with the [intcal13](https://www.doi.org/10.2458/azu_js_rc.55.16947) calibration curve data. We can load the `intcal13` table directly from Bchron.

```{r}
load(system.file('data/intcal13.rda', package = 'Bchron'))
```

For this kind of plot it is more convenient to work with the simplified `calrange` output of `c14bazAAR::calibrate()`.

```{r, include = FALSE}
Batalimo_calibrated <- Batalimo %>%
  calibrate(choices = "calrange")
```

```{r, eval = FALSE}
Batalimo_calibrated <- Batalimo %>%
  calibrate(choices = "calrange")
```

Like the `calprobdistr` option this also adds a list column to the input `c14_date_list`, but a much smaller one. For each date only the age ranges that make up the 2-sigma significance interval of the probability distribution are stored.

```{r}
Batalimo_calibrated$calrange[1:3]
```

The resulting table can also be unnested to make the list column content available in the main table.

```{r}
Batalimo_cal_range <- Batalimo_calibrated %>% tidyr::unnest(cols = c("calrange"))
```

Now we can plot the calibration curve and -- on top -- error bars with the `calrange` sequences.

```{r, warning=FALSE}
ggplot() +
  # line plot of the intcal curve
  geom_line(
    data = intcal13,
    # again we transform the age information from BP to BC
    mapping = aes(x = -V1 + 1950, y = -V2 + 1950)
  ) +
  # the errorbars are plotted on top of the curve
  geom_errorbarh(
    data = Batalimo_cal_range,
    mapping = aes(y = -c14age + 1950, xmin = -to + 1950, xmax = -from + 1950)
  ) +
  # we define the age range manually -- typically the calcurve
  # is arranged to go from the top left to the bottom right corner
  xlim(-1000, 2000) +
  ylim(2000, -1000) +
  xlab("age calBC/calAD") +
  ylab("uncalibrated age BC/AD")
```

## Spatial mapping of radiocarbon dates

Most radiocarbon dates that can be accessed with c14bazAAR have coordinate information for the respective sites where the samples were taken. Spatial maps therefore are an important form of data visualization as well.

`c14_date_list`s can directly be transformed to objects of class `sf`. `sf` objects were introduced by the R package [sf](https://CRAN.R-project.org/package=sf) which provides a tidy interface to work with spatial data in R.

```{r}
adrac_sf <- adrac %>% as.sf()
```

This tabular data structure contains the spatial point information for each date in a column *geom*, but also the initial columns of the input dataset: *data.\**

```{r, echo=FALSE}
adrac_sf %>% dplyr::select(data.labnr, data.c14age, data.c14std, geom)
```

It can be manipulated with the powerful dplyr functions. We `filter` out all dates from one particular publication (*Moga 2008*), `group` the dates `by` *site* and apply the `summarise` command to keep only one value per group. As we do not define an operation to fold the other variables in the input table, they are removed. Only the geometry column remains.

```{r}  
Moga_spatial <- adrac_sf %>%
  dplyr::filter(grepl("Moga 2008", data.shortref)) %>%
  dplyr::group_by(data.site) %>%
  dplyr::summarise(.groups = "drop")
```

### Interactive map view

The resulting `sf` object can be plotted interactively with the [mapview](https://CRAN.R-project.org/package=mapview) package. 

```{r}
# Moga_spatial %>% mapview::mapview()
```

### Static map plot

The `sf` object can also be used for a static plot -- which is useful for publications. We download some simple country border base map vector data with the [rnaturalearth](https://CRAN.R-project.org/package=rnaturalearth) R package and transform it to `sf` as well.

```{r}
countries <- rnaturalearth::ne_countries() %>% sf::st_as_sf()
```

Now we can combine the base layer and our point data to create the prototype of a static map plot. 

```{r, warning=FALSE}
ggplot() +
  # geom_sf is a special geom to handle spatial data in the sf format
  geom_sf(data = countries) +
  # the explicit mapping of variables is not necessary here, as geom_sf 
  # automatically finds the *geom* column in the input table
  geom_sf_text(data = countries, mapping = aes(label = formal_en), size = 2) +
  geom_sf(data = Moga_spatial) +
  # with geom_sf comes coord_sf to manage the underlying coordinate grid
  coord_sf(xlim = c(10, 30), ylim = c(0, 15))
```

Please feel free to open an issue [here](https://github.com/ropensci/c14bazAAR/issues) if you have questions about plotting radiocarbon dates.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_14cpalaeolithic.R, R/get_14sea.R,
%   R/get_adrac.R, R/get_agrichange.R, R/get_aida.R, R/get_austarch.R,
%   R/get_bda.R, R/get_c14data.R, R/get_calpal.R, R/get_caribbean.R,
%   R/get_context.R, R/get_eubar.R, R/get_euroevol.R, R/get_irdd.R,
%   R/get_jomon.R, R/get_katsianis.R, R/get_kiteeastafrica.R,
%   R/get_medafricarbon.R, R/get_mesorad.R, R/get_nerd.R, R/get_pacea.R,
%   R/get_palmisano.R, R/get_radon.R, R/get_radonb.R, R/get_rxpand.R,
%   R/get_sard.R
\name{get_14cpalaeolithic}
\alias{get_14cpalaeolithic}
\alias{get_14sea}
\alias{get_adrac}
\alias{get_agrichange}
\alias{get_aida}
\alias{get_austarch}
\alias{get_bda}
\alias{get_all_dates}
\alias{get_calpal}
\alias{get_caribbean}
\alias{get_context}
\alias{get_eubar}
\alias{get_euroevol}
\alias{get_irdd}
\alias{get_jomon}
\alias{get_katsianis}
\alias{get_kiteeastafrica}
\alias{get_medafricarbon}
\alias{get_mesorad}
\alias{get_nerd}
\alias{get_pacea}
\alias{get_palmisano}
\alias{get_radon}
\alias{get_radonb}
\alias{get_rxpand}
\alias{get_sard}
\title{Backend functions for data download}
\usage{
get_14cpalaeolithic(db_url = get_db_url("14cpalaeolithic"))

get_14sea(db_url = get_db_url("14sea"))

get_adrac(db_url = get_db_url("adrac"))

get_agrichange(db_url = get_db_url("agrichange"))

get_aida(db_url = get_db_url("aida"))

get_austarch(db_url = get_db_url("austarch"))

get_bda(db_url = get_db_url("bda"))

get_all_dates()

get_calpal(db_url = get_db_url("calpal"))

get_caribbean(db_url = get_db_url("caribbean"))

get_context(db_url = get_db_url("context"))

get_eubar(db_url = get_db_url("eubar"))

get_euroevol(db_url = get_db_url("euroevol"))

get_irdd(db_url = get_db_url("irdd"))

get_jomon(db_url = get_db_url("jomon"))

get_katsianis(db_url = get_db_url("katsianis"))

get_kiteeastafrica(db_url = get_db_url("kiteeastafrica"))

get_medafricarbon(db_url = get_db_url("medafricarbon"))

get_mesorad(db_url = get_db_url("mesorad"))

get_nerd(db_url = get_db_url("nerd"))

get_pacea(db_url = get_db_url("pacea"))

get_palmisano(db_url = get_db_url("palmisano"))

get_radon(db_url = get_db_url("radon"))

get_radonb(db_url = get_db_url("radonb"))

get_rxpand(db_url = get_db_url("rxpand"))

get_sard(db_url = get_db_url("sard"))
}
\arguments{
\item{db_url}{Character. URL that points to the c14 archive file. \code{c14bazAAR::get_db_url()}
fetches the URL from a reference list}
}
\description{
Backend functions to download data. See \code{?\link{get_c14data}}
for a more simple interface and further information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14_date_list_duplicates_remove.R
\name{duplicates}
\alias{duplicates}
\alias{remove_duplicates}
\alias{remove_duplicates.default}
\alias{remove_duplicates.c14_date_list}
\title{Remove duplicates in a \strong{c14_date_list}}
\usage{
remove_duplicates(
  x,
  preferences = NULL,
  supermerge = FALSE,
  log = TRUE,
  mark_only = FALSE
)

\method{remove_duplicates}{default}(
  x,
  preferences = NULL,
  supermerge = FALSE,
  log = TRUE,
  mark_only = FALSE
)

\method{remove_duplicates}{c14_date_list}(
  x,
  preferences = NULL,
  supermerge = FALSE,
  log = TRUE,
  mark_only = FALSE
)
}
\arguments{
\item{x}{an object of class c14_date_list}

\item{preferences}{character vector with the order of source databases by
which the deduping should be executed. If e.g. preferences = c("radon", "calpal")
and a certain date appears in radon and euroevol, then only the radon entry remains.
Default: NULL. With preferences = NULL all overlapping, conflicting information in
individual columns of one duplicated date is removed. See Option 2 and 3.}

\item{supermerge}{boolean. Should the duplicated datasets be merged on the column level?
Default: FALSE. See Option 3.}

\item{log}{logical. If log = TRUE, an additional column is added that contains a string
documentation of all variants of the information for one date from all conflicting
databases. Default = TRUE.}

\item{mark_only}{boolean. Should duplicates not be removed, but only indicated? Default: FALSE.}
}
\value{
an object of class c14_date_list with the additional
columns \strong{duplicate_group} or \strong{duplicate_remove_log}
}
\description{
Duplicates are found by comparison of \strong{labnr}s.
Only dates with exactly equal \strong{labnr}s are considered duplicates.
Duplicate groups are numbered (from 0) and these numbers linked to
the individual dates in a internal column \strong{duplicate_group}.
If you only want to see this grouping without removing anything use the \code{mark_only} flag.
\code{c14bazAAR::remove_duplicates()} can remove duplicates with three different strategies
according to the value of the arguments \code{preferences} and \code{supermerge}:
\enumerate{
  \item Option 1: By merging all dates in a \strong{duplicate_group}. All non-equal variables
  in the duplicate group are turned to \code{NA}. This is the default option.
  \item Option 2: By selecting individual database entries in a \strong{duplicate_group}
  according to a trust hierarchy as defined by the parameter \code{preferences}.
  In case of duplicates within one database the first occurrence in the table (top down)
  is selected. All databases not mentioned in \code{preferences} are dropped.
  \item Option 3: Like option 2, but in this case the different datasets in a
  \strong{duplicate_group} are merged column by column to
  create a superdataset with a maximum of information. The column \strong{sourcedb} is
  dropped in this case to indicate that multiple databases have been merged. Data
  citation is a lot more difficult with this option. It can be activated with \code{supermerge}.
}
The option \code{log} allows to add a new column \strong{duplicate_remove_log}
that documents the variety of values provided by all databases for this
duplicated date.
}
\examples{
library(magrittr)

test_data <- tibble::tribble(
  ~sourcedb, ~labnr,  ~c14age, ~c14std,
 "A",       "lab-1", 1100,    10,
 "A",       "lab-1", 2100,    20,
 "B",       "lab-1", 3100,    30,
 "A",       "lab-2", NA,      10,
 "B",       "lab-2", 2200,    20,
 "C",       "lab-3", 1300,    10
) \%>\% as.c14_date_list()

# remove duplicates with option 1:
test_data \%>\% remove_duplicates()

# remove duplicates with option 2:
test_data \%>\% remove_duplicates(
  preferences = c("A", "B")
)

# remove duplicates with option 3:
test_data \%>\% remove_duplicates(
  preferences = c("A", "B"),
  supermerge = TRUE
)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14bazAAR.R
\docType{package}
\name{c14bazAAR-package}
\alias{c14bazAAR}
\alias{c14bazAAR-package}
\title{c14bazAAR: Download and Prepare C14 Dates from Different Source Databases}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Query different C14 date databases and apply basic data cleaning, merging and calibration steps. Currently available databases: 14cpalaeolithic, 14sea, adrac, agrichange, aida, austarch, bda, calpal, caribbean, context, eubar, euroevol, irdd, jomon, katsianis, kiteeastafrica, medafricarbon, mesorad, nerd, pacea, palmisano, radon, radonb, rxpand, sard.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/c14bazAAR}
  \item \url{https://github.com/ropensci/c14bazAAR}
  \item Report bugs at \url{https://github.com/ropensci/c14bazAAR/issues}
}

}
\author{
\strong{Maintainer}: Clemens Schmid \email{clemens@nevrome.de} (\href{https://orcid.org/0000-0003-3448-5715}{ORCID}) [copyright holder]

Authors:
\itemize{
  \item Dirk Seidensticker \email{dirk.seidensticker@gmail.com} (\href{https://orcid.org/0000-0002-8155-7702}{ORCID})
  \item Daniel Knitter \email{knitter@geographie.uni-kiel.de} (\href{https://orcid.org/0000-0003-3014-4497}{ORCID})
  \item Martin Hinz \email{martin.hinz@iaw.unibe.ch} (\href{https://orcid.org/0000-0002-9904-6548}{ORCID})
  \item David Matzig \email{matzigdavid@gmail.com} (\href{https://orcid.org/0000-0001-7349-5401}{ORCID})
  \item Wolfgang Hamer \email{hamer@geographie.uni-kiel.de} (\href{https://orcid.org/0000-0002-5943-5020}{ORCID})
  \item Kay Schmütz \email{kschmuetz@ufg.uni-kiel.de}
}

Other contributors:
\itemize{
  \item Thomas Huet \email{thomashuet7@gmail.com} (\href{https://orcid.org/0000-0002-1112-6122}{ORCID}) [contributor]
  \item Nils Mueller-Scheessel \email{nils.mueller-scheessel@ufg.uni-kiel.de} (\href{https://orcid.org/0000-0001-7992-8722}{ORCID}) [contributor]
  \item Ben Marwick (\href{https://orcid.org/0000-0001-7879-4531}{ORCID}) [reviewer]
  \item Enrico R. Crema (\href{https://orcid.org/0000-0001-6727-5138}{ORCID}) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{deprecated_functions}
\alias{deprecated_functions}
\alias{mark_duplicates}
\alias{coordinate_precision}
\alias{finalize_country_name}
\alias{standardize_country_name}
\alias{get_emedyd}
\alias{fix_database_country_name}
\alias{classify_material}
\title{Deprecated functions}
\usage{
mark_duplicates(...)

coordinate_precision(...)

finalize_country_name(...)

standardize_country_name(...)

get_emedyd(...)

fix_database_country_name(...)

classify_material(...)
}
\arguments{
\item{...}{...}
}
\description{
Run them anyway to get some information about their replacements
or why they were removed.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14_date_list_calibrate.R
\name{calibrate}
\alias{calibrate}
\alias{calibrate.default}
\alias{calibrate.c14_date_list}
\title{Calibrate all valid dates in a \strong{c14_date_list}}
\usage{
calibrate(
  x,
  choices = c("calrange"),
  sigma = 2,
  calCurves = rep("intcal20", nrow(x)),
  ...
)

\method{calibrate}{default}(
  x,
  choices = c("calrange"),
  sigma = 2,
  calCurves = rep("intcal20", nrow(x)),
  ...
)

\method{calibrate}{c14_date_list}(
  x,
  choices = c("calrange"),
  sigma = 2,
  calCurves = rep("intcal20", nrow(x)),
  ...
)
}
\arguments{
\item{x}{an object of class c14_date_list}

\item{choices}{whether the result should include the full calibrated
probability dataframe ('calprobdistr') or the sigma range ('calrange').
Both arguments may be given at the same time.}

\item{sigma}{the desired sigma value (1,2,3) for the calibrated sigma ranges}

\item{calCurves}{a vector of values containing either intcal20, shcal20,
marine20, or normal (older calibration curves are supposed such as intcal13).
Should be the same length the number of ages supplied.
See \link[Bchron]{BchronCalibrate} for more information}

\item{...}{passed to \link[Bchron]{BchronCalibrate}}
}
\value{
an object of class c14_date_list with the additional columns
\strong{calprobdistr} or \strong{calrange} and \strong{sigma}
}
\description{
Calibrate all dates in a \strong{c14_date_list} with
\code{Bchron::BchronCalibrate()}. The function provides two different
kinds of output variables that are added as new list columns to the input
\strong{c14_date_list}: \strong{calprobdistr} and \strong{calrange}.
\strong{calrange} is accompanied by \strong{sigma}. See
\code{?Bchron::BchronCalibrate} and \code{?c14bazAAR:::hdr} for some more
information. \cr
\strong{calprobdistr}: The probability distribution of the individual date
for all ages with an individual probability >= 1e-06. For each date there's
a data.frame with the columns \strong{calage} and \strong{density}. \cr
\strong{calrange}: The contiguous ranges which cover the probability interval
requested for the individual date. For each date there's a data.frame with the
columns \strong{dens} and \strong{from} and \strong{to}.
}
\examples{
calibrate(
  example_c14_date_list,
  choices = c("calprobdistr", "calrange"),
  sigma = 1
)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_c14data.R
\name{get_c14data}
\alias{get_c14data}
\title{Download radiocarbon source databases and convert them to a \strong{c14_date_list}}
\usage{
get_c14data(databases = c())
}
\arguments{
\item{databases}{Character vector. Names of databases to be downloaded. "all" causes the download of all databases. \code{get_c14data()} prints a list of the currently available databases}
}
\description{
\code{get_c14data()} allows to download source databases and adjust their variables to conform to the
definition in the
\href{https://github.com/ropensci/c14bazAAR/blob/master/data-raw/variable_reference.csv}{variable_reference}
table. That includes renaming and arranging the variables (with \code{c14bazAAR::order_variables()})
as well as type conversion (with \code{c14bazAAR::enforce_types()}) -- so all the steps undertaken by
\code{as.c14_date_list()}. \cr
All databases require different downloading and data wrangling steps. Therefore
there's a custom getter function for each of them (see \code{?get_all_dates}). \cr

\code{get_c14data()} is a wrapper to download all dates from multiple databases and
\code{c14bazAAR::fuse()} the results.
}
\examples{

\dontrun{
 get_c14data(databases = c("adrac", "palmisano"))
  get_all_dates()}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14_date_list_country_attribution.R
\name{country_attribution}
\alias{country_attribution}
\alias{determine_country_by_coordinate}
\alias{determine_country_by_coordinate.default}
\alias{determine_country_by_coordinate.c14_date_list}
\title{Determine the country of all dates in a \strong{c14_date_list} from their coordinates}
\usage{
determine_country_by_coordinate(x, suppress_spatial_warnings = TRUE)

\method{determine_country_by_coordinate}{default}(x, suppress_spatial_warnings = TRUE)

\method{determine_country_by_coordinate}{c14_date_list}(x, suppress_spatial_warnings = TRUE)
}
\arguments{
\item{x}{an object of class c14_date_list}

\item{suppress_spatial_warnings}{suppress some spatial data messages and warnings}
}
\value{
an object of class c14_date_list with the additional column \strong{country_coord}
}
\description{
\code{c14bazAAR::determine_country_by_coordinate()} adds the column
\strong{country_coord} with standardized country attribution based on the coordinate
information for the dates.
Due to the inconsistencies in the \strong{country} column in many c14 source databases
it's often necessary to rely on the coordinate position (\strong{lat} & \strong{lon})
for country attribution information. Unfortunately not all source databases store
coordinates.
}
\examples{
library(magrittr)
example_c14_date_list \%>\%
  determine_country_by_coordinate()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{db_info_table}
\alias{db_info_table}
\title{Database lookup table}
\format{
a data.frame. Columns:
\itemize{
 \item{db: database name}
 \item{version: database version}
 \item{url_num: url number (some databases are spread over multiple files)}
 \item{url: file url where the database can be downloaded}
}
}
\description{
Lookup table for general source database information.
}
\concept{lookup_tables}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_dbs.R
\name{get_db_url}
\alias{get_db_url}
\alias{get_db_version}
\title{Get information for c14 databases}
\usage{
get_db_url(..., db_info_table = c14bazAAR::db_info_table)

get_db_version(..., db_info_table = c14bazAAR::db_info_table)
}
\arguments{
\item{...}{names of the databases}

\item{db_info_table}{db info reference table}
}
\description{
Looks for information for the c14 source databases in \link{db_info_table}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14_date_list_fuse.R
\name{fuse}
\alias{fuse}
\alias{fuse.default}
\alias{fuse.c14_date_list}
\title{Fuse multiple \strong{c14_date_list}s}
\usage{
fuse(...)

\method{fuse}{default}(...)

\method{fuse}{c14_date_list}(...)
}
\arguments{
\item{...}{objects of class c14_date_list}
}
\value{
an object of class c14_date_list
}
\description{
This function combines \strong{c14_date_list}s with
\code{dplyr::bind_rows()}. \cr
This is not a joining operation and it therefore
might introduce duplicates. See \code{c14bazAAR::mark_duplicates()}
and \code{c14bazAAR::remove_duplicates()} for a way to find and remove
them.
}
\examples{
# fuse three identical example c14_date_lists
fuse(example_c14_date_list, example_c14_date_list, example_c14_date_list)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14_date_list_convert.R
\name{as.sf}
\alias{as.sf}
\alias{as.sf.default}
\alias{as.sf.c14_date_list}
\title{Convert a \strong{c14_date_list} to a sf object}
\usage{
as.sf(x, quiet = FALSE)

\method{as.sf}{default}(x, quiet = FALSE)

\method{as.sf}{c14_date_list}(x, quiet = FALSE)
}
\arguments{
\item{x}{an object of class c14_date_list}

\item{quiet}{suppress warning about the removal of dates without coordinates}
}
\value{
an object of class sf
}
\description{
Most 14C dates have point position information in
the coordinates columns \strong{lat} and \strong{lon}. This allows
them to be converted to a spatial simple feature collection as provided
by the \code{sf} package. This simplifies for example mapping of the
dates.
}
\examples{
sf_c14 <- as.sf(example_c14_date_list)

\dontrun{
library(mapview)
mapview(sf_c14$geom)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14_date_list_basic.R
\name{c14_date_list}
\alias{c14_date_list}
\alias{as.c14_date_list}
\alias{is.c14_date_list}
\alias{format.c14_date_list}
\alias{print.c14_date_list}
\alias{plot.c14_date_list}
\title{\strong{c14_date_list}}
\usage{
as.c14_date_list(x, ...)

is.c14_date_list(x, ...)

\method{format}{c14_date_list}(x, ...)

\method{print}{c14_date_list}(x, ...)

\method{plot}{c14_date_list}(x, ...)
}
\arguments{
\item{x}{an object}

\item{...}{further arguments passed to or from other methods}
}
\description{
The \strong{c14_date_list} is the central data structure of the
\code{c14bazAAR} package. It's a tibble with set of custom methods and
variables. Please see the
\href{https://github.com/ropensci/c14bazAAR/blob/master/data-raw/variable_reference.csv}{variable_reference}
table for a description of the variables. Further available variables are ignored. \cr
If an object is of class data.frame or tibble (tbl & tbl_df), it can be
converted to an object of class \strong{c14_date_list}. The only requirement
is that it contains the essential columns \strong{c14age} and \strong{c14std}.
The \code{as} function adds the string "c14_date_list" to the classes vector
of the object and applies \code{order_variables()}, \code{enforce_types()} and
the helper function \code{clean_latlon()} to it.
}
\examples{
as.c14_date_list(data.frame(c14age = c(2000, 2500), c14std = c(30, 35)))
is.c14_date_list(5) # FALSE
is.c14_date_list(example_c14_date_list) # TRUE

print(example_c14_date_list)
plot(example_c14_date_list)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14_date_list_write_c14.R
\name{write_c14}
\alias{write_c14}
\alias{write_c14.default}
\alias{write_c14.c14_date_list}
\title{write \strong{c14_date_list}s to files}
\usage{
write_c14(x, format = c("csv"), ...)

\method{write_c14}{default}(x, format = c("csv"), ...)

\method{write_c14}{c14_date_list}(x, format = c("csv"), ...)
}
\arguments{
\item{x}{an object of class c14_date_list}

\item{format}{the output format: 'csv' (default) or 'xlsx'.
'csv' calls \code{utils::write.csv()},
'xlsx' calls \code{writexl::write_xlsx()}}

\item{...}{passed to the actual writing functions}
}
\description{
write \strong{c14_date_list}s to files
}
\examples{
csv_file <- tempfile(fileext = ".csv")
write_c14(
  example_c14_date_list,
  format = "csv",
  file = csv_file
)
\donttest{
xlsx_file <- tempfile(fileext = ".xlsx")
write_c14(
  example_c14_date_list,
  format = "xlsx",
  path = xlsx_file,
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14_date_list_order_variables.R
\name{order_variables}
\alias{order_variables}
\alias{order_variables.default}
\alias{order_variables.c14_date_list}
\title{Order the variables in a \strong{c14_date_list}}
\usage{
order_variables(x)

\method{order_variables}{default}(x)

\method{order_variables}{c14_date_list}(x)
}
\arguments{
\item{x}{an object of class c14_date_list}
}
\value{
an object of class c14_date_list
}
\description{
Arrange variables according to a defined order. This makes
sure that a \strong{c14_date_list} always appears with the same
outline. \cr
A \strong{c14_date_list} has at least the columns \strong{c14age}
and \strong{c14std}. Beyond that there's a selection of additional
variables depending on the input from the source databases, as a
result of the \code{c14bazAAR} functions or added by other data
analysis steps. This function arranges the expected variables in
a distinct, predefined order. Undefined variables are added at the
end.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c14_date_list_enforce_types.R
\name{enforce_types}
\alias{enforce_types}
\alias{enforce_types.default}
\alias{enforce_types.c14_date_list}
\title{Enforce variable types in a \strong{c14_date_list}}
\usage{
enforce_types(x, suppress_na_introduced_warnings = TRUE)

\method{enforce_types}{default}(x, suppress_na_introduced_warnings = TRUE)

\method{enforce_types}{c14_date_list}(x, suppress_na_introduced_warnings = TRUE)
}
\arguments{
\item{x}{an object of class c14_date_list}

\item{suppress_na_introduced_warnings}{suppress warnings caused by data removal in
type transformation due to wrong database entries (such as text in a number column)}
}
\value{
an object of class c14_date_list
}
\description{
Enforce variable types in a \strong{c14_date_list} and remove
everything that doesn't fit (e.g. text in a number field).
See the
\href{https://github.com/ropensci/c14bazAAR/blob/master/data-raw/variable_reference.csv}{variable_reference}
table for a documentation of the variable types.
\code{enforce_types()} is called in \code{c14bazAAR::as.c14_date_list()}.
}
\examples{
# initial situation
ex <- example_c14_date_list
class(ex$c14age)

# modify variable/column type
ex$c14age <- as.character(ex$c14age)
class(ex$c14age)

# fix type with enforce_types()
ex <- enforce_types(ex)
class(ex$c14age)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{example_c14_date_list}
\alias{example_c14_date_list}
\title{Example c14_date_list}
\format{
a c14_date_list.
See data_raw/variable_definition.csv for an explanation of
the variable meaning.
}
\description{
c14_date_list for tests and example code.
}
\concept{c14_date_lists}
