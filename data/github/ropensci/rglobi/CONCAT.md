
# rglobi 
R library to access species interaction data of http://globalbioticinteractions.org

[![R-check](https://github.com/ropensci/rglobi/workflows/R-check/badge.svg)](https://github.com/ropensci/rglobi/actions) [![Build Status](https://travis-ci.org/ropensci/rglobi.svg?branch=master)](https://travis-ci.org/ropensci/rglobi) [![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/rglobi?color=E664A4)](https://github.com/r-hub/cranlogs.app) [![cran version](https://www.r-pkg.org/badges/version/rglobi)](https://CRAN.R-project.org/package=rglobi) [![cran checks](https://cranchecks.info/badges/worst/rglobi)](https://cranchecks.info/pkgs/rglobi)


## install 
To install ```rglobi``` from [CRAN](https://CRAN.R-project.org/package=rglobi):
```R
install.packages("rglobi")
```

Or install development version:
```R
install.packages("devtools")
devtools::install_github("ropensci/rglobi")
```

## examples

```R
library(rglobi)
# find all unique prey names of Homo sapiens
prey_of("Homo sapiens")$target_taxon_name
# is a shortcut of
get_interactions_by_taxa(sourcetaxon='Homo sapiens', interactiontype='preysOn')$target_taxon_name

# list of supported interactions types
get_interaction_types()

# all known prey names and locations (latitude, longitude) where birds (Aves) preyed on rodents (Rodentia) in California
obs <- get_interactions_by_taxa(sourcetaxon = "Aves", bbox=c(-125.53344800000002,32.750323,-114.74487299999998,41.574361), targettaxon = "Rodentia", returnobservations=T)
locations <- cbind(obs$target_taxon_name, obs$latitude, obs$longitude)
```
Please see R help pages (e.g. ```?get_interactions_by_taxa``` and [vignettes](https://CRAN.R-project.org/package=rglobi) for more information.

## tests
Tests can be executed using devtools package.
```R
# workdir should be rglobi repo root directory (check with getwd())
# install dependencies 
devtools::install('.')
devtools::test()
```
This should reload the library, executes the test_that testcases and show test reports.

## documentation
roxygen2 is used to generate .Rd and NAMESPACE by running:
```R
 library(roxygen2)
 roxygenize(".")
```

Vignettes are generated using ```knitr``` and ```markdown``` packages.

## meta

Please [report any issues or bugs](https://github.com/ropensci/rglobi/issues).

This package is part of the [rOpenSci](http://ropensci.org/packages) project.

[![rOpenSci footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
Dear Reviewers:

I built this package using [R CMD build .] and checked it with command [R CMD check --as-cran rglobi_0.2.27.tar.gz] on my Ubuntu 18.04 LTS Bionic Beaver (4.15.0-159-generic) machine running R 4.1.1 (2021-08-10). Also, the checks were executed on travis-ci.com (running Ubuntu), github workflows (windows/mac/ubuntu), and winbuilder (running windows). 

IMPROVEMENT
* improve test cases to skip with message if web apis are unavailable

I've checked the mis-spelled words and confirmed that they are in fact not mis-spelled. 

Thank you for taking the time to review my submission.

-jorrit
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Introduction to the rglobi package}
-->



rglobi vignette - an R interface to the aggregated biotic interaction data of GloBI (http://globalbioticinteractions.org)
======

### About the package

`rglobi` allows easy access to the GLoBI database by wrapping the existing API.  

********************

### Quick start

#### First, install `rglobi`



```r
install.packages("devtools")
```

```
## 
## The downloaded binary packages are in
## 	/var/folders/v4/0_y78_xj1x7b72l89cwc2kq80000gp/T//RtmpF6qc1T/downloaded_packages
```

```r
devtools::install_github("ropensci/rglobi")
```

```
## Downloading github repo ropensci/rglobi@master
## Installing rglobi
## '/Library/Frameworks/R.framework/Resources/bin/R' --vanilla CMD INSTALL  \
##   '/private/var/folders/v4/0_y78_xj1x7b72l89cwc2kq80000gp/T/RtmpF6qc1T/devtools15b730380a08/ropensci-rglobi-0b49250'  \
##   --library='/Library/Frameworks/R.framework/Versions/3.1/Resources/library'  \
##   --install-tests 
## 
## Reloading installed rglobi
```


```r
library("rglobi")
```

Note that since rglobi is still pretty new, only the dev version of the library is available at this time. We hope to publish a version to CRAN once the library matures. For more information see [GitHub](https://github.com/ropensci/rglobi).

#### Find interactions involving a certain species

Determining which species interact with each other (and how and where) is a major focus of ecology.  The Global Biotic Interactions Database offers data on multiple interaction types, and the `rglobi` library offers several ways to specify and access these interactions. `get_interactions()` is the primary function used if data on a specific interaction type is required.  A focus taxon is entered (may be specified as "Genus species" or higher level (e.g., Genus, Family, Class)).


```r
hsapiens <- get_interactions(taxon = "Homo sapiens", interaction.type = "preysOn")
head(hsapiens)
```

```
##   source_taxon_external_id source_taxon_name
## 1               EOL:327955      Homo sapiens
## 2               EOL:327955      Homo sapiens
## 3               EOL:327955      Homo sapiens
## 4               EOL:327955      Homo sapiens
## 5               EOL:327955      Homo sapiens
## 6               EOL:327955      Homo sapiens
##                                                                                                                                                                                                           source_taxon_path
## 1 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 2 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 3 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 4 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 5 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 6 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
##   source_specimen_life_stage interaction_type target_taxon_external_id
## 1                         NA          preysOn               EOL:328048
## 2                         NA          preysOn               EOL:328716
## 3                         NA          preysOn                 EOL:1666
## 4                         NA          preysOn               EOL:326518
## 5                         NA          preysOn              EOL:4453589
## 6                         NA          preysOn              EOL:2815988
##                   target_taxon_name
## 1                 Nandinia binotata
## 2             Cephalophus monticola
## 3                           Manidae
## 4               Atherurus africanus
## 5 Cercopithecus nictitans nictitans
## 6                         Serpentes
##                                                                                                                                                                                                                                                                                                      target_taxon_path
## 1                                                                                                                                                                                                                              Animalia | Chordata | Mammalia | Carnivora | Nandiniidae | Nandinia | Nandinia binotata
## 2                                                                                                                                                                                                                        Animalia | Chordata | Mammalia | Artiodactyla | Bovidae | Cephalophus | Cephalophus monticola
## 3                                                                                                                                                                                                                                                                 Animalia | Chordata | Mammalia | Pholidota | Manidae
## 4                                                                                                                                                                                                                            Animalia | Chordata | Mammalia | Rodentia | Hystricidae | Atherurus | Atherurus africanus
## 5 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Cercopithecoidea | Cercopithecidae | Cercopithecinae | Cercopithecini | Cercopithecus | Cercopithecus nictitans | Cercopithecus nictitans nictitans
## 6                                                                                                                                                                                           Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Reptilia | Squamata | Serpentes
##   target_specimen_life_stage latitude longitude study_title
## 1                         NA       NA        NA          NA
## 2                         NA       NA        NA          NA
## 3                         NA       NA        NA          NA
## 4                         NA       NA        NA          NA
## 5                         NA       NA        NA          NA
## 6                         NA       NA        NA          NA
```

`get_predators_of()` and `get_prey_of()` are wrappers for `get_interactions()` that remove the need to specify interaction type


```r
hsapiens <- get_prey_of("Homo sapiens")
head(hsapiens)
```

```
##   source_taxon_external_id source_taxon_name
## 1               EOL:327955      Homo sapiens
## 2               EOL:327955      Homo sapiens
## 3               EOL:327955      Homo sapiens
## 4               EOL:327955      Homo sapiens
## 5               EOL:327955      Homo sapiens
## 6               EOL:327955      Homo sapiens
##                                                                                                                                                                                                           source_taxon_path
## 1 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 2 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 3 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 4 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 5 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
## 6 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Hominoidea | Hominidae | Homininae | Homo | Homo sapiens
##   source_specimen_life_stage interaction_type target_taxon_external_id
## 1                         NA          preysOn               EOL:328048
## 2                         NA          preysOn               EOL:328716
## 3                         NA          preysOn                 EOL:1666
## 4                         NA          preysOn               EOL:326518
## 5                         NA          preysOn              EOL:4453589
## 6                         NA          preysOn              EOL:2815988
##                   target_taxon_name
## 1                 Nandinia binotata
## 2             Cephalophus monticola
## 3                           Manidae
## 4               Atherurus africanus
## 5 Cercopithecus nictitans nictitans
## 6                         Serpentes
##                                                                                                                                                                                                                                                                                                      target_taxon_path
## 1                                                                                                                                                                                                                              Animalia | Chordata | Mammalia | Carnivora | Nandiniidae | Nandinia | Nandinia binotata
## 2                                                                                                                                                                                                                        Animalia | Chordata | Mammalia | Artiodactyla | Bovidae | Cephalophus | Cephalophus monticola
## 3                                                                                                                                                                                                                                                                 Animalia | Chordata | Mammalia | Pholidota | Manidae
## 4                                                                                                                                                                                                                            Animalia | Chordata | Mammalia | Rodentia | Hystricidae | Atherurus | Atherurus africanus
## 5 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Primates | Haplorrhini | Simiiformes | Cercopithecoidea | Cercopithecidae | Cercopithecinae | Cercopithecini | Cercopithecus | Cercopithecus nictitans | Cercopithecus nictitans nictitans
## 6                                                                                                                                                                                           Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Reptilia | Squamata | Serpentes
##   target_specimen_life_stage latitude longitude study_title
## 1                         NA       NA        NA          NA
## 2                         NA       NA        NA          NA
## 3                         NA       NA        NA          NA
## 4                         NA       NA        NA          NA
## 5                         NA       NA        NA          NA
## 6                         NA       NA        NA          NA
```

For a complete list of supported interaction types, 


```r
get_interaction_types()
```

```
##      interaction     source     target
## 1        preysOn   predator       prey
## 2   preyedUponBy       prey   predator
## 3     parasiteOf   parasite       host
## 4    hasParasite       host   parasite
## 5     pollinates pollinator      plant
## 6   pollinatedBy      plant pollinator
## 7     pathogenOf   pathogen       host
## 8    hasPathogen      plant pollinator
## 9     symbiontOf     source     target
## 10 interactsWith     source     target
```

For data sources in which type of interactions was not specified, the interaction is labeled "interacts_with".

If you wish to view all interactions instead of specific types (e.g., parasitism and predation instead of just one of the two), the `get_interactions_by_taxa()` function allows this. In addition, the function provides options to search for interactions between two groups (source and target taxa, see above table) and only find interactions in a certain region.


```r
rattus <- get_interactions_by_taxa(sourcetaxon = "Rattus")
head(rattus)
```

```
##   source_taxon_external_id source_taxon_name
## 1               EOL:328447     Rattus rattus
## 2               EOL:328447     Rattus rattus
## 3               EOL:328447     Rattus rattus
## 4               EOL:328447     Rattus rattus
## 5               EOL:328447     Rattus rattus
## 6               EOL:328447     Rattus rattus
##                                                                                                                                                                             source_taxon_path
## 1 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Rodentia | Myomorpha | Muridae | Murinae | Rattus | Rattus rattus
## 2 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Rodentia | Myomorpha | Muridae | Murinae | Rattus | Rattus rattus
## 3 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Rodentia | Myomorpha | Muridae | Murinae | Rattus | Rattus rattus
## 4 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Rodentia | Myomorpha | Muridae | Murinae | Rattus | Rattus rattus
## 5 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Rodentia | Myomorpha | Muridae | Murinae | Rattus | Rattus rattus
## 6 Animalia | Bilateria | Deuterostomia | Chordata | Vertebrata | Gnathostomata | Tetrapoda | Mammalia | Theria | Eutheria | Rodentia | Myomorpha | Muridae | Murinae | Rattus | Rattus rattus
##   source_specimen_life_stage interaction_type target_taxon_external_id
## 1                         NA    interactsWith             EOL:13340393
## 2                         NA    interactsWith             EOL:10702710
## 3                         NA    interactsWith              EOL:5004474
## 4                         NA    interactsWith              EOL:4968671
## 5                         NA    interactsWith              EOL:4965412
## 6                         NA    interactsWith              EOL:4963889
##              target_taxon_name
## 1            Mastophorus muris
## 2           Trypanosoma lewisi
## 3           Calodium hepaticum
## 4                  Raillietina
## 5 Nippostrongylus brasiliensis
## 6                     Syphacia
##                                                                                                                  target_taxon_path
## 1                                     Animalia | Nematoda | Secernentea | Spirurida | Spiruridae | Mastophorus | Mastophorus muris
## 2 Cellular organisms | Eukaryota | Euglenozoa | Kinetoplastida | Trypanosomatidae | Trypanosoma | Herpetosoma | Trypanosoma lewisi
## 3                                                                                                               Calodium hepaticum
## 4                                                Animalia | Platyhelminthes | Cestoda | Cyclophyllidea | Davaineidae | Raillietina
## 5                 Animalia | Nematoda | Secernentea | Strongylida | Heligosomidae | Nippostrongylus | Nippostrongylus brasiliensis
## 6                                                            Animalia | Nematoda | Secernentea | Ascaridida | Oxyuridae | Syphacia
##   target_specimen_life_stage latitude longitude study_title
## 1                         NA       NA        NA          NA
## 2                         NA       NA        NA          NA
## 3                         NA       NA        NA          NA
## 4                         NA       NA        NA          NA
## 5                         NA       NA        NA          NA
## 6                         NA       NA        NA          NA
```

Only a source taxa need be identified, but you can also specify a target taxon and/or geographic boundary (Coordinates must be in decimal degrees ([EPSG:4326](http://spatialreference.org/ref/epsg/wgs-84/)) and correspond to the west, south, east, and northern boundaries (i.e., min x, min y, max x, max y).  


```r
aves_crustacea_northern_hemisphere <- get_interactions_by_taxa( sourcetaxon = "Aves", targettaxon = "Crustacea", bbox=c(-180, 0, 180, 90 ))
head(aves_crustacea_northern_hemisphere)
```

```
##   source_taxon_external_id       source_taxon_name
## 1                EOL:18884                 Cepphus
## 2              EOL:1047918      Anas platyrhynchos
## 3              EOL:1049560 Recurvirostra americana
## 4              EOL:1049376  Arenaria melanocephala
## 5              EOL:1048977       Bucephala albeola
## 6              EOL:1048977       Bucephala albeola
##                                                                                           source_taxon_path
## 1                                          Animalia | Chordata | Aves | Charadriiformes | Alcidae | Cepphus
## 2                          Animalia | Chordata | Aves | Anseriformes | Anatidae | Anas | Anas platyrhynchos
## 3 Animalia | Chordata | Aves | Charadriiformes | Recurvirostridae | Recurvirostra | Recurvirostra americana
## 4           Animalia | Chordata | Aves | Charadriiformes | Scolopacidae | Arenaria | Arenaria melanocephala
## 5                      Animalia | Chordata | Aves | Anseriformes | Anatidae | Bucephala | Bucephala albeola
## 6                      Animalia | Chordata | Aves | Anseriformes | Anatidae | Bucephala | Bucephala albeola
##   source_specimen_life_stage interaction_type target_taxon_external_id
## 1                         NA          preysOn              EOL:2625033
## 2                         NA          preysOn              EOL:2620777
## 3                         NA    interactsWith              EOL:2625033
## 4                         NA          preysOn              EOL:2620777
## 5                         NA    interactsWith              EOL:2625033
## 6                         NA          preysOn              EOL:2620777
##   target_taxon_name
## 1          Copepoda
## 2        Gammaridea
## 3          Copepoda
## 4        Gammaridea
## 5          Copepoda
## 6        Gammaridea
##                                                                                                                               target_taxon_path
## 1                                              Animalia | Bilateria | Protostomia | Ecdysozoa | Arthropoda | Crustacea | Maxillopoda | Copepoda
## 2 Animalia | Bilateria | Protostomia | Ecdysozoa | Arthropoda | Crustacea | Malacostraca | Eumalacostraca | Peracarida | Amphipoda | Gammaridea
## 3                                              Animalia | Bilateria | Protostomia | Ecdysozoa | Arthropoda | Crustacea | Maxillopoda | Copepoda
## 4 Animalia | Bilateria | Protostomia | Ecdysozoa | Arthropoda | Crustacea | Malacostraca | Eumalacostraca | Peracarida | Amphipoda | Gammaridea
## 5                                              Animalia | Bilateria | Protostomia | Ecdysozoa | Arthropoda | Crustacea | Maxillopoda | Copepoda
## 6 Animalia | Bilateria | Protostomia | Ecdysozoa | Arthropoda | Crustacea | Malacostraca | Eumalacostraca | Peracarida | Amphipoda | Gammaridea
##   target_specimen_life_stage latitude longitude study_title
## 1                         NA       NA        NA          NA
## 2                         NA       NA        NA          NA
## 3                         NA       NA        NA          NA
## 4                         NA       NA        NA          NA
## 5                         NA       NA        NA          NA
## 6                         NA       NA        NA          NA
```

#### Find interactions in a geographic areas

Instead of a taxon-specific approach, users may also wish to gather information on all interactions occuring in a specific area.  For example, a group developing ecoystem models for the Gulf of Mexico may want to consider all the data from that region.  `rglobi` enables this type of search by allowing users to specify a rectangular bounding box. Coordinates must be in decimal degrees ([EPSG:4326](http://spatialreference.org/ref/epsg/wgs-84/)) and correspond to the west, south, east, and northern boundaries (i.e., min x, min y, max x, max y).  


```r
gulfinteractions <- get_interactions_in_area( bbox=c(-97.0, 17.5, -81, 31))
head(gulfinteractions)
```

```
##   source_taxon_external_id        source_taxon_name
## 1                    EOL:1                 Animalia
## 2                 EOL:1905           Actinopterygii
## 3              EOL:2598871                Crustacea
## 4               EOL:208553 Chloroscombrus chrysurus
## 5               EOL:208553 Chloroscombrus chrysurus
## 6               EOL:208553 Chloroscombrus chrysurus
##                                                                                             source_taxon_path
## 1                                                                                                    Animalia
## 2                          Biota | Animalia | Chordata | Vertebrata | Gnathostomata | Pisces | Actinopterygii
## 3                                     Animalia | Bilateria | Protostomia | Ecdysozoa | Arthropoda | Crustacea
## 4 Animalia | Chordata | Actinopterygii | Perciformes | Carangidae | Chloroscombrus | Chloroscombrus chrysurus
## 5 Animalia | Chordata | Actinopterygii | Perciformes | Carangidae | Chloroscombrus | Chloroscombrus chrysurus
## 6 Animalia | Chordata | Actinopterygii | Perciformes | Carangidae | Chloroscombrus | Chloroscombrus chrysurus
##   source_specimen_life_stage interaction_type target_taxon_external_id
## 1                         NA     preyedUponBy               EOL:208553
## 2                         NA     preyedUponBy               EOL:208553
## 3                         NA     preyedUponBy               EOL:208553
## 4                         NA          preysOn                    EOL:1
## 5                         NA          preysOn                 EOL:1905
## 6                         NA          preysOn              EOL:2598871
##          target_taxon_name
## 1 Chloroscombrus chrysurus
## 2 Chloroscombrus chrysurus
## 3 Chloroscombrus chrysurus
## 4                 Animalia
## 5           Actinopterygii
## 6                Crustacea
##                                                                                             target_taxon_path
## 1 Animalia | Chordata | Actinopterygii | Perciformes | Carangidae | Chloroscombrus | Chloroscombrus chrysurus
## 2 Animalia | Chordata | Actinopterygii | Perciformes | Carangidae | Chloroscombrus | Chloroscombrus chrysurus
## 3 Animalia | Chordata | Actinopterygii | Perciformes | Carangidae | Chloroscombrus | Chloroscombrus chrysurus
## 4                                                                                                    Animalia
## 5                          Biota | Animalia | Chordata | Vertebrata | Gnathostomata | Pisces | Actinopterygii
## 6                                     Animalia | Bilateria | Protostomia | Ecdysozoa | Arthropoda | Crustacea
##   target_specimen_life_stage latitude longitude study_title
## 1                         NA       NA        NA          NA
## 2                         NA       NA        NA          NA
## 3                         NA       NA        NA          NA
## 4                         NA       NA        NA          NA
## 5                         NA       NA        NA          NA
## 6                         NA       NA        NA          NA
```

To see all locations for which interactions have entered in GloBi, 


```r
areas <- get_interaction_areas()
head(areas)
```

```
##   Latitude Longitude
## 1 29.34695 -92.98061
## 2 29.03260 -92.28701
## 3 28.03673 -96.11108
## 4 27.62409 -95.77403
## 5 26.33146 -96.03294
## 6 30.25024 -86.13114
```
You can also restrict this search to a certain area:


```r
areas <- get_interaction_areas (bbox=c(-67.87,12.79,-57.08,23.32))
head(areas)
```

```
##      Latitude Longitude
## 1251 18.24829 -66.49989
## 1255 18.33884 -65.82572
## 1268 18.06667 -63.06667
## 1339 13.16667 -59.53333
```

#### Custom Queries using Cypher

Currently, GloBI is powered by [neo4j](http://neo4j.org), a graph database. Most results created by the R functions provided by `rglobi` are basically wrappers for commonly used cypher queries. This way, you can still use the library without having to learn Cypher, a graph query language. However, if you feel adventurous, you can execute queries using Cypher only using the function `query()`:


```r
result <- query("START taxon = node:taxons(name='Homo sapiens') MATCH taxon<-[:CLASSIFIED_AS]-predator-[:ATE|PREYS_ON]->prey-[:CLASSIFIED_AS]->preyTaxon RETURN distinct(preyTaxon.name) as `prey.taxon.name` LIMIT 10")
result
```

```
##         prey.taxon.name
## 1          Homo sapiens
## 2    Macropus bernardus
## 3             Mysticeti
## 4   Dendrolagus scottae
## 5   Engraulis japonicus
## 6         Sergia lucens
## 7       Martes melampus
## 8   Lophocebus albigena
## 9            Salmonidae
## 10 Moschus chrysogaster
```

This particular query returns the names of the first ten known prey or diet items of humans (_Homo sapiens_). More examples of GloBI specific Cypher queries can be found on the [GloBI Cypher Wiki](https://github.com/globalbioticinteractions/globalbioticinteractions/wiki/cypher). 
