
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
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Introduction to the rglobi package}
-->

```{r echo=FALSE, eval=FALSE}
knitr::opts_chunk$set(
  comment = "#>", 
  warning = FALSE, 
  message = FALSE,
  cache = TRUE
)
```

rglobi vignette - an R interface to the aggregated biotic interaction data of GloBI (http://globalbioticinteractions.org)
======

### About the package

`rglobi` allows easy access to the GLoBI database by wrapping the existing API.  

********************

### Quick start

#### First, install `rglobi`


```{r}
install.packages("devtools")
devtools::install_github("ropensci/rglobi")
```

```{r}
library("rglobi")
```

Note that since rglobi is still pretty new, only the dev version of the library is available at this time. We hope to publish a version to CRAN once the library matures. For more information see [GitHub](https://github.com/ropensci/rglobi).

#### Find interactions involving a certain species

Determining which species interact with each other (and how and where) is a major focus of ecology.  The Global Biotic Interactions Database offers data on multiple interaction types, and the `rglobi` library offers several ways to specify and access these interactions. `get_interactions()` is the primary function used if data on a specific interaction type is required.  A focus taxon is entered (may be specified as "Genus species" or higher level (e.g., Genus, Family, Class)).

```{r}
hsapiens <- get_interactions(taxon = "Homo sapiens", interaction.type = "preysOn")
head(hsapiens)
```

`get_predators_of()` and `get_prey_of()` are wrappers for `get_interactions()` that remove the need to specify interaction type

```{r}
hsapiens <- get_prey_of("Homo sapiens")
head(hsapiens)
```

For a complete list of supported interaction types, 

```{r}
get_interaction_types()
```

For data sources in which type of interactions was not specified, the interaction is labeled "interacts_with".

If you wish to view all interactions instead of specific types (e.g., parasitism and predation instead of just one of the two), the `get_interactions_by_taxa()` function allows this. In addition, the function provides options to search for interactions between two groups (source and target taxa, see above table) and only find interactions in a certain region.

```{r}
rattus <- get_interactions_by_taxa(sourcetaxon = "Rattus")
head(rattus)
```

Only a source taxa need be identified, but you can also specify a target taxon and/or geographic boundary (Coordinates must be in decimal degrees ([EPSG:4326](http://spatialreference.org/ref/epsg/wgs-84/)) and correspond to the west, south, east, and northern boundaries (i.e., min x, min y, max x, max y).  

```{r}
aves_crustacea_northern_hemisphere <- get_interactions_by_taxa( sourcetaxon = "Aves", targettaxon = "Crustacea", bbox=c(-180, 0, 180, 90 ))
head(aves_crustacea_northern_hemisphere)
```

#### Find interactions in a geographic areas

Instead of a taxon-specific approach, users may also wish to gather information on all interactions occuring in a specific area.  For example, a group developing ecoystem models for the Gulf of Mexico may want to consider all the data from that region.  `rglobi` enables this type of search by allowing users to specify a rectangular bounding box. Coordinates must be in decimal degrees ([EPSG:4326](http://spatialreference.org/ref/epsg/wgs-84/)) and correspond to the west, south, east, and northern boundaries (i.e., min x, min y, max x, max y).  

```{r}
gulfinteractions <- get_interactions_in_area( bbox=c(-97.0, 17.5, -81, 31))
head(gulfinteractions)
```

To see all locations for which interactions have entered in GloBi, 

```{r}
areas <- get_interaction_areas()
head(areas)
```
You can also restrict this search to a certain area:

```{r}
areas <- get_interaction_areas (bbox=c(-67.87,12.79,-57.08,23.32))
head(areas)
```

#### Custom Queries using Cypher

Currently, GloBI is powered by [neo4j](http://neo4j.org), a graph database. Most results created by the R functions provided by `rglobi` are basically wrappers for commonly used cypher queries. This way, you can still use the library without having to learn Cypher, a graph query language. However, if you feel adventurous, you can execute queries using Cypher only using the function `query()`:

```{r}
result <- query("START taxon = node:taxons(name='Homo sapiens') MATCH taxon<-[:CLASSIFIED_AS]-predator-[:ATE|PREYS_ON]->prey-[:CLASSIFIED_AS]->preyTaxon RETURN distinct(preyTaxon.name) as `prey.taxon.name` LIMIT 10")
result
```

This particular query returns the names of the first ten known prey or diet items of humans (_Homo sapiens_). More examples of GloBI specific Cypher queries can be found on the [GloBI Cypher Wiki](https://github.com/globalbioticinteractions/globalbioticinteractions/wiki/cypher). 
---
title: Introduction to the rglobi package
author: Jorrit Poelen
vignette: >
    %\VignetteIndexEntry{introduction}
    %\VignetteEngine{knitr::knitr}
    %\VignetteEncoding{UTF-8}
---

`rglobi` allows easy access to biotic interactions via the Global Biotic Interactions (GLoBI, https://www.globalbioticinteractions.org/) database by wrapping the existing API.  

### Quick start

#### First, install `rglobi`

For released version, use:

```r
install.packages("rglobi")
```

For most recent development version, use:

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

#### Find citations using individual interaction records

By default, `get_interactions()` returns interaction records aggregated by source/target taxa and interaction type. This means that a single result can be based on a single claim, or many claims (e.g., many individual records claiming that sea otters (_Enhydra lutris_) eat (eats, [RO:0002470](http://purl.obolibrary.org/obo/RO_0002470)) Red Rock Crabs (_Cancer productus_)). In addition, the aggregated results exclude values specific to individual claims like latitude, longitude, study_citation, study_source_citation. This is why these and other values will show up as NA by default. 

To retrieve an exhaustive list of _all_, unaggregated, interaction records, set the ```returnobservations=T```. For instance, 

```r
sea_otter_prey <- get_interactions(taxon = "Enhydra lutris", interaction.type = "preysOn", returnobservations=T)
head(sea_otter_prey)
```

results in 

```
##  source_taxon_external_id      source_taxon_name
## 1        INAT_TAXON:133061 Enhydra lutris kenyoni
## 2        INAT_TAXON:117520  Enhydra lutris nereis
## 3        INAT_TAXON:117520  Enhydra lutris nereis
## 4        INAT_TAXON:117520  Enhydra lutris nereis
## 5        INAT_TAXON:117520  Enhydra lutris nereis
## 6        INAT_TAXON:117520  Enhydra lutris nereis
##                                                                           source_taxon_path
## 1 Animalia | Chordata | Mammalia | Carnivora | Mustelidae | Enhydra | Enhydra lutris kenyoni
## 2  Animalia | Chordata | Mammalia | Carnivora | Mustelidae | Enhydra | Enhydra lutris nereis
## 3  Animalia | Chordata | Mammalia | Carnivora | Mustelidae | Enhydra | Enhydra lutris nereis
## 4  Animalia | Chordata | Mammalia | Carnivora | Mustelidae | Enhydra | Enhydra lutris nereis
## 5  Animalia | Chordata | Mammalia | Carnivora | Mustelidae | Enhydra | Enhydra lutris nereis
## 6  Animalia | Chordata | Mammalia | Carnivora | Mustelidae | Enhydra | Enhydra lutris nereis
##  source_specimen_life_stage interaction_type target_taxon_external_id
## 1                       <NA>          preysOn             GBIF:2289226
## 2                       <NA>          preysOn         INAT_TAXON:47828
## 3                       <NA>          preysOn         INAT_TAXON:47828
## 4                       <NA>          preysOn             GBIF:2285679
## 5                       <NA>          preysOn                 GBIF:137
## 6                       <NA>          preysOn         INAT_TAXON:47826
##  target_taxon_name
## 1           Octopus
## 2         Cancridae
## 3         Cancridae
## 4           Mytilus
## 5          Bivalvia
## 6  Cancer productus
##                                                                                                                           target_taxon_path
## 1                                                                       Animalia | Mollusca | Cephalopoda | Octopoda | Octopodidae | Octopus
## 2 Animalia | Arthropoda | Crustacea | Multicrustacea | Malacostraca | Eucarida | Decapoda | Pleocyemata | Brachyura | Cancroidea | Cancridae
## 3 Animalia | Arthropoda | Crustacea | Multicrustacea | Malacostraca | Eucarida | Decapoda | Pleocyemata | Brachyura | Cancroidea | Cancridae
## 4                                                                            Animalia | Mollusca | Bivalvia | Mytilida | Mytilidae | Mytilus
## 5                                                                                                             Animalia | Mollusca | Bivalvia
## 6                                                    Animalia | Arthropoda | Malacostraca | Decapoda | Cancridae | Cancer | Cancer productus
##  target_specimen_life_stage latitude longitude
## 1                       <NA> 56.93835 -135.7816
## 2                       <NA> 35.39714 -120.9534
## 3                       <NA> 35.39714 -120.9534
## 4                       <NA> 36.95229 -121.6419
## 5                       <NA> 36.90766 -121.7446
## 6                       <NA> 36.99515 -121.6459
##                                                                                                                                                       study_citation
## 1        Donna Pomeroy. 2012. Enhydra lutris kenyoni prey id Octopus. iNaturalist.org. Accessed at <https://www.inaturalist.org/observations/1181454> on 25 Jul 2021.
## 2      Marisa Agarwal. 2014. Enhydra lutris nereis prey id Cancridae. iNaturalist.org. Accessed at <https://www.inaturalist.org/observations/1130763> on 25 Jul 2021.
## 3      Marisa Agarwal. 2014. Enhydra lutris nereis prey id Cancridae. iNaturalist.org. Accessed at <https://www.inaturalist.org/observations/1130763> on 25 Jul 2021.
## 4          James Maughn. 2015. Enhydra lutris nereis prey id Mytilus. iNaturalist.org. Accessed at <https://www.inaturalist.org/observations/1964813> on 25 Jul 2021.
## 5         James Maughn. 2014. Enhydra lutris nereis prey id Bivalvia. iNaturalist.org. Accessed at <https://www.inaturalist.org/observations/1119439> on 25 Jul 2021.
## 6 Donna Pomeroy. 2014. Enhydra lutris nereis prey id Cancer productus. iNaturalist.org. Accessed at <https://www.inaturalist.org/observations/799531> on 25 Jul 2021.
##                                                                                                                        study_source_citation
## 1 http://iNaturalist.org is a place where you can record what you see in nature, meet other nature lovers, and learn about the natural world.
## 2 http://iNaturalist.org is a place where you can record what you see in nature, meet other nature lovers, and learn about the natural world.
## 3 http://iNaturalist.org is a place where you can record what you see in nature, meet other nature lovers, and learn about the natural world.
## 4 http://iNaturalist.org is a place where you can record what you see in nature, meet other nature lovers, and learn about the natural world.
## 5 http://iNaturalist.org is a place where you can record what you see in nature, meet other nature lovers, and learn about the natural world.
## 6 http://iNaturalist.org is a place where you can record what you see in nature, meet other nature lovers, and learn about the natural world.
```

After introducing the ```returnobservations=T```, the results now include values for ```study_citation``` and ```study_source_citation```. Also, if available, values for individual properties of interaction records (e.g., ```latitude```, ```longitude```) will be included. 

Note that the arguments can also be used in related methods (e.g., ```get_interactions_by_taxa(sourcetaxon="Enhydra lutris", targettaxon="Cancer productus", interactiontype="eats", returnobservations=T)```). 

#### Listing supported interaction types

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

Only a source taxa need be identified, but you can also specify a target taxon and/or geographic boundary (Coordinates must be in decimal degrees ([EPSG:4326 aka WGS 84 or WGS 1984](https://en.wikipedia.org/wiki/World_Geodetic_System)) and correspond to the west, south, east, and northern boundaries (i.e., min x, min y, max x, max y).  


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

Instead of a taxon-specific approach, users may also wish to gather information on all interactions occurring in a specific area.  For example, a group developing ecoystem models for the Gulf of Mexico may want to consider all the data from that region.  `rglobi` enables this type of search by allowing users to specify a rectangular bounding box. Coordinates must be in decimal degrees ([EPSG:4326 aka WGS 84 or WGS 1984](https://en.wikipedia.org/wiki/World_Geodetic_System)) and correspond to the west, south, east, and northern boundaries (i.e., min x, min y, max x, max y).  


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

To see all locations for which interactions have entered in GloBI, 


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

Currently, GloBI is powered by [neo4j](https://neo4j.com/), a graph database. Most results created by the R functions provided by `rglobi` are basically wrappers for commonly used cypher queries. This way, you can still use the library without having to learn Cypher, a graph query language. However, if you feel adventurous, you can execute queries using Cypher only using the function `query()`:


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

#### Pagination

By default, the amount of results are limited. If you'd like to retrieve all results, you can used pagination. For instance, to retrieve parasitic interactions using pagination, you can use:


```r
otherkeys = list("limit"=10, "skip"=0)
first_page_of_ten <- get_interactions_by_type(interactiontype = c("hasParasite"), otherkeys = otherkeys)
otherkeys = list("limit"=10, "skip"=10)
second_page_of_ten <- get_interactions_by_type(interactiontype = c("hasParasite"), otherkeys = otherkeys)
```

To exhaust all available interactions, you can keep paging results until the size of the page is less than the limit (e.g., ```nrows(interactions) < limit```).

#### Selecting Taxonomic Schemes

GloBI integrations with various different taxonomies (see [Selecting specific taxonomies](https://github.com/globalbioticinteractions/globalbioticinteractions/wiki/api#selecting-specific-taxonomies)). 

To select a specific preferred taxonomy like the [Open Tree of Life](http://opentreeoflife.org) reference taxonomy to retrieve up to 10 parasitic interactions, use:

```r
parasitic_interactions <- get_interactions_by_type(interactiontype = c("hasParasite"), otherkeys = list("taxonIdPrefix"="OTT", "limit"=10))
```

Note that if no taxonomic mapping for the preferred taxonomy is available, another taxonomy will be selected.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_interactions}
\alias{get_interactions}
\title{Get Species Interaction from GloBI}
\usage{
get_interactions(taxon = "Homo sapiens", interaction.type = "preysOn",
  ...)
}
\arguments{
\item{taxon}{canonical scientific name of source taxon (e.g. Homo sapiens)}

\item{interaction.type}{the preferred interaction type (e.g. preysOn)}

\item{...}{list of options to configure GloBI API}
}
\value{
species interactions between source and target taxa
}
\description{
Get Species Interaction from GloBI
}
\examples{
\dontrun{
get_interactions("Homo sapiens", "preysOn")
get_interactions("Insecta", "parasiteOf")
}
}
\seealso{
Other interactions: \code{\link{get_child_taxa}},
  \code{\link{get_interaction_matrix}},
  \code{\link{get_interaction_table}},
  \code{\link{get_interaction_types}},
  \code{\link{get_interactions_by_taxa}},
  \code{\link{get_interactions_by_type}},
  \code{\link{get_predators_of}}, \code{\link{get_prey_of}}
}
\concept{interactions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_prey_of}
\alias{get_prey_of}
\title{Get a List of Prey for given Predator Taxon}
\usage{
get_prey_of(taxon = "Homo sapiens", ...)
}
\arguments{
\item{taxon}{scientific name of predator taxon. Can be any taxonomic rank (e.g. Homo sapiens, Animalia)}

\item{...}{list of named options to configure GloBI API}
}
\value{
list of recorded predator-prey interactions that involve the desired predator taxon
}
\description{
Get a List of Prey for given Predator Taxon
}
\examples{
\dontrun{
get_prey_of("Homo sapiens")
get_prey_of("Primates")
}
}
\seealso{
Other interactions: \code{\link{get_child_taxa}},
  \code{\link{get_interaction_matrix}},
  \code{\link{get_interaction_table}},
  \code{\link{get_interaction_types}},
  \code{\link{get_interactions_by_taxa}},
  \code{\link{get_interactions_by_type}},
  \code{\link{get_interactions}},
  \code{\link{get_predators_of}}
}
\concept{interactions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_interaction_areas}
\alias{get_interaction_areas}
\title{Find locations at which interactions were observed}
\usage{
get_interaction_areas(bbox = NULL, read_csv = read_csv_online, ...)
}
\arguments{
\item{bbox}{Coordinates in EPSG:4326 decimal degrees defining "left, bottom,
right, top" of bounding box}

\item{read_csv}{function used to find csv associated to query url, defaulting to online query method}

\item{...}{list of named options to configure GloBI API}
}
\value{
Returns data frame of coordinates
}
\description{
Returns all locations (latitude,longitude) of interactions in data
base or area specified in arguments
}
\examples{
\dontrun{
get_interaction_areas ()
get_interaction_areas (bbox=c(-67.87,12.79,-57.08,23.32))
}
}
\seealso{
Other areas: \code{\link{get_interactions_in_area}}
}
\concept{areas}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_interactions_by_type}
\alias{get_interactions_by_type}
\title{Get Species Interactions by Interaction Type from GloBI}
\usage{
get_interactions_by_type(interactiontype = c("interactsWith"), ...)
}
\arguments{
\item{interactiontype}{the requested interaction type (e.g. preysOn)}

\item{...}{list of options to configure GloBI API}
}
\value{
species interactions given provided interaction type(s)
}
\description{
Get Species Interactions by Interaction Type from GloBI
}
\examples{
\dontrun{
get_interactions_by_type(interactiontype = c("eats", "eatenBy"))
get_interactions_by_type(interactiontype = "parasiteOf")
}
}
\seealso{
Other interactions: \code{\link{get_child_taxa}},
  \code{\link{get_interaction_matrix}},
  \code{\link{get_interaction_table}},
  \code{\link{get_interaction_types}},
  \code{\link{get_interactions_by_taxa}},
  \code{\link{get_interactions}},
  \code{\link{get_predators_of}}, \code{\link{get_prey_of}}
}
\concept{interactions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_interactions_by_taxa}
\alias{get_interactions_by_taxa}
\title{Return interactions involving specific taxa}
\usage{
get_interactions_by_taxa(sourcetaxon, targettaxon = NULL,
  interactiontype = NULL, accordingto = NULL,
  showfield = c("source_taxon_external_id", "source_taxon_name",
  "source_taxon_path", "source_specimen_life_stage", "interaction_type",
  "target_taxon_external_id", "target_taxon_name", "target_taxon_path",
  "target_specimen_life_stage", "latitude", "longitude", "study_citation",
  "study_external_id", "study_source_citation"), otherkeys = NULL,
  bbox = NULL, returnobservations = F, opts = list(),
  read_csv = read_csv_online)
}
\arguments{
\item{sourcetaxon}{Taxa of interest (consumer, predator, parasite); may be
specified as "Genus species" or higher level (e.g., Genus, Family, Class).}

\item{targettaxon}{Taxa of interest (prey, host); may be specified as "Genus
species" or higher level (e.g., Genus, Family, Class)}

\item{interactiontype}{Interaction types of interest (prey, host); may be
specified as listed by get_interaction_types()}

\item{accordingto}{Data source of interest}

\item{showfield}{Data fields of interest (e. g. source_taxon_external_id, source_taxon_name);
may be specified as listed by get_data_fields()}

\item{otherkeys}{list of key-value pairs to query any field not covered by other parameters;
keys may be specified as listed by get_data_fields()}

\item{bbox}{Coordinates in EPSG:4326 decimal degrees defining "left, bottom,
right, top" of bounding box}

\item{returnobservations}{if true, all individual observations are returned,
else only distinct relationships}

\item{opts}{list of named options to configure GloBI API}

\item{read_csv}{function used to find csv associated to query url, defaulting to online query method}
}
\value{
Returns data frame of interactions
}
\description{
Returns interactions involving specific taxa.  Secondary (target)
taxa and spatial boundaries may also be set
}
\note{
For data sources in which type of interactions were not specified, the interaction is labeled "interacts_with"
}
\examples{
\dontrun{
get_interactions_by_taxa(sourcetaxon = "Rattus")
get_interactions_by_taxa(sourcetaxon = "Aves", targettaxon = "Rattus")
get_interactions_by_taxa(sourcetaxon = "Rattus rattus",
bbox = c(-67.87,12.79,-57.08,23.32))
}
}
\seealso{
Other interactions: \code{\link{get_child_taxa}},
  \code{\link{get_interaction_matrix}},
  \code{\link{get_interaction_table}},
  \code{\link{get_interaction_types}},
  \code{\link{get_interactions_by_type}},
  \code{\link{get_interactions}},
  \code{\link{get_predators_of}}, \code{\link{get_prey_of}}
}
\concept{interactions}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_interactions_in_area}
\alias{get_interactions_in_area}
\title{Return all interactions in specified area}
\usage{
get_interactions_in_area(bbox, ...)
}
\arguments{
\item{bbox}{Coordinates in EPSG:4326 decimal degrees defining "left, bottom,
right, top" of bounding box}

\item{...}{list of named options to configure GloBI API}
}
\value{
Returns data frame of interactions
}
\description{
Returns all interactions in data base in area specified in
arguments
}
\examples{
\dontrun{
get_interactions_in_area(bbox = c(-67.87, 12.79, -57.08, 23.32))
}
}
\seealso{
Other areas: \code{\link{get_interaction_areas}}
}
\concept{areas}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_interaction_matrix}
\alias{get_interaction_matrix}
\title{Get Interaction Matrix. Constructs an interaction matrix indicating whether source taxa (rows) or target taxa (columns) are known to interact with given type.}
\usage{
get_interaction_matrix(source.taxon.names = list("Homo sapiens"),
  target.taxon.names = list("Mammalia"), interaction.type = "eats",
  opts = list(), read_csv = read_csv_online)
}
\arguments{
\item{source.taxon.names}{list of source taxon names (e.g. list('Mammalia', 'Aves', 'Ariopsis felis'))}

\item{target.taxon.names}{list of target taxon names}

\item{interaction.type}{the preferred interaction type (e.g. preysOn)}

\item{opts}{list of options to configure GloBI API}

\item{read_csv}{function used to find csv associated to query url, defaulting to online query method}
}
\value{
matrix representing species interactions between source and target taxa
}
\description{
Get Interaction Matrix. Constructs an interaction matrix indicating whether source taxa (rows) or target taxa (columns) are known to interact with given type.
}
\examples{
\dontrun{
get_interaction_matrix("Homo sapiens", "Mammalia", "interactsWith")
}
}
\seealso{
Other interactions: \code{\link{get_child_taxa}},
  \code{\link{get_interaction_table}},
  \code{\link{get_interaction_types}},
  \code{\link{get_interactions_by_taxa}},
  \code{\link{get_interactions_by_type}},
  \code{\link{get_interactions}},
  \code{\link{get_predators_of}}, \code{\link{get_prey_of}}
}
\concept{interactions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_interaction_table}
\alias{get_interaction_table}
\title{Returns all known child taxa with known interaction of specified source and target taxa on any rank.}
\usage{
get_interaction_table(source.taxon.names = list(),
  target.taxon.names = list(), interaction.type = "preysOn",
  skip = 0, limit = 100, opts = list())
}
\arguments{
\item{source.taxon.names}{list of taxon names for source}

\item{target.taxon.names}{list of taxon names for target}

\item{interaction.type}{kind of interaction}

\item{skip}{number of records skipped before including record in result table, used in pagination}

\item{limit}{maximum number of interaction to include}

\item{opts}{connection parameters and other options}
}
\value{
table of matching source, target and interaction types
}
\description{
Returns all known child taxa with known interaction of specified source and target taxa on any rank.
}
\examples{
\dontrun{
get_interaction_table(source.taxon.names = list("Aves"), target.taxon.names = list('Insecta'))
}
}
\seealso{
Other interactions: \code{\link{get_child_taxa}},
  \code{\link{get_interaction_matrix}},
  \code{\link{get_interaction_types}},
  \code{\link{get_interactions_by_taxa}},
  \code{\link{get_interactions_by_type}},
  \code{\link{get_interactions}},
  \code{\link{get_predators_of}}, \code{\link{get_prey_of}}
}
\concept{interactions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_predators_of}
\alias{get_predators_of}
\title{Get a List of Predators of a Given Prey Taxon}
\usage{
get_predators_of(taxon = "Rattus rattus", ...)
}
\arguments{
\item{taxon}{scientific name of prey taxon. Can be any taxonomic rank (e.g. Rattus rattus, Decapoda)}

\item{...}{list of named options to configure the GloBI API}
}
\value{
list of recorded prey-predator interactions that involve the desired prey taxon.
}
\description{
Get a List of Predators of a Given Prey Taxon
}
\examples{
\dontrun{
get_predators_of("Rattus rattus")
get_predators_of("Primates")
}
}
\seealso{
Other interactions: \code{\link{get_child_taxa}},
  \code{\link{get_interaction_matrix}},
  \code{\link{get_interaction_table}},
  \code{\link{get_interaction_types}},
  \code{\link{get_interactions_by_taxa}},
  \code{\link{get_interactions_by_type}},
  \code{\link{get_interactions}}, \code{\link{get_prey_of}}
}
\concept{interactions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{query}
\alias{query}
\title{Executes a Cypher Query Against GloBI's Neo4j Instance}
\usage{
query(cypherQuery, opts = list())
}
\arguments{
\item{cypherQuery}{Cypher query (see http://github.com/globalbioticinteractions/globalbioticinteractions/wiki/cypher for examples)}

\item{opts}{list of named options to configure GloBI API}
}
\value{
result of cypher query string
}
\description{
Executes a Cypher Query Against GloBI's Neo4j Instance
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_child_taxa}
\alias{get_child_taxa}
\title{Returns all known child taxa with known interaction of specified taxa and rank.}
\usage{
get_child_taxa(taxon.names, rank = "Species", skip = 0, limit = 25,
  opts = list())
}
\arguments{
\item{taxon.names}{list of taxa of which child taxa should be included.}

\item{rank}{selected taxonomic rank of child taxa}

\item{skip}{number of child taxon names to skip before returning result. May be used for pagination.}

\item{limit}{maximum number of child taxon names returned}

\item{opts}{list of options including web service configuration like "port" and "host"}
}
\value{
list of child taxon names
}
\description{
Returns all known child taxa with known interaction of specified taxa and rank.
}
\examples{
\dontrun{
get_child_taxa(list("Aves"))
}
}
\seealso{
Other interactions: \code{\link{get_interaction_matrix}},
  \code{\link{get_interaction_table}},
  \code{\link{get_interaction_types}},
  \code{\link{get_interactions_by_taxa}},
  \code{\link{get_interactions_by_type}},
  \code{\link{get_interactions}},
  \code{\link{get_predators_of}}, \code{\link{get_prey_of}}
}
\concept{interactions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_interaction_types}
\alias{get_interaction_types}
\title{List interactions identified in GloBI database}
\usage{
get_interaction_types(opts = list(), read_csv = read_csv_online)
}
\arguments{
\item{opts}{list of named options to configure GloBI API}

\item{read_csv}{function used to find csv associated to query url, defaulting to online query method}
}
\value{
Returns data frame of supported interaction types
}
\description{
Returns data frame with supported interaction types
}
\examples{
\dontrun{
get_interaction_types()
}
}
\seealso{
Other interactions: \code{\link{get_child_taxa}},
  \code{\link{get_interaction_matrix}},
  \code{\link{get_interaction_table}},
  \code{\link{get_interactions_by_taxa}},
  \code{\link{get_interactions_by_type}},
  \code{\link{get_interactions}},
  \code{\link{get_predators_of}}, \code{\link{get_prey_of}}
}
\concept{interactions}
\keyword{database}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rglobi.R
\name{get_data_fields}
\alias{get_data_fields}
\title{List data fields identified in GloBI database}
\usage{
get_data_fields(opts = list(), read_csv = read_csv_online)
}
\arguments{
\item{opts}{list of named options to configure GloBI API}

\item{read_csv}{function used to find csv associated to query url, defaulting to online query method}
}
\value{
Returns data frame of supported data fields
}
\description{
Returns data frame with supported data fields
}
\examples{
\dontrun{
get_data_fields()
}
}
\concept{data}
\keyword{database}
