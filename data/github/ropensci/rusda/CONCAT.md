rusda
=======



[![Build Status](https://api.travis-ci.org/ropensci/rusda.png)](https://travis-ci.org/ropensci/rusda)
[![codecov.io](https://codecov.io/github/ropensci/rusda/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rusda?branch=master)

## Interface to USDA databases

## Description

An interface to the web service methods provided by the United States Department of Agriculture (USDA). The Agricultural Research Service (ARS) provides a large set of databases. The current version of the package holds interfaces to the Systematic Mycology and Microbiology Laboratory (SMML), which consists of four databases: Fungus-Host Distributions, Specimens, Literature and the Nomenclature database. It provides functions for querying these databases. The main function is `associations()`, which allows searching for fungus-host combinations.

## Get rusda

From CRAN (not on CRAN yet)


```r
install.packages("rusda")
```

or from rOpenSci


```r
install.packages("devtools")
devtools::install_github("ropensci/rusda")
```

And load rusda


```r
library("rusda")
```

## Example 1
In the following case we want to search for the host of a fungus (Rosellinia ligniaria) and the fungal associations for a give host (Fagus sylvatica). From our expert knowledge we already know that these two species are associated. Let's see if the Fungus-Hosts Distributions database from the USDA confirms this hypothesis.
We first specify the input species vectors. In this case they they are only one element long (they are allowed to be longer, too).


```r
host <- "Fagus sylvatica"
fungus <- "Rosellinia aquila"
```

Then we search for the associations and look at the output. Since we are interested in "real" species names we choose the 'clean' output. We want the full possible list of associations so we choose 'synonyms=TRUE'.


```r
fungi <- associations(host, database = "both", clean = TRUE, syn_include = TRUE, spec_type = "plant", process = TRUE)
hosts <- associations(fungus, database = "both", clean = TRUE, syn_include = TRUE, spec_type = "fungus", process = TRUE)

head(fungi$association$`Fagus sylvatica`)
```

```
#> [1] "Absidia glauca"            "Acia stenodon"            
#> [3] "Acrogenospora megalospora" "Actinocladium rhodosporum"
#> [5] "Actinonema fagicola"       "Alternaria alternata"
```

```r
head(hosts$association$`Rosellinia aquila`)
```

```
#> [1] "Acer pseudoplatanus" "Acer rubrum"         "Acer sp."           
#> [4] "Alnus incana"        "Alnus rubra"         "Asclepias sp."
```

Now we want to check if our initial knowledge is correct:


```r
is.element("Rosellinia aquila", fungi$associations[[1]])
```

```
#> [1] TRUE
```

```r
is.element("Fagus sylvatica", hosts$association[[1]])
```

```
#> [1] TRUE
```

We were right. Now we can be happy and search for other associations ;-).

## Example 2
We want to know the mean number of associations for a group, e.g. the Polyporales. Lets create a species input vector with Linnean species names derived from GenBank. In a first step we might want to check how many species are deposited in the database.


```r
polyporus <- c("Polyporus_admirabilis", "Polyporus_alveoaris", "Polyporus_americanus", "Polyporus_arcularius", "Polyporus_brumalis", "Polyporus_chozeniae", "Polyporus_ciliatus", "Polyporus_corylinus", "Polyporus_craterellus", "Polyporus_dictyopus", "Polyporus_favescens", "Polyporus_fraxineus", "Polyporus_gayanus", "Polyporus_grammocephalus", "Polyporus_guianensis", "Polyporus_lepideus", "Polyporus_leprieurii", "Polyporus_leptocephalus", "Polyporus_longiporus", "Polyporus_melanopus", "Polyporus_meridionalis", "Polyporus_pinsitus", "Polyporus_pseudobetulinus", "Polyporus_radicatus", "Polyporus_rhizophilus", "Polyporus_squamosus", "Polyporus_squamulosus", "Polyporus_submelanopus", "Polyporus_subvarius", "Polyporus_tenuiculus", "Polyporus_tessellatus", "Polyporus_tricholoma", "Polyporus_tuberaster", "Polyporus_tubiformis", "Polyporus_udus", "Polyporus_umbellatus", "Polyporus_varius", "Polyporus_virgatus")

poly_meta <- meta_smml(polyporus, process = TRUE, spec_type = "fungus")
head(poly_meta)
```

```
#>                       Nomenclature Specimens Host_Fungus Literature
#> Polyporus_admirabilis            1         1           1          1
#> Polyporus_alveoaris              0         0           0          0
#> Polyporus_americanus             0         0           0          0
#> Polyporus_arcularius             1         1           1          1
#> Polyporus_brumalis               1         1           1          1
#> Polyporus_chozeniae              0         0           0          0
```

Do we need to query associations for all Polyporus species?

```
length(polyporus)                             # 38 all species
nrow(poly_meta[rowSums(poly_meta)>0,])        # 27 with data species
```

No, only 27 species are supported with data...


```r
polyporus <- rownames(poly_meta[rowSums(poly_meta) > 0, ])
polyporus_ass <- associations(polyporus, database = "both", clean = TRUE, syn_include = TRUE, spec_type = "fungus", process = TRUE)
mean(unlist(lapply(polyporus_ass[[2]], length)))
```

```
#> [1] 23.25926
```

So within the genus Polyporus the mean number of host associations is:

[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
rusda 1.0.8
==============

### IMROVEMENTS
* now the user can also provide a genus name or multiple of them. Then Linnean Species names are downloaded from the NCBI taxonomy and used as query input.
* added tests, e.g. if USDA website is available
* updated vignetterusda
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  warning = FALSE,
  message = FALSE
)
```

[![Build Status](https://api.travis-ci.org/ropensci/rusda.png)](https://travis-ci.org/ropensci/rusda)
[![codecov.io](https://codecov.io/github/ropensci/rusda/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rusda?branch=master)

## Interface to USDA databases

## Description

An interface to the web service methods provided by the United States Department of Agriculture (USDA). The Agricultural Research Service (ARS) provides a large set of databases. The current version of the package holds interfaces to the Systematic Mycology and Microbiology Laboratory (SMML), which consists of four databases: Fungus-Host Distributions, Specimens, Literature and the Nomenclature database. It provides functions for querying these databases. The main function is `associations()`, which allows searching for fungus-host combinations.
## Get rusda
From CRAN
```r
install.packages("rusda")
```
or from rOpenSci
```r
install.packages("devtools")
library("devtools")

install_github("ropensci/rusda")
```
And load rusda
```r
library("rusda")
```
## Example 1
In the following example, we want to search for the host of a fungus (Rosellinia ligniaria) and the fungal associations for a give host (Fagus sylvatica). From our expert knowledge, we already know that these two species are associated. Let us see, wheather the USDA Fungus-Hosts Distributions database confirms this knowledge. 
We first specify the input species vectors. In this case they they are only one element long (they are allowed to be longer, too).
```r
host <- "Fagus sylvatica"
fungus <- "Rosellinia aquila"
```
Then we search for the associations and look at the output. Since we are interested in "real" species names we choose the 'clean' output. We want the full possible list of associations so we choose 'synonyms=TRUE'. 
```r
fungi <- associations(x = host, database = "both", clean = TRUE, syn_include = TRUE, spec_type = "plant", process = TRUE)
hosts <- associations(x = fungus, database = "both", clean = TRUE, syn_include = TRUE, spec_type = "fungus", process = TRUE)

head(fungi$association$`Fagus sylvatica`)
head(hosts$association$`Rosellinia aquila`)
```

```
## > head(fungi$association$`Fagus sylvatica`)
## [1] "Polyporus squamosus"       "Absidia glauca"            "Acia stenodon"            
## [4] "Acrogenospora megalospora" "Actinocladium rhodosporum" "Actinonema fagicola"

## > head(hosts$association$`Rosellinia aquila`)
## [1] "Acer pseudoplatanus" "Acer rubrum"         "Acer sp."            "Alnus incana"       
## [5] "Alnus rubra"         "Asclepias sp." 
```

Now we want to check if our initial knowledge is correct:
```r
cat("Is R. aqulia a fungus growing on F. sylvatica? \n", is.element("Rosellinia aquila", pathogens$association[[1]]))
cat("Is F. sylvatica a host for R. aqulia?\n", is.element("Fagus sylvatica", hosts$association[[1]]))
```

```
## > cat("Is R. aqulia a fungus growing on F. sylvatica? \n", is.element("Rosellinia aquila",
## pathogens$association[[1]]))
## Is R. aqulia a fungus growing on F. sylvatica? 
##  TRUE
## > cat("Is F. sylvatica a host for R. aqulia?\n", is.element("Fagus sylvatica", hosts$association[[1]]))
## Is F. sylvatica a host for R. aqulia?
##  TRUE
```
Our expert knowledge is excellent.

## Example 2
We want to know the mean number of associations for a group, e.g. the Polyporales. Lets create a species input vector with Linnean species names derived from GenBank. In a first step we might want to check how many species are deposited in the database.
```r
polyporus <- c("Polyporus_admirabilis", "Polyporus_alveoaris", "Polyporus_americanus", "Polyporus_arcularius", "Polyporus_brumalis", "Polyporus_chozeniae", "Polyporus_ciliatus", "Polyporus_corylinus", "Polyporus_craterellus", "Polyporus_dictyopus", "Polyporus_favescens", "Polyporus_fraxineus", "Polyporus_gayanus", "Polyporus_grammocephalus", "Polyporus_guianensis", "Polyporus_lepideus", "Polyporus_leprieurii", "Polyporus_leptocephalus", "Polyporus_longiporus", "Polyporus_melanopus", "Polyporus_meridionalis", "Polyporus_pinsitus", "Polyporus_pseudobetulinus", "Polyporus_radicatus", "Polyporus_rhizophilus", "Polyporus_squamosus", "Polyporus_squamulosus", "Polyporus_submelanopus", "Polyporus_subvarius", "Polyporus_tenuiculus", "Polyporus_tessellatus", "Polyporus_tricholoma", "Polyporus_tuberaster", "Polyporus_tubiformis", "Polyporus_udus", "Polyporus_umbellatus", "Polyporus_varius", "Polyporus_virgatus")

poly_meta <- meta_smml(x = polyporus, spec_type = "fungus", process = TRUE)
head(poly_meta)
```
```
## > head(poly_meta)
##                      Nomenclature Specimens Host_Fungus Literature
## Polyporus_admirabilis            1         1           1          1
## Polyporus_alveoaris              0         0           0          0
## Polyporus_americanus             0         0           0          0
## Polyporus_arcularius             1         1           1          1
## Polyporus_brumalis               1         1           1          1
## Polyporus_chozeniae              0         0           0          0
```
If we are interested in associations for this group, do we need to query associations for all Polyporus species? We can check that by firest running meta_smml and pruning the species without host data.

```
length(polyporus)                             # 38 all species
nrow(poly_meta[rowSums(poly_meta)>0,])        # 27 with data species
```
No, 27 of 38 species are supported with data ...
```r
polyporus <- rownames(poly_meta[rowSums(poly_meta)>0,])
polyporus_ass <- associations(x = polyporus, database = "both", clean = TRUE, syn_include = TRUE,
spec_type = "fungus", process = TRUE)
cat("Mean of hosts: ", mean(unlist(lapply(polyporus_ass[[2]], length))))
```
So within the genus Polyporus the mean number of host associations is:
```
## > cat("Mean of hosts: ", mean(unlist(lapply(polyporus_ass[[2]], length))))
## Mean of hosts:  22.88889
```


[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "rusda"
author: "Franz-Sebastian Krah"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Interface to USDA databases

## Description
An interface to the web service methods provided by the United States Department of Agriculture (USDA). The Agricultural Research Service (ARS) provides a large set of databases. The current version of the package holds interfaces to the Systematic Mycology and Microbiology Laboratory (SMML), which consists of four databases: Fungus-Host Distributions, Specimens, Literature and the Nomenclature database. It provides functions for querying these databases. The main function is \code{associations}, which allows searching for fungus-host combinations.

## Get rusda
From CRAN
```r
install.packages("rusda")
```
or from rOpenSci
```r
install.packages("devtools")
library("devtools")

install_github("ropensci/rusda")
```
And load rusda
```r
library("rusda")
```
## Example 1
In the following example, we want to search for the host of a fungus (Rosellinia ligniaria) and the fungal associations for a give host (Fagus sylvatica). From our expert knowledge, we already know that these two species are associated. Let us see, wheather the USDA Fungus-Hosts Distributions database confirms this knowledge. 
We first specify the input species vectors. In this case they they are only one element long (they are allowed to be longer, too).
```r
host <- "Fagus sylvatica"
fungus <- "Rosellinia aquila"
```
Then we search for the associations and look at the output. Since we are interested in "real" species names we choose the 'clean' output. We want the full possible list of associations so we choose 'synonyms=TRUE'. 
```r
fungi <- associations(x = host, database = "both", clean = TRUE, syn_include = TRUE, spec_type = "plant", process = TRUE)
hosts <- associations(x = fungus, database = "both", clean = TRUE, syn_include = TRUE, spec_type = "fungus", process = TRUE)

head(fungi$association$`Fagus sylvatica`)
head(hosts$association$`Rosellinia aquila`)
```

```
## > head(fungi$association$`Fagus sylvatica`)
## [1] "Polyporus squamosus"       "Absidia glauca"            "Acia stenodon"            
## [4] "Acrogenospora megalospora" "Actinocladium rhodosporum" "Actinonema fagicola"

## > head(hosts$association$`Rosellinia aquila`)
## [1] "Acer pseudoplatanus" "Acer rubrum"         "Acer sp."            "Alnus incana"       
## [5] "Alnus rubra"         "Asclepias sp." 
```

Now we want to check if our initial knowledge is correct:
```r
cat("Is R. aqulia a fungus growing on F. sylvatica? \n", is.element("Rosellinia aquila", pathogens$association[[1]]))
cat("Is F. sylvatica a host for R. aqulia?\n", is.element("Fagus sylvatica", hosts$association[[1]]))
```

```
## > cat("Is R. aqulia a fungus growing on F. sylvatica? \n", is.element("Rosellinia aquila",
## pathogens$association[[1]]))
## Is R. aqulia a fungus growing on F. sylvatica? 
##  TRUE
## > cat("Is F. sylvatica a host for R. aqulia?\n", is.element("Fagus sylvatica", hosts$association[[1]]))
## Is F. sylvatica a host for R. aqulia?
##  TRUE
```
Our expert knowledge is excellent.

## Example 2
We want to know the mean number of associations for a group, e.g. the Polyporales. Lets create a species input vector with Linnean species names derived from GenBank. In a first step we might want to check how many species are deposited in the database.
```r
polyporus <- c("Polyporus_admirabilis", "Polyporus_alveoaris", "Polyporus_americanus", "Polyporus_arcularius", "Polyporus_brumalis", "Polyporus_chozeniae", "Polyporus_ciliatus", "Polyporus_corylinus", "Polyporus_craterellus", "Polyporus_dictyopus", "Polyporus_favescens", "Polyporus_fraxineus", "Polyporus_gayanus", "Polyporus_grammocephalus", "Polyporus_guianensis", "Polyporus_lepideus", "Polyporus_leprieurii", "Polyporus_leptocephalus", "Polyporus_longiporus", "Polyporus_melanopus", "Polyporus_meridionalis", "Polyporus_pinsitus", "Polyporus_pseudobetulinus", "Polyporus_radicatus", "Polyporus_rhizophilus", "Polyporus_squamosus", "Polyporus_squamulosus", "Polyporus_submelanopus", "Polyporus_subvarius", "Polyporus_tenuiculus", "Polyporus_tessellatus", "Polyporus_tricholoma", "Polyporus_tuberaster", "Polyporus_tubiformis", "Polyporus_udus", "Polyporus_umbellatus", "Polyporus_varius", "Polyporus_virgatus")

poly_meta <- meta_smml(x = polyporus, spec_type = "fungus", process = TRUE)
head(poly_meta)
```
```
## > head(poly_meta)
##                      Nomenclature Specimens Host_Fungus Literature
## Polyporus_admirabilis            1         1           1          1
## Polyporus_alveoaris              0         0           0          0
## Polyporus_americanus             0         0           0          0
## Polyporus_arcularius             1         1           1          1
## Polyporus_brumalis               1         1           1          1
## Polyporus_chozeniae              0         0           0          0
```
If we are interested in associations for this group, do we need to query associations for all Polyporus species? We can check that by firest running meta_smml and pruning the species without host data.

```
length(polyporus)                             # 38 all species
nrow(poly_meta[rowSums(poly_meta)>0,])        # 27 with data species
```
No, 27 of 38 species are supported with data ...
```r
polyporus <- rownames(poly_meta[rowSums(poly_meta)>0,])
polyporus_ass <- associations(x = polyporus, database = "both", clean = TRUE, syn_include = TRUE,
spec_type = "fungus", process = TRUE)
cat("Mean of hosts: ", mean(unlist(lapply(polyporus_ass[[2]], length))))
```
So within the genus Polyporus the mean number of host associations is:
```
## > cat("Mean of hosts: ", mean(unlist(lapply(polyporus_ass[[2]], length))))
## Mean of hosts:  22.88889
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/literature.R
\name{literature}
\alias{literature}
\title{Downloads literature from SMML Literature DB}
\usage{
literature(x, spec_type = c("plant", "fungus"), process = TRUE)
}
\arguments{
\item{x}{a vector of class \code{character} containing fungal or plant species names}

\item{spec_type}{a character string specifying the type of \code{spec}. Can be either 
\code{"plant"} or \code{"fungus"}}

\item{process}{logical, if \code{TRUE} downloading and extraction process is displayed

an object of class \code{list}}
}
\value{
a vector of mode \code{list} with literature entries for \code{x}
}
\description{
Searches and downloads literature entries from the SMML Literature database
}
\examples{
\dontrun{
x <- "Polyporus badius"
lit <- literature(x, process = TRUE, spec_type = "fungus")
lit
}
}
\author{
Franz-Sebastian Krah
}
\name{rusda-internal}

\alias{clean_step}
\alias{getCOND}
\alias{getHF}
\alias{getMETA}
\alias{getSYNS}
\alias{ncbiSpecies}

\title{Internal rusda Functions}
\description{Internal \pkg{rusda} functions.}
\note{
These are internal functions in \code{rusda} and are not intended to be called by the user.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/associations.R
\name{associations}
\alias{associations}
\title{Downloads associations for input species from SMML Fungus-Host DB}
\usage{
associations(x, database = c("FH", "SP", "both"),
  spec_type = c("plant", "fungus"), clean = TRUE, syn_include = TRUE,
  process = TRUE, db = "ncbi")
}
\arguments{
\item{x}{a vector of class \code{character} containing fungal or plant species names or a genus name (see Details)}

\item{database}{a character string specifying the databases that should be queried. Valid are
\code{"FH"} (Fungus-Host Distributions), \code{"SP"} (Specimens) or \code{"both"} databases}

\item{spec_type}{a character string specifying the type of \code{x}. 
Can be either \code{"plant"} or \code{"fungus"}}

\item{clean}{logical, if \code{TRUE} a cleaning step is run of the resulting associations list}

\item{syn_include}{logical, if \code{TRUE} associations for synonyms are searched and added. For a
complete synonyms list check \code{rusda::synonyms}}

\item{process}{logical, if \code{TRUE} downloading and extraction process is displayed}

\item{db}{if x is higher than species level, all species for the higher taxon are retrived using the function taxize::downstream. Here one of ITIS (itis), Catalogue of Life (col), GBIF (gbif), or NCBI (ncbi) has to be selected. NCBI is default.}
}
\value{
an object of class \code{list}.

First is synonyms, second is associations. Synonmys is a
vector of mode \code{list} with synonyms for \code{x}. Notice: This is not a
complete list of synonym data in the database. This is the list of synonyms that contain data for
the input \code{x}. For a complete synonyms list check \code{rusda::synonyms} or (if needed) for fungi R package rmycobank.

Associations is a vector of mode \code{list} of associations for \code{x}
}
\description{
Searches and downloads associations from SMML Fungus-Hosts Distributions and Specimens database
for fungus or plant species input vector
}
\details{
The Fungus-Hosts distributions database 'FH' comprises data compiled from Literature. In
the uncleaned output all kinds of unspecified substrates are documented like "submerged wood".
Cleanded data displayes Linnean names only and species names with either "subsp.","f. sp." "f.",
"var.". The Specimens database comprises entries from field collections.

If genera names are supplied, then species are derived from the NCBI taxonomy.
}
\examples{
\dontrun{
## Example for species name(s) as input
x <- "Fagus sylvatica"
pathogens <- associations(x, database = "both", clean = TRUE, syn_include = TRUE,
spec_type = "plant", process = TRUE)
x <- "Rosellinia ligniaria"
hosts <- associations(x, database = "both", clean = TRUE, syn_include = TRUE, 
spec_type = "fungus", process = TRUE)
is.element("Rosellinia ligniaria", pathogens$association[[1]])
is.element("Fagus sylvatica", hosts$association[[1]])

## Example for genus/genera name(s) as input
x <- "Zehneria"
# or
x <- c("Zehneria", "Momordica")
hosts <- associations(x, database = "both", clean = TRUE, syn_include = TRUE, 
spec_type = "plant", process = TRUE)
}
}
\author{
Franz-Sebastian Krah
}
\name{rusda-package}
\alias{rusda-package}
%\alias{rusda}
\docType{package}
\title{
Interface to USDA Databases
}
\description{An interface to the web service methods provided by the United States Department of Agriculture (USDA). The Agricultural Research Service (ARS) provides a large set of databases. The current version of the package holds interfaces to the Systematic Mycology and Microbiology Laboratory (SMML), which consists of four databases: Fungus-Host Distributions, Specimens, Literature and the Nomenclature database. It provides functions for querying these databases. The main function is \code{associations}, which allows searching for fungus-host combinations.}
\details{
\tabular{ll}{
Package: \tab rusda\cr
Type: \tab Package\cr
Version: \tab 1.0.7\cr
Date: \tab 2016-01-20\cr
}

}
\author{
Franz-Sebastian Krah \cr
Maintainer: Franz-Sebastian Krah <f.krah@mailbox.org>
}
\references{
Farr, D.F., & Rossman, A.Y. Fungal Databases, Systematic Mycology and Microbiology Laboratory, ARS, USDA

 \url{http://nt.ars-grin.gov/sbmlweb/fungi/databases.cfm}, 
 \url{http://www.usda.gov/wps/portal/usda/usdahome}

}
\keyword{ package }
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/substrate.R
\name{substrate}
\alias{substrate}
\title{Downloads substrate data from SMML Nomenclature DB}
\usage{
substrate(x, process = TRUE)
}
\arguments{
\item{x}{a vector of class \code{character} containing fungal or plant species names}

\item{process}{logical, if \code{TRUE} downloading and extraction process is displayed}
}
\value{
an object of mode \code{list} containing substrate for fungus species
}
\description{
Searches and downloads substrate data from SMML Nomenclature database
}
\details{
Don't be disappointed. Not much data there. 
But depends on the study group, so give it try.
}
\examples{
\dontrun{
x <- c("Polyporus_rhizophilus", "Polyporus_squamosus")
subs.poly <- substrate(x, process=TRUE)
subs.poly
}
}
\author{
Franz-Sebastian Krah
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_smml.R
\name{meta_smml}
\alias{meta_smml}
\title{Downloads and evaluate species presence in SMML DBs}
\usage{
meta_smml(x, spec_type = c("plant", "fungus"), process = TRUE)
}
\arguments{
\item{x}{a vector of class \code{character} containing fungal or plant species or genus names}

\item{spec_type}{a character string specifying the type of \code{x}. 
Can be either \code{"plant"} or \code{"fungus"}}

\item{process}{logical, if \code{TRUE} downloading and extraction process is displayed}
}
\value{
an object of class \code{data.frame}: presence/absence
}
\description{
Searches, downloads and evaluates presence/absence of data in the SMML databases
}
\details{
Use this function before deriving data from one of the databases in order to prune your
input species vector. With pruned species vectors the functions will run faster. This is important
if \code{x} is some hundred species long.
}
\examples{
\dontrun{
fungus.meta <- meta_smml(x = "Picea abies", process = TRUE, spec_type = "plant")
fungus.meta
hosts.meta <- meta_smml(x = "Antrodiella citrinella", process = TRUE, spec_type = "fungus")
hosts.meta
}
}
\author{
Franz-Sebastian Krah
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synonyms_smml.R
\name{synonyms_smml}
\alias{synonyms_smml}
\title{Downloads synonym data from SMML Nomenclature DB}
\usage{
synonyms_smml(x, spec_type = c("plant", "fungus"), clean = TRUE,
  process = TRUE)
}
\arguments{
\item{x}{a vector of class \code{character} containing fungal or plant species or genus names}

\item{spec_type}{a character string specifying the type of \code{x}. 
Can be either \code{"plant"} or \code{"fungus"}}

\item{clean}{logical, if \code{TRUE} a cleaning step is run of the resulting associations list}

\item{process}{logical, if \code{TRUE} downloading and extraction process is displayed}
}
\value{
an object of class \code{list} containing synonyms for \code{x}
}
\description{
Searches and downloads synonym data from SMML Nomenclature database
}
\examples{
\dontrun{
x <- "Solanum tuberosum"
synonyms_usda(x, spec_type = "plant", process = TRUE, clean = TRUE)
x <- c("Phytophthora infestans", "Polyporus badius")
synonyms_usda(x, spec_type = "fungus", process = TRUE, clean = TRUE)
}
}
\author{
Franz-Sebastian Krah
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getStudy.R
\name{getStudy}
\alias{getStudy}
\title{Downloads studies for input study IDs}
\usage{
getStudy(id, sep = "; ")
}
\arguments{
\item{id}{a single study id or a vector of study ids}

\item{sep}{seperator how to collapse output literature references}
}
\value{
an object of class \code{data.frame} with studies
}
\description{
Downloads studies for study IDs which are output of associations
}
\author{
Franz-Sebastian Krah
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getBPI.R
\name{getBPI}
\alias{getBPI}
\title{Downloads specimens records for input study BPI accession number}
\usage{
getBPI(BPI, sep = "; ")
}
\arguments{
\item{BPI}{a single study BPI or a vector of study BPIs}

\item{sep}{seperator how to collapse output literature references}
}
\value{
an object of class \code{data.frame} with studies
}
\description{
Downloads specimens records for input study BPI accession number
}
\author{
Franz-Sebastian Krah
}
