
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/Logo_FishPhyloMaker.png" alt="fish logo" width="200px" align="right"/>

# FishPhyloMaker

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/FishPhyloMaker)](https://cran.r-project.org/package=FishPhyloMaker)

[![DOI](https://zenodo.org/badge/336899540.svg)](https://zenodo.org/badge/latestdoi/336899540)

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://choosealicense.com/licenses/mit/)

`{FishPhyloMaker}` is an R package that allows to construct synthesis
phylogenies for finned-ray fishes. The package has two main functions,
`FishTaxaMaker` and `FishPhyloMaker`. The first generate a data frame
from fish species names provided by the user, checking the validity of
these names and possible synonyms by using the information contained in
[Fishbase database](http://www.fishbase.org) through the package
[rfishbase](https://CRAN.R-project.org/package=rfishbase). The output of
`FishTaxaMaker()` function is a list containing the following objects:

  - **All\_info\_fishbase**: A data frame containing the taxonomic
    classification of valid species accordingly to Fishbase;

  - **Taxon\_data\_FishPhyloMaker**: A data frame with three columns
    containing the valid scientific name of the species (s), its family
    (f) and order (o). This data frame can be used in `FishPhyloMaker()`
    to generate the phylogeny;

  - **Species\_not\_in\_Fishbase**: A character vector containing the
    names of species that was not found in Fishbase with a valid name
    data frame object containing three columns with the name of species
    (s), the Family (f) and the Order (o) of all species provided by the
    user.

Besides to help checking the validity of the names of species, its
synonyms and duplicated species, the data frame returned by
`FishTaxaMaker()` in **Taxon\_data\_FishPhyoMaker** are formatted so
that it can be directly used in the core function `FishPhyloMaker()`.
This function will use the information of the taxonomic hierarchy
contained in the data frame returned from `FishTaxaMaker()`, joint with
the information present in the the [fishtree of life
project](https://fishtreeoflife.org/) to construct the phylogenetic
tree.

The paper describing the functionalities of FishPhyloMaker is now
published in Ecological Informatics and can be assessed
[here](https://www.sciencedirect.com/science/article/pii/S1574954121002727).

# Download

A stable version of FishPhyloMaker can be installed from CRAN

``` r
install.packages("FishPhyloMaker")
```

To install the development version of this package the user must type:

``` r
# install.packages("devtools")
devtools::install_github("GabrielNakamura/FishPhyloMaker", ref = "main")
```

# Example

To run an example the user can load a data set present in the package:

``` r
library(FishPhyloMaker)
data(neotropical_comm)
data_comm <- neotropical_comm[, -c(1, 2)] # removing latitude and longitude
```

This data set comprises a community matrix with the occurrences of 59
fish species in headwater streams of Parana and Paraguai River Basins,
in Brazil. The coordinates of these streams are presented in the two
first columns of this data set.

First the user must obtain the data necessary to enter in
`FishPhyloMaker` using `FishTaxaMaker` function.

``` r
taxon_data <- FishTaxaMaker(data_comm, allow.manual.insert = TRUE)
Characidae
Characiformes
Characidae
Characiformes
Characidae
Characiformes
Loricariidae
Siluriformes
Cichlidae
Cichliformes
Crenuchidae
Characiformes
Gymnotidae
Gymnotiformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Heptapteridae
Siluriformes
Characidae
Characiformes
Loricariidae
Siluriformes
Characidae
Characiformes
```

`FishTaxaMaker` finds in [Fishbase](http://www.fishbase.org/search.php)
for the family and the order of species provided in data argument. If
any species was not find in Fishbase, the user will be asked to type the
Family and the Order of this species manually. This function can also be
useful to check possible misspelling errors in the name of species.

Finally run `FishPhyloMaker`

``` r
res_phylo <- FishPhyloMaker(data = taxon_data$Taxon_data_FishPhyloMaker,
                            insert.base.node = TRUE, 
                            return.insertions = TRUE, 
                            progress.bar = TRUE)
```

The species are inserted in a sequential procedure. Species with any
family representatives will be printed in the console jointly with a
list of Genus contained of its family contained in the three, so that
the user must choose. The user have three options:

1.  Insert **near to a specific Genus**: the user must type the name of
    this Genus;
2.  Insert **between to Genus**: the user must type the names of these
    two Genus separated by a white space;
3.  Insert **at the node that correspond to the Family** of the species
    being inserted: the user must type the name of the Family.

The output has two objects, a phylogenetic tree that can be directly
plot with the following code:

``` r
plot(res_phylo$Phylogeny, cex = 0.7)
```

And a data frame indicating at which level the species was inserted (one
of the six categories detailed above).

``` r
res_phylo$Insertions_data
```

For more details and updates see [FishPhyloMaker web
page](https://gabrielnakamura.github.io/FishPhyloMaker/)

# Next steps

  - [x] CRAN release
  - [ ] Implement user option to insert a taxonomic table and backbone
    phylogeny
  - [ ] Incorporate butterfly phylogeny
  - [ ] Automatic phylogeny plot
# FishPhyloMaker 0.2.0

* bug fixed in `FishPhyloMaker` function in family insertion
procedure.

* new argument in `FishPhyloMaker` function that allows to insert all species at node base of the family.

# FishPhyloMaker 0.1.2

* Second release of FishPhyloMaker - errors in insertions fixed

* A new function `whichFishAdd` was released

* A stable version of FishPhyloMaker is now on CRAN

* Function PD_deficit now outputs total and inserted phylogenetic information 

* Paper describing FishPhyloMaker published in [Ecological Informatics](https://www.sciencedirect.com/science/article/pii/S1574954121002727) ## Test environments

* local OS X install, R 4.1.0
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu Linux 20.04.1 LTS, R-release, GCC

## R CMD check results

There were no ERRORs or WARNINGs.

NOTES.

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Gabriel Nakamura <gabriel.nakamura.souza@gmail.com>’

  New submission

MANUAL CHECKING

* Please omit the redundant "An R Package to" from your title.
  
  I removed the redundant part of the title
  
* Is there a doi available to the reference in your description field you could add in the form <doi:10.prefix/suffix>?
  
  I added a doi to the reference
  
* Please unwrap the examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing.
(You could also replace \dontrun{} with \donttest, if it takes longer than 5 sec to be executed, but it would be preferable to have automatic checks for functions. Otherwise, you can also write some tests.)
  
  I remove  \dontrun{} and change it to \donttest{} in fishPhyloMaker.R and whichFishAdd.R functions, since the example are executed in > 5 sec. I can not remove \dontrun{} from tab_function.R since it depends on a interactive procedure with the user.  


* Possibly mis-spelled words in DESCRIPTION:
  al (19:91)
  et (19:88)
  Phylogenies (3:8)
  Rabosky (19:80)
  
The spelling is correct for all words


* Found the following (possibly) invalid URLs:
   URL: %22http://www.fishbase.org%22
     From: inst/doc/FishPhyloMaker_vignette.html
     Message: Invalid URI scheme
     
I changed the URL to [Fishbase database](https://www.fishbase.se) and now the link is correct with the proper https address

* Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1590/1982-0224-2020-0126
    From: man/neotropical_comm.Rd
    Status: Error
    Message: libcurl error code 60:
      	SSL certificate problem: unable to get local issuer certificate
      	(Status without verification: OK)

The URL is valid ---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

<img src="man/figures/Logo_FishPhyloMaker.png" alt="fish logo" width="200px" align="right"/>

# FishPhyloMaker

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/FishPhyloMaker)](https://cran.r-project.org/package=FishPhyloMaker)

[![DOI](https://zenodo.org/badge/336899540.svg)](https://zenodo.org/badge/latestdoi/336899540)

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://choosealicense.com/licenses/mit/)

`{FishPhyloMaker}` is an R package that allows to construct synthesis phylogenies for finned-ray fishes.
The package has two main functions, `FishTaxaMaker` and `FishPhyloMaker`. The first generate a data frame from fish species names provided by the user, checking the validity of these names and possible synonyms by using the information contained in [Fishbase database](http://www.fishbase.org) through the package [rfishbase](https://CRAN.R-project.org/package=rfishbase). The output of `FishTaxaMaker()` function is a list containing the following objects: 

- **All_info_fishbase**: A data frame containing the taxonomic classification of valid species accordingly to Fishbase;

- **Taxon_data_FishPhyloMaker**: A data frame with three columns containing the valid scientific name of the species
    (s), its family (f) and order (o). This data frame can be used in `FishPhyloMaker()` to generate the phylogeny;
- **Species_not_in_Fishbase**: A character vector containing the names of species that was not found      in Fishbase with a valid name
data frame object containing three columns with the name of species (s), the Family (f) and the Order (o) of all species provided by the user. 

Besides to help checking the validity of the names of species, its synonyms and duplicated species, the data frame returned by `FishTaxaMaker()` in **Taxon_data_FishPhyoMaker** are formatted so that it can be directly used in the core function `FishPhyloMaker()`. This function will use the information of the taxonomic hierarchy contained in the data frame returned from `FishTaxaMaker()`, joint with the information present in the the [fishtree of life project](https://fishtreeoflife.org/) to construct the phylogenetic tree.

The paper describing the functionalities of FishPhyloMaker is now published in Ecological Informatics and can be assessed [here](https://www.sciencedirect.com/science/article/pii/S1574954121002727).


# Download 

A stable version of FishPhyloMaker can be installed from CRAN

```{.r}
install.packages("FishPhyloMaker")
```

To install the development version of this package the user must type:

```{.r}
# install.packages("devtools")
devtools::install_github("GabrielNakamura/FishPhyloMaker", ref = "main")
```

# Example

To run an example the user can load a data set present in the package:

```{.r}
library(FishPhyloMaker)
data(neotropical_comm)
data_comm <- neotropical_comm[, -c(1, 2)] # removing latitude and longitude
```

This data set comprises a community matrix with the occurrences of 59 fish species in headwater streams of 
Parana and Paraguai River Basins, in Brazil. The coordinates of these streams are presented in the two first columns of this data set.

First the user must obtain the data necessary to enter in `FishPhyloMaker`  using `FishTaxaMaker` function.

```{.r}
taxon_data <- FishTaxaMaker(data_comm, allow.manual.insert = TRUE)
Characidae
Characiformes
Characidae
Characiformes
Characidae
Characiformes
Loricariidae
Siluriformes
Cichlidae
Cichliformes
Crenuchidae
Characiformes
Gymnotidae
Gymnotiformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Heptapteridae
Siluriformes
Characidae
Characiformes
Loricariidae
Siluriformes
Characidae
Characiformes
```

`FishTaxaMaker` finds in [Fishbase](http://www.fishbase.org/search.php) for the family and the order of species provided in data argument. If any species was not find in Fishbase, the user will be asked to type the Family and the Order of this species manually. This function can also be useful to check possible misspelling errors in the name of species.


Finally run `FishPhyloMaker`

```{.r}
res_phylo <- FishPhyloMaker(data = taxon_data$Taxon_data_FishPhyloMaker,
                            insert.base.node = TRUE, 
                            return.insertions = TRUE, 
                            progress.bar = TRUE)
```

The species are inserted in a sequential procedure. Species with any family representatives will be printed in the console jointly with a list of Genus contained of its family contained in the three, so that the user must choose. The user have three options: 

1. Insert **near to a specific Genus**: the user must type the name of this Genus;
2. Insert **between to Genus**: the user must type the names of these two Genus separated by a white space;
3. Insert **at the node that correspond to the Family** of the species being inserted: the user must type the name of the Family.

The output has two objects, a phylogenetic tree that can be directly plot with the following code:

```{.r}
plot(res_phylo$Phylogeny, cex = 0.7)
```

And a data frame indicating at which level the species was inserted (one of the six categories detailed above).

```{.r}
res_phylo$Insertions_data
```

For more details and updates see [FishPhyloMaker web page](https://gabrielnakamura.github.io/FishPhyloMaker/)

# Next steps

- [x] CRAN release
- [ ] Implement user option to insert a taxonomic table and backbone phylogeny
- [ ] Incorporate butterfly phylogeny
- [ ] Automatic phylogeny plot
---
title: "How to know which species must be added to phylogeny?"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{How to find for species that must be added?}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

In this article we will show how to find for species that must be added to the mega-tree, before running the `FishPhyloMaker` function. When there are few species that must be added, the interactive insertion process is ease to follow. However, when we are dealing with big data bases with many species, it is desirable that before to initiate the insertion process with `FishPhyloMaker` function we know which species must be inserted and at which position in phylogenetic tree. So, in order to facilitate this process the function `whichFishAdd` take all species informed by the user and provide a classification of each level in which the species be inserted when running `FishPhyloMaker` function. The code of classification follows similar classification from that presented in the insertion process:

-   **Present_in_tree** the species was already present in the original tree;
-   **Congeneric_insertion** species inserted as a sister species of the same genus presented in the tree;
-   **Family_insertion** if not found any congeneric species, the species will be inserted near to, or between genus of the same family presented in the tree. The user can also insert the species in the base of the family;
-   **Order_insertion** if not found any genus of the same family of the species that must be inserted, the species will be inserted near to a given family, between two families or at the node that corresponds to the Order of this species;


First, lets read `{FishPhyloMaker}` package and the data containing more than 650 fish species with occurrence in the Neotropical region.

```{r setup, eval=FALSE, echo=TRUE}
library(FishPhyloMaker)
data("fish_SAmerica")
```

Then, we will check for the names of the species in the data base through the function `FishTaxaMaker`, and inform manually the Family and the Order of the species that were not find in Fishbase. This procedure is better explained in other [articles](https://gabrielnakamura.github.io/FishPhyloMaker/articles/FishPhyloMaker_vignette.html)

```{r taxondata, eval=FALSE, echo=TRUE}
taxon_data <- FishTaxaMaker(data = fish_SAmerica)
Loricariidae
Siluriformes
Anostomidae 
Characiformes
Anostomidae 
Characiformes
Anostomidae 
Characiformes
Serrasalmidae
Characiformes
Auchenipteridae
Siluriformes
Auchenipteridae
Siluriformes
Auchenipteridae
Siluriformes
```

taxon_data contain the taxonomic information for all species, and can be used in `whichFishAdd` function to perform the classification of each species.

```{r addFish, eval=FALSE, echo=TRUE}
class_fish <- whichFishAdd(data = taxon_data) 
```

class_fish is a data frame that contain in insertion column an indication of which level each species must be inserted (or if the species are already on tree). This data can be used previous than `FishPhyloMaker` to find for phylogenetic position of the absent species.
---
title: "Using FishPhyloMaker to build phylogenies for finned-ray fishes"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using FishPhyloMaker to build phylogenies for finned-ray fishes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

In this article we will show how to obtain a phylogeny for a pool of species using FishPhyloMaker package. First of all, we need to install and read `{FishPhyloMaker}` package, that can be made using the following code:


```{r install_pkg, echo=TRUE, eval=FALSE}
devtools::install_github("GabrielNakamura/FishPhyloMaker", ref = "main")
```

We will use a data base that contain stream fish species from Paraná and Paraguai basins. These comprises two of the main basins in Neotropical region. This data are embedded in the package and can be easily accessed by typing:

```{r read_data, eval=FALSE, echo=TRUE} 
library(FishPhyloMaker)
data(neotropical_comm)
data_comm <- neotropical_comm[, -c(1, 2)] # removing Latitude and Longitude
```

`neotropical_comm` correspond to a data frame containing the abundance of 59 fish species collected in 20 headwater streams, being 10 in Paraná and 10 in Paraguay Basin. To obtain a phylogenetic tree for these species we first need to format these data accordingly to enter in `FishPhyloMaker` function, that correspond to a data frame with three column containing the taxonomic information of each one of the 59 species. To prepare this data we will use the function `FishTaxaMaker` that finds for each species name present in the data_comm and checks if these names are found in [Fishbase database](https://www.fishbase.se/). 
Note that not all species was found in Fishbase data, so, the user must type manually, through an interactive process, the name of the Family and the Order of these species that was not found. A message will appear in the command prompt indicating when to type the names. The names must be typed sequentially, as in the following code.

```{r taxon_data, echo=TRUE, eval=FALSE}
taxon_data <- FishTaxaMaker(data_comm, allow.manual.insert = TRUE)
Characidae
Characiformes
Characidae
Characiformes
Characidae
Characiformes
Loricariidae
Siluriformes
Characidae
Characiformes
Cichlidae
Cichliformes
Crenuchidae
Characiformes
Gymnotidae
Gymnotiformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Heptapteridae
Siluriformes
Characidae
Characiformes
Loricariidae
Siluriformes
Characidae
Characiformes
```

This function will return a data frame with three columns and the same number of species present in data_comm, and now taxon_data can be used in `FishPhyloMaker` function to generate the phylogenetic tree for these species. In this example we will set the argument return.insertions = TRUE, since we want to know at which level of hierarchical taxonomy each species was inserted in the process, the levels are explained in the [introductory article](https://gabrielnakamura.github.io/FishPhyloMaker/index.html). Like `FishTaxaMaker` function, `FishPhyloMaker` also have a interactive process of insertion in which the user must inform at which position species must be inserted in the mega-tree. There are three options: 

1. **Insert the species near to a specific genus/family**: To insert the species near to a specific genus/family,
    the user must type the name of the genus (or family) at which the species must be inserted;

2. **Insert the species between two genus/families**: To insert the species between two different genus/families, the user 
    must type the names of the two genus/families separated by a blank space;

3. **Insert at the root of the family/order**: To insert the species at the node that corresponds to the family (or order), 
    the user must type the name of the family.

For each species a list of names containing Genus names or Family names will appear at command prompt, and the user
    must decide among one of the previous option.

*NOTE* - The function can take several minutes to run depending on the number of species that must be inserted.

```{r phylo_make, eval=FALSE, echo=TRUE}
phylo_fish_streams <- FishPhyloMaker(data = taxon_data, 
                                     return.insertions = TRUE,
                                     insert.base.node = TRUE, 
                                     progress.bar = TRUE)
```

In this example I choose to insert all the species that did not presented congeneric relatives in the 
    base of its family.
    
The result of `FishPhyloMaker` is a list containing two objects, a phylogenetic tree and a data frame indicating at which level each species was inserted in the phylogeny.
We can use the result contained in phylo_fish_streams to plot the tree. We can use the following code to highlight the families in the tree.

```{r plot_phylo, eval=FALSE, echo=TRUE}
library(ggtree)
tree.PR<- phylo_fish_streams$Phylogeny

tree.PR <- ape::makeNodeLabel(tree.PR)
phylo <- tree.PR

rm.famNames <- which(table(taxon_dataPR$f) == 1) # monotipic families
names.fam <- setdiff(unique(taxon_dataPR$f), names(rm.famNames)) # removing monotipic families from the names 

for (i in 1:length(names.fam)) {
  set <- subset(taxon_dataPR, f == names.fam[i])
  phylo <- ape::makeNodeLabel(phylo, "u", nodeList = list(Fam_name = set$s))
  
  phylo$node.label[which(phylo$node.label == 
                           "Fam_name") ] <- paste(set$f[1])
}

pos.node <- unlist(lapply(names.fam, function(x){
  which(phylo$node.label == x) + length(phylo$tip.label)
}))

df.phylo <- data.frame(Fam.names = names.fam,
                       node.number = pos.node)

plot.base <- ggtree(phylo) + theme_tree2()
plot1 <- revts(plot.base) + scale_x_continuous(labels=abs)


PR.PG <- plot1 + geom_hilight(data = df.phylo, aes(node = node.number, fill = Fam.names), 
                      alpha = .6) +
  scale_fill_viridis(discrete = T, name = "Family names")
```

With this code we will obtain a phylogenetic tree like that in Figure 1. The species that are not shaded with colors corresponds to monotipic families.

![Fig. 1 - Phylogenetic tree containing 59 species of stream fishes in the Paraná and Paraguay Basins (Brazil).](../vignettes/Phylo_Parana-Paraguai.png)

---
title: "How to map species insertions in phylogenetic tree?"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{How to map species insertions in phylogenetic tree?}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

In this article we will show how we can map the insertions made in phylogenetic tree through the procedure realized in `FishPhyloMaker` function. Mapping the insertions is of special importance since allow to understand how the gaps in phylogenetic knowledge is distributed in the phylogenetic tree. For example, we can understand if the Characidae family present more phylogenetic gaps than Loricariidae family.

This can be done by using the information contained in the output of `FishPhyloMaker` function when we set the argument `return.insertions = TRUE`. This returns a data frame containing the categories of insertion for each species, these categories are:

-   **Present_in_tree** the species was already present in the original tree;
-   **Congeneric_insertion** species inserted as a sister species of the same genus presented in the tree;
-   **Congeneric_Family_level** species inserted as a sister species of the same genus presented in the tree, but that were added after a species of local pool of the same genus be inserted in the tree;
-   **Family_insertion** if not found any congeneric species, the species will be inserted near to, or between genus of the same family presented in the tree. The user can also insert the species in the base of the family;
-   **Order_insertion** if not found any genus of the same family of the species that must be inserted, the species will be inserted near to a given family, between two families or at the node that corresponds to the Order of this species;
-   **Not_inserted** if species was not inserted in any of the previous steps, it will not be inserted in the final tree;

To obtain this data frame we need first to run `FishPhyloMaker` setting the argument `return.insertions = TRUE`. We will use data of fish occurrence in the Neotropical region present in `{FishPhyloMaker}` package.

```{r setup}
library(FishPhyloMaker)
data("spp_afrotropic")
```

We need to format this data using function `FishTaxaMaker`

```{r formatData, echo=TRUE, eval=FALSE}
taxon_data <- FishTaxaMaker(data = spp_afrotropic, allow.manual.insert = TRUE)
```

With taxon_data we can run `FishPhyloMaker` to obtain the phylogeny and the data frame with all the insertions made by each species

```{r makingPhylo, echo=TRUE, eval=FALSE}
phylo_fish_Afrotropics <- FishPhyloMaker(data = taxon_data$Taxon_data_FishPhyloMaker, 
                                         return.insertions = TRUE,
                                         insert.base.node = TRUE, 
                                         progress.bar = TRUE)
```

The data frame can be extracted and the categories of insertion can be plotted in the phylogenetic tree by using the information on `phylo_fish_SAmerica$Insertions`

```{r familyNames, echo=TRUE, eval=FALSE}
library(phytools)
library(ggtree)
library(ggplot2)
insertions_org <- phylo_fish_Afrotropics$Insertions_data[match(tree$tip.label, phylo_fish_Afrotropics$Insertions_data$s), ]
p.base <- ggtree(tree, layout = "circular", size = .3)  %<+% insertions_org +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                 fontsize = 0) + #  plot the scale bar
  annotate("text", x = 4, y = 500, label = "20 myr", size = 1.5) # an attempt for add a scale bar

p.full <- p.base +
  geom_tippoint(aes(color = insertions), 
                size = .5, alpha = .8) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2)))  +
  scale_color_viridis_d(name = NULL, na.translate = F,
                        labels = c("Congeneric F", "Congeneric", "Family", "Order", "Present")) 
```

In Figure 1 we can see all the insertions made in the insertion process.

![Fig. 1 - Phylogenetic tree showing at which level each species was inserted in the tree.](../vignettes/phylo_afrotropics_insertions.png)% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_filter_rank.R
\name{filter_rank}
\alias{filter_rank}
\title{Function to download species from fishtreeoflife}
\usage{
filter_rank(order)
}
\arguments{
\item{order}{Character with names of species}
}
\value{
phylogeny
}
\description{
Function to download species from fishtreeoflife
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PD_deficit.R
\name{PD_deficit}
\alias{PD_deficit}
\title{Title Calculate the amount of phylogenetic deficit in assemblages}
\usage{
PD_deficit(phylo, data, level = "Congeneric_insertion")
}
\arguments{
\item{phylo}{Phylogenetic tree in newick format, can be an object from \code{\link{FishPhyloMaker}} function}

\item{data}{A data frame containing the classification  informing the level of insertions. This can be obtained
from \code{\link{FishPhyloMaker}} function}

\item{level}{Character indicating which level must be considered in the calculation of PD deficit.
Can be a vector with the levels ("Congeneric_insertion", "Congeneric_Family_level", "Family_insertion", "Order_insertion")
which will be considered in the calculation of phylogenetic deficit.
default is "Congeneric_insertion".}
}
\value{
A vector containing four values:
- Amount phylogenetic information present in the tree before insertions (PDintree)
- Amount of phylogenetic information inserted in the tree (PDdeficit)
- Total Phylogenetic information of the tree (PDtotal)
- A ratio calculated as PDdeficit/PDtotal (Darwinian_deficit)
}
\description{
Title Calculate the amount of phylogenetic deficit in assemblages
}
\seealso{
\code{\link{FishPhyloMaker}} for phylogeny and data frame containing the classification of insertions
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tab_function.R
\name{FishTaxaMaker}
\alias{FishTaxaMaker}
\title{Generate a list of species
Auxiliary function to obtain taxonomic classification and check the names of species
present in species pool}
\usage{
FishTaxaMaker(data, allow.manual.insert = TRUE)
}
\arguments{
\item{data}{A character vector with species names or a community matrix with species names in columns}

\item{allow.manual.insert}{Logical, if TRUE (default), the user must type the names of Family and Order of species
not found in Fishbase}
}
\value{
List with three elements.\preformatted{- A data frame containing the taxonomic classification of valid species accordingy to Fishbase

- A data frame with three columns containing the name of species (s), the Family (f) and Order (o) that can be used in
   FishPhyloMaker function
   
- A character vector containing all names of species that was not find in Fishbase
}
}
\description{
Generate a list of species
Auxiliary function to obtain taxonomic classification and check the names of species
present in species pool
}
\examples{
 \dontrun{
 data(neotropical_comm)
 data_comm <- neotropical_comm[, -c(1, 2)]
 taxon_data <- FishTaxaMaker(data_comm, allow.manual.insert = TRUE)
 Characidae
 Characiformes
 Characidae
 Characiformes
 Characidae
 Characiformes
 Loricariidae
 Siluriformes
 Characidae
 Characiformes
 Cichlidae
 Cichliformes
 Crenuchidae
 Characiformes
 Gymnotidae
 Gymnotiformes
 Loricariidae
 Siluriformes
 Loricariidae
 Siluriformes
 Loricariidae
 Siluriformes
 Loricariidae
 Siluriformes
 Heptapteridae
 Siluriformes
 Characidae
 Characiformes
 Loricariidae
 Siluriformes
 Characidae
 Characiformes
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_user_opt_printCatFamily.R
\name{print_cat_family}
\alias{print_cat_family}
\title{Internal function to help in interactive process}
\usage{
print_cat_family(print_cat, spp, order)
}
\arguments{
\item{print_cat}{Character}

\item{spp}{Character}

\item{order}{Character}
}
\value{
interactive action in console
}
\description{
Internal function to help in interactive process
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_description.R
\docType{data}
\name{taxon_data_PhyloMaker}
\alias{taxon_data_PhyloMaker}
\title{Data frame with species names needed to assemble the phylogenetic tree}
\format{
A data frame with taxonomic classification (species, family and order) of 45 species
}
\usage{
taxon_data_PhyloMaker
}
\description{
A data frame that can be directly used in FishPhyloMaker to obtain a phylogenetic tree
}
\references{
Species that make up the dataset in the paper published in Neotropical Ichthyology  \doi{10.1590/1982-0224-2020-0126}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_description.R
\docType{data}
\name{spp_afrotropic}
\alias{spp_afrotropic}
\title{List of fish species with occurrence in Afrotropical ecoregion}
\format{
A character vector with 767 species names:
}
\usage{
spp_afrotropic
}
\description{
A list of species that occur in basins of
Afrotropical ecoregion
}
\references{
\url{https://www.nature.com/articles/sdata2017141}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_treedata_modif.R
\name{treedata_modif}
\alias{treedata_modif}
\title{Internal function to show the number of species from species pool that lacks in phylogeny}
\usage{
treedata_modif(phy, data, sort = FALSE, warnings = TRUE)
}
\arguments{
\item{phy}{Phylogenetic hypothesis in newick format}

\item{data}{Data frame with species to be added in tree}

\item{sort}{Sorting species in alphabetic order}

\item{warnings}{Logical}
}
\value{
character vector
}
\description{
Internal function to show the number of species from species pool that lacks in phylogeny
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_user_opt_printCat.R
\name{print_cat}
\alias{print_cat}
\title{Internal function - auxiliary to interactive procedure}
\usage{
print_cat(print_cat, spp, family)
}
\arguments{
\item{print_cat}{Character}

\item{spp}{Character}

\item{family}{Character}
}
\value{
interactive action in console
}
\description{
Internal function - auxiliary to interactive procedure
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_description.R
\docType{data}
\name{neotropical_comm}
\alias{neotropical_comm}
\title{Abundance of stream fish species in Parana and Paraguay streams}
\format{
A data frame with 20 rows and 61 variables:
}
\source{
Article published in Neotropical Ichthyology  \doi{10.1590/1982-0224-2020-0126}
}
\usage{
neotropical_comm
}
\description{
A dataset containing the abundance of stream fish species distributed in streams of
Parana and Paraguay river Basins
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fishPhyloMaker.R
\name{FishPhyloMaker}
\alias{FishPhyloMaker}
\title{Obtaining fish phylogeny according to a local pool of species}
\usage{
FishPhyloMaker(
  data,
  insert.base.node = FALSE,
  return.insertions = TRUE,
  progress.bar = TRUE
)
}
\arguments{
\item{data}{A data frame with three columns containing the name of species (s), the Family (f) and the Order (o). This data frame can be generated
with tab_function function.}

\item{insert.base.node}{Logical argument indicating if the species must be added automatically
in the family and order (when needed) nodes. Default is FALSE}

\item{return.insertions}{Logical, if TRUE (default) the output is a list of length two containing the phylogeny
and a dataframe with a column indicating at which level each species was inserted.}

\item{progress.bar}{Logical argument. If TRUE (default) a progress bar will be shown in console.}
}
\value{
A newick object containing the phylogeny with the species in data object. If return.insertions = TRUE the output
will be a list of length two containing the newick phylogeny and a data frame equal that provided in data plus a column
indicating at which level each species was inserted in the tree.
}
\description{
Obtaining fish phylogeny according to a local pool of species
}
\examples{
\donttest{
    data("taxon_data_PhyloMaker")
    res_phylo <- FishPhyloMaker(data = taxon_data_PhyloMaker,
    insert.base.node = TRUE, 
    return.insertions = TRUE, 
    progress.bar = TRUE)
}
    

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/whichFishAdd.R
\name{whichFishAdd}
\alias{whichFishAdd}
\title{Function to inform which species must be added to the mega-tree phylogeny in the insertion process.}
\usage{
whichFishAdd(data)
}
\arguments{
\item{data}{A data frame with three column containing the name of species (s), the Family (f) and Order (o). This
can be generated with function \code{\link{FishTaxaMaker}}}
}
\value{
A data frame containing a column informing at which level the species in data must be added.
}
\description{
Function to inform which species must be added to the mega-tree phylogeny in the insertion process.
}
\details{
This function can be used  in order to known which species that must be added in the insertion process
made by \code{\link{FishPhyloMaker}}.
}
\examples{
\donttest{
    data("taxon_data_PhyloMaker")
    res_test <- whichFishAdd(data = taxon_data_PhyloMaker)
}


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_user_opt_printCat2.R
\name{print_cat2}
\alias{print_cat2}
\title{Internal function}
\usage{
print_cat2(spp)
}
\arguments{
\item{spp}{Character}
}
\value{
interactive action in console
}
\description{
Internal function
}
\keyword{internal}
