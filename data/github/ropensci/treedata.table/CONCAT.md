# treedata.table
##  An R package for manipulating phylogenetic data with _data.table_
A wrapper for data.table that enables fast manipulation of  phylogenetic trees matched to data. <img src='man/figures/logo.png' align="right" height="250" />

<!-- badges: start -->
[![codecov](https://codecov.io/gh/ropensci/treedata.table/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/treedata.table)
[![Devel version](https://img.shields.io/badge/devel%20version-0.1.0-blue.svg)](https://github.com/ropensci/treedata.table)
[![Lifecycle](https://img.shields.io/badge/lifecycle-ropensci/treedata.table-blue.svg)](https://www.tidyverse.org/lifecycle/#ropensci/treedata.table)
[![Code size](https://img.shields.io/github/languages/code-size/ropensci/treedata.table.svg)](https://github.com/ropensci/treedata.table)
[![Latest commit](https://img.shields.io/github/last-commit/ropensci/treedata.table.svg)](https://github.com/ropensci/treedata.table/commits/master)
[![rOpenSci badge](https://badges.ropensci.org/156_status.svg)](https://github.com/ropensci/onboarding/issues/367)
[![R-CMD-check](https://github.com/ropensci/treedata.table/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/treedata.table/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/ropensci/treedata.table/badge)](https://www.codefactor.io/repository/github/ropensci/treedata.table)
[![](https://cranlogs.r-pkg.org/badges/treedata.table)](https://cran.r-project.org/package=treedata.table)
<!-- badges: end -->

The [`data.table` package](https://github.com/Rdatatable/data.table) enables high-performance extended functionality for
data tables in R. `treedata.table` is a wrapper for `data.table` for phylogenetic analyses that matches a phylogeny to the
data.table, and preserves matching during `data.table` operations.


### tl;dr

`treedata.table` is a library for syncing data between a dataset and the tips of trees to enable efficient data and taxon management, as well as on-the-fly indexing and data selection in comparative analyses.


## Why use `treedata.table`?

Simultaneous processing of phylogenetic trees and data remains a computationally-intensive task. For example, processing a dataset of phylogenetic characters alongside a tree in `treedata.table` takes **90%** longer than processing the data alone in `data.table` (*Fig. 1A*). `treedata.table` provides new tools for increasing the speed and efficiency of phylogenetic data processing. Data manipulation in `treedata.table` is significantly faster than in other commonly used packages such as `base` (**>35%**), `treeplyr` (**>60%**), and `dplyr` (**>90%**). Additionally, `treedata.table` is **>400%** faster than `treeplyr` during the initial data/tree matching step (*Fig. 1B*).  

<div style="text-align:center">
<img src='man/figures/bench_TDT_Aug14.png' align="middle"width="600" />
</div>

 <font size="2"> **Fig. 1.** Results for the `treedata.table` microbenchmark during (**A**) data manipulation (`treedata.table[`) and (**B**) tree/data matching steps. We compare the performance of `treedata.table` against `treedata.table`, `base`, `treeplyr`, and `dplyr` using the `microbenchmark` R package.</font>


## Installing `treedata.table`

treedata.table can be installed from CRAN or GitHub at the present. We presently recommend installing using
[`remotes`](https://cran.r-project.org/web/packages/remotes/index.html) if using the GitHib version:
```r
install.packages("treedata.table")
remotes::install_github("ropensci/treedata.table")
library(treedata.table)
 ```


## What Can I Do With `treedata.table`?

`treedata.table` is designed with the intention of being able to efficiently manipulate trait data and
phylogenetic trees to enable comparative analyses. With the package is bundled some example data. Let's load it in and look at some common analyses.

```r
data(anolis)
td <- as.treedata.table(tree = anolis$phy, data = anolis$dat)
```

The function `as.treedata.table` converts a normal comma- or tab-delimited file to the `data.table` [format](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html). This enables a range of efficient and intuitive indexing and selection operations.

Now, we're going to use the `tdt` and `extractVector` functions in `data.table`. As an example, in our dataset is the column `SVL`, or snout-to-vent length in our anoles. We can index out this column on the fly, and run a Brownian motion analysis on this trait using the R package Geiger:

```r
tdt(td, geiger::fitContinuous(phy, extractVector(td, 'SVL'), model="BM", ncores=1))
```

We can also do efficient dropping of taxa from the analysis like so:

```r
dt <- droptreedata.table(tdObject=td, taxa=c("chamaeleonides" ,"eugenegrahami" ))
```

## Additional resources

More details about the functions implemented in `treedata.table` can be found in the different vignettes associated with the package or in our [website](https://ropensci.github.io/treedata.table/).

## Contributing

Please see our [contributing guide](CONTRIBUTING).

## Contact

Please see the package [DESCRIPTION](DESCRIPTION) for package authors.

## Citation

Please use the following citation when using the package:

*Josef Uyeda, Cristian Román-Palacios, and April Wright (2021). treedata.table: Manipulation of Matched Phylogenies and Data using 'data.table'. Available at: https://cran.r-project.org/web/packages/treedata.table/*


## License

Please see the package [LICENSE](LICENSE).


---

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Breve Introducción al Paquete treedata.table"
author: "Cristian Roman-Palacios, April Wright, Josef Uyeda"
date: "09/06/2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Breve Introducción al Paquete treedata.table}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Breve Introducción al Paquete `treedata.table`

El paquete `treedata.table` tiene como objetivo permitir a investigadores acceder y manipular datos filogenéticos usando herramientas del paquete `data.table`. `data.table` tiene diferentes funciones para manipular rápidamente datos en una forma eficiente.

El primer paso para usar `treedata.table` consiste en crear un objeto `treedata.table`. El objeto `treedata.table` empareja los tip.labels de la filogenia con una columna en el `data.frame` que contiene los nombres de los taxones. Este paso inicial permite la manipulación subsecuente y coordinada de los datos en el árbol y la matriz de caracteres dentro de `treedata.table`.

Para este tutorial vamos a usar los datos de _Anolis_ creados en `treeplyr`. Estos caracteres fueron generados aleatoriamente. Es importante resaltar tres aspectos. Primero, el árbol tiene que ser en formato `phylo` (o `multiPhylo` en caso de múltiples árboles). Segundo, la matriz de caracteres tiene que estár en formato `data.frame`. Tercero, la matriz de caracteres tiene que contener una columna con los nombres de los taxones coincidiendo con los tip.labels del árbol (o árboles). 

El objeto `treedata.table` se crea usando la función `as.treedata.table`.

```{r}
library(ape)
library(treedata.table)

# Cargamos los datos del ejemplo
data(anolis)
# Creamos el objecto treedata.table con as.treedata.table
td <- as.treedata.table(tree = anolis$phy, data = anolis$dat)
```

Podemos revisar el objeto resultante simplemente llamando el nombre del objeto en la consola. La matriz de caracteres, antes en `data.frame`, es ahora un `data.table`

```{r}
td
```

Adicionalmente, la matriz de caracteres en el nuevo formato `data.table` ha sido re-ordenado para tener las filas en el mismo orden que los tip.labels en el árbol.


```{r}
td$phy$tip.label == td$dat$tip.label
```


## Manipulando Datos

La matriz de datos en el objeto treedata.table puede ser indexada como cualquier otro objeto `data.table`. Por ejemplo, podemos hacer lo siguiente para extraer la columna con información de longitud hocico-cloaca (SVL) para cada especie.

```{r pressure, echo=TRUE}
td$dat[,'SVL']

```

También podemos usar los brackets dobles para extraer directamente la misma columna como un vector con nombres.

```{r}
td[["SVL"]]
```

El mismo resultado puede también ser logrado usando la función `extractVector`. Al igual que con los brackets dobles, el resultado de la función `extractVector` es un vector con nombres.


```{r}
extractVector(td, 'SVL')
```

Múltiples columnas pueden tambien ser extraídas usando `extractVector`.

```{r}
extractVector(td, 'SVL','ecomorph')
```

Hay un par de aspectos que son únicos a `[[.treedata.table()` y `extractVector()`. Primero, `[[.treedata.table()` tiene un argumento adicional que permite un correspondencia parcial del nombre de la columna (i.e. cuando el nombre objetivo tiene una superposicion parcial con los elementos en el objeto). Segundo, `extractVector()` puede extraer múltiples columnas y permite una evaluación no estandard (i.e. los nombres son tratados como cadenas de texto literal).


El poder real de `treedata.table` está en coindexar el árbol con la matriz de caracteres. Por ejemplo, en el siguiente comando usamos la sintaxis de `data.table` para extraer el primer representante de cada ecomorfo y retener todas las columnas. 

```{r}
 td[, head(.SD, 1), by = "ecomorph"]

```

Podemos hacer la misma operación con múltiples columnas.

```{r}
td[, head(.SD, 1), by = .(ecomorph, island)]

```

También implementamos la función `tail`.

```{r}
 td[, tail(.SD, 1), by = "ecomorph"]

```

Las columnas en `treedata.table` pueden ser operadas usando la misma sintaxis de `data.table`. En el siguiente ejemplo, los árboles solo incluirán especies distribuidas en Cuba. Este es el equivalente a filtrar usando `dplyr`. Después, una nueva columna llamada “Index” es creada en el objeto `data.table` dentro del objeto `treedata.table` con los valores de SVL+hostility. En resumen, la siguiente línea permite en forma simultánea crear una nueva columna y reducir el número de taxones en la filogenia a las especies de interés.

```{r}
td[island == "Cuba",.(Index=SVL+hostility)]
```

`treedata.table` permite aplicar funciones directamente en nuestros datos de interés. En el siguiente ejemplo, evaluamos un modelo de evolución browniano sobre los datos de SVL en nuestro set de datos. Usamos una combinación de `tdt`, `extractVector` y `geiger::fitContinuous` para correr funciones en nuestros datos, extraer un vector de caracteres y ajustar el model en cuestión, respectivamente.

```{r}

tdt(td, geiger::fitContinuous(phy, extractVector(td, 'SVL'), model="BM", ncores=1))

```

Los terminales en el árbol también pueden ser removidos fácilmente, con los cambios también reflejados sobre la matriz de caracteres. En el siguiente ejemplo, removemos dos taxones por sus nombres.


```{r}
dt <- droptreedata.table(tdObject=td, taxa=c("chamaeleonides" ,"eugenegrahami" ))
```

Revisamos si *A. chamaeleonides* y *A. eugenegrahami* aún están en el árbol.

```{r}
c("chamaeleonides" ,"eugenegrahami" ) %in% dt$phy$tip.label
```

Y podemos hacer lo mismo con la matriz de caracteres en el nuevo objeto `treedata.table`.

```{r}
c("chamaeleonides" ,"eugenegrahami" ) %in% dt$dat$X
```

Por último, el árbol y la matriz de caracteres pueden ser extraídos de el objeto `treedata.table` fácilmente usando la función ` pulltreedata.table`. 

```{r}
df <- pulltreedata.table(td, "dat")
tree <- pulltreedata.table(td, "phy")
```

La tabla

```{r}
df
```

Y el árbol

```{r}
tree
```

La misma funcionalidad explicada en este tutorial sobre objetos `phylo` aplica directamente a objetos `multiPhylo`.
---
title: "Additional Functions for Manipulating Data in treedata.table"
author: "Josef Uyeda, Cristian Roman-Palacios, April Wright"
date: "08/08/2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Additional Functions for Manipulating Data in treedata.table}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Additional functions for manipulating data

`treedata.table` includes additional functions that allow the identification of `discrete` and `continuous` characters in a given dataset. We first load the dataset:

```{r}
library(ape)
library(treedata.table)

# Load example data
data(anolis)
#Create treedata.table object with as.treedata.table
td <- as.treedata.table(tree = anolis$phy, data = anolis$dat)
```

The `detectCharacterType()` function can be used to examine whether `SVL` is `discrete` or `continuous`: 

```{r}
detectCharacterType(anolis$dat$SVL)
```

We can further examine the type of characters we have in our dataset by using the `detectAllCharacters()` function:

```{r}
detectAllCharacters(anolis$dat)
```

Summarizing this information in a table, we get:

```{r}
cbind.data.frame(character=colnames(anolis$dat),type=detectAllCharacters(anolis$dat))
```

Finally, we can use the `filterMatrix()` function to subset our dataset to only a certain type of characters. For instance, let's extract all discrete characters in the *Anolis* dataset:

```{r}
filterMatrix(anolis$dat, "discrete")
```

Two additional functions in `treedata.table` are designed to examine and modify column and row names in any dataset. For instance, we can ask if the *Anolis* dataset has column names:

```{r}
hasNames(anolis$dat, "col")
```

It does have column names. Let's just remove the column names and check if `hasNames()` can detect this change. Here's our new dataset:

```{r}
data=anolis$dat
colnames(data)<-NULL
head(data,2)
```

Let's run `hasNames()` on our new dataset:

```{r}
hasNames(data, "col")
```

Now, we can create new column names using the `forceNames()` function:

```{r}
data <- forceNames(data, "col")
```

The new dataset, with column names (n1...), looks like this:

```{r}
head(data,2)
```

And we can finally ask whether our new dataset actually have column names by running the `hasNames()` function again:


```{r}
hasNames(data, "col")
```


We can apply the same procedure for columns (`col`), rows (`row`) or both (`rowcol`).

---
title: "Partially Matching of Trait Data and Tree(s) in treedata.table"
author: "Josef Uyeda, Cristian Roman-Palacios, April Wright"
date: "08/08/2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Partially Matching of Trait Data and Tree(s) in treedata.table}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Partially matching trait data and tree(s)

The `as.treedata.table` function enables users to match a tree (or multiple trees) against a single trait database. We first load the sample dataset.

```{r}
library(ape)
library(treedata.table)

# Load example data
data(anolis)
#Create treedata.table object with as.treedata.table
td <- as.treedata.table(tree = anolis$phy, data = anolis$dat)
```


Tips that are not common between the tree (or trees) and dataset are dropped from the resulting `treedata.table` object. For instance, below I have modified the original anole phylogeny such that *A. ahli* (**ahli**) is replaced for a label that is not present in the dataset (**NAA**).

```{r}
anolis_newtip<-anolis$phy
anolis_newtip$tip.label[1]<-'NAA'
anolis_newtip
```


We then use this modified tree to fit a `treedata.table` object using the `as.treedata.table` function:


```{r}
td <- as.treedata.table(tree=anolis_newtip, data=anolis$dat)
```

Note that `as.treedata.table` drops all non-overlapping tips (**NAA** [present in the tree but not in the trait data] and **ahi** [present in the database but not in tree] in this case) and returns a `treedata.table` object with fully matching `phy` and `data` objects. 

```{r}
td
```


Fully-matching matrix and trees are also returned in `treedata.table` objects with `multiPhylo` objects in their `phy` component. See the example below.

We first construct a `multiPhylo` object that partially overlaps the original trait database by using **NAA** instead of **ahi**.

```{r}
anolis2<-anolis$phy
anolis2$tip.label[1]<-'NAA'
anolis1<-anolis$phy
anolis1$tip.label[1]<-'NAA'
trees<-list(anolis1,anolis2)
class(trees) <- "multiPhylo"
trees
```

Next, we fit the `treedata.table` object using the relevant `multiPhylo` object and the original trait database.

```{r}
td <- as.treedata.table(tree=trees, data=anolis$dat)
```

Note that 1 tip was dropped for all trees in the `multiPhylo` object and a single row was deleted from the `data.table` object in the `treedata.table` object.

```{r}
td
```

---
title: "Getting Started With The treedata.table Package"
author: "Josef Uyeda, Cristian Roman-Palacios, April Wright"
date: "08/08/2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting Started With The treedata.table Package}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Getting Started With The `treedata.table` Package

The aim of the `treedata.table` R package is to allow researchers to access and manipulate phylogenetic data using tools from the `data.table` package. `data.table` has many functions for rapidly manipulating data in a memory efficient way. 

Using the `treedata.table` package begins with creating a `treedata.table` object. The `treedata.table` matches the tip.labels of the phylogeny to a column of names in your `data.frame`. This allows you to manipulate the data, and the corresponding tree together. 

Importantly, the character matrix must must include a column with the taxa names and should be of class `data.frame`. The tree must be of class `phylo` or `multiPhylo`.

A `treedata.table` is created using the `as.treedata.table` function. Here we use the 
*Anolis* dataset from `treeplyr`. Traits in this dataset were randomly generated for a set of 100 species.

```{r}
library(ape)
library(treedata.table)

# Load example data
data(anolis)
#Create treedata.table object with as.treedata.table
td <- as.treedata.table(tree = anolis$phy, data = anolis$dat)
```

We may inspect our object by calling it by name. You will notice that your `data.frame` is now a `data.table`. A `data.table` is simply an advanced version of a `data.frame` that, among other, increase speed in data manipulation steps while simplifying syntax. 

```{r}
td
```


Furthermore, the new `data.table` has been reordered into the same order as the tip.labels of your tree.

```{r}
td$phy$tip.label == td$dat$tip.label
```


## Manipulating Data


### Coindexing

Your data table can be indexed in the same way any other `data.table` object would be. For example, if we wanted to look at our snout-vent length column, we can do that like so.

```{r pressure, echo=TRUE}
td$dat[,'SVL']

```

You can also use double bracket syntax to directly return column data as a named list.

```{r}
td[["SVL"]]
```

The same functionality can also be accomplished through the `extractVector` function. Both the double bracket syntax and the `extractVector` function will return a named vector.

```{r}
extractVector(td, 'SVL')
```

Multiple traits can also be extracted using `extractVector`.

```{r}
extractVector(td, 'SVL','ecomorph')
```

However, there's a couple aspects that are unique to `[[.treedata.table()` and `extractVector()`. First, `[[.treedata.table()` has an extra exact argument to enable partial match (i.e. when target strings and those in the `treedata.table` object match partially). Second, `extractVector()` can extract multiple columns and accepts non-standard evaluation (i.e. names are treated as string literals).

The real power in treedata.table is in co-indexing the tree and table. For example, in the below command, we use `data.table` syntax to take the first representative from each ecomorph. We retain all data columns. If you examine the tree object, you will see that it has had all the tips absent from the resultant `data.table`.

```{r}
 td[, head(.SD, 1), by = "ecomorph"]

```

We could also do the same operation with multiple columns:

```{r}
td[, head(.SD, 1), by = .(ecomorph, island)]

```

Tail is also implemented

```{r}
 td[, tail(.SD, 1), by = "ecomorph"]

```

Columns in the `treedata.table` object can also be operated on using general `data.table` syntax. In the below example, the tree is pruned to those tips that occur in Cuba. This is the `data.table` equivalent of `dplyr`'s filter. Then, a new column is created in the `data.table`, assigned the name "Index", and assigned the value of the SVL + the hostility index. This enables concurrent manipulation of the phylogeny, and the calculation of a new index for only those tips we would actually like to use.

```{r}
td[island == "Cuba",.(Index=SVL+hostility)]
```

### Running functions on `treedata.table` objects

In the below command, we extract one vector from our data.table and use `geiger`'s continuous model fitting to estimate a Brownian motion model for the data using the `tdt` function. 


```{r}

tdt(td, geiger::fitContinuous(phy, extractVector(td, 'SVL'), model="BM", ncores=1))

```

### Dropping and extracting taxa from `treedata.table` objects

We can also drop tips directly from the tree, and have those tips drop concurrently from the data.table. In the example below, we  remove two taxa by name. 

```{r}
dt <- droptreedata.table(tdObject=td, taxa=c("chamaeleonides" ,"eugenegrahami" ))
```

We can check if *A. chamaeleonides* and *A. eugenegrahami* are still in the tree

```{r}
c("chamaeleonides" ,"eugenegrahami" ) %in% dt$phy$tip.label
```

And we can do the same with the data in the `treedata.table` object

```{r}
c("chamaeleonides" ,"eugenegrahami" ) %in% dt$dat$X
```


When you're done, the data.table and tree can both be extracted from the object:

```{r}
df <- pulltreedata.table(td, "dat")
tree <- pulltreedata.table(td, "phy")
```

The table

```{r}
df
```

and the corresponding tree

```{r}
tree
```








---
title: "Working With multiPhylo Objects in treedata.table"
author: "Josef Uyeda, Cristian Roman-Palacios, April Wright"
date: "08/08/2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Working With multiPhylo Objects in treedata.table}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Working with multiphylo objects

`treedata.table` further allows the matching of multiple phylogenies (`multiPhylo`) against a single dataset (`data.frame`). Below, we modified the anole dataset to explain the extended functionality of `treedata.table` with `multiPhylo` objects. Note that all the trees in the `multiPhylo` must have exactly the same taxa.

We first load the sample dataset.

```{r}
library(ape)
library(treedata.table)

# Load example data
data(anolis)
#Create treedata.table object with as.treedata.table
td <- as.treedata.table(tree = anolis$phy, data = anolis$dat)
```


We then create a `multiPhylo` object including only two `phylo` objects. Users can provide any number of `phylo` objects within the `multiPhylo` object. However, trees can only differ in their topology. In other words, all trees must have the same tip labels. 

We also note that both the provided `multiPhylo` and `data.frame` should partially overlap

```{r}
trees<-list(anolis$phy,anolis$phy)
class(trees) <- "multiPhylo"
trees
```

Now, we create our treedata.table object by combining the trait data (`data.frame`) and the newly generated `multiPhylo` object. Note that there is only a single character matrix.

```{r}
td <- as.treedata.table(tree=trees, data=anolis$dat)
```

The resulting `td` object now returns a `multiPhylo` object under `phy`. This objectcontains only the overlapping taxa between the multiphylo objects and the input dataset.

```{r}
class(td$phy);td$phy
```

Please note that all the basic `treedata.table` functions highlighted above for `phylo` objects are still functional when `treedata.table` objects include `multiPhylo` objects.

```{r}
td[, head(.SD, 1), by = "ecomorph"]
```

Functions can also be run on any `treedata.table` object with `multiphylo` data. For instance, the following line will fit a phenogram for `SVL` on each of the trees we provided in the `multiPhylo` object.

```{r}
tdt(td, geiger::fitContinuous(phy, extractVector(td, 'SVL'), model="BM", ncores=1))
```

The output is an object of class `list` with each element corresponding to the output function of each tree in the provided `multiPhylo` object. 



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treedata.table.R
\name{[[.treedata.table}
\alias{[[.treedata.table}
\title{Function for extract a named vector from an object of class \code{treedata.table}}
\usage{
\method{[[}{treedata.table}(x, ..., exact = TRUE)
}
\arguments{
\item{x}{An object of class \code{treedata.table}}

\item{...}{Column name in class \code{character}}

\item{exact}{whether exact search should be conducted}
}
\value{
A new object of class \code{vector} with names set to labels corresponding
to tip labels in the provided \code{treedata.table} object.
}
\description{
This function extracts a named vector for any  trait from an object of class
\code{treedata.table}.
}
\examples{
data(anolis)
# With a phylo object
td <- as.treedata.table(anolis$phy, anolis$dat)
td[["SVL"]]

# With a multiPhylo object
treesFM <- list(anolis$phy, anolis$phy)
class(treesFM) <- "multiPhylo"
td <- as.treedata.table(treesFM, anolis$dat)
td[["SVL"]]
}
\seealso{
\code{\link[=data.table]{data.table()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractVector.R
\name{extractVector}
\alias{extractVector}
\title{Returning a named vector from a treedata.table object}
\usage{
extractVector(tdObject, ...)
}
\arguments{
\item{tdObject}{A treedata.table object}

\item{...}{The name of the column or columns to select.}
}
\value{
A named vector or a list of named vectors
}
\description{
Returning a named vector from a treedata.table object
}
\examples{

data(anolis)
td <- as.treedata.table(tree = anolis$phy, data = anolis$dat)
# extracts the named vector for SVL from the td object
extractVector(td, "SVL")
# extracts the named vector for SVL and ecomorph from the td object
extractVector(td, "SVL", "ecomorph")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anolis.R
\docType{data}
\name{anolis}
\alias{anolis}
\title{Anole data}
\format{
An object of class \code{list} of length 2.
}
\usage{
data(anolis)
}
\description{
Anole data for treedata.table functions. Many
of the traits (e.g. awesomeness, hostility) in this dataset
retrieved from treeplyr are based on random numbers.
}
\author{
Luke Harmon
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treedata.table.R
\name{[.treedata.table}
\alias{[.treedata.table}
\title{Function for performing data.table operations on an object of class
\code{treedata.table}}
\usage{
\method{[}{treedata.table}(x, ...)
}
\arguments{
\item{x}{An object of class \code{treedata.table}}

\item{...}{Arguments in the structure of \code{data.table} used to perform changes
on the \code{treedata.table} object}
}
\value{
A new object of class \code{treedata.table} with \verb{$dat} and \verb{$phy}
corresponding with the changes set to \verb{$dat} using
\link[data.table:data.table]{data.table}'s structure.
}
\description{
This function can be used to subset rows, select and compute on columns
\link[data.table:data.table]{data.table}.
}
\examples{

data(anolis)
anolis2 <- anolis$phy
anolis2$tip.label[1] <- "NAA"
anolis1 <- anolis$phy
anolis1$tip.label[1] <- "NAA"
trees <- list(anolis1, anolis2)
class(trees) <- "multiPhylo"
treesFM <- list(anolis$phy, anolis$phy)
class(treesFM) <- "multiPhylo"

# A phylo object that fully matches the data
td <- as.treedata.table(tree = anolis$phy, data = anolis$dat)
td <- as.treedata.table(anolis$phy, anolis$dat)
td[, SVL]
td[island == "Cuba" & ecomorph == "TG", .(ecomorph, island, SVL)]
td[, utils::head(.SD, 1), by = .(ecomorph, island)]

# A multiphylo object that fully matches the data
td <- as.treedata.table(tree = treesFM, data = anolis$dat)
td <- as.treedata.table(treesFM, anolis$dat)
td[, SVL]
td[island == "Cuba" & ecomorph == "TG", .(ecomorph, island, SVL)]
td[, utils::head(.SD, 1), by = .(ecomorph, island)]

# A phylo object that partially matches the data
td <- as.treedata.table(tree = anolis1, data = anolis$dat)
td <- as.treedata.table(anolis1, anolis$dat)
td[, SVL]
td[island == "Cuba" & ecomorph == "TG", .(ecomorph, island, SVL)]
td[, utils::head(.SD, 1), by = .(ecomorph, island)]

# A multiphylo object that partially matches the data
td <- as.treedata.table(tree = trees, data = anolis$dat)
td <- as.treedata.table(trees, anolis$dat)
td[, SVL]
td[island == "Cuba" & ecomorph == "TG", .(ecomorph, island, SVL)]
td[, utils::head(.SD, 1), by = .(ecomorph, island)]
}
\seealso{
\link[data.table:data.table]{data.table}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/droptreedata.table.R
\name{droptreedata.table}
\alias{droptreedata.table}
\title{Function dropping taxa from an object of class \code{treedata.table}}
\usage{
droptreedata.table(tdObject, taxa)
}
\arguments{
\item{tdObject}{An object of class \code{treedata.table}}

\item{taxa}{A vector class \code{character} containing all taxa that needs to be
dropped from the original \code{treedata.table} object}
}
\value{
An object of class \code{treedata.table} with matching \verb{$dat} and \verb{$phy}
elements based on whether \code{taxa} were dropped or not.
}
\description{
This function can be used to remove species from an object of class
\code{treedata.table}. The resulting \code{treedata.table} will include fully matching
\verb{$dat} and \verb{$phy} elements. The user should confirm the changes before they
are processed.
}
\examples{
data(anolis)
# With a multiphylo object in the treedata.table object
td <- as.treedata.table(anolis$phy, anolis$dat)
droptreedata.table(
  tdObject = td, taxa =
    c("chamaeleonides", "eugenegrahami")
)

# With a multiphylo object in the treedata.table object
treesFM <- list(anolis$phy, anolis$phy)
class(treesFM) <- "multiPhylo"
td <- as.treedata.table(treesFM, anolis$dat)
droptreedata.table(
  tdObject = td, taxa =
    c("chamaeleonides", "eugenegrahami")
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdt_methods.R
\name{summary.treedata.table}
\alias{summary.treedata.table}
\title{Summarizing treedata.table objects}
\usage{
\method{summary}{treedata.table}(object, ...)
}
\arguments{
\item{object}{an object of class "treedata.table"}

\item{...}{additional arguments passed to "head.treedata.table"}
}
\value{
Function tries to be smart about summarizing the data
and detecting continuous vs. discrete data, as well as whether
any data have missing data. Also returns the taxa that are
dropped from either the original tree or the original data.
}
\description{
Summarizing treedata.table objects
}
\examples{
data(anolis)
td <- as.treedata.table(anolis$phy, anolis$dat)
summary(td)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectCharacterType.R
\name{detectCharacterType}
\alias{detectCharacterType}
\title{Function to detect whether a character is continuous or discrete}
\usage{
detectCharacterType(dat, cutoff = 0.1)
}
\arguments{
\item{dat}{A vector of data}

\item{cutoff}{Cutoff value for deciding if numeric data might actually be
discrete: if nlev is the number of levels and n the length of dat, then
nlev / n should exceed cutoff, or the data will be classified as discrete}
}
\value{
Either "discrete" or "continuous"
}
\description{
This function detects whether a given vector is a continuous
(e.g., with values 2.45, 9.35, and so on) or a discrete
(e.g., with values blue, red, yellow) character.
}
\examples{
data(anolis)
detectCharacterType(anolis$dat[, 1])
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectCharacterType.R
\name{filterMatrix}
\alias{filterMatrix}
\title{Filter a character matrix, returning either all continuous or all discrete characters}
\usage{
filterMatrix(mat, returnType = "discrete")
}
\arguments{
\item{mat}{A character matrix of class data.frame}

\item{returnType}{Either discrete or continuous}
}
\value{
data.frame with only discrete (default) or continuous characters
}
\description{
This function filters a character matrix based on continuous
(e.g., with values 2.45, 9.35, and so on) or discrete characters
(e.g., with values blue, red, yellow).
}
\examples{
data(anolis)
filterMatrix(anolis$dat, "discrete")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectCharacterType.R
\name{forceNames}
\alias{forceNames}
\title{Force names for rows, columns, or both}
\usage{
forceNames(dat, nameType = "row")
}
\arguments{
\item{dat}{A vector of data}

\item{nameType, }{either:
\describe{
\item{"row"}{Rows (default)}
\item{"col"}{Columns}
\item{"rowcol"}{Both rows and columns}
}}
}
\value{
An object of type `data.frame with labeled columns, rows, or both.
}
\description{
This function creates column names (\code{colnames}), row.names (\code{row.names}),
or both in an unnamed \code{data.frame} or \code{matrix}.
}
\examples{
data(anolis)
forceNames(anolis$dat, "row")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdt_methods.R
\name{tail.treedata.table}
\alias{tail.treedata.table}
\title{Return the last part of an treedata.table object}
\usage{
\method{tail}{treedata.table}(x, ...)
}
\arguments{
\item{x}{a treedata.table object}

\item{...}{Additional arguments passed to head.data.table}
}
\description{
Return the last part of an treedata.table object
}
\examples{
data(anolis)
td <- as.treedata.table(anolis$phy, anolis$dat)
tail(td)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.treedata.table.R
\name{as.treedata.table}
\alias{as.treedata.table}
\title{Combine tree (or set of trees) and data.frame into a single treedata.table
object}
\usage{
as.treedata.table(tree, data, name_column = "detect")
}
\arguments{
\item{tree}{A tree of class \code{phylo} or multiple trees of class \code{multiPhylo}}

\item{data}{A dataset in format \code{data.frame}}

\item{name_column}{A character indicating the name of taxa in \code{data.frame}.
If set to \code{detect} (default) \verb{as treedata.table} will auto-detect this
column}
}
\value{
An object of type \code{treedata.table} containing the tree and data.table
}
\description{
This function takes as input a tree of class \code{phylo} or \code{multiPhylo} and a
\code{data.frame} and combines them into a treedata.table. If a \code{multiPhylo} is
provided, all trees must have the same tip.labels. \code{treedata.table} object is
sorted such that the rows in the data.table are matched to the tip.labels
of the phylogeny. Tip.labels on the tree must match a column of tip
names in the input data.frame. The output of this function will be a
treedata.table, which can be manipulated as a data.table.
}
\examples{

data(anolis)
anolis2 <- anolis$phy
anolis2$tip.label[1] <- "NAA"
anolis1 <- anolis$phy
anolis1$tip.label[1] <- "NAA"
trees <- list(anolis1, anolis2)
class(trees) <- "multiPhylo"
treesFM <- list(anolis$phy, anolis$phy)
class(treesFM) <- "multiPhylo"

# A phylo object that fully matches the data
td <- as.treedata.table(tree = anolis$phy, data = anolis$dat)
# A multiphylo object that fully matches the data
td <- as.treedata.table(tree = treesFM, data = anolis$dat)
# A phylo object that partially matches the data
td <- as.treedata.table(tree = anolis1, data = anolis$dat)
# A multiphylo object that partially matches the data
td <- as.treedata.table(tree = trees, data = anolis$dat)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdt_methods.R
\name{print.treedata.table}
\alias{print.treedata.table}
\title{Print method treedata.table objects}
\usage{
\method{print}{treedata.table}(x, ...)
}
\arguments{
\item{x}{an object of class "treedata.table"}

\item{...}{additional arguments passed to "head.treedata.table"}
}
\value{
Function uses prints the tree and the first lines of the
data.table object.
}
\description{
Print method treedata.table objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulltreedata.table.R
\name{pulltreedata.table}
\alias{pulltreedata.table}
\title{Returns the character matrix or phylogeny from a treedata.table object}
\usage{
pulltreedata.table(tdObject, type = c("dat", "phy"))
}
\arguments{
\item{tdObject}{A treedata.table object}

\item{type}{Whether the output of the function is a tree ('type="phylo"')
or a data.table ('type="dat"')}
}
\value{
A \code{data.table} or \code{phylo} object from the original \code{treedata.table}
object
}
\description{
Returns the character matrix or phylogeny from a treedata.table object
}
\examples{
data(anolis)
td <- as.treedata.table(anolis$phy, anolis$dat)
pulltreedata.table(td, type = "phy")
pulltreedata.table(td, type = "dat")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdt_methods.R
\name{head.treedata.table}
\alias{head.treedata.table}
\title{Return the first part of an treedata.table object}
\usage{
\method{head}{treedata.table}(x, ...)
}
\arguments{
\item{x}{a treedata.table object}

\item{...}{Additional arguments passed to head.data.table}
}
\description{
Return the first part of an treedata.table object
}
\examples{
data(anolis)
td <- as.treedata.table(anolis$phy, anolis$dat)
head(td)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectCharacterType.R
\name{detectAllCharacters}
\alias{detectAllCharacters}
\title{Apply detectCharacterType over an entire matrix}
\usage{
detectAllCharacters(mat, cutoff = 0.1)
}
\arguments{
\item{mat}{A matrix of data}

\item{cutoff}{Cutoff value for deciding if numeric data might actually be
discrete: if nlev is the number of levels and n the length of dat, then
nlev / n should exceed cutoff, or the data will be classified as discrete}
}
\value{
Vector of either "discrete" or "continuous" for each variable in
matrix
}
\description{
This function detects whether each column in a matrix is a continuous
(e.g., with values 2.45, 9.35, and so on) or a discrete character
(e.g., with values blue, red, yellow).
}
\examples{
data(anolis)
detectAllCharacters(anolis$dat)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectCharacterType.R
\name{hasNames}
\alias{hasNames}
\title{Row and column name check}
\usage{
hasNames(dat, nameType = "row")
}
\arguments{
\item{dat}{A vector of data}

\item{nameType, }{either:
\describe{
\item{"row"}{Rows (default)}
\item{"col"}{Columns}
\item{"rowcol"}{Both rows and columns}
}}
}
\value{
\code{TRUE} or \code{FALSE} indicating if the object has names (\code{columns},
\code{rows}, or
\code{both})
}
\description{
This function checks whether a given \code{data.frame} or \code{matrix} has
column names (\code{colnames}), row.names (\code{row.names}), or both.
}
\examples{
data(anolis)
hasNames(anolis$dat, "row")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdt.R
\name{tdt}
\alias{tdt}
\title{Run a function on a \code{treedata.table} object}
\usage{
tdt(tdObject, ...)
}
\arguments{
\item{tdObject}{A treedata.table object}

\item{...}{A function call.}
}
\value{
Function output for a single tree (phylo) or a list of function
outputs (one per each tree in the MultiPhylo object)
}
\description{
Run a function on a \code{treedata.table} object
}
\details{
This function allows R functions that use trees and data to be run
on\code{treedata.table} objects.
}
\examples{
data(anolis)
\donttest{

# A treedata.table object with a phylo $phy
td <- as.treedata.table(anolis$phy, anolis$dat)
tdt(td, geiger::fitContinuous(phy, extractVector(td, "SVL"),
  model = "BM", ncores = 1
))


# A treedata.table object with a multiPhylo $phy
treesFM <- list(anolis$phy, anolis$phy)
class(treesFM) <- "multiPhylo"
td <- as.treedata.table(treesFM, anolis$dat)
tdt(td, geiger::fitContinuous(phy, extractVector(td, "SVL"),
  model = "BM",
  ncores = 1
))
}
}
