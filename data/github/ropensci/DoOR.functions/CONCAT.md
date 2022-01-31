DoOR.functions
==============
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.375617.svg)](http://dx.doi.org/10.5281/zenodo.375617)
[![](https://badges.ropensci.org/34_status.svg)](https://github.com/ropensci/onboarding/issues/34)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/DoOR.functions.svg?branch=master)](https://travis-ci.org/ropensci/DoOR.functions)

R package containing the functions used to build the Database of Odor Responses. The corresponding data package can be found at [https://github.com/ropensci/DoOR.data](https://github.com/ropensci/DoOR.data).

## The DoOR Project
Find more information, precompiled R-packages and an interactive web-version of the DoOR-Database at: **[http://neuro.uni.kn/DoOR](http://neuro.uni.kn/DoOR)**

## Install
Either download a packaged version or install _via_ `devtools`:
```{r}
# install devtools
install.packages("devtools")
library(devtools)

# install DoOR.functions 2.0.1
install_github("ropensci/DoOR.functions", ref="v2.0.1")

# or install the latest version available on Github
install_github("ropensci/DoOR.functions")
```

## Publications
DoOR was first published in 2010, the OpenAccess publication is available from
[http://chemse.oxfordjournals.org/content/35/7/551](http://chemse.oxfordjournals.org/content/35/7/551 "10.1093/chemse/bjq042")

An OpenAccess publication regarding the comprehensive update to **DoOR version 2.0** is available from
A preprint of manuscript related to DoOR 2.0 can be found at bioRxiv: [http://www.nature.com/articles/srep21841](http://www.nature.com/articles/srep21841)

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# DoOR.functions v2.0.1
- added unit testing

# DoOR.functions v2.0.1
- converted function names from camelCase to snake_case
- inform about new function names upon package loading
- added warning message that data is loaded into current workspace and forced user feedback (#17)


# DoOR.functions v2.0.0
- A comprehensive update to data and functions of the DoOR project. Please see the publication for details: http://doi.org/10.1038/srep21841
- rewrote many of DoOR's core functions to improve speed and readability
- updated all documentation
- introduced InChIKeys as main chemical identifier
- added several new functions for analysis and plotting
- added vignettes

## Bugfixes

- updateDatabase\(\) calculates with wrong permutation (#2)
- calModel\\(\\) fails for studies containing only 0s (#1) 


# DoOR.functions v1.0.2
- several bugfixes


# DoOR.functions v1.0
- initial release as published in: Integrating heterogeneous odor response data into a common response model: A DoOR to the complete olfactome. ChemSenses 35, 551–63. http://doi.org/10.1093/chemse/bjq042---
title: "The Database of Odor Responses - DoOR functions package"
author: "Daniel Münch"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{The Database of Odor Responses - DoOR functions package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

DoOR consists of two R packages and both are needed for DoOR to work properly.
One package, `DoOR.data` contains all the _Drosophila_ odor responses we
gathered from labs around the world or recorded ourselves. The other package
`DoOR.functions` contains the DoOR framework for integrating heterogeneous data
sets as well as analysis tools and functions for visualization.

In this vignette we describe how to build, modify and update DoOR and introduce
some helper functions. There are two other vignettes explaining the [plotting
functions](DoOR_visualizations.html) and the [analysis tools](DoOR_tools.html)
in detail.
##

## Content
* [loading DoOR](#loading)
* [Modifying, building and updating DoOR](#building)
    * [Importing new data with `import_new_data()`](#import_new_data)
    * [Building the complete data base with `create_door_database()`](
    #create_door_database)
    * [Updating parts of the data base with `update_door_database()`](
    #update_door_database)
    * [`model_response()` and `model_response_seq()`](#model)
    * [Removing a study with `remove_study()`](#remove_study)
    * [Updating the odor information with `update_door_odorinfo()`](
    #update_door_odorinfo)
* [Helper functions](#helper)


## Loading DoOR{#loading}
The first step after starting R is to attach both packages and to load the
response data:

```{r, results='hide'}
library(DoOR.data)
library(DoOR.functions)
load_door_data(nointeraction = TRUE)
```
`load_door_data()` attaches the data from `DoOR.data`.



# Modifying, building and updating DoOR{#building}
DoOR comes with all the original data sets as well as with a pre-computed
version of the consensus matrix `door_response_matrix` where all data was
integrated using the DoOR merging algorithms (see paper for details on how the
algorithm works). The values in `door_response_matrix` are globally normalized
with values scaled `[0,1]`. `door_response_matrix_non_normalized` is a version
of the consensus data that is not globally normalized meaning that responses are
scaled `[0,1]` within each _responding unit_ (receptor, sensory neuron,
glomerulus...).


## Importing new data with `import_new_data()`{#import_new_data}
It is easy to add new response data to DoOR, we only have to take care to
provide it in the right format:

* either a .csv or a .txt file with fields separated by colons or tabs (see
`?read.table` for detailed specifications). * the filename corresponds to the
later name of the data set * if we add e.g. recordings obtained with different
methods, these should go into two data sets and thus into two different files
that we import * e.g. "Hallem.2004.EN" and "Hallem.2004.WT" are the "empty
neuron" and the "wildtype neuron" recordings from Elissa Hallem's 2004
publication * the file needs at least two columns: 1. one column named
"InChIKey" holding the InChIKey of the odorant 1. one column named after the
responding unit the recording comes from (e.g. "Or22a")

A minimal example file could look like this:
```{r, echo=FALSE, }
tmp <- Or22a[c(1,3:5), c(3,6)]
colnames(tmp)[2] <- "Or22a"
knitr::kable(tmp)
```

We can provide more chemical identifiers:
```{r, echo=FALSE, }
tmp <- Or22a[c(1,3:5), c(1:6)]
colnames(tmp)[6] <- "Or22a"
knitr::kable(tmp)
```

Any of the following will be imported:

**`Class`**
: e.g. "ester"
: the chemical class an odorant belongs to

**`Name`**
: e.g. "isopentyl acetate"

**`InChIKey`**
: e.g. "MLFHJEHSLIIPHL-UHFFFAOYSA-N" 
([details](https://en.wikipedia.org/wiki/International_Chemical_Identifier))

**`InChI`**
: e.g. "InChI=1S/C7H14O2/c1-6(2)4-5-9-7(3)8/h6H,4-5H2,1-3H3" ([details](https://en.wikipedia.org/wiki/International_Chemical_Identifier))

**`CAS`**
: e.g. "123-92-2" ([details](https://en.wikipedia.org/wiki/CAS_Registry_Number))

**`CID`**
: e.g. "31276" ([details](https://en.wikipedia.org/wiki/PubChem))

**`SMILES`**
: e.g. "C(C(C)C)COC(=O)C" 
([details](
https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system))


See `?import_new_data` for more details. We can e.g. import data also based on
CAS or CID instead of InChIKey.

#### Looking up InChIKeys If you do not know the InChIKeys of the odorants in
your data set, we recommend using the
[`webchem`](https://cran.r-project.org/package=webchem) package for automated
lookup or doing it manually _via_ <http://cactus.nci.nih.gov/chemical/structure>
or any other chemical lookup service.


## Building the complete data base with
`create_door_database()`{#create_door_database} Once we imported new data we can
use `create_door_database()` in order to rebuild both response matrices. During
the merge process some data sets might be excluded because either their overlap
with other studies is too low or the fit against other studies is too bad; these
studies will be recorded in `door_excluded_data`.


## Updating parts of the data base with
`update_door_database()`{#update_door_database} If we imported new data only for
a few receptors, we can update the data base with `update_door_database()`.
There are two ways to update the data base:

### Using the heuristic approach This is the faster way to perform a merge of
all data sets. All possible binary combinations of data sets will be merged
using 10 different fitting functions on the odorants that were measured in both
data sets. The two data sets yielding the "best merge" (i.e. lowest mean
deviations of points from the fitted function) will be merged. The process of
pairwise merges will be repeated with the "merged_data" against the remaining
data sets until all of them are included:

```{r, fig.width = 7.1, fig.height = 5.5}
update_door_database("Or92a", permutation = FALSE, plot = TRUE)
```

### Trying all permutations The more exhaustive way to update the data base is
to test all possible sequences of data set merges, calculating the mean
deviations from all original data sets and selecting the merge that produces the
lowest mean deviations. This approach works well for responding units that
contain a low number of recorded data sets. For responding units containing 5
data sets we have to calculate merges for 120 different sequences. With 6 it is
already 720 sequences and with 10 data sets we have to test > 3.6 million
different sequences.

While this can be done _via_ parallel computing, this is nothing you should try
on your home PC. For the pre-computed response matrices we performed matches
using the permutation approach for all responding units that contained a maximum
of 10 different data sets on a computing cluster. For DoOR 2.0 these are all
responding units except Or22a.

```{r, fig.width = 7.1, fig.height = 5.5}
update_door_database("Or67a", permutation = TRUE, plot = FALSE)
```


## `model_response()` and `model_response_seq()`{#model} 
`update_door_database()` and `createDatabse()` call `model_response()` and
`model_response_seq()` to perform the merges and update the different DoOR
objects. If we only want to perform a merge we can call them both directly.

### Merging using the heuristic with `model_response()` `model_response()`
returns a list containing the merged data, the names of the excluded data sets
(if any) and the names of the included data sets (if any were excluded).
```{r, fig.width = 7.1, fig.height = 5.5}
merge <- model_response(Or67a, plot = FALSE)
knitr::kable(head(merge$model.response))
```

### Merging in a specific sequence with `model_response_seq()` 
`update_door_database()` with `permutation = TRUE` calls `model_response_seq()`.
Like `model_response()` we can also call model_response_seq directly:
```{r, fig.width=5, fig.height=5.5}
SEQ <- c("Hallem.2006.EN","Kreher.2008.EN","Hallem.2006.EN")
merge <- model_response_seq(Or35a, SEQ = SEQ, plot = TRUE)
head(merge)
```


## Removing a study with `remove_study()`{#remove_study} `remove_study()` will
remove a data set from all DoOR data objects. If we import a data set that
already exists with `import_new_data()`, `remove_study()` will automatically run
before the data is imported.

```{r}
remove_study(study = "Hallem.2004.EN")
```


## Updating the odor information with
`update_door_odorinfo()`{#update_door_odorinfo} If we edit the general odor
information in `DoOR.data::odor` we need to update all other DoOR objects with
the new information. `update_door_odorinfo()` overwrites the first 5 columns of
the DoOR responding units data frames (e.g. `Or22a`), it does not add or remove
lines!



# Helper functions{#helper}
There are several small helper functions that belong to `DoOR.functions`.


## `trans_id()`{#trans_id}
Maybe **the** most important little function in DoOR. With `trans_id()` we can
translate odorant identifiers, e.g. from CAS numbers to InChIKeys or to names.
The information is taken from `DoOR.data::odor`, any `colnames(odor)` can be
used to define input or output:
```{r}
trans_id("123-92-2")
trans_id("123-92-2", to = "Name")
trans_id("carbon dioxide", from = "Name", to = "SMILES")

odorants <- c("carbon dioxide", "pentanoic acid", "water", "benzaldehyde", 
              "isopentyl acetate")
trans_id(odorants, from = "Name", to = "InChI")

```


## `reset_sfr()`{#reset_sfr} `reset_sfr()` subtracts the values of a specified
odorant from a response vector or from the whole response matrix. It is usually
used to subtract the spontaneous firing rate of an odorant, thus setting it to
zero and restoring inhibitory responses. We treat SFR like a normal odorant
during the merging process, thus it becomes > 0 if negative values exist (as all
data gets rescaled `[0,1]` before merging).

`reset_sfr()` works either on the whole `door_response_matrix`, then an odorant
InChIKey has to be specified for subtraction. Or it subtracts a value from a
response vector.

```{r}
rm_sfrReset <- reset_sfr(x = door_response_matrix, sfr = "SFR")
knitr::kable(rm_sfrReset[1:10,6:15], digits = 2)
```

```{r}
reset_sfr(x = c(1:10), sfr = 4)
```


## `door_default_values()`{#door_default_values} `door_default_values()` returns
default values for several parameters used by the DoOR functions, e.g. the
default odor identifier of the colors used in plots.

```{r}
door_default_values("ident")
door_default_values("colors")
```


## `get_responses()`{#get_responses} `get_responses()` returns the response
values of one or several odorants across individual data sets.
```{r}
odorants  <- trans_id(c("carbon dioxide", "isopentyl acetate"), from = "Name")
responses <- get_responses(odorants)
responses <- na.omit(responses)
knitr::kable(head(responses))
```


## `get_normalized_responses()`{#get_normalized_responses} 
`get_normalized_responses()` gathers responses to the specified odorants from
the door_response_matrix and resets the SFR _via_ `reset_sfr()`:
```{r}
odorants  <- trans_id(c("carbon dioxide", "isopentyl acetate"), from = "Name")
responses <- get_normalized_responses(odorants)
responses <- na.omit(responses)
knitr::kable(head(responses))
```


## `countStudies()`{#countStudies} `countStudies()` counts the number of studies
that measured a given odorant-responding unit combination.
```{r}
counts <- countStudies()
knitr::kable(counts[1:10,6:15])
```


## `export_door_data()`{#export_door_data} `export_door_data()` exports all or
selected DoOR data objects in txt or csv format.
```{r}
# export_door_data(".csv")                  	# export all data as .csv files
# export_door_data(".txt", all.data = FALSE) 	# export odorant responses data 
                                              # only as .txt files
```
---
title: "Visualizing DoOR"
author: "Daniel Münch"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Visualizing DoOR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
Apart from merging data and performing calculations on it, the `DoOR.functions`
package provides several ways to visualize the data served by the `DoOR.data`
package. Please see the [DoOR.function vignette](DoOR.functions.html) for
instructions on how to install and load both functions and data. The plotting
functions can be identified by their common prefix `dplot_`, most of them make
use of the `ggplot2` package which allows to override most, if not all of the
theming options.

## Content
* [Loading data](#loading)
* [Visualizing odorants _vs_ responding units with `dplot_response_matrix()`](
#responseMatrix)
* [Visualizing ensemble responses with `dplot_al_map()`](#ALmap)
* [Visualizing tuning curves with `dplot_tuningCurve()`](#tuningCurve)
* [Visualizing response profiles `dplot_response_profile()`](#responseProile)
* [Comparing response profiles `dplot_compare_profiles()`](#compareProfiles)
* [Visualizing odorant responses across responding units with `dplot_across_ru()
   `](#acrossReceptors)
* [Visualizing odorant responses across OSNs with `dplot_across_osns()`](
   #across_OSNs)


## Loading data{#loading}
First we need to load packages and data:
```{r, echo = T, message = TRUE, results='hide'}
#load data
library(DoOR.functions)
library(DoOR.data)
load_door_data(nointeraction = TRUE)
```

## Visualizing odorants _vs_ responding units with `dplot_response_matrix()`
{#responseMatrix}
`dplot_response_matrix()` visualizes the DoOR consensus response matrix either
as a point matrix with the size of the points relating to the response strength,
or as a heatmap, relating color to response strength. The default plot that is
used depends on the data we enter: positive values only (like from
`door_response_matrix` or `door_response_matrix_non_normalized`) which are
scaled `[0,1]` will by default be plotted as a point-matrix. If we provide data
that contains negative values (e.g. data with the spontaneous firing rate set to
0 via `reset_sfr(door_response_matrix, "SFR")`) it will be shown as color-coded
heatmap.

```{r,fig.width=7.1, fig.height=4}
dplot_response_matrix(door_response_matrix[2:50,], tag = "Name", base_size = 8)
```

If the data contains negative values, a colored response matrix will be plotted:
```{r,fig.width=7.1, fig.height=4}
dplot_response_matrix(reset_sfr(door_response_matrix, "SFR")[2:50,], 
                      tag = "Name", base_size = 8)
```

## Visualizing ensemble responses with `dplot_al_map()`{#ALmap}
With `dplot_al_map()` we can visualize the ensemble response a given odorant
elicits across DoOR responding units (receptors, sensory neurons, glomeruli, ...
) as a hypothetical antennal lobe activation pattern.
```{r, fig.width=7.1}
dplot_al_map("QSJXEFYPDANLFS-UHFFFAOYSA-N", base_size = 8)
```
If we do not know the InChIKey of the substance we are interested in, we can use
the `trans_id()` function for conversion. If we are interested in the expressed
receptor rather than the glomerulus names:
```{r, fig.width=7.1, warning=FALSE}
dplot_al_map(trans_id("benzaldehyde", from = "Name"), tag = "receptor", 
             main = "SMILES", base_size = 8)
```

If we prefer the plain activation pattern without annotations at all:
```{r, fig.width=7.1, fig.height=2}
dplot_al_map(trans_id("628-63-7", from = "CAS"), tag = "", main = "", 
             legend = FALSE, scalebar = FALSE)
```

## Visualizing tuning curves with `dplot_tuningCurve()`{#tuningCurve}
### Responding units

The set of odorants a given responding unit is responsive to can be described as
its tuning curve. With `dplot_tuningCurve()` we can easily display such a tuning
curve together with its kurtosis. Kurtosis is a measure of the shape of the
tuning curve, i.e. whether a responding unit is narrowly tuned to a few odorants
(high kurtosis), or whether it responds to many odorants (low kurtosis). The
gaussian distribution has a kurtosis of 0.

If we only specify the receptor/response unit name, the data is taken from
`door_response_matrix` and `dplot_tuningCurve()` uses `reset_sfr()` to reset
spontaneous firing rate (or any other specified odorant) to zero before
plotting.
```{r}
dplot_tuningCurve(receptor = "Or22a", base_size = 8)
```

To prevent resetting by SFR:
```{r}
dplot_tuningCurve(receptor = "Or22a", zero = "", base_size = 8)
```

We can also plot any other vector of responses, once `response.vector` is
specified, the value of `receptor` is only used for the plot title and not for
data lookup anymore:
```{r}
dplot_tuningCurve(receptor = "receptor X", response.vector = c(1:100), 
                  base_size = 8)
```

### Odorants
`dplot_tuningCurve()` can as well be used to visualize the ensemble of
responding units that is activated by a given odorant. Therefore, we specify an
odorant name instead of a response unit name:
```{r}
dplot_tuningCurve(odorant = "PGMYKACGEOXYJE-UHFFFAOYSA-N", base_size = 8)
```

We can specify the chemical identifier to plot via `odor.main`:
```{r}
dplot_tuningCurve(odorant = "PGMYKACGEOXYJE-UHFFFAOYSA-N", odor.main = "SMILES",
                  base_size = 8)
dplot_tuningCurve(odorant = "CURLTUGMZLYLDI-UHFFFAOYSA-N", odor.main = "InChI",
                  base_size = 8)
```

And finally we can control the color of the bars:
```{r}
dplot_tuningCurve(odorant = trans_id("carbon dioxide", from = "Name"), 
                  fill.odorant = "#FF0000", base_size = 8)
```

As mentioned, all of these plots are generated with the ggplot2 package which
allows to override theming:
```{r}
library(ggplot2)
dplot_tuningCurve(odorant = trans_id("carbon dioxide", from = "Name"), 
                  base_size = 8) +
  theme(panel.background = element_rect(fill = "grey", color = "magenta"))
```

## Visualizing response profiles `dplot_response_profile()`{#responseProile}
`dplot_response_profile` creates a horizontal bar plot of the response profile
of a given receptor. It displays the same data as `dplot_tuningCurve` but
focusses on the odorant identity.

Per default the response strength is displayed as bar height as well as color
code:
```{r, fig.width=5, fig.height=6}
dplot_response_profile("Gr21a.Gr63a", tag = "Name", base_size = 8)
```

We can again omit to reset to SFR:
```{r, fig.width=5, fig.height=6}
dplot_response_profile("Gr21a.Gr63a", tag = "Name", base_size = 8, zero ="")
```

And if we prefer monochrome data:
```{r, fig.width=5, fig.height=6}
dplot_response_profile("Gr21a.Gr63a", tag = "CAS", base_size = 8, 
                       colored = FALSE)
```

## Comparing response profiles `dplot_compare_profiles()`{#compareProfiles}
As the name indicates, with `dplot_compare_profiles()` we can plot two response
profiles side by side for comparison. The syntax here differs from the previous
plots as we can use it also to compare the original data sets in DoOR.


### Comparing original data sets
This time, with `x` and `y` we have to specify a whole data.frame, `by.x` and
`by.y`  take the corresponding column names that will be plotted. If `x` is not
specified, both `by.x` and `by.y` will be taken from `x`.

```{r}
dplot_compare_profiles(x = Or22a, y = Or22a, by.x = "Pelz.2006.AntEC50",
                         by.y = "Hallem.2004.EN", tag = "Name", base_size = 8)
```

We see that the measured Or22a sensory neuron responses in these two data sets
are in good accordance.

### Comparing DoOR response profiles
Next we compare two DoOR consensus response profiles, measurements of the
misexpressed receptor Or35a and recordings from the ac3B sensory neuron that
naturally expresses Or35a:
```{r, fig.width = 7.1}
dplot_compare_profiles(
  x = door_response_matrix,
  by.x = "Or35a",
  by.y = "ac3B",
  tag = "Name",
  base_size = 8
  )
```

or with "SFR" set to 0:
```{r, fig.width = 7.1}
dplot_compare_profiles(
  x = reset_sfr(door_response_matrix, "SFR"),
  by.x = "Or35a",
  by.y = "ac3B",
  tag = "Name",
  base_size = 8
  )
```

We see that the response profiles are very similar but not identical, which is
expected as the ac3B neuron expresses other receptors in addition to Or35a.

## Visualizing responses across responding units with `dplot_across_ru()`
{#acrossReceptors}
With `dplot_across_ru()` we can visualize the responses that one or several
odorants elicit across receptors / responding units.

```{r, fig.width = 6, fig.height = 6}
odors <-
  trans_id(c("pentyl acetate", "carbon dioxide", "2,3-butanedione"), 
           from = "Name")
  dplot_across_ru(odors, tag = "Name", base_size = 8)
```


## Visualizing odorant responses across OSNs with `dplot_across_osns()`
{#across_OSNs}
`dplot_across_osns()` is similar to `dplot_across_ru()` but the responding units
are sorted according to the sensory neuron they belong to. This sorting is
controlled via `door_mappings` from the `DoOR.data` package. There are two types
of plot that `dplot_across_osns()` can return. Type 2 directly resembles
`dplot_across_ru()`:


```{r, fig.width = 6, fig.height = 6}
dplot_across_osns(odors, base_size = 8, plot.type = 2)
```

In type 1 the data is split according to odorant X sensillum, the color is
assigned to the corresponding neuron A-D or X-Z if the neuron's identity is
unknown:
```{r, fig.width = 7.1, fig.height = 5}
dplot_across_osns(odors, base_size = 8, plot.type = 1)
```

As this plot gets pretty messy, we have the option to restrict plotting to
certain subsets of sensilla. We can for example only plot the responses of
antennal basiconic (ab) sensilla:

```{r, fig.width = 7.1, fig.height = 5}
dplot_across_osns(odors, base_size = 8, plot.type = 1, sub = "ab")
```

Or we plot coeloconic and trichoid sensilla:

```{r, fig.width = 7.1, fig.height = 5}
dplot_across_osns(odors, base_size = 8, plot.type = 1, sub = c("ac", "at"))
```
---
title: "DoOR analysis tools"
author: "Daniel Münch"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{DoOR analysis tools}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The `DoOR.functions` package provides tools for the analysis of the data
provided by the `DoOR.data` package.

## Content
* [Loading data](#loading)
* [Identifying the sensillum we are recording from with `identify_sensillum()`](
#identify_sensillum)
* [Finding neuron-specific odorants with `private_odorant()`](#private_odorant)
* [Mapping response data from an unknown source with `map_receptor()`](
#map_receptor)
* [Changing the response unit with `back_project()`](#back_project)
* [Quantifying tuning curves with `sparse()`](#sparse)



## Loading data{#loading}
First we need to load packages and data:
```{r, echo = TRUE, message = TRUE, results='hide'}
#load data
library(DoOR.functions)
library(DoOR.data)
load_door_data(nointeraction = TRUE)
```

## Identifying the sensillum we are recording from with `identify_sensillum()` 
{#identify_sensillum}
Imagine we perform an electrophysiological recording from a _Drosophila_
sensillum (single sensillum recording, SSR) and we are not sure what sensillum
we are recording from. In order to identify the sensillum we used several
diagnostic odorants (maybe selected using
[`private_odorant()`](#private_odorant)) and got responses from the different
sensory neurons the sensillum houses. We can now pass our recorded data to
`identify_sensillum()`.

Let's make up some simple fake data. We pretend to have recorded with three
odorants (2,3-butanedione, ethanoic acid and carbon dioxide) and we could
separate the responses of two units. Unit1 responded strongly to 2,3-butanedione
only, unit2 only responded to carbon dioxide. We create a data.frame that
contains a column called `odorants` with the InChIKeys of our test odorants, and
one column for each unit (name the colnames as you like, e.g. unit1-n or Aneuro
if you are sure about the neuron).

```{r}
recording <- data.frame(
  odorants = c(trans_id(c("BEDN", "ETAS"), "Code"),
               trans_id("carbon dioxide", "Name")),
  unit1 = c(.9,.1,.1),
  unit2 = c(0, .1, 1)
)
```
Next we feed the recording to the function:

### using correlations
```{r, fig.width=7.1, fig.height=4.5}
identify_sensillum(recording, base_size = 8)
```
Note that the function tells us that it found hits for all units in ab1 and ab5,
meaning that e.g. within the four neurons housed in the ab1 sensillum both of
our units had good matches. You can set this correlation threshold with
`min.cor`. If we increase the threshold to 0.99 only ab1 is returned as a double
match:

```{r, fig.width=7.1, fig.height=4.5, fig.show='hide'}
identify_sensillum(recording, min.cor = .99)
```

We can define the number of best hits that we want to get returned (the default
is 10):
```{r, fig.width=7.1, fig.height=4.5}
identify_sensillum(recording, nshow = 5, base_size = 8)
```
And if we know e.g. that we are recording from a basiconic sensillum we can
restrict the search to one or a few sensilla types:
```{r, fig.width=7.1, fig.height=4.5}
identify_sensillum(recording, sub = "ab", nshow = 5, base_size = 8)
identify_sensillum(recording, sub = c("ac","at"), nshow = 5, base_size = 8)
```

### using Euclidean distances
Instead of correlations we can also use the Euclidean distance as a
(dis)similarity measure:

```{r, fig.width=7.1, fig.height=4.5}
identify_sensillum(recording, method = "dist", sub = "ab", nshow = 5, 
                   base_size = 8)
```

### returning data instead of plots
We can also return the correlation/distance data instead of the plot when
setting `plot =FALSE`:
```{r}
sensillumX <-
  identify_sensillum(recording,
  method = "dist",
  sub = "ab",
  plot = FALSE)
  head(sensillumX)
```

So apparently our fake recording came from the ab1 sensillum, which was
admittedly quite obvious as we had a strong carbon dioxide response and ab1
houses the carbon dioxide receptor :)



## Finding neuron-specific odorants with `private_odorant()` {#private_odorant}
There may be several cases where we might be interested in so called 
_private odorants_, odorants that specifically activate a given receptor or 
sensory neuron. Maybe we are looking for diagnostic odorants for sensillum
identification or we want to activate a specific neuronal pathway,
`private_odorant()` returns candidate odorants for that task.

Let's say we want to specifically activate Or22a neurons:
```{r}
private_odorant("Or22a")
```

We might want to return the odorant names instead of InChiKeys:
```{r}
private_odorant("Or22a", tag = "Name")
```
So according to the function sec-amyl acetate would be a good candidate. It
activates Or22a at 0.4 (DoOR response, max is 1) while the maximum activation in
all other tested responding units (receptors, neurons, glomeruli) is 0.016, a
difference of 0.40. Sounds good, but it was tested only in 4 other responding
units, so I would rather go for ethyl hexanoate with about the same difference
but being tested in 29 other responding units.

We can also restrict the search to the sensillum the responding units of
interest is related to:
```{r}
private_odorant("Or22a", tag = "Name", sensillum = T)
```
Ethyl 2-methylbutanoate sounds like a good hit, it has the same difference to
the other units as ethyl hexanoate but hardly elicits a response at all from the
other neuron. The n of 1 is fine as there are only 2 neurons housed in the ab3
sensillum.



## Mapping response data from an unknown source with `map_receptor()`
{#map_receptor}
Similar to `identify_sensillum()`, `map_receptor()` correlates a response vector
to all responding units of the existing DoOR consensus data. Let's grab a data
set from Or22a and see where it ends up:
```{r}
data <- data.frame(odorants  = Or22a$InChIKey, responses = Or22a$Hallem.2006.EN)
data <- na.omit(data)
head(data)
map_receptor(data = data, nshow = 5)
```
This example was a bit circular as the tested data contributed to the consensus
data...



## Changing the response unit (spikes, deltaF/F, ...) with `back_project()`
{#back_project}
The DoOR consensus data is normalized to values between 0 and 1. If we want to
compare the DoOR data to one of our own recordings, it would be great to have
the DoOR data in the same unit as our owne data (e.g. spikerate). This is what
we can do with `back_project()`,  we can rescale the DoOR data to fit a given
response template.

As an example, let's take the data Hallem et al. recorded from Or22a _via_
calcium imaging and rescale the DoOR responses accordingly. The template has to
have 2 columns named `odorants` and `responses`:
```{r, fig.width = 4.5, fig.height = 4.5}
template <- data.frame(odorants  = Or22a$InChIKey,
                       responses = Or22a$Hallem.2006.EN)

bp <- back_project(template, responding.unit = "Or22a")
plot(bp$back_projected$original.data,
     bp$back_projected$back_projected.data,
     xlab = "DoOR consensus response",
     ylab = "back_projected data [spikes, Hallem.2006.EN]"
)

head(bp$back_projected)

```
All the yellow lines in the first plot represent odorant responses that were not
available in the original data set but were projected onto the fitted function
and rescaled to the units in "Hallem.2006.EN". The second plot shows the
relationship between the rescaled data and the original consensus responses.


## Quantifying tuning curves with `sparse()`{#sparse}
The width of a tuning curve, i.e. for example to how many odorants a receptor
shows strong responses, can be quantified using different sparseness measures.
We implemented kurtosis[^1] and sparseness[^2] in `sparse()`. A high kurtosis or
sparseness value indicates a narrow tuning curve (see also `dplot_tuningCurve()`
in the [DoOR visualization vignette](DoOR_visualizations.html#tuningCurve)).
While kurtosis is able to deal with negative values, sparseness can't and thus
all values need to be be transformed to absolute values first. Sparseness scales
between 0 and 1, kurtosis between -∞ and ∞, a kurtisis of 0 corresponds to the
Gaussian distribution.

```{r}
rm.SFRreset <- reset_sfr(door_response_matrix, "SFR")

sparse(x = rm.SFRreset[,"Or69a"], method = "ltk")
sparse(x = rm.SFRreset[,"Or69a"], method = "lts")

sparse(x = rm.SFRreset[,"Gr21a.Gr63a"], method = "ltk")
sparse(x = abs(rm.SFRreset[,"Gr21a.Gr63a"]), method = "lts")

```



[^1]: Willmore, B., Tolhurst, D.J., 2001. Characterizing the sparseness of
neural codes. Network 12, 255–270. dx.doi.org/10.1080/net.12.3.255.270
[^2]: Bhandawat, V., Olsen, S.R., Gouwens, N.W., Schlief, M.L., Wilson, R.I.,
2007. Sensory processing in the Drosophila antennal lobe increases reliability
and separability of ensemble odor representations. Nature neuroscience 10,
1474–82. dx.doi.org/10.1038/nn1976
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_door_odorinfo.R
\name{update_door_odorinfo}
\alias{update_door_odorinfo}
\alias{updateOdorInfo}
\title{update_door_odorinfo}
\usage{
update_door_odorinfo()
}
\description{
Update the DoOR odor data with info from \code{odor}. For the function to
work, all DoOR data has to be loaded to the current environment.
}
\examples{
# load data
load_door_data(nointeraction = TRUE)  
# modify odor
odor[1,1] <- "acid"

# run 
update_door_odorinfo()

# check that data sets have been updated
head(Or22a)


}
\author{
Daniel Münch, \email{daniel@muench.bio}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_normalized_responses.R
\name{get_normalized_responses}
\alias{get_normalized_responses}
\alias{getNormalizedResponses}
\title{Find normalised receptor responses}
\usage{
get_normalized_responses(odors, zero = door_default_values("zero"),
  response_matrix = door_default_values("door_response_matrix"), round = 3,
  na.rm = FALSE)
}
\arguments{
\item{odors}{character vector, one or more odors provided as InChIKey.}

\item{zero}{InChIKey of background that should be set to zero. The default is
"SFR", i.e. the spontaneous firing rate.}

\item{response_matrix}{a data frame, as e.g. "door_response_matrix" that is
loaded by \code{\link{model_response}}. It is also possible to create this
frame manually using \code{\link{model_response}}.}

\item{round}{numeric, round to this amount of digits, set to NA if you do not
want to round}

\item{na.rm}{logical, remove NAs?}
}
\description{
given a chemical, get normalised receptor responses from all studies in the 
database.
}
\examples{
# load data
library(DoOR.data)
data(door_response_matrix)

# define a list of odorants
odors <- c("MLFHJEHSLIIPHL-UHFFFAOYSA-N",
           "OBNCKNCVKJNDBV-UHFFFAOYSA-N",
           "IKHGUXGNUITLKF-UHFFFAOYSA-N")

# get the normalized responses for these odorants
result <- get_normalized_responses(odors, 
                                   response_matrix = door_response_matrix)

}
\seealso{
\code{\link{model_response}},\code{\link{create_door_database}}
}
\author{
Daniel Münch \email{daniel.muench@uni-konstanz.de}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplot_tuningCurve.R
\name{dplot_tuningCurve}
\alias{dplot_tuningCurve}
\title{dplot_tuningCurve}
\usage{
dplot_tuningCurve(receptor, odorant, response.vector,
  response_matrix = door_default_values("door_response_matrix"),
  odor_data = door_default_values("odor"),
  zero = door_default_values("zero"),
  fill.receptor = door_default_values("color.receptor"),
  fill.odorant = door_default_values("color.odorant"), odor.main = "Name",
  limits, base_size = 12)
}
\arguments{
\item{receptor}{character, a receptor name (one of ORS$OR)}

\item{odorant}{character, an odorant name (InChIKey)}

\item{response.vector}{numerical vector, a vector with responses, if empty 
this is taken from door_response_matrix}

\item{response_matrix}{DoOR response matrix, response vector will be taken 
from here, not needed if response.vector is given}

\item{odor_data}{data frame, contains the odorant information.}

\item{zero}{InChIKey, will be set to zero, default is SFR, ignored when data 
is provided via response.vector, set to "" if you don't want to subtract 
anything}

\item{fill.receptor}{color code, bar color for receptor tuning curve}

\item{fill.odorant}{color code, bar color for odorant tuning curve}

\item{odor.main}{the odor identifier to plot, one of colnamed(odor)}

\item{limits}{the numerical vector of length 2, y limits for the tuning curve}

\item{base_size}{numeric, the base font size for the ggplot2 plot}
}
\value{
a ggplot object
}
\description{
plot a receptor or odorant tuning curve
}
\examples{
# load data
library(DoOR.data)
data(door_response_matrix)
data(odor)

# plot a tuning curve for an odorant
dplot_tuningCurve(odorant = odor$InChIKey[2])
# or for a receptor
dplot_tuningCurve(receptor = "Or22a")

# adjust the plotting range
range <- range(reset_sfr(door_response_matrix, "SFR"), na.rm = TRUE)
dplot_tuningCurve(receptor = "Or10a", limits = range, 
                  fill.receptor = "magenta")

# plot with manual input data as receptor tuning curve
dplot_tuningCurve(receptor = "receptor X", response.vector = c(1:100))
# or as odor tuning curve
dplot_tuningCurve(odorant = "odor X", response.vector = rnorm(200))

}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reset_sfr.R
\name{reset_sfr}
\alias{reset_sfr}
\alias{resetSFR}
\title{reset SFR}
\usage{
reset_sfr(x, sfr)
}
\arguments{
\item{x}{numeric or DoOR response matrix, input values}

\item{sfr}{numeric or character, either a value to subtract if x is a vector
or an InChIKey if x is a DoOR response matrix}
}
\description{
A function for reseting SFR to zero
}
\details{
Performs a simple subtraction of the SFR value.
}
\examples{
# load data
library(DoOR.data)
data(door_response_matrix)

# create a response matrix with the SFR reset to 0
door_response_matrix_SFRreset <- reset_sfr(door_response_matrix, "SFR")

}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export_door_data.R
\name{export_door_data}
\alias{export_door_data}
\alias{exportData}
\title{export data}
\usage{
export_door_data(file.format, directory,
  odorantReceptors = door_default_values("ORs"),
  response_matrix = door_default_values("door_response_matrix"),
  responseRange = door_default_values("door_response_range"),
  unglobalNorm_RM = door_default_values("door_response_matrix_non_normalized"),
  weightGlobNorm = door_default_values("door_global_normalization_weights"),
  all.data = TRUE)
}
\arguments{
\item{file.format}{character string, the format of given file, either ".txt" 
or ".csv"}

\item{directory}{character string, naming a directory for writing. If 
missing, the exported data are saved in current working directory.}

\item{odorantReceptors}{data frame, receptor names and expressions}

\item{response_matrix}{data matrix, an global unnormalized responses matrix}

\item{responseRange}{data frame, response ranges for each study}

\item{unglobalNorm_RM}{data matrix, an unnormalized responses matrix}

\item{weightGlobNorm}{data frame, weight matrix for global normalizazion}

\item{all.data}{logical, if TRUE, export odorant response data and supported 
data "door_response_matrix", "door_response_range", 
"door_response_matrix_non_normalized", "door_response_matrix", 
"door_global_normalization_weights" and "ORs".}
}
\description{
export odor response data and supported data
}
\details{
Please load ORs from data package DoOR.data by typing (\code{data(ORs)}) 
before use.
}
\examples{
\dontrun{
# load data
library(DoOR.data)
library(DoOR.functions)
load_door_data()

# export odorant response data only
export_door_data(".txt", all.data = FALSE) 	
}
}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_model.R
\name{select_model}
\alias{select_model}
\alias{selectModel}
\title{compute the data pairwise and and selects a pair with the lowest "MD" value.}
\usage{
select_model(candidate, data_candidate, merged_data,
  overlapValues = door_default_values("overlapValues"),
  merged = door_default_values("merged"))
}
\arguments{
\item{candidate}{a character vector, contains the names of studies.}

\item{data_candidate}{a data frame, odorant response data that only contains
value columns.}

\item{merged_data}{numeric vector, merged data}

\item{overlapValues}{numeric, a criterion using to refuse a data set that has
not enough overlap value.}

\item{merged}{logical, if merged is \code{TRUE}, calculate models between
merged_data and candidate data. If \code{FALSE}, calculate models between
candidates.}
}
\description{
compute the data pairwise using function \code{\link{calculate_model}} and
selects a pair with the lowest "MD" value.
}
\details{
This function is used in \code{\link{model_response}} to select the first
pair or next data set for merging. The output is a list containing
"selected.x" and "selected.y" that specify which data plots against another,
and "best.model" that is used in function \code{\link{project_points}}.
}
\examples{
# load data
library(DoOR.data)
data(ac3B)

# split into data and header
studies <- names(ac3B)[c(7:8)]
data_candidate <- ac3B[,c(7:8)]

# rescale data
norm_data_candidate <- apply(data_candidate, 2, door_norm)

# find the best fitting model
select_model(candidate = studies, data_candidate = norm_data_candidate,
             merged = FALSE)

}
\seealso{
\code{\link{project_points}},\code{\link{model_response}}
}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplot_across_ru.R
\name{dplot_across_ru}
\alias{dplot_across_ru}
\alias{dplot_acrossReceptors}
\title{dplot_across_ru}
\usage{
dplot_across_ru(odorant,
  response_matrix = door_default_values("door_response_matrix"),
  odor_data = door_default_values("odor"),
  zero = door_default_values("zero"), tag = "Name", limits,
  base_size = 12)
}
\arguments{
\item{odorant}{character vecto, one or several InChIKeys}

\item{response_matrix}{DoOR response matrix, a DoOR response matrix as data
source}

\item{odor_data}{data frame, contains the odorant information.}

\item{zero}{character, an InChIKey of the odorant that should be set to 0}

\item{tag}{character, the chemical identifier to plot as odorant name (one of
colnames(odor))}

\item{limits}{numeric of length 2, if provided the ylim will range
accordingly}

\item{base_size}{numeric, the base font size for the ggplot2 plot}
}
\value{
a ggplot object
}
\description{
barplot of DoOR responses of a set of odorant across all responding units in
DoOR
}
\examples{
# load data
library(DoOR.data)
library(DoOR.functions)
data(odor)

# plot activation pattern of one or several odorants
dplot_across_ru(trans_id("123-92-2"), tag = "CAS")
dplot_across_ru(odor$InChIKey[4:10])

}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DoOR.functions-package.R
\docType{package}
\name{DoOR.functions.package}
\alias{DoOR.functions.package}
\alias{DoOR.functions}
\alias{DoOR.function}
\alias{DoOR.functions.package-package}
\title{DoOR Functions}
\description{
Functions package providing manipulation and application of the DoOR.
}
\details{
\tabular{ll}{ Package: \tab DoOR.functions\cr Type: \tab Package\cr Version: 
\tab 2.0.0\cr Date: \tab 2016-01-25\cr License: \tab GPL-3\cr LazyLoad: \tab 
yes\cr }

\bold{Type \code{help(package = DoOR.functions)} to see a complete list of 
datasets and functions.  Below is what you need for a quick start.}

First, load the DoOR packages, data and function package: \tabular{ll}{ 
\code{library(DoOR.functions)}: \tab \cr \code{library(DoOR.data)}: \tab \cr 
}

then, load all datasets including the precomputed response matrix: 
\tabular{ll}{ \code{load_door_data}: \tab Load all data into current active 
environment (function comes with DoOR.data) . \cr } or, load all odorant
reseponse data into a list: \tabular{ll}{ \code{\link{load2list}}: \tab Load
odorant response data only and compose them as a list.  \cr }

Try some visualizations (e.g. producing the plots from the paper): 
\tabular{ll}{ \code{\link{dplot_al_map}}: \tab response to a chemical mapped
onto an image of the antennal lobe.\cr \code{\link{dplot_compare_profiles}}:
\tab compare the results of two studies.  \cr
\code{\link{dplot_response_matrix}}: \tab Dot Plot of Odorant Responses
Across Receptors. \cr \code{\link{dplot_response_profile}}: \tab bar plot:
one receptor, all chemicals. \cr \code{\link{dplot_tuningCurve}}: \tab
pyramid diagram depicting a receptor's tuning breadth. \cr } Try some
queries: \tabular{ll}{ \code{\link{get_responses}}: \tab given a chemical,
get original responses from all studies in the database.\cr
\code{\link{get_normalized_responses}}: \tab given a chemical, get normalised
responses from all studies in the database.\cr}

In case you wish to create your own response model (e.g. because you want to 
include your own datasets): \tabular{ll}{ \code{\link{create_door_database}}:
\tab compute the complete response model for all receptors in the database
(calls \code{\link{model_response}} for all receptors). \cr
\code{\link{model_response}}: \tab run the DoOR algorithm, that merges all
measurements for one receptor. \cr }

Estimate odorant responses: \tabular{ll}{
\code{\link{estimate_missing_value}}: \tab estimate NA entries in a consensus
response data. \cr } Project the model response values back to tested values:
\tabular{ll}{ \code{\link{back_project}}: \tab project the model response
values back to tested values. \cr }

Introduce new data into DoOR and update the supported data sets: 
\tabular{ll}{ \code{\link{import_new_data}}: \tab import new data into DoOR, 
and update the weight, response range and receptor names. \cr 
\code{\link{update_door_database}}: \tab update response matrix by
calculating new consensus response data for a given receptor. \cr }

See the Vignettes and the help pages for more documentation.
}
\references{
\url{http://neuro.uni-konstanz.de/DoOR}
}
\seealso{
\code{DoOR.data}
}
\author{
C. Giovanni Galizia \cr Daniel Muench \cr Martin Strauch \cr Anja 
  Nissler \cr Shouwen Ma \cr
  
  Maintainer: Daniel Münch <daniel.muench@uni-konstanz.de>
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_studies.R
\name{count_studies}
\alias{count_studies}
\alias{countStudies}
\title{count_studies}
\usage{
count_studies(ors = door_default_values("ORs"),
  odor_data = door_default_values("odor"),
  char.columns = door_default_values("num.charColumns"),
  ident = door_default_values("ident"))
}
\arguments{
\item{ors}{data.frame containing all receptors exidting in DoOR.}

\item{odor_data}{data.frame containing information about the odorants in
DoOR.}

\item{char.columns}{number of character columns in each receptor data.frame.}

\item{ident}{odorant identifier to be used as rownames.}
}
\value{
matrix
}
\description{
returns a matrix indiating how many studies have recorded individual
receptor-odorant combinations
}
\examples{
# load some data
library(DoOR.data)
load_door_data(nointeraction = TRUE)

#run count studies and plot the result
count <- count_studies()
image(count)
head(count)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_study.R
\name{remove_study}
\alias{remove_study}
\alias{removeStudy}
\title{Remove a study from DoOR}
\usage{
remove_study(study, receptors = door_default_values("ORs"),
  responseRange = door_default_values("door_response_range"),
  weightGlobNorm = door_default_values("door_global_normalization_weights"))
}
\arguments{
\item{study}{a string containing the name of the study you want to remove
(e.g. 'Bruyne.2001.WT')}

\item{receptors}{a vector of all the receptors to be checked. Defaults to all
receptors exidting in DoOR.}

\item{responseRange}{the dataframe containing the info about the response
ranges of all studies (\code{door_response_range})}

\item{weightGlobNorm}{the dataframe containing the info about the relative
weights between receptors (\code{door_global_normalization_weights})}
}
\description{
Use this function to remove a study from the DoOR database.
\code{import_new_data.R} uses this function when it detects an existing study
during the import process (e.g. because you imported updated data).
}
\examples{
# load data
library(DoOR.data)
load_door_data(nointeraction = TRUE)

# remove Bruyne.2001.WT from DoOR
remove_study('Bruyne.2001.WT')

}
\seealso{
\code{\link{import_new_data}}
}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplot_response_matrix.R
\name{dplot_response_matrix}
\alias{dplot_response_matrix}
\alias{dplot_responseMatrix}
\title{dplot_response_matrix}
\usage{
dplot_response_matrix(data, odor_data = door_default_values("odor"),
  tag = door_default_values("tag"), colors = door_default_values("colors"),
  flip = FALSE, fix = TRUE, bw = FALSE, point = FALSE, limits,
  base_size = 12)
}
\arguments{
\item{data}{a subset of e.g. door_response_matrix}

\item{odor_data}{data frame, contains the odorant information.}

\item{tag}{the chemical identfier to plot (one of colnames(odor))}

\item{colors}{the colors to use if negative values are supplied (range of 5 
colors, 2 for negative values, 1 for 0 and 3 for positive values)}

\item{flip}{logical, if TRUE the x and y axes will be flipped}

\item{fix}{logical, whether to fix the ratio of the tiles when plotting as a 
heatmap}

\item{bw}{logical, whether to plot b&w or colored}

\item{point}{logical, if \code{TRUE} a point matrix instead of a heatmap will
be returned (the default if you supply only positive values)}

\item{limits}{the limits of the scale, will be calculated if not set}

\item{base_size}{numeric, the base font size for the ggplot2 plot}
}
\value{
a dotplot if limits[1] >= 0  or a heatmap if limits[1] < 0
}
\description{
plot DoOR responses as a point matrix
}
\examples{
# load data
library(DoOR.data)
data(door_response_matrix)

# reset the spontaneous firing rate to 0
tmp <- reset_sfr(door_response_matrix, "SFR")

# plot heatmap / coloured tiles
dplot_response_matrix(tmp[10:50,], tag = "Name", 
 limits = range(tmp, na.rm = TRUE))

# plot dotplot
dplot_response_matrix(door_response_matrix[10:50,], tag = "Name",
                        limits = range(door_response_matrix, na.rm = TRUE))

}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private_odorant.R
\name{private_odorant}
\alias{private_odorant}
\alias{privateOdorant}
\title{private_odorant}
\usage{
private_odorant(receptor, sensillum = FALSE,
  response_matrix = door_default_values("door_response_matrix"),
  door_mappings = door_default_values("door_mappings"),
  zero = door_default_values("zero"), nshow = 5, tag)
}
\arguments{
\item{receptor}{character, name of a DoOR responding unit (one of
\code{ORs$Or})}

\item{sensillum}{logical, restrict the search to the sensillum the receptor
is expressed in?}

\item{response_matrix}{DoOR response matrix, the input data to perform the
search on}

\item{door_mappings}{the data frame containing the mapping information}

\item{zero}{character, an odorant that should be set to 0}

\item{nshow}{numeric, the number of private odorants to return}

\item{tag}{character, the chemical identifier to give the odorant names in
(on of \code{colnames(odor)})}
}
\value{
a data.frame containing odorants and the response in the receptor of
  interest as well as the maximum response of the remaining receptors and
  their difference
}
\description{
return an odorant that activates the receptor of interest exclusively
}
\examples{
# load data
library(DoOR.data)

# find a private odorant for Gr21a.Gr63a (the carbon dioxide receptor)
# private_odorant("Gr21a.Gr63a", tag = "Name")

# now find an odorant that within the ab3 sensillum specifically activates
# Or22a
private_odorant("Or22a", tag = "Name", sensillum = TRUE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_missing_value.R
\name{estimate_missing_value}
\alias{estimate_missing_value}
\alias{DoOREst}
\title{Estimate the missing entries in a response data}
\usage{
estimate_missing_value(data, nodor, method = "PC")
}
\arguments{
\item{data}{a data frame or matrix, contaning the consensus response values}

\item{nodor}{a numeric value, specifying the number of the selected odors}

\item{method}{character string, specifying the method ("PC" (Pearson's
coefficient) and "Knn" (k nearest neighbors)) for estimation, the default is
"PC".}
}
\description{
Estimate the missing entries in a response data
}
\details{
A wrapper programe for using Pearson Correlation or k-nearest neighbors to
estimate the missing entries in a response matrix.
}
\examples{
\dontrun{
# load data
library(DoOR.data)
data(door_response_matrix)

# pick an example subset
subset <- door_response_matrix[1:100, 11:30]

# estimate missing values
est.data <- estimate_missing_value(data = subset, nodor = 6)
}
}
\references{
Kim, H., Golub, G. H. & Park, H., Missing value estimation for
DNA microarray gene expression data: local least squares imputation., 2005,
Bioinformatics, 21, 187-198
}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>
}
\keyword{math}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_response.R
\name{model_response}
\alias{model_response}
\alias{modelRP}
\title{Generates a model response}
\usage{
model_response(da, select.MDValue = door_default_values("select.MDValue"),
  overlapValues = door_default_values("overlapValues"),
  responseRange = door_default_values("door_response_range"),
  weightGlobNorm = door_default_values("door_global_normalization_weights"),
  glob.normalization = door_default_values("glob.normalization"),
  plot = door_default_values("plot"))
}
\arguments{
\item{da}{data frame, a selected receptor containing measured responses from
studies.}

\item{select.MDValue}{numeric, threshold on the MD for rejecting a fit.}

\item{overlapValues}{numeric, a criterion using to refuse a data set that
has not enough overlap value.}

\item{responseRange}{data frame, contains response ranges for all studies.}

\item{weightGlobNorm}{data frame, a binary data matrix, 1 indicates given
odor has been measured in given study, NA indicates NOT.}

\item{glob.normalization}{logical, default is \code{TRUE}, performs a global
normalization for the model response. Otherwise (\code{FALSE}) response
values will be given in value from 0 to 1.}

\item{plot}{logical, If \code{FALSE}, plotting is suppressed. Default is
\code{FALSE}.}
}
\description{
Runs the DoOR algorithm, that merges all measurements for one receptor into
a common response model.
}
\details{
Merging a data is processed by following: \enumerate{ \item Normalize all
response data in value [0,1].  \item Compute the correlation between studies
and selected the best pair using \code{\link{select_model}}.  \item Merge the
first pair using function \code{\link{project_points}}.  \item Add other
datasets if the correlation between the growing model response and the new
dataset is below the correlation threshold (select.MDValue). Datasets
excluded based on this criterion will be appended in a separate list.  }
}
\examples{
# load data
library(DoOR.data)
data(Or35a)
data(door_global_normalization_weights)
data(door_response_range)

# merge all existing data sets for Or35a into a consensus model response
model_response_Or35a <- model_response(Or35a, plot = TRUE)

}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_door_database.R
\name{update_door_database}
\alias{update_door_database}
\title{update response matrix}
\usage{
update_door_database(receptor, permutation = TRUE, perm,
  response_matrix_nn = door_default_values("door_response_matrix_non_normalized"),
  response_matrix = door_default_values("door_response_matrix"),
  responseRange = door_default_values("door_response_range"),
  weightGlobNorm = door_default_values("door_global_normalization_weights"),
  select.MDValue = door_default_values("select.MDValue"), strict = TRUE,
  overlapValues = door_default_values("overlapValues"),
  door_excluded_data = door_default_values("door_excluded_data"),
  plot = FALSE)
}
\arguments{
\item{receptor}{character string, name of given odorant receptor.}

\item{permutation}{logical, if TRUE, the sequence is chosen from permutation,
if FALSE, sequence is chosen by the routine process.}

\item{perm}{a matrix with one sequence of study names per row, if empty, all 
possible permutations of study names will be provided.}

\item{response_matrix_nn}{data frame, response data that has not been 
globally normalized.}

\item{response_matrix}{data frame, globally normalized response data.}

\item{responseRange}{data frame, response range of studies.}

\item{weightGlobNorm}{data frame, weight matrix for global normalization.}

\item{select.MDValue}{the minimum mean distance between studies to perfom a 
merge (used if permutation == FALSE or if permutation == TRUE AND strict ==
TRUE)}

\item{strict}{logical, if TRUE merging a permutation will be stopped once a 
single merge has a mean distance above select.MDValue (only valid if 
permutation == TRUE)}

\item{overlapValues}{minimum overlap between studies to perfom a merge}

\item{door_excluded_data}{the data frame that contains the list of excluded
data sets.}

\item{plot}{logical}
}
\description{
update the globally \code{response matrix} and the unglobally normalized 
response matrix \code{door_response_matrix_non_normalized} by introducing new
consensus response data of given receptor.
}
\details{
The merging sequence could be arranged by the routine process (using 
\code{\link{model_response}} or taking the optimized sequence that is chosen
from permutations. The mean correlation between merged responses and each
original recording will be computed for each permutation, the optimozed
sequence is with the highest correlation.
}
\examples{
\dontrun{
# load data
library(DoOR.data)
load_door_data()
# update the entry "Or67b" of data "door_response_matrix" and
# "door_response_matrix_non_normalized" with permutations.
 update_door_database(receptor="Or67b", permutation = TRUE)
}
}
\seealso{
\code{\link{model_response}},\code{\link{model_response_seq}}
}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>

Shouwen Ma <\email{daniel.muench@uni-konstanz.de}>
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_response_seq.R
\name{model_response_seq}
\alias{model_response_seq}
\alias{modelRPSEQ}
\title{model_response_seq}
\usage{
model_response_seq(data, SEQ,
  overlapValues = door_default_values("overlapValues"),
  select.MDValue = door_default_values("select.MDValue"), strict = TRUE,
  plot = FALSE)
}
\arguments{
\item{data}{data frame, odorant response data, e.g. Or22a.}

\item{SEQ}{character vector, containing the names of studies indicating given
sequence for merging data.}

\item{overlapValues}{minimum overlap between studies to perfom a merge}

\item{select.MDValue}{the minimum mean distance between studies to perfom a
merge}

\item{strict}{logical, if TRUE merging a permutation will be stopped once a
single merge has a mean distance above select.MDValue}

\item{plot}{logical}
}
\description{
generates a model response and merge data in given sequence
}
\details{
# model_response_seq.R: #################

# merges studies in a given sequence (determined by the user or by exhaustive
enumeration and choosing the optimal sequence)


# input parameters: ####################


# data  : data frame, odorant response data for a given receptor, e.g. Or22a
# SEQ 	: character vector, contains the names of studies that measured this
receptor in a specific order (the merging sequence)

# output is a numeric vector: response values
}
\examples{
# load data
library(DoOR.data)
data(Or35a)
data(door_response_range)

# specify a sequence of merging
SEQ <- c("Hallem.2006.EN","Kreher.2008.EN","Hallem.2006.EN")

# perform the merging
selected_merg <- model_response_seq(Or35a, SEQ = SEQ, plot = TRUE)

}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load2list.R
\name{load2list}
\alias{load2list}
\title{load2list}
\usage{
load2list()
}
\value{
a list
}
\description{
returns all original DoOR response data as a list
}
\examples{
# load DoOR.data
library(DoOR.data)
load_door_data(nointeraction = TRUE)

# write the data into a list
lst <- load2list()
}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_new_data.R
\name{import_new_data}
\alias{import_new_data}
\alias{importNewData}
\title{Import new data into DoOR}
\usage{
import_new_data(file.name,
  dataFormat = door_default_values("door_data_format"),
  odor_data = door_default_values("odor"),
  weightGlobNorm = door_default_values("door_global_normalization_weights"),
  responseRange = door_default_values("door_response_range"),
  receptors = door_default_values("ORs"),
  ident = door_default_values("ident"), round = 3)
}
\arguments{
\item{file.name}{character string, the name of given file that contains 
response values of one or more odorant receptors, either a .csv or .txt 
file.}

\item{dataFormat}{data frame, a data frame does not contain any response 
value but odorant information.}

\item{odor_data}{data frame, contains the odorant information.}

\item{weightGlobNorm}{data matrix, indicates whether given receptor has been 
measured by given study.}

\item{responseRange}{data frame, contains the information about response 
range of each study and how many odors have been measured in each study.}

\item{receptors}{data frame, contains the receptor and OSN names and their 
expression.}

\item{ident}{the identifier used for matching, usually the InChIKey is used.}

\item{round}{the number of digits the imported values are rounded to.}
}
\description{
Import or update new data and update
\code{door_global_normalization_weights}, \code{door_response_range},
\code{odor}, \code{ORs} and receptor data frames.
}
\details{
\code{\link{import_new_data}} is used to import new data into database. If
the data contains a new receptor or ORN, then build a new data frame for this
receptor or ORN. If the data contains a receptor that is already present in 
database, then merge the imported data into old data frame. The information 
(e.g. response range, how many receptors and odors were measured from given 
study) will be integrated into data \code{door_response_range}, \code{odor}, 
\code{ORs} and \code{door_global_normalization_weights}. If an existing study
is imported, \code{\link{remove_study}} will be run first in order to perform
an update.
}
\examples{
\dontrun{
import new data named "odorantResponses_Orx.txt" into database and update the
support data.
library(DoOR.data)
import_new_data(file.name = "odorantResponses_Orx.csv")
}

}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>

Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_responses.R
\name{get_responses}
\alias{get_responses}
\alias{getResponses}
\title{Find receptor responses}
\usage{
get_responses(odorant,
  responseRange = door_default_values("door_response_range"),
  Or.list = load2list())
}
\arguments{
\item{odorant}{a single odor provided as InChIKey}

\item{responseRange}{data frame, response ranges of studies}

\item{Or.list}{a list contains reponse data of all available receptors. It
can be loaded using \code{\link{load2list}}.}
}
\description{
given a chemical, get original receptor responses from all studies in the
database.
}
\details{
output is a data frame containing response values of given odor across
receptors from all available studies.
}
\examples{
# load data
library(DoOR.data)
load_door_data(nointeraction = TRUE)

# get raw responses for odorant MLFHJEHSLIIPHL-UHFFFAOYSA-N
responses <- get_responses(odorant = 'MLFHJEHSLIIPHL-UHFFFAOYSA-N')

}
\author{
Daniel Münch \email{daniel.muench@uni-konstanz.de}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/back_project.R
\name{back_project}
\alias{back_project}
\title{back_project}
\usage{
back_project(template.data, responding.unit,
  response_matrix = door_default_values("door_response_matrix"))
}
\arguments{
\item{template.data}{data frame, the data that provides the scale to
transform to, containing InChIKeys in a column called "odorants" and
responses in a column called "responses"}

\item{responding.unit}{character, the name of the receptor/OSN/glomerulus
which responses should be transformed}

\item{response_matrix}{DoOR response matrix, the source data is picked from
here}
}
\value{
Output of back_project is a list containing a data frame with the
  back_projected data, the original data, the data used as a template and the
  data that resulted from fitting source and template (before rescaling to
  the template scale). additionaly the parameters of the linear fit between
  the source and template response scale is returned.
}
\description{
project the model response values back into your scale of interest (spike,
deltaF/F...)
}
\details{
The process of back projection is the following: \itemize{ \item 1. rescale
both data sets to [0,1], \item 2. find the best fitting model between
"bp.data" and "cons.data" (lowest MD value), \item 3. project the consensus
data onto the fitted curv, these are now our normalized, back projected
responses \item 4. rescale all responses to the scale of the original data
via a linear fit.  }
}
\examples{
# load some data sets
data(Or22a)
data(door_response_matrix)

# create example data we are going to back project
template.data <- data.frame(odorants = Or22a$InChIKey,
                            responses = Or22a$Hallem.2004.EN)

# run back_project and plot the results
bp <- back_project(template.data, "Or22a")

plot(bp$back_projected$original.data,
     bp$back_projected$back_projected.data,
     xlab = "DoOR consensus response",
     ylab = "back_projected data [spikes, Hallem.2004.EN]"
)

}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>

Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/door_default_values.R
\name{door_default_values}
\alias{door_default_values}
\alias{default.val}
\title{default values for DoOR functions}
\usage{
door_default_values(DoOR_default)
}
\arguments{
\item{DoOR_default}{a character string, indicating which argument is to be 
returned for DoOR functions.}
}
\description{
\code{door_default_values} is used to return default values for DoOR
functions.
}
\details{
There are six categories for default value. real number, integer, logical, 
NULL, character string and character vector.
}
\examples{
# extract DoOR default values
door_default_values(DoOR_default = "select.MD")

}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplot_al_map.R
\name{dplot_al_map}
\alias{dplot_al_map}
\alias{dplot_ALmap}
\title{dplot_al_map}
\usage{
dplot_al_map(InChIKey,
  response_matrix = door_default_values("door_response_matrix"),
  odor_data = door_default_values("odor"),
  door_mappings = door_default_values("door_mappings"),
  zero = door_default_values("zero"),
  tag = door_default_values("tag.ALmap"), main = "Name",
  scalebar = door_default_values("scalebar"),
  door_AL_map = door_default_values("door_AL_map"),
  colors = door_default_values("colors"), legend = TRUE, limits,
  base_size = 12)
}
\arguments{
\item{InChIKey}{InChIKey specifying the odorant to plot}

\item{response_matrix}{the input data (e.g. door_response_matrix or
door_response_matrix_non_normalized)}

\item{odor_data}{data frame, contains the odorant information.}

\item{door_mappings}{the data frame containing the mapping information}

\item{zero}{the odorant to set to zero (defaults to "SFR")}

\item{tag}{the labels to plot on top of the glomeruli (one of the following
\code{door_mappings} columns: "receptor", "sensillum", "ORN", "glomerulus"
or "co.receptor")}

\item{main}{the title, one column of \code{odor}, defaults to "Name"}

\item{scalebar}{whether or not to add a scalebar}

\item{door_AL_map}{a list containing the AL model}

\item{colors}{a vector containing 6 color values (2 for values below 0, 1 0
value and 3 steps between 0 and 1)}

\item{legend}{logical, plot a legend?}

\item{limits}{the limits for the color scale, if empty the range of the
response matrix is taken (after setting ``zero`` to 0)}

\item{base_size}{numeric, the base font size for the ggplot plot}
}
\value{
a ggplot2 object
}
\description{
Plot an antennal lobe map with color coded odorant responses.
}
\details{
Normalized, color coded odor responses across receptors are mapped
  onto a map of the \emph{Drosophila} antennal lobe. The antennal lobe map
  was a kind gift from Veit Grabe.
}
\examples{
# load data
library(DoOR.data)
library(DoOR.functions)

# map responses on antennal lobe scheme
dplot_al_map("MLFHJEHSLIIPHL-UHFFFAOYSA-N", scalebar = FALSE)

# change colors
dplot_al_map("MLFHJEHSLIIPHL-UHFFFAOYSA-N", tag = "Ors",
   color = c("magenta", "pink", "white", "yellow", "orange", "red"))

# pass some ggplot2 theming parameters
dplot_al_map(trans_id("123-92-2"), scalebar = FALSE) +
ggplot2::theme(legend.position  = "bottom",
     panel.background = ggplot2::element_rect(fill = "grey90", color = NA)) +
ggplot2::ggtitle("responses elicited by isopentyl acetate")

# export as pdf
\dontrun{
p <- dplot_al_map(trans_id("123-92-2"))
ggplot2::ggsave("AL.response.pdf", p, width = 6, height = 2, scale = 2)
}
}
\references{
Grabe, V., Strutz, A., Baschwitz, A., Hansson, B.S., Sachse, S.,
  2014. A digital in vivo 3D atlas of the antennal lobe of Drosophila
  melanogaster. J. Comp. Neurol. n/a–n/a. doi:10.1002/cne.23697
}
\seealso{
\link{get_normalized_responses}, \pkg{ggplot2}, \pkg{grid}
}
\author{
Daniel Münch \email{daniel.muench@uni-konstanz.de}

Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trans_id.R
\name{trans_id}
\alias{trans_id}
\alias{transID}
\title{trans_id}
\usage{
trans_id(x, from = "CAS", to = "InChIKey",
  odor_data = door_default_values("odor"))
}
\arguments{
\item{x}{character vector, one or many chemical identifiers}

\item{from}{character, the type of identifier to translate from (one of the
column names of ``odor``)}

\item{to}{character, the type of identifier to translate from (one of the
column names of ``odor``)}

\item{odor_data}{the data frame containing the odor information (defaults to
``odor``).}
}
\value{
character vector of translated chemical identifiers
}
\description{
Translate chemical identifiers from one to the other.
}
\examples{
# load data
library(DoOR.data)

# transform CAS to InChIKey
trans_id("123-92-2")

# transform Name to InChIKey
trans_id("isopentyl acetate", "Name")

# transform SMILE to InChIKey
trans_id("C(C(C)C)COC(=O)C", "SMILES", "Name")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_door_database.R
\name{create_door_database}
\alias{create_door_database}
\alias{CreateDatabase}
\title{Compose a Response Matrix of All Odor Receptors}
\usage{
create_door_database(tag = door_default_values("tag"),
  select.MDValue = door_default_values("select.MDValue"),
  overlapValues = door_default_values("overlapValues"), ...)
}
\arguments{
\item{tag}{character string, format for rownames, possibilities: "InChIKey",
CAS", "CID", "Name"}

\item{select.MDValue}{a numeric, threshold on the MD, this is used to reject
studies that do not align sufficiently well to the response model}

\item{overlapValues}{numeric, a criterion using to refuse a data set that
has not enough overlap value.}

\item{...}{pass more parameters to \code{\link{model_response}}}
}
\description{
computes the complete response model for all receptors in the database (calls
\code{\link{model_response}} for all receptors). Overwrites response_matrix,
door_response_matrix_non_normalized and door_excluded_data.
}
\examples{
\dontrun{
# load DoOR data
library(DoOR.data)
load_door_data()

# create a new consensus matrix
create_door_database()
}
}
\seealso{
\code{\link{model_response}}
}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>

Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_model.R
\name{calculate_model}
\alias{calculate_model}
\alias{calModel}
\title{' select the best model function}
\usage{
calculate_model(x, y, select.MD = door_default_values("select.MD"))
}
\arguments{
\item{x, y}{data vectors from study x and y (can contain NA)}

\item{select.MD}{logical, if TRUE, only the best model function (in terms of
MD) will be returned.}
}
\description{
\code{calculate_model} is used to return the best model function that
represent the relationship between responses from study x and y.
}
\details{
\code{calculate_model} chooses the best model function from following: linear,
exponential function, sigmoid, asymptotic model with x intercept, asympototic
model with y intercept and their inverse versions. (If your are interested in
these functions please check the sources at
https://github.com/Dahaniel/DoOR.functions)
}
\examples{
# load a data set
library(DoOR.data)
data(Or35a)

# pick 2 data sets for Or35a and rescale the data [0,1]
x <- door_norm(Or35a[,6])
y <- door_norm(Or35a[,9])
# run calculate_model
calM_xy <- calculate_model(x,y, select.MD = door_default_values("select.MD"))

}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>

Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplot_across_osns.R
\name{dplot_across_osns}
\alias{dplot_across_osns}
\alias{dplot_acrossOSNs}
\title{dplot_across_osns}
\usage{
dplot_across_osns(odorants,
  response_matrix = door_default_values("door_response_matrix"),
  odor_data = door_default_values("odor"),
  door_mappings = door_default_values("door_mappings"),
  zero = door_default_values("zero"), tag = "Name", sub, plot.type = 1,
  base_size = 12)
}
\arguments{
\item{odorants}{character vector, one or several InChIKeys}

\item{response_matrix}{DOOR response matrix, contains the data to plot}

\item{odor_data}{data frame, contains the odorant information.}

\item{door_mappings}{data frame, containing the mappings of response profiles
to morphological structures.}

\item{zero}{character, InChIKey of the odorant that should be set to 0 (e.g.
SFR)}

\item{tag}{character, the chemical identifier to use in the plot, one of
\code{colnames(odor)}}

\item{sub}{character vector, specify one or several classes of sensilla the
plot should be restricted to. One or several of: "ab" "ac", "at", "ai",
"pb", "sacI", "sacII"}

\item{plot.type}{interger, 1 or 2, defining the type of plot (1: facet_grid
of odorants * sensillae, 2: facet_wrap across OSNs)}

\item{base_size}{numeric, the base font size for the ggplot2 plot}
}
\value{
a ggplot2 object
}
\description{
plot the activation patterns of one or several odorants across OSNs
}
\details{
DoOR response profiles will be selected and ordered according to the
  OSNs they are related to. Several DoOR response profiles might exist for a
  given OSN (e.g. one for the OSN itself and one for the OSNs misexpressed
  receptor protein) but only one will be shown. Which DoOR profile is mapped
  to which OSN is controlled via the "code.OSN" column in \code{DoORmapings}.
}
\examples{
# load DoOR data & functions
library(DoOR.data)
library(DoOR.functions)

# pick example odorants by name ans transform their ID to InChIKey 
odorants <- trans_id(c("1-butanol", "isopentyl acetate", "carbon dioxide", "water"), 
from = "Name", to = "InChIKey")
 
# plot                                      
dplot_across_osns(odorants)
# plot only across ac and at sensillae
dplot_across_osns(odorants, sub = c("ac", "at"))
# plot across sensory neurons
dplot_across_osns(odorants, plot.type = 2)
}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplot_response_profile.R
\name{dplot_response_profile}
\alias{dplot_response_profile}
\alias{dplot_responseProfile}
\title{dplot_response_profile}
\usage{
dplot_response_profile(receptor,
  response_matrix = door_default_values("door_response_matrix"),
  odor_data = door_default_values("odor"), tag = door_default_values("tag"),
  colored = TRUE, colors = door_default_values("colors"), limits,
  zero = door_default_values("zero"),
  scalebar = door_default_values("scalebar"), base_size = 12)
}
\arguments{
\item{receptor}{character, receptor name, any of
colnames(door_response_matrix)}

\item{response_matrix}{a DoOR door_response_matrix}

\item{odor_data}{data frame, contains the odorant information.}

\item{tag}{character, chemical identifier for annotation}

\item{colored}{logical, color code the bars according to the response value?}

\item{colors}{character vector, a vector of 5 colors (2 for values < 0, 1 
value for 0 and 3 values > 0)}

\item{limits}{numeric of length 2, the limits for the colorscale and the x 
axis, global range of data will be used if empty}

\item{zero}{character, the odorant response that is set to 0, defaults to 
"SFR"}

\item{scalebar}{logical, add or suppress scalebars}

\item{base_size}{numeric, the base font size for the ggplot2 plot}
}
\value{
ggplot2 plot
}
\description{
create a barplot of a DoOR response profile
}
\examples{
# load data
library(DoOR.data)
data(door_response_matrix)

# plot with default parameters
dplot_response_profile("Or22a", door_response_matrix)

# plot wit odorant names
dplot_response_profile("Or22a", door_response_matrix, tag = "Name")
}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/door_norm.R
\name{door_norm}
\alias{door_norm}
\alias{DoORnorm}
\title{rescale the data values from 0 to 1}
\usage{
door_norm(x)
}
\arguments{
\item{x}{a numeric vector}
}
\description{
\code{door_norm} is used to normalize the data in values from 0 to 1.
}
\examples{
# create example data
x <- rnorm(10)

# run door_norm on it
door_norm(x)

}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>

Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
\keyword{math}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplot_compare_profiles.R
\name{dplot_compare_profiles}
\alias{dplot_compare_profiles}
\alias{dplot_compareProfiles}
\title{Compare two response profiles}
\usage{
dplot_compare_profiles(x, y, by.x, by.y, tag = "Name", base_size = 12)
}
\arguments{
\item{x}{the input data frame the first response profile will be taken from}

\item{y}{the input data frame the second response profile will be taken from
(x will be taken if y is missing)}

\item{by.x}{character string, specifying a column in x}

\item{by.y}{character string, specifying a column in y}

\item{tag}{character, the chemical identifier that will be used as odorant
label.}

\item{base_size}{numeric, the base font size for the ggplot2 plot}
}
\description{
Orderdered bar plots for two studies, allowing for an easy comparison of the
two studies / response profiles'.
}
\examples{
# load data
library(DoOR.data)
library(DoOR.functions)
data(Or22a)
data(door_response_range)
data(door_response_matrix)

# compare the Hallem and the Pelz data set for Or22a
dplot_compare_profiles(x = Or22a, y = Or22a,
                         by.x = "Hallem.2006.EN",
                         by.y = "Pelz.2006.AntEC50")

# comparedata from two different sensory neurons and add odorant labels 
dplot_compare_profiles(x = cbind(door_response_matrix, InChIKey =
rownames(door_response_matrix)), y = cbind(door_response_matrix, InChIKey =
rownames(door_response_matrix)), by.x = "Or22a", by.y = "Or10a")

}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/identify_sensillum.R
\name{identify_sensillum}
\alias{identify_sensillum}
\alias{identifySensillum}
\title{identify_sensillum}
\usage{
identify_sensillum(recording,
  response_matrix = door_default_values("door_response_matrix"),
  odor_data = door_default_values("odor"),
  door_mappings = door_default_values("door_mappings"), tag = "Name",
  min.cor = 0.9, nshow = 10, method = "cor", sub, base_size = 12,
  plot = TRUE, use = "everything")
}
\arguments{
\item{recording}{data frame, a data frame with the following columns 
"odorants" containing InChIKeys of the tested odorrant, and one column 
called "unit1" etc. for each unit, containing responses (or estimates) 
scaled between 0 and 1 (see examples)}

\item{response_matrix}{DoOR response matrix, the data to compair against}

\item{odor_data}{data frame, contains the odorant information.}

\item{door_mappings}{the data frame containing the mapping information}

\item{tag}{character, the chemical identifier to use in plots, one of 
\code{colnames(odor)}}

\item{min.cor}{numeric, a minimum correlation value, the function will check 
wether there is a higher correlation for all units within a single 
sensillum}

\item{nshow}{numeric, the number of plots to nshow, plot e.g. only the top 10
matches}

\item{method}{character, the method for similarity calculations: correlation 
("cor") or Euclidean distances ("dist")}

\item{sub}{character, if you know the class of sensillum you were recording 
from you can restrict the search to this subset here ("ab", "ac", "at", 
"pb", "sac")}

\item{base_size}{numeric, the base font size of the ggplot plots}

\item{plot}{logical, if TRUE returns the plot, else returns the data frame
with the correlations/distances}

\item{use}{character, the "use" option from the \code{cor} function, "all" 
returns NA when pairs are incomplete, "na.or.complete" only uses complete 
observations to calculate correlations; see \code{\link{cor}} for details}
}
\value{
eithe& Carolin G.(†27)r a plot (gtable) with responses sorted by
  highest correlations or lowest distances, or a data frame containing all
  calculated correlations or Euclidean distances
}
\description{
correlates the result from a SSR recording of several odorants against all 
DoOR response profiles
}
\examples{
# load data
library(DoOR.data)

# create an example recording
recording <- data.frame(
   odorants = c(trans_id(c("BEDN", "ETAS"), "Code"),
   trans_id("carbon dioxide", "Name")),
   unit1 = c(.9,.1,.1),
   unit2 = c(0, .1, 1)
)

# run the identification
identify_sensillum(recording)
identify_sensillum(recording, method = "dist", nshow = 5)

}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_dataset.R
\name{get_dataset}
\alias{get_dataset}
\alias{getDataset}
\title{getDataset}
\usage{
get_dataset(study, na.rm = FALSE)
}
\arguments{
\item{study}{character, the name of the study you want to aggregate the dta
from}

\item{na.rm}{logical, whether or not to exclude odorants that were not
measured in the study}
}
\value{
returns a data frame containing all the odorant responses measured in
  \code{study}
}
\description{
aggregates original data from a given study
}
\examples{
# load data
library(DoOR.data)
load_door_data(nointeraction = TRUE)

# get all recordings from the Hallem.2004.EN data set
get_dataset("Hallem.2004.EN", na.rm = TRUE)
}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_points.R
\name{project_points}
\alias{project_points}
\alias{projectPoints}
\title{project_points}
\usage{
project_points(x, y, xylim, best.model, plot = door_default_values("plot"),
  points_cex = door_default_values("points.cex"),
  title = door_default_values("title"), ...)
}
\arguments{
\item{x, y}{numeric vectors of data values, coordinate vectors of points to
plot, the coordinates can contain \code{NA} values.}

\item{xylim}{numeric vectors, x, y limits of the plot.}

\item{best.model}{a list, containing the parameters, function, inverse
function, Leibniz's notation for distance calculation and MD value. if
missing, the best model will be generated automatically.}

\item{plot}{logical, If \code{FALSE}, plotting is suppressed. Default is
\code{FALSE}.}

\item{points_cex}{a numerical value, giving the magnification level of
symbols relative to the default size.}

\item{title}{logical, If \code{TRUE}, title is shown. Default is
\code{FALSE}.}

\item{\dots}{further graphical parameters}
}
\description{
projects data points onto the curve defined by the model function
}
\details{
For internal use in the merging process (see also
\code{\link{model_response}}). The model function is choosen by
\code{\link{calculate_model}}. \code{\link{project_points}} then projects the
data points from the datasets to be merged onto the curve defined by the
model function. It computes the closest distance from a data point to a point
on the curve by numerical optimisation.

A list with two data frames "double.observations" and "single.observations"
is returned, which give the coordinates of double observations (defined as
(x,y)) and coordinates of single observation (defined as (x,NA) or (NA,y)).
Both data frames contain seven columns: "ID" indicating the original position
of data x and y, "x", "y" indicating the coordinate of observation, "X", "Y"
indicating the coordinate of projected point on the function, "distance"
indicating the distances between \code{(xmin, f(xmin))} and all points on the
functional line, "NDR" indicating the normalized distances across all the
distance values.
}
\examples{
# load data
library(DoOR.data)
data(Or23a)

# normalize two example data sets
x <- door_norm(Or23a[,'Hallem.2004.EN'])
y <- door_norm(Or23a[,'Hallem.2006.EN'])

# find the best fitting function and project the remaining points (measured
# only in one of the data sets) onto the fit.
project_points(x = x, y = y, plot = TRUE)

}
\seealso{
\code{\link{calculate_model}}, \code{\link{optimize}},
  \code{\link{integrate}}
}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_receptor.R
\name{map_receptor}
\alias{map_receptor}
\alias{mapReceptor}
\title{map_receptor}
\usage{
map_receptor(data,
  response_matrix = door_default_values("door_response_matrix"), sub,
  threshold.p, threshold.cor, nshow)
}
\arguments{
\item{data}{data frame, containing two columns, one called "odorants" and one
"responses" providing InChIKeys and odorant responses respectively.}

\item{response_matrix}{output is a numeric vector that contains the Pearson 
Correlation Coefficient between given data and selected consensus data in}

\item{sub}{character, a subset of responding units returned response matrix}

\item{threshold.p}{numeric, a p-value threshold, only correlations below will
be returned}

\item{threshold.cor}{numeric, a correlation-coefficient threshold, only
correlations above will be returned}

\item{nshow}{numeric, if defined, only this number of results will be}
}
\description{
Identifying the source of unknown response data by correlating it agains all 
DoOR responding units.
}
\examples{
# load data
load_door_data(nointeraction = TRUE)

# pick example data
data <- data.frame(odorants  = Or22a$InChIKey,
                   responses = Or22a$Hallem.2004.EN)
data <- na.omit(data)

# find the corresponding receptor / responding unit
map_receptor(data = data)
}
\author{
Shouwen Ma <\email{shouwen.ma@uni-konstanz.de}>

Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse.R
\name{sparse}
\alias{sparse}
\title{Calculate the sparseness of a vector.}
\usage{
sparse(x, method = "ltk")
}
\arguments{
\item{x}{numerical input vector}

\item{method}{type of sparseness measure, either 'ltk' for lifetime kurtosis
or 'lts' for lifetime sparseness (see references).}
}
\description{
Sparseness can be calculated as lifetime kurtosis (LTK, Willmore and
Tolhurst, 2001) or as lifetime sparseness (LTS, Bhandawat et al., 2007).
}
\details{
LTS scales between \[0,1\] while LTK is not restricted. LTS only
  takes positive values.
}
\references{
Bhandawat, V., Olsen, S.R., Gouwens, N.W., Schlief, M.L., Wilson,
  R.I., 2007. Sensory processing in the Drosophila antennal lobe increases
  reliability and separability of ensemble odor representations. Nature
  neuroscience 10, 1474–82. doi:10.1038/nn1976

Willmore, B., Tolhurst, D.J., 2001. Characterizing the sparseness
  of neural codes. Network 12, 255–270. doi:10.1080/net.12.3.255.270
}
\author{
Daniel Münch <\email{daniel.muench@uni-konstanz.de}>
}
\keyword{kurtosis}
\keyword{sparseness,}
