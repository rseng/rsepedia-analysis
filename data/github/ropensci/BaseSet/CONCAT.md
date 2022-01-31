
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/ropensci/BaseSet.svg?branch=master)](https://travis-ci.org/ropensci/BaseSet)
[![R build
status](https://github.com/ropensci/BaseSet/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/BaseSet/actions)
[![Coverage
status](https://codecov.io/gh/ropensci/BaseSet/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/BaseSet?branch=master)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![rOpenSci](https://badges.ropensci.org/359_status.svg)](https://github.com/ropensci/software-review/issues/359)
[![CRAN
status](https://www.r-pkg.org/badges/version/BaseSet)](https://CRAN.R-project.org/package=BaseSet)
<!-- badges: end -->

# BaseSet

The goal of BaseSet is to facilitate working with sets in an efficient
way. The package implements methods to work on sets, doing intersection,
union, complementary, power sets, cartesian product and other set
operations in a tidy way.

The package supports
[classical](https://en.wikipedia.org/wiki/Set_(mathematics)) and
[fuzzy](https://en.wikipedia.org/wiki/Fuzzy_set) sets. Fuzzy sets are
similar to classical sets but there is some vagueness on the
relationship between the element and the set.

It also allows to import from several formats used in the life science
world. Like the
[GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
and the
[GAF](http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/)
or the [OBO format](http://www.obofoundry.org/) file for ontologies.

You can save information about the elements, sets and their relationship
on the object itself. For instance origin of the set, categorical or
numeric data associated with sets…

Watch BaseSet working on the [examples](#Examples) below and in the
vignettes. You can also find [related packages](#Related-packages) and
the differences with BaseSet. If you have some questions or bugs [open
an issue](https://github.com/ropensci/BaseSet/issues) (remember the
[Code of Conduct](#Code-of-Conduct))

# Installation

The package depends on some packages from Bioconductor. In order to
install some of its dependencies you’ll need first to install
`{BiocManager}`:

``` r
if (!require("BiocManager")) {
  install.packages("BiocManager")
}
```

You can install the latest version of BaseSet from
[Github](https://github.com/ropensci/BaseSet) with:

``` r
BiocManager::install("ropensci/BaseSet", 
                     dependencies = TRUE, build_vignettes = TRUE, force = TRUE)
```

# Examples

## Sets

We can create a set like this:

``` r
sets <- list(A = letters[1:5], B = c("a", "f"))
sets_analysis <- tidySet(sets)
sets_analysis
#>   elements sets fuzzy
#> 1        a    A     1
#> 2        b    A     1
#> 3        c    A     1
#> 4        d    A     1
#> 5        e    A     1
#> 6        a    B     1
#> 7        f    B     1
```

Perform typical operations like union, intersection. You can name the
resulting set or let the default name:

``` r
union(sets_analysis, sets = c("A", "B")) 
#>   elements sets fuzzy
#> 1        a  A∪B     1
#> 2        b  A∪B     1
#> 3        c  A∪B     1
#> 4        d  A∪B     1
#> 5        e  A∪B     1
#> 6        f  A∪B     1
# Or we can give a name to the new set
union(sets_analysis, sets = c("A", "B"), name = "D")
#>   elements sets fuzzy
#> 1        a    D     1
#> 2        b    D     1
#> 3        c    D     1
#> 4        d    D     1
#> 5        e    D     1
#> 6        f    D     1
# Or the intersection
intersection(sets_analysis, sets = c("A", "B"))
#>   elements sets fuzzy
#> 1        a  A∩B     1
# Keeping the other sets:
intersection(sets_analysis, sets = c("A", "B"), name = "D", keep = TRUE) 
#>   elements sets fuzzy
#> 1        a    A     1
#> 2        b    A     1
#> 3        c    A     1
#> 4        d    A     1
#> 5        e    A     1
#> 6        a    B     1
#> 7        f    B     1
#> 8        a    D     1
```

And compute size of sets among other things:

``` r
set_size(sets_analysis)
#>   sets size probability
#> 1    A    5           1
#> 2    B    2           1
```

The elements in one set not present in other:

``` r
subtract(sets_analysis, set_in = "A", not_in = "B", keep = FALSE)
#>   elements sets fuzzy
#> 1        b  A∖B     1
#> 2        c  A∖B     1
#> 3        d  A∖B     1
#> 4        e  A∖B     1
```

Or any other verb from
[dplyr](https://cran.r-project.org/package=dplyr). We can add columns,
filter, remove them and add information about the sets:

``` r
library("magrittr")
#> 
#> Attaching package: 'magrittr'
#> The following object is masked from 'package:BaseSet':
#> 
#>     subtract
set.seed(4673) # To make it reproducible in your machine
sets_enriched <- sets_analysis %>% 
  mutate(Keep = sample(c(TRUE, FALSE), 7, replace = TRUE)) %>% 
  filter(Keep == TRUE) %>% 
  select(-Keep) %>% 
  activate("sets") %>% 
  mutate(sets_origin = c("Reactome", "KEGG"))
sets_enriched
#>   elements sets fuzzy sets_origin
#> 1        a    A     1    Reactome
#> 2        b    A     1    Reactome
#> 3        c    A     1    Reactome
#> 4        d    A     1    Reactome
#> 5        e    A     1    Reactome
#> 6        f    B     1        KEGG

# Activating sets makes the verb affect only them:
elements(sets_enriched)
#>   elements
#> 1        a
#> 2        b
#> 3        c
#> 4        d
#> 5        e
#> 6        f
relations(sets_enriched)
#>   elements sets fuzzy
#> 1        a    A     1
#> 2        b    A     1
#> 3        c    A     1
#> 4        d    A     1
#> 5        e    A     1
#> 6        f    B     1
sets(sets_enriched)
#>   sets sets_origin
#> 1    A    Reactome
#> 6    B        KEGG
```

## Fuzzy sets

In [fuzzy sets](https://en.wikipedia.org/wiki/Fuzzy_set) the elements
are vaguely related to the set by a numeric value usually between 0 and
1. This implies that the association is not guaranteed.

``` r
relations <- data.frame(sets = c(rep("A", 5), "B", "B"), 
                        elements = c("a", "b", "c", "d", "e", "a", "f"),
                        fuzzy = runif(7))
fuzzy_set <- tidySet(relations)
fuzzy_set
#>   elements sets     fuzzy
#> 1        a    A 0.1837246
#> 2        b    A 0.4567009
#> 3        c    A 0.8152075
#> 4        d    A 0.5800610
#> 5        e    A 0.5724973
#> 6        a    B 0.9381182
#> 7        f    B 0.9460158
```

The equivalent operations performed on classical sets are possible with
fuzzy sets:

``` r
union(fuzzy_set, sets = c("A", "B")) 
#>   elements sets     fuzzy
#> 1        a  A∪B 0.9381182
#> 2        b  A∪B 0.4567009
#> 3        c  A∪B 0.8152075
#> 4        d  A∪B 0.5800610
#> 5        e  A∪B 0.5724973
#> 6        f  A∪B 0.9460158
# Or we can give a name to the new set
union(fuzzy_set, sets = c("A", "B"), name = "D")
#>   elements sets     fuzzy
#> 1        a    D 0.9381182
#> 2        b    D 0.4567009
#> 3        c    D 0.8152075
#> 4        d    D 0.5800610
#> 5        e    D 0.5724973
#> 6        f    D 0.9460158
# Or the intersection
intersection(fuzzy_set, sets = c("A", "B"))
#>   elements sets     fuzzy
#> 1        a  A∩B 0.1837246
# Keeping the other sets:
intersection(fuzzy_set, sets = c("A", "B"), name = "D", keep = TRUE) 
#>   elements sets     fuzzy
#> 1        a    A 0.1837246
#> 2        b    A 0.4567009
#> 3        c    A 0.8152075
#> 4        d    A 0.5800610
#> 5        e    A 0.5724973
#> 6        a    B 0.9381182
#> 7        f    B 0.9460158
#> 8        a    D 0.1837246
```

Assuming that the fuzzy value is a probability, we can calculate which
is the probability of having several elements:

``` r
# A set could be empty
set_size(fuzzy_set)
#>   sets size probability
#> 1    A    0 0.014712455
#> 2    A    1 0.120607154
#> 3    A    2 0.318386944
#> 4    A    3 0.357078627
#> 5    A    4 0.166499731
#> 6    A    5 0.022715089
#> 7    B    0 0.003340637
#> 8    B    1 0.109184679
#> 9    B    2 0.887474684
# The more probable size of the sets:
set_size(fuzzy_set) %>% 
  group_by(sets) %>% 
  filter(probability == max(probability))
#> # A tibble: 2 x 3
#> # Groups:   sets [2]
#>   sets   size probability
#>   <chr> <dbl>       <dbl>
#> 1 A         3       0.357
#> 2 B         2       0.887
# Probability of belonging to several sets:
element_size(fuzzy_set)
#>    elements size probability
#> 1         a    0  0.05051256
#> 2         a    1  0.77713204
#> 3         a    2  0.17235540
#> 4         b    0  0.54329910
#> 5         b    1  0.45670090
#> 6         c    0  0.18479253
#> 7         c    1  0.81520747
#> 8         d    0  0.41993900
#> 9         d    1  0.58006100
#> 10        e    0  0.42750268
#> 11        e    1  0.57249732
#> 12        f    0  0.05398419
#> 13        f    1  0.94601581
```

With fuzzy sets we can filter at certain levels (called alpha cut):

``` r
fuzzy_set %>% 
  filter(fuzzy > 0.5) %>% 
  activate("sets") %>% 
  mutate(sets_origin = c("Reactome", "KEGG"))
#>   elements sets     fuzzy sets_origin
#> 1        c    A 0.8152075    Reactome
#> 2        d    A 0.5800610    Reactome
#> 3        e    A 0.5724973    Reactome
#> 4        a    B 0.9381182        KEGG
#> 5        f    B 0.9460158        KEGG
```

# Related packages

There are several other packages related to sets, which partially
overlap with BaseSet functionality:

-   [sets](https://CRAN.R-project.org/package=sets)  
    Implements a more generalized approach, that can store functions or
    lists as an element of a set (while BaseSet only allows to store a
    character or factor), but it is harder to operate in a tidy/long
    way. Also the operations of intersection and union need to happen
    between two different objects, while a single TidySet object (the
    class implemented in BaseSet) can store one or thousands of sets.

-   [`{GSEABase}`](https://bioconductor.org/packages/GSEABase/)  
    Implements a class to store sets and related information, but it
    doesn’t allow to store fuzzy sets and it is also quite slow as it
    creates several classes for annotating each set.

-   [`{BiocSet}`](https://bioconductor.org/packages/BiocSet/)  
    Implements a tidy class for sets but does not handle fuzzy sets. It
    also has less functionality to operate with sets, like power sets
    and cartesian product. BiocSet was influenced by the development of
    this package.

-   [`{hierarchicalSets}`](https://CRAN.R-project.org/package=hierarchicalSets)  
    This package is focused on clustering of sets that are inside other
    sets and visualizations. However, BaseSet is focused on storing and
    manipulate sets including hierarchical sets.

-   [`{set6}`](https://cran.r-project.org/package=set6) This package
    implements different classes for different type of sets including
    fuzzy sets, conditional sets. However, it doesn’t handle information
    associated to elements, sets or relationship.

# Why this package?

On bioinformatics when looking for the impact of an experiment
enrichment methods are applied. This involves obtaining several sets of
genes from several resources and methods. Usually these curated sets of
genes are taken at face value. However, there are several resources of
sets and they [do not agree between
them](https://doi.org/10.1186/1471-2105-14-112), regardless they are
used without considering any uncertainty on sets composition.

Fuzzy theory has long studied sets whose elements have degrees of
membership and/or uncertainty. Therefore one way to improve the methods
involve using fuzzy methods and logic on this field. As I couldn’t find
any package that provided methods for this I set on creating it (after
trying to [expand](https://github.com/llrs/GSEAdv) the existing one I
knew).

This package is intended to be easy to use for someone who is working
with collections of sets but flexible about the methods and logic it can
use. To be consistent, the standard fuzzy logic is the default but it
might not be the right one for your data. Consider changing the defaults
to match with the framework the data was obtained with.

# Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
# BaseSet (development version)

# BaseSet 0.0.17

* Fixing test when missing a package. 
* Adding copyright holder.

# BaseSet 0.0.16

* Update Code of Conduct to rOpenSci template
* Fix NOTE about LazyData

# BaseSet 0.0.15

* Upgrade R version requirements
* Fix some links
* Make sure that vignettes run when the Bioconductor packages are available

# BaseSet 0.0.14

* Remove unused dependency to BiocStyle
* Fix unicode characters for windows

# BaseSet 0.0.12

# BaseSet 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
This release aims to remove the errors on CRAN checks due to bad handling of suggested packages used on tests.

## Test environments
* local R installation, R 4.0.1
* win-builder (devel, release)
* Github Actions: windows (release, oldrel), macOS (release, oldrel), ubuntu (release oldrel)

## R CMD check results

0 errors | 0 warnings | 0 note

---
output: github_document
editor_options: 
  chunk_output_type: console
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/ropensci/BaseSet.svg?branch=master)](https://travis-ci.org/ropensci/BaseSet)
[![R build status](https://github.com/ropensci/BaseSet/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/BaseSet/actions)
[![Coverage status](https://codecov.io/gh/ropensci/BaseSet/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/BaseSet?branch=master)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![rOpenSci](https://badges.ropensci.org/359_status.svg)](https://github.com/ropensci/software-review/issues/359)
[![CRAN status](https://www.r-pkg.org/badges/version/BaseSet)](https://CRAN.R-project.org/package=BaseSet)
<!-- badges: end -->

# BaseSet

The goal of BaseSet is to facilitate working with sets in an efficient way. 
The package implements methods to work on sets, doing intersection, union, complementary, power sets, cartesian product and other set operations in a tidy way. 


The package supports [classical](https://en.wikipedia.org/wiki/Set_(mathematics)) and [fuzzy](https://en.wikipedia.org/wiki/Fuzzy_set) sets. 
Fuzzy sets are similar to classical sets but there is some vagueness on the relationship between the element and the set. 


It also allows to import from several formats used in the life science world. 
Like the [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) and the [GAF](http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/) or the [OBO format](http://www.obofoundry.org/) file for ontologies.

You can save information about the elements, sets and their relationship on the object itself. 
For instance origin of the set, categorical or numeric data associated with sets...

Watch BaseSet working on the [examples](#Examples) below and in the vignettes. 
You can also find [related packages](#Related-packages) and the differences with BaseSet. 
If you have some questions or bugs [open an issue](https://github.com/ropensci/BaseSet/issues) (remember the [Code of Conduct](#Code-of-Conduct))

# Installation

The package depends on some packages from Bioconductor. In order to install some of its dependencies you'll need first to install `{BiocManager}`:

```{r dep, eval = FALSE}
if (!require("BiocManager")) {
  install.packages("BiocManager")
}
```

You can install the latest version of BaseSet from [Github](https://github.com/ropensci/BaseSet) with:

```{r eval=FALSE}
BiocManager::install("ropensci/BaseSet", 
                     dependencies = TRUE, build_vignettes = TRUE, force = TRUE)
```

 
# Examples {#Examples}

```{r include=FALSE}
library("BaseSet")
```

## Sets

We can create a set like this:

```{r TidySet}
sets <- list(A = letters[1:5], B = c("a", "f"))
sets_analysis <- tidySet(sets)
sets_analysis
```

Perform typical operations like union, intersection. You can name the resulting set or let the default name:

```{r union-intersection}
union(sets_analysis, sets = c("A", "B")) 
# Or we can give a name to the new set
union(sets_analysis, sets = c("A", "B"), name = "D")
# Or the intersection
intersection(sets_analysis, sets = c("A", "B"))
# Keeping the other sets:
intersection(sets_analysis, sets = c("A", "B"), name = "D", keep = TRUE) 
```

And compute size of sets among other things:

```{r set_size}
set_size(sets_analysis)
```

The elements in one set not present in other:

```{r subraction}
subtract(sets_analysis, set_in = "A", not_in = "B", keep = FALSE)
```

Or any other verb from [dplyr](https://cran.r-project.org/package=dplyr). We can add columns, filter, remove them and add information about the sets:

```{r dplyr}
library("magrittr")
set.seed(4673) # To make it reproducible in your machine
sets_enriched <- sets_analysis %>% 
  mutate(Keep = sample(c(TRUE, FALSE), 7, replace = TRUE)) %>% 
  filter(Keep == TRUE) %>% 
  select(-Keep) %>% 
  activate("sets") %>% 
  mutate(sets_origin = c("Reactome", "KEGG"))
sets_enriched

# Activating sets makes the verb affect only them:
elements(sets_enriched)
relations(sets_enriched)
sets(sets_enriched)
```

## Fuzzy sets

In [fuzzy sets](https://en.wikipedia.org/wiki/Fuzzy_set) the elements are vaguely related to the set by a numeric value usually between 0 and 1.
This implies that the association is not guaranteed.

```{r fuzzy}
relations <- data.frame(sets = c(rep("A", 5), "B", "B"), 
                        elements = c("a", "b", "c", "d", "e", "a", "f"),
                        fuzzy = runif(7))
fuzzy_set <- tidySet(relations)
fuzzy_set
```

The equivalent operations performed on classical sets are possible with fuzzy sets:

```{r fuzzy-operations}
union(fuzzy_set, sets = c("A", "B")) 
# Or we can give a name to the new set
union(fuzzy_set, sets = c("A", "B"), name = "D")
# Or the intersection
intersection(fuzzy_set, sets = c("A", "B"))
# Keeping the other sets:
intersection(fuzzy_set, sets = c("A", "B"), name = "D", keep = TRUE) 
```

Assuming that the fuzzy value is a probability, we can calculate which is the probability of having several elements:

```{r prob}
# A set could be empty
set_size(fuzzy_set)
# The more probable size of the sets:
set_size(fuzzy_set) %>% 
  group_by(sets) %>% 
  filter(probability == max(probability))
# Probability of belonging to several sets:
element_size(fuzzy_set)
```

With fuzzy sets we can filter at certain levels (called alpha cut):

```{r alphaCut}
fuzzy_set %>% 
  filter(fuzzy > 0.5) %>% 
  activate("sets") %>% 
  mutate(sets_origin = c("Reactome", "KEGG"))
```

# Related packages {#related}

There are several other packages related to sets, which partially overlap with BaseSet functionality:

 - [sets]( https://CRAN.R-project.org/package=sets)  
 Implements a more generalized approach, that can store functions or lists as an element of a set (while BaseSet only allows to store a character or factor), but it is harder to operate in a tidy/long way. Also the operations of intersection and union need to happen between two different objects, while a single TidySet object (the class implemented in BaseSet) can store one or thousands of sets.

 - [`{GSEABase}`](https://bioconductor.org/packages/GSEABase/)  
 Implements a class to store sets and related information, but it doesn't allow to store fuzzy sets and it is also quite slow as it creates several classes for annotating each set. 
  
 - [`{BiocSet}`](https://bioconductor.org/packages/BiocSet/)  
 Implements a tidy class for sets but does not handle fuzzy sets. It also has less functionality to operate with sets, like power sets and cartesian product. BiocSet was influenced by the development of this package. 

 - [`{hierarchicalSets}`](https://CRAN.R-project.org/package=hierarchicalSets)  
 This package is focused on clustering of sets that are inside other sets and visualizations. However, BaseSet is focused on storing and manipulate sets including hierarchical sets.
 
 - [`{set6}`](https://cran.r-project.org/package=set6)
 This package implements different classes for different type of sets including fuzzy sets, conditional sets. However, it doesn't handle information associated to elements, sets or relationship. 
 
# Why this package? {#why}

On bioinformatics when looking for the impact of an experiment enrichment methods are applied.
This involves obtaining several sets of genes from several resources and methods.
Usually these curated sets of genes are taken at face value. 
However, there are several resources of sets and they [do not agree between them](https://doi.org/10.1186/1471-2105-14-112), regardless they are used without considering any uncertainty on sets composition. 


Fuzzy theory has long studied sets whose elements have degrees of membership and/or uncertainty. 
Therefore one way to improve the methods involve using fuzzy methods and logic on this field. 
As I couldn't find any package that provided methods for this I set on creating it (after trying to [expand](https://github.com/llrs/GSEAdv) the existing one I knew).

This package is intended to be easy to use for someone who is working with collections of sets but flexible about the methods and logic it can use. 
To be consistent, the standard fuzzy logic is the default but it might not be the right one for your data. 
Consider changing the defaults to match with the framework the data was obtained with. 

# Code of Conduct {#CoC}

Please note that this package is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.
---
title: "Fuzzy sets"
abstract: >
  Describes the fuzzy sets, interpretation and how to work with them.
date: "`r format(Sys.time(), '%Y %b %d')`"
output:
  html_document:
    fig_caption: true
    code_folding: show
    self_contained: yes
    toc_float:
      collapsed: true
      toc_depth: 3
author:
- name: Lluís Revilla
  affiliation: 
    - August Pi i Sunyer Biomedical Research Institute (IDIBAPS); Liver Unit, Hospital Clinic
  email: lluis.revilla@gmail.com
vignette: >
  %\VignetteIndexEntry{Fuzzy sets}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
  %\DeclareUnicodeCharacter{2229}{$\cap$}
  %\DeclareUnicodeCharacter{222A}{$\cup$}
editor_options: 
  chunk_output_type: console
---

```{r setup, message=FALSE, warning=FALSE, include=FALSE}
knitr::opts_knit$set(root.dir = ".")
knitr::opts_chunk$set(collapse = TRUE, 
                      warning = TRUE,
                      comment = "#>")
library("BaseSet")
library("dplyr")
```


# Getting started

This vignettes supposes that you already read the "About BaseSet" vignette. 
This vignette explains what are the fuzzy sets and how to use them. 
As all methods for "normal" sets are available for fuzzy sets this vignette focuses on how to create, use them.

# What are fuzzy sets and why/when use them ?

Fuzzy sets are generalizations of classical sets where there is some vagueness, one doesn't know for sure something on this relationship between the set and the element.
This vagueness can be on the assignment to a set and/or on the membership on the set.
One way these vagueness arise is when classifying a continuous scale on a categorical scale, for example: when the temperature is 20ºC is it hot or not? If the temperature drops to 15ºC is it warm?
When does the switch happen from warm to hot?

In fuzzy set theories the step from a continuous scale to a categorical scale is performed by the [membership function](https://en.wikipedia.org/wiki/Membership_function_(mathematics)) and is called [fuzzification](https://en.wikipedia.org/wiki/Fuzzy_logic#Fuzzification). 

When there is a degree of membership and uncertainty on the membership it is considered a [type-2 fuzzy set](https://en.wikipedia.org/wiki/Type-2_fuzzy_sets_and_systems).
We can understand it with using as example a paddle ball, it can be used for tennis or paddle (membership), but until we don't test it bouncing on the court we won't be sure (assignment) if it is a paddle ball or a tennis ball. 
We could think about a ball equally good for tennis and paddle (membership) but which two people thought it is for tennis and the other for paddle. 

These voting/rating system is also a common scenario where fuzzy sets arise.
When several graders/people need to agree but compromise on a middle ground. 
One modern example of this is ratting apps where several people vote between 1 and 5 an app and the displayed value takes into consideration all the votes. 

As you can see when one does have some vagueness or uncertainty then fuzzy logic is a good choice.
There have been developed several logic and methods. 

# Creating a fuzzy set

To create a fuzzy set you need to have a column named "fuzzy" if you create it 
from a `data.frame` or have a named numeric vector if you create it from a `list`.
These values are restricted to a numeric value between 0 and 1. 
The value indicates the strength, membership, truth value (or probability) of the relationship between 
the element and the set. 

```{r fuzzy}
set.seed(4567) # To be able to have exact replicates
relations <- data.frame(sets = c(rep("A", 5), "B", "C"),
                          elements = c(letters[seq_len(6)], letters[6]),
                          fuzzy = runif(7))
fuzzy_set <- tidySet(relations)
```


# Working with fuzzy sets

We can work with fuzzy sets as we do with normal sets. 
But if you remember that at the end of the previous vignette we used an Important column, now it is already included as fuzzy. 
This allows us to use this information for union and intersection methods and other operations:

## Union

You can make a union of two sets present on the same object.

```{r union}
BaseSet::union(fuzzy_set, sets = c("A", "B"))
BaseSet::union(fuzzy_set, sets = c("A", "B"), name = "D")
```

We get a new set with all the elements on both sets.
There isn't an element in both sets A and B, so the fuzzy values here do not change.

If we wanted to use other logic we can with provide it with the FUN argument. 

```{r union_logic}

BaseSet::union(fuzzy_set, sets = c("A", "B"), FUN = function(x){sqrt(sum(x))})
```

There are several logic see for instance `?sets::fuzzy_logic`. 
You should pick the operators that fit on the framework of your data. 
Make sure that the defaults arguments of logic apply to your data obtaining process. 

## Intersection

However if we do the intersection between B and C we can see some changes on the fuzzy value:

```{r intersection}
intersection(fuzzy_set, sets = c("B", "C"), keep = FALSE)
intersection(fuzzy_set, sets = c("B", "C"), keep = FALSE, FUN = "mean")
```

Different logic on the `union()`, `intersection()`, `complement()`, and `cardinality()`, we get different fuzzy values. 
Depending on the nature of our fuzziness and the intended goal we might apply one or other rules.

## Complement

We can look for the complement of one or several sets:

```{r complement}
complement_set(fuzzy_set, sets = "A", keep = FALSE)
```

Note that the values of the complement are `1-fuzzy` but can be changed:

```{r complement_previous}
filter(fuzzy_set, sets == "A")
complement_set(fuzzy_set, sets = "A", keep = FALSE, FUN = function(x){1-x^2})
```

## Subtract

This is the equivalent of `setdiff`, but clearer:

```{r subtract}
subtract(fuzzy_set, set_in = "A", not_in = "B", keep = FALSE, name = "A-B")
# Or the opposite B-A, but using the default name:
subtract(fuzzy_set, set_in = "B", not_in = "A", keep = FALSE)
```

Note that here there is also a subtraction of the fuzzy value.

# Sizes

If we consider the fuzzy values as probabilities then the size of a set is not fixed.
To calculate the size of a given set we have `set_size()`:

```{r set_size}
set_size(fuzzy_set)
```

Or an element can be in 0 sets:

```{r element_size}
element_size(fuzzy_set)
```

In this example we can see that it is more probable that the element "a" is not 
present than the element "f" being present in one set. 

# Interpretation

Sometimes it can be a bit hard to understand what do the fuzzy sets mean on your analysis.
To better understand let's dive a bit in the interpretation with an example:

Imagine you have your experiment where you collected data from a sample of cells for each cell (our elements). 
Then you used some program to classify which type of cell it is (alpha, beta, delta, endothelial), this are our sets. 
The software returns a probability for each type it has: the higher, the more confident it is of the assignment:

```{r cells_0}
sc_classification <- data.frame(
  elements = c("D2ex_1", "D2ex_10", "D2ex_11", "D2ex_12", "D2ex_13", "D2ex_14", 
               "D2ex_15", "D2ex_16", "D2ex_17", "D2ex_18", "D2ex_1", "D2ex_10", 
               "D2ex_11", "D2ex_12", "D2ex_13", "D2ex_14", "D2ex_15", "D2ex_16",
               "D2ex_17", "D2ex_18", "D2ex_1", "D2ex_10", "D2ex_11", "D2ex_12", 
               "D2ex_13", "D2ex_14", "D2ex_15", "D2ex_16", "D2ex_17", "D2ex_18", 
               "D2ex_1", "D2ex_10", "D2ex_11", "D2ex_12", "D2ex_13", "D2ex_14", 
               "D2ex_15", "D2ex_16", "D2ex_17", "D2ex_18"), 
  sets = c("alpha", "alpha", "alpha", "alpha", "alpha", "alpha", "alpha", 
           "alpha", "alpha", "alpha", "endothel", "endothel", "endothel", 
           "endothel", "endothel", "endothel", "endothel", "endothel", 
           "endothel", "endothel", "delta", "delta", "delta", "delta", "delta", 
           "delta", "delta", "delta", "delta", "delta", "beta", "beta", "beta", 
           "beta", "beta", "beta", "beta", "beta", "beta", "beta"), 
  fuzzy = c(0.18, 0.169, 0.149, 0.192, 0.154, 0.161, 0.169, 0.197, 0.162, 0.201, 
            0.215, 0.202, 0.17, 0.227, 0.196, 0.215, 0.161, 0.195, 0.178, 
            0.23, 0.184, 0.172, 0.153, 0.191, 0.156, 0.167, 0.165, 0.184, 
            0.162, 0.194, 0.197, 0.183, 0.151, 0.208, 0.16, 0.169, 0.169, 
            0.2, 0.154, 0.208), stringsAsFactors = FALSE)
head(sc_classification)
```

Our question is **which type of cells did we have on the original sample?**  

We can easily answer this by looking at the relations that have higher confidence of the relationship for each cell.

```{r cells_classification}
sc_classification %>% 
  group_by(elements) %>% 
  filter(fuzzy == max(fuzzy)) %>% 
  group_by(sets) %>% 
  count()
```

There is a cell that can be in two cell types, because we started with 10 cells and we have 11 elements here. 
However, how likely is that a cell is placed to just a single set? 

```{r cells_subset}
scTS <- tidySet(sc_classification) # Conversion of format
sample_cells <- scTS %>% 
  element_size() %>% 
  group_by(elements) %>% 
  filter(probability == max(probability))
sample_cells
```


There must be some cell misclassification: we have `r nElements(scTS)` cells in this example but there are `r sum(sample_cells$size)` cells that the maximum probability predicted for these types is a single cell type. 
Even the cell that had two equal probabilities for two cell types is more probable to be in no one of these cell types than in any of them. 

Ideally the predicted number of cells per type and the cells with higher confidence about the type should match. 

We can also look the other way around: How good is the prediction of a cell type for each cell? 

```{r celltypes}
scTS %>% 
  set_size() %>% 
  group_by(sets) %>% 
  filter(probability == max(probability))
```

We can see that for each cell type it is probable to have at least one cell and in the endothelial cell type two cells is the most probable outcome. 
However, these probabilities are lower than the probabilities of cells being assigned a cell type.
This would mean that this method is not a good method or that the cell types are not specific enough for the cell. 

In summary, the cells that we had most probable are not those 4 cell types except in two cells were it might be. 


# Session info {.unnumbered}

```{r sessionInfo, echo=FALSE}
sessionInfo()
```
---
title: "Advanced examples"
abstract: >
  This vignette assumes you are familiar with set operations from the basic vignette.
date: "`r format(Sys.time(), '%Y %b %d')`"
output:
  html_document:
    fig_caption: true
    code_folding: show
    self_contained: yes
    toc_float:
      collapsed: true
      toc_depth: 3
author:
- name: Lluís Revilla
  affiliation: 
    - August Pi i Sunyer Biomedical Research Institute (IDIBAPS); Liver Unit, Hospital Clinic
  email: lluis.revilla@gmail.com
vignette: >
  %\VignetteIndexEntry{Advanced examples}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteDepends{GO.db}
  %\VignetteDepends{reactome.db}
  %\VignetteDepends{GSEABase}
  %\VignetteDepends{org.Hs.eg.db}
  %\DeclareUnicodeCharacter{2229}{$\cap$}
  %\DeclareUnicodeCharacter{222A}{$\cup$}
editor_options: 
  chunk_output_type: console
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, 
                      warning = TRUE,
                      comment = "#>")
run_vignette <- requireNamespace("GSEABase", quietly = TRUE) && requireNamespace("GO.db", quietly = TRUE) &&  requireNamespace("reactome.db", quietly = TRUE) &&  requireNamespace("org.Hs.eg.db", quietly = TRUE)
```


# Initial setup

To show compatibility with tidy workflows we will use magrittr pipe operator and the
dplyr verbs.

```{r setupr, message=FALSE}
library("BaseSet")
library("dplyr")
```

# Human gene ontology

We will explore the genes with assigned gene ontology terms. 
These terms describe what is the process and role of the genes.
The links are annotated with different [evidence codes](http://geneontology.org/docs/guide-go-evidence-codes/) to indicate how such annotation is supported. 

```{r prepare_GO, message=FALSE, eval=run_vignette}
# We load some libraries
library("org.Hs.eg.db")
library("GO.db")
library("ggplot2")
# Prepare the data 
h2GO_TS <- tidySet(org.Hs.egGO)
h2GO <- as.data.frame(org.Hs.egGO)
```

We can now explore if there are differences in evidence usage for each ontology in gene ontology:

```{r evidence_ontology, eval=run_vignette}
library("forcats")
h2GO %>% 
  group_by(Evidence, Ontology) %>% 
  count(name = "Freq") %>% 
  ungroup() %>% 
  mutate(Evidence = fct_reorder2(Evidence, Ontology, -Freq),
         Ontology = case_when(Ontology == "CC" ~ "Cellular Component",
                              Ontology == "MF" ~ "Molecular Function",
                              Ontology == "BP" ~ "Biological Process",
                              TRUE ~ NA_character_)) %>% 
  ggplot() +
  geom_col(aes(Evidence, Freq)) +
  facet_grid(~Ontology) + 
  theme_minimal() +
  coord_flip() +
  labs(x = element_blank(), y = element_blank(),
       title = "Evidence codes for each ontology")
```

We can see that biological process are more likely to be defined by IMP evidence code that means inferred from mutant phenotype. 
While inferred from physical interaction (IPI) is almost exclusively used to assign molecular functions. 

This graph doesn't consider that some relationships are better annotated than other:

```{r nEvidence_plot, eval=run_vignette}
h2GO_TS %>% 
  relations() %>% 
  group_by(elements, sets) %>% 
  count(sort = TRUE, name = "Annotations") %>% 
  ungroup() %>% 
  count(Annotations, sort = TRUE) %>% 
  ggplot() +
  geom_col(aes(Annotations, n)) +
  theme_minimal() +
  labs(x = "Evidence codes", y = "Annotations", 
       title = "Evidence codes for each annotation",
       subtitle = "in human") +
  scale_x_continuous(breaks = 1:7)
```

We can see that mostly all the annotations are done with a single evidence code.
So far we have explored the code that it is related to a gene but how many genes don't have any annotation?

```{r numbers, eval=run_vignette}
# Add all the genes and GO terms
h2GO_TS <- add_elements(h2GO_TS, keys(org.Hs.eg.db)) %>% 
  add_sets(grep("^GO:", keys(GO.db), value = TRUE))

sizes_element <- element_size(h2GO_TS) %>% 
    arrange(desc(size))
sum(sizes_element$size == 0)
sum(sizes_element$size != 0)

sizes_set <- set_size(h2GO_TS) %>% 
    arrange(desc(size))
sum(sizes_set$size == 0)
sum(sizes_set$size != 0)
```

So we can see that both there are more genes without annotation and more gene ontology terms without a (direct) gene annotated.

```{r plots_GO, eval=run_vignette}
sizes_element %>% 
    filter(size != 0) %>% 
    ggplot() +
    geom_histogram(aes(size), binwidth = 1) +
    theme_minimal() +
    labs(x = "# sets per element", y = "Count")

sizes_set %>% 
    filter(size != 0) %>% 
    ggplot() +
    geom_histogram(aes(size), binwidth = 1) +
    theme_minimal() +
    labs(x = "# elements per set", y = "Count")
```

As you can see on the second plot we have very large values but that are on associated on many genes:

```{r distr_sizes, eval=run_vignette}
head(sizes_set, 10)
```


## Using fuzzy values

This could radically change if we used fuzzy values. 
We could assign a fuzzy value to each evidence code given the lowest fuzzy value for the [IEA (Inferred from Electronic Annotation)](http://wiki.geneontology.org/index.php/Inferred_from_Electronic_Annotation_(IEA)) evidence.
The highest values would be for evidence codes coming from experiments or alike. 

```{r fuzzy_setup, eval=run_vignette}
nr <- h2GO_TS %>% 
  relations() %>% 
  dplyr::select(sets, Evidence) %>% 
  distinct() %>% 
  mutate(fuzzy = case_when(
    Evidence == "EXP" ~ 0.9,
    Evidence == "IDA" ~ 0.8,
    Evidence == "IPI" ~ 0.8,
    Evidence == "IMP" ~ 0.75,
    Evidence == "IGI" ~ 0.7,
    Evidence == "IEP" ~ 0.65,
    Evidence == "HEP" ~ 0.6,
    Evidence == "HDA" ~ 0.6,
    Evidence == "HMP" ~ 0.5,
    Evidence == "IBA" ~ 0.45,
    Evidence == "ISS" ~ 0.4,
    Evidence == "ISO" ~ 0.32,
    Evidence == "ISA" ~ 0.32,
    Evidence == "ISM" ~ 0.3,
    Evidence == "RCA" ~ 0.2,
    Evidence == "TAS" ~ 0.15,
    Evidence == "NAS" ~ 0.1,
    Evidence == "IC" ~ 0.02,
    Evidence == "ND" ~ 0.02,
    Evidence == "IEA" ~ 0.01,
    TRUE ~ 0.01)) %>% 
  dplyr::select(sets = "sets", elements = "Evidence", fuzzy = fuzzy)
```

We have several evidence codes for the same ontology, this would result on different fuzzy values for each relation. 
Instead, we extract this and add them as new sets and elements and add an extra column to classify what are those elements:

```{r fuzzy_setup2, eval=run_vignette}
ts <- h2GO_TS %>% 
  relations() %>% 
  dplyr::select(-Evidence) %>% 
  rbind(nr) %>% 
  tidySet() %>% 
  mutate_element(Type = ifelse(grepl("^[0-9]+$", elements), "gene", "evidence"))
```

Now we can see which gene ontologies are more supported by the evidence:

```{r cardinality, eval=run_vignette}
ts %>% 
  dplyr::filter(Type != "Gene") %>% 
  cardinality() %>% 
  arrange(desc(cardinality)) %>% 
  head()
```

Surprisingly the most supported terms are protein binding, nucleus and cytosol. 
I would expect them to be the top three terms for cellular component, biological function and molecular function. 

Calculating set sizes would be interesting but it requires computing a big number of combinations that make it last long and require many memory available.

```{r size_go, eval=run_vignette}
ts %>% 
  filter(sets %in% c("GO:0008152", "GO:0003674", "GO:0005575"),
         Type != "gene") %>% 
  set_size()
```

Unexpectedly there is few evidence for the main terms:

```{r evidence_go, eval=run_vignette}
ts %>% 
  filter(sets %in% c("GO:0008152", "GO:0003674", "GO:0005575")) %>% 
  filter(Type != "gene") 
```

In fact those terms are arbitrarily decided or inferred from electronic analysis. 

# Human pathways

Now we will repeat the same analysis with pathways:

```{r prepare_reactome, eval=run_vignette}
# We load some libraries
library("reactome.db")

# Prepare the data (is easier, there isn't any ontoogy or evidence column)
h2p <- as.data.frame(reactomeEXTID2PATHID)
colnames(h2p) <- c("sets", "elements")
# Filter only for human pathways
h2p <- h2p[grepl("^R-HSA-", h2p$sets), ]

# There are duplicate relations with different evidence codes!!: 
summary(duplicated(h2p[, c("elements", "sets")]))
h2p <- unique(h2p)
# Create a tidySet and 
h2p_TS <- tidySet(h2p) %>% 
    # Add all the genes 
    add_elements(keys(org.Hs.eg.db))
```

Now that we have everything ready we can start measuring some things...

```{r numbers_pathways, eval=run_vignette}
sizes_element <- element_size(h2p_TS) %>% 
    arrange(desc(size))
sum(sizes_element$size == 0)
sum(sizes_element$size != 0)

sizes_set <- set_size(h2p_TS) %>% 
    arrange(desc(size))
```

We can see there are more genes without pathways than genes with pathways.

```{r pathways_plots, eval=run_vignette}
sizes_element %>% 
    filter(size != 0) %>% 
    ggplot() +
    geom_histogram(aes(size), binwidth = 1) +
    scale_y_log10() +
    theme_minimal() +
    labs(x = "# sets per element", y = "Count")

sizes_set %>% 
    ggplot() +
    geom_histogram(aes(size), binwidth = 1) +
    scale_y_log10() +
    theme_minimal() +
    labs(x = "# elements per set", y = "Count")
```

As you can see on the second plot we have very large values but that are on associated on many genes:

```{r distr_sizes_pathways, eval=run_vignette}
head(sizes_set, 10)
```

# Session info {.unnumbered}

```{r sessionInfo, echo=FALSE}
sessionInfo()
```
---
title: "About BaseSet"
abstract: >
  Describes the background of the package, important functions defined in the
  package and some of the applications and usages.
date: "`r format(Sys.time(), '%Y %b %d')`"
output:
  html_document:
    fig_caption: true
    code_folding: show
    self_contained: yes
    toc_float:
      collapsed: true
      toc_depth: 3
author:
- name: Lluís Revilla
  affiliation: 
    - August Pi i Sunyer Biomedical Research Institute (IDIBAPS); Liver Unit, Hospital Clinic
  email: lluis.revilla@gmail.com
vignette: >
  %\VignetteIndexEntry{About BaseSet}
  %\VignetteEncoding{UTF-8}  
  %\VignetteEngine{knitr::rmarkdown}
  %\DeclareUnicodeCharacter{2229}{$\cap$}
  %\DeclareUnicodeCharacter{222A}{$\cup$}
editor_options: 
  chunk_output_type: console
---
```{r setup, message=FALSE, warning=FALSE, include=FALSE}
knitr::opts_knit$set(root.dir = ".")
knitr::opts_chunk$set(collapse = TRUE, 
                      warning = TRUE,
                      comment = "#>")
```
# Getting started

This vignette explains how to use the methods available in this package.


# The TidySet class

This is a basic example which shows you how to create a `TidySet` object, to store associations between genes and sets:

```{r from_list, message=FALSE}
library("BaseSet")
gene_lists <- list(
    geneset1 = c("A", "B"),
    geneset2 = c("B", "C", "D")
)
tidy_set <- tidySet(gene_lists)
tidy_set
```

This is then stored internally in three slots `relations`, `elements`, and `sets` slots.

If you have more information for each element or set it can be added:

```{r metadata, message=FALSE}
gene_data <- data.frame(
    stat1     = c( 1,   2,   3,   4 ),
    info1     = c("a", "b", "c", "d")
)

tidy_set <- add_column(tidy_set, "elements", gene_data)
set_data <- data.frame(
    Group     = c(      100,        200 ),
    Colum     = c(     "abc",      "def")
)
tidy_set <- add_column(tidy_set, "sets", set_data)
tidy_set
```

This data is stored in one of the three slots, which can be directly accessed using their getter methods:

```{r getters}
relations(tidy_set)
elements(tidy_set)
sets(tidy_set)
```

You can add as much information as you want, with the only restriction for a "fuzzy" column for the `relations`. See the Fuzzy sets vignette.

# Creating a TidySet

As you can see it is possible to create a TidySet from a list and a data.frame, but it is also possible from a matrix:

```{r tidySet_matrix}
m <- matrix(c(0, 0, 1, 1, 1, 1, 0, 1, 0), ncol = 3, nrow =3,  
               dimnames = list(letters[1:3], LETTERS[1:3]))
m
tidy_set <- tidySet(m)
```

Or they can be created from a GeneSet and GeneSetCollection objects. Additionally it has several function to read files related to sets like the OBO files (`getOBO`) and GAF (`getGAF`)

# Converting to other formats

It is possible to extract the gene sets as a `list`, for use with functions such as `lapply`.

```{r as.list}
as.list(tidy_set)
```

Or if you need to apply some network methods and you need a matrix, you can create it with `incidence`:

```{r incidence}
incidence(tidy_set)
```

# Operations with sets

To work with sets several methods are provided. In general you can provide a new name for the resulting set of the operation, but if you don't one will be automatically provided using `naming`. All methods work with fuzzy and non-fuzzy sets

## Union

You can make a union of two sets present on the same object.

```{r union}
BaseSet::union(tidy_set, sets = c("C", "B"), name = "D")
```

## Intersection


```{r intersection}
intersection(tidy_set, sets = c("A", "B"), name = "D", keep = TRUE)
```

The keep argument used here is if you want to keep all the other previous sets:

```{r intersection2}
intersection(tidy_set, sets = c("A", "B"), name = "D", keep = FALSE)
```

## Complement

We can look for the complement of one or several sets:

```{r complement}
complement_set(tidy_set, sets = c("A", "B"))
```

Observe that we haven't provided a name for the resulting set but we can provide one if we prefer to

```{r complement2}
complement_set(tidy_set, sets = c("A", "B"), name = "F")
```

## Subtract

This is the equivalent of `setdiff`, but clearer:

```{r subtract}
out <- subtract(tidy_set, set_in = "A", not_in = "B", name = "A-B")
out
name_sets(out)
subtract(tidy_set, set_in = "B", not_in = "A", keep = FALSE)
```

See that in the first case there isn't any element present in B not in set A, but the new set is stored. 
In the second use case we focus just on the elements that are present on B but not in A.

# Additional information

The number of unique elements and sets can be obtained using the `nElements` and `nSets` methods.

```{r n}
nElements(tidy_set)
nSets(tidy_set)
nRelations(tidy_set)
```

The size of each gene set can be obtained using the `set_size` method.

```{r set_size}
set_size(tidy_set, "A")
```

Conversely, the number of sets associated with each gene is returned by the 
`element_size` function.

```{r element_size}
element_size(tidy_set)
```

The identifiers of elements and sets can be inspected and renamed using `name_elements` and 

```{r name}
name_elements(tidy_set)
name_elements(tidy_set) <- paste0("Gene", seq_len(nElements(tidy_set)))
name_elements(tidy_set)
name_sets(tidy_set)
name_sets(tidy_set) <- paste0("Geneset", seq_len(nSets(tidy_set)))
name_sets(tidy_set)
```


# Using `dplyr` verbs

You can also use `mutate`, `filter` and other `dplyr` verbs with TidySets (with the only exception being `group_by`), but you usually need to activate which three slots you want to affect with `activate`:

```{r tidyverse}
library("dplyr")
m_TS <- tidy_set %>% 
  activate("relations") %>% 
  mutate(Important = runif(nRelations(tidy_set)))
m_TS
```

You can use activate to select what are the verbs modifying:

```{r deactivate}
set_modified <- m_TS %>% 
  activate("elements") %>% 
  mutate(Pathway = if_else(elements %in% c("Gene1", "Gene2"), 
                           "pathway1", 
                           "pathway2"))
set_modified
set_modified %>% 
  deactivate() %>% # To apply a filter independently of where it is
  filter(Pathway == "pathway1")
```


If you think you need group_by usually this would mean that you need a new set. You can create a new one with `group`. 
If you want to use `group_by` to group some elements then you need to create a 
new set:

```{r group}
# A new group of those elements in pathway1 and with Important == 1
set_modified %>% 
  deactivate() %>% 
  group(name = "new", Pathway == "pathway1")
```
```{r group2}
set_modified %>% 
  group("pathway1", elements %in% c("Gene1", "Gene2"))
```

After grouping or mutating sometimes we might be interested in moving a column describing something to other places. We can do by this with:

```{r moving}
elements(set_modified)
out <- move_to(set_modified, "elements", "relations", "Pathway")
relations(out)
```

# Session info {.unnumbered}

```{r sessionInfo, echo=FALSE}
sessionInfo()
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arrange.R
\name{arrange.TidySet}
\alias{arrange.TidySet}
\alias{arrange_set}
\alias{arrange_element}
\alias{arrange_relation}
\title{Arrange the order of a TidySet}
\usage{
\method{arrange}{TidySet}(.data, ...)

arrange_set(.data, ...)

arrange_element(.data, ...)

arrange_relation(.data, ...)
}
\arguments{
\item{.data}{The TidySet object}

\item{...}{Comma separated list of variables names or expressions
integer column position to be used to reorder the TidySet.}
}
\value{
A TidySet object
}
\description{
Use arrange to extract the columns of a TidySet object. You can use activate
with filter or use the specific function. The S3 method filters using all
the information on the TidySet.
}
\examples{
relations <- data.frame(
    sets = c(rep("A", 5), "B", rep("A2", 5), "B2"),
    elements = rep(letters[seq_len(6)], 2),
    fuzzy = runif(12)
)
a <- tidySet(relations)
a <- mutate_element(a,
    type = c(rep("Gene", 4), rep("lncRNA", 2))
)

b <- arrange(a, desc(type))
elements(b)
b <- arrange_element(a, elements)
elements(b)
# Arrange sets
arrange_set(a, sets)
}
\seealso{
dplyr \code{\link[dplyr]{arrange}} and \code{\link{activate}}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/cardinality.R
\name{cardinality}
\alias{cardinality}
\alias{cardinality,TidySet-method}
\title{Cardinality or membership of sets}
\usage{
cardinality(object, sets = NULL, ...)

\S4method{cardinality}{TidySet}(object, sets, FUN = "sum", ...)
}
\arguments{
\item{object}{A TidySet object.}

\item{sets}{Character vector with the name of the sets.}

\item{...}{Other arguments passed to \code{FUN}.}

\item{FUN}{Function that returns a single numeric value given a vector of
fuzzy values.}
}
\description{
Calculates the membership of sets according to the logic defined in FUN.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Cardinality of sets
}}

\examples{
rel <- list(A = letters[1:3], B = letters[1:2])
TS <- tidySet(rel)
cardinality(TS, "A")
}
\seealso{
\code{\link[=size]{size()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy-set.R
\name{tidySet}
\alias{tidySet}
\alias{tidySet.data.frame}
\alias{tidySet.list}
\alias{tidySet.matrix}
\alias{tidySet.Go3AnnDbBimap}
\title{Create a TidySet object}
\usage{
tidySet(relations)

\method{tidySet}{data.frame}(relations)

\method{tidySet}{list}(relations)

\method{tidySet}{matrix}(relations)

\method{tidySet}{Go3AnnDbBimap}(relations)
}
\arguments{
\item{relations}{An object to be coerced to a TidySet.}
}
\value{
A TidySet object.
}
\description{
These functions help to create a \code{TidySet} object from
\code{data.frame}, \code{list}, \code{matrix}, and \code{GO3AnnDbBimap}.
They can create both fuzzy and standard sets.
}
\details{
Elements or sets without any relation are not shown when printed.
}
\section{Methods (by class)}{
\itemize{
\item \code{data.frame}: Given the relations in a data.frame

\item \code{list}: Convert to a TidySet from a list

\item \code{matrix}: Convert an incidence matrix into a TidySet

\item \code{Go3AnnDbBimap}: Convert Go3AnnDbBimap into a TidySet object.
}}

\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b"),
    elements = letters[seq_len(6)]
)
tidySet(relations)
relations2 <- data.frame(
    sets = c(rep("A", 5), "B"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
tidySet(relations2)
# A
x <- list("A" = letters[1:5], "B" = LETTERS[3:7])
tidySet(x)
# A fuzzy set taken encoded as a list
A <- runif(5)
names(A) <- letters[1:5]
B <- runif(5)
names(B) <- letters[3:7]
relations <- list(A, B)
tidySet(relations)
# Will error
# x <- list("A" = letters[1:5], "B" = LETTERS[3:7], "c" = runif(5))
# a <- tidySet(x) # Only characters or factors are allowed as elements.
M <- matrix(c(1, 0.5, 1, 0), ncol = 2,
            dimnames = list(c("A", "B"), c("a", "b")))
tidySet(M)
}
\seealso{
\code{\link{TidySet-class}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/move_to.R
\name{move_to}
\alias{move_to}
\alias{move_to,TidySet,characterORfactor,characterORfactor,character-method}
\title{Move columns between slots}
\usage{
move_to(object, from, to, columns)

\S4method{move_to}{TidySet,characterORfactor,characterORfactor,character}(object, from, to, columns)
}
\arguments{
\item{object}{A TidySet object.}

\item{from}{The name of the slot where the content is.}

\item{to}{The name of the slot to move the content.}

\item{columns}{The name of the columns that should be moved.}
}
\value{
A TidySet object where the content is moved from one slot to other.
}
\description{
Moves information from one slot to other slots.
For instance from the sets to the relations.
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,from = characterORfactor,to = characterORfactor,columns = character}: Move columns
}}

\examples{
x <- list("A" = c("a" = 0.1, "b" = 0.5), "B" = c("a" = 0.2, "b" = 1))
TS <- tidySet(x)
TS <- mutate_element(TS, b = runif(2))
TS2 <- move_to(TS, from = "elements", to = "relations", "b")
# Note that apparently we haven't changed anything:
TS2
}
\seealso{
Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complement.R
\name{complement}
\alias{complement}
\title{Complement TidySet}
\usage{
complement(.data, ...)
}
\arguments{
\item{.data}{The TidySet object}

\item{...}{Other arguments passed to either \code{\link{complement_set}} or
\code{\link{complement_element}}.}
}
\value{
A TidySet object
}
\description{
Use complement to find elements or sets the TidySet object. You can use
activate with complement or use the specific function. You must specify if
you want the complements of sets or elements.
}
\examples{
rel <- data.frame(
    sets = c("A", "A", "B", "B", "C", "C"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(rel)
TS \%>\%
    activate("elements") \%>\%
    complement("a")
TS \%>\%
    activate("elements") \%>\%
    complement("a", "C_a", keep = FALSE)
TS \%>\%
    activate("set") \%>\%
    complement("A")
TS \%>\%
    activate("set") \%>\%
    complement("A", keep = FALSE)
TS \%>\%
    activate("set") \%>\%
    complement("A", FUN = function(x){abs(x - 0.2)}, keep = FALSE)
}
\seealso{
\code{\link{activate}}

Other complements: 
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{subtract}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{complements}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/add_column.R
\name{add_column}
\alias{add_column}
\alias{add_column,TidySet,character-method}
\title{Add column}
\usage{
add_column(object, slot, columns)

\S4method{add_column}{TidySet,character}(object, slot, columns)
}
\arguments{
\item{object}{A TidySet object.}

\item{slot}{A TidySet slot.}

\item{columns}{The columns to add.}
}
\value{
A \code{TidySet} object.
}
\description{
Add column to a slot of the TidySet object.
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,slot = character}: Add a column to any slot
}}

\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
add_column(TS, "relations", data.frame(well = c(
    "GOOD", "BAD", "WORSE",
    "UGLY", "FOE", "HEY"
)))
}
\seealso{
\code{\link{rename_set}}

Other column: 
\code{\link{remove_column}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{column}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/elements.R
\name{elements}
\alias{elements}
\alias{elements<-}
\alias{elements,TidySet-method}
\alias{elements<-,TidySet-method}
\alias{replace_elements}
\alias{nElements,TidySet-method}
\title{Elements of the TidySet}
\usage{
elements(object)

elements(object) <- value

\S4method{elements}{TidySet}(object)

\S4method{elements}{TidySet}(object) <- value

replace_elements(object, value)

\S4method{nElements}{TidySet}(object)
}
\arguments{
\item{object}{A TidySet object.}

\item{value}{Modification of the elements.}
}
\value{
A \code{data.frame} with information about the elements
}
\description{
Given TidySet retrieve the elements or substitute them.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Retrieve the elements

\item \code{TidySet}: Modify the elements

\item \code{TidySet}: Return the number of elements
}}

\examples{
TS <- tidySet(list(A = letters[1:5], B = letters[2:10]))
elements(TS)
elements(TS) <- data.frame(elements = letters[10:1])
TS2 <- replace_elements(TS, data.frame(elements = letters[1:11]))
nElements(TS)
nElements(TS2)
}
\seealso{
\code{\link{nElements}}

Other slots: 
\code{\link{relations}()},
\code{\link{sets}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{slots}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/relations.R
\name{is.fuzzy}
\alias{is.fuzzy}
\alias{is.fuzzy,TidySet-method}
\title{Check if a TidySet is fuzzy.}
\usage{
is.fuzzy(object)

\S4method{is.fuzzy}{TidySet}(object)
}
\arguments{
\item{object}{Object to be coerced or tested.}
}
\value{
A logical value.
}
\description{
Check if there are fuzzy sets. A fuzzy set is a set where the relationship
between elements is given by a probability (or uncertainty).
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Check if it is fuzzy
}}

\examples{
TS <- tidySet(list(A = letters[1:5], B = letters[2:10]))
is.fuzzy(TS)
}
\seealso{
Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nested.R
\name{is_nested}
\alias{is_nested}
\alias{is_nested.TidySet}
\title{Are some sets as elements of other sets?}
\usage{
is_nested(object)

\method{is_nested}{TidySet}(object)
}
\arguments{
\item{object}{A TidySet object.}
}
\value{
A logical value: TRUE if there are some sets included as elements of
others.
}
\description{
Check if some elements are also sets of others. This is also known as
hierarchical sets.
}
\examples{
relations <- list(A = letters[1:3], B = c(letters[4:5]))
TS <- tidySet(relations)
is_nested(TS)
TS2 <- add_relation(TS, data.frame(elements = "A", sets = "B"))
# Note that A is both a set and an element of B
TS2
is_nested(TS2)
}
\seealso{
adjacency

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/rename.R
\name{rename_elements}
\alias{rename_elements}
\alias{rename_elements,TidySet-method}
\title{Rename elements}
\usage{
rename_elements(object, old, new)

\S4method{rename_elements}{TidySet}(object, old, new)
}
\arguments{
\item{object}{A TidySet object.}

\item{old}{A character vector of to be renamed.}

\item{new}{A character vector of with new names.}
}
\value{
A \code{TidySet} object.
}
\description{
Change the default names of sets and elements.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Rename elements
}}

\examples{
x <- list("A" = letters[1:5], "B" = letters[3:7])
TS <- tidySet(x)
name_elements(TS)
TS2 <- rename_elements(TS, "a", "first")
name_elements(TS2)
}
\seealso{
\code{\link{name_elements}}

Other renames: 
\code{\link{rename_set}()}

Other names: 
\code{\link{name_elements<-}()},
\code{\link{name_elements}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{rename_set}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{names}
\concept{renames}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/activate.R, R/deactivate.R
\name{activate}
\alias{activate}
\alias{active}
\alias{deactivate}
\title{Determine the context of subsequent manipulations.}
\usage{
activate(.data, what)

active(.data)

deactivate(.data)
}
\arguments{
\item{.data}{A \code{TidySet} object.}

\item{what}{Either "elements", "sets" or "relations"}
}
\value{
A \code{TidySet} object.
}
\description{
Functions to help to perform some action to just some type of data: elements,
sets or relations.
\code{activate}: To table the focus of future manipulations: elements, sets
or relations.
\code{active}: To check the focus on the \code{TidySet}.
\code{deactivate}: To remove the focus on a specific \code{TidySet}-
}
\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b", rep("a2", 5), "b2"),
    elements = rep(letters[seq_len(6)], 2),
    fuzzy = runif(12)
)
a <- tidySet(relations)
elements(a) <- cbind(elements(a),
    type = c(rep("Gene", 4), rep("lncRNA", 2))
)
# Filter in the whole TidySet
filter(a, elements == "a")
filter(a, elements == "a", type == "Gene")
# Equivalent to filter_elements
filter_element(a, type == "Gene")
a <- activate(a, "elements")
active(a)
filter(a, type == "Gene")
a <- deactivate(a)
active(a)
filter(a, type == "Gene")
}
\seealso{
Other methods: 
\code{\link{TidySet-class}},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filter.TidySet}
\alias{filter.TidySet}
\alias{filter_set}
\alias{filter_element}
\alias{filter_relation}
\title{Filter TidySet}
\usage{
\method{filter}{TidySet}(.data, ...)

filter_set(.data, ...)

filter_element(.data, ...)

filter_relation(.data, ...)
}
\arguments{
\item{.data}{The TidySet object.}

\item{...}{The logical predicates in terms of the variables of the sets.}
}
\value{
A TidySet object.
}
\description{
Use filter to subset the TidySet object. You can use activate with filter or
use the specific function. The S3 method filters using all the information
on the TidySet.
}
\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b", rep("a2", 5), "b2"),
    elements = rep(letters[seq_len(6)], 2),
    fuzzy = runif(12),
    type = c(rep("Gene", 4), rep("lncRNA", 2))
)
TS <- tidySet(relations)
TS <- move_to(TS, from = "relations", to = "elements", column = "type")
filter(TS, elements == "a")
# Equivalent to filter_relation
filter(TS, elements == "a", sets == "a")
filter_relation(TS, elements == "a", sets == "a")
# Filter element
filter_element(TS, type == "Gene")
# Filter sets and by property of elements simultaneously
filter(TS, sets == "b", type == "lncRNA")
# Filter sets
filter_set(TS, sets == "b")
}
\seealso{
dplyr \code{\link[dplyr]{filter}} and \code{\link{activate}}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMT.R
\name{getGMT}
\alias{getGMT}
\title{Import GMT (Gene Matrix Transposed) files}
\usage{
getGMT(con, sep = "\\t", ...)
}
\arguments{
\item{con}{File name of the GMT file.}

\item{sep}{GMT file field separator, by default tabs.}

\item{...}{Other arguments passed to \code{readLines}.}
}
\value{
A TidySet object.
}
\description{
The GMT (Gene Matrix Transposed) file format is a tab delimited file format
that describes groups of genes. In this format, each row represents a group.
Each group is described by a name, a description, and the genes in it.
}
\examples{
gmtFile <- system.file(
    package = "BaseSet", "extdata",
    "hallmark.gene.symbol.gmt"
)
gs <- getGMT(gmtFile)
nRelations(gs)
nElements(gs)
nSets(gs)
}
\references{
The file format is defined by the Broad Institute \href{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}{here}
}
\seealso{
Other IO functions: 
\code{\link{getGAF}()},
\code{\link{getOBO}()}
}
\concept{IO functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_frame.R
\name{df2TS}
\alias{df2TS}
\title{The opposite of as.data.frame}
\usage{
df2TS(.data = NULL, df)
}
\arguments{
\item{.data}{The original TidySet}

\item{df}{The flattened data.frame}
}
\value{
A TidySet object
}
\description{
Convert a data.frame to a TidySet by first using the relations.
It requires the original TidySet in order to convert it back to resemble
the position of the columns.
}
\seealso{
\code{\link[=tidySet.data.frame]{tidySet.data.frame()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list.R
\name{as.list.TidySet}
\alias{as.list.TidySet}
\title{Convert to list}
\usage{
\method{as.list}{TidySet}(x, ...)
}
\arguments{
\item{x}{A TidySet object to be coerced to a list.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}
}
\value{
A list.
}
\description{
Converts a TidySet to a list.
}
\examples{
r <- data.frame(sets = c("A", "A", "A", "B", "C"),
             elements = c(letters[1:3], letters[2:3]),
             fuzzy = runif(5),
             info = rep_len(c("important", "very important"), 5))
TS <- tidySet(r)
TS
as.list(TS)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/complement.R
\name{complement_set}
\alias{complement_set}
\alias{complement_set,TidySet,characterORfactor-method}
\title{Complement of a set}
\usage{
complement_set(object, sets, ...)

\S4method{complement_set}{TidySet,characterORfactor}(
  object,
  sets,
  name = NULL,
  FUN = NULL,
  keep = TRUE,
  keep_relations = keep,
  keep_elements = keep,
  keep_sets = keep
)
}
\arguments{
\item{object}{A TidySet object.}

\item{sets}{The name of the set to look for the complement.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}

\item{name}{Name of the new set. By default it adds a "C".}

\item{FUN}{A function to be applied when performing the union.
The standard union is the "max" function, but you can provide any other
function that given a numeric vector returns a single number.}

\item{keep}{Logical value to keep all the other sets.}

\item{keep_relations}{A logical value if you wan to keep old relations.}

\item{keep_elements}{A logical value if you wan to keep old elements.}

\item{keep_sets}{A logical value if you wan to keep old sets.}
}
\value{
A \code{TidySet} object.
}
\description{
Return the complement for a set
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,sets = characterORfactor}: Complement of the sets.
}}

\examples{
relations <- data.frame(
    sets = c("A", "A", "B", "B", "C", "C"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
complement_set(TS, "A")
}
\seealso{
\code{\link{filter}}

Other complements: 
\code{\link{complement_element}()},
\code{\link{complement}()},
\code{\link{subtract}()}

Other methods that create new sets: 
\code{\link{complement_element}()},
\code{\link{intersection}()},
\code{\link{subtract}()},
\code{\link{union}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{complements}
\concept{methods}
\concept{methods that create new sets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/intersection.R
\name{intersection}
\alias{intersection}
\alias{intersect}
\alias{intersection,TidySet,character-method}
\title{Intersection of two or more sets}
\usage{
intersection(object, sets, ...)

\S4method{intersection}{TidySet,character}(
  object,
  sets,
  name = NULL,
  FUN = "min",
  keep = FALSE,
  keep_relations = keep,
  keep_elements = keep,
  keep_sets = keep,
  ...
)
}
\arguments{
\item{object}{A TidySet object.}

\item{sets}{The character of sets to be intersect.}

\item{...}{Other named arguments passed to \code{FUN}.}

\item{name}{The name of the new set. By defaults joins the sets with an
\ifelse{latex}{\out{$\cup$}}{\ifelse{html}{\out{&cup;}}{}}.}

\item{FUN}{A function to be applied when performing the union.
The standard intersection is the "min" function, but you can provide any
other function that given a numeric vector returns a single number.}

\item{keep}{A logical value if you want to keep originals sets.}

\item{keep_relations}{A logical value if you wan to keep old relations.}

\item{keep_elements}{A logical value if you wan to keep old elements.}

\item{keep_sets}{A logical value if you wan to keep old sets.}
}
\value{
A \code{TidySet} object.
}
\description{
Given a TidySet creates a new set with the elements on the both of them
following the logic defined on FUN.
}
\details{
#' The default uses the \code{min} function following the \href{https://en.wikipedia.org/wiki/Fuzzy_set_operations}{standard fuzzy definition}, but it can be
changed.
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,sets = character}: Applies the standard intersection
}}

\examples{
rel <- data.frame(
    sets = c(rep("A", 5), "B"),
    elements = c("a", "b", "c", "d", "f", "f")
)
TS <- tidySet(rel)
intersection(TS, c("A", "B")) # Default Name
intersection(TS, c("A", "B"), "C") # Set the name
# Fuzzy set
rel <- data.frame(
    sets = c(rep("A", 5), "B"),
    elements = c("a", "b", "c", "d", "f", "f"),
    fuzzy = runif(6)
)
TS2 <- tidySet(rel)
intersection(TS2, c("A", "B"), "C")
intersection(TS2, c("A", "B"), "C", FUN = function(x){max(sqrt(x))})
}
\seealso{
Other methods that create new sets: 
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{subtract}()},
\code{\link{union}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{methods that create new sets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/obo.R
\name{getGAF}
\alias{getGAF}
\title{Read a GAF file}
\usage{
getGAF(x)
}
\arguments{
\item{x}{A file in GAF format}
}
\value{
A TidySet object
}
\description{
Read a GO Annotation File (GAF) formatted file
}
\examples{
gafFile <- system.file(
    package = "BaseSet", "extdata",
    "go_human_rna_valid_subset.gaf"
)
gs <- getGAF(gafFile)
head(gs)
}
\references{
The format is defined \href{http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/}{here}.
}
\seealso{
Other IO functions: 
\code{\link{getGMT}()},
\code{\link{getOBO}()}
}
\concept{IO functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{show,TidySet-method}
\alias{show,TidySet-method}
\title{Method to show the TidySet object}
\usage{
\S4method{show}{TidySet}(object)
}
\arguments{
\item{object}{A TidySet}
}
\value{
A table with the information of the relationships.
}
\description{
Prints the resulting table of a TidySet object. Does not shown elements or
sets without any relationship (empty sets). To see them use \code{\link[=sets]{sets()}} or
\code{\link[=elements]{elements()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.R
\name{select.TidySet}
\alias{select.TidySet}
\alias{select_set}
\alias{select_element}
\alias{select_relation}
\title{select from a TidySet}
\usage{
\method{select}{TidySet}(.data, ...)

select_set(.data, ...)

select_element(.data, ...)

select_relation(.data, ...)
}
\arguments{
\item{.data}{The TidySet object}

\item{...}{The name of the columns you want to keep, remove or rename.}
}
\value{
A TidySet object
}
\description{
Use select to extract the columns of a TidySet object. You can use activate
with filter or use the specific function. The S3 method filters using all
the information on the TidySet.
}
\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b", rep("a2", 5), "b2"),
    elements = rep(letters[seq_len(6)], 2),
    fuzzy = runif(12)
)
a <- tidySet(relations)
a <- mutate_element(a,
    type = c(rep("Gene", 4), rep("lncRNA", 2))
)
a <- mutate_set(a, Group = c("UFM", "UAB", "UPF", "MIT"))
b <- select(a, -type)
elements(b)
b <- select_element(a, elements)
elements(b)
# Select sets
select_set(a, sets)
}
\seealso{
dplyr \code{\link[dplyr]{select}} and \code{\link{activate}}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/names.R
\name{name_sets}
\alias{name_sets}
\alias{name_sets,TidySet-method}
\alias{name_sets<-,TidySet,characterORfactor-method}
\title{Name sets}
\usage{
name_sets(object)

\S4method{name_sets}{TidySet}(object)

\S4method{name_sets}{TidySet,characterORfactor}(object) <- value
}
\arguments{
\item{object}{A TidySet object.}

\item{value}{A character with the new names for the sets.}
}
\value{
A \code{TidySet} object.
}
\description{
Retrieve the name of the sets.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Name sets

\item \code{object = TidySet,value = characterORfactor}: Rename sets
}}

\examples{
relations <- data.frame(
    sets = c(rep("A", 5), "B"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
name_sets(TS)
}
\seealso{
Other names: 
\code{\link{name_elements<-}()},
\code{\link{name_elements}()},
\code{\link{name_sets<-}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{names}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basesets-package.R
\docType{package}
\name{BaseSet-package}
\alias{BaseSet}
\alias{BaseSet-package}
\title{BaseSet: Working with Sets the Tidy Way}
\description{
Implements a class and methods to work with sets,
    doing intersection, union, complementary sets, power sets, cartesian
    product and other set operations in a "tidy" way. These set operations
    are available for both classical sets and fuzzy sets. Import sets from
    several formats or from other several data structures.
}
\details{
It provides a class \code{\link{TidySet}} with methods to do operations with sets.
}
\examples{
set <- list("A" = letters[1:5], "B" = letters[4:7])
TS <- tidySet(set)
cardinality(TS)
intersection(TS, c("A", "B"))
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/BaseSet}
  \item \url{https://docs.ropensci.org/BaseSet/}
  \item Report bugs at \url{https://github.com/ropensci/BaseSet/issues}
}

}
\author{
\strong{Maintainer}: Lluís Revilla Sancho \email{lluis.revilla@gmail.com} (\href{https://orcid.org/0000-0001-9747-2570}{ORCID}) [copyright holder]

Other contributors:
\itemize{
  \item Zebulun Arendsee [reviewer]
  \item Jennifer Chang [reviewer]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/length.R
\name{union_probability}
\alias{union_probability}
\alias{length_probability}
\title{Calculates the probability of a single length}
\usage{
union_probability(p)

length_probability(p, size)
}
\arguments{
\item{p}{Numeric vector of probabilities.}

\item{size}{Integer value of the size of the selected values.}
}
\value{
A numeric value of the probability of the given size.
}
\description{
Creates all the possibilities and then add them up.
\code{union_probability} Assumes independence between the probabilities to
calculate the final size.
}
\examples{
length_probability(c(0.5, 0.75), 2)
length_probability(c(0.5, 0.75, 0.66), 1)
length_probability(c(0.5, 0.1, 0.3, 0.5, 0.25, 0.23), 2)
union_probability(c(0.5, 0.1, 0.3))
}
\seealso{
\code{\link[=multiply_probabilities]{multiply_probabilities()}} and \code{\link[=length_set]{length_set()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/remove.R
\name{remove_relation}
\alias{remove_relation}
\alias{remove_relation,TidySet,characterORfactor,characterORfactor-method}
\title{Remove a relation}
\usage{
remove_relation(object, elements, sets, ...)

\S4method{remove_relation}{TidySet,characterORfactor,characterORfactor}(object, elements, sets)
}
\arguments{
\item{object}{A TidySet object}

\item{elements}{The elements of the sets.}

\item{sets}{The name of the new set.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}
}
\value{
A \code{TidySet} object.
}
\description{
Given a TidySet removes relations between elements and sets
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,elements = characterORfactor,sets = characterORfactor}: Removes a relation between elements and sets.
}}

\examples{
relations <- data.frame(
    sets = c(rep("A", 5), "B"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
remove_relation(TS, "A", "a")
}
\seealso{
Other remove functions: 
\code{\link{remove_element}()},
\code{\link{remove_set}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{remove functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_frame.R
\name{as.data.frame.TidySet}
\alias{as.data.frame.TidySet}
\title{Transforms a TidySet to a data.frame}
\usage{
\method{as.data.frame}{TidySet}(x, ...)
}
\arguments{
\item{x}{The \code{TidySet} object.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}
}
\value{
A \code{data.frame} table.
}
\description{
Flattens the three slots to a single big table
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R
\name{add_sets}
\alias{add_sets}
\title{Add sets to a TidySet}
\usage{
add_sets(object, sets, ...)
}
\arguments{
\item{object}{A \code{\link{TidySet}} object}

\item{sets}{A character vector of sets to be added.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}
}
\value{
A \code{\link{TidySet}} object with the new sets.
}
\description{
Functions to add sets. If the sets are new they are added,
otherwise they are omitted.
}
\note{
\code{add_sets} doesn't set up any other information about the sets.
Remember to add/modify them if needed with \code{\link{mutate}} or \code{\link{mutate_set}}
}
\examples{
x <- list("a" = letters[1:5], "b" = LETTERS[3:7])
a <- tidySet(x)
b <- add_sets(a, "fg")
sets(b)
}
\seealso{
Other add_*: 
\code{\link{add_elements}()},
\code{\link{add_relations}()}
}
\concept{add_*}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/set.R
\name{sets}
\alias{sets}
\alias{sets<-}
\alias{sets,TidySet-method}
\alias{sets<-,TidySet-method}
\alias{replace_sets}
\alias{nSets,TidySet-method}
\title{Sets of the TidySet}
\usage{
sets(object)

sets(object) <- value

\S4method{sets}{TidySet}(object)

\S4method{sets}{TidySet}(object) <- value

replace_sets(object, value)

\S4method{nSets}{TidySet}(object)
}
\arguments{
\item{object}{A \code{SetCollection} object.}

\item{value}{Modification of the sets.}
}
\value{
A \code{data.frame} with information from the sets.
}
\description{
Given TidySet retrieve the sets or substitute them.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Retrieve the sets information

\item \code{TidySet}: Modify the sets information

\item \code{TidySet}: Return the number of sets
}}

\examples{
TS <- tidySet(list(A = letters[1:5], B = letters[2:10]))
sets(TS)
sets(TS) <- data.frame(sets = c("B", "A"))
TS2 <- replace_sets(TS, data.frame(sets = c("A", "B", "C")))
sets(TS2)
nSets(TS)
nSets(TS2)
}
\seealso{
\code{\link{nSets}}

Other slots: 
\code{\link{elements}()},
\code{\link{relations}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{slots}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/relations.R
\name{relations}
\alias{relations}
\alias{relations<-}
\alias{relations,TidySet-method}
\alias{replace_relations}
\alias{relations<-,TidySet-method}
\alias{nRelations,TidySet-method}
\title{Relations of the TidySet}
\usage{
relations(object)

relations(object) <- value

\S4method{relations}{TidySet}(object)

replace_relations(object, value)

\S4method{relations}{TidySet}(object) <- value

\S4method{nRelations}{TidySet}(object)
}
\arguments{
\item{object}{Object to be coerced or tested.}

\item{value}{Modification of the relations.}
}
\value{
A \code{data.frame} with information about the relations between
elements and sets.
}
\description{
Given TidySet retrieve the relations or substitute them.
\code{\link{TidySet}} object
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Retrieve the relations

\item \code{TidySet}: Modify the relations

\item \code{TidySet}: Return the number of unique relations
}}

\examples{
TS <- tidySet(list(A = letters[1:2], B = letters[5:7]))
relations(TS)
}
\seealso{
\code{\link{nRelations}}

Other slots: 
\code{\link{elements}()},
\code{\link{sets}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{slots}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/union.R
\name{union}
\alias{union}
\alias{union.TidySet}
\title{Join sets}
\usage{
union(object, ...)

\method{union}{TidySet}(
  object,
  sets,
  name = NULL,
  FUN = "max",
  keep = FALSE,
  keep_relations = keep,
  keep_elements = keep,
  keep_sets = keep,
  ...
)
}
\arguments{
\item{object}{A TidySet object.}

\item{...}{Other named arguments passed to \code{FUN}.}

\item{sets}{The name of the sets to be used.}

\item{name}{The name of the new set. By defaults joins the sets with an
\ifelse{latex}{\out{$\cap$}}{\ifelse{html}{\out{&cap;}}{}}.}

\item{FUN}{A function to be applied when performing the union.
The standard union is the "max" function, but you can provide any other
function that given a numeric vector returns a single number.}

\item{keep}{A logical value if you want to keep.}

\item{keep_relations}{A logical value if you wan to keep old relations.}

\item{keep_elements}{A logical value if you wan to keep old elements.}

\item{keep_sets}{A logical value if you wan to keep old sets.}
}
\value{
A \code{TidySet} object.
}
\description{
Given a TidySet merges several sets into the new one using the logic
defined on FUN.
}
\details{
The default uses the \code{max} function following the \href{https://en.wikipedia.org/wiki/Fuzzy_set_operations}{standard fuzzy definition}, but it can be
changed. See examples below.
}
\examples{
# Classical set
rel <- data.frame(
    sets = c(rep("A", 5), "B", "B"),
    elements = c(letters[seq_len(6)], "a")
)
TS <- tidySet(rel)
union(TS, c("B", "A"))
# Fuzzy set
rel <- data.frame(
    sets = c(rep("A", 5), "B", "B"),
    elements = c(letters[seq_len(6)], "a"),
    fuzzy = runif(7)
)
TS2 <- tidySet(rel)
# Standard default logic
union(TS2, c("B", "A"), "C")
# Probability logic
union(TS2, c("B", "A"), "C", FUN = union_probability)
}
\seealso{
\code{\link[=union_probability]{union_probability()}}

Other methods that create new sets: 
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{intersection}()},
\code{\link{subtract}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()}
}
\concept{methods}
\concept{methods that create new sets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/remove_column.R
\name{remove_column}
\alias{remove_column}
\alias{remove_column,TidySet,character,character-method}
\title{Remove column}
\usage{
remove_column(object, slot, column_names)

\S4method{remove_column}{TidySet,character,character}(object, slot, column_names)
}
\arguments{
\item{object}{A TidySet object.}

\item{slot}{A TidySet slot.}

\item{column_names}{The name of the columns.}
}
\value{
A \code{TidySet} object.
}
\description{
Removes column from a slot of the TidySet object.
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,slot = character,column_names = character}: Remove columns to any slot
}}

\examples{
x <- data.frame(sets = c(rep("A", 5), rep("B", 5)),
                elements = c(letters[1:5], letters[3:7]),
                extra = sample(c("YES", "NO"), 10, replace = TRUE))
TS <- tidySet(x)
TS
remove_column(TS, "relations", "extra")
}
\seealso{
\code{\link{rename_set}}

Other column: 
\code{\link{add_column}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{column}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{nSets}
\alias{nSets}
\title{Number of sets}
\usage{
nSets(object)
}
\arguments{
\item{object}{Object to be coerced or tested.}
}
\value{
The number of sets present.
}
\description{
Check the number of sets of the TidySet
}
\examples{
TS <- tidySet(list(A = letters[1:2], B = letters[5:7]))
nSets(TS)
}
\seealso{
Other count functions: 
\code{\link{nElements}()},
\code{\link{nRelations}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{count functions}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutate.R
\name{mutate.TidySet}
\alias{mutate.TidySet}
\alias{mutate_set}
\alias{mutate_element}
\alias{mutate_relation}
\title{Mutate}
\usage{
\method{mutate}{TidySet}(.data, ...)

mutate_set(.data, ...)

mutate_element(.data, ...)

mutate_relation(.data, ...)
}
\arguments{
\item{.data}{The TidySet object.}

\item{...}{The logical predicates in terms of the variables of the sets.}
}
\value{
A TidySet object
}
\description{
Use mutate to alter the TidySet object. You can use activate with mutate or
use the specific function. The S3 method filters using all the information
on the TidySet.
}
\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b", rep("a2", 5), "b2"),
    elements = rep(letters[seq_len(6)], 2),
    fuzzy = runif(12)
)
a <- tidySet(relations)
a <- mutate_element(a, Type = c(rep("Gene", 4), rep("lncRNA", 2)))
a
b <- mutate_relation(a, Type = sample(c("PPI", "PF", "MP"), 12,
    replace = TRUE
))
}
\seealso{
\code{\link[dplyr]{mutate}} and \code{\link{activate}}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{nRelations}
\alias{nRelations}
\title{Number of relations}
\usage{
nRelations(object)
}
\arguments{
\item{object}{Object to be coerced or tested.}
}
\value{
A numeric value with the number of the relations.
}
\description{
Check the number of relations of the TidySet.
}
\examples{
TS <- tidySet(list(A = letters[1:2], B = letters[5:7]))
nRelations(TS)
}
\seealso{
Other count functions: 
\code{\link{nElements}()},
\code{\link{nSets}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{count functions}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/names.R
\name{name_elements}
\alias{name_elements}
\alias{name_elements,TidySet-method}
\alias{name_elements<-,TidySet,characterORfactor-method}
\title{Name elements}
\usage{
name_elements(object)

\S4method{name_elements}{TidySet}(object)

\S4method{name_elements}{TidySet,characterORfactor}(object) <- value
}
\arguments{
\item{object}{A TidySet object.}

\item{value}{A character with the new names for the elements.}
}
\value{
A \code{TidySet} object.
}
\description{
Retrieve the name of the elements.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Name elements

\item \code{object = TidySet,value = characterORfactor}: Rename elements
}}

\examples{
relations <- data.frame(
    sets = c(rep("A", 5), "B"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
name_elements(TS)
}
\seealso{
Other names: 
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()}
}
\concept{names}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_set.R
\name{power_set}
\alias{power_set}
\title{Create the power set}
\usage{
power_set(object, set, name, ...)
}
\arguments{
\item{object}{A TidySet object.}

\item{set}{The name of the set to be used for the power set}

\item{name}{The root name of the new set.}

\item{...}{Other arguments passed down if possible. Currently ignored.}
}
\value{
A TidySet object with the new set
}
\description{
Create the power set
}
\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b"),
    elements = letters[seq_len(6)]
)
TS <- tidySet(relations)
power_set(TS, "a", name = "power_set")
}
\seealso{
Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/size.R
\name{size}
\alias{size}
\title{Size}
\usage{
size(object, ...)
}
\arguments{
\item{object}{A TidySet object}

\item{...}{Character vector with the name of elements or sets you want to
calculate the size of.}
}
\value{
The size of the elements or sets. If there is no active slot or it
is the relations slot returns the TidySet object with a warning.
}
\description{
Calculate the size of the elements or sets, using the fuzzy values as
probabilities. First it must have active either sets or elements.
}
\examples{
rel <- data.frame(
    sets = c(rep("A", 5), "B", "C"),
    elements = c(letters[seq_len(6)], letters[6])
)
TS <- tidySet(rel)
TS <- activate(TS, "elements")
size(TS)
TS <- activate(TS, "sets")
size(TS)
# With fuzzy sets
relations <- data.frame(
    sets = c(rep("A", 5), "B", "C"),
    elements = c(letters[seq_len(6)], letters[6]),
    fuzzy = runif(7)
)
TS <- tidySet(relations)
TS <- activate(TS, "elements")
size(TS)
TS <- activate(TS, "sets")
size(TS)
}
\seealso{
A related concept \code{\link[=cardinality]{cardinality()}}. It is calculated using
\code{\link[=length_set]{length_set()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R
\name{add_relations}
\alias{add_relations}
\title{Add relations to a TidySet}
\usage{
add_relations(object, elements, sets, fuzzy, ...)
}
\arguments{
\item{object}{A \code{\link{TidySet}} object}

\item{elements}{A character vector of the elements.}

\item{sets}{A character vector of sets to be added.}

\item{fuzzy}{The strength of the membership.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}
}
\value{
A \code{\link{TidySet}} object with the new relations.
}
\description{
Adds new relations to existing or new sets and elements.
If the sets or elements do not exist they are added.
}
\note{
\code{add_relations} doesn't set up any other information about the
relationship.
Remember to add/modify them if needed with \code{\link{mutate}} or \code{\link{mutate_relation}}
}
\examples{
x <- list("a" = letters[1:5], "b" = LETTERS[3:7])
a <- tidySet(x)
add_relations(a, elements = c("a", "b", "g"), sets = "d")
add_relations(a, elements = c("a", "b"), sets = c("d", "g"))
add_relations(a, elements = c("a", "b"), sets = c("d", "g"), fuzzy = 0.5)
add_relations(a,
    elements = c("a", "b"), sets = c("d", "g"),
    fuzzy = c(0.5, 0.7)
)
}
\seealso{
\code{\link{add_relation}} to add relations with new sets or/and
new elements.

Other add_*: 
\code{\link{add_elements}()},
\code{\link{add_sets}()}
}
\concept{add_*}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/droplevels.R
\name{droplevels.TidySet}
\alias{droplevels.TidySet}
\title{Drop unused elements and sets}
\usage{
\method{droplevels}{TidySet}(x, elements = TRUE, sets = TRUE, relations = TRUE, ...)
}
\arguments{
\item{x}{A TidySet object.}

\item{elements}{Logical value: Should elements be dropped?}

\item{sets}{Logical value: Should sets be dropped?}

\item{relations}{Logical value: Should sets be dropped?}

\item{...}{Other arguments, currently ignored.}
}
\value{
A TidySet object.
}
\description{
Drop elements and sets without any relation.
}
\examples{
rel <- list(A = letters[1:3], B = character())
TS <- tidySet(rel)
TS
sets(TS)
TS2 <- droplevels(TS)
TS2
sets(TS2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/naming.R
\name{naming}
\alias{naming}
\title{Name an operation}
\usage{
naming(
  start = NULL,
  sets1,
  middle = NULL,
  sets2 = NULL,
  collapse_symbol = "union"
)
}
\arguments{
\item{start, middle}{Character used as a start symbol or to divide
\code{sets1} and \code{sets2}.}

\item{sets1, sets2}{Character of sets}

\item{collapse_symbol}{Name of the symbol that joins the sets on
\code{sets1} and \code{sets2}.}
}
\value{
A character vector combining the sets
}
\description{
Helps setting up the name of an operation.
}
\examples{
naming(sets1 = c("a", "b"))
naming(sets1 = "a", middle = "union", sets2 = "b")
naming(sets1 = "a", middle = "intersection", sets2 = c("b", "c"))
naming(sets1 = "a", middle = "intersection", sets2 = c("b", "c"))
naming(
    start = "complement", sets1 = "a", middle = "intersection",
    sets2 = c("b", "c"), collapse_symbol = "intersection"
)
}
\seealso{
\code{\link{set_symbols}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/length.R
\name{multiply_probabilities}
\alias{multiply_probabilities}
\alias{independent_probabilities}
\title{Probability of a vector of probabilities}
\usage{
multiply_probabilities(p, i)

independent_probabilities(p, i)
}
\arguments{
\item{p}{Numeric vector of probabilities.}

\item{i}{Numeric integer index of the complementary probability.}
}
\value{
A numeric value of the probability.
}
\description{
Calculates the probability that all probabilities happened simultaneously.
\code{independent_probabilities} just multiply the probabilities of the index passed.
}
\examples{
multiply_probabilities(c(0.5, 0.1, 0.3, 0.5, 0.25, 0.23), c(1, 3))
independent_probabilities(c(0.5, 0.1, 0.3, 0.5, 0.25, 0.23), c(1, 3))
}
\seealso{
\code{\link[=length_probability]{length_probability()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\value{
The result of the left-hand side (lhs) is passed to the right-hand
side (rhs).
}
\description{
See \code{magrittr::\link[magrittr]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/subtract.R
\name{subtract}
\alias{subtract}
\alias{subtract,TidySet,characterORfactor,characterORfactor-method}
\title{Subtract}
\usage{
subtract(object, set_in, not_in, ...)

\S4method{subtract}{TidySet,characterORfactor,characterORfactor}(
  object,
  set_in,
  not_in,
  name = NULL,
  keep = TRUE,
  keep_relations = keep,
  keep_elements = keep,
  keep_sets = keep
)
}
\arguments{
\item{object}{A TidySet object.}

\item{set_in}{Name of the sets where the elements should be present.}

\item{not_in}{Name of the sets where the elements should not be present.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}

\item{name}{Name of the new set. By default it adds a "C".}

\item{keep}{Logical value to keep all the other sets.}

\item{keep_relations}{A logical value if you wan to keep old relations.}

\item{keep_elements}{A logical value if you wan to keep old elements.}

\item{keep_sets}{A logical value if you wan to keep old sets.}
}
\value{
A \code{TidySet} object.
}
\description{
Elements in a set not present in the other set. Equivalent to
\code{\link{setdiff}}.
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,set_in = characterORfactor,not_in = characterORfactor}: Elements present in sets but not in other sets
}}

\examples{
relations <- data.frame(
    sets = c("A", "A", "B", "B", "C", "C"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
subtract(TS, "A", "B")
subtract(TS, "A", "B", keep = FALSE)
}
\seealso{
\code{\link{setdiff}}

Other complements: 
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()}

Other methods that create new sets: 
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{intersection}()},
\code{\link{union}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{union}()}
}
\concept{complements}
\concept{methods}
\concept{methods that create new sets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneSetCollection.R
\name{tidy}
\alias{tidy}
\alias{tidy.GeneSetCollection}
\alias{tidy.GeneSet}
\title{Convert GSEABase classes to a TidySet}
\usage{
tidy(object)

\method{tidy}{GeneSetCollection}(object)

\method{tidy}{GeneSet}(object)
}
\arguments{
\item{object}{A GeneSetCollection or a GeneSet derived object}
}
\value{
A TidySet object
}
\description{
Convert GSEABase classes to a TidySet
}
\section{Methods (by class)}{
\itemize{
\item \code{GeneSetCollection}: Converts to a tidySet given a GeneSetCollection

\item \code{GeneSet}: Converts to a tidySet given a GeneSet
}}

\examples{
# Needs GSEABase pacakge from Bioconductor
if (requireNamespace("GSEABase", quietly = TRUE)) {
    library("GSEABase")
    gs <- GeneSet()
    gs
    tidy(gs)
    fl <- system.file("extdata", "Broad.xml", package="GSEABase")
    gs2 <- getBroadSets(fl) # actually, a list of two gene sets
    TS <- tidy(gs2)
    dim(TS)
    sets(TS)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arrange.R, R/filter.R, R/group_by.R,
%   R/mutate.R, R/pull.R, R/select.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{arrange}
\alias{filter}
\alias{group_by}
\alias{mutate}
\alias{pull}
\alias{select}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{dplyr}{\code{\link[dplyr]{arrange}}, \code{\link[dplyr]{filter}}, \code{\link[dplyr]{group_by}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{pull}}, \code{\link[dplyr]{select}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/rename.R
\name{rename_set}
\alias{rename_set}
\alias{rename_set,TidySet-method}
\title{Rename sets}
\usage{
rename_set(object, old, new)

\S4method{rename_set}{TidySet}(object, old, new)
}
\arguments{
\item{object}{A TidySet object.}

\item{old}{A character vector of to be renamed.}

\item{new}{A character vector of with new names.}
}
\value{
A \code{TidySet} object.
}
\description{
Change the default names of sets and elements.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Rename sets
}}

\examples{
x <- list("A" = letters[1:5], "B" = letters[3:7])
TS <- tidySet(x)
name_sets(TS)
TS2 <- rename_set(TS, "A", "C")
name_sets(TS2)
}
\seealso{
\code{\link{name_sets}}

Other renames: 
\code{\link{rename_elements}()}

Other names: 
\code{\link{name_elements<-}()},
\code{\link{name_elements}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{rename_elements}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{names}
\concept{renames}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{name_elements<-}
\alias{name_elements<-}
\title{Rename elements}
\usage{
name_elements(object) <- value
}
\arguments{
\item{object}{A TidySet object.}

\item{value}{A character with the new names for the elements.}
}
\value{
A \code{TidySet} object.
}
\description{
Rename elements.
}
\examples{
relations <- data.frame(
    sets = c(rep("A", 5), "B"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
TS
name_elements(TS) <- letters[1:6]
}
\seealso{
\code{\link{rename_elements}}

Other names: 
\code{\link{name_elements}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{names}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/length.R
\name{set_size}
\alias{set_size}
\alias{set_size,TidySet-method}
\title{Calculates the size of a set}
\usage{
set_size(object, sets = NULL)

\S4method{set_size}{TidySet}(object, sets = NULL)
}
\arguments{
\item{object}{A TidySet object.}

\item{sets}{The sets from which the length is calculated.}
}
\value{
A list with the size of the set or the probability of having that
size.
}
\description{
Assuming that the fuzzy values are probabilities,
calculates the probability of being of different sizes for a given set.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Calculates the size of a set using \code{\link[=length_set]{length_set()}}
}}

\examples{
relations <- data.frame(
    sets = c(rep("A", 5), "B", "C"),
    elements = c(letters[seq_len(6)], letters[6]),
    fuzzy = runif(7)
)
a <- tidySet(relations)
set_size(a)
}
\seealso{
cardinality

Other sizes: 
\code{\link{element_size}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{sizes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{name_sets<-}
\alias{name_sets<-}
\title{Rename sets}
\usage{
name_sets(object) <- value
}
\arguments{
\item{object}{A TidySet object.}

\item{value}{A character with the new names for the sets.}
}
\value{
A \code{TidySet} object.
}
\description{
Rename sets.
}
\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
TS
name_sets(TS) <- LETTERS[1:2]
}
\seealso{
\code{\link{rename_set}}

Other names: 
\code{\link{name_elements<-}()},
\code{\link{name_elements}()},
\code{\link{name_sets}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{names}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjacency.R
\name{adjacency}
\alias{adjacency}
\alias{adjacency_element}
\alias{adjacency_set}
\alias{adjacency.TidySet}
\title{Adjacency}
\usage{
\method{adjacency}{TidySet}(object)

adjacency_element(object)

adjacency_set(object)

\method{adjacency}{TidySet}(object)
}
\arguments{
\item{object}{A TidySet object}
}
\value{
A square matrix, 1 if two nodes are connected, 0 otherwise.
}
\description{
Are two elements connected ?
}
\examples{
x <- list("SET1" = letters[1:5], "SET2" = LETTERS[3:7])
a <- tidySet(x)
adjacency_element(a)
adjacency_set(a)
}
\seealso{
\code{\link{incidence}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/length.R
\name{length_set}
\alias{length_set}
\title{Calculates the probability}
\usage{
length_set(probability)
}
\arguments{
\item{probability}{A numeric vector of probabilities.}
}
\value{
A vector with the probability of each set.
}
\description{
Given several probabilities it looks for how probable is to have a vector of
each length
}
\examples{
length_set(c(0.5, 0.1, 0.3, 0.5, 0.25, 0.23))
}
\seealso{
\code{\link[=length_probability]{length_probability()}} to calculate the probability of a specific
length.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_by.R
\name{group_by.TidySet}
\alias{group_by.TidySet}
\title{group_by TidySet}
\usage{
\method{group_by}{TidySet}(.data, ...)
}
\arguments{
\item{.data}{The TidySet object}

\item{...}{The logical predicates in terms of the variables of the sets}
}
\value{
A grouped data.frame (See The dplyr help page)
}
\description{
Use group_by to group the TidySet object. You can use activate with
group_by or with the whole data.
}
\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b", rep("a2", 5), "b2"),
    elements = rep(letters[seq_len(6)], 2),
    fuzzy = runif(12)
)
a <- tidySet(relations)
elements(a) <- cbind(elements(a),
    type = c(rep("Gene", 4), rep("lncRNA", 2))
)
group_by(a, elements)
}
\seealso{
dplyr \code{\link[dplyr]{group_by}} and \code{\link{activate}}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/complement.R
\name{complement_element}
\alias{complement_element}
\alias{complement_element,TidySet,characterORfactor-method}
\title{Complement of elements}
\usage{
complement_element(object, elements, ...)

\S4method{complement_element}{TidySet,characterORfactor}(
  object,
  elements,
  name = NULL,
  FUN = NULL,
  keep = TRUE,
  keep_relations = keep,
  keep_elements = keep,
  keep_sets = keep
)
}
\arguments{
\item{object}{A TidySet object.}

\item{elements}{The set to look for the complement.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}

\item{name}{Name of the new set. By default it adds a "C".}

\item{FUN}{A function to be applied when performing the union.
The standard union is the "max" function, but you can provide any other
function that given a numeric vector returns a single number.}

\item{keep}{Logical value to keep all the other sets.}

\item{keep_relations}{A logical value if you wan to keep old relations.}

\item{keep_elements}{A logical value if you wan to keep old elements.}

\item{keep_sets}{A logical value if you wan to keep old sets.}
}
\value{
A \code{TidySet} object.
}
\description{
Return the objects without the elements listed
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,elements = characterORfactor}: Complement of the elements.
}}

\examples{
relations <- data.frame(
    sets = c("A", "A", "B", "B", "C", "C"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
complement_element(TS, "a", "C_a")
complement_element(TS, "a", "C_a", keep = FALSE)
}
\seealso{
Other complements: 
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{subtract}()}

Other methods that create new sets: 
\code{\link{complement_set}()},
\code{\link{intersection}()},
\code{\link{subtract}()},
\code{\link{union}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{complements}
\concept{methods}
\concept{methods that create new sets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/independent.R
\name{independent}
\alias{independent}
\title{Independence of the sets}
\usage{
independent(object, sets)
}
\arguments{
\item{object}{A \code{\link{TidySet}} object.}

\item{sets}{A character vector with the names of the sets to analyze.}
}
\value{
A logical value indicating if the sets are independent (TRUE) or not.
}
\description{
Checks if the elements of the sets are present in more than one set.
}
\examples{
x <- list("A" = letters[1:5], "B" = letters[3:7], "C" = letters[6:10])
TS <- tidySet(x)
independent(TS)
independent(TS, c("A", "B"))
independent(TS, c("A", "C"))
independent(TS, c("B", "C"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{TidySet-class}
\alias{TidySet-class}
\alias{TidySet}
\title{A tidy class to represent a set}
\description{
A set is a group of unique elements it can be either a fuzzy set, where the
relationship is between 0 or 1 or nominal.
}
\details{
When printed if an element or a set do not have any relationship is not
shown.
They can be created from lists, matrices or data.frames. Check \code{\link[=tidySet]{tidySet()}}
constructor for more information.
}
\section{Slots}{

\describe{
\item{\code{relations}}{A data.frame with elements and the sets were they belong.}

\item{\code{elements}}{A data.frame of unique elements and related information.}

\item{\code{sets}}{A data.frame of unique sets and related information.}
}}

\examples{
x <- list("A" = letters[1:5], "B" = LETTERS[3:7])
a <- tidySet(x)
a
x <- list("A" = letters[1:5], "B" = character())
b <- tidySet(x)
b
name_sets(b)
}
\seealso{
\link{tidySet}

Other methods: 
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/add_relation.R
\name{add_relation}
\alias{add_relation}
\alias{add_relation,TidySet,data.frame-method}
\title{Add relations}
\usage{
add_relation(object, relations, ...)

\S4method{add_relation}{TidySet,data.frame}(object, relations)
}
\arguments{
\item{object}{A TidySet object}

\item{relations}{A data.frame object}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}
}
\value{
A \code{TidySet} object.
}
\description{
Given a TidySet adds new relations between elements and sets.
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,relations = data.frame}: Adds relations
}}

\examples{
relations <- data.frame(
    sets = c(rep("A", 5), "B"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
relations <- data.frame(
    sets = c(rep("A2", 5), "B2"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6),
    new = runif(6)
)
add_relation(TS, relations)
}
\seealso{
Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{add functions}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pull.R
\name{pull.TidySet}
\alias{pull.TidySet}
\alias{pull_set}
\alias{pull_element}
\alias{pull_relation}
\title{Pull from a TidySet}
\usage{
\method{pull}{TidySet}(.data, var = -1, name = NULL, ...)

pull_set(.data, var = -1, name = NULL, ...)

pull_element(.data, var = -1, name = NULL, ...)

pull_relation(.data, var = -1, name = NULL, ...)
}
\arguments{
\item{.data}{The TidySet object}

\item{var}{The literal variable name, a positive integer or a negative
integer column position.}

\item{name}{Column used to name the output.}

\item{...}{Currently not used.}
}
\value{
A TidySet object
}
\description{
Use pull to extract the columns of a TidySet object. You can use activate
with filter or use the specific function. The S3 method filters using all
the information on the TidySet.
}
\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b", rep("a2", 5), "b2"),
    elements = rep(letters[seq_len(6)], 2),
    fuzzy = runif(12)
)
a <- tidySet(relations)
a <- mutate_element(a, type = c(rep("Gene", 4), rep("lncRNA", 2)))
pull(a, type)
# Equivalent to pull_relation
b <- activate(a, "relations")
pull_relation(b, elements)
pull_element(b, elements)
# Filter element
pull_element(a, type)
# Filter sets
pull_set(a, sets)
}
\seealso{
dplyr \code{\link[dplyr]{pull}} and \code{\link{activate}}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cartesian.R
\name{cartesian}
\alias{cartesian}
\alias{cartesian.TidySet}
\title{Create the cartesian product of two sets}
\usage{
cartesian(object, set1, set2, name = NULL, ...)

\method{cartesian}{TidySet}(
  object,
  set1,
  set2,
  name = NULL,
  keep = TRUE,
  keep_relations = keep,
  keep_elements = keep,
  keep_sets = keep,
  ...
)
}
\arguments{
\item{object}{A TidySet object.}

\item{set1, set2}{The name of the sets to be used for the cartesian product}

\item{name}{The name of the new set.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}

\item{keep}{A logical value if you want to keep.}

\item{keep_relations}{A logical value if you wan to keep old relations.}

\item{keep_elements}{A logical value if you wan to keep old elements.}

\item{keep_sets}{A logical value if you wan to keep old sets.}
}
\value{
A TidySet object with the new set
}
\description{
Given two sets creates new sets with one element of each set
}
\examples{
relations <- data.frame(
    sets = c(rep("a", 5), "b"),
    elements = letters[seq_len(6)]
)
TS <- tidySet(relations)
cartesian(TS, "a", "b")
}
\seealso{
Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/obo.R
\name{getOBO}
\alias{getOBO}
\title{Read an OBO file}
\usage{
getOBO(x)
}
\arguments{
\item{x}{Path to a file in OBO format.}
}
\value{
A TidySet object.
}
\description{
Read an Open Biological and Biomedical Ontologies (OBO) formatted file
}
\examples{
oboFile <- system.file(
    package = "BaseSet", "extdata",
    "go-basic_subset.obo"
)
gs <- getOBO(oboFile)
head(gs)
}
\references{
The format is described \href{https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html}{here}
}
\seealso{
Other IO functions: 
\code{\link{getGAF}()},
\code{\link{getGMT}()}
}
\concept{IO functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/length.R
\name{element_size}
\alias{element_size}
\alias{element_size,TidySet-method}
\title{Calculates the size of the elements}
\usage{
element_size(object, elements = NULL)

\S4method{element_size}{TidySet}(object, elements = NULL)
}
\arguments{
\item{object}{A TidySet object.}

\item{elements}{The element from which the length is calculated.}
}
\value{
A list with the size of the elements or the probability of having
that size.
}
\description{
Assuming that the fuzzy values are probabilities, calculates the probability
of being of different sizes for a given set.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Calculates the number of sets an element appears
with \code{\link[=length_set]{length_set()}}
}}

\examples{
relations <- data.frame(
    sets = c(rep("A", 5), "B", "C"),
    elements = c(letters[seq_len(6)], letters[6]),
    fuzzy = runif(7)
)
a <- tidySet(relations)
element_size(a)
}
\seealso{
cardinality

Other sizes: 
\code{\link{set_size}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{sizes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operations.R
\name{fapply}
\alias{fapply}
\title{Apply to fuzzy}
\usage{
fapply(relations, FUN, ...)
}
\arguments{
\item{relations}{A data.frame or similar with fuzzy, sets and elements
columns.}

\item{FUN}{A function to perform on the fuzzy numbers.}

\item{...}{Other named arguments passed to \code{FUN}.}
}
\value{
A modified TidySet object
}
\description{
Simplify and returns unique results of the object.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/incidence.R
\name{incidence}
\alias{incidence}
\alias{incidence,TidySet-method}
\title{Incidence}
\usage{
incidence(object)

\S4method{incidence}{TidySet}(object)
}
\arguments{
\item{object}{Object to be coerced or tested.}
}
\value{
A matrix with elements in rows and sets in columns where the values
indicate the relationship between the element and the set.
}
\description{
Check which elements are in which sets.
}
\section{Methods (by class)}{
\itemize{
\item \code{TidySet}: Incidence of the TidySet
}}

\examples{
x <- list("a" = letters[1:5], "b" = LETTERS[3:7])
a <- tidySet(x)
incidence(a)
}
\seealso{
\code{\link{adjacency}}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R
\name{add_elements}
\alias{add_elements}
\title{Add elements to a TidySet}
\usage{
add_elements(object, elements, ...)
}
\arguments{
\item{object}{A \code{\link{TidySet}} object}

\item{elements}{A character vector of the elements.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}
}
\value{
A \code{\link{TidySet}} object with the new elements.
}
\description{
Functions to add elements. If the elements are new they are added,
otherwise they are omitted.
}
\note{
\code{add_element} doesn't set up any other information about the elements.
Remember to add/modify them if needed with \code{\link{mutate}} or \code{\link{mutate_element}}
}
\examples{
x <- list("a" = letters[1:5], "b" = LETTERS[3:7])
a <- tidySet(x)
b <- add_elements(a, "fg")
elements(b)
}
\seealso{
Other add_*: 
\code{\link{add_relations}()},
\code{\link{add_sets}()}
}
\concept{add_*}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{data}
\name{set_symbols}
\alias{set_symbols}
\title{A subset of symbols related to sets}
\format{
An object of class \code{character} of length 16.
}
\usage{
set_symbols
}
\description{
Name and symbol of operations related to sets, including intersection
and union among others:
}
\examples{
set_symbols
}
\references{
\url{https://www.fileformat.info/info/unicode/category/Sm/list.htm}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/remove.R
\name{remove_element}
\alias{remove_element}
\alias{remove_element,TidySet,characterORfactor-method}
\title{Remove elements}
\usage{
remove_element(object, elements, ...)

\S4method{remove_element}{TidySet,characterORfactor}(object, elements)
}
\arguments{
\item{object}{A TidySet object.}

\item{elements}{The elements to be removed.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}
}
\value{
A \code{TidySet} object.
}
\description{
Given a TidySet remove elements and the related relations and if
required also the sets.
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,elements = characterORfactor}: Removes everything related to an element
}}

\examples{
relations <- data.frame(
    sets = c(rep("A", 5), "B"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
remove_element(TS, "c")
}
\seealso{
Other remove functions: 
\code{\link{remove_relation}()},
\code{\link{remove_set}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{remove functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group.R
\name{group}
\alias{group}
\alias{group.TidySet}
\title{Create a new set from existing elements}
\usage{
group(object, name, ...)

\method{group}{TidySet}(object, name, ...)
}
\arguments{
\item{object}{A TidySet object.}

\item{name}{The name of the new set.}

\item{...}{A logical condition to subset some elements.}
}
\value{
A TidySet object with the new set.
}
\description{
It allows to create a new set given some condition. If no element meet the
condition an empty set is created.
}
\examples{
x <- list("A" = c("a" = 0.1, "b" = 0.5), "B" = c("a" = 0.2, "b" = 1))
TS <- tidySet(x)
TS1 <- group(TS, "C", fuzzy < 0.5)
TS1
sets(TS1)
TS2 <- group(TS, "D", fuzzy < 0)
sets(TS2)
r <- data.frame(
    sets = c(rep("A", 5), "B", rep("A2", 5), "B2"),
    elements = rep(letters[seq_len(6)], 2),
    fuzzy = runif(12),
    type = c(rep("Gene", 2), rep("Protein", 2), rep("lncRNA", 2))
)
TS3 <- tidySet(r)
group(TS3, "D", sets \%in\% c("A", "A2"))
}
\seealso{
Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/remove.R
\name{remove_set}
\alias{remove_set}
\alias{remove_set,TidySet,characterORfactor-method}
\title{Remove sets}
\usage{
remove_set(object, sets, ...)

\S4method{remove_set}{TidySet,characterORfactor}(object, sets)
}
\arguments{
\item{object}{A TidySet object.}

\item{sets}{The sets to be removed.}

\item{...}{Placeholder for other arguments that could be passed to the
method. Currently not used.}
}
\value{
A \code{TidySet} object.
}
\description{
Given a TidySet remove sets and the related relations and if
required also the elements
}
\section{Methods (by class)}{
\itemize{
\item \code{object = TidySet,sets = characterORfactor}: Removes everything related to a set
}}

\examples{
relations <- data.frame(
    sets = c("A", "A", "B", "B", "C", "C"),
    elements = letters[seq_len(6)],
    fuzzy = runif(6)
)
TS <- tidySet(relations)
remove_set(TS, "B")
}
\seealso{
Other remove functions: 
\code{\link{remove_element}()},
\code{\link{remove_relation}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nElements}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{methods}
\concept{remove functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{nElements}
\alias{nElements}
\title{Number of elements}
\usage{
nElements(object)
}
\arguments{
\item{object}{Object to be coerced or tested.}
}
\value{
A numeric value with the number of elements.
}
\description{
Check the number of elements of the TidySet.
}
\examples{
TS <- tidySet(list(A = letters[1:2], B = letters[5:7]))
nElements(TS)
}
\seealso{
Other count functions: 
\code{\link{nRelations}()},
\code{\link{nSets}()}

Other methods: 
\code{\link{TidySet-class}},
\code{\link{activate}()},
\code{\link{add_column}()},
\code{\link{add_relation}()},
\code{\link{arrange.TidySet}()},
\code{\link{cartesian}()},
\code{\link{complement_element}()},
\code{\link{complement_set}()},
\code{\link{complement}()},
\code{\link{element_size}()},
\code{\link{elements}()},
\code{\link{filter.TidySet}()},
\code{\link{group_by.TidySet}()},
\code{\link{group}()},
\code{\link{incidence}()},
\code{\link{intersection}()},
\code{\link{is.fuzzy}()},
\code{\link{is_nested}()},
\code{\link{move_to}()},
\code{\link{mutate.TidySet}()},
\code{\link{nRelations}()},
\code{\link{nSets}()},
\code{\link{name_elements<-}()},
\code{\link{name_sets<-}()},
\code{\link{name_sets}()},
\code{\link{power_set}()},
\code{\link{pull.TidySet}()},
\code{\link{relations}()},
\code{\link{remove_column}()},
\code{\link{remove_element}()},
\code{\link{remove_relation}()},
\code{\link{remove_set}()},
\code{\link{rename_elements}()},
\code{\link{rename_set}()},
\code{\link{select.TidySet}()},
\code{\link{set_size}()},
\code{\link{sets}()},
\code{\link{subtract}()},
\code{\link{union}()}
}
\concept{count functions}
\concept{methods}
