---
title: 'phylogram: an R package for phylogenetic analysis with nested lists'
tags:
  - R
  - phylogeny
  - dendrogram
  - evolutionary tree
authors:
  - name: Shaun P. Wilkinson
    orcid: 0000-0002-7332-7931
    affiliation: 1
  - name: Simon K. Davy
    orcid: 0000-0003-3584-5356
    affiliation: 1
affiliations:
 - name: School of Biological Sciences, Victoria University of Wellington, P.O. Box 600, Wellington, New Zealand.
   index: 1
date: 21 June 2018
bibliography: paper.bib
---

# Summary

The R environment continues to gain popularity as a platform for
bioinformatic analysis, due to the reproducible 
code-based workflow and the many powerful analytical tools 
available in a suite of open-source packages such as ``ape`` 
[@Paradis2004], ``phangorn`` [@Schliep2011] and ``Phytools`` [@Revell2012]. 
These packages typically employ a tree structure known as the 
"phylo" object, whose primary element is an integer matrix 
with one row for each edge in the graph, 
and two columns giving the indices of the connecting nodes. 
This is a versatile and memory-efficient structure suitable 
for most applications encountered by evolutionary biologists, 
and hence a comprehensive array of tools has been developed for 
editing, analyzing and visualizing trees in this format. 

An alternative tree structure is the "dendrogram" object, 
whose nodes consist of deeply nested lists.
While less memory-efficient than matrix-based trees, 
a useful feature of this representation is its modularity, 
whereby the sub tree of a tree is itself a tree - 
a dendrogram within a dendrogram. 
This means that dendrograms are subsettable in the same 
way that standard lists are, 
which facilitates intuitive command-line tree manipulation. 
An especially powerful feature of this object type is that tree-editing 
operations can be carried out recursively 
using fast inbuilt functions in the "apply" family such as `dendrapply` 
and `lapply`. There is also a large and growing number of resources for 
manipulating and plotting dendrograms in contributed packages such as 
``dendextend`` [@Galili2015], and
hence bi-directional conversion between "dendrogram"
and "phylo" class objects would expand the range of tools available for 
both object types. 

Here, we introduce ``phylogram``, an R package for developing 
phylogenies as deeply-nested lists, converting trees between 
list- and matrix-type objects, importing and exporting trees 
to and from parenthetic text, and editing/manipulating dendrogram objects.
``phylogram`` is available from GitHub (<https://github.com/ropensci/phylogram>)
and CRAN (<https://CRAN.R-project.org/package=phylogram>), 
and version 2.1 of the package is archived to Zenodo (<http://dx.doi.org/10.5281/zenodo.1293634>).
A full reference manual with worked examples can be found at <https://cran.r-project.org/web/packages/phylogram/vignettes/phylogram-vignette.html>. 
Bug reports and other feedback are welcomed and can be directed to the GitHub issues
page at <http://github.com/ropensci/phylogram/issues>, 
or the ``phylogram`` google group at <https://groups.google.com/group/phylogram>.


# Acknowledgements

This software was developed with funding from a Rutherford Foundation Postdoctoral 
Research Fellowship from the Royal Society of New Zealand. The authors declare no 
competing interests.


# References
# phylogram

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00790/status.svg)](https://doi.org/10.21105/joss.00790)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/phylogram)](https://cran.r-project.org/package=phylogram)
[![](http://cranlogs.r-pkg.org/badges/grand-total/phylogram)](https://cran.r-project.org/package=phylogram)
[![](https://badges.ropensci.org/212_status.svg)](https://github.com/ropensci/onboarding/issues/212)
[![Build Status](https://travis-ci.org/ropensci/phylogram.svg?branch=master)](https://travis-ci.org/ropensci/phylogram)
[![codecov](https://codecov.io/github/ropensci/phylogram/branch/master/graphs/badge.svg)](https://codecov.io/github/ropensci/phylogram)
[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1293634.svg)](https://doi.org/10.5281/zenodo.1293634)
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--7332--7931-brightgreen.svg)](http://orcid.org/0000-0002-7332-7931)

--------------------------------------------------------------------------------

The **phylogram** R package is a tool for for developing 
phylogenetic trees as deeply-nested lists known as "dendrogram" objects. 
It provides functions for conversion between "dendrogram" and 
"phylo" class objects, as well as several tools for command-line tree 
manipulation and import/export via Newick parenthetic text.
This improves accessibility to the comprehensive range of object-specific 
analytical and tree-visualization functions found across a wide array of 
bioinformatic R packages.


### Installation

To download **phylogram** from CRAN and load the package, run

```R
install.packages("phylogram")
library("phylogram")
```

To download the latest development version from GitHub, run:

```R
devtools::install_github("ropensci/phylogram", build_vignettes = TRUE) 
library("phylogram")
```


### Example: reading and writing trees

Consider the simple example of a tree with three members named 
"A", "B" and "C", where "B" and "C" are more closely related
to eachother than they are to "A". 
An unweighted Newick string for this tree would be "(A,(B,C));".
This text can be imported as a 
dendrogram object using the `read.dendrogram` function 
as follows:

```R
library("phylogram")
newick <- "(A,(B,C));"
x <- read.dendrogram(text = newick)
plot(x)
```

The following command writes the object back to the console in 
Newick format without edge weights:

```R
write.dendrogram(x, edges = FALSE)
```

The syntax is similar when reading and writing text files, 
except that the `text` argument is replaced by `file` and a 
valid file path is passed to the function.

To convert the dendrogram to a "phylo" object, run

```R
y <- as.phylo(x)
```

These and more examples are available in the package vignette.
To view the vignette, run `vignette("phylogram-vignette")` or access it 
directly from [CRAN](https://CRAN.R-project.org/package=phylogram).


### Help

An overview of the package with links to the function documentation can be found by running

```R
?phylogram
```

If you experience a problem using this package please
either raise it as an issue on [GitHub](http://github.com/ropensci/phylogram/issues) 
or post it on the **phylogram** [google group](https://groups.google.com/group/phylogram).


### Acknowledgements

This software was developed at 
[Victoria University of Wellington](http://www.victoria.ac.nz/) 
with funding from a Rutherford Foundation Postdoctoral Research Fellowship 
award from the Royal Society of New Zealand.


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
#phylogram 2.1.0

* Addresses performance issues when converting between "dendrogram" and "phylo" objects
* More examples and detail added to vignette
* `as.cladogram` replaces `ultrametricize`
* `read.dendrogram` now wraps `ape::read.tree` and converts to dendrogram via `as.dendrogram.phylo`
* Thanks very much to the ROpenSci editors and reviewers for improvements and suggestions

#phylogram 2.0.1

Patch release addressing issue where as.phylo was masked from ape.

#phylogram 2.0.0

Major release retaining the newick parsing and tree manipulation 
functions but migrating the k-mer counting and clustering functions 
to the **kmer** package. 

Package no longer needs compilation as all Rcpp functions 
are migrated. 

All dependencies are eliminated with the exception of **ape**.


# phylogram 1.0.0

Released on CRAN 2017-06-12.
#phylogram version 2.1.0

This minor release addresses some performance issues and features additional
details/examples in the vignette.

## Test environments

 * local ubuntu 16.04.4 x86_64-pc-linux-gnu; R version 3.4.4
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R 3.5.0
 * winbuilder R Under development (2018-06-21 r74929)

## R CMD check results

There were no ERRORs WARNINGs. There was one NOTE:

Possibly mis-spelled words in DESCRIPTION:
  Newick (20:40)
  Paradis (18:10)
  al (18:21)
  et (18:18)
  phylo (17:6)
  
These have been checked and are all fine.

## Downstream dependencies

kmer:  OK
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->

# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/phylogram/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/phylogram.git`
* Make sure to track progress upstream (i.e., on our version of `phylogram` at `ropensci/phylogram`) by doing `git remote add upstream https://github.com/ropensci/phylogram.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/phylogram`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
title: "phylogram: an R package for phylogenetic analysis with dendrograms"
author: "Shaun P. Wilkinson *^1,2^* and Simon K. Davy *^1^*"
output:
  html_document:
    keep_md: true
bibliography: phylogram.bib
vignette: >
  %\VignetteIndexEntry{Introduction to the phylogram package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE, message = FALSE, warning = FALSE}
#knitr::opts_chunk$set(out.width='750px', dpi=200)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", dpi=500, out.width='500px',
                      fig.path = 'figures/')
```

--------------------------------------------------------------------------------

*^1^* School of Biological Sciences, Victoria University of Wellington, P.O. Box 600, Wellington, New Zealand.

*^2^* Correspondence author. E-mail: shaunpwilkinson\@gmail.com

\  

## Abstract
The **phylogram** R package is a tool for for developing 
phylogenetic trees as deeply-nested lists known as "dendrogram" objects. 
It provides functions for conversion between dendrograms and 
"phylo" class objects, as well as several tools for command-line tree 
manipulation and import/export via Newick parenthetic text.
This improves accessibility to the comprehensive range of object-specific 
analytical and tree-visualization functions found across a wide array of 
bioinformatic R packages.
The **phylogram** package is released under the GPL-3 license, 
and is available for download from 
CRAN <https://CRAN.R-project.org/package=phylogram> and GitHub 
<https://github.com/ropensci/phylogram>.  

\  

## Introduction
The R environment continues to gain popularity as a platform for
bioinformatic analysis, due to the reproducible 
code-based workflow and the many powerful analytical tools 
available in a suite of open-source packages such as **ape** 
[@Paradis2004], **phangorn** [@Schliep2011] and **Phytools** [@Revell2012]. 
These packages typically employ a tree structure known as the 
"phylo" object, whose primary element is an integer matrix 
with one row for each edge in the graph, 
and two columns giving the indices of the connecting nodes. 
This is a versatile and memory-efficient representation suitable 
for most applications encountered by evolutionary biologists, 
and hence a comprehensive array of tools has been developed for 
parsing, manipulating, analyzing and visualizing trees in matrix format. 

An alternative tree structure is the "dendrogram" object, generated 
using the `as.dendrogram` function in the **stats** package [@RCoreTeam2015].
Rather than a matrix of edges, a dendrogram is a hierarchical list. 
These 'lists of lists' can be deeply nested, with the limit depending on 
the C stack size (settable *via* `options("expressions")`).
While less memory-efficient than matrix-based tree objects, 
a useful feature of the nested-list representation is its modularity, 
whereby the sub tree of a tree is itself a tree - 
a dendrogram within a dendrogram. 
This means that dendrogram objects are subsettable in the same 
way that standard lists are, which in addition to the inbuilt
editing functions such as `cut` and `merge`, 
facilitates intuitive command-line tree manipulation. 
An especially powerful feature of this object type is that tree 
editing operations can be carried out recursively 
using fast inbuilt functions in the "apply" family such as `dendrapply` 
and `lapply`. 

Each node of a dendrogram object has the following mandatory attributes: 

* "height" the position of the node along the vertical axis 
   (assuming the graph is orientated vertically)
* "midpoint" the horizontal distance of the node from the left-most member 
   of the sub tree (where the horizontal distance between adjacent leaves is 1 unit)
* "members" the number of terminal leaf nodes belonging to the node 
* "class" all nodes have the class attribute "dendrogram" 

Rather than lists, terminal leaf nodes are length-1 integer vectors whose
values correspond to the indices of the members in the set. 
The "members" attributes of leaf nodes is always 1, 
the "midpoint" attribute is 0, and they have two additional attributes:    

* "leaf" TRUE for terminal nodes (NULL otherwise)
* "label" an optional character string giving the name of the taxon or group

Aside from those listed above, users may attach other objects 
as attributes to the dendrogram nodes. For example, "label" attributes
can be attached to inner nodes, and users can specify plotting parameters
for each node by setting the attributes "nodePar" and "edgePar".
Any number of non-standard objects can also be attached as 
node attributes, which can be useful for storing additional metadata 
such as bootstrap values and taxonomic information.

The flexibility, modularity and intuitive structure of dendrogram 
objects are appealing to many users, particularly where highly dynamic 
tree structures are required for applications such as 
machine learning clustering and classification. 
There is also a large and growing number of resources for manipulating and 
plotting dendrograms in contributed packages such as 
**dendextend** [@Galili2015], and
hence functions enabling bi-directional conversion between "dendrogram"
and "phylo" class objects would expand the range of tools available for 
both object types. 
While conversion is currently possible using the "hclust" 
object as an intermediary, this object type does not support 
non-zero leaf node heights and hence is limited to ultrametric 
trees only.

\  

## The phylogram package
Here, we introduce **phylogram**, an R package for structuring 
evolutionary trees as deeply-nested lists and transforming trees between 
list- and matrix-type objects. 
The package also contains functions for importing and exporting dendrogram 
objects to and from parenthetic text, as well as several functions for 
manipulating trees in nested-list format.
These functions are detailed below with examples of their utility.

\  

### Importing and exporting trees
The Newick parenthetic text format 
[a.k.a. the New Hampshire format; @Felsenstein1986] 
is a universal phylogenetic tree 
representation that is compatible with most tree-editing software. 
The **phylogram** function `read.dendrogram` wraps the the Newick 
parser `read.tree` from the **ape** package [@Paradis2004], 
and converts the intermediate "phylo" object to a dendrogram.
This function supports weighted edges, 
labels with special meta-characters 
(enclosed in single quotation marks), comments 
(enclosed in square brackets; ignored by the parser), 
multifuricating nodes, and both rooted and unrooted trees.
Inner-node labels are also supported, and are 
attached as "label" attributes to non-leaf nodes.
Objects of class "dendrogram" can also be exported as 
Newick-style parenthetic text using the function 
`write.dendrogram`.

\  

#### Example 1: Import and export a tree from a Newick string
Consider the simple example of a tree with three members named 
"A", "B" and "C", where "B" and "C" are more closely related
to each other than either is to "A". 
An unweighted Newick string for this tree would be *(A,(B,C));*
This can be imported as a 
dendrogram object using the `read.dendrogram` function 
as follows:

```{r, fig.keep='none'}
library(phylogram)
x <- read.dendrogram(text = "(A,(B,C));")
plot(x, yaxt = "n")
```

```{r, echo=FALSE, fig.height=1.5, fig.width=3, fig.align='left', dpi=500}
par(mar = c(1, 1, 1, 1))
plot(x, yaxt = "n")
```

\ \ \ \ **Figure 1:** A simple dendrogram with three terminal leaf nodes

\  

The following command writes the object back to the console in 
Newick format without edge weights:

```{r}
write.dendrogram(x, edges = FALSE)
```

The syntax is similar when reading and writing text files, 
except that the `text` argument is replaced by `file`, 
and a valid file path is passed to the function. 

\  

### Converting tree objects
Dendrograms can be converted to "phylo" 
objects and *vice versa* using the `as.phylo.dendrogram` and 
`as.dendrogram.phylo` methods. 
Unlike functions that employ an "hclust" object as an
intermediary (e.g. `as.dendrogram(as.hclust(phy))`),
these methods retain all weighted edges and do not require trees 
to be ultrametric. 
This facilitates access to the comprehensive range of functions 
that are specific for either dendrograms or "phylo" objects in 
contributed packages such as **dendextend** [@Galili2015] and 
**ape** [@Paradis2004], respectively.
Note that other packages may employ the same function names, 
and hence the method dispatched may depend on the order 
in the which packages are loaded.
For this reason it may be safer to use the full function call 
(e.g. `phylogram::as.phylo.dendrogram(x)` and 
`phylogram::as.dendrogram.phylo(x)`) when using these methods.

\  

#### Example 2: Convert a "phylo" object to a dendrogram
A common application requiring conversion between
"phylo" and "dendrogram" objects involves plotting tanglegrams to
visualize incongruence between two phylogenetic trees. 
The **dendextend** package [@Galili2015] features the function `tanglegram` 
for versatile plotting of two distinct trees and 
indicating the discordant nodes using a series of non-parallel edges. 
However this function does not currently support non-ultrametric 
"phylo" objects.
In this example, two weighted neighbor-joining trees are generated 
from the left and right sections of the woodmouse alignment from the **ape**
package [@Paradis2004], and converted to dendrograms for visual comparison using `dendextend::tanglegram`.

```{r, message = FALSE, fig.height=4, fig.width=10, fig.align='left', out.width= '1000px'}
library(ape)
data(woodmouse)
## generate distance matrices for each section of the alignment
dist1 <- dist.dna(woodmouse[, 1:482])
dist2 <- dist.dna(woodmouse[, 483:965])
## build neighbor-joining trees
phy1 <- nj(dist1)
phy2 <- nj(dist2)
## root with No0912S as outgroup
phy1 <- root(phy1, "No0912S")
phy2 <- root(phy2, "No0912S")
## convert phylo objects to dendrograms
dnd1 <- as.dendrogram(phy1)
dnd2 <- as.dendrogram(phy2)
## rearrange in ladderized fashion
dnd1 <- ladder(dnd1)
dnd2 <- ladder(dnd2)
## plot the tanglegram
dndlist <- dendextend::dendlist(dnd1, dnd2)
dendextend::tanglegram(dndlist, fast = TRUE, margin_inner = 5)
```

\ \ \ \ **Figure 2:** Tanglegram showing incongruence between the left- and right-hand sections of the woodmouse alignment.

\  

### Tree editing/manipulation
The **phylogram** package features 
several additional functions to facilitate some of the more common 
manipulation operations.
Leaf nodes and internal branching nodes can be removed 
using the function `prune`, which identifies and 
recursively deletes nodes based on pattern 
matching of "label" attributes.
This is slower than the **ape** function `drop.tip`, but offers
the benefits of versatile string matching using regular expressions,
and the ability to remove inner nodes (and by extension all of 
their subnodes) that feature matching "label" attributes.
To aid visualization, the function `ladder` rearranges
the tree, sorting nodes by the number of members (analogous to the
`ladderize` function in the **ape** package). Another function 
aiding in tree visualization is `as.cladogram`, which
resets the "height" attributes of all terminal leaf nodes to zero 
and progressively resets the heights of the inner nodes
by single incremental units in a bottom-up fashion.
The function `reposition` scales the heights of all nodes in
a tree by a given constant (passed *via* the argument `shift`), 
and features the option to reset all node heights so that height of 
the farthest terminal leaf node from the root is zero (by specifying 
`shift = "reset"`). 
The function `remidpoint` recursively corrects all "midpoint", 
"members" and "leaf" attributes following manual editing of a tree or
while converting a nested list to a "dendrogram" object.

\  

#### Example 3: Building and manipulating dendrograms

The simple three-leaf dendrogram in Figure 1 can be created 
manually as follows:

```{r, fig.keep='none'}
x <- list(1, list(2, 3))
## attach "leaf" and "label" attributes to leaf nodes
attr(x[[1]], "leaf") <- TRUE
attr(x[[2]][[1]], "leaf") <- attr(x[[2]][[2]], "leaf") <- TRUE
attr(x[[1]], "label") <- "A"
attr(x[[2]][[1]], "label") <- "B"
attr(x[[2]][[2]], "label") <- "C"
## set "height" attributes for all nodes
attr(x, "height") <- 2
attr(x[[1]], "height") <- 0
attr(x[[2]], "height") <- 1
attr(x[[2]][[1]], "height") <- attr(x[[2]][[2]], "height") <- 0
## set "midpoints" attributes for all nodes
attr(x, "midpoint") <- 0.75
attr(x[[1]], "midpoint") <- 0
attr(x[[2]], "midpoint") <- 0.5
attr(x[[2]][[1]], "midpoint") <- attr(x[[2]][[2]], "midpoint") <- 0
## set "members" attributes for all nodes
attr(x, "members") <- 3
attr(x[[1]], "members") <- 1
attr(x[[2]], "members") <- 2
attr(x[[2]][[1]], "members") <- attr(x[[2]][[2]], "members") <- 1
## set class as "dendrogram" 
## Note that setting the class for the root node
## automatically sets the class of all nested subnodes
class(x) <- "dendrogram"
x
```

As demonstrated above, manually setting attributes on dendrogram
objects can be rather tedious, motivating the development of functions 
to automate the generation and manipulation of these tree 
structures.

This simple tree can be recreated more succinctly using the 
**phylogram** package functions as follows:

```{r}
x <- list(1, list(2, 3))
## recursively set class, midpoint, members and leaf attributes
x <- remidpoint(x)
## set incremental height attributes
x <- as.cladogram(x)
## set label attributes using dendrapply
set_label <- function(node){
  if(is.leaf(node)) attr(node, "label") <- LETTERS[node]
  return(node)
}
x <- dendrapply(x, set_label)
x
```

Similarly, dendrogram objects can be subset using either the `prune` 
function or standard list-subsetting syntax, again with the help of 
utility functions to recursively reset node attributes. 
The following code demonstrates one option for rearranging the tree with 
species A and B as sister taxa and C as the ancestor: 

```{r, fig.keep='none'}
## isolate root node (species C)
ancestor <- prune(x, pattern = "C", keep = TRUE) 
## alternative option using subset operator
ancestor <- x[[2]][[2]]
## create subtree without species C
subtree <- prune(x, pattern = "C")
## graft subtree onto root
x <- list(ancestor, subtree)
## set attributes as above
x <- as.cladogram(remidpoint(x))
## plot dendrogram
plot(x, yaxt = "n")
```

```{r, echo=FALSE, fig.height=1.5, fig.width=3, fig.align='left', dpi=500}
par(mar = c(1, 1, 1, 1))
plot(x, yaxt = "n")
```

\ \ \ \ **Figure 3:** Rearranged dendrogram with species C ancestral to A and B

\  

### Tree visualization 
Publication-quality trees can be generated from dendrogram objects 
using the **stats** plotting function `plot.dendrogram`, and the extensive
plotting functions available in dendrogram-enhancing packages such as 
**circlize** [@Gu2014] and **dendextend** [@Galili2015].
The latter also offers the facility to convert dendrograms to "ggdend" objects, 
for which many powerful 'grammar of graphics' plotting functions are available in 
the **ggplot2** [@Wickham2009] and **ggdendro** [@deVries2016] packages. 
Moreover, there are several advanced plotting options for "phylo" objects in 
the **ape** package [@Paradis2004], as well as the
Bioconductor package **ggtree** [@Guangchuang2017].
Given the extensive tree visualization options already available, 
we do not include any additional plotting functions in the **phylogram** package. 

\  

### Summary
The **phylogram** package offers a dendrogram parser for phylogenetic trees, several new tree-editing functions,
and a bridge between the "dendrogram"" and "phylo" object types that 
improves accessibility to the comprehensive number of object-specific 
functions found across a suite of contributed packages.
Future versions of the package will aim to further expand the range of input formats and object types available, 
thereby helping to integrate the wide variety of 
phylogenetic applications implemented in the R programming language.
This software is still under active development, and will
continue to be upgraded and expanded as new applications arise.
Bug reports and other feedback are welcomed and can be directed to the GitHub issues
page at <http://github.com/ropensci/phylogram/issues>, 
or the **phylogram** google group at <https://groups.google.com/group/phylogram>.

\  

## Acknowledgements 
This software was developed with funding from a Rutherford Foundation Postdoctoral 
Research Fellowship from the Royal Society of New Zealand. The authors declare no 
competing interests.

\  

## Author Contributions
Both authors conceived and designed the software. 
SPW wrote the package functions and documentation.
Both authors wrote the manuscript and gave final approval for publication.

\  

## References 
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylogram.R
\docType{package}
\name{phylogram}
\alias{phylogram}
\alias{phylogram-package}
\title{Dendrograms for evolutionary analysis.}
\description{
The \strong{phylogram} R package is a tool for for developing
phylogenetic trees as deeply-nested lists known as "dendrogram" objects.
It provides functions for conversion between "dendrogram" and
"phylo" class objects, as well as several tools for command-line tree
manipulation and import/export via Newick parenthetic text.
This improves accessibility to the comprehensive range of object-specific
analytical and tree-visualization functions found across a wide array of
bioinformatic R packages.
}
\section{Functions}{

A brief description of the primary \pkg{phylogram} functions are
  provided with links to their help pages below.
}

\section{File import/export}{

\itemize{
\item \code{\link{read.dendrogram}} reads a Newick
  parenthetic text string from a file or text connection
  and creates an object of class \code{"dendrogram"}
\item \code{\link{write.dendrogram}} outputs an object of class
  \code{"dendrogram"} to a text string or file in Newick
  format
}
}

\section{Object conversion}{

\itemize{
\item \code{\link{as.phylo.dendrogram}} converts a dendrogram to
  an object of class "phylo"
  \code{"dendrogram"}
\item \code{\link{as.dendrogram.phylo}} converts a "phylo" object
  to a dendrogram
}
}

\section{Tree editing and manipulation}{

\itemize{
\item \code{\link{prune}} remove branches from a \code{dendrogram} object
  based on regular expression pattern matching
\item \code{\link{ladder}} reorders the branches of a \code{dendrogram}
  object to aid visualization
\item \code{\link{remidpoint}} recursively sets "midpoint" and "members"
  attributes for a nested list/\code{dendrogram} object
\item \code{\link{reposition}} shifts a \code{dendrogram} object up or
  down (or sideways if plotted horizontally)
\item \code{\link{as.cladogram}} modifies the "height" attributes of the
  nodes such that all leaves terminate at zero
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ladder.R
\name{ladder}
\alias{ladder}
\title{Reorder tree branches in ladderized pattern.}
\usage{
ladder(x, decreasing = FALSE)
}
\arguments{
\item{x}{an object of class \code{"dendrogram"}.}

\item{decreasing}{logical indicating whether the tree should be
ladderized upwards or downwards. Defaults to FALSE (downwards).}
}
\value{
Returns an object of class \code{dendrogram}.
}
\description{
This function ladderizes the branches of a \code{dendrogram} object
  to aid in visual interpretation.
}
\details{
This function is the \code{dendrogram} analogue of the
  \code{\link[ape]{ladderize}} function in the \code{\link[ape]{ape}}
  package (Paradis et al 2004, 2012).
}
\examples{
  x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
  plot(x, horiz = TRUE)
  x <- ladder(x, decreasing = TRUE)
  plot(x, horiz = TRUE)
}
\references{
Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
  and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.

  Paradis E (2012) Analysis of Phylogenetics and Evolution with R
  (Second Edition). Springer, New York.
}
\seealso{
The \code{\link[ape]{ladderize}} function in the
  \code{\link[ape]{ape}} package performs a similar operation for objects
  of class \code{"phylo"}.
}
\author{
Shaun Wilkinson
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asdendrogram.R
\name{as.dendrogram.phylo}
\alias{as.dendrogram.phylo}
\title{Convert a "phylo" object to a dendrogram.}
\usage{
\method{as.dendrogram}{phylo}(object, ...)
}
\arguments{
\item{object}{an object of class "phylo".}

\item{...}{further arguments to be passed between methods.}
}
\value{
an object of class "dendrogram".
}
\description{
This function converts a "phylo" object (Paradis et al 2004)
  to a dendrogram.
}
\examples{
  newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  x <- read.dendrogram(text = newick)
  y <- as.phylo(x)
  z <- as.dendrogram(y)
  identical(x, z)
}
\references{
Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
  and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.

  Paradis E (2008) Definition of Formats for Coding Phylogenetic Trees in R.
  \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}

  Paradis E (2012) Analysis of Phylogenetics and Evolution with R
  (Second Edition). Springer, New York.
}
\author{
Shaun Wilkinson
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{as.cladogram}
\alias{as.cladogram}
\title{Apply unweighted branch lengths.}
\usage{
as.cladogram(x)
}
\arguments{
\item{x}{an object of class \code{"dendrogram"}.}
}
\value{
an object of class \code{"dendrogram"}.
}
\description{
This function sets the 'height' attributes of all leaf nodes to zero and
  progressively resets the heights of the inner nodes by single incremental
  units in a bottom-up fashion.
}
\examples{
  x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
  plot(x, horiz = TRUE)
  x <- as.cladogram(x)
  plot(x, horiz = TRUE)
}
\author{
Shaun Wilkinson
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asphylo.R
\name{as.phylo.dendrogram}
\alias{as.phylo.dendrogram}
\title{Convert a dendrogram to a "phylo" object.}
\usage{
\method{as.phylo}{dendrogram}(x, ...)
}
\arguments{
\item{x}{a dendrogram.}

\item{...}{further arguments to be passed between methods.}
}
\value{
an object of class "phylo".
}
\description{
This function converts a dendrogram into an object of class "phylo"
  (see Paradis et al 2004).
}
\examples{
  newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  x <- read.dendrogram(text = newick)
  y <- as.phylo(x)
  z <- as.dendrogram(y)
}
\references{
Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
  and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.

  Paradis E (2008) Definition of Formats for Coding Phylogenetic Trees in R.
  \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}

  Paradis E (2012) Analysis of Phylogenetics and Evolution with R
  (Second Edition). Springer, New York.
}
\author{
Shaun Wilkinson
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write.R
\name{write.dendrogram}
\alias{write.dendrogram}
\title{Write a dendrogram object to parenthetic text.}
\usage{
write.dendrogram(x, file = "", append = FALSE, edges = TRUE, ...)
}
\arguments{
\item{x}{an object of class \code{"dendrogram"}.}

\item{file}{a character string naming a file or connection to write the
output to. If no file path is specified or \code{file = ""} the result
is printed to the console.}

\item{append}{logical indicating whether the output should be
appended to the file. If \code{append = FALSE} the contents of the
file will be overwritten (the default setting).}

\item{edges}{logical indicating whether edge weights should be
included in the output string.}

\item{...}{further arguments to be passed to \code{format}. Used to
specify the numbering style of the edge weights (if edges = TRUE).}
}
\description{
This function exports a dendrogram object as a Newick/New Hampshire
  text string.
}
\examples{
  newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  x <- read.dendrogram(text = newick)
  write.dendrogram(x, edges = TRUE)
}
\references{
Felsenstein J (1986) The Newick tree format.
  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}

  Olsen G (1990) Interpretation of the "Newick's 8:45" tree format standard.
  \url{http://evolution.genetics.washington.edu/phylip/newick_doc.html}
}
\seealso{
\code{\link{read.dendrogram}} to parse a \code{"dendrogram"}
  object from a text file.
  The \code{\link[ape]{write.tree}} function in the \code{\link[ape]{ape}}
  package performs a similar operation for \code{"phylo"}
  and \code{"multiPhylo"} objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prune.R
\name{prune}
\alias{prune}
\title{Remove tree nodes by regular expression pattern matching.}
\usage{
prune(tree, pattern, keep = FALSE, ...)
}
\arguments{
\item{tree}{an object of class \code{"dendrogram"}.}

\item{pattern}{a regular expression.}

\item{keep}{logical indicating whether the nodes whose labels match
the regular expression provided in "pattern" should be
kept and the remainder discarded. Defaults to FALSE.
Note that nodes without "label" attributes are ignored.}

\item{...}{further arguments to be passed to \code{grepl} and \code{gsub}.}
}
\value{
Returns an object of class \code{"dendrogram"}.
}
\description{
\code{"prune"} takes an object of class \code{"dendrogram"} and
  removes all branches whose branch labels match a given regular
  expression.
}
\details{
This function recursively tests the "label" attribute of each
  dendrogram node (including non-leaf inner nodes if applicable) for
  the specified pattern, removing those that register a positive hit.
  Note that positive matching inner nodes are removed along with all of
  their sub-nodes, regardless of whether the "label" attributes of the
  sub-nodes match the pattern.
}
\examples{
  x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
  plot(x, horiz = TRUE)
  x <- prune(x, pattern = "^A$")
  plot(x, horiz = TRUE)
}
\seealso{
The \code{\link[ape]{drop.tip}} function in the
  \code{\link[ape]{ape}} package performs a similar operation for objects
  of class \code{"phylo"}. See \code{\link{regex}} for help with
  compiling regular expressions.
}
\author{
Shaun Wilkinson
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{reposition}
\alias{reposition}
\title{Reset dendrogram height attributes.}
\usage{
reposition(x, shift = "reset")
}
\arguments{
\item{x}{an object of class \code{"dendrogram"}.}

\item{shift}{either the character string "reset" (shift the graph so that
the height of the farthest leaf from the root is zero), or a numeric value
giving the amount to shift the graph along the primary axis.}
}
\value{
Returns an object of class \code{"dendrogram"}.
}
\description{
\code{reposition} is a helper function used for manually creating
  \code{"dendrogram"} objects from nested lists. The function
  recursively reassigns the 'height' attributes at each node by
  a given constant.
}
\examples{
  x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
  plot(x, horiz = TRUE)
  x <- reposition(x)
  plot(x, horiz = TRUE)
}
\author{
Shaun Wilkinson
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{remidpoint}
\alias{remidpoint}
\title{Set dendrogram attributes for a nested list.}
\usage{
remidpoint(x)
}
\arguments{
\item{x}{a nested list, possibly of class \code{"dendrogram"}}
}
\value{
returns a nested list, or an object of class \code{"dendrogram"}
  depending on the class of the input object.
}
\description{
\code{remidpoint} is a helper function used for manually creating
  \code{"dendrogram"} objects from nested lists. The function
  recursively assigns the necessary 'midpoint', 'members', and
  'leaf' attributes at each node.
}
\examples{
  ## manually create a small dendrogram with three members, A, B and C
  x <- list("A", list("B", "C"))
  attr(x[[1]], "leaf") <- TRUE
  attr(x[[2]][[1]], "leaf") <- TRUE
  attr(x[[2]][[2]], "leaf") <- TRUE
  attr(x[[1]], "label") <- "A"
  attr(x[[2]][[1]], "label") <- "B"
  attr(x[[2]][[2]], "label") <- "C"
  attr(x, "height") <- 1
  attr(x[[1]], "height") <- 0
  attr(x[[2]], "height") <- 0.5
  attr(x[[2]][[1]], "height") <- 0
  attr(x[[2]][[2]], "height") <- 0
  x <- remidpoint(x)
  class(x) <- "dendrogram"
  plot(x, horiz = TRUE)
}
\author{
Shaun Wilkinson
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read.dendrogram}
\alias{read.dendrogram}
\title{Read a dendrogram from parenthetic text.}
\usage{
read.dendrogram(file = "", text = NULL, ...)
}
\arguments{
\item{file}{character string giving a valid path to the file from
which to read the data.}

\item{text}{optional character string in lieu of a "file" argument.
If a text argument is provided instead of a file path, the data
are read via a text connection.}

\item{...}{further arguments to be passed to
\code{\link[ape]{read.tree}} (which may then be passed on to
\code{scan}).}
}
\value{
an object of class \code{"dendrogram"}.
}
\description{
This function wraps the \code{\link[ape]{read.tree}} parser from the
  \code{\link[ape]{ape}} package to read a phylogenetic tree from
  parenthetic text in the Newick/New Hampshire format, and
  converts it to object of class "dendrogram".
}
\examples{
  x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
  plot(x, horiz = TRUE)
}
\references{
Felsenstein J (1986) The Newick tree format.
  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}

  Olsen G (1990) Interpretation of the "Newick's 8:45" tree format standard.
  \url{http://evolution.genetics.washington.edu/phylip/newick_doc.html}

  Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
  and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.

  Paradis E (2008) Definition of Formats for Coding Phylogenetic Trees in R.
  \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}

  Paradis E (2012) Analysis of Phylogenetics and Evolution with R
  (Second Edition). Springer, New York.
}
\seealso{
\code{\link{write.dendrogram}} writes an object of
  class \code{"dendrogram"} to a Newick text string.
  The \code{\link[ape]{read.tree}} function in the
  \code{\link[ape]{ape}} package parses objects
  of class \code{"phylo"} and \code{"multiPhylo"}.
}
\author{
Shaun Wilkinson
}
