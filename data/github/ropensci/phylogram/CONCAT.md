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
