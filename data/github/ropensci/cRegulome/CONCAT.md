[![Build Status](https://travis-ci.org/ropensci/cRegulome.svg?branch=master)](https://travis-ci.org/ropensci/cRegulome)
[![codecov](https://codecov.io/gh/ropensci/cRegulome/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/cRegulome)
[![Build status](https://ci.appveyor.com/api/projects/status/gcmojtcsyt7rcwtk?svg=true)](https://ci.appveyor.com/project/ropensci/cregulome)
[![](https://badges.ropensci.org/149_status.svg)](https://github.com/ropensci/onboarding/issues/149)  
[![CRAN version](https://img.shields.io/badge/CRAN-v0.3.0-blue.svg)](https://CRAN.R-project.org/package=cRegulome) 
![downloads](https://cranlogs.r-pkg.org/badges/grand-total/cRegulome)  

# cRegulome  
## Overview  
Transcription factors and microRNAs are important for regulating the gene
expression in normal physiology and pathological conditions. Many
bioinformatics tools were built to predict and identify transcription
factors and microRNA targets and their role in development of diseases
including cancers. The availability of public access high-throughput data
allowed for data-driven predictions and discoveries.
Here, we build on some of these tools and integrative analyses and provide a
tool to access, manage and visualize data from open source databases.
cRegulome provides a programmatic access to the regulome (microRNA and
transcription factor) correlations with target genes in cancer. The package
obtains a local instance of 
[Cistrome Cancer](http://cistrome.org/CistromeCancer/) and 
[miRCancerdb](https://mahshaaban.shinyapps.io/miRCancerdb/) databases and
provides classes and methods to interact with and visualize the correlation
data.  

## What is cRegulome used for?  
cRegulome provides programmatic access to regulome-gene correlation data in 
cancer from different data sources. Researches who are interested in studying 
the role of microRNAs and transcription factors in cancer can use this package 
to construct a small or large scale queries to answer different questions:  

* Which microRNAs and/or transcription factors are associated with a particular
set of genes?  
* What different regulation patterns a microRNA or a transcription factor can 
take in different types of cancer?  
* For a given set of regulatory elements, which genes are likely to be 
regulated by these elements in a certain type of cancer?  

In addition, cRegulome can be used with other R packages like `igraph` to 
study the co-regulation networks in different types of cancer.  
    

## Getting started  
To get starting with cRegulome we show a very quick example. We first start
by downloading a small test database file, make a simple query and convert
the output to a cRegulome object to print and visualize.  

```r
# install the package from CRAN
install.packages('cRegulome')
```

```r
# install the development version from github
devtools::install_github('ropensci/cRegulome')

# install the development version and build vignette from github 
devtools::install_github('ropensci/cRegulome', build_vignettes = TRUE)
```

```{r load_libraries}
# load required libraries
library(cRegulome)
library(RSQLite)
library(ggplot2)
```

```r
if(!file.exists('cRegulome.db')) {
    get_db(test = TRUE)
}

# connect to the db file
conn <- dbConnect(SQLite(), 'cRegulome.db')
```

Or access the same test set file from the package directly  

```r
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- dbConnect(SQLite(), fl)
```

```r
# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = 'hsa-let-7g',
               study = 'STES',
               min_abs_cor = .3,
               max_num = 5)

# make a cmicroRNA object   
ob <- cmicroRNA(dat)
```

```r
# print object
ob
```

## Documentation

```r
browseVignettes("cRegulome")
```
Alternatively, the vingettes can be found online, [case_study](http://rpubs.com/MahShaaban/cRegulome1) 
and [using_cRegulome](http://rpubs.com/MahShaaban/cRegulome2).

## Citation  

```r
citation("cRegulome")
```

[![](http://www.ropensci.org/public_images/github_footer.png)](http://ropensci.org)

# cRegulome 0.3.2
    
    - Restored on CRAN
    
# cRegulome 0.3.1
    
    - Fixed a CRAN test failure due to changes in a dep pkg

# cRegulome 0.3.0
    
    - Bug fix: since 0.2.0 the argument targets_only did not work properly.
    The bug is fixed and tested in this release.
    - Added the option directed to cor_igraph which allows for constructing
    a directed graph when desired.
    - Added the option to limit the query output of get_tf and get_mir to 
    a predifined set of genes.
  
# cRegulome 0.2.0
    
    - Reduced code dependencies
    - Improved code performance

# cRegulome 0.1.1

    - fix installing in default library tree
    
# cRegulome 0.1.0

    - cRegulome v0.1.0 (2018-02-08) Approved by rOpenSci
    - On CRAN

# cRegulome 0.99.0

    - cRegulome v0.99.0 (2017-09-06) Submit to rOpenSci

# All contributions are most welcomed  

I'd be glad to recieve any comments or ideas to help the package forward.  

## Bugs  

* To report a bug please use the [issue](https://github.com/MahShaaban/cRegulome/issues) page on github  

## Code contributions  

* Fork this repo to your github account  
* Clone the repo to your machine and make changes  
* Push to your account  
* Submit a [pull](https://github.com/MahShaaban/cRegulome/pulls) request at this repo  
    
## Email: [mahmoud.s.fahmy@students.kasralainy.egu.eg](mahmoud.s.fahmy@students.kasralainy.egu.eg)  

## Thank you for contributing!  ## Test environments
* local OS X install, R 3.4.3
* ubuntu 12.04 (on travis-ci), R 3.4.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

## Submissions

* First submission v0.1.0
* Updated the title and description fields
* Updated the description with quoted software names
* Default 'get_db' to write to a temporary directory

---

## Fixed

* Fix installing in default library tree

--- 

## v0.2.0

* Submission v0.2.0
* Reduced code dependencies
* Improved code performance

## v0.3.0

* Submission v0.3.0
* A bug fix
* Minor changes to get_tf, get_mir and cor_igraph
