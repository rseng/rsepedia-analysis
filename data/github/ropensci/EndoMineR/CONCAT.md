
[![Build
Status](https://travis-ci.org/sebastiz/EndoMineR.svg?branch=master)](https://travis-ci.org/sebastiz/EndoMineR)
[![ropensci](https://badges.ropensci.org/153_status.svg)](https://github.com/ropensci/onboarding/issues/153)
[![Coverage
status](https://codecov.io/gh/sebastiz/EndoMineR/branch/master/graph/badge.svg)](https://codecov.io/github/sebastiz/EndoMineR?branch=master)

<!-- README.md is generated from README.Rmd. Please edit that file -->

    ## here() starts at /home/rstudio/EndoMineR

This package has undergone a major revision to make it much more user
friendly. THe documentation has been updated to reflect this. I am
always happy to hear of any feedback, positive and negative.

## **Aims of EndoMineR**

The goal of EndoMineR is to extract as much information as possible from
free or semi-structured endoscopy reports and their associated pathology
specimens. A full tutorial can be found
[here](https://docs.ropensci.org/EndoMineR/articles/EndoMineR.html)

## Installation

You can install EndoMineR from github with:

``` r
# install.packages("devtools")
devtools::install_github("ropenSci/EndoMineR")
```

If you dont have access to github, then download the zip and change the
working dirctory to the place you have downloaded it, then do

``` r
setwd("C:/Users/Desktop/")

#On windows you cand cd to change the directory or us pushd to create a temporary directory indtead of cd and then setwd to the temporary directory
unzip("EndoMineR.zip")
file.rename("EndoMineR.zip-master", "EndoMineR.zip")
shell("R CMD build EndoMineR.zip")

#Then install the resulting tarball with:

install.packages("EndoMineR_0.2.0.9000.tar.gz", repos = NULL)
```

### How to contribute

Contributions to this project are most welcome. There are just a few
small guidelines you need to follow.

#### Submitting a patch

It’s generally best to start by opening a new issue describing the bug
or feature you’re intending to fix. Even if you think it’s relatively
minor, it’s helpful to know what people are working on. Mention in the
initial issue that you are planning to work on that bug or feature so
that it can be assigned to you.

Follow the normal process of forking the project, and setup a new branch
to work in. It’s important that each group of changes be done in
separate branches in order to ensure that a pull request only includes
the commits related to that bug or feature.

The best way to ensure your code is properly formatted is to use lint.
Various packages in R provide this.

Any significant changes should almost always be accompanied by tests.
The project already has good test coverage, so look at some of the
existing tests if you’re unsure how to go about it.

Do your best to have well-formed commit messages for each change. This
provides consistency throughout the project, and ensures that commit
messages are able to be formatted properly by various git tools.

Finally, push the commits to your fork and submit a pull request.
Please, remember to rebase properly in order to maintain a clean, linear
git
history.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# EndoMineR 2.0.0.9000

* EndoMineR has evolved both as a project and as a code base
  It is now fully Git compliant allowing real collaboration
  The functions have been hugely simplified. The text preparation for example is not just one function 'textPrep' which means the user can concentrate more on the higher order function.
  The whole package has now becaome more lexicon-centric meaning that the term mapping, spell correction etc is all guided by the prepackaged lexicons. The lexicons act as a dictionary for the package to look up (eg 'RFA' and 'HALO' is mapped to 'Radiofrequency ablation'). The lexions are likely to grow and further work is underway to expand different sub-speciality lexicons useing the United Medical Language System.
  Several further modules are being developed for subspecialties. The Barrett's functions have grown as have the polyp functions.
  Further work is underwy for inflammatory bowel disease and ERCP.
  The data visualisation functions have become more streamlined so that all graphs adhere to publication ready specifications
  A template project has been included so the user can see how to run a script from question to abstract creation. The template project also using standardised project infrastructure so the user has a head start in using best practice for project development.
  And there is so much more to come....



# Contributing to `EndoMineR`

Thank you for any and all contributions! Following these guidelines will help streamline the process of contributing and make sure that we're all on the same page. While we ask that you read this guide and follow it to the best of your abilities, we welcome contributions from all, regardless of your level of experience.



# Types of contributions 

All contributions are welcome, even just a comment or a slap on the back to say well done. Examples of contributions include:
  
- Identify areas for future development ([open an Issue](https://github.com/sebastiz/EndoMineR/issues))
- Identify issues/bugs ([open an Issue] (https://github.com/sebastiz/EndoMineR/issues))
- Write tutorials/vignettes ([open a Pull Request](https://github.com/sebastiz/EndoMineR/pulls) to contribute to the ones here, or make your own elsewhere and send us a link)
- Add functionality ([open a Pull Request](https://github.com/sebastiz/EndoMineR/pulls))
- Fix bugs ([open a Pull Request](https://github.com/sebastiz/EndoMineR/pulls))

# New to GitHub?

Getting ready to make your first contribution? Here are a couple of tutorials you may wish to check out:
  
  - [Tutorial for first-timers](https://github.com/Roshanjossey/first-contributions)
- [How to contribute (in-depth lessons)](https://egghead.io/series/how-to-contribute-to-an-open-source-project-on-github)
- [GitHub on setup](https://help.github.com/articles/set-up-git)
- [GitHub on pull requests](https://help.github.com/articles/using-pull-requests/).)


# How to contribute code

- Fork the repository
- Clone the repository from GitHub to your computer e.g,. `git clone https://github.com/ropensci/EndoMineR.git`
- Make sure to track progress upstream (i.e., on our version of `EndoMineR` at `ropensci/EndoMineR`)
- `git remote add upstream https://github.com/ropensci/EndoMineR.git`
- Before making changes make sure to pull changes in from upstream with `git pull upstream`
- Make your changes
- For changes beyond minor typos, add an item to NEWS.md describing the changes and add yourself to the DESCRIPTION file as a contributor
- Push to your GitHub account


# Code formatting

- In general follow the convention of <http://r-pkgs.had.co.nz/r.html#style> (snake_case functions and argument names, etc.)
- Where there is conflict, default to the style of `EndoMineR`
- Use explicit package imports (i.e. package_name::package_function) and avoid @import if at all possible---
title: 'EndoMineR for the extraction of endoscopic and associated pathology data from medical reports'
tags:
  - example
  - tags
  - for the paper
authors:
  - name: Sebastian S Zeki
    orcid: 0000-0003-1673-2663
    affiliation: "1"

affiliations:
  - name: Department of Gastroenterology, St Thomas' Hospital, Westminster Bridge Bridge Road, London SE1 7EH
    index: 1
date: 25th April 2018
bibliography: paper.bib
---

# Summary


Medical data is increasingly kept in an electronic format worldwide [@Bretthauer2016Reporting]. This serves many purposes including more efficient storage, distribution and accessibility of patient-focussed data. As important is the ability to analyse healthcare data for to optimize resource deployment and usage.  The tools for the analysis are often statistical and rely on the provision of ‘clean’ datasets before this can be done. ‘Cleaning’ a dataset is often the most difficult aspect of any data analysis and involves the provision of meaningful and well-formatted data so that the interpretation of the analysis is not subject to the doubts of the data quality. 

The British Society of Gastroenterology recommends that all endoscopic data is kept in an electronic format particularly to facilitate audit and maintain standards through the Global Rating Scale (GRS) [@Stebbing2011quality]. The endoscopic dataset is however only part of the patient’s story as many aspects of a patient’s gastroenterological care depend on the results of histopathological analysis of tissue taken during the examination. Pathology results are often available many days after the endoscopic result and usually stored in a separate data repository, although this may change with the arrival of an encompassing electronic patient record. 
Regardless of the method of storage, it is often difficult to associate the particular  histopathological result with an endoscopic result. Further, even if the two data sets can be merged, a problem occurs in the isolation of various parts of each report such that each part can be individually analysed.  Examples include the isolation of who the endoscopist was or the presence of dysplasia within a histopathology report. This is all the more difficult if the report is unstructured or partially structured free text. 

However if this can be done then many downstream analyses which benefit individual patients as well as the department, can be automated and include more complex analyses to determine follow-up regimes or endoscopic –pathologic lesion recognition performance.

The EndoMineR package provides a comprehensive way to extract information from natural language endoscopy ann pathology reports as well as merging the two datasets so that pathology specimens are relevant to the endoscopy they came from. Furthermore the package also provides functions for the following types of analysis of endoscopic and pathological datasets:

 + 1. Patient surveillance. Examples including determining when patients should return for surveillance and who is overdue.
 + 2. Patient tracking. -Examples include determining the length of time since the last endoscopy, as well as aggregate functions such as finding how many endoscopies of a certain type have been done and predicting future burden.
 + 3. Patient flow - determining the kinds of endoscopies an individual patient may get over time eg for ablation of Barrett's oesophagus.
 + 4. Quality of endoscopy and pathology reporting- Determining whether endoscopy quality is being maintained using some of the Global Rating scale metrics. Also making sure the pathology reports are complete.
 + 5. Diagnostic yield. Examples include determination of detection of dysplasia and cancer by endoscopist as a measure of lesion quality.

 It is the purpose of the package to create a unified process for merging of endoscopy reports with their associated pathology reports and to allow the extraction and tidying of commonly need data. Furthermore the package has methods for the analysis of the data in areas that are commonly required for high quality endoscopic services. This includes methods to track patients who need endoscopic surveillance, methods to determine endoscopic quality and disease detection rates. Also included are methods to assess patient flow through different types of endoscopy and to predict future usage of certain endoscopic techniques.
 
The package is in the process of having each analysis function validated and functions some validation has been submitted in abstract form to gastroenterological societies. 


# References
# Contributing to EndoMineR

This outlines how to propose a change to EndoMineR. For more detailed
info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib).

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  New code should follow the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the current
development version header describing the changes made followed by your GitHub
username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to
abide by its terms.

### See tidyverse [development contributing guide](https://rstd.io/tidy-contrib) for further details.
