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
- initial release as published in: Integrating heterogeneous odor response data into a common response model: A DoOR to the complete olfactome. ChemSenses 35, 551–63. http://doi.org/10.1093/chemse/bjq042