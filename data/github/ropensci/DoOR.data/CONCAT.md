DoOR.data
=========
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.375618.svg)](http://dx.doi.org/10.5281/zenodo.375618)
[![](https://badges.ropensci.org/35_status.svg)](https://github.com/ropensci/onboarding/issues/35)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/DoOR.data.svg?branch=master)](https://travis-ci.org/ropensci/DoOR.data)


R package containing the data for the Database of Odor Responses. The corresponding functions package can be found at [https://github.com/ropensci/DoOR.functions](https://github.com/ropensci/DoOR.functions).

## The DoOR Project
Find more information, precompiled R-packages and a interactive web-version of the DoOR-Database at: **[http://neuro.uni.kn/DoOR](http://neuro.uni.kn/DoOR)**


## Install
Either download a packaged version or install _via_ `devtools`:
```{r}
# install devtools
install.packages("devtools")
library(devtools)

# install DoOR.data 2.0.1
install_github("ropensci/DoOR.data", ref="v2.0.1")

# or install the latest version available on Github
install_github("ropensci/DoOR.data")
```

## Publications
DoOR was first published in 2010, the OpenAccess publication is available from
[http://chemse.oxfordjournals.org/content/35/7/551](http://chemse.oxfordjournals.org/content/35/7/551 "10.1093/chemse/bjq042")

An OpenAccess publication regarding the comprehensive update to **DoOR version 2.0** is available from
A preprint of manuscript related to DoOR 2.0 can be found at bioRxiv: [http://www.nature.com/articles/srep21841](http://www.nature.com/articles/srep21841)

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# DoOR.data v2.0.2
- added unit testing


# DoOR.data v2.0.1
- replaced '.' with '\_' in data names (#19)
- converted function names from camelCase to snake_case
- added warning message that data is loaded into current workspace and forced user feedback (#17)
- wrote documentation for datasets (#15)
- revert to store data as text rather than binary to allow change tracking (#14)

# DoOR.data v2.0.0
- A comprehensive update to data and functions of the DoOR project. Please see the publication for details: http://doi.org/10.1038/srep21841
- added several new data from literature and updated existing datasets with original data from authors
- introduced InChIKeys as main chemical identifier

## Bugfixes
- renamed ab3B \> Or85b in DoOR.mappings (#7)
- set SFR == 0 for all Ca²⁺ imaging studies (#4)
- set SFR == 0 for all studies that subtracted SFR but did not report it (#3)


# DoOR.data v1.0.2
- several bugfixes

# DoOR.data v1.0
- initial release as published in: Integrating heterogeneous odor response data into a common response model: A DoOR to the complete olfactome. ChemSenses 35, 551–63. http://doi.org/10.1093/chemse/bjq042
