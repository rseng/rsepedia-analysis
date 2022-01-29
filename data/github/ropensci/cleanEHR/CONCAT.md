[![Build Status](https://travis-ci.org/ropensci/cleanEHR.svg?branch=master)](https://travis-ci.org/ropensci/cleanEHR)
[![codecov](https://codecov.io/gh/ropensci/cleanEHR/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/cleanEHR)
[![AUR](https://img.shields.io/aur/license/yaourt.svg)]()
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/cleanEHR)](https://cran.r-project.org/package=cleanEHR)
[![Downloads](http://cranlogs.r-pkg.org/badges/cleanEHR?color=brightgreen)](http://www.r-pkg.org/pkg/cleanEHR)
[![DOI](https://zenodo.org/badge/50511003.svg)](https://zenodo.org/badge/latestdoi/50511003)
[![](https://badges.ropensci.org/145_status.svg)](https://github.com/ropensci/onboarding/issues/102)


cleanEHR is an electronic health care record (EHR) data cleaning and
processing platform, which works with the Critical Care Health Informatics
Collaborative's data set. The purpose of the project is to enable researchers
to answer clinical questions that are important to patients, but which are
normally too difficult because data is unstandardised, siloed, and
inaccessible. 

Since 2014 data from the critical care units at Cambridge, Guys/Kings/St
Thomas', Imperial, Oxford, and University College London has been extracted and
stored securely in a standardised format. 

These data are crucially needed by healthcare professionals for the delivery
and continuity of care; by administrators for audit, planning and service
improvement; and by academic and industry researchers for the translation of
scientific progress into patient benefit. Through this process, CC-HIC can
improve patient outcomes, reduce the costs of care, and accelerate the pace of
translational health research. 

The physical database is held at the [UCL
IDHS](http://www.ucl.ac.uk/isd/itforslms/services/handling-sens-data/tech-soln)
within the Information Services Division of University College London (UCL).
UCL manage and ensure that the database and the surrounding governance
structures are appropriate for holding identifiable and sensitive NHS data. The
safe haven is compliant to NHS Information Governance Toolkit Level 2 and
operates to the ISO 27001, the Safe Haven already holds identifiable, sensitive
NHS data for secondary purpose. The Critical Care HIC management group will
maintain oversight and be the point of contact for researchers wishing to
access data.

The security put in place to ensure the safety of this resource inevitably
creates challenges for the researchers. We have therefore created a three sided
tool to make CCHIC _research ready_.

- this shared code library
- an anonymised development data set 
- a virtual machine for simulating work within the safe haven

You request access to the anonymised toy dataset from
[here](https://form.jotformeu.com/80383154562355)

## Required packages
* R (>= 3.1.0),
* XML,
* data.table,
* yaml,
* pander,
* Rcpp,
* methods


## How to install the R package
### From CRAN to get the last stable version

```
install.packages("cleanEHR")
```

### From Github to install the latest development version.
```
install.packages("devtools")
devtools::install_github("CC-HIC/cleanEHR")
```
## Vignette

* Introduction to CCHIC critical care data [here](https://docs.ropensci.org/cleanEHR/articles/cchic_overview.html)
* Data cleaning and wrangling with cleanEHR [here](https://docs.ropensci.org/cleanEHR/articles/data_clean.html)


## How to contribute
The cleanEHR package is currently under development. We wish you will find our
software tools useful to your research project. If you have any question about
the code or the data, please just raise an issue ticket on Github or contact us
via email (s.shi@ucl.ac.uk). As cleanEHR is an open source and community
project, we also welcome contributions of any kind. We would like to make
cleanEHR the platform of collaboration critical care data analysis. If you have
any idea or comment, please feel free to raise the issue ticket [here](https://github.com/ropensci/cleanEHR/issues). We would
also encourage everyone to integrate their code into the cleanEHR repository to
benefit the entire community. Please feel free to contact us when you wish to
contribute your code.


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Critical care data processing tools"
tags:
  - critical care
  - electronic health record
  - data manipulation
  - data pipeline
  - data cleaning
  - anonymisation
authors:
  - name: Sinan Shi
    orcid: 0000-0001-6935-0429
    affiliation: 1
  - name: David Pérez-Suárez
    orcid: 0000-0003-0784-6909
    affiliation: 1
  - name: Steve Harris
    orcid: 0000-0002-4982-1374
    affiliation: 2
  - name: Niall MacCallum
    orcid: 0000-0001-6168-1506
    affiliation: 2
  - name: David Brealey
    orcid: 0000-0002-1982-3379
    affiliation: 2
  - name: Mervyn Singer
    orcid: 0000-0002-1042-6350
    affiliation: 2
  - name: James Hetherington
    orcid: 0000-0001-6993-0319
    affiliation: 1
affiliations:
  - name: Research Software Development Group, Research IT Services, University College London
    index: 1
  - name: University College London Hospitals
    index: 2
date: 16 December 2016
bibliography: paper.bib
---

# Summary

cleanEHR [@cleanEHR] is a data cleaning and wrangling platform which works with
the Critical Care Health Informatics Collaborative (CCHIC) database.  CCHIC
collects and gathers high resolution longitudinal patient record from critical
care units at Cambridge, Guys/Kings/St Thomas', Imperial, Oxford, UCL
Hospitals. 

The increased adoption of high resolution longitudinal EHRs has created novel
opportunities for researchers, clinicians and data scientists to access large,
enriched patient databases [@icnarc] [@mimic]. The purpose of cleanEHR is to
enable researchers to answer clinical questions that are important to patients.
cleanEHR is a solution to address various data reliability and accessibility
problems as well. It provides a platform that enables data manipulation,
transformation, reduction, cleaning and validation with a friendly user
interface which empowers non-programmers to conduct basic data analysis by
simply writing a human-readable configuration file.

# High resolution longitudinal EHR: CCHIC

CCHIC database has in total collected 22,628 admissions (18,074 unique
patients) from 2014 to 2017. It contains 119 million data points (mean 6626
data points per patient). The recruited patients have an age range from 18 to
116 years old. Physiological, laboratory, drugs and nursing information are the
longitudinal data recorded during a patient's stay of the ICU.  The full list
of longitudinal data collected by CCHIC is listed below.

![List of CCHIC longitudinal data fields](graph/item_ref_time.png)

![Selected data fields of an admission](graph/episode_graph.pdf)

# Data cleaning and wrangling

Data of this kind, though provides vast information, often faces two
main issues, a) low data quality, b) low accessibility due to the complexity.
We proposed a workflow, which has been incorporated in the cleanEHR package, to address
the these issues. The highlight of this workflow includes the following,

* A table structure (ccTable) for data manipulation.
* Configuration file for researchers without technical knowledge to select and clean 
the data. The data cleaning includes various filters and data interpolation (impute)
function.

For detail description of the functions and examples, please see the manual and
the vignettes of cleanEHR [@cleanEHR]

![An example of filtering abnormal heart rate values by range](graph/range_filter.png)

# Reference
