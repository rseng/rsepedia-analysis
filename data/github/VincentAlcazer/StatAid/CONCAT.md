---
title: 'StatAid: An R package with a graphical user interface for data analysis'
authors:
- affiliation: 1, 2
  name: Vincent Alcazer
  orcid: 0000-0003-1843-6286
date: "20 October 2020"
bibliography: paper.bib
tags:
- R
- Data analysis
- Medicine
- Science
- Survival analysis
affiliations:
- index: 1
  name: Cancer Research Center of Lyon, INSERM U1052, Lyon, FRANCE
- index: 2
  name: Hospices Civils de Lyon, Lyon, FRANCE

---

# Summary

Data analysis is a crucial step in every research project in life science. Every clinician or researcher is one day faced with the need to perform statistical analysis. However, few free accessible solutions exist to date and most of the relevant software programs need a paid license. R is a free language which allows one to perform statistical analysis [@RCoreTeam:2017].
While the R environment is very powerful, its learning curve can be very steep in the beginning, especially for people with no previous coding skills or those with less time to learn them. A graphical user interface has already been provided as an independent package, but its features are limited for medical and applied life science studies and its usage remains difficult and unintuitive for new users [@Fox:2005]. Other free software programs exist, such as iNZight or Jamovi. However, while providing solutions with multiple features such as variable recoding, they do not guide the user through the analysis and can lack some key features such as time-dependent outcome analysis.

`StatAid` is a free open-source software provided as an R package which allows clinicians and researchers to perform statistical analysis through an intuitive graphical interface. It has been developed using the Shiny package [@Chang:2020], while Golem was used for package compilation and deployment [@Guyader:2020].

This software guides the user through all the steps of data analysis and includes multiple features such as:

- Exploratory data analysis: distribution, count, missing-values and outliers check

- Descriptive analysis, simple comparative analysis and publication ready 'table 1' output 

- Publication-ready graph customization 

- Paired data analysis (e.g. for repeated measures and matched case-control studies)

- Univariate analysis and models for continuous and categorical outcome: correlation, linear and logistic regression 

- Univariate analysis and models for time-dependent outcome: Kaplan-Meier curves and Cox regression 

- Multivariate analysis and models for continuous, categorical and time-dependent outcomes

Its user-friendly interface can guide clinicians or researchers, even those without previous software experience, through statistical analysis. In addition to a local version, a ready-to-use online version [(https://vincentalcazer.shinyapps.io/StatAid/)](https://vincentalcazer.shinyapps.io/StatAid/) is also available, providing access to `StatAid` everywhere, even on hospital/research center computers where external software installation is not allowed.
 

# Statement of need 

`StatAid` is an R package designed to fit the needs of every-day statistical analysis in science. This package provides all the tools necessary to perform data analysis in an intuitive, ready-to-use graphical interface. Users without coding skills or the availability of softwares with paid-licenses can easily perform all the steps of a good statistical analysis, from data-exploration/quality check to multivariate modeling. 

Compared with other similar free software, `StatAid` has been designed to quickly produce publication-ready graphs and tables by guiding the user through their data analysis and providing multiple graph customization options. By limiting the number of choices and integrating different checks and variable controls, `StatAid` helps the user prevent bad test use or bad graph choice. Besides, as an evolving software, `StatAid` also provides the possibility for users to request the implementation of additional features or to contribute to software development. Its open-source aspect can also be seen as a security for people working with sensitive data (e.g. data from clinical trials/patients). The online version of `StatAid` enables it to be accessable everywhere, even on computers with restrictive policies for software installation.

`StatAid` was designed to be used by clinicians, researchers, students, and any person wanting to perform statistical analysis with no prior coding skills. Primarily designed for medical/life science data analysis, `StatAid` can also easily be extended to other fields.



# References

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to StatAid

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

`StatAid` is a free open-source software provided as an R package
allowing clinicians and researchers to perform statistical analysis
through an intuitive graphical interface. It has been developed with the
R software, using the [Shiny package](https://shiny.rstudio.com/).
[Golem](https://github.com/ThinkR-open/golem) has been used for package
compilation and deployment.

The software guides the users through the steps of a good data analysis,
including multiple features such as:

<ul>

<li>

Exploratory data analysis: distribution, count, missing-values and
outliers check

</li>

<li>

Descriptive analysis, simple comparative analysis and publication ready
‘table 1’ output

</li>

<li>

Publication-ready graph customization

</li>

<li>

Paired data analysis (matched case-control studies, repeated measures)

</li>

<li>

Univariate analysis and models for continuous and categorical outcome:
Correlation, linear and logistic regression

</li>

<li>

Univariate analysis and models for time-dependent outcome: Kaplan-Meier
curves and cox regression

</li>

<li>

Multivariate analysis and models for continuous, categorical and
time-dependent outcomes

</li>

</ul>

# Getting started

## Online version

StatAid has a ready-to-use online version available at
<https://vincentalcazer.shinyapps.io/StatAid/>.

## Local version

You can install the development version from
[GitHub](https://github.com/VincentAlcazer/StatAid) else by cloning the
repository or directly by downloading the package in R:

``` r

install.packages("remotes")
remotes::install_github("VincentAlcazer/StatAid")

StatAid::run_app()
```

## Quick-start user guide

If you are not familiar with StatAid or just want to have an overview of
the different possibilities, you can check the [StatAid’s quick-start
user
guide](https://github.com/VincentAlcazer/StatAid/blob/master/STATAID_QUICK_START_USER_GUIDE.pdf)

## Citing StatAid

If you found StatAid useful and used it for your research, please cite
the [paper published in the Journal of Open Source
Software.](https://joss.theoj.org/papers/10.21105/joss.02630)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02630/status.svg)](https://doi.org/10.21105/joss.02630)

# Troubleshooting and contribution

All troubleshooting and contributions can be found on the [Github
page.](https://github.com/VincentAlcazer/StatAid/issues)

## Bug report

If you encounter any problem with the software or find a bug, please
report it on GitHub:

  - Create a [new
    issue](https://github.com/VincentAlcazer/StatAid/issues) on the
    Github page
  - Try to describe the problem/bug with reproductible steps

## Feature request

To ask for new feature implementation/current feature enhancemenet:

  - Create a [new
    issue](https://github.com/VincentAlcazer/StatAid/issues) on the
    Github page
  - Briefly describe the research question you want to answer and the
    type of data you have
  - If possible: provide pictures of the graph you would like to make or
    references from the paper you saw the analysis in.

## Contribution proposal

Contributions to new features or code enhancement are welcomed by
creating a new [pull
request.](https://github.com/VincentAlcazer/StatAid/pulls)
# StatAid 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
