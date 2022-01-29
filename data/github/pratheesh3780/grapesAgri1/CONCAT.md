# `grapesAgri1`: Collection of shiny applications for data analysis in Agriculture-Part 1 <img src="man/figures/logo.PNG" align="right" alt="logo" width="173" height = "200" style = "border: none; float: right;">
### General R-shiny based Analysis Platform Empowered by Statistics in Agriculture part-1 (grapesAgri1)

## R-Package for Data Analysis in Agriculture.
###### Version : 0.1.0; Copyright (C) 2021-2025: [Kerala Agricultural University](https://www.kaugrapes.com); License: [GPL-3](https://www.r-project.org/Licenses/) 

##### *Gopinath, P. P.<sup>1</sup>, Parsad, R.<sup>2</sup>, Joseph, B.<sup>1</sup>, Adarsh, V.S.<sup>3</sup>*

1.  Department of Agricultural Statistics, College of Agriculture, Vellayani, Kerala Agricultural Univesity.
2.  ICAR-Indian Agricultural Statistics Research Institute,
    New Delhi.
3.  Department of Agricultural Statistics, BCKV, West Bengal

---

![CRAN/METACRAN](https://img.shields.io/cran/v/grapesAgri1?style=for-the-badge)
![GitHub](https://img.shields.io/github/license/pratheesh3780/grapesAgri1)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version-last-release/grapesAgri1)](https://cran.r-project.org/package=grapesAgri1)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4923220.svg)](https://doi.org/10.5281/zenodo.4923220)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/pratheesh3780/grapesAgri1)
![GitHub R package version](https://img.shields.io/github/r-package/v/pratheesh3780/grapesAgri1)
![GitHub language count](https://img.shields.io/github/languages/count/pratheesh3780/grapesAgri1)
![Libraries.io dependency status for GitHub repo](https://img.shields.io/librariesio/github/pratheesh3780/grapesAgri1)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/pratheesh3780/grapesAgri1)
[![R-CMD-check](https://github.com/pratheesh3780/grapesAgri1/workflows/R-CMD-check/badge.svg)](https://github.com/pratheesh3780/grapesAgri1/actions)
[![](https://cranlogs.r-pkg.org/badges/grapesAgri1)](https://cran.r-project.org/package=grapesAgri1)
[![Build Status](https://www.travis-ci.com/pratheesh3780/grapesAgri1.svg?branch=master)](https://www.travis-ci.com/pratheesh3780/grapesAgri1)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03437/status.svg)](https://doi.org/10.21105/joss.03437)

---

## Introduction
<div align="justify">Agricultural experiments demands a wide range of statistical tools for analysis, which includes from Exploratory analysis, Design of experiments and Statistical genetics. Majority of the Agricultural scientists prefer graphical user interface for performing analysis . As R-shiny provides a platform to create interactive user interface, we have utilized it to produce interactive interfaces for commonly used analysis tools in Agrciultural experiments. grapesAgri1(General R-shiny based Analysis Platform Empowered by Statistics for data analysis in Agriculture-part1) is a collection of shiny based applications for some basic statistical analysis commonly used in agricultural research. It can be utilised by scientific community who prefers an interactive user interface. On using the functions in this package a Graphical User Interface will pop up. Apps Works by simple upload of files in CSV format. Results can be downloaded in HTML format. Plots and Graphs can be generated, which is also downloadable as .png file.</div>

## Installation
The package can be installed from CRAN as follows:

``` r
# Install from CRAN
install.packages('grapesAgri1', dependencies=TRUE)
```

The development version can be installed from github as follows:

``` r
# Install grapesAgri1 development version from Github using the code below:
if (!require('devtools')) install.packages('devtools')
devtools::install_github("pratheesh3780/grapesAgri1")
```

## usage
```r
grapesAgri1::descApp() # descriptive Statistics and Visualization 
grapesAgri1::corrApp() # Correlation Analysis
grapesAgri1::ttApp() # Compare Means
grapesAgri1::crdApp() # Completely Randomized Design
grapesAgri1::layoutApp() # Field layout of experiments
grapesAgri1::rbdApp() # Randomized Block Design 
```
## Apps included in the package

|Sl. No.| App Title | Function to call |Remark |
|:-----:| :----------- | :-----------:|:----------------|
|1|Descriptive Statistics and Visualization   | descApp()      |Summary Statistics, Summary Statistics by Group, Box plot, Histogram, Q-Q plot and Shapiro-Wilk's test|
|2|Correlation Analysis   | corrApp()      | Simple correlation, Correlation Matrix, correlogram and scatterplot|
|3|Compare Means: Small samle tests  | ttApp()      | One sample t-test, Two sample unpaired t-test, paired t-test, Two sample unpaired Welch t-test, F test, Box plot, Paired Plot|
|4|Completely Randomized Design  | crdApp()      |One-way Analysis of variance (equal and unequal replication), Multiple comparison test, boxplot and barchart with confidence interval|
|5|Field layout of experiments | layoutApp()      |Field layout of following designs can be obtained: Completely Randomized Design (CRD), Randomized Complete Block Design (RCBD), Split-plot design, Strip-plot design, Augmented RCBD|
|6|Randomized Block Design  | rbdApp()      |Two-way Analysis of variance, Multiple comparison test, boxplot and barchart with confidence interval|

## Further Reading
To know more about analysis tools included in the package see the following links
1. [Design and Analysis of experiments](http://apps.iasri.res.in/ebook/EBADAT/2-Basic%20Statistical%20Techniques/9-Fundamentals%20Of%20Designsf.pdf).
2. [Correlation Analysis](http://apps.iasri.res.in/ebook/EBADAT/2-Basic%20Statistical%20Techniques/6-Correlation_and_regression.pdf).
3. [Hypothesis Testing](http://apps.iasri.res.in/ebook/EBADAT/2-Basic%20Statistical%20Techniques/4-TEST%20OF%20HYPOTHESIS.pdf).
4. [Descriptive statistics and Exploratory data analysis](http://apps.iasri.res.in/ebook/EBADAT/2-Basic%20Statistical%20Techniques/1-Descriptive%20Statistics.pdf).
5. [Test for significance](http://apps.iasri.res.in/ebook/EBADAT/2-Basic%20Statistical%20Techniques/5-Tests%20of%20Significance-Seema.pdf).

## Glimpse to grapesAgri1 in Action!
It is very user friendly. Just upload your file in CSV format.

Note: we apologize that in grapesAgri1 version 1.0.0 in CRAN you may not be able to download model data set in crdApp()and rbdApp(). Issue will be cleared in the version 1.1.0 releasing by next month. You can instead download from github where the issue is resolved. 

Thank You

See below for some random images of GUI of grapesAgri1

![](man/figures/Corr.png) 

![](man/figures/Corr1.png)

![](man/figures/crd.png)  

![](man/figures/crd2.png)

![](man/figures/crd3.PNG)

![](man/figures/crd4.png)

![](man/figures/desc.png)

![](man/figures/layout.png)

![](man/figures/rbd.png)

![](man/figures/rbd1.png)

![](man/figures/rbd2.png)

# Community guidelines

Report Issues:

-   Questions, feedback, bug reports: please open an issue in the [issue tracker of the project](https://github.com/pratheesh3780/grapesAgri1/issues).

Contribution to the software:

-   Please open an issue in the issue tracker of the project that describes the changes you would like to make to the software and open a pull request with the changes. The description of the pull request must reference the corresponding issue.

# grapesAgri1 1.1.0

* Grey scale plots added in RCBD and CRD analysis  
* Bar charts with letter grouping added in RCBD and CRD  
* Download Files issue resolved  

# Contributing to grapesAgri1

This outlines how to propose a change to grapesAgri1. 
For more detailed info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib). 

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file. 
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("pratheesh3780/grapesAgri1", fork = TRUE)`.

*   Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.

### Code style

*   New code should follow the tidyverse [style guide](https://style.tidyverse.org). 
    You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles, but please don't restyle code that has nothing to do with your PR.  

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  

## Code of Conduct

Please note that the grapesAgri1 project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
---
title: 'grapesAgri1: Collection of Shiny Apps for Data Analysis in Agriculture'
tags:
  - R
  - Data analysis in Agriculture
  - shiny app
  - Design of experiments
  - Compare means
  - Field Layout
  - Correlation Analysis
  - Descriptive Statistics
  
authors:
  - name: Pratheesh P. Gopinath
    orcid: 0000-0003-3290-0436
    affiliation: "1"
  - name: Rajender Parsad
    affiliation: "2"
  - name: Brigit Joseph
    affiliation: "1"
  - name: Adarsh V. S.
    affiliation: "3"
    
affiliations:
 - name: Department of Agricultural Statistics, College of Agriculture, Vellayani, Kerala Agricultural University, Kerala, India.
   index: 1
 - name: ICAR-Indian Agricultural Statistics Research Institute, New Delhi, India.
   index: 2
 - name: Department of Agricultural Statistics, BCKV, West Bengal, India
   index: 3
   
date: June 11, 2021 
bibliography: paper.bib

---

# Summary

Agricultural experiments demand a wide range of statistical tools for analysis, which includes exploratory analysis, design of experiments, and statistical genetics. It is a challenge for scientists and students to find a suitable platform for data analysis and publish the research outputs in quality journals. Most of the software available for data analysis are proprietary or lack a simple user interface, for example SAS&reg; is available in ICAR (Indian Council of Agricultural Research) for data analysis, though it is a highly advanced statistical analysis platform, and its complexity holds back students and researchers from using it. Some web applications like WASP (https://ccari.res.in/waspnew.html) and OPSTAT (http://14.139.232.166/opstat/) used by the agricultural research community are user friendly but these applications don't provide options to generate plots and graphs.

The open source programming language R and associated ecosystem of packages, provides an excellent platform for data analysis but as of yet, is not heavily utilised by researchers in agricultural disciplines. Insufficient programming and computational knowledge are the primary challenges for agricultural researchers using R for analysis, as well as a preferences for researchers in agriculture to prefer a graphical user interface.

Efforts should therefore be made to develop a high quality, reliable open-source platform with a simple interactive user interface for data analysis in agriculture. Utilising the features of `shiny` package in R, we have developed a collection of `shiny` apps for agricultural research called `grapesAgri1` (General R `shiny` based Analysis Platform Empowered by Statistics for data analysis in Agriculture-part1). `grapesAgri1` is hosted on the web (http://www.kaugrapes.com), downloadable as a standalone application, and now released as an R package.

By calling the functions in the `grapesAgri1` package, a graphical user interface will open. Apps are self-explanatory and work by uploading files in CSV format. Results can be downloaded in a HTML format. Plots and graphs can be generated, which are also downloadable as .png files.

# Statement of need

India has one of the world's largest agricultural education systems. According to the Indian Council of Agricultural Research (ICAR), the main regulator of agricultural education in India, there are 63 State Agricultural Universities, 4 deemed universities, and 3 Central Agricultural Universities. These institutions enroll annually about 15,000 students in undergraduate programs and over 7,000 students in post graduate programs and more than 2000 at PhD level. At any point, there are over 75,000 students studying in these institutions. Research activities are performed actively in all these universities. `grapesAgri1` not only serves as a platform for data analysis but also can be used as a teaching tool in agricultural statistics. `grapesAgri1` includes some basic statistical tools which were covered in the syllabi of undergraduate programs as well as in post graduate programs.

# Information for Users

`grapesAgri1` is a collection of six `shiny` applications. Details of six applications are given below. Model dataset for testing can be downloaded from the main window of the application. Details for preparing CSV files are also included in the main window. Users just need to click on the browse button and upload the CSV file for analysis.

## Apps included in the package

|Sl. No.| App Title | Function to call |Utility|
|:-----:| :----------- | :-----------:|:----------------|
|1|Descriptive Statistics and Visualization   | descApp()      |Summary statistics, summary statistics by group, box plot, histogram, Q-Q plot and Shapiro-Wilk's test|
|2|Correlation Analysis   | corrApp()      | Simple correlation, correlation Matrix, correlogram and scatterplot|
|3|Compare Means: Small samle tests  | ttApp()      | One sample t-test, two sample unpaired t-test, paired t-test, two sample unpaired Welch t-test, F test, box plot, paired Plot|
|4|Completely Randomized Design  | crdApp()      |One-way ANOVA (equal or unequal replications), multiple comparison test, boxplot and barchart with confidence interval|
|5|Field layout of experiments | layoutApp()      |Field layouts of following designs can be obtained: completely randomized design (CRD), randomized complete block design (RCBD), split-plot design, strip-plot design, augmented RCBD|
|6|Randomized Block Design  | rbdApp()      |Two-way ANOVA, multiple comparison test, boxplot and barchart with confidence interval|

The package can be installed from CRAN as follows:

``` r
# Install from CRAN
install.packages('grapesAgri1', dependencies=TRUE)
```

The development version can be installed from GitHub as follows:

``` r
# Install grapesAgri1 development version from Github using the code below:
if (!require('devtools')) install.packages('devtools')
devtools::install_github("pratheesh3780/grapesAgri1")
```

## Package dependencies and details of functions used 

DescApp() function uses `descr` and `stby` functions of `summarytools` package [@Dominic_Comtois_2021] to calculate summary statistics and summary statistics by group. `knitr` [@Yihui_Xie_2021] and `kableExtra`[@Hao_zhu_2021] packages were used to produce HTML tables. `shapiro.test`, `qqnorm` and `qqline` functions of `stats` package were used for the Test of Homogeneity of variance and obtaining Q-Q plot. `hist` and `boxplot` of package `graphics` were used to obtain histogram and boxplot respectively. `ggqqplot` of package `ggpubr `[@ggpubr_2020] is also used to plot Q-Q plot in the app.

CorrApp() function uses `cor.test` to calculate correlation. Correlation matrix is calculated using `rcorr` function in `Hmisc` package [@hmisc_2021]. Correlogram is obtained using `corrplot` function in `corrplot`[@corrplot2021] package.

ttApp()function uses `t.test` function to calculate t statistic. Descriptive statistics were calculated using `stat.desc` function of `pastecs` package. `var.test` function is used for F-test. `ggboxplot` function of `ggpubr` [@ggpubr_2020] package is used to draw boxplot. Paired plot is obtained using `paired` function of package `PairedData`[@paired_data2018].

crdApp() uses `anova` function of `stats` package to obtain one-way ANOVA. `LSD.test`,`duncan.test` and `HSD.test` functions of `agricolae` [@agricolae_2020] package is used for multiple comparison test like LSD,DMRT and Tukey respectively. `ggboxplot` function of `ggpubr` [@ggpubr_2020] package is used for boxplot. `ggplot` function of `ggplot2`[@ggplot_2016] is used for barchart with confidence interval.

layoutApp() uses `design.crd`, `design.rcbd`, `design.dau`, `design.strip`, `design.split` functions of package `agricolae` [@agricolae_2020] to generate random layout of designs. Field layout were plotted using `desplot` function in `desplot` package [@desplot_2020].

rbdApp() uses `anova` function of `stats` package to obtain two-way ANOVA. `LSD.test`,`duncan.test` and `HSD.test` functions of `agricolae` package [@agricolae_2020] is used for multiple comparison test like LSD,DMRT and Tukey respectively. `ggboxplot` function of `ggpubr` package [@ggpubr_2020] is used for boxplot. `ggplot` function of `ggplot2` [@ggplot_2016] is used for barchart with confidence interval.

# Usage

``` r
grapesAgri1::descApp() # descriptive Statistics and Visualization 
grapesAgri1::corrApp() # Correlation Analysis
grapesAgri1::ttApp() # Compare Means
grapesAgri1::crdApp() # Completely Randomized Design
grapesAgri1::layoutApp() # Field layout of experiments
grapesAgri1::rbdApp() # Randomized Block Design 
```
# Acknowledgements

We wish to thank Kerala Agricultural University and Regional Agricultural Reaserch Station, Vellayani, Thiruvananthapuram for the financial support through the revolving fund scheme and observation trial.

# References
