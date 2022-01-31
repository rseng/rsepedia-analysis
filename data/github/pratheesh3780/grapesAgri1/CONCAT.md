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
---
title: "Design Layout"
output: html_document
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
library(rmarkdown)
library(agricolae)
library(dplyr)
library(ggplot2)
library(magrittr)
library(desplot)
library(knitr)
```

Details of the Experiment is given below.

```{r , echo=FALSE}


  if(input$design == 'crd'){
    if(input$rep_crd>0 && input$submit1>0){
      Treatments= input$trt_crd
      Replication=input$rep_crd
      Experimental.units=Treatments*Replication
      det=cbind(Treatments,Replication,Experimental.units)
      det=as.data.frame(det)
      kable(det,row.names = FALSE) %>%
      kable_styling() %>%
      kable_paper("hover", full_width = F)
    }
  }

    if(input$design == 'rbd'){
      if(input$rep_rbd>0 && input$submit2>0){
        Treatments= input$trt_rbd
        Blocks=input$rep_rbd
        Experimental.units=Treatments*Blocks
        det=cbind(Treatments,Blocks,Experimental.units)
        det=as.data.frame(det)
        kable(det,row.names = FALSE) %>%
      kable_styling() %>%
      kable_paper("hover", full_width = F)
      }
    }
    
    if(input$design == 'aug'){
      if(input$check_aug>0&&input$trt_aug&&input$rep_aug>0 && input$submit3>0){
        Check =input$check_aug
        Treatments= input$trt_aug
        Blocks=input$rep_aug
        Plots =Check*Blocks+Treatments
        det=cbind(Check,Treatments,Blocks,Plots)
        det=as.data.frame(det)
        kable(det,row.names = FALSE) %>%
      kable_styling() %>%
      kable_paper("hover", full_width = F)
      }
    }
  
    if(input$design == 'split'){
      if(input$maintrt_split>0&&input$subtrt_split&&input$rep_split>0 &&
         input$submit4>0){
        Main_plot_trt =input$maintrt_split
        Sub_plot_trt =input$subtrt_split
        Replications =input$rep_split
        No_of_Plots =Main_plot_trt*Sub_plot_trt*Replications
        det=cbind(Main_plot_trt,Sub_plot_trt,Replications,No_of_Plots)
        det=as.data.frame(det)
        kable(det,row.names = FALSE) %>%
      kable_styling() %>%
      kable_paper("hover", full_width = F)
      }
    }

    if(input$design == 'strip'){
      if(input$maintrt_strip>0&&input$subtrt_strip&&input$rep_strip>0 &&
         input$submit5>0){
        Main_plot_A =input$maintrt_strip
        Main_plot_B =input$subtrt_strip
        Replications =input$rep_strip
        No_of_Plots =Main_plot_A*Main_plot_B*Replications
        det=cbind(Main_plot_A,Main_plot_B,Replications,No_of_Plots)
        det=as.data.frame(det)
        kable(det,row.names = FALSE) %>%
      kable_styling() %>%
      kable_paper("hover", full_width = F)
      }
    }

```


Field layout of your experiment

```{r echo=FALSE, error=FALSE, message=FALSE, warning=FALSE, include=TRUE,warn.conflicts=FALSE}

 if (input$design == "crd") {
  if (input$trt_crd > 0) {
    if (input$rep_crd > 0) {
      if (input$submit1 > 0) {
        t = input$trt_crd
        r = input$rep_crd
        s = input$size
        rand = input$rand
        trtname = sprintf("T%d", 1:t)
        outdesign = agricolae::design.crd(trtname, r = r, seed = rand, serie = 0, kinds = "Super-Duper", randomization = TRUE)
        CRD = outdesign$book
        CRD = CRD[order(CRD$r), ]
        CRD$col = CRD$r
        CRD$row = rep(1:t, r)
        desplot::desplot(
          form = trtname ~ col + row, data = CRD, text = trtname, out1 = col,
          out2 = row, out2.gpar = list(col = "black", lwd = 3),
          cex = s, main = "Layout of Completely Randomized Design", show.key = FALSE
        )
      }
    }
  }
}
if (input$design == "rbd") {
  if (input$trt_rbd > 0) {
    if (input$rep_rbd > 0) {
      if (input$submit2 > 0) {
        t = input$trt_rbd
        rep = input$rep_rbd
        s = input$size_rbd
        trtname = sprintf("T%d", 1:t)
        n = t * rep
        outdesign = agricolae::design.rcbd(trt = trtname, r = rep, seed = input$rand2, serie = 0, kinds = "Super-Duper", randomization = TRUE)
        RBD = outdesign$book
        RBD = RBD[order(RBD$block), ]
        RBD$blocks = as.numeric(RBD$block)
        RBD$plots = rep(1:t, rep)
        desplot::desplot(
          form = block ~ blocks + plots, data = RBD, text = trtname, out1 = blocks,
          out2 = plots, out2.gpar = list(col = "green"),
          cex = s, main = "Layout of Randomized Block Design", show.key = TRUE
        )
      }
    }
  }
}
  
  if(input$design == 'aug'){
    if(input$check_aug> 0){
      if(input$trt_aug> 0){
        if(input$rep_aug> 0){
            if(input$submit3>0){
              NC=input$check_aug # no. of check
              blk=input$rep_aug # no: of Blocks
              trt=input$trt_aug # no.of trt
              s=input$size_aug
              T1=sprintf("C%d", 1:NC) # checks
              T2=sprintf("T%d", 1:trt) # treatments
              outdesign=agricolae::design.dau(T1,T2, r=blk,seed=input$rand3,serie=0,randomization = TRUE)
              aug= outdesign$book
              aug = aug[order(aug$block),]
              aug$blocks = as.numeric(aug$block)
              x=aug %>%
                group_by(block) %>%
                mutate(row = row_number())
              aug1=as.data.frame(cbind(aug,row=x$row))
              desplot::desplot(form=block~ blocks+row, data=aug1, text=trt, out1=blocks, out2=row,
                               cex=s, main="Layout of Augmented Randomized Block Design",show.key=FALSE)
              
            }
        }
      }
    }
  }
  
 if(input$design == 'strip'){
    if(input$maintrt_strip> 0){
      if(input$subtrt_strip> 0){
        if(input$rep_strip> 0){
            if(input$submit5>0){
              a=input$maintrt_strip
              b=input$subtrt_strip
              rep=input$rep_strip
              s=input$size_strip
              main1=sprintf("A%d", 1:a)
              main2=sprintf("B%d", 1:b)
              outdesign =agricolae::design.strip(main1,main2,r=rep,serie=0,seed=input$rand5,kinds="Super-Duper",randomization=TRUE)
              strip= outdesign$book
              strip$block2=strip$block
              v1 = do.call(paste, as.data.frame(t(apply(strip[4:5], 1, sort))))
              strip$row =  match(v1, unique(v1))
              v2 = do.call(paste, as.data.frame(t(apply(strip[2:3], 1, sort))))
                 strip$column =  match(v2, unique(v2))
              desplot::desplot(form=main1~ row+column, data=strip,text=main2, out1=column,
                               out2=row,out2.gpar = list(col = "#a83232"),
                               cex=s, main="Layout of Strip-Plot Design",show.key=TRUE)
              
            }
        }
      }
    }
  }
  
 if(input$design == 'split'){
    if(input$maintrt_split> 0){
      if(input$subtrt_split> 0){
        if(input$rep_split> 0){
            if(input$submit4>0){
              a=input$maintrt_split
              b=input$subtrt_split
              rep=input$rep_split
              s=input$size_split
              main=sprintf("A%d", 1:a)
              sub=sprintf("b%d", 1:b)
              outdesign =agricolae::design.split(main,sub,r=rep,serie=0,seed=input$rand4,kinds="Super-Duper",randomization=TRUE)
              split = outdesign$book
              split = split[order(split$plots),]
              split$plots = split$plots
              split$plots=as.numeric(split$plots)
              split$splots=as.numeric(split$splots)
              desplot::desplot(form=main~ plots+splots, data=split,col=main, text=sub, out1=plots, out2=splots,
                               cex=s, main="Layout of Split-Plot Design",show.key=TRUE)
              
            }
        }
      }
    }
  }
```


```{r , echo=FALSE}
if (input$design == "crd") {
  if (input$table_butt1 > 0) {
    t = input$trt_crd
    r = input$rep_crd
    n = t * r
    trtname = sprintf("T%d", 1:t)
    outdesign = agricolae::design.crd(trtname, r = r, seed = input$rand, serie = 0, kinds = "Super-Duper", randomization = TRUE)
    CRD = outdesign$book
    CRD = CRD[order(CRD$r), ]
    CRD$exp = rep(1:n)
    final = as.data.frame(cbind(Experimental.Unit = CRD$exp, Treatments = as.character(CRD$trtname)))
    kable(final, caption = "Layout of CRD", row.names = FALSE)
  }
}

if (input$design == "rbd") {
  if (input$table_butt2 > 0) {
    t = input$trt_rbd
    rep = input$rep_rbd
    trtname = sprintf("T%d", 1:t)
    n = t * rep
    outdesign = agricolae::design.rcbd(trt = trtname, r = rep, seed = input$rand2, serie = 0, kinds = "Super-Duper", randomization = TRUE)
    RBD = outdesign$book
    RBD = RBD[order(RBD$block), ]
    RBD$plot = rep(1:t, rep)
    final = as.data.frame(cbind(Block = RBD$block, Plot = RBD$plot, Treatment = as.character(RBD$trtname)))
    kable(final, caption = "Layout of RBD", row.names = FALSE)
  }
}


if (input$design == "aug") {
  if (input$table_butt3 > 0) {
    NC = input$check_aug # no. of check
    blk = input$rep_aug # no: of Blocks
    trt = input$trt_aug # no.of trt
    s = input$size_aug
    T1 = sprintf("C%d", 1:NC) # checks
    T2 = sprintf("T%d", 1:trt) # treatments
    outdesign = agricolae::design.dau(T1, T2, r = blk, seed = input$rand3, serie = 0, randomization = TRUE)
    aug = outdesign$book
    aug = aug[order(aug$block), ]
    x1 = aug %>%
      group_by(block) %>%
      mutate(Plot = row_number())
    final = as.data.frame(cbind(Block = x1$block, Plot = x1$Plot, Treatment = as.character(x1$trt)))
    kable(final, caption = "Layout of Augmented Design", row.names = FALSE)
  }
}


if (input$design == "split") {
  if (input$table_butt4 > 0) {
    a = input$maintrt_split
    b = input$subtrt_split
    rep = input$rep_split
    s = input$size_split
    main = sprintf("A%d", 1:a)
    sub = sprintf("b%d", 1:b)
    outdesign = agricolae::design.split(main, sub, r = rep, serie = 0, seed = input$rand4, kinds = "Super-Duper", randomization = TRUE)
    split = outdesign$book
    split = split[order(split$plots), ]
    final = as.data.frame(cbind(
      Main_plot = split$plots, Sub_plot = split$splots,
      Replication = split$block, Main_treatment = as.character(split$main), Sub_treatment = as.character(split$sub)
    ))
    kable(final, caption = "Layout of Split-Plot design", row.names = FALSE)
  }
}

if (input$design == "strip") {
  if (input$table_butt5 > 0) {
    a = input$maintrt_strip
    b = input$subtrt_strip
    rep = input$rep_strip
    s = input$size_strip
    main1 = sprintf("A%d", 1:a)
    main2 = sprintf("B%d", 1:b)
    outdesign = agricolae::design.strip(main1, main2, r = rep, serie = 0, seed = input$rand5, kinds = "Super-Duper", randomization = TRUE)
    strip = outdesign$book
    strip$block2 = strip$block
    v1 = do.call(paste, as.data.frame(t(apply(strip[4:5], 1, sort))))
    strip$row = match(v1, unique(v1))
    v2 = do.call(paste, as.data.frame(t(apply(strip[2:3], 1, sort))))
    strip$column = match(v2, unique(v2))
    final = as.data.frame(cbind(
      Replication = strip$block, Horizontal_row_No. = strip$row,
      Treatment_A = as.character(strip$main1),
      Vertical_row_No. = strip$column,
      Treatment_B = as.character(strip$main2)
    ))
    kable(final, caption = "Layout of Strip plot design", row.names = FALSE)
  }
}
  
```
**Package: grapesAgri1, Version;1.0.0**
  
  
  
---
title: "Two way ANOVA (RBD)"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(shiny)
library(ggplot2)
library(Hmisc)
library(agricolae)
library(shinyWidgets)
library(shinycssloaders)
library(dplyr)
library(gtools)
library(knitr)
library(kableExtra)
library(magrittr)
```
## RESULT
RBD ANALYSIS

```{r, echo = FALSE}
csvfile <- reactive({
  csvfile <- input$file1
  if (is.null(csvfile)) {
    return(NULL)
  }
  dt <- read.csv(csvfile$datapath, header = input$header, sep = ",")
  dt
})

if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  r <- as.numeric(input$rep)
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  replication <- d[, input$Replication]
  treatment <- factor(treatment)
  replication <- factor(replication)
  anvaTable <- lm(response ~ treatment + replication)
  result <- as.data.frame(anova(anvaTable))
  out <- agricolae::LSD.test(csvfile()[, input$yield], csvfile()[, input$treatment], result[3, 1], result[3, 3])
  trtmeans <- out$means
   trtmeans <- trtmeans[ gtools::mixedsort(row.names(trtmeans)), ]
  colnames(trtmeans)[1] <- "Treatment_means"
  drops <- c("r", "Q25", "Q50", "Q75")
  trtmeans <- trtmeans[, !(names(trtmeans) %in% drops)]
  kable(trtmeans, caption = "Treatment mean and other statistics", digits = 3, align = "c", row.names = TRUE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
}
tags$br()
################ anova table

if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  r <- as.numeric(input$rep)
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  replication <- d[, input$Replication]
  treatment <- factor(treatment)
  replication <- factor(replication)
  anvaTable <- lm(response ~ treatment + replication)
  result <- as.data.frame(anova(anvaTable))
  SoV <- c("Treatment", "Replication", "Error")
  final <- cbind(SoV, result)
  kable(final, caption = "ANOVA TABLE", digits = 3, align = "c", row.names = FALSE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
}
tags$br()
############################################## SEM
if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  r <- as.numeric(input$rep)
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  replication <- d[, input$Replication]
  treatment <- factor(treatment)
  replication <- factor(replication)
  anvaTable <- lm(response ~ treatment + replication)
  result <- as.data.frame(anova(anvaTable))
  out <- LSD.test(d[, input$yield], d[, input$treatment], result[3, 1], result[3, 3])
  colnam <- c("MSE", "SE(d)", "SE(m)", "CV")
  stat <- out$statistics
  MSE <- stat[1, 1]
  SED <- sqrt((2 * MSE) / r)
  SEM <- sqrt(MSE / r)
  CV <- stat[1, 4]
  Result <- cbind(MSE, SED, SEM, CV)
  colnames(Result) <- colnam
  kable(Result, caption = "SEM & Other statistics", digits = 3, align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
}
tags$br()
if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  r <- as.numeric(input$rep)
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  replication <- d[, input$Replication]
  treatment <- factor(treatment)
  replication <- factor(replication)
  anvaTable <- lm(response ~ treatment + replication)
  result <- as.data.frame(anova(anvaTable))
  if (result[1, 5] <= 0.05) {
    cat("Since the P-value in ANOVA table is < 0.05,\nthere is a significant difference between atleast a pair of treatments,\nso multiple comparison is required to identify best treatment(s)")
  }
  if (result[1, 5] > 0.05) {
    cat("Treatment means are not significantly different")
  }
}
############################### CD and measures
tags$br()
if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  r <- as.numeric(input$rep)
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  replication <- d[, input$Replication]
  treatment <- factor(treatment)
  replication <- factor(replication)
  anvaTable <- lm(response ~ treatment + replication)
  result <- as.data.frame(anova(anvaTable))
  if (result[1, 5] > 0.05) {
    return()
  }

  else if (input$req == "lsd") {
    d <- as.data.frame(csvfile())
    r <- as.numeric(input$rep)
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    replication <- d[, input$Replication]
    treatment <- factor(treatment)
    replication <- factor(replication)
    anvaTable <- lm(response ~ treatment + replication)
    result <- as.data.frame(anova(anvaTable))
    out <- agricolae::LSD.test(d[, input$yield], d[, input$treatment], result[3, 1], result[3, 3])
    colnam <- c("MSE", "SE(d)", "SE(m)", "CD", "t value", "CV")
    stat <- out$statistics
    MSE <- stat[1, 1]
    SED <- sqrt((2 * MSE) / r)
    SEM <- sqrt(MSE / r)
    CD <- stat[1, 6]
    t <- stat[1, 5]
    CV <- stat[1, 4]
    Result <- cbind(MSE, SED, SEM, CD, t, CV)
    colnames(Result) <- colnam
    kable(Result, caption = "LSD test", digits = 3, align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }
  else if (input$req == "dmrt") {
    d <- as.data.frame(csvfile())
    r <- as.numeric(input$rep)
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    replication <- d[, input$Replication]
    treatment <- factor(treatment)
    replication <- factor(replication)
    anvaTable <- lm(response ~ treatment + replication)
    result <- as.data.frame(anova(anvaTable))
    out <- agricolae::duncan.test(d[, input$yield], d[, input$treatment], result[3, 1], result[3, 3], alpha = 0.05, group = TRUE, main = NULL, console = FALSE)
    Result <- out$duncan
    kable(Result, caption = "DMRT", digits = 3, align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }
  else if (input$req == "tukey") {
    d <- as.data.frame(csvfile())
    r <- as.numeric(input$rep)
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    replication <- d[, input$Replication]
    treatment <- factor(treatment)
    replication <- factor(replication)
    anvaTable <- lm(response ~ treatment + replication)
    result <- as.data.frame(anova(anvaTable))
    out <- agricolae::HSD.test(d[, input$yield], d[, input$treatment], result[3, 1], result[3, 3], alpha = 0.05, group = TRUE, main = NULL, unbalanced = FALSE, console = FALSE)
    Result <- out$statistics
    kable(Result, caption = "Tukey's HSD", digits = 3, align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }
}
tags$br()
########################################### grouping
if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  r <- as.numeric(input$rep)
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  replication <- d[, input$Replication]
  treatment <- factor(treatment)
  replication <- factor(replication)
  anvaTable <- lm(response ~ treatment + replication)
  result <- as.data.frame(anova(anvaTable))
  if (result[1, 5] > 0.05) {
    return()
  }

  else if (input$req == "lsd") {
    d <- as.data.frame(csvfile())
    r <- as.numeric(input$rep)
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    replication <- d[, input$Replication]
    treatment <- factor(treatment)
    replication <- factor(replication)
    anvaTable <- lm(response ~ treatment + replication)
    result <- as.data.frame(anova(anvaTable))
    out <- agricolae::LSD.test(d[, input$yield], d[, input$treatment], result[3, 1], result[3, 3])
    outgroup <- out$groups
    colnames(outgroup) <- c("trt_mean", "grouping")
    kable(outgroup, caption = "Treatment Grouping", align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }
  else if (input$req == "dmrt") {
    d <- as.data.frame(csvfile())
    r <- as.numeric(input$rep)
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    replication <- d[, input$Replication]
    treatment <- factor(treatment)
    replication <- factor(replication)
    anvaTable <- lm(response ~ treatment + replication)
    result <- as.data.frame(anova(anvaTable))
    out <- agricolae::duncan.test(d[, input$yield], d[, input$treatment], result[3, 1], result[3, 3], alpha = 0.05, group = TRUE, main = NULL, console = FALSE)
    outgroup <- out$groups
    colnames(outgroup) <- c("trt_mean", "grouping")
    kable(outgroup, caption = "Treatment Grouping", align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }
  else if (input$req == "tukey") {
    d <- as.data.frame(csvfile())
    r <- as.numeric(input$rep)
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    replication <- d[, input$Replication]
    treatment <- factor(treatment)
    replication <- factor(replication)
    anvaTable <- lm(response ~ treatment + replication)
    result <- as.data.frame(anova(anvaTable))
    out <- agricolae::HSD.test(d[, input$yield], d[, input$treatment], result[3, 1], result[3, 3], alpha = 0.05, group = TRUE, main = NULL, unbalanced = FALSE, console = FALSE)
    outgroup <- out$groups
    colnames(outgroup) <- c("trt_mean", "grouping")
    kable(outgroup, caption = "Treatment Grouping", align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }
}
############
if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  r <- as.numeric(input$rep)
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  replication <- d[, input$Replication]
  treatment <- factor(treatment)
  replication <- factor(replication)
  anvaTable <- lm(response ~ treatment + replication)
  result <- as.data.frame(anova(anvaTable))
  if (result[1, 5] > 0.05) {
    return("")
  }
  if (input$req == "tukey" || input$req == "dmrt" || input$req == "lsd") {
    cat("Treatments with same letters are not significantly different")
  }
}

tags$br()
tags$br()
```
Report generated from: **Package: grapesAgri1, Version 1.0.0**
---
title: "Completely Randomized Design"
output: html_document
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
library(knitr)
library(Hmisc)
library(agricolae)
library(dplyr)
library(magrittr)
library(kableExtra)
library(gtools)
```
## RESULT
CRD ANALYSIS

```{r, echo = FALSE}
csvfile <- reactive({
  csvfile <- input$file1
  if (is.null(csvfile)) {
    return(NULL)
  }
  dt <- read.csv(csvfile$datapath, header = input$header, sep = ",", check.names = FALSE)
  dt
})
################################ trt means
if (input$submit > 0) {
  input$reload
  Sys.sleep(2)
  d <- as.data.frame(csvfile())
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  treatment <- factor(treatment)
  anvaTable <- lm(response ~ treatment)
  result <- as.data.frame(anova(anvaTable))
  out <- LSD.test(csvfile()[, input$yield], csvfile()[, input$treatment], result[2, 1], result[2, 3])
  trtmeans <- out$means
  trtmeans <- trtmeans[ gtools::mixedsort(row.names(trtmeans)), ]
  colnames(trtmeans)[1] <- "Treatment_means"
  drops <- c("r", "Q25", "Q50", "Q75")
  result <- trtmeans[, !(names(trtmeans) %in% drops)]
  kable(result, digits = 3, caption = "Treatment mean and other statistics", align = "c", row.names = TRUE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
}
################ finish
tags$br()
################################ ANOVA TABLE
if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  treatment <- factor(treatment)
  anvaTable <- lm(response ~ treatment)
  result <- as.data.frame(anova(anvaTable))
  SoV <- c("Treatment", "Error")
  final <- cbind(SoV, result)
  kable(final, digits = 3, caption = "ANOVA TABLE", align = "c", row.names = FALSE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
}
############################# finish
tags$br()
###################### SEM
if (input$submit > 0) {
  if (input$filerepli == "equal") {
    d <- as.data.frame(csvfile())
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    treatment <- factor(treatment)
    anvaTable <- lm(response ~ treatment)
    result <- as.data.frame(anova(anvaTable))
    out <- LSD.test(csvfile()[, input$yield], csvfile()[, input$treatment], result[2, 1], result[2, 3])
    colnam <- c("MSE", "SE(d)", "SE(m)", "CV(%)")
    stat <- out$statistics
    repl <- out$means
    r <- repl[1, 3]
    MSE <- stat[1, 1]
    SED <- sqrt((2 * MSE) / r)
    SEM <- sqrt(MSE / r)
    CV <- stat[1, 4]
    Result <- cbind(MSE, SED, SEM, CV)
    colnames(Result) <- colnam
    kable(Result, digits = 3, caption = "Other important Statistics", align = "c", row.names = FALSE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }
}
#################################### finish
tags$br()
########################### TEXT INFERENCE
if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  treatment <- factor(treatment)
  anvaTable <- lm(response ~ treatment)
  result <- as.data.frame(anova(anvaTable))
  if (result[1, 5] <= 0.05) {
    cat("Since the P-value in ANOVA table is < 0.05,\nthere is a significant difference between atleast a pair of treatments,\nso multiple comparison is required to identify best treatment(s)")
  }
  if (result[1, 5] > 0.05) {
    cat("Treatment means are not significantly different")
  }
}
#########################
tags$br()
################################### Multi
if (input$submit > 0) {
  if (input$filerepli == "equal") {
    d <- as.data.frame(csvfile())
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    treatment <- factor(treatment)
    anvaTable <- lm(response ~ treatment)
    result <- as.data.frame(anova(anvaTable))
    if (result[1, 5] > 0.05) {
      return("Multiple comparison test is not performed")
    }

    else if (input$req == "lsd") {
      d <- as.data.frame(csvfile())
      t <- as.numeric(input$trt)
      response <- d[, input$yield]
      treatment <- d[, input$treatment]
      treatment <- factor(treatment)
      anvaTable <- lm(response ~ treatment)
      result <- as.data.frame(anova(anvaTable))
      out <- LSD.test(d[, input$yield], d[, input$treatment], result[2, 1], result[2, 3])
      result1 <- out$statistics
      kable(result1, caption = "LSD test", digits = 3, align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
    }


    else if (input$req == "dmrt") {
      d <- as.data.frame(csvfile())
      t <- as.numeric(input$trt)
      response <- d[, input$yield]
      treatment <- d[, input$treatment]
      treatment <- factor(treatment)
      anvaTable <- lm(response ~ treatment)
      result <- as.data.frame(anova(anvaTable))
      out <- duncan.test(d[, input$yield], d[, input$treatment], result[2, 1], result[2, 3], alpha = 0.05, group = TRUE, main = NULL, console = FALSE)
      result1 <- out$duncan
      kable(result1, caption = "DMRT", digits = 3, align = "c", row.names = FALSE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
    }


    else if (input$req == "tukey") {
      d <- as.data.frame(csvfile())
      t <- as.numeric(input$trt)
      response <- d[, input$yield]
      treatment <- d[, input$treatment]
      treatment <- factor(treatment)
      anvaTable <- lm(response ~ treatment)
      result <- as.data.frame(anova(anvaTable))
      out <- HSD.test(d[, input$yield], d[, input$treatment], result[2, 1], result[2, 3], alpha = 0.05, group = TRUE, main = NULL, unbalanced = FALSE, console = FALSE)
      result1 <- out$statistics
      kable(result1, caption = "Tukey", digits = 3, align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
    }
  }
}
########################## finish
tags$br()
##################### CD matrix
if (input$submit > 0) {
  if (input$filerepli == "unequal") {
    d <- as.data.frame(csvfile())
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    treatment <- factor(treatment)
    anvaTable <- lm(response ~ treatment)
    result <- as.data.frame(anova(anvaTable))
    if (result[1, 5] > 0.05) {
      return("Multiple comparison test is not performed")
    }
    else if (input$req == "lsd") {
      d <- as.data.frame(csvfile())
      t <- as.numeric(input$trt)
      response <- d[, input$yield]
      treatment <- d[, input$treatment]
      treatment <- factor(treatment)
      anvaTable <- lm(response ~ treatment)
      result <- as.data.frame(anova(anvaTable))
      count <- table(d[, input$treatment]) # count the number of replications of treatment
      t <- (result[1, 1] + 1) # no.of treatments
      repli <- as.data.frame(count)
      reciproc <- 1 / repli[, 2] # 1/ri
      npw <- choose(t, 2) # number of pairwise combinations
      sumres <- as.vector(apply(combn(reciproc, 2), 2, sum)) # all pairwise sum of reciprocals (1/ri+1/rj)
      ems <- as.vector(replicate(npw, (result[2, 3])))
      SE_D <- sqrt((ems * sumres)) # standard error of difference
      tvalue <- qt(0.975, result[2, 1]) # tvalue
      vect <- replicate(npw, tvalue) # vector of t value
      CD <- vect * SE_D # critical difference
      means <- aggregate(response, list(treatment), mean)
      std <- aggregate(response, list(treatment), sd)
      finalmean <- cbind(means, std[, 2])
      rownames(finalmean) <- NULL
      colnam <- c("Treatment", "mean", "std")
      colnames(finalmean) <- colnam
      b <- matrix(0, t, t)
      b[upper.tri(b, diag = FALSE)] <- CD
      name <- finalmean$Treatment
      colnames(b) <- name
      row.names(b) <- name
      kable(b, caption = "LSD test Matrix of CD values", digits = 3, align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
    }
  }
}
###################
tags$br()
################################### group
if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  treatment <- factor(treatment)
  anvaTable <- lm(response ~ treatment)
  result <- as.data.frame(anova(anvaTable))
  if (result[1, 5] > 0.05) {
    return("")
  }

  else if (input$req == "lsd") {
    d <- as.data.frame(csvfile())
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    treatment <- factor(treatment)
    anvaTable <- lm(response ~ treatment)
    result <- as.data.frame(anova(anvaTable))
    out <- LSD.test(d[, input$yield], d[, input$treatment], result[2, 1], result[2, 3])
    outgroup <- out$groups
    colnames(outgroup) <- c("trt_mean", "grouping")
    kable(outgroup, caption = "Treatment Grouping", align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }

  else if (input$req == "dmrt") {
    d <- as.data.frame(csvfile())
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    treatment <- factor(treatment)
    anvaTable <- lm(response ~ treatment)
    result <- as.data.frame(anova(anvaTable))
    out <- duncan.test(d[, input$yield], d[, input$treatment], result[2, 1], result[2, 3], alpha = 0.05, group = TRUE, main = NULL, console = FALSE)
    outgroup <- out$groups
    colnames(outgroup) <- c("trt_mean", "grouping")
    kable(outgroup, caption = "Treatment Grouping", align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }

  else if (input$req == "tukey") {
    d <- as.data.frame(csvfile())
    t <- as.numeric(input$trt)
    response <- d[, input$yield]
    treatment <- d[, input$treatment]
    treatment <- factor(treatment)
    anvaTable <- lm(response ~ treatment)
    result <- as.data.frame(anova(anvaTable))
    out <- HSD.test(d[, input$yield], d[, input$treatment], result[2, 1], result[2, 3], alpha = 0.05, group = TRUE, main = NULL, unbalanced = FALSE, console = FALSE)
    outgroup <- out$groups
    colnames(outgroup) <- c("trt_mean", "grouping")
    outgroup
    kable(outgroup, caption = "Treatment Grouping", align = "c") %>% kable_styling() %>% kable_paper("hover", full_width = F)
  }
}
tags$br()
############
if (input$submit > 0) {
  d <- as.data.frame(csvfile())
  t <- as.numeric(input$trt)
  response <- d[, input$yield]
  treatment <- d[, input$treatment]
  treatment <- factor(treatment)
  anvaTable <- lm(response ~ treatment)
  result <- as.data.frame(anova(anvaTable))
  if (result[1, 5] > 0.05) {
    return("")
  }
  if (input$req == "tukey" || input$req == "dmrt" || input$req == "lsd") {
    cat("Treatments with same letters are not significantly different")
  }
}

tags$br()
tags$br()
```
report generated from:
**Package: grapesAgri1, Version 1.0.0**


---
title: "Summary Statistics"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE,error = FALSE, warning = FALSE)
library(shiny)
library(pastecs)
library(ggpubr)
library(rmarkdown)
library(knitr)
library(kableExtra)
library(magrittr)
library(summarytools)
library(dplyr)
```
```{r echo=FALSE, error=FALSE, message=FALSE, warning=FALSE, include=TRUE,warn.conflicts=FALSE}

csvfile = reactive({
  csvfile = input$file1
  if (is.null(csvfile)) {
    return(NULL)
  }
  dt = read.csv(csvfile$datapath, header = input$header, sep = ",",check.names = FALSE)
})
# output
if (input$req == "summary") {
  if (input$submit4 > 0) {
    y = subset(csvfile(), select = input$var)
    final =
      summarytools::descr(y) %>%
      summarytools::tb(order = 3) %>%
      knitr::kable(digits = 2, caption = "Summary Statistics") %>%
      kableExtra::kable_styling("bordered", full_width = F) %>%
      kableExtra::collapse_rows(columns = 1, valign = "top")
    final
  }
}

if (input$req == "sumbygrp") {
  if (input$submit5 > 0) {
    y1 = subset(csvfile(), select = input$var)
    y2 = subset(csvfile(), select = input$group)
    final =
      summarytools::stby(y1, y2, descr) %>%
      summarytools::tb(order = 1) %>%
      knitr::kable(digits = 2, caption = "Summary Statistics by Group") %>%
      kableExtra::kable_styling("bordered", full_width = F) %>%
      kableExtra::collapse_rows(columns = 1, valign = "top")
    final
  }
}
```


**package: grapesAgri1, Version 1.0.0**



---
title: "Compare Means: Small Sample Test"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(rmarkdown)
library(pastecs)
library(dplyr)
library(reshape2)
library(knitr)
library(magrittr)
library(kableExtra)

```
## RESULT

```{r, echo = FALSE}
csvfile <- reactive({
    csvfile <- input$file1
    if (is.null(csvfile)){return(NULL)}
    dt <- read.csv(csvfile$datapath, header=input$header, sep=",",check.names = FALSE)
    dt
  })

################## summary one sample ttest
if(input$req1=='ottest'){
      if(input$submit3 > 0){
      grp1<- subset(csvfile(),select=input$dvar)
      final<-grp1
      res <- stat.desc(final)
      result<-as.data.frame(res)
      rownames(result)[rownames(result) == "nbr.val"] <- "Number of Obs."
      rownames(result)[rownames(result) == "nbr.null"] <- "null values"
      rownames(result)[rownames(result) == "nbr.na"] <- "NA"
      round(result,3)
      kable(result,caption = "Summary Statistics",digits=3) %>% kable_styling() %>% kable_paper("hover", full_width = F)
      }
    }
################# summary
if(input$req1 == 'ttest'||input$req1 == 'wttest'||input$req1=='ftest'||input$req1=='pttest'){
      if(input$submit1 > 0||input$submit2 > 0||input$submit4 > 0||input$submit5 > 0){
      grp1<- subset(csvfile(),select=input$dvar)
      grp2<- subset(csvfile(),select=input$ivar)
      final<-cbind(grp1,grp2)
      res <- pastecs::stat.desc(final)
      result<-as.data.frame(res)
      rownames(result)[rownames(result) == "nbr.val"] <- "Number of Obs."
      rownames(result)[rownames(result) == "nbr.null"] <- "null values"
      rownames(result)[rownames(result) == "nbr.na"] <- "NA"
      round(result,3)
      kable(result,caption = "Summary Statistics",digits=3) %>% kable_styling() %>% kable_paper("hover", full_width = F)
      }
    }
########## two sample ttest
 if(input$req1 == 'ttest'){
      if(input$submit1 > 0){
        x<-as.vector(csvfile()[,input$dvar])
        y<-as.vector(csvfile()[,input$ivar])
        t<-t.test(x, y, alternative = input$alt,
                  var.equal = TRUE,paired = FALSE,na.omit=TRUE)
        t_value<-round(t$statistic,3)
        df<-t$parameter
        Pvalue<-round(t$p.value,6)
        alt.Hypothesis<-t$alternative
        result<-cbind(t_value,df,Pvalue,alt.Hypothesis)
        result<-as.data.frame(result)
         kable(result,caption = "Two sample unpaired t-test",row.names = FALSE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
      }
 }
h3("")
########## note 1
  if(input$req1 == 'ttest'){
      if(input$submit1 > 0){
        x<-as.vector(csvfile()[,input$dvar])
        y<-as.vector(csvfile()[,input$ivar])
        t<-t.test(x, y, alternative = input$alt,
                  var.equal = TRUE,paired = FALSE, na.omit=TRUE)

        if(t$p.value<=0.05){
        
        cat("Since P-value is < 0.05 we reject the null hypothesis at 5% level of significance.\n(here null hypothesis is: Population mean of group A = Population mean of group B)\n ")
       
        }
       else if(t$p.value>0.05){
       
        cat("Since P-value is > 0.05 we don't have enough evidence to reject the null hypothesis at 5% level of significance.\n(here null hypothesis is: Population mean of group A = Population mean of group B)\n")
        
         }
      }
  }
#############################################welch ttest
 if(input$req1 == 'wttest'){
      if(input$submit2 > 0){
        x<-as.vector(csvfile()[,input$dvar])
        y<-as.vector(csvfile()[,input$ivar])
        t<-t.test(x, y, alternative = input$alt,
                  var.equal = FALSE,paired = FALSE,na.omit=TRUE)
        t_value<-round(t$statistic,3)
        df<-round(t$parameter,2)
        Pvalue<-round(t$p.value,6)
        alt.Hypothesis<-t$alternative
        result<-cbind(t_value,df,Pvalue,alt.Hypothesis)
        result<-as.data.frame(result)
        result
        kable(result,caption = "Welch unpaired t-test",row.names = FALSE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
      }
 }
h3("")
########### note2
 if(input$req1 == 'wttest'){
      if(input$submit2 > 0){
        x<-as.vector(csvfile()[,input$dvar])
        y<-as.vector(csvfile()[,input$ivar])
        t<-t.test(x, y, alternative = input$alt,
                  var.equal = FALSE,paired = FALSE, na.omit=TRUE)

        if(t$p.value<=0.05){
      cat("Since P-value is < 0.05 we reject the null hypothesis at 5% level of significance.\n(here null hypothesis is: Population mean of group A = Population mean of group B)\n ")
     
        }
        else if(t$p.value>0.05){
         cat("Since P-value is > 0.05 we don't have enough evidence to reject the null hypothesis at 5% level of significance.\n(here null hypothesis is: Population mean of group A = Population mean of group B)\n")
        
        }
      }
    }
############### one sample ttest
if(input$req1 == 'ottest'){
      if(input$submit3 > 0){
        x<-as.vector(csvfile()[,input$dvar])
        t<-t.test(x, y=NULL, alternative = input$alt,
                  mu = input$mu, na.omit=TRUE)
        t_value<-round(t$statistic,3)
        df<-round(t$parameter,2)
        Pvalue<-round(t$p.value,6)
        alt.Hypothesis<-t$alternative
        result<-cbind(t_value,df,Pvalue,alt.Hypothesis)
        result<-as.data.frame(result)
        kable(result,caption = "One sample t-test",row.names = FALSE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
      }
}
h3("")

############### note 3
if(input$req1 == 'ottest'){
      if(input$submit3 > 0){
        x<-as.vector(csvfile()[,input$dvar])
        t<-t.test(x, y=NULL, alternative = input$alt,
                  mu = input$mu, na.omit=TRUE)
        mu<-as.numeric(input$mu)

        if(t$p.value<=0.05){
          cat("Since P-value is < 0.05.\nwe can reject the null hypothesis at 5 percent level of significance.\n")
        }
        else if(t$p.value>0.05){
          cat("Since P-value is > 0.05.\nwe don't have enough evidence to reject null hypothesis at 5 percent level of significance.\n")
        }
      }
    }
############### F test
if(input$req1 == 'ftest'){
      if(input$submit4 > 0){
        x<-as.vector(csvfile()[,input$dvar])
        y<-as.vector(csvfile()[,input$ivar])
        f<-var.test(x, y, ratio = 1,
                 alternative = input$alt)
        F_value<-round(f$statistic,3)
        df<-as.data.frame(f$parameter)
        Pvalue<-round(f$p.value,6)
        alt.Hypothesis<-f$alternative
        result<-cbind(F_value,t(df),Pvalue,alt.Hypothesis)
        result<-as.data.frame(result)
         kable(result,caption = "F test for Homogenity of Variance",row.names = FALSE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
      }
}
h3("")
##### note 4
if(input$req1 == 'ftest'){
      if(input$submit4 > 0){
        x<-as.vector(csvfile()[,input$dvar])
        y<-as.vector(csvfile()[,input$ivar])
        f<-var.test(x, y, ratio = 1,
                    alternative = input$alt)

        if(f$p.value<=0.05){
          cat("Since P-value is < 0.05 we can reject the null hypothesis at 5 percent level of significance.\n(Here null hypothesis: Variances are Homogenous)\n")
        }
        else if(f$p.value>0.05){
          cat("Since P-value is > 0.05 we don't have enough evidence to reject null hypothesis at 5 percent level of significance. Variance is homogenous.\n(Here null hypothesis: Variances are Homogenous)\n")
          
        }
      }
    }
############################# Paired t test
if(input$req1 == 'pttest'){
      if(input$submit5 > 0){
      grp1<- subset(csvfile(),select=input$dvar)
      grp2<- subset(csvfile(),select=input$ivar)
      final<-cbind(grp1,grp2)
      ttest <- t.test(final[,1], final[,2], paired = TRUE, alternative = input$alt)
      tvalue<-round(ttest$statistic,3)
      df<-ttest$parameter
      pvalue<-round(ttest$p.value,6)
      meandiff<-round(ttest$estimate,3)
      result1<-cbind(tvalue,df,pvalue,meandiff)
      row.names(result1)<-NULL
       kable(result1,caption = "Paired t test",row.names = FALSE) %>% kable_styling() %>% kable_paper("hover", full_width = F)
      }
}
h3("")
 if(input$req1 == 'pttest'){
      if(input$submit5 > 0){
        x<-as.vector(csvfile()[,input$dvar])
        y<-as.vector(csvfile()[,input$ivar])
        t<-t.test(x, y, alternative = input$alt,
                  var.equal = FALSE,paired = TRUE, na.omit=TRUE)

        if(t$p.value<=0.05){
          cat("Since P-value is < 0.05 we reject the null hypothesis at 5% level of significance.\n(here null hypothesis is: Population mean of group A = Population mean of group B) \n")
          
        }
        else if(t$p.value>0.05){
         cat("Since P-value is > 0.05 we don't have enough evidence to reject the null hypothesis at 5% level of significance.\n(here null hypothesis is: Population mean of group A = Population mean of group B)\n")
          
        }
      }
    }

```
  
  
    
**Package: grapesAgri1, Version: 1.0.0**



---
title: "Correlation"
output: html_document
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
library(knitr)
library(corrplot)
library(Hmisc)
library(magrittr)
library(kableExtra)
```
## RESULT
CORRELATION ANALYSIS 

```{r, echo = FALSE}
csvfile = reactive({
  csvfile = input$file1
  if (is.null(csvfile)) {
    return(NULL)
  }
  dt = read.csv(csvfile$datapath, header = input$header, sep = ",", check.names = FALSE)
  dt
})

# output
if (input$req1 == "correlation") {
  if (input$submit > 0) {
    a = as.vector(csvfile()[, input$dvar])
    y = as.vector(csvfile()[, input$ivar])
    x = cor.test(a, y,
      method = input$req, conf.level = as.numeric(input$ci),
      alternative = input$alt, exact = FALSE
    )
    t_value = round(x$statistic, 3)
    correlation = round(x$estimate, 3)
    df = x$parameter
    pvalue = round(x$p.value, 3)
    alt.Hypothesis = x$alternative
    result = cbind(correlation, t_value, df, pvalue, alt.Hypothesis)
    nam = x$method
    rownames(result) = nam
    result = as.data.frame(result)
    kable(result, caption = "Correlation Analysis", row.names = FALSE) %>%
      kable_styling() %>%
      kable_paper("hover", full_width = F)
  }
}

if (input$req1 == "correlation") {
  if (input$submit > 0) {
    if (input$req == "pearson") {
      a = as.vector(csvfile()[, input$dvar])
      y = as.vector(csvfile()[, input$ivar])
      x = cor.test(a, y,
        method = input$req, conf.level = as.numeric(input$ci),
        alternative = input$alt, exact = FALSE
      )
      ci = x$conf.int
      ci_nw = melt(ci, value.name = "Lower Limit and Upper limit")
      kable(round(ci_nw, 3), caption = "Confidence Interval", row.names = FALSE) %>%
        kable_styling() %>%
        kable_paper("hover", full_width = F)
    }
  }
}



if (input$req1 == "corrmat") {
  if (input$submit2 > 0) {
    x = as.data.frame(csvfile()[, input$selvar])
    cormat = rcorr(as.matrix(x), type = input$req)
    R = round(cormat$r, 3)
    p = cormat$P
    ## Define notions for significance levels; spacing is important.
    mystars = ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    ")))
    Rnew = matrix(paste(R, mystars, sep = ""), ncol = ncol(x))
    diag(Rnew) = paste(diag(R), " ", sep = "")
    row.names(Rnew) = names(x)
    colnames(Rnew) = names(x)
    kable(Rnew, caption = "Correlation Matrix") %>%
      kable_styling() %>%
      kable_paper("hover", full_width = F)
  }
}
tags$br()
if (input$req1 == "corrmat") {
  if (input$submit2 > 0) {
    cat("*** Correlation is significant at 0.001 level (two tailed) \n** Correlation is significant at 0.01 level (two tailed)\n* Correlation is significant at 0.05 level (two tailed)")
  }
}
tags$br()

if (input$req1 == "corrmat") {
  if (input$submit2 > 0) {
    x = as.data.frame(csvfile()[, input$selvar])
    cormat = rcorr(as.matrix(x), type = input$req)
    correlmat1 = cormat$P
    row.names(correlmat1) = names(x)
    kable(round(correlmat1, 3), caption = "Matrix of P-values") %>%
      kable_styling() %>%
      kable_paper("hover", full_width = F)
  }
}



h3("")



```








**package: grapesAgri1, Version; 1.0.0**



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ttest.R
\name{ttApp}
\alias{ttApp}
\title{t-test and Paired t-test}
\usage{
ttApp()
}
\value{
Nothing
}
\description{
ttApp() function opens up an interactive shiny app which will allow
user to easily perform one sample t-test, unpaired two sample t-test,
unpaired two sample Welch t-test, paired t-test, test for homogeneity of variance (F-test),
and obtain plots like boxplot and paired plot by uploading CSV file.
}
\details{
This app uses \code{t.test} function to calculate t statistic.Descriptive statistics
were calculated using \code{stat.desc} function of \code{pastecs} package.
\code{var.test} function is used for F-test.\code{ggboxplot} function
of \code{ggpubr} package is used to draw boxplot. Paired plot is obtained
using \code{paired} function of package \code{PairedData}.
}
\examples{
if (interactive()) {
  ttApp()
}
}
\references{
\insertRef{R_2021}{grapesAgri1}

\insertRef{shiny_2021}{grapesAgri1}

\insertRef{sw_2021}{grapesAgri1}

\insertRef{dplyr_2021}{grapesAgri1}

\insertRef{ggpubr_2020}{grapesAgri1}

\insertRef{past_2018}{grapesAgri1}

\insertRef{gupta1985statistical}{grapesAgri1}

\insertRef{paired_data2018}{grapesAgri1}
}
\keyword{Testing}
\keyword{Welch}
\keyword{and}
\keyword{boxplot}
\keyword{homogeneity}
\keyword{of}
\keyword{one}
\keyword{paired}
\keyword{plot}
\keyword{sample}
\keyword{t-test}
\keyword{two}
\keyword{unpaired}
\keyword{variance}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr.R
\name{corrApp}
\alias{corrApp}
\title{Correlation Analysis}
\usage{
corrApp()
}
\value{
Nothing
}
\description{
corrApp() function opens up an interactive shiny app which will
allow the user to easily calculate Simple correlation, Correlation Matrix and obtain plots
like correlogram and scatterplot by uploading CSV file.
}
\details{
This app uses \code{cor.test} to calculate correlation. Correlation matrix
is calculated using \code{rcorr} function in \code{Hmisc} package. Correlogram
is obtained using \code{corrplot} function in \code{corrplot} package.
}
\examples{
if (interactive()) {
  corrApp()
}
}
\references{
\insertRef{corrplot2021}{grapesAgri1}

\insertRef{Hmisc_2021}{grapesAgri1}

\insertRef{R_2021}{grapesAgri1}

\insertRef{shiny_2021}{grapesAgri1}

\insertRef{sw_2021}{grapesAgri1}

\insertRef{ggplot_2016}{grapesAgri1}

\insertRef{gupta1985statistical}{grapesAgri1}
}
\keyword{Correlation}
\keyword{Correlogram}
\keyword{Matrix}
\keyword{Scatter}
\keyword{plot}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crd.R
\name{crdApp}
\alias{crdApp}
\title{Completely Randomized Design}
\usage{
crdApp()
}
\value{
Nothing
}
\description{
crdApp() function opens up an interactive shiny app which will allow
the user to perform analysis of completely randomized design with
equal or unequal replications. Multiple comparison tests like LSD,
DMRT and Tukey can be performed. Box-pot and Bar-chart with confidence interval
can be plotted. All these can be achieved by uploading CSV file.
}
\details{
This app uses \code{anova} function of \code{stats} package to
obtain one-way ANOVA.\code{LSD.test},\code{duncan.test} and
\code{HSD.test} functions of \code{agricolae} package is used for
multiple comparison test like LSD,DMRT and Tukey respectively.
\code{ggboxplot} function of \code{ggpubr} package is used for
boxplot.'\code{ggplot} function of \code{ggplot2} is used for
barchart with confidence interval.
}
\examples{
if (interactive()) {
  crdApp()
}
}
\references{
\insertRef{R_2021}{grapesAgri1}

\insertRef{shiny_2021}{grapesAgri1}

\insertRef{sw_2021}{grapesAgri1}

\insertRef{dplyr_2021}{grapesAgri1}

\insertRef{ggpubr_2020}{grapesAgri1}

\insertRef{ggplot_2016}{grapesAgri1}

\insertRef{gupta1985statistical}{grapesAgri1}

\insertRef{tukey1977exploratory}{grapesAgri1}

\insertRef{hmisc_2021}{grapesAgri1}

\insertRef{agricolae_2020}{grapesAgri1}

\insertRef{rcol_2014}{grapesAgri1}

\insertRef{shinycss_2020}{grapesAgri1}

\insertRef{das1979design}{grapesAgri1}
}
\keyword{ANOVA}
\keyword{Barchart}
\keyword{Box}
\keyword{Completely}
\keyword{DMRT}
\keyword{Design}
\keyword{Interval}
\keyword{Multiple}
\keyword{One-way}
\keyword{Randomized}
\keyword{Tests}
\keyword{comparison}
\keyword{confidence}
\keyword{plot}
\keyword{with}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout.R
\name{layoutApp}
\alias{layoutApp}
\title{Field Layout of Designs}
\usage{
layoutApp()
}
\value{
Nothing
}
\description{
layoutApp() function opens up an interactive shiny app which will allow
the user to create field layout of Completely Randomized Design (CRD),
Randomized Complete Block Design (RCBD), Split-plot design, Strip-plot design
and Augmented Randomized complete block design. Layout generated are
random. Field layout in table format can also be prepared for recording
observations from the field. Results can be downloaded in HTML format
}
\details{
This app uses \code{design.crd}, \code{design.rcbd}, \code{design.dau},
\code{design.strip}, \code{design.split} functions of package
\code{agricolae} to generate random layout of designs. Field layout
were plotted using \code{desplot} function in \code{desplot} package.
}
\examples{
if (interactive()) {
  layoutApp()
}
}
\references{
\insertRef{R_2021}{grapesAgri1}

\insertRef{shiny_2021}{grapesAgri1}

\insertRef{sw_2021}{grapesAgri1}

\insertRef{shinycss_2020}{grapesAgri1}

\insertRef{dplyr_2021}{grapesAgri1}

\insertRef{agricolae_2020}{grapesAgri1}

\insertRef{desplot_2020}{grapesAgri1}

\insertRef{magi_2020}{grapesAgri1}

\insertRef{Yihui_Xie_2021}{grapesAgri1}

\insertRef{gupta1985statistical}{grapesAgri1}

\insertRef{das1979design}{grapesAgri1}
}
\keyword{Augmented}
\keyword{Completely}
\keyword{Design}
\keyword{RCBD}
\keyword{Randomized}
\keyword{Split-plot}
\keyword{Strip-plot}
\keyword{block}
\keyword{complete}
\keyword{design}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/desc.R
\name{descApp}
\alias{descApp}
\title{Descriptive statistics and Visualization}
\usage{
descApp()
}
\value{
Nothing
}
\description{
descApp() function opens up an interactive shiny app which will allow
the user to easily calculate Summary Statistics, Summary Statistics by Group, Box plot,
Histogram, Q-Q plot and Shapiro-Wilk's test by uploading CSV file.
}
\details{
This app uses \code{descr} and \code{stby} functions of \code{summarytools}
package (Dominic Comtois, 2021) to calculate summary statistics and
summary statistics by group. \code{knitr} (Yihui Xie,2021) and \code{kableExtra}(Hao Zhu,2021) packages
were used to produce HTML tables. \code{shapiro.test}, \code{qqnorm} and \code{qqline} functions of
\code{stats}  package were used for Test of Homogeneity of variance and obtaining
Q-Q plot. \code{hist} and \code{boxplot} of package \code{graphics} were used
to obtain histogram and boxplot respectively. \code{ggqqplot} of package \code{ggpubr} (Alboukadel Kassambara,2020)
is also used to plot Q-Q plot in the app.
}
\examples{
if (interactive()) {
  descApp()
}
}
\references{
\insertRef{Dominic_Comtois_2021}{grapesAgri1}

\insertRef{Hao_zhu_2021}{grapesAgri1}

\insertRef{Yihui_Xie_2021}{grapesAgri1}

\insertRef{R_2021}{grapesAgri1}

\insertRef{shiny_2021}{grapesAgri1}

\insertRef{sw_2021}{grapesAgri1}

\insertRef{dplyr_2021}{grapesAgri1}

\insertRef{ggpubr_2020}{grapesAgri1}

\insertRef{past_2018}{grapesAgri1}

\insertRef{magi_2020}{grapesAgri1}

\insertRef{gridG_2020}{grapesAgri1}

\insertRef{gupta1985statistical}{grapesAgri1}

\insertRef{tukey1977exploratory}{grapesAgri1}

\insertRef{ggplot_2016}{grapesAgri1}
}
\keyword{box}
\keyword{by}
\keyword{descriptive}
\keyword{group}
\keyword{histogram}
\keyword{plot}
\keyword{q-q}
\keyword{statistics}
\keyword{summary}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbd.R
\name{rbdApp}
\alias{rbdApp}
\title{Randomized Block Design}
\usage{
rbdApp()
}
\value{
Nothing
}
\description{
rbdApp() function opens up an interactive shiny app which will allow
the user to perform analysis of randomized Block design.
Multiple comparison tests like LSD,DMRT and Tukey can be performed.
Box-pot and Bar-chart with confidence interval
can be plotted. All these can be achieved by uploading CSV file.
}
\details{
This app uses \code{anova} function of \code{stats} package to
obtain two-way ANOVA.\code{LSD.test},\code{duncan.test} and
\code{HSD.test} functions of \code{agricolae} package is used for
multiple comparison test like LSD,DMRT and Tukey respectively.
\code{ggboxplot} function of \code{ggpubr} package is used for
boxplot.'\code{ggplot} function of \code{ggplot2} is used for
barchart with confidence interval.
}
\examples{
if (interactive()) {
  rbdApp()
}
}
\references{
\insertRef{R_2021}{grapesAgri1}

\insertRef{shiny_2021}{grapesAgri1}

\insertRef{sw_2021}{grapesAgri1}

\insertRef{dplyr_2021}{grapesAgri1}

\insertRef{ggpubr_2020}{grapesAgri1}

\insertRef{ggplot_2016}{grapesAgri1}

\insertRef{gupta1985statistical}{grapesAgri1}

\insertRef{tukey1977exploratory}{grapesAgri1}

\insertRef{hmisc_2021}{grapesAgri1}

\insertRef{agricolae_2020}{grapesAgri1}

\insertRef{rcol_2014}{grapesAgri1}

\insertRef{shinycss_2020}{grapesAgri1}

\insertRef{das1979design}{grapesAgri1}
}
\keyword{ANOVA}
\keyword{Barchart}
\keyword{Block}
\keyword{Box}
\keyword{DMRT}
\keyword{Design}
\keyword{Interval}
\keyword{Multiple}
\keyword{Randomized}
\keyword{Tests}
\keyword{Two-way}
\keyword{comparison}
\keyword{confidence}
\keyword{plot}
\keyword{with}
