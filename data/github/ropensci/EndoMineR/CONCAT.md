
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


---
output: github_document
---


[![Build Status](https://travis-ci.org/sebastiz/EndoMineR.svg?branch=master)](https://travis-ci.org/sebastiz/EndoMineR)
[![ropensci](https://badges.ropensci.org/153_status.svg)](https://github.com/ropensci/onboarding/issues/153)
[![Coverage status](https://codecov.io/gh/sebastiz/EndoMineR/branch/master/graph/badge.svg)](https://codecov.io/github/sebastiz/EndoMineR?branch=master)


<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
library(pander)
library(EndoMineR)
library(here)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```


This package has undergone a major revision to make it much more user friendly. THe documentation has been updated to reflect this. I am always happy to hear of any feedback, positive and negative.

## **Aims of EndoMineR**

The goal of EndoMineR is to extract as much information as possible from free or semi-structured endoscopy reports and their associated pathology specimens.  A full tutorial can be found [here](https://docs.ropensci.org/EndoMineR/articles/EndoMineR.html)

## Installation

You can install EndoMineR from github with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropenSci/EndoMineR")
```

If you dont have access to github, then download the zip and change the working dirctory to the place you have downloaded it, then do

```{r gh-installation2, eval = FALSE}
setwd("C:/Users/Desktop/")

#On windows you cand cd to change the directory or us pushd to create a temporary directory indtead of cd and then setwd to the temporary directory
unzip("EndoMineR.zip")
file.rename("EndoMineR.zip-master", "EndoMineR.zip")
shell("R CMD build EndoMineR.zip")

#Then install the resulting tarball with:

install.packages("EndoMineR_0.2.0.9000.tar.gz", repos = NULL)
```


### How to contribute

Contributions to this project are most welcome. There are just a few small guidelines you need to follow.

#### Submitting a patch

It's generally best to start by opening a new issue describing the bug or feature you're intending to fix. Even if you think it's relatively minor, it's helpful to know what people are working on. Mention in the initial issue that you are planning to work on that bug or feature so that it can be assigned to you.

Follow the normal process of forking the project, and setup a new branch to work in. It's important that each group of changes be done in separate branches in order to ensure that a pull request only includes the commits related to that bug or feature.

The best way to ensure your code is properly formatted is to use lint. Various packages in R provide this.

Any significant changes should almost always be accompanied by tests. The project already has good test coverage, so look at some of the existing tests if you're unsure how to go about it. 

Do your best to have well-formed commit messages for each change. This provides consistency throughout the project, and ensures that commit messages are able to be formatted properly by various git tools.

Finally, push the commits to your fork and submit a pull request. Please, remember to rebase properly in order to maintain a clean, linear git history.


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)



---
title: "Analysis"
author: "Sebastian Zeki"
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    toc_depth: 5
    css: style.css
vignette: >
  %\VignetteIndexEntry{Analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
library(pander)
library(EndoMineR)
library(dplyr)
knitr::opts_chunk$set(echo = TRUE)

```

```{r global_options, include=FALSE}
knitr::opts_chunk$set( warning=FALSE, message=FALSE)
```

<br>
A fundamental design principle of EndoMineR was that it should address the important categories of questions we all have in gastroenterology, and endoscopy in particular. These questions roughly fall into the following: surveillance, quality and also operational questions (eg patient flow through endoscopy).



##**1. Surveillance functions**

<br>
Surveillance tracking is difficult because it relies on assessment at several timepoint and then deciding on the next examination based on a ruleset. A basic question is often: 'How good are our surveillance programmes?" which really means "How good are we at making sure patients come back in a timely way for their endoscopy after a polyp removal or for Barrett's surveillance", for example?

<br>

```{r fig.width=12, fig.height=8,fig.align='center',echo=FALSE}
knitr::include_graphics("img/EndoMineR_Surveillance.svg")
```

<br>

Surveillance relates to the timing of a test relative to other tests or all tests done for a patient. To do this, the EndoMineR surveillance functions simply order the endoscopies by patient and date, and extract the date the first test was done, as well as the last test (of the same type) and the difference in timing between each test, always grouped by patient. 

<br>
As all these functions are simply looking at the date of the test, they can take a raw dataset, as long as a date column is present and use that, rather than have a lot of pre-processing steps. Of course, the pre-processing steps explained in the EndoMineR vignette (mainly using the **textPrep** function) are recommended however as then the user will be able to perform any other additional analyses if needed.

<br>
The basic surveillance functions are simple but are the most used. **SurveilTimeByRow** will extract the time difference between each individual endoscopy for an individual patient. This is useful to see how adherent the surveillance endoscopy is to guidelines. 

<br>

**SurveilLastTest** simply extracts the last and first test respectively for each patient so you can assess how long the patient has been surveilled for. This is likely to come in useful for future iterations of EndoMineR as patient Theographs are developed (work in progress). 

<br>

```{r exampleSurveillanceTimeByRow, eval = TRUE}
em1<-SurveilTimeByRow(Myendo,'HospitalNumber','Dateofprocedure')
```

<br>

```{r exampleSurveillanceTimeByRowtbl, echo = FALSE}
pander(head(data.frame(em1[2],em1[ncol(em1)]),5))

```

<br>

```{r exampleSurveilLastTest, echo = TRUE}
em3<-SurveilLastTest(Myendo,'HospitalNumber','Dateofprocedure')
```

<br>

```{r exampleSurveilLastTesttbl, echo = FALSE}
pander(head(data.frame(em3[2],em3[5]),5))
```

<br>

```{r exampleSurveilFirstTest, echo = TRUE}
em4<-SurveilFirstTest(Myendo,'HospitalNumber','Dateofprocedure')
```

<br>

```{r exampleSurveilFirstTesttbl, echo = FALSE}
pander(head(data.frame(em4[2],em4[5]),5))
```

<br>

Of course we may also want to know how many tests have been done over a time period and this is provided by the function **HowManyTests**

<br>

This function will return the number of tests by day, month and year so they can be easily graphed according to what you want.

<br>

```{r exampleHowManyTests, echo = TRUE}
how<-HowManyOverTime(Myendo,'Indications','Dateofprocedure','Surv')
```

<br>

```{r exampleHowManyTeststbl, echo = FALSE}
pander(head(data.frame(how),5))
```



## **2. Assessment of quality functions**

<br>

Quality is measured in a variety of ways. For endoscopy it is measured according to the adherence to a) standards for endoscopic documentation as well as b) detection of certain pathological conditions such as dysplasia (best summarised as lesion recognition)


<br>

### **a) Documentation Quality**

<br>

As regards adherence to documentation for example, a generic function is provided that will look up the presence of words presented in a list in a target column. It will then output the proportion of reports that have these words, as well as a barchart to show what proportion of the endoscopies showed these words. The list can be comprised of terms that should be mentioned in a report.

<br>

_Input_
                      
```{r exampleListLookup, eval = TRUE,echo=FALSE}
panderOptions('table.split.table', Inf)
pander(head(data.frame(Myendo[2:3],Myendo[13])))
```

<br>

In this example we are looking for the words Barrett's and coeliac as perhaps we have chosen the macroscopic recognition of these features to denote what an endoscopist should always describe in the endoscopy report:

<br>

```{r exampleListLookup2, eval = TRUE}
library(tm)
myNotableWords <- c("barrett", "coeliac")
ListLookup(Myendo,'Findings',myNotableWords)
```

<br>

```{r exampleListLookup3, echo=FALSE}
#pander::panderOptions('table.split.table', Inf)
#pander(head(tt))
```

<br>

So we can see that the terms are present in the minority of reports across endoscopists, so perhaps we can look into this further..

<br>

### **b) Endoscopic Quality**

#### **Sedation Usage**

<br>

Another measure of quality is the assessment of those factors that are recorded at endoscopy such as degree of sedation used etc. Rather than provide a function for each metric, again a generic function is provided that uses any quantifiable metric and plots it against the endoscopist. This function returns a list with two elements- the plot and the table:

<br>

```{r exampleEndoscChopperMeds,echo=TRUE}
#We have to attach the output of EndoscMeds to the original dataframe
MyendoNew<-cbind(EndoscMeds(Myendo$Medications),Myendo)
#Average Fentanyl use by endoscopist:
Mytable<-MetricByEndoscopist(MyendoNew,'Endoscopist','Fent')
```

<br>

```{r exampleEndoscChopperMedstbl,echo=FALSE}

pander(head(data.frame(MyendoNew$Endoscopist,MyendoNew$Fent),10))
```

<br>
<br>



<br>


## **4.Patient flow functions**

<br>

### **Sankey plots**

<br>

We often like to get an overview of how patients are flowing through a system overall. This can give a nice visual representation of whether which patients diverge from the normal flow through a system so we can study them further. There are two ways to look at this. Sankey plots give good timepoint by timepoint representation of flow. This really works with more than one type of event at each timepoint. 

<br>

For example, if we have a dataset with events such as 'radiofrequency ablation' and 'endoscopic mucosal resection' or 'nothing' we can use the Sankey plot to determine the order of events over a large patient population. You choose the column in the dataframe that describes the Procedure type ("EMR","RFA","nothing" in this case)

<br>

```{r exampleSurveySankey, eval = FALSE}
#how<-SurveySankey(Myendo,"ProcedurePerformed")
```

<br>



```{r fig.width=12, fig.height=8,fig.align='center',echo=FALSE,out.width = "100%"}
knitr::include_graphics("img/EndoMineR_Sankey.svg")
```

<br>


### **Circos plots**

<br>

We may need something even more aggregated. Perhaps we want to see the overall number of patients that go from one event to another regardless of which timepoint it is at. To do this we can use a circos plot, which makes use of the circlize library, as follows:

```{r examplePatientFlow_CircosPlots, eval = FALSE}
#flow<-PatientFlow_CircosPlots(v,"Date.y","pHospitalNum","ProcedurePerformed")
```

<br>


```{r fig.width=12, fig.height=8,fig.align='center',echo=FALSE,out.width = "60%"}
knitr::include_graphics("img/EndoMineR_Circos.svg")
```

<br>---
title: "Polyps"
author: "Sebastian Zeki"
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    toc_depth: 5
    css: style.css
vignette: >
  %\VignetteIndexEntry{Polyps}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
library(pander)
library(EndoMineR)
knitr::opts_chunk$set(echo = TRUE)

```

```{r global_options, include=FALSE}
knitr::opts_chunk$set( warning=FALSE, message=FALSE)
```


##**Pathology Detection Quality**

<br>

###**Polyp and sample location**

<br>

A difficult area is the assessment of endoscopic quality by looking at the pathology processed from an endoscopy. This package is excellent at dealing with this kind of question because of its ability to merge the datasets together:

A particularly well developed area to look at is that of the [Global Rating Scale](http://www.healthcareimprovementscotland.org/our_work/governance_and_assurance/endoscopy_accreditation/grs_reports.aspx) for assessing the quality of colonoscopy. One of the metrics- the adenoma detection rate assesses the number of colonoscopies where at least one adenoma was detected.

One function is provided to produce a table that gives the number of adenomas, adenocarcinomas and hyperplastic polyps (also as a ration to adenomas) by endoscopist therefore immediately fulfilling the GRS requirement for the ADR as well as providing further metrics alongside.

Having made sure the text is prepared (using the **textPrep** function) and merged with Pathology reports (which have also been textPrep'd) we can go straight ahead and use the function **GRS_Type_Assess_By_Unit** on one of our example datasets called vColon (the format of the dataset can be seen in the Data vignette. Here we will use the columns "ProcedurePerformed","Endoscopist", "Diagnosis" and  "Histology" to calculate the adenoma detection rate, as well as the rate of detection of adenocarcinoma, high grade dysplasia, low grade dysplasia, serrated adenomas and hyperplastic adenomas.


```{r exampleGRS}
  data(vColon)
  nn<-GRS_Type_Assess_By_Unit(
  vColon, "ProcedurePerformed",
  "Endoscopist", "Diagnosis", "Histology")
```

```{r exampleGRStbl,echo=FALSE}
   pander(nn)
```

<br>
 
This of course can then be graphed using the pre-formatted function **EndoBasicGraph** 
```{r exampleGRSGraph,warning=FALSE,message=FALSE}
EndoBasicGraph(nn, "Endoscopist", "Adenocarcinoma")
```
---
title: "Barrett's Oesophagus"
author: "Sebastian Zeki"
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    css: style.css
vignette: >
  %\VignetteIndexEntry{Barretts}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



```{r setup, include=FALSE}

library(knitr)
library(EndoMineR)
library(pander)
library(prettydoc)
```

```{r global_options, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
```

<br>

## **Specific diseases - Barrett's oesophagus**

<br>

One particular disease that lends itself well to analytics, particularly as it is part of a surveillance programme, is the premalignant oesophageal condition Barrett's oesophagus. This is characterised by the growth of cells (called columnar lined epithelium) in the oesophagus. These cells usually occupy the lower part of the oesophagus as a continuous sheet from the top of the stomach to varying lengths up the oesophagus. 

This condition requires endoscopic surveillance and the timing of this depends on the prior endoscopic features (namely the length of the Barretts segment as measured by the Prague score- explained below) and the pathological stage at that endoscopy (which for non-dysplastic samples, since the revised 2013 UK guidelines, means the presence or absence of intestinal metaplasia). This can be seen in the image below (from Fitzgerald RC, et al. Gut 2013;0:1–36. doi:10.1136/gutjnl-2013-305372)

<br>

```{r ,fig.align='center',echo=FALSE,fig.width=18, fig.height=12}
knitr::include_graphics("img/BarrettsGuide.png")
```

<br>

## **1. Pre-processing Barrett's samples**

<br>

Such a dataset needs some processing prior to the analysis so for this we can turn to a specific set of function for Barrett's oesophagus itself. 

<br>

### **a) Prague score**

<br>

Firstly we need to extract the length of the Barrett’s segment. This is known as the [Prague score](https://services.nhslothian.scot/endoscopyunit/InformationForClinicalStaff/Documents/PRAGUE%20CRITERIA.pdf) and is made up of the length from the top of the gastric folds (just below the gastro-oesophageal junction) to the top of the circumferential extent of the Barrett's segment (C). In addition the maximal extent is from the top of the gastric folds to the top of the tongues of Barrett's segment (M). This gives an overall score such as C1M2.

<br>

After filtering for endoscopic indication (eg “Surveillance-Barrett’s”- this is stored in the 'Indication' column in our data set) the aim of the following function is to extract a C and M stage (Prague score) for Barrett’s samples. This is done using a regular expression where C and M stages are explicitly mentioned in the free text. Specifically it extracts the Prague score. This is usually mentoned in the 'Findings' column in our dataset but obviously the user can define which column should be searched.

```{r exampleBarretts_PragueScore,echo = TRUE,message=FALSE, warning=FALSE}
v<-Barretts_PragueScore(Myendo,'Findings','OGDReportWhole')
```

```{r exampleBarretts_PragueScore2,echo = FALSE,message=FALSE, warning=FALSE}
panderOptions('table.split.table', Inf)
pander(v[23:27,(ncol(v)-4):ncol(v)])
```

### **b) Worst pathological stage**

<br>

We also need to extract the worst pathological stage for a sample, and if non-dysplastic, determine whether the sample has intestinal metaplasia or not. This is done using 'degredation' so that it will look for the worst overall grade in the histology specimen and if not found it will look for the next worst and so on. 

It looks per report not per biopsy (it is more common for histopathology reports to contain the worst overall grade rather than individual biopsy grades).


```{r exampleBarretts_PathStage, echo = TRUE,message=FALSE, warning=FALSE}
#The histology column is the one we are interested in:
Mypath$b <- Barretts_PathStage(Mypath, "Histology")
```

```{r exampleBarretts_PathStage2, echo = TRUE,message=FALSE, warning=FALSE,echo=FALSE}
panderOptions('table.split.table', Inf)
pander(Mypath[2:3,(ncol(Mypath)-4):ncol(Mypath)])
```



### **c)Follow-up groups**

<br>

Having done these pre-processing steps, the follow-up group to which the last endoscopy belongs (rather than the patient as their biopsy results or Barrett's segment length and therefore their follow-up timing, may fluctuate over time) can be determined. 

The follow-up timing, as explained in the the original guideline flowchart above, depends on the length of the Barrett's segment and the presence of intestinal metaplasia (a type of columnar lined epithelium). If abnormal cells (dysplasia) are present the there is a different follow-up regime which we won't concern ourselves with at the moment. 

The timing of follow-up is done with the function **Barretts_FUType.** This relies on the previous functions called **Barretts_PathStage**  and **Barretts_PragueScore** having been run. The Barretts_FUType function will tell you which follow up Rule the patient should be on so that the timing of the next endoscopy can be determined. As these functions usually go together a wrapper function called **BarrettAll** is also provided.

<br>

```{r exampleBarretts_FUType, echo = TRUE,message=FALSE, warning=FALSE}
#Create the merged dataset
v<-Endomerge2(Myendo,"Dateofprocedure","HospitalNumber",Mypath,"Dateofprocedure","HospitalNumber")
#Find the worst pathological grade for that endoscopy
v$IMorNoIM <- Barretts_PathStage(v, "Histology")
#Find the Prague score for that endoscopy
b1<-Barretts_PragueScore(v, "Findings", "OGDReportWhole")
#Get the follow-up type for that endoscopy
b1$FU_Type<-Barretts_FUType(b1,"CStage","MStage","IMorNoIM")

```


```{r exampleBarretts_FUType2, echo = FALSE,message=FALSE, warning=FALSE,echo=FALSE}
panderOptions('table.split.table', Inf)
pander(b1[23:27,(ncol(b1)-4):ncol(b1)])
```





## **2.Quality assessment in Barrett's surveillance**

<br>

Many of the aspects of generic quality monitoring apply to Barrett's as well. For example, you may want to make sure that the endoscopy reports adhere to guidance about what should be in an endoscopy report. This can be assessed using the generic function **ListLookUp**




### **Quality of perfomance of Barrett's surveillance endoscopies as just by tissue sampling**

Some things are specific to Barrett's oesophagus however. One of the essential requirements to demonstrate adequate sampling of Barrett's oesophagus during endoscopy is that the endoscopist should adhere to the 'Seattle protocol' for biopsies which is to take  4 equally spaced biopsies at 2cm intervals in the circumferential part of the oesophagus. Because the macroscopic description of the pathological specimen tells us how many samples are taken overall (and rarely how many at each level but this is usually not the case for a variety of reasons) we can determine the shortfall in the number of biopsies taken, per endoscopist. Again pre-processing the Barrett's samples is pre-requisite. The Number of biopsies and their size should also be extracted using the histopathology functions.

```{r exampleBarrettsQuality_AnalysisBiopsyNumber,echo = FALSE,message=FALSE, warning=FALSE}
 # The number of average number of biopsies is then calculated and
 # compared to the average Prague C score so that those who are taking
 # too few biopsies can be determined
#Lets just use the Surveillance endoscopies:
b1<-b1[grepl("[Ss]urv",b1$Indications),]
b1$NumBx<-HistolNumbOfBx(b1$Macroscopicdescription,'specimen')
b1$BxSize <- HistolBxSize(b1$Macroscopicdescription)
b2<-BarrettsBxQual(b1,'Date.x','HospitalNumber',
                                    'Endoscopist')
```
 
 
```{r exampleBarrettsQuality_AnalysisBiopsyNumbertbl, echo = FALSE,message=FALSE, warning=FALSE,echo=FALSE}

#panderOptions('table.split.table', Inf)
panderOptions('table.split.table', Inf)
pander(b2)
```
 <br>
 
 This function will again return a dataframe with the number of biopsies taken that is outside of the number that should have been taken for a certain Prague score.
---
title: "Data"
author: Sebastian Zeki
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    toc_depth: 5
    css: style.css
vignette: >
  %\VignetteIndexEntry{Data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```{r setup, include=FALSE}
library(pander)
library(EndoMineR)
knitr::opts_chunk$set(echo = TRUE)
```

```{r global_options, include=FALSE}
knitr::opts_chunk$set( warning=FALSE, message=FALSE)
```

##Overview

It is envisaged that different users will start at different points of their data preparation. THis section is intended to explain the fake data I have created so the type of data used for the examples can be better understood.

There are several data files used. These are detailed below

## Gastroscopy

####  Raw datasets:

#####TheOGDReportFinal 
A dataset containing fake upper GI endoscopy reports. The report field is provided as a whole report without any fields having been already extracted. There are 2000 rows

#####PathDataFrameFinal 
A dataset containing fake upper GI pathology reports. The report field is provided as a whole report without any fields having been already extracted. There are 2000 rows

####Pre-extracted datasets:

#####Myendo 
This has been extracted using the Extractor method as follows from the raw text within Mypath:

```{r MyendoExtract, eval = FALSE}
mywords <- c("OGDReportWhole","HospitalNumber","PatientName",
             "GeneralPractitioner","Dateofprocedure","Endoscopist",
             "Secondendoscopist","Medications","Instrument","ExtentofExam",
             "Indications","ProcedurePerformed","Findings")
Extractor(TheOGDReportFinal,"OGDReportWhole",mywords)
```

####Mypath 

This has been extracted using the Extractor method as follows from the raw text within Mypath:

```{r MypathExtract, eval = FALSE}
mywords<-c("HospitalNumber","PatientName","DOB","GeneralPractitioner",
           "Dateofprocedure","ClinicalDetails","Macroscopicdescription",
           "Histology","Diagnosis")
Extractor(PathDataFrameFinal,"PathReportWhole",mywords)
```

The original dataset has also been added as "PathReportWhole",

##Colonoscopy

###Raw datasets

####ColonFinal

A dataset containing fake lower GI endoscopy reports. The report field is provided as a whole report without any fields having been already extracted. There are 2000 rows

####PathDataFrameFinalColon

A dataset containing fake lower GI pathology reports. The report field is provided as a whole report without any fields having been already extracted. There are 2000 rows.
---
title: "EndoMineR"
author: "Sebastian Zeki `r Sys.Date()`"
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    tod_depth: 4
    css: style.css

vignette: >
  %\VignetteIndexEntry{EndoMineR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
library(pander)
library(EndoMineR)
library(here)
knitr::opts_chunk$set(echo = TRUE)
```

```{r global_options, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
```
## **Aims of EndoMineR**

The goal of EndoMineR is to extract as much information as possible from free or semi-structured endoscopy reports and their associated pathology specimens. 

Gastroenterology now has many standards against which practice is measured although many reporting systems do not include the reporting capability to give anything more than basic analysis. Much of the data is locked in semi-structured text. However the nature of semi-structured text means that data can be extracted in a standardised way- it just requires more manipulation. 

This package provides that manipulation so that complex endoscopic-pathological analyses, in line with recognised standards for these analyses, can be done.






##**How is the package divided?**

<br>

```{r fig.width=8, fig.height=8,fig.align='center',out.width = "100%",echo=FALSE}
knitr::include_graphics("img/EndoMineRBasic.jpg")
```

<br>

The package is basically divided into  parts. How all the parts are connected in shown in the figure above. 
The import of the raw data is left up to the user with the overall aim being that all the data is present in one dataframe. 


## **1. Extraction and cleaning**
If using raw reports, then once you have imported the data into your R environment, you can go ahead and use the **textPrep** function. If you have pre-segregated data, then use the **EndoPaste** function (explained at the end of this section). 

Once the data is ready you can use the **textPrep** function. This function is really a wrapper for a number of different functions that work on your data. 

#### **i. Dusting the data** 
Firstly the data is cleaned up so that extra newlines, spaces, unnecessary punctuation etc is dealt with. Although you dont need to know the details of this, if you want to look under the hood you can see that this is part of the **ColumnCleanUp** function.

#### **ii. Spell checking** 
The textPrep also implements a spell checker (the **spellCheck** function). This really checks the spelling of gastroenterology terms that are present in the in-built lexicons. So, if for example the report contains the terms 'Radiafrequency ablashion' then this will be corrected using Radiofrequency ablation. This is one function that acts to standardise the text in the report so that downstream analyses can be robust. 

#### **iii.Negative phrase removal** 
The text is then passed along to the **NegativeRemove** function. This will remove any phrases that contain negative sentences indicating a non-positive finding. For example it would remove 'There is no evidence of dysplasia here'. This makes text extraction analyses, which often report the detection of a positive finding, much more accurate. If you wish to include negative phrases however, you can. This function is part of the parameters you can switch on and off for the **textPrep** function.

#### **iv. Term mapping** 
The next step is to perform term mapping. This means that variations of a term are all mapped to a single common term. For example 'RFA' and 'HALO' may be both mapped to 'Radiofrequency ablation'. This is performed using the using the **DictionaryInPlaceReplace** function. This function looks through all of the lexicons included in this package to perform this. The lexicons are all manually created and consist of key-value pairs which therefore map a key which is the term variant) to a value (which is the standardised term that should be used).

#### **Sv. egregating the data** 
  Finally the text is ready to separate into columns. The basic **Extractor** function will take your data with the list of terms you have supplied that act as the column boundaries and separate your data accordingly. It is up to the user to define the list of boundary keywords. This list is made up of the words that will be used to split the document. It is very common for individual departments in both gastroenterology and pathology to use semi-structured reporting so that the same headers are used between patient reports. The Extractor then does the splitting for each pair of words, dumps the text between the delimiter in a new column and names that column with the first keyword in the pair with some cleaning up and then returns the new dataframe.



### **What if my data is already in columns?**

Often the raw data is derived from a basic query to some endoscopy software. The output is likely to be a spreadsheet with some of the data (eg Endoscopist) already placed in separate columns. The report still needs to be tidied up in order to maximise what we can get from the data and especially the free text. In order to do the term mapping and negative removal etc it is therefore necessary for EndoMineR to merge all of the columns together, for each endoscopy, but it keeps the column headers as delimiters for use later. This is the case with the top left box in the first diagram.

**EndoPaste** is provided as an optional function to get your data in the right shape to allow EndoMineR to process the data properly in this situation

**EndoPaste** outputs two things: the data, as well as a list of delimiters (basically your column headings). The delimiters can then be used in the **textPrep** function. This would work as follows:

```{r exampleEndoPaste,echo=TRUE}

#An example dataset
testList<-structure(list(PatientName = c("Tom Hardy", "Elma Fudd", "Bingo Man"), 
                         HospitalNumber = c("H55435", "Y3425345", "Z343424"), Text = c("All bad. Not good", 
"Serious issues", "from a land far away")), class = "data.frame", row.names = c(NA, -3L))
myReadyDataset<-EndoPaste(testList)


```


The dataframe can be obtained from myReadyDataset[1] and the delimiters from myReadyDataset[2]


<br>
<br>
<br>

### **Can I have an example please?**

As an example, here we use an sample dataset (which has not had separate columns selected already, so no need to use EndoPaste) as the input. This contains a synthetic dataset provided as part of the EndoMineR package.

```{r exampleExtractor,echo=TRUE}
PathDataFrameFinalColon2<-PathDataFrameFinalColon
names(PathDataFrameFinalColon2)<-"PathReportWhole"
pander(head(PathDataFrameFinalColon2,1))
```

```{r exampleExtractortbl,echo=FALSE}
pander(head(PathDataFrameFinalColon2,1))
```
<br>

We can then define the list of delimiters that will split this text into separate columns, title the columns according to the delimiters and return a dataframe. each column simply contains the text between the delimiters that the user has defined. These columns are then ready for the more refined cleaning provided by subesquent functions.

<br>


```{r exampleExtractor2,echo=TRUE}
library(EndoMineR)
mywords<-c("Hospital Number","Patient Name:","DOB:","General Practitioner:",
"Date received:","Clinical Details:","Macroscopic description:",
"Histology:","Diagnosis:")

PathDataFrameFinalColon3<-Extractor(PathDataFrameFinalColon2$PathReportWhole,mywords)

```

```{r exampleExtractor3,echo=FALSE}

panderOptions('table.split.table', Inf)
pander(head(PathDataFrameFinalColon3[,1:9],1))
```

The **Extractor** function is embedded within **textPrep** so you may never have to use it directly, but you will always have to submit a list of delimiters so that **textPrep** can use the **Extractor** to do its segregation:

```{r exampletextPrep,echo=FALSE}
#Submit delimiters
mywords<-c("Hospital Number","Patient Name:","DOB:","General Practitioner:",
"Date received:","Clinical Details:","Macroscopic description:",
"Histology:","Diagnosis:")
CleanResults<-textPrep(PathDataFrameFinal$PathReportWhole,mywords)
```

**textPrep** takes the optional parameters NegEx which tells the **textPrep** to remove negative phrases like 'There is no pathology here' from the text.





There will always be a certain amount of data cleaning that only the end user can do before data can be extracted. There is some cleaning that is common to many endoscopy reports and so functions have been provided for this. An abvious function is the cleaning of the endoscopist's name. This is done with the function **EndoscEndoscopist**. This removes titles and tidies up the names so that there aren't duplicate names (eg "Dr. Sebastian Zeki" and "Sebastian Zeki"). This is applied to any endoscopy column where the Endoscopist name has been isolated into its own column.

```{r exampleEndoscEndoscopist2,echo = TRUE}
EndoscEndoscopist(Myendo$Endoscopist[2:6])
```



## **2. Data linkage**

Endoscopy data may be linked with other types of data. The most common associated dataset is pathology data from biopsies etc taken at endoscopy. This pathology data should be processed in exactly the same way as the endoscopy data- namely with **textPrep** (with or without **EndoPaste**). 

The resulting pathology dataset should then be merged with the endoscopy dataset using **Endomerge2** which will merge all rows with the same hospital number and do a fuzzy match (up to 7 days after the endoscopy) with pathology so the right endoscopy event is associated with the right pathology report.An example of merging the included datasets Mypath and Myendo is given here:

```{r exampleEndomerge2,echo=TRUE}
v<-Endomerge2(Myendo,'Dateofprocedure','HospitalNumber',Mypath,'Dateofprocedure','HospitalNumber')
```

## **3. Deriving new data from what you have**

Once the text has been separated in to the columns of your choosing you are ready to see what we can extract. Functions are provided to allow quite complex extractions at both a general level and also for specific diseases.




#### **i) Medication**

The extraction of medication type and dose is important for lots of analyses in endoscopy. This is provided with the function **EndoscMeds**. The function currently extracts Fentanyl, Pethidine, Midazolam and Propofol doses into a separate column and reformats them as numeric columns so further calculations can be done. It outputs a dataframe so you will need to re-bind this output with the original dataframe for further analyses. This is shown below:


```{r exampleEndoCleaningFuncMed, echo = TRUE}
MyendoMeds<-cbind(EndoscMeds(Myendo$Medications), Myendo)
```

```{r exampleEndoCleaningFuncMedtbl, echo = FALSE}
pander(head(MyendoMeds[1:4],5))

```
#### **ii) Endosccopy Event extraction**

The **EndoscopyEvent** function will extract any event that has been performed at the endoscopy and dump it in a new column. It does this by looking in pairs of sentences and therefore wraps around a more basic function called **EndoscopyPairs_TwoSentence** which you dont have to directly use. This allows us to get the site that the event (usually a therapy) happened on. A example sentence might be 'There was a bleeding gastric ulcer. A clip was applied to it' We can only extract stomach:clip by reference to the two sentences.


#### **iii) Histology biopsy number extraction**

There is also a lot of information we can extract into a separate column from the histopathology information.
We can derive the number of biopsies taken (usually specified in the macroscopic description of a sample) using the function **HistolNumbOfBx**. 

In order to extract the numbers, the limit of what has to be extracted has to be set as part of the regex so that the function takes whatever word limits the selection.It collects everything from the regex [0-9]{1,2}.{0,3} to whatever the string boundary is.

For example, a report might say:

<br>

```{r exampleHistolNumbOfBx1, echo = TRUE}
sg<-data.frame(Mypath$HospitalNumber,Mypath$PatientName,Mypath$Macroscopicdescription)
```

```{r exampleHistolNumbOfBx1tbl, echo = FALSE}
pander(head(sg,5))
```
<br>

Based on this, the word that limits the number you are interested in is 'specimen' so the function and it's output is:

<br>

```{r exampleHistolNumbOfBx2, echo = TRUE}
Mypath$NumbOfBx<-HistolNumbOfBx(Mypath$Macroscopicdescription,'specimen')
sh<-data.frame(Mypath$HospitalNumber,Mypath$PatientName,Mypath$NumbOfBx)
```

```{r exampleHistolNumbOfBx2tbl, echo = FALSE}
pander(head(sh,5))
```

#### **iv) Histology biopsy size extraction**

We might also be interested in the size of the biopsy taken. A further function called **HistolBxSize** is provided for this. This is also derived from the macroscopic description of the specimen

```{r exampleHistolBxSize1, echo = TRUE}
Mypath$BxSize<-HistolBxSize(Mypath$Macroscopicdescription)
sh<-data.frame(Mypath$HospitalNumber,Mypath$PatientName,Mypath$BxSize)
```

```{r exampleHistolBxSize1tbl, echo = FALSE}

pander(head(sh,5))
```
#### **v) Histology type and site extraction**

A final function is also provided called **HistolTypeAndSite**. This is particularly useful when trying to find out which biopsies came from which site. 
The output is provided as a site:specimen type pair. An alternative output is also provided which groups locations (eg the gastro-oesophageal junction is seen as part of the oesophagus). 

##**And this is just the beginning...**

What is described above are the building blocks for starting a more complex analysis of the Endoscopic and Pathological data. 

Generic analyses, such as figuring out surveillance intervals, can be determined using the appropriate functions in the Analysis module.

The dataset can also be fed in to more complex functions as are described in the Barrett's and Polyp modules

The package also provides generic data visualisation tools to assess Patient flow, amonst other visualisations. All of these can be found in the associated vignettes.
---
title: "Analysis"
author: "Sebastian Zeki"
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    toc_depth: 5
    css: style.css
vignette: >
  %\VignetteIndexEntry{Analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
library(pander)
library(EndoMineR)
library(dplyr)
knitr::opts_chunk$set(echo = TRUE)

```

```{r global_options, include=FALSE}
knitr::opts_chunk$set( warning=FALSE, message=FALSE)
```

<br>
A fundamental design principle of EndoMineR was that it should address the important categories of questions we all have in gastroenterology, and endoscopy in particular. These questions roughly fall into the following: surveillance, quality and also operational questions (eg patient flow through endoscopy).



##**1. Surveillance functions**

<br>
Surveillance tracking is difficult because it relies on assessment at several timepoint and then deciding on the next examination based on a ruleset. A basic question is often: 'How good are our surveillance programmes?" which really means "How good are we at making sure patients come back in a timely way for their endoscopy after a polyp removal or for Barrett's surveillance", for example?

<br>

```{r fig.width=12, fig.height=8,fig.align='center',echo=FALSE}
knitr::include_graphics("img/EndoMineR_Surveillance.svg")
```

<br>

Surveillance relates to the timing of a test relative to other tests or all tests done for a patient. To do this, the EndoMineR surveillance functions simply order the endoscopies by patient and date, and extract the date the first test was done, as well as the last test (of the same type) and the difference in timing between each test, always grouped by patient. 

<br>
As all these functions are simply looking at the date of the test, they can take a raw dataset, as long as a date column is present and use that, rather than have a lot of pre-processing steps. Of course, the pre-processing steps explained in the EndoMineR vignette (mainly using the **textPrep** function) are recommended however as then the user will be able to perform any other additional analyses if needed.

<br>
The basic surveillance functions are simple but are the most used. **SurveilTimeByRow** will extract the time difference between each individual endoscopy for an individual patient. This is useful to see how adherent the surveillance endoscopy is to guidelines. 

<br>

**SurveilLastTest** simply extracts the last and first test respectively for each patient so you can assess how long the patient has been surveilled for. This is likely to come in useful for future iterations of EndoMineR as patient Theographs are developed (work in progress). 

<br>

```{r exampleSurveillanceTimeByRow, eval = TRUE}
em1<-SurveilTimeByRow(Myendo,'HospitalNumber','Dateofprocedure')
```

<br>

```{r exampleSurveillanceTimeByRowtbl, echo = FALSE}
pander(head(data.frame(em1[2],em1[ncol(em1)]),5))

```

<br>

```{r exampleSurveilLastTest, echo = TRUE}
em3<-SurveilLastTest(Myendo,'HospitalNumber','Dateofprocedure')
```

<br>

```{r exampleSurveilLastTesttbl, echo = FALSE}
pander(head(data.frame(em3[2],em3[5]),5))
```

<br>

```{r exampleSurveilFirstTest, echo = TRUE}
em4<-SurveilFirstTest(Myendo,'HospitalNumber','Dateofprocedure')
```

<br>

```{r exampleSurveilFirstTesttbl, echo = FALSE}
pander(head(data.frame(em4[2],em4[5]),5))
```

<br>

Of course we may also want to know how many tests have been done over a time period and this is provided by the function **HowManyTests**

<br>

This function will return the number of tests by day, month and year so they can be easily graphed according to what you want.

<br>

```{r exampleHowManyTests, echo = TRUE}
how<-HowManyOverTime(Myendo,'Indications','Dateofprocedure','Surv')
```

<br>

```{r exampleHowManyTeststbl, echo = FALSE}
pander(head(data.frame(how),5))
```



## **2. Assessment of quality functions**

<br>

Quality is measured in a variety of ways. For endoscopy it is measured according to the adherence to a) standards for endoscopic documentation as well as b) detection of certain pathological conditions such as dysplasia (best summarised as lesion recognition)


<br>

### **a) Documentation Quality**

<br>

As regards adherence to documentation for example, a generic function is provided that will look up the presence of words presented in a list in a target column. It will then output the proportion of reports that have these words, as well as a barchart to show what proportion of the endoscopies showed these words. The list can be comprised of terms that should be mentioned in a report.

<br>

_Input_
                      
```{r exampleListLookup, eval = TRUE,echo=FALSE}
panderOptions('table.split.table', Inf)
pander(head(data.frame(Myendo[2:3],Myendo[13])))
```

<br>

In this example we are looking for the words Barrett's and coeliac as perhaps we have chosen the macroscopic recognition of these features to denote what an endoscopist should always describe in the endoscopy report:

<br>

```{r exampleListLookup2, eval = TRUE}
library(tm)
myNotableWords <- c("barrett", "coeliac")
ListLookup(Myendo,'Findings',myNotableWords)
```

<br>

```{r exampleListLookup3, echo=FALSE}
#pander::panderOptions('table.split.table', Inf)
#pander(head(tt))
```

<br>

So we can see that the terms are present in the minority of reports across endoscopists, so perhaps we can look into this further..

<br>

### **b) Endoscopic Quality**

#### **Sedation Usage**

<br>

Another measure of quality is the assessment of those factors that are recorded at endoscopy such as degree of sedation used etc. Rather than provide a function for each metric, again a generic function is provided that uses any quantifiable metric and plots it against the endoscopist. This function returns a list with two elements- the plot and the table:

<br>

```{r exampleEndoscChopperMeds,echo=TRUE}
#We have to attach the output of EndoscMeds to the original dataframe
MyendoNew<-cbind(EndoscMeds(Myendo$Medications),Myendo)
#Average Fentanyl use by endoscopist:
Mytable<-MetricByEndoscopist(MyendoNew,'Endoscopist','Fent')
```

<br>

```{r exampleEndoscChopperMedstbl,echo=FALSE}

pander(head(data.frame(MyendoNew$Endoscopist,MyendoNew$Fent),10))
```

<br>
<br>



<br>


## **4.Patient flow functions**

<br>

### **Sankey plots**

<br>

We often like to get an overview of how patients are flowing through a system overall. This can give a nice visual representation of whether which patients diverge from the normal flow through a system so we can study them further. There are two ways to look at this. Sankey plots give good timepoint by timepoint representation of flow. This really works with more than one type of event at each timepoint. 

<br>

For example, if we have a dataset with events such as 'radiofrequency ablation' and 'endoscopic mucosal resection' or 'nothing' we can use the Sankey plot to determine the order of events over a large patient population. You choose the column in the dataframe that describes the Procedure type ("EMR","RFA","nothing" in this case)

<br>

```{r exampleSurveySankey, eval = FALSE}
#how<-SurveySankey(Myendo,"ProcedurePerformed")
```

<br>



```{r fig.width=12, fig.height=8,fig.align='center',echo=FALSE,out.width = "100%"}
knitr::include_graphics("img/EndoMineR_Sankey.svg")
```

<br>


### **Circos plots**

<br>

We may need something even more aggregated. Perhaps we want to see the overall number of patients that go from one event to another regardless of which timepoint it is at. To do this we can use a circos plot, which makes use of the circlize library, as follows:

```{r examplePatientFlow_CircosPlots, eval = FALSE}
#flow<-PatientFlow_CircosPlots(v,"Date.y","pHospitalNum","ProcedurePerformed")
```

<br>


```{r fig.width=12, fig.height=8,fig.align='center',echo=FALSE,out.width = "60%"}
knitr::include_graphics("img/EndoMineR_Circos.svg")
```

<br>---
title: "Polyps"
author: "Sebastian Zeki"
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    toc_depth: 5
    css: style.css
vignette: >
  %\VignetteIndexEntry{Polyps}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
library(pander)
library(EndoMineR)
knitr::opts_chunk$set(echo = TRUE)

```

```{r global_options, include=FALSE}
knitr::opts_chunk$set( warning=FALSE, message=FALSE)
```


##**Pathology Detection Quality**

<br>

###**Polyp and sample location**

<br>

A difficult area is the assessment of endoscopic quality by looking at the pathology processed from an endoscopy. This package is excellent at dealing with this kind of question because of its ability to merge the datasets together:

A particularly well developed area to look at is that of the [Global Rating Scale](http://www.healthcareimprovementscotland.org/our_work/governance_and_assurance/endoscopy_accreditation/grs_reports.aspx) for assessing the quality of colonoscopy. One of the metrics- the adenoma detection rate assesses the number of colonoscopies where at least one adenoma was detected.

One function is provided to produce a table that gives the number of adenomas, adenocarcinomas and hyperplastic polyps (also as a ration to adenomas) by endoscopist therefore immediately fulfilling the GRS requirement for the ADR as well as providing further metrics alongside.

Having made sure the text is prepared (using the **textPrep** function) and merged with Pathology reports (which have also been textPrep'd) we can go straight ahead and use the function **GRS_Type_Assess_By_Unit** on one of our example datasets called vColon (the format of the dataset can be seen in the Data vignette. Here we will use the columns "ProcedurePerformed","Endoscopist", "Diagnosis" and  "Histology" to calculate the adenoma detection rate, as well as the rate of detection of adenocarcinoma, high grade dysplasia, low grade dysplasia, serrated adenomas and hyperplastic adenomas.


```{r exampleGRS}
  data(vColon)
  nn<-GRS_Type_Assess_By_Unit(
  vColon, "ProcedurePerformed",
  "Endoscopist", "Diagnosis", "Histology")
```

```{r exampleGRStbl,echo=FALSE}
   pander(nn)
```

<br>
 
This of course can then be graphed using the pre-formatted function **EndoBasicGraph** 
```{r exampleGRSGraph,warning=FALSE,message=FALSE}
EndoBasicGraph(nn, "Endoscopist", "Adenocarcinoma")
```
---
title: "Barrett's Oesophagus"
author: "Sebastian Zeki"
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    css: style.css
vignette: >
  %\VignetteIndexEntry{Barretts}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



```{r setup, include=FALSE}

library(knitr)
library(EndoMineR)
library(pander)
library(prettydoc)
```

```{r global_options, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
```

<br>

## **Specific diseases - Barrett's oesophagus**

<br>

One particular disease that lends itself well to analytics, particularly as it is part of a surveillance programme, is the premalignant oesophageal condition Barrett's oesophagus. This is characterised by the growth of cells (called columnar lined epithelium) in the oesophagus. These cells usually occupy the lower part of the oesophagus as a continuous sheet from the top of the stomach to varying lengths up the oesophagus. 

This condition requires endoscopic surveillance and the timing of this depends on the prior endoscopic features (namely the length of the Barretts segment as measured by the Prague score- explained below) and the pathological stage at that endoscopy (which for non-dysplastic samples, since the revised 2013 UK guidelines, means the presence or absence of intestinal metaplasia). This can be seen in the image below (from Fitzgerald RC, et al. Gut 2013;0:1–36. doi:10.1136/gutjnl-2013-305372)

<br>

```{r ,fig.align='center',echo=FALSE,fig.width=18, fig.height=12}
knitr::include_graphics("img/BarrettsGuide.png")
```

<br>

## **1. Pre-processing Barrett's samples**

<br>

Such a dataset needs some processing prior to the analysis so for this we can turn to a specific set of function for Barrett's oesophagus itself. 

<br>

### **a) Prague score**

<br>

Firstly we need to extract the length of the Barrett’s segment. This is known as the [Prague score](https://services.nhslothian.scot/endoscopyunit/InformationForClinicalStaff/Documents/PRAGUE%20CRITERIA.pdf) and is made up of the length from the top of the gastric folds (just below the gastro-oesophageal junction) to the top of the circumferential extent of the Barrett's segment (C). In addition the maximal extent is from the top of the gastric folds to the top of the tongues of Barrett's segment (M). This gives an overall score such as C1M2.

<br>

After filtering for endoscopic indication (eg “Surveillance-Barrett’s”- this is stored in the 'Indication' column in our data set) the aim of the following function is to extract a C and M stage (Prague score) for Barrett’s samples. This is done using a regular expression where C and M stages are explicitly mentioned in the free text. Specifically it extracts the Prague score. This is usually mentoned in the 'Findings' column in our dataset but obviously the user can define which column should be searched.

```{r exampleBarretts_PragueScore,echo = TRUE,message=FALSE, warning=FALSE}
v<-Barretts_PragueScore(Myendo,'Findings','OGDReportWhole')
```

```{r exampleBarretts_PragueScore2,echo = FALSE,message=FALSE, warning=FALSE}
panderOptions('table.split.table', Inf)
pander(v[23:27,(ncol(v)-4):ncol(v)])
```

### **b) Worst pathological stage**

<br>

We also need to extract the worst pathological stage for a sample, and if non-dysplastic, determine whether the sample has intestinal metaplasia or not. This is done using 'degredation' so that it will look for the worst overall grade in the histology specimen and if not found it will look for the next worst and so on. 

It looks per report not per biopsy (it is more common for histopathology reports to contain the worst overall grade rather than individual biopsy grades).


```{r exampleBarretts_PathStage, echo = TRUE,message=FALSE, warning=FALSE}
#The histology column is the one we are interested in:
Mypath$b <- Barretts_PathStage(Mypath, "Histology")
```

```{r exampleBarretts_PathStage2, echo = TRUE,message=FALSE, warning=FALSE,echo=FALSE}
panderOptions('table.split.table', Inf)
pander(Mypath[2:3,(ncol(Mypath)-4):ncol(Mypath)])
```



### **c)Follow-up groups**

<br>

Having done these pre-processing steps, the follow-up group to which the last endoscopy belongs (rather than the patient as their biopsy results or Barrett's segment length and therefore their follow-up timing, may fluctuate over time) can be determined. 

The follow-up timing, as explained in the the original guideline flowchart above, depends on the length of the Barrett's segment and the presence of intestinal metaplasia (a type of columnar lined epithelium). If abnormal cells (dysplasia) are present the there is a different follow-up regime which we won't concern ourselves with at the moment. 

The timing of follow-up is done with the function **Barretts_FUType.** This relies on the previous functions called **Barretts_PathStage**  and **Barretts_PragueScore** having been run. The Barretts_FUType function will tell you which follow up Rule the patient should be on so that the timing of the next endoscopy can be determined. As these functions usually go together a wrapper function called **BarrettAll** is also provided.

<br>

```{r exampleBarretts_FUType, echo = TRUE,message=FALSE, warning=FALSE}
#Create the merged dataset
v<-Endomerge2(Myendo,"Dateofprocedure","HospitalNumber",Mypath,"Dateofprocedure","HospitalNumber")
#Find the worst pathological grade for that endoscopy
v$IMorNoIM <- Barretts_PathStage(v, "Histology")
#Find the Prague score for that endoscopy
b1<-Barretts_PragueScore(v, "Findings", "OGDReportWhole")
#Get the follow-up type for that endoscopy
b1$FU_Type<-Barretts_FUType(b1,"CStage","MStage","IMorNoIM")

```


```{r exampleBarretts_FUType2, echo = FALSE,message=FALSE, warning=FALSE,echo=FALSE}
panderOptions('table.split.table', Inf)
pander(b1[23:27,(ncol(b1)-4):ncol(b1)])
```





## **2.Quality assessment in Barrett's surveillance**

<br>

Many of the aspects of generic quality monitoring apply to Barrett's as well. For example, you may want to make sure that the endoscopy reports adhere to guidance about what should be in an endoscopy report. This can be assessed using the generic function **ListLookUp**




### **Quality of perfomance of Barrett's surveillance endoscopies as just by tissue sampling**

Some things are specific to Barrett's oesophagus however. One of the essential requirements to demonstrate adequate sampling of Barrett's oesophagus during endoscopy is that the endoscopist should adhere to the 'Seattle protocol' for biopsies which is to take  4 equally spaced biopsies at 2cm intervals in the circumferential part of the oesophagus. Because the macroscopic description of the pathological specimen tells us how many samples are taken overall (and rarely how many at each level but this is usually not the case for a variety of reasons) we can determine the shortfall in the number of biopsies taken, per endoscopist. Again pre-processing the Barrett's samples is pre-requisite. The Number of biopsies and their size should also be extracted using the histopathology functions.

```{r exampleBarrettsQuality_AnalysisBiopsyNumber,echo = FALSE,message=FALSE, warning=FALSE}
 # The number of average number of biopsies is then calculated and
 # compared to the average Prague C score so that those who are taking
 # too few biopsies can be determined
#Lets just use the Surveillance endoscopies:
b1<-b1[grepl("[Ss]urv",b1$Indications),]
b1$NumBx<-HistolNumbOfBx(b1$Macroscopicdescription,'specimen')
b1$BxSize <- HistolBxSize(b1$Macroscopicdescription)
b2<-BarrettsBxQual(b1,'Date.x','HospitalNumber',
                                    'Endoscopist')
```
 
 
```{r exampleBarrettsQuality_AnalysisBiopsyNumbertbl, echo = FALSE,message=FALSE, warning=FALSE,echo=FALSE}

#panderOptions('table.split.table', Inf)
panderOptions('table.split.table', Inf)
pander(b2)
```
 <br>
 
 This function will again return a dataframe with the number of biopsies taken that is outside of the number that should have been taken for a certain Prague score.
---
title: "Data"
author: Sebastian Zeki
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    toc_depth: 5
    css: style.css
vignette: >
  %\VignetteIndexEntry{Data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```{r setup, include=FALSE}
library(pander)
library(EndoMineR)
knitr::opts_chunk$set(echo = TRUE)
```

```{r global_options, include=FALSE}
knitr::opts_chunk$set( warning=FALSE, message=FALSE)
```

##Overview

It is envisaged that different users will start at different points of their data preparation. THis section is intended to explain the fake data I have created so the type of data used for the examples can be better understood.

There are several data files used. These are detailed below

## Gastroscopy

####  Raw datasets:

#####TheOGDReportFinal 
A dataset containing fake upper GI endoscopy reports. The report field is provided as a whole report without any fields having been already extracted. There are 2000 rows

#####PathDataFrameFinal 
A dataset containing fake upper GI pathology reports. The report field is provided as a whole report without any fields having been already extracted. There are 2000 rows

####Pre-extracted datasets:

#####Myendo 
This has been extracted using the Extractor method as follows from the raw text within Mypath:

```{r MyendoExtract, eval = FALSE}
mywords <- c("OGDReportWhole","HospitalNumber","PatientName",
             "GeneralPractitioner","Dateofprocedure","Endoscopist",
             "Secondendoscopist","Medications","Instrument","ExtentofExam",
             "Indications","ProcedurePerformed","Findings")
Extractor(TheOGDReportFinal,"OGDReportWhole",mywords)
```

####Mypath 

This has been extracted using the Extractor method as follows from the raw text within Mypath:

```{r MypathExtract, eval = FALSE}
mywords<-c("HospitalNumber","PatientName","DOB","GeneralPractitioner",
           "Dateofprocedure","ClinicalDetails","Macroscopicdescription",
           "Histology","Diagnosis")
Extractor(PathDataFrameFinal,"PathReportWhole",mywords)
```

The original dataset has also been added as "PathReportWhole",

##Colonoscopy

###Raw datasets

####ColonFinal

A dataset containing fake lower GI endoscopy reports. The report field is provided as a whole report without any fields having been already extracted. There are 2000 rows

####PathDataFrameFinalColon

A dataset containing fake lower GI pathology reports. The report field is provided as a whole report without any fields having been already extracted. There are 2000 rows.
---
title: "EndoMineR"
author: "Sebastian Zeki `r Sys.Date()`"
output:
  rmarkdown::html_document:
    theme: lumen
    toc: true
    tod_depth: 4
    css: style.css

vignette: >
  %\VignetteIndexEntry{EndoMineR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
library(pander)
library(EndoMineR)
library(here)
knitr::opts_chunk$set(echo = TRUE)
```

```{r global_options, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
```
## **Aims of EndoMineR**

The goal of EndoMineR is to extract as much information as possible from free or semi-structured endoscopy reports and their associated pathology specimens. 

Gastroenterology now has many standards against which practice is measured although many reporting systems do not include the reporting capability to give anything more than basic analysis. Much of the data is locked in semi-structured text. However the nature of semi-structured text means that data can be extracted in a standardised way- it just requires more manipulation. 

This package provides that manipulation so that complex endoscopic-pathological analyses, in line with recognised standards for these analyses, can be done.






##**How is the package divided?**

<br>

```{r fig.width=8, fig.height=8,fig.align='center',out.width = "100%",echo=FALSE}
knitr::include_graphics("img/EndoMineRBasic.jpg")
```

<br>

The package is basically divided into  parts. How all the parts are connected in shown in the figure above. 
The import of the raw data is left up to the user with the overall aim being that all the data is present in one dataframe. 


## **1. Extraction and cleaning**
If using raw reports, then once you have imported the data into your R environment, you can go ahead and use the **textPrep** function. If you have pre-segregated data, then use the **EndoPaste** function (explained at the end of this section). 

Once the data is ready you can use the **textPrep** function. This function is really a wrapper for a number of different functions that work on your data. 

#### **i. Dusting the data** 
Firstly the data is cleaned up so that extra newlines, spaces, unnecessary punctuation etc is dealt with. Although you dont need to know the details of this, if you want to look under the hood you can see that this is part of the **ColumnCleanUp** function.

#### **ii. Spell checking** 
The textPrep also implements a spell checker (the **spellCheck** function). This really checks the spelling of gastroenterology terms that are present in the in-built lexicons. So, if for example the report contains the terms 'Radiafrequency ablashion' then this will be corrected using Radiofrequency ablation. This is one function that acts to standardise the text in the report so that downstream analyses can be robust. 

#### **iii.Negative phrase removal** 
The text is then passed along to the **NegativeRemove** function. This will remove any phrases that contain negative sentences indicating a non-positive finding. For example it would remove 'There is no evidence of dysplasia here'. This makes text extraction analyses, which often report the detection of a positive finding, much more accurate. If you wish to include negative phrases however, you can. This function is part of the parameters you can switch on and off for the **textPrep** function.

#### **iv. Term mapping** 
The next step is to perform term mapping. This means that variations of a term are all mapped to a single common term. For example 'RFA' and 'HALO' may be both mapped to 'Radiofrequency ablation'. This is performed using the using the **DictionaryInPlaceReplace** function. This function looks through all of the lexicons included in this package to perform this. The lexicons are all manually created and consist of key-value pairs which therefore map a key which is the term variant) to a value (which is the standardised term that should be used).

#### **Sv. egregating the data** 
  Finally the text is ready to separate into columns. The basic **Extractor** function will take your data with the list of terms you have supplied that act as the column boundaries and separate your data accordingly. It is up to the user to define the list of boundary keywords. This list is made up of the words that will be used to split the document. It is very common for individual departments in both gastroenterology and pathology to use semi-structured reporting so that the same headers are used between patient reports. The Extractor then does the splitting for each pair of words, dumps the text between the delimiter in a new column and names that column with the first keyword in the pair with some cleaning up and then returns the new dataframe.



### **What if my data is already in columns?**

Often the raw data is derived from a basic query to some endoscopy software. The output is likely to be a spreadsheet with some of the data (eg Endoscopist) already placed in separate columns. The report still needs to be tidied up in order to maximise what we can get from the data and especially the free text. In order to do the term mapping and negative removal etc it is therefore necessary for EndoMineR to merge all of the columns together, for each endoscopy, but it keeps the column headers as delimiters for use later. This is the case with the top left box in the first diagram.

**EndoPaste** is provided as an optional function to get your data in the right shape to allow EndoMineR to process the data properly in this situation

**EndoPaste** outputs two things: the data, as well as a list of delimiters (basically your column headings). The delimiters can then be used in the **textPrep** function. This would work as follows:

```{r exampleEndoPaste,echo=TRUE}

#An example dataset
testList<-structure(list(PatientName = c("Tom Hardy", "Elma Fudd", "Bingo Man"), 
                         HospitalNumber = c("H55435", "Y3425345", "Z343424"), Text = c("All bad. Not good", 
"Serious issues", "from a land far away")), class = "data.frame", row.names = c(NA, -3L))
myReadyDataset<-EndoPaste(testList)


```


The dataframe can be obtained from myReadyDataset[1] and the delimiters from myReadyDataset[2]


<br>
<br>
<br>

### **Can I have an example please?**

As an example, here we use an sample dataset (which has not had separate columns selected already, so no need to use EndoPaste) as the input. This contains a synthetic dataset provided as part of the EndoMineR package.

```{r exampleExtractor,echo=TRUE}
PathDataFrameFinalColon2<-PathDataFrameFinalColon
names(PathDataFrameFinalColon2)<-"PathReportWhole"
pander(head(PathDataFrameFinalColon2,1))
```

```{r exampleExtractortbl,echo=FALSE}
pander(head(PathDataFrameFinalColon2,1))
```
<br>

We can then define the list of delimiters that will split this text into separate columns, title the columns according to the delimiters and return a dataframe. each column simply contains the text between the delimiters that the user has defined. These columns are then ready for the more refined cleaning provided by subesquent functions.

<br>


```{r exampleExtractor2,echo=TRUE}
library(EndoMineR)
mywords<-c("Hospital Number","Patient Name:","DOB:","General Practitioner:",
"Date received:","Clinical Details:","Macroscopic description:",
"Histology:","Diagnosis:")

PathDataFrameFinalColon3<-Extractor(PathDataFrameFinalColon2$PathReportWhole,mywords)

```

```{r exampleExtractor3,echo=FALSE}

panderOptions('table.split.table', Inf)
pander(head(PathDataFrameFinalColon3[,1:9],1))
```

The **Extractor** function is embedded within **textPrep** so you may never have to use it directly, but you will always have to submit a list of delimiters so that **textPrep** can use the **Extractor** to do its segregation:

```{r exampletextPrep,echo=FALSE}
#Submit delimiters
mywords<-c("Hospital Number","Patient Name:","DOB:","General Practitioner:",
"Date received:","Clinical Details:","Macroscopic description:",
"Histology:","Diagnosis:")
CleanResults<-textPrep(PathDataFrameFinal$PathReportWhole,mywords)
```

**textPrep** takes the optional parameters NegEx which tells the **textPrep** to remove negative phrases like 'There is no pathology here' from the text.





There will always be a certain amount of data cleaning that only the end user can do before data can be extracted. There is some cleaning that is common to many endoscopy reports and so functions have been provided for this. An abvious function is the cleaning of the endoscopist's name. This is done with the function **EndoscEndoscopist**. This removes titles and tidies up the names so that there aren't duplicate names (eg "Dr. Sebastian Zeki" and "Sebastian Zeki"). This is applied to any endoscopy column where the Endoscopist name has been isolated into its own column.

```{r exampleEndoscEndoscopist2,echo = TRUE}
EndoscEndoscopist(Myendo$Endoscopist[2:6])
```



## **2. Data linkage**

Endoscopy data may be linked with other types of data. The most common associated dataset is pathology data from biopsies etc taken at endoscopy. This pathology data should be processed in exactly the same way as the endoscopy data- namely with **textPrep** (with or without **EndoPaste**). 

The resulting pathology dataset should then be merged with the endoscopy dataset using **Endomerge2** which will merge all rows with the same hospital number and do a fuzzy match (up to 7 days after the endoscopy) with pathology so the right endoscopy event is associated with the right pathology report.An example of merging the included datasets Mypath and Myendo is given here:

```{r exampleEndomerge2,echo=TRUE}
v<-Endomerge2(Myendo,'Dateofprocedure','HospitalNumber',Mypath,'Dateofprocedure','HospitalNumber')
```

## **3. Deriving new data from what you have**

Once the text has been separated in to the columns of your choosing you are ready to see what we can extract. Functions are provided to allow quite complex extractions at both a general level and also for specific diseases.




#### **i) Medication**

The extraction of medication type and dose is important for lots of analyses in endoscopy. This is provided with the function **EndoscMeds**. The function currently extracts Fentanyl, Pethidine, Midazolam and Propofol doses into a separate column and reformats them as numeric columns so further calculations can be done. It outputs a dataframe so you will need to re-bind this output with the original dataframe for further analyses. This is shown below:


```{r exampleEndoCleaningFuncMed, echo = TRUE}
MyendoMeds<-cbind(EndoscMeds(Myendo$Medications), Myendo)
```

```{r exampleEndoCleaningFuncMedtbl, echo = FALSE}
pander(head(MyendoMeds[1:4],5))

```
#### **ii) Endosccopy Event extraction**

The **EndoscopyEvent** function will extract any event that has been performed at the endoscopy and dump it in a new column. It does this by looking in pairs of sentences and therefore wraps around a more basic function called **EndoscopyPairs_TwoSentence** which you dont have to directly use. This allows us to get the site that the event (usually a therapy) happened on. A example sentence might be 'There was a bleeding gastric ulcer. A clip was applied to it' We can only extract stomach:clip by reference to the two sentences.


#### **iii) Histology biopsy number extraction**

There is also a lot of information we can extract into a separate column from the histopathology information.
We can derive the number of biopsies taken (usually specified in the macroscopic description of a sample) using the function **HistolNumbOfBx**. 

In order to extract the numbers, the limit of what has to be extracted has to be set as part of the regex so that the function takes whatever word limits the selection.It collects everything from the regex [0-9]{1,2}.{0,3} to whatever the string boundary is.

For example, a report might say:

<br>

```{r exampleHistolNumbOfBx1, echo = TRUE}
sg<-data.frame(Mypath$HospitalNumber,Mypath$PatientName,Mypath$Macroscopicdescription)
```

```{r exampleHistolNumbOfBx1tbl, echo = FALSE}
pander(head(sg,5))
```
<br>

Based on this, the word that limits the number you are interested in is 'specimen' so the function and it's output is:

<br>

```{r exampleHistolNumbOfBx2, echo = TRUE}
Mypath$NumbOfBx<-HistolNumbOfBx(Mypath$Macroscopicdescription,'specimen')
sh<-data.frame(Mypath$HospitalNumber,Mypath$PatientName,Mypath$NumbOfBx)
```

```{r exampleHistolNumbOfBx2tbl, echo = FALSE}
pander(head(sh,5))
```

#### **iv) Histology biopsy size extraction**

We might also be interested in the size of the biopsy taken. A further function called **HistolBxSize** is provided for this. This is also derived from the macroscopic description of the specimen

```{r exampleHistolBxSize1, echo = TRUE}
Mypath$BxSize<-HistolBxSize(Mypath$Macroscopicdescription)
sh<-data.frame(Mypath$HospitalNumber,Mypath$PatientName,Mypath$BxSize)
```

```{r exampleHistolBxSize1tbl, echo = FALSE}

pander(head(sh,5))
```
#### **v) Histology type and site extraction**

A final function is also provided called **HistolTypeAndSite**. This is particularly useful when trying to find out which biopsies came from which site. 
The output is provided as a site:specimen type pair. An alternative output is also provided which groups locations (eg the gastro-oesophageal junction is seen as part of the oesophagus). 

##**And this is just the beginning...**

What is described above are the building blocks for starting a more complex analysis of the Endoscopic and Pathological data. 

Generic analyses, such as figuring out surveillance intervals, can be determined using the appropriate functions in the Analysis module.

The dataset can also be fed in to more complex functions as are described in the Barrett's and Polyp modules

The package also provides generic data visualisation tools to assess Patient flow, amonst other visualisations. All of these can be found in the associated vignettes.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoMineR.R
\name{SurveilFirstTest}
\alias{SurveilFirstTest}
\title{Extracts the first test only per patient}
\usage{
SurveilFirstTest(dataframe, HospNum_Id, Endo_ResultPerformed)
}
\arguments{
\item{dataframe}{dataframe}

\item{HospNum_Id}{Patient ID}

\item{Endo_ResultPerformed}{Date of the Endoscopy}
}
\description{
Extracts the first test only per patient and returns a new dataframe listing the
patientID and the first test done
}
\examples{
dd <- SurveilFirstTest(
  Myendo, "HospitalNumber",
  "Dateofprocedure"
)
}
\seealso{
Other Basic Analysis - Surveillance Functions: 
\code{\link{HowManyOverTime}()},
\code{\link{SurveilLastTest}()},
\code{\link{SurveilTimeByRow}()},
\code{\link{TimeToStatus}()}
}
\concept{Basic Analysis - Surveillance Functions}
\keyword{Surveillance}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graphics.R
\name{scale_fill_Publication}
\alias{scale_fill_Publication}
\title{Set the fills for all the ggplots}
\usage{
scale_fill_Publication()
}
\description{
This standardises the fills for any ggplot plot produced.
If you do use it, like all ggplots it can be extended using the 
"+" to add whatever else is necessary
}
\examples{
# None needed
}
\seealso{
Other Data Presentation helpers: 
\code{\link{EndoBasicGraph}()},
\code{\link{scale_colour_Publication}()},
\code{\link{theme_Publication}()}
}
\concept{Data Presentation helpers}
\keyword{ggplot}
\keyword{themes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{NegativeRemoveWrapper}
\alias{NegativeRemoveWrapper}
\title{Wrapper for Negative Remove}
\usage{
NegativeRemoveWrapper(inputText)
}
\arguments{
\item{inputText}{the text to remove Negatives from}
}
\value{
This returns a column within a dataframe. This should be changed to a 
character vector eventually
}
\description{
This performs negative removal on a per sentance basis
}
\examples{
# Build a character vector and then
# incorporate into a dataframe
anexample<-c("There is no evidence of polyp here",
"Although the prep was poor,there was no adenoma found",
"The colon was basically inflammed, but no polyp was seen",
"The Barrett's segment was not biopsied",
"The C0M7 stretch of Barrett's was flat")
anexample<-data.frame(anexample)
names(anexample)<-"Thecol"
# Run the function on the dataframe and it should get rid of sentences (and
# parts of sentences) with negative parts in them.
#hh<-NegativeRemoveWrapper(anexample$Thecol)
}
\seealso{
Other NLP - Text Cleaning and Extraction: 
\code{\link{ColumnCleanUp}()},
\code{\link{DictionaryInPlaceReplace}()},
\code{\link{Extractor}()},
\code{\link{NegativeRemove}()},
\code{\link{textPrep}()}
}
\concept{NLP - Text Cleaning and Extraction}
\keyword{Negative}
\keyword{Sentences}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Module_Polyps.R
\name{GRS_Type_Assess_By_Unit}
\alias{GRS_Type_Assess_By_Unit}
\title{Create GRS metrics by endoscopist (X-ref with pathology)}
\usage{
GRS_Type_Assess_By_Unit(dataframe, ProcPerformed, Endo_Endoscopist, Dx, Histol)
}
\arguments{
\item{dataframe}{The dataframe}

\item{ProcPerformed}{The column containing the Procedure type performed}

\item{Endo_Endoscopist}{column containing the Endoscopist name}

\item{Dx}{The column with the Histological diagnosis}

\item{Histol}{The column with the Histology text in it}
}
\description{
This extracts the polyps types from the data
(for colonoscopy and flexible sigmoidosscopy data)
and outputs the adenoma,adenocarcinoma and
hyperplastic detection rate by endoscopist as well
as overall number of colonoscopies.
This will be extended to other GRS outputs in the future.
}
\examples{
nn <- GRS_Type_Assess_By_Unit(
  vColon, "ProcedurePerformed",
  "Endoscopist", "Diagnosis", "Original.y"
)
}
\concept{Disease Specific Analysis - Polyp functions}
\keyword{Withdrawal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{vColon}
\alias{vColon}
\title{Fake Lower GI Endoscopy Set including Pathology}
\format{
A data frame with 2000 rows and 26 variables:
\describe{
  \item{pHospitalNum}{The HospitalNum, in text}
  \item{PatientName.x}{The PatientName, in text}
  \item{GeneralPractitioner.x}{The GeneralPractitioner report, in text}
  \item{Date.x}{The Date, in date}
  \item{Endoscopist}{The Endoscopist report, in text}
  \item{Secondendoscopist}{The Secondendoscopist report, in text}
  \item{Medications}{The Medications report, in text}
  \item{Instrument}{The Instrument report, in text}
  \item{ExtentofExam}{The ExtentofExam report, in text}
  \item{Indications}{The Indications report, in text}
  \item{ProcedurePerformed}{The ProcedurePerformed report, in text}
  \item{Findings}{The Findings report, in text}
  \item{EndoscopicDiagnosis}{The EndoscopicDiagnosis report, in text}
  \item{Original.x}{The Original endosocpy report, in text}
  \item{eHospitalNum}{The HospitalNum, in text}
  \item{PatientName.y}{The PatientName, in text}
  \item{DOB}{The DOB, in date}
  \item{GeneralPractitioner.y}{The GeneralPractitioner report, in text}
  \item{Date.y}{The Date.y , in date}
  \item{ClinicalDetails}{The ClinicalDetails report, in text}
  \item{Natureofspecimen}{The Natureofspecimen report, in text}
  \item{Macroscopicdescription}{The Macroscopicdescription report, in text}
  \item{Histology}{The Histology report, in text}
  \item{Diagnosis}{The Diagnosis report, in text}
  \item{Original.y}{The whole report, in text}
  \item{Days}{Days, in numbers}
}
}
\usage{
vColon
}
\description{
A dataset containing fake lower GI endoscopy reports and pathology reports all
pre-extracted
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graphics.R
\name{SurveySankey}
\alias{SurveySankey}
\title{Create a Sankey plot for patient flow}
\usage{
SurveySankey(dfw, ProcPerformedColumn, PatientID)
}
\arguments{
\item{dfw}{the dataframe extracted using the standard cleanup scripts}

\item{ProcPerformedColumn}{the column containing the test like P
rocPerformed for example}

\item{PatientID}{the column containing the patients unique identifier
eg hostpital number}
}
\description{
The purpose of the function is to provide a Sankey plot 
which allows the analyst to see the proportion
of patients moving from one state (in this case type of Procedure) to
another. This allows us to see for example how many EMRs are done after
RFA.
}
\examples{
names(Myendo)[names(Myendo) == "HospitalNumber"] <- "PatientID"
gg <- SurveySankey(Myendo, "ProcedurePerformed", "PatientID")

}
\seealso{
Other Patient Flow functions: 
\code{\link{PatientFlowIndividual}()}
}
\concept{Patient Flow functions}
\keyword{Sankey}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_EndoPathExtractors.R
\name{EndoscEndoscopist}
\alias{EndoscEndoscopist}
\title{Clean endoscopist column}
\usage{
EndoscEndoscopist(EndoscopistColumn)
}
\arguments{
\item{EndoscopistColumn}{The endoscopy text column}
}
\value{
This returns a character vector
}
\description{
If an endoscopist column is part of the dataset once the extractor
function has been used this cleans the endoscopist column from the report.
It gets rid of titles
It gets rid of common entries that are not needed.
It should be used after the textPrep function
}
\examples{
Myendo$Endoscopist <- EndoscEndoscopist(Myendo$Endoscopist)
}
\seealso{
Other Endoscopy specific cleaning functions: 
\code{\link{EndoscInstrument}()},
\code{\link{EndoscMeds}()},
\code{\link{EndoscopyEvent}()}
}
\concept{Endoscopy specific cleaning functions}
\keyword{Endoscopist}
\keyword{extraction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_EndoPathExtractors.R
\name{HistolTypeAndSite}
\alias{HistolTypeAndSite}
\title{Extract the site a specimen was removed from as well as the type}
\usage{
HistolTypeAndSite(inputString1, inputString2, procedureString)
}
\arguments{
\item{inputString1}{The first column to look in}

\item{inputString2}{The second column to look in}

\item{procedureString}{The column with the procedure in it}
}
\value{
a list with two columns, one is the type and site and the other
is the index to be used for OPCS4 coding later if needed.
}
\description{
This needs some blurb to be written. Used in the OPCS4 coding
}
\examples{
Myendo2<-Endomerge2(Myendo,'Dateofprocedure','HospitalNumber',
Mypath,'Dateofprocedure','HospitalNumber')
PathSiteAndType <- HistolTypeAndSite(Myendo2$PathReportWhole,
Myendo2$Macroscopicdescription, Myendo2$ProcedurePerformed)
}
\seealso{
Other Histology specific cleaning functions: 
\code{\link{HistolBxSize}()},
\code{\link{HistolNumbOfBx}()}
}
\concept{Histology specific cleaning functions}
\keyword{Find}
\keyword{and}
\keyword{replace}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PathDataFrameFinalColon}
\alias{PathDataFrameFinalColon}
\title{Fake Lower GI Pathology Set}
\format{
A data frame with 2000 rows and 1 variables:
\describe{
  \item{PathReportWhole}{The whole report, in text}
}
}
\usage{
PathDataFrameFinalColon
}
\description{
A dataset containing fake pathology reports for lower GI endoscopy tissue specimens.
The report field is provided as a whole report
without any fields having been already extracted
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{LocationListUpper}
\alias{LocationListUpper}
\title{Use list of standard locations for upper GI endoscopy}
\usage{
LocationListUpper()
}
\description{
The is a list of standard locations at endoscopy that is used in the
extraction of the site of biopsies/EMRs and potentially in functions looking at the site of a 
therapeutic event. It just returns the list in the function.
}
\seealso{
Other NLP - Lexicons: 
\code{\link{BiopsyIndex}()},
\code{\link{EventList}()},
\code{\link{GISymptomsList}()},
\code{\link{HistolType}()},
\code{\link{LocationListLower}()},
\code{\link{LocationListUniversal}()},
\code{\link{LocationList}()},
\code{\link{RFACath}()},
\code{\link{WordsToNumbers}()}
}
\concept{NLP - Lexicons}
\keyword{Location}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Module_Barretts.R
\name{Barretts_PragueScore}
\alias{Barretts_PragueScore}
\title{Extract the Prague score}
\usage{
Barretts_PragueScore(dataframe, EndoReportColumn, EndoReportColumn2)
}
\arguments{
\item{dataframe}{dataframe with column of interest}

\item{EndoReportColumn}{column of interest}

\item{EndoReportColumn2}{second column of interest}
}
\description{
The aim is to extract a C and M stage (Prague score) for Barrett's samples.
This is done using a regex where C and M stages are explicitly mentioned in
the free text
Specfically it extracts the Prague score
}
\examples{
# The example takes the endoscopy demo dataset and searches the
# Findings column (which contains endoscopy free text about the
# procedure itself). It then extracts the Prague score if relevant. I
# find it easiest to use this on a Barrett's subset of data rather than
# a dump of all endoscopies but of course this is a permissible dataset
# too


aa <- Barretts_PragueScore(Myendo, "Findings", "OGDReportWhole")
}
\seealso{
Other Disease Specific Analysis - Barretts Data: 
\code{\link{BarrettsAll}()},
\code{\link{BarrettsBxQual}()},
\code{\link{BarrettsParisEMR}()},
\code{\link{Barretts_FUType}()},
\code{\link{Barretts_PathStage}()}
}
\concept{Disease Specific Analysis - Barretts Data}
\keyword{Prague}
\keyword{score}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{spellCheck}
\alias{spellCheck}
\title{Find and Replace}
\usage{
spellCheck(pattern, replacement, x, fixed = FALSE)
}
\arguments{
\item{pattern}{the pattern to look for}

\item{replacement}{the pattern replaceme with}

\item{x}{the target string}

\item{fixed}{whether the pattern is regex or not. Default not.}
}
\value{
This returns a character vector
}
\description{
This is a helper function for finding and replacing from lexicons
like the event list. The lexicons are all named lists where the name
is the text to replace and the value what it should be replaced with
It uses fuzzy find and replace to account for spelling errors
}
\examples{
L <- tolower(stringr::str_split(HistolType(),"\\\\|"))
}
\concept{NLP - Text Cleaning and Extraction
inputText<-TheOGDReportFinal$OGDReportWhole
inputText<-Reduce(function(x, nm) spellCheck(nm, L[[nm]], x), init = inputText, names(L))}
\keyword{Find}
\keyword{and}
\keyword{replace}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Module_Eosinophil.R
\name{Eosinophilics}
\alias{Eosinophilics}
\title{Extract the Prague score}
\usage{
Eosinophilics(dataframe, findings, histol, IndicationsFroExamination)
}
\arguments{
\item{dataframe}{dataframe with column of interest}

\item{findings}{column of interest}

\item{histol}{second column of interest}

\item{IndicationsFroExamination}{second column of interest}
}
\description{
The aim is to extract a C and M stage (Prague score) for Barrett's samples.
This is done using a regex where C and M stages are explicitly mentioned in
the free text
Specfically it extracts the Prague score
}
\examples{
# Firstly relevant columns are extrapolated from the
# Mypath demo dataset. These functions are all part of Histology data
# cleaning as part of the package.
# Mypath demo dataset. These functions are all part of Histology data
# cleaning as part of the package.
v <- Mypath
v$NumBx <- HistolNumbOfBx(v$Macroscopicdescription, "specimen")
v$BxSize <- HistolBxSize(v$Macroscopicdescription)
# The histology is then merged with the Endoscopy dataset. The merge occurs
# according to date and Hospital number
v <- Endomerge2(
  Myendo, "Dateofprocedure", "HospitalNumber", v, "Dateofprocedure",
  "HospitalNumber"
)


aa <- Eosinophilics(v, "Findings", "Histology","Indications")
}
\concept{Disease Specific Analysis - Eosinophilic Data}
\keyword{EoE}
\keyword{Eosinophilic}
\keyword{Oesophagitis}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_EndoPathExtractors.R
\name{dev_ExtrapolateOPCS4Prep}
\alias{dev_ExtrapolateOPCS4Prep}
\title{OPCS-4 Coding}
\usage{
dev_ExtrapolateOPCS4Prep(dataframe, Procedure, PathSite, Event, extentofexam)
}
\arguments{
\item{dataframe}{the dataframe}

\item{Procedure}{The Procedure column}

\item{PathSite}{The column containing the Pathology site}

\item{Event}{the EVENT column}

\item{extentofexam}{the furthest point reached in the examination}
}
\description{
This function extracts the OPCS-4 codes for all Barrett's procedures
It should take the OPCS-4 from the EVENT and perhaps also using extent
depending on how the coding is done. The EVENT column will need to 
extract multiple findings
The hope is that the OPCS-4 column will then map from the EVENT column. This returns a nested list 
column with the procedure, furthest path site and event performed
}
\examples{
# Need to run the HistolTypeSite and EndoscopyEvent functions first here
# SelfOGD_Dunn$OPCS4w<-ExtrapolateOPCS4Prep(SelfOGD_Dunn,"PROCEDUREPERFORMED",
# "PathSite","EndoscopyEvent")
}
\keyword{OPCS-4}
\keyword{codes}
\keyword{extraction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{EventList}
\alias{EventList}
\title{Use list of endoscopic events and procedures}
\usage{
EventList()
}
\description{
This function returns all the conversions from common version of events to
a standardised event list, much like the Location standardisation function
This does not include EMR as this is
extracted from the pathology so is part of pathology type.
}
\examples{
# unique(unlist(EventList(), use.names = FALSE))
}
\seealso{
Other NLP - Lexicons: 
\code{\link{BiopsyIndex}()},
\code{\link{GISymptomsList}()},
\code{\link{HistolType}()},
\code{\link{LocationListLower}()},
\code{\link{LocationListUniversal}()},
\code{\link{LocationListUpper}()},
\code{\link{LocationList}()},
\code{\link{RFACath}()},
\code{\link{WordsToNumbers}()}
}
\concept{NLP - Lexicons}
\keyword{Event}
\keyword{extraction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{EntityPairs_OneSentence}
\alias{EntityPairs_OneSentence}
\title{See if words from two lists co-exist within a sentence}
\usage{
EntityPairs_OneSentence(inputText, list1, list2)
}
\arguments{
\item{inputText}{The relevant pathology text column}

\item{list1}{First list to refer to}

\item{list2}{The second list to look for}
}
\description{
See if words from two lists co-exist within a sentence. Eg site and tissue type.
This function only looks in one sentence for the two terms. If you suspect the terms may
occur in adjacent sentences then use the EntityPairs_TwoSentence function.
}
\examples{
# tbb<-EntityPairs_OneSentence(Mypath$Histology,HistolType(),LocationList())
}
\seealso{
Other Basic Column mutators: 
\code{\link{EntityPairs_TwoSentence}()},
\code{\link{ExtrapolatefromDictionary}()},
\code{\link{ListLookup}()},
\code{\link{MyImgLibrary}()}
}
\concept{Basic Column mutators}
\keyword{PathPairLookup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{EndoPaste}
\alias{EndoPaste}
\title{Paste endoscopy and histology results into one}
\usage{
EndoPaste(x)
}
\arguments{
\item{x}{the dataframe}
}
\value{
This returns a list with a dataframe containing one column of the merged text
and a character vector which is the delimiter list for when the textPrep function is used
}
\description{
As spreadsheets are likely to be submitted with pre-segregated data as appears from 
endoscopy software output, these should be remerged prior to cleaning. This function
takes the column headers and places it before each text so that the original
full text is recreated. It will use the column headers as the delimiter. This should 
be used before textPrep as the textPrep function takes a character vector (ie the whole
report and not a segregated one) only
}
\examples{
testList<-structure(list(PatientName = c("Tom Hardy", "Elma Fudd", "Bingo Man"
), HospitalNumber = c("H55435", "Y3425345", "Z343424"), Text = c("All bad. Not good", 
"Serious issues", "from a land far away")), class = "data.frame", row.names = c(NA, -3L))
EndoPaste(testList)
}
\concept{NLP - Text merging:}
\keyword{Merge}
\keyword{columns}
\keyword{dataframe}
\keyword{into}
\keyword{one}
\keyword{text}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{MyImgLibrary}
\alias{MyImgLibrary}
\title{Clean html endoscopic images}
\usage{
MyImgLibrary(file, delim, location)
}
\arguments{
\item{file}{The html report to extract (the html will have all the images references in it)}

\item{delim}{The phrase that separates individual endoscopies}

\item{location}{The folder containing the actual images}
}
\description{
This is used to pick and clean endoscopic images from html exports so they can be prepared
before being linked to pathology and endoscopy reports
}
\examples{
# MyImgLibrary("~/Images Captured with Proc Data Audit_Findings1.html",
#                         "procedureperformed","~/")
}
\seealso{
Other Basic Column mutators: 
\code{\link{EntityPairs_OneSentence}()},
\code{\link{EntityPairs_TwoSentence}()},
\code{\link{ExtrapolatefromDictionary}()},
\code{\link{ListLookup}()}
}
\concept{Basic Column mutators}
\keyword{Image}
\keyword{extraction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoFakeData.R
\name{Endomerge2}
\alias{Endomerge2}
\title{Merge endoscopy and histology data.}
\usage{
Endomerge2(x, EndoDate, EndoHospNumber, y, PathDate, PathHospNumber)
}
\arguments{
\item{x}{Endoscopy dataframe}

\item{EndoDate}{The date the endoscopy was performed}

\item{EndoHospNumber}{The unique hospital number in the endoscopy dataset}

\item{y}{Histopathology dataframe}

\item{PathDate}{The date the endoscopy was performed}

\item{PathHospNumber}{The unique hospital number in the endoscopy dataset}
}
\description{
This takes the endoscopy dataset date
performed and the hospital number column
and merges with the equivalent column in the pathology dataset. This is
merged within a 7 day time frame as pathology is often reported after
endoscopic
}
\examples{
v <- Endomerge2(
  Myendo, "Dateofprocedure", "HospitalNumber",
  Mypath, "Dateofprocedure", "HospitalNumber"
)
}
\keyword{and}
\keyword{endoscopy}
\keyword{histology}
\keyword{merge}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoMineR1.R
\docType{package}
\name{EndoMineR}
\alias{EndoMineR}
\title{EndoMineR: A package for analysis of endoscopic and related pathology}
\description{
The goal of EndoMineR is to extract as much information as possible from
endoscopy reports and their associated pathology specimens. The package
is intended for use by gastroenterologists, pathologists and anyone
interested in the analysis of endoscopic and ppathological datasets
Gastroenterology now has many standards against which practice is measured
although many reporting systems do not include the reporting capability to
give anything more than basic analysis. Much of the data is locked in
semi-structured text.However the nature of semi-structured text means that
data can be extracted in a standardised way- it just requires more
manipulation. This package provides that manipulation so that complex
endoscopic-pathological analyses, in line with recognised standards for
these analyses, can be done.The package is basically in three parts/
}
\details{
\itemize{
\item The extraction- This is really when the data is provided as full text
reports. You may already have the data in a spreadsheet in which case
this part isn't necessary.

\item Cleaning- These are a group of functions
that allow the user to extract and clean data commonly found in endoscopic
and pathology reports. The cleaning functions usually remove common typos or
extraneous information and do some reformatting.

\item  Analyses- The analyses provide graphing function as well as analyses
according to the cornerstone questions in gastroenterology- namely
surveillance, patient tracking, quality of endoscopy and pathology
reporting and diagnostic yield questions.
}

To learn more about EndoMineR, start with the vignettes:
`browseVignettes(package = "EndoMineR")`
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{LocationListUniversal}
\alias{LocationListUniversal}
\title{Use list of standard locations for upper GI endoscopy}
\usage{
LocationListUniversal()
}
\description{
The is a list of standard locations at endoscopy that is used in the
extraction of the site of biopsies/EMRs and potentially in functions looking at the site of a
therapeutic event. It just returns the list in the function
}
\seealso{
Other NLP - Lexicons: 
\code{\link{BiopsyIndex}()},
\code{\link{EventList}()},
\code{\link{GISymptomsList}()},
\code{\link{HistolType}()},
\code{\link{LocationListLower}()},
\code{\link{LocationListUpper}()},
\code{\link{LocationList}()},
\code{\link{RFACath}()},
\code{\link{WordsToNumbers}()}
}
\concept{NLP - Lexicons}
\keyword{Location}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Module_Barretts.R
\name{BarrettsParisEMR}
\alias{BarrettsParisEMR}
\title{Run the Paris classification versus worst histopath grade for Barrett's}
\usage{
BarrettsParisEMR(Column, Column2)
}
\arguments{
\item{Column}{Endoscopy report field of interest as a string vector}

\item{Column2}{Another endoscopy report field of interest as a string vector}
}
\value{
a string vector
}
\description{
This creates a column of Paris grade for all samples where this is mentioned.
}
\examples{
# 
Myendo$EMR<-BarrettsParisEMR(Myendo$ProcedurePerformed,Myendo$Findings)
}
\seealso{
Other Disease Specific Analysis - Barretts Data: 
\code{\link{BarrettsAll}()},
\code{\link{BarrettsBxQual}()},
\code{\link{Barretts_FUType}()},
\code{\link{Barretts_PathStage}()},
\code{\link{Barretts_PragueScore}()}
}
\concept{Disease Specific Analysis - Barretts Data}
\keyword{Does}
\keyword{data}
\keyword{something}
\keyword{with}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Module_Barretts.R
\name{BarrettsBxQual}
\alias{BarrettsBxQual}
\title{Get the number of Barrett's biopsies taken}
\usage{
BarrettsBxQual(dataframe, Endo_ResultPerformed, PatientID, Endoscopist)
}
\arguments{
\item{dataframe}{dataframe}

\item{Endo_ResultPerformed}{Date of the Endoscopy}

\item{PatientID}{Patient's unique identifier}

\item{Endoscopist}{name of the column with the Endoscopist names}
}
\description{
This function gets the number of biopsies taken per 
endoscopy and compares it to the
Prague score for that endoscopy.Endoscopists should be taking a certain
number of biopsies given the length of a Barrett's segment so it
should be straightforward to detect a shortfall in the number
of biopsies being taken. The output is the shortfall per endoscopist
}
\examples{
# Firstly relevant columns are extrapolated from the
# Mypath demo dataset. These functions are all part of Histology data
# cleaning as part of the package.
Mypath$NumBx <- HistolNumbOfBx(Mypath$Macroscopicdescription, "specimen")
Mypath$BxSize <- HistolBxSize(Mypath$Macroscopicdescription)

# The histology is then merged with the Endoscopy dataset. The merge occurs
# according to date and Hospital number
v <- Endomerge2(
  Myendo, "Dateofprocedure", "HospitalNumber", Mypath, "Dateofprocedure",
  "HospitalNumber"
)

# The function relies on the other Barrett's functions being run as well:
b1 <- Barretts_PragueScore(v, "Findings")
b1$PathStage <- Barretts_PathStage(b1, "Histology")

# The follow-up group depends on the histology and the Prague score for a
# patient so it takes the processed Barrett's data and then looks in the
# Findings column for permutations of the Prague score.
b1$FU_Type <- Barretts_FUType(b1, "CStage", "MStage", "PathStage")


colnames(b1)[colnames(b1) == "pHospitalNum"] <- "HospitalNumber"
# The number of average number of biopsies is then calculated and
# compared to the average Prague C score so that those who are taking
# too few biopsies can be determined
hh <- BarrettsBxQual(
  b1, "Date.x", "HospitalNumber",
  "Endoscopist"
)
rm(v)
}
\seealso{
Other Disease Specific Analysis - Barretts Data: 
\code{\link{BarrettsAll}()},
\code{\link{BarrettsParisEMR}()},
\code{\link{Barretts_FUType}()},
\code{\link{Barretts_PathStage}()},
\code{\link{Barretts_PragueScore}()}
}
\concept{Disease Specific Analysis - Barretts Data}
\keyword{Does}
\keyword{data}
\keyword{something}
\keyword{with}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ColonFinal}
\alias{ColonFinal}
\title{Fake Lower GI Endoscopy Set}
\format{
A data frame with 2000 rows and 1 variables:
\describe{
  \item{OGDReportWhole}{The whole report, in text}
}
}
\usage{
ColonFinal
}
\description{
A dataset containing fake lower GI endoscopy reports. The report field is provided as a whole report
without any fields having been already extracted
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{WordsToNumbers}
\alias{WordsToNumbers}
\title{Convetr words to numbers especially for the histopathology text}
\usage{
WordsToNumbers()
}
\description{
This function converts words to numbers.
}
\seealso{
Other NLP - Lexicons: 
\code{\link{BiopsyIndex}()},
\code{\link{EventList}()},
\code{\link{GISymptomsList}()},
\code{\link{HistolType}()},
\code{\link{LocationListLower}()},
\code{\link{LocationListUniversal}()},
\code{\link{LocationListUpper}()},
\code{\link{LocationList}()},
\code{\link{RFACath}()}
}
\concept{NLP - Lexicons}
\keyword{Event}
\keyword{extraction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_EndoPathExtractors.R
\name{HistolNumbOfBx}
\alias{HistolNumbOfBx}
\title{Extract the number of biopsies taken from the histology report}
\usage{
HistolNumbOfBx(inputString, regString)
}
\arguments{
\item{inputString}{The input text to process}

\item{regString}{The keyword to remove and to stop at in the regex}
}
\description{
This extracts the number of biopsies taken from the pathology report.
This is usually from the Macroscopic description column.
It collects everything from the regex [0-9]{1,2}.{0,3}
to whatever the string boundary is (z).
}
\examples{
qq <- HistolNumbOfBx(Mypath$Macroscopicdescription, "specimen")
}
\seealso{
Other Histology specific cleaning functions: 
\code{\link{HistolBxSize}()},
\code{\link{HistolTypeAndSite}()}
}
\concept{Histology specific cleaning functions}
\keyword{Biopsy}
\keyword{number}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graphics.R
\name{PatientFlowIndividual}
\alias{PatientFlowIndividual}
\title{Create a plot over time of patient categorical findings as a line chart}
\usage{
PatientFlowIndividual(
  theframe,
  EndoReportColumn,
  myNotableWords,
  DateofProcedure,
  PatientID
)
}
\arguments{
\item{theframe}{dataframe}

\item{EndoReportColumn}{the column containing the date of the procedure}

\item{myNotableWords}{The terms from a column with categorical variables}

\item{DateofProcedure}{Column with the date of the procedure}

\item{PatientID}{Column with the patient's unique identifier}
}
\description{
This plots the findings at endoscopy (or pathology) over time for individual
patients. An example might be with worst pathological grade on biopsy for 
Barrett's oesophagus over time
}
\examples{
# This function builds chart of categorical outcomes for individal patients over time
# It allows a two dimensional visualisation of patient progress. A perfect example is 
# visualising the Barrett's progression for patients on surveillance and then
# therapy if dysplasia develops and highlighting recurrence if it happens
# Barretts_df <- BarrettsAll(Myendo, "Findings", "OGDReportWhole", Mypath, "Histology")
# myNotableWords<-c("No_IM","IM","LGD","HGD","T1a","IGD","SM1","SM2")
# PatientFlowIndividual(Barretts_df,"IMorNoIM",myNotableWords,DateofProcedure,"HospitalNumber")
# Once the function is run you should always call dev.off()
}
\seealso{
Other Patient Flow functions: 
\code{\link{SurveySankey}()}
}
\concept{Patient Flow functions}
\keyword{Patient}
\keyword{flow}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graphics.R
\name{scale_colour_Publication}
\alias{scale_colour_Publication}
\title{Set the colour theme for all the ggplots}
\usage{
scale_colour_Publication()
}
\description{
This standardises the colours for any ggplot plot produced.
If you do use it, like all ggplots it can be extended using the 
"+" to add whatever else is necessary
}
\examples{
# None needed
}
\seealso{
Other Data Presentation helpers: 
\code{\link{EndoBasicGraph}()},
\code{\link{scale_fill_Publication}()},
\code{\link{theme_Publication}()}
}
\concept{Data Presentation helpers}
\keyword{ggplot}
\keyword{themes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graphics.R
\name{theme_Publication}
\alias{theme_Publication}
\title{Set the publication theme for all the ggplots}
\usage{
theme_Publication(base_size = 14, base_family = "Helvetica")
}
\arguments{
\item{base_size}{the base size}

\item{base_family}{the base family}
}
\description{
This standardises the theme for any ggplot plot produced.
If you do use it, like all ggplots it can be extended using the 
"+" to add whatever else is necessary
}
\examples{
# None needed
}
\seealso{
Other Data Presentation helpers: 
\code{\link{EndoBasicGraph}()},
\code{\link{scale_colour_Publication}()},
\code{\link{scale_fill_Publication}()}
}
\concept{Data Presentation helpers}
\keyword{ggplot}
\keyword{themes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{EntityPairs_TwoSentence}
\alias{EntityPairs_TwoSentence}
\title{Look for relationships between site and event}
\usage{
EntityPairs_TwoSentence(inputString, list1, list2)
}
\arguments{
\item{inputString}{The relevant pathology text column}

\item{list1}{The intial list to assess}

\item{list2}{The other list to look for}
}
\description{
This is used to look for relationships between site and event especially for endoscopy events
where sentences such as 'The stomach polyp was large. It was removed with a snare' ie the therapy
and the site are in two different locations.
}
\examples{
# tbb<-EntityPairs_TwoSentence(Myendo$Findings,EventList(),HistolType())
}
\seealso{
Other Basic Column mutators: 
\code{\link{EntityPairs_OneSentence}()},
\code{\link{ExtrapolatefromDictionary}()},
\code{\link{ListLookup}()},
\code{\link{MyImgLibrary}()}
}
\concept{Basic Column mutators}
\keyword{Find}
\keyword{and}
\keyword{replace}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Module_Barretts.R
\name{BarrettsAll}
\alias{BarrettsAll}
\title{Run all the basic Barrett's functions}
\usage{
BarrettsAll(
  Endodataframe,
  EndoReportColumn,
  EndoReportColumn2,
  Pathdataframe,
  PathColumn
)
}
\arguments{
\item{Endodataframe}{endoscopy dataframe of interest}

\item{EndoReportColumn}{Endoscopy report field of interest as a string vector}

\item{EndoReportColumn2}{Second endoscopy report field of interest as a string vector}

\item{Pathdataframe}{pathology dataframe of interest}

\item{PathColumn}{Pathology report field of interest as a string vector}
}
\value{
Newdf
}
\description{
Function to encapsulate all the Barrett's functions together. This includes the Prague
score and the worst pathological grade and then feeds both of these things into
the follow up function. The output is a dataframe with all the original data as
well as the new columns that have been created.
}
\examples{
Barretts_df <- BarrettsAll(Myendo, "Findings", "OGDReportWhole", Mypath, "Histology")
}
\seealso{
Other Disease Specific Analysis - Barretts Data: 
\code{\link{BarrettsBxQual}()},
\code{\link{BarrettsParisEMR}()},
\code{\link{Barretts_FUType}()},
\code{\link{Barretts_PathStage}()},
\code{\link{Barretts_PragueScore}()}
}
\concept{Disease Specific Analysis - Barretts Data}
\keyword{Does}
\keyword{data}
\keyword{something}
\keyword{with}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{LocationListLower}
\alias{LocationListLower}
\title{Use list of standard locations for lower GI endoscopy}
\usage{
LocationListLower()
}
\description{
The is a list of standard locations at endoscopy that is used in the
extraction of the site of biopsies/EMRs and potentially in functions looking at the site of a
therapeutic event. It just returns the list in the function
}
\seealso{
Other NLP - Lexicons: 
\code{\link{BiopsyIndex}()},
\code{\link{EventList}()},
\code{\link{GISymptomsList}()},
\code{\link{HistolType}()},
\code{\link{LocationListUniversal}()},
\code{\link{LocationListUpper}()},
\code{\link{LocationList}()},
\code{\link{RFACath}()},
\code{\link{WordsToNumbers}()}
}
\concept{NLP - Lexicons}
\keyword{Location}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{Extractor}
\alias{Extractor}
\title{Extract columns from the raw text}
\usage{
Extractor(inputString, delim)
}
\arguments{
\item{inputString}{the column to extract from}

\item{delim}{the vector of words that will be used as the boundaries to
extract against}
}
\description{
This is the main extractor for the Endoscopy and Histology report.
This relies on the user creating a list of words representing the
subheadings. The list is then fed to the
Extractor so that it acts as the beginning and the end of the
regex used to split the text. Whatever has been specified in the list
is used as a column header. Column headers don't tolerate special characters
like : or ? and / and don't allow numbers as the start character so these
have to be dealt with in the text before processing
}
\examples{
# As column names cant start with a number, one of the dividing
# words has to be converted
# A list of dividing words (which will also act as column names)
# is then constructed
mywords<-c("Hospital Number","Patient Name:","DOB:","General Practitioner:",
"Date received:","Clinical Details:","Macroscopic description:",
"Histology:","Diagnosis:")
Mypath2<-Extractor(PathDataFrameFinal$PathReportWhole,mywords)
}
\seealso{
Other NLP - Text Cleaning and Extraction: 
\code{\link{ColumnCleanUp}()},
\code{\link{DictionaryInPlaceReplace}()},
\code{\link{NegativeRemoveWrapper}()},
\code{\link{NegativeRemove}()},
\code{\link{textPrep}()}
}
\concept{NLP - Text Cleaning and Extraction}
\keyword{Extraction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{GISymptomsList}
\alias{GISymptomsList}
\title{Index of GI symptoms}
\usage{
GISymptomsList()
}
\description{
This function returns all the common GI symptoms. 
They are simply listed as is without grouping
or mapping. They have been derived from a manual 
list with synonyms derived from the UMLS Methatharus
using the browser.
}
\seealso{
Other NLP - Lexicons: 
\code{\link{BiopsyIndex}()},
\code{\link{EventList}()},
\code{\link{HistolType}()},
\code{\link{LocationListLower}()},
\code{\link{LocationListUniversal}()},
\code{\link{LocationListUpper}()},
\code{\link{LocationList}()},
\code{\link{RFACath}()},
\code{\link{WordsToNumbers}()}
}
\concept{NLP - Lexicons}
\keyword{Event}
\keyword{extraction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{RFACath}
\alias{RFACath}
\title{Use list of catheters used in radiofrequency ablation}
\usage{
RFACath()
}
\description{
The takes a list of catheters used in radiofrequency ablation.
}
\seealso{
Other NLP - Lexicons: 
\code{\link{BiopsyIndex}()},
\code{\link{EventList}()},
\code{\link{GISymptomsList}()},
\code{\link{HistolType}()},
\code{\link{LocationListLower}()},
\code{\link{LocationListUniversal}()},
\code{\link{LocationListUpper}()},
\code{\link{LocationList}()},
\code{\link{WordsToNumbers}()}
}
\concept{NLP - Lexicons}
\keyword{RFA}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{BiopsyIndex}
\alias{BiopsyIndex}
\title{Index biopsy locations}
\usage{
BiopsyIndex()
}
\description{
This function returns all the conversions from common version of events to
a standardised event list, much like the Location standardidastion function
This does not include EMR as this is 
extracted from the pathology so is part of pathology type. It is used for
automated OPCS-4 coding.
}
\seealso{
Other NLP - Lexicons: 
\code{\link{EventList}()},
\code{\link{GISymptomsList}()},
\code{\link{HistolType}()},
\code{\link{LocationListLower}()},
\code{\link{LocationListUniversal}()},
\code{\link{LocationListUpper}()},
\code{\link{LocationList}()},
\code{\link{RFACath}()},
\code{\link{WordsToNumbers}()}
}
\concept{NLP - Lexicons}
\keyword{Event}
\keyword{extraction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoMineR.R
\name{CategoricalByEndoscopist}
\alias{CategoricalByEndoscopist}
\title{Group anything by Endoscopist and returns the table}
\usage{
CategoricalByEndoscopist(ProportionColumn, EndoscopistColumn)
}
\arguments{
\item{ProportionColumn}{The column (categorical data) of interest}

\item{EndoscopistColumn}{The endoscopist column}
}
\description{
This creates a proportion table for categorical variables by endoscopist
It of course relies on a Endoscopist column being present
}
\examples{
# The function plots any numeric metric by endoscopist
# Mypath demo dataset. These functions are all part of Histology data
# cleaning as part of the package.
v <- Mypath
v$NumBx <- HistolNumbOfBx(Mypath$Macroscopicdescription, "specimen")
v$BxSize <- HistolBxSize(v$Macroscopicdescription)
# The histology is then merged with the Endoscopy dataset. The merge occurs
# according to date and Hospital number
v <- Endomerge2(
  Myendo, "Dateofprocedure", "HospitalNumber", v, "Dateofprocedure",
  "HospitalNumber"
)
# The function relies on the other Barrett's functions being run as well:
v$IMorNoIM <- Barretts_PathStage(v, "Histology")
colnames(v)[colnames(v) == "pHospitalNum"] <- "HospitalNumber"
# The function takes the column with the extracted worst grade of
# histopathology and returns the proportion of each finding (ie
# proportion with low grade dysplasia, high grade etc.) for each
# endoscopist
kk <- CategoricalByEndoscopist(v$IMorNoIM, v$Endoscopist)
rm(Myendo)
}
\seealso{
Other Grouping by endoscopist: 
\code{\link{MetricByEndoscopist}()}
}
\concept{Grouping by endoscopist}
\keyword{Endoscopist}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graphics.R
\name{PatientFlow_CircosPlots}
\alias{PatientFlow_CircosPlots}
\title{Create a Circos plot for patient flow}
\usage{
PatientFlow_CircosPlots(
  dataframe,
  Endo_ResultPerformed,
  HospNum_Id,
  ProcPerformed
)
}
\arguments{
\item{dataframe}{dataframe}

\item{Endo_ResultPerformed}{the column containing the date of the procedure}

\item{HospNum_Id}{Column with the patient's unique hospital number}

\item{ProcPerformed}{The procedure that you want to plot (eg EMR,
radiofrequency ablation for Barrett's but can be
any dscription of a procedure you desire)}
}
\description{
This allows us to look at the overall flow from one
type of procedure to another using circos plots. A good example of it's 
use might be to see how patients move from one state (e.g. having an
EMR), to another state (e.g. undergoing RFA)
}
\examples{
# This function builds a circos plot which gives a more aggregated
# overview of how patients flow from one state to another than the
# SurveySankey function
# Build a list of procedures
Event <- list(
  x1 = "Therapeutic- Dilatation",
  x2 = "Other-", x3 = "Surveillance",
  x4 = "APC", x5 = "Therapeutic- RFA TTS",
  x5 = "Therapeutic- RFA 90",
  x6 = "Therapeutic- EMR", x7 = "Therapeutic- RFA 360"
)
EndoEvent <- replicate(2000, sample(Event, 1, replace = FALSE))
# Merge the list with the Myendo dataframe
fff <- unlist(EndoEvent)
fff <- data.frame(fff)
names(fff) <- "col1"
Myendo$EndoEvent<-fff$col1
names(Myendo)[names(Myendo) == "HospitalNumber"] <- "PatientID"
names(Myendo)[names(Myendo) == "fff$col1"] <- "EndoEvent"
# Myendo$EndoEvent<-as.character(Myendo$EndoEvent)
# Run the function using the procedure information (the date of the
# procedure, the Event type and the individual patient IDs)
hh <- PatientFlow_CircosPlots(Myendo, "Dateofprocedure", "PatientID", "EndoEvent")
rm(Myendo)
rm(EndoEvent)
}
\keyword{@family}
\keyword{Circos}
\keyword{Flow}
\keyword{Patient}
\keyword{functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{LocationList}
\alias{LocationList}
\title{Use list of upper and lower GI standard locations}
\usage{
LocationList()
}
\description{
The is a list of standard locations at endoscopy. It used for the site of biopsies/EMRs 
and potentially in functions looking at the site of a 
therapeutic event. It just returns the list in the function.
}
\seealso{
Other NLP - Lexicons: 
\code{\link{BiopsyIndex}()},
\code{\link{EventList}()},
\code{\link{GISymptomsList}()},
\code{\link{HistolType}()},
\code{\link{LocationListLower}()},
\code{\link{LocationListUniversal}()},
\code{\link{LocationListUpper}()},
\code{\link{RFACath}()},
\code{\link{WordsToNumbers}()}
}
\concept{NLP - Lexicons}
\keyword{Location}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_EndoPathExtractors.R
\name{EndoscopyEvent}
\alias{EndoscopyEvent}
\title{Extract the endoscopic event.}
\usage{
EndoscopyEvent(dataframe, EventColumn1, Procedure, Macroscopic, Histology)
}
\arguments{
\item{dataframe}{datafrane of interest}

\item{EventColumn1}{The relevant endoscopt free text column describing the findings}

\item{Procedure}{Column saying which procedure was performed}

\item{Macroscopic}{Column describing all the macroscopic specimens}

\item{Histology}{Column with free text histology (usually microscopic histology)}
}
\value{
This returns a character vector
}
\description{
This extracts the endoscopic event. It looks for the event 
term and then looks in the event sentence as well as the one above to see if
the location is listed. It only looks within the endoscopy fields. If tissue is taken
then this will be extracted with the HistolTypeAndSite function rather than being 
listed as a result as this is cleaner and more robust.
}
\examples{
# Myendo$EndoscopyEvent<-EndoscopyEvent(Myendo,"Findings",
# "ProcedurePerformed","MACROSCOPICALDESCRIPTION","HISTOLOGY")
}
\seealso{
Other Endoscopy specific cleaning functions: 
\code{\link{EndoscEndoscopist}()},
\code{\link{EndoscInstrument}()},
\code{\link{EndoscMeds}()}
}
\concept{Endoscopy specific cleaning functions}
\keyword{Find}
\keyword{and}
\keyword{replace}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoMineR.R
\name{HowManyOverTime}
\alias{HowManyOverTime}
\title{Number of tests done per month and year by indication}
\usage{
HowManyOverTime(dataframe, Indication, Endo_ResultPerformed, StringToSearch)
}
\arguments{
\item{dataframe}{dataframe}

\item{Indication}{Indication column}

\item{Endo_ResultPerformed}{column containing date the Endoscopy was
performed}

\item{StringToSearch}{The string in the Indication to search for}
}
\description{
Get an overall idea of how many endoscopies have been done for an indication
by year and month. This is a more involved version of
SurveilCapacity function. It takes string for
the Indication for the test
}
\details{
This returns a list which contains a plot (number of tests for that
indication over time and a table with the same information broken down
by month and year).
}
\examples{
# This takes the dataframe MyEndo (part of the package examples) and looks in
# the column which holds the test indication (in this example it is called
# 'Indication' The date of the procedure column(which can be date format or
# POSIX format) is also necessary.  Finally the string which indicates the text
# indication needs to be inpoutted. In this case we are looking for all endoscopies done
# where the indication is surveillance (so searching on 'Surv' will do fine).
# If you want all the tests then put '.*' instead of Surv
rm(list = ls(all = TRUE))
ff <- HowManyOverTime(Myendo, "Indications", "Dateofprocedure", ".*")
}
\seealso{
Other Basic Analysis - Surveillance Functions: 
\code{\link{SurveilFirstTest}()},
\code{\link{SurveilLastTest}()},
\code{\link{SurveilTimeByRow}()},
\code{\link{TimeToStatus}()}
}
\concept{Basic Analysis - Surveillance Functions}
\keyword{Tests}
\keyword{number}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{DictionaryInPlaceReplace}
\alias{DictionaryInPlaceReplace}
\title{Dictionary In Place Replace}
\usage{
DictionaryInPlaceReplace(inputString, list)
}
\arguments{
\item{inputString}{the input string (ie the full medical report)}

\item{list}{The replacing list}
}
\value{
This returns a character vector
}
\description{
This maps terms in the text and replaces them with the 
standardised term (mapped in the lexicon file) within the text.
It is used within the textPrep function.
}
\examples{
inputText<-DictionaryInPlaceReplace(TheOGDReportFinal$OGDReportWhole,LocationList())
}
\seealso{
Other NLP - Text Cleaning and Extraction: 
\code{\link{ColumnCleanUp}()},
\code{\link{Extractor}()},
\code{\link{NegativeRemoveWrapper}()},
\code{\link{NegativeRemove}()},
\code{\link{textPrep}()}
}
\concept{NLP - Text Cleaning and Extraction}
\keyword{Replace}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoMineR.R
\name{TimeToStatus}
\alias{TimeToStatus}
\title{Extract the time to an event}
\usage{
TimeToStatus(dataframe, HospNum, EVENT, indicatorEvent, endEvent)
}
\arguments{
\item{dataframe}{The dataframe}

\item{HospNum}{The Hospital Number column}

\item{EVENT}{The column that contains the outcome of choice}

\item{indicatorEvent}{The name of the start event (can be a regular expression)}

\item{endEvent}{The name of the endpoint (can be a regular expression)}
}
\description{
This function selects patients who have had a start event and an end
event of the users choosing so you can determine things like how long
it takes to get a certain outcome. For example, how long does it take to
get a patient into a fully squamous oesophagus after Barrett's ablation
for dysplasia?
}
\examples{
# Firstly relevant columns are extrapolated from the
# Mypath demo dataset. These functions are all part of Histology data
# cleaning as part of the package.
v <- Mypath
v$NumBx <- HistolNumbOfBx(v$Macroscopicdescription, "specimen")
v$BxSize <- HistolBxSize(v$Macroscopicdescription)

# The histology is then merged with the Endoscopy dataset. The merge occurs
# according to date and Hospital number
v <- Endomerge2(
  Myendo, "Dateofprocedure", "HospitalNumber", v, "Dateofprocedure",
  "HospitalNumber"
)

# The function relies on the other Barrett's functions being run as well:
b1 <- Barretts_PragueScore(v, "Findings")
b1$IMorNoIM <- Barretts_PathStage(b1, "Histology")
colnames(b1)[colnames(b1) == "pHospitalNum"] <- "HospitalNumber"

# The function groups the procedures by patient and gives
# all the procedures between
# the indicatorEvent amd the procedure just after the endpoint.
# Eg if the start is RFA and the
# endpoint is biopsies then it will give all RFA procedures and
# the first biopsy procedure

b1$EndoscopyEvent <- EndoscopyEvent(
  b1, "Findings", "ProcedurePerformed",
  "Macroscopicdescription", "Histology"
)
nn <- TimeToStatus(b1, "eHospitalNum", "EndoscopyEvent", "rfa", "dilat")
rm(v)
}
\seealso{
Other Basic Analysis - Surveillance Functions: 
\code{\link{HowManyOverTime}()},
\code{\link{SurveilFirstTest}()},
\code{\link{SurveilLastTest}()},
\code{\link{SurveilTimeByRow}()}
}
\concept{Basic Analysis - Surveillance Functions}
\keyword{ourcome}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{ExtrapolatefromDictionary}
\alias{ExtrapolatefromDictionary}
\title{Extrapolate from Dictionary}
\usage{
ExtrapolatefromDictionary(inputString, list)
}
\arguments{
\item{inputString}{The text string to process}

\item{list}{of words to iterate through}
}
\description{
Provides term mapping and extraction in one.
Standardises any term according to a mapping lexicon provided and then
extracts the term. This is
different to the DictionaryInPlaceReplace in that it provides a new column
with the extracted terms as opposed to changing it in place
}
\examples{
#Firstly we extract histology from the raw report
# The function then standardises the histology terms through a series of
# regular expressions and then extracts the type of tissue 
Mypath$Tissue<-suppressWarnings(
suppressMessages(
ExtrapolatefromDictionary(Mypath$Histology,HistolType()
)
)
)
rm(MypathExtraction)
}
\seealso{
Other Basic Column mutators: 
\code{\link{EntityPairs_OneSentence}()},
\code{\link{EntityPairs_TwoSentence}()},
\code{\link{ListLookup}()},
\code{\link{MyImgLibrary}()}
}
\concept{Basic Column mutators}
\keyword{Withdrawal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TheOGDReportFinal}
\alias{TheOGDReportFinal}
\title{Fake Upper GI Endoscopy Set}
\format{
A data frame with 2000 rows and 1 variables:
\describe{
  \item{OGDReportWhole}{The whole report, in text}
}
}
\usage{
TheOGDReportFinal
}
\description{
A dataset containing fake endoscopy reports. The report field is provided as a whole report
without any fields having been already extracted
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Module_Barretts.R
\name{Barretts_FUType}
\alias{Barretts_FUType}
\title{Determine the Follow up group}
\usage{
Barretts_FUType(dataframe, CStage, MStage, IMorNoIM)
}
\arguments{
\item{dataframe}{the dataframe(which has to have been processed by the
Barretts_PathStage function first to get IMorNoIM and the Barretts_PragueScore
to get the C and M stage if available),}

\item{CStage}{CStage column}

\item{MStage}{MStage column}

\item{IMorNoIM}{IMorNoIM column}
}
\description{
This determines the follow up rule a patient should fit in to (according to
the British Society for Gastroenterology guidance on Barrett's oesophagus)
Specfically it combines the presence of intestinal metaplasia with
Prague score so the follow-up group can be determined. It relies on the
presence of a Prague score. It should be run after
Barretts_PathStage which looks for the worst stage of a
specimen and which will determine the presence or absence of intestinal
metaplasia if the sample is non-dysplastic. Because reports often do not record
a Prague score a more pragmatic approach as been to assess the M stage and if
this is not present then to use the C stage extrapolated using the
Barretts_Prague function
}
\examples{
# Firstly relevant columns are extrapolated from the
# Mypath demo dataset. These functions are all part of Histology data
# cleaning as part of the package.
# Mypath demo dataset. These functions are all part of Histology data
# cleaning as part of the package.
v <- Mypath
v$NumBx <- HistolNumbOfBx(v$Macroscopicdescription, "specimen")
v$BxSize <- HistolBxSize(v$Macroscopicdescription)
# The histology is then merged with the Endoscopy dataset. The merge occurs
# according to date and Hospital number
v <- Endomerge2(
  Myendo, "Dateofprocedure", "HospitalNumber", v, "Dateofprocedure",
  "HospitalNumber"
)
# The function relies on the other Barrett's functions being run as well:
v$IMorNoIM <- Barretts_PathStage(v, "Histology")
v <- Barretts_PragueScore(v, "Findings")

# The follow-up group depends on the histology and the Prague score for a
# patient so it takes the processed Barrett's data and then looks in the
# Findings column for permutations of the Prague score.
v$FU_Type <- Barretts_FUType(v, "CStage", "MStage", "IMorNoIM")
rm(v)
}
\seealso{
Other Disease Specific Analysis - Barretts Data: 
\code{\link{BarrettsAll}()},
\code{\link{BarrettsBxQual}()},
\code{\link{BarrettsParisEMR}()},
\code{\link{Barretts_PathStage}()},
\code{\link{Barretts_PragueScore}()}
}
\concept{Disease Specific Analysis - Barretts Data}
\keyword{Follow-Up}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Mypath}
\alias{Mypath}
\title{Fake Pathology report}
\format{
A data frame with 2000 rows and 1 variables:
\describe{
  \item{PathReportWhole}{The whole report, in text}
  \item{HospitalNumber}{Hospital Number, in text}
  \item{PatientName}{Patient Name, in text}
  \item{DOB}{Date of Birth, in text}
  \item{GeneralPractitioner}{General Practitioner, in text}
  \item{Dateofprocedure}{Date of the procedure, as date}
  \item{ClinicalDetails}{Clinical Details, in text}
  \item{Macroscopicdescription}{Macroscopic description of the report, in text}
  \item{Histology}{Histology, in text}
  \item{Diagnosis}{Diagnosis, in text}
}
}
\usage{
Mypath
}
\description{
A dataset containing fake pathology reports.
The report field is derived from the whole report as follows:
Mypath<-PathDataFrameFinalColon
HistolTree<-list('Hospital Number','Patient Name','DOB:','General Practitioner:',
'Date of procedure:','Clinical Details:','Macroscopic description:','Histology:','Diagnosis:','')
for(i in 1:(length(HistolTree)-1)) {
 Mypath<-Extractor(Mypath,'PathReportWhole',as.character(HistolTree[i]),
 as.character(HistolTree[i+1]),as.character(HistolTree[i]))
}
Mypath$Dateofprocedure<-as.Date(Mypath$Dateofprocedure)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{ColumnCleanUp}
\alias{ColumnCleanUp}
\title{Tidy up messy columns}
\usage{
ColumnCleanUp(vector)
}
\arguments{
\item{vector}{column of interest}
}
\value{
This returns a character vector
}
\description{
This does a general clean up of whitespace,
semi-colons,full stops at the start
of lines and converts end sentence full stops to new lines.
}
\examples{
ii<-ColumnCleanUp(Myendo$Findings)
}
\seealso{
Other NLP - Text Cleaning and Extraction: 
\code{\link{DictionaryInPlaceReplace}()},
\code{\link{Extractor}()},
\code{\link{NegativeRemoveWrapper}()},
\code{\link{NegativeRemove}()},
\code{\link{textPrep}()}
}
\concept{NLP - Text Cleaning and Extraction}
\keyword{Cleaner}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReportBuilder.R
\name{sanity}
\alias{sanity}
\title{Create a basic consort diagram from dataframes}
\usage{
sanity(pathName)
}
\arguments{
\item{pathName}{The string in the Indication to search for}
}
\description{
This function creates a consort diagram using 
diagrammeR by assessing all of the dataframes in your script
and populating each box in the consort diagram with the 
number of rows in each dataframe as well as how the dataframes are linked
together. The user just provides a pathname for the script
}
\examples{
#pathName<-paste0(here::here(),"/inst/TemplateProject/munge/PreProcessing.R")
#sanity(pathName)
# This creates a consort diagram from any R script (not Rmd). It
# basically tells you how all the dataframes are related and how many
# rows each dataframe has so you can see if any data has been lost
# on the way.
}
\keyword{consort}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_EndoPathExtractors.R
\name{HistolBxSize}
\alias{HistolBxSize}
\title{Determine the largest biopsy size from the histology report}
\usage{
HistolBxSize(MacroColumn)
}
\arguments{
\item{MacroColumn}{Macdescrip}
}
\description{
This extracts the biopsy size from the report. If there are multiple
biopsies it will extract the overall size of each one (size is calculated
usually in cubic mm from the three dimensions provided). This will result
in row duplication.
}
\details{
This is usually from the Macroscopic description column.
}
\examples{
rr <- HistolBxSize(Mypath$Macroscopicdescription)
}
\seealso{
Other Histology specific cleaning functions: 
\code{\link{HistolNumbOfBx}()},
\code{\link{HistolTypeAndSite}()}
}
\concept{Histology specific cleaning functions}
\keyword{biopsy}
\keyword{size}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Myendo}
\alias{Myendo}
\title{Fake Endoscopies}
\format{
A data frame with 2000 rows and 1 variables:
\describe{
  \item{OGDReportWhole}{The whole report, in text}
  \item{HospitalNumber}{Hospital Number, in text}
  \item{PatientName}{Patient Name, in text}
  \item{GeneralPractitioner}{General Practitioner, in text}
  \item{Dateofprocedure}{Date of the procedure, as date}
  \item{Endoscopist}{Endoscopist, in text}
  \item{Secondendoscopist}{Secondendoscopist, in text}
  \item{Medications}{Medications, in text}
  \item{Instrument}{Instrument, in text}
  \item{ExtentofExam}{ExtentofExam, in text}
  \item{Indications}{Indications, in text}
  \item{ProcedurePerformed}{Procedure Performed, in text}
  \item{Findings}{Endoscopic findings, in text}
}
}
\usage{
Myendo
}
\description{
A dataset containing fake endoscopy reports. The report fields have already been
The report field is derived from the whole report as follows:
Myendo<-TheOGDReportFinal
Myendo$OGDReportWhole<-gsub('2nd Endoscopist:','Second endoscopist:',Myendo$OGDReportWhole)
EndoscTree<-list('Hospital Number:','Patient Name:','General Practitioner:',
'Date of procedure:','Endoscopist:','Second endoscopist:','Medications',
'Instrument','Extent of Exam:','Indications:','Procedure Performed:','Findings:',
'Endoscopic Diagnosis:')
for(i in 1:(length(EndoscTree)-1)) {
 Myendo<-Extractor(Myendo,'OGDReportWhole',as.character(EndoscTree[i]),
 as.character(EndoscTree[i+1]),as.character(EndoscTree[i]))
}
Myendo$Dateofprocedure<-as.Date(Myendo$Dateofprocedure)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Module_Barretts.R
\name{Barretts_PathStage}
\alias{Barretts_PathStage}
\title{Get the worst pathological stage for Barrett's}
\usage{
Barretts_PathStage(dataframe, PathColumn)
}
\arguments{
\item{dataframe}{dataframe with column of interest}

\item{PathColumn}{column of interest}
}
\description{
This extracts the pathological stage from the histopathology specimen. It is
done using 'degradation' so that it will look for the worst overall grade
in the histology specimen and if not found it will look for the next worst
and so on. It looks per report not per biopsy (it is more common
for histopathology reports to contain the worst overall grade
rather than individual biopsy grades).
Specfically it extracts the histopathology worst grade within the specimen
FOr the sake of accuracy this should alwats be used after the HistolDx function
and this removes negative sentences such as 'there is no dysplasia'.
This current function should be used on the column derived from HistolDx
which is called Dx_Simplified
}
\examples{
# Firstly relevant columns are extrapolated from the
# Mypath demo dataset. These functions are all part of Histology data
# cleaning as part of the package.
# The function then takes the Histology column from the merged data set (v).
# It extracts the worst histological grade for a specimen
b <- Barretts_PathStage(Mypath, "Histology")
rm(v)
}
\seealso{
Other Disease Specific Analysis - Barretts Data: 
\code{\link{BarrettsAll}()},
\code{\link{BarrettsBxQual}()},
\code{\link{BarrettsParisEMR}()},
\code{\link{Barretts_FUType}()},
\code{\link{Barretts_PragueScore}()}
}
\concept{Disease Specific Analysis - Barretts Data}
\keyword{Pathology}
\keyword{extraction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{textPrep}
\alias{textPrep}
\title{Combine all the text cleaning and extraction functions into one}
\usage{
textPrep(inputText, delim)
}
\arguments{
\item{inputText}{The relevant pathology text columns}

\item{delim}{the delimitors so the extractor can be used}
}
\value{
This returns a string vector.
}
\description{
This function prepares the data by cleaning 
punctuation, checking spelling against the lexicons, mapping terms
according to the lexicons and lower casing everything. 
It contains several of the other functions
in the package for ease of use.
}
\examples{
mywords<-c("Hospital Number","Patient Name:","DOB:","General Practitioner:",
"Date received:","Clinical Details:","Macroscopic description:",
"Histology:","Diagnosis:")
CleanResults<-textPrep(PathDataFrameFinal$PathReportWhole,mywords)
}
\seealso{
Other NLP - Text Cleaning and Extraction: 
\code{\link{ColumnCleanUp}()},
\code{\link{DictionaryInPlaceReplace}()},
\code{\link{Extractor}()},
\code{\link{NegativeRemoveWrapper}()},
\code{\link{NegativeRemove}()}
}
\concept{NLP - Text Cleaning and Extraction}
\keyword{cleaning}
\keyword{text}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PathDataFrameFinal}
\alias{PathDataFrameFinal}
\title{Fake Upper GI Pathology Set}
\format{
A data frame with 2000 rows and 1 variables:
\describe{
  \item{PathReportWhole}{The whole report, in text}
}
}
\usage{
PathDataFrameFinal
}
\description{
A dataset containing fake pathology reports for upper GI endoscopy tissue specimens.
The report field is provided as a whole report
without any fields having been already extracted
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoMineR.R
\name{SurveilLastTest}
\alias{SurveilLastTest}
\title{Extract the last test done by a patient only}
\usage{
SurveilLastTest(dataframe, HospNum_Id, Endo_ResultPerformed)
}
\arguments{
\item{dataframe}{dataframe}

\item{HospNum_Id}{Patient ID}

\item{Endo_ResultPerformed}{Date of the Endoscopy}
}
\description{
This extracts the last test only per patient and returns a new dataframe listing the
patientID and the last test done
}
\examples{
cc <- SurveilLastTest(Myendo, "HospitalNumber", "Dateofprocedure")
}
\seealso{
Other Basic Analysis - Surveillance Functions: 
\code{\link{HowManyOverTime}()},
\code{\link{SurveilFirstTest}()},
\code{\link{SurveilTimeByRow}()},
\code{\link{TimeToStatus}()}
}
\concept{Basic Analysis - Surveillance Functions}
\keyword{Surveillance}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoMineR.R
\name{MetricByEndoscopist}
\alias{MetricByEndoscopist}
\title{Plot a metric by endoscopist}
\usage{
MetricByEndoscopist(dataframe, Column, EndoscopistColumn)
}
\arguments{
\item{dataframe}{The dataframe}

\item{Column}{The column (numeric data) of interest}

\item{EndoscopistColumn}{The endoscopist column}
}
\description{
This takes any of the numerical metrics in the dataset and plots it by
endoscopist.
It of course relies on a Endoscopist column being present
}
\examples{
#The function gives a table with any numeric
# metric by endoscopist
# In this example we tabulate medication by
# endoscopist
# Lets bind the output of EndoscMeds to the main dataframe so we
# have a complete dataframe with all the meds extracted
MyendoNew<-cbind(EndoscMeds(Myendo$Medications),Myendo)

# Now lets look at the fentanly use per Endoscopist:
kk<-MetricByEndoscopist(MyendoNew,'Endoscopist','Fent')
#EndoBasicGraph(MyendoNew, "Endoscopist", "Fent") #run this
#if you want to see the graph
rm(Myendo)
}
\seealso{
Other Grouping by endoscopist: 
\code{\link{CategoricalByEndoscopist}()}
}
\concept{Grouping by endoscopist}
\keyword{Endoscopist}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_EndoPathExtractors.R
\name{EndoscMeds}
\alias{EndoscMeds}
\title{Clean medication column}
\usage{
EndoscMeds(MedColumn)
}
\arguments{
\item{MedColumn}{column of interest as a string vector}
}
\value{
This returns a dataframe
}
\description{
This cleans medication column from the report assuming such a column exists.
It gets rid of common entries that are not needed. It also splits the
medication into fentanyl and midazolam numeric doses for use. 
It should be used after the textPrep function.
}
\examples{
MyendoNew <- cbind(EndoscMeds(Myendo$Medications), Myendo)
}
\seealso{
Other Endoscopy specific cleaning functions: 
\code{\link{EndoscEndoscopist}()},
\code{\link{EndoscInstrument}()},
\code{\link{EndoscopyEvent}()}
}
\concept{Endoscopy specific cleaning functions}
\keyword{Endoscopy}
\keyword{medications}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoMineR.R
\name{ListLookup}
\alias{ListLookup}
\title{Extract from report, using words from a list}
\usage{
ListLookup(theframe, EndoReportColumn, myNotableWords)
}
\arguments{
\item{theframe}{the dataframe,}

\item{EndoReportColumn}{the column of interest,}

\item{myNotableWords}{list of words you are interested in}
}
\description{
The aim here is simply to
produce a document term matrix to get the frequency
of all the words, then extract the words you are
interested in with tofind then find which reports
have those words. Then find what proportion of the reports
have those terms.
}
\examples{
# The function relies on defined a list of
# words you are interested in and then choosing the column you are
# interested in looking in for these words. This can be for histopathology
# free text columns or endoscopic. In this example it is for endoscopic
# columns
myNotableWords <- c("arrett", "oeliac")
jj <- ListLookup(Myendo, "Findings", myNotableWords)
}
\seealso{
Other Basic Column mutators: 
\code{\link{EntityPairs_OneSentence}()},
\code{\link{EntityPairs_TwoSentence}()},
\code{\link{ExtrapolatefromDictionary}()},
\code{\link{MyImgLibrary}()}
}
\concept{Basic Column mutators}
\keyword{Lookup}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lexicons.R
\name{HistolType}
\alias{HistolType}
\title{Use list of pathology types}
\usage{
HistolType()
}
\description{
This standardizes terms to describe the pathology tissue type being examined
}
\seealso{
Other NLP - Lexicons: 
\code{\link{BiopsyIndex}()},
\code{\link{EventList}()},
\code{\link{GISymptomsList}()},
\code{\link{LocationListLower}()},
\code{\link{LocationListUniversal}()},
\code{\link{LocationListUpper}()},
\code{\link{LocationList}()},
\code{\link{RFACath}()},
\code{\link{WordsToNumbers}()}
}
\concept{NLP - Lexicons}
\keyword{Pathology}
\keyword{type}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graphics.R
\name{EndoBasicGraph}
\alias{EndoBasicGraph}
\title{Basic graph creation using the template specified in theme_Publication.}
\usage{
EndoBasicGraph(dataframe, xdata, number)
}
\arguments{
\item{dataframe}{dataframe}

\item{xdata}{The x column}

\item{number}{The numeric column}
}
\value{
Myplot This is the final plot

Myplot
}
\description{
This creates a basic graph using the template specified in theme_Publication.
It takes a numeric column and plots it against any non-numeric x axis in a ggplot
}
\examples{
# This function plots numeric y vs non-numeric x
# Get some numeric columns e.g. number of biopsies and size
Mypath$Size <- HistolBxSize(Mypath$Macroscopicdescription)
Mypath$NumBx <- HistolNumbOfBx(Mypath$Macroscopicdescription, "specimen")
Mypath2 <- Mypath[, c("NumBx", "Size")]
EndoBasicGraph(Mypath, "Size", "NumBx")
}
\seealso{
Other Data Presentation helpers: 
\code{\link{scale_colour_Publication}()},
\code{\link{scale_fill_Publication}()},
\code{\link{theme_Publication}()}
}
\concept{Data Presentation helpers}
\keyword{Time}
\keyword{plots}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_CleanUp.R
\name{NegativeRemove}
\alias{NegativeRemove}
\title{Remove negative and normal sentences}
\usage{
NegativeRemove(inputText)
}
\arguments{
\item{inputText}{column of interest}
}
\value{
This returns a column within a dataframe. THis should be changed to a 
character vector eventually
}
\description{
Extraction of the negative sentences so that normal findings can be
removed and not counted when searching for true diseases. eg remove
'No evidence of candidal infection' so it doesn't get included if
looking for candidal infections. It is used by default as part of
the textPrep function but can be turned off as an optional parameter
}
\examples{
# Build a character vector and then
# incorporate into a dataframe
anexample<-c("There is no evidence of polyp here",
"Although the prep was poor,there was no adenoma found",
"The colon was basically inflammed, but no polyp was seen",
"The Barrett's segment was not biopsied",
"The C0M7 stretch of Barrett's was flat")
anexample<-data.frame(anexample)
names(anexample)<-"Thecol"
# Run the function on the dataframe and it should get rid of sentences (and
# parts of sentences) with negative parts in them.
hh<-NegativeRemove(anexample$Thecol)
}
\seealso{
Other NLP - Text Cleaning and Extraction: 
\code{\link{ColumnCleanUp}()},
\code{\link{DictionaryInPlaceReplace}()},
\code{\link{Extractor}()},
\code{\link{NegativeRemoveWrapper}()},
\code{\link{textPrep}()}
}
\concept{NLP - Text Cleaning and Extraction}
\keyword{Negative}
\keyword{Sentences}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EndoMineR.R
\name{SurveilTimeByRow}
\alias{SurveilTimeByRow}
\title{Extract the time difference between each test in days}
\usage{
SurveilTimeByRow(dataframe, HospNum_Id, Endo_ResultPerformed)
}
\arguments{
\item{dataframe}{dataframe,}

\item{HospNum_Id}{Patient ID}

\item{Endo_ResultPerformed}{Date of the Endoscopy}
}
\description{
This determines the time difference between each test for a patient in days
It returns the time since the first and the last study as a new dataframe.
}
\examples{
aa <- SurveilTimeByRow(
  Myendo, "HospitalNumber",
  "Dateofprocedure"
)
}
\seealso{
Other Basic Analysis - Surveillance Functions: 
\code{\link{HowManyOverTime}()},
\code{\link{SurveilFirstTest}()},
\code{\link{SurveilLastTest}()},
\code{\link{TimeToStatus}()}
}
\concept{Basic Analysis - Surveillance Functions}
\keyword{Surveillance}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Module_IBD.R
\name{IBD_Scores}
\alias{IBD_Scores}
\title{Cleans medication column if present}
\usage{
IBD_Scores(inputColumn1)
}
\arguments{
\item{inputColumn1}{column of interest as a string vector}
}
\value{
This returns a dataframe with all the scores in it
}
\description{
This extracts all of the relevant IBD scores where present
from the medical text.
}
\examples{
 # Example to be provided
}
\keyword{IBD}
\keyword{scores}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic_EndoPathExtractors.R
\name{EndoscInstrument}
\alias{EndoscInstrument}
\title{Clean instrument column}
\usage{
EndoscInstrument(EndoInstrument)
}
\arguments{
\item{EndoInstrument}{column of interest}
}
\value{
This returns a character vector
}
\description{
This cleans the Instument column from the report assuming such a column exists
(where instrument usually refers to the endoscope number being used.)
It gets rid of common entries that are not needed.
It should be used after the textPrep function.
Note this is possibly going to be deprecated in the next version 
as the endoscope coding used here is not widely used.
}
\examples{
Myendo$Instrument <- EndoscInstrument(Myendo$Instrument)
}
\seealso{
Other Endoscopy specific cleaning functions: 
\code{\link{EndoscEndoscopist}()},
\code{\link{EndoscMeds}()},
\code{\link{EndoscopyEvent}()}
}
\concept{Endoscopy specific cleaning functions}
\keyword{Instrument}
