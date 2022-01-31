---
title: "Bioassays: a new package in R for analyzing cellular assays"

tags:
  - R
  - Multi-well plate
  - bioassays

authors:
  - name: Anwar Azad Palakkan
    orcid: 0000-0002-1693-1219
    affiliation: 1

affiliations:
  - name: Deanery of Biomedical Sciences, University of Edinburgh, UK
    index: 1

date: 11 June 2020

bibliography: paper.bib
---

# Summary

  Experimental techniques that acquire data from multi-well plates are very common in life science research; today they account for half of all assays used in drug development and toxicity testing. These experiments include colorimetric, fluorometric and luminometric assays for cytotoxicity / cell viability [@Markossian:2004], bacterial assays [@Burton:2007; @Pant:2016] and immunoassays [@Lin:2015; @Song:2015]. Such assays feature in small-scale experiments and also in high-throughput screens. High-throughput imaging and spectrofluorometric quantification-based methods are increasingly used in the pharmaceutical industry and in academic research aiming to screen large libraries (e.g. potential drug compounds, RNAi, CRISPR/Cas9) or to quantify rare biological events.

Plate-reading equipment is usually supplied with accompanying analytic software, such as  Gen5 (BioTek), Multi-Mode Analysis Software (Molecular Devices), and Multiwell-Analyzer (Multi Channel Systems). There are, however, drawbacks to its use. Usually, the software is bound to the data coming from a specific piece of equipment, it is seldom possible to use the software on multiple computers without making additional license purchases. This is a problem in academic settings in which a single machine is shared and used intensively by several research groups. Also, proprietary software has closed code that cannot be checked or understood in detail by the  user and cannot be expanded or adapted for analyses not envisaged by the equipment manufacturer.

In response to these limitations, some open-source software packages have been developed to serve particular unusual applications such as single-cell migration [@Masuzzo:2017] or differential scanning fluorimetry [@Wang:2012]. Open software with more general capabilities like DRfit [@Hofmann:2019] that run on JAVA platforms has been developed, but they are not ideal for high throughput screening. Many biologists are familiar with the statistical language, R, and R libraries have been created to assist with analysis of data from plate readers. Examples include ‘platetools’ (https://cran.r-project.org/package=platetools), ‘plater’ (https://cran.r-project.org/package=plater) and ‘phenoScreen’ (https://rpubs.com/Swarchal/phenoScreen), each of which addresses particular functions. ‘Plater’, for example, is focused on data formatting, ‘phenoScreen’ on data visualization, and ‘platetools’ on both data visualizing and formatting. The packages are useful, but can be difficult for a beginner to grasp, partly because their documentation seems not to have been written with beginners in mind.

Here, we introduce a new open-source package, 'bioassays', designed to provide a wide range of functions relevant to multi-well plate assays. ‘Bioassays’ is an R package freely available on CRAN [@R:2014], and has supports for formatting, visualizing and analyzing multi-well plate assays. It has functions for handling outliers, for handling multiple data sets with separate blanks, for estimating values from standard curves, for summarizing data and for doing statistical analysis. Moreover, we have provided examples (bioassays-examples) in the vignettes, which will be very helpful for beginners to grasp how this package can be used. Because this is a R package, more experienced users will be able add to its functionality using other R packages.


# Bioassays
Bioassays can handle data from any of the standard multi-well plate format: 6, 12, 24, 96 or 384 well plates. A prerequisite for the package is the need for both input data (\autoref{figure 1}) and metadata (\autoref{figure 2}) to be in comma-separated variable (csv) format; most plate readers can export data as a .csv file. Metadata need  “row” and “col” columns to indicate the location of well.  The package has functions for data extraction (extract_filename), formatting (data2plateformat, plate2df, matrix96 and plate_metadata), visualization (heatplate) and data analysis (reduceblank, estimate, dfsummary, pvalue).

The function extract_filename is useful for extracting specific information from file names, such as compound name, plate number etc., which can be used for automated data analysis. This function provides a very easy way to pass data of this type into the analysis by simply editing file names.

The function data2plateformat converts the data (eg: readings from a 96 well plate) to an appropriate matrix format with suitable row and column names. The plate2df function can convert such data into a data frame with ‘row’ and ‘col’ columns to indicate well position. The function matrix96 can convert both character and numeric columns to a matrix. This function can automatically determine the type of multi-well from which the data is coming, and it have provision to convert negative and NA values to 0, if needed. The plate_metadata function combines plate-specific information (such as compound used, standard concentration, dilution of samples, etc.), to produce unique plate metadata.

For data visualization, the ‘heatplate’ function can be used. This can plot both heat plots (\autoref{figure 3}) and categorical plots (\autoref{figure 4}) automatically, depending on the data type.  The ‘heatplate’ function also has provisions for adjusting the displayed well size, for visual esthetics. Any multi-well plate data can be structured for the plots (function 'heatplate') by using functions 'data2plateformat', 'plate2df', 'matrix96' alone or in combinations (depending on the input data).

The function ‘reduceblank’ can help to reduce blank values from the readings. This function can handle separate blanks for different datasets. The function ‘estimate’ can be used to estimate unknown variable (eg: concentration) based on standard curve. The function ‘dfsummary’ is really versatile for handling multiple data sets. It can group samples and summarize data sets separately. It has additional controls for handling outliers and omitting unwanted data sets. Function ‘pvalue’ can be used to test significance by t-test. It has provisions for asserting control group and level of significance.

# Availability and installation

Bioassays is supported on Windows and macOS. The package can be installed using install.packages command in R. The source code, vignette, datasets and detailed examples on how to use the package are available on CRAN ( https://CRAN.R-project.org/package=bioassays) and GitHub (https://github.com/anwarbio/bioassays).

# Acknowledgment
The author acknowledges Prof Jamie Davies for suggestions and correcting the article. Support from lab members (Jamie Davies Group, Centre for Discovery Brain Sciences, University of Ediburgh) are acknowledged.

# Figures
![Input data format from a 96 well plate reading.\label{figure 1}](figure1.png)

![Metafile data format.\label{figure 2}](figure2.png)

![Heat map of 384 well (normalized values).\label{figure 3}](figure3.png)

![Categorical plot of 384 well plate.\label{figure 4}](figure4.png)

# References
# Contributing

Thanks for your interest in contributing to bioassays.
Here are some guidelines to help make it easier to merge your Pull Request:

* For potentially large changes, please open an Issue first to discuss
* Please follow the [Hadley style guide][style]
* Run `devtools::test()` to run the tests
* (Optional) Add new test(s) in `tests/testthat/`

If you're new to submitting Pull Requests, please read the section [Contribute
to other projects][contribute] in the tutorial [A quick introduction to version
control with Git and GitHub][git-tutorial].

## More about this repository

For the most part, I tried to follow the guidelines from [R packages][r-pkg] by
[Hadley Wickham][hadley]. The unit tests are performed with [testthat][], the
documentation is built with [roxygen2][], and the online package documentation
is created with [rmarkdown](https://github.com/rstudio/rmarkdown). Continuous integration testing is performed for
for macOS by [Travis CI](https://travis-ci.com/github/anwarbio/bioassays/builds), and for Windows
by [AppVeyor](https://ci.appveyor.com/project/anwarbio/bioassays).

The repository contains the files `LICENSE` and `LICENSE.md` to both adhere to
[R package conventions for defining the license][r-exts-licensing] and also to
make the license clear in a more conventional manner (suggestions for
improvement welcome). Directories are standard for R packages
as described in the manual [Writing R Extensions][r-exts].

## Release checklist (for maintainers)

* Bump version with [scripts/bump-version.R](scripts/bump-version.R)
* Update [NEWS.md](NEWS.md): Check `git log` and make sure to reference GitHub
Issues/PRs
* Run [scripts/build.sh](scripts/build.sh) to confirm tests pass locally
* Test on [rhub][]:
    * Have to validate email first with `rhub::validate_email()`. Copy-paste
    token from email into R console.
    * Check on Ubuntu with `rhub::check_on_ubuntu()`
    * Check on macOS with `rhub::check_on_macos()`
    * Check on Windows with `rhub::check_on_windows()`
* Test on [winbuilder][]:
    * Check with R release with `devtools::check_win_release()`
    * Check with R devel with `devtools::check_win_devel()`
* Update [cran-comments.md](cran-comments.md)
* Commit with `git commit -am "Bump version: x.x.x.9xxx -> x.x.x and re-build
docs."`
* Push with `git push origin master` and wait for CI builds to pass
* Tag with `git tag -a vx.x.x`. Summarize [NEWS.md](NEWS.md) entry into bullet
points. Run ` git tag -l -n9` for past examples. Push with `git push origin
--tags`.
* Make a release. On GitHub, go to Releases -> Tags -> Edit release notes. Name
the release "workflowr x.x.x" and copy-paste the Markdown entry from
[NEWS.md](NEWS.md).
* Build tarball with `R CMD build .` and upload to [CRAN submission
site][cran-submit]. You will receive an email to request confirmation, then an
email confirming the package was submitted, and then an email with the test
results. Once it is accepted to CRAN, monitor the [check results][check-results]
for any surprise errors. Also, these builds are when the binaries are built for
Windows and macOS, so they aren't available until they are finished. You will
receive an email once all the Windows binaries are available for download
(devel, release, oldrel).

[appveyor]: https://ci.appveyor.com
[check-results]: https://cran.r-project.org/web/checks/check_results_workflowr.html
[circleci]: https://circleci.com
[Codecov]: https://codecov.io/
[contribute]: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004668#sec011
[covr]: https://github.com/jimhester/covr
[cran-submit]: https://cran.r-project.org/submit.html
[foghorn]: https://cran.r-project.org/package=foghorn
[git-tutorial]: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004668
[hadley]: http://hadley.nz/
[pkgdown]: https://github.com/r-lib/pkgdown
[pt]: https://rstudio.github.io/rstudio-extensions/rstudio_project_templates.html
[r-exts]: https://cran.r-project.org/doc/manuals/R-exts.html
[r-exts-licensing]: https://cran.r-project.org/doc/manuals/R-exts.html#Licensing
[r-pkg]: http://r-pkgs.had.co.nz/
[rhub]: https://r-hub.github.io/rhub/
[roxygen2]: https://github.com/klutometis/roxygen
[style]: http://adv-r.had.co.nz/Style.html
[testthat]: https://github.com/hadley/testthat
[travis]: https://travis-ci.org/
[winbuilder]: https://win-builder.r-project.org/

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Bioassays

<!-- badges: start -->

<!-- badges: end -->

'Bioassays' is an R package, designed to provide a wide range of functions relevant to multi-well plate assays. It can handle data from any of the standard multi-well plate format: 6, 12, 24, 96 or 384 well plates. ‘Bioassays’ can help in formatting, visualizing and analyzing multi-well plate data. It has functions for handling outliers, for handling multiple data sets with separate blanks, for estimating values from standard curves, for summarizing data and for doing statistical analysis. Moreover, it is strongly documented in a manner designed to be easy for even beginners to grasp.
## Installation

You can install the released version of bioassays from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bioassays")
```
The latest development version can be installed from github with the remotes package
``` r
library(remotes)
install_github("anwarbio/bioassays")
```
## Example

Detail’s of various funtions in this package is provided in
‘bioassays-vignette’. Examples on how to use this package is provided
in ‘bioassays-example’ in Vignette.

## License: GPL-3.0
https://github.com/anwarbio/bioassays/blob/master/LICENSE

## Contribute

We love your input! Users may request new features by opening a [GitHub Issue](https://github.com/anwarbio/bioassays/issues), or may contribute their own additions and improvements via a pull request. Similarly, if you run into problems while using this package, or require technical support, do not hesitate to request support through a [GitHub Issue](https://github.com/anwarbio/bioassays/issues). If you use ‘Bioassays’ in your work and would like to further collaborate, I would be more than willing to discuss it over email or [GitHub Issue](https://github.com/anwarbio/bioassays/issues). When you submit code changes, your submissions are understood to be under the same [GPL-3.0](https://github.com/anwarbio/bioassays/blob/master/LICENSE) that covers this project.

An incomplete list of possible improvements:
* Include support for qPCR data analysis
* Include provisions for user template
* Include capability for plotting data

## version 1.0.0
- OCT 2020
---

### NEWS.md setup

- added NEWS.md creation with newsmd

---
title: "Cellular assays using ‘bioassays’ package in R"
author: "Anwar Azad Palakkan, Jamie Davies"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{bioassays-examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(bioassays)
```


# Libraries 

The following packages will be useful. If they are not installed, please install them using **install.packages()**.

```{r libraries, message=FALSE, warning=FALSE, echo=TRUE, eval = FALSE}
library(tcltk)# for selecting the folder for analysis
library(dplyr)
library(ggplot2)# for plotting graphs
library(reshape2)
library(nplr)# for the standard curve fitting
```


# Introduction

In a cell culture lab various cellular assays are performed. The package
"bioassays" helps to analyse the results of these experiments performed in
multiwell plates. The functions in this package can be used to summarise data
from any multiwell plate, and by incorporating them in a loop several plates can
be analyzed automatically. Two examples are shown in this article. All the csv
files used for the examples are available in the `exdata` folder.

# Example 1: Analysing result from a 96 well plate

## Input data

To set up a folder as working directory

```{r directory,echo=TRUE, eval = FALSE}
path1<-tk_choose.dir(getwd(), "Choose the folder for Analysis")
# A window will popup and ask you to select the folder containg the data  
setwd(path1)
```

## Read files

Files can be read using the code example below. The files need to be in .csv
format. In this example file names reflect the type of the assay used. For
example  *"L_DIFF_P3_72HRS.csv"* :*L* is assay code (Lactate assay), *DIFF*
differentiation (type of cells used), *p3*: Plate 3 (each plate represent
specific compounds tested), *72HRS* 72hrs treatment with compound. This
information will help to summarise the results.

```{r files}
filelist <-list("L_HEPG2_P3_72HRS.csv","L_HEPG2_P3_24HRS.csv") # list of files
fno<-1 # file number (in fileslist) that is going to be analyzed
result <- data.frame(stringsAsFactors= FALSE) ## An empty dataframe to dump result
zzz <- data.frame(stringsAsFactors= FALSE) ## An empty dataframe to dump result
```
 
To read the first file

```{r readcsv}
filename<-extract_filename(filelist[fno])[1]
filename
nickname<-extract_filename(filelist[fno], split="_",end=".csv",remove="",sep="")[2]
nickname
```

```{r readcsv2, eval=FALSE}
rawdata<-read.csv(filename,stringsAsFactors = FALSE, strip.white = TRUE, 
                  na.strings = c("NA",""),header = TRUE,skip=1)
head(rawdata)
```

```{r readcsv3, echo=FALSE}
data(rawdata96)
rawdata<-rawdata96
head(rawdata)
```

Reading the metadata file

```{r metafile, eval=FALSE}
metadata<-read.csv("metafile.csv",stringsAsFactors = FALSE,strip.white = TRUE, 
                  na.strings = c("NA",""),header = TRUE)
head(metadata)
```

```{r metafile2, echo=FALSE}
data(metafile96)
metadata<-metafile96
head(metadata)
```

## Rearranging the data

96 well plates were used for the assay, so it is assumed that data is arrayed in
rows A to G and columns 1 to 12 of a 96 well plate. The `data2plateformat`
function is used to label them correctly.

```{r rearrange-plate-reading}
rawdata<-data2plateformat(rawdata,platetype = 96)
head(rawdata)
```

## To convert data into a data.frame.

The `data2plateformat` function uses the column and row names of the `rawdata`
object to coerce it into a data.frame. 

```{r plate2df OD}
OD_df<- plate2df(rawdata)
head(OD_df)
```

## Create a graphical overview

```{r heatmap, normalization, echo = TRUE, eval=TRUE, fig.width=4, fig.height=3.5}
data<-matrix96(OD_df,"value",rm="TRUE")
heatplate(data,"Plate 1", size=5)
```

## Filling metadata file using plate specific details

The example function given below determines the compound, concentration,type and
dilution from the file name.

```{r joining metatadata}
plate_info<-function(file,i){
  
  file<-file[1]
  plate<- extract_filename(file,split = "_",end = ".csv", remove = " ",
                           sep=" ")[5]
  
  if(plate == "P2"){
    compound<-"CyclosporinA"   # Concentration of cyclosporinA used for experiment
    concentration<-c(0,1,5,10,15,20,25,50) # Concentration of cyclosporinA used for experiment
    type<-c("S1","S2","S3","S4","S5","S6","S7","S8") # sample names of corresponding concentration
    dilution<-5
    plate_meta<-list(compound=compound,concentration=round(concentration,2),
                     type=type,dilution=dilution)
  }
  
  
  if(plate == "P3"){
    compound<-"Taxol"
    concentration<-c(0,0.0125,.025,.05,0.1,1,5,10) 
    type<-c("S1","S2","S3","S4","S5","S6","S7","S8")
    dilution<-5
    plate_meta<-list(compound=compound,concentration=round(concentration,2),
                     type=type, dilution=dilution)
  }
  
  
  if(plate =="p4"){
    compound<-c("Cisplatin")
    concentration<-c(0,0.5,2,4,8,16,64,"") 
    type<-c("S1","S2","S3","S4","S5","S6","S7","") 
    dilution <- 5
    plate_meta<-list(compound=compound,concentration=round(concentration,2),
                     type=type,dilution=dilution)
  }
  return(plate_meta)
}

plate_details<-plate_info(filelist,1)
plate_details
```

These 'plate details' can be used to add the metadata with the `plate_metadata` 
function.

```{r metadata compiling}
metadata1<-plate_metadata(plate_details,metadata,mergeby="type")
head(metadata1)
```

For joining the metadata and platelayout dplyr's `inner_join` function can be
used:

```{r metafile-joining, warning=FALSE}
data_DF<- dplyr::inner_join(OD_df,metadata1,by=c("row","col","position"))
assign(paste("data",sep="_",nickname),data_DF) # create a copy of data with plate name
head(data_DF)
```

## Sorting blank wells and reducing blanks
 
Blank values can be subtracted from the measurements with the `reduceblank`
function.

```{r reduce-blank}
data_DF<-reduceblank(data_DF, x_vector =c("All"),
                     blank_vector = c("Blank"), "value")
head(data_DF)
assign(paste("Blkmin",sep="_",nickname),data_DF) # create a copy of data with plate name
```

## Plotting a standard curve

To filter standards

```{r standard-curve-calculation, warning=FALSE}
std<- dplyr::filter(data_DF, data_DF$id=="STD")  
std<- aggregate(std$blankminus ~ std$concentration, FUN = mean )
colnames (std) <-c("con", "OD")
head(std)
```

Calculations for standard curve : nonparametric logistic regression curve

```{r nprc graph, fig.width= 3, fig.height=3,warning=FALSE}

fit1<-nplr::nplr(std$con,std$OD,npars=3,useLog = FALSE) 
#npars = 3 for 3 parametric regression curve

#for graph
x1 <- nplr::getX(fit1); y1 <- nplr::getY(fit1)
x2 <- nplr::getXcurve(fit1); y2 <- nplr::getYcurve(fit1)
plot(x1, y1, pch=15, cex=1, col="red", xlab="Concentration",
     ylab="Mean OD", main=paste("Standard Curve: ", nickname), cex.main=1)
lines(x2, y2, lwd=3, col="seagreen4")
```


To evaluate nonparametric logistic regression fitting

```{r evaluate linear fitting}
params<-nplr::getPar(fit1)$params
nplr::getGoodness(fit1)
```


To estimate the values based on logistic regression fitting

```{r as a fxn estimate1, warning= FALSE}
estimated_nplr<-estimate(data_DF,colname="blankminus",fitformula=fit1,method="nplr")
head(estimated_nplr)
```


Calculations for standard curve : linear regression curve

```{r linear-regression-plot, fig.width= 3, fig.height=3}

fit2<-stats::lm(formula = con ~ OD,data = std)
ggplot2::ggplot(std, ggplot2::aes(x=OD,y=con))+
ggplot2::ggtitle(paste("Std Curve:", nickname))+
ggplot2::geom_point(color="red",size=2)+
ggplot2::geom_line(data = ggplot2::fortify(fit2),ggplot2::aes(x=OD,y=.fitted),
          colour="seagreen4",size=1)+
ggplot2::theme_bw()
```


To evaluate the fit of the linear model


```{r linear-fit-summary, warning=FALSE}
conpred<-estimate(std,colname="OD",fitformula=fit2,method="linear")
compare<-conpred[,c(1,3)]
corAccuracy<-cor(compare)[1,2]
corAccuracy ## Correlation accuracy of regression line
summary(fit2)

```


To estimate the values based on linear curve

```{r as a fxn estimate3, warning= FALSE}
estimated_lr<-estimate(data_DF,colname="blankminus",fitformula=fit2,
                       method="linear")
head(estimated_lr)
```


For multiply estimated by dilution

```{r multiply by dilution1, warning=FALSE}
estimated_lr$estimated2 <- estimated_lr$estimated * estimated_lr$dilution
head(estimated_lr)
```

## Summarise the data
 
 For summarising the "estimated_lr" based on "id" and "type". 

```{r plate-summary, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, out.width=50}

result<-dfsummary(estimated_lr,"estimated2",c("id","type"),
        c("STD","Blank"),"plate1", rm="FALSE",
        param=c(strict="FALSE",cutoff=40,n=12))

result
```

## Statistical test

For the t test (S1 as control)

```{r t test 1}
pval<-pvalue(result, control="S1", sigval=0.05)
head(pval)

```

          
\newpage

# Example 2: Analysing result from a 384 well plate.

This example data contain result of dose response of few drugs (drug1, drug2,
drug3, drug4) at three concentrations (C1,C2,C3) from two different cell lines
(hepg2 and huh7)

## Input data.

Reading the spectrophotometer readings from a csv file

```{r 384 rawdata, eval=FALSE}
rawdata2<-read.csv("384.csv",stringsAsFactors = FALSE,strip.white = TRUE, 
                  na.strings = c("NA",""),header = TRUE,skip=1)
dim(rawdata2)
head(rawdata2)

```

```{r 384 rawdata2, echo= FALSE }
data(rawdata384)
rawdata2<-rawdata384
dim(rawdata2)
head(rawdata2)
```

Reading the metadata file

```{r 384 metadata, eval=FALSE}
metadata2<-read.csv("metafile_384_plate.csv",stringsAsFactors = FALSE,
                    strip.white = TRUE, na.strings = c("NA",""), header = TRUE)
head(metadata2)

```

```{r 384 metadata2, echo=FALSE}
data(metafile384)
metadata2<-metafile384
head(metadata2)
```


## Rearranging the data.

Renaming  rows and columns of the `rawdata2` object

```{r 384 rearrange}
rawdata2<-data2plateformat(rawdata2,platetype = 384)
head(rawdata2)

```

Converting the data into a data.frame.

```{r 384 plate2df}
OD_df2 <- plate2df(rawdata2)
head(OD_df2)
```

Create a graphical overview of the plate.

```{r 384overview1, echo = TRUE, eval=TRUE, fig.width=4, fig.height=3.5}
data2<-matrix96(OD_df2,"value",rm="TRUE")
heatplate(data2,"Plate 384", size=1.5)
```

Joining the `metadata2` and `OD_df2` objects

```{r 384 innerjoin}
data_DF2<- dplyr::inner_join(OD_df2,metadata2,by=c("row","col","position"))
head(data_DF2)
```

Display the cell categories:

```{r 384overview2, echo = TRUE, eval=TRUE, fig.width=4, fig.height=3.5}
data3<-matrix96(data_DF2,"cell",rm="TRUE")
heatplate(data3,"Plate 384", size=2)
```

Discplay the compounds

```{r 384overview3, echo = TRUE, eval=TRUE, fig.width=4, fig.height=3.5}
data4<-matrix96(data_DF2,"compound",rm="TRUE")
heatplate(data4,"Plate 384", size=2)
```


## Sorting blank wells and reducing blanks.

The data contain separate blanks wells for drug1, drug2, drug3 and drug4
(blank1, blank2,blank3 and blank 4 respectively) that can be used for background
correction.

```{r 384 blank}
data_blk<-reduceblank(data_DF2, 
x_vector=c("drug1","drug2","drug3","drug4"),
blank_vector = c("blank1","blank2","blank3","blank4"), "value")
dim(data_blk)
```

```{r 384 blank2}
head(data_blk)
```

For summarising the result in the order cell, compound,concentration,type and by
omitting blanks.

```{r 384 summary 1}
result2<-dfsummary(data_blk,"blankminus",
        c("cell","compound","concentration","type"),
        c("blank1","blank2","blank3","blank4"),
        nickname="384well", 
        rm="FALSE",param=c(strict="FALSE",cutoff=40,n=12))
head (result2)
dim (result2)

```

The following example summarises the results in the following order: cell,
compound and concentration, and by omitting blanks (all blanks are marked as "B"
in concentration), drug 2 and huh7.

```{r 384 summary 2}
result3<-dfsummary(data_blk,"blankminus",
        c("cell","compound","concentration"),
        c("B","drug2","huh7"),
        nickname="", 
        rm="FALSE",param=c(strict="FALSE",cutoff=40,n=12))
head (result3)
dim (result3)

```

## Statistical testing

The following command will perform a t-test:

```{r 384 ttest, echo=TRUE}
pvalue<-pvalue(result3,"C3",sigval=0.05)
pvalue

```
---
title: "Vignette of package bioassays"
author: "Anwar Azad Palakkan, Jamie Davies"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{bioassays-vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(bioassays)
```

# Introduction

In a cell culture lab various cellular assays are performed. The package
"bioassays" will help to analyse the results of these experiments performed in
multiwell plates. The usage of various functions in the "bioassays" package is
provided in this article.

The functions in this package can be used to summarise data from any multiwell
plate, and by incorporating them in a loop several plates can be analyzed
automatically. Two examples are also provided in the article.

# prerequisite

The output reading from the instrument (eg.spectrophotometer) should be in a
matrix format. An example data file (csv format) is shown below. If the data is
in .xls/.xlsx format `read_excel` function from the `readxl` package can be used.


```{r rawdata}
data(rawdata96)
head(rawdata96)
```

A metadata file is needed for the whole experiment. The `row` and `col` columns
must be present in the metadata file to indicate the location of well. An
example is given below.

```{r metadata}
data(metafile96)
head(metafile96)
```

# Functions

## 1. Function: extract_filename

`extract_filename` helps to extract information from the file name. syntax is
`extract_filename(filename,split = " ",end = ".csv", remove = " ", sep="-")`.
**filename** is the file name. **split** is the portions at which the name has
to be split (default is space " "). **end** is the extension of file name that
need to be removed (default is ".csv"). **remove** is the portion from the file
name that need to be omitted after splitting (default is space " "). **sep** add
a symbol between separate sections, default is "-".

This function is useful for extracting specific information from file names,
like compound name, plate number etc, to provide appropriate analysis.

For e.g.

```{r extract_fname eg1}
extract_filename("L HEPG2 P3 72HRS.csv")
```

```{r extract_fname eg2}
extract_filename("L HEPG2 P3 72HRS.csv", split=" ",end=".csv",remove="L",sep="")
```

## 2. Function: rmodd_summary

`rmodd_summary` helps to remove outliers and summarises the values from a given
set of values. The syntax is `rmodd_summary(x, rm = "FALSE", strict= "FALSE",
cutoff=80,n=3)`. **x** is a numeric vector. **rm = TRUE** if want to remove
outliers. If **strict = FALSE** values above/below 1.5 IQR are omitted (outliers
omitted). If **strict = TRUE** more aggresive outlier removal is used to bring
the %cv below the **cutoff**. **n** is the minimum number of samples you need
per group if more aggressive outlier removal is used.

For example:

```{r rmodd_summary eg1}
x<- c(1.01,0.98,0.6,0.54,0.6,0.6,0.4,3)
rmodd_summary(x, rm = "FALSE", strict= "FALSE", cutoff=80,n=3)
```

```{r ermodd_summary eg2}
rmodd_summary(x, rm = "TRUE", strict= "FALSE", cutoff=80,n=3)
```

```{r rmodd_summary eg3}
rmodd_summary(x, rm = "TRUE", strict= "TRUE", cutoff=20,n=5)
```


## 3. Function: data2plateformat

`data2plateformat` converts the data (eg: readings from a 96 well plate) to an
appropriate matrix format. Syntax is `data2plateformat(data, platetype = 96)`.
**data** is the data to be formatted. **platetype** is the plate from which the
data is coming. It can take 6, 12, 24, 96 or 384 values to represent the
corresponding multiwell plate type.

For e.g. To rename columns and rows of `rawdata96` to right format.


```{r data2plateformat eg2 }
rawdata<-data2plateformat(rawdata96,platetype = 96)
head(rawdata)
```

## 4. Function: plate2df

`plate2df` coereces a matrix with values from multiwell plates into data.frame.
The function uses column names and row names of `datamatrix` (2D data of a
mutliwell plate) and generates a dataframe with row, col (column) and position
indices. The `value` column represents the corresponding value in the
`datamarix`.

Syntax is `plate2df(datamatrix)`. **datamatrix** is the data in matrix format. 

For eg.

```{r plate2df eg1}
OD_df <- plate2df(rawdata)
head(OD_df)
```

## 5. Function: matrix96

`matrix96` helps to convert a data.frame into a matrix. The syntax is
`matrix96(dataframe,column,rm="FALSE")`. **dataframe** is the data.frame to be
formatted. The data.frame must contain `row` and `col` columns. **column** is the
name of column that needs to be converted into a matrix. If **rm= "TRUE"** then
-ve and NA are set to 0.


For e.g.

```{r matrix96 eg2}
matrix96(OD_df,"value")
```

```{r matrix96 eg3}
matrix96(OD_df,"position")
```

## 6. Function: plate_metadata

`plate_metadata` combines the plate-specific information (like compound used,
standard concentration, dilution of samples, etc) and metadata into  unique plate
metadata. The syntax is `plate_metadata(plate_details,
metadata,mergeby="type")`. **plate details** is the plate specific information
that needs to be added to metadata. **metadata** is the metadata for the whole
experiment. **mergeby** is the column that is common to both metadata and
`plate_meta` (this column will be used to merge the information).

For eg. An incomplete meta data

```{r plate_metadata eg 2 }
head(metafile96)
```

Plate specific details are.

```{r plate_metadata eg 3, echo=TRUE}
plate_details <- list("compound" = "Taxol",
                "concentration" = c(0.00,0.01,0.02,0.05,0.10,1.00,5.00,10.00),
                "type" = c("S1","S2","S3","S4","S5","S6","S7","S8"),
                "dilution" = 1)
```


Using the plate-specific info, the metadata can be filled by calling 
`plate_metadata` function.


```{r plate_metadata eg 5, message=FALSE, warning=FALSE}
plate_meta<-plate_metadata(plate_details,metafile96,mergeby="type")
head(plate_meta)
```

To join both plate_meta and OD_df, the `dplyr::inner_join` function can be used.


```{r plate_metadata eg 6, message=FALSE, warning=FALSE}
data_DF<- dplyr::inner_join(OD_df,plate_meta,by=c("row","col","position"))

head(data_DF)
```


## 7. Function: heatplate

`heatplate` creates a heatmap representation of a multiwell plate. The syntax is `heatplate(datamatrix,name,size=7.5)`. **datamatrix** is the data in matrix
format. An easy way to create this is by calling the `matrix96` function, as 
explained before. **name** is the name to be given for heatmap, 
**size** is the size of each well in the heatmap (default is 7.5). 

This function will produce a heatmap of normalized values if the `variable` is 
numeric. If it is a factor, it will simple provide a coloured  categorical plot.

eg 1. Categorical plot

```{r heatplate eg:1, warning= FALSE}
datamatrix<-matrix96(metafile96,"id")
datamatrix
```


```{r heatplate eg:2, warning= FALSE, fig.width=4, fig.height=3.5}
heatplate(datamatrix,"Plate 1", size=5)
```

eg 2. Heatmap

```{r heatplate eg:3, warning=FALSE}
rawdata<-data2plateformat(rawdata96,platetype = 96)
OD_df<- plate2df(rawdata)
data<-matrix96(OD_df,"value")
data
```


```{r heatplate eg:4, warning=FALSE, fig.width=4, fig.height=3.5}
heatplate(data,"Plate 1", size=5)
```


## 8. Function: reduceblank

`reduceblank` helps to reduce blank values from the readings.

The syntax is `reduceblank (dataframe,x_vector,blank_vector,y)`. **dataframe**
is the data. **x_vector** are the entries for which the blank has to be
subtracted If all entries have to subtracted, then use "All". x_vector should be
a vector eg: c("drug1","drug2",drug3" etc). **blank_vector** is a vector of
blank names whose value has to be subtracted eg:
c("blank1","blank2","blank3","blank4")). This function will subtract the first
blank vector element from first x_vector element and so on. **y** is the column
name where the action will take place. y should be numeric. The results will
appear as a new column named `blankminus`.

For eg.

```{r reduceblank eg1, warning=FALSE}
data_DF<-reduceblank(data_DF, x_vector =c("All"),
                     blank_vector = c("Blank"), "value")
head(data_DF)

```

## 9. Function: estimate

`estimate` estimates an unknown variable (eg: concentration) based on the
standard curve. Syntax is
`estimate(data=dataframe,colname="blankminus",fitformula=fit,
methord="linear/nplr")`. **data** is the dataframe which needs to be evaluated.
**colname** is the column name for which the values has to be estimated.
**fitformula** represents the model to be fit **methord** specifies if linear or
nonparametric logistic curve was used for the fitformula.

For eg: data_DF is a dataframe for which the concentration has to be estimated based on the value of blankminus.

For filtering the `standards`

```{r standards}
std<- dplyr::filter(data_DF, data_DF$id=="STD")  
std<- aggregate(std$blankminus ~ std$concentration, FUN = mean )
colnames (std) <-c("con", "OD")
head(std)
```

To fit a standard curve:

fit1 is the 3 parameter logistic curve model and fit2 is the linear regression
model. The appropriate one for your experiment can be used.

```{r fitmodels}
fit2<-stats::lm(formula = con ~ OD,data = std)# linear model
fit1<-nplr::nplr(std$con,std$OD,npars=3,useLog = FALSE)#  nplr, 3 parameter model
```

For estimating the concentration using linear model

```{r nplr estimating nplr, message=FALSE, warning=FALSE}
estimated<-estimate(data_DF,colname="blankminus",fitformula=fit2,method="linear")
head(estimated)
```

For estimating the concentration using nplr methord

```{r nplr estimating linear, message=FALSE, warning=FALSE}
estimated2<-estimate(data_DF,colname="blankminus",fitformula=fit1,method="nplr")
head(estimated2)
```

## 10. Function: dfsummary

`dfsummary()` summarizes the data.frame (based on a column). It has additional
controls to group samples and to omit variables not needed. The syntax is
`dfsummary(dataframe,y,grp_vector,rm_vector,nickname,rm="FALSE",param)`.
**dataframe** is the data. **y** is the numeric variable (column name) that has
to be summarized. **grp_vector** is a vector of column names, determining which
samples will be grouped. The order of the elements in grp_vector determines the
order of grouping. **rm_vector** is the vector of items to be omitted before
summarizing. **nickname** is the name for the output data.frame. Set
**rm="FALSE"** if outliers should not be removed. To apply more stringent
thresholds for removing outliers, parameters can be provided in the  **param**
vector. **param** has to be entered in the format
`c(strict="TRUE",cutoff=40,n=12)`. For details please refer to the
`rmodd_summary `function.

In the following example, the data has to be summarized based on the "type"
column. "estimated" values are summarized. Samples are grouped as per "id".
"STD" and "Blank" values need to be omitted. Outliers are retained (rm="FALSE").
The nickname for the plate is "plate1".

```{r summary 3, message=FALSE, warning=FALSE, out.width=40}
result<-dfsummary(estimated,"estimated",c("id","type"),
        c("STD","Blank"),"plate1", rm="FALSE",
        param=c(strict="FALSE",cutoff=40,n=12))

result
```

## 11. Function: pvalue

The `pvalue()` function applies a t-test to the result dataframe. The syntax is
`pvalue(dataframe,control,sigval)`. **dataframe** is the result of dfsummary.
**control** is the group that is considered as control, **sigval** is the pvalue
cutoff (a value below this is considered as significant).

For example:

```{r pvalue eg1, warning=FALSE, message = FALSE}
pval<-pvalue(result, control="S8", sigval=0.05)
head(pval)
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metafile96.R
\docType{data}
\name{metafile96}
\alias{metafile96}
\title{metadata of 96 well plate.}
\format{
A data frame with 96 rows and 7 variables.
\describe{
\item{row}{Row number of multi well well plate.}
\item{col}{Column number of multi well plate.}
\item{position}{Well position address of multi well plate.}
\item{id}{Type of sample. 'STD' represent standards, 'sample' represent samples.}
\item{type}{Type of sample. 'STD1','STD2','STD3' etc represent different standards, 'S1','S2','S3' represent different samples.}
\item{concentration}{concentration of different standards (mg/ml).}
\item{dilution}{dilution of samples used for assay.}
}
}
\source{
{User generated metadata of 96 well plate.}
}
\usage{
metafile96
}
\description{
A dataset containing metadata.
}
\keyword{datasets}
\name{plate_metadata}
\alias{plate_metadata}

\title{Combining Plate Specific Information with Metadata}
\description{
plate_metadata combine the plate specific information (like compounds used, standard concentration, dilution of samples, etc) and metadata, to produce a plate specific metadata.
}
\usage{
plate_metadata (plate_details, metadata, mergeby = "type")
}

\arguments{
  \item{plate_details}{plate specific information that need to be added to metadata}
  \item{metadata}{metadata for whole experiment}
  \item{mergeby}{column that is common to both metadata and plate_meta (as a string in "")}
}
\details{
plate_details need to be in a list format. Metadata should have a 'row' and 'col' columns representing the row and column names of the corresponding multi well plate.
}
\value{
 A dataframe. Each element of 'plate_details' will appear as a new column to the left of 'metadata'

}

\author{
A.A Palakkan
}

\examples{
## loading data
data(metafile96)
plate_details <- list("compound" = "Taxol",
                "concentration" = c(0.00,0.01,0.02,0.05,0.10,1.00,5.00,10.00),
                "type" = c("S1","S2","S3","S4","S5","S6","S7","S8"),
                "dilution" = 1)

## eg:1 filling metadata96 using plate_details
plate_meta<-plate_metadata(plate_details,metafile96,mergeby="type")
head(plate_meta)


}

\keyword{ manip }
\name{plate2df}
\alias{plate2df}

\title{Format Matrix Type 2D Data of Multi well Plate as Dataframe }

\description{
This function uses column names and row names of  'datamatrix' (2D data of a mutli well plate) and generate a dataframe with row, col (column) and position indices. The 'value' column represent corresponding value in the 'datamarix'.
}
\usage{
plate2df(datamatrix)
}

\arguments{
  \item{datamatrix}{datamatrix is the 2D data of a mutli well plate. Usually the result of \code{\link{data2plateformat}}:}
}

\value{
A dataframe with 4 columns. Number of rows is equal to the number of wells (plate type of 'datamatrix'). The columns represent
\item{row }{Row number of the entry}
\item{col }{Column number of the entry}
\item{position }{Position (Row+column number) of the entry}
\item{value}{Individual entries in the 'datamatrix'}

}

\author{
A.A Palakkan
}

\examples{
## loading data
data(rawdata24,rawdata96,rawdata384)

## eg:1 spectrophotometer reading from 24 well plate in dataframe format
datamatrix<- data2plateformat(rawdata24, platetype = 24)
head(plate2df(datamatrix))

## eg:2 spectrophotometer reading from 96 well plate in dataframe format
datamatrix<- data2plateformat(rawdata96, platetype = 96)
head(plate2df(datamatrix))

## eg:3 spectrophotometer reading from 384 well plate in dataframe format
datamatrix<- data2plateformat(rawdata384, platetype = 384)
head(plate2df(datamatrix))

}

\keyword{ manip }

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metafile384.R
\docType{data}
\name{metafile384}
\alias{metafile384}
\title{metadata of 384 well plate.}
\format{
A data frame with 384 rows and 8 variables.
\describe{
\item{row}{Row number of multi well well plate.}
\item{col}{Column number of multi well plate.}
\item{position}{Well position address of multi well plate.}
\item{cell}{Type of cells used for the assay.}
\item{compound}{Different drugs (drug1,drug2,etc) used for the assay.}
\item{concentration}{'C1','C2','C3' etc represent different concentration used for the same compound. 'B' represent blank wells}
\item{type}{'treated' and 'untreated' shows if the wells had received pretreatment (example:inhibitors) or not. 'Blank1','Blank2','Blank3' etc represent separate blanks for different drugs.}
\item{dilution}{dilution of samples used for the assay.}
}
}
\source{
{User generated metadata of the 384 well plate.}
}
\usage{
metafile384
}
\description{
A dataset containing metadata.
}
\keyword{datasets}
\name{heatplate}
\alias{heatplate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Heatmap of multi well plate}
\description{
This function generate a heatmap (for numeric vector) or categorical plot (character vector) of multi well plate}
\usage{
heatplate(datamatrix, name, size = 7.5)
}

\arguments{
  \item{datamatrix}{data in matrix format. An easy way to create this is by calling \code{\link{matrix96}}}
  \item{name}{name to be given for the heatmap}
  \item{size}{plot size for each well in the heatmap (default is 7.5)}
}
\details{
Heat map can be generated for any multi well plate data in matrix format (datamatrix). The columns and rows of datamatrix should be labelled appropriately using \code{\link{matrix96}}. A heatplot is generated if datamatrix is numeric, but a categorical plot is generated if datamatrix is a character matrix.
}
\value{
A graphical plot.

}

\author{
A.A Palakkan
}

\examples{
## loading data
data(metafile96, rawdata96,rawdata384)
rawdata96 <- data2plateformat(rawdata96,platetype = 96)
rawdata384 <- data2plateformat(rawdata384,platetype = 384)

## eg:1 heat map of rawdata96
data<-matrix96(plate2df(rawdata96),"value")
heatplate(data,"Plate 1", size=5)

## eg:2 heat map of rawdata96 can also be called as
heatplate(as.matrix(rawdata96),"Plate 1", size=5)

## eg:3 heat map of rawdata384
heatplate(as.matrix(rawdata384),"Plate 1", size=2)

## eg:4 catagorical map of metafile96 (column:id)
data<-matrix96(metafile96,"id")
heatplate(data,"Plate 1", size=5)

}

\keyword{ hplot }
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bioassays.R
\docType{package}
\name{bioassays}
\alias{bioassays}
\title{Bioassays: A Package for Analyzing Multi Well Plate Bioassay's.}
\description{
The bioassays package provides three categories of important functions for extracting, formatting, plotting and analyzing.
}
\section{Extract functions}{

The function that help to extract information from file name is \cr
\code{\link{extract_filename}}'
}

\section{Format functions}{

The function to format data are \cr
\code{\link{data2plateformat}}\cr
\code{\link{plate2df}}\cr
\code{\link{matrix96}}\cr
\code{\link{plate_metadata}}
}

\section{Plot functions}{

The function to plot is\cr
\code{\link{heatplate}}
}

\section{Analysis functions}{

The function to analyse data are\cr
\code{\link{reduceblank}}\cr
\code{\link{estimate}}\cr
\code{\link{dfsummary}}\cr
\code{\link{pvalue}}
}


\name{dfsummary}
\alias{dfsummary}

\title{Summarize a Dataframe After Grouping Samples}
\description{
This function summarize the dataframe (based on a column). It has additional controls to group samples and to omit variables not needed.
}
\usage{ dfsummary(dataframe, y, grp_vector, rm_vector, nickname, rm="FALSE", param)
}

\arguments{
  \item{dataframe}{data in dataframe format}
  \item{y}{column name whose values has to be summarized (column elements need to be numeric)}
  \item{grp_vector}{a character vector of column names whose order indicate the order of grouping}
  \item{rm_vector}{a character vector of items that need to be omitted before summarizing}
  \item{nickname}{label name for the entries in output dataframe}
  \item{rm}{rm = “FALSE” if outliers not to be removed. rm = “TRUE” If outliers to be removed}
\item{param}{a vector of parameters for more stringent outlier removal. param has to be entered in the format c(strict, cutoff, n). For details please refer \code{\link{rmodd_summary}}:.}
}
\details{
This function first remove 'rm_vector' elements from the 'dataframe'. Samples are grouped (each level of a 'grp_vector' element as separate group) and sorted (based on 'grp_vector' elements order). column 'y' is then summarized for each group (please refer \code{\link{rmodd_summary}}: for details).
}

\value{
A dataframe. First columns are named as grp_vector elements. Followed by a 'label' column (element is 'nickname').This 'label' column will be useful when analyzing multiple plates. Summary statistics of 'y' appear as columns: N (number of samples/group), Mean (average/group), SD (standard deviation/group) and  CV (percentage cv/group).

}

\author{
A.A Palakkan

}

\examples{
## loading data
data(metafile384, rawdata384)
rawdata<-plate2df(data2plateformat(rawdata384,platetype = 384))
data_DF2<- dplyr::inner_join(rawdata,metafile384,by=c("row","col","position"))

## eg:1 summarising the 'value' after grouping samples and omitting blanks.
  # grouping order cell, compound, concentration and type.

result2 <- dfsummary(data_DF2,y = "value",
           grp_vector = c("cell","compound","concentration","type"),
           rm_vector = c("blank1","blank2","blank3","blank4"),
           nickname = "384well",
           rm = "FALSE",param = c(strict="FALSE",cutoff=40,n=12))

}

\keyword{ arith }

\name{data2plateformat}
\alias{data2plateformat}

\title{Renaming column and Row of Multiwell Data to Match Plate Format}

\description{Convert the data (example: readings from mutli well plate) to appropriate plate format by renaming column and rownames.}

\usage{data2plateformat(data, platetype = 96)}

\arguments{
  \item{data}{Matrix data to be formatted.}
  \item{platetype}{Plate from which the data is coming. It can take 6, 12, 24, 96 and 384 values to represent the corresponding multi well plate.}
  }

\details{This function will label the columns and rows correctly to match the plate format, and discard the extras. For example, if the 'data' is coming from a a '96' well plate ('platetype'), the function will rename rows as A to H and columns as 1 to 12. Extra columns and rows of 'data' is discarded.}

\value{A data frame with columns and rows matching (label and numbers) the mutli well plate format.}

\author{
A.A Palakkan
}

\examples{
## loading data
data(rawdata24,rawdata96,rawdata384)

## eg:1 spectrophotometer reading from 24 well plate
data2plateformat(rawdata24, platetype = 24)

## eg:2 spectrophotometer reading from 96 well plate
data2plateformat(rawdata96, platetype = 96)

## eg:3 spectrophotometer reading from 384 well plate
data2plateformat(rawdata384, platetype = 384)

}

\keyword{ manip }

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rawdata384.R
\docType{data}
\name{rawdata384}
\alias{rawdata384}
\title{rawdata of 384 well plate.}
\format{
A data frame with 16 rows and 25 columns: First column shows 'row name'.
}
\source{
{Spectrophotometer output reading (OD) of 384 well plate.}
}
\usage{
rawdata384
}
\description{
A dataset of spectrophotometer readings (OD).
}
\keyword{datasets}
\name{reduceblank}

\alias{reduceblank}

\title{Reduce Blank Values}

\description{This function can reduce 'blank' value from readings. Can handle separate blanks for separate groups in the dataframe}

\usage{
reduceblank(dataframe, x_vector, blank_vector, y)
}

\arguments{
  \item{dataframe}{Data in the form of dataframe.}
  \item{x_vector}{A character vector of groups/entries for which the blank has to be reduced.}
  \item{blank_vector}{A character vector of blank names whose value has to be reduced.}
  \item{y}{Name of the column (column should be numeric in nature) whose values has be reduced.}
}
\details{
 This function will reduce the first blank vector element from first x_vector element and so on.}

\value{
A dataframe with a new column 'blankminus' (result of the blankminus function) added to the right.
}

\author{
A.A Palakkan
}

\examples{
## loading data
data(metafile384, rawdata384)
rawdata<-plate2df(data2plateformat(rawdata384,platetype = 384))
data_DF2<- dplyr::inner_join(rawdata,metafile384,by=c("row","col","position"))

## eg:1 reduce blanks of data_DF2.
  # reduce seperate blanks (mean of blank wells) for drug1, drug2, drug3 and drug4.
  #blanks are blank1, blank2, blank3 and blank4 respectively for different drugs.

data_blk<-reduceblank(data_DF2,
          x_vector=c("drug1","drug2","drug3","drug4"),
          blank_vector = c("blank1","blank2","blank3","blank4"), "value")


}

\keyword{ manip}

\name{extract_filename}
\alias{extract_filename}

\title{Extract Information From File Name}

\description{
This function split a string (file name) as per the requirement of the user. It is useful to extract  informations like compound name, plate number etc from the file name.}
\usage{extract_filename(filename,split = " ",end = ".csv", remove = " ", sep = "-")}

\arguments{
  \item{filename}{name of the file (string).}
  \item{split}{regular expressions at which filename has to be split to create different sections.}
  \item{end}{extension (end portion) of filename that need to be removed.}
  \item{remove}{section that need to be omitted after splitting the filename.}
  \item{sep}{symbol to be added to separate sections (obtained after splitting) before combining. default is ”-"}

  }

\value{
A character vector. First element is the unsplit 'filename'. Second element is the processed 'filename'.Other elements are different sections after splitting the 'filename'.

}

\author{
Anwar Azad Palakkan
}

\examples{
extract_filename("L HEPG2 P3 72HRS.csv")
extract_filename("L HEPG2 P3 72HRS.csv", split=" ",end=".csv",remove="L",sep="")


}

\keyword{ character }

\name{matrix96}
\alias{matrix96}

\title{Formatting Long Dataframe in to a Matrix Layout of Multi well Plate}
\description{
This function format a long dataframe (with col and row columns) in to a multiwell plate matrix layout.
}
\usage{
matrix96 (dataframe, column, rm = "FALSE")
}

\arguments{
  \item{dataframe}{dataframe to be formatted}
  \item{column}{name of column (as a string in "") that need be converted as a matrix}
  \item{rm}{If rm = “TRUE” then -ve and NA are assigned as 0}

}
\details{
The 'dataframe' to be formatted should have a 'col' and 'row' columns representing the column and rowname of the corresponding multiwell plate.
}
\value{
A matrix data with row and column names corresponding to multiwell plate

}

\author{
A.A Palakkan
}

\examples{
## loading data
data(rawdata96, metafile96, metafile384)
rawdata<- data2plateformat(rawdata96, platetype = 96)
rawdata<- plate2df(rawdata)

## eg:1 rawdata to matrix format (column: value)
matrix96(rawdata,"value")

## eg:2 metafile96 to matrix format (column: id)
matrix96(metafile96,"id")

## eg:3 metafile384 to matrix format (column: cell)
matrix96(metafile384,"cell")


}

\keyword{manip }

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rawdata24.R
\docType{data}
\name{rawdata24}
\alias{rawdata24}
\title{rawdata of 24 well plate.}
\format{
A data frame with 4 rows and 7 columns: First column shows 'row name'. .
}
\source{
{Spectrophotometer output reading (OD) of 24 well plate.}
}
\usage{
rawdata24
}
\description{
A dataset of spectrophotometer readings (OD).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rawdata96.R
\docType{data}
\name{rawdata96}
\alias{rawdata96}
\title{rawdata of 96 well plate.}
\format{
A data frame with 8 rows and 13 columns: First column shows 'row name'.
}
\source{
{Spectrophotometer output reading (OD) of 96 well plate.}
}
\usage{
rawdata96
}
\description{
A dataset of spectrophotometer readings (OD).
}
\keyword{datasets}
\name{estimate}
\alias{estimate}

\title{Estimating Variable Based on Standard Curve}

\description{ This function will estimate the unknown variable (example: concentration) based on a standard curve.}

\usage{ estimate (data, colname = "blankminus", fitformula = fiteq, method = "linear/nplr")
}

\arguments{
  \item{data}{data in dataframe format}
  \item{colname}{column name whose values has to be estimated}
  \item{fitformula}{formula used for fitting standard curve}
  \item{method}{method = "linear" if standard curve is linear in nature. method = "nplr" if standard curve is nonparametric logistic curve.}

}
\details{
For linear standard curve 'fitformula' need to generated using \code{\link[stats]{lm}}.
For nonparametric logistic curve 'fitformula' need to generated using \code{\link[nplr]{nplr}}.
}
\value{
A dataframe with estimated values added to right as a new column "estimated".

}

\author{
A.A Palakkan
}

\examples{
## loading data
data(data_DF1)

## Filtering standards
std<- dplyr::filter(data_DF1, data_DF1$id=="STD")
std <- aggregate(std$blankminus ~ std$concentration, FUN = mean )
colnames (std) <-c("con", "OD")

## 3-parametric regression curve fitting
fit1<-nplr::nplr(std$con,std$OD,npars=3,useLog = FALSE)

## Linear regression curve fitting
fit2<- stats::lm(formula = con ~ OD,data = std)

## Estimating the 'blankminus'
## eg:1 Based on nonparametric logistic regression fitting
estimated_nplr <- estimate(data_DF1,colname = "blankminus",fitformula = fit1,method = "nplr")

## eg:2 Based on linear regression fitting
estimated_lr<-estimate(data_DF1,colname="blankminus",fitformula=fit2,method="linear")

}

\keyword{ math }

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_DF1.R
\docType{data}
\name{data_DF1}
\alias{data_DF1}
\title{Data of 96 well plate.}
\format{
A data frame with 96 rows and 10 variables
\describe{
\item{row}{Row number of multi well well plate.}
\item{col}{Column number of multi well plate.}
\item{position}{Well position address of multi well plate.}
\item{value}{Spectrophotometer reading (OD).}
\item{id}{Type of sample.'STD' represent standards and 'sample' represent samples.}
\item{type}{Type of sample.'STD1','STD2','STD3' etc represent different standards.'S1','S2','S3' represent different samples.}
\item{dilution}{dilution of samples used for the assay.}
\item{concentration}{Concentration of respective standards.}
\item{compound}{Compound used for the assay.}
\item{blankminus}{Blank reduced OD (value - mean(blank))}
}
}
\source{
{User generated dataframe of the 96 well plate.}
}
\usage{
data_DF1
}
\description{
A complete dataset containing both metadata and spectrophotometer reading.
}
\keyword{datasets}
\name{rmodd_summary}
\alias{rmodd_summary}
\title{Summarise a Numerical Vector with Control on Outlier Removal}
\description{
Summarise a numerical vector with control on how the outliers has to be treated.
}
\usage{
rmodd_summary(x, rm = "FALSE", strict = "FALSE", cutoff = 80,n = 3)
}

\arguments{
  \item{x}{ numerical vector}
  \item{rm}{ if rm = "TRUE" outliers are omitted. If rm = "FALSE" all elements in the vector are considered for summarising}
  \item{strict}{ if strict = "FALSE" outliers are omitted based on IQR rule. If strict = "TRUE" more aggressive outlier omitting method is used to bring CV below a cutoff value}
  \item{cutoff}{ cv cutoff value for the aggressive outlier removal}
  \item{n}{minimum number of samples needed}

}
\details{
In IQR rule (ie when strict = "FALSE") those values above 'Q3 + 1.5 IQR' and those below 'Q1 - 1.5 IQR' is considered as outlier.
For the aggressive outlier removal (ie when strict = "TRUE") those values above 90th percentile and  below 10th percentile are removed consecutively till the cv fall below the 'cutoff' or only the minimum number of samples is leftover (whichever happens first halt the loop).
}
\value{
  A  numeric vector of length 5 with the elements representing

 \item{mean }{the average of samples}
 \item{median }{the median of samples}
\item{n }{number of samples}
\item{sd }{standard deviation of samples}
\item{cv }{percentage cv of samples}
}

\author{
A.A Palakkan
}

\examples{
## data set x
x <- c(1.01,0.98,0.6,0.54,0.6,0.6,0.4,3)

## summarising without removing outliers
rmodd_summary(x, rm = "FALSE", strict= "FALSE", cutoff=80, n=3)

## summarising after removing outliers (IQR methord)
rmodd_summary(x, rm = "TRUE", strict= "FALSE", cutoff=20, n=5)

## summarising after removing outliers (Stringent to reduce cv)
rmodd_summary(x, rm = "TRUE", strict= "TRUE", cutoff=20, n=5)

}

\keyword{ arith }

\name{pvalue}
\alias{pvalue}

\title{t-Test on Summary Dataframe}

\description{
This function calculate the significance (t-test) within groups of 'dataframe'.
}
\usage{
pvalue (dataframe, control, sigval)
}

\arguments{
  \item{dataframe}{a summary dataframe of \code{\link{dfsummary}}: output}
  \item{control}{control group name}
    \item{sigval}{pvalue cutoff for significance}
}
\details{
The 'dataframe' should be having similar format of \code{\link{dfsummary}} output. 'control' should be an element from the column just before 'label'. 'N', 'Mean', 'SD' and 'CV' columns in the 'dataframe' are used for calculating p value by t-test (one to one t-test with 'control' in that group). significant if pvalue is < 'sigval'. Different groups in 'dataframe' are evaluated separately (columns before label is used for grouping).
}

\value{
A dataframe. New columns named 'pvalue' (p values of t-test.If the value is less than 0.001, then appear as "< 0.001") and 'significance' (yes if pvalue less than 'sigval') are attached to the left. }


\author{
A.A Palakkan
}

\examples{
## loading data
data(metafile384, rawdata384)
rawdata<-plate2df(data2plateformat(rawdata384,platetype = 384))
data_DF2<- dplyr::inner_join(rawdata,metafile384,by=c("row","col","position"))
result3 <- dfsummary(data_DF2,y = "value",
           grp_vector = c("cell","compound","concentration"),
           rm_vector = c("B", "drug2", "huh7"),
           nickname = "",
           rm = "FALSE", param = c(strict = "FALSE", cutoff = 40,n = 12))

## eg:1 t-test on result3.
pvalue(result3,"C3",sigval=0.05)



}

\keyword{ htest }
