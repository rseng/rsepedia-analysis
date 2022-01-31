---
title: Akmedoids R package for generating directionally-homogeneous clusters of longitudinal data sets
date: "09 January 2020"
bibliography: paper.bib
affiliations:
- index: 1
  name: Crime and Well-being Big Data Centre, Manchester Metropolitan University
authors:
- affiliation: 1
  name: Monsuru Adepeju
  orcid: 0000-0002-9006-4934
- affiliation: 1
  name: Samuel Langton
  orcid: 0000-0002-1322-1553
- affiliation: 1
  name: Jon Bannister
  orcid: 0000-0002-1350-510X
tags:
- Anchored k-medoids
- k-means
- crime
- longitudinal clustering
- long-term trends
---

# Statement of Need

In social and behavioural sciences, longitudinal clustering is widely used for identifying groups of individual trends that correspond to certain developmental processes over time. Whilst popular clustering techniques, such as k-means, are suited for identifying spherical clusters [@GenoFali:2010; @Curman:2015], there has been little attempt to modify such methods to identify alternative forms of cluster, such as those that represent linear growth over time (i.e. directionally-homogeneous clusters). To address this shortcoming, we introduce `Anchored k-medoids`, a package referred to as `Ak-medoids`, which implements a medoid-based expectation maximisation (MEM) procedure within a classical k-means clustering framework.  The package includes functions to assist in the manipulation of longitudinal data sets prior to the clustering procedure, and the visualisation of solutions post-procedure. The potential application areas of `Ak-medoids` include criminology, transport, epidemiology and brain imaging.

# Design and implementation

Previous studies have taken advantage of the various functional characteristics of longitudinal data in order to extract theoretically or empirically interesting clusters of subjects. Examples include using the Fourier basis [@Tarpey:2003] or the coefficients of the B-spline derivative estimates [@Boor:1978; @Schumaker:2007] which anchor clustering routines to better capture a presumed developmental process. Here, we develop an `Anchored k-medoids` (`Akmedoids`) clustering package, which employs the ordinary least square (OLS) trend lines of subjects, and a bespoke expectation-maximisation procedure, specifically to capture long-term linear growth. In criminology, identifying such slow-changing trends helps to unravel place-based characteristics that drive crime-related events, such as street gun and homicide, across a geographical space [@Griffith:2004]. To date, explorations of these trends have deployed existing techniques, namely k-means [@Curman:2015; @Andresen:2017] and group-based trajectory modelling [@Weisburd:2004; @Chavez:2009; @Bannister:2017], which are suited for spherical clusters [@GenoFali:2010]. The sensitivity of such techniques to short-term fluctuations and outliers in longitudinal datasets makes it more difficult to extract clusters based on the underlying long-term trends. `Akmdeoids` is tailored for such a scenario.
The main clustering function in the `Akmedoids` package  implements a medoid-based expectation maximisation (MEM) procedure by integrating certain key modifications into the classical k-means routine. First, it approximates longitudinal trajectories using OLS regression and second, anchors the initialisation process with medoid observations. It then deploys the medoid observations as new anchors for each iteration of the expectation-maximisation procedure [@Celeux:1992], until convergence. In a similar fashion to classical k-means, the routine relies on distance-based similarity between vectors of observations and is scale invariant. This implementation ensures that the impact of short-term fluctuations and outliers are minimised. The final groupings are augmented with the raw trajectories, and visualised, in order to provide a clearer delineation of the long-term linear trends of subject trajectories. Given an `l` number of iterations, the computational complexity of the clustering routine is the same as that of a classical k-means algorithm, i.e. `O(lkn)`, where `k` is the specified number of clusters and `n`, the number of individual trajectories. The optimal number of clusters for a given data may be determined using the average silhouette [@Rousseeuw:1987] or the Calinski and Harabasz criterion [@Calinski:1974]. A full demonstration is provided in the package vignette of how to deploy `Akmedoids` to examine long-term relative exposure to crime in `R` [@rManual:2020]. We encourage the use of the package outside of criminology.

# Clustering and cluster representations

The main clustering function of akmedoids is `akclustr`. The function captures directionally homogeneous clusters within any given longitudinal dataset using the procedure detailed above. For crime inequality studies, the package includes the `props` function for converting the absolute (or rate) measures of individual trajectories into a relative measures over time. The package includes the `print_akstats` and `plot_akstats` functions to generate the properties of clusters and visualize the clusters, respectively. In particular, the `plot_akstats` function draws from the `ggplot2` library [@ggplot:2016] in order to visualize the resulting clusters in either a line or an areal-stacked graph format.

# Acknowledgment

We gratefully acknowledge the Economic and Social Research Council (ESRC), who funded the Understanding Inequalities project (Grant Reference ES/P009301/1) through which this research was conducted.

# References

<!-- README.md is generated from README.Rmd. Please edit that file -->

akmedoids
=========

<!-- badges: start -->
<!-- badges: end -->

An R package for analyzing and clustering longitudinal data

### Description

The `akmedoids` package advances the clustering of longitudinal datasets
in order to identify clusters of trajectories with similar long-term
linear trends over time, providing an improved cluster identification as
compared with the classic kmeans algorithm. The package also includes a
set of functions for addressing common data issues, such as missing
entries and outliers, prior to conducting advance longitudinal data
analysis. One of the key objectives of this package is to facilitate
easy replication of a recent paper which examined small area inequality
in the crime drop (Adepeju et al. 2020). Many of the functions provided
in the `akmedoids` package may be applied to longitudinal data in
general.

**For more information and usability, check out details on
[CRAN](https://cran.r-project.org/web/packages/akmedoids/index.html).**

### Installation from `CRAN`

From an R console, type:

    #install.packages("akmedoids")
    library(akmedoids)

    #Other libraries
    library(tidyr)
    library(ggplot2)
    #> Warning: package 'ggplot2' was built under R version 4.0.3
    library(reshape)
    #> Warning: package 'reshape' was built under R version 4.0.3
    library(readr)

To install the development version of the package, type
`remotes::install_github("MAnalytics/akmedoids")`. Please, report any
installation problems in the issues.

### Example usage:

Given a longitudinal datasets, the following is an example of how
`akmedoids` could be used to extract clusters of trajectories with
similar long-term trends over time. We will use a simulated dataset
(named `simulated.rda`) stored in the `data/` directory.

### Generating artificial dataset

Simulating data set which comprised of three clusters with distinct mean
directional change over time. Each group contains 50 trajectories.

    dir.create("input") # create a folder
    #> Warning in dir.create("input"): 'input' already exists

    #function for creating longitudinal noise
    noise_fn <- function(x=3, time){
      rnorm(length(time), mean=0, x)}

    #function for simulating a trajectory group
    sim_group <- function(gr_baseline, sd, time){
      intcp_errors <- rgamma(1, shape=2, scale=sd) #intercept error
      mean_traj = gr_baseline + intcp_errors
      traj = mean_traj + noise_fn(intcp_errors, time)
    }

    #time steps
    t_steps <- c(0:20)

    #increasing group
    i_gr <- NULL
    for(i in seq_len(50)){
      i_gr <- rbind(i_gr,
                    sim_group(gr_baseline=(0.5*t_steps),
                              sd=1, time=t_steps))
    }

    #stable group
    s_gr <- NULL
    for(i in seq_len(50)){
      s_gr <- rbind(s_gr,
                    sim_group(gr_baseline=rep(3,length(t_steps)),
                              sd=1, time=t_steps))
    }

    #decreasing group
    d_gr <- NULL
    for(i in seq_len(50)){
      d_gr <- rbind(d_gr,
                    sim_group(gr_baseline=(10 - (0.5*t_steps)),
                              sd=1, time=t_steps))
    }

    #combine groups
    simulated <- data.frame(rbind(i_gr, s_gr, d_gr))

    #add group label
    simulated <- data.frame(cbind(ID=1:nrow(simulated), simulated))

    colnames(simulated) <- c("ID", 1:(ncol(simulated)-1))

    #save data set
    ##simulated = readr::write_csv(simulated, "input/example-simulated.csv")

### Visualising artificial dataset


    #import already save simulated data
    Import_simulated = read_csv(file="./input/example-simulated.csv")

    #preview the data
    head(Import_simulated)
    #> # A tibble: 6 x 22
    #>      ID    `1`   `2`    `3`    `4`    `5`   `6`   `7`   `8`   `9`  `10`  `11`
    #>   <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    #> 1     1 -0.759  1.90  3.02  -0.355  1.63   4.24  4.52  7.93  5.18  6.84  6.09
    #> 2     2  5.52   3.22  4.49   4.61  11.2    5.90 15.7   7.76  7.46  9.76  6.61
    #> 3     3  1.52   1.09  2.05   1.89   0.728  1.27  4.76  3.75  4.58  5.81  5.39
    #> 4     4  3.24   2.41 -0.169  5.59   1.29   3.42  4.59  7.29  5.23  7.25  6.90
    #> 5     5  0.950  2.58  2.58   1.39   3.22   3.28  3.40  6.21  5.71  3.69  7.39
    #> 6     6  2.33   7.88  2.15   1.89   7.38   3.57  7.33 10.4   5.90  5.12 10.3 
    #> # ... with 10 more variables: `12` <dbl>, `13` <dbl>, `14` <dbl>, `15` <dbl>,
    #> #   `16` <dbl>, `17` <dbl>, `18` <dbl>, `19` <dbl>, `20` <dbl>, `21` <dbl>


    #convert wide-format into long
    simulated_long <- melt(t(Import_simulated), id.vars=c("ID"))

     simulated_long <- simulated_long %>%
      dplyr::filter(X1!="ID") %>%
      dplyr::rename(Time=X1, ID=X2)
    # 
    simulated_long <- data.frame(cbind(simulated_long,
                                       Groups= c(rep("Increasing", 50*21),
                                       rep("Stable", 50*21),
                                       rep("Decreasing", 50*21))))

    #re-order levels
    simulated_long$Time <- factor(simulated_long$Time, 
                                    levels = c(1:21))

    ggplot(simulated_long, aes(x = Time, y = value, group=ID, color=Groups)) +
      geom_point(size=0.5) + 
      geom_line() +
      scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
      theme_light()

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

### Performing clustering using `akmedoids`

Performing clustering analysis using `akmedoids` package.


    output <- akclustr(Import_simulated, id_field=TRUE, verbose = FALSE, k=c(3,12), crit = "Silhouette",
                      quality_plot=TRUE)
    #> `geom_smooth()` using formula 'y ~ x'

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />


    #Remark: The entire output can be printed out to the console by typing `output`

    #To check the optimal cluster size, type:
    output$optimal_k
    #> [1] 3.000046

    #Show quality plot
    #plot(output$QltyPlot)

    #user may decide the examine the plot post clustering....

### Documentation

From an R console type `??akmedoids` for help on the package. The
package page on CRAN is
[here](https://cran.r-project.org/web/packages/akmedoids/index.html),
package reference manual is
[here](https://cran.r-project.org/web/packages/akmedoids/akmedoids.pdf),
package vignette is
[here](https://cran.r-project.org/web/packages/akmedoids/vignettes/akmedoids-vignette.html).

### Support and Contributions:

For support and bug reports send an email to:
<a href="mailto:monsuur2010@yahoo.com" class="email">monsuur2010@yahoo.com</a>
or open an issue [here](https://github.com/MAnalytics/akmedoids/issues).
Code contributions to akmedoids are also very welcome.

### References:

Rousseeuw, P. J. 1987. “Silhouettes: A Graphical Aid to the
Interpretation and Validation of Cluster Analysis.” Journal of
Computational and Applied Mathematics, no. 20: 53–6.
[link](https://www.bibsonomy.org/bibtex/bc0f62c7895f91c787354d03f23da976)

Caliński, T., and J. Harabasz. 1974. “A Dendrite Method for Cluster
Analysis.” Communications in Statistics-Theory and Methods, 3(1): 1–27.
[link](https://www.tandfonline.com/doi/abs/10.1080/03610927408827101)

Adepeju, M., Langton, S. and Bannister, J. 2020. Anchored k-medoids: a
novel adaptation of k-means further refined to measure instability in
the exposure to crime. Journal of Computational Social Science,
(revised).
---
title: "NEWS.md"
authors: "geoMADE"
date: "12th April 2021"
output: html_document
---


Akmedoids' package updated (Version: v1.3.0)

Updates:
1. Added 4 new sample datasets, and used in the function examples. The datasets include `clustr`, `simulated`, `TO1Risk`, `traj_w_space`.
2. Added detailed descriptions of fields in each dataset. See `format` section of each dataset documentation.
3. Change the citation `Adepeju et al. 2019` to `Adepeju et al. 2021` based on the newly published article upon which the `akmedoids` package is based.
4. Modified the names of some functions, e.g. `akmedoids.clust` changed to `akclustr`,  `statPrint` changed to `print_akstats`, `population` changed to `popl`.
5. Added two new functions, namely; `remove_rows_n`, `plot_akstats` (see the documentation for details)
6. Ensured that the names of functions are in lower cases, e.g. `dataImputation` changed to `data_imputation`, `outlierDetect` changed to `outlier_detect`, etc.
7. Modify example in 'elbow_point' function to reduce run time


Your faithfully.
Monsuru.


## Test environments
* local OS Windows, R 3.5.3
* Debian Linux (R-devel), R 3.5.3
* win-builder (devel and release), R 3.5.3 


## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

##

'Akmedoids' package updated 30th March 2021 (New version: v1.3.0)

### Updates:

'Akmedoids' package updated (Version: v1.3.0)

### Updates:
1. Added 4 new sample datasets, and used in the function examples. The datasets include `clustr`, `simulated`, `TO1Risk`, `traj_w_space`.
2. Added detailed descriptions of fields in each dataset. See `format` section of each dataset documentation.
3. Change the citation `Adepeju et al. 2019` to `Adepeju et al. 2021` based on the newly published article upon which the `akmedoids` package is based.
4. Modified the names of some functions, e.g. `akmedoids.clust` changed to `akclustr`,  `statPrint` changed to `print_akstats`, `population` changed to `popl`.
5. Added two new functions, namely; `remove_rows_n`, `plot_akstats` (see the documentation for details)
6. Ensured that the names of functions are in lower cases, e.g. `dataImputation` changed to `data_imputation`, `outlierDetect` changed to `outlier_detect`, etc.
7. Modify example in 'elbow_point' function to reduce run time

Your faithfully.
Monsuru.


---
#output: github_document
output: md_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# akmedoids

<!-- badges: start -->
<!-- badges: end -->

An R package for analyzing and clustering longitudinal data

### Description

The `akmedoids` package advances the clustering of longitudinal datasets in order to identify clusters of trajectories with similar long-term linear trends over time, providing an improved cluster identification as compared with the classic kmeans algorithm. The package also includes a set of functions for addressing common data issues, such as missing entries and outliers, prior to conducting advance longitudinal data analysis. One of the key objectives of this package is to facilitate easy replication of a recent paper which examined small area inequality in the crime drop (Adepeju et al. 2020). Many of the functions provided in the `akmedoids` package may be applied to longitudinal data in general. 

**For more information and usability, check out details on [CRAN](https://cran.r-project.org/web/packages/akmedoids/index.html).**


### Installation from `CRAN`

From an R console, type:

```{r, echo=TRUE, message=FALSE, eval=TRUE}
#install.packages("akmedoids")
library(akmedoids)

#Other libraries
library(tidyr)
library(ggplot2)
library(reshape)
library(readr)
```
To install the development version of the package, type `remotes::install_github("MAnalytics/akmedoids")`. Please, report any installation problems in the issues.

### Example usage:

Given a longitudinal datasets, the following is an example of how `akmedoids` could be used to extract clusters of trajectories with similar long-term trends over time. We will use a simulated dataset (named `simulated.rda`) stored in the `data/` directory.


### Generating artificial dataset

Simulating data set which comprised of three clusters with distinct mean directional change over time. Each group contains 50 trajectories.


```{r, echo=TRUE, message=FALSE, eval=TRUE}
dir.create("input") # create a folder

#function for creating longitudinal noise
noise_fn <- function(x=3, time){
  rnorm(length(time), mean=0, x)}

#function for simulating a trajectory group
sim_group <- function(gr_baseline, sd, time){
  intcp_errors <- rgamma(1, shape=2, scale=sd) #intercept error
  mean_traj = gr_baseline + intcp_errors
  traj = mean_traj + noise_fn(intcp_errors, time)
}

#time steps
t_steps <- c(0:20)

#increasing group
i_gr <- NULL
for(i in seq_len(50)){
  i_gr <- rbind(i_gr,
                sim_group(gr_baseline=(0.5*t_steps),
                          sd=1, time=t_steps))
}

#stable group
s_gr <- NULL
for(i in seq_len(50)){
  s_gr <- rbind(s_gr,
                sim_group(gr_baseline=rep(3,length(t_steps)),
                          sd=1, time=t_steps))
}

#decreasing group
d_gr <- NULL
for(i in seq_len(50)){
  d_gr <- rbind(d_gr,
                sim_group(gr_baseline=(10 - (0.5*t_steps)),
                          sd=1, time=t_steps))
}

#combine groups
simulated <- data.frame(rbind(i_gr, s_gr, d_gr))

#add group label
simulated <- data.frame(cbind(ID=1:nrow(simulated), simulated))

colnames(simulated) <- c("ID", 1:(ncol(simulated)-1))

#save data set
##simulated = readr::write_csv(simulated, "input/example-simulated.csv")

```


### Visualising artificial dataset

```{r, echo=TRUE, message=FALSE, eval=TRUE}

#import already save simulated data
Import_simulated = read_csv(file="./input/example-simulated.csv")

#preview the data
head(Import_simulated)


#convert wide-format into long
simulated_long <- melt(t(Import_simulated), id.vars=c("ID"))

 simulated_long <- simulated_long %>%
  dplyr::filter(X1!="ID") %>%
  dplyr::rename(Time=X1, ID=X2)
# 
simulated_long <- data.frame(cbind(simulated_long,
                                   Groups= c(rep("Increasing", 50*21),
                                   rep("Stable", 50*21),
                                   rep("Decreasing", 50*21))))

#re-order levels
simulated_long$Time <- factor(simulated_long$Time, 
                                levels = c(1:21))

ggplot(simulated_long, aes(x = Time, y = value, group=ID, color=Groups)) +
  geom_point(size=0.5) + 
  geom_line() +
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  theme_light()

```

### Performing clustering using `akmedoids`

Performing clustering analysis using `akmedoids` package.

```{r, echo=TRUE, message=TRUE, eval=TRUE}

output <- akclustr(Import_simulated, id_field=TRUE, verbose = FALSE, k=c(3,12), crit = "Silhouette",
                  quality_plot=TRUE)

#Remark: The entire output can be printed out to the console by typing `output`

#To check the optimal cluster size, type:
output$optimal_k

#Show quality plot
#plot(output$QltyPlot)

#user may decide the examine the plot post clustering....
```

### Documentation

From an R console type `??akmedoids` for help on the package. The package page on CRAN is [here](https://cran.r-project.org/web/packages/akmedoids/index.html), package reference manual is [here](https://cran.r-project.org/web/packages/akmedoids/akmedoids.pdf), package vignette is [here](https://cran.r-project.org/web/packages/akmedoids/vignettes/akmedoids-vignette.html). 

### Support and Contributions:

For support and bug reports send an email to: monsuur2010@yahoo.com or open an issue [here](https://github.com/MAnalytics/akmedoids/issues). Code contributions to akmedoids are also very welcome.

### References:

Rousseeuw, P. J. 1987. “Silhouettes: A Graphical Aid to the Interpretation and Validation of Cluster Analysis.” Journal of Computational and Applied Mathematics, no. 20: 53–6. [link](https://www.bibsonomy.org/bibtex/bc0f62c7895f91c787354d03f23da976)

Caliński, T., and J. Harabasz. 1974. “A Dendrite Method for Cluster Analysis.” Communications in Statistics-Theory and Methods, 3(1): 1–27. [link](https://www.tandfonline.com/doi/abs/10.1080/03610927408827101)

Adepeju, M., Langton, S. and Bannister, J. 2020. Anchored k-medoids: a novel adaptation of k-means further refined to measure instability in the exposure to crime. Journal of Computational Social Science, (revised).

---
title: "A guide to measuring long-term inequality in the exposure to crime at micro-area levels using `'Akmedoids'` package"

author: |
  | `Authors:`
  | `Adepeju, M., Langton, S., and Bannister, J.`
  | `Big Data Centre, Manchester Metropolitan University, Manchester, M15 6BH`
  
date: |
  | `Date:`
  | ``r Sys.Date()``

output:
  rmarkdown::html_vignette
  
#dev: png
#output:
  #word_document: default
  #always_allow_html: yes
#  pdf_document: default
always_allow_html: yes
#fig_caption: yes
bibliography: references.bib

abstract: The `'akmedoids'` package advances a set of R-functions for longitudinal clustering of long-term trajectories and determines the optimal solution based on the `Caliński-Harabasz` criterion (Caliński and Harabasz 1974). The package also includes a set of functions for addressing common data issues, such as missing entries and outliers, prior to conducting advance longitudinal data analysis. One of the key objectives of this package is to facilitate easy replication of a recent paper which examined small area inequality in the crime drop (see Adepeju et al. 2021). This document is created to provide a guide towards accomplishing this objective. Many of the functions provided in the `akmedoids` package may be applied to longitudinal data in general.
  
vignette: >
  %\VignetteIndexEntry{A guide to measuring long-term inequality in the exposure to crime at micro-area levels using 'Akmedoids' package}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

<style type="text/css">

h1.title {
  font-size: 26px;
  line-height: 130%;
  color: Black;
  text-align: center;
}

h2.subtitle {
  font-size: 13px;
  line-height: 120%;
  color: Black;
  text-align: center;
}

h4.author { /* Header 4 - and the author and data headers use this too  */
  font-size: 17px;
  font-family: "Arial";
  color: Black;
  text-align: center;
}
h4.date { /* Header 4 - and the author and data headers use this too  */
  font-size: 17px;
  font-family: "Arial", Times, serif;
  color: Black;
  text-align: center;
}

h4.abstract { /* Header 4 - and the author and data headers use this too  */
  font-size: 10px;
  font-family: "Arial", Times, serif;
  color: black;
  text-align: center;
}

h4.institute{ /* Header 4 - and the author and data headers use this too  */
  font-size: 10px;
  font-family: "Arial", Times, serif;
  color: black;
  text-align: center;
}

body, td {
   font-size: 14px;
}
code.r{
  font-size: 13px;
}
pre {
  font-size: 13px
}
h1 { /* Header 1 */
  font-size: 16px;
  color: DarkBlue;
}
h2 { /* Header 2 */
    font-size: 16px;
  color: DarkBlue;
}
h3 { /* Header 3 */
  font-size: 15px;
  font-family: "Times New Roman", Times, serif;
  color: DarkBlue;

</style>

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r functions, include=FALSE}
# A function for captioning and referencing images
fig <- local({
    i <- 0
    ref <- list()
    list(
        cap=function(refName, text) {
            i <<- i + 1
            ref[[refName]] <<- i
            paste("Figure ", i, ": ", text, sep="")
        },
        ref=function(refName) {
            ref[[refName]]
        })
})
```



# Introduction

Longitudinal clustering analysis is ubiquitous in social and behavioral sciences for investigating the developmental processes of a phenomenon over time. Examples of the commonly used techniques in these areas include group-based trajectory modeling (GBTM) and the non-parametric kmeans method. Whilst kmeans has a number of benefits over GBTM, such as more relaxed statistical assumptions, generic implementations render it more sensitive to outliers and short-term fluctuations, which minimises its ability to identify long-term linear trends in data. In crime and place research, for example, the identification of such `long-term` linear trends may help to develop some theoretical understanding of criminal victimization within a geographical space [@Weisburd2004; @Griffith2004]. In order to address this sensitivity problem, we advance a novel technique named `anchored kmedoids` (`'akmedoids'`) which implements three key modifications to the existing longitudinal `kmeans` approach. First, it approximates trajectories using ordinary least square regression (`OLS`) and second, `anchors` the initialisation process with median observations. It then deploys the `medoids` observations as new anchors for each iteration of the expectation-maximization procedure [@Celeux1992]. These modifications ensure that the impacts of short-term fluctuations and outliers are minimized. By linking the final groupings back to the original trajectories, a clearer delineation of the long-term linear trends of trajectories are obtained.

We facilitate the easy use of `akmedoids` through an open-source package using `R`. We encourage the use of the package outside of criminology, should it be appropriate. Before outlining the main `clustering` functions, we demonstrate the use of a few `data manipulation` functions that assist in data preparation. The worked demonstration uses a small example dataset which should allow users to get a clear understanding of the operation of each function. 

```{r, eval=TRUE, echo=FALSE, include=FALSE}
#install.packages("kableExtra")
require(knitr)
library(gdtools)
library(kableExtra)
library(clusterCrit)
library(dplyr)
```

```{r, echo=FALSE, include=FALSE}
col1 <- c("1", "2","3","4", "5", "6")
col2 <- c("`data_imputation`","`rates`", "`props`", "`outlier_detect`","`w_spaces`", "`remove_rows_n`")
col3 <- c("Data imputation for longitudinal data", "Conversion of 'counts' to 'rates'", "Conversion of 'counts' (or 'rates') to 'Proportion'", "Outlier detection and replacement","Whitespace removal", "Incomplete rows removal")
col4 <- c("Calculates any missing entries (`NA`, `Inf`, `null`) in a longitudinal data, according to a specified method","Calculates rates from observed 'counts' and its associated denominator data", "Converts 'counts' or 'rates' observation to 'proportion'", "Identifies outlier observations in the data, and replace or remove them","Removes all the leading and trailing whitespaces in a longitudinal data", "Removes rows which contain 'NA' and 'inf' entries")
tble <- data.frame(col1, col2, col3, col4)
tble <- tble
```

```{r table1, results='asis', echo=FALSE, tidy.opts=list(width.cutoff=50)}
knitr::kable(tble, caption = "Table 1. `Data manipulation` functions", col.names = c("SN","Function","Title","Description")) %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T, border_right = T) %>%
  column_spec(2, width = "8em", background = "white") %>%
  column_spec(3, width = "12em", background = "white") %>%
  column_spec(4, width = "16em", background = "white")#%>%
  #row_spec(3:5, bold = T, color = "white", background = "#D7261E")
```

# 1. Data manipulation

Table 1 shows the main data manipulation functions and their descriptions. These functions help to address common data issues prior to analysis, as well as basic data manipulation tasks such as converting longitudinal data from `count` to `proportion` measures (as per the crime inequality paper where `akmedoids` was first implemented). In order to demonstrate the utility of these functions, we provide a simulated dataset `traj` which can be called by typing `traj` in `R` console after loading the `akmedoids` library.


## (i) `"data_imputation"` function

Calculates any missing entries in a data, according to a chosen method. This function recognizes three kinds of data entries as missing. These are `NA`, `Inf`, `null`, and an option of whether or not to consider `0`'s as missing values. The function provides a replacement option for the missing entries using two methods. First, an `arithmetic` method which uses the `mean`, `minimum` or `maximum` value of the corresponding rows or columns of the missing values. Second, a `regression` method which uses OLS  regression lines to estimate the missing values. Using the regression method, only the missing data points derive values from the regression line while the remaining (observed) data points retain their original values. The function terminates if there is any trajectory with only one observation in it. Using the `'traj'` dataset, we demonstrate how the `'regression'` method estimates missing values. 


```{r, eval=FALSE}
#installing the `akmedoids` packages
install.packages("devtools")
devtools::install_github("manalytics/akmedoids")

```

```{r, eval=TRUE}

#loading the package
library(akmedoids)
```

```{r, eval=TRUE}

#import and preview the first 6 rows of 'traj' object
data(traj)

head(traj)

#no. of rows
nrow(traj) 

#no. of columns
ncol(traj) 

```

The first column of the `traj` object is the `id` (unique) field. In many applications, it is necessary to preserve the `id` column in order to allow linking of outputs to other external datasets, such as spatial location data. Most of the functions of the `akmedoids` provides an option to recognise the first column of an input dataset as the unique field. The `data_imputation` function can be used to imput the missing data point of `traj` object as follows:


```{r, eval=TRUE}
imp_traj <- data_imputation(traj, id_field = TRUE, method = 2, 
               replace_with = 1, fill_zeros = FALSE)

imp_traj <- imp_traj$CompleteData
#viewing the first 6 rows
head(imp_traj)
```

```{r figs1, echo=FALSE, fig.width=5,fig.height=6,fig.align="center", fig.cap=fig$cap("figs1", "Data imputation with regression")}

par(mar=c(2,2,2,2)+0.1)
par(adj = 0)
par(mfrow=c(6,2))
dev.new()
dat <- as.data.frame(traj)
t_name <- as.vector(traj[,1])
dat <- dat[,2:ncol(dat)]
#if(k==nrow(dat)){
  #}
#head(dat)
for(k in seq_len(nrow(dat))){ #k<-2
  y <- suppressWarnings(as.numeric(as.character(dat[k,])))
  x <- seq_len(length(y))
  known <- data.frame(x, y)
  known_1 <- data.frame(known[is.na(known[,2])|is.infinite(known[,2]),])  #
  known_2 <- data.frame(known[!is.na(known[,2])&!is.infinite(known[,2]),])
  #train the available data using linear regression
  model.lm <- lm(y ~ x, data = known_2)
  # Use predict the y value for the removed data
  newY <- predict(model.lm, newdata = data.frame(x = known_1[,1]))
   l_pred <- predict(model.lm, newdata = data.frame(1:9)) #line
  #add to the original data.
  dat[k, known_1[,1]] <- newY
  #Add the predicted points to the original data
  #dev.new()
  #plot(1:10, col=2)
  #options(rgl.useNULL = TRUE)
  plot (known$x, known$y, type="o", main=paste("traj_id:",t_name[k], sep=" "), font.main = 1)
  if(!length(newY)==0){#plot only if it has elements
  lines(l_pred, lty="dotted", col="red", lwd=2)
  }
  points(known_1[,1], newY, col = "red")
}
#point legend
plot_colors <- c("black","red")
text <- c("Observed points", "Predicted points")
#options(rgl.useNULL = TRUE)
par(xpd=TRUE)
legend("center",legend = text, text.width = max(sapply(text, strwidth)),
       col=plot_colors, pch = 1, cex=1, horiz = FALSE)
par(xpd=FALSE)

#line legend
plot_colors <- c("black","red")
text <- c("line joining observed points", "regression line predicting missing points")
#options(rgl.useNULL = TRUE)
plot.new()
par(xpd=TRUE)
legend("center",legend = text, text.width = max(sapply(text, strwidth)),
       col=plot_colors, lwd=1, cex=1, lty=c(1,2), horiz = FALSE)
par(xpd=FALSE)
```


The argument `method = 2` refers to the `regression` technique, while the argument `replace_with = 1` refers to the `linear` option (currently the only available option). Figure `r fig$ref("figs1")` is a graphical illustration of how this method approximates the missing values of the `traj` object.

### Estimating the population data using the `'data_imputation'` function

Obtaining the denominator information (e.g. population estimates to normalize counts) of local areas within a city for non-census years is problematic in longitudinal studies. This challenge poses a significant drawback to the accurate estimation of various measures, such as crime rates and population-at-risk of an infectious disease. Assuming a limited amount of denominator information is available, an alternative way of obtaining the missing data points is to interpolate and/or extrapolate the missing population information using the available data points. The `data_imputation` function can be used to perform this task. 

The key step towards using the function for this purpose is to create a matrix, containing both the available fields and the missing fields arranged in their appropriate order. All the entries of the missing fields can be filled with either `NA` or `null`. Below is a demonstration of this task with a sample population dataset with only two available data fields. The corresponding `input` matrix is constructed as shown.

```{r, eval=TRUE}

#import population data
data(popl)

#preview the data
head(popl)

nrow(popl) #no. of rows

ncol(popl) #no. of columns
```


The corresponding `input` dataset is prepared as follows and saved as `population2`:


```{r, echo=FALSE}
#create a matrix of the same rows and column as the `traj` data
pop <- as.data.frame(matrix(0, nrow(popl), ncol(traj)))
colnames(pop) <- names(traj) 
pop[,1] <- as.vector(as.character(popl[,1]))
pop[,4] <- as.vector(as.character(popl[,2]))
pop[,8] <- as.vector(as.character(popl[,3]))
list_ <- c(2, 3, 5, 6, 7, 9, 10)
for(u_ in seq_len(length(list_))){ #u_<-1
  pop[,list_[u_]] <- "NA"
}

head(pop)

population2 <- pop
```


The missing values are estimated as follows using the `regression` method of the `data_imputation` function: 


```{r, eval=TRUE}

pop_imp_result <- data_imputation(population2, id_field = TRUE, method = 2, 
               replace_with = 1, fill_zeros = FALSE)

pop_imp_result <- pop_imp_result$CompleteData

#viewing the first 6 rows
head(pop_imp_result)

```


Given that there are only two data points in each row, the `regression` method will simply generate the missing values by fitting a straight line to the available data points. The higher the number of available data points in any trajectory the better the estimation of the missing points. Figure `r fig$ref("figs1")` illustrates this estimation process. 


## (ii) `"rates"` function

Given a longitudinal data ($m\times n$) and its associated denominator data ($s\times n$), the `'rates'` function converts the longitudinal data to 'rates' measures (e.g. counts per 100 residents). Both the longitudinal and the denominator data may contain different number of rows, but need to have the same number of columns, and must include the `id` (unique) field as their first column. The rows do not have to be sorted in any particular order. The rate measures (i.e. the output) will contain only rows whose `id's` match from both datasets. We demonstrate the utility of this function using the `imp_traj` object (above) and the estimated population data ('`pop_imp_result`'). 


```{r, eval=TRUE}

#example of estimation of 'crimes per 200 residents'
crime_per_200_people <- rates(imp_traj, denomin=pop_imp_result, id_field=TRUE, 
                              multiplier = 200)

#view the full output
crime_per_200_people <- crime_per_200_people$rates_estimates

#check the number of rows
nrow(crime_per_200_people)

```

From the output, it can be observed that the number of rows of the output data is 9. This implies that only 9 `location_ids` match between the two datasets. The unmatched `ids` are ignored. **Note**: the calculation of `rates` often returns outputs with some of the cell entries having `Inf` and `NA` values, due to calculation errors and character values in the data. We therefore recommend that users re-run the `data_imputation` function after generating `rates` measures, especially for a large data matrix.



## (iii) `"props"` function


Given a longitudinal data, the `props` function converts each data point (i.e. entry in each cell) to the proportion of the sum of their corresponding column. Using the `crime_per_200_people` estimated above, we can derive the `proportion of crime per 200 people` for each entry as follows:  


```{r, eval=TRUE}

#Proportions of crimes per 200 residents
prop_crime_per200_people <- props(crime_per_200_people, id_field = TRUE, scale = 1, digits=2)

#view the full output
prop_crime_per200_people


#A quick check that sum of each column of proportion measures adds up to 1.  
colSums(prop_crime_per200_people[,2:ncol(prop_crime_per200_people)])


```

In line with the demonstration in Adepeju et al. (2021), we will use these `proportion` measures to demonstrate the main clustering function of this package. 


## (iv) `"outlier_detect"` function
This function is aimed at allowing users to identify any outlier observations in their longitudinal data, and replace or remove them accordingly. The first step towards identifying outliers in any data is to visualize the data. A user can then decide a cut-off value for isolating the outliers. The `outlier_detect` function provides two options for doing this: (`i`) a `quantile` method, which isolates any observations with values higher than a specified quantile of the data values distribution, and (`ii`) a `manual` method, in which a user specifies the cut-off value. The '`replace_with`' argument is used to determine whether an outlier value should be replaced with the mean value of the row or the mean value of the column in which the outlier is located. The user also has the option to simply remove the trajectory that contains an outlier value. In deciding whether a trajectory contains outlier or not, the `count` argument allows the user to set an horizontal threshold (i.e. number of outlier values that must occur in a trajectory) in order for the trajectory to be considered as having outlier observations. Below, we demonstrate the utility of the `outlier_detect` function using the `imp_traj` data above. 

## (v) `"w_spaces"` function
Given a matrix suspected to contain whitespaces, this function removes all the whitespaces and returns a cleaned data. ’Whitespaces’ are white characters often introduced during data entry, for instance by wrongly pressing the spacebar. For example, neither " A" nor "A " equates "A" because of the whitespaces that exist in them. They can also result from systematic errors in data recording devices.

## (vi) `"remove_rows_n"` function
This function removes any rows in which an 'NA' or an 'Inf' entry is found.


```{r figs2, echo=TRUE, fig.width=6,fig.height=3,fig.align="center", fig.cap=fig$cap("figs2", "Identifying outliers")}

#Plotting the data using ggplot library
library(ggplot2)
#library(reshape2)

#converting the wide data format into stacked format for plotting
#doing it manually instead of using 'melt' function from 'reshape2'

#imp_traj_long <- melt(imp_traj, id="location_ids") 

coln <- colnames(imp_traj)[2:length(colnames(imp_traj))]
code_ <- rep(imp_traj$location_ids, ncol(imp_traj)-1)
d_bind <- NULL
  for(v in seq_len(ncol(imp_traj)-1)){
    d_bind <- c(d_bind, as.numeric(imp_traj[,(v+1)]))
  }

code <- data.frame(location_ids=as.character(code_))
variable <- data.frame(variable=as.character(rep(coln,
                        each=length(imp_traj$location_ids))))
value=data.frame(value = as.numeric(d_bind))

imp_traj_long <- bind_cols(code, variable,value) 
  
#view the first 6 rows
head(imp_traj_long)

#plot function
p <-  ggplot(imp_traj_long, aes(x=variable, y=value,
            group=location_ids, color=location_ids)) + 
            geom_point() + 
            geom_line()

#options(rgl.useNULL = TRUE)
print(p)

```

Figure `r fig$ref("figs2")`  is the output of the above plot function.  

Based on Figure `r fig$ref("figs2")` if we assume that observations of `x2001`, `x2007` and `x2008` of trajectory id `E01004806` are outliers, we can set the `threshold` argument as `20`. In this case, setting `count=1` will suffice as the trajectory is clearly separable from the rest of the trajectories.  


```{r figs3, echo=TRUE, fig.width=6,fig.height=3,fig.align="center", fig.cap=fig$cap("figs3", "Replacing outliers with mean observation")}

imp_traj_New <- outlier_detect(imp_traj, id_field = TRUE, method = 2, 
                              threshold = 20, count = 1, replace_with = 2)

imp_traj_New <- imp_traj_New$Outliers_Replaced 

#options(rgl.useNULL = TRUE)
print(imp_traj_New)

#imp_traj_New_long <- melt(imp_traj_New, id="location_ids") 

coln <- colnames(imp_traj_New)[2:length(colnames(imp_traj_New))]
code_ <- rep(imp_traj_New$location_ids, ncol(imp_traj_New)-1)

d_bind <- NULL
  for(v in seq_len(ncol(imp_traj_New)-1)){
    d_bind <- c(d_bind, as.numeric(imp_traj_New[,(v+1)]))
  }

code <- data.frame(location_ids=as.character(code_))
variable <- data.frame(variable=as.character(rep(coln,
              each=length(imp_traj_New$location_ids))))
value=data.frame(value = as.numeric(d_bind))

imp_traj_New_long <- bind_cols(code, variable,value)

#plot function
#options(rgl.useNULL = TRUE)
p <-  ggplot(imp_traj_New_long, aes(x=variable, y=value,
            group=location_ids, color=location_ids)) + 
            geom_point() + 
            geom_line()

#options(rgl.useNULL = TRUE)
print(p)

```

Setting `replace_with = 2`, that is to replace the outlier points with the 'mean of the row observations', the function generates outputs re-plotted in Figure `r fig$ref("figs3")`.

## (vii) 'Other' functions

Please see the `akmedoids` user manual for the remaining `data manipulation` functions. 


# 2. Data Clustering

Table 2 shows the two main functions required to carry out the longitudinal clustering and generate the descriptive statistics of the resulting groups. The relevant functions are `akclustr` and `print_akstats`. The `akclustr` function clusters trajectories according to the similarities of their long-term trends, while the `print_akstats` function extracts descriptive and change statistics for each of the clusters. The former also generates `quality` plots for the best cluster solution.

The long-term trends of trajectories are defined in terms of a set of OLS regression lines. This allows the clustering function to classify the final groupings in terms of their slopes as `rising`, `stable`, and `falling`. The key benefits of this implementation is that it allows the clustering process to ignore the short-term fluctuations of actual trajectories and focus on their long-term linear trends. Adepeju and colleagues (2021) applied this technique in crime concentration research for measuring long-term inequalities in the exposure to crime at find-grained spatial scales. 

```{r figs4, echo=FALSE, fig.cap=fig$cap("figs4", paste("Long-time linear trends of relative (`proportion`, `p`) crime exposure. Three inequality trends: trajectory i1: crime exposure is falling faster, i2, crime exposure is falling at the same rate, and i3, crime exposure is falling slower or increasing, relatively to the citywide trend. (Source:", "Adepeju et al. 2021)", sep=" ")), out.width = '60%', fig.align="center"} 
knitr::include_graphics("inequality.png")
```

Their implementation was informed by the conceptual (`inequality`) framework shown in Figure `r fig$ref("figs4")`. That said, `akmedoids` can be deployed on any measure (counts, rates) and is not limited to criminology, but rather, any field where the aim is to cluster longitudinal data based on long-term trajectories. By mapping the resulting trend lines grouping to the original trajectories, various performance statistics can be generated. 

In addition to the use of trend lines, the `akmedoids` makes two other modifications to the expectation-maximisation clustering routines [@Celeux1992]. First, the `akmedoids` implements an anchored median-based initialisation strategy for the clustering to begin. The purpose behind this step is to give the algorithm a theoretically-driven starting point and try and ensure that heterogenous trend slopes end up in different clusters (@Khan2004; @Steinley2007). Second, instead of recomputing centroids based on the mean distances between each trajectory trend lines and the cluster centers, the median of each cluster is selected and then used as the next centroid. This then becomes the new anchor for the current iteration of the expectation-maximisation step [@Celeux1992]. This strategy is implemented in order to minimize the impact of outliers. The iteration then continues until an objective function is maximised.  


```{r, echo=FALSE, include=FALSE}

col1 <- c("1", "2", "3")
col2 <- c("`akclustr`","`print_akstats`", "`plot_akstats`")
col3 <- c("`Anchored k-medoids clustering`","`Descriptive (Change) statistics of clusters`", "`Plots of cluster groups`")
col4 <- c("Clusters trajectories into a `k` number of groups according to the similarities in their long-term trend and determines the best solution based on the Silhouette width measure or the Calinski-Harabasz criterion","Generates the descriptive and change statistics of groups, and also plots the groups performances", "Generates different plots of cluster groups")
tble2 <- data.frame(col1, col2, col3, col4)
tble2 <- tble2

```

```{r table2, results='asis', echo=FALSE, tidy.opts=list(width.cutoff=50)}
knitr::kable(tble2, caption = "Table 2. `Data clustering` functions", col.names = c("SN","Function","Title","Description")) %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T, border_right = T) %>%
  column_spec(2, width = "8em", background = "white") %>%
  column_spec(3, width = "12em", background = "white") %>%
  column_spec(4, width = "16em", background = "white")#%>%
  #row_spec(3:5, bold = T, color = "white", background = "#D7261E")

```


```{r figs5, echo=TRUE, fig.width=6,fig.height=3,fig.align="center", fig.cap=fig$cap("figs5",  "Trajectory of crime proportions over time")}

#Visualizing the proportion data

#view the first few rows
head(prop_crime_per200_people)

#prop_crime_per200_people_melt <- melt(prop_crime_per200_people, id="location_ids") 
coln <- colnames(prop_crime_per200_people)[2:length(colnames(prop_crime_per200_people))]
code_ <- rep(prop_crime_per200_people$location_ids, ncol(prop_crime_per200_people)-1)
d_bind <- NULL
  for(v in seq_len(ncol(prop_crime_per200_people)-1)){
    d_bind <- c(d_bind, prop_crime_per200_people[,(v+1)])
  }

prop_crime_per200_people_melt <- data.frame(cbind(location_ids=as.character(code_), variable =
                        rep(coln,
                        each=length(prop_crime_per200_people$location_ids)), value=d_bind))

#plot function
#options(rgl.useNULL = TRUE)
p <-  ggplot(prop_crime_per200_people_melt, aes(x=variable, y=value,
            group=location_ids, color=location_ids)) + 
            geom_point() + 
            geom_line()

#options(rgl.useNULL = TRUE)
print(p)

```

In the following sections, we provide a worked example of clustering with `akclustr` function using the `prop_crime_per200_people` object. The function will generate cluster solution over a set of `k` values, determine the optimal value of `k`. The `print_akstats` function will then be applied to generate the descriptive summary and the change statistics of the clusters. The `prop_crime_per200_people` object is plotted in `r fig$ref("figs5")`.  


## (i) `akclustr` function

***Dataset***: 

Each trajectory in Figure `r fig$ref("figs5")` represents the proportion of crimes per 200 residents in each location over time. The goal is to first extract the inequality trend lines such as in Figure `r fig$ref("figs4")` and then cluster them according to the similarity of their slopes. For the `akclustr` function, a user sets the `k` value which may be an integer or a vector of length two specifying the minimum and maximum numbers of clusters to loop through. In the latter case, the `akclustr` function employs either the `Silhouette` (@Rousseeuw1987)) or the `Calinski_Harabasz` score (@Calinski1974;  @Genolini2010) to determine the best cluster solution. In other words, it determines the `k` value that optimizes the specified criterion. The `verbose` argument can be used to control the processing messages. The function is executed as follows:

```{r, echo=TRUE, include=TRUE}

#clustering
akObj <- akclustr(prop_crime_per200_people, id_field = TRUE, 
                                method = "linear", k = c(3,8), crit = "Calinski_Harabasz", verbose=TRUE)

```

In order to preview all the variables of the `quality_plot` object, type: 

```{r, echo=TRUE, message=TRUE, eval=TRUE}

names(akObj)

```

* The description of these variables are as follow:

  + `traj` - returns the input data set used for the clustering.

  + `id_field` - indicates whether the input data set included the id field.

  + `solutions` - the list of cluster solutions by `k` values.

  + `qualitycriterion` - the quality criterion specified.

  + `optimal_k` - the optimal value of `k` as determined by the quality criterion.

  + `qualityCrit.List` - the estimated quality of cluster solutions by `k` values.

  + `qltyplot` - the plot of `qualityCrit.List`, with a red vertical line to indicate the optimal value of `k`.


***Accessing the optimal solution****

The `qualityCrit.List` can be viewed graphically by setting the `quality_plot` argument as `TRUE`. Also, the plot may still be accessed after clustering by printing the variable `akObj$qltyplot`. 


```{r figs6, echo=FALSE, fig.cap=fig$cap("figs6", "Clustering performance at different values of k"), out.width = '80%', fig.align="center"} 
knitr::include_graphics("caliHara.png")
```


From (Figure `r fig$ref("figs6")`), the best value of `k` is the highest at `k=5`, and therefore determined as the best solution. It is recommended that the determination based on either of the quality criteria should be used complementarily with users judgment in relation to the problem at hand. 

Given a value of `k`, the group membership (labels) of its cluster solution can be extracted by entering `k'= k - 2`) into the variable `akObj$solutions[[k']]`. E.g.  

```{r, echo=TRUE, include=TRUE}

#5-group clusters
akObj$solutions[[3]] #for `k=5` solution

```


Also, note that the indexes of the group memberships correspond to that of the trajectory object (`prop_crime_per200_people`) inputted into the function. That is, the membership labels, `"D"`, `"A"`, `"A"`, `....` are the group membership of the trajectories `"E01012628"`,`"E01004768"`,`"E01004803"`,`...` of the object `prop_crime_per200_people`.

## (ii) `print_akstats` function:

The properties (i.e. the descriptive and change statistics) of a cluster solutions (i.e. solution for any value of `k`) such as in `k = 5` above can be generated by  using the special `print` function `print_akstats`. The print function takes as input the `akobject` class, e.g. `akObj`. The descriptive statistics shows the `group memberships` and their `performances` in terms of their shares of the `proportion` measure captured over time. The change statistics shows the information regarding the direction variances of the groups in relation to reference direction. In trajectory clustering analysis, the resulting groups are often re-classified into larger classes based on the slopes, such as `Decreasing`, `Stable`, or `Increasing` classes (@Weisburd2004; Andresen et al. 2017). The slope of a group is the angle made by the medoid of the group relative to a reference line ($R$). The `reference` argument is specified as `1`, `2` or `3`, representing the `mean trajectory`, `medoid trajectory`, or a `horizontal line with slope = 0`, respectively. Let $\vartheta_{1}$ and $\vartheta_{n}$ represent the angular deviations of the group medoids with the lowest slope (negative) and highest (positive) slopes, respectively. 


```{r figs7, echo=FALSE, fig.cap=fig$cap("figs7", "Quantile sub-divisions of most-diverging groups (n_quant=4)"), out.width = '80%', fig.align="center"} 

knitr::include_graphics("Nquant.png")

```

If we sub-divide each of these slopes into a specified number of equal intervals (quantiles), the specific interval within which each group medoid falls can be determined. This specification is made using the `n_quant` argument. Figure `r fig$ref("figs7")` illustrates the quantiles sub-divisions for `n_quant = 4`.

In addition to the slope composition of trajectories found in each group, the quantile location of each group medoid can be used to further categorize the groups into larger classes. We refer users to the package `user manual` for more details about these parameters. Using the current example, the function can be ran as follows: 


```{r, echo=TRUE, include=TRUE}

#Specifying the optimal solution, output$optimal_k (i.e. `k = 5`) and using `stacked` type graph
prpties = print_akstats(akObj, k = 5, show_plots = FALSE)

prpties
```

```{r figs8, echo=FALSE, eval=FALSE, include=FALSE, fig.cap=fig$cap("figs8","group memberships"), out.width = '85%', fig.align="center"} 

knitr::include_graphics("traj_perfm.png")

```

## (iii) `plot_akstats` function:

The above printouts represent the properties (i.e. the descriptive and change properties) of the clusters. Note: the `show_plots` argument of `print_akstats` function, if set as `TRUE`, will produce the plot of group trajectories, representing the group directional change over time. However, the `plot_akstats` has been designed to generate different performance plots of the groups. See below:

(a) ***Group trajectories (directional change over time)***

```{r, echo=TRUE, include=TRUE, fig.width=5,fig.height=5,fig.align="center"}

  #options(rgl.useNULL = TRUE)
  plot_akstats(akObj, k = 5, type="lines", y_scaling="fixed")

```

(b) ***Proportional change of groups change over time***

```{r, echo=TRUE, include=TRUE, fig.width=5,fig.height=5,fig.align="center"}
  #options(rgl.useNULL = TRUE)
  plot_akstats(akObj, k = 5, reference = 1, n_quant = 4, type="stacked")

```

In the context of the long-term inequality study, broad conclusions can be made from both the statistical properties and the plots regarding relative crime exposure in the area represented by each group or class (Adepeju et al. 2021). For example, whilst relative crime exposure have declined in 33.3% (groups `A` and `B`) of the study area, the relative crime exposure have risen in 44.4% (groups `D` and `E`) of the area. The relative crime exposure can be said to be `stable` in 22.2% (group `C`) of the area, based on its close proximity to the reference line. The medoid of the group falls within the $1^{st}(+ve)$ quantile (see Figure `r fig$ref("figs8")`). In essence, we determine that groups `A` and `B` belong to the `Decreasing` class, while groups `D` and `E` belong to the `Increasing` class. 

It is important to state that this proposed classification method is  simply advisory; you may devise a different approach or interpretation depending on your research questions and data.


```{r figs9, echo=FALSE, eval=FALSE, include=FALSE, fig.cap=fig$cap("figs9", "group quality over time"), out.width = '60%', fig.align="center"} 

knitr::include_graphics("traj_perfm2.png")

```

By changing the argument `type="lines"` to `type="stacked"`, a `quality plot` is generated instead (see Figure `r fig$ref("figs9")`). Note that these plots make use of functions within the `ggplot2` library [@Wickham2016]. For a more customized visualization, we recommend that users deploy the `ggplot2` library directly. 


```{r, echo=FALSE, include=FALSE}

col1 <- c("1", "2","3","4","5","6", "7","8","9","10")
col2 <- c("`group`", "`n`", "`n(%)`", "`%Prop.time1`", "`%Prop.timeT`", "`Change`", "`%Change`", "`%+ve Traj.`", "`%-ve Traj.`", "`Qtl:1st-4th`")
col3 <- c("`group membershp`", "`size (no.of.trajectories.)`", "`% size`", "`% proportion of obs. at time 1 (2001)`", "`proportion of obs. at time T (2009)`", "`absolute change in proportion between time1 and timeT`", "`% change in proportion between time 1 and time T`", "`% of trajectories with positive slopes`", "`% of trajectories with negative slopes`", "`Position of a group medoid in the quantile subdivisions`")
tble3 <- data.frame(col1, col2, col3)
tble3 <- tble3

```

```{r table3, results='asis', echo=FALSE, tidy.opts=list(width.cutoff=50)}

knitr::kable(tble3, caption = "Table 3. field description of clustering outputs", col.names = c("SN","field","Description")) %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T, border_right = T) %>%
  column_spec(2, width = "8em", background = "white") %>%
  column_spec(3, width = "12em", background = "white") #%>%
  #row_spec(3:5, bold = T, color = "white", background = "#D7261E")

```


# Conclusion

The `akmedoids` package has been developed in order to aid the replication of a place-based crime inequality investigation conducted in @Adepeju2021. Meanwhile, the utility of the functions in this package are not limited to criminology, but rather can be applicable to longitudinal datasets more generally. This package is being updated on a regular basis to add more functionalities to the existing `functions` and add new functions to carry out other longitudinal data analysis. 

We encourage users to report any bugs encountered while using the package so that they can be fixed immediately. Welcome contributions to this package which will be acknowledged accordingly. 

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_akstats.R
\name{print_akstats}
\alias{print_akstats}
\title{Descriptive (Change) statistics}
\usage{
print_akstats(ak_object, k = 3, reference = 1, n_quant = 4, show_plots = FALSE)
}
\arguments{
\item{ak_object}{An output of \code{\link{akclustr}} function.
The object contains individual trajectories and their cluster
solution(s) at the specified values of \code{k}. Also, includes
the optimal value of \code{k} based on the criterion specified.
at (different) values of \code{k} the \code{traj}.}

\item{k}{[integer] \code{k} cluster to generate its solution.}

\item{reference}{[numeric] Specifying the reference line from
which the direction of each group is measured. Options are:
\code{1}: slope of mean trajectory, \code{2}: slope of medoid
trajectory, \code{3}: slope of a horizontal line
(i.e. slope = 0). Default: \code{1}.}

\item{n_quant}{[numeric] Number of equal intervals (quantiles)
to create between the reference line \code{(R)} and the medoids
\code{(M)} of the most-diverging groups of both sides of
\code{(R)}. Default is \code{4} - meaning quartile subdivisions
on each side of \code{(R)}. In this scenario, the function
returns the quartile in which the medoid of each group falls.
This result can be used to further categorize the groups into
'classes'. For example, groups that fall within the \code{1st}
quartile may be classified as 'Stable' groups (Adepeju et al. 2021).}

\item{show_plots}{[TRUE or FALSE] Provides the trajectory group
plot. Please, see \code{plot_akstats} function for more
plot options.
Defaults \code{FALSE}}
}
\value{
A plot showing group membership or sizes (proportion)
and statistics.
}
\description{
This function perform two tasks:
(i) it generate the descriptive and change statistics
of groups, particularly suited for the outputs form
the \code{\link{akclustr}} function, and
(ii) generates the plots of the groups (performances).
}
\details{
Generates the plot of trajectory groupings.
Given an \code{ak_object} class (from
the \code{akclustr} function), this function show the
plots of cluster groups. The plot component draws from
\code{plot_akstats} function.
}
\examples{

data(traj)

trajectry <- data_imputation(traj, id_field = TRUE, method = 1,
replace_with = 1, fill_zeros = FALSE)

print(trajectry$CompleteData)

trajectry <- props(trajectry$CompleteData, id_field = TRUE)

aksolution <- akclustr(trajectry, id_field = TRUE,
    method = "linear", k = c(3,5), crit='Calinski_Harabasz')

print_akstats(aksolution, k = 4, show_plots=FALSE)

}
\references{
\code{1}. Adepeju, M. et al. (2021). Anchored k-medoids:
A novel adaptation of k-medoids further refined to measure
inequality in the exposure to crime across micro places,
doi: 10.1007/s42001-021-00103-1.

\code{2}. Wickham H. (2016). Elegant graphics for
Data Analysis. Spring-Verlag New York (2016).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_rows_n.R
\name{remove_rows_n}
\alias{remove_rows_n}
\title{Removes rows that contain 'NA' and/or 'Inf' entries}
\usage{
remove_rows_n(traj, id_field=TRUE, remove=1)
}
\arguments{
\item{traj}{[data.frame (numeric)]: longitudinal data.
Each row represents an individual trajectory (of observations).
The columns show the observations at consecutive time points.}

\item{id_field}{[numeric or character] Whether the first column
of the \code{traj} is a unique (\code{id}) field.
Default: \code{FALSE}. If \code{TRUE} the function recognises
the second column as the first time step.}

\item{remove}{[integer] Type of missing entries to remove.
\code{1} for 'NA', \code{2} for 'Inf', and \code{3} for both.
Default:\code{1}.}
}
\value{
A matrix with complete observations
}
\description{
This function removes any rows in which an 'NA'
or an 'Inf' entry is found. The function is also able to
remove records with 'Inf' entries, distinguishing it from
the popular 'na.omit()' function in R.
}
\details{
Given a matrix (or a dataframe) containing an 'NA' or
an 'Inf' entry, the function returns only rows with
complete observations.
}
\examples{

data(traj)

remove_rows_n(traj, id_field=TRUE, remove=3)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_akstats.R
\name{plot_akstats}
\alias{plot_akstats}
\title{Plot of cluster groups.}
\usage{
plot_akstats(
  ak_object,
  k = 3,
  reference = 1,
  n_quant = 4,
  type = "lines",
  y_scaling = "fixed"
)
}
\arguments{
\item{ak_object}{An output of \code{\link{akclustr}} function.
The object contains individual trajectories and their cluster
solution(s) at the specified values of \code{k}. Also, includes
the optimal value of \code{k} based on the criterion specified.
at (different) values of \code{k} the \code{traj}.}

\item{k}{[integer] \code{k} cluster to generate its solution.}

\item{reference}{[numeric] Specifying the reference line from
which the direction of each group is measured. Options are:
\code{1}: slope of mean trajectory, \code{2}: slope of medoid
trajectory, \code{3}: slope of a horizontal line
(i.e. slope = 0). Default: \code{1}.}

\item{n_quant}{[numeric] Number of equal intervals (quantiles)
to create between the reference line \code{(R)} and the medoids
\code{(M)} of the most-diverging groups of both sides of
\code{(R)}. Default is \code{4} - meaning quartile subdivisions
on each side of \code{(R)}. In this scenario, the function
returns the quartile in which the medoid of each group falls.
This result can be used to further categorize the groups into
'classes'. For example, groups that fall within the \code{1st}
quartile may be classified as 'Stable' groups (Adepeju et al. 2021).}

\item{type}{[character] plot type. Available options are:
\code{"lines"} and \code{"stacked"}.}

\item{y_scaling}{[character] works only if \code{type="lines"}.
\code{y_scaling} set the vertical scales of the cluster panels.
Options are: \code{"fixed"}: uses uniform scale for all panels,
\code{"free"}: uses variable scales for panels.}
}
\value{
A plot showing group membership or sizes (proportion)
and statistics.
}
\description{
Takes the 'ak_object' from the
\code{'akclustr'} as input and produce either the 'line' plot
or 'stacked' histogram.
}
\details{
Generates the plots of cluster groups - same plots
generated by the \code{'show_plots'} argument of \code{print_akstats}.
The function draw from the functionalities of the
\code{ggplot2} library.
For a more customized visualisation, we recommend that users
deploy \code{ggplot2} directly (\code{Wickham H. (2016)}).
}
\examples{

data(traj)

trajectry <- data_imputation(traj, id_field = TRUE, method = 1,
replace_with = 1, fill_zeros = FALSE)

print(trajectry$CompleteData)

trajectry <- props(trajectry$CompleteData, id_field = TRUE)

aksolution <- akclustr(trajectry, id_field = TRUE,
method = "linear", k = c(3,5), crit='Calinski_Harabasz')

plot_akstats(aksolution, k = 4, type="lines",
y_scaling="fixed")

plot_akstats(aksolution, k = 4, reference = 1,
n_quant = 4, type="stacked")

}
\references{
\code{1}. Adepeju, M. et al. (2021). Anchored k-medoids:
A novel adaptation of k-medoids further refined to measure
inequality in the exposure to crime across micro places,
doi: 10.1007/s42001-021-00103-1.

\code{2}. Wickham H. (2016). Elegant graphics for
Data Analysis. Spring-Verlag New York (2016).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/props.R
\name{props}
\alias{props}
\title{Conversion of counts (or rates) to 'Proportion'}
\usage{
props(traj, id_field = TRUE, scale = 1, digits = 4)
}
\arguments{
\item{traj}{[matrix (numeric)]: longitudinal data. Each row
represents an individual trajectory (of observations). The
columns show the observations at consecutive time points.}

\item{id_field}{[numeric or character] Whether the first
column of the \code{traj} is a unique (\code{id}) field.
Default: \code{FALSE}. If \code{TRUE} the function recognizes
the second column as the first time step.}

\item{scale}{[numeric] To scale the 'proportion' measures.
Default: \code{1}}

\item{digits}{[numeric] Specifying number of digits to
approximate the output to. Default: \code{4}.}
}
\value{
A dataframe of proportion measures
}
\description{
This function converts counts or rates to proportions.
}
\details{
Given a matrix of observations (counts or rates), this
function converts each observation to a proportion equivalent to
the sum of each column. In other words, each observation is divided
by the sum of the column where it is located, i.e.
\code{prop = [a cell value] / sum[corresponding column]}
}
\examples{

trajectry <- data_imputation(traj, id_field = TRUE, method = 2,
replace_with = 1, fill_zeros = FALSE) #filling the missing values

trajectry <- props(trajectry$CompleteData, id_field = TRUE,
scale=1, digits=4)

print(trajectry)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{clustr}
\alias{clustr}
\title{Sample labels of cluster groups}
\format{
A dataframe containing one variable:
\itemize{
  \item label: alphabetical label by clusters
    }
}
\usage{
clustr
}
\description{
A dataframe of alphabetical labels representing
the optimal solution of `traj` dataset based on `akClust`
function
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TO1Risk}
\alias{TO1Risk}
\title{Time-at-risk for the Adjudicated Toronto
Youth Data (Sample 1)}
\format{
A dataframe with the following variables:
  \itemize{
  \item 8: Time-at-risk per year at age 8
  \item 9: Time-at-risk per year at age 9
  \item 10: Time-at-risk per year at age 10
  \item 11: Time-at-risk per year at age 11
  \item 12: Time-at-risk per year at age 12
  \item 13: Time-at-risk per year at age 13
  \item 14: Time-at-risk per year at age 14
  \item 15: Time-at-risk per year at age 15
  \item 16: Time-at-risk per year at age 16
  \item 17: Time-at-risk per year at age 17
  \item 18: Time-at-risk per year at age 18
  \item 19: Time-at-risk per year at age 19
  \item 20: Time-at-risk per year at age 20
  \item 21: Time-at-risk per year at age 21
  \item 22: Time-at-risk per year at age 22
  \item 23: Time-at-risk per year at age 23
  \item 24: Time-at-risk per year at age 24
  \item 25: Time-at-risk per year at age 25
  \item 26: Time-at-risk per year at age 26
  \item 27: Time-at-risk per year at age 27
  \item 28: Time-at-risk per year at age 28
  \item 29: Time-at-risk per year at age 29
  \item 30: Time-at-risk per year at age 30
  \item 31: Time-at-risk per year at age 31
  \item 32: Time-at-risk per year at age 32
  \item 33: Time-at-risk per year at age 33
  \item 34: Time-at-risk per year at age 34
  \item 35: Time-at-risk per year at age 35
  \item 36: Time-at-risk per year at age 36
  \item 37: Time-at-risk per year at age 37
  \item 38: Time-at-risk per year at age 38
    }
}
\usage{
TO1Risk
}
\description{
Real-life time-at-risk per year for 378 individuals
from the age of 8 to 38 in the Toronto, Ontario, Canada.
The data is obtained through the R package `crimCV`.
For further information, please see: Nielsen, J. (2018)
crimCV: Group-Based Modelling of Longitudinal Data.
R package version 0.9.6.
URL https://CRAN.R-project.org/package=crimCV.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alpha_label.R
\name{alpha_label}
\alias{alpha_label}
\title{Numerics ids to alphabetical ids}
\usage{
alpha_label(x)
}
\arguments{
\item{x}{A vector of numeric ids}
}
\value{
A vector of alphabetical ids.
}
\description{
Function to transform a list of numeric ids
to alphabetic ids
}
\details{
Given a vector of numeric cluster ids,
`alpha_label` converts each id to its corresponding alphabets.
It combines alphabets for ids greater than 26.
}
\examples{

data(TO1Risk)

set.seed(1000)
#pick 4 random clusters
center <- TO1Risk[runif(4,1,nrow(TO1Risk)), ]

#Assigning each individual to nearest centre
numeric_Labels <- kml::affectIndivC(TO1Risk, center)

mode(numeric_Labels)

#transform numeric cluster labels to alphabets
alphab_Labels <- alpha_label(numeric_Labels)

mode(alphab_Labels)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elbow_point.R
\name{elbow_point}
\alias{elbow_point}
\title{Determine the elbow point on a curve}
\usage{
elbow_point(x, y)
}
\arguments{
\item{x}{vector of x coordinates of points on the curve}

\item{y}{vector of y coordinates of points on the curve}
}
\value{
indicate the optimal k value determined by
the elbow point point.
}
\description{
Given a list of x, y coordinates on a curve,
function determines the elbow point of the curve.
}
\details{
highlight the maximum curvature to identify the
elbow point (credit: 'github.com/agentlans')
}
\examples{

# Generate some curve
x <- runif(100, min=-2, max=3)
y <- -exp(-x) * (1+rnorm(100)/3)
plot(x, y)
#Plot elbow points
abline(v=elbow_point(x,y)$y, col="blue", pch=20, cex=3)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rates.R
\name{rates}
\alias{rates}
\title{Conversion of counts to rates}
\usage{
rates(traj, denomin, id_field, multiplier)
}
\arguments{
\item{traj}{[matrix (numeric)] longitudinal (e.g.
observed count) data (\code{m x n}). Each row represents an
individual trajectory (of observations). The columns show the
observations at consecutive time steps.}

\item{denomin}{[matrix (numeric)] longitudinal (denominator)
data of the same column as `traj` (\code{n}).}

\item{id_field}{[numeric or character] Default is \code{TRUE}.
The first column of both the `traj` and the `denomin` object
must be the unique (\code{id}) field. If \code{FALSE}, the function
will terminate. The assumption is that columns of both the
\code{traj} and \code{denominat} corresponds. That is, column2,
column3, ... represent time points 2, 3, ..., respectively, in
each object.}

\item{multiplier}{[numeric] A quantify by which to the ratio
\code{traj/denomin} is expressed. Default is \code{100}.}
}
\value{
An object which comprised of four output variables, namely:
(i) `$common_ids` - individual ids present in both
`traj` (trajectory data) and `denomin` (denominator data);
(ii) `$ids_unique_to_traj_data` - individual ids unique to
trajectory data (i.e. not present in the denominator data);
(iii) `$ids_unique_to_denom_data` - individual ids unique
to denominator data (i.e. not present in the trajectory data);
(iv) `` - a dataframe of rates estimates. Note: only the individual
ids in `$rates_estimates` are used in the `rates` estimation.
}
\description{
Calculates rates from 'observed' count and a
denominator data
}
\examples{

traj2 <- data_imputation(traj, id_field = TRUE, method = 2,
replace_with = 1, fill_zeros = FALSE)

pop <- popl #read denominator data

pop2 <- as.data.frame(matrix(0, nrow(popl), ncol(traj)))

colnames(pop2) <- names(traj2$CompleteData)

pop2[,1] <- as.vector(as.character(pop[,1]))

pop2[,4] <- as.vector(as.character(pop[,2]))

pop2[,8] <- as.vector(as.character(pop[,3]))

list_ <- c(2, 3, 5, 6, 7, 9, 10) #vector of missing years

#fill the missing fields with 'NA'
for(u_ in seq_len(length(list_))){
    pop2[,list_[u_]] <- "NA"
}

#estimate missing fields
pop_imp_result <- data_imputation(pop2, id_field = TRUE, method = 2,
replace_with = 1, fill_zeros = FALSE)

#calculate rates i.e. crimes per 200 population
crime_rates <- rates(traj2$CompleteData, denomin=pop_imp_result$CompleteData,
id_field=TRUE, multiplier = 200)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_imputation.R
\name{data_imputation}
\alias{data_imputation}
\title{Data imputation for longitudinal data}
\usage{
data_imputation(traj, id_field = FALSE, method = 2,
replace_with = 1, fill_zeros = FALSE, verbose=TRUE)
}
\arguments{
\item{traj}{[\code{matrix (numeric)}]: longitudinal data. Each
row represents an individual trajectory (of observations). The
columns show the observations at consecutive time points.}

\item{id_field}{[numeric or character] Whether the first column
of the \code{traj} is a unique (\code{id}) field. Default:
\code{FALSE}. If \code{TRUE} the function recognises the second
column as the first time step.}

\item{method}{[an integer] indicating a method for calculating
the missing values. Options are: \code{'1'}: \code{arithmetic}
 method, and \code{'2'}: \code{regression} method. The default
 is \code{'1'}: \code{arithmetic} method}

\item{replace_with}{[an integer from 1 to 6] indicating the technique,
based on a specified \code{method}, for calculating the missing entries.
\code{'1'}: \code{arithmetic} method, \code{replace_with} options are:
\code{'1'}: Mean value of the corresponding column;
\code{'2'}: Minimum value of corresponding column; \code{'3'}:
Maximum value of corresponding column;
\code{'4'}: Mean value of corresponding row; \code{'5'}:
Minimum value of corresponding row,
or \code{'6'}: Maximum value of corresponding row. For \code{'2'}:
regression method:
the available option for the \code{replace_with} is: \code{'1'}:
\code{linear}.
The regression method fits a linear regression line to a trajectory
with missing entry(s)
and estimates the missing data values from the regression line.
Note: only the missing data points derive their new values from the
regression line
while the rest of the data points retain their original values. The
function terminates if there are
trajectories with only one observation. The default is \code{'1'}: Mean
value of the corresponding column}

\item{fill_zeros}{[TRUE or FALSE] whether to consider zeros \code{0}
as missing values when \code{2: regression} method is used. The default
is \code{FALSE}.}

\item{verbose}{to suppress printing output messages (to the console).
Default: \code{TRUE}.}
}
\value{
A data.frame with missing values (\code{NA}, \code{Inf},
\code{null}) imputed according to the a specified technique.
}
\description{
This function fills any missing entries (\code{NA},
\code{Inf}, \code{null}) in a matrix or dataframe, according to
a specified method. By default, \code{'0'} is considered a value.
}
\details{
Given a matrix or data.frame with some missing values
indicated by (\code{NA}, \code{Inf}, \code{null}), this function
impute the missing value by using either an estimation from the
corresponding rows or columns, or to use a regression method to
estimate the missing values.
}
\examples{

# Using the example 'traj' datasets

imp_data <- data_imputation(traj, id_field = TRUE, method = 2,
replace_with = 1,
fill_zeros = FALSE, verbose=FALSE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outlier_detect.R
\name{outlier_detect}
\alias{outlier_detect}
\title{Outlier detection and replacement}
\usage{
outlier_detect(traj, id_field = FALSE, method = 1, threshold = 0.95,
count = 1, replace_with = 1, verbose=TRUE)
}
\arguments{
\item{traj}{[matrix (numeric)]: longitudinal data. Each row
represents an individual trajectory (of observations). The columns
show the observations at consecutive time points.}

\item{id_field}{[numeric or character] Whether the first column
of the \code{traj} is a unique (\code{id}) field.
Default: \code{FALSE}. If \code{TRUE} the function recognizes
the second column as the first time step.}

\item{method}{[integer (numeric)] indicating the method for
identifying the outlier. Options are: \code{'1'}: quantile method
(\code{default}), and \code{'2'}: manual method. The \code{manual}
method requires a user-defined value.}

\item{threshold}{[numeric] A cut-off value for outliers. If the
\code{method} parameter is set as \code{'1'}:quantile, the \code{threshold}
should be a numeric vector of probability between \code{[0,1]}, whilst if
the \code{method} is set as \code{'2'}: \code{manual}, the
\code{threshold} could be any numeric vector.}

\item{count}{[integer (numeric)] indicating the number of observations
(in a trajectory) that must exceed the \code{threshold} in order for the
trajectory to be considered an \code{outlier}. Default is \code{1}.}

\item{replace_with}{[integer (numeric)] indicating the technique to
use for calculating a replacement for an outlier observation. The remaining
observations on the row or the column in which the outlier observation is
located are used to calculate the replacement.
The replacement options are: \code{'1'}: Mean value of the column,
\code{'2'}: Mean value of the row and \code{'3'}: remove the row
(trajectory) completely from the data. Default value is the
\code{'1'} option.}

\item{verbose}{to suppress output messages (to the console).
Default: \code{TRUE}.}
}
\value{
A dataframe with outlier observations replaced or removed.
}
\description{
This function identifies outlier observations
in the trajectories, and allows users to replace the observations
or remove trajectories entirely.
}
\details{
Given a matrix, this function identifies outliers that
exceed the threshold and replaces the outliers with an estimate
calculated using the other observations either the rows or the columns
in which the outlier observation is located. Option is also provided to
remove the trajectories (containing the outlier) from the data.
}
\examples{

data(traj)

trajectry <- data_imputation(traj, id_field=TRUE, method = 1,
   replace_with = 1, verbose=FALSE)

trajectry <- props(trajectry$CompleteData, id_field=TRUE)

outp <- outlier_detect(trajectry, id_field = TRUE, method = 1,
threshold = 0.95, count = 1, replace_with = 1, verbose=TRUE)

outp <- outlier_detect(trajectry, id_field = TRUE, method = 2, threshold = 15,
  count = 4, replace_with = 3, verbose=TRUE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{popl}
\alias{popl}
\title{Simulated population data.}
\format{
A dataframe with the following variables:
\itemize{
  \item location_id: Character id of sample census unit
  at which the population is obtained.
  \item census_2003: Population estimates at the sample
  locations for the census year 2003.
  \item census_2007: Population estimates at the sample
  locations for the census year 2007.
    }
}
\usage{
popl
}
\description{
Sample simulated population data to be used as
the denominator variable against `traj` dataset. Contains
data for two consecutive census years
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{simulated}
\alias{simulated}
\title{Simulated longitudinal dataset}
\format{
A dataframe with the following variables:
\itemize{
  \item X1: Values across locations at time step 1
  \item X2: Values across locations at time step 2
  \item X3: Values across locations at time step 3
  \item X4: Values across locations at time step 4
  \item X5: Values across locations at time step 5
  \item X6: Values across locations at time step 6
  \item X7: Values across locations at time step 7
  \item X8: Values across locations at time step 8
  \item X9: Values across locations at time step 9
  \item X10: Values across locations at time step 10
  \item X11: Values across locations at time step 11
  \item X12: Values across locations at time step 12
  \item X13: Values across locations at time step 13
  \item X14: Values across locations at time step 14
  \item X15: Values across locations at time step 15
  \item X16: Values across locations at time step 16
  \item X17: Values across locations at time step 17
  \item X18: Values across locations at time step 18
  \item X19: Values across locations at time step 19
  \item X20: Values across locations at time step 20
  \item X21: Values across locations at time step 21
    }
}
\usage{
simulated
}
\description{
Contains simulated trajectories belonging to one of the
three pre-defined groups, namely (a) decreasing, (b) stable
and (c) increasing groups.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{traj}
\alias{traj}
\title{Sample longitudinal dataset}
\format{
A dataframe with the following variables:
  \itemize{
  \item location_ids: Character id of sample locations
  at which values are obtained.
  \item X2001: Values at time step 1 (i.e. year 2001)
  \item X2002: Values at time step 2 (i.e. year 2002)
  \item X2003: Values at time step 3 (i.e. year 2003)
  \item X2004: Values at time step 4 (i.e. year 2004)
  \item X2005: Values at time step 5 (i.e. year 2005)
  \item X2006: Values at time step 6 (i.e. year 2006)
  \item X2007: Values at time step 7 (i.e. year 2007)
  \item X2008: Values at time step 8 (i.e. year 2008)
  \item X2009: Values at time step 9 (i.e. year 2009)
    }
}
\usage{
traj
}
\description{
Simulated longitudinal datasets containing
trajectories with missing values
(\code{NA}, \code{Inf}, \code{null})
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/w_spaces.R
\name{w_spaces}
\alias{w_spaces}
\title{Whitespaces removal}
\usage{
w_spaces(traj, remove="Both", verbose=TRUE)
}
\arguments{
\item{traj}{[matrix (numeric)]: longitudinal data. Each row
represents an individual trajectory (of observations).
The columns show the observations at consecutive time points.}

\item{remove}{[character]: Type of whitespace to remove.
That is, "Left" (leading), (2) "Right" (trailing), or "Both"
(both leading and trailing whitespaces). Default: "Both".}

\item{verbose}{to suppress output messages (to the console).
Default: \code{TRUE}.}
}
\value{
A matrix with all whitespaces (if any) removed.
}
\description{
This function removes all the leading and the
trailing whitespaces in data
}
\details{
Given a matrix suspected to contain whitespaces,
this function removes the type of the whitespaces specified and returns a
cleaned data. ’Whitespaces’ are white characters often
introduced into data during data entry, for instance by
wrongly pressing the spacebar. For example, neither " A"
nor "A " is the same as "A" because of the whitespaces that
exist in them. They can also result from systematic
errors in data recording devices.
}
\examples{

data(traj_w_spaces)

w_spaces(traj_w_spaces, remove="Both", verbose=TRUE)

}
\references{
\url{https://en.wikipedia.org/wiki/Whitespace_character}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/akclustr.R
\name{akclustr}
\alias{akclustr}
\title{Anchored k-medoids clustering}
\usage{
akclustr(traj, id_field = FALSE, method = "linear",
k = c(3,6), crit="Silhouette", verbose = TRUE, quality_plot=FALSE)
}
\arguments{
\item{traj}{[matrix (numeric)]: longitudinal data. Each row represents
an individual trajectory (of observations). The columns show the
observations at consecutive time steps.}

\item{id_field}{[numeric or character] Whether the first column of the
\code{traj} is a unique (\code{id}) field. Default: \code{FALSE}.
If \code{TRUE} the function recognizes the second column as the first
time points.}

\item{method}{[character] The parametric initialization strategy.
Currently, the only available method is a \code{linear} method, set as
\code{"linear"}. This uses the time-dependent linear regression lines
and the resulting groups are order in the order on increasing slopes.}

\item{k}{[integer or vector (numeric)] either an exact integer number
of clusters, or a vector of length two specifying the minimum and
maximum numbers of clusters to be examined from which the best
solution will be determined. In either case, the minimum number
of clusters is \code{3}. The default is \code{c(3,6)}.}

\item{crit}{[character] a string specifying the type of the criterion
to use for assessing the quality of the cluster solutions, when
\code{k} is a vector of two values (as above). Default:
\code{crit="Silhouette"}, use the average Silhouette width
(\code{Rousseeuw P. J. 1987}). Using the \code{"Silhouette"} criterion,
the optimal value of \code{k} can be determined as the elbow point of
the curve. Other valid criterion is the "Calinski_Harabasz"
(\code{Caliński T. & Harabasz J. 1974}) in which the maximum score
represents the point of optimality. Having determined the optimal
\code{k}, the function can then be re-run, using the exact (optimal)
value of \code{k}.}

\item{verbose}{to suppress output messages (to the console)
during clustering. Default: \code{TRUE}.}

\item{quality_plot}{Whether to show plot of quality criteria across
different values of \code{k}. Default: \code{FALSE}.}
}
\value{
generates an \code{akobject} consisting of the
cluster solutions at the specified values of \code{k}. Also,
the graphical plot of the quality scores of the cluster
solutions.
}
\description{
Given a list of trajectories and a functional method,
this function clusters the trajectories into a \code{k} number of
groups. If a vector of two numbers is given, the function determines
the best solution from those options based on the Caliński-Harabasz
criterion.
}
\details{
This function works by first approximating the trajectories
based on the chosen parametric forms (e.g. linear), and then partitions
the original trajectories based on the form groupings, in similar
fashion to k-means clustering \code{(Genolini et al. 2015)}. The key
distinction of \code{akmedoids} compared with existing longitudinal
approaches is that both the initial starting points as well as the
subsequent cluster centers (as the iteration progresses) are based
the selection of observations (medoids) as oppose to centroids.
}
\examples{

data(traj)

trajectry <- data_imputation(traj, id_field = TRUE, method = 2,
replace_with = 1, fill_zeros = FALSE)

trajectry <- props(trajectry$CompleteData, id_field = TRUE)

print(trajectry)

output <- akclustr(trajectry, id_field = TRUE,
method = "linear", k = c(3,7), crit='Calinski_Harabasz',
verbose = FALSE, quality_plot=FALSE)

print(output)

}
\references{
\code{1}. Genolini, C. et al. (2015) kml and kml3d:
R Packages to Cluster Longitudinal Data. Journal of Statistical
Software, 65(4), 1-34. URL http://www.jstatsoft.org/v65/i04/.

\code{2}. Rousseeuw P. J. (1987) Silhouettes: A graphical aid
to the interpretation and validation of cluster analysis.
J. Comput. Appl. Math 20:53–65.

\code{3}. Caliński T, Harabasz J (1974) A dendrite method for
cluster analysis. Commun. Stat. 3:1-27.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{traj_w_spaces}
\alias{traj_w_spaces}
\title{Sample longitudinal dataset containing whitespaces}
\format{
A dataframe with the following variables:
  \itemize{
  \item location_ids: Character id of sample locations
  at which values are obtained.
  \item X2001: Values at time step 1 (i.e. year 2001)
  \item X2002: Values at time step 2 (i.e. year 2002)
  \item X2003: Values at time step 3 (i.e. year 2003)
  \item X2004: Values at time step 4 (i.e. year 2004)
  \item X2005: Values at time step 5 (i.e. year 2005)
  \item X2006: Values at time step 6 (i.e. year 2006)
  \item X2007: Values at time step 7 (i.e. year 2007)
  \item X2008: Values at time step 8 (i.e. year 2008)
  \item X2009: Values at time step 9 (i.e. year 2009)
    }
}
\usage{
traj_w_spaces
}
\description{
Longitudinal dataset with both trailing and leading
whitespaces. For example, there is a trailing whitespace
at cell [3, 6], while there is a leading whitespace
at cell [9, 4].
}
\keyword{datasets}
