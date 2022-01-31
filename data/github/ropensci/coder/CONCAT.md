
# coder <img src="man/figures/logo.png" align="right"/>

[![R build
status](https://github.com/ropensci/coder/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/coder/actions)
[![codecov](https://codecov.io/gh/ropensci/coder/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/coder)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/65975808.svg)](https://zenodo.org/badge/latestdoi/65975808)
[![status](https://joss.theoj.org/papers/10.21105/joss.02916/status.svg)](https://joss.theoj.org/papers/10.21105/joss.02916)
[![CRAN
status](https://www.r-pkg.org/badges/version/coder)](https://CRAN.R-project.org/package=coder)
![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/coder)

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Aim of the package

The goal of `{coder}` is to classify items from one dataset, using codes
from a secondary source with classification schemes based on regular
expressions and weighted indices.

## Installation

You can install the released version of coder from
[CRAN](https://CRAN.R-project.org) with:

``` r
# install.packages("coder")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("eribul/coder")
```

## Typical use case

-   Determining comorbidities before clinical trials
-   Discovering adverse events after surgery

**Patient data:** The initial rationale for the package was to classify
patient data based on medical coding. A typical use case would consider
patients from a medical/administrative data base, as identified by some
patient id and possibly with some associated date of interest (date of
diagnoses/treatment/intervention/hospitalization/rehabilitation). This
data source could be for example an administrative hospital register or
a national quality register.

**Codify:** The primary source could then be linked to a secondary
(possibly larger) data base including the same patients with
corresponding id:s and some coded patient data. This could be a national
patient register with medical codes from the International
Classification of Diseases *(ICD)* with corresponding dates of hospital
visits/admission/discharge, or a medical prescription register with
codes from the Anatomic Therapeutic Chemical *(ATC)* classification
system with dates of medical prescription/dispatch/usage. A time window
could be specified relating the date of the primary source (i. e. the
date of a primary total hip arthroplasty; THA), to dates from the
secondary source (i.e. the date of a medical prescription). ATC codes
associated with medical prescriptions during one year prior to THA,
could thus be identified and used as a measure of comorbidity. Another
time window of 90 days after THA, might instead be used to identify
adverse events after surgery.

**Classify:** To work with medical/chemical codes directly might be
cumbersome, since those classifications tend to be massive and therefore
hard to interpret. It is thus common to use data aggregation as proposed
by some classification or combined index from the literature. This could
be the *Charlson* or *Elixhauser* comorbidity indices based on
ICD-codes, or the *RxRisk V* classification based on ATC-codes. Each of
those tools appear with different code versions (ICD-8, ICD-9, ICD-9-CM,
ICD-10, ICD-10-CA, ICD-10-SE, ICD-10-CM et cetera) and with different
codes recognized as relevant comorbidities (the Charlson index proposed
by Charlson et al, Deyo et al, Romano et al. Quan et al. et cetera).
Using a third object (in addition to the primary and secondary patient
data sets) helps to formalize and structure the use of such
classifications. This is implemented in the `coder` package by
`classcodes` objects based on regular expressions (often with several
alternative versions). Those `classcodes` objects could be prepared by
the user, although a number of default `classcodes` are also included in
the package (table below).

**Index:** Now, instead of working with tens of thousands of individual
ICD-codes, each patient might be recognized to have none or some
familiar comorbidity such as hypertension, cancer or dementia. This
granularity might be too fine-grained still, wherefore an even simpler
index score might be searched for. Such scores/indices/weighted sums
have been proposed as well and exist in many versions for each of the
standard classifications. Some are simple counts, some are weighted
sums, and some accounts for some inherited hierarchy (such that
ICD-codes for diabetes with and without complications might be
recognized in the same patient, although the un-complicated version
might be masked by the complicated version in the index).

**Conditions:** Some further complexity might appear if some codes are
only supposed to be recognized based on certain conditions. Patients
with THA for example might have an adverse event after surgery if a
certain ICD-code is recorded as the main diagnose at a later hospital
visit, although the same code could be ignored if recorded only as a
secondary diagnosis.

**To summarize:** The coder package takes three objects: (1) a data
frame/table/tibble with id and possible dates from a primary source; (2)
coded data from a secondary source with the same id and possibly
different dates and; (3) a `classcodes` object, either a default one
from the package, or as specified by the user. The outcome is then: (i)
codes associated with each element from (1) identified from (2),
possibly limited to a relevant time window; (ii) a broader
categorization of the relevant codes as prescribed by (3), and; (iii) a
summarized index score based on the relevant categories from (3).

(i-iii) corresponds to the output from functions `codify()`,
`classify()` and `index()`, which could be chained explicitly as
`codify() %>% classify() %>% index()`, or implicitly by the
`categorize()` function.

## Usage

Assume we have some patients with surgery at specified dates:

``` r
library(coder)
ex_people
#> # A tibble: 100 x 2
#>    name              surgery   
#>    <chr>             <date>    
#>  1 Chen, Trevor      2021-01-04
#>  2 Graves, Acineth   2020-09-26
#>  3 Trujillo, Yanelly 2020-09-13
#>  4 Simpson, Kenneth  2020-12-16
#>  5 Chin, Nelson      2020-11-29
#>  6 Le, Christina     2020-07-03
#>  7 Kang, Xuan        2020-10-05
#>  8 Shuemaker, Lauren 2020-07-04
#>  9 Boucher, Teresa   2020-12-10
#> 10 Le, Soraiya       2020-11-14
#> # ... with 90 more rows
```

Those patients (among others) were also recorded in a national patient
register with date of hospital admissions and diagnoses codes coded by
the International Classification of Diseases (ICD) version 10:

``` r
ex_icd10
#> # A tibble: 2,376 x 4
#>    name                 admission  icd10 hdia 
#>    <chr>                <date>     <chr> <lgl>
#>  1 Tran, Kenneth        2020-07-18 S134A FALSE
#>  2 Tran, Kenneth        2021-01-01 W3319 FALSE
#>  3 Tran, Kenneth        2020-12-11 Y0262 TRUE 
#>  4 Tran, Kenneth        2020-11-03 X0488 FALSE
#>  5 Sommerville, Dominic 2020-12-23 V8104 FALSE
#>  6 Sommerville, Dominic 2020-08-03 B853  FALSE
#>  7 Sommerville, Dominic 2020-12-18 Q174  FALSE
#>  8 Sommerville, Dominic 2020-08-08 A227  FALSE
#>  9 Sommerville, Dominic 2020-12-13 H702  FALSE
#> 10 Sommerville, Dominic 2020-04-06 X6051 TRUE 
#> # ... with 2,366 more rows
```

Using those two data sets, as well as a classification scheme
(`classcodes` object; see below), we can easily identify all Charlson
comorbidities for each patient:

``` r
ch <- 
  categorize(
    ex_people,                  # patients of interest 
    codedata = ex_icd10,        # Medical codes from national patient register
    cc = "charlson",            # Calculate Charlson comorbidity
    id = "name", code = "icd10" # Specify column names
  )
#> Classification based on: icd10

ch
#> # A tibble: 100 x 25
#>    name  surgery    myocardial.infa~ congestive.hear~ peripheral.vasc~
#>    <chr> <date>     <lgl>            <lgl>            <lgl>           
#>  1 Chen~ 2021-01-04 FALSE            FALSE            FALSE           
#>  2 Grav~ 2020-09-26 FALSE            FALSE            FALSE           
#>  3 Truj~ 2020-09-13 FALSE            FALSE            FALSE           
#>  4 Simp~ 2020-12-16 FALSE            FALSE            FALSE           
#>  5 Chin~ 2020-11-29 FALSE            FALSE            FALSE           
#>  6 Le, ~ 2020-07-03 FALSE            FALSE            FALSE           
#>  7 Kang~ 2020-10-05 FALSE            FALSE            FALSE           
#>  8 Shue~ 2020-07-04 FALSE            FALSE            FALSE           
#>  9 Bouc~ 2020-12-10 FALSE            FALSE            TRUE            
#> 10 Le, ~ 2020-11-14 FALSE            FALSE            FALSE           
#> # ... with 90 more rows, and 20 more variables: cerebrovascular.disease <lgl>,
#> #   dementia <lgl>, chronic.pulmonary.disease <lgl>, rheumatic.disease <lgl>,
#> #   peptic.ulcer.disease <lgl>, mild.liver.disease <lgl>,
#> #   diabetes.without.complication <lgl>, hemiplegia.or.paraplegia <lgl>,
#> #   renal.disease <lgl>, diabetes.complication <lgl>, malignancy <lgl>,
#> #   moderate.or.severe.liver.disease <lgl>, metastatic.solid.tumor <lgl>,
#> #   AIDS.HIV <lgl>, charlson <dbl>, deyo_ramano <dbl>, dhoore <dbl>,
#> #   ghali <dbl>, quan_original <dbl>, quan_updated <dbl>
```

How many patients were diagnosed with malignancy?

``` r
sum(ch$malignancy)
#> [1] 5
```

What is the distribution of the combined comorbidity index for each
patient?

``` r
barplot(table(ch$charlson))
```

<img src="man/figures/READMEunnamed-chunk-5-1.png" width="100%" />

There are many versions of the Charlson comorbidity index, which might
be controlled by the `index` argument. We might also be interested only
in diagnoses from 90 days before surgery as specified with an argument
list `codify_args`as passed to `codify()`:

``` r
ch <- 
  categorize(
    ex_people, codedata = ex_icd10, cc = "charlson", id = "name", code = "icd10",
    
    # Additional arguments
    index       = c("quan_original", "quan_updated"), # Indices
    codify_args = list(
      date      = "surgery",   # Name of column with index dates
      code_date = "admission", # Name of column with code dates
      days      = c(-90, -1)   # Time window
    )
  )
#> Classification based on: icd10
```

Number of malignancies during this period?

``` r
sum(ch$malignancy, na.rm = TRUE)
#> [1] 3
```

Distribution of the index as proposed by Quan et al 2011 during the 90
day period:

``` r
barplot(table(ch$quan_updated))
```

<img src="man/figures/READMEunnamed-chunk-8-1.png" width="100%" />

## Classification schemes

Classification schemes (`classcodes` objects, see
`vignette("classcodes")`) are based on regular expressions for
computational speed (see `vignette("Interpret_regular_expressions")`),
but their content can be summarized and visualized for clarity.
Arbitrary `classcodes` objects can also be specified by the user.

The package includes default `classcodes` for medical patient data based
on the international classification of diseases version 8, 9 and 10
(ICD-8/9/10), as well as the Anatomical Therapeutic Chemical
Classification System (ATC) for medical prescription data.

Default `classcades` are listed in the table. Each classification
(classcodes column) can be based on several code systems (regex column)
and have several alternative weighted indices (indices column). Those
might be combined freely.

``` r
coder::all_classcodes()
#> # A tibble: 7 x 3
#>   classcodes   regex                             indices                        
#>   <chr>        <chr>                             <chr>                          
#> 1 charlson     icd10, icd9cm_deyo, icd9cm_enhan~ "charlson, deyo_ramano, dhoore~
#> 2 cps          icd10                             "only_ordinary"                
#> 3 elixhauser   icd10, icd10_short, icd9cm, icd9~ "sum_all, sum_all_ahrq, walrav~
#> 4 hip_ae       icd10, kva, icd10_fracture        ""                             
#> 5 hip_ae_hail~ icd10, kva                        ""                             
#> 6 knee_ae      icd10, kva                        ""                             
#> 7 rxriskv      atc_pratt, atc_caughey, atc_garl~ "pratt, sum_all"
```

# Relation to other packages

`coder` uses `data.table` as a backend to increase computational speed
for large datasets. There are some R packages with a narrow focus on
Charlson and Elixhauser co-morbidity based on ICD-codes
([icd](https://CRAN.R-project.org/package=icd),
[comorbidity](https://CRAN.R-project.org/package=comorbidity),
[medicalrisk](https://CRAN.R-project.org/package=medicalrisk),
[comorbidities.icd10](https://github.com/gforge/comorbidities.icd10),
[icdcoder](https://github.com/wtcooper/icdcoder)). The `coder` package
includes similar functionalities but has a wider scope.

# Code of conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
# coder 0.13.1

Release candidate submitted to JOSS.

# coder 0.5.5

* Include Walraven score for Elixhauser comorbidity index

# coder 0.5.3.9000

* Include method summary.classcodes (#39)


# coder 0.5.2.9000

- Internal changes only.


---
title: 'coder: An R package for code-based item classification and categorization'
tags:
  - R
  - epidemiology
  - statistics
  - Administrative data
authors:
  - name: Erik Bülow
    orcid: 0000-0002-9973-456X
    affiliation: "1, 2"
affiliations:
 - name: The Swedish Arthroplasty Register, Registercentrum Västra Götaland, Gothenburg, Sweden
   index: 1
 - name: Department of Orthopaedics, Institute of Clinical Sciences, Sahlgrenska Academy, University of Gothenburg, Gothenburg, Sweden
   index: 2
date: 16 December 2020
bibliography: paper.bib
---

# Summary

The `coder` package lets researchers classify and categorize coded item data using pre-specified classification schemes based on regular expressions. Default classification schemes are included for commonly used medical and clinical classifications. The package implementation aim for high performance and the package can be used with large data sets.


# Medical coding and classifications

Registry based research and the use of real world evidence (RWE) and data (RWD) have gained popularity over the last years [@Sherman2016], both as an epidemiological research tool, and for monitoring post market safety and adverse events due to regulatory decisions. Data from administrative, clinical and medical registries are often coded based on standardized classifications for diagnostics, procedures/interventions, medications/medical devices and health status/functioning.

Codes and classifications are maintained and developed by several international bodies, such as The World Health Organization [(WHO)](https://www.who.int/classifications/), [SNOMED International](snomed.org), and the Nordic Medico-Statistical Committee [(NOMESCO)](http://nowbase.org/). 


# Challanges 

Common classifications such as the International Classification of Diseases (ICD) or the Anatomical Therapeutic Chemical Classification System (ATC) entails thousands of codes which are hard to use and interpret in applied research. This is often solved by an abstraction layer combining individual codes into broader categories, sometimes further simplified by a single index value based on a weighted sum of individual categories [@Charlson1987; @Elixhauser1998; @Quan2005; @Sloan2003; @Pratt2018]. 
 
 
# Statement of Need

Large and long-standing national databases often contain millions of entries and span several Gigabytes (GB) in size. This leads to high computational burden and a time-consuming data managing process, a cumbersome but necessary prerequisite before any relevant analysis can be performed. There are several R-packages with a deliberate focus on comorbidity data coded by ICD and summarized by the Charlson or Elixhauser comorbidity indices ([icd](https://jackwasey.github.io/icd), [comorbidity](https://ellessenne.github.io/comorbidity/) [@Gasparini2018] and [medicalrisk](https://github.com/patrickmdnet/medicalrisk)). The `coder` package includes such capabilities as well, but takes a more general approach to deterministic item classification and categorization.


# The coder package

[coder](https://docs.ropensci.org/coder/) is an R package with a scope to combine items (i.e. patients) with generic code sets, and to classify and categorize such data based on generic classification schemes defined by regular expressions. It is easy to combine different classifications (such as multiple versions of ICD, ATC or NOMESCO codes), with different classification schemes (such as Charlson, Elixhauser, RxRisk V or for example local definitions of adverse events after total hip arthroplasty) and different weighted indices based on those classifications. The package includes default classification schemes for all those settings, as well as an infrastructure to implement and visualize custom classification schemes. Additional functions simplify identification of codes and events within limited time frames, such as comorbidity during one year before surgery or adverse events within 30 days after. `coder` can also be used in tandem with [decoder](https://cancercentrum.bitbucket.io/decoder/), a package facilitating interpretation of individual codes.

`coder` has been optimized for speed and large data sets using reference semantics from [data.table](https://rdatatable.gitlab.io/data.table/), matrix-based computations and code profiling. The prevalence of large datasets makes it difficult to use parallel computing however, since the limit of available random-access memory (RAM) often implies a more serious bottleneck, which limits the possibility to manifold data sets for multiple cores.

`coder` has been used in ongoing, as well as in previously published research [@Bulow2017; @Bulow2019; @Bulow2020; @Cnudde2017; @Cnudde2018a; @Berg2018; @Jawad2019; @Wojtowicz2019; @Hansson2020; @Nemes2018a].


# References
---
output: github_document
---

# coder <img src="man/figures/logo.png" align="right"/>

[![R build status](https://github.com/ropensci/coder/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/coder/actions) [![codecov](https://codecov.io/gh/ropensci/coder/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/coder) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![DOI](https://zenodo.org/badge/65975808.svg)](https://zenodo.org/badge/latestdoi/65975808) [![status](https://joss.theoj.org/papers/10.21105/joss.02916/status.svg)](https://joss.theoj.org/papers/10.21105/joss.02916)
[![CRAN status](https://www.r-pkg.org/badges/version/coder)](https://CRAN.R-project.org/package=coder)
![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/coder)

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README",
  out.width = "100%"
) 
```

## Aim of the package

The goal of `{coder}` is to classify items from one dataset, using codes from a secondary source with classification schemes based on regular expressions and weighted indices.

## Installation

You can install the released version of coder from [CRAN](https://CRAN.R-project.org) with:

``` r
# install.packages("coder")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("eribul/coder")
```

## Typical use case

-   Determining comorbidities before clinical trials
-   Discovering adverse events after surgery

**Patient data:** The initial rationale for the package was to classify patient data based on medical coding. A typical use case would consider patients from a medical/administrative data base, as identified by some patient id and possibly with some associated date of interest (date of diagnoses/treatment/intervention/hospitalization/rehabilitation). This data source could be for example an administrative hospital register or a national quality register.

**Codify:** The primary source could then be linked to a secondary (possibly larger) data base including the same patients with corresponding id:s and some coded patient data. This could be a national patient register with medical codes from the International Classification of Diseases *(ICD)* with corresponding dates of hospital visits/admission/discharge, or a medical prescription register with codes from the Anatomic Therapeutic Chemical *(ATC)* classification system with dates of medical prescription/dispatch/usage. A time window could be specified relating the date of the primary source (i. e. the date of a primary total hip arthroplasty; THA), to dates from the secondary source (i.e. the date of a medical prescription). ATC codes associated with medical prescriptions during one year prior to THA, could thus be identified and used as a measure of comorbidity. Another time window of 90 days after THA, might instead be used to identify adverse events after surgery.

**Classify:** To work with medical/chemical codes directly might be cumbersome, since those classifications tend to be massive and therefore hard to interpret. It is thus common to use data aggregation as proposed by some classification or combined index from the literature. This could be the *Charlson* or *Elixhauser* comorbidity indices based on ICD-codes, or the *RxRisk V* classification based on ATC-codes. Each of those tools appear with different code versions (ICD-8, ICD-9, ICD-9-CM, ICD-10, ICD-10-CA, ICD-10-SE, ICD-10-CM et cetera) and with different codes recognized as relevant comorbidities (the Charlson index proposed by Charlson et al, Deyo et al, Romano et al. Quan et al. et cetera). Using a third object (in addition to the primary and secondary patient data sets) helps to formalize and structure the use of such classifications. This is implemented in the `coder` package by `classcodes` objects based on regular expressions (often with several alternative versions). Those `classcodes` objects could be prepared by the user, although a number of default `classcodes` are also included in the package (table below).

**Index:** Now, instead of working with tens of thousands of individual ICD-codes, each patient might be recognized to have none or some familiar comorbidity such as hypertension, cancer or dementia. This granularity might be too fine-grained still, wherefore an even simpler index score might be searched for. Such scores/indices/weighted sums have been proposed as well and exist in many versions for each of the standard classifications. Some are simple counts, some are weighted sums, and some accounts for some inherited hierarchy (such that ICD-codes for diabetes with and without complications might be recognized in the same patient, although the un-complicated version might be masked by the complicated version in the index).

**Conditions:** Some further complexity might appear if some codes are only supposed to be recognized based on certain conditions. Patients with THA for example might have an adverse event after surgery if a certain ICD-code is recorded as the main diagnose at a later hospital visit, although the same code could be ignored if recorded only as a secondary diagnosis.

**To summarize:** The coder package takes three objects: (1) a data frame/table/tibble with id and possible dates from a primary source; (2) coded data from a secondary source with the same id and possibly different dates and; (3) a `classcodes` object, either a default one from the package, or as specified by the user. The outcome is then: (i) codes associated with each element from (1) identified from (2), possibly limited to a relevant time window; (ii) a broader categorization of the relevant codes as prescribed by (3), and; (iii) a summarized index score based on the relevant categories from (3).

(i-iii) corresponds to the output from functions `codify()`, `classify()` and `index()`, which could be chained explicitly as `codify() %>% classify() %>% index()`, or implicitly by the `categorize()` function.

## Usage

Assume we have some patients with surgery at specified dates:

```{r}
library(coder)
ex_people
```

Those patients (among others) were also recorded in a national patient register with date of hospital admissions and diagnoses codes coded by the International Classification of Diseases (ICD) version 10:

```{r}
ex_icd10
```

Using those two data sets, as well as a classification scheme (`classcodes` object; see below), we can easily identify all Charlson comorbidities for each patient:

```{r}
ch <- 
  categorize(
    ex_people,                  # patients of interest 
    codedata = ex_icd10,        # Medical codes from national patient register
    cc = "charlson",            # Calculate Charlson comorbidity
    id = "name", code = "icd10" # Specify column names
  )

ch
```

How many patients were diagnosed with malignancy?

```{r}
sum(ch$malignancy)
```

What is the distribution of the combined comorbidity index for each patient?

```{r}
barplot(table(ch$charlson))
```

There are many versions of the Charlson comorbidity index, which might be controlled by the `index` argument. We might also be interested only in diagnoses from 90 days before surgery as specified with an argument list `codify_args`as passed to `codify()`:

```{r}
ch <- 
  categorize(
    ex_people, codedata = ex_icd10, cc = "charlson", id = "name", code = "icd10",
    
    # Additional arguments
    index       = c("quan_original", "quan_updated"), # Indices
    codify_args = list(
      date      = "surgery",   # Name of column with index dates
      code_date = "admission", # Name of column with code dates
      days      = c(-90, -1)   # Time window
    )
  )
```

Number of malignancies during this period?

```{r}
sum(ch$malignancy, na.rm = TRUE)
```

Distribution of the index as proposed by Quan et al 2011 during the 90 day period:

```{r}
barplot(table(ch$quan_updated))
```

## Classification schemes

Classification schemes (`classcodes` objects, see `vignette("classcodes")`) are based on regular expressions for computational speed (see `vignette("Interpret_regular_expressions")`), but their content can be summarized and visualized for clarity. Arbitrary `classcodes` objects can also be specified by the user.

The package includes default `classcodes` for medical patient data based on the international classification of diseases version 8, 9 and 10 (ICD-8/9/10), as well as the Anatomical Therapeutic Chemical Classification System (ATC) for medical prescription data.

Default `classcades` are listed in the table. Each classification (classcodes column) can be based on several code systems (regex column) and have several alternative weighted indices (indices column). Those might be combined freely.

```{r}
coder::all_classcodes()
```

# Relation to other packages

`coder` uses `data.table` as a backend to increase computational speed for large datasets. There are some R packages with a narrow focus on Charlson and Elixhauser co-morbidity based on ICD-codes ([icd](https://CRAN.R-project.org/package=icd), [comorbidity](https://CRAN.R-project.org/package=comorbidity), [medicalrisk](https://CRAN.R-project.org/package=medicalrisk), [comorbidities.icd10](https://github.com/gforge/comorbidities.icd10), [icdcoder](https://github.com/wtcooper/icdcoder)). The `coder` package includes similar functionalities but has a wider scope.

# Code of conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
---
title: "Review"
author: "Erik Bülow"
date: "11/20/2020"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Thank you!

First of all, I would like to express my gratitude, both for the impressive work performed for the review, but also for letting me take the time to address the raised concerns! I have now (finally) been able to (hopefully) make all the changes and improvements as requested. I will address each comment in a point-by point manner below!

# Zabore

## Review Comments

> I have had to classify patients based on data that contained a variety of types of ICD codes in the past and it was incredibly cumbersome, so a package that automates this process for some standard use cases and with customization possible is a great contribution!

Thank you very much for your kind words and appreciation!

### README

> The statement of need could be fleshed out more. As a biostatistician working in the medical field, I was able to get a sense of what the package is meant to do since I'm familiar with ICD codes and co-morbidity indices such as Charlson. But some more details about the problem this solves, and who it is meant to help, could be usefully added. I believe that the package can take a dataset where patients have various codes, such as ICD codes, with associated dates, and it then uses regular expressions to identify specific codes of interest for creating, for example, the Charlson comorbidity index. Is that right? I understand that you are trying to be vague since the package could be used generally for any application, but I find that in not being specific about the medical application at all times it is harder to understand the core functionality. A key to clarity of purpose may lie in terminology, and clearly defining the different pieces of information this package can bring together.

**A:** I agree and have expanded this section of the README.

> There are currently some errors printed in the output on the README. Those errors were not replicated when I ran the code locally. But the README needs to be updated without errors.

**A:** Thanks for noticing! This has now been fixed in later versions!

### coder vignette

> I would have preferred that the trivial example from the vignette coder actually be a simple example of classifying patients into groups based on ICD codes, rather than an unrelated example based on car makes and models. You could also show in this vignette how the ex_carbrands classcodes object was created. I found myself wondering that as I read the vignette.

**A:** I agree! The vignettes have been extensively rewritten and reorganized. Examples with car data have been removed. `vignette("coder")` now includes examples based on the RxRisk V classification. Additional examples with Charlson and Elixhauser comorbidity are also included in other vignettes.

> Under the "Interpretation" section, the listed names are mostly not in line with the actual output printed above. I assume you changed the contents of the example datasets and did not change the corresponding vignette text so this should be updated.

**A:** The vignette has been completely rewritten wherefore I hope this is no longer an issue.

> What was the index based on in this example? i.e. how did the function categorize() know that the index should be a sum of the matched cars owned by each person? In the documentation for the function categorize it says that a value of NULL would lead to the use of all available indices, but what does that mean exactly? Available from where?

**A:** The vignette is rewritten and the documentation for the `index` argument used by `categorize()` has been updated:

    Argument passed to index(). A character vector of names of columns with index weights from the corresponding classcodes object (as supplied by the cc argument). See attr(cc, "indices") for available options. Set to FALSE if no index should be calculated. If NULL, the default, all available indices (from attr(cc, "indices")) are provided.

> In fact, the example on the README page is much more relevant and could easily be expanded to demonstrate the functionality in place of the car brand example. Maybe both are not needed.

**A:** I agree, and have done so!

### classcodes vignette

> I liked the output of the visualize() function, but seeing an example of an ICD code that would be identified by the given regular expression could be helpful to the less experienced reader. Presumably many users are familiar with ICD codes generally, but for the less experienced user this would be valuable.

**A:** I agree! I have introduced a new `vignette("Interpret_regular_expressions")` where I have also presented a new `summary.classcodes()` method as well as two functions `coodebook()`and `codebooks()`. Using those, it is much easier to get an overview of all codes being recognized by the relevant regular expressions.

> I found the classcodes vignette a little hard to follow and found myself wondering throughout, "when would I use this?" and wishing for a concrete example that was used throughout demonstrating: this is what I want to do, and here is how I use these functions to do it.

A: I have re-organized all vignettes and have tried to include some concrete use cases.

## Comorbidity vignette

> The link to the Elixhauser co-morbidity index is broken

A: Thank you! I have changed the link to a formal citation of the relevant paper.

> In the Risk Rx V example you mention creating dataframe ex_atc in which patients could have zero, one, or several codes prescribed. However, the way you create the data each patient in fact has at least one code. Either alter the example to match (preferable to demonstrate how patients without codes are handled) or alter the text.

**A:** Good point! I have updated to only include 90 out of the 100 patients. This (as well as the other example data sets) are now presented in `vignette("ex_data")`. The data set `ex_atc` has also been included in the package. The code to generate the sample is therefore no longer presented (but is available in `data-raw/example_data.R`

> You mention surgery dates farther down in the vignette, but at the beginning only refer to a vague event date. This example would be improved by making it more specific. Who are the patients, what are the dates, what are we trying to do, and why?

**A:** I agree and have changed "event" to "surgery".

> I noticed that when the default categorization was used, a nice message was printed in the R console stating what the classification was based on. This is a great use of messages. However, when a different classification was specified using cc_args this message was no longer printed. This could be standardized as I think the messaging is nice, and especially important when the user is changing it from some default.

**A:** You are correct that the message is only displayed if relying on the default setting. This was inspired by the `join`-functions from the dplyr package, as described for the `by` argument in the manual:

    `If NULL, the default, *_join() will perform a natural join, using all variables in common across x and y. A message lists the variables so that you can check they're correct; suppress the message by supplying by explicitly.`.

This design pattern is also described and recommended in chap 15 of the Tidyverse design manual ([https://design.tidyverse.org/def-inform.html).](https://design.tidyverse.org/def-inform.html).)

I like this approach since I think the message is less needed if the user explicitly specified the setting. I have nevertheless borrowed the text above for the manual to be more transparent considering the underlying motif for this behavior.

### Documentation/Examples

> Both the DESCRIPTION file and coder help file link only to the GitHub repo, even though this package does have a website. Make the link to the website easier to find since it's the best place for new users to start.

**A:** Thank you for noticing. This has been updated!

> Good use of multiple points of entry. I found descriptions of the package through the help file, the README, the website.

**A: Thank you :-)**

> The help file for classcodes mentions the "triad of case data, code data, and classification data." I found this to be a real "aha!" moment for understanding your package but this was the first place I saw this particular description. I think you could usefully focus more on this three-piece idea throughout the README, help files, and vignettes, with consistent terminology used throughout. Also, perhaps function and argument names could be altered/unified to reflect the trinity of information sources. I didn't find the names particularly intuitive as they are and was constantly having to reference them again when I was testing the functions.

**A:** Thanks for the suggestion! This is now (I hope) better explained in the package README, DESCRIPTION, main vignette and in the manual page for `categorize()`.

> I'd like to see examples that utilize every option available in each function, rather than only simple examples. I often find myself wanting to use a non-default option in a function and not finding any example that helps me understand how to specify the option. Having simple examples is essential, but also including a more complex example that utilizes more than just the basic function arguments will really help users with more complex problems or more advanced users. For example, the example given in the classcodes help file isn't very useful. The object you pass to as.classcodes() in the example is already a classcodes object. How would I use the argument hierarchy? I don't know because you haven't given an example.

**A:** I agree and have added examples which now should cover the use of all arguments etc.

> Many Roxygen comments went beyond my 80 character line, making readability difficult. You could use carriage returns to improve readability of comments.

**A:** I agree and have changed accordingly.

> Code is formatted well, and in line with tidy style guidelines.

**A:** Thank you very much!

> First, the documentation for filter_dates should note that there is a default lower date limit. I had to go to the documentation of dates_within before I could understand the default behavior of filter_dates. Second, why is the default lower limit set to 1970-01-01? Does that have something to do with the coding schemes available within the package?

**A:** Those functions have been removed due to suggestions from @drgtwo (see below).

> No example was included in the help file for set_classcodes No example was included in the help file for as.keyvalue.classcodes

**A:** All examples have been rewritten (see comment above). They should now cover all functions and arguments.

> I see that codify throws an error if the values specified in the date argument aren't of class "Date". However I didn't see anywhere in the documentation information for the user specifying that dates need to be formatted specially. Typically when I read clinical data into R I end up with date columns that are formatted as character and so I would need to convert them to use these functions, so it might be worthwhile to include this detail somewhere.

**A:** Good point! This has now been specified!

### Tests

> There's an issue with the test-as.keyvalue file. It is currently saved as "test-as.keyvalueR" rather than "test-as.keyvalue.R". Also, this test file has an error. I believe instead of a \<- as.keyvalue.classcodes(elixhauser, "icd10se") you should have a \<- decoder::as.keyvalue(elixhauser, "icd10se"). And similarly within the test.

**A:** Thank you for noticing! I have changed the file name and the code!

> There are at least two tests named "multiplication works" that are unrelated to multiplication: test-all_classcodes, test-as.keyvalue

**A:** Thank you! I have now removed all `context()`-lines from the test files in accordance with the latest version of the testthat package.

### Functionality/performance

> The use of 365 days to indicate a year is not exact. Depending on the year, some dates could be missed due to leap years. I see that using windows in units of days is an easy coding solution here, but perhaps some more complex options for when time windows are specified in months or years, rather than days, could be supplied as well.

**A:** I agree in theory, although I think this is a rather complex issue. Rules for leap years are different depending on whether the year is divisible by 4, 100 and 400 for example. Adding leap seconds makes it even more difficult. Also the week numbers are treated differently in different parts of the world (a year might have 53 weeks in some countries but only 52 in others). Further complexity might consider different time zones (including some regions with 30 minutes intervals, summer times etc). The length of the year is also different in different religions (although the christian/Gregorian calender might of course be the norm :-). Altogether, I think more complicated time windows are better handled by dedicated date packages, which could be used explicitly by the user. A better approximation for the number of days of a year might be 365.241 (as advocated by some), but in this setting it does not really matter since time points were only recorded as dates (not time of day). I have therefore left this example as is.

> I thought it was strange that codify did not return the same variable name for the date of interest as specified in data. Using your example data, ex_people has a date called event but when I combine it with ex_icd10 using codify the resulting data frame has a date variable named date rather than event. I would prefer to maintain the original variable name. This behavior does not occur in other functions; for example, categorize returns the original date variable name even when arguments are passed to codify.

**A:** I agree! This was a bug which is now fixed!

> I expected categorize to be able to take the results of a call to codify as input. This could be a useful alternative specification option, and is being done internally anyhow

**A:** This is a good point, wherefore `categorize` has now been made S3-generic with a dedicated method for objects of class `codified` (a new class introduced for the output from `codify()`)

> I didn't have a large dataset of the type that might be used in this setting to realistically test the claim that the performance is fast.

**A:** Thank you nevertheless for the throughout review regarding all other aspects :-)

### Package check notes:

> I got the following notes about example time and file size when I ran check on my local Windows machine:

         installed size is 37.0Mb
         sub-directories of 1Mb or more:
           Meta   8.0Mb
           R      3.0Mb
           data   3.0Mb
           doc   10.0Mb
           help   5.0Mb
           html   2.0Mb
       Examples with CPU (user + system) or elapsed time > 5s
                  user system elapsed
       codebooks  1.39   2.53   57.08
       codebook   1.39   2.41   59.13
       categorize 0.23   1.05   22.28
       classify   0.13   0.28    7.01

**A:** I was not able to reproduce this, neither locally, nor by the GitHub actions set up for Windows, Mac or Ubunto. When I check the file sizes of the whole repo (excluding Git-related files) I have approximately 6 MB in total. The R folder for example is 94 KB. Wen I build the package locally I get a 349 KB source archive. I therefore leave this as is for now.

# drgtwo

> `coder` takes a set of patients, links them to diagnosis codes, and links those to healthcare indices. From my brief experience in healthcare I know this is an important use case with many applications. While other packages (which the author cites perfectly in the JOSS manuscript) include tools for calculating , the ability to match many diagnosis codes at once, and to calculate multiple indices in the same way, make the coder package an important contribution to the field. Great job!
>
> The package follows the ROpenSci guidelines, is thoroughly tested and has chosen a high-performance approach. I've installed and checked it, and outside of the one bug described below (duplicate names to `codify()`) I found it to work for me on use cases.
>
> I have some recommendations for specific changes below, but I'm of course I'm open to dialogue.

**A:** Thank you very much! I really appreciate your throughout review and all the time and energy invested in it!

#### **Documentation on a "real" use case**

> One thing that's missing from the documentation (for functions, the README, and the vignettes) is an example of what someone would do with the *output* of the package. The example in the "coder" vignette isn't healthcare-related, instead focused on a toy example involving patient cars. Perhaps replace it with a vignette that works through a "real" example, and shows at least one result from it that would be similar to what your next typical step would be. Besides being more interesting, it serves as documentation for the output: if you created a histogram of Charlson indices, it brings the user's attention to the importance of that column.

**A:** I agree! I have removed the toy example with car brands and have tried to include more realistic examples in the README and vignettes. Some histograms (or at least bar plots) with Charlson (and Rx Risk V) indices are also provided.

> You could add something from the result of your `ex_people` and `ex_icd10`, but that join has only a single positive result (one patient with peripheral vascular disease). Of course real healthcare data is necessarily private, but you could instead consider taking a small sample of rows and columns from the [SynPUF](https://www.cms.gov/Research-Statistics-Data-and-Systems/Downloadable-Public-Use-Files/SynPUFs/DE_Syn_PUF) data, which is synthetic emergency room data including admissions and diagnoses. After using coder to determine what diagnoses occurred after an emergency room visit, the vignette could get one or two results from the data (e.g. showing the average Charlson comorbidity index of emergency room admissions over time, creating a histogram of them, or showing the most common diagnoses within the window). This would help communicate why it's helpful to annotate a dataset in this way.
>
> (You don't need to use a sample from SynPUF if you don't want to; you could also just construct the simulated data to have more joined diagnoses).

**A:** Thank you for really good suggestions! Unfortunately, I was not able to include the SynPUF data but I did find a previsously exported data set `icd.data::uranium_pathology` which I have now modified as part of the internal example data sets `ex_icd10` and `ex_people`. I think those data sets should be more realistic.

> The README is solid, but does jumps immediately into examples of simulated data, without discussing why someone might want to join diagnoses in the previous year. This doesn't have to take much text; it could be a short bulleted list of use cases like "Discovering adverse events after surgery"/"Determining comorbidities before clinical trials." It's likely that users already know their use case, but an example lets them recognize it and think "this package is for me!"

**A:** I agree and have substantially expanded the README. The suggested sentences are now also included under the "Typical use case" section.

> A second piece of advice with the README is to start with an example that doesn't include a date column before showing the one that does, to ease the user into the use of the `categorize` function with relatively few arguments.

**A:** Thanks for the suggestion, which has now been included in the README!

> The comorbidities vignette and the JOSS paper are very well done in terms of giving the appropriate level of background and describing the use case.

**A:** Thank you very much! I appreciate it!

#### **Duplicate names to codify**

> If there are duplicate names in the `data` passed to `codify()`, it returns a data.table error that isn't informative toward fixing the problem. (`categorize()` does catch this with "Non-unique ids!" but not `codify()`).

    people_doubled <- rbind(ex_people, ex_people)
    codify(people_doubled, ex_icd10, id = "name", date = "event", days = c(-365, 0))

A: Thank you for noticing! The message should now be the same as for `categorize()`.

> More importantly, don't there exist use cases for `categorize` where there are multiple events for the same patient, with different dates? Examples could include adverse events after starting multiple lines of therapy, or comorbidities before multiple diagnoses. In those cases, doesn't it make sense to return one row for each event, even if there are multiple for a patient? Should the check only error out when there are duplicate name/date **pairs**?

**A:** I agree that such a feature is relevant. The problem, however, is that unit data is matched to code data based on the index variable and that I cannot perform such matching based on the date column (which would be a non-equi-join, as allowed for some `data.table` operations but not in `merge` which is currently used). Although this would be possible after some factoring of internal functions, I think it is currently better to perform such operations using standard functionality outside the package, such as with `x %>% group_by(y) %>% codify(...)` for `dplyr` or `x[, codify(...), by = y]` with `data.table`.

#### **as.codedata**

> I think the `as.codedata()` approach can be improved to make the package more understandable and usable. Some issues:
>
> -   By convention `as.X` functions in R return an object of class X, but this returns a data.table.
>
> -   `codify()` describes the second argument as "output from as.codedata", but the function still works if given a data frame, data.table, or tibble.
>
> -   By default, `as.codedata()` filters out dates in the future and dates before 1970. I assume this is meant to remove bad data, but isn't it better to leave such data quality filters to the user? As it is, the user must go through a few pages of documentation (codify/categorize -\> as.codedata -\> dates_within) to learn about this behavior. And in any case where there's a date window, the extreme date values won't affect the coding anyway.
>
> It looks to me like the main reason for `as.codedata()` is to speed up the function by making it a data.table and setting keys. But you could do this within `codify()` as well; the only advantage this provides is if you run many codings with different ids/dates (or different arguments) while keeping the code data the same. I've done some benchmarking and it looks like the improvements become visible (in tens of milliseconds) when there are around million coded events.
>
> Do we expect it to be common for users to run the package with millions of coding events, where the codedata stays the same while the input events change, and in an environment where fractions of a second matter? Is this common enough to be worth imposing extra instructions on every user of the package?

> My recommendation is
>
> -   Don't export `as.codedata`, and instead do the preprocessing/checking of codedata within the codify function instead of suggesting that the user use `as.codedata`

A: I agree and have now integrated `as.codedata()` into `codify()` as suggested!

> -   In the documentation for `codify` and `classify`, as well as your documentation examples, describe the codedata input as "a table with columns `id`, `code`, and optionally `date`."

**A:** The documentation has been improved. I have changed the requirements that the codedata must have columns named `id`, `code` and `code_date`. Instead, the names of corresponding columns should now be specified as arguments in relevant functions.

> -   If you're very confident that performance when keeping the codedata the same and trying many different datasets is important, you could add a function called `classcodes_prepare` or `prepare_classcodes` that does the conversion to data.table and sets the indices, and could describe that in the `details` section of the documentation. But I'd want to understand why that's a typical use case.

**A:** See previous comment (`as.codedata()` removed in favor of a re-factored `codify()`).

> Relatedly, I recommend removing (or at least making internal) `dates_within()` and `filter_dates()`. Their purpose (applying a filter on dates with some defaults) has no relationship to the rest of the package, and is something the user can do themselves with tools they're accustomed to (base R, data.table, or dplyr).

**A:** I agree and have removed those functions. They were previously included since I was unaware of `data.table::IDate` and that `data.table::between` is optimized for date comparisons. The problem with normal date comparisons is that those can be really slow for big data.

> #### **regex\_ in column names when tech_names = TRUE**
>
> The output of `categorize()` on a table returns columns with spaces in their names. This isn't well set up for additional analysis, since it makes it difficult to do any kind of programming with them, including using data.table to filter for one diagnosis or to aggregate the percentage of patients (perhaps within each group) that have a condition. It's nice for displaying the names in a table, but is it a common use case to display individual patients in a table (as opposed to aggregated statistics?)
>
> It seems like the `tech_names` argument is designed to fix this, but it leaves prefixes like `charlson_regex_` on every column name, which will need to be removed for meaningful downstream analysis. How about removing the `charlson_regex_`, or at least the `regex_`, in these cases? (Indeed, is there a reason that the `charlson` classcodes object itself has to have the `regex_` prefixes? It already has an attribute `regexprs` that includes those column names). Besides which, perhaps consider leaving `tech_names` to default to TRUE for the reasons described above.

**A:** Good point! I have made several changes:

-   `classcodes` object no longer have column prefixes `(reg|ind)ex_`.

-   I have introduced a new `print.classcodes()` method for a better default display of classcodes where regex and indices are identified by a heading and not by column names prefixes

-   `categorize()` has a new argument `check.names` (same as `data.frame`/`data.table`). This argument is `TRUE` by default, making the column names syntactically correct (using dots instead of spaces). The original names (possibly with spaces) are received by `check.names = FALSE`, which might sometimes be useful.

The reason for the long names implied by `tech_names` is that `categorize` is sometimes used multiple times, for example to enhance a data set with both comorbidity and adverse events. To group such variable names by common and descriptive prefixes might then be useful (due to tab completion etc). But the `(reg|ind)ex_` part is no longer included and the need to use those longer names has probably decreased by the use of the `check.names` argument.

> #### **tibbles and data.tables**
>
> Your examples like `ex_people` are tibbles, but when `categorize()` or `codify()` is passed a tibble, it returns a data.table. This would be a surprising behavior for people using these packages within a tidyverse workflow. I think data.table is a terrific package, but there's not a reason to surprise users with the data type if they're not accustomed to it. (And the fact that the example datasets are tibbles rather than data.frames or data.tables adds to the inconsistency a bit).
>
> I recommend ending the functions with something like
>
>     # Where data was the argument passed in, and ret is what's about to be returned
>     if (tibble::is_tibble(data)) {
>       ret <- tibble::as_tibble(ret)
>     }
>
> This would mean that it returns a data.table when it's passed a data.frame or data.table, but a tibble if and only if it's passed a tibble. Admittedly, this requires adding an import for tibble (which perhaps is why it wasn't done), but since tibble is imported by 800 CRAN packages (including dplyr + ggplot2, each depended on by \~2000 packages) it's a fairly low-impact dependency. This also doesn't strike me as a utility package that will frequently be installed in production systems; it's a scientific package that would typically used with other data analysis tools. I think there are some useful thoughts on tibble dependencies [here](https://www.tidyverse.org/blog/2019/05/itdepends/).

**A:** I agree and have made the following changes:

-   included tibble as a dependency

-   `all_classcodes()` now returns a tibble.

-   `categorize()` has been re-factored into an S3-generic and returns data sets of the same class as the input (data.table/data.frame/tibble).

-   `as.classcodes()` returns as tibble with an additional class attribute.

-   `codify()` is often used to return large data sets (several millions of rows) which should be used in a following step by `classify()`. I therefore think that the `data.table` format is preferred here. I have re-factored the function into S3-methods, however, to treat input as data.frame/data.table/tibble in a better way. I have also implemented a `print.codified()` method, which prints the first `n` rows as a tibble (possible to override with `print(..., n = NULL)`, which will print the object as is [a data.table]). I think this might be a good compromise for most users. I have also clearly stated that the preview is simply a preview as part of the output from the print method.

-   `classify()` returns a matrix for efficiency, since this object is always logical/boolean. I think this should be kept as is. The two methods `as.data.frame.classified()` and `as.data.table.classified()` should make it relatively simple to convert the output if desired (including to tibbles inhereting from data.frame).

-   `summary.classcodes()` returns a tibble, which is also printed as such through `print.summary.classcodes()`.

> Relatedly (though less important), the example datasets don't print as tibbles by default. If you follow the instructions in `usethis::use_tibble()`, you could support printing it as a tibble even when the tibble/dplyr packages aren't loaded. The additional advantage of this is that you could get rid of most of the uses of `head()` in the README, making your examples more concise and focused on your use case.

**A: I agree and have changed accordingly!**

> #### **Naming**
>
> -   `index` and especially `visualize` are very generic names for very specific functions, and doesn't give any hints about what they're used for. How about `visualize_classcodes`?

**A:** The name `visualize` was actually a half-baked S3-method for `generics::visualize()` with the aim to "Visualize a data set or object". It has now been changed into an actual S3-method to fit with this generic.

I am less confident that the `index()` function is that specific since its meaning will differ depending on different `classcodes` object. The terminology is borrowed from the typical use case of calculating a comorbidity index, which was the initial motif for developing the package ( "Charlson comorbidity index", etc). Similar verbs would be for example "aggregate", "count" or "summarize" but those would both be uninformative and collide with other common functions. `index()` on the other hand is not used in any widely used package as from what I have seen. I guess `calculate_index()` or similar would work as well but it just seems longer without any convincing reason :-)

> -   An alternative for function naming is to have a common prefix for functions, e.g., `coder_classify`, `coder_categorize`, `coder_index`, `coder_codify`, `coder_visualize`. This has both the advantage of ensuring it doesn't overlap with other packages and making it easy to find codify-related functions with autocomplete. But that's just a suggestion.

**A:** Thank you for the suggestion! I agree that this convention can be useful. I do think it is almost as simple to use autocomplete with the `coder::` prefix instead of `coder_`, however. The difference is just one additional key stroke :-) and I would prefer shorter names without this redundancy `coder::coder_xxx`.\
This naming convention is also discussed in the Tidyverse design guide, section 5.2 (<https://design.tidyverse.org/function-names.html#function-families>) with the following paragraph:

    Not sure about common prefixes for a package. Works well for stringr (esp. with stringi), forcats, xml2, and rvest. But there’s only a limited number of short prefixes and I think it would break down if every package did it.

> I agree with Noam that coder isn't an ideal package name, if only because it makes the online resources a bit harder for users to find. Try Googling "coder", "R coder", or "R coder github"! But if it's too late to change the name, I don't consider it a dealbreaker.

**A:** Thank yo for the suggestion. I do understand the possible downsides of the name and I partially agree. On the other hand, try Googling the "sister" package "CRAN decoder", which does find the `decoder` package easily. Consider also some popular RStudio packages such as: `generics`, `baguette`, `glue`, `haven`, `parsnip`, `tune`, `reciepies`, `dials`, `workflows` and `yardstick`, none of which are very Googlable on their own. I thus hope that a later CRAN release will make Googling simpler :-)
---
title: "coder"
author: "Erik Bulow"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    df_print: tibble
vignette: >
  %\VignetteIndexEntry{coder}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: references.bib
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The `coder` package simplifies unit classification based on external code data. this is a generic aim that might be hard to grasp without further concretization. In this vignette, I will first explain the overall design principles, and then exemplify the concept with a typical use case involving patients with total hip arthroplasty (THA) and their pre-surgery comorbidity. Note, however, that the package is not limited to patient data or medical settings.

```{r}
library(coder)
```

# Triad of objects

Functions of the package relies on a triad of objects:

1.  Case data with unit id:s and possible dates of interest
2.  External code data for corresponding units in (1) and with optional dates of interest and
3.  A classification scheme ('classcodes' object) with regular expressions to identify and categorize relevant codes from (2).

It is easy to introduce new classification schemes ('classcodes' objects) or to use default schemes included in the package (see `vignette("classcodes")`).

# Triad of functions

There are three important functions to control the intended work flow of the package:

i.  `codify()` will merge object (1) and (2) for a coded data set of the intended format. If optional dates are specified, those will be used to construct time windows in order to filter out only the important dates (i.e. comorbidity during one year before surgery or adverse events 90 days after).
ii. `classify()` will then use the coded data and classify it using the `classcodes` object (3) (i.e. to code comorbidity data by the Charlson or Elixhauser comorbidity classifications).
iii. `index()` is a third optional step to summarize the individual `classcodes` categories to a (possibly weighted) index sum for each coded item (i.e. to calculate the Charlson comorbidity index for each patient).

Those steps could be performed explicitly as `codify() %>% classify() %>% index()` or implicitly by the main function `categorize()` combining all steps automatically.

# Use case

A typical use case of the `coder` package would consider patient data and comorbidity as described in the package [readme](https://docs.ropensci.org/coder/).

The concept of comorbidity is often attributed to Feinstein [-@Feinstein1970]:

> [T]he term co-morbidity will refer to any distinct additional clinical entity that has existed or that may occur during the clinical course of a patient who has the index disease under study.

Let's consider a group of patients with THA, as identified from a national quality register, which might be large in size. Assume we are interested in those patients' pre-surgery comorbidity, which is not captured by the quality register itself. Instead, this data might be codified in a secondary source, such as a national patient register containing all hospital visits and admissions during several years, both before and after the THA-surgery. Each hospital visit/admission might be recorded with one or several medical codes, for example using the International classification of diseases version 10 (ICD-10). Similarly, a medical prescription register might hold records of prescribed drugs with their corresponding codes from the Anatomic therapeutic chemical classification (ATC) system.

Thus, combining the primary and secondary data sets (objects 1-2 above) using some unique patient id, and a possible time window (i.e. to only consider comorbidity as recorded during one year before the THA), is a first step to identify patient comorbidity. This step is performed by the `codify()` function in step (i) above.

We have now gathered all the relevant codes for each patient. Common classifications (i.e. ICD-10 and ATC) are wast, however, including tens of thousands of medical/chemical codes, which are cumbersome and impractical to use directly. It is therefore common to categorize such codes into broader categories (i.e. by the Charlson, Elixhauser or RxRisk V classifications as below). Such classification could be a simple code matching problem using a look-up table. This is generally a slow, cumbersome and error-prone process, however. I therefore recommend to use regular expression for a compact code representation, as well as a computationally faster procedure. This is implemented in the `classify()` function from step (ii) above.

We have now reduced the data from tens of thousands of codes to perhaps 10-50 combined categories. This might be sufficient in some cases, although further simplifications might also be needed. It is thus common to simplify comorbidity into a single number, an index score, as the sum of individual comorbidities, possible weighted to differentiate more serious conditions from more trivial. Different weights might be of relevance under different circumstances or in different fields. This is implemented by the `index()` function in step (iii) above.

# Charlson and Elixhauser

The Charlson [-@Charlson1987] and Elixhauser [-@Elixhauser1998] comorbidity indices are two examples used in medical research. Each index consist of several medical conditions, possibly summarized by a (weighted) index. Each condition is defined by a set of medical codes [@Quan2005]. Different versions of the International Classification of Diseases (ICD) codes are often used.

The `coder` package provides substantial functionality for both Charlson and Elixhauser, although we will not focus on those indices here (but see examples in `vignette("classcodes")`). Several other R packages have functions for Charlson and Elixhauser:

-   [icd (CRAN)](https://CRAN.R-project.org/package=icd)
-   [comorbidity (CRAN)](https://CRAN.R-project.org/package=comorbidity)
-   [medicalrisk (CRAN)](https://CRAN.R-project.org/package=medicalrisk)
-   [comorbidities.icd10 (GitHub)](https://github.com/gforge/comorbidities.icd10)
-   [icdcoder (GitHub)](https://github.com/wtcooper/icdcoder)

`icd` and `comorbidity` are both good packages well suited for their purpose based on effective implementations. `medicalrisk` can be used with ICD-9-CM codes but is not up-to-date with the latest version of ICD-10. `comorbidities.icd10` and `icdcoder` are not actively developed or maintained.

One advantage with the `coder` package is the great flexibility for combining different sets of codes (ICD-8, ICD-9, ICD-9-CM and ICD-10 et cetera), with different weighted indices.

# Risk Rx V

Another advantage of the `coder` package is the inclusion of additional classifications (see `?all_classcodes()`), such as the pharmacy-based case-mix instrument Rx Risk V [@Sloan2003]. We will use this classification in an example. This classification, in contrast to Charlson and Elixhauser, relies on medical prescription data codified by the Anatomic Therapeutic Chemical classification system (ATC).

As for all classcodes objects in the package, additional information and references are found in the object documentation (`?rxriskv`).

# Concrete example

Let us consider the hypothetical setting above using some example data (`ex_peopple` and `ex_atc`) as described in `vignette("ex_data")`.

# Default categorization

A first attempt to calculate the Rx Risk V score for each patient:

```{r}
default <- categorize(
    ex_people, codedata = ex_atc, cc = rxriskv, id = "name", code = "atc")
default
```

The first two columns are identical to `ex_people`. Additional columns indicate whether patients had any of the individual comorbidities identified by Rx Risk V. Patients without any medical prescriptions have `NA` values (which might be substituted by `FALSE`). The last columns contain summarized index values (weighted sums of individual comorbidities). Let's summarize the distribution of a weighted index according to `pratt` [@Pratt2018]:

```{r}
hist2 <- function(x) {
  hist(x$pratt, main = NULL, xlab = "RxRisk V", col = "lightblue")
}
hist2(default)
```

# Specified time-window

Some prescriptions might have been filed long before surgery, or even after. Those codes are less relevant for comorbidities present at surgery. We can limit the categorization to a time window of one year (365 days) prior to surgery. This is done internally by the `codify()` function, hence by specifying a list of arguments passed to this function:

```{r}
codify_args <- 
  list(date = "surgery", code_date = "prescription", days = c(-365, -1))

ct <- 
  categorize(
    ex_people, 
    codedata    = ex_atc, 
    cc          = rxriskv, 
    id          = "name", 
    code        = "atc", 
    codify_args = codify_args
  )
  
hist2(ct)
```

# Alternative classification

Comorbidities are identified from ATC codes captured by regular expression (see `vignette("classcodes")` and `vignette("Intrpret_regular_expressions")`). Codes identified by `atc_pratt` are used by default. Let's use an alternative version adopted from Caughy [-@Caughey2010] as specified by an argument passed by the `cc_args` argument.

```{r}
hist2(
  categorize(
    ex_people, 
    codedata      = ex_atc, 
    cc            = rxriskv, 
    id            = "name", 
    code          = "atc",
    codify_args   = codify_args,
    cc_args       = list(regex = "caughey")
  )
)
```

# Specified index

We did not specify how to calculate the weighted index sum above, wherefore all available indices were provided by default. We might go back to Pratt's classification scheme (`atc_pratt`) and only calculate the corresponding index `pratt`. Let´s also perform the three computational steps explicitly instead of using the combining `categorize()` function and tabulate the result

```{r}
codify(
  ex_people, 
  ex_atc, 
  id        = "name", 
  code      = "atc",  
  date      = "surgery", 
  code_date = "prescription",
  days      = c(-365, -1)
) %>% 
  classify(rxriskv) %>% 
  index("pratt") %>% 
  table()
```

# Dirty code data

Let's assume that our code data is not as clean as simulated above.

```{r}
s <- function(x) sample(x, 1e3, replace = TRUE)

ex_atc$code <- 
  paste0(
    s(letters), s(0:9), s(letters), s(c(".", "-", "?")), 
    ex_atc$atc, s(letters), s(0:9)
  )

ex_atc

sum(
  categorize(
    ex_people, 
    codedata = ex_atc, 
    cc       = rxriskv, 
    id       = "name",
    code     = "code"
  )$pratt,
  na.rm      = TRUE
)
```

Thus, no codes are recognized (every one got index = 0). By default, codes are only recognized if found immediate in its corresponding column. This can be controlled by arguments `start` and `stop` specified via `cc_args`. We can also ignore all non alphanumeric characters by setting `alnum = TRUE` as passed to `codify()` by argument `codify_args`.

```{r}
hist2(
  categorize(
    ex_people, 
    codedata = ex_atc, 
    cc       = rxriskv, 
    id       = "name",
    code     = "code",
    cc_args  = list(
      start  = FALSE, 
      stop   = FALSE
    ),
    codify_args = list(
      alnum = TRUE
    )
  )
)
```

# Bibliography
---
title: "Example data"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Example data}
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
library(coder)
```

This vignette contains some example data used in the other vignettes.

# Patients

`ex_people` contains 100 patients (with random names from the [`randomNames`](https://centerforassessment.github.io/randomNames/) package) who received total hip arthroplasty (THA) surgery at given (random) dates (`surgery` column). This data represent a sample from a national quality register.

See also `?ex_people`.

```{r}
ex_people
```

# Diagnoses data

We are interested in comorbidity for the patients above and have collected some synthesized diagnostics data (`ex_icd10`) from a national patient register (we can at least assume that for now). Patients have one entry for every combination of recorded diagnoses codes according to the International classification of diseases version 10, `icd10`, and corresponding dates of hospital `admission`s for which those codes were recorded. (Column `hdia` is `TRUE` for main diagnoses and `FALSE` for underlying/less relevant codes).

See also `?ex_icd10`.

```{r}
ex_icd10
```

# Medical data

Assume we have some external code data from a national prescription register. Such register would likely cover additional patients but let's just consider a small sample with ATC codes for patients above, such that each patient can have zero, one, or several codes prescribed at different dates.

```{r}
ex_atc

```
---
title: "Interpret regular expressions"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Interpret regular expressions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: references.bib
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(coder)
```

Classcodes objects (as described in `vignette("classcodes")`) use regular expressions to classify/categorize individual codes into groups (i.e. comorbidity conditions). Those regular expressions might be hard to interpret on their own. Several methods are therefore available to aid such interpretation of the classcodes objects.

# `visualize()`

A graphical representation of a classcodes object is created by `visualize()`. It will be showed in the default web browser (requires an Internet connection; not available within this vignette).

```{r, eval = FALSE}
visualize(charlson)
```

Visualization of all groups (comorbidity conditions) simultaneously might lead to complex figures. We can focus on a specific group (comorbidity) by the `group` argument. How is `r charlson$group[1]` codified by `regex_icd9cm_deyo`?

```{r, eval = FALSE}
visualize(charlson, "myocardial infarction", regex    = "icd9cm_deyo")
```

```{r, echo = FALSE}
knitr::include_graphics("regexp_charlson_ci_icd9.png")

```

Hence, all ICD-9 codes starting with `41` followed by either `0` or `2` will be recognized as myocardial infarction according to `icd9cm_deyo`. The corresponding regular expression for ICD-10 is:

```{r, eval = FALSE}
visualize(charlson, "myocardial infarction", regex = "icd10")
```

```{r, echo = FALSE}
knitr::include_graphics("regexp_charlson_ci_icd10.png")

```

Such codes should start with `I2` followed by either `1`, `2` or `52`. The vertical bar `|` (in the regular expression of the heading) indicates a logical "or". See `?regex` for more details on how to use regular expressions in R (Perl-like versions are currently not allowed).

# `summary()`

An alternative representation is to list all relevant codes identified by each regular expression. This is implemented by the `summary()` method for classcodes objects. Note, however, that the regular expressions are stand alone in each classcodes object. Hence, there are no static look-up-tables to map individual codes to each group. We therefore need to specify a code list/dictionary of all possible codes to be recognized by those regular expressions. Then `summary()` will categorize those and display the result. Common code lists are found in the [decoder](https://cancercentrum.bitbucket.io/decoder/) package and are accessed automatically through the `coding` argument to `summary()`. Hence, there is a "keyvalue" object `icd10cm` with all ICD-10-CM codes in `{decoder}:`

```{r}
head(decoder::icd10cm)

```

We can use this code list to identify all codes recognized by `charlson` with its default classification based on "icd10". The printed result (see `?print.summary.classcodes`) is a tibble with each group and a comma separated code list.

```{r}
s <- summary(charlson, coding = "icd10cm")
s
```

A list with all code vectors (to use for programmatic purposes) is also returned (invisible) and accessed by `s$codes_vct`.

Now, compare the result above with the output based on a different code list, namely ICD-10-SE, the Swedish version of ICD-10, instead of ICD-10-CM:

```{r}
summary(charlson, coding = "icd10se")
```

There are some noticeable differences. AIDS/HIV for example has only one code deemed clinically relevant in the USA (thus included in the CM-version of ICD-10), although there are 22 different codes potentially used in the Swedish national patient register. There are additional differences concerning the fifth code position (digits in ICD-10-CM and characters in ICD-10-SE). Those mark national modifications to the original ICD-10 codes, which has only 4 positions (one character and three digits). For this example, the `charlson$icd10` column was based on ICD-10-CM [@Quan2005]. The comparison above thus highlights potential differences when using this classification in a setting based on another classification (such as with data from the Swedish national patient register).

If we are interested in another code version, for example as specified by ICD-9-CM [@Deyo1992] , this can be specified by the `regex`-argument passed by the `cc_args` argument to the `set_classcodes` function. Simultaneously, the `coding` argument is set to `icd9cmd` to match the regular expressions to the disease part of ICD-9-CM classification.

```{r}
summary(
  charlson, coding = "icd9cmd",
  cc_args = list(regex = "icd9cm_deyo")
)
```

# 

# `codebook()`

Even with individual codes summarized, those might still be hard to interpret on their own. The [decoder](https://cancercentrum.bitbucket.io/decoder/) package can help to translate codes to readable names/description. This is facilitated by the `codebook()` function in the `{coder}` package.

The main purpose is to export an Excel-file (if path specified by argument `file`). The output is otherwise a list, including both a summary table (described above) and a tibble with "all_codes" explaining the meaning of each code.

We can compare the codes recognized as AIDS/HIV by either ICD-10-CM or ICD-10-SE:

```{r}

cm <- codebook(charlson, "icd10cm")$all_codes
cm[cm$group == "AIDS/HIV", ]

se <- codebook(charlson, "icd10se")$all_codes
se[se$group == "AIDS/HIV", ]

```

# `codebooks()`

Several codebooks can be combined (exported to a single Excel-file) by the function `codebooks()` (note the plural s). This is difficult to illustrate in a vignette but examples are provided in `?codebooks`

# Bibliography
---
title: "Classcodes"
output: 
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Classcodes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: references.bib
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(coder)
```

# Motivating example

Let's consider some example data (`ex_peopple` and `ex_icd10`) from `vignette("ex_data")`.

Let's categorize those patients by their Charlson comorbidity:

```{r}
categorize(ex_people, codedata = ex_icd10, cc = charlson, id = "name", code = "icd10")
```

Here, `charlson` (as supplied by the `cc` argument) is a "classcodes" object containing a classification scheme. This is the specification of how to match `ex_icd10$icd10` to each condition recognized by the Charlson comorbidity classification. It is based on regular expressions (see `?regex`).

# Default classcodes

There are `r nrow(all_classcodes())` default "classcodes" objects in the package (`classcodes` column below). Each of them might have several versions of regular expressions (column `regex`) and weighted indices (column `indices`):

```{r}
all_classcodes()
```

# classcodes object

Each of those classcodes objects are documented (see for example `?charlson`). Those objects are basically tibbles (data frames) with some additional attributes:

```{r}
charlson
```

Columns have pre-specified names and/or content:

-   `group`: short descriptive names of all groups to classify by (i.e. medical conditions/comorbidities in the Charlson case)
-   `description:` (optional) details describing each group
-   regular expressions identifying each group (see `vignette("Interpret_regular_expressions")` for details and `?charlson` for concrete examples). Multiple versions might be used if combined with different code sets (i.e. ICD-9 versus ICD-10) or as suggested by different sources/authors. (Column names are arbitrary but identified by `attr(., "regexprs")` and specified by argument `regex` in `as.classcodes()`).
-   numeric vectors used as weights when calculating index sums based on all (or a subset of) individual groups. (Column names are arbitrary but identified by `attr(., "indices")` and specified by argument `indices` in `as.classcodes()`.)
-   `condition`: (optional) conditional classification (not used with `charlson` but see example below).

In the example above, we did not specify which version of the regular expressions to use. We see from the printed output above (or by `attr(charlson, "regexprs")`), that the first regular expression is "icd10". This will be used by default. We have ICD-10 codes recorded in our code data set (`ex_icd10$icd10`). We might therefore use either "icd10" or the alternative "icd10_rcs". Other versions might be relevant if the medical data is coded by other codes (such as earlier versions of ICD). We will show below how to alter this setting in practice.

# Hierarchy

Some classcodes objects have an additional class attribute "hierarchy", controlling hierarchical groups where only one of possibly several groups should be used in weighted index sums. The classcodes object for the Elixhauser comorbidity classification has this property:\

```{r}
print(elixhauser, n = 0) # preview 0 rows but present the attributes
```

This means that patients who have both metastatic cancer and solid tumors should be recognized as such if classified. If such patient are assigned an aggregated index score, however, only the largest score is used (in this case for a metastatic cancer as superior to a solid tumor). The same is true for patients diagnosed with both uncomplicated and complicated diabetes.

Consider a patient Alice with some diagnoses:

```{r}
pat <- tibble::tibble(id = "Alice")
diags <- c("C01", "C801", "E1010", "E1021")
decoder::decode(diags, decoder::icd10cm)
```

According to Elixhauser, poor Alice has both a solid tumor and a metastatic cancer, as well as diabetes both with and without complications. The (unweighted) index "sum_all", however will not equal 4 but 2, since metastatic cancer and diabetes with complications subsume solid tumors and diabetes without complications.

```{r}
icd10 <- tibble::tibble(id = "Alice", icd10 = diags)
x <- categorize(pat, codedata = icd10, cc = elixhauser, 
                id = "id", code = "icd10", index = "sum_all", check.names = FALSE)
t(x)
```

# Conditions

Consider Alice once more. Suppose she got a THA and had some surgical procedure codes recorded at hospital visits either before, during or after her index surgery. Those codes are recorded by the Nomesco classification of surgical procedures (also known as KVA codes in Swedish). Here, "post_op" indicates whether the code was recorded after surgery or not. This information is not always accessible by pure date stamps (if so, the approach illustrated in `vignette("coder")` could be used instead).

```{r}

nomesco <- 
  tibble::tibble(
    id      = "Alice",
    kva     = c("AA01", "NFC01"),
    post_op = c(TRUE, FALSE)
  )
```

Thus, the "post_op" column is a Boolean/logical vector with a name recognized from the "condition" column in `hip_ae`, a classcodes object used to identify adverse events after THA (the use of `set_classcodes()` is further explained below and is used here since `hip_ae` includes codes for both ICD and NOMESCO/KVA).

```{r}
set_classcodes(hip_ae, regex = "kva")

```

A code from `nomesco$kva` will only be recognized as an adverse events if 1) the code is matched by the relevant regular expression, and 2) the extra condition (from `nomesco$post_op`) is `TRUE.`

We need to specify that codes are based on regular expressions matching NOMESCO codes. We do this by the `regex` argument passed to `set_classcodes()` by the `cc_args` argument.

In the data set (`nomesco`), "AA01" was recorded after surgery but does not indicate a potential adverse event. "NFC01" is a potential adverse event but was recorded already before surgery. Therefore, no adverse event will be recognized in this case.

```{r}
categorize(pat, codedata = nomesco, cc = hip_ae, id = "id", code = "kva",
           cc_args = list(regex = "kva"))

```

# Use classcodes objects

Most functions do not use the classcodes object themselves, but a modified version passed through `set_classcodes()`. This function can be called directly but is more often invoked by arguments passed by the `cc_args` argument used in other functions (as in the example above).

## Explicit use of `set_classcodes()`

We might use `set_classcodes()` to prepare a classification scheme according to the Charlson comorbidity index based on ICD-8 [@Brusselaers2017]. Assume that such codes might be found in character strings with leading prefixes or in the middle of a more verbatim description. This is controlled by setting the argument `start = FALSE`, meaning that the identified ICD-8 codes do not need to appear in the beginning of the character string. We might assume, however, that there is no more information after the code (as specified by `stop = TRUE`). We can also use some more specific and unique group names as specified by `tech_names`.

```{r}
charlson_icd8 <- 
  set_classcodes(
    "charlson",
    regex      = "icd8_brusselaers", # Version based on ICD-8
    start      = FALSE, # Codes do not have to occur in the beginning of a vector
    stop       = TRUE,  # Code vector must end with the specified codes
    tech_names = TRUE   # Use long but unique and descriptive variable names
  )
```

The resulting object has only one version of regular expressions (`icd8_brusselaers` as specified). Each regular expression is suffixed with `$` (due to `stop = TRUE`). Group names might seem cumbersome but this will help to distinguish column names added by `categorize()` if this function is run repeatedly with different classcodes (i.e. if we calculate both Charlson and Elixhauser indices for the same patients). The original `charlson` object had `r nrow(charlson)` rows, but `charlson_icd8` has only `r nrow(charlson_icd8)`, since not all groups are used in this version.

```{r}
charlson_icd8
```

Note that all index columns remain in the tibble. It is thus possible to combine any categorization with any index, although some combinations might be preferred (such as `regex_icd9cm_deyo` combined with `index_deyo_ramano`).

We can now use `charlson_icd8` for classification:

```{r}
classify(410, charlson_icd8)
```

The ICD-8 code `410`is recognized as (only) myocardial infarction.

## Implicit use of `set_classcodes()`

Instead of pre-specifying the `charlson_icd8`, a similar result is achieved by:

```{r}
classify(
  410,
  "charlson",
  cc_args = list(
    regex      = "icd8_brusselaers", 
    start      = FALSE, 
    stop       = TRUE,
    tech_names = TRUE
  )
)
```

# Bibliography
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{copybig}
\alias{copybig}
\title{Decide if large objects should be copied}
\usage{
copybig(x, .copy = NA)
}
\arguments{
\item{x}{object (potentially of large size)}

\item{.copy}{Should the object be copied internally by \code{\link[data.table:copy]{data.table::copy()}}?
\code{NA} (by default) means that objects smaller than 1 GB are copied.
If the size is larger, the argument must be set explicitly. Set \code{TRUE}
to make copies regardless of object size. This is recommended if enough RAM
is available. If set to \code{FALSE}, calculations might be carried out
but the object will be changed by reference.
IMPORTANT! This might lead to undesired consequences and should only be used
if absolutely necessary!}
}
\value{
Either \code{x} unchanged, or a fresh copy of \code{x}.
}
\description{
Decide if large objects should be copied
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classcodes.R
\name{print.classcodes}
\alias{print.classcodes}
\title{Print classcodes object}
\usage{
\method{print}{classcodes}(x, n = NULL, ...)
}
\arguments{
\item{x}{object of type classcodes}

\item{n}{number of rows to preview (\code{n = 0} is allowed)}

\item{...}{arguments passed to print method for tibble}
}
\description{
Print classcodes object
}
\examples{
# Default printing
elixhauser

# Print attributes data but no data preview
print(elixhauser, n = 0)

# Print all rows
print(elixhauser, n = 31)
}
\seealso{
Other classcodes: 
\code{\link{all_classcodes}()},
\code{\link{as.data.frame.classified}()},
\code{\link{classcodes}},
\code{\link{codebook}()},
\code{\link{print.classified}()},
\code{\link{set_classcodes}()},
\code{\link{summary.classcodes}()},
\code{\link{visualize.classcodes}()}
}
\concept{classcodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_for_datasets.R
\docType{data}
\name{ex_icd10}
\alias{ex_icd10}
\title{Example data for random codes assigned to random people}
\format{
Data frames with 1,000 rows and 4 variables:
\describe{
\item{id}{Random names corresponding to column \code{name} in dataset
\code{ex_people}}
\item{date}{random dates corresponding to registered (comorbidity) codes}
\item{code}{ICD-10 codes from the \code{uranium_pathology}
dataset in the \code{icd.data} package by Jack Wasey originating from the
United States Transuranium and Uranium Registries,
published in the public domain.}
\item{hdia}{boolean marker if corresponding code is the main diagnose of
the hospital visit (randomly assigned to 10 percent of the codes)}
}
}
\source{
https://github.com/jackwasey/icd.data
https://ustur.wsu.edu/about-us/
}
\usage{
ex_icd10
}
\description{
Example data for fictive ICD-10-diagnoses to use for testing and
in examples.
}
\seealso{
Other example data: 
\code{\link{ex_atc}},
\code{\link{ex_people}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_for_datasets.R
\docType{data}
\name{ex_atc}
\alias{ex_atc}
\title{Example data for random ATC codes}
\format{
Data frames with 100 rows and 2 variables:
\describe{
\item{name}{random person names}
\item{atc}{Random codes from the Anatomic Therapeutic Chemical
classification (ATC) system.}
\item{prescription}{random dates of prescription of medications with
corresponding ATC codes}
}
}
\usage{
ex_atc
}
\description{
Example data for fictive people to use for testing and in examples.
}
\seealso{
Other example data: 
\code{\link{ex_icd10}},
\code{\link{ex_people}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_for_datasets.R
\docType{data}
\name{cps}
\alias{cps}
\title{Classcodes for the comorbidity-polypharmacy score (CPS) based on ICD-10 codes}
\format{
A data frame with 2 rows and 2 variables:
\describe{
\item{group}{comorbidity groups, either "ordinary" for most ICD-10-codes or
"special" for codes beginning with "UA", "UB" and "UP"}
\item{icd10}{regular expressions identifying ICD-10 codes of each
group}
\item{only_ordinary}{index weights, 1 for ordinary and 0 for special}
}
}
\usage{
cps
}
\description{
Classcodes for the comorbidity-polypharmacy score (CPS) based on ICD-10 codes
}
\references{
Stawicki, Stanislaw P., et al.
"Comorbidity polypharmacy score and its clinical utility: A pragmatic
practitioner's perspective." Journal of emergencies, trauma, and shock 8.4
(2015): 224.
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4626940/}
}
\seealso{
Other default classcodes: 
\code{\link{ae}},
\code{\link{charlson}},
\code{\link{elixhauser}},
\code{\link{hip_ae_hailer}},
\code{\link{rxriskv}}
}
\concept{default classcodes}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify.R
\name{as.data.frame.classified}
\alias{as.data.frame.classified}
\alias{as.data.table.classified}
\alias{as.matrix.classified}
\title{Convert output from classify() to matrix/data.frame/data.table}
\usage{
\method{as.data.frame}{classified}(x, ...)

\method{as.data.table}{classified}(x, ...)

\method{as.matrix}{classified}(x, ...)
}
\arguments{
\item{x}{output from \code{\link[=classify]{classify()}}}

\item{...}{ignored}
}
\value{
data frame/data table with:
\itemize{
\item first column named as "id" column specified as input
to \code{\link[=classify]{classify()}} and with data from \code{row.names(x)}
\item all columns from \code{classified}
\item no row names
}

or simply the input matrix without additional attributes
}
\description{
Convert output from classify() to matrix/data.frame/data.table
}
\examples{
x <- classify(c("C80", "I20", "unvalid_code"), "elixhauser")

as.matrix(x)[, 1:3]
as.data.frame(x)[, 1:3]
data.table::as.data.table(x)[, 1:3]

# `as_tibble()` works automatically due to internal use of `as.data.frame()`.
tibble::as_tibble(x)
}
\seealso{
Other classcodes: 
\code{\link{all_classcodes}()},
\code{\link{classcodes}},
\code{\link{codebook}()},
\code{\link{print.classcodes}()},
\code{\link{print.classified}()},
\code{\link{set_classcodes}()},
\code{\link{summary.classcodes}()},
\code{\link{visualize.classcodes}()}
}
\concept{classcodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classcodes.R
\name{classcodes}
\alias{classcodes}
\alias{as.classcodes}
\alias{as.classcodes.classcodes}
\alias{as.classcodes.data.frame}
\alias{is.classcodes}
\title{Classcodes methods}
\usage{
as.classcodes(x, ...)

\method{as.classcodes}{classcodes}(
  x,
  ...,
  regex = attr(x, "regexpr"),
  indices = attr(x, "indices"),
  hierarchy = attr(x, "hierarchy")
)

\method{as.classcodes}{data.frame}(
  x,
  ...,
  regex = NULL,
  indices = NULL,
  hierarchy = attr(x, "hierarchy"),
  .name = NULL
)

is.classcodes(x)
}
\arguments{
\item{x}{data frame with columns described in the details section.
Alternatively a \code{classcodes} object to be modified.}

\item{...}{arguments passed between methods#'}

\item{regex, indices}{character vector with names of columns in \code{x} containing
regular expressions/indices.}

\item{hierarchy}{named list of pairwise group names to appear as superior and
subordinate for indices.
To be used for indexing when the subordinate class is redundant
(see the details section of \code{\link{elixhauser}} for an example).}

\item{.name}{used internally for name dispatch}
}
\value{
Object of class \code{classcodes} (inheriting from data frame)
with additional attributes:
\itemize{
\item \verb{code:} the coding used (for example "icd10", or "ATC").
\code{NULL} for unknown/arbitrary coding.
\item \verb{regexprs:} name of columns with regular expressions
(as specified by the \code{regex}argument)
\item \verb{indices:} name of columns with (optional) index weights
(as specified by the \code{indices}argument)
\item \verb{hierarchy:} list as specified by the \code{hierarchy} argument.
\item \verb{name:} name as specified by the \code{.name} argument.
}
}
\description{
\code{classcodes} are classification schemes based on regular expression stored in
data frames. These are essential to the package and constitute the third
part of the triad of case data, code data and a classification scheme.
}
\details{
A classcodes object is a data frame with mandatory columns:
\itemize{
\item \code{group}: unique and non missing class names
\item At least one column with regular expressions
(\link{regex} without Perl-like versions) defining class
membership. Those columns can have arbitrary names
(as specified by the \code{regex} argument).
Occurrences of non unique regular expressions will lead to the same class
having multiple names. This is accepted but will raise a warning.
Classes do not have to be disjunct.
}

The object can have additional optional columns:
\itemize{
\item \code{description}: description of each category
\item \code{condition}: a class might have conditions additional to what
is expressed by the regular expressions.
If so, these should be specified as quoted
expressions that can be evaluated within the data frame used by
\code{\link[=classify]{classify()}}
\item weights for each class used by
\code{\link[=index]{index()}}. Could be more than one and could have arbitrary names
(as specified by the \code{indices}argument).
}
}
\examples{
# The Elixhauser comorbidity classification is already a classcodes object
is.classcodes(coder::elixhauser)

# Strip its class attributes to use in examples
df <- as.data.frame(coder::elixhauser)

# Specify which columns store regular expressions and indices
# (assume no hierarchy)
elix <-
  as.classcodes(
    df,
    regex     = c("icd10", "icd10_short", "icd9cm", "icd9cm_ahrqweb", "icd9cm_enhanced"),
    indices   = c("sum_all", "sum_all_ahrq", "walraven",
                "sid29", "sid30", "ahrq_mort", "ahrq_readm"),
    hierarchy = NULL
  )
elix

# Specify hierarchy for patients with different types of cancer and diabetes
# See `?elixhauser` for details
as.classcodes(
  elix,
  hierarchy = list(
    cancer   = c("metastatic cancer", "solid tumor"),
    diabetes = c("diabetes complicated", "diabetes uncomplicated")
  )
)

# Several checks are performed to not allow any erroneous classcodes object
\dontrun{
  as.classcodes(iris)
  as.classcodes(iris, regex = "Species")
}
}
\seealso{
\code{vignette("classcodes")}
\code{vignette("Interpret_regular_expressions")}
The package have several default classcodes included, see \code{\link[=all_classcodes]{all_classcodes()}}.

Other classcodes: 
\code{\link{all_classcodes}()},
\code{\link{as.data.frame.classified}()},
\code{\link{codebook}()},
\code{\link{print.classcodes}()},
\code{\link{print.classified}()},
\code{\link{set_classcodes}()},
\code{\link{summary.classcodes}()},
\code{\link{visualize.classcodes}()}

Other classcodes: 
\code{\link{all_classcodes}()},
\code{\link{as.data.frame.classified}()},
\code{\link{codebook}()},
\code{\link{print.classcodes}()},
\code{\link{print.classified}()},
\code{\link{set_classcodes}()},
\code{\link{summary.classcodes}()},
\code{\link{visualize.classcodes}()}
}
\concept{classcodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_for_datasets.R
\docType{data}
\name{ex_people}
\alias{ex_people}
\title{Example data for random people}
\format{
Data frames with 100 rows and 2 variables:
\describe{
\item{name}{random person names}
\item{surgery}{random dates for a relevant event}
}
}
\usage{
ex_people
}
\description{
Example data for fictive people to use for testing and in examples.
}
\seealso{
Other example data: 
\code{\link{ex_atc}},
\code{\link{ex_icd10}}
}
\concept{example data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_for_datasets.R
\docType{data}
\name{charlson}
\alias{charlson}
\title{Classcodes for Charlson comorbidity based on ICD-codes}
\format{
A data frame with 17 rows and 8 variables:
\itemize{
\item \verb{group:} comorbidity groups
\item \verb{description:} Verbal description of codes as described by
Deyo et al. (1992).
\item \verb{icd10:} regular expressions identifying ICD-10 codes of each
group as decoded from Quan et al. 2005. Note that this classification was
not originally used with all weights! To simply use this classification
table with weights other than \code{quan_original} and \code{quan_updated}
might therefore lead to different results than originally intended for each
index.
\item \verb{icd9cm_deyo:}Codes from table 1 column "Deyo's ICD-9-CM"
in Quan et al. (2005).
Procedure code 38.48 for peripheral vascular disease ignored.
\item \verb{icd9cm_enhanced:} Codes from table 1 column "Enhanced ICD-9-CM"
in Quan et al. (2005).
\item \verb{icd10_rcs:} Codification by Armitage (2010).
Note that Peptic ulcer disease is not included.
All liver diseases (including mild) are included in
"moderate or severe liver disease".
All diabetes is included in "diabetes complication"
\item \verb{icd8_brusselaers:} Back translated version from ICD-10 to
ICD-8 by Brusselaers et al. (2017).
"Moderate and severe liver disease" contains all liver disease and
"diabetes complication" contains all diabetes.
\item \verb{icd9_brusselaers:} Back translated version from ICD-10 to
ICD-9 by Brusselaers et al. (2017).
"Moderate and severe liver disease" contains all liver disease and
"diabetes complication" contains all diabetes.
\item \verb{charlson:} original weights as suggested by Charlson et al.
(1987)*
\item \verb{deyo_ramano:} weights suggested by Deyo and Romano*
\item \verb{dhoore:} weights suggested by D'Hoore*
\item \verb{ghali:} weights suggested by Ghali*
\item \verb{quan_original:} weights suggested by Quan (2005)
\item \verb{quan_updated:} weights suggested by Quan (2011)
}
\itemize{
\item Weights decoded from Yurkovich et al. (2015).
}
}
\usage{
charlson
}
\description{
Classcodes for Charlson comorbidity based on ICD-codes
}
\references{
Armitage, J. N., & van der Meulen, J. H. (2010).
Identifying co-morbidity in surgical patients using administrative data
with the Royal College of Surgeons Charlson Score.
British Journal of Surgery, 97(5), 772–781.
\doi{10.1002/bjs.6930}

Brusselaers N, Lagergren J. (2017)
The Charlson Comorbidity Index in Registry-based Research.
Methods Inf Med 2017;56:401–6. \doi{10.3414/ME17-01-0051}.

Deyo, R. A., Cherkin, D. C., & Ciol, M. A. (1992).
Adapting a clinical comorbidity index for use with ICD-9-CM
administrative databases.
Journal of Clinical Epidemiology, 45(6), 613–619.
\doi{10.1016/0895-4356(92)90133-8}

Quan Hude et al. (2005). Coding algorithms for defining
comorbidities in ICD-9-CM and ICD-10 administrative data.
Medical care, 1130-1139.
\url{https://www.jstor.org/stable/3768193}

Yurkovich, M., Avina-Zubieta, J. A., Thomas, J., Gorenchtein, M., & Lacaille,
D. (2015). A systematic review identifies valid comorbidity indices derived
from administrative health data.
Journal of clinical epidemiology, 68(1), 3-14.
}
\seealso{
Other default classcodes: 
\code{\link{ae}},
\code{\link{cps}},
\code{\link{elixhauser}},
\code{\link{hip_ae_hailer}},
\code{\link{rxriskv}}
}
\concept{default classcodes}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codify.R
\name{codify}
\alias{codify}
\alias{codify.data.frame}
\alias{codify.data.table}
\alias{print.codified}
\title{Codify case data with external code data (within specified time frames)}
\usage{
codify(x, codedata, ..., id, code, date = NULL, code_date = NULL, days = NULL)

\method{codify}{data.frame}(x, ..., id, date = NULL, days = NULL)

\method{codify}{data.table}(
  x,
  codedata,
  ...,
  id,
  code,
  date = NULL,
  code_date = NULL,
  days = NULL,
  alnum = FALSE,
  .copy = NA
)

\method{print}{codified}(x, ..., n = 10)
}
\arguments{
\item{x}{data set with mandatory character id column
(identified by argument \code{id = "<col_name>"}),
and optional \code{\link{Date}}  of interest
(identified by argument \code{date = "<col_name>"}).
Alternatively, the output from \code{\link[=codify]{codify()}}}

\item{codedata}{additional data with columns
including case id (\code{character}), code and an optional date (\link{Date}) for
each code. An optional column \code{condition} might distinguish codes/dates
with certain characteristics (see example).}

\item{...}{arguments passed between methods}

\item{id, code, date, code_date}{column names with case id
(\code{character} from \code{x} and \code{codedata}), \code{code} (from \code{x}) and
optional date (\link{Date} from \code{x}) and
\code{code_date} (\link{Date} from \code{codedata}).}

\item{days}{numeric vector of length two with lower and upper bound for range
of relevant days relative to \code{date}. See "Relevant period".}

\item{alnum}{Should codes be cleaned from all non alphanumeric characters?}

\item{.copy}{Should the object be copied internally by \code{\link[data.table:copy]{data.table::copy()}}?
\code{NA} (by default) means that objects smaller than 1 GB are copied.
If the size is larger, the argument must be set explicitly. Set \code{TRUE}
to make copies regardless of object size. This is recommended if enough RAM
is available. If set to \code{FALSE}, calculations might be carried out
but the object will be changed by reference.
IMPORTANT! This might lead to undesired consequences and should only be used
if absolutely necessary!}

\item{n}{number of rows to preview as tibble.
The output is technically a \link[data.table:data.table]{data.table::data.table}, which might be an
unusual format to look at. Use \code{n = NULL} to print the object as is.}
}
\value{
Object of class \code{codified} (inheriting from \link[data.table:data.table]{data.table::data.table}).
Essentially \code{x} with additional columns:
\verb{code, code_date}: left joined from \code{codedata} or \code{NA}
if no match within period. \code{in_period}: Boolean indicator if the case
had at least one code within the specified period.

The output has one row for each combination of "id" from \code{x} and
"code" from \code{codedata}. Rows from \code{x} might be repeated
accordingly.
}
\description{
This is the first step of \code{codify() \%>\% classify() \%>\% index()}.
The function combines case data from one data set with related code data from
a second source, possibly limited to codes valid at certain time points
relative to case dates.
}
\section{Relevant period}{

Some examples for argument \code{days}:
\itemize{
\item \code{c(-365, -1)}: window of one year prior to the \code{date}
column of \code{x}. Useful for patient comorbidity.
\item \code{c(1, 30)}: window of 30 days after \code{date}.
Useful for adverse events after a surgical procedure.
\item \code{c(-Inf, Inf)}: no limitation on non-missing dates.
\item \code{NULL}: no time limitation at all.
}
}

\examples{
# Codify all patients from `ex_people` with their ICD-10 codes from `ex_icd10`
x <- codify(ex_people, ex_icd10, id = "name", code = "icd10")
x

# Only consider codes if recorded at hospital admissions within one year prior
# to surgery
codify(
  ex_people,
  ex_icd10,
  id        = "name",
  code      = "icd10",
  date      = "surgery",
  code_date = "admission",
  days      = c(-365, 0)   # admission during one year before surgery
)

# Only consider codes if recorded after surgery
codify(
  ex_people,
  ex_icd10,
  id        = "name",
  code      = "icd10",
  date      = "surgery",
  code_date = "admission",
  days      = c(1, Inf)     # admission any time after surgery
)


# Dirty code data ---------------------------------------------------------

# Assume that codes contain unwanted "dirty" characters
# Those could for example be a dot used by ICD-10 (i.e. X12.3 instead of X123)
dirt <- c(strsplit(c("!#\%&/()=?`,.-_"), split = ""), recursive = TRUE)
rdirt <- function(x) sample(x, nrow(ex_icd10), replace = TRUE)
sub <- function(i) substr(ex_icd10$icd10, i, i)
ex_icd10$icd10 <-
  paste0(
    rdirt(dirt), sub(1),
    rdirt(dirt), sub(2),
    rdirt(dirt), sub(3),
    rdirt(dirt), sub(4),
    rdirt(dirt), sub(5)
  )
head(ex_icd10)

# Use `alnum = TRUE` to ignore non alphanumeric characters
codify(ex_people, ex_icd10, id = "name", code = "icd10", alnum = TRUE)



# Big data ----------------------------------------------------------------

# If `data` or `codedata` are large compared to available
# Random Access Memory (RAM) it might not be possible to make internal copies
# of those objects. Setting `.copy = FALSE` might help to overcome such problems

# If no copies are made internally, however, the input objects (if data tables)
# would change in the global environment
x2 <- data.table::as.data.table(ex_icd10)
head(x2) # Look at the "icd10" column (with dirty data)

# Use `alnum = TRUE` combined with `.copy = FALSE`
codify(ex_people, x2, id = "name", code = "icd10", alnum = TRUE, .copy = FALSE)

# Even though no explicit assignment was specified
# (neither for the output of codify(), nor to explicitly alter `x2`,
# the `x2` object has changed (look at the "icd10" column!):
head(x2)

# Hence, the `.copy` argument should only be used if necessary
# and if so, with caution!


# print.codify() ----------------------------------------------------------

x # Preview first 10 rows as a tibble
print(x, n = 20) # Preview first 20 rows as a tibble
print(x, n = NULL) # Print as data.table (ignoring the 'classified' class)
}
\seealso{
Other verbs: 
\code{\link{categorize}()},
\code{\link{classify}()},
\code{\link{index_fun}}
}
\concept{verbs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{visualize.classcodes}
\alias{visualize.classcodes}
\title{Visualize classification scheme in web browser}
\usage{
\method{visualize}{classcodes}(x, group = NULL, show = TRUE, ...)
}
\arguments{
\item{x}{\link{classcodes} object or name of such object included in the package
(see \code{\link[=all_classcodes]{all_classcodes()}}).}

\item{group}{names (as character vector) of groups to visualize
(subset of \code{rownames(x)}). (All groups if \code{NULL}.)}

\item{show}{should a visualization be shown in the default web browser.
Set to \code{FALSE} to just retrieve a URL for later use.}

\item{...}{
  Arguments passed on to \code{\link[=set_classcodes]{set_classcodes}}
  \describe{
    \item{\code{regex}}{name of column with regular expressions to use for
classification.
\code{NULL} (default) uses \code{attr(obj, "regexpr")[1]}.}
  }}
}
\value{
URL to website with visualization (invisible)
}
\description{
Groups from a \code{classcodes} object are visualized by their regular expressions
in the default web browser.
The visualization does not give any details on group names, conditions or
weights but might be useful both for understanding of a classification scheme
in use, and during the creation and debugging of such.
}
\examples{

# The default behavior is to open a visualization in the default web browser
\dontrun{

 # How is depression classified according to Elixhauser?
 visualize("elixhauser", "depression")

 # Compare the two diabetes groups according to Charlson
 visualize("charlson",
   c("diabetes without complication", "diabetes complication"))

 # Is this different from the "Royal College of Surgeons classification?
 # Yes, there is only one group for diabetes
 visualize("charlson",
   c("diabetes without complication", "diabetes complication"),
   regex = "rcs"
 )

 # Show all groups from Charlson
 visualize("charlson")

 # It is also possible to visualize an arbitrary regular expression
 # from a character string
 visualize("I2([12]|52)")
}

 # The URL is always returned (invisable) but the visual display can
 # also be omitted
url <- visualize("hip_ae", show = FALSE)
url
}
\seealso{
Other classcodes: 
\code{\link{all_classcodes}()},
\code{\link{as.data.frame.classified}()},
\code{\link{classcodes}},
\code{\link{codebook}()},
\code{\link{print.classcodes}()},
\code{\link{print.classified}()},
\code{\link{set_classcodes}()},
\code{\link{summary.classcodes}()}
}
\concept{classcodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{clean_text}
\alias{clean_text}
\title{Make clean text with only lowercase alphanumeric characters and "_"}
\usage{
clean_text(x_name, x)
}
\arguments{
\item{x_name}{Name of object to use as prefix}

\item{x}{character vector}
}
\value{
character vector of the same length as \code{x}
}
\description{
Make clean text with only lowercase alphanumeric characters and "_"
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_classcodes.R
\name{set_classcodes}
\alias{set_classcodes}
\title{Set classcodes object}
\usage{
set_classcodes(
  cc,
  classified = NULL,
  regex = NULL,
  start = TRUE,
  stop = FALSE,
  tech_names = NULL
)
}
\arguments{
\item{cc}{\code{\link{classcodes}} object (or name of a default object from
\code{\link[=all_classcodes]{all_classcodes()}}).}

\item{classified}{object that classcodes could be inherited from}

\item{regex}{name of column with regular expressions to use for
classification.
\code{NULL} (default) uses \code{attr(obj, "regexpr")[1]}.}

\item{start, stop}{should codes start/end with the specified regular
expressions? If \code{TRUE}, column "regex" is prefixed/suffixed
by \verb{^/$}.}

\item{tech_names}{should technical column names be used? If \code{FALSE},
colnames are taken directly from group names of \code{cc}, if \code{TRUE},
these are changed to more technical names avoiding special characters and
are prefixed by the name of the classification scheme.
\code{NULL} (by default) preserves previous names if \code{cc} is inherited from
\code{classified} (fall backs to \code{FALSE} if not already set).}
}
\value{
\code{\link{classcodes}} object.
}
\description{
Prepare a \code{classcodes}object by specifying the regular expressions
to use for classification.
}
\examples{
# Prepare a classcodes object for the Charlson comorbidity classification
# based on the default regular expressions
set_classcodes(charlson)   # by object
set_classcodes("charlson") # by name

# Same as above but based on regular expressions for ICD-8 (see `?charlson`)
set_classcodes(charlson, regex = "icd8_brusselaers")

# Only recognize codes if no other characters are found after the relevant codes
# Hence if the code vector stops with the code
set_classcodes(charlson, stop = TRUE)

# Accept code vectors with strings which do not necessarily start with the code.
# This is useful if the code might appear in the middle of a longer character
# string or if a common prefix is used for all codes.
set_classcodes(charlson, start = FALSE)

# Use technical names to clearly describe the origin of each group.
# Note that the `cc` argument must be specified by a character string
# since this name is used as part of the column names
x <- set_classcodes("charlson", tech_names = TRUE)
x$group
}
\seealso{
Other classcodes: 
\code{\link{all_classcodes}()},
\code{\link{as.data.frame.classified}()},
\code{\link{classcodes}},
\code{\link{codebook}()},
\code{\link{print.classcodes}()},
\code{\link{print.classified}()},
\code{\link{summary.classcodes}()},
\code{\link{visualize.classcodes}()}
}
\concept{classcodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codebook.R
\name{codebook}
\alias{codebook}
\alias{print.codebook}
\alias{codebooks}
\title{codebook(s) for classcodes object}
\usage{
codebook(object, coding, ..., file = NULL)

\method{print}{codebook}(x, ...)

codebooks(..., file = NULL)
}
\arguments{
\item{object}{classcodes object}

\item{coding}{either a vector with codes from the original classification,
or a name (character vector of length one) of a keyvalue object
from package "decoder" (for example "icd10cm" or "atc")}

\item{...}{Additional arguments for each function:
\itemize{
\item \code{codebook()}: arguments passed to \code{\link[=summary.classcodes]{summary.classcodes()}}
\item \code{codebooks()}: multiple named outputs from \code{\link[=codebook]{codebook()}}
\item \code{print.codebook()}: arguments passed to \code{tibble:::print.tbl()}
}}

\item{file}{name/path to Excel file for data export}

\item{x}{output from \code{codebook()}}
}
\value{
Functions are primarily called for their side effects (exporting data to
Excel or printing to screen). In addition:
\itemize{
\item \code{codebook()}returns list of data frames describing relationship
between groups and individual codes
\item \code{codebooks()} returns a concatenated list with output from \code{codebook()}.
Only one 'README' object is kept however and renamed as such.
\item \code{print.codebook()}returns \code{x} (invisible)
}
}
\description{
\code{\link[=summary.classcodes]{summary.classcodes()}} and \code{\link[=visualize.classcodes]{visualize.classcodes()}} are used to
summarize/visualize classcodes in R. A codebook, on the other hand,
is an exported summary
saved in an Excel spreadsheet to use in collaboration with non R-users.
Several codebooks might be combined into a single Excel document with
several sheets (one for each codebook).
}
\examples{
# codebook() --------------------------------------------------------------
\dontrun{
# Export codebook (to temporary file) with all codes identified by the
# Elixhauser comorbidity classification based on ICD-10-CM
codebook(elixhauser, "icd10cm", file = tempfile("codebook", fileext = ".xlsx"))

# All codes from ICD-9-CM Disease part used by Elixhauser enhanced version
codebook(elixhauser, "icd9cmd",
  cc_args = list(regex = "icd9cm_enhanced",
  file = tempfile("codebook", fileext = ".xlsx"))
)

# The codebook returns a list with three objects.
# Access a dictionary table with translates of each code to text:
codebook(charlson, "icd10cm")$all_codes
}

# print.codebook() --------------------------------------------------------

# If argument `file` is unspecified, a preview of each sheet of the codebook is
# printed to the screen
(cb <- codebook(charlson, "icd10cm"))

# The preview can be modified by arguments to the print-method
print(cb, n = 20)


# codebooks() -------------------------------------------------------------

# Combine codebooks based on different versions of the regular expressions
# and export to a single (temporary) Excel file
c1 <- codebook(elixhauser, "icd10cm")
c2 <- codebook(elixhauser, "icd9cmd",
  cc_args = list(regex = "icd9cm_enhanced")
  )

codebooks(
  elix_icd10 = c1, elix_icd9cm = c2,
  file = tempfile("codebooks", fileext = ".xlsx")
)

}
\seealso{
Other classcodes: 
\code{\link{all_classcodes}()},
\code{\link{as.data.frame.classified}()},
\code{\link{classcodes}},
\code{\link{print.classcodes}()},
\code{\link{print.classified}()},
\code{\link{set_classcodes}()},
\code{\link{summary.classcodes}()},
\code{\link{visualize.classcodes}()}
}
\concept{classcodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.classcodes.R
\name{summary.classcodes}
\alias{summary.classcodes}
\alias{print.summary.classcodes}
\title{Summarizing a classcodes object}
\usage{
\method{summary}{classcodes}(object, coding, ..., cc_args = list())

\method{print}{summary.classcodes}(x, ...)
}
\arguments{
\item{object}{classcodes object}

\item{coding}{either a vector with codes from the original classification,
or a name (character vector of length one) of a keyvalue object
from package "decoder" (for example "icd10cm" or "atc")}

\item{...}{\itemize{
\item \code{summary.classcodes()}: ignored
\item \code{print.summary.classcodes()}: arguments passed to \code{tibble:::print.tbl()}
}}

\item{cc_args}{List of named arguments passed to \code{\link[=set_classcodes]{set_classcodes()}}}

\item{x}{output from \code{summary.classcodes()}}
}
\value{
Methods primarily called for their side effects (printing to the screen) but
with additional invisible objects returned:
\itemize{
\item \code{summary.classcodes()}: list with input arguments \code{object} and \code{coding}
unchanged, as well as a data frame (\code{summary}) with columns for groups
identified (\code{group}); the number of codes to be recognized for each group
(\code{n}) and individual codes within each group (\code{codes}).
\item \code{print.summary.classcodes()}: argument \code{x} unchanged
}
}
\description{
Classification schemes are formalized by regular expressions within the
classcodes objects. These are computationally effective but sometimes hard to
interpret. Use this function to list all codes identified for each
group.
}
\examples{

# summary.classcodes() ----------------------------------------------------

# Summarize all ICD-10-CM codes identified by the Elixhauser
# comorbidity classification
# See `?decoder::icd10cm` for details
summary(elixhauser, coding = "icd10cm")

# Is there a difference if instead considering the Swedish ICD-10-SE?
# See `?decoder::icd10se` for details
summary(elixhauser, coding = "icd10se")

# Which ICD-9-CM diagnostics codes are recognized by Charlson according to
# Brusselears et al. 2017 (see `?charlson`)
summary(
  charlson, coding = "icd9cmd",
  cc_args = list(regex = "icd9_brusselaers")
)


# print.summary.classcodes() ----------------------------------------------

# Print all 31 lines of the summarized Elixhauser classcodes object
print(
  summary(elixhauser, coding = "icd10cm"),
  n = 31
)

}
\seealso{
Other classcodes: 
\code{\link{all_classcodes}()},
\code{\link{as.data.frame.classified}()},
\code{\link{classcodes}},
\code{\link{codebook}()},
\code{\link{print.classcodes}()},
\code{\link{print.classified}()},
\code{\link{set_classcodes}()},
\code{\link{visualize.classcodes}()}
}
\concept{classcodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_for_datasets.R
\docType{data}
\name{ae}
\alias{ae}
\alias{knee_ae}
\alias{hip_ae}
\title{Classcodes for adverse events after knee and hip arthroplasty}
\format{
Data frame with 3 columns:
\describe{
\item{group}{Different types of adverse events (see reference section)}
\item{icd10}{regular expressions identifying ICD-10 codes for each
group}
\item{icd10_fracture}{regular expressions for fracture patients.
Essentially the same as \code{regex} but with some additional codes for
group "DM1 other"}
\item{kva}{regular expressions identifying KVA codes}
\item{condition}{special conditions are used, see below.}
}

An object of class \code{classcodes} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 7 rows and 5 columns.
}
\source{
Knee (p. 83): \url{http://www.myknee.se/pdf/SVK-2016_1.1.pdf}.

Hip (p. 162): \url{https://registercentrum.blob.core.windows.net/shpr/r/Arsrapport_2018_Hoftprotes_final_web-rJgg8LvkOB.pdf}
}
\usage{
knee_ae

hip_ae
}
\description{
ICD-10 group names are prefixed by two letters as given by the references.
Two groups (DB and DM) are split into two due to different conditions.
}
\section{Hip fractures}{

Adverse events (AE) codes for hip fractures are based on codes for elective
cases but with some additional codes for DM 1 (N300, N308, N309 and N390).
}

\section{Conditions}{

Special conditions apply to all categories.
Those require non-standard modifications
of the classcodes data prior to categorization.

\describe{
\item{hbdia1_hdia}{\code{TRUE} if the code was
given as any type of diagnose during hospital visit for index operation,
or as main diagnose for later visits, otherwise \code{FALSE}}
\item{late_hdia}{\code{TRUE} if the code was
given as main diagnose at a later visit after the index operation,
otherwise \code{FALSE}}
\item{post_op}{\code{TRUE} if the code was
given at a later visit after the index operation, otherwise \code{FALSE}}
}
}

\references{
Magneli M, Unbeck M, Rogmark C, Rolfson O, Hommel A, Samuelsson B, et al.
Validation of adverse events after hip arthroplasty:
a Swedish multi-centre cohort study.
BMJ Open. 2019 Mar 7;9(3):e023773.
Available from: \url{https://pubmed.ncbi.nlm.nih.gov/30850403/}
}
\seealso{
hip_ae_hailer

Other default classcodes: 
\code{\link{charlson}},
\code{\link{cps}},
\code{\link{elixhauser}},
\code{\link{hip_ae_hailer}},
\code{\link{rxriskv}}
}
\concept{default classcodes}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index_fun.R
\name{index_fun}
\alias{index_fun}
\alias{index}
\alias{index.data.frame}
\alias{index.matrix}
\title{Calculate index based on classification scheme}
\usage{
index(classified, ...)

\method{index}{data.frame}(classified, ...)

\method{index}{matrix}(classified, index = NULL, cc = NULL, ...)
}
\arguments{
\item{classified}{output from \code{\link[=classify]{classify()}}}

\item{...}{used internally}

\item{index}{name of column with 'weights' from corresponding
\code{\link{classcodes}} object. Can be \code{NULL} if the index is just a unweighted
count of all identified groups.}

\item{cc}{\code{\link{classcodes}} object. Can be \code{NULL} if a \code{classcodes} object is
already available as an attribute of \code{classified} (which is often the case)
and/or if \code{index = NULL}.}
}
\value{
Named numeric index vector with names corresponding to
\code{rownames(classified)}
}
\description{
This is the third step of \code{codify() \%>\% classify() \%>\% index()}.
The function takes classified case data and calculates
(weighted) index sums as specified by weights from a \code{classcodes} object.
}
\details{
Index weights for subordinate hierarchical classes
(as identified by \code{attr(cc, "hierarchy")}) are excluded in presence of
superior classes if index specified with argument \code{index}.
}
\examples{

# Prepare some codified data with ICD-10 codes during 1 year (365 days)
# before surgery
x <-
  codify(
    ex_people,
    ex_icd10,
    id        = "name",
    code      = "icd10",
    date      = "surgery",
    days      = c(-365, 0),
    code_date = "admission"
  )

# Classify those patients by the Charlson comorbidity indices
cl <- classify(x, "charlson")

# Calculate (weighted) index values
head(index(cl))                  # Un-weighted sum/no of conditions for each patient
head(index(cl, "quan_original")) # Weighted index (Quan et al. 2005; see `?charlson`)
head(index(cl, "quan_updated"))  # Weighted index (Quan et al. 2011; see `?charlson`)

# Tabulate index for all patients.
# As expected, most patients are healthy and have index = 0/NA,
# where NA indicates no recorded hospital visits
# found in `ex_icd10` during codification.
# In practice, those patients might be assumed to have 0 comorbidity as well.
table(index(cl, "quan_original"), useNA = "always")

# If `cl` is a matrix without additional attributes (as imposed by `codify()`)
# an explicit classcodes object must be specified by the `cc` argument
cl2 <- as.matrix(cl)
head(index(cl2, cc = "charlson"))
}
\seealso{
Other verbs: 
\code{\link{categorize}()},
\code{\link{classify}()},
\code{\link{codify}()}
}
\concept{verbs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.keyvalue.R
\name{as.keyvalue.classcodes}
\alias{as.keyvalue.classcodes}
\title{Make keyvalue object from classcodes object}
\usage{
\method{as.keyvalue}{classcodes}(x, coding, cc_args = list())
}
\arguments{
\item{x}{classcodes object}

\item{coding}{either a vector with codes from the original classification,
or a name (character vector of length one) of a keyvalue object
from package "decoder" (for example "icd10cm" or "atc")}

\item{cc_args}{List of named arguments passed to \code{\link[=set_classcodes]{set_classcodes()}}}
}
\value{
Object of class \code{keyvalue} where \code{key} is the subset of codes from
\code{object$key}identified by the regular expression from \code{x} and where
\code{value} is the corresponding \code{x$group}. Hence, note that the original
\code{object$value} is not used in the output.
}
\description{
S3-method for generic \code{\link[decoder:keyvalue]{decoder::as.keyvalue()}}
}
\examples{
# List all codes with corresponding classes as recognized by the Elixhauser
# comorbidity classification according to the Swedish version of the
# international classification of diseases version 10 (ICD-10-SE)
head(decoder::as.keyvalue(elixhauser, "icd10se"))

# Similar but with the American ICD-10-CM instead
# Note that the `value` column is similar as above
# (with names from `x$group`) and not
# from `object$value`
head(decoder::as.keyvalue(elixhauser, "icd10cm"))

# Codes identified by regular expressions based on ICD-9-CM and found in
# the Swedish version of ICD-9 used within the national cancer register
# (thus, a subset of the whole classification).
head(
  decoder::as.keyvalue(
    elixhauser, "icd9",
    cc_args = list(regex = "icd9cm")
  )
)
}
\concept{helper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexports.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\alias{visualize}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{generics}{\code{\link[generics]{visualize}}}

  \item{tibble}{\code{\link[tibble:reexports]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_for_datasets.R
\docType{data}
\name{hip_ae_hailer}
\alias{hip_ae_hailer}
\title{Classcodes for infection and dislocation after hip arthroplasty}
\format{
Data frame with 3 columns:
\describe{
\item{group}{Infection or dislocation}
\item{icd10}{regular expressions based on ICD-10}
\item{kva}{regular expressions based on NOMESCO/KVA codes}
}
}
\usage{
hip_ae_hailer
}
\description{
Classcodes for infection and dislocation after hip arthroplasty
}
\seealso{
ae

Other default classcodes: 
\code{\link{ae}},
\code{\link{charlson}},
\code{\link{cps}},
\code{\link{elixhauser}},
\code{\link{rxriskv}}
}
\concept{default classcodes}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{cols}
\alias{cols}
\title{Return all columns from x with names matching "find"}
\usage{
cols(find, x)
}
\arguments{
\item{find}{character vector with names to match}

\item{x}{matrix}
}
\description{
Return all columns from x with names matching "find"
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify.R
\name{print.classified}
\alias{print.classified}
\title{Printing classified data}
\usage{
\method{print}{classified}(x, ...)
}
\arguments{
\item{x}{output from \code{\link[=classify]{classify()}}}

\item{...}{additional arguments passed to printing method for a \code{tibble}.
\code{n} is the number of rows to preview.
Set \code{n = NULL} to disable the \code{tibble}
preview and print the object as is (a matrix).}
}
\description{
Preview first \code{n} rows as tibble
}
\examples{
# Preview all output
classify(c("C80", "I20", "unvalid_code"), "elixhauser")

# Preview only the first row
print(classify(c("C80", "I20", "unvalid_code"), "elixhauser"), n = 1)

# Print object as is (matrix)
print(classify(c("C80", "I20", "unvalid_code"), "elixhauser"), n = NULL)
}
\seealso{
Other classcodes: 
\code{\link{all_classcodes}()},
\code{\link{as.data.frame.classified}()},
\code{\link{classcodes}},
\code{\link{codebook}()},
\code{\link{print.classcodes}()},
\code{\link{set_classcodes}()},
\code{\link{summary.classcodes}()},
\code{\link{visualize.classcodes}()}
}
\concept{classcodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify.R
\name{classify}
\alias{classify}
\alias{classify.default}
\alias{classify.codified}
\alias{classify.data.frame}
\alias{classify.data.table}
\title{Classify codified data}
\usage{
classify(codified, cc, ..., cc_args = list())

\method{classify}{default}(codified, cc, ..., cc_args = list())

\method{classify}{codified}(codified, ...)

\method{classify}{data.frame}(codified, ...)

\method{classify}{data.table}(codified, cc, ..., id, code, cc_args = list())
}
\arguments{
\item{codified}{output from \code{\link[=codify]{codify()}}}

\item{cc}{\code{\link{classcodes}} object (or name of a default object from
\code{\link[=all_classcodes]{all_classcodes()}}).}

\item{...}{arguments passed between methods}

\item{cc_args}{List with named arguments passed to
\code{\link[=set_classcodes]{set_classcodes()}}}

\item{code, id}{name of code/id columns (in \code{codified}).}
}
\value{
Object of class "classified". Inheriting from a Boolean matrix with
one row for each element/row of \code{codified}
and columns for each class with corresponding class names (according to the
\code{\link{classcodes}} object). Note, however, that \code{\link[=print.classified]{print.classified()}} preview
this output as a tibble.
}
\description{
This is the second step of \code{codify() \%>\% classify() \%>\% index()}.
Hence, the function takes a codified data set and classify each case based on
relevant codes as identified by the classification scheme provided by a
\code{classcodes} object.
}
\examples{


# classify.default() ------------------------------------------------------

# Classify individual ICD10-codes by Elixhauser
classify(c("C80", "I20", "unvalid_code"), "elixhauser")



# classify.codified() -----------------------------------------------------

# Prepare some codified data with ICD-10 codes during 1 year (365 days)
# before surgery
x <-
  codify(
    ex_people,
    ex_icd10,
    id        = "name",
    code      = "icd10",
    date      = "surgery",
    days      = c(-365, 0),
    code_date = "admission"
  )

# Classify those patients by the Charlson and Elixhasuer comorbidity indices
classify(x, "charlson")        # classcodes object by name ...
classify(x, coder::elixhauser) # ... or by the object itself


# -- start/stop --
# Assume that a prefix "ICD-10 = " is used for all codes and that some
# additional numbers are added to the end
x$icd10 <- paste0("ICD-10 = ", x$icd10)

# Set start = FALSE to identify codes which are not necessarily found in the
# beginning of the string
classify(x, "charlson", cc_args = list(start = FALSE))


# -- regex --
# Use a different version of Charlson (as formulated by regular expressions
# according to the Royal College of Surgeons (RCS) by passing arguments to
# `set_classcodes()` using the `cc_args` argument
y <-
  classify(
    x,
    "charlson",
    cc_args = list(regex = "icd10_rcs")
  )


# -- tech_names --
# Assume that we want to compare the results using the default ICD-10
# formulations (from Quan et al. 2005) and the RCS version and that the result
# should be put into the same data frame. We can use `tech_names = TRUE`
# to distinguish variables with otherwise similar names
cc <- list(tech_names = TRUE) # Prepare sommon settings
compare <-
  merge(
  classify(x, "charlson", cc_args = cc),
  classify(x, "charlson", cc_args = c(cc, regex = "icd10_rcs"))
)
names(compare) # long but informative and distinguishable column names



# classify.data.frame() / classify.data.table() ------------------------

# Assume that `x` is a data.frame/data.table without additional attributes
# from `codify()` ...
xdf <- as.data.frame(x)
xdt <- data.table::as.data.table(x)

# ... then the `id` and `code` columns must be specified explicitly
classify(xdf, "charlson", id = "name", code = "icd10")
classify(xdt, "charlson", id = "name", code = "icd10")
}
\seealso{
\code{\link[=as.data.frame.classified]{as.data.frame.classified()}}, \code{\link[=as.data.table.classified]{as.data.table.classified()}} and
\code{\link[=as.matrix.classified]{as.matrix.classified()}}, \code{\link[=print.classified]{print.classified()}}

Other verbs: 
\code{\link{categorize}()},
\code{\link{codify}()},
\code{\link{index_fun}}
}
\concept{verbs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/categorize.R
\name{categorize}
\alias{categorize}
\alias{categorize.data.frame}
\alias{categorize.tbl_df}
\alias{categorize.data.table}
\alias{categorize.codified}
\title{Categorize cases based on external data and classification scheme}
\usage{
categorize(x, ...)

\method{categorize}{data.frame}(x, ...)

\method{categorize}{tbl_df}(x, ...)

\method{categorize}{data.table}(x, ..., codedata, id, code, codify_args = list())

\method{categorize}{codified}(
  x,
  ...,
  cc,
  index = NULL,
  cc_args = list(),
  check.names = TRUE,
  .data_cols = NULL
)
}
\arguments{
\item{x}{data set with mandatory character id column
(identified by argument \code{id = "<col_name>"}),
and optional \code{\link{Date}}  of interest
(identified by argument \code{date = "<col_name>"}).
Alternatively, the output from \code{\link[=codify]{codify()}}}

\item{...}{arguments passed between methods}

\item{codedata}{external code data with mandatory character id column
(identified by \code{id = "<col_name>"}),
code column (identified by argument \code{code = "<col_name>"})
and optional \code{\link{Date}} column
(identified by \code{codify_args = list(code_date = "<col_name>")}).}

\item{id}{name of unique character id column found in
both \code{x}and \code{codedata}.
(where it must not be unique).}

\item{code}{name of code column in \code{codedata}.}

\item{codify_args}{Lists of named arguments passed to \code{\link[=codify]{codify()}}}

\item{cc}{\code{\link{classcodes}} object (or name of a default object from
\code{\link[=all_classcodes]{all_classcodes()}}).}

\item{index}{Argument passed to \code{\link[=index]{index()}}.
A character vector of names of columns with index weights from the
corresponding classcodes object (as supplied by the \code{cc}argument).
See \code{attr(cc, "indices")} for available options.
Set to \code{FALSE} if no index should be calculated.
If \code{NULL}, the default, all available indices (from \code{attr(cc, "indices")})
are provided.}

\item{cc_args}{List with named arguments passed to
\code{\link[=set_classcodes]{set_classcodes()}}}

\item{check.names}{Column names are based on \code{cc$group}, which might include
spaces. Those names are changed to syntactically correct names by
\code{check.names = TRUE}. Syntactically invalid, but grammatically correct
names might be preferred for presentation of the data as achieved by
\code{check.names = FALSE}. Alternatively, if \code{categorize} is called repeatedly,
longer informative names might be created by
\code{cc_args = list(tech_names = TRUE)}.}

\item{.data_cols}{used internally}
}
\value{
Object of the same class as \code{x} with additional logical columns
indicating membership of groups identified by the
\code{classcodes} object (the \code{cc} argument).
Numeric indices are also included if requested by the \code{index} argument.
}
\description{
This is the main function of the package, which relies of a triad of objects:
(1) \code{data} with unit id:s and possible dates of interest;
(2) \code{codedata} for corresponding
units and with optional dates of interest and;
(3) a classification scheme (\code{\link{classcodes}} object; \code{cc}) with regular
expressions to identify and categorize relevant codes.
The function combines the three underlying steps performed by
\code{\link[=codify]{codify()}}, \code{\link[=classify]{classify()}} and \code{\link[=index]{index()}}.
Relevant arguments are passed to those functions by
\code{codify_args} and \code{cc_args}.
}
\examples{
# For some patient data (ex_people) and related hospital visit code data
# with ICD 10-codes (ex_icd10), add the Elixhauser comorbidity
# conditions based on all registered ICD10-codes
categorize(
   x            = ex_people,
   codedata     = ex_icd10,
   cc           = "elixhauser",
   id           = "name",
   code         = "icd10"
)


# Add Charlson categories and two versions of a calculated index
# ("quan_original" and "quan_updated").
categorize(
   x            = ex_people,
   codedata     = ex_icd10,
   cc           = "charlson",
   id           = "name",
   code         = "icd10",
   index        = c("quan_original", "quan_updated")
)


# Only include recent hospital visits within 30 days before surgery,
categorize(
   x            = ex_people,
   codedata     = ex_icd10,
   cc           = "charlson",
   id           = "name",
   code         = "icd10",
   index        = c("quan_original", "quan_updated"),
   codify_args  = list(
      date      = "surgery",
      days      = c(-30, -1),
      code_date = "admission"
   )
)



# Multiple versions -------------------------------------------------------

# We can compare categorization by according to Quan et al. (2005); "icd10",
# and Armitage et al. (2010); "icd10_rcs" (see `?charlson`)
# Note the use of `tech_names = TRUE` to distinguish the column names from the
# two versions.

# We first specify some common settings ...
ind <- c("quan_original", "quan_updated")
cd  <- list(date = "surgery", days = c(-30, -1), code_date = "admission")

# ... we then categorize once with "icd10" as the default regular expression ...
categorize(
   x            = ex_people,
   codedata     = ex_icd10,
   cc           = "charlson",
   id           = "name",
   code         = "icd10",
   index        = ind,
   codify_args  = cd,
   cc_args      = list(tech_names = TRUE)
) \%>\%

# .. and once more with `regex = "icd10_rcs"`
categorize(
   codedata     = ex_icd10,
   cc           = "charlson",
   id           = "name",
   code         = "icd10",
   index        = ind,
   codify_args  = cd,
   cc_args      = list(regex = "icd10_rcs", tech_names = TRUE)
)



# column names ------------------------------------------------------------

# Default column names are based on row names from corresponding classcodes
# object but are modified to be syntactically correct.
default <-
   categorize(ex_people, codedata = ex_icd10, cc = "elixhauser",
              id = "name", code = "icd10")

# Set `check.names = FALSE` to retain original names:
original <-
  categorize(
    ex_people, codedata = ex_icd10, cc = "elixhauser",
    id = "name", code = "icd10",
    check.names = FALSE
   )

# Or use `tech_names = TRUE` for informative but long names (use case above)
tech <-
  categorize(ex_people, codedata = ex_icd10, cc = "elixhauser",
    id = "name", code = "icd10",
    cc_args = list(tech_names = TRUE)
  )

# Compare
tibble::tibble(names(default), names(original), names(tech))
}
\seealso{
Other verbs: 
\code{\link{classify}()},
\code{\link{codify}()},
\code{\link{index_fun}}
}
\concept{verbs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coder-package.R
\docType{package}
\name{coder-package}
\alias{coder}
\alias{coder-package}
\title{coder: Deterministic Categorization of Items Based on External Code Data}
\description{
Fast categorization of items based on external code data identified by 
  regular expressions. A typical use case considers patient with medically coded 
  data, such as codes from the International Classification of Diseases ('ICD') or 
  the Anatomic Therapeutic Chemical ('ATC') classification system. 
  Functions of the package relies on a triad of objects: (1) case data with unit 
  id:s and possible dates of interest; (2) external code data for corresponding 
  units in (1) and with optional dates of interest and; (3) a classification 
  scheme ('classcodes' object) with regular expressions to identify and 
  categorize relevant codes from (2). 
  It is easy to introduce new classification schemes ('classcodes' objects) or  
  to use default schemes included in the package. Use cases includes patient 
  categorization based on 'comorbidity indices' such as 'Charlson', 'Elixhauser', 
  'RxRisk V', or the 'comorbidity-polypharmacy' score (CPS), as well as adverse 
  events after hip and knee replacement surgery.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/coder/}
  \item Report bugs at \url{https://github.com/ropensci/coder/issues}
}

}
\author{
\strong{Maintainer}: Erik Bulow \email{eriklgb@gmail.com} (\href{https://orcid.org/0000-0002-9973-456X}{ORCID})

Other contributors:
\itemize{
  \item Emely C Zabore (Emily reviewed the package (v. 0.12.1) for rOpenSci, see <https://github.com/ropensci/software-review/issues/381>) [reviewer]
  \item David Robinson (David reviewed the package (v. 0.12.1) for rOpenSci, see <https://github.com/ropensci/software-review/issues/381>) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_for_datasets.R
\docType{data}
\name{rxriskv}
\alias{rxriskv}
\title{Classcodes for RxRisk V based on ATC codes}
\format{
Data frames with 46 rows and 6 variables:
\describe{
\item{group}{medical condition}
\item{pratt}{ATC codes from table 1 in Pratt et al. 2018
(ignoring PBS item codes and extra conditions).}
\item{garland}{Modified version by Anne
Garland to resemble medical use in Sweden 2016 (Unpublished).}
\item{caughey}{From appendix 1 in Caughey et al. 2010}
\item{pratt}{Mortality weights from table 1 in Pratt et al. 2018}
\item{sum_all}{Unweighted count of all conditions.}
}
}
\usage{
rxriskv
}
\description{
Note that desired implementation might differ over time and by country.
}
\references{
Caughey GE, Roughead EE, Vitry AI, McDermott RA, Shakib S, Gilbert AL.
Comorbidity in the elderly with diabetes:
Identification of areas of potential treatment conflicts.
Diabetes Res Clin Pract 2010;87:385–93.
\doi{10.1016/j.diabres.2009.10.019}.

Pratt NL, Kerr M, Barratt JD, Kemp-Casey A, Kalisch Ellett LM,
Ramsay E, et al.
The validity of the Rx-Risk Comorbidity Index using medicines mapped to
the Anatomical Therapeutic Chemical (ATC) Classification System.
BMJ Open 2018;8.
}
\seealso{
Other default classcodes: 
\code{\link{ae}},
\code{\link{charlson}},
\code{\link{cps}},
\code{\link{elixhauser}},
\code{\link{hip_ae_hailer}}
}
\concept{default classcodes}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_classcodes.R
\name{all_classcodes}
\alias{all_classcodes}
\title{Summary data for all default classcodes object in the package}
\usage{
all_classcodes()
}
\value{
\link[tibble:tibble]{tibble::tibble} with columns describing all default classcodes
objects from the package.
}
\description{
Tabulate object names and list all related versions of implemented regular
expressions
and index weights.
}
\examples{
all_classcodes()
}
\seealso{
Other classcodes: 
\code{\link{as.data.frame.classified}()},
\code{\link{classcodes}},
\code{\link{codebook}()},
\code{\link{print.classcodes}()},
\code{\link{print.classified}()},
\code{\link{set_classcodes}()},
\code{\link{summary.classcodes}()},
\code{\link{visualize.classcodes}()}
}
\concept{classcodes}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_for_datasets.R
\docType{data}
\name{elixhauser}
\alias{elixhauser}
\title{Classcodes for Elixhauser based on ICD-codes}
\format{
A data frame with 31 rows and 8 variables:
\describe{
\item{group}{comorbidity groups}
\item{icd10}{regular expressions identifying ICD-10 codes of
each group. Corresponds to column 'ICD-10' in table 2
of Quan et al. (2005).}
\item{icd10_short}{regular expressions identifying only the first
three characters of ICD-10 codes of each group. This alternative version
was added only to use in emergency when only the first three digits are
available. It is not an official version and we do not recommend to
use it!}
\item{icd9cm}{Corresponds to column 'Elixhauser's Original ICD-9-CM'
in table 2 of Quan et al. (2005).}
\item{icd9cm_ahrqweb}{Corresponds to column
'Elixhauser AHRQ-Web ICD-9-CM' in table 2 of Quan et al. (2005).}
\item{icd9cm_enhanced}{Corresponds to column 'Enhanced ICD-9-CM'
in table 2 of Quan et al. (2005).}
\item{sum_all}{all weights = 1 (thus no weights)}
\item{sum_all_ahrq}{as \code{sum_all} excluding "cardiac arrhythmia.
Compare to \code{icd9cm_ahrqweb} which does not
consider this condition.}
\item{walraven}{weights suggested by Walraven et al. (2009)}
\item{sid29}{weights suggested by Thompson et al. (2015)
based on all conditions except cardiac arrhythmia}
\item{sid30}{weights suggested by Thompson et al. (2015)
based on all conditions}
\item{ahrq_mort}{weights for in-hospital mortality suggested by
Moore et al. (2017)}
\item{ahrq_readm}{weights for readmissions suggested by
Moore et al. (2017)}
}
}
\usage{
elixhauser
}
\description{
Solid tumors are subordinate to metastatic cancer. A patient with both
conditions will still be classified as such but a possible (weighted)
index value will only account for metastatic cancer. The same is true for
"diabetes uncomplicated", which is subordinate of "diabetes complicated".
See Elixhauser et al. (1998).
}
\details{
Note that "DRG screen" as proposed in table 1 of Elixhauser et al. (1998)
is not handled by the coder package. This should instead be considered as
an additional pre- or post-processing step!
}
\references{
Elixhauser A, Steiner C, Harris DR, Coffey RM (1998).
Comorbidity Measures for Use with Administrative Data.
Med Care. 1998;36(1):8–27.

Moore, B. J., White, S., Washington, R., Coenen, N., & Elixhauser, A. (2017).
Identifying Increased Risk of Readmission and In-hospital Mortality Using
Hospital Administrative Data.
Medical Care, 55(7), 698–705.
\doi{10.1097/MLR.0000000000000735}

Quan Hude et al. (2005). Coding algorithms for defining
comorbidities in ICD-9-CM and ICD-10 administrative data.
Medical care, 1130-1139.
\url{https://www.jstor.org/stable/3768193}

Thompson, N. R., Fan, Y., Dalton, J. E., Jehi, L., Rosenbaum, B. P.,
Vadera, S., & Griffith, S. D. (2015).
A new Elixhauser-based comorbidity summary measure to predict in-hospital
mortality. Med Care, 53(4), 374–379.
\doi{10.1097/MLR.0000000000000326}

Walraven, C. Van, Austin, P. C., Jennings, A., Quan, H., Alan, J., Walraven,
C. Van, … Jennings, A. (2009).
A Modification of the Elixhauser Comorbidity Measures Into a Point System
for Hospital Death Using Administrative Data.
Medical Care, 47(6), 626–633.
}
\seealso{
Other default classcodes: 
\code{\link{ae}},
\code{\link{charlson}},
\code{\link{cps}},
\code{\link{hip_ae_hailer}},
\code{\link{rxriskv}}
}
\concept{default classcodes}
\keyword{datasets}
