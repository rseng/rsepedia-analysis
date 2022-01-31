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
---
title: Data Quality Report
nep: `r ccd@nepisodes`
np: `r max(ccd@infotb$pid, na.rm=TRUE)`
site: `r s <- unique(ccd@infotb$site_id); s[s!="NA"]`
admt: `r min(ccd@infotb$t_admission, na.rm=TRUE)`
disct: `r max(ccd@infotb$t_discharge, na.rm=TRUE)`
dbupdate: `r max(ccd@infotb$parse_time)`
fulldb: `r dbfull`

---

    
```{r, echo=FALSE, warning=FALSE}
item2d <- c("h_rate", "bp_m_a", "bp_sys_a", "lactate_abg", "bilirubin",
            "fluid_balance_h", "pao2_abg")

new.config <- function(items) {
    conf <- list()
    for (i in items) conf[[stname2code(i)]] <- list()
    conf
}

cctb <- create.cctable(ccd, conf=new.config(item2d), freq=1)$torigin

suppressWarnings(library(pander))
panderOptions("digits", 2) 
panderOptions("table.split.table", Inf)
demg <- sql.demographic.table(ccd)
```



# Data Summary
```{r, echo=FALSE}
dp <- total.data.point(ccd)
np <- max(ccd@infotb$pid, na.rm=TRUE)
```

This database contains **`r ccd@nepisodes`** episode data from 
**`r length(s[s!="NA"])`** sites. Based on the NHS numbers and
PAS numbers, we can identify **`r np`** unique patients, among which the earliest
admission is `r min(ccd@infotb$t_admission, na.rm=TRUE)` and the latest discharge time is
`r max(ccd@infotb$t_discharge, na.rm=TRUE)`. There are `r dp `
total data points found in the current database, which makes in average
`r round(dp/np)` per unique patient. 


## Site reference
```{r, echo=FALSE}
pander(as.data.frame(site.info()[, 1:3], style="rmarkdown"))
```

## The original XML files and parse information
```{r, echo=FALSE, results="asis"}
fs <- file.summary(ccd)
pander(as.data.frame(fs[, c("File", "Number of Episode", "Sites"), with=FALSE]))
```

```{r fig.width=15, fig.height=10, echo=FALSE}
xml.file.duration.plot(ccd)
```
```{r fig.width=15, fig.height=10, echo=FALSE}
xml.site.duration.plot(ccd)
```
\newpage

```{r, echo=FALSE, results="asis"}
table1(demg, c("SEX", "ETHNIC", "LOCA", "DIS", "HDIS"))
```
# Completeness
## Demographic data completeness
```{r, echo=FALSE, warning=FALSE, results="asis"}
demographic.data.completeness(demg)
```
## Sample period of time-wise data
This is the section to show the completeness of some key time-wise data, e.g.
physiological data. The completeness of data is measured by sample period. 
The sample period $P$ is defined as the ratio of the total
admission hour $T$ to the number of valid hours $D$ where data can be found.
$$P = \frac{T}{D}$$




```{r, echo=FALSE, warning=FALSE, results="asis"}
samplerate2d(cctb)
```



\newpage
# Data Distribution
## Demographic
```{r, echo=FALSE, results="asis", warning=FALSE}
demg.distribution(demg, c("HCM", "apache_score"))
```


## Physiological & Drugs
```{r, echo=FALSE, results="asis", warning=FALSE}
physio.distribution(cctb, item2d)
```

---
title: "Introduction to CCHIC critical care data"
author: David Perez Suarez & Sinan Shi
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Introduction to CCHIC critical care date}
  \usepackage[utf8]{inputenc}
---

CCHIC data strongly impacts the design of cleanEHR. Complex, heterogeneous, and
high resolution longitudinal data is a trend of EHR analysis, which takes
advantage of the ever more sophisticated statistical techniques and growing
computation capability. CCHIC is a representative example database of such
kind. It records 263 fields including 154 time-varying fields of patients
during their stay in intensive care units across five NHS trusts in England.
The recorded variables include patient demographics, time of admission and
discharge, survival status, diagnosis, physiology, laboratory, nursing and
drug. The latest database contains 22,628 admissions from 2014 to 2016, with
about 119 million data points (~6k per patient).

The anonymised data subset can be obtained [here](https://form.jotformeu.com/drstevok/cchic-end-user-license---cleanEHR). 
The selected number of researchers can get the access to the identifiable data 
[UCL IDHS](http://www.ucl.ac.uk/isd/itforslms/services/handling-sens-data/tech-soln).

### Episode and ICU stay timeline
We introduced the concept of 'episode' as the fundamental entity of EHRs, which
comprises all the data being recorded during the ICU stay. Each episode also
contains the demographic information of the patient, ward transferring origin
and destination within a hospital and diagnosis information. It allows us to 
link episode data from a single patient across the entire multi-centre database. 
The key dates and times are recorded as follow. 
<img src="cchic_timeline.png" width=480 height=360 />
```{r, message=FALSE, warning=FALSE}
library(cleanEHR)
data("sample_ccd")

# Extract all non-longitudinal data (demographic, time, survival status, diagnosis)
dt <- ccd_demographic_table(ccd, dtype=TRUE)
```
`ccd_demographic_table` function returns a table of all non-time-varying fields alongside with 
several derived fields -- the fields that are not directly recorded in the original data. 
Each row of the table is a unique admission, and every column is a non-time-varying data field. 

* `pid`: unique patient ID derived from NHS number or PAS number.
* `AGE`: date of birth - unit admission time
* `spell`: see *Spell*

```{r}
print(dt[1:3, ])
```

Data missing are caused by many reasons. We have to understand that in such a large database, 
the data quality varies. In the anonymised dataset, data can be missing due to security reason. 
Missing data are either NULL or "NA" depending on the field data type. 

### Demographic data
Every patient in England has a unique NHS number and PAS (Patient Admission
System) number. These can be used to identify a unique patient. Other
demographic information can also be found in the CCHIC dataset such as age, sex,
GP code, postcode and so on. Most of the demographic data will be removed,
pseudonimised, or modified in the anonymised dataset to protect the patient
privacy. 


### Spell
Ward transferring within or beyond ICUs are counted as different episodes respectively. 
In some research, one may be interested in looking at the sickness development and the treatment 
history beyond each ICU stays. A spell includes number of episodes from a
unique patient which occured in a user defined period One can link episodes by
spell ID. 
```{r}
head(ccd_unique_spell(ccd, duration=1)[, c("episode_id", "spell")])
```

### Diagnosis data
Instead of using free text, we adopted [ICNARC diagnosis](https://www.icnarc.org/Our-Audit/Audits/Cmp/Resources/Icm-Icnarc-Coding-Method) 
code system to record the diagnosis. The full ICNARC code is a five digit number separated 
with dots. From left to right each digit represents a higher category of diagnosis. Due 
to the privacy concerns, in the anonysmised dataset, last two digits will be removed. 
One may use the function `icnarc2diagnosis` to look up the diagnosis code.
```{r}
icnarc2diagnosis("1.1")
icnarc2diagnosis("1.1.4")
icnarc2diagnosis("1.1.4.27.1")
```

### Longitudinal data
CCHIC measures 154 longitudinal data. The full list of longitudinal data are shown below. 
<img src="item_ref.png" width=480 height=360 />

We can also easily plot the data from a single selected admission. 
```{r, fig.width=10, fig.height=11, out.width='700px', results='hide', message=FALSE, warning=FALSE}
plot_episode(ccd@episodes[[7]], c("h_rate",  "bilirubin", "fluid_balance_d"))
```
---
title: "Data cleaning and wrangling with cleanEHR"
author: David Perez Suarez & Sinan Shi
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 3
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Data cleaning and wrangling with cleanEHR}
  \usepackage[utf8]{inputenc}

---


# Preparation
[data.table](https://CRAN.R-project.org/package=data.table) package is the backbone of cleanEHR package. 
You can find in the above link some useful information and tutorial, if you are not familiar with
`data.table`.

### Load data
```{r}
library(cleanEHR)
data("sample_ccd")
```


### Inspect individual episode 
There are 263 fields which covers patient demographics, physiology, laboratory,
and medication information. Each field has 2 labels, NHIC code and short name.
There is a function `lookup.items()` to look up the fields you need.
`lookup.items()` function is case insensitive and allows fuzzy search.
```
# searching for heart rate
lookup.items('heart') # fuzzy search

+-------------------+--------------+--------------+--------+-------------+
|     NHIC.Code     |  Short.Name  |  Long.Name   |  Unit  |  Data.type  |
+===================+==============+==============+========+=============+
| NIHR_HIC_ICU_0108 |    h_rate    |  Heart rate  |  bpm   |   numeric   |
+-------------------+--------------+--------------+--------+-------------+
| NIHR_HIC_ICU_0109 |   h_rhythm   | Heart rhythm |  N/A   |    list     |
+-------------------+--------------+--------------+--------+-------------+

```


# Non-longitudinal Data 
`ccd_demographic_table()` can generate a `data.table` that contains all the
non-longitudinal variables. A demonstration of how to do some work on a subset
of data.  
```{r, fig.width=10, fig.height=6, out.width='600px', results='hide', message=FALSE, warning=FALSE}
# contains all the 1D fields i.e. non-longitudinal
tb1 <- ccd_demographic_table(ccd)

# filter out all dead patient. (All patients are dead in the dataset.)
tb1 <- tb1[DIS=="D"]

# subset variables we want (ARSD = Advanced respiratory support days,
# apache_prob = APACHE II probability)
tb <- tb1[, c("SEX", "ARSD", "apache_prob"), with=FALSE]
tb <- tb[!is.na(apache_prob)]

# plot
library(ggplot2)
ggplot(tb, aes(x=apache_prob, y=ARSD, color=SEX)) + geom_point()

```


# Longitudinal data 
## Longitudinal table structure: `ccTable`

To deal with longitudinal data, we need to first to transform it into a table format. 
cleanEHR provides a refclass `ccTable`. There are several key components in the `ccTable`
structure. 

* `torigin`: the `ccRecord` dataset will be converted into a table format, where each row is a data point and each column is a field and pivoted by `time`, `site`, and `eipisode_id`. 
* `tclean`: Same structure like the `torigin` but the values are modified with the cleaning process. 
* filters: `filter_range`, `filter_categories`, `filter_nodata`, `filter_missingness`.
* imputation: filling missing data. 

### Create a `cctable`
First we need to prepare a simple YAML configuration file. YAML is a human freindly 
data serialization standard, see [YAML](http://yaml.org/). 
The first level item is
CCHIC code, see `lookup.items()`.  We suggest users to write the short name and
long name (dataItem) to avoid confusion, though the both names will not be
taken into account in the process. We selected three items, heart rate
(longitudinal), Systolic arterial blood pressure (longitudinal), and sex
(non-longitudinal).

```{r}
# To prepare a YAML configuration file like this. You write the following text
# in a YAML file. 
conf <- "
NIHR_HIC_ICU_0108:
  shortName: hrate
NIHR_HIC_ICU_0112:
  shortName: bp_sys_a
  dataItem: Systolic Arterial blood pressure - Art BPSystolic Arterial blood pressure
NIHR_HIC_ICU_0093:
   shortName: sex
"
library(yaml)
conf <- yaml.load(conf)
```

```{r}
# conf is the full path of the YAML configuration.
tb <- create_cctable(ccd, conf, freq=1)
print(tb$torigin) # the table
```
In this table we can find the following columns, 

* time: number of hours from the unit admission. Since we set the `freq`=1, the cadence between rows is always 1 hour. 
* site, episode_id: combine these two columns will give you a unique admission.
* fields: three selected fields. 
* extra fields: depending on the variable we choose, some extra information are given. 


### Get the mean heart rate of each patient. 
```{r}
tb$tclean[, mean(NIHR_HIC_ICU_0108, na.rm=TRUE), by=c("site", "episode_id")]
```

# Data cleaning with `ccTable`
## Numerical range filter
The numerical range filter can only be applied on variables. 
We envisaged three different cases for the numerical ranges -- values that are impossible, e.g. 
negative heart rate; possible but unlikely, e.g. heart rate of 200; within a normal range. The 
filter will label all these scenarios using "red", "amber", "green" respectively. The definition 
of these ranges can be configured by users based on their judgement and the purpose of research. 
Note, from "red" to "green", the next range must be a subset of the previous range.

In the following section, we would like to apply a range filter to heart rate by modifying the previous 
YAML configuration file. 

```{r}
conf <- "NIHR_HIC_ICU_0108:
  shortName: h_rate
  dataItem: Heart rate
  range:
    labels:
      red: (0, 300)
      amber: (11, 150]
      green: (50, 100]
    apply: drop_entry
NIHR_HIC_ICU_0112:
  shortName: bp_sys_a
  dataItem: Systolic Arterial blood pressure - Art BPSystolic Arterial blood pressure
NIHR_HIC_ICU_0093:
   shortName: sex
   category:
      M: male
      F: female
      m: male
      f: female
"
conf <- yaml.load(conf)
```

```{r}
tb <- create_cctable(ccd, conf, freq=1)
tb$filter_range("amber") # chose only the entry with amber
tb$apply_filters() # apply the filter to the clean table
```
Now let's see the effect on the cleaned data `tclean`

```{r, fig.width=12, fig.height=12, out.width='700px', results='hide', message=FALSE, warning=FALSE}
cptb <- rbind(cbind(tb$torigin, data="origin"), 
              cbind(tb$tclean, data="clean"))

ggplot(cptb, aes(x=time, y=NIHR_HIC_ICU_0108, color=data)) + 
  geom_point(size=1.5) + facet_wrap(~episode_id, scales="free_x")
```

In the case of changing the filter range from amber to green, 
```{r}
#tb$reset() # reset the all the filters first.
tb$filter_range("green")
tb$apply_filters()
```
```{r, fig.width=12, fig.height=12, out.width='700px', results='hide', message=FALSE, warning=FALSE}
cptb <- rbind(cbind(tb$torigin, data="origin"), 
              cbind(tb$tclean, data="clean"))

ggplot(cptb, aes(x=time, y=NIHR_HIC_ICU_0108, color=data)) + 
  geom_point(size=1.5) + facet_wrap(~episode_id, scales="free_x")
```


## Categorical data filter
The purpose of categorical data filter is to remove the unexpected categorical data. 
We can extend the previous configuration file as such,  
```{r}
conf <- "NIHR_HIC_ICU_0108:
  shortName: h_rate
  dataItem: Heart rate
NIHR_HIC_ICU_0112:
  shortName: bp_sys_a
  dataItem: Systolic Arterial blood pressure - Art BPSystolic Arterial blood pressure
NIHR_HIC_ICU_0093:
   shortName: sex
   category:
    levels:
      M: male
      F: female
      m: male
      f: female
    apply: drop_entry
"
conf <- yaml.load(conf)

# Try to modify the original data
tb$torigin$NIHR_HIC_ICU_0093[1] <- "ERROR"

tb$reload_conf(conf)  # change configuration file
tb$filter_categories() 
tb$apply_filters() 
```
There is one error gender introduced in the sex field. After the filtering 
process, the error entry is substitute by NA. 
```{r}
unique(tb$torigin$NIHR_HIC_ICU_0093)
unique(tb$tclean$NIHR_HIC_ICU_0093)
```


## Missingness filter
In some cases, we wish to exclude episodes where the data is too scarce. There are 
three components in the missingness filter. In the following example, we arbitrarily
name the filter "daily". We gave 24 hours interval and 70% accepting rate. It is to say
in any 24 hours interval, if the heart rate missing rate is higher than 30%, we will 
exclude the entire episode. Note, the unit `labels: daily: 24` number of rows instead of 
hours. It represents 24 hours because the cadence of the `ccTable` is 1 hour.

```{r}
conf <- "NIHR_HIC_ICU_0108:
  shortName: h_rate
  dataItem: Heart rate
  missingness:
    labels:
      daily: 24
    accept_2d:
      daily: 70
    apply: drop_episode
NIHR_HIC_ICU_0112:
  shortName: bp_sys_a
  dataItem: Systolic Arterial blood pressure - Art BPSystolic Arterial blood pressure
NIHR_HIC_ICU_0093:
   shortName: sex
"
conf <- yaml.load(conf)

tb$reload_conf(conf)  # change configuration file
tb$filter_missingness() 
tb$apply_filters()

# episodes in the original data table 
unique(paste(tb$torigin$site, tb$torigin$episode_id))
# episodes in the cleaned data table
unique(paste(tb$tclean$site, tb$tclean$episode_id))
```

## Nodata filter
Similarly, we can setup the no data filter as following. Here it means, 
drop the entire episode if no hear rate data is found. 

```yaml
NIHR_HIC_ICU_0108:
  shortName: h_rate
  dataItem: Heart rate
  no_data:
    apply: drop_episode
NIHR_HIC_ICU_0112:
  shortName: bp_sys_a
  dataItem: Systolic Arterial blood pressure - Art BPSystolic Arterial blood pressure
NIHR_HIC_ICU_0093:
   shortName: sex
```

## Run all filters together
To wrap up, we can put all the above stated filter configurations together in the 
YAML file and run the filter together. 

```{r}
conf <- "NIHR_HIC_ICU_0108:
  shortName: h_rate
  dataItem: Heart rate
  range:
    labels:
      red: (0, 300)
      amber: (11, 150]
      green: (50, 100]
    apply: drop_entry
  missingness:
    labels:
      daily: 24
    accept_2d:
      daily: 70
    apply: drop_episode
  nodata:
    apply: drop_episode
NIHR_HIC_ICU_0112:
  shortName: bp_sys_a
  dataItem: Systolic Arterial blood pressure - Art BPSystolic Arterial blood pressure
NIHR_HIC_ICU_0093:
   shortName: sex
   category:
    levels:
      M: male
      F: female
      m: male
      f: female
    apply: drop_entry
"
conf <- yaml.load(conf)

# Method 1
tb <- create_cctable(ccd, conf, freq=1)
tb$filter_range("amber")
tb$filter_missingness()
tb$filter_nodata()
tb$filter_categories()
tb$apply_filters() 

tb$reset() # reset

# Method 2
#tb$clean()
```


## Imputation

We provide the `impute()` to interpolate the missing data. For each missing value, 
the interpolation will be only based on the nearby values which are specified by 
`lead` and `lag` arguments. `lead` suggests the number of previous values and `lag` suggests 
the number of later values. The corresponding time will be related to the 
`freq` you set for the `ccTable`, e.g. `lead: 2` means previous 4 hours when `freq=0.5`.
One can also set the `fun` to determine the interpolation function.
The imputation step usually should be carried out after filtering, otherwise
imputation will take values that are out of range its into account.

One needs to be always careful when impute the data to make the best 
trade-off between usefulness and correctness. The interpolation methods should 
be carried out wisely based on the characteristics of the variable. A good 
overview of how to deal with the missing data can be found 
[(Salagodo et al. 2016)](https://link.springer.com/content/pdf/10.1007%2F978-3-319-43742-2_13.pdf)
If you are not sure about the characteristics of the variable, we would 
suggest you to keep the window size small and use `median` as the interpolation 
function. 


```{r}
# Initialise the simulated ccRecord
hr <- c(rep(80, 10), rep(NA, 10), rep(90, 10), NA, NA, rep(90, 10), rep(NA, 10), 180, NA, NA, 
        rep(90, 10), 180, NA, 0, NA, NA, rep(60, 10))
# hr <- hr + runif(length(hr)) * 15 # adding noise if needed. 
data <- data.frame(time=as.numeric(seq(hr)), item2d=hr)
rec <- ccRecord()+new.episode(list(NIHR_HIC_ICU_0108=data))

# Prepare the plotting function
library(data.table)
plot_imputation <- function() {
cptb <- data.table(episode_id=as.integer(tb$torigin$episode_id), 
                   time=tb$torigin$time, origin=tb$torigin$NIHR_HIC_ICU_0108, 
                   clean=tb$tclean$NIHR_HIC_ICU_0108)

ggplot(cptb, aes(x=time)) + 
  geom_point(size=5, shape=16, aes(y=origin), colour="red") + 
  geom_point(size=2, aes(y=clean)) + 
  geom_line(aes(y=clean)) + 
  scale_x_continuous(minor_breaks = seq(length(hr)))+ 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5), 
        panel.grid.major = element_line(colour="grey", size=0.5))
}



```
**Example 1**: median interpolation with a window [-2, 2]
```{r, fig.width=12, fig.height=12, out.width='700px', results='hide', message=FALSE, warning=FALSE}
# mock the configuration YAML
conf <- "NIHR_HIC_ICU_0108:
  shortName: h_rate
  dataItem: Heart rate
  missingness:
    impute:
      lead: 2 # 2 previous values 
      lag: 2  # 2 later values
      fun: median # missing value filled by the median of 2 previous and 2 later values. 
  nodata:
    apply: drop_episode
"

conf <- yaml.load(conf)
tb <- create_cctable(rec, conf, freq=1)
tb$imputation()
plot_imputation()
```

**Example 2**: increase the window size to [-10, 10]
We can increase the window size to fill more data, 
```{r, fig.width=12, fig.height=12, out.width='700px', results='hide', message=FALSE, warning=FALSE}
conf <- "NIHR_HIC_ICU_0108:
  shortName: h_rate
  dataItem: Heart rate
  missingness:
    impute:
      lead: 10
      lag: 10
      fun: median
  nodata:
    apply: drop_episode
"

rec <- ccRecord()+new.episode(list(NIHR_HIC_ICU_0108=data))

conf <- yaml.load(conf)
tb <- create_cctable(rec, conf, freq=1)
tb$imputation()
plot_imputation()
```

**Example 3**: use `mean` as the interpolation function. 
```{r, fig.width=12, fig.height=12, out.width='700px', results='hide', message=FALSE, warning=FALSE}
conf <- "NIHR_HIC_ICU_0108:
  shortName: h_rate
  dataItem: Heart rate
  missingness:
    impute:
      lead: 10
      lag: 10
      fun: mean
  nodata:
    apply: drop_episode
"

rec <- ccRecord()+new.episode(list(NIHR_HIC_ICU_0108=data))

conf <- yaml.load(conf)
tb <- create_cctable(rec, conf, freq=1)
tb$imputation()
plot_imputation()

```

**Advanced Example**: Use a self-defined function. 
```{r, fig.width=12, fig.height=12, out.width='700px', results='hide', message=FALSE, warning=FALSE}
conf <- "NIHR_HIC_ICU_0108:
  shortName: h_rate
  dataItem: Heart rate
  missingness:
    impute:
      lead: 40
      lag: 40
      fun: myfun
  nodata:
    apply: drop_episode
"
# Define my own interpolation function. 
# We use piecewise polynomial interpolation spline here for 
# the demonstration purpose. 
myfun <- function(x) {
    return(splinefun(x)(ceiling(length(x)/2)))
}


conf <- yaml.load(conf)
tb <- create_cctable(rec, conf, freq=1)
tb$imputation()
plot_imputation()
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data.checklist}
\alias{data.checklist}
\title{This a reference table of NHIC data items.}
\description{
This a reference table of NHIC data items.
}
\author{
Sinan Shi \email{s.shi@ucl.ac.uk}
}
\keyword{data}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{methods}
\name{+,ccRecord,ccRecord-method}
\alias{+,ccRecord,ccRecord-method}
\title{Combine two ccRecord objects}
\usage{
\S4method{+}{ccRecord,ccRecord}(e1, e2)
}
\arguments{
\item{e1}{ccRecord-class}

\item{e2}{ccRecord-class}
}
\value{
ccRecord-class
}
\description{
Combine two ccRecord objects
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cchic_xml.R
\name{xml2Data}
\alias{xml2Data}
\title{Convert the XML file to ccRecord}
\usage{
xml2Data(file, select.episode = NULL, quiet = TRUE, xml = NULL,
  file_origin = "NA", parse_time = Sys.time())
}
\arguments{
\item{file}{character string. The path of XML file.}

\item{select.episode}{integer vector. Load only a selected number of episodes. 
It is NULL by default which loads all the episodes in a file.}

\item{quiet}{logical. Switch on/off the progress bar.}

\item{xml}{XML object. Usually not needed.}

\item{file_origin}{character string. The XML file name. The file name will be 
extracted automatically while argument xml is NULL.}

\item{parse_time}{POSIXct. By default is the time of the execution of this function.}
}
\value{
ccRecord-class
}
\description{
Convert the XML file to ccRecord. For more details, see ccRecord-class.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{methods}
\name{+,ccRecord,list-method}
\alias{+,ccRecord,list-method}
\title{Adding a list of ccEpisode to ccRecord}
\usage{
\S4method{+}{ccRecord,list}(e1, e2)
}
\arguments{
\item{e1}{ccRecord}

\item{e2}{a list of ccEpisode objects}
}
\value{
ccRecord
}
\description{
Adding a list of one or multiple ccEpisode objects to a
ccRecord object, the information table (infotb) will be updated automatically.
It is the more efficient way to add multiple ccEpisode objects.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ccd}
\alias{ccd}
\title{Synthetic example dataset}
\description{
This dataset has the same data strucutre that of the CCHIC data, though
the data are synthetic.
}
\keyword{data}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{stname2code}
\alias{stname2code}
\title{Convert short names to NHIC codes}
\usage{
stname2code(stname)
}
\arguments{
\item{stname}{character short names of data item h_rate}
}
\value{
NIHC code character such as NIHR_HIC_ICU_0108
}
\description{
Convert short names to NHIC codes
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{is.demographic}
\alias{is.demographic}
\title{Check if the item NHIC code or short name belongs to the demographic
category.}
\usage{
is.demographic(item_name)
}
\arguments{
\item{item_name}{character the NHIC code or the short name}
}
\value{
logical
}
\description{
Check if the item NHIC code or short name belongs to the demographic
category.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{ccTable_apply_filters}
\alias{ccTable_apply_filters}
\title{Apply all the setup filters.}
\arguments{
\item{warnings}{logical value to indicate more or less messages with an 
default value TRUE.}
}
\description{
Once filters are applied, the processed data will be stored in tclean. Note, 
running filtering function before apply_filters is necessary. This function 
will have no effect on tclean if no filter is ran prior.
Filters will decide to preserve or remove particular entries or episodes.
}
\examples{
\dontrun{
tb <- create_cctable(ccd, conf, 1)
tb$range_filter() 
tb$apply_filter() # apply only the range filter ragardless of the conf. 
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{methods}
\name{[[,ccRecord-method}
\alias{[[,ccRecord-method}
\title{Subsetting a ccRecord object and return a list of ccEpisode objects.}
\usage{
\S4method{[[}{ccRecord}(x, i)
}
\arguments{
\item{x}{ccRecord-class}

\item{i}{integer vector}
}
\description{
Subsetting a ccRecord object and return a list of ccEpisode objects.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demographics.R
\name{ccd_demographic_table}
\alias{ccd_demographic_table}
\title{Create the demographic tables, which includes all non-time-varying variables.}
\usage{
ccd_demographic_table(record, dtype = TRUE)
}
\arguments{
\item{record}{ccRecord-class}

\item{dtype}{logical column will be type aware, else all in character.}
}
\description{
The data type of each column is in its corresponding data 
type.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cchic_xml.R
\name{xmlLoad}
\alias{xmlLoad}
\title{load xml clinical data}
\usage{
xmlLoad(file)
}
\arguments{
\item{file}{character string. The path of the XML file.}
}
\value{
the root of the xml data.
}
\description{
load xml clinical data
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{total.data.point}
\alias{total.data.point}
\title{Return total data point of the ccRecord object.}
\usage{
total.data.point(ccd)
}
\arguments{
\item{ccd}{ccRecord-class}
}
\description{
Return total data point of the ccRecord object.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{create_cctable}
\alias{create_cctable}
\title{Create a ccTable object}
\usage{
create_cctable(rec, conf = NULL, freq = 1)
}
\arguments{
\item{rec}{ccRecord}

\item{conf}{either the path of YAML configuration file or the configuration}

\item{freq}{a unique sampling frequency (in hours) for all variables. e.g. if freq is set to 
1, each row in ccTable will represent a record of one hour.}
}
\value{
ccTable
}
\description{
Re-arrange the ccRecord object to table format where each column stands 
for a variable and each row a record data point. The number of rows will 
depend on the sampling frequency set in this function. If the original data
has a higher recording frequency than the set frequency (freq), the closest 
data point will be taken. It is suggested the `freq` should not be set 
lower than the maximum recording frequency in the original dataset.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{is.physiology}
\alias{is.physiology}
\title{Check if the item NHIC code or short name belongs to the physiology
category.}
\usage{
is.physiology(item_name)
}
\arguments{
\item{item_name}{character the NHIC code or the short name}
}
\value{
logical
}
\description{
Check if the item NHIC code or short name belongs to the physiology
category.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ITEM_REF}
\alias{ITEM_REF}
\title{Field reference table}
\description{
Field reference table
}
\keyword{data}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{getEpisodePeriod}
\alias{getEpisodePeriod}
\title{Get the length of stay based on the first and the last data point.}
\usage{
getEpisodePeriod(e, unit = "hours")
}
\arguments{
\item{e}{ccEpisode object.}

\item{unit}{character string.  Units in which the results are desired. Can be abbreviated.}
}
\value{
length of stay
}
\description{
Get the length of stay based on the first and the last data point.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deltaTime.R
\name{deltaTime}
\alias{deltaTime}
\title{Convert calendar date-time to the time difference comparing to the ICU
admission time.}
\usage{
deltaTime(record, pseudotime = FALSE, units = "hours", tdiff = FALSE)
}
\arguments{
\item{record}{ccRecord}

\item{pseudotime}{logical If pseudotime is set to be TRUE, then the
admission and discharge time will be set as the earliest and latest data stamp
in the record.}

\item{units}{units of delta time, which can be "hours", "mins", "days".}

\item{tdiff}{if false the delta time will be written in numeric format.}
}
\description{
Convert calendar date-time to the time difference comparing to the ICU
admission time.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{class}
\name{ccEpisode-class}
\alias{ccEpisode}
\alias{ccEpisode-class}
\title{The S4 class which holds data of a single episode.}
\description{
The S4 class which holds data of a single episode.
}
\section{Fields}{

\describe{
\item{\code{site_id}}{character string. Site ID, if presented, otherwise "NA".}

\item{\code{episode_id}}{character string. Episode ID, if presented, otherwise "NA".}

\item{\code{nhs_number}}{character string. NHS number, if presented, otherwise "NA".}

\item{\code{pas_number}}{character string. PAS number, if presented, otherwise "NA".}

\item{\code{parse_file}}{character string. The source XML file. If the source is not a file then "NA".}

\item{\code{t_admission}}{POSIXct. Time of Admission to the ICU, if presented, otherwise NA.}

\item{\code{t_discharge}}{POSIXct. Time of discharge of the ICU, if presented, otherwise NA.}

\item{\code{parse_time}}{POSIXct. Parse time.}

\item{\code{data}}{A list which holds all the data of this episode which is indexed by NIHIC code.}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{ccTable_clean}
\alias{ccTable_clean}
\title{Apply all the filters}
\description{
All the filters in configuration will be applied to create the 
clean dataset. The filters include range, categories, missingness, 
no_data.
}
\examples{
\dontrun{
tb <- create_cctable(ccd, conf, 1)
tb$clean()
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{icnarc_table}
\alias{icnarc_table}
\title{ICNARC diagnosis reference table}
\description{
ICNARC diagnosis reference table
}
\references{
\url{https://www.icnarc.org/Our-Audit/Audits/Cmp/Resources/Icm-Icnarc-Coding-Method}
}
\keyword{data}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{demographic.data.completeness}
\alias{demographic.data.completeness}
\title{Create a demographic completeness table (in pander table)}
\usage{
demographic.data.completeness(demg, names = NULL, return.data = FALSE)
}
\arguments{
\item{demg}{data.table the demographic data table created by ccd_demographic_table()}

\item{names}{short name of selected items}

\item{return.data}{logical return the table if TRUE}
}
\description{
Create a demographic completeness table (in pander table)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.categorical.R
\name{ccTable_filter_categories}
\alias{ccTable_filter_categories}
\title{Categorical data filter}
\description{
Categorical variables only allow a set of values to appear in the variable
. Due to various reasons, a categorical variable may contain values that are not 
standard. The allowed values can be set in the YAML configuration while initialising 
the ccTable (see ccTable-class, create_cctable).  
In the following example, we can see how to set up the categorical filter 
for the variable dead_icu (NIHR_HIC_ICU_0097) which only allows its value to 
be A, D, E.
}
\examples{
\dontrun{
# Example for categorical filter setup in the YAML configuration
NIHR_HIC_ICU_0097:
 category:
   levels:
     A: Alive
     D: Dead
     E: Alive - not discharged
   apply: drop_entry

# Run the filter on ccTable ct
ct$filter_categories() # run the filter
ct$apply_filters()     # apply the filter and create the clean table
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{is.drugs}
\alias{is.drugs}
\title{Check if the item NHIC code or short name belongs to the drugs 
category.}
\usage{
is.drugs(item_name)
}
\arguments{
\item{item_name}{character the NHIC code or the short name}
}
\value{
logical
}
\description{
Check if the item NHIC code or short name belongs to the drugs 
category.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{which.classification}
\alias{which.classification}
\title{Identify the classification - classification1}
\usage{
which.classification(item_name)
}
\arguments{
\item{item_name}{NHIC code or the short name}
}
\value{
character the item classification
}
\description{
Identify the classification of a given item code or short
name. Classification1 has 5 labels: 
[1] "Demographic", [2] "Physiology" 
[3] "Drugs" [4] "Nursing_other" [5] "Laboratory"
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{xml.file.duration.plot}
\alias{xml.file.duration.plot}
\title{plot the duration of XML files.}
\usage{
xml.file.duration.plot(ccd)
}
\arguments{
\item{ccd}{ccRecord-class}
}
\description{
plot the duration of XML files.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{stname2longname}
\alias{stname2longname}
\title{Convert short names to long names.}
\usage{
stname2longname(stname)
}
\arguments{
\item{stname}{character short names of data item h_rate}
}
\value{
longname character such as "heart rate"
}
\description{
Convert short names to long names.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{ccTable_reset}
\alias{ccTable_reset}
\title{Reset the ccTable}
\description{
Restore the object to its initial status. All the filters, quality and the 
cleaned table will be removed.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{table1}
\alias{table1}
\title{Produce the item specified table one.}
\usage{
table1(demg, names, return.data = FALSE)
}
\arguments{
\item{demg}{demographic table created by ccd_demographic_table()}

\item{names}{character string. Short names of data items, e.g. h_rate.}

\item{return.data}{logical, FALSE: printing the pander table, TRUE: return the table but not print out the pander table.}
}
\value{
if return.data is TRUE, return data.table
}
\description{
Produce the item specified table one.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.missingness.R
\name{ccTable_filter_nodata}
\alias{ccTable_filter_nodata}
\title{No data filter}
\description{
Remove the episode when a particular field is not presented.
It need to be set up in the YAML configuration file.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\name{for_each_episode}
\alias{for_each_episode}
\title{loop over all episodes of a ccRecord object}
\usage{
for_each_episode(record, fun)
}
\arguments{
\item{record}{ccRecord}

\item{fun}{function}
}
\description{
loop over all episodes of a ccRecord object
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{code2stname}
\alias{code2stname}
\title{Convert NHIC codes to the short names}
\usage{
code2stname(code)
}
\arguments{
\item{code}{character NIHC code, e.g. NIHR_HIC_ICU_0108}
}
\value{
shortname character e.g. h_rate
}
\description{
Convert NHIC codes to the short names
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{file.summary}
\alias{file.summary}
\title{Produce a file summary table}
\usage{
file.summary(ccd)
}
\arguments{
\item{ccd}{ccRecord-class}
}
\value{
data.table
}
\description{
Produce a file summary table
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{ccTable_create_cctable}
\alias{ccTable_create_cctable}
\title{Create a ccTable object}
\description{
This is a member function of ccTable-class. Using create_cctable is a safer and 
easier way to create the ccTable. See create_cctable.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{as.number}
\alias{as.number}
\title{Convert standard IDs to numbers (character) which can be used for indexing.}
\usage{
as.number(obj)
}
\arguments{
\item{obj}{a StdId object.}
}
\description{
Convert standard IDs to numbers (character) which can be used for indexing.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\name{new.episode}
\alias{new.episode}
\title{Create a new episode}
\usage{
new.episode(lt = list(), parse_file = "NA", parse_time = as.POSIXct(NA))
}
\arguments{
\item{lt}{is a list}

\item{parse_file}{the file location from which the episode comes from.}

\item{parse_time}{the parse date and time of the episode.}
}
\value{
ccEpisode object
}
\description{
create a new ccEpisode object by given the episode data as a
list. The list should be organised in data items and indexed with NIHC code,
e.g. NIHR_HIC_ICU_0108.
}
\examples{
eps <- list()
eps[["NIHR_HIC_ICU_0018"]] <- data.frame(time=seq(10), rep(70, 10))
new.episode(eps)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{samplerate2d}
\alias{samplerate2d}
\title{Produce a pander table of sample rate of longitudinal data.}
\usage{
samplerate2d(cctb)
}
\arguments{
\item{cctb}{ccTable-class, see create.cctable().}
}
\description{
Produce a pander table of sample rate of longitudinal data.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{long2stname}
\alias{long2stname}
\title{Convert long names to short names.}
\usage{
long2stname(l)
}
\arguments{
\item{l}{long name such as "heart rate"}
}
\value{
short name character such as "h_rate"
}
\description{
Convert long names to short names.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{reallocateTimeRecord}
\alias{reallocateTimeRecord}
\title{Propagate a numerical delta time interval record.}
\usage{
reallocateTimeRecord(record, delta = 0.5)
}
\arguments{
\item{record}{ccRecord}

\item{delta}{time frequency in hours}
}
\description{
Propagate a numerical delta time interval record.
}
\details{
when discharge time and admission time are missing, the latest  and
the earliest data time stamp will be used instead.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{methods}
\name{plot_episode,ccEpisode,missing-method}
\alias{plot_episode,ccEpisode,missing-method}
\title{Episode chart default fields}
\usage{
\S4method{plot_episode}{ccEpisode,missing}(r)
}
\arguments{
\item{r}{ccEpisode-class}
}
\description{
Episode chart default fields
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{extract_info}
\alias{extract_info}
\title{Extract information from data.checklist}
\usage{
extract_info()
}
\value{
list of time [data.frame(id, idt)], meta [data.frame(id, idmeta)], 
        nontime [numeric], MAX_NUM_NHIC
}
\description{
Extract information from data.checklist
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.range.R
\name{ccTable_filter_range}
\alias{ccTable_filter_range}
\title{Numerical range filter}
\arguments{
\item{select}{the range label - "red", "amber", "green"
If I give "yellow to select, it means I only want the values which is
labeled as "yellow" to be in the clean table.}
}
\description{
Range filter can only be applied on numerical fields. 
For those fields which requires a range filter to be applied, 
one needs to set a series ranges from the broadest to the narrowest in
the YAML configuration. We can set three levels (labels) of ranges, red, amber, and
green.  It is also OK to set only one range instead of three. 
The range filter will first assign a label to every data entry.
}
\details{
The range in the YAML configuration file can be (l, h), [l, h], (l, h], [h, l) 
standing for close, open and half open intervals.
}
\examples{
\dontrun{
# YAML example
NIHR_HIC_ICU_0108:
   shortName: h_rate
   dataItem: Heart rate
   range:
    labels:
     red: (0, 300)     # broader
     amber: (11, 170) 
     green: (60, 100)  # narrower
    apply: drop_entry
# apply range filter on ccTable ct
ct$filter_range("yellow")
ct$apply_filters
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{ccd_select_table}
\alias{ccd_select_table}
\title{Create the table for ccTable from ccRecord}
\usage{
ccd_select_table(record, items_opt = NULL, items_obg = NULL, freq,
  return_list = FALSE)
}
\arguments{
\item{record}{ccRecord}

\item{items_opt}{character vectors. Items (HIC code) selected in item_opt are optional items, which will be automatically 
filled when item is missing.}

\item{items_obg}{obligatory items that is obligatory; Any episode that does not contain
item in this vector will be removed.}

\item{freq}{numeric cadence in hour.}

\item{return_list}{logical if TRUE return as a list.}
}
\value{
data.table
}
\description{
Create the table for ccTable from ccRecord
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{icnarc2diagnosis}
\alias{icnarc2diagnosis}
\title{Convert ICNARC codes to diagnosis (text)}
\usage{
icnarc2diagnosis(icnarc, surgery = TRUE, levels = NULL)
}
\arguments{
\item{icnarc}{the ICNARC code, e.g. 1.1.1.1.1}

\item{surgery}{T/F with or without surgical information}

\item{levels}{category level, from [1 - 5]. TODO level 4.}
}
\value{
character ICNARC diagnosis
}
\description{
NOTE: There are still ~600 code missing. see issue #133
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{methods}
\name{+,ccRecord,ccEpisode-method}
\alias{+,ccRecord,ccEpisode-method}
\title{Adding one ccEpisode object to a ccRecord}
\usage{
\S4method{+}{ccRecord,ccEpisode}(e1, e2)
}
\arguments{
\item{e1}{ccRecord-class}

\item{e2}{ccEpisode-class}
}
\value{
ccRecord-class
}
\description{
Adding one ccEpisode object to a ccRecord
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{xmlTime2POSIX}
\alias{xmlTime2POSIX}
\title{Convert time from xml to POSIX format.}
\usage{
xmlTime2POSIX(xml.time, allow = FALSE)
}
\arguments{
\item{xml.time}{character. Time in XML format such as 2014-02-01T03:00:00}

\item{allow}{logical. Wrong format will be accepted when \code{allow} is set
to be TRUE and NA will be the return value, otherwise return error. 
It is useful while dealing with pseudonymous data where the time format is
not presented correctly.}
}
\description{
Convert the XML time The XML time format to POSIXct.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{create2dclean}
\alias{create2dclean}
\title{Clean table - low memory}
\usage{
create2dclean(record, config, freq = 1, nchunks = 1)
}
\arguments{
\item{record}{ccRecord}

\item{config}{the full path of the YAML configuration file}

\item{freq}{table cadence}

\item{nchunks}{integer number. The larger the nchunks the less memory
requirement.}
}
\value{
A cleaned 2d wide table
}
\description{
The cleaning process is specified by the YAML configuration. All the filters
presented in the configuration will be applied. It returns only the cleaned
data. However all the data quality information will be lost. This function
is useful when the memory is not sufficiently enough to hold all the
information.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\docType{class}
\name{ccTable-class}
\alias{ccTable}
\alias{ccTable-class}
\title{Process the EHR data in table format}
\description{
ccRecord data are re-arranged into tables where the columns stands for 
data fields (e.g. heart rate, blood pressure) and the rows stands for 
each data record within a unique cadence. See ccTable_create_cctable.
ccTable is the data processing platform. It stores both original data 
and processed data alongside with the process details. It also contains 
various commonly used data filters.
}
\section{Fields}{

\describe{
\item{\code{record}}{ccRecord.}

\item{\code{conf}}{the YAML style configuration.}

\item{\code{torigin}}{the original data table.}

\item{\code{tclean}}{the data table after cleaning processes.}

\item{\code{dfilter}}{list contains data filtering information.}

\item{\code{dquality}}{list contains data quality.}

\item{\code{summary}}{list}

\item{\code{base_cadence}}{the base cadence is specified in hours}
}}
\section{Methods}{

\describe{
\item{\code{apply_filters(warnings = TRUE)}}{Apply all filters specified in the configuration to update the clean
table (tclean)}

\item{\code{create_table(freq)}}{Create a table contains the selected items in the conf with a given
frequency (in hour)}

\item{\code{export_csv(file = NULL)}}{Export the cleaned table to a CSV file.}

\item{\code{filter_categories()}}{Check individual entries if they are the in the categories specified
in conf.}

\item{\code{filter_missingness(recount = FALSE)}}{filter out the where missingness is too low.}

\item{\code{filter_nodata()}}{Exclude episodes when no data is presented in certain fields}

\item{\code{imputation()}}{Filling missing data to a time series data by performing a given imputation
method on a selected window period nearby the missing data.}

\item{\code{reload_conf(conf)}}{reload yaml configuration.}
}}
\examples{
rec <- ccRecord()
cctable <- create_cctable(rec, freq=1)
cctable <- cctable$clean()
#table <- cctable$tclean 
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{is.laboratory}
\alias{is.laboratory}
\title{Check if the item NHIC code or short name belongs to the Laboratory 
category.}
\usage{
is.laboratory(item_name)
}
\arguments{
\item{item_name}{character the NHIC code or the short name}
}
\value{
logical
}
\description{
Check if the item NHIC code or short name belongs to the Laboratory 
category.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{physio.distribution}
\alias{physio.distribution}
\title{Plot the physiological data distribution.}
\usage{
physio.distribution(cctb, names)
}
\arguments{
\item{cctb}{ccTable-class, see create.cctable().}

\item{names}{character vector of short names of numerical demographic data.}
}
\description{
Plot the physiological data distribution.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{ccTable_export_csv}
\alias{ccTable_export_csv}
\title{Export the clean table as a CSV file}
\arguments{
\item{file}{the full path of the output CSV file.}
}
\description{
Export tclean as a CSV file.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.missingness.R
\name{ccTable_filter_missingness}
\alias{ccTable_filter_missingness}
\title{Data missing filter}
\arguments{
\item{recount}{logical value. Recount the missingness if TRUE.}
}
\description{
Deal with data when insufficient data points are supported. There are 
two key items to be set in the YAML configuration file. 
1) labels -- time interval. 2) accept_2d -- the accept present ratio. 
So if we set the labels is 24, and accept_2d is 70. It means we accept 
all the missing rate that is lower than 30% every 24 data points.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{methods}
\name{+,ccRecord,NULL-method}
\alias{+,ccRecord,NULL-method}
\title{Adding nothing to a ccRecord object and return the original ccRecord}
\usage{
\S4method{+}{ccRecord,`NULL`}(e1, e2)
}
\arguments{
\item{e1}{ccRecord-class}

\item{e2}{NULL}
}
\description{
Adding nothing to a ccRecord object and return the original ccRecord
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{data.quality.report.brc}
\alias{data.quality.report.brc}
\title{Create the data quality report}
\usage{
data.quality.report.brc(ccd, pdf = TRUE, brc = NULL, path = NULL)
}
\arguments{
\item{ccd}{ccRecord}

\item{pdf}{logical create the pdf version of the DQ report, 
otherwise stay in markdown format}

\item{brc}{character BRC names which can be Cambridge, GSTT, Imperial,
Oxford, and UCLH.}

\item{path}{report export path}
}
\description{
Create a detailed data quality report, including file summary, site 
summary, data completeness, and density plot. The result can be found 
in {path}/report/data_quality_report.{pdf}/{md}. Using this function, 
one can also create a site/trust specified report, see the argument "site". 
You need to make sure that you have the right to write into the {work_dir}.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{xml.site.duration.plot}
\alias{xml.site.duration.plot}
\title{Plot the XML duration in terms of sites.}
\usage{
xml.site.duration.plot(ccd)
}
\arguments{
\item{ccd}{ccRecord-class}
}
\description{
Plot the XML duration in terms of sites.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{whichIsCode}
\alias{whichIsCode}
\title{give id number from NHIC code like "NIHR_HIC_ICU_xxxx"}
\usage{
whichIsCode(nhic)
}
\arguments{
\item{nhic}{NHIC code}
}
\description{
give id number from NHIC code like "NIHR_HIC_ICU_xxxx"
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccTable.R
\name{ccTable_reload_conf}
\alias{ccTable_reload_conf}
\title{Reload the YAML configuration file}
\arguments{
\item{conf}{full path of the YAML configuration file or the parsed config list.}
}
\description{
Note, this function will also reset all the operations and 
remove the tclean.
}
\examples{
\dontrun{
tb$reload_conf("REF.yaml")
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cchic_xml.R
\name{getXmlepisode}
\alias{getXmlepisode}
\title{get the episode data from xml}
\usage{
getXmlepisode(xml.root, id)
}
\arguments{
\item{xml.root}{root of xml data returned by xmlLoad()}

\item{id}{integer}
}
\description{
get the episode data from xml
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{methods}
\name{[,ccRecord,ANY-method}
\alias{[,ccRecord,ANY-method}
\title{Create a subset of ccRecord object from the original one via specifying the row number of episodes.}
\usage{
\S4method{[}{ccRecord,ANY}(x, i)
}
\arguments{
\item{x}{ccRecord-class}

\item{i}{integer vector}
}
\description{
Create a subset of ccRecord object from the original one via specifying the row number of episodes.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{methods}
\name{[,ccRecord,character-method}
\alias{[,ccRecord,character-method}
\title{Create a ccRecord subsetting via selected sites.}
\usage{
\S4method{[}{ccRecord,character}(x, i)
}
\arguments{
\item{x}{ccRecord-class}

\item{i}{character vector which contains site_ids, e.g. c("Q70", "Q70W")}
}
\description{
Create a ccRecord subsetting via selected sites.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{lookup.items}
\alias{lookup.items}
\title{Lookup items information by keywords}
\usage{
lookup.items(keyword, style = "grid")
}
\arguments{
\item{keyword}{character e.g. "h_rate", "heart", "108".}

\item{style}{character, the style of the table output which can be "simple",
"rmarkdown", and "grid"}
}
\value{
character the short names of the selected items.
}
\description{
This function tries to match keywords in short names, long names and NHIC code. 
The matched items will be displayed.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdid.R
\name{StdId}
\alias{StdId}
\title{S4 class to hold standard IDs such as "NIHR_HIC_ICU_0001"}
\usage{
StdId(text)

StdId(text)
}
\arguments{
\item{text}{NIHC code which should be in a format like NIHR_HIC_ICU_xxxx}
}
\description{
S4 class to hold standard IDs such as "NIHR_HIC_ICU_0001"

constructor of StdId class
}
\section{Slots}{

\describe{
\item{\code{ids}}{single or multiple characters}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demographics.R
\name{ccd_unique_spell}
\alias{ccd_unique_spell}
\title{find the unique spell ID.}
\usage{
ccd_unique_spell(rec, duration = 2)
}
\arguments{
\item{rec}{ccRecord-class}

\item{duration}{integer hours}
}
\value{
data.table contains spell id.
}
\description{
find the unique spell ID.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{class}
\name{ccRecord-class}
\alias{ccRecord}
\alias{ccRecord-class}
\title{The S4 class which holds all the CCHIC patient record - served as a database.}
\description{
ccRecord is a class to hold the raw episode data parsed directly 
from XML or CSV files.
}
\section{Fields}{

\describe{
\item{\code{nepisodes}}{is an integer number indicates the total number of episode
the record is holding.}

\item{\code{dmgtb}}{a data.table containing all the demographic information of each
episode, including site_id, NHS number, PAS number, admission date/time,
and discharge date/time. This field is usually left empty.}

\item{\code{infotb}}{a data.table holding the parsing information of each episode such as the
parsing time and from which file it parsed from.}

\item{\code{episdoes}}{a list of ccEpisode objects.}
}}
\examples{
heart_rate <- data.frame(seq(10), rep(70, 10)) # NIHR_HIC_ICU_0108
site_id <- "Q70" #  NIHR_HIC_ICU_0002
episode_id <- "0000001" # NIHR_HIC_ICU_0005

# Create a new episode 
ep <- new.episode(list(NIHR_HIC_ICU_0108=heart_rate, 
                         NIHR_HIC_ICU_0002=site_id, 
                         NIHR_HIC_ICU_0005=episode_id)) 

# modifying records 
rec <- ccRecord() # a new record 
rec <- rec + ep # adding a new episode to the record
rec <- rec + NULL # adding nothing to the record
rec <- rec + rec # adding a record to a record
# Adding a list of episodes 
rec <- ccRecord()
ep1 <- new.episode()
ep2 <- new.episode()
eps.list <- list(ep1, ep2)
new.rec <- rec + eps.list
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demographics.R
\name{lenstay}
\alias{lenstay}
\title{Calculate the length of stay in the ICU.}
\usage{
lenstay(demg, units = "hours")
}
\arguments{
\item{demg}{data.table the demographic table which should at least contain
column DAICU and DDICU}

\item{units}{character The unit of lenstay column, by default the output will be in hours}
}
\value{
data.table It is the original data.table with lenstay column (in 
difftime) appended.
}
\description{
Calculate the length of stay in the ICU and append it to the original demographic
table.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.range.R
\name{inrange}
\alias{inrange}
\title{Check if the values of a vector v is in the given ranges.}
\usage{
inrange(v, range)
}
\arguments{
\item{v}{vector numeric}

\item{range}{A string contains the numeric ranges in a form such as (low,
up) for open range and [low, up] for close range. Multiple
ranges should be separated by semi-columns which is equivalent to logical
OR e.g. (low1, up1); (low2, up2)}
}
\description{
Check if the values of a vector v is in the given ranges.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{site.info}
\alias{site.info}
\title{Produce a site id reference table.}
\usage{
site.info()
}
\value{
data.frame
}
\description{
Produce a site id reference table.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demographics.R
\name{ccd_demographic_spell}
\alias{ccd_demographic_spell}
\title{Create demographic table with spell IDs}
\usage{
ccd_demographic_spell(rec, duration = 2)
}
\arguments{
\item{rec}{ccRecord}

\item{duration}{the maximum hours of transition period}
}
\value{
data.table demographic table with spell ID in column spell
}
\description{
same output like ccd_demographic_table but in 
addition with a spell ID.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{demg.distribution}
\alias{demg.distribution}
\title{demg.distribution
Create a plot of the distribution of numerical demographic data.}
\usage{
demg.distribution(demg, names)
}
\arguments{
\item{demg}{ccRecord or demographic table created by ccd_demographic_table()}

\item{names}{character vector of short names of numerical demographic data.}
}
\description{
demg.distribution
Create a plot of the distribution of numerical demographic data.
}
\examples{
\dontrun{tdemg.distribution(ccd, "HCM")}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\name{plot_episode}
\alias{plot_episode}
\title{Individual episode chart}
\usage{
plot_episode(r, v)
}
\arguments{
\item{r}{ccEpisode-class}

\item{v}{short name of longitudinal data. While v is not given, the chart 
will only display h_rate, spo2, bilirubin, platelets, pao2_fio2, gcs_total.}
}
\value{
a table of selected vars of an episode
}
\description{
Create an individual episode chart for its diagnosis, drugs and physiological
variables. Diagnosis and drugs are always included, while the user can
select other longitudinal data.
}
\examples{
\dontrun{
plot_episode(ccd@episodes[[1]]) # plot first episode with default variables. 
plot_episode(ccd@episodes[[1]], "h_rate") # plot first episode heart rate
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.quality.report.R
\name{data.quality.report}
\alias{data.quality.report}
\title{Create the data quality report}
\usage{
data.quality.report(ccd, site = NULL, file = NULL, pdf = TRUE,
  out = "report")
}
\arguments{
\item{ccd}{ccRecord}

\item{site}{a vector of the site ids for the site specified report.}

\item{file}{character a list of XML file origins.}

\item{pdf}{logical create the pdf version of the DQ report, 
otherwise stay in markdown format}

\item{out}{character output path}
}
\description{
Create a detailed data quality report, including file summary, site 
summary, data completeness, and density plot. The result can be found 
in {work_dir}/report/data_quality_report.{pdf}/{md}. Using this function, 
one can also create a site/trust specified report, see the argument "site". 
You need to make sure that you have the right to write into the {work_dir}.
}
\examples{
\dontrun{data.quality.report(ccd, c("Q70", "C90"))}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccRecord.R
\docType{methods}
\name{plot_episode,ccEpisode,character-method}
\alias{plot_episode,ccEpisode,character-method}
\title{Episode chart}
\usage{
\S4method{plot_episode}{ccEpisode,character}(r, v)
}
\arguments{
\item{r}{ccEpisode-class}

\item{v}{character}
}
\description{
Episode chart
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cchic_xml.R
\name{extract_file_origin}
\alias{extract_file_origin}
\title{Extract the original file name from a path and file removing
all the suffixes.}
\usage{
extract_file_origin(pathfile, removestr = ".xml")
}
\arguments{
\item{pathfile}{a particular file name which may have a suffix}

\item{removestr}{last bit from the original filename}
}
\value{
string
}
\description{
Extract the original file name from a path and file removing
all the suffixes.
}

