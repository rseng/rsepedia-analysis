kmdata
====

This R<sup>[[1]](#1)</sup> package contains a database of 304 reconstructed,
patient-level clinical trial datasets on multiple survival endpoints,
including overall and progression-free survival.

The data have been extracted from Kaplan-Meier (KM) curves reported in 153
oncology Phase III clinical trial publications.

Studies have been identified through a [PubMed](https://pubmed.ncbi.nlm.nih.gov/)
search of clinical trials in breast, lung cancer, prostate, and colorectal
cancer published between 2014 and 2016.

For each trial that met the our search criteria, we extracted and curated
study-level information.

All reported KM survival curves were digitized with the software
[DigitizeIt](https://www.digitizeit.de).

The digitized KM survival curves were used to estimate (possibly censored)
patient-level event times using the Guyot-algorithm<sup>[[2]](#2)</sup>.

### References
<a id="1">[1]</a>
R-Core-Team.
A Language and Environment for Statistical Computing.
R Found Stat Comput.
2018;2: https://www.r-project.org.

<a id="2">[2]</a>
Guyot P, Ades AE, Ouwens MJNM, Welton NJ.
Enhanced secondary analysis of survival data: Reconstructing the data from
published Kaplan-Meier survival curves.
_BMC Med Res Methodol_.
2012;12. doi:10.1186/1471-2288-12-9

---

### Installation

```r
# install.packages('devtools')
devtools::install_github('raredd/kmdata', build_vignettes = TRUE)
```

### Getting started

View a list of the data available in the `kmdata` package:

```r
data(package = 'kmdata')
```

Or generate individual patient data (IPD) from digitized survival curve data
using the method described by Guyot:

```r
library('survival')
sf <- survfit(Surv(time, status) ~ 1, cancer)
aa <- summary(sf, times = 0:8 * 100)

pt <- ipd(sf$time, sf$surv, aa$time, aa$n.risk)
```

Compare:

```r
sf
````

```
# Call: survfit(formula = Surv(time, status) ~ 1, data = cancer)
# 
#       n  events  median 0.95LCL 0.95UCL 
#     228     165     310     285     363 
```


```r
survfit(Surv(time, event) ~ 1, pt)
```

```
# Call: survfit(formula = Surv(time, event) ~ 1, data = pt)
# 
#       n  events  median 0.95LCL 0.95UCL 
#     228     163     310     285     363 
```

### Additional resources

See the `kmdata` [intro vignette](vignettes/kmdata-intro.Rmd) or more
information about [generating IPD](vignettes/kmdata-ipd.Rmd) from digitized
survival curves.
---
title: "kmdata - intro"
output:
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{kmdata - intro}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r, include=FALSE}
library('knitr')

opts_chunk$set(
  cache = FALSE, fig.align = 'center', dev = 'png',
  fig.width=9, fig.height=7, echo = TRUE, message = FALSE
)
```

## Data sets available

```{r}
library('kmdata')

data(package = 'kmdata')
```

This will show a list of the data sets available. Data objects can be
references by `ATTENTION_2A` or `` `TRIBE(2)_2A` `` (with backticks) for data
names with special characters.

The name consists of the study short-name and figure identifier:

```{r}
data <- ls('package:kmdata', pattern = '^[A-Z]')
cbind(
  name = data,
  study = gsub('_.*', '', data),
  figure = gsub('.*_', '', data)
)[1:5, ]
```

For example, study "ATTENTION" has two figures: "2A" and "2B." Study "ACTSCC"
has only one figure: "2A." If data sets sourced from multi-panel figures, the
name will look similar to study "ATTENTION" with figure IDs "2A" and "2B."

All data sets are listed in `kmdata_key` along with some useful metadata for
each including the journal and publication identifiers, outcomes and study
arms, quality of the re-capitulated data, and other information.

## Working with the data

Each data set contains the same format for consistency:

```{r, echo=FALSE}
knitr::kable(
  data.frame(
    time = 'time-to-event (in units)', event = 'event indicator (0/1)',
    arm = 'treatment arm identifier (e.g., arm-1 vs arm-2)'
  )
)
```

The time unit, event type, and treatment arms can be found in the help page
for each data set, e.g., `?ACT1_2A`. Additionally, the data objects contain
metadata stored as attributes:

```{r}
head(ACT1_2A)

attr(ACT1_2A, 'event')

attributes(ACT1_2A)[-(1:3)]
```

Data may be examined and plotted using the built-in functions `summary` and
`kmplot`.

```{r, fig.show='hold', message = TRUE}
summary(ACT1_2A)

kmplot(ACT1_2A)
```

## Selecting data

The `kmdata` package contains a function, `select_kmdata`, to easily search
and filter data sets which share common features. Any of the columns in
`kmdata_key` may be used to filter.

For example, if we wanted a list of lung cancer data sets with overall
survival (OS) in months with fewer than 500 patients reporting at least a
1.2 hazard ratio for treatment compared to a reference arm, we can use the
following:

```{r}
select_kmdata(
  Cancer %in% 'Lung' &
    Outcome %in% 'OS' &
    Units %in% 'months' &
    ReportedSampleSize < 500 &
    HazardRatio >= 1.2,
  return = 'name'
)
```

By default, `select_kmdata` returns only the names of the data sets for
reference individually (i.e., `select_kmdata(..., return = 'name')`), but
it can also return the matching rows of `kmdata_key` or the matching data
sets as a list.

```{r, fig.height=11, fig.width=12, results='hide'}
key <- select_kmdata(
  Cancer %in% 'Lung' &
    Outcome %in% 'OS' &
    Units %in% 'months' &
    ReportedSampleSize < 500 &
    HazardRatio >= 1.2,
  return = 'key'
)

dat <- select_kmdata(
  Cancer %in% 'Lung' &
    Outcome %in% 'OS' &
    Units %in% 'months' &
    ReportedSampleSize < 500 &
    HazardRatio >= 1.2,
  return = 'data'
)

par(mfrow = n2mfrow(length(dat)))
for (dd in dat)
  kmplot(dd)
```

## Data quality

Each figure and data set contains a quality score which represents how well
the re-capitulated agrees with the original publication. Scores range from
0 (worst) to 100% (best) and are an aggregation of four metrics: hazard ratio,
total events, median time-to-event, and number at-risk.

Each metric is score from 0 (worst) to 3 (best); the maximum score per figure
may vary with the metrics reported in the original publication. For example,
if only one was reported, the maximum score is 3/3.

A score of 3 points is given per metric per figure if the re-capitulated
metric is no more than 5% different than the published, 2 points are given
if the metric is 5-10% different, 1 point for 10-20%, and 0 points for more
than 20% different.

|  % difference from publication |  Quality points per metric |
|-------------------------------:|---------------------------:|
|                            0-5 |                          3 |
|                           5-10 |                          2 |
|                          10-20 |                          1 |
|                        &gt; 20 |                          0 |

---

## References

The publications and figures available in this package are listed below by
first author.

<details>
<summary>Click to expand</summary>
```{r, echo=FALSE}
cit <- system.file('docs', 'Citations_final.xlsx', package = 'kmdata')
cit <- as.data.frame(readxl::read_excel(cit, skip = 1L))

cit <- within(cit, {
  Title    <- gsub('^.*?\\.\\s+|\\.\\s+[A-z ]+\\d{4};.*$', '', Reference)
  Author   <- gsub('^([^.]+\\.)|.', '\\1', Reference)
  PubData  <- gsub('([A-z ]+\\s+\\d{4};.*)\\.$|.', '\\1', Reference)
  Journal  <- gsub('(.*?)\\d{4};|.', '\\1', PubData)
  Year     <- gsub('(\\d{4});|.', '\\1', PubData)
  Location <- gsub('^.*?\\d{4};\\s*', '\\1', PubData)
})[, c('PMID', 'Author', 'Journal', 'Year', 'Title', 'Location')]
cit <- cit[order(cit$Author), ]
rownames(cit) <- NULL

knitr::kable(cit, format = 'markdown', caption = 'List of publications.')
```
</details>
---
title: "kmdata - ipd"
output:
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{kmdata - ipd}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r, include=FALSE}
library('knitr')

opts_chunk$set(
  cache = FALSE, fig.align = 'center', dev = 'png',
  fig.width=9, fig.height=7, echo = TRUE, message = FALSE
)
```

## Generating IPD

In addition to the recapitulated study data, the `kmdata` package contains
the algorithm to re-create patient-level time-to-event data from digitized
survival curves.

Digitized curves are simply the x- and y-coordinates of the Kaplan-Meier
curves which can be generated several ways including the \code{R} package
**`digitize`** or proprietary software such as
[DigitizeIt](https://digitizeit.soft112.com/).

```{r}
library('kmdata')

xy <- system.file(
  'etc', 'init', 'Checkmate_067_S3A_Nivolumab.csv', package = 'kmdata'
)
ar <- system.file(
  'etc', 'init', 'Checkmate_067_At_Risk.csv', package = 'kmdata'
)

dd <- read.csv(xy)
aa <- read.csv(ar)
aa <- aa[aa$nrisk > 0, ]

ipd <- ipd(dd$T, dd$S, aa$trisk, aa$nrisk, arm = 'Nivo')
kmplot(ipd, conf.int = FALSE, median = FALSE)

points(S ~ T, dd, type = 's', col = 'red')

legend(
  'topright', col = 1:2, lty = 1, bty = 'n', lwd = 2,
  legend = c('Re-capitulated', 'Truth (no censoring\ndata available)')
)
```

## References

Guyot, Patricia, et al. Enhanced secondary analysis of survival data:
reconstructing the data from published Kaplan-Meier survival curves.
_BMC Medical Research Methodology_ 2012, **12**:9.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MARQUEE_2A}
\alias{MARQUEE_2A}
\title{MARQUEE, figure 2A}
\format{
A data frame of 1,048 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (erlotinib_placebo, erlotinib_tivantinib) \cr
}
}
\source{
Scagliotti G, von Pawel J, Novello S, et al. Phase III Multinational,
Randomized, Double-Blind, Placebo-Controlled Study of Tivantinib (ARQ
197) Plus Erlotinib Versus Erlotinib Alone in Previously Treated
Patients With Locally Advanced or Metastatic Nonsquamous
Non-Small-Cell Lung Cancer. J Clin Oncol 2015; 33: 2667–74.
}
\usage{
MARQUEE_2A
}
\description{
Kaplan-Meier digitized data from MARQUEE, figure 2A (PMID 26169611). A reported sample size of 1,048 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(MARQUEE_2A)

kmplot(MARQUEE_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00673049_2A}
\alias{NCT00673049_2A}
\title{NCT00673049, figure 2A}
\format{
A data frame of 583 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, figitumumab) \cr
}
}
\source{
Scagliotti GV, Bondarenko I, Blackhall F, et al. Randomized, phase
III trial of figitumumab in combination with erlotinib versus
erlotinib alone in patients with nonadenocarcinoma nonsmall-cell lung
cancer. Ann Oncol 2015; 26: 497–504.
}
\usage{
NCT00673049_2A
}
\description{
Kaplan-Meier digitized data from NCT00673049, figure 2A (PMID 25395283). A reported sample size of 583 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(NCT00673049_2A)

kmplot(NCT00673049_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CALGB90202_2B}
\alias{CALGB90202_2B}
\title{CALGB90202, figure 2B}
\format{
A data frame of 645 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, za) \cr
}
}
\source{
Smith MR, Halabi S, Ryan CJ, et al. Randomized controlled trial of
early zoledronic acid in men with castration-sensitive prostate
cancer and bone metastases: results of CALGB 90202 (alliance). J Clin
Oncol 2014; 32: 1143–50.
}
\usage{
CALGB90202_2B
}
\description{
Kaplan-Meier digitized data from CALGB90202, figure 2B (PMID 24590644). A reported sample size of 645 for a primary endpoint of time_SRE in prostate cancer.
}
\examples{
summary(CALGB90202_2B)

kmplot(CALGB90202_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AIO0207(a)_3C}
\alias{AIO0207(a)_3C}
\title{AIO0207(a), figure 3C}
\format{
A data frame of 314 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (5fu_bev, no_tx) \cr
}
}
\source{
Hegewisch-Becker S, Graeven U, Lerchenmüller CA, et al. Maintenance
strategies after first-line oxaliplatin plus fluoropyrimidine plus
bevacizumab for patients with metastatic colorectal cancer (AIO
0207): a randomised, non-inferiority, open-label, phase 3 trial.
Lancet Oncol 2015; 16: 1355–69.
}
\usage{
`AIO0207(a)_3C`
}
\description{
Kaplan-Meier digitized data from AIO0207(a), figure 3C (PMID 26361971). A reported sample size of 472 for a primary endpoint of time_tx_failure in colorectal cancer.
}
\examples{
summary(`AIO0207(a)_3C`)

kmplot(`AIO0207(a)_3C`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ACT1_3A}
\alias{ACT1_3A}
\title{ACT1, figure 3A}
\format{
A data frame of 637 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (amrubicin, topotecan) \cr
}
}
\source{
von Pawel J, Jotte R, Spigel DR, et al. Randomized phase III trial of
amrubicin versus topotecan as second-line treatment for patients with
small-cell lung cancer. J Clin Oncol 2014; 32: 4012–9.
}
\usage{
ACT1_3A
}
\description{
Kaplan-Meier digitized data from ACT1, figure 3A (PMID 25385727). A reported sample size of 637 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ACT1_3A)

kmplot(ACT1_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BREC_2A}
\alias{BREC_2A}
\title{BREC, figure 2A}
\format{
A data frame of 277 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, experimental) \cr
}
}
\source{
Moran T, Wei J, Cobo M, et al. Two biomarker-directed randomized
trials in European and Chinese patients with nonsmall-cell lung
cancer: the BRCA1-RAP80 Expression Customization (BREC) studies. Ann
Oncol 2014; 25: 2147–55.
}
\usage{
BREC_2A
}
\description{
Kaplan-Meier digitized data from BREC, figure 2A (PMID 25164908). A reported sample size of 279 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(BREC_2A)

kmplot(BREC_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CAIRO3_2D}
\alias{CAIRO3_2D}
\title{CAIRO3, figure 2D}
\format{
A data frame of 557 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (maintenance, observation) \cr
}
}
\source{
Simkens LHJ, van Tinteren H, May A, et al. Maintenance treatment with
capecitabine and bevacizumab in metastatic colorectal cancer
(CAIRO3): a phase 3 randomised controlled trial of the Dutch
Colorectal Cancer Group. Lancet 2015; 385: 1843–52.
}
\usage{
CAIRO3_2D
}
\description{
Kaplan-Meier digitized data from CAIRO3, figure 2D (PMID 25862517). A reported sample size of 558 for a primary endpoint of PFS2 in colorectal cancer.
}
\examples{
summary(CAIRO3_2D)

kmplot(CAIRO3_2D)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{KEYNOTE024_1A}
\alias{KEYNOTE024_1A}
\title{KEYNOTE024, figure 1A}
\format{
A data frame of 305 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, pembrolizumab) \cr
}
}
\source{
Reck M, Rodríguez-Abreu D, Robinson AG, et al. Pembrolizumab versus
Chemotherapy for PD-L1-Positive Non-Small-Cell Lung Cancer. N Engl J
Med 2016; 375: 1823–33.
}
\usage{
KEYNOTE024_1A
}
\description{
Kaplan-Meier digitized data from KEYNOTE024, figure 1A (PMID 27718847). A reported sample size of 305 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(KEYNOTE024_1A)

kmplot(KEYNOTE024_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00938652_2A}
\alias{NCT00938652_2A}
\title{NCT00938652, figure 2A}
\format{
A data frame of 519 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (gc, gci) \cr
}
}
\source{
O’Shaughnessy J, Schwartzberg L, Danso MA, et al. Phase III study of
iniparib plus gemcitabine and carboplatin versus gemcitabine and
carboplatin in patients with metastatic triple-negative breast
cancer. J Clin Oncol 2014; 32: 3840–7.
}
\usage{
NCT00938652_2A
}
\description{
Kaplan-Meier digitized data from NCT00938652, figure 2A (PMID 25349301). A reported sample size of 519 for a primary endpoint of OS/PFS in breast cancer.
}
\examples{
summary(NCT00938652_2A)

kmplot(NCT00938652_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BEBYP_1}
\alias{BEBYP_1}
\title{BEBYP, figure 1}
\format{
A data frame of 184 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev_ct, ct) \cr
}
}
\source{
Masi G, Salvatore L, Boni L, et al. Continuation or reintroduction of
bevacizumab beyond progression to first-line therapy in metastatic
colorectal cancer: final results of the randomized BEBYP trial. Ann
Oncol 2015; 26: 724–30.
}
\usage{
BEBYP_1
}
\description{
Kaplan-Meier digitized data from BEBYP, figure 1 (PMID 25600568). A reported sample size of 185 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(BEBYP_1)

kmplot(BEBYP_1)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CHHiP(b)_2A}
\alias{CHHiP(b)_2A}
\title{CHHiP(b), figure 2A}
\format{
A data frame of 2,142 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab RFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (57gy, 74gy) \cr
}
}
\source{
Dearnaley D, Syndikus I, Mossop H, et al. Conventional versus
hypofractionated high-dose intensity-modulated radiotherapy for
prostate cancer: 5-year outcomes of the randomised, non-inferiority,
phase 3 CHHiP trial. Lancet Oncol 2016; 17: 1047–60.
}
\usage{
`CHHiP(b)_2A`
}
\description{
Kaplan-Meier digitized data from CHHiP(b), figure 2A (PMID 27339115). A reported sample size of 3,216 for a primary endpoint of bRFS in prostate cancer.
}
\examples{
summary(`CHHiP(b)_2A`)

kmplot(`CHHiP(b)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{INSPIRE_2A}
\alias{INSPIRE_2A}
\title{INSPIRE, figure 2A}
\format{
A data frame of 633 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ct, necitumumab_ct) \cr
}
}
\source{
Paz-Ares L, Mezger J, Ciuleanu TE, et al. Necitumumab plus pemetrexed
and cisplatin as first-line therapy in patients with stage IV
non-squamous non-small-cell lung cancer (INSPIRE): an open-label,
randomised, controlled phase 3 study. Lancet Oncol 2015; 16: 328–37.
}
\usage{
INSPIRE_2A
}
\description{
Kaplan-Meier digitized data from INSPIRE, figure 2A (PMID 25701171). A reported sample size of 633 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(INSPIRE_2A)

kmplot(INSPIRE_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{FIRE3_2A}
\alias{FIRE3_2A}
\title{FIRE3, figure 2A}
\format{
A data frame of 592 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfiri_bevacizumab, folfiri_cetuximab) \cr
}
}
\source{
Heinemann V, von Weikersthal LF, Decker T, et al. FOLFIRI plus
cetuximab versus FOLFIRI plus bevacizumab as first-line treatment for
patients with metastatic colorectal cancer (FIRE-3): a randomised,
open-label, phase 3 trial. Lancet Oncol 2014; 15: 1065–75.
}
\usage{
FIRE3_2A
}
\description{
Kaplan-Meier digitized data from FIRE3, figure 2A (PMID 25088940). A reported sample size of 592 for a primary endpoint of objective_response in colorectal cancer.
}
\examples{
summary(FIRE3_2A)

kmplot(FIRE3_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PROFILE1014_2A}
\alias{PROFILE1014_2A}
\title{PROFILE1014, figure 2A}
\format{
A data frame of 343 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, crizotinib) \cr
}
}
\source{
Solomon BJ, Cappuzzo F, Felip E, et al. Intracranial Efficacy of
Crizotinib Versus Chemotherapy in Patients With Advanced ALK-Positive
Non–Small-Cell Lung Cancer: Results From PROFILE 1014. J Clin Oncol
2016; 34: 2858-67.
}
\usage{
PROFILE1014_2A
}
\description{
Kaplan-Meier digitized data from PROFILE1014, figure 2A (PMID 27022118). A reported sample size of 343 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(PROFILE1014_2A)

kmplot(PROFILE1014_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ERCC1_2D}
\alias{ERCC1_2D}
\title{ERCC1, figure 2D}
\format{
A data frame of 173 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cisplatin, p1) \cr
}
}
\source{
Lee SM, Falzon M, Blackhall F, et al. Randomized Prospective
Biomarker Trial of ERCC1 for Comparing Platinum and Nonplatinum
Therapy in Advanced Non-Small-Cell Lung Cancer: ERCC1 Trial (ET). J
Clin Oncol 2017; 35: 402–11.
}
\usage{
ERCC1_2D
}
\description{
Kaplan-Meier digitized data from ERCC1, figure 2D (PMID 27893326). A reported sample size of 648 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ERCC1_2D)

kmplot(ERCC1_2D)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NSABPR04_3B}
\alias{NSABPR04_3B}
\title{NSABPR04, figure 3B}
\format{
A data frame of 1,284 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab LRR event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (no_oxali, oxali) \cr
}
}
\source{
Allegra CJ, Yothers G, O’Connell MJ, et al. Neoadjuvant 5-FU or
Capecitabine Plus Radiation With or Without Oxaliplatin in Rectal
Cancer Patients: A Phase III Randomized Clinical Trial. J Natl Cancer
Inst 2015; 107. DOI:10.1093/jnci/djv248.
}
\usage{
NSABPR04_3B
}
\description{
Kaplan-Meier digitized data from NSABPR04, figure 3B (PMID 26374429). A reported sample size of 1,608 for a primary endpoint of Local-regional tumor control in colorectal cancer.
}
\examples{
summary(NSABPR04_3B)

kmplot(NSABPR04_3B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUMELung1_2A}
\alias{LUMELung1_2A}
\title{LUMELung1, figure 2A}
\format{
A data frame of 1,134 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (docetaxel_nintedanib, docetaxel_placebo) \cr
}
}
\source{
Reck M, Kaiser R, Mellemgaard A, et al. Docetaxel plus nintedanib
versus docetaxel plus placebo in patients with previously treated
non-small-cell lung cancer (LUME-Lung 1): a phase 3, double-blind,
randomised controlled trial. Lancet Oncol 2014; 15: 143–55.
}
\usage{
LUMELung1_2A
}
\description{
Kaplan-Meier digitized data from LUMELung1, figure 2A (PMID 24411639). A reported sample size of 655 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(LUMELung1_2A)

kmplot(LUMELung1_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{E1199(c)_2C}
\alias{E1199(c)_2C}
\title{E1199(c), figure 2C}
\format{
A data frame of 2,483 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d1, p3) \cr
}
}
\source{
Sparano JA, Zhao F, Martino S, et al. Long-Term Follow-Up of the
E1199 Phase III Trial Evaluating the Role of Taxane and Schedule in
Operable Breast Cancer. J Clin Oncol 2015; 33: 2353–60.
}
\usage{
`E1199(c)_2C`
}
\description{
Kaplan-Meier digitized data from E1199(c), figure 2C (PMID 26077235). A reported sample size of 4,954 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`E1199(c)_2C`)

kmplot(`E1199(c)_2C`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00337103_2A}
\alias{NCT00337103_2A}
\title{NCT00337103, figure 2A}
\format{
A data frame of 1,102 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cap, eribulin) \cr
}
}
\source{
Kaufman PA, Awada A, Twelves C, et al. Phase III open-label
randomized study of eribulin mesylate versus capecitabine in patients
with locally advanced or metastatic breast cancer previously treated
with an anthracycline and a taxane. J Clin Oncol 2015; 33: 594–601.
}
\usage{
NCT00337103_2A
}
\description{
Kaplan-Meier digitized data from NCT00337103, figure 2A (PMID 25605862). A reported sample size of 1,102 for a primary endpoint of OS/PFS in breast cancer.
}
\examples{
summary(NCT00337103_2A)

kmplot(NCT00337103_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NTR1527_4}
\alias{NTR1527_4}
\title{NTR1527, figure 4}
\format{
A data frame of 495 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, thoracic rt) \cr
}
}
\source{
Slotman BJ, van Tinteren H, Praag JO, et al. Use of thoracic
radiotherapy for extensive stage small-cell lung cancer: a phase 3
randomised controlled trial. Lancet 2015; 385: 36–42.
}
\usage{
NTR1527_4
}
\description{
Kaplan-Meier digitized data from NTR1527, figure 4 (PMID 25230595). A reported sample size of 498 for a primary endpoint of OS-1yr in lung cancer.
}
\examples{
summary(NTR1527_4)

kmplot(NTR1527_4)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SELECTBC_2A}
\alias{SELECTBC_2A}
\title{SELECTBC, figure 2A}
\format{
A data frame of 592 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (s_1, taxane) \cr
}
}
\source{
Takashima T, Mukai H, Hara F, et al. Taxanes versus S-1 as the
first-line chemotherapy for metastatic breast cancer (SELECT BC): an
open-label, non-inferiority, randomised phase 3 trial. Lancet Oncol
2016; 17: 90–8.
}
\usage{
SELECTBC_2A
}
\description{
Kaplan-Meier digitized data from SELECTBC, figure 2A (PMID 26617202). A reported sample size of 618 for a primary endpoint of OS in breast cancer.
}
\examples{
summary(SELECTBC_2A)

kmplot(SELECTBC_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{COUAA202_2}
\alias{COUAA202_2}
\title{COUAA202, figure 2}
\format{
A data frame of 1,088 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (abiraterone, placebo) \cr
}
}
\source{
Ryan CJ, Smith MR, Fizazi K, et al. Abiraterone acetate plus
prednisone versus placebo plus prednisone in chemotherapy-naive men
with metastatic castration-resistant prostate cancer (COU-AA-302):
final overall survival analysis of a randomised, double-blind,
placebo-controlled phase 3 study. Lancet Oncol 2015; 16: 152–60.
}
\usage{
COUAA202_2
}
\description{
Kaplan-Meier digitized data from COUAA202, figure 2 (PMID 25601341). A reported sample size of 1,088 for a primary endpoint of rPFS/OS in prostate cancer.
}
\examples{
summary(COUAA202_2)

kmplot(COUAA202_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NSABPB35_2}
\alias{NSABPB35_2}
\title{NSABPB35, figure 2}
\format{
A data frame of 3,077 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab EFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (anastrazole, tamoxifen) \cr
}
}
\source{
Margolese RG, Cecchini RS, Julian TB, et al. Anastrozole versus
tamoxifen in postmenopausal women with ductal carcinoma in situ
undergoing lumpectomy plus radiotherapy (NSABP B-35): a randomised,
double-blind, phase 3 clinical trial. Lancet 2016; 387: 849–56.
}
\usage{
NSABPB35_2
}
\description{
Kaplan-Meier digitized data from NSABPB35, figure 2 (PMID 26686957). A reported sample size of 3,104 for a primary endpoint of BCFS in breast cancer.
}
\examples{
summary(NSABPB35_2)

kmplot(NSABPB35_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NeoALTTO(a)_2A}
\alias{NeoALTTO(a)_2A}
\title{NeoALTTO(a), figure 2A}
\format{
A data frame of 303 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab EFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l, t) \cr
}
}
\source{
de Azambuja E, Holmes AP, Piccart-Gebhart M, et al. Lapatinib with
trastuzumab for HER2-positive early breast cancer (NeoALTTO):
survival outcomes of a randomised, open-label, multicentre, phase 3
trial and their association with pathological complete response.
Lancet Oncol 2014; 15: 1137–46.
}
\usage{
`NeoALTTO(a)_2A`
}
\description{
Kaplan-Meier digitized data from NeoALTTO(a), figure 2A (PMID 25130998). A reported sample size of 455 for a primary endpoint of pCR in breast cancer.
}
\examples{
summary(`NeoALTTO(a)_2A`)

kmplot(`NeoALTTO(a)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BEACON_2A}
\alias{BEACON_2A}
\title{BEACON, figure 2A}
\format{
A data frame of 852 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, etirinotecan_pegol) \cr
}
}
\source{
Perez EA, Awada A, O’Shaughnessy J, et al. Etirinotecan pegol
(NKTR-102) versus treatment of physician’s choice in women with
advanced breast cancer previously treated with an anthracycline, a
taxane, and capecitabine (BEACON): a randomised, open-label,
multicentre, phase 3 trial. Lancet Oncol 2015; 16: 1556–68.
}
\usage{
BEACON_2A
}
\description{
Kaplan-Meier digitized data from BEACON, figure 2A (PMID 26482278). A reported sample size of 852 for a primary endpoint of OS in breast cancer.
}
\examples{
summary(BEACON_2A)

kmplot(BEACON_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CONCUR_2A}
\alias{CONCUR_2A}
\title{CONCUR, figure 2A}
\format{
A data frame of 204 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, regorafenib) \cr
}
}
\source{
Li J, Qin S, Xu R, et al. Regorafenib plus best supportive care
versus placebo plus best supportive care in Asian patients with
previously treated metastatic colorectal cancer (CONCUR): a
randomised, double-blind, placebo-controlled, phase 3 trial. Lancet
Oncol 2015; 16: 619–29.
}
\usage{
CONCUR_2A
}
\description{
Kaplan-Meier digitized data from CONCUR, figure 2A (PMID 25981818). A reported sample size of 243 for a primary endpoint of OS in colorectal cancer.
}
\examples{
summary(CONCUR_2A)

kmplot(CONCUR_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AURA3_1A}
\alias{AURA3_1A}
\title{AURA3, figure 1A}
\format{
A data frame of 419 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (osimertinib, platinum-pemetrexed) \cr
}
}
\source{
Mok TS, Wu Y-L, Ahn M-J, et al. Osimertinib or Platinum-Pemetrexed in
EGFR T790M-Positive Lung Cancer. N Engl J Med 2017; 376: 629–40.
}
\usage{
AURA3_1A
}
\description{
Kaplan-Meier digitized data from AURA3, figure 1A (PMID 27959700). A reported sample size of 419 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(AURA3_1A)

kmplot(AURA3_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{INT0142_2A}
\alias{INT0142_2A}
\title{INT0142, figure 2A}
\format{
A data frame of 337 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (tamoxifen, tamoxifen_ofs) \cr
}
}
\source{
Tevaarwerk AJ, Wang M, Zhao F, et al. Phase III comparison of
tamoxifen versus tamoxifen plus ovarian function suppression in
premenopausal women with node-negative, hormone receptor-positive
breast cancer (E-3193, INT-0142): a trial of the Eastern Cooperative
Oncology Group. J Clin Oncol 2014; 32: 3948–58.
}
\usage{
INT0142_2A
}
\description{
Kaplan-Meier digitized data from INT0142, figure 2A (PMID 25349302). A reported sample size of 345 for a primary endpoint of DFS/OS in breast cancer.
}
\examples{
summary(INT0142_2A)

kmplot(INT0142_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ExteNET_2A}
\alias{ExteNET_2A}
\title{ExteNET, figure 2A}
\format{
A data frame of 2,840 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (neratinib, placebo) \cr
}
}
\source{
Chan A, Delaloge S, Holmes FA, et al. Neratinib after
trastuzumab-based adjuvant therapy in patients with HER2-positive
breast cancer (ExteNET): a multicentre, randomised, double-blind,
placebo-controlled, phase 3 trial. Lancet Oncol 2016; 17: 367–77.
}
\usage{
ExteNET_2A
}
\description{
Kaplan-Meier digitized data from ExteNET, figure 2A (PMID 26874901). A reported sample size of 2,840 for a primary endpoint of iDFS in breast cancer.
}
\examples{
summary(ExteNET_2A)

kmplot(ExteNET_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{EAGLE_1}
\alias{EAGLE_1}
\title{EAGLE, figure 1}
\format{
A data frame of 368 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev10, bev5) \cr
}
}
\source{
Iwamoto S, Takahashi T, Tamagawa H, et al. FOLFIRI plus bevacizumab
as second-line therapy in patients with metastatic colorectal cancer
after first-line bevacizumab plus oxaliplatin-based therapy: the
randomized phase III EAGLE study. Ann Oncol 2015; 26: 1427–33.
}
\usage{
EAGLE_1
}
\description{
Kaplan-Meier digitized data from EAGLE, figure 1 (PMID 25908603). A reported sample size of 369 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(EAGLE_1)

kmplot(EAGLE_1)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CAIRO3_2A}
\alias{CAIRO3_2A}
\title{CAIRO3, figure 2A}
\format{
A data frame of 557 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (maintenance, observation) \cr
}
}
\source{
Simkens LHJ, van Tinteren H, May A, et al. Maintenance treatment with
capecitabine and bevacizumab in metastatic colorectal cancer
(CAIRO3): a phase 3 randomised controlled trial of the Dutch
Colorectal Cancer Group. Lancet 2015; 385: 1843–52.
}
\usage{
CAIRO3_2A
}
\description{
Kaplan-Meier digitized data from CAIRO3, figure 2A (PMID 25862517). A reported sample size of 558 for a primary endpoint of PFS2 in colorectal cancer.
}
\examples{
summary(CAIRO3_2A)

kmplot(CAIRO3_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MA31_3}
\alias{MA31_3}
\title{MA31, figure 3}
\format{
A data frame of 652 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ltax_l, ttax_t) \cr
}
}
\source{
Gelmon KA, Boyle FM, Kaufman B, et al. Lapatinib or Trastuzumab Plus
Taxane Therapy for Human Epidermal Growth Factor Receptor 2-Positive
Advanced Breast Cancer: Final Results of NCIC CTG MA.31. J Clin Oncol
2015; 33: 1574–83.
}
\usage{
MA31_3
}
\description{
Kaplan-Meier digitized data from MA31, figure 3 (PMID 25779558). A reported sample size of 652 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(MA31_3)

kmplot(MA31_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TRIBE(2)_2B}
\alias{TRIBE(2)_2B}
\title{TRIBE(2), figure 2B}
\format{
A data frame of 508 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfiri_bevacizumab, folfoxiri_bevacizumab) \cr
}
}
\source{
Loupakis F, Cremolini C, Masi G, et al. Initial Therapy with
FOLFOXIRI and Bevacizumab for Metastatic Colorectal Cancer. N Engl J
Med 2014; 371: 1609-18.
}
\usage{
`TRIBE(2)_2B`
}
\description{
Kaplan-Meier digitized data from TRIBE(2), figure 2B (PMID 25337750). A reported sample size of 508 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(`TRIBE(2)_2B`)

kmplot(`TRIBE(2)_2B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CALGB40502_2B}
\alias{CALGB40502_2B}
\title{CALGB40502, figure 2B}
\format{
A data frame of 542 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (nab_paclitaxel, p1) \cr
}
}
\source{
Rugo HS, Barry WT, Moreno-Aspitia A, et al. Randomized Phase III
Trial of Paclitaxel Once Per Week Compared With Nanoparticle
Albumin-Bound Nab-Paclitaxel Once Per Week or Ixabepilone With
Bevacizumab As First-Line Chemotherapy for Locally Recurrent or
Metastatic Breast Cancer: CALGB 40502/NCCTG N063H (Alliance). J Clin
Oncol 2015; 33: 2361–9.
}
\usage{
CALGB40502_2B
}
\description{
Kaplan-Meier digitized data from CALGB40502, figure 2B (PMID 26056183). A reported sample size of 799 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CALGB40502_2B)

kmplot(CALGB40502_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GIM2_2B}
\alias{GIM2_2B}
\title{GIM2, figure 2B}
\format{
A data frame of 2,091 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ec_p, fec_p) \cr
}
}
\source{
Del Mastro L, De Placido S, Bruzzi P, et al. Fluorouracil and
dose-dense chemotherapy in adjuvant treatment of patients with
early-stage breast cancer: an open-label, 2 × 2 factorial, randomised
phase 3 trial. Lancet 2015; 385: 1863–72.
}
\usage{
GIM2_2B
}
\description{
Kaplan-Meier digitized data from GIM2, figure 2B (PMID 25740286). A reported sample size of 2,091 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(GIM2_2B)

kmplot(GIM2_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{METlung_2A}
\alias{METlung_2A}
\title{METlung, figure 2A}
\format{
A data frame of 499 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (onaturzumab_erlotinib, placebo_erlotinib) \cr
}
}
\source{
Spigel DR, Edelman MJ, O’Byrne K, et al. Results From the Phase III
Randomized Trial of Onartuzumab Plus Erlotinib Versus Erlotinib in
Previously Treated Stage IIIB or IV Non-Small-Cell Lung Cancer:
METLung. J Clin Oncol 2017; 35: 412–20.
}
\usage{
METlung_2A
}
\description{
Kaplan-Meier digitized data from METlung, figure 2A (PMID 27937096). A reported sample size of 499 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(METlung_2A)

kmplot(METlung_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TRIBE_2}
\alias{TRIBE_2}
\title{TRIBE, figure 2}
\format{
A data frame of 508 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfiri_bev, folfoxiri_bev) \cr
}
}
\source{
Cremolini C, Loupakis F, Antoniotti C, et al. FOLFOXIRI plus
bevacizumab versus FOLFIRI plus bevacizumab as first-line treatment
of patients with metastatic colorectal cancer: updated overall
survival and molecular subgroup analyses of the open-label, phase 3
TRIBE study. Lancet Oncol 2015; 16: 1306–15.
}
\usage{
TRIBE_2
}
\description{
Kaplan-Meier digitized data from TRIBE, figure 2 (PMID 26338525). A reported sample size of 508 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(TRIBE_2)

kmplot(TRIBE_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{OMEGA_3A}
\alias{OMEGA_3A}
\title{OMEGA, figure 3A}
\format{
A data frame of 78 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cap, peg_doxorubicin) \cr
}
}
\source{
Smorenburg CH, de Groot SM, van Leeuwen-Stok AE, et al. A randomized
phase III study comparing pegylated liposomal doxorubicin with
capecitabine as first-line chemotherapy in elderly patients with
metastatic breast cancer: results of the OMEGA study of the Dutch
Breast Cancer Research Group BOOG. Ann Oncol 2014; 25: 599–605.
}
\usage{
OMEGA_3A
}
\description{
Kaplan-Meier digitized data from OMEGA, figure 3A (PMID 24504445). A reported sample size of 78 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(OMEGA_3A)

kmplot(OMEGA_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AZURE_4}
\alias{AZURE_4}
\title{AZURE, figure 4}
\format{
A data frame of 3,359 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, zoledronic_acid) \cr
}
}
\source{
Coleman R, Cameron D, Dodwell D, et al. Adjuvant zoledronic acid in
patients with early breast cancer: final efficacy analysis of the
AZURE (BIG 01/04) randomised open-label phase 3 trial. Lancet Oncol
2014; 15: 997–1006.
}
\usage{
AZURE_4
}
\description{
Kaplan-Meier digitized data from AZURE, figure 4 (PMID 25035292). A reported sample size of 3,360 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(AZURE_4)

kmplot(AZURE_4)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TRAPEZE_1A}
\alias{TRAPEZE_1A}
\title{TRAPEZE, figure 1A}
\format{
A data frame of 757 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (no_za, za) \cr
}
}
\source{
James ND, Pirrie SJ, Pope AM, et al. Clinical Outcomes and Survival
Following Treatment of Metastatic Castrate-Refractory Prostate Cancer
With Docetaxel Alone or With Strontium-89, Zoledronic Acid, or Both:
The TRAPEZE Randomized Clinical Trial. JAMA Oncol 2016; 2: 493–9.
}
\usage{
TRAPEZE_1A
}
\description{
Kaplan-Meier digitized data from TRAPEZE, figure 1A (PMID 26794729). A reported sample size of 757 for a primary endpoint of cPFS,cost-effectiveness in prostate cancer.
}
\examples{
summary(TRAPEZE_1A)

kmplot(TRAPEZE_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MAGRIT_2A}
\alias{MAGRIT_2A}
\title{MAGRIT, figure 2A}
\format{
A data frame of 2,272 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (mage_a3, placebo) \cr
}
}
\source{
Vansteenkiste JF, Cho BC, Vanakesa T, et al. Efficacy of the MAGE-A3
cancer immunotherapeutic as adjuvant therapy in patients with
resected MAGE-A3-positive non-small-cell lung cancer (MAGRIT): a
randomised, double-blind, placebo-controlled, phase 3 trial. Lancet
Oncol 2016; 17: 822–35.
}
\usage{
MAGRIT_2A
}
\description{
Kaplan-Meier digitized data from MAGRIT, figure 2A (PMID 27132212). A reported sample size of 2,312 for a primary endpoint of DFS(3 primary in diff pop) in lung cancer.
}
\examples{
summary(MAGRIT_2A)

kmplot(MAGRIT_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{DART0105_2A}
\alias{DART0105_2A}
\title{DART0105, figure 2A}
\format{
A data frame of 343 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ltad, stad) \cr
}
}
\source{
Zapatero A, Guerrero A, Maldonado X, et al. High-dose radiotherapy
with short-term or long-term androgen deprivation in localised
prostate cancer (DART01/05 GICOR): a randomised, controlled, phase 3
trial. Lancet Oncol 2015; 16: 320–7.
}
\usage{
DART0105_2A
}
\description{
Kaplan-Meier digitized data from DART0105, figure 2A (PMID 25702876). A reported sample size of 178 for a primary endpoint of bDFS in prostate cancer.
}
\examples{
summary(DART0105_2A)

kmplot(DART0105_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NITRO_2A}
\alias{NITRO_2A}
\title{NITRO, figure 2A}
\format{
A data frame of 372 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, nitro_chemo) \cr
}
}
\source{
Davidson A, Veillard A-S, Tognela A, et al. A phase III randomized
trial of adding topical nitroglycerin to first-line chemotherapy for
advanced nonsmall-cell lung cancer: the Australasian lung cancer
trials group NITRO trial. Ann Oncol 2015; 26: 2280–6.
}
\usage{
NITRO_2A
}
\description{
Kaplan-Meier digitized data from NITRO, figure 2A (PMID 26347110). A reported sample size of 372 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(NITRO_2A)

kmplot(NITRO_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{IMPRESS_2}
\alias{IMPRESS_2}
\title{IMPRESS, figure 2}
\format{
A data frame of 265 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (gefitinib, placebo) \cr
}
}
\source{
Soria J-C, Wu Y-L, Nakagawa K, et al. Gefitinib plus chemotherapy
versus placebo plus chemotherapy in EGFR-mutation-positive
non-small-cell lung cancer after progression on first-line gefitinib
(IMPRESS): a phase 3 randomised trial. Lancet Oncol 2015; 16: 990–8.
}
\usage{
IMPRESS_2
}
\description{
Kaplan-Meier digitized data from IMPRESS, figure 2 (PMID 26159065). A reported sample size of 265 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(IMPRESS_2)

kmplot(IMPRESS_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PROCTOR-SCRIPT_1A}
\alias{PROCTOR-SCRIPT_1A}
\title{PROCTOR-SCRIPT, figure 1A}
\format{
A data frame of 437 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, observation) \cr
}
}
\source{
Breugom AJ, van Gijn W, Muller EW, et al. Adjuvant chemotherapy for
rectal cancer patients treated with preoperative (chemo)radiotherapy
and total mesorectal excision: a Dutch Colorectal Cancer Group (DCCG)
randomized phase III trial. Ann Oncol 2015; 26: 696–701.
}
\usage{
`PROCTOR-SCRIPT_1A`
}
\description{
Kaplan-Meier digitized data from PROCTOR-SCRIPT, figure 1A (PMID 25480874). A reported sample size of 437 for a primary endpoint of OS in colorectal cancer.
}
\examples{
summary(`PROCTOR-SCRIPT_1A`)

kmplot(`PROCTOR-SCRIPT_1A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ACT2_2A}
\alias{ACT2_2A}
\title{ACT2, figure 2A}
\format{
A data frame of 71 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (wt_b, wt_be) \cr
}
}
\source{
Hagman H, Frödin J-E, Berglund Å, et al. A randomized study of
KRAS-guided maintenance therapy with bevacizumab, erlotinib or
metronomic capecitabine after first-line induction treatment of
metastatic colorectal cancer: the Nordic ACT2 trial. Ann Oncol 2016;
27: 140–7.
}
\usage{
ACT2_2A
}
\description{
Kaplan-Meier digitized data from ACT2, figure 2A (PMID 26483047). A reported sample size of 233 for a primary endpoint of PFS-3mo in colorectal cancer.
}
\examples{
summary(ACT2_2A)

kmplot(ACT2_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SIRFLOX_3}
\alias{SIRFLOX_3}
\title{SIRFLOX, figure 3}
\format{
A data frame of 530 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (mfolfox6, mfolfox6_rt) \cr
}
}
\source{
van Hazel GA, Heinemann V, Sharma NK, et al. SIRFLOX: Randomized
Phase III Trial Comparing First-Line mFOLFOX6 (Plus or Minus
Bevacizumab) Versus mFOLFOX6 (Plus or Minus Bevacizumab) Plus
Selective Internal Radiation Therapy in Patients With Metastatic
Colorectal Cancer. J Clin Oncol 2016; 34: 1723–31.
}
\usage{
SIRFLOX_3
}
\description{
Kaplan-Meier digitized data from SIRFLOX, figure 3 (PMID 26903575). A reported sample size of 530 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(SIRFLOX_3)

kmplot(SIRFLOX_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00294996_1A}
\alias{NCT00294996_1A}
\title{NCT00294996, figure 1A}
\format{
A data frame of 363 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (mpt, pt) \cr
}
}
\source{
Baselga J, Manikhas A, Cortés J, et al. Phase III trial of
nonpegylated liposomal doxorubicin in combination with trastuzumab
and paclitaxel in HER2-positive metastatic breast cancer. Ann Oncol
2014; 25: 592–8.
}
\usage{
NCT00294996_1A
}
\description{
Kaplan-Meier digitized data from NCT00294996, figure 1A (PMID 24401928). A reported sample size of 181 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(NCT00294996_1A)

kmplot(NCT00294996_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ENESTg1_2B}
\alias{ENESTg1_2B}
\title{ENESTg1, figure 2B}
\format{
A data frame of 644 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (imatinib, nilotinib) \cr
}
}
\source{
Blay J-Y, Shen L, Kang Y-K, et al. Nilotinib versus imatinib as
first-line therapy for patients with unresectable or metastatic
gastrointestinal stromal tumours (ENESTg1): a randomised phase 3
trial. Lancet Oncol 2015; 16: 550–60.
}
\usage{
ENESTg1_2B
}
\description{
Kaplan-Meier digitized data from ENESTg1, figure 2B (PMID 25882987). A reported sample size of 647 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(ENESTg1_2B)

kmplot(ENESTg1_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MA17R_1B}
\alias{MA17R_1B}
\title{MA17R, figure 1B}
\format{
A data frame of 1,918 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l, placebo) \cr
}
}
\source{
Goss PE, Ingle JN, Pritchard KI, et al. Extending Aromatase-Inhibitor
Adjuvant Therapy to 10 Years. N Engl J Med 2016; 375: 209–19.
}
\usage{
MA17R_1B
}
\description{
Kaplan-Meier digitized data from MA17R, figure 1B (PMID 27264120). A reported sample size of 1,918 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(MA17R_1B)

kmplot(MA17R_1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AMAROS_2A}
\alias{AMAROS_2A}
\title{AMAROS, figure 2A}
\format{
A data frame of 1,425 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (dissection, radiotherapy) \cr
}
}
\source{
Donker M, van Tienhoven G, Straver ME, et al. Radiotherapy or surgery
of the axilla after a positive sentinel node in breast cancer (EORTC
10981-22023 AMAROS): a randomised, multicentre, open-label, phase 3
non-inferiority trial. Lancet Oncol 2014; 15: 1303–10.
}
\usage{
AMAROS_2A
}
\description{
Kaplan-Meier digitized data from AMAROS, figure 2A (PMID 25439688). A reported sample size of 4,823 for a primary endpoint of AxRecurrence-5yr in breast cancer.
}
\examples{
summary(AMAROS_2A)

kmplot(AMAROS_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCICBR26_2B}
\alias{NCICBR26_2B}
\title{NCICBR26, figure 2B}
\format{
A data frame of 720 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (dacomitinib, placebo) \cr
}
}
\source{
Ellis PM, Shepherd FA, Millward M, et al. Dacomitinib compared with
placebo in pretreated patients with advanced or metastatic
non-small-cell lung cancer (NCIC CTG BR.26): a double-blind,
randomised, phase 3 trial. Lancet Oncol 2014; 15: 1379–88.
}
\usage{
NCICBR26_2B
}
\description{
Kaplan-Meier digitized data from NCICBR26, figure 2B (PMID 25439692). A reported sample size of 720 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(NCICBR26_2B)

kmplot(NCICBR26_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PALOMA2_1A}
\alias{PALOMA2_1A}
\title{PALOMA2, figure 1A}
\format{
A data frame of 666 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (palbociclib_letrozole, placebo_letrozole) \cr
}
}
\source{
Finn RS, Martin M, Rugo HS, et al. Palbociclib and Letrozole in
Advanced Breast Cancer. N Engl J Med 2016; 375: 1925–36.
}
\usage{
PALOMA2_1A
}
\description{
Kaplan-Meier digitized data from PALOMA2, figure 1A (PMID 27959613). A reported sample size of 666 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(PALOMA2_1A)

kmplot(PALOMA2_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{R98_4}
\alias{R98_4}
\title{R98, figure 4}
\format{
A data frame of 357 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (5fulv, 5fulv_cpt11) \cr
}
}
\source{
Delbaldo C, Ychou M, Zawadi A, et al. Postoperative irinotecan in
resected stage II-III rectal cancer: final analysis of the French R98
Intergroup trial†. Ann Oncol 2015; 26: 1208–15.
}
\usage{
R98_4
}
\description{
Kaplan-Meier digitized data from R98, figure 4 (PMID 25739671). A reported sample size of 357 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(R98_4)

kmplot(R98_4)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{JCOG0509_2B}
\alias{JCOG0509_2B}
\title{JCOG0509, figure 2B}
\format{
A data frame of 284 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (amrubicin_cisplatin, irinotecan_cisplatin) \cr
}
}
\source{
Satouchi M, Kotani Y, Shibata T, et al. Phase III study comparing
amrubicin plus cisplatin with irinotecan plus cisplatin in the
treatment of extensive-disease small-cell lung cancer: JCOG 0509. J
Clin Oncol 2014; 32: 1262–8.
}
\usage{
JCOG0509_2B
}
\description{
Kaplan-Meier digitized data from JCOG0509, figure 2B (PMID 24638015). A reported sample size of 284 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(JCOG0509_2B)

kmplot(JCOG0509_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CHHiP(b)_2B}
\alias{CHHiP(b)_2B}
\title{CHHiP(b), figure 2B}
\format{
A data frame of 2,142 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (57gy, 74gy) \cr
}
}
\source{
Dearnaley D, Syndikus I, Mossop H, et al. Conventional versus
hypofractionated high-dose intensity-modulated radiotherapy for
prostate cancer: 5-year outcomes of the randomised, non-inferiority,
phase 3 CHHiP trial. Lancet Oncol 2016; 17: 1047–60.
}
\usage{
`CHHiP(b)_2B`
}
\description{
Kaplan-Meier digitized data from CHHiP(b), figure 2B (PMID 27339115). A reported sample size of 3,216 for a primary endpoint of bRFS in prostate cancer.
}
\examples{
summary(`CHHiP(b)_2B`)

kmplot(`CHHiP(b)_2B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MARIANNE(a)_2A}
\alias{MARIANNE(a)_2A}
\title{MARIANNE(a), figure 2A}
\format{
A data frame of 732 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (t-dm1, trastuzumab_taxane) \cr
}
}
\source{
Perez EA, Barrios C, Eiermann W, et al. Trastuzumab Emtansine With or
Without Pertuzumab Versus Trastuzumab Plus Taxane for Human Epidermal
Growth Factor Receptor 2-Positive, Advanced Breast Cancer: Primary
Results From the Phase III MARIANNE Study. J Clin Oncol 2017; 35:
141–8.
}
\usage{
`MARIANNE(a)_2A`
}
\description{
Kaplan-Meier digitized data from MARIANNE(a), figure 2A (PMID 28056202). A reported sample size of 1,095 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(`MARIANNE(a)_2A`)

kmplot(`MARIANNE(a)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PROSE_2A}
\alias{PROSE_2A}
\title{PROSE, figure 2A}
\format{
A data frame of 263 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, erlotinib) \cr
}
}
\source{
Gregorc V, Novello S, Lazzari C, et al. Predictive value of a
proteomic signature in patients with non-small-cell lung cancer
treated with second-line erlotinib or chemotherapy (PROSE): a
biomarker-stratified, randomised phase 3 trial. Lancet Oncol 2014;
15: 713–21.
}
\usage{
PROSE_2A
}
\description{
Kaplan-Meier digitized data from PROSE, figure 2A (PMID 24831979). A reported sample size of 142 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(PROSE_2A)

kmplot(PROSE_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BEYOND_3A}
\alias{BEYOND_3A}
\title{BEYOND, figure 3A}
\format{
A data frame of 276 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev_cp, placebo_cp) \cr
}
}
\source{
Zhou C, Wu Y-L, Chen G, et al. BEYOND: A Randomized, Double-Blind,
Placebo-Controlled, Multicenter, Phase III Study of First-Line
Carboplatin/Paclitaxel Plus Bevacizumab or Placebo in Chinese
Patients With Advanced or Recurrent Nonsquamous Non-Small-Cell Lung
Cancer. J Clin Oncol 2015; 33: 2197–204.
}
\usage{
BEYOND_3A
}
\description{
Kaplan-Meier digitized data from BEYOND, figure 3A (PMID 26014294). A reported sample size of 276 for a primary endpoint of PFS in prostate cancer.
}
\examples{
summary(BEYOND_3A)

kmplot(BEYOND_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUXLung3-6_2A}
\alias{LUXLung3-6_2A}
\title{LUXLung3-6, figure 2A}
\format{
A data frame of 345 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (afatinib, pemetrexed-cisplatin) \cr
}
}
\source{
Yang JC-H, Wu Y-L, Schuler M, et al. Afatinib versus cisplatin-based
chemotherapy for EGFR mutation-positive lung adenocarcinoma (LUX-Lung
3 and LUX-Lung 6): analysis of overall survival data from two
randomised, phase 3 trials. Lancet Oncol 2015; 16: 141–51.
}
\usage{
`LUXLung3-6_2A`
}
\description{
Kaplan-Meier digitized data from LUXLung3-6, figure 2A (PMID 25589191). A reported sample size of 345,364 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(`LUXLung3-6_2A`)

kmplot(`LUXLung3-6_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PROCLAIM_2A}
\alias{PROCLAIM_2A}
\title{PROCLAIM, figure 2A}
\format{
A data frame of 598 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (eto_cis, pem_cis) \cr
}
}
\source{
Senan S, Brade A, Wang L-H, et al. PROCLAIM: Randomized Phase III
Trial of Pemetrexed-Cisplatin or Etoposide-Cisplatin Plus Thoracic
Radiation Therapy Followed by Consolidation Chemotherapy in Locally
Advanced Nonsquamous Non-Small-Cell Lung Cancer. J Clin Oncol 2016;
34: 953–62.
}
\usage{
PROCLAIM_2A
}
\description{
Kaplan-Meier digitized data from PROCLAIM, figure 2A (PMID 26811519). A reported sample size of 598 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(PROCLAIM_2A)

kmplot(PROCLAIM_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{KEYNOTE024_2}
\alias{KEYNOTE024_2}
\title{KEYNOTE024, figure 2}
\format{
A data frame of 305 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, pembrolizumab) \cr
}
}
\source{
Reck M, Rodríguez-Abreu D, Robinson AG, et al. Pembrolizumab versus
Chemotherapy for PD-L1-Positive Non-Small-Cell Lung Cancer. N Engl J
Med 2016; 375: 1823–33.
}
\usage{
KEYNOTE024_2
}
\description{
Kaplan-Meier digitized data from KEYNOTE024, figure 2 (PMID 27718847). A reported sample size of 305 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(KEYNOTE024_2)

kmplot(KEYNOTE024_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Chronicle_2B}
\alias{Chronicle_2B}
\title{Chronicle, figure 2B}
\format{
A data frame of 113 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (capecitabine_oxaliplatin, follow_up_only) \cr
}
}
\source{
Glynne-Jones R, Counsell N, Quirke P, et al. Chronicle: results of a
randomised phase III trial in locally advanced rectal cancer after
neoadjuvant chemoradiation randomising postoperative adjuvant
capecitabine plus oxaliplatin (XELOX) versus control. Ann Oncol 2014;
25: 1356–62.
}
\usage{
Chronicle_2B
}
\description{
Kaplan-Meier digitized data from Chronicle, figure 2B (PMID 24718885). A reported sample size of 113 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(Chronicle_2B)

kmplot(Chronicle_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUXLung5_1A}
\alias{LUXLung5_1A}
\title{LUXLung5, figure 1A}
\format{
A data frame of 202 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (afatinib_paclitaxel, chemo) \cr
}
}
\source{
Schuler M, Yang JC-H, Park K, et al. Afatinib beyond progression in
patients with non-small-cell lung cancer following chemotherapy,
erlotinib/gefitinib and afatinib: phase III randomized LUX-Lung 5
trial. Ann Oncol 2016; 27: 417–23.
}
\usage{
LUXLung5_1A
}
\description{
Kaplan-Meier digitized data from LUXLung5, figure 1A (PMID 26646759). A reported sample size of 220 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(LUXLung5_1A)

kmplot(LUXLung5_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Checkmate057_1A}
\alias{Checkmate057_1A}
\title{Checkmate057, figure 1A}
\format{
A data frame of 582 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d1, nivolumab) \cr
}
}
\source{
Borghaei H, Paz-Ares L, Horn L, et al. Nivolumab versus Docetaxel in
Advanced Nonsquamous Non-Small-Cell Lung Cancer. N Engl J Med 2015;
373: 1627–39.
}
\usage{
Checkmate057_1A
}
\description{
Kaplan-Meier digitized data from Checkmate057, figure 1A (PMID 26412456). A reported sample size of 582 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(Checkmate057_1A)

kmplot(Checkmate057_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{OPTIMOX3_3B}
\alias{OPTIMOX3_3B}
\title{OPTIMOX3, figure 3B}
\format{
A data frame of 452 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev, bev_erlotinib) \cr
}
}
\source{
Tournigand C, Chibaudel B, Samson B, et al. Bevacizumab with or
without erlotinib as maintenance therapy in patients with metastatic
colorectal cancer (GERCOR DREAM; OPTIMOX3): a randomised, open-label,
phase 3 trial. Lancet Oncol 2015; 16: 1493–505.
}
\usage{
OPTIMOX3_3B
}
\description{
Kaplan-Meier digitized data from OPTIMOX3, figure 3B (PMID 26474518). A reported sample size of 700 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(OPTIMOX3_3B)

kmplot(OPTIMOX3_3B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{JCOG0509_2A}
\alias{JCOG0509_2A}
\title{JCOG0509, figure 2A}
\format{
A data frame of 284 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (amrubicin_cisplatin, irinotecan_cisplatin) \cr
}
}
\source{
Satouchi M, Kotani Y, Shibata T, et al. Phase III study comparing
amrubicin plus cisplatin with irinotecan plus cisplatin in the
treatment of extensive-disease small-cell lung cancer: JCOG 0509. J
Clin Oncol 2014; 32: 1262–8.
}
\usage{
JCOG0509_2A
}
\description{
Kaplan-Meier digitized data from JCOG0509, figure 2A (PMID 24638015). A reported sample size of 284 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(JCOG0509_2A)

kmplot(JCOG0509_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PALOMA3_1A}
\alias{PALOMA3_1A}
\title{PALOMA3, figure 1A}
\format{
A data frame of 521 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (palbociclib_fulvestrant, placebo_fulvestrant) \cr
}
}
\source{
Turner NC, Ro J, André F, et al. Palbociclib in
Hormone-Receptor-Positive Advanced Breast Cancer. N Engl J Med 2015;
373: 209–19.
}
\usage{
PALOMA3_1A
}
\description{
Kaplan-Meier digitized data from PALOMA3, figure 1A (PMID 26030518). A reported sample size of 521 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(PALOMA3_1A)

kmplot(PALOMA3_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ITACa_2A}
\alias{ITACa_2A}
\title{ITACa, figure 2A}
\format{
A data frame of 370 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ct, ct_b) \cr
}
}
\source{
Passardi A, Nanni O, Tassinari D, et al. Effectiveness of bevacizumab
added to standard chemotherapy in metastatic colorectal cancer: final
results for first-line treatment from the ITACa randomized clinical
trial. Ann Oncol 2015; 26: 1201–7.
}
\usage{
ITACa_2A
}
\description{
Kaplan-Meier digitized data from ITACa, figure 2A (PMID 25735317). A reported sample size of 376 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(ITACa_2A)

kmplot(ITACa_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NSABPR04_5A}
\alias{NSABPR04_5A}
\title{NSABPR04, figure 5A}
\format{
A data frame of 1,567 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (5fu, cape) \cr
}
}
\source{
Allegra CJ, Yothers G, O’Connell MJ, et al. Neoadjuvant 5-FU or
Capecitabine Plus Radiation With or Without Oxaliplatin in Rectal
Cancer Patients: A Phase III Randomized Clinical Trial. J Natl Cancer
Inst 2015; 107. DOI:10.1093/jnci/djv248.
}
\usage{
NSABPR04_5A
}
\description{
Kaplan-Meier digitized data from NSABPR04, figure 5A (PMID 26374429). A reported sample size of 1,608 for a primary endpoint of Local-regional tumor control in colorectal cancer.
}
\examples{
summary(NSABPR04_5A)

kmplot(NSABPR04_5A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ALTTO(c)_2A}
\alias{ALTTO(c)_2A}
\title{ALTTO(c), figure 2A}
\format{
A data frame of 4,197 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l, t) \cr
}
}
\source{
Piccart-Gebhart M, Holmes E, Baselga J, et al. Adjuvant Lapatinib and
Trastuzumab for Early Human Epidermal Growth Factor Receptor
2-Positive Breast Cancer: Results From the Randomized Phase III
Adjuvant Lapatinib and/or Trastuzumab Treatment Optimization Trial. J
Clin Oncol 2016; 34: 1034–42.
}
\usage{
`ALTTO(c)_2A`
}
\description{
Kaplan-Meier digitized data from ALTTO(c), figure 2A (PMID 26598744). A reported sample size of 8,381 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`ALTTO(c)_2A`)

kmplot(`ALTTO(c)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ASPECCT_2A}
\alias{ASPECCT_2A}
\title{ASPECCT, figure 2A}
\format{
A data frame of 999 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cetuximab, panitumumab) \cr
}
}
\source{
Price TJ, Peeters M, Kim TW, et al. Panitumumab versus cetuximab in
patients with chemotherapy-refractory wild-type KRAS exon 2
metastatic colorectal cancer (ASPECCT): a randomised, multicentre,
open-label, non-inferiority phase 3 study. Lancet Oncol 2014; 15:
569–79.
}
\usage{
ASPECCT_2A
}
\description{
Kaplan-Meier digitized data from ASPECCT, figure 2A (PMID 24739896). A reported sample size of 1,010 for a primary endpoint of OS in colorectal cancer.
}
\examples{
summary(ASPECCT_2A)

kmplot(ASPECCT_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TRIBE(2)_2A}
\alias{TRIBE(2)_2A}
\title{TRIBE(2), figure 2A}
\format{
A data frame of 508 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfiri_bevacizumab, folfoxiri_bevacizumab) \cr
}
}
\source{
Mayer RJ, Van Cutsem E, Falcone A, et al. Randomized trial of TAS-102
for refractory metastatic colorectal cancer. N Engl J Med 2015; 372:
1909–19.
}
\usage{
`TRIBE(2)_2A`
}
\description{
Kaplan-Meier digitized data from TRIBE(2), figure 2A (PMID 25970050). A reported sample size of 508 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(`TRIBE(2)_2A`)

kmplot(`TRIBE(2)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUXLung6(2)_2A}
\alias{LUXLung6(2)_2A}
\title{LUXLung6(2), figure 2A}
\format{
A data frame of 364 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (afatinib, gemcitabine_cisplatin) \cr
}
}
\source{
Wu Y-L, Zhou C, Hu C-P, et al. Afatinib versus cisplatin plus
gemcitabine for first-line treatment of Asian patients with advanced
non-small-cell lung cancer harbouring EGFR mutations (LUX-Lung 6): an
open-label, randomised phase 3 trial. Lancet Oncol 2014; 15: 213–22.
}
\usage{
`LUXLung6(2)_2A`
}
\description{
Kaplan-Meier digitized data from LUXLung6(2), figure 2A (PMID 24439929). A reported sample size of 364 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(`LUXLung6(2)_2A`)

kmplot(`LUXLung6(2)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ALTTO(b)_2A}
\alias{ALTTO(b)_2A}
\title{ALTTO(b), figure 2A}
\format{
A data frame of 4,188 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (t, t_then_l) \cr
}
}
\source{
Piccart-Gebhart M, Holmes E, Baselga J, et al. Adjuvant Lapatinib and
Trastuzumab for Early Human Epidermal Growth Factor Receptor
2-Positive Breast Cancer: Results From the Randomized Phase III
Adjuvant Lapatinib and/or Trastuzumab Treatment Optimization Trial. J
Clin Oncol 2016; 34: 1034–42.
}
\usage{
`ALTTO(b)_2A`
}
\description{
Kaplan-Meier digitized data from ALTTO(b), figure 2A (PMID 26598744). A reported sample size of 8,381 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`ALTTO(b)_2A`)

kmplot(`ALTTO(b)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SAKK_2B}
\alias{SAKK_2B}
\title{SAKK, figure 2B}
\format{
A data frame of 232 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo_surgery, crt_surgery) \cr
}
}
\source{
Pless M, Stupp R, Ris H-B, et al. Induction chemoradiation in stage
IIIA/N2 non-small-cell lung cancer: a phase 3 randomised trial.
Lancet 2015; 386: 1049–56.
}
\usage{
SAKK_2B
}
\description{
Kaplan-Meier digitized data from SAKK, figure 2B (PMID 26275735). A reported sample size of 232 for a primary endpoint of EFS in lung cancer.
}
\examples{
summary(SAKK_2B)

kmplot(SAKK_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NSABPR04_3A}
\alias{NSABPR04_3A}
\title{NSABPR04, figure 3A}
\format{
A data frame of 1,567 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab LRR event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (5fu, cape) \cr
}
}
\source{
Allegra CJ, Yothers G, O’Connell MJ, et al. Neoadjuvant 5-FU or
Capecitabine Plus Radiation With or Without Oxaliplatin in Rectal
Cancer Patients: A Phase III Randomized Clinical Trial. J Natl Cancer
Inst 2015; 107. DOI:10.1093/jnci/djv248.
}
\usage{
NSABPR04_3A
}
\description{
Kaplan-Meier digitized data from NSABPR04, figure 3A (PMID 26374429). A reported sample size of 1,608 for a primary endpoint of Local-regional tumor control in colorectal cancer.
}
\examples{
summary(NSABPR04_3A)

kmplot(NSABPR04_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ELDA_1B}
\alias{ELDA_1B}
\title{ELDA, figure 1B}
\format{
A data frame of 299 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cmf, wdocetaxel) \cr
}
}
\source{
Perrone F, Nuzzo F, Di Rella F, et al. Weekly docetaxel versus CMF as
adjuvant chemotherapy for older women with early breast cancer: final
results of the randomized phase III ELDA trial. Ann Oncol 2015; 26:
675–82.
}
\usage{
ELDA_1B
}
\description{
Kaplan-Meier digitized data from ELDA, figure 1B (PMID 25488686). A reported sample size of 302 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(ELDA_1B)

kmplot(ELDA_1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{WJOG5108L_2B}
\alias{WJOG5108L_2B}
\title{WJOG5108L, figure 2B}
\format{
A data frame of 559 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (erlotinib, gefitinib) \cr
}
}
\source{
Urata Y, Katakami N, Morita S, et al. Randomized Phase III Study
Comparing Gefitinib With Erlotinib in Patients With Previously
Treated Advanced Lung Adenocarcinoma: WJOG 5108L. J Clin Oncol 2016;
34: 3248–57.
}
\usage{
WJOG5108L_2B
}
\description{
Kaplan-Meier digitized data from WJOG5108L, figure 2B (PMID 27022112). A reported sample size of 561 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(WJOG5108L_2B)

kmplot(WJOG5108L_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{select_kmdata}
\alias{select_kmdata}
\title{Select data sets}
\usage{
select_kmdata(..., return = c("name", "key", "data"))
}
\arguments{
\item{...}{an expression to be evaluated within \code{\link{kmdata_key}}
such as \code{ReportedSampleSize < 100}; see examples}

\item{return}{type of object to return; one of \code{"name"} (default) for
the study names that match the criteria, \code{"key"} for the matching
rows of \code{kmdata_key}, or \code{"data"} for a list of data frames for
each match}
}
\description{
Select publication data sets based on study characteristics including
outcome, sample size, treatment arms, journal, disease, etc.
}
\examples{
names(kmdata_key)
select_kmdata(ReportedSampleSize < 100)
select_kmdata(grepl('folfiri', Arms) & Outcome == 'OS')
select_kmdata(ReportedSampleSize < 100 |
  Cancer \%in\% c('Lung/Colorectal', 'Prostate'))

## get a list of the data sets
l <- select_kmdata(ReportedSampleSize < 100, return = 'data')
par(mfrow = n2mfrow(length(l)))
sapply(l, kmplot)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NSABPB40_3B}
\alias{NSABPB40_3B}
\title{NSABPB40, figure 3B}
\format{
A data frame of 1,186 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev, no_bev) \cr
}
}
\source{
Bear HD, Tang G, Rastogi P, et al. Neoadjuvant plus adjuvant
bevacizumab in early breast cancer (NSABP B-40 [NRG Oncology]):
secondary outcomes of a phase 3, randomised controlled trial. Lancet
Oncol 2015; 16: 1037–48.
}
\usage{
NSABPB40_3B
}
\description{
Kaplan-Meier digitized data from NSABPB40, figure 3B (PMID 26272770). A reported sample size of 1,206 for a primary endpoint of pCR in breast cancer.
}
\examples{
summary(NSABPB40_3B)

kmplot(NSABPB40_3B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MAPS_2C}
\alias{MAPS_2C}
\title{MAPS, figure 2C}
\format{
A data frame of 448 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (pc, pcb) \cr
}
}
\source{
Zalcman G, Mazieres J, Margery J, et al. Bevacizumab for newly
diagnosed pleural mesothelioma in the Mesothelioma Avastin Cisplatin
Pemetrexed Study (MAPS): a randomised, controlled, open-label, phase
3 trial. Lancet 2016; 387: 1405–14.
}
\usage{
MAPS_2C
}
\description{
Kaplan-Meier digitized data from MAPS, figure 2C (PMID 26719230). A reported sample size of 448 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(MAPS_2C)

kmplot(MAPS_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUXBreast1_2}
\alias{LUXBreast1_2}
\title{LUXBreast1, figure 2}
\format{
A data frame of 508 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (afatinib, t) \cr
}
}
\source{
Harbeck N, Huang C-S, Hurvitz S, et al. Afatinib plus vinorelbine
versus trastuzumab plus vinorelbine in patients with
HER2-overexpressing metastatic breast cancer who had progressed on
one previous trastuzumab treatment (LUX-Breast 1): an open-label,
randomised, phase 3 trial. Lancet Oncol 2016; 17: 357–66.
}
\usage{
LUXBreast1_2
}
\description{
Kaplan-Meier digitized data from LUXBreast1, figure 2 (PMID 26822398). A reported sample size of 508 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(LUXBreast1_2)

kmplot(LUXBreast1_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{FALCON_2}
\alias{FALCON_2}
\title{FALCON, figure 2}
\format{
A data frame of 462 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (anastrazole, fulvestrant) \cr
}
}
\source{
Robertson JFR, Bondarenko IM, Trishkina E, et al. Fulvestrant 500 mg
versus anastrozole 1 mg for hormone receptor-positive advanced breast
cancer (FALCON): an international, randomised, double-blind, phase 3
trial. Lancet 2016; 388: 2997–3005.
}
\usage{
FALCON_2
}
\description{
Kaplan-Meier digitized data from FALCON, figure 2 (PMID 27908454). A reported sample size of 462 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(FALCON_2)

kmplot(FALCON_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Checkmate057_1C}
\alias{Checkmate057_1C}
\title{Checkmate057, figure 1C}
\format{
A data frame of 582 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d1, nivolumab) \cr
}
}
\source{
Borghaei H, Paz-Ares L, Horn L, et al. Nivolumab versus Docetaxel in
Advanced Nonsquamous Non-Small-Cell Lung Cancer. N Engl J Med 2015;
373: 1627–39.
}
\usage{
Checkmate057_1C
}
\description{
Kaplan-Meier digitized data from Checkmate057, figure 1C (PMID 26412456). A reported sample size of 582 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(Checkmate057_1C)

kmplot(Checkmate057_1C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{E1199(b)_2A}
\alias{E1199(b)_2A}
\title{E1199(b), figure 2A}
\format{
A data frame of 2,484 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d3, p3) \cr
}
}
\source{
Sparano JA, Zhao F, Martino S, et al. Long-Term Follow-Up of the
E1199 Phase III Trial Evaluating the Role of Taxane and Schedule in
Operable Breast Cancer. J Clin Oncol 2015; 33: 2353–60.
}
\usage{
`E1199(b)_2A`
}
\description{
Kaplan-Meier digitized data from E1199(b), figure 2A (PMID 26077235). A reported sample size of 4,954 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`E1199(b)_2A`)

kmplot(`E1199(b)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUXBreast1_4}
\alias{LUXBreast1_4}
\title{LUXBreast1, figure 4}
\format{
A data frame of 508 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (afatinib, t) \cr
}
}
\source{
Harbeck N, Huang C-S, Hurvitz S, et al. Afatinib plus vinorelbine
versus trastuzumab plus vinorelbine in patients with
HER2-overexpressing metastatic breast cancer who had progressed on
one previous trastuzumab treatment (LUX-Breast 1): an open-label,
randomised, phase 3 trial. Lancet Oncol 2016; 17: 357–66.
}
\usage{
LUXBreast1_4
}
\description{
Kaplan-Meier digitized data from LUXBreast1, figure 4 (PMID 26822398). A reported sample size of 508 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(LUXBreast1_4)

kmplot(LUXBreast1_4)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TCOG0701_1B}
\alias{TCOG0701_1B}
\title{TCOG0701, figure 1B}
\format{
A data frame of 596 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (docetaxel_cis, s1_cis) \cr
}
}
\source{
Kubota K, Sakai H, Katakami N, et al. A randomized phase III trial of
oral S-1 plus cisplatin versus docetaxel plus cisplatin in Japanese
patients with advanced non-small-cell lung cancer: TCOG0701 CATS
trial. Ann Oncol 2015; 26: 1401–8.
}
\usage{
TCOG0701_1B
}
\description{
Kaplan-Meier digitized data from TCOG0701, figure 1B (PMID 25908605). A reported sample size of 608 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(TCOG0701_1B)

kmplot(TCOG0701_1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SAKK_2A}
\alias{SAKK_2A}
\title{SAKK, figure 2A}
\format{
A data frame of 232 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab EFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo_surgery, crt_surgery) \cr
}
}
\source{
Pless M, Stupp R, Ris H-B, et al. Induction chemoradiation in stage
IIIA/N2 non-small-cell lung cancer: a phase 3 randomised trial.
Lancet 2015; 386: 1049–56.
}
\usage{
SAKK_2A
}
\description{
Kaplan-Meier digitized data from SAKK, figure 2A (PMID 26275735). A reported sample size of 232 for a primary endpoint of EFS in lung cancer.
}
\examples{
summary(SAKK_2A)

kmplot(SAKK_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{E1199(b)_2C}
\alias{E1199(b)_2C}
\title{E1199(b), figure 2C}
\format{
A data frame of 2,484 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d3, p3) \cr
}
}
\source{
Sparano JA, Zhao F, Martino S, et al. Long-Term Follow-Up of the
E1199 Phase III Trial Evaluating the Role of Taxane and Schedule in
Operable Breast Cancer. J Clin Oncol 2015; 33: 2353–60.
}
\usage{
`E1199(b)_2C`
}
\description{
Kaplan-Meier digitized data from E1199(b), figure 2C (PMID 26077235). A reported sample size of 4,954 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`E1199(b)_2C`)

kmplot(`E1199(b)_2C`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MINDACT_2D}
\alias{MINDACT_2D}
\title{MINDACT, figure 2D}
\format{
A data frame of 690 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, no_chemotherapy) \cr
}
}
\source{
Cardoso F, van’t Veer LJ, Bogaerts J, et al. 70-Gene Signature as an
Aid to Treatment Decisions in Early-Stage Breast Cancer. N Engl J Med
2016; 375: 717–29.
}
\usage{
MINDACT_2D
}
\description{
Kaplan-Meier digitized data from MINDACT, figure 2D (PMID 27557300). A reported sample size of 1,550 for a primary endpoint of DMFS in breast cancer.
}
\examples{
summary(MINDACT_2D)

kmplot(MINDACT_2D)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NOAH_2B}
\alias{NOAH_2B}
\title{NOAH, figure 2B}
\format{
A data frame of 235 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, chemotherapy_trastuzumab) \cr
}
}
\source{
Gianni L, Eiermann W, Semiglazov V, et al. Neoadjuvant and adjuvant
trastuzumab in patients with HER2-positive locally advanced breast
cancer (NOAH): follow-up of a randomised controlled superiority trial
with a parallel HER2-negative cohort. Lancet Oncol 2014; 15: 640–7.
}
\usage{
NOAH_2B
}
\description{
Kaplan-Meier digitized data from NOAH, figure 2B (PMID 24657003). A reported sample size of 235 for a primary endpoint of EFS in breast cancer.
}
\examples{
summary(NOAH_2B)

kmplot(NOAH_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCICBR26_2A}
\alias{NCICBR26_2A}
\title{NCICBR26, figure 2A}
\format{
A data frame of 720 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (dacomitinib, placebo) \cr
}
}
\source{
Ellis PM, Shepherd FA, Millward M, et al. Dacomitinib compared with
placebo in pretreated patients with advanced or metastatic
non-small-cell lung cancer (NCIC CTG BR.26): a double-blind,
randomised, phase 3 trial. Lancet Oncol 2014; 15: 1379–88.
}
\usage{
NCICBR26_2A
}
\description{
Kaplan-Meier digitized data from NCICBR26, figure 2A (PMID 25439692). A reported sample size of 720 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(NCICBR26_2A)

kmplot(NCICBR26_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00614393(a)_3B}
\alias{NCT00614393(a)_3B}
\title{NCT00614393(a), figure 3B}
\format{
A data frame of 227 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (arma, armc) \cr
}
}
\source{
Sclafani F, Kim TY, Cunningham D, et al. A Randomized Phase II/III
Study of Dalotuzumab in Combination With Cetuximab and Irinotecan in
Chemorefractory, KRAS Wild-Type, Metastatic Colorectal Cancer. J Natl
Cancer Inst 2015; 107: djv258.
}
\usage{
`NCT00614393(a)_3B`
}
\description{
Kaplan-Meier digitized data from NCT00614393(a), figure 3B (PMID 26405092). A reported sample size of 344 for a primary endpoint of PFS,OS in colorectal cancer.
}
\examples{
summary(`NCT00614393(a)_3B`)

kmplot(`NCT00614393(a)_3B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ipd.R
\name{ipd}
\alias{ipd}
\title{IPD}
\usage{
ipd(
  time,
  prob,
  t.atrisk,
  n.atrisk,
  total = NA,
  arm = "arm",
  tau = max(t.atrisk)
)
}
\arguments{
\item{time, prob}{vectors of survival time and probabilities}

\item{t.atrisk, n.atrisk}{vectors of at-risk times and number at-risk}

\item{total}{(optional) total number of events}

\item{arm}{(optional) arm label}

\item{tau}{maximum follow-up time}
}
\description{
Generate individual patient data (IPD) from digitized survival curve data
using the method described by Guyot (2012).
}
\examples{
\dontrun{
xy <- system.file('etc', 'init', 'Checkmate 067_S3A_Nivolumab.csv',
                  package = 'kmdata')
ar <- system.file('etc', 'init', 'At Risk.csv', package = 'kmdata')

dd <- read.csv(xy)
aa <- read.csv(ar)
aa <- aa[aa$nrisk > 0, ]

## use full at-risk table data
ipd1 <- ipd(dd$T, dd$S, aa$trisk, aa$nrisk, arm = 'Nivo')
kmplot(ipd1, xaxis.at = 0:4 * 10, xlim = c(0, 45))


## only using 1/2 of at-risk table
idx <- replace(seq_along(aa$trisk) \%\% 2 == 1, nrow(aa), TRUE)
ipd2 <- ipd(dd$T, dd$S, aa$trisk[idx], aa$nrisk[idx], arm = 'Nivo')
kmplot(ipd2, xaxis.at = 0:4 * 10, xlim = c(0, 45))


## only using 1/4 of at-risk table
idx <- replace(seq_along(aa$trisk) \%\% 4 == 1, nrow(aa), TRUE)
ipd3 <- ipd(dd$T, dd$S, aa$trisk[idx], aa$nrisk[idx], arm = 'Nivo')
kmplot(ipd3, xaxis.at = 0:4 * 10, xlim = c(0, 45))


## only using 1/8 of at-risk table
idx <- replace(seq_along(aa$trisk) \%\% 8 == 1, nrow(aa), TRUE)
ipd4 <- ipd(dd$T, dd$S, aa$trisk[idx], aa$nrisk[idx], arm = 'Nivo')
kmplot(ipd4, xaxis.at = 0:4 * 10, xlim = c(0, 45))
}

}
\references{
Guyot, Patricia, et al. Enhanced secondary analysis of survival data:
reconstructing the data from published Kaplan-Meier survival curves.
\emph{BMC Medical Research Methodology} 2012, \strong{12}:9.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GIM2_2A}
\alias{GIM2_2A}
\title{GIM2, figure 2A}
\format{
A data frame of 2,091 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ec_p, fec_p) \cr
}
}
\source{
Del Mastro L, De Placido S, Bruzzi P, et al. Fluorouracil and
dose-dense chemotherapy in adjuvant treatment of patients with
early-stage breast cancer: an open-label, 2 × 2 factorial, randomised
phase 3 trial. Lancet 2015; 385: 1863–72.
}
\usage{
GIM2_2A
}
\description{
Kaplan-Meier digitized data from GIM2, figure 2A (PMID 25740286). A reported sample size of 2,091 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(GIM2_2A)

kmplot(GIM2_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CALGB40502_2A}
\alias{CALGB40502_2A}
\title{CALGB40502, figure 2A}
\format{
A data frame of 487 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ixabepilone, p1) \cr
}
}
\source{
Rugo HS, Barry WT, Moreno-Aspitia A, et al. Randomized Phase III
Trial of Paclitaxel Once Per Week Compared With Nanoparticle
Albumin-Bound Nab-Paclitaxel Once Per Week or Ixabepilone With
Bevacizumab As First-Line Chemotherapy for Locally Recurrent or
Metastatic Breast Cancer: CALGB 40502/NCCTG N063H (Alliance). J Clin
Oncol 2015; 33: 2361–9.
}
\usage{
CALGB40502_2A
}
\description{
Kaplan-Meier digitized data from CALGB40502, figure 2A (PMID 26056183). A reported sample size of 799 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CALGB40502_2A)

kmplot(CALGB40502_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TURANDOT_2B}
\alias{TURANDOT_2B}
\title{TURANDOT, figure 2B}
\format{
A data frame of 564 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev_capecitabine, bev_paclitaxel) \cr
}
}
\source{
Zielinski C, Láng I, Inbar M, et al. Bevacizumab plus paclitaxel
versus bevacizumab plus capecitabine as first-line treatment for
HER2-negative metastatic breast cancer (TURANDOT): primary endpoint
results of a randomised, open-label, non-inferiority, phase 3 trial.
Lancet Oncol 2016; 17: 1230–9.
}
\usage{
TURANDOT_2B
}
\description{
Kaplan-Meier digitized data from TURANDOT, figure 2B (PMID 27501767). A reported sample size of 564 for a primary endpoint of OS in breast cancer.
}
\examples{
summary(TURANDOT_2B)

kmplot(TURANDOT_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RAISE_2B}
\alias{RAISE_2B}
\title{RAISE, figure 2B}
\format{
A data frame of 1,072 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo_folfiri, ramucirumab_folfiri) \cr
}
}
\source{
Tabernero J, Yoshino T, Cohn AL, et al. Ramucirumab versus placebo in
combination with second-line FOLFIRI in patients with metastatic
colorectal carcinoma that progressed during or after first-line
therapy with bevacizumab, oxaliplatin, and a fluoropyrimidine
(RAISE): a randomised, double-blind, multicentre, phase 3 study.
Lancet Oncol 2015; 16: 499–508.
}
\usage{
RAISE_2B
}
\description{
Kaplan-Meier digitized data from RAISE, figure 2B (PMID 25877855). A reported sample size of 1,072 for a primary endpoint of OS in colorectal cancer.
}
\examples{
summary(RAISE_2B)

kmplot(RAISE_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MINDACT_2F}
\alias{MINDACT_2F}
\title{MINDACT, figure 2F}
\format{
A data frame of 690 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, no_chemotherapy) \cr
}
}
\source{
Cardoso F, van’t Veer LJ, Bogaerts J, et al. 70-Gene Signature as an
Aid to Treatment Decisions in Early-Stage Breast Cancer. N Engl J Med
2016; 375: 717–29.
}
\usage{
MINDACT_2F
}
\description{
Kaplan-Meier digitized data from MINDACT, figure 2F (PMID 27557300). A reported sample size of 1,550 for a primary endpoint of DMFS in breast cancer.
}
\examples{
summary(MINDACT_2F)

kmplot(MINDACT_2F)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CA184043_4}
\alias{CA184043_4}
\title{CA184043, figure 4}
\format{
A data frame of 799 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ipilimumab, placebo) \cr
}
}
\source{
Kwon ED, Drake CG, Scher HI, et al. Ipilimumab versus placebo after
radiotherapy in patients with metastatic castration-resistant
prostate cancer that had progressed after docetaxel chemotherapy
(CA184-043): a multicentre, randomised, double-blind, phase 3 trial.
Lancet Oncol 2014; 15: 700–12.
}
\usage{
CA184043_4
}
\description{
Kaplan-Meier digitized data from CA184043, figure 4 (PMID 24831977). A reported sample size of 799 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(CA184043_4)

kmplot(CA184043_4)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NeoALTTO(b)_2A}
\alias{NeoALTTO(b)_2A}
\title{NeoALTTO(b), figure 2A}
\format{
A data frame of 301 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab EFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (lapatinib_trastuzumab, t) \cr
}
}
\source{
de Azambuja E, Holmes AP, Piccart-Gebhart M, et al. Lapatinib with
trastuzumab for HER2-positive early breast cancer (NeoALTTO):
survival outcomes of a randomised, open-label, multicentre, phase 3
trial and their association with pathological complete response.
Lancet Oncol 2014; 15: 1137–46.
}
\usage{
`NeoALTTO(b)_2A`
}
\description{
Kaplan-Meier digitized data from NeoALTTO(b), figure 2A (PMID 25130998). A reported sample size of 455 for a primary endpoint of pCR in breast cancer.
}
\examples{
summary(`NeoALTTO(b)_2A`)

kmplot(`NeoALTTO(b)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CA097375_2}
\alias{CA097375_2}
\title{CA097375, figure 2}
\format{
A data frame of 499 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (a, b) \cr
}
}
\source{
Love RR, Laudico AV, Van Dinh N, et al. Timing of adjuvant surgical
oophorectomy in the menstrual cycle and disease-free and overall
survival in premenopausal women with operable breast cancer. J Natl
Cancer Inst 2015; 107: djv064.
}
\usage{
CA097375_2
}
\description{
Kaplan-Meier digitized data from CA097375, figure 2 (PMID 25794890). A reported sample size of 740 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(CA097375_2)

kmplot(CA097375_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{COMET1_2B}
\alias{COMET1_2B}
\title{COMET1, figure 2B}
\format{
A data frame of 1,028 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cabozantinib, prednisone) \cr
}
}
\source{
Smith M, De Bono J, Sternberg C, et al. Phase III Study of
Cabozantinib in Previously Treated Metastatic Castration-Resistant
Prostate Cancer: COMET-1. J Clin Oncol 2016; 34: 3005–13.
}
\usage{
COMET1_2B
}
\description{
Kaplan-Meier digitized data from COMET1, figure 2B (PMID 27400947). A reported sample size of 1,028 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(COMET1_2B)

kmplot(COMET1_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SOFT_2A}
\alias{SOFT_2A}
\title{SOFT, figure 2A}
\format{
A data frame of 2,033 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (tamoxifen, tamoxifen_os) \cr
}
}
\source{
Francis PA, Regan MM, Fleming GF, et al. Adjuvant ovarian suppression
in premenopausal breast cancer. N Engl J Med 2015; 372: 436–46.
}
\usage{
SOFT_2A
}
\description{
Kaplan-Meier digitized data from SOFT, figure 2A (PMID 25495490). A reported sample size of 3,066 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(SOFT_2A)

kmplot(SOFT_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{FFCD2001_S1B}
\alias{FFCD2001_S1B}
\title{FFCD2001, figure S1B}
\format{
A data frame of 282 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (fu, iri) \cr
}
}
\source{
Aparicio T, Lavau-Denes S, Phelip JM, et al. Randomized phase III
trial in elderly patients comparing LV5FU2 with or without irinotecan
for first-line treatment of metastatic colorectal cancer (FFCD
2001-02). Ann Oncol 2016; 27: 121–7.
}
\usage{
FFCD2001_S1B
}
\description{
Kaplan-Meier digitized data from FFCD2001, figure S1B (PMID 26487578). A reported sample size of 212 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(FFCD2001_S1B)

kmplot(FFCD2001_S1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ENSURE_3A}
\alias{ENSURE_3A}
\title{ENSURE, figure 3A}
\format{
A data frame of 217 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (erlotinib, gp) \cr
}
}
\source{
Wu Y-L, Zhou C, Liam C-K, et al. First-line erlotinib versus
gemcitabine/cisplatin in patients with advanced EGFR
mutation-positive non-small-cell lung cancer: analyses from the phase
III, randomized, open-label, ENSURE study. Ann Oncol 2015; 26:
1883–9.
}
\usage{
ENSURE_3A
}
\description{
Kaplan-Meier digitized data from ENSURE, figure 3A (PMID 26105600). A reported sample size of 217 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(ENSURE_3A)

kmplot(ENSURE_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{WJOG5208L_2A}
\alias{WJOG5208L_2A}
\title{WJOG5208L, figure 2A}
\format{
A data frame of 349 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cisplatin_docetaxel, nedaplatin_docetaxel) \cr
}
}
\source{
Shukuya T, Yamanaka T, Seto T, et al. Nedaplatin plus docetaxel
versus cisplatin plus docetaxel for advanced or relapsed squamous
cell carcinoma of the lung (WJOG5208L): a randomised, open-label,
phase 3 trial. Lancet Oncol 2015; 16: 1630–8.
}
\usage{
WJOG5208L_2A
}
\description{
Kaplan-Meier digitized data from WJOG5208L, figure 2A (PMID 26522337). A reported sample size of 355 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(WJOG5208L_2A)

kmplot(WJOG5208L_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AMAROS_2B}
\alias{AMAROS_2B}
\title{AMAROS, figure 2B}
\format{
A data frame of 1,425 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (dissection, radiotherapy) \cr
}
}
\source{
Donker M, van Tienhoven G, Straver ME, et al. Radiotherapy or surgery
of the axilla after a positive sentinel node in breast cancer (EORTC
10981-22023 AMAROS): a randomised, multicentre, open-label, phase 3
non-inferiority trial. Lancet Oncol 2014; 15: 1303–10.
}
\usage{
AMAROS_2B
}
\description{
Kaplan-Meier digitized data from AMAROS, figure 2B (PMID 25439688). A reported sample size of 4,823 for a primary endpoint of AxRecurrence-5yr in breast cancer.
}
\examples{
summary(AMAROS_2B)

kmplot(AMAROS_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ENSURE_2A}
\alias{ENSURE_2A}
\title{ENSURE, figure 2A}
\format{
A data frame of 217 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (erlotinib, gp) \cr
}
}
\source{
Wu Y-L, Zhou C, Liam C-K, et al. First-line erlotinib versus
gemcitabine/cisplatin in patients with advanced EGFR
mutation-positive non-small-cell lung cancer: analyses from the phase
III, randomized, open-label, ENSURE study. Ann Oncol 2015; 26:
1883–9.
}
\usage{
ENSURE_2A
}
\description{
Kaplan-Meier digitized data from ENSURE, figure 2A (PMID 26105600). A reported sample size of 217 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(ENSURE_2A)

kmplot(ENSURE_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GETUG12_2}
\alias{GETUG12_2}
\title{GETUG12, figure 2}
\format{
A data frame of 413 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab RFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (adt, adt_de) \cr
}
}
\source{
Fizazi K, Faivre L, Lesaunier F, et al. Androgen deprivation therapy
plus docetaxel and estramustine versus androgen deprivation therapy
alone for high-risk localised prostate cancer (GETUG 12): a phase 3
randomised controlled trial. Lancet Oncol 2015; 16: 787–94.
}
\usage{
GETUG12_2
}
\description{
Kaplan-Meier digitized data from GETUG12, figure 2 (PMID 26028518). A reported sample size of 207 for a primary endpoint of RFS in prostate cancer.
}
\examples{
summary(GETUG12_2)

kmplot(GETUG12_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{REVEL_3}
\alias{REVEL_3}
\title{REVEL, figure 3}
\format{
A data frame of 1,253 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo_docetaxel, ramucirumab_docetaxel) \cr
}
}
\source{
Garon EB, Ciuleanu T-E, Arrieta O, et al. Ramucirumab plus docetaxel
versus placebo plus docetaxel for second-line treatment of stage IV
non-small-cell lung cancer after disease progression on
platinum-based therapy (REVEL): a multicentre, double-blind,
randomised phase 3 trial. Lancet 2014; 384: 665–73.
}
\usage{
REVEL_3
}
\description{
Kaplan-Meier digitized data from REVEL, figure 3 (PMID 24933332). A reported sample size of 1,825 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(REVEL_3)

kmplot(REVEL_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{VANTAGE014_2}
\alias{VANTAGE014_2}
\title{VANTAGE014, figure 2}
\format{
A data frame of 661 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in weeks) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, vorinostat) \cr
}
}
\source{
Krug LM, Kindler HL, Calvert H, et al. Vorinostat in patients with
advanced malignant pleural mesothelioma who have progressed on
previous chemotherapy (VANTAGE-014): a phase 3, double-blind,
randomised, placebo-controlled trial. Lancet Oncol 2015; 16: 447–56.
}
\usage{
VANTAGE014_2
}
\description{
Kaplan-Meier digitized data from VANTAGE014, figure 2 (PMID 25800891). A reported sample size of 661 for a primary endpoint of OS/safety in lung cancer.
}
\examples{
summary(VANTAGE014_2)

kmplot(VANTAGE014_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BOLERO2_1}
\alias{BOLERO2_1}
\title{BOLERO2, figure 1}
\format{
A data frame of 724 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (eve_exe, pbo_exe) \cr
}
}
\source{
Piccart M, Hortobagyi GN, Campone M, et al. Everolimus plus
exemestane for hormone-receptor-positive, human epidermal growth
factor receptor-2-negative advanced breast cancer: overall survival
results from BOLERO-2†. Ann Oncol 2014; 25: 2357–62.
}
\usage{
BOLERO2_1
}
\description{
Kaplan-Meier digitized data from BOLERO2, figure 1 (PMID 25231953). A reported sample size of 720 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(BOLERO2_1)

kmplot(BOLERO2_1)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{QUASAR2_2B}
\alias{QUASAR2_2B}
\title{QUASAR2, figure 2B}
\format{
A data frame of 1,939 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cap, cap_bev) \cr
}
}
\source{
Kerr RS, Love S, Segelov E, et al. Adjuvant capecitabine plus
bevacizumab versus capecitabine alone in patients with colorectal
cancer (QUASAR 2): an open-label, randomised phase 3 trial. Lancet
Oncol 2016; 17: 1543–57.
}
\usage{
QUASAR2_2B
}
\description{
Kaplan-Meier digitized data from QUASAR2, figure 2B (PMID 27660192). A reported sample size of 1,952 for a primary endpoint of DFS3 in colorectal cancer.
}
\examples{
summary(QUASAR2_2B)

kmplot(QUASAR2_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUXLung5_1B}
\alias{LUXLung5_1B}
\title{LUXLung5, figure 1B}
\format{
A data frame of 202 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (afatinib_paclitaxel, chemo) \cr
}
}
\source{
Schuler M, Yang JC-H, Park K, et al. Afatinib beyond progression in
patients with non-small-cell lung cancer following chemotherapy,
erlotinib/gefitinib and afatinib: phase III randomized LUX-Lung 5
trial. Ann Oncol 2016; 27: 417–23.
}
\usage{
LUXLung5_1B
}
\description{
Kaplan-Meier digitized data from LUXLung5, figure 1B (PMID 26646759). A reported sample size of 220 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(LUXLung5_1B)

kmplot(LUXLung5_1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GECESTRO-APBI_3}
\alias{GECESTRO-APBI_3}
\title{GECESTRO-APBI, figure 3}
\format{
A data frame of 1,184 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (apbi, wbi) \cr
}
}
\source{
Strnad V, Ott OJ, Hildebrandt G, et al. 5-year results of accelerated
partial breast irradiation using sole interstitial multicatheter
brachytherapy versus whole-breast irradiation with boost after
breast-conserving surgery for low-risk invasive and in-situ carcinoma
of the female breast: a randomised, phase 3, non-inferiority trial.
Lancet 2016; 387: 229–38.
}
\usage{
`GECESTRO-APBI_3`
}
\description{
Kaplan-Meier digitized data from GECESTRO-APBI, figure 3 (PMID 26494415). A reported sample size of 551 for a primary endpoint of LR in breast cancer.
}
\examples{
summary(`GECESTRO-APBI_3`)

kmplot(`GECESTRO-APBI_3`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00614393(a)_3A}
\alias{NCT00614393(a)_3A}
\title{NCT00614393(a), figure 3A}
\format{
A data frame of 227 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (arma, armc) \cr
}
}
\source{
Sclafani F, Kim TY, Cunningham D, et al. A Randomized Phase II/III
Study of Dalotuzumab in Combination With Cetuximab and Irinotecan in
Chemorefractory, KRAS Wild-Type, Metastatic Colorectal Cancer. J Natl
Cancer Inst 2015; 107: djv258.
}
\usage{
`NCT00614393(a)_3A`
}
\description{
Kaplan-Meier digitized data from NCT00614393(a), figure 3A (PMID 26405092). A reported sample size of 344 for a primary endpoint of PFS,OS in colorectal cancer.
}
\examples{
summary(`NCT00614393(a)_3A`)

kmplot(`NCT00614393(a)_3A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{DELTA_2B}
\alias{DELTA_2B}
\title{DELTA, figure 2B}
\format{
A data frame of 301 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d1, erlotinib) \cr
}
}
\source{
Kawaguchi T, Ando M, Asami K, et al. Randomized phase III trial of
erlotinib versus docetaxel as second- or third-line therapy in
patients with advanced non-small-cell lung cancer: Docetaxel and
Erlotinib Lung Cancer Trial (DELTA). J Clin Oncol 2014; 32: 1902–8.
}
\usage{
DELTA_2B
}
\description{
Kaplan-Meier digitized data from DELTA, figure 2B (PMID 24841974). A reported sample size of 301 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(DELTA_2B)

kmplot(DELTA_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SQUIRE_2B}
\alias{SQUIRE_2B}
\title{SQUIRE, figure 2B}
\format{
A data frame of 1,093 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (gem_cis, necitumumab_gem_cis) \cr
}
}
\source{
Thatcher N, Hirsch FR, Luft AV, et al. Necitumumab plus gemcitabine
and cisplatin versus gemcitabine and cisplatin alone as first-line
therapy in patients with stage IV squamous non-small-cell lung cancer
(SQUIRE): an open-label, randomised, controlled phase 3 trial. Lancet
Oncol 2015; 16: 763–74.
}
\usage{
SQUIRE_2B
}
\description{
Kaplan-Meier digitized data from SQUIRE, figure 2B (PMID 26045340). A reported sample size of 1,093 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(SQUIRE_2B)

kmplot(SQUIRE_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PROCTOR-SCRIPT_1B}
\alias{PROCTOR-SCRIPT_1B}
\title{PROCTOR-SCRIPT, figure 1B}
\format{
A data frame of 437 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, observation) \cr
}
}
\source{
Breugom AJ, van Gijn W, Muller EW, et al. Adjuvant chemotherapy for
rectal cancer patients treated with preoperative (chemo)radiotherapy
and total mesorectal excision: a Dutch Colorectal Cancer Group (DCCG)
randomized phase III trial. Ann Oncol 2015; 26: 696–701.
}
\usage{
`PROCTOR-SCRIPT_1B`
}
\description{
Kaplan-Meier digitized data from PROCTOR-SCRIPT, figure 1B (PMID 25480874). A reported sample size of 437 for a primary endpoint of OS in colorectal cancer.
}
\examples{
summary(`PROCTOR-SCRIPT_1B`)

kmplot(`PROCTOR-SCRIPT_1B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ACTSCC_2A}
\alias{ACTSCC_2A}
\title{ACTSCC, figure 2A}
\format{
A data frame of 1,518 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (s_1, uft_lv) \cr
}
}
\source{
Yoshida M, Ishiguro M, Ikejiri K, et al. S-1 as adjuvant chemotherapy
for stage III colon cancer: a randomized phase III study (ACTS-CC
trial). Ann Oncol 2014; 25: 1743–9.
}
\usage{
ACTSCC_2A
}
\description{
Kaplan-Meier digitized data from ACTSCC, figure 2A (PMID 24942277). A reported sample size of 1,518 for a primary endpoint of DFS_3yr in colorectal cancer.
}
\examples{
summary(ACTSCC_2A)

kmplot(ACTSCC_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RTOG9804_2D}
\alias{RTOG9804_2D}
\title{RTOG9804, figure 2D}
\format{
A data frame of 585 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (observation, rt) \cr
}
}
\source{
McCormick B, Winter K, Hudis C, et al. RTOG 9804: a prospective
randomized trial for good-risk ductal carcinoma in situ comparing
radiotherapy with observation. J Clin Oncol 2015; 33: 709–15.
}
\usage{
RTOG9804_2D
}
\description{
Kaplan-Meier digitized data from RTOG9804, figure 2D (PMID 25605856). A reported sample size of 632 for a primary endpoint of ILF in breast cancer.
}
\examples{
summary(RTOG9804_2D)

kmplot(RTOG9804_2D)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BEBYP_3}
\alias{BEBYP_3}
\title{BEBYP, figure 3}
\format{
A data frame of 184 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev_ct, ct) \cr
}
}
\source{
Masi G, Salvatore L, Boni L, et al. Continuation or reintroduction of
bevacizumab beyond progression to first-line therapy in metastatic
colorectal cancer: final results of the randomized BEBYP trial. Ann
Oncol 2015; 26: 724–30.
}
\usage{
BEBYP_3
}
\description{
Kaplan-Meier digitized data from BEBYP, figure 3 (PMID 25600568). A reported sample size of 185 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(BEBYP_3)

kmplot(BEBYP_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{QUASAR2_2A}
\alias{QUASAR2_2A}
\title{QUASAR2, figure 2A}
\format{
A data frame of 1,933 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cap, cap_bev) \cr
}
}
\source{
Kerr RS, Love S, Segelov E, et al. Adjuvant capecitabine plus
bevacizumab versus capecitabine alone in patients with colorectal
cancer (QUASAR 2): an open-label, randomised phase 3 trial. Lancet
Oncol 2016; 17: 1543–57.
}
\usage{
QUASAR2_2A
}
\description{
Kaplan-Meier digitized data from QUASAR2, figure 2A (PMID 27660192). A reported sample size of 1,952 for a primary endpoint of DFS3 in colorectal cancer.
}
\examples{
summary(QUASAR2_2A)

kmplot(QUASAR2_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CALGB90202_2C}
\alias{CALGB90202_2C}
\title{CALGB90202, figure 2C}
\format{
A data frame of 645 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, za) \cr
}
}
\source{
Smith MR, Halabi S, Ryan CJ, et al. Randomized controlled trial of
early zoledronic acid in men with castration-sensitive prostate
cancer and bone metastases: results of CALGB 90202 (alliance). J Clin
Oncol 2014; 32: 1143–50.
}
\usage{
CALGB90202_2C
}
\description{
Kaplan-Meier digitized data from CALGB90202, figure 2C (PMID 24590644). A reported sample size of 645 for a primary endpoint of time_SRE in prostate cancer.
}
\examples{
summary(CALGB90202_2C)

kmplot(CALGB90202_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TANIA_2}
\alias{TANIA_2}
\title{TANIA, figure 2}
\format{
A data frame of 494 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chem_bev, chemo) \cr
}
}
\source{
von Minckwitz G, Puglisi F, Cortes J, et al. Bevacizumab plus
chemotherapy versus chemotherapy alone as second-line treatment for
patients with HER2-negative locally recurrent or metastatic breast
cancer after first-line treatment with bevacizumab plus chemotherapy
(TANIA): an open-label, randomised phase 3 trial. Lancet Oncol 2014;
15: 1269–78.
}
\usage{
TANIA_2
}
\description{
Kaplan-Meier digitized data from TANIA, figure 2 (PMID 25273342). A reported sample size of 494 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(TANIA_2)

kmplot(TANIA_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{INT0142_2B}
\alias{INT0142_2B}
\title{INT0142, figure 2B}
\format{
A data frame of 337 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (tamoxifen, tamoxifen_ofs) \cr
}
}
\source{
Tevaarwerk AJ, Wang M, Zhao F, et al. Phase III comparison of
tamoxifen versus tamoxifen plus ovarian function suppression in
premenopausal women with node-negative, hormone receptor-positive
breast cancer (E-3193, INT-0142): a trial of the Eastern Cooperative
Oncology Group. J Clin Oncol 2014; 32: 3948–58.
}
\usage{
INT0142_2B
}
\description{
Kaplan-Meier digitized data from INT0142, figure 2B (PMID 25349302). A reported sample size of 345 for a primary endpoint of DFS/OS in breast cancer.
}
\examples{
summary(INT0142_2B)

kmplot(INT0142_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LEA_2A}
\alias{LEA_2A}
\title{LEA, figure 2A}
\format{
A data frame of 374 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (et, et_b) \cr
}
}
\source{
Martín M, Loibl S, von Minckwitz G, et al. Phase III trial evaluating
the addition of bevacizumab to endocrine therapy as first-line
treatment for advanced breast cancer: the letrozole/fulvestrant and
avastin (LEA) study. J Clin Oncol 2015; 33: 1045–52.
}
\usage{
LEA_2A
}
\description{
Kaplan-Meier digitized data from LEA, figure 2A (PMID 25691671). A reported sample size of 374 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(LEA_2A)

kmplot(LEA_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ACT2_2D}
\alias{ACT2_2D}
\title{ACT2, figure 2D}
\format{
A data frame of 67 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (mut_b, mut_c) \cr
}
}
\source{
Hagman H, Frödin J-E, Berglund Å, et al. A randomized study of
KRAS-guided maintenance therapy with bevacizumab, erlotinib or
metronomic capecitabine after first-line induction treatment of
metastatic colorectal cancer: the Nordic ACT2 trial. Ann Oncol 2016;
27: 140–7.
}
\usage{
ACT2_2D
}
\description{
Kaplan-Meier digitized data from ACT2, figure 2D (PMID 26483047). A reported sample size of 233 for a primary endpoint of PFS-3mo in colorectal cancer.
}
\examples{
summary(ACT2_2D)

kmplot(ACT2_2D)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ELM-PC4_2A}
\alias{ELM-PC4_2A}
\title{ELM-PC4, figure 2A}
\format{
A data frame of 1,560 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (orteronel, placebo) \cr
}
}
\source{
Saad F, Fizazi K, Jinga V, et al. Orteronel plus prednisone in
patients with chemotherapy-naive metastatic castration-resistant
prostate cancer (ELM-PC 4): a double-blind, multicentre, phase 3,
randomised, placebo-controlled trial. Lancet Oncol 2015; 16: 338–48.
}
\usage{
`ELM-PC4_2A`
}
\description{
Kaplan-Meier digitized data from ELM-PC4, figure 2A (PMID 25701170). A reported sample size of 1,560 for a primary endpoint of rPFS/OS in prostate cancer.
}
\examples{
summary(`ELM-PC4_2A`)

kmplot(`ELM-PC4_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{METlung_3A}
\alias{METlung_3A}
\title{METlung, figure 3A}
\format{
A data frame of 499 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (onaturzumab_erlotinib, placebo_erlotinib) \cr
}
}
\source{
Spigel DR, Edelman MJ, O’Byrne K, et al. Results From the Phase III
Randomized Trial of Onartuzumab Plus Erlotinib Versus Erlotinib in
Previously Treated Stage IIIB or IV Non-Small-Cell Lung Cancer:
METLung. J Clin Oncol 2017; 35: 412–20.
}
\usage{
METlung_3A
}
\description{
Kaplan-Meier digitized data from METlung, figure 3A (PMID 27937096). A reported sample size of 499 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(METlung_3A)

kmplot(METlung_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CONCEPT_1A}
\alias{CONCEPT_1A}
\title{CONCEPT, figure 1A}
\format{
A data frame of 139 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab TTF event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (co, io) \cr
}
}
\source{
Hochster HS, Grothey A, Hart L, et al. Improved time to treatment
failure with an intermittent oxaliplatin strategy: results of
CONcePT. Ann Oncol 2014; 25: 1172–8.
}
\usage{
CONCEPT_1A
}
\description{
Kaplan-Meier digitized data from CONCEPT, figure 1A (PMID 24608198). A reported sample size of 140 for a primary endpoint of TTF in colorectal cancer.
}
\examples{
summary(CONCEPT_1A)

kmplot(CONCEPT_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ROSE_2A}
\alias{ROSE_2A}
\title{ROSE, figure 2A}
\format{
A data frame of 1,144 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (pbo_doc, ram_doc) \cr
}
}
\source{
Mackey JR, Ramos-Vazquez M, Lipatov O, et al. Primary results of
ROSE/TRIO-12, a randomized placebo-controlled phase III trial
evaluating the addition of ramucirumab to first-line docetaxel
chemotherapy in metastatic breast cancer. J Clin Oncol 2015; 33:
141–8.
}
\usage{
ROSE_2A
}
\description{
Kaplan-Meier digitized data from ROSE, figure 2A (PMID 25185099). A reported sample size of 1,144 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(ROSE_2A)

kmplot(ROSE_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NSABPB35_3}
\alias{NSABPB35_3}
\title{NSABPB35, figure 3}
\format{
A data frame of 3,077 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (anastrazole, tamoxifen) \cr
}
}
\source{
Margolese RG, Cecchini RS, Julian TB, et al. Anastrozole versus
tamoxifen in postmenopausal women with ductal carcinoma in situ
undergoing lumpectomy plus radiotherapy (NSABP B-35): a randomised,
double-blind, phase 3 clinical trial. Lancet 2016; 387: 849–56.
}
\usage{
NSABPB35_3
}
\description{
Kaplan-Meier digitized data from NSABPB35, figure 3 (PMID 26686957). A reported sample size of 3,104 for a primary endpoint of BCFS in breast cancer.
}
\examples{
summary(NSABPB35_3)

kmplot(NSABPB35_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RADIANT4_2C}
\alias{RADIANT4_2C}
\title{RADIANT4, figure 2C}
\format{
A data frame of 302 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (everolimus, placebo) \cr
}
}
\source{
Yao JC, Fazio N, Singh S, et al. Everolimus for the treatment of
advanced, non-functional neuroendocrine tumours of the lung or
gastrointestinal tract (RADIANT-4): a randomised, placebo-controlled,
phase 3 study. Lancet 2016; 387: 968–77.
}
\usage{
RADIANT4_2C
}
\description{
Kaplan-Meier digitized data from RADIANT4, figure 2C (PMID 26703889). A reported sample size of 302 for a primary endpoint of PFS in lung/colorectal cancer.
}
\examples{
summary(RADIANT4_2C)

kmplot(RADIANT4_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{QUARTZ_2A}
\alias{QUARTZ_2A}
\title{QUARTZ, figure 2A}
\format{
A data frame of 538 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in weeks) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (osc, osc_wbrt) \cr
}
}
\source{
Mulvenna P, Nankivell M, Barton R, et al. Dexamethasone and
supportive care with or without whole brain radiotherapy in treating
patients with non-small cell lung cancer with brain metastases
unsuitable for resection or stereotactic radiotherapy (QUARTZ):
results from a phase 3, non-inferiority, randomised trial. Lancet
2016; 388: 2004–14.
}
\usage{
QUARTZ_2A
}
\description{
Kaplan-Meier digitized data from QUARTZ, figure 2A (PMID 27604504). A reported sample size of 538 for a primary endpoint of QALYs in lung cancer.
}
\examples{
summary(QUARTZ_2A)

kmplot(QUARTZ_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PETACC8_3A}
\alias{PETACC8_3A}
\title{PETACC8, figure 3A}
\format{
A data frame of 1,602 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfox4, folfox4_cetuximab) \cr
}
}
\source{
Taieb J, Tabernero J, Mini E, et al. Oxaliplatin, fluorouracil, and
leucovorin with or without cetuximab in patients with resected stage
III colon cancer (PETACC-8): an open-label, randomised phase 3 trial.
Lancet Oncol 2014; 15: 862–73.
}
\usage{
PETACC8_3A
}
\description{
Kaplan-Meier digitized data from PETACC8, figure 3A (PMID 24928083). A reported sample size of 2,559 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(PETACC8_3A)

kmplot(PETACC8_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MINDACT_2C}
\alias{MINDACT_2C}
\title{MINDACT, figure 2C}
\format{
A data frame of 1,497 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, no_chemotherapy) \cr
}
}
\source{
Cardoso F, van’t Veer LJ, Bogaerts J, et al. 70-Gene Signature as an
Aid to Treatment Decisions in Early-Stage Breast Cancer. N Engl J Med
2016; 375: 717–29.
}
\usage{
MINDACT_2C
}
\description{
Kaplan-Meier digitized data from MINDACT, figure 2C (PMID 27557300). A reported sample size of 1,550 for a primary endpoint of DMFS in breast cancer.
}
\examples{
summary(MINDACT_2C)

kmplot(MINDACT_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00294996_1B}
\alias{NCT00294996_1B}
\title{NCT00294996, figure 1B}
\format{
A data frame of 363 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (mpt, pt) \cr
}
}
\source{
Baselga J, Manikhas A, Cortés J, et al. Phase III trial of
nonpegylated liposomal doxorubicin in combination with trastuzumab
and paclitaxel in HER2-positive metastatic breast cancer. Ann Oncol
2014; 25: 592–8.
}
\usage{
NCT00294996_1B
}
\description{
Kaplan-Meier digitized data from NCT00294996, figure 1B (PMID 24401928). A reported sample size of 181 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(NCT00294996_1B)

kmplot(NCT00294996_1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00596830_2A}
\alias{NCT00596830_2A}
\title{NCT00596830, figure 2A}
\format{
A data frame of 681 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, figitumumab) \cr
}
}
\source{
Langer CJ, Novello S, Park K, et al. Randomized, phase III trial of
first-line figitumumab in combination with paclitaxel and carboplatin
versus paclitaxel and carboplatin alone in patients with advanced
non-small-cell lung cancer. J Clin Oncol 2014; 32: 2059–66.
}
\usage{
NCT00596830_2A
}
\description{
Kaplan-Meier digitized data from NCT00596830, figure 2A (PMID 24888810). A reported sample size of 681 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(NCT00596830_2A)

kmplot(NCT00596830_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CEREBEL_2B}
\alias{CEREBEL_2B}
\title{CEREBEL, figure 2B}
\format{
A data frame of 540 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (lapatinib_cape, trastuzumab_cape) \cr
}
}
\source{
Pivot X, Manikhas A, Żurawski B, et al. CEREBEL (EGF111438): A Phase
III, Randomized, Open-Label Study of Lapatinib Plus Capecitabine
Versus Trastuzumab Plus Capecitabine in Patients With Human Epidermal
Growth Factor Receptor 2-Positive Metastatic Breast Cancer. J Clin
Oncol 2015; 33: 1564–73.
}
\usage{
CEREBEL_2B
}
\description{
Kaplan-Meier digitized data from CEREBEL, figure 2B (PMID 25605838). A reported sample size of 540 for a primary endpoint of Incidence_CNSmets in breast cancer.
}
\examples{
summary(CEREBEL_2B)

kmplot(CEREBEL_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{kmplot}
\alias{kmplot}
\alias{plot.kmdata}
\title{kmplot}
\usage{
kmplot(
  x,
  ...,
  relevel = FALSE,
  col = NULL,
  plot = TRUE,
  xlim = NULL,
  ylim = NULL,
  xaxis.at = NULL,
  xlab = NULL,
  ylab = NULL,
  lr_test = TRUE,
  test_details = TRUE,
  median = TRUE,
  atrisk = TRUE,
  mark.time = TRUE,
  title = NULL,
  legend = TRUE
)

\method{plot}{kmdata}(
  x,
  ...,
  relevel = FALSE,
  col = NULL,
  plot = TRUE,
  xlim = NULL,
  ylim = NULL,
  xaxis.at = NULL,
  xlab = NULL,
  ylab = NULL,
  lr_test = TRUE,
  test_details = TRUE,
  median = TRUE,
  atrisk = TRUE,
  mark.time = TRUE,
  title = NULL,
  legend = TRUE
)
}
\arguments{
\item{x}{a data set of class \code{"kmdata"}, i.e., one of the data sets
in the \pkg{kmdata} package; alternatively, a data frame with columns
"time", "event", and "arm"}

\item{...}{additional arguments passed to \code{\link[survival]{plot.survfit}}
or further to \code{\link{plot.default}} or \code{\link{par}}}

\item{relevel}{logical; if \code{TRUE}, group order is reversed; by default,
the arms are in alphabetical order which may not be desired for some
placebo or control arms}

\item{col}{a vector of colors for the survival curves (recycled for at-risk
table, medians, and legend)}

\item{plot}{logical; if \code{TRUE}, a KM figure is drawn}

\item{xlim, ylim}{x- and y-axis limits}

\item{xaxis.at}{x-axis positions of ticks and at-risk table}

\item{xlab, ylab}{the x- and y-axis labels}

\item{lr_test}{logical; if \code{TRUE}, log-rank test is shown in upper
right corner of figure}

\item{test_details}{logical; if \code{TRUE}, test statistic and degrees
of freedom are added with p-value for log-rank test}

\item{median}{logical; if \code{TRUE}, the medians for each curve is added}

\item{atrisk}{logical; if \code{TRUE}, an at-risk table is drawn below plot}

\item{mark.time}{passed to \code{\link[survival]{plot.survfit}}}

\item{title}{optional title for plot; default is \code{attr(data, 'title')}}

\item{legend}{logical; if \code{TRUE}, a legend for the arms, total events,
and hazard ratios is added}
}
\description{
Make a Kaplan-Meier plot of a data set.
}
\examples{
kmplot(ATTENTION_2A)

kmplot(
  ATTENTION_2A,
  relevel = TRUE, median = FALSE,
  col = 3:4, xaxis.at = 0:4 * 6
)

## or equivalently with ?survival::plot.survfit
s <- survfit(Surv(time, event) ~ arm, ATTENTION_2A)
plot(s, col = 1:2, mark.time = TRUE)


## kmplot can be used generically given the proper data structure
dat <- data.frame(time = aml$time, event = aml$status, arm = aml$x)
kmplot(dat)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{IMELDA_2B}
\alias{IMELDA_2B}
\title{IMELDA, figure 2B}
\format{
A data frame of 185 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bevacizumab, bevacizumab_capecitabine) \cr
}
}
\source{
Gligorov J, Doval D, Bines J, et al. Maintenance capecitabine and
bevacizumab versus bevacizumab alone after initial first-line
bevacizumab and docetaxel for patients with HER2-negative metastatic
breast cancer (IMELDA): a randomised, open-label, phase 3 trial.
Lancet Oncol 2014; 15: 1351–60.
}
\usage{
IMELDA_2B
}
\description{
Kaplan-Meier digitized data from IMELDA, figure 2B (PMID 25273343). A reported sample size of 284 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(IMELDA_2B)

kmplot(IMELDA_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BOLERO1_2B}
\alias{BOLERO1_2B}
\title{BOLERO1, figure 2B}
\format{
A data frame of 719 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (everolimus, placebo) \cr
}
}
\source{
Hurvitz SA, Andre F, Jiang Z, et al. Combination of everolimus with
trastuzumab plus paclitaxel as first-line treatment for patients with
HER2-positive advanced breast cancer (BOLERO-1): a phase 3,
randomised, double-blind, multicentre trial. Lancet Oncol 2015; 16:
816–29.
}
\usage{
BOLERO1_2B
}
\description{
Kaplan-Meier digitized data from BOLERO1, figure 2B (PMID 26092818). A reported sample size of 719 for a primary endpoint of PFS(overall & HR- population) in breast cancer.
}
\examples{
summary(BOLERO1_2B)

kmplot(BOLERO1_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BOLERO3_2}
\alias{BOLERO3_2}
\title{BOLERO3, figure 2}
\format{
A data frame of 569 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in weeks) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (eve_exe, placebo) \cr
}
}
\source{
André F, O’Regan R, Ozguroglu M, et al. Everolimus for women with
trastuzumab-resistant, HER2-positive, advanced breast cancer
(BOLERO-3): a randomised, double-blind, placebo-controlled phase 3
trial. Lancet Oncol 2014; 15: 580–91.
}
\usage{
BOLERO3_2
}
\description{
Kaplan-Meier digitized data from BOLERO3, figure 2 (PMID 24742739). A reported sample size of 569 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(BOLERO3_2)

kmplot(BOLERO3_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00673049_2B}
\alias{NCT00673049_2B}
\title{NCT00673049, figure 2B}
\format{
A data frame of 583 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, figitumumab) \cr
}
}
\source{
Scagliotti GV, Bondarenko I, Blackhall F, et al. Randomized, phase
III trial of figitumumab in combination with erlotinib versus
erlotinib alone in patients with nonadenocarcinoma nonsmall-cell lung
cancer. Ann Oncol 2015; 26: 497–504.
}
\usage{
NCT00673049_2B
}
\description{
Kaplan-Meier digitized data from NCT00673049, figure 2B (PMID 25395283). A reported sample size of 583 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(NCT00673049_2B)

kmplot(NCT00673049_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{X20050181_2B}
\alias{X20050181_2B}
\title{20050181, figure 2B}
\format{
A data frame of 486 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfiri, panitumumab_folfiri) \cr
}
}
\source{
Peeters M, Price TJ, Cervantes A, et al. Final results from a
randomized phase 3 study of FOLFIRI {+/-} panitumumab for second-line
treatment of metastatic colorectal cancer. Ann Oncol 2014; 25:
107–16.
}
\usage{
X20050181_2B
}
\description{
Kaplan-Meier digitized data from 20050181, figure 2B (PMID 24356622). A reported sample size of 1,186 for a primary endpoint of PFS/OS in colorectal cancer.
}
\examples{
summary(X20050181_2B)

kmplot(X20050181_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GEICAM2003_2B}
\alias{GEICAM2003_2B}
\title{GEICAM2003, figure 2B}
\format{
A data frame of 1,384 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ec_t, et_x) \cr
}
}
\source{
Martín M, Ruiz Simón A, Ruiz Borrego M, et al. Epirubicin Plus
Cyclophosphamide Followed by Docetaxel Versus Epirubicin Plus
Docetaxel Followed by Capecitabine As Adjuvant Therapy for
Node-Positive Early Breast Cancer: Results From the GEICAM/2003-10
Study. J Clin Oncol 2015; 33: 3788–95.
}
\usage{
GEICAM2003_2B
}
\description{
Kaplan-Meier digitized data from GEICAM2003, figure 2B (PMID 26416999). A reported sample size of 1,384 for a primary endpoint of iDFS in breast cancer.
}
\examples{
summary(GEICAM2003_2B)

kmplot(GEICAM2003_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{JCOG0605_2B}
\alias{JCOG0605_2B}
\title{JCOG0605, figure 2B}
\format{
A data frame of 180 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (combination_chemotherapy, topotecan) \cr
}
}
\source{
Goto K, Ohe Y, Shibata T, et al. Combined chemotherapy with
cisplatin, etoposide, and irinotecan versus topotecan alone as
second-line treatment for patients with sensitive relapsed small-cell
lung cancer (JCOG0605): a multicentre, open-label, randomised phase 3
trial. Lancet Oncol 2016; 17: 1147–57.
}
\usage{
JCOG0605_2B
}
\description{
Kaplan-Meier digitized data from JCOG0605, figure 2B (PMID 27312053). A reported sample size of 180 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(JCOG0605_2B)

kmplot(JCOG0605_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{X20050181_1B}
\alias{X20050181_1B}
\title{20050181, figure 1B}
\format{
A data frame of 486 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfiri, panitumumab_folfiri) \cr
}
}
\source{
Peeters M, Price TJ, Cervantes A, et al. Final results from a
randomized phase 3 study of FOLFIRI {+/-} panitumumab for second-line
treatment of metastatic colorectal cancer. Ann Oncol 2014; 25:
107–16.
}
\usage{
X20050181_1B
}
\description{
Kaplan-Meier digitized data from 20050181, figure 1B (PMID 24356622). A reported sample size of 1,186 for a primary endpoint of PFS/OS in colorectal cancer.
}
\examples{
summary(X20050181_1B)

kmplot(X20050181_1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BREC_2C}
\alias{BREC_2C}
\title{BREC, figure 2C}
\format{
A data frame of 124 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, experimental) \cr
}
}
\source{
Moran T, Wei J, Cobo M, et al. Two biomarker-directed randomized
trials in European and Chinese patients with nonsmall-cell lung
cancer: the BRCA1-RAP80 Expression Customization (BREC) studies. Ann
Oncol 2014; 25: 2147–55.
}
\usage{
BREC_2C
}
\description{
Kaplan-Meier digitized data from BREC, figure 2C (PMID 25164908). A reported sample size of 279 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(BREC_2C)

kmplot(BREC_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CA184-095_2A}
\alias{CA184-095_2A}
\title{CA184-095, figure 2A}
\format{
A data frame of 602 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ipilimumab, placebo) \cr
}
}
\source{
Beer TM, Kwon ED, Drake CG, et al. Randomized, Double-Blind, Phase
III Trial of Ipilimumab Versus Placebo in Asymptomatic or Minimally
Symptomatic Patients With Metastatic Chemotherapy-Naive
Castration-Resistant Prostate Cancer. J Clin Oncol 2017; 35: 40–7.
}
\usage{
`CA184-095_2A`
}
\description{
Kaplan-Meier digitized data from CA184-095, figure 2A (PMID 28034081). A reported sample size of 602 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(`CA184-095_2A`)

kmplot(`CA184-095_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-demo.R
\docType{data}
\name{kmdata_demo}
\alias{kmdata_demo}
\title{\code{kmdata} demographics data}
\format{
A data frame of 331 observations and 12 variables:

\tabular{lll}{

\tab \code{Publication} \tab publication identifier \cr
\tab \code{Therapy} \tab therapy received \cr
\tab \code{Description} \tab additional information about \code{Therapy}
  where appropriate \cr
\tab \code{Median Age} \tab the median age by publication/therapy \cr
\tab \code{Age Range} \tab the age range by publication/therapy \cr
\tab \code{Sex:Males} \tab the number of males by publication/therapy \cr
\tab \code{Sex:Females} \tab the number of females by publication/therapy \cr
\tab \code{Race:White} \tab the number of White/Caucasian by
  publication/therapy \cr
\tab \code{Race:Other} \tab the number of non-White/Caucasian by
  publication/therapy \cr
\tab \code{ECOG-PS 0 or 1} \tab the number of patients with an ECOG
  performance score of 0 or 1 by publication/therapy \cr
\tab \code{ECOG=PS >=2} \tab the number of patients with an ECOG
  performance score of 2+ by publication/therapy \cr
\tab \code{ECOG-PS unknown} \tab the number of patients with an ECOG
  performance score of unknown by publication/therapy \cr

}
}
\usage{
kmdata_demo
}
\description{
Demographics of data sets in \code{kmdata} package by treatment arm. Missing
values are shown when either the publication did not report or the number
could not be inferred from available information.
}
\seealso{
\code{\link{kmdata_key}}; \code{\link{select_kmdata}};
\code{\link{summary.kmdata}}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{X20050181_1A}
\alias{X20050181_1A}
\title{20050181, figure 1A}
\format{
A data frame of 597 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfiri, panitumumab_folfiri) \cr
}
}
\source{
Peeters M, Price TJ, Cervantes A, et al. Final results from a
randomized phase 3 study of FOLFIRI {+/-} panitumumab for second-line
treatment of metastatic colorectal cancer. Ann Oncol 2014; 25:
107–16.
}
\usage{
X20050181_1A
}
\description{
Kaplan-Meier digitized data from 20050181, figure 1A (PMID 24356622). A reported sample size of 1,186 for a primary endpoint of PFS/OS in colorectal cancer.
}
\examples{
summary(X20050181_1A)

kmplot(X20050181_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AIO0207(b)_3C}
\alias{AIO0207(b)_3C}
\title{AIO0207(b), figure 3C}
\format{
A data frame of 312 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (5fu_bev, bev) \cr
}
}
\source{
Hegewisch-Becker S, Graeven U, Lerchenmüller CA, et al. Maintenance
strategies after first-line oxaliplatin plus fluoropyrimidine plus
bevacizumab for patients with metastatic colorectal cancer (AIO
0207): a randomised, non-inferiority, open-label, phase 3 trial.
Lancet Oncol 2015; 16: 1355–69.
}
\usage{
`AIO0207(b)_3C`
}
\description{
Kaplan-Meier digitized data from AIO0207(b), figure 3C (PMID 26361971). A reported sample size of 472 for a primary endpoint of time_tx_failure in colorectal cancer.
}
\examples{
summary(`AIO0207(b)_3C`)

kmplot(`AIO0207(b)_3C`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GIM2_2C}
\alias{GIM2_2C}
\title{GIM2, figure 2C}
\format{
A data frame of 2,003 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (q2, q3) \cr
}
}
\source{
Del Mastro L, De Placido S, Bruzzi P, et al. Fluorouracil and
dose-dense chemotherapy in adjuvant treatment of patients with
early-stage breast cancer: an open-label, 2 × 2 factorial, randomised
phase 3 trial. Lancet 2015; 385: 1863–72.
}
\usage{
GIM2_2C
}
\description{
Kaplan-Meier digitized data from GIM2, figure 2C (PMID 25740286). A reported sample size of 2,091 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(GIM2_2C)

kmplot(GIM2_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CALGB40502_2D}
\alias{CALGB40502_2D}
\title{CALGB40502, figure 2D}
\format{
A data frame of 542 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (nab_paclitaxel, p1) \cr
}
}
\source{
Rugo HS, Barry WT, Moreno-Aspitia A, et al. Randomized Phase III
Trial of Paclitaxel Once Per Week Compared With Nanoparticle
Albumin-Bound Nab-Paclitaxel Once Per Week or Ixabepilone With
Bevacizumab As First-Line Chemotherapy for Locally Recurrent or
Metastatic Breast Cancer: CALGB 40502/NCCTG N063H (Alliance). J Clin
Oncol 2015; 33: 2361–9.
}
\usage{
CALGB40502_2D
}
\description{
Kaplan-Meier digitized data from CALGB40502, figure 2D (PMID 26056183). A reported sample size of 799 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CALGB40502_2D)

kmplot(CALGB40502_2D)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NeoALTTO(b)_2D}
\alias{NeoALTTO(b)_2D}
\title{NeoALTTO(b), figure 2D}
\format{
A data frame of 301 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (lapatinib_trastuzumab, t) \cr
}
}
\source{
de Azambuja E, Holmes AP, Piccart-Gebhart M, et al. Lapatinib with
trastuzumab for HER2-positive early breast cancer (NeoALTTO):
survival outcomes of a randomised, open-label, multicentre, phase 3
trial and their association with pathological complete response.
Lancet Oncol 2014; 15: 1137–46.
}
\usage{
`NeoALTTO(b)_2D`
}
\description{
Kaplan-Meier digitized data from NeoALTTO(b), figure 2D (PMID 25130998). A reported sample size of 455 for a primary endpoint of pCR in breast cancer.
}
\examples{
summary(`NeoALTTO(b)_2D`)

kmplot(`NeoALTTO(b)_2D`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GEICAM2003_2A}
\alias{GEICAM2003_2A}
\title{GEICAM2003, figure 2A}
\format{
A data frame of 1,384 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ec_t, et_x) \cr
}
}
\source{
Martín M, Ruiz Simón A, Ruiz Borrego M, et al. Epirubicin Plus
Cyclophosphamide Followed by Docetaxel Versus Epirubicin Plus
Docetaxel Followed by Capecitabine As Adjuvant Therapy for
Node-Positive Early Breast Cancer: Results From the GEICAM/2003-10
Study. J Clin Oncol 2015; 33: 3788–95.
}
\usage{
GEICAM2003_2A
}
\description{
Kaplan-Meier digitized data from GEICAM2003, figure 2A (PMID 26416999). A reported sample size of 1,384 for a primary endpoint of iDFS in breast cancer.
}
\examples{
summary(GEICAM2003_2A)

kmplot(GEICAM2003_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ROSE_2C}
\alias{ROSE_2C}
\title{ROSE, figure 2C}
\format{
A data frame of 1,144 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (pbo_doc, ram_doc) \cr
}
}
\source{
Mackey JR, Ramos-Vazquez M, Lipatov O, et al. Primary results of
ROSE/TRIO-12, a randomized placebo-controlled phase III trial
evaluating the addition of ramucirumab to first-line docetaxel
chemotherapy in metastatic breast cancer. J Clin Oncol 2015; 33:
141–8.
}
\usage{
ROSE_2C
}
\description{
Kaplan-Meier digitized data from ROSE, figure 2C (PMID 25185099). A reported sample size of 1,144 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(ROSE_2C)

kmplot(ROSE_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ALTTO(a)_2B}
\alias{ALTTO(a)_2B}
\title{ALTTO(a), figure 2B}
\format{
A data frame of 4,190 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l_t, t) \cr
}
}
\source{
Piccart-Gebhart M, Holmes E, Baselga J, et al. Adjuvant Lapatinib and
Trastuzumab for Early Human Epidermal Growth Factor Receptor
2-Positive Breast Cancer: Results From the Randomized Phase III
Adjuvant Lapatinib and/or Trastuzumab Treatment Optimization Trial. J
Clin Oncol 2016; 34: 1034–42.
}
\usage{
`ALTTO(a)_2B`
}
\description{
Kaplan-Meier digitized data from ALTTO(a), figure 2B (PMID 26598744). A reported sample size of 8,381 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`ALTTO(a)_2B`)

kmplot(`ALTTO(a)_2B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PREVAIL_1B}
\alias{PREVAIL_1B}
\title{PREVAIL, figure 1B}
\format{
A data frame of 1,717 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (enzalutamide, placebo) \cr
}
}
\source{
Beer TM, Armstrong AJ, Rathkopf DE, et al. Enzalutamide in metastatic
prostate cancer before chemotherapy. N Engl J Med 2014; 371: 424–33.
}
\usage{
PREVAIL_1B
}
\description{
Kaplan-Meier digitized data from PREVAIL, figure 1B (PMID 24881730). A reported sample size of 1,717 for a primary endpoint of PFS/OS in prostate cancer.
}
\examples{
summary(PREVAIL_1B)

kmplot(PREVAIL_1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00596830_2B}
\alias{NCT00596830_2B}
\title{NCT00596830, figure 2B}
\format{
A data frame of 681 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, figitumumab) \cr
}
}
\source{
Langer CJ, Novello S, Park K, et al. Randomized, phase III trial of
first-line figitumumab in combination with paclitaxel and carboplatin
versus paclitaxel and carboplatin alone in patients with advanced
non-small-cell lung cancer. J Clin Oncol 2014; 32: 2059–66.
}
\usage{
NCT00596830_2B
}
\description{
Kaplan-Meier digitized data from NCT00596830, figure 2B (PMID 24888810). A reported sample size of 681 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(NCT00596830_2B)

kmplot(NCT00596830_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUMELung1_3C}
\alias{LUMELung1_3C}
\title{LUMELung1, figure 3C}
\format{
A data frame of 1,314 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (docetaxel_nintedanib, docetaxel_placebo) \cr
}
}
\source{
Reck M, Kaiser R, Mellemgaard A, et al. Docetaxel plus nintedanib
versus docetaxel plus placebo in patients with previously treated
non-small-cell lung cancer (LUME-Lung 1): a phase 3, double-blind,
randomised controlled trial. Lancet Oncol 2014; 15: 143–55.
}
\usage{
LUMELung1_3C
}
\description{
Kaplan-Meier digitized data from LUMELung1, figure 3C (PMID 24411639). A reported sample size of 655 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(LUMELung1_3C)

kmplot(LUMELung1_3C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{FFCD2001_S1A}
\alias{FFCD2001_S1A}
\title{FFCD2001, figure S1A}
\format{
A data frame of 282 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (fu, iri) \cr
}
}
\source{
Aparicio T, Lavau-Denes S, Phelip JM, et al. Randomized phase III
trial in elderly patients comparing LV5FU2 with or without irinotecan
for first-line treatment of metastatic colorectal cancer (FFCD
2001-02). Ann Oncol 2016; 27: 121–7.
}
\usage{
FFCD2001_S1A
}
\description{
Kaplan-Meier digitized data from FFCD2001, figure S1A (PMID 26487578). A reported sample size of 212 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(FFCD2001_S1A)

kmplot(FFCD2001_S1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{COMET1_2A}
\alias{COMET1_2A}
\title{COMET1, figure 2A}
\format{
A data frame of 1,028 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cabozantinib, prednisone) \cr
}
}
\source{
Smith M, De Bono J, Sternberg C, et al. Phase III Study of
Cabozantinib in Previously Treated Metastatic Castration-Resistant
Prostate Cancer: COMET-1. J Clin Oncol 2016; 34: 3005–13.
}
\usage{
COMET1_2A
}
\description{
Kaplan-Meier digitized data from COMET1, figure 2A (PMID 27400947). A reported sample size of 1,028 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(COMET1_2A)

kmplot(COMET1_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MARIANNE(b)_2A}
\alias{MARIANNE(b)_2A}
\title{MARIANNE(b), figure 2A}
\format{
A data frame of 728 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (t-dm1_pertuzumab, trastuzumab_taxane) \cr
}
}
\source{
Perez EA, Barrios C, Eiermann W, et al. Trastuzumab Emtansine With or
Without Pertuzumab Versus Trastuzumab Plus Taxane for Human Epidermal
Growth Factor Receptor 2-Positive, Advanced Breast Cancer: Primary
Results From the Phase III MARIANNE Study. J Clin Oncol 2017; 35:
141–8.
}
\usage{
`MARIANNE(b)_2A`
}
\description{
Kaplan-Meier digitized data from MARIANNE(b), figure 2A (PMID 28056202). A reported sample size of 1,095 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(`MARIANNE(b)_2A`)

kmplot(`MARIANNE(b)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ALTTO(c)_2B}
\alias{ALTTO(c)_2B}
\title{ALTTO(c), figure 2B}
\format{
A data frame of 4,197 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l, t) \cr
}
}
\source{
Piccart-Gebhart M, Holmes E, Baselga J, et al. Adjuvant Lapatinib and
Trastuzumab for Early Human Epidermal Growth Factor Receptor
2-Positive Breast Cancer: Results From the Randomized Phase III
Adjuvant Lapatinib and/or Trastuzumab Treatment Optimization Trial. J
Clin Oncol 2016; 34: 1034–42.
}
\usage{
`ALTTO(c)_2B`
}
\description{
Kaplan-Meier digitized data from ALTTO(c), figure 2B (PMID 26598744). A reported sample size of 8,381 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`ALTTO(c)_2B`)

kmplot(`ALTTO(c)_2B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ACT1_2A}
\alias{ACT1_2A}
\title{ACT1, figure 2A}
\format{
A data frame of 637 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (amrubicin, topotecan) \cr
}
}
\source{
von Pawel J, Jotte R, Spigel DR, et al. Randomized phase III trial of
amrubicin versus topotecan as second-line treatment for patients with
small-cell lung cancer. J Clin Oncol 2014; 32: 4012–9.
}
\usage{
ACT1_2A
}
\description{
Kaplan-Meier digitized data from ACT1, figure 2A (PMID 25385727). A reported sample size of 637 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ACT1_2A)

kmplot(ACT1_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MARQUEE_2B}
\alias{MARQUEE_2B}
\title{MARQUEE, figure 2B}
\format{
A data frame of 1,048 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (erlotinib_placebo, erlotinib_tivantinib) \cr
}
}
\source{
Scagliotti G, von Pawel J, Novello S, et al. Phase III Multinational,
Randomized, Double-Blind, Placebo-Controlled Study of Tivantinib (ARQ
197) Plus Erlotinib Versus Erlotinib Alone in Previously Treated
Patients With Locally Advanced or Metastatic Nonsquamous
Non-Small-Cell Lung Cancer. J Clin Oncol 2015; 33: 2667–74.
}
\usage{
MARQUEE_2B
}
\description{
Kaplan-Meier digitized data from MARQUEE, figure 2B (PMID 26169611). A reported sample size of 1,048 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(MARQUEE_2B)

kmplot(MARQUEE_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ERCC1_2B}
\alias{ERCC1_2B}
\title{ERCC1, figure 2B}
\format{
A data frame of 468 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cisplatin_pemetrexed, p1) \cr
}
}
\source{
Lee SM, Falzon M, Blackhall F, et al. Randomized Prospective
Biomarker Trial of ERCC1 for Comparing Platinum and Nonplatinum
Therapy in Advanced Non-Small-Cell Lung Cancer: ERCC1 Trial (ET). J
Clin Oncol 2017; 35: 402–11.
}
\usage{
ERCC1_2B
}
\description{
Kaplan-Meier digitized data from ERCC1, figure 2B (PMID 27893326). A reported sample size of 648 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ERCC1_2B)

kmplot(ERCC1_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CAIRO3_2B}
\alias{CAIRO3_2B}
\title{CAIRO3, figure 2B}
\format{
A data frame of 557 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (maintenance, observation) \cr
}
}
\source{
Simkens LHJ, van Tinteren H, May A, et al. Maintenance treatment with
capecitabine and bevacizumab in metastatic colorectal cancer
(CAIRO3): a phase 3 randomised controlled trial of the Dutch
Colorectal Cancer Group. Lancet 2015; 385: 1843–52.
}
\usage{
CAIRO3_2B
}
\description{
Kaplan-Meier digitized data from CAIRO3, figure 2B (PMID 25862517). A reported sample size of 558 for a primary endpoint of PFS2 in colorectal cancer.
}
\examples{
summary(CAIRO3_2B)

kmplot(CAIRO3_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{FIRE3_2B}
\alias{FIRE3_2B}
\title{FIRE3, figure 2B}
\format{
A data frame of 592 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfiri_bevacizumab, folfiri_cetuximab) \cr
}
}
\source{
Heinemann V, von Weikersthal LF, Decker T, et al. FOLFIRI plus
cetuximab versus FOLFIRI plus bevacizumab as first-line treatment for
patients with metastatic colorectal cancer (FIRE-3): a randomised,
open-label, phase 3 trial. Lancet Oncol 2014; 15: 1065–75.
}
\usage{
FIRE3_2B
}
\description{
Kaplan-Meier digitized data from FIRE3, figure 2B (PMID 25088940). A reported sample size of 592 for a primary endpoint of objective_response in colorectal cancer.
}
\examples{
summary(FIRE3_2B)

kmplot(FIRE3_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PROFILE1014_1A}
\alias{PROFILE1014_1A}
\title{PROFILE1014, figure 1A}
\format{
A data frame of 343 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, crizotinib) \cr
}
}
\source{
Solomon BJ, Mok T, Kim D-W, et al. First-line crizotinib versus
chemotherapy in ALK-positive lung cancer. N Engl J Med 2014; 371:
2167–77.
}
\usage{
PROFILE1014_1A
}
\description{
Kaplan-Meier digitized data from PROFILE1014, figure 1A (PMID 25470694). A reported sample size of 343 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(PROFILE1014_1A)

kmplot(PROFILE1014_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MRC-RT01_2A}
\alias{MRC-RT01_2A}
\title{MRC-RT01, figure 2A}
\format{
A data frame of 836 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, escalated) \cr
}
}
\source{
Dearnaley DP, Jovic G, Syndikus I, et al. Escalated-dose versus
control-dose conformal radiotherapy for prostate cancer: long-term
results from the MRC RT01 randomised controlled trial. Lancet Oncol
2014; 15: 464–73.
}
\usage{
`MRC-RT01_2A`
}
\description{
Kaplan-Meier digitized data from MRC-RT01, figure 2A (PMID 24581940). A reported sample size of 862 for a primary endpoint of bPFS/OS in prostate cancer.
}
\examples{
summary(`MRC-RT01_2A`)

kmplot(`MRC-RT01_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MA31_2}
\alias{MA31_2}
\title{MA31, figure 2}
\format{
A data frame of 652 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ltax_l, ttax_t) \cr
}
}
\source{
Gelmon KA, Boyle FM, Kaufman B, et al. Lapatinib or Trastuzumab Plus
Taxane Therapy for Human Epidermal Growth Factor Receptor 2-Positive
Advanced Breast Cancer: Final Results of NCIC CTG MA.31. J Clin Oncol
2015; 33: 1574–83.
}
\usage{
MA31_2
}
\description{
Kaplan-Meier digitized data from MA31, figure 2 (PMID 25779558). A reported sample size of 652 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(MA31_2)

kmplot(MA31_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MIROX_2B}
\alias{MIROX_2B}
\title{MIROX, figure 2B}
\format{
A data frame of 282 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfox4, folfox7_folfiri) \cr
}
}
\source{
Hebbar M, Chibaudel B, André T, et al. FOLFOX4 versus sequential
dose-dense FOLFOX7 followed by FOLFIRI in patients with resectable
metastatic colorectal cancer (MIROX): a pragmatic approach to
chemotherapy timing with perioperative or postoperative chemotherapy
from an open-label, randomized phase III trial. Ann Oncol 2015; 26:
340–7.
}
\usage{
MIROX_2B
}
\description{
Kaplan-Meier digitized data from MIROX, figure 2B (PMID 25403578). A reported sample size of 284 for a primary endpoint of DFS-2yr in colorectal cancer.
}
\examples{
summary(MIROX_2B)

kmplot(MIROX_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT01234311_2B}
\alias{NCT01234311_2B}
\title{NCT01234311, figure 2B}
\format{
A data frame of 1,158 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, tasq) \cr
}
}
\source{
Sternberg C, Armstrong A, Pili R, et al. Randomized, Double-Blind,
Placebo-Controlled Phase III Study of Tasquinimod in Men With
Metastatic Castration-Resistant Prostate Cancer. J Clin Oncol 2016;
34: 2636–43.
}
\usage{
NCT01234311_2B
}
\description{
Kaplan-Meier digitized data from NCT01234311, figure 2B (PMID 27298414). A reported sample size of 1,245 for a primary endpoint of rPFS in prostate cancer.
}
\examples{
summary(NCT01234311_2B)

kmplot(NCT01234311_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MRC-RT01_2C}
\alias{MRC-RT01_2C}
\title{MRC-RT01, figure 2C}
\format{
A data frame of 835 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, escalated) \cr
}
}
\source{
Dearnaley DP, Jovic G, Syndikus I, et al. Escalated-dose versus
control-dose conformal radiotherapy for prostate cancer: long-term
results from the MRC RT01 randomised controlled trial. Lancet Oncol
2014; 15: 464–73.
}
\usage{
`MRC-RT01_2C`
}
\description{
Kaplan-Meier digitized data from MRC-RT01, figure 2C (PMID 24581940). A reported sample size of 862 for a primary endpoint of bPFS/OS in prostate cancer.
}
\examples{
summary(`MRC-RT01_2C`)

kmplot(`MRC-RT01_2C`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{X2005B3030_2A}
\alias{X2005B3030_2A}
\title{2005B3030, figure 2A}
\format{
A data frame of 156 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, pci) \cr
}
}
\source{
Li N, Zeng Z-F, Wang S-Y, et al. Randomized phase III trial of
prophylactic cranial irradiation versus observation in patients with
fully resected stage IIIA-N2 nonsmall-cell lung cancer and high risk
of cerebral metastases after adjuvant chemotherapy. Ann Oncol 2015;
26: 504–9.
}
\usage{
X2005B3030_2A
}
\description{
Kaplan-Meier digitized data from 2005B3030, figure 2A (PMID 25515658). A reported sample size of 156 for a primary endpoint of DFS in lung cancer.
}
\examples{
summary(X2005B3030_2A)

kmplot(X2005B3030_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ENESTg1_2A}
\alias{ENESTg1_2A}
\title{ENESTg1, figure 2A}
\format{
A data frame of 644 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (imatinib, nilotinib) \cr
}
}
\source{
Blay J-Y, Shen L, Kang Y-K, et al. Nilotinib versus imatinib as
first-line therapy for patients with unresectable or metastatic
gastrointestinal stromal tumours (ENESTg1): a randomised phase 3
trial. Lancet Oncol 2015; 16: 550–60.
}
\usage{
ENESTg1_2A
}
\description{
Kaplan-Meier digitized data from ENESTg1, figure 2A (PMID 25882987). A reported sample size of 647 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(ENESTg1_2A)

kmplot(ENESTg1_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MA17R_1A}
\alias{MA17R_1A}
\title{MA17R, figure 1A}
\format{
A data frame of 1,918 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l, placebo) \cr
}
}
\source{
Goss PE, Ingle JN, Pritchard KI, et al. Extending Aromatase-Inhibitor
Adjuvant Therapy to 10 Years. N Engl J Med 2016; 375: 209–19.
}
\usage{
MA17R_1A
}
\description{
Kaplan-Meier digitized data from MA17R, figure 1A (PMID 27264120). A reported sample size of 1,918 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(MA17R_1A)

kmplot(MA17R_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GIM_3}
\alias{GIM_3}
\title{GIM, figure 3}
\format{
A data frame of 281 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, lhrha) \cr
}
}
\source{
Lambertini M, Boni L, Michelotti A, et al. Ovarian Suppression With
Triptorelin During Adjuvant Breast Cancer Chemotherapy and Long-term
Ovarian Function, Pregnancies, and Disease-Free Survival: A
Randomized Clinical Trial. JAMA 2015; 314: 2632–40.
}
\usage{
GIM_3
}
\description{
Kaplan-Meier digitized data from GIM, figure 3 (PMID 26720025). A reported sample size of 281 for a primary endpoint of Incidence_chemo-induced menopause in breast cancer.
}
\examples{
summary(GIM_3)

kmplot(GIM_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{HYPRO_2B}
\alias{HYPRO_2B}
\title{HYPRO, figure 2B}
\format{
A data frame of 804 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (conventional fractionated, hypofractionated) \cr
}
}
\source{
Incrocci L, Wortel RC, Alemayehu WG, et al. Hypofractionated versus
conventionally fractionated radiotherapy for patients with localised
prostate cancer (HYPRO): final efficacy results from a randomised,
multicentre, open-label, phase 3 trial. Lancet Oncol 2016; 17:
1061–9.
}
\usage{
HYPRO_2B
}
\description{
Kaplan-Meier digitized data from HYPRO, figure 2B (PMID 27339116). A reported sample size of 804 for a primary endpoint of RFS in prostate cancer.
}
\examples{
summary(HYPRO_2B)

kmplot(HYPRO_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GIM2_2D}
\alias{GIM2_2D}
\title{GIM2, figure 2D}
\format{
A data frame of 2,003 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (q2, q3) \cr
}
}
\source{
Del Mastro L, De Placido S, Bruzzi P, et al. Fluorouracil and
dose-dense chemotherapy in adjuvant treatment of patients with
early-stage breast cancer: an open-label, 2 × 2 factorial, randomised
phase 3 trial. Lancet 2015; 385: 1863–72.
}
\usage{
GIM2_2D
}
\description{
Kaplan-Meier digitized data from GIM2, figure 2D (PMID 25740286). A reported sample size of 2,091 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(GIM2_2D)

kmplot(GIM2_2D)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CA184-095_3}
\alias{CA184-095_3}
\title{CA184-095, figure 3}
\format{
A data frame of 602 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ipilimumab, placebo) \cr
}
}
\source{
Beer TM, Kwon ED, Drake CG, et al. Randomized, Double-Blind, Phase
III Trial of Ipilimumab Versus Placebo in Asymptomatic or Minimally
Symptomatic Patients With Metastatic Chemotherapy-Naive
Castration-Resistant Prostate Cancer. J Clin Oncol 2017; 35: 40–7.
}
\usage{
`CA184-095_3`
}
\description{
Kaplan-Meier digitized data from CA184-095, figure 3 (PMID 28034081). A reported sample size of 602 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(`CA184-095_3`)

kmplot(`CA184-095_3`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{IBCSG2200_2A}
\alias{IBCSG2200_2A}
\title{IBCSG2200, figure 2A}
\format{
A data frame of 1,081 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cm, no cm) \cr
}
}
\source{
Colleoni M, Gray KP, Gelber S, et al. Low-Dose Oral Cyclophosphamide
and Methotrexate Maintenance for Hormone Receptor-Negative Early
Breast Cancer: International Breast Cancer Study Group Trial 22-00. J
Clin Oncol 2016; 34: 3400–8.
}
\usage{
IBCSG2200_2A
}
\description{
Kaplan-Meier digitized data from IBCSG2200, figure 2A (PMID 27325862). A reported sample size of 1,086 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(IBCSG2200_2A)

kmplot(IBCSG2200_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUXLung8_2A}
\alias{LUXLung8_2A}
\title{LUXLung8, figure 2A}
\format{
A data frame of 795 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (afatinib, erlotinib) \cr
}
}
\source{
Soria J-C, Felip E, Cobo M, et al. Afatinib versus erlotinib as
second-line treatment of patients with advanced squamous cell
carcinoma of the lung (LUX-Lung 8): an open-label randomised
controlled phase 3 trial. Lancet Oncol 2015; 16: 897–907.
}
\usage{
LUXLung8_2A
}
\description{
Kaplan-Meier digitized data from LUXLung8, figure 2A (PMID 26156651). A reported sample size of 795 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(LUXLung8_2A)

kmplot(LUXLung8_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NeoALTTO(a)_2D}
\alias{NeoALTTO(a)_2D}
\title{NeoALTTO(a), figure 2D}
\format{
A data frame of 303 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l, t) \cr
}
}
\source{
de Azambuja E, Holmes AP, Piccart-Gebhart M, et al. Lapatinib with
trastuzumab for HER2-positive early breast cancer (NeoALTTO):
survival outcomes of a randomised, open-label, multicentre, phase 3
trial and their association with pathological complete response.
Lancet Oncol 2014; 15: 1137–46.
}
\usage{
`NeoALTTO(a)_2D`
}
\description{
Kaplan-Meier digitized data from NeoALTTO(a), figure 2D (PMID 25130998). A reported sample size of 455 for a primary endpoint of pCR in breast cancer.
}
\examples{
summary(`NeoALTTO(a)_2D`)

kmplot(`NeoALTTO(a)_2D`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RTOG9804_2C}
\alias{RTOG9804_2C}
\title{RTOG9804, figure 2C}
\format{
A data frame of 585 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (observation, rt) \cr
}
}
\source{
McCormick B, Winter K, Hudis C, et al. RTOG 9804: a prospective
randomized trial for good-risk ductal carcinoma in situ comparing
radiotherapy with observation. J Clin Oncol 2015; 33: 709–15.
}
\usage{
RTOG9804_2C
}
\description{
Kaplan-Meier digitized data from RTOG9804, figure 2C (PMID 25605856). A reported sample size of 632 for a primary endpoint of ILF in breast cancer.
}
\examples{
summary(RTOG9804_2C)

kmplot(RTOG9804_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{R98_2}
\alias{R98_2}
\title{R98, figure 2}
\format{
A data frame of 357 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (5fulv, 5fulv_cpt11) \cr
}
}
\source{
Delbaldo C, Ychou M, Zawadi A, et al. Postoperative irinotecan in
resected stage II-III rectal cancer: final analysis of the French R98
Intergroup trial†. Ann Oncol 2015; 26: 1208–15.
}
\usage{
R98_2
}
\description{
Kaplan-Meier digitized data from R98, figure 2 (PMID 25739671). A reported sample size of 357 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(R98_2)

kmplot(R98_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PROCLAIM_2B}
\alias{PROCLAIM_2B}
\title{PROCLAIM, figure 2B}
\format{
A data frame of 598 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (eto_cis, pem_cis) \cr
}
}
\source{
Senan S, Brade A, Wang L-H, et al. PROCLAIM: Randomized Phase III
Trial of Pemetrexed-Cisplatin or Etoposide-Cisplatin Plus Thoracic
Radiation Therapy Followed by Consolidation Chemotherapy in Locally
Advanced Nonsquamous Non-Small-Cell Lung Cancer. J Clin Oncol 2016;
34: 953–62.
}
\usage{
PROCLAIM_2B
}
\description{
Kaplan-Meier digitized data from PROCLAIM, figure 2B (PMID 26811519). A reported sample size of 598 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(PROCLAIM_2B)

kmplot(PROCLAIM_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{WJOG5208L_2B}
\alias{WJOG5208L_2B}
\title{WJOG5208L, figure 2B}
\format{
A data frame of 349 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cisplatin_docetaxel, nedaplatin_docetaxel) \cr
}
}
\source{
Shukuya T, Yamanaka T, Seto T, et al. Nedaplatin plus docetaxel
versus cisplatin plus docetaxel for advanced or relapsed squamous
cell carcinoma of the lung (WJOG5208L): a randomised, open-label,
phase 3 trial. Lancet Oncol 2015; 16: 1630–8.
}
\usage{
WJOG5208L_2B
}
\description{
Kaplan-Meier digitized data from WJOG5208L, figure 2B (PMID 26522337). A reported sample size of 355 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(WJOG5208L_2B)

kmplot(WJOG5208L_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RTOG0617_2B}
\alias{RTOG0617_2B}
\title{RTOG0617, figure 2B}
\format{
A data frame of 465 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cetuximab, no_cetuximab) \cr
}
}
\source{
Bradley JD, Paulus R, Komaki R, et al. Standard-dose versus high-dose
conformal radiotherapy with concurrent and consolidation carboplatin
plus paclitaxel with or without cetuximab for patients with stage
IIIA or IIIB non-small-cell lung cancer (RTOG 0617): a randomised,
two-by-two factorial phase 3 study. Lancet Oncol 2015; 16: 187–99.
}
\usage{
RTOG0617_2B
}
\description{
Kaplan-Meier digitized data from RTOG0617, figure 2B (PMID 25601342). A reported sample size of 166 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(RTOG0617_2B)

kmplot(RTOG0617_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ELM-PC5_2B}
\alias{ELM-PC5_2B}
\title{ELM-PC5, figure 2B}
\format{
A data frame of 1,099 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (orteronel, placebo) \cr
}
}
\source{
Fizazi K, Jones R, Oudard S, et al. Phase III, randomized,
double-blind, multicenter trial comparing orteronel (TAK-700) plus
prednisone with placebo plus prednisone in patients with metastatic
castration-resistant prostate cancer that has progressed during or
after docetaxel-based therapy: ELM-PC 5. J Clin Oncol 2015; 33:
723–31.
}
\usage{
`ELM-PC5_2B`
}
\description{
Kaplan-Meier digitized data from ELM-PC5, figure 2B (PMID 25624429). A reported sample size of 1,099 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(`ELM-PC5_2B`)

kmplot(`ELM-PC5_2B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Chronicle_2A}
\alias{Chronicle_2A}
\title{Chronicle, figure 2A}
\format{
A data frame of 113 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (capecitabine_oxaliplatin, follow_up_only) \cr
}
}
\source{
Glynne-Jones R, Counsell N, Quirke P, et al. Chronicle: results of a
randomised phase III trial in locally advanced rectal cancer after
neoadjuvant chemoradiation randomising postoperative adjuvant
capecitabine plus oxaliplatin (XELOX) versus control. Ann Oncol 2014;
25: 1356–62.
}
\usage{
Chronicle_2A
}
\description{
Kaplan-Meier digitized data from Chronicle, figure 2A (PMID 24718885). A reported sample size of 113 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(Chronicle_2A)

kmplot(Chronicle_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AIO0207(a)_3B}
\alias{AIO0207(a)_3B}
\title{AIO0207(a), figure 3B}
\format{
A data frame of 301 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (5fu_bev, no_tx) \cr
}
}
\source{
Hegewisch-Becker S, Graeven U, Lerchenmüller CA, et al. Maintenance
strategies after first-line oxaliplatin plus fluoropyrimidine plus
bevacizumab for patients with metastatic colorectal cancer (AIO
0207): a randomised, non-inferiority, open-label, phase 3 trial.
Lancet Oncol 2015; 16: 1355–69.
}
\usage{
`AIO0207(a)_3B`
}
\description{
Kaplan-Meier digitized data from AIO0207(a), figure 3B (PMID 26361971). A reported sample size of 472 for a primary endpoint of time_tx_failure in colorectal cancer.
}
\examples{
summary(`AIO0207(a)_3B`)

kmplot(`AIO0207(a)_3B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{DART0105_2B}
\alias{DART0105_2B}
\title{DART0105, figure 2B}
\format{
A data frame of 347 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ltad, stad) \cr
}
}
\source{
Zapatero A, Guerrero A, Maldonado X, et al. High-dose radiotherapy
with short-term or long-term androgen deprivation in localised
prostate cancer (DART01/05 GICOR): a randomised, controlled, phase 3
trial. Lancet Oncol 2015; 16: 320–7.
}
\usage{
DART0105_2B
}
\description{
Kaplan-Meier digitized data from DART0105, figure 2B (PMID 25702876). A reported sample size of 178 for a primary endpoint of bDFS in prostate cancer.
}
\examples{
summary(DART0105_2B)

kmplot(DART0105_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{JFMC33-0502_2B}
\alias{JFMC33-0502_2B}
\title{JFMC33-0502, figure 2B}
\format{
A data frame of 1,060 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, study) \cr
}
}
\source{
Sadahiro S, Tsuchiya T, Sasaki K, et al. Randomized phase III trial
of treatment duration for oral uracil and tegafur plus leucovorin as
adjuvant chemotherapy for patients with stage IIB/III colon cancer:
final results of JFMC33-0502. Ann Oncol 2015; 26: 2274–80.
}
\usage{
`JFMC33-0502_2B`
}
\description{
Kaplan-Meier digitized data from JFMC33-0502, figure 2B (PMID 26347106). A reported sample size of 1,071 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(`JFMC33-0502_2B`)

kmplot(`JFMC33-0502_2B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CA184043_2A}
\alias{CA184043_2A}
\title{CA184043, figure 2A}
\format{
A data frame of 799 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ipilimumab, placebo) \cr
}
}
\source{
Kwon ED, Drake CG, Scher HI, et al. Ipilimumab versus placebo after
radiotherapy in patients with metastatic castration-resistant
prostate cancer that had progressed after docetaxel chemotherapy
(CA184-043): a multicentre, randomised, double-blind, phase 3 trial.
Lancet Oncol 2014; 15: 700–12.
}
\usage{
CA184043_2A
}
\description{
Kaplan-Meier digitized data from CA184043, figure 2A (PMID 24831977). A reported sample size of 799 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(CA184043_2A)

kmplot(CA184043_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CALGB40502_2C}
\alias{CALGB40502_2C}
\title{CALGB40502, figure 2C}
\format{
A data frame of 487 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ixabepilone, p1) \cr
}
}
\source{
Rugo HS, Barry WT, Moreno-Aspitia A, et al. Randomized Phase III
Trial of Paclitaxel Once Per Week Compared With Nanoparticle
Albumin-Bound Nab-Paclitaxel Once Per Week or Ixabepilone With
Bevacizumab As First-Line Chemotherapy for Locally Recurrent or
Metastatic Breast Cancer: CALGB 40502/NCCTG N063H (Alliance). J Clin
Oncol 2015; 33: 2361–9.
}
\usage{
CALGB40502_2C
}
\description{
Kaplan-Meier digitized data from CALGB40502, figure 2C (PMID 26056183). A reported sample size of 799 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CALGB40502_2C)

kmplot(CALGB40502_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{INSPIRE_2B}
\alias{INSPIRE_2B}
\title{INSPIRE, figure 2B}
\format{
A data frame of 633 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ct, necitumumab_ct) \cr
}
}
\source{
Paz-Ares L, Mezger J, Ciuleanu TE, et al. Necitumumab plus pemetrexed
and cisplatin as first-line therapy in patients with stage IV
non-squamous non-small-cell lung cancer (INSPIRE): an open-label,
randomised, controlled phase 3 study. Lancet Oncol 2015; 16: 328–37.
}
\usage{
INSPIRE_2B
}
\description{
Kaplan-Meier digitized data from INSPIRE, figure 2B (PMID 25701171). A reported sample size of 633 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(INSPIRE_2B)

kmplot(INSPIRE_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ESPATUE_3}
\alias{ESPATUE_3}
\title{ESPATUE, figure 3}
\format{
A data frame of 161 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (arma, armb) \cr
}
}
\source{
Eberhardt WEE, Pöttgen C, Gauler TC, et al. Phase III Study of
Surgery Versus Definitive Concurrent Chemoradiotherapy Boost in
Patients With Resectable Stage IIIA(N2) and Selected IIIB
Non-Small-Cell Lung Cancer After Induction Chemotherapy and
Concurrent Chemoradiotherapy (ESPATUE). J Clin Oncol 2015; 33:
4194–201.
}
\usage{
ESPATUE_3
}
\description{
Kaplan-Meier digitized data from ESPATUE, figure 3 (PMID 26527789). A reported sample size of 246 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ESPATUE_3)

kmplot(ESPATUE_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Checkmate017_1}
\alias{Checkmate017_1}
\title{Checkmate017, figure 1}
\format{
A data frame of 272 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d1, nivolumab) \cr
}
}
\source{
Brahmer J, Reckamp KL, Baas P, et al. Nivolumab versus Docetaxel in
Advanced Squamous-Cell Non-Small-Cell Lung Cancer. N Engl J Med 2015;
373: 123–35.
}
\usage{
Checkmate017_1
}
\description{
Kaplan-Meier digitized data from Checkmate017, figure 1 (PMID 26028407). A reported sample size of 272 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(Checkmate017_1)

kmplot(Checkmate017_1)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUXLung8_3A}
\alias{LUXLung8_3A}
\title{LUXLung8, figure 3A}
\format{
A data frame of 795 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (afatinib, erlotinib) \cr
}
}
\source{
Soria J-C, Felip E, Cobo M, et al. Afatinib versus erlotinib as
second-line treatment of patients with advanced squamous cell
carcinoma of the lung (LUX-Lung 8): an open-label randomised
controlled phase 3 trial. Lancet Oncol 2015; 16: 897–907.
}
\usage{
LUXLung8_3A
}
\description{
Kaplan-Meier digitized data from LUXLung8, figure 3A (PMID 26156651). A reported sample size of 795 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(LUXLung8_3A)

kmplot(LUXLung8_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MIROX_2A}
\alias{MIROX_2A}
\title{MIROX, figure 2A}
\format{
A data frame of 282 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfox4, folfox7_folfiri) \cr
}
}
\source{
Hebbar M, Chibaudel B, André T, et al. FOLFOX4 versus sequential
dose-dense FOLFOX7 followed by FOLFIRI in patients with resectable
metastatic colorectal cancer (MIROX): a pragmatic approach to
chemotherapy timing with perioperative or postoperative chemotherapy
from an open-label, randomized phase III trial. Ann Oncol 2015; 26:
340–7.
}
\usage{
MIROX_2A
}
\description{
Kaplan-Meier digitized data from MIROX, figure 2A (PMID 25403578). A reported sample size of 284 for a primary endpoint of DFS-2yr in colorectal cancer.
}
\examples{
summary(MIROX_2A)

kmplot(MIROX_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00337103_2B}
\alias{NCT00337103_2B}
\title{NCT00337103, figure 2B}
\format{
A data frame of 1,102 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cap, eribulin) \cr
}
}
\source{
Kaufman PA, Awada A, Twelves C, et al. Phase III open-label
randomized study of eribulin mesylate versus capecitabine in patients
with locally advanced or metastatic breast cancer previously treated
with an anthracycline and a taxane. J Clin Oncol 2015; 33: 594–601.
}
\usage{
NCT00337103_2B
}
\description{
Kaplan-Meier digitized data from NCT00337103, figure 2B (PMID 25605862). A reported sample size of 1,102 for a primary endpoint of OS/PFS in breast cancer.
}
\examples{
summary(NCT00337103_2B)

kmplot(NCT00337103_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-demo.R
\docType{data}
\name{kmdata_key}
\alias{kmdata_key}
\title{\code{kmdata} key}
\format{
A data frame of 304 observations and 29 variables:

\tabular{lll}{

\tab \code{name} \tab a unique data set identifier; use this to access data
  sets for analysis \cr
\tab \code{Journal} \tab the publishing journal \cr
\tab \code{PubMedID} \tab the publication PubMed identifier \cr
\tab \code{TrialRegistration} \tab trial registration identifier, if
  applicable; for example, the ClinicalTrials.gov identifier (NCT number)
  or other \cr
\tab \code{Title} \tab publication title \cr
\tab \code{Year} \tab publication year \cr
\tab \code{ClinicalTrial} \tab unique clinical trial identifier; note that
  \code{ClinicalTrial} and \code{Figure} are joined to create the unique
  trial-figure identifier, \code{name} \cr
\tab \code{Figure} \tab figure identifier; note that \code{ClinicalTrial}
  and \code{Figure} are joined to create the unique trial-figure
  identifier, \code{name} \cr
\tab \code{Cancer} \tab study cancer type \cr
\tab \code{Subgroups} \tab subgroups studied, if any \cr
\tab \code{ReportedSampleSize} \tab overall study sample size \cr
\tab \code{RandomizationRatio} \tab the ratio of patients randomized to
  each study arm \cr
\tab \code{RandomizationType} \tab treatment characteristic by which
  patients were randomized, e.g., patients randomized to receive different
  dose levels, durations or timing of therapy, or therapy agent \cr
\tab \code{TrialDesign} \tab trial design, i.e., superiority,
  inferiority, or both \cr
\tab \code{Metastatic} \tab indicator of metastatic or non-metastatic
  cancer \cr
\tab \code{InterventionClass} \tab the class of drug/treatment used in
  the intervention/treatment arm \cr
\tab \code{Outcome} \tab the outcome or event for each data set, e.g., OS,
  PFS, DFS, etc. \cr
\tab \code{PrimaryStudyResults} \tab the primary outcome results, either
  positive, negative, or mixed (i.e., co-primary endpoints with both
  positive and negative outcomes) \cr
\tab \code{Units} \tab the time units for \code{Outcome}, e.g., days,
  weeks, months, years \cr
\tab \code{HazardRatio} \tab the reported/published hazard ratio for
\code{Outcome} \cr
\tab \code{Arms} \tab the two treatment arms for each data set \cr
\tab \code{NumberofArms} \tab the number of treatment arms for the entire
  study (note that data sets in this package have two arms each) \cr
\tab \code{qsAtRisk} \tab quality score for at-risk tables (0-3) \cr
\tab \code{qsTotalEvents} \tab quality score total events (0-3) \cr
\tab \code{qsHazardRatio} \tab quality score for hazard ratio (0-3) \cr
\tab \code{qsMedian} \tab quality score for median time-to-event (0-3) \cr
\tab \code{QualityScore} \tab aggregate quality score for each figure
  over the four metrics (hazard ratio, total events, median time-to-event,
  number at-risk) comparing published to re-capitulated data \cr
\tab \code{QualityMax} \tab maximum quality score possible for each figure;
  each quality metric can be score from 0 (worst) to 3 (best); the maximum
  score possible for each figure depends on the number of metrics reported
  by the publication \cr
\tab \code{QualityPercent} \tab a quality score percentage for each figure,
  from 0 (worst) to 100 (best): \code{QualityScore / QualityMax * 100} \cr

}
}
\usage{
kmdata_key
}
\description{
Description of data sets available in \code{kmdata} package.
}
\seealso{
\code{\link{kmdata_demo}}; \code{\link{select_kmdata}};
\code{\link{summary.kmdata}}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{HYPRO_2A}
\alias{HYPRO_2A}
\title{HYPRO, figure 2A}
\format{
A data frame of 804 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab RFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (conventional fractionated, hypofractionated) \cr
}
}
\source{
Incrocci L, Wortel RC, Alemayehu WG, et al. Hypofractionated versus
conventionally fractionated radiotherapy for patients with localised
prostate cancer (HYPRO): final efficacy results from a randomised,
multicentre, open-label, phase 3 trial. Lancet Oncol 2016; 17:
1061–9.
}
\usage{
HYPRO_2A
}
\description{
Kaplan-Meier digitized data from HYPRO, figure 2A (PMID 27339116). A reported sample size of 804 for a primary endpoint of RFS in prostate cancer.
}
\examples{
summary(HYPRO_2A)

kmplot(HYPRO_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SQUIRE_2A}
\alias{SQUIRE_2A}
\title{SQUIRE, figure 2A}
\format{
A data frame of 1,093 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (gem_cis, necitumumab_gem_cis) \cr
}
}
\source{
Thatcher N, Hirsch FR, Luft AV, et al. Necitumumab plus gemcitabine
and cisplatin versus gemcitabine and cisplatin alone as first-line
therapy in patients with stage IV squamous non-small-cell lung cancer
(SQUIRE): an open-label, randomised, controlled phase 3 trial. Lancet
Oncol 2015; 16: 763–74.
}
\usage{
SQUIRE_2A
}
\description{
Kaplan-Meier digitized data from SQUIRE, figure 2A (PMID 26045340). A reported sample size of 1,093 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(SQUIRE_2A)

kmplot(SQUIRE_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RAISE_2A}
\alias{RAISE_2A}
\title{RAISE, figure 2A}
\format{
A data frame of 1,072 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo_folfiri, ramucirumab_folfiri) \cr
}
}
\source{
Tabernero J, Yoshino T, Cohn AL, et al. Ramucirumab versus placebo in
combination with second-line FOLFIRI in patients with metastatic
colorectal carcinoma that progressed during or after first-line
therapy with bevacizumab, oxaliplatin, and a fluoropyrimidine
(RAISE): a randomised, double-blind, multicentre, phase 3 study.
Lancet Oncol 2015; 16: 499–508.
}
\usage{
RAISE_2A
}
\description{
Kaplan-Meier digitized data from RAISE, figure 2A (PMID 25877855). A reported sample size of 1,072 for a primary endpoint of OS in colorectal cancer.
}
\examples{
summary(RAISE_2A)

kmplot(RAISE_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CHHiP(a)_2A}
\alias{CHHiP(a)_2A}
\title{CHHiP(a), figure 2A}
\format{
A data frame of 2,139 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab RFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (60gy, 74gy) \cr
}
}
\source{
Dearnaley D, Syndikus I, Mossop H, et al. Conventional versus
hypofractionated high-dose intensity-modulated radiotherapy for
prostate cancer: 5-year outcomes of the randomised, non-inferiority,
phase 3 CHHiP trial. Lancet Oncol 2016; 17: 1047–60.
}
\usage{
`CHHiP(a)_2A`
}
\description{
Kaplan-Meier digitized data from CHHiP(a), figure 2A (PMID 27339115). A reported sample size of 3,216 for a primary endpoint of bRFS in prostate cancer.
}
\examples{
summary(`CHHiP(a)_2A`)

kmplot(`CHHiP(a)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{REVEL_2}
\alias{REVEL_2}
\title{REVEL, figure 2}
\format{
A data frame of 1,253 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo_docetaxel, ramucirumab_docetaxel) \cr
}
}
\source{
Garon EB, Ciuleanu T-E, Arrieta O, et al. Ramucirumab plus docetaxel
versus placebo plus docetaxel for second-line treatment of stage IV
non-small-cell lung cancer after disease progression on
platinum-based therapy (REVEL): a multicentre, double-blind,
randomised phase 3 trial. Lancet 2014; 384: 665–73.
}
\usage{
REVEL_2
}
\description{
Kaplan-Meier digitized data from REVEL, figure 2 (PMID 24933332). A reported sample size of 1,825 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(REVEL_2)

kmplot(REVEL_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{E1199(a)_2A}
\alias{E1199(a)_2A}
\title{E1199(a), figure 2A}
\format{
A data frame of 2,481 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (p1, p3) \cr
}
}
\source{
Sparano JA, Zhao F, Martino S, et al. Long-Term Follow-Up of the
E1199 Phase III Trial Evaluating the Role of Taxane and Schedule in
Operable Breast Cancer. J Clin Oncol 2015; 33: 2353–60.
}
\usage{
`E1199(a)_2A`
}
\description{
Kaplan-Meier digitized data from E1199(a), figure 2A (PMID 26077235). A reported sample size of 4,954 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`E1199(a)_2A`)

kmplot(`E1199(a)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ESPATUE_2}
\alias{ESPATUE_2}
\title{ESPATUE, figure 2}
\format{
A data frame of 161 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (arma, armb) \cr
}
}
\source{
Eberhardt WEE, Pöttgen C, Gauler TC, et al. Phase III Study of
Surgery Versus Definitive Concurrent Chemoradiotherapy Boost in
Patients With Resectable Stage IIIA(N2) and Selected IIIB
Non-Small-Cell Lung Cancer After Induction Chemotherapy and
Concurrent Chemoradiotherapy (ESPATUE). J Clin Oncol 2015; 33:
4194–201.
}
\usage{
ESPATUE_2
}
\description{
Kaplan-Meier digitized data from ESPATUE, figure 2 (PMID 26527789). A reported sample size of 246 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ESPATUE_2)

kmplot(ESPATUE_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{WSG-ARA_3}
\alias{WSG-ARA_3}
\title{WSG-ARA, figure 3}
\format{
A data frame of 1,170 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab EFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (da_negative, da_positive) \cr
}
}
\source{
Nitz U, Gluz O, Zuna I, et al. Final results from the prospective
phase III WSG-ARA trial: impact of adjuvant darbepoetin alfa on
event-free survival in early breast cancer. Ann Oncol 2014; 25:
75–80.
}
\usage{
`WSG-ARA_3`
}
\description{
Kaplan-Meier digitized data from WSG-ARA, figure 3 (PMID 24356620). A reported sample size of 1,234 for a primary endpoint of EFS in breast cancer.
}
\examples{
summary(`WSG-ARA_3`)

kmplot(`WSG-ARA_3`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NO16968_2C}
\alias{NO16968_2C}
\title{NO16968, figure 2C}
\format{
A data frame of 1,886 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (fu_fa, xelox) \cr
}
}
\source{
Schmoll H-J, Tabernero J, Maroun J, et al. Capecitabine Plus
Oxaliplatin Compared With Fluorouracil/Folinic Acid As Adjuvant
Therapy for Stage III Colon Cancer: Final Results of the NO16968
Randomized Controlled Phase III Trial. J Clin Oncol 2015; 33:
3733–40.
}
\usage{
NO16968_2C
}
\description{
Kaplan-Meier digitized data from NO16968, figure 2C (PMID 26324362). A reported sample size of 1,886 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(NO16968_2C)

kmplot(NO16968_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{FFCD2001_S1C}
\alias{FFCD2001_S1C}
\title{FFCD2001, figure S1C}
\format{
A data frame of 282 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (classic, simplified) \cr
}
}
\source{
Aparicio T, Lavau-Denes S, Phelip JM, et al. Randomized phase III
trial in elderly patients comparing LV5FU2 with or without irinotecan
for first-line treatment of metastatic colorectal cancer (FFCD
2001-02). Ann Oncol 2016; 27: 121–7.
}
\usage{
FFCD2001_S1C
}
\description{
Kaplan-Meier digitized data from FFCD2001, figure S1C (PMID 26487578). A reported sample size of 212 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(FFCD2001_S1C)

kmplot(FFCD2001_S1C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GETUG16_2}
\alias{GETUG16_2}
\title{GETUG16, figure 2}
\format{
A data frame of 730 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (rt, rt_goserelin) \cr
}
}
\source{
Carrie C, Hasbini A, de Laroche G, et al. Salvage radiotherapy with
or without short-term hormone therapy for rising prostate-specific
antigen concentration after radical prostatectomy (GETUG-AFU 16): a
randomised, multicentre, open-label phase 3 trial. Lancet Oncol 2016;
17: 747–56.
}
\usage{
GETUG16_2
}
\description{
Kaplan-Meier digitized data from GETUG16, figure 2 (PMID 27160475). A reported sample size of 743 for a primary endpoint of PFS in prostate cancer.
}
\examples{
summary(GETUG16_2)

kmplot(GETUG16_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AIO0207(b)_3B}
\alias{AIO0207(b)_3B}
\title{AIO0207(b), figure 3B}
\format{
A data frame of 300 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (5fu_bev, bev) \cr
}
}
\source{
Hegewisch-Becker S, Graeven U, Lerchenmüller CA, et al. Maintenance
strategies after first-line oxaliplatin plus fluoropyrimidine plus
bevacizumab for patients with metastatic colorectal cancer (AIO
0207): a randomised, non-inferiority, open-label, phase 3 trial.
Lancet Oncol 2015; 16: 1355–69.
}
\usage{
`AIO0207(b)_3B`
}
\description{
Kaplan-Meier digitized data from AIO0207(b), figure 3B (PMID 26361971). A reported sample size of 472 for a primary endpoint of time_tx_failure in colorectal cancer.
}
\examples{
summary(`AIO0207(b)_3B`)

kmplot(`AIO0207(b)_3B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RADIANT4_2A}
\alias{RADIANT4_2A}
\title{RADIANT4, figure 2A}
\format{
A data frame of 302 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (everolimus, placebo) \cr
}
}
\source{
Yao JC, Fazio N, Singh S, et al. Everolimus for the treatment of
advanced, non-functional neuroendocrine tumours of the lung or
gastrointestinal tract (RADIANT-4): a randomised, placebo-controlled,
phase 3 study. Lancet 2016; 387: 968–77.
}
\usage{
RADIANT4_2A
}
\description{
Kaplan-Meier digitized data from RADIANT4, figure 2A (PMID 26703889). A reported sample size of 302 for a primary endpoint of PFS in lung/colorectal cancer.
}
\examples{
summary(RADIANT4_2A)

kmplot(RADIANT4_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RADIANT_2A}
\alias{RADIANT_2A}
\title{RADIANT, figure 2A}
\format{
A data frame of 973 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (erlotinib, placebo) \cr
}
}
\source{
Kelly K, Altorki NK, Eberhardt WEE, et al. Adjuvant Erlotinib Versus
Placebo in Patients With Stage IB-IIIA Non-Small-Cell Lung Cancer
(RADIANT): A Randomized, Double-Blind, Phase III Trial. J Clin Oncol
2015; 33: 4007–14.
}
\usage{
RADIANT_2A
}
\description{
Kaplan-Meier digitized data from RADIANT, figure 2A (PMID 26324372). A reported sample size of 973 for a primary endpoint of DFS in lung cancer.
}
\examples{
summary(RADIANT_2A)

kmplot(RADIANT_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SAKK4106_3}
\alias{SAKK4106_3}
\title{SAKK4106, figure 3}
\format{
A data frame of 262 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev, no_bev) \cr
}
}
\source{
Koeberle D, Betticher DC, von Moos R, et al. Bevacizumab continuation
versus no continuation after first-line chemotherapy plus bevacizumab
in patients with metastatic colorectal cancer: a randomized phase III
non-inferiority trial (SAKK 41/06). Ann Oncol 2015; 26: 709–14.
}
\usage{
SAKK4106_3
}
\description{
Kaplan-Meier digitized data from SAKK4106, figure 3 (PMID 25605741). A reported sample size of 262 for a primary endpoint of TTP in colorectal cancer.
}
\examples{
summary(SAKK4106_3)

kmplot(SAKK4106_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ELM-PC4_3A}
\alias{ELM-PC4_3A}
\title{ELM-PC4, figure 3A}
\format{
A data frame of 1,560 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (orteronel, placebo) \cr
}
}
\source{
Saad F, Fizazi K, Jinga V, et al. Orteronel plus prednisone in
patients with chemotherapy-naive metastatic castration-resistant
prostate cancer (ELM-PC 4): a double-blind, multicentre, phase 3,
randomised, placebo-controlled trial. Lancet Oncol 2015; 16: 338–48.
}
\usage{
`ELM-PC4_3A`
}
\description{
Kaplan-Meier digitized data from ELM-PC4, figure 3A (PMID 25701170). A reported sample size of 1,560 for a primary endpoint of rPFS/OS in prostate cancer.
}
\examples{
summary(`ELM-PC4_3A`)

kmplot(`ELM-PC4_3A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CONCUR_3A}
\alias{CONCUR_3A}
\title{CONCUR, figure 3A}
\format{
A data frame of 204 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, regorafenib) \cr
}
}
\source{
Li J, Qin S, Xu R, et al. Regorafenib plus best supportive care
versus placebo plus best supportive care in Asian patients with
previously treated metastatic colorectal cancer (CONCUR): a
randomised, double-blind, placebo-controlled, phase 3 trial. Lancet
Oncol 2015; 16: 619–29.
}
\usage{
CONCUR_3A
}
\description{
Kaplan-Meier digitized data from CONCUR, figure 3A (PMID 25981818). A reported sample size of 243 for a primary endpoint of OS in colorectal cancer.
}
\examples{
summary(CONCUR_3A)

kmplot(CONCUR_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{WSGAGO_2A}
\alias{WSGAGO_2A}
\title{WSGAGO, figure 2A}
\format{
A data frame of 1,930 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab EFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, ed_doc) \cr
}
}
\source{
Nitz U, Gluz O, Huober J, et al. Final analysis of the prospective
WSG-AGO EC-Doc versus FEC phase III trial in intermediate-risk (pN1)
early breast cancer: efficacy and predictive value of Ki67
expression. Ann Oncol 2014; 25: 1551–7.
}
\usage{
WSGAGO_2A
}
\description{
Kaplan-Meier digitized data from WSGAGO, figure 2A (PMID 24827128). A reported sample size of 2,011 for a primary endpoint of EFS in breast cancer.
}
\examples{
summary(WSGAGO_2A)

kmplot(WSGAGO_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmdata-package.R
\docType{package}
\name{kmdata-package}
\alias{kmdata-package}
\alias{kmdata}
\title{kmdata}
\description{
Re-constructed Kaplan-Meier data from publications gathered on
\href{https://pubmed.ncbi.nlm.nih.gov/}{PubMed}.
}
\details{
Data sets were gathered from PubMed and limited to phase 3, randomized
trials completed between January 1st 2014 and December 31st 2016 and
published by JAMA Oncology, New England Journal of Medicine, Lancet,
Lancet Oncology, Journal of Clinical Oncology, Annals of Oncology, Journal
of the American Medical Association, or Journal of the National Cancer
Institute.

The final list includs 263 figures from 152 publications. The data spans
four cancers: colorectal, lung, prostate, and breast. All data sets feature
time and event indicators as well as treatment arm (two arms per data set).

Additional information for each study is included in attributes of the data
objects or in the key data, \code{kmdata_key}. The key includes publication
identifiers, journal, title, outcome, sample size, etc.

This package also includes utilities to generate individual patient data
(IPD) from digitized survival curve data using the method described by
Guyot (2012); see \code{\link{ipd}}.
}
\examples{
## all data sets included
data(package = 'kmdata')


## basic usage
summary(ATTENTION_2A)
kmplot(ATTENTION_2B)


## list of data sets estimating PFS
select_kmdata(Outcome \%in\% 'PFS')


## list of studies in breast cancer with fewer than 200 patients
l <- select_kmdata(Cancer \%in\% 'Breast' & ReportedSampleSize < 200, return = 'data')
par(mfrow = n2mfrow(length(l)))
sapply(l, kmplot)

}
\seealso{
\code{\link{kmdata_key}} for a list of the data sets available with some
additional information about each including data source, outcomes, disease
type, results, and quality of each re-capitulated data set.

\code{\link[=summary.kmdata]{summary}} for a method to summarize data sets

\code{\link{kmplot}} to plot Kaplan-Meier curves for each data set

\code{\link{select_kmdata}} to select studies based on characteristics
listed in the \code{\link{kmdata_key}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CALGB40503_2B}
\alias{CALGB40503_2B}
\title{CALGB40503, figure 2B}
\format{
A data frame of 343 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l, l_b) \cr
}
}
\source{
Dickler MN, Barry WT, Cirrincione CT, et al. Phase III Trial
Evaluating Letrozole As First-Line Endocrine Therapy With or Without
Bevacizumab for the Treatment of Postmenopausal Women With Hormone
Receptor-Positive Advanced-Stage Breast Cancer: CALGB 40503
(Alliance). J Clin Oncol 2016; 34: 2602–9.
}
\usage{
CALGB40503_2B
}
\description{
Kaplan-Meier digitized data from CALGB40503, figure 2B (PMID 27138575). A reported sample size of 350 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CALGB40503_2B)

kmplot(CALGB40503_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ATTENTION_2A}
\alias{ATTENTION_2A}
\title{ATTENTION, figure 2A}
\format{
A data frame of 307 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, tivantinib) \cr
}
}
\source{
Yoshioka H, Azuma K, Yamamoto N, et al. A randomized, double-blind,
placebo-controlled, phase III trial of erlotinib with or without a
c-Met inhibitor tivantinib (ARQ 197) in Asian patients with
previously treated stage IIIB/IV nonsquamous nonsmall-cell lung
cancer harboring wild-type epidermal growth factor receptor
(ATTENTION study). Ann Oncol 2015; 26: 2066–72.
}
\usage{
ATTENTION_2A
}
\description{
Kaplan-Meier digitized data from ATTENTION, figure 2A (PMID 26153496). A reported sample size of 460 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ATTENTION_2A)

kmplot(ATTENTION_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ALTTO(b)_2B}
\alias{ALTTO(b)_2B}
\title{ALTTO(b), figure 2B}
\format{
A data frame of 4,188 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (t, t_then_l) \cr
}
}
\source{
Piccart-Gebhart M, Holmes E, Baselga J, et al. Adjuvant Lapatinib and
Trastuzumab for Early Human Epidermal Growth Factor Receptor
2-Positive Breast Cancer: Results From the Randomized Phase III
Adjuvant Lapatinib and/or Trastuzumab Treatment Optimization Trial. J
Clin Oncol 2016; 34: 1034–42.
}
\usage{
`ALTTO(b)_2B`
}
\description{
Kaplan-Meier digitized data from ALTTO(b), figure 2B (PMID 26598744). A reported sample size of 8,381 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`ALTTO(b)_2B`)

kmplot(`ALTTO(b)_2B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RTOG0415_2}
\alias{RTOG0415_2}
\title{RTOG0415, figure 2}
\format{
A data frame of 1,092 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (70gy, 73.8gy) \cr
}
}
\source{
Lee WR, Dignam JJ, Amin MB, et al. Randomized Phase III
Noninferiority Study Comparing Two Radiotherapy Fractionation
Schedules in Patients With Low-Risk Prostate Cancer. J Clin Oncol
2016; 34: 2325–32.
}
\usage{
RTOG0415_2
}
\description{
Kaplan-Meier digitized data from RTOG0415, figure 2 (PMID 27044935). A reported sample size of 1,115 for a primary endpoint of DFS in prostate cancer.
}
\examples{
summary(RTOG0415_2)

kmplot(RTOG0415_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TH3RESA_2A}
\alias{TH3RESA_2A}
\title{TH3RESA, figure 2A}
\format{
A data frame of 602 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (physicians_choice, trastuzumab_emtansine) \cr
}
}
\source{
Krop IE, Kim S-B, González-Martín A, et al. Trastuzumab emtansine
versus treatment of physician’s choice for pretreated HER2-positive
advanced breast cancer (TH3RESA): a randomised, open-label, phase 3
trial. Lancet Oncol 2014; 15: 689–99.
}
\usage{
TH3RESA_2A
}
\description{
Kaplan-Meier digitized data from TH3RESA, figure 2A (PMID 24793816). A reported sample size of 602 for a primary endpoint of PFS/OS in breast cancer.
}
\examples{
summary(TH3RESA_2A)

kmplot(TH3RESA_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NO16968_2A}
\alias{NO16968_2A}
\title{NO16968, figure 2A}
\format{
A data frame of 1,886 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (fu_fa, xelox) \cr
}
}
\source{
Schmoll H-J, Tabernero J, Maroun J, et al. Capecitabine Plus
Oxaliplatin Compared With Fluorouracil/Folinic Acid As Adjuvant
Therapy for Stage III Colon Cancer: Final Results of the NO16968
Randomized Controlled Phase III Trial. J Clin Oncol 2015; 33:
3733–40.
}
\usage{
NO16968_2A
}
\description{
Kaplan-Meier digitized data from NO16968, figure 2A (PMID 26324362). A reported sample size of 1,886 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(NO16968_2A)

kmplot(NO16968_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{OPTIMOX3_3A}
\alias{OPTIMOX3_3A}
\title{OPTIMOX3, figure 3A}
\format{
A data frame of 452 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev, bev_erlotinib) \cr
}
}
\source{
Tournigand C, Chibaudel B, Samson B, et al. Bevacizumab with or
without erlotinib as maintenance therapy in patients with metastatic
colorectal cancer (GERCOR DREAM; OPTIMOX3): a randomised, open-label,
phase 3 trial. Lancet Oncol 2015; 16: 1493–505.
}
\usage{
OPTIMOX3_3A
}
\description{
Kaplan-Meier digitized data from OPTIMOX3, figure 3A (PMID 26474518). A reported sample size of 700 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(OPTIMOX3_3A)

kmplot(OPTIMOX3_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TRAPEZE_1B}
\alias{TRAPEZE_1B}
\title{TRAPEZE, figure 1B}
\format{
A data frame of 757 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (no_sr89, sr89) \cr
}
}
\source{
James ND, Pirrie SJ, Pope AM, et al. Clinical Outcomes and Survival
Following Treatment of Metastatic Castrate-Refractory Prostate Cancer
With Docetaxel Alone or With Strontium-89, Zoledronic Acid, or Both:
The TRAPEZE Randomized Clinical Trial. JAMA Oncol 2016; 2: 493–9.
}
\usage{
TRAPEZE_1B
}
\description{
Kaplan-Meier digitized data from TRAPEZE, figure 1B (PMID 26794729). A reported sample size of 757 for a primary endpoint of cPFS,cost-effectiveness in prostate cancer.
}
\examples{
summary(TRAPEZE_1B)

kmplot(TRAPEZE_1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PROFILE1014_1B}
\alias{PROFILE1014_1B}
\title{PROFILE1014, figure 1B}
\format{
A data frame of 343 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, crizotinib) \cr
}
}
\source{
Solomon BJ, Mok T, Kim D-W, et al. First-line crizotinib versus
chemotherapy in ALK-positive lung cancer. N Engl J Med 2014; 371:
2167–77.
}
\usage{
PROFILE1014_1B
}
\description{
Kaplan-Meier digitized data from PROFILE1014, figure 1B (PMID 25470694). A reported sample size of 343 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(PROFILE1014_1B)

kmplot(PROFILE1014_1B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NSABPB40_3A}
\alias{NSABPB40_3A}
\title{NSABPB40, figure 3A}
\format{
A data frame of 1,184 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev, no_bev) \cr
}
}
\source{
Bear HD, Tang G, Rastogi P, et al. Neoadjuvant plus adjuvant
bevacizumab in early breast cancer (NSABP B-40 [NRG Oncology]):
secondary outcomes of a phase 3, randomised controlled trial. Lancet
Oncol 2015; 16: 1037–48.
}
\usage{
NSABPB40_3A
}
\description{
Kaplan-Meier digitized data from NSABPB40, figure 3A (PMID 26272770). A reported sample size of 1,206 for a primary endpoint of pCR in breast cancer.
}
\examples{
summary(NSABPB40_3A)

kmplot(NSABPB40_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CLEOPATRA_2A}
\alias{CLEOPATRA_2A}
\title{CLEOPATRA, figure 2A}
\format{
A data frame of 808 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, pertuzumab) \cr
}
}
\source{
Swain SM, Baselga J, Kim S-B, et al. Pertuzumab, trastuzumab, and
docetaxel in HER2-positive metastatic breast cancer. N Engl J Med
2015; 372: 724–34.
}
\usage{
CLEOPATRA_2A
}
\description{
Kaplan-Meier digitized data from CLEOPATRA, figure 2A (PMID 25693012). A reported sample size of 808 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CLEOPATRA_2A)

kmplot(CLEOPATRA_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BEYOND_2A}
\alias{BEYOND_2A}
\title{BEYOND, figure 2A}
\format{
A data frame of 276 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev_cp, placebo_cp) \cr
}
}
\source{
Zhou C, Wu Y-L, Chen G, et al. BEYOND: A Randomized, Double-Blind,
Placebo-Controlled, Multicenter, Phase III Study of First-Line
Carboplatin/Paclitaxel Plus Bevacizumab or Placebo in Chinese
Patients With Advanced or Recurrent Nonsquamous Non-Small-Cell Lung
Cancer. J Clin Oncol 2015; 33: 2197–204.
}
\usage{
BEYOND_2A
}
\description{
Kaplan-Meier digitized data from BEYOND, figure 2A (PMID 26014294). A reported sample size of 276 for a primary endpoint of PFS in prostate cancer.
}
\examples{
summary(BEYOND_2A)

kmplot(BEYOND_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CBCSG006_2B}
\alias{CBCSG006_2B}
\title{CBCSG006, figure 2B}
\format{
A data frame of 236 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cis_gemcitabine, paclitaxel_gemcitabine) \cr
}
}
\source{
Hu X-C, Zhang J, Xu B-H, et al. Cisplatin plus gemcitabine versus
paclitaxel plus gemcitabine as first-line therapy for metastatic
triple-negative breast cancer (CBCSG006): a randomised, open-label,
multicentre, phase 3 trial. Lancet Oncol 2015; 16: 436–46.
}
\usage{
CBCSG006_2B
}
\description{
Kaplan-Meier digitized data from CBCSG006, figure 2B (PMID 25795409). A reported sample size of 240 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CBCSG006_2B)

kmplot(CBCSG006_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{DELTA_2A}
\alias{DELTA_2A}
\title{DELTA, figure 2A}
\format{
A data frame of 301 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d1, erlotinib) \cr
}
}
\source{
Kawaguchi T, Ando M, Asami K, et al. Randomized phase III trial of
erlotinib versus docetaxel as second- or third-line therapy in
patients with advanced non-small-cell lung cancer: Docetaxel and
Erlotinib Lung Cancer Trial (DELTA). J Clin Oncol 2014; 32: 1902–8.
}
\usage{
DELTA_2A
}
\description{
Kaplan-Meier digitized data from DELTA, figure 2A (PMID 24841974). A reported sample size of 301 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(DELTA_2A)

kmplot(DELTA_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TH3RESA_4}
\alias{TH3RESA_4}
\title{TH3RESA, figure 4}
\format{
A data frame of 602 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (physicians_choice, trastuzumab_emtansine) \cr
}
}
\source{
Krop IE, Kim S-B, González-Martín A, et al. Trastuzumab emtansine
versus treatment of physician’s choice for pretreated HER2-positive
advanced breast cancer (TH3RESA): a randomised, open-label, phase 3
trial. Lancet Oncol 2014; 15: 689–99.
}
\usage{
TH3RESA_4
}
\description{
Kaplan-Meier digitized data from TH3RESA, figure 4 (PMID 24793816). A reported sample size of 602 for a primary endpoint of PFS/OS in breast cancer.
}
\examples{
summary(TH3RESA_4)

kmplot(TH3RESA_4)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RTOG0617_2A}
\alias{RTOG0617_2A}
\title{RTOG0617, figure 2A}
\format{
A data frame of 424 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (60gy, 74gy) \cr
}
}
\source{
Bradley JD, Paulus R, Komaki R, et al. Standard-dose versus high-dose
conformal radiotherapy with concurrent and consolidation carboplatin
plus paclitaxel with or without cetuximab for patients with stage
IIIA or IIIB non-small-cell lung cancer (RTOG 0617): a randomised,
two-by-two factorial phase 3 study. Lancet Oncol 2015; 16: 187–99.
}
\usage{
RTOG0617_2A
}
\description{
Kaplan-Meier digitized data from RTOG0617, figure 2A (PMID 25601342). A reported sample size of 166 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(RTOG0617_2A)

kmplot(RTOG0617_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{WJOG5108L_2A}
\alias{WJOG5108L_2A}
\title{WJOG5108L, figure 2A}
\format{
A data frame of 559 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (erlotinib, gefitinib) \cr
}
}
\source{
Urata Y, Katakami N, Morita S, et al. Randomized Phase III Study
Comparing Gefitinib With Erlotinib in Patients With Previously
Treated Advanced Lung Adenocarcinoma: WJOG 5108L. J Clin Oncol 2016;
34: 3248–57.
}
\usage{
WJOG5108L_2A
}
\description{
Kaplan-Meier digitized data from WJOG5108L, figure 2A (PMID 27022112). A reported sample size of 561 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(WJOG5108L_2A)

kmplot(WJOG5108L_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{summary.kmdata}
\alias{summary.kmdata}
\title{Summarize \code{kmdata}}
\usage{
\method{summary}{kmdata}(object, message = TRUE, ...)
}
\arguments{
\item{object}{an object of class \code{"kmdata"}}

\item{message}{logical; if \code{TRUE} (default), a description of the
study will be printed with the summary}

\item{...}{ignored}
}
\description{
Print a summary of a \code{kmdata} object.
}
\examples{
summary(ATTENTION_2A)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ITACa_2B}
\alias{ITACa_2B}
\title{ITACa, figure 2B}
\format{
A data frame of 370 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (ct, ct_b) \cr
}
}
\source{
Passardi A, Nanni O, Tassinari D, et al. Effectiveness of bevacizumab
added to standard chemotherapy in metastatic colorectal cancer: final
results for first-line treatment from the ITACa randomized clinical
trial. Ann Oncol 2015; 26: 1201–7.
}
\usage{
ITACa_2B
}
\description{
Kaplan-Meier digitized data from ITACa, figure 2B (PMID 25735317). A reported sample size of 376 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(ITACa_2B)

kmplot(ITACa_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00614393(b)_3B}
\alias{NCT00614393(b)_3B}
\title{NCT00614393(b), figure 3B}
\format{
A data frame of 228 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (armb, armc) \cr
}
}
\source{
Sclafani F, Kim TY, Cunningham D, et al. A Randomized Phase II/III
Study of Dalotuzumab in Combination With Cetuximab and Irinotecan in
Chemorefractory, KRAS Wild-Type, Metastatic Colorectal Cancer. J Natl
Cancer Inst 2015; 107: djv258.
}
\usage{
`NCT00614393(b)_3B`
}
\description{
Kaplan-Meier digitized data from NCT00614393(b), figure 3B (PMID 26405092). A reported sample size of 344 for a primary endpoint of PFS,OS in colorectal cancer.
}
\examples{
summary(`NCT00614393(b)_3B`)

kmplot(`NCT00614393(b)_3B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GECESTRO-APBI_4}
\alias{GECESTRO-APBI_4}
\title{GECESTRO-APBI, figure 4}
\format{
A data frame of 1,184 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (apbi, wbi) \cr
}
}
\source{
Strnad V, Ott OJ, Hildebrandt G, et al. 5-year results of accelerated
partial breast irradiation using sole interstitial multicatheter
brachytherapy versus whole-breast irradiation with boost after
breast-conserving surgery for low-risk invasive and in-situ carcinoma
of the female breast: a randomised, phase 3, non-inferiority trial.
Lancet 2016; 387: 229–38.
}
\usage{
`GECESTRO-APBI_4`
}
\description{
Kaplan-Meier digitized data from GECESTRO-APBI, figure 4 (PMID 26494415). A reported sample size of 551 for a primary endpoint of LR in breast cancer.
}
\examples{
summary(`GECESTRO-APBI_4`)

kmplot(`GECESTRO-APBI_4`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NTR1527_2}
\alias{NTR1527_2}
\title{NTR1527, figure 2}
\format{
A data frame of 495 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, thoracic rt) \cr
}
}
\source{
Slotman BJ, van Tinteren H, Praag JO, et al. Use of thoracic
radiotherapy for extensive stage small-cell lung cancer: a phase 3
randomised controlled trial. Lancet 2015; 385: 36–42.
}
\usage{
NTR1527_2
}
\description{
Kaplan-Meier digitized data from NTR1527, figure 2 (PMID 25230595). A reported sample size of 498 for a primary endpoint of OS-1yr in lung cancer.
}
\examples{
summary(NTR1527_2)

kmplot(NTR1527_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CBCSG006_2A}
\alias{CBCSG006_2A}
\title{CBCSG006, figure 2A}
\format{
A data frame of 236 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cis_gemcitabine, paclitaxel_gemcitabine) \cr
}
}
\source{
Hu X-C, Zhang J, Xu B-H, et al. Cisplatin plus gemcitabine versus
paclitaxel plus gemcitabine as first-line therapy for metastatic
triple-negative breast cancer (CBCSG006): a randomised, open-label,
multicentre, phase 3 trial. Lancet Oncol 2015; 16: 436–46.
}
\usage{
CBCSG006_2A
}
\description{
Kaplan-Meier digitized data from CBCSG006, figure 2A (PMID 25795409). A reported sample size of 240 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CBCSG006_2A)

kmplot(CBCSG006_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{VANTAGE014_4}
\alias{VANTAGE014_4}
\title{VANTAGE014, figure 4}
\format{
A data frame of 661 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in weeks) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, vorinostat) \cr
}
}
\source{
Krug LM, Kindler HL, Calvert H, et al. Vorinostat in patients with
advanced malignant pleural mesothelioma who have progressed on
previous chemotherapy (VANTAGE-014): a phase 3, double-blind,
randomised, placebo-controlled trial. Lancet Oncol 2015; 16: 447–56.
}
\usage{
VANTAGE014_4
}
\description{
Kaplan-Meier digitized data from VANTAGE014, figure 4 (PMID 25800891). A reported sample size of 661 for a primary endpoint of OS/safety in lung cancer.
}
\examples{
summary(VANTAGE014_4)

kmplot(VANTAGE014_4)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ALTTO(a)_2A}
\alias{ALTTO(a)_2A}
\title{ALTTO(a), figure 2A}
\format{
A data frame of 4,190 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l_t, t) \cr
}
}
\source{
Piccart-Gebhart M, Holmes E, Baselga J, et al. Adjuvant Lapatinib and
Trastuzumab for Early Human Epidermal Growth Factor Receptor
2-Positive Breast Cancer: Results From the Randomized Phase III
Adjuvant Lapatinib and/or Trastuzumab Treatment Optimization Trial. J
Clin Oncol 2016; 34: 1034–42.
}
\usage{
`ALTTO(a)_2A`
}
\description{
Kaplan-Meier digitized data from ALTTO(a), figure 2A (PMID 26598744). A reported sample size of 8,381 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`ALTTO(a)_2A`)

kmplot(`ALTTO(a)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NOAH_2A}
\alias{NOAH_2A}
\title{NOAH, figure 2A}
\format{
A data frame of 235 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab EFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, chemotherapy_trastuzumab) \cr
}
}
\source{
Gianni L, Eiermann W, Semiglazov V, et al. Neoadjuvant and adjuvant
trastuzumab in patients with HER2-positive locally advanced breast
cancer (NOAH): follow-up of a randomised controlled superiority trial
with a parallel HER2-negative cohort. Lancet Oncol 2014; 15: 640–7.
}
\usage{
NOAH_2A
}
\description{
Kaplan-Meier digitized data from NOAH, figure 2A (PMID 24657003). A reported sample size of 235 for a primary endpoint of EFS in breast cancer.
}
\examples{
summary(NOAH_2A)

kmplot(NOAH_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT01234311_2A}
\alias{NCT01234311_2A}
\title{NCT01234311, figure 2A}
\format{
A data frame of 577 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, tasq) \cr
}
}
\source{
Sternberg C, Armstrong A, Pili R, et al. Randomized, Double-Blind,
Placebo-Controlled Phase III Study of Tasquinimod in Men With
Metastatic Castration-Resistant Prostate Cancer. J Clin Oncol 2016;
34: 2636–43.
}
\usage{
NCT01234311_2A
}
\description{
Kaplan-Meier digitized data from NCT01234311, figure 2A (PMID 27298414). A reported sample size of 1,245 for a primary endpoint of rPFS in prostate cancer.
}
\examples{
summary(NCT01234311_2A)

kmplot(NCT01234311_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CLEOPATRA_3A}
\alias{CLEOPATRA_3A}
\title{CLEOPATRA, figure 3A}
\format{
A data frame of 808 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, pertuzumab) \cr
}
}
\source{
Swain SM, Baselga J, Kim S-B, et al. Pertuzumab, trastuzumab, and
docetaxel in HER2-positive metastatic breast cancer. N Engl J Med
2015; 372: 724–34.
}
\usage{
CLEOPATRA_3A
}
\description{
Kaplan-Meier digitized data from CLEOPATRA, figure 3A (PMID 25693012). A reported sample size of 808 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CLEOPATRA_3A)

kmplot(CLEOPATRA_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RECOURSE_1A}
\alias{RECOURSE_1A}
\title{RECOURSE, figure 1A}
\format{
A data frame of 800 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, tas-102) \cr
}
}
\source{
Mayer RJ, Van Cutsem E, Falcone A, et al. Randomized trial of TAS-102
for refractory metastatic colorectal cancer. N Engl J Med 2015; 372:
1909–19.
}
\usage{
RECOURSE_1A
}
\description{
Kaplan-Meier digitized data from RECOURSE, figure 1A (PMID 25970050). A reported sample size of 800 for a primary endpoint of OS in breast cancer.
}
\examples{
summary(RECOURSE_1A)

kmplot(RECOURSE_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{OMEGA_3B}
\alias{OMEGA_3B}
\title{OMEGA, figure 3B}
\format{
A data frame of 78 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cap, peg_doxorubicin) \cr
}
}
\source{
Smorenburg CH, de Groot SM, van Leeuwen-Stok AE, et al. A randomized
phase III study comparing pegylated liposomal doxorubicin with
capecitabine as first-line chemotherapy in elderly patients with
metastatic breast cancer: results of the OMEGA study of the Dutch
Breast Cancer Research Group BOOG. Ann Oncol 2014; 25: 599–605.
}
\usage{
OMEGA_3B
}
\description{
Kaplan-Meier digitized data from OMEGA, figure 3B (PMID 24504445). A reported sample size of 78 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(OMEGA_3B)

kmplot(OMEGA_3B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Checkmate017_2B}
\alias{Checkmate017_2B}
\title{Checkmate017, figure 2B}
\format{
A data frame of 272 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d1, nivolumab) \cr
}
}
\source{
Brahmer J, Reckamp KL, Baas P, et al. Nivolumab versus Docetaxel in
Advanced Squamous-Cell Non-Small-Cell Lung Cancer. N Engl J Med 2015;
373: 123–35.
}
\usage{
Checkmate017_2B
}
\description{
Kaplan-Meier digitized data from Checkmate017, figure 2B (PMID 26028407). A reported sample size of 272 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(Checkmate017_2B)

kmplot(Checkmate017_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ELM-PC5_2A}
\alias{ELM-PC5_2A}
\title{ELM-PC5, figure 2A}
\format{
A data frame of 1,099 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (orteronel, placebo) \cr
}
}
\source{
Fizazi K, Jones R, Oudard S, et al. Phase III, randomized,
double-blind, multicenter trial comparing orteronel (TAK-700) plus
prednisone with placebo plus prednisone in patients with metastatic
castration-resistant prostate cancer that has progressed during or
after docetaxel-based therapy: ELM-PC 5. J Clin Oncol 2015; 33:
723–31.
}
\usage{
`ELM-PC5_2A`
}
\description{
Kaplan-Meier digitized data from ELM-PC5, figure 2A (PMID 25624429). A reported sample size of 1,099 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(`ELM-PC5_2A`)

kmplot(`ELM-PC5_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CALGB40503_2A}
\alias{CALGB40503_2A}
\title{CALGB40503, figure 2A}
\format{
A data frame of 343 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (l, l_b) \cr
}
}
\source{
Dickler MN, Barry WT, Cirrincione CT, et al. Phase III Trial
Evaluating Letrozole As First-Line Endocrine Therapy With or Without
Bevacizumab for the Treatment of Postmenopausal Women With Hormone
Receptor-Positive Advanced-Stage Breast Cancer: CALGB 40503
(Alliance). J Clin Oncol 2016; 34: 2602–9.
}
\usage{
CALGB40503_2A
}
\description{
Kaplan-Meier digitized data from CALGB40503, figure 2A (PMID 27138575). A reported sample size of 350 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(CALGB40503_2A)

kmplot(CALGB40503_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{EORTC22881_2}
\alias{EORTC22881_2}
\title{EORTC22881, figure 2}
\format{
A data frame of 5,318 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (boost, no boost) \cr
}
}
\source{
Bartelink H, Maingon P, Poortmans P, et al. Whole-breast irradiation
with or without a boost for patients treated with breast-conserving
surgery for early breast cancer: 20-year follow-up of a randomised
phase 3 trial. Lancet Oncol 2015; 16: 47–56.
}
\usage{
EORTC22881_2
}
\description{
Kaplan-Meier digitized data from EORTC22881, figure 2 (PMID 25500422). A reported sample size of 2,657 for a primary endpoint of OS in breast cancer.
}
\examples{
summary(EORTC22881_2)

kmplot(EORTC22881_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{JCOG0605_2A}
\alias{JCOG0605_2A}
\title{JCOG0605, figure 2A}
\format{
A data frame of 180 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (combination_chemotherapy, topotecan) \cr
}
}
\source{
Goto K, Ohe Y, Shibata T, et al. Combined chemotherapy with
cisplatin, etoposide, and irinotecan versus topotecan alone as
second-line treatment for patients with sensitive relapsed small-cell
lung cancer (JCOG0605): a multicentre, open-label, randomised phase 3
trial. Lancet Oncol 2016; 17: 1147–57.
}
\usage{
JCOG0605_2A
}
\description{
Kaplan-Meier digitized data from JCOG0605, figure 2A (PMID 27312053). A reported sample size of 180 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(JCOG0605_2A)

kmplot(JCOG0605_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CEREBEL_2A}
\alias{CEREBEL_2A}
\title{CEREBEL, figure 2A}
\format{
A data frame of 540 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (lapatinib_cape, trastuzumab_cape) \cr
}
}
\source{
Pivot X, Manikhas A, Żurawski B, et al. CEREBEL (EGF111438): A Phase
III, Randomized, Open-Label Study of Lapatinib Plus Capecitabine
Versus Trastuzumab Plus Capecitabine in Patients With Human Epidermal
Growth Factor Receptor 2-Positive Metastatic Breast Cancer. J Clin
Oncol 2015; 33: 1564–73.
}
\usage{
CEREBEL_2A
}
\description{
Kaplan-Meier digitized data from CEREBEL, figure 2A (PMID 25605838). A reported sample size of 540 for a primary endpoint of Incidence_CNSmets in breast cancer.
}
\examples{
summary(CEREBEL_2A)

kmplot(CEREBEL_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RECOURSE_2A}
\alias{RECOURSE_2A}
\title{RECOURSE, figure 2A}
\format{
A data frame of 800 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, tas102) \cr
}
}
\source{
Mayer RJ, Van Cutsem E, Falcone A, et al. Randomized trial of TAS-102
for refractory metastatic colorectal cancer. N Engl J Med 2015; 372:
1909–19.
}
\usage{
RECOURSE_2A
}
\description{
Kaplan-Meier digitized data from RECOURSE, figure 2A (PMID 25970050). A reported sample size of 800 for a primary endpoint of OS in breast cancer.
}
\examples{
summary(RECOURSE_2A)

kmplot(RECOURSE_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SELECTBC_2B}
\alias{SELECTBC_2B}
\title{SELECTBC, figure 2B}
\format{
A data frame of 592 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab TTF event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (s1, taxane) \cr
}
}
\source{
Takashima T, Mukai H, Hara F, et al. Taxanes versus S-1 as the
first-line chemotherapy for metastatic breast cancer (SELECT BC): an
open-label, non-inferiority, randomised phase 3 trial. Lancet Oncol
2016; 17: 90–8.
}
\usage{
SELECTBC_2B
}
\description{
Kaplan-Meier digitized data from SELECTBC, figure 2B (PMID 26617202). A reported sample size of 618 for a primary endpoint of OS in breast cancer.
}
\examples{
summary(SELECTBC_2B)

kmplot(SELECTBC_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PALOMA3_2A}
\alias{PALOMA3_2A}
\title{PALOMA3, figure 2A}
\format{
A data frame of 521 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (fulvestrant_palbociclib, fulvestrant_placebo) \cr
}
}
\source{
Cristofanilli M, Turner NC, Bondarenko I, et al. Fulvestrant plus
palbociclib versus fulvestrant plus placebo for treatment of
hormone-receptor-positive, HER2-negative metastatic breast cancer
that progressed on previous endocrine therapy (PALOMA-3): final
analysis of the multicentre, double-blind, phase 3 randomised
controlled trial. Lancet Oncol 2016; 17: 425–39.
}
\usage{
PALOMA3_2A
}
\description{
Kaplan-Meier digitized data from PALOMA3, figure 2A (PMID 26947331). A reported sample size of 521 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(PALOMA3_2A)

kmplot(PALOMA3_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TROG0306_2}
\alias{TROG0306_2}
\title{TROG0306, figure 2}
\format{
A data frame of 293 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (delayed_adt, immediate_adt) \cr
}
}
\source{
Duchesne GM, Woo HH, Bassett JK, et al. Timing of
androgen-deprivation therapy in patients with prostate cancer with a
rising PSA (TROG 03.06 and VCOG PR 01-03 [TOAD]): a randomised,
multicentre, non-blinded, phase 3 trial. Lancet Oncol 2016; 17:
727–37.
}
\usage{
TROG0306_2
}
\description{
Kaplan-Meier digitized data from TROG0306, figure 2 (PMID 27155740). A reported sample size of 293 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(TROG0306_2)

kmplot(TROG0306_2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00614393(b)_3A}
\alias{NCT00614393(b)_3A}
\title{NCT00614393(b), figure 3A}
\format{
A data frame of 228 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (armb, armc) \cr
}
}
\source{
Sclafani F, Kim TY, Cunningham D, et al. A Randomized Phase II/III
Study of Dalotuzumab in Combination With Cetuximab and Irinotecan in
Chemorefractory, KRAS Wild-Type, Metastatic Colorectal Cancer. J Natl
Cancer Inst 2015; 107: djv258.
}
\usage{
`NCT00614393(b)_3A`
}
\description{
Kaplan-Meier digitized data from NCT00614393(b), figure 3A (PMID 26405092). A reported sample size of 344 for a primary endpoint of PFS,OS in colorectal cancer.
}
\examples{
summary(`NCT00614393(b)_3A`)

kmplot(`NCT00614393(b)_3A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MAINSAIL_2B}
\alias{MAINSAIL_2B}
\title{MAINSAIL, figure 2B}
\format{
A data frame of 1,059 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (lenolidomide, placebo) \cr
}
}
\source{
Petrylak DP, Vogelzang NJ, Budnik N, et al. Docetaxel and prednisone
with or without lenalidomide in chemotherapy-naive patients with
metastatic castration-resistant prostate cancer (MAINSAIL): a
randomised, double-blind, placebo-controlled phase 3 trial. Lancet
Oncol 2015; 16: 417–25.
}
\usage{
MAINSAIL_2B
}
\description{
Kaplan-Meier digitized data from MAINSAIL, figure 2B (PMID 25743937). A reported sample size of 1,059 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(MAINSAIL_2B)

kmplot(MAINSAIL_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PETACC8_2A}
\alias{PETACC8_2A}
\title{PETACC8, figure 2A}
\format{
A data frame of 1,602 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfox4, folfox4_cetuximab) \cr
}
}
\source{
Taieb J, Tabernero J, Mini E, et al. Oxaliplatin, fluorouracil, and
leucovorin with or without cetuximab in patients with resected stage
III colon cancer (PETACC-8): an open-label, randomised phase 3 trial.
Lancet Oncol 2014; 15: 862–73.
}
\usage{
PETACC8_2A
}
\description{
Kaplan-Meier digitized data from PETACC8, figure 2A (PMID 24928083). A reported sample size of 2,559 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(PETACC8_2A)

kmplot(PETACC8_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BEACON_2B}
\alias{BEACON_2B}
\title{BEACON, figure 2B}
\format{
A data frame of 852 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, etirinotecan_pegol) \cr
}
}
\source{
Perez EA, Awada A, O’Shaughnessy J, et al. Etirinotecan pegol
(NKTR-102) versus treatment of physician’s choice in women with
advanced breast cancer previously treated with an anthracycline, a
taxane, and capecitabine (BEACON): a randomised, open-label,
multicentre, phase 3 trial. Lancet Oncol 2015; 16: 1556–68.
}
\usage{
BEACON_2B
}
\description{
Kaplan-Meier digitized data from BEACON, figure 2B (PMID 26482278). A reported sample size of 852 for a primary endpoint of OS in breast cancer.
}
\examples{
summary(BEACON_2B)

kmplot(BEACON_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ASPECCT_3A}
\alias{ASPECCT_3A}
\title{ASPECCT, figure 3A}
\format{
A data frame of 999 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cetuximab, panitumumab) \cr
}
}
\source{
Price TJ, Peeters M, Kim TW, et al. Panitumumab versus cetuximab in
patients with chemotherapy-refractory wild-type KRAS exon 2
metastatic colorectal cancer (ASPECCT): a randomised, multicentre,
open-label, non-inferiority phase 3 study. Lancet Oncol 2014; 15:
569–79.
}
\usage{
ASPECCT_3A
}
\description{
Kaplan-Meier digitized data from ASPECCT, figure 3A (PMID 24739896). A reported sample size of 1,010 for a primary endpoint of OS in colorectal cancer.
}
\examples{
summary(ASPECCT_3A)

kmplot(ASPECCT_3A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SMX1132531_1}
\alias{SMX1132531_1}
\title{SMX1132531, figure 1}
\format{
A data frame of 98 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (srs, upfront_chemotherapy) \cr
}
}
\source{
Lim SH, Lee JY, Lee M-Y, et al. A randomized phase III trial of
stereotactic radiosurgery (SRS) versus observation for patients with
asymptomatic cerebral oligo-metastases in non-small-cell lung cancer.
Ann Oncol 2015; 26: 762–8.
}
\usage{
SMX1132531_1
}
\description{
Kaplan-Meier digitized data from SMX1132531, figure 1 (PMID 25538174). A reported sample size of 105 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(SMX1132531_1)

kmplot(SMX1132531_1)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NITRO_2B}
\alias{NITRO_2B}
\title{NITRO, figure 2B}
\format{
A data frame of 372 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, nitro_chemo) \cr
}
}
\source{
Davidson A, Veillard A-S, Tognela A, et al. A phase III randomized
trial of adding topical nitroglycerin to first-line chemotherapy for
advanced nonsmall-cell lung cancer: the Australasian lung cancer
trials group NITRO trial. Ann Oncol 2015; 26: 2280–6.
}
\usage{
NITRO_2B
}
\description{
Kaplan-Meier digitized data from NITRO, figure 2B (PMID 26347110). A reported sample size of 372 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(NITRO_2B)

kmplot(NITRO_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MINDACT_2E}
\alias{MINDACT_2E}
\title{MINDACT, figure 2E}
\format{
A data frame of 690 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, no_chemotherapy) \cr
}
}
\source{
Cardoso F, van’t Veer LJ, Bogaerts J, et al. 70-Gene Signature as an
Aid to Treatment Decisions in Early-Stage Breast Cancer. N Engl J Med
2016; 375: 717–29.
}
\usage{
MINDACT_2E
}
\description{
Kaplan-Meier digitized data from MINDACT, figure 2E (PMID 27557300). A reported sample size of 1,550 for a primary endpoint of DMFS in breast cancer.
}
\examples{
summary(MINDACT_2E)

kmplot(MINDACT_2E)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{JFMC33-0502_2A}
\alias{JFMC33-0502_2A}
\title{JFMC33-0502, figure 2A}
\format{
A data frame of 1,060 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, study) \cr
}
}
\source{
Sadahiro S, Tsuchiya T, Sasaki K, et al. Randomized phase III trial
of treatment duration for oral uracil and tegafur plus leucovorin as
adjuvant chemotherapy for patients with stage IIB/III colon cancer:
final results of JFMC33-0502. Ann Oncol 2015; 26: 2274–80.
}
\usage{
`JFMC33-0502_2A`
}
\description{
Kaplan-Meier digitized data from JFMC33-0502, figure 2A (PMID 26347106). A reported sample size of 1,071 for a primary endpoint of DFS in colorectal cancer.
}
\examples{
summary(`JFMC33-0502_2A`)

kmplot(`JFMC33-0502_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PEGASE07_2A}
\alias{PEGASE07_2A}
\title{PEGASE07, figure 2A}
\format{
A data frame of 174 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (arma, armb) \cr
}
}
\source{
Gonçalves A, Pierga J-Y, Ferrero J-M, et al. UNICANCER-PEGASE 07
study: a randomized phase III trial evaluating postoperative
docetaxel-5FU regimen after neoadjuvant dose-intense chemotherapy for
treatment of inflammatory breast cancer. Ann Oncol 2015; 26: 1692–7.
}
\usage{
PEGASE07_2A
}
\description{
Kaplan-Meier digitized data from PEGASE07, figure 2A (PMID 25943350). A reported sample size of 174 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(PEGASE07_2A)

kmplot(PEGASE07_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MAGRIT_2C}
\alias{MAGRIT_2C}
\title{MAGRIT, figure 2C}
\format{
A data frame of 2,272 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (mage_a3, placebo) \cr
}
}
\source{
Vansteenkiste JF, Cho BC, Vanakesa T, et al. Efficacy of the MAGE-A3
cancer immunotherapeutic as adjuvant therapy in patients with
resected MAGE-A3-positive non-small-cell lung cancer (MAGRIT): a
randomised, double-blind, placebo-controlled, phase 3 trial. Lancet
Oncol 2016; 17: 822–35.
}
\usage{
MAGRIT_2C
}
\description{
Kaplan-Meier digitized data from MAGRIT, figure 2C (PMID 27132212). A reported sample size of 2,312 for a primary endpoint of DFS(3 primary in diff pop) in lung cancer.
}
\examples{
summary(MAGRIT_2C)

kmplot(MAGRIT_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LEA_2C}
\alias{LEA_2C}
\title{LEA, figure 2C}
\format{
A data frame of 374 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (et, et_b) \cr
}
}
\source{
Martín M, Loibl S, von Minckwitz G, et al. Phase III trial evaluating
the addition of bevacizumab to endocrine therapy as first-line
treatment for advanced breast cancer: the letrozole/fulvestrant and
avastin (LEA) study. J Clin Oncol 2015; 33: 1045–52.
}
\usage{
LEA_2C
}
\description{
Kaplan-Meier digitized data from LEA, figure 2C (PMID 25691671). A reported sample size of 374 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(LEA_2C)

kmplot(LEA_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ERCC1_2A}
\alias{ERCC1_2A}
\title{ERCC1, figure 2A}
\format{
A data frame of 464 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cisplatin, p1) \cr
}
}
\source{
Lee SM, Falzon M, Blackhall F, et al. Randomized Prospective
Biomarker Trial of ERCC1 for Comparing Platinum and Nonplatinum
Therapy in Advanced Non-Small-Cell Lung Cancer: ERCC1 Trial (ET). J
Clin Oncol 2017; 35: 402–11.
}
\usage{
ERCC1_2A
}
\description{
Kaplan-Meier digitized data from ERCC1, figure 2A (PMID 27893326). A reported sample size of 648 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ERCC1_2A)

kmplot(ERCC1_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{FFCD2001_S1D}
\alias{FFCD2001_S1D}
\title{FFCD2001, figure S1D}
\format{
A data frame of 282 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (classic, simplified) \cr
}
}
\source{
Aparicio T, Lavau-Denes S, Phelip JM, et al. Randomized phase III
trial in elderly patients comparing LV5FU2 with or without irinotecan
for first-line treatment of metastatic colorectal cancer (FFCD
2001-02). Ann Oncol 2016; 27: 121–7.
}
\usage{
FFCD2001_S1D
}
\description{
Kaplan-Meier digitized data from FFCD2001, figure S1D (PMID 26487578). A reported sample size of 212 for a primary endpoint of PFS in colorectal cancer.
}
\examples{
summary(FFCD2001_S1D)

kmplot(FFCD2001_S1D)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AZURE_2A}
\alias{AZURE_2A}
\title{AZURE, figure 2A}
\format{
A data frame of 3,359 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, zoledronic_acid) \cr
}
}
\source{
Coleman R, Cameron D, Dodwell D, et al. Adjuvant zoledronic acid in
patients with early breast cancer: final efficacy analysis of the
AZURE (BIG 01/04) randomised open-label phase 3 trial. Lancet Oncol
2014; 15: 997–1006.
}
\usage{
AZURE_2A
}
\description{
Kaplan-Meier digitized data from AZURE, figure 2A (PMID 25035292). A reported sample size of 3,360 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(AZURE_2A)

kmplot(AZURE_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{X2005B3030_2B}
\alias{X2005B3030_2B}
\title{2005B3030, figure 2B}
\format{
A data frame of 156 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (control, pci) \cr
}
}
\source{
Li N, Zeng Z-F, Wang S-Y, et al. Randomized phase III trial of
prophylactic cranial irradiation versus observation in patients with
fully resected stage IIIA-N2 nonsmall-cell lung cancer and high risk
of cerebral metastases after adjuvant chemotherapy. Ann Oncol 2015;
26: 504–9.
}
\usage{
X2005B3030_2B
}
\description{
Kaplan-Meier digitized data from 2005B3030, figure 2B (PMID 25515658). A reported sample size of 156 for a primary endpoint of DFS in lung cancer.
}
\examples{
summary(X2005B3030_2B)

kmplot(X2005B3030_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{X20050181_2A}
\alias{X20050181_2A}
\title{20050181, figure 2A}
\format{
A data frame of 597 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (folfiri, panitumumab_folfiri) \cr
}
}
\source{
Peeters M, Price TJ, Cervantes A, et al. Final results from a
randomized phase 3 study of FOLFIRI {+/-} panitumumab for second-line
treatment of metastatic colorectal cancer. Ann Oncol 2014; 25:
107–16.
}
\usage{
X20050181_2A
}
\description{
Kaplan-Meier digitized data from 20050181, figure 2A (PMID 24356622). A reported sample size of 1,186 for a primary endpoint of PFS/OS in colorectal cancer.
}
\examples{
summary(X20050181_2A)

kmplot(X20050181_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ACT2_2B}
\alias{ACT2_2B}
\title{ACT2, figure 2B}
\format{
A data frame of 67 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (mut_b, mut_c) \cr
}
}
\source{
Hagman H, Frödin J-E, Berglund Å, et al. A randomized study of
KRAS-guided maintenance therapy with bevacizumab, erlotinib or
metronomic capecitabine after first-line induction treatment of
metastatic colorectal cancer: the Nordic ACT2 trial. Ann Oncol 2016;
27: 140–7.
}
\usage{
ACT2_2B
}
\description{
Kaplan-Meier digitized data from ACT2, figure 2B (PMID 26483047). A reported sample size of 233 for a primary endpoint of PFS-3mo in colorectal cancer.
}
\examples{
summary(ACT2_2B)

kmplot(ACT2_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUXLung3-6_2B}
\alias{LUXLung3-6_2B}
\title{LUXLung3-6, figure 2B}
\format{
A data frame of 364 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (afatinib, gemicitabine-cisplatin) \cr
}
}
\source{
Yang JC-H, Wu Y-L, Schuler M, et al. Afatinib versus cisplatin-based
chemotherapy for EGFR mutation-positive lung adenocarcinoma (LUX-Lung
3 and LUX-Lung 6): analysis of overall survival data from two
randomised, phase 3 trials. Lancet Oncol 2015; 16: 141–51.
}
\usage{
`LUXLung3-6_2B`
}
\description{
Kaplan-Meier digitized data from LUXLung3-6, figure 2B (PMID 25589191). A reported sample size of 345,364 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(`LUXLung3-6_2B`)

kmplot(`LUXLung3-6_2B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PEGASE07_2B}
\alias{PEGASE07_2B}
\title{PEGASE07, figure 2B}
\format{
A data frame of 174 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (arma, armb) \cr
}
}
\source{
Gonçalves A, Pierga J-Y, Ferrero J-M, et al. UNICANCER-PEGASE 07
study: a randomized phase III trial evaluating postoperative
docetaxel-5FU regimen after neoadjuvant dose-intense chemotherapy for
treatment of inflammatory breast cancer. Ann Oncol 2015; 26: 1692–7.
}
\usage{
PEGASE07_2B
}
\description{
Kaplan-Meier digitized data from PEGASE07, figure 2B (PMID 25943350). A reported sample size of 174 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(PEGASE07_2B)

kmplot(PEGASE07_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MAPS_2A}
\alias{MAPS_2A}
\title{MAPS, figure 2A}
\format{
A data frame of 448 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (pc, pcb) \cr
}
}
\source{
Zalcman G, Mazieres J, Margery J, et al. Bevacizumab for newly
diagnosed pleural mesothelioma in the Mesothelioma Avastin Cisplatin
Pemetrexed Study (MAPS): a randomised, controlled, open-label, phase
3 trial. Lancet 2016; 387: 1405–14.
}
\usage{
MAPS_2A
}
\description{
Kaplan-Meier digitized data from MAPS, figure 2A (PMID 26719230). A reported sample size of 448 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(MAPS_2A)

kmplot(MAPS_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TEXT-SOFT(2)_2A}
\alias{TEXT-SOFT(2)_2A}
\title{TEXT-SOFT(2), figure 2A}
\format{
A data frame of 4,690 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (exemestane_os, tamoxifen_os) \cr
}
}
\source{
Pagani O, Regan MM, Walley BA, et al. Adjuvant exemestane with
ovarian suppression in premenopausal breast cancer. N Engl J Med
2014; 371: 107–18.
}
\usage{
`TEXT-SOFT(2)_2A`
}
\description{
Kaplan-Meier digitized data from TEXT-SOFT(2), figure 2A (PMID 24881463). A reported sample size of 4,690 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`TEXT-SOFT(2)_2A`)

kmplot(`TEXT-SOFT(2)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PREVAIL_1A}
\alias{PREVAIL_1A}
\title{PREVAIL, figure 1A}
\format{
A data frame of 1,633 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab RPFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (enzalutamide, placebo) \cr
}
}
\source{
Beer TM, Armstrong AJ, Rathkopf DE, et al. Enzalutamide in metastatic
prostate cancer before chemotherapy. N Engl J Med 2014; 371: 424–33.
}
\usage{
PREVAIL_1A
}
\description{
Kaplan-Meier digitized data from PREVAIL, figure 1A (PMID 24881730). A reported sample size of 1,717 for a primary endpoint of PFS/OS in prostate cancer.
}
\examples{
summary(PREVAIL_1A)

kmplot(PREVAIL_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NSABPR04_5B}
\alias{NSABPR04_5B}
\title{NSABPR04, figure 5B}
\format{
A data frame of 1,284 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (no_oxali, oxali) \cr
}
}
\source{
Allegra CJ, Yothers G, O’Connell MJ, et al. Neoadjuvant 5-FU or
Capecitabine Plus Radiation With or Without Oxaliplatin in Rectal
Cancer Patients: A Phase III Randomized Clinical Trial. J Natl Cancer
Inst 2015; 107. DOI:10.1093/jnci/djv248.
}
\usage{
NSABPR04_5B
}
\description{
Kaplan-Meier digitized data from NSABPR04, figure 5B (PMID 26374429). A reported sample size of 1,608 for a primary endpoint of Local-regional tumor control in colorectal cancer.
}
\examples{
summary(NSABPR04_5B)

kmplot(NSABPR04_5B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MONALEESA2_1}
\alias{MONALEESA2_1}
\title{MONALEESA2, figure 1}
\format{
A data frame of 668 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, ribociclib_letrozole) \cr
}
}
\source{
Hortobagyi GN, Stemmer SM, Burris HA, et al. Ribociclib as First-Line
Therapy for HR-Positive, Advanced Breast Cancer. N Engl J Med 2016;
375: 1738–48.
}
\usage{
MONALEESA2_1
}
\description{
Kaplan-Meier digitized data from MONALEESA2, figure 1 (PMID 27717303). A reported sample size of 668 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(MONALEESA2_1)

kmplot(MONALEESA2_1)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TCOG0701_1A}
\alias{TCOG0701_1A}
\title{TCOG0701, figure 1A}
\format{
A data frame of 596 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (docetaxel_cis, s1_cis) \cr
}
}
\source{
Kubota K, Sakai H, Katakami N, et al. A randomized phase III trial of
oral S-1 plus cisplatin versus docetaxel plus cisplatin in Japanese
patients with advanced non-small-cell lung cancer: TCOG0701 CATS
trial. Ann Oncol 2015; 26: 1401–8.
}
\usage{
TCOG0701_1A
}
\description{
Kaplan-Meier digitized data from TCOG0701, figure 1A (PMID 25908605). A reported sample size of 608 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(TCOG0701_1A)

kmplot(TCOG0701_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ELDA_1A}
\alias{ELDA_1A}
\title{ELDA, figure 1A}
\format{
A data frame of 299 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cmf, wdocetaxel) \cr
}
}
\source{
Perrone F, Nuzzo F, Di Rella F, et al. Weekly docetaxel versus CMF as
adjuvant chemotherapy for older women with early breast cancer: final
results of the randomized phase III ELDA trial. Ann Oncol 2015; 26:
675–82.
}
\usage{
ELDA_1A
}
\description{
Kaplan-Meier digitized data from ELDA, figure 1A (PMID 25488686). A reported sample size of 302 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(ELDA_1A)

kmplot(ELDA_1A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{E1199(a)_2C}
\alias{E1199(a)_2C}
\title{E1199(a), figure 2C}
\format{
A data frame of 2,481 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (p1, p3) \cr
}
}
\source{
Sparano JA, Zhao F, Martino S, et al. Long-Term Follow-Up of the
E1199 Phase III Trial Evaluating the Role of Taxane and Schedule in
Operable Breast Cancer. J Clin Oncol 2015; 33: 2353–60.
}
\usage{
`E1199(a)_2C`
}
\description{
Kaplan-Meier digitized data from E1199(a), figure 2C (PMID 26077235). A reported sample size of 4,954 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`E1199(a)_2C`)

kmplot(`E1199(a)_2C`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CHHiP(a)_2B}
\alias{CHHiP(a)_2B}
\title{CHHiP(a), figure 2B}
\format{
A data frame of 2,139 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (60gy, 74gy) \cr
}
}
\source{
Dearnaley D, Syndikus I, Mossop H, et al. Conventional versus
hypofractionated high-dose intensity-modulated radiotherapy for
prostate cancer: 5-year outcomes of the randomised, non-inferiority,
phase 3 CHHiP trial. Lancet Oncol 2016; 17: 1047–60.
}
\usage{
`CHHiP(a)_2B`
}
\description{
Kaplan-Meier digitized data from CHHiP(a), figure 2B (PMID 27339115). A reported sample size of 3,216 for a primary endpoint of bRFS in prostate cancer.
}
\examples{
summary(`CHHiP(a)_2B`)

kmplot(`CHHiP(a)_2B`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{IMELDA_2A}
\alias{IMELDA_2A}
\title{IMELDA, figure 2A}
\format{
A data frame of 185 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bevacizumab, bevacizumab_capecitabine) \cr
}
}
\source{
Gligorov J, Doval D, Bines J, et al. Maintenance capecitabine and
bevacizumab versus bevacizumab alone after initial first-line
bevacizumab and docetaxel for patients with HER2-negative metastatic
breast cancer (IMELDA): a randomised, open-label, phase 3 trial.
Lancet Oncol 2014; 15: 1351–60.
}
\usage{
IMELDA_2A
}
\description{
Kaplan-Meier digitized data from IMELDA, figure 2A (PMID 25273343). A reported sample size of 284 for a primary endpoint of PFS in breast cancer.
}
\examples{
summary(IMELDA_2A)

kmplot(IMELDA_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ACT2_2C}
\alias{ACT2_2C}
\title{ACT2, figure 2C}
\format{
A data frame of 71 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (wt_b, wt_be) \cr
}
}
\source{
Hagman H, Frödin J-E, Berglund Å, et al. A randomized study of
KRAS-guided maintenance therapy with bevacizumab, erlotinib or
metronomic capecitabine after first-line induction treatment of
metastatic colorectal cancer: the Nordic ACT2 trial. Ann Oncol 2016;
27: 140–7.
}
\usage{
ACT2_2C
}
\description{
Kaplan-Meier digitized data from ACT2, figure 2C (PMID 26483047). A reported sample size of 233 for a primary endpoint of PFS-3mo in colorectal cancer.
}
\examples{
summary(ACT2_2C)

kmplot(ACT2_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NCT00938652_2C}
\alias{NCT00938652_2C}
\title{NCT00938652, figure 2C}
\format{
A data frame of 519 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (gc, gci) \cr
}
}
\source{
O’Shaughnessy J, Schwartzberg L, Danso MA, et al. Phase III study of
iniparib plus gemcitabine and carboplatin versus gemcitabine and
carboplatin in patients with metastatic triple-negative breast
cancer. J Clin Oncol 2014; 32: 3840–7.
}
\usage{
NCT00938652_2C
}
\description{
Kaplan-Meier digitized data from NCT00938652, figure 2C (PMID 25349301). A reported sample size of 519 for a primary endpoint of OS/PFS in breast cancer.
}
\examples{
summary(NCT00938652_2C)

kmplot(NCT00938652_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MAINSAIL_2A}
\alias{MAINSAIL_2A}
\title{MAINSAIL, figure 2A}
\format{
A data frame of 1,059 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (lenolidomide, placebo) \cr
}
}
\source{
Petrylak DP, Vogelzang NJ, Budnik N, et al. Docetaxel and prednisone
with or without lenalidomide in chemotherapy-naive patients with
metastatic castration-resistant prostate cancer (MAINSAIL): a
randomised, double-blind, placebo-controlled phase 3 trial. Lancet
Oncol 2015; 16: 417–25.
}
\usage{
MAINSAIL_2A
}
\description{
Kaplan-Meier digitized data from MAINSAIL, figure 2A (PMID 25743937). A reported sample size of 1,059 for a primary endpoint of OS in prostate cancer.
}
\examples{
summary(MAINSAIL_2A)

kmplot(MAINSAIL_2A)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ATTENTION_2B}
\alias{ATTENTION_2B}
\title{ATTENTION, figure 2B}
\format{
A data frame of 307 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (placebo, tivantinib) \cr
}
}
\source{
Yoshioka H, Azuma K, Yamamoto N, et al. A randomized, double-blind,
placebo-controlled, phase III trial of erlotinib with or without a
c-Met inhibitor tivantinib (ARQ 197) in Asian patients with
previously treated stage IIIB/IV nonsquamous nonsmall-cell lung
cancer harboring wild-type epidermal growth factor receptor
(ATTENTION study). Ann Oncol 2015; 26: 2066–72.
}
\usage{
ATTENTION_2B
}
\description{
Kaplan-Meier digitized data from ATTENTION, figure 2B (PMID 26153496). A reported sample size of 460 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ATTENTION_2B)

kmplot(ATTENTION_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PROSE_2B}
\alias{PROSE_2B}
\title{PROSE, figure 2B}
\format{
A data frame of 263 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (chemo, erlotinib) \cr
}
}
\source{
Gregorc V, Novello S, Lazzari C, et al. Predictive value of a
proteomic signature in patients with non-small-cell lung cancer
treated with second-line erlotinib or chemotherapy (PROSE): a
biomarker-stratified, randomised phase 3 trial. Lancet Oncol 2014;
15: 713–21.
}
\usage{
PROSE_2B
}
\description{
Kaplan-Meier digitized data from PROSE, figure 2B (PMID 24831979). A reported sample size of 142 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(PROSE_2B)

kmplot(PROSE_2B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TEXT-SOFT(2)_2D}
\alias{TEXT-SOFT(2)_2D}
\title{TEXT-SOFT(2), figure 2D}
\format{
A data frame of 4,690 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (exemestane_os, tamoxifen_os) \cr
}
}
\source{
Pagani O, Regan MM, Walley BA, et al. Adjuvant exemestane with
ovarian suppression in premenopausal breast cancer. N Engl J Med
2014; 371: 107–18.
}
\usage{
`TEXT-SOFT(2)_2D`
}
\description{
Kaplan-Meier digitized data from TEXT-SOFT(2), figure 2D (PMID 24881463). A reported sample size of 4,690 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`TEXT-SOFT(2)_2D`)

kmplot(`TEXT-SOFT(2)_2D`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ERCC1_2C}
\alias{ERCC1_2C}
\title{ERCC1, figure 2C}
\format{
A data frame of 168 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (cisplatin, paclitaxel_gemcitabine) \cr
}
}
\source{
Lee SM, Falzon M, Blackhall F, et al. Randomized Prospective
Biomarker Trial of ERCC1 for Comparing Platinum and Nonplatinum
Therapy in Advanced Non-Small-Cell Lung Cancer: ERCC1 Trial (ET). J
Clin Oncol 2017; 35: 402–11.
}
\usage{
ERCC1_2C
}
\description{
Kaplan-Meier digitized data from ERCC1, figure 2C (PMID 27893326). A reported sample size of 648 for a primary endpoint of OS in lung cancer.
}
\examples{
summary(ERCC1_2C)

kmplot(ERCC1_2C)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{OPTIMAL_1}
\alias{OPTIMAL_1}
\title{OPTIMAL, figure 1}
\format{
A data frame of 154 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (erlotinib, gem_carboplatin) \cr
}
}
\source{
Zhou C, Wu YL, Chen G, et al. Final overall survival results from a
randomised, phase III study of erlotinib versus chemotherapy as
first-line treatment of EGFR mutation-positive advanced
non-small-cell lung cancer (OPTIMAL, CTONG-0802). Ann Oncol 2015; 26:
1877–83.
}
\usage{
OPTIMAL_1
}
\description{
Kaplan-Meier digitized data from OPTIMAL, figure 1 (PMID 26141208). A reported sample size of 165 for a primary endpoint of PFS in lung cancer.
}
\examples{
summary(OPTIMAL_1)

kmplot(OPTIMAL_1)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SMART_3}
\alias{SMART_3}
\title{SMART, figure 3}
\format{
A data frame of 203 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in days) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (deferred_rt, immediate_rt) \cr
}
}
\source{
Clive AO, Taylor H, Dobson L, et al. Prophylactic radiotherapy for
the prevention of procedure-tract metastases after surgical and
large-bore pleural procedures in malignant pleural mesothelioma
(SMART): a multicentre, open-label, phase 3, randomised controlled
trial. Lancet Oncol 2016; 17: 1094–104.
}
\usage{
SMART_3
}
\description{
Kaplan-Meier digitized data from SMART, figure 3 (PMID 27345639). A reported sample size of 203 for a primary endpoint of PTM-incidence in lung cancer.
}
\examples{
summary(SMART_3)

kmplot(SMART_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{E1199(c)_2A}
\alias{E1199(c)_2A}
\title{E1199(c), figure 2A}
\format{
A data frame of 2,483 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab DFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (d1, p3) \cr
}
}
\source{
Sparano JA, Zhao F, Martino S, et al. Long-Term Follow-Up of the
E1199 Phase III Trial Evaluating the Role of Taxane and Schedule in
Operable Breast Cancer. J Clin Oncol 2015; 33: 2353–60.
}
\usage{
`E1199(c)_2A`
}
\description{
Kaplan-Meier digitized data from E1199(c), figure 2A (PMID 26077235). A reported sample size of 4,954 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(`E1199(c)_2A`)

kmplot(`E1199(c)_2A`)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SAKK4106_1}
\alias{SAKK4106_1}
\title{SAKK4106, figure 1}
\format{
A data frame of 262 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev, no_bev) \cr
}
}
\source{
Koeberle D, Betticher DC, von Moos R, et al. Bevacizumab continuation
versus no continuation after first-line chemotherapy plus bevacizumab
in patients with metastatic colorectal cancer: a randomized phase III
non-inferiority trial (SAKK 41/06). Ann Oncol 2015; 26: 709–14.
}
\usage{
SAKK4106_1
}
\description{
Kaplan-Meier digitized data from SAKK4106, figure 1 (PMID 25605741). A reported sample size of 262 for a primary endpoint of TTP in colorectal cancer.
}
\examples{
summary(SAKK4106_1)

kmplot(SAKK4106_1)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CA097375_3}
\alias{CA097375_3}
\title{CA097375, figure 3}
\format{
A data frame of 499 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in years) \cr
\tab \code{event} \tab OS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (a, b) \cr
}
}
\source{
Love RR, Laudico AV, Van Dinh N, et al. Timing of adjuvant surgical
oophorectomy in the menstrual cycle and disease-free and overall
survival in premenopausal women with operable breast cancer. J Natl
Cancer Inst 2015; 107: djv064.
}
\usage{
CA097375_3
}
\description{
Kaplan-Meier digitized data from CA097375, figure 3 (PMID 25794890). A reported sample size of 740 for a primary endpoint of DFS in breast cancer.
}
\examples{
summary(CA097375_3)

kmplot(CA097375_3)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TURANDOT_3}
\alias{TURANDOT_3}
\title{TURANDOT, figure 3}
\format{
A data frame of 564 observations and 3 variables:
\tabular{lll}{
\tab \code{time} \tab event time (in months) \cr
\tab \code{event} \tab PFS event indicator (\code{0}: no event, \code{1}: event) \cr
\tab \code{arm} \tab treatment arms (bev_capecitabine, bev_paclitaxel) \cr
}
}
\source{
Zielinski C, Láng I, Inbar M, et al. Bevacizumab plus paclitaxel
versus bevacizumab plus capecitabine as first-line treatment for
HER2-negative metastatic breast cancer (TURANDOT): primary endpoint
results of a randomised, open-label, non-inferiority, phase 3 trial.
Lancet Oncol 2016; 17: 1230–9.
}
\usage{
TURANDOT_3
}
\description{
Kaplan-Meier digitized data from TURANDOT, figure 3 (PMID 27501767). A reported sample size of 564 for a primary endpoint of OS in breast cancer.
}
\examples{
summary(TURANDOT_3)

kmplot(TURANDOT_3)
}
\keyword{datasets}
