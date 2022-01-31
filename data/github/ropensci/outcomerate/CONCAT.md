
<!-- badges: start -->
[![Coverage
status](https://codecov.io/gh/ropensci/outcomerate/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/outcomerate?branch=master)
[![Travis build
status](https://travis-ci.org/ropensci/outcomerate.svg?branch=master)](https://travis-ci.org/ropensci/outcomerate)
[![Ropensci
status](https://badges.ropensci.org/213_status.svg)](https://github.com/ropensci/onboarding/issues/213)
[![CRAN
status](https://www.r-pkg.org/badges/version/outcomerate)](https://CRAN.R-project.org/package=outcomerate)
[![R build status](https://github.com/ropensci/outcomerate/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/outcomerate/actions)
<!-- badges: end -->

# outcomerate

`outcomerate` is a lightweight R package that implements the standard
outcome rates for surveys, as defined in the [Standard
Definitions](https://www.aapor.org/Standards-Ethics/Standard-Definitions-\(1\).aspx)
of the American Association of Public Opinion Research (AAPOR).

Although the mathematical formulas are straightforward, it can get
tedious and repetitive calculating all the rates by hand, especially for
sub-groups of your study. The formulas are similar to one another and so
it is also dangerously easy to make a clerical mistake. The
`outcomerate` package simplifies the analytically workflow by defining
all formulas as a collection of functions.

## Installation

Install the package from CRAN:

``` r
install.packages("outcomerate")
```

Alternatively, install the latest development version via github:

``` r
#install.packages("devtools")
devtools::install_github("ropensci/outcomerate")
```

## Example

Let’s say you try to survey 12 people. After finishing the fieldwork,
you tabulate all your attempts into a table of disposition outcomes:

| code | disposition           | n |
| :--- | :-------------------- | -: |
| I    | Complete interview    | 4 |
| P    | Partial interview     | 2 |
| R    | Refusal and break-off | 1 |
| NC   | Non-contact           | 1 |
| O    | Other                 | 1 |
| UH   | Unknown if household  | 1 |
| NE   | Known ineligible      | 1 |
| UO   | Unknown, other        | 1 |

Using this table, you may wish to report some of the common survey
outcome rates, such as:

  - **Response Rate:** The proportion of your sample that results in an
    interview.
  - **Cooperation Rate:** The proportion of people contacted who
    participate in your survey.
  - **Refusal Rate:** The proportion of your sample that refused to
    participate.
  - **Contact Rate:** The proportion of sampled cases where you manage
    to reach the respondent.
  - **Location Rate:** The proportion of cases (say, in an establishment
    survey) that you manage to locate.

Most of these rates come under a number of variants, having definitions
that are standardized by AAPOR. The `outcomerate` function lets your
calculate these rates seamlessly:

``` r
# load package
library(outcomerate)

# set counts per disposition code (needs to be a named vector)
freq <- c(I = 4, P = 2, R = 1, NC = 1, O = 1, UH = 1, UO = 1, NE = 1)

# calculate rates, assuming 90% of unknown cases are elligble
outcomerate(freq, e = eligibility_rate(freq))
#>   RR1   RR2   RR3   RR4   RR5   RR6 COOP1 COOP2 COOP3 COOP4  REF1  REF2  REF3 
#> 0.364 0.545 0.370 0.556 0.444 0.667 0.500 0.750 0.571 0.857 0.091 0.093 0.111 
#>  CON1  CON2  CON3  LOC1  LOC2 
#> 0.727 0.741 0.889 0.818 0.833
```

Dispositions do not always come in a tabulated format. Survey analysts
often work with microdata directly, where each row represents an
interview. The `outcomerate` package allows you to obtain rates using
such a format as well:

``` r
# define a vector of dispositions
x <- c("I", "P", "I", "UO", "R", "I", "NC", "I", "O", "P", "UH")

# calculate desired rates
outcomerate(x, rate = c("RR2", "CON1"))
#>  RR2 CON1 
#> 0.55 0.73

# obtain a weighted rate
w <- c(rep(1.3, 6), rep(2.5, 5))
outcomerate(x, weight = w, rate = c("RR2", "CON1"))
#>  RR2w CON1w 
#>  0.50  0.69
```

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# outcomerate (development version)

* Minor updates to vignettes to reflect changes in the tidyverse.

# outcomerate 1.0.1

#### Documentation

* Added CITATION details to the package
* Added documentation as `pkgdown` site

# outcomerate 1.0.0

#### New Features

* `eligibility_rate()` function added to estimate the proportion of eligible cases from the unknowns, based on the known ineligibles (`NE`'s).

#### Improvements

* Refactoring of code based on ROpenSci peer review feedback.
* Added S3 method for factors.
* Addition of many more unit tests.
* Additional of more helpful error messages.

#### Breaking Changes

* `weight` argument no longer accepts scalar inputs.
* If weights are provided, the output labels are renamed in the form 'RR2w' instead of "RR2"
* If `rate = NULL` in the function parameters, the default behavior will be to return all possible rates given the other parameters specified.
* Disposition codes now accept "NE" for known ineligibles. Within `outcomerate()`, these are largely ignored, but are used by `eligibility_rate()` to estimate `e`

#### Documentation

* Added documentation for the (internal) `fmat` formula matrix object
* Added documentation on the `middleearth` toy dataset

# outcomerate 0.0.0.9000

* Created `outcomerate` package
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
## Test environments
* local OS X install, R 3.4.4
* ubuntu 14.04 (on travis-ci), R 3.4.4
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
---
output: github_document
editor_options: 
  chunk_output_type: console
---
 <!-- badges: start -->
[![R build status](https://github.com/ropensci/outcomerate/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/outcomerate/actions) [![Coverage status](https://codecov.io/gh/ropensci/outcomerate/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/outcomerate?branch=master) [![Travis build status](https://travis-ci.org/ropensci/outcomerate.svg?branch=master)](https://travis-ci.org/ropensci/outcomerate) [![Ropensci status](https://badges.ropensci.org/213_status.svg)](https://github.com/ropensci/onboarding/issues/213) [![CRAN status](https://www.r-pkg.org/badges/version/outcomerate)](https://CRAN.R-project.org/package=outcomerate)
<!-- badges: end -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# outcomerate

`outcomerate` is a lightweight R package that implements the standard outcome rates for surveys, as defined in the [Standard Definitions](https://www.aapor.org/Standards-Ethics/Standard-Definitions-(1).aspx) of the American Association of Public Opinion Research (AAPOR).

Although the mathematical formulas are straightforward, it can get tedious and repetitive calculating all the rates by hand, especially for sub-groups of your study. The formulas are similar to one another and so it is also dangerously easy to make a clerical mistake. The `outcomerate` package simplifies the analytically workflow by defining all formulas as a collection of functions.


## Installation


Install the package from CRAN:

``` r
install.packages("outcomerate")
```

Alternatively, install the latest development version via github:

``` r
#install.packages("devtools")
devtools::install_github("ropensci/outcomerate")
```

## Example

Let's say you try to survey 12 people. After finishing the fieldwork, you tabulate all your attempts into a table of disposition outcomes:

```{r, echo=FALSE, message=FALSE}
library(dplyr)
library(forcats)
library(knitr)
options(digits = 2)


x <- c("I", "P", "R", "NC", "O", "UH",  "I", "NE", "I", "I", "P", "UO")
data.frame(code = x) %>%
  mutate(across(code, fct_inorder)) %>%
  mutate(disposition = fct_recode(x, `Complete interview` = "I", 
                                  `Partial interview` = "P",
                                  `Refusal and break-off` = "R",
                                  `Non-contact` = "NC",
                                  `Other` = "O",
                                  `Unknown if household` = "UH",
                                  `Unknown, other` = "UO",
                                  `Known ineligible` = "NE")) %>%
  count(code, disposition) %>%
  kable()
```

Using this table, you may wish to report some of the common survey outcome rates, such as:

* __Response Rate:__ The proportion of your sample that results in an interview.
* __Cooperation Rate:__ The proportion of people contacted who participate in your survey.
* __Refusal Rate:__ The proportion of your sample that refused to participate.
* __Contact Rate:__ The proportion of sampled cases where you manage to reach the respondent.
* __Location Rate:__ The proportion of cases (say, in an establishment survey) that you manage to locate.

Most of these rates come under a number of variants, having definitions that are standardized by AAPOR. The `outcomerate` function lets your calculate these rates seamlessly:

```{r example}
# load package
library(outcomerate)

# set counts per disposition code (needs to be a named vector)
freq <- c(I = 4, P = 2, R = 1, NC = 1, O = 1, UH = 1, UO = 1, NE = 1)

# calculate rates, assuming 90% of unknown cases are eligible
outcomerate(freq, e = eligibility_rate(freq))
```

Dispositions do not always come in a tabulated format. Survey analysts often work with microdata directly, where each row represents an interview. The `outcomerate` package allows you to obtain rates using such a format as well:

```{r example2}
# define a vector of dispositions
x <- c("I", "P", "I", "UO", "R", "I", "NC", "I", "O", "P", "UH")

# calculate desired rates
outcomerate(x, rate = c("RR2", "CON1"))

# obtain a weighted rate
w <- c(rep(1.3, 6), rep(2.5, 5))
outcomerate(x, weight = w, rate = c("RR2", "CON1"))
```

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: "Intro to outcomerate"
author: "Rafael Pilliard Hellwig"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Intro to outcomerate}
  %\VignetteEncoding{UTF-8}
  %\VignetteDepends{dplyr, ggplot2, tidyr, knitr, stringr}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  fig.width = 7,
  fig.height = 4,
  collapse = TRUE,
  comment = "#>"
)
options(dplyr.summarise.inform = FALSE)
```



This vignette demonstrates the basic applications of the `outcomerate` package in R. I will draw on the popular `tidyverse` family of packages for the analysis.

```{r, message=FALSE}
# load packages
library(outcomerate)
library(dplyr)
library(tidyr)
library(knitr)
```

To keep things lighthearted, I will use a toy dataset named `middleearth`. The data consists of 1691 rows, each representing an attempt to interview a member of middle earth. Not all elements in the sample resulted in a completed interview, however. Some cases could not be located, others were located but no one was available, some individuals were found but refused to participate, etc. These particular 'dispositions' can be summarized from the `code` variable in the data:

```{r, message=FALSE}
# load dataset
data(middleearth)

# tabulate frequency table of outcomes
kable(count(middleearth, code, outcome))
```

```{r include = FALSE}
attach(as.list(table(middleearth$code)))
```


It is common for survey practitioners to report a number of outcome rates. These rates give an indication as to the quality of the field work. For example, you may want to know the response rate: the proportion of all cases from our intended sample that actually resulted in an interview. 

How might we go about calculating this?

When we inspect our disposition codes, it become apparent that there could be several ways to do this. For example, you may start by using the total number of complete cases (`r I`) and diving this by the number of observations in the data, `r I` / `r nrow(middleearth)` = `r round(I / nrow(middleearth), 2)`. But what about partially completed interviews? If you include those, you would get a rate of  (`r I` + `r P`) / `r nrow(middleearth)` =  `r round((I + P) / nrow(middleearth), 2)`.

It turns out that there are a lot of ways to calculate such outcome rates. Unless we specify exactly what we mean by "response rate", it is easy for claims regarding survey quality to become opaque, lacking comparability with other surveys. For this reason, the American Association for Public Opinion Research (AAPOR) has published a [set of standardized definitions](https://www.aapor.org/Standards-Ethics/Standard-Definitions-(1).aspx) for practitioners. The guide has no fewer than 6 different variants of the 'response rate.' In the our example, the rates we calculated would match to AAPOR's "Response Rate 1" and "Response Rate 2":


$$
\textrm{RR1} = \frac{\textrm{I}}{\textrm{(I + P) + R + O + NC + (UO + UH)}} 
= \frac{`r I`}{(`r I` + `r P`) + `r R` + `r O` + `r NC` + (`r UO` + 0)} 
= `r round(outcomerate(middleearth$code, rate = "RR1"), 2)`
$$

$$
\textrm{RR2} = \frac{\textrm{(I + P)}}{\textrm{(I + P) + R + O + NC + (UO + UH)}} 
= \frac{(`r I` + `r P`)}{(`r I` + `r P`) + `r R` + `r O` + `r NC` + (`r UO` + 0)} 
= `r round(outcomerate(middleearth$code, rate = "RR2"), 2)`
$$

What's more, the guide has multiple definitions for contact rates, refusal rates, and cooperation rates, and weighted rates. It can easily become tedious to look all these up and calculate them by hand. The `outcomerate` package makes it easier by giving all rates (and more) in one go:


```{r}
disp_counts <- c(I = 760, P = 339, R = 59, NC = 288, O = 1, UO = 173, NE = 71) 

e <- eligibility_rate(disp_counts)
outcomerate(disp_counts, e = e)
```

Each of these rates has a precise definition (see `?outcomerate` for details). As we can see, `RR1` and `RR2` match our earlier calculations. In the example, I needed to specify the parameter `e`, the estimated proportion of unknown cases unknowns (`UO`) that were eligible. The `eligibility_rate()` offers a default way to calculate this, but others may be appropriate.

If we had wanted just to return the two rates from above, we could specify this:

```{r}
outcomerate(disp_counts, rate = c("RR1", "RR2"))
```


## More Advanced Uses

In certain situations, you may want to calculate outcome rates based on a vector of codes, rather than a table of frequency counts. It is just as easy to obtain rates this way using `outcomerate`:

```{r}
# print the head of the dataset
head(middleearth)

# calculate rates using codes; should be same result as before
outcomerate(middleearth$code, e = e)
```

Why might we prefer this input format, when it is just as easy to specify the counts?

Well, if we want to calculate outcome rates by some other covariate, we typically need to go back to the original data. For example, here we use `dplyr` and `tidyr` to calculate outcome rates of interest by race:


```{r}
# create a small wrapper function
get_rates <- function(x, ...){
  rlist <- c("RR1", "RR2", "COOP1", "COOP2", "CON1", "REF1", "LOC1")
  as.data.frame(as.list(outcomerate(x, rate = rlist, e = e, ...)))
}

# calculate rates by group
middleearth %>%
  group_by(race) %>%
  summarise(n     = n(),
            Nhat  = sum(svywt),
            rates = list(get_rates(code))) %>%
  unnest(cols = c(rates)) %>%
  kable(digits = 2, caption = "Outcome Rates by Race")
```


### Weighted Outcome Rates

In certain situations, we also wish to produce _weighted_ outcome rates, using the survey weights that are provided in the data. This is easy to do with one additional parameter:

```{r}
# calculate weighted rates by group
middleearth %>%
  group_by(region) %>%
  summarise(n     = n(),
            Nhat  = sum(svywt),
            rates = list(get_rates(code, weight = svywt))) %>%
  unnest(cols = c(rates)) %>%
  kable(digits = 2, caption = "Weighted Outcome Rates by Region")
```

Compare this to the equivalent unweighted estimates, and you see that the results are not the same.

```{r, echo=FALSE}
# calculate weighted rates by group
middleearth %>%
  group_by(region) %>%
  summarise(n     = n(),
            Nhat  = sum(svywt),
            rates = list(get_rates(code))) %>%
  unnest(cols = c(rates)) %>%
  kable(digits = 2, caption = "Unweighted Outcome Rates by Region")
```



### By Date

Lastly, another useful application of grouped analysis is to calculate the rates by date. This allows you to monitor the quality day by day and notice if performance starts to change over time.


```{r}
library(ggplot2)
library(stringr)

# day-by-day quality monitoring
middleearth %>%
  group_by(day) %>%
  summarise(rates = list(get_rates(code))) %>%
  unnest(cols = c(rates)) %>%
  gather(rate, value, -day) %>%
  mutate(type = str_sub(rate, start = -9, end = -2)) %>%
  ggplot(aes(x = day, y = value, colour = rate)) +
  geom_line(size = 1) +
  facet_wrap(~type) +
  labs(title = "Outcome Rates Over Time")
```

In this example, we can see that the contact rate (`CON`) and response rate (`RR`) start to degrade in quality towards day 30. If fieldwork was still continuing, this could be something to look into and attempt to explain and/or redress.


## Variance Estimation

To estimate the errors from estimates generated by `outcomerate()`, the simplest approach is to use the normal approximation. Since outcome rates are nothing more than proportions (or nearly so), their standard error is given by $SE(p) = \sqrt{(p(1-p))/n}$.

```{r}
# first, calculate the outcome rates
(res <- outcomerate(middleearth$code))

# estimate standard errors using the Normal approximation for proportions 
se <- sapply(res, function(p) sqrt((p * (1 - p)) / nrow(middleearth)))
```

With the standard error in hand, we can then construct frequentist confidence intervals:

```{r}
# calculate 95% confidence intervals
rbind(res - (se * 1.96), res + (se * 1.96))
```

Weighted variance estimation in complex surveys require different procedures that go beyond the scope of this vignette. We recommend using `svycontrast()` from the `survey` package to obtain design-based errors that account for elements such as clustering and stratification. Bootstrapping primary sampling units (PSUs) may also be an appropriate method depending on the design at hand. 
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{middleearth}
\alias{middleearth}
\title{middleearth Dataset}
\description{
\code{middlearth} is a toy dataset consisting of 1691 fake survey interviews
conducted in J.R.R. Tolkien's fictional world of Middle Earth.
}
\details{
Variables contained in the data:
\itemize{
\item \strong{code:} one of the outcome codes {I, P, R, NC, O, UH, UO, UH, UO, NE}
\item \strong{outcome:} A human-interpretable label for the \code{code} variable
\item \strong{researcher}: An identifier for the researcher conducting the interview
\item \strong{region}: The region of the respondent (one of five)
\item \strong{Q1}: A hypothetical binary research question posed to respondents
\item \strong{Q2}: A hypothetical continuous scale question posed to respondents
\item \strong{day}: The day the interview took place (1 being the first day of fieldwork)
\item \strong{race}: The race of the respondent in middle earth (Dwarf, Elf, Hobbit, Man, or Wizard)
\item \strong{svywt}: The survey weight (inverse probability of selection)
}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{fmat}
\alias{fmat}
\title{outcomerate Formula Matrix (Internal Data)}
\description{
The \code{fmat} object is the internal dataset used by the \code{outcomerate} package.
It holds all definitions for the outcome rates. With the exception of location
rates, these are taken from the AAPOR Standard Definitions (2016).
}
\details{
The data is a 3-dimensional binary array consisting of:
\itemize{
\item outcome: codes {I, P, R, NC, O, UH, UO, eUH, eUO, NE}
\item rate: the shorthand name for the rate (e.g. RR1)
\item side: numerator (NUM) and denominator (DEN)
}

Given these three dimensions, each outcome rate can be defined as a rational
number (i.e. a fraction) consisting of a summation of frequencies of
outcome codes (where the matrix entries are nonzero).

The input parameters given by the user are {I, P, R, NC, O, UH, UO} and
the parameter 'e'. The parameter e is multiplied by {UH, UO} internally so as to
produce {eUH, eUO}.

The reason for this implementation is:

a) It conforms to a DRY (don't repeat yourself) philosophy by
holding all definitions in one place. These definitions can be used as upstream
inputs to functions/test suites requiring them.

b) It makes it easier to use intermediate steps in the formula calculations.
For instance, it may be of use to a researchers to want to obtain the
numerator/denominators of calculations, instead of only the output.

c) it makes it easy to compare the output

d) It is easier to maintain
}
\examples{
fmat <- outcomerate:::fmat

# Print the dimensions
dimnames(fmat)

# Say we want to know the defintion of Response Rate 2, RR2. We see
# below that the numerator (NUM) column is defined by the entries with a 1,
# or (I + P). Likewise, the denominator (DEN) is defined as
# (I + P + R + NC + O + UH + UO)
fmat[, "RR2", ]


# To use linear algebra, we define a zero-one numerator matrix 'N'
# and a zero-one denominator matrix 'D'. Our count of disposition codes
# is given here manually as 'x' (in the same order as N and D).
N = fmat[ , , 1]
D = fmat[ , , 2]
x <- c(I = 5, P = 2, R = 1, NC = 7, O = 3,
      UH = 4, UO = 8,  NE = 1, eUH = 3, eUO = 6)

# Return all rates
(x \%*\% N) / (x \%*\% D)


# The same thing can be achieved with the apply family of functions
numden <- apply(x * fmat, 2:3, sum)
numden[, 1] / numden[, 2]
}
\references{
\url{https://www.aapor.org/Standards-Ethics/Standard-Definitions-(1).aspx}
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eligibility_rate.R
\name{eligibility_rate}
\alias{eligibility_rate}
\title{Survey Eligibility Rate}
\usage{
eligibility_rate(x, weight = NULL)
}
\arguments{
\item{x}{a character vector of disposition outcomes (I, P, R, NC, O, UH, UO,
U, or NE). Alternatively, a named vector/table of (weighted) disposition
counts.}

\item{weight}{an optional numeric vector that specifies the weight of each
element in 'x' if x is a character vector. If none is provided (the
default), an unweighted estimate is returned.}
}
\description{
Provides an estimate for the proportion of cases of unknown eligibility
that are eligible, as described by \insertCite{vdk}{outcomerate}. The
rate is typically (but not necessarily) calculated on the screener data
or other sources depending on the type of survey, and approaches to
calculating 'e' may therefore differ from one survey to the next.
}
\details{
The present implementation follows the default used in the Excel-based \emph{AAPOR
Outcome Rate Calculator (Version 4.0, May, 2016)} on the basis of
known ineligibles being coded as "NE".

The eligibility rate (ELR) is defined as
\itemize{
\item ELR = (I + P + R + NC + O) / (I + P + R + NC + O + NE)
}
}
\examples{
# load the outcomerate package
library(outcomerate)

# Create a vector of survey dispositions
#
# I  = Complete interview
# P  = Partial interview
# R  = Refusal and break-off
# NC = Non-contact
# O  = Other
# UH = Unknown if household/occupied housing unit
# UO = Unknown, other
# NE = Not eligible
x <- c("I", "P", "I", "NE", "NC", "UH", "I", "R", "UO", "I", "O", "P", "I")

# calculate all rates, assume 80\% of unknown cases are elligble
eligibility_rate(x)

# calculate weighted rates
w <- runif(13, 0, 5)
eligibility_rate(x, weight = w)

# alternatively, provide input as counts
freq <- c(I = 6, P = 2, NC = 3, NE = 1)
eligibility_rate(freq)

}
\references{
\insertRef{aapor}{outcomerate} \insertAllCited
}
\seealso{
\link{outcomerate}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outcomerate.R
\name{outcomerate}
\alias{outcomerate}
\title{AAPOR Survey Outcome Rates}
\usage{
outcomerate(x, e = NULL, rate = NULL, weight = NULL,
  return_nd = FALSE)
}
\arguments{
\item{x}{a character vector of disposition outcomes (I, P, R, NC, O, UH, or
UO). Alternatively, a named vector/table of (weighted) disposition counts.}

\item{e}{a scalar number that specifies the eligibility rate (the estimated
proportion of unknown cases which are eligible). A default method
of calculating 'e' is provided by \code{\link[=eligibility_rate]{eligibility_rate()}}.}

\item{rate}{an optional character vector specifying the rates to be
calculated. If set to NA (the default), all rates are returned.}

\item{weight}{an optional numeric vector that specifies the weight of each
element in 'x' if x is a character vector or factor. If none is provided (the
default), an unweighted estimate is returned.}

\item{return_nd}{a logical to switch to having the function return the
numerator and denominator instead of the rate. Defaults to FALSE.}
}
\description{
Provides standardized outcome rates for surveys, primarily as defined by the
\href{http://www.aapor.org/}{American Association for Public Opinion Research (AAPOR)}. Details can be
found in the Standard Definitions manual \insertCite{aapor}{outcomerate}.
}
\details{
Survey and public opinion research often categorizes interview attempts of
of a survey according to a set of outcome codes as follows:
\itemize{
\item I  = Complete interview
\item P  = Partial interview
\item R  = Refusal and break-off
\item NC = Non-contact
\item O  = Other
\item UH = Unknown if household/occupied housing unit
\item UO = Unknown, other
\item NE = Known ineligible
}

These high-level classes are used to calculate outcome rates that
provide some measure of quality over the fieldwork. These outcome rates
are defined here as follows:

\strong{AAPOR Response Rate}

The proportion of your intended sample that participate in the survey.
\itemize{
\item RR1 = I / ((I + P) + (R + NC + O) + (UH + UO))
\item RR2 = (I + P) / ((I + P) + (R + NC + O) + (UH + UO))
\item RR3 = I / ((I + P) + (R + NC + O) + e(UH + UO))
\item RR4 = (I + P) / ((I + P) + (R + NC + O) + e(UH + UO))
\item RR5 = I / ((I + P) + (R + NC + O))
\item RR6 = (I + P) / ((I + P) + (R + NC + O))
}

\strong{AAPOR Cooperation Rates}

The proportion of contacted respondents who participate in the survey.
\itemize{
\item COOP1 = I / ((I + P) + R + O)
\item COOP2 = (I + P) / ((I + P) + R + O)
\item COOP3 = I / ((I + P) + R)
\item COOP4 = (I + P) / ((I + P) + R)
}

\strong{AAPOR Refusal Rates}

The proportion of the sample that refuses to participate in the survey.
\itemize{
\item REF1 = R / ((I + P) + (R + NC + O) + (UH + UO))
\item REF2 = R / ((I + P) + (R + NC + O) + e(UH + UO))
\item REF3 = R / ((I + P) + (R + NC + O))
}

\strong{AAPOR Contact Rates}

The proportion of the sample that is successfully contacted for
an interview (whether they chose to participate or not).
\itemize{
\item CON1 = ((I + P) + (R + O)) / ((I + P) + (R + NC + O) + (UH+ UO))
\item CON2 = ((I + P) + (R + O)) / ((I + P) + (R + NC + O) + e(UH + UO))
\item CON3 = ((I + P) + (R + O)) / ((I + P) + (R + NC + O))
}

\strong{Location Rate}

The proportion of cases that could be located for an interview.

The location rate is not defined in AAPOR's Standards, but can be found in
\insertCite{vdk}{outcomerate}. Note: depending on how the
located cases are encoded, this may or may not be the correct formula.
\itemize{
\item LOC1 = ((I + P) + (R + O + NC)) / ((I + P) + (R + NC + O) + (UH + UO))
\item LOC2 = ((I + P) + (R + O + NC)) / ((I + P) + (R + NC + O) + e(UH + UO))
}
}
\examples{
# load the outcomerate package
library(outcomerate)

# Create a vector of survey dispositions
#
# I  = Complete interview
# P  = Partial interview
# R  = Refusal and break-off
# NC = Non-contact
# O  = Other
# UH = Unknown if household/occupied housing unit
# UO = Unknown, other
# NE = Known ineligible
x <- c("I", "P", "I", "NC", "UH", "I", "R", "NE",
      "UO", "I", "O", "P", "I")

# calculate all rates
elr <- eligibility_rate(x)
outcomerate(x, e = elr)

# return only one rate
outcomerate(x, rate = "COOP1")

# calculate weighted rates
w <- runif(length(x), 0, 5)
outcomerate(x, e = elr, weight = w)

# alternatively, provide input as counts
freq <- c(I = 6, P = 2, NC = 3, R = 1)
outcomerate(freq, e = elr)
}
\references{
\insertAllCited \insertRef{aapor}{outcomerate}
}
