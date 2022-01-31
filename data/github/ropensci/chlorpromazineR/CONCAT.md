# chlorpromazineR: Convert antipsychotic doses to chlorpromazine equivalents

[![cran checks](https://cranchecks.info/badges/summary/chlorpromazineR)](https://cran.r-project.org/web/checks/check_results_chlorpromazineR.html) [![R-CMD-check](https://github.com/ropensci/chlorpromazineR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/chlorpromazineR/actions) [![Coverage status](https://codecov.io/gh/ropensci/chlorpromazineR/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/chlorpromazineR?branch=master) [![DOI](https://zenodo.org/badge/175675220.svg)](https://zenodo.org/badge/latestdoi/175675220) [![](https://badges.ropensci.org/307_status.svg)](https://github.com/ropensci/software-review/issues/307/) 


Studies investigating or controlling for the impact of antipsychotic medications often need to quantify the amount of medication to which an individual is or has been exposed. As different antipsychotics have different potencies, the task is more complicated than using each medication’s daily dosage in milligrams, for example. `chlorpromazineR` is an R package to calculate dose equivalents for common oral and injectable antipsychotic medications based on conversion factors from the published literature. We do not propose to suggest which conversion factors are appropriate to use, or how to interpret the converted data. All users should also refer to the papers from which the conversion factor data originates to determine whether the use of such data is appropriate for their study.

We hope that this package is of use to scientists who do clinical research involving antipsychotic medications. Specifically, the goals of this package are:

* to improve transparency and consistency in calculating chlorpromazine equivalents,
* to reduce human error and improve accuracy,
* and to simplify workflows for large datasets, as from chart reviews of electronic health records.

For further details and usage, please see the [walkthrough vignette](https://htmlpreview.github.io/?https://github.com/ropensci/chlorpromazineR/blob/master/doc/walkthrough.html).

This results from this package should be double checked for accuracy when used in production. We welcome feedback--please contact via eb@ericebrown.com or file an issue.

## Installation

The CRAN release version (recommended) can be installed via the command: `install.packages("chlorpromazineR")`.

The development version of this package can be installed via the command: `devtools::install_github("ropensci/chlorpromazineR")`.

## Usage

Once installed, load package with `library(chlorpromazineR)`. The package's main conversion functions are documented and usage and examples can be seen with `help(to_cpz)` and `help(to_ap)`. It is strongly recommended to read the articles from which the keys are derived, and to verify that the program's output is producing results as expected. This package facilitates bulk conversion, but should be verified to ensure it produces the results as expected.

### Online calculator

An online calculator using this package is available as a [shiny app](http://ap.eebc.ca/).

### Convert data to chlorpromazine equivalents

    participant_ID <- c("P01", "P02", "P03", "P04")
    age <- c(42, 29, 30, 60) # not used in calculation
    antipsychotic <- c("olanzapine", "olanzapine", "quetiapine", "ziprasidone")
    dose <- c(10, 12.5, 300, 60)
    example_oral <- data.frame(participant_ID, age, antipsychotic, dose, 
                               stringsAsFactors = FALSE)
    to_cpz(example_oral, ap_label = "antipsychotic", dose_label = "dose", 
           route = "oral")

## Disclaimer

This package is not for clinical use. The authors assume no liability. All work must be verified independently. Use at own risk.

## Licence

Copyright (C) 2019-2021 Eric E. Brown. This software is licensed under the GPL-3.

## Citation

If you use this package in your scientific paper, please cite this package and the original papers from which the conversion factors are derived. The references can be viewed by using the built-in help function, e.g. `help(gardner2010)`, and as listed below. In addition, a spreadsheet-based tool facilitating dose equivalence conversion has been published by Leucht et al. (2020).

    Davis, J. (1974). Dose equivalence of the anti-psychotic drugs.
    Journal of Psychiatric Research, 11, 65-69.
    <https://doi.org/10.1016/0022-3956(74)90071-5>

    Gardner, D. M., Murphy, A. L., O’Donnell, H., Centorrino, F., & 
    Baldessarini, R. J. (2010). International consensus study of 
    antipsychotic dosing. The American Journal of Psychiatry, 167(6),
    686–693. <https://doi.org/10.1176/appi.ajp.2009.09060802>

    Leucht, S., Samara, M., Heres, S., & Davis, J. M. (2016). Dose
    Equivalents for Antipsychotic Drugs: The DDD Method. Schizophrenia
    Bulletin, 42(suppl_1), S90–S94. <https://doi.org/10.1093/schbul/sbv167>
    
    Leucht, S., Crippa, A., Siafis, S., Patel, M., Orsini, N. & Davis, J. M. 
    (2020). Dose-Response Meta-Analysis of Antipsychotic Drugs for Acute 
    Schizophrenia. American Journal of Psychiatry. 117(4).
    <https://doi.org/10.1176/appi.ajp.2019.19010034>

    Woods, S. (2003). Chlorpromazine Equivalent Doses for the Newer
    Atypical Antipsychotics. Journal of Clinical Psychiatry. 64(6).
    663-667. <https://doi.org/10.4088/JCP.v64n0607>
    
    

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# chlorpromazineR 0.1.2

This is the first CRAN release of the package. The package has just passed peer-review through rOpenSci. It is published in a fully functional state.

# chlorpromazineR 0.1.3

## Major changes

Added new vignette with examples from each of the keys (references).

## Minor changes

Added additional tests and code cleanup.

# chlorpromazineR 0.2.0

## Major changes

Added new key: leucht2020 (<https://doi.org/10.1176/appi.ajp.2019.19010034>)

### Minor changes

Added additional tests to compare results to spreadsheet by Leucht et al. with similar purpose.
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
# CONTRIBUTING #

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
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the chlorpromazineR project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, an use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
title: "Using chlorpromazineR to calculate chlorpromazine-equivalent doses"
author: "Eric Brown, Parita Shah, Julia Kim"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using chlorpromazineR to calculate chlorpromazine-equivalent doses}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(chlorpromazineR)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Background

Studies investigating or controlling for the impact of antipsychotic medications often need to quantify the amount of medication to which an individual is or has been exposed. As different antipsychotics have different potencies, the task is more complicated than using each medication's daily dosage in milligrams, for example. Further complicating the matter is that antipsychotic medications have different formulations including oral medications, and injectable formulations that are either short- or long-acting (also called depots). A commonly used strategy to account for the differing potencies is to calculate the equivalent dose of each medication in terms of a reference medication.

Chlorpromazine (CPZ) is an antipsychotic medication most often used for this purpose. A CPZ-equivalent dose is typically defined as a dose of antipsychotic that is comparable to 100 mg of CPZ.^[Davis, J. M. (1974). Dose equivalence of the antipsychotic drugs. Journal of Psychiatric Research, 1, 65–69.] The total daily dose of a medication expressed in milligrams of CPZ per day is the daily-dose equivalent, and is commonly utilized in both clinical and research settings.^[Danivas, V., & Venkatasubramanian, G. (2013). Current perspectives on chlorpromazine equivalents: Comparing apples and oranges! Indian Journal of Psychiatry, 55(2), 207–208. <https://doi.org/10.4103/0019-5545.111475>] Both strategies allow the dose of one medication, such as quetiapine, to be expressed in equivalents of chlorpromazine.

There are various approaches to defining dose equivalents, including the "defined daily doses" (Methodology WCCfDS, 2013), and Delphie methods.^[Dalkey, N. (1969). The Delphi Method: An Experimental Study of Group Opinion | RAND. Retrieved from <https://www.rand.org/pubs/research_memoranda/RM5888.html>] Defined daily doses, which are produced by the World Health Organization, provide an estimate of the average maintenance antipsychotic daily dose required for a 70-kg adult with psychosis.^[Leucht, S., Samara, M., Heres, S., & Davis, J. M. (2016). Dose Equivalents for Antipsychotic Drugs: The DDD Method. Schizophrenia Bulletin, 42(suppl_1), S90–S94. <https://doi.org/10.1093/schbul/sbv167>] On the other hand, the Delphi method is an approach to reach a consensus to antipsychotic dose equivalents based on a two-stage, questionnaire-based survey of experts.^[Hasson, F., Keeney, S., & McKenna, H. (2000). Research guidelines for the Delphi survey technique. Journal of Advanced Nursing, 32(4), 1008–1015.] Using the doses defined through these methods, equivalency ratio for each antipsychotic can be calculated, with CPZ as a reference medication (i.e. x mg/day of antipsychotic A is equal to 1 mg/day CPZ). 

While there is no single standard method to defining dose equivalence, the International Consensus Study by Gardner et al. (2010), which used the Delphi method,^[Gardner, D. M., Murphy, A. L., O’Donnell, H., Centorrino, F., & Baldessarini, R. J. (2010). International consensus study of antipsychotic dosing. The American Journal of Psychiatry, 167(6), 686–693. <https://doi.org/10.1176/appi.ajp.2009.09060802>] is a widely used resource to calculate equivalency ratios.

## chlorpromazineR

We present here an easy-to-use R package to calculate CPZ daily dose equivalents for common oral and injectable antipsychotic medications. We do not propose to suggest which conversion factors are appropriate to use, or how to interpret the converted data. All users should also refer to the papers from which the conversion factor data originates to determine whether the use of such data is appropriate for their study.

With this package, we hope:

* to improve transparency and consistency in calculating CPZ equivalents,
* to reduce human error and improve accuracy,
* and to simplify workflows for large datasets, as from chart reviews of electronic health records.

We use the data from the International Consensus Study (Gardner et al. 2010) as the default conversion "key". We also include keys generated from the most widely cited papers on CPZ equivalence. The user can also create their own keys (and contribute them to this package), which are `list` objects as described below.

While the default key is `gardner2010`, derived from the 2010 paper by Gardner et al, no equivalency was provided for short acting vs. oral chlorpromazine, so `gardner2010` only includes oral and long acting injectable conversion factors. However, the SAI equivalencies are available in terms of chlorpromazine SAI equivalents, so to use that data knowing this limitation, the key `gardner2010_withsai` is available if needed. The SAI values calculated could then be converted using an external reference value).

### Input data structure

The data is assumed to be stored in a `data.frame` where each row represents a participant, and with the required variables stored in columns. The two variables always required are the antipsychotic names and the dose to convert. The generic name of the medication must be used. There are differing conventions for drug name abbreviations, and different brand names around the world, so brand names are not built into the package.

For oral medications or short-acting injectables, the dose is converted directly, i.e. not accounting for dosing frequencies, which must be done separately. The exception is for long-acting injectables (LAIs), which converts to daily equivalent oral CPZ dose by dividing by the frequency of administration in days. Thus when using LAIs, a third variable, the administration frequency in days, must be available, unless the stored doses have already been divided by the administration frequency.

### Key data structure

The conversion factors are stored in "keys". A key is a `list` object containing 3 named sub-`lists`: `oral`, `sai` and `lai`. 

```{r}
names(gardner2010)
```

Each sub-`list` is itself a named list where the names are the antipsychotic names, and the values are the conversion factors. 

```{r}
names(gardner2010$oral)[1:5]
```

The conversion factor is what a given dose needs to be multiplied by to get the CPZ-equivalent.

```{r}
gardner2010$oral$aripiprazole
```

For example, with aripiprazole, the value is 20, meaning that the CPZ-equivalent dose of aripiprazole 5 mg is `5 * 20 = 100`.

## Tutorial

### Example for one route only

To demonstrate with a simple example using only oral medications, consider the following data:

```{r}
participant_ID <- c("P01", "P02", "P03", "P04")
age <- c(42, 29, 30, 60)
antipsychotic <- c("olanzapine", "olanzapine", "quetiapine", "ziprasidone")
dose <- c(10, 12.5, 300, 60)
example_oral <- data.frame(participant_ID, age, antipsychotic, dose, 
                           stringsAsFactors = FALSE)

example_oral
```

In this case we have `example_oral`, a `data.frame`, with one row per participant, and a variable for the antipsychotic name, labeled `antipsychotic` and a variable for the dose, named `dose`.

To calculate CPZ-equivalent doses, use the `to_cpz()` function:

```{r}
to_cpz(example_oral, ap_label = "antipsychotic", dose_label = "dose", 
       route = "oral")
```

Notice that 2 new columns appear, which are named based on the default arguments in `to_cpz()`. The CPZ-equivalent doses are in the column labeled `cpz_eq`.

### Example for data with oral and injectable routes

The default key from the study by Gardner et al. included dosing equivalents for oral medications, short acting injectables and long-acting injectables. If the original data includes more than one of these routes, the `route` argument can be set to "mixed".

In this case, more information is required to calculate the equivalent doses. Specifically, a variable to indicate the route ("oral", "sai", or "lai") is required, and if long-acting injectables are used, then the injection frequency is also required.* 

*If the doses have already manually been divided by the frequency, then setting the `q` parameter to `1` will forgo the need to specify a variable storing the frequency information.

Here is some example data:

```{r}
example_lai <- example_oral
example_lai$participant_ID <- c("P05", "P06", "P07", "P08")
example_lai$antipsychotic <- c("zuclopenthixol decanoate", 
                               "zuclopenthixol decanoate", 
                               "perphenazine enanthate", "fluspirilene")
example_lai$q <- c(14, 21, 14, 14)
example_lai$dose <- c(200, 200, 50, 6)

example_sai <- example_oral
example_sai$participant_ID <- c("P09", "P10", "P11", "P12")
example_sai$antipsychotic <- c("Chlorpromazine HCl", "Loxapine HCl", 
                               "Olanzapine tartrate", "Olanzapine tartrate")
example_sai$dose <- c(100, 50, 10, 5)

example_oral$q <- NA
example_sai$q <- NA

example_oral$route <- "oral"
example_lai$route <- "lai"
example_sai$route <- "sai"

example_mixed <- rbind(example_oral, example_sai, example_lai)

example_mixed
```

The function `to_cpz()` is used similarly as above, with the extra variable names (route_label, q) specified:

```{r}
to_cpz(example_mixed, ap_label = "antipsychotic", dose_label = "dose", 
       route = "mixed", route_label = "route", q_label = "q", key=gardner2010_withsai) 
```

### Participants on multiple medications

If the included participants are on multiple medications, then the name and dose of the antipsychotic should be specified as separate variables in the `data.frame`. In this case, `to_cpz()` can be used multiple times, each time specifying the appropriate variables. Run it first specifying the name and dose of the first medication, and then run it again specifying the name and dose of the second medication. Be sure to change the name of the `eq_label` where the result will be stored.

## Checking for inclusion and matching of antipsychotic names

The software depends on an exact match (except for case) of the generic names of antipsychotics. To check if there are rows that do not match, `check_ap()` can be used. Here is an example of data where all the antipsychotic names match the key, showing that 0 is returned (the number of mismatches):

```{r}
check_ap(example_oral, ap_label = "antipsychotic", route = "oral", key=gardner2010)
```

When there are mismatches, they are listed, and the number is returned:

```{r}
example_sai <- example_sai
check_ap(example_sai, ap_label = "antipsychotic", route = "sai", key=gardner2010)
```

As with `to_cpz()`, mixed routes can be used with `check_ap()`:

```{r}
example_mixed_bad <- example_mixed
example_mixed_bad[1,3] <- "olanzzzzapine"
example_mixed_bad[4,3] <- "Seroquel"
example_mixed_bad[5,3] <- "chlorpromazineHCl"
example_mixed_bad[9,3] <- "zuclo decanoate"
check_ap(example_mixed_bad, ap_label = "antipsychotic", 
         route = "mixed", route_label = "route", key=gardner2010_withsai)
```

While the default `gardner2010` key uses the full antipsychotic name for SAI and LAI formulations (e.g. Chlorpromazine HCl) in order to be consistent with the publication and reduce the risk of errors due to ambiguity, the use can modify the key to "trim" the antipsychotic names to just 1 word (e.g. chlorpromazine). 

```{r}
names(gardner2010$lai)
```

To trim to 1-word antipsychotic names, we have made the `trim_key()` function, which takes a key and generates a new key from it with the 1-word antipsychotic labels. If you want to use the `gardner2010` key but want 1-word antipsychotic names, then you can use `trim_key` to generate that new key, which you can then use in `to_cpz()`. For example:

```{r}
trimmed_gardner <- trim_key(gardner2010)

names(trimmed_gardner$lai)
```

## Using other keys

The above examples have used the default CPZ conversion key, `gardner2010`, based on data extracted from Gardner et al. 2010. However, not all antipsychotic medications were included in that study. This package also includes a key extracted from the data included in the Leucht et al. 2016 publication based on the World Health Organization's Defined Daily Doses. To explicitly change the CPZ key, use the key parameter when calling `to_cpz()`.

```{r}
to_cpz(example_oral, ap_label = "antipsychotic", dose_label = "dose", 
       route = "oral", key = leucht2016)
```

However, the authors do not recommend using the DDD values when more scientific values are available. Therefore, it can be desirable to create a new key with data from the `gardner2010` key (or your own created key) combined with the `leucht2016` key, using data from `gardner2010` key where possible. 

For this purpose, we have created the `add_key()` function, which takes a preferred base key from which ALL available values are taken, and adds only the missing antipsychotic conversion values from the second key to be added. Further, there is a `trim` parameter, which invokes the `trim_key()` function on both keys, in the case as in `leucht2016` where only 1-word antipsychotic names are provided.

```{r}
gardner_leucht <- add_key(base = gardner2010, add = leucht2016, trim = TRUE)

names(gardner_leucht$sai)
```

Let's change a value in the `example_oral` data so that there is an antipsychotic that is not listed in the `gardner2010` key: asenapine.

```{r}
other_oral <- example_oral
other_oral[2,3] <- "asenapine"
other_oral[2,4] <- 10

check_ap(other_oral, ap_label = "antipsychotic", route = "oral")
```

We can see that asenapine is not in the default key. We could use the `leucht2016` key:

```{r}
to_cpz(other_oral, ap_label = "antipsychotic", dose_label = "dose", 
       route = "oral", key = gardner_leucht)
```

Note that it used the `gardner2010` values where they were available (olanzapine, quetiapine, ziprasidone) and the `leucht2016` key for asenapine.

To draw from values from additional keys (e.g. user-created keys), the `add_key()` function could be used iteratively. The first argument of `add_key()`, `base`, is always the key with values that take priority.

## List of included keys

Use the built in help function to see the full reference for any key, e.g. `help(davis1974)`. The data may need to be loaded prior to using the help function (e.g. load via `chlorpromazineR::davis1974`, then use `help(davis1974)`.

| Key name | Reference | Description |
|----------|-----------|-------------|
|davis1974 |<https://doi.org/10.1016/0022-3956(74)90071-5>|Classic paper on chlorpromazine equivalent doses.|
|gardner2010|<https://doi.org/10.1176/appi.ajp.2009.09060802>|Consensus study.|
|gardner2010_withsai|<https://doi.org/10.1176/appi.ajp.2009.09060802>|The included SAI data result in equivalent to SAI (parenteral) chlorpromazine, not oral.|
|leucht2016|<https://doi.org/10.1093/schbul/sbv167>|WHO DDD method.|
|leucht2020|<https://doi.org/10.1176/appi.ajp.2019.19010034>|Dose-response meta-analysis.|
|woods2003|<https://doi.org/10.4088/JCP.v64n0607>||

## Converting to equivalents of an arbitrary antipsychotic

Historically, CPZ has been used as the reference antipsychotic. However, using the same key data, any arbitrary reference antipsychotic could be used. Similar to `to_cpz()`, `to_ap()` can perform this conversion. The desired reference antipsychotic can be specified, e.g. with the parameters `convert_to_ap = "olanzapine", convert_to_route = "oral"`.

```{r}
to_ap(other_oral, convert_to_ap = "olanzapine", convert_to_route = "oral", 
      ap_label = "antipsychotic", dose_label = "dose", 
      route = "oral", key = gardner_leucht)
```

---
title: "Examples produced by chlorpromazineR"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Examples produced by chlorpromazineR}
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
library(chlorpromazineR)
```

For a walkthrough on this package's functions and how to use them, see the walkthrough vignette. The purpose of this vignette is to exhibit the antipsychotics included by the conversion keys and their dose equivalents.

## Gardner 2010

For the reference, see `help(gardner2010)`.

```{r}
data_gardner_oral_names <- c("Amisulpride", "Aripiprazole", "Benperidol", "Chlorpromazine",
                             "Clopenthixol", "Clorprothixene", "Clotiapine", "Clozapine",
                             "Droperidol", "Flupenthixol", "Fluphenazine", "Haloperidol",
                             "Levomepromazine", "Loxapine", "Mesoridazine", 
                             "Methotrimeprazine", "Molindone", "Olanzapine", "Oxypertine",
                             "Paliperidone", "Pericyazine", "Perphenazine", "Pimozide",
                             "Prochlorperazine", "Quetiapine", "Remoxipride", "Risperidone",
                             "Sertindole", "Sulpiride", "Thioridazine", "Thiothixene",
                             "Trifluoperazine", "Trifluperidol", "Triflupromazine",
                             "Ziprasidone", "Zotepine", "Zuclopenthixol")

data_gardner_oral_median <- c(700, 30, 5, 600, 60, 500, 100, 400, 10, 10, 12, 10, 400, 
                              60, 300, 300, 100, 20, 240, 9, 50, 30, 8, 88, 750, 212, 6, 
                              20, 800, 500, 30, 20, 2, 100, 160, 300, 50)

data_gardner_oral <- data.frame(ap = data_gardner_oral_names, 
                                dose = data_gardner_oral_median)

```
## Oral route

### Equivalent to olanzapine 20 mg (CPZ 600 mg)
```{r}
to_ap(data_gardner_oral, convert_to_ap = "olanzapine", ap_label = "ap", 
      dose_label = "dose", route = "oral")

```
### Equivalent to chlorpromazine 100 mg 
```{r}
data_gardner_oral_median_cpz100 <- data_gardner_oral_median / 6
data_gardner_oral_cpz100 <- data.frame(ap = data_gardner_oral_names,
                                       dose=data_gardner_oral_median_cpz100)

to_ap(data_gardner_oral_cpz100, convert_to_ap = "olanzapine", 
      ap_label = "ap", dose_label = "dose", route = "oral")

```

## Short-acting injectables

```{r}
data_gardner_sai_names <- c("Chlorpromazine HCl", "Clotiapine injectable",
                            "Fluphenazine HCl", "Haloperidol lactate",
                            "Loxapine HCl", "Mesoridazine besylate",
                            "Olanzapine tartrate", "Perphenazine USP",
                            "Prochlorperazine mesylate", "Promazine HCl",
                            "Trifluoperazine HCl", "Triflupromazine HCl",
                            "Ziprasidone mesylate", "Zuclopenthixol acetate")

data_gardner_sai_median <- c(100, 40, 5, 5, 25, 100, 10, 10, 22, 100, 
                             5, 60, 20, 50)

data_gardner_sai <- data.frame(ap = data_gardner_sai_names, 
                               dose = data_gardner_sai_median)


to_cpz(data_gardner_sai, key=gardner2010_withsai, ap_label = "ap", 
      dose_label = "dose", route = "sai")

```
### Equivalent to haloperidol 5 mg IM
```{r}
to_ap(data_gardner_sai, key=gardner2010_withsai, 
      convert_to_ap = "haloperidol lactate", 
      ap_label = "ap", dose_label = "dose", route = "sai",
      convert_to_route = "sai")
```
    
## Long-acting injectables

```{r}
data_gardner_lai_names <- c("Clopenthixol decanoate", "Flupenthixol decanoate", 
                            "Fluphenazine decanoate", "Fluphenazine enanthate", 
                            "Fluspirilene", "Haloperidol decanoate", 
                            "Perphenazine enanthate", "Pipotiazine palmitate", 
                            "Risperidone microspheres", "Zuclopenthixol decanoate")

data_gardner_lai_median <- c(300, 40, 25, 25, 6, 150, 100, 100, 50, 200)

data_gardner_lai_q <- c(14, 14, 14, 14, 7, 28, 14, 28, 14, 14)

data_gardner_lai <- data.frame(ap = data_gardner_lai_names,
                               dose = data_gardner_lai_median,
                               q = data_gardner_lai_q)

to_cpz(data_gardner_lai, key=gardner2010, ap_label = "ap", 
       dose_label = "dose", route = "lai", q_label = "q")

```

## Davis 1974

For the reference, see `help(davis1974)`.

```{r}
data_davis_names <- c("Chlorpromazine", "Triflupromazine", "Thioridazine", "Prochlorperazine",
                      "Perphenazine", "Fluphenazine", "Trifluoperazine", "Acetophenazine", 
                      "Carphenazine", "Butaperazine", "Mesoridazine", "Piperacetazine", 
                      "Haloperidol", "Chlorprothixene", "Thiothixene")

data_davis_doses <- c(100, 28.4, 95.3, 14.3, 8.9, 1.2, 2.8, 23.5, 24.3, 8.9, 55.3, 10.5, 1.6, 
                     43.9, 5.2)

data_davis_oral <- data.frame(ap = data_davis_names, 
                              dose = data_davis_doses)
```

### Oral

```{r}
to_cpz(data_davis_oral, ap_label = "ap", 
      dose_label = "dose", route = "oral", key=davis1974)
```
### Short acting injectable

The Davis ket converts parenteral (SAI) to oral chlorpromazine equivalents on the basis of the statement in the text that oral is assumed to be 3x the potency of oral. 

```{r}
to_cpz(data_davis_oral, ap_label = "ap", 
      dose_label = "dose", route = "sai", key=davis1974)
```


## Leucht 2016

For the reference, see `help(leucht2016)`.

```{r}

leucht_names <- c("Acepromazine", "Acetophenazine", "Amisulpride", "Aripiprazole", 
                  "Asenapine", "Benperidol", "Bromperidol", "Butaperazine", "Cariprazine",
                  "Chlorproethazine", "Chlorpromazine", "Chlorprothixene", "Clopenthixol",
                  "Clotiapine", "Clozapine", "Cyamemazine", "Dixyrazine", "Droperidol",
                  "Fluanisone", "Flupentixol", "Fluphenazine", "Fluspirilene", "Haloperidol",
                  "Iloperidone", "Levomepromazine", "Levosulpiride", "Loxapine", "Lurasidone",
                  "Melperone", "Mesoridazine", "Molindone", "Moperone", "Mosapramine",
                  "Olanzapine", "Oxypertine", "Paliperidone", "Penfluridol", "Perazine",
                  "Periciazine", "Perphenazine", "Pimozide", "Pipamperone", "Pipotiazine",
                  "Prochlorperazine", "Promazine", "Prothipendyl", "Quetiapine", "Remoxipride",
                  "Risperidone", "Sertindole", "Sulpiride", "Sultopride", "Thiopropazate",
                  "Thioproperazine", "Thioridazine", "Tiapride", "Tiotixene", 
                  "Trifluoperazine", "Trifluperidol", "Triflupromazine", "Veralipride",
                  "Ziprasidone", "Zotepine", "Zuclopenthixol")

leucht_DDD_oral <- c(100, 50, 400, 15, 20, 1.5, 10, 10, NA, NA, 300, 300, 100, 80, 300, NA, 
                     50, NA, NA, 6, 10, NA, 8, NA, 300, 400, 100, 60, 300, 200, 50, 20, NA, 
                     10, 120, 6, 6, 100, 50, 30, 4, 200, 10, 100, 300, 240, 400, 300, 5, 16, 
                     800, 1200, 60, 75, 300, 400, 30, 20, 2, 100, NA, 80, 200, 30)

leucht_DDD_sai <- c(50, NA, NA, 15, NA, NA, 10, NA, NA, NA, 100, 50, 100, 80, 300, NA, 30, 
                    2.5, NA, NA, NA, NA, 8, NA, 100, NA, NA, NA, 300, 200, NA, 20, NA, 10, 
                    NA, NA, NA, 100, 20, 10, NA, NA, NA, 50, 100, 240, NA, 300, NA, NA, 
                    800, NA, NA, 20, NA, 400, NA, 8, NA, 100, NA, 40, NA, 30)

leucht_DDD_lai <- c(NA, NA, NA, NA, NA, NA, 3.3, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                    NA, NA, 4, 1, 0.7, 3.3, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 10.0, 
                    NA, 2.5, NA, NA, NA, 7.0, NA, NA, 5, NA, NA, NA, NA, NA, 2.7, NA, NA, 
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 15)



data_leucht_DDD_oral <- data.frame(ap = leucht_names, 
                              dose = leucht_DDD_oral)

data_leucht_DDD_sai <- data.frame(ap = leucht_names, 
                              dose = leucht_DDD_sai)

# pretend that all are given q 14 days
data_leucht_DDD_lai <- data.frame(ap = leucht_names, 
                                  dose = (leucht_DDD_lai*14), 
                                  q = rep(14, 64))

data_leucht_DDD_oral <- data_leucht_DDD_oral[!is.na(data_leucht_DDD_oral$dose),]
data_leucht_DDD_sai <- data_leucht_DDD_sai[!is.na(data_leucht_DDD_sai$dose),]
data_leucht_DDD_lai <- data_leucht_DDD_lai[!is.na(data_leucht_DDD_lai$dose),]


to_ap(data_leucht_DDD_oral, ap_label = "ap", dose_label = "dose", 
      route = "oral", key=leucht2016, convert_to_ap = "olanzapine")

to_ap(data_leucht_DDD_sai, ap_label = "ap", dose_label = "dose", 
      route = "sai", key=leucht2016, convert_to_ap = "olanzapine", 
      convert_to_route = "sai")

to_ap(data_leucht_DDD_lai, ap_label = "ap", dose_label = "dose", 
      route = "lai", key=leucht2016, convert_to_ap = "olanzapine", q = "q")

```

## Woods 2003

For the reference, see `help(woods2003)`.

```{r}

woods_names <- c("haloperidol", "risperidone", "olanzapine",
                 "quetiapine", "ziprasidone", "aripiprazole")

woods_doses <- c(2, 2, 5, 75, 60, 7.5)

woods_oral <- data.frame(ap = woods_names, 
                         dose = woods_doses)

to_ap(woods_oral, route="oral", ap_label="ap", 
       dose="dose", key=woods2003, 
      convert_to_ap = "olanzapine")

```
---
title: "Using chlorpromazineR to calculate chlorpromazine-equivalent doses"
author: "Eric Brown, Parita Shah, Julia Kim"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using chlorpromazineR to calculate chlorpromazine-equivalent doses}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(chlorpromazineR)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Background

Studies investigating or controlling for the impact of antipsychotic medications often need to quantify the amount of medication to which an individual is or has been exposed. As different antipsychotics have different potencies, the task is more complicated than using each medication's daily dosage in milligrams, for example. Further complicating the matter is that antipsychotic medications have different formulations including oral medications, and injectable formulations that are either short- or long-acting (also called depots). A commonly used strategy to account for the differing potencies is to calculate the equivalent dose of each medication in terms of a reference medication.

Chlorpromazine (CPZ) is an antipsychotic medication most often used for this purpose. A CPZ-equivalent dose is typically defined as a dose of antipsychotic that is comparable to 100 mg of CPZ.^[Davis, J. M. (1974). Dose equivalence of the antipsychotic drugs. Journal of Psychiatric Research, 1, 65–69.] The total daily dose of a medication expressed in milligrams of CPZ per day is the daily-dose equivalent, and is commonly utilized in both clinical and research settings.^[Danivas, V., & Venkatasubramanian, G. (2013). Current perspectives on chlorpromazine equivalents: Comparing apples and oranges! Indian Journal of Psychiatry, 55(2), 207–208. <https://doi.org/10.4103/0019-5545.111475>] Both strategies allow the dose of one medication, such as quetiapine, to be expressed in equivalents of chlorpromazine.

There are various approaches to defining dose equivalents, including the "defined daily doses" (Methodology WCCfDS, 2013), and Delphie methods.^[Dalkey, N. (1969). The Delphi Method: An Experimental Study of Group Opinion | RAND. Retrieved from <https://www.rand.org/pubs/research_memoranda/RM5888.html>] Defined daily doses, which are produced by the World Health Organization, provide an estimate of the average maintenance antipsychotic daily dose required for a 70-kg adult with psychosis.^[Leucht, S., Samara, M., Heres, S., & Davis, J. M. (2016). Dose Equivalents for Antipsychotic Drugs: The DDD Method. Schizophrenia Bulletin, 42(suppl_1), S90–S94. <https://doi.org/10.1093/schbul/sbv167>] On the other hand, the Delphi method is an approach to reach a consensus to antipsychotic dose equivalents based on a two-stage, questionnaire-based survey of experts.^[Hasson, F., Keeney, S., & McKenna, H. (2000). Research guidelines for the Delphi survey technique. Journal of Advanced Nursing, 32(4), 1008–1015.] Using the doses defined through these methods, equivalency ratio for each antipsychotic can be calculated, with CPZ as a reference medication (i.e. x mg/day of antipsychotic A is equal to 1 mg/day CPZ). 

While there is no single standard method to defining dose equivalence, the International Consensus Study by Gardner et al. (2010), which used the Delphi method,^[Gardner, D. M., Murphy, A. L., O’Donnell, H., Centorrino, F., & Baldessarini, R. J. (2010). International consensus study of antipsychotic dosing. The American Journal of Psychiatry, 167(6), 686–693. <https://doi.org/10.1176/appi.ajp.2009.09060802>] is a widely used resource to calculate equivalency ratios.

## chlorpromazineR

We present here an easy-to-use R package to calculate CPZ daily dose equivalents for common oral and injectable antipsychotic medications. We do not propose to suggest which conversion factors are appropriate to use, or how to interpret the converted data. All users should also refer to the papers from which the conversion factor data originates to determine whether the use of such data is appropriate for their study.

With this package, we hope:

* to improve transparency and consistency in calculating CPZ equivalents,
* to reduce human error and improve accuracy,
* and to simplify workflows for large datasets, as from chart reviews of electronic health records.

We use the data from the International Consensus Study (Gardner et al. 2010) as the default conversion "key". We also include keys generated from the most widely cited papers on CPZ equivalence. The user can also create their own keys (and contribute them to this package), which are `list` objects as described below.

While the default key is `gardner2010`, derived from the 2010 paper by Gardner et al, no equivalency was provided for short acting vs. oral chlorpromazine, so `gardner2010` only includes oral and long acting injectable conversion factors. However, the SAI equivalencies are available in terms of chlorpromazine SAI equivalents, so to use that data knowing this limitation, the key `gardner2010_withsai` is available if needed. The SAI values calculated could then be converted using an external reference value).

### Input data structure

The data is assumed to be stored in a `data.frame` where each row represents a participant, and with the required variables stored in columns. The two variables always required are the antipsychotic names and the dose to convert. The generic name of the medication must be used. There are differing conventions for drug name abbreviations, and different brand names around the world, so brand names are not built into the package.

For oral medications or short-acting injectables, the dose is converted directly, i.e. not accounting for dosing frequencies, which must be done separately. The exception is for long-acting injectables (LAIs), which converts to daily equivalent oral CPZ dose by dividing by the frequency of administration in days. Thus when using LAIs, a third variable, the administration frequency in days, must be available, unless the stored doses have already been divided by the administration frequency.

### Key data structure

The conversion factors are stored in "keys". A key is a `list` object containing 3 named sub-`lists`: `oral`, `sai` and `lai`. 

```{r}
names(gardner2010)
```

Each sub-`list` is itself a named list where the names are the antipsychotic names, and the values are the conversion factors. 

```{r}
names(gardner2010$oral)[1:5]
```

The conversion factor is what a given dose needs to be multiplied by to get the CPZ-equivalent.

```{r}
gardner2010$oral$aripiprazole
```

For example, with aripiprazole, the value is 20, meaning that the CPZ-equivalent dose of aripiprazole 5 mg is `5 * 20 = 100`.

## Tutorial

### Example for one route only

To demonstrate with a simple example using only oral medications, consider the following data:

```{r}
participant_ID <- c("P01", "P02", "P03", "P04")
age <- c(42, 29, 30, 60)
antipsychotic <- c("olanzapine", "olanzapine", "quetiapine", "ziprasidone")
dose <- c(10, 12.5, 300, 60)
example_oral <- data.frame(participant_ID, age, antipsychotic, dose, 
                           stringsAsFactors = FALSE)

example_oral
```

In this case we have `example_oral`, a `data.frame`, with one row per participant, and a variable for the antipsychotic name, labeled `antipsychotic` and a variable for the dose, named `dose`.

To calculate CPZ-equivalent doses, use the `to_cpz()` function:

```{r}
to_cpz(example_oral, ap_label = "antipsychotic", dose_label = "dose", 
       route = "oral")
```

Notice that 2 new columns appear, which are named based on the default arguments in `to_cpz()`. The CPZ-equivalent doses are in the column labeled `cpz_eq`.

### Example for data with oral and injectable routes

The default key from the study by Gardner et al. included dosing equivalents for oral medications, short acting injectables and long-acting injectables. If the original data includes more than one of these routes, the `route` argument can be set to "mixed".

In this case, more information is required to calculate the equivalent doses. Specifically, a variable to indicate the route ("oral", "sai", or "lai") is required, and if long-acting injectables are used, then the injection frequency is also required.* 

*If the doses have already manually been divided by the frequency, then setting the `q` parameter to `1` will forgo the need to specify a variable storing the frequency information.

Here is some example data:

```{r}
example_lai <- example_oral
example_lai$participant_ID <- c("P05", "P06", "P07", "P08")
example_lai$antipsychotic <- c("zuclopenthixol decanoate", 
                               "zuclopenthixol decanoate", 
                               "perphenazine enanthate", "fluspirilene")
example_lai$q <- c(14, 21, 14, 14)
example_lai$dose <- c(200, 200, 50, 6)

example_sai <- example_oral
example_sai$participant_ID <- c("P09", "P10", "P11", "P12")
example_sai$antipsychotic <- c("Chlorpromazine HCl", "Loxapine HCl", 
                               "Olanzapine tartrate", "Olanzapine tartrate")
example_sai$dose <- c(100, 50, 10, 5)

example_oral$q <- NA
example_sai$q <- NA

example_oral$route <- "oral"
example_lai$route <- "lai"
example_sai$route <- "sai"

example_mixed <- rbind(example_oral, example_sai, example_lai)

example_mixed
```

The function `to_cpz()` is used similarly as above, with the extra variable names (route_label, q) specified:

```{r}
to_cpz(example_mixed, ap_label = "antipsychotic", dose_label = "dose", 
       route = "mixed", route_label = "route", q_label = "q", key=gardner2010_withsai) 
```

### Participants on multiple medications

If the included participants are on multiple medications, then the name and dose of the antipsychotic should be specified as separate variables in the `data.frame`. In this case, `to_cpz()` can be used multiple times, each time specifying the appropriate variables. Run it first specifying the name and dose of the first medication, and then run it again specifying the name and dose of the second medication. Be sure to change the name of the `eq_label` where the result will be stored.

## Checking for inclusion and matching of antipsychotic names

The software depends on an exact match (except for case) of the generic names of antipsychotics. To check if there are rows that do not match, `check_ap()` can be used. Here is an example of data where all the antipsychotic names match the key, showing that 0 is returned (the number of mismatches):

```{r}
check_ap(example_oral, ap_label = "antipsychotic", route = "oral", key=gardner2010)
```

When there are mismatches, they are listed, and the number is returned:

```{r}
example_sai <- example_sai
check_ap(example_sai, ap_label = "antipsychotic", route = "sai", key=gardner2010)
```

As with `to_cpz()`, mixed routes can be used with `check_ap()`:

```{r}
example_mixed_bad <- example_mixed
example_mixed_bad[1,3] <- "olanzzzzapine"
example_mixed_bad[4,3] <- "Seroquel"
example_mixed_bad[5,3] <- "chlorpromazineHCl"
example_mixed_bad[9,3] <- "zuclo decanoate"
check_ap(example_mixed_bad, ap_label = "antipsychotic", 
         route = "mixed", route_label = "route", key=gardner2010_withsai)
```

While the default `gardner2010` key uses the full antipsychotic name for SAI and LAI formulations (e.g. Chlorpromazine HCl) in order to be consistent with the publication and reduce the risk of errors due to ambiguity, the use can modify the key to "trim" the antipsychotic names to just 1 word (e.g. chlorpromazine). 

```{r}
names(gardner2010$lai)
```

To trim to 1-word antipsychotic names, we have made the `trim_key()` function, which takes a key and generates a new key from it with the 1-word antipsychotic labels. If you want to use the `gardner2010` key but want 1-word antipsychotic names, then you can use `trim_key` to generate that new key, which you can then use in `to_cpz()`. For example:

```{r}
trimmed_gardner <- trim_key(gardner2010)

names(trimmed_gardner$lai)
```

## Using other keys

The above examples have used the default CPZ conversion key, `gardner2010`, based on data extracted from Gardner et al. 2010. However, not all antipsychotic medications were included in that study. This package also includes a key extracted from the data included in the Leucht et al. 2016 publication based on the World Health Organization's Defined Daily Doses. To explicitly change the CPZ key, use the key parameter when calling `to_cpz()`.

```{r}
to_cpz(example_oral, ap_label = "antipsychotic", dose_label = "dose", 
       route = "oral", key = leucht2016)
```

However, the authors do not recommend using the DDD values when more scientific values are available. Therefore, it can be desirable to create a new key with data from the `gardner2010` key (or your own created key) combined with the `leucht2016` key, using data from `gardner2010` key where possible. 

For this purpose, we have created the `add_key()` function, which takes a preferred base key from which ALL available values are taken, and adds only the missing antipsychotic conversion values from the second key to be added. Further, there is a `trim` parameter, which invokes the `trim_key()` function on both keys, in the case as in `leucht2016` where only 1-word antipsychotic names are provided.

```{r}
gardner_leucht <- add_key(base = gardner2010, add = leucht2016, trim = TRUE)

names(gardner_leucht$sai)
```

Let's change a value in the `example_oral` data so that there is an antipsychotic that is not listed in the `gardner2010` key: asenapine.

```{r}
other_oral <- example_oral
other_oral[2,3] <- "asenapine"
other_oral[2,4] <- 10

check_ap(other_oral, ap_label = "antipsychotic", route = "oral")
```

We can see that asenapine is not in the default key. We could use the `leucht2016` key:

```{r}
to_cpz(other_oral, ap_label = "antipsychotic", dose_label = "dose", 
       route = "oral", key = gardner_leucht)
```

Note that it used the `gardner2010` values where they were available (olanzapine, quetiapine, ziprasidone) and the `leucht2016` key for asenapine.

To draw from values from additional keys (e.g. user-created keys), the `add_key()` function could be used iteratively. The first argument of `add_key()`, `base`, is always the key with values that take priority.

## List of included keys

Use the built in help function to see the full reference for any key, e.g. `help(davis1974)`. The data may need to be loaded prior to using the help function (e.g. load via `chlorpromazineR::davis1974`, then use `help(davis1974)`.

| Key name | Reference | Description |
|----------|-----------|-------------|
|davis1974 |<https://doi.org/10.1016/0022-3956(74)90071-5>|Classic paper on chlorpromazine equivalent doses.|
|gardner2010|<https://doi.org/10.1176/appi.ajp.2009.09060802>|Consensus study.|
|gardner2010_withsai|<https://doi.org/10.1176/appi.ajp.2009.09060802>|The included SAI data result in equivalent to SAI (parenteral) chlorpromazine, not oral.|
|leucht2016|<https://doi.org/10.1093/schbul/sbv167>|WHO DDD method.|
|leucht2020|<https://doi.org/10.1176/appi.ajp.2019.19010034>|Dose-response meta-analysis.|
|woods2003|<https://doi.org/10.4088/JCP.v64n0607>||

## Converting to equivalents of an arbitrary antipsychotic

Historically, CPZ has been used as the reference antipsychotic. However, using the same key data, any arbitrary reference antipsychotic could be used. Similar to `to_cpz()`, `to_ap()` can perform this conversion. The desired reference antipsychotic can be specified, e.g. with the parameters `convert_to_ap = "olanzapine", convert_to_route = "oral"`.

```{r}
to_ap(other_oral, convert_to_ap = "olanzapine", convert_to_route = "oral", 
      ap_label = "antipsychotic", dose_label = "dose", 
      route = "oral", key = gardner_leucht)
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_functions.R
\name{check_key}
\alias{check_key}
\title{Check whether a conversion key is the expected format}
\usage{
check_key(key)
}
\arguments{
\item{key}{the key to check}
}
\value{
TRUE if the key is valid, otherwise a error is thrown.
}
\description{
chlorpromazineR uses conversion factors stored in a named list of 3 named
lists. This verifies that the key is in a usable format, which can be helpful
when creating custom keys or modifying included keys.
}
\examples{
check_key(gardner2010)
}
\seealso{
Other key functions: 
\code{\link{add_key}()},
\code{\link{trim_key}()}
}
\concept{key functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_functions.R
\name{add_key}
\alias{add_key}
\title{Combine 2 keys with base key taking precedence}
\usage{
add_key(base, added, trim, verbose = TRUE)
}
\arguments{
\item{base}{the base key}

\item{added}{the key from which other antipsychotics are found to add}

\item{trim}{TRUE to use trim_key on both the base and added key, needed when
one does not use the full names (e.g. leucht2016).}

\item{verbose}{If TRUE, added antipsychotic names will be shown in a message}
}
\value{
a merged key
}
\description{
Use this to combine 2 keys by using the whole "base" key, and adding any
antipsychotics from the "added" key that are not in the "base" key.
}
\examples{
add_key(gardner2010, leucht2016, trim = TRUE)
}
\seealso{
Other key functions: 
\code{\link{check_key}()},
\code{\link{trim_key}()}
}
\concept{key functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_data.R
\docType{data}
\name{gardner2010_withsai}
\alias{gardner2010_withsai}
\title{Chlorpromazine equivalent key from Gardner et al. 2010 data}
\format{
A named list of 3 named lists (1 for each route) and each sub-list
contains the conversion factors for each antipsychotic. The 3 top-level lists
are named `oral`, `sai`, and `lai` (route), and the lists they contain have
names corresponding to the antipsychotic, e.g. `olanzapine`.
}
\source{
Gardner, D. M., Murphy, A. L., O’Donnell, H., Centorrino, F., &
Baldessarini, R. J. (2010). International consensus study of antipsychotic
dosing. The American Journal of #' Psychiatry, 167(6), 686–693.
<https://doi.org/10.1176/appi.ajp.2009.09060802>
}
\usage{
gardner2010_withsai
}
\description{
A list of antipsychotics and their chlorpromazine equivalent doses, generated
from the following file included with the package:
system.file("extdata", "gardner2010.csv", package="chlorpromazineR").
}
\details{
The SAI equivalents produced by this key are equivalent to chlorpromazine SAI
not oral. They could be manually converted.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_data.R
\docType{data}
\name{leucht2020}
\alias{leucht2020}
\title{Antipsychotic equivalent key from Leucht et al. 2020 data}
\format{
A named list of 3 named lists (1 for each route) and each sub-list
contains the conversion factors for each antipsychotic. The 3 top-level lists
are named `oral`, `sai`, and `lai` (route), and the lists they contain have
names corresponding to the antipsychotic, e.g. `olanzapine`.
}
\source{
Leucht, S., Crippa, A., Siafis, S., Patel, M., Orsini, N. & Davis,
J. M. (2020). Dose-Response Meta-Analysis of Antipsychotic Drugs for Acute
Schizophrenia. American Journal of Psychiatry. 117(4).
<https://doi.org/10.1176/appi.ajp.2019.19010034>
}
\usage{
leucht2020
}
\description{
A list of antipsychotics and their olanzapine-equivalent doses, generated
from the following file included with the package:
system.file("extdata", "leucht2020.csv", package="chlorpromazineR").
}
\details{
This reference does not include chlorpromazine, so the conversion from
leucht2016 is implied (i.e. chlorpromazine = olanzapine * 30).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_functions.R
\name{trim_key}
\alias{trim_key}
\title{Modify the names in a conversion key to only include the first word}
\usage{
trim_key(key)
}
\arguments{
\item{key}{the key to trim}
}
\value{
the key that was trimmed (a named list of 3 named lists)
}
\description{
For parenteral (sai) and long-acting/depot (lai) antipsychotics, the name 
consists of the usual generic name (such as haloperidol) and a second word
describing the formulation (e.g. haloperidol decanoate). Since to_cpz() and
add_key() require exact matches to work properly, removing the second word
may be required, but should be done with care as it can add ambiguity (e.g.
fluphenazine enanthate and decanoate).
}
\examples{
trim_key(gardner2010)
}
\seealso{
Other key functions: 
\code{\link{add_key}()},
\code{\link{check_key}()}
}
\concept{key functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_data.R
\docType{data}
\name{leucht2016}
\alias{leucht2016}
\title{Chlorpromazine equivalent key from Leucht et al. 2016 data}
\format{
A named list of 3 named lists (1 for each route) and each sub-list
contains the conversion factors for each antipsychotic. The 3 top-level lists
are named `oral`, `sai`, and `lai` (route), and the lists they contain have
names corresponding to the antipsychotic, e.g. `olanzapine`.
}
\source{
Leucht, S., Samara, M., Heres, S., & Davis, J. M. (2016). Dose
Equivalents for Antipsychotic Drugs: The DDD Method. Schizophrenia Bulletin,
42(suppl_1), S90–S94. <https://doi.org/10.1093/schbul/sbv167>
}
\usage{
leucht2016
}
\description{
A list of antipsychotics and their chlorpromazine equivalent doses, generated
from the following file included with the package:
system.file("extdata", "leucht2016.csv", package="chlorpromazineR").
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_data.R
\docType{data}
\name{woods2003}
\alias{woods2003}
\title{Chlorpromazine equivalent key from Woods 2003 data}
\format{
A named list of 3 named lists (1 for each route) and each sub-list
contains the conversion factors for each antipsychotic. The 3 top-level lists
are named `oral`, `sai`, and `lai` (route), and the lists they contain have
names corresponding to the antipsychotic, e.g. `olanzapine`.
}
\source{
Scott Woods (2003). Chlorpromazine Equivalent Doses for the Newer
Atypical Antipsychotics. Journal of Clinical Psychiatry. 64(6). 663-667.
<https://doi.org/10.4088/JCP.v64n0607>
}
\usage{
woods2003
}
\description{
A list of antipsychotics and their chlorpromazine equivalent doses, generated
from the following file included with the package:
system.file("extdata", "woods2003.csv", package="chlorpromazineR").
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{check_ap}
\alias{check_ap}
\title{Checks whether antipsychotic names are in the key}
\usage{
check_ap(
  input_data,
  key = chlorpromazineR::gardner2010,
  ap_label,
  route,
  route_label
)
}
\arguments{
\item{input_data}{data.frame with antipsychotic name and dose data}

\item{key}{source of the conversion factors--defaults to Gardner et al. 2010}

\item{ap_label}{column in x that stores antipsychotic name}

\item{route}{options include "oral", "sai", "lai" or "mixed"}

\item{route_label}{if "mixed" route is specified, provide the column that
stores the route information}
}
\value{
number of antipsychotic names in x[,ap_label] that don't match key
}
\description{
Provided a data.frame, x, this checks that the antipsychotic names stored in
the x's variable ap_label are present in the key.
}
\examples{
participant_ID <- c("P01", "P02", "P03", "P04")
age <- c(42, 29, 30, 60) # not used in calculation, just shows other data
                         # can exist in the data.frame
antipsychotic <- c("olanzapine", "olanzapine", "quetiapine", "ziprasidone")
dose <- c(10, 12.5, 300, 60)
example_oral <- data.frame(participant_ID, age, antipsychotic, dose, 
                           stringsAsFactors = FALSE)
check_ap(example_oral, ap_label = "antipsychotic", route = "oral", 
         key = gardner2010)
}
\concept{checking functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_data.R
\docType{data}
\name{gardner2010}
\alias{gardner2010}
\title{Chlorpromazine equivalent key from Gardner et al. 2010 data}
\format{
A named list of 3 named lists (1 for each route) and each sub-list
contains the conversion factors for each antipsychotic. The 3 top-level lists
are named `oral`, `sai`, and `lai` (route), and the lists they contain have
names corresponding to the antipsychotic, e.g. `olanzapine`.
}
\source{
Gardner, D. M., Murphy, A. L., O’Donnell, H., Centorrino, F., &
Baldessarini, R. J. (2010). International consensus study of antipsychotic
dosing. The American Journal of Psychiatry, 167(6), 686–693.
<https://doi.org/10.1176/appi.ajp.2009.09060802>
}
\usage{
gardner2010
}
\description{
A list of antipsychotics and their chlorpromazine equivalent doses, generated
from the following file included with the package:
system.file("extdata", "gardner2010.csv", package="chlorpromazineR").
}
\details{
The SAI data is not included in this key, because the original study did not
specify a conversion factor from chlorpromazine LAI to oral. The alternative
key gardner2010_withsai can be used, which includes the SAI data, but the
chlorpromazine equivalent doses produced are equivalent to chlorpromazine SAI
not chlorpromazine oral. They could be manually converted (e.g. by
multiplying the SAI doses by 3 per equivalence noted by Davis 1974
<https://doi.org/10.1016/0022-3956(74)90071-5>)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{to_cpz}
\alias{to_cpz}
\title{Calculates chlorpromazine-equivalent doses}
\usage{
to_cpz(
  input_data,
  ap_label,
  dose_label,
  route = "oral",
  key = chlorpromazineR::gardner2010,
  eq_label = "cpz_eq",
  factor_label = "cpz_conv_factor",
  route_label = NULL,
  q_label = NULL
)
}
\arguments{
\item{input_data}{data.frame with antipsychotic name and dose data}

\item{ap_label}{column in x that stores antipsychotic name}

\item{dose_label}{column in x that stores dose}

\item{route}{options include "oral", "sai", "lai" or "mixed"}

\item{key}{source of the conversion factors--defaults to Gardner et al. 2010}

\item{eq_label}{the name of the column to be created, to save the 
calculated CPZ-equivalent dose}

\item{factor_label}{the name of the column to be created to store the
conversion factors}

\item{route_label}{if "mixed" route is specified, provide the column that
stores the route information}

\item{q_label}{if long-acting injectable doses are included, provide the 
column that stores the injection frequency (days), or only if the doses have
already been divided, set q_label = 1.}
}
\value{
data.frame with new variables storing conversion factor and 
CPZ-equivalent doses
}
\description{
Given a data.frame containing doses of antipsychotics to_cpz() converts the
doses into the equivalent chlorpromazine (CPZ) doses, using the conversion
factors specified in the key.
}
\details{
The default key is gardner2010 which has data for both oral and long-acting
antipsychotic medications. See help(gardner2010) for the source reference.
}
\examples{
participant_ID <- c("P01", "P02", "P03", "P04")
age <- c(42, 29, 30, 60)
antipsychotic <- c("olanzapine", "olanzapine", "quetiapine", "ziprasidone")
dose <- c(10, 12.5, 300, 60)
example_oral <- data.frame(participant_ID, age, antipsychotic, dose, 
                           stringsAsFactors = FALSE)
to_cpz(example_oral, ap_label = "antipsychotic", dose_label = "dose", 
       route = "oral")
}
\seealso{
Other conversion functions: 
\code{\link{to_ap}()}
}
\concept{conversion functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{to_ap}
\alias{to_ap}
\title{Calculates equivalent doses}
\usage{
to_ap(
  input_data,
  convert_to_ap = "olanzapine",
  convert_to_route = "oral",
  ap_label,
  dose_label,
  route = "oral",
  key = chlorpromazineR::gardner2010,
  cpz_eq_label = "cpz_eq",
  ref_eq_label = "ap_eq",
  factor_label = "cpz_conv_factor",
  route_label = NULL,
  q_label = NULL
)
}
\arguments{
\item{input_data}{data.frame with antipsychotic name and dose data}

\item{convert_to_ap}{name of desired reference antipsychotic}

\item{convert_to_route}{the route of the desired reference antipsychotic}

\item{ap_label}{column in x that stores antipsychotic name}

\item{dose_label}{column in x that stores dose}

\item{route}{options include "oral", "sai", "lai" or "mixed"}

\item{key}{source of the conversion factors--defaults to Gardner et al. 2010}

\item{cpz_eq_label}{the name of the column to be created, to save the
calculated CPZ-equivalent dose}

\item{ref_eq_label}{the name of the column to be created to save the doses in
terms of the specified reference antipsychotic (in convert_to_ap)}

\item{factor_label}{the name of the column to be created to store the
conversion factors}

\item{route_label}{if "mixed" route is specified, provide the column that
stores the route information}

\item{q_label}{if long-acting injectable doses are included, provide the
column that stores the injection frequency (days), or only if the doses have
already been divided, set q_label = 1.}
}
\value{
data.frame with new variables storing conversion factor and
CPZ-equivalent doses
}
\description{
As in to_cpz(), to_ap() converts doses of antipsychotics into equivalent
doses to a reference antipsychotic. Whereas in to_cpz() the reference
antipsychotic is chlorpromazine (CPZ), to_ap() converts to equivalents of an
arbitrary antipsychotic specified as a string to convert_to_ap. Conversion
factors are specified in the key.
}
\examples{
participant_ID <- c("P01", "P02", "P03", "P04")
age <- c(42, 29, 30, 60) # not used in calculation, just shows other data
                         # can exist in the data.frame
antipsychotic <- c("olanzapine", "olanzapine", "quetiapine", "ziprasidone")
dose <- c(10, 12.5, 300, 60)
example_oral <- data.frame(participant_ID, age, antipsychotic, dose,
                           stringsAsFactors = FALSE)
to_ap(example_oral, convert_to_ap="olanzapine", convert_to_route="oral",
      ap_label = "antipsychotic", dose_label = "dose", route = "oral")
}
\seealso{
Other conversion functions: 
\code{\link{to_cpz}()}
}
\concept{conversion functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_data.R
\docType{data}
\name{davis1974}
\alias{davis1974}
\title{Chlorpromazine equivalent key from Davis 1974 data}
\format{
A named list of 3 named lists (1 for each route) and each sub-list
contains the conversion factors for each antipsychotic. The 3 top-level lists
are named `oral`, `sai`, and `lai` (route), and the lists they contain have
names corresponding to the antipsychotic, e.g. `olanzapine`.
}
\source{
John Davis (1974). Dose equivalence of the anti-psychotic drugs.
Journal of Psychiatric Research, 11, 65-69.
<https://doi.org/10.1016/0022-3956(74)90071-5>
}
\usage{
davis1974
}
\description{
A list of antipsychotics and their chlorpromazine equivalent doses, generated
from the following file included with the package:
system.file("extdata", "davis1974.csv", package="chlorpromazineR").
}
\keyword{datasets}
