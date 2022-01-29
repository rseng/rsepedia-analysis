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
