# rodev <img src='man/figures/logo.png' align="right" height="134.5" />

[![Travis build status](https://travis-ci.com/ropensci/rodev.svg?branch=master)](https://travis-ci.com/ropensci/rodev) [![Coverage status](https://codecov.io/gh/ropensci/rodev/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/rodev?branch=master) [![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

The goal of rodev is to help rOpenSci package developers with common tasks, and to promote best practices like the use of status badges across the entire suite.

## Installation

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/rodev")
```

## Active project

All functions will work on the active project, by default the current directory. It is not possible to change this at the moment.

## Functions

Refer to the [reference](https://docs.ropensci.org/rodev/reference/index.html).

See the [issue tracker](https://github.com/ropensci/rodev/issues) for feature suggestions.

## Meta

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.

# rodev 0.0.1.900

* Functions related to rOpenSci branded pkgdown were removed, no longer needed cf https://ropensci.org/technotes/2019/06/07/ropensci-docs/

* `add_ro_desc()` and `add_ro_footer()` now error since they were used for now defunct recommendations.

* Added a `NEWS.md` file to track changes to the package.
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
make sure someone from the team agrees that itâ€™s a problem. If youâ€™ve found a
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

Please note that the rodev project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

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

Please note that the {{{ package }}} project is released with a
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
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_docs_link.R
\name{use_docs_link}
\alias{use_docs_link}
\title{Add rOpenSci docs URL to DESCRIPTION}
\usage{
use_docs_link()
}
\description{
Add rOpenSci docs URL to DESCRIPTION
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_review_badge.R
\name{use_review_badge}
\alias{use_review_badge}
\title{Use rOpenSci review badge}
\usage{
use_review_badge(onboarding_issue_id)
}
\arguments{
\item{onboarding_issue_id}{Issue number of the onboarding thread of your package}
}
\value{
It will help you add the rOpenSci review badge to your README.
}
\description{
Use rOpenSci review badge
}
\examples{
\dontrun{
use_review_badge(24)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collaboration_templates.R
\name{use_ro_github}
\alias{use_ro_github}
\title{Add collaboration files in .github}
\usage{
use_ro_github(
  package_name = desc::desc_get_field("Package", file = usethis::proj_get())
)
}
\arguments{
\item{package_name}{Package name}
}
\description{
Add issue, PR templates and CONTRIBUTING.md file
in .github
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_ro_fnd.R
\name{add_ro_fnd}
\alias{add_ro_fnd}
\title{Add rOpenSci as funder to DESCRIPTION}
\usage{
add_ro_fnd(path = usethis::proj_get())
}
\arguments{
\item{path}{path to package}
}
\description{
Add rOpenSci as funder to DESCRIPTION
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{repostatus_badges}
\alias{repostatus_badges}
\title{Names and Markdown code of repostatus.org badges}
\format{
A data.frame with columns status (lowercase) and md_code.
}
\source{
\url{http://www.repostatus.org/}
}
\description{
Names, Markdown code, href and src of repostatus.org badges
as retrieved from \url{https://github.com/jantman/repostatus.org/tree/master/badges/latest}
before this package release. See also \url{http://www.repostatus.org/}
The data is distributed with a CC-BY SA license. Refer to \url{http://www.repostatus.org/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_cchecks_badge.R
\name{use_cchecks_badge}
\alias{use_cchecks_badge}
\title{Use CRAN checks badge}
\usage{
use_cchecks_badge()
}
\value{
It will help you add the CRAN checks badge to your README.
}
\description{
Use CRAN checks badge
}
\examples{
\dontrun{
use_cchecks_badge()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_repostatus_badge.R
\name{use_repostatus_badge}
\alias{use_repostatus_badge}
\title{Use repostatus.org Badge}
\usage{
use_repostatus_badge(status)
}
\arguments{
\item{status}{current status of the project, cf details.}
}
\value{
It will help you add the repostatus.org badge to your README.
}
\description{
Use repostatus.org Badge
}
\details{
Possible statuses are \code{rodev::repostatus_badges$status},
\itemize{
\item abandoned
\item active
\item concept
\item inactive
\item moved
\item suspended
\item unsupported
\item wip
}
For more details refer to \url{https://www.repostatus.org}.
}
\examples{
\dontrun{
use_repostatus_badge("wip")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_ro_desc.R
\name{add_ro_desc}
\alias{add_ro_desc}
\title{Add rOpenSci review mention to DESCRIPTION (no longer recommended)}
\usage{
add_ro_desc(path = usethis::proj_get())
}
\arguments{
\item{path}{path to package}
}
\description{
No longer recommended!
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_gitbook.R
\name{read_pkg_guide}
\alias{read_pkg_guide}
\title{Browse packaging guidelines}
\usage{
read_pkg_guide()
}
\description{
Browse packaging guidelines
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_ro_footer.R
\name{use_ro_footer}
\alias{use_ro_footer}
\title{Use rOpenSci footer}
\usage{
use_ro_footer()
}
\description{
No longer recommended!
}
