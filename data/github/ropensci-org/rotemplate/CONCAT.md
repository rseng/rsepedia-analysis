# rotemplate <a href='https://docs.ropensci.org/rotemplate'><img src='man/figures/logo.png' align="right" height="134.5" /></a>

rotemplate provides a custom pkgdown template for rOpenSci packages. We
use this to render sites at `https://docs.ropensci.org`. Please don’t
use it for your own package if it’s not an rOpenSci package (i.e. only
use it for packages listed on <https://ropensci.org/packages/>).

Inspired by [tidytemplate](https://github.com/tidyverse/tidytemplate/)
and [lockedatapkg](https://github.com/lockedatapublished/lockedatapkg).

## How to use `rotemplate`

Documentation rOpenSci packages will automatically be generated from
your master branch and published to <https://docs.ropensci.org>. You
don’t have to do anything to make this work. If you want to test your
site locally use this:

``` r
library(rotemplate)
#install.packages("yourpkg")
rotemplate::build_ropensci_docs("path/to/yourpkg")
```

Everything else can be configured as usual via the `_pkgdown.yml` file
as described in the pkgdown documentation.

If your website is not deploying or you run into another problem, please
open an issue in the [ropensci/docs](https://github.com/ropensci/docs)
repository.

### Mathjax

if you want to use Mathjax you’ll need to specify it in the `pkgdown`
config file like so:

``` yaml
template:
  params:
    mathjax: true
```

## Example sites

-   [`cyphr`](https://docs.ropensci.org/cyphr/)

-   [`drake`](https://docs.ropensci.org/drake/)

-   [`riem`](https://docs.ropensci.org/riem/)

-   [`ropenaq`](https://docs.ropensci.org/ropenaq/)

-   [`rotl`](https://docs.ropensci.org/rotl/)

-   [`stplanr`](https://docs.ropensci.org/stplanr/)

-   [`visdat`](http://visdat.njtierney.com/)

-   [`magick`](https://docs.ropensci.org/magick/)
# rotemplate 0.0.1.9001

* Adapted the template to BS5 usage (pkgdown 2.0.0).

# rotemplate 0.0.1.9000

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

Please note that the rotemplate project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
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
