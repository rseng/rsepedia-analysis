
<!-- README.md is generated from README.Rmd. Please edit that file -->

# skater

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/skater)](https://CRAN.R-project.org/package=skater)

[![biorXiv](https://img.shields.io/badge/biorXiv-10.1101%2F2021.07.21.453083-red)](https://www.biorxiv.org/content/10.1101/2021.07.21.453083v1)

[![DOI](https://zenodo.org/badge/339462170.svg)](https://zenodo.org/badge/latestdoi/339462170)

[![R-CMD-check-stable](https://github.com/signaturescience/skater/workflows/R-CMD-check-stable/badge.svg)](https://github.com/signaturescience/skater/actions)

[![R-CMD-check-dev](https://github.com/signaturescience/skater/workflows/R-CMD-check-dev/badge.svg)](https://github.com/signaturescience/skater/actions)

<!-- badges: end -->

**S**NP-based **K**inship **A**nalysis **T**esting and **E**valuation:
miscellaneous **R** data analysis utilties.

## Installation

Install stable release from CRAN:

``` r
install.packages("skater")
```

Install development version from GitHub:

``` r
remotes::install_github("signaturescience/skater", build_vignettes=TRUE)
```

## Usage

The [“Basic Usage”
vignette](https://signaturescience.github.io/skater/articles/basic_usage.html)
steps through the primary functionality of the `skater` package:

``` r
vignette("basic_usage", package = "skater")
```

Full documentation: <https://signaturescience.github.io/skater/>.
# skater 0.1.1

- Updated vignette to handle future changes in tidyr (thanks @DavisVaughan, #56)
- Updated DESCRIPTION with `URL` and `BugReports` links.

# skater 0.1.0

- Initial release
## Test environments

- Local MacOS install, R 4.0.4
- R hub
    - Fedora Linux, R-devel
    - Ubuntu Linux 20.04.1 LTS, R-release
    - Windows Server 2008 R2 SP1, R-devel

## R CMD check results

- Local `R CMD check`: Status OK, 0 errors, 0 warnings, 0 notes
- R hub: 
    - NOTE, New submission
    - NOTE possibly mis-spelled words in description: IBD, et, al, benchmarking, polymorphism. IBD is defined in the DESCRIPTION as "identical by descent"; benchmarking and polymorphism are spelled correctly; and "et al." is used in a reference before linking to the doi with `<doi:...>`.

## Revisions after initial CRAN inspection

- Added more detailed description about package functionality in DESCRIPTION.
- Defined acronyms in DESCRIPTION.
- Added Signature Science, LLC as `cph` in DESCRIPTION Authors.
- Added reference to Description field of DESCRIPTION in the form: `authors (year) <doi:...>` with reference to preprint describing methods.
- Better explanation of identical by descent (IBD) segment to kinship coefficient math in function documentation.
- Stopped exporting two internal functions, removed examples, clarified documentation.
- Added a return value for `plot_pedigree()` (called for side effects).
- Updated exported functions to ensure `@return` `\value` notes class of the output value and what it means.
