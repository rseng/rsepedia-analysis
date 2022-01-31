<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R build
status](https://github.com/ropenscilabs/allcontributors/workflows/R-CMD-check/badge.svg)](https://github.com/ropenscilabs/allcontributors/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropenscilabs/allcontributors/branch/master/graph/badge.svg)](https://codecov.io/gh/ropenscilabs/allcontributors)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/allcontributors)](https://cran.r-project.org/web/packages/allcontributors)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/allcontributors?color=orange)](https://cran.r-project.org/package=allcontributors)
<!-- badges: end -->

An alternative implementation in R of the original
[`all-contributors`](https://allcontributors.org/) to acknowledge all
contributors in your ‘README’ (or elsewhere). The original is intended
to help acknowledge *all* contributions including those beyond the
contents of an actual repository, such as community or other or
less-tangible organisational contributions. This version only
acknowledges tangible contributions to a repository, but automates that
task to a single function call, in the hope that such simplicity will
spur greater usage. In short: This package can’t do everything the
original does, but it makes what it does much easier.

## Why then?

The original [`all-contributors`](https://allcontributors.org/) is
primarily a bot which responds to commit messages such as
`add @user for <contribution>`, where `<contribution>` is one of the
[recognized types](https://allcontributors.org/docs/en/emoji-key). As
said above, the relative advantage of that original system lies
primarily in the diversity of contribution types able to be
acknowledged, with each type for a given user appearing as a
corresponding [emoji](https://allcontributors.org/docs/en/emoji-key)
below their github avatar as listed on the README. In comparison, this R
package:

1.  Works automatically, by calling `add_contributors()` at any time to
    add or update contributor acknowledgements.
2.  Works locally without any bot integration
3.  Can add contributors to any file, not just the main README
4.  Offers a variety of formats for listing contributors:
    1.  divided into sections by types of contributions, or as a single
        section
    2.  presented as full grids (like [the
        original](https://github.com/all-contributors/all-contributors/blob/master/README.md#contributors-)),
        numbered lists of github user names only, or single text strings
        of comma-separated names.

## Installation

The package is now on CRAN (as of 2nd Dec 2020), so can be installed
with,

``` r
install.packages ("allcontributors")
```

Alternatively, a development version can be installed from remote
repository host systems using any one of the following options:

``` r
# install.packages("remotes")
remotes::install_git("https://git.sr.ht/~ropenscilabs/allcontributors")
remotes::install_bitbucket("mpadge/allcontributors")
remotes::install_gitlab("mpadge/allcontributors")
remotes::install_github("mpadge/allcontributors")
```

The package can then be loaded the usual way:

``` r
library (allcontributors)
```

## Usage

The primary function of the package,
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html),
adds a table of all contributors by default to the main `README.md` file
(and `README.Rmd` if that exists). Tables or lists can be added to other
files by specifying the `files` argument of that function. The
appearance of the contributors table is determined by several parameters
in that function, including:

1.  `type` For the type of contributions to include (code, contributors
    who open issues, contributors who discuss issues).
2.  `num_sections` For whether to present contributors in 1, 2, or 3
    distinct sections, dependent upon which `type`s of contributions are
    to be acknowledged.
3.  `format` Determining whether contributors are presented in a grid
    with associated avatars of each contributor, as in [the
    original](https://github.com/all-contributors/all-contributors/blob/master/README.md#contributors-),
    an enumerated list of github user names only, or a single text
    string of comma-separated names.

Contribution data are obtained by querying the github API, for which a
local key should be set as an environmental variable containing the name
`"GITHUB"` (either via `Sys.setenv()`, or as an equivalent entry in a
file `~/.Renviron`).

If the main `README` file(s) contains a markdown section entitled
`"Contributors"`, the
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html)
function will add a table of contributors there, otherwise it will be
appended to the end of the document(s). If you wish your contributors
table to be somewhere other than at the end of the `README` file(s),
start by adding an empty `"## Contributors` section to the file(s) and
the function will insert the table at that point.

Any time you wish to update your contributor list, simply re-run the
`add_contributors()` function. There’s even an `open_issue` parameter
that will automatically open or update a github issue on your repository
so that contributors will be pinged about them being added to your list
of contributors.

The data used to construct the contributions table can also be extracted
without writing to the `README` file(s) with the function
[`get_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/get_contributors.html):

``` r
get_contributors(org = "ropenscilabs", repo = "allcontributors")
#> ★  Extracting code contributors
#> ✔ Extracted code contributors
#> ★  Extracting github issue contributors
#> ✔ Extracted github issue contributors
#>       logins contributions
#> 1     mpadge           142
#> 2     maelle            NA
#> 3 shamindras            NA
#>                                                                                           avatar
#> 1                                            https://avatars.githubusercontent.com/u/6697851?v=4
#> 2 https://avatars.githubusercontent.com/u/8360597?u=144e03ae2bbe8a69318cb0c6c3f647e25aec6763&v=4
#> 3 https://avatars.githubusercontent.com/u/7627188?u=d05fb551796e6ce6db64ae43cd8ce48a0217ef85&v=4
#>            type
#> 1          code
#> 2 issue_authors
#> 3 issue_authors
```

## Updating Contributor Acknowledgements

“Contributors” sections of files will be automatically updated to
reflect any new contributions by simply calling
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html).
If your contributors have not changed then your lists of
acknowledgements will not be changed. The
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html)
function has an additional parameter which may be set to
`force_update = TRUE` to force lists to be updated regardless of whether
contributions have changed. This can be used to change the formats of
acknowledgements at any time. If anything goes wrong, the easiest way to
replace a contributions section is to simply delete the old ones from
all files, and call
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html)
again.

## More Information

The package has a [single
vignette](https://ropenscilabs.github.io/allcontributors/articles/allcontributors.html)
which visually demonstrates the various formats in which an
“allcontributors” section can be presented.

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the
[`allcontributors`
package](https://github.com/ropenscilabs/allcontributors) following the
[all-contributors](https://allcontributors.org) specification.
Contributions of any kind are welcome!

### Code

<table>
<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://avatars.githubusercontent.com/u/6697851?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=mpadge">mpadge</a>
</td>
</tr>
</table>

### Issues

<table>
<tr>
<td align="center">
<a href="https://github.com/maelle">
<img src="https://avatars.githubusercontent.com/u/8360597?u=144e03ae2bbe8a69318cb0c6c3f647e25aec6763&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Amaelle">maelle</a>
</td>
<td align="center">
<a href="https://github.com/shamindras">
<img src="https://avatars.githubusercontent.com/u/7627188?u=d05fb551796e6ce6db64ae43cd8ce48a0217ef85&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Ashamindras">shamindras</a>
</td>
</tr>
</table>
<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
# v 0.0.2.00x

Minor changes:

- Add new `exclude_issues` parameter to main functions
# CRAN notes for allcontributors_0.0.3 submission

This submission rectifies the single NOTE previously generated on some CRAN systems.

The submission generates no notes or warnings on:

* Ubuntu 20.04: R-oldrelease, R-release
* Windows: R-oldrelease, R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)

---
title: "allcontributors"
output:
  rmarkdown::html_vignette:
    self_contained: no

  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = TRUE,
  message = TRUE,
  width = 120,
  comment = "#>",
  fig.retina = 2,
  fig.path = "README-"
)
```

<!-- badges: start -->

[![R build
status](https://github.com/ropenscilabs/allcontributors/workflows/R-CMD-check/badge.svg)](https://github.com/ropenscilabs/allcontributors/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropenscilabs/allcontributors/branch/master/graph/badge.svg)](https://codecov.io/gh/ropenscilabs/allcontributors)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/allcontributors)](https://cran.r-project.org/web/packages/allcontributors) 
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/allcontributors?color=orange)](https://cran.r-project.org/package=allcontributors)
<!-- badges: end -->

An alternative implementation in R of the original
[`all-contributors`](https://allcontributors.org/) to acknowledge all
contributors in your 'README' (or elsewhere). The original is intended to help
acknowledge *all* contributions including those beyond the contents of an
actual repository, such as community or other or less-tangible organisational
contributions. This version only acknowledges tangible contributions to
a repository, but automates that task to a single function call, in the hope
that such simplicity will spur greater usage. In short: This package can't do
everything the original does, but it makes what it does much easier.

## Why then?

The original [`all-contributors`](https://allcontributors.org/) is primarily
a bot which responds to commit messages such as `add @user for <contribution>`,
where `<contribution>` is one of the [recognized
types](https://allcontributors.org/docs/en/emoji-key). As said above, the
relative advantage of that original system lies primarily in the diversity of
contribution types able to be acknowledged, with each type for a given user
appearing as a corresponding
[emoji](https://allcontributors.org/docs/en/emoji-key) below their github
avatar as listed on the README. In comparison, this R package:

1. Works automatically, by calling `add_contributors()` at any time to add or
   update contributor acknowledgements.
2. Works locally without any bot integration
3. Can add contributors to any file, not just the main README
4. Offers a variety of formats for listing contributors:
   (i) divided into sections by types of contributions, or as a single section
   (ii) presented as full grids (like [the
        original](https://github.com/all-contributors/all-contributors/blob/master/README.md#contributors-)),
        numbered lists of github user names only, or single text strings of
        comma-separated names.

## Installation

The package is now on CRAN (as of 2nd Dec 2020), so can be installed with,
```{r cran-install, eval = FALSE}
install.packages ("allcontributors")
```

Alternatively, a development version can be installed from remote repository
host systems using any one of the following options:

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_git("https://git.sr.ht/~ropenscilabs/allcontributors")
remotes::install_bitbucket("mpadge/allcontributors")
remotes::install_gitlab("mpadge/allcontributors")
remotes::install_github("mpadge/allcontributors")
```

The package can then be loaded the usual way:
```{r load, echo = FALSE, message = FALSE}
devtools::load_all (".", export_all = FALSE)
```
```{r load-fakey, eval = FALSE}
library (allcontributors)
```

## Usage

The primary function of the package,
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html),
adds a table of all contributors by default to the main `README.md` file (and
`README.Rmd` if that exists). Tables or lists can be added to other files by
specifying the `files` argument of that function. The appearance of the
contributors table is determined by several parameters in that function,
including:

1. `type` For the type of contributions to include (code, contributors who open
   issues, contributors who discuss issues).
2. `num_sections` For whether to present contributors in 1, 2, or 3 distinct
   sections, dependent upon which `type`s of contributions are to be
   acknowledged.
3. `format` Determining whether contributors are presented in a grid with
   associated avatars of each contributor, as in [the
   original](https://github.com/all-contributors/all-contributors/blob/master/README.md#contributors-), 
   an enumerated list of github user names only, or a single text string of
   comma-separated names.

Contribution data are obtained by querying the github API, for which a local
key should be set as an environmental variable containing the name `"GITHUB"`
(either via `Sys.setenv()`, or as an equivalent entry in a file `~/.Renviron`).

If the main `README` file(s) contains a markdown section entitled
`"Contributors"`, the
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html)
function will add a table of contributors there, otherwise it will be appended
to the end of the document(s). If you wish your contributors table to be
somewhere other than at the end of the `README` file(s), start by adding an
empty `"## Contributors` section to the file(s) and the function will insert
the table at that point.

Any time you wish to update your contributor list, simply re-run the
`add_contributors()` function. There's even an `open_issue` parameter that will
automatically open or update a github issue on your repository so that
contributors will be pinged about them being added to your list of
contributors.

The data used to construct the contributions table can also be extracted
without writing to the `README` file(s) with the function
[`get_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/get_contributors.html):

```{r get_contributors}
get_contributors(org = "ropenscilabs", repo = "allcontributors")
```

## Updating Contributor Acknowledgements

"Contributors" sections of files will be automatically updated to reflect any new
contributions by simply calling
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html).
If your contributors have not changed then your lists of acknowledgements will
not be changed. The
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html)
function has an additional parameter which may be set to `force_update = TRUE`
to force lists to be updated regardless of whether contributions have changed.
This can be used to change the formats of acknowledgements at any time. If
anything goes wrong, the easiest way to replace a contributions section is to
simply delete the old ones from all files, and call
[`add_contributors()`](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html)
again.

## More Information

The package has a [single
vignette](https://ropenscilabs.github.io/allcontributors/articles/allcontributors.html)
which visually demonstrates the various formats in which an "allcontributors"
section can be presented.


## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.

## Contributors




<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the [`allcontributors` package](https://github.com/ropenscilabs/allcontributors) following the [all-contributors](https://allcontributors.org) specification. Contributions of any kind are welcome!

### Code

<table>

<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://avatars.githubusercontent.com/u/6697851?v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=mpadge">mpadge</a>
</td>
</tr>

</table>


### Issues

<table>

<tr>
<td align="center">
<a href="https://github.com/maelle">
<img src="https://avatars.githubusercontent.com/u/8360597?u=144e03ae2bbe8a69318cb0c6c3f647e25aec6763&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Amaelle">maelle</a>
</td>
<td align="center">
<a href="https://github.com/shamindras">
<img src="https://avatars.githubusercontent.com/u/7627188?u=d05fb551796e6ce6db64ae43cd8ce48a0217ef85&v=4" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Ashamindras">shamindras</a>
</td>
</tr>

</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
---
title: "allcontributors"
author: "Mark Padgham"
date: "`r Sys.Date()`"
output: 
    html_document:
        toc: true
        toc_float: true
        number_sections: false
        theme: flatly
vignette: >
  %\VignetteIndexEntry{allcontributors}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r pkg-load, echo = FALSE, message = FALSE}
library (allcontributors)
```

The main functionality of the [`allcontributors`
package](https://github.com/ropenscilabs/allcontributors) is described in the
main [`README`](https://ropenscilabs.github.io/allcontributors/). This vignette
provides a visual reference for the various options available for formatting
contributors.


## Default Grid Format

The following represents the default format of contributors divided into three
sections ("Code", "Issue Authors", and "Issue Contributors"), with each section
formatted as a grid with seven columns (determined by the [`ncols`
parameter](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html)).
The images ("Avatars") are hyperlinked to the main github pages of each
contributor, and the github names below them are linked to the contributions
made to the package by each contributor. The following uses dummy avatars
simply to reduce the compiled size of this vignette, and also uses dummy names
for all except the first. (The names are dummy only in the sense of being
entirely generic, although they actually do belong to real people - click to
find out.)



### Code

<table>

<tr>
<td align="center">
<a href="https://github.com/mpadge">
<img src="https://github.com/identicons/0.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=mpadge">mpadge</a>
</td>
<td align="center">
<a href="https://github.com/this-person">
<img src="https://github.com/identicons/1.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=this-person">this-person</a>
</td>
<td align="center">
<a href="https://github.com/that-person">
<img src="https://github.com/identicons/2.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=that-person">that-person</a>
</td>
<td align="center">
<a href="https://github.com/somebody">
<img src="https://github.com/identicons/3.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=somebody">somebody</a>
</td>
<td align="center">
<a href="https://github.com//somebody-else">
<img src="https://github.com/identicons/4.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=somebody-else">somebody-else</a>
</td>
<td align="center">
<a href="https://github.com/them">
<img src="https://github.com/identicons/5.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=them">them</a>
</td>
<td align="center">
<a href="https://github.com/others">
<img src="https://github.com/identicons/6.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=others">others</a>
</td>
</tr>

</table>


### Issue Authors

<table>

<tr>
<td align="center">
<a href="https://github.com/nobody">
<img src="https://github.com/identicons/7.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Anobody">nobody</a>
</td>
<td align="center">
<a href="https://github.com/somebody">
<img src="https://github.com/identicons/8.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Asomebody">somebody</a>
</td>
<td align="center">
<a href="https://github.com/anybody">
<img src="https://github.com/identicons/9.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Aanybody">anybody</a>
</td>
<td align="center">
<a href="https://github.com/nope">
<img src="https://github.com/identicons/10.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Anope">nope</a>
</td>
<td align="center">
<a href="https://github.com/yep">
<img src="https://github.com/identicons/11.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Ayep">yep</a>
</td>
<td align="center">
<a href="https://github.com/maybe">
<img src="https://github.com/identicons/12.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Amaybe">maybe</a>
</td>
<td align="center">
<a href="https://github.com/doubtful">
<img src="https://github.com/identicons/13.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Adoubtful">doubtful</a>
</td>
</tr>

</table>


### Issue Contributors

<table>

<tr>

<td align="center">
<a href="https://github.com/here">
<img src="https://github.com/identicons/14.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Ahere">here</a>
</td>
<td align="center">
<a href="https://github.com/there">
<img src="https://github.com/identicons/15.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Athere">there</a>
</td>
<td align="center">
<a href="https://github.com/anywhere">
<img src="https://github.com/identicons/16.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Aanywhere">anywhere</a>
</td>
<td align="center">
<a href="https://github.com/somewhere">
<img src="https://github.com/identicons/17.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Asomewhere">somewhere</a>
</td>
<td align="center">
<a href="https://github.com/nowhere">
<img src="https://github.com/identicons/18.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Anowhere">nowhere</a>
</td>
<td align="center">
<a href="https://github.com/sometime">
<img src="https://github.com/identicons/19.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Asometime">sometime</a>
</td>
<td align="center">
<a href="https://github.com/later">
<img src="https://github.com/identicons/20.png" width="100px;" alt=""/>
</a><br>
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Alater">later</a>
</td>

</tr>

</table>

---

## Section Organisation

The default output shown above has three sections of "Code", "Issue Authors"
and "Issue Contributors". The organisation of these sections can be controlled
by the parameters `num_sections`, `type`, and `section_names`.

The `type` parameter enables sections to be removed by reducing them from the
default three referred to as `code`, `issues` (for those who open issues), and
`discussion` (for those who contribute to issues). For example, passing `type =
"code"` will only acknowledge direct contributions to code, while ignoring all
those who contributed to issues only.

The `num_sections` argument is provided for convenience, primarily in order to
allow default formats to have either one, two, or three sections. Specifying
`num_sections = 2` will by default collapse the "Issue Authors" and "Issue
Contributors" sections into a single section named "Issues". (This section title
may be renamed with the `section_names` parameter.) 

## List Format 

The `format` parameter of the [`add_contributors()`
function](https://ropenscilabs.github.io/allcontributors/reference/add_contributors.html)
accepts the three options of "grid", "list", or "text." With the three default
section titles as shown above, the "list" option gives output that looks like this:


### Code

<ol>

<li><a href="https://github.com/ropenscilabs/allcontributors/commits?author=mpadge">mpadge</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/commits?author=this-person">this-person</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/commits?author=that-person">that-person</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/commits?author=somebody">somebody</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/commits?author=somebody-else">somebody-else</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/commits?author=them">them</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/commits?author=others">others</a></li>

</ol>


### Issue Authors

<ol>

<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Anobody">nobody</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Asomebody">somebody</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Aanybody">anybody</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Anope">nope</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Ayep">yep</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Amaybe">maybe</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Adoubtful">doubtful</a></li>

</ol>


### Issue Contributors

<ol>

<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Ahere">here</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Athere">there</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Aanywhere">anywhere</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Asomewhere">somewhere</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Anowhere">nowhere</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Asometime">sometime</a></li>
<li><a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Alater">later</a></li>

</ol>

## Text Format 

Finally, the text format enables contributors to be acknowledged as a single
lines of text.


### Code


<a href="https://github.com/ropenscilabs/allcontributors/commits?author=mpadge">mpadge</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=this-person">this-person</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=that-person">that-person</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=somebody">somebody</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=somebody-else">somebody-else</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=them">them</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=others">others</a> 



### Issue Authors


<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Anobody">nobody</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Asomebody">somebody</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Aanybody">anybody</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Anope">nope</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Ayep">yep</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Amaybe">maybe</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Adoubtful">doubtful</a> 


### Issue Contributors


<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Ahere">here</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Athere">there</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Aanywhere">anywhere</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Asomewhere">somewhere</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Anowhere">nowhere</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Asometime">sometime</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Alater">later</a>

The shortest possible way of acknowledging your contributors would be like this:
```{r short, eval = FALSE}
add_contributors (num_sectons = 1, format = "text")
```
which would in this case convert the above into the single list of,

<a href="https://github.com/ropenscilabs/allcontributors/commits?author=mpadge">mpadge</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=this-person">this-person</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=that-person">that-person</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=somebody">somebody</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=somebody-else">somebody-else</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=them">them</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/commits?author=others">others</a>,
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Anobody">nobody</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Asomebody">somebody</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Aanybody">anybody</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Anope">nope</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Ayep">yep</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Amaybe">maybe</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+author%3Adoubtful">doubtful</a>,
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Ahere">here</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Athere">there</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Aanywhere">anywhere</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Asomewhere">somewhere</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Anowhere">nowhere</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Asometime">sometime</a>, 
<a href="https://github.com/ropenscilabs/allcontributors/issues?q=is%3Aissue+commenter%3Alater">later</a>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{get_contributors}
\alias{get_contributors}
\title{get_contributors}
\usage{
get_contributors(
  org,
  repo,
  type = c("code", "issues", "discussion"),
  exclude_issues = NULL,
  alphabetical = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{org}{Github organisation name for repository}

\item{repo}{Repository within \code{org} for which contributors are to be
extracted}

\item{type}{Type of contributions to include: 'code' for direct code
contributions (including documentation), 'issues' to recognise contributors
who open issues, and 'discussion' for contributing to discussions within
issues. Discussion contributions are only from individuals not present in
either 'issues' or 'code'; and 'issues' contributions are only from
individuals not present in 'code'.}

\item{exclude_issues}{Numbers of any issues (or pull requests) to be excluded
from lists of contributors.}

\item{alphabetical}{If \code{TRUE}, order contributors alphabetically, otherwise
order by decreasing numbers of contributions.}

\item{quiet}{If \code{FALSE}, display progress information on screen.}
}
\description{
Get all contributors to a repository, including those who contribute to code,
open issues, and contribute to discussions in issues.
}
\examples{
\dontrun{
get_contributors (org = "ropenscilabs", repo = "allcontributors")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{get_gh_code_contributors}
\alias{get_gh_code_contributors}
\title{get_gh_code_contributors}
\usage{
get_gh_code_contributors(org, repo, alphabetical = FALSE)
}
\arguments{
\item{org}{Github organisation name for repository}

\item{repo}{Repository within \code{org} for which contributors are to be
extracted}

\item{alphabetical}{If \code{TRUE}, order contributors alphabetically, otherwise
order by decreasing numbers of contributions.}
}
\value{
A \code{data.frame} of two columns of contributor (name, login)
}
\description{
Get list of all code contributors to the code of a repository
}
\examples{
\dontrun{
get_gh_code_contributors (org = "ropenscilabs", repo = "allcontributors")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{get_gh_issue_people}
\alias{get_gh_issue_people}
\title{get_gh_issue_people}
\usage{
get_gh_issue_people(org, repo, exclude_issues = NULL)
}
\arguments{
\item{org}{Github organisation name for repository}

\item{repo}{Repository within \code{org} for which contributors are to be
extracted}

\item{exclude_issues}{Numbers of any issues (or pull requests) to be excluded
from lists of contributors.}
}
\value{
List of (authors, contributors), each as character vector of github
login names.
}
\description{
Extract lists of (1) all authors of, and (2) all contributors to, all github
issues for nominated repository
}
\examples{
\dontrun{
get_gh_issue_people (org = "ropenscilabs", repo = "allcontributors")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{get_gh_contrib_issue}
\alias{get_gh_contrib_issue}
\title{get_gh_contrib_issue}
\usage{
get_gh_contrib_issue(org, repo)
}
\arguments{
\item{org}{Github organisation name for repository}

\item{repo}{Repository within \code{org} for which contributors are to be
extracted}
}
\value{
Character vector of github logins for all contributors listed in
current issue, or empty character string if there no issue named "All
Contributors".
}
\description{
Extract contributors currently listed on an "All Contributions" issue in a
github repository.
}
\examples{
\dontrun{
get_gh_contrib_issue (org = "ropenscilabs", repo = "allcontributors")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allcontributor-package.R
\docType{package}
\name{allcontributors-package}
\alias{allcontributors}
\alias{allcontributors-package}
\title{allcontributors: Acknowledge all Contributors to a Project}
\description{
Acknowledge all contributors to a project via a single function call. The function appends to a 'README' or other specified file(s) a table with names of all individuals who contributed via code or repository issues. The package also includes several additional functions to extract and quantify contributions to any repository.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropenscilabs/allcontributors}
  \item Report bugs at \url{https://github.com/ropenscilabs/allcontributors/issues}
}

}
\author{
\strong{Maintainer}: Mark Padgham \email{mark.padgham@email.com}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{get_gh_issue_titles}
\alias{get_gh_issue_titles}
\title{get_gh_issue_titles}
\usage{
get_gh_issue_titles(org, repo)
}
\arguments{
\item{org}{Github organisation name for repository}

\item{repo}{Repository within \code{org} for which contributors are to be
extracted}
}
\value{
\code{data.frame} with one column of issue numbers, and one column of
issue titles.
}
\description{
Extract titles and numbers of all issues associated with a nominated
repository
}
\examples{
\dontrun{
get_gh_issue_titles (org = "ropenscilabs", repo = "allcontributors")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add-contributors.R
\name{add_contributors}
\alias{add_contributors}
\title{add_contributors}
\usage{
add_contributors(
  repo = ".",
  ncols = 7,
  files = c("README.Rmd", "README.md"),
  type = c("code", "issues", "discussion"),
  exclude_issues = NULL,
  num_sections = 3,
  section_names = c("Code", "Issue Authors", "Issue Contributors"),
  format = "grid",
  alphabetical = FALSE,
  open_issue = FALSE,
  force_update = FALSE
)
}
\arguments{
\item{repo}{Location of repository for which contributions are to be
extracted. This must be a git project with a github remote.}

\item{ncols}{Number of columns for contributors in 'README'}

\item{files}{Names of files in which to add contributors}

\item{type}{Type of contributions to include: 'code' for direct code
contributions (including documentation), 'issues' to recognise contributors
who open issues, and 'discussion' for contributing to discussions within
issues. Discussion contributions are only from individuals not present in
either 'issues' or 'code'; and 'issues' contributions are only from
individuals not present in 'code'.}

\item{exclude_issues}{Numbers of any issues (or pull requests) to be excluded
from lists of contributors.}

\item{num_sections}{Number of sections in which to divide contributors:
\itemize{
\item{1} All contributions within single section regardless of \code{type}
\item{2} Contributions divided between a single section for \code{code} and a
second section for all other issue-related contributions.
\item{3} Contributions divided into single sections for each of the three
\code{type} arguments.
}}

\item{section_names}{Names of the sections to appear on the nominated
\code{files}.}

\item{format}{One of ("grid", "list", "text") to control display of
contributors as
\itemize{
\item{1} "grid" for a rectangular grid, with each contributor represented by
their github avatar, with avatar linked to contributor profile and name
linked to repository contributions.
\item{2} "list" for a more condensed list with github user names only
and no avatars, one contributor per line linked to issue contributions.
\item{3} "text" for a single line of text containing comma-separated github
user names linked to issue contributions.
}}

\item{alphabetical}{If \code{TRUE}, order contributors alphabetically, otherwise
order by decreasing numbers of contributions.}

\item{open_issue}{If \code{TRUE}, open or edit an issue on github in order to
notify all contributors that they've been added to your \code{README} (see Note).}

\item{force_update}{If \code{TRUE}, update the specified files even if
contributions have not changed.}
}
\value{
Named list of logical values indicating whether files of given names
were updated or not is returned invisibly (that is, only if explicitly
assigned to a return value).
}
\description{
Add contributors to README.Rmd
}
\note{
Opening an issue on github requires the github command-line interface
to be locally installed. See \url{https://cli.github.com/}.
}
\examples{
# The following code extracts the contributors from the git repository
# associated with current working directory and writes them to a file.
\dontrun{
f <- tempfile (fileext = ".Rmd")
writeLines ("", f) # blank file in tempdir()
add_contributors (repo = ".", files = f)
}
}
