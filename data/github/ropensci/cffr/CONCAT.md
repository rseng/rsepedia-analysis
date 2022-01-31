
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cffr <a href='https://docs.ropensci.org/cffr/'><img src="man/figures/logo.png" align="right" height="139"/></a>

<!-- badges: start -->

[![CRAN-status](https://www.r-pkg.org/badges/version/cffr)](https://CRAN.R-project.org/package=cffr)
[![CRAN-results](https://cranchecks.info/badges/worst/cffr)](https://cran.r-project.org/web/checks/check_results_cffr.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/cffr?color=blue)](https://cran.r-project.org/package=cffr)
[![R-CMD-check](https://github.com/ropensci/cffr/actions/workflows/check-full.yaml/badge.svg)](https://github.com/ropensci/cffr/actions/workflows/check-full.yaml)
[![codecov](https://codecov.io/gh/ropensci/cffr/branch/main/graph/badge.svg?token=YRO3XL8RWK)](https://app.codecov.io/gh/ropensci/cffr)
[![r-universe](https://ropensci.r-universe.dev/badges/cffr)](https://ropensci.r-universe.dev/)
[![CITATION-cff](https://github.com/ropensci/cffr/actions/workflows/cff-validator.yml/badge.svg)](https://github.com/ropensci/cffr/actions/workflows/cff-validator.yml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03900/status.svg)](https://doi.org/10.21105/joss.03900)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub code size in
bytes](https://img.shields.io/github/languages/code-size/ropensci/cffr)
[![peer-review](https://badges.ropensci.org/463_status.svg)](https://github.com/ropensci/software-review/issues/463)

<!-- badges: end -->

**cffr** provides utilities to generate, parse, modify and validate
`CITATION.cff` files automatically for **R** packages, as well as tools
and examples for working with .cff more generally.

## What is a `CITATION.cff` file?

[Citation File Format (CFF](https://citation-file-format.github.io/))
(Druskat et al. [2021](#ref-druskat_citation_2021)) (v1.2.0) are plain
text files with human- and machine-readable citation information for
software (and datasets). Code developers can include them in their
repositories to let others know how to correctly cite their software.

This format is becoming popular within the software citation ecosystem.
Recently
[GitHub](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files),
[Zenodo](https://twitter.com/ZENODO_ORG/status/1420357001490706442) and
[Zotero](https://twitter.com/zotero/status/1420515377390530560) have
included full support of this citation format (Druskat
[2021](#ref-druskat_stephan_making_2021)). GitHub support is of special
interest:

<img src="vignettes/tweet-1.png" title="GitHub-link" alt="GitHub-link" width="400" style="display: block; margin: auto;" />

*— Nat Friedman (@natfriedman) [July 27,
2021](https://twitter.com/natfriedman/status/1420122675813441540?ref_src=twsrc%5Etfw)*

See [Enhanced support for citations on
GitHub](https://github.blog/2021-08-19-enhanced-support-citations-github/)
(Smith [2021](#ref-smith2021)) for more info.

### Related projects

[The CodeMeta Project](https://codemeta.github.io/) (Jones et al.
[2017](#ref-jones2017)) creates a concept vocabulary that can be used to
standardize the exchange of software metadata across repositories and
organizations. One of the many uses of a `codemeta.json` file (created
following the standards defined on The CodeMeta Project) is to provide
citation metadata such as title, authors, publication year, and venue
(Fenner [2021](#ref-fenner2021)). The packages
[**codemeta**](https://github.com/cboettig/codemeta)/
[**codemetar**](https://github.com/ropensci/codemetar) allows to
generate `codemeta.json` files from R packages metadata.

## The cffr package

**cffr** maximizes the data extraction by using both the `DESCRIPTION`
file and the `CITATION` file (if present) of your package. Note that
**cffr** works best if your package pass `R CMD
check/devtools::check()`.

As per 2022-01-26 there are at least 92 repos on GitHub using **cffr**.
[Check them out
here](https://github.com/search?l=&o=desc&q=cffr+extension%3Acff+filename%3ACITATION&s=indexed&type=Code).

### Installation

Install **cffr** from [CRAN](https://CRAN.R-project.org/package=cffr):

``` r
install.packages("cffr")
```

You can install the developing version of **cffr** with:

``` r
devtools::install_github("ropensci/cffr")
```

Alternatively, you can install **cffr** using the
[r-universe](https://ropensci.r-universe.dev/ui#builds):

``` r

# Enable this universe
options(repos = c(
  ropensci = "https://ropensci.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

# Install some packages
install.packages("cffr")
```

### Example

By default most often from within your package folder you’ll simply run
`cff_write()`, that creates a `cff` object, write it on a `CITATION.cff`
file and validates it on a single command:

``` r

library(cffr)

# For in-development packages
cff_write()
#>
#> CITATION.cff generated
#>
#> cff_validate results-----
#> Congratulations! This .cff file is valid
```

However, **cffr** provides also custom print methods and mechanisms that
allows you to customize the `CITATION.cff` and integrate them in your
workflows.

This is a basic example which shows you how to create a `cff` object
(see `?cff` for more info). In this case, we are creating a `cff` object
from the metadata of the **rmarkdown** package:

``` r
library(cffr)

# Example with an installed package
test <- cff_create("rmarkdown")
```

<details>

<summary><code>CITATION.cff</code> for <strong>rmarkdown</strong>
</summary>

    cff-version: 1.2.0
    message: 'To cite package "rmarkdown" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'rmarkdown: Dynamic Documents for R'
    version: '2.11'
    abstract: Convert R Markdown documents into a variety of formats.
    authors:
    - family-names: Allaire
      given-names: JJ
      email: jj@rstudio.com
    - family-names: Xie
      given-names: Yihui
      email: xie@yihui.name
      orcid: https://orcid.org/0000-0003-0645-5666
    - family-names: McPherson
      given-names: Jonathan
      email: jonathan@rstudio.com
    - family-names: Luraschi
      given-names: Javier
      email: javier@rstudio.com
    - family-names: Ushey
      given-names: Kevin
      email: kevin@rstudio.com
    - family-names: Atkins
      given-names: Aron
      email: aron@rstudio.com
    - family-names: Wickham
      given-names: Hadley
      email: hadley@rstudio.com
    - family-names: Cheng
      given-names: Joe
      email: joe@rstudio.com
    - family-names: Chang
      given-names: Winston
      email: winston@rstudio.com
    - family-names: Iannone
      given-names: Richard
      email: rich@rstudio.com
      orcid: https://orcid.org/0000-0003-3925-190X
    preferred-citation:
      type: manual
      title: 'rmarkdown: Dynamic Documents for R'
      authors:
      - family-names: Allaire
        given-names: JJ
        email: jj@rstudio.com
      - family-names: Xie
        given-names: Yihui
        email: xie@yihui.name
        orcid: https://orcid.org/0000-0003-0645-5666
      - family-names: McPherson
        given-names: Jonathan
        email: jonathan@rstudio.com
      - family-names: Luraschi
        given-names: Javier
        email: javier@rstudio.com
      - family-names: Ushey
        given-names: Kevin
        email: kevin@rstudio.com
      - family-names: Atkins
        given-names: Aron
        email: aron@rstudio.com
      - family-names: Wickham
        given-names: Hadley
        email: hadley@rstudio.com
      - family-names: Cheng
        given-names: Joe
        email: joe@rstudio.com
      - family-names: Chang
        given-names: Winston
        email: winston@rstudio.com
      - family-names: Iannone
        given-names: Richard
        email: rich@rstudio.com
        orcid: https://orcid.org/0000-0003-3925-190X
      year: '2021'
      notes: R package version 2.11
      url: https://github.com/rstudio/rmarkdown
    repository: https://CRAN.R-project.org/package=rmarkdown
    repository-code: https://github.com/rstudio/rmarkdown
    url: https://pkgs.rstudio.com/rmarkdown/
    date-released: '2021-09-14'
    contact:
    - family-names: Xie
      given-names: Yihui
      email: xie@yihui.name
      orcid: https://orcid.org/0000-0003-0645-5666
    keywords:
    - literate-programming
    - markdown
    - pandoc
    - r
    - r-package
    - rmarkdown
    references:
    - type: book
      title: 'R Markdown: The Definitive Guide'
      authors:
      - family-names: Xie
        given-names: Yihui
      - family-names: Allaire
        given-names: J.J.
      - family-names: Grolemund
        given-names: Garrett
      publisher:
        name: Chapman and Hall/CRC
        address: Boca Raton, Florida
      year: '2018'
      notes: ISBN 9781138359338
      url: https://bookdown.org/yihui/rmarkdown
    - type: book
      title: R Markdown Cookbook
      authors:
      - family-names: Xie
        given-names: Yihui
      - family-names: Dervieux
        given-names: Christophe
      - family-names: Riederer
        given-names: Emily
      publisher:
        name: Chapman and Hall/CRC
        address: Boca Raton, Florida
      year: '2020'
      notes: ISBN 9780367563837
      url: https://bookdown.org/yihui/rmarkdown-cookbook
    - type: software
      title: 'R: A Language and Environment for Statistical Computing'
      notes: Depends
      authors:
      - name: R Core Team
      location:
        name: Vienna, Austria
      year: '2022'
      url: https://www.R-project.org/
      institution:
        name: R Foundation for Statistical Computing
      version: '>= 3.0'
    - type: software
      title: tools
      abstract: 'R: A Language and Environment for Statistical Computing'
      notes: Imports
      authors:
      - name: R Core Team
      location:
        name: Vienna, Austria
      year: '2022'
      url: https://www.R-project.org/
      institution:
        name: R Foundation for Statistical Computing
    - type: software
      title: utils
      abstract: 'R: A Language and Environment for Statistical Computing'
      notes: Imports
      authors:
      - name: R Core Team
      location:
        name: Vienna, Austria
      year: '2022'
      url: https://www.R-project.org/
      institution:
        name: R Foundation for Statistical Computing
    - type: software
      title: knitr
      abstract: 'knitr: A General-Purpose Package for Dynamic Report Generation in R'
      notes: Imports
      authors:
      - family-names: Xie
        given-names: Yihui
        email: xie@yihui.name
        orcid: https://orcid.org/0000-0003-0645-5666
      year: '2022'
      url: https://yihui.org/knitr/
      version: '>= 1.22'
    - type: software
      title: yaml
      abstract: 'yaml: Methods to Convert R Data to YAML and Back'
      notes: Imports
      authors:
      - family-names: Stephens
        given-names: Jeremy
      - family-names: Simonov
        given-names: Kirill
      - family-names: Xie
        given-names: Yihui
      - family-names: Dong
        given-names: Zhuoer
      - family-names: Wickham
        given-names: Hadley
      - family-names: Horner
        given-names: Jeffrey
      - name: reikoch
      - family-names: Beasley
        given-names: Will
      - family-names: O'Connor
        given-names: Brendan
      - family-names: Warnes
        given-names: Gregory R.
      - family-names: Quinn
        given-names: Michael
      - family-names: Kamvar
        given-names: Zhian N.
      year: '2022'
      url: https://CRAN.R-project.org/package=yaml
      version: '>= 2.1.19'
    - type: software
      title: htmltools
      abstract: 'htmltools: Tools for HTML'
      notes: Imports
      authors:
      - family-names: Cheng
        given-names: Joe
        email: joe@rstudio.com
      - family-names: Sievert
        given-names: Carson
        email: carson@rstudio.com
        orcid: https://orcid.org/0000-0002-4958-2844
      - family-names: Schloerke
        given-names: Barret
        email: barret@rstudio.com
        orcid: https://orcid.org/0000-0001-9986-114X
      - family-names: Chang
        given-names: Winston
        email: winston@rstudio.com
        orcid: https://orcid.org/0000-0002-1576-2126
      - family-names: Xie
        given-names: Yihui
        email: yihui@rstudio.com
      - family-names: Allen
        given-names: Jeff
        email: jeff@rstudio.com
      year: '2022'
      url: https://github.com/rstudio/htmltools
      version: '>= 0.3.5'
    - type: software
      title: evaluate
      abstract: 'evaluate: Parsing and Evaluation Tools that Provide More Details than
        the Default'
      notes: Imports
      authors:
      - family-names: Wickham
        given-names: Hadley
      - family-names: Xie
        given-names: Yihui
        email: xie@yihui.name
        orcid: https://orcid.org/0000-0003-0645-5666
      year: '2022'
      url: https://github.com/r-lib/evaluate
      version: '>= 0.13'
    - type: software
      title: jsonlite
      abstract: 'jsonlite: A Simple and Robust JSON Parser and Generator for R'
      notes: Imports
      authors:
      - family-names: Ooms
        given-names: Jeroen
        email: jeroen@berkeley.edu
        orcid: https://orcid.org/0000-0002-4035-0289
      year: '2022'
    - type: software
      title: tinytex
      abstract: 'tinytex: Helper Functions to Install and Maintain TeX Live, and Compile
        LaTeX Documents'
      notes: Imports
      authors:
      - family-names: Xie
        given-names: Yihui
        email: xie@yihui.name
        orcid: https://orcid.org/0000-0003-0645-5666
      year: '2022'
      url: https://github.com/yihui/tinytex
      version: '>= 0.31'
    - type: software
      title: xfun
      abstract: 'xfun: Supporting Functions for Packages Maintained by ''Yihui Xie'''
      notes: Imports
      authors:
      - family-names: Xie
        given-names: Yihui
        email: xie@yihui.name
        orcid: https://orcid.org/0000-0003-0645-5666
      year: '2022'
      url: https://github.com/yihui/xfun
      version: '>= 0.21'
    - type: software
      title: jquerylib
      abstract: 'jquerylib: Obtain ''jQuery'' as an HTML Dependency Object'
      notes: Imports
      authors:
      - family-names: Sievert
        given-names: Carson
        email: carson@rstudio.com
        orcid: https://orcid.org/0000-0002-4958-2844
      - family-names: Cheng
        given-names: Joe
        email: joe@rstudio.com
      year: '2022'
    - type: software
      title: methods
      abstract: 'R: A Language and Environment for Statistical Computing'
      notes: Imports
      authors:
      - name: R Core Team
      location:
        name: Vienna, Austria
      year: '2022'
      url: https://www.R-project.org/
      institution:
        name: R Foundation for Statistical Computing
    - type: software
      title: stringr
      abstract: 'stringr: Simple, Consistent Wrappers for Common String Operations'
      notes: Imports
      authors:
      - family-names: Wickham
        given-names: Hadley
        email: hadley@rstudio.com
      year: '2022'
      version: '>= 1.2.0'
    - type: software
      title: testthat
      abstract: 'testthat: Unit Testing for R'
      notes: Suggests
      authors:
      - family-names: Wickham
        given-names: Hadley
        email: hadley@rstudio.com
      year: '2022'
      version: '>= 3.0.0'
    - type: software
      title: digest
      abstract: 'digest: Create Compact Hash Digests of R Objects'
      notes: Suggests
      authors:
      - family-names: Lucas
        given-names: Dirk Eddelbuettel with contributions by Antoine
        email: edd@debian.org
      - family-names: Tuszynski
        given-names: Jarek
      - family-names: Bengtsson
        given-names: Henrik
      - family-names: Urbanek
        given-names: Simon
      - family-names: Frasca
        given-names: Mario
      - family-names: Lewis
        given-names: Bryan
      - family-names: Stokely
        given-names: Murray
      - family-names: Muehleisen
        given-names: Hannes
      - family-names: Murdoch
        given-names: Duncan
      - family-names: Hester
        given-names: Jim
      - family-names: Wu
        given-names: Wush
      - family-names: Kou
        given-names: Qiang
      - family-names: Onkelinx
        given-names: Thierry
      - family-names: Lang
        given-names: Michel
      - family-names: Simko
        given-names: Viliam
      - family-names: Hornik
        given-names: Kurt
      - family-names: Neal
        given-names: Radford
      - family-names: Bell
        given-names: Kendon
      - family-names: de Queljoe
        given-names: Matthew
      - family-names: Suruceanu
        given-names: Ion
      - family-names: Denney
        given-names: Bill
      - family-names: Schumacher
        given-names: Dirk
      - family-names: Chang.
        given-names: and Winston
      year: '2022'
    - type: software
      title: vctrs
      abstract: 'vctrs: Vector Helpers'
      notes: Suggests
      authors:
      - family-names: Wickham
        given-names: Hadley
        email: hadley@rstudio.com
      - family-names: Henry
        given-names: Lionel
        email: lionel@rstudio.com
      - family-names: Vaughan
        given-names: Davis
        email: davis@rstudio.com
      year: '2022'
      url: https://vctrs.r-lib.org/
    - type: software
      title: tibble
      abstract: 'tibble: Simple Data Frames'
      notes: Suggests
      authors:
      - family-names: Müller
        given-names: Kirill
        email: krlmlr+r@mailbox.org
      - family-names: Wickham
        given-names: Hadley
        email: hadley@rstudio.com
      year: '2022'
    - type: software
      title: fs
      abstract: 'fs: Cross-Platform File System Operations Based on ''libuv'''
      notes: Suggests
      authors:
      - family-names: Hester
        given-names: Jim
      - family-names: Wickham
        given-names: Hadley
        email: hadley@rstudio.com
      - family-names: Csárdi
        given-names: Gábor
        email: csardi.gabor@gmail.com
      year: '2022'
    - type: software
      title: withr
      abstract: 'withr: Run Code ''With'' Temporarily Modified Global State'
      notes: Suggests
      authors:
      - family-names: Hester
        given-names: Jim
      - family-names: Henry
        given-names: Lionel
        email: lionel@rstudio.com
      - family-names: Müller
        given-names: Kirill
        email: krlmlr+r@mailbox.org
      - family-names: Ushey
        given-names: Kevin
        email: kevinushey@gmail.com
      - family-names: Wickham
        given-names: Hadley
        email: hadley@rstudio.com
      - family-names: Chang
        given-names: Winston
      year: '2022'
      version: '>= 2.4.2'
    - type: software
      title: bslib
      abstract: 'bslib: Custom ''Bootstrap'' ''Sass'' Themes for ''shiny'' and ''rmarkdown'''
      notes: Suggests
      authors:
      - family-names: Sievert
        given-names: Carson
        email: carson@rstudio.com
        orcid: https://orcid.org/0000-0002-4958-2844
      - family-names: Cheng
        given-names: Joe
        email: joe@rstudio.com
      year: '2022'
      version: '>= 0.2.5.1'
    - type: software
      title: sass
      abstract: 'sass: Syntactically Awesome Style Sheets (''Sass'')'
      notes: Suggests
      authors:
      - family-names: Cheng
        given-names: Joe
        email: joe@rstudio.com
      - family-names: Mastny
        given-names: Timothy
        email: tim.mastny@gmail.com
      - family-names: Iannone
        given-names: Richard
        email: rich@rstudio.com
        orcid: https://orcid.org/0000-0003-3925-190X
      - family-names: Schloerke
        given-names: Barret
        email: barret@rstudio.com
        orcid: https://orcid.org/0000-0001-9986-114X
      - family-names: Sievert
        given-names: Carson
        email: carson@rstudio.com
        orcid: https://orcid.org/0000-0002-4958-2844
      year: '2022'
      url: https://github.com/rstudio/sass
      version: '>= 0.4.0'

</details>

<p>

We can validate the result using `cff_validate()`:

``` r

cff_validate(test)
#> 
#> cff_validate results-----
#> Congratulations! This cff object is valid
```

Check the [docs](https://docs.ropensci.org/cffr/reference/index.html)
and `vignette("cffr", package = "cffr")` to learn how to work with `cff`
objects.

### Keep your `CITATION.cff` file up-to-date

#### GitHub Actions

The easiest way for keeping you `CITATION.cff` file up-to-date is using
GitHub Actions. Use `cff_gha_update()`function to install a GitHub
Action that would update your `CITATION.cff` file on the following
events:

  - When you publish a new release of the package on your GitHub repo.
  - Each time that you modify your DESCRIPTION or inst/CITATION files.
  - The action can be run also manually.

<!-- end list -->

``` r
cff_gha_update()

#> Installing update-citation-cff.yaml on './.github/workflows'
#> Adding .github to .Rbuildignore
```

See the example workflow file
[here](https://github.com/ropensci/cffr/blob/main/.github/workflows/update-citation-cff.yaml).

#### Git pre-commit hook [![Experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

You can also use a [git pre-commit
hook](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks#_committing_workflow_hooks):

> The `pre-commit` hook is run first, before you even type in a commit
> message. It’s used to inspect the snapshot that’s about to be
> committed, to see if you’ve forgotten something, to make sure tests
> run, or to examine whatever you need to inspect in the code. Exiting
> non-zero from this hook aborts the commit, although you can bypass it
> with `git commit --no-verify`.

A specific pre-commit hook can be installed with
`cff_git_hook_install()`. If you want to use a pre-commit hook, please
make sure you have the **testthat** package installed.

### Learn more

Check the following articles to learn more about **cffr**:

  - [cffr: Create a CITATION.cff File for your R
    Package](https://ropensci.org/blog/2021/11/23/cffr/)
  - [How I Test cffr on (about) 2,000 Packages using GitHub Actions and
    R-universe](https://ropensci.org/blog/2021/11/23/how-i-test-cffr/)

## Related packages

  - [**citation**](https://github.com/pik-piam/citation/): The
    development version (at the time of this writing) includes a new
    function `r2cff` that creates a `CITATION.cff` file (v1.1.0) using
    the information of your `DESCRIPTION` file. It also provide minimal
    validity checks.
  - [**handlr**](https://github.com/ropensci/handlr): Tool for
    converting among citation formats, including `*.cff` files. At the
    time of this writing only CFF v1.1.0 was supported (see
    [\#24](https://github.com/ropensci/handlr/issues/24)).
  - [**codemeta**](https://github.com/cboettig/codemeta)/
    [**codemetar**](https://github.com/ropensci/codemetar) provides
    similar solutions for creating `codemeta.json` file, another format
    for storing and sharing software metadata.

## Citation

Hernangómez D (2021). “cffr: Generate Citation File Format Metadata for
R Packages.” *Journal of Open Source Software*, *6*(67), 3900. doi:
10.21105/joss.03900 (URL: <https://doi.org/10.21105/joss.03900>), \<URL:
<https://doi.org/10.21105/joss.03900>\>.

A BibTeX entry for LaTeX users is

    @Article{hernangomez2021,
      doi = {10.21105/joss.03900},
      url = {https://doi.org/10.21105/joss.03900},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {67},
      pages = {3900},
      author = {Diego Hernangómez},
      title = {cffr: Generate Citation File Format Metadata for R Packages},
      journal = {Journal of Open Source Software},
    }

You can also use the [citation provided by
GitHub](https://github.com/ropensci/cffr), that is generated from the
information of a `CITATION.cff` created with **cffr**. See [About
CITATION
files](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)
for more info.

## References

<div id="refs">

<div id="ref-druskat_stephan_making_2021">

Druskat, Stephan. 2021. “Making Software Citation Easi(er) - the
Citation File Format and Its Integrations.”
<https://doi.org/10.5281/zenodo.5529914>.

</div>

<div id="ref-druskat_citation_2021">

Druskat, Stephan, Jurriaan H. Spaaks, Neil Chue Hong, Robert Haines,
James Baker, Spencer Bliven, Egon Willighagen, David Pérez-Suárez, and
Alexander Konovalov. 2021. “Citation File Format.”
<https://doi.org/10.5281/zenodo.5171937>.

</div>

<div id="ref-fenner2021">

Fenner, Martin. 2021. “We Need Your Feedback: Aligning the CodeMeta
Vocabulary for Scientific Software with Schema.org.”
<https://doi.org/10.5438/a49j-x692>.

</div>

<div id="ref-jones2017">

Jones, Matthew B, Carl Boettiger, Abby Cabunoc Mayes, Arfon Smith, Peter
Slaughter, Kyle Niemeyer, Yolanda Gil, et al. 2017. *CodeMeta: An
Exchange Schema for Software Metadata*. KNB Data Repository.
<https://doi.org/10.5063/SCHEMA/CODEMETA-2.0>.

</div>

<div id="ref-smith2021">

Smith, Arfon. 2021. “Enhanced Support for Citations on GitHub.”
<https://github.blog/2021-08-19-enhanced-support-citations-github/>.

</div>

</div>

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# cffr 0.2.1

-   GitHub Action now runs only on `master` or `main`branch.

-   Better handling of references

# cffr 0.2.0

-   Now **cffr** extracts also information of the package dependencies and adds
    the main citation of the dependencies to the `references` field, using
    `citation(auto = TRUE)`.

    -   New `dependencies` parameter on `cff_create()` and `cff_write()`.

-   Other improvements on `cff_parse_citation():`

    -   `cff_parse_citation()` extracts more information of authors, based on
        the fields provided on the DESCRIPTION file.

    -   `cff_parse_citation()` does a better job extracting information from
        `bibentry()` /BibTeX and mapping it to `preferred-citation/references`
        fields of CFF.

-   Add new functions for working with git pre-commit hooks
    [![Experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental):

    -   `cff_git_hook_install()`
    -   `cff_git_hook_remove()`

-   New BibTeX functions:

    -   `cff_extract_to_bibtex()`
    -   `cff_to_bibtex()`
    -   `cff_parse_person_bibtex()`
    -   `write_bib()`

-   Add a new dependency: `lifecycle`.

# cffr 0.1.1

-   Accepted on JOSS
    [![DOI](https://joss.theoj.org/papers/10.21105/joss.03900/status.svg)](https://doi.org/10.21105/joss.03900)
-   Include `pages` on `cff_parse_citation()` .
-   New `gh_keywords` parameter on `cff_create()` /`cff_write()`. If `TRUE`, it
    would include GitHub repo topics as keywords.

# cffr 0.1.0

-   First CRAN release

# cffr 0.0.2

-   `cffr` is part now of rOpenSci.
-   Update on docs and README.
-   Add fuzzy match on `keys` parameter.
-   New dataset: `cran_to_spdx`.
-   Add DOI <https://doi.org/10.5281/zenodo.5509766>
-   Citation of installed packages extracted using `citation().`
-   Auto-generating `preferred-citation` key from DESCRIPTION.
-   Rename `cff_schema_definitions_reference()` to
    `cff_schema_definitions_refs()`.
-   "repository" key is supported.
-   Added vignette: `vignette("crosswalk", package = "cffr")`.
-   Add support to Bioconductor packages.
-   New function: `cff_gha_update()`.

# cffr 0.0.1

-   First stable release
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

Please note that the cffr project is released with a
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
# Test on local installation

**This folder is `.Rbuildignored`**.

This test validates the `cff` parsing for \>1500 packages:

-   Core packages of every [CRAN Task
    Views](https://cran.r-project.org/web/views/) and their dependencies.
-   All the packages available in the [rOpenSci
    r-universe](https://ropensci.r-universe.dev/) and their dependencies.
-   All the packages of r-forge, r-lib and r-studio: lists extracted from
    https://r-universe.dev/organizations/.

This test is deployed in [GitHub
Actions](https://github.com/ropensci/cffr/actions/workflows/test-ci.yaml) and
the results are uploaded as an artifact. We use here Windows and MacOS binaries
for speeding up the process.

As the installations differs across users and machines, the snapshot testing is
expected to fail on a normal run. However, the snapshots are quite useful for
extensive tests and debugging, as well as for capturing corner cases. For that
reason, these tests are no run in **CRAN** or in the regular package development
workflow.

However, the test can be run locally with

``` r
# Load package
devtools::load_all()

# Run the tests
testthat::test_dir("tests/testthat/test_ci")
```
# Parse date

    cff-version: 1.2.0
    message: 'To cite package "rgeos" in publications use:'
    type: software
    license: GPL-2.0-or-later
    title: 'rgeos: Interface to Geometry Engine - Open Source (''GEOS'')'
    version: 0.5-7
    abstract: 'Interface to Geometry Engine - Open Source (''GEOS'') using the C ''API''
      for topology operations on geometries. Please note that ''rgeos'' will be retired
      by the end of 2023, plan transition to sf functions using ''GEOS'' at your earliest
      convenience. The ''GEOS'' library is external to the package, and, when installing
      the package from source, must be correctly installed first. Windows and Mac Intel
      OS X binaries are provided on ''CRAN''. (''rgeos'' >= 0.5-1): Up to and including
      ''GEOS'' 3.7.1, topological operations succeeded with some invalid geometries for
      which the same operations fail from and including ''GEOS'' 3.7.2. The ''checkValidity=''
      argument defaults and structure have been changed, from default FALSE to integer
      default ''0L'' for ''GEOS'' < 3.7.2 (no check), ''1L'' ''GEOS'' >= 3.7.2 (check
      and warn). A value of ''2L'' is also provided that may be used, assigned globally
      using ''set_RGEOS_CheckValidity(2L)'', or locally using the ''checkValidity=2L''
      argument, to attempt zero-width buffer repair if invalid geometries are found. The
      previous default (FALSE, now ''0L'') is fastest and used for ''GEOS'' < 3.7.2, but
      will not warn users of possible problems before the failure of topological operations
      that previously succeeded. From ''GEOS'' 3.8.0, repair of geometries may also be
      attempted using ''gMakeValid()'', which may, however, return a collection of geometries
      of different types.'
    authors:
    - family-names: Bivand
      given-names: Roger
      email: Roger.Bivand@nhh.no
      orcid: https://orcid.org/0000-0003-2392-6140
    - family-names: Rundel
      given-names: Colin
    preferred-citation:
      type: manual
      title: 'rgeos: Interface to Geometry Engine - Open Source (''GEOS'')'
      authors:
      - family-names: Bivand
        given-names: Roger
        email: Roger.Bivand@nhh.no
        orcid: https://orcid.org/0000-0003-2392-6140
      - family-names: Rundel
        given-names: Colin
      version: 0.5-7
      abstract: 'Interface to Geometry Engine - Open Source (''GEOS'') using the C ''API''
        for topology operations on geometries. Please note that ''rgeos'' will be retired
        by the end of 2023, plan transition to sf functions using ''GEOS'' at your earliest
        convenience. The ''GEOS'' library is external to the package, and, when installing
        the package from source, must be correctly installed first. Windows and Mac Intel
        OS X binaries are provided on ''CRAN''. (''rgeos'' >= 0.5-1): Up to and including
        ''GEOS'' 3.7.1, topological operations succeeded with some invalid geometries
        for which the same operations fail from and including ''GEOS'' 3.7.2. The ''checkValidity=''
        argument defaults and structure have been changed, from default FALSE to integer
        default ''0L'' for ''GEOS'' < 3.7.2 (no check), ''1L'' ''GEOS'' >= 3.7.2 (check
        and warn). A value of ''2L'' is also provided that may be used, assigned globally
        using ''set_RGEOS_CheckValidity(2L)'', or locally using the ''checkValidity=2L''
        argument, to attempt zero-width buffer repair if invalid geometries are found.
        The previous default (FALSE, now ''0L'') is fastest and used for ''GEOS'' < 3.7.2,
        but will not warn users of possible problems before the failure of topological
        operations that previously succeeded. From ''GEOS'' 3.8.0, repair of geometries
        may also be attempted using ''gMakeValid()'', which may, however, return a collection
        of geometries of different types.'
      repository: https://CRAN.R-project.org/package=rgeos
      repository-code: https://r-forge.r-project.org/projects/rgeos/
      url: https://trac.osgeo.org/geos/
      identifiers:
      - type: url
        value: http://rgeos.r-forge.r-project.org/index.html
      date-released: '2020-09-07'
      contact:
      - family-names: Bivand
        given-names: Roger
        email: Roger.Bivand@nhh.no
        orcid: https://orcid.org/0000-0003-2392-6140
      license: GPL-2.0-or-later
      year: '2020'
    repository: https://CRAN.R-project.org/package=rgeos
    repository-code: https://r-forge.r-project.org/projects/rgeos/
    url: https://trac.osgeo.org/geos/
    date-released: '2020-09-07'
    contact:
    - family-names: Bivand
      given-names: Roger
      email: Roger.Bivand@nhh.no
      orcid: https://orcid.org/0000-0003-2392-6140
    identifiers:
    - type: url
      value: http://rgeos.r-forge.r-project.org/index.html

# Parse date in another format

    cff-version: 1.2.0
    message: 'To cite package "basicdescdate" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'basicdescdate: A Basic Description with Date'
    version: 0.1.6
    abstract: A very basic description. Should parse without problems. I have a Date
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    preferred-citation:
      type: manual
      title: 'basicdescdate: A Basic Description with Date'
      authors:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      version: 0.1.6
      abstract: A very basic description. Should parse without problems. I have a Date
      repository-code: https://github.com/basic/package
      url: https://basic.github.io/package
      date-released: '1999-01-01'
      contact:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      license: GPL-3.0-only
      year: '1999'
    repository-code: https://github.com/basic/package
    url: https://basic.github.io/package
    date-released: '1999-01-01'
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com

# Parsing many urls

    cff-version: 1.2.0
    message: 'To cite package "manyurls" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'manyurls: A lot of urls'
    version: 0.1.6
    abstract: This package has many urls. Specifically, 1 Bug Reports and 6 URLs. Expected
      is to have 1 repository-code, 1 url and 3 URLs, since there is 1 duplicate and 1
      invalid url.
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    preferred-citation:
      type: manual
      title: 'manyurls: A lot of urls'
      authors:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      version: 0.1.6
      abstract: This package has many urls. Specifically, 1 Bug Reports and 6 URLs. Expected
        is to have 1 repository-code, 1 url and 3 URLs, since there is 1 duplicate and
        1 invalid url.
      repository-code: https://github.com/test/package
      url: https://test.github.io/package/
      identifiers:
      - type: url
        value: https://r-forge.r-project.org/projects/test/
      - type: url
        value: http://google.ru
      - type: url
        value: https://gitlab.com/r-packages/behaviorchange
      contact:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      license: GPL-3.0-only
      year: '2022'
    repository-code: https://github.com/test/package
    url: https://test.github.io/package/
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    identifiers:
    - type: url
      value: https://r-forge.r-project.org/projects/test/
    - type: url
      value: http://google.ru
    - type: url
      value: https://gitlab.com/r-packages/behaviorchange

# Parsing Gitlab

    cff-version: 1.2.0
    message: 'To cite package "codemetar" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
    version: 0.1.6
    abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
      metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
      to generate, parse, and modify 'codemeta.json' files automatically for R packages,
      as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more generally.
    authors:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    - family-names: Salmon
      given-names: Maëlle
      orcid: https://orcid.org/0000-0002-2815-0399
    preferred-citation:
      type: manual
      title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
      authors:
      - family-names: Boettiger
        given-names: Carl
        email: cboettig@gmail.com
        orcid: https://orcid.org/0000-0002-1642-628X
      - family-names: Salmon
        given-names: Maëlle
        orcid: https://orcid.org/0000-0002-2815-0399
      version: 0.1.6
      abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
        metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
        to generate, parse, and modify 'codemeta.json' files automatically for R packages,
        as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more
        generally.
      repository: https://CRAN.R-project.org/package=codemetar
      repository-code: https://gitlab.com/ninijay/methoden
      url: https://ropensci.github.io/codemetar
      contact:
      - family-names: Boettiger
        given-names: Carl
        email: cboettig@gmail.com
        orcid: https://orcid.org/0000-0002-1642-628X
      keywords:
      - metadata
      - codemeta
      - ropensci
      - citation
      - credit
      - linked-data
      license: GPL-3.0-only
      year: '2022'
    repository: https://CRAN.R-project.org/package=codemetar
    repository-code: https://gitlab.com/ninijay/methoden
    url: https://ropensci.github.io/codemetar
    contact:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    keywords:
    - metadata
    - codemeta
    - ropensci
    - citation
    - credit
    - linked-data

# Parsing many persons

    cff-version: 1.2.0
    message: 'To cite package "manypersons" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'manypersons: A lot of persons'
    version: 0.1.6
    abstract: Overkill desc with many persons. Try this
    authors:
    - family-names: Hernangómez
      given-names: Diego
      email: fake@gmail.com
      orcid: https://orcid.org/0000-0001-8457-4658
    - family-names: Doe
      given-names: Joe
      affiliation: This One
      country: ES
    - family-names: Doe
      given-names: Pepe
      email: fake@gmail.com
    - name: I am an entity
      date-end: '2020-01-01'
    preferred-citation:
      type: manual
      title: 'manypersons: A lot of persons'
      authors:
      - family-names: Hernangómez
        given-names: Diego
        email: fake@gmail.com
        orcid: https://orcid.org/0000-0001-8457-4658
      - family-names: Doe
        given-names: Joe
        affiliation: This One
        country: ES
      - family-names: Doe
        given-names: Pepe
        email: fake@gmail.com
      - name: I am an entity
        date-end: '2020-01-01'
      version: 0.1.6
      abstract: Overkill desc with many persons. Try this
      repository-code: https://github.com/many/persons
      url: https://many.github.io/persons
      contact:
      - family-names: Hernangómez
        given-names: Diego
        email: fake@gmail.com
        orcid: https://orcid.org/0000-0001-8457-4658
      - name: I am an entity
        date-end: '2020-01-01'
      keywords:
      - metadata
      - cffr
      - ropensci
      - citation
      - credit
      - linked-data
      - one
      - two
      license: GPL-3.0-only
      year: '2022'
    repository-code: https://github.com/many/persons
    url: https://many.github.io/persons
    contact:
    - family-names: Hernangómez
      given-names: Diego
      email: fake@gmail.com
      orcid: https://orcid.org/0000-0001-8457-4658
    - name: I am an entity
      date-end: '2020-01-01'
    keywords:
    - metadata
    - cffr
    - ropensci
    - citation
    - credit
    - linked-data
    - one
    - two

# Parsing wrong urls

    cff-version: 1.2.0
    message: 'To cite package "wrongurls" in publications use:'
    type: software
    license: MIT
    title: 'wrongurls: Generate CodeMeta Metadata for R Packages'
    version: 0.1.0
    abstract: Codemeta defines a 'JSON-LD' format for describing software metadata. This
      package provides utilities to generate, parse, and modify codemeta.jsonld files
      automatically for R packages.
    authors:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    preferred-citation:
      type: manual
      title: 'wrongurls: Generate CodeMeta Metadata for R Packages'
      authors:
      - family-names: Boettiger
        given-names: Carl
        email: cboettig@gmail.com
        orcid: https://orcid.org/0000-0002-1642-628X
      version: 0.1.0
      abstract: Codemeta defines a 'JSON-LD' format for describing software metadata.
        This package provides utilities to generate, parse, and modify codemeta.jsonld
        files automatically for R packages.
      url: https://httpbin.org/status/404
      identifiers:
      - type: url
        value: https://httpbin.org/status/429
      - type: url
        value: https://www.github.es/ropensci/codemeta
      contact:
      - family-names: Boettiger
        given-names: Carl
        email: cboettig@gmail.com
        orcid: https://orcid.org/0000-0002-1642-628X
      keywords:
      - metadata
      - codemeta
      - ropensci
      - citation
      - credit
      - linked-data
      license: MIT
      year: '2022'
    url: https://httpbin.org/status/404
    contact:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    keywords:
    - metadata
    - codemeta
    - ropensci
    - citation
    - credit
    - linked-data
    identifiers:
    - type: url
      value: https://httpbin.org/status/429
    - type: url
      value: https://www.github.es/ropensci/codemeta

# Parsing two maintainers

    cff-version: 1.2.0
    message: 'To cite package "codemetar" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
    version: 0.1.6
    abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
      metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
      to generate, parse, and modify 'codemeta.json' files automatically for R packages,
      as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more generally.
    authors:
    - family-names: Ok
      given-names: John
      email: email@email.edu
    - family-names: Doe
      given-names: Jane
      email: email2@email.edu
    - family-names: Doo
      given-names: Jane
    preferred-citation:
      type: manual
      title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
      authors:
      - family-names: Ok
        given-names: John
        email: email@email.edu
      - family-names: Doe
        given-names: Jane
        email: email2@email.edu
      - family-names: Doo
        given-names: Jane
      version: 0.1.6
      abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
        metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
        to generate, parse, and modify 'codemeta.json' files automatically for R packages,
        as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more
        generally.
      repository: https://CRAN.R-project.org/package=codemetar
      repository-code: https://github.com/ropensci/codemetar
      url: https://ropensci.github.io/codemetar
      contact:
      - family-names: Ok
        given-names: John
        email: email@email.edu
      - family-names: Doe
        given-names: Jane
        email: email2@email.edu
      license: GPL-3.0-only
      year: '2022'
    repository: https://CRAN.R-project.org/package=codemetar
    repository-code: https://github.com/ropensci/codemetar
    url: https://ropensci.github.io/codemetar
    contact:
    - family-names: Ok
      given-names: John
      email: email@email.edu
    - family-names: Doe
      given-names: Jane
      email: email2@email.edu

# Parsing r-universe

    cff-version: 1.2.0
    message: 'To cite package "codemetar" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
    version: 0.3.2
    abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
      metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
      to generate, parse, and modify 'codemeta.json' files automatically for R packages,
      as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more generally.
    authors:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    - family-names: Salmon
      given-names: Maëlle
      orcid: https://orcid.org/0000-0002-2815-0399
    preferred-citation:
      type: manual
      title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
      authors:
      - family-names: Boettiger
        given-names: Carl
        email: cboettig@gmail.com
        orcid: https://orcid.org/0000-0002-1642-628X
      - family-names: Salmon
        given-names: Maëlle
        orcid: https://orcid.org/0000-0002-2815-0399
      version: 0.3.2
      abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
        metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
        to generate, parse, and modify 'codemeta.json' files automatically for R packages,
        as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more
        generally.
      repository: https://ropensci.r-universe.dev
      repository-code: https://github.com/ropensci/codemetar
      url: https://docs.ropensci.org/codemetar/
      contact:
      - family-names: Boettiger
        given-names: Carl
        email: cboettig@gmail.com
        orcid: https://orcid.org/0000-0002-1642-628X
      keywords:
      - metadata
      - codemeta
      - ropensci
      - citation
      - credit
      - linked-data
      license: GPL-3.0-only
      year: '2022'
    repository: https://ropensci.r-universe.dev
    repository-code: https://github.com/ropensci/codemetar
    url: https://docs.ropensci.org/codemetar/
    contact:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    keywords:
    - metadata
    - codemeta
    - ropensci
    - citation
    - credit
    - linked-data

# Parsing Bioconductor

    cff-version: 1.2.0
    message: 'To cite package "GenomicRanges" in publications use:'
    type: software
    license: Artistic-2.0
    title: 'GenomicRanges: Representation and manipulation of genomic intervals'
    version: 1.44.0
    abstract: The ability to efficiently represent and manipulate genomic annotations
      and alignments is playing a central role when it comes to analyzing high-throughput
      sequencing data (a.k.a. NGS data). The GenomicRanges package defines general purpose
      containers for storing and manipulating genomic intervals and variables defined
      along a genome. More specialized containers for representing and manipulating short
      alignments against a reference genome, or a matrix-like summarization of an experiment,
      are defined in the GenomicAlignments and SummarizedExperiment packages, respectively.
      Both packages build on top of the GenomicRanges infrastructure.
    authors:
    - name: Bioconductor Package Maintainer
      email: maintainer@bioconductor.org
    - family-names: Aboyoun
      given-names: P.
    - family-names: Pagès
      given-names: H.
    - family-names: Lawrence
      given-names: M.
    preferred-citation:
      type: manual
      title: 'GenomicRanges: Representation and manipulation of genomic intervals'
      authors:
      - name: Bioconductor Package Maintainer
        email: maintainer@bioconductor.org
      - family-names: Aboyoun
        given-names: P.
      - family-names: Pagès
        given-names: H.
      - family-names: Lawrence
        given-names: M.
      version: 1.44.0
      abstract: The ability to efficiently represent and manipulate genomic annotations
        and alignments is playing a central role when it comes to analyzing high-throughput
        sequencing data (a.k.a. NGS data). The GenomicRanges package defines general purpose
        containers for storing and manipulating genomic intervals and variables defined
        along a genome. More specialized containers for representing and manipulating
        short alignments against a reference genome, or a matrix-like summarization of
        an experiment, are defined in the GenomicAlignments and SummarizedExperiment packages,
        respectively. Both packages build on top of the GenomicRanges infrastructure.
      repository: https://bioconductor.org/
      repository-code: https://github.com/Bioconductor/GenomicRanges
      url: https://bioconductor.org/packages/GenomicRanges
      date-released: '2021-05-19'
      contact:
      - name: Bioconductor Package Maintainer
        email: maintainer@bioconductor.org
      license: Artistic-2.0
      year: '2021'
    repository: https://bioconductor.org/
    repository-code: https://github.com/Bioconductor/GenomicRanges
    url: https://bioconductor.org/packages/GenomicRanges
    date-released: '2021-05-19'
    contact:
    - name: Bioconductor Package Maintainer
      email: maintainer@bioconductor.org

# Search package on CRAN

    cff-version: 1.2.0
    message: 'To cite package "ggplot2" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'ggplot2: A Basic Description'
    version: 0.1.6
    abstract: A very basic description. Should parse without problems.
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    preferred-citation:
      type: manual
      title: 'ggplot2: A Basic Description'
      authors:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      version: 0.1.6
      abstract: A very basic description. Should parse without problems.
      repository: https://CRAN.R-project.org/package=ggplot2
      repository-code: https://github.com/basic/package
      url: https://basic.github.io/package
      contact:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      license: GPL-3.0-only
      year: '2022'
    repository: https://CRAN.R-project.org/package=ggplot2
    repository-code: https://github.com/basic/package
    url: https://basic.github.io/package
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com

# Keys snapshot

    Code
      cff_schema_keys(FALSE)
    Output
       [1] "cff-version"         "message"             "type"               
       [4] "license"             "title"               "version"            
       [7] "doi"                 "abstract"            "authors"            
      [10] "preferred-citation"  "repository"          "repository-artifact"
      [13] "repository-code"     "url"                 "date-released"      
      [16] "contact"             "keywords"            "references"         
      [19] "commit"              "identifiers"         "license-url"        

# Test full with CITATION and (option = author)

    cff-version: 1.2.0
    message: 'To cite package "rgeos" in publications use:'
    type: software
    license: GPL-2.0-or-later
    title: 'rgeos: Interface to Geometry Engine - Open Source (''GEOS'')'
    version: 0.5-7
    abstract: 'Interface to Geometry Engine - Open Source (''GEOS'') using the C ''API''
      for topology operations on geometries. Please note that ''rgeos'' will be retired
      by the end of 2023, plan transition to sf functions using ''GEOS'' at your earliest
      convenience. The ''GEOS'' library is external to the package, and, when installing
      the package from source, must be correctly installed first. Windows and Mac Intel
      OS X binaries are provided on ''CRAN''. (''rgeos'' >= 0.5-1): Up to and including
      ''GEOS'' 3.7.1, topological operations succeeded with some invalid geometries for
      which the same operations fail from and including ''GEOS'' 3.7.2. The ''checkValidity=''
      argument defaults and structure have been changed, from default FALSE to integer
      default ''0L'' for ''GEOS'' < 3.7.2 (no check), ''1L'' ''GEOS'' >= 3.7.2 (check
      and warn). A value of ''2L'' is also provided that may be used, assigned globally
      using ''set_RGEOS_CheckValidity(2L)'', or locally using the ''checkValidity=2L''
      argument, to attempt zero-width buffer repair if invalid geometries are found. The
      previous default (FALSE, now ''0L'') is fastest and used for ''GEOS'' < 3.7.2, but
      will not warn users of possible problems before the failure of topological operations
      that previously succeeded. From ''GEOS'' 3.8.0, repair of geometries may also be
      attempted using ''gMakeValid()'', which may, however, return a collection of geometries
      of different types.'
    authors:
    - family-names: Bivand
      given-names: Roger
      email: Roger.Bivand@nhh.no
      orcid: https://orcid.org/0000-0003-2392-6140
    - family-names: Rundel
      given-names: Colin
    preferred-citation:
      type: manual
      title: 'rgeos: Interface to Geometry Engine - Open Source (''GEOS'')'
      authors:
      - family-names: Bivand
        given-names: Roger
        email: Roger.Bivand@nhh.no
        orcid: https://orcid.org/0000-0003-2392-6140
      - family-names: Rundel
        given-names: Colin
      version: 0.5-7
      abstract: 'Interface to Geometry Engine - Open Source (''GEOS'') using the C ''API''
        for topology operations on geometries. Please note that ''rgeos'' will be retired
        by the end of 2023, plan transition to sf functions using ''GEOS'' at your earliest
        convenience. The ''GEOS'' library is external to the package, and, when installing
        the package from source, must be correctly installed first. Windows and Mac Intel
        OS X binaries are provided on ''CRAN''. (''rgeos'' >= 0.5-1): Up to and including
        ''GEOS'' 3.7.1, topological operations succeeded with some invalid geometries
        for which the same operations fail from and including ''GEOS'' 3.7.2. The ''checkValidity=''
        argument defaults and structure have been changed, from default FALSE to integer
        default ''0L'' for ''GEOS'' < 3.7.2 (no check), ''1L'' ''GEOS'' >= 3.7.2 (check
        and warn). A value of ''2L'' is also provided that may be used, assigned globally
        using ''set_RGEOS_CheckValidity(2L)'', or locally using the ''checkValidity=2L''
        argument, to attempt zero-width buffer repair if invalid geometries are found.
        The previous default (FALSE, now ''0L'') is fastest and used for ''GEOS'' < 3.7.2,
        but will not warn users of possible problems before the failure of topological
        operations that previously succeeded. From ''GEOS'' 3.8.0, repair of geometries
        may also be attempted using ''gMakeValid()'', which may, however, return a collection
        of geometries of different types.'
      repository: https://CRAN.R-project.org/package=rgeos
      repository-code: https://r-forge.r-project.org/projects/rgeos/
      url: https://trac.osgeo.org/geos/
      identifiers:
      - type: url
        value: http://rgeos.r-forge.r-project.org/index.html
      date-released: '2020-09-07'
      contact:
      - family-names: Bivand
        given-names: Roger
        email: Roger.Bivand@nhh.no
        orcid: https://orcid.org/0000-0003-2392-6140
      license: GPL-2.0-or-later
      year: '2020'
    repository: https://CRAN.R-project.org/package=rgeos
    repository-code: https://r-forge.r-project.org/projects/rgeos/
    url: https://trac.osgeo.org/geos/
    date-released: '2020-09-07'
    contact:
    - family-names: Bivand
      given-names: Roger
      email: Roger.Bivand@nhh.no
      orcid: https://orcid.org/0000-0003-2392-6140
    references:
    - type: manual
      title: 'rgeos: Interface to Geometry Engine - Open Source (''GEOS'')'
      authors:
      - family-names: Bivand
        given-names: Roger
        email: Roger.Bivand@nhh.no
        orcid: https://orcid.org/0000-0003-2392-6140
      - family-names: Rundel
        given-names: Colin
      year: '2020'
      notes: R package version 0.5-7
      url: https://CRAN.R-project.org/package=rgeos
    - type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org
    identifiers:
    - type: url
      value: http://rgeos.r-forge.r-project.org/index.html

# Parsed several citations

    - type: manual
      title: 'rgeos: Interface to Geometry Engine - Open Source (''GEOS'')'
      authors:
      - family-names: Bivand
        given-names: Roger
        email: Roger.Bivand@nhh.no
        orcid: https://orcid.org/0000-0003-2392-6140
      - family-names: Rundel
        given-names: Colin
      year: '2020'
      notes: R package version 0.5-7
      url: https://CRAN.R-project.org/package=rgeos
    - type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

# Add wrong field to citation

    cff-version: 1.2.0
    message: If you use this software, please cite it using these metadata.
    title: My Research Software
    authors:
    - family-names: Doe
      given-names: John
    preferred-citation:
      type: manual
      title: My Research Software
      authors:
      - family-names: Doe
        given-names: John
      year: '2022'
    references:
    - type: manual
      title: favoritefood is not valid on cff schema
      authors:
      - family-names: Smith
        given-names: Jane

# Fix wrong orcid

    cff-version: 1.2.0
    message: If you use this software, please cite it using these metadata.
    title: My Research Software
    authors:
    - family-names: Doe
      given-names: John
    preferred-citation:
      type: manual
      title: My Research Software
      authors:
      - family-names: Doe
        given-names: John
      year: '2022'
    references:
    - type: manual
      title: Wrong orcid fixed by cffr
      authors:
      - family-names: Smith
        given-names: Jane
        orcid: https://orcid.org/0000-0000-0000-306X

# Several identifiers and duplicates

    cff-version: 1.2.0
    message: If you use this software, please cite it using these metadata.
    title: My Research Software
    authors:
    - family-names: Doe
      given-names: John
    preferred-citation:
      type: manual
      title: My Research Software
      authors:
      - family-names: Doe
        given-names: John
      year: '2022'
    references:
    - type: manual
      title: A Language and Environment for Statistical Computing
      authors:
      - name: R Core Team
      year: '2022'
      url: https://www.R-project.org/
      doi: 10.5281/zenodo.5366600
      identifiers:
      - type: doi
        value: 10.5281/zenodo.5366601
      - type: doi
        value: 10.5281/zenodo.5366602
      - type: url
        value: https://google.com/

# Test keywords and urls

    cff-version: 1.2.0
    message: If you use this software, please cite it using these metadata.
    title: My Research Software
    authors:
    - family-names: Doe
      given-names: John
    preferred-citation:
      type: manual
      title: My Research Software
      authors:
      - family-names: Doe
        given-names: John
      year: '2022'
    references:
    - type: manual
      title: A Language and Environment for Statistical Computing
      authors:
      - name: R Core Team
      year: '2022'
      url: https://www.R-project.org/
      keywords:
      - Some
      - random keywords
      - in
      - here
      identifiers:
      - type: url
        value: https://google.com/

# Parse persons on CITATION

    type: manual
    title: A Language and Environment for Statistical Computing
    authors:
    - name: R Core Team
    year: '2021'
    contact:
    - family-names: name
      given-names: A
    - family-names: contact
      given-names: A
    conference:
      name: A conference
    database-provider:
      name: Database provider
    editors:
    - family-names: editor
      given-names: A
    - name: Ben and Jerry
    editors-series:
    - family-names: editor series
      given-names: An
    - name: Another
    publisher:
      name: A publisher
      address: A location
    recipients:
    - family-names: recipient
      given-names: A
    senders:
    - name: A Sender
    - family-names: Sender
      given-names: Another
    translators:
    - family-names: one
      given-names: Translator
    - family-names: two
      given-names: Translator

# Test inputs

    type: book
    title: Test
    authors:
    - family-names: Jean
      given-names: Billy
    year: '2021'
    publisher:
      name: Random House

# Merge all DESCRIPTION files with CITATION_basic

    cff-version: 1.2.0
    message: 'To cite package "basicdesc" in publications use:'
    type: software
    title: 'basicdesc: A Basic Description'
    version: 0.1.6
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    abstract: A very basic description. Should parse without problems.
    repository-code: https://github.com/basic/package
    url: https://basic.github.io/package
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "basicdesc" in publications use:'
    type: software
    title: 'basicdesc: A Basic Description'
    version: 0.1.6
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    abstract: A very basic description. Should parse without problems.
    repository-code: https://github.com/basic/package
    url: https://basic.github.io/package
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "basicdescdate" in publications use:'
    type: software
    title: 'basicdescdate: A Basic Description with Date'
    version: 0.1.6
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    abstract: A very basic description. Should parse without problems. I have a Date
    repository-code: https://github.com/basic/package
    url: https://basic.github.io/package
    date-released: '1999-01-01'
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "GenomicRanges" in publications use:'
    type: software
    title: 'GenomicRanges: Representation and manipulation of genomic intervals'
    version: 1.44.0
    authors:
    - name: Bioconductor Package Maintainer
      email: maintainer@bioconductor.org
    - family-names: Aboyoun
      given-names: P.
    - family-names: Pagès
      given-names: H.
    - family-names: Lawrence
      given-names: M.
    abstract: The ability to efficiently represent and manipulate genomic annotations
      and alignments is playing a central role when it comes to analyzing high-throughput
      sequencing data (a.k.a. NGS data). The GenomicRanges package defines general purpose
      containers for storing and manipulating genomic intervals and variables defined
      along a genome. More specialized containers for representing and manipulating short
      alignments against a reference genome, or a matrix-like summarization of an experiment,
      are defined in the GenomicAlignments and SummarizedExperiment packages, respectively.
      Both packages build on top of the GenomicRanges infrastructure.
    repository: https://bioconductor.org/
    repository-code: https://github.com/Bioconductor/GenomicRanges
    url: https://bioconductor.org/packages/GenomicRanges
    date-released: '2021-05-19'
    contact:
    - name: Bioconductor Package Maintainer
      email: maintainer@bioconductor.org
    license: Artistic-2.0
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "codemetar" in publications use:'
    type: software
    title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
    version: 0.1.6
    authors:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    - family-names: Salmon
      given-names: Maëlle
      orcid: https://orcid.org/0000-0002-2815-0399
    abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
      metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
      to generate, parse, and modify 'codemeta.json' files automatically for R packages,
      as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more generally.
    repository: https://CRAN.R-project.org/package=codemetar
    repository-code: https://gitlab.com/ninijay/methoden
    url: https://ropensci.github.io/codemetar
    contact:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    keywords:
    - metadata
    - codemeta
    - ropensci
    - citation
    - credit
    - linked-data
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "manypersons" in publications use:'
    type: software
    title: 'manypersons: A lot of persons'
    version: 0.1.6
    authors:
    - family-names: Hernangómez
      given-names: Diego
      email: fake@gmail.com
      orcid: https://orcid.org/0000-0001-8457-4658
    - family-names: Doe
      given-names: Joe
      affiliation: This One
      country: ES
    - family-names: Doe
      given-names: Pepe
      email: fake@gmail.com
    - name: I am an entity
      date-end: '2020-01-01'
    abstract: Overkill desc with many persons. Try this
    repository-code: https://github.com/many/persons
    url: https://many.github.io/persons
    contact:
    - family-names: Hernangómez
      given-names: Diego
      email: fake@gmail.com
      orcid: https://orcid.org/0000-0001-8457-4658
    - name: I am an entity
      date-end: '2020-01-01'
    keywords:
    - metadata
    - cffr
    - ropensci
    - citation
    - credit
    - linked-data
    - one
    - two
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "manyurls" in publications use:'
    type: software
    title: 'manyurls: A lot of urls'
    version: 0.1.6
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    abstract: This package has many urls. Specifically, 1 Bug Reports and 6 URLs. Expected
      is to have 1 repository-code, 1 url and 3 URLs, since there is 1 duplicate and 1
      invalid url.
    repository-code: https://github.com/test/package
    url: https://test.github.io/package/
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    identifiers:
    - type: url
      value: https://r-forge.r-project.org/projects/test/
    - type: url
      value: http://google.ru
    - type: url
      value: https://gitlab.com/r-packages/behaviorchange
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "codemetar" in publications use:'
    type: software
    title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
    version: 0.1.6
    authors:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    - family-names: Salmon
      given-names: Maëlle
      orcid: https://orcid.org/0000-0002-2815-0399
    abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
      metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
      to generate, parse, and modify 'codemeta.json' files automatically for R packages,
      as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more generally.
    repository: https://CRAN.R-project.org/package=codemetar
    contact:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    keywords:
    - metadata
    - codemeta
    - ropensci
    - citation
    - credit
    - linked-data
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "idonthavencoding" in publications use:'
    type: software
    title: 'idonthavencoding: No Encoding'
    version: 0.1.6
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
      orcid: https://orcid.org/0000-0000-0000-001X
    abstract: I don't have Encoding field. Should parse without problems. Maybe some warning,
      but because the missing Encoding field.
    repository-code: https://github.com/idonthavencoding/package
    url: https://idonthavencoding.github.io/package
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
      orcid: https://orcid.org/0000-0000-0000-001X
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "codemetar" in publications use:'
    type: software
    title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
    version: 0.3.2
    authors:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    - family-names: Salmon
      given-names: Maëlle
      orcid: https://orcid.org/0000-0002-2815-0399
    abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
      metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
      to generate, parse, and modify 'codemeta.json' files automatically for R packages,
      as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more generally.
    repository: https://ropensci.r-universe.dev
    repository-code: https://github.com/ropensci/codemetar
    url: https://docs.ropensci.org/codemetar/
    contact:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    keywords:
    - metadata
    - codemeta
    - ropensci
    - citation
    - credit
    - linked-data
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "rgeos" in publications use:'
    type: software
    title: 'rgeos: Interface to Geometry Engine - Open Source (''GEOS'')'
    version: 0.5-7
    authors:
    - family-names: Bivand
      given-names: Roger
      email: Roger.Bivand@nhh.no
      orcid: https://orcid.org/0000-0003-2392-6140
    - family-names: Rundel
      given-names: Colin
    abstract: 'Interface to Geometry Engine - Open Source (''GEOS'') using the C ''API''
      for topology operations on geometries. Please note that ''rgeos'' will be retired
      by the end of 2023, plan transition to sf functions using ''GEOS'' at your earliest
      convenience. The ''GEOS'' library is external to the package, and, when installing
      the package from source, must be correctly installed first. Windows and Mac Intel
      OS X binaries are provided on ''CRAN''. (''rgeos'' >= 0.5-1): Up to and including
      ''GEOS'' 3.7.1, topological operations succeeded with some invalid geometries for
      which the same operations fail from and including ''GEOS'' 3.7.2. The ''checkValidity=''
      argument defaults and structure have been changed, from default FALSE to integer
      default ''0L'' for ''GEOS'' < 3.7.2 (no check), ''1L'' ''GEOS'' >= 3.7.2 (check
      and warn). A value of ''2L'' is also provided that may be used, assigned globally
      using ''set_RGEOS_CheckValidity(2L)'', or locally using the ''checkValidity=2L''
      argument, to attempt zero-width buffer repair if invalid geometries are found. The
      previous default (FALSE, now ''0L'') is fastest and used for ''GEOS'' < 3.7.2, but
      will not warn users of possible problems before the failure of topological operations
      that previously succeeded. From ''GEOS'' 3.8.0, repair of geometries may also be
      attempted using ''gMakeValid()'', which may, however, return a collection of geometries
      of different types.'
    repository: https://CRAN.R-project.org/package=rgeos
    repository-code: https://r-forge.r-project.org/projects/rgeos/
    url: https://trac.osgeo.org/geos/
    date-released: '2020-09-07'
    contact:
    - family-names: Bivand
      given-names: Roger
      email: Roger.Bivand@nhh.no
      orcid: https://orcid.org/0000-0003-2392-6140
    license: GPL-2.0-or-later
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    identifiers:
    - type: url
      value: http://rgeos.r-forge.r-project.org/index.html
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "surveillance" in publications use:'
    type: software
    title: 'surveillance: Temporal and Spatio-Temporal Modeling and Monitoring of Epidemic
      Phenomena'
    version: 1.19.1
    authors:
    - family-names: Höhle
      given-names: Michael
      email: hoehle@math.su.se
      orcid: https://orcid.org/0000-0002-0423-6702
    - family-names: Meyer
      given-names: Sebastian
      email: seb.meyer@fau.de
      orcid: https://orcid.org/0000-0002-1791-9449
    - family-names: Paul
      given-names: Michaela
    abstract: Statistical methods for the modeling and monitoring of time series of counts,
      proportions and categorical data, as well as for the modeling of continuous-time
      point processes of epidemic phenomena. The monitoring methods focus on aberration
      detection in count data time series from public health surveillance of communicable
      diseases, but applications could just as well originate from environmetrics, reliability
      engineering, econometrics, or social sciences. The package implements many typical
      outbreak detection procedures such as the (improved) Farrington algorithm, or the
      negative binomial GLR-CUSUM method of Höhle and Paul (2008) <doi:10.1016/j.csda.2008.02.015>.
      A novel CUSUM approach combining logistic and multinomial logistic modeling is also
      included. The package contains several real-world data sets, the ability to simulate
      outbreak data, and to visualize the results of the monitoring in a temporal, spatial
      or spatio-temporal fashion. A recent overview of the available monitoring procedures
      is given by Salmon et al. (2016) <doi:10.18637/jss.v070.i10>. For the retrospective
      analysis of epidemic spread, the package provides three endemic-epidemic modeling
      frameworks with tools for visualization, likelihood inference, and simulation. hhh4()
      estimates models for (multivariate) count time series following Paul and Held (2011)
      <doi:10.1002/sim.4177> and Meyer and Held (2014) <doi:10.1214/14-AOAS743>. twinSIR()
      models the susceptible-infectious-recovered (SIR) event history of a fixed population,
      e.g, epidemics across farms or networks, as a multivariate point process as proposed
      by Höhle (2009) <doi:10.1002/bimj.200900050>. twinstim() estimates self-exciting
      point process models for a spatio-temporal point pattern of infective events, e.g.,
      time-stamped geo-referenced surveillance data, as proposed by Meyer et al. (2012)
      <doi:10.1111/j.1541-0420.2011.01684.x>. A recent overview of the implemented space-time
      modeling frameworks for epidemic phenomena is given by Meyer et al. (2017) <doi:10.18637/jss.v077.i11>.
    repository: https://CRAN.R-project.org/package=surveillance
    url: https://surveillance.R-Forge.R-project.org/
    date-released: '2021-03-30'
    contact:
    - family-names: Meyer
      given-names: Sebastian
      email: seb.meyer@fau.de
      orcid: https://orcid.org/0000-0002-1791-9449
    license: GPL-2.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "codemetar" in publications use:'
    type: software
    title: 'codemetar: Generate ''CodeMeta'' Metadata for R Packages'
    version: 0.1.6
    authors:
    - family-names: Ok
      given-names: John
      email: email@email.edu
    - family-names: Doe
      given-names: Jane
      email: email2@email.edu
    - family-names: Doo
      given-names: Jane
    abstract: The 'Codemeta' Project defines a 'JSON-LD' format for describing software
      metadata, as detailed at <https://codemeta.github.io>. This package provides utilities
      to generate, parse, and modify 'codemeta.json' files automatically for R packages,
      as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more generally.
    repository: https://CRAN.R-project.org/package=codemetar
    repository-code: https://github.com/ropensci/codemetar
    url: https://ropensci.github.io/codemetar
    contact:
    - family-names: Ok
      given-names: John
      email: email@email.edu
    - family-names: Doe
      given-names: Jane
      email: email2@email.edu
    license: GPL-3.0-only
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

---

    cff-version: 1.2.0
    message: 'To cite package "wrongurls" in publications use:'
    type: software
    title: 'wrongurls: Generate CodeMeta Metadata for R Packages'
    version: 0.1.0
    authors:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    abstract: Codemeta defines a 'JSON-LD' format for describing software metadata. This
      package provides utilities to generate, parse, and modify codemeta.jsonld files
      automatically for R packages.
    url: https://httpbin.org/status/404
    contact:
    - family-names: Boettiger
      given-names: Carl
      email: cboettig@gmail.com
      orcid: https://orcid.org/0000-0002-1642-628X
    keywords:
    - metadata
    - codemeta
    - ropensci
    - citation
    - credit
    - linked-data
    license: MIT
    doi: 10.1111/2041-210X.12469
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    identifiers:
    - type: url
      value: https://httpbin.org/status/429
    - type: url
      value: https://www.github.es/ropensci/codemeta
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org

# preferred-citation-book

    @Book{bueler:2021,
      title = {PETSc for Partial Differential Equations: Numerical Solutions in C and Python},
      author = {Ed Bueler},
      year = {2021},
      publisher = {SIAM Press},
      address = {Philadelphia},
      isbn = {978111976304},
      url = {https://github.com/bueler/p4pdes},
    }

# preferred-citation-conference-paper-2

    @InProceedings{gamblinlegendre:2016,
      title = {The Spack Package Manager: Bringing Order to HPC Software Chaos},
      author = {Todd Gamblin and Matthew LeGendre and Michael R. Collette and Gregory L. Lee and Adam Moody and Bronis R. {de Supinski} and Scott Futral},
      year = {2016},
      month = {dec},
      booktitle = {Proceedings of Supercomputing 2015},
      publisher = {ACM/IEEE},
      address = {Austin, Texas, USA},
      series = {Supercomputing 2015 (SC'15)},
      doi = {10.1244/2886907.2879400},
    }

# preferred-citation-conference-paper-missing

    @InProceedings{gamblinlegendre:2015,
      title = {The Spack Package Manager: Bringing Order to HPC Software Chaos},
      author = {Todd Gamblin and Matthew LeGendre and Michael R. Collette and Gregory L. Lee and Adam Moody and Bronis R. {de Supinski} and Scott Futral},
      year = {2015},
      month = {nov},
      booktitle = {Proceedings of Supercomputing 2015},
      doi = {10.1244/2886907.2879400},
    }

# preferred-citation-conference-paper

    @InProceedings{rampinfreire:2016,
      title = {ReproZip: Computational Reproducibility With Ease},
      author = {Rémi Rampin and Juliana Freire and Fernando Chirigati and Dennis Shasha},
      year = {2016},
      month = {jun},
      booktitle = {Proceedings of the 2016 ACM SIGMOD International Conference on Management of Data (SIGMOD)},
      publisher = {ACM},
      address = {San Francisco, US},
      series = {SIGMOD '16},
      pages = {2085--2088},
      doi = {10.1145/2882903.2899401},
      date = {2016-06-26},
    }

# preferred-citation-manual

    @Manual{hernang?mez:2021,
      title = {cffr: Generate Citation File Format ('cff') Metadata for R Packages},
      author = {Diego Hernangómez},
      year = {2021},
      doi = {10.5281/zenodo.5509766},
      url = {https://dieghernan.github.io/cffr/},
    }

# preferred-citation-no-month

    @Article{lisa:2021,
      title = {My awesome research software},
      author = {Mona Lisa},
      year = {2021},
      journal = {Journal Title},
    }

# preferred-citation-no-vol

    @Article{hartmannwong,
      title = {An image-based data-driven analysis of cellular architecture in a developing tissue},
      author = {Jonas Hartmann and Mie Wong and Elisa Gallo and Darren Gilmour},
      year = {2020},
      journal = {eLife},
      pages = {888},
      doi = {10.7554/eLife.55913},
      url = {https://elifesciences.org/articles/55913},
    }

# preferred-citation-pamphlet

    @Booklet{haines:2021,
      title = {Why and how to use CFF},
      author = {Robert Haines},
      year = {2021},
      month = {sep},
      doi = {10.5281/zenodo.1184077},
      url = {https://github.com/citation-file-format/ruby-cff},
    }

# preferred-citation-report-no-institution

    @TechReport{haines:2021,
      title = {The benefits of using CFF files},
      author = {Robert Haines},
      year = {2021},
      month = {sep},
      doi = {10.5281/zenodo.1184077},
      url = {https://github.com/citation-file-format/ruby-cff},
      institution = {The University of Manchester},
    }

# preferred-citation-report

    @TechReport{haines:2021,
      title = {The benefits of using CFF files},
      author = {Robert Haines},
      year = {2021},
      month = {sep},
      doi = {10.5281/zenodo.1184077},
      url = {https://github.com/citation-file-format/ruby-cff},
      institution = {The CFF Institute},
    }

# preferred-citation-unpublished

    @Unpublished{haines:2021,
      title = {Why and how to use CFF},
      author = {Robert Haines},
      year = {2021},
      month = {sep},
      doi = {10.5281/zenodo.1184077},
      url = {https://github.com/citation-file-format/ruby-cff},
      note = {Self-published by the author.},
    }

# reprozip

    @Proceedings{rampinfreire:2016,
      title = {ReproZip: Computational Reproducibility With Ease},
      author = {Remi Rampin and Juliana Freire and Fernando Chirigati and Dennis Shasha},
      year = {2016},
      month = {jun},
      address = {San Francisco, US},
      series = {SIGMOD '16},
      doi = {10.1145/2882903.2899401},
      abstract = {We present ReproZip, the recommended packaging tool for the SIGMOD Reproducibility Review. ReproZip was designed to simplify the process of making an existing computational experiment reproducible across platforms, even when the experiment was put together without reproducibility in mind. The tool creates a self-contained package for an experiment by automatically tracking and identifying all its required dependencies. The researcher can share the package with others, who can then use ReproZip to unpack the experiment, reproduce the findings on their favorite operating system, as well as modify the original experiment for reuse in new research, all with little effort. The demo will consist of examples of non-trivial experiments, showing how these can be packed in a Linux machine and reproduced on different machines and operating systems. Demo visitors will also be able to pack and reproduce their own experiments.},
      date = {2016-06-26},
    }

# smith-et-al

    @Article{smithkatz:2016,
      title = {Software citation principles},
      author = {A. M. Smith and D. S. Katz and K. E. Niemeyer and {FORCE11 Software Citation Working Group}},
      year = {2016},
      month = {sep},
      journal = {PeerJ Computer Science},
      volume = {2},
      number = {123},
      pages = {e86},
      doi = {10.7717/peerj-cs.86},
    }

# tidyverse-joss-paper

    @Article{wickham:2019,
      title = {Welcome to the Tidyverse},
      author = {Hadley Wickham},
      year = {2019},
      month = {nov},
      journal = {Journal of Open Source Software},
      volume = {4},
      number = {43},
      pages = {1686},
      doi = {10.21105/joss.01686},
    }

# Creating cff from packages encoded in latin1

    cff-version: 1.2.0
    message: 'To cite package "surveillance" in publications use:'
    type: software
    license: GPL-2.0-only
    title: 'surveillance: Temporal and Spatio-Temporal Modeling and Monitoring of Epidemic
      Phenomena'
    version: 1.19.1
    abstract: Statistical methods for the modeling and monitoring of time series of counts,
      proportions and categorical data, as well as for the modeling of continuous-time
      point processes of epidemic phenomena. The monitoring methods focus on aberration
      detection in count data time series from public health surveillance of communicable
      diseases, but applications could just as well originate from environmetrics, reliability
      engineering, econometrics, or social sciences. The package implements many typical
      outbreak detection procedures such as the (improved) Farrington algorithm, or the
      negative binomial GLR-CUSUM method of Höhle and Paul (2008) <doi:10.1016/j.csda.2008.02.015>.
      A novel CUSUM approach combining logistic and multinomial logistic modeling is also
      included. The package contains several real-world data sets, the ability to simulate
      outbreak data, and to visualize the results of the monitoring in a temporal, spatial
      or spatio-temporal fashion. A recent overview of the available monitoring procedures
      is given by Salmon et al. (2016) <doi:10.18637/jss.v070.i10>. For the retrospective
      analysis of epidemic spread, the package provides three endemic-epidemic modeling
      frameworks with tools for visualization, likelihood inference, and simulation. hhh4()
      estimates models for (multivariate) count time series following Paul and Held (2011)
      <doi:10.1002/sim.4177> and Meyer and Held (2014) <doi:10.1214/14-AOAS743>. twinSIR()
      models the susceptible-infectious-recovered (SIR) event history of a fixed population,
      e.g, epidemics across farms or networks, as a multivariate point process as proposed
      by Höhle (2009) <doi:10.1002/bimj.200900050>. twinstim() estimates self-exciting
      point process models for a spatio-temporal point pattern of infective events, e.g.,
      time-stamped geo-referenced surveillance data, as proposed by Meyer et al. (2012)
      <doi:10.1111/j.1541-0420.2011.01684.x>. A recent overview of the implemented space-time
      modeling frameworks for epidemic phenomena is given by Meyer et al. (2017) <doi:10.18637/jss.v077.i11>.
    authors:
    - family-names: Höhle
      given-names: Michael
      email: hoehle@math.su.se
      orcid: https://orcid.org/0000-0002-0423-6702
    - family-names: Meyer
      given-names: Sebastian
      email: seb.meyer@fau.de
      orcid: https://orcid.org/0000-0002-1791-9449
    - family-names: Paul
      given-names: Michaela
    preferred-citation:
      type: manual
      title: 'surveillance: Temporal and Spatio-Temporal Modeling and Monitoring of Epidemic
        Phenomena'
      authors:
      - family-names: Höhle
        given-names: Michael
        email: hoehle@math.su.se
        orcid: https://orcid.org/0000-0002-0423-6702
      - family-names: Meyer
        given-names: Sebastian
        email: seb.meyer@fau.de
        orcid: https://orcid.org/0000-0002-1791-9449
      - family-names: Paul
        given-names: Michaela
      version: 1.19.1
      abstract: Statistical methods for the modeling and monitoring of time series of
        counts, proportions and categorical data, as well as for the modeling of continuous-time
        point processes of epidemic phenomena. The monitoring methods focus on aberration
        detection in count data time series from public health surveillance of communicable
        diseases, but applications could just as well originate from environmetrics, reliability
        engineering, econometrics, or social sciences. The package implements many typical
        outbreak detection procedures such as the (improved) Farrington algorithm, or
        the negative binomial GLR-CUSUM method of Höhle and Paul (2008) <doi:10.1016/j.csda.2008.02.015>.
        A novel CUSUM approach combining logistic and multinomial logistic modeling is
        also included. The package contains several real-world data sets, the ability
        to simulate outbreak data, and to visualize the results of the monitoring in a
        temporal, spatial or spatio-temporal fashion. A recent overview of the available
        monitoring procedures is given by Salmon et al. (2016) <doi:10.18637/jss.v070.i10>.
        For the retrospective analysis of epidemic spread, the package provides three
        endemic-epidemic modeling frameworks with tools for visualization, likelihood
        inference, and simulation. hhh4() estimates models for (multivariate) count time
        series following Paul and Held (2011) <doi:10.1002/sim.4177> and Meyer and Held
        (2014) <doi:10.1214/14-AOAS743>. twinSIR() models the susceptible-infectious-recovered
        (SIR) event history of a fixed population, e.g, epidemics across farms or networks,
        as a multivariate point process as proposed by Höhle (2009) <doi:10.1002/bimj.200900050>.
        twinstim() estimates self-exciting point process models for a spatio-temporal
        point pattern of infective events, e.g., time-stamped geo-referenced surveillance
        data, as proposed by Meyer et al. (2012) <doi:10.1111/j.1541-0420.2011.01684.x>.
        A recent overview of the implemented space-time modeling frameworks for epidemic
        phenomena is given by Meyer et al. (2017) <doi:10.18637/jss.v077.i11>.
      repository: https://CRAN.R-project.org/package=surveillance
      url: https://surveillance.R-Forge.R-project.org/
      date-released: '2021-03-30'
      contact:
      - family-names: Meyer
        given-names: Sebastian
        email: seb.meyer@fau.de
        orcid: https://orcid.org/0000-0002-1791-9449
      license: GPL-2.0-only
      year: '2021'
    repository: https://CRAN.R-project.org/package=surveillance
    url: https://surveillance.R-Forge.R-project.org/
    date-released: '2021-03-30'
    contact:
    - family-names: Meyer
      given-names: Sebastian
      email: seb.meyer@fau.de
      orcid: https://orcid.org/0000-0002-1791-9449
    references:
    - type: article
      title: 'Monitoring Count Time Series in R: Aberration Detection in Public Health
        Surveillance'
      authors:
      - family-names: Salmon
        given-names: Maëlle
      - family-names: Schumacher
        given-names: Dirk
      - family-names: Höhle
        given-names: Michael
      journal: Journal of Statistical Software
      year: '2016'
      volume: '70'
      issue: '10'
      doi: 10.18637/jss.v070.i10
      start: '1'
      end: '35'
    - type: article
      title: Spatio-Temporal Analysis of Epidemic Phenomena Using the R Package surveillance
      authors:
      - family-names: Meyer
        given-names: Sebastian
      - family-names: Held
        given-names: Leonhard
      - family-names: Höhle
        given-names: Michael
      journal: Journal of Statistical Software
      year: '2017'
      volume: '77'
      issue: '11'
      doi: 10.18637/jss.v077.i11
      start: '1'
      end: '55'

# Auto generate preferred citations

    cff-version: 1.2.0
    message: 'To cite package "rgeos" in publications use:'
    type: software
    license: GPL-2.0-or-later
    title: 'rgeos: Interface to Geometry Engine - Open Source (''GEOS'')'
    version: 0.5-7
    abstract: 'Interface to Geometry Engine - Open Source (''GEOS'') using the C ''API''
      for topology operations on geometries. Please note that ''rgeos'' will be retired
      by the end of 2023, plan transition to sf functions using ''GEOS'' at your earliest
      convenience. The ''GEOS'' library is external to the package, and, when installing
      the package from source, must be correctly installed first. Windows and Mac Intel
      OS X binaries are provided on ''CRAN''. (''rgeos'' >= 0.5-1): Up to and including
      ''GEOS'' 3.7.1, topological operations succeeded with some invalid geometries for
      which the same operations fail from and including ''GEOS'' 3.7.2. The ''checkValidity=''
      argument defaults and structure have been changed, from default FALSE to integer
      default ''0L'' for ''GEOS'' < 3.7.2 (no check), ''1L'' ''GEOS'' >= 3.7.2 (check
      and warn). A value of ''2L'' is also provided that may be used, assigned globally
      using ''set_RGEOS_CheckValidity(2L)'', or locally using the ''checkValidity=2L''
      argument, to attempt zero-width buffer repair if invalid geometries are found. The
      previous default (FALSE, now ''0L'') is fastest and used for ''GEOS'' < 3.7.2, but
      will not warn users of possible problems before the failure of topological operations
      that previously succeeded. From ''GEOS'' 3.8.0, repair of geometries may also be
      attempted using ''gMakeValid()'', which may, however, return a collection of geometries
      of different types.'
    authors:
    - family-names: Bivand
      given-names: Roger
      email: Roger.Bivand@nhh.no
      orcid: https://orcid.org/0000-0003-2392-6140
    - family-names: Rundel
      given-names: Colin
    preferred-citation:
      type: manual
      title: 'rgeos: Interface to Geometry Engine - Open Source (''GEOS'')'
      authors:
      - family-names: Bivand
        given-names: Roger
        email: Roger.Bivand@nhh.no
        orcid: https://orcid.org/0000-0003-2392-6140
      - family-names: Rundel
        given-names: Colin
      version: 0.5-7
      abstract: 'Interface to Geometry Engine - Open Source (''GEOS'') using the C ''API''
        for topology operations on geometries. Please note that ''rgeos'' will be retired
        by the end of 2023, plan transition to sf functions using ''GEOS'' at your earliest
        convenience. The ''GEOS'' library is external to the package, and, when installing
        the package from source, must be correctly installed first. Windows and Mac Intel
        OS X binaries are provided on ''CRAN''. (''rgeos'' >= 0.5-1): Up to and including
        ''GEOS'' 3.7.1, topological operations succeeded with some invalid geometries
        for which the same operations fail from and including ''GEOS'' 3.7.2. The ''checkValidity=''
        argument defaults and structure have been changed, from default FALSE to integer
        default ''0L'' for ''GEOS'' < 3.7.2 (no check), ''1L'' ''GEOS'' >= 3.7.2 (check
        and warn). A value of ''2L'' is also provided that may be used, assigned globally
        using ''set_RGEOS_CheckValidity(2L)'', or locally using the ''checkValidity=2L''
        argument, to attempt zero-width buffer repair if invalid geometries are found.
        The previous default (FALSE, now ''0L'') is fastest and used for ''GEOS'' < 3.7.2,
        but will not warn users of possible problems before the failure of topological
        operations that previously succeeded. From ''GEOS'' 3.8.0, repair of geometries
        may also be attempted using ''gMakeValid()'', which may, however, return a collection
        of geometries of different types.'
      repository: https://CRAN.R-project.org/package=rgeos
      repository-code: https://r-forge.r-project.org/projects/rgeos/
      url: https://trac.osgeo.org/geos/
      identifiers:
      - type: url
        value: http://rgeos.r-forge.r-project.org/index.html
      date-released: '2020-09-07'
      contact:
      - family-names: Bivand
        given-names: Roger
        email: Roger.Bivand@nhh.no
        orcid: https://orcid.org/0000-0003-2392-6140
      license: GPL-2.0-or-later
      year: '2020'
    repository: https://CRAN.R-project.org/package=rgeos
    repository-code: https://r-forge.r-project.org/projects/rgeos/
    url: https://trac.osgeo.org/geos/
    date-released: '2020-09-07'
    contact:
    - family-names: Bivand
      given-names: Roger
      email: Roger.Bivand@nhh.no
      orcid: https://orcid.org/0000-0003-2392-6140
    identifiers:
    - type: url
      value: http://rgeos.r-forge.r-project.org/index.html

---

    cff-version: 1.2.0
    message: 'To cite package "basicdescdate" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'basicdescdate: A Basic Description with Date'
    version: 0.1.6
    abstract: A very basic description. Should parse without problems. I have a Date
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    preferred-citation:
      type: manual
      title: 'basicdescdate: A Basic Description with Date'
      authors:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      version: 0.1.6
      abstract: A very basic description. Should parse without problems. I have a Date
      repository-code: https://github.com/basic/package
      url: https://basic.github.io/package
      date-released: '1999-01-01'
      contact:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      license: GPL-3.0-only
      year: '1999'
    repository-code: https://github.com/basic/package
    url: https://basic.github.io/package
    date-released: '1999-01-01'
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com

# Fuzzy match on cff_create

    
    
    ## Fuzzy match on cff_create 
    
    cff-version: 1.2.0
    message: If you use this software, please cite it using these metadata.
    title: Modifying a 'cff' object
    version: 0.0.1
    authors:
    - family-names: Doe
      given-names: John
    preferred-citation:
      type: manual
      title: My Research Software
      authors:
      - family-names: Doe
        given-names: John
      year: '2022'
    repository: https://github.com/user/repo
    url: https://ropensci.org/
    
    ---

# Article

    type: article
    title: Literate Programming
    authors:
    - name: R Core Team
    journal: The Computer Journal
    year: '1984'
    volume: '27'
    issue: '2'
    month: '1'
    notes: Example modified for testing purposes
    start: '97'
    end: '111'

# Book

    type: book
    title: The LaTeX Companion
    authors:
    - family-names: Mittelbach
      given-names: Frank
    - family-names: Gossens
      given-names: Michel
    - family-names: Braams
      given-names: Johannes
    - family-names: Carlisle
      given-names: David
    - family-names: Rowley
      given-names: Chris
    editors:
    - name: Barnes and Noble
    publisher:
      name: Addison-Wesley Professional
      address: Santa Monica
    year: '2004'
    volume: '3'
    issue: '7'
    edition: Fourth
    month: '8'
    notes: Example modified for testing purposes
    keywords:
    - Two
    - keyword
    collection-title: The LateX Books

# Booklet

    type: pamphlet
    title: Java Booklet
    authors:
    - family-names: Mustermann
      given-names: Max
    medium: Internet
    location:
      name: Stuttgart
    month: '2'
    year: '2016'
    notes: Example modified from Jabref

# Conference

    type: conference-paper
    title: On Notions of Information Transfer in VLSI Circuits
    authors:
    - family-names: Oaho
      given-names: Alfred V.
    - family-names: Ullman
      given-names: Jeffrey D.
    - family-names: Yannakakis
      given-names: Mihalis
    collection-title: Proc. Fifteenth Annual ACM STOC
    year: '1983'
    editors:
    - family-names: Oz
      given-names: Wizard V.
    - family-names: Yannakakis
      given-names: Mihalis
    volume: '41'
    issue: '17'
    location:
      name: Boston
    publisher:
      name: Academic Press
    notes: Example modified for testing purposes
    start: '133'
    end: '139'
    conference:
      name: All ACM Conferences
      address: Boston
    institution:
      name: ACM

# InBook

    type: book
    title: A Framework for Freeness Analysis
    authors:
    - family-names: King
      given-names: A.
    editors:
    - family-names: Tick
      given-names: E
    - family-names: Succi
      given-names: G
    section: 7, 14
    publisher:
      name: Kluwer Academic Publishers
      address: Dordrecht
    year: '1994'
    volume: '27'
    issue: '2'
    edition: Second
    month: '1'
    notes: Example modified for testing purposes
    start: '137'
    end: '149'
    collection-title: Implementations of Logic Programming Systems

# InCollection

    type: generic
    title: Knowledge-Based Methods for WSD
    authors:
    - family-names: Mihalcea
      given-names: Rada
    collection-title: 'Word Sense Disambiguation: Algorithms and Applications'
    publisher:
      name: Springer
      address: 107--132
    year: '2006'
    editors:
    - family-names: Agirre
      given-names: Eneko
    - family-names: Edmonds
      given-names: Philip
    volume: '23'
    issue: '3'
    section: '1,2,3'
    edition: Third
    month: '8'
    notes: A note
    start: '24'
    end: '57'

# InProceedings

    type: conference-paper
    title: On Notions of Information Transfer in VLSI Circuits
    authors:
    - family-names: Oaho
      given-names: Alfred V.
    - family-names: Ullman
      given-names: Jeffrey D.
    - family-names: Yannakakis
      given-names: Mihalis
    collection-title: Proc. Fifteenth Annual ACM STOC
    year: '1983'
    editors:
    - family-names: Oz
      given-names: Wizard V.
    - family-names: Yannakakis
      given-names: Mihalis
    volume: '41'
    issue: '17'
    location:
      name: Boston
    publisher:
      name: Academic Press
    notes: Example modified for testing purposes
    start: '133'
    end: '139'
    conference:
      name: All ACM Conferences
      address: Boston
    institution:
      name: ACM

# Manual

    type: manual
    title: A Language and Environment for Statistical Computing
    authors:
    - name: R Core Team
    location:
      name: Vienna, Austria
    edition: Fourth
    month: '8'
    year: '2021'
    notes: Example modified for testing purposes
    institution:
      name: R Foundation for Statistical Computing

# MastersThesis

    type: thesis
    title: An examination of keystroke dynamics for continuous user authentication
    authors:
    - family-names: Alsolami
      given-names: Eesa
    year: '2012'
    month: '8'
    notes: Example modified for testing purposes
    institution:
      name: Queensland University of Technology
      address: Queensland, NZ
    thesis-type: Master's Thesis

# Misc

    type: generic
    title: A Language and Environment for Statistical Computing
    authors:
    - name: R Core Team
    medium: CD-ROM
    month: '1'
    year: '2021'
    notes: A note

# PhdThesis

    type: thesis
    title: An examination of keystroke dynamics for continuous user authentication
    authors:
    - family-names: Alsolami
      given-names: Eesa
    year: '2012'
    month: '8'
    notes: Example modified for testing purposes
    institution:
      name: Queensland University of Technology
      address: Queensland, NZ
    thesis-type: PhD Thesis

# Proceedings

    type: proceedings
    title: Proc. Fifteenth Annual STOC
    authors:
    - name: anonymous
    year: '1983'
    editors:
    - family-names: Oz
      given-names: Wizard V.
    - family-names: Yannakakis
      given-names: Mihalis
    volume: '1'
    issue: '17'
    location:
      name: Boston, US
    month: '8'
    publisher:
      name: Academic Press
    notes: Example modified for testing purposes
    conference:
      name: All ACM Conferences
      address: Boston, US
    institution:
      name: The OX Association for Computing Machinery

# TechReport

    type: report
    title: Naive tools for studying compilation histories
    authors:
    - family-names: Jadud
      given-names: Matthew C.
    - family-names: Fincher
      given-names: Sally A.
    institution:
      name: University of Kent Canterbury
      address: Computing Laboratory, University of Kent, Canterbury, Kent, CT2 7NF
    year: '2003'
    issue: 3-03
    month: '3'
    notes: Example modified for testing purposes

# Unpublished

    type: unpublished
    title: Demonstratives
    authors:
    - family-names: Kaplan
      given-names: D.
    notes: Unpublished manuscript, UCLA
    year: '1977'
    month: '8'

# Test entry without author

    type: proceedings
    title: Proceedings of the 6th European Conference on Computer Systems
    authors:
    - name: anonymous
    editors:
    - family-names: Berbers
      given-names: Yolande
    - family-names: Zwaenepoel
      given-names: Willy
    collection-title: Proceedings of the 6th European Conference on Computer Systems
    publisher:
      name: ACM
    month: '4'
    year: '2006'
    isbn: 1-59593-322-02

# Test entry without author but has a key

    type: generic
    title: Proceedings of the 6th European Conference on Computer Systems
    authors:
    - name: anonymous
    collection-title: Proceedings of the 6th European Conference on Computer Systems
    publisher:
      name: ACM
    month: '4'
    year: '2006'
    isbn: 1-59593-322-02

# Test entry without author and key

    type: generic
    title: Proceedings of the 6th European Conference on Computer Systems
    authors:
    - name: anonymous
    collection-title: Proceedings of the 6th European Conference on Computer Systems
    publisher:
      name: ACM
    month: '4'
    year: '2006'
    isbn: 1-59593-322-02

# Skip misc without title

    cff-version: 1.2.0
    message: If you use this software, please cite it using these metadata.
    title: My Research Software
    authors:
    - family-names: Doe
      given-names: John
    preferred-citation:
      type: manual
      title: My Research Software
      authors:
      - family-names: Doe
        given-names: John
      year: '2022'

# Skip misc without title, not skipping the good one

    cff-version: 1.2.0
    message: If you use this software, please cite it using these metadata.
    title: My Research Software
    authors:
    - family-names: Doe
      given-names: John
    preferred-citation:
      type: manual
      title: My Research Software
      authors:
      - family-names: Doe
        given-names: John
      year: '2022'
    references:
    - type: generic
      title: 'rromeo: An R Client for SHERPA/RoMEO API'
      authors:
      - family-names: Grenié
        given-names: Matthias
      - family-names: Gruson
        given-names: Hugo
      year: '2019'
      url: https://CRAN.R-project.org/package=rromeo

# Check extended BibLatex Fields

    type: article
    title: Computation of methodology hyphen independent ionic solvation free energies
      from molecular simulations
    authors:
    - family-names: Kastenholz
      given-names: M. A.
    - family-names: Hünenbergerb
      given-names: Philippe H.
    journal: J. Chem. Phys.
    year: '2006'
    notes: Example modified for testing purposes
    date-published: '2006-03-15'
    filename: a_file.pdf
    issue-title: Semantic 3D Media and Content
    translators:
    - family-names: Wicksteed
      given-names: P. H.
    - family-names: Cornford
      given-names: F. M.
    date-accessed: '2006-10-01'
    pages: '528'
    abstract: The computation of ionic solvation free energies from atomistic simulations
      is a surprisingly difficult problem that has found no satisfactory solution for
      more than 15 years.
    doi: 10.1063/1.2172593
    isbn: 0-816-52066-6
    issn: 0097-8493
    url: http://www.ctan.org
    start: '55'
    end: '65'
    month: '3'

# Parse CITATION_basic

    
    To cite RNeXML in publications, please use:
    
    Boettiger C, Chamberlain S, Vos R, Lapp H (2016). "RNeXML: A Package
    for Reading and Writing Richly Annotated Phylogenetic, Character, and
    Trait Data in R." _Methods in Ecology and Evolution_, *7*, 352-357.
    doi: 10.1111/2041-210X.12469 (URL:
    https://doi.org/10.1111/2041-210X.12469).
    
    A BibTeX entry for LaTeX users is
    
      @Article{,
        title = {{RNeXML}: {A} Package for Reading and Writing Richly Annotated Phylogenetic, Character, and Trait Data in {R}},
        journal = {Methods in Ecology and Evolution},
        author = {Carl Boettiger and Scott Chamberlain and Rutger Vos and Hilmar Lapp},
        year = {2016},
        volume = {7},
        pages = {352--357},
        doi = {10.1111/2041-210X.12469},
      }
    
      H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
      Springer-Verlag New York, 2016.
    
    A BibTeX entry for LaTeX users is
    
      @Book{,
        author = {Hadley Wickham},
        title = {ggplot2: Elegant Graphics for Data Analysis},
        publisher = {Springer-Verlag New York},
        year = {2016},
        isbn = {978-3-319-24277-4},
        url = {https://ggplot2.tidyverse.org},
      }
    

# Parse CITATION with no encoding

    
    To cite RNeXML in publications, please use:
    
    Boettiger C, Chamberlain S, Vos R, Lapp H (2016). "RNeXML: A Package
    for Reading and Writing Richly Annotated Phylogenetic, Character, and
    Trait Data in R." _Methods in Ecology and Evolution_, *7*, 352-357.
    doi: 10.1111/2041-210X.12469 (URL:
    https://doi.org/10.1111/2041-210X.12469).
    
    A BibTeX entry for LaTeX users is
    
      @Article{,
        title = {{RNeXML}: {A} Package for Reading and Writing Richly Annotated Phylogenetic, Character, and Trait Data in {R}},
        journal = {Methods in Ecology and Evolution},
        author = {Carl Boettiger and Scott Chamberlain and Rutger Vos and Hilmar Lapp},
        year = {2016},
        volume = {7},
        pages = {352--357},
        doi = {10.1111/2041-210X.12469},
      }
    
      H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
      Springer-Verlag New York, 2016.
    
    A BibTeX entry for LaTeX users is
    
      @Book{,
        author = {Hadley Wickham},
        title = {ggplot2: Elegant Graphics for Data Analysis},
        publisher = {Springer-Verlag New York},
        year = {2016},
        isbn = {978-3-319-24277-4},
        url = {https://ggplot2.tidyverse.org},
      }
    

# Parse CITATION_auto

    
      Roger Bivand and Colin Rundel (2020). rgeos: Interface to Geometry
      Engine - Open Source ('GEOS'). R package version 0.5-7.
      https://CRAN.R-project.org/package=rgeos
    
    A BibTeX entry for LaTeX users is
    
      @Manual{,
        title = {rgeos: Interface to Geometry Engine - Open Source ('GEOS')},
        author = {Roger Bivand and Colin Rundel},
        year = {2020},
        note = {R package version 0.5-7},
        url = {https://CRAN.R-project.org/package=rgeos},
      }
    
    To cite RNeXML in publications, please use:
    
    Boettiger C, Chamberlain S, Vos R, Lapp H (2016). "RNeXML: A Package
    for Reading and Writing Richly Annotated Phylogenetic, Character, and
    Trait Data in R." _Methods in Ecology and Evolution_, *7*, 352-357.
    doi: 10.1111/2041-210X.12469 (URL:
    https://doi.org/10.1111/2041-210X.12469).
    
    A BibTeX entry for LaTeX users is
    
      @Article{,
        title = {{RNeXML}: {A} Package for Reading and Writing Richly Annotated Phylogenetic, Character, and Trait Data in {R}},
        journal = {Methods in Ecology and Evolution},
        author = {Carl Boettiger and Scott Chamberlain and Rutger Vos and Hilmar Lapp},
        year = {2016},
        volume = {7},
        pages = {352--357},
        doi = {10.1111/2041-210X.12469},
      }
    
      H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
      Springer-Verlag New York, 2016.
    
    A BibTeX entry for LaTeX users is
    
      @Book{,
        author = {Hadley Wickham},
        title = {ggplot2: Elegant Graphics for Data Analysis},
        publisher = {Springer-Verlag New York},
        year = {2016},
        isbn = {978-3-319-24277-4},
        url = {https://ggplot2.tidyverse.org},
      }
    

# Parse CITATION_rmarkdown

    
    To cite the 'rmarkdown' package in publications, please use:
    
      (2022). rmarkdown: A Basic Description. R package version 0.1.6. URL
      https://rmarkdown.rstudio.com.
    
    A BibTeX entry for LaTeX users is
    
      @Manual{,
        title = {rmarkdown: A Basic Description},
        year = {2022},
        note = {R package version 0.1.6},
        url = {https://github.com/basic/package},
      }
    
      Yihui Xie and J.J. Allaire and Garrett Grolemund (2018). R Markdown:
      The Definitive Guide. Chapman and Hall/CRC. ISBN 9781138359338. URL
      https://bookdown.org/yihui/rmarkdown.
    
    A BibTeX entry for LaTeX users is
    
      @Book{,
        title = {R Markdown: The Definitive Guide},
        author = {Yihui Xie and J.J. Allaire and Garrett Grolemund},
        publisher = {Chapman and Hall/CRC},
        address = {Boca Raton, Florida},
        year = {2018},
        note = {ISBN 9781138359338},
        url = {https://bookdown.org/yihui/rmarkdown},
      }
    
      Yihui Xie and Christophe Dervieux and Emily Riederer (2020). R
      Markdown Cookbook. Chapman and Hall/CRC. ISBN 9780367563837. URL
      https://bookdown.org/yihui/rmarkdown-cookbook.
    
    A BibTeX entry for LaTeX users is
    
      @Book{,
        title = {R Markdown Cookbook},
        author = {Yihui Xie and Christophe Dervieux and Emily Riederer},
        publisher = {Chapman and Hall/CRC},
        address = {Boca Raton, Florida},
        year = {2020},
        note = {ISBN 9780367563837},
        url = {https://bookdown.org/yihui/rmarkdown-cookbook},
      }
    

# Article to bibtex

    @Article{knuth:1984,
      author = {{R Core Team}},
      title = {Literate Programming},
      journal = {The Computer Journal},
      year = {1984},
      volume = {27},
      number = {2},
      pages = {97--111},
      month = {January},
      keywords = {Some},
      keywords = {simple},
      keywords = {keywords},
    }

---

    @Article{rcoreteam:1984,
      title = {Literate Programming},
      author = {{R Core Team}},
      year = {1984},
      month = {jan},
      journal = {The Computer Journal},
      volume = {27},
      number = {2},
      pages = {97--111},
      keywords = {Some, simple, keywords},
    }

# Book to bibtex

    @Book{latex:companion,
      author = {Frank Mittelbach and Michel Gossens and Johannes Braams and David Carlisle and Chris Rowley},
      editor = {{{Barnes} and {Noble}}},
      title = {The LaTeX Companion},
      publisher = {Addison-Wesley Professional},
      year = {2004},
      volume = {3},
      number = {7},
      series = {The LateX Books},
      address = {Santa Monica},
      edition = {Fourth},
      month = {August},
      note = {Example modified for testing purposes},
      keywords = {Two, keyword},
    }

---

    @Book{mittelbachgossens:2004,
      title = {The LaTeX Companion},
      author = {Frank Mittelbach and Michel Gossens and Johannes Braams and David Carlisle and Chris Rowley},
      year = {2004},
      month = {aug},
      publisher = {Addison-Wesley Professional},
      address = {Santa Monica},
      editor = {{Barnes and Noble}},
      series = {The LateX Books},
      volume = {3},
      number = {7},
      note = {Example modified for testing purposes},
      edition = {Fourth},
      keywords = {Two, keyword},
    }

# Booklet to bibtex

    @Booklet{Mustermann2016,
      title = {Java Booklet},
      author = {Max Mustermann},
      howpublished = {Internet},
      address = {Stuttgart},
      month = {feb},
      year = {2016},
      note = {Example modified from Jabref},
      keywords = {java},
    }

---

    @Booklet{mustermann:2016,
      title = {Java Booklet},
      author = {Max Mustermann},
      year = {2016},
      month = {feb},
      address = {Stuttgart},
      note = {Example modified from Jabref},
      howpublished = {Internet},
    }

# InBook to bibtex with pages

    @InBook{,
      year = {2003},
      month = {oct},
      pages = {175--196},
      title = {Architectural Mismatch Tolerance},
      chapter = {Tolerances and Other Notes},
      author = {R. {de Lemos} and C. Gacek and A. Romanovsky},
      url = {http://www.cs.kent.ac.uk/pubs/2003/1773},
      publication_type = {inbook},
      submission_id = {12884_1074884456},
      isbn = {3-540-40727-8},
      editor = {A. Lalanda},
      edition = {Fifth},
      publisher = {Springer},
      volume = {2677},
      number = {234},
      address = {Lozoya},
      series = {Lecture Notes in Computer Science},
      type = {Architecting Dependable Systems},
    }

---

    @InBook{delemosgacek:2003,
      title = {Architectural Mismatch Tolerance},
      author = {R. {de Lemos} and C. Gacek and A. Romanovsky},
      year = {2003},
      month = {oct},
      publisher = {Springer},
      address = {Lozoya},
      editor = {A. Lalanda},
      series = {Lecture Notes in Computer Science},
      volume = {2677},
      number = {234},
      pages = {175--196},
      isbn = {3-540-40727-8},
      url = {http://www.cs.kent.ac.uk/pubs/2003/1773},
      chapter = {Tolerances and Other Notes},
      edition = {Fifth},
    }

# InCollection to bibtex

    @InCollection{,
      author = {Klaus Abels},
      title = {Who Gives a Damn about Minimizers in Questions?},
      booktitle = {Proceedings from Semantics and Linguistic Theory {XIII}},
      publisher = {Cornell University},
      year = {2003},
      editor = {Robert B. Young and Yuping Zhou},
      pages = {1--18},
      address = {Ithaca, New York},
      topic = {interrogatives;nl-semantics;polarity;},
    }

---

    @InCollection{abels:2003,
      title = {Who Gives a Damn about Minimizers in Questions?},
      author = {Klaus Abels},
      year = {2003},
      booktitle = {Proceedings from Semantics and Linguistic Theory XIII},
      publisher = {Cornell University},
      address = {Ithaca, New York},
      editor = {Robert B. Young and Yuping Zhou},
      pages = {1--18},
    }

# InProceedings to bibtex

    @InProceedings{,
      author = {John Aberdeen and Samuel Bayer and Sasha Caskey and Laurie Damianos and Alan Goldschen and Lynette Hirschman and Dan Loehr and Hugo Trapper},
      title = {Implementing Practical Dialogue Systems with the
                        {DARPA} Communicator Architecture},
      booktitle = {Proceedings of the {IJCAI}-99 Workshop on
                        Knowledge and Reasoning in Practical Dialogue Systems},
      year = {1999},
      editor = {Jan Alexandersson},
      pages = {81--86},
      series = {A Series},
      organization = {IJCAI},
      publisher = {International Joint Conference on Artificial Intelligence},
      address = {Murray Hill, New Jersey},
    }

---

    @InProceedings{aberdeenbayer:1999,
      title = {Implementing Practical Dialogue Systems with the DARPA Communicator Architecture},
      author = {John Aberdeen and Samuel Bayer and Sasha Caskey and Laurie Damianos and Alan Goldschen and Lynette Hirschman and Dan Loehr and Hugo Trapper},
      year = {1999},
      booktitle = {Proceedings of the IJCAI-99 Workshop on Knowledge and Reasoning in Practical Dialogue Systems},
      publisher = {International Joint Conference on Artificial Intelligence},
      address = {Murray Hill, New Jersey},
      editor = {Jan Alexandersson},
      series = {A Series},
      pages = {81--86},
      organization = {IJCAI},
    }

# Manual to bibtex

    @Manual{,
      author = {Gerhard Allwein and Dave Barker-Plummer and Jon Barwise and John Etchemendy},
      title = {{LPL} Software Manual},
      publisher = {{CSLI} Publications},
      year = {1999},
      address = {Stanford, California},
      howpublished = {CD-Rom},
    }

---

    @Manual{allweinbarker?plummer:1999,
      title = {LPL Software Manual},
      author = {Gerhard Allwein and Dave Barker-Plummer and Jon Barwise and John Etchemendy},
      year = {1999},
      publisher = {CSLI Publications},
      address = {Stanford, California},
      howpublished = {CD-Rom},
    }

# MastersThesis to bibtex

    @MastersThesis{,
      author = {Murat Bayraktar},
      title = {Computer-Aided Analysis of {E}nglish Punctuation on a
                        Parsed Corpus: The Special Case of Comma},
      school = {Department of Computer Engineering and Information
                        Science, Bilkent University, Turkey},
      address = {Ankara, Turkey},
      year = {1996},
      note = {Forthcoming},
    }

---

    @MastersThesis{bayraktar:1996,
      title = {Computer-Aided Analysis of English Punctuation on a Parsed Corpus: The Special Case of Comma},
      author = {Murat Bayraktar},
      year = {1996},
      note = {Forthcoming},
      school = {Department of Computer Engineering and Information Science, Bilkent University, Turkey},
    }

# PhdThesis to bibtex

    @PhdThesis{,
      author = {David I. Beaver},
      title = {Presupposition and Assertion in Dynamic Semantics},
      school = {Centre for Cognitive Science, University of Edinburgh},
      year = {1995},
      type = {Ph.D. Dissertation},
      address = {Edinburgh},
    }

---

    @PhdThesis{beaver:1995,
      title = {Presupposition and Assertion in Dynamic Semantics},
      author = {David I. Beaver},
      year = {1995},
      school = {Centre for Cognitive Science, University of Edinburgh},
    }

# Proceedings to bibtex

    @Proceedings{,
      title = {An Abductive Framework for Negation in Disjunctive
                        Logic Programming},
      organization = {{JELIA}'96},
      year = {1996},
      editor = {Jose Julio Alferes and Luis Moniz Pereira and Eva Orlowska},
      publisher = {Springer-Verlag},
      address = {Berlin},
      missinginfo = {pages},
    }

---

    @Proceedings{anonymous:1996,
      title = {An Abductive Framework for Negation in Disjunctive Logic Programming},
      year = {1996},
      publisher = {Springer-Verlag},
      address = {Berlin},
      editor = {Jose Julio Alferes and Luis Moniz Pereira and Eva Orlowska},
      organization = {JELIA'96},
    }

# TechReport to bibtex

    @TechReport{,
      author = {John M. Aronis},
      title = {Implementing Inheritance on the Connection Machine},
      institution = {Intelligent Systems Program, University of Pittsburgh},
      number = {ISP 93-1},
      year = {1993},
      address = {Pittsburgh, PA 15260},
    }

---

    @TechReport{aronis:1993,
      title = {Implementing Inheritance on the Connection Machine},
      author = {John M. Aronis},
      year = {1993},
      address = {Pittsburgh, PA 15260},
      number = {ISP 93-1},
      institution = {Intelligent Systems Program, University of Pittsburgh},
    }

# Unpublished to bibtex

    @Unpublished{,
      author = {John M. Aronis and Foster J. Provost},
      title = {Efficiently Constructing Relational Features from Background},
      year = {1959},
      note = {Unpublished MS, Computer Science Department, University of Pittsburgh.},
      missinginfo = {Date is  guess.},
    }

---

    @Unpublished{aronisprovost:1959,
      title = {Efficiently Constructing Relational Features from Background},
      author = {John M. Aronis and Foster J. Provost},
      year = {1959},
      note = {Unpublished MS, Computer Science Department, University of Pittsburgh.},
    }

# particle names

    type: book
    title: A Handbook for Scholars
    authors:
    - family-names: Leunen
      given-names: Mary-Claire
      name-particle: van
    - family-names: Davis
      given-names: Sammy
      name-suffix: Jr.
    year: '1979'
    publisher:
      name: Knopf

---

    @Book{vanleunendavisjr?:1979,
      title = {A Handbook for Scholars},
      author = {Mary-Claire {van Leunen} and Sammy {Davis Jr.}},
      year = {1979},
      publisher = {Knopf},
    }

# From plain cff with pref citation

    @Manual{doe:2022,
      title = {My Research Software},
      author = {John Doe},
      year = {2022},
      month = {mar},
      editor = {A name},
    }

# From plain cff

    @Misc{doe,
      title = {My Research Software},
      author = {John Doe},
    }

# From file

    @InBook{vanderrealpersoniventityprojectteamconferenceentity:2017,
      title = {Book Title},
      author = {One Truly {van der Real Person IV} and {Entity Project Team Conference entity}},
      year = {2017},
      month = {mar},
      journal = {PeerJ},
      publisher = {Entity Project Team Conference entity},
      address = {22 Acacia Avenue, Citationburgh, Renfrewshire, GB},
      editor = {One Truly {van der Real Person IV} and {Entity Project Team Conference entity}},
      series = {Collection Title},
      volume = {2},
      number = {123},
      pages = {123--456},
      doi = {10.5281/zenodo.1003150},
      isbn = {978-1-89183-044-0},
      issn = {1234-543X},
      url = {http://j.mp},
      note = {A field for general notes about the reference, usable in other formats such as BibTeX.},
      chapter = {Chapter 2 - "Reference keys"},
      edition = {2nd edition},
      howpublished = {Hardcover book},
      institution = {Entity Project Team Conference entity},
      keywords = {Software, Citation},
      abstract = {Description of the book.},
      date = {2017-10-31},
      file = {book.zip},
      issuetitle = {Special Issue on Software Citation},
      pagetotal = {765},
      urldate = {2017-10-31},
      version = {0.0.1423-BETA},
      translator = {van der Real Person, IV, One Truly and {Entity Project Team Conference entity}},
    }

# Test anonymous

    @Booklet{anonymous,
      title = {A booklet},
    }

---

    @Manual{anonymous,
      title = {A manual},
    }

---

    @Misc{anonymous,
      title = {A misc},
    }

---

    @Proceedings{anonymous:1984,
      title = {proceedings},
      year = {1984},
    }

# Fallback month

    @Article{,
      title = {An Article},
      author = {John Doe},
      journal = {El Adelantado de Segovia},
      year = {1678},
      date = {1678-04-23},
    }

---

    @Article{doe:1678,
      title = {An Article},
      author = {John Doe},
      year = {1678},
      month = {apr},
      journal = {El Adelantado de Segovia},
      date = {1678-04-23},
    }

# Test BibLateX entry

    @Article{,
      author = {M. A. Kastenholz and Philippe H. Hünenbergerb},
      title = {Computation of methodology hyphen independent ionic solvation
                      free energies from molecular simulations},
      journal = {J. Chem. Phys.},
      year = {2006},
      note = {Example modified for testing purposes},
      pages = {55--65},
      date = {2006-03-15},
      file = {a_file.pdf},
      issuetitle = {Semantic {3D} Media and Content},
      translator = {Wicksteed, P. H. and {The translator factory}},
      urldate = {2006-10-01},
      pagetotal = {528},
      abstract = {The computation of ionic solvation free energies from
                      atomistic simulations is a surprisingly difficult problem that
                      has found no satisfactory solution for more than 15 years.},
      doi = {10.1063/1.2172593},
      isbn = {0-816-52066-6},
      issn = {0097-8493},
      url = {http://www.ctan.org},
    }

---

    @Article{kastenholzh?nenbergerb:2006,
      title = {Computation of methodology hyphen independent ionic solvation free energies from molecular simulations},
      author = {M. A. Kastenholz and Philippe H. Hünenbergerb},
      year = {2006},
      month = {mar},
      journal = {J. Chem. Phys.},
      pages = {55--65},
      doi = {10.1063/1.2172593},
      isbn = {0-816-52066-6},
      issn = {0097-8493},
      url = {http://www.ctan.org},
      note = {Example modified for testing purposes},
      abstract = {The computation of ionic solvation free energies from atomistic simulations is a surprisingly difficult problem that has found no satisfactory solution for more than 15 years.},
      date = {2006-03-15},
      file = {a_file.pdf},
      issuetitle = {Semantic 3D Media and Content},
      pagetotal = {528},
      urldate = {2006-10-01},
      translator = {Wicksteed, P. H. and {The translator factory}},
    }

# Add new keys

    cff-version: 1.2.0
    message: This overwrites fields
    type: software
    license: GPL-3.0-only
    title: 'basicdesc: A Basic Description'
    version: 0.1.6
    abstract: New abstract
    authors:
    - family-names: Nadie
      given-names: Don
    preferred-citation:
      type: manual
      title: 'basicdesc: A Basic Description'
      authors:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      version: 0.1.6
      abstract: A very basic description. Should parse without problems.
      repository-code: https://github.com/basic/package
      url: https://basic.github.io/package
      contact:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      license: GPL-3.0-only
      year: '2022'
    repository-code: https://github.com/basic/package
    url: https://basic.github.io/package
    date-released: '1900-01-01'
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    keywords:
    - A
    - new
    - list
    - of
    - keywords

# Append keys

    cff-version: 1.2.0
    message: 'To cite package "basicdesc" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'basicdesc: A Basic Description'
    version: 0.1.6
    abstract: A very basic description. Should parse without problems.
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    preferred-citation:
      type: manual
      title: 'basicdesc: A Basic Description'
      authors:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      version: 0.1.6
      abstract: A very basic description. Should parse without problems.
      repository-code: https://github.com/basic/package
      url: https://basic.github.io/package
      contact:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      license: GPL-3.0-only
      year: '2022'
    repository-code: https://github.com/basic/package
    url: https://basic.github.io/package
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com

---

    cff-version: 1.2.0
    message: 'To cite package "basicdesc" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'basicdesc: A Basic Description'
    version: 0.1.6
    abstract: A very basic description. Should parse without problems.
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    - family-names: author
      given-names: New
      website: https://stackoverflow.com/
      country: IT
    preferred-citation:
      type: manual
      title: 'basicdesc: A Basic Description'
      authors:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      version: 0.1.6
      abstract: A very basic description. Should parse without problems.
      repository-code: https://github.com/basic/package
      url: https://basic.github.io/package
      contact:
      - family-names: Basic
        given-names: Marc
        email: marcbasic@gmail.com
      license: GPL-3.0-only
      year: '2022'
    repository-code: https://github.com/basic/package
    url: https://basic.github.io/package
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com

# Try not writing

    @Misc{,
      title = {My title},
      author = {Fran Pérez},
    }

---

    @Misc{,
      title = {My title},
      author = {Fran P{\'e}rez},
    }

# Walk trough full lifecycle

    
    
    ## Read object 
    
    cff-version: 1.2.0
    message: If you use this software, please cite it as below.
    abstract: This is an awesome piece of research software!
    authors:
    - family-names: Real Person
      given-names: One Truly
      name-particle: van der
      name-suffix: IV
      alias: Citey
      affiliation: Excellent University, Niceplace, Arcadia
      address: 22 Acacia Avenue
      city: Citationburgh
      region: Renfrewshire
      post-code: C13 7X7
      country: GB
      orcid: https://orcid.org/0000-0001-2345-6789
      email: project@entity.com
      tel: +44(0)141-323 4567
      fax: +44(0)141-323 45678
      website: https://www.entity-project-team.io
    - name: Entity Project Team Conference entity
      address: 22 Acacia Avenue
      city: Citationburgh
      region: Renfrewshire
      post-code: C13 7X7
      country: GB
      orcid: https://orcid.org/0000-0001-2345-6789
      email: project@entity.com
      tel: +44(0)141-323 4567
      fax: +44(0)141-323 45678
      website: https://www.entity-project-team.io
      date-start: '2017-01-01'
      date-end: '2017-01-31'
      location: The team garage
    commit: 156a04c74a8a79d40c5d705cddf9d36735feab4d
    contact:
    - family-names: Real Person
      given-names: One Truly
      name-particle: van der
      name-suffix: IV
      alias: Citey
      affiliation: Excellent University, Niceplace, Arcadia
      address: 22 Acacia Avenue
      city: Citationburgh
      region: Renfrewshire
      post-code: C13 7X7
      country: GB
      orcid: https://orcid.org/0000-0001-2345-6789
      email: project@entity.com
      tel: +44(0)141-323 4567
      fax: +44(0)141-323 45678
      website: https://www.entity-project-team.io
    - name: Entity Project Team Conference entity
      address: 22 Acacia Avenue
      city: Citationburgh
      region: Renfrewshire
      post-code: C13 7X7
      country: GB
      orcid: https://orcid.org/0000-0001-2345-6789
      email: project@entity.com
      tel: +44(0)141-323 4567
      fax: +44(0)141-323 45678
      website: https://www.entity-project-team.io
      date-start: '2017-01-01'
      date-end: '2017-01-31'
      location: The team garage
    date-released: '2017-12-11'
    doi: 10.5281/zenodo.1003150
    identifiers:
    - type: doi
      value: 10.5281/zenodo.1003150
    - type: swh
      value: swh:1:rel:99f6850374dc6597af01bd0ee1d3fc0699301b9f
    - type: url
      value: https://example.com
    - type: other
      value: other-schema://abcd.1234.efgh.5678
    keywords:
    - One
    - Two
    - Three
    - '4'
    license: CC-BY-SA-4.0
    license-url: https://spdx.org/licenses/CC-BY-SA-4.0.html#licenseText
    repository: https://www.example.com/foo/?bar=baz&inga=42&quux
    repository-code: http://foo.com/blah_(wikipedia)_blah#cite-1
    repository-artifact: https://files.pythonhosted.org/packages/0a/84/10507b69a07768bc16981184b4d147a0fc84b71fbf35c03bafc8dcced8e1/cffconvert-1.3.3.tar.gz
    title: Citation File Format 1.0.0
    type: software
    url: http://userid:password@example.com:8080/
    version: 1.0.0
    preferred-citation:
      type: book
      title: Book Title
      abbreviation: Abbr
      abstract: Description of the book.
      collection-doi: 10.5281/zenodo.1003150
      collection-title: Collection Title
      collection-type: Collection Type
      commit: 156a04c74a8a79d40c5d705cddf9d36735feab4d
      copyright: 2017 Stephan Druskat
      data-type: Data Type
      database: Database
      date-accessed: '2017-10-31'
      date-downloaded: '2017-10-31'
      date-released: '2017-10-31'
      date-published: '2017-10-31'
      department: Department
      doi: 10.5281/zenodo.1003150
      edition: 2nd edition
      end: '456'
      entry: Chapter 9
      filename: book.zip
      format: Printed book
      identifiers:
      - type: doi
        value: 10.5281/zenodo.1003150
      - type: swh
        value: swh:1:rel:99f6850374dc6597af01bd0ee1d3fc0699301b9f
      - type: url
        value: https://example.com
      - type: other
        value: other-schema://abcd.1234.efgh.5678
      isbn: 978-1-89183-044-0
      issn: 1234-543X
      issue: '123'
      issue-date: December
      issue-title: Special Issue on Software Citation
      journal: PeerJ
      keywords:
      - Software
      - Citation
      languages:
      - aaa
      - zu
      license: Apache-2.0
      license-url: https://spdx.org/licenses/Apache-2.0.html#licenseText
      loc-start: '14'
      loc-end: '54'
      medium: hardcover book
      month: '3'
      nihmsid: Don't know what format a NIHMSID is in
      notes: A field for general notes about the reference, usable in other formats such
        as BibTeX.
      number: A general-purpose field for accession numbers, cf. the specifications for
        examples.
      number-volumes: '7'
      pages: '765'
      patent-states:
      - Germany
      - ROI
      - 'but also for example US states, such as:'
      - IL
      - RI
      pmcid: PMC1234567
      repository: http://code.google.com/events/#&product=browser
      repository-code: http://142.42.1.1:8080/
      repository-artifact: https://files.pythonhosted.org/packages/0a/84/10507b69a07768bc16981184b4d147a0fc84b71fbf35c03bafc8dcced8e1/cffconvert-1.3.3.tar.gz
      scope: Cite this book if you want to reference the general concepts implemented
        in Citation File Format 1.0.0.
      section: Chapter 2 - "Reference keys"
      status: advance-online
      start: '123'
      thesis-type: Doctoral dissertation
      url: http://j.mp
      version: 0.0.1423-BETA
      volume: '2'
      volume-title: Advances in Software Citation
      year: '2017'
      year-original: '2012'
      conference:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      authors:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      contact:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      database-provider:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      editors:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      editors-series:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      institution:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      location:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      publisher:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      recipients:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      senders:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      translators:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
    references:
    - type: book
      title: Book Title
      abbreviation: Abbr
      abstract: Description of the book.
      collection-doi: 10.5281/zenodo.1003150
      collection-title: Collection Title
      collection-type: Collection Type
      commit: 156a04c74a8a79d40c5d705cddf9d36735feab4d
      copyright: 2017 Stephan Druskat
      data-type: Data Type
      database: Database
      date-accessed: '2017-10-31'
      date-downloaded: '2017-10-31'
      date-released: '2017-10-31'
      date-published: '2017-10-31'
      department: Department
      doi: 10.5281/zenodo.1003150
      edition: 2nd edition
      end: '123'
      entry: Chapter 9
      filename: book.zip
      format: Printed book
      identifiers:
      - type: doi
        value: 10.5281/zenodo.1003150
      - type: swh
        value: swh:1:rel:99f6850374dc6597af01bd0ee1d3fc0699301b9f
      - type: url
        value: https://example.com
      - type: other
        value: other-schema://abcd.1234.efgh.5678
      isbn: 978-1-89183-044-0
      issn: 1234-543X
      issue: '123'
      issue-date: December
      issue-title: Special Issue on Software Citation
      journal: PeerJ
      keywords:
      - Software
      - Citation
      languages:
      - aaa
      - zu
      license: Apache-2.0
      license-url: https://spdx.org/licenses/Apache-2.0.html#licenseText
      loc-start: '14'
      loc-end: '54'
      medium: hardcover book
      month: '3'
      nihmsid: Don't know what format a NIHMSID is in
      notes: A field for general notes about the reference, usable in other formats such
        as BibTeX.
      number: A general-purpose field for accession numbers, cf. the specifications for
        examples.
      number-volumes: '7'
      pages: '765'
      patent-states:
      - Germany
      - ROI
      - 'but also for example US states, such as:'
      - IL
      - RI
      pmcid: PMC1234567
      repository: http://code.google.com/events/#&product=browser
      repository-code: http://142.42.1.1:8080/
      repository-artifact: https://files.pythonhosted.org/packages/0a/84/10507b69a07768bc16981184b4d147a0fc84b71fbf35c03bafc8dcced8e1/cffconvert-1.3.3.tar.gz
      scope: Cite this book if you want to reference the general concepts implemented
        in Citation File Format 1.0.0.
      section: Chapter 2 - "Reference keys"
      status: advance-online
      start: '123'
      thesis-type: Doctoral dissertation
      url: http://j.mp
      version: 0.0.1423-BETA
      volume: '2'
      volume-title: Advances in Software Citation
      year: '2017'
      year-original: '2012'
      conference:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      authors:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      contact:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      database-provider:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      editors:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      editors-series:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      institution:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      location:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      publisher:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      recipients:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      senders:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      translators:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
    
    ---

---

    
    
    ## Modify object 
    
    cff-version: 1.2.0
    message: If you use this software, please cite it as below.
    type: software
    license: CC-BY-SA-4.0
    title: A new title
    version: 1.0.0
    doi: 10.5281/zenodo.1003150
    abstract: This is an awesome piece of research software!
    authors:
    - family-names: Real Person
      given-names: One Truly
      name-particle: van der
      name-suffix: IV
      alias: Citey
      affiliation: Excellent University, Niceplace, Arcadia
      address: 22 Acacia Avenue
      city: Citationburgh
      region: Renfrewshire
      post-code: C13 7X7
      country: GB
      orcid: https://orcid.org/0000-0001-2345-6789
      email: project@entity.com
      tel: +44(0)141-323 4567
      fax: +44(0)141-323 45678
      website: https://www.entity-project-team.io
    - name: Entity Project Team Conference entity
      address: 22 Acacia Avenue
      city: Citationburgh
      region: Renfrewshire
      post-code: C13 7X7
      country: GB
      orcid: https://orcid.org/0000-0001-2345-6789
      email: project@entity.com
      tel: +44(0)141-323 4567
      fax: +44(0)141-323 45678
      website: https://www.entity-project-team.io
      date-start: '2017-01-01'
      date-end: '2017-01-31'
      location: The team garage
    preferred-citation:
      type: book
      title: Book Title
      abbreviation: Abbr
      abstract: Description of the book.
      collection-doi: 10.5281/zenodo.1003150
      collection-title: Collection Title
      collection-type: Collection Type
      commit: 156a04c74a8a79d40c5d705cddf9d36735feab4d
      copyright: 2017 Stephan Druskat
      data-type: Data Type
      database: Database
      date-accessed: '2017-10-31'
      date-downloaded: '2017-10-31'
      date-released: '2017-10-31'
      date-published: '2017-10-31'
      department: Department
      doi: 10.5281/zenodo.1003150
      edition: 2nd edition
      end: '456'
      entry: Chapter 9
      filename: book.zip
      format: Printed book
      identifiers:
      - type: doi
        value: 10.5281/zenodo.1003150
      - type: swh
        value: swh:1:rel:99f6850374dc6597af01bd0ee1d3fc0699301b9f
      - type: url
        value: https://example.com
      - type: other
        value: other-schema://abcd.1234.efgh.5678
      isbn: 978-1-89183-044-0
      issn: 1234-543X
      issue: '123'
      issue-date: December
      issue-title: Special Issue on Software Citation
      journal: PeerJ
      keywords:
      - Software
      - Citation
      languages:
      - aaa
      - zu
      license: Apache-2.0
      license-url: https://spdx.org/licenses/Apache-2.0.html#licenseText
      loc-start: '14'
      loc-end: '54'
      medium: hardcover book
      month: '3'
      nihmsid: Don't know what format a NIHMSID is in
      notes: A field for general notes about the reference, usable in other formats such
        as BibTeX.
      number: A general-purpose field for accession numbers, cf. the specifications for
        examples.
      number-volumes: '7'
      pages: '765'
      patent-states:
      - Germany
      - ROI
      - 'but also for example US states, such as:'
      - IL
      - RI
      pmcid: PMC1234567
      repository: http://code.google.com/events/#&product=browser
      repository-code: http://142.42.1.1:8080/
      repository-artifact: https://files.pythonhosted.org/packages/0a/84/10507b69a07768bc16981184b4d147a0fc84b71fbf35c03bafc8dcced8e1/cffconvert-1.3.3.tar.gz
      scope: Cite this book if you want to reference the general concepts implemented
        in Citation File Format 1.0.0.
      section: Chapter 2 - "Reference keys"
      status: advance-online
      start: '123'
      thesis-type: Doctoral dissertation
      url: http://j.mp
      version: 0.0.1423-BETA
      volume: '2'
      volume-title: Advances in Software Citation
      year: '2017'
      year-original: '2012'
      conference:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      authors:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      contact:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      database-provider:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      editors:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      editors-series:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      institution:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      location:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      publisher:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      recipients:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      senders:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      translators:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
    repository: https://www.example.com/foo/?bar=baz&inga=42&quux
    repository-artifact: https://files.pythonhosted.org/packages/0a/84/10507b69a07768bc16981184b4d147a0fc84b71fbf35c03bafc8dcced8e1/cffconvert-1.3.3.tar.gz
    repository-code: http://foo.com/blah_(wikipedia)_blah#cite-1
    url: http://userid:password@example.com:8080/
    date-released: '2017-12-11'
    contact:
    - family-names: Real Person
      given-names: One Truly
      name-particle: van der
      name-suffix: IV
      alias: Citey
      affiliation: Excellent University, Niceplace, Arcadia
      address: 22 Acacia Avenue
      city: Citationburgh
      region: Renfrewshire
      post-code: C13 7X7
      country: GB
      orcid: https://orcid.org/0000-0001-2345-6789
      email: project@entity.com
      tel: +44(0)141-323 4567
      fax: +44(0)141-323 45678
      website: https://www.entity-project-team.io
    - name: Entity Project Team Conference entity
      address: 22 Acacia Avenue
      city: Citationburgh
      region: Renfrewshire
      post-code: C13 7X7
      country: GB
      orcid: https://orcid.org/0000-0001-2345-6789
      email: project@entity.com
      tel: +44(0)141-323 4567
      fax: +44(0)141-323 45678
      website: https://www.entity-project-team.io
      date-start: '2017-01-01'
      date-end: '2017-01-31'
      location: The team garage
    keywords:
    - One
    - Two
    - Three
    - '4'
    references:
    - type: book
      title: Book Title
      abbreviation: Abbr
      abstract: Description of the book.
      collection-doi: 10.5281/zenodo.1003150
      collection-title: Collection Title
      collection-type: Collection Type
      commit: 156a04c74a8a79d40c5d705cddf9d36735feab4d
      copyright: 2017 Stephan Druskat
      data-type: Data Type
      database: Database
      date-accessed: '2017-10-31'
      date-downloaded: '2017-10-31'
      date-released: '2017-10-31'
      date-published: '2017-10-31'
      department: Department
      doi: 10.5281/zenodo.1003150
      edition: 2nd edition
      end: '123'
      entry: Chapter 9
      filename: book.zip
      format: Printed book
      identifiers:
      - type: doi
        value: 10.5281/zenodo.1003150
      - type: swh
        value: swh:1:rel:99f6850374dc6597af01bd0ee1d3fc0699301b9f
      - type: url
        value: https://example.com
      - type: other
        value: other-schema://abcd.1234.efgh.5678
      isbn: 978-1-89183-044-0
      issn: 1234-543X
      issue: '123'
      issue-date: December
      issue-title: Special Issue on Software Citation
      journal: PeerJ
      keywords:
      - Software
      - Citation
      languages:
      - aaa
      - zu
      license: Apache-2.0
      license-url: https://spdx.org/licenses/Apache-2.0.html#licenseText
      loc-start: '14'
      loc-end: '54'
      medium: hardcover book
      month: '3'
      nihmsid: Don't know what format a NIHMSID is in
      notes: A field for general notes about the reference, usable in other formats such
        as BibTeX.
      number: A general-purpose field for accession numbers, cf. the specifications for
        examples.
      number-volumes: '7'
      pages: '765'
      patent-states:
      - Germany
      - ROI
      - 'but also for example US states, such as:'
      - IL
      - RI
      pmcid: PMC1234567
      repository: http://code.google.com/events/#&product=browser
      repository-code: http://142.42.1.1:8080/
      repository-artifact: https://files.pythonhosted.org/packages/0a/84/10507b69a07768bc16981184b4d147a0fc84b71fbf35c03bafc8dcced8e1/cffconvert-1.3.3.tar.gz
      scope: Cite this book if you want to reference the general concepts implemented
        in Citation File Format 1.0.0.
      section: Chapter 2 - "Reference keys"
      status: advance-online
      start: '123'
      thesis-type: Doctoral dissertation
      url: http://j.mp
      version: 0.0.1423-BETA
      volume: '2'
      volume-title: Advances in Software Citation
      year: '2017'
      year-original: '2012'
      conference:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      authors:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      contact:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      database-provider:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      editors:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      editors-series:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      institution:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      location:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      publisher:
        name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      recipients:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      senders:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
      translators:
      - family-names: Real Person
        given-names: One Truly
        name-particle: van der
        name-suffix: IV
        alias: Citey
        affiliation: Excellent University, Niceplace, Arcadia
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
      - name: Entity Project Team Conference entity
        address: 22 Acacia Avenue
        city: Citationburgh
        region: Renfrewshire
        post-code: C13 7X7
        country: GB
        orcid: https://orcid.org/0000-0001-2345-6789
        email: project@entity.com
        tel: +44(0)141-323 4567
        fax: +44(0)141-323 45678
        website: https://www.entity-project-team.io
        date-start: '2017-01-01'
        date-end: '2017-01-31'
        location: The team garage
    commit: 156a04c74a8a79d40c5d705cddf9d36735feab4d
    identifiers:
    - type: doi
      value: 10.5281/zenodo.1003150
    - type: swh
      value: swh:1:rel:99f6850374dc6597af01bd0ee1d3fc0699301b9f
    - type: url
      value: https://example.com
    - type: other
      value: other-schema://abcd.1234.efgh.5678
    license-url: https://spdx.org/licenses/CC-BY-SA-4.0.html#licenseText
    
    ---

# Fuzzy matching of keys on cff

    
    
    ## Fuzzy keys 
    
    title: a
    cff-version: 1.2.0
    version: '200'
    authors:
    - family-names: a
      given-names: b
    message: Fix my keys
    
    ---

# Test in mock package

    cff-version: 1.2.0
    message: 'To cite package "manyurls" in publications use:'
    type: software
    license: GPL-3.0-only
    title: 'manyurls: A lot of urls'
    version: 0.1.6
    doi: 10.1111/2041-210X.12469
    abstract: This package has many urls. Specifically, 1 Bug Reports and 6 URLs. Expected
      is to have 1 repository-code, 1 url and 3 URLs, since there is 1 duplicate and 1
      invalid url.
    authors:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    preferred-citation:
      type: article
      title: 'RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic,
        Character, and Trait Data in R'
      authors:
      - family-names: Boettiger
        given-names: Carl
      - family-names: Chamberlain
        given-names: Scott
      - family-names: Vos
        given-names: Rutger
      - family-names: Lapp
        given-names: Hilmar
      journal: Methods in Ecology and Evolution
      year: '2016'
      volume: '7'
      doi: 10.1111/2041-210X.12469
      start: '352'
      end: '357'
    repository-code: https://github.com/test/package
    url: https://test.github.io/package/
    contact:
    - family-names: Basic
      given-names: Marc
      email: marcbasic@gmail.com
    references:
    - type: book
      title: 'ggplot2: Elegant Graphics for Data Analysis'
      authors:
      - family-names: Wickham
        given-names: Hadley
      publisher:
        name: Springer-Verlag New York
      year: '2016'
      isbn: 978-3-319-24277-4
      url: https://ggplot2.tidyverse.org
    identifiers:
    - type: url
      value: https://r-forge.r-project.org/projects/test/
    - type: url
      value: http://google.ru
    - type: url
      value: https://gitlab.com/r-packages/behaviorchange

# Validate error CITATION.cff

                field                          message
    1            data        has additional properties
    2  data.authors.0                 no schemas match
    3        data.doi referenced schema does not match
    4 data.keywords.0                is the wrong type
    5    data.license referenced schema does not match
    6        data.url referenced schema does not match

# Parse one person

    family-names: person
    given-names: one

# Parse several persons

    - family-names: person
      given-names: one
    - family-names: human
      given-names: another
    - family-names: more
      given-names: and one

# Parse bibtex persons

    family-names: Wright
    given-names: Frank Edwin
    name-suffix: III

---

    - family-names: person
      given-names: A
    - name: another
    - family-names: one
      given-names: Another

# Parse bibtex persons with masks

    - name: Elephant
    - name: Castle

---

    name: Elephant and Castle

---

    - name: Elephant and Castle
    - name: this
    - name: Ltd.

---

    - name: Elephant and Castle
    - name: this AND Ltd.

# Test first von last

     given family 
      "AA"   "BB" 

---

    family 
      "AA" 

---

     given family 
      "AA"   "bb" 

---

    family 
      "aa" 

---

     given    von family 
      "AA"   "bb"   "CC" 

---

         given        von     family 
          "AA" "bb CC dd"       "EE" 

---

      given     von  family 
    "AA bB"    "cc"    "dd" 

---

         given        von     family 
    "AA \\BBb"       "cc"       "dd" 

---

      given     von  family 
    "AA bb"    "cc"    "DD" 

---

      given     von  family 
       "AA"    "bb" "cc DD" 

---

      given  family 
    "AA bb"    "CC" 

# Testing with random names First von Last

            given           von        family 
           "Jean"          "de" "La Fontaine" 

---

               given           family 
             "Diego" "Hernandez Sanz" 

---

            given        family 
    "Juan Manuel"  "Miramontes" 

---

                  given              family 
          "Juan Manuel" "Miramontes Garcia" 

---

            given           von        family 
    "Juan Manuel"         "van"       "Halen" 

---

                   given               family 
                 "Bosco" "de la Cruz y Ochoa" 

# Test von Last, First

     given    von family 
      "AA"   "bb"   "CC" 

---

     given    von family 
      "aa"   "bb"   "CC" 

---

         given        von     family 
          "AA" "bb CC dd"       "EE" 

---

     given family 
      "AA"   "bb" 

---

    family 
      "BB" 

# Test von Last, First with brackets, etc

      given     von  family 
      "Ana"    "de" "Armas" 

---

         given     family 
         "Ana" "de Armas" 

---

                 given             family 
                 "Ana" "de Armas, Aguero" 

---

                 given             family 
           "Ana Maria" "de Armas, Aguero" 

# Test von Last, Jr,  First

     given    von family     jr 
      "AA"   "bb"   "CC"   "XX" 

---

     given family 
      "AA"   "BB" 

---

    family     jr 
      "BB"   "AA" 

# Test von Last, Jr,  First with masking

      given  family      jr 
    "Sammy" "Davis"    "Jr" 

---

            given        family            jr 
          "Sammy"  "Davis, and" "Jr, another" 

# Rest of cases

                              family 
    "David, and, Jr, another, Sammy" 

# tames da beast

             von       family 
    "jean de la"   "fontaine" 

---

         given        von     family 
        "Jean"    "de la" "fontaine" 

---

         given        von     family 
     "Jean de"       "la" "fontaine" 

---

                 von           family 
              "jean" "de la fontaine" 

---

           given       family 
    "Jean de la"   "fontaine" 

---

           given       family 
    "Jean De La"   "Fontaine" 

---

             von       family 
    "jean De la"   "Fontaine" 

---

            given           von        family 
           "Jean"          "de" "La Fontaine" 

---

             von       family 
    "jean de la"   "fontaine" 

---

         given        von     family 
        "Jean"    "de la" "fontaine" 

---

               given           family 
              "Jean" "De La Fontaine" 

---

         given        von     family 
        "Jean"    "De la" "Fontaine" 

---

            given           von        family 
           "Jean"          "de" "La Fontaine" 

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
---
output: github_document
bibliography: inst/REFERENCES.bib
link-citations: yes
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# cffr <a href='https://docs.ropensci.org/cffr/'><img src="man/figures/logo.png" align="right" height="139"/></a>

<!-- badges: start -->

[![CRAN-status](https://www.r-pkg.org/badges/version/cffr)](https://CRAN.R-project.org/package=cffr)
[![CRAN-results](https://cranchecks.info/badges/worst/cffr)](https://cran.r-project.org/web/checks/check_results_cffr.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/cffr?color=blue)](https://cran.r-project.org/package=cffr)
[![R-CMD-check](https://github.com/ropensci/cffr/actions/workflows/check-full.yaml/badge.svg)](https://github.com/ropensci/cffr/actions/workflows/check-full.yaml)
[![codecov](https://codecov.io/gh/ropensci/cffr/branch/main/graph/badge.svg?token=YRO3XL8RWK)](https://app.codecov.io/gh/ropensci/cffr)
[![r-universe](https://ropensci.r-universe.dev/badges/cffr)](https://ropensci.r-universe.dev/)
[![CITATION-cff](https://github.com/ropensci/cffr/actions/workflows/cff-validator.yml/badge.svg)](https://github.com/ropensci/cffr/actions/workflows/cff-validator.yml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03900/status.svg)](https://doi.org/10.21105/joss.03900)
[![Project Status: Active - The project has reached a stable, usable state and
is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub code size in
bytes](https://img.shields.io/github/languages/code-size/ropensci/cffr)
[![peer-review](https://badges.ropensci.org/463_status.svg)](https://github.com/ropensci/software-review/issues/463)

<!-- badges: end -->

**cffr** provides utilities to generate, parse, modify and validate
`CITATION.cff` files automatically for **R** packages, as well as tools and
examples for working with .cff more generally.

## What is a `CITATION.cff` file?

[Citation File Format (CFF](https://citation-file-format.github.io/))
[@druskat_citation_2021] (v1.2.0) are plain text files with human- and
machine-readable citation information for software (and datasets). Code
developers can include them in their repositories to let others know how to
correctly cite their software.

This format is becoming popular within the software citation ecosystem. Recently
[GitHub](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files),
[Zenodo](https://twitter.com/ZENODO_ORG/status/1420357001490706442) and
[Zotero](https://twitter.com/zotero/status/1420515377390530560) have included
full support of this citation format [@druskat_stephan_making_2021]. GitHub
support is of special interest:

```{r echo=FALSE, out.width="400", fig.align='center', fig.alt="GitHub-link"}
knitr::include_graphics("vignettes/tweet-1.png")
```

*--- Nat Friedman (\@natfriedman) [July 27,
2021](https://twitter.com/natfriedman/status/1420122675813441540?ref_src=twsrc%5Etfw)*

See [Enhanced support for citations on
GitHub](https://github.blog/2021-08-19-enhanced-support-citations-github/)
[@smith2021] for more info.

### Related projects

[The CodeMeta Project](https://codemeta.github.io/) [@jones2017] creates a
concept vocabulary that can be used to standardize the exchange of software
metadata across repositories and organizations. One of the many uses of a
`codemeta.json` file (created following the standards defined on The CodeMeta
Project) is to provide citation metadata such as title, authors, publication
year, and venue [@fenner2021]. The packages
[**codemeta**](https://github.com/cboettig/codemeta)/
[**codemetar**](https://github.com/ropensci/codemetar) allows to generate
`codemeta.json` files from R packages metadata.

## The cffr package

**cffr** maximizes the data extraction by using both the `DESCRIPTION` file and
the `CITATION` file (if present) of your package. Note that **cffr** works best
if your package pass 
`R CMD check/devtools::check()`.

```{r count_cffr, echo=FALSE, results='asis'}
cat("\n")
today <- Sys.Date()
# Try get the count of GitHub repos here
token <- (Sys.getenv(c("GITHUB_PAT", "GITHUB_TOKEN")))
token <- token[!token %in% c(NA, NULL, "")][1]
ghtoken <- paste("token", token)
tmpfile <- tempfile(fileext = ".json")
# Get numbers of repos
api_url <- "https://api.github.com/search/code?q=cffr+extension:cff+filename:CITATION"
res <- tryCatch(download.file(api_url,
  tmpfile,
  quiet = TRUE,
  headers = c(Authorization = ghtoken)
),
warning = function(e) {
  return(TRUE)
},
error = function(e) {
  return(TRUE)
}
)
# If not successful
if (isTRUE(res)) {
  cat(paste0(
    "\n", "See [some projects already using **cffr**]",
    "(https://github.com/search?l=&o=desc&q=cffr+extension%3Acff+filename%3ACITATION&s=indexed&type=Code)",
    "."
  ))
} else {
  nreps <- as.integer(jsonlite::read_json(tmpfile)$total_count)
  cat(paste0(
    "As per ", today, " there are at least ", nreps, " repos on GitHub using **cffr**. ",
    "[Check them out here]",
    "(https://github.com/search?l=&o=desc&q=cffr+extension%3Acff+filename%3ACITATION&s=indexed&type=Code)."
  ))
}
cat("\n")
```

### Installation

Install **cffr** from [CRAN](https://CRAN.R-project.org/package=cffr):

```{r, eval=FALSE}
install.packages("cffr")
```

You can install the developing version of **cffr** with:

```{r, eval=FALSE}
devtools::install_github("ropensci/cffr")
```

Alternatively, you can install **cffr** using the
[r-universe](https://ropensci.r-universe.dev/ui#builds):

```{r, eval=FALSE}

# Enable this universe
options(repos = c(
  ropensci = "https://ropensci.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

# Install some packages
install.packages("cffr")
```

### Example

By default most often from within your package folder you'll simply run
`cff_write()`, that creates a `cff` object, write it on a `CITATION.cff` file
and validates it on a single command:

```{r, eval=FALSE}

library(cffr)

# For in-development packages
cff_write()
#>
#> CITATION.cff generated
#>
#> cff_validate results-----
#> Congratulations! This .cff file is valid
```

However, **cffr** provides also custom print methods and mechanisms that allows
you to customize the `CITATION.cff` and integrate them in your workflows.

This is a basic example which shows you how to create a `cff` object (see `?cff`
for more info). In this case, we are creating a `cff` object from the metadata
of the **rmarkdown** package:

```{r }
library(cffr)

# Example with an installed package
test <- cff_create("rmarkdown")
```

<details>
<summary><code>CITATION.cff</code> for <strong>rmarkdown</strong>
</summary>

```{r, echo=FALSE, comment=""}

test
```

</details>

<p>

We can validate the result using `cff_validate()`:

```{r }

cff_validate(test)
```

Check the [docs](https://docs.ropensci.org/cffr/reference/index.html) and
`vignette("cffr", package = "cffr")` to learn how to work with `cff` objects.

### Keep your `CITATION.cff` file up-to-date

#### GitHub Actions

The easiest way for keeping you `CITATION.cff` file up-to-date is using GitHub
Actions. Use `cff_gha_update()`function to install a GitHub Action that would
update your `CITATION.cff` file on the following events:

-   When you publish a new release of the package on your GitHub repo.
-   Each time that you modify your DESCRIPTION or inst/CITATION files.
-   The action can be run also manually.

```{r, eval=FALSE}
cff_gha_update()

#> Installing update-citation-cff.yaml on './.github/workflows'
#> Adding .github to .Rbuildignore
```

See the example workflow file
[here](https://github.com/ropensci/cffr/blob/main/.github/workflows/update-citation-cff.yaml).

#### Git pre-commit hook [![Experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

You can also use a [git pre-commit
hook](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks#_committing_workflow_hooks):

> The `pre-commit` hook is run first, before you even type in a commit message.
> It's used to inspect the snapshot that's about to be committed, to see if
> you've forgotten something, to make sure tests run, or to examine whatever you
> need to inspect in the code. Exiting non-zero from this hook aborts the
> commit, although you can bypass it with `git commit --no-verify`.

A specific pre-commit hook can be installed with `cff_git_hook_install()`. If
you want to use a pre-commit hook, please make sure you have the **testthat**
package installed.

### Learn more

Check the following articles to learn more about **cffr**:

-   [cffr: Create a CITATION.cff File for your R
    Package](https://ropensci.org/blog/2021/11/23/cffr/)
-   [How I Test cffr on (about) 2,000 Packages using GitHub Actions and
    R-universe](https://ropensci.org/blog/2021/11/23/how-i-test-cffr/)

## Related packages

-   [**citation**](https://github.com/pik-piam/citation/): The development
    version (at the time of this writing) includes a new function `r2cff` that
    creates a `CITATION.cff` file (v1.1.0) using the information of your
    `DESCRIPTION` file. It also provide minimal validity checks.
-   [**handlr**](https://github.com/ropensci/handlr): Tool for converting among
    citation formats, including `*.cff` files. At the time of this writing only
    CFF v1.1.0 was supported (see
    [#24](https://github.com/ropensci/handlr/issues/24)).
-   [**codemeta**](https://github.com/cboettig/codemeta)/
    [**codemetar**](https://github.com/ropensci/codemetar) provides similar
    solutions for creating `codemeta.json` file, another format for storing and
    sharing software metadata.

## Citation

```{r echo=FALSE, results='asis'}
print(citation("cffr")[1], bibtex = FALSE)
```

A BibTeX entry for LaTeX users is

```{r echo=FALSE, comment=''}
toBibtex(citation("cffr")[1])
```

You can also use the [citation provided by
GitHub](https://github.com/ropensci/cffr), that is generated from the
information of a `CITATION.cff` created with **cffr**. See [About CITATION
files](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)
for more info.

## References

::: {#refs}
:::

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "\\BibTeX\\ and CFF"
subtitle: A potential crosswalk
bibliography: REFERENCES.bib
author: Diego Hernangómez
abstract: >-
  This article presents a crosswalk between \BibTeX\ and Citation File Format 
  [@druskat_citation_2021], as it is performed by the cffr package 
  [@hernangomez2021].
link-citations: yes
documentclass: article
urlcolor: blue
linkcolor: brown
header-includes:
  \usepackage{fvextra}
  \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
output:
  pdf_document:
    latex_engine: pdflatex
    keep_tex: true
    includes:
      in_header: preamble.tex
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "",
  tidy = "styler"
)

options(width = 60)
```
---
title: "BibTeX and CFF"
subtitle: A potential crosswalk
bibliography: REFERENCES.bib
author: Diego Hernangómez
description: >-
  This article presents a crosswalk between BibTeX and Citation File Format 
  [@druskat_citation_2021], as it is performed by the cffr package 
  [@hernangomez2021].
abstract: >-
  This article presents a crosswalk between BibTeX and Citation File Format 
  [@druskat_citation_2021], as it is performed by the cffr package 
  [@hernangomez2021]. Several crosswalk models specific for each BibTeX entry
  type [@patashnik1988] are proposed. The article also provide examples using real
  BibTeX entries and tips for developers that would like to implement the crosswalk
  on different programming languages.
link-citations: yes
documentclass: article
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{BibTeX and CFF}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)

options(width = 60)
```

## Citation

Please cite this article as:

```{r cit, echo=FALSE, results='asis'}
thisart <- bibentry("article",
  title = "{BibTeX} and {CFF}, a potential crosswalk",
  key = "hernangomez2022",
  author = "Diego Hernangómez",
  journal = "The {cffr} package",
  year = 2022,
  volume = "Vignettes",
)
cat("  \n")
thisart
```

A BibTeX entry for LaTeX users:

``` bibtex
@article{hernangomez2022,
    title        = {{BibTeX} and {CFF}, a potential crosswalk},
    author       = {Diego Hernangómez},
    year         = 2022,
    journal      = {The {cffr} package},
    volume       = {Vignettes}
}
```

## BibTeX and R

[BibTeX](https://en.wikipedia.org/wiki/BibTeX) is a well-known format for
storing references created by [Oren
Patashnik](https://en.wikipedia.org/wiki/Oren_Patashnik "Oren Patashnik") and
[Leslie Lamport](https://en.wikipedia.org/wiki/Leslie_Lamport "Leslie Lamport")
back in 1985. BibTeX that may be reused by another software, like
[LaTeX](https://en.wikipedia.org/wiki/LaTeX), for adding references to a work.
An example structure of a BibTeX entry would be:

``` bibtex
@book{einstein1921,
    title        = {Relativity: The Special and the General Theory},
    author       = {Einstein, A.},
    year         = 1920,
    publisher    = {Henry Holt and Company},
    address      = {London, United Kingdom},
    isbn         = 9781587340925
}
```

On this case, the entry (identified as `einstein1921`) would refer to a book.
This entry then can be used on a document and include references to it.

On **R** [@R_2021], we can replicate this structure using the `bibentry()` and
`toBibtex()` functions:

```{r bibentry, comment="#>"}

entry <- bibentry("book",
  key = "einstein1921",
  title = "Relativity: The Special and the General Theory",
  author = person("A.", "Einstein"),
  year = 1920,
  publisher = "Henry Holt and Company",
  address = "London, United Kingdom",
  isbn = 9781587340925,
)

toBibtex(entry)
```

The final results of the entry as a text string would be parsed as[^1]:

[^1]: By default R Pandoc would generate the cite on the Chicago author-date
    format [@rmarkdowncookbook2020]

```{r echo=FALSE, results='asis'}

entry
```

## BibTeX definitions

@patashnik1988 provides a comprehensive explanation of the BibTeX formats. We
can distinguish between **Entries** and **Fields**.

### Entries {#entries}

Each entry type defines a different type of work. The 14 entry types defined on
BibTeX[^2] are:

[^2]: Other implementations similar to BibTeX, as
    [BibLaTeX](https://www.ctan.org/pkg/biblatex), expand the definitions of
    entries including other types as **online**, **software** or **dataset**. On
    BibTeX these entries should be reclassified to **misc**.

-   **\@article**: An article from a journal or magazine.
-   **\@book**: A book with an explicit publisher.
-   **\@booklet**: A work that is printed and bound, but without a named
    publisher or sponsoring institution.
-   **\@conference**: The same as **\@inproceedings**, included for Scribe
    compatibility.
-   **\@inbook**: A part of a book, which may be a chapter (or section or
    whatever) and/or a range of pages.
-   **\@incollection**: A part of a book having its own title.
-   **\@inproceedings**: An article in a conference proceedings.
-   **\@manual**: Technical documentation.
-   **\@mastersthesis**: A Master's thesis.
-   **\@misc**: Use this type when nothing else fits.
-   **\@phdthesis**: A PhD thesis.
-   **\@proceedings**: The proceedings of a conference.
-   **\@techreport**: A report published by a school or other institution,
    usually numbered within a series.
-   **\@unpublished**: A document having an author and title, but not formally
    published.

Regarding the entries, `bibentry()` **R** function does not implement
**\@conference** . However, we can replace that key by **\@inproceedings** given
that the definition is identical.

### Fields

As in the case of Entries, @patashnik1988 provides also a definition for each of
the possible standard BibTeX fields[^3]. An entry can include other fields that
would be ignored on the raw implementation of BibTeX:

[^3]: As in the case of the entries, other implementations based on BibTeX may
    recognize additional fields.

-   **address**: Usually the address of the **publisher** or other of
    **institution**.
-   **annote**: An annotation. It is not used by the standard bibliography
    styles, but may be used by others that produce an annotated bibliography.
-   **author**: The name(s) of the author(s), in the format described in the
    LaTeX book [@lamport86latex].
-   **booktitle**: Title of a book, part of which is being cited. For **\@book**
    entries, use the **title** field instead.
-   **chapter**: A chapter (or section or whatever) number.
-   **crossref**: The database key of the entry being cross referenced.
-   **edition**: The edition of a **\@book** - for example, "Second". This
    should be an ordinal, and should have the first letter capitalized, the
    standard styles convert to lower case when necessary.
-   **editor**: Name(s) of editor(s), typed as indicated in the LaTeX book
    [@lamport86latex]. If there is also an **author** field, then the editor
    field gives the editor of the book or collection in which the reference
    appears.
-   **howpublished**: How something strange has been published. The first word
    should be capitalized.
-   **institution**: The sponsoring institution of a technical report.
-   **journal**: A journal name.
-   **key**: Used for alphabetizing, cross referencing, and creating a label
    when the **author** information is missing.
-   **month**: The month in which the work was published or, for an unpublished
    work, in which it was written. You should use the standard three-letter
    abbreviation, as described in Appendix B.1.3 of the LaTeX book
    [@lamport86latex] (i.e. `jan, feb, mar`).
-   **note**: Any additional information that can help the reader. The first
    word should be capitalized.
-   **number**: The number of a journal, magazine, technical report, or of a
    work in a series. An issue of a journal or magazine is usually identified by
    its **volume**: and number; the organization that issues a technical report
    usually gives it a number; and sometimes books are given numbers in a named
    series.
-   **organization**: The organization that sponsors a **\@conference** or that
    publishes a manual.
-   **pages**: One or more page numbers or range of numbers, such as `42--111`
    or `7,41,73--97` or `43+`.
-   **publisher**: The publisher's name.
-   **school**: The name of the school where a thesis was written.
-   **series**: The name of a series or set of books. When citing an entire
    book, the **title** field gives its title and an optional **series** field
    gives the name of a series or multi-volume set in which the book is
    published.
-   **title**: The work's title.
-   **type**: The type of a technical report---for example, "Research Note".
-   **volume**: The volume of a journal or multivolume book.
-   **year**: The year of publication or, for an unpublished work, the year it
    was written. Generally it should consist of four numerals, such as `1984`.

There is a strict relation between Entries and Fields on BibTeX. Depending on
the type of entries, some fields are required while others are optional or even
ignored. On the following table, required field are flagged as **\*** and
optional fields are flagged as **-**. Fields on parenthesis **()** denotes that
there are some degree of flexibility on the requirement of the field, see
@patashnik1988 for more information.

```{r entry_fields1, echo=FALSE}

bibtex_field_entry <- read.csv(system.file("extdata/bibtex_field_entry.csv",
  package = "cffr"
),
sep = ","
)

t1 <- bibtex_field_entry[, c(1:7)]

knitr::kable(t1,
  col.names = gsub("\\.", ",", names(t1)),
  align = c("l", rep("c", 6)),
  caption = "BibTeX, required fields by entry"
)
```

```{r entry_fields2, echo=FALSE}
t2 <- bibtex_field_entry[, c(1, 8:13)]

knitr::kable(t2,
  col.names = gsub("\\.", ",", names(t2)),
  align = c("l", rep("c", 6)),
  caption = "(cont) BibTeX, required fields by entry"
)
```

It can be observed that just a subset of fields is required in any of the
Entries. For example, **title**, **year** and **author** are either required or
optional on almost every entry, while **crossref**, **annote** or **key** are
never required.

## Citation File Format

[Citation File Format (CFF](https://citation-file-format.github.io/))
[@druskat_citation_2021] are plain text files with human- and machine-readable
citation information for software (and datasets). Among the [valid keys of
CFF](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#valid-keys)
there are two keys, `preferred-citation` and `references` of special interest
for citing and referring to related works[^4]:

[^4]: See [Guide to Citation File Format schema version
    1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#preferred-citation).

-   **`preferred-citation`**: A reference to another work that should be cited
    instead of the software or dataset itself.

-   **`references`**: Reference(s) to other creative works. Similar to a list of
    references in a paper, references of the software or dataset may include
    other software (dependencies), or other research products that the software
    or dataset builds on, but not work describing the software or dataset.

These two keys are expected to be
[`definition.reference`](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#definitionsreference)
objects, therefore they may contain the following keys:

```{r refkeys, echo=FALSE, message=FALSE, warning=FALSE, results='asis'}
library(cffr)

# Fill with whites
l <- c(cff_schema_definitions_refs(), rep("", 4))


refkeys <- matrix(l, ncol = 5, byrow = TRUE)

knitr::kable(refkeys,
  caption = "Valid keys on CFF `definition-reference` objects"
)
```

These keys are the equivalent to the fields of BibTeX (see [Fields]), with the
exception of the key **type**. On CFF, this key defines the type of work[^5],
therefore this is the equivalent to the BibTeX entries (see
[Entries](#entries)).

[^5]: See a complete list of possible values of CFF type on the [Guide to
    Citation File Format schema version
    1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#definitionsreferencetype).

## Proposed crosswalk

The **cffr** package [@hernangomez2021] provides utilities from converting
BibTeX entries (via the **R** base function `bibentry()`) to CFF files and
vice-versa. This section describes how the conversion between both formats have
been implemented. This crosswalk is based partially on
@Haines_Ruby_CFF_Library_2021[^6].

[^6]: Note that this software performs only the conversion from CFF to BibTeX,
    however **cffr** can perform the conversion in both directions.

On the following two section I present an overview of the proposed mapping
between the Entries and Fields of BibTeX and the CFF keys. After this initial
mapping, I propose further transformations to improve the compatibility between
both systems using different [Entry Models].

### Entry/Type crosswalk

For converting general BibTeX entries to CFF types, the following crosswalk is
proposed:

| BibTeX Entry        | Value of CFF key: type | Notes              |
|---------------------|------------------------|--------------------|
| **\@article**       | article                |                    |
| **\@book**          | book                   |                    |
| **\@booklet**       | pamphlet               |                    |
| **\@conference**    | conference-paper       |                    |
| **\@inbook**        | book                   | See [Entry Models] |
| **\@incollection**  | generic                | See [Entry Models] |
| **\@inproceedings** | conference-paper       |                    |
| **\@manual**        | manual                 |                    |
| **\@mastersthesis** | thesis                 | See [Entry Models] |
| **\@misc**          | generic                |                    |
| **\@phdthesis**     | thesis                 | See [Entry Models] |
| **\@proceedings**   | proceedings            |                    |
| **\@techreport**    | report                 |                    |
| **\@unpublished**   | unpublished            |                    |

: Entry/Type crosswalk: From BibTeX to CFF

Also, given that CFF provides with a [wide range of allowed
values](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#definitionsreferencetype)
on type, the following conversion would be performed from CFF to BibTeX:

+---------------------+--------------------------+----------------------------+
| Value of CFF key:   | BibTeX Entry             | Notes                      |
| type                |                          |                            |
+=====================+==========================+============================+
| book                | **\@book** **/           | See [Entry Models]         |
|                     | \@inbook**               |                            |
+---------------------+--------------------------+----------------------------+
| conference          | **\@inproceedings**      |                            |
+---------------------+--------------------------+----------------------------+
| conference-paper    | **\@inproceedings**      |                            |
+---------------------+--------------------------+----------------------------+
| magazine-article    | **\@article**            |                            |
+---------------------+--------------------------+----------------------------+
| manual              | **\@manual**             |                            |
+---------------------+--------------------------+----------------------------+
| newspaper-article   | **\@article**            |                            |
+---------------------+--------------------------+----------------------------+
| pamphlet            | **\@booklet**            |                            |
+---------------------+--------------------------+----------------------------+
| proceedings         | **\@proceedings**        |                            |
+---------------------+--------------------------+----------------------------+
| report              | **\@techreport**         |                            |
+---------------------+--------------------------+----------------------------+
| thesis              | **\@mastersthesis /      | See [Entry Models]         |
|                     | \@phdthesis**            |                            |
+---------------------+--------------------------+----------------------------+
| unpublished         | **\@unpublished**        |                            |
+---------------------+--------------------------+----------------------------+
| generic             | **\@misc /               | Under specific conditions, |
|                     | \@incollection**         | see [Entry Models]         |
+---------------------+--------------------------+----------------------------+
| \<any other value>  | **\@misc**               |                            |
+---------------------+--------------------------+----------------------------+

: Entry/Type crosswalk: From CFF to BibTeX

### Fields/Key crosswalk

There is a large degree of similarity between the definition and names of some
BibTeX fields and CFF keys[^7]. On the following cases, the equivalence is
almost straightforward:

[^7]: To avoid errors on the interpretation, names in **bold** correspond to
    BibTeX fields while the same number on [underscore]{.ul} would refer to CFF.

| BibTeX Field     | CFF key                   |
|------------------|---------------------------|
| **address**      | See [Entry Models]        |
| ***annote***     | \-                        |
| **author**       | [authors]{.ul}            |
| **booktitle**    | [collection-title]{.ul}   |
| **chapter**      | [section]{.ul}            |
| ***crossref***   | \-                        |
| **edition**      | [edition]{.ul}            |
| **editor**       | [editors]{.ul}            |
| **howpublished** | [medium]{.ul}             |
| **institution**  | See [Entry Models]        |
| **journal**      | [journal]{.ul}            |
| ***key***        | \-                        |
| **month**        | [month]{.ul}              |
| **note**         | [notes]{.ul}              |
| **number**       | [issue]{.ul}              |
| **organization** | See [Entry Models]        |
| **pages**        | [start]{.ul} & [end]{.ul} |
| **publisher**    | [publisher]{.ul}          |
| **school**       | See [Entry Models]        |
| **series**       | See [Entry Models]        |
| **title**        | [title]{.ul}              |
| ***type***       | \-                        |
| **volume**       | [volume]{.ul}             |
| **year**         | [year]{.ul}               |

: BibTeX - CFF Field/Key crosswalk

Additionally, there are other additional CFF keys that have a correspondence
with BibLaTeX fields. We propose also to include these fields on the
crosswalk[^8], although they are not part of the core BibTeX fields definition.

[^8]: See @biblatexcheatsheet for a preview of the accepted BibLaTeX fields.

| BibLaTeX Field | CFF key               |
|----------------|-----------------------|
| **abstract**   | [abstract]{.ul}       |
| **date**       | [date-published]{.ul} |
| **doi**        | [doi]{.ul}            |
| **file**       | [filename]{.ul}       |
| **isbn**       | [isbn]{.ul}           |
| **issn**       | [issn]{.ul}           |
| **issuetitle** | [issue-title]{.ul}    |
| **pagetotal**  | [pages]{.ul}          |
| **translator** | [translators]{.ul}    |
| **url**        | [url]{.ul}            |
| **urldate**    | [date-accessed]{.ul}  |
| **version**    | [version]{.ul}        |

: BibLaTeX - CFF Field/Key crosswalk

## Entry Models

This section presents the specific mapping proposed for each of the BibTeX
entries, providing further information on how each field is treated. Examples
are adapted from the [xampl.bib](https://tug.org/texmf-docs/bibtex/xampl.bib)
file provided with the bibtex package [@patashnik].

### article

The crosswalk of **\@article** does not require any special treatment.

+------------------+----------------------+-----------------------------------+
| BibTeX           | CFF                  | Note                              |
+==================+======================+===================================+
| **\@article**    | [type: article]{.ul} | When converting CFF to BibTeX,    |
|                  |                      | [type magazine-article]{.ul} and  |
|                  |                      | [newspaper-article]{.ul} are      |
|                  |                      | converted to **\@article**.       |
+------------------+----------------------+-----------------------------------+
| **author\***     | [authors]{.ul}       |                                   |
+------------------+----------------------+-----------------------------------+
| **title\***      | [title]{.ul}         |                                   |
+------------------+----------------------+-----------------------------------+
| **journal\***    | [journal]{.ul}       |                                   |
+------------------+----------------------+-----------------------------------+
| **year\***       | [year]{.ul}          |                                   |
+------------------+----------------------+-----------------------------------+
| **volume**       | [volume]{.ul}        |                                   |
+------------------+----------------------+-----------------------------------+
| **number**       | [issue]{.ul}         |                                   |
+------------------+----------------------+-----------------------------------+
| **pages**        | [start]{.ul} and     | Separated by `--`, i.e,           |
|                  | [end]{.ul}           |                                   |
|                  |                      | **pages** = {3--5}                |
|                  |                      |                                   |
|                  |                      | would be parsed as                |
|                  |                      |                                   |
|                  |                      | [start]{.ul}: 3                   |
|                  |                      |                                   |
|                  |                      | [end]{.ul}: 5                     |
+------------------+----------------------+-----------------------------------+
| **month**        | [month]{.ul}         | As a fallback, **month** could be |
|                  |                      | extracted also from **date**      |
|                  |                      | (BibLaTeX field)/                 |
|                  |                      | [date-published]{.ul}             |
+------------------+----------------------+-----------------------------------+
| **note**         | [notes]{.ul}         |                                   |
+------------------+----------------------+-----------------------------------+

: **\@article** Model

**Examples**

[*BibTeX entry*]{.ul}

``` bibtex
@article{article-full,
    title        = {The Gnats and Gnus Document Preparation System},
    author       = {Leslie A. Aamport},
    year         = 1986,
    month        = jul,
    journal      = {{G-Animal's} Journal},
    volume       = 41,
    number       = 7,
    pages        = {73+},
    note         = {This is a full ARTICLE entry}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE}

bib <- bibentry("article",
  title        = "The Gnats and Gnus Document Preparation System",
  author       = "Leslie A. Aamport",
  year         = 1986,
  month        = "jul",
  journal      = "{G-Animal's} Journal",
  volume       = 41,
  number       = 7,
  pages        = "73+",
  note         = "This is a full ARTICLE entry"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### book/inbook

In terms of field required on BibTeX, the only difference between **\@book** and
**\@inbook** is that the latter requires also a **chapter** or **pages**, while
for **\@book** these fields are not even optional. So we propose here to
identify an **\@inbook** on CFF as a [book]{.ul} with [section]{.ul} and
[start]{.ul}-[end]{.ul} fields (CFF).

Another specificity is that **series** field is mapped to
[collection-title]{.ul} and **address** is mapped as the [address]{.ul} of the
[publisher]{.ul} (CFF).

+------------------------+--------------------------+--------------------------+
| BibTeX                 | CFF                      | Note                     |
+========================+==========================+==========================+
| **\@book**             | [type: book]{.ul}        |                          |
+------------------------+--------------------------+--------------------------+
| **\@inbook**           | [type: book]{.ul}        | For identifying an       |
|                        |                          | **\@inbook** in CFF,     |
|                        |                          | assess if [section]{.ul} |
|                        |                          | or                       |
|                        |                          | [start]{.ul}-[end]{.ul}  |
|                        |                          | information is available |
+------------------------+--------------------------+--------------------------+
| **author\***           | [authors]{.ul}           |                          |
+------------------------+--------------------------+--------------------------+
| **editor\***           | [editors]{.ul}           |                          |
+------------------------+--------------------------+--------------------------+
| **title\***            | [title]{.ul}             |                          |
+------------------------+--------------------------+--------------------------+
| **publisher\***        | [publisher]{.ul}         |                          |
+------------------------+--------------------------+--------------------------+
| **year\***             | [year]{.ul}              |                          |
+------------------------+--------------------------+--------------------------+
| **chapter\***          | [section]{.ul}           | Only required on         |
|                        |                          | **\@inbook**             |
+------------------------+--------------------------+--------------------------+
| **pages\***            | [start]{.ul} and         | Only required on         |
|                        | [end]{.ul}               | **\@inbook**             |
+------------------------+--------------------------+--------------------------+
| **volume**             | [volume]{.ul}            |                          |
+------------------------+--------------------------+--------------------------+
| **number**             | [issue]{.ul}             |                          |
+------------------------+--------------------------+--------------------------+
| **series**             | [collection-title]{.ul}  |                          |
+------------------------+--------------------------+--------------------------+
| **address**            | [address]{.ul} property  | As a fallback, the field |
|                        | of [publisher]{.ul}      | [location]{.ul} can be   |
|                        |                          | used                     |
+------------------------+--------------------------+--------------------------+
| **edition**            | [edition]{.ul}           |                          |
+------------------------+--------------------------+--------------------------+
| **month**              | [month]{.ul}             | See **Note** on          |
|                        |                          | [article]                |
+------------------------+--------------------------+--------------------------+
| **note**               | [notes]{.ul}             |                          |
+------------------------+--------------------------+--------------------------+

: **\@book/\@inbook** Model

**Examples: book**

[*BibTeX entry*]{.ul}

``` bibtex
@book{book-full,
    title        = {Seminumerical Algorithms},
    author       = {Donald E. Knuth},
    year         = 1981,
    month        = 10,
    publisher    = {Addison-Wesley},
    address      = {Reading, Massachusetts},
    series       = {The Art of Computer Programming},
    volume       = 2,
    note         = {This is a full BOOK entry},
    edition      = {Second}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE}

bib <- bibentry("Book",
  title = "Seminumerical Algorithms",
  author = "Donald E. Knuth",
  year = 1981,
  month = 10,
  publisher = "Addison-Wesley",
  address = "Reading, Massachusetts",
  series = "The Art of Computer Programming",
  volume = 2,
  note = "This is a full BOOK entry",
  edition = "Second"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r,  echo=FALSE}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

**Examples: inbook**

[*BibTeX entry*]{.ul}

``` bibtex
@inbook{inbook-full,
    title        = {Fundamental Algorithms},
    author       = {Donald E. Knuth},
    year         = 1973,
    month        = 10,
    publisher    = {Addison-Wesley},
    address      = {Reading, Massachusetts},
    series       = {The Art of Computer Programming},
    volume       = 1,
    pages        = {10--119},
    note         = {This is a full INBOOK entry},
    edition      = {Second},
    type         = {Section},
    chapter      = {1.2}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE,}

bib <- bibentry("inbook",
  title        = "Fundamental Algorithms",
  author       = "Donald E. Knuth",
  year         = 1973,
  month        = 10,
  publisher    = "Addison-Wesley",
  address      = "Reading, Massachusetts",
  series       = "The Art of Computer Programming",
  volume       = 1,
  pages        = "10--119",
  note         = "This is a full INBOOK entry",
  edition      = "Second",
  type         = "Section",
  chapter      = "1.2"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE,}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### booklet

In **\@booklet** **address** is mapped to [location]{.ul}.

+-----------------------------+------------------+----------------------------+
| BibTeX                      | CFF              | Note                       |
+=============================+==================+============================+
| **\@booklet**               | [type:           |                            |
|                             | pamphlet]{.ul}   |                            |
+-----------------------------+------------------+----------------------------+
| **title\***                 | [title]{.ul}     |                            |
+-----------------------------+------------------+----------------------------+
| **author\***                | [authors]{.ul}   |                            |
+-----------------------------+------------------+----------------------------+
| **howpublished**            | [medium]{.ul}    |                            |
+-----------------------------+------------------+----------------------------+
| **address**                 | [location]{.ul}  |                            |
+-----------------------------+------------------+----------------------------+
| **month**                   | [month]{.ul}     | See **Note** on [article]  |
+-----------------------------+------------------+----------------------------+
| **year**                    | [year]{.ul}      | As a fallback, **year**    |
|                             |                  | could be extracted also    |
|                             |                  | from **date** (BibLaTeX    |
|                             |                  | field)/                    |
|                             |                  | [date-published]{.ul}      |
+-----------------------------+------------------+----------------------------+
| **note**                    | [notes]{.ul}     |                            |
+-----------------------------+------------------+----------------------------+

: **\@booklet** Model

**Examples**

[*BibTeX entry*]{.ul}

``` bibtex
@booklet{booklet-full,
    title        = {The Programming of Computer Art},
    author       = {Jill C. Knvth},
    date         = {1988-03-14},
    month        = feb,
    address      = {Stanford, California},
    note         = {This is a full BOOKLET entry},
    howpublished = {Vernier Art Center}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE, }

bib <- bibentry("booklet",
  title = "The Programming of Computer Art",
  author = "Jill C. Knvth",
  date = "1988-03-14",
  month = "feb",
  address = "Stanford, California",
  note = "This is a full BOOKLET entry",
  howpublished = "Vernier Art Center"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE, }
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### conference/inproceedings

Note that in this case, **organization** is mapped to [institution]{.ul}, as
BibTeX does not prescribe the use of **institution** on these entries.

+-------------------------+---------------------------+-----------------------+
| BibTeX                  | CFF                       | Note                  |
+=========================+===========================+=======================+
| **\@conference /        | [type:                    | CFF entries with      |
| \@inproceedings**       | conference-paper]{.ul}    | [type]{.ul}           |
|                         |                           | =conference are       |
|                         |                           | mapped back to        |
|                         |                           | **\@inproceedings**.  |
+-------------------------+---------------------------+-----------------------+
| **author\***            | [authors]{.ul}            |                       |
+-------------------------+---------------------------+-----------------------+
| **title\***             | [title]{.ul}              |                       |
+-------------------------+---------------------------+-----------------------+
| **booktitle\***         | [collection-title]{.ul}   |                       |
+-------------------------+---------------------------+-----------------------+
| **year\***              | [year]{.ul}               |                       |
+-------------------------+---------------------------+-----------------------+
| **editor**              | [editors]{.ul}            |                       |
+-------------------------+---------------------------+-----------------------+
| **volume**              | [volume]{.ul}             |                       |
+-------------------------+---------------------------+-----------------------+
| **number**              | [issue]{.ul}              |                       |
+-------------------------+---------------------------+-----------------------+
| **series**              | [conference]{.ul}         |                       |
+-------------------------+---------------------------+-----------------------+
| **pages**               | [start]{.ul} and          | See **Note** on       |
|                         | [end]{.ul}                | [article]             |
+-------------------------+---------------------------+-----------------------+
| **address**             | [location]{.ul}           | As a fallback,        |
|                         |                           | [address]{.ul}        |
|                         |                           | property of           |
|                         |                           | [conference]{.ul} can |
|                         |                           | be used               |
+-------------------------+---------------------------+-----------------------+
| **month**               | [month]{.ul}              | See **Note** on       |
|                         |                           | [article]             |
+-------------------------+---------------------------+-----------------------+
| **organization**        | [institution]{.ul}        |                       |
+-------------------------+---------------------------+-----------------------+
| **publisher**           | [publisher]{.ul}          |                       |
+-------------------------+---------------------------+-----------------------+
| **note**                | [notes]{.ul}              |                       |
+-------------------------+---------------------------+-----------------------+

: **\@conference/\@inproceedings** Model

**Examples**

[*BibTeX entry*]{.ul}

``` bibtex
@inproceedings{inproceedings-full,
    title        = {On Notions of Information Transfer in {VLSI} Circuits},
    author       = {Alfred V. Oaho and Jeffrey D. Ullman and Mihalis Yannakakis},
    year         = 1983,
    month        = mar,
    booktitle    = {Proc. Fifteenth Annual ACM Symposium on the Theory of Computing},
    publisher    = {Academic Press},
    address      = {Boston},
    series       = {All ACM Conferences},
    number       = 17,
    pages        = {133--139},
    editor       = {Wizard V. Oz and Mihalis Yannakakis},
    organization = {The OX Association for Computing Machinery}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE,}

bib <- bibentry("inproceedings",
  title        = "On Notions of Information Transfer in {VLSI} Circuits",
  author       = "Alfred V. Oaho and Jeffrey D. Ullman and Mihalis Yannakakis",
  year         = 1983,
  month        = "mar",
  booktitle    = "Proc. Fifteenth Annual ACM Symposium on the Theory of Computing",
  publisher    = "Academic Press",
  address      = "Boston",
  series       = "All ACM Conferences",
  number       = 17,
  pages        = "133--139",
  editor       = "Wizard V. Oz and Mihalis Yannakakis",
  organization = "The OX Association for Computing Machinery"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE,}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### incollection

As **booktitle** is a required field, we propose to map that field to
[collection-title]{.ul} and the [type]{.ul} to [generic]{.ul}. Therefore, an
**\@incollection** is a [type: generic]{.ul} with a [collection-title]{.ul} key.

+-------------------------+-------------------------+-------------------------+
| **BibTeX**              | CFF                     | Note                    |
+=========================+=========================+=========================+
| **\@incollection**      | [type: generic]{.ul}    | Including a             |
|                         |                         | [collection-title]{.ul} |
|                         |                         | value                   |
+-------------------------+-------------------------+-------------------------+
| **author\***            | [authors]{.ul}          |                         |
+-------------------------+-------------------------+-------------------------+
| **title\***             | [title]{.ul}            |                         |
+-------------------------+-------------------------+-------------------------+
| **booktitle\***         | [collection-title]{.ul} |                         |
+-------------------------+-------------------------+-------------------------+
| **publisher\***         | [publisher]{.ul}        |                         |
+-------------------------+-------------------------+-------------------------+
| **year\***              | [year]{.ul}             |                         |
+-------------------------+-------------------------+-------------------------+
| **editor**              | [editors]{.ul}          |                         |
+-------------------------+-------------------------+-------------------------+
| **volume**              | [volume]{.ul}           |                         |
+-------------------------+-------------------------+-------------------------+
| **number**              | [issue]{.ul}            |                         |
+-------------------------+-------------------------+-------------------------+
| **series**              | [series]{.ul}           |                         |
+-------------------------+-------------------------+-------------------------+
| **type**                | \-                      | Ignored                 |
+-------------------------+-------------------------+-------------------------+
| **chapter**             | [section]{.ul}          |                         |
+-------------------------+-------------------------+-------------------------+
| **pages**               | [start]{.ul} and        | See **Note** on         |
|                         | [end]{.ul}              | [article]               |
+-------------------------+-------------------------+-------------------------+
| **address**             | [address]{.ul} property | See **Note** on         |
|                         | of [publisher]{.ul}     | [book/inbook]           |
+-------------------------+-------------------------+-------------------------+
| **edition**             | [edition]{.ul}          |                         |
+-------------------------+-------------------------+-------------------------+
| **month**               | [month]{.ul}            | See **Note** on         |
|                         |                         | [article]               |
+-------------------------+-------------------------+-------------------------+
| **note**                | [notes]{.ul}            |                         |
+-------------------------+-------------------------+-------------------------+

: **\@incollection** Model

**Examples**

[*BibTeX entry*]{.ul}

``` bibtex
@incollection{incollection-full,
    title        = {Semigroups of Recurrences},
    author       = {Daniel D. Lincoll},
    year         = 1977,
    month        = sep,
    booktitle    = {High Speed Computer and Algorithm Organization},
    publisher    = {Academic Press},
    address      = {New York},
    series       = {Fast Computers},
    number       = 23,
    pages        = {179--183},
    note         = {This is a full INCOLLECTION entry},
    editor       = {David J. Lipcoll and D. H. Lawrie and A. H. Sameh},
    chapter      = 3,
    type         = {Part},
    edition      = {Third}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE,}

bib <- bibentry("incollection",
  title        = "Semigroups of Recurrences",
  author       = "Daniel D. Lincoll",
  year         = 1977,
  month        = "sep",
  booktitle    = "High Speed Computer and Algorithm Organization",
  publisher    = "Academic Press",
  address      = "New York",
  series       = "Fast Computers",
  number       = 23,
  pages        = "179--183",
  note         = "This is a full INCOLLECTION entry",
  editor       = "David J. Lipcoll and D. H. Lawrie and A. H. Sameh",
  chapter      = 3,
  type         = "Part",
  edition      = "Third"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE,}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### manual

As in the case of [conference/inproceedings], **organization** is mapped to
[institution]{.ul}.

+--------------------------+-----------------------------+---------------------+
| **BibTeX**               | CFF                         | Note                |
+==========================+=============================+=====================+
| **\@manual**             | [type: manual]{.ul}         |                     |
+--------------------------+-----------------------------+---------------------+
| **title\***              | [title]{.ul}                |                     |
+--------------------------+-----------------------------+---------------------+
| **author**               | [authors]{.ul}              |                     |
+--------------------------+-----------------------------+---------------------+
| **organization**         | [institution]{.ul}          |                     |
+--------------------------+-----------------------------+---------------------+
| **address**              | [location]{.ul}             |                     |
+--------------------------+-----------------------------+---------------------+
| **edition**              | [edition]{.ul}              |                     |
+--------------------------+-----------------------------+---------------------+
| **month**                | [month]{.ul}                | See **Note** on     |
|                          |                             | [article]           |
+--------------------------+-----------------------------+---------------------+
| **year**                 | [year]{.ul}                 | See **Note** on     |
|                          |                             | [booklet]           |
+--------------------------+-----------------------------+---------------------+
| **note**                 | [notes]{.ul}                |                     |
+--------------------------+-----------------------------+---------------------+

: **\@manual** Model

**Examples**

[*BibTeX entry*]{.ul}

Note that **month** can't be parsed to a single integer in the range `1--12` as
required on CFF, so it is not parsed to avoid validation errors.

``` bibtex
@manual{manual-full,
  title        = {The Definitive Computer Manual},
    author       = {Larry Manmaker},
    year         = 1986,
    month        = {apr-may},
    address      = {Silicon Valley},
    note         = {This is a full MANUAL entry},
    organization = {Chips-R-Us},
    edition      = {Silver}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE,}

bib <- bibentry("Manual",
  title        = "The Definitive Computer Manual",
  author       = "Larry Manmaker",
  year         = 1986,
  month        = "apr-may",
  address      = "Silicon Valley",
  note         = "This is a full MANUAL entry",
  organization = "Chips-R-Us",
  edition      = "Silver"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE,}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### mastersthesis/phdthesis

In terms of field required on BibTeX, it is identical for both
**\@mastersthesis** and **\@phdthesis.**

We propose here to identify each type of thesis using the field
[thesis-type]{.ul} (CFF). So if [thesis-type]{.ul} contains a [regex
pattern](https://regex101.com/r/mBWfbs/1) `(?i)(phd)` it would be recognized as
**\@phdthesis**.

+---------------------------+----------------------+---------------------------+
| BibTeX                    | CFF                  | Note                      |
+===========================+======================+===========================+
| **\@mastersthesis**       | [type: thesis]{.ul}  | Use also                  |
|                           |                      | [thesis-type]{.ul} for    |
|                           |                      | identifying the thesis    |
|                           |                      | type.                     |
+---------------------------+----------------------+---------------------------+
| **\@phdthesis**           | [type: thesis]{.ul}  | Use also                  |
|                           |                      | [thesis-type]{.ul} for    |
|                           |                      | identifying the thesis    |
|                           |                      | type.                     |
+---------------------------+----------------------+---------------------------+
| **author\***              | [authors]{.ul}       |                           |
+---------------------------+----------------------+---------------------------+
| **title\***               | [title]{.ul}         |                           |
+---------------------------+----------------------+---------------------------+
| **school\***              | [institution]{.ul}   |                           |
+---------------------------+----------------------+---------------------------+
| **year\***                | [year]{.ul}          |                           |
+---------------------------+----------------------+---------------------------+
| **type**                  |                      |                           |
+---------------------------+----------------------+---------------------------+
| **address**               | [address]{.ul}       | See **Note** on           |
|                           | property of          | [book/inbook]             |
|                           | [institution]{.ul}   |                           |
+---------------------------+----------------------+---------------------------+
| **month**                 | [month]{.ul}         | See **Note** on [article] |
+---------------------------+----------------------+---------------------------+
| **note**                  | [notes]{.ul}         |                           |
+---------------------------+----------------------+---------------------------+

: **\@mastersthesis/phdthesis** Model

**Examples: mastersthesis**

[*BibTeX entry*]{.ul}

``` bibtex
@mastersthesis{mastersthesis-full,
    title        = {Mastering Thesis Writing},
    author       = {Edouard Masterly},
    year         = 1988,
    month        = jun,
    address      = {English Department},
    note         = {This is a full MASTERSTHESIS entry},
    school       = {Stanford University},
    type         = {Master's project}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE}

bib <- bibentry("mastersthesis",
  title        = "Mastering Thesis Writing",
  author       = "Edouard Masterly",
  year         = 1988,
  month        = "jun",
  address      = "English Department",
  note         = "This is a full MASTERSTHESIS entry",
  school       = "Stanford University",
  type         = "Master's project"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r,  echo=FALSE}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

**Examples: phdthesis**

[*BibTeX entry*]{.ul}

``` bibtex
@phdthesis{phdthesis-full,
    title        = {Fighting Fire with Fire: Festooning {F}rench Phrases},
    author       = {F. Phidias Phony-Baloney},
    year         = 1988,
    month        = jun,
    address      = {Department of French},
    note         = {This is a full PHDTHESIS entry},
    school       = {Fanstord University},
    type         = {{PhD} Dissertation}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE,}

bib <- bibentry("phdthesis",
  title = "Fighting Fire with Fire: Festooning {F}rench Phrases",
  author = "F. Phidias Phony-Baloney",
  year = 1988,
  month = "jun",
  address = "Department of French",
  note = "This is a full PHDTHESIS entry",
  school = "Fanstord University",
  type = "{PhD} Dissertation"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE,}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### misc

The crosswalk of **\@misc** does not require any special treatment. This entry
does not require any field.

Note als that it is mapped to [type: generic]{.ul} as [incollection], but in
this case **booktitle** is not even an option, so the proposed definition should
cover both **\@misc** and **\@incollection** without problems.

+--------------------------+-----------------------------+---------------------+
| **BibTeX**               | CFF                         | Note                |
+==========================+=============================+=====================+
| **\@misc**               | [type: generic]{.ul}        |                     |
+--------------------------+-----------------------------+---------------------+
| **author**               | [authors]{.ul}              |                     |
+--------------------------+-----------------------------+---------------------+
| **title**                | [title]{.ul}                |                     |
+--------------------------+-----------------------------+---------------------+
| **howpublished**         | [medium]{.ul}               |                     |
+--------------------------+-----------------------------+---------------------+
| **month**                | [month]{.ul}                | See **Note** on     |
|                          |                             | [article]           |
+--------------------------+-----------------------------+---------------------+
| **year**                 | [year]{.ul}                 | See **Note** on     |
|                          |                             | [booklet]           |
+--------------------------+-----------------------------+---------------------+
| **note**                 | [notes]{.ul}                |                     |
+--------------------------+-----------------------------+---------------------+

: **\@misc** Model

**Examples**

[*BibTeX entry*]{.ul}

``` bibtex
@misc{misc-full,
    title        = {Handing out random pamphlets in airports},
    author       = {Joe-Bob Missilany},
    year         = 1984,
    month        = oct,
    note         = {This is a full MISC entry},
    howpublished = {Handed out at O'Hare}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE,}

bib <- bibentry("Misc",
  title        = "Handing out random pamphlets in airports",
  year         = 1984,
  month        = "oct",
  note         = "This is a MISC entry",
  howpublished = "Handed out at O'Hare"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE,}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### proceedings

The proposed model is similar to [conference/inproceedings]. Note that
**\@proceedings** does not prescribe a **author** field. On this cases, as
[authors]{.ul} is required on CFF, we would use *anonymous*[^9] when converting
to CFF and omit it on the conversion back to CFF.

[^9]: As proposed on [*How to deal with unknown individual
    authors?*](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#how-to-deal-with-unknown-individual-authors),
    **(Guide to Citation File Format schema version 1.2.0)**

+-------------------------+---------------------------+-----------------------+
| BibTeX                  | CFF                       | Note                  |
+=========================+===========================+=======================+
| **\@proceedings**       | [type: proceedings]{.ul}  |                       |
+-------------------------+---------------------------+-----------------------+
| **title\***             | [title]{.ul}              |                       |
+-------------------------+---------------------------+-----------------------+
| **year\***              | [year]{.ul}               |                       |
+-------------------------+---------------------------+-----------------------+
| **editor**              | [editors]{.ul}            |                       |
+-------------------------+---------------------------+-----------------------+
| **volume**              | [volume]{.ul}             |                       |
+-------------------------+---------------------------+-----------------------+
| **number**              | [issue]{.ul}              |                       |
+-------------------------+---------------------------+-----------------------+
| **series**              | [conference]{.ul}         |                       |
+-------------------------+---------------------------+-----------------------+
| **address**             | [location]{.ul}           | As a fallback,        |
|                         |                           | [address]{.ul}        |
|                         |                           | property of           |
|                         |                           | [conference]{.ul} can |
|                         |                           | be used               |
+-------------------------+---------------------------+-----------------------+
| **month**               | [month]{.ul}              | See **Note** on       |
|                         |                           | [article]             |
+-------------------------+---------------------------+-----------------------+
| **organization**        | [institution]{.ul}        |                       |
+-------------------------+---------------------------+-----------------------+
| **publisher**           | [publisher]{.ul}          |                       |
+-------------------------+---------------------------+-----------------------+
| **note**                | [notes]{.ul}              |                       |
+-------------------------+---------------------------+-----------------------+

: **\@conference/\@inproceedings** Model

**Examples**

[*BibTeX entry*]{.ul}

``` bibtex
@proceedings{proceedings-full,
    title        = {Proc. Fifteenth Annual ACM Symposium on the Theory of Computing},
    year         = 1983,
    month        = mar,
    publisher    = {Academic Press},
    address      = {Boston},
    series       = {All ACM Conferences},
    number       = 17,
    note         = {This is a full PROCEEDINGS entry},
    editor       = {Wizard V. Oz and Mihalis Yannakakis},
    organization = {The OX Association for Computing Machinery}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE,}

bib <- bibentry("proceedings",
  title        = "Proc. Fifteenth Annual ACM Symposium on the Theory of Computing",
  year         = 1983,
  month        = "mar",
  publisher    = "Academic Press",
  address      = "Boston",
  series       = "All ACM Conferences",
  number       = 17,
  note         = "This is a full PROCEEDINGS entry",
  editor       = "Wizard V. Oz and Mihalis Yannakakis",
  organization = "The OX Association for Computing Machinery"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE,}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### techreport

+--------------------------+-----------------------------+---------------------+
| **BibTeX**               | CFF                         | Note                |
+==========================+=============================+=====================+
| **\@techreport**         | [type: report]{.ul}         |                     |
+--------------------------+-----------------------------+---------------------+
| **author\***             | [authors]{.ul}              |                     |
+--------------------------+-----------------------------+---------------------+
| **title\***              | [title]{.ul}                |                     |
+--------------------------+-----------------------------+---------------------+
| **institution\***        | [institution]{.ul}          |                     |
+--------------------------+-----------------------------+---------------------+
| **year\***               | [year]{.ul}                 |                     |
+--------------------------+-----------------------------+---------------------+
| **type**                 | \-                          | Ignored             |
+--------------------------+-----------------------------+---------------------+
| **number**               | [issue]{.ul}                |                     |
+--------------------------+-----------------------------+---------------------+
| **address**              | [address]{.ul} property of  | See **Note** on     |
|                          | [institution]{.ul}          | [book/inbook]       |
+--------------------------+-----------------------------+---------------------+
| **month**                | [month]{.ul}                | See **Note** on     |
|                          |                             | [article]           |
+--------------------------+-----------------------------+---------------------+
| **note**                 | [notes]{.ul}                |                     |
+--------------------------+-----------------------------+---------------------+

: **\@techreport** Model

**Examples**

[*BibTeX entry*]{.ul}

``` bibtex
@techreport{techreport-full,
    title        = {A Sorting Algorithm},
    author       = {Tom Terrific},
    year         = 1988,
    month        = oct,
    address      = {Computer Science Department, Fanstord, California},
    number       = 7,
    note         = {This is a full TECHREPORT entry},
    institution  = {Fanstord University},
    type         = {Wishful Research Result}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE,}

bib <- bibentry("techreport",
  title = "A Sorting Algorithm",
  author = "Tom Terrific",
  year = 1988,
  month = "oct",
  address = "Computer Science Department, Fanstord, California",
  number = 7,
  note = "This is a full TECHREPORT entry",
  institution = "Fanstord University",
  type = "Wishful Research Result"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE,}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

### unpublished

+--------------------------+-----------------------------+---------------------+
| **BibTeX**               | CFF                         | Note                |
+==========================+=============================+=====================+
| **\@unpublished**        | [type: unpublished]{.ul}    |                     |
+--------------------------+-----------------------------+---------------------+
| **author\***             | [authors]{.ul}              |                     |
+--------------------------+-----------------------------+---------------------+
| **title\***              | [title]{.ul}                |                     |
+--------------------------+-----------------------------+---------------------+
| **note\***               | [notes]{.ul}                |                     |
+--------------------------+-----------------------------+---------------------+
| **month**                | [month]{.ul}                | See **Note** on     |
|                          |                             | [article]           |
+--------------------------+-----------------------------+---------------------+
| **year**                 | [year]{.ul}                 | See **Note** on     |
|                          |                             | [booklet]           |
+--------------------------+-----------------------------+---------------------+

: **\@unpublished** Model

**Examples**

[*BibTeX entry*]{.ul}

``` bibtex
@unpublished{unpublished-minimal,
    title        = {Lower Bounds for Wishful Research Results},
    author       = {Ulrich Underwood and Ned Net and Paul Pot},
    note         = {Talk at Fanstord University (this is a minimal UNPUBLISHED entry)}
}
```

[*CFF entry*]{.ul}

```{r echo=FALSE,}

bib <- bibentry("unpublished",
  title        = "Lower Bounds for Wishful Research Results",
  author       = "Ulrich Underwood and Ned Net and Paul Pot",
  note         = "Talk at Fanstord University (this is a minimal UNPUBLISHED entry)"
)

cff_parse_citation(bib)
```

[*From CFF to BibTeX*]{.ul}

```{r echo=FALSE,}
toBibtex(cff_to_bibtex(cff_parse_citation(bib)))
```

## References
---
title: "cffr: Generate Citation File Format Metadata for R Packages"
subtitle: "JOSS paper"
description: >
  Paper published on The Journal of Open Source Software.
tags:
  - R
  - cff
  - citation
  - credit
  - metadata
author: Diego Hernangómez
date: 09 November 2021
output: rmarkdown::html_vignette
bibliography: REFERENCES.bib
link-citations: yes
vignette: >
  %\VignetteIndexEntry{cffr: Generate Citation File Format Metadata for R Packages}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03900/status.svg)](https://doi.org/10.21105/joss.03900)

## Summary

The Citation File Format project [@druskat_citation_2021] defines a standardized
format for providing software or datasets citation metadata in plaintext files
that are easy to read by both humans and machines.

This metadata format is being adopted by GitHub as the primary format for its
built-in citation support [@github_about_citation]. Other leading archives for
scientific software, including Zenodo and Zotero [@druskat_stephan_making_2021],
have included as well support for CITATION.cff files in their GitHub
integration.

The cffr package provides utilities to generate and validate these CITATION.cff
files automatically for R [@R_2021] packages by parsing the DESCRIPTION file and
the native R citation file. The package also includes utilities and examples for
parsing components as persons and additional citations, as well as several
vignettes which illustrate both the basic usage of the package as well as some
more technical details about the metadata extraction process.

## Statement of need

Citation of research software on research project is often omitted [@salmon2021]
. Among many reasons why software is not cited, one is the lack of a clear
citation information from package developers.

Some of the main reasons for citing software used on research are:

1.  **Reproducibility**: Software and their versions are important information
    to include in any research project. It helps peers to understand and
    reproduce effectively the results of any work. Including versions is also
    crucial as a way of recording the context of your manuscript when software
    changes.
2.  **Developer Credit:** On the context of Free and Open Source Software
    (FOSS), many of the software developers themselves are also researches.
    Receive credit for software development shouldn't be different from the
    credit received on other formats, as books or articles.

CITATION.cff files provides a clear citation rules for software. The format is
easily readable by humans and also can be parsed by appropriate software. The
adoption of GitHub of this format sends a strong message that research software
is something worthy of citation, and therefore deserves credit.

The cffr package allow R software developers to create CITATION.cff files from
the metadata already included on the package. Additionally, the package also
include validation tools via the jsonvalidate package [@jsonvalidate2021], that
allow developers to assess the validity of the file created using the latest CFF
schema.json.

## Acknowledgements

I would like to thank [Carl
Boettiger](https://ropensci.org/author/carl-boettiger/), [Maëlle
Salmon](https://ropensci.org/author/ma%C3%ABlle-salmon/) and the rest of
contributors of the [codemetar](https://docs.ropensci.org/codemetar/) package.
This package was the primary inspiration for developing cffr and shares a common
goal of increasing awareness on the efforts of software developers.

I would like also to thank [João Martins](https://zambujo.github.io/) and [Scott
Chamberlain](https://ropensci.org//author/scott-chamberlain/) for thorough
reviews, that helps improving the package and the documentation as well as
[Emily Riederer](https://emilyriederer.netlify.app/) for handling the [review
process](https://github.com/ropensci/software-review/issues/463).

## Citation

Hernangómez, D., (2021). cffr: Generate Citation File Format Metadata for R
Packages. Journal of Open Source Software, 6(67), 3900,
<https://doi.org/10.21105/joss.03900>

``` bibtex
@article{hernangomez2021,
  doi = {10.21105/joss.03900},
  url = {https://doi.org/10.21105/joss.03900},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {67},
  pages = {3900},
  author = {Diego Hernangómez},
  title = {cffr: Generate Citation File Format Metadata for R Packages},
  journal = {Journal of Open Source Software}
}
```

## References
---
title: "From R to CFF"
subtitle: "Crosswalk"
description: >
  A comprehenshive description of the internal mappings performed by `cffr`.
author: Diego Hernangómez
bibliography: REFERENCES.bib
link-citations: yes
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{From R to CFF}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(cffr)
```

The goal of this vignette is to provide an explicit map between the metadata
fields used by **cffr** and each one of the valid keys of the [Citation File
Format schema version
1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#valid-keys).

## Summary {#summary}

We summarize here the fields that **cffr** can parse and the original source of
information for each one of them. The details on each key are presented on the
next section of the document. The assessment of fields are based on the [Guide
to Citation File Format schema version
1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#valid-keys)
[@druskat_citation_2021].

```{r summary , echo=FALSE}

keys <- cff_schema_keys(sorted = TRUE)
origin <- vector(length = length(keys))
origin[keys == "cff-version"] <- "parameter on function"
origin[keys == "type"] <- "Fixed value: 'software'"
origin[keys == "identifiers"] <- "DESCRIPTION/CITATION files"
origin[keys == "references"] <- "DESCRIPTION/CITATION files"

origin[keys %in% c(
  "message",
  "title",
  "version",
  "authors",
  "abstract",
  "repository",
  "repository-code",
  "url",
  "date-released",
  "contact",
  "keywords",
  "license"
)] <- "DESCRIPTION file"

origin[keys %in% c(
  "doi",
  "preferred-citation"
)] <- "CITATION file"


origin[origin == FALSE] <- "Not parsed by cffr"

df <- data.frame(
  key = paste0("<a href='#", keys, "'>", keys, "</a>"),
  source = origin
)


knitr::kable(df, escape = FALSE)
```

## Details

### abstract

This key is extracted from the "Description" field of the DESCRIPTION file.

<details>

<summary>

<strong>Example</strong>

</summary>

```{r abstract}

library(cffr)

# Create cffr for yaml

cff_obj <- cff_create("rmarkdown")

# Get DESCRIPTION of rmarkdown to check

pkg <- desc::desc(file.path(find.package("rmarkdown"), "DESCRIPTION"))

cat(cff_obj$abstract)

cat(pkg$get("Description"))
```

</details>

[Back to summary](#summary).

### authors

This key is parsed from the "Authors" or "Authors\@R" field of the DESCRIPTION
file. Only persons with the role "aut" or "cre" are considered.

<details>

<summary>

<strong>Example</strong>

</summary>

```{r authors}

cff_obj <- cff_create("rmarkdown")
pkg <- desc::desc(file.path(find.package("rmarkdown"), "DESCRIPTION"))

cff_obj$authors

authors <- pkg$get_authors()

authors[vapply(authors, function(x) {
  "aut" %in% x$role || "cre" %in% x$role
}, logical(1))]
```

</details>

[Back to summary](#summary).

### cff-version

This key can be set via the parameters of the `cff_create()`/`cff_write()`
functions:

<details>

<summary>

<strong>Example</strong>

</summary>

```{r cffversion}

cff_objv110 <- cff_create("jsonlite", cff_version = "v1.1.0")

cat(cff_objv110$`cff-version`)
```

</details>

[Back to summary](#summary).

### commit

This key is not extracted from the metadata of the package. See the description
on the [Guide to CFF schema
v1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#commit).

> -   **description**: The commit hash or revision number of the software
>     version.
>
> -   **usage**:<br><br>
>
>     ``` yaml
>     commit: 1ff847d81f29c45a3a1a5ce73d38e45c2f319bba
>
>     commit: "Revision: 8612"
>     ```

[Back to summary](#summary).

### contact

This key is parsed from the "Authors" or "Authors\@R" field of the DESCRIPTION
file. Only persons with the role "cre" (i.e, the maintainer(s)) are considered.

<details>

<summary>

<strong>Example</strong>

</summary>

```{r contact}

cff_obj <- cff_create("rmarkdown")
pkg <- desc::desc(file.path(find.package("rmarkdown"), "DESCRIPTION"))

cff_obj$contact

pkg$get_author()
```

</details>

[Back to summary](#summary).

### date-released

This key is parsed from the "Date" field or, if not present, from the
"Date/Publication" field that is present on packages built on CRAN.

<details>

<summary>

<strong>Example</strong>

</summary>

```{r date-released}

# From an installed package

cff_obj <- cff_create("rmarkdown")
pkg <- desc::desc(file.path(find.package("rmarkdown"), "DESCRIPTION"))


cat(pkg$get("Date/Publication"))


cat(cff_obj$`date-released`)



# A DESCRIPTION file without a Date
nodate <- system.file("examples/DESCRIPTION_basic", package = "cffr")
tmp <- tempfile("DESCRIPTION")

# Create a temporary file
file.copy(nodate, tmp)


pkgnodate <- desc::desc(tmp)
cffnodate <- cff_create(tmp)

# Won't appear
cat(cffnodate$`date-released`)

pkgnodate

# Adding a Date

desc::desc_set("Date", "1999-01-01", file = tmp)

cat(cff_create(tmp)$`date-released`)
```

</details>

[Back to summary](#summary).

### doi {#doi}

This key is parsed from the "doi" field of the
[preferred-citation](#preferred-citation) object.

<details>

<summary>

<strong>Example</strong>

</summary>

```{r doi}

cff_doi <- cff_create("cffr")

cat(cff_doi$doi)

cat(cff_doi$`preferred-citation`$doi)
```

</details>

[Back to summary](#summary).

### identifiers

This key includes all the possible identifiers of the package:

-   From the DESCRIPTION field, it includes all the urls not included in
    [url](#url) or [repository-code](#repository-code).

-   From the CITATION file, it includes all the dois not included in [doi](#doi)
    and the identifiers (if any) not included in the "identifiers" key of
    [preferred-citation](#preferred-citation).

<details>

<summary>

<strong>Example</strong>

</summary>

```{r identifiers}
file <- system.file("examples/DESCRIPTION_many_urls", package = "cffr")

pkg <- desc::desc(file)

cat(pkg$get_urls())

cat(cff_create(file)$url)

cat(cff_create(file)$`repository-code`)

cff_create(file)$identifiers
```

</details>

[Back to summary](#summary).

### keywords

This key is extracted from the DESCRIPTION file. The keywords should appear in
the DESCRIPTION as:

    ...
    X-schema.org-keywords: keyword1, keyword2, keyword3

<details>

<summary>

<strong>Example</strong>

</summary>

```{r keyword}

# A DESCRIPTION file without keywords
nokeywords <- system.file("examples/DESCRIPTION_basic", package = "cffr")
tmp2 <- tempfile("DESCRIPTION")

# Create a temporary file
file.copy(nokeywords, tmp2)


pkgnokeywords <- desc::desc(tmp2)
cffnokeywords <- cff_create(tmp2)

# Won't appear
cat(cffnokeywords$keywords)

pkgnokeywords

# Adding Keywords

desc::desc_set("X-schema.org-keywords", "keyword1, keyword2, keyword3", file = tmp2)

cat(cff_create(tmp2)$keywords)
```

</details>

Additionally, if the source code of the package is hosted on GitHub, **cffr**
can retrieve the topics of your repo via the [GitHub
API](https://docs.github.com/en/rest) and include those topics as keywords. This
option is controlled via the `gh_keywords` parameter:

<details>

<summary>

<strong>Example</strong>

</summary>

```{r ghkeyword}

# Get cff object from jsonvalidate

jsonval <- cff_create("jsonvalidate")

# Keywords are retrieved from the GitHub repo

jsonval

# Check keywords
jsonval$keywords

# The repo
jsonval$`repository-code`
```

</details>

[Back to summary](#summary).

### license

This key is extracted from the "License" field of the DESCRIPTION file.

<details>

<summary>

<strong>Example</strong>

</summary>

```{r license}

cff_obj <- cff_create("yaml")

cat(cff_obj$license)

pkg <- desc::desc(file.path(find.package("yaml"), "DESCRIPTION"))

cat(pkg$get("License"))
```

</details>

[Back to summary](#summary).

### license-url

This key is not extracted from the metadata of the package. See the description
on the [Guide to CFF schema
v1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#license-url).

> -   **description**: The URL of the license text under which the software or
>     dataset is licensed (only for non-standard licenses not included in the
>     [SPDX License
>     List](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#definitionslicense-enum)).
> -   **usage**:<br><br>
>     `yaml     license-url: "https://obscure-licenses.com?id=1234"`

[Back to summary](#summary).

### message

This key is extracted from the DESCRIPTION field, specifically as:

```{r eval=FALSE}
msg <- paste0(
  'To cite package "',
  "NAME_OF_THE_PACKAGE",
  '" in publications use:'
)
```

<details>

<summary>

<strong>Example</strong>

</summary>

```{r message}

cat(cff_create("jsonlite")$message)
```

</details>

[Back to summary](#summary).

### preferred-citation {#preferred-citation}

This key is extracted from the CITATION file. If several references are
provided, it would select the first citation as the "preferred-citation" and the
rest of them as [references](#references).

<details>

<summary>

<strong>Example</strong>

</summary>

```{r preferred-citation}

cffobj <- cff_create("rmarkdown")

cffobj$`preferred-citation`

citation("rmarkdown")[1]
```

</details>

[Back to summary](#summary).

### references {#references}

This key is extracted from the CITATION file if several references are provided.
The first citation is considered as the
[preferred-citation](#preferred-citation) and the rest of them as "references".
It also extracts the package dependencies and adds those to this fields using
`citation(auto = TRUE)` on each dependency.

<details>

<summary>

<strong>Example</strong>

</summary>

```{r references}

cffobj <- cff_create("rmarkdown")

cffobj$references

citation("rmarkdown")[-1]
```

</details>

[Back to summary](#summary).

### repository

This key is extracted from the "Repository" field of the DESCRIPTION file.
Usually, this field is auto-populated when a package is hosted on a repo (like
CRAN or the [r-universe](https://r-universe.dev/)). For packages without this
field on the DESCRIPTION (that is the typical case for an in-development
package), **cffr** would try to search the package on any of the default
repositories specified on `options("repos")`.

In the case of [Bioconductor](https://bioconductor.org/) packages, those are
identified if a
["biocViews"](http://contributions.bioconductor.org/description.html#biocviews)
is present on the DESCRIPTION file.

If **cffr** detects that the package is available on CRAN, it would return the
canonical url form of the package (i.e.
<https://CRAN.R-project.org/package=jsonlite>).

<details>

<summary>

<strong>Example</strong>

</summary>

```{r repository}

# Installed package

inst <- cff_create("jsonlite")

cat(inst$repository)

# Demo file downloaded from the r-universe

runiv <- system.file("examples/DESCRIPTION_r_universe", package = "cffr")
runiv_cff <- cff_create(runiv)

cat(runiv_cff$repository)

desc::desc(runiv)$get("Repository")

# For in development package

norepo <- system.file("examples/DESCRIPTION_basic", package = "cffr")

# No repo
norepo_cff <- cff_create(norepo)

cat(norepo_cff[["repository"]])

# Change the name to a known package on CRAN: ggplot2

tmp <- tempfile("DESCRIPTION")
file.copy(norepo, tmp)


# Change name
desc::desc_set("Package", "ggplot2", file = tmp)

cat(cff_create(tmp)[["repository"]])

# Show what happens if another repo is set

# Save original config
orig_options <- options()
getOption("repos")


# Set new repos
options(repos = c(
  tidyverse = "https://tidyverse.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

# Load again the library
# Repos are evaluated on load
unloadNamespace("cffr")
library(cffr)


cat(cff_create(tmp)[["repository"]])

# Now it is the tidyverse repo, due to our new config!

# Reset original config
options(orig_options)
getOption("repos")
```

</details>

[Back to summary](#summary).

### repository-artifact

This key is not extracted from the metadata of the package. See the description
on the [Guide to CFF schema
v1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#repository-artifact).

> -   **description**: The URL of the work in a build artifact/binary repository
>     (when the work is software).
>
> -   **usage**:<br><br>
>
>     ``` yaml
>     repository-artifact: "https://search.maven.org/artifact/org.corpus-tools/cff-maven-plugin/0.4.0/maven-plugin"
>     ```

[Back to summary](#summary).

### repository-code {#repository-code}

This key is extracted from the "BugReports" or "URL" fields on the DESCRIPTION
file. **cffr** tries to identify the url of the source on the following
repositories:

-   [GitHub](https://github.com/).
-   [GitLab](https://about.gitlab.com/).
-   [R-Forge](https://r-forge.r-project.org/).
-   [Bitbucket](https://bitbucket.org/).

<details>

<summary>

<strong>Example</strong>

</summary>

```{r repository-code}

# Installed package on GitHub

cff_create("jsonlite")$`repository-code`



# GitLab

gitlab <- system.file("examples/DESCRIPTION_gitlab", package = "cffr")

cat(cff_create(gitlab)$`repository-code`)


# Check

desc::desc(gitlab)
```

</details>

[Back to summary](#summary).

### title

This key is extracted from the "Description" field of the DESCRIPTION file.

```{r eval=FALSE}
title <- paste0(
  "NAME_OF_THE_PACKAGE",
  ": ",
  "TITLE_OF_THE_PACKAGE"
)
```

<details>

<summary>

<strong>Example</strong>

</summary>

```{r title}

# Installed package

cat(cff_create("testthat")$title)
```

</details>

[Back to summary](#summary).

### type

Fixed value equal to "software". The other possible value is "dataset". See the
description on the [Guide to CFF schema
v1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md#type).

[Back to summary](#summary).

### url {#url}

This key is extracted from the "BugReports" or "URL" fields on the DESCRIPTION
file. It corresponds to the first url that is different to
[repository-code](#repository-code).

<details>

<summary>

<strong>Example</strong>

</summary>

```{r url}

# Many urls
manyurls <- system.file("examples/DESCRIPTION_many_urls", package = "cffr")

cat(cff_create(manyurls)$url)

# Check

desc::desc(manyurls)
```

</details>

[Back to summary](#summary).

### version

This key is extracted from the "Version" field on the DESCRIPTION file.

```{r version}

# Should be (>= 3.0.0)
cat(cff_create("testthat")$version)
```

[Back to summary](#summary).

## References
---
title: "Manipulating Citations with cffr"
description: >
  Learn how to modify `cff` objects.
output: 
  rmarkdown::html_vignette:
    toc: true
bibliography: REFERENCES.bib
link-citations: yes
vignette: >
  %\VignetteIndexEntry{Manipulating Citations with cffr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = TRUE
)

library(cffr)
```

**cffr** is a tool whose target audience are **R** package developers. The main
goal of **cffr** is to create a `CITATION.cff` file using the metadata
information of the following files:

-   Your `DESCRIPTION` file.
-   If available, the citation information located in `inst/CITATION`.

## What is a `CITATION.cff` file?

[Citation File Format (CFF](https://citation-file-format.github.io/))
[@druskat_citation_2021] (v1.2.0) are plain text files with human- and
machine-readable citation information for software (and datasets). Code
developers can include them in their repositories to let others know how to
correctly cite their software.

This format is becoming popular within the software citation ecosystem. Recently
[GitHub](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files),
[Zenodo](https://twitter.com/ZENODO_ORG/status/1420357001490706442) and
[Zotero](https://twitter.com/zotero/status/1420515377390530560) have included
full support of this citation format [@druskat_stephan_making_2021].

GitHub support is of special interest:

```{r echo=FALSE, out.width="400", fig.align='center', fig.alt="GitHub-link"}
knitr::include_graphics("tweet-1.png")
```

*--- Nat Friedman (\@natfriedman) [July 27,
2021](https://twitter.com/natfriedman/status/1420122675813441540?ref_src=twsrc%5Etfw)*

See [Customize your repository/About CITATION
files](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)
for more info.

## Creating a `CITATION.cff` file for my R package

With **cffr** creating a `CITATION.cff` file is quite straightforward. You just
need to run `cff_write()`:

```{r setup, eval=FALSE}

library(cffr)

cff_write()

# You are done!
```

Under the hood, `cff_write()` performs the following tasks:

-   It extracts the metadata using `cff_create()`.
-   Writes a `CITATION.cff` file using `yaml::write_yaml()`.
-   Validates the result using `cff_validate()`.

Congratulations! Now you have a full `CITATION.cff` file for your **R** package.

## Modifying your `CITATION.cff` file

You can easily customize the `cff` object (a custom class of **cffr**) using the
parsers provided in the package, as well as making use of the `keys` parameter.

We would create a `cff` object using `cff()` (for example purposes only) and we
would add or modify contents of it.

### Adding new fields

```{r newfields}

newobject <- cff_create(cff())

# For modifying your auto-generated object, run this line instead:
# newoobject <- cff_create()

newobject
```

The valid keys of the [Citation File Format schema version
1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md)
can be displayed with `cff_schema_keys()`:

```{r validkeys}

cff_schema_keys()
```

In this case, we are going to add `url`, `version` and `repository`. We would
also overwrite the `title` key. We just need to create a list and pass it to the
`keys` argument of `cff_create()`:

```{r modify}

newkeys <- list(
  "url" = "https://ropensci.org/",
  "version" = "0.0.1",
  "repository" = "https://github.com/user/repo",
  # If the field is already present, it would be overridden
  title = "Modifying a 'cff' object"
)

modobject <- cff_create(newobject, keys = newkeys)

modobject

# Validate against the schema

cff_validate(modobject)
```

### Parsing persons and citations

**cffr** provides two functions that parse `person` objects and `bibentry`
objects (See `?person` and `?bibentry`). These objects are included in the
**utils** package and are a core part of the metadata of any **R** package.

Following the previous example, we are going to add a new author first. For
doing that, we need first to extract the current author of the package and
append the parsed person:

```{r includeauthor}

# Valid person keys

cff_schema_definitions_person()

# Create the person

chiquito <- person("Gregorio",
  "Sánchez Fernández",
  email = "fake@email2.com",
  comment = c(
    alias = "Chiquito de la Calzada",
    city = "Malaga",
    country = "ES",
    ORCID = "0000-0000-0000-0001"
  )
)

chiquito

# Parse it
chiquito_parsed <- cff_parse_person(chiquito)
chiquito_parsed


# Append to previous authors

# Needs to be append as a list
newauthors <- c(modobject$authors, list(chiquito_parsed))
newauthors

newauthorobject <- cff_create(modobject, keys = list(authors = newauthors))

newauthorobject

cff_validate(newauthorobject)
```

Now, we may want to add `references` to our data. **cffr** supports two types of
references:

-   References created with `bibentry()`
-   References extracted from packages using `citation()`

On the following example, we would add two references, one of each type:

```{r parsingcits}
# Valid reference keys

cff_schema_definitions_refs()

# Auto parsed from another R package
base_r <- cff_parse_citation(citation("base"))

base_r

# Create with bibentry

bib <- bibentry("Book",
  title = "This is a book",
  author = "Lisa Lee",
  year = 1980,
  publisher = "McGraw Hill",
  volume = 2
)
bib

# Now parse it

bookparsed <- cff_parse_citation(bib)

bookparsed
```

Now the process is similar to the example with `person`: we append both
references (as lists) and add them to our object:

```{r references}

refkeys <- list(references = c(list(base_r), list(bookparsed)))

refkeys

finalobject <- cff_create(newauthorobject, keys = refkeys)

finalobject

cff_validate(finalobject)
```

### Create your modified `CITATION.cff` file

The results can be written with `cff_write()`:

```{r write}

# For example
tmp <- tempfile(fileext = ".cff")

see_res <- cff_write(finalobject, outfile = tmp)

see_res
```

And finally we can read our created `CITATION.cff` file using `cff()`:

```{r read}

reading <- cff(tmp)

reading
```

Note that `cff_write()` also has the `keys` param, so the workflow can be
simplified as:

```{r}

allkeys <- list(
  "url" = "https://ropensci.org/",
  "version" = "0.0.1",
  "repository" = "https://github.com/user/repo",
  # If the field is already present, it would be overridden
  title = "Modifying a 'cff' object",
  authors = newauthors,
  references = c(list(base_r), list(bookparsed))
)

tmp2 <- tempfile(fileext = ".cff")

res <- cff_write(cff(), outfile = tmp2, keys = allkeys)

res
```

## References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_extract_to_bibtex.R
\name{cff_extract_to_bibtex}
\alias{cff_extract_to_bibtex}
\title{Create BibTeX entries from a package}
\usage{
cff_extract_to_bibtex(x, what = "preferred")
}
\arguments{
\item{x}{The source that would be used for generating
the \code{\link{cff}} object. It could be:
\itemize{
\item A missing value. That would retrieve the DESCRIPTION
file on your in-development package.
\item An existing \code{\link{cff}} object,
\item The name of an installed package (\code{"jsonlite"}), or
\item Path to a DESCRIPTION file (\code{"*/DESCRIPTION*"}).
}}

\item{what}{Fields to extract. The value could be:
\itemize{
\item \code{preferred}: This would create a single entry with the main citation
info of the package.
\item \code{references}: Extract all the entries on \code{references}.
\item \code{all}: A combination of the previous two options. This would extract
both the preferred citation info and the references.
}}
}
\value{
A \code{bibentry} object or a list of \code{bibentry} objects. This could
be parsed to BibTeX using \code{\link[=toBibtex]{toBibtex()}}
}
\description{
Extract the information of a package to BibTeX. This is done by creating a
\code{cff} object with \code{cff_create()} and extracting the corresponding entries
with \code{cff_to_bibtex()}.
}
\examples{
\donttest{

jsonvalidate <- cff_extract_to_bibtex("jsonvalidate")

jsonvalidate

toBibtex(jsonvalidate)

lite <- cff_extract_to_bibtex("jsonlite", "references")

lite

toBibtex(lite)
}
}
\seealso{
Other bibtex: 
\code{\link{cff_to_bibtex}()},
\code{\link{encoded_utf_to_latex}()},
\code{\link{write_bib}()}
}
\concept{bibtex}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-bibtex.R
\name{encoded_utf_to_latex}
\alias{encoded_utf_to_latex}
\title{Encode UTF-8 text to LaTeX}
\usage{
encoded_utf_to_latex(x)
}
\arguments{
\item{x}{A string, possibly encoded in UTF-8 encoding system.}
}
\value{
A string with the corresponding transformations.
}
\description{
Transform a UTF-8 string to LaTeX special characters.
}
\details{
This is a variation of \code{\link[tools:encoded]{tools::encoded_text_to_latex()}} performing some
additional replacements to increase compatibility.
}
\examples{
\dontshow{if (getRversion() >= "4.0.0") (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Full range of supported characters on R
library(tools)

range <- 1:511

ascii_table <- data.frame(
  dec = range,
  utf8 = intToUtf8(range, multiple = TRUE)
)

# Add latex using base approach
ascii_table$latex_base <- encoded_text_to_latex(ascii_table$utf8,
  encoding = "UTF-8"
)

# With cffr
ascii_table$latex_cffr <- encoded_utf_to_latex(ascii_table$utf8)

ascii_table
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link[tools:encoded]{tools::encoded_text_to_latex()}}

Other bibtex: 
\code{\link{cff_extract_to_bibtex}()},
\code{\link{cff_to_bibtex}()},
\code{\link{write_bib}()}
}
\concept{bibtex}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff.R
\name{cff}
\alias{cff}
\alias{as.cff}
\title{\code{cff} objects}
\usage{
cff(path, ...)

as.cff(x)
}
\arguments{
\item{path}{The path to a \code{CITATION.cff} file.}

\item{...}{Named arguments to be used for creating a \code{\link{cff}} object. See
\strong{Details}.}

\item{x}{a character string for the \code{\link{as.cff}} default method}
}
\value{
A \code{cff} object. Under the hood, a \code{cff} object is a regular \code{\link{list}}
object with a special \code{\link[=print]{print()}} method.
}
\description{
A class and utility methods for reading, creating and holding CFF
information.
}
\details{
This object can be manipulated using \code{\link[=cff_create]{cff_create()}}.

\strong{Note that} this function reads \code{CITATION.cff} files. If you want to
create \code{cff} objects from DESCRIPTION files use \code{\link[=cff_create]{cff_create()}}.

If no additional \code{...} parameters are supplied (the default behavior),
a minimal valid \code{cff} object is created. Valid parameters are those
specified on \code{\link[=cff_schema_keys]{cff_schema_keys()}}:\tabular{l}{
   \strong{valid cff keys} \cr
   cff-version \cr
   message \cr
   type \cr
   license \cr
   title \cr
   version \cr
   doi \cr
   abstract \cr
   authors \cr
   preferred-citation \cr
   repository \cr
   repository-artifact \cr
   repository-code \cr
   url \cr
   date-released \cr
   contact \cr
   keywords \cr
   references \cr
   commit \cr
   identifiers \cr
   license-url \cr
}
}
\examples{

# Blank cff
cff()

# From file
cff(system.file("examples/CITATION_basic.cff",
  package = "cffr"
))

# Use custom params
test <- cff(
  title = "Manipulating files",
  keywords = c("A", "new", "list", "of", "keywords"),
  authors = list(cff_parse_person("New author"))
)
test
\donttest{
# Would fail
cff_validate(test)


# Modify with cff_create
new <- cff_create(test, keys = list(
  "cff-version" = "1.2.0",
  message = "A blank file"
))
new

# Would pass
cff_validate(new)
}


# Convert a list to "cff" object
cffobj <- as.cff(list(
  "cff-version" = "1.2.0",
  title = "Manipulating files"
))

class(cffobj)

# Nice display thanks to yaml package
cffobj
}
\seealso{
Other core functions: 
\code{\link{cff_create}()},
\code{\link{cff_validate}()},
\code{\link{cff_write}()}
}
\concept{core functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_parse_person.R
\name{cff_parse_person}
\alias{cff_parse_person}
\alias{cff_parse_person_bibtex}
\title{Parse a person to \code{cff}}
\usage{
cff_parse_person(person)

cff_parse_person_bibtex(person)
}
\arguments{
\item{person}{A \code{person} object created with \code{\link[=person]{person()}} or a character string.
See \strong{Details}.}
}
\value{
A \code{\link{cff}} object ready to be used on \code{\link[=cff_create]{cff_create()}}.
}
\description{
Parse a person or string to a valid format for a \code{CITATION.cff} file. This
is a helper function designed to help on adding or replacing the
auto-generated authors of the package.
}
\details{
The \code{person} parameter of the function could be:
\itemize{
\item For \code{cff_parse_person()}: A \code{person} object or a character coercible to
\code{person}. See \code{\link[=person]{person()}} for details.
\item For \code{cff_parse_person_bibtex()}: A string with the definition of an author
or several authors, using the standard BibTeX notation. See Markey (2007)
for a full explanation.
}

See \strong{Examples} for more information.
}
\examples{
# Parse a person object

cff_parse_person(person(
  given = "First",
  family = "Author",
  role = c("aut", "cre"),
  email = "first.last@example.com",
  comment = c(
    ORCID = "0000-0001-8457-4658",
    affiliation = "An affiliation"
  )
))

# Parse a string

cff_parse_person("Julio Iglesias <fake@email.com>")

# Several persons
persons <- c(person("Clark", "Kent"), person("Lois", "Lane"))

cff_parse_person(persons)

# Or you can use BibTeX style if you prefer

x <- "Frank Sinatra and Dean Martin and Davis, Jr., Sammy and Joey Bishop"

cff_parse_person_bibtex(x)

cff_parse_person_bibtex("Herbert von Karajan")
}
\references{
\itemize{
\item Patashnik, Oren. "BIBTEXTING" February 1988.
\url{https://osl.ugr.es/CTAN/biblio/bibtex/base/btxdoc.pdf}.
\item Markey, Nicolas. "Tame the BeaST."
\emph{The B to X of BibTeX, Version 1.4} (October 2007).
\url{https://osl.ugr.es/CTAN/info/bibtex/tamethebeast/ttb_en.pdf}.
}
}
\seealso{
\code{\link[=cff_create]{cff_create()}}, \code{vignette("cffr", "cffr")}, \code{\link[utils:person]{utils::person()}}

Other parsers: 
\code{\link{cff_parse_citation}()}
}
\concept{parsers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_to_bibtex.R
\name{cff_to_bibtex}
\alias{cff_to_bibtex}
\title{Create a BibTeX entry from a CITATION file or a \code{cff} object}
\usage{
cff_to_bibtex(x)
}
\arguments{
\item{x}{The source that would be used for generating
the \code{\link{cff}} object. It could be:
\itemize{
\item An existing \code{\link{cff}} object,
\item A CITATION.cff file.
}}
}
\value{
A \code{bibentry} object that can be parsed to BibTeX format with
\code{\link[=toBibtex]{toBibtex()}}
}
\description{
Creates a \code{bibentry} object (\code{\link[=bibentry]{bibentry()}}) from a \code{cff} object
}
\examples{
\donttest{

# From a cff object
package <- cff_create("rmarkdown")

obj <- cff_to_bibtex(package)

obj

toBibtex(obj)
}
}
\references{
\itemize{
\item Patashnik, Oren. "BIBTEXTING" February 1988.
\url{https://osl.ugr.es/CTAN/biblio/bibtex/base/btxdoc.pdf}.
\item Haines, R., & The Ruby Citation File Format Developers. (2021).
\emph{Ruby CFF Library (Version 0.9.0)} (Computer software).
\doi{10.5281/zenodo.1184077}.
}
}
\seealso{
\code{\link[=cff_parse_citation]{cff_parse_citation()}}, \code{\link[=bibentry]{bibentry()}}, \code{\link[=toBibtex]{toBibtex()}}

Other bibtex: 
\code{\link{cff_extract_to_bibtex}()},
\code{\link{encoded_utf_to_latex}()},
\code{\link{write_bib}()}
}
\concept{bibtex}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_write.R
\name{cff_write}
\alias{cff_write}
\title{Write a \code{CITATION.cff} file}
\usage{
cff_write(
  x,
  outfile = "CITATION.cff",
  keys = list(),
  cff_version = "1.2.0",
  gh_keywords = TRUE,
  dependencies = TRUE,
  validate = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{x}{The source that would be used for generating
the \code{CITATION.cff} file. It could be:
\itemize{
\item A missing value. That would retrieve the DESCRIPTION
file on your in-development package.
\item A \code{\link{cff}} object,
\item The name of an installed package (\code{"jsonlite"}), or
\item Path to a DESCRIPTION file (\code{"*/DESCRIPTION*"}).
}}

\item{outfile}{The name and path of the \code{CITATION.cff} to be created.}

\item{keys}{List of additional keys to add to the \code{\link{cff}} object. See
\code{\link[=cff_create]{cff_create()}} for details and examples.}

\item{cff_version}{The Citation File Format schema version that the
\code{CITATION.cff} file adheres to for providing the citation metadata.}

\item{gh_keywords}{Logical \code{TRUE/FALSE}. If the package is hosted on
GitHub, would you like to add the repo topics as keywords?}

\item{dependencies}{Logical \code{TRUE/FALSE}. Would you like to add the
of your package to the \code{reference} key?}

\item{validate}{Logical \code{TRUE/FALSE}. Should the new file be validated using
\code{\link[=cff_validate]{cff_validate()}}?}

\item{verbose}{Logical \code{TRUE/FALSE}. On \code{TRUE} the function would display
informative messages.}
}
\value{
A \code{CITATION.cff} file and an (invisible) \code{\link{cff}} object.
}
\description{
\strong{This is the core function of the package and likely to be the only one
you would need when developing a package}.

This function writes out a \code{CITATION.cff} file for a given package. This
function is basically a wrapper around \code{\link[=cff_create]{cff_create()}} to both create the
\code{\link{cff}} object and writes it out to a YAML-formatted file in one command.
}
\details{
When creating and writing a \code{CITATION.cff} for the first time, the function
adds "CITATION.cff" to ".Rbuildignore".
}
\examples{
\donttest{
tmpfile <- tempfile(fileext = ".cff")
cff_obj <- cff_write("jsonlite", outfile = tmpfile)

cff_obj

# Force clean-up
file.remove(tmpfile)
}
}
\seealso{
\href{https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md}{Guide to Citation File Format schema version 1.2.0}.

Other core functions: 
\code{\link{cff_create}()},
\code{\link{cff_validate}()},
\code{\link{cff}()}
}
\concept{core functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cffr-package.R
\docType{package}
\name{cffr-package}
\alias{cffr}
\alias{cffr-package}
\title{cffr: Generate Citation File Format ('cff') Metadata for R Packages}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

The Citation File Format version 1.2.0 <doi:10.5281/zenodo.5171937> is a human and machine readable file format which provides citation metadata for software. This package provides core utilities to generate and validate this metadata.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/cffr/}
  \item \url{https://github.com/ropensci/cffr}
  \item Report bugs at \url{https://github.com/ropensci/cffr/issues}
}

}
\author{
\strong{Maintainer}: Diego Hernangómez \email{diego.hernangomezherrero@gmail.com} (\href{https://orcid.org/0000-0001-8457-4658}{ORCID}) [copyright holder]

Other contributors:
\itemize{
  \item João Martins (\href{https://orcid.org/0000-0001-7961-4280}{ORCID}) [reviewer]
  \item Scott Chamberlain (\href{https://orcid.org/0000-0003-1444-9135}{ORCID}) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cran_to_spdx}
\alias{cran_to_spdx}
\title{Mapping between \code{License} fields and SPDX}
\format{
A data frame with 91 rows and 2 variables:
\itemize{
\item LICENSE: A valid \code{License} string on CRAN.
\item SPDX. A valid SPDX License Identifier.
}
}
\source{
\url{https://spdx.org/licenses/}
}
\usage{
cran_to_spdx
}
\description{
A dataset containing the mapping between the \code{License} strings observed
on CRAN packages and its (approximate) match on the
\href{https://spdx.org/licenses/}{SPDX License List}.
}
\examples{

data("cran_to_spdx")

head(cran_to_spdx, 20)
}
\seealso{
\emph{Writing R Extensions}, \href{https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Licensing}{Licensing section}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_validate.R
\name{cff_validate}
\alias{cff_validate}
\title{Validate a \code{CITATION.cff} file or a \code{\link{cff}} object}
\usage{
cff_validate(x = "CITATION.cff", verbose = TRUE)
}
\arguments{
\item{x}{This is expected to be either a \code{\link{cff}} object created
with \code{\link[=cff_create]{cff_create()}} or the path to a \code{CITATION.cff} file to be validated.}

\item{verbose}{Logical \code{TRUE/FALSE}. On \code{TRUE} the function would display
informative messages.}
}
\value{
A message indicating the result of the validation and an invisible
value \code{TRUE/FALSE}.
}
\description{
Validate a \code{CITATION.cff} file or a \code{\link{cff}} object created with
\code{\link[=cff_create]{cff_create()}} using the corresponding validation
\href{https://github.com/citation-file-format/citation-file-format/blob/main/schema.json}{schema.json}.
}
\examples{
\donttest{
# Full .cff example
cff_validate(system.file("examples/CITATION_complete.cff", package = "cffr"))

# Validate a cffr object
cffr <- cff_create("jsonlite")
class(cffr)
cff_validate(cffr)
}
\dontrun{
# .cff with errors
cff_validate(system.file("examples/CITATION_error.cff", package = "cffr"))
# If a CITATION file (note that is not .cff) it throws an error
cff_validate(system.file("CITATION", package = "cffr"))
}
}
\seealso{
\href{https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md}{Guide to Citation File Format schema version 1.2.0}.

Other core functions: 
\code{\link{cff_create}()},
\code{\link{cff_write}()},
\code{\link{cff}()}
}
\concept{core functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_git_hook.R
\name{cff_git_hook}
\alias{cff_git_hook}
\alias{cff_git_hook_install}
\alias{cff_git_hook_remove}
\title{Use a git pre-commit hook \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}}
\usage{
cff_git_hook_install()

cff_git_hook_remove()
}
\value{
Invisible. This function is called for its side effects.
}
\description{
Install a
\href{https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks#_committing_workflow_hooks}{pre-commit hook}
that remembers you to update your \code{CITATION.cff} file.
}
\details{
This function would install a pre-commit hook using
\code{\link[usethis:use_git_hook]{usethis::use_git_hook()}}.

A pre-commit hook is a script that identifies  simple issues before
submission to code review. This pre-commit hook would warn you if any of the
following conditions are met:
\itemize{
\item You included in a commit your DESCRIPTION or inst/CITATION file, you
are not including your \code{CITATION.cff} and the \code{CITATION.cff} file is
"older" than any of your DESCRIPTION or inst/CITATION file, or
\item You have updated your \code{CITATION.cff} but you are not including it on
your commit.
}
}
\section{A word of caution}{
The pre-commit hook may prevent you to commit if you are not updating your
\code{CITATION.cff}. However, the mechanism of detection is not perfect and would
be triggered also even if you have tried to update your \code{CITATION.cff} file.

This is typically the case when you have updated your DESCRIPTION or
inst/CITATION files but those changes doesn't make a change on your
\code{CITATION.cff} file (i.e. you are including new dependencies).

In those cases, you can override the check running \verb{git commit --no-verify}
on the Terminal tab. If you are using
RStudio you can run also this command from a R script by selecting that
line and sending it to the Terminal using:
\itemize{
\item \code{Ctrl+Alt+Enter} (Windows & Linux), or
\item \code{Cmd+Option+Return} (Mac).
}
}

\section{Removing the git pre-commit hook}{
You can remove the pre-commit hook by running \code{cff_git_hook_remove()}.
}

\examples{
\dontrun{
cff_git_hook_install()
}

}
\seealso{
\code{\link[usethis:use_git_hook]{usethis::use_git_hook()}}, \code{\link[usethis:use_git]{usethis::use_git()}}

Other git: 
\code{\link{cff_gha_update}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_bib.R
\name{write_bib}
\alias{write_bib}
\title{Create a .bib file}
\usage{
write_bib(x, file = NULL, append = FALSE, verbose = TRUE, ascii = FALSE)
}
\arguments{
\item{x}{A \code{bibentry} object created with:
\itemize{
\item \code{\link[=cff_extract_to_bibtex]{cff_extract_to_bibtex()}}, \code{\link[=cff_to_bibtex]{cff_to_bibtex()}}
\item \code{\link[=citation]{citation()}} or \code{\link[=bibentry]{bibentry()}}
}}

\item{file}{Name of the file. If \code{NULL} it would display the lines to be
written.}

\item{append}{Whether to append the entries to an existing file or not.}

\item{verbose}{Display informative messages}

\item{ascii}{Whether to write the entries using ASCII characters only or not.}
}
\description{
Creates a .bib file from a \code{bibentry} object(s)
}
\details{
For security reasons, if the file already exists the function would create
a backup copy on the same directory.
}
\examples{

bib <- bibentry("Misc",
  title = "My title",
  author = "Fran Pérez"
)

write_bib(bib)

write_bib(bib, ascii = TRUE)
}
\seealso{
\code{\link[knitr:write_bib]{knitr::write_bib()}} and the following packages:
\itemize{
\item \href{https://github.com/ropensci/bibtex}{bibtex} package.
\item \href{https://github.com/ropensci/RefManageR}{RefManageR} package.
\item \href{https://github.com/GeoBosh/rbibutils/}{rbibutils}
}

Other bibtex: 
\code{\link{cff_extract_to_bibtex}()},
\code{\link{cff_to_bibtex}()},
\code{\link{encoded_utf_to_latex}()}
}
\concept{bibtex}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_gha_update.R
\name{cff_gha_update}
\alias{cff_gha_update}
\title{Install a cffr GitHub Action}
\usage{
cff_gha_update(path = ".", overwrite = FALSE)
}
\arguments{
\item{path}{Project directory}

\item{overwrite}{If already present, do you want to overwrite your action?}
}
\value{
Invisible, this function is called by its side effects.
}
\description{
This function would install a GitHub Action on your repo. The action
will update your \code{CITATION.cff} when any of these events occur:
\itemize{
\item You publish a new release of the package.
\item Your DESCRIPTION or inst/CITATION are modified.
\item The action can be run also manually.
}
}
\details{
Triggers on your action can be modified, see
\href{https://docs.github.com/en/actions/learn-github-actions/events-that-trigger-workflows}{Events that trigger workflows}.
}
\examples{
\dontrun{
cff_gha_update()
}
}
\seealso{
Other git: 
\code{\link{cff_git_hook}}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_create.R
\name{cff_create}
\alias{cff_create}
\title{Create \code{cff} object}
\usage{
cff_create(
  x,
  keys = list(),
  cff_version = "1.2.0",
  gh_keywords = TRUE,
  dependencies = TRUE
)
}
\arguments{
\item{x}{The source that would be used for generating
the \code{\link{cff}} object. It could be:
\itemize{
\item A missing value. That would retrieve the DESCRIPTION
file on your in-development package.
\item An existing \code{\link{cff}} object,
\item The name of an installed package (\code{"jsonlite"}), or
\item Path to a DESCRIPTION file (\code{"*/DESCRIPTION*"}).
}}

\item{keys}{List of additional keys to add to the \code{\link{cff}} object. See
\strong{Details}.}

\item{cff_version}{The Citation File Format schema version that the
\code{CITATION.cff} file adheres to for providing the citation metadata.}

\item{gh_keywords}{Logical \code{TRUE/FALSE}. If the package is hosted on
GitHub, would you like to add the repo topics as keywords?}

\item{dependencies}{Logical \code{TRUE/FALSE}. Would you like to add the
of your package to the \code{reference} key?}
}
\value{
A \code{\link{cff}} list object.
}
\description{
Create a \code{\link{cff}} object from a given source for further manipulation.
Similar to \code{\link[=cff_write]{cff_write()}}, but returns a object rather than writing
directly to a file. See \strong{Examples}.
}
\details{
It is possible to add additional keys not detected by \code{\link[=cff_create]{cff_create()}} using
the \code{keys} argument. A list of valid keys can be retrieved with
\code{\link[=cff_schema_keys]{cff_schema_keys()}}.

Please refer to
\href{https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md}{Guide to Citation File Format schema version 1.2.0}.
for additional details.

If \code{x} is a path to a DESCRIPTION file or \code{inst/CITATION}, is not present on
your package, \strong{cffr} would auto-generate a \code{preferred-citation} key using
the information provided on that file. On
}
\examples{
\donttest{
# Installed package
cff_create("jsonlite")

# Demo file
demo_file <- system.file("examples/DESCRIPTION_basic", package = "cffr")
cff_create(demo_file)

# Add additional keys

newkeys <- list(
  message = "This overwrites fields",
  abstract = "New abstract",
  keywords = c("A", "new", "list", "of", "keywords"),
  authors = list(cff_parse_person("New author"))
)

cff_create(demo_file, keys = newkeys)

# Update a field on a list - i,e: authors, contacts, etc.
# We are adding a new contact here

old <- cff_create(demo_file)

new_contact <- append(
  old$contact,
  list(
    cff_parse_person(person(
      given = "I am",
      family = "New Contact"
    ))
  )
)


cff_create(demo_file, keys = list("contact" = new_contact))
}
}
\seealso{
\href{https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md}{Guide to Citation File Format schema version 1.2.0}.

\code{vignette("cffr", "cffr")}

Other core functions: 
\code{\link{cff_validate}()},
\code{\link{cff_write}()},
\code{\link{cff}()}
}
\concept{core functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cff_parse_citation.R
\name{cff_parse_citation}
\alias{cff_parse_citation}
\title{Parse a \code{bibentry} to \code{cff}}
\usage{
cff_parse_citation(bib)
}
\arguments{
\item{bib}{A \code{bibentry} object, either created with \code{\link[=bibentry]{bibentry()}}
(preferred) or \code{\link[=citEntry]{citEntry()}}.}
}
\value{
A \code{\link{cff}} object ready to be used on \code{\link[=cff_create]{cff_create()}}.
}
\description{
Parse a \code{bibentry} object to a valid format for a \code{CITATION.cff} file.
}
\details{
This is a helper function designed to help on adding or
replacing the auto-generated authors of the package. See \strong{Examples}.

This function tries to adapt a \code{bibentry} object (generated with \code{\link[=bibentry]{bibentry()}}
or \code{\link[=citEntry]{citEntry()}}) to the CFF standard.
\subsection{Entry types considered}{
\itemize{
\item \strong{Article}, \strong{Book}, \strong{Booklet}, \strong{InBook}, \strong{InCollection},
\strong{InProceedings}, \strong{Manual}, \strong{MastersThesis}, \strong{Misc}, \strong{PhDThesis},
\strong{Proceedings}, \strong{TechReport}, \strong{Unpublished}. See \code{\link[=bibentry]{bibentry()}}
for more information.
}

Note that \strong{Conference} is not implemented in
\code{\link[=bibentry]{bibentry()}}, however is equivalent to \strong{InProceedings} (Patashnik (1988)).
}

\subsection{Fields considered}{
\itemize{
\item \strong{address}, \strong{author}, \strong{booktitle}, \strong{chapter}, \strong{edition},
\strong{editor}, \strong{howpublished}, \strong{institution},  \strong{journal}, \strong{key},
\strong{month}, \strong{note}, \strong{number}, \strong{organization}, \strong{pages},
\strong{publisher}, \strong{school}, \strong{series}, \strong{title}, \strong{type}, \strong{year}.
}

\strong{annote} and \strong{crossref} fields are ignored.
}
}
\examples{
\donttest{
bib <- citation("base")
bib


# To cff
bib_to_cff <- cff_parse_citation(bib)
bib_to_cff

# Create the object
new_cff <- cff()

full <- cff_create(new_cff, keys = list("preferred-citation" = bib_to_cff))

full
# Validate
cff_validate(full)

# Several citations

cff_parse_citation(citation("rmarkdown"))
}
}
\references{
\itemize{
\item Patashnik, Oren. "BIBTEXTING" February 1988.
\url{https://osl.ugr.es/CTAN/biblio/bibtex/base/btxdoc.pdf}.
\item Haines, R., & The Ruby Citation File Format Developers. (2021).
\emph{Ruby CFF Library (Version 0.9.0)} (Computer software).
\doi{10.5281/zenodo.1184077}.
}
}
\seealso{
\code{\link[=cff_create]{cff_create()}}, \code{vignette("cffr", "cffr")}, \code{\link[=bibentry]{bibentry()}}

Other parsers: 
\code{\link{cff_parse_person}()}
}
\concept{parsers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-schema.R
\name{cff_schema}
\alias{cff_schema}
\alias{cff_schema_keys}
\alias{cff_schema_keys_license}
\alias{cff_schema_definitions_person}
\alias{cff_schema_definitions_entity}
\alias{cff_schema_definitions_refs}
\title{Schema utils}
\source{
\href{https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md}{Guide to Citation File Format schema version 1.2.0}.
}
\usage{
cff_schema_keys(sorted = FALSE)

cff_schema_keys_license()

cff_schema_definitions_person()

cff_schema_definitions_entity()

cff_schema_definitions_refs()
}
\arguments{
\item{sorted}{Logical \code{TRUE/FALSE}. Should the keys be arranged
alphabetically?}
}
\value{
A vector of characters with the names of the valid keys to be used on a
Citation File Format version 1.2.0
}
\description{
Helper functions with the valid values of different fields, according to the
\href{https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md}{Citation File Format schema version 1.2.0}.
\itemize{
\item \code{\link[=cff_schema_keys]{cff_schema_keys()}} provides the valid high-level keys of the Citation
File Format.
\item \code{\link[=cff_schema_keys_license]{cff_schema_keys_license()}} provides the valid
\href{https://spdx.dev/ids/}{SPDX license identifier(s)} to be used on the
\code{CITATION.cff} file.
\item \code{\link[=cff_schema_definitions_person]{cff_schema_definitions_person()}} and \code{\link[=cff_schema_definitions_entity]{cff_schema_definitions_entity()}}
returns the valid fields to be included when defining a
person or entity.
\item \code{\link[=cff_schema_definitions_refs]{cff_schema_definitions_refs()}} provides the valid
keys to be used on the \code{preferred-citation} and \code{references} keys.
}
}
\examples{

cff_schema_keys(sorted = TRUE)

# Valid Licenses keys
head(cff_schema_keys_license(), 20)

cff_schema_definitions_person()

cff_schema_definitions_entity()

cff_schema_definitions_refs()
}
\concept{schema}
