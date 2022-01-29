
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
