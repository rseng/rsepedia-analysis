---
title: 'Generating CodeMeta Metadata for R Packages'
tags:
 - metadata
 - codemeta
 - ropensci
 - citation
 - credit
authors:
 - name: Carl Boettiger
   orcid: 0000-0002-1642-628X
   affiliation: 1
affiliations:
 - name: University of California, Berkeley
   index: 1
date: 2017-06-29
bibliography: paper.bib
---

# Summary

The CodeMeta project defines a JSON-LD format [@jsonld] for describing software
metadata, based largely on `schema.org` terms.
This metadata format is being adopted by many leading archives for scientific software, including DataCite,
Zenodo, and DataONE to address many of the needs identified in the NIH report on the need for a
"Software Discovery Index" [@SDI].
Many common software metadata formats have been mapped into CodeMeta by means of a crosswalk table [@codemeta], also implemented in this package.
The `codemetar` package provides utilities to generate and validate these `codemeta.json`
files automatically for R packages by parsing the DESCRIPTION file
and other common locations for R metadata.
The package also includes utilities and examples for parsing and working with existing codemeta files,
and includes several vignettes which illustrate both the basic usage of the package as well as some more advanced applications.

# References
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
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/

<!-- README.md is generated from README.Rmd. Please edit that file -->

# codemetar

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![R build
status](https://github.com/ropensci/codemetar/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/codemetar/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/ropensci/codemetar/master.svg)](https://codecov.io/github/ropensci/codemetar?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/codemetar)](https://cran.r-project.org/package=codemetar)
[![](http://badges.ropensci.org/130_status.svg)](https://github.com/ropensci/software-review/issues/130)
[![DOI](https://zenodo.org/badge/86626030.svg)](https://zenodo.org/badge/latestdoi/86626030)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/codemetar)](https://CRAN.R-project.org/package=codemetar)

The [codemeta](https://cran.r-project.org/package=codemeta) package
provides a more minimalist approach to generating codemeta based only on
DESCRIPTION and CITATION files, while `codemetar` provides additional
abilities to detect metadata from README and GitHub sources, and provide
more user feedback, suggestions, and messaging.

**Why codemetar?** The ‘Codemeta’ Project defines a ‘JSON-LD’ format for
describing software metadata, as detailed at
<https://codemeta.github.io>. This package provides utilities to
**generate, parse, and modify codemeta.jsonld files automatically for R
packages**, as well as tools and examples for **working with codemeta
json-ld more generally**.

It has three main goals:

-   Quickly **generate a valid codemeta.json file from any valid R
    package**. To do so, we automatically extract as much metadata as
    possible using the DESCRIPTION file, as well as extracting metadata
    from other common best-practices such as the presence of Travis and
    other badges in README, etc.
-   Facilitate the addition of further metadata fields into a
    codemeta.json file, as well as general manipulation of codemeta
    files.
-   Support the ability to crosswalk between terms used in other
    metadata standards, as identified by the Codemeta Project Community,
    see <https://codemeta.github.io/crosswalk/>

## Why create a codemeta.json for your package?

**Why bother creating a codemeta.json for your package?** R packages
encode lots of metadata in the `DESCRIPTION` file, `README`, and other
places, telling users and developers about the package purpose, authors,
license, dependencies, and other information that facilitates discovery,
adoption, and credit for your software. Unfortunately, because each
software language records this metadata in a different format, that
information is hard for search engines, software repositories, and other
developers to find and integrate.

By generating a `codemeta.json` file, you turn your metadata into a
format that can easily crosswalk between metadata in many other software
languages. CodeMeta is built on [schema.org](https://schema.org) a
simple [structured
data](https://developers.google.com/search/docs/advanced/structured-data/intro-structured-data)
format developed by major search engines like Google and Bing to improve
discoverability in search. CodeMeta is also understood by significant
software archiving efforts such as [Software
Heritage](https://www.softwareheritage.org/) Project, which seeks to
permanently archive all open source software.

For more general information about the CodeMeta Project for defining
software metadata, see <https://codemeta.github.io>. In particular, new
users might want to start with the [User
Guide](https://codemeta.github.io/user-guide/), while those looking to
learn more about JSON-LD and consuming existing codemeta files should
see the [Developer Guide](https://codemeta.github.io/developer-guide/).

## Create a codemeta.json in one function call

`codemetar` can take the path to the source package root to glean as
much information as possible.

``` r
codemetar::write_codemeta()
```

    … Getting CRAN metadata from RStudio CRAN mirror
    ✓ Got CRAN metadata!
    … Asking README URL from GitHub API
    ✓ Got README URL!
    … Asking README URL from GitHub API
    ✓ Got README URL!
    … Getting repo topics from GitHub API
    ✓ Got repo topics!

``` r
library("magrittr")
"../../codemeta.json" %>%
  details::details(summary = "codemetar's codemeta.json",
                   lang = "json")
```

<details closed>
<summary>
<span title="Click to Expand"> codemetar’s codemeta.json </span>
</summary>

``` json
{
  "@context": "https://doi.org/10.5063/schema/codemeta-2.0",
  "@type": "SoftwareSourceCode",
  "identifier": "codemetar",
  "description": "The 'Codemeta' Project defines a 'JSON-LD' format for describing software metadata, as detailed at <https://codemeta.github.io>. This package provides utilities to generate, parse, and modify 'codemeta.json' files automatically for R packages, as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more generally.",
  "name": "codemetar: Generate 'CodeMeta' Metadata for R Packages",
  "relatedLink": ["https://docs.ropensci.org/codemetar/", "https://CRAN.R-project.org/package=codemetar"],
  "codeRepository": "https://github.com/ropensci/codemetar",
  "issueTracker": "https://github.com/ropensci/codemetar/issues",
  "license": "https://spdx.org/licenses/GPL-3.0",
  "version": "0.3.3",
  "programmingLanguage": {
    "@type": "ComputerLanguage",
    "name": "R",
    "url": "https://r-project.org"
  },
  "runtimePlatform": "R version 4.1.0 (2021-05-18)",
  "provider": {
    "@id": "https://cran.r-project.org",
    "@type": "Organization",
    "name": "Comprehensive R Archive Network (CRAN)",
    "url": "https://cran.r-project.org"
  },
  "author": [
    {
      "@type": "Person",
      "givenName": "Carl",
      "familyName": "Boettiger",
      "email": "cboettig@gmail.com",
      "@id": "https://orcid.org/0000-0002-1642-628X"
    },
    {
      "@type": "Person",
      "givenName": "Maëlle",
      "familyName": "Salmon",
      "@id": "https://orcid.org/0000-0002-2815-0399"
    }
  ],
  "contributor": [
    {
      "@type": "Person",
      "givenName": "Anna",
      "familyName": "Krystalli",
      "@id": "https://orcid.org/0000-0002-2378-4915"
    },
    {
      "@type": "Person",
      "givenName": "Maëlle",
      "familyName": "Salmon",
      "@id": "https://orcid.org/0000-0002-2815-0399"
    },
    {
      "@type": "Person",
      "givenName": "Katrin",
      "familyName": "Leinweber",
      "@id": "https://orcid.org/0000-0001-5135-5758"
    },
    {
      "@type": "Person",
      "givenName": "Noam",
      "familyName": "Ross",
      "@id": "https://orcid.org/0000-0002-2136-0000"
    },
    {
      "@type": "Person",
      "givenName": "Arfon",
      "familyName": "Smith"
    },
    {
      "@type": "Person",
      "givenName": "Jeroen",
      "familyName": "Ooms",
      "@id": "https://orcid.org/0000-0002-4035-0289"
    },
    {
      "@type": "Person",
      "givenName": "Sebastian",
      "familyName": "Meyer",
      "@id": "https://orcid.org/0000-0002-1791-9449"
    },
    {
      "@type": "Person",
      "givenName": "Michael",
      "familyName": "Rustler",
      "@id": "https://orcid.org/0000-0003-0647-7726"
    },
    {
      "@type": "Person",
      "givenName": "Hauke",
      "familyName": "Sonnenberg",
      "@id": "https://orcid.org/0000-0001-9134-2871"
    },
    {
      "@type": "Person",
      "givenName": "Sebastian",
      "familyName": "Kreutzer",
      "@id": "https://orcid.org/0000-0002-0734-2199"
    },
    {
      "@type": "Person",
      "givenName": "Thierry",
      "familyName": "Onkelinx",
      "@id": "https://orcid.org/0000-0001-8804-4216"
    }
  ],
  "copyrightHolder": [
    {
      "@type": "Person",
      "givenName": "Carl",
      "familyName": "Boettiger",
      "email": "cboettig@gmail.com",
      "@id": "https://orcid.org/0000-0002-1642-628X"
    }
  ],
  "funder": [
    {
      "@type": "Organization",
      "name": "rOpenSci"
    }
  ],
  "maintainer": [
    {
      "@type": "Person",
      "givenName": "Carl",
      "familyName": "Boettiger",
      "email": "cboettig@gmail.com",
      "@id": "https://orcid.org/0000-0002-1642-628X"
    }
  ],
  "softwareSuggestions": [
    {
      "@type": "SoftwareApplication",
      "identifier": "withr",
      "name": "withr",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=withr"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "covr",
      "name": "covr",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=covr"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "details",
      "name": "details",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=details"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "dplyr",
      "name": "dplyr",
      "version": ">= 0.7.0",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=dplyr"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "jsonld",
      "name": "jsonld",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=jsonld"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "jsonvalidate",
      "name": "jsonvalidate",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=jsonvalidate"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "knitr",
      "name": "knitr",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=knitr"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "printr",
      "name": "printr",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=printr"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "rmarkdown",
      "name": "rmarkdown",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=rmarkdown"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "testthat",
      "name": "testthat",
      "version": ">= 3.0.0",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=testthat"
    },
    {
      "@type": "SoftwareApplication",
      "identifier": "usethis",
      "name": "usethis",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=usethis"
    }
  ],
  "softwareRequirements": {
    "1": {
      "@type": "SoftwareApplication",
      "identifier": "R",
      "name": "R",
      "version": ">= 3.2.0"
    },
    "2": {
      "@type": "SoftwareApplication",
      "identifier": "commonmark",
      "name": "commonmark",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=commonmark"
    },
    "3": {
      "@type": "SoftwareApplication",
      "identifier": "crul",
      "name": "crul",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=crul"
    },
    "4": {
      "@type": "SoftwareApplication",
      "identifier": "desc",
      "name": "desc",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=desc"
    },
    "5": {
      "@type": "SoftwareApplication",
      "identifier": "gert",
      "name": "gert",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=gert"
    },
    "6": {
      "@type": "SoftwareApplication",
      "identifier": "gh",
      "name": "gh",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=gh"
    },
    "7": {
      "@type": "SoftwareApplication",
      "identifier": "jsonlite",
      "name": "jsonlite",
      "version": ">= 1.6",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=jsonlite"
    },
    "8": {
      "@type": "SoftwareApplication",
      "identifier": "magrittr",
      "name": "magrittr",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=magrittr"
    },
    "9": {
      "@type": "SoftwareApplication",
      "identifier": "memoise",
      "name": "memoise",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=memoise"
    },
    "10": {
      "@type": "SoftwareApplication",
      "identifier": "methods",
      "name": "methods"
    },
    "11": {
      "@type": "SoftwareApplication",
      "identifier": "pingr",
      "name": "pingr",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=pingr"
    },
    "12": {
      "@type": "SoftwareApplication",
      "identifier": "purrr",
      "name": "purrr",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=purrr"
    },
    "13": {
      "@type": "SoftwareApplication",
      "identifier": "remotes",
      "name": "remotes",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=remotes"
    },
    "14": {
      "@type": "SoftwareApplication",
      "identifier": "sessioninfo",
      "name": "sessioninfo",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=sessioninfo"
    },
    "15": {
      "@type": "SoftwareApplication",
      "identifier": "stats",
      "name": "stats"
    },
    "16": {
      "@type": "SoftwareApplication",
      "identifier": "urltools",
      "name": "urltools",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=urltools"
    },
    "17": {
      "@type": "SoftwareApplication",
      "identifier": "xml2",
      "name": "xml2",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=xml2"
    },
    "18": {
      "@type": "SoftwareApplication",
      "identifier": "cli",
      "name": "cli",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=cli"
    },
    "19": {
      "@type": "SoftwareApplication",
      "identifier": "codemeta",
      "name": "codemeta",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=codemeta"
    },
    "SystemRequirements": null
  },
  "isPartOf": "https://ropensci.org",
  "keywords": ["metadata", "codemeta", "ropensci", "citation", "credit", "linked-data", "json-ld", "r", "rstats", "r-package", "peer-reviewed"],
  "fileSize": "NAKB",
  "releaseNotes": "https://github.com/ropensci/codemetar/blob/master/NEWS.md",
  "readme": "https://github.com/ropensci/codemetar/blob/master/README.md",
  "contIntegration": ["https://github.com/ropensci/codemetar/actions", "https://codecov.io/github/ropensci/codemetar?branch=master"],
  "developmentStatus": "https://www.repostatus.org/",
  "review": {
    "@type": "Review",
    "url": "https://github.com/ropensci/software-review/issues/130",
    "provider": "https://ropensci.org"
  }
}
```

</details>

<br>

By default most often from within your package folder you’ll simply run
`codemetar::write_codemeta()`.

## Keep codemeta.json up-to-date

**How to keep codemeta.json up-to-date?** In particular, how to keep it
up to date with `DESCRIPTION`? `codemetar` itself no longer supports
automatic sync, but there are quite a few methods available out there.
Choose one that fits well into your workflow!

-   You could rely on `devtools::release()` since it will ask you
    whether you updated codemeta.json when such a file exists.

-   You could use a git pre-commit hook that prevents a commit from
    being done if DESCRIPTION is newer than codemeta.json.

    -   You can use the [precommit
        package](https://github.com/lorenzwalthert/precommit) in which
        there’s a “codemeta-description-updated” hook.

    -   If that’s your only pre-commit hook (i.e. you don’t have one
        created by e.g. `usethis::use_readme_rmd()`), then you can
        create it using

``` r
script = readLines(system.file("templates", "description-codemetajson-pre-commit.sh", package = "codemetar"))
usethis::use_git_hook("pre-commit",
                     script = script)
```

-   You could use GitHub actions. Refer to GitHub actions docs
    <https://github.com/features/actions>, and to the example workflow
    provided in this package (type
    `system.file("templates", "codemeta-github-actions.yml", package = "codemetar")`).
    You can use the `cm-skip` keyword in your commit message if you
    don’t want this to run on a specific commit. The example workflow
    provided is setup to only run when a push is made to the master
    branch. This setup is designed for if you’re using a [git
    flow](https://nvie.com/posts/a-successful-git-branching-model/#the-main-branches)
    setup where the master branch is only committed and pushed to via
    pull requests. After each PR merge (and the completion of this
    GitHub action), your master branch will always be up to date and so
    long as you don’t make manual changes to the codemeta.json file, you
    won’t have merge conflicts.

Alternatively, you can have GitHub actions route run `codemetar` on each
commit. If you do this you should try to remember to run `git pull`
before making any new changes on your local project. However, if you
forgot to pull and already committed new changes, fret not, you can use
([`git pull --rebase`](https://stackoverflow.com/questions/18930527/difference-between-git-pull-and-git-pull-rebase/38139843#38139843))
to rewind you local changes on top of the current upstream `HEAD`.

<details closed>
<summary>
<span title="Click to Expand"> click here to see the workflow </span>
</summary>

``` yaml
on:
  push:
    branches: master
    paths:
      - DESCRIPTION
      - .github/workflows/main.yml

name: Render codemeta
jobs:
  render:
    name: Render codemeta
    runs-on: macOS-latest
    if: "!contains(github.event.head_commit.message, 'cm-skip')"
    steps:
      - uses: actions/checkout@v1
      - uses: r-lib/actions/setup-r@v1
      - name: Install codemetar
        run: Rscript -e 'install.packages("codemetar")'
      - name: Render codemeta
        run: Rscript -e 'codemetar::write_codemeta()'
      - name: Commit results
        run: |
          git commit codemeta.json -m 'Re-build codemeta.json' || echo "No changes to commit"
          git push https://${{github.actor}}:${{secrets.GITHUB_TOKEN}}@github.com/${{github.repository}}.git HEAD:${{ github.ref }} || echo "No changes to commit"
```

</details>

<br>

## How to improve your package’s codemeta.json?

The best way to ensure `codemeta.json` is as complete as possible is to
set metadata in all the usual places, and then if needed add more
metadata.

To ensure you have metadata in the usual places, you can run
`codemetar::give_opinions()`.

### Usual terms in DESCRIPTION

-   Fill `BugReports` and `URL`.

-   Using the `Authors@R` notation allows a much richer specification of
    author roles, correct parsing of given vs family names, and email
    addresses.

In the current implementation, developers may specify an ORCID url for
an author in the optional `comment` field of `Authors@R`, e.g.

    Authors@R: c(person(given = "Carl",
                 family = "Boettiger",
                 role = c("aut", "cre", "cph"),
                 email = "cboettig@gmail.com",
                 comment = c(ORCID = "0000-0002-1642-628X")))

which will allow `codemetar` to associate an identifier with the person.
This is clearly something of a hack since R’s `person` object lacks an
explicit notion of `id`, and may be frowned upon.

### Usual terms in the README

In the README, you can use badges for continuous integration, repo
development status (repostatus.org or lifecycle.org), provider
([e.g. for CRAN](https://docs.r-hub.io/#badges)).

### GitHub repo topics

If your package source is hosted on GitHub and there’s a way for
codemetar to determine that (URL in DESCRIPTION, or git remote URL)
codemetar will use [GitHub repo
topics](https://docs.github.com/en/github/administering-a-repository/classifying-your-repository-with-topics)
as keywords in codemeta.json. If you also set keywords in DESCRIPTION
(see next section), codemetar will merge the two lists.

### Set even more terms via DESCRIPTION

In general, setting metadata via the places stated earlier is the best
solution because that metadata is used by other tools (e.g. the URLs in
DESCRIPTION can help the package users, not only codemetar).

The DESCRIPTION file is the natural place to specify any metadata for an
R package. The `codemetar` package can detect certain additional terms
in the [CodeMeta context](https://codemeta.github.io/terms/). Almost any
additional codemeta field can be added to and read from the DESCRIPTION
into a `codemeta.json` file (see `codemetar:::additional_codemeta_terms`
for a list).

CRAN requires that you prefix any additional such terms to indicate the
use of `schema.org` explicitly, e.g. `keywords` would be specified in a
DESCRIPTION file as:

    X-schema.org-keywords: metadata, codemeta, ropensci, citation, credit, linked-data

Where applicable, these will override values otherwise guessed from the
source repository. Use comma-separated lists to separate multiple values
to a property, e.g. keywords.

See the
[DESCRIPTION](https://github.com/ropensci/codemetar/blob/master/DESCRIPTION)
file of the `codemetar` package for an example.

### Set the branch that codemetar references

There are a number of places that codemetar will reference a github
branch if your code is hosted on github (e.g. for release notes, readme,
etc.). By default, codemetar will use the name “master” but you can
change that to whatever your default branch is by setting the option
“codemeta_branch” (e.g. `options(codemeta_branch = "main")` before
calling `write_codemeta()` to use the branch named “main” as the default
branch).

## Installation and usage requirements

You can install the latest version from CRAN using:

``` r
install.packages("codemetar")
```

You can also install the development version of `codemetar` from GitHub
with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/codemetar")
```

For optimal results you need a good internet connection.

The package queries

-   `utils::available.packages()` for CRAN and Bioconductor packages;

-   GitHub API via the [`gh` package](https://github.com/r-lib/gh), if
    it finds a GitHub repo URL in DESCRIPTION or as git remote. GitHub
    API is queried to find the [preferred
    README](https://developer.github.com/v3/repos/contents/#get-the-readme),
    and the [repo
    topics](https://developer.github.com/v3/repos/#list-all-topics-for-a-repository).
    If you use codemetar for many packages having a
    [GITHUB_PAT](https://github.com/r-lib/gh#environment-variables) is
    better;

-   [R-hub sysreqs API](https://docs.r-hub.io/#sysreqs) to parse
    SystemRequirements.

If your machine is offline, a more minimal codemeta.json will be
created. If your internet connection is poor or there are firewalls, the
codemeta creation might indefinitely hang.

## Going further

Check out all the [codemetar
man](https://docs.ropensci.org/codemetar/articles/index.html) for
tutorials on other cool stuff you can do with codemeta and json-ld.

A new feature is the creation of a minimal schemaorg.json for insertion
on your website’s webpage for Search Engine Optimization, when the
`write_minimeta` argument of `write_codemeta()` is `TRUE`.

You could e.g. use the code below in a chunk in README.Rmd with
`results="asis"`.

``` r
glue::glue('<script type="application/ld+json">
      {glue::glue_collapse(readLines("schemaorg.json"), sep = "\n")}
    </script>')
```

Refer to [Google
documentation](https://developers.google.com/search/docs) for more
guidance.

<script type="application/ld+json">
      {
  "@context": "https://schema.org",
  "type": "SoftwareSourceCode",
  "author": [
    {
      "id": "https://orcid.org/0000-0002-2815-0399"
    },
    {
      "id": "https://orcid.org/0000-0002-1642-628X"
    }
  ],
  "codeRepository": "https://github.com/ropensci/codemetar",
  "contributor": [
    {
      "id": "https://orcid.org/0000-0002-2378-4915",
      "type": "Person",
      "familyName": "Krystalli",
      "givenName": "Anna"
    },
    {
      "id": "https://orcid.org/0000-0002-2815-0399",
      "type": "Person",
      "familyName": "Salmon",
      "givenName": "Maëlle"
    },
    {
      "id": "https://orcid.org/0000-0001-5135-5758",
      "type": "Person",
      "familyName": "Leinweber",
      "givenName": "Katrin"
    },
    {
      "id": "https://orcid.org/0000-0002-2136-0000",
      "type": "Person",
      "familyName": "Ross",
      "givenName": "Noam"
    },
    {
      "type": "Person",
      "familyName": "Smith",
      "givenName": "Arfon"
    },
    {
      "id": "https://orcid.org/0000-0002-4035-0289",
      "type": "Person",
      "familyName": "Ooms",
      "givenName": "Jeroen"
    },
    {
      "id": "https://orcid.org/0000-0002-1791-9449",
      "type": "Person",
      "familyName": "Meyer",
      "givenName": "Sebastian"
    },
    {
      "id": "https://orcid.org/0000-0003-0647-7726",
      "type": "Person",
      "familyName": "Rustler",
      "givenName": "Michael"
    },
    {
      "id": "https://orcid.org/0000-0001-9134-2871",
      "type": "Person",
      "familyName": "Sonnenberg",
      "givenName": "Hauke"
    },
    {
      "id": "https://orcid.org/0000-0002-0734-2199",
      "type": "Person",
      "familyName": "Kreutzer",
      "givenName": "Sebastian"
    }
  ],
  "copyrightHolder": {
    "id": "https://orcid.org/0000-0002-1642-628X",
    "type": "Person",
    "email": "cboettig@gmail.com",
    "familyName": "Boettiger",
    "givenName": "Carl"
  },
  "description": "The 'Codemeta' Project defines a 'JSON-LD' format\n    for describing software metadata, as detailed at\n    <https://codemeta.github.io>. This package provides utilities to\n    generate, parse, and modify 'codemeta.json' files automatically for R\n    packages, as well as tools and examples for working with\n    'codemeta.json' 'JSON-LD' more generally.",
  "funder": {
    "type": "Organization",
    "name": "rOpenSci"
  },
  "license": "https://spdx.org/licenses/GPL-3.0",
  "name": "codemetar: Generate 'CodeMeta' Metadata for R Packages",
  "programmingLanguage": {
    "type": "ComputerLanguage",
    "name": "R",
    "url": "https://r-project.org"
  },
  "provider": {
    "id": "https://cran.r-project.org",
    "type": "Organization",
    "name": "Comprehensive R Archive Network (CRAN)",
    "url": "https://cran.r-project.org"
  },
  "runtimePlatform": "R version 3.6.1 (2019-07-05)",
  "version": "0.1.8.9000"
}
    </script>

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# codemetar (development version)

# codmetar 0.3.3

* Now uses the existing `codemeta` package on CRAN for all core functionality.
* This package provides more advanced metadata detection than the minimalist `codemeta`
* Return to CRAN


# codemetar 0.3.2

* Skip on CRAN a test relying on internet connection.
* Bug fix: Fix regression in Windows path handling (#315, @ms609).

# codemeta 0.3.1

* CRAN Solaris issues

# codemetar 0.3.0

* `write_codemetar()` can now be called from anywhere within a package directory structure. (#305, @mpadge)
* Breaking change: `write_codemeta()` writes the JSON file at `path` relative to `pkg`, not the current directory. (#303, @ThierryO)
* Added documentation for changing the default branch (#302, @jonkeane)
* Breaking change: Relatedly, it is no longer possible to use `write_codemeta()` on an installed packages, in that case one would have to use `create_codemeta()` together with `jsonlite::write_json()`.
* Use R.Version() instead of R.version to allow mocking in tests.
* Bug fix: now able to parse a README where badges are in a table with non badges links.
* Bug fix: `guess_fileSize()` properly handles `.Rbuildignore` (#299, @ThierryO).
* Bug fix: `create_codemetar()` handles minimal packages (#298, @ThierryO).

# codemetar 0.1.9 2020-07-16

## Deprecation

* The use_git_hook argument of write_codemeta() has been deprecated. Solutions for keeping DESCRIPTION and codemeta.json in sync are available in the docs.

* drops `crosswalk`, [#288]

## Enhancements

* Docs were improved to make a better case for codemetar.

* Changes in the way codeRepository is guessed. codemetar can now recognize an URL from GitHub, GitLab, Bitbucket, R-Forge among several URLs in DESCRIPTION, to assign it to codeRepository. If no URL in DESCRIPTION is from any of these providers, `guess_github()` is called.

* Adds documentation of internet needs and verbosity to steps downloading information from the web (#270, @Bisaloo)

* New argument `write_minimeta` for `write_codemeta()` indicating whether to also create the file schemaorg.json that  corresponds to the metadata Google would validate, to be inserted to a webpage for SEO. It is saved as "schemaorg.json" alongside `path` (by default, "codemeta.json"). This functionality requires the `jsonld` package (listed under `Suggests`).

* Updated the GitHub action template to only run on pushes to the master branch and added an explanation of how that works to the readme. (@jonkeane)

## Bug fixes

* Fix for detecting rOpenSci review badge (@sckott, #236)

* Fix extraction of ORCID when composite comment (@billy34, #231)

* Fix bug in crosswalking (#243)

* Bug fix: the codeRepository is updated if there's any URL in DESCRIPTION.

* Bug fix: the README information is now updated by codemeta_readme(). Previously if e.g. a developmentStatus had been set previously, it was never updated.

## Internals

* Code cleaning following the book Martin, Robert C. Clean code: a handbook of agile software craftsmanship. Pearson Education, 2009. (@hsonne, #201, #202, #204, #205, #206, #207, #209, #210, #211, #212, #216, #218, #219, #220, #221).

* Use of re-usable Rmd pieces for the README, intro vignette and man pages to reduce copy-pasting.

# codemetar 0.1.8 2019-05

* address internet timeout issues
* tidy source code
* update test suite to reflect newly available metadata.
* `write_codemeta()` and `create_codemeta()`: `use_filesize = FALSE` is now the default and the estimation of the file size does not leave any more unwanted files behind [PR #239](https://github.com/ropensci/codemetar/pull/239). Furthermore, the way the file size is calculated changed: Before we used the size of the package built with `pkgbuild::build()`, which took rather long. Now the size is calculated based on the source files minus files excluded via  
 `.Rbuildignore` (if such a file exists).
* `write_codemeta()`: the default of argument `use_git_hook` is now `FALSE` to avoid an 
unwanted alteration of the user's git environment [issue #240](https://github.com/ropensci/codemetar/issues/240).
* Package dependency to 'pkgbuild' has been dropped.
* `write_codemeta()` does not crash anymore if the `CITATION` file contains a line `citation(auto = meta)` [Issue #238](https://github.com/ropensci/codemetar/issues/238).


# codemetar 0.1.7 2018-12

* `jsonld` that is used only in `crosswalk()` and `codemeta_validate()` is now an optional dependency (Suggests rather than Imports).

* The CodeRepository URL is now cleaned a bit (removing direct link to the README).

* `write_codemeta()` gains a new argument `use_git_hook` to make the creation of a DESCRIPTION/codemeta.json git pre-commit hook optional.

* `create_codemeta()` and `write_codemeta()` gain a new argument `use_filesize` to make the building of the package to get its size optional.

* Encoding bug fixed in `extract_badges()`.

* `extract_badges()` now uses `commonmark` and `xml2` instead of only regular expressions.

* `pkgbuild` is now used directly instead of `devtools`.

* `give_opinion()` now recognizes lifecycle badges, not only repostatus.org badges.

# codemetar 0.1.6 2018-04

* New functions

    * extract_badges for extracting information from all badges in a Markdown file.
    
    * give_opinion giving opiniated advice about package metadata
    
* Changes to the create_codemeta output

    * relatedLink field now include provider URL and URL(s) from DESCRIPTION that are not the code repository
    
    * maintainer is now a list allowing for several maintainers since e.g. the BioConductor a4 package has two maintainers.
    
    * if more than one CI service among Travis, Appveyor and Circle CI are used and shown via a README badge they're all added to the contIntegration field. URLs from codecov and coveralls badges are also added to the contIntegration field.
    
    * repo status inferred from the README now 1) is an URL instead of a word 2) recognizes either repostatus.org or Tidyverse lifecycle badges.
    
    * if present, priority is given to the Repository and BugReports fields of DESCRIPTION for filling the codeRepository and issueTracker fields of codemeta.json (which means working on a fork won't change these).
    
    * ability to parse all CRAN-allowed MARC roles.
    
    * if there is a badge for an rOpenSci onboarding review and the review issue is closed, basic review metadata is added to codemeta.json
    
    * For dependencies, if the provider guessed is CRAN or BioConductor, their canonic CRAN/BioConductor URL is added to codemeta.json as sameAs, unless there's a GitHub repo mentioned for them in Remotes in DESCRIPTION, in which case sameAs is that GitHub repo.
    
    * CRAN is now correctly translated as "Comprehensive R Archive Network"
    
    * If codeRepository is guessed to be a GitHub repo (via the URL field of DESCRIPTION or via git remote URL), the repo topics are queried via GitHub API V3 and added to the keywords (in combination with keywords stored in the X-schema.org-keywords field of DESCRIPTION)
    
    * SystemRequirements are now parsed using https://sysreqs.r-hub.io/, outputting URLs then stored in softwareRequirements

* Help to remind to update codemeta.json regularly: Writing codemeta.json for the first time adds a git pre-commit hook and suggests adding a release question for devtools::release.

* Internal changes

    * Now uses desc to parse DESCRIPTION files.

    * Package license changed to GPL because of code borrowed from usethis
    
    * Uses crul instead of httr and uses crul to check some URLs.
    
    * write_codemeta only uses Rbuildignore and a pre-commit git hook if the function is called from a package folder directly and with the path argument equal to "codemeta.json"
    
    * The calls to available.packages() for guess_provider now happen inside memoised functions.
    
    * codemeta_readme function.

# codemetar 0.1.5 2018-03-21

* Default to DOI-based schema. (previous CN issues now resolved)

# codemetar 0.1.4 2018-02-12

* Allow vignettes to gracefully handle network timeout errors that
  may occur on CRAN's Windows build server.

# codemetar 0.1.3 2018-02-08

* CRAN release
* Switch to <https://purl.org/> based URIs for the JSON-LD 
  Context file instead of a DOI, due to frequent failure
  of content negotiation on DataCite servers
  ([#34](https://github.com/ropensci/codemetar/issues/34))
* bugfix UTF-8 characters in CITATION files 
  ([#44](https://github.com/ropensci/codemetar/issues/44))
* bugfix to git URLs
* Use `https` on ORCID `@id` URIs

# codemetar 0.1.2

* JOSS release

# codemetar 0.1.1

* Post onboarding release

# codemetar 0.1.0

* Added a `NEWS.md` file to track changes to the package.





Contributing Guidelines
=======================

Repository structure
--------------------

This repository is structured as a standard R package
following the conventions outlined in the [Writing R
extensions](http://cran.r-project.org/doc/manuals/R-exts.html) manual.
A few additional files are provided that are not part of the built
R package and are listed in `.Rbuildignore`, such as `.travis.yml`,
which is used for continuous testing and integration.


Code
----

All code for this package is found in `R/`, (except compiled source
code, if used, which is in `/src`).  All functions should be thoroughly
documented with `roxygen2` notation; see Documentation. Code should
conform to our [Style guide](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md)

Testing
-------

Any new feature or bug-fix should include a unit-test demonstrating the
change.  Unit tests follow the `testthat` framework with files in
`tests/testthat`.  Please make sure that the testing suite passes
before issuing a pull request.  This can be done by running `check()`
from the `devtools` package, which will also check for consistent
documentation, etc.


This package uses the [travis](https://github.com/craigcitro/r-travis)
continuous testing mechanism for R to ensure that the test suite is run
on each push to Github.  An icon at the top of the README.md indicates
whether or not the tests are currently passing. 

This package also uses
[codecov.io](https://codecov.io/) to
measure test coverage.  While not all code can be covered by automated 
tests (in particular, functions involving user prompts), try to avoid
decreasing coverage by writing unit tests for any contributed code. 
Codecov.io will flag PRs that decrease coverage. 


Documentation
-------------

All of the function documentation is generated automatically.
Please do not edit any of the documentation files in `man/`
or the `NAMESPACE`.  Instead, construct the appropriate
[roxygen2](https://github.com/klutometis/roxygen) documentation in the
function files in `R/` themselves.  The documentation is then generated
by running the `document()` function from the `devtools` package.  Please
consult the [Advanced R programming](http://adv-r.had.co.nz/) guide if
this workflow is unfamiliar to you.  Note that functions should include
examples in the documentation. Please use `\dontrun` for examples that
take more than a few seconds to execute or require an internet connection.

Likewise, the README.md file in the base directory should not be edited
directly.  This file is created automatically from code that runs the
examples shown, helping to ensure that they are functioning as advertised
and consistent with the package README vignette.  Instead, edit the
`README.Rmd` source file in `manuscripts` and run `make` to build
the README.

General Development Goals & Guidelines
---------------------------------------

1. Not having too many high-level functions, 
2. Using sensible defaults, (driven by use cases).
3. Docs should point advanced users to the lower-level API when they need special cases.  
4. Maintain a consistent user-facing API.   
Dear CRAN team,

Bugfix as described in NEWS.md
Resolves previous CRAN check problems.


## Test environments

* local Ubuntu install
* GitHub Actions (macOS-latest, R release; Windows latest, R release, Windows latest, Ubuntu R-devel) 
* win-builder
* R-hub

## R CMD check results

0 errors | 0 warnings | 1 note

```
New submission
  
  Package was archived on CRAN
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2021-05-15 as check problems were not
      corrected in time.
```
Parsing codemeta data
================


**NB: sources moved to vignettes**

``` r
library(jsonld)
library(jsonlite)
library(magrittr)
library(codemetar)
library(tidyverse)
library(printr)
```

``` r
write_codemeta("codemetar", "codemeta.json")
```

Digest input with a frame:

``` r
frame <- system.file("schema/frame_schema.json", package="codemetar")

meta <- 
  jsonld_frame("codemeta.json", frame) %>%
  fromJSON(FALSE) %>% getElement("@graph") %>% getElement(1)
```

Construct a citation

``` r
authors <- 
lapply(meta$author, 
       function(author) 
         person(given = author$given, 
                family = author$family, 
                email = author$email,
                role = "aut"))
year <- meta$datePublished
if(is.null(year)) 
  year <- format(Sys.Date(), "%Y")
bibitem <- 
 bibentry(
     bibtype = "Manual",
     title = meta$name,
     author = authors,
     year = year,
     note = paste0("R package version ", meta$version),
     url = meta$URL,
     key = meta$identifier
   )

cat(format(bibitem, "bibtex"))
```

    ## @Manual{codemetar,
    ##   title = {codemetar: Generate CodeMeta Metadata for R Packages},
    ##   author = {Carl Boettiger},
    ##   year = {2017},
    ##   note = {R package version 0.1.0},
    ## }

``` r
bibitem
```

    ## Boettiger C (2017). _codemetar: Generate CodeMeta Metadata for R
    ## Packages_. R package version 0.1.0.

Parsing the ropensci corpus
---------------------------

Frame, expanding any referenced nodes

``` r
corpus <- 
    jsonld_frame("ropensci.json", frame) %>%
    fromJSON(simplifyVector = FALSE) %>%
    getElement("@graph") 
```

Some basics:

``` r
## deal with nulls explicitly by starting with map
pkgs <- map(corpus, "name") %>% compact() %>% as.character()

# keep only those with package identifiers (names)
keep <- map_lgl(corpus, ~ length(.x$identifier) > 0)
corpus <- corpus[keep]

## now we can just do
all_pkgs <- map_chr(corpus, "name")
head(all_pkgs)
```

    ## [1] "AntWeb: programmatic interface to the AntWeb"                                
    ## [2] "aRxiv: Interface to the arXiv API"                                           
    ## [3] "chromer: Interface to Chromosome Counts Database API"                        
    ## [4] "ckanr: Client for the Comprehensive Knowledge Archive Network ('CKAN') 'API'"
    ## [5] "dashboard: A package status dashboard"                                       
    ## [6] "ggit: Git Graphics"

``` r
## 60 unique maintainers
map_chr(corpus, c("maintainer", "familyName")) %>% unique() %>% length()
```

    ## [1] 61

``` r
## Mostly Scott
map_chr(corpus, c("maintainer", "familyName")) %>% 
  as_tibble() %>%
  group_by(value) %>%
  tally(sort=TRUE)
```

| value        |    n|
|:-------------|----:|
| Chamberlain  |  105|
| Ooms         |   12|
| Mullen       |    8|
| Ram          |    8|
| Boettiger    |    6|
| Salmon       |    5|
| FitzJohn     |    4|
| Hart         |    2|
| Leeper       |    2|
| Marwick      |    2|
| Müller       |    2|
| Padgham      |    2|
| South        |    2|
| Varela       |    2|
| Vitolo       |    2|
| Arnold       |    1|
| Attali       |    1|
| Banbury      |    1|
| Becker       |    1|
| Bengtsson    |    1|
| Braginsky    |    1|
| Broman       |    1|
| Bryan        |    1|
| Dallas       |    1|
| de Queiroz   |    1|
| Drost        |    1|
| Fischetti    |    1|
| Ghahraman    |    1|
| Goring       |    1|
| hackathoners |    1|
| Harrison     |    1|
| Hughes       |    1|
| Jahn         |    1|
| Jones        |    1|
| Keyes        |    1|
| Krah         |    1|
| Lehtomaki    |    1|
| Lovelace     |    1|
| Lundstrom    |    1|
| McGlinn      |    1|
| McVey        |    1|
| Meissner     |    1|
| Michonneau   |    1|
| Moroz        |    1|
| Otegui       |    1|
| Pardo        |    1|
| Pennell      |    1|
| Poelen       |    1|
| Robinson     |    1|
| Ross         |    1|
| Rowlingson   |    1|
| Scott        |    1|
| Seers        |    1|
| Shotwell     |    1|
| Sievert      |    1|
| Sparks       |    1|
| Stachelek    |    1|
| Szöcs        |    1|
| Widgren      |    1|
| Wiggin       |    1|
| Winter       |    1|

``` r
## number of co-authors ... 
map_int(corpus, function(r) length(r$author)) %>% 
  as_tibble() %>%
  group_by(value) %>%
  tally(sort=TRUE)
```

|  value|    n|
|------:|----:|
|      1|  146|
|      2|   30|
|      3|   17|
|      4|    8|
|      5|    5|
|      7|    3|
|     13|    1|

``` r
## Contributors isn't used as much...
map_int(corpus, function(r) length(r$contributor)) %>% 
  as_tibble() %>%
  group_by(value) %>%
  tally(sort=TRUE)
```

|  value|    n|
|------:|----:|
|      0|  178|
|      2|   13|
|      4|    9|
|      3|    7|
|      5|    1|
|      6|    1|
|      8|    1|

Numbers (n) of packages with a total of (value) dependencies:

``` r
map_int(corpus, function(r) length(r$softwareRequirements))  %>% 
  as_tibble() %>%
  group_by(value) %>%
  tally(sort=TRUE)
```

|  value|    n|
|------:|----:|
|      4|   39|
|      5|   35|
|      2|   25|
|      3|   25|
|      7|   19|
|      6|   16|
|      8|   13|
|      9|    8|
|     12|    7|
|     10|    6|
|     11|    6|
|     13|    3|
|      0|    2|
|     14|    1|
|     17|    1|
|     18|    1|
|     21|    1|
|     22|    1|
|     23|    1|

which dependencies are used most frequently?

``` r
corpus %>%
map_df(function(x){
  ## single, unboxed dep
  if("name" %in% names(x$softwareRequirements))
    dep <- x$name
  else if("name" %in% names(x$softwareRequirements[[1]]))
    dep <- map_chr(x$softwareRequirements, "name")
  else { ## No requirementsß
    dep <- NA
  }
  
  tibble(identifier = x$identifier, dep = dep)
}) -> dep_df


dep_df %>%
group_by(dep) %>% 
  tally(sort = TRUE)
```

| dep                                                              |    n|
|:-----------------------------------------------------------------|----:|
| jsonlite                                                         |   99|
| httr                                                             |   92|
| R                                                                |   66|
| tibble                                                           |   46|
| dplyr                                                            |   43|
| methods                                                          |   37|
| xml2                                                             |   37|
| data.table                                                       |   35|
| utils                                                            |   35|
| crul                                                             |   31|
| plyr                                                             |   29|
| XML                                                              |   25|
| magrittr                                                         |   24|
| sp                                                               |   22|
| stringr                                                          |   21|
| curl                                                             |   18|
| ggplot2                                                          |   18|
| lazyeval                                                         |   17|
| stats                                                            |   17|
| lubridate                                                        |   14|
| R6                                                               |   14|
| rappdirs                                                         |   13|
| assertthat                                                       |   12|
| digest                                                           |   12|
| RCurl                                                            |   12|
| readr                                                            |   11|
| rgdal                                                            |   10|
| whisker                                                          |   10|
| scales                                                           |    9|
| ape                                                              |    8|
| raster                                                           |    8|
| tidyr                                                            |    8|
| Rcpp                                                             |    7|
| reshape2                                                         |    7|
| rvest                                                            |    7|
| rgeos                                                            |    6|
| V8                                                               |    6|
| hoardr                                                           |    5|
| rjson                                                            |    5|
| taxize                                                           |    5|
| tools                                                            |    5|
| git2r                                                            |    4|
| maps                                                             |    4|
| oai                                                              |    4|
| openssl                                                          |    4|
| R(&gt;=3.2.1)                                                    |    4|
| solrium                                                          |    4|
| urltools                                                         |    4|
| foreach                                                          |    3|
| knitr                                                            |    3|
| leaflet                                                          |    3|
| maptools                                                         |    3|
| memoise                                                          |    3|
| mime                                                             |    3|
| pdftools                                                         |    3|
| purrr                                                            |    3|
| RColorBrewer                                                     |    3|
| rgbif                                                            |    3|
| rmarkdown                                                        |    3|
| shiny                                                            |    3|
| spocc                                                            |    3|
| stringi                                                          |    3|
| uuid                                                             |    3|
| wicket                                                           |    3|
| yaml                                                             |    3|
| base64enc                                                        |    2|
| bibtex                                                           |    2|
| Biostrings                                                       |    2|
| crayon                                                           |    2|
| devtools                                                         |    2|
| downloader                                                       |    2|
| fauxpas                                                          |    2|
| gdata                                                            |    2|
| gistr                                                            |    2|
| graphics                                                         |    2|
| grid                                                             |    2|
| htmltools                                                        |    2|
| htmlwidgets                                                      |    2|
| httpcode                                                         |    2|
| igraph                                                           |    2|
| jqr                                                              |    2|
| MASS                                                             |    2|
| miniUI                                                           |    2|
| ncdf4                                                            |    2|
| png                                                              |    2|
| R.cache                                                          |    2|
| R.utils                                                          |    2|
| rcrossref                                                        |    2|
| rentrez                                                          |    2|
| reshape                                                          |    2|
| rmapshaper                                                       |    2|
| rplos                                                            |    2|
| rvertnet                                                         |    2|
| shinyjs                                                          |    2|
| storr                                                            |    2|
| tm                                                               |    2|
| NA                                                               |    2|
| analogue                                                         |    1|
| antiword: Extract Text from Microsoft Word Documents             |    1|
| apipkgen: Package Generator for HTTP API Wrapper Packages        |    1|
| appl: Approximate POMDP Planning Software                        |    1|
| aRxiv                                                            |    1|
| binman                                                           |    1|
| Biobase                                                          |    1|
| BiocGenerics                                                     |    1|
| biomaRt                                                          |    1|
| bold                                                             |    1|
| caTools                                                          |    1|
| ckanr                                                            |    1|
| cld2: Google's Compact Language Detector 2                       |    1|
| countrycode                                                      |    1|
| cranlogs                                                         |    1|
| crminer                                                          |    1|
| crosstalk                                                        |    1|
| DBI                                                              |    1|
| dirdf: Extracts Metadata from Directory and File Names           |    1|
| doParallel                                                       |    1|
| DT(&gt;=0.1)                                                     |    1|
| elastic                                                          |    1|
| EML                                                              |    1|
| fastmatch                                                        |    1|
| foreign                                                          |    1|
| functionMap                                                      |    1|
| genderdata: Historical Datasets for Predicting Gender from Names |    1|
| GenomeInfoDb                                                     |    1|
| GenomicFeatures                                                  |    1|
| GenomicRanges(&gt;=1.23.24)                                      |    1|
| geoaxe                                                           |    1|
| geojson                                                          |    1|
| geojsonrewind: Fix 'GeoJSON' Winding Direction                   |    1|
| geonames                                                         |    1|
| geoops: 'GeoJSON' Manipulation Operations                        |    1|
| geosphere                                                        |    1|
| getPass                                                          |    1|
| ggm                                                              |    1|
| ggmap                                                            |    1|
| ggthemes                                                         |    1|
| graphql                                                          |    1|
| grDevices                                                        |    1|
| gridExtra                                                        |    1|
| gtools                                                           |    1|
| hash                                                             |    1|
| hexbin                                                           |    1|
| historydata: Data Sets for Historians                            |    1|
| Hmisc                                                            |    1|
| httpuv                                                           |    1|
| IRanges                                                          |    1|
| isdparser                                                        |    1|
| jsonvalidate                                                     |    1|
| jsonvalidate: Validate 'JSON'                                    |    1|
| leafletR                                                         |    1|
| loggr                                                            |    1|
| mapproj                                                          |    1|
| markdown                                                         |    1|
| Matrix                                                           |    1|
| memisc                                                           |    1|
| miniUI(&gt;=0.1.1)                                               |    1|
| nabor                                                            |    1|
| natserv                                                          |    1|
| openxlsx                                                         |    1|
| osmar                                                            |    1|
| outliers                                                         |    1|
| pdftools: Text Extraction and Rendering of PDF Documents         |    1|
| phytools                                                         |    1|
| plotly                                                           |    1|
| plumber                                                          |    1|
| progress                                                         |    1|
| protolite                                                        |    1|
| qlcMatrix                                                        |    1|
| RApiSerialize                                                    |    1|
| rapport                                                          |    1|
| rbhl                                                             |    1|
| rbison                                                           |    1|
| rebird                                                           |    1|
| redland                                                          |    1|
| redux                                                            |    1|
| remotes                                                          |    1|
| ridigbio                                                         |    1|
| ritis                                                            |    1|
| rJava                                                            |    1|
| RJSONIO                                                          |    1|
| rlist                                                            |    1|
| Rmpfr                                                            |    1|
| RMySQL                                                           |    1|
| rncl                                                             |    1|
| rnoaa                                                            |    1|
| rnrfa                                                            |    1|
| rotl                                                             |    1|
| rowr                                                             |    1|
| RPostgreSQL                                                      |    1|
| rredis                                                           |    1|
| rredlist                                                         |    1|
| RSQLite                                                          |    1|
| rstudioapi(&gt;=0.5)                                             |    1|
| rtracklayer                                                      |    1|
| rworldmap                                                        |    1|
| rzmq: R Bindings for ZeroMQ                                      |    1|
| S4Vectors                                                        |    1|
| scrapeR                                                          |    1|
| selectr                                                          |    1|
| sf                                                               |    1|
| shiny(&gt;=0.13.2)                                               |    1|
| snow                                                             |    1|
| SnowballC                                                        |    1|
| spatstat                                                         |    1|
| SSOAP                                                            |    1|
| stringdist                                                       |    1|
| sys                                                              |    1|
| tabulizerjars                                                    |    1|
| testthat                                                         |    1|
| tif: Text Interchange Format                                     |    1|
| USAboundariesData: Datasets for the 'USAboundaries' package      |    1|
| VariantAnnotation                                                |    1|
| viridisLite                                                      |    1|
| wdman(&gt;=0.2.2)                                                |    1|
| wellknown                                                        |    1|
| wicket: Utilities to Handle WKT Spatial Data                     |    1|
| WikidataR                                                        |    1|
| wikitaxa                                                         |    1|
| withr                                                            |    1|
| worrms                                                           |    1|
| xslt: XSLT 1.0 Transformations                                   |    1|
| zoo                                                              |    1|

Alternate approach using a frame, gets all Depends and suggests (really all `SoftwareApplication` types mentioned)

``` r
dep_frame <- '{
  "@context": "https://raw.githubusercontent.com/codemeta/codemeta/master/codemeta.jsonld",
  "@explicit": "true",
  "name": {}
}'
jsonld_frame("ropensci.json", dep_frame) %>% 
  fromJSON() %>% 
  getElement("@graph") %>%
  filter(type == "SoftwareApplication") %>%
  group_by(name) %>% 
  tally(sort = TRUE)
```

| name                        |    n|
|:----------------------------|----:|
| testthat                    |  168|
| knitr                       |  122|
| jsonlite                    |  105|
| httr                        |   96|
| roxygen2                    |   92|
| R                           |   72|
| rmarkdown                   |   68|
| covr                        |   52|
| dplyr                       |   49|
| tibble                      |   48|
| xml2                        |   41|
| methods                     |   38|
| utils                       |   37|
| data.table                  |   36|
| ggplot2                     |   36|
| crul                        |   33|
| plyr                        |   32|
| magrittr                    |   28|
| sp                          |   26|
| XML                         |   25|
| curl                        |   21|
| stringr                     |   21|
| lazyeval                    |   18|
| stats                       |   18|
| lubridate                   |   16|
| R6                          |   14|
| readr                       |   14|
| rgdal                       |   14|
| rappdirs                    |   13|
| assertthat                  |   12|
| devtools                    |   12|
| digest                      |   12|
| raster                      |   12|
| RCurl                       |   12|
| scales                      |   12|
| Rcpp                        |   11|
| whisker                     |   11|
| leaflet                     |   10|
| rgeos                       |   10|
| taxize                      |   10|
| tidyr                       |   10|
| reshape2                    |    9|
| ape                         |    8|
| maps                        |    8|
| V8                          |    8|
| maptools                    |    7|
| purrr                       |    7|
| rvest                       |    7|
| pdftools                    |    6|
| rgbif                       |    6|
| shiny                       |    6|
| ggmap                       |    5|
| git2r                       |    5|
| hoardr                      |    5|
| ncdf4                       |    5|
| png                         |    5|
| rjson                       |    5|
| tools                       |    5|
| oai                         |    4|
| openssl                     |    4|
| R(&gt;=3.2.1)               |    4|
| rcrossref                   |    4|
| RSQLite                     |    4|
| sf                          |    4|
| solrium                     |    4|
| urltools                    |    4|
| uuid                        |    4|
| yaml                        |    4|
| DBI                         |    3|
| fauxpas                     |    3|
| foreach                     |    3|
| gdata                       |    3|
| gistr                       |    3|
| graphics                    |    3|
| lintr                       |    3|
| MASS                        |    3|
| memoise                     |    3|
| mime                        |    3|
| miniUI                      |    3|
| R.utils                     |    3|
| RColorBrewer                |    3|
| rentrez                     |    3|
| rmapshaper                  |    3|
| rvertnet                    |    3|
| rworldmap                   |    3|
| spocc                       |    3|
| stringi                     |    3|
| wicket                      |    3|
| base64enc                   |    2|
| bibtex                      |    2|
| Biostrings                  |    2|
| broom                       |    2|
| crayon                      |    2|
| downloader                  |    2|
| elastic                     |    2|
| geiger                      |    2|
| getPass                     |    2|
| GGally                      |    2|
| ggthemes                    |    2|
| grDevices                   |    2|
| grid                        |    2|
| gridExtra                   |    2|
| htmltools                   |    2|
| htmlwidgets                 |    2|
| httpcode                    |    2|
| igraph                      |    2|
| jqr                         |    2|
| jsonvalidate                |    2|
| listviewer                  |    2|
| mapproj                     |    2|
| Matrix                      |    2|
| phylobase                   |    2|
| phytools                    |    2|
| R.cache                     |    2|
| RcppRedis                   |    2|
| readxl                      |    2|
| remotes                     |    2|
| reshape                     |    2|
| rplos                       |    2|
| shinyjs                     |    2|
| storr                       |    2|
| sys                         |    2|
| tm                          |    2|
| viridis                     |    2|
| webp                        |    2|
| zoo                         |    2|
| akima                       |    1|
| analogue                    |    1|
| aRxiv                       |    1|
| binman                      |    1|
| Biobase                     |    1|
| BiocGenerics                |    1|
| biomaRt                     |    1|
| bold                        |    1|
| Cairo                       |    1|
| caTools                     |    1|
| ckanr                       |    1|
| corrplot                    |    1|
| countrycode                 |    1|
| cranlogs                    |    1|
| crminer                     |    1|
| crosstalk                   |    1|
| dendextend                  |    1|
| doParallel                  |    1|
| dplyr(&gt;=0.3.0.2)         |    1|
| DT(&gt;=0.1)                |    1|
| EML                         |    1|
| etseed                      |    1|
| fastmatch                   |    1|
| fields                      |    1|
| forecast                    |    1|
| foreign                     |    1|
| fulltext                    |    1|
| functionMap                 |    1|
| genderdata                  |    1|
| GenomeInfoDb                |    1|
| GenomicFeatures             |    1|
| GenomicRanges(&gt;=1.23.24) |    1|
| geoaxe                      |    1|
| geojson                     |    1|
| geojsonio                   |    1|
| geojsonlint                 |    1|
| geonames                    |    1|
| geosphere                   |    1|
| ggalt                       |    1|
| ggm                         |    1|
| graphql                     |    1|
| GSODR                       |    1|
| gtools                      |    1|
| hash                        |    1|
| hexbin                      |    1|
| historydata                 |    1|
| Hmisc                       |    1|
| httpuv                      |    1|
| IRanges                     |    1|
| IRdisplay                   |    1|
| isdparser                   |    1|
| janeaustenr                 |    1|
| jpeg                        |    1|
| knitcitations               |    1|
| leafletR                    |    1|
| loggr                       |    1|
| magick                      |    1|
| mapdata                     |    1|
| markdown                    |    1|
| MCMCglmm                    |    1|
| memisc                      |    1|
| miniUI(&gt;=0.1.1)          |    1|
| mongolite                   |    1|
| nabor                       |    1|
| natserv                     |    1|
| openair                     |    1|
| openxlsx                    |    1|
| osmar                       |    1|
| outliers                    |    1|
| pander                      |    1|
| parallel                    |    1|
| plot3D                      |    1|
| plotKML                     |    1|
| plotly                      |    1|
| plumber                     |    1|
| progress                    |    1|
| protolite                   |    1|
| purrrlyr                    |    1|
| qlcMatrix                   |    1|
| RApiSerialize               |    1|
| rapport                     |    1|
| rbhl                        |    1|
| rbison                      |    1|
| rcdk                        |    1|
| Rcompression                |    1|
| readtext                    |    1|
| rebird                      |    1|
| RedisAPI                    |    1|
| redland                     |    1|
| redux                       |    1|
| reeack                      |    1|
| rfigshare                   |    1|
| ridigbio                    |    1|
| rinat                       |    1|
| ritis                       |    1|
| rJava                       |    1|
| RJSONIO                     |    1|
| rlist                       |    1|
| Rmpfr                       |    1|
| RMySQL                      |    1|
| rnaturalearthdata           |    1|
| rnaturalearthhires          |    1|
| rncl                        |    1|
| RNeXML                      |    1|
| rnoaa                       |    1|
| rnrfa                       |    1|
| ropenaq                     |    1|
| rotl                        |    1|
| rowr                        |    1|
| RPostgreSQL                 |    1|
| rrdf                        |    1|
| rredis                      |    1|
| rredlist                    |    1|
| rrlite                      |    1|
| RSclient                    |    1|
| RSelenium                   |    1|
| Rserve                      |    1|
| rstudioapi(&gt;=0.5)        |    1|
| rsvg                        |    1|
| rtracklayer                 |    1|
| RUnit                       |    1|
| S4Vectors                   |    1|
| sangerseqR                  |    1|
| scrapeR                     |    1|
| selectr                     |    1|
| seqinr                      |    1|
| shiny(&gt;=0.13.2)          |    1|
| snow                        |    1|
| SnowballC                   |    1|
| sofa                        |    1|
| spacetime                   |    1|
| spatstat                    |    1|
| SSOAP                       |    1|
| stringdist                  |    1|
| Suggests:testthat           |    1|
| Sxslt                       |    1|
| tabulizerjars               |    1|
| testthat(&gt;=0.7)          |    1|
| tidytext                    |    1|
| tidyverse                   |    1|
| tiff                        |    1|
| tmap                        |    1|
| USAboundaries               |    1|
| USAboundariesData           |    1|
| VariantAnnotation           |    1|
| vegan                       |    1|
| viridisLite                 |    1|
| wdman(&gt;=0.2.2)           |    1|
| weathermetrics              |    1|
| webmockr                    |    1|
| webshot                     |    1|
| wellknown                   |    1|
| WikidataR                   |    1|
| wikitaxa                    |    1|
| withr                       |    1|
| wordcloud2                  |    1|
| worrms                      |    1|
| XMLSchema                   |    1|
| xtable                      |    1|
| xts                         |    1|

``` r
#  summarise(count(name))
```

<!-- README.md is generated from README.Rmd. Please edit that file -->

# fakepackage

[![Travis build status](https://travis-ci.org/fakeorg/fakepackage.svg?branch=master)](https://travis-ci.org/fakeorg/fakepackage)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/fakeorg/fakepackage?branch=master&svg=true)](https://ci.appveyor.com/project/fakeorg/fakepackage)
[![Coverage status](https://codecov.io/gh/fakeorg/fakepackage/branch/master/graph/badge.svg)](https://codecov.io/github/fakeorg/fakepackage?branch=master)
[![CRAN status](http://www.r-pkg.org/badges/version/fakepackage)](https://cran.r-project.org/package=fakepackage)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![](https://badges.ropensci.org/24_status.svg)](https://github.com/ropensci/onboarding/issues/24)

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![Travis-CI Build Status](https://travis-ci.org/codemeta/codemetar.svg?branch=master)](https://travis-ci.org/codemeta/codemetar) [![Coverage Status](https://img.shields.io/codecov/c/github/codemeta/codemetar/master.svg)](https://codecov.io/github/codemeta/codemetar?branch=master)[![R build status](https://github.com/ropensci/codemetar/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/codemetar/actions) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/codemetar)](https://cran.r-project.org/package=codemetar)

<!-- README.md is generated from README.Rmd. Please edit that file -->
codemetar
=========

The goal of codemetar is to generate the JSON-LD file, `codemeta.json` containing software metadata describing an R package

Installation
------------

You can install codemetar from github with:

``` r
# install.packages("devtools")
devtools::install_github("codemeta/codemetar")
```

``` r
library("codemetar")
```

Example
-------

This is a basic example which shows you how to generate a `codemeta.json` for an R package (e.g. for `testthat`):

``` r
write_codemeta("testthat")
```

`codemetar` can take the path to the package root instead. This may allow `codemetar` to glean some additional information that is not available from the description file alone.

``` r
write_codemeta(".")
```

    {
      "@context": "https://raw.githubusercontent.com/codemeta/codemeta/master/codemeta.jsonld",
      "@type": "SoftwareSourceCode",
      "identifier": "codemetar",
      "description": "Codemeta defines a 'JSON-LD' format for describing software metadata.\n    This package provides utilities to generate, parse, and modify codemeta.jsonld\n    files automatically for R packages.",
      "name": "codemetar: Generate CodeMeta Metadata for R Packages",
      "codeRepository": "https://github.com/codemeta/codemetar",
      "issueTracker": "https://github.com/codemeta/codemetar/issues",
      "license": "https://spdx.org/licenses/MIT",
      "version": "0.1.0",
      "programmingLanguage": {
        "@type": "ComputerLanguage",
        "name": "R",
        "version": "3.4.0",
        "url": "https://r-project.org"
      },
      "runtimePlatform": "R version 3.4.0 (2017-04-21)",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Central R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "author": [
        {
          "@type": "Person",
          "givenName": "Carl",
          "familyName": "Boettiger",
          "email": "cboettig@gmail.com",
          "@id": "http://orcid.org/0000-0002-1642-628X"
        }
      ],
      "copyrightHolder": [
        {
          "@type": "Person",
          "givenName": "Carl",
          "familyName": "Boettiger",
          "email": "cboettig@gmail.com",
          "@id": "http://orcid.org/0000-0002-1642-628X"
        }
      ],
      "maintainer": {
        "@type": "Person",
        "givenName": "Carl",
        "familyName": "Boettiger",
        "email": "cboettig@gmail.com",
        "@id": "http://orcid.org/0000-0002-1642-628X"
      },
      "softwareSuggestions": [
        {
          "@type": "SoftwareApplication",
          "name": "testthat",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "jsonvalidate",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "covr",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "knitr",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "rmarkdown",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "httr",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "magrittr",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "readr",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "xml2",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        }
      ],
      "softwareRequirements": [
        {
          "@type": "SoftwareApplication",
          "name": "jsonlite",
          "version": "1.3",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "jsonld",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "git2r",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "devtools",
          "provider": {
            "@id": "https://cran.r-project.org",
            "@type": "Organization",
            "name": "Central R Archive Network (CRAN)",
            "url": "https://cran.r-project.org"
          }
        },
        {
          "@type": "SoftwareApplication",
          "name": "methods"
        },
        {
          "@type": "SoftwareApplication",
          "name": "R",
          "version": "3.0.0"
        }
      ],
      "contIntegration": "https://travis-ci.org/codemeta/codemetar",
      "developmentStatus": "wip",
      "releaseNotes": "https://github.com/codemeta/codemetar/blob/master/README.md",
      "readme": "https://github.com/codemeta/codemetar/blob/master/README.md",
      "fileSize": "119.263KB"
    }

Enriching CodeMeta metadata
---------------------------

The best way to ensure `codemeta.json` is as complete as possible is to begin by making full use of the fields that can be set in an R package DESCRIPTION file, such as `BugReports` and `URL`. Using the `Authors@R` notation allows a much richer specification of author roles, correct parsing of given vs family names, and email addresses.

In the current implementation, developers may specify an ORCID url for an author in the optional `comment` field of `Authors@R`, e.g.

    Authors@R: person("Carl", "Boettiger", role=c("aut", "cre", "cph"), email="cboettig@gmail.com", comment="http://orcid.org/0000-0002-1642-628X")

which will allow `codemetar` to associate an identifier with the person. This is clearly something of a hack since R's `person` object lacks an explicit notion of `id`, and may be frowned upon.

Additional metadata can be added by creating and manipulating a `codemeta` list in R:

``` r
cm <- create_codemeta(".")
cm$keywords <- list("metadata", "ropensci")
write_codemeta(cm)
```


<!-- README.md is generated from README.Rmd. Please edit that file -->
codemetar
=========

The goal of codemetar is to generate the JSON-LD file, `codemeta.json` containing software metadata describing an R package

Installation
------------

You can install codemetar from github with:

``` r
# install.packages("devtools")
devtools::install_github("codemeta/codemetar")
```

``` r
library("codemetar")
```

``` r
cm <- create_codemeta(".")
cm$keywords <- list("metadata", "ropensci")
write_codemeta(cm)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Travis-CI Build Status](https://travis-ci.org/ropensci/codemetar.svg?branch=master)](https://travis-ci.org/ropensci/codemetar) [![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/csawpip238vvbd72/branch/master?svg=true)](https://ci.appveyor.com/project/cboettig/codemetar/branch/master) [![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/codemetar/master.svg)](https://codecov.io/github/ropensci/codemetar?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/codemetar)](https://cran.r-project.org/package=codemetar) [![](http://badges.ropensci.org/130_status.svg)](https://github.com/ropensci/onboarding/issues/130) [![DOI](https://zenodo.org/badge/86626030.svg)](https://zenodo.org/badge/latestdoi/86626030)

<!-- README.md is generated from README.Rmd. Please edit that file -->
codemetar
=========

The goal of codemetar is to generate the JSON-LD file, `codemeta.json` containing software metadata describing an R package. For more general information about the CodeMeta Project for defining software metadata, see <https://codemeta.github.io>. In particular, new users might want to start with the [User Guide](https://codemeta.github.io/user-guide/), while those looking to learn more about JSON-LD and consuming existing codemeta files should see the [Developer Guide](https://codemeta.github.io/developer-guide/).

Installation
------------

You can install the latest version from CRAN using:

``` r
install.packages("codemetar")
```

You can also install the development version of `codemetar` from github with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/codemetar")
```

``` r
library("codemetar")
```

Example
-------

This is a basic example which shows you how to generate a `codemeta.json` for an R package (e.g. for `testthat`):

``` r
write_codemeta("testthat")
```

`codemetar` can take the path to the package root instead. This may allow `codemetar` to glean some additional information that is not available from the description file alone.

``` r
write_codemeta(".")
```

Which creates a file looking like so (first 10 lines; see full [codemeta.json here](https://github.com/codemeta/codemetar/blob/master/codemeta.json)):

    {
      "@context": [
        "http://purl.org/codemeta/2.0",
        "http://schema.org"
      ],
      "@type": "SoftwareSourceCode",
      "identifier": "codemetar",
      "description": "The 'Codemeta' Project defines a 'JSON-LD' format for describing\n  software metadata, as detailed at <https://codemeta.github.io>. This package\n  provides utilities to generate, parse, and modify 'codemeta.json' files \n  automatically for R packages, as well as tools and examples for working with\n  'codemeta.json' 'JSON-LD' more generally.",
      "name": "codemetar: Generate 'CodeMeta' Metadata for R Packages",
      "issueTracker": "https://github.com/ropensci/codemetar/issues",

Modifying or enriching CodeMeta metadata
----------------------------------------

The best way to ensure `codemeta.json` is as complete as possible is to begin by making full use of the fields that can be set in an R package DESCRIPTION file, such as `BugReports` and `URL`. Using the `Authors@R` notation allows a much richer specification of author roles, correct parsing of given vs family names, and email addresses.

In the current implementation, developers may specify an ORCID url for an author in the optional `comment` field of `Authors@R`, e.g.

    Authors@R: person("Carl", "Boettiger", role=c("aut", "cre", "cph"), email="cboettig@gmail.com", comment="http://orcid.org/0000-0002-1642-628X")

which will allow `codemetar` to associate an identifier with the person. If the package is hosted on CRAN, including the ORCiD in this way will cause an ORCiD logo and link to the ORCiD page to be added to the package CRAN webpage.

### Using the DESCRIPTION file

The DESCRIPTION file is the natural place to specify any metadata for an R package. The `codemetar` package can detect certain additional terms in the [CodeMeta context](https://codemeta.github.io/terms). Almost any additional codemeta field (see `codemetar:::additional_codemeta_terms` for a list) and can be added to and read from the DESCRIPTION into a `codemeta.json` file.

CRAN requires that you prefix any additional such terms to indicate the use of `schema.org` explicitly, e.g. `keywords` would be specified in a DESCRIPTION file as:

    X-schema.org-keywords: metadata, codemeta, ropensci, citation, credit, linked-data

Where applicable, these will override values otherwise guessed from the source repository. Use comma-separated lists to separate multiple values to a property, e.g. keywords.

See the [DESCRIPTION](https://github.com/codemeta/codemetar/blob/master/DESCRIPTION) file of the `codemetar` package for an example.

Going further
-------------

Check out all the [codemetar vignettes](https://docs.ropensci.org/codemetar/articles/index.html) for tutorials on other cool stuff you can do with codemeta and json-ld.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# fakepackage

[![Travis build
status](https://travis-ci.org/fakeorg/fakepackage.svg?branch=master)](https://travis-ci.org/fakeorg/fakepackage)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/fakeorg/fakepackage?branch=master&svg=true)](https://ci.appveyor.com/project/fakeorg/fakepackage)
[![Coverage
status](https://codecov.io/gh/fakeorg/fakepackage/branch/master/graph/badge.svg)](https://codecov.io/github/fakeorg/fakepackage?branch=master)
[![CRAN
status](http://www.r-pkg.org/badges/version/fakepackage)](https://cran.r-project.org/package=fakepackage)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

 [![Travis-CI Build Status](https://travis-ci.org/ropensci/codemetar.svg?branch=master)](https://travis-ci.org/ropensci/codemetar) [![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/csawpip238vvbd72/branch/master?svg=true)](https://ci.appveyor.com/project/cboettig/codemetar/branch/master) [![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/codemetar/master.svg)](https://codecov.io/github/ropensci/codemetar?branch=master)  [![](http://badges.ropensci.org/130_status.svg)](https://github.com/ropensci/onboarding/issues/130) [![DOI](https://zenodo.org/badge/86626030.svg)](https://zenodo.org/badge/latestdoi/86626030)

<!-- README.md is generated from README.Rmd. Please edit that file -->
codemetar
=========

The goal of codemetar is to generate the JSON-LD file, `codemeta.json` containing software metadata describing an R package. For more general information about the CodeMeta Project for defining software metadata, see <https://codemeta.github.io>. In particular, new users might want to start with the [User Guide](https://codemeta.github.io/user-guide/), while those looking to learn more about JSON-LD and consuming existing codemeta files should see the [Developer Guide](https://codemeta.github.io/developer-guide/).

Installation
------------

You can install the latest version from CRAN using:

``` r
install.packages("codemetar")
```

You can also install the development version of `codemetar` from github with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/codemetar")
```

``` r
library("codemetar")
```

Example
-------

This is a basic example which shows you how to generate a `codemeta.json` for an R package (e.g. for `testthat`):

``` r
write_codemeta("testthat")
```

`codemetar` can take the path to the package root instead. This may allow `codemetar` to glean some additional information that is not available from the description file alone.

``` r
write_codemeta(".")
```

Which creates a file looking like so (first 10 lines; see full [codemeta.json here](https://github.com/codemeta/codemetar/blob/master/codemeta.json)):

    {
      "@context": [
        "http://purl.org/codemeta/2.0",
        "http://schema.org"
      ],
      "@type": "SoftwareSourceCode",
      "identifier": "codemetar",
      "description": "The 'Codemeta' Project defines a 'JSON-LD' format for describing\n  software metadata, as detailed at <https://codemeta.github.io>. This package\n  provides utilities to generate, parse, and modify 'codemeta.json' files \n  automatically for R packages, as well as tools and examples for working with\n  'codemeta.json' 'JSON-LD' more generally.",
      "name": "codemetar: Generate 'CodeMeta' Metadata for R Packages",
      "issueTracker": "https://github.com/ropensci/codemetar/issues",

Modifying or enriching CodeMeta metadata
----------------------------------------

The best way to ensure `codemeta.json` is as complete as possible is to begin by making full use of the fields that can be set in an R package DESCRIPTION file, such as `BugReports` and `URL`. Using the `Authors@R` notation allows a much richer specification of author roles, correct parsing of given vs family names, and email addresses.

In the current implementation, developers may specify an ORCID url for an author in the optional `comment` field of `Authors@R`, e.g.

    Authors@R: person("Carl", "Boettiger", role=c("aut", "cre", "cph"), email="cboettig@gmail.com", comment="http://orcid.org/0000-0002-1642-628X")

which will allow `codemetar` to associate an identifier with the person. If the package is hosted on CRAN, including the ORCiD in this way will cause an ORCiD logo and link to the ORCiD page to be added to the package CRAN webpage.

### Using the DESCRIPTION file

The DESCRIPTION file is the natural place to specify any metadata for an R package. The `codemetar` package can detect certain additional terms in the [CodeMeta context](https://codemeta.github.io/terms). Almost any additional codemeta field (see `codemetar:::additional_codemeta_terms` for a list) and can be added to and read from the DESCRIPTION into a `codemeta.json` file.

CRAN requires that you prefix any additional such terms to indicate the use of `schema.org` explicitly, e.g. `keywords` would be specified in a DESCRIPTION file as:

    X-schema.org-keywords: metadata, codemeta, ropensci, citation, credit, linked-data

Where applicable, these will override values otherwise guessed from the source repository. Use comma-separated lists to separate multiple values to a property, e.g. keywords.

See the [DESCRIPTION](https://github.com/codemeta/codemetar/blob/master/DESCRIPTION) file of the `codemetar` package for an example.

Going further
-------------

Check out all the [codemetar vignettes](https://codemeta.github.io/codemetar/articles/index.html) for tutorials on other cool stuff you can do with codemeta and json-ld.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# fakepackage

[![Travis build status](https://travis-ci.org/fakeorg/fakepackage.svg?branch=master)](https://travis-ci.org/fakeorg/fakepackage)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/fakeorg/fakepackage?branch=master&svg=true)](https://ci.appveyor.com/project/fakeorg/fakepackage)
[![Coverage status](https://codecov.io/gh/fakeorg/fakepackage/branch/master/graph/badge.svg)](https://codecov.io/github/fakeorg/fakepackage?branch=master)
[![CRAN status](http://www.r-pkg.org/badges/version/fakepackage)](https://cran.r-project.org/package=fakepackage)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![](https://badges.ropensci.org/9999999999_status.svg)](https://github.com/ropensci/onboarding/issues/9999999999)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# fakepackage

[![Travis build
status](https://travis-ci.org/fakeorg/fakepackage.svg?branch=master)](https://travis-ci.org/fakeorg/fakepackage)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/fakeorg/fakepackage?branch=master&svg=true)](https://ci.appveyor.com/project/fakeorg/fakepackage)
[![Coverage
status](https://codecov.io/gh/fakeorg/fakepackage/branch/master/graph/badge.svg)](https://codecov.io/github/fakeorg/fakepackage?branch=master)
[![CRAN
status](http://www.r-pkg.org/badges/version/fakepackage)](https://cran.r-project.org/package=fakepackage)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# fakepackage

[![Travis build status](https://travis-ci.org/fakeorg/fakepackage.svg?branch=master)](https://travis-ci.org/fakeorg/fakepackage)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/fakeorg/fakepackage?branch=master&svg=true)](https://ci.appveyor.com/project/fakeorg/fakepackage)
[![Coverage status](https://codecov.io/gh/fakeorg/fakepackage/branch/master/graph/badge.svg)](https://codecov.io/github/fakeorg/fakepackage?branch=master)
[![CRAN status](http://www.r-pkg.org/badges/version/fakepackage)](https://cran.r-project.org/package=fakepackage)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<table class="table"><thead><tr class="header">
<th align="left">Release</th>
<th align="left">Usage</th>
<th align="left">Development</th>
</tr></thead>
<tbody>
<tr class="odd">
<td align="left"><a href="https://doi.org/"><img src="http://joss.theoj.org/papers.svg" alt="JOSS"></a></td>
<td align="left"><a href="https://www.gnu.org/licenses/gpl-3.0.en.html"><img src="https://img.shields.io/badge/licence-GPL--3-blue.svg" alt="Licence"></a></td>
<td align="left"><a href="https://ci.appveyor.com/project/fakeorg/fakepackage"><img src="https://ci.appveyor.com/api/projects/status/4ypc9xnmqt70j94e?svg=true&amp;branch=master" alt="AppVeyor"></a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://github.com/fakeorg/onboarding/issues/156"><img src="https://badges.fakeorg.org/156_status.svg" alt="fakeorg"></a></td>
<td align="left"><a href="https://cran.r-project.org/"><img src="https://img.shields.io/badge/R%3E%3D-3.2.0-blue.svg" alt="Project Status: Active – The project has reached a stable, usable state and is being actively developed."></a></td>
<td align="left"><a href="https://travis-ci.org/fakeorg/fakepackage"><img src="https://travis-ci.org/fakeorg/fakepackage.svg?branch=master" alt="Travis"></a></td>
</tr>
<tr class="odd">
<td align="left"><a href="http://cran.r-project.org/package=fakepackage"><img src="http://www.r-pkg.org/badges/version/fakepackage" alt="CRAN"></a></td>
<td align="left"><a href="http://cran.rstudio.com/package=fakepackage"><img src="http://cranlogs.r-pkg.org/badges/fakepackage" alt="downloads"></a></td>
<td align="left"><a href="https://codecov.io/github/fakeorg/fakepackage?branch=master"><img src="https://codecov.io/github/fakeorg/fakepackage/coverage.svg?branch=master" alt="Codecov"></a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://zenodo.org/badge/latestdoi/82609103"><img src="https://zenodo.org/badge/82609103.svg" alt="Zenodo"></a></td>
<td align="left"></td>
<td align="left"></td>
</tr>
</tbody>
</table>
<br>


[![fakeorg_footer](http://fakeorg.org/public_images/github_footer.png)](https://fakeorg.org)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# fakepackage

[![Travis build
status](https://travis-ci.org/fakeorg/fakepackage.svg?branch=master)](https://travis-ci.org/fakeorg/fakepackage)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/fakeorg/fakepackage?branch=master&svg=true)](https://ci.appveyor.com/project/fakeorg/fakepackage)
[![Coverage
status](https://codecov.io/gh/fakeorg/fakepackage/branch/master/graph/badge.svg)](https://codecov.io/github/fakeorg/fakepackage?branch=master)
[![CRAN
status](http://www.r-pkg.org/badges/version/fakepackage)](https://cran.r-project.org/package=fakepackage)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
---
output: github_document
---


<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  comment = "",
  fig.path = "tools/README-"
)
```

# codemetar

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![R build status](https://github.com/ropensci/codemetar/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/codemetar/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/codemetar/master.svg)](https://codecov.io/github/ropensci/codemetar?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/codemetar)](https://cran.r-project.org/package=codemetar)
[![](http://badges.ropensci.org/130_status.svg)](https://github.com/ropensci/software-review/issues/130)
[![DOI](https://zenodo.org/badge/86626030.svg)](https://zenodo.org/badge/latestdoi/86626030)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/codemetar)](https://CRAN.R-project.org/package=codemetar)


The [codemeta](https://cran.r-project.org/package=codemeta) package provides a more minimalist approach to generating codemeta based only on DESCRIPTION and CITATION files,
while `codemetar` provides additional abilities to detect metadata from README and GitHub sources, and provide
more user feedback, suggestions, and messaging.



```{r child='man/rmdhunks/goals.Rmd'} 
```

## Why create a codemeta.json for your package?


```{r child='man/rmdhunks/whybother.Rmd'} 
```


## Create a codemeta.json in one function call

```{r child='man/rmdhunks/example.Rmd'} 
```
## Keep codemeta.json up-to-date

```{r child='man/rmdhunks/uptodate.Rmd'} 
```

## How to improve your package's codemeta.json?

```{r child='man/rmdhunks/improvecodemetadata.Rmd'} 
```

## Installation and usage requirements

```{r child='man/rmdhunks/install-instructions.Rmd'} 
```


```{r child='man/rmdhunks/internet.Rmd'} 
```


## Going further

Check out all the [codemetar man](https://docs.ropensci.org/codemetar/articles/index.html) for tutorials on other cool stuff you can do with codemeta and json-ld.  

```{r child='man/rmdhunks/minimeta.Rmd'} 
```

```{r addscript, results="asis", echo=FALSE}
glue::glue('<script type="application/ld+json">
      {glue::glue_collapse(readLines(file.path("inst", "schemaorg.json")), sep = "\n")}
    </script>')
```


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: "Translating between schema using JSON-LD"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Translating between various software metaata formats with JSON-LD in codemetar}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



```{r include=FALSE}
knitr::opts_chunk$set(comment="")
if(grepl("windows", tolower(Sys.info()[["sysname"]])))
  knitr::opts_chunk$set(comment="", error =TRUE)
```

```{r message=FALSE}
library("codemetar")
library("magrittr")
library("jsonlite")
library("jsonld")
```

## JSON-LD transforms: Expansion and Compaction

One of the central motivations of JSON-LD is making it easy to translate between different representations of what are fundamentally the same data types. Doing so uses the two core algorithms of JSON-LD: *expansion* and *compaction*, as [this excellent short video by JSON-LD creator Manu Sporny](https://www.youtube.com/watch?v=Tm3fD89dqRE) describes.

Here's how we would use JSON-LD (from R) to translate between the two examples of JSON data from different providers as shown in the video.  First, the JSON from the original provider:

```{r}
ex <-
'{
"@context":{
  "shouter": "http://schema.org/name",
  "txt": "http://schema.org/commentText"
},
"shouter": "Jim",
"txt": "Hello World!"
}'
```

Next, we need the context of the second data provider.  This will let us translate the JSON format used by provider one ("Shouttr") to the second ("BigHash"):

```{r}
bighash_context <- 
'{
"@context":{
  "user": "http://schema.org/name",
  "comment": "http://schema.org/commentText"
}
}'
```

With this in place, we simply expand the original JSON and then compact using the new context:

```{r}
jsonld_expand(ex) %>%
  jsonld_compact(context = bighash_context)
```

## Crosswalking

The CodeMeta Crosswalk table seeks to accomplish a very similar goal.  The crosswalk table provides a human-readable mapping of different software metadata providers into the codemeta context (an extension of schema.org).  For instance, we'll read in some data from GitHub:


### GitHub

Here we crosswalk the JSON data returned as the repository information from the GitHub API: 

```{r eval = FALSE}
repo_info <- gh::gh("/repos/:owner/:repo", owner = "ropensci", repo = "EML")
```

```{r include = FALSE}
## Actually, we'll use a chached copy and not eval the above chunk to avoid a `gh` dependency:
repo_info <- read_json(system.file("examples/github_format.json", package = "codemetar"))
```

Let's just take a look at what the returned json data looks like:

```{r}
repo_info %>% toJSON()
```

We can crosswalk this information into codemeta just by supplying the column name to the `crosswalk` function.  This performs the same expansion of the metadata in the GitHub context, followed by compaction into the codemeta context.  Note that terms not recognized/included in the codemeta context will be dropped:

```{r}
github_meta <- crosswalk(repo_info, "GitHub")
github_meta
```

We can verify that the result is a valid codemeta document:

```{r}
codemeta_validate(github_meta)
```


## Transforming into other column schema


The above transform showed the process of going from plain JSON data into the codemeta standard serialization.  Similarly, we can crosswalk into any of the other columns in the crosswalk table.  To do so, the `crosswalk` function will first expand any of the recognized properties into the codemeta JSON-LD context, just as above.  Unrecognized properties are dropped, since there is no consensus context into which we can expand them.  Second, the expanded terms are then compacted down into the new context (Zenodo in this case.)  This time,
any terms that are not part of the codemeta context are kept, but not compacted, since they still have meaningful contexts (that is, full [URI](https://en.wikipedia.org/wiki/Uniform_Resource_Identifier)s, e.g. URLs) that can be associated with them:

```{r}
crosswalk(repo_info, "GitHub", "Zenodo") %>%
drop_context()
```

Thus terms that still have a uncompacted prefix like `schema:` or `codemeta:` reflect properties that we could successfully extract from the input data, but do not have corresponding properties in the Zenodo context.  This is the standard behavior of the 
compaction algorithm: unrecognized fields are not dropped, but are also not compacted,
making them accessible only if referenced explicitly.  


## NodeJS example

NodeJS uses a `package.json` format that is very similar to a simple `codemeta.json` file, though it is not Linked Data as it does not declare a context.  Here we crosswalk an example `package.json` file into proper `codemeta` standard.  

```{r}
package.json <- read_json(
"https://raw.githubusercontent.com/heroku/node-js-sample/master/package.json")
package.json
```


```{r}
crosswalk(package.json, "NodeJS")
```

Note that while nested structures per se pose no special problem, the compaction/expansion paradigm lacks a mechanism to capture differences in nesting between schema.  For instance,
in `codemeta` (i.e. in schema.org), a `codeRepository` is expected to be a URL, while NodeJS `package.json` permits it to be another object node with sub-properties `type` and `url`.  There is no way in JSON-LD transforms or context definitions to indicate that the `url` sub-property in the NodeJS case, e.g. `codeRepository.url` maps to schema's `codeRepository`. 
(This same limitation is also true of the 2-dimensional table structure of the crosswalk itself, though it is important to keep in mind that this 1:1 mapping requirement is not unique to the the `.csv` representation but also inherent in JSON-LD contexts.)  

Consequently, a thorough translation between formats that do not provide their own JSON-LD contexts will ultimately require more manual implementation, which would be expressed within a particular programming language (e.g. R) rather than in the generic algorithms of JSON-LD available in many common programming languages.  




---
draft: true
---


```{r message=FALSE}
library(httr)
library(magrittr)
library(jsonld)
library(jsonlite)

```


```{r}
date_parts <- function(parts)
  do.call(paste,
          list(unlist(parts),
               collapse = "-"))

author_parse <- function(authors){
  lapply(authors, function(a){
    list("@type" = "Person",
      givenName = a$given,
         familyName = a$family
  #       affiliation = list("@type" = "Organization", name = a$affiliation[[1]]$name) # handle if empty
         )
  })
}
```

```{r}
csl_to_schema <- function(doi, as = "list", silent=FALSE){
  message(doi)
csl <-
  content(
    GET(doi, add_headers(Accept="application/vnd.citationstyles.csl+json")),
    as = "parsed", type = "application/json")

schema <- 
list( "@type" = "ScholarlyArticle",
  "@id" = paste0("https://doi.org/", csl$DOI),
  "name" = csl$title,
  
  "dateCreated" = date_parts(csl$created$`date-parts`), 
  "isPartOf" = list(
    "@type" = "PublicationIssue",
    "issueNumber" = csl$issue,
    "datePublished" = date_parts(csl$`published-print`$`date-parts`), 
    "isPartOf" = list(
      "@type"= c("Periodical", "PublicationVolume"),
      "volumeNumber" = csl$volume,
      "name"= csl$`container-title`, # or csl$`journal-title`
      "publisher" = csl$publisher
    )
  ), 
  "author" = author_parse(csl$author),
  "pagination" = csl$page,
  license = csl$license[[1]]$URL) 
 
  switch(as, 
         "list" = return(schema),
         "json" = return(toJSON(schema, auto_unbox=TRUE, pretty=TRUE)),
         "file" = return(write_json(schema, "csl.json", auto_unbox=TRUE, pretty=TRUE))
         )
  
}
csl_to_schema(doi = "https://doi.org/10.1002/ece3.2314", as="file")


```



```{r}
c(
"https://doi.org/10.1093/biosci/bix025",
"https://doi.org/10.1002/ece3.2314",
"https://doi.org/10.1890/15-0236",
"https://doi.org/10.5334/jors.bu",
"https://doi.org/10.1111/2041-210X.12469",
"https://doi.org/10.1145/2723872.2723882",
"https://doi.org/10.1098/rspb.2014.1631",
"https://doi.org/10.1098/rspb.2013.1372",
"https://doi.org/10.1007/s12080-013-0192-6",
"https://doi.org/10.1038/493157a",
"https://doi.org/10.1098/rspb.2012.2085",
"https://doi.org/10.1111/j.1095-8649.2012.03464.x",
"https://doi.org/10.1111/j.2041-210X.2012.00247.x",
"https://doi.org/10.1098/rsif.2012.0125",
"https://doi.org/10.1111/j.1558-5646.2012.01619.x", 
"https://doi.org/10.1111/j.1558-5646.2011.01574.x",
"https://doi.org/10.1016/j.tpb.2009.10.003",
"https://doi.org/10.1086/508600"
) %>% lapply(csl_to_schema, as="list") -> out


out %>%
  write_json(schema, "csl.json", auto_unbox=TRUE, pretty=TRUE)
	 



```







```{r}
"https://doi.org/10.1093/biosci/bix025" %>%
  GET(add_headers(Accept="application/rdf+xml")) %>%
  content(as = "parsed", type = "application/xml") %>%
  write_xml("rdf.xml")
```
```{r}
library("gh")
library("jsonlite")
library("jsonld")
library("codemetar")
library("magrittr")
library("readr")
library("dplyr")
library("purrr")
```

```{r message=FALSE}
crosswalk <- function(x,
                      schema,
                      codemeta_context = 
             "https://doi.org/10.5063/schema/codemeta-2.0"
                      ){
  
  df <- crosswalk_table(schema)
  context <- crosswalk_context(df, codemeta_context)
  crosswalk_transform(x, 
                      crosswalk_context = context,
                      codemeta_context = codemeta_context)
  
}

crosswalk_table <- function(schema){
  df <-
    read_csv("https://github.com/codemeta/codemeta/raw/master/crosswalk.csv",
             col_types = cols(.default = "c"))
  df <- df[c("Property", schema)]
  df[!is.na(df[[schema]]),] ## stats::na.omit
}

crosswalk_context <- 
  function(df, codemeta_context = 
             "https://doi.org/10.5063/schema/codemeta-2.0"){
    
  context <- jsonlite::read_json(codemeta_context)
  context[[1]][["id"]] <- NULL ## avoid collisions with @id
  properties <- names(context[[1]])  
    
  ## replace property names according to crosswalk data.frame
  for(i in 1:dim(df)[1]){
   properties[properties == df[i, 1][[1]] ]  <- df[i,2]
  }
  names(context[[1]]) <- properties
  context
}

crosswalk_transform <- function(x, 
                                crosswalk_context = NULL,
                                codemeta_context = 
             "https://doi.org/10.5063/schema/codemeta-2.0"){
  
  ## FIXME consider case of x is `json` class, and x is file-path / url
  
  if(is.null(x$`@context`) && is.null(crosswalk_context)){
    stop("crosswalk_context required if input json has no '@context' property")
  }
  if(is.null(x$`@context`))
    x$`@context` <- crosswalk_context[[1]]

  y <- toJSON(x, auto_unbox = TRUE, pretty = TRUE)
  y <- jsonld_expand(y)
  y <- jsonld_compact(y, context = codemeta_context)
  y
}

```

Example: 

```{r}
owner <- "ropensci"
repo <- "EML"
r <- gh("/repos/:owner/:repo", owner = owner, repo = repo)
crosswalk(r, "GitHub")
```


Add Contributor information. Get contributor list, then contributor names

```{r}
contribs <- gh("/repos/:owner/:repo/contributors", owner = owner, repo = repo)

author_threshold <- 2
authors <- contribs %>% map_int("contributions") >= author_threshold

## Sort authors by contribution number.  Seems contribs already does this?
contribs[authors] %>%
  map_chr("login") %>%
  map(function(user) gh("/users/:user", user = user)) %>%
  map_chr("name")

```


Other mappable info: user:company, user:blog, user:email (API doesn't return email)

```{r}
gh("/users/:user", user = "cboettig")
```

Use the preview version of the API to get license information

```{r}

r <- gh("/repos/:owner/:repo", owner = owner, repo = repo, .send_headers = c(ACCEPT = "application/vnd.github.drax-preview+json"))

```
---
title: "CodeMeta for rOpenSci packages"
author: "Carl Boettiger"
date: "5/13/2017"
output: github_document
---

```{r include=FALSE}
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
```

All rOpenSci packages via GitHub API:

```{r}
# remotes::install_github("r-pkgs/gh")
library("gh")

repos <- gh("/users/:username/repos", username = "ropensci", .limit = Inf)
repo_names <- vapply(repos, "[[", "", "name")
#not_a_pkg <- grepl("[-_]", repo_names)
#pkgs <- repo_names[!not_a_pkg]
writeLines(paste("ropensci", repo_names, sep="/"), "ropensci_gh.txt")

## Clone all, check for which ones have DESCRIPTION files

```

All rOpenSci packages via ropkgs: 

```{r}
#remotes::install_github("ropensci/ropkgs")
library("ropkgs")
out <- ro_pkgs()
good <- out$packages$status == "good"
installable <- out$packages$installable

pkgs <- gsub("https://github.com/", "", out$packages$url)[installable & good]

```

```{r}
#pkgs <- pkgs[1:3]
dir.create("pkg_src"); 
setwd("pkg_src")
for(p in pkgs){
  system(paste0("git clone https://github.com/", p, " ", p))
}
```

```{r}
setwd("pkg_src") # setwd doesn't persist across chunks
cm <- lapply(pkgs, function(p){
  message(p)
  codemetar::create_codemeta(p)
})
```


```{r}
codemetar::write_codemeta(cm, path="ropensci.json")
```



---
title: "sparql.Rmd"
author: "Carl Boettiger"
date: "8/10/2017"
output: html_document
---

```{r message=FALSE}
library("jsonld")
library("jsonlite")
library("redland")
library("tidyverse")
library("codemetar")
```

# Practice examples

We'll start by generating a codemeta json-ld file, and converting this to RDF using the `jsonld` library.  The official `jsonld` serialization uses `nquads` format for RDF, (it is not clear if any others are supported at this time)

```{r}
write_codemeta("codemetar")
jsonld::jsonld_to_rdf("codemeta.json") %>% 
  writeLines("codemeta.nquads")
```

Now we can use the `redland` library to parse and query the resulting RDF.  Note that the we must tell `rdflib` what parser and mime type to use, matching the names and `q 1.0` recognized mimeTypes as listed here: <http://librdf.org/raptor/api/raptor-formats-types-by-parser.html>.

```{r}
world <- new("World")
storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
model <- new("Model", world=world, storage, options="")
parser <- new("Parser", world, name="nquads", mimeType="text/x-nquads")
parseFileIntoModel(parser, world, "codemeta.nquads", model)
```

We now have a `model` object loaded into memory.  We can serialize this into different formats, e.g. rdf+xml:

```{r}
serializer <- new("Serializer", world, name="rdfxml", mimeType="application/rdf+xml")
status <- serializeToFile(serializer, world, model, "codemeta.xml")
```


against which we can make SPARQL queries:


```{r}
queryString <- '
PREFIX schema: <http://schema.org>
SELECT ?name
WHERE { 
  ?a schema:name ?name .
}'
query <- new("Query", world, queryString, base_uri=NULL, query_language="sparql", query_uri=NULL)
queryResult <- executeQuery(query, model)
result <- getNextResult(queryResult)

result
```





--------


```{r}
xmlrdf <- 
'<?xml version="1.0"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
     xmlns:dc="http://purl.org/dc/elements/1.1/">
  <rdf:Description rdf:about="http://www.dajobe.org/">
    <dc:title>Dave Becketts Home Page</dc:title>
    <rdf:type>Webpage</rdf:type>
    <dc:creator>Dave Beckett</dc:creator>
    <dc:address>
      <dc:state>California</dc:state>
    </dc:address>
    <dc:description>The generic home page of Dave Beckett.</dc:description>
  </rdf:Description> 
</rdf:RDF>
'
writeLines(xmlrdf, "ex.rdf")
```

```{r}
#f <- system.file("extdata", "dc.rdf", package="redland")
#cat(readLines(f), sep="\n")
f <- "ex.rdf"
parseFileIntoModel(parser, world, f, model)

```

```{r}
queryString <- '
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
SELECT *
WHERE { 
  ?a dc:creator ?creator ;
     dc:title ?title ;
     dc:description ?description .
  FILTER (?title = "Dave Beckett\'s Home Page") .

}'
query <- new("Query", world, queryString, base_uri=NULL, query_language="sparql", query_uri=NULL)
queryResult <- executeQuery(query, model)
result <- getNextResult(queryResult)

result
```


------

# rOpenSci corpus



```{r}
#nquads <- jsonld_to_rdf("ropensci.json") 
#xml <- jsonld_to_rdf("ropensci.json", options = list(format="application/rdf+xml")) 
#writeLines(nquads, "ropensci.nquads")
```

```{r}
world <- new("World")
storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
model <- new("Model", world, storage, options="")
# Create the default "rdfxml" parser
parser <- new("Parser", world)
parseFileIntoModel(parser, world, "ropensci.rdf", model)

```


```{r}
queryString <- '
PREFIX schema: <http://schema.org>
SELECT ?name
WHERE { 
  ?a schema:name ?name .
}'
query <- new("Query", world, queryString, base_uri=NULL, query_language="sparql", query_uri=NULL)
queryResult <- executeQuery(query, model)
result <- getNextResult(queryResult)

result
```
---
title: "Intro to codemetar, Codemeta creator for R packages"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Intro to codemetar, Codemeta creator for R packages}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r include=FALSE}
Sys.setenv("ON_CRAN" = "true")
knitr::opts_chunk$set(comment="")
if(grepl("windows", tolower(Sys.info()[["sysname"]])))
  knitr::opts_chunk$set(comment="", error =TRUE)
```



# codemetar: generate codemeta metadata for R packages

```{r child='../man/rmdhunks/goals.Rmd'} 
```


## Why create a codemeta.json for your package?

```{r child='../man/rmdhunks/whybother.Rmd'} 
```

## Installation and usage requirements

```{r child='../man/rmdhunks/install-instructions.Rmd'} 
```


```{r child='../man/rmdhunks/internet.Rmd'} 
```

## Create a codemeta.json in one function call

```{r child='../man/rmdhunks/example.Rmd'} 
```

## Keep codemeta.json up-to-date

```{r child='../man/rmdhunks/uptodate.Rmd'} 
```

## A brief intro to common terms we'll use:

- [Linked data](https://en.wikipedia.org/wiki/Linked_data): We often use different words to mean the same thing. And sometimes the same word to mean different things.  Linked data seeks to address this issue by using [URIs](https://en.wikipedia.org/wiki/Uniform_Resource_Identifier) (i.e. URLs) to make this explicit.

- context: No one likes typing out long URLs all the time.  So instead, the *context* of a JSON-LD file (`"@context"` element) gives us the context for the terms we use, that is, the root URL.  Usually schema.org but domain specific ones also (eg codemeta)

- [Schema.org](https://schema.org/): A major initiative led by Google and other search engines to define a simple and widely used context to link data on the web through a catalogue of standard metadata fields

- [The CodeMeta Project](https://codemeta.github.io/): an academic led community initiative to formalise the metadata fields included in typical software metadata records and introduce important fields that did not have clear equivalents. The codemeta [crosswalk](https://codemeta.github.io/crosswalk/) provides an explicit map between the metadata fields used by a broad range of software repositories, registries and archives

- [JSON-LD](https://json-ld.org): While 'linked data' can be represented in many different formats, these have consistently proven a bit tricky to use, either for consumers or developers or both.  JSON-LD provides a simple adaptation of the JSON format, which has proven much more popular with both audiences, that allows it to express (most) linked-data concepts.  It is now the format of choice for expressing linked data by Google and many others.  Any JSON-LD file is valid JSON, and any JSON file can be treated as JSON-LD.  

- [codemetar](https://codemeta.github.io/codemetar/): The CodeMeta Project has created tools in several languages to implement the CodeMeta Crosswalk (using JSON-LD) and help extract software metadata into `codemeta.json` records.  `codemetar` is one such tool, focused on R and R packages.

## How to improve your package's codemeta.json?

```{r child='../man/rmdhunks/improvecodemetadata.Rmd'} 
```

## Going further

Check out all the [codemetar vignettes](https://docs.ropensci.org/codemetar/articles/index.html) for tutorials on other cool stuff you can do with codemeta and json-ld.  

```{r child='../man/rmdhunks/minimeta.Rmd'} 
```
---
title: "Parsing Codemeta Data"
subtitle: "Illustrate the kind of information we can discover by parsing collections of codemeta documents programmatically"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
---


```{r include=FALSE}
knitr::opts_chunk$set(comment="", warning = FALSE)
if(grepl("windows", tolower(Sys.info()[["sysname"]])))
  knitr::opts_chunk$set(comment="", error =TRUE)
```


Here we illustrate some example use cases that involve parsing codemeta data.  

```{r message=FALSE}
library(jsonld)
library(jsonlite)
library(magrittr)
library(codemetar)
library(purrr)
library(dplyr)
library(printr)
library(tibble)
```

We start with a simple example from the `codemeta.json` file of `codemetar` itself.  First, we'll just generate a copy of the codemeta record for the package:

```{r, eval=FALSE}
write_codemeta("codemetar", "codemeta.json")
```

```{r include=FALSE}

# Fix, the vignette is not parsing due to an error when writing:

codemeta <- create_codemeta("codemetar")

jsonlite::write_json(codemeta, "codemeta.json", pretty = TRUE, auto_unbox = TRUE)
  

```


We then digest this input using a JSON-LD "frame." While not strictly necessary, this helps ensure the data matches the format we expect, even if the original file had errors or missing data.  See the vignette "Validating in JSON-LD" in this package and the official [JSON-LD docs](https://json-ld.org/spec/latest/json-ld-framing/) for details).  The `codemetar` package includes a reasonably explicit frame to get us started:

```{r}
frame <- system.file("schema/frame_schema.json", package="codemetar")

meta <- 
  jsonld_frame("codemeta.json", frame) %>%
  fromJSON(FALSE) %>% getElement("@graph") %>% getElement(1)
```

Construct a citation

```{r}
authors <- 
lapply(meta$author, 
       function(author) 
         person(given = author$given, 
                family = author$family, 
                email = author$email,
                role = "aut"))
year <- meta$datePublished
if(is.null(year)) 
  year <- format(Sys.Date(), "%Y")
bibitem <- 
 bibentry(
     bibtype = "Manual",
     title = meta$name,
     author = authors,
     year = year,
     note = paste0("R package version ", meta$version),
     url = meta$URL,
     key = meta$identifier
   )

cat(format(bibitem, "bibtex"))

bibitem
```


## Parsing the ropensci corpus

The ropensci corpus consists of a list of codemeta files for all packages provided by the rOpenSci project, <ropensci.org>.  This provides a good test-case for how a large collection of codemeta files can be manipulated to help us get a better picture of the corpus.  

```{r}
download.file("https://github.com/codemeta/codemetar/raw/master/inst/notebook/ropensci.json",
              "ropensci.json")
```


As before, it is helpful, though not essential, to start off by framing the input data.

```{r}
frame <- system.file("schema/frame_schema.json", package="codemetar")

corpus <- 
    jsonld_frame("ropensci.json", frame) %>%
    fromJSON(simplifyVector = FALSE) %>%
    getElement("@graph") 
```

We're now ready to start exploring.  As usual, functions from `purrr` prove very useful for iterating through large JSON files.  First, we look at some basic summary data:

```{r}
## deal with nulls explicitly by starting with map
pkgs <- map(corpus, "name") %>% compact() %>% as.character()

# keep only those with package identifiers (names)
keep <- map_lgl(corpus, ~ length(.x$identifier) > 0)
corpus <- corpus[keep]

## now we can just do
all_pkgs <- map_chr(corpus, "name")
head(all_pkgs)
```

```{r}
## 60 unique maintainers
map_chr(corpus, c("maintainer", "familyName")) %>% unique() %>% length()

## Mostly Scott
map_chr(corpus, c("maintainer", "familyName")) %>% 
  as_tibble() %>%
  group_by(value) %>%
  tally(sort=TRUE)
```


```{r}
## number of co-authors ... 
map_int(corpus, function(r) length(r$author)) %>% 
  as_tibble() %>%
  group_by(value) %>%
  tally(sort=TRUE)
```

```{r}
## Contributors isn't used as much...
map_int(corpus, function(r) length(r$contributor)) %>% 
  as_tibble() %>%
  group_by(value) %>%
  tally(sort=TRUE)

```

Numbers (n) of packages with a total of (value) dependencies:

```{r}
map_int(corpus, function(r) length(r$softwareRequirements))  %>% 
  as_tibble() %>%
  group_by(value) %>%
  tally(sort=TRUE)
```

which dependencies are used most frequently?

```{r}
corpus %>%
map_df(function(x){
  ## single, unboxed dep
  if("name" %in% names(x$softwareRequirements))
    dep <- x$name
  else if("name" %in% names(x$softwareRequirements[[1]]))
    dep <- map_chr(x$softwareRequirements, "name")
  else { ## No requirementsß
    dep <- NA
  }
  
  tibble(identifier = x$identifier, dep = dep)
}) -> dep_df


dep_df %>%
group_by(dep) %>% 
  tally(sort = TRUE)

```

Alternate approach using a frame instead of `purrr` functions for subsetting the data.  Note that this gets all Depends and suggests (really all `SoftwareApplication` types mentioned)

```{r}
dep_frame <- '{
  "@context": "https://raw.githubusercontent.com/codemeta/codemeta/master/codemeta.jsonld",
  "@explicit": "true",
  "name": {}
}'
jsonld_frame("ropensci.json", dep_frame) %>% 
  fromJSON() %>% 
  getElement("@graph") %>%
  filter(type == "SoftwareApplication") %>%
  group_by(name) %>% 
  tally(sort = TRUE)
  
#  summarise(count(name))
```

```{r include = FALSE}
unlink("ropensci.json")
unlink("codemeta.json")
```
---
title: "Validation in JSON-LD"
subtitle: "Validating and consuming in JSON-LD in software development"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
---

## Introduction

Schema validation is a useful and important concept to the distribution of metadata in formats such as XML and JSON, in which the standard-provider creates a schema (specified in an XML-schema, XSD, for XML documents, or [json-schema](http://json-schema.org/) for JSON documents).  Schemas allow us to go beyond the basic notation of making sure a file is simply valid XML or valid JSON, a requirement just to be read in by any parser.  By detailing how the metadata must be structured, what elements must, can, and may not be included, and what data types may be used for those elements, schema help developers consuming the data to anticipate these details and thus build applications which know how to process them.  For the data creator, validation is a convenient way to catch data input errors and ensure a consistent data structure.  

Because schema validation must ensure predictable behavior without knowledge of what any specific application is going to do with the data, it tends to be very strict.  A simple application may not care if certain fields are missing or if integers are mistaken for characters, while to another application these differences could lead it to throw fatal errors.  

The approach of JSON-LD is less prescriptive.  JSON-LD uses the notion of "framing" to let each application specify how it expects it data to be structured.  JSON frames allow each developer consuming the data to handle many of the same issues that schema validation have previously assured.   Readers should consult the [official json-ld framing](https://json-ld.org/spec/latest/json-ld-framing/) documentation for details on this approach.  




```{r message=FALSE}
library(jsonld)
library(jsonlite)
library(magrittr)
library(codemetar)
```


```{r include=FALSE}
knitr::opts_chunk$set(comment="")
if(grepl("windows", tolower(Sys.info()[["sysname"]])))
  knitr::opts_chunk$set(comment="", error =TRUE)
```



## A motivating example:

Consider the following codemeta document:

```{r}
codemeta <- 
'
{
  "@context": "https://doi.org/10.5063/schema/codemeta-2.0",
  "@type": "SoftwareSourceCode",
  "name": "codemetar: Generate CodeMeta Metadata for R Packages",
  "datePublished":" 2017-05-20",
  "version": 1.0,
  "author": [
    {
      "@type": "Person",
      "givenName": "Carl",
      "familyName": "Boettiger",
      "email": "cboettig@gmail.com",
      "@id": "http://orcid.org/0000-0002-1642-628X"
    }],
  "maintainer":  {"@id": "http://orcid.org/0000-0002-1642-628X"}
}
'
```


```{r}
jsonld::jsonld_compact(codemeta, "https://doi.org/10.5063/schema/codemeta-2.0")
```

## Framing: subsetting data

By default, frames return all the input data, while our application may only be interested in some subset.  Often it is sufficient just to ignore these additional terms: in the example above it's just as easy for our application to work with author elements whether or not we have dropped other elements such as `meta$version`.  To restrict a frame to returning only the nodes we explicitly mention, we can use the keyword `@explicit`:



```{r}
frame <- '{
  "@context": "https://doi.org/10.5063/schema/codemeta-2.0",
  "@explicit": "true",
  "@type": "Person",
  "givenName": {},
  "familyName": {}
}'

meta <- 
  jsonld_frame(codemeta, frame)  %>%
  fromJSON(codemeta, simplifyVector = FALSE) %>% 
  getElement("@graph") 

meta[[1]]
```

Note that this has only returned the requested fields in the graph (along with the `@id` and `@type`, which are always included if provided, since they may be required to interpret the data properly). This frame extracts the `givenName` and `familyName` of any `Person` node it finds, regardless of where it occurs, while omitting the rest of the data.  Note that since the frame requests these elements at the top level, they are returned as such, with each match a separate entry in the `@graph`.  Our example has only one person in `meta[[1]]`, had we more matches they would appear in `meta[[2]]`, etc.  Note these returns are un-ordered.


## Framing: expanding node references

The same underlying data can often be expressed in different ways, particularly when dealing with nested data.  Framing can be of great help here to reshape the data into the structure required by the application.  For instance, it would be natural to access the `email` of the `maintainer` in the same manner we did the author, but this fails for our example as `maintainer` is defined only by reference to an ID:

```{r error=TRUE}
meta <- fromJSON(codemeta, simplifyVector = FALSE) 
paste("For complaints, email", meta$maintainer$email)
```

We can confirm that `maintainer` is just an ID:

```{r}
meta$maintainer
```


We can use a frame with the special directive `"@embed": "@always"` to say that we want the full maintainer information embedded an not just referred to by id alone.  Then we can subset `maintainer` just like we do author.  

```{r}
frame <- '{
  "@context": "https://doi.org/10.5063/schema/codemeta-2.0",
  "@embed": "@always"
}'


meta <- 
  jsonld_frame(codemeta, frame) %>%
  fromJSON(codemeta, simplifyVector = FALSE) %>% 
  getElement("@graph") %>% getElement(1)
```



Now we can do


```{r }
paste("For complaints, email", meta$maintainer$email)
```

and see that `email` has been successfully returned from the matching ID under author data.


## Handling unexpected types

JSON-LD routines will simply refuse to compact data if the type differs from what the context expects.  Here is a sample data file
that declares that the `buildInstructions` are included as text, which differs from the `context` file which explicitly states that 
`buildInstructions` should be a URL:

```{r}
codemeta <- 
'
{
  "@context": "https://purl.org/codemeta/2.0",
  "@type": "SoftwareSourceCode",
  "name": "codemetar: Generate CodeMeta Metadata for R Packages",
  "buildInstructions": { 
      "@value": "Just install this package using devtools::install_github", 
      "@type": "Text"
  }
}
'
```

When we perform a framing or compaction operation, `buildInstructions` gets de-referenced to `codemeta:buildInstructions`, because it does not match the context.  This means that if our application asks for `meta$buildInstructions`:

```{r}
meta <-
  jsonld_frame(codemeta, '{"@context": "https://purl.org/codemeta/2.0"}') %>% 
  fromJSON(codemeta, simplifyVector = FALSE) %>%  
  getElement("@graph") %>% getElement(1)    

## above is same as compacting:
#jsonld_compact(codemeta, "https://doi.org/10.5063/schema/codemeta-2.0") %>% 
#  fromJSON(codemeta, simplifyVector = FALSE)

meta$buildInstructions

```

We just get `NULL`, rather than some unexpected type of object (e.g. a string that is not a URL.)  Note that the data is not lost, but simply not dereferenced:


```{r}
names(meta)
meta["codemeta:buildInstructions"]
```


Note that this behavior only happens because the data declared the `"@type": "Text"` explicitly.  JSON-LD algorithms only believe what they are told about type and only look for consistency in declared types.  If you give text but declare it as a `"@type": "URL"`, or don't declare the type at all, JSON-LD algorithms won't know anything is amiss and the property will be compacted as usual.  
`codemetar` can take the path to the source package root to glean as much information as possible.

```{r echo=TRUE, eval=FALSE}
codemetar::write_codemeta()
```

```{r echo=FALSE, eval = identical(Sys.getenv("NOT_CRAN"), "true")}
pkg <- "../.."
codemetar::write_codemeta(pkg = pkg)
```


```{r eval = identical(Sys.getenv("NOT_CRAN"), "true")}
library("magrittr")
"../../codemeta.json" %>%
  details::details(summary = "codemetar's codemeta.json",
                   lang = "json")
```

By default most often from within your package folder you'll simply run `codemetar::write_codemeta()`.

```{r echo = FALSE, results='hide', eval = identical(Sys.getenv("NOT_CRAN"), "true")}
file.remove(file.path(pkg, "codemeta.json"))
```
A new feature is the creation of a minimal schemaorg.json for insertion on your website's webpage for Search Engine Optimization, 
when the `write_minimeta` argument of `write_codemeta()` is `TRUE`.

You could e.g. use the code below in a chunk in README.Rmd with `results="asis"`.

```r
glue::glue('<script type="application/ld+json">
      {glue::glue_collapse(readLines("schemaorg.json"), sep = "\n")}
    </script>')
```

Refer to [Google documentation](https://developers.google.com/search/docs) for more guidance.
**How to keep codemeta.json up-to-date?** In particular, how to keep it up to date with `DESCRIPTION`? `codemetar` itself no longer supports automatic sync, but there are quite a few methods available out there. Choose one that fits well into your workflow!

* You could rely on `devtools::release()` since it will ask you whether you updated codemeta.json when such a file exists.

* You could use a git pre-commit hook that prevents a commit from being done if DESCRIPTION is newer than codemeta.json.

    * You can use the [precommit package](https://github.com/lorenzwalthert/precommit) in which there's a "codemeta-description-updated" hook.
      
    * If that's your only pre-commit hook (i.e. you don't have one created by e.g. `usethis::use_readme_rmd()`), then you can create it using
    
```r
script = readLines(system.file("templates", "description-codemetajson-pre-commit.sh", package = "codemetar"))
usethis::use_git_hook("pre-commit",
                     script = script)
```

* You could use GitHub actions. Refer to GitHub actions docs <https://github.com/features/actions>, and to the example workflow provided in this package (type `system.file("templates", "codemeta-github-actions.yml", package = "codemetar")`). You can use the `cm-skip` keyword in your commit message if you don't want this to run on a specific commit. 
The example workflow provided is setup to only run when a push is made to the master branch. This setup is designed for if you're using a [git flow](https://nvie.com/posts/a-successful-git-branching-model/#the-main-branches) setup where the master branch is only committed and pushed to via pull requests. After each PR merge (and the completion of this GitHub action), your master branch will always be up to date and so long as you don't make manual changes to the codemeta.json file, you won't have merge conflicts.

Alternatively, you can have GitHub actions route run `codemetar` on each commit. If you do this you should try to remember to run `git pull` before making any new changes on your local project. However, if you forgot to pull and already committed new changes, fret not, you can use ([`git pull --rebase`](https://stackoverflow.com/questions/18930527/difference-between-git-pull-and-git-pull-rebase/38139843#38139843)) to rewind you local changes on top of the current upstream `HEAD`. 

```{r, echo = FALSE}
details::details(system.file("templates", "codemeta-github-actions.yml", package = "codemetar"), 
                 summary = "click here to see the workflow",
                 lang = "yaml")
```
For optimal results you need a good internet connection.

The package queries

* `utils::available.packages()` for CRAN and Bioconductor packages;

* GitHub API via the [`gh` package](https://github.com/r-lib/gh), if it finds a GitHub repo URL in DESCRIPTION or as git remote. GitHub API is queried to find the [preferred README](https://developer.github.com/v3/repos/contents/#get-the-readme), and the [repo topics](https://developer.github.com/v3/repos/#list-all-topics-for-a-repository). If you use codemetar for many packages having a [GITHUB_PAT](https://github.com/r-lib/gh#environment-variables) is better;

* [R-hub sysreqs API](https://docs.r-hub.io/#sysreqs) to parse SystemRequirements.

If your machine is offline, a more minimal codemeta.json will be created. If your internet connection is poor or there are firewalls, the codemeta creation might indefinitely hang.
You can install the latest version from CRAN using:

```{r cran-installation, eval = FALSE}
install.packages("codemetar")
```

You can also install the development version of `codemetar` from GitHub with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/codemetar")
```
**Why codemetar?** The 'Codemeta' Project defines a 'JSON-LD' format for describing software metadata, as detailed at <https://codemeta.github.io>.
This package provides utilities to **generate, parse, and modify codemeta.jsonld files automatically for R packages**, as well as tools and examples for **working with codemeta json-ld more generally**.

It has three main goals:

- Quickly **generate a valid codemeta.json file from any valid R package**. To do so, we automatically extract as much metadata as possible using the DESCRIPTION file, as well as extracting metadata from other common best-practices such as the presence of Travis and other badges in README, etc.
- Facilitate the addition of further metadata fields into a codemeta.json file, as well as general manipulation of codemeta files.
- Support the ability to crosswalk between terms used in other metadata standards, as identified by the Codemeta Project Community, see <https://codemeta.github.io/crosswalk/>
The best way to ensure `codemeta.json` is as complete as possible is to set metadata in all the usual places, and then if needed add more metadata.

To ensure you have metadata in the usual places, you can run `codemetar::give_opinions()`.

### Usual terms in DESCRIPTION

* Fill `BugReports` and `URL`. 

* Using the `Authors@R` notation allows a much richer specification of author roles, correct parsing of given vs family names, and email addresses. 

In the current implementation, developers may specify an ORCID url for an author in the optional `comment` field of `Authors@R`, e.g.

```
Authors@R: c(person(given = "Carl",
             family = "Boettiger",
             role = c("aut", "cre", "cph"),
             email = "cboettig@gmail.com",
             comment = c(ORCID = "0000-0002-1642-628X")))
```

which will allow `codemetar` to associate an identifier with the person.  This is clearly something of a hack since R's `person` object lacks an explicit notion of `id`, and may be frowned upon.

### Usual terms in the README

In the README, you can use badges for continuous integration, repo development status (repostatus.org or lifecycle.org), provider ([e.g. for CRAN](https://docs.r-hub.io/#badges)).

### GitHub repo topics

If your package source is hosted on GitHub and there's a way for codemetar to determine that (URL in DESCRIPTION, or git remote URL) codemetar will use [GitHub repo topics](https://docs.github.com/en/github/administering-a-repository/classifying-your-repository-with-topics) as keywords in codemeta.json. If you also set keywords in DESCRIPTION (see next section), codemetar will merge the two lists.

### Set even more terms via DESCRIPTION

In general, setting metadata via the places stated earlier is the best solution because that metadata is used by other tools (e.g. the URLs in DESCRIPTION can help the package users, not only codemetar).

The DESCRIPTION file is the natural place to specify any metadata for an R package.  The `codemetar` package can detect certain additional terms in the [CodeMeta context](https://codemeta.github.io/terms/).  Almost any additional codemeta field can be added to and read from the DESCRIPTION into a `codemeta.json` file (see `codemetar:::additional_codemeta_terms` for a list).  

CRAN requires that you prefix any additional such terms to indicate the use of `schema.org` explicitly, e.g. `keywords` would be specified in a DESCRIPTION file as:

```
X-schema.org-keywords: metadata, codemeta, ropensci, citation, credit, linked-data
```

Where applicable, these will override values otherwise guessed from the source repository.  Use comma-separated lists to separate multiple values to a property, e.g. keywords.  

See the [DESCRIPTION](https://github.com/ropensci/codemetar/blob/master/DESCRIPTION) file of the `codemetar` package for an example.  

### Set the branch that codemetar references

There are a number of places that codemetar will reference a github branch if your code is hosted on github (e.g. for release notes, readme, etc.). By default, codemetar will use the name "master" but you can change that to whatever your default branch is by setting the option "codemeta_branch" (e.g. `options(codemeta_branch = "main")` before calling `write_codemeta()` to use the branch named "main" as the default branch).
**Why bother creating a codemeta.json for your package?** R packages encode lots of metadata in the `DESCRIPTION` file, `README`, and other places, telling users and developers about the package purpose, authors, license, dependencies, and other information that facilitates discovery, adoption, and credit for your software.  Unfortunately, because each software language records this metadata in a different format, that information is hard for search engines, software repositories, and other developers to find and integrate.

By generating a `codemeta.json` file, you turn your metadata into a format that can easily crosswalk between metadata in many other software languages.  CodeMeta is built on [schema.org](https://schema.org) a simple [structured data](https://developers.google.com/search/docs/advanced/structured-data/intro-structured-data) format developed by major search engines like Google and Bing to improve discoverability in search.  CodeMeta is also understood by significant software archiving efforts such as [Software Heritage](https://www.softwareheritage.org/) Project, which seeks to permanently archive all open source software.

For more general information about the CodeMeta Project for defining software metadata, see <https://codemeta.github.io>.  In particular, new users might want to start with the [User Guide](https://codemeta.github.io/user-guide/), while those looking to learn more about JSON-LD and consuming existing codemeta files should see the [Developer Guide](https://codemeta.github.io/developer-guide/).   
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/give_opinions.R
\name{give_opinions}
\alias{give_opinions}
\title{Function giving opinions about a package}
\usage{
give_opinions(pkg_path = getwd(), verbose = FALSE)
}
\arguments{
\item{pkg_path}{Path to the package root}

\item{verbose}{Whether to print message related to internet download progress.}
}
\value{
A data.frame of opinions
}
\description{
Function giving opinions about a package
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_badges.R
\name{extract_badges}
\alias{extract_badges}
\title{Extract all badges from Markdown file}
\usage{
extract_badges(path)
}
\arguments{
\item{path}{Path to Markdown file}
}
\value{
A data.frame with for each badge its text, link and link to
its image.
}
\description{
Extract all badges from Markdown file
}
\examples{
\dontrun{
extract_badges(system.file("examples/README_fakepackage.md", package="codemetar"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codemetar-package.R
\docType{package}
\name{codemetar-package}
\alias{codemetar}
\alias{codemetar-package}
\title{codemetar: generate codemeta metadata for R packages}
\description{
The 'Codemeta' Project defines a 'JSON-LD' format for describing software metadata, as detailed at <https://codemeta.github.io>. This package provides utilities to generate, parse, and modify 'codemeta.json' files automatically for R packages, as well as tools and examples for working with 'codemeta.json' 'JSON-LD' more generally.
}
\details{
\strong{Why bother creating a codemeta.json for your package?} R packages
encode lots of metadata in the \code{DESCRIPTION} file, \code{README}, and other
places, telling users and developers about the package purpose, authors,
license, dependencies, and other information that facilitates discovery,
adoption, and credit for your software. Unfortunately, because each
software language records this metadata in a different format, that
information is hard for search engines, software repositories, and other
developers to find and integrate.

By generating a \code{codemeta.json} file, you turn your metadata into a
format that can easily crosswalk between metadata in many other software
languages. CodeMeta is built on \href{https://schema.org}{schema.org} a
simple \href{https://developers.google.com/search/docs/advanced/structured-data/intro-structured-data}{structured data}
format developed by major search engines like Google and Bing to improve
discoverability in search. CodeMeta is also understood by significant
software archiving efforts such as \href{https://www.softwareheritage.org/}{Software Heritage} Project, which seeks to
permanently archive all open source software.

For more general information about the CodeMeta Project for defining
software metadata, see \url{https://codemeta.github.io}. In particular, new
users might want to start with the \href{https://codemeta.github.io/user-guide/}{User Guide}, while those looking to
learn more about JSON-LD and consuming existing codemeta files should
see the \href{https://codemeta.github.io/developer-guide/}{Developer Guide}.

\strong{Why codemetar?} The ‘Codemeta’ Project defines a ‘JSON-LD’ format for
describing software metadata, as detailed at
\url{https://codemeta.github.io}. This package provides utilities to
\strong{generate, parse, and modify codemeta.jsonld files automatically for R
packages}, as well as tools and examples for \strong{working with codemeta
json-ld more generally}.

It has three main goals:
\itemize{
\item Quickly \strong{generate a valid codemeta.json file from any valid R
package}. To do so, we automatically extract as much metadata as
possible using the DESCRIPTION file, as well as extracting metadata
from other common best-practices such as the presence of Travis and
other badges in README, etc.
\item Facilitate the addition of further metadata fields into a
codemeta.json file, as well as general manipulation of codemeta
files.
\item Support the ability to crosswalk between terms used in other
metadata standards, as identified by the Codemeta Project Community,
see \url{https://codemeta.github.io/crosswalk/}
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/codemetar}
  \item \url{https://docs.ropensci.org/codemetar/}
  \item Report bugs at \url{https://github.com/ropensci/codemetar/issues}
}

}
\author{
\strong{Maintainer}: Carl Boettiger \email{cboettig@gmail.com} (\href{https://orcid.org/0000-0002-1642-628X}{ORCID}) [copyright holder]

Authors:
\itemize{
  \item Maëlle Salmon (\href{https://orcid.org/0000-0002-2815-0399}{ORCID}) [contributor]
}

Other contributors:
\itemize{
  \item Anna Krystalli (\href{https://orcid.org/0000-0002-2378-4915}{ORCID}) [reviewer, contributor]
  \item Toph Allen (\href{https://orcid.org/0000-0003-4580-091X}{ORCID}) [reviewer]
  \item rOpenSci (https://ropensci.org/) [funder]
  \item Katrin Leinweber (\href{https://orcid.org/0000-0001-5135-5758}{ORCID}) [contributor]
  \item Noam Ross (\href{https://orcid.org/0000-0002-2136-0000}{ORCID}) [contributor]
  \item Arfon Smith [contributor]
  \item Jeroen Ooms (\href{https://orcid.org/0000-0002-4035-0289}{ORCID}) [contributor]
  \item Sebastian Meyer (\href{https://orcid.org/0000-0002-1791-9449}{ORCID}) [contributor]
  \item Michael Rustler (\href{https://orcid.org/0000-0003-0647-7726}{ORCID}) [contributor]
  \item Hauke Sonnenberg (\href{https://orcid.org/0000-0001-9134-2871}{ORCID}) [contributor]
  \item Sebastian Kreutzer (\href{https://orcid.org/0000-0002-0734-2199}{ORCID}) [contributor]
  \item Thierry Onkelinx (\href{https://orcid.org/0000-0001-8804-4216}{ORCID}) [contributor]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_codemeta.R
\name{write_codemeta}
\alias{write_codemeta}
\title{write_codemeta}
\usage{
write_codemeta(
  pkg = ".",
  path = "codemeta.json",
  root = ".",
  id = NULL,
  use_filesize = TRUE,
  force_update = getOption("codemeta_force_update", TRUE),
  use_git_hook = NULL,
  verbose = TRUE,
  write_minimeta = FALSE,
  ...
)
}
\arguments{
\item{pkg}{package path to package root, or package name, or description file
(character), or a codemeta object (list)}

\item{path}{file name of the output, leave at default "codemeta.json"}

\item{root}{if pkg is a codemeta object, optionally give the path to package
root. Default guess is current dir.}

\item{id}{identifier for the package, e.g. a DOI (or other resolvable URL)}

\item{use_filesize}{whether to try to estimating and adding a filesize by using
\code{base::file.size()}. Files in \code{.Rbuildignore} are ignored.}

\item{force_update}{Update guessed fields even if they are defined in an
existing codemeta.json file}

\item{use_git_hook}{Deprecated argument.}

\item{verbose}{Whether to print messages indicating opinions e.g. when
DESCRIPTION has no URL. -- See \code{\link{give_opinions}};
and indicating the progress of internet downloads.}

\item{write_minimeta}{whether to also create the file schemaorg.json that
corresponds to the metadata Google would validate, to be inserted to a
webpage for SEO. It is saved as "inst/schemaorg.json" alongside \code{path} (by
default, "codemeta.json").}

\item{...}{additional arguments to \code{\link{write_json}}}
}
\value{
writes out the codemeta.json file, and schemaorg.json if \code{write_codemeta}
is \code{TRUE}.
}
\description{
write out a codemeta.json file for a given package.  This function is
basically a wrapper around create_codemeta() to both create the codemeta
object and write it out to a JSON-LD-formatted file in one command. It can
also be used simply to write out to JSON-LD any existing object created with
\code{create_codemeta()}.
}
\details{
\strong{Why bother creating a codemeta.json for your package?} R packages
encode lots of metadata in the \code{DESCRIPTION} file, \code{README}, and other
places, telling users and developers about the package purpose, authors,
license, dependencies, and other information that facilitates discovery,
adoption, and credit for your software. Unfortunately, because each
software language records this metadata in a different format, that
information is hard for search engines, software repositories, and other
developers to find and integrate.

By generating a \code{codemeta.json} file, you turn your metadata into a
format that can easily crosswalk between metadata in many other software
languages. CodeMeta is built on \href{https://schema.org}{schema.org} a
simple \href{https://developers.google.com/search/docs/advanced/structured-data/intro-structured-data}{structured data}
format developed by major search engines like Google and Bing to improve
discoverability in search. CodeMeta is also understood by significant
software archiving efforts such as \href{https://www.softwareheritage.org/}{Software Heritage} Project, which seeks to
permanently archive all open source software.

For more general information about the CodeMeta Project for defining
software metadata, see \url{https://codemeta.github.io}. In particular, new
users might want to start with the \href{https://codemeta.github.io/user-guide/}{User Guide}, while those looking to
learn more about JSON-LD and consuming existing codemeta files should
see the \href{https://codemeta.github.io/developer-guide/}{Developer Guide}.

\strong{How to keep codemeta.json up-to-date?} In particular, how to keep it
up to date with \code{DESCRIPTION}? \code{codemetar} itself no longer supports
automatic sync, but there are quite a few methods available out there.
Choose one that fits well into your workflow!
\itemize{
\item You could rely on \code{devtools::release()} since it will ask you
whether you updated codemeta.json when such a file exists.
\item You could use a git pre-commit hook that prevents a commit from
being done if DESCRIPTION is newer than codemeta.json.
\itemize{
\item You can use the \href{https://github.com/lorenzwalthert/precommit}{precommit package} in which
there’s a “codemeta-description-updated” hook.
\item If that’s your only pre-commit hook (i.e. you don’t have one
created by e.g. \code{usethis::use_readme_rmd()}), then you can
create it using
}
}\if{html}{\out{<div class="sourceCode r">}}\preformatted{script = readLines(system.file("templates", "description-codemetajson-pre-commit.sh", package = "codemetar"))
usethis::use_git_hook("pre-commit",
                     script = script)
}\if{html}{\out{</div>}}
\itemize{
\item You could use GitHub actions. Refer to GitHub actions docs
\url{https://github.com/features/actions}, and to the example workflow
provided in this package (type
\code{system.file("templates", "codemeta-github-actions.yml", package = "codemetar")}).
You can use the \code{cm-skip} keyword in your commit message if you
don’t want this to run on a specific commit. The example workflow
provided is setup to only run when a push is made to the master
branch. This setup is designed for if you’re using a \href{https://nvie.com/posts/a-successful-git-branching-model/#the-main-branches}{git flow}
setup where the master branch is only committed and pushed to via
pull requests. After each PR merge (and the completion of this
GitHub action), your master branch will always be up to date and so
long as you don’t make manual changes to the codemeta.json file, you
won’t have merge conflicts.
}

Alternatively, you can have GitHub actions route run \code{codemetar} on each
commit. If you do this you should try to remember to run \verb{git pull}
before making any new changes on your local project. However, if you
forgot to pull and already committed new changes, fret not, you can use
(\href{https://stackoverflow.com/questions/18930527/difference-between-git-pull-and-git-pull-rebase/38139843#38139843}{\verb{git pull --rebase}})
to rewind you local changes on top of the current upstream \code{HEAD}.\if{html}{\out{
<details closed>
<summary>
<span title="Click to Expand"> click here to see the workflow </span>
</summary>
}}
\if{html}{\out{<div class="sourceCode yaml">}}\preformatted{on:
  push:
    branches: master
    paths:
      - DESCRIPTION
      - .github/workflows/main.yml

name: Render codemeta
jobs:
  render:
    name: Render codemeta
    runs-on: macOS-latest
    if: "!contains(github.event.head_commit.message, 'cm-skip')"
    steps:
      - uses: actions/checkout@v1
      - uses: r-lib/actions/setup-r@v1
      - name: Install codemetar
        run: Rscript -e 'install.packages("codemetar")'
      - name: Render codemeta
        run: Rscript -e 'codemetar::write_codemeta()'
      - name: Commit results
        run: |
          git commit codemeta.json -m 'Re-build codemeta.json' || echo "No changes to commit"
          git push https://$\{\{github.actor\}\}:$\{\{secrets.GITHUB_TOKEN\}\}@github.com/$\{\{github.repository\}\}.git HEAD:$\{\{ github.ref \}\} || echo "No changes to commit"
}\if{html}{\out{</div>}}\if{html}{\out{
</details>
}}
\if{html}{\out{
<br>
}}
}
\section{Technical details}{

If pkg is a codemeta object, the function will attempt to update any
fields it can guess (i.e. from the DESCRIPTION file), overwriting any
existing data in that block. In this case, the package root directory
should be the current working directory.

When creating and writing a codemeta.json for the first time, the function
adds "codemeta.json" to .Rbuildignore.
}

\examples{
\dontrun{
# from anywhere in the package source directory
write_codemeta()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \verb{magrittr::[\\\%>\\\%][magrittr::\\\%>\\\%]} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_codemeta.R
\name{create_codemeta}
\alias{create_codemeta}
\title{create_codemeta}
\usage{
create_codemeta(
  pkg = ".",
  root = ".",
  id = NULL,
  use_filesize = FALSE,
  force_update = getOption("codemeta_force_update", TRUE),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{pkg}{package path to package root, or package name, or description file
(character), or a codemeta object (list)}

\item{root}{if pkg is a codemeta object, optionally give the path to package
root. Default guess is current dir.}

\item{id}{identifier for the package, e.g. a DOI (or other resolvable URL)}

\item{use_filesize}{whether to try to estimating and adding a filesize by using
\code{base::file.size()}. Files in \code{.Rbuildignore} are ignored.}

\item{force_update}{Update guessed fields even if they are defined in an
existing codemeta.json file}

\item{verbose}{Whether to print messages indicating opinions e.g. when
DESCRIPTION has no URL. -- See \code{\link{give_opinions}};
and indicating the progress of internet downloads.}

\item{...}{additional arguments to \code{\link{write_json}}}
}
\value{
a codemeta list object
}
\description{
create a codemeta list object in R for further manipulation. Similar
to \code{\link[=write_codemeta]{write_codemeta()}}, but returns an R list object rather
than writing directly to a file.  See examples.
}
\examples{
\donttest{
path <- system.file("", package="codemeta")
cm <- create_codemeta(path)
cm$keywords <- list("metadata", "ropensci")
}
}
