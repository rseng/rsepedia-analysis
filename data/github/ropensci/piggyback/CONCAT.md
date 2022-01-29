---
title: 'Managing Larger Data on a GitHub Repository'
tags:
 - data
 - GitHub
authors:
 - name: Carl Boettiger
   orcid: 0000-0002-1642-628X
   affiliation: 1
affiliations:
 - name: University of California, Berkeley
   index: 1
date: 2018-05-30
bibliography: paper.bib
---

# Piggyback: Working with larger data in GitHub

GitHub has become a central component for preserving and sharing software-driven
analysis in academic research [@Ram2013].  As scientists adopt this workflow, 
a desire to manage data associated with the analysis in the same manner soon emerges.
While small data can easily be committed to GitHub repositories along-side source
code and analysis scripts, files larger than 50 MB cannot. Existing work-arounds introduce
significant complexity and break the ease of sharing [@Boettiger2018]. 

This package provides a simple work-around by allowing larger
(up to 2 GB) data files to piggyback on a repository as assets attached to individual
GitHub releases. `piggyback` provides a workflow similar to Git LFS [@GitLFS], in 
which data files can be tracked by type and pushed and pulled to GitHub with dedicated
commands. These files are not handled by git in any way, but instead are
uploaded, downloaded, or edited directly by calls through the GitHub API [@API3]. These
data files can be versioned manually by creating different releases.  This approach
works equally well with public or private repositories.  Data can be uploaded
and downloaded programmatically from scripts. No authentication is required to
download data from public repositories.

## Examples

As long as a repository has at least one release, users can upload a set of specified
files from the current repository to that release by simply passing the file names to
`pb_upload()`.  Specify individual files to download using `pb_download()`, or use no
arguments to download all data files attached to the latest release.  Alternatively,
users can track files by a given pattern: for instance, `pb_track("*.csv")` will track
all `*.csv` files in the repository.  Then use `pb_upload(pb_track())` to upload all
currently tracked files.  `piggyback` compares timestamps to avoid unnecessary transfer.
The `piggyback` package looks for the same `GITHUB_TOKEN` environmental variable for 
authentication that is used across GitHub APIs. Details are provided in an introductory
vignette [@Boettiger2018b].

# References

<!-- README.md is generated from README.Rmd. Please edit that file -->

# piggyback <img src="man/figures/logo.svg" align="right" alt="" width="120" />

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/ropensci/piggyback/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/piggyback/actions)
[![Coverage
status](https://codecov.io/gh/ropensci/piggyback/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/piggyback?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/piggyback)](https://cran.r-project.org/package=piggyback)
[![Peer Review
Status](https://badges.ropensci.org/220_status.svg)](https://github.com/ropensci/software-review/issues/220)
[![DOI](https://zenodo.org/badge/132979724.svg)](https://zenodo.org/badge/latestdoi/132979724)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00971/status.svg)](https://doi.org/10.21105/joss.00971)
<!-- badges: end -->

Because larger (&gt; 50 MB) data files cannot easily be committed to
git, a different approach is required to manage data associated with an
analysis in a GitHub repository. This package provides a simple
work-around by allowing larger ([up to 2 GB per
file](https://docs.github.com/en/github/managing-large-files/distributing-large-binaries))
data files to piggyback on a repository as assets attached to individual
GitHub releases. These files are not handled by git in any way, but
instead are uploaded, downloaded, or edited directly by calls through
the GitHub API. These data files can be versioned manually by creating
different releases. This approach works equally well with public or
private repositories. Data can be uploaded and downloaded
programmatically from scripts. No authentication is required to download
data from public repositories.

## Installation

Install from CRAN via

``` r
install.packages("piggyback")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/piggyback")
```

## Quickstart

See the [piggyback
vignette](https://docs.ropensci.org/piggyback/articles/intro.html) for
details on authentication and additional package functionality.

Piggyback can download data attached to a release on any repository:

``` r
library(piggyback)
pb_download("iris.tsv.gz", repo = "cboettig/piggyback-tests", dest = tempdir())
#> Warning in pb_download("iris.tsv.gz", repo = "cboettig/piggyback-tests", :
#> file(s) iris.tsv.gz not found in repo cboettig/piggyback-tests
```

Downloading from private repos or uploading to any repo requires
authentication, so be sure to set a `GITHUB_TOKEN` (or `GITHUB_PAT`)
environmental variable, or include the `.token` argument. Omit the file
name to download all attached objects. Omit the repository name to
default to the current repository. See [introductory
vignette](https://docs.ropensci.org/piggyback/articles/intro.html) or
function documentation for details.

We can also upload data to any existing release (defaults to `latest`):

``` r
## We'll need some example data first.
## Pro tip: compress your tabular data to save space & speed upload/downloads
readr::write_tsv(mtcars, "mtcars.tsv.gz")

pb_upload("mtcars.tsv.gz", repo = "cboettig/piggyback-tests")
```

## Git LFS and other alternatives

`piggyback` acts like a poor soul’s [Git
LFS](https://git-lfs.github.com/). Git LFS is not only expensive, it
also [breaks GitHub’s collaborative
model](https://angryfrenchman.org/github-s-large-file-storage-is-no-panacea-for-open-source-quite-the-opposite-12c0e16a9a91)
– basically if someone wants to submit a PR with a simple edit to your
docs, they cannot fork your repository since that would otherwise count
against your Git LFS storage. Unlike Git LFS, `piggyback` doesn’t take
over your standard `git` client, it just perches comfortably on the
shoulders of your existing GitHub API. Data can be versioned by
`piggyback`, but relative to `git LFS` versioning is less strict:
uploads can be set as a new version or allowed to overwrite previously
uploaded data.

## But what will GitHub think of this?

[GitHub
documentation](https://docs.github.com/en/github/managing-large-files/distributing-large-binaries)
at the time of writing endorses the use of attachments to releases as a
solution for distributing large files as part of your project:

![](man/figures/github-policy.png)

Of course, it will be up to GitHub to decide if this use of release
attachments is acceptable in the long term.

<!--
 When GitHub first came online, it was questioned whether committing binary objects and data to GitHub was acceptable or an abuse of a *source code* repository.  GitHub has since clearly embraced a inclusive notion of "repository" for containing far more than pure source.  I believe attaching data that is essential to replicating an analysis and within the 2 GB file limits enforced by GitHub to be in the same spirit of this inclusive notion, but GitHub may decide otherwise. 
 -->

Also see our [vignette comparing
alternatives](https://docs.ropensci.org/piggyback/articles/alternatives.html).

------------------------------------------------------------------------

Please note that this project is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By participating in
this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# piggyback 0.1.1

* switch to gh::gh_token() for token management.  Still supports the same env var approach, but also compatible with `gitcreds` and other use.
* resolve issue in `pb_upload()` when creating a new tag in the process, previously data would be attached to the previously `latest` tag instead of the newly created one. 
* resolve issue in `pb_download()` where httr would report a 401 status even after data successfully downloads. 

# piggyback 0.1.0

* address remaining authentication issue in changes to GitHub API (on pb_upload()) [#47]
* Use flat file structure on upload/download instead of encoding path [#48]
* improve performance via more aggressive memoising of `pb_info()` calls, inceasing default `piggyback_cache_duration` to 10 minutes [#46]
* Resolve bug introduced by API changes that would stop creation of tags on repos with default branch called `main` or without previous releases [#48]


# piggyback 0.0.12

* address issues in authentication due to changes in GitHub API (#37)

# piggyback 0.0.11 2020-02-25

* `guess_repo()` now infers a remote when there are multiple associated with the repo. The "upstream" (preferred) or "origin" repo is selected if either exists, otherwise the function errors and asks the user to explicitly specify a repo (#31).
* `release_info()` now works properly when there are no existing releases, which enables the usage of `pb_new_release()` on repos without a release (#29).
* Fix error on `pb_info()` under certain cases which resulted in `Error in a[[1]] : subscript out of bounds`, (#36)
* Fix CRAN unit-test on deleting file

# piggyback 0.0.10 2018-02-06

* Improve interface regarding `overwrite` behavior in `pb_upload()` (#25)
* Bugfixes for errors introduced in 0.0.9: 
   - Access all assets on a release instead of first 30.  This could break upload and download. (#23, #24)
   - Uploading of directory paths could cause download errors in `pb_download()`. (#24, #26)

# piggyback 0.0.9, 2019-01-08

* Enable re-upload and deletion of partially uploaded files (#19)

# piggyback 0.0.8, 2018-10-06

* Updates to documentation, streamlining tests
* remove dependency on `utils::askYesNo` which is only available in R >= 3.5.0

# piggyback 0.0.7, 2018-09-30

* Initial release to CRAN

--------------------------------------------

# piggyback 0.0.6, 2018-09-21

* bugfix for migrating unit test

# piggyback 0.0.6, 2018-09-21

* bugfix for migrating unit test, JOSS submission

# piggyback 0.0.5, 2018-09-21

* initial Onboarding to rOpenSci

# piggyback 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
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
Dear CRAN Maintainers,

This release fixes a few minor bugs due to changes in the API.  See NEWS.  

- Carl


