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


---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  eval = TRUE,
  collapse = TRUE,
  message = FALSE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# piggyback <img src="man/figures/logo.svg" align="right" alt="" width="120" />

  <!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/ropensci/piggyback/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/piggyback/actions)
[![Coverage status](https://codecov.io/gh/ropensci/piggyback/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/piggyback?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/piggyback)](https://cran.r-project.org/package=piggyback)
[![Peer Review Status](https://badges.ropensci.org/220_status.svg)](https://github.com/ropensci/software-review/issues/220)
[![DOI](https://zenodo.org/badge/132979724.svg)](https://zenodo.org/badge/latestdoi/132979724)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00971/status.svg)](https://doi.org/10.21105/joss.00971)
  <!-- badges: end -->


Because larger (> 50 MB) data files cannot easily be committed to git, a different approach is required to manage data associated with an analysis in a GitHub repository.  This package provides a simple work-around by allowing larger ([up to 2 GB per file](https://docs.github.com/en/github/managing-large-files/distributing-large-binaries)) data files to piggyback on a repository as assets attached to individual GitHub releases.  These files are not handled by git in any way, but instead are uploaded, downloaded, or edited directly by calls through the GitHub API. These data files can be versioned manually by creating different releases.  This approach works equally well with public or private repositories.  Data can be uploaded and downloaded programmatically from scripts. No authentication is required to download data from public repositories.




## Installation


Install from CRAN via

``` r
install.packages("piggyback")
```

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/piggyback")
```

## Quickstart

See the [piggyback vignette](https://docs.ropensci.org/piggyback/articles/intro.html) for details on authentication and additional package functionality. 

Piggyback can download data attached to a release on any repository: 

```{r results="hide"}
library(piggyback)
pb_download("iris.tsv.gz", repo = "cboettig/piggyback-tests", dest = tempdir())
```


Downloading from private repos or uploading to any repo requires authentication, so be sure to set a `GITHUB_TOKEN` (or `GITHUB_PAT`) environmental variable, or include the `.token` argument.  Omit the file name to download all attached objects. Omit the repository name to default to the current repository.  See [introductory vignette](https://docs.ropensci.org/piggyback/articles/intro.html) or function documentation for details.  

We can also upload data to any existing release (defaults to `latest`):

```{r eval=FALSE}
## We'll need some example data first.
## Pro tip: compress your tabular data to save space & speed upload/downloads
readr::write_tsv(mtcars, "mtcars.tsv.gz")

pb_upload("mtcars.tsv.gz", repo = "cboettig/piggyback-tests")
```




## Git LFS and other alternatives


`piggyback` acts like a poor soul's [Git LFS](https://git-lfs.github.com/). Git LFS is not only expensive, it also [breaks GitHub's collaborative model](https://angryfrenchman.org/github-s-large-file-storage-is-no-panacea-for-open-source-quite-the-opposite-12c0e16a9a91) -- basically if someone wants to submit a PR with a simple edit to your docs, they cannot fork your repository since that would otherwise count against your Git LFS storage.   Unlike Git LFS, `piggyback` doesn't take over your standard `git` client, it just perches comfortably on the shoulders of your existing GitHub API.  Data can be versioned by `piggyback`, but relative to `git LFS` versioning is less strict: uploads can be set as a new version or allowed to overwrite previously uploaded data.  

## But what will GitHub think of this?

[GitHub documentation](https://docs.github.com/en/github/managing-large-files/distributing-large-binaries) at the time of writing endorses the use of attachments to releases as a solution for distributing large files as part of your project:

![](man/figures/github-policy.png)


Of course, it will be up to GitHub to decide if this use of release attachments is acceptable in the long term. 

<!--
 When GitHub first came online, it was questioned whether committing binary objects and data to GitHub was acceptable or an abuse of a *source code* repository.  GitHub has since clearly embraced a inclusive notion of "repository" for containing far more than pure source.  I believe attaching data that is essential to replicating an analysis and within the 2 GB file limits enforced by GitHub to be in the same spirit of this inclusive notion, but GitHub may decide otherwise. 
 -->

Also see our [vignette comparing alternatives](https://docs.ropensci.org/piggyback/articles/alternatives.html).

----

Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/).
By participating in this project you agree to abide by its terms.

```{r include=FALSE}
unlink("*.gz")
codemeta::write_codemeta()
```


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Piggyback Data atop your GitHub Repository!"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{piggyback}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  results="hide",
  eval = Sys.getenv("CBOETTIG_TOKEN", FALSE)
)

Sys.setenv(piggyback_cache_duration=0)

```



# Why `piggyback`?

`piggyback` grew out of the needs of students both in my classroom and in my research group, who frequently need to work with data files somewhat larger than one can conveniently manage by committing directly to GitHub.  As we frequently want to share and run code that depends on >50MB data files on each of our own machines, on continuous integration (i.e. [travis](https://travis-ci.org)), and on larger computational servers, data sharing quickly becomes a bottleneck. 

[GitHub allows](https://docs.github.com/en/github/managing-large-files/distributing-large-binaries) repositories to attach files of up to 2 GB each to releases as a way to distribute large files associated with the project source code.  There is no limit on the number of files or bandwidth to deliver them.  

## Installation

Install the latest release from CRAN using:

``` r
install.packages("piggyback")
```

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/piggyback")
```

## Authentication

No authentication is required to download data from *public* GitHub repositories using `piggyback`. Nevertheless, `piggyback` recommends setting a token when possible to avoid rate limits. To upload data to any repository, or to download data from *private* repositories, you will need to authenticate first. 

To do so, add your [GitHub Token](https://github.com/settings/tokens/new?scopes=repo,gist&description=R:GITHUB_PAT) to an environmental variable, e.g. in a `.Renviron` file in your home directory or project directory (any private place you won't upload), see `usethis::edit_r_environ()`.  For one-off use you can also set your token from the R console using:

```r
Sys.setenv(GITHUB_TOKEN="xxxxxx")
```

But try to avoid putting `Sys.setenv()` in  any R scripts -- remember, the goal here is to avoid writing your private token in any file that might be shared, even privately.   For more help setting up a GitHub token, for the first time, see `usethis::browse_github_pat()`.


## Downloading data

Download the latest version or a specific version of the data:

```r
library(piggyback)
```

```r
pb_download("iris2.tsv.gz", 
            repo = "cboettig/piggyback-tests",
            tag = "v0.0.1",
            dest = tempdir())
```

**Note**: Whenever you are working from a location inside a git repository corresponding to your GitHub repo, you can simply omit the `repo` argument and it will be detected automatically. Likewise, if you omit the release `tag`, `the`pb_download` will simply pull data from most recent release (`latest`).  Third, you can omit `tempdir()` if you are using an RStudio Project (`.Rproj` file) in your repository, and then the download location will be relative to Project root.  `tempdir()` is used throughout the examples only to meet CRAN policies and is unlikely to be the choice you actually want here.  


Lastly, simply omit the file name to download all assets connected with a given release.  

```r
pb_download(repo = "cboettig/piggyback-tests",
            tag = "v0.0.1",
            dest = tempdir())
```  

These defaults mean that in most cases, it is sufficient to simply call `pb_download()` without additional arguments to pull in any data associated with a project on a GitHub repo that is too large to commit to git directly. 

`pb_download()` will skip the download of any file that already exists locally if the timestamp on the local copy is more recent than the timestamp on the GitHub copy.  `pb_download()` also includes arguments to control the timestamp behavior, progress bar, whether existing files should be overwritten, or if any particular files should not be downloaded.  See function documentation for details.  


Sometimes it is preferable to have a URL from which the data can be read in directly, rather than downloading the data to a local file.  For example, such a URL can be embedded directly into another R script, avoiding any dependence on `piggyback` (provided the repository is already public.)  To get a list of URLs rather than actually downloading the files, use `pb_download_url()`:

```r
pb_download_url("data/mtcars.tsv.gz", 
                repo = "cboettig/piggyback-tests",  
                tag = "v0.0.1") 
```



## Uploading data

If your GitHub repository doesn't have any [releases](https://docs.github.com/en/github/administering-a-repository/managing-releases-in-a-repository) yet, `piggyback` will help you quickly create one.  Create new releases to manage multiple versions of a given data file. While you can create releases as often as you like, making a new release is by no means necessary each time you upload a file.  If maintaining old versions of the data is not useful, you can stick with a single release and upload all of your data there.  

```r
pb_new_release("cboettig/piggyback-tests", "v0.0.2")
```

Once we have at least one release available, we are ready to upload.  By default, `pb_upload` will attach data to the latest release.  

```r
## We'll need some example data first.
## Pro tip: compress your tabular data to save space & speed upload/downloads
readr::write_tsv(mtcars, "mtcars.tsv.gz")

pb_upload("mtcars.tsv.gz", 
          repo = "cboettig/piggyback-tests", 
          tag = "v0.0.1")
```

Like `pb_download()`, `pb_upload()` will overwrite any file of the same name already attached to the release file by default, unless the timestamp the previously uploaded version is more recent.  You can toggle these settings with `overwrite=FALSE` and `use_timestamps=FALSE`.  


## Additional convenience functions

List all files currently piggybacking on a given release.  Omit the `tag` to see files on all releases.  


```r
pb_list(repo = "cboettig/piggyback-tests", 
        tag = "v0.0.1")
```

Delete a file from a release:

```r
pb_delete(file = "mtcars.tsv.gz", 
          repo = "cboettig/piggyback-tests", 
          tag = "v0.0.1")
```

Note that this is irreversible unless you have a copy of the data elsewhere. 

## git-style tracking

`piggyback` can be used in a Git-LFS-like manner by tracking all files that match a particular pattern, typically a file extension such as `*.tif` or `*.tar.gz` frequently found on large binary data files associated with a project but too big to commit to git. Similarly, specific directories for data files can be tracked.  `pb_track()` function takes such patterns and stores them in into a hidden config file, `.pbattributes` (just like `.gitattributes` in Git LFS, which you can also edit manually).    

```r
pb_track(c("*.tsv.gz", "*.tif", "*.zip"))
pb_track("data/*")
```

Adding a pattern with `pb_track()` will also automatically add that pattern to `.gitignore`, since these data files will be piggybacking on top of the repo rather than being version managed by `git`.  You probably will want to check in the `.pbattributes` file to version control, just as you would a `.gitattributes` or `.gitignore`. 

Once you have tracked certain file types, it is easy to push all such files up to GitHub by piping `pb_track() %>% pb_upload()`.  `pb_track()` just returns file paths to all matching files.  As usual, this can upload to a specific repository and tag or merely to the defaults.  

```r
library(magrittr)
pb_track() %>% pb_upload(repo = "cboettig/piggyback-tests", tag = "v0.0.1")
```

Similarly, you can download all current data assets of the latest or specified release by using `pb_download()` with no arguments.


## Caching
 

To reduce API calls to GitHub, piggyback caches most calls with a timeout of 1 second by default.  This avoids repeating identical requests to update it's internal record of the repository data (releases, assets, timestamps, etc) during programmatic use.  You can increase or decrease this delay by setting the environmental variable in seconds, e.g. `Sys.setenv("piggyback_cache_duration"=10)` for a longer delay or `Sys.setenv("piggyback_cache_duration"=0)` to disable caching. 


## Path names

GitHub assets attached to a release do not support file paths, and will convert most special characters (`#`, `%`, etc) to `.` or throw an error (e.g. for file names containing `$`, `@`, `/`).  To preserve path information on uploading data, `piggyback` uses relative paths (relative to the working directory, or for `pb_push()` and `pb_pull`, relative to the project directory, see `here::here()`) in data file names, and encodes the system path delimiter as `.2f` (`%2f` is the HTML encoding of a literal `/`, but `%` cannot be used in asset names).  `piggyback` functions will always show and use the decoded file names, e.g. `data/mtcars.csv`, but you'll see `data.2fmtcars.csv` if you look at the release attachment on GitHub.


## A Note on GitHub Releases vs Data Archiving

`piggyback` is not intended as a data archiving solution.  Importantly, bear in mind that there is nothing special about multiple "versions" in releases, as far as data assets uploaded by `piggyback` are concerned.  The data files `piggyback` attaches to a Release can be deleted or modified at any time -- creating a new release to store data assets is the functional equivalent of just creating new directories `v0.1`, `v0.2` to store your data.  (GitHub Releases are always pinned to a particular `git` tag, so the code/git-managed contents associated with repo are more immutable, but remember our data assets just piggyback on top of the repo).  

Permanent, published data should always be archived in a proper data repository with a DOI, such as [zenodo.org](https://zenodo.org). Zenodo can freely archive public research data files up to 50 GB in size, and data is strictly versioned (once released, a DOI always refers to the same version of the data, new releases are given new DOIs). `piggyback` is meant only to lower the friction of working with data during the research process.  (e.g. provide data accessible to collaborators or continuous integration systems during research process, including for private repositories.)



## What will GitHub think of this?

[GitHub documentation](https://docs.github.com/en/github/managing-large-files/distributing-large-binaries) at the time of writing endorses the use of attachments to releases as a solution for distributing large files as part of your project:

![](https://github.com/ropensci/piggyback/raw/83776863b34bb1c9962154608a5af41867a0622f/man/figures/github-policy.png)


Of course, it will be up to GitHub to decide if this use of release attachments is acceptable in the long term. 





---
title: "Piggyback comparison to alternatives"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{alternatives}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
## `piggyback` vs the alternatives

There are many alternatives to `piggyback`, and after considerable experience I haven't found any that ticked all the boxes for me:

- [x] Free storage
- [x] Can be integrated into private code / private workflows
- [x] Simple and practical to deploy on continuous integration
- [x] Works well with private data
- [x] Minimal configuration



### Git LFS

Git LFS provides the closest user experience to what I was going for. It stands out above all other alternatives for providing both the *best authentication* experience (relying directly on any of the standard `git` authentication mechanisms such as https, ssh keys, app integration), and it provides the most legitimate version control of the data.  However, there are many show-stoppers to using Git LFS for me.

- GitHub pricing & resulting problems for GitHub's fork /  PR model.  [Described eloquently here](https://angryfrenchman.org/github-s-large-file-storage-is-no-panacea-for-open-source-quite-the-opposite-12c0e16a9a91).  Basically, despite generous rates and free data options everywhere else, GitHub's LFS storage and bandwidth not only cost a lot, but also make it impossible to have public forks and pull request for your repository.  Technically this is a problem only for GitHub's LFS (since it stems from the pricing rules); and can be avoided by using LFS on GitLab or other platform, as [Jim Hester has described](https://github.com/jimhester/test-glfs/).  Still, this proved [unsuccessful for me](https://github.com/jimhester/test-glfs/issues/2), and still faces the other big issue with `git-lfs`:

- Overwrites `git` itself.  Git LFS is just *too* integrated into `git` -- it replaces your authentic `git` engine with `git-lfs`, such that the identical `git` command can have different behaviors on a machine with `git-lfs` installed vs just plain `git`.  Maybe fine for a professional team that is "all in" on `git-lfs`, but is a constant source of pitfalls when working with students and moving between machines that all have only authentic `git` installed.  The difficulties with supporting pull requests etc are also related to this -- in some sense, once you have a `git-lfs` repository, you're really using an entirely new version control system that isn't going to be 100% compatible with the nearly-ubiquitous authentic `git`.

### Amazon S3

Amazon S3 is perhaps the most universal and most obvious go-to place for online-available public and private data storage.  The 5 GB/mo free tier is nice and the pricing is very reasonable and only very incremental after that.  It is easily the most industry-standard solution, and still probably the best way to go in many cases.  It is probably the most scalable solution for very large data, and the only such that has built in support/integration to larger query services like Apache Spark / `sparklyr`.  It falls short of my own use case though in the authentication area. I require students create a GitHub account for my courses and my lab group.  I don't like requiring such third-party accounts, but this one is fundamental to our daily use in classroom and in research, and most of them will continue using the service afterwards.  I particularly don't like having people create complex accounts that they might not even use much in the class or afterwards, just to deal with some pesky minor issue of some data file that is just a little too big for GitHub.

Amazon's authentication is also much more complex than GitHub's passwords or tokens, as is the process of uploading and downloading data from S3 (though the `aws.s3` R package is rather nice remedy here, it doesn't conform to the same user API as the `aws-cli` (python) tool, leaving some odd quirks and patterns that don't match standard Linux commands.)  Together, these make it significantly more difficult to deploy as a quick solution for moving private data around with private repositories.

### Scientific repositories with private storage

For scientific research purposes, this would be my ideal solution.  Encouraging researchers to submit data to a repository at the time of publication is always a challenge, since doing so inevitably involves time & effort and the immediate benefit to the researcher is relatively minimal.  If uploading the data to a repository served an immediate practical purpose of facilitating collaboration, backing up and possibly versioning data, etc, during the research process itself rather than after all is said and done, it would be much more compelling.  Several repositories permit sharing of private data, at least up to some threshold, including DataONE and figshare.  Unfortunately, at this time, I have found the interfaces and R tooling for these too limited or cumbersome for everyday use.

### [datastorr](https://github.com/traitecoevo/datastorr)

The `piggyback` approach is partly inspired by the strategy used in the `datastorr` package, which also uploads data to GitHub releases.  `datastorr` envisions a rather different workflow around this storage strategy, based on the concept of an R "data package" rather than the Git LFS.  I am not a fan of the "data package" approach in general -- I think data should be stored in a platform agnostic way, not as `.Rdata` files, and I often want to first download my data to disk and read it with dedicated functions, not load it "auto-magically" as a package.  This latter issue is particularly important when the data files are larger than what can conveniently fit into working memory, and is better accessed as a database (e.g. SQLite for tabular data, postgis spatial data, etc).

In terms of practical implementation, `datastorr` also creates a new release every time the data file is updated, rather than letting you overwrite files.  In principle `piggyback` will let you version data this way as well, simply create a new release first using `pb_new_release(tag="v2")` or whatever tag you like.  I have not opted for this workflow since in reality, versioning data with releases this way is technically equivalent to creating a new folder for each new version of the data and storing that -- unlike true git commits, release assets such as `datastorr` creates can be easily deleted or overwritten.  I still believe permanent versioned archives like Zenodo should be used for long-term versioned distribution.  Meanwhile, for day-to-day use I often want to overwrite data files with their most recent versions.  (In my case these 'data' files are most often created from upstream data and/or other possibly-long-running code, and are tracked for convenience.  As such they often change as a result of continued work on the upstream processing code.  Perhaps this is not the case for many users and more attention should be paid to versioning.)


### Sharding on GitHub

Another creative solution (hack), at least for some file types, is to break large files into multiple smaller files, and commit those to one or many GitHub repositories.  While [sharding](https://en.wikipedia.org/wiki/Shard_(database_architecture)) is sometimes a legitimate strategy, it has many obvious practical disadvantages and limitations.





% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pb_download_url.R
\name{pb_download_url}
\alias{pb_download_url}
\title{Get the download url of a given file}
\usage{
pb_download_url(
  file = NULL,
  repo = guess_repo(),
  tag = "latest",
  .token = get_token()
)
}
\arguments{
\item{file}{name or vector of names of files to be downloaded. If \code{NULL},
all assets attached to the release will be downloaded.}

\item{repo}{Repository name in format "owner/repo". Will guess the current
repo if not specified.}

\item{tag}{tag for the GitHub release to which this data should be attached.}

\item{.token}{GitHub authentication token, see \verb{[gh::gh_token()]}}
}
\value{
the URL to download a file
}
\description{
Returns the URL download for a public file. This can be useful when writing
scripts that may want to download the file directly without introducing any
dependency on \code{piggyback} or authentication steps.
}
\examples{
\dontrun{

pb_download_url("iris.tsv.xz",
                repo = "cboettig/piggyback-tests",
                tag = "v0.0.1")

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pb_download.R
\name{pb_download}
\alias{pb_download}
\title{Download data from an existing release}
\usage{
pb_download(
  file = NULL,
  dest = ".",
  repo = guess_repo(),
  tag = "latest",
  overwrite = TRUE,
  ignore = "manifest.json",
  use_timestamps = TRUE,
  show_progress = TRUE,
  .token = get_token()
)
}
\arguments{
\item{file}{name or vector of names of files to be downloaded. If \code{NULL},
all assets attached to the release will be downloaded.}

\item{dest}{name of vector of names of where file should be downloaded.
Can be a directory or a list of filenames the same length as \code{file}
vector. Any directories in the path provided must already exist.}

\item{repo}{Repository name in format "owner/repo". Will guess the current
repo if not specified.}

\item{tag}{tag for the GitHub release to which this data should be attached.}

\item{overwrite}{Should any local files of the same name be overwritten?
default \code{TRUE}.}

\item{ignore}{a list of files to ignore (if downloading "all" because
\code{file=NULL}).}

\item{use_timestamps}{DEPRECATED.}

\item{show_progress}{logical, show a progress bar be shown for uploading?
Defaults to \code{TRUE}.}

\item{.token}{GitHub authentication token, see \verb{[gh::gh_token()]}}
}
\description{
Download data from an existing release
}
\examples{
\dontrun{
 ## Download a specific file.
 ## (dest can be omitted when run inside and R project)
 piggyback::pb_download("iris.tsv.gz",
                        repo = "cboettig/piggyback-tests",
                        dest = tempdir())
}
\dontrun{
 ## Download all files
 piggyback::pb_download(repo = "cboettig/piggyback-tests",
                        dest = tempdir())

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pb_new_release.R
\name{pb_new_release}
\alias{pb_new_release}
\title{Create a new release on GitHub repo}
\usage{
pb_new_release(
  repo = guess_repo(),
  tag,
  commit = NULL,
  name = tag,
  body = "Data release",
  draft = FALSE,
  prerelease = FALSE,
  .token = get_token()
)
}
\arguments{
\item{repo}{Repository name in format "owner/repo". Will guess
the current repo if not specified.}

\item{tag}{tag to create for this release}

\item{commit}{Specifies the commit-ish value that
determines where the Git tag is created from.
Can be any branch or commit SHA. Unused if the
git tag already exists. Default: the repository's
default branch (usually \code{master}).}

\item{name}{The name of the release. Defaults to tag.}

\item{body}{Text describing the contents of the tag.
default text is "Data release".}

\item{draft}{default \code{FALSE}. Set to \code{TRUE} to create
a draft (unpublished) release.}

\item{prerelease}{default \code{FALSE}. Set to \code{TRUE} to
identify the release as a pre-release.}

\item{.token}{GitHub authentication token, see \verb{[gh::gh_token()]}}
}
\description{
Create a new release on GitHub repo
}
\examples{
\dontrun{
pb_new_release("cboettig/piggyback-tests", "v0.0.5")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pb_list.R
\name{pb_list}
\alias{pb_list}
\title{List all assets attached to a release}
\usage{
pb_list(
  repo = guess_repo(),
  tag = NULL,
  ignore = "manifest.json",
  .token = get_token()
)
}
\arguments{
\item{repo}{Repository name in format "owner/repo". Will guess the current
repo if not specified.}

\item{tag}{which release tag do we want information for? If \code{NULL} (default),
will return a table for all available release tags.}

\item{ignore}{a list of files to ignore (if downloading "all" because
\code{file=NULL}).}

\item{.token}{GitHub authentication token, see \verb{[gh::gh_token()]}}
}
\value{
a data.frame of release asset names, (normalized to local paths), release tag,
timestamp, owner, and repo.
}
\description{
List all assets attached to a release
}
\details{
To preserve path information, local path delimiters are converted to \verb{.2f}
when files are uploaded as assets.  Listing will display the local filename,
with asset names converting the \verb{.2f} escape code back to the system delimiter.
}
\examples{
\dontrun{
pb_list("cboettig/piggyback-tests")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pb_delete.R
\name{pb_delete}
\alias{pb_delete}
\title{Delete an asset attached to a release}
\usage{
pb_delete(
  file = NULL,
  repo = guess_repo(),
  tag = "latest",
  .token = get_token()
)
}
\arguments{
\item{file}{file(s) to be deleted from the release. If \code{NULL} (default
when argument is omitted), function will delete all attachments to the release.
delete}

\item{repo}{Repository name in format "owner/repo". Will guess the current
repo if not specified.}

\item{tag}{tag for the GitHub release to which this data should be attached.}

\item{.token}{GitHub authentication token, see \verb{[gh::gh_token()]}}
}
\value{
\code{TRUE} (invisibly) if a file is found and deleted.
Otherwise, returns \code{NULL} (invisibly) if no file matching the name was found.
}
\description{
Delete an asset attached to a release
}
\examples{
\dontrun{
readr::write_tsv(mtcars, "mtcars.tsv.gz")
## Upload
pb_upload("mtcars.tsv.gz",
          repo = "cboettig/piggyback-tests",
           overwrite = TRUE)
pb_delete("mtcars.tsv.gz",
          repo = "cboettig/piggyback-tests",
          tag = "v0.0.1")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pb_upload.R
\name{pb_upload}
\alias{pb_upload}
\title{Upload data to an existing release}
\usage{
pb_upload(
  file,
  repo = guess_repo(),
  tag = "latest",
  name = NULL,
  overwrite = "use_timestamps",
  use_timestamps = NULL,
  show_progress = TRUE,
  .token = get_token(),
  dir = "."
)
}
\arguments{
\item{file}{path to file to be uploaded}

\item{repo}{Repository name in format "owner/repo". Will guess the current
repo if not specified.}

\item{tag}{tag for the GitHub release to which this data should be attached.}

\item{name}{name for uploaded file. If not provided will use the basename of
\code{file} (i.e. filename without directory)}

\item{overwrite}{overwrite any existing file with the same name already
attached to the on release? Default behavior is based on timestamps,
only overwriting those files which are older.}

\item{use_timestamps}{DEPRECATED.}

\item{show_progress}{logical, show a progress bar be shown for uploading?
Defaults to \code{TRUE}.}

\item{.token}{GitHub authentication token, see \verb{[gh::gh_token()]}}

\item{dir}{directory relative to which file names should be based.}
}
\description{
NOTE: you must first create a release if one does not already exists.
}
\examples{
\dontrun{
# Needs your real token to run

readr::write_tsv(mtcars,"mtcars.tsv.xz")
pb_upload("mtcars.tsv.xz", "cboettig/piggyback-tests")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/piggyback.R
\docType{package}
\name{piggyback-package}
\alias{piggyback}
\alias{piggyback-package}
\title{piggyback: Managing Larger Data on a GitHub Repository}
\description{
Because larger (> 50 MB) data files cannot easily be committed to git,
a different approach is required to manage data associated with an analysis in a
GitHub repository.  This package provides a simple work-around by allowing larger
(up to 2 GB) data files to piggyback on a repository as assets attached to individual
GitHub releases.  These files are not handled by git in any way, but instead are
uploaded, downloaded, or edited directly by calls through the GitHub API. These
data files can be versioned manually by creating different releases.  This approach
works equally well with public or private repositories.  Data can be uploaded
and downloaded programmatically from scripts. No authentication is required to
download data from public repositories.
}
\details{
It has two main modes or workflows:
\itemize{
\item \code{\link[=pb_upload]{pb_upload()}} / \code{\link[=pb_download]{pb_download()}}:  Upload and download individual files to/from
the desired release of the specified repository
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/piggyback}
  \item Report bugs at \url{https://github.com/ropensci/piggyback/issues}
}

}
\author{
\strong{Maintainer}: Carl Boettiger \email{cboettig@gmail.com} (\href{https://orcid.org/0000-0002-1642-628X}{ORCID}) [copyright holder]

Other contributors:
\itemize{
  \item Mark Padgham (\href{https://orcid.org/0000-0003-2172-5265}{ORCID}) [contributor]
  \item Jeffrey O Hanson (\href{https://orcid.org/0000-0002-4716-6134}{ORCID}) [contributor]
  \item Kevin Kuo (\href{https://orcid.org/0000-0001-7803-7901}{ORCID}) [contributor]
}

}
