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
[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/onekp.svg?branch=master)](https://travis-ci.org/ropensci/onekp)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/onekp/master.svg)](https://codecov.io/github/ropensci/onekp?branch=master)

# onekp

The [1000 Plants initiative
(1KP)](https://sites.google.com/a/ualberta.ca/onekp/) provides the
transcriptome sequences to over 1000 plants from diverse lineages. `onekp`
allows researchers in plant genomics and transcriptomics to access this dataset
through a simple R interface. The metadata for each transcriptome project is
scraped from the 1KP project website. This metadata includes the species,
tissue, and research group for each sequence sample. `onekp` leverages the
taxonomy program `taxizedb`, a local database version of `taxize` package, to
allow filtering of the metadata by taxonomic group (entered as either a taxon
name or NCBI ID). The raw nucleotide or translated peptide sequence can then be
downloaded for the full, or filtered, table of transcriptome projects. 

## Alternatives to `onekp`

The data may also be accessed directly through CyVerse (previously iPlant).
CyVerse efficiently distributes data using the iRODS data system. This approach
is preferable for high-throughput cases or in where iRODS is already in play.
Further, accessing data straight from the source at CyVerse is more stable than
scraping it from project website. However, the `onekp` R package is generally
easier to use (no iRODS dependency or CyVerse API) and offers powerful
filtering solutions. 

## Contact info

1KP staff

 * [Gane Ka-Shu Wong](https://sites.google.com/a/ualberta.ca/professor-gane-ka-shu-wong/) - Principal investigator

 * [Michael Deyholos](mkdeyholos@gmail.com) - Alberta co-investigator

 * [Yong Zhang](zhangy@genomics.org.cn) - Shenzhen co-investigator

 * [Eric Carpenter](ejc@ualberta.ca) - Database manager

R package maintainer

 * [Zebulun Arendsee](arendsee@iastate.edu)


## Installation

`onekp` is on CRAN, but currently is a little out of date. So for now it is
better to install through github. 


```r
library(devtools)
install_github('ropensci/onekp')
```

## Examples

Retrieve the protein and gene transcript FASTA files for two 1KP transcriptomes: 


```r
onekp <- retrieve_onekp()
seqs <- filter_by_code(onekp, c('URDJ', 'ROAP'))
download_peptides(seqs, 'oneKP/pep')
download_nucleotides(seqs, 'oneKP/nuc')
```

This will create the following directory:

```
oneKP
 ├── nuc 
 │   ├── ROAP.fna
 │   └── URDJ.fna
 └── pep
     ├── ROAP.faa
     └── URDJ.faa
```

`onekp` can also filter by species names, taxon ids, or clade.


```r
# filter by species name
filter_by_species(onekp, 'Pinus radiata')

# filter by species NCBI taxon ID
filter_by_species(onekp, 3347)

# filter by clade name scientific name (get all data for the Brassicaceae family)
filter_by_clade(onekp, 'Brassicaceae')

# filter by clade NCBI taxon ID
filter_by_clade(onekp, 3700)
```

So to get the protein sequences for all species in Brassicaceae:


```r
onekp <- retrieve_onekp()
seqs <- filter_by_clade(onekp, 'Brassicaceae')
download_peptides(seqs, 'oneKP/pep')
download_nucleotides(seqs, 'oneKP/nuc')
```

# Funding

Development of this R package was supported by the National Science Foundation under Grant No. IOS 1546858.

# Contributing

We welcome any contributions!

 * By participating in this project you agree to abide by the terms outlined in
   the [Contributor Code of Conduct](CONDUCT.md).

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# onekp 0.2.1

 First rOpenSci version

 * Transfer ownership to rOpenSci

 * Update README
 
 * Add .github files

 * Add reviewers to DESCRIPTION

# onekp 0.2.0

 The goal of this version is to address all issues raised by the rOpenSci
 editor (Scott Chamberlain) and my two reviewers: Jessica Minnier and Zachary
 Foster.

## Updates based on editor comments

 * Added `URL` tag to DESCRIPTION

 * Added `BugReports` tag to DESCRIPTION

 * Wrapped all lines longer than 80 characters

 * Removed unused imports `readr` and `tidyr`

 * Split test functions into smaller test files

 * Add `skip_on_cran` and `skip_on_travis` to `retrieve_onekp` test

 * Save metadata as package data: `data(onekp)` 

 * Add minimum version to `taxizedb`

## Updates based on Zachary Foster's review

 * Add space between equal signs in argument lists

 * Change package name from `oneKP` to `onekp`

 * Let user set the download directory (defaulting to a tempdir).

 * Re-export the `magrittr` `%>%` operator

 * Update README

 * Add CONDUCT.md

 * Add R `devel` test to travis

 * More package-level documentation

 * Do not export print.OneKP function

 * Add documentation details to `retrieve_onekp`

 * Update intro vignette

 * Do not run examples that access the internet (call `retrieve_onekp`)

## Updates based on Jessica Minnier's review  

 There was some overlap between the two reviews, here I just list updates made that are unique to this review.

 * Clarify `filter_by_code` and `filter_by_species` examples

 * Cleanup the `download` documentation

 * Add `absolute` parameter to `download_*` functions. If TRUE, this returns
   the absolute path to the retrieved data.

 * Add handling for downloading nothing: if everything is filtered out, just
   return an empty character vector.

 * Add test coverage for corner cases

   - filtering by clade with a clade that is not in the NCBI common tree raises
     an error.

   - filtering with valid input, but where no rows are selected, returns
     a OneKP object with 0 rows (no error).

   - download an empty object returns an empty vector of filenames (no error)

# onekp 0.1.0

 * Initial release
## Test environments
* local linux-86_64-arch, R 3.4.3
* ubuntu 12.04 (on travis-ci), R 3.4.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
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

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/onekp/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/onekp.git`
* Make sure to track progress upstream (i.e., on our version of `onekp` at `ropensci/onekp`) by doing `git remote add upstream https://github.com/ropensci/onekp.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/onekp`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
