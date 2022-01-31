
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- devtools::rmarkdown::render("README.Rmd") -->

<!-- Rscript -e "library(knitr); knit('README.Rmd')" -->

# Install and run programs, outside of R, inside of R <img src="logo.png" height="200" align="right"/>

[![Build
Status](https://travis-ci.org/ropensci/outsider.svg?branch=master)](https://travis-ci.org/ropensci/outsider)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider)
[![Coverage
Status](https://coveralls.io/repos/github/ropensci/outsider/badge.svg?branch=master)](https://coveralls.io/github/ropensci/outsider?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3615177.svg)](https://doi.org/10.5281/zenodo.3615177)
[![ropensci](https://badges.ropensci.org/282_status.svg)](https://github.com/ropensci/software-review/issues/282)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02038/status.svg)](https://doi.org/10.21105/joss.02038)
[![CRAN
downloads](http://cranlogs.r-pkg.org/badges/grand-total/outsider)](https://CRAN.R-project.org/package=outsider)

> The Outsider is always unhappy, but he is an agent that ensures the
> happiness for millions of ‘Insiders’.<br><br> *[The Outsider,
> Wilson, 1956](https://en.wikipedia.org/wiki/The_Outsider_\(Colin_Wilson\)).*

<br> Integrating external programs into a deployable, R workflow can be
challenging. Although there are many useful functions and packages (e.g.
`base::system()`) for calling code and software from alternative
languages, these approaches require users to independently install
dependant software and may not work across platforms. `outsider` aims to
make this easier by allowing users to install, run and control programs
*outside of R* across all operating systems.

It’s like [whalebrew](https://github.com/whalebrew/whalebrew) but
exclusively for R.

**For more detailed information, check out the [`outsider`
website](https://docs.ropensci.org/outsider/articles/outsider.html)**

## Installation

To install the development version of the package …

``` r
remotes::install_github('ropensci/outsider')
```

Additionally, you will also need to install **Docker desktop**. To
install Docker visit the Docker website and follow the instructions for
your operating system: [Install
Docker](https://www.docker.com/products/docker-desktop).

### Compatibility

Tested and functioning on Linux, Mac OS and Windows. (For some older
versions of Windows, the legacy [Docker
Toolbox](https://docs.docker.com/toolbox/toolbox_install_windows/) may
be required instead of Docker Desktop.)

## Quick example

``` r
library(outsider)
#> ----------------
#> outsider v 0.1.1
#> ----------------
#> - Security notice: be sure of which modules you install
# outsider modules are hosted on GitHub and other code-sharing sites
# this repo is a demonstration outsider module
# it contains a function for printing 'Hello World!' in Ubuntu 18.04
repo <- 'dombennett/om..hello.world'
module_install(repo = repo, force = TRUE)

# look up the help files for the module
module_help(repo = repo)

# import the 'hello_world' function
hello_world <- module_import(fname = 'hello_world', repo = repo)

# run the imported function
hello_world()
#> Hello world!
#> ------------
#> DISTRIB_ID=Ubuntu
#> DISTRIB_RELEASE=18.04
#> DISTRIB_CODENAME=bionic
#> DISTRIB_DESCRIPTION="Ubuntu 18.04.1 LTS"
```

## Available external programs

Modules available on GitHub since 12:08 30 May 2020 (CEST)

● astral

● beast

● PyRate

● RAxML

● pasta

…. Plus, at least, 10 more\!

For more details, see the [available modules
table](https://docs.ropensci.org/outsider/articles/available.html)

### Real-World Example: Aligning biological sequences

Installing and running a multiple sequence alignment program
([mafft](https://mafft.cbrc.jp/alignment/software/)).

![](https://raw.githubusercontent.com/ropensci/outsider/master/other/alignment_example.gif)

(See [“Evolutionary tree
pipeline”](https://docs.ropensci.org/outsider/articles/phylogenetic_pipeline.html)
for running this program yourself.)

### Not finding a module you need?

Try raising an issue to request someone make a module, [Raise an
Issue](https://github.com/ropensci/outsider/issues/new).

Otherwise, why not make it yourself? Check out the
[`outsider.devtools`](https://github.com/ropensci/outsider.devtools)
package.

## Security notice :rotating\_light:

There is a risk that `outsider` modules may be malicious. Modules make
use of the program Docker which allows any program to be efficiently
deployed by wrapping the program’s code, along with everything that
program requires to run e.g. operating system, dependent libraries, into
an executable container.

While this is useful for providing users with whichever programs they
require, there is a potential security risk if, along with the desired
program and dependencies, malicious software is also shipped.

A well-known malicious example of Docker container exploitation is in
cryptocurrency mining. A container may ship with a cryptocurrency mining
software that would make use of your computer’s resources while you ran
you the module.

To minimise any security risks **Be sure of which modules you install on
your machine.** Whenever installing a new module, `outsider` will alert
you to the potential security risks. Before installing a new module, ask
yourself:

  - Is this module from a well-known developer?
  - How many others are using this module?

Consider checking the stats on the module’s GitHub page (e.g. number of
stars/watchers) or looking-up the details of the developer (e.g. email
forums, twitter, academic profile).

Additionally, you may wish to check the Dockerfile of the module. Does
it install programs from well-known code repositories (e.g. `apt-get`)?
Or is it running lines of code from unknown/untrackable URL sources?

## How does it work?

`outsider` makes use of the program [docker](https://www.docker.com/)
which allows users to create small, deployable machines, called Docker
images. The advantage of these images is that they can be run on any
machine that has Docker installed, regardless of operating system. The
`outsider` package makes external programs available in R by
facilitating the interaction between Docker and the R console through
**outsider modules**. These modules consist of two parts: a Dockerfile
that describes the Docker image that contains the external program and
an R package for interacting with the Docker image. Upon installing and
running a module through `outsider`, a Docker image is launched and the
R code of the module is used to interact with the external program.
Anyone can create a module. They are hosted on
[GitHub](https://github.com/) as well as other code-sharing sites and
can be searched for and downloaded through
`outsider`.

![outsider\_outline](https://raw.githubusercontent.com/ropensci/outsider/master/other/outline.png)

## Outsider CI statuses

*Statuses of package building checks and tests, run
monthly.*

| Repo                                                                      | Linux ([Travis CI](https://travis-ci.org/))                                                                                                 | Windows 10 ([Appveyor](https://www.appveyor.com/))                                                                                                                                                 |
| ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`outsider.base`](https://github.com/ropensci/outsider.base)              | [![Build Status](https://travis-ci.org/ropensci/outsider.base.svg?branch=master)](https://travis-ci.org/ropensci/outsider.base)             | [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider.base?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider.base)             |
| [`outsider`](https://github.com/ropensci/outsider)                        | [![Build Status](https://travis-ci.org/ropensci/outsider.svg?branch=master)](https://travis-ci.org/ropensci/outsider)                       | [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider)                       |
| [`outsider.devtools`](https://github.com/ropensci/outsider.devtools)      | [![Build Status](https://travis-ci.org/ropensci/outsider.devtools.svg?branch=master)](https://travis-ci.org/ropensci/outsider.devtools)     | [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider.devtools?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider.devtools)     |
| [Outsider Test suites](https://github.com/ropensci/outsider-testsuites)\* | [![Build Status](https://travis-ci.org/ropensci/outsider-testsuites.svg?branch=master)](https://travis-ci.org/ropensci/outsider-testsuites) | [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider-testsuites?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider-testsuites) |

*\*Mock pipelines to test the interaction of all the packages.*

## Version

Released version 0.1, see
[NEWS](https://github.com/ropensci/outsider/blob/master/NEWS.md).

## Citation

Bennett et al., (2020). outsider: Install and run programs, outside of
R, inside of R. Journal of Open Source Software, 5(45), 2038,
<https://doi.org/10.21105/joss.02038>

## Maintainer

[Dom
Bennett](https://github.com/DomBennett)

-----

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# outsider 0.1.1

## ROpenSci accepted version

* Improved documentation
* Test suites
* Better security warnings/information

# outsider 0.1.0

## Post-review version of outsider

* Split into three packages: `outsider.base`, `outsider` and `outsider.devtools`
* More `module_` functions
* Modules hostable on more code-sharing sites and available offline
* More security notices
* SSH commands to remote host

# outsider 0.0.0

## Pre-review version of outsider

* Install outsider modules
* Search outsider modules on GitHub
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
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
# Contributing

You are very welcome to help out in the development of `outsider`!

If you have any ideas for future features then please create an issue.
Otherwise, if you have the guile, time and inspriation dig deep into the code
and fix/create, then please fork and send a pull request!

## Outsider structure

The `outsider` universe is split between four different entities:

* [`outsider.base`](https://github.com/ropensci/outsider.base): low-level
interface between Docker and modules.
* [`outsider`](https://github.com/ropensci/outsider): the user-friendly
functions for running and installing modules.
* [`outsider.devtools`](https://github.com/ropensci/outsider.devtools):
toolkit for making it easier to develop modules.
* ["Test suites"](https://github.com/ropensci/outsider-testsuites):
pipelines connecting modules to test holistic functionality.

If you wish to raise an issue or create a pull-request, please first
determine which part of the `outsider` structure is most relevant. For example,
if the problem seems to be Docker-related it is likely `outsider.base`. If it
seems to be module discoverability it is likely `outsider`. (If in doubt
select `outsider` as default.)

## Areas for possible contribution

### More modules!

If you think that you can code up a module that does not yet exist in `outsider`
module form -- then please it give a go!

A whole other [package](https://github.com/ropensci/outsider.devtools) and
[website](https://docs.ropensci.org/outsider.devtools/) exists for helping
uses contribute modules.

### Documentation

Any help with documentation to do with the development of modules would be
greatly appreciated. In particular, more examples of module designs/structures
along with descriptions. Additionally, any more tips and tricks or contributions
to the FAQs would be equally welcomed.

### Scalability

The idea of `outsider` was to allow incorporation of code outside of R into R.
But some of this integration can come at a computational cost: code that may
have otherwise been run outside of an R analysis on a powerful desktop or
a HPC facility could be run on a desktop instead. One solution so far has been
to enable SSH'ing to remote servers. Any ideas on improving this or
documentation on setting up hosting services (AWS, Google Cloud, Azure ... etc.)
that can interact with `outsider` would be greatly welcomed.

### New test suite

`outsider` aims to make it easier to string pipelines together in R that make
use of multiple external programs. To ensure that this is continually
functional, as the project and its packages develop, a GitHub repo
called "outsider-testsuites" is used to run a series of pipelines ("suites")
that run domain-specific pipelines. For example, "suite_1" runs a biological
analysis for inferring an evolutionary tree of pineapple-like plants.

If you have created a series of `outsider` modules for a specific domain and
think they would make for a given test suite, then please add a new pipeline to
the
["outsider-testsuites"](https://github.com/ropensci/outsider-testsuites)!

## How to contribute to an R package

To contribute you will need a GitHub account and to have knowledge of
R (plus, depending on where you contribute, knowledge of bash and Docker).
You can then create a fork of the repo in your own GitHub account
and download the repository to your local machine. `devtools` is recommended.

```r
devtools::install_github('[your account]/outsider')
```

All new functions must be tested. For every new file in `R/`, a new test file
must be created in `tests/testthat/`. To test the package and make sure it
meets CRAN guidelines use `devtools`. 

```r
devtools::test()
devtools::check_cran()
```

For help, refer to Hadley Wickham's book, [R packages](http://r-pkgs.had.co.nz/).

## Style guide

`outsider` is developed for submission to ROpenSci. This means the package and
its code should meet ROpenSci style and standards. For example, function
names should be all lowercase, separated by underscores and the last word
should, ideally, be a verb.

```
# e.g.
species_ids_retrieve()  # good
sppIDs()                # not great
sp.IDS_part2()          # really bad
sigNXTprt.p()           # awful
```

It is best to make functions small, with specific names. Feel free to break up code into multiple separate files (e.g. tools,
helper functions, stages ...). For more details and better explanations refer to the ROpenSci [guide](https://devguide.ropensci.org/building.html).

## API tokens

The repository search features can make use of private tokens which make
searching faster and more reliable. In the case of GitLab, however, an
"access token" is a requirement for using the API.

A token is set up by generating one by logging into your code-sharing service
and then saving the generated token to an R environment file.

To create a token for your code-sharing site, visit either
https://github.com/settings/tokens for GitHub or
https://gitlab.com/profile/personal_access_tokens for GitLab. The generated
tokens should allow searching of all public repositories. Then add the token to
your R environment using `usethis::edit_r_environ()`. The `.Renviron` file
should contain the lines:

```
GITHUB_PAT=[your private token]
GITLAB_PAT=[your private token]
```

## Code of Conduct

Please note that the 'outsider' project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
---
title: "outsider: Install and run programs, outside of R, inside of R"
tags:
  - R
  - Docker
  - GitHub
  - GitLab
  - BitBucket
  - Reproducibility
  - Workflow
  - Integration
authors:
 - name: Dominic J. Bennett
   orcid: 0000-0003-2722-1359
   affiliation: "1, 2"
 - name: Hannes Hettling
   orcid: 0000-0003-4144-2238
   affiliation: "3"
 - name: Daniele Silvestro
   orcid: 0000-0003-0100-0961
   affiliation: "1, 2"
 - name: Rutger Vos
   orcid: 0000-0001-9254-7318
   affiliation: "3"
 - name: Alexandre Antonelli
   orcid: 0000-0003-1842-9297
   affiliation: "1, 2, 4"
affiliations:
 - name: Gothenburg Global Biodiversity Centre, Box 461, SE-405 30 Gothenburg, Sweden
   index: 1
 - name: Department of Biological and Environmental Sciences, University of Gothenburg, Box 461, SE-405 30 Gothenburg, Sweden
   index: 2
 - name: Naturalis Biodiversity Center, P.O. Box 9517, 2300 RA Leiden, The Netherlands
   index: 3
 - name: Royal Botanic Gardens, Kew, TW9 3AE, Richmond, Surrey, UK
   index: 4
date: 29 January 2020
bibliography: paper.bib

---

# Statement of need

Enable integration of R and non-R code and programs to facilitate reproducible workflows.

# Summary

In many areas of research, product development and software engineering, analytical pipelines – workflows connecting output from multiple software – are key for processing and running tests on data.
They can provide results in a consistent, modular and transparent manner.
Pipelines also make it easier to demonstrate the reproducibility of one’s research as well as enabling analyses that update as new data are added.
Not all analyses, however, can necessarily be run or coded in one’s favoured programming language as different parts of an analysis may require external software or packages.
Integrating a variety of programs and software can lead to issues of portability (additional software may not run across all operating systems) and versioning errors (differing arguments across additional software versions).
For the ideal pipeline, it should be possible to install and run any command-line software, within the main programming language of the pipeline, without concern for software versions or operating system.
R [@cran] is one of the most popular computer languages amongst researchers, and many packages exist for calling programs and code from non-R sources (e.g. `sys` [@sys] for shell commands, `reticulate` [@reticulate] for `python` and `rJava` [@rJava] for `Java`).
To our knowledge, however, no R package exists with the ability to launch external programs originating from *any* UNIX command-line source.

The `outsider` packages work through `docker` [@docker] – a service that, through OS-level virtualization, enables deployment of isolated software "containers" – and a code-sharing service, e.g. GitHub [@github], to allow a user to install and run, in theory, any external, command-line program or package, on any of the major operating systems (Windows, Linux, OSX).

## How it works

`outsider` packages provide an interface to install and run *outsider modules*.
These modules are hostable on GitHub [@github], GitLab [@gitlab] and/or BitBucket [@bitbucket] and consist of two parts: a (barebones) R package and a Dockerfile.
The Dockerfile details the installation process for an external program contained within a Docker image, while the R package comprises functions and documentation for interacting with the external program via a Docker container.
For many programs, Dockerfiles are readily available online and require minor changes to adapt for `outsider`.
By default, a module’s R code simply passes command-line arguments through Docker.
After installation, a module’s functions can then be imported and launched using `outsider` functions.
Upon running a module’s code, `outsider` code will first launch a Docker container of the image as described by the module’s Dockerfile.
`outsider` then facilitates the communication between the module’s R code and the Docker container that hosts the external program (developers of modules have the choice of determining default behaviours for handling generated files).
`outsider` modules thus wrap external command-line programs into R functions in a convenient manner.
`outsider` functions allow users to look up available modules and determine build statuses (i.e. whether the package is passing its online tests) before installing.


At time of writing, `outsider` modules for some of the most popular bioinformatics tools have been developed: BLAST [@blast], MAFFT [@mafft], *BEAST [@beast], RAxML [@raxml], bamm [@bamm], PyRate [@pyrate].
(See the `outsider` website for an up-to-date and [complete list](https://docs.ropensci.org/outsider/articles/available.html)).
All that is required to run these modules is R and Docker.
Docker Desktop [@docker_desktop] can be installed for all operating systems but for older versions of OSX and Windows the legacy "Docker Toolbox" [@docker_toolbox] may instead need to be installed.
(Note, users may need to create an account with Docker-Hub to install Docker.)

![An outline of the outsider module ecosystem.](https://raw.githubusercontent.com/ropensci/outsider/master/other/outline.png)


### Code structure

The code-base that allows for the installation, execution and development of `outsider` modules is held across three different R packages. For end-users of modules, however, only the `outsider` module is required. For those who wish to develop their own modules, the `outsider.devtools` package provides helper functions for doing so. In addition, there is a test suites repository that hosts mock analysis pipelines that initiate several modules in sequence to test the interaction of all the packages.

* [`outsider`](https://github.com/ropensci/outsider): The main package for installing, importing and running **outsider modules** [@outsider_z].
* [`outsider.base`](https://github.com/ropensci/outsider.base): The package for low-level interaction between module R code and Docker containers (not user-facing) [@outsiderbase_z].
* [`outsider.devtools`](https://github.com/ropensci/outsider.devtools): The development tools package for facilitating the creation of new modules [@outsiderdev_z].
* ["outsider-testuites"](https://github.com/ropensci/outsider-testsuites): A repository hosting a series of "test" pipelines for ensuring modules can be successfully strung together to form R/non-R workflows [@outsiderts_z].

![How the outsider packages interact](https://raw.githubusercontent.com/ropensci/outsider.devtools/master/other/package_structures.png)

------

# Examples

## Saying hello from Ubuntu

By hosting a Docker container, `outsider` can run any UNIX-based, external command-line program.
To demonstrate this process we can say "hello world" via a container hosting the [Ubuntu operating system](https://en.wikipedia.org/wiki/Ubuntu).
In this short example, we will install a small `outsider module` – `om..hello.world` – that installs a local copy of the latest version of Ubuntu and contains a function for saying hello using the command `echo`.

```r
library(outsider)
# outsider modules are hosted on GitHub
# this repo is a demonstration of an outsider module
# it contains a function for printing 'Hello World!'
repo <- 'dombennett/om..hello.world'
module_install(repo = repo)

# look up the help files for the module
module_help(repo = repo)

# import the 'hello_world' function
hello_world <- module_import(fname = 'hello_world', repo = repo)

# run the imported function
hello_world()
#> Hello world!
#> ------------
#> DISTRIB_ID=Ubuntu
#> DISTRIB_RELEASE=18.04
#> DISTRIB_CODENAME=bionic
#> DISTRIB_DESCRIPTION="Ubuntu 18.04.1 LTS"
```

## A basic bioinformatic pipeline

To better demonstrate the power of the `outsider` package, we will run a simple bioinformatic pipeline that downloads a file of biological sequence data (borrowed from [@genomeworkshop]) and aligns the separate strands of DNA using the [multiple sequence alignment](https://en.wikipedia.org/wiki/Multiple_sequence_alignment) program MAFFT [@mafft].
Note that we can pass arguments to an `outsider` module, such as `mafft` in the example below, using separate R arguments for each command-line argument.

```r
library(outsider)
repo <- 'dombennett/om..mafft'
module_install(repo = repo)
mafft <- module_import(fname = 'mafft', repo = repo)

# some example file
download.file('https://molb7621.github.io/workshop/_downloads/sample.fa',
              'sample.fa')

# run maft with --auto and write results to alignment.fa
mafft(arglist = c('--auto', 'sample.fa', '>', 'alignment.fa'))

# view alignment
cat(readLines('alignment.fa'), sep = '\n')
#> >derice
#> -actgactagctagctaactg
#> >sanka
#> -gcatcgtagctagctacgat
#> >junior
#> catcgatcgtacgtacg-tag
#> >yul
#> -atcgatcgatcgtacgatcg
```

**For more detailed and up-to-date examples and tutorials, see the `outsider` GitHub page [@outsider_gh].**

# Availability

`outsider` (and its sister packages) are open-source software made available under the MIT licence allowing reuse of the software with limited constraints.
It is aimed that all packages will be made available through CRAN [@cran], e.g. `install.package("outsider")`.
Currently, all are available from GitHub source code repositories using the `remotes` package, e.g. `remotes::install_github("ropensci/outsider")`

# Funding

This package has been developed as part of the supersmartR project [@supersmartR] which has received funding through A.A. (from the Swedish Research Council [B0569601], the Swedish Foundation for Strategic Research and a Wallenberg Academy Fellowship) and through D.S. (from the Swedish Research Council [2015-04748]).

# Acknowledgements

The authors wish to thank [Daniel Nüst](https://github.com/nuest) and [Bob Rudis](https://github.com/hrbrmstr) for taking the time to review the package and providing valuable ideas and suggestions.
Additionally, we would like to thank [Anna Krystalli](https://github.com/annakrystalli) for overseeing the review process and providing useful advice in turn.
Finally, we would like to thank all the people behind [ROpenSci](https://ropensci.org/) and [JOSS](https://joss.theoj.org/) for making this project and its review possible.

# References

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

<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>



---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- devtools::rmarkdown::render("README.Rmd") -->
<!-- Rscript -e "library(knitr); knit('README.Rmd')" -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# Install and run programs, outside of R, inside of R <img src="logo.png" height="200" align="right"/>

[![Build Status](https://travis-ci.org/ropensci/outsider.svg?branch=master)](https://travis-ci.org/ropensci/outsider) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider) [![Coverage Status](https://coveralls.io/repos/github/ropensci/outsider/badge.svg?branch=master)](https://coveralls.io/github/ropensci/outsider?branch=master) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3615177.svg)](https://doi.org/10.5281/zenodo.3615177) [![ropensci](https://badges.ropensci.org/282_status.svg)](https://github.com/ropensci/software-review/issues/282) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02038/status.svg)](https://doi.org/10.21105/joss.02038) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/outsider)](https://CRAN.R-project.org/package=outsider)

> The Outsider is always unhappy, but he is an agent that ensures the happiness for millions of 'Insiders'.<br><br>
*[The Outsider, Wilson, 1956](https://en.wikipedia.org/wiki/The_Outsider_(Colin_Wilson)).*

<br>
Integrating external programs into a deployable, R workflow can be challenging. Although there are many useful functions and packages (e.g. `base::system()`) for calling code and software from alternative languages, these approaches require users to independently install dependant software and may not work across platforms. `outsider` aims to make this easier by allowing users to install, run and control programs *outside of R* across all operating systems.

It's like [whalebrew](https://github.com/whalebrew/whalebrew) but exclusively
for R.

**For more detailed information, check out the [`outsider` website](https://docs.ropensci.org/outsider/articles/outsider.html)**

## Installation

To install the development version of the package ...

```{r outsider-install, eval=FALSE, include=TRUE}
remotes::install_github('ropensci/outsider')
```

Additionally, you will also need to install **Docker desktop**. To install
Docker visit the Docker website and follow the instructions for your operating
system: [Install Docker](https://www.docker.com/products/docker-desktop).

### Compatibility

Tested and functioning on Linux, Mac OS and Windows. (For some older versions of
Windows, the legacy
[Docker Toolbox](https://docs.docker.com/toolbox/toolbox_install_windows/) may
be required instead of Docker Desktop.)

## Quick example

```{r example-setup, eval=TRUE, include=FALSE}
outsider::module_uninstall(repo = 'dombennett/om..hello.world')
```

```{r outsider-example}
library(outsider)
# outsider modules are hosted on GitHub and other code-sharing sites
# this repo is a demonstration outsider module
# it contains a function for printing 'Hello World!' in Ubuntu 18.04
repo <- 'dombennett/om..hello.world'
module_install(repo = repo, force = TRUE)

# look up the help files for the module
module_help(repo = repo)

# import the 'hello_world' function
hello_world <- module_import(fname = 'hello_world', repo = repo)

# run the imported function
hello_world()
```

## Available external programs

```{r, available, echo=FALSE, results='asis'}
knitr::opts_chunk$set(
  comment = ""
)
nmax <- 5
avlbl <- suppressWarnings(outsider::module_details(service = 'github'))
prgms <- unique(avlbl[['program']])
prgms <- sort(prgms[prgms != ''])
prgms <- prgms[prgms != 'hello world']
n <- length(prgms)
if (n >= nmax) {
  prgms <- sample(x = prgms, size = nmax)
} else {
  prgms <- sample(x = prgms)
}
time_date <- as.character(format(Sys.time(), '%H:%M %d %B %Y (%Z)'))
cat('Modules available on GitHub since ', crayon::bold(time_date), '\n\n')
for (prgm in prgms) {
  cli::cat_bullet(prgm, '\n')
}
if (n > nmax) {
  cat('.... Plus, at least, ', crayon::bold(n - nmax), ' more!\n\n')
}
```

For more details, see the [available modules table](https://docs.ropensci.org/outsider/articles/available.html)

### Real-World Example: Aligning biological sequences

Installing and running a multiple sequence alignment program 
([mafft](https://mafft.cbrc.jp/alignment/software/)).

![](https://raw.githubusercontent.com/ropensci/outsider/master/other/alignment_example.gif)

(See ["Evolutionary tree pipeline"](https://docs.ropensci.org/outsider/articles/phylogenetic_pipeline.html) for running this program yourself.)

### Not finding a module you need?

Try raising an issue to request someone make a module,
[Raise an Issue](https://github.com/ropensci/outsider/issues/new).

Otherwise, why not make it yourself? Check out the [`outsider.devtools`](https://github.com/ropensci/outsider.devtools) package.

## Security notice :rotating_light:

There is a risk that `outsider` modules may be malicious. Modules make use of
the program Docker which allows any program to be efficiently deployed by
wrapping the program's code, along with everything that program requires to run
e.g. operating system, dependent libraries, into an executable container.

While this is useful for providing users with whichever programs they require,
there is a potential security risk if, along with the desired program and
dependencies, malicious software is also shipped.

A well-known malicious example of Docker container exploitation is in
cryptocurrency mining. A container may ship with a cryptocurrency mining
software that would make use of your computer's resources while you ran you the
module.

To minimise any security risks **Be sure of which modules you install on your
machine.** Whenever installing a new module, `outsider` will alert you to the
potential security risks. Before installing a new module, ask yourself:

* Is this module from a well-known developer?
* How many others are using this module?

Consider checking the stats on the module's GitHub page (e.g. number of
stars/watchers) or looking-up the details of the developer (e.g. email forums,
twitter, academic profile).

Additionally, you may wish to check the Dockerfile of the module. Does it
install programs from well-known code repositories (e.g. `apt-get`)? Or is it
running lines of code from unknown/untrackable URL sources?

## How does it work?

`outsider` makes use of the program [docker](https://www.docker.com/) which
allows users to create small, deployable machines, called Docker images.
The advantage of these images is that they can be run on any machine that has
Docker installed, regardless of operating system. The `outsider` package makes
external programs available in R by facilitating the interaction between Docker
and the R console through **outsider modules**. These modules consist of two
parts: a Dockerfile that describes the Docker image that contains the external
program and an R package for interacting with the Docker image. Upon installing
and running a module through `outsider`, a Docker image is launched and the R
code of the module is used to interact with the external program. Anyone can
create a module. They are hosted on [GitHub](https://github.com/) as well as
other code-sharing sites and can be searched for and downloaded through
`outsider`.

![outsider_outline](https://raw.githubusercontent.com/ropensci/outsider/master/other/outline.png)

## Outsider CI statuses

*Statuses of package building checks and tests, run monthly.*

|Repo|Linux ([Travis CI](https://travis-ci.org/))|Windows 10 ([Appveyor](https://www.appveyor.com/))|
|---|---|---|
|[`outsider.base`](https://github.com/ropensci/outsider.base)|[![Build Status](https://travis-ci.org/ropensci/outsider.base.svg?branch=master)](https://travis-ci.org/ropensci/outsider.base)|[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider.base?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider.base)|
|[`outsider`](https://github.com/ropensci/outsider)|[![Build Status](https://travis-ci.org/ropensci/outsider.svg?branch=master)](https://travis-ci.org/ropensci/outsider)|[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider)|
|[`outsider.devtools`](https://github.com/ropensci/outsider.devtools)|[![Build Status](https://travis-ci.org/ropensci/outsider.devtools.svg?branch=master)](https://travis-ci.org/ropensci/outsider.devtools)|[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider.devtools?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider.devtools)|
|[Outsider Test suites](https://github.com/ropensci/outsider-testsuites)*|[![Build Status](https://travis-ci.org/ropensci/outsider-testsuites.svg?branch=master)](https://travis-ci.org/ropensci/outsider-testsuites)|[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider-testsuites?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider-testsuites)|

*&ast;Mock pipelines to test the interaction of all the packages.*

## Version

Released version 0.1, see [NEWS](https://github.com/ropensci/outsider/blob/master/NEWS.md).

## Citation

Bennett et al., (2020). outsider: Install and run programs, outside of R, inside of R. Journal of Open Source Software, 5(45), 2038, https://doi.org/10.21105/joss.02038

## Maintainer

[Dom Bennett](https://github.com/DomBennett)

---

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Available modules"
output: rmarkdown::html_vignette
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, comment = '')
```

### Modules on GitHub

```{r available, echo=FALSE, results='asis'}
module_info <- suppressWarnings(outsider::module_details(service = 'github'))
module_info <- module_info[module_info[['program']] != '', ]
# convert to links
module_info[['program']] <- paste0('[', module_info[['program']], '](',
                                   module_info[['url']], ')')
module_info[['url']] <- NULL
module_info[['repo']] <- NULL
module_info[['updated_at']] <- format(module_info[['updated_at']], '%d %b %y')
colnames(module_info) <- c('Program', 'Details', 'Versions', 'Updated',
                           'Stars')
time_date <- as.character(format(Sys.time(), '%H:%M, %d %B %Y (%Z)'))
cat('Last updated: ', crayon::bold(time_date), '\n')
knitr::kable(as.data.frame(module_info), format = "html")
```---
title: "Pipeline to construct evolutionary trees"
date: "2020-01-20"
output: rmarkdown::html_vignette
---



The whole point of `outsider` is to bring lots of programs together into a
single place to allow the development of analysis pipelines. Let's demonstrate
this process by creating a simple pipeline for constructing an evolutionary
tree. We will generate a tree in three steps: obtain orthologous DNA sequences,
align these sequences and then estimate a tree.

In our example, we will keep it simple and fast and generate a tree of
[CytB](https://en.wikipedia.org/wiki/Cytochrome_b) sequences in the primate
genus of Night Monkeys, **Aotus**.

# phylotaR

To get us started we need DNA sequences that represent the same gene region,
i.e. sequences that are orthologous. We can obtain our sequences quickly using
the example data from the R package
[`phylotaR`](https://github.com/ropensci/phylotaR).

```r
# phylotaR is currently not avaialble on CRAN but can be install from Github
remotes::install_github("ropensci/phylotaR")
```

`phylotaR` provides a pipeline for identifying orthologous sequences from a
reputable online biological sequence database. The package comes with a few
pre-calculated example sequences. We can extract sequence clusters that
represent the gene region [cytb](https://en.wikipedia.org/wiki/Cytochrome_b)
with the following script.


```r
aotus <- NULL
```

```r
library(phylotaR)
# Example data
data("aotus")
# Generate summary of identified clusters
smmry <- summary(aotus)
# Extract cluster with 'cytb' in feature name
cytb <- smmry$ID[which(grepl('cytb', smmry$Feature))[1]]
cytb <- drop_clstrs(aotus, cytb)
# Reduce cluster to just one sequence per taxon
cytb <- drop_by_rank(cytb, n = 1)
# Get taxonomic IDs for taxa in cluster
txids <- get_txids(cytb, cid = cytb@cids[[1]])
# Convert IDs to species names
sp_nms <- get_tx_slot(cytb, txids, slt_nm = 'scnm')
sp_nms <- sub(' ', '_', sp_nms)
# Write out
write_sqs(cytb, 'aotus_cytb.fasta', sq_nm = sp_nms, sid = names(sp_nms))
# What do the first 50 lines of the file look like?
cat(paste0(readLines('aotus_cytb.fasta', n = 50), collapse = '\n'))
```

```
## >Aotus_nancymaae
## atgacttctccccgcaaaacacacccactaacaaagatcattaacgaatcattcattgatctacccacaccacccaacat
## ttcctcctgatgaaattttggctcactcttaggcatttgcctaattattcaaatcaccaccggcctgttcttagctatac
## actacacaccagatacctcaaccgccttctcctccgtcgcccatatcacccgagacgtcaactatggctgaataattcgc
## tacatacatgccaacggtgcttccatattcttcgtatgcctttttctccacattggtcgaggactttactatggatcctt
## tctttctctgaagacttgaaatatcggtaccatcctactacttacaaccatagccacagcattcataggctatgttcttc
## catgaggccaaatatcattctgaggggctacagtaattacaaatcttttatcagccattccctatatcggatctgacctt
## gtacaatgaatttgaggtggcttctcagtagataaagccactctcacacgattctttacttttcactttatcttaccctt
## cattatcgcagccctagcaactatccatctattatttctgcatgaaacaggatcaagtaacccatcaggaataacatctg
## accccgacaaaatcacatttcacccctattatacagctaaagacattctaggattaatctttcttctcttatccctaata
## agcctaaccctatttatacccgaccttttaaccgacccagataattatacactggctaatcccctcaacactccacccca
## catcaagccagagtgatattttctatttgcatacgcaatcctacgatctatccctaataaacttggaggagttctagccc
## tagtactttctattttaattctaatagttatccctatactacatctctccaaacaacaaagcataatatttcgacccatc
## actcaaattctattctgaactctagtagctgacctactaactctcacatgaattggaggccaaccggttgaatacccctt
## cgtaaccattggccaaaccgcatccattacatacttcttcattattattatcctaatgcccctttccgcctcaatcgaaa
## atatattacttaaatgataa
## 
## >Aotus_azarai
## atgacctccccccgcaaaacacacccactagcaaagattattaacgaatcattcatcgatctccccacaccatccaacat
## ttcctcttgatgaaattttggctcactcttaggcatttgcctaatcattcaaatcaccaccggcctgttcttagctatac
## attacacaccagatacctcaactgccttctcctccgtcgctcatatcacccgagacgttaactatggctgaataattcgc
## tatatacatgccaacggcgcttccatattcttcgtatgcctttttctccatattggccgaggactttactatggatcttt
## cctttttctgaagacttgaaatatcggtattatcctactacttacaaccatagccacagcattcataggctatgttcttc
## catgaggccaaatatcattctgaggggccacagtaattacaaaccttctatcagctatcccctatatcgggtctgacctt
## gtacaatgaatttgaggtggcttctcagtagataaagccactctcacacgattctttacttttcactttatcttaccctt
## tattatcgcagccctagcaactattcacctcttatttctacatgaaacaggatcaagcaacccatcaggaataacatctg
## accccgacaaagtcacattccacccctattatacagctaaggatattctaggattaatctttcttctcttatccctaata
## agcctaaccctatttatacccgaccttctaacggacccagataattatacactagccaaccctctcaacaccccgcctca
## cattaagccagagtggtattttctatttgcatacgcaattctacgatctatccctaataaacttggaggagtactagccc
## tagtactttctatcctaatcttaatggctatccccgtactacatttctccaaacagcaaagtataatatttcgacccatt
## actcaaattctattctgagctctggtagctgacctactaactctcacgtgaattggaggtcaaccagttgagtacccctt
## cgtaaccattggccagaccgcatccattatatacttcttcattattattaccctaataccccttttcgccttaattgaaa
## ataaattacttaaatgatag
## 
## >Aotus_griseimembra
## atgacttctccccgcaaaacacacccactaacaaagatcattaacgaatcattcattgatctacccacaccacccaacat
## ttcctcctgatgaaattttggctcactcttaggcatttgcctaattattcaaatcaccaccggcctgttcttagctatac
## actacacaccagatacctcaaccgccttctcctccgtcgcccatatcacccgagacgtcaactatggctgaataattcgc
## tacatacatgccaacggtgcttccatattcttcgtatgcctttttctccacattggtcgaggactttactatggatcctt
## cctttctctgaagacttgaaatatcggtaccatcctactacttacaaccatagccacagcattcataggctatgttcttc
## catgaggccaaatatcattctgaggggctacagtaattacaaatcttatatcagccattccctatatcggatctgacctt
## gtacaatgaatttgaggtggcttctcagtagataaagccactctcacacgattctttacttttcactttatcttaccctt
## cattatcgcagccctagcaactatccatctattatttctgcatgaaacaggatcaagtaacccatcaggaataacatctg
## accccgacaaaatcacatttcacccctattatacagctaaagacattctagggttaatctttcttctcttatccctaata
## agcctaaccctatttatacccgaccttttaaccgacccagataattatacactggctaatcccctcaacactccacccca
## catcaagccagagtgatattttctatttgcatacgcaattctacgatctatccctaataaacttggaggagttctagccc
## tagtactttctattttaattctaatagttatccctatactacatctctccaaacaacaaagcataatatttcgacccatc
## actcaaattctattctgaactctagtagctgacctactaactctcacatgaattggaggccaaccggttgaatacccctt
## cgtaaccattggccaaaccgcatccattacatacttcttcattattattatcctaatgcccctttccgcctcaatcgaaa
## atatattacttaaatgataa
```

```r
cat('...')
```

```
## ...
```

# Alignment

Now we have written our sequences to a text-based `.fasta` file, we need to
align them! To do this, we can use
[`mafft`](https://mafft.cbrc.jp/alignment/software/). We first need to install
the module and then import the key function, `mafft`. After doing this, we can
call the `mafft` program using the same arguments that we would if it were
calling the program via command-line. For simplicity, we can run our alignment
using the '--auto' parameter.


```r
library(outsider)
```

```
## ----------------
## outsider v 0.1.0
## ----------------
## - Security notice: be sure of which modules you install
```

```r
# Install module
repo <- 'dombennett/om..mafft'
module_install(repo = repo, force = TRUE)
```

```r
# Import the function
mafft <- module_import(fname = 'mafft', repo = repo)
# Run our program with normal mafft arguments
# Note: all arguments must be separate character variables.
mafft(arglist = c('--auto', 'aotus_cytb.fasta', '>', 'alignment.fasta'))
```

```
## outputhat23=16
## treein = 0
## compacttree = 0
## stacksize: 8192 kb
## generating a scoring matrix for nucleotide (dist=200) ... done
## All-to-all alignment.
##     0 / 8

-
    1 / 8

|
    2 / 8

-
    3 / 8

|
    4 / 8

-
    5 / 8

|
    6 / 8

-
tbfast-pair (nuc) Version 7.407
## alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
## 0 thread(s)
## 
## outputhat23=16
## Loading 'hat3.seed' ... 
## done.
## 
|
Writing hat3 for iterative refinement
## generating a scoring matrix for nucleotide (dist=200) ... done
## Gap Penalty = -1.53, +0.00, +0.00
## tbutree = 1, compacttree = 0
## Constructing a UPGMA tree ... 
## 
    0 / 8
## done.
## 
## 
-
Progressive alignment ... 
## 
STEP     1 /7 
|

STEP     2 /7 
-

STEP     3 /7 
|

STEP     4 /7 
-

STEP     5 /7 
|

STEP     6 /7 
-

STEP     7 /7 
|
## done.
## tbfast (nuc) Version 7.407
## alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
## 1 thread(s)
## 
## 
-
minimumweight = 0.000010
## autosubalignment = 0.000000
## nthread = 0
## randomseed = 0
## blosum 62 / kimura 200
## poffset = 0
## niter = 16
## sueff_global = 0.100000
## nadd = 16
## 
|
Loading 'hat3' ... done.
## generating a scoring matrix for nucleotide (dist=200) ... done
## 
-
## 
    0 / 8
## Segment   1/  1    1-1141
## 
|
STEP 001-001-0 
-
 identical.   
STEP 001-001-1 
|
 identical.   
STEP 001-002-0 
-
 identical.   
STEP 001-002-1 
|
 identical.   
STEP 001-003-0 
-
 identical.   
STEP 001-003-1 
|
 identical.   
STEP 001-004-0 
-
 identical.   
STEP 001-004-1 
|
 identical.   
STEP 001-005-0 
-
 identical.   
STEP 001-005-1 
|
 identical.   
STEP 001-006-0 
-
 identical.   
STEP 001-006-1 
|
 identical.   
STEP 001-007-1 
-
 identical.   
STEP 002-007-1 
|
 identical.   
STEP 002-006-0 
-
 identical.   
STEP 002-006-1 
|
 identical.   
## Converged.
## 
## done
## 
-
dvtditr (nuc) Version 7.407
## alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
## 0 thread(s)
## 
## 
|
## Strategy:
##  L-INS-i (Probably most accurate, very slow)
##  Iterative refinement method (<16) with LOCAL pairwise alignment information
## 
## If unsure which option to use, try 'mafft --auto input > output'.
## For more information, see 'mafft --help', 'mafft --man' and the mafft page.
## 
## The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
## It tends to insert more gaps into gap-rich regions than previous versions.
## To disable this change, add the --leavegappyregion option.
## 
## 
```

# Phylogeny

We can repeat the same process as we did for `mafft` but instead for
[`RAxML`](https://github.com/stamatak/standard-RAxML) -- a command-line program
for estimating evolutionary trees from an alignment of DNA sequences using a
model of molecular evolution (GTRGAMMA) with maximum likelihood.


```r
library(outsider)
# Install
repo <- 'dombennett/om..raxml'
module_install(repo = repo, force = TRUE)
```

```r
# Import
raxml <- module_import(fname = 'raxml', repo = repo)
# Run
raxml(arglist = c('-m', 'GTRGAMMA', '-s', 'alignment.fasta', '-p', '1234',
                  '-n', 'aotus_cytb', '-T', '2'))
```

```
##
## RAxML can't, parse the alignment file as phylip file 
## it will now try to parse it as FASTA file
## 
## 
## 
## Using BFGS method to optimize GTR rate parameters, to disable this specify "--no-bfgs" 
## 
## 
## This is the RAxML Master Pthread
## 
## This is RAxML Worker Pthread Number: 1
## 
## 
## This is RAxML version 8.2.12 released by Alexandros Stamatakis on May 2018.
## 
## With greatly appreciated code contributions by:
## Andre Aberer      (HITS)
## Simon Berger      (HITS)
## Alexey Kozlov     (HITS)
## Kassian Kobert    (HITS)
## David Dao         (KIT and HITS)
## Sarah Lutteropp   (KIT and HITS)
## Nick Pattengale   (Sandia)
## Wayne Pfeiffer    (SDSC)
## Akifumi S. Tanabe (NRIFS)
## Charlie Taylor    (UF)
## 
## 
## Alignment has 74 distinct alignment patterns
## 
## Proportion of gaps and completely undetermined characters in this alignment: 0.00%
## 
## RAxML rapid hill-climbing mode
## 
## Using 1 distinct models/data partitions with joint branch length optimization
## 
## 
## Executing 1 inferences on the original alignment using 1 distinct randomized MP trees
## 
## All free model parameters will be estimated by RAxML
## GAMMA model of rate heterogeneity, ML estimate of alpha-parameter
## 
## GAMMA Model parameters will be estimated up to an accuracy of 0.1000000000 Log Likelihood units
## 
## Partition: 0
## Alignment Patterns: 74
## Name: No Name Provided
## DataType: DNA
## Substitution Matrix: GTR
## 
## 
## 
## 
## RAxML was called as follows:
## 
## raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -s alignment.fasta -p 1234 -n aotus_cytb -T 2 
## 
## 
## Partition: 0 with name: No Name Provided
## Base frequencies: 0.289 0.285 0.120 0.306 
## 
## Inference[0]: Time 0.056735 GAMMA-based likelihood -2591.962827, best rearrangement setting 5
## 
## 
## Conducting final model optimizations on all 1 trees under GAMMA-based models ....
## 
## Inference[0] final GAMMA-based Likelihood: -2591.962819 tree written to file /working_dir/RAxML_result.aotus_cytb
## 
## 
## Starting final GAMMA-based thorough Optimization on tree 0 likelihood -2591.962819 .... 
## 
## Final GAMMA-based Score of best tree -2591.962819
## 
## Program execution info written to /working_dir/RAxML_info.
-
aotus_cytb
## Best-scoring ML tree written to: /working_dir/RAxML_bestTree.aotus_cytb
## 
## Overall execution time: 0.095162 secs or 0.000026 hours or 0.000001 days
## 
## 
```

# Visualisation

Now, let's check out our tree! We can use the R package
["Analysis of Phylogenetics and Evolution" or `ape`](http://ape-package.ird.fr/)
to do this.


```r
library(ape)
tree <- read.tree('RAxML_bestTree.aotus_cytb')
plot(tree)
```

![plot of chunk visualise](../figure/visualise-1.png)


---
title: "Install and run programs, outside of R, inside of R"
output: html_document
---



Many data science, statistical and academic projects require running models on
external software and then performing subsequent analyses of the results in R 
e.g. for project specific testing and visualisation. `outsider` aims to make
this process simpler by first enabling non-R software to be run from within R
and, second, by making it easier to install external program.

## How it works

The `outsider` package acts as an interface between the R environment and
external programs that are hosted on virtual machines. A virtual machine is
hosted on a user's computer but acts like an external computer with its own
operating system. These virtual machines are run through the program [Docker](https://www.docker.com/). So long as a computer is running Docker,
then any of these virtual machines can be downloaded and run without any
installation process. Docker runs on multiple operating systems including
Windows, OSX and Linux.

For every external program provided through `outsider`, a virtual machine, or
"Docker image", needs to be described and specific R code -- for launching and
interacting with the program -- is required. This Docker image and R code are
provided through **outsider modules** that are hosted on
[GitHub](https://github.com/).

Users can install any of the available `outsider` modules. With two commands,
a user can install and import an external program for calling within R:
`module_install()` and `module_import()`.

![outsider_outline](https://raw.githubusercontent.com/ropensci/outsider/master/other/outline.png)

---

## Installing `outsider`

Before you can make the most of `outsider` you will need to install and start
running Docker. Follow the installation instructions for your specific operating
system, ["Install Docker"](https://www.docker.com/products/docker-desktop).

> For some operating systems, "Docker Desktop" is not available. If that is the
case, try ["Docker Toolbox"](https://docs.docker.com/toolbox/). This is a
legacy Docker for older operating systems. It has similar functionality but
requires a virtual machine and has greater computational overhead.

With Docker installed, you then can install `outsider` via GitHub.


```r
library(remotes)
install_remotes('ropensci/outsider')
```

### Non-R dependencies

* For Windows users you may need to install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/).
* For Mac users you may need to install and setup 
"Xcode command-line tools", by running the command `xcode-select --install` in
the terminal.

## Finding and installing modules

To see what modules are available you can see the "[available modules page](https://docs.ropensci.org/outsider/articles/available.html)".
Alternatively, for the latest available information you can search for modules
using the `module_details()` function.


```r
library(outsider)
# repo = NULL will search for ALL available modules
#  (this may take a long time, depends on internet connection and remote server)
print(module_details(repo = 'dombennett/om..mafft'))
```

```
## # A tibble: 1 x 7
##   repo         program details                         versions updated_at          watchers_count url               
##   <chr>        <chr>   <chr>                           <chr>    <dttm>                       <int> <chr>             
## 1 dombennett/… mafft   Multiple alignment program for… latest   2020-01-16 10:36:00              0 https://github.co…
```

To install a module, all that is required is to provide the repo name to the
function `module_install()`.


```r
library(outsider)
module_install(repo = 'dombennett/om..mafft', force = TRUE)
```

> **What is `repo`?** The `repo` is the unique name for a GitHub repository that
hosts an `outsider` module. It consists of two parts: a GitHub username and a
project name. Given its uniqueness, all modules are referred to by their `repo`.

To confirm the module is installed on a computer, it might be useful to use
`module_installed()`. This function returns a table of all installed modules.


```r
library(outsider)
print(module_installed())
```

```
## # A tibble: 3 x 7
##   package         image                 tag    program      url                              image_created image_id  
##   <fct>           <chr>                 <chr>  <chr>        <chr>                            <chr>         <chr>     
## 1 om..hello.world dombennett/om_hello_… latest hello world  https://github.com/DomBennett/o… 8 months ago  acdff0a24…
## 2 om..mafft       dombennett/om_mafft   latest mafft        https://github.com/DomBennett/o… 13 months ago 97170a5f7…
## 3 om..partitionf… dombennett/om_partit… <NA>   PartitionFi… https://github.com/DomBennett/o… <NA>          <NA>
```

## Importing

All modules contain functions for interacting with the external program that
they host. To see these functions we can use `module_help()` to look up the help
documents.


```r
library(outsider)
# the whole module
module_help(repo = 'dombennett/om..mafft')
# specific function of a module (if known)
module_help(repo = 'dombennett/om..mafft', fname = 'mafft')
```

Once a function name is known of a particular module, the function can be
imported with `module_import()`.


```r
library(outsider)
mafft <- module_import(fname = 'mafft', repo = 'dombennett/om..mafft')
print(is(mafft))
```

```
## [1] "function"         "OptionalFunction" "PossibleMethod"
```

> **What is `mafft`?** "mafft" is a multiple alignment tool for for biological
sequences. Note, a user can use *any* name they wish for the function when it
is imported. For example, `mafftymcmafftface <- module_import( ...` would work
equally well.

## Commands

The imported functions from modules act like portals to the external programs
the modules host. To run a command, a user needs to use the function name and
give arguments corresponding to the arguments of the external program. For
example, on command-line to list the help information for `mafft`, we would
write `mafft --help`. With `outsider` we can do run `mafft('--help')`.

For a more complicated example, we could launch a small analysis with `mafft`
like so.


```r
library(outsider)
mafft <- module_import(fname = 'mafft', repo = 'dombennett/om..mafft')
mafft(arglist = c('--auto', 'input_sequences.fasta', '>',
                  'output_alignment.fasta'))
```

>**Why the spaces between arguments?** All the arguments of the external,
command-line program must be provided as separated characters. This helps
`outsider` parse the elements.

### How does that last line work?

```
mafft(arglist = c('--auto', 'input_sequences.fasta', '>',
'output_alignment.fasta'))
```
describes in R how to call the MAFFT program via command-line/terminal. It is
equivalent to `mafft --auto input_sequences.fasta > output_alignment.fasta` if
we were to call the program via command-line/terminal. How do we know how to
structure the program arguments? In the case of MAFFT we can look-up the
arguments on their website, [mafft.cbrc.jp](mafft.cbrc.jp/). But often for
command-line programs we can call for help with `-h` or `--help`. For MAFFT at
the command-line, we could run `mafft --help` or with `outsider` we can run:


```r
mafft(arglist = '--help')
```

```
## 
------------------------------------------------------------------------------
##   MAFFT v7.407 (2018/Jul/23)
## 
-
  https://mafft.cbrc.jp/alignment/software/
##   MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)
## ------------------------------------------------------------------------------
## High speed:
##   % mafft in > out
##   % mafft --retree 1 in > out (fast)
## 
## High accuracy (for <~200 sequences x <~2,000 aa/nt):
##   % mafft --maxiterate 1000 --localpair  in > out (% linsi in > out is also ok)
##   % mafft --maxiterate 1000 --genafpair  in > out (% einsi in > out)
##   % mafft --maxiterate 1000 --globalpair in > out (% ginsi in > out)
## 
## If unsure which option to use:
##   % mafft --auto in > out
## 
## --op # :         Gap opening penalty, default: 1.53
## --ep # :         Offset (works like gap extension penalty), default: 0.0
## --maxiterate # : Maximum number of iterative refinement, default: 0
## --clustalout :   Output: clustal format, default: fasta
## --reorder :      Outorder: aligned, default: input order
## --quiet :        Do not report progress
## --thread # :     Number of threads (if unsure, --thread -1)
## 
```

The help page returned tells us how to structure the arguments:

```
[options] [input_file] > [output_file]
```

Where the options (e.g. alignment method, number of threads) are always
indicated first with `--` and the input and output files are indicated second
with the `>`.

## Uninstalling

Clean up your computer by removing unwanted modules with `module_uninstall()`.

```r
library(outsider)
module_uninstall(repo = 'dombennett/om..mafft')
```

## Building your own module

Unfortunately, `outsider`'s utility is limited by the number of available
modules. Fortunately, it is very easy to create and upload your own module.
The package comes with a range of helper functions for minimising the amount of
coding for a module developer. If you know how to install an external program on
your own computer which you would like  would like to run it through `outsider`
and you have some experience with GitHub, then explore the
["outsider.devtools"](https://github.com/ropensci/outsider.devtools)
package.
---
title: "Run commands on a server"
output: rmarkdown::html_vignette
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, comment = '')
```

`outsider` allows users to launch commands on a remote machine. This can be
handy if the user wishes to run their external programs on more powerful
machines. All that is required is that the remote server has
[Docker](https://www.docker.com/) and can be connected to via
[ssh](https://en.wikipedia.org/wiki/Secure_Shell). To then run `outsider`
module programs on the remote machine, an `ssh` connection must be set up
using the `ssh_setup` function and the [`ssh`](https://github.com/ropensci/ssh)
package. Once set-up, all external programs will be run on the external
machine, not the local, with all input and output files being automatically sent
and retrieved by `outsider`.

```{r eval=FALSE}
# Libs
library(ssh)
library(outsider)

# SSH
session <- ssh_connect(host = "username@ip.address:port")
ssh_setup(session)

# Run
repo <- 'dombennett/om.hello.world'
if (!is_module_installed(repo = repo)) {
  module_install(repo = repo)
}
module_import(repo = repo, fname = 'hello_world')
# hello world will be run on the remote machine
hello_world()

# Clean up
ssh_teardown()
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/travis.R
\name{travis_build_status}
\alias{travis_build_status}
\title{Check Travis build status}
\usage{
travis_build_status(repo)
}
\arguments{
\item{repo}{repo}
}
\value{
Logical
}
\description{
Is build passing on travis? Returns either TRUE or FALSE.
}
\details{
For GitHub-based repositories only.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{github_search}
\alias{github_search}
\title{Search for outsider modules in GitHub}
\usage{
github_search()
}
\value{
data.frame
}
\description{
Returns GitHub API item results for outsider module search.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install.R
\name{module_import}
\alias{module_import}
\title{Import functions from a module}
\usage{
module_import(fname, repo)
}
\arguments{
\item{fname}{Function name to import}

\item{repo}{Module repo}
}
\value{
Function
}
\description{
Import specific functions from an outsider module to the
Global Environment.
}
\examples{
library(outsider)
if (is_outsider_ready()) {
  # simplest repo
  repo <- 'dombennett/om..hello.world'

  # is module_installed?
  if (is_module_installed(repo = repo)) {

    # get help for package
    module_help(repo = repo)
    
    # list functions available
    module_functions(repo = repo)
    
    # import
    hello_world <- module_import(fname = 'hello_world', repo = repo)
    
    # get help for function
    module_help(repo = repo, fname = 'hello_world')
    # also works
    ?hello_world

    # run function
    hello_world()
    
    # change verbosity settings
    
    # print nothing to console
    verbosity_set(show_program = FALSE, show_docker = FALSE)
    hello_world()
    
    # print everything to console
    verbosity_set(show_program = TRUE, show_docker = TRUE)
    hello_world()
    
    # write program output to a file
    log_file <- tempfile()
    verbosity_set(show_program = log_file, show_docker = FALSE)
    
    hello_world()
    
    (readLines(con = log_file))
    
    # Clean up
    file.remove(log_file)
  }
}
}
\seealso{
Other public: 
\code{\link{is_module_installed}()},
\code{\link{module_details}()},
\code{\link{module_functions}()},
\code{\link{module_help}()},
\code{\link{module_installed}()},
\code{\link{module_install}()},
\code{\link{module_search}()},
\code{\link{module_uninstall}()}
}
\concept{public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{yaml_fetch}
\alias{yaml_fetch}
\title{Safely fetch om.yaml}
\usage{
yaml_fetch(url)
}
\arguments{
\item{url}{URL to repo}
}
\value{
list
}
\description{
Return list of 'program' and 'details'.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gitlab.R
\name{gitlab_tags}
\alias{gitlab_tags}
\title{Module tags from GitLab}
\usage{
gitlab_tags(repo_ids)
}
\arguments{
\item{repo_ids}{Character vector of outsider module repositories ids.}
}
\value{
tbl_df
}
\description{
Return tbl_df of module tags for a list of outsider
modules hosted on gitlab.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install.R
\name{module_functions}
\alias{module_functions}
\title{List the functions associated with a module}
\usage{
module_functions(repo)
}
\arguments{
\item{repo}{Module repo}
}
\value{
character
}
\description{
Return a vector of functions that can be imported from the
module.
}
\examples{
library(outsider)
if (is_outsider_ready()) {
  # simplest repo
  repo <- 'dombennett/om..hello.world'

  # is module_installed?
  if (is_module_installed(repo = repo)) {

    # get help for package
    module_help(repo = repo)
    
    # list functions available
    module_functions(repo = repo)
    
    # import
    hello_world <- module_import(fname = 'hello_world', repo = repo)
    
    # get help for function
    module_help(repo = repo, fname = 'hello_world')
    # also works
    ?hello_world

    # run function
    hello_world()
    
    # change verbosity settings
    
    # print nothing to console
    verbosity_set(show_program = FALSE, show_docker = FALSE)
    hello_world()
    
    # print everything to console
    verbosity_set(show_program = TRUE, show_docker = TRUE)
    hello_world()
    
    # write program output to a file
    log_file <- tempfile()
    verbosity_set(show_program = log_file, show_docker = FALSE)
    
    hello_world()
    
    (readLines(con = log_file))
    
    # Clean up
    file.remove(log_file)
  }
}
}
\seealso{
Other public: 
\code{\link{is_module_installed}()},
\code{\link{module_details}()},
\code{\link{module_help}()},
\code{\link{module_import}()},
\code{\link{module_installed}()},
\code{\link{module_install}()},
\code{\link{module_search}()},
\code{\link{module_uninstall}()}
}
\concept{public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{module_details}
\alias{module_details}
\title{Look up details on module(s)}
\usage{
module_details(repo = NULL, service = c("github", "bitbucket", "gitlab"))
}
\arguments{
\item{repo}{Vector of one or more outsider module repositories, default NULL.}

\item{service}{Code-sharing service, e.g. GitHub}
}
\value{
tbl_df
}
\description{
Return a tbl_df of information for outsider module(s) for a
given code-sharing service. If \code{repo} is NULL, will return details on
all available modules.
}
\details{
Module details in tibble format include: repository name
(user/repo), last time repo was updated, number of watchers (or stars in the
case of GitLab), url to web presence, names of tagged versions.
}
\examples{
library(outsider)
# return table of ALL available modules on GitHub
# NOT RUN - takes too long
\dontrun{
  (available_modules <- module_search())
}

# look-up specific modules
repo <- 'dombennett/om..goodbye.world'
(suppressWarnings(module_details(repo = repo))) # no module exists, expect warning
repo <- 'dombennett/om..hello.world'
(module_details(repo = repo))
}
\seealso{
Other public: 
\code{\link{is_module_installed}()},
\code{\link{module_functions}()},
\code{\link{module_help}()},
\code{\link{module_import}()},
\code{\link{module_installed}()},
\code{\link{module_install}()},
\code{\link{module_search}()},
\code{\link{module_uninstall}()}
}
\concept{public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install.R
\name{is_module_installed}
\alias{is_module_installed}
\title{Is module installed?}
\usage{
is_module_installed(repo)
}
\arguments{
\item{repo}{Module repo}
}
\value{
Logical
}
\description{
Check if a module is installed on your system.
}
\details{
Searches for \code{repo} among installed outsider modules. Returns
TRUE if found, else FALSE.
}
\examples{
library(outsider)
# NOT RUN (too slow for automated testing)
\dontrun{
  if (is_outsider_ready()) {
    # simplest repo
    repo <- 'dombennett/om..hello.world'
    # install
    module_install(repo = repo, force = TRUE, update = 'never')
    # is module_installed?
    (is_module_installed(repo = repo))
    # uninstall
    module_uninstall(repo)
  }
}
}
\seealso{
Other public: 
\code{\link{module_details}()},
\code{\link{module_functions}()},
\code{\link{module_help}()},
\code{\link{module_import}()},
\code{\link{module_installed}()},
\code{\link{module_install}()},
\code{\link{module_search}()},
\code{\link{module_uninstall}()}
}
\concept{public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install.R
\name{module_uninstall}
\alias{module_uninstall}
\title{Uninstall and remove a module}
\usage{
module_uninstall(repo)
}
\arguments{
\item{repo}{Module repo}
}
\value{
Logical
}
\description{
Uninstall outsider module and removes it from your docker
}
\details{
If program is successfully removed from your system, TRUE is
returned else FALSE.
}
\examples{
library(outsider)
# NOT RUN (too slow for automated testing)
\dontrun{
  if (is_outsider_ready()) {
    # simplest repo
    repo <- 'dombennett/om..hello.world'
    # install
    module_install(repo = repo, force = TRUE, update = 'never')
    # is module_installed?
    (is_module_installed(repo = repo))
    # uninstall
    module_uninstall(repo)
  }
}
}
\seealso{
Other public: 
\code{\link{is_module_installed}()},
\code{\link{module_details}()},
\code{\link{module_functions}()},
\code{\link{module_help}()},
\code{\link{module_import}()},
\code{\link{module_installed}()},
\code{\link{module_install}()},
\code{\link{module_search}()}
}
\concept{public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{github_tags}
\alias{github_tags}
\title{Module tags from GitHub}
\usage{
github_tags(repos)
}
\arguments{
\item{repos}{Character vector of outsider module repositories.}
}
\value{
tbl_df
}
\description{
Return tbl_df of module tags for a list of outsider
modules hosted on GitHub.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bitbucket.R
\name{bitbucket_repo_search}
\alias{bitbucket_repo_search}
\title{Search for repository}
\usage{
bitbucket_repo_search(repo)
}
\arguments{
\item{repo}{bitbucket repo}
}
\value{
data.frame
}
\description{
Return bitbucket API item for specific repository.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{module_search}
\alias{module_search}
\title{Search for available outsider modules}
\usage{
module_search(service = c("github", "gitlab"))
}
\arguments{
\item{service}{Code-sharing service, e.g. GitHub}
}
\value{
Character vector
}
\description{
Return a list of available outsider modules. (Not possible for
BitBucket.)
}
\details{
Note: To search GitLab an access token is required. To create one:
1. Visit the personal access tokens section of your GitLab profile
\url{https://about.gitlab.com/}
2. Create a new token with api scope
3. Save the generated token to .Renviron
(try \code{usethis::edit_r_environ()}) with the line
"\code{GITLAB_PAT=your access token}"

For increased search relaiability, a token can be created for GitHub as well.
Visit \url{https://github.com/settings/tokens} to create a token and save
it to .Renviron as "GITHUB_PAT".
}
\examples{
library(outsider)
# return table of ALL available modules on GitHub
# NOT RUN - takes too long
\dontrun{
  (available_modules <- module_search())
}

# look-up specific modules
repo <- 'dombennett/om..goodbye.world'
(suppressWarnings(module_details(repo = repo))) # no module exists, expect warning
repo <- 'dombennett/om..hello.world'
(module_details(repo = repo))
}
\seealso{
Other public: 
\code{\link{is_module_installed}()},
\code{\link{module_details}()},
\code{\link{module_functions}()},
\code{\link{module_help}()},
\code{\link{module_import}()},
\code{\link{module_installed}()},
\code{\link{module_install}()},
\code{\link{module_uninstall}()}
}
\concept{public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bitbucket.R
\name{bitbucket_search}
\alias{bitbucket_search}
\title{Search for outsider modules in bitbucket}
\usage{
bitbucket_search(...)
}
\arguments{
\item{...}{Arguments}
}
\value{
data.frame
}
\description{
Returns bitbucket API item results for outsider module search.
}
\details{
Function is NOT available. This is a stub for when BitBucket API
updates.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repo.R
\name{pkgnm_guess}
\alias{pkgnm_guess}
\title{Guess package name}
\usage{
pkgnm_guess(repo, call_error = TRUE)
}
\arguments{
\item{repo}{Repository (e.g. "username/repo") or package name associated
with module}

\item{call_error}{Call error if no package found? Default, TRUE.}
}
\value{
character(1)
}
\description{
Return package name from a repo name.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gitlab.R
\name{gitlab_search}
\alias{gitlab_search}
\title{Search for outsider modules in GitLab}
\usage{
gitlab_search()
}
\value{
data.frame
}
\description{
Returns GitLab API item results for outsider module search.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outsider.R
\docType{package}
\name{outsider}
\alias{outsider}
\title{outsider: Install and run programs, outside of R, inside of R}
\description{
The outsider package facilitates the installation and running of external
software by interfacing with docker (\url{https://www.docker.com/}).
External software are contained within mini-R-packages, called "outsider
modules" and can be installed directly to a user's computer through online
code-sharing services such as GitHub (\url{https://github.com/}). The
outsider package comes with a series of functions for identifying,
installing and importing these outsider modules.
}
\details{
For more information visit the outsider website
(\url{https://docs.ropensci.org/outsider/}).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{yaml_read}
\alias{yaml_read}
\title{Module YAML information}
\usage{
yaml_read(repos, service = c("github", "gitlab", "bitbucket"))
}
\arguments{
\item{repos}{Character vector of outsider module repositories on GitHub,
GitLab or BitBucket.}

\item{service}{Code-sharing service. Character.}
}
\value{
tbl_df
}
\description{
Return tbl_df of all YAML information of given outsider
module repos.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outsider.R
\name{is_outsider_ready}
\alias{is_outsider_ready}
\title{Is outsider ready to run?}
\usage{
is_outsider_ready()
}
\value{
logical
}
\description{
Return TRUE if modules can be installed and run. Provides
helpful messages on what may be missing.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install.R
\name{module_help}
\alias{module_help}
\title{Get help for outsider modules}
\usage{
module_help(repo, fname = NULL)
}
\arguments{
\item{repo}{Module repo}

\item{fname}{Function name}
}
\description{
Look up help files for specific outsider module functions or
whole modules.
}
\examples{
library(outsider)
if (is_outsider_ready()) {
  # simplest repo
  repo <- 'dombennett/om..hello.world'

  # is module_installed?
  if (is_module_installed(repo = repo)) {

    # get help for package
    module_help(repo = repo)
    
    # list functions available
    module_functions(repo = repo)
    
    # import
    hello_world <- module_import(fname = 'hello_world', repo = repo)
    
    # get help for function
    module_help(repo = repo, fname = 'hello_world')
    # also works
    ?hello_world

    # run function
    hello_world()
    
    # change verbosity settings
    
    # print nothing to console
    verbosity_set(show_program = FALSE, show_docker = FALSE)
    hello_world()
    
    # print everything to console
    verbosity_set(show_program = TRUE, show_docker = TRUE)
    hello_world()
    
    # write program output to a file
    log_file <- tempfile()
    verbosity_set(show_program = log_file, show_docker = FALSE)
    
    hello_world()
    
    (readLines(con = log_file))
    
    # Clean up
    file.remove(log_file)
  }
}
}
\seealso{
Other public: 
\code{\link{is_module_installed}()},
\code{\link{module_details}()},
\code{\link{module_functions}()},
\code{\link{module_import}()},
\code{\link{module_installed}()},
\code{\link{module_install}()},
\code{\link{module_search}()},
\code{\link{module_uninstall}()}
}
\concept{public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gitlab.R
\name{gitlab_repo_search}
\alias{gitlab_repo_search}
\title{Search for repository}
\usage{
gitlab_repo_search(repo)
}
\arguments{
\item{repo}{gitlab repo}
}
\value{
data.frame
}
\description{
Return gitlab API item for specific repository.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bitbucket.R
\name{bitbucket_tags}
\alias{bitbucket_tags}
\title{Module tags from bitbucket}
\usage{
bitbucket_tags(repos)
}
\arguments{
\item{repos}{Character vector of outsider module repositories.}
}
\value{
tbl_df
}
\description{
Return tbl_df of module tags for a list of outsider
modules hosted on bitbucket.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/github.R
\name{github_repo_search}
\alias{github_repo_search}
\title{Search for repository}
\usage{
github_repo_search(repo)
}
\arguments{
\item{repo}{GitHub repo}
}
\value{
data.frame
}
\description{
Return GitHub API item for specific repository.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install.R
\name{module_install}
\alias{module_install}
\title{Install an outsider module}
\usage{
module_install(
  repo = NULL,
  url = NULL,
  filepath = NULL,
  git = NULL,
  service = c("github", "bitbucket", "gitlab"),
  tag = "latest",
  manual = FALSE,
  verbose = FALSE,
  force = FALSE,
  update = c("default", "ask", "always", "never")
)
}
\arguments{
\item{repo}{Module repo, character.}

\item{url}{URL to downloadable compressed (zip, tar or bzipped/gzipped)
folder of a module, character.}

\item{filepath}{Filepath to uncompressed directory of module, character.}

\item{git}{URL to git repository}

\item{service}{Code-sharing service. Character.}

\item{tag}{Module version, default latest. Character.}

\item{manual}{Build the docker image? Default FALSE. Logical.}

\item{verbose}{Be verbose? Default FALSE.}

\item{force}{Ignore warnings and install anyway? Default FALSE.}

\item{update}{Update dependent R packages?}
}
\value{
Logical
}
\description{
Install a module through multiple different methods: via a code
sharing site such as GitHub, a URL, a git repository or local filepath.
The function will first install the R package and then build the Docker
image. Docker image version is determined by "tag". To avoid pulling
the image from DockerHub set "manual" to TRUE.
}
\details{
All installation options depend on the installation functions of
\code{remotes}. E.g. GitHub packages are installed with
\code{\link[remotes]{install_github}}. See these functions for more details
on the R package installation process.
}
\examples{
library(outsider)
# NOT RUN (too slow for automated testing)
\dontrun{
  if (is_outsider_ready()) {
    # simplest repo
    repo <- 'dombennett/om..hello.world'
    # install
    module_install(repo = repo, force = TRUE, update = 'never')
    # is module_installed?
    (is_module_installed(repo = repo))
    # uninstall
    module_uninstall(repo)
  }
}
}
\seealso{
Other public: 
\code{\link{is_module_installed}()},
\code{\link{module_details}()},
\code{\link{module_functions}()},
\code{\link{module_help}()},
\code{\link{module_import}()},
\code{\link{module_installed}()},
\code{\link{module_search}()},
\code{\link{module_uninstall}()}
}
\concept{public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install.R
\name{module_installed}
\alias{module_installed}
\title{Which outsider modules are installed?}
\usage{
module_installed()
}
\value{
tbl_df
}
\description{
Returns tbl_df of details for all outsider modules
installed on the user's computer.
}
\examples{
library(outsider)
# check what modules are installed
if (is_outsider_ready()) {
  (module_installed())
}
}
\seealso{
Other public: 
\code{\link{is_module_installed}()},
\code{\link{module_details}()},
\code{\link{module_functions}()},
\code{\link{module_help}()},
\code{\link{module_import}()},
\code{\link{module_install}()},
\code{\link{module_search}()},
\code{\link{module_uninstall}()}
}
\concept{public}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outsider.R
\name{ssh_teardown}
\alias{ssh_teardown}
\title{Teardown SSH}
\usage{
ssh_teardown()
}
\value{
logical
}
\description{
Disconnect from a remote host and stop commands being
transferred.
}
\examples{
library(outsider)

# To forward all Docker commands to a remote host:
# 1. Gain ssh access to a remote host
# 2. Ensure Docker is running on the remote machine
# 3. Supply the IP address and authentication args to ssh::ssh_connect
ip_address <- NULL

if (!is.null(ip_address)) {
  # Create an ssh session
  session <- ssh::ssh_connect(host = ip_address)
  
  # Setup the session for running outsider
  ssh_setup(session)
}

# After setup, run outsider as normal

# simplest repo
repo <- 'dombennett/om..hello.world'

if (is_module_installed(repo = repo)) {
  # import
  hello_world <- module_import(fname = 'hello_world', repo = repo)
  
  # run function
  hello_world()
}

# Always ensure to disconnect after a session
if (!is.null(ip_address)) {
  ssh_teardown()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install.R
\name{user_warn}
\alias{user_warn}
\title{Warn users on the dangers of outsider modules}
\usage{
user_warn(pkgnm)
}
\arguments{
\item{pkgnm}{Package name}
}
\value{
Logical
}
\description{
Warn users on the dangers of installing an outsider module
whose origin is potentially unknown.
}
\details{
Prints additional info to screen based on module YAML file.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outsider.R
\name{ssh_setup}
\alias{ssh_setup}
\title{Setup SSH}
\usage{
ssh_setup(session)
}
\arguments{
\item{session}{\code{ssh} session, see \code{\link[ssh]{ssh_connect}}}
}
\value{
logical
}
\description{
Send all outsider commands to an external host. Provide an
\code{ssh} session to this function and all subsequent commands will be run
on the host rather than the local machine. When finished it is always good
practice to disconnect from the remote host by running \code{ssh_teardown}.
It is required that the remote host has Docker running.
}
\examples{
library(outsider)

# To forward all Docker commands to a remote host:
# 1. Gain ssh access to a remote host
# 2. Ensure Docker is running on the remote machine
# 3. Supply the IP address and authentication args to ssh::ssh_connect
ip_address <- NULL

if (!is.null(ip_address)) {
  # Create an ssh session
  session <- ssh::ssh_connect(host = ip_address)
  
  # Setup the session for running outsider
  ssh_setup(session)
}

# After setup, run outsider as normal

# simplest repo
repo <- 'dombennett/om..hello.world'

if (is_module_installed(repo = repo)) {
  # import
  hello_world <- module_import(fname = 'hello_world', repo = repo)
  
  # run function
  hello_world()
}

# Always ensure to disconnect after a session
if (!is.null(ip_address)) {
  ssh_teardown()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outsider.R
\name{verbosity_set}
\alias{verbosity_set}
\title{Set the verbosity of modules}
\usage{
verbosity_set(show_program = TRUE, show_docker = FALSE)
}
\arguments{
\item{show_program}{Show external program messages? Default TRUE.}

\item{show_docker}{Show docker messages? Default FALSE.}
}
\value{
data.frame
}
\description{
Control console messages of running outsider modules. Allow
either the external program messages to run, the Docker messages or both.
}
\details{
For more control see \code{\link[outsider.base]{log_set}}
}
\examples{
library(outsider)
# NOT RUN (too slow for automated testing)
\dontrun{
  if (is_outsider_ready()) {
    # simplest repo
    repo <- 'dombennett/om..hello.world'
    # install
    module_install(repo = repo, force = TRUE, update = 'never')
    # is module_installed?
    (is_module_installed(repo = repo))
    # uninstall
    module_uninstall(repo)
  }
}
}
