
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



