
Development tools package for `outsider` <img src="https://raw.githubusercontent.com/ropensci/outsider.devtools/master/other/logo_devtools.png" height="200" align="right"/>
----

[![Build Status](https://travis-ci.org/ropensci/outsider.devtools.svg?branch=master)](https://travis-ci.org/ropensci/outsider.base) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/outsider.devtools?branch=master&svg=true)](https://ci.appveyor.com/project/DomBennett/outsider-devtools) [![Coverage Status](https://coveralls.io/repos/github/ropensci/outsider.devtools/badge.svg?branch=master)](https://coveralls.io/github/ropensci/outsider.devtools?branch=master) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3615074.svg)](https://doi.org/10.5281/zenodo.3615074)

## Build [`outsider`](https://github.com/ropensci/outsider#readme) modules!

This package aims to make life easier for  those who wish to build their own
`outsider` modules. In just a few function calls: build a module skeleton,
check the file structures, test the module and upload it online!
Top banana! :banana:

Acquaint yourself better with all these steps by reading up on:

* ["Basic: Module Building"](https://docs.ropensci.org/outsider.devtools/articles/basic.html)
* ["Intermediate: Module Building"](https://docs.ropensci.org/outsider.devtools/articles/intermediate.html)
* ["Advanced: Module Building"](https://docs.ropensci.org/outsider.devtools/articles/advanced.html)

Happy building! :wrench:

## Install

Install via GitHub ....

```r
# install.packages("remotes")
remotes::install_github("ropensci/outsider.devtools")
```

In addition to installing `outsider.devtools`, the above code will also install
the key dependency packages
[`outsider.base`](https://github.com/ropensci/outsider.base#readme) and 
[`outsider`](https://github.com/ropensci/outsider#readme). (Read up on
[`remotes`](https://github.com/r-lib/remotes))

## Quick guide

Build an [`outsider`](https://github.com/ropensci/outsider#readme) module to run [`echo`](https://en.wikipedia.org/wiki/Echo_(command)) via the [Linux distribution Ubuntu](https://en.wikipedia.org/wiki/Ubuntu) in just a few function calls.

```r
# make my own quick package
library(outsider.devtools)

# construct a skeleton file structure for the module
module_path <- module_skeleton(program_name = 'echo', flpth = getwd())

# check the file structure
module_check(flpth = module_path)

# look-up key identifying names: R package name, Docker image name
module_identities(flpth = module_path)

# build the R package and Docker image
module_build(flpth = module_path, tag = 'latest')

# test the module
module_test(flpth = module_path)
```

![](https://raw.githubusercontent.com/ropensci/outsider.devtools/master/other/build_example.gif)

Visit the webpage ["The Basics"](https://docs.ropensci.org/outsider.devtools/articles/basic.html) to find out more.

## How do the `outsider` packages interact?

![](https://raw.githubusercontent.com/ropensci/outsider.devtools/master/other/package_structures.png)

## Citation

Bennett et al., (2020). outsider: Install and run programs, outside of R, inside of R. Journal of Open Source Software, 5(45), 2038, https://doi.org/10.21105/joss.02038

## Maintainer

[Dom Bennett](https://github.com/dombennett/)

---

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# outsider 0.1.1

## Accepted ROpenSci version

* Improved documentation
* Reexports from base
* Improved skeleton
* Better tests

# outsider 0.1.0

## Initial version

* Develop modules
---
output: github_document
---
<!--
The README should be used to describe the program. It acts like the homepage of
your module.

Edit README.Rmd not README.md. The .Rmd file can be knitted to parse real-code
examples and show their output in the .md file.

To knit, use devtools::build_readme() or outsider.devtools::build()

Edit the template to describe your program: how to install, import and run;
run exemplary, small demonstrations; present key arguments; provide links and
references to the program that the module wraps.

Learn more about markdown and Rmarkdown:
https://daringfireball.net/projects/markdown/syntax
https://rmarkdown.rstudio.com/
-->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# Run `%program_name%` with `outsider` in R

[![Build Status](https://travis-ci.org/%repo%.svg?branch=master)](https://travis-ci.org/%repo%)

> Summary line of what the program does.


<!-- Install information -->
## Install and look up help

```{r install-snippet, eval=FALSE, include=TRUE}
library(outsider)
module_install(repo = "%repo%")
module_help(repo = "%repo%")
```

<!-- Detailed examples -->
## A detailed example

<!-- Note: set eval=TRUE to run example and show output -->
```{r detailed-example, eval=FALSE, include=TRUE}
library(outsider)
%cmd% <- module_import(fname = '%cmd%', repo = "%repo%")
%cmd%(arglist = c('--a', '-more=10', '-complicated=a', 'example.txt'))
```

<!-- Remove module after running above example -->
```{r uninstall-snippet, eval=FALSE, include=FALSE}
module_uninstall(repo = "%repo%")
```

### Key arguments

Describe key arguments and show program output. E.g. in table format.

|Argument|Usage|Description|
|--------|-----|-----------|
|a       |--a  |Runs "a" method|
|more    |-more=#  |Run # more|
|complicated|-complicated=[a-z]|Run [a-z]|
|Output file|Final argument|Where to storeoutput|

## Links

Find out more by visiting the
[%program_name%'s homepage](www.external_program.org).

## Please cite

* Smith, J. et al. (2020) %program_name% reference. *Journal of Outsider
Modules*
* Bennett et al. (2020). outsider: Install and run programs, outside of R,
inside of R. *Journal of Open Source Software*, In review


<!-- Footer -->
---

<img align="left" width="120" height="125" src="https://raw.githubusercontent.com/ropensci/outsider/master/logo.png">

**An `outsider` module**

Learn more at [outsider website](https://docs.ropensci.org/outsider/).
Want to build your own module? Check out [`outsider.devtools` website](https://docs.ropensci.org/outsider.devtools/).
---
title: "Frequently Asked Questions"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


## Why is my module not returned by `outsider::module_search`?

`outsider::module_search` makes it easier for users to locate and install
modules. It can search for modules across both GitHub and GitLab (but not
BitBucket). For both services, the repo names must begin `om..` and there must
be an `om.yml` in the top-level directory.


Additionally, GitHub repos must have **"outsider-module"** at the beginning of
their repo descriptions.


Try comparing and contrasting your set-up with these repos:

* [Working example on GitHub](https://github.com/DomBennett/om..hello.world)
* [Working example on GitLab](https://gitlab.com/DomBennett/om..hello.world)

## Can a module be an R package?

There is nothing preventing you develop an R package while enabling the
package to also act as an `outsider` module. You do not need to call the package
`om..`, all that is required is the `om.yml` and a functioning Dockerfile.

## How do I delete an old image?

Docker maintains images and containers. In the advanced vignette we showed you
how to stop and remove containers. To delete unwanted images:

```bash
# list images
docker image ls
# delete an image
docker image rm [IMAGE ID]
```
---
title: "Advanced: Building a module"
output: html_document
---



In [intermediate](https://docs.ropensci.org/outsider.devtools/articles/intermediate.html)
we demonstrated how a complex, command-line program with input and output files
(that we partly developed) can be cast as an `outsider` module.

On this page, we will develop a module for a typical command-line program with
lots of arguments.

## [`RAxML`](https://cme.h-its.org/exelixis/web/software/raxml/)

For this walkthrough we will create a module for "RAxML" a program
that can generate evolutionary trees using a maximum-likelihood. The exact
functioning (or use even!) of the program is not relevant, but it demonstrates
well the typical program which `outsider` aims to bring into the R environment:
lots of arguments, complex algorithm developed in an alternative language.


To demonstrate the program's complexity, here is the command's syntax:

```
raxmlHPC[-SSE3|-AVX|-PTHREADS|-PTHREADS-SSE3|-PTHREADS-AVX|-HYBRID|-HYBRID-SSE3|HYBRID-AVX]
      -s sequenceFileName -n outputFileName -m substitutionModel
      [-a weightFileName] [-A secondaryStructureSubstModel]
      [-b bootstrapRandomNumberSeed] [-B wcCriterionThreshold]
      [-c numberOfCategories] [-C] [-d] [-D]
      [-e likelihoodEpsilon] [-E excludeFileName]
      [-f a|A|b|B|c|C|d|D|e|E|F|g|G|h|H|i|I|j|J|k|m|n|N|o|p|P|q|r|R|s|S|t|T|u|v|V|w|W|x|y] [-F]
      [-g groupingFileName] [-G placementThreshold] [-h] [-H]
      [-i initialRearrangementSetting] [-I autoFC|autoMR|autoMRE|autoMRE_IGN]
      [-j] [-J MR|MR_DROP|MRE|STRICT|STRICT_DROP|T_<PERCENT>] [-k] [-K] 
      [-L MR|MRE|T_<PERCENT>] [-M]
      [-o outGroupName1[,outGroupName2[,...]]][-O]
      [-p parsimonyRandomSeed] [-P proteinModel]
      [-q multipleModelFileName] [-r binaryConstraintTree]
      [-R binaryModelParamFile] [-S secondaryStructureFile] [-t userStartingTree]
      [-T numberOfThreads] [-u] [-U] [-v] [-V] [-w outputDirectory] [-W slidingWindowSize]
      [-x rapidBootstrapRandomNumberSeed] [-X] [-y] [-Y quartetGroupingFileName|ancestralSequenceCandidatesFileName]
      [-z multipleTreesFile] [-#|-N numberOfRuns|autoFC|autoMR|autoMRE|autoMRE_IGN]
      [--mesquite][--silent][--no-seq-check][--no-bfgs]
      [--asc-corr=stamatakis|felsenstein|lewis]
      [--flag-check][--auto-prot=ml|bic|aic|aicc]
      [--epa-keep-placements=number][--epa-accumulated-threshold=threshold]
      [--epa-prob-threshold=threshold]
      [--JC69][--K80][--HKY85]
      [--bootstop-perms=number]
      [--quartets-without-replacement]
      [---without-replacement]
      [--print-identical-sequences]
```

How are we going to code a module that is able to account for all these
arguments?

# Walkthrough

## Building

As before, let's just jump right in with the module skeleton.


```r
library(outsider.devtools)
flpth <- module_skeleton(repo_user = 'dombennett', program_name = 'raxml',
                         docker_user = 'dombennett', flpth = tempdir(),
                         service = 'github')
# folder name where module is stored
print(basename(flpth))
```

```
## [1] "om..raxml"
```

Next, the Dockerfile. This is a much more complicated installation process.
First we require the apt programs `wget`, `make` and `gcc` to download the
program and compile it. Then we need to download and compile the program.
The majority of the Dockerfile contents can be found in the installation
instructions of the
[RAxML GitHub site](https://github.com/stamatak/standard-RAxML#readme).


```r
dockerfile_text <- "
FROM ubuntu:latest
RUN apt-get update && apt-get install -y \
    wget make gcc
# 8.2 for this demonstration
RUN wget https://github.com/stamatak/standard-RAxML/archive/v8.2.12.tar.gz && \
    tar zxvf v8.2.12.tar.gz && \
    rm v8.2.12.tar.gz && \
    mv standard-RAxML-8.2.12 raxml
RUN cd /raxml && make -f Makefile.SSE3.PTHREADS.gcc && \
    rm *.o && cp raxmlHPC-PTHREADS-SSE3 /usr/bin/.
RUN mkdir /working_dir
WORKDIR /working_dir
"
# write to latest Dockerfile
cat(dockerfile_text, file = file.path(flpth, 'inst', 'dockerfiles', 'latest',
                                      'Dockerfile'))
```

Dockerfile defined. What should our R function be?


```r
#' @name raxml
#' @title raxml
#' @description Run raxml
#' @param arglist Arguments to raxml provided as a character vector
#' @param outdir Filepath to where all output files should be returned.
#' @example examples/example.R
#' @export
raxml <- function(arglist = arglist_get(...), outdir = getwd()) {
  files_to_send <- filestosend_get(arglist)
  arglist <- arglist_parse(arglist)
  otsdr <- outsider_init(pkgnm = 'om..raxml', cmd = 'raxmlHPC-PTHREADS-SSE3',
                         wd = outdir, files_to_send = files_to_send,
                         arglist = arglist)
  run(otsdr)
}
```

In this function, input and output files are intermixed in the arguments passed
to RAxML. To determine what is a file that should be sent to the container from
a normal argument, we have `filestosend_get`. This function goes through all the
provided arguments to check for any files (it tests whether they are valid file
paths) and returns them as a character vector to `files_to_send`. Then we have
`arglist_parse`, this converts all given arguments to their basename so that
absolute filepaths are dropped. All input files will be passed to the containers
`working_dir/` and therefore no parental file path are required. This function
has additional functions for dropping irrelevant arguments in the context of an
`outsider` module (e.g. no need to specify an output directory at the 
command-line, this would need to be done at the initiation of the `outsider`
object.)

> **What's with this `outdir`-thingy?** The RAxML function returns a series of
files. To ensure they are all returned in a convenient location, we have
created an additional argument to the `raxml()` function called `outdir`. It can
be good practice to break-up information that needs to be processed by the
module function into things like `arglist` and `outdir`. In some instances, you
may be calling a program where additional information can be provided outside of
remit of the main program being emulated. For example in our last example we had
`input_file` and `output_file`. Other examples include, number of threads,
memory allocation for the program, etc.

OK! So that's what the function looks like. Let's write it to file.


```r
function_text <- "
#' @name raxml
#' @title raxml
#' @description Run raxml
#' @param arglist Arguments to raxml provided as a character vector
#' @param outdir Filepath to where all output files should be returned.
#' @example examples/example.R
#' @export
raxml <- function(arglist = arglist_get(...), outdir = getwd()) {
  files_to_send <- filestosend_get(arglist)
  arglist <- arglist_parse(arglist)
  otsdr <- outsider_init(pkgnm = 'om..raxml', cmd = 'raxmlHPC-PTHREADS-SSE3',
                         wd = outdir, files_to_send = files_to_send,
                         arglist = arglist)
  run(otsdr)
}
"
# write to R/functions.R
cat(function_text, file = file.path(flpth, 'R', 'functions.R'))
```

Let's check the setup.

Ok, let's build, test and upload the module! (The test should pass because by
default, `-h` is used in the example and this is a valid argument for RAxML.)


```r
module_build(flpth = flpth)
```

```
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::document() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```
## Updating om..raxml documentation
```

```
## First time using roxygen2. Upgrading automatically...
```

```
## Updating roxygen version in /private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..raxml/DESCRIPTION
```

```
## Writing NAMESPACE
```

```
## Loading om..raxml
```

```
## Writing NAMESPACE
## Writing raxml.Rd
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::install() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##   
   checking for file ‘/private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..raxml/DESCRIPTION’ ...
  
✓  checking for file ‘/private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..raxml/DESCRIPTION’
## 
  
─  preparing ‘om..raxml’:
## 
  
   checking DESCRIPTION meta-information ...
  
✓  checking DESCRIPTION meta-information
## 
  
─  checking for LF line-endings in source and make files and shell scripts
## 
  
─  checking for empty or unneeded directories
## 
  
─  building ‘om..raxml_0.0.1.tar.gz’
## 
  
   
## 
Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
##   /var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T//RtmpquGDgt/om..raxml_0.0.1.tar.gz --install-tests 
## * installing to library ‘/Library/Frameworks/R.framework/Versions/3.6/Resources/library’
## 
-
* installing *source* package ‘om..raxml’ ...
## 
|
** using staged installation
## 
-
** R
## 
|
** inst
## 
-
** byte-compile and prepare package for lazy loading
## 
|

-
** help
## 
|
*** installing help indices
## 
-
** building package indices
## 
|
** testing if installed package can be loaded from temporary location
## 
-
** testing if installed package can be loaded from final location
## 
|
** testing if installed package keeps a record of temporary installation path
## 
-
* DONE (om..raxml)
## 
|

-

 
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running docker_build()
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Command:
## docker build -t dombennett/om_raxml:latest /Library/Frameworks/R.framework/Versions/3.6/Resources/library/om..raxml/dockerfiles/latest
## .....................................................................................................................
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::build_readme() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```
## Building om..raxml readme
```

```r
module_test(flpth = flpth)
```

```
## Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
##   /private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..raxml --install-tests 
## * installing to library ‘/Library/Frameworks/R.framework/Versions/3.6/Resources/library’
## 
-
* installing *source* package ‘om..raxml’ ...
## 
|
** using staged installation
## 
-
** R
## 
|
** inst
## 
-
** byte-compile and prepare package for lazy loading
## 
|

-
** help
## 
|
*** installing help indices
## 
-
** building package indices
## 
|
** testing if installed package can be loaded from temporary location
## 
-
** testing if installed package can be loaded from final location
## 
|
** testing if installed package keeps a record of temporary installation path
## 
-
* DONE (om..raxml)
## 
|

-

 
```

```
## No Docker image found for 'om..raxml' -- attempting to pull/build image with tag 'latest'
```

```
## Removing package from '/Library/Frameworks/R.framework/Versions/3.6/Resources/library'
## (as 'lib' is unspecified)
```

```
## The gloating goat got away! The module is not working....
## But keep on forming! You're doing correctly!
```

And then upload ....


```r
module_upload(flpth = flpth, code_sharing = TRUE, dockerhub = TRUE)
```

### Continuous Integration (GitHub only)

OK. So you have a module that works and passes its tests. How do we inform the
world that it works? The answer: [Travis-CI](https://travis-ci.org/)

Travis is a free service that monitors your GitHub for updates and then runs
a series of tests on a remote server to determine whether it passes its tests or
not. The result of these tests can be added to a repo's README in form of an
image badge: e.g. "Tests Passing".

To set-up travis for your repo, visit the website https://travis-ci.org/, sign
in with your GitHub account and then activate the service your repo of choice,
e.g. `om..raxml`.

Then we need to provide the test instructions for Travis in the form of a YAML
file called `.travis.yml` and upload this to our GitHub repo. The `.travis.yml`
can be created for `outsider` modules so:


```r
module_travis(flpth = flpth)
```

Then simply update the repo on GitHub!

## Taking the module forward

After you have built a module that passes its tests and is available for
download, there are a few extra things to consider:

* **Update the `om.yml`** Add a more detailed description
* **Update the README** The default README is very bland. It provides a template
for what sorts of things a good README should contain: installation, how to run,
where to get more help. Remember to only edit the `README.Rmd`. To update the
`README.md` for GitHub, you need to "knit" the `.Rmd` file by running
`devtools::build_readme()` -- this is the function `module_build()` runs.
* **Ensure discoverability** `outsider::module_search` can find `outsider`
module repos on GitHub if they have and `om.yml`, have their name in the form
`om..` and have the phrase `outsider-module: ` in their description.

---

**Delete it all**


```r
module_uninstall(repo = 'om..raxml')
unlink(x = flpth, recursive = TRUE, force = TRUE)
```

---

# Tips and tricks

* The best way to learn how to build your own module is to look at how others
have created modules for programs you are familiar with. Check out the list of
already available modules at this
[page](https://docs.ropensci.org/outsider/articles/available.html).
* If you find issues with the docker image build step, you can try "logging
into" a container and running the installation steps to test out what works.
At the terminal, use


```bash
# run an ubuntu container
docker run -t -d --name test ubuntu:latest
# iteractively run bash
docker exec -it test bash
# when finished, list all running containers and stop/rm
docker ps -a
docker stop test
docker rm test
```
---
title: "Basic: Module Building"
output: html_document
---



`outsider` allows users to install and run programs on their own computer
without the need to leave the R environment. If there is a command-line program
that is not available through `outsider` you can create your own module!

If you are able to install the program on your own machine and have some
experience with R packages and GitHub then you should be readily able to
create one.

On this page, we will outline what an `outsider` module is and then provide
a simple walkthrough for creating the simplest of simple modules: a command line
program that prints whatever text the user provides it.

To follow this guide you will need the following:

* A GitHub account (sign-up [here](https://github.com/))
* A Docker-Hub account (sign-up [here](https://hub.docker.com/))

---

# Basics

## Module structure

At its heart, an `outsider` module is just an R package with a `Dockerfile`.
All modules have the following basic file/folder structure:

```
# Core files and folders of a module #
- DESCRIPTION
- R/
    - functions.R
    - import.R
- inst/
    - dockerfiles/
        - latest/
            - Dockerfile
    - om.yaml
- README.md
- README.Rmd
- .travis.yml
```

### R files and folders

The R package is determined by the `DESCRIPTION` file and the `R/` folder.
The first of these describes the package details (package name, author,
dependencies, etc.) and the second contains the R code that makes up the
package.
The `R/` folder can have any number of scripts with whatever names a developer
chooses (skeleton default `functions.R` and `import.R`).

### Docker files and folders

The `outsider` package and its modules depend on Docker to run. Docker is a
service that is able to host isolated software packages, termed "containers",
that act _like_ virtual machines but, thanks to
[OS-level virtualisation](https://en.wikipedia.org/wiki/OS-level_virtualization),
are in fact much more lightweight. By using a Docker container a user is able to
run specific code, often designed for different operating systems, on any
machine that has Docker installed. A Docker container is created by
running a Docker "image" where the image acts as a description of the code
necessary for the container to function (base operating system, required
programs, etc.). And these images are described by a `Dockerfile`.


|Docker term|Description|
| --------- |:----------|
|Container|Isolated environment where applications are hosted and launched|
|Image|A file that acts as the "blueprint" for containers|
|Dockerfile|Text-based file describing the steps (layers) that define an image|


In an `outsider` module, the `dockerfiles/` folder in the module contains the
`Dockerfile` that describes how to install the external command-line program
that the developer wishes to run through R. `dockerfiles/` can have multiple
versions of a `Dockerfile` but every module must have a `latest/`.

> **Why "latest"?** Every Docker image is tagged with a version name/number.
For `outsider` modules all images must have a "latest" tag which acts as the
default tag. Additionally, a developer may add any number of additional tagged
versions of their program (e.g.legacy version of a program) to their module by
creating new Dockerfiles in separate directories under `dockerfiles/`.

For more information on Docker, see:

* [Official Docker Docs](https://docs.docker.com/)
* [Docker glossary](https://docs.docker.com/glossary/)
* [Learn Docker Basics](https://docker-curriculum.com/)
* [Wikipedia Page](https://en.wikipedia.org/wiki/Docker_(software))

### GitHub files and folders

To make modules discoverable on GitHub, all modules require an `om.yml`.
This file has two elements (program and details) encoded in the
[YAML format](https://en.wikipedia.org/wiki/YAML). In addition, the module
has a `README.md` file that provides the text describing the module on the
module's GitHub (or other code-sharing site) homepage. `.md` format is like a
simplified HTML. It can be generated by hand or -- by default with
`outsider.devtools` -- it can be rendered from a `.Rmd` file, in this case
`README.Rmd`. The advantage of the `.Rmd`-approach is that R code chunks are
parsed and their output is then included in the `.md` version of the file. This
allows users to better understand the code as they can directly see the output.


Finally, a `.travis.yml` provides instructions detailing how the module should
be tested on [Travis-CI](https://en.wikipedia.org/wiki/Travis_CI). By default,
the `outsider.devtools` helper functions create a `.travis.yml` file that tests
by installing the module and then running the examples of the exported R
functions.

---

# Walkthrough

We will walk you through how to create your own `outsider` module that will
simply print (through `echo`) any text provided. This process comes in a
series of steps:

i. Generate core files and folders
ii. Create Docker image
iii. Document and build the R package
iv. Try the module
v. Upload to GitHub and Docker Hub

> **What's "Docker Hub"?** Docker Hub is service that hosts Docker images. It
can simplify Dockerfile creation as image layers can be sourced from
pre-existing images available via Docker Hub. Additionally, by sharing your
images generated from your module Dockerfile on Docker Hub, you will be
speeding up the installation step for end-users.

## Generate the files and folders

As displayed above, we need to generate the core files and folders of a module.
This process can be easily performed using the `module_skeleton()` function.
This function takes a few details about the developer and the program and then
generates all the necessary core file structures.

The necessary information required for the module to run are our GitHub and
Docker Hub usernames plus the name for the program we wish to provide as a
module (which is `echo` -- a UNIX command for printing). In these code snippets,
the usernames are of the `outsider` maintainer. In order for these examples to
work for you, you will need to change "dombennett" to your own usernames.
(Note, your usernames may differ for GitHub and Docker Hub).


```r
library(outsider.devtools)
# the file location of where the module is saved is returned to "module_path"
module_path <- module_skeleton(repo_user = 'dombennett', program_name = 'echo',
                               docker_user = 'dombennett', flpth = tempdir(),
                               full_name = 'D.J. Bennett',
                               email = 'dominic.john.bennett@gmail.com',
                               service = 'github')
# folder name where module is stored
print(basename(module_path))
```

```
## [1] "om..echo"
```

The above code will create an `outsider` module with the module and directory
name `om..echo` at the file location `module_path`.
After running the above code, you should take a minute to inspect the generated
files. In particular you should look at the `DESCRIPTION` file, the `Dockerfile`
and `functions.R`.

**File tree `module_path` (before build)**


```
## ├── DESCRIPTION
## ├── R
## │   ├── functions.R
## │   └── import.R
## ├── README.Rmd
## ├── examples
## │   └── example.R
## └── inst
##     ├── dockerfiles
##     │   └── latest
##     │       └── Dockerfile
##     └── om.yml
```

> **Why `om..echo`?** All modules must start with "om.." in order for them to
be discovered on GitHub. This is not a requirement for the functioning of the
module, it just allows `outsider::module_search()` to find them.

At this stage, we would then edit the `inst/dockerfiles/latest/Dockerfile` and
the `R/functions.R` to work for our chosen external program. But because `echo`
-- the program we wish to port through `outsider` -- is so simple we don't
actually have to make any changes to these files. By default, the `Dockerfile`
is based on Ubuntu which ships with `echo` and our starting function in
`functions.R` is based around running `echo`: it parses the arguments, creates
an outsider object, and then launches the object.


```r
# The echo function in om..echo
echo <- function(...) {
  # convert the ... into an argument list
  arglist <- arglist_get(...)
  # create an outsider object: describe the arguments and program
  otsdr <- outsider_init(pkgnm = 'om..echo',
                         cmd = 'echo', arglist = arglist)
  # run the command
  run(otsdr)
}
```

> **What's `...`?** In function calls in R, `...` indicate that any number of
arguments can be provided to a function. The `arglist_get` function will take
the `...` and convert them into a character vector that can be parsed.
`outsider` recommends module functions make use of this feature so that any
number of arguments can be passed to external programs. Additionally, this has
the advantage that the developer would then not need to document all the
arguments of the external program. For many external programs there may be
hundreds of arguments, all of which are likely to be already documented;
viewable through commands like `-h` or `--help`. For programs with few arguments
or where the execution of the external program requires additional setting-up
(input files, environment settings, etc.) or  where argument definitions cannot
be displayed with `-h` or `--help` then it is best to provide richer
documentation at the R level. In these instances, the specific arguments to
the external program can either be separate R arguments, `function(arg1, arg2)`,
or they can be a character vector along with other arguments,
`function(arglist=arglist_get(...), input=NULL, memmory="1GB")`.

## Building the R package

With the skeleton set-up, we can now install the R package using
`module_build()`. At this stage we only want to build the package components,
not the Docker image, so we can set all build options to TRUE except for the
`build_image` option.


```r
module_build(flpth = module_path, build_documents = TRUE, build_package = TRUE,
             build_image = FALSE, build_readme = TRUE)
```

```
## Updating om..echo documentation
```

```
## First time using roxygen2. Upgrading automatically...
```

```
## Updating roxygen version in /private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..echo/DESCRIPTION
```

```
## Loading om..echo
```

```
## Building om..echo readme
```

```
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::document() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Writing NAMESPACE
## Writing NAMESPACE
## Writing echo.Rd
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::install() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##   
   checking for file ‘/private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..echo/DESCRIPTION’ ...
  
✓  checking for file ‘/private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..echo/DESCRIPTION’
## 
  
─  preparing ‘om..echo’:
## 
  
   checking DESCRIPTION meta-information ...
  
✓  checking DESCRIPTION meta-information
## 
  
─  checking for LF line-endings in source and make files and shell scripts
## 
  
─  checking for empty or unneeded directories
## 
  
─  building ‘om..echo_0.0.1.tar.gz’
## 
  
   
## 
Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
##   /var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T//RtmpquGDgt/om..echo_0.0.1.tar.gz --install-tests 
## 
-
* installing to library ‘/Library/Frameworks/R.framework/Versions/3.6/Resources/library’
## 
|
* installing *source* package ‘om..echo’ ...
## 
-
** using staged installation
## 
|
** R
## 
-
** inst
## 
|
** byte-compile and prepare package for lazy loading
## 
-

|
** help
## 
-
*** installing help indices
## 
|
** building package indices
## 
-
** testing if installed package can be loaded from temporary location
## 
|
** testing if installed package can be loaded from final location
## 
-
** testing if installed package keeps a record of temporary installation path
## 
|
* DONE (om..echo)
## 
-

|

 
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::build_readme() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

> **What does `build_documents = TRUE` do?** In addition to the core files
indicated above, an R package also requires R documentation files that are
stored in `man/` -- these provide the `?[function]` utility. The above
function generates these files via
[`roxygen`](https://github.com/klutometis/roxygen) comments and
tags located in the R scripts that make up the package, i.e. all comments that
begin `#'`. For more information, look up
["object documentation in R"](http://r-pkgs.had.co.nz/man.html).

> **What does `build_readme = TRUE` do?** This builds the `README.md` file
from the `README.Rmd` file. The two files will contain the same text and images,
but for any R code snippets in the `.Rmd`, these are first run and their output
is captured and placed in the `.md` file. Remember that the `README.md` acts
as the homepage to the module on GitHub.

**File tree `module_path` (after build)**


```
## ├── DESCRIPTION
## ├── NAMESPACE
## ├── R
## │   ├── functions.R
## │   └── import.R
## ├── README.Rmd
## ├── README.html
## ├── README.md
## ├── examples
## │   └── example.R
## ├── inst
## │   ├── dockerfiles
## │   │   └── latest
## │   │       └── Dockerfile
## │   └── om.yml
## └── man
##     └── echo.Rd
```

## Creating/Building the Docker image

> **Dockerfile commands** Dockerfiles are series of instructions for
constructing a "containerised" machine called a Docker image (functionally it's
a bit like a virtual machine, but more lightweight). Each command begins
with a capitalised instruction followed by arguments. The most common
instruction would be "RUN", this executes command-line code within the Docker
image system. For example, `RUN echo "command!"` would pass "command!" to the
program "echo". All Dockerfiles begin with a FROM instruction, which pulls a
Docker image on which to build your own image. For example, many images are
built upon the Linux operating system Ubuntu in which case the first line of the
Dockerfile would be "FROM ubuntu:latest". This first line would then download
the latest Ubuntu Docker image and all subsequent "RUN" instructions would be
running in Ubuntu command-line. For far more detailed information on
Dockerfiles, see the
[Docker docs](https://docs.docker.com/engine/docker-overview/).

In our `om..echo` we have a dockerfiles folder that contains a Dockerfile
describing the Docker image for our `echo` program. Our "latest" Dockerfile
contains the instructions to pull the Docker image of the latest Ubuntu release,
to create a folder called "working_dir" and then set this new folder as the
"WORKDIR".

```
# Example host distro
FROM ubuntu:latest

# Install program using RUN lines

# outsider *requires* working_dir
RUN mkdir /working_dir
WORKDIR /working_dir
```

> **What's the WORKDIR?** The WORKDIR sets the working directory when a command
is passed to the Docker image. All `outsider` modules require this to be
"working_dir" as it allows `outsider` functions to know where to transfer files
to and from the container.

Using the Dockerfile that was created with the skeleton, we can now build our
Docker image using `module_build()`. (In practice, we'd combine steps 2 and 3
by calling `module_build` and setting all the build arguments to TRUE.) Docker
will build the image from the Dockerfile and store it along with all the other
images that are available on your machine, you do not need to worry about where
the Docker image is stored. With an existing image that is associated with the
module, when the module's code is called a new container is created from the
image and commands are passed to the container from R.


```r
module_build(flpth = module_path, tag = 'latest', build_image = TRUE,
             build_documents = FALSE, build_package = FALSE,
             build_readme = FALSE)
```

```
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running docker_build()
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Command:
## docker build -t dombennett/om_echo:latest /Library/Frameworks/R.framework/Versions/3.6/Resources/library/om..echo/dockerfiles/latest
## .....................................................................................................................
```

> **What is a tag?** A Docker tag is akin to the version number of the image.
By default, if no tag is provided,  Docker will use 'latest'.

> **Note**: The functions try to be helpful by printing to screen the exact
Docker commands of the tasks being performed. This is to give you a better idea
of what is happening and to allow you to recreate the operations via a terminal.

## Try the module.

With a Docker image and an installed R package we are now ready to try out the
module before we upload it online. The below command will look up the associated
image for module, create a new container, run the command and then shut-down the
container. (If no associated image can be found, it will attempt to pull one
from Docker Hub.)


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
# the repo always refers to the future github repo
echo <- module_import('echo', repo = 'dombennett/om..echo')
echo('hello world!')
```

## Upload to GitHub and Docker Hub

After we have played with the module and ensured it works as we would hope, we
can upload it to our GitHub and Docker accounts so that others may download it.

To upload to GitHub we must first create a version of the repository online.
Visit your GitHub account and then create a new repository by clicking
`Repositories > New`. Ensure to name the online version with the same name as
your module, i.e. `om..echo`. (No initial steps are required to upload to
Docker-Hub.)


```r
# based on the module details, the function will determine which code-sharing
# service to upload to
module_upload(flpth = module_path, code_sharing = TRUE, dockerhub = TRUE)
```

---

**Delete it all**

Don't want `om..echo` on your computer? Delete it so ...


```r
# to delete the Docker image and uninstall the R package
module_uninstall(repo = 'dombennett/om..echo')
```

```
## Removing 'om..echo'
```

```
## Removing package from '/Library/Frameworks/R.framework/Versions/3.6/Resources/library'
## (as 'lib' is unspecified)
```

```r
# to delete the repo folder
unlink(x = module_path, recursive = TRUE, force = TRUE)
```


---

## Next-up: [Intermediate](https://docs.ropensci.org/outsider.devtools/articles/intermediate.html)
---
title: "Intermediate: Module Building"
output: html_document
---



In [basics](https://docs.ropensci.org/outsider.devtools/articles/intermediate.html)
we introduced the essential elements that make up a module and then built a
super simple module.

On this page, we will build a module for a slightly more complex program while
introducing some more functionalities of `outsider`.

## `figlet`

For this walkthrough, we will be creating a slightly more sophisticated module
than `om..echo` by building a module for [`figlet`](http://www.figlet.org/),
a program that takes textual input and returns ASCII art.

The program itself only takes text input and returns text output. We, however,
will use `outsider` to create a module that can either take an input file or
text and then return ASCII art to either the console or to an output file.

Not all command-line programs follow the same standards (different input/output
file arguments, different help pages ...), so `outsider` comes with a series of
functions to provide options for developers for how best they want to
encapsulate the program in module form.

(For reference, this module does already exist in the wild:
[dombennett/om..figlet](https://github.com/DomBennett/om..figlet))

# Walkthrough

## Building

To get started, we will first need our module skeleton.


```r
library(outsider.devtools)
flpth <- module_skeleton(repo_user = 'dombennett', program_name = 'figlet',
                         docker_user = 'dombennett', flpth = tempdir(),
                         service = 'github')
# folder name where module is stored
print(basename(flpth))
```

```
## [1] "om..figlet"
```

Next we will need to develop the Dockerfile! `figlet` is easy enough to install
on Ubuntu Linux, we can just use the "apt-get" system to install it. (Here I'm
editing the file from within R to make the process reproducible, you can edit
the file directly with a text editor.)


```r
dockerfile_text <- "
# select Ubuntu Linux image
FROM ubuntu:latest
# Use the 'apt-get' system to install figlet
RUN apt-get update && apt-get install -y figlet
# set the working_dir/
RUN mkdir /working_dir
WORKDIR /working_dir
"
# write to latest Dockerfile
cat(dockerfile_text, file = file.path(flpth, 'inst', 'dockerfiles', 'latest',
                                      'Dockerfile'))
```

Now, we have our Dockerfile set-up and our module skeleton, we must develop the
R code. The skeleton code for our `figlet`-module is currently:


```r
#' @name figlet
#' @title figlet
#' @description Run figlet
#' @param ... Arguments
#' @example /examples/example.R
#' @export
figlet <- function(...) {
  # convert the ... into a argument list
  arglist <- arglist_get(...)
  # create an outsider object: describe the arguments and program
  otsdr <- outsider_init(pkgnm = 'om..figlet', cmd = 'figlet',
                         arglist = arglist)
  # run the command
  run(otsdr)
}
```

The above function would work fine. It would pass any arguments from the user
in R to the `figlet` program. But because `figlet` does not come with input/output arguments, we will need to hardcode these into our function and
use some if statements and bash to create this functionality. Because this will
require calling several commands via the terminal, we can create a shell script
that we can then pass through `outsider` to our container.


```r
function_text <- "
#' @name figlet
#' @title figlet
#' @description Run figlet with input and output files.
#' @param arglist Arguments for figlet
#' @param input_file Text input file
#' @param output_file ASCII art output file
#' @details If no input_file, will use arglist argument.
#' If no output_file, will print to console.
#' @example /examples/example.R
#' @export
figlet <- function(arglist = arglist_get(...), input_file = NULL,
                   output_file = NULL) {
  wd <- NULL
  # construct shell script from arglist
  arglist <- c('figlet', arglist)
  if (!is.null(input_file)) {
      # cat the input_file contents to figlet
      # (basename is used because on the container,
      #  filepaths cannot be used.)
      arglist <- c('cat', basename(input_file), '|', arglist)
  }
  if (!is.null(output_file)) {
      # write out the results of figlet to output_file
      arglist <- c(arglist, '>', basename(output_file))
      # determine where returned files should be sent
      wd <- dirpath_get(output_file)
  }
  # write arglist to temp file
  script <- file.path(tempdir(), 'script.sh')
  on.exit(file.remove(script))
  # ensure script is written in binary format
  script_cnntn <- file(script, 'wb')
  cmds <- paste(arglist, collapse = ' ')
  # debug print
  print(cmds)
  write(x = cmds, file = script_cnntn)
  close(script_cnntn)
  # initialise outsider container by specifying the command,
  # the arguments, the files to be sent, and the directory to where
  # returned files should be sent
  otsdr <- outsider_init(pkgnm = 'om..figlet', cmd = 'bash',
                         arglist = 'script.sh', wd = wd,
                         files_to_send = c(script, input_file))
  # run the command
  run(otsdr)
}
"
# write to R/functions.R
cat(function_text, file = file.path(flpth, 'R', 'functions.R'))
```

> **What does `arglist` do?** This function takes R objects and converts them
into character strings that then be passed to a command-line program. It can
take as many arguments as you want. For example, with `figlet` you may want to
convert the value of `a` to ASCII art. In this way, you would simply be able to
pass `a` to the `figlet` function rather than spell it out again.

> **What's with the `#'@`?** These are the roxygen tags. They determine the
documentation of the function. See
["Object documentation"](http://r-pkgs.had.co.nz/man.html) for more details.

Ok, let's build the module!


```r
module_build(flpth = flpth)
```

```
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::document() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```
## Updating om..figlet documentation
```

```
## First time using roxygen2. Upgrading automatically...
```

```
## Updating roxygen version in /private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..figlet/DESCRIPTION
```

```
## Writing NAMESPACE
```

```
## Loading om..figlet
```

```
## Writing NAMESPACE
## Writing figlet.Rd
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::install() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
##   
   checking for file ‘/private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..figlet/DESCRIPTION’ ...
  
✓  checking for file ‘/private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..figlet/DESCRIPTION’
## 
  
─  preparing ‘om..figlet’:
##    checking DESCRIPTION meta-information ...
  
✓  checking DESCRIPTION meta-information
## 
  
─  checking for LF line-endings in source and make files and shell scripts
## 
  
─  checking for empty or unneeded directories
## 
  
─  building ‘om..figlet_0.0.1.tar.gz’
## 
  
   
## 
Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
##   /var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T//RtmpquGDgt/om..figlet_0.0.1.tar.gz --install-tests 
## 
-
* installing to library ‘/Library/Frameworks/R.framework/Versions/3.6/Resources/library’
## 
|
* installing *source* package ‘om..figlet’ ...
## 
-
** using staged installation
## 
|
** R
## 
-
** inst
## 
|
** byte-compile and prepare package for lazy loading
## 
-

|
** help
## 
-
*** installing help indices
## 
|
** building package indices
## 
-
** testing if installed package can be loaded from temporary location
## 
|
** testing if installed package can be loaded from final location
## 
-
** testing if installed package keeps a record of temporary installation path
## 
|
* DONE (om..figlet)
## 
-

|

 
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running docker_build()
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Command:
## docker build -t dombennett/om_figlet:latest /Library/Frameworks/R.framework/Versions/3.6/Resources/library/om..figlet/dockerfiles/latest
## .....................................................................................................................
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::build_readme() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```
## Building om..figlet readme
```

Does it work?


```r
library(outsider)
figlet <- module_import('figlet', repo = 'dombennett/om..figlet')
# without files
figlet(arglist = 'hello!')
```

```
## [1] "figlet hello!"
```

```r
# with input file
input_file <- file.path(tempdir(), 'input_figlet.txt')
cat('This is from the input text file!', file = input_file)
figlet(arglist = '', input_file = input_file)
```

```
## [1] "cat input_figlet.txt | figlet "
```

```r
# with output file
output_file <- file.path(tempdir(), 'output_figlet.txt')
figlet(arglist = 'Into the output file', output_file = output_file)
```

```
## [1] "figlet Into the output file > output_figlet.txt"
```

```r
cat(readLines(con = output_file), sep = '\n')
```

```
##  ___       _          _   _                        _               _   
## |_ _|_ __ | |_ ___   | |_| |__   ___    ___  _   _| |_ _ __  _   _| |_ 
##  | || '_ \| __/ _ \  | __| '_ \ / _ \  / _ \| | | | __| '_ \| | | | __|
##  | || | | | || (_) | | |_| | | |  __/ | (_) | |_| | |_| |_) | |_| | |_ 
## |___|_| |_|\__\___/   \__|_| |_|\___|  \___/ \__,_|\__| .__/ \__,_|\__|
##                                                       |_|              
##   __ _ _      
##  / _(_) | ___ 
## | |_| | |/ _ \
## |  _| | |  __/
## |_| |_|_|\___|
## 
```

```r
# from input file to output file... with font block!
figlet(arglist = c('-f', 'block'), input_file = input_file,
       output_file = output_file)
```

```
## [1] "cat input_figlet.txt | figlet -f block > output_figlet.txt"
```

```r
cat(readLines(con = output_file), sep = '\n')
```

```
##                                                       
## _|_|_|_|_|  _|        _|                _|            
##     _|      _|_|_|          _|_|_|            _|_|_|  
##     _|      _|    _|  _|  _|_|          _|  _|_|      
##     _|      _|    _|  _|      _|_|      _|      _|_|  
##     _|      _|    _|  _|  _|_|_|        _|  _|_|_|    
##                                                       
##                                                       
##                                               
##     _|_|                                      
##   _|      _|  _|_|    _|_|    _|_|_|  _|_|    
## _|_|_|_|  _|_|      _|    _|  _|    _|    _|  
##   _|      _|        _|    _|  _|    _|    _|  
##   _|      _|          _|_|    _|    _|    _|  
##                                               
##                                               
##                                                                               
##   _|      _|                      _|                                  _|      
## _|_|_|_|  _|_|_|      _|_|            _|_|_|    _|_|_|    _|    _|  _|_|_|_|  
##   _|      _|    _|  _|_|_|_|      _|  _|    _|  _|    _|  _|    _|    _|      
##   _|      _|    _|  _|            _|  _|    _|  _|    _|  _|    _|    _|      
##     _|_|  _|    _|    _|_|_|      _|  _|    _|  _|_|_|      _|_|_|      _|_|  
##                                                 _|                            
##                                                 _|                            
##                                                                             
##   _|                            _|              _|_|  _|  _|            _|  
## _|_|_|_|    _|_|    _|    _|  _|_|_|_|        _|          _|    _|_|    _|  
##   _|      _|_|_|_|    _|_|      _|          _|_|_|_|  _|  _|  _|_|_|_|  _|  
##   _|      _|        _|    _|    _|            _|      _|  _|  _|            
##     _|_|    _|_|_|  _|    _|      _|_|        _|      _|  _|    _|_|_|  _|  
##                                                                             
## 
```

```r
# clean-up
file.remove(c(input_file, output_file))
```

```
## [1] TRUE TRUE
```

## Check, test and deploy

Ok! So we have a functioning module. But how do we know it's a functioning
module? To check our module structure, configuration and functioning we have a 
few helpful functions.


```r
# make sure the folder structure is correct
module_check(flpth = flpth)
```

```
## DESCRIPTION found ✓
## R folder with files found ✓
## inst found ✓
## inst/om.yml found ✓
## inst/dockerfiles found ✓
## inst/dockerfiles/latest with one Dockerfile found ✓
```

```r
# check are the names of the module components
module_identities(flpth = flpth)
```

```
## R package name: 'om..figlet'
## URL: 'https://github.com/dombennett/om..figlet'
## Docker images: 'dombennett/om_figlet:latest'
```

```r
# test that the module works
# module_test(flpth = flpth)
```

Running `module_test` we will find that there is an issue with the module.
`module_test` checks the functioning of the code by running the example; our
current example tries to call `-h` with `figlet` -- which is
an invalid argument. To correct this we can simply updated our example.


```r
example_text <- "
library(outsider)
figlet <- module_import('figlet', repo = 'dombennett/om..figlet')
figlet('hello!')
"
# write to examples/example.R
cat(example_text, file = file.path(flpth, 'examples', 'example.R'))
# re-build package
module_build(flpth = flpth, build_image = FALSE, verbose = FALSE)
```

```
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::document() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```
## Updating om..figlet documentation
```

```
## Writing NAMESPACE
```

```
## Loading om..figlet
```

```
## Writing NAMESPACE
## Writing figlet.Rd
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::install() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```
## 
## Attaching package: 'om..figlet'
```

```
## The following object is masked _by_ '.GlobalEnv':
## 
##     figlet
```

```
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
## Running devtools::build_readme() ...
## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

```
## Building om..figlet readme
```

Now re-running the test ...


```r
module_test(flpth = flpth)
```

```
## Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
##   /private/var/folders/x9/m8kwpxps2v93xk52zqzm5lkh0000gp/T/RtmpquGDgt/om..figlet --install-tests 
## * installing to library ‘/Library/Frameworks/R.framework/Versions/3.6/Resources/library’
## 
-
* installing *source* package ‘om..figlet’ ...
## 
|
** using staged installation
## 
-
** R
## 
|
** inst
## 
-
** byte-compile and prepare package for lazy loading
## 
|

-
** help
## 
|
*** installing help indices
## 
-
** building package indices
## 
|
** testing if installed package can be loaded from temporary location
## 
-
** testing if installed package can be loaded from final location
## 
|
** testing if installed package keeps a record of temporary installation path
## 
-
* DONE (om..figlet)
## 
|

-

 
```

```
## 
## Attaching package: 'om..figlet'
```

```
## The following object is masked _by_ '.GlobalEnv':
## 
##     figlet
```

```
## [1] "figlet hello!"
```

```
## Removing package from '/Library/Frameworks/R.framework/Versions/3.6/Resources/library'
## (as 'lib' is unspecified)
```

```
## Yikes! The module works! You are super-excellent!
```


Nice! The package passes. That means we can upload the package to GitHub and
Docker-Hub using the same functions as we saw in "basic".


```r
module_upload(flpth = flpth, code_sharing = TRUE, dockerhub = TRUE)
```

---

**Delete it all**


```r
module_uninstall(repo = 'om..figlet')
unlink(x = flpth, recursive = TRUE, force = TRUE)
```

---

## Next-up: [Advanced](https://docs.ropensci.org/outsider.devtools/articles/advanced.html)


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/git.R
\name{remote_git_exists}
\alias{remote_git_exists}
\title{Does a remote git repo exist}
\usage{
remote_git_exists(url)
}
\arguments{
\item{url}{URL to check}
}
\value{
Logical
}
\description{
Return TRUE if request returns a 200 status code.
}
\details{
Private repositories will not be discovered. Doesn't work for
bitbucket.
}
\seealso{
Other git: 
\code{\link{git_upload}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/git.R
\name{git_upload}
\alias{git_upload}
\title{Upload module to code-sharing service}
\usage{
git_upload(flpth, username, service = c("github", "gitlab", "bitbucket"))
}
\arguments{
\item{flpth}{File path to module.}

\item{username}{Username for code-sharing service.}

\item{service}{Code-sharing service}
}
\value{
Logical
}
\description{
Upload module to a git-based code-sharing service. Initiate
a git repo, add core module files, commit and push to remote.
}
\details{
Remote URL is determined to be: code sharing URL + username + R
package name.
Git must be configured on a user's system before this function can be run.
}
\examples{
# Git must be configured before function can be run.
# To configure, provide your name and email via the command-line, e.g.
# git config --global user.name "John Doe"
# git config --global user.email johndoe@example.com
# To check your username and email, try:
# git config --list

# NOT RUN
# # construct a simple module
# module_path <- module_skeleton(program_name = 'echo', flpth = getwd())
# module_check(flpth = module_path)
# module_identities(flpth = module_path)
# # after these steps all files are built, to upload to git use:
# git_upload(flpth = module_path, username = '[YOUR USERNAME]',
#            service = '[SERVICE]')
}
\seealso{
Other git: 
\code{\link{remote_git_exists}()}
}
\concept{git}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{string_replace}
\alias{string_replace}
\title{Replace patterns in a string}
\usage{
string_replace(string, patterns, values)
}
\arguments{
\item{string}{Text}

\item{patterns}{Patterns to replace with values}

\item{values}{Values to be put in place}
}
\value{
character
}
\description{
For a given character string, replace patterns with values.
}
\seealso{
Other utils: 
\code{\link{description_get}()},
\code{\link{examples_get}()},
\code{\link{file_create}()},
\code{\link{pkgdetails_get}()},
\code{\link{pkgnm_get}()},
\code{\link{tags_get}()},
\code{\link{templates_get}()},
\code{\link{yaml_get}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docker.R
\name{docker_push}
\alias{docker_push}
\title{Push a Docker image to Docker Hub}
\usage{
docker_push(username, img, tag, verbose)
}
\arguments{
\item{username}{Login username for Docker Hub.}

\item{img}{Image name}

\item{tag}{Tag, e.g. 'latest'}

\item{verbose}{Be verbose? T/F}
}
\value{
Logical
}
\description{
Push a Docker image to Docker Hub. Requires a user to have login
details for Docker Hub, \url{https://hub.docker.com/}.
Returns TRUE if no errors are raised.
}
\examples{
library(outsider.devtools)

# NOT RUN
# # create image from om..echo Docker file
# url <- paste0('https://raw.githubusercontent.com/DomBennett/om..hello.world/',
#               'master/inst/dockerfiles/latest/Dockerfile')
# docker_build(img = 'test_docker_push', tag = 'latest', verbose = TRUE,
#              url_or_path = url)
# # push to your Docker Hub
# docker_push(username = '[YOUR USERNAME]', img = 'test_docker_push',
#             tag = 'latest', verbose = TRUE)
}
\seealso{
Other docker: 
\code{\link{docker_build}()}
}
\concept{docker}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{description_get}
\alias{description_get}
\title{Return DESCRIPTION}
\usage{
description_get(flpth)
}
\arguments{
\item{flpth}{Path to module}
}
\value{
named list
}
\description{
Return contents of DESCRIPTION
}
\seealso{
Other utils: 
\code{\link{examples_get}()},
\code{\link{file_create}()},
\code{\link{pkgdetails_get}()},
\code{\link{pkgnm_get}()},
\code{\link{string_replace}()},
\code{\link{tags_get}()},
\code{\link{templates_get}()},
\code{\link{yaml_get}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{examples_test}
\alias{examples_test}
\title{Run each example of an outsider module}
\usage{
examples_test(flpth)
}
\arguments{
\item{flpth}{File path to location of module}
}
\value{
logical
}
\description{
Return TRUE if all of the outsider module functions successfully
run.
}
\seealso{
Other testing: 
\code{\link{test}()}
}
\concept{testing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{yaml_get}
\alias{yaml_get}
\title{Return om.yml}
\usage{
yaml_get(flpth)
}
\arguments{
\item{flpth}{Path to module}
}
\value{
list
}
\description{
Return contents of om.yml
}
\seealso{
Other utils: 
\code{\link{description_get}()},
\code{\link{examples_get}()},
\code{\link{file_create}()},
\code{\link{pkgdetails_get}()},
\code{\link{pkgnm_get}()},
\code{\link{string_replace}()},
\code{\link{tags_get}()},
\code{\link{templates_get}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/praise.R
\name{celebrate}
\alias{celebrate}
\title{Say nice things}
\usage{
celebrate()
}
\description{
Celebrate the passing of module tests.
}
\seealso{
Other praise: 
\code{\link{comfort}()}
}
\concept{praise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/praise.R
\name{comfort}
\alias{comfort}
\title{Say consoling things}
\usage{
comfort()
}
\description{
Commiserate the failing of module tests.
}
\details{
(Idea borrowed from \code{devtools}.)
}
\seealso{
Other praise: 
\code{\link{celebrate}()}
}
\concept{praise}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outsider.devtools.R
\docType{package}
\name{outsider.devtools}
\alias{outsider.devtools}
\title{outsider.devtools: Build 'outsider' Modules}
\description{
Tools, resources and information for making it easier to build your own
outsider modules.
}
\details{
For more information visit the outsider website
(\url{https://docs.ropensci.org/outsider.devtools/}).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{templates_get}
\alias{templates_get}
\title{Retrieve template files}
\usage{
templates_get()
}
\value{
character vector
}
\description{
Return template files for an outsider module.
}
\seealso{
Other utils: 
\code{\link{description_get}()},
\code{\link{examples_get}()},
\code{\link{file_create}()},
\code{\link{pkgdetails_get}()},
\code{\link{pkgnm_get}()},
\code{\link{string_replace}()},
\code{\link{tags_get}()},
\code{\link{yaml_get}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{module_identities}
\alias{module_identities}
\title{Return identities for a module}
\usage{
module_identities(flpth = getwd())
}
\arguments{
\item{flpth}{File path to location of module}
}
\value{
Logical
}
\description{
Returns a list of the identities (GitHub repo, Package name,
Docker images) for an outsider module. Works for modules in development.
Requires module to have a file path.
}
\examples{
library(outsider)

# build file structure for an example module
module_path <- module_skeleton(program_name = "goldenhind",
                               repo_user = "drake_on_github",
                               docker_user = "drake_on_docker",
                               full_name = 'Sir Francis Drake',
                               email = 'f.drake@goldenhind.gov.uk',
                               service = 'github',
                               flpth = tempdir())
# new path created 
(module_path)
# check the generated names and links
module_identities(flpth = module_path)
# check the files are in the right locations
module_check(flpth = module_path)
# deliberately break: delete a folder and check again
unlink(x = file.path(module_path, 'inst'), recursive = TRUE, force = TRUE)
module_check(flpth = module_path)

# clean-up
unlink(x = module_path, recursive = TRUE, force = TRUE)
}
\seealso{
Other build: 
\code{\link{module_build}()},
\code{\link{module_check}()},
\code{\link{module_skeleton}()},
\code{\link{module_test}()},
\code{\link{module_travis}()},
\code{\link{module_upload}()}
}
\concept{build}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{module_skeleton}
\alias{module_skeleton}
\title{Generate a skeleton for a module}
\usage{
module_skeleton(
  program_name,
  repo_user = NULL,
  docker_user = NULL,
  flpth = getwd(),
  module_name = NULL,
  cmd = program_name,
  full_name = NULL,
  email = NULL,
  service = c("github", "gitlab", "bitbucket"),
  overwrite = FALSE
)
}
\arguments{
\item{program_name}{Name of the command-line program.}

\item{repo_user}{Developer's username for code sharing service. If NULL, no
code sharing site information is added.}

\item{docker_user}{Developer's username for Docker. If NULL, no docker 
information is added.}

\item{flpth}{File path to location of where module will be created, default
current working directory.}

\item{module_name}{Name of the module, if NULL rendered as
"om..[program_name]"}

\item{cmd}{Command-line call for program, default [program_name]}

\item{full_name}{Your full name (for authorship)}

\item{email}{Your email (for authorship)}

\item{service}{Code-sharing site.}

\item{overwrite}{Automatically overwrite pre-existing files? If FALSE,
user is queried whether to overwrite for each pre-existing file.}
}
\value{
Character
}
\description{
Create all the base files and folders to kickstart the
development of a new outsider module. Returns file path to new module.
}
\details{
If \code{full_name} and \code{email} are provided, then new lines
are added to DESCRIPTION specifying the author and maintainer of the package.
}
\examples{
library(outsider)

# build file structure for an example module
module_path <- module_skeleton(program_name = "goldenhind",
                               repo_user = "drake_on_github",
                               docker_user = "drake_on_docker",
                               full_name = 'Sir Francis Drake',
                               email = 'f.drake@goldenhind.gov.uk',
                               service = 'github',
                               flpth = tempdir())
# new path created 
(module_path)
# check the generated names and links
module_identities(flpth = module_path)
# check the files are in the right locations
module_check(flpth = module_path)
# deliberately break: delete a folder and check again
unlink(x = file.path(module_path, 'inst'), recursive = TRUE, force = TRUE)
module_check(flpth = module_path)

# clean-up
unlink(x = module_path, recursive = TRUE, force = TRUE)
}
\seealso{
Other build: 
\code{\link{module_build}()},
\code{\link{module_check}()},
\code{\link{module_identities}()},
\code{\link{module_test}()},
\code{\link{module_travis}()},
\code{\link{module_upload}()}
}
\concept{build}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{module_travis}
\alias{module_travis}
\title{Generate Travis-CI file (GitHub only)}
\usage{
module_travis(flpth = getwd())
}
\arguments{
\item{flpth}{Directory in which to create .travis.yml}
}
\value{
Logical
}
\description{
Write .travis.yml to working directory.
}
\details{
Validated outsider modules must have a .travis.yml in their
repository. These \code{.travis.yml} are created with \link{module_skeleton}
but can also be generated using this function.
}
\examples{
library(outsider)

# a skeleton already comes with a .travis.yml
module_path <- module_skeleton(program_name = "goldenhind",
                               repo_user = "drake_on_github",
                               docker_user = "drake_on_docker",
                               service = 'github',
                               flpth = tempdir())
(file.exists(file.path(module_path, '.travis.yml')))
# but if it were deleted, needs updating or mistakenly edited ...
file.remove(file.path(module_path, '.travis.yml'))
# a new one can generated
module_travis(flpth = module_path)
}
\seealso{
Other build: 
\code{\link{module_build}()},
\code{\link{module_check}()},
\code{\link{module_identities}()},
\code{\link{module_skeleton}()},
\code{\link{module_test}()},
\code{\link{module_upload}()}
}
\concept{build}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{pkgdetails_get}
\alias{pkgdetails_get}
\title{Read the details of module's R package}
\usage{
pkgdetails_get(flpth)
}
\arguments{
\item{flpth}{Path to package}
}
\value{
list
}
\description{
Return a list of all package details based on a package's
DESCRIPTION file plus its \code{om.yml}.
}
\seealso{
Other utils: 
\code{\link{description_get}()},
\code{\link{examples_get}()},
\code{\link{file_create}()},
\code{\link{pkgnm_get}()},
\code{\link{string_replace}()},
\code{\link{tags_get}()},
\code{\link{templates_get}()},
\code{\link{yaml_get}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{module_test}
\alias{module_test}
\title{Test an outsider module}
\usage{
module_test(flpth = getwd(), verbose = FALSE, pull = FALSE)
}
\arguments{
\item{flpth}{File path to location of module}

\item{verbose}{Print docker and program info to console}

\item{pull}{Pull image from Docker Hub? T/F}
}
\value{
Logical
}
\description{
Ensure an outsider module builds, imports correctly and all
its functions successfully complete.
}
\details{
Success or fail, the module is uninstalled from the machine after
the test is run.
}
\examples{
library(outsider)
repo <- 'dombennett/om..hello.world'
# installs, checks functions, runs function examples, uninstalls
# module must already be uploaded to GitHub and Docker Hub
# NOT RUN
# module_test(repo = repo, verbose = TRUE)
}
\seealso{
Other build: 
\code{\link{module_build}()},
\code{\link{module_check}()},
\code{\link{module_identities}()},
\code{\link{module_skeleton}()},
\code{\link{module_travis}()},
\code{\link{module_upload}()}
}
\concept{build}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{examples_get}
\alias{examples_get}
\title{Determine examples}
\usage{
examples_get(flpth)
}
\arguments{
\item{flpth}{Path to module}
}
\value{
character
}
\description{
Determine example files for a module
}
\seealso{
Other utils: 
\code{\link{description_get}()},
\code{\link{file_create}()},
\code{\link{pkgdetails_get}()},
\code{\link{pkgnm_get}()},
\code{\link{string_replace}()},
\code{\link{tags_get}()},
\code{\link{templates_get}()},
\code{\link{yaml_get}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{module_check}
\alias{module_check}
\title{Check names and structure of a module}
\usage{
module_check(flpth = getwd())
}
\arguments{
\item{flpth}{File path to location of module}
}
\value{
Logical
}
\description{
Returns TRUE if all the names and structure of a minimal viable
outsider module are correct.
}
\examples{
library(outsider)

# NOT RUN
# # build a skeleton package
# module_path <- module_skeleton(program_name = 'echo', flpth = getwd())
# # check the file structure
# module_check(flpth = module_path)
# # look-up key identifying names: R package name, Docker image name
# module_identities(flpth = module_path)
# # build the R package and Docker image
# module_build(flpth = module_path, tag = 'latest')
# # test the module
# module_test(flpth = module_path)
# # clean-up
# unlink(x = module_path, recursive = TRUE)
}
\seealso{
Other build: 
\code{\link{module_build}()},
\code{\link{module_identities}()},
\code{\link{module_skeleton}()},
\code{\link{module_test}()},
\code{\link{module_travis}()},
\code{\link{module_upload}()}
}
\concept{build}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_from_base.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{arglist_get}
\alias{arglist_parse}
\alias{dirpath_get}
\alias{filestosend_get}
\alias{wd_get}
\alias{docker_img_ls}
\alias{image_install}
\alias{is_docker_available}
\alias{outsider_init}
\alias{run}
\alias{is_installed}
\alias{meta_get}
\alias{modules_list}
\alias{pkg_install}
\alias{uninstall}
\alias{cat_line}
\alias{char}
\alias{func}
\alias{log_set}
\alias{stat}
\alias{server_connect}
\alias{server_disconnect}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{outsider.base}{\code{\link[outsider.base]{arglist_get}}, \code{\link[outsider.base]{arglist_parse}}, \code{\link[outsider.base]{cat_line}}, \code{\link[outsider.base]{char}}, \code{\link[outsider.base]{dirpath_get}}, \code{\link[outsider.base]{docker_img_ls}}, \code{\link[outsider.base]{filestosend_get}}, \code{\link[outsider.base]{func}}, \code{\link[outsider.base]{image_install}}, \code{\link[outsider.base]{is_docker_available}}, \code{\link[outsider.base]{is_installed}}, \code{\link[outsider.base]{log_set}}, \code{\link[outsider.base]{meta_get}}, \code{\link[outsider.base]{modules_list}}, \code{\link[outsider.base]{outsider_init}}, \code{\link[outsider.base]{pkg_install}}, \code{\link[outsider.base]{run}}, \code{\link[outsider.base]{server_connect}}, \code{\link[outsider.base]{server_disconnect}}, \code{\link[outsider.base]{stat}}, \code{\link[outsider.base]{uninstall}}, \code{\link[outsider.base]{wd_get}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{tags_get}
\alias{tags_get}
\title{Determine tags}
\usage{
tags_get(flpth)
}
\arguments{
\item{flpth}{Path to module}
}
\value{
character
}
\description{
Determine module Docker tags
}
\seealso{
Other utils: 
\code{\link{description_get}()},
\code{\link{examples_get}()},
\code{\link{file_create}()},
\code{\link{pkgdetails_get}()},
\code{\link{pkgnm_get}()},
\code{\link{string_replace}()},
\code{\link{templates_get}()},
\code{\link{yaml_get}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{module_upload}
\alias{module_upload}
\title{Upload a module to code sharing site and DockerHub}
\usage{
module_upload(
  flpth = getwd(),
  code_sharing = TRUE,
  dockerhub = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{flpth}{File path to location of module}

\item{code_sharing}{Upload to code sharing service?}

\item{dockerhub}{Upload to DockerHub?}

\item{verbose}{Print docker and program info to console}
}
\value{
Logical
}
\description{
Look up usernames and other information contained in
"om.yml" to upload module to a code sharing site (github, gitlab or
bitbucket) and/or DockerHub.
}
\details{
This function runs \code{\link{git_upload}} and
\code{\link{docker_push}}.
}
\examples{
library(outsider)

# NOT RUN
# # build a simple module
# module_path <- module_skeleton(program_name = 'echo', flpth = getwd(),
#                                repo_user = '[YOUR USERNAME]',
#                                docker_user = '[YOUR USERNAME]',
#                                service = 'github')
# module_build(module_path)
# # upload to your GitHub and your DockerHub
# module_upload(flpth = module_path)
}
\seealso{
Other build: 
\code{\link{module_build}()},
\code{\link{module_check}()},
\code{\link{module_identities}()},
\code{\link{module_skeleton}()},
\code{\link{module_test}()},
\code{\link{module_travis}()}
}
\concept{build}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{module_build}
\alias{module_build}
\title{Build a module}
\usage{
module_build(
  flpth = getwd(),
  tag = "latest",
  build_documents = TRUE,
  build_package = TRUE,
  build_image = TRUE,
  build_readme = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{flpth}{File path to location of module.}

\item{tag}{Docker tag, e.g. latest.}

\item{build_documents}{Build R documentation? T/F}

\item{build_package}{Build R package? T/F}

\item{build_image}{Build Docker image? T/F}

\item{build_readme}{Build README.md? T/F}

\item{verbose}{Be verbose? T/F}
}
\value{
Logical
}
\description{
Do
}
\examples{
library(outsider)

# NOT RUN
# # build a skeleton package
# module_path <- module_skeleton(program_name = 'echo', flpth = getwd())
# # check the file structure
# module_check(flpth = module_path)
# # look-up key identifying names: R package name, Docker image name
# module_identities(flpth = module_path)
# # build the R package and Docker image
# module_build(flpth = module_path, tag = 'latest')
# # test the module
# module_test(flpth = module_path)
# # clean-up
# unlink(x = module_path, recursive = TRUE)
}
\seealso{
Other build: 
\code{\link{module_check}()},
\code{\link{module_identities}()},
\code{\link{module_skeleton}()},
\code{\link{module_test}()},
\code{\link{module_travis}()},
\code{\link{module_upload}()}
}
\concept{build}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{pkgnm_get}
\alias{pkgnm_get}
\title{Determine R package name}
\usage{
pkgnm_get(flpth)
}
\arguments{
\item{flpth}{Path to module}
}
\value{
character(1)
}
\description{
Determine R package name from file path
}
\seealso{
Other utils: 
\code{\link{description_get}()},
\code{\link{examples_get}()},
\code{\link{file_create}()},
\code{\link{pkgdetails_get}()},
\code{\link{string_replace}()},
\code{\link{tags_get}()},
\code{\link{templates_get}()},
\code{\link{yaml_get}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docker.R
\name{docker_build}
\alias{docker_build}
\title{Build a Docker image}
\usage{
docker_build(img, tag, url_or_path, verbose)
}
\arguments{
\item{img}{Image name}

\item{tag}{Tag, e.g. 'latest'}

\item{url_or_path}{URL or file path to Dockerfile}

\item{verbose}{Be verbose? T/F}
}
\value{
Logical
}
\description{
Build a Docker image through a system call. Returns TRUE if
no errors are raised.
}
\examples{
library(outsider.devtools)

# # NOT RUN
# 
# # simplest possible Dockerfile
# df_text <- "
# FROM ubuntu:latest
# "
# 
# # create dir to host Dockerfile
# flpth <- file.path(tempdir(), 'test_docker_build')
# if(!dir.exists(flpth)) {
#   dir.create(flpth)
# }
# 
# # write to file
# write(x = df_text, file = file.path(flpth, 'Dockerfile'))
# 
# # run docker_build from flpth
# docker_build(img = 'test_docker_build', tag = 'latest', verbose = TRUE,
#              url_or_path = flpth)
}
\seealso{
Other docker: 
\code{\link{docker_push}()}
}
\concept{docker}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{file_create}
\alias{file_create}
\title{Create file}
\usage{
file_create(x, flpth, overwrite)
}
\arguments{
\item{x}{Text for writing to file}

\item{flpth}{File path to be created}

\item{overwrite}{Overwrite pre-existing file? Logical.}
}
\description{
Write x to a filepath. Forces creation of directories.
}
\seealso{
Other utils: 
\code{\link{description_get}()},
\code{\link{examples_get}()},
\code{\link{pkgdetails_get}()},
\code{\link{pkgnm_get}()},
\code{\link{string_replace}()},
\code{\link{tags_get}()},
\code{\link{templates_get}()},
\code{\link{yaml_get}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{test}
\alias{test}
\title{Test a module}
\usage{
test(flpth, pull = FALSE)
}
\arguments{
\item{flpth}{File path to location of module}

\item{pull}{Pull Docker image from Docker Hub? T/F}
}
\value{
logical
}
\description{
Test an outsider module by making sure it installs,
imports and its examples run correctly. Raises an error if a test fails.
}
\seealso{
Other testing: 
\code{\link{examples_test}()}
}
\concept{testing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docker.R
\name{docker_login}
\alias{docker_login}
\title{Login to Docker}
\usage{
docker_login(username)
}
\arguments{
\item{username}{Username for Docker-Hub.}
}
\value{
Logical
}
\description{
Login to docker using username. User is prompted to provide
password.
Returns TRUE if no errors are raised.
}
\seealso{
Other docker-private: 
\code{\link{docker_cmd}()}
}
\concept{docker-private}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docker.R
\name{docker_cmd}
\alias{docker_cmd}
\title{Send commands to Docker}
\usage{
docker_cmd(args, std_out = TRUE, std_err = TRUE)
}
\arguments{
\item{args}{Vector of arguments for "docker" command. E.g. '--help'.}

\item{std_out}{Logical or file path}

\item{std_err}{Logical or file path}
}
\value{
Logical
}
\description{
Safely send commands to Docker. All commands are also echo'd.
}
\seealso{
Other docker-private: 
\code{\link{docker_login}()}
}
\concept{docker-private}
