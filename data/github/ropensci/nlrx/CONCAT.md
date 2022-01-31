
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nlrx <img src="man/figures/logo.png" align="right" width="150" />

<!-- old badges: [![Build Status](https://travis-ci.org/ropensci/nlrx.svg?branch=master)](https://travis-ci.org/ropensci/nlrx)[
![Build status](https://ci.appveyor.com/api/projects/status/swsstjxxjnkyuoh9/branch/master?svg=true)](https://ci.appveyor.com/project/marcosci/nlrx/branch/master) -->

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/nlrx/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/nlrx/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/nlrx/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/nlrx)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![CRAN
status](https://www.r-pkg.org/badges/version/nlrx)](https://cran.r-project.org/package=nlrx)
[![](http://cranlogs.r-pkg.org/badges/grand-total/nlrx)](https://cran.r-project.org/package=nlrx)
[![ropensci](https://badges.ropensci.org/262_status.svg)](https://github.com/ropensci/software-review/issues/262)
[![DOI:10.1111/2041-210X.13286](https://zenodo.org/badge/DOI/10.1111/2041-210X.13286.svg)](https://doi.org/10.1111/2041-210X.13286)
<!-- badges: end -->

The nlrx package provides tools to setup and execute NetLogo simulations
from R. NetLogo is a free, open-source and cross-platform modelling
environment for simulating natural and social phenomena. NetLogo
focusses on implementation of agent-based and spatially explicit
simulation models, although system dynamics models are supported as
well. NetLogo is developed and maintained at the Center for Connected
Learning and Computer-Based Modeling, Northwestern University, Evanston,
IL. More details on NetLogo itself are available online: [NetLogo online
documentation](https://ccl.northwestern.edu/netlogo/docs/)

NetLogo comes with the built-in experiment tool [Behavior
Space](https://ccl.northwestern.edu/netlogo/docs/behaviorspace.html)
that allows to setup and execute model simulations with different
settings and parameter variations and to collect model output. This
experiment tool can be executed via command line in combination with an
XML file that contains the experiment specifications, such as runtime,
variables, output measurements, stop conditions, and more. One
limitation of Behavior Space is, that it only supports full-factorial
parameter designs, which may not be appropriate for complex model
analyses. Furthermore, Behavior Space experiment specifications are
stored within the NetLogo file and are not easily accessible from R.
However, in many cases it is useful to store such specifications along
with the model output and analyses results in order to enable fully
reproducible model analyses.

The nlrx package utilizes the commandline functionality of Behavior
Space to execute NetLogo simulations directly from R. Instead of
defining experiments within NetLogo Behavior Space, experiments are
defined in R using the class objects of the nlrx package. These class
objects hold all the information that is needed to run these experiments
remotely from R, such as path to NetLogo installation folder, path to
the model file and the experiment specifications itself. nlrx provides
useful helper functions to generate parameter input matrices from
parameter range definitions that cover a wide range of parameter
exploration approaches. By storing all relevant information on
simulation experiments, including the output of the model simulations in
one class object, experiments can be easily stored and shared.

In summary, the nlrx package uses a similar structure as NetLogos
Behavior Space but offers more flexibility and additional tools for
running reproducible complex model analyses directly from R.

## Publication

Further information on the package functionality and detailed code
examples can be found in our accompanying publication: Salecker J,
Sciaini M, Meyer KM, Wiegand K. The nlrx r package: A next-generation
framework for reproducible NetLogo model analyses. Methods Ecol Evol.
2019;2041-210X. <https://doi.org/10.1111/2041-210X.13286>.

Get citation information for `nlrx` in R doing `citation(package =
'nlrx')`.

## Prerequirements

### NetLogo

In order to use the nlrx package, [NetLogo](http://netlogoweb.org/)
(\>=5.3.1) needs to be installed on the system that is used to execute
model simulations (local/remote). For remote execution, NetLogo needs to
be installed on remote machines as well. The nlrx package provides a
utility function (`download_netlogo()`) that can be used to download and
unzip (only unix systems) a specified NetLogo version to a local folder.
For windows machines, the downloaded file needs to be executed in order
to install NetLogo on the local system. If you are running MacOS, please
use the Linux tar.gz version of NetLogo (either from the NetLogo
Homepage or by using the `download_netlogo()` function). The dmg version
from the NetLogo homepage is not compatible with nlrx.

### Java

Because NetLogo is executed in a Java virtual machine, Java needs to be
installed on the local/remote system as well. We recommend the [Oracle
Java SE Development
Kit 8](https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html)
or the [openjdk](https://github.com/ojdkbuild/ojdkbuild). While the nlrx
package might work without setting the Java system path explicitly, we
recommend to make sure that JAVA\_HOME points to the correct Java
installation of the system.

## Installation

You can install the released version of nlrx from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("nlrx")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/nlrx")
```

## Get started

General information that is needed to run NetLogo simulations remotely,
such as path to the NetLogo installation folder is stored within a `nl`
class object. Nested within this `nl` class are the classes `experiment`
and `simdesign`. The `experiment` class stores all experiment
specifications. After attaching a valid experiment, a `simdesign` class
object can be attached to the `nl` class object, by using one of the
simdesign helper functions. These helper functions create different
parameter input matrices from the experiment variable definitions that
can then be executed by the `run_nl_one()` and `run_nl_all()` functions.
The nested design allows to store everything related to the experiment
within one R object. Additionally, different simdesign helper functions
can be applied to the same `nl` object in order to repeat the same
experiment with different parameter exploration methods (simdesigns).

### Step by step application example

The “Wolf Sheep Predation” model from the NetLogo models library is used
to present a basic example on how to setup and run NetLogo model
simulations from R.

#### Step 1: Create a nl object:

The nl object holds all information on the NetLogo version, a path to
the NetLogo directory with the defined version, a path to the model
file, and the desired memory for the java virtual machine. Depending on
the operation system, paths to NetLogo and the model need to be
adjusted.

``` r
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

#### Step 2: Attach an experiment

The experiment object is organized in a similar fashion as NetLogo
Behavior Space experiments. It contains all information that is needed
to generate a simulation parameter matrix and to execute the NetLogo
simulations. Details on the specific slots of the experiment class can
be found in the package documentation (`?experiment`) and the “Advanced
configuration” vignette.

``` r
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                             'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))
```

#### Step 3: Attach a simulation design

While the experiment defines the variables and specifications of the
model, the simulation design creates a parameter input table based on
these model specifications and the chosen simulation design method. nlrx
provides a bunch of different simulation designs, such as
full-factorial, latin-hypercube, sobol, morris and eFast (see “Simdesign
Examples” vignette for more information on simdesigns). All simdesign
helper functions need a properly defined nl object with a valid
experiment design. Each simdesign helper also allows to define a number
of random seeds that are randomly generated and can be used to execute
repeated simulations of the same parameter matrix with different
random-seeds (see “Advanced configuration” vignette for more information
on random-seed and repetition management). A simulation design is
attached to a nl object by using one of the simdesign helper functions:

``` r
nl@simdesign <- simdesign_lhs(nl=nl,
                               samples=100,
                               nseeds=3,
                               precision=3)
```

#### Step 4: Run simulations

All information that is needed to run the simulations is now stored
within the nl object. The `run_nl_one()` function allows to run one
specific simulation from the siminput parameter table. The
`run_nl_all()` function runs a loop over all simseeds and rows of the
parameter input table siminput. The loops are constructed in a way that
allows easy parallelisation, either locally or on remote HPC machines
(see “Advanced configuration” vignette for more information on
parallelisation). Before running your simulations you might want to
check your current nl object setup. `eval_variables_constants(nl)`
evaluates if the defined variables and constants are correctly defined
and are consistent with the attached model. `print(nl)` prints a
complete summary of the provided nl object including checkmarks that
might help to indicate potential problems.

``` r
# Evaluate nl object:
eval_variables_constants(nl)
print(nl)

# Run all simulations (loop over all siminputrows and simseeds)
results <- run_nl_all(nl)
```

#### Step 5: Attach results to nl and run analysis

nlrx provides method specific analysis functions for each simulation
design. Depending on the chosen design, the function reports a tibble
with aggregated results or sensitivity indices. In order to run the
analyze\_nl function, the simulation output has to be attached to the nl
object first. The simdesign class within the nl object provides a slot
for attaching output results (simoutput). An output results tibble can
be attached to this slot by using the simdesign setter function
`setsim(nl, "simoutput")`. After attaching the simulation results, these
can also be written to the defined outpath of the experiment object.

``` r
# Attach results to nl object:
setsim(nl, "simoutput") <- results

# Write output to outpath of experiment within nl
write_simoutput(nl)

# Do further analysis:
analyze_nl(nl)
```

## Meta

  - Please [report any issues or
    bugs](https://github.com/ropensci/nlrx/issues/new).
  - License: GPL3
  - Get citation information for `nlrx` in R doing `citation(package =
    'nlrx')`
  - We are very open to contributions - if you are interested check
    [Contributing](https://github.com/ropensci/nlrx/blob/master/CONTRIBUTING.md).
      - Please note that this project is released with a [Contributor
        Code of
        Conduct](https://github.com/ropensci/nlrx/blob/master/CODE_OF_CONDUCT.md).
        By participating in this project you agree to abide by its
        terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

# nlrx 0.4.3

## Functionality

* added test_nlrx() function to check functionality of the package
* added support to download_netLogo() for NetLogo version 6.2.0

## Bugfixes

* changed timestamp for nldoc function from lubridate to base date
* fixed an error in the function parser of the nldoc procedure
* fixed bug in calculation of number of computed runs in print_nl()


# nlrx 0.4.2

## Functionality

* added option to run_nl_one that allows to store results as rds files
* added eval_simoutput option to check for missing combinations of siminputrow and random-seeds
* added support for progressr progress bars for run_nl_all function (details see further notes vignette) and removed the silent parameter of the run_nl_all function

## Bugfixes
* hotfix for another dependency on external files in nldoc roxygen examples
* small bugfix in analyze_morris: A warning is now thrown if NA are present in the simulation data
* bugfix in random seed generator
* bugfix for sobol simulation design when sobolorder is higher than the available number of variables
* analyze_nl now prints a warning if missing combinations were detected in the simulation output
* updated testdata
* user rights for temporary sh scripts are now set correctly

# nlrx 0.4.1
* fixed dependency on external file source in nldoc automated tests
* these files are now included in the package
* added link to documentation website in description

# nlrx 0.4.0

* Added new simdesigns simdesign_ABCmcmc_Marjoram, simdesign_ABCmcmc_Marjoram_original and simdesign_ABCmcmc_Wegmann to perform approximate bayesian computation
* Added print function for nl objects
* Added dependencies: crayon, EasyABC
* Added pandoc to system requirements
* Added additional pandoc_available() check for nldoc function
* Added support to download_netLogo() for NetLogo version 6.1.1
* Added new vignette showing an example for approximate bayesian computation with nlrx
* Updated "Sensitivity Analyses with nlrx"" vignette
* Updated "Advanced Configuration" vignette
* Updated package tests


# nlrx 0.3.0

* Added support for self-defined evaluation functions to optimization functions simdesign_GenAlg and simdesign_GenSA.
* Added support to simdesign_simple() for models without any GUI parameters
* Added support to download_netLogo() for NetLogo version 6.1.0
* Added new vignette showing an example for Sensitivity Analysis with nlrx.
* Added new vignette showing an example for Optimization with nlrx.
* Updated "Advanced Configuration" vignette.
* Updated citation information of the package.
* Hotfix for unnest_simoutput(). In the previous package version, under some circumstances an error occured due to NA data.
* Corrected spelling errors in some vignettes and documentation files.

# nlrx 0.2.0

* nl_to_raster() hotfix

# nlrx 0.1.0

* First release to CRAN.
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
# CONTRIBUTING #

### Please contribute!

We love collaboration.
In fact, this was one of the main ideas to start the work on this package.
Of course, you are also welcome to suggest general improvements to the package structure or to whatsoever.
We appreciate any contribution and collaboration.

### Bugs?

* Submit an issue on the Issues page [here](https://github.com/ropensci/nlrx/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/ropensci/nlrx.git`
* Make sure to track progress upstream (i.e., on our version of `nlrx` at `ropensci/nlrx`) by doing `git remote add upstream https://github.com/ropensci/nlrx.git`.
Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new branch)
* If you alter package functionality at all (e.g., the code itself, not just documentation) please do write some tests to cover the new functionality
* Push up to your account
* Submit a pull request to home base at `ropensci/nlrx`

### Questions? Get in touch: [nlrx@mailbox.org](mailto:nlrx@mailbox.org)

### Thanks for contributing!
# Resubmission

This is a resubmission. In this version I have:

* Updated the broken URL for the lifecycle badge in the README.md 


# Original (Re-) Submission

nlrx version 0.4.3

## CRAN check error fixes
* updated maintainer and maintainer email address (package was archived because of an auto-response from the previous email address)

## Changes in version 0.4.3

#### Functionality
* added test_nlrx() function to check functionality of the package
* added support to download_netLogo() for NetLogo version 6.2.0

#### Bugfixes
* changed timestamp for nldoc function from lubridate to base date
* fixed an error in the function parser of the nldoc procedure
* fixed bug in calculation of number of computed runs in print_nl()

## Test environments
* win-builder (release and devel)
* Windows 10, R 4.1.0
* macOS 11.3 (Big Sur), R 4.1.0
* macOS 10.14 (Mojave), R 4.1.0
* macOS 10.15 (Catalina), R 4.1.0 (github actions)
* Windows Server 2019 x64, R 4.1.0 (github actions)
* Ubuntu 20.04, R 4.1.0 (github actions)

## R CMD check results

0 errors | 0 warnings | 0 note

## Reverse dependencies

There are currently no reverse dependencies.
---
output:
  github_document:
    html_preview: false
editor_options: 
  chunk_output_type: console
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# nlrx <img src="man/figures/logo.png" align="right" width="150" />

<!-- old badges: [![Build Status](https://travis-ci.org/ropensci/nlrx.svg?branch=master)](https://travis-ci.org/ropensci/nlrx)[
![Build status](https://ci.appveyor.com/api/projects/status/swsstjxxjnkyuoh9/branch/master?svg=true)](https://ci.appveyor.com/project/marcosci/nlrx/branch/master) -->

<!-- badges: start -->
[![R build status](https://github.com/ropensci/nlrx/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/nlrx/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/nlrx/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/nlrx)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN status](https://www.r-pkg.org/badges/version/nlrx)](https://cran.r-project.org/package=nlrx)
[![](http://cranlogs.r-pkg.org/badges/grand-total/nlrx)](https://cran.r-project.org/package=nlrx)
[![ropensci](https://badges.ropensci.org/262_status.svg)](https://github.com/ropensci/software-review/issues/262)
[![DOI:10.1111/2041-210X.13286](https://zenodo.org/badge/DOI/10.1111/2041-210X.13286.svg)](https://doi.org/10.1111/2041-210X.13286)
<!-- badges: end -->

The nlrx package provides tools to setup and execute NetLogo simulations from R.
NetLogo is a free, open-source and cross-platform modelling environment for simulating natural and social phenomena.
NetLogo focusses on implementation of agent-based and spatially explicit simulation models, although system dynamics models are supported as well.
NetLogo is developed and maintained at the Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.
More details on NetLogo itself are available online: [NetLogo online documentation](https://ccl.northwestern.edu/netlogo/docs/)

NetLogo comes with the built-in experiment tool [Behavior Space](https://ccl.northwestern.edu/netlogo/docs/behaviorspace.html) that allows to setup and execute model simulations with different settings and parameter variations and to collect model output. This experiment tool can be executed via command line in combination with an XML file that contains the experiment specifications, such as runtime, variables, output measurements, stop conditions, and more. One limitation of Behavior Space is, that it only supports full-factorial parameter designs, which may not be appropriate for complex model analyses. Furthermore, Behavior Space experiment specifications are stored within the NetLogo file and are not easily accessible from R. However, in many cases it is useful to store such specifications along with the model output and analyses results in order to enable fully reproducible model analyses.

The nlrx package utilizes the commandline functionality of Behavior Space to execute NetLogo simulations directly from R.
Instead of defining experiments within NetLogo Behavior Space, experiments are defined in R using the class objects of the nlrx package.
These class objects hold all the information that is needed to run these experiments remotely from R, such as path to NetLogo installation folder, path to the model file and the experiment specifications itself. nlrx provides useful helper functions to generate parameter input matrices from parameter range definitions that cover a wide range of parameter exploration approaches. By storing all relevant information on simulation experiments, including the output of the model simulations in one class object, experiments can be easily stored and shared.

In summary, the nlrx package uses a similar structure as NetLogos Behavior Space but offers more flexibility and additional tools for running reproducible complex model analyses directly from R.

## Publication
Further information on the package functionality and detailed code examples can be found in our accompanying publication:
Salecker J, Sciaini M, Meyer KM, Wiegand K. The nlrx r package: A next-generation framework for reproducible NetLogo model analyses. Methods Ecol Evol. 2019;2041-210X. [https://doi.org/10.1111/2041-210X.13286](https://doi.org/10.1111/2041-210X.13286).

Get citation information for `nlrx` in R doing `citation(package = 'nlrx')`.

## Prerequirements

### NetLogo
In order to use the nlrx package, [NetLogo](http://netlogoweb.org/) (>=5.3.1) needs to be installed on the system that is used to execute model simulations (local/remote). For remote execution, NetLogo needs to be installed on remote machines as well. The nlrx package provides a utility function (`download_netlogo()`) that can be used to download and unzip (only unix systems) a specified NetLogo version to a local folder. For windows machines, the downloaded file needs to be executed in order to install NetLogo on the local system. If you are running MacOS, please use the Linux tar.gz version of NetLogo (either from the NetLogo Homepage or by using the `download_netlogo()` function). The dmg version from the NetLogo homepage is not compatible with nlrx.

### Java
Because NetLogo is executed in a Java virtual machine, Java needs to be installed on the local/remote system as well. We recommend the [Oracle Java SE Development Kit 8](https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) or the [openjdk](https://github.com/ojdkbuild/ojdkbuild).
While the nlrx package might work without setting the Java system path explicitly, we recommend to make sure that JAVA_HOME points to the correct Java installation of the system.

## Installation

You can install the released version of nlrx from
[CRAN](https://CRAN.R-project.org) with:
``` r
install.packages("nlrx")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/nlrx")
```

## Get started

General information that is needed to run NetLogo simulations remotely, such as path to the NetLogo installation folder is stored within a `nl` class object.
Nested within this `nl` class are the classes `experiment` and `simdesign`. The `experiment` class stores all experiment specifications. After attaching a valid experiment, a `simdesign` class object can be attached to the `nl` class object, by using one of the simdesign helper functions. These helper functions create different parameter input matrices from the experiment variable definitions that can then be executed by the `run_nl_one()` and `run_nl_all()` functions. The nested design allows to store everything related to the experiment within one R object. Additionally, different simdesign helper functions can be applied to the same `nl` object in order to repeat the same experiment with different parameter exploration methods (simdesigns).

### Step by step application example

The "Wolf Sheep Predation" model from the NetLogo models library is used to present a basic example on how to setup and run NetLogo model simulations from R.

#### Step 1: Create a nl object:

The nl object holds all information on the NetLogo version, a path to the NetLogo directory with the defined version, a path to the model file, and the desired memory for the java virtual machine. Depending on the operation system, paths to NetLogo and the model need to be adjusted.

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

#### Step 2: Attach an experiment

The experiment object is organized in a similar fashion as NetLogo Behavior Space experiments.
It contains all information that is needed to generate a simulation parameter matrix and to execute the NetLogo simulations.
Details on the specific slots of the experiment class can be found in the package documentation (`?experiment`) and the "Advanced configuration" vignette.

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                             'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))
```

#### Step 3: Attach a simulation design

While the experiment defines the variables and specifications of the model, the simulation design creates a parameter input table based on these model specifications and the chosen simulation design method.
nlrx provides a bunch of different simulation designs, such as full-factorial, latin-hypercube, sobol, morris and eFast (see "Simdesign Examples" vignette for more information on simdesigns).
All simdesign helper functions need a properly defined nl object with a valid experiment design. Each simdesign helper also allows to define a number of random seeds that are randomly generated and can be used to execute repeated simulations of the same parameter matrix with different random-seeds (see "Advanced configuration" vignette for more information on random-seed and repetition management).
A simulation design is attached to a nl object by using one of the simdesign helper functions:

```{r eval=FALSE}
nl@simdesign <- simdesign_lhs(nl=nl,
                               samples=100,
                               nseeds=3,
                               precision=3)
```

#### Step 4: Run simulations

All information that is needed to run the simulations is now stored within the nl object.
The `run_nl_one()` function allows to run one specific simulation from the siminput parameter table.
The `run_nl_all()` function runs a loop over all simseeds and rows of the parameter input table siminput.
The loops are constructed in a way that allows easy parallelisation, either locally or on remote HPC machines (see "Advanced configuration" vignette for more information on parallelisation).
Before running your simulations you might want to check your current nl object setup. `eval_variables_constants(nl)` evaluates if the defined variables and constants are correctly defined and are consistent with the attached model. `print(nl)` prints a complete summary of the provided nl object including checkmarks that might help to indicate potential problems.


```{r eval=FALSE}
# Evaluate nl object:
eval_variables_constants(nl)
print(nl)

# Run all simulations (loop over all siminputrows and simseeds)
results <- run_nl_all(nl)
```

#### Step 5: Attach results to nl and run analysis

nlrx provides method specific analysis functions for each simulation design.
Depending on the chosen design, the function reports a tibble with aggregated results or sensitivity indices.
In order to run the analyze_nl function, the simulation output has to be attached to the nl object first.
The simdesign class within the nl object provides a slot for attaching output results (simoutput). An output results tibble can be attached to this slot by using the simdesign setter function `setsim(nl, "simoutput")`.
After attaching the simulation results, these can also be written to the defined outpath of the experiment object.

```{r eval=FALSE}
# Attach results to nl object:
setsim(nl, "simoutput") <- results

# Write output to outpath of experiment within nl
write_simoutput(nl)

# Do further analysis:
analyze_nl(nl)
``` 

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/nlrx/issues/new).
* License: GPL3
* Get citation information for `nlrx` in R doing `citation(package = 'nlrx')`
* We are very open to contributions - if you are interested check [Contributing](https://github.com/ropensci/nlrx/blob/master/CONTRIBUTING.md).
    * Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/nlrx/blob/master/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Get Started"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Get Started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# nlrx

The nlrx package provides tools to setup and execute NetLogo simulations from R.
NetLogo is a free, open-source and cross-platform modelling environment for simulating natural and social phenomena.
NetLogo focusses on implementation of agent-based and spatially explicit simulation models, although system dynamics models are supported as well.
NetLogo is developed and maintained at the Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.
More details on NetLogo itself are available online: [NetLogo online documentation](https://ccl.northwestern.edu/netlogo/docs/)

NetLogo comes with the built-in experiment tool [Behavior Space](https://ccl.northwestern.edu/netlogo/docs/behaviorspace.html) that allows to setup and execute model simulations with different settings and parameter variations and to collect model output. This experiment tool can be executed via command line in combination with an XML file that contains the experiment specifications, such as runtime, variables, output measurements, stop conditions, and more. One limitation of Behavior Space is, that it only supports full-factorial parameter designs, which may not be appropriate for complex model analyses. Furthermore, Behavior Space experiment specifications are stored within the NetLogo file and are not easily accessible from R. However, in many cases it is useful to store such specifications along with the model output and analyses results in order to enable fully reproducible model analyses.

The nlrx package utilizes the commandline functionality of Behavior Space to execute NetLogo simulations directly from R.
Instead of defining experiments within NetLogo Behavior Space, experiments are defined in R using the class objects of the nlrx package.
These class objects hold all the information that is needed to run these experiments remotely from R, such as path to NetLogo installation folder, path to the model file and the experiment specifications itself. nlrx provides useful helper functions to generate parameter input matrices from parameter range definitions that cover a wide range of parameter exploration approaches. By storing all relevant information on simulation experiments, including the output of the model simulations in one class object, experiments can be easily stored and shared.

In summary, the nlrx package uses a similar structure as NetLogos Behavior Space but offers more flexibility and additional tools for running reproducible complex model analyses directly from R.

## Publication
Further information on the package functionality and detailed code examples can be found in our accompanying publication:
Salecker J, Sciaini M, Meyer KM, Wiegand K. The nlrx r package: A next-generation framework for reproducible NetLogo model analyses. Methods Ecol Evol. 2019;2041-210X. [https://doi.org/10.1111/2041-210X.13286](https://doi.org/10.1111/2041-210X.13286).

Get citation information for `nlrx` in R doing `citation(package = 'nlrx')`.

## Prerequirements

### NetLogo
In order to use the nlrx package, [NetLogo](http://netlogoweb.org/) (>=5.3.1) needs to be installed on the system that is used to execute model simulations (local/remote). For remote execution, NetLogo needs to be installed on remote machines as well. The nlrx package provides a utility function (`download_netlogo()`) that can be used to download and unzip (only unix systems) a specified NetLogo version to a local folder. For windows machines, the downloaded file needs to be executed in order to install NetLogo on the local system. If you are running MacOS, please use the Linux tar.gz version of NetLogo (either from the NetLogo Homepage or by using the `download_netlogo()` function). The dmg version from the NetLogo homepage is not compatible with nlrx.

### Java
Because NetLogo is executed in a Java virtual machine, Java needs to be installed on the local/remote system as well. We recommend the [Oracle Java SE Development Kit 8](https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) or the [openjdk](https://github.com/ojdkbuild/ojdkbuild).
While the nlrx package might work without setting the Java system path explicitly, we recommend to make sure that JAVA_HOME points to the correct Java installation of the system.

## Installation

You can install the released version of nlrx from
[CRAN](https://CRAN.R-project.org) with:
``` r
install.packages("nlrx")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/nlrx")
```

## Get started

General information that is needed to run NetLogo simulations remotely, such as path to the NetLogo installation folder is stored within a `nl` class object.
Nested within this `nl` class are the classes `experiment` and `simdesign`. The `experiment` class stores all experiment specifications. After attaching a valid experiment, a `simdesign` class object can be attached to the `nl` class object, by using one of the simdesign helper functions. These helper functions create different parameter input matrices from the experiment variable definitions that can then be executed by the `run_nl_one()` and `run_nl_all()` functions. The nested design allows to store everything related to the experiment within one R object. Additionally, different simdesign helper functions can be applied to the same `nl` object in order to repeat the same experiment with different parameter exploration methods (simdesigns).

### Step by step application example

The "Wolf Sheep Predation" model from the NetLogo models library is used to present a basic example on how to setup and run NetLogo model simulations from R.

#### Step 1: Create a nl object:

The nl object holds all information on the NetLogo version, a path to the NetLogo directory with the defined version, a path to the model file, and the desired memory for the java virtual machine. Depending on the operation system, paths to NetLogo and the model need to be adjusted.

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

#### Step 2: Attach an experiment

The experiment object is organized in a similar fashion as NetLogo Behavior Space experiments.
It contains all information that is needed to generate a simulation parameter matrix and to execute the NetLogo simulations.
Details on the specific slots of the experiment class can be found in the package documentation (`?experiment`) and the "Advanced configuration" vignette.

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                             'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))
```

#### Step 3: Attach a simulation design

While the experiment defines the variables and specifications of the model, the simulation design creates a parameter input table based on these model specifications and the chosen simulation design method.
nlrx provides a bunch of different simulation designs, such as full-factorial, latin-hypercube, sobol, morris and eFast (see "Simdesign Examples" vignette for more information on simdesigns).
All simdesign helper functions need a properly defined nl object with a valid experiment design. Each simdesign helper also allows to define a number of random seeds that are randomly generated and can be used to execute repeated simulations of the same parameter matrix with different random-seeds (see "Advanced configuration" vignette for more information on random-seed and repetition management).
A simulation design is attached to a nl object by using one of the simdesign helper functions:

```{r eval=FALSE}
nl@simdesign <- simdesign_lhs(nl=nl,
                               samples=100,
                               nseeds=3,
                               precision=3)
```

#### Step 4: Run simulations

All information that is needed to run the simulations is now stored within the nl object.
The `run_nl_one()` function allows to run one specific simulation from the siminput parameter table.
The `run_nl_all()` function runs a loop over all simseeds and rows of the parameter input table siminput.
The loops are constructed in a way that allows easy parallelisation, either locally or on remote HPC machines (see "Advanced configuration" vignette for more information on parallelisation).
Before running your simulations you might want to check your current nl object setup. `eval_variables_constants(nl)` evaluates if the defined variables and constants are correctly defined and are consistent with the attached model. `print(nl)` prints a complete summary of the provided nl object including checkmarks that might help to indicate potential problems.


```{r eval=FALSE}
# Evaluate nl object:
eval_variables_constants(nl)
print(nl)

# Run all simulations (loop over all siminputrows and simseeds)
results <- run_nl_all(nl)
```

#### Step 5: Attach results to nl and run analysis

nlrx provides method specific analysis functions for each simulation design.
Depending on the chosen design, the function reports a tibble with aggregated results or sensitivity indices.
In order to run the analyze_nl function, the simulation output has to be attached to the nl object first.
The simdesign class within the nl object provides a slot for attaching output results (simoutput). An output results tibble can be attached to this slot by using the simdesign setter function `setsim(nl, "simoutput")`.
After attaching the simulation results, these can also be written to the defined outpath of the experiment object.

```{r eval=FALSE}
# Attach results to nl object:
setsim(nl, "simoutput") <- results

# Write output to outpath of experiment within nl
write_simoutput(nl)

# Do further analysis:
analyze_nl(nl)
``` 

## Complete code example 

This code block contains the above code example for easy copy-pasting into your R script:

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

# Setup nl object
nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

# Attach experiment
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                             'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

# Attach simdesign
nl@simdesign <- simdesign_lhs(nl=nl,
                               samples=100,
                               nseeds=3,
                               precision=3)

# Evaluate nl object:
eval_variables_constants(nl)
print(nl)

# Run all simulations (loop over all siminputrows and simseeds)
results <- run_nl_all(nl)

# Attach results to nl object:
setsim(nl, "simoutput") <- results

# Write output to outpath of experiment within nl
write_simoutput(nl)

# Do further analysis:
analyze_nl(nl)
```

---
title: "Sensitivity Analysis"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Sensitivity Analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


# Sensitivity Analysis with nlrx

Different types of sensitivity analyses can be conducted using the nlrx package.
To perform a local sensitivity analysis, we recommend using the `simdesign_distinct()` to specify local changes of parameters. Afterwards the proportion of output change can be easily calculated from the simulation results.
The nlrx package also provides simdesign helper functions to conduct more sophisticated methods such as Morris Elementary Effects Screening (`simdesign_morris()`), Sobol variance decomposition (`simdesign_sobol()`, `simdesign_sobol2007()`, `simdesign_soboljansen()`) and Extended Fourier amplitude sensitivity test (`simdesign_eFAST`). Additionally, output of Latin Hypercube Sampling designs (`simdesign_lhs()`) can be used to calculate parameter effects based on Partial (rank) correlation coefficients or Standardised (rank) regression coefficients.

In this vignette, we present an example of the Morris Elementary Effects screening. Other sensitivity analyses simdesigns work in a quite similar way. Details on the specific methods can be found in the corresponding simdesign help pages and the documentation of the [sensitivity package](https://cran.r-project.org/package=sensitivity). The second example shows how Latin Hypercube Sampling can be used to calculate Partial (rank) correlation coefficients and Standardised (rank) regression coefficients.


## Example 1: Morris elementary effects screening

Here we present a simple example for running a Morris Sensitivity Analysis with nlrx.
We use the Wolf Sheep Predation model from the models library for this example.

#### Step 1: Create a nl object:

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

#### Step 2: Attach an experiment

In this example, we want to calculate sensitivity for 3 outputs (number of sheep, number of wolves, number of grass patches).
We vary all numeric model parameters to estimate their sensitivity on the three defined output metrics.
Thus, we define parameter ranges and distribution functions for all our numeric model parameters.
We set the runtime of the model to 500 ticks and measure our metrics on each tick (`tickmetrics = "true"`).
However, for calculation of sensitivity indices, we only want to consider the last 200 ticks. Thus, we set evalticks to `seq(300,500)`.

```{r eval=FALSE}
nl@experiment <- experiment(expname = "wolf-sheep-morris",
                            outpath = outpath,
                            repetition = 1,   
                            tickmetrics = "true",
                            idsetup = "setup",  
                            idgo = "go",        
                            runtime = 500,
                            evalticks = seq(300,500),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list("initial-number-sheep" = list(min=50, max=150, step=10, qfun="qunif"),
                                             "initial-number-wolves" = list(min=50, max=150, step=10, qfun="qunif"),
                                             "grass-regrowth-time" = list(min=0, max=100, step=10, qfun="qunif"),
                                             "sheep-gain-from-food" = list(min=0, max=50, step=10, qfun="qunif"),
                                             "wolf-gain-from-food" = list(min=0, max=100, step=10, qfun="qunif"),
                                             "sheep-reproduce" = list(min=0, max=20, step=5, qfun="qunif"),
                                             "wolf-reproduce" = list(min=0, max=20, step=5, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "show-energy?" = "false"))
```

#### Step 3: Attach a simulation design

We use the `simdesgin_morris()` function to attach a Morris Sensitivity Analysis design.
The `morrislevels` parameter sets the number of different values for each parameter (sampling density).
The `morrisr` paramater sets the number of repeated samplings (sampling size).
The `morrisgridjump` parameter sets the number of levels that are increased/decreased for computing the elementary effects. Morris recommendation is to set this value to `levels / 2`.
We can increase the `nseeds` parameter in order to perform multiple runs of the same parameter matrix with different random seeds.
The variation between those repetitions is an indicator of the stochasticity effects within the model.
More information on the Morris specific parameters can be found in the description of the morris function in the sensitivity package (`?morris`).

```{r eval=FALSE}
nl@simdesign <- simdesign_morris(nl=nl,
                                 morristype="oat",
                                 morrislevels=4,
                                 morrisr=1000,
                                 morrisgridjump=2,
                                 nseeds=5)
```

#### Step 4: Run simulations

To execute the simulations, we can use the function `run_nl_all()`.
Sensitivity analyses typically have many runs that need to be simulated, thus we recommend to parallelize model runs by adjusting the future plan (more details on parallelization can be found in the "Advanced configuration" vignette).

```{r eval=FALSE}
library(future)
plan(multisession)
results <- run_nl_all(nl)
```

#### Step 5: Investigate output

First, we need to attach the results to the nl object.

```{r eval=FALSE}
setsim(nl, "simoutput") <- results
saveRDS(nl, file.path(nl@experiment@outpath, "morris.rds"))
```

After results have been attached, we can use the `analyze_nl()` function to calculate morris sensetivity indices.

```{r eval=FALSE}
morris <- analyze_nl(nl)
```


## Example 2: Latin Hypercube Sampling

Here we perform a Latin Hypercube Sampling to calculate Partial (rank) correlation coefficients and Standardised (rank) regression coefficients.


#### Step 1: Create a nl object:

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

#### Step 2: Attach an experiment

In this example, we want to calculate sensitivity for 3 outputs (number of sheep, number of wolves, number of grass patches).
We vary all numeric model parameters to estimate their sensitivity on the three defined output metrics.
Thus, we define parameter ranges and distribution functions for all our numeric model parameters.
We set the runtime of the model to 500 ticks and measure our metrics on each tick (`evalticks = "true"`).

```{r eval=FALSE}
nl@experiment <- experiment(expname = "wolf-sheep-morris",
                            outpath = outpath,
                            repetition = 1,   
                            tickmetrics = "true",
                            idsetup = "setup",  
                            idgo = "go",        
                            runtime = 500,
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list("initial-number-sheep" = list(min=50, max=150, step=10, qfun="qunif"),
                                             "initial-number-wolves" = list(min=50, max=150, step=10, qfun="qunif"),
                                             "grass-regrowth-time" = list(min=0, max=100, step=10, qfun="qunif"),
                                             "sheep-gain-from-food" = list(min=0, max=50, step=10, qfun="qunif"),
                                             "wolf-gain-from-food" = list(min=0, max=100, step=10, qfun="qunif"),
                                             "sheep-reproduce" = list(min=0, max=20, step=5, qfun="qunif"),
                                             "wolf-reproduce" = list(min=0, max=20, step=5, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "show-energy?" = "false"))
```

#### Step 3: Attach a simulation design

Here we want to run a Latin Hypercube Sampling, thus we use the `simdesign_lhs()` function.

```{r eval=FALSE}
nl@simdesign <- simdesign_lhs(nl, samples=500, nseeds=1, precision=3)
```

#### Step 4: Run simulations

To execute the simulations, we can use the function `run_nl_all()`.
Sensitivity analyses typically have many runs that need to be simulated, thus we recommend to parallelize model runs by adjusting the future plan (more details on parallelization can be found in the "Advanced configuration" vignette).

```{r eval=FALSE}
library(future)
plan(multisession)
results <- run_nl_all(nl, split=10)
```

#### Step 5: Investigate output

First, we need to attach the results to the nl object.

```{r eval=FALSE}
setsim(nl, "simoutput") <- results
saveRDS(nl, file.path(nl@experiment@outpath, "lhs.rds"))
```

After results have been attached, we need to post-process our data to run the `pcc` and `src` function of the sensitivity package.
We first take our parameter matrix (`siminput`) and select only columns with variable parameters and drop all other columns.
We also need to rename the columns because `pcc` and `src` do not support special characters (-) in column names.

Our simulation results are measured for each tick, thus we first need to aggregate our output. Here we just calculate the mean and standard deviation of outputs for each random-seed and siminputrow combination. Afterwards, we drop the random seed and siminputrow columns and rename the columns to remove special characters.

Finally, we use both datasets to run the `pcc` and `src` functions. These functions can only compute coefficients for one output at a time. Thus, we nested the function call inside a `purrr::map()` function that iterates over the column names of our output tibble.

```{r eval=FALSE}
library(tidyverse)
input <- getsim(nl, "siminput") %>%    # Take input parameter matrix
  dplyr::select(names(getexp(nl, "variables"))) %>%  # Select variable parameters only
  dplyr::rename_all(~str_replace_all(., c("-" = "_", "\\s+" = "_"))) # Remove - and space characters.

output <- getsim(nl, "simoutput") %>%   # Take simulation output
  dplyr::group_by(`random-seed`, siminputrow) %>% # Group by random seed and siminputrow
  dplyr::summarise_at(getexp(nl, "metrics"), list(mean=mean, sd=sd)) %>% # Aggregate output
  dplyr::ungroup() %>%  # Ungroup
  dplyr::select(-`random-seed`, -siminputrow) %>%  # Only select metrics
  dplyr::rename_all(~str_replace_all(., c("-" = "_", "\\s+" = "_", "\\[" = "_", "\\]" = "_", "=" = ""))) # Remove - and space characters.

# Perform pcc and src for each output separately (map)
pcc.result <- purrr::map(names(output), function(x) sensitivity::pcc(X=input, y=output[,x], nboot = 100, rank = FALSE)) 
src.result <- purrr::map(names(output), function(x) sensitivity::src(X=input, y=output[,x], nboot = 100, rank = FALSE)) 
```

The results are reported as a nested list, where each outer element represents one of the calculated model outputs. The inner list items represent the different outputs from the `pcc` and `src` functions.

We can for example look at the `pcc` results of one specific output by using the basic `plot` function:

```{r eval=FALSE}
plot(pcc.result[[1]])
```

We can also extract all the data to a tidy data format and create nice plots with the ggplot package:

```{r eval=FALSE}
pcc.result.tidy <- purrr::map_dfr(seq_along(pcc.result), function(x) {
  pcc.result[[x]]$PCC %>% 
    tibble::rownames_to_column(var="parameter") %>% 
    dplyr::mutate(metric = names(output)[x])
})

ggplot(pcc.result.tidy) +
  coord_flip() +
  facet_wrap(~metric) +
  geom_point(aes(x=parameter, y=original, color=metric)) +
  geom_errorbar(aes(x=parameter, ymin=`min. c.i.`, ymax=`max. c.i.`, color=metric), width=0.1)

src.result.tidy <- purrr::map_dfr(seq_along(src.result), function(x) {
  src.result[[x]]$SRC %>% 
    tibble::rownames_to_column(var="parameter") %>% 
    dplyr::mutate(metric = names(output)[x])
})

ggplot(src.result.tidy) +
  coord_flip() +
  facet_wrap(~metric) +
  geom_point(aes(x=parameter, y=original, color=metric)) +
  geom_errorbar(aes(x=parameter, ymin=`min. c.i.`, ymax=`max. c.i.`, color=metric), width=0.1)
```
---
title: "Advanced configuration"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Advanced configuration}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Variables and constants definitions

Correctly defining variables within the experiment class object is crucial for creating simdesigns.
The implemented simdesigns have different requirements for variable definitions:

   Simdesign       | Variable requirements              |  data type 
------------------ | ---------------------------------- | -----------
simdesign_simple   | only constants are used            | any
simdesign_distinct | values (need to have equal length) | any
simdesign_ff       | values, or min, max, step (values is prioritized) | any
simdesign_lhs      | min, max, qfun                     | numeric
simdesign_sobol    | min, max, qfun                     | numeric
simdesign_sobol2007 | min, max, qfun                     | numeric
simdesign_soboljansen | min, max, qfun                     | numeric
simdesign_morris   | min, max, qfun                     | numeric
simdesign_eFast    | min, max, qfun                     | numeric
simdesign_genSA    | min, max                          | numeric
simdesign_genAlg    | min, max                         | numeric
simdesign_ABCmcmc_Marjoram | min, max, qfun            | numeric
simdesign_ABCmcmc_Marjoram_original | min, max, qfun   | numeric
simdesign_ABCmcmc_Wegmann | min, max, qfun             | numeric


Additionally, please note the following restrictions in order to define variables and constants correctly:

* Categorical variable values are currently only allowed for simdesign_simple, simdesign_distinct and simdesign_ff.
* Variable values that should be recognized by NetLogo as strings need to be nested inside escaped quotes (e.g. `"\\"string\\""`).
* Variable values that should be recognized by NetLogo as logical need to be entered as strings (e.g. `"false"`).
* It is not allowed to list the same variable in the variables and constants list. 
* NetLogo model parameters that are not listed in any of these two lists will be set with their default value from the NetLogo model interface.
* simdesign_simple() requires at least one defined constant within the constants list. If your model does not have any globals (GUI and code), please create a dummy global (either create a global widget on the GUI or add a dummy variable to the globals[] code section) for your model and put it in the constants list with an appropriate value.

A complete list of all valid NetLogo parameters can be loaded by commiting a nl object with a valid modelpath to the function `report_model_parameters()`. This function reads all GUI elements of the NetLogo model that can be set by nlrx.
After attaching an experiment to an nl object, validity of defined experiment variables and constants can be checked by commiting the nl object to the function `eval_variables_constants()`. The function will report detailed warnings or error messages, if definitions of variables or constants are invalid.

## Print functions

`print(nl)` can be used with any nl class object. The function will print a formatted overview of class contents to the console.
The function will also print a short summary checklist that may be helpful for debugging certain issues.
Depending on the simdesign, the summary table also prints the estimated number of runs.


## Capturing progress of model simulations

The `run_nl_all()` function supports the [progressr](https://github.com/HenrikBengtsson/progressr) framework for capturing progress of simulations. Following this logic, the function itself will be silent unless it is wrapped within a `progressr::with_progress()` call. By using different handlers, the layout of the progress bar is different. For example, installing the package [progress](https://github.com/r-lib/progress) and using the handler `progressr::handlers("progress")` will additionally print the current row and seed of the parameter matrix. For example, to run a nl parameter matrix with progress bar one can do:

```{r progress, eval=FALSE}
progressr::handlers("progress")
results <- progressr::with_progress(run_nl_all(nl))
```

The `run_nl_dyn()` function might provide progress output depending on the chosen method (for example ABC offers a progress bar). However, dynamic experiments are difficult to track because it is not clear how long the complete experiment will take from the beginning.

The `run_nl_one()` does not report any progress because it only executes one simulation.

In addition, NetLogo print commands are redirected to the R console. Thus, print commands can be used within the NetLogo model code to display the current progress of simulations in the R console. Another possibility is, to define a print reporter procedure in the experiment slot "idfinal" that is executed at the end of each simulation.

Capturing output from multiple processes in parallelized environments to one R console is not straightforward. If such functionality is needed, we suggest to write the current progress to an output file directly from NetLogo (for example using the idrunnum functionality of nlrx, see section "Notes on self-written output"). These output files can then be monitored to capture the progress of the parallelized model executions. 

## Handling NetLogo runtime errors

Usually, runtime errors of NetLogo model simulations are printed to the R console and the current execution of `run_nl_all()`, `run_nl_one()` or `run_nl_dyn()` breaks. However, it can happen that a model simulation freezes due to Java runtime errors. Unfortunately it is not possible to terminate the Java virtual machine or print an error message to the console after such a runtime error occurred.
The current R session and the freezed Java Virtual Machine need to be terminated manually.
Thus, NetLogo models should be debugged in NetLogo prior to execution of large experiments with nlrx.
Capturing progress of model simulations (see section "Capturing progress of model simulations") might help in debugging runtime errors that freeze the Java Virtual Machine.


## Self-written NetLogo output

The experiment provides a slot called "idrunnum".
This slot can be used to transfer the current nlrx experiment name, random seed and runnumber (siminputrow) to NetLogo.
To use this functionality, a string input field widget needs to be created on the GUI of your NetLogo model.
The name of this widget can be entered into the "idrunnum" field of the experiment.
During simulations, the value of this widget is automatically updated with a generated string that contains the current nlrx experiment name, random seed and siminputrow ("expname_seed_siminputrow").
For self-written output In NetLogo, we suggest to include this global variable which allows referencing the self-written output files to the collected output of the nlrx simulations in R.

## Temporary files management

nlrx uses temporary files to store experiment xml files, commandline batch files to start NetLogo simulations and csv output files.
These temporary files are stored in the default temporary files directory of the R session.
By default, these files are deleted after each simulation run. However, if it is needed to look at this files, automatic deletion of temporary files can be disabled by setting the corresponding cleanup parameters in the `run_nl` functions (cleanup.csv, cleanup.xml, cleanup.bat function parameters).

On unix systems, it can happen that system processes delete files in the default temporary files folder. Thus, we recommend to reassign the temporary files folder for the R-session to another folder. The R-package [unixtools](https://www.rforge.net/unixtools/) provides a function to reassign the temporary files folder for the current R-session:
```{r eval=FALSE}
install.packages('unixtools', repos = 'http://www.rforge.net/')
unixtools::set.tempdir("<path-to-temp-dir>")
```


## random-seed and repetitions management

The experiment provides a slot called "repetition" which allows to run multiple simulations of one parameterization.
This is only useful if you manually generate a new random-seed during the setup of your model.
By default, the NetLogo random-seed is set by the simdesign that is attached to your nl object.
If your model does not reset the random seed manually, the seed will always be the same for each repetition.

However, the concept of nlrx is based on sensitivity analyses. Here, you may want to exclude stochasticity from your output and instead do multiple sensitivity analyses with the same parameter matrix but different random seeds. You can then observe the effect of stochasticity on the level of your final output, the sensitivity indices. Thus we suggest to set the experiment repetition to 1 and instead use the nseeds variable of the desired simdesign to run multiple simulations with different random seeds.

In summary, if you set the random-seed of your NetLogo model manually, you can increase the repitition of the experiment to run several simulations with equal parameterization and different random seeds.
Otherwise, set the experiment repetition to 1 and increase the nseeds variable of your desired simdesign.

## Runtime and measurements

The runtime of NetLogo model simulations can be controlled with two slots of the experiment class:

* runtime - an integer that defines the maximum number of ticks that are simulated
* stopcond - a NetLogo reporter that reports true/false. Simulations are automatically stopped if TRUE is reported. Defaults to NA_character_ which means, no stop condition is applied.

Runtime can be set to 0 or NA_integer_ to allow for simulations without a pre-defined maximum runtime.
However, this should only be done in combination with a proper stop condition (stopcond) or with NetLogo models that have a built-in stop condition. Otherwise, simulations might get stuck in endless loops.

Two slots of the experiment class further define when measurements are taken:

* tickmetrics - defines if measurements are taken at the end of the simulation (final tick) or on each tick
* evalticks - only applied when tickmetrics = "true"; defines a vector of ticks for which output metrics are reported. Set to NA_integer_ to report all collected output.

Depending on the evalticks definition, it might happen, that a simulation stops before any output has been collected.
In such cases, output is still reported but all metrics that could not be collected for any defined evalticks will be filled up with NA.

Four slots of the experiment class further define which measurements are taken:

* metrics - vector of valid NetLogo reporters that are used to collect data (e.g. c("count turtles"))
* metrics.turtles - a list with named vectors of strings defining valid turtles-own variables that are taken as output measurements (e.g. list("turtles" = c("who", "pxcor", "pycor", "color"))
* metrics.patches - vector of valid patches-own variables that are used to collect patches data (e.g. c("pxcor", "pycor", "pcolor"))
* metrics.links - a list with named vectors of strings defining valid links-own variables that are taken as output measurements (e.g. list("links" = c("[who] of end1", "[who] of end2")))

Although the metrics slot accepts any valid NetLogo reporter, such as "count patches", reporter strings can become quite long and confusing. We suggest to create NetLogo reporter procedures for complex reporters in order to get a nice and clean results data frame.
For example, the NetLogo reporter "count patches with [pcolor = green]" could be written as a NetLogo reporter function:
```{r eval=FALSE}
to-report green.patches
  report count patches with [pcolor = green]
end
```
In your nlrx experiment metrics field you can then enter "green.patches" which is way more intuitive then "count patches with [pcolor = green]".

## NetLogo extensions

Usually, all NetLogo extensions that are shipped with NetLogo should also work with nlrx.
However, depending on the system it can happen that NetLogo extensions are not found properly. To solve such problems we advise to put your `.nlogo model` file in the `app/models` subdirectory of the NetLogo installation path.
A special case is the NetLogo r-extension because it needs to be stopped manually in between model runs.
To achieve that, simply put the `r:stop` command in the `idfinal` slot of your experiment: `idfinal = "r:stop"`.


## Parallelisation and the future concept

The run_nl_all function uses the `future_map_dfr()` function from the [furrr package](https://github.com/DavisVaughan/furrr). The simulations are executed in a nested loop where the outer loop iterates over the random seeds of your simdesign, and the inner loop iterates over the rows of the siminput parameter matrix of your simdesign. These loops can be executed in parallel by setting up an appropriate plan from the [future package](https://github.com/HenrikBengtsson/future). See examples below for more details on parallelisation on local machines and remote HPC clusters.

### Parallelisation on local machine

Model simulations can be distributed to each logical processor on the local machine in parallel. The future package provides two options for parallelization, explicit futures and implicit futures which are executed in the background and do not block the console while running. 

Running parallel simulations with an explicit future command:
```{r eval=FALSE}
library(future)
plan(multiprocess)
results <- run_nl_all(nl = nl)
```

For running parallel simulations with an implicit future command we need to define the type of parallelisation for each level of the nested `furrr::future_map_dfr()` function individually. Because we want to parallelize the actual simulations, we need to use the multiprocess plan on the inner level.
Note, that we use the assignment operator for implicit futures (`%<-%`):
```{r eval=FALSE}
library(future)
plan(list(sequential, multiprocess))
results %<-% run_nl_all(nl = nl)
```

In cases, where the number of random seeds is lower than the available processor cores, parallelisation may not be completely efficient. To allow efficient parallelisation, even for a small number of random seeds the split parameter of the `run_nl_all()` function can be used to split the parameter matrix into smaller chunks, which can be distributed to separate processor cores. For example, a simulation with 1000 runs (rows of the siminput matrix) and 2 random seeds should be distributed to 8 processor cores. By default, the parallelisation loop would consist of two jobs (one for each random seed) with 1000 simulation runs each. This experiment would only utilize 2 of the 8 available processor cores. By setting the split parameter to 4, we increase the total number of jobs from 2 to 8 (2 random-seeds * 4 parameter matrix chunks). Each job runs 1/4th of the parameter input matrix (250 rows) using one of the 2 defined random seeds.
```{r eval=FALSE}
library(future)
plan(multisession)
results <- run_nl_all(nl = nl, split = 4)
```



### Parallelisation on remote HPC cluster

This option requires access to a remote HPC cluster.
This example gives you some guidance and examples for sending jobs from an R session on your local machine to an HPC running slurm. Details might be different depending on the HPC setup.

Some settings, such as ssh access and slurm templates need to be defined to access remote HPC clusters from your local R session. Please check out this [detailed HPC setup manual](https://r-spatialecology.github.io/gwdg_hpc_guide/) for examples on required settings for an HPC running slurm.

In order to run NetLogo models on remote HPCs, required software needs to be installed on the remote system as well: java, Netlogo, R, nlrx and further required R-packages (future/clusterMQ, ...). Of course the NetLogo model files need to be available on the remote machine as well.

#### Using the future framework

For sending jobs to the remote HPC under the future framework, we need to install and load additional R packages and adjust the future plan before executing `run_nl_all()`.
You need to define a path to your ssh key, a server adress for the HPC, a user name for the HPC and a path to the slurm template file on the HPC.
Please also make sure that the `nlpath`, `modelpath` and `outpath` variables within your nl object point to locations on the filesystem of the HPC and not your local filesystem.

```{r eval=FALSE}
# Load required packages
library(future)
library(future.batchtools)
library(debugme)
Sys.setenv(DEBUGME='batchtools')
library(batchtools)

# Define the path to the ssh key for your machine:
options(future.makeNodePSOCK.rshopts = c("-i", "/patch/to/id_rsa"))
# Define server and your login credentials for the remote HPC:
login <- tweak(remote, workers="server.HPC.de", user="username")

# Define plan for future environment:
bsub <- tweak(batchtools_slurm, template = "slurm.tmpl", # define name of slurm tmeplate on HPC filesystem
              resources = list(job.name = "jobname", # define jobname
                               log.file = "jobname.log", # define logfile name
                               queue = "medium",  # define HPC queue
                               service = "normal", # define HPC service
                               walltime = "00:30", # define walltime
                               n_jobs = "1",   # define number of processes per job
                               mem_cpu = "4000") # define memory per cpu   

# Load HPC plan:
plan(list(login,
          bsub,
          multisession))

# Execute simulations
results <- run_nl_all(nl = nl)
```

#### Using the clusterMQ framework

The clusterMQ framework is somewhat different from the future framework. However, in our experience it worked more reliable in combination with a slurm HPC.
For installation of clustermq and .Rprofile settings see also the [detailed HPC setup manual](https://r-spatialecology.github.io/gwdg_hpc_guide/).

clusterMQ does not directly support parallelisation of the nested `furrr::future_map_dfr()` loops of the `run_nl_all()` function. We need to define our own parallel simulation function, using the `run_nl_one()` function of the nlrx package:

```{r eval=FALSE}
library(clustermq)

# First, we set the total number of jobs for the HPC
# In this example we run each simulation as an individual job (recommended).
# Thus to calculate the number of jobs we just multiply the number of parameterizations of the simdesign with the number of random seeds.
# If you want to group several runs into the same job you can adjust this line and choose a lower number.
# However, the number must be chosen that nrow(nl@simdesign@siminput)/njobs results in an integer value
njobs <- nrow(nl@simdesign@siminput) * length(nl@simdesign@simseeds)

# Second, we generate vectors for looping trough model runs.
# We generate a vector for simpinputrows by repeating the sequence of parameterisations for each seed.
# Then, we generate a vector of random-seeds by repeating each seed for n times, where n is the number of siminputrows.
siminputrows <- rep(seq(1:nrow(nl@simdesign@siminput)), length(nl@simdesign@simseeds))
rndseeds <- rep(nl@simdesign@simseeds, each=nrow(nl@simdesign@siminput))

# Third, we define our simulation function
# Please adjust the path to the temporary file directory
simfun <- function(nl, siminputrow, rndseed, writeRDS=FALSE)
{
  unixtools::set.tempdir("/hpath/to/temp/dir")
  library(nlrx)
  res <- run_nl_one(nl = nl, siminputrow = siminputrow, seed = rndseed, writeRDS = TRUE)
  return(res)
}

# Fourth, use the Q function from the clustermq package to run the jobs on the HPC:
# The q function loops through our siminputrows and rndseeds vectors.
# The function creates njobs jobs and distributes corresponding chunks of the input vectors to these jobs for executing the simulation function simfun.
# As constants we provide our nl object and the function parameter writeRDS. 
# If write RDS is true, an *.rds file will be written on the HPC after each jobs has finished.
# This can be useful to gain results from completed runs after a crash has occured.
results <- clustermq::Q(fun = simfun, 
                        siminputrow = siminputrows,
                        rndseed = rndseeds,
                        const = list(nl = nl,
                                     writeRDS = TRUE),
                        export = list(), 
                        seed = 42, 
                        n_jobs = njobs, 
                        template = list(job_name = "jobname", # define jobname
                                        log.file = "jobname.log", # define logfile name
                                        queue = "medium",  # define HPC queue
                                        service = "normal", # define HPC service
                                        walltime = "00:30:00", # define walltime
                                        mem_cpu = "4000")) # define memory per cpu   

# The Q function reports the individual results of each job as a list
# Thus, we convert the list results to tibble format:
results <- dplyr::bind_rows(results)
```
---
title: "Approximate Bayesian Computation (ABC)"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Approximate Bayesian Computation (ABC)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


# Approximate bayesian computation (ABC) with nlrx

Approximate bayesian computation (ABC) algorithms have been increasingly used for calibration of agent-based simulation models.
The nlrx package provides different algorithms from the [EasyABC package](https://cran.r-project.org/package=EasyABC).
These algorithms can be used by attaching the corresponding simdesigns (`simdesign_ABCmcmc_Marjoram()`, `simdesign_ABCmcmc_Marjoram_original()`, `simdesign_ABCmcmc_Wegmann()`). Example 1 shows the process of how to use ABC with nlrx.
Additionally, Latin Hypercube Sampling output can be used to calculate parameter distributions based on rejection sampling and local linear regression.
Example 2 shows, how the `simdesign_lhs()` can be used in combination with the [abc package](https://cran.r-project.org/package=abc).

## Example 1: Approximate bayesian computation with Monte-Carlo Markov-Chain

Here we present one example for the widely used Marjoram algorithm which combines ABC with a Markov-Chain Monte-Carlo parameter sampling scheme.
However, the other two ABCmcmc simdesigns work in a very similar way except for the parameter definitions within the simdesigns (see respective documentation pages for help).

We use the Wolf Sheep Predation model from the models library to show a basic example of the calibration workflow.

#### Step 1: Create a nl object:

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

#### Step 2: Attach an experiment

Because we want to apply a calibration algorithm, we need to define proper variable ranges.
The algorithm is allowed to change the values of these parameters within these ranges in order to reach our specified target values.
Possible choices for the distribution type are qunif, qnorm, qlnorm and qexp.
Each model run is evaluated for each specified metric within the metrics slot. When the simdesign is attached (see step 3) we will define target values for each metric. Thus, we should only enter metrics and reporters that should be used for calibrating the model.

For this simple example, we just want to check if we can find a parameterization that leads to a specific number of wolves and sheep after a runtime of 100 ticks. Thus, we define `count sheep` and `count wolves` as metrics, set the runtime to 100 and set tickmetrics to `false` (because we only want to measure the last simulation tick).

If more than one tick would be measured, the algorithm automatically calculates the mean value of the selected reporter over all measured ticks. If you wish to apply other functions to aggregate temporal information into one value, you can use a self-defined post-processing function when attaching the simdesign (see step 3).

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="false",
                            idsetup="setup",
                            idgo="go",
                            runtime=100,
                            metrics=c("count sheep", "count wolves"),
                            variables = list("sheep-gain-from-food" = list(min=2, max=6, qfun="qunif"),
                                             "wolf-gain-from-food" = list(min=10, max=30, qfun="qunif")),
                            constants = list('initial-number-sheep' = 100,
                                             'initial-number-wolves' = 50,
                                             "grass-regrowth-time" = 30,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "model-version" = "\"sheep-wolves-grass\"",
                                             "show-energy?" = "false"))

```

#### Step 3: Attach a simulation design

We use the `simdesign_ABCmcmc_Marjoram()` function to attach a calibration simdesign.
The `summary_stat_target` represents the vector of our target values that we want to reach. These values corresponds to the defined metrics of the experiment and should have the same length and order. `n_rec` defines the number of samples and `n_calibration` defines the number of calibration runs. If this value is too low, a subscript out of bounds error message might appear. If `use_seed` is set to `TRUE`, the algorithm will automatically use a newly generated seed for each model run. If it is set to false, a user-specified seed (set in the `run_nl_dyn()` function) will be used instead. The `progress_bar` gives you expected runtime information during the execution of the algorithm. The nseeds command allows to generate a vector of random-seeds that may be used for setting model seeds.

```{r eval=FALSE}
nl@simdesign <- simdesign_ABCmcmc_Marjoram(nl=nl,
                                           summary_stat_target = c(100, 80),
                                           n_rec = 100, 
                                           n_calibration=200,
                                           use_seed = TRUE,
                                           progress_bar = TRUE,
                                           nseeds = 1)
```


These are the most important simdesign parameters, but there are many more to fine control the behavior of the algorithm. Check out the simdesign help page for more information. As already mentioned before, it is also possible to define a custom post-processing function. To apply this function for model output, the function name needs to be entered as `postpro_function` within the simdesign. For example, we might want to use the maximum value of some measured metrics to calibrate our model. Or maybe we want to run some tests and use the test statistics as calibration criterion. In such cases, we can define a function that accepts the nl object (with simulation results attached) as function input and returns a vector of numerics. This vector needs to represent the values that we defined as `sumary_stat_target` and thus should have the same length and order. Below is an example of a custom post-processing function that calculates the maximum value of selected metrics over all simulation ticks.

```{r eval=FALSE}
post <- function(nl){
  res <- getsim(nl, "simoutput") %>% 
    dplyr::select(getexp(nl, "metrics")) %>% 
    dplyr::summarise_each(list(max=max))
  return(as.numeric(res))
}

nl@simdesign <- simdesign_ABCmcmc_Marjoram(nl=nl,
                                           postpro_function = post,
                                           summary_stat_target = c(100, 80),
                                           n_rec = 100, 
                                           n_calibration=200,
                                           use_seed = TRUE,
                                           progress_bar = TRUE,
                                           nseeds = 1)
```

#### Step 4: Run simulations

For calibration simdesigns, the `run_nl_dyn()` function lets you execute the simulations.
There are some notable differences between `run_nl_all()` and `run_nl_dyn()`.
First, because parameterizations depend of results from previous runs, `run_nl_dyn()` can not be parallelized.
Second, the procedure does not automatically loop over created random seeds of the simdesign.
If you want to repeat the same algorithm several times, just embed the `run_nl_dyn()` function in any kind of loop and iterate through the `nl@simdesign@simseeds` vector. We set the `use_seed` parameter of the simdesign to true, however we still need to define a random-seed in run_nl_dyn although it will be overwritten by the EasyABC functions.

```{r eval=FALSE}
results <- run_nl_dyn(nl, seed = nl@simdesign@simseeds[1])
```

#### Step 5: Investigate output

The output is reported as nested tibble which can be attached to the nl object.
There are many possible ways to inspect the simulation output of the `ABCmcmc` functions.
Below you find some guidance on how to summarize the output for calculating parameter statistics, sampling distributions, sampling density and exporting the best parameter combination.

```{r eval=FALSE}
setsim(nl, "simoutput") <- results
saveRDS(nl, file.path(nl@experiment@outpath, "ABCmcmc.rds"))

## Calculate descriptive statistics
getsim(nl, "simoutput") %>% # get simulation results from nl object
  dplyr::select(param) %>% # select param column
  tidyr::unnest(cols=param) %>%  # unnest param column
  dplyr::summarise_each(list(min=min, max=max, mean=mean, median=median)) %>% # calculate statistics
  tidyr::gather(parameter, value) %>% # convert to long format
  tidyr::separate(parameter, into = c("parameter", "stat"), sep = "_") %>% # seperate parameter name and statistic
  tidyr::spread(stat, value) # convert back to wide format

## Plot histogram of parameter sampling distribution:
getsim(nl, "simoutput") %>% # get simulation results from nl object
  dplyr::select(param) %>% # select param column
  tidyr::unnest(cols=param) %>%  # unnest param column
  tidyr::gather(parameter, value) %>% # convert to long format
  ggplot2::ggplot() + # plot histogram with a facet for each parameter
  ggplot2::facet_wrap(~parameter, scales="free") +
  ggplot2::geom_histogram(ggplot2::aes(x=value), bins = 40)

## Plot density of parameter sampling distribution:
getsim(nl, "simoutput") %>% # get simulation results from nl object
  dplyr::select(param) %>% # select param column
  tidyr::unnest(cols=param) %>% # unnest param column
  tidyr::gather(parameter, value) %>% # convert to long format
  ggplot2::ggplot() + # plot density with a facet for each parameter
  ggplot2::facet_wrap(~parameter, scales="free") +
  ggplot2::geom_density(ggplot2::aes(x=value, fill=parameter))

## Get best parameter combinations and corresponding function values
getsim(nl, "simoutput") %>%  # get simulation results from nl object
  dplyr::select(dist,epsilon) %>%  # select dist and epsilon columns
  tidyr::unnest(cols=c(dist,epsilon)) %>%  # unnest dist and epsilon columns
  dplyr::mutate(runID=dplyr::row_number()) %>% # add row ID column
  dplyr::filter(dist == epsilon) %>% # only keep runs with dist=epsilon
  dplyr::left_join(getsim(nl, "simoutput") %>% # join parameter values of best runs
                     dplyr::select(param) %>%
                     tidyr::unnest(cols=param) %>% 
                     dplyr::mutate(runID=dplyr::row_number())) %>% 
  dplyr::left_join(getsim(nl, "simoutput") %>% # join output values best runs
                     dplyr::select(stats) %>%
                     tidyr::unnest(cols=stats) %>% 
                     dplyr::mutate(runID=dplyr::row_number())) %>% 
  dplyr::select(runID, dist, epsilon, dplyr::everything()) # update order of columns

## Analyse mcmc using coda summary and plot functions:
summary(coda::as.mcmc(getsim(nl, "simoutput") %>%
                        dplyr::select(param) %>%
                        tidyr::unnest(cols=param)), quantiles =c(0.05,0.95,0.5))

plot(coda::as.mcmc(getsim(nl, "simoutput") %>%
                        dplyr::select(param) %>%
                        tidyr::unnest(cols=param)))

```

## Example 2: Using Latin Hypercube Sampling for Approximate bayesian computation


#### Step 1: Create a nl object:

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

#### Step 2: Attach an experiment

We want to run a Latin Hypercube sampling, thus we need to define proper variable ranges.
We also need to define our output metrics. These metrics are also used for the rejection sampling later on.

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="false",
                            idsetup="setup",
                            idgo="go",
                            runtime=100,
                            metrics=c("count sheep", "count wolves"),
                            variables = list("sheep-gain-from-food" = list(min=2, max=6, qfun="qunif"),
                                             "wolf-gain-from-food" = list(min=10, max=30, qfun="qunif")),
                            constants = list('initial-number-sheep' = 100,
                                             'initial-number-wolves' = 50,
                                             "grass-regrowth-time" = 30,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "model-version" = "\"sheep-wolves-grass\"",
                                             "show-energy?" = "false"))

```

#### Step 3: Attach a simulation design

We use the `simdesign_lhs()` helper function to generate a Latin Hypercube Sampling with 500 samples

```{r eval=FALSE}
nl@simdesign <- simdesign_lhs(nl, 
                              samples=500, 
                              nseeds=1, 
                              precision=3)
```

#### Step 4: Run simulations

We can simply use `run_nl_all()` to execute our simulations.

```{r eval=FALSE}
results <- run_nl_all(nl)
```

#### Step 5: Investigate output

We first attach the output results to our nl object and store a copy of the nl object on disk.

```{r eval=FALSE}
setsim(nl, "simoutput") <- results
saveRDS(nl, file.path(nl@experiment@outpath, "ABClhs.rds"))
```

For post-processing, we need a tibble with input parameter distributions. This can be easily extracted from the `siminput` slot of the `simdesign` by selecting only the columns with variable (non-constant) parameters. Next, we need corresponding outputs for these parameters. In this example we just take the measured output tibble (`simoutput`) and only select the columns with our metrics. Of course, you can also perform additional post-processing of these outputs if desired.
Third, we need to define expected values for our outputs. In our case, we just assume that `count sheep` should have a value of 100, whereas `count wolves` should have a value of 80.

```{r eval=FALSE}
input <- getsim(nl, "siminput") %>% 
  dplyr::select(names(getexp(nl, "variables")))
output <- getsim(nl, "simoutput") %>% 
  dplyr::select(getexp(nl, "metrics"))
target <- c("count sheep"=100, "count wolves"=80)
```

We use the `abc` function of the [abc package](https://cran.r-project.org/package=abc) to perform the rejection sampling. For this example, we perform both algorithms provided by this function ("rejection" and "loclinear").

```{r eval=FALSE}
results.abc.reject <- abc::abc(target=target, 
                        param=input,
                        sumstat=output,
                        tol=0.3, 
                        method="rejection")

results.abc.loclin <- abc::abc(target=target, 
                               param=input,
                               sumstat=output,
                               tol=0.3, 
                               method="loclinear")
```

Finally, we might want to compare the accepted parameter distributions of both algorithms and the initial distribution of the Latin Hypercube sampling.
Thus, we reformat the results to a tidy data format and attach the initial parameter distributions. This dataset can now be used for displaying the parameter distributions with ggplot.

```{r eval=FALSE}
results.abc.all <- tibble::as_tibble(results.abc.reject$unadj.values) %>% # results from rejection method
  tidyr::gather(parameter, value) %>% 
  dplyr::mutate(method="rejection") %>% 
  dplyr::bind_rows(tibble::as_tibble(results.abc.loclin$adj.values) %>% # results from local linear regression method
                     tidyr::gather(parameter, value) %>% 
                     dplyr::mutate(method="loclinear")) %>% 
  dplyr::bind_rows(input %>%                # initial parameter distribution (lhs)
                     tidyr::gather(parameter, value) %>% 
                     dplyr::mutate(method="lhs"))

ggplot2::ggplot(results.abc.all) +
  ggplot2::facet_grid(method~parameter, scales="free") +
  ggplot2::geom_histogram(ggplot2::aes(x=value))

ggplot2::ggplot(results.abc.all) +
  ggplot2::facet_wrap(~parameter, scales="free") +
  ggplot2::geom_density(ggplot2::aes(x=value, fill=method), alpha=0.1)

```

---
title: "Optimization"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Optimization}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


## Optimization with nlrx

Here we present two simple examples for running an optimization algorithm on a NetLogo model with nlrx.
In our example, we use the Simulated Annealing simdesign (`simdesign_GenSA()`).
However, except for the parameter definitions in the simdesign function and the output of the function, the genetic algorithm optimization (`simdesign_GenAlg()`) works in the same way.

We use the Wolf Sheep Predation model from the models library to show a basic example of the optimization workflow.
Example 1 shows, how a NetLogo reporter can be used as a fitness criterion for optimization. Example 2 uses a self-defined evaluation function that calculates landscape metrics that are then used as fitness criterion.

## Example 1: NetLogo reporter as fitness criterion

#### Step 1: Create a nl object:

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

#### Step 2: Attach an experiment

Because we want to apply an optimization algorithm, we need to define proper variable ranges.
The algorithm is allowed to change the values of these parameters within these ranges in order to minimize our fitness criterion.
In this example we want to use a reporter from the metrics slot for evaluating our model runs.
Here we want to find a parameterization that leads to the maximum number of wolfs after 50 ticks.
Because the algorithm automatically searches for minimum values, we add `"1 / count wolves"` to the metrics vector in order to find the maximum number of wolves.

It is also important to think about the settings for tickmetrics, runtime and evalticks.
Because we only want to consider the last tick of the simulation, we set tickmetrics to "false" and runtime to 50. If more than one tick would be measured, the algorithm automatically calculates the mean value of the selected reporter. If you wish to apply other functions to aggregate temporal information into one value, you can use a self-defined evaluation function (see Example 2).

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep-GenSA1",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="false",
                            idsetup="setup",
                            idgo="go",
                            runtime=50,
                            metrics=c("(1 / count wolves)"),
                            variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                             'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

```

#### Step 3: Attach a simulation design

We use the `simdesgin_GenSA()` function to attach a Simulated Annealing simdesign.
We select the evaluation criterion (`evalcrit`) by assigning the position of the reporter that we want to evaluate within the metrics vector of the experiment.
In our case, there is only one reporter in the metrics vector thus we set evalcrit to use the first reporter (`evalcrit = 1`).
The control parameter allows us to provide additional parameters for the GenSA function (see ?GenSA for details). For demonstration purposes, we set the maximum number of iterations to 20.

```{r eval=FALSE}
nl@simdesign <- simdesign_GenSA(nl, 
                                evalcrit = 1, 
                                nseeds = 1, 
                                control=list(maxit = 20))
```

#### Step 4: Run simulations

For optimization simdesign, the `run_nl_dyn()` function lets you execute the simulations.
There are some notable differences between `run_nl_all()` and `run_nl_dyn()`.
First, because parameterizations depend of results from previous runs, `run_nl_dyn()` can not be parallelized.
Second, the procedure does not automatically loop over created random seeds of the simdesign.
If you want to repeat the same algorithm several times, just embed the `run_nl_dyn()` function in any kind of loop and iterate through the `nl@simdesign@simseeds` vector.
Third, the output of `run_nl_dyn()` is reported as objects from the specific optimization procedures and not in tibble format. In order to attach these results to the nl object, the output needs to be converted to tibble format first. However, attaching optimization results to the nl does not enable any further post-processing functions of the nlrx package and is only relevant for storing results together with the nl object. This design decision was made in order to allow application of the method specific summary functions to the results of the optimization.

```{r eval=FALSE}
results <- run_nl_dyn(nl, seed = nl@simdesign@simseeds[1])
```

#### Step 5: Investigate output

The output list of the Simulated Annealing procedure contains four elements:
`value` reports the minimum final value of the evaluation criterion.
`par` reports the parameter settings of the final parameterisation in the same order as defined in the experiment of the nl object.
`trace.mat` gives you detailed information on the optimization process over all iterations.
`counts` indicates how often the optimization procedure was executed in total.

```{r eval=FALSE}
results
```

In order to store our results together with the nl object we need to attach the results to the nl object first.
As explained above, we need to enframe the results as a tibble.

```{r eval=FALSE}
setsim(nl, "simoutput") <- tibble::enframe(results)
saveRDS(nl, file.path(nl@experiment@outpath, "genSA_1.rds"))

```



## Example 2: Evaluation function as fitness criterion

#### Step 1: Create a nl object:

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

#### Step 2: Attach an experiment

Because we want to apply an optimization algorithm, we need to define proper variable ranges.
The algorithm is allowed to change the values of these parameters within these ranges in order to minimize our fitness criterion.
In this example we want to use a self-defined evaluation function to calculate a fitness criterion.
Thus, we add the patch coordinates and patch color (as a patch class indicator) to the `metrics.patches` vector.
We want to use spatial data to calculate the landscape edge density index of the final tick and find a parameterization that leads to the edge density.
Because we only want to consider the last tick of the simulation, we set tickmetrics to "false" and runtime to 50. 

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep-GenSA2",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="false",
                            idsetup="setup",
                            idgo="go",
                            runtime=50,
                            metrics.patches = c("pxcor", "pycor", "pcolor"),
                            variables = list('initial-number-sheep' = list(min=50, max=150),
                                             'initial-number-wolves' = list(min=50, max=150)),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

```

#### Step 3: Attach a simulation design

We use the `simdesgin_GenSA()` function to attach a Simulated Annealing simdesign.
Because we want to post-process our simulation results, we need to define an evaluation function.
The evaluation function needs to accept the nl object as input and must return a single numeric value.
First we load the package landscapemetrics. We then convert the spatial data to a raster format and calculate the landscape edge density index.
Finally, we report only the index value of the resulting tibble.

```{r eval=FALSE}
critfun <- function(nl) {
  library(landscapemetrics)
  res_spat <- nl_to_raster(nl)
  res_spat_raster <- res_spat$spatial.raster[[1]]
  lsm <- lsm_l_ed(res_spat_raster)
  crit <- lsm$value
  return(crit)
}
```

In the `simdesign_GenSA()` function we now provide our evaluation function (`critfun`) as evaluation criterion (`evalcrit`).
The control parameter allows us to provide additional parameters for the GenSA function (see ?GenSA for details). For demonstration purposes, we set the maximum number of iterations to 20.

```{r eval=FALSE}
nl@simdesign <- simdesign_GenSA(nl, 
                                evalcrit = critfun, 
                                nseeds = 1, 
                                control=list(maxit = 20))
```

#### Step 4: Run simulations

For optimization simdesign, the `run_nl_dyn()` function lets you execute the simulations.
There are some notable differences between `run_nl_all()` and `run_nl_dyn()`.
First, because parameterizations depend of results from previous runs, `run_nl_dyn()` can not be parallelized.
Second, the procedure does not automatically loop over created random seeds of the simdesign.
If you want to repeat the same algorithm several times, just embed the `run_nl_dyn()` function in any kind of loop and iterate through the `nl@simdesign@simseeds` vector.
Third, the output of `run_nl_dyn()` is reported as objects from the specific optimization procedures and not in tibble format. In order to attach these results to the nl object, the output needs to be converted to tibble format first. However, attaching optimization results to the nl does not enable any further post-processing functions of the nlrx package and is only relevant for storing results together with the nl object. This design decision was made in order to allow application of the method specific summary functions to the results of the optimization.

```{r eval=FALSE}
results <- run_nl_dyn(nl, seed = nl@simdesign@simseeds[1])
```

#### Step 5: Investigate output

The output list of the Simulated Annealing procedure contains four elements:
`value` reports the minimum final value of the evaluation criterion.
`par` reports the parameter settings of the final parameterisation in the same order as defined in the experiment of the nl object.
`trace.mat` gives you detailed information on the optimization process over all iterations.
`counts` indicates how often the optimization procedure was executed in total.

```{r eval=FALSE}
results
```

In order to store our results together with the nl object we need to attach the results to the nl object first.
As explained above, we need to enframe the results as a tibble.

```{r eval=FALSE}
setsim(nl, "simoutput") <- tibble::enframe(results)
saveRDS(nl, file.path(nl@experiment@outpath, "genSA_2.rds"))

```






---
title: "Simdesign examples"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Simdesign examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

nlrx provides a variety of different simdesigns that generate parameter input matrices from valid experiment variable definitions.
These simdesigns represent different techniques of parameter space explorations.
This is not a technical guide on parameter space exploration but rather a collection of code examples on how to setup these different approaches with nlrx.
Details on the sensitivity analysis methods can be found in the documentation of the [sensitivity package](https://CRAN.R-project.org/package=sensitivity).
Details on the optimization methods can be found in the documentation of the [genalg package](https://CRAN.R-project.org/package=genalg) and the [GenSA package](https://CRAN.R-project.org/package=GenSA). Details on the approximate bayesian computation (ABC) methods can be found in the documentation of the [EasyABC package](https://cran.r-project.org/package=EasyABC)

## Simdesign examples for the Wolf Sheep Predation model

The following section provides valid experiment setups for all supported simdesigns using the Wolf Sheep Model from the NetLogo models library.

First we set up a nl object. We use this nl object for all simdesigns:

```{r eval=FALSE}
library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
```

## Simple simulation using constants (simdesign_simple)

The simple simdesign only uses defined constants and reports a parameter matrix with only one parameterization.
To setup a simple simdesign, no variables have to be defined. The simple simdesign can be used to test a specific parameterisation of the model.
It is also useful to for creating animated output of one specific parameterisation (see "Capturing Spatial NetLogo Output" vignette).

`simdesign_simple` is most useful, if you want to imitate pressing the _go_ button in your model in NetLogo - just a single run, with a
defined parameterisation.

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath="C:/out/",
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            idfinal=NA_character_,
                            idrunnum=NA_character_,
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list(),
                            constants = list("initial-number-sheep" = 20,
                                             "initial-number-wolves" = 20,
                                             "model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

nl@simdesign <- simdesign_simple(nl=nl,
                                 nseeds=3)
```

## Distinct parameter combinations (simdesign_distinct)

The distinct simdesign can be used to run distinct parameter combinations. To setup a distinct simdesign, vectors of values need to be defined for each variable. These vectors must have the same number of elements across all variables. The first simulation run consist of all 1st elements of these variable vectors; the second run uses all 2nd values, and so on.

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath="C:/out/",
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            idfinal=NA_character_,
                            idrunnum=NA_character_,
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(values=c(10, 20, 30, 40)),
                                             'initial-number-wolves' = list(values=c(30, 40, 50, 60))),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

nl@simdesign <- simdesign_distinct(nl=nl,
                                   nseeds=3)
```

## Full-factorial simulation (simdesign_ff)

The full factorial simdesign creates a full-factorial parameter matrix with all possible combinations of parameter values.
To setup a full-factorial simdesign, vectors of values need to be defined for each variable. Alternatively, a sequence can be defined by setting min, max and step. However, if both (values and min, max, step) are defined, the values vector is prioritized.

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath="C:/out/",
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            idfinal=NA_character_,
                            idrunnum=NA_character_,
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(values=c(10, 20, 30, 40)),
                                             'initial-number-wolves' = list(min=0, max=50, step=10)),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

nl@simdesign <- simdesign_ff(nl=nl,
                             nseeds=3)
```

## Latin Hypercube Sampling (simdesign_lhs)

The latin hypercube simdesign creates a Latin Hypercube sampling parameter matrix.
The method can be used to generate a near-random sample of parameter values from the defined parameter distributions.
More Details on Latin Hypercube Sampling can be found in [McKay 1979](https://www.tandfonline.com/doi/abs/10.1080/00401706.1979.10489755). 
nlrx uses the [lhs](https://CRAN.R-project.org/package=lhs/index.html) package to generate the Latin Hypercube parameter matrix.
To setup a latin hypercube sampling simdesign, variable distributions need to be defined (min, max, qfun).

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath="C:/out/",
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            idfinal=NA_character_,
                            idrunnum=NA_character_,
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                             'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

nl@simdesign <- simdesign_lhs(nl=nl,
                               samples=100,
                               nseeds=3,
                               precision=3)

```




## Sensitivity Analyses (simdesign_sobol, _sobol2007, _soboljansen, _morris, _eFast)

Sensitivity analyses are useful to estimate the importance of model parameters and to scan the parameter space in an efficient way.
nlrx uses the [sensitivity](https://CRAN.R-project.org/package=sensitivity/index.html) package to setup sensitivity analysis parameter matrices.
All supported sensitivity analysis simdesigns can be used to calculate sensitivity indices for each parameter-output combination. These indices can be calculated by using the `analyze_nl()` function after attaching the simulation results to the nl object. 
To setup sensitivity analysis simdesigns, variable distributions (min, max, qfun) need to be defined.

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath="C:/out/",
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            idfinal=NA_character_,
                            idrunnum=NA_character_,
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                             'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

nl@simdesign <- simdesign_lhs(nl=nl,
                               samples=100,
                               nseeds=3,
                               precision=3)


nl@simdesign <- simdesign_sobol(nl=nl,
                                 samples=200,
                                 sobolorder=2,
                                 sobolnboot=20,
                                 sobolconf=0.95,
                                 nseeds=3,
                                 precision=3)

nl@simdesign <- simdesign_sobol2007(nl=nl,
                                     samples=200,
                                     sobolnboot=20,
                                     sobolconf=0.95,
                                     nseeds=3,
                                     precision=3)

nl@simdesign <- simdesign_soboljansen(nl=nl,
                                       samples=200,
                                       sobolnboot=20,
                                       sobolconf=0.95,
                                       nseeds=3,
                                       precision=3)


nl@simdesign <- simdesign_morris(nl=nl,
                                  morristype="oat",
                                  morrislevels=4,
                                  morrisr=100,
                                  morrisgridjump=2,
                                  nseeds=3)

nl@simdesign <- simdesign_eFast(nl=nl,
                                 samples=100,
                                 nseeds=3)
```

## Optimization techniques (simdesign_GenSA, _GenAlg)

Optimization techniques are a powerful tool to search the parameter space for specific solutions.
Both approaches try to minimize a specified model output reporter by systematically (genetic algorithm, utilizing the [genalg](https://CRAN.R-project.org/package=genalg/index.html) package) or randomly (simulated annealing, utilizing the [genSA](https://CRAN.R-project.org/package=GenSA/index.html) package) changing the model parameters within the allowed ranges.
To setup optimization simdesigns, variable ranges (min, max) need to be defined.
Optimization simdesigns can only be executed using the `run_nl_dyn()` function - instead of `run_nl_all()` or `run_nl_one()`.

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath="C:/out/",
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            idfinal=NA_character_,
                            idrunnum=NA_character_,
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(min=50, max=150),
                                             'initial-number-wolves' = list(min=50, max=150)),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

nl@simdesign <- simdesign_GenAlg(nl=nl, 
                                 popSize = 200, 
                                 iters = 100, 
                                 evalcrit = 1,
                                 elitism = NA, 
                                 mutationChance = NA, 
                                 nseeds = 1)

nl@simdesign <- simdesign_GenSA(nl=nl,
                                par=NULL,
                                evalcrit=1,
                                control=list(max.time = 600),
                                nseeds=1)

```

## Calibration simdesigns (Approximate Bayesian Computation (ABC))

Approximate bayesian computation (ABC) algorithms have been increasingly used for calibration of agent-based simulation models.
The nlrx package provides different algorithms from the [EasyABC package](https://cran.r-project.org/package=EasyABC).
These algorithms can be used by attaching the corresponding simdesigns (`simdesign_ABCmcmc_Marjoram()`, `simdesign_ABCmcmc_Marjoram_original()`, `simdesign_ABCmcmc_Wegmann()`). 

```{r eval=FALSE}
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath="C:/out/",
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            idfinal=NA_character_,
                            idrunnum=NA_character_,
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            variables = list('initial-number-sheep' = list(min=50, max=150),
                                             'initial-number-wolves' = list(min=50, max=150)),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

nl@simdesign <- simdesign_ABCmcmc_Marjoram(nl=nl,
                                           summary_stat_target = c(100, 80),
                                           n_rec = 100, 
                                           n_calibration=200,
                                           use_seed = TRUE,
                                           progress_bar = TRUE,
                                           nseeds = 1)

nl@simdesign <- simdesign_ABCmcmc_Marjoram_original(nl=nl,
                                           summary_stat_target = c(100, 80),
                                           n_rec = 10, 
                                           use_seed = TRUE,
                                           progress_bar = TRUE,
                                           nseeds = 1)

nl@simdesign <- simdesign_ABCmcmc_Wegmann(nl=nl,
                                          summary_stat_target = c(100, 80),
                                          n_rec = 10, 
                                          n_calibration=200,
                                          use_seed = TRUE,
                                          progress_bar = TRUE,
                                          nseeds = 1)

```
---
title: "Spatial Output"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
editor_options: 
  chunk_output_type: console
vignette: >
  %\VignetteIndexEntry{Spatial Output}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Gathering spatial output from NetLogo model simulations

nlrx is able to gather spatial output from your NetLogo simulations.
The experiment class object provides slots for measuring turtles, patches and link variables:

* metrics.turtles - a list with named vectors of strings defining valid turtles-own variables that are taken as output measurements (e.g. `list("turtles" = c("who", "pxcor", "pycor", "color")`)
* metrics.patches - vector of valid patches-own variables that are used to collect patches data (e.g. `c("pxcor", "pycor", "pcolor")`)
* metrics.links - a list with named vectors of strings defining valid links-own variables that are taken as output measurements (e.g. `list("links" = c("[who] of end1", "[who] of end2"))`)

Basically, you can enter any variable of your model that is listed in turtles-own, patches-own or links-own, however if you add variables that contain strings, these strings must not contain any whitespaces or the output data will not be parsed correctly.

If your model has agent variables that only exist for one specific breed (breed-own), measuring those variables for all turtles would result in a runtime error. Thus, different vectors of metrics can be provided for each specific breed (e.g. `metrics.turtles = list("breed_A" = c("who", "pxcor", "pycor", "var_of_breed_A"), "breed_B" = c("who", "pxcor", "pycor", "var_of_breed_B"))`). This also works for links breeds.

If an experiment contains any agent metrics, the results of these metrics will be nested inside the output results tibble. In order to utilize these data, they need to be postprocessed. nlrx provides two types of postprocessing functions:

* `unnest_simoutput()` - This function may be used to unnest the data within the output results tibble. It reports a long tibble containing all parameters, variables and agent metrics. This format is easily subsettable and is suited to produce plots using the ggplot package, or create animations with gganimate.
* `nl_to_graph()`, `nl_to_raster()`, `nl_to_points()` - These functions create spatial objects from the nested data within the output results tibble. Please note, that these functions have requirements on measured spatial variables:
    * `nl_to_graph()` - Reports igraph objects and needs at least turtle who numbers and who numbers of link ends (`metrics.turtles = list("turtles" = c("who")), metrics.links = list("links" = c("[who] of end1", "[who] of end2"))`). Additional turtle and link variables will be stored as properties of the igraph nodes and vertices.
    * `nl_to_raster()` - Reports raster objects and needs patch coordinates and at least one patch variable (`metrics.patches = c("pxcor", "pycor", "pcolor")`). If several patch variables are provided, a raster stack is created with rasters for each patch variable.
    * `nl_to_points()` - Reports spatial point objects and needs at least turtle coordinates, either pxcor/pycor or xcor/ycor (`metrics.turtles = list("turtles" = c("xcor", "ycor", "who", "color"))`). Additional turtle variables will be stored as properties of the spatial points.


## Application Example 1: Wolf Sheep Predation

We use the Wolf Sheep Model from the NetLogo models library to capture metrics of patches and turtles. We measure coordinates, who numbers and the breed of turtles, and coordinates of patches and the corresponding pcolor on each tick. We define our experiment accordingly and run the simulations:

```{r eval=FALSE}
library(nlrx)
library(ggplot2)
library(gganimate) # devtools::install_github('thomasp85/gganimate') - if you have troubles installing gganimate, you most likely also need to install gifski as system dependency

# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

# Define nl object
nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

# Define experiment
nl@experiment <- experiment(expname = "nlrx_spatial",
                            outpath = outpath,
                            repetition = 1,      
                            tickmetrics = "true",
                            idsetup = "setup",   
                            idgo = "go",         
                            runtime = 100,
                            metrics = c("count sheep","count wolves"),
                            metrics.turtles = list("turtles" = c("who", "pxcor", "pycor")),
                            metrics.patches = c("pxcor", "pycor", "pcolor"),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             'initial-number-sheep' = 100,
                                             'initial-number-wolves' = 50,
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false")
                            )

# Attach simdesign simple using only constants
nl@simdesign <- simdesign_simple(nl=nl,
                                 nseeds=1)

# Run simulations and store output in results
results <- run_nl_all(nl = nl)
```

This experiment will run for 100 ticks (runtime) and collects all metrics, metrics.turtles and metrics.patches on each tick (evalticks).
Thus, executing `run_nl_all()` will report a tibble containing all metrics, metrics.turtles and metrics.patches. However, because the spatial metrics contain more than one value, these datasets are stored as lists inside the output tibble.
These lists already contain all measured agent metrics and can for example be used to analyze distributions of these variables for specific agent groups.


#### Postprocessing example 1: `unnest_simoutput()`

We use the `unnest_simoutput()` function to create a large tibble format which we can use for plotting. In order to use this function, the simulation output need to be attached to the nl object first (otherwise a warning will appear). The simdesign class within the nl object provides a slot for attaching output results (simoutput). An output results tibble can be attached to this slot by using the simdesign setter function `setsim(nl, "simoutput")`.

```{r eval=FALSE}
# Attach results to nl object:
setsim(nl, "simoutput") <- results

# Report spatial data:
results_unnest <- unnest_simoutput(nl)
``` 

The spatial tibble output from `unnest_simoutput()` can for example be used to plot maps for different ticks of the model simulation.
Here is an example to create a facet plot using spatial simulation data of every 10th simulation tick:

```{r eval=FALSE}

# Split tibble into turtles and patches tibbles and select each 10th step:
results_unnest_turtles <- results_unnest %>% 
  dplyr::filter(agent=="turtles") %>% 
  dplyr::filter(`[step]` %in% seq(10,80,10))
results_unnest_patches <- results_unnest %>% 
  dplyr::filter(agent=="patches") %>% 
  dplyr::filter(`[step]` %in% seq(10,80,10))

# Create facet plot:
ggplot() +
  facet_wrap(~`[step]`, ncol=4) +
  coord_equal() +
  geom_tile(data=results_unnest_patches, aes(x=pxcor, y=pycor, fill=factor(pcolor))) +
  geom_point(data=results_unnest_turtles, aes(x = pxcor, y = pycor, color = breed), size=1) +
  scale_fill_manual(breaks=c("35", "55"), values = c("35" = "#D9AF6B", "55" = "#68855C")) +
  scale_color_manual(breaks=c("sheep", "wolves"), values = c("sheep" = "beige", "wolves" = "black")) +
  guides(fill=guide_legend(title="LandCover")) +
  theme_minimal() +
  ggtitle("Output maps of each 10th simulation tick")

``` 

<img src="wolfsheep_world.png" align="center" width="100%" />

Using the gganimate package(https://github.com/thomasp85/gganimate), it is even possible to generate animated plots from this spatial data tibble.
Here is an example for a plot that has been generated by running the above experiment and postprocessing the data with `unnest_simoutput()`.


```{r eval=FALSE}

# Split tibble into turtles and patches tibbles:
results_unnest_turtles <- results_unnest %>% 
  dplyr::filter(agent == "turtles")
results_unnest_patches <- results_unnest %>% 
  dplyr::filter(agent == "patches")

# Create an animated plot, using the step column as animation variable
p1 <- ggplot() +
  geom_tile(data=results_unnest_patches, aes(x=pxcor, y=pycor, fill=factor(pcolor))) +
  geom_point(data=results_unnest_turtles, aes(x = pxcor, y = pycor, group=who, color = breed), size=2) +
  scale_fill_manual(breaks=c("35", "55"), values = c("35" = "#D9AF6B", "55" = "#68855C")) +
  scale_color_manual(breaks=c("sheep", "wolves"), values = c("sheep" = "beige", "wolves" = "black")) +
  guides(fill=guide_legend(title="LandCover")) +
  transition_time(`[step]`) +
  coord_equal() +
  labs(title = 'Step: {frame_time}') +
  theme_void()

# Animate the plot and use 1 frame for each step of the model simulations
gganimate::animate(p1, nframes = length(unique(results_unnest_patches$`[step]`)), width=400, height=400, fps=4)
anim_save("wolfsheep_world.gif")

``` 

<center>
<img src="wolfsheep_world.gif" align="center" width="50%" />
</center>

#### Postprocessing example 2: `nl_to_raster()` and `nl_to_points()`

A second option to postprocess spatial data is to convert patches data to raster objects and turtles data to spatial point objects. This can be handy if further spatial analyses are carried out, such as calculating landscape metrics (e.g. with R-package [landscapemetrics](https://github.com/r-spatialecology/landscapemetrics)).
In this example, we use the nl object with attached simulation output of the previous example to generate raster and spatial point objects.

```{r eval=FALSE}
## Create raster and point objects from patches and turtles data:
library(raster)
library(sf)
nlraster <- nl_to_raster(nl)
nlpoints <- nl_to_points(nl, coords = "px")

## Plot raster and turtles of tick n:
n <- 1
plot(nlraster$spatial.raster[[n]], col=c("35" = "#D9AF6B", "55" = "#68855C"))
plot(nlpoints$spatial.turtles[[n]]["breed"], add=TRUE, pch=16, col=c("sheep" = "beige", "wolves" = "black"))

``` 


## Application example 2: Giant component

We use the Giant Component model from the NetLogo models library to capture metrics of turtles and links. 
The `nl_to_graph()` function generates igraph objects from measured turtles and links data.
It produces an igraph object for each row of the simoutput results tibble. Thus, it reports one igraph network for each combination of random-seed, siminputrow and step.

In order to generate igraph objects some metrics are mandatory:
The metrics.turtles slot of the experiment must contain "who" numbers (see example experiment).
Additional turtle metrics will be stored as properties of the igraph vertices.
The metrics.links slot of the experiment must contain "who" numbers of link end1 and end2 (see example experiment).
Additional link metrics will be stored as properties of the igraph edges.

For our application example, we are mainly interested in the final network structure. We want to measure who numbers and color of turtles, and who numbers of end1 and end2 of links on the final tick. We define our experiment accordingly and run the simulations:

```{r eval=FALSE}

library(nlrx)
library(igraph)

# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.4")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Networks/Giant Component.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.4")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Networks/Giant Component.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.4",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

nl@experiment <- experiment(expname="networks",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="false",
                            idsetup="setup",
                            idgo="go",
                            runtime=50,
                            metrics.turtles = list("turtles" = c("who", "color")),
                            metrics.links = list("links" = c("[who] of end1", "[who] of end2")),
                            constants = list("num-nodes" = 80,
                                             "layout?" = "true"))

nl@simdesign <- simdesign_simple(nl, 1)
results <- run_nl_all(nl)

```

In order to execute the `nl_to_graph()` function, the output results need to be attached to the nl object.
Afterwards we can create the igraph object with `nl_to_graph()`:

```{r eval=FALSE}
# Attach results to nl object:
setsim(nl, "simoutput") <- results
# Create igraph:
nl.graph <- nl_to_graph(nl)
``` 

The igraph objects are attached to the output results tibble in the newly generated column "spatial.links".
In our case, we only have one row of results and thus only one igraph object.
We can now extract this object and do some plotting, using the igraph plotting function:

```{r eval=FALSE}
## Extract graph of tick 1:
nl.graph.i <- nl.graph$spatial.links[[1]]
## Set vertex colors by measured color variable:
vcols <- c("7" = "grey", "15" = "red")
V(nl.graph.i)$color <- vcols[as.character(V(nl.graph.i)$color)]
## Set edge colors by measured link breed:
ecols <- c("links" = "black")
E(nl.graph.i)$color <- ecols[E(nl.graph.i)$breed]

## Plot:
plot.igraph(nl.graph.i, vertex.size=8, vertex.label=NA, edge.arrow.size=0.2)
``` 

We can also calculate network metrics using the igraph package:

```{r eval=FALSE}
## Extract graph of tick 1:
nl.graph.i <- nl.graph$spatial.links[[1]]

## Vertex and edge betweenness centrality
betweenness(nl.graph.i)

## Clusters
clusters(nl.graph.i)

``` 
---
title: "Capturing NetLogo output manually"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
editor_options: 
  chunk_output_type: console
vignette: >
  %\VignetteIndexEntry{Capturing NetLogo output manually}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Capturing output manually

While nlrx provides the `metrics`, `metrics.turtles`, `metrics.patches` and `metrics.links` slots of the `exerpiment` class to capture global and agent related output from NetLogo models, this might not be sufficient under certain circumstances. For example, the `metrics` slot cannot track and collect nested list output or similar complex data structures such as NetLogo arrays or matrices. Another example would be spatial output such as shapefiles or grids, generated with the GIS extension. Additionally, you might already have implemented complex routines in your model which write model output to disk.

The main question for this vignette to answer is: How can you link such self-written output to your nlrx experiment. Here we show you one basic example how this works. In this example we write ascii raster files using the GIS extension in NetLogo. However, the same workflow can of course be applied to other types of output, for example text files written with `file-type` primitives or the csv extension.


#### Step 1: Creating an ID widget on your model interface

In order to link self-written output to our nlrx simulations within the R session, we need to transfer the current nlrx siminputrow and seed to our NetLogo model. That way, when executing `run_nl_all()`, within each simulation the NetLogo model exactly "knows" which parameter row (siminputrow) is currently simulated and the correpsonding random seed. To transfer this information we just need to create a string input widget on our model interface. In the example below, we created such a widget at the bottom of the Wolf Sheep model called `nlrx_id`.

<center>
<img src="idrunnum_gui.png" align="center" width="75%" />
</center>

#### Step 2: Define idrunnum of the experiment

The idrunnum field of the experiment is a built-in feature that allows to transfer the current experiment name, siminputrow and random seed to a defined NetLogo gui parameter widget.
In our case we want to use our newly created `nlrx_id` string input field, so we just set `idrunnum = nlrx_id`.

```{r eval=FALSE}
# Attach experiment
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            idrunnum="nlrx_id",
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("count sheep"),
                            variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                             'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))
```


#### Step 3: Use nlrx_id within NetLogo model to tag self-written output

We can now use the `nlrx_id` field within the NetLogo model to tag our self-written output.
One way to do this is to tag the filenames of self-written output. Here is an example where raster grids are created and written at the end of each model run:

```{netlogo eval=FALSE}
to write_final_landscape
  let filename (word nlrx_id "_" but-first (word (1000 + ticks)) ".asc")
  if (file-exists? filename) [file-delete filename]
  gis:store-dataset final_landscape filename  
end
```

Here, we first create our filename by using the string from the `nlrx_id` field (which contains the current experiment name, siminputrow and random seed divided by an underscore) and add the current tick (containing leading zeros by adding 1000 to the tick count and removing the first number afterwards).

Of course you could also use different approaches to utilize the information within the `nlrx_id` field for tagging your output.


#### Step 4: Linking self-written output to nlrx collected output

Finally, we want to read in our self-written raster files, calculate some landscape metrics and add those to the results table that we received from the `run_nl_all()` function.

We basically have two types of output now: The tibble from our model executions, and a folder with self-written output (here defined as a subfolder of our modelpath).

```{r eval=FALSE}
## Output from nlrx simulations:
results <- run_nl_all(nl)
## Self-writte output directory:
ascdir <- file.path(dirname(nl@modelpath), "output")
``` 

We can now loop over the files in that folder, read the content and use the `strplit` function on the filename to identify the experiment name, the siminputrow, the random seed and the tick count. As an example application, we then calculate some landscapemetrics using the `landscapemetrics` package. Finally we put everything into a results tibble:

```{r eval=FALSE}
results.lsm <- purrr::map_dfr(list.files(ascdir, pattern = "asc", full.names = TRUE), function(x) {
  x.split <- strsplit(x, "_")[[1]]
  x.tick <- as.numeric(strsplit(x.split[[length(x.split)]], "\\.")[[1]][[1]])
  x.siminputrow <- as.numeric(x.split[[length(x.split) - 1]][[1]])
  x.seed <- as.numeric(x.split[[length(x.split) - 2]][[1]])
  x.raster <- raster(x)
  
  ## netlogo asc files use NaN as default nodata value in the asc file header
  ## this leads to problems when reading the raster because it sets ´zeros to NA
  ## here we set NAs back to zertos manually:
  x.raster <- reclassify(x.raster, cbind(NA, 0))
  
  ## Calculate landscape metrics:
  metrics <- c("lsm_l_ed", "lsm_l_shdi", "lsm_l_lsi", "lsm_l_lpi", "lsm_l_area_mn")
  x.metrics <- landscapemetrics::calculate_lsm(x.raster, what=metrics) %>% 
    dplyr::select(metric, value) %>% 
    tidyr::pivot_wider(names_from=metric, values_from = value)
  
  x.final <- tibble::tibble(siminputrow = x.siminputrow,
                            `[step]` = x.tick,
                            `random-seed` = x.seed)
  
  x.final <- cbind(x.final, x.metrics)
  return(x.final)
})
```

Now we have the same identifier columns in the nlrx reported output and the self-written output tibble which allows us to join both tibbles together:

```{r eval=FALSE}
## Combine results with lsm and store:
results <- results %>% left_join(results.lsm, by = c("siminputrow", "random-seed", "[step]"))

## Attach output to nl object:
setsim(nl, "simoutput") <- results
saveRDS(nl, file = file.path(outpath, "my_final_nl_object.rds"))
```
---
title: "Publication record"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Publication record}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Publication record

On this page we aim to collect research papers that have used and cited the nlrx package. We will update this database from time to time. If your article is missing, please drop us an [email](mailto:jan.salecker@posteo.de) or file an [issue](https://github.com/ropensci/nlrx/issues/new).

## Share your nlrx use cases

We are always interested in all kind of use cases where nlrx is applied. We therefore encourage you to share your tips and tricks for nlrx so that the nlrx community can benefit from your experience. A nice way to share your use cases is the [ropensci community board](https://discuss.ropensci.org/c/usecases/10). Your use case will then also be featured via the ropensci twitter account. Of course you can also get in touch with us if you want to share your work, e.g. on our documentation and examples pages.

## Using nlrx in your research

If you use the nlrx package in your research, please cite the nlrx package as:

Salecker J, Sciaini M, Meyer KM, Wiegand K. The nlrx r package: A next-generation framework for reproducible NetLogo model analyses. Methods Ecol Evol. 2019;2041-210X. [https://doi.org/10.1111/2041-210X.13286](https://doi.org/10.1111/2041-210X.13286).

Citation information for `nlrx` can also be accessed in R doing `citation(package = 'nlrx')`.


## Publication record

<script src="//bibbase.org/show?bib=https://raw.githubusercontent.com/ropensci/nlrx/master/vignettes/articles/nlrx_papers.bib&jsonp=1&sort=author_short"></script>

---
title: "Interfacing a variety of different NetLogo models"
author: "Marco Sciaini"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Interfacing a variety of different NetLogo models}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Here, we present examples on how to interface a variety of different NetLogo models.
NetLogo model interfaces are very diverse across models and might contain spatial elements of different types (patches, different groups of turtles, links). With our interfacing examples, we try to cover a wide range of different models to show how **nlrx** can be utilized for various applications.

Our examples show how spatial output from NetLogo models can be used to create animated gif files, displaying patches and turtles properties over time. We explicitly show how to setup, postprocess and plot spatial output for very different model types. Although, the following examples cover a wide range of output types, it is likely that different experiment measurements or postprocessing needs to be done for specific and complex models. If you have problems to interface a specific type of model, feel free to send one of the authors of nlrx an email!

First, we need to load a bunch of packages, mainly for the subsetting of the spatial data
and the plotting:

```{r eval = FALSE}
## Load packages
library(nlrx)
library(tidyverse)
library(ggplot2)
library(gganimate)
library(cartography) 
library(rcartocolor)
library(ggthemes) 
```

Second, we need to define the path variables pointing to the NetLogo installation.
Depending on the operation system, the NetLogo and output paths need to be adjusted.

```{r eval = FALSE}
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.4")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.4")
outpath <- file.path("/home/out")
```

... and now are ready to interface some of the models of the NetLogo model library. 



## Diffusion model

The diffusion model is a simple example of how to plot moving turtles and patches.
Turtles move over the lattice and generate heat on cells they move over. The accumulated heat also diffuses to neighboring cells.

In order to transfer the spatial information from NetLogo to R, we measure the who number and coordinates of all turtles (metrics.turtles) and coordinates, patch color and amount of heat of all patches (metrics.patches).

After converting the spatial data to tibble format (`unnest_simoutput()`), we split the spatial data into a turtles tibble and a patches tibble using the agent column of the returned tibble from `unnest_simoutput()`.
The plot consist of a geom_tile layer, which displays the amount of heat in each cell and a geom_point layer that displays the movement of the turtles. In order to render animated gifs with gganimate it is very important to set the group property of the turtles layer correctly (group=who). 

```{r  eval = FALSE}
## Step1: Create a nl obejct:
nl <- nl(nlversion = "6.0.4",
         nlpath = netlogopath,
         modelpath = file.path(netlogopath, "app/models/Sample Models/Art/Diffusion Graphics.nlogo"),
         jvmmem = 1024)

## Step2: Add Experiment
nl@experiment <- experiment(expname = "diffusion",
                            outpath = outpath,
                            repetition = 1,   
                            tickmetrics = "true",
                            idsetup = "setup",  
                            idgo = "go",        
                            runtime = 500,
                            evalticks = seq(1,500),
                            constants = list("num-turtles" = 20,
                                             "diffusion-rate" = 1,
                                             "turtle-heat" = 139,
                                             "turtle-speed" = 1.0,
                                             "wander?" = "TRUE"),
                            metrics.turtles = list("turtles" = c("who", "xcor", "ycor")),
                            metrics.patches = c("pxcor", "pycor", "pcolor", "heat"))

# Evaluate if variables and constants are valid:
eval_variables_constants(nl)

## Step3: Add a Simulation Design
nl@simdesign <- simdesign_simple(nl = nl,
                                 nseeds = 1)


# Step4: Run simulations:
results <- run_nl_all(nl = nl)

## Postprocessing:
## Step1: Attach results to nl:
setsim(nl, "simoutput") <- results


# Prepare data for plotting
nl_spatial <- unnest_simoutput(nl)
n <- nl@experiment@runtime

turtles <- nl_spatial %>% dplyr::filter(agent == "turtles" & `[step]` < n) %>% dplyr::select(xcor, ycor, `[step]`, who)
patches <- nl_spatial %>% dplyr::filter(agent == "patches"& `[step]` < n) %>% dplyr::select(pxcor, pycor, pcolor, `[step]`)

## Plot animation:
p1 <- ggplot() +
  geom_tile(data = patches, aes(x=pxcor, y = pycor, fill=pcolor)) +
  geom_point(data = turtles, aes(x = xcor, y = ycor, group = who), size=2) +
  scale_fill_carto_c(palette = "Prism") +
  transition_time(`[step]`) +
  guides(fill="none", color="none") +
  coord_equal() +
  theme_void() 

gganimate::animate(p1, nframes = n, width=800, height=800, duration=6)
```

<center>
<img src="https://github.com/nldoc/nlrx_playground/raw/master/gif/diffusion.gif" align="center" width="100%" />
</center>

## Fire model

Similar to the diffusion example, the fire model displays patches and moving turtles.
However, the model contains two turtle breeds which need to be displayed differently.

Thus, we also measure the turtle breed (metrics.turtles).
Before plotting we split up the turtle data into two tibbles, one for each breed. This allows us to create two different geom_point layers with different aesthetics.

```{r  eval = FALSE}
## Step1: Create a nl obejct:
nl <- nl(nlversion = "6.0.4",
         nlpath = netlogopath,
         modelpath = file.path(netlogopath, "app/models/Sample Models/Earth Science/Fire.nlogo"),
         jvmmem = 1024)

## Step2: Add Experiment
nl@experiment <- experiment(expname = "fire",
                            outpath = outpath,
                            repetition = 1,    
                            tickmetrics = "true",
                            idsetup = "setup", 
                            idgo = "go",       
                            runtime = 500,
                            evalticks = seq(1,500),
                            metrics = c("count patches"),
                            metrics.turtles = list("turtles" = c("who", "pxcor", "pycor", "breed", "color")),
                            metrics.patches = c("pxcor", "pycor", "pcolor"),
                            constants = list('density' = 62)
)

# Evaluate if variables and constants are valid:
eval_variables_constants(nl)

## Step3: Add a Simulation Design
nl@simdesign <- simdesign_simple(nl = nl,
                                 nseeds = 1)


# Step4: Run simulations:
results <- run_nl_all(nl = nl)

## Postprocessing:
## Step1: Attach results to nl:
setsim(nl, "simoutput") <- results

nl_spatial <- unnest_simoutput(nl)
n <- max(nl_spatial$`[step]`)
         
embers <- nl_spatial %>% dplyr::filter(breed == "embers" & agent == "turtles" & `[step]` < n) %>% dplyr::select(pxcor, pycor, `[step]`, color, who)
fires <- nl_spatial %>% dplyr::filter(breed == "fires" & agent == "turtles" & `[step]` < n) %>% dplyr::select(pxcor, pycor, `[step]`, color, who)
patches <- nl_spatial %>% dplyr::filter(agent == "patches" & `[step]` < n) %>% dplyr::select(pxcor, pycor, pcolor, `[step]`)

# make same space in the memory
rm(nl)
rm(results)
rm(nl_spatial)
gc()

#----------------------------------
## Plot animation:
p1 <- ggplot(embers) +
  geom_tile(data = patches, aes(x=pxcor, y=pycor, fill=factor(pcolor))) +
  geom_point(data = embers, aes(x = pxcor, y = pycor, color = color, group = who), size=2) +
  scale_color_gradientn(colors = rev(cartography::carto.pal("orange.pal", n1 = 8))) +
  geom_point(data = fires, aes(x = pxcor, y = pycor, color = color, group = who), size=3) +
  scale_fill_manual(values = c("0" = "gray24", "55" = "#83B37D", "11.4" = "#B59A89")) +
  transition_time(`[step]`) +
  guides(fill="none", color="none") +
  coord_equal() +
  theme_void() 
  
g1 <- gganimate::animate(p1, nframes = n, width=800, height=800, duration = 6)
```

<center>
<img src="https://github.com/nldoc/nlrx_playground/raw/master/gif/fire.gif" align="center" width="100%" />
</center>



## Ants model

The Ants model is another typical example for displaying patches and moving turtles. 
Ants move across the lattice in order to gather food from food sources and returning food to their nest in the center.
Ants carrying food leave a trail of chemicals on cells they move over. These chemicals also diffuse to neighboring cells over time. Displaying the patches output from the model is not straightforward because three different groups of patches exist: food sources, ants nest and all other patches which are coloured depending on the amount of chemicals that the patch contains.

In order to transfer the spatial information from NetLogo to R, we measure the who number, coordinates and breed for all turtles (metrics.turtles) and coordinates, patch color, amount of chemicals, presence of food and the food source number of all patches (metrics.patches).

After converting the spatial data to tibble format (`unnest_simoutput()`), we split the patches information in three objects in order to display the three groups of patches:

* chem_dat: Contains chemical distribution of all patches. The amount of chemicals of each patch is plotted by using a geom_tile layer that acts as base layer of our plot.
* food_dat: Contains patches coordinates of patches that contain food sources (food > 0). The food sources are plotted by adding a square shaped geom_point layer on top of the chemicals distribution layer.
* nest_dat: Contains patches coordinates of the nest patch. The nest is plotted by adding a large circle geom_point layer on top.

The turtle data is plotted by adding another geom_text layer, using a Yen sign (￥) that displays the movement of the turtles. In order to render animated gifs with gganimate it is very important to set the group property of the layer correctly (group=who).

The plot need to be stored in an object before it can be rendered with the `gganimate` function.
In some cases, it is important to set the number of frames explicitly to the number of steps of your simulation output.
Otherwise, gganimate may interpolate between measured steps, which may lead to incorrect results.

```{r eval = FALSE}
## Step1: Create a nl obejct:
nl <- nl(nlversion = "6.0.4",
         nlpath = netlogopath,
         modelpath = file.path(netlogopath, "app/models/Sample Models/Biology/Ants.nlogo"),
         jvmmem = 1024)

## Step2: Add Experiment
nl@experiment <- experiment(expname = "ants",
                            outpath = outpath,
                            repetition = 1,      
                            tickmetrics = "true",
                            idsetup = "setup",   
                            idgo = "go",         
                            runtime = 1000,
                            evalticks = seq(1,1000),
                            metrics.turtles = list("turtles" = c("who", "pxcor", "pycor", "breed")),
                            metrics.patches = c("pxcor", "pycor", "pcolor", "chemical", "food", "food-source-number"),
                            constants = list("population" = 125,
                                             'diffusion-rate' = 50,
                                             'evaporation-rate' = 10))

## Step3: Add a Simulation Design
nl@simdesign <- simdesign_simple(nl = nl,
                                 nseeds = 1)


## Step4: Run simulations:
results <- run_nl_all(nl = nl)

## Step5: Attach results to nl and reformat spatial data with get_nl_spatial()
setsim(nl, "simoutput") <- results
nl_spatial <- unnest_simoutput(nl)

## Step6: Prepare data for plotting
# Extract infromation on food sources and select maximum step as simulation length:
food_dat <- nl_spatial %>% dplyr::filter(food > 0) %>% dplyr::select(pxcor, pycor, `[step]`)
nmax <- max(food_dat$`[step]`)
food_dat <- food_dat %>% dplyr::filter(`[step]` %in% 1:nmax)
# Extract information on chemicals and apply minimum treshhold for coloring:
chem_dat <- nl_spatial %>% dplyr::filter(`[step]` %in% 1:nmax) %>% dplyr::select(pxcor, pycor, chemical, `[step]`)
chem_dat$chemical <- ifelse(chem_dat$chemical <= 0.2, NA, chem_dat$chemical)
# Extract information on turtle positions:
turt_dat <- nl_spatial %>% dplyr::filter(!is.na(who)) %>% dplyr::filter(`[step]` %in% 1:nmax) %>% dplyr::select(pxcor, pycor, who, `[step]`)
# Create a new data frame to overlay the nest position (in this case the center of the world 0,0)
nest_dat <- food_dat
nest_dat$pxcor <- 0
nest_dat$pycor <- 0

## Step7: Plotting
p1 <- ggplot(food_dat) +
  geom_tile(data=chem_dat, aes(x=pxcor, y=pycor, fill=sqrt(chemical))) +
  geom_point(data=food_dat, aes(x=pxcor, y=pycor), color="black", shape=15, size=4.5, alpha=1) +
  geom_text(data=turt_dat, aes(x=pxcor, y=pycor, group=who, color=as.numeric(who)), size=5, alpha=1, label="￥") +
  geom_point(data=nest_dat, aes(x=pxcor, y=pycor), color="brown", fill="white", size=30, stroke=2, shape=21) +
  scale_fill_viridis_c(direction=-1, option="magma", na.value = "white") +
  scale_color_gradient_tableau(palette="Orange") +
  transition_time(`[step]`) +
  guides(fill="none", color="none") +
  coord_equal() +
  labs(title = 'Step: {frame_time}') +
  theme_void()

## Step8: Animate the plot and use 1 frame for each step of the model simulations
gganimate::animate(p1, nframes = nmax, width=800, height=800, fps=10)
```

<center>
<img src="https://github.com/nldoc/nlrx_playground/raw/master/gif/ants.gif" align="center" width="100%" />
</center>



## Waves model

The waves model contains no patch information. Three types of turtle agents are present in the model: waves, edges and a driver. However, these agents are not defined by using breeds but by using turtles-own variables (driver?, edge?).

In addition to these two variables we measure the who number, turtle coordinates and turtle color.
Again, before plotting we split the spatial data tibble from `unnest_simoutput()` into three different tibbles, one for each agent group. These agent groups are then plotted by using geom_point layers with different aesthetics.

Instead of a simple simdesign with only one parameterization, this example uses a distinct simdesign with 4 different friction values. We can use the facet functionality of ggplot to visualize all four model runs at once.

```{r  eval = FALSE}

## The models library filepath pointing to the waves model contains an ampersand which needs to be escaped on Windows but not on Linux.
## Choose the filepath according to your OS:
## Windows:
modelpath <- "app/models/Sample Models/Chemistry ^& Physics/Waves/Wave Machine.nlogo"
## Linux:
modelpath <- "app/models/Sample Models/Chemistry & Physics/Waves/Wave Machine.nlogo"

## Step1: Create a nl obejct:
nl <- nl(nlversion = "6.0.4",
         nlpath = netlogopath,
         modelpath = file.path(netlogopath, modelpath),
         jvmmem = 1024)

## Step2: Add Experiment
nl@experiment <- experiment(expname = "waves",
                            outpath = outpath,
                            repetition = 1,      
                            tickmetrics = "true",
                            idsetup = "setup",  
                            idgo = "go",
                            runtime = 100,
                            evalticks = seq(1,100),
                            variables = list('friction' = list(values = c(5,25,50,90))),
                            constants = list("stiffness" = 20),
                            metrics.turtles = list("turtles" = c("who", "xcor", "ycor", "driver?", "edge?", "color")))

## Step3: Add a Simulation Design
nl@simdesign <- simdesign_distinct(nl = nl,
                                   nseeds = 1)


# Step4: Run simulations:
results <- run_nl_all(nl = nl)

## Postprocessing:
## Step1: Attach results to nl:
setsim(nl, "simoutput") <- results


# Prepare data for plotting
nl_spatial <- unnest_simoutput(nl)

nl_spatial$friction <- factor(nl_spatial$friction)
levels(nl_spatial$friction) <- c("Friction = 5",
                                 "Friction = 25",
                                 "Friction = 50",
                                 "Friction = 90")

n <- nl@experiment@runtime

waves  <- nl_spatial %>% dplyr::filter(`driver?` == "false" & `edge?` == "false" & `[step]` < n) %>% dplyr::select(xcor, ycor, `[step]`, who, color, friction)
egde   <- nl_spatial %>% dplyr::filter(`driver?` == "false" & `edge?` == "true" & `[step]` < n) %>% dplyr::select(xcor, ycor, `[step]`, who, friction)
driver <- nl_spatial %>% dplyr::filter(`driver?` == "true" & `edge?` == "false" & `[step]` < n) %>% dplyr::select(xcor, ycor, `[step]`, who, friction)

p1 <- ggplot(waves) +
  geom_point(data = waves, aes(x = xcor, y = ycor, group = who, color = color), size=2) +
  geom_point(data = driver, aes(x = xcor, y = ycor, group = who), size=2, color = "grey") +
  geom_point(data = egde, aes(x = xcor, y = ycor, group = who), size=2, color = "black") +
  facet_wrap(~friction) +
  transition_time(`[step]`) +
  guides(color="none") +
  coord_equal() +
  theme_void() 

gganimate::animate(p1, width=800, height=800, duration = 6)
```

<center>
<img src="https://github.com/nldoc/nlrx_playground/raw/master/gif/waves.gif" align="center" width="100%" />
</center>

## Flocking Model

The flocking model describes bird flocking formation using very simplified movement rules.
For this model, headings (=viewing angles) of agents are crucial.
The NetLogo heading system is measured in degree, ranging from 0 to 359. An agent with heading = 0 is pointing straight to the top border of the lattice.
In order to transfer the NetLogo heading system to the ggplot and gganimate angle system, some conversions need to be done. The type of conversion depends on the geom that is used for plotting.
We present two examples for heading conversion:

(a) geom_text can be used to plot arrows pointing into the heading of the agent. geom_text also uses degree angles, however an angle of 0 means text is pointing to the right edge of the plotting area. Furthermore, angles are counter-clockwise, whereas NetLogo headings are defined in clockwise direction. In our case we use a right-pointing ASCII-arrow as geom_text label. Thus, we need to shift all NetLogo headings by 90 degree to the right and invert the direction by multiplying all headings with -1.

(b) geom_spoke can be used to plot arrows pointing into the heading of the agent. geom_spoke uses radians angles. Similar to geom_text angles, an angle of 0 is pointing to the right and angles are defined in counter-clockwise direction. Thus, we need to shift NetLogo headings by 90 degree to the right, inverse the direction by multiplying all headings with -1 and finally transform the degree angles to radians angles.

```{r eval = FALSE}
## Step1: Create a nl obejct:
nl <- nl(nlversion = "6.0.4",
         nlpath = netlogopath,
         modelpath = file.path(netlogopath, "app/models/Sample Models/Biology/Flocking.nlogo"),
         jvmmem = 1024)

## Step2: Add Experiment
nl@experiment <- experiment(expname = "flocking",
                            outpath = outpath,
                            repetition = 1,      
                            tickmetrics = "true",
                            idsetup = "setup",   
                            idgo = "go",         
                            runtime = 300,
                            evalticks = seq(1,300),
                            constants = list("population" = 100,
                                             "vision" = 5,
                                             "minimum-separation" = 1, 
                                             "max-align-turn" = 4,
                                             "max-cohere-turn" = 4,
                                             "max-separate-turn" = 4),
                            metrics.turtles = list("turtles" = c("who", "xcor", "ycor", "heading", "color")))

# Evaluate if variables and constants are valid:
eval_variables_constants(nl)

## Step3: Add a Simulation Design
nl@simdesign <- simdesign_simple(nl = nl,
                                 nseeds = 1)


# Step4: Run simulations:
results <- run_nl_all(nl = nl)

setsim(nl, "simoutput") <- results

nl_spatial <- unnest_simoutput(nl)

## Calculate angles for plotting:
# (a) convert NetLogo degree headings (clockwise with 0 = top) to geom_text degree angle (counter-clockwise with 0 = right)
# (b) convert geom_text degree angle to geom_spoke radians angle
nl_spatial <- nl_spatial %>% 
  dplyr::select(`[step]`, who, xcor, ycor, heading, color) %>% 
  mutate(heading_text = ((heading * -1) + 90)) %>% 
  mutate(heading_radians = ((heading_text * pi) / 180))

# Plot with geom_point and geom_spoke
p1 <- ggplot(nl_spatial, aes(x=xcor, y=ycor)) +
  geom_point(aes(color=color, group=who)) +
  geom_spoke(aes(color=color, group=who, angle=heading_radians), arrow=arrow(length=unit(0.2, "inches")), radius=1) +
  scale_color_viridis_c() +
  guides(color=FALSE) +
  transition_time(`[step]`) +
  coord_equal() +
  theme_void() 

gganimate::animate(p1, nframes=max(nl_spatial$`[step]`), width=600, height=600, fps=20)
```

<center>
<img src="https://github.com/nldoc/nlrx_playground/raw/master/gif/flocking_spoke.gif" align="center" width="75%" />
</center>

```{r eval = FALSE}
# Plot with geom_text
p2 <- ggplot(nl_spatial, aes(x=xcor, y=ycor)) +
  geom_text(aes(color=color, group=who, angle=heading_text), label="→", size=5) +
  scale_color_viridis_c() +
  guides(color=FALSE) +
  transition_time(`[step]`) +
  coord_equal() +
  theme_void() 

gganimate::animate(p2, nframes=max(nl_spatial$`[step]`), width=600, height=600, fps=20)

```

<center>
<img src="https://github.com/nldoc/nlrx_playground/raw/master/gif/flocking_text.gif" align="center" width="75%" />
</center>

## Preferential Attachment

The preferential attachment model from the NetLogo models library is a simple model of a growing network.
New nodes spawn each time step and connect to the already existing network.
The network is realized by connecting turtles with links.
This example shows how link metrics can be measured, using the metrics.links slot of the experiment.
Link positions in NetLogo are defined by start and end points (end1, end2) using who numbers of the turtles they connect.
Thus, on each tick we measure end1 and end2 of each link (metrics.links) and the who number, xcor, ycor and color of each node turtle (metrics.turtles).

#### Visualize network using ggplot

Again, after attaching the simulation results and postprocessing with `unnest_simoutput()` we split the tibble into a turtles tibble and a links tibble by subsetting the group column.

Finally, we need to reference the who number start and end points of the links to the actual coordinates of the corresponding turtle on each step. This can be easily achieved with two left_join calls, one for end1 who numbers and one for end2 who numbers.

```{r eval = FALSE}
## Step1: Create a nl obejct:
nl <- nl(nlversion = "6.0.4",
         nlpath = netlogopath,
         modelpath = file.path(netlogopath, "app/models/Sample Models/Networks/Preferential Attachment.nlogo"),
         jvmmem = 1024)

## Step2: Add Experiment
nl@experiment <- experiment(expname = "networks",
                            outpath = outpath,
                            repetition = 1,    
                            tickmetrics = "true",
                            idsetup = "setup", 
                            idgo = "go",       
                            runtime = 200,
                            evalticks = seq(1,200),
                            constants = list("layout?" = TRUE),
                            metrics.turtles = list("turtles" = c("who", 
                                                                 "xcor",                            
                                                                 "ycor",
                                                                 "color")),
                            metrics.links = list("links" = c("[who] of end1","[who] of end2")))

## Step3: Add simdesign
nl@simdesign <- simdesign_simple(nl=nl, nseeds = 1)

## Run simulation:
results <- run_nl_one(nl = nl,
                      seed = getsim(nl, "simseeds")[1],
                      siminputrow = 1)

## Attach results to nl
setsim(nl, "simoutput") <- results

## Postprocess spatial metrics with get_nl_spatial:
nl_spatial <- unnest_simoutput(nl)

## Subset nl_spatial using the group column:
nl_links <- nl_spatial %>%
  dplyr::filter(agent == "links") %>% 
  dplyr::select(-who, -xcor, -ycor, -color)

nl_turtles <- nl_spatial %>%
  dplyr::filter(agent == "turtles") %>% 
  dplyr::select(`[step]`, who, xcor, ycor, color)

## Reference who numbers of link start (end1) and end (end2) points to actual coordinates:
nl_links <- nl_links %>% 
  dplyr::left_join(nl_turtles, by=c("[step]"="[step]","end1" = "who")) %>% 
  dplyr::left_join(nl_turtles, by=c("[step]"="[step]","end2" = "who"))

## Plot:
p1 <- ggplot() +
  geom_point(data=nl_turtles, aes(x = xcor, y = ycor, group=who), color="red", size=2) +
  geom_segment(data=nl_links, aes(x = xcor.x, y = ycor.x, xend = xcor.y, yend = ycor.y), size=0.5) +
  transition_time(`[step]`) +
  coord_equal() +
  theme_void()

gganimate::animate(p1, nframes = max(nl_turtles$`[step]`), width=400, height=400, fps=8)

```

<center>
<img src="https://github.com/nldoc/nlrx_playground/raw/master/gif/networks.gif" align="center" width="75%" />
</center>

#### Visualize network using igraph package

It is also possible to use the [igraph](https://igraph.org/r/) package to visualize a NetLogo network.
We use the `nl_to_graph()` function from the nlrx package to convert our spatial data into an igraph object.
Here we want to display the network for only one simulation step, thus we filter our spatial data and select only the last simulation step. Again, we split the data into a turtles and a links tibble and rename the columns.

After conversion to an igraph object we can use the power of network analysis functions from the igraph package to do further analyses. In this example we perform a community cluster algorithm to color the nodes based on their community membership within the network.

```{r eval = FALSE}

## Load igraph package and convert spatial tibbles to igraph:
library(igraph)
nw.all <- nl_to_graph(nl)

## Select the last step and extract list element to get igraph object
nw <- nw.all %>% dplyr::filter(`[step]` == max(`[step]`))
nw <- nw$spatial.links[[1]]

## Perform community cluster algorithm and set community dependend colors:
com <- walktrap.community(nw)
V(nw)$community <- com$membership
rain <- rainbow(14, alpha=.5)
V(nw)$color <- rain[V(nw)$community]

## Plot:
plot(nw,
     vertex.label=NA,
     vertex.size=4,
     edge.curved=0,
     layout=layout_with_fr(nw, niter = 3000))

```

<center>
<img src="networks_igraph.png" align="center" width="75%" />
</center>
---
title: "NetLogo documentation with nlrx"
author: "Jan Salecker"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{NetLogo documentation with nlrx}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# nldoc

The `nldoc` function within nlrx can be used to create NetLogo documentations using R.
The function searches NetLogo model files for Markdown headers.
The documentation can be created in different formats (html, pdf, docx) and styles.
There are also utility functions, such as creation of a procedure network (`nldoc_network()`).

## Example

#### Model documentation

In order to create NetLogo model documentations with nlrx, documentation tags need to be added to the NetLogo model.
These tags are very similar to roxygen documentation tags and are called noxygen tags for the purpose of this package.
Noxygen tags are organized in three main groups:

* General model information
    * `@model Defines the title for the model documentation
    * `@author Defines author of the model (multiple author tags can be used to define several authors)
* Global definitions
    * `@global Defines the name of a global definitions (e.g. globals, patches-owns, breeds, ...)
    * `@details Further description of the global definition (multi line: each further line needs to start with details tag)
    * `@code TRUE/FALSE if following code should be included in the documentation
* Procedures
    * `@procedure Defines the name of a model procedure
    * `@param Defines a parameter that needs to be provided for the function call
    * `@return Description of the return value (in case of to-report functions)
    * `@details Further description of the model procedure (multi line: each further line needs to start with details tag)
    * `@code TRUE/FALSE if following code should be included in the documentation

For example, the start of a NetLogo model file with noxygen code could look like this:

```{r eval=FALSE}
;`@model Wolf Sheep Predation (NetLogo models library)
;`@author Uri Wilenski

;`@global Subfiles
;`@details We included a subilfe to the wolf sheep model to demonstrate that nls files are supported by nldoc
;`@code TRUE
__includes["Wolf Sheep Predation Extra.nls"]

;`@global Global variables
;`@details There is only one global variable that stores the max number of sheep allowed
;`@code TRUE
globals [ max-sheep ]  ; dont let sheep population grow too large

;`@global Breeds
;`@details There are two breeds: sheep and wolves
;`@code TRUE
breed [ sheep a-sheep ]  ; sheep is its own plural, so we use "a-sheep" as the singular.
breed [ wolves wolf ]

;`@global Agent properties
;`@details Sheep and wolves have a energy variable. Patches have a countdown variable to store the current state of the grass regrowth countdown.
;`@code TRUE
turtles-own [ energy ]       ; both wolves and sheep have energy
patches-own [ countdown ]

;`@procedure Setup
;`@details The setup procedure first resets the model.
;`@details Depending on the chosen model version, grass patches are initialized.
;`@details Finally, wolves and sheep are created.
;`@code FALSE
to setup
  ...
end

```

One, or more model files containing noxygen tags, can be used to call the `nldoc` function of the nlrx package.
The procedure creates a markdown file and renders the documentation in the specified output.
As an example, we added noxygen tags to the Wolf Sheep model from the NetLogo models library.
These modified model files are hosted on github. We use these files to render a documentation in html.
Example output: [nldoc documentation](nlrx_nldocumentation.html)

```{r eval=FALSE}
# Load nlrx:
library(nlrx)
outpath <- tempdir()

# List model files (.nls subfiles are also supported)
modelfiles <- c("https://raw.githubusercontent.com/nldoc/nldoc_pg/master/WSP.nlogo",
                "https://raw.githubusercontent.com/nldoc/nldoc_pg/master/WSP.nls")



## Create documentation:
# Themes: "journal", "cerulean", "flatly", "readable", "spacelab", "united", "cosmo"
# output_format: "pdf "html" "docx"
nldoc(modelfiles = modelfiles,
      infotab=TRUE,
      gui=TRUE,
      bs=TRUE,
      outpath = outpath,
      output_format = "html",
      number_sections = TRUE,
      theme = "cosmo",
      date = date(),
      toc = TRUE)
```

#### Model procedure network

Another useful function of the nlrx package is the creation of model procedure graphs (`nldoc_network()`).
This function can be used for any NetLogo model, even if the code does not contain noxygen tags.

```{r eval=FALSE}
# Load nlrx:
library(nlrx)

# List model files (.nls subfiles are also supported)
modelfiles <- c("https://raw.githubusercontent.com/nldoc/nldoc_pg/master/WSP.nlogo",
                "https://raw.githubusercontent.com/nldoc/nldoc_pg/master/WSP.nls")


## Determine the function network with nldoc:
nw <- nldoc_network(modelfiles)

## Determine communities within the network and plot using Igraph package:
library(igraph)
com <- walktrap.community(nw)
V(nw)$community <- com$membership
rain <- rainbow(14, alpha=.5)
V(nw)$color <- rain[V(nw)$community]

plot(nw,
     edge.arrow.size=1,
     vertex.label.color="black",
     vertex.label.dist=2.5,
     vertex.size=10,
     edge.curved=0,
     vertex.label.cex=1.5,
     layout=layout_with_fr(nw, niter = 2000))

## Interactive plot using igraph::tkplot
tkplot(nw, layout=layout_with_fr(nw, niter = 2000))
```

<center>
<img src="nldoc_network.png" align="center" width="75%" />
</center>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_simple}
\alias{nl_simple}
\title{Wolf Sheep model sample data: simdesign simple}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_simple
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a simple simulation design.
It also contains the model outputs within the simdesign object.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl.R
\name{util_call_nl}
\alias{util_call_nl}
\title{Setup and execute NetLogo via command line}
\usage{
util_call_nl(nl, xmlfile, outfile, batchfile)
}
\arguments{
\item{nl}{nl object}

\item{xmlfile}{file location of the experiment xml file}

\item{outfile}{file location for output results}

\item{batchfile}{file location of system specific batch file to call NetLogo
via command line}
}
\description{
Setup and execute NetLogo via command line
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nldoc_table_bs.R
\name{nldoc_table_bs}
\alias{nldoc_table_bs}
\title{Read NetLogo behavior space experiments from files}
\usage{
nldoc_table_bs(modelfiles)
}
\arguments{
\item{modelfiles}{vector of filepaths to model files}
}
\value{
list containing NetLogo behavior space experiments
}
\description{
Read NetLogo behavior space experiments from files
}
\details{
The procedure reads text from the provided model files and reports a list of behavior space experiments.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl.R
\name{util_create_sim_XML}
\alias{util_create_sim_XML}
\title{Create a temporary behavior space xml file to setup NetLogo via command line}
\usage{
util_create_sim_XML(nl, seed, siminputrow, xmlfile)
}
\arguments{
\item{nl}{nl object}

\item{seed}{random-seed for NetLogo simulation}

\item{siminputrow}{row id of the simulation input tibble of the simdesign
within the provided nl object}

\item{xmlfile}{filepath where the xml file is stored}
}
\description{
Create a temporary behavior space xml file to setup NetLogo via
command line
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_ABCmcmc_Marjoram}
\alias{simdesign_ABCmcmc_Marjoram}
\title{Add an Approximate Bayesian Computation (Monte-Carlo Markov-Chain) simdesign using the Majoram algorithm to a nl object}
\usage{
simdesign_ABCmcmc_Marjoram(
  nl,
  postpro_function = NULL,
  summary_stat_target,
  prior_test = NULL,
  n_rec,
  n_between_sampling = 10,
  n_cluster = 1,
  use_seed = FALSE,
  dist_weights = NULL,
  n_calibration = 10000,
  tolerance_quantile = 0.01,
  proposal_phi = 1,
  seed_count = 0,
  progress_bar = FALSE,
  nseeds
)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{postpro_function}{default is NULL. Allows to provide a function that is called to post-process the output Tibble of the NetLogo simulations. The function must accept the nl object with attached results as input argument. The function must return a one-dimensional vector of output metrics that corresponds in leght and order to the specified summary_stat_target.}

\item{summary_stat_target}{a vector of target values in the same order as the defined metrics of the experiment}

\item{prior_test}{a string expressing the constraints between model parameters. This expression will be evaluated as a logical expression, you can use all the logical operators including "<", ">", ... Each parameter should be designated with "X1", "X2", ... in the same order as in the prior definition. Set to NULL to disable.}

\item{n_rec}{Number of samples along the MCMC}

\item{n_between_sampling}{a positive integer equal to the desired spacing between sampled points along the MCMC.}

\item{n_cluster}{number of cores to parallelize simulations. Due to the design of the EasyABC parallelization it is currently not possible to use this feature with cores > 1.}

\item{use_seed}{if TRUE, seeds will be automatically created for each new model run}

\item{dist_weights}{a vector containing the weights to apply to the distance between the computed and the targeted statistics. These weights can be used to give more importance to a summary statistisc for example. The weights will be normalized before applying them. Set to NULL to disable.}

\item{n_calibration}{a positive integer. This is the number of simulations performed during the calibration step. Default value is 10000.}

\item{tolerance_quantile}{a positive number between 0 and 1 (strictly). This is the percentage of simulations retained during the calibration step to determine the tolerance threshold to be used during the MCMC. Default value is 0.01.}

\item{proposal_phi}{a positive number. This is a scaling factor defining the range of MCMC jumps. Default value is 1.}

\item{seed_count}{a positive integer, the initial seed value provided to the function model (if use_seed=TRUE). This value is incremented by 1 at each call of the function model.}

\item{progress_bar}{logical, FALSE by default. If TRUE, ABC_mcmc will output a bar of progression with the estimated remaining computing time. Option not available with multiple cores.}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add an Approximate Bayesian Computation (Monte-Carlo Markov-Chain) simdesign using the Majoram algorithm to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min, max and a shape of the distribution (qunif, qnorm, qlnorm, qexp)(e.g. list(min=1, max=4, qfun="qunif")).

The function uses the EasyABC package to set up the ABC_mcmc function.
For details on the ABC_mcmc function parameters see ?EasyABC::ABC_mcmc
Finally, the function reports a simdesign object.

Approximate Bayesian Computation simdesigns can only be executed using the \link[nlrx]{run_nl_dyn} function instead of \link[nlrx]{run_nl_all} or \link[nlrx]{run_nl_one}.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_lhs

# Attach the simdesign to the nl object
nl@simdesign <- simdesign_ABCmcmc_Marjoram(nl = nl,
                                            summary_stat_target = c(100, 80),
                                            n_rec = 100,
                                            n_between_sampling = 10,
                                            n_cluster = 1,
                                            use_seed = FALSE,
                                            n_calibration = 10000,
                                            tolerance_quantile = 0.01,
                                            proposal_phi = 1,
                                            progress_bar = FALSE,
                                            nseeds = 1)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl_dyn.R
\name{util_run_nl_dyn_GenAlg}
\alias{util_run_nl_dyn_GenAlg}
\title{Genetic Algorithm call simulations function}
\usage{
util_run_nl_dyn_GenAlg(nl, seed, cleanup.csv, cleanup.xml, cleanup.bat)
}
\arguments{
\item{nl}{nl object}

\item{seed}{current model seed}

\item{cleanup.csv}{TRUE/FALSE, if TRUE temporary created csv output files will be deleted after gathering results.}

\item{cleanup.xml}{TRUE/FALSE, if TRUE temporary created xml output files will be deleted after gathering results.}

\item{cleanup.bat}{TRUE/FALSE, if TRUE temporary created bat/sh output files will be deleted after gathering results.}
}
\description{
Genetic Algorithm call simulations function
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nldoc_parse_modelcode.R
\name{nldoc_parse_modelcode}
\alias{nldoc_parse_modelcode}
\title{Parse model code}
\usage{
nldoc_parse_modelcode(nlogocode)
}
\arguments{
\item{nlogocode}{vector of netlogo code strings}
}
\value{
tibble with structured netlogo code
}
\description{
Parse model code
}
\details{
The procedure searches for noxygen commands within the NetLogo code.
This information is used to structure the netlogo code strings in a tibble.
Additionally, tibbles with gui and behavior space information are created.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_nl.R
\name{import_nl}
\alias{import_nl}
\title{Import NetLogo Experiment}
\usage{
import_nl(tarfile, targetdir, new_session = FALSE)
}
\arguments{
\item{tarfile}{Path to tarfile that contains files to run NetLogo experiment}

\item{targetdir}{Path to folder where the experiments gets extracted}

\item{new_session}{If TRUE, opens a new RStudio Session with an Rproj}
}
\value{
The status value returned by the external command, invisibly.
}
\description{
Import NetLogo Experiment from export_nl
}
\details{
Imports NetLogo experiments that were saved with \code{export_nl}.
If the folder comes with an .Rproj file (which is recommended because
relative paths enhance the reproducability of your analysis),
\code{import_nl} opens this project and loads the nl object in your R environment.
}
\examples{
\dontrun{

infile <- "/home/user/test.zip"
targetdirectory <- "/home/user/test"
import_nl(infile, targetdirectory)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_nl.R
\name{analyze_sobol2007}
\alias{analyze_sobol2007}
\title{Analyze NetLogo simulation output of simdesign sobol2007}
\usage{
analyze_sobol2007(nl, metrics, funs)
}
\arguments{
\item{nl}{nl object}

\item{metrics}{vector of strings defining metric columns for evaluation. Defaults to metrics of the experiment within the nl object}

\item{funs}{list with the summary metrics for the sensitivity results}
}
\description{
Analyze NetLogo simulation output of simdesign sobol2007
}
\details{
The function calculates sobol sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl_dyn.R
\name{util_run_nl_dyn_GenSA_fn}
\alias{util_run_nl_dyn_GenSA_fn}
\title{Simulated Annealing run simulation function}
\usage{
util_run_nl_dyn_GenSA_fn(
  param,
  nl,
  evalcrit,
  seed,
  cleanup.csv,
  cleanup.xml,
  cleanup.bat
)
}
\arguments{
\item{param}{vector of model parameters, generated by GenSA function}

\item{nl}{nl object}

\item{evalcrit}{evaluation criterion for simulated annealing}

\item{seed}{current model seed}

\item{cleanup.csv}{TRUE/FALSE, if TRUE temporary created csv output files will be deleted after gathering results.}

\item{cleanup.xml}{TRUE/FALSE, if TRUE temporary created xml output files will be deleted after gathering results.}

\item{cleanup.bat}{TRUE/FALSE, if TRUE temporary created bat/sh output files will be deleted after gathering results.}
}
\description{
Simulated Annealing run simulation function
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_distinct}
\alias{simdesign_distinct}
\title{Add a distinct simdesign to a nl object}
\usage{
simdesign_distinct(nl, nseeds)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add a distinct simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.
The distinct simdesign allows to create a parameter matrix with distinct parameterisations.

Variables in the experiment variable list need to provide a vector of distinct values (e.g. list(values=c(1,2,3,4)).
All vectors of values must have the same length across variables.

The distinct simdesign then creates one simulation run for all first elements of these values vectors,
one run for all second items, and so on.
With this function, multiple distinct simulations can be run at once.
Finally, the function reports a simdesign object.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_distinct
nl@simdesign <- simdesign_distinct(nl = nl, nseeds = 3)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nldoc_write_nldoc.R
\name{nldoc_write_nldoc}
\alias{nldoc_write_nldoc}
\title{Write NetLogo documentation}
\usage{
nldoc_write_nldoc(
  noxygen,
  noxygen_it,
  noxygen_gui,
  noxygen_bs,
  outpath,
  output_format,
  number_sections,
  theme,
  date,
  toc
)
}
\arguments{
\item{noxygen}{list with parsed and processed noxygen tags from NetLogo model code}

\item{noxygen_it}{list with parsed and processed infotab strings}

\item{noxygen_gui}{list with parsed and processed noxygen tags from NetLogo GUI elements}

\item{noxygen_bs}{list with parsed and processed noxygen tags from NetLogo behavior space experiments}

\item{outpath}{Path to folder where rendered documentation should be created}

\item{output_format}{either "html", "pdf" or "docx"}

\item{theme}{markdown theme, supported themes are "journal", "cerulean", "flatly", "readable", "spacelab", "united", "cosmo"}

\item{date}{date that is printed in the documentation header}

\item{toc}{TRUE/FALSE, if TRUE the documentation contains a table of contents - only for html and pdf output format}

\item{number_ections}{TRUE/FALSE, if TRUE sections in the documentation will be numbered}
}
\value{
list containing NetLogo GUI elements
}
\description{
Write NetLogo documentation
}
\details{
The procedure uses parsed and processed noxygen tags to create and render a markdown documentation in the specified format.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_netlogo.R
\name{download_netlogo}
\alias{download_netlogo}
\title{Download NetLogo}
\usage{
download_netlogo(to, version, extract = FALSE)
}
\arguments{
\item{to}{Path to folder where the downloaded file is saved.}

\item{version}{Character string naming which NetLogo Version to download (see Details)}

\item{extract}{TRUE/FALSE, if TRUE downloaded archive is extracted to subfolder of \code{to} (only unix)}
}
\description{
Auxiliary function to download NetLogo
}
\details{
Supported Versions for Download and Usage (parameter \code{version}):
\itemize{
\item "6.2.0" = NetLogo Version 6.2.0
\item "6.1.1" = NetLogo Version 6.1.1
\item "6.1.0" = NetLogo Version 6.1.0
\item "6.0.4" = NetLogo Version 6.0.4
\item "6.0.3" = NetLogo Version 6.0.3
\item "6.0.2" = NetLogo Version 6.0.2
\item "6.0.1" = NetLogo Version 6.0.1
\item "6.0" = NetLogo Version 6.0
\item "5.3.1" = NetLogo Version 5.3.1
}
}
\examples{
\dontrun{
dlpath <- tempdir()  # adjust path to your needs
try(download_netlogo(dlpath, "6.0.3"))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_eval.R
\name{util_eval_variables_distinct}
\alias{util_eval_variables_distinct}
\title{Evaluate variables list of an experiment object for distinct simdesign}
\usage{
util_eval_variables_distinct(nl)
}
\arguments{
\item{nl}{nl object}
}
\description{
Evaluate variables list of an experiment object for distinct simdesign
}
\details{
util_eval_variables_distinct checks if the variables list of an experiment within a nl object has enough information to create a \link[nlrx]{simdesign_distinct}.
It reports an error message if at least one variable does not have a vector of distinct values or if there is a mismatch in length of these values vectors.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_setget.R
\name{setexp<-}
\alias{setexp<-}
\alias{setexp}
\title{Setter function to set a variable of an experiment object}
\usage{
setexp(nl, var) <- value
}
\arguments{
\item{nl}{nl object}

\item{var}{valid experiment variable string}

\item{value}{valid value for the specified variable}
}
\description{
Setter function to set a variable of an experiment object
}
\examples{

# Example for Wolf Sheep Predation model from NetLogo models library:
nl <- nl(nlversion = "6.0.3",
nlpath = "/home/user/NetLogo 6.0.3/",
modelpath = "/home/user/NetLogo 6.0.3/app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo",
jvmmem = 1024)

# Set experiment name
setexp(nl, "expname") <- "experimentName"

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_GenSA}
\alias{simdesign_GenSA}
\title{Add a Simulated Annealing simdesign to a nl object}
\usage{
simdesign_GenSA(nl, par = NULL, evalcrit = 1, control = list(), nseeds = 1)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{par}{optional vector of start values for each parameter defined in variables of experiment}

\item{evalcrit}{position of evaluation criterion within defined NetLogo metrics of nl experiment or a function that reports a single numeric value}

\item{control}{list with further arguments passed to the GenSA function (see ?GenSA for details)}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add a Simulated Annealing simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min and max (e.g. list(min=1, max=4)).

The GenSA simdesign generates a simulated Annealing experiment within the defined min and max parameter boundaries
that are defined in the variables field of the experiment object within the nl object.

The evalcrit reporter defines the evaluation criterion for the simulated annealing procedure.
There are two options to evaluate the fitness value of each iteration of the algorithm:
\enumerate{
\item Use a reporter that is defined within the experiment metrics vector.
You can just enter the position of that metric within the experiment metrics vector (e.g. 1 would use the first defined metric of the experiment to evaluate each iteration).
The algorithm automatically calculates the mean value of this reporter if evalticks is defined to measure multiple ticks during each simulation.
You can define a function that post-processes NetLogo output to calculate an evaluation value. This function must accept the nl object as input and return one single numeric value.
The nl object that is then provided to the evaluation function will have results of the current iteration attached. The results can be accessed via the simoutput slot of the simdesign.
You can pass this function to evalcrit. It is then applied to the output of each iteration.
}

The function uses the GenSA package to set up a Simulated Annealing function.
For details on the GenSA function parameters see ?GenSA
Finally, the function reports a simdesign object.

Simulated Annealing simdesigns can only be executed using the \link[nlrx]{run_nl_dyn} function instead of \link[nlrx]{run_nl_all} or \link[nlrx]{run_nl_one}.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_lhs

# Example 1: Using a metric from the experiment metrics vector for evaluation:
nl@simdesign <- simdesign_GenSA(nl=nl,
                                 par=NULL,
                                 evalcrit=1,
                                 control=list(max.time = 600),
                                 nseeds=1)


# Example 2: Using a self-defined evaluation function
# For demonstration we define a simple function that calculates
# the maximum value of count sheep output.
critfun <- function(nl) {
results <- nl@simdesign@simoutput
crit <- as.integer(max(results$`count sheep`))
return(crit)
}

nl@simdesign <- simdesign_GenSA(nl=nl,
                                 par=NULL,
                                 evalcrit=critfun,
                                 control=list(max.time = 600),
                                 nseeds=1)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_eval.R
\name{util_eval_variables_op}
\alias{util_eval_variables_op}
\title{Evaluate variables list of an experiment object for optimization simdesigns}
\usage{
util_eval_variables_op(nl)
}
\arguments{
\item{nl}{nl object}
}
\description{
Evaluate variables list of an experiment object for optimization simdesigns
}
\details{
util_eval_variables_op checks if the variables list of an experiment within a nl object has enough information to create an optimization simdesign.
It reports an error message if at least one variable does not have a defined range (min, max).
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nldoc_network.R
\name{nldoc_network}
\alias{nldoc_network}
\title{Create NetLogo procedure network}
\usage{
nldoc_network(modelfiles)
}
\arguments{
\item{modelfiles}{vector of filepaths to model files}
}
\value{
network of model procedures (igraph)
}
\description{
Create NetLogo procedure network
}
\details{
Reads model code from the provided model files.
The procedure identifies NetLogo procedures and searches for procedure calls within the code.
From this data, an igraph network is created and returned.
This network can be used to plot the model procedure network and identify model components.
}
\examples{
\dontrun{

# List model files (.nls subfiles are also supported)
modelfiles <- c("https://raw.githubusercontent.com/nldoc/nldoc_pg/master/WSP.nlogo",
                "https://raw.githubusercontent.com/nldoc/nldoc_pg/master/WSP.nls")

# Determine the function network:
nw <- nldoc_network(modelfiles)

# Determine communities within the network and plot using Igraph package:
library(igraph)
com <- walktrap.community(nw)
V(nw)$community <- com$membership
rain <- rainbow(14, alpha=.5)
V(nw)$color <- rain[V(nw)$community]

plot(nw,
     edge.arrow.size=.2,
     vertex.label.color="black",
     vertex.label.dist=1,
     vertex.size=5,
     edge.curved=0,
     vertex.label.cex=.5,
     layout=layout_with_fr(nw, niter = 2000))

# Interactive plot using igraph::tkplot
tkplot(nw, layout=layout_with_fr(nw, niter = 2000))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unnest_simoutput.R
\name{unnest_simoutput}
\alias{unnest_simoutput}
\title{Get spatial data from metrics.turtles and metrics.patches output}
\usage{
unnest_simoutput(nl)
}
\arguments{
\item{nl}{nl object}
}
\value{
tibble with spatial data objects
}
\description{
Turn results from NetLogo in spatial data objects
}
\details{
Unnests output from run_nl into long format.
}
\examples{

# To unnest data a nl object containing spatial output data is needed.
# For this example, we load a nl object from test data.

nl <- nl_spatial
unnest_simoutput(nl)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_nl.R
\name{analyze_nl}
\alias{analyze_nl}
\title{Analyze NetLogo simulation output}
\usage{
analyze_nl(nl, metrics = getexp(nl, "metrics"), funs = list(mean = mean))
}
\arguments{
\item{nl}{nl object}

\item{metrics}{vector of strings defining metric columns for evaluation. Defaults to metrics of the experiment within the nl object}

\item{funs}{list with the summary metrics for the sensitivity results}
}
\value{
analysis summary tibble
}
\description{
Analyze NetLogo simulation output
}
\details{
The analyze_nl function runs basic analyses on NetLogo simulation output.
In order to execute this function, simulation output needs to be attached to the simdesign first with \code{setsim(nl, "output") <- results}.

analyze_nl calls different post-processing analysis functions, depending on the specified method in the simdesign object of the nl object.

\strong{The following simdesign are currently supported:}

\link[nlrx]{simdesign_ff}

Calls \link[nlrx]{analyze_ff}.
The function calculates aggregated output metrics by dropping random seeds and aggregating values with the provided functions.

\link[nlrx]{simdesign_lhs}

Calls \link[nlrx]{analyze_lhs}.
The function calculates aggregated output metrics by dropping random seeds and aggregating values with the provided functions.

\link[nlrx]{simdesign_sobol}

Calls \link[nlrx]{analyze_sobol}.
The function calculates sobol sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.

\link[nlrx]{simdesign_sobol2007}

Calls \link[nlrx]{analyze_sobol2007}.
The function calculates sobol sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.

\link[nlrx]{simdesign_soboljansen}

Calls \link[nlrx]{analyze_soboljansen}.
The function calculates sobol sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.

\link[nlrx]{simdesign_morris}

Calls \link[nlrx]{analyze_morris}.
The function calculates morris sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.

\link[nlrx]{simdesign_eFast}

Calls \link[nlrx]{analyze_eFast}.
The function calculates eFast sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.

\strong{For the following simdesign no postprocessing analysis function has been implemented yet:}

\link[nlrx]{simdesign_simple},
\link[nlrx]{simdesign_distinct},
\link[nlrx]{simdesign_GenSA},
\link[nlrx]{simdesign_GenAlg}
}
\examples{

# Load nl object including output data from testdata
nl <- nl_sobol

# Define aggregation measurements:
myfuns <- list(mean=mean, sd=sd, min=min, max=max)

# Calculate sensitivity indices:
analyze_nl(nl, funs = myfuns)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_eval.R
\name{util_eval_constants}
\alias{util_eval_constants}
\title{Evaluate if constants list of an experiment object is empty}
\usage{
util_eval_constants(nl)
}
\arguments{
\item{nl}{nl object
util_eval_constants checks if the constants list of an experiment within a nl object is empty.
It reports an error message if no constants were defined. This evaluation is only done for \link[nlrx]{simdesign_simple}.}
}
\description{
Evaluate if constants list of an experiment object is empty
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_spatial}
\alias{nl_spatial}
\title{Wolf Sheep model sample data: spatial}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_spatial
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a simple simulation design.
For these simulations, spatial agent output has been collected.
metrics.turtles was set up to measure turtle coordinates.
metrics.patches was set up to measure colors of cells.
It also contains the model outputs within the simdesign object.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_ABCmcmc_Wegmann}
\alias{simdesign_ABCmcmc_Wegmann}
\title{Add an Approximate Bayesian Computation (Monte-Carlo Markov-Chain) simdesign using the Wegmann algorithm to a nl object}
\usage{
simdesign_ABCmcmc_Wegmann(
  nl,
  postpro_function = NULL,
  summary_stat_target,
  prior_test = NULL,
  n_rec,
  n_between_sampling = 10,
  n_cluster = 1,
  use_seed = FALSE,
  dist_weights = NULL,
  n_calibration = 10000,
  tolerance_quantile = 0.01,
  proposal_phi = 1,
  numcomp = 0,
  seed_count = 0,
  progress_bar = FALSE,
  nseeds
)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{postpro_function}{default is NULL. Allows to provide a function that is called to post-process the output Tibble of the NetLogo simulations. The function must accept the nl object with attached results as input argument. The function must return a one-dimensional vector of output metrics that corresponds in leght and order to the specified summary_stat_target.}

\item{summary_stat_target}{a vector of target values in the same order as the defined metrics of the experiment}

\item{prior_test}{a string expressing the constraints between model parameters. This expression will be evaluated as a logical expression, you can use all the logical operators including "<", ">", ... Each parameter should be designated with "X1", "X2", ... in the same order as in the prior definition. Set to NULL to disable.}

\item{n_rec}{Number of samples along the MCMC}

\item{n_between_sampling}{a positive integer equal to the desired spacing between sampled points along the MCMC.}

\item{n_cluster}{number of cores to parallelize simulations. Due to the design of the EasyABC parallelization it is currently not possible to use this feature with cores > 1.}

\item{use_seed}{if TRUE, seeds will be automatically created for each new model run}

\item{dist_weights}{a vector containing the weights to apply to the distance between the computed and the targeted statistics. These weights can be used to give more importance to a summary statistisc for example. The weights will be normalized before applying them. Set to NULL to disable.}

\item{n_calibration}{a positive integer. This is the number of simulations performed during the calibration step. Default value is 10000.}

\item{tolerance_quantile}{a positive number between 0 and 1 (strictly). This is the percentage of simulations retained during the calibration step to determine the tolerance threshold to be used during the MCMC. Default value is 0.01.}

\item{proposal_phi}{a positive number. This is a scaling factor defining the range of MCMC jumps. Default value is 1.}

\item{numcomp}{a positive integer. This is the number of components to be used for PLS transformations. Default value is 0 which encodes that this number is equal to the number of summary statistics.}

\item{seed_count}{a positive integer, the initial seed value provided to the function model (if use_seed=TRUE). This value is incremented by 1 at each call of the function model.}

\item{progress_bar}{logical, FALSE by default. If TRUE, ABC_mcmc will output a bar of progression with the estimated remaining computing time. Option not available with multiple cores.}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add an Approximate Bayesian Computation (Monte-Carlo Markov-Chain) simdesign using the Wegmann algorithm to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min, max and a shape of the distribution (qunif, qnorm, qlnorm, qexp)(e.g. list(min=1, max=4, qfun="qunif")).

The function uses the EasyABC package to set up the ABC_mcmc function.
For details on the ABC_mcmc function parameters see ?EasyABC::ABC_mcmc
Finally, the function reports a simdesign object.

Approximate Bayesian Computation simdesigns can only be executed using the \link[nlrx]{run_nl_dyn} function instead of \link[nlrx]{run_nl_all} or \link[nlrx]{run_nl_one}.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_lhs

# Attach the simdesign to the nl object
nl@simdesign <- simdesign_ABCmcmc_Wegmann(nl = nl,
                                            summary_stat_target = c(100, 80),
                                            n_rec = 100,
                                            n_between_sampling = 10,
                                            n_cluster = 1,
                                            use_seed = FALSE,
                                            n_calibration = 10000,
                                            tolerance_quantile = 0.01,
                                            proposal_phi = 1,
                                            progress_bar = FALSE,
                                            nseeds = 1)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_lhs}
\alias{nl_lhs}
\title{Wolf Sheep model sample data: simdesign lhs}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_lhs
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a latin hypercube simulation design.
It also contains the model outputs within the simdesign object.
Further analysis output can be derived by submitting the dataset to analyze_nl().
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nldoc_read_nlogo.R
\name{nldoc_read_nlogo}
\alias{nldoc_read_nlogo}
\title{Read NetLogo model code from files}
\usage{
nldoc_read_nlogo(modelfiles)
}
\arguments{
\item{modelfiles}{vector of filepaths to model files}
}
\value{
vector of strings containing NetLogo model code
}
\description{
Read NetLogo model code from files
}
\details{
The procedure reads text from the provided model files and reports the code as a vector of strings.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_gensa}
\alias{nl_gensa}
\title{Wolf Sheep model sample data: gensa}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_gensa
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a simulated annealing simulation design.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_sobol2007}
\alias{simdesign_sobol2007}
\title{Add a sobol2007 simdesign to a nl object}
\usage{
simdesign_sobol2007(nl, samples, sobolnboot, sobolconf, nseeds, precision)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{samples}{number of samples for the sobol sensitivity analysis}

\item{sobolnboot}{number of bootstrap replicates of the sobol sensitivity analysis}

\item{sobolconf}{the confidence level for bootstrap confidence intervals}

\item{nseeds}{number of seeds for this simulation design}

\item{precision}{number of digits for the decimal fraction of parameter values}
}
\value{
simdesign S4 class object
}
\description{
Add a sobol2007 simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min, max and qfun (e.g. list(min=1, max=4, qfun="qunif")).

The sobol2007 simdesign uses the sensitivity package to set up a sobol2007 sensitivity analysis, including a simobject of class sobol and a input tibble for simulations.
For details on method specific sensitivity analysis function parameters see ?sobol2007
Finally, the function reports a simdesign object.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_sobol2007
nl@simdesign <- simdesign_sobol2007(nl=nl,
samples=1000,
sobolnboot=100,
sobolconf=0.95,
nseeds=3,
precision=3)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_eFast}
\alias{nl_eFast}
\title{Wolf Sheep model sample data: simdesign eFast}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_eFast
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run an eFast simulation design.
It also contains the model outputs within the simdesign object.
Further analysis output can be derived by submitting the dataset to analyze_nl().
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_ff}
\alias{nl_ff}
\title{Wolf Sheep model sample data: simdesign ff}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_ff
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a full-factorial simulation design.
It also contains the model outputs within the simdesign object.
Further analysis output can be derived by submitting the dataset to analyze_nl().
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_nl.R
\name{util_print.summary}
\alias{util_print.summary}
\title{Print nl object summary}
\usage{
util_print.summary(x, ...)
}
\arguments{
\item{x}{nl object}

\item{...}{further arguments passed to or from other methods}
}
\description{
Print nl object summary
}
\details{
Print summary of nl object and embedded experiment and simdesign objects in a nice formatted overview to console
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_nl.R
\name{analyze_eFast}
\alias{analyze_eFast}
\title{Analyze NetLogo simulation output of simdesign eFast}
\usage{
analyze_eFast(nl, metrics, funs)
}
\arguments{
\item{nl}{nl object}

\item{metrics}{vector of strings defining metric columns for evaluation. Defaults to metrics of the experiment within the nl object}

\item{funs}{list with the summary metrics for the sensitivity results}
}
\description{
Analyze NetLogo simulation output of simdesign eFast
}
\details{
The function calculates eFast sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_morris}
\alias{nl_morris}
\title{Wolf Sheep model sample data: simdesign morris}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_morris
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a morris simulation design.
It also contains the model outputs within the simdesign object.
Further analysis output can be derived by submitting the dataset to analyze_nl().
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nldoc_table_gui.R
\name{nldoc_table_gui}
\alias{nldoc_table_gui}
\title{Read NetLogo GUI elements from files}
\usage{
nldoc_table_gui(modelfiles)
}
\arguments{
\item{modelfiles}{vector of filepaths to model files}
}
\value{
list containing NetLogo GUI elements
}
\description{
Read NetLogo GUI elements from files
}
\details{
The procedure reads text from the provided model files and reports a list of GUI elements.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_nl.R
\name{util_print.simdesign}
\alias{util_print.simdesign}
\title{Print simdesign object content}
\usage{
util_print.simdesign(x, ...)
}
\arguments{
\item{x}{simdesign object}

\item{...}{further arguments passed to or from other methods}
}
\description{
Print simdesign object content
}
\details{
Print content of simdesign object in a nice formatted overview to console
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nldoc.R
\name{nldoc}
\alias{nldoc}
\title{Create NetLogo documentation}
\usage{
nldoc(
  modelfiles,
  outpath,
  infotab = TRUE,
  gui = TRUE,
  bs = TRUE,
  output_format = "html",
  number_sections = TRUE,
  theme = "journal",
  date = as.Date(Sys.time()),
  toc = TRUE
)
}
\arguments{
\item{modelfiles}{vector of filepaths to model files}

\item{outpath}{Path to folder where rendered documentation should be created}

\item{infotab}{TRUE/FALSE, if TRUE infotab section will be included in the documentation}

\item{gui}{TRUE/FALSE, if TRUE a table with GUI elements from the model will be included in the documentation}

\item{bs}{TRUE/FALSE, if TRUE a table with behavior space experiments will be included in the documentation}

\item{output_format}{either "html", "pdf" or "docx"}

\item{number_sections}{TRUE/FALSE, if TRUE sections in the documentation will be numbered}

\item{theme}{markdown theme, supported themes are "journal", "cerulean", "flatly", "readable", "spacelab", "united", "cosmo"}

\item{date}{date that is printed in the documentation header}

\item{toc}{TRUE/FALSE, if TRUE the documentation contains a table of contents - only for html and pdf output format}
}
\description{
Create NetLogo documentation
}
\details{
nldoc reads model code from the provided model files.
The code is then split into several groups (code, gui, behavior space).
The procedures then finds noxygen commands within the NetLogo model code.
For a complete list of noxygen commands type ?nldoc
These commands are translated into a markdown documentation file.
If needed, tables of gui elements and behavior space experiments are added to the markdown file.
Finally, the document is rendered in the specified format.
}
\examples{
\dontrun{

# List model files (.nls subfiles are also supported)
modelfiles <- c("https://raw.githubusercontent.com/nldoc/nldoc_pg/master/WSP.nlogo",
                "https://raw.githubusercontent.com/nldoc/nldoc_pg/master/WSP.nls")

# Define output directory:
outdir <- tempdir()  # adjust path to your needs

# Create documentation:
nldoc(modelfiles = modelfiles,
      infotab=TRUE,
      gui=TRUE,
      bs=TRUE,
      outpath = outdir,
      output_format = "html",
      number_sections = TRUE,
      theme = "cosmo",
      date = date(),
      toc = TRUE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl.R
\name{util_read_write_batch}
\alias{util_read_write_batch}
\title{Write a modified batchfile that executes NetLogo}
\usage{
util_read_write_batch(nl)
}
\arguments{
\item{nl}{nl object}
}
\description{
Write a modified batchfile that executes NetLogo
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_ABCmcmc_Marjoram_original}
\alias{simdesign_ABCmcmc_Marjoram_original}
\title{Add an Approximate Bayesian Computation (Monte-Carlo Markov-Chain) simdesign using the Majoram Original algorithm to a nl object}
\usage{
simdesign_ABCmcmc_Marjoram_original(
  nl,
  postpro_function = NULL,
  summary_stat_target,
  prior_test = NULL,
  n_rec,
  n_between_sampling = 10,
  n_cluster = 1,
  use_seed = FALSE,
  dist_weights = NULL,
  dist_max = 0,
  tab_normalization = summary_stat_target,
  proposal_range = vector(mode = "numeric", length = length(getexp(nl, "variables"))),
  seed_count = 0,
  progress_bar = FALSE,
  nseeds
)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{postpro_function}{default is NULL. Allows to provide a function that is called to post-process the output Tibble of the NetLogo simulations. The function must accept the nl object with attached results as input argument. The function must return a one-dimensional vector of output metrics that corresponds in leght and order to the specified summary_stat_target.}

\item{summary_stat_target}{a vector of target values in the same order as the defined metrics of the experiment}

\item{prior_test}{a string expressing the constraints between model parameters. This expression will be evaluated as a logical expression, you can use all the logical operators including "<", ">", ... Each parameter should be designated with "X1", "X2", ... in the same order as in the prior definition. Set to NULL to disable.}

\item{n_rec}{Number of samples along the MCMC}

\item{n_between_sampling}{a positive integer equal to the desired spacing between sampled points along the MCMC.}

\item{n_cluster}{number of cores to parallelize simulations. Due to the design of the EasyABC parallelization it is currently not possible to use this feature with cores > 1.}

\item{use_seed}{if TRUE, seeds will be automatically created for each new model run}

\item{dist_weights}{a vector containing the weights to apply to the distance between the computed and the targeted statistics. These weights can be used to give more importance to a summary statistisc for example. The weights will be normalized before applying them. Set to NULL to disable.}

\item{dist_max}{a positive number. This is the tolerance threshold used during the MCMC. If not provided by the user, it is automatically computed as half the distance between the first simulation and the target summary statistics and a warning is printed.}

\item{tab_normalization}{a vector of the same length as summary_stat_target. Each element contains a positive number by which each summary statistics must be divided before the computation of the Euclidean distance between simulations and data. If not provided by the user, the simulated summary statistics are divided by the target summary statistics and a warning is printed.}

\item{proposal_range}{a vector of the same length as the number of model parameters, used when method is "Marjoram_original". Each element contains a positive number defining the range of MCMC jumps for each model parameter. If not provided by the user, a default value is used for each parameter and a warning is printed. The default value is 1/50 of the prior range for uniform distributions, 1/20 of the standard deviation of the prior distribution for normal distributions, 1/20 * exp ( sigma * sigma}

\item{seed_count}{a positive integer, the initial seed value provided to the function model (if use_seed=TRUE). This value is incremented by 1 at each call of the function model.}

\item{progress_bar}{logical, FALSE by default. If TRUE, ABC_mcmc will output a bar of progression with the estimated remaining computing time. Option not available with multiple cores.}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add an Approximate Bayesian Computation (Monte-Carlo Markov-Chain) simdesign using the Majoram Original algorithm to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min, max and a shape of the distribution (qunif, qnorm, qlnorm, qexp)(e.g. list(min=1, max=4, qfun="qunif")).

The function uses the EasyABC package to set up the ABC_mcmc function.
For details on the ABC_mcmc function parameters see ?EasyABC::ABC_mcmc
Finally, the function reports a simdesign object.

Approximate Bayesian Computation simdesigns can only be executed using the \link[nlrx]{run_nl_dyn} function instead of \link[nlrx]{run_nl_all} or \link[nlrx]{run_nl_one}.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_lhs

# Attach the simdesign to the nl object
nl@simdesign <- simdesign_ABCmcmc_Marjoram_original(nl = nl,
                                                     summary_stat_target = c(100, 80),
                                                     n_rec = 100,
                                                     n_between_sampling = 10,
                                                     nseeds = 1)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl_dyn.R
\name{util_run_nl_dyn_GenSA}
\alias{util_run_nl_dyn_GenSA}
\title{Simulated Annealing call simulations function}
\usage{
util_run_nl_dyn_GenSA(nl, seed, cleanup.csv, cleanup.xml, cleanup.bat)
}
\arguments{
\item{nl}{nl object}

\item{seed}{current model seed}

\item{cleanup.csv}{TRUE/FALSE, if TRUE temporary created csv output files will be deleted after gathering results.}

\item{cleanup.xml}{TRUE/FALSE, if TRUE temporary created xml output files will be deleted after gathering results.}

\item{cleanup.bat}{TRUE/FALSE, if TRUE temporary created bat/sh output files will be deleted after gathering results.}
}
\description{
Simulated Annealing call simulations function
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_sobol}
\alias{simdesign_sobol}
\title{Add a sobol simdesign to a nl object}
\usage{
simdesign_sobol(
  nl,
  samples,
  sobolorder,
  sobolnboot,
  sobolconf,
  nseeds,
  precision
)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{samples}{number of samples for the sobol sensitivity analysis}

\item{sobolorder}{order of interactions of the sobol sensitivity analysis}

\item{sobolnboot}{number of bootstrap replicates of the sobol sensitivity analysis}

\item{sobolconf}{the confidence level for bootstrap confidence intervals}

\item{nseeds}{number of seeds for this simulation design}

\item{precision}{number of digits for the decimal fraction of parameter values}
}
\value{
simdesign S4 class object
}
\description{
Add a sobol simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min, max and qfun (e.g. list(min=1, max=4, qfun="qunif")).

The sobol simdesign uses the sensitivity package to set up a sobol sensitivity analysis, including a simobject of class sobol and a input tibble for simulations.
For details on method specific sensitivity analysis function parameters see ?sobol
Finally, the function reports a simdesign object.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_sobol
nl@simdesign <- simdesign_sobol(nl=nl,
samples=1000,
sobolorder=2,
sobolnboot=100,
sobolconf=0.95,
nseeds=3,
precision=3)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_setget.R
\name{getsim}
\alias{getsim}
\title{Getter function to get a variable of a simdesign object}
\usage{
getsim(nl, var)
}
\arguments{
\item{nl}{nl object}

\item{var}{valid simdesign variable string}
}
\description{
Getter function to get a variable of a simdesign object
}
\examples{
# Example for Wolf Sheep Predation model from NetLogo models library:
nl <- nl(nlversion = "6.0.3",
nlpath = "/home/user/NetLogo 6.0.3/",
modelpath = "/home/user/NetLogo 6.0.3/app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo",
jvmmem = 1024)

# Set simulation seeds
setsim(nl, "simseeds") <- c(123, 456, 789)

# Set simulation seeds
getsim(nl, "simseeds")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constr.R
\name{nl}
\alias{nl}
\title{Construct a new nl object}
\usage{
nl(
  nlversion = "6.0.2",
  nlpath = character(),
  modelpath = character(),
  jvmmem = 1024,
  experiment = methods::new("experiment"),
  simdesign = methods::new("simdesign"),
  ...
)
}
\arguments{
\item{nlversion}{A character string defining the NetLogo version that is used}

\item{nlpath}{Path to the NetLogo main directory matching the defined version}

\item{modelpath}{Path to the NetLogo model file (*.nlogo) that is used for simulations}

\item{jvmmem}{Java virtual machine memory capacity in megabytes}

\item{experiment}{Holds a experiment S4 class object}

\item{simdesign}{Holds a simdesign S4 class object}

\item{...}{...}
}
\value{
nl S4 class object
}
\description{
Construct a new nl object
}
\details{
nl objects are the main class objects used in the nlrx package.
These objects store all information that is needed to run NetLogo simulations.
nl objects are initialized with basic information on Netlogo and the model.

After setting up the nl object, an experiment needs to be attached.
The experiment class stores all information related to the NetLogo simulation experiment, such as runtime,
variables, constants, measurements, and more.

After attaching an experiment, different simdesign helper functions can be used to attach a simdesign to the nl object.
The simdesign helper functions use the variable definitions from the experiment within the nl object to generate a parameter tibble for simulations.
}
\examples{
# Example for Wolf Sheep Predation model from NetLogo models library:
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_eval.R
\name{util_eval_variables}
\alias{util_eval_variables}
\title{Evaluate if variables list of an experiment object is empty}
\usage{
util_eval_variables(nl)
}
\arguments{
\item{nl}{nl object}
}
\description{
Evaluate if variables list of an experiment object is empty
}
\details{
util_eval_variables checks if the variables list of an experiment within a nl object is empty.
It reports an error message if no variables were defined.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_morris}
\alias{simdesign_morris}
\title{Add a morris elementary effects simdesign to a nl object}
\usage{
simdesign_morris(nl, morristype, morrislevels, morrisr, morrisgridjump, nseeds)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{morristype}{morris design type}

\item{morrislevels}{number of parameter levels}

\item{morrisr}{morris r value}

\item{morrisgridjump}{morris grid jump value}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add a morris elementary effects simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min, max and qfun (e.g. list(min=1, max=4, qfun="qunif")).

The morris simdesign uses the sensitivity package to set up a morris elementary effects sensitivity analysis, including a simobject of class morris and a input tibble for simulations.
For details on method specific sensitivity analysis function parameters see ?morris
Finally, the function reports a simdesign object.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_morris
nl@simdesign <- simdesign_morris(nl=nl,
                                  morristype="oat",
                                  morrislevels=4,
                                  morrisr=20,
                                  morrisgridjump=2,
                                  nseeds=3)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_setget.R
\name{getexp}
\alias{getexp}
\title{Getter function to get a variable of an experiment object}
\usage{
getexp(nl, var)
}
\arguments{
\item{nl}{nl object}

\item{var}{valid experiment variable string}
}
\description{
Getter function to get a variable of an experiment object
}
\examples{
# Example for Wolf Sheep Predation model from NetLogo models library:
nl <- nl(nlversion = "6.0.3",
nlpath = "/home/user/NetLogo 6.0.3/",
modelpath = "/home/user/NetLogo 6.0.3/app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo",
jvmmem = 1024)

# Set experiment name
setexp(nl, "expname") <- "experimentName"

# Get experiment name
getexp(nl, "experiment")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_stuff.R
\name{util_create_lhs}
\alias{util_create_lhs}
\title{Identify and report the current OS}
\usage{
util_create_lhs(input, samples, precision)
}
\arguments{
\item{input}{list with variables and value ranges}

\item{samples}{number of lhs samples}

\item{precision}{number of digits for the decimal fraction of parameter
values}
}
\description{
Identify and report the current OS
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_setget.R
\name{setsim<-}
\alias{setsim<-}
\alias{setsim}
\title{Setter function to set a variable of a simdesign object}
\usage{
setsim(nl, var) <- value
}
\arguments{
\item{nl}{nl object}

\item{var}{valid simdesign variable string}

\item{value}{valid value for the specified variable}
}
\description{
Setter function to set a variable of a simdesign object
}
\examples{
# Example for Wolf Sheep Predation model from NetLogo models library:
nl <- nl(nlversion = "6.0.3",
nlpath = "/home/user/NetLogo 6.0.3/",
modelpath = "/home/user/NetLogo 6.0.3/app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo",
jvmmem = 1024)

# Set simulation seeds
setsim(nl, "simseeds") <- c(123, 456, 789)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_lhs}
\alias{simdesign_lhs}
\title{Add a latin-hypercube simdesign to a nl object}
\usage{
simdesign_lhs(nl, samples, nseeds, precision)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{samples}{number of samples for the latin hypercube}

\item{nseeds}{number of seeds for this simulation design}

\item{precision}{number of digits for the decimal fraction of parameter values}
}
\value{
simdesign S4 class object
}
\description{
Add a latin-hypercube simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min, max and qfun (e.g. list(min=1, max=4, qfun="qunif")).

The latin hypercube simdesign creates a parameter matrix based on these defined distributions.
Finally, the function reports a simdesign object.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_lhs
nl@simdesign <- simdesign_lhs(nl=nl,
                               samples=100,
                               nseeds=3,
                               precision=3)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_distinct}
\alias{nl_distinct}
\title{Wolf Sheep model sample data: simdesign distinct}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_distinct
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a distinct simulation design.
It also contains the model outputs within the simdesign object.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_eval.R
\name{util_eval_variables_ff}
\alias{util_eval_variables_ff}
\title{Evaluate variables list of an experiment object for full-factorial simdesign}
\usage{
util_eval_variables_ff(nl)
}
\arguments{
\item{nl}{nl object}
}
\description{
Evaluate variables list of an experiment object for full-factorial simdesign
}
\details{
util_eval_variables_ff checks if the variables list of an experiment within a nl object has enough information to create a \link[nlrx]{simdesign_ff}.
It reports an error message if at least one variable does not have a defined sequence (min, max, step) or a vector of distinct values.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval_variables_constants.R
\name{eval_variables_constants}
\alias{eval_variables_constants}
\title{Evaluate variable validity}
\usage{
eval_variables_constants(nl)
}
\arguments{
\item{nl}{nl object}
}
\description{
Evaluate variables and constants defined in experiment
}
\details{
This function checks if the variables and constants that are defined in the
experiment are valid.
It loads the model code of the NetLogo model and checks if these variables
and constants really exist.
In case of nonvalid entries, the function throws an error message, indicating
which variables and constants are not valid.
Please note, that this function might fail if the supported modelpath does
not point to an existing nlogo file.
This might for example happen, if the modelpath is set up for a remote
cluster execution.
}
\examples{
\dontrun{
nl <- nl_lhs
eval_variables_constants(nl)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_eFast}
\alias{simdesign_eFast}
\title{Add an eFast simdesign to a nl object}
\usage{
simdesign_eFast(nl, samples, nseeds)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{samples}{number of samples for the eFast sensitivity analysis}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add an eFast simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min, max and qfun (e.g. list(min=1, max=4, qfun="qunif")).

The eFast simdesign uses the sensitivity package to set up a fast99 elementary effects sensitivity analysis, including a simobject of class fast99 and a input tibble for simulations.
For details on method specific sensitivity analysis function parameters see ?fast99
Finally, the function reports a simdesign object.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_eFast
nl@simdesign <- simdesign_eFast(nl=nl,
                                 samples=100,
                                 nseeds=1)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report_model_parameters.R
\name{report_model_parameters}
\alias{report_model_parameters}
\title{Report globals from a NetLogo model that is defined within a nl object}
\usage{
report_model_parameters(nl)
}
\arguments{
\item{nl}{nl object with a defined modelpath that points to a NetLogo model
(*.nlogo)}
}
\description{
Report globals from a NetLogo model that is defined within a nl
object
}
\details{
The function reads the NetLogo model file that is defined within the nl object
and reports all global parameters that are defined as widget elements on
the GUI of the NetLogo model.
Only globals that are found by this function are valid globals that can be
entered into the variables or constants vector of an experiment object.
}
\examples{
\dontrun{
nl <- nl_lhs
report_model_parameters(nl)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_GenAlg}
\alias{simdesign_GenAlg}
\title{Add a Genetic Algorithm simdesign to a nl object}
\usage{
simdesign_GenAlg(
  nl,
  popSize = 200,
  iters = 100,
  evalcrit = 1,
  elitism = NA,
  mutationChance = NA,
  nseeds = 1
)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{popSize}{population Size parameter for genetic algorithm}

\item{iters}{number of iterations for genetic algorithm function}

\item{evalcrit}{position of evaluation criterion within defined NetLogo metrics of nl experiment or a function that reports a single numeric value}

\item{elitism}{elitism rate of genetic algorithm function}

\item{mutationChance}{mutation rate of genetic algorithm function}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add a Genetic Algorithm simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min and max (e.g. list(min=1, max=4)).

The GenAlg simdesign generates a Genetic Algorithm experiment within the defined min and max parameter boundaries
that are defined in the variables field of the experiment object within the nl object.

The evalcrit reporter defines the evaluation criterion for the Genetic algorithm procedure.
There are two options to evaluate the fitness value of each iteration of the algorithm:
\enumerate{
\item Use a reporter that is defined within the experiment metrics vector.
You can just enter the position of that metric within the experiment metrics vector (e.g. 1 would use the first defined metric of the experiment to evaluate each iteration).
The algorithm automatically calculates the mean value of this reporter if evalticks is defined to measure multiple ticks during each simulation.
\item Use a self-defined evaluation function
You can define a function that post-processes NetLogo output to calculate an evaluation value. This function must accept the nl object as input and return one single numeric value.
The nl object that is then provided to the evaluation function will have results of the current iteration attached. The results can be accessed via the simoutput slot of the simdesign.
You can pass this function to evalcrit. It is then applied to the output of each iteration.
}

The function uses the genalg package to set up a Genetic Algorithm function.
For details on the genalg function parameters see ?genalg::rbga
Finally, the function reports a simdesign object.

Genetic Algorithm simdesigns can only be executed using the \link[nlrx]{run_nl_dyn} function instead of \link[nlrx]{run_nl_all} or \link[nlrx]{run_nl_one}.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_lhs

# Example 1: Using a metric from the experiment metrics vector for evaluation:
nl@simdesign <- simdesign_GenAlg(nl=nl,
                                  evalcrit=1,
                                  nseeds=1)

# Example 2: Using a self-defined evaluation function
# For demonstration we define a simple function that calculates
# the maximum value of count sheep output.
critfun <- function(nl) {
results <- nl@simdesign@simoutput
crit <- as.integer(max(results$`count sheep`))
return(crit)
}

nl@simdesign <- simdesign_GenAlg(nl=nl,
                                  evalcrit=critfun,
                                  nseeds=1)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_setget.R
\name{getnl}
\alias{getnl}
\title{Getter function to get a variable of a nl object}
\usage{
getnl(nl, var)
}
\arguments{
\item{nl}{nl object}

\item{var}{valid nl variable string}
}
\description{
Getter function to get a variable of a nl object
}
\examples{

# Example for Wolf Sheep Predation model from NetLogo models library:
nl <- nl(nlversion = "6.0.3",
nlpath = "/home/user/NetLogo 6.0.3/",
modelpath = "/home/user/NetLogo 6.0.3/app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo",
jvmmem = 1024)

# get NetLogo version
getnl(nl, "nlversion")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_nl.R
\name{print.nl}
\alias{print.nl}
\title{Print content of nl object}
\usage{
\method{print}{nl}(x, ...)
}
\arguments{
\item{x}{nl object to print}

\item{...}{further arguments passed to or from other methods}
}
\description{
Print content of nl object and embedded experiment and simdesign objects to console
}
\details{
Print content of the provided nl object in a readable format.
}
\examples{

print(nl_lhs)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_stuff.R
\name{util_generate_seeds}
\alias{util_generate_seeds}
\title{Generate a vector of random seeds}
\usage{
util_generate_seeds(nseeds)
}
\arguments{
\item{nseeds}{desired length of the random seeds vector}
}
\description{
Generate a vector of random seeds
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_soboljansen}
\alias{simdesign_soboljansen}
\title{Add a soboljansen simdesign to a nl object}
\usage{
simdesign_soboljansen(nl, samples, sobolnboot, sobolconf, nseeds, precision)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{samples}{number of samples for the sobol sensitivity analysis}

\item{sobolnboot}{number of bootstrap replicates of the sobol sensitivity analysis}

\item{sobolconf}{the confidence level for bootstrap confidence intervals}

\item{nseeds}{number of seeds for this simulation design}

\item{precision}{number of digits for the decimal fraction of parameter values}
}
\value{
simdesign S4 class object
}
\description{
Add a soboljansen simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a numeric distribution with min, max and qfun (e.g. list(min=1, max=4, qfun="qunif")).

The soboljansen simdesign uses the sensitivity package to set up a soboljansen sensitivity analysis, including a simobject of class sobol and a input tibble for simulations.
For details on method specific sensitivity analysis function parameters see ?soboljansen
Finally, the function reports a simdesign object.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_soboljansen
nl@simdesign <- simdesign_soboljansen(nl=nl,
samples=1000,
sobolnboot=100,
sobolconf=0.95,
nseeds=3,
precision=3)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_nl.R
\name{util_print.nl}
\alias{util_print.nl}
\title{Print nl object content}
\usage{
util_print.nl(x, ...)
}
\arguments{
\item{x}{nl object}

\item{...}{further arguments passed to or from other methods}
}
\description{
Print nl object content
}
\details{
Print content of nl object in a nice formatted overview to console
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlrx-package.R
\docType{package}
\name{nlrx-package}
\alias{nlrx}
\alias{nlrx-package}
\title{nlrx: A package for running NetLogo simulations from R.}
\description{
The nlrx package provides tools to setup NetLogo simulations in R.
It uses a similar structure as NetLogos Behavior Space but offers more flexibility and additional tools for running sensitivity analyses.
}
\details{
Get started:

General information that is needed to run NetLogo simulations remotely, such as path to the NetLogo installation folder is stored within a \code{nl} class object.
Nested within this \code{nl} class are the classes \code{experiment} and \code{simdesign}. The \code{experiment} class stores all experiment specifications.
After attaching a valid experiment, a \code{simdesign} class object can be attached to the \code{nl} class object, by using one of the simdesign helper functions.
These helper functions create different parameter input matrices from the experiment variable definitions that can then be executed by the \code{run_nl_one()} and \code{run_nl_all()} functions.
The nested design allows to store everything related to the experiment within one R object.
Additionally, different simdesign helper functions can be applied to the same \code{nl} object in order to repeat the same experiment with different parameter exploration methods (simdesigns).

Step by step application example

The "Wolf Sheep Predation" model from the NetLogo models library is used to present a basic example on how to setup and run NetLogo model simulations from R.

Step 1: Create a nl object:

The nl object holds all information on the NetLogo version, a path to the NetLogo directory with the defined version, a path to the model file, and the desired memory for the java virtual machine.
Depending on the operation system, paths to NetLogo and the model need to be adjusted.\preformatted{library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.3")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo")
outpath <- file.path("/home/out")
nl <- nl(nlversion = "6.0.3",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)
}

Step 2: Attach an experiment

The experiment object is organized in a similar fashion as NetLogo Behavior Space experiments.
It contains all information that is needed to generate a simulation parameter matrix and to execute the NetLogo simulations.
Details on the specific slots of the experiment class can be found in the package documentation (\code{?experiment}) and the "Advanced configuration" vignette.\preformatted{nl@experiment <- experiment(expname="wolf-sheep",
                             outpath=outpath,
                             repetition=1,
                             tickmetrics="true",
                             idsetup="setup",
                             idgo="go",
                             runtime=50,
                             evalticks=seq(40,50),
                             metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                             variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                              'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                             constants = list("model-version" = "\\"sheep-wolves-grass\\"",
                                              "grass-regrowth-time" = 30,
                                              "sheep-gain-from-food" = 4,
                                              "wolf-gain-from-food" = 20,
                                              "sheep-reproduce" = 4,
                                              "wolf-reproduce" = 5,
                                              "show-energy?" = "false"))
}

Step 3: Attach a simulation design

While the experiment defines the variables and specifications of the model, the simulation design creates a parameter input table based on these model specifications and the chosen simulation design method.
nlrx provides a bunch of different simulation designs, such as full-factorial, latin-hypercube, sobol, morris and eFast (see "Simdesign Examples" vignette for more information on simdesigns).
All simdesign helper functions need a properly defined nl object with a valid experiment design.
Each simdesign helper also allows to define a number of random seeds that are randomly generated and can be used to execute repeated simulations of the same parameter matrix with different random-seeds (see "Advanced configuration" vignette for more information on random-seed and repetition management).
A simulation design is attached to a nl object by using one of the simdesign helper functions:\preformatted{nl@simdesign <- simdesign_lhs(nl=nl,
                               samples=100,
                               nseeds=3,
                               precision=3)
}

Step 4: Run simulations

All information that is needed to run the simulations is now stored within the nl object.
The \code{run_nl_one()} function allows to run one specific simulation from the siminput parameter table.
The \code{run_nl_all()} function runs a loop over all simseeds and rows of the parameter input table siminput.
The loops are constructed in a way that allows easy parallelisation, either locally or on remote HPC machines (see "Advanced configuration" vignette for more information on parallelisation).\preformatted{results <- run_nl_all(nl = nl)
}

Step 5: Attach results to nl and run analysis

nlrx provides method specific analysis functions for each simulation design.
Depending on the chosen design, the function reports a tibble with aggregated results or sensitivity indices.
In order to run the \code{analyze_nl()} function, the simulation output has to be attached to the nl object first.
After attaching the simulation results, these can also be written to the defined outpath of the experiment object.\preformatted{# Attach results to nl object:
setsim(nl, "simoutput") <- results
# Write output to outpath of experiment within nl
write_simoutput(nl)
# Do further analysis:
analyze_nl(nl)
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/nlrx/}
  \item \url{https://github.com/ropensci/nlrx/}
  \item Report bugs at \url{https://github.com/ropensci/nlrx/issues/}
}

}
\author{
Jan Salecker \email{jan.salecker@posteo.de}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_sobol}
\alias{nl_sobol}
\title{Wolf Sheep model sample data: simdesign sobol}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_sobol
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a sobol simulation design.
It also contains the model outputs within the simdesign object.
Further analysis output can be derived by submitting the dataset to analyze_nl().
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_nl.R
\name{analyze_simple}
\alias{analyze_simple}
\title{Analyze NetLogo simulation output of simdesign simple}
\usage{
analyze_simple(nl, metrics, funs)
}
\arguments{
\item{nl}{nl object}

\item{metrics}{vector of strings defining metric columns for evaluation. Defaults to metrics of the experiment within the nl object}

\item{funs}{list with the summary metrics for the sensitivity results}
}
\description{
Analyze NetLogo simulation output of simdesign simple
}
\details{
The simdesign_simple analysis functions is not yet supported and will only print a warning message.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nl_to_points.R
\name{nl_to_points}
\alias{nl_to_points}
\title{Get spatial data from metrics.turtles output}
\usage{
nl_to_points(nl, coords)
}
\arguments{
\item{nl}{nl object}

\item{coords}{nl object}
}
\value{
tibble with spatial data objects
}
\description{
Turn turtle metrics from NetLogo in spatial data objects
}
\details{
Converts measured metrics.turtles into spatial sf point objects.
In order to so, a pair of turtle coordinates needs to be measured.
Any additional metrics will be stored as properties of the spatial points.
Because turtle coordinates in NetLogo can be measured in two formats,
pxcor/pycor or xcor/ycor coordinates, the type of coordinate used for
transformation to spatial objects need to be defined, using the parameter
coords: "px" for pxcor/pycor coordinates, "x" for xcor/ycor coordinates.

In order to use this function, simulation results need to be attached to
the nl object first.
}
\examples{

# Load nl object (with spatial output data already attached) from test data:
nl <- nl_spatial

# Convert turtle metrics (pxcor/pycor) to spatial point objects:
results.sf <- nl_to_points(nl, coords="px")


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl.R
\name{util_gather_results}
\alias{util_gather_results}
\title{Load output file from simulations}
\usage{
util_gather_results(nl, outfile, seed, siminputrow)
}
\arguments{
\item{nl}{nl object}

\item{outfile}{file location of output results}

\item{seed}{model random-seed}

\item{siminputrow}{current row of siminput tibble}
}
\description{
Load output file from simulations
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nl_to_graph.R
\name{nl_to_graph}
\alias{nl_to_graph}
\title{Generate igraph objects from measured turtles and links metrics}
\usage{
nl_to_graph(nl)
}
\arguments{
\item{nl}{nl object

Generate igraph objects from measured turtles and links metrics.
Output has to be attached to the simdesign first with simoutput(nl) <- results
I graph objects are created automatically for each combination of random-seed, siminputrow and step.
An additional column with igraph objects is attached to the original output results tibble of the nl object.
In order to generate igraph objects some metrics are mandatory:
The metrics.turtles slot of the experiment must contain "who" numbers (see example experiment).
Additional turtle metrics will be stored as properties of the igraph vertices.
The metrics.links slot of the experiment must contain "who" numbers of link end1 and end2 (see example experiment).
Additional link metrics will be stored as properties of the igraph edges.}
}
\description{
Generate igraph objects from measured turtles and links metrics
}
\examples{
\dontrun{
## Example running the Giant Component model from the NetLogo models library:
library(nlrx)
library(igraph)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.0.4")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Networks/Giant Component.nlogo")
outpath <- file.path("C:/out")
# Unix default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("/home/NetLogo 6.0.4")
modelpath <- file.path(netlogopath, "app/models/Sample Models/Networks/Giant Component.nlogo")
outpath <- file.path("/home/out")

nl <- nl(nlversion = "6.0.4",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

nl@experiment <- experiment(expname="networks",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="false",
                            idsetup="setup",
                            idgo="go",
                            runtime=50,
                            metrics.turtles = list("turtles" = c("who", "color")),
                            metrics.links = list("links" = c("[who] of end1", "[who] of end2")),
                            constants = list("num-nodes" = 80,
                                             "layout?" = "true"))

nl@simdesign <- simdesign_simple(nl, 1)
nl@simdesign@simoutput <- run_nl_all(nl)
nl.graph <- nl_to_graph(nl)

## Extract graph of tick 1:
nl.graph.i <- nl.graph$spatial.links[[1]]
## Set vertex colors by measured color variable:
vcols <- c("7" = "grey", "15" = "red")
V(nl.graph.i)$color <- vcols[as.character(V(nl.graph.i)$color)]
## Set edge colors by measured link breed:
ecols <- c("links" = "black")
E(nl.graph.i)$color <- ecols[E(nl.graph.i)$breed]

## Plot:
plot.igraph(nl.graph.i, vertex.size=8, vertex.label=NA, edge.arrow.size=0.2)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_simoutput.R
\name{write_simoutput}
\alias{write_simoutput}
\title{Write attached NetLogo simulation output to file}
\usage{
write_simoutput(nl, outpath = NA)
}
\arguments{
\item{nl}{nl object}

\item{outpath}{optional path to directory where output is written

Write NetLogo simulation output to a csv file in the directory outpath of the nl object
Output has to be attached to the simdesign first with simoutput(nl) <- results
The outpath argument can be optionally used to write output to a different directory than the defined outpath of the nl object.}
}
\description{
Write attached NetLogo simulation output to file
}
\examples{

# Load nl object including output data from testdata
nl <- nl_lhs

# Write output to outpath directory
write_simoutput(nl, outpath=tempdir())

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_stuff.R
\name{util_get_os}
\alias{util_get_os}
\title{Identify and report the current OS}
\usage{
util_get_os()
}
\description{
Identify and report the current OS
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_eval.R
\name{util_eval_simdesign}
\alias{util_eval_simdesign}
\title{Evaluate all slots of a simdesign object}
\usage{
util_eval_simdesign(nl)
}
\arguments{
\item{nl}{nl object}
}
\description{
Evaluate all slots of a simdesign object
}
\details{
util_eval_simdesign checks if the information stored within the simdesign slots are valid.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_nlrx.R
\name{test_nlrx}
\alias{test_nlrx}
\title{Test if nlrx runs on the local system}
\usage{
test_nlrx(nlpath, nlversion)
}
\arguments{
\item{nlpath}{Provide a path to a netlogo folder}

\item{nlversion}{Matching version string of the provided NetLogo folder (e.g. "6.1.1")}
}
\description{
Runs a short test if nlrx runs on the local system
}
\details{
Runs a short test if nlrx runs on the local system. Reports TRUE if successful!
}
\examples{
\dontrun{
test_nlrx(nlpath="/Users/xyz/netlogo/NetLogo 6.1.1", nlversion="6.1.1")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl_dyn.R
\name{util_run_nl_dyn_GenAlg_fn}
\alias{util_run_nl_dyn_GenAlg_fn}
\alias{util_run_nl_dyn_ABCmcmc_fn}
\title{Genetic Algorithm run simulation function}
\usage{
util_run_nl_dyn_GenAlg_fn(
  param,
  nl,
  evalcrit,
  seed,
  cleanup.csv,
  cleanup.xml,
  cleanup.bat
)

util_run_nl_dyn_ABCmcmc_fn(param)
}
\arguments{
\item{param}{vector of model parameters passed from ABC_mcmc function. If use_seeds = TRUE, the first element of this vector is a random seed}

\item{nl}{nl object}

\item{evalcrit}{evaluation criterion for simulated annealing}

\item{seed}{current model seed}

\item{cleanup.csv}{TRUE/FALSE, if TRUE temporary created csv output files will be deleted after gathering results.}

\item{cleanup.xml}{TRUE/FALSE, if TRUE temporary created xml output files will be deleted after gathering results.}

\item{cleanup.bat}{TRUE/FALSE, if TRUE temporary created bat/sh output files will be deleted after gathering results.}
}
\description{
Genetic Algorithm run simulation function

Genetic Algorithm run simulation function
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export_nl.R
\name{export_nl}
\alias{export_nl}
\title{Export NetLogo Experiment}
\usage{
export_nl(nl, path = dirname(getnl(nl, "modelpath")), tarfile)
}
\arguments{
\item{nl}{nl object}

\item{path}{Path to folder that contains files to run NetLogo experiment}

\item{tarfile}{Path to folder where the experiments gets stored as zip file}
}
\value{
The status value returned by the external command, invisibly.
}
\description{
Export NetLogo Experiment as zip file
}
\details{
Exports your folder that contains data and scripts for NetLogo + nlrx analyses
as a zip file. Furthermore, \code{export_nl} takes your nl object and saves it in the
zipped folder. This enables other person to run the same experiment as you.
}
\examples{
\dontrun{

# Load nl object from testdata:
nl <- nl_lhs
path <- getwd() # adjust path to your needs, path should point to a directory with model data
outfile <- tempfile(fileext = ".zip") # adjust file path to your needs
export_nl(nl, path = path, tarfile = outfile)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_nl.R
\name{analyze_lhs}
\alias{analyze_lhs}
\title{Analyze NetLogo simulation output of simdesign latin-hypercube}
\usage{
analyze_lhs(nl, metrics, funs)
}
\arguments{
\item{nl}{nl object}

\item{metrics}{vector of strings defining metric columns for evaluation. Defaults to metrics of the experiment within the nl object}

\item{funs}{list with the summary metrics for the sensitivity results}
}
\description{
Analyze NetLogo simulation output of simdesign latin-hypercube
}
\details{
The function calculates aggregated output metrics by dropping random seeds and aggregating values with the provided functions.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_eval.R
\name{util_eval_variables_sa}
\alias{util_eval_variables_sa}
\title{Evaluate variables list of an experiment object for sensitivity analysis simdesigns}
\usage{
util_eval_variables_sa(nl)
}
\arguments{
\item{nl}{nl object}
}
\description{
Evaluate variables list of an experiment object for sensitivity analysis simdesigns
}
\details{
util_eval_variables_sa checks if the variables list of an experiment within a nl object has enough information to create a sensitivity analysis simdesign.
It reports an error message if at least one variable does not have a defined distribution (min, max, qfun).
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_setget.R
\name{setnl<-}
\alias{setnl<-}
\alias{setnl}
\title{Setter function to set a variable of a nl object}
\usage{
setnl(nl, var) <- value
}
\arguments{
\item{nl}{nl object}

\item{var}{valid nl variable string}

\item{value}{valid value for the specified variable}
}
\description{
Setter function to set a variable of a nl object
}
\examples{
# Example for Wolf Sheep Predation model from NetLogo models library:
nl <- nl(
nlpath = "/home/user/NetLogo 6.0.3/",
modelpath = "/home/user/NetLogo 6.0.3/app/models/Sample Models/Biology/Wolf Sheep Predation.nlogo",
jvmmem = 1024)

# set NetLogo version
setnl(nl, "nlversion") <- "6.0.3"

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_nl.R
\name{analyze_soboljansen}
\alias{analyze_soboljansen}
\title{Analyze NetLogo simulation output of simdesign soboljansen}
\usage{
analyze_soboljansen(nl, metrics, funs)
}
\arguments{
\item{nl}{nl object}

\item{metrics}{vector of strings defining metric columns for evaluation. Defaults to metrics of the experiment within the nl object}

\item{funs}{list with the summary metrics for the sensitivity results}
}
\description{
Analyze NetLogo simulation output of simdesign soboljansen
}
\details{
The function calculates sobol sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_ff}
\alias{simdesign_ff}
\title{Add a full-factorial simdesign to a nl object}
\usage{
simdesign_ff(nl, nseeds)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add a full-factorial simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.

Variables in the experiment variable list need to provide a vector of distinct values (e.g. list(values=c(1,2,3,4)).
Or a sequence definition with min, max and step (e.g. list=(min=1, max=4, step=1)).
If both (values and sequence) are defined, the full-factorial design gives priority to the values.

The full-factorial simdesign uses these defined parameter ranges within the nl object.
A full-factorial matrix of all parameter combinations is created as input tibble for the simdesign.
Finally, the function reports a simdesign object.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_ff
nl@simdesign <- simdesign_ff(nl = nl, nseeds = 3)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_eval.R
\name{util_eval_experiment}
\alias{util_eval_experiment}
\title{Evaluate all slots of an experiment object}
\usage{
util_eval_experiment(nl)
}
\arguments{
\item{nl}{nl object}
}
\description{
Evaluate all slots of an experiment object
}
\details{
util_eval_experiment checks if the information stored within the experiment slots are valid.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_nl.R
\name{util_print.experiment}
\alias{util_print.experiment}
\title{Print experiment object content}
\usage{
util_print.experiment(x, ...)
}
\arguments{
\item{x}{experiment object}

\item{...}{further arguments passed to or from other methods}
}
\description{
Print experiment object content
}
\details{
Print content of experiment object in a nice formatted overview to console
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl.R
\name{util_cleanup}
\alias{util_cleanup}
\title{Delete temporary files}
\usage{
util_cleanup(
  nl,
  cleanup.csv = TRUE,
  cleanup.xml = TRUE,
  cleanup.bat = TRUE,
  cleanup.files
)
}
\arguments{
\item{nl}{nl object}

\item{cleanup.csv}{TRUE/FALSE, if TRUE temporary created csv output files will be deleted after gathering results.}

\item{cleanup.xml}{TRUE/FALSE, if TRUE temporary created xml output files will be deleted after gathering results.}

\item{cleanup.bat}{TRUE/FALSE, if TRUE temporary created bat/sh output files will be deleted after gathering results.}

\item{cleanup.files}{vector with paths to temporary created files (csv, xml, bat)}
}
\description{
Delete temporary files
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_nl.R
\name{run_nl_all}
\alias{run_nl_all}
\title{Execute all NetLogo simulations from a nl object}
\usage{
run_nl_all(
  nl,
  split = 1,
  cleanup.csv = TRUE,
  cleanup.xml = TRUE,
  cleanup.bat = TRUE
)
}
\arguments{
\item{nl}{nl object}

\item{split}{number of parts the job should be split into}

\item{cleanup.csv}{TRUE/FALSE, if TRUE temporary created csv output files will be deleted after gathering results.}

\item{cleanup.xml}{TRUE/FALSE, if TRUE temporary created xml output files will be deleted after gathering results.}

\item{cleanup.bat}{TRUE/FALSE, if TRUE temporary created bat/sh output files will be deleted after gathering results.}
}
\value{
tibble with simulation output results
}
\description{
Execute all NetLogo simulations from a nl object with a defined experiment and simdesign
}
\details{
run_nl_all executes all simulations of the specified NetLogo model within the provided nl object.
The function loops over all random seeds and all rows of the siminput table of the simdesign of nl.
The loops are created by calling \link[furrr]{future_map_dfr}, which allows running the function either locally or on remote HPC machines.
The logical cleanup variables can be set to FALSE to preserve temporary generated output files (e.g. for debugging).
cleanup.csv deletes/keeps the temporary generated model output files from each run.
cleanup.xml deletes/keeps the temporary generated experiment xml files from each run.
cleanup.bat deletes/keeps the temporary generated batch/sh commandline files from each run.

When using run_nl_all in a parallelized environment (e.g. by setting up a future plan using the future package),
the outer loop of this function (random seeds) creates jobs that are distributed to available cores of the current machine.
The inner loop (siminputrows) distributes simulation tasks to these cores.
However, it might be advantageous to split up large jobs into smaller jobs for example to reduce the total runtime of each job.
This can be done using the split parameter. If split is > 1 the siminput matrix is split into smaller parts.
Jobs are created for each combination of part and random seed.
If the split parameter is set such that the siminput matrix can not be splitted into equal parts, the procedure will stop and throw an error message.
\subsection{Debugging "Temporary simulation output file not found" error message:}{

Whenever this error message appears it means that the simulation did not produce any output.
Two main reasons can lead to this problem, either the simulation did not even start or the simulation crashed during runtime.
Both can happen for several reasons and here are some hints for debugging this:
\enumerate{
\item Missing software:
Make sure that java is installed and available from the terminal (java -version).
Make sure that NetLogo is installed and available from the terminal.
\item Wrong path definitions:
Make sure your nlpath points to a folder containing NetLogo.
Make sure your modelpath points to a *.nlogo model file.
Make sure that the nlversion within your nl object matches the NetLogo version of your nlpath.
Use the convenience function of nlrx for checking your nl object (print(nl), eval_variables_constants(nl)).
\item Temporary files cleanup:
Due to automatic temp file cleanup on unix systems temporary output might be deleted.
Try reassigning the default temp folder for this R session (the unixtools package has a neat function).
\item NetLogo runtime crashes:
It can happen that your NetLogo model started but failed to produce output because of a NetLogo runtime error.
Make sure your model is working correctly or track progress using print statements.
Sometimes the java virtual machine crashes due to memory constraints.
}
}
}
\examples{
\dontrun{

# Load nl object from test data:
nl <- nl_lhs

# Execute all simulations from an nl object with properly attached simdesign.
results <- run_nl_all(nl)

# Run in parallel on local machine:
library(future)
plan(multisession)
results <- run_nl_all(nl)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_nl.R
\name{run_nl_one}
\alias{run_nl_one}
\title{Execute one NetLogo simulation from a nl object}
\usage{
run_nl_one(
  nl,
  seed,
  siminputrow,
  cleanup.csv = TRUE,
  cleanup.xml = TRUE,
  cleanup.bat = TRUE,
  writeRDS = FALSE
)
}
\arguments{
\item{nl}{nl object}

\item{seed}{a random seed for the NetLogo simulation}

\item{siminputrow}{rownumber of the input tibble within the attached simdesign object that should be executed}

\item{cleanup.csv}{TRUE/FALSE, if TRUE temporary created csv output files will be deleted after gathering results.}

\item{cleanup.xml}{TRUE/FALSE, if TRUE temporary created xml output files will be deleted after gathering results.}

\item{cleanup.bat}{TRUE/FALSE, if TRUE temporary created bat/sh output files will be deleted after gathering results.}

\item{writeRDS}{TRUE/FALSE, if TRUE an rds file with the simulation results will be written to the defined outpath folder of the experiment within the nl object.}
}
\value{
tibble with simulation output results
}
\description{
Execute one NetLogo simulation from a nl object with a defined experiment and simdesign
}
\details{
run_nl_one executes one simulation of the specified NetLogo model within the provided nl object.
The random seed is set within the NetLogo model to control stochasticity.
The siminputrow number defines which row of the input data tibble within the simdesign object of the provided nl object is executed.
The logical cleanup variables can be set to FALSE to preserve temporary generated output files (e.g. for debugging).
cleanup.csv deletes/keeps the temporary generated model output files from each run.
cleanup.xml deletes/keeps the temporary generated experiment xml files from each run.
cleanup.bat deletes/keeps the temporary generated batch/sh commandline files from each run.

This function can be used to run single simulations of a NetLogo model.
}
\examples{
\dontrun{

# Load nl object from test data:
nl <- nl_lhs

# Run one simulation:
results <- run_nl_one(nl = nl,
                      seed = getsim(nl, "simseeds")[1],
                      siminputrow = 1)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_runnl_dyn.R
\name{util_run_nl_dyn_ABCmcmc}
\alias{util_run_nl_dyn_ABCmcmc}
\title{ABCmcmc call simulations function}
\usage{
util_run_nl_dyn_ABCmcmc(
  nl,
  seed,
  cleanup.csv = TRUE,
  cleanup.xml = TRUE,
  cleanup.bat = TRUE
)
}
\arguments{
\item{nl}{nl object}

\item{seed}{current model seed}

\item{cleanup.csv}{TRUE/FALSE, if TRUE temporary created csv output files will be deleted after gathering results.}

\item{cleanup.xml}{TRUE/FALSE, if TRUE temporary created xml output files will be deleted after gathering results.}

\item{cleanup.bat}{TRUE/FALSE, if TRUE temporary created bat/sh output files will be deleted after gathering results.}
}
\description{
ABCmcmc call simulations function
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constr.R
\name{simdesign}
\alias{simdesign}
\title{Construct a new simdesign object}
\usage{
simdesign(
  simmethod = character(),
  siminput = tibble::tibble(),
  simobject = list(),
  simseeds = NA_integer_,
  simoutput = tibble::tibble(),
  ...
)
}
\arguments{
\item{simmethod}{character string defining the method of the simulation design}

\item{siminput}{tibble providing input parameterisations for the NetLogo model (cols=parameter, rows=runs)}

\item{simobject}{used for some methods to store additional information (sobol, morris, eFast)}

\item{simseeds}{a vector or model random seeds}

\item{simoutput}{tibble containing model results}

\item{...}{...}
}
\value{
simdesign S4 class object
}
\description{
Construct a new simdesign object
}
\details{
The simulation design class holds information on the input parameter design of model simulations.
It also stores information that is needed to run method specific analysis functions.
The simseeds can be used to run all model simulations that are defined within the siminput tibble several times with changing random-seeds.
While it is possible to add simdesign directly with this function, we suggest to use our simdesign_helper functions.
A simulation design can be attached to a nl object by using one of these simdesign_helper functions on an already defined \link[nlrx]{nl}
object with a valid \link[nlrx]{experiment}.
All simdesign helpers use the defined constants and variables of the experiment to create the siminput tibble.
NetLogo parameters that are not defined in constants or variables will be set with their default value from the NetLogo interface.

Currently, following simdesign_helper functions are provided:

\link[nlrx]{simdesign_simple}

The simple simdesign only uses defined constants and reports a parameter matrix with only one parameterization.
To setup a simple simdesign, no variables have to be defined.

\link[nlrx]{simdesign_distinct}

The distinct simdesign can be used to run distinct parameter combinations.
To setup a distinct simdesign, vectors of values need to be defined for each variable.
These vectors must have the same number of elements across all variables.
The first simulation run consist of all 1st elements of these variable vectors; the second run uses all 2nd values, and so on.

\link[nlrx]{simdesign_ff}

The full factorial simdesign creates a full-factorial parameter matrix with all possible combinations of parameter values.
To setup a full-factorial simdesign, vectors of values need to be defined for each variable.
Alternatively, a sequence can be defined by setting min, max and step.
However, if both (values and min, max, step) are defined, the values vector is prioritized.

\link[nlrx]{simdesign_lhs}

The latin hypercube simdesign creates a Latin Hypercube sampling parameter matrix.
The method can be used to generate a near-random sample of parameter values from the defined parameter distributions.
More Details on Latin Hypercube Sampling can be found in \href{https://www.tandfonline.com/doi/abs/10.1080/00401706.1979.10489755}{McKay 1979}.
nlrx uses the \href{https://CRAN.R-project.org/package=lhs/index.html}{lhs} package to generate the Latin Hypercube parameter matrix.
To setup a latin hypercube sampling simdesign, variable distributions need to be defined (min, max, qfun).

Sensitivity Analyses: \link[nlrx]{simdesign_sobol}, \link[nlrx]{simdesign_sobol2007}, \link[nlrx]{simdesign_soboljansen}, \link[nlrx]{simdesign_morris}, \link[nlrx]{simdesign_eFast}

Sensitivity analyses are useful to estimate the importance of model parameters and to scan the parameter space in an efficient way.
nlrx uses the \href{https://CRAN.R-project.org/package=sensitivity/index.html}{sensitivity} package to setup sensitivity analysis parameter matrices.
All supported sensitivity analysis simdesigns can be used to calculate sensitivity indices for each parameter-output combination.
These indices can be calculated by using the \link[nlrx]{analyze_nl} function after attaching the simulation results to the nl object.
To setup sensitivity analysis simdesigns, variable distributions (min, max, qfun) need to be defined.

Optimization techniques: \link[nlrx]{simdesign_GenSA}, \link[nlrx]{simdesign_GenAlg}

Optimization techniques are a powerful tool to search the parameter space for specific solutions.
Both approaches try to minimize a specified model output reporter by systematically (genetic algorithm, utilizing the \href{https://CRAN.R-project.org/package=genalg/index.html}{genalg} package) or randomly (simulated annealing, utilizing the \href{https://CRAN.R-project.org/package=GenSA/index.html}{genSA} package) changing the model parameters within the allowed ranges.
To setup optimization simdesigns, variable ranges (min, max) need to be defined.
Optimization simdesigns can only be executed using the \link[nlrx]{run_nl_dyn} function instead of \link[nlrx]{run_nl_all} or \link[nlrx]{run_nl_one}.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load nl objects from test data.

# Simdesign examples for Wolf Sheep Predation model from NetLogo models library:

nl <- nl_simple
nl@simdesign <- simdesign_simple(nl = nl,
                                 nseeds = 3)

nl <- nl_distinct
nl@simdesign <- simdesign_distinct(nl = nl,
                                   nseeds = 3)

nl <- nl_ff
nl@simdesign <- simdesign_ff(nl = nl,
                             nseeds = 3)

nl <- nl_lhs
nl@simdesign <- simdesign_lhs(nl=nl,
                              samples=100,
                              nseeds=3,
                              precision=3)

nl <- nl_sobol
nl@simdesign <- simdesign_sobol(nl=nl,
                                samples=200,
                                sobolorder=2,
                                sobolnboot=20,
                                sobolconf=0.95,
                                nseeds=3,
                                precision=3)

nl <- nl_sobol2007
nl@simdesign <- simdesign_sobol2007(nl=nl,
                                    samples=200,
                                    sobolnboot=20,
                                    sobolconf=0.95,
                                    nseeds=3,
                                    precision=3)

nl <- nl_soboljansen
nl@simdesign <- simdesign_soboljansen(nl=nl,
                                      samples=200,
                                      sobolnboot=20,
                                      sobolconf=0.95,
                                      nseeds=3,
                                      precision=3)

nl <- nl_morris
nl@simdesign <- simdesign_morris(nl=nl,
                                 morristype="oat",
                                 morrislevels=4,
                                 morrisr=100,
                                 morrisgridjump=2,
                                 nseeds=3)

nl <- nl_eFast
nl@simdesign <- simdesign_eFast(nl=nl,
                                samples=100,
                                nseeds=3)

nl <- nl_lhs
nl@simdesign <- simdesign_GenAlg(nl=nl,
                                 popSize = 200,
                                 iters = 100,
                                 evalcrit = 1,
                                 elitism = NA,
                                 mutationChance = NA,
                                 nseeds = 1)

nl <- nl_lhs
nl@simdesign <- simdesign_GenSA(nl=nl,
                                par=NULL,
                                evalcrit=1,
                                control=list(max.time = 600),
                                nseeds=1)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_soboljansen}
\alias{nl_soboljansen}
\title{Wolf Sheep model sample data: simdesign soboljansen}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_soboljansen
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a soboljansen simulation design.
It also contains the model outputs within the simdesign object.
Further analysis output can be derived by submitting the dataset to analyze_nl().
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_nl.R
\name{analyze_morris}
\alias{analyze_morris}
\title{Analyze NetLogo simulation output of simdesign morris}
\usage{
analyze_morris(nl, metrics, funs)
}
\arguments{
\item{nl}{nl object}

\item{metrics}{vector of strings defining metric columns for evaluation. Defaults to metrics of the experiment within the nl object}

\item{funs}{list with the summary metrics for the sensitivity results}
}
\description{
Analyze NetLogo simulation output of simdesign morris
}
\details{
The function calculates morris sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_nl.R
\name{analyze_ff}
\alias{analyze_ff}
\title{Analyze NetLogo simulation output of simdesign full-factorial}
\usage{
analyze_ff(nl, metrics, funs)
}
\arguments{
\item{nl}{nl object}

\item{metrics}{vector of strings defining metric columns for evaluation. Defaults to metrics of the experiment within the nl object}

\item{funs}{list with the summary metrics for the sensitivity results}
}
\description{
Analyze NetLogo simulation output of simdesign full-factorial
}
\details{
The function calculates aggregated output metrics by dropping random seeds and aggregating values with the provided functions.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdesign_helper.R
\name{simdesign_simple}
\alias{simdesign_simple}
\title{Add a simple simdesign to a nl object}
\usage{
simdesign_simple(nl, nseeds)
}
\arguments{
\item{nl}{nl object with a defined experiment}

\item{nseeds}{number of seeds for this simulation design}
}
\value{
simdesign S4 class object
}
\description{
Add a simple simdesign to a nl object
}
\details{
This function creates a simdesign S4 class which can be added to a nl object.
The simple simdesign only uses model parameters that are defined in the constants field of the experiment object within the nl object.
Thus, the resulting input tibble of the simdesign has only one run with constant parameterisations.
This can be useful to run one simulation with a specific parameterset.
Finally, the function reports a simdesign object.
}
\examples{

# To attach a simdesign, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

nl <- nl_simple
nl@simdesign <- simdesign_simple(nl = nl, nseeds = 3)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_init.R
\name{init-classes}
\alias{init-classes}
\alias{.initClasses}
\title{initClasses}
\usage{
.initClasses()
}
\description{
empty
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_nl.R
\name{run_nl_dyn}
\alias{run_nl_dyn}
\title{Execute NetLogo simulation without pregenerated parametersets}
\usage{
run_nl_dyn(
  nl,
  seed,
  cleanup.csv = TRUE,
  cleanup.xml = TRUE,
  cleanup.bat = TRUE
)
}
\arguments{
\item{nl}{nl object}

\item{seed}{a random seed for the NetLogo simulation}

\item{cleanup.csv}{TRUE/FALSE, if TRUE temporary created csv output files will be deleted after gathering results.}

\item{cleanup.xml}{TRUE/FALSE, if TRUE temporary created xml output files will be deleted after gathering results.}

\item{cleanup.bat}{TRUE/FALSE, if TRUE temporary created bat/sh output files will be deleted after gathering results.}
}
\value{
simulation output results can be tibble, list, ...
}
\description{
Execute NetLogo simulation from a nl object with a defined experiment and simdesign but no pregenerated input parametersets
}
\details{
run_nl_dyn can be used for simdesigns where no predefined parametersets exist.
This is the case for dynamic designs, such as Simulated Annealing and Genetic Algorithms, where parametersets are dynamically generated, based on the output of previous simulations.
The logical cleanup variables can be set to FALSE to preserve temporary generated output files (e.g. for debugging).
cleanup.csv deletes/keeps the temporary generated model output files from each run.
cleanup.xml deletes/keeps the temporary generated experiment xml files from each run.
cleanup.bat deletes/keeps the temporary generated batch/sh commandline files from each run.
}
\examples{
\dontrun{

# Load nl object form test data:
nl <- nl_lhs

# Add genalg simdesign:
nl@simdesign <- simdesign_GenAlg(nl=nl,
                                  popSize = 200,
                                  iters = 100,
                                  evalcrit = 1,
                                  nseeds = 1)

# Run simulations:
results <- run_nl_dyn(nl)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constr.R
\name{experiment}
\alias{experiment}
\title{Construct a new experiment object}
\usage{
experiment(
  expname = "defaultexp",
  outpath = NA_character_,
  repetition = 1,
  tickmetrics = "true",
  idsetup = "setup",
  idgo = "go",
  idfinal = NA_character_,
  idrunnum = NA_character_,
  runtime = 1,
  evalticks = NA_integer_,
  stopcond = NA_character_,
  metrics = c("count turtles"),
  metrics.turtles = list(),
  metrics.patches = NA_character_,
  metrics.links = list(),
  variables = list(),
  constants = list(),
  ...
)
}
\arguments{
\item{expname}{A character string defining the name of the experiment, no whitespaces allowed}

\item{outpath}{Path to a directory where experiment output will be stored}

\item{repetition}{A number which gives the number of repetitions for each row of the simulation design input tibble}

\item{tickmetrics}{Character string "true" runs defined metrics on each simulation tick. "false" runs metrics only after simulation is finished}

\item{idsetup}{character string or vector of character strings, defining the name of the NetLogo setup procedure}

\item{idgo}{character string or vector of character strings, defining the name of the NetLogo go procedure}

\item{idfinal}{character string or vector of character strings, defining the name of NetLogo procedures that should be run after the last tick}

\item{idrunnum}{character string, defining the name of a NetLogo global that should be used to parse the current siminputrow during model executions which can then be used for self-written output.}

\item{runtime}{number of model ticks that should be run for each simulation}

\item{evalticks}{vector of tick numbers defining when measurements are taken. NA_integer_ to measure each tick}

\item{stopcond}{a NetLogo reporter that reports TRUE/FALSE. If it reports TRUE the current simulation run is stopped (e.g. "not any? turtles")}

\item{metrics}{vector of strings defining valid NetLogo reporters that are taken as output measurements (e.g. c("count turtles", "count patches"))}

\item{metrics.turtles}{a list with named vectors of strings defining valid turtles-own variables that are taken as output measurements (e.g. list("turtles" = c("who", "pxcor", "pycor", "color"))}

\item{metrics.patches}{vector of strings defining valid patches-own variables that are taken as output measurements (e.g. c("pxcor", "pycor", "pcolor"))}

\item{metrics.links}{a list with named vectors of strings defining valid links-own variables that are taken as output measurements (e.g. list("links" = c("end1", "end2")))}

\item{variables}{a nested list of variables that are changed within a simulation design. The name of each sublist item has to be a valid global of the defined NetLogo model. Depending on the desired simdesign each list item consist of a vector of values, a min value, a max value, a step value and a qfun (e.g. list("paramA" = list(values=c(0, 0.5, 1), min=0, max=1, step=0.1, qfun="qunif")))}

\item{constants}{a list of constants that are kept constant within a simulation design. The name of each list item has to be a valid global of the defined NetLogo model (e.g. list("pNUM" = 12, "pLOGIC"="TRUE", "pSTRING"="\"default\""))}

\item{...}{...}
}
\value{
experiment S4 class object
}
\description{
Construct a new experiment object
}
\details{
The experiment class stores all information related to the NetLogo simulation experiment.
The class holds all information that is typically entered into NetLogo Behavior Space experiments.
When setting up an experiment, it is usually attached to an already defined \link[nlrx]{nl} object (see examples).
After attaching an experiment, different simdesign helper functions can be used to attach a simdesign to the nl object \link[nlrx]{simdesign}.
The simdesign helper functions use the variable definitions from the experiment within the nl object to generate a parameter tibble for simulations.

\strong{The following class slots are obligatory to run an experiment:}

\emph{repetition}

In cases, where the random seed is controlled by nlrx simdesigns, repitition should be set to one as random seeds would not differ between simulations.
In cases, where the random seed is set within the NetLogo model, repitition can be increased to repeat the same parameterisation with different random seeds.

\emph{tickmetrics}

If "true", the defined output reporters are collected on each simulation tick that is defined in evalticks. If "false" measurements are taken only on the last tick.

\emph{idsetup, idgo}

These two class slots accept strings, or vectors of strings, defining NetLogo model procedures that should be executed for model setup (idsetup) and model execution (idgo).

\emph{runtime}

Defines the maximum number of simulation ticks that are executed.
Use 0 or NA_integer_ to execute simulations without predefined number of ticks.
Warning: This is only recommended in combination with a stop condition (see slot stopcond), or NetLogo models with a built-in stop condition.
Otherwise, simulations might get stuck in endless loops.

\strong{Depending on the simdesign, the following slots may be obligatory:}

\emph{metrics}

A vector of valid netlogo reporters that defines which measurements are taken.
The optimization simdesigns need at least one defined metrics reporter for fitness calculation of the optimization algorithm.

\emph{constants, variables}

These slots accept lists with NetLogo parameters that should be varied within a simdesign (variables) or should be kept constant (constants) for each simulation.
Any NetLogo parameter that is not entered in at least one of these two lists will be set up as constant with the default value from the NetLogo interface.
It is not possible to enter a NetLogo parameter in both lists (a warning message will appear when a simdesign is attached to such an experiment).
All simdesigns except \link[nlrx]{simdesign_simple} need defined variables for setting up a parameter matrix.
Variables can be defined as distinct values, value distributions or range with increment.
The information that is needed, depends on the chosen simdesign (details on variable definition requirements can be found in the helpfiles of each simdesign helper function).

\strong{All remaining slots are optional:}

\emph{expname}

A character string defining the name of the experiment, useful for documentation purposes. The string must not contain any whitespaces.

\emph{outpath}

A valid path to an existing directory. The directory is used by the \link[nlrx]{write_simoutput} function to store attached simulation results to disk in csv format.

\emph{idfinal}

A character string or vector of strings defining NetLogo procedures that are executed at the end of each simulation (e.g. cleanup or self-written output procedures).

\emph{idrunnum}

This slot can be used to transfer the current nlrx experiment name, random seed and runnumber (siminputrow) to NetLogo.
To use this functionality, a string input field widget needs to be created on the GUI of your NetLogo model.
The name of this widget can be entered into the "idrunnum" field of the experiment.
During simulations, the value of this widget is automatically updated with a generated string that contains the current nlrx experiment name, random seed and siminputrow ("expname_seed_siminputrow").
For self-written output In NetLogo, we suggest to include this global variable which allows referencing the self-written output files to the collected output of the nlrx simulations in R.

\emph{evalticks}

Only applied if tickmetrics = TRUE.
Evalticks may contain a vector of integers, defining the ticks for which the defined metrics will be measured.
Set evalticks to NA_integer_ to measure on every tick.

\emph{stopcond}

The stopcond slot can be used to define a stop condition by providing a string with valid NetLogo code that reports either true or false.
Each simulation will be stopped automatically, once the reporter reports true.

\emph{metrics.patches}

The metrics.patches slot accepts a vector of valid patches-own variables of the NetLogo model.
These patch variables are measured in addition to the defined metrics.
Results of these metrics will be nested inside the output results tibble of the simulations.
Please note that NetLogo models may contain a huge number of patches and output measurements of agent variables on each tick may need a lot of ressources.

\emph{metrics.turtles}

The metrics.turtles slot accepts a list with named vectors of valid turtle breed metrics.
Each name of a vector in this list defines a specified breed of the NetLogo model, whereas the vector defines the variables that are measured for this breed.
For example metrics.turtles = list("sheep"=c("color"), "wolves"=c("who")) - would measure the color of each sheep and the who number of each wolf agent.
To measure <turtles-own> variables for all turtles, use "turtles" = c(...).
Be aware, that NetLogo will produce runtime errors if you measure <breed-own> variables for agents that do not belong to this breed.
Please note that NetLogo models may contain a huge number of turtles and output measurements of agent variables on each tick may need a lot of ressources.

\emph{metrics.links}

The metrics.links slot accepts a list with named vectors of valid link breed metrics.
Each name of a vector in this list defines a specified link breed of the NetLogo model, whereas the vector defines the variables that are measured for this link breed.
For example metrics.links = list("linktype-a"=c("end1"), "linktype-b"=c("end2")) - would measure the start agent of each linktype-a link and the end agent of each linktype-b link.
To measure <links-own> variables for all links, use "links" = c(...).
Be aware, that NetLogo will produce runtime errors if you measure <link-breed-own> variables for agents that do not belong to this breed.
Please note that NetLogo models may contain a huge number of turtles and output measurements of agent variables on each tick may need a lot of ressources.
}
\examples{

# To attach an experiment, a nl object needs to be created first (see ?nl).
# For this example, we load a nl object from test data.

# Example for Wolf Sheep Predation model from NetLogo models library:
nl <- nl_simple
nl@experiment <- experiment(expname="wolf-sheep",
                             outpath="C:/out/",
                             repetition=1,
                             tickmetrics="true",
                             idsetup="setup",
                             idgo="go",
                             idfinal=NA_character_,
                             idrunnum=NA_character_,
                             runtime=50,
                             evalticks=seq(40,50),
                             stopcond="not any? turtles",
                             metrics=c("count sheep",
                                       "count wolves",
                                     "count patches with [pcolor = green]"),
                             metrics.turtles=list("turtles" = c("who",
                                               "pxcor",
                                               "pycor",
                                               "color")),
                             metrics.patches=c("pxcor", "pycor", "pcolor"),
                             variables = list('initial-number-sheep' =
                             list(min=50, max=150, step=10, qfun="qunif"),
                                              'initial-number-wolves' =
                             list(min=50, max=150, step=10, qfun="qunif")),
                             constants = list("model-version" =
                                              "\"sheep-wolves-grass\"",
                                              "grass-regrowth-time" = 30,
                                              "sheep-gain-from-food" = 4,
                                              "wolf-gain-from-food" = 20,
                                              "sheep-reproduce" = 4,
                                              "wolf-reproduce" = 5,
                                              "show-energy?" = "false"))


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{nl_sobol2007}
\alias{nl_sobol2007}
\title{Wolf Sheep model sample data: simdesign sobol2007}
\format{
nl object with defined experiment, simdesign and model output
}
\usage{
nl_sobol2007
}
\description{
The dataset contains a complete nl object.
The nl object has been used to setup and run a sobol2007 simulation design.
It also contains the model outputs within the simdesign object.
Further analysis output can be derived by submitting the dataset to analyze_nl().
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval_simoutput.R
\name{eval_simoutput}
\alias{eval_simoutput}
\title{Evaluate input/output integrity}
\usage{
eval_simoutput(nl)
}
\arguments{
\item{nl}{nl object with attached simulation output}
}
\description{
Evaluate input/output integrity
}
\details{
This function checks if the attached simulation output in the simoutput slot of the simdesign,
corresponds to the defined siminput matrix.

Warning messages are thrown if data is missing in the simoutput tibble.
Additionally, missing combinations of siminputrow and random seed for which no data was found can be reported as tibble.
Such a tibble can then be used directly to rerun missing combinations conveniently (see examples below)
}
\examples{
\dontrun{
# Check eval_simoutput for testdata nl_lhs:
nl <- nl_lhs
eval_simoutput(nl)

# Now remove one row of simoutput and check output:
nl <- nl_lhs
nl@simdesign@simoutput <- nl@simdesign@simoutput[-1,]
check <- eval_simoutput(nl)
check

# Rerun missing combinations within check tibble:
rerun <- purrr::map_dfr(seq(nrow(check)), function(x) {
  res <- run_nl_one(nl, siminputrow=check$siminputrow[x], seed=check$seed[x])
    return(res)
    }) \%>\%
      dplyr::bind_rows(., nl@simdesign@simoutput)


}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nl_to_raster.R
\name{nl_to_raster}
\alias{nl_to_raster}
\title{Get spatial data from metrics.patches output}
\usage{
nl_to_raster(nl)
}
\arguments{
\item{nl}{nl object}
}
\value{
tibble with spatial data objects
}
\description{
Turn patch metrics from NetLogo in spatial data objects
}
\details{
Converts measured metrics.patches into spatial raster objects.
In order to so, a patch coordinates needs to be measured (pxcor/pycor).
For each additional patch metric, a raster will be created using the
measured coordinates as x and y and the additional metric as z dimension.
In case of multiple measured metrics, a raster stack with one raster
for each metric will be reported.

In order to use this function, simulation results need to be attached to
the nl object first.
}
\examples{

# Load nl object (with spatial output data already attached) from test data:
nl <- nl_spatial

# Convert patch metrics to spatial raster objects:
results.raster <- nl_to_raster(nl)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nldoc_find_procedure_calls.R
\name{nldoc_find_procedure_calls}
\alias{nldoc_find_procedure_calls}
\title{Determine procedure calls}
\usage{
nldoc_find_procedure_calls(nlogocode)
}
\arguments{
\item{nlogocode}{vector of netlogo code strings}
}
\value{
tibble with procedure names and procedure calls
}
\description{
Determine procedure calls
}
\details{
The procedure searches netlogo code for procedure definitions and calls.
The information is stored within a tibble that can be further processed.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze_nl.R
\name{analyze_sobol}
\alias{analyze_sobol}
\title{Analyze NetLogo simulation output of simdesign sobol}
\usage{
analyze_sobol(nl, metrics, funs)
}
\arguments{
\item{nl}{nl object}

\item{metrics}{vector of strings defining metric columns for evaluation. Defaults to metrics of the experiment within the nl object}

\item{funs}{list with the summary metrics for the sensitivity results}
}
\description{
Analyze NetLogo simulation output of simdesign sobol
}
\details{
The function calculates sobol sensitivity indices from the output results using the \link[sensitivity]{sensitivity} package.
}
\keyword{internal}
