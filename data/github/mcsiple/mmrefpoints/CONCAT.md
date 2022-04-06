
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mmrefpoints: Projecting long-term marine mammal abundance with bycatch

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://github.com/mcsiple/mmrefpoints/issues)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6336011.svg)](https://doi.org/10.5281/zenodo.6336011)

<!-- badges: end -->

`mmrefpoints` is an R package that generates marine mammal population
projections based on starting abundance, life history, and bycatch
rates, based on the BALEEN II population dynamics model.

## Authors

Margaret C. Siple  
André E. Punt  
Tessa B. Francis  
Philip S. Hammond  
Dennis Heinemann  
Kristy J. Long  
Jeffrey E. Moore  
Maritza Sepúlveda  
Randall R. Reeves  
Guðjón Már Sigurðsson  
Gísli Víkingsson  
Paul R. Wade  
Rob Williams  
Alexandre N. Zerbini

## Contents

-   [Need](#need)
-   [Details](#details)
-   [Installation](#installation)
-   [Contributing](#contributing)
-   [References](#references) <!-- end toc -->

## Need

Stakeholders involved in the management of marine mammal bycatch in
marine fisheries need tools to simulate the effects of management
decisions on marine mammal populations. Population models are a key part
of this process. This package contains the tools to simulate marine
mammal populations and an app that shows model outputs in a
user-friendly way.

## Details

This R package contains the functions used in the Marine Mammal Bycatch
Impacts Exploration Tool (MMBIET), a Shiny app built by Margaret Siple,
André Punt, and the Ocean Modeling Forum’s [Marine Mammal Bycatch
Working
Group](https://oceanmodelingforum.org/working-groups/marine-mammal-bycatch-working-group/).

The functions in this package, and the Shiny app, are intended to be
used in cases where data on bycatch and/or population status are sparse
or unavailable.

Our target audience is stakeholders interested in projecting marine
mammal populations to examine the impacts of bycatch. Those code could
also be used as a teaching tool, or for anyone who is more familiar with
R than FORTRAN and wants to use some components of the BALEEN II model
(Punt 1999).

## Installation

This package can be downloaded directly from GitHub:

    devtools::install_github("mcsiple/mmrefpoints")

-   NOTE: For Linux users, if you run into an error about the `magick`
    package, try installing `magick` first using the instructions
    [here](https://cran.r-project.org/web/packages/magick/vignettes/intro.html).

## Contributing [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)

### Community guidelines

We would like this package to be sustainable in the long term and
welcome contributions. If you are interested in contributing, please
check out our [Contribution
Guide](https://github.com/mcsiple/mmrefpoints/blob/master/CONTRIBUTING.md).

Bugs and enhancements are tracked through GitHub issues. If you have a
bug to report, there is a [bug report
template](https://github.com/mcsiple/mmrefpoints/blob/master/.github/ISSUE_TEMPLATE/bug_report.md)
to help maximize the benefit of your report for everyone. The same goes
for requests for
[enhancements](https://github.com/mcsiple/mmrefpoints/blob/master/.github/ISSUE_TEMPLATE/feature_request.md).

## Accessing the MMBIET Shiny app

The functions in this package can also be accessed through the Shiny app
for this project, which is located online
[here](https://msiple.shinyapps.io/mmrefpoints/). The app provides an
easy way to explore outcomes and print out a report with inputs and
outputs.

The mmBIET Shiny app can also be accessed through the R package:

``` r
library(mmrefpoints)
run_app()
```

<img src="https://github.com/mcsiple/ltbycatch/blob/master/docs/screenshot1.png" alt="screenshot1" width="400">

## Functionality

The foundation of this package is an age-structured population
projection model with bycatch mortality. Key functions in this package:

| Function      | Purpose                                                         |
|:--------------|:----------------------------------------------------------------|
| dynamics()    | Generate a single trajectory for marine mammal population size  |
| projections() | Generate several trajectories for marine mammal population size |

To create a single projection for a marine mammal population, use the
`dynamics()` function:

``` r
x <- mmrefpoints::dynamics(lh.params = list(S0 = 0.944, S1plus = 0.99, 
                           K1plus = 9000, AgeMat = 17, z = 2.39, nages = 25, lambdaMax = 1.04),
                           InitDepl = 0.6, 
                           ConstantCatch = NA, 
                           ConstantF = rep(0.01, times = 100), 
                           nyears = 100)
plot(1:100, x$TotalPop, type = 'l', xlab = "Year", ylab = "Population size")
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

Variation in the model is introduced through variation in bycatch
mortality over time and uncertainty in the estimate of starting
abundance.

``` r
x <- mmrefpoints::projections(
  NOut = 1,
  ConstantBycatch = list(Catch = 50, CV = 0.7),
  InitDepl = 0.6,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, nages = 25, z = 2.39, lambdaMax = 1.04
  ),
  nyears = 100, obs_CV = 0.1
)

# One trajectory with bycatch uncertainty and an observation CV:
plot(x$trajectories, type = "l", xlab = "Year", ylab = "Population size")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Projections shown in the app are based on simulation parameters provided
by the user. They include a “high”, “medium”, and “low” bycatch level
based on a user-determined range.

``` r
x_lo <- mmrefpoints::projections(
  NOut = 100,
  ConstantBycatch = list(Catch = 0, CV = 0),
  InitDepl = 0.6,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, nages = 25, z = 2.39, lambdaMax = 1.04
  ),
  nyears = 100, obs_CV = 0.1
)

x_med <- mmrefpoints::projections(
  NOut = 100,
  ConstantBycatch = list(Catch = 50, CV = 0.7),
  InitDepl = 0.6,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, nages = 25, z = 2.39, lambdaMax = 1.04
  ),
  nyears = 100, obs_CV = 0.1
)

x_hi <- mmrefpoints::projections(
  NOut = 100,
  ConstantBycatch = list(Catch = 200, CV = 0.7),
  InitDepl = 0.6,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, nages = 25, z = 2.39, lambdaMax = 1.04
  ),
  nyears = 100, obs_CV = 0.1
)

mmrefpoints::plot_proj(high = x_hi,
                       med = x_med,
                       low = x_lo,
                       years.plot = 100,
                       ylims = c(0, 9000),
                       K1plus = 9000)
#> Warning: Removed 38 row(s) containing missing values (geom_path).
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

## References

Punt, A. E. 1999. Annex R: A full description of the standard Baleen II
model and some variants thereof. Division of Marine Research, CSIRO
Marine Laboratories, Hobart, Australia. Available from
<https://duwamish.lib.washington.edu/uwnetid/illiad.dll?Action=10&Form=75&Value=1651729>
(accessed August 7, 2018).

## How to cite

To cite this package or the MMBIET Shiny app, please use the following
citation:

> Margaret C. Siple, André E. Punt, Tessa B. Francis, Philip S. Hammond,
> Dennis Heinemann, Kristy J. Long, Jeffrey E. Moore, Randall R. Reeves,
> Sepúlveda Maritza, Guðjón Már Sigurðsson, Gísli Víkingsson, Paul R.
> Wade, Rob Williams and Alexandre N. Zerbini (NA). mmrefpoints: Project
> Marine Mammal Populations and Calculate Reference Points. R package
> version 1.0.1. url: <https://github.com/mcsiple/mmrefpoints> doi:
> 10.5281/zenodo.4758401

NOTE that if you want to cite all versions of the software, you can use
the doi
[10.5281/zenodo.4758401](https://zenodo.org/record/5949332#.Yf1infXMI-Q).
When additional releases happen, there will be a doi for each new
release as well.
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at [INSERT CONTACT
METHOD]. All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/
code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
---
title: 'mmrefpoints: Projecting long-term marine mammal abundance with bycatch'
tags:
  - Shiny
  - R
  - marine mammal bycatch
authors:
  - name: Margaret C. Siple^[*]
    orcid: 0000-0002-4260-9856
    affiliation: 1
  - name: André E. Punt
    orcid: 0000-0001-8489-2488
    affiliation: 2
  - name: Tessa B. Francis
    orcid: 0000-0002-3383-5392
    affiliation: 3
  - name: Philip S. Hammond
    orcid: 0000-0002-2381-8302
    affiliation: 4
  - name: Dennis Heinemann
    orcid: 0000-0002-1434-2445
    affiliation: 5
  - name: Kristy J. Long
    orcid: 0000-0001-6970-0935
    affiliation: 6
  - name: Jeffrey Moore
    orcid: 0000-0003-3715-7442
    affiliation: 7
  - name: Maritza Sepúlveda
    orcid: 0000-0002-1403-176X
    affiliation: 8
  - name: Randall R. Reeves
    orcid: 0000-0002-6512-6507
    affiliation: 9
  - name: Guðjón Már Sigurðsson
    orcid: 0000-0001-9390-6693
    affiliation: 10
  - name: Gísli Víkingsson
    orcid: 0000-0002-4501-193X
    affiliation: 10
  - name: Paul R. Wade
    orcid: 0000-0003-2428-9323
    affiliation: 11
  - name: Rob Williams
    orcid: 0000-0001-7496-453X
    affiliation: 12
  - name: Alexandre N. Zerbini
    orcid: 0000-0002-9776-6605
    affiliation: 13
affiliations:
 - name: Resource Assessment and Conservation Engineering Division, Alaska Fisheries Science Center, National Oceanic and Atmospheric Administration, Seattle, WA, 98115, USA
   index: 1
 - name: School of Aquatic and Fishery Sciences, University of Washington, 1122 NE Boat St, Seattle, WA 98115
   index: 2
 - name: Puget Sound Institute, University of Washington Tacoma, 326 East D Street, Tacoma, WA 98421, USA
   index: 3
 - name: Sea Mammal Research Unit, Scottish Oceans Institute, University of St Andrews, Fife KY16 8LB, UK
   index: 4
 - name: U.S. Marine Mammal Commission, 4340 East-West Hwy, Rm 700, Bethesda, MD 20814, USA
   index: 5
 - name: Office of Protected Resources, NOAA's National Marine Fisheries Service, Silver Spring, MD 20910, USA
   index: 6
 - name: Protected Resources Division, NOAA SWFSC, La Jolla, CA 92037, USA
   index: 7
 - name: Facultad de Ciencias, Universidad de Valparaíso, Gran Bretaña 1111, Playa Ancha, Valparaíso, Chile
   index: 8
 - name: Okapi Wildlife Associates, Hudson, Quebec, Canada
   index: 9
 - name: Marine and Freshwater Research Institute, Skúlagata 4, 121, Reykjavík, Iceland
   index: 10
 - name: Marine Mammal Laboratory, NOAA AFSC, Seattle, WA, 98115-6349, USA
   index: 11
 - name: Oceans Initiative, 117 E. Louisa Street No. 135. Seattle, WA 98102, USA
   index: 12
 - name: Cascadia Research Collective, 218 ½ 4th Ave W, Olympia, WA, 98501, USA
   index: 13
 - name: Marine Ecology and Telemetry Research, 2468 Camp McKenzie Tr NW, Seabeck, WA, 98380, USA
   index: 14
   
date: 8 March 2022
bibliography: paper.bib
---

## Statement of Need
Fisheries bycatch is one of the top threats to marine mammals worldwide. Data on bycatch and marine mammal abundance are required to make management decisions, yet many marine mammal populations exist in locations where these data are sparse. Population models offer a way to explore the long-term impacts of different bycatch scenarios even when data are sparse or absent. Our modeling tool, `mmrefpoints`, uses very basic information about marine mammal life history, estimates of abundance and bycatch, and population models to project long-term (e.g., 20-100 yr) outcomes of different bycatch rates. These long-term outcomes are the basis of several reference points used for marine mammal management including the Potential Biological Removal (PBR) approach [@wade_calculating_1998]. The goal is to make complex population models accessible to managers and stakeholders who need them, and to students who are learning about risk-based approaches to managing marine mammal populations.

## Summary
This tool provides a way for managers and other stakeholders to explore bycatch scenarios, based on simple information about marine mammal life history and rough estimates of abundance and bycatch. The tool consists of an R package [@r2021] and a Shiny application [@shiny2021]. The primary machinery in the package is an age-structured population dynamics model, which is used to model future population size based on the current population size, life history traits, and bycatch rates. The package also contains tools for calculating performance metrics, as well as the U.S. reference point for bycatch (Potential Biological Removal or PBR) and a solver that estimates the maximum bycatch rate that will meet management objectives. For users who would prefer to see outputs in an interactive user interface, there is a Shiny app called the Marine Mammal Bycatch Impacts Exploration Tool (MMBIET) that shows projections, explains and calculates reference points, and creates a report summarizing inputs and outputs. The app is hosted publicly on the web via the [shinyapps.io server](https://msiple.shinyapps.io/mmrefpoints/), or it can be run locally by installing the package and running the `run_app()` function in the R or RStudio console:

```{r}
remotes::install_github("mcsiple/mmrefpoints")
library(mmrefpoints)
run_app()
```

## Population model
The population model is a simplified version of the model described by @breiwick_population_1984 and @punt_a._e._annex_1999. The model is a single-sex, age-structured population model. The number of calves or pups born each year is density-dependent, based on the number of mature adults, the pregnancy rate at pre-exploitation equilibrium, maximum theoretical fecundity rate, degree of compensation, and the total abundance relative to carrying capacity K. 
User inputs determine the parameters for calf survival, adult survival, age at maturity, and plus group age. Default parameters are based on those used by @punt_conserving_2018, with additional default parameters for other species drawn from literature values [@arso_civil_variations_2019; @olafsdottir_growth_2003; @speakman_mark-recapture_2010]. Because of the variation in pinniped life history parameters, parameter estimates for several pinniped species are provided within the Shiny app. These life history parameters are taken from @butterworth_effects_1995, @delong_age-_2017, @dillingham_improved_2016, @hastings_sex-_2012, and @moore_unpublished_2019.

The full model description including equations is contained in the “Model description” vignette and the “About the Model” tab of the app.

## Intended use
### `mmrefpoints` and the MMBIET app are intended to be used for the following:

* Exploring outcomes for bycatch rates in populations with little or no information
* Calculating reference points like PBR to estimate maximum allowable bycatch for marine mammal populations

### The MMBIET app is *not* intended for the following:

* Permitting specific management actions regarding bycatch. We refer users to @hammond_estimating_2021, @moore_estimating_2021, and @wade_best_2021 for guidance on developing a management program for marine mammal bycatch including monitoring.
* Calculating PBR for marine mammal stocks that already have a stock assessment. If reference points have already been calculated for the stock, those should be used.
* Fitting population models to data (we direct readers to other tools like [rSPAMM](https://github.com/NorskRegnesentral/rSPAMM) for this type of need)

## Acknowledgements
The authors would like to thank several pilot testers for reviewing a beta version of the MMBIET Shiny app, and Christine Stawitz and Jay Barlow for reviewing an earlier version of the `mmrefpoints` R package and the app.

## Ongoing projects using MMBIET
At the time of this submission, three papers have cited `mmrefpoints` and/or MMBIET:

Hammond, P. S., Francis, T. B., Heinemann, D., Long, K. J., Moore, J. E., Punt, A. E., Reeves, R. R., Sepúlveda, M., Sigurðsson, G. M., Siple, M. C., Víkingsson, G., Wade, P. R., Williams, R., & Zerbini, A. N. (2021). Estimating the Abundance of Marine Mammal Populations. Frontiers in Marine Science, 8, 1316. [https://doi.org/10.3389/fmars.2021.735770](https://doi.org/10.3389/fmars.2021.735770)

Moore, J. E., Heinemann, D., Francis, T. B., Hammond, P. S., Long, K. J., Punt, A. E., Reeves, R. R., Sepúlveda, M., Sigurðsson, G. M., Siple, M. C., Víkingsson, G. A., Wade, P. R., Williams, R., & Zerbini, A. N. (2021). Estimating Bycatch Mortality for Marine Mammals: Concepts and Best Practices. Frontiers in Marine Science, 8, 1793. [https://doi.org/10.3389/fmars.2021.752356](https://doi.org/10.3389/fmars.2021.752356)

Wade, P. R., Long, K. J., Francis, T. B., Punt, A. E., Hammond, P. S., Heinemann, D., Moore, J. E., Reeves, R. R., Sepúlveda, M., Sullaway, G., Sigurðsson, G. M., Siple, M. C., Víkingsson, G. A., Williams, R., & Zerbini, A. N. (2021). Best Practices for Assessing and Managing Bycatch of Marine Mammals. Frontiers in Marine Science, 8, 1566. [https://doi.org/10.3389/fmars.2021.757330](https://doi.org/10.3389/fmars.2021.757330)


## Financial support
Support for this project is provided by the Ocean Modeling Forum and the Lenfest Ocean Program. M. Siple was supported by the Ocean Modeling Forum and a James S. McDonnell Postdoctoral Fellowship.

## References

# mmrefpoints 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
# Contributing

## Welcome
Thank you for your interest in contributing to `mmrefpoints`! We hope to maintain this package as a resource for a broad set of stakeholders and modelers.

## Table of Contents

[Short Links](#short-links)

[How Can I Contribute?](#how-can-i-contribute)

[Style Guides](#style-guides)


## Short Links
* View the rest of the Ocean Modeling Forum Marine Mammal Bycatch Working Group products [here](https://oceanmodelingforum.org/working-groups/marine-mammal-bycatch-working-group/)
* Report bugs [here](https://github.com/mcsiple/mmrefpoints/issues)
* Contact Margaret at margaret (dot) siple (at) noaa.gov

## How Can I Contribute?
 * [Reporting Bugs](#reporting-bugs)
 * [Suggesting Enhancements](#suggesting-enhancements)
 * [Pull Requests](#pull-requests)
 
## Style Guides
* [Write an issue on GitHub](#issue-style-guides)
* [Code style guides](#code-style-guides)

## Reporting Bugs 
Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/). If you find a bug, please make your report as detailed as possible. Fill out the [required bug report template](https://github.com/mcsiple/mmrefpoints/blob/master/.github/ISSUE_TEMPLATE/bug_report.md); the information it asks for helps us resolve issues faster.

**Note:** If you find a **Closed** issue that looks like it is the same thing that you're experiencing, open a new issue and include a link to the original issue in the body of your new one (you can tag issues by a hashtag (#) followed by the issue number).

## Suggesting Enhancements 
:pencil: We are always interested in how to make this package more useful and accessible. If you have an enhancement to suggest, please submit a [Feature Request](https://github.com/mcsiple/mmrefpoints/blob/master/.github/ISSUE_TEMPLATE/feature_request.md). 

### How do I submit a good enhancement suggestion?
Enhancement suggestions are tracked as [GitHub issues](https://guides.github.com/features/issues/). Create an issue on the mmrefpoints repository and include the following:

* **Use a clear and descriptive title** for the issue to identify the suggestion.
* **Provide a step-by-step description of the suggested enhancement** in as many details as possible.
* **Provide specific examples to demonstrate the steps**. Include copy/pasteable snippets which you use in those examples, as [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the current behavior** and **explain which behavior you expected to see instead** and why.
* **Include screenshots and animated GIFs** which help you demonstrate the steps or point out the part of Atom which the suggestion is related to. You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on macOS and Windows, and [this tool](https://github.com/colinkeenan/silentcast) or [this tool](https://github.com/GNOME/byzanz) on Linux.
* **Explain why this enhancement would be useful** to most mmrefpoints users and stakeholders.
* **Specify what version of mmrefpoints you are using** and whether you are using the Shiny app in the browser, the local Shiny app within the mmrefpoints package, or the functions in mmrefpoints by themselves.
* **Specify the name and version of the OS you're using.**


## Pull requests
If you have a code suggestion in hand, please submit a pull request. 

## Write an issue on GitHub
### Git commit messages
(many of these are lifted from the Atom [style guide](https://github.com/atom/atom/blob/master/CONTRIBUTING.md) for commit messages)

* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* When only changing documentation, include `[ci skip]` in the commit title
* Consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :racehorse: `:racehorse:` when improving performance
    * :memo: `:memo:` when writing docs
    * :penguin: `:penguin:` when fixing something on Linux
    * :apple: `:apple:` when fixing something on macOS
    * :checkered_flag: `:checkered_flag:` when fixing something on Windows
    * :bug: `:bug:` when fixing a bug
    * :fire: `:fire:` when removing code or files
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :arrow_up: `:arrow_up:` when upgrading dependencies
    * :arrow_down: `:arrow_down:` when downgrading dependencies

## Code style guides
We have adapted the following coding style guidelines from simulation developers with established guidelines, like [ss3sim](https://github.com/ss3sim). 

You will notice that some of the functionality of mmrefpoints may not follow all the advice in this style guide. I am working on streamlining much of this code (especially the Shiny code :construction_worker:), so in the meantime we request that you do as we style guide, not as we do. 

### Coding style
* Guides
  * Follow the [tidyverse style guide](https://style.tidyverse.org)
  * References on R package development from Hadley's book [R Packages](http://r-pkgs.had.co.nz)
  * [Advanced R](http://adv-r.had.co.nz)
  * [CRAN page](http://cran.r-project.org/doc/manuals/r-release/R-exts.html)
  * [Semantic versioning](https://semver.org/)
* Code rules
  * Use `vapply` or `lapply` instead of `sapply` because `sapply` is not consistent in the what class of object it returns (this may not yet be implemented in the package)
  * Use `[` rather than `subset` because `subset` will lead to unassigned objects, i.e., R will not know where the column vector name in the `subset` call comes from
  * Use `<-` rather than `=` to assign objects
  * Use `seq_len(NROW(x))` or `seq_along(x)` rather than `1:NROW(x)` or `1:length(x)` (not yet implemented in the package)
  * [Importing functions from other packages in ss3sim functions](Importing-functions-from-other-packages)

### Importing functions from other packages
(these guidelines are lifted and adapted from the ss3sim [developer wiki](https://github.com/ss3sim/ss3sim/wiki/developers), which contains a lot of useful resources for developing code)

What if you need to use a function from another package? Normally you might attach all the functions from a package with a call to `library()`. This is a bad practice in the context of writing an R package. We shouldn't be attaching other packages in the user's environment. In the context of writing a package we will instead import that function for our use. The package structure ensures that the user has the package installed.

For more details, see http://r-pkgs.had.co.nz/namespace.html

Say you want to use the function `mle2` from the bbmle package.

Here are the steps to use the function within your package:

* Add the package to the list of `Imports` in the `DESCRIPTION` file. If it's possible that functionality has changed in the past for the function then specify a `>=` version number for the package. You can find the version you have installed in many ways. One way is to type `help(package = "bbmle")` within R. You can also look at your `sessionInfo()` in R.
* In the function Rogxyen documentation include `@importFrom bbmle mle2`. If there are multiple functions to import than list them as `@importFrom packagename function1 function2` etc.
* Call the function as you would normally, e.g., `mle2(x, ...)`

The above steps outline one of two preferred notations. Within mmrefpoints we also use the `::` notation to be explicit about what package a function comes from. E.g., `shiny.i18n::update_lang()` clearly tells readers of the code that the `update_lang()` function comes from the shiny.i18n package without readers having to navigate to the top of the code to read the roxygen documentation. This notation also requires that packages be listed in the `Imports` section of the `DESCRIPTION` file. We support either notation, but we have largely switched to this one as time goes on.

Occasionally your code might call so many functions from another package that you just want all the functions available. In that case use `@import bbmle` at the top. I have done this in places in mmrefpoints. But, it is best practices to import specific functions that we need because it saves time, which users appreciate. It also makes the code more explicit.


## Code of conduct
Contributors to mmrefpoints follow the [Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/) code of conduct. For answers to common questions about this code of conduct, see https://www.contributor-covenant.org/faq.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is. Please check whether you encountered this bug:

- [ ] while using the Shiny app in the browser (https://msiple.shinyapps.io/mmrefpoints/)
- [ ] while using the Shiny app locally (i.e., using `run_app()` in R on your computer)
- [ ] while using the `mmrefpoints` functions locally

**To Reproduce**
Steps to reproduce the behavior:

1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**

 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Information about your version of R**

 - R version (e.g., 4.1.2, get this by running `R.Version()` in the console)
 - Session info (`sessionInfo()` in the console)
 - How you are accessing R (RStudio and version, R GUI)

**Additional context**
Add any other context about the problem here.
---
output: github_document
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

# mmrefpoints: Projecting long-term marine mammal abundance with bycatch

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://github.com/mcsiple/mmrefpoints/issues)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6336011.svg)](https://doi.org/10.5281/zenodo.6336011)

<!-- badges: end -->

`mmrefpoints` is an R package that generates marine mammal population projections based on starting abundance, life history, and bycatch rates, based on the BALEEN II population dynamics model.

## Authors
Margaret C. Siple  
André E. Punt  
Tessa B. Francis  
Philip S. Hammond  
Dennis Heinemann  
Kristy J. Long  
Jeffrey E. Moore  
Maritza Sepúlveda  
Randall R. Reeves  
Guðjón Már Sigurðsson  
Gísli Víkingsson  
Paul R. Wade  
Rob Williams  
Alexandre N. Zerbini  


## Contents
-   [Need](#need)
-   [Details](#details)
-   [Installation](#installation)
-   [Contributing](#contributing)
-   [References](#references)
<!-- end toc -->

## Need
Stakeholders involved in the management of marine mammal bycatch in marine fisheries need tools to simulate the effects of management decisions on marine mammal populations. Population models are a key part of this process. This package contains the tools to simulate marine mammal populations and an app that shows model outputs in a user-friendly way.

## Details
This R package contains the functions used in the Marine Mammal Bycatch Impacts Exploration Tool (MMBIET), a Shiny app built by Margaret Siple, André Punt, and the Ocean Modeling Forum's [Marine Mammal Bycatch Working Group](https://oceanmodelingforum.org/working-groups/marine-mammal-bycatch-working-group/). 

The functions in this package, and the Shiny app, are intended to be used in cases where data on bycatch and/or population status are sparse or unavailable. 

Our target audience is stakeholders interested in projecting marine mammal populations to examine the impacts of bycatch. Those code could also be used as a teaching tool, or for anyone who is more familiar with R than FORTRAN and wants to use some components of the BALEEN II model (Punt 1999).

## Installation
This package can be downloaded directly from GitHub:

```
devtools::install_github("mcsiple/mmrefpoints")
```

* NOTE: For Linux users, if you run into an error about the `magick` package, try installing `magick` first using the instructions [here](https://cran.r-project.org/web/packages/magick/vignettes/intro.html).

## Contributing [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)
### Community guidelines
We would like this package to be sustainable in the long term and welcome contributions. If you are interested in contributing, please check out our [Contribution Guide](https://github.com/mcsiple/mmrefpoints/blob/master/CONTRIBUTING.md). 

Bugs and enhancements are tracked through GitHub issues. If you have a bug to report, there is a [bug report template](https://github.com/mcsiple/mmrefpoints/blob/master/.github/ISSUE_TEMPLATE/bug_report.md) to help maximize the benefit of your report for everyone. The same goes for requests for [enhancements](https://github.com/mcsiple/mmrefpoints/blob/master/.github/ISSUE_TEMPLATE/feature_request.md). 





## Accessing the MMBIET Shiny app
The functions in this package can also be accessed through the Shiny app for this project, which is located online [here](https://msiple.shinyapps.io/mmrefpoints/). The app provides an easy way to explore outcomes and print out a report with inputs and outputs. 

The mmBIET Shiny app can also be accessed through the R package:
```{r eval=FALSE}
library(mmrefpoints)
run_app()
```

<img src="https://github.com/mcsiple/ltbycatch/blob/master/docs/screenshot1.png" alt="screenshot1" width="400">

## Functionality

The foundation of this package is an age-structured population projection model with bycatch mortality. Key functions in this package:

```{r eval=TRUE,echo=FALSE, results='asis'}
x <- data.frame("Function" = c("dynamics()",
                               "projections()"),
                "Purpose" = c("Generate a single trajectory for marine mammal population size", 
                              "Generate several trajectories for marine mammal population size"))

knitr::kable(x)
```



To create a single projection for a marine mammal population, use the `dynamics()` function:
```{r message=FALSE}
x <- mmrefpoints::dynamics(lh.params = list(S0 = 0.944, S1plus = 0.99, 
                           K1plus = 9000, AgeMat = 17, z = 2.39, nages = 25, lambdaMax = 1.04),
                           InitDepl = 0.6, 
                           ConstantCatch = NA, 
                           ConstantF = rep(0.01, times = 100), 
                           nyears = 100)
plot(1:100, x$TotalPop, type = 'l', xlab = "Year", ylab = "Population size")
```

Variation in the model is introduced through variation in bycatch mortality over time and uncertainty in the estimate of starting abundance.

```{r message=FALSE}
x <- mmrefpoints::projections(
  NOut = 1,
  ConstantBycatch = list(Catch = 50, CV = 0.7),
  InitDepl = 0.6,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, nages = 25, z = 2.39, lambdaMax = 1.04
  ),
  nyears = 100, obs_CV = 0.1
)

# One trajectory with bycatch uncertainty and an observation CV:
plot(x$trajectories, type = "l", xlab = "Year", ylab = "Population size")
```

Projections shown in the app are based on simulation parameters provided by the user. They include a "high", "medium", and "low" bycatch level based on a user-determined range. 

```{r message=FALSE}
x_lo <- mmrefpoints::projections(
  NOut = 100,
  ConstantBycatch = list(Catch = 0, CV = 0),
  InitDepl = 0.6,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, nages = 25, z = 2.39, lambdaMax = 1.04
  ),
  nyears = 100, obs_CV = 0.1
)

x_med <- mmrefpoints::projections(
  NOut = 100,
  ConstantBycatch = list(Catch = 50, CV = 0.7),
  InitDepl = 0.6,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, nages = 25, z = 2.39, lambdaMax = 1.04
  ),
  nyears = 100, obs_CV = 0.1
)

x_hi <- mmrefpoints::projections(
  NOut = 100,
  ConstantBycatch = list(Catch = 200, CV = 0.7),
  InitDepl = 0.6,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, nages = 25, z = 2.39, lambdaMax = 1.04
  ),
  nyears = 100, obs_CV = 0.1
)

mmrefpoints::plot_proj(high = x_hi,
                       med = x_med,
                       low = x_lo,
                       years.plot = 100,
                       ylims = c(0, 9000),
                       K1plus = 9000)
```

## References
Punt, A. E. 1999. Annex R: A full description of the standard Baleen II model and some variants thereof. Division of Marine Research, CSIRO Marine Laboratories, Hobart, Australia. Available from https://duwamish.lib.washington.edu/uwnetid/illiad.dll?Action=10&Form=75&Value=1651729 (accessed August 7, 2018).

## How to cite
To cite this package or the MMBIET Shiny app, please use the following citation:

> Margaret C. Siple, André E. Punt, Tessa B. Francis, Philip S. Hammond, Dennis Heinemann, Kristy J. Long, Jeffrey E. Moore,
  Randall R. Reeves, Sepúlveda Maritza, Guðjón Már Sigurðsson, Gísli Víkingsson, Paul R. Wade, Rob Williams and Alexandre N.
  Zerbini (NA). mmrefpoints: Project Marine Mammal Populations and Calculate Reference Points. R package version 1.0.1.
  url: <https://github.com/mcsiple/mmrefpoints> doi: 10.5281/zenodo.4758401
  
NOTE that if you want to cite all versions of the software, you can use the doi [10.5281/zenodo.4758401](https://zenodo.org/record/5949332#.Yf1infXMI-Q). When additional releases happen, there will be a doi for each new release as well. ---
title: "About the projection model"
author: "Siple MC and Punt AE"
date: 'Last updated: June 12, 2020'
output:
  html_document:
    highlight: tango
    theme: flatly
  pdf_document: default
  word_document: default
---

<!-- NOTE: in order to include this html in the shiny app, delete everything in the knitted html before <body> and after </body> -->

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# devtools::install_github("haozhu233/kableExtra")
 library(kableExtra)
```


### Population model

Population projections are generated using a single-sex age-structured model. In this projection model, the number of calves or pups born each year is density dependent, with the extent of density dependence being a function of the number of mature adults $\tilde N$, the fecundity (pregnancy rate) at pre-exploitation equilibrium $f_0$, the maximum theoretical fecundity rate $f_{max}$, the degree of compensation $z$, and the abundance of individuals aged 1+ ($N_{y+1}^{1+}$) relative to carrying capacity $K^{1+}$.

The number of individuals age 1 and older is a function of calf/pup survival, $S_0$, and the survival rate of adults, $S_{1+}$, and removals due to bycatch mortality, $C_y$ (all symbols are defined in Table 1).


$$
\tag{Eq 1}
N_{y+1,a} = 
\begin{cases}
\tilde N_{y+1}\bigg(f_0 + (f_{max}-f_0)\bigg[1-\bigg(\frac{N_{y+1}^{1+}}{K^{1+}}\bigg)^z\bigg]\bigg) &   a=0\\ 
N_{y,0}S_0 &    a=1 \\
(N_{y,a-1}-C_{y,a-1})S_{1+}  & 2\leq a<x \\
(N_{y,x-1}-C_{y,x-1})S_{1+} + (N_{y,x}-C_{y,x})S_{1+} & a=x \\
\end{cases}
$$

where $N_{y,a}$ is the number of animals of age $a$ at the start of year $y$. The number of animals of age $a$ dying due to bycatch mortality during year $y$, $C_{y,a}$, is removed uniformly from the 1+ component of the population (calves and pups are assumed not to die due to bycatch), i.e.:

$$
\tag{Eq 2}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
\frac{C_y N_{y,a}}{N_y^{1+}} & a>0 \\
\end{cases}
$$

Given this assumption, our model will not adequately characterize cases where bycatch mortality occurs predominantly among calves/pups. 

If the user specifies a constant rate of bycatch mortality, individuals are removed from the population according to a bycatch mortality rate $E_y$ and vulnerability (which is assumed constant through time, and uniform on age 1+ individuals):

$$
\tag{Eq 3}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
N_{y,a} E_{y} & a>0 \\
\end{cases}
$$
If bycatch is specified by the user as a number of individuals per year, annual bycatch rates $E_y$ are lognormally distributed with mean $log(C_{user})$ and standard deviation $\sigma_E$ of $\sqrt{log(CV_C^2 + 1)}$, where $C_{user}$ is the bycatch in numbers per year defined by the user and $CV_C$ is the bycatch CV defined by the user. If bycatch is specified as a constant rate, $E_y$ is beta-distributed with mean $E_{user}$ and standard deviation $E_{user} \cdot CV_E$.




**Table 1.** Symbols included in the projection model.
 
```{r echo = F}
Symbol=c(#User-specified parameters
         "$S_0$",
         "$S_{1+}$",
         "$x$",
         "$\\lambda_{max}$",
          "$a_p$",
         "$E_y$",
         # Parameters derived from user-specified parameters
         "$f_0$",
         "$f_{max}$",
         "$K^{1+}$",
         "$z$",
         # Derived variables
         "$\\tilde N_{y,a}$", # need double dash for using in kableExtra
         "$N_{y,a}$",
         "$N_y^{1+}$",
         "$C_{y,a}$")

Description=c(#User-specified parameters
              "Pup or calf survival",
              "Survival of individuals aged 1+",
              "Plus-group age",
              "Maximum steady rate of increase (population growth rate)",
              "Age at first parturition",
              "Bycatch mortality rate in year $y$ (specified by the user or computed from the user-specified total bycatch mortality)", 
              # Parameters derived from user-specified parameters
              "Fecundity (pregnancy rate) at pre-exploitation equilibrium",
              "Maximum theoretical fecundity (pregnancy rate)",
              "Carrying capacity in terms of 1+ component of the population",
              "Degree of compensation",
              # Derived variables
              "Number of mature animals of age $a$ at the start of year $y$",
              "Number of animals of age $a$ at the start of year $y$",
              "Number of animals aged 1 and older at the start of year $y$",
              "Mortality due to bycatch of animals of age $a$ in year $y$")
Code <- c(#User-specified parameters
         "S0",
         "S1plus",
         "nages",
         "lambdaMax",
         "AgeMat + 1",
         "ConstantF",
         # Parameters derived from user-specified parameters
         "Fec0",
         "fmax",
         "K1plus",
         "z",
         # Derived variables
         "--", # need double dash for using in kableExtra
         "--",
         "--",
         "--")

Symbol_definitions <- data.frame(cbind(Symbol,
                                      Description, Code))
Symbol_definitions$Code <- cell_spec(Symbol_definitions$Code, font_size = 12, extra_css = "font-family:'Monaco', sans-serif;")
  
knitr::kable(format="html", #html  - change to pandoc for ms word 
             Symbol_definitions, 
             escape = FALSE) %>% 
            kable_styling('bordered') %>%
            pack_rows('User-specified parameters', 1, 6, label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
            pack_rows('Parameters derived from user-specified parameters', 7, 10, 
                      label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
            pack_rows('Derived variables', 11, 14, 
                      label_row_css = "font-style: italic;border-bottom: 1px solid black;")
```

### Performance measures
There are three primary performance measures included in the outputs:

-	the probability that the abundance of age 1+ animals exceeds the Maximum Net Productivity Level (MNPL) after 50 and 100 years;
-	the abundance of age 1+ animals compared to 1+ carrying capacity after 10, 20, and 50 years; and
-	the abundance of age 1+ animals compared to the abundance of age 1+ animals under a no-bycatch scenario after 10, 20, and 50 years.

These performance measures are all related to population recovery. Wade (1998) identifies a performance ‘standard’ of a 95% probability of recovery to MNPL after 100 years for a population initially at 30% of its carrying capacity, with an MNPL of 0.5K. MNPL is the lower bound for the U.S. Marine Mammal Protection Act Optimum Sustainable Population (OSP), which is defined as the “number of any animals which will result in the maximum productivity of the population or the species, keeping in mind the carrying capacity of the habitat and the health of the ecosystem of which they form a constituent element” (16 USCS § 1362 (9)). In the US management scheme, marine mammal stocks are considered depleted when they are below OSP (16 USCS § 1362 (1B)).

Single values reported in tables are the medians from the number of projections.


### Solving for bycatch rate

The app can solve for the bycatch level that would result in recovery to a given proportion of K in a specific amount of time. The rebuilding probability is calculated as the proportion of simulations in which the population has recovered to the rebuilding goal by the rebuilding year. The “Solve for bycatch” tab uses root finding (the function `uniroot()` in our code) to find the bycatch mortality rate that minimizes the difference between the desired recovery probability and the recovery probability under bycatch mortality rate $E$.


### Maximum theoretical fecundity
The maximum theoretical fecundity rate is calculated based on the population size when the population is growing at its maximum population growth rate $\lambda_{max}$. Parts of this derivation can be found in Butterworth and Punt (1992) and Breiwick et al. (1984). Starting with the population model for the female component of the population:

$$
\tag{Eq 4}
N_{t+1} = N_t(1-E)S_{1+}+N_{t-a_p}f q_f S_0 (S_{1+})^{a_p -2}(1-E)^{a_p-1-a_r}
$$

where $a_p$ is the age at first parturition (assumed to be one year after the age at maturity), $f$ is the fecundity (sometimes expressed as the product of pregnancy rate and pup survival), $q_f$ is the proportion of calves/pups that are female and $a_r$ is the age at which the mammals first suffer any bycatch mortality (in our case, this is equal to 1 year). When the population is growing at its fastest rate, fishing is zero, and $N_{t+1}=N_t \lambda_{max}$. In our case $q_f = 1$ because we are modeling all adults, instead of just mature females. Thus, the above equation becomes:

$$
\tag{Eq 5}
\lambda_{max}N_t=N_t S_{1+}+f_{max} N_{t-a_p-1} \lambda _{max}^{a_p-1} S_0 S_{1+}^{(a_p-2)}
$$
Solving for $f_{max}$ gives the maximum theoretical fecundity as a function of $\lambda_{max}$, survival, and age at first parturition:

$$
\tag{Eq. 6}
f_{max} = \frac {\lambda_{max}^{a_p -1} - \lambda_{max}^{{a_p -2}}}  {S_0 {S_{1+}^{a_p-2}}}
$$
This is referred to as $p$ in Butterworth and Punt (1992). 

### Maximum net productivity level
To calculate the maximum net productivity level ($MNPL$) given $z$, we first calculate the sustainable yield $C$  as a function of bycatch mortality rate $E$.

$$
\tag{Eq. 7}
C = E \tilde B(E)\tilde P(E)
$$
where $\tilde B (E)$ is the normalized (to numbers at unfished equilibrium) recruitment when the bycatch rate is fixed at $E$ and $\tilde P (E)$ is the equilibrium number of "recruited" (age 1+) animals per calf/pup when the bycatch mortality rate is fixed at $E$. The normalized recruitment $\tilde B(E)$ is calculated as follows:

$$
\tag{Eq. 8}
\tilde B(E) = \bigg(1 - \frac{1-f_0 \tilde N(E)} {Af_0\tilde N (E)}\bigg)^{1/z}  \bigg(\frac{\tilde B(0)\tilde P(0)}{\tilde P(E)}\bigg)
$$
where $f_0 = \frac{1}{\tilde N(0)}$,  $\tilde N(E)$ is the number of animals at the age of first parturition and older (i.e., reproducing animals) per recruit at no-bycatch-mortality equilibrium, and $A$ is the Pella-Tomlinson resilience parameter ($A=\frac{f_{max}-f_0}{f_0}$ ; Punt 1999). $\tilde B(0)$ is assumed to be equal to 1, because all calculations are per-recruit.


The number of reproducing animals per recruit at exploitation rate $\boldsymbol{E}$ is the sum of adult animals per recruit $\tilde N_a$ from the age at first parturition $a_p$ to the plus group age $x$:

$$
\tag{Eq. 9}
\tilde N(E) = \sum_{a=a_p}^{x} \tilde N_{a}(E)
$$

We solve for the bycatch mortality rate at which $\tilde C$ is maximized (i.e., where $\frac{dC}{dE}$ is zero). This is $MSYR$, which is analogous to $F_{MSY}$ in fisheries.

The number of 1+ animals per recruit at bycatch mortality rate $E$, $\tilde P(E)$ is defined as:
$$
\tag{Eq. 10}
\tilde P(E)=\sum_{a=1}^{x} \tilde N_{a}(E)
$$

where $\tilde N_{0,a}(E)$ is the numbers per recruit at each age $a$ given a stable age structure:

$$
\tag{Eq. 11}
\tilde N_{0,a}(E) = 
\begin{cases}
1 &   a=0 \\
S_0[S_{1+}(1-E)]^{a-1} &    1\leq a<x \\
\frac{S_0[S_{1+}(1-E)]^{x-1}}{1-[S_{1+}(1-E)]} &    a=x \\
\end{cases}
$$ 

$MSYR$ is the value of $E$ for which the derivative of $C$ with respect to $E$ is zero, which we determined through numerical differentiation:

$$
\tag{Eq. 12}
\frac{dC}{dE} = \frac {C(E+0.001) - C(E-0.001)} {0.002}
$$

Then the relative abundance that corresponds to $MSYR$, $MNPL$, is determined by calculating the total 1+ population size at $MSYR$ relative to the equilibrium 1+ population size with no bycatch mortality:

$$
\tag{Eq. 13}
MNPL = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)\tilde B(0)} = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)} 
$$

### Parameterization
We assume that the population starts with a stable age structure in year 1 of the projection period (Eq. 11). The numbers at age at the start of projections correspond to a constant bycatch mortality rate $E$, which is calculated by solving the following equation for $E$: 

$$
\tag{Eq. 14}
\frac{\tilde B(E)\tilde P(E)}{\tilde P(0)}= \frac{N_0^{1+}}{K^{1+}}
$$

The initial depletion $\frac{N_0^{1+}}{K^{1+}}$ is based on the history of human-caused mortality for the population, which is provided by the user.

For cases where observation error is given for the initial population size, the starting abundance is drawn from a lognormal distribution with a standard deviation proportional to the observation CV.

#### Life history types
Each marine mammal life history type presented as a default option in this app corresponds to a unique value of calf/pup survival $S_0$, adult survival $S_{1+}$, age at first parturition, $a_p$, and intrinsic rate of population growth $\lambda_{max}$. These values are presented in Table 2. For computation purposes, we assumed that the plus group age $x$ is two years past the age at maturity ($x=a_p+1$). 

#### Compensation
We solve for the value of the degree of compensation $z$ that corresponds to the value of MNPL provided by the user. This involves solving the equation   $\tilde P(E^*) \tilde B(E^*) = MSYL * \tilde P(0)$ for $z$  where $E^*$ depends on $z$ as outlined above.



**Table 2.**
```{r echo = F}
lh.sources <- data.frame(stringsAsFactors=FALSE,
                           Type = c("Bowhead whale", "Bottlenose dolphin", "Humpback whale",
                                    "Phocid seal", "Fur seal", "Sea lion",
                                    "Porpoise", "Minke whale",
                                    "False killer whale/killer whale", "Pilot whale", "Right whale"),
                 Representative = c("Balaena mysticetus", "Tursiops truncatus",
                                    "Megaptera novaeangliae", "Phoca vitulina",
                                    "Arctocephalus pusillus pusillus",
                                    "Zalophus californianus", "Phocoena phocoena",
                                    "Balaenoptera bonaerensis", "Orcinus orca",
                                    "Globicephala macrorhynchus", "Eubalaena glacialis"),
                             S0 = c(0.944, 0.865, 0.9, 0.802, 0.77, 0.83, 0.8096, 0.84216,
                                    0.84744, 0.85008, 0.85536),
                         S1plus = c(0.99, 0.951, 0.95, 0.92, 0.88, 0.95, 0.92, 0.957, 0.963,
                                    0.966, 0.972),
                         AgeMat = c(17, 6, 10, 6, 3, 4, 3, 7, 9, 9, 8),
                   #PlusGroupAge = c(25, 10, 15, 8, 10, 5, 7, 9, 11, 11, 10),
                         Source = c("Punt et al. (2018) and references therein",
                                    "Punt et al. (2018) and references therein,
                                    except juvenile survival which is set to 0.88($S_{1+}$)", "Punt et al. (2018),
                                    Arso Civil et al. (2019), Speakman et al. (2010)",
                                    "Punt et al. (2018); Hastings et al. (2012)",
                                    "Punt et al. (2018); Butterworth et al. (1995)",
                                    "Punt et al. (2018); De Long et al. (2017)",
                                    "Moore unpublished (2019); G. Vikingsson pers. comm.; Olafsdottir et al. (2003)",
                                    "Moore unpublished (2019)",
                                    "Moore unpublished (2019)",
                                    "Moore unpublished (2019)",
                                    "Moore unpublished (2019)")
              )


x <- knitr::kable(format="html",
             col.names = c('Type','Representative','$S_0$','$S_{1+}$','Age at maturity','Source'),
             lh.sources, 
             escape = FALSE) %>% 
            kable_styling('bordered')
x <- column_spec(x, column = 2,italic = TRUE)
x
```


<!-- ### Citations -->
<!-- Arso Civil, M., Cheney, B., Quick, N.J., Islas-Villanueva, V., Graves, J.A., Janik, V.M., et al. (2019). Variations in age- and sex-specific survival rates help explain population trend in a discrete marine mammal population. Ecology and Evolution, 9, 533–544. -->

<!-- Breiwick, J.M., Eberhardt, L.L. & Braham, H.W. (1984). Population Dynamics of Western Arctic Bowhead Whales ( Balaena mysticetus ). Canadian Journal of Fisheries and Aquatic Sciences, 41, 484–496. -->

<!-- Butterworth, D.S. & Punt, A.E. (1992). The Scientific Committee “...Agreed That The MSY Rate Would Most Likely Lie Between 1 and 4%” - But Which MSY Rate? Reports of the International Whaling Commission, 42, 583–591. -->

<!-- Butterworth, D.S., Punt, A.E., Oosthuizen, W.H., Wickens, P.A. (1995) The effects of future consumption by the Cape fur seal on catches and catch rates of the Cape hakes. 3. Modelling the dynamics of the Cape fur seal Arctocephalus pusillus pusillus. South African Journal of Marine Science, 16, 161–183. -->

<!-- DeLong, R.L., Melin, S.R., Laake, J.L., Morris, P., Orr, A.J. & Harris, J.D. (2017). Age- and sex-specific survival of California sea lions (Zalophus californianus) at San Miguel Island, California. Marine Mammal Science, 33, 1097–1125. -->

<!-- Dillingham, P.W., Moore, J.E., Fletcher, D., Cortés, E., Curtis, K.A., James, K.C., et al. (2016). Improved estimation of intrinsic growth rmax for long-lived species: integrating matrix models and allometry. Ecological Applications, 26, 322–333. -->

<!-- Moore 2019. Unpublished estimate following the methods of Dillingham et al. 2016. -->

<!-- Ólafsdóttir, D., Víkingsson, G. A., Halldórsson, S. D., & Sigurjónsson, J. (2003). Growth and reproduction in harbour porpoises (Phocoena phocoena) in Icelandic waters. NAMMCO Scientific Publications, 5, 195–210. -->

<!-- Punt, A. E. (1999). Annex R: A full description of the standard Baleen II model and some variants thereof. Division of Marine Research, CSIRO Marine Laboratories, Hobart, Australia. -->

<!-- Punt, A.E., Moreno, P., Brandon, J.R. & Mathews, M.A. (2018). Conserving and recovering vulnerable marine species: a comprehensive evaluation of the US approach for marine mammals. ICES Journal of Marine Science, 75, 1813–1831. -->
---
title: 'Sobre el modelo de proyección'
author: 'Siple MC and Punt AE'
date: 'última actualización: October 6, 2020'
output:
  html_document:
    highlight: tango
    theme: flatly
  word_document: default
---

<!-- NOTE: in order to include this html in the shiny app, delete everything in the knitted html before <body> and after </body> -->

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# devtools::install_github("haozhu233/kableExtra")
library(kableExtra)
```


### Modelo de población

Las proyecciones poblacionales se generan utilizando un modelo estructurado por edad de un solo sexo. En este modelo de proyección, el número de crías nacidas cada año es denso dependiente, con la extensión de denso dependencia como función del número de adultos maduros $\tilde N$, la fecundidad (tasa de embarazo) en equilibrio previo a la explotación $f_0$, la tasa de fecundidad teórica máxima $f_{max}$, el grado de compensación $z$, y  la abundancia de individuos de edades 1+ ($N_{y+1}^{1+}$) en relación a la capacidad de carga $K^{1+}$.

El número de individuos de 1 año de edad depende de la supervivencia de crías $S_0$, de la tasa de supervivencia de los adultos, $S_{1+}$ y de remociones debido a la mortalidad por captura incidental $C_y$ (todos los símbolos están definidos en la Tabla 1).


$$
\tag{Eq 1}
N_{y+1,a} = 
\begin{cases}
\tilde N_{y+1}\bigg(f_0 + (f_{max}-f_0)\bigg[1-\bigg(\frac{N_{y+1}^{1+}}{K^{1+}}\bigg)^z\bigg]\bigg) &   a=0\\ 
N_{y,0}S_0 &    a=1 \\
(N_{y,a-1}-C_{y,a-1})S_{1+}  & 2\leq a<x \\
(N_{y,x-1}-C_{y,x-1})S_{1+} + (N_{y,x}-C_{y,x})S_{1+} & a=x \\
\end{cases}
$$

dónde $N_{y,a}$ es  el número de animales de edad $a$ a principios del año $y$. El número de animales de edad $a$ muriendo  producto de las capturas incidentales durante el año $y$, $C_{y,a}$, es removida uniformemente desde el componente 1+ de la población (se supone que las crías no mueren debido a la captura incidental), es decir:

$$
\tag{Eq 2}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
\frac{C_y N_{y,a}}{N_y^{1+}} & a>0 \\
\end{cases}
$$

Dado este supuesto, nuestro modelo no caracterizará adecuadamente los escenarios si la mortalidad por captura incidental ocurre predominantemente en las crías.

Si el usuario especifica una tasa constante de mortalidad por captura incidental, los individuos se eliminan de la población de acuerdo a una tasa de mortalidad por captura incidental $E_y$ y vulnerabilidad (que se supone constante a través del tiempo, y uniforme en individuos de edad 1+):

$$
\tag{Eq 3}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
N_{y,a} E_{y} & a>0 \\
\end{cases}
$$

**Tabla 1.**  Símbolos incluidos en el modelo de proyección.

```{r echo = F}
Símbolo=c(#User-specified parameters
  "$S_0$",
  "$S_{1+}$",
  "$x$",
  "$\\lambda_{max}$",
  "$a_p$",
  "$E_y$",
  # Parameters derived from user-specified parameters
  "$f_0$",
  "$f_{max}$",
  "$K^{1+}$",
  "$z$",
  # Derived variables
  "$\\tilde N_{y,a}$", # need double dash for using in kableExtra
  "$N_{y,a}$",
  "$N_y^{1+}$",
  "$C_{y,a}$")

Descripción=c(#User-specified parameters
  
  "Supervivencia de cachorro o ballenato",
  "Supervivencia de individuos de edad 1+",
  "Edad del grupo plus",
  "Máxima tasa constante de aumento (tasa de crecimiento de la población)",
  "Edad al primer parto",
  "Tasa de mortalidad por captura incidental en el año $y$ (especificado por el usuario o calculado a partir de la mortalidad por captura incidental total especificada por el usuario)", 
  # Parameters derived from user-specified parameters
  "Fecundidad (tasa de preñez) en equilibrio previo a la explotación",
  "Fecundidad teórica máxima (tasa de preñez)",
  "Capacidad de carga en términos de 1+ componente de la población",
  "Grado de compensación",
  # Derived variables
  "Número de animales maduros de edad $a$ a principios del año $y$",
  "Número de animales de edad $a$ a principios del año $y$",
  "Número de animales de edad 1 y mayores a principios del año $y$",
  "Mortalidad por captura incidental de animales de edad $a$ en año $y$")

Symbol_definitions <- data.frame(cbind(Símbolo,
                                       Descripción))
knitr::kable(format="html", #html  - change to pandoc for ms word 
             Symbol_definitions, 
             escape = FALSE) %>% 
  kable_styling('bordered') %>%
  pack_rows('Parámetros especificados por el usuario', 1, 6,
            label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
  pack_rows('Parámetros derivados de los parámetros especificados por el usuario', 7, 10,
            label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
  pack_rows('Parámetros derivados', 11, 14, 
            label_row_css = "font-style: italic;border-bottom: 1px solid black;")
```

### Medidas de desempeño
Hay tres medidas principales de rendimiento incluidas en los resultados:

- la probabilidad de que la abundancia de animales de edad 1+ exceda el Nivel de Productividad Neta Máxima (MNPL por sus siglas en inglés) después de 50 y 100 años;
- la abundancia de animales de edad 1+ en comparación con la capacidad de carga 1+ después de 10, 20 y 50 años; y
- la abundancia de animales de edad 1+ en comparación con la abundancia de animales de edad 1+ en un escenario sin captura incidental después de 10, 20 y 50 años.

Estas medidas de rendimiento están relacionadas con la recuperación de la población. Wade (1998) identifica un 'estándar' de rendimiento de un 95% de probabilidad de recuperación a MNPL después de 100 años para una población inicialmente al 30% de su capacidad de carga, con un MNPL de 0.5K. MNPL es el límite inferior para la Población Óptima Sostenible (OSP de su sigla en inglés), que se define como el "número de animales que dará como resultado la productividad máxima de la población o la especie, teniendo en cuenta la capacidad de carga del hábitat y la salud del ecosistema del cual forman un elemento constituyente ”(16 USCS § 1362 (9)). En el esquema de manejo de EE. UU., las poblaciones de mamíferos marinos se consideran agotadas cuando están por debajo de OSP (16 USCS § 1362 (1B)).

Los valores individuales informados en las tablas son las medianas del número de proyecciones.


### Resolviendo la tasa de captura incidental

Resolvimos el nivel de captura incidental que resultaría en la recuperación de una proporción dada de K en un período específico de tiempo. La probabilidad de reconstrucción se calcula como la proporción de simulaciones en las que la población se ha recuperado a la meta de reconstrucción para el año de reconstrucción. La pestaña "Resolver para captura incidental" utiliza la búsqueda de raíz (la función `uniroot()` en nuestro código) para encontrar la tasa de mortalidad de captura incidental que minimiza la diferencia entre la probabilidad de recuperación deseada y la probabilidad de recuperación bajo la tasa de mortalidad de captura incidental $E$.


### Fecundidad teórica máxima

La tasa de fecundidad teórica máxima se calcula en función del tamaño de la población cuando la población está creciendo a su tasa máxima de crecimiento de la población $\lambda_{max}$. Se pueden encontrar partes de esta derivación en Butterworth y Punt (1992) y Breiwick et al. (1984) Comenzando con el modelo de población para el componente femenino de la población:

$$
\tag{Eq 4}
N_{t+1} = N_t(1-E)S_{1+}+N_{t-a_p}f q_f S_0 (S_{1+})^{a_p -2}(1-E)^{a_p-1-a_r}
$$

dónde $a_p$ es la edad en el primer parto (se supone que es un año después de la edad de madurez sexual), $f$ es la fecundidad (a veces expresada como el producto de la tasa de preñez y la supervivencia de las crías), $q_f$ es la proporción de crias que son hembras y $a_r$ es la edad a la que los animales sufren por primera vez alguna mortalidad por captura incidental (en nuestro caso, esto es igual a 1 año). Cuando la población crece a su ritmo más rápido, la pesca es cero, y $N_{t+1}=N_t \lambda_{max}$. En nuestro caso $q_f = 1$ porque estamos modelando a todos los adultos, en lugar de solo hembras maduras. Por lo tanto, la ecuación anterior se convierte en:

$$
\tag{Eq 5}
\lambda_{max}N_t=N_t S_{1+}+f_{max} N_{t-a_p-1} \lambda _{max}^{a_p-1} S_0 S_{1+}^{(a_p-2)}
$$
Resolviendo para $f_{max}$ dada la máxima fecundidad teórica en función de $\lambda_{max}$, supervivencia y edad en el primer parto:

$$
\tag{Eq. 6}
f_{max} = \frac {\lambda_{max}^{a_p -1} - \lambda_{max}^{{a_p -2}}}  {S_0 {S_{1+}^{a_p-2}}}
$$
Esto se conoce como $p$ en Butterworth y Punt (1992).


### Máximo nivel de productividad neta

Para calcular el nivel máximo de productividad neta ($MNPL$) dato $z$, primero calculamos el rendimiento sostenible por recluta $C$ en función de la tasa de mortalidad por captura incidental $E$.

$$
\tag{Eq. 7}
C = E \tilde B(E)\tilde P(E)
$$
dónde $\tilde B (E)$ es el reclutamiento normalizado cuando la tasa de captura incidental se fija en $E$ y $\tilde P (E)$ es el número de equilibrio de animales "reclutados" (edad 1+) por cría cuando la tasa de mortalidad por captura incidental se fija en $E$. El reclutamiento normalizado $\tilde B(E)$ se calcula de la siguiente manera:

$$
\tag{Eq. 8}
\tilde B(E) = \bigg(1 - \frac{1-f_0 \tilde N(E)} {Af_0\tilde N (E)}\bigg)^{1/z}  \bigg(\frac{\tilde B(0)\tilde P(0)}{\tilde P(E)}\bigg)
$$
dondé $f_0 = \frac{1}{\tilde N(0)}$,  $\tilde N(E)$ es el número de animales a la edad del primer parto y mayores (es decir, animales reproductores) por recluta al equilibrio sin mortalidad incidental, y $A$ es el parámetro de resiliencia de Pella-Tomlinson (($A=\frac{f_{max}-f_0}{f_0}$ ; Punt 1999). $\tilde B(0)$ asumiendo que es igual a 1, porque todos los cálculos son por recluta.

El número de animales reproductores por recluta a la tasa de explotación $\boldsymbol{E}$ es la suma de animales adultos por recluta $\tilde N_a$ desde la edad del primer parto $a_p$ a la edad del grupo plus $x$:

$$
\tag{Eq. 9}
\tilde N(E) = \sum_{a=a_p}^{x} \tilde N_{a}(E)
$$

Resolvemos la tasa de mortalidad por captura incidental a la que $\tilde C$ es maximizada (es decir, donde $\frac{dC}{dE}$ es cero) Esto es $MSYR$, que es análogo a $F_{MSY}$ en pesquerías.

El número de 1+ animales por recluta a la tasa de mortalidad por captura incidental $E$, $\tilde P(E)$ se define como:
$$
\tag{Eq. 10}
\tilde P(E)=\sum_{a=1}^{x} \tilde N_{a}(E)
$$

dónde $\tilde N_{0,a}(E)$ son los números por recluta en cada edad $a$  dada una estructura de edad estable:

$$
\tag{Eq. 11}
\tilde N_{0,a}(E) = 
\begin{cases}
1 &   a=0 \\
S_0[S_{1+}(1-E)]^{a-1} &    1\leq a<x \\
\frac{S_0[S_{1+}(1-E)]^{x-1}}{1-[S_{1+}(1-E)]} &    a=x \\
\end{cases}
$$ 

$MSYR$ es el valor de $E$ para el cual la derivada de $C$ con respecto a $E$ es cero, que determinamos mediante diferenciación numérica:

$$
\tag{Eq. 12}
\frac{dC}{dE} = \frac {C(E+0.001) - C(E-0.001)} {0.002}
$$

Entonces, la abundancia relativa que corresponde a $MSYR$, $MNPL$, se determina calculando el tamaño total de la población 1+ en $MSYR$ en relación con el tamaño de la población de equilibrio 1+ sin mortalidad por captura incidental:

$$
\tag{Eq. 13}
MNPL = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)\tilde B(0)} = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)} 
$$

### Parametrización
Supongamos que la población comienza con una estructura de edad estable en el año 1 del período de proyección (Ec. 11). Los números a la edad al comienzo de las proyecciones corresponden a una tasa constante de mortalidad por captura incidental $E$, que se calcula resolviendo la siguiente ecuación para $E$: 

$$
\tag{Eq. 14}
\frac{\tilde B(E)\tilde P(E)}{\tilde P(0)}= \frac{N_0^{1+}}{K^{1+}}
$$

El agotamiento inicial $\frac{N_0^{1+}}{K^{1+}}$  se basa en la historia de mortalidad causada por humanos para la población, que es proporcionada por el usuario.

Para los casos en que se da un error de observación para la población, la abundancia inicial se extrae de una distribución log normal con una desviación estándar proporcional al CV de observación.

#### Tipos de historia de vida
Cada tipo de historia de vida de mamíferos marinos presentado como una opción predeterminada en esta aplicación corresponde a un valor único de supervivencia de las crías $S_0$, sobrevivencia adulta $S_{1+}$, edad al primer parto, $a_p$, y tasa intrínseca de crecimiento de la población $\lambda_{max}$. Estos valores se presentan en la Tabla 2. Para fines de cálculo, supusimos que la edad del grupo más $x$ es dos años después de la edad de madurez ($x=a_p+1$).

#### Compensación
Resolvemos el valor del grado de compensación $z$ que corresponde al valor de MNPL proporcionado por el usuario. Esto implica resolver la ecuación $\tilde P(E^*) \tilde B(E^*) = MSYL * \tilde P(0)$ para $z$  donde $E^*$ depende de $z$  como se describe arriba.



**Tabla 2.**
```{r echo = F}
lh.sources <- data.frame(stringsAsFactors=FALSE,
                         Tipo = c("Ballena de Groenlandia o ballena boreal",
                                  "Delfín nariz de botella",
                                  "Ballena jorobada",
                                  "Foca común",
                                  "Lobo fino o de dos pelos",
                                  "Lobo común o de un pelo",
                                  "Marsopa común",
                                  "Ballena minke antártica",
                                  "Orca",
                                  "Calderón de aletas cortas",
                                  "Ballena franca"),
                         Representante = c("Balaena mysticetus",
                                           "Tursiops truncatus",
                                           "Megaptera novaeangliae",
                                           "Phoca vitulina",
                                           "Arctocephalus pusillus pusillus",
                                           "Zalophus californianus",
                                           "Phocoena phocoena",
                                           "Balaenoptera bonaerensis",
                                           "Orcinus orca",
                                           "Globicephala macrorhynchus",
                                           "Eubalaena glacialis"),
                         S0 = c(0.944, 0.865, 0.9, 0.802, 0.77, 0.83, 0.8096, 0.84216,
                                0.84744, 0.85008, 0.85536),
                         S1plus = c(0.99, 0.951, 0.95, 0.92, 0.88, 0.95, 0.92, 0.957, 0.963,
                                    0.966, 0.972),
                         AgeMat = c(17, 6, 10, 6, 3, 4, 3, 7, 9, 9, 8),
                         #PlusGroupAge = c(25, 10, 15, 8, 10, 5, 7, 9, 11, 11, 10),
                         Source = c("Punt et al. (2018) and references therein",
                                    "Punt et al. (2018) and references therein,
                                    except juvenile survival which is set to 0.88($S_{1+}$)", "Punt et al. (2018),
                                    Arso Civil et al. (2019), Speakman et al. (2010)",
                                    "Punt et al. (2018); Hastings et al. (2012)",
                                    "Punt et al. (2018); Butterworth et al. (1995)",
                                    "Punt et al. (2018); De Long et al. (2017)",
                                    "Moore unpublished (2019); G. Vikingsson pers. comm.; Olafsdottir et al. (2003)",
                                    "Moore unpublished (2019)",
                                    "Moore unpublished (2019)",
                                    "Moore unpublished (2019)",
                                    "Moore unpublished (2019)")
)


x <- knitr::kable(format="html",
                  col.names = c('Tipo','Representante','$S_0$','$S_{1+}$','Edad de madurez','Referencias'),
                  lh.sources, 
                  escape = FALSE) %>% 
  kable_styling('bordered')
x <- column_spec(x, column = 2,italic = TRUE)
x
```


<!-- ### Citations -->
<!-- Arso Civil, M., Cheney, B., Quick, N.J., Islas-Villanueva, V., Graves, J.A., Janik, V.M., et al. (2019). Variations in age- and sex-specific survival rates help explain population trend in a discrete marine mammal population. Ecology and Evolution, 9, 533–544. -->

<!-- Breiwick, J.M., Eberhardt, L.L. & Braham, H.W. (1984). Population Dynamics of Western Arctic Bowhead Whales ( Balaena mysticetus ). Canadian Journal of Fisheries and Aquatic Sciences, 41, 484–496. -->

<!-- Butterworth, D.S. & Punt, A.E. (1992). The Scientific Committee “...Agreed That The MSY Rate Would Most Likely Lie Between 1 and 4%” - But Which MSY Rate? Reports of the International Whaling Commission, 42, 583–591. -->

<!-- Butterworth, D.S., Punt, A.E., Oosthuizen, W.H., Wickens, P.A. (1995) The effects of future consumption by the Cape fur seal on catches and catch rates of the Cape hakes. 3. Modelling the dynamics of the Cape fur seal Arctocephalus pusillus pusillus. South African Journal of Marine Science, 16, 161–183. -->

<!-- DeLong, R.L., Melin, S.R., Laake, J.L., Morris, P., Orr, A.J. & Harris, J.D. (2017). Age- and sex-specific survival of California sea lions (Zalophus californianus) at San Miguel Island, California. Marine Mammal Science, 33, 1097–1125. -->

<!-- Dillingham, P.W., Moore, J.E., Fletcher, D., Cortés, E., Curtis, K.A., James, K.C., et al. (2016). Improved estimation of intrinsic growth rmax for long-lived species: integrating matrix models and allometry. Ecological Applications, 26, 322–333. -->

<!-- Moore 2019. Unpublished estimate following the methods of Dillingham et al. 2016. -->

<!-- Ólafsdóttir, D., Víkingsson, G. A., Halldórsson, S. D., & Sigurjónsson, J. (2003). Growth and reproduction in harbour porpoises (Phocoena phocoena) in Icelandic waters. NAMMCO Scientific Publications, 5, 195–210. -->

<!-- Punt, A. E. (1999). Annex R: A full description of the standard Baleen II model and some variants thereof. Division of Marine Research, CSIRO Marine Laboratories, Hobart, Australia. -->

<!-- Punt, A.E., Moreno, P., Brandon, J.R. & Mathews, M.A. (2018). Conserving and recovering vulnerable marine species: a comprehensive evaluation of the US approach for marine mammals. ICES Journal of Marine Science, 75, 1813–1831. -->
---
title: 'Sur le Modèle'
author: 'Siple MC and Punt AE'
date: 'Dernière mise à jour: October 6, 2020'
output:
  html_document:
    highlight: tango
    theme: flatly
  word_document: default
---

<!-- NOTE: in order to include this html in the shiny app, delete everything in the knitted html before <body> and after </body> -->

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# devtools::install_github("haozhu233/kableExtra")
library(kableExtra)
```


### Modèle de population

Les projections de la population sont générées à l'aide d'un modèle unisexe structuré en âge. Dans ce modèle de projection, le nombre de baleineaux ou de nouveau-nés chaque année dépend de la densité, l'intensité de la densité dépendance étant fonction du nombre d'adultes matures $\tilde N$, de la fécondité (taux de gestation) à l'équilibre de pré-exploitation $f_0$, du taux maximum de fécondité théorique  $f_{max}$, du degré de compensation $z$, et de l'abondance des individus d'âges 1+ ($N_{y+1}^{1+}$) par rapport à la capacité de charge $K^{1+}$.

Le nombre d'individus d'âge 1 est fonction de la survie des baleineaux/nouveaux nés $S_0$, et du taux de survie des adultes, $S_{1+}$ et des prélèvements dus à la mortalité des prises accessoires $C_y$ (tous les symboles sont définis dans le Tableau 1).

$$
\tag{Eq 1}
N_{y+1,a} = 
\begin{cases}
\tilde N_{y+1}\bigg(f_0 + (f_{max}-f_0)\bigg[1-\bigg(\frac{N_{y+1}^{1+}}{K^{1+}}\bigg)^z\bigg]\bigg) &   a=0\\ 
N_{y,0}S_0 &    a=1 \\
(N_{y,a-1}-C_{y,a-1})S_{1+}  & 2\leq a<x \\
(N_{y,x-1}-C_{y,x-1})S_{1+} + (N_{y,x}-C_{y,x})S_{1+} & a=x \\
\end{cases}
$$

où $N_{y,a}$ est le nombre d'animaux  d'âge $a$ au début de l'année $y$. Le nombre d'animaux d'âge $a$ mourant en raison de la mortalité due aux prises accessoires au cours de l'année $y$, $C_{y,a}$, est retiré uniformément de la catégorie d'âge 1+ de la population (les veaux et les petits sont supposés ne pas mourir en raison des prises accessoires), c'est-à-dire:

$$
\tag{Eq 2}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
\frac{C_y N_{y,a}}{N_y^{1+}} & a>0 \\
\end{cases}
$$

Compte tenu de cette hypothèse, notre modèle ne caractérisera pas correctement les scénarios où la mortalité due aux prises accessoires se produit principalement chez les baleineaux/nouveaux-nés.

Si l'utilisateur spécifie un taux constant de mortalité par prise accessoire, les individus sont retirés de la population en fonction d'un taux de mortalité par prise accessoire $E_y$ et de la vulnérabilité (qui est supposée constante dans le temps, et uniforme sur les individus d'âge 1+) :

$$
\tag{Eq 3}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
N_{y,a} E_{y} & a>0 \\
\end{cases}
$$

**Tableau 1.** Symboles inclus dans le modèle de projection.

```{r echo = F}
Symbole=c(#User-specified parameters
  "$S_0$",
  "$S_{1+}$",
  "$x$",
  "$\\lambda_{max}$",
  "$a_p$",
  "$E_y$",
  # Parameters derived from user-specified parameters
  "$f_0$",
  "$f_{max}$",
  "$K^{1+}$",
  "$z$",
  # Derived variables
  "$\\tilde N_{y,a}$", # need double dash for using in kableExtra
  "$N_{y,a}$",
  "$N_y^{1+}$",
  "$C_{y,a}$")
Description=c(#User-specified parameters
  
  "Survie du nouveau-né ou du baleineau",
  "Survie des individus âgés de 1 an et plus",
  "Âge du groupe-plus",
  "Taux maximum d'accroissement constant (taux de croissance de la population)",
  "L'âge à la première mise bas",
  "Taux de mortalité par prises accessoires durant l'année $y$ (spécifié par l'utilisateur ou calculé à partir de la mortalité totale des prises accessoires spécifiée par l'utilisateur)", 
  # Parameters derived from user-specified parameters
  "Fécondité (taux de grossesse) à l'équilibre pré-exploitation",
  "Fécondité maximale théorique (taux de grossesse)",
  "Capacité de charge exprimée en termes de composante 1+ de la population",
  "Degré de compensation",
  # Derived variables
  "Nombre d'animaux matures d'âge $a$ au début de l'année $y$",
  "Nombre d'animaux d'âge $a$ au début de l'année $y$",
  "Nombre d'animaux d'âge 1 et plus au début de l'année $y$",
  "Mortalité due aux prises accessoires d'animaux d’âges $a$ durant l'année $y$")
Symbol_definitions <- data.frame(cbind(Symbole,
                                       Description))
knitr::kable(format="html", #html  - change to pandoc for ms word 
             Symbol_definitions, 
             escape = FALSE) %>% 
  kable_styling('bordered') %>%
  pack_rows("Paramètres spécifiés par l'utilisateur", 1, 6, label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
  pack_rows("Paramètres dérivés des paramètres spécifiés par l'utilisateur", 7, 10, label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
  pack_rows("Variables dérivées", 11, 14, label_row_css = "font-style: italic;border-bottom: 1px solid black;")
```

### Mesures de performance
Trois principales mesures de performance sont incluses dans les résultats :

- la probabilité que l'abondance des animaux d'âge 1+ dépasse le Niveau de Productivité Net Maximal (NPNM) après 50 et 100 ans ;
- l'abondance des animaux d'âge 1+ par rapport à la capacité de charge des animaux d'âge 1+ après 10, 20 et 50 ans; et
-	l'abondance des animaux d'âge 1+ par rapport à l'abondance des animaux d'âge 1+ dans un scénario de non prise accidentelle après 10, 20 et 50 ans.

Ces mesures de performance sont toutes associées à la reconstitution de la population. Wade (1998) a établi une "norme" de performance d'une probabilité de 95 % de rétablissement de la NPNM après 100 ans pour une population se trouvant initialement à 30 % de sa capacité de charge, avec une NPNM de 0,5 K. La NPNM est la limite inférieure de la population optimale durable (POD), qui est définie comme le "nombre d'animaux qui permettra d'obtenir la productivité maximale de la population ou de l'espèce, en tenant compte de la capacité de charge de l'habitat et de la santé de l'écosystème dont ils sont un élément constitutif" (16 USCS § 1362 (9)). Dans le système de gestion américain, les stocks de mammifères marins sont considérés comme épuisés lorsqu'ils sont inférieurs à la POD (16 USCS § 1362 (1B)).

Les valeurs uniques indiquées dans les tableaux sont les médianes du nombre de projections.


### Résoudre le problème du taux de prises accessoires

Nous avons résolu le problème du niveau de prises accessoires qui aboutirait à la reconstitution d'une proportion donnée de K dans un intervalle de temps donné. La probabilité de reconstitution est calculée comme étant la proportion de simulations dans lesquelles la population a atteint l'objectif de reconstitution au cours de l'année de reconstitution. L'onglet "Résoudre le problème des prises accessoires" utilise l'algorithme de recherche d'un zéro d'une fonction (la fonction `uniroot()` dans notre code) pour trouver le taux de mortalité des prises accessoires qui minimise la différence entre la probabilité de reconstitution souhaitée et la probabilité de reconstitution avec un taux de mortalité des prises accessoires de $E$.


### Fécondité théorique maximale
Le taux de fécondité théorique maximal est calculé sur la base de la taille de la population à son taux maximal de croissance démographique $\lambda_{max}$. Des parties de cette dérivation mathematique se trouvent dans Butterworth et Punt (1992) et Breiwick et al. (1984). Commençons par le modèle de population pour la composante féminine de la population :
$$
\tag{Eq 4}
N_{t+1} = N_t(1-E)S_{1+}+N_{t-a_p}f q_f S_0 (S_{1+})^{a_p -2}(1-E)^{a_p-1-a_r}
$$

où $a_p$ est l'âge à la première mise bas (supposée être un an après l'âge à maturité), $f$ est la fécondité (parfois exprimée comme le produit du taux de grossesse et de la survie des nouveau-nés), $q_f$ est la proportion de baleineaux/nouveaux-nés qui sont des femelles et $a_r$ est l'âge auquel les mammifères subissent pour la première fois une mortalité due à des prises accessoires (dans notre cas, cela est égal à 1 an). Lorsque la population croît à son rythme le plus rapide, la pêche est nulle et $N_{t+1}=N_t \lambda_{max}$. Dans notre cas, $q_f = 1$ parce que nous modélisons tous les adultes, et non uniquement femelles matures. Ainsi, l'équation ci-dessus devient :
$$
\tag{Eq 5}
\lambda_{max}N_t=N_t S_{1+}+f_{max} N_{t-a_p-1} \lambda _{max}^{a_p-1} S_0 S_{1+}^{(a_p-2)}
$$
Le calcul de $f_{max}$ permet d’obtenir la fécondité maximale théorique comme une fonction de $\lambda_{max}$, de la survie et de l'âge à la première mise bas :

$$
\tag{Eq. 6}
f_{max} = \frac {\lambda_{max}^{a_p -1} - \lambda_{max}^{{a_p -2}}}  {S_0 {S_{1+}^{a_p-2}}}
$$
Ceci est désigné sous le nom de $p$ dans Butterworth et Punt (1992).

### Niveau de productivité net maximal
Pour calculer le niveau de productivité net maximale ($NPNM$) sachany $z$, nous calculons d'abord le rendement durable par recrue $C$ en fonction du taux de mortalité des prises accessoires $E$.

$$
\tag{Eq. 7}
C = E \tilde B(E)\tilde P(E)
$$
où $\tilde B (E)$ est le recrutement standardisé lorsque le taux de prises accessoires est fixé à $E$ et $\tilde P (E)$ est le nombre d'animaux "recrutés" (âge 1+) par baleineaux/nouveaux nés lorsque le taux de mortalité des prises accessoires est fixé à $E$. Le recrutement standardisé $\tilde B(E)$ est calculé de la manière suivante :

$$
\tag{Eq. 8}
\tilde B(E) = \bigg(1 - \frac{1-f_0 \tilde N(E)} {Af_0\tilde N (E)}\bigg)^{1/z}  \bigg(\frac{\tilde B(0)\tilde P(0)}{\tilde P(E)}\bigg)
$$
où $f_0 = \frac{1}{\tilde N(0)}$, $\tilde N(E)$ est le nombre d'animaux à l'âge de la première mise bas et plus âgés (c'est-à-dire les animaux reproducteurs) par recrue à l'équilibre de la mortalité sans prise accessoire, et $A$ est le paramètre de résilience de Pella-Tomlinson ($A=\frac{f_{max}-f_0}{f_0}$ ; Punt 1999). $\tilde B(0)$ est considéré comme égal à 1, car tous les calculs sont effectués par recrue.

Le nombre d'animaux reproducteurs par recrue, pour un taux d'exploitation $E$, est égal à la somme des animaux adultes par recrue $\tilde N_a$ de l'âge à la première mise bas $a_p$ jusqu’ à l'âge du groupe plus $x$ :

$$
\tag{Eq. 9}
\tilde N(E) = \sum_{a=a_p}^{x} E \tilde N_{0,a}(E)
$$

Nous déterminons le taux de mortalité des prises accessoires auquel $\tilde C$ est maximal (c'est-à-dire où $\frac{dC}{dE}$ est égal à zéro). C'est le RMDT, qui est analogue à $F_{MSY}$ dans les pêcheries.

Le nombre d'animaux 1+ par recrue, pour un taux de mortalité des prises accessoires $\boldsymbol{E}$, $\tilde P(E)$ est défini de la façon suivante:
$$
\tag{Eq. 10}
\tilde P(E)=\sum_{a=1}^{x} \tilde N_{0,a}(E)
$$

où $\tilde N_{0,a}(E)$ est le nombre par recrue à chaque âge $a$ pour une structure d'âge stable donnée:

$$
\tag{Eq. 11}
\tilde N_{0,a}(E) = 
\begin{cases}
1 &   a=0 \\
S_0[S_{1+}(1-E)]^{a-1} &    1\leq a<x \\
\frac{S_0[S_{1+}(1-E)]^{x-1}}{1-[S_{1+}(1-E)]} &    a=x \\
\end{cases}
$$ 

$MSYR$ est la valeur de $E$ pour laquelle la dérivée de $C$ par rapport à $E$ est égale à zéro, déterminée par différenciation numérique :

$$
\tag{Eq. 12}
\frac{dC}{dE} = \frac {C(E+0.001) - C(E-0.001)} {0.002}
$$

Ensuite, l'abondance relative qui correspond au $MSYR$, $MNPL$, est déterminée en calculant la taille totale de la population 1+ au $RMDT$ par rapport à la taille de la population 1+ à l'équilibre, sans mortalité due aux prises accessoires:

$$
\tag{Eq. 13}
MNPL = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)\tilde B(0)} = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)} 
$$

### Paramétrisation
Nous considérons que la population commence avec une structure d'âge stable la première année de la période de projection (équation 11). Les valeurs à l'âge au début des projections correspondent à un taux de mortalité par prise accessoire constant $E$, qui est calculé en résolvant l'équation suivante pour $E$ :

$$
\tag{Eq. 14}
\frac{\tilde B(E)\tilde P(E)}{\tilde P(0)}= \frac{N_0^{1+}}{K^{1+}}
$$

L'épuisement initial $\frac{N_0^{1+}}{K^{1+}}$ est basé sur l'historique de la mortalité causée par l'homme au sein de la population, qui est fourni par l'utilisateur.

Dans les cas où l'erreur d'observation est donnée pour la population, l'abondance de départ est déterminée à partir d'une distribution lognormale avec un écart-type proportionnel au CV d'observation.

#### Types d'histoire de vie
Chaque type d'histoire de vie de mammifère marin présenté comme option par défaut dans cette application correspond à une valeur unique de la survie du baleineau/nouveau né $S_0$, de la survie de l'adulte $S_{1+}$, de l'âge à la première mise bas, $a_p$ et du taux intrinsèque de croissance de la population $\lambda_{max}$. Ces valeurs sont présentées dans le tableau 2. Pour des raisons de calcul, nous avons supposé que l'âge du groupe plus $x$ est supérieur de deux ans à l'âge à la maturité ($x=a_p+1$). 

#### Compensation

Nous déterminons la valeur du degré de compensation $z$ qui correspond à la valeur du NPNM fourni par l'utilisateur. Cela implique de résoudre l'équation $\tilde P(E^*) \tilde B(E^*) = MSYL * \tilde P(0)$ pour $z$ où $E^*$ dépend de $z$ comme indiqué ci-dessus.



**Tableau 2.**
```{r echo = F}
lh.sources <- data.frame(stringsAsFactors=FALSE,
                         Type = c("Baleine boréale",
                                  "Grand dauphin",
                                  "Baleine à bosse",
                                  "Phoque commun",
                                  "Otarie à fourrure",
                                  "Lion de mer",
                                  "Marsouin",
                                  "Petit rorqual",
                                  "Fausse orque/orque",
                                  "Globicéphales",
                                  "Baleine franche"),
                         Representant = c("Balaena mysticetus",
                                          "Tursiops truncatus",
                                          "Megaptera novaeangliae",
                                          "Phoca vitulina",
                                          "Arctocephalus pusillus pusillus",
                                          "Zalophus californianus",
                                          "Phocoena phocoena",
                                          "Balaenoptera bonaeresis",
                                          "Orcinus orca",
                                          "Globicephala macrorhynchus",
                                          "Eubalaena glacialis"),
                         S0 = c(0.944, 0.865, 0.9, 0.802, 0.77, 0.83, 0.8096, 0.84216,
                                0.84744, 0.85008, 0.85536),
                         S1plus = c(0.99, 0.951, 0.95, 0.92, 0.88, 0.95, 0.92, 0.957, 0.963,
                                    0.966, 0.972),
                         AgeMat = c(17, 6, 10, 6, 3, 4, 3, 7, 9, 9, 8),
                         #PlusGroupAge = c(25, 10, 15, 8, 10, 5, 7, 9, 11, 11, 10),
                         Source = c("Punt et al. (2018) and references therein",
                                    "Punt et al. (2018) and references therein,
                                    except juvenile survival which is set to 0.88($S_{1+}$)", "Punt et al. (2018),
                                    Arso Civil et al. (2019), Speakman et al. (2010)",
                                    "Punt et al. (2018); Hastings et al. (2012)",
                                    "Punt et al. (2018); Butterworth et al. (1995)",
                                    "Punt et al. (2018); De Long et al. (2017)",
                                    "Moore unpublished (2019); G. Vikingsson pers. comm.; Olafsdottir et al. (2003)",
                                    "Moore unpublished (2019)",
                                    "Moore unpublished (2019)",
                                    "Moore unpublished (2019)",
                                    "Moore unpublished (2019)")
)
x <- knitr::kable(format="html",
                  col.names = c('Type','Representant','$S_0$','$S_{1+}$','Âge à maturité','Référence'),
                  lh.sources, 
                  escape = FALSE) %>% 
  kable_styling('bordered')
x <- column_spec(x, column = 2,italic = TRUE)
x
```


<!-- ### Citations -->
<!-- Arso Civil, M., Cheney, B., Quick, N.J., Islas-Villanueva, V., Graves, J.A., Janik, V.M., et al. (2019). Variations in age- and sex-specific survival rates help explain population trend in a discrete marine mammal population. Ecology and Evolution, 9, 533–544. -->

<!-- Breiwick, J.M., Eberhardt, L.L. & Braham, H.W. (1984). Population Dynamics of Western Arctic Bowhead Whales ( Balaena mysticetus ). Canadian Journal of Fisheries and Aquatic Sciences, 41, 484–496. -->

<!-- Butterworth, D.S. & Punt, A.E. (1992). The Scientific Committee “...Agreed That The MSY Rate Would Most Likely Lie Between 1 and 4%” - But Which MSY Rate? Reports of the International Whaling Commission, 42, 583–591. -->

<!-- Butterworth, D.S., Punt, A.E., Oosthuizen, W.H., Wickens, P.A. (1995) The effects of future consumption by the Cape fur seal on catches and catch rates of the Cape hakes. 3. Modelling the dynamics of the Cape fur seal Arctocephalus pusillus pusillus. South African Journal of Marine Science, 16, 161–183. -->

<!-- DeLong, R.L., Melin, S.R., Laake, J.L., Morris, P., Orr, A.J. & Harris, J.D. (2017). Age- and sex-specific survival of California sea lions (Zalophus californianus) at San Miguel Island, California. Marine Mammal Science, 33, 1097–1125. -->

<!-- Dillingham, P.W., Moore, J.E., Fletcher, D., Cortés, E., Curtis, K.A., James, K.C., et al. (2016). Improved estimation of intrinsic growth rmax for long-lived species: integrating matrix models and allometry. Ecological Applications, 26, 322–333. -->

<!-- Moore 2019. Unpublished estimate following the methods of Dillingham et al. 2016. -->

<!-- Ólafsdóttir, D., Víkingsson, G. A., Halldórsson, S. D., & Sigurjónsson, J. (2003). Growth and reproduction in harbour porpoises (Phocoena phocoena) in Icelandic waters. NAMMCO Scientific Publications, 5, 195–210. -->

<!-- Punt, A. E. (1999). Annex R: A full description of the standard Baleen II model and some variants thereof. Division of Marine Research, CSIRO Marine Laboratories, Hobart, Australia. -->

<!-- Punt, A.E., Moreno, P., Brandon, J.R. & Mathews, M.A. (2018). Conserving and recovering vulnerable marine species: a comprehensive evaluation of the US approach for marine mammals. ICES Journal of Marine Science, 75, 1813–1831. -->
---
title: "Bycatch impacts exploration tool: report"
output: html_document
params:
  lh.params: NA
  performance: NA
  pbr: NA
---

This report was generated using the [marine mammal bycatch impacts exploration tool](https://msiple.shinyapps.io/mammaltool/). The output reflects parameters that the user has entered for their marine mammal population of interest. These values are entered in the  **Advanced** and **PBR & PBR calculator** tabs. 

### Life history parameters used

```{r echo=FALSE,  warning=FALSE, message=FALSE}
# The `params` object is available in the document.
options(knitr.table.format = "html") 
library(magrittr)

lhdf <- as.data.frame(params$lh.params)

lhdf <- lhdf %>%
  dplyr::select(S0,S1plus,AgeMat,lambdaMax,z) %>%
  dplyr::rename("$S_0$" = S0,
         "$S_{1+}$" = S1plus,
         "Age at maturity" = AgeMat,
         "$\\lambda_{max}$" = lambdaMax)

knitr::kable(format="html", # make table showing all the parameters
             lhdf) %>% 
            kableExtra::kable_styling()
```

### Performance measures
These are the median outcomes for each of the bycatch and depletion levels. 

```{r echo=FALSE,  warning=FALSE, message=FALSE}
# The `params` object is available in the document.
pm <-  params$performance 

knitr::kable(format="html",pm) %>%
  kableExtra::kable_styling(full_width = F,position="center")

```

### Estimated reference points

The following table is based on values that the user has entered in the **PBR** tab of the app. 

If $N_{MIN}$ and $N_{BEST}$ are the same, the user has not specified a survey CV of abundance. If several values are missing, the user may have run some projections without entering anything in the **PBR** tab. 

**Note:** For this table to give the correct values, the user has to specify a current population size and an annual bycatch rate in the **Advanced** tab. Otherwise, the numbers below will not be based on a real population size. 


```{r echo=FALSE,  warning=FALSE, message=FALSE}
# The `params` object is available in the document.
PBR <- params$pbr %>%
  dplyr::mutate(Parameter = dplyr::recode(Parameter, 
                            Nbest = "$N_{BEST}$",
                            Nmin = "$N_{MIN}$",
                            Rmax = "$R_{MAX}$",
                            Fr = "$F_R$"))
PBR %>% 
  dplyr::mutate_if(is.numeric, list(~as.character(., signif(., 2)))) %>%
  knitr::kable(format="html") %>%
  kableExtra::kable_styling(full_width = F,position="center")

```---
title: "Description du modèle"
output:
  rmarkdown::html_vignette:
  toc: true
bibliography: refs.bib
vignette: >
  %\VignetteIndexEntry{Model description (French)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(kableExtra)
```

### Modèle de population

Les projections de la population sont générées à l'aide d'un modèle unisexe structuré en âge. Dans ce modèle de projection, le nombre de baleineaux ou de nouveau-nés chaque année dépend de la densité, l'intensité de la densité dépendance étant fonction du nombre d'adultes matures $\tilde N$, de la fécondité (taux de gestation) à l'équilibre de pré-exploitation $f_0$, du taux maximum de fécondité théorique  $f_{max}$, du degré de compensation $z$, et de l'abondance des individus d'âges 1+ ($N_{y+1}^{1+}$) par rapport à la capacité de charge $K^{1+}$.

Le nombre d'individus d'âge 1 est fonction de la survie des baleineaux/nouveaux nés $S_0$, et du taux de survie des adultes, $S_{1+}$ et des prélèvements dus à la mortalité des prises accessoires $C_y$ (tous les symboles sont définis dans le Tableau 1).

$$
\tag{Eq 1}
N_{y+1,a} = 
\begin{cases}
\tilde N_{y+1}\bigg(f_0 + (f_{max}-f_0)\bigg[1-\bigg(\frac{N_{y+1}^{1+}}{K^{1+}}\bigg)^z\bigg]\bigg) &   a=0\\ 
N_{y,0}S_0 &    a=1 \\
(N_{y,a-1}-C_{y,a-1})S_{1+}  & 2\leq a<x \\
(N_{y,x-1}-C_{y,x-1})S_{1+} + (N_{y,x}-C_{y,x})S_{1+} & a=x \\
\end{cases}
$$

où $N_{y,a}$ est le nombre d'animaux  d'âge $a$ au début de l'année $y$. Le nombre d'animaux d'âge $a$ mourant en raison de la mortalité due aux prises accessoires au cours de l'année $y$, $C_{y,a}$, est retiré uniformément de la catégorie d'âge 1+ de la population (les veaux et les petits sont supposés ne pas mourir en raison des prises accessoires), c'est-à-dire:

$$
\tag{Eq 2}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
\frac{C_y N_{y,a}}{N_y^{1+}} & a>0 \\
\end{cases}
$$

Compte tenu de cette hypothèse, notre modèle ne caractérisera pas correctement les scénarios où la mortalité due aux prises accessoires se produit principalement chez les baleineaux/nouveaux-nés.

Si l'utilisateur spécifie un taux constant de mortalité par prise accessoire, les individus sont retirés de la population en fonction d'un taux de mortalité par prise accessoire $E_y$ et de la vulnérabilité (qui est supposée constante dans le temps, et uniforme sur les individus d'âge 1+) :

$$
\tag{Eq 3}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
N_{y,a} E_{y} & a>0 \\
\end{cases}
$$

**Tableau 1.** Symboles inclus dans le modèle de projection.

```{r echo = F}
Symbole=c(#User-specified parameters
  "$S_0$",
  "$S_{1+}$",
  "$x$",
  "$\\lambda_{max}$",
  "$a_p$",
  "$E_y$",
  # Parameters derived from user-specified parameters
  "$f_0$",
  "$f_{max}$",
  "$K^{1+}$",
  "$z$",
  # Derived variables
  "$\\tilde N_{y,a}$", # need double dash for using in kableExtra
  "$N_{y,a}$",
  "$N_y^{1+}$",
  "$C_{y,a}$")
Description=c(#User-specified parameters
  
  "Survie du nouveau-né ou du baleineau",
  "Survie des individus âgés de 1 an et plus",
  "Âge du groupe-plus",
  "Taux maximum d'accroissement constant (taux de croissance de la population)",
  "L'âge à la première mise bas",
  "Taux de mortalité par prises accessoires durant l'année $y$ (spécifié par l'utilisateur ou calculé à partir de la mortalité totale des prises accessoires spécifiée par l'utilisateur)", 
  # Parameters derived from user-specified parameters
  "Fécondité (taux de grossesse) à l'équilibre pré-exploitation",
  "Fécondité maximale théorique (taux de grossesse)",
  "Capacité de charge exprimée en termes de composante 1+ de la population",
  "Degré de compensation",
  # Derived variables
  "Nombre d'animaux matures d'âge $a$ au début de l'année $y$",
  "Nombre d'animaux d'âge $a$ au début de l'année $y$",
  "Nombre d'animaux d'âge 1 et plus au début de l'année $y$",
  "Mortalité due aux prises accessoires d'animaux d’âges $a$ durant l'année $y$")
Symbol_definitions <- data.frame(cbind(Symbole,
                                       Description))
knitr::kable(format="html", #html  - change to pandoc for ms word 
             Symbol_definitions, 
             escape = FALSE) %>% 
  kable_styling('bordered') %>%
  pack_rows("Paramètres spécifiés par l'utilisateur", 1, 6, label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
  pack_rows("Paramètres dérivés des paramètres spécifiés par l'utilisateur", 7, 10, label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
  pack_rows("Variables dérivées", 11, 14, label_row_css = "font-style: italic;border-bottom: 1px solid black;")
```

### Mesures de performance
Trois principales mesures de performance sont incluses dans les résultats :

- la probabilité que l'abondance des animaux d'âge 1+ dépasse le Niveau de Productivité Net Maximal (NPNM) après 50 et 100 ans ;
- l'abondance des animaux d'âge 1+ par rapport à la capacité de charge des animaux d'âge 1+ après 10, 20 et 50 ans; et
-	l'abondance des animaux d'âge 1+ par rapport à l'abondance des animaux d'âge 1+ dans un scénario de non prise accidentelle après 10, 20 et 50 ans.

Ces mesures de performance sont toutes associées à la reconstitution de la population. Wade (1998) a établi une "norme" de performance d'une probabilité de 95 % de rétablissement de la NPNM après 100 ans pour une population se trouvant initialement à 30 % de sa capacité de charge, avec une NPNM de 0,5 K. La NPNM est la limite inférieure de la population optimale durable (POD), qui est définie comme le "nombre d'animaux qui permettra d'obtenir la productivité maximale de la population ou de l'espèce, en tenant compte de la capacité de charge de l'habitat et de la santé de l'écosystème dont ils sont un élément constitutif" (16 USCS § 1362 (9)). Dans le système de gestion américain, les stocks de mammifères marins sont considérés comme épuisés lorsqu'ils sont inférieurs à la POD (16 USCS § 1362 (1B)).

Les valeurs uniques indiquées dans les tableaux sont les médianes du nombre de projections.


### Résoudre le problème du taux de prises accessoires

Nous avons résolu le problème du niveau de prises accessoires qui aboutirait à la reconstitution d'une proportion donnée de K dans un intervalle de temps donné. La probabilité de reconstitution est calculée comme étant la proportion de simulations dans lesquelles la population a atteint l'objectif de reconstitution au cours de l'année de reconstitution. L'onglet "Résoudre le problème des prises accessoires" utilise l'algorithme de recherche d'un zéro d'une fonction (la fonction `uniroot()` dans notre code) pour trouver le taux de mortalité des prises accessoires qui minimise la différence entre la probabilité de reconstitution souhaitée et la probabilité de reconstitution avec un taux de mortalité des prises accessoires de $E$.


### Fécondité théorique maximale
Le taux de fécondité théorique maximal est calculé sur la base de la taille de la population à son taux maximal de croissance démographique $\lambda_{max}$. Des parties de cette dérivation mathematique se trouvent dans Butterworth et Punt (1992) et Breiwick et al. (1984). Commençons par le modèle de population pour la composante féminine de la population :
$$
\tag{Eq 4}
N_{t+1} = N_t(1-E)S_{1+}+N_{t-a_p}f q_f S_0 (S_{1+})^{a_p -2}(1-E)^{a_p-1-a_r}
$$

où $a_p$ est l'âge à la première mise bas (supposée être un an après l'âge à maturité), $f$ est la fécondité (parfois exprimée comme le produit du taux de grossesse et de la survie des nouveau-nés), $q_f$ est la proportion de baleineaux/nouveaux-nés qui sont des femelles et $a_r$ est l'âge auquel les mammifères subissent pour la première fois une mortalité due à des prises accessoires (dans notre cas, cela est égal à 1 an). Lorsque la population croît à son rythme le plus rapide, la pêche est nulle et $N_{t+1}=N_t \lambda_{max}$. Dans notre cas, $q_f = 1$ parce que nous modélisons tous les adultes, et non uniquement femelles matures. Ainsi, l'équation ci-dessus devient :
$$
\tag{Eq 5}
\lambda_{max}N_t=N_t S_{1+}+f_{max} N_{t-a_p-1} \lambda _{max}^{a_p-1} S_0 S_{1+}^{(a_p-2)}
$$
Le calcul de $f_{max}$ permet d’obtenir la fécondité maximale théorique comme une fonction de $\lambda_{max}$, de la survie et de l'âge à la première mise bas :

$$
\tag{Eq. 6}
f_{max} = \frac {\lambda_{max}^{a_p -1} - \lambda_{max}^{{a_p -2}}}  {S_0 {S_{1+}^{a_p-2}}}
$$
Ceci est désigné sous le nom de $p$ dans Butterworth et Punt (1992).

### Niveau de productivité net maximal
Pour calculer le niveau de productivité net maximale ($NPNM$) sachany $z$, nous calculons d'abord le rendement durable par recrue $C$ en fonction du taux de mortalité des prises accessoires $E$.

$$
\tag{Eq. 7}
C = E \tilde B(E)\tilde P(E)
$$
où $\tilde B (E)$ est le recrutement standardisé lorsque le taux de prises accessoires est fixé à $E$ et $\tilde P (E)$ est le nombre d'animaux "recrutés" (âge 1+) par baleineaux/nouveaux nés lorsque le taux de mortalité des prises accessoires est fixé à $E$. Le recrutement standardisé $\tilde B(E)$ est calculé de la manière suivante :

$$
\tag{Eq. 8}
\tilde B(E) = \bigg(1 - \frac{1-f_0 \tilde N(E)} {Af_0\tilde N (E)}\bigg)^{1/z}  \bigg(\frac{\tilde B(0)\tilde P(0)}{\tilde P(E)}\bigg)
$$
où $f_0 = \frac{1}{\tilde N(0)}$, $\tilde N(E)$ est le nombre d'animaux à l'âge de la première mise bas et plus âgés (c'est-à-dire les animaux reproducteurs) par recrue à l'équilibre de la mortalité sans prise accessoire, et $A$ est le paramètre de résilience de Pella-Tomlinson ($A=\frac{f_{max}-f_0}{f_0}$ ; Punt 1999). $\tilde B(0)$ est considéré comme égal à 1, car tous les calculs sont effectués par recrue.

Le nombre d'animaux reproducteurs par recrue, pour un taux d'exploitation $E$, est égal à la somme des animaux adultes par recrue $\tilde N_a$ de l'âge à la première mise bas $a_p$ jusqu’ à l'âge du groupe plus $x$ :

$$
\tag{Eq. 9}
\tilde N(E) = \sum_{a=a_p}^{x} E \tilde N_{0,a}(E)
$$

Nous déterminons le taux de mortalité des prises accessoires auquel $\tilde C$ est maximal (c'est-à-dire où $\frac{dC}{dE}$ est égal à zéro). C'est le RMDT, qui est analogue à $F_{MSY}$ dans les pêcheries.

Le nombre d'animaux 1+ par recrue, pour un taux de mortalité des prises accessoires $\boldsymbol{E}$, $\tilde P(E)$ est défini de la façon suivante:
$$
\tag{Eq. 10}
\tilde P(E)=\sum_{a=1}^{x} \tilde N_{0,a}(E)
$$

où $\tilde N_{0,a}(E)$ est le nombre par recrue à chaque âge $a$ pour une structure d'âge stable donnée:

$$
\tag{Eq. 11}
\tilde N_{0,a}(E) = 
\begin{cases}
1 &   a=0 \\
S_0[S_{1+}(1-E)]^{a-1} &    1\leq a<x \\
\frac{S_0[S_{1+}(1-E)]^{x-1}}{1-[S_{1+}(1-E)]} &    a=x \\
\end{cases}
$$ 

$MSYR$ est la valeur de $E$ pour laquelle la dérivée de $C$ par rapport à $E$ est égale à zéro, déterminée par différenciation numérique :

$$
\tag{Eq. 12}
\frac{dC}{dE} = \frac {C(E+0.001) - C(E-0.001)} {0.002}
$$

Ensuite, l'abondance relative qui correspond au $MSYR$, $MNPL$, est déterminée en calculant la taille totale de la population 1+ au $RMDT$ par rapport à la taille de la population 1+ à l'équilibre, sans mortalité due aux prises accessoires:

$$
\tag{Eq. 13}
MNPL = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)\tilde B(0)} = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)} 
$$

### Paramétrisation
Nous considérons que la population commence avec une structure d'âge stable la première année de la période de projection (équation 11). Les valeurs à l'âge au début des projections correspondent à un taux de mortalité par prise accessoire constant $E$, qui est calculé en résolvant l'équation suivante pour $E$ :

$$
\tag{Eq. 14}
\frac{\tilde B(E)\tilde P(E)}{\tilde P(0)}= \frac{N_0^{1+}}{K^{1+}}
$$

L'épuisement initial $\frac{N_0^{1+}}{K^{1+}}$ est basé sur l'historique de la mortalité causée par l'homme au sein de la population, qui est fourni par l'utilisateur.

Dans les cas où l'erreur d'observation est donnée pour la population, l'abondance de départ est déterminée à partir d'une distribution lognormale avec un écart-type proportionnel au CV d'observation.

#### Types d'histoire de vie
Chaque type d'histoire de vie de mammifère marin présenté comme option par défaut dans cette application correspond à une valeur unique de la survie du baleineau/nouveau né $S_0$, de la survie de l'adulte $S_{1+}$, de l'âge à la première mise bas, $a_p$ et du taux intrinsèque de croissance de la population $\lambda_{max}$. Ces valeurs sont présentées dans le tableau 2. Pour des raisons de calcul, nous avons supposé que l'âge du groupe plus $x$ est supérieur de deux ans à l'âge à la maturité ($x=a_p+1$). 

#### Compensation

Nous déterminons la valeur du degré de compensation $z$ qui correspond à la valeur du NPNM fourni par l'utilisateur. Cela implique de résoudre l'équation $\tilde P(E^*) \tilde B(E^*) = MSYL * \tilde P(0)$ pour $z$ où $E^*$ dépend de $z$ comme indiqué ci-dessus.



**Tableau 2.**
```{r echo = F}
lh.sources <- data.frame(stringsAsFactors=FALSE,
                           Type = c("Bowhead whale", "Bottlenose dolphin", "Humpback whale",
                                    "Phocid seal", "Fur seal", "Sea lion",
                                    "Porpoise", "Minke whale",
                                    "False killer whale/killer whale", "Pilot whale", "Right whale"),
                 Representative = c("Balaena mysticetus", "Tursiops truncatus",
                                    "Megaptera novaeangliae", "Phoca vitulina",
                                    "Arctocephalus pusillus pusillus",
                                    "Zalophus californianus", "Phocoena phocoena",
                                    "Balaenoptera bonaerensis", "Orcinus orca",
                                    "Globicephala macrorhynchus", "Eubalaena glacialis"),
                             S0 = c(0.944, 0.865, 0.9, 0.802, 0.77, 0.83, 0.8096, 0.84216,
                                    0.84744, 0.85008, 0.85536),
                         S1plus = c(0.99, 0.951, 0.95, 0.92, 0.88, 0.95, 0.92, 0.957, 0.963,
                                    0.966, 0.972),
                         AgeMat = c(17, 6, 10, 6, 3, 4, 3, 7, 9, 9, 8),
                   #PlusGroupAge = c(25, 10, 15, 8, 10, 5, 7, 9, 11, 11, 10),
                         Source = c("@punt_conserving_2018 and references therein",
                                    "@punt_conserving_2018 and references therein,
                                    except juvenile survival which is set to 0.88($S_{1+}$)", "@punt_conserving_2018,
                                    @arso_civil_variations_2019, @speakman_mark-recapture_2010",
                                    "@punt_conserving_2018; @hastings_sex-_2012",
                                    "@punt_conserving_2018; @butterworth_effects_1995",
                                    "@punt_conserving_2018; @delong_age-_2017",
                                    "@moore_unpublished_2019; G. Vikingsson pers. comm.; olafsdottir_growth_2003",
                                    "@moore_unpublished_2019",
                                    "@moore_unpublished_2019",
                                    "@moore_unpublished_2019",
                                    "@moore_unpublished_2019")
              )


x <- knitr::kable(format="html",
             col.names = c('Type','Representant','$S_0$','$S_{1+}$','Âge à maturité','Référence'),
             lh.sources, 
             escape = FALSE) %>% 
            kable_styling('bordered')
x <- column_spec(x, column = 2,italic = TRUE)
x
```


<!-- ### Citations -->
<!-- Arso Civil, M., Cheney, B., Quick, N.J., Islas-Villanueva, V., Graves, J.A., Janik, V.M., et al. (2019). Variations in age- and sex-specific survival rates help explain population trend in a discrete marine mammal population. Ecology and Evolution, 9, 533–544. -->

<!-- Breiwick, J.M., Eberhardt, L.L. & Braham, H.W. (1984). Population Dynamics of Western Arctic Bowhead Whales ( Balaena mysticetus ). Canadian Journal of Fisheries and Aquatic Sciences, 41, 484–496. -->

<!-- Butterworth, D.S. & Punt, A.E. (1992). The Scientific Committee “...Agreed That The MSY Rate Would Most Likely Lie Between 1 and 4%” - But Which MSY Rate? Reports of the International Whaling Commission, 42, 583–591. -->

<!-- Butterworth, D.S., Punt, A.E., Oosthuizen, W.H., Wickens, P.A. (1995) The effects of future consumption by the Cape fur seal on catches and catch rates of the Cape hakes. 3. Modelling the dynamics of the Cape fur seal Arctocephalus pusillus pusillus. South African Journal of Marine Science, 16, 161–183. -->

<!-- DeLong, R.L., Melin, S.R., Laake, J.L., Morris, P., Orr, A.J. & Harris, J.D. (2017). Age- and sex-specific survival of California sea lions (Zalophus californianus) at San Miguel Island, California. Marine Mammal Science, 33, 1097–1125. -->

<!-- Dillingham, P.W., Moore, J.E., Fletcher, D., Cortés, E., Curtis, K.A., James, K.C., et al. (2016). Improved estimation of intrinsic growth rmax for long-lived species: integrating matrix models and allometry. Ecological Applications, 26, 322–333. -->

<!-- Moore 2019. Unpublished estimate following the methods of Dillingham et al. 2016. -->

<!-- Ólafsdóttir, D., Víkingsson, G. A., Halldórsson, S. D., & Sigurjónsson, J. (2003). Growth and reproduction in harbour porpoises (Phocoena phocoena) in Icelandic waters. NAMMCO Scientific Publications, 5, 195–210. -->

<!-- Punt, A. E. (1999). Annex R: A full description of the standard Baleen II model and some variants thereof. Division of Marine Research, CSIRO Marine Laboratories, Hobart, Australia. -->

<!-- Punt, A.E., Moreno, P., Brandon, J.R. & Mathews, M.A. (2018). Conserving and recovering vulnerable marine species: a comprehensive evaluation of the US approach for marine mammals. ICES Journal of Marine Science, 75, 1813–1831. -->---
title: "Description del Modelo"
output:
  rmarkdown::html_vignette:
  toc: true
bibliography: refs.bib
vignette: >
  %\VignetteIndexEntry{Model description (Spanish)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(kableExtra)
```

### Modelo de población

Las proyecciones poblacionales se generan utilizando un modelo estructurado por edad de un solo sexo. En este modelo de proyección, el número de crías nacidas cada año es denso dependiente, con la extensión de denso dependencia como función del número de adultos maduros $\tilde N$, la fecundidad (tasa de embarazo) en equilibrio previo a la explotación $f_0$, la tasa de fecundidad teórica máxima $f_{max}$, el grado de compensación $z$, y  la abundancia de individuos de edades 1+ ($N_{y+1}^{1+}$) en relación a la capacidad de carga $K^{1+}$.

El número de individuos de 1 año de edad depende de la supervivencia de crías $S_0$, de la tasa de supervivencia de los adultos, $S_{1+}$ y de remociones debido a la mortalidad por captura incidental $C_y$ (todos los símbolos están definidos en la Tabla 1).


$$
\tag{Eq 1}
N_{y+1,a} = 
\begin{cases}
\tilde N_{y+1}\bigg(f_0 + (f_{max}-f_0)\bigg[1-\bigg(\frac{N_{y+1}^{1+}}{K^{1+}}\bigg)^z\bigg]\bigg) &   a=0\\ 
N_{y,0}S_0 &    a=1 \\
(N_{y,a-1}-C_{y,a-1})S_{1+}  & 2\leq a<x \\
(N_{y,x-1}-C_{y,x-1})S_{1+} + (N_{y,x}-C_{y,x})S_{1+} & a=x \\
\end{cases}
$$

dónde $N_{y,a}$ es  el número de animales de edad $a$ a principios del año $y$. El número de animales de edad $a$ muriendo  producto de las capturas incidentales durante el año $y$, $C_{y,a}$, es removida uniformemente desde el componente 1+ de la población (se supone que las crías no mueren debido a la captura incidental), es decir:

$$
\tag{Eq 2}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
\frac{C_y N_{y,a}}{N_y^{1+}} & a>0 \\
\end{cases}
$$

Dado este supuesto, nuestro modelo no caracterizará adecuadamente los escenarios si la mortalidad por captura incidental ocurre predominantemente en las crías.

Si el usuario especifica una tasa constante de mortalidad por captura incidental, los individuos se eliminan de la población de acuerdo a una tasa de mortalidad por captura incidental $E_y$ y vulnerabilidad (que se supone constante a través del tiempo, y uniforme en individuos de edad 1+):

$$
\tag{Eq 3}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
N_{y,a} E_{y} & a>0 \\
\end{cases}
$$

**Tabla 1.**  Símbolos incluidos en el modelo de proyección.

```{r echo = F}
Símbolo=c(#User-specified parameters
  "$S_0$",
  "$S_{1+}$",
  "$x$",
  "$\\lambda_{max}$",
  "$a_p$",
  "$E_y$",
  # Parameters derived from user-specified parameters
  "$f_0$",
  "$f_{max}$",
  "$K^{1+}$",
  "$z$",
  # Derived variables
  "$\\tilde N_{y,a}$", # need double dash for using in kableExtra
  "$N_{y,a}$",
  "$N_y^{1+}$",
  "$C_{y,a}$")

Descripción=c(#User-specified parameters
  
  "Supervivencia de cachorro o ballenato",
  "Supervivencia de individuos de edad 1+",
  "Edad del grupo plus",
  "Máxima tasa constante de aumento (tasa de crecimiento de la población)",
  "Edad al primer parto",
  "Tasa de mortalidad por captura incidental en el año $y$ (especificado por el usuario o calculado a partir de la mortalidad por captura incidental total especificada por el usuario)", 
  # Parameters derived from user-specified parameters
  "Fecundidad (tasa de preñez) en equilibrio previo a la explotación",
  "Fecundidad teórica máxima (tasa de preñez)",
  "Capacidad de carga en términos de 1+ componente de la población",
  "Grado de compensación",
  # Derived variables
  "Número de animales maduros de edad $a$ a principios del año $y$",
  "Número de animales de edad $a$ a principios del año $y$",
  "Número de animales de edad 1 y mayores a principios del año $y$",
  "Mortalidad por captura incidental de animales de edad $a$ en año $y$")

Symbol_definitions <- data.frame(cbind(Símbolo,
                                       Descripción))
knitr::kable(format="html", #html  - change to pandoc for ms word 
             Symbol_definitions, 
             escape = FALSE) %>% 
  kable_styling('bordered') %>%
  pack_rows('Parámetros especificados por el usuario', 1, 6,
            label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
  pack_rows('Parámetros derivados de los parámetros especificados por el usuario', 7, 10,
            label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
  pack_rows('Parámetros derivados', 11, 14, 
            label_row_css = "font-style: italic;border-bottom: 1px solid black;")
```

### Medidas de desempeño
Hay tres medidas principales de rendimiento incluidas en los resultados:

- la probabilidad de que la abundancia de animales de edad 1+ exceda el Nivel de Productividad Neta Máxima (MNPL por sus siglas en inglés) después de 50 y 100 años;
- la abundancia de animales de edad 1+ en comparación con la capacidad de carga 1+ después de 10, 20 y 50 años; y
- la abundancia de animales de edad 1+ en comparación con la abundancia de animales de edad 1+ en un escenario sin captura incidental después de 10, 20 y 50 años.

Estas medidas de rendimiento están relacionadas con la recuperación de la población. Wade (1998) identifica un 'estándar' de rendimiento de un 95% de probabilidad de recuperación a MNPL después de 100 años para una población inicialmente al 30% de su capacidad de carga, con un MNPL de 0.5K. MNPL es el límite inferior para la Población Óptima Sostenible (OSP de su sigla en inglés), que se define como el "número de animales que dará como resultado la productividad máxima de la población o la especie, teniendo en cuenta la capacidad de carga del hábitat y la salud del ecosistema del cual forman un elemento constituyente ”(16 USCS § 1362 (9)). En el esquema de manejo de EE. UU., las poblaciones de mamíferos marinos se consideran agotadas cuando están por debajo de OSP (16 USCS § 1362 (1B)).

Los valores individuales informados en las tablas son las medianas del número de proyecciones.


### Resolviendo la tasa de captura incidental

Resolvimos el nivel de captura incidental que resultaría en la recuperación de una proporción dada de K en un período específico de tiempo. La probabilidad de reconstrucción se calcula como la proporción de simulaciones en las que la población se ha recuperado a la meta de reconstrucción para el año de reconstrucción. La pestaña "Resolver para captura incidental" utiliza la búsqueda de raíz (la función `uniroot()` en nuestro código) para encontrar la tasa de mortalidad de captura incidental que minimiza la diferencia entre la probabilidad de recuperación deseada y la probabilidad de recuperación bajo la tasa de mortalidad de captura incidental $E$.


### Fecundidad teórica máxima

La tasa de fecundidad teórica máxima se calcula en función del tamaño de la población cuando la población está creciendo a su tasa máxima de crecimiento de la población $\lambda_{max}$. Se pueden encontrar partes de esta derivación en Butterworth y Punt (1992) y Breiwick et al. (1984) Comenzando con el modelo de población para el componente femenino de la población:

$$
\tag{Eq 4}
N_{t+1} = N_t(1-E)S_{1+}+N_{t-a_p}f q_f S_0 (S_{1+})^{a_p -2}(1-E)^{a_p-1-a_r}
$$

dónde $a_p$ es la edad en el primer parto (se supone que es un año después de la edad de madurez sexual), $f$ es la fecundidad (a veces expresada como el producto de la tasa de preñez y la supervivencia de las crías), $q_f$ es la proporción de crias que son hembras y $a_r$ es la edad a la que los animales sufren por primera vez alguna mortalidad por captura incidental (en nuestro caso, esto es igual a 1 año). Cuando la población crece a su ritmo más rápido, la pesca es cero, y $N_{t+1}=N_t \lambda_{max}$. En nuestro caso $q_f = 1$ porque estamos modelando a todos los adultos, en lugar de solo hembras maduras. Por lo tanto, la ecuación anterior se convierte en:

$$
\tag{Eq 5}
\lambda_{max}N_t=N_t S_{1+}+f_{max} N_{t-a_p-1} \lambda _{max}^{a_p-1} S_0 S_{1+}^{(a_p-2)}
$$
Resolviendo para $f_{max}$ dada la máxima fecundidad teórica en función de $\lambda_{max}$, supervivencia y edad en el primer parto:

$$
\tag{Eq. 6}
f_{max} = \frac {\lambda_{max}^{a_p -1} - \lambda_{max}^{{a_p -2}}}  {S_0 {S_{1+}^{a_p-2}}}
$$
Esto se conoce como $p$ en Butterworth y Punt (1992).


### Máximo nivel de productividad neta

Para calcular el nivel máximo de productividad neta ($MNPL$) dato $z$, primero calculamos el rendimiento sostenible por recluta $C$ en función de la tasa de mortalidad por captura incidental $E$.

$$
\tag{Eq. 7}
C = E \tilde B(E)\tilde P(E)
$$
dónde $\tilde B (E)$ es el reclutamiento normalizado cuando la tasa de captura incidental se fija en $E$ y $\tilde P (E)$ es el número de equilibrio de animales "reclutados" (edad 1+) por cría cuando la tasa de mortalidad por captura incidental se fija en $E$. El reclutamiento normalizado $\tilde B(E)$ se calcula de la siguiente manera:

$$
\tag{Eq. 8}
\tilde B(E) = \bigg(1 - \frac{1-f_0 \tilde N(E)} {Af_0\tilde N (E)}\bigg)^{1/z}  \bigg(\frac{\tilde B(0)\tilde P(0)}{\tilde P(E)}\bigg)
$$
dondé $f_0 = \frac{1}{\tilde N(0)}$,  $\tilde N(E)$ es el número de animales a la edad del primer parto y mayores (es decir, animales reproductores) por recluta al equilibrio sin mortalidad incidental, y $A$ es el parámetro de resiliencia de Pella-Tomlinson (($A=\frac{f_{max}-f_0}{f_0}$ ; Punt 1999). $\tilde B(0)$ asumiendo que es igual a 1, porque todos los cálculos son por recluta.

El número de animales reproductores por recluta a la tasa de explotación $\boldsymbol{E}$ es la suma de animales adultos por recluta $\tilde N_a$ desde la edad del primer parto $a_p$ a la edad del grupo plus $x$:

$$
\tag{Eq. 9}
\tilde N(E) = \sum_{a=a_p}^{x} \tilde N_{a}(E)
$$

Resolvemos la tasa de mortalidad por captura incidental a la que $\tilde C$ es maximizada (es decir, donde $\frac{dC}{dE}$ es cero) Esto es $MSYR$, que es análogo a $F_{MSY}$ en pesquerías.

El número de 1+ animales por recluta a la tasa de mortalidad por captura incidental $E$, $\tilde P(E)$ se define como:
$$
\tag{Eq. 10}
\tilde P(E)=\sum_{a=1}^{x} \tilde N_{a}(E)
$$

dónde $\tilde N_{0,a}(E)$ son los números por recluta en cada edad $a$  dada una estructura de edad estable:

$$
\tag{Eq. 11}
\tilde N_{0,a}(E) = 
\begin{cases}
1 &   a=0 \\
S_0[S_{1+}(1-E)]^{a-1} &    1\leq a<x \\
\frac{S_0[S_{1+}(1-E)]^{x-1}}{1-[S_{1+}(1-E)]} &    a=x \\
\end{cases}
$$ 

$MSYR$ es el valor de $E$ para el cual la derivada de $C$ con respecto a $E$ es cero, que determinamos mediante diferenciación numérica:

$$
\tag{Eq. 12}
\frac{dC}{dE} = \frac {C(E+0.001) - C(E-0.001)} {0.002}
$$

Entonces, la abundancia relativa que corresponde a $MSYR$, $MNPL$, se determina calculando el tamaño total de la población 1+ en $MSYR$ en relación con el tamaño de la población de equilibrio 1+ sin mortalidad por captura incidental:

$$
\tag{Eq. 13}
MNPL = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)\tilde B(0)} = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)} 
$$

### Parametrización
Supongamos que la población comienza con una estructura de edad estable en el año 1 del período de proyección (Ec. 11). Los números a la edad al comienzo de las proyecciones corresponden a una tasa constante de mortalidad por captura incidental $E$, que se calcula resolviendo la siguiente ecuación para $E$: 

$$
\tag{Eq. 14}
\frac{\tilde B(E)\tilde P(E)}{\tilde P(0)}= \frac{N_0^{1+}}{K^{1+}}
$$

El agotamiento inicial $\frac{N_0^{1+}}{K^{1+}}$  se basa en la historia de mortalidad causada por humanos para la población, que es proporcionada por el usuario.

Para los casos en que se da un error de observación para la población, la abundancia inicial se extrae de una distribución log normal con una desviación estándar proporcional al CV de observación.

#### Tipos de historia de vida
Cada tipo de historia de vida de mamíferos marinos presentado como una opción predeterminada en esta aplicación corresponde a un valor único de supervivencia de las crías $S_0$, sobrevivencia adulta $S_{1+}$, edad al primer parto, $a_p$, y tasa intrínseca de crecimiento de la población $\lambda_{max}$. Estos valores se presentan en la Tabla 2. Para fines de cálculo, supusimos que la edad del grupo más $x$ es dos años después de la edad de madurez ($x=a_p+1$).

#### Compensación
Resolvemos el valor del grado de compensación $z$ que corresponde al valor de MNPL proporcionado por el usuario. Esto implica resolver la ecuación $\tilde P(E^*) \tilde B(E^*) = MSYL * \tilde P(0)$ para $z$  donde $E^*$ depende de $z$  como se describe arriba.



**Tabla 2.**
```{r echo = F}
lh.sources <- data.frame(stringsAsFactors=FALSE,
                           Tipo = c("Ballena de Groenlandia o ballena boreal",
                                  "Delfín nariz de botella",
                                  "Ballena jorobada",
                                  "Foca común",
                                  "Lobo fino o de dos pelos",
                                  "Lobo común o de un pelo",
                                  "Marsopa común",
                                  "Ballena minke antártica",
                                  "Orca",
                                  "Calderón de aletas cortas",
                                  "Ballena franca"),
                         Representante = c("Balaena mysticetus",
                                           "Tursiops truncatus",
                                           "Megaptera novaeangliae",
                                           "Phoca vitulina",
                                           "Arctocephalus pusillus pusillus",
                                           "Zalophus californianus",
                                           "Phocoena phocoena",
                                           "Balaenoptera bonaerensis",
                                           "Orcinus orca",
                                           "Globicephala macrorhynchus",
                                           "Eubalaena glacialis"),
                             S0 = c(0.944, 0.865, 0.9, 0.802, 0.77, 0.83, 0.8096, 0.84216,
                                    0.84744, 0.85008, 0.85536),
                         S1plus = c(0.99, 0.951, 0.95, 0.92, 0.88, 0.95, 0.92, 0.957, 0.963,
                                    0.966, 0.972),
                         AgeMat = c(17, 6, 10, 6, 3, 4, 3, 7, 9, 9, 8),
                   #PlusGroupAge = c(25, 10, 15, 8, 10, 5, 7, 9, 11, 11, 10),
                         Source = c("@punt_conserving_2018 and references therein",
                                    "@punt_conserving_2018 and references therein,
                                    except juvenile survival which is set to 0.88($S_{1+}$)", "@punt_conserving_2018,
                                    @arso_civil_variations_2019, @speakman_mark-recapture_2010",
                                    "@punt_conserving_2018; @hastings_sex-_2012",
                                    "@punt_conserving_2018; @butterworth_effects_1995",
                                    "@punt_conserving_2018; @delong_age-_2017",
                                    "@moore_unpublished_2019; G. Vikingsson pers. comm.; olafsdottir_growth_2003",
                                    "@moore_unpublished_2019",
                                    "@moore_unpublished_2019",
                                    "@moore_unpublished_2019",
                                    "@moore_unpublished_2019")
              )


x <- knitr::kable(format="html",
             col.names = c('Tipo','Representante','$S_0$','$S_{1+}$','Edad de madurez','Referencias'),
             lh.sources, 
             escape = FALSE) %>% 
            kable_styling('bordered')
x <- column_spec(x, column = 2,italic = TRUE)
x
```


<!-- ### Citations -->
<!-- Arso Civil, M., Cheney, B., Quick, N.J., Islas-Villanueva, V., Graves, J.A., Janik, V.M., et al. (2019). Variations in age- and sex-specific survival rates help explain population trend in a discrete marine mammal population. Ecology and Evolution, 9, 533–544. -->

<!-- Breiwick, J.M., Eberhardt, L.L. & Braham, H.W. (1984). Population Dynamics of Western Arctic Bowhead Whales ( Balaena mysticetus ). Canadian Journal of Fisheries and Aquatic Sciences, 41, 484–496. -->

<!-- Butterworth, D.S. & Punt, A.E. (1992). The Scientific Committee “...Agreed That The MSY Rate Would Most Likely Lie Between 1 and 4%” - But Which MSY Rate? Reports of the International Whaling Commission, 42, 583–591. -->

<!-- Butterworth, D.S., Punt, A.E., Oosthuizen, W.H., Wickens, P.A. (1995) The effects of future consumption by the Cape fur seal on catches and catch rates of the Cape hakes. 3. Modelling the dynamics of the Cape fur seal Arctocephalus pusillus pusillus. South African Journal of Marine Science, 16, 161–183. -->

<!-- DeLong, R.L., Melin, S.R., Laake, J.L., Morris, P., Orr, A.J. & Harris, J.D. (2017). Age- and sex-specific survival of California sea lions (Zalophus californianus) at San Miguel Island, California. Marine Mammal Science, 33, 1097–1125. -->

<!-- Dillingham, P.W., Moore, J.E., Fletcher, D., Cortés, E., Curtis, K.A., James, K.C., et al. (2016). Improved estimation of intrinsic growth rmax for long-lived species: integrating matrix models and allometry. Ecological Applications, 26, 322–333. -->

<!-- Moore 2019. Unpublished estimate following the methods of Dillingham et al. 2016. -->

<!-- Ólafsdóttir, D., Víkingsson, G. A., Halldórsson, S. D., & Sigurjónsson, J. (2003). Growth and reproduction in harbour porpoises (Phocoena phocoena) in Icelandic waters. NAMMCO Scientific Publications, 5, 195–210. -->

<!-- Punt, A. E. (1999). Annex R: A full description of the standard Baleen II model and some variants thereof. Division of Marine Research, CSIRO Marine Laboratories, Hobart, Australia. -->

<!-- Punt, A.E., Moreno, P., Brandon, J.R. & Mathews, M.A. (2018). Conserving and recovering vulnerable marine species: a comprehensive evaluation of the US approach for marine mammals. ICES Journal of Marine Science, 75, 1813–1831. -->---
title: "Model description"
output:
  rmarkdown::html_vignette:
  toc: true
bibliography: refs.bib
vignette: >
  %\VignetteIndexEntry{Model description}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(kableExtra)
```

### Population model

Population projections are generated using a single-sex age-structured model. In this projection model, the number of calves or pups born each year is density dependent, with the extent of density dependence being a function of the number of mature adults $\tilde N$, the fecundity (pregnancy rate) at pre-exploitation equilibrium $f_0$, the maximum theoretical fecundity rate $f_{max}$, the degree of compensation $z$, and the abundance of individuals aged 1+ ($N_{y+1}^{1+}$) relative to carrying capacity $K^{1+}$.

The number of individuals age 1 and older is a function of calf/pup survival, $S_0$, and the survival rate of adults, $S_{1+}$, and removals due to bycatch mortality, $C_y$ (all symbols are defined in Table 1).


$$
\tag{Eq 1}
N_{y+1,a} = 
\begin{cases}
\tilde N_{y+1}\bigg(f_0 + (f_{max}-f_0)\bigg[1-\bigg(\frac{N_{y+1}^{1+}}{K^{1+}}\bigg)^z\bigg]\bigg) &   a=0\\ 
N_{y,0}S_0 &    a=1 \\
(N_{y,a-1}-C_{y,a-1})S_{1+}  & 2\leq a<x \\
(N_{y,x-1}-C_{y,x-1})S_{1+} + (N_{y,x}-C_{y,x})S_{1+} & a=x \\
\end{cases}
$$

where $N_{y,a}$ is the number of animals of age $a$ at the start of year $y$. The number of animals of age $a$ dying due to bycatch mortality during year $y$, $C_{y,a}$, is removed uniformly from the 1+ component of the population (calves and pups are assumed not to die due to bycatch), i.e.:

$$
\tag{Eq 2}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
\frac{C_y N_{y,a}}{N_y^{1+}} & a>0 \\
\end{cases}
$$

Given this assumption, our model will not adequately characterize cases where bycatch mortality occurs predominantly among calves/pups. 

If the user specifies a constant rate of bycatch mortality, individuals are removed from the population according to a bycatch mortality rate $E_y$ and vulnerability (which is assumed constant through time, and uniform on age 1+ individuals):

$$
\tag{Eq 3}
C_{y,a}=
\begin{cases}
0 &   a=0 \\ 
N_{y,a} E_{y} & a>0 \\
\end{cases}
$$
If bycatch is specified by the user as a number of individuals per year, annual bycatch rates $E_y$ are lognormally distributed with mean $log(C_{user})$ and standard deviation $\sigma_E$ of $\sqrt{log(CV_C^2 + 1)}$, where $C_{user}$ is the bycatch in numbers per year defined by the user and $CV_C$ is the bycatch CV defined by the user. If bycatch is specified as a constant rate, $E_y$ is beta-distributed with mean $E_{user}$ and standard deviation $E_{user} \cdot CV_E$.


**Table 1.** Symbols included in the projection model.
 
```{r echo = F}
Symbol=c(#User-specified parameters
         "$S_0$",
         "$S_{1+}$",
         "$x$",
         "$\\lambda_{max}$",
          "$a_p$",
         "$E_y$",
         # Parameters derived from user-specified parameters
         "$f_0$",
         "$f_{max}$",
         "$K^{1+}$",
         "$z$",
         # Derived variables
         "$\\tilde N_{y,a}$", # need double dash for using in kableExtra
         "$N_{y,a}$",
         "$N_y^{1+}$",
         "$C_{y,a}$")

Description=c(#User-specified parameters
              "Pup or calf survival",
              "Survival of individuals aged 1+",
              "Plus-group age",
              "Maximum steady rate of increase (population growth rate)",
              "Age at first parturition",
              "Bycatch mortality rate in year $y$ (specified by the user or computed from the user-specified total bycatch mortality)", 
              # Parameters derived from user-specified parameters
              "Fecundity (pregnancy rate) at pre-exploitation equilibrium",
              "Maximum theoretical fecundity (pregnancy rate)",
              "Carrying capacity in terms of 1+ component of the population",
              "Degree of compensation",
              # Derived variables
              "Number of mature animals of age $a$ at the start of year $y$",
              "Number of animals of age $a$ at the start of year $y$",
              "Number of animals aged 1 and older at the start of year $y$",
              "Mortality due to bycatch of animals of age $a$ in year $y$")

Symbol_definitions <- data.frame(cbind(Symbol,
                                      Description))
knitr::kable(format="html", #html  - change to pandoc for ms word 
             Symbol_definitions, 
             escape = FALSE) %>% 
            kable_styling('bordered') %>%
            pack_rows('User-specified parameters', 1, 6, label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
            pack_rows('Parameters derived from user-specified parameters', 7, 10, label_row_css = "font-style: italic;border-bottom: 1px solid black;") %>%
            pack_rows('Derived variables', 11, 14, label_row_css = "font-style: italic;border-bottom: 1px solid black;")
```

### Performance measures
There are three primary performance measures included in the outputs:

-	the probability that the abundance of age 1+ animals exceeds the Maximum Net Productivity Level (MNPL) after 50 and 100 years;
-	the abundance of age 1+ animals compared to 1+ carrying capacity after 10, 20, and 50 years; and
-	the abundance of age 1+ animals compared to the abundance of age 1+ animals under a no-bycatch scenario after 10, 20, and 50 years.

These performance measures are all related to population recovery. Wade (1998) identifies a performance ‘standard’ of a 95% probability of recovery to MNPL after 100 years for a population initially at 30% of its carrying capacity, with an MNPL of 0.5K. MNPL is the lower bound for the U.S. Marine Mammal Protection Act Optimum Sustainable Population (OSP), which is defined as the “number of any animals which will result in the maximum productivity of the population or the species, keeping in mind the carrying capacity of the habitat and the health of the ecosystem of which they form a constituent element” (16 USCS § 1362 (9)). In the US management scheme, marine mammal stocks are considered depleted when they are below OSP (16 USCS § 1362 (1B)).

Single values reported in tables are the medians from the number of projections.


### Solving for bycatch rate

The app can solve for the bycatch level that would result in recovery to a given proportion of K in a specific amount of time. The rebuilding probability is calculated as the proportion of simulations in which the population has recovered to the rebuilding goal by the rebuilding year. The “Solve for bycatch” tab uses root finding (the function `uniroot()` in our code) to find the bycatch mortality rate that minimizes the difference between the desired recovery probability and the recovery probability under bycatch mortality rate $E$.


### Maximum theoretical fecundity
The maximum theoretical fecundity rate is calculated based on the population size when the population is growing at its maximum population growth rate $\lambda_{max}$. Parts of this derivation can be found in @butterworth_scientific_1992 and @breiwick_population_1984. Starting with the population model for the female component of the population:

$$
\tag{Eq 4}
N_{t+1} = N_t(1-E)S_{1+}+N_{t-a_p}f q_f S_0 (S_{1+})^{a_p -2}(1-E)^{a_p-1-a_r}
$$

where $a_p$ is the age at first parturition (assumed to be one year after the age at maturity), $f$ is the fecundity (sometimes expressed as the product of pregnancy rate and pup survival), $q_f$ is the proportion of calves/pups that are female and $a_r$ is the age at which the mammals first suffer any bycatch mortality (in our case, this is equal to 1 year). When the population is growing at its fastest rate, fishing is zero, and $N_{t+1}=N_t \lambda_{max}$. In our case $q_f = 1$ because we are modeling all adults, instead of just mature females. Thus, the above equation becomes:

$$
\tag{Eq 5}
\lambda_{max}N_t=N_t S_{1+}+f_{max} N_{t-a_p-1} \lambda _{max}^{a_p-1} S_0 S_{1+}^{(a_p-2)}
$$
Solving for $f_{max}$ gives the maximum theoretical fecundity as a function of $\lambda_{max}$, survival, and age at first parturition:

$$
\tag{Eq. 6}
f_{max} = \frac {\lambda_{max}^{a_p -1} - \lambda_{max}^{{a_p -2}}}  {S_0 {S_{1+}^{a_p-2}}}
$$
This is referred to as $p$ in @butterworth_scientific_1992. 

### Maximum net productivity level
To calculate the maximum net productivity level ($MNPL$) given $z$, we first calculate the sustainable yield $C$  as a function of bycatch mortality rate $E$.

$$
\tag{Eq. 7}
C = E \tilde B(E)\tilde P(E)
$$
where $\tilde B (E)$ is the normalized (to numbers at unfished equilibrium) recruitment when the bycatch rate is fixed at $E$ and $\tilde P (E)$ is the equilibrium number of "recruited" (age 1+) animals per calf/pup when the bycatch mortality rate is fixed at $E$. The normalized recruitment $\tilde B(E)$ is calculated as follows:

$$
\tag{Eq. 8}
\tilde B(E) = \bigg(1 - \frac{1-f_0 \tilde N(E)} {Af_0\tilde N (E)}\bigg)^{1/z}  \bigg(\frac{\tilde B(0)\tilde P(0)}{\tilde P(E)}\bigg)
$$
where $f_0 = \frac{1}{\tilde N(0)}$,  $\tilde N(E)$ is the number of animals at the age of first parturition and older (i.e., reproducing animals) per recruit at no-bycatch-mortality equilibrium, and $A$ is the Pella-Tomlinson resilience parameter ($A=\frac{f_{max}-f_0}{f_0}$ ; @punt_a._e._annex_1999). $\tilde B(0)$ is assumed to be equal to 1, because all calculations are per-recruit.


The number of reproducing animals per recruit at exploitation rate $\boldsymbol{E}$ is the sum of adult animals per recruit $\tilde N_a$ from the age at first parturition $a_p$ to the plus group age $x$:

$$
\tag{Eq. 9}
\tilde N(E) = \sum_{a=a_p}^{x} \tilde N_{a}(E)
$$

We solve for the bycatch mortality rate at which $\tilde C$ is maximized (i.e., where $\frac{dC}{dE}$ is zero). This is $MSYR$, which is analogous to $F_{MSY}$ in fisheries.

The number of 1+ animals per recruit at bycatch mortality rate $E$, $\tilde P(E)$ is defined as:
$$
\tag{Eq. 10}
\tilde P(E)=\sum_{a=1}^{x} \tilde N_{a}(E)
$$

where $\tilde N_{0,a}(E)$ is the numbers per recruit at each age $a$ given a stable age structure:

$$
\tag{Eq. 11}
\tilde N_{0,a}(E) = 
\begin{cases}
1 &   a=0 \\
S_0[S_{1+}(1-E)]^{a-1} &    1\leq a<x \\
\frac{S_0[S_{1+}(1-E)]^{x-1}}{1-[S_{1+}(1-E)]} &    a=x \\
\end{cases}
$$ 

$MSYR$ is the value of $E$ for which the derivative of $C$ with respect to $E$ is zero, which we determined through numerical differentiation:

$$
\tag{Eq. 12}
\frac{dC}{dE} = \frac {C(E+0.001) - C(E-0.001)} {0.002}
$$

Then the relative abundance that corresponds to $MSYR$, $MNPL$, is determined by calculating the total 1+ population size at $MSYR$ relative to the equilibrium 1+ population size with no bycatch mortality:

$$
\tag{Eq. 13}
MNPL = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)\tilde B(0)} = \frac{\tilde P(E=MSYR)\tilde B(E=MSYR)}{\tilde P(0)} 
$$

### Parameterization
We assume that the population starts with a stable age structure in year 1 of the projection period (Eq. 11). The numbers at age at the start of projections correspond to a constant bycatch mortality rate $E$, which is calculated by solving the following equation for $E$: 

$$
\tag{Eq. 14}
\frac{\tilde B(E)\tilde P(E)}{\tilde P(0)}= \frac{N_0^{1+}}{K^{1+}}
$$

The initial depletion $\frac{N_0^{1+}}{K^{1+}}$ is based on the history of human-caused mortality for the population, which is provided by the user.

For cases where observation error is given for the initial population size, the starting abundance is drawn from a lognormal distribution with a standard deviation proportional to the observation CV.

#### Life history types
Each marine mammal life history type presented as a default option in this app corresponds to a unique value of calf/pup survival $S_0$, adult survival $S_{1+}$, age at first parturition, $a_p$, and intrinsic rate of population growth $\lambda_{max}$. These values are presented in Table 2. For computation purposes, we assumed that the plus group age $x$ is two years past the age at maturity ($x=a_p+1$). 

#### Compensation
We solve for the value of the degree of compensation $z$ that corresponds to the value of MNPL provided by the user. This involves solving the equation   $\tilde P(E^*) \tilde B(E^*) = MSYL * \tilde P(0)$ for $z$  where $E^*$ depends on $z$ as outlined above.



**Table 2.**
```{r echo = F}
lh.sources <- data.frame(stringsAsFactors=FALSE,
                           Type = c("Bowhead whale", "Bottlenose dolphin", "Humpback whale",
                                    "Phocid seal", "Fur seal", "Sea lion",
                                    "Porpoise", "Minke whale",
                                    "False killer whale/killer whale", "Pilot whale", "Right whale"),
                 Representative = c("Balaena mysticetus", "Tursiops truncatus",
                                    "Megaptera novaeangliae", "Phoca vitulina",
                                    "Arctocephalus pusillus pusillus",
                                    "Zalophus californianus", "Phocoena phocoena",
                                    "Balaenoptera bonaerensis", "Orcinus orca",
                                    "Globicephala macrorhynchus", "Eubalaena glacialis"),
                             S0 = c(0.944, 0.865, 0.9, 0.802, 0.77, 0.83, 0.8096, 0.84216,
                                    0.84744, 0.85008, 0.85536),
                         S1plus = c(0.99, 0.951, 0.95, 0.92, 0.88, 0.95, 0.92, 0.957, 0.963,
                                    0.966, 0.972),
                         AgeMat = c(17, 6, 10, 6, 3, 4, 3, 7, 9, 9, 8),
                   #PlusGroupAge = c(25, 10, 15, 8, 10, 5, 7, 9, 11, 11, 10),
                         Source = c("@punt_conserving_2018 and references therein",
                                    "@punt_conserving_2018 and references therein,
                                    except juvenile survival which is set to 0.88($S_{1+}$)", "@punt_conserving_2018,
                                    @arso_civil_variations_2019, @speakman_mark-recapture_2010",
                                    "@punt_conserving_2018; @hastings_sex-_2012",
                                    "@punt_conserving_2018; @butterworth_effects_1995",
                                    "@punt_conserving_2018; @delong_age-_2017",
                                    "@moore_unpublished_2019; G. Vikingsson pers. comm.; olafsdottir_growth_2003",
                                    "@moore_unpublished_2019",
                                    "@moore_unpublished_2019",
                                    "@moore_unpublished_2019",
                                    "@moore_unpublished_2019")
              )


x <- knitr::kable(format="html",
             col.names = c('Type','Representative','$S_0$','$S_{1+}$','Age at maturity','Source'),
             lh.sources, 
             escape = FALSE) %>% 
            kable_styling('bordered')
x <- column_spec(x, column = 2,italic = TRUE)
x
```


<!-- ### Citations -->
<!-- Arso Civil, M., Cheney, B., Quick, N.J., Islas-Villanueva, V., Graves, J.A., Janik, V.M., et al. (2019). Variations in age- and sex-specific survival rates help explain population trend in a discrete marine mammal population. Ecology and Evolution, 9, 533–544. -->

<!-- Breiwick, J.M., Eberhardt, L.L. & Braham, H.W. (1984). Population Dynamics of Western Arctic Bowhead Whales ( Balaena mysticetus ). Canadian Journal of Fisheries and Aquatic Sciences, 41, 484–496. -->

<!-- Butterworth, D.S. & Punt, A.E. (1992). The Scientific Committee “...Agreed That The MSY Rate Would Most Likely Lie Between 1 and 4%” - But Which MSY Rate? Reports of the International Whaling Commission, 42, 583–591. -->

<!-- Butterworth, D.S., Punt, A.E., Oosthuizen, W.H., Wickens, P.A. (1995) The effects of future consumption by the Cape fur seal on catches and catch rates of the Cape hakes. 3. Modelling the dynamics of the Cape fur seal Arctocephalus pusillus pusillus. South African Journal of Marine Science, 16, 161–183. -->

<!-- DeLong, R.L., Melin, S.R., Laake, J.L., Morris, P., Orr, A.J. & Harris, J.D. (2017). Age- and sex-specific survival of California sea lions (Zalophus californianus) at San Miguel Island, California. Marine Mammal Science, 33, 1097–1125. -->

<!-- Dillingham, P.W., Moore, J.E., Fletcher, D., Cortés, E., Curtis, K.A., James, K.C., et al. (2016). Improved estimation of intrinsic growth rmax for long-lived species: integrating matrix models and allometry. Ecological Applications, 26, 322–333. -->

<!-- Moore 2019. Unpublished estimate following the methods of Dillingham et al. 2016. -->

<!-- Ólafsdóttir, D., Víkingsson, G. A., Halldórsson, S. D., & Sigurjónsson, J. (2003). Growth and reproduction in harbour porpoises (Phocoena phocoena) in Icelandic waters. NAMMCO Scientific Publications, 5, 195–210. -->

<!-- Punt, A. E. (1999). Annex R: A full description of the standard Baleen II model and some variants thereof. Division of Marine Research, CSIRO Marine Laboratories, Hobart, Australia. -->

<!-- Punt, A.E., Moreno, P., Brandon, J.R. & Mathews, M.A. (2018). Conserving and recovering vulnerable marine species: a comprehensive evaluation of the US approach for marine mammals. ICES Journal of Marine Science, 75, 1813–1831. -->% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11_get_dz.R
\name{get_dz}
\alias{get_dz}
\title{Derivative of z (function to minimize)}
\usage{
get_dz(z, MNPL, lh.params)
}
\arguments{
\item{z}{the degree of compensation}

\item{MNPL}{desired maximum net productivity level}

\item{lh.params}{a list of life history parameters (juvenile survival S0, adult survival S1plus, age at maturity AgeMat, plus group age, max theoretical fecundity fmax, and maximum steady rate of increase (population growth rate) lambdaMax)}
}
\value{
The difference between the MNPL associated with the value of z that the user defined and the MNPL that the user has defined.
}
\description{
Derivative of z (function to minimize)
}
\examples{
get_dz(z = 2.39, MNPL = 0.5, 
lh.params = list(S0 = 0.944, S1plus = 0.99, 
AgeMat = 17, nages = 19, 
lambdaMax = 1.04, K1plus = 9000))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dat-data.R
\docType{data}
\name{dat}
\alias{dat}
\title{Life history parameters for pinnipeds}
\format{
A data frame with 11 rows and 11 variables:
\describe{
  \item{Species}{Species}
  \item{Parameter}{Life history parameter}
  \item{Value}{Parameter value}
  \item{Range}{Parameter range, if a range is reported in the paper instead of a single value}
  \item{Citation}{Source for parameter value or range}
  \item{doi}{doi for source}
  \item{Location}{Location where data were collected}
  \item{Time.period}{Time period studied in the paper}
  \item{Sex.specific.}{If only one sex was included in the paper}
  \item{Notes}{Messy notes about data collection and source}
  ...
}
}
\source{
Compiled from literature by M Siple and coauthors
}
\usage{
dat
}
\description{
A table with life history information for pinnipeds.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07_ce.R
\name{ce}
\alias{ce}
\title{Calculate normalized sustainable yield}
\usage{
ce(
  S0 = NA,
  S1plus = NA,
  AgeMat = NA,
  nages = NA,
  z = NA,
  E = NA,
  A = NA,
  P0 = NA,
  N0 = NA
)
}
\arguments{
\item{S0}{Calf/pup survival, a numeric value between 0 and 1}

\item{S1plus}{1+ survival rate for animals age 1 year and older, a numeric value between 0 and 1}

\item{AgeMat}{Age at maturity (= age at first parturition - 1). Must be less than \code{nages}}

\item{nages}{"maximum" age, treated as the plus group age. The plus group age can be set equal to the age at maturity +2 years without losing accuracy.}

\item{z}{degree of compensation}

\item{E}{bycatch mortality rate (applies to 1+ numbers)}

\item{A}{the Pella-Tomlinson resilience parameter ((fmax - f0)/f0)}

\item{P0}{unfished number-per-recruit - 1+ adults}

\item{N0}{unfished numbers-per-recruit - mature adults}
}
\value{
a single value of normalized yield for exploitation rate E
}
\description{
This function calculates the normalized sustainable yield, which is used to find MNPL (the population size at which productivity is maximized).
}
\examples{
# Set parameters
S0.w = 0.5; S1plus.w = 0.944; nages.w = 25; AgeMat.w = 18 
# Get number of individuals per recruit in terms of mature individuals (N0.w)
NPROut <- npr(S0 = S0.w, S1plus = S1plus.w, nages = nages.w, AgeMat = AgeMat.w, E = 0)

N0 <- NPROut$npr # mature numbers per recruit
# Get number of individuals per recruit in terms of individuals aged 1+ (P0.w)
P0 <- NPROut$P1r # 1+ nums per recruit

ce(S0 = S0.w, S1plus = S1plus.w, 
nages = nages.w, 
AgeMat = AgeMat.w, 
E=0.01, z=2.39,A=2, N0 = N0, P0 = P0)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/13_pop_vs_yield.R
\name{pop_vs_yield}
\alias{pop_vs_yield}
\title{Yield or productivity curve}
\usage{
pop_vs_yield(
  z.vec = c(1, 2.39, 5.99),
  lh.params = lh.params1,
  add.legend = FALSE,
  ggp = TRUE,
  linecolor = "#123f5a",
  lang = "en"
)
}
\arguments{
\item{z.vec}{a vector of z values}

\item{lh.params}{a list of life history parameters}

\item{add.legend}{logical; whether or not to add a legend}

\item{ggp}{logical; whether to plot in ggplot and plot a single yield curve; set to \code{FALSE} for base R plot and multiple yield curves.}

\item{linecolor}{color of yield curve line}

\item{lang}{language selected by the user (character)}
}
\value{
a plot of population size vs. 'yield'
}
\description{
relative population size (x) vs. Sustainable yield (y)
}
\examples{
pop_vs_yield(z.vec = c( 1.1, 2.5, 5.99),
lh.params = list(
  S0 = 0.944, S1plus = 0.99, AgeMat = 17,
  nages = 19, lambdaMax = 1.04
), ggp = FALSE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_inv_logit.R
\name{inv_logit}
\alias{inv_logit}
\title{Get inverse logit of a value}
\usage{
inv_logit(x)
}
\arguments{
\item{x}{A number}
}
\value{
the inverse logit of x.
}
\description{
The inverse logit of x
}
\examples{
inv_logit(0)
inv_logit(1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/23_prob_rebuilt_goal.R
\name{prob_rebuilt_goal}
\alias{prob_rebuilt_goal}
\title{Calculate the probability of reaching a rebuilding goal}
\usage{
prob_rebuilt_goal(traj, goal = 2000, rebuild.yr = 100)
}
\arguments{
\item{traj}{a matrix of trajectories, with rows=nsims and cols=nyears}

\item{goal}{rebuilding goal as an absolute number of animals}

\item{rebuild.yr}{the year by which you want pop to be recovered (calculate probability that pop will recover to MNPL by rebuild.yr)}
}
\value{
- the probability that the stock will be rebuilt to the goal population size by year \code{rebuilt.yr}
}
\description{
function that calculates the probability of recovery (in terms of 1+ numbers) to MNPL (or other population size) at \code{rebuilt.yr} years
}
\examples{
lh.params <- list(
  S0 = 0.944, S1plus = 0.99, 
  AgeMat = 17,
  nages = 19, z = 2.39, 
  lambdaMax = 1.04, K1plus = 9000
)
Projection <- projections(
  NOut = 100,
  ConstantBycatch = list(Catch = 40, CV = 01),
  InitDepl = 0.3,
  lh.params = lh.params,
  nyears = 100, obs_CV = 0.2
)
prob_rebuilt_goal(traj = Projection$trajectories, goal = 3500, rebuild.yr = 25)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/19_ggradar.R
\name{ggradar}
\alias{ggradar}
\title{Radar plot}
\usage{
ggradar(
  plot.data,
  axis.labels = colnames(plot.data)[-1],
  grid.label.size = 7,
  axis.label.size = 8,
  plot.legend = if (nrow(plot.data) > 1) TRUE else FALSE,
  legend.text.size = grid.label.size,
  palette.vec = c("#D53E4F", "#FC8D59", "#FEE08B", "#E6F598", "#99D594", "#3288BD")
)
}
\arguments{
\item{plot.data}{a dataframe with performance measures as columns and management scenarios as rows}

\item{axis.labels}{a vector of names for the performance measures}

\item{grid.label.size}{a numeric value for grid label size}

\item{axis.label.size}{a numeric value for axis label size}

\item{plot.legend}{logical; whether or not to plot a legend to the right of the plot}

\item{legend.text.size}{numeric value for size of legend text}

\item{palette.vec}{a vector of colors to use for the different scenarios (each row = 1 color)}
}
\description{
Radar plot
}
\details{
I modified the `ggradar()` function slightly to do nice things that I like, like use a custom color palette, including different line types, etc.

Since this code was originally written, ggradar has becomes its own standalone package. For more information and for the most current version of the function, see Ricardo Bion's \href{https://github.com/ricardo-bion/ggradar/}{GitHub}
}
\examples{
plot.data <- data.frame(
  bycatch = factor(c(
    "Higher end of bycatch range",
    "Midpoint of bycatch range",
    "Lower end of bycatch range"
  ),
  ,levels = c(
  "Higher end of bycatch range",
  "Midpoint of bycatch range",
  "Lower end of bycatch range"
  )),
  prebuild50 = c(0, 0, 1),
  prebuild100 = c(0, 0, 1),
  abundrel10 = c(0, 0, 0.29),
  abundrel20 = c(0, 0, 0.34),
  abundrel50 = c(0, 0, 0.53)
)

ggradar(
  plot.data = plot.data,
  axis.label.size = 4,
  grid.label.size = 4,
  palette.vec = c("#99D594","#E6F598", "#FEE08B") 
)

}
\author{
Ricardo Bion, modified slightly by Margaret Siple
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/12_calc_z.R
\name{calc_z}
\alias{calc_z}
\title{Degree of compensation}
\usage{
calc_z(MNPL_in, lh.params_in)
}
\arguments{
\item{MNPL_in}{User-specified value for MNPL as a proportion of K (between 0 and 1)}

\item{lh.params_in}{a list of life history parameters. Must contain S0, S1plus, nages, AgeMat, lambdaMax and z.}
}
\value{
the value of z corresponding to the user-defined value of MNPL.
}
\description{
Calculate the parameter z, the degree of compensation
}
\details{
Helper function for calculating z when user specifies MNPL
}
\examples{
calc_z(MNPL_in = 0.5,
lh.params_in = list(S0 = 0.944, S1plus = 0.99, AgeMat = 17, nages = 19,
fmax = 0.29, lambdaMax = 1.04, K1plus = 9000))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/25_abund_rel_vec.R
\name{abund_rel_vec}
\alias{abund_rel_vec}
\title{Calculate relative abundance}
\usage{
abund_rel_vec(traj.list, K = NA, years.vec = 10)
}
\arguments{
\item{traj.list}{a list of simulation outputs (high.list,med.list,low.lis,zero.list) each of these lists is a list of three, one output at each starting depletion (high, med, low depletion)}

\item{K}{carrying capacity}

\item{years.vec}{vector of years to check abundances at. If length = 1, the fn returns relative abundance at that year}
}
\value{
A vector of relative abundances sorted by bycatch (high, med, low, zero) within depletion level (low, med, high).
}
\description{
This function returns a vector of relative abundances based on a list of projection outputs. It is primarily for building the performance table.
}
\examples{
parms <- list(S0 = 0.944, S1plus = 0.99, K1plus = 9000, AgeMat = 18, 
              nages = 25, z = 2.39, lambdaMax = 1.02)
nyears <- 50
initdepl.vec <- c(0.2, 0.5, 0.9)
high.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 25, CV = 0.3),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
med.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 12, CV = 0.3),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
low.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 2, CV = 0.3),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
zero.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 0, CV = 0),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
traj.list <- list(
  high.list.const,
  med.list.const,
  low.list.const,
  zero.list.const
)

abund_rel_vec(traj.list = traj.list, K = parms$K1plus, years.vec = 10)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/14_plot_yield_curve.R
\name{plot_yield_curve}
\alias{plot_yield_curve}
\title{Yield curve for Shiny}
\usage{
plot_yield_curve(lh.params, z, MNPL_in, lang = "en")
}
\arguments{
\item{lh.params}{a list of life history parameters: S0, S1plus, AgeMat, nages, lambdaMax, K1plus}

\item{z}{degree of compensation. If this function is used outside the Shiny app, z is calculated from MNPL.}

\item{MNPL_in}{Maximum Net Productivity Level (MNPL) defined as the greatest net annual increment in population numbers or biomass resulting from additions to the population due to reproduction and/or growth less losses due to natural mortality. If the function is used outside Shiny, it will calculate z from this value. In Shiny, only z is used for productivity.}

\item{lang}{language selected by the user (character)}
}
\value{
a ggplot object showing depletion (1+ population size relative to K) vs. production. In fisheries this is a yield curve; in marine mammal management it shows where the productivity level is highest, i.e., MNPL.
}
\description{
Yield curve for Shiny
}
\examples{
plot_yield_curve(
  lh.params = list(
    S0 = 0.944, S1plus = 0.99, AgeMat = 17,
    nages = 19, lambdaMax = 1.02, K1plus = 9000
  ),
  MNPL_in = 0.5, z = NA, lang = "en"
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/15_plot_bycatch_guesses.R
\name{plot_bycatch_guesses}
\alias{plot_bycatch_guesses}
\title{Plot distribution of bycatch values}
\usage{
plot_bycatch_guesses(
  highval,
  medval,
  lowval,
  cv,
  set_size = 18,
  color.palette = c("forestgreen", "orange", "red"),
  lang = "en"
)
}
\arguments{
\item{highval}{The high end of the user-defined range of bycatch values}

\item{medval}{The middle of the user-defined bycatch range}

\item{lowval}{The low end of the user-defined range}

\item{cv}{CV of bycatch mortality or bycatch mortality rate}

\item{set_size}{base size to pass to \code{ggplot}}

\item{color.palette}{A vector of colors to represent bycatch levels, from low to high end of range}

\item{lang}{Language selected by the user (character, 2 letters)}
}
\value{
A \code{ggplot2} grob showing the distributions of bycatch values based on what the user entered.
}
\description{
Plot user-defined bycatch values (uses \code{geom_linerange()} from \code{ggplot2})
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/17_plot_proj.R
\name{plot_proj}
\alias{plot_proj}
\title{Plot projections}
\usage{
plot_proj(
  high,
  med,
  low,
  years.plot = 50,
  ylims,
  spaghetti = FALSE,
  K1plus = 9000,
  InitDepl = 0.8,
  color.palette = c("#7bbcb0", "#3a7c89", "#123f5a"),
  lang = "en"
)
}
\arguments{
\item{high}{a list containing output from \code{projections()} (including a matrix of simulation trajectories) - corresponding to a high bycatch level (this is at the high end of the range defined by the user in the Shiny app)}

\item{med}{a list containing output from \code{projections()} (including a matrix of simulation trajectories) - corresponding to a medium bycatch level (this is at the median of the low and high end of the range defined by the user)}

\item{low}{a list containing output from \code{projections()} (including a matrix of simulation trajectories) - corresponding to a high bycatch level (this is at the high end of the range defined by the user)}

\item{years.plot}{number of years to plot on the x axis}

\item{ylims}{y-limits of projection plot}

\item{spaghetti}{either FALSE or a number, where the number is how many simulations to show from the same scenario}

\item{K1plus}{carrying capacity in terms of age 1+ individuals}

\item{InitDepl}{initial depletion level (1+ population size relative to K). Must be between 0 and 1.}

\item{color.palette}{a vector of three colors to use for low, medium and high bycatch rates}

\item{lang}{language selected by the user (character)}
}
\value{
A plot of 50 percent and 90 percent confidence intervals of population projections (if \code{spaghetti == FALSE}) or a spaghetti plot (if \code{is.numeric(spaghetti)}),  from \code{Projections()}.
}
\description{
plots outputs from several projections that result from a Projections() call.
}
\examples{
parms <- list(
  S0 = 0.944, S1plus = 0.99,
  K1plus = 9000, 
  AgeMat = 18, nages = 20,
  z = 2.39, lambdaMax = 1.02
)
initdepl <- 0.5
high.simple <- projections(
  NOut = 50,
  ConstantBycatch = list(
    Catch = 100,
    CV = 0.3
  ),
  InitDepl = initdepl,
  lh.params = parms,
  nyears = 100
)
med.simple <- projections(
  NOut = 50,
  ConstantBycatch = list(
    Catch = 50,
    CV = 0.3
  ),
  InitDepl = initdepl,
  lh.params = parms,
  nyears = 100
)
low.simple <- projections(
  NOut = 50,
  ConstantBycatch = list(
    Catch = 10,
    CV = 0.3
  ),
  InitDepl = initdepl,
  lh.params = parms,
  nyears = 100
)

x <- plot_proj(
  high = high.simple,
  med = med.simple,
  low = low.simple,
  years.plot = 50,
  ylims = c(0, parms$K1plus), InitDepl = initdepl,
  K1plus = parms$K1plus
)

x
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/26_make_ptable.R
\name{make_ptable}
\alias{make_ptable}
\title{Make a performance table. Trajectories must be based on 100-year simulations.}
\usage{
make_ptable(traj.list, depletion, mnpl = NA)
}
\arguments{
\item{traj.list}{A list of trajectories from \code{Projections()}}

\item{depletion}{A vector of starting depletions (abundance relative to carrying capacity)}

\item{mnpl}{Max net productivity level (MNPL) defined by user. If MNPL is specified, that value is used. If the user does not specify MNPL, it is calculated from the life history parameters.}
}
\value{
A dataframe containing performance metrics.
}
\description{
Make a performance table. Trajectories must be based on 100-year simulations.
}
\examples{
parms <- list(S0 = 0.944, S1plus = 0.99, K1plus = 9000, AgeMat = 18, nages = 20,
              z = 2.39, lambdaMax = 1.02)
initdepl.vec <- c(0.2, 0.5, 0.9)
nyears <- 100
high.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 25, CV = 0.3),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
med.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 12, CV = 0.3),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
low.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 2, CV = 0.3),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
zero.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 0, CV = 0),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
traj.list <- list(
  high.list.const,
  med.list.const,
  low.list.const,
  zero.list.const
)
make_ptable(traj.list = traj.list, depletion = initdepl.vec, mnpl = 0.5)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/22_plot_pinnipeds.R
\name{plot_pinnipeds}
\alias{plot_pinnipeds}
\title{Plot all the pinniped life history parameters}
\usage{
plot_pinnipeds(dat, central = FALSE, set_size = 14)
}
\arguments{
\item{dat}{dat a dataframe with all the parameters (from Pinniped parameters - Data.csv)}

\item{central}{Whether or not to plot a line showing the central tendency of all the parameters}

\item{set_size}{Base size for plot}
}
\value{
a plotly interactive plot showing life history parameters for pinnipeds.
}
\description{
Make the table on the documentation tab of the app, which shows the range of pinniped life history parameters
}
\examples{
plot_pinnipeds(dat = dat)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/20_dynamics.R
\name{dynamics}
\alias{dynamics}
\title{Generate one marine mammal population trajectory}
\usage{
dynamics(lh.params, InitDepl, ConstantCatch = NA, ConstantF = NA, nyears)
}
\arguments{
\item{lh.params}{A list containing life history parameters: \cr\cr
\code{S0} Calf/pup survival, a numeric value between 0 and 1 \cr\cr
\code{S1plus} Survival for animals age 1 year and older, a numeric value between 0 and 1 \cr\cr
\code{K1plus} The pre-exploitation population size of individuals aged 1 and older.  If this value is unavailable, it can be approximated by using the initial depletion and the estimate of current abundance \cr\cr
\code{AgeMat} Age at maturity in years (assumed to be age at first parturition - 1) \cr\cr
\code{nages} "Maximum" age, treated as the plus group age. The plus group age can be set equal to the age at maturity +2 years without losing accuracy. Must be greater than \code{AgeMat}.\cr\cr
\code{z} The degree of compensation.  The default value is \code{z = 2.39}.\cr\cr
\code{lambdaMax} Maximum steady rate of increase (population growth rate)}

\item{InitDepl}{Starting depletion level}

\item{ConstantCatch}{Total bycatch each year, expressed as a vector of length \code{nyears}}

\item{ConstantF}{vector (length = \code{nyears}) rate of bycatch each year}

\item{nyears}{Number of years to project}
}
\value{
A list containing a matrix \code{N} of numbers at age (dimensions \code{nyears} (rows) x \code{nages} (columns)) and one vector \code{TotalPop} (a vector of length \code{nyears}), containing the number of age 1+ individuals in the population.
}
\description{
This function generates one trajectory for a marine mammal population, starting at a user-specified depletion level \code{InitDepl}.
}
\details{
The population model is a single-sex age-structured model in which the number of calves or pups born each year is density dependent, with the extent of density dependence a function of the number of mature adults \eqn{\tildeN}, the fecundity (pregnancy rate) at pre-exploitation equilibrium \eqn{f_0}, the maximum theoretical fecundity rate fmax, the degree of compensation \eqn{z}, and the abundance of individuals aged 1+ \eqn{N_{y+1}^{1+}} relative to carrying capacity \eqn{K^{1+}}. This function can be used alone but is intended to be used with \code{Projections()} to generate multiple simulations. NOTE: Either \code{ConstantCatch} or \code{ConstantF} can be specified, but not both.
}
\examples{
# Generate a time series of abundance for a bowhead whale
dynamics(lh.params = list(S0 = 0.944, S1plus = 0.99, 
K1plus = 9000, AgeMat = 17,nages = 25,
z = 2.39, lambdaMax = 1.04),
 InitDepl = 0.6, ConstantCatch = NA, ConstantF = rep(0.01, times = 100), 
 nyears = 100)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/28_traj_list_to_df.R
\name{traj_list_to_df}
\alias{traj_list_to_df}
\title{Turn a nested list of projections into a single dataframe}
\usage{
traj_list_to_df(x)
}
\arguments{
\item{x}{A nested list of projections, with depletion levels nested within bycatch levels.}
}
\value{
a dataframe containing projection outputs, initial depletion levels, and bycatch levels.
}
\description{
This function is for any case where the user wants to pull the projection outputs out of their nested list format. In the Shiny app, it is used to generate a simple table of raw outputs in case users want to create their own plots.
}
\examples{
 parms <- list(S0 = 0.944, S1plus = 0.99, K1plus = 9000, AgeMat = 18,
nages = 25, z = 2.39, lambdaMax = 1.02)
nyears <- 50
initdepl.vec <- c(0.2, 0.5, 0.9)
high.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 25, CV = 0.3),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
med.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 12, CV = 0.3),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
low.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 2, CV = 0.3),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
zero.list.const <- lapply(
  X = initdepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantBycatch = list(Catch = 0, CV = 0),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.1
    )
  }
)
traj.list <- list(
  high.list.const,
  med.list.const,
  low.list.const,
  zero.list.const
)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/27_rebuild_by_x.R
\name{rebuild_by_x}
\alias{rebuild_by_x}
\title{Calculate the bycatch rate needed to reach a performance goal}
\usage{
rebuild_by_x(
  needf.start,
  init.depl.w,
  goal.w,
  desired.prob.w,
  when.w,
  lh.params.w,
  fixed.cv.catch.w
)
}
\arguments{
\item{needf.start}{Starting guess for the bycatch mortality rate needed to recover the population}

\item{init.depl.w}{Initial depletion (a fraction)}

\item{goal.w}{Population size goal (number of whales or pinnipeds) to rebuild to, expressed as a whole number}

\item{desired.prob.w}{What probability you want (e.g., 0.75 probability that the population will rebuild to X in Y years)}

\item{when.w}{The year Y when rebuilding is desired by (number of years into future; current year = 0)}

\item{lh.params.w}{A list of life history parameters used by \code{projections()}}

\item{fixed.cv.catch.w}{The CV of the bycatch rate - should be fixed}
}
\value{
The bycatch rate that will result in the specified rebuilding goal. If Shiny is running, it will return a list containing the bycatch rate \code{f} and a matrix of guesses that \code{optim()} has searched through to find the solution.
}
\description{
Takes a performance goal (what level you want to rebuild to) and a time window (how long you want that to take) and calculates what the bycatch rate needs to be
}
\examples{
rebuild_by_x(
  needf.start = 0.001,
  init.depl.w = 0.5, goal.w = 4500,
  desired.prob.w = 0.8, when.w = 100,
  lh.params.w = list(
    S0 = 0.944, S1plus = 0.99,
    AgeMat = 17, nages = 19,
    z = 2.39, lambdaMax = 1.04, K1plus = 9000
  ),
  fixed.cv.catch.w = 0
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/09_find_msyr.R
\name{find_msyr}
\alias{find_msyr}
\title{calculate MSYR}
\usage{
find_msyr(E.start, lh.params, fmax)
}
\arguments{
\item{E.start}{a starting guess for bycatch mortality rate that will result in FMSY (numeric value)}

\item{lh.params}{a list of life history parameters (S0, S1plus, nages, AgeMat, lmabdaMax, K1plus, and z)}

\item{fmax}{Max theoretical fecundity (numeric)}
}
\description{
one of the functions needed for getting MNPL. MSYR is the same as FMSY.
}
\examples{
# Set parameters
S0.w = 0.5; S1plus.w = 0.944; nages.w = 25; K1plus.w = 9000; AgeMat.w = 18 
InitDepl.w = 0.9; z.w = 2.39; lambdaMax.w = 1.04
lh.params = list(S0 = S0.w,S1plus = S1plus.w, nages = nages.w,
AgeMat = AgeMat.w,K1plus=9000,z=z.w,lambdaMax = lambdaMax.w) 
# Get number of individuals per recruit in terms of mature individuals (\eqn{N0.w})
NPROut <- npr(S0 = S0.w, S1plus = S1plus.w, nages = nages.w, AgeMat = AgeMat.w, E = 0)
N0 <- NPROut$npr # mature numbers per recruit
# Get number of individuals per recruit in terms of individuals aged 1+ (\eqn{P0.w})
P0 <- NPROut$P1r # 1+ nums per recruit

fmax <- getfecmax(lambdaMax = lambdaMax.w, S0 = S0.w, S1plus = S1plus.w, AgeMat = AgeMat.w)
Fec0 <- 1.0 / N0
A <- (fmax - Fec0) / Fec0
fmsy<-find_msyr(E.start=0.01, lh.params=lh.params, fmax=fmax)
cat("fmsy =",fmsy,"\n")
results <- matrix(0,ncol=2,nrow=101)
results[,1] <- c(0:100)*fmsy/50
for (II in 1:101)
 results[II,2] <- ce(S0 = S0.w, S1plus = S1plus.w, nages = nages.w, AgeMat = AgeMat.w,z = z.w,
                     E=results[II,1], A = A, P0 = P0, N0 = N0)
plot(results[,1],results[,2],xlab="f",ylab="ce",type="l",yaxs="i",ylim=c(0,max(results[,2])*1.2))
abline(v=fmsy)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/21_projections.R
\name{projections}
\alias{projections}
\title{Run simulations of a marine mammal population with bycatch mortality}
\usage{
projections(
  NOut,
  ConstantBycatch = list(Catch = NA, CV = NA),
  ConstantRateBycatch = list(Rate = NA, CV = NA),
  InitDepl,
  lh.params,
  nyears,
  obs_CV = 0
)
}
\arguments{
\item{NOut}{Number of simulations}

\item{ConstantBycatch}{Mean and CV of number of animals killed as bycatch per year (assumed lognormal)}

\item{ConstantRateBycatch}{Mean and CV of bycatch rate (assumed normal)}

\item{InitDepl}{Initial depletion. If obs_CV>0, this is the mean depletion.}

\item{lh.params}{Life history parameters as a list. The list must include S0, S1plus, K1plus, AgeMat, nages, z, and lambdaMax.}

\item{nyears}{Number of years to do projections}

\item{obs_CV}{Observation CV. Default to 1 for simple projections}
}
\value{
A list of outputs from simulations: \cr\cr
* \code{params} contains parameter values for each trajectory as a matrix; \cr
* \code{trajectories} contains simulation outputs as a matrix; \cr
* \code{fishing.rates} contain the bycatch rates for each year in each simulation as a matrix; \cr
* \code{InitDepl} returns the initial depletion for the projections; \cr
* \code{ConstantBycatch} provides Catch (total individuals killed in bycatch events per year) and CV of Catch (if the user has specified bycatch as a constant number); \cr
* \code{ConstantRateBycatch} contains Bycatch Rate (additional mortality from bycatch each year) and CV of ByCatch rate. Other parameters are the same as in the \code{dynamics()} function.
}
\description{
Generates several projections, stochasticity is in the number of catches from year to year
}
\examples{
projections(
  NOut = 3, ConstantRateBycatch = list(Rate = 0.01, CV = 0.3),
  InitDepl = 0.8,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, 
    nages = 20, z = 2.39, 
    lambdaMax = 1.02
  ),
  nyears = 50, obs_CV = 0
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_get_rf.R
\name{get_rf}
\alias{get_rf}
\title{Get relative recruitment at bycatch mortality rate E}
\usage{
get_rf(E_in, S0, S1plus, nages, AgeMat, z, A, P0, N0)
}
\arguments{
\item{E_in}{Bycatch mortality rate, a numeric value between 0 and 1}

\item{S0}{calf/pup survival, a numeric value between 0 and 1}

\item{S1plus}{adult survival, a numeric value between 0 and 1}

\item{nages}{number of age classes-- including plus group age-- in years}

\item{AgeMat}{age at maturity in years, must be less than nages}

\item{z}{degree of compensation, in the app calculated from the value of MNPL defined by the user}

\item{A}{Pella-Tomlinson resilience parameter (see Punt 1999; Annex R).
A = (FecMax - Fec0) / Fec0 (num)}

\item{P0}{unfished 1+ numbers per recruit, \eqn{\tildeP(0)}}

\item{N0}{unfished mature numbers per recruit, \eqn{\tildeN(0)}}
}
\value{
recruitment given exploitation rate \emph{E} - this value is multiplied by the initial abundance \eqn{N_{init}} to get initial nums at age ( a vector)
}
\description{
calculates recruitment at bycatch mortality rate \emph{E}, relative to recruitment for a no-bycatch scenario
}
\examples{
S0 = 0.944; S1plus = 0.99; nages = 25; AgeMat = 17; z = 2.39; lambdaMax = 1.04;
NPROut <- npr(S0 = S0, S1plus = S1plus, nages = nages, AgeMat = AgeMat, E = 0)
N0 <- NPROut$npr # mature numbers per recruit

P0 <- NPROut$P1r # 1+ nums per recruit
Fec0 <- 1.0 / N0
FecMax <- getfecmax(S1plus = S1plus, S0 = S0, 
AgeMat = AgeMat, lambdaMax = lambdaMax)
A <- (FecMax - Fec0) / Fec0

get_rf(E_in = 0.01, S0 = S0, S1plus = S1plus, 
nages = nages, AgeMat = AgeMat, z = z, A = A,P0 = P0,N0 = N0)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_logit.R
\name{logit}
\alias{logit}
\title{Logit transform}
\usage{
logit(p)
}
\arguments{
\item{p}{Value to logit-transform}
}
\description{
Take the logit transform of a value \emph{p}
}
\examples{
x <- logit(p = 0.01)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/08_get_diff.R
\name{get_diff}
\alias{get_diff}
\title{Calculate the derivative of C in getMNPL()}
\usage{
get_diff(
  logit.E,
  S0 = S0.w,
  S1plus = S1plus.w,
  AgeMat = AgeMat.w,
  nages = nages.w,
  A = A.w,
  z = z.w
)
}
\arguments{
\item{logit.E}{logit transform of bycatch mortality}

\item{S0}{Calf/pup survival, a numeric value between 0 and 1}

\item{S1plus}{adult survival, a numeric value between 0 and 1}

\item{AgeMat}{Age at maturity (= age at first parturition - 1)}

\item{nages}{"maximum" age, treated as the plus group age. The plus group age can be set equal to the age at maturity +2 years without losing accuracy.}

\item{A}{the Pella-Tomlinson resilience parameter ((fmax - f0)/f0)}

\item{z}{Degree of compensation (also known as the Pella-Tomlinson parameter)}
}
\description{
this is the function to minimize, but needs to be restricted to that yield > 0
All life history params are as above
}
\examples{
S0.w = 0.5; S1plus.w = 0.944; nages.w = 25; AgeMat.w = 18 
InitDepl.w = 0.9; z.w = 2.39; lambdaMax.w = 1.04
# Get A parameter
NPROut <- npr(S0 = S0.w, S1plus = S1plus.w, nages = nages.w, AgeMat = AgeMat.w, E = 0)
N0 <- NPROut$npr # mature numbers per recruit
Fec0 <- 1.0 / N0
fmax <- getfecmax(lambdaMax = lambdaMax.w, S0 = S0.w, S1plus = S1plus.w, AgeMat = AgeMat.w)
A.w <- (fmax - Fec0) / Fec0
# Get number of individuals per recruit in terms of mature individuals (\eqn{N0.w})
get_diff(
 logit.E = logit(0.01), S0 = S0.w, S1plus = S1plus.w, nages = nages.w, A = A.w, AgeMat = AgeMat.w, 
  z = 2.39
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_app.R
\name{run_app}
\alias{run_app}
\title{Run the Shiny Application}
\usage{
run_app(...)
}
\arguments{
\item{...}{A series of options to be used inside the app.}
}
\description{
Run the Shiny Application
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10_get_mnpl.R
\name{get_mnpl}
\alias{get_mnpl}
\title{Calculate Maximum Net Productivity Level (MNPL)}
\usage{
get_mnpl(E.start = 0.001, lh.params)
}
\arguments{
\item{E.start}{a starting guess for bycatch mortality rate that will result in MSYR (equivalent mathematically to FMSY for fisheries). A numeric value between 0 and 1.}

\item{lh.params}{a list of life history parameters}
}
\description{
Calculate Maximum Net Productivity Level (MNPL)
}
\examples{
get_mnpl(E.start = 0.001, 
lh.params = list(S0 = 0.944, S1plus = 0.99, AgeMat = 17, nages = 19,
z = 2.39, lambdaMax = 1.04, 
K1plus = 9000))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_add_trans.R
\name{add_trans}
\alias{add_trans}
\title{Add transparency to a color}
\usage{
add_trans(color, trans)
}
\arguments{
\item{color}{input color}

\item{trans}{degree of transparency - a numeric values between 0 and 255 where 255 is fully visible)}
}
\value{
a hex code for the transparent version of \code{color}
}
\description{
Add transparency to a color
}
\details{
This function adds transparency to a color. It is originally by Tim Essington. Define transparency with an integer between 0 and 255, 0 being fully transparent and 255 being fully visible. Works with either color and trans a vector of equal length, or one of the two of length 1.
}
\examples{
x <- add_trans("red", 100)
print(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/24_abund_rel.R
\name{abund_rel}
\alias{abund_rel}
\title{Calculate relative abundance}
\usage{
abund_rel(traj, zero.traj, K = NA, years.vec = c(10, 20, 50), fulldist = TRUE)
}
\arguments{
\item{traj}{a matrix of trajectories, with rows=nsims and cols=nyears. Can be produced from \code{projections()}}

\item{zero.traj}{a matrix of trajectories, same as traj but unfished. CHECK to make sure both start at the same simulation, in the app code.}

\item{K}{if calculating abundance relative to K, put in K.}

\item{years.vec}{a vector of years to check abundances at. E.g., how many years after the start of projections do you want to know abundance?}

\item{fulldist}{logical saying whether to return the full distribution of relative abundances or the median of all the relative abundances.}
}
\value{
a vector of abundance relative to K or zero-exploitation where nrows=length of years.vec (10 years after, 20 years, etc.) and ncol=number of simulations.
}
\description{
Calculates the expected abundance relative to zero exploitation, years.vec years after projections start
}
\examples{
parms <- list(
  S0 = 0.944, S1plus = 0.99,
  K1plus = 9000,
  AgeMat = 18, nages = 20,
  z = 2.39, lambdaMax = 1.02
)
initdepl <- 0.5
traj <- projections(
  NOut = 50,
  ConstantBycatch = list(
    Catch = 50,
    CV = 0.3
  ),
  InitDepl = initdepl,
  lh.params = parms,
  nyears = 100
)$trajectories
traj0 <- projections(
  NOut = 50,
  ConstantBycatch = list(
    Catch = 0,
    CV = 0
  ),
  InitDepl = initdepl,
  lh.params = parms,
  nyears = 100
)$trajectories
abund_rel(traj = traj, zero.traj = traj0, K = parms$K1plus)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/18_multiplot_proj.R
\name{multiplot_proj}
\alias{multiplot_proj}
\title{multiplot.proj}
\usage{
multiplot_proj(
  high.d1,
  med.d1,
  low.d1,
  high.d2,
  med.d2,
  low.d2,
  high.d3,
  med.d3,
  low.d3,
  spaghetti = FALSE,
  years.to.plot = 50,
  color.palette = c("forestgreen", "darkorange", "red"),
  lang = "en"
)
}
\arguments{
\item{high.d1}{a list containing output from \code{projections()} (including a matrix of simulation trajectories) - this corresponds to the highest bycatch level (thus "high") and the depletion level ( \code{d1} indicates the lowest starting depletion level)}

\item{med.d1}{a list containing the middle bycatch value and lowest starting depletion}

\item{low.d1}{a list containing the lowest bycatch value and lowest starting depletion}

\item{high.d2}{a list containing the highest bycatch value and middle starting depletion}

\item{med.d2}{a list containing the middle bycatch value and middle starting depletion}

\item{low.d2}{a list containing the lowest bycatch value and middle starting depletion}

\item{high.d3}{a list containing the highest bycatch value and highest starting depletion}

\item{med.d3}{a list containing the middle bycatch value and highest starting depletion}

\item{low.d3}{a list containing the lowest bycatch value and highest starting depletion}

\item{spaghetti}{either FALSE or a number, where the number is how many simulations to show from the same scenario}

\item{years.to.plot}{How many years to plot on the x axis}

\item{color.palette}{A vector of length 3 giving the color values for low, medium, and high bycatch mortality or bycatch mortality rate}

\item{lang}{language to use. "en" = English; "es" = Spanish; "fr" = French.}
}
\value{
A plot of 50 percent and 90 percent confidence intervals of population projections if \code{spaghetti == FALSE} or a spaghetti plot with n individual projections if \code{spaghetti == n },  from \code{projections()}.
}
\description{
plots outputs from several projections that result from \code{projections()}, with multiple depletion levels.
}
\examples{
parms <- list(
  S0 = 0.944, S1plus = 0.99,
  K1plus = 9000,
  AgeMat = 18, nages = 20,
  z = 2.39, lambdaMax = 1.02
)
InitDepl.vec <- c(0.1, 0.5, 0.9)
BycatchCV <- 0.2
nyears <- 100

high.list <- lapply(
  X = InitDepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantRateBycatch = list(Rate = 0.3, CV = BycatchCV),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.2
    )
  }
)

med.list <- lapply(
  X = InitDepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantRateBycatch = list(Rate = 0.02, CV = BycatchCV),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.2
    )
  }
)
low.list <- lapply(
  X = InitDepl.vec,
  function(x) {
    projections(
      NOut = 50,
      ConstantRateBycatch = list(Rate = 0.001, CV = BycatchCV),
      InitDepl = x,
      lh.params = parms,
      nyears = nyears,
      obs_CV = 0.2
    )
  }
)
multiplot_proj(
  high.d1 = high.list[[1]], # d1 is the lowest depletion
  med.d1 = med.list[[1]],
  low.d1 = low.list[[1]],
  high.d2 = high.list[[2]],
  med.d2 = med.list[[2]],
  low.d2 = low.list[[2]],
  high.d3 = high.list[[3]],
  med.d3 = med.list[[3]],
  low.d3 = low.list[[3]],
  years.to.plot = nyears
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/28_extract_df.R
\name{extract_df}
\alias{extract_df}
\title{Extract results from `projections()` function as a dataframe}
\usage{
extract_df(x)
}
\arguments{
\item{x}{Outputs from a call to projections(), which are a named list}
}
\value{
a dataframe with population trajectories, initial depletion, and bycatch levels that can be printed as a table.
}
\description{
Unnest results from projections() and turn them into a dataframe that you can bind together
}
\examples{
x <- projections(
  NOut = 3, ConstantRateBycatch = list(Rate = 0.01, CV = 0.3),
  InitDepl = 0.8,
  lh.params = list(
    S0 = 0.944, S1plus = 0.99,
    K1plus = 9000, AgeMat = 18, nages = 20, z = 2.39, lambdaMax = 1.02
  ),
  nyears = 50, obs_CV = 0
)
extract_df(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lh-data.R
\docType{data}
\name{lh}
\alias{lh}
\title{Life history parameters for marine mammals.}
\format{
A data frame with 11 rows and 11 variables:
\describe{
  \item{Type}{Life history type}
  \item{Code}{a quick reference for that life history type}
  \item{Representative}{The species representing that life history type}
  \item{S0}{Calf/pup survival}
  \item{S1plus}{Adult survival}
  \item{AgeMat}{Age at maturity}
  \item{nages}{Plus group age; referred to as nages in code}
  \item{fmax}{Maximum theoretical fecundity}
  \item{z}{Degree of compensation}
  \item{lambdaMax}{Maximum theoretical population growth rate}
  \item{K1plus}{Carrying capacity in terms of the 1+ component of the population}
  ...
}
}
\source{
Compiled from literature by M. Siple and coauthors
}
\usage{
lh
}
\description{
A table with the default parameter values for the marine mammals included in this app.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_getfecmax.R
\name{getfecmax}
\alias{getfecmax}
\title{Calculate maximum theoretical fecundity rate}
\usage{
getfecmax(S0, lambdaMax, S1plus, AgeMat)
}
\arguments{
\item{S0}{calf/pup survival, a numeric value between 0 and 1}

\item{lambdaMax}{maximum population growth rate (must exceed 0; default value is 1.04 for cetaceans and 1.12 for pinnipeds)}

\item{S1plus}{survival of age 1+ individuals, a numeric value between 0 and 1}

\item{AgeMat}{age at maturity in years (must be equal to or less than \code{nages})}
}
\value{
a numeric value for maximum theoretical fecundity.
}
\description{
Calculate maximum theoretical fecundity rate \emph{fmax}
}
\details{
Parts of this derivation can be found in Breiwick et al. (1984) and Butterworth and Punt (1992).
\emph{Important}: when applying this calculation, use age at maturity ( \code{AgeMat} )
}
\examples{
# This is fmax, the maximum theoretical fecundity
x <- getfecmax(lambdaMax = 1.04, S0 = 0.944, S1plus = 0.99, AgeMat = 17)
# This is fec0, the fecundity when there is no bycatch mortality (only M)
unpr <- npr(S0 = 0.944, S1plus = 0.99, AgeMat = 17, nages = 10000, E = 0)
print(x)
print(1 / unpr$npr)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_npr.R
\name{npr}
\alias{npr}
\title{Calculate numbers per recruit}
\usage{
npr(S0, S1plus, nages, AgeMat, E = 0)
}
\arguments{
\item{S0}{Calf/pup survival, a numeric value between 0 and 1}

\item{S1plus}{Adult survival, a numeric value between 0 and 1}

\item{nages}{Plus group age in years}

\item{AgeMat}{Age at maturity in years (must be equal to or less than nages)}

\item{E}{Bycatch mortality rate, a numeric value between 0 and 1}
}
\value{
A list of numbers per recruit (\code{npr}), 1+ numbers per recruit (\code{P1r}), and numbers at age per recruit (\code{nvec})
}
\description{
Calculate numbers-per-recruit as a function of the bycatch rate assuming that 1+ animals are subject to bycatch  \emph{E}.
}
\examples{
(unpr <- npr(S0 = 0.9, S1plus = 0.9, 
AgeMat = 11, nages = 13, E = 0)) # unfished nums per recruit
(nprf <- npr(S0 = 0.9, S1plus = 0.9, 
AgeMat = 11, nages = 13, E = 0.8)) # nums per recruit at bycatch rate E
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_get_f.R
\name{get_f}
\alias{get_f}
\title{Get bycatch mortality rate given depletion}
\usage{
get_f(
  f.start = NA,
  S0.w = NA,
  S1plus.w = NA,
  nages.w = NA,
  AgeMat.w = NA,
  InitDepl.w = NA,
  z.w = NA,
  lambdaMax.w = NA,
  N0.w = NA,
  P0.w = NA,
  Check = FALSE
)
}
\arguments{
\item{f.start}{an initial guess for the bycatch mortality rate E. The default value is E = 0.5}

\item{S0.w}{Calf/pup survival, a numeric value between 0 and 1. (Note: the 'w' suffix indicates that \eqn{z} is in the wrapper function, and is used inside the function by \code{optim})}

\item{S1plus.w}{survival rate for animals age 1 year and older, a numeric value between 0 and 1}

\item{nages.w}{"maximum" age, treated as the plus group age. The plus group age can be set equal to the age at maturity +2 years without losing accuracy.}

\item{AgeMat.w}{Age at maturity in years (assumed to be age at first parturition - 1)}

\item{InitDepl.w}{The depletion level to solve for as a proportion of carrying capacity; a numeric value between 0 and 1. This is equivalent to 'starting depletion'.}

\item{z.w}{The degree of compensation. The default value is \code{z = 2.39}.}

\item{lambdaMax.w}{The maximum intrinsic growth rate}

\item{N0.w}{unfished numbers per recruit in terms of mature individuals individuals}

\item{P0.w}{Number of individuals per recruit in terms of individuals aged 1+}

\item{Check}{logical; if \code{TRUE}, prints the value to be minimized (should be zero)}
}
\value{
The bycatch rate that would lead to a depletion level of \code{InitDepl.w} -- a single value.
}
\description{
This function solves for the bycatch mortality rate \eqn{E} that gives a pre-specified depletion level \code{InitDepl.w}. It is used within the \code{projections()} function to calculate the stable age distribution at which to start the projections.
}
\examples{
# Set parameters
S0.w = 0.5; S1plus.w = 0.944; 
nages.w = 25; AgeMat.w = 18 
InitDepl.w = 0.9; z.w = 2.39; lambdaMax.w = 1.04
# Get number of individuals per recruit in terms of mature individuals (N0.w)
NPROut <- npr(S0 = S0.w, S1plus = S1plus.w, 
nages = nages.w, AgeMat = AgeMat.w, E = 0)
N0 <- NPROut$npr # mature numbers per recruit
# Get number of individuals per recruit in terms of individuals aged 1+ (P0.w)
P0 <- NPROut$P1r # 1+ nums per recruit

# Get bycatch mortality rate for the initial depletion defined above
get_f(f.start = 0.5, 
S0.w = S0.w, S1plus.w = S1plus.w, nages.w = nages.w, AgeMat.w = AgeMat.w,
InitDepl.w = InitDepl.w, z.w = z.w, lambdaMax.w = lambdaMax.w, 
N0.w = N0, P0.w = P0, Check = FALSE)
}
