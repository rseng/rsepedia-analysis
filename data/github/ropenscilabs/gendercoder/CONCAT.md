
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gendercoder

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/gendercoder)](https://CRAN.R-project.org/package=gendercoder)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropenscilabs/gendercoder/workflows/R-CMD-check/badge.svg)](https://github.com/ropenscilabs/gendercoder/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropenscilabs/gendercoder/branch/master/graph/badge.svg)](https://codecov.io/gh/ropenscilabs/gendercoder?branch=master)
<!-- badges: end -->

The goal of gendercoder is to allow simple re-coding of free-text gender
responses. This is intended to permit representation of gender
diversity, while managing troll-responses and the workload implications
of manual coding.

## Installation

This package is not on CRAN. To use this package please run the
following code:

``` r
devtools::install_github("ropenscilabs/gendercoder")
library(gendercoder)
```

## Basic use

The gendercoder package permits the efficient re-coding of free-text
gender responses within a tidyverse pipeline. It contains two built-in
English output dictionaries, a default `manylevels_en` dictionary which
corrects spelling and standardises terms while maintaining the diversity
of responses and a `fewlevels_en` dictionary which contains fewer gender
categories, “man”, “woman”, “boy”, “girl”, and “sex and gender diverse”.

The core function, `gender_recode()`, takes 3 arguments,

-   `gender` the vector of free-text gender,

-   `dictionary` the preferred dictionary, and

-   `retain_unmatched` a logical indicating whether original values
    should be carried over if there is no match.

Basic usage is demonstrated below.

``` r
library(gendercoder)

tibble(gender = c("male", "MALE", "mle", "I am male", "femail", "female", "enby")) %>% 
  mutate(manylevels_gender  = recode_gender(gender, dictionary = manylevels_en, retain_unmatched = TRUE),
         fewlevels_gender = recode_gender(gender, dictionary = fewlevels_en, retain_unmatched = FALSE)
  )
#> Results not matched from the dictionary have been filled with the user inputted values
#> # A tibble: 7 x 3
#>   gender    manylevels_gender fewlevels_gender      
#>   <chr>     <chr>             <chr>                 
#> 1 male      man               man                   
#> 2 MALE      man               man                   
#> 3 mle       man               man                   
#> 4 I am male I am male         <NA>                  
#> 5 femail    woman             woman                 
#> 6 female    woman             woman                 
#> 7 enby      non-binary        sex and gender diverse
```

The package does not need to be used as part of a tidyverse pipeline:

``` r
df <- tibble(gender = c("male", "MALE", "mle", "I am male", "femail", "female", "enby")) 

df$manylevels_gender <- recode_gender(df$gender, dictionary = manylevels_en)
df
#> # A tibble: 7 x 2
#>   gender    manylevels_gender
#>   <chr>     <chr>            
#> 1 male      man              
#> 2 MALE      man              
#> 3 mle       man              
#> 4 I am male <NA>             
#> 5 femail    woman            
#> 6 female    woman            
#> 7 enby      non-binary
```

## Contributing to this package

This package is a reflection of cultural context of the package
contributors. We acknowledge that understandings of gender are bound by
both culture and time and are continually changing. As such, we welcome
issues and pull requests to make the package more inclusive, more
reflective of current understandings of gender inclusive languages
and/or suitable for a broader range of cultural contexts. We
particularly welcome addition of non-English dictionaries or of other
gender-diverse responses to the manylevels\_en and fewlevels\_en
dictionaries.

## Citation Information

Please cite this package as:

Jennifer Beaudry, Emily Kothe, Felix Singleton Thorn, Rhydwyn McGuire,
Nicholas Tierney and Mathew Ling (2020). gendercoder: Recodes Sex/Gender
Descriptions into a Standard Set. R package version 0.0.0.9000.
<https://github.com/ropenscilabs/gendercoder>

## Acknowledgement of Country

We acknowledge the Wurundjeri people of the Kulin Nation as the
custodians of the land on which this package was developed and pay
respects to elders past, present and future.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at emily.kothe@deakin.edu.au. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
## Test environments
* macOS X Catalina (gh-actions), R 4.0.5
* ubuntu 20.04 (gh-actions), R 4.0.5
* ubuntu 20.04 (gh-actions), dev 2021-04-05 r8014
* windows server 19 (gh-actions), R 4.0.5

## R CMD check results
Note: This is a new submission.

## Downstream dependencies

There are currently no downstream dependencies for this package
