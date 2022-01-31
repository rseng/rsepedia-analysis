
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
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
library(dplyr)
```
# gendercoder
<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/gendercoder)](https://CRAN.R-project.org/package=gendercoder)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropenscilabs/gendercoder/workflows/R-CMD-check/badge.svg)](https://github.com/ropenscilabs/gendercoder/actions)
[![Codecov test coverage](https://codecov.io/gh/ropenscilabs/gendercoder/branch/master/graph/badge.svg)](https://codecov.io/gh/ropenscilabs/gendercoder?branch=master)
<!-- badges: end -->


The goal of gendercoder is to allow simple re-coding of free-text gender 
responses. This is intended to permit representation of gender diversity, 
while managing troll-responses and the workload implications of manual coding. 

## Installation

This package is not on CRAN. To use this package please run the following code:

``` r
devtools::install_github("ropenscilabs/gendercoder")
library(gendercoder)
```
## Basic use

The gendercoder package permits the efficient re-coding of free-text gender 
responses within a tidyverse pipeline. It contains two built-in English output 
dictionaries, a default `manylevels_en` dictionary which corrects spelling and 
standardises terms while maintaining the diversity of responses and a 
`fewlevels_en` dictionary which contains fewer gender categories, "man", "woman", 
"boy", "girl", and "sex and gender diverse". 

The core function, `gender_recode()`, takes 3 arguments, 

- `gender` the vector of free-text gender,

- `dictionary` the preferred dictionary, and

- `retain_unmatched` a logical indicating whether original values should be carried over if 
there is no match. 

Basic usage is demonstrated below. 

```{r}
library(gendercoder)

tibble(gender = c("male", "MALE", "mle", "I am male", "femail", "female", "enby")) %>% 
  mutate(manylevels_gender  = recode_gender(gender, dictionary = manylevels_en, retain_unmatched = TRUE),
         fewlevels_gender = recode_gender(gender, dictionary = fewlevels_en, retain_unmatched = FALSE)
  )
  
```

The package does not need to be used as part of a tidyverse pipeline:

```{r}
df <- tibble(gender = c("male", "MALE", "mle", "I am male", "femail", "female", "enby")) 

df$manylevels_gender <- recode_gender(df$gender, dictionary = manylevels_en)
df
```


## Contributing to this package

This package is a reflection of cultural context of the package contributors. 
We acknowledge that understandings of gender are bound by both culture and time 
and are continually changing. As such, we welcome issues and pull requests to 
make the package more inclusive, more reflective of current understandings of 
gender inclusive languages and/or suitable for a broader range of cultural 
contexts. We particularly welcome addition of non-English dictionaries or of 
other gender-diverse responses to the manylevels_en and fewlevels_en dictionaries.

## Citation Information

Please cite this package as:

Jennifer Beaudry, Emily Kothe, Felix Singleton Thorn, Rhydwyn McGuire, Nicholas Tierney and Mathew Ling (2020).
  gendercoder: Recodes Sex/Gender Descriptions into a Standard Set. R package version 0.0.0.9000.
  https://github.com/ropenscilabs/gendercoder

## Acknowledgement of Country

We acknowledge the Wurundjeri people of the Kulin Nation as the custodians of 
the land on which this package was developed and pay respects to elders past, 
present and future.
---
title: "Adding to the dictionary"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Adding to the dictionary}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup, echo = FALSE, message=FALSE}
df <-  dplyr::tibble(gender = c("male", "enby", "womn", "mlae", "mann", "frau", "femme", 
                         "homme", "nin"))
```
## Outline

While the `gendercoder` dictionaries aim to be as comprehensive as possible, it
is inevitable that new typos and variations will occur in wild data. Moreover, 
at present, the dictionaries are limited to data the authors have had access to 
which has been collected in English. As such, if you are collecting data, you 
will at some point want to add to or create your own dictionaries (and if so, 
we strongly encourage contributions either as [a pull request via github](https://github.com/ropenscilabs/gendercoder), or by [raising an issue](https://github.com/ropenscilabs/gendercoder/issues/new) so the team can 
help).

## Adding to the dictionary

Let's say I have free-text gender data, but some of it is not in English.

```{r example-data, message=FALSE}
library(gendercoder)
library(dplyr)
df
```

I can create a new dictionary by creating a named vector, where the names are 
the raw, uncoded values, and the values are the desired outputs. This can then 
be used as the dictionary in the `recode_gender()` function.

```{r new-dictionary}
new_dictionary <- c(
  mann = "man", 
  frau = "woman", 
  femme = "woman", 
  homme = "man", 
  nin = "man")

df %>% 
  mutate(recoded_gender = recode_gender(gender, 
                                        dictionary = new_dictionary, 
                                        retain_unmatched = TRUE))

```
However, as you can see using just this new dictionary leaves a number of 
responses uncoded that the built-in dictionaries could handle. As the 
dictionaries are just vectors, we can simply concatenate these to use both at 
the same time. 

We can do this in-line...

```{r inline-augmentation}
df %>% 
  mutate(recoded_gender = recode_gender(gender, 
                                        dictionary = c(manylevels_en, new_dictionary), 
                                        retain_unmatched = TRUE))
```
Or otherwise we can create a new dictionary and call that later, useful if you 
might want to save an augmented dictionary for later use or for contributing to 
the package.
```{r stepped-augmentation}
manylevels_plus <-  c(manylevels_en, new_dictionary)

df %>% 
  mutate(recoded_gender = recode_gender(gender, 
                                        dictionary = manylevels_plus, 
                                        retain_unmatched = TRUE))
```

## Making it official

Let's say you are happy with your `manylevels_plus` dictionary and think it should be 
part of the `manylevels_en` dictionary in the package. All you need to do is [fork the 
gendercoder repo](https://docs.github.com/en/get-started/quickstart/fork-a-repo), 
[clone it to your local device](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository-from-github/cloning-a-repository), and then rename 
your vector and use the `usethis::use_data()` function to overwrite the `manylevels_en` 
dictionary as shown below. 

```{r replacing-dictionaries, eval=FALSE}
manylevels_en <-  manylevels_plus
usethis::use_data(manylevels_en, overwrite = TRUE)
```

Once you've pushed the changes to your fork, you can [make a pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork). Please tell us 
what you're adding so we know what to look out for and how to test it.  

---
title: "Introduction to gendercoder"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to gendercoder}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The goal of gendercoder is to allow simple recoding of free-text gender responses.

## Installation

This package is not on CRAN. To use this package please run the following code:

``` r
devtools::install_github("ropenscilabs/gendercoder")
library(gendercoder)
```

## Why would we do this?

Researchers who collect self-reported demographic data from respondents 
occasionally collect gender using a free-text response option. This has the 
advantage of respecting the gender diversity of respondents without prompting 
users and potentially including misleading responses. However, this presents a 
challenge to researchers in that some inconsistencies in typography and spelling 
create a larger set of responses than would be required to fully capture the 
demographic characteristics of the sample. 

For example, male participants may provide free-text responses as "male", "man", 
"mail", "mael". Non-binary participants may provide responses as "nonbinary", 
"enby", "non-binary", "non binary"

Manually coding of such free-text responses this is often not feasible with 
larger datasets. `gendercoder()` uses dictionaries of common 
misspellings to re-code free-text responses into a consistent set of 
responses. The small number of responses not automatically re-coded by 
gendercoder() can then be feasibly manually recoded. 

## Motivating example

`gendercoder()` includes a sample dataset with actual free-text 
responses to the question "What is your gender?" from a number of studies of 
English-speaking participants. The sample dataset includes responses from 7756
participants. Naive coding identifies 103 unique responses to this item. 


```{r message = FALSE, warning= FALSE}
library(gendercoder)
library(dplyr)

sample %>% 
  group_by(Gender) %>% 
  summarise(count = n()) %>% 
  arrange(-count) %>% 
  knitr::kable(caption = "Summary of gender categories before coding")
```

Recoding using the `gender_coder()` function classifies all but 28 responses 
into pre-defined response categories. 

```{r}
sample %>% 
  head(10) %>% 
  mutate(recoded_gender = recode_gender(gender = Gender, dictionary = manylevels_en)) %>% 
  knitr::kable(caption = "The manylevels_en dictionary applied to `head(sample)`")

sample %>% 
  mutate(recoded_gender = recode_gender(gender = Gender, dictionary = manylevels_en)) %>% 
  filter(!is.na(recoded_gender)) %>% 
  group_by(recoded_gender) %>% 
  summarise(count = n()) %>% 
  arrange(-count) %>% 
  knitr::kable(caption = "Summary of gender categories after use of the *manylevels_en* dictionary")

```

In this dataset unclassified responses are a mix of unusual responses and 
apparent response errors (e.g. numbers and symbols). While some of these are
genuinely missing (i.e. Gender = 40), other could be manually recoded, or added
to a custom dictionary. 

```{r}

sample %>% 
  mutate(recoded_gender = recode_gender(gender = Gender, dictionary = manylevels_en)) %>% 
  filter(is.na(recoded_gender)) %>% 
  knitr::kable(caption = "All responses not classified by the built-in dictionary")

```

## Options within the function

### dictionary

The package provides two built-in dictionaries. The use of these is controlled 
using the `dictionary` argument. The first `dictionary = manylevels_en` provides 
corrects spelling and standardises terms while maintaining the diversity of responses. 
This is the default dictionary for gendercoder() as it preserves as much gender
diversity as possible. 

However in some cases you may wish to collapse gender into a smaller set of 
categories by using the `fewlevels_en` dictionary (`dictionary = fewlevels_en`). This 
dictionary contains fewer gender categories, "man", "woman", 
"boy", "girl", and "sex and gender diverse". 

The "man" category includes all participants who indicate that they are

- male
- trans male  (including female to male transgender respondents)
- cis male

The "woman" category includes all participants who indicate that they are

- female
- trans female (including male to female transgender respondents)
- cis female

The "sex and gender diverse" category includes all participants who indicate 
that they are

- agender
- androgynous
- intersex
- non-binary
- gender-queer

```{r} 
sample %>% 
  head(10) %>% 
  mutate(recoded_gender = recode_gender(gender = Gender, dictionary = fewlevels_en)) %>% 
  knitr::kable(caption = "The fewlevels_en dictionary applied to `head(sample)`")

sample %>% 
  mutate(recoded_gender = recode_gender(gender = Gender, dictionary = fewlevels_en)) %>% 
  group_by(recoded_gender) %>% 
  summarise(count = n()) %>% 
  arrange(-count) %>% 
  knitr::kable(caption = "Summary of gender categories after use of the *fewlevels_en* dictionary")
```

You can also specify a custom dictionary to replace or supplement the built-in 
dictionary. The custom dictionary should be a list in the following 
format. 

```{r}
# name of the vector element is the user input value and the vector element is the 
# replacement value corresponding to that name as a lower case string.
custom_dictionary <- c(
  masculino = "man",
  hombre = "man",
  mujer = "woman",
  femenina = "woman"
)

str(custom_dictionary)
```

Custom dictionaries can be used in place of a built-in dictionary or can 
supplement the built-in dictionary by providing a vector of vectors to the 
dictionary argument. Where the lists contain duplicated elements, the last 
version of the duplicated value will be used for recoding. This allows you to 
use the built-in dictionary but change the coding of one or more responses from 
that dictionary. Here the addition of Spanish terms allows for recoding of 11 
previously uncoded responses.

```{r}

sample %>% 
  mutate(recoded_gender = recode_gender(gender = Gender, 
                                        dictionary = c(fewlevels_en, 
                                                       custom_dictionary))) %>% 
  group_by(recoded_gender) %>% 
  summarise(count = n()) %>% 
  arrange(-count) %>% 
  knitr::kable(caption = "Summary of gender categories after use of the combined dictionaries")


sample %>% 
  mutate(recoded_gender = recode_gender(gender = Gender, 
                                        dictionary = c(fewlevels_en, 
                                                       custom_dictionary))) %>% 
  filter(is.na(recoded_gender)) %>% 
  knitr::kable(caption = "All responses not classified by the combined dictionaries")
```

### retain_unmatched

The `retain_unmatched` argument is used to determine the handling for recoding of values
not contained in the dictionary. By default, unmatched values are coded as NA. 
`retain_unmatched = TRUE` will fill unmatched responses with the participant provided 
response. 

```{r}
sample %>% 
  mutate(recoded_gender = recode_gender(gender = Gender, 
                                        dictionary = c(fewlevels_en, 
                                                       custom_dictionary),
                                        retain_unmatched = TRUE)) %>% 
  group_by(recoded_gender) %>% 
  summarise(count = n()) %>% 
  arrange(-count) %>% 
  knitr::kable(caption = "Summary of gender categories after use of the combined dictionary and `retain_unmatched = TRUE`")
```

# A disclaimer on handling gender responses

This package attempts to remove typographical errors from free text gender data.
The defaults that we used are specific to our context and the time at which the 
package was developed and your data or context may be different. 

We offer two built-in dictionaries, manylevels_en and fewlevels_en. Both are necessarily 
opinionated about how gender descriptors collapse into categories. 

However, as these are culturally specific, they may not be suitable for your data. 
In particular the fewlevels_en option makes opinionated choices about some responses that we want to acknowledge are potentially problematic. Specifically,

- In 'fewlevels_en' coding intersex responses are recoded as 'sex and gender
diverse'
- In 'fewlevels_en' responses where people indicate they are trans and
indicate their identified gender are recoded as the identified gender
(e.g. 'Male to Female' is recoded as 'woman'). We wish to acknowledge
that this may not reflect how some individuals would classify
themselves when given these categories and in some contexts may make
systematic errors. The manylevels_en coding dictionary attempts to avoid these
issues as much as possible - however users can provide a custom
dictionary to add to or overwrite our coding decisions if they feel
this is more appropriate. We welcome people to update the built-in dictionary 
where desired responses are missing.
- In both dictionaries, we assume that typographical features such as spacing
are not relevant to recoding the gender response (e.g. we assume that
"genderqueer" and "gender queer" are equivalent). This is unlikely to be 
true for all contexts. 


The 'manylevels_en' coding separates out those who identify as trans
female/male or cis female/male into separate categories it should not
be assumed that all people who describe as male/female are cis, if you
are assessing trans status we recommend a two part question see:

Bauer, Greta & Braimoh, Jessica & Scheim, Ayden & Dharma, Christoffer.
(2017). Transgender-inclusive measures of sex/gender for population surveys:
Mixed-methods evaluation and recommendations. PLoS ONE. 12.

# Contributing to this package

This package is a reflection of cultural context of the package contributors. 
We acknowledge that understandings of gender are bound by both culture and time 
and are continually changing. As such, we welcome issues and pull requests to 
make the package more inclusive, more reflective of current understandings of 
gender inclusive languages and/or suitable for a broader range of cultural 
contexts. We particularly welcome addition of non-English dictionaries or of 
other gender-diverse responses to the manylevels_en and fewlevels_en dictionaries.

# Acknowledgement of Country

We acknowledge the Wurundjeri people of the Kulin Nation as the custodians of 
the land on which this package was developed and pay respects to elders past, 
present and future.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recode_gender.R
\name{recode_gender}
\alias{recode_gender}
\title{recode_gender}
\usage{
recode_gender(
  gender = gender,
  dictionary = gendercoder::manylevels_en,
  retain_unmatched = FALSE
)
}
\arguments{
\item{gender}{a character vector of gender responses for recoding}

\item{dictionary}{a list that the contains gender responses and their
replacement values. A built-in dictionary \code{manylevels_en} is used by
default if an alternative dictionary is not supplied.}

\item{retain_unmatched}{logical indicating if gender responses that are not found in
dictionary should be filled with the uncleaned values during recoding}
}
\value{
a character vector of recoded genders
}
\description{
\code{recode_gender} matches uncleaned gender responses to cleaned list using
an built-in or custom dictionary.
}
\examples{

\dontrun{

df <- data.frame(
  stringsAsFactors = FALSE,
  gender = c("male", "MALE", "mle", "I am male", "femail", "female", "enby"),
  age = c(34L, 37L, 77L, 52L, 68L, 67L, 83L)
)

df \%>\% mutate(recoded_gender = recode_gender(gender,
  dictionary = manylevels_en,
  retain_unmatched = TRUE
))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gendercodeR.R
\docType{data}
\name{fewlevels_en}
\alias{fewlevels_en}
\title{fewlevels_en}
\description{
A English dictionary for the recode_gender function that has fewer levels
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gendercodeR.R
\docType{data}
\name{manylevels_en}
\alias{manylevels_en}
\title{manylevels_en}
\description{
A English dictionary for the recode_gender function that has many levels
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gendercodeR.R
\docType{data}
\name{sample}
\alias{sample}
\title{sample}
\description{
A sample data.frame of free text gender in English for testing and demonstration
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gendercodeR.R
\docType{package}
\name{gendercoder}
\alias{gendercoder}
\title{gendercoder: A Package for Recoding Freetext Gender Data}
\description{
Provides dictionaries and  a function \code{recode_gender}
to allow for easy automatic coding of common variations in free text
responses to the question \bold{"What is your gender?"}
}
