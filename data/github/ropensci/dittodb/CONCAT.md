# dittodb
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/dittodb)](https://CRAN.R-project.org/package=dittodb)
[![macOS](https://github.com/ropensci/dittodb/workflows/check-macOS/badge.svg)](https://github.com/ropensci/dittodb/actions?workflow=check-macOS)
[![Linux](https://github.com/ropensci/dittodb/workflows/check-linux/badge.svg)](https://github.com/ropensci/dittodb/actions?workflow=check-linux)
[![Windows](https://github.com/ropensci/dittodb/workflows/check-windows/badge.svg)](https://github.com/ropensci/dittodb/actions?workflow=check-windows)
[![Codecov test coverage](https://codecov.io/gh/ropensci/dittodb/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/dittodb?branch=main)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

some new content

{dittodb} is a package that makes testing against databases easy. When writing code that relies on interactions with databases, testing has been difficult without recreating test databases in your continuous integration (aka CI) environment, or resorting to using SQLite databases instead of the database engines you have in production. Both have their downsides: recreating database infrastructure is slow, error prone, and hard to iterate with. Using SQLite works well, right up until you use a feature (like [a full outer join](https://www.sqlite.org/omitted.html)) or has [quirks](https://www.sqlite.org/quirks.html) that might differ from your production database. {dittodb} solves this by recording database interactions, saving them as mocks, and then replaying them seamlessly during testing. This means that if you can get a query from your database, you can record the response and reliably reproduce that response in tests.

{dittodb} is heavily inspired by [{httptest}](https://CRAN.R-project.org/package=httptest), if you've used {httptest} before, you'll find many of the interactions similar.

## A quick example
Say we have a database with some [{nycflights}](https://CRAN.R-project.org/package=nycflights13) data in it and we are writing functions that query this data that we want to test. 

For example, we have the simple function that retrieves one airline:

```r
get_an_airline <- function(con) {
  return(dbGetQuery(con, "SELECT carrier, name FROM airlines LIMIT 1"))
}
```

But we want to make sure that this function returns what we expect. To do this, we first record the response we get from the production database:

## {.tabset}

### RMariaDB
```r
start_db_capturing()

con <- DBI::dbConnect(
  RMariaDB::MariaDB(),
  dbname = "nycflights"
)

get_an_airline(con)
DBI::dbDisconnect(con)

stop_db_capturing()
```

### RPostgres
```r
start_db_capturing()

con <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname = "nycflights"
)

get_an_airline(con)
DBI::dbDisconnect(con)

stop_db_capturing()
```

### RSQLite
```r
start_db_capturing()

con <- DBI::dbConnect(RSQLite::SQLite())

get_an_airline(con)
DBI::dbDisconnect(con)

stop_db_capturing()
```

## {.tabset}

This will run the query from `get_an_airline()`, and save the response in a mock directory and file. Then, when we are testing, we can use the following:


### RMariaDB
```r
with_mock_db({
  con <- DBI::dbConnect(
    RMariaDB::MariaDB(),
    dbname = "nycflights"
  )
  
  test_that("We get one airline", {
    one_airline <- get_an_airline()
    expect_s3_class(one_airline, "data.frame")
    expect_equal(nrow(one_airline), 1)
    expect_equal(one_airline$carrier, "9E")
    expect_equal(one_airline$name, "Endeavor Air Inc.")
  })
})
```

### RPostgres
```r
with_mock_db({
  con <- DBI::dbConnect(
    RPostgres::Postgres(),
    dbname = "nycflights"
  )
  
  test_that("We get one airline", {
    one_airline <- get_an_airline()
    expect_s3_class(one_airline, "data.frame")
    expect_equal(nrow(one_airline), 1)
    expect_equal(one_airline$carrier, "9E")
    expect_equal(one_airline$name, "Endeavor Air Inc.")
  })
})
```

### RSQLite
```r
with_mock_db({
  con <- DBI::dbConnect(RSQLite::SQLite())
  
  test_that("We get one airline", {
    one_airline <- get_an_airline()
    expect_s3_class(one_airline, "data.frame")
    expect_equal(nrow(one_airline), 1)
    expect_equal(one_airline$carrier, "9E")
    expect_equal(one_airline$name, "Endeavor Air Inc.")
  })
})
```

##

All without having to ever set a database up on Travis or GitHub Actions 🎉

## Installation
Currently, {dittodb} is on CRAN (The Comprehensive R Archive Network), so you can install it with `install.packages("dittodb")`. 

### Installing a development version

If you would like to use the development version, you can install from GitHub with: `remotes::install_github("ropensci/dittodb")`

_Note_ You may need to add `@main` at the end if you are using a version of {remotes} prior to 2.2.0. Alternatively, you can use `remotes::install_git()` directly: `remotes::install_git("https://github.com/ropensci/dittodb.git")`

## Setup a package to use {dittodb}
Use the function `dittodb::use_dittodb()` to easily get started using {dittodb}. It will add {dittodb} to `Suggests` in the `DESCRIPTION` file and add `library(dittodb)` to `tests/testthat/helper.R`.

## Development
There is extensive information about developing {dittodb} in the vignette [Developing {dittodb}](https://dittodb.jonkeane.com/articles/developing-dittodb.html, please read that before trying to make changes to {dittodb} or running any of the scripts provided in the `db-setup` directory.

In order to test {dittodb} recording functionality locally or on continuous integration, it is helpful to have databases with test data available. This can be accomplished using the scripts in the `db-setup` directory. By default, {dittodb} does not run any tests that require database infrastructure locally.

To get local databases, the easiest way is to use docker and run either the `db-setup/local-mariadb-docker-setup.sh` or `db-setup/local-postgres-docker-setup.sh` which will pull a docker image and set up a test database with the user and passwords that the {dittodb} tests are expecting (and will stop and remove the docker images if they are present). 

On continuous integration, (using GitHub Actions) these scripts in the `db-setup` directory are used to set up these test databases so we can run integration tests (predominantly in the file `tests/testthat/test-dbi-generic-integration.R`).

## Code of Conduct

Please note that the {dittodb} project is released with a [Contributor Code of Conduct](https://dittodb.jonkeane.com/CODE_OF_CONDUCT). By contributing to this project, you agree to abide by its terms.
# dittodb (development version)

# dittodb 0.1.3
* Minor CRAN update that makes vignette execution conditional when `Suggests` packages are not available.

# dittodb 0.1.2
## New features
* Experimental support for [`expect_sql()`] to check if a specific SQL statement is sent in a test without needing a fixture. Useful for when you only want or need to check that a specific query was sent and you don't need to check any code after that. This feature is experimental, so might change in a subsequent release. 

## Bug fixes and test improvements 
* ODBC connections that only specify a dsn now use the dsn as the path (@klmr, #132). 
* Compatibility for the forthcoming {testthat} 3e.
* Test changes for {dbplyr} (@hadley, #134).
* Internal updates for changes in an upcoming {dbplyr} release.

# dittodb 0.1.1
* Minor CRAN update that makes tests and examples conditional when `Suggests` packages are not available.

# dittodb 0.1.0 
* Initial release with functionality for recording and playing back database fixtures from a number of DBI-based drivers ({RSQLite}, {RPostgres}, {RMariaDB}, {RPostgreSQL})
* Thanks to @maelle for PR#12
* `nycflights13_create_sql()` now always uses {DBI}
* bug fixes to cope with {dbplyr}'s unique table name functions and quoting

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

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
* Focusing on what is best not just for us as individuals, but for the
  overall community

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

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
[jkeane@gmail.com].
All complaints will be reviewed and investigated promptly and fairly.

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

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

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
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# Contributing to {dittodb}

This outlines how to propose a change to {dittodb}. For detailed information 
about developing {dittodb} see the [developing {dittodb} vignette](articles/developing-dittodb.html).


### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.
We love getting help with this!

* **Yes, please:** you edit a roxygen comment in a `.R` file under the `R/` directory.
* **No, thanks:** you edit an `.Rd` file under the `man/` directory.

### Prerequisites

In order to not waste effort and time, before you make a substantial pull request,
you should file an issue and make sure someone from the team agrees that it’s a 
problem. If you’ve found a bug, create an associated issue and illustrate the 
bug with a minimal [reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a new git branch for each pull request (PR).
*  Look at the GitHub Actions build status before and after making changes.
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html),
for documentation.
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.
*  New code should generally follow the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  

### Code of Conduct

Please note that the {dittodb} project is released with a  [Contributor Code of Conduct](CODE_OF_CONDUCT.md). 
By contributing to this project, you agree to abide by its terms.
This submission fixes an Error in the check process when `Suggests` packages are unavailable (i.e. RPostgres).

## Test environments
* local R installation, R 4.0.2
* ubuntu 18.04 (on GitHub actions), R 3.3-4.0
* macOS (on GitHub actions), R 3.6, 4.0, devel
* windows (on GitHub actions), R 4.0

## R CMD check results

0 errors | 0 warnings | 0 note

## revdepcheck results

We checked 1 reverse dependency, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
Please briefly describe your problem and what output you expect. If you have a more general question, you might find a quicker answer on <https://community.rstudio.com/> or <https://stackoverflow.com/>.

Please include a minimal reproducible example (also called a reprex). If you've never heard of a [reprex](https://reprex.tidyverse.org/) before, there's a great overview at <https://www.tidyverse.org/help/#reprex>.

---

### Brief description of the problem

### The kind of database backend you are trying to test

```r
# insert reprex here
```
---
title: 'Dittodb: A Test Environment for Database Requests'
authors:
- affiliation: 1
  name: Jonathan Keane
  orcid: 0000-0001-7087-9776
- affiliation: 2
  name: Mauricio Vargas
  orcid: 0000-0003-1017-7574
date: "8 August 2020"
output: pdf_document
bibliography: REFERENCES.bib
tags:
- R
- SQL
- testing
- reproducibility
affiliations:
- index: 1
  name: Socure, Inc.
- index: 2
  name: Pontifical Catholic University of Chile
---

# Summary

`dittodb::` is an R package that makes testing against databases easy, which is heavily inspired by `httptest::` (@webmockr), and follows the same philosophy to make the interaction similar.

When writing code that relies on interactions with databases, testing has been difficult without recreating test databases in your Continuous integration (CI) environment, or resorting to using SQLite databases instead of the database engines you have in production.

On the one hand, recreating database infrastructure is slow, error prone, and hard to iterate with. On the other, SQLite works well right up until you use a feature (i.e. a full outer join) or has a differenyt syntax your production database. `dittodb::` solves this by recording database interactions, saving them as mocks, and then replaying them seamlessly during testing.

With our software you can get a query from your database, record the response and reliably reproduce that response in your workflow no matter if in the context of an isolated script or in tests created with `testthat::` (@testthat).

To provide a seamless experience between using a real database connection and using the mocked version of the database `dittodb::` uses some features of R that are pretty uncommon.

In order to record fixtures while using a real database connection, we use `::trace` from base R (@base) to add code that inspects the queries, defines unique hashes, and saves the results so that they can be used later. This functionality should generally be used to see what interactions a piece of code to be tested is having with a database and either use or edit and use the fixtures it produces in testing.

When using fixtures (i.e. with a mocked database), we use some internals to mock the `::dbConnect` from `DBI::` (@DBI) and replace the true connection with a special mock connection class from `dittodb::` itself (i.e. DBIMockConnection, though there are specific sub-classes for some drivers, such as DBIMockRPostgresConnection). `dittodb::` relies on standard S4 method dispatch to find the appropriate fixture for queries being run during testing.

`pointblank::` (@pointblank), which tests against public databases, provides a good example that shows the benefits of our software. Running formal tests against it can fail because of internet connection issues or anything not related to the software or the database itself, but thanks to dittodb, all the database related tests in `pointblank::`  are mocked, therefore it prevents timeout errors or similar that would break the workflow in gh-actions and would return a false negative for the check results on gh-actions.

# References
# Making testing databases easy with dittodb

Testing is critical for good package development. R has a robust testing ecosystem, and pacakges like httptest+webmockr that ease HTTP testing. However, it can be difficult to test functionality that connects with external databases. dittodb solves this by providing an easy way to test functions that connect to databases.

Using dittodb means that one doesn't need to maintain dbs on testing infrastructure, can run tests quickly without the latency of connecting to a db, and can develop functions that interact with a db without direct access to that same db.

While developing dittodb I learned a few lessons about testing and package development that I will share during the talk: why is CI important and what makes a good test.


# Video

Hello I'm Jonathan Keane — while I was trained as a linguist — I have been working as a data scientist in a number of capacities for at least the past 5 years.

Early on in my career I fell down the rabbit hole of test driven development and have been using it ever since. Testing has not only saved me from costly mistakes or bugs I never would have caught without it, it has also made it easier to write better code.

This passion for testing led me to developing the package {dittodb}. {dittodb} makes testing database interactions fast, simple, and fun. It uses the DBI interface to make it easy to set up a mocked database driven by static fixtures so tests are quick and can run anywhere. Though there are packages that make testing HTTP-based APIs easy, like {httptest} and {webmockr}, until {dittodb} there wasn't an equivalent package that made testing database interactions simple. 

Through my talk, people will have an understanding of how to test databse interactions using {dittodb} as well as some more general testing-related ideas about what makes a good test and how testing can improve your R code as well your life as a developer.


I love building and contributing to tools that help people interact with data and empower them to be able to analyze and model their data better, faster, and easier.

https://docs.google.com/forms/d/e/1FAIpQLSduWOeixz5qtw5i4IEhOTAGqFOwmDy-2NjYYlFyu2Dk9WJewA/viewform
Reviewer comments:

# @helenmiller16's comments

## Vignette comments
Thank you for pointing out about the special getting started link in pkgdown. I've updated the vignette name and it now appears on the front page of the pkgdown website.

I've also included a bit more information in the nycflights13 data, a link to the package, as well as references to the functions to set the data up in a database (or with SQLite)

## Documentation and examples comments
I've added some more documentation as well as a few more exported functions (with documentation). I also added a bit more detail as well. 

I've also added to the examples and made sure that they are runable. One thing I have kept is many of the `\dontrun{}` wrapping around the examples. The reason for keeping that is that many of the tests will interact with either the filesystem (for recording fixtures) or change some aspects of options (for setting mock paths). Neither of which I wanted to do automatically. I did, however, add to the examples so that they can be copied and pasted and run without anything else.

As for links to httptest, Linking directly in the Rd causes a CMD check error (below) when httptest isn't available, and I didn't think that adding it as a suggest was worth it, I have made sure there are a few links to the CRAN page for httptest, however.

```
N  checking Rd cross-references (3.6s)
   Package unavailable to check Rd xrefs: ‘httptest’
```

## Community guidelines comments
I've add a contributing guide.

## Other comments
Starting and stopping mocks. Both reviewers mentioned wanting an additional way to start and stop mocking besides `with_mock_db()`. I've added `start_mock_db()` and `stop_mock_db()` to turn mocking on and off. This way people can step through their tests if they want (and this addresses some of the other comments as well). 


## Organization 
I've reorganized some of the methods to be more thematic and all of the methods that work with the various `DBI` methods (e.g. `DBI::dbQueries-Results.R` for query and result focused methods).

## Debugging
Thank you for the suggestion, I've added `set_dittodb_debug_level()` and included some documentation about how to use it and what is printed with debugging turned on.

## Use of `.Deprecated()`
I've removed `start_capturing()` entirely, but am using `.Deprecated()` for the deprecation of `.db_mock_paths()`

## `redact_columns()` bug
Thank you for catching the bug in `redact_columns()`. I intended for `grep` to work there and possibly even match multiple columns, but it would not have done that without a fix. I've not fixed that and added to the documentation how to use regex to specify more than one column at a time.

## Developer setup
Developer setup, I've added a whole vignette devoted to how to develop `dittodb` and to set up testing databases. I've also reorganized some of the setup scripts to make it easier to see which ones to use, how to turn on the tests for those when you're ready. Finally, I've added a method for changing the ports that are used by the test databases. 


# @etiennebr's comments

## Installation instructions
I've added a new vignette that describes the process to set up development databases, various configuration options (when, why, and how to turn on or off tests locally). I've also clarified and described more about what the database setup scripts are for (and warned against running ones that are intended for / needed for CI compared to ones that are intended for running locally)

## Improved vignettes
I've added a bit more to the getting started vignette to show an iterative process of the record - test process I use with dittodb. I've also added a developer vignette to help getting started with development.

## Community Guidelines
I've added a contributing guidelines document, thanks for noticing this oversight.

## Automated Tests
I've added more description of what tests are not run by default locally (which makes coverage lower like you noted). To get full coverage one has to run all of the database backends during tests and turn on those tests with environment variables (both to make local development easier when one is not messing with some backends as well as to make sure these tests are not run on CRAN). The CI setup turns these on (and sets up the databases) and the coverage figure comes from there.

Thanks for noticing I wasn't recommending the correct `.sh` scripts in my very brief discussion in the readme. I've expanded that discussion in the new development vignette and I renamed some of those scripts to hopefully be more user-friendly and obvious.


## Workflow
Thank you for this suggestion. I thought about this a lot and one of the things I did try out implementing was this kind of dual-use class/connector that would use mocks if available and then use the database if not. It turned out that it actually ended up getting pretty confusing when it was using which connector under the hood, especially when there were different failures on CI compared to locally. I also wanted to avoid needing to make any operational code changes to start using `dittodb` (with the current architecture, the only time the `dittodb` classes are used or needed is when tests are run).

I did make a number of changes that I think address the main struggles here: there is now a `start_mock_db()` and `stop_mock_db()` that make it easier to step through tests line-by-line without needing to use `browser()` or breakpoints. I've added a bit more to the getting started vignette to show an iterative process of the record - test process I use with dittodb

## wrapper for capturing requests 
Thank you for this suggestion, I've implemented this with `capture_db_requests({})` which can be wrapped around any code that you want to record interactions for.

## Data format
Yes, I've actually thought about this a little bit already and you note the pretty clear tension involved in this. Ideally, the fixtures would be plain text and easy to see if things are changing in gitdiffs. Ideally something like a CSV would be great, and I tried that initially, but there's some type-loss in CSVs that require some kind of sidecar metadata file. I've also looked in to using something like [csvy](https://cran.r-project.org/web/packages/csvy/index.html) which has a metadata header that would work well, but if that is a dependency it also means adding `data.table` (and it's dependencies) as a dependency to `dittodb`. 

I think a better approach might be to replicate some of what csvy is doing inside of `dittodb` so we don't take on another dependency and we get the nice result of having CSV (like) fixture files. Collaborators have also wanted to include the ability to use other kinds of serializers for cases when they want to save a lot of data. This is something I don't necessarily want to fully support out of the box / make super easy because I think like you mention having small, digestable, readable fixtures is important to making good tests, it would be nice of us to provide a hook that someone could use if they wanted to to make their own serializing and deserializing functions as well.

I have been tracking thoughts and ideas here: https://github.com/jonkeane/dittodb/issues/61 I agree that this is something that would be really good to resolve, but since it's not a blocker to functionality right now, think it would be good to add this in as a feature in a later release.

### About lintr complaining about some of the fixtures
Yeah, that's an issue. I've added an issue to investigate the ways to address this.
https://github.com/jonkeane/dittodb/issues/94 Since it's not a blocker, I'm going to save this enhancement for a future release, but thank you for pointing it out and the suggestion!

## I wondered why is `.db_mock_path()` prefixed with a `.`?
I've changed this. It was originally based off of httptest's `.mockPaths()` which
took inspiration (and interface) from `.libPaths()`. At this point I think the lineage is too tenuous to make it useful to have the name the way it is, so I've removed the dot entirely. Thanks for commenting on this, a great example of something only making sense through the developer/creator's eyes that would be mind boggling for most people using it.

## `start_capturing()` and `stop_capturing()`
Thank you for the suggestion to use `.Deprecated()` or remove them. Like you mentioned, this is still in development and they have been soft-deprecated for a while so I just removed them. Though I did use `.Deprecated()` when deprecating `.db_mock_paths()` above.

## What's the relation between rstudio/dbtest and jonkeane/dittodb? Maybe it
doesn't even have to be mentioned?
There isn't any relation, I had originally called `dittodb` `dbtest` (after `httptest`) but @maelle pointed out that `rstudio/dbtest` already existed (though was not on CRAN). There isn't really any other connection and they are looking to solve very different problems. I can add a section about that if you think it's important.

## Fixtures vs. snapshots
I had originally conceived of this as making it easier to understand (I know that beginners to testing can find the multiplicity of terms fixture, mock, stud, harness, etc. confusing). I thought using snapshot in the getting started would be helpful, but I agree with you that it wasn't used particularly well. What I've done is used `snapshot` in the very beginning to describe it abstractly and then specifically name it fixture and use fixture throughout elsewhere. 
## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.2 (2020-06-22) |
|os       |macOS Catalina 10.15.6       |
|system   |x86_64, darwin17.0           |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |America/Chicago              |
|date     |2020-10-07                   |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|dittodb |0.1.1 |0.1.2 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
title: "Recording queries with {dittodb} for travelling"
author: "Mauricio Vargas"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Recording queries with {dittodb} for travelling}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r recording, include = FALSE, eval = FALSE}
library(dplyr)
library(dbplyr)

con_psql <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname = "nycflights",
  host = "127.0.0.1",
  user = getOption("dittodb.test.user"),
  password = getOption("dittodb.test.pw")
)
DBI::dbSendStatement(con_psql, "CREATE DATABASE travelling")
DBI::dbDisconnect(con_psql)


con_psql <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname = "travelling",
  host = "127.0.0.1",
  user = getOption("dittodb.test.user"),
  password = getOption("dittodb.test.pw")
)
nycflights13_create_sql(con_psql)
DBI::dbDisconnect(con_psql)

start_db_capturing(path = "./")
con_psql <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname = "travelling",
  host = "127.0.0.1",
  user = getOption("dittodb.test.user"),
  password = getOption("dittodb.test.pw")
)

tbl(con_psql, "flights") %>%
  filter(!is.na(tailnum)) %>%
  filter(arr_delay >= 180) %>%
  select(tailnum) %>%
  distinct() %>%
  collect()

dbDisconnect(con_psql)
stop_db_capturing()
```

```{r setup, include=FALSE}
library(dittodb)

# set the mockPaths for this vignette
db_mock_paths("travelling")

has_postgres <- check_for_pkg("RPostgres", func = warning)
has_dbplyr <- check_for_pkg("dbplyr", func = warning)
has_dplyr <- check_for_pkg("dplyr", func = warning)
can_eval <- has_postgres & has_dbplyr & has_dplyr

knitr::opts_chunk$set(eval = TRUE, message = FALSE, warning = FALSE)
```

# Scope

The present consists in mocking the connection to a real PostgreSQL server that
contains a database version of the {nycflights13} dataset (among other
databases). See [the {nycflights13} vignette](nycflights.html) for
more information about this database.

This example is for you if you ever wondered how to use scripts that you use
at the office when you are at home or travelling. Or how to continue developing
these scripts while you don't have an internet connection.

Many of us have to use databases that are only accessible from a local network. 
The package {dittodb} provides `with_mock_db()` that wraps the code and makes it
possible to run outside the office (or even with no internet access at all!).

# Recording queries

Suppose we are asked to analyze the flights to only show flights with planes 
that have been delayed at least 3 hours.

One would find all the flights that have been delayed by over 3 hours, and then 
only grab the distinct tail numbers. The only consideration would be to filter 
those flights with missing tail number or those will be treated as a single 
plane.

We could run the following code to get that data with a direct connection to the 
database (i.e. at the office):
```{r, error=TRUE, eval=FALSE}
library(dplyr)
library(dbplyr)

con_psql <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname = "travelling",
  host = "127.0.0.1",
  user = "m.ciccone"
)

tbl(con_psql, "flights") %>%
  filter(!is.na(tailnum)) %>%
  filter(arr_delay >= 180) %>%
  select(tailnum) %>%
  distinct()
```

However, this won't work if we can't connect to our database server. And since 
`postgres.server` is an alias to an IP only accessible from the local network at 
our office, we couldn't run this code and get a result elsewhere. But what if we 
wanted to continue work on this analysis on the train home?

*Important:* This example is using phony authentication. Please never write your 
passwords in scripts, use your `.Rprofile`, an environment variable, or some other 
more secure method instead.

One option would be saving a CSV or TXT file of the data manually, and then
manually reading it in to our R session. But this has a number of drawbacks: we
have to mentally keep track of where each query is from, save it to the right
file, read it in to the right place, etc. We also have to maintain a separate
system or code path for reading in the saved files. {dittodb} can take care of
all of this for us in the background, allowing us to record the results of the
necessary queries, and playing them back when those same queries are called
without a connection to the database.

While we are able to connect to the database (i.e. when we are at the office) we
can save the results returned by queries with code like the following (by
calling `start_db_capturing()` before the connection and the code that executes the
queries and then `stop_db_capturing()` at the end):

```{r, eval=FALSE}
library(dittodb)

start_db_capturing()

con_psql <- DBI::dbConnect(
    RPostgres::Postgres(),
    dbname = "dittodb",
    host = "postgres.server",
    user = "m.ciccone"
  )

flights_delayed <- tbl(con_psql, "flights") %>%
  filter(!is.na(tailnum)) %>%
  filter(arr_delay >= 180) %>%
  select(tailnum) %>%
  distinct() %>%
  collect()

flights_delayed

dbDisconnect(con_psql)

stop_db_capturing()
```

```{r cooking show trick, echo=FALSE, eval=can_eval}
library(dplyr)
library(dbplyr)

# this is the same code that is echoed below, but used here to show output that
# the chunk above would produce if it were able to connect
with_mock_db({
  con_psql <- DBI::dbConnect(
    RPostgres::Postgres(),
    dbname = "travelling",
    host = "127.0.0.1",
    user = "m.ciccone"
  )

  flights_delayed_from_mock <- tbl(con_psql, "flights") %>%
    filter(!is.na(tailnum)) %>%
    filter(arr_delay >= 180) %>%
    select(tailnum) %>%
    distinct() %>%
    collect()

  flights_delayed_from_mock
})

# `dbDisconnect` returns TRUE
TRUE
```

# Reproducing query results

If there was a success capturing one or more queries, then we are able to
replicate the result connected to a different network or even without internet
access:

```{r, eval=can_eval}
with_mock_db({
  con_psql <- DBI::dbConnect(
    RPostgres::Postgres(),
    dbname = "travelling",
    host = "127.0.0.1",
    user = "m.ciccone"
  )

  flights_delayed_from_mock <- tbl(con_psql, "flights") %>%
    filter(!is.na(tailnum)) %>%
    filter(arr_delay >= 180) %>%
    select(tailnum) %>%
    distinct() %>%
    collect()

  flights_delayed_from_mock
})
```

One thing to note is that when using `dbplyr`, we need to be a bit careful that 
we wrap the entire interaction in with the database objects in `with_mock_db` if
we are taking advantage of `dbplyr`'s lazy evaluation (which is by default) and
use `collect()` to return the results when you want them recorded. Because 
`dbplyr` waits until the last possible second to request the data, if you don't 
have a `collect()` call (or a call the will implicitly send the query) there 
won't be a query called, and {dittodb} won't see be able to record the response 
from that query.
---
title: "Getting Started with {dittodb}"
author: "Jonathan Keane"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting Started with {dittodb}}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r, include = FALSE}
library(testthat)
library(dittodb)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

db_mock_paths("dittodb")
```

{dittodb} is designed to make it easy and fun to test functions that interact with a database. It works by looking for mock responses for each query you send while you run your tests and will seamlessly pretend that those mocks were provided by the database connection without needing a connection at all.

```{r setup db, include = FALSE, eval = FALSE}
# setup the DB used in the rest of the vignette

con <- dbConnect(
  RMariaDB::MariaDB(),
  dbname = "nycflights",
  host = "127.0.0.1",
  username = "travis",
  password = ""
)

nycflights13_create_sql(con, schema = "nycflights13")

# record interactions

start_db_capturing()
out <- mean_delays("day")
out <- mean_delays("month")
stop_db_capturing()
```

## Our function `mean_delays()` {.tabset}
To get started, imagine that we are working on a package that queries a database that consists of the [nycflights13 data](nycflights.html). We have the following function which takes a column to aggregate by and returns a dataframe with that column and the mean delay for groups based on the values in the column name given.

### RMariaDB
```{r mean_delays}
library(DBI)

mean_delays <- function(group_col) {
  con <- dbConnect(
    RMariaDB::MariaDB(),
    dbname = "nycflights"
  )
  on.exit(dbDisconnect(con))

  query <- glue::glue(
    "SELECT {group_col}, AVG(arr_delay) as mean_delay from nycflights13.flights ",
    "WHERE arr_delay > 0 GROUP BY {group_col}"
  )

  return(dbGetQuery(con, query))
}
```

### RPostgres
```{r mean_delays_rpostgres, eval = FALSE}
library(DBI)

mean_delays <- function(group_col) {
  con <- dbConnect(
    RPostgres::Postgres(),
    dbname = "nycflights"
  )
  on.exit(dbDisconnect(con))

  query <- glue::glue(
    "SELECT {group_col}, AVG(arr_delay) as mean_delay from nycflights13.flights ",
    "WHERE arr_delay > 0 GROUP BY {group_col}"
  )

  return(dbGetQuery(con, query))
}
```

### RSQLite
```{r mean_delays_rsqlite, eval = FALSE}
library(DBI)

mean_delays <- function(group_col) {
  con <- dbConnect(RSQLite::SQLite())
  on.exit(dbDisconnect(con))

  query <- glue::glue(
    "SELECT {group_col}, AVG(arr_delay) as mean_delay from nycflights13.flights ",
    "WHERE arr_delay > 0 GROUP BY {group_col}"
  )

  return(dbGetQuery(con, query))
}
```

## {.unlisted .unnumbered}

If we give it the column `"month"`, we get the following dataframe:

```{r month, eval = FALSE}
mean_delays("month")
```

```{r cooking_show, echo = FALSE}
with_mock_db(mean_delays("month"))
```

Great, now that we have our function we want to test it to make sure it is operating as expected. Normally, we could write something like:

```{r tests_1, eval = FALSE}
library(testthat)

test_that("mean_delays()", {
  out <- mean_delays("month")
  expect_named(out, c("month", "mean_delay"))
  expect_equal(dim(out), c(12, 2))
})
```

And this works just fine if we only ever run your tests locally, but if we want to [run our tests with a Continuous Integration system](http://r-pkgs.had.co.nz/check.html#check) (and yes, we want to do that!), this won't work without first setting up our production database of flights. For our tests, we don't actually need to connect to the database and get new data (and, in fact, that would make some tests fail erroneously suddenly if the underlying changed). Instead, what we want is to take a snapshot of what happens when running the test code, and then be able to use that snapshot when we run tests later. These snapshots are frequently called fixtures (though you might hear people use other names like stubs or mocks).

# Recording fixtures

We can record fixtures of the database interactions with the commands `start_db_capturing()`, run the functions we want to record, and then stop recording with `stop_db_capturing()`.

```{r recording, eval = FALSE}
start_db_capturing()
out <- mean_delays("month")
stop_db_capturing()
```

This will write a new folder (by default in `./tests/testthat/`) with the name of the database (here: `nycflights`) and then write one file with the name `SELECT-e53189.R` which is the fixture for this example. This `SELECT-*` file contains the data that was received from the database for use in tests.

# `with_mock_db()`

Now that we have a fixture, we can use that fixture by wrapping our call that includes a database interaction with the function `with_mock_db()`. This will look for fixtures and use those.

```{r with_mock_1}
with_mock_db(
  mean_delays("month")
)
```

So, now we can write our tests like:

```{r tests_2}
library(testthat)
library(dittodb)

with_mock_db(
  test_that("mean_delays()", {
    out <- mean_delays("month")
    expect_named(out, c("month", "mean_delay"))
    expect_equal(dim(out), c(12, 2))
  })
)
```

When placed inside of `with_mock_db(...)` a call to `mean_delays("month")` will return what we saved as our fixture _as if it had actually connected to the database_ without needing the database to be installed, reachable, operational, or to exist at all anywhere.

If we wanted to test that a day-based aggregation works, we can! Although we will have to make a new fixture. First we would run the following interactively:

```{r recording_days, eval = FALSE}
start_db_capturing()
out <- mean_delays("day")
stop_db_capturing()
```

This will create a new file (`SELECT-16d120.R`) which contains the response when aggregating by day. dittodb saves each database interaction with a hash of the query that is sent, so that a number of different responses from a database can be saved and the correct one will be used when called inside of `with_mock_db(...)`. So now, we could write our new test with:

```{r tests_day_2}
with_mock_db(
  test_that("mean_delays()", {
    out <- mean_delays("day")
    expect_named(out, c("day", "mean_delay"))
    expect_equal(dim(out), c(31, 2))
  })
)
```

# Getting setup to use {dittodb}
Use the function `dittodb::use_dittodb()` to easily get started using {dittodb}. 
It will add {dittodb} to `Suggests` in the `DESCRIPTION` file and add `library(dittodb)`
to `tests/testthat/helper.R`.

# Things to be careful about

There are a few things to be careful about when using dittodb. 

## When to call `dbConnect()`

Always call `dbConnect()` inside of `with_mock_db(...)`. You can make as many calls as you want to the mock database inside of a `with_mock_db(...)`, but you should always make sure that you connect to the database inside of and not outside of `with_mock_db(...)`. This is because when you "connect" to the mock database, a few variables are set that tell dittodb where to look for mocks. It's less important (though still a good idea) to call `dbDisconnect()` inside of `with_mock_db()`. This is also true when recording fixtures with `start_db_recording()`, you should start the recording and then call `dbConnect()`.

## Query size

Recording fixtures saves the whole query to disk in a relatively inefficient way (from a data storage perspective), so be careful with what you save. And you'll want to not save extremely large results if at all possible. This is also a best-practice for writing tests anyway: you should have mocks that are as minimal as possible to test the functionality you need to. Minimal mocks make it easier to change things that aren't relevant to the test (you don't have to change the way data is represented if it's not important to what you're testing) and it makes your tests run faster.


# Advanced uses

There are a number of advanced features that might be useful. However they take a bit of configuration to use.

## Specify a new path

You can control where mocks are read from (when you're using `with_mock_db(...)`) as well as where they are written to (when using `start_db_capturing()`). To do this, use the function `db_mock_paths()`.

You can see what paths are being used by calling `db_mock_paths()` with no arguments. dittodb will look for mocks in each path starting with the first one. When recording mocks, dittodb always uses the first path that is returned by `db_mock_paths()`.

You can add a new path by calling `db_mock_paths("new/path/here")` which will add the path provided to the top of the list of paths to use.

## Redacting

Sometimes (much? most? of the time!) there is sensitive data in your database that you don't actually want to put into your test fixtures. {dittodb} allows you to specify columns that should always be redacted by specifying them like so: 
```
start_db_capturing(redact_columns = c("sensitive_column", "other_sensitive_column"))
```

This will always redact the columns "sensitive_column" and "other_sensitive_column" every time a query is recorded that includes either. The redactor replaces every value in the column with a standard value (for example "[redacted]" for characters, `9` for numerics, `1988-10-11T17:00:00` for date times) see `redact_columns()` for more information.

## You, too, can write a fixture!

When we use `start_db_recording()` to record fixtures, we are creating what some people call fixtures (though other terms for these abound). These are files that are used during testing to represent and provide some data or state necessary to execute the test. In the case of dittodb, these files contain the data that dittodb uses when it pretends to be a live database. During recording, each query that is sent to the database gets a unique identifier (the first 6 digits of the hash of the query) and when the response is received, that response is saved to a file with the first SQL (Structured Query Language) verb (e.g. `SELECT`), a dash, and the hash using the `dput()` function. This lets you craft a fixture that tests exactly what you need to without having extraneous rows or columns that might not be relevant.

You can save our own responses for queries by getting figuring out the hash (the easiest way to do this now is to write the test that you want to create a fixture for, run it and see the error message that looks something like "Couldn't find the file nycflights/SELECT-16d120.R in any of the mock directories." and use the filename from there.) and then saving the dataframe that you want the test to use with the command `dput(df, file = "nycflights/SELECT-16d120.R", control = c("all", "hexNumeric"))` (if the dataframe you want to save is `df` and we are using the path we saw in the example error message). And you've created your own fixture!

You can also take the approach of recording fixtures and then editing them manually to pare them down. The workflow for that would be something like:

```{r editing, eval = FALSE}
# read in the recorded fixture
df_fixt <- source("nycflights/SELECT-16d120.R", keep.source = FALSE)$value

# filter out anything after february and all days after the 9th of the month
df_fixt <- dplyr::filter(df_fixt, month <= 2 & day < 10)

# save the fixture for use in tests
dput(df_fixt, file = "nycflights/SELECT-16d120.R", control = c("all", "hexNumeric"))
```
---
title: "Developing {dittodb}"
author: "Jonathan Keane"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Developing {dittodb}}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(testthat)
library(dittodb)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

db_mock_paths("dittodb")
```

We welcome contributions from anyone, no matter how small or trivial. Documentation additions or typo fixes are especially welcome. For larger, more-functional changes, either see if there is an issue open on GitHub already, or open one with a proposal of the change you would like to make. Please also see our [code of conduct](../docs/CODE_OF_CONDUCT.html) and [contributing guide](../docs/CONTRIBUTING.html).

Developing {dittodb} requires is a bit more complication than developing other R packages for a few reasons: 

1. setting up all of the databases to fully test recording is complicated (which is in some ways the exact reason {dittodb} exists, so you don't have to do this!)
2. some of the mechanisms that make {dittodb} work aren't commonly used in other R packages.

## Setting up databases
In order to fully test that {dittodb} works, we aim to have full coverage and test as many database backends as possible for both recording and using as a mocked database. To do this on continuous integration (CI for short) can be finicky to get working (and on the CI front, we did it once, so that you can use {dittodb} and you won't have to setup your own database backend just to run tests!). Frankly, even doing this set up locally on a second computer can be a pain! We include in the repository a few scripts that make it (relatively) easy to setup testing database backends, as well as some scripts that we use to setup database backends on GitHub Actions.

### What we test
We currently test against the following database backends with GitHub Actions for CI: 

* Postgres (with drivers: [RPostgres](https://CRAN.R-project.org/package=RPostgres), [RPostgreSQL](https://CRAN.R-project.org/package=RPostgreSQL), and [odbc](https://CRAN.R-project.org/package=odbc))
* MariaDB (with driver: [RMariaDB](https://CRAN.R-project.org/package=RMariaDB))
* SQLite (with driver: [RSQLite](https://CRAN.R-project.org/package=RSQLite))

### How to setup test databases locally
All of these (with the exception of SQLite) are tested in the test file `test-dbi-generic-integration.R`. However, tests for each database are only run if specific environment variables are set that trigger them. The reason for this is so that it is easy to test locally without needing to setup databases, but we are covered by these tests being run on GitHub Actions. If you would like to run these tests locally, you can set the following environment variables and run tests as usual (e.g. with `R CMD check`, `devtools::check()`, `devtools::test()`)

* if `DITTODB_ENABLE_PG_TESTS` is `TRUE`, then Postgres-based tests will be run
* if `DITTODB_ENABLE_MARIA_TESTS` is `TRUE`, then MariaDB-based tests will be run

There are a few scripts included in the `db-setup` folder that are helpful for setting up databases. For local tests, we highly recommend using the docker scripts:

* `db-setup/local-mariadb-docker-setup.sh` which starts (or stops and then starts if it's already running) a docker container, installs MariaDB in that container (running on the default port 3306), and loads the correct test user and test data into the database for running tests.
* `db-setup/local-postgres-docker-setup.sh` which starts (or stops and then starts if it's already running) a docker container, installs Postgres in that container (running on the default port 5432), and loads the correct test user and test data into the database for running tests.

If you've already got databases running on the default ports (3306 for MariaDB and 5432 for Postgres) and you want to use the docker scripts, we recommend that you change the ports that docker is using for any databases you're already running. You can use the `DITTODB_MARIA_TEST_PORT` and `DITTODB_PG_TEST_PORT` environment variables to change which port {dittodb} uses to connect to the test databases. The docker scripts above will use these environment variables to map ports if they are set (and exported) for convenience. One thing to note: during {dittodb} tests, if some database drivers attempt to connect to not-running or on-the-wrong-port database backends, they can segfault instead of erroring with a more informative error. If you see this, the first thing to check is that the port variables are being set correctly and that the database backend is up and running normally.

Both of these utilize a few SQL (Structured Query Language) scripts for their respective backends. These might be useful if you're manually adding the test data into a database you already have running, but if you're using the docker scripts above, you shouldn't need to use them at all. 

* `db-setup/[mariadb|postgres]-reset.sql` creates the database `nycflights` and test users (dropping them if they already exist so they are fresh).
* `db-setup/[mariadb|postgres]-nycflights.sql` creates the necessary tables in the `nycflights` database for use in testing.
* `db-setup/populate-dbs.sh` uses the above scripts to populate the databases on GitHub Actions.

### ☠️ What not to run ☠️
The other scripts (e.g. `db-setup/[mariadb|postgres]-brew.sh` and `db-setup/[mariadb|postgres]-docker-container-only.sh`) are only intended for use on GitHub Actions and should not be run locally. They include commands that will remove files necessary to reset database setups that allow for tests to be run. Running them locally will delete files that you might care about.

## Some of the tricky bits that {dittodb} uses
In order to provide a seamless experience between using a real database connection and using the mocked version of the database {dittodb} uses some features of R that are pretty uncommon. This is not intended to be a comprehensive description of {dittodb}'s architecture, but a few things that are uncommon or a little strange.

### Recording
In order to record fixtures while using a real database connection, we use `base::trace()` to add code that inspects the queries (to define unique hashes) and saves the results so that they can be used later. This tracing only happens when using the `start_db_capturing()` functions and should generally not be used during testing by packages that use {dittodb}. Rather, this functionality should generally be used to see what interactions a piece of code to be tested is having with a database and either use or edit and use the fixtures it produces in testing.

### Using a mocked database
When using fixtures (i.e. with a mocked database), we use some internals to mock the `DBI::dbConnect()` function and replace the true connection with a special mock connection class from {dittodb} (`DBIMockConnection`, though there are specific sub-classes for some drivers, e.g. `DBIMockRPostgresConnection`). Then {dittodb} relies on standard S4 method dispatch to find the appropriate fixture for queries being run during testing. 
---
title: "nycflights13 data"
author: "Mauricio Vargas and Jonathan Keane"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{nycflights13 data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
library(dittodb)

# set the mockPaths for this vignette
db_mock_paths("nycflights13")

knitr::opts_chunk$set(eval = TRUE, message = FALSE, warning = FALSE)
```

{dittodb} uses the [{nycflights13}](https://CRAN.R-project.org/package=nycflights13) dataset for testing and example purposes.

# Exploring {nycflights13}
The {nycflights13} dataset contains airline on-time data for all flights 
departing NYC in 2013. It also includes useful metadata on airlines, airports,
weather, and planes.
    
Have a look to the database schema:

![{nycflights13} relational diagram.](relational-nycflights.svg)

# {nycflights13} test database
{dittodb} comes with a small subset of {nycflights13} to be used in testing and 
examples. To access it, use the convenience function `nycflights_sqlite()` which 
will return an `RSQLite` connection the the `nycflights.sqlite` file included 
with {dittodb}. Alternatively, you can connect to this file with 
`system.file("nycflights.sqlite", package = "dittodb")`.


# Adding {nycflights13} data to a database
{dittodb} has a few functions that make loading {nycflights13} data into a 
database easier. `nycflights13_create_sql(con, schema = "nycflights")` will 
write the {nycflights13} data to the database connect to with `con` and write it 
to the schema `nycflights`.

To quickly set up an SQLite version `nycflights13_create_sql()` will create an 
in-memory SQLite database with the {nycflights13} data.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{set_dittodb_debug_level}
\alias{set_dittodb_debug_level}
\title{Set {dittodb}'s debug level}
\usage{
set_dittodb_debug_level(level)
}
\arguments{
\item{level}{a numeric, the level to set to (e.g. 1)}
}
\value{
the level, invisibly
}
\description{
It can be helpful to see what's going on by increasing {dittodb}'s verbosity
which will show what's going on under the hood (e.g. what queries are being
requested, from where). This sets the option \code{dittodb.debug} to the value
given in the \code{level} argument. The option can be set directly with
\code{options(dittodb.debug = n)} as well.
}
\details{
The \code{level} argument is a numeric, where 0 is the default and (relatively)
silent. The higher the level, the more verbose {dittodb} will be.

Currently, {dittodb} only has one level of debugging (any value 1 or
greater), but more might be used in the future.
}
\examples{
set_dittodb_debug_level(1)
set_dittodb_debug_level(0)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_for_pkg}
\alias{check_for_pkg}
\title{Check if a package is installed}
\usage{
check_for_pkg(package, func = stop)
}
\arguments{
\item{package}{the name of the package to check for}

\item{func}{what should this check call if the package is not installed?
This can be any function, but \code{stop}, \code{warning}, \code{skip}, etc. are likely
candidates (default: \code{stop})}
}
\value{
\code{TRUE} if the package is installed, \code{FALSE} if it is not (invisibly)
}
\description{
Uses \code{requireNamespace()} to check if a package is already installed and
provides options for issuing an error, warning, etc. in case the package is
not installed.
}
\details{
It is only exported for use in examples.
}
\examples{
check_for_pkg("DBI")
check_for_pkg("no-such-package", func = message)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dittodb_debug_level}
\alias{dittodb_debug_level}
\title{Get the dittodb debug level and evaluate if it is above a level}
\usage{
dittodb_debug_level(level)
}
\arguments{
\item{level}{the level to test against (greater than or equal to)}
}
\value{
logical
}
\description{
Get the dittodb debug level and evaluate if it is above a level
}
\examples{
dittodb_debug_level(0)
dittodb_debug_level(2)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/capture-requests.R
\name{get_redactor}
\alias{get_redactor}
\title{Get the current redactor}
\usage{
get_redactor()
}
\value{
the current list of columns to redact
}
\description{
This function should generally not be used, but must be exported for the
query recording function to work properly
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paths.R
\name{make_path}
\alias{make_path}
\title{make a mock path}
\usage{
make_path(path, type, hash)
}
\arguments{
\item{path}{the path to look in}

\item{type}{what type of query is it? (e.g. \code{SELECT}, \code{INSERT})}

\item{hash}{the hash of the query}
}
\value{
a constructed path to a mock
}
\description{
make a mock path
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/driver-specific-connections.R
\docType{class}
\name{driver-specifc-mock-connections}
\alias{driver-specifc-mock-connections}
\alias{DBIMockSQLiteConnection-class}
\alias{DBIMockRPostgreSQLConnection-class}
\alias{DBIMockRPostgresConnection-class}
\alias{DBIMockMariaDBConnection-class}
\title{Driver-specific mock classes}
\description{
Each of the drivers that are supported have their own mock connection class.
They all inherit from \code{DBIMockConnection} as well as their own driver's
connection class. Each is only really available if the corresponding package
is installed.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{db_path_sanitize}
\alias{db_path_sanitize}
\title{Switch illegal characters for legal ones}
\usage{
db_path_sanitize(filename, replacement = "_")
}
\arguments{
\item{filename}{the file or folder to sanitize}

\item{replacement}{what should the illegal character(s) be replaced with?
(default: "_")}
}
\value{
the sanitized string
}
\description{
Inspired by the \href{https://CRAN.R-project.org/package=fs}{fs} package's
\code{path_sanitize} function
}
\examples{
db_path_sanitize('this:string"has?issues')
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paths.R
\name{clean_statement}
\alias{clean_statement}
\title{Clean a statement string}
\usage{
clean_statement(string)
}
\arguments{
\item{string}{an SQL statement to clean}
}
\value{
the SQL statement stripped of extraneous bits
}
\description{
SQL statement strings sometimes have characters and specifications that don't
change the meaning or are determined at query time. To avoid this, before
hashing a statement we clean/strip these from the statement
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nycflights13-sql.R
\name{nycflights13_create_sql}
\alias{nycflights13_create_sql}
\title{Create a standardised database for testing}
\usage{
nycflights13_create_sql(con, schema = "", ...)
}
\arguments{
\item{con}{an SQL connection (i.e a PostgreSQL connection)}

\item{schema}{schema to write the tables ("", or no schema by default)}

\item{...}{additional parameters to connect to a database}
}
\value{
the connection given in \code{con} invisibly, generally called for the
side effects of writing to the database
}
\description{
Using the connection given in \code{con}, create a database including a few tables
from the \href{https://CRAN.R-project.org/package=nycflights13}{\code{nycflights13}} dataset.
}
\examples{
\donttest{
if (check_for_pkg("RSQLite", message)) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

  nycflights13_create_sql(con)

  DBI::dbGetQuery(
    con,
    "SELECT year, month, day, carrier, flight, tailnum FROM flights LIMIT 10"
  )

  DBI::dbDisconnect(con)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connection.R, R/dbExistsTable.R,
%   R/dbListTables-Fields.R, R/dbQueries-Results.R, R/dbMisc.R, R/quote.R,
%   R/transactions.R
\docType{class}
\name{mock-db-methods}
\alias{mock-db-methods}
\alias{DBIMockConnection-class}
\alias{dbDisconnect,DBIMockConnection-method}
\alias{dbExistsTable,DBIMockConnection,character-method}
\alias{dbExistsTable,DBIMockConnection,Id-method}
\alias{dbListTables,DBIMockConnection-method}
\alias{dbListFields,DBIMockConnection,character-method}
\alias{dbListFields,DBIMockConnection,Id-method}
\alias{dbListFields,DBIMockConnection,ANY-method}
\alias{DBIMockResult-class}
\alias{dbSendQuery,DBIMockConnection,character-method}
\alias{dbSendQuery,DBIMockConnection,SQL-method}
\alias{dbSendStatement,DBIMockConnection,character-method}
\alias{dbFetch,DBIMockResult-method}
\alias{fetch,DBIMockResult,ANY-method}
\alias{fetch,DBIMockResult-method}
\alias{fetch,DBIMockResult,missing-method}
\alias{dbClearResult,DBIMockResult-method}
\alias{dbHasCompleted,DBIMockResult-method}
\alias{dbGetQuery,DBIMockRPostgreSQLConnection,character-method}
\alias{dbGetRowsAffected,DBIMockResult-method}
\alias{dbGetInfo,DBIMockConnection-method}
\alias{dbWriteTable,DBIMockConnection,character,data.frame-method}
\alias{dbRemoveTable,DBIMockConnection,character-method}
\alias{dbColumnInfo,DBIMockResult-method}
\alias{dbGetInfo,DBIMockResult-method}
\alias{dbQuoteIdentifier,DBIMockRPostgresConnection,character-method}
\alias{dbQuoteIdentifier,DBIMockRPostgresConnection,SQL-method}
\alias{dbQuoteString,DBIMockRPostgresConnection,character-method}
\alias{dbQuoteString,DBIMockRPostgresConnection,SQL-method}
\alias{dbQuoteString,DBIMockMariaDBConnection,character-method}
\alias{dbQuoteString,DBIMockMariaDBConnection,SQL-method}
\alias{dbBegin,DBIMockConnection-method}
\alias{dbCommit,DBIMockConnection-method}
\alias{dbRollback,DBIMockConnection-method}
\title{Methods for interacting with DB mocks instead of an actual database}
\usage{
\S4method{dbDisconnect}{DBIMockConnection}(conn, ...)

\S4method{dbExistsTable}{DBIMockConnection,character}(conn, name, ...)

\S4method{dbExistsTable}{DBIMockConnection,Id}(conn, name, ...)

\S4method{dbListTables}{DBIMockConnection}(conn, ...)

\S4method{dbListFields}{DBIMockConnection,character}(conn, name, ...)

\S4method{dbListFields}{DBIMockConnection,Id}(conn, name, ...)

\S4method{dbListFields}{DBIMockConnection,ANY}(conn, name, ...)

\S4method{dbSendQuery}{DBIMockConnection,character}(conn, statement, ...)

\S4method{dbSendQuery}{DBIMockConnection,SQL}(conn, statement, ...)

\S4method{dbSendStatement}{DBIMockConnection,character}(conn, statement, ...)

\S4method{dbFetch}{DBIMockResult}(res, n = -1, ...)

\S4method{fetch}{DBIMockResult,ANY}(res, n = -1, ...)

\S4method{fetch}{DBIMockResult,missing}(res, n = -1, ...)

\S4method{dbClearResult}{DBIMockResult}(res, n, ...)

\S4method{dbHasCompleted}{DBIMockResult}(res, ...)

\S4method{dbGetQuery}{DBIMockRPostgreSQLConnection,character}(conn, statement, ...)

\S4method{dbGetRowsAffected}{DBIMockResult}(res, ...)

\S4method{dbGetInfo}{DBIMockConnection}(dbObj, ...)

\S4method{dbWriteTable}{DBIMockConnection,character,data.frame}(conn, name, value, ...)

\S4method{dbRemoveTable}{DBIMockConnection,character}(conn, name, ...)

\S4method{dbColumnInfo}{DBIMockResult}(res, ...)

\S4method{dbGetInfo}{DBIMockResult}(dbObj, ...)

\S4method{dbQuoteIdentifier}{DBIMockRPostgresConnection,character}(conn, x, ...)

\S4method{dbQuoteIdentifier}{DBIMockRPostgresConnection,SQL}(conn, x, ...)

\S4method{dbQuoteString}{DBIMockRPostgresConnection,character}(conn, x, ...)

\S4method{dbQuoteString}{DBIMockRPostgresConnection,SQL}(conn, x, ...)

\S4method{dbQuoteString}{DBIMockMariaDBConnection,character}(conn, x, ...)

\S4method{dbQuoteString}{DBIMockMariaDBConnection,SQL}(conn, x, ...)

\S4method{dbBegin}{DBIMockConnection}(conn, ..., name = NULL)

\S4method{dbCommit}{DBIMockConnection}(conn, ..., name = NULL)

\S4method{dbRollback}{DBIMockConnection}(conn, ..., name = NULL)
}
\arguments{
\item{conn}{a database connection (for dispatch with these methods, it should
be of class \code{DBIMockConnection})}

\item{...}{arguments passed on inside of the methods}

\item{name}{name of the table (for \code{\link{dbListFields}}, \code{\link{dbWriteTable}},
\code{\link{dbRemoveTable}})}

\item{statement}{an SQL statement to execute}

\item{res}{a result object (for dispatch with these methods, it should be of
class \code{DBIMockResult})}

\item{n}{number of results to fetch (ignored)}

\item{dbObj}{a database object (a connection, result, etc.) for use in
\code{\link{dbGetInfo}}}

\item{value}{a value (generally a \code{data.frame}) for use in \code{\link{dbWriteTable}}}

\item{x}{a name to quote (for \code{\link{dbQuoteIdentifier}})}
}
\description{
Various methods (\code{dbSendQuery}, \code{dbFetchQuery}) that are mocks of the
\href{https://CRAN.R-project.org/package=DBI}{DBI} methods of the same name.
Instead of actually interacting with a database, they read in mock responses
and the code proceeds after that. These aren't used directly, but are part of
how {dittodb} works.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expect-sql.R
\name{expect_sql}
\alias{expect_sql}
\title{Detect if a specific SQL statement is sent}
\usage{
expect_sql(object, regexp = NULL, ...)
}
\arguments{
\item{object}{the expression to evaluate}

\item{regexp}{the statement to match}

\item{...}{arguments passed to \code{\link[testthat:expect_error]{testthat::expect_error()}}}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
}
\details{
Sometimes all you need to check is if a specific SQL statement has been sent
and you don't care about retrieving the results.

This works by raising an error that contains the statement that is sent to the
database as well as the location of the result. Currently, \code{expect_sql()} only
works with \code{\link[DBI:dbSendQuery]{DBI::dbSendQuery()}} (and most implementations of \code{\link[DBI:dbGetQuery]{DBI::dbGetQuery()}}
which call \code{\link[DBI:dbSendQuery]{DBI::dbSendQuery()}} internally).

\emph{Note:} this function is experimental and will likely evolve over time. Please
be prepared that new releases might break backwards compatibility.
}
\examples{
if (check_for_pkg("RSQLite", message)) {
  with_mock_db({
    con <- dbConnect(RSQLite::SQLite(), dbname = "not_a_db")

    expect_sql(
      dbGetQuery(con, "SELECT carrier, name FROM airlines LIMIT 3"),
      "SELECT carrier, name FROM airlines LIMIT 3"
    )
  })
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock-paths.R
\name{with_mock_path}
\alias{with_mock_path}
\title{Run the DBI queries in an alternate mock directory}
\usage{
with_mock_path(path, expr, replace = FALSE)
}
\arguments{
\item{path}{the alternate directory}

\item{expr}{the expression to execute}

\item{replace}{logical, should the path replace the current mock paths
(\code{TRUE}) or should they be appended (to the beginning) of the current mock
paths (default, \code{FALSE})}
}
\value{
nothing, called to execute the expression(s) in \code{expr}
}
\description{
When testing with dittodb, wrap your tests in \code{with_mock_path({})} to use the
database fixtures located in other directories. {dittodb} will look for
fixtures in the directory specified by the user, which can be a temporary
or permanent location.
}
\examples{
# Only run if RSQLite and testthat are available
if (check_for_pkg("RSQLite", message) & check_for_pkg("testthat", message)) {
  with_mock_path(
    system.file("nycflight_mocks", package = "dittodb"),
    with_mock_db({
      con <- DBI::dbConnect(
        RSQLite::SQLite(),
        dbname = "nycflights"
      )

      one_airline <- dbGetQuery(
        con,
        "SELECT carrier, name FROM airlines LIMIT 1"
      )
      testthat::test_that("We get one airline", {
        testthat::expect_s3_class(one_airline, "data.frame")
        testthat::expect_equal(nrow(one_airline), 1)
        testthat::expect_equal(one_airline$carrier, "9E")
        testthat::expect_equal(one_airline$name, "Endeavor Air Inc.")
      })
      one_airline
    })
  )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{sanitize_table_id}
\alias{sanitize_table_id}
\title{Sanitize Table id}
\usage{
sanitize_table_id(id, ...)
}
\arguments{
\item{id}{the table identifier (an \code{Id}, a vector of strings, or a string)}

\item{...}{additional arguments (to allow for things like \code{schema_name} that
\code{odbc} uses.)}
}
\value{
the first word in the statement
}
\description{
Tables are identified and specified with a large number of ways across
drivers. For the purposes of {dittodb}, the details are less important since
we almost always just want a flat representation (\emph{ie} for filenames). This
takes the various formats and returns a string with various elements
separated by dots.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock-db.R
\name{mockdb}
\alias{mockdb}
\alias{with_mock_db}
\alias{start_mock_db}
\alias{stop_mock_db}
\title{Run DBI queries against a mocked database}
\usage{
with_mock_db(expr)

start_mock_db()

stop_mock_db()
}
\arguments{
\item{expr}{the expression to execute}
}
\value{
nothing
}
\description{
Wrap a chunk of code in \code{with_mock_db()} to use mocked databases that will
use fixtures instead of connecting to a real database. Alternatively, you can
start and stop using a mocked database with \code{start_mock_db()} and
\code{stop_mock_db()} respectively.to execute the whole thing without needing to
remember to stop the mocking. When testing with {dittodb}, it will look for
fixtures in all entries of \code{\link{db_mock_paths}}.
}
\details{
You only need to use one approach: either use \code{start_mock_db()} to start
using mocks and then \code{stop_mock_db()} to stop or use \code{with_mock_db()} wrapped
around the code you want to execute against the mocked database. You don't
need to (and should not) use both at the same time. Generally
\code{with_mock_db()} is preferred because it is slightly safer and you don't have
to remember to \code{stop_mock_db()} when you're done. However, it is easier to
step through tests interactively using \code{start_mock_db()}/\code{stop_mock_db()}.

Connections should be made after \code{start_mock_db()} if you're using that
function or they should be made inside of \code{with_mock_db()} if you're using
that function because {dittodb} uses the database name (given in \code{dbname} or
\code{Database} argument of \code{\link{dbConnect}} depending on the driver) to separate
different fixtures. For ODBC connections with only a dsn provided, the dsn is
used for this directory.
}
\examples{
# Add the mocks included with dittodb to the db_mock_paths to use them below
db_mock_paths(system.file("nycflight_mocks", package = "dittodb"), last = TRUE)

if (check_for_pkg("RSQLite", message) & check_for_pkg("testthat", message)) {
  # using  `with_mock_db()`
  with_mock_db({
    con <- dbConnect(
      RSQLite::SQLite(),
      dbname = "nycflights"
    )

    testthat::test_that("We get one airline", {
      one_airline <- dbGetQuery(
        con,
        "SELECT carrier, name FROM airlines LIMIT 1"
      )
      testthat::expect_s3_class(one_airline, "data.frame")
      testthat::expect_equal(nrow(one_airline), 1)
      testthat::expect_equal(one_airline$carrier, "9E")
      testthat::expect_equal(one_airline$name, "Endeavor Air Inc.")
    })

    dbDisconnect(con)
  })

  # using `start_mock_db()` and `stop_mock_db()`
  start_mock_db()
    con <- dbConnect(
      RSQLite::SQLite(),
      dbname = "nycflights"
    )

  testthat::test_that("We get one airline", {
    one_airline <- dbGetQuery(
      con,
      "SELECT carrier, name FROM airlines LIMIT 1"
    )
    testthat::expect_s3_class(one_airline, "data.frame")
    testthat::expect_equal(nrow(one_airline), 1)
    testthat::expect_equal(one_airline$carrier, "9E")
    testthat::expect_equal(one_airline$name, "Endeavor Air Inc.")
  })

  dbDisconnect(con)
  stop_mock_db()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_type}
\alias{get_type}
\title{Get the type of an SQL statement}
\usage{
get_type(statement)
}
\arguments{
\item{statement}{the statement to extract the first word from}
}
\value{
the first word in the statement
}
\description{
Get the type of an SQL statement
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_dbname}
\alias{get_dbname}
\title{Get the \code{dbname} from a connection call}
\usage{
get_dbname(dots)
}
\arguments{
\item{dots}{from the argument being passed to the connection}
}
\value{
the name, sanitized if needed
}
\description{
Get the \code{dbname} from a connection call
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/redact.R
\name{redact_columns}
\alias{redact_columns}
\title{Redact columns from a dataframe with the default redactors}
\usage{
redact_columns(data, columns, ignore.case = TRUE, ...)
}
\arguments{
\item{data}{a dataframe to redact}

\item{columns}{character, the columns to redact}

\item{ignore.case}{should case be ignored? (default: \code{TRUE})}

\item{...}{additional options to pass on to \code{grep()} when matching the column
names}
}
\value{
data, with the columns specified in \code{columns} duly redacted
}
\description{
This function redacts the columns specified in \code{columns} in the data given in
\code{data} using {dittodb}'s standard redactors.
}
\details{
The column names given in the \code{columns} argument are treated as regular
expressions, however they always have \code{^} and \code{$} added to the beginning and
end of the strings. So if you would like to match any column that starts with
the string \code{sensitive} (e.g. \code{sensitive_name}, \code{sensitive_date}) you could
use \verb{"sensitive.*} and this would catch all of those columns (though it would
not catch a column called \code{most_sensitive_name}).

The standard redactors replace all values in the column with the following
values based on the columns type:
\itemize{
\item integer -- \code{9L}
\item numeric -- \code{9}
\item character -- \code{"[redacted]"}
\item \code{POSIXct} (date times) -- \code{as.POSIXct("1988-10-11T17:00:00", tz = tzone)}
}
}
\examples{
if (check_for_pkg("nycflights13", message)) {
  small_flights <- head(nycflights13::flights)

  # with no columns specified, redacting does nothing
  redact_columns(small_flights, columns = NULL)

  # integer
  redact_columns(small_flights, columns = c("arr_time"))

  # numeric
  redact_columns(small_flights, columns = c("arr_delay"))

  # characters
  redact_columns(small_flights, columns = c("origin", "dest"))

  # datetiems
  redact_columns(small_flights, columns = c("time_hour"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nycflights13-sql.R
\name{nycflights_sqlite}
\alias{nycflights_sqlite}
\title{An SQLite connection to a subset of nycflights13}
\usage{
nycflights_sqlite()
}
\value{
an RSQLiteConnection
}
\description{
Included with {dittodb} is a small subset of
\href{https://CRAN.R-project.org/package=nycflights13}{\code{nycflights13}}
prepopulated into a \code{sqlite} database.
}
\details{
This database is helpful for getting to know {dittodb} and running example
code. It contains a small subset of the data in nycflights13: namely only the
flights and planes that had a destination of ORD or MDW (the codes for the
two major airports in Chicago) in February of 2013. The airports table has
also been limited to only the New York and Chicago area airports.
}
\examples{
if (check_for_pkg("RSQLite", message)) {
  con <- nycflights_sqlite()

  DBI::dbGetQuery(con, "SELECT flight, tailnum, origin, dest FROM flights LIMIT 10")
  DBI::dbGetQuery(con, "SELECT faa, name, lat, lon, alt, tz FROM airports")

  DBI::dbDisconnect(con)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dittodb-env.R
\docType{data}
\name{.dittodb_env}
\alias{.dittodb_env}
\title{an environment for dittodb storing state}
\format{
An object of class \code{environment} of length 1.
}
\usage{
.dittodb_env
}
\description{
an environment for dittodb storing state
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use-dittodb.R
\name{use_dittodb}
\alias{use_dittodb}
\title{Use {dittodb} in your tests}
\usage{
use_dittodb(path = ".")
}
\arguments{
\item{path}{character path to the package}
}
\value{
Nothing: called for file system side effects.
}
\description{
If you would like to use {dittodb} in your package, and you are already using
\href{https://CRAN.R-project.org/package=testthat}{{testthat}}, use this function to
add {dittodb} to Suggests in the package DESCRIPTION and loads it in
\code{tests/testthat/helper.R}. Call it once when you're setting up a new package
test suite.
}
\details{
This function should be called with the path to your package source as the
\code{path} argument. The function is idempotent: if {dittodb} is already added to
these files, no additional changes will be made.

It will:
\itemize{
\item add {dittodb} to the \code{Suggests} field of the DESCRIPTION file in the
current working directory
\item add \code{library(dittodb)} to the file \code{tests/testthat/helper.R} (creating it
if it doesn't already exist)
}
}
\examples{
\dontrun{
use_dittodb()
use_dittodb("/path/to/package")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/capture-requests.R
\name{capture_requests}
\alias{capture_requests}
\alias{start_db_capturing}
\alias{stop_db_capturing}
\alias{capture_db_requests}
\title{Capture and record database transactions and save them as mocks}
\usage{
start_db_capturing(path, redact_columns = NULL)

stop_db_capturing()

capture_db_requests(expr, path, redact_columns = NULL)
}
\arguments{
\item{path}{the path to record mocks (default if missing: the first path in
\code{db_mock_paths()}.}

\item{redact_columns}{a character vector of columns to redact. Any column
that matches an entry will be redacted with a standard value for the column
type (e.g. characters will be replaced with "[redacted]")}

\item{expr}{an expression to evaluate while capturing requests (for
\code{capture_db_requests()})}
}
\value{
\code{NULL} (invisibily)
}
\description{
When creating database fixtures, it can sometimes be helpful to record
the responses from the database for use in crafting tests.
}
\details{
You can start capturing with \code{start_db_capturing()} and end it with
\code{stop_db_capturing()}. All queries run against a database will be executed like
normal, but their responses will be saved to the mock path given, so that if
you use the same queries later inside of a \code{\link{with_mock_db}} block, the
database functions will return as if they had been run against the database.

Alternatively, you can wrap the code that you are trying to capture in the
function \code{capture_db_requests({...})} this does the same thing as
\code{start_db_capturing()} and \code{stop_db_capturing()} but without needing to
remember to stop the recording.

You can redact certain columns using the \code{redact_columns} argument. This will
replace the values in the column with a generic redacted version. This works
by always passing the data being saved through \code{\link{redact_columns}}.

\emph{note} You should always call \code{\link[DBI:dbConnect]{DBI::dbConnect}} inside of the capturing
block. When you connect to the database, dittodb sets up the mocks for the
specific database you're connecting to when you call \code{\link[DBI:dbConnect]{DBI::dbConnect}}.
}
\examples{
\donttest{
if (check_for_pkg("RSQLite", message)) {
  # Temporary files for examples
  nycflights_path <- tempfile()

  con <- nycflights13_create_sqlite(location = nycflights_path)
  dbDisconnect(con)

  start_db_capturing()
  con <- dbConnect(RSQLite::SQLite(), nycflights_path)

  df_1 <- dbGetQuery(con, "SELECT * FROM airlines LIMIT 1")
  res <- dbSendQuery(con, "SELECT * FROM airlines LIMIT 2")
  df_2 <- dbFetch(res)
  dbClearResult(res)

  dbDisconnect(con)
  stop_db_capturing()

  start_db_capturing(redact_columns = "carrier")
  con <- dbConnect(RSQLite::SQLite(), nycflights_path)

  df_3 <- dbGetQuery(con, "SELECT * FROM airlines LIMIT 3")

  dbDisconnect(con)
  stop_db_capturing()

  with_mock_db({
    con <- dbConnect(RSQLite::SQLite(), nycflights_path)

    # the result from df1 above
    print(dbGetQuery(con, "SELECT * FROM airlines LIMIT 1"))

    # the result from df3 above
    print(dbGetQuery(con, "SELECT * FROM airlines LIMIT 3"))
  })
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nycflights13-sql.R
\name{nycflights13_create_sqlite}
\alias{nycflights13_create_sqlite}
\title{Create an in-memory SQLite database for testing}
\usage{
nycflights13_create_sqlite(location = ":memory:", ...)
}
\arguments{
\item{location}{where to store the database}

\item{...}{additional parameters to connect to a database (most are passed on
to \code{\link{nycflights13_create_sql}})}
}
\value{
RSQLiteConnection
}
\description{
Create an in-memory SQLite database for testing
}
\examples{
\donttest{
if (check_for_pkg("RSQLite", message)) {
  con <- nycflights13_create_sqlite()

  DBI::dbGetQuery(
    con,
    "SELECT year, month, day, carrier, flight, tailnum FROM flights LIMIT 10"
  )

  DBI::dbDisconnect(con)
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paths.R
\name{hash}
\alias{hash}
\title{Make a (short) hash from a string}
\usage{
hash(string, n = 6)
}
\arguments{
\item{string}{the string to hash}

\item{n}{how long should the hash be? (default: 6)}
}
\value{
a hash for the string
}
\description{
Make a (short) hash from a string
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock-paths.R
\name{mockPaths}
\alias{mockPaths}
\alias{db_mock_paths}
\alias{.db_mock_paths}
\title{Set an alternate directory for mock API fixtures}
\usage{
db_mock_paths(new, last = FALSE, replace = FALSE)

.db_mock_paths(new)
}
\arguments{
\item{new}{either a character vector of path(s) to add, or \code{NULL} to reset
to the default.}

\item{last}{a logical, should the new path given be added to the end of the
list of paths? (default: \code{FALSE})}

\item{replace}{logical, should the path replace the current mock paths
(\code{TRUE}) or should they be appended (to the beginning) of the current mock
paths (default, \code{FALSE})}
}
\value{
If \code{new} is omitted, the function returns the current search paths, a
character vector. If \code{new} is provided, the updated value will be returned
invisibly.
}
\description{
By default, \code{with_mock_api} will look for mocks relative to the current
working directory (or the test directory). If you want to look in other
places, you can call \code{db_mock_paths} to add directories to the search path.
}
\details{
It works like \code{\link[base:libPaths]{base::.libPaths()}}: any directories you specify will be added
to the list and searched first. The default directory will be searched last.
Only unique values are kept: if you provide a path that is already found in
\code{db_mock_paths}, the result effectively moves that path to the first
position.

When you are capturing fixtures (e.g. with \code{\link{start_db_capturing}}), the first
path is used as the path to save the fixtures in. For this reason, you may
want to set the \code{last} argument to \code{TRUE} if you want to read from a
directory but don't want to write to it.

For finer-grained control, or to completely override the defaults or any
additions made by calls to \code{db_mock_paths(...)}, you can set the option
"dittodb.mock.paths". If the option "dittodb.mock.paths" is set it will be
used instead of any paths set with \code{db_mock_paths(...)} or even inside of
\code{with_mock_path()}

This function is similar to \code{.mockPaths()} from
\href{https://CRAN.R-project.org/package=httptest}{httptest}

The function \code{.db_mock_paths} is the same as \code{db_mock_paths} although it is
deprecated and should not be used.
}
\examples{
# reset mock paths to default
db_mock_paths(NULL)

identical(db_mock_paths(), c("tests/testthat/", "."))
db_mock_paths("/var/somewhere/else")
identical(db_mock_paths(), c("/var/somewhere/else", "tests/testthat/", "."))
db_mock_paths(NULL)
identical(db_mock_paths(), c("tests/testthat/", "."))
db_mock_paths("/var/somewhere/else", last = TRUE)
identical(db_mock_paths(), c("tests/testthat/", ".", "/var/somewhere/else"))
}
\keyword{internal}
