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

*Wow, no problems at all. :)**Wow, no problems at all. :)*