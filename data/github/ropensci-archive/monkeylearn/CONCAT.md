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
monkeylearn
===========

  [![Project Status: Abandoned – Initial development has started, but there has not yet been a stable, usable release; the project has been abandoned and the author(s) do not intend on continuing development.](https://www.repostatus.org/badges/latest/abandoned.svg)](https://www.repostatus.org/#abandoned)
[![](https://badges.ropensci.org/45_status.svg)](https://github.com/ropensci/onboarding/issues/45)

This package has been archived.

This R package is an interface to the [MonkeyLearn API](http://docs.monkeylearn.com/article/api-reference/). MonkeyLearn is a Machine Learning platform on the cloud that allows software companies and developers to easily extract actionable data from text. :monkey:
# monkeylearn 0.2.0

* New functions `monkey_classify()` and `monkey_extract()` that:
    * Accept as input both a vector and a dataframe and named column
    * Always return a tibble explicitly relating each input to its classification, allowing for the removal of the MD5 hash
    * Have an `unnest` flag to unnest the output (turn 1 row per input into 1 row per output)
    * Have a `.keep_all` flag to retain other columns if input is a dataframe
    * Coerce `NULL` values and empty vectors returned from MonkeyLearn to `NA`s
    * Include inputs that could not be processed as `NA`s in the output
    * Message the first 20 indices of inputs that are not sent to the API (these now include `NA` and `NULL` values as well as empty strings)
    * Message the currently processing batch

* Bug fixes and improvements to `monkeylearn_classify()` and `monkeylearn_extract()`
    * `monkeylearn_classify()` can now accept `params`
    * Fix to messaging when unable to connect to MonkeyLearn API
    * Default texts per request is set to 200 now (the recommended number), rather than 20
    * Addition of message suggesting that users switch to newer functions

* Implementation of `ratelimitr`. Creation and documentation of two environment variables allowing smarter rate handling when querying the MonkeyLearn API.

* Creation of `pkgdown` website

* Programmatic test coverage to re-use common tests for multiple circumstances.

* Use of a `cowsay` monkey when verbose=TRUE.


# monkeylearn 0.1.3

* Better states the dependency on tibble, it is tibble >= 1.2.

* Better handles blank text in input, outputs an empty tibble and a warning if the request is only blank, and a message if only parts of the request are blank.


# monkeylearn 0.1.2

* Disables HTTP2 for now because of a bug for Windows users. Fix by Jeroen Ooms.

# monkeylearn 0.1.1

* Added a `NEWS.md` file to track changes to the package.



## Test environments
* local x86_64-w64-mingw32/x64 install, R 3.3.1
* Ubuntu 12.04 (on Travis CI), R devel, release and oldrel
* Windows on Appveyor CI (stable, patched and oldrel)

## R CMD check results

0 errors | 0 warnings | 0 note

## Release summary


* New functions `monkey_classify()` and `monkey_extract()` that:
    * Accept as input both a vector and a dataframe and named column
    * Always return a tibble explicitly relating each input to its classification, allowing for the removal of the MD5 hash
    * Have an `unnest` flag to unnest the output (turn 1 row per input into 1 row per output)
    * Have a `.keep_all` flag to retain other columns if input is a dataframe
    * Coerce `NULL` values and empty vectors returned from MonkeyLearn to `NA`s
    * Include inputs that could not be processed as `NA`s in the output
    * Message the first 20 indices of inputs that are not sent to the API (these now include `NA` and `NULL` values as well as empty strings)
    * Message the currently processing batch

* Bug fixes and improvements to `monkeylearn_classify()` and `monkeylearn_extract()`
    * `monkeylearn_classify()` can now accept `params`
    * Fix to messaging when unable to connect to MonkeyLearn API
    * Default texts per request is set to 200 now (the recommended number), rather than 20
    * Addition of message suggesting that users switch to newer functions

* Implementation of `ratelimitr`. Creation and documentation of two environment variables allowing smarter rate handling when querying the MonkeyLearn API.

* Creation of `pkgdown` website

* Programmatic test coverage to re-use common tests for multiple circumstances.

* Use of a `cowsay` monkey when verbose=TRUE.


---
