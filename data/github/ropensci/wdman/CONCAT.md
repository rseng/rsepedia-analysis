# wdman

| CRAN version       | Travis build status   | Appveyor build status   | Coverage |
| :-------------: |:-------------:|:-------------:|:-------------:|
| [![](http://www.r-pkg.org/badges/version/wdman)](https://CRAN.R-project.org/package=wdman) | [![Build Status](https://travis-ci.org/ropensci/binman.svg?branch=master)](https://travis-ci.org/ropensci/wdman) | [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/o8q2d6gdm9su5mcy?svg=true)](https://ci.appveyor.com/project/juyeongkim/wdman) | [![codecov](https://codecov.io/gh/ropensci/wdman/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/wdman)|


## Introduction

`wdman` (Webdriver Manager) is an R package that allows the user to manage the downloading/running of third party binaries relating to the webdriver/selenium projects. The package was inspired by a similar node package [webdriver-manager](https://www.npmjs.com/package/webdriver-manager).

The checking/downloading of binaries is handled by the [`binman`](https://github.com/ropensci/binman) package, and the running of the binaries as processes is handled by the [`subprocess`](https://github.com/lbartnik/subprocess) package.

The `wdman` package currently manages the following binaries:

* [Selenium standalone binary](http://selenium-release.storage.googleapis.com/index.html)
* [chromedriver](https://chromedriver.storage.googleapis.com/index.html)
* [PhantomJS binary](http://phantomjs.org/download.html)
* [geckodriver](https://github.com/mozilla/geckodriver/releases)
* [iedriver](https://github.com/SeleniumHQ/selenium/wiki/InternetExplorerDriver)

Associated with the above are five functions to download/manage the binaries:

* `selenium(...)`
* `chrome(...)`
* `phantomjs(...)`
* `gecko(...)`
* `iedriver(...)`


## Installation

You can install `wdman` from GitHub with:

```R
# install.packages("remotes")
remotes::install_github("ropensci/wdman")
```

The package can also be installed from CRAN:

```R
install.packages("wdman")
```


## Example

As an example, we show how one would run the Selenium standalone binary as a process:

### Running the Selenium binary

The binary takes a port argument which defaults to `port = 4567L`. There are a number of optional arguments to use a particular version of the binaries related to browsers selenium may control. By default, the `selenium` function will look to use the latest version of each. 

```R
selServ <- selenium(verbose = FALSE)
selServ$process

## PROCESS 'file50e6163b37b8.sh', running, pid 21289.
```

The `selenium` function returns a list of functions and a handle representing the running process.

The returned `output`, `error` and `log` functions give access to the stdout/stderr pipes and the cumulative stdout/stderr messages respectively.

```R
selServ$log()

## $stderr
## [1] "13:25:51.744 INFO [GridLauncherV3.parse] - Selenium server version: 4.0.0-alpha-2, revision: f148142cf8"         
## [2] "13:25:52.174 INFO [GridLauncherV3.lambda$buildLaunchers$3] - Launching a standalone Selenium Server on port 4567"
## [3] "13:25:54.018 INFO [WebDriverServlet.<init>] - Initialising WebDriverServlet"                                     
## [4] "13:25:54.539 INFO [SeleniumServer.boot] - Selenium Server is up and running on port 4567"                        

## $stdout
## character(0)
```

The `stop` function sends a signal that terminates the process:

```R
selServ$stop()

## TRUE
```

### Available browsers

By default, the `selenium` function includes paths to chromedriver/geckodriver/phantomjs so that the Chrome/Firefox and PhantomJS browsers are available respectively. All versions (chromever, geckover etc) are given as `"latest"`. If the user passes a value of `NULL` for any driver, it will be excluded.

On Windows operating systems, the option to included the Internet Explorer driver is also given. This is set to `iedrver = NULL` so not ran by default. Set it to `iedrver = "latest"` or a specific version string to include it on your Windows.


## Further details

For further details, please see [the package vignette](https://docs.ropensci.org/wdman/articles/basics.html).

---

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# wdman 0.2.5

* Replaced subprocess (archived in CRAN) with [processx](https://processx.r-lib.org/).

# wdman 0.2.4

* Fixed tests for CRAN re-submission.

# wdman 0.2.3

* Fixed issue with checking for JAVA (thanks Michal Libura #15)

# wdman 0.2.2

* Moved unix based systems to write pipes to file.
* Fixed issue with shell escaping paths.

# wdman 0.2.1

* Added a read_pipes internal function for windows drivers
* Fixed an issue with Windows and blocking pipes. A batch file is now ran
  with stdout/stderr piped to file.

# wdman 0.2.0

* Added verbose arguments to the driver functions
* Import semver for parsing semantic versions
* Added basic vignette on operation.
* Set default PhantomJS version to 2.1.1 (2.5.0-beta runs old ghostdriver
  currently).
* Added a check argument to all driver functions.
* Added tests and refactored code.

# wdman 0.1.5

* Use binman::sem_ver for versioning.
* Remove the v in gecko versions.
* Improve logging in selenium function.

# wdman 0.1.4

* Added tests for driver functions.

# wdman 0.1.3

* Added selenium function.
* Added iedriver function.

# wdman 0.1.2

* Added phantomjs function.

# wdman 0.1.1

* Added chrome function.
* Added gecko function.

# wdman 0.1.0

* Added a `NEWS.md` file to track changes to the package.



## Test environments
* local OS X install, R 3.6.2
* ubuntu 14.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
