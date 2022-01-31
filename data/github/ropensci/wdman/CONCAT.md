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
---
title: "Basics"
author: "John D Harrison"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 3
vignette: >
  %\VignetteIndexEntry{Basics}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The goal of this vignette is to describe the basic functionality of the `wdman` package.


## Introduction

`wdman` (Webdriver Manager) is an R package that allows the user to manage the downloading/running of third party binaries relating to the webdriver/selenium projects. The package was inspired by a similar node package [webdriver-manager](https://www.npmjs.com/package/webdriver-manager).

The checking/downloading of binaries is handled by the [binman](https://github.com/ropensci/binman) package and the running of the binaries as processes is handled by the [processx](https://processx.r-lib.org/) package.

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

The driver functions take a number of common arguments (verbose, check, retcommand) which we describe: 

### Verbosity

Each of the driver functions has a `verbose` argument which controls message output to the user. If `verbose = TRUE` then messages are relayed to the user to inform them when drivers are checked/downloaded/ran. The default value for the driver functions is `TRUE`.

```R
selServ <- selenium(verbose = TRUE)

## checking Selenium Server versions:

## BEGIN: PREDOWNLOAD

## BEGIN: DOWNLOAD

## BEGIN: POSTDOWNLOAD

## checking chromedriver versions:

## BEGIN: PREDOWNLOAD

## BEGIN: DOWNLOAD

## BEGIN: POSTDOWNLOAD

## checking geckodriver versions:

## BEGIN: PREDOWNLOAD

## BEGIN: DOWNLOAD

## BEGIN: POSTDOWNLOAD

## checking phantomjs versions:

## BEGIN: PREDOWNLOAD

## BEGIN: DOWNLOAD

## BEGIN: POSTDOWNLOAD

selServ$stop()

## TRUE
```

versus's

```R
selServ <- selenium(verbose = FALSE)
selServ$stop()

## TRUE
```

### Check for updates 

Each driver function has a `check` argument. If `check= TRUE` the function will liaise with the driver repository for any updates. If new
driver versions are available these will be downloaded. The [binman](https://github.com/ropensci/binman) package is used for this purpose.

### Command line output

For diagnostic purposes each driver function has a `retcommand` argument. If `retcommand = TRUE` the command that would have been launched as a process is instead returned as a string. As an example:

```R
selCommand <- selenium(retcommand = TRUE, verbose = FALSE, check = FALSE)
selCommand

## [1] "/usr/bin/java -Dwebdriver.chrome.driver='/Users/jkim/Library/Application Support/binman_chromedriver/mac64/80.0.3987.16/chromedriver' -Dwebdriver.gecko.driver='/Users/jkim/Library/Application Support/binman_geckodriver/macos/0.26.0/geckodriver' -Dphantomjs.binary.path='/Users/jkim/Library/Application Support/binman_phantomjs/macosx/2.1.1/phantomjs-2.1.1-macosx/bin/phantomjs' -jar '/Users/jkim/Library/Application Support/binman_seleniumserver/generic/4.0.0-alpha-2/selenium-server-standalone-4.0.0-alpha-2.jar' -port 4567"
```

```R
chromeCommand <- chrome(retcommand = TRUE, verbose = FALSE, check = FALSE)
chromeCommand

## [1] "/Users/jkim/Library/Application Support/binman_chromedriver/mac64/80.0.3987.16/chromedriver --port=4567 --url-base=wd/hub --verbose"
```


## Selenium Standalone

The `selenium` function manages the Selenium Standalone binary. It can check for updates at http://selenium-release.storage.googleapis.com/index.html and run the resulting binaries as processes.

### Running the Selenium binary

The binary takes a port argument which defaults to `port = 4567L`. There are a number of optional arguments to use a particular version of the binaries related to browsers selenium may control. By default the `selenium` function will look to use the latest version of each.

```R
selServ <- selenium(verbose = FALSE, check = FALSE)
selServ$process

## PROCESS 'file50e6163b37b8.sh', running, pid 21289.
```

The selenium function returns a list of functions and a handle representing the running process.

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

By default the `selenium` function includes paths to chromedriver/geckodriver/ phantomjs so that the Chrome/Firefox and PhantomJS browsers are available respectively. All versions (chromever, geckover etc) are given as "latest". If the user passes a value of NULL for any driver, it will be excluded.

On Windows operating systems the option to included the Internet Explorer driver is also given. This is set to `iedrver = NULL` so not ran by default. Set it to `iedrver = "latest"` or a specific version string to include it on your Windows.


## Chrome Driver

The `chrome` function manages the Chrome Driver binary. It can check for updates at https://chromedriver.storage.googleapis.com/index.html
and run the resulting binaries as processes.

**The `chrome` function runs the Chrome Driver binary as a standalone process. It takes a default `port` argument `port = 4567L`. Users can then connect directly to the chrome driver to drive a chrome browser.**

Similarly to the `selenium` function, the `chrome` function returns a list of four functions and a handle to the underlying running process.

```R
cDrv <- chrome(verbose = FALSE, check = FALSE)
cDrv$process

## PROCESS 'file534c4e940dd8.sh', running, pid 21386.
```

```R
cDrv$log()

## $stderr
## character(0)

## $stdout
## [1] "Starting ChromeDriver 80.0.3987.16 (320f6526c1632ad4f205ebce69b99a062ed78647-refs/branch-heads/3987@{#185}) on port 4567"
## [2] "Only local connections are allowed."                                                                                     
## [3] "Please protect ports used by ChromeDriver and related test frameworks to prevent access by malicious code."  
```

```R
cDrv$stop()

## TRUE
```


## PhantomJS

The `phantomjs` function manages the PhantomJS binary. It can check for updates at https://bitbucket.org/ariya/phantomjs/downloads and run the resulting binaries as processes.

**The `phantomjs` function runs the PhantomJS binary as a standalone process in webdriver mode. It takes a default `port` argument `port = 4567L`. Users can then connect directly to the "ghostdriver" to drive a PhantomJS browser. Currently the default `version` is set to `version = "2.1.1"`. At the time of writing `2.5.0-beta` has been released. It currently does not have an up-to-date version of ghostdriver associated with it. For this reason it will be unstable/unpredictable to use it in webdriver mode. **

Similarly to the `selenium` function, the `phantomjs` function returns a list of four functions and a handle to the underlying running process.

```R
pjsDrv <- phantomjs(verbose = FALSE, check = FALSE)
pjsDrv$process

## PROCESS 'file5394b74d790.sh', running, pid 21443.
```

```R
pjsDrv$log()

## $stderr
## character(0)

## $stdout
## [1] "[INFO  - 2020-01-31T21:32:04.538Z] GhostDriver - Main - running on port 4567"
```

```R
pjsDrv$stop()

## TRUE
```


## Gecko Driver

The `gecko` function manages the Gecko Driver binary. It can check for updates at https://github.com/mozilla/geckodriver/releases and run the resulting binaries as processes.

**The `gecko` function runs the Gecko Driver binary as a standalone process. It takes a default `port` argument `port = 4567L`. Users can then connect directly to the gecko driver to drive a firefox browser. Currently the default `version` is set to `version = "2.1.1"`.**

**A very IMPORTANT point to note is that geckodriver implements the W3C webdriver protocol which as at the time of writing is not finalised. Currently packages such as RSelenium implement the JSONwireprotocol which whilst similar expects different return from the underlying driver.**

**The geckodriver implementation like the W3C webdriver specification is incomplete at this point in time.**

Similarly to the `selenium` function, the `gecko` function returns a list of four functions and a handle to the underlying running process.

```R
gDrv <- gecko(verbose = FALSE, check = FALSE)
gDrv$process

## PROCESS 'file53946017eccb.sh', running, pid 21458.
```

```R
gDrv$log()

## $stderr
## character(0)

## $stdout
## character(0)
```

```R
gDrv$stop()

## TRUE
```


## IE Driver

The `iedriver` function manages the Internet Explorer Driver binary. It can check for updates at http://selenium-release.storage.googleapis.com/index.html
and run the resulting binaries as processes (the iedriver is distributed currently with the Selenium standalone binary amongst other files). 

**The `chrome` function runs the Chrome Driver binary as a standalone process. It takes a default `port` argument `port = 4567L`. Users can then connect directly to the chrome driver to drive a chrome browser.**

**Please note that additional settings are required to drive an Internet Explorer browser. Security settings and zoom level need to be set correctly in the browser. The author of this document needed to set a registry entry (for ie 11). This is outlined at https://github.com/SeleniumHQ/selenium/wiki/InternetExplorerDriver in the required configuration section.**

Similarly to the `selenium` function the `gecko` function returns a list of four functions and a handle to the underlying running process.

```R
ieDrv <- iedriver(verbose = FALSE, check = FALSE)
ieDrv$process

## Process Handle
## command   : C:\Users\john\AppData\Local\binman\binman_iedriverserver\win64\3.0.0\IEDriverServer.exe /port=4567 /log-level=FATAL /log-file=C:\Users\john\AppData\Local\Temp\RtmpqSdw94\file5247395f2a.txt
## system id : 7484
## state     : running
```

```R
ieDrv$log()

## $stderr
## character(0)
## 
## $stdout
## [1] "Started InternetExplorerDriver server (64-bit)"                                          
## [2] "3.0.0.0"                                                                                 
## [3] "Listening on port 4567"                                                                  
## [4] "Log level is set to FATAL"                                                               
## [5] "Log file is set to C:\\Users\\john\\AppData\\Local\\Temp\\RtmpqSdw94\\file5247395f2a.txt"
## [6] "Only local connections are allowed"
```

```R
ieDrv$stop()

## [1] TRUE
```


## Issues and problems

If you experience issues or problems running one of the drivers/functions, please try running the command in a terminal on your OS initially. You can access the command to run by using the `retcommand` argument in each of the main package functions. If you continue to have problems, consider posting an issue at https://github.com/ropensci/wdman/issues
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chrome.R
\name{chrome}
\alias{chrome}
\title{Start chrome driver}
\usage{
chrome(
  port = 4567L,
  version = "latest",
  path = "wd/hub",
  check = TRUE,
  verbose = TRUE,
  retcommand = FALSE,
  ...
)
}
\arguments{
\item{port}{Port to run on}

\item{version}{what version of chromedriver to run. Default = "latest"
which runs the most recent version. To see other version currently
sourced run binman::list_versions("chromedriver")}

\item{path}{base URL path prefix for commands, e.g. wd/hub}

\item{check}{If TRUE check the versions of chromedriver available. If
new versions are available they will be downloaded.}

\item{verbose}{If TRUE, include status messages (if any)}

\item{retcommand}{If TRUE return only the command that would be passed
to \code{\link[processx]{process}}}

\item{...}{pass additional options to the driver}
}
\value{
Returns a list with named elements \code{process}, \code{output},
    \code{error}, \code{stop}, and \code{log}.
    \code{process} is the object from calling \code{\link[processx]{process}}.
    \code{output} and \code{error} are the functions reading the latest
    messages from "stdout" and "stderr" since the last call whereas \code{log}
    is the function that reads all messages.
    Lastly, \code{stop} call the \code{kill} method in
    \code{\link[processx]{process}} to the kill the \code{process}.
}
\description{
Start chrome driver
}
\examples{
\dontrun{
cDrv <- chrome()
cDrv$output()
cDrv$stop()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iedriver.R
\name{iedriver}
\alias{iedriver}
\title{Start IE driver server}
\usage{
iedriver(
  port = 4567L,
  version = "latest",
  check = TRUE,
  loglevel = c("FATAL", "TRACE", "DEBUG", "INFO", "WARN", "ERROR"),
  verbose = TRUE,
  retcommand = FALSE,
  ...
)
}
\arguments{
\item{port}{Port to run on}

\item{version}{what version of IE driver server to run. Default = "latest"
which runs the most recent version. To see other version currently
sourced run binman::list_versions("iedriverserver")}

\item{check}{If TRUE check the versions of IE driver available. If
new versions are available they will be downloaded.}

\item{loglevel}{Specifies the log level used by the server. Valid values
are: TRACE, DEBUG, INFO, WARN, ERROR, and FATAL. Defaults to FATAL
if not specified.}

\item{verbose}{If TRUE, include status messages (if any)}

\item{retcommand}{If TRUE return only the command that would be passed
to \code{\link[processx]{process}}}

\item{...}{pass additional options to the driver}
}
\value{
Returns a list with named elements \code{process}, \code{output},
    \code{error}, \code{stop}, and \code{log}.
    \code{process} is the object from calling \code{\link[processx]{process}}.
    \code{output} and \code{error} are the functions reading the latest
    messages from "stdout" and "stderr" since the last call whereas \code{log}
    is the function that reads all messages.
    Lastly, \code{stop} call the \code{kill} method in
    \code{\link[processx]{process}} to the kill the \code{process}.
}
\description{
Start IE driver server
}
\examples{
\dontrun{
ieDrv <- iedriver()
ieDrv$output()
ieDrv$stop()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wdman-package.r
\docType{package}
\name{wdman}
\alias{wdman}
\title{wdman.}
\description{
wdman.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phantom.R
\name{phantomjs}
\alias{phantomjs}
\title{Start phantomjs}
\usage{
phantomjs(
  port = 4567L,
  version = "2.1.1",
  check = TRUE,
  loglevel = c("INFO", "ERROR", "WARN", "DEBUG"),
  verbose = TRUE,
  retcommand = FALSE,
  ...
)
}
\arguments{
\item{port}{Port to run on}

\item{version}{what version of phantomjs to run. Default = "2.2.1"
which runs the most recent stable version. To see other version currently
sourced run binman::list_versions("phantomjs")}

\item{check}{If TRUE check the versions of phantomjs available. If
new versions are available they will be downloaded.}

\item{loglevel}{Set phantomjs log level [values: fatal, error,
warn, info, config, debug, trace]}

\item{verbose}{If TRUE, include status messages (if any)}

\item{retcommand}{If TRUE return only the command that would be passed
to \code{\link[processx]{process}}}

\item{...}{pass additional options to the driver}
}
\value{
Returns a list with named elements \code{process}, \code{output},
    \code{error}, \code{stop}, and \code{log}.
    \code{process} is the object from calling \code{\link[processx]{process}}.
    \code{output} and \code{error} are the functions reading the latest
    messages from "stdout" and "stderr" since the last call whereas \code{log}
    is the function that reads all messages.
    Lastly, \code{stop} call the \code{kill} method in
    \code{\link[processx]{process}} to the kill the \code{process}.
}
\description{
Start phantomjs in webdriver mode
}
\examples{
\dontrun{
pjs <- phantomjs()
pjs$output()
pjs$stop()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gecko.R
\name{gecko}
\alias{gecko}
\title{Start gecko driver}
\usage{
gecko(
  port = 4567L,
  version = "latest",
  check = TRUE,
  loglevel = c("info", "fatal", "error", "warn", "config", "debug", "trace"),
  verbose = TRUE,
  retcommand = FALSE,
  ...
)
}
\arguments{
\item{port}{Port to run on}

\item{version}{what version of geckodriver to run. Default = "latest"
which runs the most recent version. To see other version currently
sourced run binman::list_versions("geckodriver")}

\item{check}{If TRUE check the versions of geckodriver available. If
new versions are available they will be downloaded.}

\item{loglevel}{Set Gecko log level [values: fatal, error,
warn, info, config, debug, trace]}

\item{verbose}{If TRUE, include status messages (if any)}

\item{retcommand}{If TRUE return only the command that would be passed
to \code{\link[processx]{process}}}

\item{...}{pass additional options to the driver}
}
\value{
Returns a list with named elements \code{process}, \code{output},
    \code{error}, \code{stop}, and \code{log}.
    \code{process} is the object from calling \code{\link[processx]{process}}.
    \code{output} and \code{error} are the functions reading the latest
    messages from "stdout" and "stderr" since the last call whereas \code{log}
    is the function that reads all messages.
    Lastly, \code{stop} call the \code{kill} method in
    \code{\link[processx]{process}} to the kill the \code{process}.
}
\description{
Start gecko driver
}
\examples{
\dontrun{
gDrv <- gecko()
gDrv$output()
gDrv$stop()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selenium.R
\name{selenium}
\alias{selenium}
\title{Start Selenium Server}
\usage{
selenium(
  port = 4567L,
  version = "latest",
  chromever = "latest",
  geckover = "latest",
  iedrver = NULL,
  phantomver = "2.1.1",
  check = TRUE,
  verbose = TRUE,
  retcommand = FALSE,
  ...
)
}
\arguments{
\item{port}{Port to run on}

\item{version}{what version of Selenium Server to run. Default = "latest"
which runs the most recent version. To see other version currently
sourced run binman::list_versions("seleniumserver")}

\item{chromever}{what version of Chrome driver to run. Default = "latest"
which runs the most recent version. To see other version currently
sourced run binman::list_versions("chromedriver"), A value of NULL
excludes adding the chrome browser to Selenium Server.}

\item{geckover}{what version of Gecko driver to run. Default = "latest"
which runs the most recent version. To see other version currently
sourced run binman::list_versions("geckodriver"), A value of NULL
excludes adding the firefox browser to Selenium Server.}

\item{iedrver}{what version of IEDriverServer to run. Default = "latest"
which runs the most recent version. To see other version currently
sourced run binman::list_versions("iedriverserver"), A value of NULL
excludes adding the internet explorer browser to Selenium Server.
NOTE this functionality is Windows OS only.}

\item{phantomver}{what version of PhantomJS to run. Default = "2.2.1"
which runs the most recent stable version. To see other version
currently
sourced run binman::list_versions("phantomjs"), A value of NULL
excludes adding the PhantomJS headless browser to Selenium Server.}

\item{check}{If TRUE check the versions of selenium available and the
versions of associated drivers (chromever, geckover, phantomver,
iedrver). If new versions are available they will be downloaded.}

\item{verbose}{If TRUE, include status messages (if any)}

\item{retcommand}{If TRUE return only the command that would be passed
to \code{\link[processx]{process}}}

\item{...}{pass additional options to the driver}
}
\value{
Returns a list with named elements \code{process}, \code{output},
    \code{error}, \code{stop}, and \code{log}.
    \code{process} is the object from calling \code{\link[processx]{process}}.
    \code{output} and \code{error} are the functions reading the latest
    messages from "stdout" and "stderr" since the last call whereas \code{log}
    is the function that reads all messages.
    Lastly, \code{stop} call the \code{kill} method in
    \code{\link[processx]{process}} to the kill the \code{process}.
}
\description{
Start Selenium Server
}
\examples{
\dontrun{
selServ <- selenium()
selServ$output()
selServ$stop()
}
}
