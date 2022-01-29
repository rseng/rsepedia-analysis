RSelenium
================

[![Build Status](https://travis-ci.org/ropensci/RSelenium.svg?branch=master)](https://travis-ci.org/ropensci/RSelenium)
[![codecov](https://codecov.io/gh/ropensci/RSelenium/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/RSelenium)
[![](http://www.r-pkg.org/badges/version/RSelenium)](https://CRAN.R-project.org/package=RSelenium)
![](http://cranlogs.r-pkg.org/badges/RSelenium?color=yellow)
![](http://cranlogs.r-pkg.org/badges/grand-total/RSelenium?color=yellowgreen)


This is a set of R Bindings for Selenium 2.0 Remote WebDriver, which you can download from http://selenium-release.storage.googleapis.com/index.html. This binding will not work with the 1.0 version of Selenium.


## Install 

To install `RSelenium` from CRAN, run:

```R
install.packages("RSelenium")
```

To install the development version from GitHub, run:

```R
# install.packages("devtools")
devtools::install_github("ropensci/RSelenium")
```

To get started using `RSelenium` you can look at the introduction vignette located in `/doc/basics.html` once `RSelenium` is installed or run

```R
vignette("basics", package = "RSelenium")
```

or the basic vignette can be viewed [here](http://docs.ropensci.org/RSelenium/articles/basics.html).

There is a second vignette dealing with running RSelenium on different browsers/OS locally and remotely which can be viewed at [Driving OS/Browsers Local and Remote](http://docs.ropensci.org/RSelenium/articles/saucelabs.html). Finally, you can read all about running a headless browser or running a normal browser on a headless server at [Headless Browsing](http://docs.ropensci.org/RSelenium/articles/headless.html).

### Summary of Vignettes

1. [Basics](http://docs.ropensci.org/RSelenium/articles/basics.html)
1. [Driving OS/Browsers Local and Remote](http://docs.ropensci.org/RSelenium/articles/saucelabs.html)
1. [Testing Shiny Apps](http://docs.ropensci.org/RSelenium/articles/shinytesting.html)
1. [Headless Browsing](http://docs.ropensci.org/RSelenium/articles/headless.html)
1. [Docker](http://docs.ropensci.org/RSelenium/articles/docker.html)
1. [Internet Explorer](http://docs.ropensci.org/RSelenium/articles/internetexplorer.html)
1. [Orange County R Users Group Webinar](http://docs.ropensci.org/RSelenium/articles/webinar.html)


## Test Shiny Apps

Use `RSelenium` to test your Shiny Apps. Read the introductory tutorial [here](http://docs.ropensci.org/RSelenium/articles/shinytesting.html).


## Use [Sauce Labs](https://saucelabs.com/) and [BrowserStack](https://www.browserstack.com/)

### Sauce Labs

```R
user <- "rselenium0"
pass <- "*******************************"
port <- 80
ip <- paste0(user, ':', pass, "@ondemand.saucelabs.com")
browser <- "firefox"
version <- "25"
platform <- "OS X 10.9"
extraCapabilities <- list(
  name = "Test RSelenium",
  username = user,
  accessKey = pass
)

remDr <- remoteDriver$new(
  remoteServerAddr = ip,
  port = port,
  browserName = browser,
  version = version,
  platform = platform,
  extraCapabilities = extraCapabilities
)
```

### BrowserStack

```R
user <- "johnharrison" 
pass <- "*******************"
port <- 80
ip <- paste0(user, ':', pass, "@hub.browserstack.com")
extraCapabilities <- list(
  "browser" = "IE",
  "browser_version" = "7.0",
  "os" = "Windows",
  "os_version" = "XP",
  "browserstack.debug" = "true"
)

remDr <- remoteDriver$new(
  remoteServerAddr = ip,
  port = port,
  extraCapabilities = extraCapabilities
)
```


## Related Work

* [seleniumPipes](https://github.com/johndharrison/seleniumPipes): A lightweight implementation of the w3c webdriver specification. It has been built utilising `xml2`, `httr` and `magrittr` so provides an alternative for users who are familiar with piping.
* [webdriver](https://github.com/rstudio/webdriver): A client for the 'WebDriver API'. It allows driving a (probably headless) web browser, and can be used to test web applications, including `Shiny` apps. In theory it works with any 'WebDriver' implementation, but it was only tested with 'PhantomJS'.
* [rwebdriver](https://github.com/crubba/Rwebdriver): R bindings to the Webdriver API
* [rdom](https://github.com/cpsievert/rdom): Render and parse the DOM from R via phantomjs.


## License

The RSelenium package is licensed under the [AGPLv3](https://www.r-project.org/Licenses/AGPL-3). The help files are licensed under the creative commons attribution, non-commercial, share-alike license [CC-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/).

As a summary, the AGPLv3 license requires, attribution, include copyright and license in copies of the software, state changes if you modify the code, and disclose all source code. Details are in the COPYING file.

---

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# RSelenium 1.7.6
* No functional changes in this version (need to re-submit to CRAN for being archived)
* Fixed typos in vignettes and documentation
* Styled the package with `styler` package following the tidyverse formatting rules

# RSelenium 1.7.5
* Fix switchToWindow issue in fiefox (#143)
* Add a tutorial to allow running RSelenium Tests in Internet Explorer (thanks @zappingseb #193)
* Updated vignettes and documentation

# RSelenium 1.7.4
* `executeScript` now passes a dummy argument
* Defunct `phantom()` function
* Updated unit tests and test environment
* Updated vignettes and documentation

# RSelenium 1.7.3
* Address issue with user/pass credentials being exposed using SauceLabs (thanks @jstockwin #131)
* Cache packages on TRAVIS to reduce runtime (thanks @jstockwin #132)

# RSelenium 1.7.2
* Fixed issue where rsDriver client when failing to open didn't catch error
* Correctly pass the check argument in rsDriver to wdman (thanks @bourdieu #123)

# RSelenium 1.7.1
* Fixed issue where rsDriver was not passing additional arguments via ...
* Fixed issue with rsDriver and Win/Firefox
* serverURL field in remoteDriver class is now set in initialize method

# RSelenium 1.7.0
* Basic vignette update with appendix on using rsDriver
* Print method added for environment returned by rsDriver
* Default PhantomJS version switched to 2.1.1 (2.5.0-beta has old
  version of ghostdriver)

# RSelenium 1.6.6
* phantom is marked as deprecated. To drive PhantomJS via selenium use the
  rsDriver function. To drive directly use wdman::phantomjs

# RSelenium 1.6.5
* checkForServer and startServer are now defunct. rsDriver is marked as a
  dual replacement. Docker is recommended to run a selenium server/browser.

# RSelenium 1.6.4
* Add a rsDriver function to return a Selenium/webdriver server and a 
  browser client.

# RSelenium 1.6.3
* Return a selected value with the selectTag method.

# RSelenium 1.6.1
* Added a selectTag method to the webElement class see #108.
* RSelenium Basics vignette was updated/revised.

# RSelenium 1.6.0
* Moved http package from RCurl to httr see #106.
* Removed dependence on rjson. httr incorporates jsonlite.
* Import base64_decode from openssl.
* Fixed issue with attributes.Selenium not firing error see #109

# RSelenium 1.5.1
* Added a path argument to the remoteDriver class.

# RSelenium 1.4.9
* Fix .DollarNames to correct issues running under recent RStudio version.

# RSelenium 1.4.8
* Added tests for executeScript
* Fixed issue in executeScript/executeAsyncScript with returning nested
    web elements

# RSelenium 1.4.7
* Code tidied up
* statCodes added as an internal data.frame
* tidy up imports. importFrom instead of import

# RSelenium 1.4.6
* Replace calls to cat with message when error

# RSelenium 1.4.5
* Use canonical form for referring to r-project

# RSelenium 1.4.4
* Deprecate startServer and checkForServer (look at processx to manage process)
* Use message rather than print (thanks Dean Attali #88) in checkForServer. Fix typo in startServer (thanks Charles Thompson #85)
* Copy startServer and checkForServer to examples/serverUtils 

# RSelenium 1.4.3
* Moved testing to TRAVIS
* Switch to rjson from RJSONIO as issue with RJSONIO and TRAVIS/covr
* Ported api tests to TRAVIS

# RSelenium 1.4.2
* Add vignette on RSelenium and Docker containers.

# RSelenium 1.4.1
* Add option to pass arguments to JVM in startServer.
* In startServer look for multiple copies of selenium binary in selDIR 
* Make renaming selenium binary optional in checkForServer
* Add option to download beta releases in checkForServer

# RSelenium 1.4.0
* startServer utility function now returns a list of function; getpid returns the process id of the
  started server, the stop function stops the started server using the process id. Thanks to  
  Dan Tenenbaum #67 and Toby Dylan Hocking #72

# RSelenium 1.3.7
* Add fix for multiple/Beta JARS in checkForServer (Thanks Dean Attali #79)
* Update reference for Selenium download (Thanks @mnel)

# RSelenium 1.3.6
* Allow passing of system2 arguments in startServer utility function

# RSelenium 1.3.4
* Fix custom path not being passed correctly to phantom utility function.
* Allowing passing of commandline arguments via utility function startServer.

# RSelenium 1.3.3
* Add utility function makeFirefoxProfile (Thanks Shan Huang #24)
* Fix phantom utility function for OSX (Thanks Carson Sievert #25)

# RSelenium 1.3.2
* Methods now fail with errors if the server returns an error related status code. Summary and Detail of the error are outputted as well as the associated java class.
* Add a phantom utility function to enable driving of phantomjs in webdriver mode independent of Selenium Server.
* Fixed file paths in startServer for windows (Thanks @mnel #22)

# RSelenium 1.3.0
* Add the content from OC-RUG webinar as a vignette.
* Update the Driving OS/Browsers local and remote vignette.

# RSelenium 1.2.5
* Update reference classes to use `@field` and inline docstrings for methods
* Allow partial string matching on the `using` argument of the findElement and findElements method from the remoteDriver class.
* Allow partial string matching on the `using` argument of the findChildElement and findChildElements method from the webElement class.

# RSelenium 1.2.4
* Add getLogtypes() and log(type) methods to remoteDriver class
* Fix getFirefoxProfile so useBase = TRUE works under windows.
* Add additional support for encoding (thanks to Nicola Logrillo issue #16)
* Add file argument to screenshot method in remoteDriver class to allow writing screenshot to file
* Add a getChromeProfile utility function.

# RSelenium 1.2.3
* Add option to display screenshot in viewer panel if using RStudio
## Test environments
* local OS X install, R 3.6.2
* ubuntu 14.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
# Testing `RSelenium`

These tests are converted from the Python tests in the Selenium project. The tests use a set of HTML documents that can be sourced using.

```sh
svn checkout https://github.com/SeleniumHQ/selenium/trunk/common/src/web
```

The tests assume these HTML documents are available and served locally. To serve the files, we use a Docker image [redsadic/docker-http-server](https://hub.docker.com/r/redsadic/docker-http-server/). This image runs the node application http-server exposing the `/public` directory at port 8080. We map the public directory on the container to the web directory above on the Host (We assume the docker commands are issued from the parent folder containing the web directory - a legacy of using the same calls on TRAVIS):

```sh
docker run -d -p 3000:8080 --name http-server -v $(pwd)/web:/public redsadic/docker-http-server&
```

Next, we run a Docker image containing the standalone Selenium server and a chrome browser:

```sh
docker run -d -p 127.0.0.1:4444:4444 -v /dev/shm:/dev/shm --link http-server selenium/standalone-chrome:2.53.1
```

or a debug version with VNC exposed on port 5901 of the host:

```sh
docker run -d -p 5901:5900 -p 127.0.0.1:4444:4444 -v /dev/shm:/dev/shm --link http-server selenium/standalone-chrome-debug:2.53.1
```

The two Docker containers are linked so the Selenium server will be able to access the http server on its port 8080 and referencing the http server as "http-server"

```
http-server:8080/*.html
```

Normally, on the test machine, docker containers are stopped and removed prior to testing:

```sh
docker stop $(docker ps -q)
docker rm $(docker ps -aq)
```

## System variables

For CRAN and TRAVIS compatibility, two environmental variables are looked for: `NOT_CRAN` and `TRAVIS`. If the tests are being run locally with the above setup, you can set these environmental variables `= "true"` or set them in R:

```R
Sys.setenv("NOT_CRAN" = "true")
Sys.setenv("TRAVIS" = "true")
```
