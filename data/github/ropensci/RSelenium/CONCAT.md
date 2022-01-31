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
---
title: "Testing Shiny Apps"
output:
  html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Testing Shiny Apps with RSelenium}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


## Introduction

The goal of this vignette is to give a basic overview of how one might approach "testing" a shiny app. [Shiny](http://www.rstudio.com/shiny/) is a new package from [RStudio](http://www.rstudio.com/) that makes it dramatically easier to build interactive web applications with R. Shiny Uses a reactive programming model and has built-in widgets derived from the [Bootstrap](http://getbootstrap.com/javascript/) front-end framework. In this vignette we will looking at writing unit tests for a simple shiny wep app. The testing package we will use is [testthat](https://github.com/hadley/testthat) which has a short introduction [here](https://journal.R-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf). I am using `testthat` version 0.8. The version on cran is version 0.7.1 and may give trouble for tests where I manipulate the test environment. You can install 0.8 from github `devtools::install_github("testthat", "hadley")`

Each section will be an introduction to an idea in testing shiny apps with Selenium, and point to more detailed explanation in other vignettes.

## Some Thoughts on Testing

### Why Test?

When faced with testing for the first time the natural reaction is to think what now? what do i test? how much/many tests do I write?

Tests need to do something useful to survive. Automated tests should help the team to make the next move by providing justified confidence a bug has been fixed, confirming refactored code still works as intended, or demonstrating that new features have been successfully implemented. There should be sufficient tests - neither more nor less: more increase the support burden, fewer leave us open to unpleasant surprises in production.

One way to create our tests is to take the view of the user. What does the user want to do?

They want to see this particular graph of a given data set. How do they do that? They select various options and input various choices. From this list of actions we can create an outline of our code for the test. 

For each method, we need to work out how to implement it in code. How could an automated test select the sliderInput bar? Do alternative ways exist? An understanding of HTML, CSS, and JavaScript will help you if you plan to use browser automation tools. All the visible elements of a web application are reflected in the Document Object Model (DOM) in HTML, and they can be addressed in various ways. Some simple examples of interacting with the DOM using `RSelenium` are given in the `Rselenium-basic` vignette.

### Vary the Tests

Having static tests can lead to problems. Introducing variance into the tests can help pick up unexpected errors. This can be achieved by introducing an element of randomness into automatic inputs or randomizing order of selection etc.

### Vary the Browsers/OS

It can help to test against a variety of browsers and operating systems. `RSelenium` can interact with services like [sauceLabs](http://saucelabs.com/). `sauceLabs` allows one to choose the browser or operating system or the version of the selenium server to use. You can test with iOS/Android/Windows/Mac/Linux and browsers like firefox/chrome/ie/opera/safari. This can be very useful to test how your app works on a range of platforms. More detailed information and examples can be seen on the sauceLabs vignette.

### Record the Tests

RSelenium has the ability to take screenshots of the browser at a particular point in time. On failure of a test a screenshot can be useful to understand what happened. If you interface RSelenium with `sauceLabs` you get screenshots and videos automatically. See the sauceLabs vignette for further details.

### Test for Fixes

Lots of bugs are discovered by means other than automated testing - they might be reported by users, for example. Once these bugs are fixed, the fixes must be tested. The tests must establish whether the problem has been fixed and, where practical, show that the root cause has been addressed. Since we want to make sure the bug doesn't resurface unnoticed in future releases, having automated tests for the bug seems sensible.


## The Shiny Test App

### Introduction

The shiny test app is composed of various widgets from the shiny package (0.8.0.99 at time of writing). We have also included the `ggplot2` library as output for one of the charts adapted from a discussion on [stackoverflow](http://stackoverflow.com/questions/11687739/two-legends-based-on-different-datasets-with-ggplot2). The app includes examples of some of the controls included with the `shiny` package namely `selectInput`, `numericInput`, `dateRangeInput` and a `sliderInput`. These controls are used to produce output rendered using `renderPrint`, `renderPlot(base)`  , `renderPlot(ggplot2)` and `renderDataTable`.

The app can be viewed if you have `shiny` installed. 

```R
require(shiny)
runApp(paste0(find.package("RSelenium"), "/apps/shinytestapp"), port = 6012)
```

An image of the app using `RSelenium` on a windows 8.1 machine running firefox 26.0

<h6 align = center>shinytestapp on win 8.1 firefox 26.0</h6>

<img src="https://res.cloudinary.com/johndharrison/image/upload/v1497012339/RSelenium/shinytesting/shinytestapp.png"  title = "shinytestapp on win 8.1 firefox 26.0"  width = '100%'/>

The image was generated using `RSelenium` and the following code.

```R
user <- "rselenium0"
pass <- "***************************"
port <- 80
ip <- paste0(user, ':', pass, "@ondemand.saucelabs.com")
browser <- "firefox"
version <- "26"
platform <- "Windows 8.1"
extraCapabilities <- list(name = "shinytestapp screenshot", username = user, accessKey = pass)

remDr <- remoteDriver$new(remoteServerAddr = ip, port = port, browserName = browser
                          , version = version, platform = platform
                          , extraCapabilities = extraCapabilities)
remDr$open()
remDr$navigate("http://spark.rstudio.com/johnharrison/shinytestapp/")
webElems <- remDr$findElements("css selector", "#ctrlSelect input")
lapply(webElems, function(x){x$clickElement()})
scr <- remDr$screenshot(display = TRUE)
```

### Observations

From the screenshot we retrieved from the remote Driver there are some interesting observations to make. Note that the `selectInput` and `numericInput` boxes are sticking out. This is occurring because the sidePanel is given a bootstrap span of 3. This is however fluid. The resolution on the remote machine is low so the pixel count on the span 3 is also low. On a local machine with high resolution (Nothing amazing just a laptop) we did not observe the `selectInput` and `numericInput` boxes sticking out. 

We could have run with a higher resolution by passing the additional `screen-resolution` parameter to `sauceLabs`. 

```R
extraCapabilities <- list(name = "shinytestapp screenshot", username = user
                          , accessKey = pass, "screen-resolution" = "1280x1024")
```

<h6 align = center>shinytestapp on win 8.1 firefox 26.0 res 1280x1024</h6>

<img src="https://res.cloudinary.com/johndharrison/image/upload/v1497012340/RSelenium/shinytesting/STA-highres.png" title = "shinytestapp on win 8.1 firefox 26.0 res 1280x1024"  width = '100%'/>

We can see things look a bit better but the `data-table` search box is a bit compacted.

### Inputs and Outputs

The app is designed to show testing of the basic shiny components. It is a bit contrived so testing it may not be as natural as testing a live working app. The outputs (charts and tables) are designed to sit side by side if possible with a maximum of 2 on a "row" then drop down to the next "row". We can test to see if this is happening by checking the posistionof elements. We will investigate this later. 

## Basic Tests

### Basic Functionality

The first test we will look at implementing will be basic connection to the app. Typically we would make a request for the page and then observe what status code was returned. Selenium doesn't currently give the html status code of a navigation request so instead we will check if the title of the web page is correct. Our `Shiny Test App` has a title of "Shiny Test App" so we will check for this.

We create a `test/` directory in our `Shiny Test App` folder. The first set of tests will be basic so we create a file `test-basic.r`. In this file we have the following code to start with:

```R
context("basic")

library(RSelenium)
library(testthat)

remDr <- remoteDriver()
remDr$open(silent = TRUE)
appURL <- "http://127.0.0.1:6012"

test_that("can connect to app", {  
  remDr$navigate(appURL)
  appTitle <- remDr$getTitle()[[1]]
  expect_equal(appTitle, "Shiny Test App")  
})

remDr$close()
```

We have a context of "basic" for the tests in this file. The test "can connect to app" simply navigates to the app URL and attempts to get the page title. If the page title is "Shiny Test App" the test is deemed successful. For testing purposes we assume the app is running locally. The easiest way to do this is open a second R session and issue the command:

```R
runApp(paste0(find.package("RSelenium"), "/apps/shinytestapp"), port = 6012)
```

The second R session will listen for connection on port 6012 and return the `Shiny Test App`. If we ran this basic test we would expect the following output:

```R
test_dir(paste0(find.package("RSelenium"), "/apps/shinytestapp/tests/"), filter = 'basic', reporter = "Tap")
```

```
[1] "Connecting to remote server"
1..1
# Context basic 
ok 1 can connect to app 

```

So running the test we observe that we can successfully "connect" to the `Shiny Test App`. What other functionality can we add to our "basic" test context. We can check that the controls and the tabs are present. We can add these tests to our `test-basic.r` file. 

```R
test_that("controls are present", {  
  webElems <- remDr$findElements("css selector", "#ctrlSelect label")
  appCtrlLabels <- sapply(webElems, function(x){x$getElementText()})
  expect_equal(appCtrlLabels[[1]], "Select controls required:")  
  expect_equal(appCtrlLabels[[2]], "selectInput")  
  expect_equal(appCtrlLabels[[3]], "numericInput")  
  expect_equal(appCtrlLabels[[4]], "dateRangeInput")  
  expect_equal(appCtrlLabels[[5]], "sliderInput")  
})

test_that("tabs are present", {  
  webElems <- remDr$findElements("css selector", ".nav a")
  appTabLabels <- sapply(webElems, function(x){x$getElementText()})
  expect_equal(appTabLabels[[1]], "Plots")  
  expect_equal(appTabLabels[[2]], "About")  
})
```

When we rerun our basic test we should hopefully now see that it is checking for the presence of
the controls and the tabs.

```R
test_dir(paste0(find.package("RSelenium"), "/apps/shinytestapp/tests/"), filter = 'basic', reporter = "Tap")
```

```
[1] "Connecting to remote server"
1..8
# Context basic 
ok 1 can connect to app 
ok 2 controls are present 
ok 3 controls are present 
ok 4 controls are present 
ok 5 controls are present 
ok 6 controls are present 
ok 7 tabs are present 
ok 8 tabs are present 
```

That concludes our basic test of the `Shiny Test App` functionality. Next we look at testing the input controls.


## Testing the Controls

Our first test of the controls will be the functioning of the checkbox. We open a new file in the test directory of our `Shiny Test App` and give it the name `test-checkbox.r`. We also give it a context of `controls`.

```R
context("controls")

library(RSelenium)
library(testthat)

remDr <- remoteDriver()
remDr$open(silent = TRUE)
sysDetails <- remDr$getStatus()
browser <- remDr$sessionInfo$browserName
appURL <- "http://127.0.0.1:6012"

test_that("can select/deselect checkbox 1", {  
  remDr$navigate(appURL)
  webElem <- remDr$findElement("css selector", "#ctrlSelect1")
  initState <- webElem$isElementSelected()[[1]]
  # check if we can select/deselect
  if(browser == "internet explorer"){
    webElem$sendKeysToElement(list(key = "space"))
  }else{
    webElem$clickElement()
  }
  changeState <- webElem$isElementSelected()[[1]]
  expect_is(initState, "logical")  
  expect_is(changeState, "logical")  
  expect_false(initState == changeState)  
})

remDr$close()
```

In this case I am informed there maybe issues with `Internet Explorer`. Usually one would select the element for the checkbox and click it. In the case of `Internet Explorer` it maybe necessary to pass a `space` key to the element instead. Otherwise the test is straightforward. We check the initial state of the checkbox. We click the checkbox or send a keypress of space to it. We check the changed state of the checkbox. If the initial state is different to the changed state the test is deemed a success. For good measure we also check that the initial and changed states are of class "logical". We add code for the other 3 checkboxes. We can check our test as follows:

```R
test_dir(paste0(find.package("RSelenium"), "/apps/shinytestapp/tests/"), reporter = "Tap", filter = "checkbox")
```

```
[1] "Connecting to remote server"
1..12
# Context controls 
ok 1 can select/deselect checkbox 1 
ok 2 can select/deselect checkbox 1 
ok 3 can select/deselect checkbox 1 
ok 4 can select/deselect checkbox 2 
ok 5 can select/deselect checkbox 2 
ok 6 can select/deselect checkbox 2 
ok 7 can select/deselect checkbox 3 
ok 8 can select/deselect checkbox 3 
ok 9 can select/deselect checkbox 3 
ok 10 can select/deselect checkbox 4 
ok 11 can select/deselect checkbox 4 
ok 12 can select/deselect checkbox 4 
```

We filter here on "checkbox" to only select this test file to run. If you watch the test running it will filter through the checkbox control checking each checkbox is functioning. The `checkboxGroupInput` drives the required controls which has id `reqcontrols`. Each of these controls is one of the building blocks of shiny and we will add a test for each.

### Testing the selectInput

We write a simple test for the `selectInput`. It tests the options presented and the label of the control. We isolate the code in a separate file `test-selectinput.r` in the test folder of our `Shiny Test App`. It also then selects an element from the options at random. It is tested whether the output changes or not.

```R
test_that("selectInput dataSet correct", {  
  remDr$navigate(appURL)
  webElem <- remDr$findElement("css selector", "#ctrlSelect1")
  initState <- webElem$isElementSelected()[[1]]
  if(!initState){
    # select the checkbox
    if(browser == "internet explorer"){
      webElem$sendKeysToElement(list(key = "space"))
    } else {
      webElem$clickElement()
    }
  }
  
  webElem <- remDr$findElement("css selector", "#reqcontrols #dataset")
  # check the available datasets
  childElems <- webElem$findChildElements("css selector", "[value]")
  appDataSets <- sapply(childElems, function(x){x$getElementAttribute("value")})
  expect_true(all(c("rock", "pressure", "cars") %in% appDataSets))
})

test_that("selectInput label correct", {
  webElem <- remDr$findElement("css selector", "#reqcontrols label[for = 'dataset']")
  expect_output(webElem$getElementText()[[1]], "Choose a dataset:")
}
)


test_that("selectInput selection invokes change", {
  webElem <- remDr$findElement("css selector", "#reqcontrols #dataset")
  childElems <- webElem$findChildElements("css selector", "[value]")
  ceState <- sapply(childElems, function(x){x$isElementSelected()})
  newState <- sample(seq_along(ceState)[!unlist(ceState)], 1)
  
  outElem <- remDr$findElement("css selector", "#summary")
  initOutput <- outElem$getElementText()[[1]]
  
  # change dataset 
  childElems[[newState]]$clickElement()
  outElem <- remDr$findElement("css selector", "#summary")  
  changeOutput <- outElem$getElementText()[[1]]
  
  expect_false(initOutput == changeOutput)
}
)
```

Running the `selectInput` test we get:

```R
test_dir(paste0(find.package("RSelenium"), "/apps/shinytestapp/tests/"), reporter = "Tap", filter = "selectinput")
```

```
[1] "Connecting to remote server"
1..3
# Context controls 
ok 1 selectInput dataSet correct 
ok 2 selectInput label correct 
ok 3 selectInput selection invokes change 
```

Note we set `remDr$setImplicitWaitTimeout(3000)` in this test so that we get a 3 second limit to find an element. 

### Testing the "numericInput"

The ideas behind testing the numericInput are similar to testing the selectInput. We test the label. We then test a random value between the allowable limits of the numericInput and check that the output changes. Finally a character string "test" is sent to the element and the appropriate error message on the output is checked. The final test can be adjusted to suit whatever bespoke error display etc is in your app. The test code is in the tests folder of the `Shiny Test App` in a file named `test-numericinput.r`. Again `remDr$setImplicitWaitTimeout(3000)` is called to give some leeway for element loading. Some commented out code indicates other methods one could deal with checking for element existence. Additional detail on timing races in Selenium can be found [here](http://www.bizalgo.com/2012/01/14/timing-races-selenium-2-implicit-waits-explicit-waits/).

### Testing the "dateRangeInput"

The test on the dateRangeInput compose of two tests. We test the label and we test the two input dates. We choose two random dates from the set of allowable dates. The output is tested for change after the two dates ave been set. `remDr$setImplicitWaitTimeout(3000)` is set in the test to allow for race conditions on elements. The test code is in the tests folder of the `Shiny Test App` in a file named `test-daterangeinput.r`.

### Testing the "sliderInput" 

For the sliderInput we test the label and we test changing the controls. The test code is in the tests folder of the `Shiny Test App` in a file named `test-sliderinput.r`. The label is tested in a similar fashion as the other controls. The second test needs a bit of explaining. There are a number of ways we could interact with the slider control to change its values. Some of the easiest ways would be to execute javascript with `Shiny.onInputChange("range", [2000, 10000])` or 
`Shiny.shinyapp.sendInput({range: [6222, 9333]})`. Both these methods would currently work. The Shiny server side would get the new values however the UI would show no change. The underlying sliderInput control is a `jslider`. Normally one can interact with the `jslider` thru calls similar to `$(".selector").slider("value", p1, p2)` as outlined [here](http://egorkhmelev.github.io/jslider/). We will use mouse movements and the `buttondown` `buttonup` methods of the remoteDriver class. __Note that one may have problems forming the test in this manner, see for example [here](http://stackoverflow.com/questions/19922578/understanding-of-cannot-perform-native-interaction-could-not-load-native-event)__. However it is useful to illustrate mouse and keyboard interactions in `RSelenium`.

We get the attributes of the slider initially. We then get the dimension of the slider

```R
webElem <- remDr$findElement("css selector", "#reqcontrols input#range + .jslider")
sliderDim <- webElem$getElementSize()
```

This gives us the pixel width of the slider as it currently stands. This will be different across machines. We generate some random values for the two slider points and then we calculate roughly how many pixels we need to move the sliders.

```R
remDr$mouseMoveToLocation(webElement = webElems[[x]])
remDr$buttondown()
remDr$mouseMoveToLocation(x = as.integer(pxToMoveSldr[x]), y = -1L)#, webElement = webElems[[x]])
remDr$buttonup()
```

The above code moves to the slider element. Pushes the left button down. Moves the mouse on the x axis in the direction calculated then releases the left mouse button. The output of the related data-table before and after the change is recorded and the test should result in the before and after not being equal.

It is interesting to note that during initial writing of this vignette a new version of firefox 27.0.1 was released. As expected native events did not work under version 2.39 of selenium server and this updated version of firefox. Subsequently our test as formulated above would fail. There is an option to pass a list `rsel.opt` for use with some of the tests. Using this we can set `nativeEvents = FALSE` and the test above will pass again. When your tests fail it is not necessarily bad. This failure indicates a problem with your test setup rather then your app however.

```R
testsel <- test_env()
with(testsel, rsel.opt <- list(nativeEvents = FALSE))
test_dir(paste0(find.package("RSelenium"), "/apps/shinytestapp/tests/"), reporter = "Tap", filter = "slider", env = testsel)
```


## Testing the Output

Finally for this simple example we will look at testing the output. The test code is in the tests folder of the `Shiny Test App` in a file named `test-output.r`. The outputs should line up side by side with a maximum of 2 on a line. We can check the position of the outputs. Our first test will check whether the four outputs line up in a grid. This test will fail on low resolution setups which we will observe latter. We can check the headers on the outputs. The two chart plots are base64 encoded images which we can check in the HTML source. We can check the headers on the outputs. Finally we can check the controls on the datatable.

The first test use the `getElementLocation` method of the `webElement` class to find the location in pixels of the output objects.

```R
webElems <- remDr$findElements("css selector", "#reqplots .span5")
out <- sapply(webElems, function(x){x$getElementLocation()})
```

The 1st and 2nd and the 3rd and 4th objects should share rows. The 1st and 3rd and the 2nd and 4th should share a column. This test will fail as the resolution of the app decreases and the output objects get compacted. The second test checks output labels in a similar fashion to other test.

The third test checks whether the chart output are base 64 encoded png. The final test selects the data-table output and randomly selects a column from carat or price. It then checks whether the ordering functions when the column header is clicked.

Finally running all tests with a "summary" reporter we would hope to get:

```R
test_dir(paste0(find.package("RSelenium"), "/apps/shinytestapp/tests/"))
```

```
basic : [1] "Connecting to remote server"
........
controls : [1] "Connecting to remote server"
............
controls : [1] "Connecting to remote server"
..
controls : [1] "Connecting to remote server"
...
outputs : [1] "Connecting to remote server"
.......
controls : [1] "Connecting to remote server"
...
controls : [1] "Connecting to remote server"
..
```


## Further Tests

* Test across multiple browsers and OS. See the saucelabs testing vignette
* Longitudinal type test. Record access times for various components of your app across time. See the RBMproxy testing vignette.
* Analysis current page load times. See the RBMproxy vignette
---
title: "OCRUG Webinar"
author: "John Harrison"
date: "22 MAY 2014"
output:
  html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{OCRUG Webinar}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Yesterday I gave a webinar on the [Rselenium](https://github.com/johndharrison/RSelenium) package to the OC-RUG. When I published RSelenium on CRAN Ray DiGiacomo the president of OC-RUG was kind enough to invite me to present RSelenium and that webinar was held on the 21st of May. I would like to thank OC-RUG and Ray in particular for the opportunity.

I will use this blog to present the video of the webinar and also for a bit of bookkeeeping. Firstly I would like to return to two of the questions asked in the Q&A.


## Q&A Question 2 and 3

### How Would You Save a Screenshot to File Using RSelenium?

In the seminar I suggested that the Selenium webdriver returned the png in a  base 64 encoded format. I tried to save the output to file using `writeBin` and remarked that I probably needed to decode the output. So here is what I should have done:

```R
b64out <- remDr$screenshot()
writeBin(base64Decode(b64out, "raw"), 'nameoffile.png')
```

### How Would You Use a Chrome Profile?

In the seminar I remarked that the google chrome browser had hundreds of startup options. They are actually termed [command line switches](http://peter.sh/experiments/chromium-command-line-switches). The two I was after were `--user-data-dir` and `--profile-directory`. It would then have been sufficient to use:

```R
args = list(
  paste0('--user-data-dir=',dataDir),
  paste0('--profile-directory=',profileDir)
)
cprof <- list(chromeOptions = list(args = args))
remDr <- remoteDriver(
  browserName = 'chrome',
  extraCapabilities = cprof
)
```

### Added to RSelenium

I added the above functionality to RSelenium. So now the `screenshot` method of the `remoteDriver` class has an optional `file` argument. If `display` is `FALSE` and a `file` argument is given `RSelenium` will now save the `png` to file.

Secondly there is now a `getChromeProfile` utility function in a similar vein to the existing getFirefoxProfile function. The function is documented in the RSelenium package. Its basic usage is as follows:

```R
# example from windows using a profile directory "Profile 1"
dataDir <- "C:\\Users\\john\\AppData\\Local\\Google\\Chrome\\User Data"
cprof <- getChromeProfile(dataDir, "Profile 1")
remDr <- remoteDriver(browserName = "chrome", extraCapabilities = cprof)
```

So it takes two arguments: a users chrome data directory and the profile directory within that data directory from which you want to run.

I will update the version of `RSelenium` on CRAN in the coming weeks. If you would like to run the dev version with added features you can install it with:
  
```R
devtools::install_github("johndharrison/RSelenium")
```

and remember to file any issues you may have. to github so `RSelenium` can improve as a package :).


## The Webinar Screencast

<iframe width="653" height="480" src="https://www.youtube.com/embed/ic65SWRWrKA" frameborder="0" allowfullscreen></iframe>


### SauceLabs Test Results

At the end of the webinar the simple test on the shinytestapp was ran against multiple OS and browsers using Sauce Labs as an external provider for the selenium server and the browsers. These test results can be seen bu navigating to the RSelenium Sauce Labs page [https://saucelabs.com/u/rselenium0](https://saucelabs.com/u/rselenium0) and searching for `OCRUG TEST 1` in the search session name panel.


## The Directors Cut

I have shot an extended version of the webinar with the fabled missing `javascript` example which follows:

<iframe width="653" height="480" src="https://www.youtube.com/embed/zcUGla8EjjY" frameborder="0" allowfullscreen></iframe>


## The Slides and R Scripts

And finally I have posted the slides for the webinar to github. Please find them at [http://johndharrison.github.io/RSOCRUG](http://johndharrison.github.io/RSOCRUG). There are also some R scripts contained in the master branch of the github project: https://github.com/johndharrison/RSOCRUG


## Conclusion

Hopefully that ties up any loose ends. It was great fun preparing the webinar. I enjoyed giving it and I'm gratefully to everyone who attended.
---
title: "Driving OS/Browsers Local and Remote"
output:
  html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Driving locally and remotely with RSelenium}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

The goal of this vignette is to give a basic overview of how one might approach using RSelenium with combinations of operating systems (OS) and browsers both locally and remotely.


## RSelenium with Local Fully-fledged Browsers

### Firefox

The default browser for `RSelenium` is firefox. When a `remoteDriver` class is instantiated using default options for example `remdr <- remoteDriver()` then the browser listed is firefox.

```R
remdr <- remoteDriver()
remDr
```

```
## $remoteServerAddr
## [1] "localhost"
## 
## $port
## [1] 4444
## 
## $browserName
## [1] "firefox"
## 
## $version
## [1] ""
## 
## $platform
## [1] "ANY"
## 
## $javascript
## [1] TRUE
## 
## $autoClose
## [1] FALSE
## 
## $nativeEvents
## [1] TRUE
## 
## $extraCapabilities
## list()
```

Other browsers can be driven using `RSelenium`. We shall split these browsers into three groups. Full-fledged browsers, headless browsers and mobile browsers.

The standalone selenium jar has the ability to drive other full-fledged browsers such as chrome, internet explorer, safari and opera. First we shall look at how to drive chrome using `RSelenium`

### Chrome

Firstly we note that chrome in this instance can be considered as having three parts. There is the browser itself ("chrome"), the language bindings provided by the Selenium project ("the driver") and an executable downloaded from the Chromium project which acts as a bridge between "chrome" and the "driver". This executable is called "chromedriver". We need to have a "chromedriver" running. First we need to locate one. The download directory for chromedriver is currently located at http://chromedriver.storage.googleapis.com/index.html. In this example we shall look at running chrome on a windows platform so we will download the windows chromedriver. The most uptodate version of chromedriver at the time of writing was 2.9. In the notes this supports chrome v31-34. We are running chrome 33 so this is fine.

```
----------ChromeDriver v2.9 (2014-01-31)----------
Supports Chrome v31-34
```

We download the appropriate [file](http://chromedriver.storage.googleapis.com/2.9/chromedriver_win32.zip) for windows and extract the .exe to our Documents folder. The .exe can be placed where the user pleases but it must be in the system path. In this case we placed in the Documents folder namely C:\Users\john\Documents. This directory was added to the system path. 

We assume that a selenium server is also running if not one can be started using `RSelenium::startServer()`. Now we are done. A chrome browser can be controlled as follows:

```R
require(RSelenium)
# RSelenium::startServer() # if needed
remDr <- remoteDriver(browserName = "chrome")
remDr$open()

head(remDr$sessionInfo)
```

```
## $platform
## [1] "WIN8"
## 
## $acceptSslCerts
## [1] TRUE
## 
## $javascriptEnabled
## [1] TRUE
## 
## $browserName
## [1] "chrome"
## 
## $chrome
## $chrome$userDataDir
## [1] "C:\\Users\\john\\AppData\\Local\\Temp\\scoped_dir24584_12002"
## 
## 
## $rotatable
## [1] FALSE
```

### Internet Explorer

Similarly to the chrome browser you do not need to run an installer before using the InternetExplorerDriver, though some configuration is required. The standalone server executable must be downloaded from the Downloads page and placed in your PATH. Again we need to download this executable and place it in our path. At the time of writing the internet explorer .exe is included with the main standalone server [here](http://selenium-release.storage.googleapis.com/index.html). The current release is 2.40. The system I am running is 64 bit so we download the [64bit version](http://selenium-release.storage.googleapis.com/2.40/IEDriverServer_x64_2.40.0.zip). For simplicity we again place this in our Documents directory namely C:\Users\john\Documents. This directory is already in the system path from running the chrome example. If  you want to place the internet explorer .exe in another folder add this folder to your system path. To control internet explorer as a browser is now as simple as:

```R
require(RSelenium)
# RSelenium::startServer() # if needed
remDr <- remoteDriver(browserName = "internet explorer")
remDr$open()

head(remDr$sessionInfo, 7)
```

```
## $platform
## [1] "WINDOWS"
## 
## $javascriptEnabled
## [1] TRUE
## 
## $elementScrollBehavior
## [1] 0
## 
## $ignoreZoomSetting
## [1] FALSE
## 
## $enablePersistentHover
## [1] TRUE
## 
## $ie.ensureCleanSession
## [1] FALSE
## 
## $browserName
## [1] "internet explorer"
```

### Safari

Currently Apple have discontinued developement of safari for windows. The latest version for windows was 5.1.7 available [here](http://filehippo.com/download_safari/). Starting with Selenium 2.30.0, the SafariDriver comes bundled with the Selenium server so nothing other then having safari installed should be required. For the purposes of this vignette I downloaded and installed safari 5.1.7 on a windows 8.1 system. 
Once installed controlling safari was as easy as:

```R
require(RSelenium)
# RSelenium::startServer() # if needed
remDr <- remoteDriver(browserName = "safari")
remDr$open()
head(remDr$sessionInfo)
```

```
## $platform
## [1] "WINDOWS"
## 
## $cssSelectorsEnabled
## [1] TRUE
## 
## $javascriptEnabled
## [1] TRUE
## 
## $secureSsl
## [1] TRUE
## 
## $browserName
## [1] "safari"
## 
## $webdriver.remote.sessionid
## [1] "a18da818-5160-47c4-8e88-7e95605c5cab"
```

### Opera

Opera is currently not supported for versions newer then 12. http://code.google.com/p/selenium/wiki/OperaDriver gives details on the current status. For the purposes of this vignette I downloaded and installed the 64 bit version of 12.16 located [here](http://www.opera.com/download/guide/?os=windows&ver=12.16&local=y). I had some issues getting this to run YMMV however. The first issue was getting the following error message

```R
require(RSelenium)
# RSelenium::startServer() # if needed
remDr <- remoteDriver(browserName = "opera")
remDr$open()

cat(remDr$value$message)
```

```
## Could not find a platform that supports bundled launchers, please set it manually
## Build info: version: '2.40.0', revision: 'fbe29a9', time: '2014-02-19 20:54:28'
## System info: host: 'johnlt', ip: '192.168.58.1', os.name: 'Windows 8', os.arch: 'amd64', os.version: '6.2', java.version: '1.7.0_45'
## Driver info: driver.version: OperaDriver
```

So not surprising win 8.1 was first unveiled and released as a public beta in June 2013 and on July 4, 2013, Opera 12.16 was released being the last current version of opera supported by selenium. OperaLauncherRunner.java located [here](https://github.com/operasoftware/operadriver/blob/master/src/com/opera/core/systems/runner/launcher/OperaLauncherRunner.java) does not currently cater for the WIN8 enum returned by Platform.getCurrent().

The solution is to start the selenium server with additional parameters passed to java (RSelenium::startServer doesn't pass arguments to java atm)

```batch
java -Dos.name=windows -jar selenium-server-standalone.jar
```

This is actually using `java -D` which is used to set a system property. The system property we set is `os.name`. This is nothing to do with selenium-server and the appearance of `Dos` is a coincidence not related to DOS. 

Now we get a new issue. 

```R
require(RSelenium)
# RSelenium::startServer() # if needed
remDr <- remoteDriver(browserName = "opera")
remDr$open()

remDr$value$message
```

```
## [1] "Unable to find executable for product Opera Desktop"
```

So in this case we refer to the operadriver [wiki](http://code.google.com/p/selenium/wiki/OperaDriver). It states that the OperaDriver server expects you to have Opera or Opera Next installed in the default location for each system which for windows is %PROGRAMFILES%\Opera\opera.exe or more precisely C:\Program Files\opera\opera.exe. As I have installed the 64bit version I need to tell opera driver where to look

```R
require(RSelenium)
# RSelenium::startServer() # if needed
remDr <- remoteDriver(
  browserName = "opera",
  extraCapabilities = list("opera.binary" = "C:\\Program Files\\Opera x64\\opera.exe")
)
remDr$open()

head(remDr$sessionInfo)
```

```
## $platform
## [1] "ANY"
## 
## $opera.binary
## [1] "C:\\Program Files\\Opera x64\\opera.exe"
## 
## $javascriptEnabled
## [1] TRUE
## 
## $opera.host
## [1] "127.0.0.1"
## 
## $browserName
## [1] "opera"
## 
## $opera.idle
## [1] FALSE
```

A few small issues. I suspect if you were running win 7 or lower and the 32 bit version of opera 12.16 it would probably run out of the box.


## RSelenium with Local Headless Browsers

Next we shall look at running what is known as headless browsers. Usually a browser can do three things 

1. For given url, download the html page (or any other content apart from html)
1. Render the content into dom, eg executing javascript inside the script tag. and the executed result will be reflected on the browsers dom.
1. Render the dom into visualised content.

A headless browser handles items 1 and 2 but doesn't carryout 3. This means it doesn't display anything. All pages etc. are in memory rather then displayed to the user. The result of this is that headless browsers should perform faster then their full-fledged competitors which could be welcome news to speed up testing.

### phantomjs

The first headless browser we shall look at is `phantomjs`. Firstly download the relevant zip file for your OS from [here](http://phantomjs.org/download.html). We are using windows so we downloaded [ phantomjs-2.1.1-windows.zip](https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-windows.zip). It is sufficient to place the location of the directory containing `phantomjs.exe` in your path. In this case we probably could have just extracted `phantomjs.exe` to the Documents folder where chromedriver etc current reside.

However I extracted it to the desktop keeping the contents of the zip. The reasoning behind this was that phantomjs is driven by selenium using [ghostdriver](https://github.com/detro/ghostdriver). At some point the version of ghostdriver phantomjs uses will be upgraded and will accept calls from an unexposed method `phantomExecute` of the RSelenium `remoteDriver` class. There are interesting scripts contained in the phantomjs /example directory like netsniff.js which captures network traffic in HAR format. When the `phantomExecute` method is exposed these scripts will be useful. So I added the location of the .exe to my path namely the directory C:\Users\john\Desktop\phantomjs-1.9.7-windows. Once your operating system can find `phantomjs.exe` or the equivalent driving a phantomjs browser is as easy as:

```R
require(RSelenium)
# RSelenium::startServer() # if needed
remDr <- remoteDriver(browserName = "phantomjs")
remDr$open()

head(remDr$sessionInfo)
```

```
## $platform
## [1] "XP"
## 
## $acceptSslCerts
## [1] FALSE
## 
## $javascriptEnabled
## [1] TRUE
## 
## $browserName
## [1] "phantomjs"
## 
## $rotatable
## [1] FALSE
## 
## $driverVersion
## [1] "1.1.0"
```

We can take a screenshot even thou the browser is headless:

```R
remDr$navigate("http://www.google.com/ncr")
remDr$screenshot(display = TRUE)
```

<h6 align = center>RSelenium on win 8.1 driving phantomjs</h6>

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/tmpScreenShot.png"  title = "RSelenium on win 8.1 driving phantomjs" style = "margin-left: auto;margin-right: auto; display: block;"  width = '100%'/>

PhantomJS is excellent. It has only recently as of version 1.8 had Ghost Driver integration and hopefully its importance will increase further.

### HtmlUnit

The original headless browser for selenium was `htmlunit`. 


## RSelenium with Local Mobile Browsers

Configuring your local machine to use mobile browsers can be slightly tricky. If you are having difficulty setting up on your particular OS you may want to skip this section.

### Android

The first mobile browser we will look at driving is Android. The selenium project had android drivers in the selenium project. The current state of these drivers is listed [here](http://code.google.com/p/selenium/wiki/AndroidDriver#REMOVED_FROM_THE_PROJECT).

As can be noted driving android using the selenium server has been deprecated in favour of the [selendroid project](http://selendroid.io/mobileWeb.html). Once selendroid has been setup this means that rather than running the selenium standalone jar as a server we will be running an equivalent selendroid jar to drive our browser on our real or emulated android phone. More on this later for now we will look at setting selendroid up.

#### Setting up Selendroid

The selendroid project has a page on [getting started](http://appium.io/docs/en/about-appium/getting-started/).

##### JDK

There are a couple of things to note **Java SDK (minimum 1.6) must be installed**. Most likely JRE is installed on your system. To check look for the directory C:\Program Files\Java\jdk1.7.0_51 or something similar. You can also check the version of java and it should indicate `Java(TM) SE Runtime Environment`. 

```batch
java -version
```

```
java version "1.7.0_51"
Java(TM) SE Runtime Environment (build 1.7.0_51-b13)
Java HotSpot(TM) 64-Bit Server VM (build 24.51-b03, mixed mode)
```

If you need the JDK you can install from [here](http://www.oracle.com/technetwork/java/javase/downloads/java-se-jre-7-download-432155.html). Once the JDK is installed the environment variable JAVA_HOME should be set.

```batch
echo %JAVA_HOME%
```

```
C:\Program Files\Java\jdk1.7.0_51
```

##### SDK

Another Development kit needs to be installed [SDK](http://developer.android.com/sdk/index.html) in this case. From the link provided you can download the ADT Bundle for Windows. I downloaded the bundle and extracted the zip to the Desktop this resulted in a directory `C:\Users\john\Desktop\adt-bundle-windows-x86_64-20131030`. There is a guide to setup the SDK [here](http://projects.spring.io/spring-android/#quick-start). The environment variable `ANDROID_HOME` needs to be set. This should be set to the /sdk/ directory in the extracted bundle.

```batch
echo %ANDROID_HOME%
```

```
C:\Users\john\Desktop\adt-bundle-windows-x86_64-20131030\sdk
```

After setting `ANDROID_HOME` `%ANDROID_HOME%\tools;%ANDROID_HOME%\platform-tools` needs to be added to the system path also. Typing android in a command prompt should bring up the following:

<h6 align = center>Android SDK on win 8.1</h6>

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/androidsdk.png"  title = "Android SDK on win 8.1" style = "margin-left: auto;margin-right: auto; display: block;"  width = '100%'/>

From the tools make sure all are installed. From the latest android release Android 4.4.2 (API 19) in this case make sure all are installed. From the Extras folder, select the checkbox for the Android Support Library and make sure it is installed. You will also want to install the Intel Hardware Accelerated Execution Manager. Instructions on how to do so are [here](http://software.intel.com/en-us/android/articles/installation-instructions-for-intel-hardware-accelerated-execution-manager-windows). Basically checking the box Intel x86 Emulator Accelerator (HAXM) and "installing" will download it to `%ANDROID_HOME%/extras/intel`. In this folder is an exe `IntelHaxm.exe` which should be ran to finish the install.

##### Creating an Android Virtual Device

Next we need to emulate a phone. The alternative is to use a hardware phone running the Android OS. Refer to [here](http://selendroid.io/setup.html#androidDevices) if you are running a hardware phone but it should be as simple as connecting it to the machine running selendroid via usb. We shall instead create an Android Virtual Device (Avd). The easiest way to do this is by typing `android avd` into a command console.

<h6 align = center>Android AVD on win 8.1</h6>

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/androidavd.png"  title = "Android AVD on win 8.1" style = "margin-left: auto;margin-right: auto; display: block;"  width = '100%'/>

You can see I have created an avd already named "my_avd". You will need to create one by clicking new. The details for the "my_avd" are shown here

<h6 align = center>Android AVD my_avd on win 8.1</h6>

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/androidmy_avd.png"  title = "Android AVD my_avd on win 8.1" style = "margin-left: auto;margin-right: auto; display: block;"  width = '100%'/>

Click ok and an avd should have been created. You can start it using the panel. Click the newly created avd then click start. It will take a few moments but a panel containing a virtual phone will hopefully boot up and eventually get to the phone screen.

<h6 align = center>Android AVD my_avd on win 8.1</h6>

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/androidrunavd.png"  title = "Android AVD my_avd on win 8.1" style = "margin-left: auto;margin-right: auto; display: block;"  width = '100%'/>

If you have got to this point it is most likely you will now be able to drive this phone with selendroid. Take a few moments to play with your virtual android phone ;).

##### Selendroid

Finally, we need to download the selendroid driver. The selendroid [homepage](http://selendroid.io/) has a link at the bottom to the most recent [jar](https://github.com/selendroid/selendroid/releases/download/0.8.0/selendroid-standalone-0.8.0-with-dependencies.jar).
This jar should be saved on the local computer. I saved the jar in the Documents folder namely `C:\Users\john\Documents`. To begin running tests you need to run the jar with java
in a command prompt.

```batch
java -jar selendroid-standalone-0.8.0-with-dependencies.jar
```

##### RSelenium and Selendroid

Finally we are ready to use `RSelenium` to control an android browser.

```R
require(RSelenium)
remDr <- remoteDriver(browserName = "android")
remDr$open()
remDr$navigate("http://www.google.com/ncr")
```

<h6 align = center>Android with RSelenium on win 8.1</h6>

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/androidRSelenium.png"  title = "Android with RSelenium on win 8.1" style = "margin-left: auto;margin-right: auto; display: block;"  width = '100%'/>

### iOS

Testing iOS requires a Mac. The setup is similar to selendroid. [Appium](http://appium.io/docs/en/about-appium/getting-started/?lang=en) and [ios-driver](http://ios-driver.github.io/ios-driver/?page=setup) can be used. An SDK is required namely iPhone simulator SDK and a virtual phone is ran in a similar fashion to Selendroid.


## RSelenium with Remote Browsers and External Sites

Setting up multiple OS/browsers combinations locally is not always the best use of ones time. It is an interesting exercise to implement for example an Android platform locally but the overhead of having multiple systems and browsers quickly overcomes the utility. There are professional service providers who maintain a suite of OS/browsers that can be utilised for testing. Companies such as [https://saucelabs.com/](Sauce Labs) and [http://www.browserstack.com](Browser Stack) offer free automated testing to open source projects. In this vignette we will demonstrate remote testing using Sauce Labs. 

### Setting up Sauce Labs

We assume in this vignette that you are setting up Sauce Labs for an open source project. Firstly you should register your project [here](https://saucelabs.com/opensauce). On the [account page](https://saucelabs.com/account) you will find the access key for your account that is in a similar format to `49953c74-5c46-4ff9-b584-cf31a4c71809`. Using Sauce Labs is pretty straightforward. You need to tell it what OS/Browser combination you would like. A list of possible setups can be viewed [here](https://saucelabs.com/platforms). As an example lets suppose we wished to run google chrome version 33 on OSX version mavericks. 

```R
require(RSelenium)
user <- "rselenium0" # Your Sauce Labs username
pass <- "49953c74-5c46-4ff9-b584-cf31a4c71809" # Your Sauce Labs access key 
port <- 80
ip <- paste0(user, ':', pass, "@ondemand.saucelabs.com")
rdBrowser <- "chrome"
version <- "33"
platform <- "OS X 10.9"
extraCapabilities <- list(
  name = "RSelenium OS/Browsers vignette first example",
  username = user,
  accessKey = pass,
  tags = list("RSelenium-vignette", "OS/Browsers-vignette")
)
remDr <- remoteDriver$new(
  remoteServerAddr = ip,
  port = port,
  browserName = rdBrowser,
  version = version,
  platform = platform,
  extraCapabilities = extraCapabilities
)
```

We state the browser and OS we require (chrome 33/ Mac OSX 10.9). The user and password are used to form an appropriate ip address for our remote server ([http://rselenium0:49953c74-5c46-4ff9-b584-cf31a4c71809@ondemand.saucelabs.com](http://rselenium0:49953c74-5c46-4ff9-b584-cf31a4c71809@ondemand.saucelabs.com) in this case). They are also passed as `extraCapabilities` to the remote Selenium server.

We give our test a `name` and any additional `tags` we wish that are passed to the remote Selenium server. Details of the name and tags are given [here](https://wiki.saucelabs.com/display/DOCS/Test+Configuration+Options). They are used to annotate our tests.

#### Basic Example

As a basic first example we will run a script using the mavericks/ chrome 33 combination. We run the following commands:

```R
testScript <- function(remDr){
  remDr$open()
  remDr$navigate("http://www.google.com/ncr")
  Sys.sleep(2)
  # highlight the query box
  remDr$findElement("name", "q")$highlightElement()
  Sys.sleep(2)
  # goto rproject
  remDr$navigate("http://www.r-project.org")
  # go Back
  remDr$goBack()
  # go Forward
  remDr$goForward()
  Sys.sleep(2)
  webElems <- remDr$findElements("css selector", "frame")
  # highlight the frames
  lapply(webElems, function(x){x$highlightElement()})
  
  remDr$close()
}

testScript(remDr)
```

And that's it. We have ran our first remote test using Sauce Labs. The results of the test can be viewed [https://saucelabs.com/tests/ae22f859de8746f9bfedad2f49c1c329](https://saucelabs.com/tests/ae22f859de8746f9bfedad2f49c1c329). I think you will agree its a nice setup. We have access to screenshots of all the commands we issued and a video (screencast) of the test run. We can view the selenium server logs and the medadata associated with our test.

### Testing Multiple OS/Browsers

We can easily extend the simple test we ran for multiple OS/Browser combinations. The browser and platform variables need to be assigned the combinations we require.

```R
osBrowser <- list(
  "OS X 10.9" = list(
    browser = list("safari", "firefox"),
    version = list('7', '28')
  ),
  "Windows 8" = list(
    browser = list("chrome", "firefox", "internet explorer"),
    version = list('33', '28', '10')
  ),
  "Linux" = list(
    browser = list("chrome", "firefox", "opera"),
    version = list('33', '28', '12')
  )
)
lapply(seq_along(osBrowser), function(x) {
  platform <- names(osBrowser)[x]
  lapply(seq_along(osBrowser[[x]]$browser), function(y){
    rdBrowser <- osBrowser[[x]]$browser[[y]]
    version <- osBrowser[[x]]$version[[y]]
    remDr <- remoteDriver$new(
      remoteServerAddr = ip,
      port = port,
      browserName = rdBrowser,
      version = version,
      platform = platform,
      extraCapabilities = extraCapabilities
    )
    testScript(remDr)
  })
})
```

To view the results you can goto the `RSelenium` project page on Sauce Labs [https://saucelabs.com/u/rselenium](https://saucelabs.com/u/rselenium). Listed here are all the tests ran on the `RSelenium` package. A partial search by name `Browsers vignette first example` will give the results of this test. There are a few repeats of the first Mavericks/ chrome 33 test where I tuned the script. 

<h6 align = center>Sauce Labs results for simple test script</h6>

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/saucelabsOSBrowsers.png"  title = "Sauce Labs results for simple test script" style = "margin-left: auto;margin-right: auto; display: block;"  width = '100%'/>

## RSelenium with Remote Browsers and Local Sites

Testing external webpages and websites across a range of operating systems and browsers can be achieved using Sauce Labs as we observed in the last section. Often however especially in a development phase of a project we either do not have or do not want an external version of our webpage/website/webapp. A good example would be our `shinytestapp`. Lets open a new R session and run our testapp. 

```R
require(shiny)
runApp(file.path(find.package("RSelenium"), "apps", "shinytestapp"), port = 3000)
```

To access our app we would require the ip address `http://localhost:3000/`. How do we access this from a remote webdriver?

### Sauce Connect

Thankfully Sauce Labs have a solution for this known as [Sauce Connect](https://saucelabs.com/docs/connect). Sauce Connect is a secure tunneling app which allows you to execute tests securely when testing behind firewalls via a secure connection between Sauce Labs client cloud and your environment. This allows you to drive an external Browser and have it interact with a local webpage/website/webapp. 

#### Setting up Sauce Connect

Firstly you need to download the Sauce Connect zip for the operating system you are using to run your tests on. This machine I will be testing from is running windows 8.1 so I download the [windows](https://saucelabs.com/downloads/sc-4.0-latest-win32.zip) zip. I unzipped Sauce Connect to the Documents folder so it is now accessible at `C:\Users\john\Documents\sc-4.1-win32`. From a windows command prompt we navigate to the Sauce Connect bin directory and run the .exe file supplying our Sauce Labs user name and access key. 

```batch
sc.exe -u rselenium0 -k 49953c74-5c46-4
```

```
ff9-b584-cf31a4c71809
```

### Basic Example

We opened our shinytestapp on port 3000 because Sauce Connect only supports a set number of ports. All ports can be used but for this you need a locally-defined domain name (which can be set in your hosts file) rather than localhost. This is simple to do but for the purposes of this vignette we shall connect to `http://localhost:3000`. Again to start with we shall use Mavericks with Chrome 33.

```R
require(RSelenium)
user <- "rselenium0" # Your Sauce Labs username
pass <- "49953c74-5c46-4ff9-b584-cf31a4c71809" # Your Sauce Labs access key 
```
~~port <- 80~~<br>
~~ip <- paste0(user, ':', pass, "@ondemand.saucelabs.com")~~

```
port <- 4445
ip <- paste0(user, ':', pass, "@localhost")
rdBrowser <- "firefox"
version <- "26"
platform <- "Linux"
extraCapabilities <- list(
  name = "RSelenium OS/Browsers vignette second example",
  username = user,
  accessKey = pass,
  tags = list(
    "RSelenium-vignette",
    "OS/Browsers-vignette",
    "Example 2"
  )
)
remDr <- remoteDriver$new(
  remoteServerAddr = ip,
  port = port,
  browserName = rdBrowser,
  version = version,
  platform = platform,
  extraCapabilities = extraCapabilities
)
```

Everything is as before the exception is that when we ask to browse to a localhost address Sauce Connect will intervene. Also we connect to Sauce Labs through Sauce Connect on `localhost:4445` by default instead of `ondemand.saucelabs.com:80`.

```R
localScript <- function(remDr){
  remDr$open()
  remDr$setImplicitWaitTimeout(2000) # wait for elements for 2 seconds
  remDr$navigate("http://localhost:3000")
  Sys.sleep(2)
  # highlight the labels
  webElems <- remDr$findElements("css selector", "#ctrlSelect span")
  lapply(webElems, function(x) {x$highlightElement()})
  Sys.sleep(2)
  appIds <- c("summary", "distPlot", "ggPlot", "dttable")
  # Click each checkbox and check for its output
  lapply(seq_along(webElems), function(x) {
    if(!webElems[[x]]$isElementSelected()[[1]]) {
      webElems[[x]]$clickElement()
      # test for its output
      out <- remDr$findElement("id", appIds[x])
      out$highlightElement()
    }
  })
  
  remDr$close()
}

localScript(remDr)
```
---
title: "Docker"
output:
  html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Using RSelenium with Docker}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


## Introduction

The goal of this document to illustrate the use of docker and the Selenium project. The selenium project has released a number of containers to run various tools in the suite. 

### Why Docker?

Running a docker container standardises the build across OS's and removes many of the issues
user may have relating to JAVA/browser version/selenium version etc.


## Preparing Docker on Windows

The following outlines preparing a windows 10 home installation to run docker containers.

### Installing Docker

Depending on your version of windows there are currently two ways of running docker.
Details are given at [Docker](https://docs.docker.com/engine/installation/windows/). For the version of windows Being used in this vignette (win 10 home) the older [Docker toolbox](https://www.docker.com/products/docker-toolbox) is what we are going to use to run containers.

Download the `.exe` for Docker toolbox and begin the install. Follow the instructions on the [toolbox install windows](https://docs.docker.com/toolbox/toolbox_install_windows/) page.

Clicking on the Docker Quickstart terminal link should eventually give you a terminal that resembles the following:

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203296/windowsToolbox.png" title = "Docker toolbox on Win 10 Home" width = '100%'/>


## Preparing Docker on Linux

The following outlines preparing a ubuntu 16.04 installation to run docker containers.

### Installing Docker

For a general guide for various Linux distros see [Docker](https://docs.docker.com/engine/installation/linux/)

#### Ubuntu 16.04

Run system update firstly

```sh
sudo apt-get update
```

Add the official Docker GPG key  to your system

```sh
sudo apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys 58118E89F3A912897C070ADBF76221572C52609D
```


Add the Docker repository to APT sources:

```sh
echo "deb https://apt.dockerproject.org/repo ubuntu-xenial main" | sudo tee /etc/apt/sources.list.d/docker.list

```

Update the package database:

```sh
sudo apt-get update
```

Now you should be able to install Docker

```sh
sudo apt-get install -y docker-engine
```

Start the docker daemon.

```sh
sudo service docker start
```

You can check that the Docker is working by running a test container:

For now run as `sudo`

```sh
sudo docker run hello-world
```

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/ubuntuHelloWorld.png" title = "Docker hello-world Ubuntu 16.04" width = '100%'/>


#### Running without sudo

If you want to run `docker` without `sudo` you can create a docker group and add the appropriate 
user to that group see [Docker](https://docs.docker.com/engine/installation/linux/ubuntulinux/#create-a-docker-group).


## Using Selenium  Docker Images

The official docker repository can be viewed [here](https://hub.docker.com/u/selenium/) at Docker Hub. 

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/SeleniumDockerHub.png" title = "Official Selenium repository at Docker Hub" width = '100%'/>

### What's in an Image?

Lets examine the [`selenium/node-firefox`](https://hub.docker.com/r/selenium/node-firefox/) image. The first thing to note is the `Tags`. 

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/node-firefox.png" title = "Selenium node-firefox image tags" width = '100%'/>

Clicking this gives us the various tagged versions of this image. Clicking Dockerfile allows us
to view the source behind the image:


```
FROM selenium/node-base:2.53.0 
MAINTAINER Selenium <selenium-developers@googlegroups.com> 

USER root #========= 
# Firefox 
#========= 
ENV FIREFOX_VERSION 45.0.2 
```

We can see that this image will run selenium server and contain firefox version 45.0.2. 

We will be interested not in the node-firefox image but the [selenium/standalone-firefox](https://hub.docker.com/r/selenium/standalone-firefox/) image. Looking at its Dockerfile we see however that it is built on the first image we looked at.

### Pulling an Image

#### Windows 10 - Docker Toolbox

To get the image we are interested in we first click our quickstart link to open a terminal. 

We then issue the following command:

```sh
docker pull selenium/standalone-firefox:2.53.0
```

Notice we appended a tag version so we are asking for the `2.53.0` version of this image.

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203296/Win10PullFirefox.png" title = "Official Selenium repository at Docker Hub" width = '100%'/>

#### Ubuntu 16.04 

Similar to windows we use the `docker pull` command. If you did not create and add your user to a `Docker` group you will need to `sudo`

```sh
sudo docker pull selenium/standalone-firefox:2.53.0
```

### Start Your Servers!

When we run an image we refer to that instance of the running image as a container. Importantly for Selenium we are dealing with a server so we require a port. We will have a port on the container and a port on the host machine (the machine that is running the container). The only complexity will involve mapping the host port and container port. 

By default the Selenium images expose port 4444 and use this port for serving. This is not a problem for us and we can map as we see fit on the host machine.

#### Windows 10 Running an Image

Again we run docker and this time use the [`run`](https://docs.docker.com/engine/reference/commandline/run/) command. 

```sh
docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.0
```

For illustration here we are mapping the container port 4444 to the host port 4445 so when we want to refer to our server on the host machine we will use 4445 as the port. 

We also run the [`ps`](https://docs.docker.com/engine/reference/commandline/ps/) command after running to show our current running containers.

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/SeleniumRunWin10.png" title = "Running standalone-firefox:2.53.0 Win 10" width = '100%'/>

So we can see an instance of standalone-firefox:2.53.0 is running. It has `id=ac61567a8f06` and is named `peaceful_aryabhata`.


#### Ubuntu 16.04 Running an Image

For ubuntu we also issue a `run` command followed by a listing of running containers with `ps`.

Again run as `sudo` if necessary:

```sh
sudo docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.0
sudo docker ps
```

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203293/SeleniumRunUbuntu.png" title = "Running standalone-firefox:2.53.0 Win 10" width = '100%'/>

The `Ubuntu` container has `id=5475923fc4dc` and is named `peaceful_payne`.

A tip if you don't want to strain your eyes is to use the `--format` argument:

```sh
sudo docker ps --format 'table {{.Names}}\t{{.Image}}\t{{.ID}}'
```

```
NAMES               IMAGE                                CONTAINER ID
peaceful_payne      selenium/standalone-firefox:2.53.0   5475923fc4dc
```

### RSelenium and Running Containers

Now that we have a container running a selenium server and its own firefox we can call it using `RSelenium`

#### Utilising the Firefox Standalone Container in Windows 10

With Windows 10 we have Docker toolbox which creates a virtual machine on windows to run our containers. So the virtual machine is the host rather than Win 10. We need the ip address of this virtual machine. This is the only added complication. We can get the ip address of the virtual machine using `docker-machine.exe`:

```sh
docker-machine ip
```

```
192.168.99.100
```

So as we can see the virtual machine has `ip = 192.168.99.100` and as we noted previously we mapped port 4444 to port 4445 on the host (the virtual machine in this case). So we should have a Selenium server accessible on port 4445 of the virtual machine. Lets try it:

```R
library(RSelenium)
remDr <- remoteDriver(
  remoteServerAddr = "192.168.99.100",
  port = 4445L
)
remDr$open()
```

```
## [1] "Connecting to remote server"
## $applicationCacheEnabled
## [1] TRUE
## 
## $rotatable
## [1] FALSE
## 
## $handlesAlerts
## [1] TRUE
## 
## $databaseEnabled
## [1] TRUE
## 
## $version
## [1] "45.0.2"
## 
## $platform
## [1] "LINUX"
## 
## $nativeEvents
## [1] FALSE
## 
## $acceptSslCerts
## [1] TRUE
## 
## $webdriver.remote.sessionid
## [1] "bfe1c24d-8b17-4e64-a062-cba14e99472b"
## 
## $webStorageEnabled
## [1] TRUE
## 
## $locationContextEnabled
## [1] TRUE
## 
## $browserName
## [1] "firefox"
## 
## $takesScreenshot
## [1] TRUE
## 
## $javascriptEnabled
## [1] TRUE
## 
## $cssSelectorsEnabled
## [1] TRUE
## 
## $id
## [1] "bfe1c24d-8b17-4e64-a062-cba14e99472b"
```

```R
remDr$navigate("http://www.google.com/ncr")
remDr$getTitle()
```

```
## [[1]]
## [1] "Google"
```

So as expected we have version `45.0.2` of firefox running in `Linux`. We note that no browser popped up and we don't see any signs in the `MINGW64` console. This is expected the browser is headless. If you are using `RStudio` you can use the `Viewer` pane to get a screenshot `remDr$screenshot(display = TRUE)`. Later we will see how we can run a debug version that uses `VNC` to allow us to interact with the containers browser.

#### Utilising the Firefox Standalone Container in Ubuntu 16.04

On our ubuntu 16.04 box life is simpler. The host in this case is the machine itself. The ip address is therefore `localhost` the default for `RSelenium`:

```R
library(RSelenium)
remDr <- remoteDriver(port = 4445L)
remDr$open()
```

```
## [1] "Connecting to remote server"
## $applicationCacheEnabled
## [1] TRUE
## 
## $rotatable
## [1] FALSE
## 
## $handlesAlerts
## [1] TRUE
## 
## $databaseEnabled
## [1] TRUE
## 
## $version
## [1] "45.0.2"
## 
## $platform
## [1] "LINUX"
## 
## $nativeEvents
## [1] FALSE
## 
## $acceptSslCerts
## [1] TRUE
## 
## $webdriver.remote.sessionid
## [1] "644c353a-34b8-4bb4-bcff-746df5a30af8"
## 
## $webStorageEnabled
## [1] TRUE
## 
## $locationContextEnabled
## [1] TRUE
## 
## $browserName
## [1] "firefox"
## 
## $takesScreenshot
## [1] TRUE
## 
## $javascriptEnabled
## [1] TRUE
## 
## $cssSelectorsEnabled
## [1] TRUE
## 
## $id
## [1] "644c353a-34b8-4bb4-bcff-746df5a30af8"
```

```R
remDr$navigate("http://www.google.com/ncr")
remDr$getTitle()
```

```
## [[1]]
## [1] "Google"

```

Again we see that we have `Firefox` 45.0.2 running on `Linux`. So running docker containers allows us to run a standard setup OS independent and we don't need to worry about machine specifics.

## Debugging Using VNC

Often it is useful to be able to interact with the browser spawned by Selenium. In this regard the `Selenium` project provides debug versions of their images which include `VNC`. This means with a VNC viewer we can remote logon to the container and interact with the browser.

### Remote Logging/Debugging with Windows

First lets stop the running container. To stop the running container we need its ID. We then pass this id to the [`stop`](https://docs.docker.com/engine/reference/commandline/stop/) command:

```sh
docker ps -q
```

```
ac61567a8f06
```

```sh
docker stop $(docker ps -q)
```

```
ac61567a8f06
```

```sh
docker ps
```

```
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS              PORTS               NAMES
```

So we have stopped our container and will have no issues with blocked ports. Next we need to [`pull`](https://docs.docker.com/engine/reference/commandline/pull/) the [`Firefox standalone debug`](https://hub.docker.com/r/selenium/standalone-firefox-debug/) image. Again this has various tags and we will choose the `2.53.0` one:

```sh
docker pull selenium/standalone-firefox-debug:2.53.0
```

```
2.53.0: Pulling from selenium/standalone-firefox-debug

0be59000882d: Already exists
f20b6f990572: Already exists
53662c966c9f: Already exists
a3ed95caeb02: Already exists
0e449738cbb6: Already exists
63921592acdf: Already exists
8553f3252cbc: Already exists
dde80f7d7068: Already exists
232a7d285d74: Already exists
284e6949cb6f: Already exists
85909c621e51: Already exists
db03d98be095: Already exists
67898bba08c8: Already exists
23edc5cca893: Already exists
5ceef68785db: Already exists
8cbb6b75dad5: Already exists
d18b13cc4e51: Pull complete
4a3e92e6c855: Pull complete
1f5786a4b2ff: Download complete
3a29027ac709: Download complete
a70ecef6edbc: Download complete
10672504dfa1: Pull complete
1f5786a4b2ff: Pull complete
3a29027ac709: Pull complete
a70ecef6edbc: Pull complete
Digest: sha256:9b03c047b68b4e1c43f6d09edd952c1200766408104f76fa055e04028390491b
Status: Downloaded newer image for selenium/standalone-firefox-debug:2.53.0
```

Next we start the container as before. This time we have two ports to map. The Selenium port as before and also now the VNC port. The Selenium image expose the VNC port on `5900` by default and we map this port to a port on the host. We will choose `5901` just to illustrate host v container:

```sh
docker run -d -p 4445:4444 -p 5901:5900 selenium/standalone-firefox-debug:2.53.0
```

```
ee0c0bb8b711723e653ccc26219e314a37c28d2027f939046adcfc90a4beb645
```

```sh
docker-machine ip
```

```
192.168.99.100
```

```sh
docker ps --format 'table {{.Names}}\t{{.Image}}\t{{.ID}}'
```

```
NAMES               IMAGE                                      CONTAINER ID
focused_tesla       selenium/standalone-firefox-debug:2.53.0   ee0c0bb8b711
```

Again we note that the container is ran on a virtual machine not on our Windows 10 box so we expose the ip of said virtual machine. Next we download a VNC viewer if we don't have one. 

[TightVNC](http://www.tightvnc.com/download.php) has worked for me so I use it here.

Our Selenium server is running so lets use RSelenium again to browse to `Google`:

```R
library(RSelenium)
remDr <- remoteDriver(
  remoteServerAddr = "192.168.99.100",
  port = 4445L
)
remDr$open(silent = TRUE)
remDr$navigate("http://www.google.com/ncr")
remDr$getTitle()
remDr$screenshot(display = TRUE)
```

Now we open `TightVNC` and use the VM ip with appropriate port for VNC:

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/TightVNC.png" title = "Accessing Firefox debug with Tight VNC" width = '100%'/>

You will be asked for a password which is `secret`. This can be seen by reading the images dockerfile:

```
RUN apt-get update -qqy \
  && apt-get -qqy install \
    x11vnc \
  && rm -rf /var/lib/apt/lists/* \
  && mkdir -p ~/.vnc \
&& x11vnc -storepasswd secret ~/.vnc/passwd
```
Enter the password and you should have access to the container browser. Navigate to another website;

```R
remDr$navigate("http://www.bbc.com")
remDr$getTitle()
```

```
## [[1]]
## [1] "BBC - Homepage"
```

and you should observe the action on the containers browser. Finally stop your containers:

```sh
docker stop $(docker ps -q)
```

```
ee0c0bb8b711
```

### Remote Logging/Debugging with Linux

Firstly, we stop any running containers which may block our ports:

```sh
docker stop $(docker ps -q)
```

```
5475923fc4dc
```

```sh
docker ps
```

```
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS              PORTS               NAMES
```

Like in windows we `pull` the `2.53.0` version of the `Firefox` standalone debug image (again using `sudo` if necessary:

```sh
sudo docker pull selenium/standalone-firefox-debug:2.53.0
```

```
2.53.0: Pulling from selenium/standalone-firefox-debug

0be59000882d: Already exists 
f20b6f990572: Already exists 
53662c966c9f: Already exists 
a3ed95caeb02: Already exists 
0e449738cbb6: Already exists 
63921592acdf: Already exists 
8553f3252cbc: Already exists 
dde80f7d7068: Already exists 
232a7d285d74: Already exists 
284e6949cb6f: Already exists 
85909c621e51: Already exists 
db03d98be095: Already exists 
67898bba08c8: Already exists 
23edc5cca893: Already exists 
5ceef68785db: Already exists 
8cbb6b75dad5: Already exists 
d18b13cc4e51: Pull complete 
4a3e92e6c855: Pull complete 
10672504dfa1: Downloading 7.915 MB/16.87 MB
1f5786a4b2ff: Downloading 5.487 MB/5.528 MB
3a29027ac709: Download complete 
a70ecef6edbc: Download complete 
```

Once the image is downloaded we can start it as we did with Windows. We map Selenium to port 4445 on the host and VNC to port 5901. As noted before the host in the case of Ubuntu 16.04 is the machine itself so `localhost` is the relevant ip:

```sh
sudo docker run -d -p 4445:4444 -p 5901:5900 selenium/standalone-firefox-debug:2.53.0
```

```
dd36811562e83e66ab145efcedca2825a621741f10d8c81a5fb3fb5ba3019032
```

```sh
sudo docker ps --format 'table {{.Names}}\t{{.Image}}\t{{.ID}}'
```

```
NAMES               IMAGE                                      CONTAINER ID
thirsty_swanson     selenium/standalone-firefox-debug:2.53.0   dd36811562e8
```

So we can see that Selenium server is running. Like in Windows we need a VNC viewer. You can use `TightVNC` as with windows and a guide to running it is [here](https://www.digitalocean.com/community/tutorials/how-to-install-and-configure-vnc-on-ubuntu-16-04). We will use `Vinagre` however. A guide to install `Vinagre` is [here](https://www.howtoinstall.co/en/ubuntu/xenial/vinagre). Basically we just need to `apt-get` it:

```sh
sudo apt-get update
sudo apt-get install vinagre
```

Again, we navigate to Google:

```R
library(RSelenium)
remDr <- remoteDriver(port = 4445L)
remDr$open(silent = TRUE)
remDr$navigate("http://www.google.com/ncr")
remDr$getTitle()
```

```
## [[1]]
## [1] "Google"
```

Now open Vinagre and input the correct ip/port (127.0.0.1:5901 in this case):

<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203295/VinagreVNC.png" title = "Accessing Firefox debug with Tight VNC" width = '100%'/>

Again it will ask for a password which as before is `secret`.

You should now be able to interact with the containers browser. Again navigate to another website and observe the containers browser:

```R
remDr$navigate("http://www.bbc.com")
remDr$getTitle()
```

```
## [[1]]
## [1] "BBC - Homepage"
```

Finally close the session:

```
remDr$close()
```

And close down any running containers:

```sh
sudo docker stop $(sudo docker ps -q)
```

```
dd36811562e8
```
---
title: "Internet Explorer"
author: "Sebastian Wolf"
output:
  html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Internet Explorer}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


## Introduction

This tutorial shall show you creating a setup that allows you to test web apps using Selenium Server + a connection to Microsoft Internet Explorer. It contains the most important tricks in Microsoft Windows you'll need to perform. Additionally some extra information is given on how to change default methods like clicking to run stable in Internet Explorer.


## Windows Registry setup

### Admin rights

You will need administrator rights to perform all steps in this chapter

### Edit Registry Main

To allow the Internet Explorer Selenium connection there are certain settings in the Windows Registry that need to be changed.

Open registry by `regedit` command on Windows

Create the Key:

```
HKEY_LOCAL_MACHINE\SOFTWARE\Wow6432Node\Microsoft\Internet Explorer\Main\FeatureControl\FEATURE_BFCACHE
```

Please note that the `FEATURE_BFCACHE` subkey may or may not be present, and should be created if it is not present.

**Important**: Inside this key, create a `DWORD` value named `iexplore.exe` with the value of `0`.

### Edit Registry User

Create the Key:

```
HKEY_LOCAL_MACHINE \SOFTW ARE\Microsoft\Internet Explorer\Main\FeatureControl\FEATURE_BFCACHE 
```

Please note that the `FEATURE_BFCACHE` subkey may or may not be present, and should be created if it is not present.

**Important**: Inside this key, create a `DWORD` value named `iexplore.exe` with the value of `0`.

### Allow Window Navigation

Go to:

```
HKEY_CURRENT_USER \Software
\Microsoft \Internet Explorer \Main
```

Inside this key (Main) , create a `DWORD` value named `TabProcGrowth` with the value of `0`.


## Selenium Server

To use Internet Explorer there is sadly just one way to use a Selenium Server which is running it via the Java Binary as explained in the **Basics** vignette of this package.


## Selenium Driver

For Internet Explorer please download the 32-bit version of the SeleniumDriver. The 64 bit version still has trouble inserting text and can make your user interface testing really slow.

Please have the `IEDriverServer.exe` in your PATH variable. You can simply do this by using 

```R
ie_driver_folder <- "C:\\Selenium"
Sys.setenv(PATH = paste(ie_driver_folder, Sys.getenv("PATH"), sep = ";"))
```

if you copied the "IEDriverServer.exe" to `C:\Selenium`


## Initialization for a Selenium Driver object

### Set extra capabilities

For initialization of a Selenium Driver object in with Internet Explorer there are certain extra settings to be made. The first one are the extraCapabilities. Not all of these are needed, but those 

- `ie.forceCreateProcessApi=FALSE` and `InternetExplorerDriver.INTRODUCE_FLAKINESS_BY_IGNORING_SECURITY_DOMAIN=TRUE` are needed to basically access the Internet Explorer from Selenium.
- `InternetExplorerDriver.IGNORE_ZOOM_SETTING=TRUE` will allow you to start Internet Explorer Selenium Driver Sessions even if you did not set the zoom to "100%" before starting your session.
- `requireWindowFocus=TRUE` allows more native browser interactions.
- `enablePersistentHover=FALSE` allows you to hover and focus elements.

So please define a list like that:

```R
extraCapabilities <- list(
  ie.forceCreateProcessApi = FALSE,
  InternetExplorerDriver.INTRODUCE_FLAKINESS_BY_IGNORING_SECURITY_DOMAIN = TRUE,
  InternetExplorerDriver.IGNORE_ZOOM_SETTING = TRUE,
  requireWindowFocus = TRUE,
  enablePersistentHover = FALSE
)
```

### Start the driver

To navigate to your first page you can now start the remoteDriver. Please note that Internet Explorer will now open and connect to the local Selenium Server. You need to have the following:

```R
remDr <- remoteDriver(
  browserName = "internet explorer",
  extraCapabilities = extraCapabilities
)
remDr$open()
remDr$setImplicitWaitTimeout(as.numeric(10000))

url <- "https://docs.ropensci.org/RSelenium"
remDr$navigate(url)
```

We use a global definition of the `remDr` element as there exists no possibility to create two parallel sessions using Internet Explorer. Additionally it is necessary to set a really long implicit wait timeout due to basic troubles Internet Explorer might have running multiple tests in a row.

### Initialization for more reproducible tests

For reproducibility reasons we noticed that in Internet Explorer you either want to always maximize the screen or set it to a fixed size. Additionally always move the Window to the top left corner of your screen. This is mainly important for checking images that you want to compare against other images created by your web app.

```R
remDr$navigate(url)
remDr$maxWindowSize()
remDr$setWindowSize(1936, 1056)
remDr$setWindowPosition(0, 0)
```

### Additional functionalities for testing shiny

Shiny may sometimes run inside `iframes`. In Internet Explorer it might be hard to get into those. Therefore in testing shiny using Internet Explorer we recommend adding a boolean variable called `in_shiny` to your sessionInfo.

```R
remDr$sessionInfo$in_shiny <- FALSE
```

This variable can be used to check if you are running inside the shiny app already or not. You do not want do go into an iframe inside the shiny app, if you are already inside the shiny app.

So after starting a Selenium Session maybe do the following:

Navigate to the mainframe

```R
remDr$sessionInfo$in_shiny <- FALSE
object$switchToFrame(NULL)
object$setImplicitWaitTimeout(1000)
```

Navigate into the first iframe if an iframe is there.

```R
iframe_found <- TRUE

if (length(remDr$findElements("tag", "iframe")) == 0 || remDr$sessionInfo$in_shiny) {
  iframe_found <- FALSE
  remDr$sessionInfo$in_shiny <- TRUE
} else {
  remDr$sessionInfo$in_shiny <- TRUE
  remDr$switchToFrame(remDr$findElements("tag", "iframe")[[1]])
}
```


## Interacting with the page

### Clicking

As simple as it might seem, during a lot of test runs using Internet Explorer for Web testing with Selenium, we found that clicking might have some hurdles. Instead of using the basic `click` functionality of Selenium we recommend either

1. Move the mouse to the element and click

```R
web_element <- remDr$findElements("tag", "a")[[1]]
remDr$mouseMoveToLocation(
 x = round(web_element$getElementSize()$width / 3),
 y = round(web_element_selector$getElementSize()$height / 3),
 webElement = web_element
)
web_element$clickElement()
```

2. Click by using javascript

```R
remDr$executeScript("arguments[0].click();", list(web_element))
```

### Entering Text in a input text field

For entering a text into a Text box in Internet Explorer we highly recommend to first set the value of the text box. Afterwards clean it and then send the character string to the textbox to type it in.

```R
web_element <- remDr$findElements("css selector", "input[type='text']")[[1]]
text_to_type = "My input text"
remDr$executeScript(
  paste0("arguments[0].setAttribute('value','", text_to_type, "');"),
  list(web_element)
)

web_element$clearElement()
web_element$sendKeysToElement(list(text_to_type))
```

### Checking a checkbox

It may seem simple, but it is one of the hardest parts using Selenium to check a checkbox. In Internet Explorer there is just one way to make it save and always happen.

**Important** You are never allowed to not have the cursor on the screen where Internet Explorer is running. You need to have the Internet Explorer Window focused.

Please always get the checkboxed focus by executing Javascript code using Selenium and afterwards click just the `input` element of this checkbox.

```R
checkboxes <- remDr$findElements("class name", "checkbox")
remDr$executeScript(
  "arguments[0].focus();",
  list(checkboxes[[1]]$findChildElements("tag", "input")[[1]])
)
checkboxes[[1]]$findChildElements("tag", "input")[[1]]$clickElement()
```
---
title: "Basics"
output:
  html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{RSelenium Basics}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


## Introduction

The goal of RSelenium is to make it easy to connect to a Selenium Server/Remote Selenium Server from within R. RSelenium provides R bindings for the Selenium Webdriver API. [Selenium](http://docs.seleniumhq.org/) is a project focused on automating web browsers. RSelenium allows you to carry out unit testing and regression testing on your webapps and webpages across a range of browser/OS combinations. This allows us to integrate from within R testing and manipulation of popular projects such as [Shiny Apps](http://www.rstudio.com/shiny/).


## Connecting to a Selenium Server

### What is a Selenium Server?

Selenium Server is a standalone java program which allows you to run HTML test suites in a range of different browsers, plus extra options like reporting.

You may, or may not, need to run a Selenium Server, depending on how you intend to use Selenium-WebDriver (`RSelenium`).

### Do I Need to Run a Selenium Server?

If you intend to drive a browser on the same machine that RSelenium is running on, you will need to have a Selenium Server running on that machine.

### How Do I Get the Selenium Server Standalone Binary?

You can download the latest Selenium Server binary manually [here](http://selenium-release.storage.googleapis.com/index.html). Look for `selenium-server-standalone-x.xx.x.jar`.

### How Do I Run the Selenium Server?

There are three ways to run a Selenium Server:

#### Docker

**The recommended way** to run a Selenium Server is by running a [Docker](https://www.docker.com/) container.

Run a server for example using Docker:

```sh
docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.1
```

Use a debug image with a VNC viewer if you wish to view the browser:

```sh
docker run -d -p 5901:5900 -p 127.0.0.1:4445:4444 --link http-server selenium/standalone-firefox-debug:2.53.1
```

There is a separate vignette which covers the using `RSelenium` with Docker see `vignette("docker", package = "RSelenium")` or [here](https://docs.ropensci.org/RSelenium/articles/docker.html).

#### `rsDriver`

For users who are not familiar with Docker, there is now a function `rsDriver` which will manage the binaries needed for running a Selenium Server. This provides a wrapper around the `wdman::selenium` function. For additional options and more control, see the [wdman](https://docs.ropensci.org/wdman/) project. Examples using the `rsDriver` function are given in the appendix. ***Please submit any issues with running binaries to the wdman project***

#### Java Binary

Alternatively, you can run the binary manually. Open a console in your OS and navigate to where the binary is located and run:

```sh
java -jar selenium-server-standalone-x.xx.x.jar
```

By default, the Selenium Server listens for connections on port 4444.

***Note for Mac OSX:*** The default port 4444 is sometimes used by other programs such as kerberos. In our examples, we use port 4445 in respect of this and for consistency with [the Docker vignette](https://docs.ropensci.org/RSelenium/articles/docker.html).

### How Do I Connect to a Running Server?

`RSelenium` has a main reference class named `remoteDriver`. To connect to a server, you need to instantiate a new `remoteDriver` with appropriate options.

```R
library(RSelenium)
remDr <- remoteDriver(
  remoteServerAddr = "localhost",
  port = 4445L,
  browserName = "firefox"
)
```

***Note for Windows:*** If you are using Docker toolbox, your remote server address will not be localhost. You need to use the ip address of the VM that is running docker.

For example:

```sh
docker-machine ip
```

```
## 192.168.99.100
```

***Note:*** See [the Docker vignette](https://docs.ropensci.org/RSelenium/articles/docker.html) for further details. The newer Docker for windows however should be accessible on the localhost.

It would have been sufficient to call `remDr <- remoteDriver(port = 4445L)`, but the options where explicitly listed to show how one may connect to an arbitrary ip/port/browser etc. More detail maybe found on [the SauceLabs vignette](https://docs.ropensci.org/RSelenium/articles/saucelabs.html).

To connect to the server, use the `open` method:

```R
remDr$open()
```

`remDr` should now have a connection to the Selenium Server. You can query the status of the remote server using the `getStatus` method:

```R
remDr$getStatus()
```

```
## $build
## $build$version
## [1] "2.53.1"
## 
## $build$revision
## [1] "a36b8b1"
## 
## $build$time
## [1] "2016-06-30 17:37:03"
## 
## 
## $os
## $os$name
## [1] "Linux"
## 
## $os$arch
## [1] "amd64"
## 
## $os$version
## [1] "4.4.0-47-generic"
## 
## 
## $java
## $java$version
## [1] "1.8.0_91"
```


## Navigating Using RSelenium

### Basic Navigation

To start with, we navigate to a url:

```R
remDr$navigate("http://www.google.com/ncr")
```

Then, we navigate to a second page:

```R
remDr$navigate("http://www.bbc.co.uk")
remDr$getCurrentUrl()
```

```
## [[1]]
## [1] "http://www.bbc.co.uk/"
```

We can go backwards and forwards using the methods `goBack` and `goForward`.

```R
remDr$goBack()
remDr$getCurrentUrl()
```

```
## [[1]]
## [1] "https://www.google.com/"
```

```R
remDr$goForward()
remDr$getCurrentUrl()
```

```
## [[1]]
## [1] "http://www.bbc.co.uk/"
```

To refresh the current page, you can use the `refresh` method:

```R
remDr$refresh()
```


## Accessing Elements in the DOM

The DOM stands for *the Document Object Model*. It is a cross-platform and language-independent convention for representing and interacting with objects in `HTML`, `XHTML` and `XML` documents. Interacting with the DOM will be very important for us with Selenium, and the webDriver provides a number of methods in which to do this.

A basic HTML page is:

```html
<!DOCTYPE html>
<html>
<body>

<h1>My First Heading</h1>

<p>My first paragraph.</p>

</body>
</html>
```

The query box on the front page of `http://www.google.com` has html code `<input id=..... name="q" ...</input>` associated with it. The full html associated with the input tag is:

```html
<input spellcheck="false" dir="ltr" style="border: medium none; padding: 0px; margin: 0px; height: auto; width: 100%; background: transparent url(&quot;data:image/gif;base64,R0lGODlhAQABAID/AMDAwAAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw%3D%3D&quot;) repeat scroll 0% 0%; position: absolute; z-index: 6; left: 0px; outline: medium none;" aria-autocomplete="both" role="combobox" aria-haspopup="false" class="gsfi" id="lst-ib" maxlength="2048" name="q" autocomplete="off" title="Search" value="" aria-label="Search" type="text">
```

***NOTE:*** The above HTML is very liable to change however the input node has had an attribute name = q for sometime so we can mostly rely on this.

### Search by Name

To find this element in the DOM, a number of methods can be used. We can search by the name:

```R
remDr$navigate("http://www.google.com/ncr")
webElem <- remDr$findElement(using = "name", value = "q")
webElem$getElementAttribute("name")
```

```
## [[1]]
## [1] "q"
```

```R
webElem$getElementAttribute("class")
```

```
## [[1]]
## [1] "gsfi lst-d-f"
```

```R
webElem$getElementAttribute("id")
```

```
## [[1]]
## [1] "lst-ib"
```

### Search by ID

In HTML, the ID of a DOM element should be unique, so this is usually a good locator to use. As noted above, the Google ID of the query box may change (the one we see is "lst-ib", so we use that). You may see an alternative ID.

```R
remDr$navigate("http://www.google.com/ncr")
webElem <- remDr$findElement(using = "id", value = "lst-ib")
```

#### Highlight an Element

An element that is visible in the DOM can usually be highlighted. Using our `webElem` which is an object of class `webElement`, we can use the associated `highlightElement` method to visually indicate to us we have the correct element. Try it:

```R
webElem$highlightElement()
```

You should see the query box flashing black and yellow to indicate visually the DOM element you have selected.

### Search by Class

We can also search by the class name:

```R
remDr$navigate("http://www.google.com/ncr")
webElem <- remDr$findElement(using = "class", "gsfi")
webElem$getElementAttribute("class")
```

```
## [[1]]
## [1] "gsfi lst-d-f"
```

```R
webElem$getElementAttribute("type")
```

```
## [[1]]
## [1] "text"
```

***NOTE:*** the class is listed as "gsfi lst-d-f", and we searched using "gsfi". This is an example of a compound class. The class name selector can only be used for elements with a single class ("gsfi lst-d-f" indicates two classes "gsfi" and "lst-d-f"). For more complicated select queries, we therefore use CSS or xpath instead.

### Search Using CSS Selectors

To replicate our name search using css selectors, we could use:

```R
webElem <- remDr$findElement(using = "css", "input[name='q']")
# OR
webElem2 <- remDr$findElement(using = "css", "[name='q']")
```

We can see we get the same element (using the `compareElements` method of the `webElement` class):

```R
webElem$compareElements(webElem2)
```

```
## [[1]]
## [1] TRUE
```

and to search on ID using the CSS Selectors (again, the ID you see maybe different):

```R
webElem <- remDr$findElement(using = "css", "input#lst-ib")
webElem$getElementAttribute("name")
```

```
## [[1]]
## [1] "q"
```

and class:

```R
webElem <- remDr$findElement(using = "css", "[class = 'gsfi lst-d-f']")
```

***NOTE:*** no issue with compound class names with CSS

A good example of searching using css-selectors is given [here](http://saucelabs.com/resources/selenium/css-selectors).

### Search Using XPath

The final option is to search using XPath. Normally, one would use XPath by default when searching or CSS. Both are the go-to options. 

XPath using ID:

```R
webElem <- remDr$findElement(using = "xpath", "//input[@id = 'lst-ib']")
```

Xpath using class:

```R
webElem <- remDr$findElement(using = "xpath", "//input[@class = 'gsfi lst-d-f']")
```

***NOTE:*** with XPath, we have no issues using a compound class name.


## Sending Events to Elements

To illustrate how to interact with elements, we will again use the `http://www.google.com/ncr` as an example.

### Sending Text to Elements

Suppose we would like to search for "R cran" on google. We would need to find the element for the query box and send the appropriate text to it. We can do this using the `sendKeysToElement` method for the `webElement` class.

```R
remDr$navigate("http://www.google.com/ncr")
webElem <- remDr$findElement(using = "css", "[name = 'q']")
webElem$sendKeysToElement(list("R Cran"))
```

### Sending Key Presses to Elements

We should see that the text "R Cran" has now been entered into the query box.

How do we press enter. We can simply send the enter key to query box. The enter key would be denoted as `"\uE007"`(its UTF-8 code). So we could use:

```R
remDr$navigate("http://www.google.com/ncr")
webElem <- remDr$findElement(using = "css", "[name = 'q']")
webElem$sendKeysToElement(list("R Cran", "\uE007"))
```

It is not very easy to remember UTF-8 codes for appropriate keys, so a mapping has been provided in `RSelenium`. `?selkeys` will bring up a help page explaining the mapping. The UTF-8 codes have been mapped to easy to remember names.

To use `selkeys`, we would send the following:

```R
remDr$navigate("http://www.google.com/ncr")
webElem <- remDr$findElement(using = "css", "[name = 'q']")
webElem$sendKeysToElement(list("R Cran", key = "enter"))
```

Typing `selKeys` into the console will bring up the list of mappings.

### Sending Mouse Events to Elements

For this example, we will go back to the Google front page and search for "R Cran". then, we will click the link for "The Comprehensive R Archive Network".

```R
remDr$navigate("http://www.google.com/ncr")
webElem <- remDr$findElement(using = "css", "[name = 'q']")
webElem$sendKeysToElement(list("R Cran", key = "enter"))
```

The header for each link is contained in a `<h3 class = "r">` tag. We will access the "h3" headers first. It will be succinct to find these elements using "css selectors".

```R
webElems <- remDr$findElements(using = "css selector", "h3.r")
resHeaders <- unlist(lapply(webElems, function(x) {x$getElementText()}))
resHeaders
```

```
## [1] "The Comprehensive R Archive Network"                                   
## [2] "Download R-3.3.2 for Windows. The R-project for statistical ... - CRAN"
## [3] "About Microsoft R Open: The Enhanced R Distribution . MRAN"            
## [4] "R (programming language) - Wikipedia"                                  
## [5] "R-Cran - StatLib - Carnegie Mellon University"                         
## [6] "Submitting your first package to CRAN, my experience | R-bloggers"     
## [7] "Debian -- Package Search Results -- r-cran"                            
## [8] "It's crantastic!"                                                      
## [9] "METACRAN"                                                              
## [10] "CRAN - Package PopGenReport"   
```

***NOTE:*** this is how the headers were presented at time of writing. Class names etc. are liable to change.

We can see that the first link is the one we want, but in case Google's search results change, we refer to it as

```R
webElem <- webElems[[which(resHeaders == "The Comprehensive R Archive Network")]]
```

How do we click the link? We can use the `clickElement` method:

```R
webElem$clickElement()
remDr$getCurrentUrl()
```

```
## [[1]]
## [1] "https://cran.r-project.org/"
```

```R
remDr$getTitle()
```

```
## [[1]]
## [1] "The Comprehensive R Archive Network"
```


## Injecting JavaScript

Sometimes it is necessary to interact with the current url using JavaScript. This maybe necessary to call bespoke methods or to have more control over the page for example by adding the `JQuery` library to the page if it is missing.

`RSelenium` has two methods we can use to execute JavaScript namely `executeScript` and `executeAsyncScript` from the `remoteDriver` class. We return to the Google front page to investigate these methods.

### Injecting JavaScript Synchronously

Returning to the Google front page, we can find the element for the "Google" image. The image has `id = "hplogo"`, and we can use this in an XPath or search by ID to select the element. In this case, we use "css selectors":

```R
remDr$navigate("http://www.google.com/ncr")
webElem <- remDr$findElement("css", "img#hplogo")
```

Is the image visible? Clearly, it is, but we can check using JavaScript:

```R
script <- "return document.getElementById('hplogo').hidden;"
remDr$executeScript(script, args = list())
```

```
## [[1]]
## [1] FALSE
```

Great! So the image is not hidden indicated by the `FALSE`. We can hide it executing some simple JavaScript.

```R
script <- "document.getElementById('hplogo').hidden = true; 
           return document.getElementById('hplogo').hidden;"
remDr$executeScript(script, args = list())
```

```
## [[1]]
## [1] TRUE

```

So now the image is hidden. We used an element here given by `id = "hplogo"`. We had to use the JavaScript function `getElementById` to access it. It would be nicer if we could have used `webElem` which we had specified earlier.

If we pass a `webElement` object as an argument to either `executeScript` or `executeAsyncScript`, `RSelenium` will pass it in an appropriate fashion.

```R
script <- "arguments[0].hidden = false; return arguments[0].hidden;"
remDr$executeScript(script, args = list(webElem))
```

```
## [[1]]
## [1] FALSE
```

Notice that how we passed the `webElement` object to the method `executeScript`. The script argument defines the script to execute in the form of a function body. The value returned by that function will be returned to the client. The function will be invoked with the provided arguments. If the function returns an element, this will be returned as an object of `webElement` class:

```R
script <- "return document.getElementsByName('q');"
test <- remDr$executeScript(script, args = list())
test[[1]]
```

```
## [1] "remoteDriver fields"
## $remoteServerAddr
## [1] "localhost"
## 
## $port
## [1] 4445
## 
## $browserName
## [1] "firefox"
## 
## $version
## [1] ""
## 
## $platform
## [1] "ANY"
## 
## $javascript
## [1] TRUE
## 
## $nativeEvents
## [1] TRUE
## 
## $extraCapabilities
## list()
## 
## [1] "webElement fields"
## $elementId
## [1] "21"
```

```R
class(test[[1]])
```

```
## [1] "webElement"
## attr(,"package")
## [1] "RSelenium"
```

Try to highlight the element as before:

```R
test[[1]]$highlightElement()
```

### Injecting JavaScript Asynchronously

I will briefly touch on async versus sync calls here. Firstly, we need to set an appropriate asynchronous timeout (that is longer than the async operation we are likely to carryout, but it will ensure that we will eventually error out in case of an issue).

```R
remDr$navigate("http://www.google.com/ncr")
remDr$setAsyncScriptTimeout(10000)
```

Observe:

```R
webElem <- remDr$findElement("css", "img#hplogo")
script <- "
cb = arguments[arguments.length -1];
webElem = arguments[0];
setTimeout(function(){webElem.hidden = true; cb('DONE');},5000);"
remDr$executeAsyncScript(script, args = list(webElem))
```

Versus

```R
remDr$navigate("http://www.google.com/ncr")
webElem <- remDr$findElement("css", "img#hplogo")
script <- "
webElem = arguments[0];
setTimeout(function(){webElem.hidden = true;},5000);
return 'DONE';
"
remDr$executeScript(script, args = list(webElem))
```

The async version waits until the callback (defined as the last argument `arguments[arguments.length - 1]` as JavaScript is zero-indexed) is called whereas the sync version returns straight away. In both cases, the Google logo disappears after five seconds.


## Frames and Windows

In the context of a web browser, a frame is a part of a web page or browser window which displays content independent of its container, with the ability to load content independently.

### Frames in Selenium

We will demonstrate interacting with frames by way of example. [The Comprehensive R Archive Network](https://CRAN.R-project.org/) conveniently contains frames so we shall use `RSelenium` to interact with it.

Assume that we have a remoteDriver opened:

```R
remDr$navigate("https://CRAN.r-project.org/")
XML::htmlParse(remDr$getPageSource()[[1]])
```

```html
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Frameset//EN" "http://www.w3.org/TR/html4/frameset.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title>The Comprehensive R Archive Network</title>
<meta content="text/html; charset=utf-8" http-equiv="content-type">
<link type="image/x-icon" href="favicon.ico" rel="icon">
<link type="image/x-icon" href="favicon.ico" rel="shortcut icon">
<link href="R.css" type="text/css" rel="stylesheet">
</head>
<frameset style="border: none;" cols="1*, 4*">
<frameset rows="120, 1*">
<frame frameborder="0" name="logo" src="logo.html">
<frame frameborder="0" name="contents" src="navbar.html">
</frameset>
<frame frameborder="0" name="banner" src="banner.shtml">
<noframes>
&lt;h1&gt;The Comprehensive R Archive Network&lt;/h1&gt;

Your browser seems not to support frames,
here is the &lt;A href="navbar.html"&gt;contents page&lt;/A&gt; of CRAN.
</noframes>
</frameset>
</html>
```

We can see that the content is contained in three frames, and we don't appear to have access to the content within a frame, but in the browser, we see all the content:

```R
remDr$maxWindowSize()
remDr$screenshot(display = TRUE)
```

<h6 align = center>RProject front page</h6>
<img src="https://res.cloudinary.com/rselenium/image/upload/v1537203294/RProject.png" title = "RProject front page on linux firefox 26.0" width = '100%'/>

To access the content, we have to switch to a frame using the `switchToFrame` method of the `remoteDriver` class.

```R
webElems <- remDr$findElements(using = "tag name", "frame")
# webElems <- remDr$findElements(value = "//frame") # using xpath
# webElems <- remDr$findElements("css", value = "frame") # using css

sapply(webElems, function(x){x$getElementAttribute("src")})
```

```
## [[1]]
## [1] "https://cran.r-project.org/logo.html"
## 
## [[2]]
## [1] "https://cran.r-project.org/navbar.html"
## 
## [[3]]
## [1] "https://cran.r-project.org/banner.shtml"
```

```R
remDr$switchToFrame(webElems[[2]])
XML::htmlParse(remDr$getPageSource()[[1]])
```

```html
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title>R Contents</title>
<meta content="text/html; charset=utf-8" http-equiv="content-type">
<link href="R.css" type="text/css" rel="stylesheet">
</head>
<body>

<em class="navigation">CRAN</em><br><a target="banner" href="mirrors.html">Mirrors</a><br><a target="banner" href="//www.R-project.org/news.html">What's new?</a><br><a target="banner" href="web/views/">Task Views</a><br><a target="banner" href="search.html">Search</a><br><!--<a href="pkg_submit.html" target="_top">Submit</a><BR>--><p>
<em class="navigation">About R</em><br><a target="_top" href="//www.R-project.org/">R Homepage</a><br><a target="_top" href="http://journal.R-project.org/">The R Journal</a>

</p>
<p>
<em class="navigation">Software</em><br><a target="banner" href="sources.html">R Sources</a><br><a target="banner" href="bin/">R Binaries</a><br><a target="banner" href="web/packages/">Packages</a><br><a target="banner" href="other-software.html">Other</a>

</p>
<p>
<em class="navigation">Documentation</em><br><a target="banner" href="manuals.html">Manuals</a><br><a target="banner" href="faqs.html">FAQs</a><br><a target="banner" href="other-docs.html">Contributed</a><br></p>
</body>
</html>
```

Now, we see the source code of the navigation side panel. Notice that how we used a `webElement` in the method `switchToFrame`. To further demonstrate, we are now "in" this frame. Let's get all the "href" attributes:

```R
webElems <- remDr$findElements(using = "css", "[href]")
unlist(sapply(webElems, function(x) {x$getElementAttribute("href")}))
```

```
## [1] "https://cran.r-project.org/R.css"              
## [2] "https://cran.r-project.org/mirrors.html"       
## [3] "https://www.r-project.org/news.html"           
## [4] "https://cran.r-project.org/web/views/"         
## [5] "https://cran.r-project.org/search.html"        
## [6] "https://www.r-project.org/"                    
## [7] "http://journal.r-project.org/"                 
## [8] "https://cran.r-project.org/sources.html"       
## [9] "https://cran.r-project.org/bin/"               
## [10] "https://cran.r-project.org/web/packages/"      
## [11] "https://cran.r-project.org/other-software.html"
## [12] "https://cran.r-project.org/manuals.html"       
## [13] "https://cran.r-project.org/faqs.html"          
## [14] "https://cran.r-project.org/other-docs.html"    
```

Notice that if we pass a `NULL` value to the method `switchToFrame`, we move back to the default view.

```R
remDr$switchToFrame(NULL)
XML::htmlParse(remDr$getPageSource()[[1]])
```

```html
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Frameset//EN" "http://www.w3.org/TR/html4/frameset.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title>The Comprehensive R Archive Network</title>
<meta content="text/html; charset=utf-8" http-equiv="content-type">
<link type="image/x-icon" href="favicon.ico" rel="icon">
<link type="image/x-icon" href="favicon.ico" rel="shortcut icon">
<link href="R.css" type="text/css" rel="stylesheet">
</head>
<frameset style="border: none;" cols="1*, 4*">
<frameset rows="120, 1*">
<frame frameborder="0" name="logo" src="logo.html">
<frame frameborder="0" name="contents" src="navbar.html">
</frameset>
<frame frameborder="0" name="banner" src="banner.shtml">
<noframes>
&lt;h1&gt;The Comprehensive R Archive Network&lt;/h1&gt;

Your browser seems not to support frames,
here is the &lt;A href="navbar.html"&gt;contents page&lt;/A&gt; of CRAN.
</noframes>
</frameset>
<body><canvas id="fxdriver-screenshot-canvas" style="display: none;" width="1360" height="559"></canvas></body>
</html>
```

Finally we can switch to the main panel using a name

```R
remDr$switchToFrame("banner")
XML::htmlParse(remDr$getPageSource()[[1]])
```

```html
<!DOCTYPE html PUBLIC "-//IETF//DTD HTML//EN">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title>The Comprehensive R Archive Network</title>
<link href="R.css" type="text/css" rel="stylesheet">
</head>
<body>

<h1>The Comprehensive R Archive Network</h1>

<div align="center">
<table border="1" width="80%"><tbody>
<tr>
<td>
<h3>Download and Install R</h3>
      Precompiled binary distributions of the base system and
      contributed packages, <strong>Windows and
      Mac</strong> users most likely want one of these versions of R:
	<ul>
<li><a href="bin/linux/">Download R for Linux</a></li>
	  <li>
<a href="bin/macosx/">Download R for (Mac) OS X</a>
	  </li>
<li><a href="bin/windows/">Download R for Windows</a></li>
	</ul>
    R is part of many Linux distributions, you should check with your Linux package management system in addition to the link above.
    </td>
  </tr>
<tr>
<td>
................
................
................
<h3>Questions About R</h3>
      <ul><li>
      If you have questions about R like how to download and install
      the software, or what the license terms are,
      please read our <a href="faqs.html">answers to frequently asked
      questions</a> before you send an email.
	</li></ul>
</td>
  </tr>
</tbody></table>
</div>

<h2>What are R and CRAN?</h2>

<p> R is GNU, a freely available language and environment for
statistical computing and graphics which provides a wide variety of
statistical and graphical techniques: linear and nonlinear modelling,
statistical tests, time series analysis, classification, clustering,
etc. Please consult the <a target="_top" href="https://www.R-project.org/">R project homepage</a> for further information.
</p>

<p> CRAN is a network of ftp and web servers around the world that
store identical, up-to-date, versions of code and documentation for
R. Please use the CRAN <a href="mirrors.html">mirror</a> nearest to you to minimize network
load.
</p>

<h2 id="submitting">Submitting to CRAN </h2>

<p>
To submit a package to CRAN,
check that your submission meets the
<a href="https://CRAN.R-project.org/web/packages/policies.html">CRAN
  Repository Policy</a> and then use the
<a href="https://xmpalantir.wu.ac.at/cransubmit/">web form</a>.
</p>

<p>
If this fails, upload to
<a target="_blank" href="ftp://CRAN.R-project.org/incoming/">ftp://CRAN.R-project.org/incoming/</a>
and send an email to
<a href="mailto:CRAN@R-project.org">CRAN@R-project.org</a> following the policy.
Please do not attach submissions to emails, because this will clutter up
the mailboxes of half a dozen people.
</p>

<p>
Note that we generally do not accept submissions of precompiled
binaries due to security reasons. All binary distribution listed above
are compiled by selected maintainers, who are in charge for all
binaries of their platform, respectively.
</p>

<p>
</p>
<hr>
<!--#if expr="$CRAN_HOST" --><!--#echo  encoding="none" var="CRAN_HOST"--><br><!--#endif -->
</body>
</html>
```

### Windows in Selenium

The easiest way to illustrate Windows in RSelenium is again by way of example. We will use the CRAN website.
First, we select the "download R"" element in the main frame. 

```R
remDr$navigate("https://cran.r-project.org/")
remDr$switchToFrame("banner")
webElems <- remDr$findElements("partial link text", "Download R")

sapply(webElems, function(x) x$getElementText())
```

```
## [[1]]
## [1] "Download R for Linux"
## 
## [[2]]
## [1] "Download R for (Mac) OS X"
## 
## [[3]]
## [1] "Download R for Windows"
```

We now send a selection of key presses to the first element to open the link it points to in a new window. If you did it manually you would move the mouse to the element right click on the link press the down arrow key twice then press enter. We will do the same

```R
loc <- webElems[[1]]$getElementLocation()
loc[c('x','y')]
```

```
## $x
## [1] 158
## 
## $y
## [1] 132
```

```R
remDr$mouseMoveToLocation(webElement = webElems[[1]]) # move mouse to the element we selected
remDr$click(2) # 2 indicates click the right mouse button
remDr$sendKeysToActiveElement(
  list(key = 'down_arrow', key = 'down_arrow', key = 'enter')
)
```

Notice now that a new windows has opened on the remote browser.

```R
remDr$getCurrentWindowHandle()
```

```
## [[1]]
## [1] "{573d17e5-b95a-41b9-a65f-04092b6a804b}"
```

```R
remDr$getWindowHandles()
```

```
## [[1]]
## [1] "{4896393a-c215-4976-b4ca-030d6b75b67d}"
## 
## [[2]]
## [1] "{69c00f18-d3a7-44d7-a236-c6b5e6c264ff}"
```

```R
remDr$getTitle()
```

```
## [[1]]
## [1] "The Comprehensive R Archive Network"
```

```R
currWin <- remDr$getCurrentWindowHandle()
allWins <- unlist(remDr$getWindowHandles())
otherWindow <- allWins[!allWins %in% currWin[[1]]]
remDr$switchToWindow(otherWindow)
remDr$getTitle()
```

```
## [[1]]
## [1] "Index of /bin/linux"
```

So using the code above one can observe how to switch between different windows on the remote browser.


## Appendix

### Server Functions

In previous versions of `RSelenium`, there were two functions `checkForServer` and `startServer` which were used to download and start a Selenium binary. These functions are defunct as users had issue with using them cross-platform. They are still in the RSelenium package and can be accessed in the serverUtils directory `file.path(find.package("RSelenium"), "examples/serverUtils")`

```R
file.path(find.package("RSelenium"), "examples/serverUtils")
```

```
## [1] "/home/john/R/x86_64-pc-linux-gnu-library/3.3/RSelenium/examples/serverUtils"
```

They may work for you depending on your platform and setup, but are not supported.

### rsDriver

The `rsDriver` function is a wrapper for the `selenium` function from the [wdman](https://docs.ropensci.org/wdman/) package. It allows the user to manage the binaries used to run a Selenium Server. It returns an environment containing a `client` and a `server`. 

By default, it runs a Chrome browser. Other browsers such as Firefox, PhantomJS, and Internet Explorer can be selected using the `browser` argument.

The default port a Selenium Server is run on using the `rsDriver` function is 4567L.

```R
rD <- rsDriver(verbose = FALSE)
rD
```

```
## $client
##   browserName                                   id
## 1      chrome 1670dcc6-3c97-4717-bf28-4d1fc8eea2c1
## 
## $server
## Process Handle
## command   : /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Dwebdriver.chrome.driver=/home/john/.local/share/binman_chromedriver/linux64/2.27/chromedriver -Dwebdriver.gecko.driver=/home/john/.local/share/binman_geckodriver/linux64/0.13.0/geckodriver -Dphantomjs.binary.path=/home/john/.local/share/binman_phantomjs/linux64/2.1.1/phantomjs-2.1.1-linux-x86_64/bin/phantomjs -jar /home/john/.local/share/binman_seleniumserver/generic/3.0.1/selenium-server-standalone-3.0.1.jar -port 4567
## system id : 120670
## state     : running
```

Assign the client to a new variable and drive the client to navigate to a page:

```R
remDr <- rD$client

remDr$navigate("http://www.r-project.org")
remDr$getTitle()
```

```
## [[1]]
## [1] "R: The R Project for Statistical Computing"
```

Check the server logs:

```R
rD$server$log()
```

```
#> $stderr
#>  [1] "10:29:11.536 INFO - Selenium build info: version: '3.0.1', revision: '1969d75'"                                                                                          
#>  [2] "10:29:11.537 INFO - Launching a standalone Selenium Server"                                                                                                              
#>  [3] "2017-01-18 10:29:11.551:INFO::main: Logging initialized @179ms"                                                                                                          
#>  [4] "10:29:11.594 INFO - Driver provider org.openqa.selenium.ie.InternetExplorerDriver registration is skipped:"                                                              
#>  [5] " registration capabilities Capabilities [{ensureCleanSession=true, browserName=internet explorer, version=, platform=WINDOWS}] does not match the current platform LINUX"
#>  [6] "10:29:11.595 INFO - Driver provider org.openqa.selenium.edge.EdgeDriver registration is skipped:"                                                                        
#>  [7] " registration capabilities Capabilities [{browserName=MicrosoftEdge, version=, platform=WINDOWS}] does not match the current platform LINUX"                             
#>  [8] "10:29:11.595 INFO - Driver class not found: com.opera.core.systems.OperaDriver"                                                                                          
#>  [9] "10:29:11.595 INFO - Driver provider com.opera.core.systems.OperaDriver registration is skipped:"                                                                         
#> [10] "Unable to create new instances on this machine."                                                                                                                         
#> [11] "10:29:11.595 INFO - Driver class not found: com.opera.core.systems.OperaDriver"                                                                                          
#> [12] "10:29:11.595 INFO - Driver provider com.opera.core.systems.OperaDriver is not registered"                                                                                
#> [13] "10:29:11.596 INFO - Driver provider org.openqa.selenium.safari.SafariDriver registration is skipped:"                                                                    
#> [14] " registration capabilities Capabilities [{browserName=safari, version=, platform=MAC}] does not match the current platform LINUX"                                        
#> [15] "2017-01-18 10:29:11.628:INFO:osjs.Server:main: jetty-9.2.15.v20160210"                                                                                                   
#> [16] "2017-01-18 10:29:11.647:INFO:osjsh.ContextHandler:main: Started o.s.j.s.ServletContextHandler@2ef5e5e3{/,null,AVAILABLE}"                                                
#> [17] "2017-01-18 10:29:11.658:INFO:osjs.ServerConnector:main: Started ServerConnector@724af044{HTTP/1.1}{0.0.0.0:4567}"                                                        
#> [18] "2017-01-18 10:29:11.658:INFO:osjs.Server:main: Started @285ms"                                                                                                           
#> [19] "10:29:11.658 INFO - Selenium Server is up and running"                                                                                                                   
#> [20] "10:29:12.295 INFO - SessionCleaner initialized with insideBrowserTimeout 0 and clientGoneTimeout 1800000 polling every 180000"                                           
#> [21] "10:29:12.319 INFO - Executing: [new session: Capabilities [{nativeEvents=true, browserName=chrome, javascriptEnabled=true, version=, platform=ANY}]])"                   
#> [22] "10:29:12.331 INFO - Creating a new session for Capabilities [{nativeEvents=true, browserName=chrome, javascriptEnabled=true, version=, platform=ANY}]"                   
#> [23] "Starting ChromeDriver 2.27.440175 (9bc1d90b8bfa4dd181fbbf769a5eb5e575574320) on port 12090"                                                                              
#> [24] "Only local connections are allowed."                                                                                                                                     
#> [25] "10:29:12.556 INFO - Attempting bi-dialect session, assuming Postel's Law holds true on the remote end"                                                                   
#> [26] "10:29:13.115 INFO - Detected dialect: OSS"                                                                                                                               
#> [27] "10:29:13.144 INFO - Done: [new session: Capabilities [{nativeEvents=true, browserName=chrome, javascriptEnabled=true, version=, platform=ANY}]]"                         
#> [28] "10:29:13.172 INFO - Executing: [get: http://www.r-project.org])"                                                                                                         
#> [29] "10:29:14.340 INFO - Done: [get: http://www.r-project.org]"                                                                                                               
#> [30] "10:29:14.356 INFO - Executing: [get title])"                                                                                                                             
#> [31] "10:29:14.363 INFO - Done: [get title]"                                                                                                                                   
#> 
#> $stdout
#> character(0)
```

Close the client and the server:

```R
remDr$close()
rD$server$stop()
```

```
## [1] TRUE
```

Check the server status:

```R
rD$server$process
```

```
## Process Handle
## command   : /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -Dwebdriver.chrome.driver=/home/john/.local/share/binman_chromedriver/linux64/2.27/chromedriver -Dwebdriver.gecko.driver=/home/john/.local/share/binman_geckodriver/linux64/0.13.0/geckodriver -Dphantomjs.binary.path=/home/john/.local/share/binman_phantomjs/linux64/2.1.1/phantomjs-2.1.1-linux-x86_64/bin/phantomjs -jar /home/john/.local/share/binman_seleniumserver/generic/3.0.1/selenium-server-standalone-3.0.1.jar -port 4567
## system id : 120670
## state     : terminated
```

For further detail and any issues you may have with the `rsDriver` function, please see the [wdman](https://docs.ropensci.org/wdman/) project.
---
title: "Headless Browsing"
output:
  html_vignette:
    toc: yes
vignette: >
  %\VignetteIndexEntry{Headless Browsing with RSelenium}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## PhantomJS

`PhantomJS` is a headless WebKit scriptable with a JavaScript API. It has fast and native support for various web standards: DOM handling, CSS selector, JSON, Canvas, and SVG. `RSelenium` can drive `PhantomJS` using two methods: directly or via the standalone Selenium Server. 

### Driving PhantomJS Directly

The `PhantomJS` binary can be driven directly with `RSelenium`. `PhantomJS` needs to be started in webdriver mode then `RSelenium` can **communicate with it directly without the need for Selenium Server.** The command line options for `PhantomJS` are outlined at [http://phantomjs.org/api/command-line.html](http://phantomjs.org/api/command-line.html). We note that it is necessary to start `PhantomJS` with the `--webdriver` option and an optional IP/port. `RSelenium` as of `v1.3.2` has a utility function `phantom` that will handle starting the `PhantomJS` binary in webdriver mode by default on port 4444. So to drive `PhantomJS` sans Selenium Server can be done as follows:

```R
require(RSelenium)
pJS <- phantom()
Sys.sleep(5) # give the binary a moment
remDr <- remoteDriver(browserName = 'phantomjs')
remDr$open()
remDr$navigate("http://www.google.com/ncr")
remDr$getTitle()[[1]] # [1] "Google"
remDr$close
pJS$stop() # close the PhantomJS process, note we dont call remDr$closeServer()
```

### Driving PhantomJS Using Selenium Server

For completeness we outline the process of opening a `PhantomJS` browser using selenium server. It is assumed that the `PhantomJS` binary is in the users path.

```R
require(RSelenium)
RSelenium::startServer()
remDr <- remoteDriver(browserName = "phantomjs")
remDr$open()
remDr$navigate("http://www.google.com/ncr")
remDr$close()
remDr$closeServer()
```

#### Providing the PhantomJS Path

It may not be possible for a user to have the `PhantomJS` binary in their path. In this case
a user may pass the path of the `PhantomJS` binary to Selenium Server:

```R
require(RSelenium)
RSelenium::startServer()
eCap <- list(phantomjs.binary.path = "C:/Users/john/Desktop/phantomjs.exe")
remDr <- remoteDriver(browserName = "phantomjs", extraCapabilities = eCap)
remDr$open()
....
```

So in the above example I suppose the `PhantomJS` binary has been moved to my Desktop which we assume is not in my path. An extra capability `phantomjs.binary.path` detailed [https://github.com/detro/ghostdriver](https://github.com/detro/ghostdriver) can be used to provide the path to `PhantomJS` to Selenium Server.

### Additional PhantomJS Capabilities

#### Setting a User Agent

A user agent can be set using the `phantomjs.page.settings.userAgent` capability. 

```R
pJS <- phantom()
Sys.sleep(5)
remDr <- remoteDriver(browserName = "phantomjs")
remDr$open()
remDr$navigate("http://www.whatsmyuseragent.com/")
remDr$findElement("id", "userAgent")$getElementText()[[1]]
```

```
## [1] "Your User Agent String is:\nMozilla/5.0 (Unknown; Linux x86_64)
## AppleWebKit/534.34 (KHTML, like Gecko) PhantomJS/1.9.7 Safari/534.34"
```

```R
remDr$close()
eCap <- list(
  phantomjs.page.settings.userAgent = "Mozilla/5.0 (Windows NT 6.1; WOW64; rv:29.0) Gecko/20120101 Firefox/29.0"
)
remDr <- remoteDriver(browserName = "phantomjs", extraCapabilities = eCap)
remDr$open()
remDr$navigate("http://www.whatsmyuseragent.com/")
remDr$findElement("id", "userAgent")$getElementText()[[1]]
```

```
## [1] "Your User Agent String is:\nMozilla/5.0 (Windows NT 6.1; WOW64; rv:29.0)
## Gecko/20120101 Firefox/29.0"
```

```R
remDr$close()
pJS$stop()
```

The https://github.com/ariya/phantomjs/wiki/API-Reference-WebPage#webpage-settings
In the above example it can be seen that the default useragent identifies us as `PhantomJS`. Some web content maybe inaccessible or blocked for `PhantomJS` users. Here we demonstrate changing our user agent so the website sees us as `Firefox 29.0`.

#### Other Possible Options

The general form of specifying [PhantomJS internal page objects](https://github.com/ariya/phantomjs/wiki/API-Reference-WebPage#webpage-settings) take the form `phantomjs.page.settings.SETTING = VALUE` where `SETTING` is the appropriate PhantomJS internal page object.
As an example we inhibit the loading of inline images:

```R
require(RSelenium)
pJS <- phantom()
Sys.sleep(5)
eCap <- list(phantomjs.page.settings.loadImages = FALSE)
remDr <- remoteDriver(browserName = "phantomjs", extraCapabilities = eCap)
remDr$open()
remDr$navigate("http://www.google.com/ncr")
remDr$screenshot(display = TRUE)
remDr$close()
pJS$stop()
```

We can see that the images are not loaded:
<img src="https://res.cloudinary.com/johndharrison/image/upload/v1497012298/RSelenium/headless/Screenshot_2014-06-03_14.32.18.png"  title = "PhantomJS loadImages = FALSE" style = "margin-left: auto;margin-right: auto; display: block;"   width = '100%'/>

## X Virtual Frame Buffer

For the discussion on `xvfb` and the related VPS, I refer you to [this blog entry](http://johndharrison.blogspot.com/2014/03/rstudioshiny-server-on-digital-ocean.html). How to setup a VPS with rstudio server and shiny server etc. is outlined. 



### Setup

The VPS i am connecting to has an ip of `128.199.255.233`. I have rstudio server running on port 8787. On the remote server we observe

```R
library(RSelenium)
RSelenium::startServer()
Sys.which('phantomjs')
```

```
                 phantomjs
"/usr/local/bin/phantomjs"
```

```R
Sys.which('firefox')
```

```
firefox
     ""
```

```R
Sys.which('chrome')
```

```
chrome
    ""
```

So  we have started a selenium server running on (default) port 4444. Firefox and google chrome are not currently installed on this remote machine. Lets install firefox first. On the remote VPS we run 

```sh
sudo apt-get install firefox
```

Now checking in the remote rstudio 

```R
Sys.which('firefox')
```

```
##            firefox 
## "/usr/bin/firefox" 
```

If we try now to connect to the remote server and open firefox:

```R
remDr <- remoteDriver(remoteServerAddr = "128.199.255.233")
remDr$open()
```

```
## [1] "Connecting to remote server"
## Error:    Summary: UnknownError
##    Detail: An unknown server-side error occurred while processing the command.
##    class: org.openqa.selenium.WebDriverException
```

We can see the problem if we try to run firefox in the remote shell:

<img src="https://res.cloudinary.com/johndharrison/image/upload/v1497012290/RSelenium/headless/firefox.png"  title = "PhantomJS loadImages = FALSE" style = "margin-left: auto;margin-right: auto; display: block;"  width = '100%'/>

Firefox is install but there is no display on our headless VPS. We can use [xvfb](http://www.x.org/archive/X11R7.7/doc/man/man1/Xvfb.1.xhtml) to provide a virtual display for our browser to run in. 

```
Xvfb :0 -screen 0 1024x768x24 2>&1 >/dev/null &
export DISPLAY=:0
nohup xvfb-run java -jar selenium-server-standalone.jar > selenium.log &
```


## PhantomJS API Examples

The `phantomExecute` method of the `remoteDriver` class allows the user to interact with the `PhantomJS` API. Currently the method only works for direct calls to `PhantomJS` using the `phantom` utility function. Driving `PhantomJS` through the `Selenium` Server and calling the `phantomExecute` method currently doesn't function and is an open issue (in the ghostDriver project). In the following sections we outline examples of using the `PhantomJS` API.

### Interacting with the Console

The `PhantomJS` API implements a number of callbacks which can be defined. [onLoadFinished](http://phantomjs.org/api/webpage/handler/on-load-finished.html) is one such callback. This callback is invoked when the page finishes the loading. It may accept a single argument indicating the pages status: `success` if no network errors occurred, otherwise `fail`.

We give a simple example of writing to the console log when a page is loaded. 

```R
library(RSelenium)
pJS <- phantom()
remDr <- remoteDriver(browserName = "phantom")
remDr$open()
result <- remDr$phantomExecute("var page = this;
                                page.onLoadFinished = function(status) {
                                var url = page.url;
                                console.log(\"Status:  \" + status);
                                console.log(\"Loaded:  \" + url);
                               };")
remDr$navigate("http://www.google.com/ncr")
```

```
## Status:  success
## Loaded:  http://www.google.com/
```

```R
remDr$navigate("http://www.bbc.co.uk")
```

```
## Status:  success
## Loaded:  http://www.bbc.co.uk/
```

```R
remDr$navigate("http://www.bbc.com")
```

```
## Status:  success
## Loaded:  http://www.bbc.com/
```

```R
pJS$stop()
```

It can be seen that the callback persists across page calls.

### PhantomJS Writing to File

The next example demonstrates writing to file from `PhantomJS`. Once again the `onLoadFinished` callback is utilised. In this example the html source of the page that is navigated to is downloaded to `output.htm` relative to `getwd()`. An example is given of using `phantom.exit()` to close `PhantomJS` from the API.

```R
library(RSelenium)
pJS <- phantom()
remDr <- remoteDriver(browserName = "phantom")
remDr$open()
result <- remDr$phantomExecute("var page = this;
                                var fs = require(\"fs\");
                                page.onLoadFinished = function(status) {
                                var file = fs.open(\"output.htm\", \"w\");
                                file.write(page.content);
                                file.close();
                                phantom.exit();
                               };")
remDr$navigate("http://www.google.com/ncr")
htmlParse("output.htm")['//title/text()'][[1]]
```

```
## Google
```

```R
pJS$stop()
```

### Injecting a Library into PhantomJS

Next we look at [includeJs](http://phantomjs.org/api/webpage/method/include-js.html).

This includes an external script from the specified url (usually a remote location) on the page and executes the callback upon completion. The library we shall include is `JQuery` using the google CDN. Now any page we call with `PhantomJS` will have the `JQuery` library loaded after the page has finished loading.

```R
library(RSelenium)
pJS <- phantom()
remDr <- remoteDriver(browserName = "phantom")
remDr$open()
remDr$navigate("http://www.google.com/ncr")
# check if the JQuery library is loaded
remDr$executeScript("return window.jQuery == undefined;")[[1]]
# TRUE is returned indicating JQuery is not present
result <- remDr$phantomExecute("var page = this;
                                page.onLoadFinished = function(status) {
                                 var url = page.url;
                                 var jURL = 'http://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js';
                                 console.log(\"Status:  \" + status);
                                 console.log(\"Loaded:  \" + url);
                                 page.includeJs(jURL, function() {console.log(\"Loaded jQuery!\");})
                                };"
                               )
remDr$navigate("http://www.google.com/ncr")
```

```
## Status:  success
## Loaded:  http://www.google.com/
## Loaded jQuery!
```

```R
remDr$executeScript("return window.jQuery == undefined;")[[1]]
# FALSE is returned indicating that JQuery is present
webElem <- remDr$executeScript("return $(\"[name='q']\").get(0);")[[1]]
webElem$sendKeysToElement(list("PhantomJS was here"))
remDr$screenshot(display = TRUE)
pJS$stop()
```

<img src="https://res.cloudinary.com/johndharrison/image/upload/v1497012295/RSelenium/headless/googlePhantomJS.png"  title = "PhantomJS with Jquery injected" style = "margin-left: auto;margin-right: auto; display: block;"   width = '100%'/>


### Starting a PhantomJS Web Server

`PhantomJS` has the ability to act as a [Web Server](http://phantomjs.org/api/webserver/). Here we demonstrate setting `PhantomJS` up as a web server on the localhost on port `8080`. When a user browses to `http://localhost:8080` they are returned a list of the current blog titles on [http://www.r-bloggers.com](http://www.r-bloggers.com). The `Jquery` library is also injected to aid extraction of the blog titles.

```R
pJS <- phantom()
remDr <- remoteDriver(browserName = "phantom")
remDr$open()
"
var server = require('webserver').create();
server.listen(8080, function (request, response) {
  var page = new WebPage();
  page.open('http://www.r-bloggers.com/', function (status) {
    var jURL = 'http://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js';
    page.includeJs(jURL, function() {
      console.log(\"Loaded jQuery!\");
      var blogs = page.evaluate(function () {
        res = $('#mainwrapper .post a[rel=\"bookmark\"]');
        return res.map(function(){return this.innerHTML}).toArray().join('\\n');
      });
      response.statusCode = 200;
      response.write('Current blogs on r-bloggers:\\n');
      response.write(blogs);
      response.write('\\n');
      response.close();
      page.close();
    });
  });
});" -> wsScript

remDr$phantomExecute(wsScript)

head(readLines("http://localhost:8080/"))
```

```
## Loaded jQuery!
## [1] "Current blogs on r-bloggers:"                        "Specifying complicated groups of time series in hts"
## [3] "Creating Inset Map with ggplot2"                     "R and Vertica"                                      
## [5] "RGolf: NGSL Scrabble"                                "European talks. June-July 2014"
```

<img src="https://res.cloudinary.com/johndharrison/image/upload/v1497012296/RSelenium/headless/phantomWebserver.png"  title = "PhantomJS as a Web Server" style = "margin-left: auto;margin-right: auto; display: block;"   width = '100%'/>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selKeys-data.R
\docType{data}
\name{selKeys}
\alias{selKeys}
\title{Selenium key mappings}
\format{A named list. The names are the descriptions of the keys. The
   values are the "UTF-8" character representations.}
\source{
https://github.com/SeleniumHQ/selenium/wiki/JsonWireProtocol#sessionsessionidelementidvalue
}
\usage{
selKeys
}
\description{
This data set contains a list of selenium key mappings.
selKeys is used when a sendKeys variable is needed.
sendKeys is defined as a list.
If an entry is needed from selKeys it is denoted by key.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{checkForServer}
\alias{checkForServer}
\title{Check for Server binary}
\usage{
checkForServer(dir = NULL, update = FALSE, rename = TRUE, beta = FALSE)
}
\arguments{
\item{dir}{A directory in which the binary is to be placed.}

\item{update}{A boolean indicating whether to update the binary if it
is present.}

\item{rename}{A boolean indicating whether to rename to
"selenium-server-standalone.jar".}

\item{beta}{A boolean indicating whether to include beta releases.}
}
\description{
Defunct. Please use \code{\link{rsDriver}}
}
\details{
\code{checkForServer}
A utility function to check if the Selenium Server stanalone binary is
   present.
}
\section{Detail}{
 The downloads for the Selenium project can be found at
   http://selenium-release.storage.googleapis.com/index.html. This
   convenient function downloads the standalone server and places it in
   the RSelenium package directory bin folder by default.
}

\examples{
\dontrun{
checkForServer()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RSelenium.R
\name{RSelenium-package}
\alias{RSelenium-package}
\alias{RSelenium}
\title{An R client for Selenium Remote Webdriver}
\description{
These are R bindings for the WebDriver API in Selenium 2.
They use the JsonWireProtocol defined at
https://github.com/SeleniumHQ/selenium/wiki/JsonWireProtocol
to communicate with a Selenium RemoteWebDriver Server.
}
\references{
http://seleniumhq.org/projects/webdriver/
}
\author{
John Harrison
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{getFirefoxProfile}
\alias{getFirefoxProfile}
\title{Get Firefox profile.}
\usage{
getFirefoxProfile(profDir, useBase = FALSE)
}
\arguments{
\item{profDir}{The directory in which the firefox profile resides}

\item{useBase}{Logical indicating whether to attempt to use zip from
utils package. Maybe easier for Windows users.}
}
\description{
\code{getFirefoxProfile}
A utility function to get a firefox profile.
}
\section{Detail}{
 A firefox profile directory is zipped and base64
   encoded. It can then be passed to the selenium server as a required
   capability with key firefox_profile
}

\examples{
\dontrun{
fprof <- getFirefoxProfile("~/.mozilla/firefox/9qlj1ofd.testprofile")
remDr <- remoteDriver(extraCapabilities = fprof)
remDr$open()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{phantom}
\alias{phantom}
\title{Start a phantomjs binary in webdriver mode.}
\usage{
phantom(pjs_cmd = "", port = 4444L, extras = "", ...)
}
\arguments{
\item{pjs_cmd}{The name, full or partial path of a phantomjs
executable. This is optional only state if the executable is not in
your path.}

\item{port}{An integer giving the port on which phantomjs will listen.
Defaults to 4444. format [[<IP>:]<PORT>]}

\item{extras}{An optional character vector: see 'Details'.}

\item{...}{Arguments to pass to \code{\link{system2}}}
}
\description{
Defunct. Please use \code{\link{rsDriver}} or \code{\link[wdman]{phantomjs}}
}
\details{
\code{phantom}
A utility function to control a phantomjs binary in webdriver mode.
}
\section{Detail}{
 phantom() is used to start a phantomjs binary in
   webdriver mode. This can be used to drive a phantomjs binary on a
   machine without selenium server. Argument extras can be used to
   specify optional extra command line arguments see
   \url{http://phantomjs.org/api/command-line.html}
}

\section{Value}{
 phantom() returns a list with two functions:
\describe{
\item{getPID}{returns the process id of the phantomjs binary running in
   webdriver mode.}
\item{stop}{terminates the phantomjs binary running in webdriver mode
   using \code{\link{pskill}}}
}
}

\examples{
\dontrun{
pJS <- phantom()
# note we are running here without a selenium server phantomjs is
# listening on port 4444
# in webdriver mode
remDr <- remoteDriver(browserName = "phantomjs")
remDr$open()
remDr$navigate("http://www.google.com/ncr")
remDr$screenshot(display = TRUE)
webElem <- remDr$findElement("name", "q")
webElem$sendKeysToElement(list("HELLO WORLD"))
remDr$screenshot(display = TRUE)
remDr$close()
# note remDr$closeServer() is not called here. We stop the phantomjs
# binary using
pJS$stop()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{getChromeProfile}
\alias{getChromeProfile}
\title{Get Chrome profile.}
\usage{
getChromeProfile(dataDir, profileDir)
}
\arguments{
\item{dataDir}{Specifies the user data directory, which is where the
browser will look for all of its state.}

\item{profileDir}{Selects directory of profile to associate with the
first browser launched.}
}
\description{
\code{getChromeProfile}
A utility function to get a Chrome profile.
}
\section{Detail}{
 A chrome profile directory is passed as an extraCapability.
The data dir has a number of default locations
\describe{
\item{Windows XP}{
Google Chrome: C:/Documents and Settings/\%USERNAME\%/Local Settings/Application Data/Google/Chrome/User Data
}
\item{Windows 8 or 7 or Vista}{
Google Chrome: C:/Users/\%USERNAME\%/AppData/Local/Google/Chrome/User Data
}
\item{Mac OS X}{
Google Chrome: ~/Library/Application Support/Google/Chrome
}
\item{Linux}{
Google Chrome: ~/.config/google-chrome
}
}
The profile directory is contained in the user directory and by default
is named "Default"
}

\examples{
\dontrun{
# example from windows using a profile directory "Profile 1"
cprof <- getChromeProfile(
  "C:\\\\Users\\\\john\\\\AppData\\\\Local\\\\Google\\\\Chrome\\\\User Data",
  "Profile 1"
)
remDr <- remoteDriver(browserName = "chrome", extraCapabilities = cprof)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{startServer}
\alias{startServer}
\title{Start the standalone server.}
\usage{
startServer(dir = NULL, args = NULL, javaargs = NULL, log = TRUE, ...)
}
\arguments{
\item{dir}{A directory in which the binary is to be placed.}

\item{args}{Additional arguments to be passed to Selenium Server.}

\item{javaargs}{arguments passed to JVM as opposed to the Selenium
Server jar.}

\item{log}{Logical value indicating whether to write a log file to the
directory containing the Selenium Server binary.}

\item{...}{arguments passed \code{\link{system2}}. Unix defaults
wait = FALSE, stdout = FALSE, stderr = FALSE. Windows defaults
wait = FALSE, invisible = TRUE.}
}
\value{
Returns a list containing two functions. The 'getpid' function
   returns the process id of the started Selenium binary. The 'stop'
   function stops the started Selenium server using the process id.
}
\description{
Defunct. Please use \code{\link{rsDriver}}
}
\details{
\code{startServer}
A utility function to start the standalone server. Return two functions
   see values.
}
\section{Detail}{
 By default the binary is assumed to be in the
   RSelenium package /bin directory. The log argument is for convenience.
   Setting it to FALSE and stipulating
   args = c("-log /user/etc/somePath/somefile.log") allows a custom
   location. Using log = TRUE sets the location to a file named
   sellog.txt in the directory containing the Selenium Server binary.
}

\examples{
\dontrun{
selServ <- startServer()
# example of commandline passing
selServ <- startServer(
  args = c("-port 4455"),
  log = FALSE,
  invisible = FALSE
)
remDr <- remoteDriver(browserName = "chrome", port = 4455)
remDr$open()
# get the process id of the selenium binary
selServ$getpid()
# stop the selenium binary
selServ$stop()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{makeFirefoxProfile}
\alias{makeFirefoxProfile}
\title{Make Firefox profile.}
\usage{
makeFirefoxProfile(opts)
}
\arguments{
\item{opts}{option list of firefox}
}
\description{
\code{makeFirefoxProfile}
A utility function to make a firefox profile.
}
\note{
Windows doesn't come with command-line zip capability.
   Installing rtools
\url{https://CRAN.R-project.org/bin/windows/Rtools/index.html} is a
   straightforward way to gain this capability.
}
\section{Detail}{
 A firefox profile directory is zipped and base64
   encoded. It can then be passed
   to the selenium server as a required capability with key
   firefox_profile
}

\examples{
\dontrun{
fprof <- makeFirefoxProfile(list(browser.download.dir = "D:/temp"))
remDr <- remoteDriver(extraCapabilities = fprof)
remDr$open()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsDriver.R
\name{rsDriver}
\alias{rsDriver}
\title{Start a selenium server and browser}
\usage{
rsDriver(
  port = 4567L,
  browser = c("chrome", "firefox", "phantomjs", "internet explorer"),
  version = "latest",
  chromever = "latest",
  geckover = "latest",
  iedrver = NULL,
  phantomver = "2.1.1",
  verbose = TRUE,
  check = TRUE,
  ...
)
}
\arguments{
\item{port}{Port to run on}

\item{browser}{Which browser to start}

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

\item{phantomver}{what version of PhantomJS to run. Default = "2.1.1"
which runs the most recent stable version. To see other version currently
sourced run binman::list_versions("phantomjs"), A value of NULL
excludes adding the PhantomJS headless browser to Selenium Server.}

\item{verbose}{If TRUE, include status messages (if any)}

\item{check}{If TRUE check the versions of selenium available and the
versions of associated drivers (chromever, geckover, phantomver,
iedrver). If new versions are available they will be downloaded.}

\item{...}{Additional arguments to pass to \code{\link{remoteDriver}}}
}
\value{
A list containing a server and a client. The server is the object
returned by \code{\link[wdman]{selenium}} and the client is an object of class
\code{\link{remoteDriver}}
}
\description{
Start a selenium server and browser
}
\details{
This function is a wrapper around \code{\link[wdman]{selenium}}.
    It provides a "shim" for the current issue running firefox on
    Windows. For a more detailed set of functions for running binaries
    relating to the Selenium/webdriver project see the
    \code{\link[wdman]{wdman}} package. Both the client and server
    are closed using a registered finalizer.
}
\examples{
\dontrun{
# start a chrome browser
rD <- rsDriver()
remDr <- rD[["client"]]
remDr$navigate("http://www.google.com/ncr")
remDr$navigate("http://www.bbc.com")
remDr$close()
# stop the selenium server
rD[["server"]]$stop()

# if user forgets to stop server it will be garbage collected.
rD <- rsDriver()
rm(rD)
gc(rD)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remoteDriver.R
\docType{class}
\name{remoteDriver-class}
\alias{remoteDriver-class}
\alias{remoteDriver}
\title{CLASS remoteDriver}
\description{
remoteDriver Class uses the JsonWireProtocol to communicate with the
   Selenium Server. If an error occurs while executing the command then
   the server sends back an HTTP error code with a JSON encoded response
   that indicates the precise Response Error Code. The remoteDriver
   class inherits from the \code{errorHandler} class. If no error
   occurred, then the subroutine called will return the value sent back
   from the server (if a return value was sent).
   So a rule of thumb while invoking methods on the driver is if the
   method did not return a status greater then zero when called, then
   you can safely assume the command was successful even if nothing was
   returned by the method.
}
\details{
remoteDriver is a generator object. To define a new remoteDriver class
   method `new` is called. The slots (default value) that are user
   defined are:
     remoteServerAddr(localhost), port(4444), browserName(firefox),
     version(""), platform(ANY),
     javascript(TRUE). See examples for more information on use.
}
\section{Fields}{

\describe{
\item{\code{remoteServerAddr}}{Object of class \code{"character"}, giving the
ip of the remote server. Defaults to localhost}

\item{\code{port}}{Object of class \code{"numeric"}, the port of the remote
server on which to connect}

\item{\code{browserName}}{Object of class \code{"character"}. The name of the
browser being used; should be one of {chrome|firefox|htmlunit|
internet explorer|iphone}.}

\item{\code{path}}{base URL path prefix for commands on the remote server.
Defaults to "/wd/hub"}

\item{\code{version}}{Object of class \code{"character"}. The browser version,
or the empty string if unknown.}

\item{\code{platform}}{Object of class \code{"character"}. A key specifying
which platform the browser is running on. This value should be one
of {WINDOWS|XP|VISTA|MAC|LINUX|UNIX}. When requesting a new session,
the client may specify ANY to indicate any available platform may be
used.}

\item{\code{javascript}}{Object of class \code{"logical"}. Whether the session
supports executing user supplied JavaScript in the context of the
current page.}

\item{\code{nativeEvents}}{Object of class \code{"logical"}. Whether the
session supports native events. n WebDriver advanced user
interactions are provided by either simulating the Javascript events
directly (i.e. synthetic events) or by letting the browser generate
the Javascript events (i.e. native events). Native events simulate
the user interactions better.}

\item{\code{serverURL}}{Object of class \code{"character"}. Url of the remote
server which JSON requests are sent to.}

\item{\code{sessionInfo}}{Object of class \code{"list"}. A list containing
information on sessions.}
}}

\section{Methods}{

\describe{
\item{\code{acceptAlert()}}{Accepts the currently displayed alert dialog.  Usually, this is
equivalent to clicking the 'OK' button in the dialog.}

\item{\code{addCookie(
  name,
  value,
  path = "/",
  domain = NULL,
  httpOnly = NULL,
  expiry = NULL,
  secure = FALSE
)}}{Set a cookie on the domain. The inputs are required apart from
those with default values.}

\item{\code{buttondown(buttonId = 0)}}{Click and hold the given mouse button (at the coordinates set by
the last moveto command). Note that the next mouse-related command
that should follow is buttondown . Any other mouse command (such
as click or another call to buttondown) will yield undefined
behaviour. buttonId - any one of 'LEFT'/0 'MIDDLE'/1 'RIGHT'/2.
Defaults to 'LEFT'}

\item{\code{buttonup(buttonId = 0)}}{Releases the mouse button previously held (where the mouse is
currently at). Must be called once for every buttondown command
issued. See the note in click and buttondown about implications of
out-of-order commands. buttonId - any one of 'LEFT'/0 'MIDDLE'/1
'RIGHT'/2. Defaults to 'LEFT'}

\item{\code{click(buttonId = 0)}}{Click any mouse button (at the coordinates set by the last
mouseMoveToLocation() command). buttonId - any one of 'LEFT'/0
'MIDDLE'/1 'RIGHT'/2. Defaults to 'LEFT'}

\item{\code{close()}}{Close the current session.}

\item{\code{closeServer()}}{Closes the server in practice terminating the process. This is
useful for linux systems. On windows the java binary operates as a
separate shell which the user can terminate.}

\item{\code{closeWindow()}}{Close the current window.}

\item{\code{deleteAllCookies()}}{Delete all cookies visible to the current page.}

\item{\code{deleteCookieNamed(name)}}{Delete the cookie with the given name. This command will be a
no-op if there is no such cookie visible to the current page.}

\item{\code{dismissAlert()}}{Dismisses the currently displayed alert dialog. For confirm() and
prompt() dialogs, this is equivalent to clicking the 'Cancel'
button. For alert() dialogs, this is equivalent to clicking the
'OK' button.}

\item{\code{doubleclick(buttonId = 0)}}{Double-Click any mouse button (at the coordinates set by the last
mouseMoveToLocation() command). buttonId - any one of 'LEFT'/0
'MIDDLE'/1 'RIGHT'/2. Defaults to 'LEFT'}

\item{\code{executeAsyncScript(script, args = list())}}{Inject a snippet of JavaScript into the page for execution in the
context of the currently selected frame. The executed script is
assumed to be asynchronous and must signal that is done by
invoking the provided callback, which is always provided as the
final argument to the function. The value to this callback will be
returned to the client. Asynchronous script commands may not span
page loads. If an unload event is fired while waiting for a script
result, an error should be returned to the client. }

\item{\code{executeScript(script, args = list(""))}}{Inject a snippet of JavaScript into the page for execution in the
context of the currently selected frame. The executed script is
assumed to be synchronous and the result of evaluating the script
is returned to the client. The script argument defines the script
to execute in the form of a function body. The value returned by
that function will be returned to the client. The function will be
invoked with the provided args array and the values may be
accessed via the arguments object in the order specified.
Arguments may be any JSON-primitive, array, or JSON object. JSON
objects that define a WebElement reference will be converted to
the corresponding DOM element. Likewise, any WebElements in the
script result will be returned to the client as WebElement JSON
objects.}

\item{\code{findElement(
  using = c("xpath", "css selector", "id", "name", "tag name", "class name",
    "link text", "partial link text"),
  value
)}}{Search for an element on the page, starting from the document
root. The located element will be returned as an object of
webElement class.The inputs are:
\describe{
  \item{\code{using}:}{Locator scheme to use to search the
    element, available schemes: Defaults to 'xpath'. Partial
    string matching is accepted.
    \describe{
      \item{"class name" :}{Returns an element whose class name
        contains the search value; compound class names are not
        permitted.}
      \item{"css selector" :}{Returns an element matching a CSS
        selector.}
      \item{"id" :}{Returns an element whose ID attribute
        matches the search value.}
      \item{"name" :}{Returns an element whose NAME attribute
        matches the search value.}
      \item{"link text" :}{Returns an anchor element whose
        visible text matches the search value.}
      \item{"partial link text" :}{Returns an anchor element
        whose visible text partially matches the search value.}
      \item{"tag name" :}{Returns an element whose tag name
        matches the search value.}
      \item{"xpath" :}{Returns an element matching an XPath
        expression.}
    }
  }
  \item{\code{value}:}{The search target. See examples.}
}}

\item{\code{findElements(
  using = c("xpath", "css selector", "id", "name", "tag name", "class name",
    "link text", "partial link text"),
  value
)}}{Search for multiple elements on the page, starting from the
document root. The located elements will be returned as an list of
objects of class WebElement. The inputs are:
\describe{
  \item{\code{using}:}{Locator scheme to use to search the
  element, available schemes: {"class name", "css selector",
  "id", "name", "link text", "partial link text",
  "tag name", "xpath" }. Defaults to 'xpath'. Partial string
  matching is accepted. See the findElement method for details}
  \item{\code{value}:}{The search target. See examples.}
}}

\item{\code{getActiveElement()}}{Get the element on the page that currently has focus. The located
element will be returned as a WebElement id.}

\item{\code{getAlertText()}}{Gets the text of the currently displayed JavaScript alert(),
confirm() or prompt() dialog.}

\item{\code{getAllCookies()}}{Retrieve all cookies visible to the current page. Each cookie
will be returned as a list with the following name and value types:
\describe{
  \item{\code{name}:}{character}
  \item{\code{value}:}{character}
  \item{\code{path}:}{character}
  \item{\code{domain}:}{character}
  \item{\code{secure}:}{logical}
}}

\item{\code{getCurrentUrl()}}{Retrieve the url of the current page.}

\item{\code{getCurrentWindowHandle()}}{Retrieve the current window handle.}

\item{\code{getLogTypes()}}{Get available log types. Common log types include 'client' = Logs
from the client, 'driver' = Logs from the webdriver, 'browser' =
Logs from the browser, 'server' = Logs from the server. Other log
types, for instance, for performance logging may also be
available. phantomjs for example returns a har log type which is a
single-entry log, with the HAR (HTTP Archive) of the current
webpage, since the first load (it's cleared at every unload event)}

\item{\code{getPageSource(...)}}{Get the current page source.}

\item{\code{getSessions()}}{Returns a list of the currently active sessions. Each session
will be returned as a list containing amongst other items:
\describe{
  \item{\code{id}:}{The session ID}
  \item{\code{capabilities}:}{An object describing session's
    capabilities}
}}

\item{\code{getStatus()}}{Query the server's current status. All server implementations
should return two basic objects describing the server's current
platform and when the server was built.}

\item{\code{getTitle(url)}}{Get the current page title.}

\item{\code{getWindowHandles()}}{Retrieve the list of window handles used in the session.}

\item{\code{getWindowPosition(windowId = "current")}}{Retrieve the window position. `windowid` is optional (default is
'current' window). Can pass an appropriate `handle`}

\item{\code{getWindowSize(windowId = "current")}}{Retrieve the window size. `windowid` is optional (default is
'current' window). Can pass an appropriate `handle`}

\item{\code{goBack()}}{Equivalent to hitting the back button on the browser.}

\item{\code{goForward()}}{Equivalent to hitting the forward button on the browser.}

\item{\code{log(type)}}{Get the log for a given log type. Log buffer is reset after each
request.
\describe{
  \item{\code{type}:}{The log type. Typically 'client', 'driver',
    'browser', 'server'}
}}

\item{\code{maxWindowSize(winHand = "current")}}{Set the size of the browser window to maximum. The windows handle
is optional. If not specified the current window in focus is used.}

\item{\code{mouseMoveToLocation(x = NA_integer_, y = NA_integer_, webElement = NULL)}}{Move the mouse by an offset of the specified element. If no
element is specified, the move is relative to the current mouse
cursor. If an element is provided but no offset, the mouse will be
moved to the center of the element. If the element is not visible,
it will be scrolled into view.}

\item{\code{navigate(url)}}{Navigate to a given url.}

\item{\code{open(silent = FALSE)}}{Send a request to the remote server to instantiate the browser.}

\item{\code{phantomExecute(script, args = list())}}{This API allows you to send a string of JavaScript via 'script',
written for PhantomJS, and be interpreted within the context of a
WebDriver Page. In other words, for the given script then this
variable is initialized to be the current Page. See
\url{https://github.com/ariya/phantomjs/wiki/API-Reference-WebPage}
and the example in this help file. NOTE: Calling the PhantomJS API
currently only works when PhantomJS is driven directly via
\code{\link{phantom}}}

\item{\code{quit()}}{Delete the session & close open browsers.}

\item{\code{refresh()}}{Reload the current page.}

\item{\code{screenshot(display = FALSE, useViewer = TRUE, file = NULL)}}{Take a screenshot of the current page. The screenshot is returned
as a base64 encoded PNG. If display is TRUE the screenshot is
displayed locally. If useViewer is TRUE and RStudio is in use the
screenshot is displayed in the RStudio viewer panel. If file is
not NULL and display = FALSE the screenshot is written to the file
denoted by file.}

\item{\code{sendKeysToActiveElement(sendKeys)}}{Send a sequence of key strokes to the active element. This
command is similar to the send keys command in every aspect except
the implicit termination: The modifiers are not released at the
end of the call. Rather, the state of the modifier keys is kept
between calls, so mouse interactions can be performed while
modifier keys are depressed. The key strokes are sent as a list.
Plain text is enter as an unnamed element of the list. Keyboard
entries are defined in `selKeys` and should be listed with name
`key`. See the examples.}

\item{\code{sendKeysToAlert(sendKeys)}}{Sends keystrokes to a JavaScript prompt() or alert() dialog.
The key strokes are sent as a list. Plain text is enter as an
unnamed element of the list. Keyboard entries are defined in
`selKeys` and should be listed with name `key`. See the examples.}

\item{\code{setAsyncScriptTimeout(milliseconds = 10000)}}{Set the amount of time, in milliseconds, that asynchronous
scripts executed by execute_async_script() are permitted to run
before they are aborted and a |Timeout| error is returned to the
client.}

\item{\code{setImplicitWaitTimeout(milliseconds = 10000)}}{Set the amount of time the driver should wait when searching for
elements. When searching for a single element, the driver will poll
the page until an element is found or the timeout expires,
whichever occurs first. When searching for multiple elements, the
driver should poll the page until at least one element is found or
the timeout expires, at which point it will return an empty list.
If this method is never called, the driver will default to an
implicit wait of 0ms.}

\item{\code{setTimeout(type = "page load", milliseconds = 10000)}}{Configure the amount of time that a particular type of operation
can execute for before they are aborted and a |Timeout| error is
returned to the client.
\describe{
  \item{\code{type}:}{The type of operation to set the timeout
    for. Valid values are: "script" for script timeouts,
    "implicit" for modifying the implicit wait timeout and
    "page load" for setting a page load timeout. Defaults to
    "page load" }
  \item{\code{milliseconds}:}{The amount of time, in
    milliseconds, that time-limited commands are permitted to run.
    Defaults to 10000 milliseconds. }
}}

\item{\code{setWindowPosition(x, y, winHand = "current")}}{Set the position (on screen) where you want your browser to be
displayed. The windows handle is optional. If not specified the
current window in focus is used.}

\item{\code{setWindowSize(width, height, winHand = "current")}}{Set the size of the browser window. The windows handle is
optional. If not specified the current window in focus is used.}

\item{\code{switchToFrame(Id)}}{Change focus to another frame on the page. Id can be
string|number|null|WebElement Object. If the Id is null, the
server should switch to the page's default content.}

\item{\code{switchToWindow(windowId)}}{Change focus to another window. The window to change focus to may
be specified by its server assigned window handle, or by the value
of its name attribute.}
}}

\examples{
\dontrun{
# start the server if one isnt running
startServer()

# use default server initialisation values
remDr <- remoteDriver$new()

# send request to server to initialise session
remDr$open()

# navigate to R home page
remDr$navigate("http://www.r-project.org")

# navigate to www.bbc.co.uk notice the need for http://
remDr$navigate("http://www.bbc.co.uk")

# go backwards and forwards
remDr$goBack()

remDr$goForward()

remDr$goBack()

# Examine the page source
frontPage <- remDr$getPageSource()

# The R homepage contains frames
webElem <- remDr$findElements(value = "//frame")
sapply(webElem, function(x) {
  x$getElementAttribute("name")
})

# The homepage contains 3 frames: logo, contents and banner
# switch to the `contents` frame
webElem <- remDr$findElement(using = "name", value = "contents")
remDr$switchToFrame(webElem$elementId)

# re-examine the page source

contentPage <- remDr$getPageSource()
identical(contentPage, frontPage) # false we hope!!

# Find the link for the search page on R homepage. Use xpath as default.
webElem <- remDr$findElement(value = '//a[@href = "search.html"]')
webElem$getElementAttribute("href")
# http://www.r-project.org/search.html

# click the search link
webElem$clickElement()

# FILL OUT A GOOGLE SEARCH FORM
remDr$navigate("http://www.google.com")

# show different methods of accessing DOM components

webElem1 <- remDr$findElement(using = "name", value = "q")
webElem2 <- remDr$findElement(
  using = "id",
  value = webElem1$getElementAttribute("id")[[1]]
)
webElem3 <- remDr$findElement(
  using = "xpath",
  value = '//input[@name = "q"]'
)

# Enter some text in the search box

webElem1$sendKeysToElement(list("RSelenium was here"))

# clear the text previously entered

webElem1$clearElement()

# show an example of sending a key press
webElem1$sendKeysToElement(list("R", key = "enter"))

# Collate the results for the `R` search
googLinkText <- remDr$findElements(value = "//h3[@class = 'r']")
linkHeading <- sapply(googLinkText, function(x) x$getElementText())
googLinkDesc <- remDr$findElements(value = "//div[@class = 's']")
linkDescription <- sapply(googLinkDesc, function(x) x$getElementText())
googLinkHref <- remDr$findElements(value = "//h3[@class = 'r']/a")
linkHref <- sapply(
  googLinkHref,
  function(x) x$getElementAttribute("href")
)

data.frame(
  heading = linkHeading,
  description = linkDescription, href = linkHref
)

# Example of javascript call
remDr$executeScript("return arguments[0] + arguments[1];", args = 1:2)
# Example of javascript async call
jsscript <-
  "arguments[arguments.length - 1](arguments[0] + arguments[1]);"
remDr$executeAsyncScript(jsscript, args = 1:2)

# EXAMPLE INJECTING INTO PHANTOMJS using phantomExecute
require(RSelenium)
pJS <- wdman::phantomjs(port = 4932L)
remDr <- remoteDriver(browserName = "phantomjs", port = 4932L)
remDr$open(silent = TRUE)
remDr$navigate("http://ariya.github.com/js/random/")
# returns a set of random numbers
remDr$findElement("id", "numbers")$getElementText()[[1]]
#  # now try injecting a new Math,random function
result <- remDr$phantomExecute("var page = this;
                               page.onInitialized = function () {
                               page.evaluate(function () {
                               Math.random = function() {return 42/100}
                               })
                               }", list())
remDr$navigate("http://ariya.github.com/js/random/")
# Math.random returns our custom function
remDr$findElement("id", "numbers")$getElementText()[[1]]
remDr$close()
pJS$stop()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errorHandler.R
\docType{class}
\name{errorHandler-class}
\alias{errorHandler-class}
\alias{errorHandler}
\title{CLASS errorHandler}
\description{
class to handle errors
}
\details{
This class is an internal class used by remoteDriver and webElement. It
   describes how drivers may respond. With a wide range of browsers etc
   the response can be variable.
}
\section{Fields}{

\describe{
\item{\code{statusCodes}}{A list with status codes and their descriptions.}

\item{\code{status}}{A status code summarizing the result of the command. A
non-zero value indicates that the command failed. A value of one is
not a failure but may  indicate a problem.}

\item{\code{statusclass}}{Class associated with the java library underlying
the server. For Example: org.openqa.selenium.remote.Response}

\item{\code{sessionid}}{An opaque handle used by the server to determine where
to route session-specific commands. This ID should be included in
all future session-commands in place of the :sessionId path segment
variable.}

\item{\code{hcode}}{A list}

\item{\code{value}}{A list containing detailed information regarding possible
errors:
\describe{
  \item{\code{message}:}{A descriptive message for the command
    failure.}
  \item{\code{screen}:}{string   (Optional) If included, a
    screenshot of the current page as a base64 encoded string.}
  \item{\code{class}:}{string   (Optional) If included, specifies
    the fully qualified class name for the exception that was thrown
    when the command failed.}
  \item{\code{stackTrace}:}{array   (Optional) If included,
    specifies an array of JSON objects describing the stack trace
    for the exception that was thrown when the command failed. The
    zeroth element of the array represents the top of the stack.}
}}

\item{\code{responseheader}}{There are two levels of error handling specified
   by the wire protocol: invalid requests and failed commands.
   Invalid Requests will probably be indicted by a status of 1.

   All invalid requests should result in the server returning a 4xx HTTP
   response. The response Content-Type should be set to text/plain and
   the message body should be a descriptive error message. The
   categories of invalid requests are as follows:
   \describe{
     \item{\code{Unknown Commands}:}{
       If the server receives a command request whose path is not mapped
       to a resource in the REST service, it should respond with a 404
       Not Found message.
     }
     \item{\code{Unimplemented Commands}:}{
       Every server implementing the WebDriver wire protocol must
       respond to every defined command. If an individual command has
       not been implemented on the server, the server should respond
       with a 501 Not Implemented error message. Note this is the only
       error in the Invalid Request category that does not return a 4xx
       status code.
     }
     \item{\code{Variable Resource Not Found}:}{
       If a request path maps to a variable resource, but that resource
       does not exist, then the server should respond with a 404 Not
       Found. For example, if ID my-session is not a valid session ID
       on the server, and a command is sent to GET /session/my-session
       HTTP/1.1, then the server should gracefully return a 404.
     }
     \item{\code{Invalid Command Method}:}{
       If a request path maps to a valid resource, but that resource
       does not respond to the request method, the server should
       respond with a 405 Method Not Allowed. The response must include
       an Allows header with a list of the allowed methods for the
       requested resource.
     }
     \item{\code{Missing Command Parameters}:}{
       If a POST/PUT command maps to a resource that expects a set of
       JSON parameters, and the response body does not include one of
       those parameters, the server should respond with a 400 Bad
       Request. The response body should list the missing parameters.
     }
   }}

\item{\code{debugheader}}{Not currently implemented}
}}

\section{Methods}{

\describe{
\item{\code{checkStatus(resContent)}}{An internal method to check the status returned by the server. If
status indicates an error an appropriate error message is thrown.}

\item{\code{errorDetails(type = "value")}}{Return error details. Type can one of c("value", "class",
"status")}

\item{\code{obscureUrlPassword(url)}}{Replaces the username and password of url with ****}

\item{\code{queryRD(ipAddr, method = "GET", qdata = NULL)}}{A method to communicate with the remote server implementing the
JSON wire protocol.}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/webElement.R
\docType{class}
\name{webElement-class}
\alias{webElement-class}
\alias{webElement}
\title{CLASS webElement}
\description{
Selenium Webdriver represents all the HTML elements as WebElements.
   This class provides a mechanism to represent them as objects &
   perform various actions on the related elements. Typically, the
   findElement method in \code{\link{remoteDriver}} returns an object
   of class webElement.
}
\details{
webElement is a generator object. To define a new webElement class
   method `new` is called.  When a webElement class is created an
   elementId should be given. Each webElement inherits from a
   remoteDriver. webElement is not usually called by the end-user.
}
\section{Fields}{

\describe{
\item{\code{elementId}}{Object of class \code{"character"}, giving a character
representation of the element id.}
}}

\section{Methods}{

\describe{
\item{\code{clearElement()}}{Clear a TEXTAREA or text INPUT element's value.}

\item{\code{clickElement()}}{Click the element.}

\item{\code{compareElements(otherElem)}}{Test if the current webElement and an other web element refer to
the same DOM element.}

\item{\code{describeElement()}}{Describe the identified element.}

\item{\code{findChildElement(
  using = c("xpath", "css selector", "id", "name", "tag name", "class name",
    "link text", "partial link text"),
  value
)}}{Search for an element on the page, starting from the node defined
by the parent webElement. The located element will be returned as
an object of webElement class.
The inputs are:
\describe{
  \item{\code{using}:}{Locator scheme to use to search the
    element, available schemes: {"class name", "css selector",
    "id", "name", "link text", "partial link text",
    "tag name", "xpath" }. Defaults to 'xpath'. Partial string
    matching is accepted.}
  \item{\code{value}:}{The search target. See examples.}
}}

\item{\code{findChildElements(
  using = c("xpath", "css selector", "id", "name", "tag name", "class name",
    "link text", "partial link text"),
  value
)}}{Search for multiple elements on the page, starting from the node
defined by the parent webElement. The located elements will be
returned as an list of objects of class WebElement.
The inputs are:
\describe{
  \item{\code{using}:}{Locator scheme to use to search the
    element, available schemes: {"class name", "css selector",
    "id", "name", "link text", "partial link text",
    "tag name", "xpath" }. Defaults to 'xpath'.
    Partial string matching is accepted.}
  \item{\code{value}:}{The search target. See examples.}
}}

\item{\code{getElementAttribute(attrName)}}{Get the value of an element's attribute. See examples.}

\item{\code{getElementLocation()}}{Determine an element's location on the page. The point (0, 0)
refers to the upper-left corner of the page.}

\item{\code{getElementLocationInView()}}{Determine an element's location on the screen once it has been
scrolled into view.
Note: This is considered an internal command and should only be
used to determine an element's location for correctly generating
native events.}

\item{\code{getElementSize()}}{Determine an element's size in pixels. The size will be returned
with width and height properties.}

\item{\code{getElementTagName()}}{Query for an element's tag name.}

\item{\code{getElementText()}}{Get the innerText of the element.}

\item{\code{getElementValueOfCssProperty(propName)}}{Query the value of an element's computed CSS property. The CSS
property to query should be specified using the CSS property name,
not the JavaScript property name (e.g. background-color instead of
backgroundColor).}

\item{\code{highlightElement(wait = 75/1000)}}{Utility function to highlight current Element. Wait denotes the
time in seconds between style changes on element.}

\item{\code{isElementDisplayed()}}{Determine if an element is currently displayed.}

\item{\code{isElementEnabled()}}{Determine if an element is currently enabled. Obviously to enable
an element just preform a click on it.}

\item{\code{isElementSelected()}}{Determine if an OPTION element, or an INPUT element of type
checkbox or radiobutton is currently selected.}

\item{\code{selectTag()}}{Utility function to return options from a select DOM node. The
option nodes are returned as webElements. The option text and the
value of the option attribute 'value' and whether the option is
selected are returned also. If this
method is called on a webElement that is not a select DOM node an
error will result.}

\item{\code{sendKeysToElement(sendKeys)}}{Send a sequence of key strokes to an element. The key strokes are
sent as a list. Plain text is enter as an unnamed element of the
list. Keyboard entries are defined in `selKeys` and should be
listed with name `key`. See the examples.}

\item{\code{setElementAttribute(attributeName, value)}}{Utility function to set an elements attributes.}

\item{\code{submitElement()}}{Submit a FORM element. The submit command may also be applied to
any element that is a descendant of a FORM element.}
}}

