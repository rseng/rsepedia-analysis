# binman

<!-- badges: start -->
[![CRAN version](http://www.r-pkg.org/badges/version/binman)](https://cran.r-project.org/package=binman)
[![R build status](https://github.com/ropensci/binman/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/binman/actions)
[![codecov](https://codecov.io/gh/ropensci/binman/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/binman)
<!-- badges: end -->


Tools and functions for managing the download of binary files. Binary repositories are defined in YAML format. Defining new pre-download, download and post-download templates allow additional repositories to be added.

## Installation

You can install `binman` from GitHub with:

```R
# install.packages("remotes")
remotes::install_github("ropensci/binman")
```

## Usage Examples

### GitHub Assets

The following is an example of using `binman` to get the GitHub assets from a project. The project is https://github.com/lightbody/browsermob-proxy/releases . When a new version is released a zipped binary is added as an "asset". A JSON representation of the project releases is available at https://api.github.com/repos/lightbody/browsermob-proxy/releases. `binman` needs a YAML file to specify how to parse this projects assets:

```yaml
name: binman-bmproxy
predlfunction:
  "binman::predl_github_assets":
    url: https://api.github.com/repos/lightbody/browsermob-proxy/releases
    platform:
    - generic
    history: 3
    appname: "binman_bmproxy"
    platformregex: browsermob-proxy
dlfunction:
  "binman::download_files": []
postdlfunction:
  "binman::unziptar_dlfiles": []
```
The file can be accessed at:

```R
ymlfile <- system.file("examples", "yaml", "bmproxy.yml", package = "binman")

```

Downloading the three most recent releases can the be done using:

```R
process_yaml(ymlfile)
```

with resulting directory structure (We omit files for brevity):

#### LINUX

```sh
john@ubuntu:~$ tree -d /home/john/.local/share/binman_bmproxy
/home/john/.local/share/binman_bmproxy
└── generic
    ├── browsermob-proxy-2.1.0
    │   └── browsermob-proxy-2.1.0
    │       ├── bin
    │       │   └── conf
    │       ├── lib
    │       └── ssl-support
    ├── browsermob-proxy-2.1.1
    │   └── browsermob-proxy-2.1.1
    │       ├── bin
    │       │   └── conf
    │       ├── lib
    │       └── ssl-support
    └── browsermob-proxy-2.1.2
        └── browsermob-proxy-2.1.2
            ├── bin
            │   └── conf
            ├── lib
            └── ssl-support

19 directories
```

#### WINDOWS

```cmd
C:\Users\john>tree C:\Users\john\AppData\Local\binman\binman_bmproxy
Folder PATH listing
Volume serial number is 7CC8-BD03
C:\USERS\JOHN\APPDATA\LOCAL\BINMAN\BINMAN_BMPROXY
└───generic
    ├───browsermob-proxy-2.1.0
    │   └───browsermob-proxy-2.1.0
    │       ├───bin
    │       │   └───conf
    │       ├───lib
    │       └───ssl-support
    ├───browsermob-proxy-2.1.1
    │   └───browsermob-proxy-2.1.1
    │       ├───bin
    │       │   └───conf
    │       ├───lib
    │       └───ssl-support
    └───browsermob-proxy-2.1.2
        └───browsermob-proxy-2.1.2
            ├───bin
            │   └───conf
            ├───lib
            └───ssl-support
```

#### MACOSX

```sh
DE529:~ admin$ tree -d /Users/admin/Library/Application\ Support/binman_bmproxy
/Users/admin/Library/Application\ Support/binman_bmproxy
└── generic
    ├── browsermob-proxy-2.1.0
    │   └── browsermob-proxy-2.1.0
    │       ├── bin
    │       │   └── conf
    │       ├── lib
    │       └── ssl-support
    ├── browsermob-proxy-2.1.1
    │   └── browsermob-proxy-2.1.1
    │       ├── bin
    │       │   └── conf
    │       ├── lib
    │       └── ssl-support
    └── browsermob-proxy-2.1.2
        └── browsermob-proxy-2.1.2
            ├── bin
            │   └── conf
            ├── lib
            └── ssl-support
19 directories
```
# binman 0.1.2

* Fixed a bug in an assertion function (thanks @vjcitn #6)
* Addressed CRAN warnings and notes

# binman 0.1.1

* Fixed tests for CRAN re-submission.

# binman 0.1.0

* Add verbose argument to process_yaml.
* Moved semantic versioning to separate package

# binman 0.0.8

* Add semantic version parsing.

# binman 0.0.6

* Added a vignette on basic functionality of binman.

# binman 0.0.5

* "binman_" is now prepended to the root directory of an app to prevent 
  possible name clashes.
* Add bzip2 file format to unziptar_dlfiles

# binman 0.0.4

* Added list_versions function, docs and test to list application version
  by platform.
* Added rm_version function, docs and test to remove application versions
  by platform.
* Added rm_platform function, docs and test to remove application
  by platform.
* Added app_dir function, docs and test to return an applications root
  directory.

# binman 0.0.3

* Added Apveyor CI for windows
* Update post download zip function to handle tar and chmod

# binman 0.0.2

* Added github assets pre download template.
* Added bitbucket downloads pre download template.

# binman 0.0.1

* Added download_files, assign_directory functions.
* Added google storage pre download template.

# binman 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.



## Test environments
* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
---
title: "binman: Basics"
author: "John D Harrison"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 3
vignette: >
  %\VignetteIndexEntry{binman: Basics}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The goal of this vignette is to describe the basic functionality of the 
`binman` package.

## Introduction

`binman` (Binary Manager) is an R package that allows the user to manage
the downloading of third party binaries. The downloading of binaries is
split into three parts: the pre-download, the download and the 
post-download. Each part of the download is controlled by an appropriate 
function: **predlfunction**, **dlfunction** and **postdlfunction** 
respectively.

### Pre-Download Function

The pre-download function should return up-to its final point a named list 
of data.frames. The name of each data.frame corresponds to a "platform". 
Each individual data.frame should be composed of three columns namely
"version", "url" and "file". This named list of data.frames can then be 
passed to the `assign_directory` function. 

```
mypredl_fun <- function(...){
# code to return named list of data.frames
  assign_directory(mynamedlist, "myappname")
}
```

#### How directories are assigned

Directories are assigned using the `rappdirs` package. See `rappdirs::user_data_dir`
for more OS specific information.

### Download Function

The download function should take as its primary input the output given by 
`assign_directory` that is a named list of data.frames with column names
"version", "url", "file", "dir" and "exists". 

The download function should process the list downloading the "url" to the
given directory "dir" assigning it the file-name "file". Whether the 
download is attempted may be dictated on the value of the "exists" 
variable (For example an overwrite argument in the download function).

The download function should return a single data.frame with columns 
"platform", "file" and "processed". The "platform" column indicates what 
platform the file downloaded relates to. The "file" column gives the full
path to the downloaded file. The "processed" column indicates whether the
file was  downloaded/processed.


```
mydl_fun <- function(dllist, ...){
# dllist is output from a predlfunction with final output from assign_directory
# CODE here to perform download  
  data.frame(
    platform = platform, 
    file = file,
    processed = processed,
    stringsAsFactors = FALSE
  )
}
```
### Post-Download Function

The post-download function takes as its primary input the output given by 
a "download function" that is a data.frame with column names
"platform", "file" and "processed".

The purpose of the post-download function is to process the downloaded 
files after the download. This could involve changing the file properties 
such as the mode, unzip/untar a downloaded zip/tar file etc.

### Application YAML file

The three functions mentioned previously are combined together to define a
process for downloading the binary file. The user stipulates the required 
functions and arguments by way of a YAML format file. For example:

```
name: superduperapp
predlfunction:
  "superduperapp::predl_superduper":
    appname: "superduperapp"
dlfunction:
  "superduperapp::download_superduper": []
postdlfunction:
  "superduperapp::postdl_superduper": []
```

The defined YAML file is then processed by the `process_yaml` function. 
The YAML file needs to define the three required functions with any 
arguments and also define a name.

## Browser Mob Proxy example

The following is an example of using binman to get the github assets from a project. The project is https://github.com/lightbody/browsermob-proxy/releases . When a new version is released a zipped binary is added as an "asset". A JSON representation of the project releases is available at https://api.github.com/repos/lightbody/browsermob-proxy/releases. 
We shall breakdown the process into its three functional parts.

### BMP Pre-Download

Firstly we note we have JSON data so we parse the JSON from the URL:

```
bmpURL <- "https://api.github.com/repos/lightbody/browsermob-proxy/releases"
bmpData <- jsonlite::fromJSON(bmpURL)
```
This gives us a list of 10 releases. A github release (https://help.github.com/articles/creating-releases/) should have a version 
number. 

```
version <- bmpData[["tag_name"]]
> version
 [1] "browsermob-proxy-2.1.2"        "browsermob-proxy-2.1.1"       
 [3] "browsermob-proxy-2.1.0"        "browsermob-proxy-2.1.0-beta-6"
 [5] "browsermob-proxy-2.1.0-beta-5" "browsermob-proxy-2.1.0-beta-4"
 [7] "browsermob-proxy-2.1.0-beta-3" "browsermob-proxy-2.1.0-beta-2"
 [9] "browsermob-proxy-2.1.0-beta-1" "browsermob-proxy-2.0.0"       
```
This is our version data but we may want to tidy it up so:

```
versionregex <- c("browsermob-proxy-(.*)$", "\\1")
version <- gsub(versionregex[1], versionregex[2], version)
> version
 [1] "2.1.2"        "2.1.1"        "2.1.0"        "2.1.0-beta-6" "2.1.0-beta-5"
 [6] "2.1.0-beta-4" "2.1.0-beta-3" "2.1.0-beta-2" "2.1.0-beta-1" "2.0.0"       
```

Each release has "assets" associated with it (a zip file). 

```
assets <- bmpData[["assets"]]
> length(assets)
[1] 10
```
The assets will contain the other information we need for our pre-download function.

From each assets item we would like the file-name:

```
bmpFiles <- vapply(assets, "[[", character(1), "name")
> bmpFiles
 [1] "browsermob-proxy-2.1.2-bin.zip"        "browsermob-proxy-2.1.1-bin.zip"       
 [3] "browsermob-proxy-2.1.0-bin.zip"        "browsermob-proxy-2.1.0-beta-6-bin.zip"
 [5] "browsermob-proxy-2.1.0-beta-5-bin.zip" "browsermob-proxy-2.1.0-beta-4-bin.zip"
 [7] "browsermob-proxy-2.1.0-beta-3-bin.zip" "browsermob-proxy-2.1.0-beta-2-bin.zip"
 [9] "browsermob-proxy-2.1.0-beta-1-bin.zip" "browsermob-proxy-2.0.0-bin.zip"       
```

and the URL for that file:

```
bmpURLs <- vapply(assets, "[[", character(1), "browser_download_url")
> bmpURLs
 [1] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.2/browsermob-proxy-2.1.2-bin.zip"              
 [2] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.1/browsermob-proxy-2.1.1-bin.zip"              
 [3] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0/browsermob-proxy-2.1.0-bin.zip"              
 [4] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0-beta-6/browsermob-proxy-2.1.0-beta-6-bin.zip"
 [5] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0-beta-5/browsermob-proxy-2.1.0-beta-5-bin.zip"
 [6] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0-beta-4/browsermob-proxy-2.1.0-beta-4-bin.zip"
 [7] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0-beta-3/browsermob-proxy-2.1.0-beta-3-bin.zip"
 [8] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0-beta-2/browsermob-proxy-2.1.0-beta-2-bin.zip"
 [9] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0-beta-1/browsermob-proxy-2.1.0-beta-1-bin.zip"
[10] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.0.0/browsermob-proxy-2.0.0-bin.zip" 
```

In the case of BMP we have a single platform which we could denote "generic".
So we could pass the following:

```
dllist <- list(
  "generic" = data.frame(version = version,
                         url = bmpURLs,
                         file = bmpFiles,
                         stringsAsFactors = FALSE)
)
```

Then finally we pass this list to the `assign_directory` function:

```
predlOut <- assign_directory(dllist, "bmpApp")
> str(predlOut, max.level = 2)
List of 1
 $ generic:'data.frame':	10 obs. of  5 variables:
  ..$ version: chr [1:10] "2.1.2" "2.1.1" "2.1.0" "2.1.0-beta-6" ...
  ..$ url    : chr [1:10] "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.2/browsermob-proxy-2.1.2-bin.zip" "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.1/browsermob-proxy-2.1.1-bin.zip" "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0/browsermob-proxy-2.1.0-bin.zip" "https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0-beta-6/browsermob-proxy-2.1.0-beta-6-bin"| __truncated__ ...
  ..$ file   : chr [1:10] "browsermob-proxy-2.1.2-bin.zip" "browsermob-proxy-2.1.1-bin.zip" "browsermob-proxy-2.1.0-bin.zip" "browsermob-proxy-2.1.0-beta-6-bin.zip" ...
  ..$ dir    :List of 10
  ..$ exists : logi [1:10] FALSE FALSE FALSE FALSE FALSE FALSE ...

```

In the `binman` package there is a function template for github asset type
downloads:

```
dllist <- predl_github_assets(
  url = "https://api.github.com/repos/lightbody/browsermob-proxy/releases",
  history = 3L, 
  platform = "generic",
  appname = "bmproxy", 
  platformregex = "browsermob-proxy",
  versionregex = c("browsermob-proxy-(.*)$", "\\1")
)

```
There are some extra arguments for generality but it basically performs 
the operations we outlined above.

### BMP Download

The download for BMP is straightforward:

```
> dllist[["generic"]][1,]
  version
1   2.1.2
                                                                                                                    url
1 https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.2/browsermob-proxy-2.1.2-bin.zip
                            file                                                  dir
1 browsermob-proxy-2.1.2-bin.zip /home/john/.local/share/binman_bmproxy/generic/2.1.2
  exists
1  FALSE
```

For each platform (there is only one "generic" platform in this case) we 
want to loop over the data.frame downloading the data at "url" to the "file"
in the "dir". If it already "exists" we may not perform the download.

In `binman` there is a simple template function for this:

```
dlfiles <- download_files(dllist)
```

So we simply need to pass the list from the pre download function to this
download function.


### BMP Post Download

The files we receive from the Browser mob Proxy project are zipped. As 
a post processing operation we would like to unzip these files. The result
of unzipping the files is a directory structure. In this case we do not 
change the mode of the ultimate binary which is contained in this 
directory structure.

Again `binman` has a simple template function to unzip/untar a downloaded
file so we pass our output from the download to this:

```
dlres <- unziptar_dlfiles(dlfiles)
```

### BMP YAML

We can incorporate our three functional calls into a YAML file for the 
Browser Mob Proxy application as follows:

```
name: bmproxy
predlfunction:
  "binman::predl_github_assets":
    url: "https://api.github.com/repos/lightbody/browsermob-proxy/releases"
    platform:
    - generic
    history: 3
    appname: "bmproxy"
    platformregex:
    - "browsermob-proxy"
    versionregex:
    - "browsermob-proxy-(.*)$"
    - "\\1"
dlfunction:
  "binman::download_files": []
postdlfunction:
  "binman::unziptar_dlfiles": []
```

We have saved this file in the `binman` package:

```
ymlfile <- system.file("examples", "yaml", "bmproxy.yml", package="binman")
```

We can now use the `process-yaml` function to download the BMP binaries:

```
> process_yaml(ymlfile)
BEGIN: PREDOWNLOAD
BEGIN: DOWNLOAD
Creating directory: /home/john/.local/share/binman_bmproxy/generic/2.1.2
Downloading binary: https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.2/browsermob-proxy-2.1.2-...

Creating directory: /home/john/.local/share/binman_bmproxy/generic/2.1.1
Downloading binary: https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.1/browsermob-proxy-2.1.1-...

Creating directory: /home/john/.local/share/binman_bmproxy/generic/2.1.0
Downloading binary: https://github.com/lightbody/browsermob-proxy/releases/download/browsermob-proxy-2.1.0/browsermob-proxy-2.1.0-...

BEGIN: POSTDOWNLOAD
```

Looking at the resulting directory structure (this is for Linux the location
of directories will be OS dependent):

```
john@ubuntu:~$ tree -d /home/john/.local/share/binman_bmproxy
/home/john/.local/share/binman_bmproxy
└── generic
    ├── 2.1.0
    │   └── browsermob-proxy-2.1.0
    │       ├── bin
    │       │   └── conf
    │       ├── lib
    │       └── ssl-support
    ├── 2.1.1
    │   └── browsermob-proxy-2.1.1
    │       ├── bin
    │       │   └── conf
    │       ├── lib
    │       └── ssl-support
    └── 2.1.2
        └── browsermob-proxy-2.1.2
            ├── bin
            │   └── conf
            ├── lib
            └── ssl-support
```

We see that `binman` has downloaded and unzipped the binaries to 
appropriate directories.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predl_templates.R
\name{predl_google_storage}
\alias{predl_google_storage}
\title{Pre-Download Google Storage}
\usage{
predl_google_storage(
  url,
  platform,
  history,
  appname,
  fileregex = "\\\\.zip$",
  platformregex = platform,
  versionregex = c(paste0("(.*)/.*", fileregex), "\\\\1")
)
}
\arguments{
\item{url}{A url giving the JSON bucket listings for a project. For
example: http://chromedriver.storage.googleapis.com/index.html
lists the chromedriver files but
https://www.googleapis.com/storage/v1/b/chromedriver/o/ is the
JSON listings for the project.}

\item{platform}{A character vector of platform names}

\item{history}{The maximum number of files to get for a platform}

\item{appname}{Name of the app}

\item{fileregex}{A filter for files}

\item{platformregex}{A filter for platforms. Defaults to the platform
names.}

\item{versionregex}{A regex for retrieving the version.}
}
\value{
A named list of data.frames. The name indicates the
    platform. The data.frame should contain the version, url and file
    to be processed. Used as input for \code{\link{download_files}} or
    an equivalent.
}
\description{
Pre-Download Google Storage template function
}
\examples{
\dontrun{
gsdata <- system.file("testdata", "test_googstor.json",
  package = "binman"
)
platform <- c("linux64", "win32", "mac64")
gsdllist <- predl_google_storage(
  url = gsdata, platform, history = 5L,
  appname = "binman_chromedriver"
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binman_utils.R
\name{rm_version}
\alias{rm_version}
\title{Remove application version}
\usage{
rm_version(appname, platform, version = c("ALL"))
}
\arguments{
\item{appname}{A character string giving the name of the application}

\item{platform}{A character string indicating the platform.}

\item{version}{A character vector of versions to remove. Defaults to
"ALL"}
}
\value{
Returns a logical vector indicating whether the removal of
    version was successful. Return is invisible.
}
\description{
Remove application version for a given platform
}
\examples{
\dontrun{
appdir <- app_dir(appname, FALSE)
platforms <- LETTERS[1:4]
versions <- LETTERS[5:7]
mkdirs <- file.path(appdir, outer(platforms, versions, file.path))
chk <- vapply(mkdirs, dir.create, logical(1), recursive = TRUE)
appver <- list_versions(appname)
rm_version(appname, platforms[2], versions[1:2])
unlink(appdir, recursive = TRUE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predl_templates.R
\name{predl_bitbucket_downloads}
\alias{predl_bitbucket_downloads}
\title{Pre download bitbucket downloads}
\usage{
predl_bitbucket_downloads(
  url,
  platform,
  history,
  appname,
  platformregex = platform,
  versionregex = "\\\\d+(?:\\\\.\\\\d+)+"
)
}
\arguments{
\item{url}{A url giving the bitbucket download JSON for a project. As
an example https://bitbucket.org/ariya/phantomjs/downloads the
phantomjs project has an asset JSON available at
https://api.bitbucket.org/2.0/repositories/ariya/phantomjs/downloads?pagelen=100}

\item{platform}{A character vector of platform names}

\item{history}{The maximum number of files to get for a platform}

\item{appname}{Name of the app}

\item{platformregex}{A filter for platforms. Defaults to the platform}

\item{versionregex}{A regex for retrieving the version.}
}
\value{
A named list of data.frames. The name indicates the
    platform. The data.frame should contain the version, url and file
    to be processed. Used as input for \code{\link{download_files}} or
    an equivalent.
}
\description{
Pre download bitbucket downloads template function
}
\examples{
\dontrun{
bbdata <- system.file("testdata", "test_bitbucketdl.json",
  package = "binman"
)
platform <- c("linux64", "linux32", "windows", "macosx")
platformregex <- c("linux-x86_64", "linux-i686", "windows", "macosx")
bbdllist <-
  predl_bitbucket_downloads(
    url = bbdata, platform, history = 3L,
    appname = "binman_chromedriver",
    platformregex
  )
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binman-package.r
\docType{package}
\name{binman}
\alias{binman}
\title{binman}
\description{
A Binary Download Manager.
}
\details{
Tools and functions for managing the download of binary files.
Binary repositories are defined in 'YAML' format. Defining new pre-download,
download and post-download templates allow additional repositories to be
added.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_files.R
\name{download_files}
\alias{download_files}
\title{Download binaries}
\usage{
download_files(dllist, overwrite = FALSE)
}
\arguments{
\item{dllist}{A named list of data.frames. The data.frame should
contain the version, url and file to be processed, the directory to
download the file to and whether the file already exists.}

\item{overwrite}{Overwrite existing binaries. Default value of FALSE}
}
\value{
A data.frame indicating whether a file was
    downloaded for a platform.
}
\description{
Download binaries from repository
}
\examples{
\dontrun{
trdata <- system.file("testdata", "test_dlres.Rdata", package = "binman")
tldata <- system.file("testdata", "test_dllist.Rdata", package = "binman")
load(trdata)
load(tldata)
dllist <- assign_directory(test_dllist, "myapp")
testthat::with_mock(
  `httr::GET` = function(...) {
    test_llres
  },
  `base::dir.create` = function(...) {
    TRUE
  },
  dlfiles <- download_files(dllist)
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{reset_version}
\alias{reset_version}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{semver}{\code{\link[semver]{reset_version}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{set_version}
\alias{set_version}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{semver}{\code{\link[semver]{set_version}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binman_utils.R
\name{rm_platform}
\alias{rm_platform}
\title{Remove application platform}
\usage{
rm_platform(appname, platform = c("ALL"))
}
\arguments{
\item{appname}{A character string giving the name of the application}

\item{platform}{A character vector indicating the platform to remove.
Defaults to "ALL"}
}
\value{
Returns a logical vector indicating whether the removal of
    platform was successful. Return is invisible.
}
\description{
Remove application files/directories for a given platform
}
\examples{
\dontrun{
appdir <- app_dir(appname, FALSE)
platforms <- LETTERS[1:4]
versions <- LETTERS[5:7]
mkdirs <- file.path(appdir, outer(platforms, versions, file.path))
chk <- vapply(mkdirs, dir.create, logical(1), recursive = TRUE)
appver <- list_versions(appname)
rm_platform(appname, platforms[2:3])
unlink(appdir, recursive = TRUE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postdl_templates.R
\name{unziptar_dlfiles}
\alias{unziptar_dlfiles}
\title{Unzip/Untar downloaded files}
\usage{
unziptar_dlfiles(dlfiles, chmod = FALSE)
}
\arguments{
\item{dlfiles}{A data.frame of files by platform and indicating
whether they were processed}

\item{chmod}{change the mode of the unarchived file/files to "755" so
they are executable on unix like systems.}
}
\value{
Returns a list of character vectors indicating files
    processed
}
\description{
Unzip/Untar downloaded files. Keeps the original zip file
}
\examples{
\dontrun{
ymlfile <- system.file("exdata", "sampleapp.yml", package = "binman")
trdata <- system.file("testdata", "test_dlres.Rdata", package = "binman")
load(trdata)
testthat::with_mock(
  `httr::GET` = function(...) {
    test_llres
  },
  `base::dir.create` = function(...) {
    TRUE
  },
  `utils::unzip` = function(zipfile, ...) {
    zipfile
  },
  procyml <- process_yaml(ymlfile)
)
procyml
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{render_version}
\alias{render_version}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{semver}{\code{\link[semver]{render_version}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign_directory.R
\name{assign_directory}
\alias{assign_directory}
\title{Assign directory}
\usage{
assign_directory(dllist, appname)
}
\arguments{
\item{dllist}{A named list of data.frames. The name indicates the
platform. The data.frame should contain the version, url and file
to be processed.}

\item{appname}{Name to give the app}
}
\value{
A named list of data.frames. The data.frame should contain the
    version, url and file to be processed, the directory to download
    the file to and whether the file already exists.
}
\description{
Assign directory to download list
}
\examples{
\dontrun{
tdata <- system.file("testdata", "test_dllist.Rdata", package = "binman")
load(tdata)
assign_directory(test_dllist, "myapp")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postdl_templates.R
\name{noproc_dlfiles}
\alias{noproc_dlfiles}
\title{Do not post process}
\usage{
noproc_dlfiles(dlfiles)
}
\arguments{
\item{dlfiles}{A data.frame of files by platform and indicating
whether they were processed}
}
\value{
Returns a list of character vectors indicating files
    processed
}
\description{
Do not post process dlfiles
}
\examples{
\dontrun{
ymlfile <- system.file("exdata", "sampleapp4.yml", package = "binman")
trdata <- system.file("testdata", "test_dlres.Rdata", package = "binman")
load(trdata)
testthat::with_mock(
  `httr::GET` = function(...) {
    test_llres
  },
  `base::dir.create` = function(...) {
    TRUE
  },
  procyml <- process_yaml(ymlfile)
)
procyml
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{parse_version}
\alias{parse_version}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{semver}{\code{\link[semver]{parse_version}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binman_utils.R
\name{app_dir}
\alias{app_dir}
\title{Get application directory}
\usage{
app_dir(appname, check = TRUE)
}
\arguments{
\item{appname}{A character string giving the name of the application}

\item{check}{check whether the app given by appname exists or not.}
}
\value{
A character string giving the path of the directory
}
\description{
Get application directory
}
\examples{
\dontrun{
appdir <- app_dir("superduperapp", FALSE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binman_utils.R
\name{list_versions}
\alias{list_versions}
\title{List app versions}
\usage{
list_versions(appname, platform = c("ALL"))
}
\arguments{
\item{appname}{A character string giving the name of the application}

\item{platform}{A character vector of platforms to list. Defaults to
"ALL"}
}
\value{
A list of platforms with version directories
}
\description{
List app versions by platform
}
\examples{
\dontrun{
appdir <- app_dir("superduperapp", FALSE)
platforms <- LETTERS[1:4]
versions <- LETTERS[5:7]
mkdirs <- file.path(appdir, outer(platforms, versions, file.path))
chk <- vapply(mkdirs, dir.create, logical(1), recursive = TRUE)
expect_true(all(chk))
res <- list_versions("superduperapp")
unlink(appdir, recursive = TRUE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predl_templates.R
\name{predl_github_assets}
\alias{predl_github_assets}
\title{Pre download Github assets}
\usage{
predl_github_assets(
  url,
  platform,
  history,
  appname,
  platformregex = platform,
  versionregex = c("", "")
)
}
\arguments{
\item{url}{A url giving the github asset JSON for a project. As an
example https://github.com/mozilla/geckodriver/releases the
geckodriver project has an asset JSON available at
https://api.github.com/repos/mozilla/geckodriver/releases}

\item{platform}{A character vector of platform names}

\item{history}{The maximum number of files to get for a platform}

\item{appname}{Name of the app}

\item{platformregex}{A filter for platforms. Defaults to the platform}

\item{versionregex}{A regex for retrieving the version.}
}
\value{
A named list of data.frames. The name indicates the
    platform. The data.frame should contain the version, url and file
    to be processed. Used as input for \code{\link{download_files}} or
    an equivalent.
}
\description{
Pre download Github assets template function
}
\examples{
\dontrun{
gadata <- system.file("testdata", "test_gitassets.json",
  package = "binman"
)
platform <- c("linux64", "win64", "macos")
gadllist <- predl_github_assets(
  url = gadata, platform, history = 3L,
  appname = "binman_chromedriver"
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_yaml.R
\name{process_yaml}
\alias{process_yaml}
\title{Process a yaml file}
\usage{
process_yaml(ymlfile, verbose = TRUE)
}
\arguments{
\item{ymlfile}{A file in a YAML format defining the pre-download/
download and post download functions together with their arguments.}

\item{verbose}{If TRUE, include status messages (if any)}
}
\value{
A list of files processed (downloaded and post processed)
}
\description{
Process a yaml file. The file defines the pre-download function,
    the download function and the post download function.
}
\examples{
\dontrun{
ymlfile <- system.file("exdata", "sampleapp.yml", package = "binman")
trdata <- system.file("testdata", "test_dlres.Rdata", package = "binman")
load(trdata)
testthat::with_mock(
  `httr::GET` = function(...) {
    test_llres
  },
  `base::dir.create` = function(...) {
    TRUE
  },
  `utils::unzip` = function(zipfile, ...) {
    zipfile
  },
  procyml <- process_yaml(ymlfile)
)
procyml
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{increment_version}
\alias{increment_version}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{semver}{\code{\link[semver]{increment_version}}}
}}

