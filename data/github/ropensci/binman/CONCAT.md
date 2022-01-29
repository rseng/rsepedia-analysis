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
