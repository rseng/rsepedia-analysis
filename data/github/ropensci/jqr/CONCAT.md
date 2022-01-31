jqr
=======




[![R-CMD-check](https://github.com/ropensci/jqr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/jqr/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/jqr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/jqr)
[![cran checks](https://cranchecks.info/badges/worst/jqr)](https://cranchecks.info/pkgs/jqr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/jqr?color=0DA6CD)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/jqr)](https://cran.r-project.org/package=jqr)

R interface to jq, a JSON processor http://stedolan.github.io/jq/

`jqr` makes it easy to process large amounts of json without having to
convert from json to R, or without using regular expressions.  This
means that the eventual loading into R can be quicker.

 - Introduction vignette at <https://cran.r-project.org/package=jqr>

## Quickstart Tutorial

The `jq` command line examples from the [jq tutorial](https://stedolan.github.io/jq/tutorial/) work exactly the same in R!


```r
library(curl)
library(jqr)
curl('https://api.github.com/repos/ropensci/jqr/commits?per_page=5') %>%
  jq('.[] | {message: .commit.message, name: .commit.committer.name}')
#> [
#>     {
#>         "message": "update cran comments and codemeta.json",
#>         "name": "Scott Chamberlain"
#>     },
#>     {
#>         "message": "change license to \"jqr authors\"",
#>         "name": "Scott Chamberlain"
#>     },
#>     {
#>         "message": "fix ci",
#>         "name": "Jeroen Ooms"
#>     },
#>     {
#>         "message": "Windows: update to libjq 1.6",
#>         "name": "Jeroen Ooms"
#>     },
#>     {
#>         "message": "Small tweak for autobrew",
#>         "name": "Jeroen Ooms"
#>     }
#> ]
```

Try running some of the [other examples](https://stedolan.github.io/jq/tutorial/).

## Installation

Binary packages for __OS-X__ or __Windows__ can be installed directly from CRAN:

```r
install.packages("jqr")
```

Installation from source on Linux or OSX requires [`libjq`](https://stedolan.github.io/jq/). On __Ubuntu 14.04 and 16.04 lower__ use [libjq-dev](https://launchpad.net/~cran/+archive/ubuntu/jq) from Launchpad:

```
sudo add-apt-repository -y ppa:cran/jq
sudo apt-get update -q
sudo apt-get install -y libjq-dev
```

More __recent Debian or Ubuntu__ install [libjq-dev](https://packages.debian.org/testing/libjq-dev) directly from Universe:

```
sudo apt-get install -y libjq-dev
```

On __Fedora__ we need [jq-devel](https://apps.fedoraproject.org/packages/jq-devel):

```
sudo yum install jq-devel
````

On __CentOS / RHEL__ we install [jq-devel](https://apps.fedoraproject.org/packages/jq-devel) via EPEL:

```
sudo yum install epel-release
sudo yum install jq-devel
```

On __OS-X__ use [jq](https://github.com/Homebrew/homebrew-core/blob/master/Formula/jq.rb) from Homebrew:

```
brew install jq
```

On __Solaris__ we can have [libjq_dev](https://www.opencsw.org/packages/libjq_dev) from [OpenCSW](https://www.opencsw.org/):
```
pkgadd -d http://get.opencsw.org/now
/opt/csw/bin/pkgutil -U
/opt/csw/bin/pkgutil -y -i libjq_dev
```


```r
library(jqr)
```

## Interfaces

### low level

There's a low level interface in which you can execute `jq` code just as you would on the command line:


```r
str <- '[{
    "foo": 1,
    "bar": 2
  },
  {
    "foo": 3,
    "bar": 4
  },
  {
    "foo": 5,
    "bar": 6
}]'
```


```r
jq(str, ".[]")
#> [
#>     {
#>         "foo": 1,
#>         "bar": 2
#>     },
#>     {
#>         "foo": 3,
#>         "bar": 4
#>     },
#>     {
#>         "foo": 5,
#>         "bar": 6
#>     }
#> ]
```


```r
jq(str, "[.[] | {name: .foo} | keys]")
#> [
#>     [
#>         "name"
#>     ],
#>     [
#>         "name"
#>     ],
#>     [
#>         "name"
#>     ]
#> ]
```

Note that we print the output to look like a valid JSON object to make it
easier to look at. However, it's a simple character string or vector of strings.
A trick you can do is to wrap your jq program in brackets like `[.[]]` instead
of `.[]`, e.g.,


```r
jq(str, ".[]") %>% unclass
#> [1] "{\"foo\":1,\"bar\":2}" "{\"foo\":3,\"bar\":4}" "{\"foo\":5,\"bar\":6}"
# vs.
jq(str, "[.[]]") %>% unclass
#> [1] "[{\"foo\":1,\"bar\":2},{\"foo\":3,\"bar\":4},{\"foo\":5,\"bar\":6}]"
```

Combine many jq arguments - they are internally combined with a pipe ` | `

(note how these are identical)


```r
jq(str, ".[] | {name: .foo} | keys")
#> [
#>     [
#>         "name"
#>     ],
#>     [
#>         "name"
#>     ],
#>     [
#>         "name"
#>     ]
#> ]
jq(str, ".[]", "{name: .foo}", "keys")
#> [
#>     [
#>         "name"
#>     ],
#>     [
#>         "name"
#>     ],
#>     [
#>         "name"
#>     ]
#> ]
```

Also accepts many JSON inputs now


```r
jq("[123, 456]   [77, 88, 99]", ".[]")
#> [
#>     123,
#>     456,
#>     77,
#>     88,
#>     99
#> ]
jq('{"foo": 77} {"bar": 45}', ".[]")
#> [
#>     77,
#>     45
#> ]
jq('[{"foo": 77, "stuff": "things"}] [{"bar": 45}] [{"n": 5}]', ".[] | keys")
#> [
#>     [
#>         "foo",
#>         "stuff"
#>     ],
#>     [
#>         "bar"
#>     ],
#>     [
#>         "n"
#>     ]
#> ]

# if you have jsons in a vector
jsons <- c('[{"foo": 77, "stuff": "things"}]', '[{"bar": 45}]', '[{"n": 5}]')
jq(paste0(jsons, collapse = " "), ".[]")
#> [
#>     {
#>         "foo": 77,
#>         "stuff": "things"
#>     },
#>     {
#>         "bar": 45
#>     },
#>     {
#>         "n": 5
#>     }
#> ]
```


### high level

The other is higher level, and uses a suite of functions to construct queries. Queries are constucted, then excuted internally with `jq()` after the last piped command.

You don't have to use pipes though. See examples below.

Examples:

Index


```r
x <- '[{"message": "hello", "name": "jenn"}, {"message": "world", "name": "beth"}]'
x %>% index()
#> [
#>     {
#>         "message": "hello",
#>         "name": "jenn"
#>     },
#>     {
#>         "message": "world",
#>         "name": "beth"
#>     }
#> ]
```

Sort


```r
'[8,3,null,6]' %>% sortj
#> [
#>     null,
#>     3,
#>     6,
#>     8
#> ]
```

reverse order


```r
'[1,2,3,4]' %>% reverse
#> [
#>     4,
#>     3,
#>     2,
#>     1
#> ]
```

Show the query to be used using `peek()`


```r
'[1,2,3,4]' %>% reverse %>% peek
#> <jq query>
#>   query: reverse
```

#### get multiple outputs for array w/ > 1 element


```r
x <- '{"user":"stedolan","titles":["JQ Primer", "More JQ"]}'
jq(x, '{user, title: .titles[]}')
#> [
#>     {
#>         "user": "stedolan",
#>         "title": "JQ Primer"
#>     },
#>     {
#>         "user": "stedolan",
#>         "title": "More JQ"
#>     }
#> ]
x %>% index()
#> [
#>     "stedolan",
#>     [
#>         "JQ Primer",
#>         "More JQ"
#>     ]
#> ]
x %>% build_object(user, title = `.titles[]`)
#> [
#>     {
#>         "user": "stedolan",
#>         "title": "JQ Primer"
#>     },
#>     {
#>         "user": "stedolan",
#>         "title": "More JQ"
#>     }
#> ]
jq(x, '{user, title: .titles[]}') %>% jsonlite::toJSON() %>% jsonlite::validate()
#> [1] TRUE
```

#### string operations

join


```r
'["a","b,c,d","e"]' %>% join
#> "a, b,c,d, e"
'["a","b,c,d","e"]' %>% join(`;`)
#> "a; b,c,d; e"
```

ltrimstr


```r
'["fo", "foo", "barfoo", "foobar", "afoo"]' %>% index() %>% ltrimstr(foo)
#> [
#>     "fo",
#>     "",
#>     "barfoo",
#>     "bar",
#>     "afoo"
#> ]
```

rtrimstr


```r
'["fo", "foo", "barfoo", "foobar", "foob"]' %>% index() %>% rtrimstr(foo)
#> [
#>     "fo",
#>     "",
#>     "bar",
#>     "foobar",
#>     "foob"
#> ]
```

startswith


```r
'["fo", "foo", "barfoo", "foobar", "barfoob"]' %>% index %>% startswith(foo)
#> [
#>     false,
#>     true,
#>     false,
#>     true,
#>     false
#> ]
'["fo", "foo"] ["barfoo", "foobar", "barfoob"]' %>% index %>% startswith(foo)
#> [
#>     false,
#>     true,
#>     false,
#>     true,
#>     false
#> ]
```

endswith


```r
'["fo", "foo", "barfoo", "foobar", "barfoob"]' %>% index %>% endswith(foo)
#> [
#>     false,
#>     true,
#>     true,
#>     false,
#>     false
#> ]
```

tojson, fromjson, tostring


```r
'[1, "foo", ["foo"]]' %>% index
#> [
#>     1,
#>     "foo",
#>     [
#>         "foo"
#>     ]
#> ]
'[1, "foo", ["foo"]]' %>% index %>% tostring
#> [
#>     "1",
#>     "foo",
#>     "[\"foo\"]"
#> ]
'[1, "foo", ["foo"]]' %>% index %>% tojson
#> [
#>     "1",
#>     "\"foo\"",
#>     "[\"foo\"]"
#> ]
'[1, "foo", ["foo"]]' %>% index %>% tojson %>% fromjson
#> [
#>     1,
#>     "foo",
#>     [
#>         "foo"
#>     ]
#> ]
```

contains


```r
'"foobar"' %>% contains("bar")
#> true
```

unique


```r
'[1,2,5,3,5,3,1,3]' %>% uniquej
#> [
#>     1,
#>     2,
#>     3,
#>     5
#> ]
```


#### filter

With filtering via `select()` you can use various operators, like `==`,
`&&`, `||`. We translate these internally for you to what `jq` wants
to see (`==`, `and`, `or`).

Simple, one condition


```r
'{"foo": 4, "bar": 7}' %>% select(.foo == 4)
#> {
#>     "foo": 4,
#>     "bar": 7
#> }
```

More complicated. Combine more than one condition; combine each individual
filtering task in parentheses


```r
x <- '{"foo": 4, "bar": 2} {"foo": 5, "bar": 4} {"foo": 8, "bar": 12}'
x %>% select((.foo < 6) && (.bar > 3))
#> {
#>     "foo": 5,
#>     "bar": 4
#> }
x %>% select((.foo < 6) || (.bar > 3))
#> [
#>     {
#>         "foo": 4,
#>         "bar": 2
#>     },
#>     {
#>         "foo": 5,
#>         "bar": 4
#>     },
#>     {
#>         "foo": 8,
#>         "bar": 12
#>     }
#> ]
```

#### types

get type information for each element


```r
'[0, false, [], {}, null, "hello"]' %>% types
#> [
#>     "number",
#>     "boolean",
#>     "array",
#>     "object",
#>     "null",
#>     "string"
#> ]
'[0, false, [], {}, null, "hello", true, [1,2,3]]' %>% types
#> [
#>     "number",
#>     "boolean",
#>     "array",
#>     "object",
#>     "null",
#>     "string",
#>     "boolean",
#>     "array"
#> ]
```

select elements by type


```r
'[0, false, [], {}, null, "hello"]' %>% index() %>% type(booleans)
#> false
```

#### key operations

get keys


```r
str <- '{"foo": 5, "bar": 7}'
str %>% keys()
#> [
#>     "bar",
#>     "foo"
#> ]
```

delete by key name


```r
str %>% del(bar)
#> {
#>     "foo": 5
#> }
```

check for key existence


```r
str3 <- '[[0,1], ["a","b","c"]]'
str3 %>% haskey(2)
#> [
#>     false,
#>     true
#> ]
str3 %>% haskey(1,2)
#> [
#>     true,
#>     false,
#>     true,
#>     true
#> ]
```

Build an object, selecting variables by name, and rename


```r
'{"foo": 5, "bar": 7}' %>% build_object(a = .foo)
#> {
#>     "a": 5
#> }
```

More complicated `build_object()`, using the included dataset `commits`


```r
commits %>%
  index() %>%
  build_object(sha = .sha, name = .commit.committer.name)
#> [
#>     {
#>         "sha": [
#>             "110e009996e1359d25b8e99e71f83b96e5870790"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "7b6a018dff623a4f13f6bcd52c7c56d9b4a4165f"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "a50e548cc5313c187483bc8fb1b95e1798e8ef65"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "4b258f7d31b34ff5d45fba431169e7fd4c995283"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "d1cb8ee0ad3ddf03a37394bfa899cfd3ddd007c5"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     }
#> ]
```

#### Maths


```r
'{"a": 7}' %>%  do(.a + 1)
#> 8
'{"a": [1,2], "b": [3,4]}' %>%  do(.a + .b)
#> [
#>     1,
#>     2,
#>     3,
#>     4
#> ]
'{"a": [1,2], "b": [3,4]}' %>%  do(.a - .b)
#> [
#>     1,
#>     2
#> ]
'{"a": 3}' %>%  do(4 - .a)
#> 1
'["xml", "yaml", "json"]' %>%  do('. - ["xml", "yaml"]')
#> ". - [\"xml\", \"yaml\"]"
'5' %>%  do(10 / . * 3)
#> 6
```

comparisons


```r
'[5,4,2,7]' %>% index() %>% do(. < 4)
#> [
#>     false,
#>     false,
#>     true,
#>     false
#> ]
'[5,4,2,7]' %>% index() %>% do(. > 4)
#> [
#>     true,
#>     false,
#>     false,
#>     true
#> ]
'[5,4,2,7]' %>% index() %>% do(. <= 4)
#> [
#>     false,
#>     true,
#>     true,
#>     false
#> ]
'[5,4,2,7]' %>% index() %>% do(. >= 4)
#> [
#>     true,
#>     true,
#>     false,
#>     true
#> ]
'[5,4,2,7]' %>% index() %>% do(. == 4)
#> [
#>     false,
#>     true,
#>     false,
#>     false
#> ]
'[5,4,2,7]' %>% index() %>% do(. != 4)
#> [
#>     true,
#>     false,
#>     true,
#>     true
#> ]
```

length


```r
'[[1,2], "string", {"a":2}, null]' %>% index %>% lengthj
#> [
#>     2,
#>     6,
#>     1,
#>     0
#> ]
```

sqrt


```r
'9' %>% sqrtj
#> 3
```

floor


```r
'3.14159' %>% floorj
#> 3
```

find minimum


```r
'[5,4,2,7]' %>% minj
#> 2
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% minj
#> {
#>     "foo": 2,
#>     "bar": 3
#> }
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% minj(foo)
#> {
#>     "foo": 1,
#>     "bar": 14
#> }
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% minj(bar)
#> {
#>     "foo": 2,
#>     "bar": 3
#> }
```

find maximum


```r
'[5,4,2,7]' %>% maxj
#> 7
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% maxj
#> {
#>     "foo": 1,
#>     "bar": 14
#> }
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% maxj(foo)
#> {
#>     "foo": 2,
#>     "bar": 3
#> }
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% maxj(bar)
#> {
#>     "foo": 1,
#>     "bar": 14
#> }
```

#### Combine into valid JSON

`jq` sometimes creates pieces of JSON that are valid in themselves, but together are not.
`combine()` is a way to make valid JSON.

This outputs a few pieces of JSON


```r
(x <- commits %>%
  index() %>%
  build_object(sha = .sha, name = .commit.committer.name))
#> [
#>     {
#>         "sha": [
#>             "110e009996e1359d25b8e99e71f83b96e5870790"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "7b6a018dff623a4f13f6bcd52c7c56d9b4a4165f"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "a50e548cc5313c187483bc8fb1b95e1798e8ef65"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "4b258f7d31b34ff5d45fba431169e7fd4c995283"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "d1cb8ee0ad3ddf03a37394bfa899cfd3ddd007c5"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     }
#> ]
```

Use `combine()` to put them together.


```r
combine(x)
#> [
#>     {
#>         "sha": [
#>             "110e009996e1359d25b8e99e71f83b96e5870790"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "7b6a018dff623a4f13f6bcd52c7c56d9b4a4165f"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "a50e548cc5313c187483bc8fb1b95e1798e8ef65"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "4b258f7d31b34ff5d45fba431169e7fd4c995283"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     },
#>     {
#>         "sha": [
#>             "d1cb8ee0ad3ddf03a37394bfa899cfd3ddd007c5"
#>         ],
#>         "name": [
#>             "Nicolas Williams"
#>         ]
#>     }
#> ]
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/jqr/issues).
* License: MIT
* Get citation information for `jqr` in R doing `citation(package = 'jqr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://www.ropensci.org/public_images/github_footer.png)](https://ropensci.org)
jqr 1.2.2
=========

* fix a failing test

jqr 1.2.1
=========

* Windows: update to libjq 1.6

jqr 1.2.0
=========

### MINOR IMPROVEMENTS

* fix to internal method `pipeline_on_exit()` for compatibility with the upcoming magrittr v2.0 (#82) thanks @lionel- !
* `jq()` errors better now when any NA class passed to it (#78)

jqr 1.1.0
=========

### NEW FEATURES

* All functions now support connection objects (file paths, urls) as input types for streaming. see new methods `jqr_feed` and `jqr_new` (#55)
* fix `jq()` to be able to accept `json` objects as input (#62)
* gains new functions `build_array()`/`build_array_()` and `build_object()`/`build_object_()` for building arrays and objects, respectively. and `select()` changes to only do filtering (instead of also doing construction) to match jq behavior (#66) (#67)
* `jq_flags()` gains parameters `stream` and `seq`

### MINOR IMPROVEMENTS

* fixed `strncpy` call in `src/jqr.c` (#75)

### BUG FIXES

* fix `select()` to handle operators other than `=` (#65) (#67)

jqr 1.0.0
=========

* Unbundle jq: the libjq library and headers are now available on all major platforms.
  See https://stedolan.github.io/jq/download/ for details. (#59)
* Removed a few authors due to n longer including jq in package
* No longer linking to BH and Rcpp. No longer using/importing Rcpp
* Use `R_registerRoutines` and `R_useDynamicSymbols` as required for 
packages with compiled code. (#57)
* Internal dataset changed name from "githubcommits" to "commits"
* Multiple JSON inputs now supported (see #53)

jqr 0.2.4
=========

### Fixes for CRAN

* Fixed the ASAN/valgrind problem (use after free) in jqr.cpp

jqr 0.2.3
=========

### Fixes for CRAN

* Fixes in v0.2.2 actually applied in this version

jqr 0.2.2
=========

### Fixes for CRAN

* Backport ec7c3cf (https://git.io/vwCFS)
* Backport eb2fc1d (https://git.io/vwCF5)
* Port 15c4a7f (https://git.io/vw1vM)

jqr 0.2.0
=========

## NEW FEATURES

* Released to CRAN.
## Test environments

* local OS X install, R 4.0.5 patched
* ubuntu 16.04 (on Github Actions), R 4.0.5
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 5 downstream dependencies. No errors were found. Summary at <https://github.com/ropensci/jqr/blob/master/revdep/README.md>.

---

This version updates the libjq version for Windows.

Thanks!
Scott Chamberlain
<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/jqr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/jqr.git`
* Make sure to track progress upstream (i.e., on our version of `jqr` at `ropensci/jqr`) by doing `git remote add upstream https://github.com/ropensci/jqr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs
* Push up to your account
* Submit a pull request to home base at `ropensci/jqr`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and the below code block and proceed :) -->

```

```
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.3 Patched (2020-11-02 r79396) |
|os       |macOS Catalina 10.15.7                      |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2020-11-12                                  |

# Dependencies

|package |old   |new   |Δ  |
|:-------|:-----|:-----|:--|
|jqr     |1.1.0 |1.2.0 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*jqr
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```


[![R-CMD-check](https://github.com/ropensci/jqr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/jqr/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/jqr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/jqr)
[![cran checks](https://cranchecks.info/badges/worst/jqr)](https://cranchecks.info/pkgs/jqr)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/jqr?color=0DA6CD)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/jqr)](https://cran.r-project.org/package=jqr)

R interface to jq, a JSON processor http://stedolan.github.io/jq/

`jqr` makes it easy to process large amounts of json without having to
convert from json to R, or without using regular expressions.  This
means that the eventual loading into R can be quicker.

 - Introduction vignette at <https://cran.r-project.org/package=jqr>

## Quickstart Tutorial

The `jq` command line examples from the [jq tutorial](https://stedolan.github.io/jq/tutorial/) work exactly the same in R! 

```{r}
library(curl)
library(jqr)
curl('https://api.github.com/repos/ropensci/jqr/commits?per_page=5') %>%
  jq('.[] | {message: .commit.message, name: .commit.committer.name}')
```

Try running some of the [other examples](https://stedolan.github.io/jq/tutorial/).

## Installation

Binary packages for __OS-X__ or __Windows__ can be installed directly from CRAN:

```r
install.packages("jqr")
```

Installation from source on Linux or OSX requires [`libjq`](https://stedolan.github.io/jq/). On __Ubuntu 14.04 and 16.04 lower__ use [libjq-dev](https://launchpad.net/~cran/+archive/ubuntu/jq) from Launchpad:

```
sudo add-apt-repository -y ppa:cran/jq
sudo apt-get update -q
sudo apt-get install -y libjq-dev
```

More __recent Debian or Ubuntu__ install [libjq-dev](https://packages.debian.org/testing/libjq-dev) directly from Universe:

```
sudo apt-get install -y libjq-dev
```

On __Fedora__ we need [jq-devel](https://apps.fedoraproject.org/packages/jq-devel):

```
sudo yum install jq-devel
````

On __CentOS / RHEL__ we install [jq-devel](https://apps.fedoraproject.org/packages/jq-devel) via EPEL:

```
sudo yum install epel-release
sudo yum install jq-devel
```

On __OS-X__ use [jq](https://github.com/Homebrew/homebrew-core/blob/master/Formula/jq.rb) from Homebrew:

```
brew install jq
```

On __Solaris__ we can have [libjq_dev](https://www.opencsw.org/packages/libjq_dev) from [OpenCSW](https://www.opencsw.org/):
```
pkgadd -d http://get.opencsw.org/now
/opt/csw/bin/pkgutil -U
/opt/csw/bin/pkgutil -y -i libjq_dev
```

```{r}
library(jqr)
```

## Interfaces

### low level

There's a low level interface in which you can execute `jq` code just as you would on the command line:

```{r}
str <- '[{
    "foo": 1,
    "bar": 2
  },
  {
    "foo": 3,
    "bar": 4
  },
  {
    "foo": 5,
    "bar": 6
}]'
```

```{r}
jq(str, ".[]")
```

```{r}
jq(str, "[.[] | {name: .foo} | keys]")
```

Note that we print the output to look like a valid JSON object to make it
easier to look at. However, it's a simple character string or vector of strings.
A trick you can do is to wrap your jq program in brackets like `[.[]]` instead
of `.[]`, e.g.,

```{r}
jq(str, ".[]") %>% unclass
# vs.
jq(str, "[.[]]") %>% unclass
```

Combine many jq arguments - they are internally combined with a pipe ` | `

(note how these are identical)

```{r}
jq(str, ".[] | {name: .foo} | keys")
jq(str, ".[]", "{name: .foo}", "keys")
```

Also accepts many JSON inputs now

```{r}
jq("[123, 456]   [77, 88, 99]", ".[]")
jq('{"foo": 77} {"bar": 45}', ".[]")
jq('[{"foo": 77, "stuff": "things"}] [{"bar": 45}] [{"n": 5}]', ".[] | keys")

# if you have jsons in a vector
jsons <- c('[{"foo": 77, "stuff": "things"}]', '[{"bar": 45}]', '[{"n": 5}]')
jq(paste0(jsons, collapse = " "), ".[]")
```


### high level

The other is higher level, and uses a suite of functions to construct queries. Queries are constucted, then excuted internally with `jq()` after the last piped command.

You don't have to use pipes though. See examples below.

Examples:

Index

```{r}
x <- '[{"message": "hello", "name": "jenn"}, {"message": "world", "name": "beth"}]'
x %>% index()
```

Sort

```{r}
'[8,3,null,6]' %>% sortj
```

reverse order

```{r}
'[1,2,3,4]' %>% reverse
```

Show the query to be used using `peek()`

```{r}
'[1,2,3,4]' %>% reverse %>% peek
```

#### get multiple outputs for array w/ > 1 element

```{r}
x <- '{"user":"stedolan","titles":["JQ Primer", "More JQ"]}'
jq(x, '{user, title: .titles[]}')
x %>% index()
x %>% build_object(user, title = `.titles[]`)
jq(x, '{user, title: .titles[]}') %>% jsonlite::toJSON() %>% jsonlite::validate()
```

#### string operations

join

```{r}
'["a","b,c,d","e"]' %>% join
'["a","b,c,d","e"]' %>% join(`;`)
```

ltrimstr

```{r}
'["fo", "foo", "barfoo", "foobar", "afoo"]' %>% index() %>% ltrimstr(foo)
```

rtrimstr

```{r}
'["fo", "foo", "barfoo", "foobar", "foob"]' %>% index() %>% rtrimstr(foo)
```

startswith

```{r}
'["fo", "foo", "barfoo", "foobar", "barfoob"]' %>% index %>% startswith(foo)
'["fo", "foo"] ["barfoo", "foobar", "barfoob"]' %>% index %>% startswith(foo)
```

endswith

```{r}
'["fo", "foo", "barfoo", "foobar", "barfoob"]' %>% index %>% endswith(foo)
```

tojson, fromjson, tostring

```{r}
'[1, "foo", ["foo"]]' %>% index
'[1, "foo", ["foo"]]' %>% index %>% tostring
'[1, "foo", ["foo"]]' %>% index %>% tojson
'[1, "foo", ["foo"]]' %>% index %>% tojson %>% fromjson
```

contains

```{r}
'"foobar"' %>% contains("bar")
```

unique

```{r}
'[1,2,5,3,5,3,1,3]' %>% uniquej
```


#### filter

With filtering via `select()` you can use various operators, like `==`, 
`&&`, `||`. We translate these internally for you to what `jq` wants 
to see (`==`, `and`, `or`).

Simple, one condition

```{r}
'{"foo": 4, "bar": 7}' %>% select(.foo == 4)
```

More complicated. Combine more than one condition; combine each individual
filtering task in parentheses

```{r}
x <- '{"foo": 4, "bar": 2} {"foo": 5, "bar": 4} {"foo": 8, "bar": 12}'
x %>% select((.foo < 6) && (.bar > 3))
x %>% select((.foo < 6) || (.bar > 3))
```

#### types

get type information for each element

```{r}
'[0, false, [], {}, null, "hello"]' %>% types
'[0, false, [], {}, null, "hello", true, [1,2,3]]' %>% types
```

select elements by type

```{r}
'[0, false, [], {}, null, "hello"]' %>% index() %>% type(booleans)
```

#### key operations

get keys

```{r}
str <- '{"foo": 5, "bar": 7}'
str %>% keys()
```

delete by key name

```{r}
str %>% del(bar)
```

check for key existence

```{r}
str3 <- '[[0,1], ["a","b","c"]]'
str3 %>% haskey(2)
str3 %>% haskey(1,2)
```

Build an object, selecting variables by name, and rename

```{r}
'{"foo": 5, "bar": 7}' %>% build_object(a = .foo)
```

More complicated `build_object()`, using the included dataset `commits`

```{r}
commits %>%
  index() %>%
  build_object(sha = .sha, name = .commit.committer.name)
```

#### Maths

```{r}
'{"a": 7}' %>%  do(.a + 1)
'{"a": [1,2], "b": [3,4]}' %>%  do(.a + .b)
'{"a": [1,2], "b": [3,4]}' %>%  do(.a - .b)
'{"a": 3}' %>%  do(4 - .a)
'["xml", "yaml", "json"]' %>%  do('. - ["xml", "yaml"]')
'5' %>%  do(10 / . * 3)
```

comparisons

```{r}
'[5,4,2,7]' %>% index() %>% do(. < 4)
'[5,4,2,7]' %>% index() %>% do(. > 4)
'[5,4,2,7]' %>% index() %>% do(. <= 4)
'[5,4,2,7]' %>% index() %>% do(. >= 4)
'[5,4,2,7]' %>% index() %>% do(. == 4)
'[5,4,2,7]' %>% index() %>% do(. != 4)
```

length

```{r}
'[[1,2], "string", {"a":2}, null]' %>% index %>% lengthj
```

sqrt

```{r}
'9' %>% sqrtj
```

floor

```{r}
'3.14159' %>% floorj
```

find minimum

```{r}
'[5,4,2,7]' %>% minj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% minj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% minj(foo)
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% minj(bar)
```

find maximum

```{r}
'[5,4,2,7]' %>% maxj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% maxj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% maxj(foo)
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% maxj(bar)
```

#### Combine into valid JSON

`jq` sometimes creates pieces of JSON that are valid in themselves, but together are not.
`combine()` is a way to make valid JSON.

This outputs a few pieces of JSON

```{r}
(x <- commits %>%
  index() %>%
  build_object(sha = .sha, name = .commit.committer.name))
```

Use `combine()` to put them together.

```{r}
combine(x)
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/jqr/issues).
* License: MIT
* Get citation information for `jqr` in R doing `citation(package = 'jqr')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![rofooter](https://www.ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "jqr introduction"
author: "Scott Chamberlain, Rich FitzJohn, Jeroen Ooms, Stefan Milton Bache"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{jqr introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
	comment = "#>",
	collapse = TRUE,
	warning = FALSE,
	message = FALSE
)
```

Introduction to jqr
===================

[jq](https://stedolan.github.io/jq/) is a _lightweight and flexible command-line JSON processor_, written in C. It's super fast, and very flexible. `jq` gives you the ability to index into, parse, and do calculations on JSON data. You can cut up and filter JSON data. You can change JSON key names and values. `jq` lets you do conditionals and comparisons, and write your own custom functions to operate on JSON data.

You can convert JSON into an R list or other R data structure, and proceed with data parsing, but why not do your JSON parsing on the actual JSON if it's easy enough?  That's where `jq` comes in. Doing your data manipulations on the actual JSON makes it easy to pass data to downstream processes that expect JSON.

If you already familiar with `jq` by using it on the command line you can use the exact same commands with `jqr`. If you've never used `jq`, `jqr` makes `jq` easy to learn with a domain specific language - and you can learn the actual `jq` syntax as you go and apply it on the command line outside of R.

## NSE vs. SE

Many functions in `jqr` have NSE (non-standard evaluation) as well as SE (standard evaluation) versions, where the NSE version for sorting an array is `sortj()` whereas the SE version is `sortj_()`. Some functions only have one version, and behave under SE rules. 

When you pass JSON into a function as the first parameter (like `ad('["a","b","c"]')`) rather than piping it in (like `'["a","b","c"]' %>% ad`), `jq()` is not executed. Rather you get back an object of class `jqr` that holds the data you passed in and the query. To execute the query on the data, run `jq()`, e.g., like `jq(ad('["a","b","c"]'))` or `ad('["a","b","c"]') %>% jq()`.

When piping JSON to DSL functions `jq()` is executed on the last DSL function used.

## jqr API

There's low and high level (or DSL [domain specific language]) interfaces in `jqr`.

### jqr low level interface

The low level and high level interfaces are unified via the function `jq()`. You can access the low leve interface by using `jq()` directly, passing a JSON string as the first parameter, the program (query) as the second, and [the flags](https://stedolan.github.io/jq/manual/#Invokingjq) as the third (by default no flags are passed).

For example, a JSON string could be `'{"a": 7, "b": 4}'`, and the program could be `.`, resulting in

```
{
    "a": 7,
    "b": 4
}
```

The program passed is exactly the same as you'd pass on the command line. Because this is a  simple replication of the command line in R, there is a higher level interface, or DSL, to make it easier to use `jq`. Nonetheless, the low level interface is important as some `jq` veterans may not want to deal with a DSL, and you may need to drop down to the low level interface if the DSL doesn't work for some reason.

### jqr DSL

The `jqr` DSL uses a suite of functions to construct queries that are executed internally with `jq()` after the last piped command. We use some logic to determine whether the function call is the last in a series of pipes, and if so, we run `jq()` on the JSON string and program/query passed.

You don't have to use pipes - they are optional. Though they do make things easier in that you can build up queries easily, just as you would with `jq`, or any other tools, on the command line.

* Execute jq
  * `jq` - execute jq
* Utility functions
  * `peek` - peek at query, without running it
  * `string` - give back character string
  * `combine` - combine pieces into proper JSON
* Identity functions
  * `dot` - takes its input and produces it unchanged as output.
  * `dotstr` - produces value at the key 'foo'
  * `index` - index to all elements, or elements by name or number
  * `indexif` - same as above, but shouldn't fail when not found
* Operations on keys, or by keys
  * `keys` - takes no input, and retrieves keys
  * `haskey` - checks if a json string has a key or keys
  * `del` - deletes provided keys
* Maths
  * `do` - arbitrary math operations
  * `lengthj` - length
  * `sqrtj` - square root
  * `floorj` - returns the floor of its numeric input
  * `minj` - minimum element of input
  * `maxj` - maximum element of input
  * `add` - adds strings or numbers together
  * `map` - for any filter X, run X for each element of input array
* Manipulation operations
  * `join` - join strings on given separator
  * `splitj` - split string on separator argument
  * `ltrimstr` - remove given prefix string, if it starts with it
  * `rtrimstr` - remove given suffix string, if it starts with it
  * `startswith` - logical output, test if input start with foo
  * `endswith` - logical output, test if input ends with foo
  * `indices` - array with numeric indices where foo occurs in inputs
  * `tojson` - dump values to JSON
  * `fromjson` - parse JSON into values
  * `tostring` - convert to string
  * `tonumber` - convert to number
  * `contains` - logical output, determine if foo is in the input
  * `uniquej` - output unique set
  * `group` - groups the elements having the same .foo field into separate arrays
* Sort
  * `sortj` - sort an array
  * `reverse` - reverse sort an array
* Types
  * `type` - select elements by type
  * `types` - get type for each element
* Functions
  * `funs` - Define and use functions
* Variables
  * `vars` - Define variables to use later
* Recursion
  * `recurse` - Search through a recursive structure - extract data from all levels
* Paths
  * `paths` - Outputs paths to all the elements in its input
* Range
  * `rangej` - Produce range of numbers
* Format strings
  * `at` - Format strings and escaping

## Load jqr

```{r}
library("jqr")
```

## Utility functions

Peek

```{r}
'{"a": 7}' %>% do(.a + 1) %>% peek
'[8,3,null,6]' %>% sortj %>% peek
```

String

```{r}
'{"a": 7}' %>% do(.a + 1) %>% string
'[8,3,null,6]' %>% sortj %>% string
```

Combine

```{r}
x <- '{"foo": 5, "bar": 7}' %>% select(a = .foo)
combine(x)
```

## index

```{r}
x <- '[{"message": "hello", "name": "jenn"}, {"message": "world", "name": "beth"}]'
x %>% index()
```

## sort

Note the function name is `sortj` to avoid collision with `base::sort`. In addition, a 
number of other functions in this package that conflict with base R functions have a 
`j` on the end.

```{r}
'[8,3,null,6]' %>% sortj
```

sort in reverse order

```{r}
'[1,2,3,4]' %>% reverse
```

## join

```{r}
'["a","b,c,d","e"]' %>% join
'["a","b,c,d","e"]' %>% join(`;`)
```

## starts- and ends-with

```{r}
'["fo", "foo", "barfoo", "foobar", "barfoob"]' %>% index %>% endswith(foo)
```

```{r}
'["fo", "foo", "barfoo", "foobar", "barfoob"]' %>% index %>% startswith(foo)
```

## contains

```{r}
'"foobar"' %>% contains("bar")
```

## unique

```{r}
'[1,2,5,3,5,3,1,3]' %>% uniquej
```

## data types

Get type information for each element

```{r}
'[0, false, [], {}, null, "hello"]' %>% types
'[0, false, [], {}, null, "hello", true, [1,2,3]]' %>% types
```

Select elements by type

```{r}
'[0, false, [], {}, null, "hello"]' %>% index() %>% type(booleans)
```

## keys

Get keys

```{r}
str <- '{"foo": 5, "bar": 7}'
str %>% keys()
```

Delete by key name

```{r}
str %>% del(bar)
str %>% del(foo)
```

Check for key existence

```{r}
str3 <- '[[0,1], ["a","b","c"]]'
str3 %>% haskey(2)
str3 %>% haskey(1,2)
```

## select 

Select variables by name, and rename

```{r}
'{"foo": 5, "bar": 7}' %>% select(a = .foo)
```

More complicated `select()`, using the included dataset `commits`

```{r}
commits %>%
  index() %>%
  build_object(sha = .sha, name = .commit.committer.name)
```

## maths 

Maths comparisons

```{r}
'[5,4,2,7]' %>% index() %>% do(. < 4)
'[5,4,2,7]' %>% index() %>% do(. > 4)
'[5,4,2,7]' %>% index() %>% do(. <= 4)
'[5,4,2,7]' %>% index() %>% do(. >= 4)
'[5,4,2,7]' %>% index() %>% do(. == 4)
'[5,4,2,7]' %>% index() %>% do(. != 4)
```

## sqrt

```{r}
'9' %>% sqrtj
```

## floor

```{r}
'3.14159' %>% floorj
```

## find minimum

```{r}
'[5,4,2,7]' %>% minj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% minj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% minj(foo)
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% minj(bar)
```

## find maximum

```{r}
'[5,4,2,7]' %>% maxj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% maxj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% maxj(foo)
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' %>% maxj(bar)
```

## connections

### files

```{r}
tmp <- tempfile()
writeLines(c("[123, 456]", "[77, 88, 99]", "[41]"), tmp)
jq(file(tmp), ".[]")
```

### urls

```{r eval=FALSE}
x <- 'http://jeroen.github.io/data/diamonds.json'
jq(url(x), "select(.carat > 3.5)")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paths.R
\name{paths}
\alias{paths}
\title{Outputs paths to all the elements in its input}
\usage{
paths(.data)
}
\arguments{
\item{.data}{input}
}
\description{
Outputs paths to all the elements in its input
}
\examples{
'[1,[[],{"a":2}]]' \%>\% paths
'[{"name":"JSON", "good":true}, {"name":"XML", "good":false}]' \%>\% paths
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logicaltests.R
\name{logicaltests}
\alias{logicaltests}
\alias{allj}
\alias{anyj}
\title{Logical tests}
\usage{
allj(.data)

anyj(.data)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}
}
\description{
Logical tests
}
\examples{
# any
'[true, false]' \%>\% anyj
'[false, false]' \%>\% anyj
'[]' \%>\% anyj

# all
'[true, false]' \%>\% allj
'[true, true]' \%>\% allj
'[]' \%>\% allj

## many JSON inputs
'[true, false] [true, true] [false, false]' \%>\% anyj
'[true, false] [true, true] [false, false]' \%>\% allj
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/types.R
\name{types}
\alias{types}
\alias{type}
\alias{type_}
\title{Types and related functions}
\usage{
types(.data)

type(.data, ...)

type_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
Types and related functions
}
\examples{
# get type information for each element
jq('[0, false, [], {}, null, "hello"]', 'map(type)')
'[0, false, [], {}, null, "hello"]' \%>\% types
'[0, false, [], {}, null, "hello", true, [1,2,3]]' \%>\% types

# select elements by type
jq('[0, false, [], {}, null, "hello"]', '.[] | numbers,booleans')
'[0, false, [], {}, null, "hello"]' \%>\% index() \%>\% type(booleans)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jqr-package.R
\docType{package}
\name{jqr}
\alias{jqr}
\alias{jqr-package}
\title{jqr}
\description{
An R client for the C library jq
}
\section{Low-level}{

Low level interface, in which you can execute `jq` code just as you
would on the command line. Available via \code{\link{jq}}
}

\section{High-level DSL}{

High-level, uses a suite of functions to construct queries. Queries
are constucted, then excuted internally with \code{\link{jq}}
}

\section{Pipes}{

The high level DSL supports piping, though you don't have to use
pipes.
}

\section{NSE and SE}{

Most DSL functions have NSE (non-standard evaluation) and SE
(standard evaluation) versions, which make \code{jqr} easy to use
for interactive use as well as programming.
}

\section{jq version}{

We link to \code{jq} through the installed version on your system,
so the version can vary. Run \code{jq --version} to get your jq version
}

\section{indexing}{

note that \code{jq} indexing starts at \code{0}, whereas R indexing
starts at \code{1}. So when you want the first thing in an array using
\code{jq}, for example, you want \code{0}, not \code{1}
}

\section{output data format}{

Note that with both the low level interface and the high level DSL, we
print the output to look like a valid JSON object to make it easier to
look at. However, it's important to know that the output is really just a
simple character string or vector of strings - it's just the print function
that pretty prints it and makes it look like a single JSON object. What jq
is giving you often is a stream of valid JSON objects, each one of which is
valid, but altogether are not valid. However, a trick you can do is to
wrap your jq program in brackets like \code{[.[]]} instead of \code{.[]}
to give a single JSON object

Related to above, you can use the function provided \code{\link{string}}
with the high level DSL to get back a character string instead of
pretty printed version
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vars.R
\name{vars}
\alias{vars}
\alias{vars_}
\title{Variables}
\usage{
vars(.data, ...)

vars_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
Variables
}
\examples{
x <- '{
 "posts": [
   {"title": "Frist psot", "author": "anon"},
   {"title": "A well-written article", "author": "person1"}
 ],
 "realnames": {
   "anon": "Anonymous Coward",
   "person1": "Person McPherson"
 }
}'

x \%>\% dotstr(posts[])
x \%>\% dotstr(posts[]) \%>\% string
x \%>\% vars(realnames = names) \%>\% dotstr(posts[]) \%>\%
   build_object(title, author = "$names[.author]")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jq.R
\name{jq}
\alias{jq}
\alias{jq.jqr}
\alias{jq.character}
\alias{jq.json}
\alias{jq.connection}
\title{Execute a query with jq}
\usage{
jq(x, ...)

\method{jq}{jqr}(x, ...)

\method{jq}{character}(x, ..., flags = jq_flags())

\method{jq}{json}(x, ..., flags = jq_flags())

\method{jq}{connection}(x, ..., flags = jq_flags(), out = NULL)
}
\arguments{
\item{x}{\code{json} object or character string with json data. this can
be one or more valid json objects}

\item{...}{character specification of jq query. Each element in code{...}
will be combined with " | ", which is convenient for long queries.}

\item{flags}{See \code{\link{jq_flags}}}

\item{out}{a filename, callback function, connection object to stream output.
Set to `NULL` to buffer all output and return a character vector.}
}
\description{
\code{jq} is meant to work with the high level interface in this package.
\code{jq} also provides access to the low level interface in which you can
use jq query strings just as you would on the command line. Output gets
class of json, and pretty prints to the console for easier viewing.
\code{jqr} doesn't do pretty printing.
}
\examples{
'{"a": 7}' \%>\%  do(.a + 1)
'[8,3,null,6]' \%>\% sortj

x <- '[{"message": "hello", "name": "jenn"},
  {"message": "world", "name": "beth"}]'
jq(index(x))

jq('{"a": 7, "b": 4}', 'keys')
jq('[8,3,null,6]', 'sort')

# many json inputs
jq(c("[123, 456]", "[77, 88, 99]", "[41]"), ".[]")
# Stream from connection
tmp <- tempfile()
writeLines(c("[123, 456]", "[77, 88, 99]", "[41]"), tmp)
jq(file(tmp), ".[]")

\dontrun{
# from a url
x <- 'http://jeroen.github.io/data/diamonds.json'
jq(url(x), ".[]")

# from a file
file <- file.path(tempdir(), "diamonds_nd.json")
download.file(x, destfile = file)
jq(file(file), ".carat")
jq(file(file), "select(.carat > 1.5)")
jq(file(file), 'select(.carat > 4 and .cut == "Fair")')
}
}
\seealso{
\code{\link{peek}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/at.R
\name{at}
\alias{at}
\alias{at_}
\title{Format strings and escaping}
\usage{
at(.data, ...)

at_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
Format strings and escaping
}
\examples{
x <- '{"user":"stedolan","titles":["JQ Primer", "More JQ"]}'
x \%>\% at(base64) \%>\% peek
x \%>\% at(base64)
x \%>\% index() \%>\% at(base64)

y <- '["fo", "foo", "barfoo", "foobar", "barfoob"]'
y \%>\% index() \%>\% at(base64)

## prepare for shell use
y \%>\% index() \%>\% at(sh)

## rendered as csv with double quotes
z <- '[1, 2, 3, "a"]'
z \%>\% at(csv)

## rendered as csv with double quotes
z \%>\% index()
z \%>\% index() \%>\% at(text)

## \% encode for URI's
#### DOESNT WORK --------------------------

## html escape
#### DOESNT WORK --------------------------

## serialize to json
#### DOESNT WORK --------------------------
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/range.R
\name{rangej}
\alias{rangej}
\title{Produce range of numbers}
\usage{
rangej(x, array = FALSE)
}
\arguments{
\item{x}{Input, single number or number range.}

\item{array}{(logical) Create array. Default: \code{FALSE}}
}
\description{
Produce range of numbers
}
\examples{
2:4 \%>\% rangej
2:1000 \%>\% rangej
1 \%>\% rangej
4 \%>\% rangej
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jqr-package.R
\docType{data}
\name{commits}
\alias{commits}
\title{GitHub Commits Data}
\format{
A character string of json github commits data for the jq repo.
}
\description{
GitHub Commits Data
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index.R
\name{index}
\alias{index}
\alias{index_}
\alias{indexif}
\alias{indexif_}
\alias{dotindex}
\alias{dotindex_}
\title{index and related functions}
\usage{
index(.data, ...)

index_(.data, ..., .dots)

indexif(.data, ...)

indexif_(.data, ..., .dots)

dotindex(.data, ...)

dotindex_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
index and related functions
}
\details{
\itemize{
 \item \code{index}/\code{index_} - queries like: \code{.[]}, \code{.[0]},
 \code{.[1:5]},
 \code{.["foo"]}
 \item \code{indexif}/\code{indexif_} - queries like: \code{.["foo"]?}
 \item \code{dotindex}/\code{dotindex_} - queries like: \code{.[].foo},
 \code{.[].foo.bar}
}
}
\examples{
str <- '[{"name":"JSON", "good":true}, {"name":"XML", "good":false}]'
str \%>\% index
'{"name":"JSON", "good":true}' \%>\% indexif(name)
'{"name":"JSON", "good":true}' \%>\% indexif(good)
'{"name":"JSON", "good":true}' \%>\% indexif(that)
'{"a": 1, "b": 1}' \%>\% index
'[]' \%>\% index
'[{"name":"JSON", "good":true}, {"name":"XML", "good":false}]' \%>\% index(0)
'["a","b","c","d","e"]' \%>\% index(2)
'["a","b","c","d","e"]' \%>\% index('2:4')
'["a","b","c","d","e"]' \%>\% index('2:5')
'["a","b","c","d","e"]' \%>\% index(':3')
'["a","b","c","d","e"]' \%>\% index('-2:')

str \%>\% index \%>\% select(bad = .name)

'[{"name":"JSON", "good":true}, {"name":"XML", "good":false}]' \%>\%
  dotindex(name)
'[{"name":"JSON", "good":true}, {"name":"XML", "good":false}]' \%>\%
  dotindex(good)
'[{"name":"JSON", "good":{"foo":5}}, {"name":"XML", "good":{"foo":6}}]' \%>\%
  dotindex(good)
'[{"name":"JSON", "good":{"foo":5}}, {"name":"XML", "good":{"foo":6}}]' \%>\%
  dotindex(good.foo)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/keys.R
\name{keys}
\alias{keys}
\alias{del}
\alias{del_}
\alias{haskey}
\alias{haskey_}
\title{Operations on keys, or by keys}
\usage{
keys(.data)

del(.data, ...)

del_(.data, ..., .dots)

haskey(.data, ...)

haskey_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
\code{keys} takes no input, and retrieves keys. \code{del} deletes
provided keys. \code{haskey} checks if a json string has a key, or the
input array has an element at the given index.
}
\examples{
# get keys
str <- '{"foo": 5, "bar": 7}'
jq(str, "keys")
str \%>\% keys()

# delete by key name
jq(str, "del(.bar)")
str \%>\% del(bar)

# check for key existence
str3 <- '[[0,1], ["a","b","c"]]'
jq(str3, "map(has(2))")
str3 \%>\% haskey(2)
jq(str3, "map(has(1,2))")
str3 \%>\% haskey(1,2)

## many JSON inputs
'{"foo": 5, "bar": 7} {"hello": 5, "world": 7}' \%>\% keys
'{"foo": 5, "bar": 7} {"hello": 5, "bar": 7}' \%>\% del(bar)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flags.R
\name{jq_flags}
\alias{jq_flags}
\alias{flags}
\title{Flags for use with jq}
\usage{
jq_flags(
  pretty = FALSE,
  ascii = FALSE,
  color = FALSE,
  sorted = FALSE,
  stream = FALSE,
  seq = FALSE
)

flags(
  .data,
  pretty = FALSE,
  ascii = FALSE,
  color = FALSE,
  sorted = FALSE,
  stream = FALSE,
  seq = FALSE
)
}
\arguments{
\item{pretty}{Pretty print the json (different to jsonlite's
pretty printing).}

\item{ascii}{Force jq to produce pure ASCII output with non-ASCII
characters replaced by equivalent escape sequences.}

\item{color}{Add ANSI escape sequences for coloured output}

\item{sorted}{Output fields of each object with keys in sorted order}

\item{stream}{Parse the input in streaming fashion, outputing arrays
of path and leaf values like \code{jq --stream} command line.}

\item{seq}{Use the application/json-seq MIME type scheme for separating
JSON like the \code{jq --seq} command line.}

\item{.data}{A \code{jqr} object.}
}
\description{
The \code{flags} function is provided for the high-level DSL
approach, whereas the \code{jq_flags} function is used to provide
the low-level \code{jq} with the appropriate flags.
}
\examples{
'{"a": 7, "z":0, "b": 4}' \%>\% flags(sorted = TRUE)
'{"a": 7, "z":0, "b": 4}' \%>\% dot \%>\% flags(sorted = TRUE)
jq('{"a": 7, "z":0, "b": 4}', ".") \%>\% flags(sorted = TRUE)
jq('{"a": 7, "z":0, "b": 4}', ".", flags = jq_flags(sorted = TRUE))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/string.R
\name{string}
\alias{string}
\title{Give back a character string}
\usage{
string(.data)
}
\arguments{
\item{.data}{(list) input, using higher level interface}
}
\description{
Give back a character string
}
\examples{
'{"a": 7}' \%>\% do(.a + 1) \%>\% string
'[8,3,null,6]' \%>\% sortj \%>\% string
}
\seealso{
\code{\link{peek}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dot.R
\name{dot}
\alias{dot}
\alias{dot_}
\alias{dotstr}
\alias{dotstr_}
\title{dot and related functions}
\usage{
dot(.data)

dot_(.data, dots = ".")

dotstr(.data, ...)

dotstr_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{dots}{dots}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}
}
\description{
dot and related functions
}
\examples{
str <- '[{"name":"JSON", "good":true}, {"name":"XML", "good":false}]'
str \%>\% dot
str \%>\% index \%>\% dotstr(name)
'{"foo": 5, "bar": 8}' \%>\% dot
'{"foo": 5, "bar": 8}' \%>\% dotstr(foo)
'{"foo": {"bar": 8}}' \%>\% dotstr(foo.bar)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine.R
\name{combine}
\alias{combine}
\title{Combine json pieces}
\usage{
combine(x)
}
\arguments{
\item{x}{Input, of class json}
}
\description{
Combine json pieces
}
\examples{
x <- '{"foo": 5, "bar": 7}' \%>\% select(a = .foo)
combine(x)

(x <- commits \%>\% index() \%>\%
 select(sha = .sha, name = .commit.committer.name))
combine(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peek.R
\name{peek}
\alias{peek}
\title{Peek at a query}
\usage{
peek(.data)
}
\arguments{
\item{.data}{(list) input, using higher level interface}
}
\description{
Prints the query resulting from \code{jq} all in one character string just
as you would execute it on the command line. Output gets class of json,
and pretty prints to the console for easier viewing.
}
\examples{
'{"a": 7}' \%>\% do(.a + 1) \%>\% peek
'[8,3,null,6]' \%>\% sortj \%>\% peek
}
\seealso{
\code{\link{jq}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
Pipe operator
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jqr.R
\name{jqr_new}
\alias{jqr_new}
\alias{jqr_feed}
\title{JQ Streaming API}
\usage{
jqr_new(query, flags = jq_flags())

jqr_feed(jqr_program, json, unlist = TRUE, finalize = FALSE)
}
\arguments{
\item{query}{string with a valid jq program}

\item{flags}{See \code{\link{jq_flags}}}

\item{jqr_program}{object returned by [jqr_new]}

\item{json}{character vector with json data. If the JSON object is incomplete, you
must set `finalize` to `FALSE` otherwise you get an error.}

\item{unlist}{if `TRUE` returns a single character vector with all output for each
each string in `json` input}

\item{finalize}{completes the parsing and verifies that the JSON string is valid. Set
this to `TRUE` when feeding the final piece of data.}
}
\description{
Low level JQ API. First create a program using a `query` and `flags` and then
feed pieces of data.
}
\examples{
program <- jqr_new(".[]")
jqr_feed(program, c("[123, 456]", "[77, 88, 99]"))
jqr_feed(program, c("[41, 234]"))
jqr_feed(program, "", finalize = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manip.R
\name{manip}
\alias{manip}
\alias{join}
\alias{join_}
\alias{splitj}
\alias{splitj_}
\alias{ltrimstr}
\alias{ltrimstr_}
\alias{rtrimstr}
\alias{rtrimstr_}
\alias{startswith}
\alias{startswith_}
\alias{endswith}
\alias{endswith_}
\alias{index_loc}
\alias{index_loc_}
\alias{rindex_loc}
\alias{rindex_loc_}
\alias{indices}
\alias{indices_}
\alias{tojson}
\alias{fromjson}
\alias{tostring}
\alias{tonumber}
\alias{contains}
\alias{contains_}
\alias{uniquej}
\alias{uniquej_}
\alias{group}
\alias{group_}
\title{Manipulation operations}
\usage{
join(.data, ...)

join_(.data, ..., .dots)

splitj(.data, ...)

splitj_(.data, ..., .dots)

ltrimstr(.data, ...)

ltrimstr_(.data, ..., .dots)

rtrimstr(.data, ...)

rtrimstr_(.data, ..., .dots)

startswith(.data, ...)

startswith_(.data, ..., .dots)

endswith(.data, ...)

endswith_(.data, ..., .dots)

index_loc(.data, ...)

index_loc_(.data, ..., .dots)

rindex_loc(.data, ...)

rindex_loc_(.data, ..., .dots)

indices(.data, ...)

indices_(.data, ..., .dots)

tojson(.data)

fromjson(.data)

tostring(.data)

tonumber(.data)

contains(.data, ...)

contains_(.data, ..., .dots)

uniquej(.data, ...)

uniquej_(.data, ..., .dots)

group(.data, ...)

group_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
Manipulation operations
}
\examples{
# join
str <- '["a","b,c,d","e"]'
jq(str, 'join(", ")')
str \%>\% join
str \%>\% join(`;`)
str \%>\% join(`yep`)
## many JSON inputs
'["a","b,c,d","e"] ["a","f,e,f"]' \%>\% join(`---`)

# split
jq('"a, b,c,d, e"', 'split(", ")')

# ltrimstr
jq('["fo", "foo", "barfoo", "foobar", "afoo"]', '[.[]|ltrimstr("foo")]')
'["fo", "foo", "barfoo", "foobar", "afoo"]' \%>\% index() \%>\% ltrimstr(foo)

# rtrimstr
jq('["fo", "foo", "barfoo", "foobar", "foob"]', '[.[]|rtrimstr("foo")]')
'["fo", "foo", "barfoo", "foobar", "foob"]' \%>\% index() \%>\% rtrimstr(foo)

# startswith
str <- '["fo", "foo", "barfoo", "foobar", "barfoob"]'
jq(str, '[.[]|startswith("foo")]')
str \%>\% index \%>\% startswith(foo)
## many JSON inputs
'["fo", "foo"] ["barfoo", "foobar", "barfoob"]' \%>\% index \%>\% startswith(foo)

# endswith
jq(str, '[.[]|endswith("foo")]')
str \%>\% index \%>\% endswith(foo)
str \%>\% index \%>\% endswith_("foo")
str \%>\% index \%>\% endswith(bar)
str \%>\% index \%>\% endswith_("bar")
## many JSON inputs
'["fo", "foo"] ["barfoo", "foobar", "barfoob"]' \%>\% index \%>\% endswith(foo)

# get index (location) of a character
## input has to be quoted
str <- '"a,b, cd, efg, hijk"'
str \%>\% index_loc(", ")
str \%>\% index_loc(",")
str \%>\% index_loc("j")
str \%>\% rindex_loc(", ")
str \%>\% indices(", ")

# tojson, fromjson, tostring, tonumber
'[1, "foo", ["foo"]]' \%>\% index \%>\% tostring
'[1, "1"]' \%>\% index \%>\% tonumber
'[1, "foo", ["foo"]]' \%>\% index \%>\% tojson
'[1, "foo", ["foo"]]' \%>\% index \%>\% tojson \%>\% fromjson

# contains
'"foobar"' \%>\% contains("bar")
'["foobar", "foobaz", "blarp"]' \%>\% contains(`["baz", "bar"]`)
'["foobar", "foobaz", "blarp"]' \%>\% contains(`["bazzzzz", "bar"]`)
str <- '{"foo": 12, "bar":[1,2,{"barp":12, "blip":13}]}'
str \%>\% contains(`{foo: 12, bar: [{barp: 12}]}`)
str \%>\% contains(`{foo: 12, bar: [{barp: 15}]}`)

# unique
'[1,2,5,3,5,3,1,3]' \%>\% uniquej
str <- '[{"foo": 1, "bar": 2}, {"foo": 1, "bar": 3}, {"foo": 4, "bar": 5}]'
str \%>\% uniquej(foo)
str \%>\% uniquej_("foo")
'["chunky", "bacon", "kitten", "cicada", "asparagus"]' \%>\% uniquej(length)

# group
x <- '[{"foo":1, "bar":10}, {"foo":3, "bar":100}, {"foo":1, "bar":1}]'
x \%>\% group(foo)
x \%>\% group_("foo")
}
\seealso{
\code{\link{add}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.R
\name{select}
\alias{select}
\alias{select_}
\title{Select - filtering}
\usage{
select(.data, ...)

select_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
The function \code{select(foo)} produces its input unchanged if 
\code{foo} returns TRUE for that input, and produces no output otherwise
}
\note{
this function has changed what it does dramatically. we were 
using this function for object construction, which is now done with 
\code{\link{build_object}}
}
\examples{
jq('[1,5,3,0,7]', 'map(select(. >= 2))')
'[1,5,3,0,7]' \%>\% map(select(. >= 2)) 


'{"foo": 4, "bar": 7}' \%>\% select(.foo == 4)
'{"foo": 5, "bar": 7} {"foo": 4, "bar": 7}' \%>\% select(.foo == 4)
'[{"foo": 5, "bar": 7}, {"foo": 4, "bar": 7}]' \%>\% index() \%>\% 
  select(.foo == 4)
'{"foo": 4, "bar": 7} {"foo": 5, "bar": 7} {"foo": 8, "bar": 7}' \%>\% 
  select(.foo < 6)

x <- '{"foo": 4, "bar": 2} {"foo": 5, "bar": 4} {"foo": 8, "bar": 12}'
jq(x, 'select((.foo < 6) and (.bar > 3))')
jq(x, 'select((.foo < 6) or (.bar > 3))')
x \%>\% select((.foo < 6) && (.bar > 3))
x \%>\% select((.foo < 6) || (.bar > 3))

x <- '[{"foo": 5, "bar": 7}, {"foo": 4, "bar": 7}, {"foo": 4, "bar": 9}]'
jq(x, '.[] | select(.foo == 4) | {user: .bar}')
x \%>\% index() \%>\% select(.foo == 4) \%>\% build_object(user = .bar)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructors.R
\name{build}
\alias{build}
\alias{build_array}
\alias{build_array_}
\alias{build_object}
\alias{build_object_}
\title{Build arrays and objects}
\usage{
build_array(.data, ...)

build_array_(.data, ..., .dots)

build_object(.data, ...)

build_object_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
Build arrays and objects
}
\examples{
## BUILD ARRAYS
x <- '{"user":"stedolan", "projects": ["jq", "wikiflow"]}' 
jq(x, "[.user, .projects[]]")
x \%>\% build_array(.user, .projects[])

jq('[1, 2, 3]', '[ .[] | . * 2]')
'[1, 2, 3]' \%>\% build_array(.[] | . * 2)


## BUILD OBJECTS
'{"foo": 5, "bar": 7}' \%>\% build_object(a = .foo) \%>\% peek
'{"foo": 5, "bar": 7}' \%>\% build_object(a = .foo)

# using json dataset, just first element
x <- commits \%>\% index(0)
x \%>\%
   build_object(message = .commit.message, name = .commit.committer.name)
x \%>\% build_object(sha = .commit.tree.sha, author = .author.login)

# using json dataset, all elements
x <- index(commits)
x \%>\% build_object(message = .commit.message, name = .commit.committer.name)
x \%>\% build_object(sha = .sha, name = .commit.committer.name)

# many JSON inputs
'{"foo": 5, "bar": 7} {"foo": 50, "bar": 7} {"foo": 500, "bar": 7}' \%>\%
  build_object(hello = .foo)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recurse.R
\name{recurse}
\alias{recurse}
\alias{recurse_}
\title{Search through a recursive structure - extract data from all levels}
\usage{
recurse(.data, ...)

recurse_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
Search through a recursive structure - extract data from all levels
}
\examples{
x <- '{"name": "/", "children": [
  {"name": "/bin", "children": [
    {"name": "/bin/ls", "children": []},
    {"name": "/bin/sh", "children": []}]},
  {"name": "/home", "children": [
    {"name": "/home/stephen", "children": [
      {"name": "/home/stephen/jq", "children": []}]}]}]}'
x \%>\% recurse(.children[]) \%>\% build_object(name)
x \%>\% recurse(.children[]) \%>\% build_object(name) \%>\% string
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs.R
\name{funs}
\alias{funs}
\title{Define and use functions}
\usage{
funs(.data, fxn, action)
}
\arguments{
\item{.data}{input}

\item{fxn}{A function definition, without \code{def} (added internally)}

\item{action}{What to do with the function on the data}
}
\description{
Define and use functions
}
\examples{
jq("[1,2,10,20]", 'def increment: . + 1; map(increment)')
"[1,2,10,20]" \%>\% funs('increment: . + 1', 'map(increment)')
"[1,2,10,20]" \%>\% funs('increment: . / 100', 'map(increment)')
"[1,2,10,20]" \%>\% funs('increment: . / 100', 'map(increment)')
'[[1,2],[10,20]]' \%>\% funs('addvalue(f): f as $x | map(. + $x)', 'addvalue(.[0])')
"[1,2]" \%>\% funs('f(a;b;c;d;e;f): [a+1,b,c,d,e,f]', 'f(.[0];.[1];.[0];.[0];.[0];.[0])')
"[1,2,3,4]" \%>\% funs('fac: if . == 1 then 1 else . * (. - 1 | fac) end', '[.[] | fac]')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sort.R
\name{sortj}
\alias{sortj}
\alias{sortj_}
\alias{reverse}
\title{Sort and related}
\usage{
sortj(.data, ...)

sortj_(.data, ..., .dots)

reverse(.data)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
Sort and related
}
\examples{
# sort
'[8,3,null,6]' \%>\% sortj
'[{"foo":4, "bar":10}, {"foo":3, "bar":100}, {"foo":2, "bar":1}]' \%>\%
  sortj(foo)

# reverse order
'[1,2,3,4]' \%>\% reverse

# many JSON inputs
'[{"foo":7}, {"foo":4}] [{"foo":300}, {"foo":1}] [{"foo":2}, {"foo":1}]' \%>\%
  sortj(foo)

'[1,2,3,4] [10,20,30,40] [100,200,300,4000]' \%>\%
  reverse
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maths.R
\name{maths}
\alias{maths}
\alias{do}
\alias{do_}
\alias{lengthj}
\alias{sqrtj}
\alias{floorj}
\alias{minj}
\alias{minj_}
\alias{maxj}
\alias{maxj_}
\alias{ad}
\alias{map}
\alias{map_}
\title{Math operations}
\usage{
do(.data, ...)

do_(.data, ..., .dots)

lengthj(.data)

sqrtj(.data)

floorj(.data)

minj(.data, ...)

minj_(.data, ..., .dots)

maxj(.data, ...)

maxj_(.data, ..., .dots)

ad(.data)

map(.data, ...)

map_(.data, ..., .dots)
}
\arguments{
\item{.data}{input. This can be JSON input, or an object of class
\code{jqr} that has JSON and query params combined, which is passed
from function to function when using the jqr DSL.}

\item{...}{Comma separated list of unquoted variable names}

\item{.dots}{Used to work around non-standard evaluation}

\item{dots}{dots}
}
\description{
Math operations
}
\examples{
# do math
jq('{"a": 7}', '.a + 1')
# adding null gives back same result
jq('{"a": 7}', '.a + null')
jq('{"a": 7}', '.a += 1')
'{"a": 7}' \%>\%  do(.a + 1)
# '{"a": 7}' \%>\%  do(.a += 1) # this doesn't work quite yet
'{"a": [1,2], "b": [3,4]}' \%>\%  do(.a + .b)
'{"a": [1,2], "b": [3,4]}' \%>\%  do(.a - .b)
'{"a": 3}' \%>\%  do(4 - .a)
'["xml", "yaml", "json"]' \%>\%  do('. - ["xml", "yaml"]')
'5' \%>\%  do(10 / . * 3)
## many JSON inputs
'{"a": [1,2], "b": [3,4]} {"a": [1,5], "b": [3,10]}' \%>\%  do(.a + .b)

# comparisons
'[5,4,2,7]' \%>\% index() \%>\% do(. < 4)
'[5,4,2,7]' \%>\% index() \%>\% do(. > 4)
'[5,4,2,7]' \%>\% index() \%>\% do(. <= 4)
'[5,4,2,7]' \%>\% index() \%>\% do(. >= 4)
'[5,4,2,7]' \%>\% index() \%>\% do(. == 4)
'[5,4,2,7]' \%>\% index() \%>\% do(. != 4)
## many JSON inputs
'[5,4,2,7] [4,3,200,0.1]' \%>\% index() \%>\% do(. < 4)

# length
'[[1,2], "string", {"a":2}, null]' \%>\% index \%>\% lengthj

# sqrt
'9' \%>\% sqrtj
## many JSON inputs
'9 4 5' \%>\% sqrtj

# floor
'3.14159' \%>\% floorj
## many JSON inputs
'3.14159 30.14 45.9' \%>\% floorj

# find minimum
'[5,4,2,7]' \%>\% minj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' \%>\% minj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' \%>\% minj(foo)
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' \%>\% minj(bar)
## many JSON inputs
'[{"foo":1}, {"foo":14}] [{"foo":2}, {"foo":3}]' \%>\% minj(foo)

# find maximum
'[5,4,2,7]' \%>\% maxj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' \%>\% maxj
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' \%>\% maxj(foo)
'[{"foo":1, "bar":14}, {"foo":2, "bar":3}]' \%>\% maxj(bar)
## many JSON inputs
'[{"foo":1}, {"foo":14}] [{"foo":2}, {"foo":3}]' \%>\% maxj(foo)

# increment values
## requires special \% operators, they get escaped internally
'{"foo": 1}' \%>\% do(.foo \%+=\% 1)
'{"foo": 1}' \%>\% do(.foo \%-=\% 1)
'{"foo": 1}' \%>\% do(.foo \%*=\% 4)
'{"foo": 1}' \%>\% do(.foo \%/=\% 10)
'{"foo": 1}' \%>\% do(.foo \%//=\% 10)
### fix me - \%= doesn't work
# '{"foo": 1}' \%>\% do(.foo \%\%=\% 10)
## many JSON inputs
'{"foo": 1} {"foo": 2} {"foo": 3}' \%>\% do(.foo \%+=\% 1)

# add
'["a","b","c"]' \%>\% ad
'[1, 2, 3]' \%>\% ad
'[]' \%>\% ad
## many JSON inputs
'["a","b","c"] ["d","e","f"]' \%>\% ad

# map
## as far as I know, this only works with numbers, thus it's
## in the maths section
'[1, 2, 3]' \%>\% map(.+1)
'[1, 2, 3]' \%>\% map(./1)
'[1, 2, 3]' \%>\% map(.*4)
# many JSON inputs
'[1, 2, 3] [100, 200, 300] [1000, 2000, 30000]' \%>\% map(.+1)
}
