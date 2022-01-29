patentsview
================

> An R client to the PatentsView API

[![](http://badges.ropensci.org/112_status.svg)](https://github.com/ropensci/software-review/issues/112)
[![R-CMD-check](https://github.com/ropensci/patentsview/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/patentsview/actions)
[![CRAN
version](http://www.r-pkg.org/badges/version/patentsview)](https://cran.r-project.org/package=patentsview)

## Installation

You can get the stable version from CRAN:

``` r
install.packages("patentsview")
```

Or the development version from GitHub:

``` r
if (!"devtools" %in% rownames(installed.packages())) 
  install.packages("devtools")

devtools::install_github("ropensci/patentsview")
```

## Basic usage

The [PatentsView API](https://patentsview.org/apis/api-endpoints)
provides an interface to a disambiguated version of USPTO. The
`patentsview` R package provides one main function, `search_pv()`, to
make it easy to interact with the API:

``` r
library(patentsview)

search_pv(query = '{"_gte":{"patent_date":"2007-01-01"}}')
#> $data
#> #### A list with a single data frame on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame': 25 obs. of  3 variables:
#>   ..$ patent_id    : chr [1:25] "10000000" ...
#>   ..$ patent_number: chr [1:25] "10000000" ...
#>   ..$ patent_title : chr [1:25] "Coherent LADAR using intra-pixel quadrature "..
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_patent_count = 100,000
```

## Learning more

Head over to the package’s
[webpage](https://docs.ropensci.org/patentsview/index.html) for more
info, including:

-   A [getting started
    vignette](https://docs.ropensci.org/patentsview/articles/getting-started.html)
    for first-time users. The package was also introduced in an
    [rOpenSci blog
    post](https://ropensci.org/blog/2017/09/19/patentsview/).
-   An in-depth tutorial on [writing
    queries](https://docs.ropensci.org/patentsview/articles/writing-queries.html)
-   A list of [basic
    examples](https://docs.ropensci.org/patentsview/articles/examples.html)
-   Two examples of data applications (e.g., a brief analysis of the
    [top
    assignees](https://docs.ropensci.org/patentsview/articles/top-assignees.html)
    in the field of databases)
[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# patentsview 0.3.0 (2021-09-03)

#### Misc

* The package is now using the new HTTPS endpoints (#17)
* The list of queryable fields was updated
* `with_qfuns()` now find objects in the calling environment (@jcheng5, #20)
* Vignettes are being pre-computed (#23)
* An issue was fixed where query strings weren't being properly URL-encoded (#24)
* Adhoc logic was added to handle API throttling

# patentsview 0.2.2 (2019-01-23)

#### Misc

* Vignettes removed from package so that CRAN builds don't fail when API is down

# patentsview 0.2.1 (2018-03-05)

#### Misc

* Examples that hit the API were wrapped in `\dontrun{}` so CRAN doesn't request fixes to package when API is down

# patentsview 0.2.0 (2018-02-08)

#### New features

* `cast_pv_data()` function added to convert the data types of the data returned by `search_pv()`
* Additional fields added to the API (e.g., fields starting with `forprior_`, `examiner_`)

#### Misc

* Additional error handler added for the locations endpoint (@mustberuss, #11)
* `error_browser` option has been deprecated

# patentsview 0.1.0 (2017-05-01)

#### New functions

* `search_pv` added to send requests to the PatentsView API
* `qry_funs` list added with functions to help users write queries
* `get_fields` and `get_endpoints` added to quickly get possible field names and endpoints, respectively
* `unnest_pv_data` added to unnest the data frames in the returned data
## Test environments

* Windows Server 2019 on Github, R 4.1.1
* Windows Server 2008 on R-hub builder, R-devel
* Ubuntu 20.04.3 on Github, R 4.1.1
* Ubuntu 20.04.3 on Github, R-devel
* Debian Linux on R-hub builder, R-devel
* OS X 10.13 locally, R 3.6.2

## R CMD check results

There were no ERRORs or WARNINGs.
Notes about these vignettes:

* They aren't included in the R package by nature of the fact that the vignettes directory is included in .Rbuildignore. Instead, I just point to the site in the R documentation.
* They mostly follow the pattern described here: https://ropensci.org/blog/2019/12/08/precompute-vignettes/
  - Files ending in .Rmd.orig are the original, totally unrendered RMarkdown files
  - Files ending in .Rmd are the half-rendered RMarkdown files. Basically we pre-compute the vignettes by rendering the .Rmd.orig files to .Rmd files (see the half_render_vignettes() for details on how this happens). This creates .Rmd docs where the R outputs are already computed and embedded in the .Rmd files. These .Rmd files get further rendered during ROpenSci's automated site builds.
  - The exception to following the pattern described in the ROpenSci post is that, when the pkgdown site is built locally (and all of the artifacts from that build placed in docs/), we don't rely on the rendered vignettes that end in .Rmd. If we did, we wouldn't be able to include HTML artifacts/dependencies easily in the site. Instead, we rename files at build time so that pkgdown ends up using the Rmd.orig files instead of the .Rmd files when it builds the site. See build_site_locally() for details on how this happens.
canvas-toBlob.js
================

canvas-toBlob.js implements the standard HTML5 [`canvas.toBlob()`][1] and
`canvas.toBlobHD()` methods in browsers that do not natively support it. canvas-toBlob.js
requires `Blob` support to function, which is not present in all browsers. [Blob.js][2]
is a cross-browser `Blob` implementation that solves this.

Supported browsers
------------------

canvas-toBlob.js has [the same browser support as FileSaver.js][3].

![Tracking image](https://in.getclicky.com/212712ns.gif)

  [1]: http://www.w3.org/TR/html5/the-canvas-element.html
  [2]: https://github.com/eligrey/Blob.js
  [3]: https://github.com/eligrey/FileSaver.js#supported-browsersFileSaver.js
============

FileSaver.js implements the HTML5 W3C `saveAs()` FileSaver interface in browsers that do
not natively support it. There is a [FileSaver.js demo][1] that demonstrates saving
various media types.

FileSaver.js is the solution to saving files on the client-side, and is perfect for
webapps that need to generate files, or for saving sensitive information that shouldn't be
sent to an external server.

Looking for `canvas.toBlob()` for saving canvases? Check out
[canvas-toBlob.js][2] for a cross-browser implementation.

Supported browsers
------------------

| Browser        | Constructs as | Filenames    | Max Blob Size | Dependencies |
| -------------- | ------------- | ------------ | ------------- | ------------ |
| Firefox 20+    | Blob          | Yes          | 800 MiB       | None         |
| Firefox < 20   | data: URI     | No           | n/a           | [Blob.js](https://github.com/eligrey/Blob.js) |
| Chrome         | Blob          | Yes          | [500 MiB][3]  | None         |
| Chrome for Android | Blob      | Yes          | [500 MiB][3]  | None         |
| IE 10+         | Blob          | Yes          | 600 MiB       | None         |
| Opera 15+      | Blob          | Yes          | 500 MiB       | None         |
| Opera < 15     | data: URI     | No           | n/a           | [Blob.js](https://github.com/eligrey/Blob.js) |
| Safari 6.1+*   | Blob          | No           | ?             | None         |
| Safari < 6     | data: URI     | No           | n/a           | [Blob.js](https://github.com/eligrey/Blob.js) |

Feature detection is possible:

```js
try {
    var isFileSaverSupported = !!new Blob;
} catch (e) {}
```

### IE < 10

It is possible to save text files in IE < 10 without Flash-based polyfills.
See [ChenWenBrian and koffsyrup's `saveTextAs()`](https://github.com/koffsyrup/FileSaver.js#examples) for more details.

### Safari 6.1+

Blobs may be opened instead of saved sometimes—you may have to direct your Safari users to manually
press <kbd>⌘</kbd>+<kbd>S</kbd> to save the file after it is opened. Using the `application/octet-stream` MIME type to force downloads [can cause issues in Safari](https://github.com/eligrey/FileSaver.js/issues/12#issuecomment-47247096).

### iOS

saveAs must be run within a user interaction event such as onTouchDown or onClick; setTimeout will prevent saveAs from triggering. Due to restrictions in iOS saveAs opens in a new window instead of downloading, if you want this fixed please [tell Apple](https://bugs.webkit.org/show_bug.cgi?id=102914) how this bug is affecting you.

Syntax
------

```js
FileSaver saveAs(Blob data, DOMString filename, optional Boolean disableAutoBOM)
```

Pass `true` for `disableAutoBOM` if you don't want FileSaver.js to automatically provide Unicode text encoding hints (see: [byte order mark](https://en.wikipedia.org/wiki/Byte_order_mark)).

Examples
--------

### Saving text

```js
var blob = new Blob(["Hello, world!"], {type: "text/plain;charset=utf-8"});
saveAs(blob, "hello world.txt");
```

The standard W3C File API [`Blob`][4] interface is not available in all browsers.
[Blob.js][5] is a cross-browser `Blob` implementation that solves this.

### Saving a canvas

```js
var canvas = document.getElementById("my-canvas"), ctx = canvas.getContext("2d");
// draw to canvas...
canvas.toBlob(function(blob) {
    saveAs(blob, "pretty image.png");
});
```

Note: The standard HTML5 `canvas.toBlob()` method is not available in all browsers.
[canvas-toBlob.js][6] is a cross-browser `canvas.toBlob()` that polyfills this.


![Tracking image](https://in.getclicky.com/212712ns.gif)

  [1]: http://eligrey.com/demos/FileSaver.js/
  [2]: https://github.com/eligrey/canvas-toBlob.js
  [3]: https://code.google.com/p/chromium/issues/detail?id=375297
  [4]: https://developer.mozilla.org/en-US/docs/DOM/Blob
  [5]: https://github.com/eligrey/Blob.js
  [6]: https://github.com/eligrey/canvas-toBlob.js

Contributing
------------

The `FileSaver.js` distribution file is compiled with Uglify.js like so:

```bash
uglifyjs FileSaver.js --mangle --comments /@source/ > FileSaver.min.js
```

Please make sure you build a production version before submitting a pull request.

Bower Installation
------------------

Please see the [this repo](http://github.com/Teleborder/FileSaver.js) for a bower-compatible fork of FileSaver.js, available under the package name `file-saver.js`.
Blob.js
==============

Blob.js implements the W3C [`Blob`][1] interface in browsers that do
not natively support it.

![Tracking image](https://in.getclicky.com/212712ns.gif)

  [1]: https://developer.mozilla.org/en-US/docs/Web/API/Blob
