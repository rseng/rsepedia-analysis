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
---
title: "patentsview"
output: github_document
---

```{r, echo = FALSE, message = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

> An R client to the PatentsView API

[![](http://badges.ropensci.org/112_status.svg)](https://github.com/ropensci/software-review/issues/112)
[![R-CMD-check](https://github.com/ropensci/patentsview/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/patentsview/actions)
[![CRAN version](http://www.r-pkg.org/badges/version/patentsview)](https://cran.r-project.org/package=patentsview)

## Installation

You can get the stable version from CRAN:

```{r eval = FALSE}
install.packages("patentsview")
```

Or the development version from GitHub:

```{r eval = FALSE}
if (!"devtools" %in% rownames(installed.packages())) 
  install.packages("devtools")

devtools::install_github("ropensci/patentsview")
```

## Basic usage

The [PatentsView API](https://patentsview.org/apis/api-endpoints) provides an interface to a disambiguated version of USPTO. The `patentsview` R package provides one main function, `search_pv()`, to make it easy to interact with the API:

```{r}
library(patentsview)

search_pv(query = '{"_gte":{"patent_date":"2007-01-01"}}')
```

## Learning more

Head over to the package's [webpage](https://docs.ropensci.org/patentsview/index.html) for more info, including:

* A [getting started vignette](https://docs.ropensci.org/patentsview/articles/getting-started.html) for first-time users. The package was also introduced in an [rOpenSci blog post](https://ropensci.org/blog/2017/09/19/patentsview/).
* An in-depth tutorial on [writing queries](https://docs.ropensci.org/patentsview/articles/writing-queries.html)
* A list of [basic examples](https://docs.ropensci.org/patentsview/articles/examples.html)
* Two examples of data applications (e.g., a brief analysis of the [top assignees](https://docs.ropensci.org/patentsview/articles/top-assignees.html) in the field of databases)
---
title: "patentsview"
output: github_document
---

```{r, echo = FALSE, message = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

<br>

## Installation

You can get the stable version from CRAN:

```{r eval = FALSE}
install.packages("patentsview")
```

Or the development version from GitHub:

```{r eval = FALSE}
if (!"devtools" %in% rownames(installed.packages())) 
  install.packages("devtools")

devtools::install_github("ropensci/patentsview")
```

## Basic usage

The [PatentsView API](http://www.patentsview.org/api/doc.html) provides an interface to a disambiguated version of USPTO. The `patentsview` R package provides one main function, `search_pv()`, to make it easy to interact with the API:

```{r}
library(patentsview)

search_pv(query = '{"_gte":{"patent_date":"2007-01-01"}}')
```

## Learning more

Check out the:

* [Getting started vignette](articles/getting-started.html) or [rOpenSci blog post](https://ropensci.org/blog/blog/2017/09/19/patentsview) if you are a first-time user
* In-depth tutorial on [writing queries](articles/writing-queries.html)
* List of [basic examples](articles/examples.html)
* Examples of data applications (e.g., a brief analysis of the [top assignees](articles/top-assignees.html) in the field of databases)
---
title: "Top assignees"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Top assignees}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



The following is a quick analysis of the top organizations patenting in the field of databases.

1. The first step is to download the relevant data fields from the PatentsView API:


```r
library(patentsview)
library(dplyr)
library(highcharter)
library(DT)
library(knitr)

# We first need to write a query. Our query will look for "database" in either 
# the patent title or abstract...Note, this isn't a terribly good way to ID our 
# patents, but it will work for the purpose of demonstration. Users who are 
# interested in writing higher-quality queries could consult the large body of 
# research that has been done in field of patent document retrieval.
query <- with_qfuns(
  or(
    text_phrase(patent_abstract = "database"),
    text_phrase(patent_title = "database")
  )
)

query
#> {"_or":[{"_text_phrase":{"patent_abstract":"database"}},{"_text_phrase":{"patent_title":"database"}}]}

# Create a list of the fields we'll need for the analysis
fields <- c(
  "patent_number", "assignee_organization",
  "patent_num_cited_by_us_patents", "app_date", "patent_date",
  "assignee_total_num_patents"
)

# Send an HTTP request to the PatentsView API to get the data
pv_out <- search_pv(query, fields = fields, all_pages = TRUE)
```

2. Now let's identify who the top assignees are based on how many patents they have in our data set. We'll also calculate how many total patents these assignees have and what fraction of their total patents relate to databases.


```r
# Unnest the data frames that are stored in the assignee list column
dl <- unnest_pv_data(pv_out$data, "patent_number")
dl
#> List of 3
#>  $ assignees   :'data.frame':	67574 obs. of  4 variables:
#>   ..$ patent_number             : chr [1:67574] "10000911" ...
#>   ..$ assignee_organization     : chr [1:67574] "DOOSAN INFRACORE CO., LTD." ...
#>   ..$ assignee_total_num_patents: chr [1:67574] "189" ...
#>   ..$ assignee_key_id           : chr [1:67574] "10886" ...
#>  $ applications:'data.frame':	66655 obs. of  3 variables:
#>   ..$ patent_number: chr [1:66655] "10000911" ...
#>   ..$ app_date     : chr [1:66655] "2014-12-05" ...
#>   ..$ app_id       : chr [1:66655] "15/101707" ...
#>  $ patents     :'data.frame':	66655 obs. of  3 variables:
#>   ..$ patent_number                 : chr [1:66655] "10000911" ...
#>   ..$ patent_num_cited_by_us_patents: chr [1:66655] "0" ...
#>   ..$ patent_date                   : chr [1:66655] "2018-06-19" ...

# Create a data frame with the top 75 assignees:
top_asgns <-
  dl$assignees %>%
    filter(!is.na(assignee_organization)) %>% # some patents are assigned to an inventor (not an org)
    mutate(ttl_pats = as.numeric(assignee_total_num_patents)) %>%
    group_by(assignee_organization, ttl_pats) %>% 
    summarise(db_pats = n()) %>% 
    mutate(frac_db_pats = round(db_pats / ttl_pats, 3)) %>%
    ungroup() %>%
    select(c(1, 3, 2, 4)) %>%
    arrange(desc(db_pats)) %>%
    slice(1:75)

# Create datatable
datatable(
  data = top_asgns,
  rownames = FALSE,
  colnames = c(
    "Assignee", "DB patents","Total patents", "DB patents / total patents"
  ),
  caption = htmltools::tags$caption(
    style = 'caption-side: top; text-align: left; font-style: italic;',
    "Table 1: Top assignees in 'databases'"
  ),
  options = list(pageLength = 10)
)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

<br>

IBM is far and away the biggest player in the field. However, we can see that Oracle and Salesforce.com are relatively more interested in this area, as indicated by the fraction of their patents that relate to databases.

3. Let's see how these assignees' level of investment in databases has changed over time.


```r
# Create a data frame with patent counts by application year for each assignee
data <- 
  top_asgns %>%
    select(-contains("pats")) %>%
    slice(1:5) %>%
    inner_join(dl$assignees) %>%
    inner_join(dl$applications) %>%
    mutate(app_yr = as.numeric(substr(app_date, 1, 4))) %>%
    group_by(assignee_organization, app_yr) %>%
    count() 

# Plot the data using highchartr:
hchart(
  data, "line", 
  hcaes(x = app_yr, y = n, group = assignee_organization)
) %>%
  hc_plotOptions(series = list(marker = list(enabled = FALSE))) %>%
  hc_xAxis(title = list(text = "Application year")) %>%
  hc_yAxis(title = list(text = "Patents")) %>%
  hc_title(text = "Top five assignees in 'databases'") %>%
  hc_subtitle(text = "Yearly patent applications over time")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

It's hard to see any clear trends in this graph. What is clear is that the top assignees have all been patenting in the field for many years.

4. Finally, let's see how the organizations compare in terms of their citation rates. First, we'll need to normalize the raw citation counts by publication year, so that older patents don't have an unfair advantage over younger patents (i.e., because they have had a longer time to accumulate citations).


```r
# Write a ranking function that will be used to rank patents by their citation counts
percent_rank2 <- function(x)
  (rank(x, ties.method = "average", na.last = "keep") - 1) / (sum(!is.na(x)) - 1)

# Create a data frame with normalized citation rates and stats from Step 2
asng_p_dat <-
  dl$patents %>%
    mutate(patent_yr = substr(patent_date, 1, 4)) %>%
    group_by(patent_yr) %>%
    mutate(perc_cite = percent_rank2(patent_num_cited_by_us_patents)) %>%
    inner_join(dl$assignees) %>%
    group_by(assignee_organization) %>%
    summarise(mean_perc = mean(perc_cite)) %>%
    inner_join(top_asgns) %>%
    arrange(desc(ttl_pats)) %>%
    filter(!is.na(assignee_organization)) %>%
    slice(1:20) %>%
    mutate(color = "#f1c40f") %>%
    as.data.frame()

kable(head(asng_p_dat), row.names = FALSE)
```



|assignee_organization                       | mean_perc| db_pats| ttl_pats| frac_db_pats|color   |
|:-------------------------------------------|---------:|-------:|--------:|------------:|:-------|
|International Business Machines Corporation | 0.4725093|    6284|   144597|        0.043|#f1c40f |
|Samsung Electronics Co., Ltd                | 0.4301865|     384|    98404|        0.004|#f1c40f |
|Canon Kabushiki Kaisha                      | 0.4299559|     217|    78447|        0.003|#f1c40f |
|SONY CORPORATION                            | 0.4214531|     429|    55766|        0.008|#f1c40f |
|Kabushiki Kaisha Toshiba                    | 0.4460499|     246|    52964|        0.005|#f1c40f |
|General Electric Company                    | 0.4908506|     305|    48821|        0.006|#f1c40f |

Now let's visualize the data. Each assignee will be represented by a point/bubble. The x-value of the point will represent the total number of patents the assignee has published in the field of databases (on a log scale), while the y-value will represent its average normalized citation rate. The size of the bubble will be proportional to the percent of the assignee's patents that relate to databases.


```r
# Adapted from http://jkunst.com/highcharter/showcase.html
hchart(
  asng_p_dat, "scatter", 
  hcaes(x = db_pats, y = mean_perc, size = frac_db_pats, 
        group = assignee_organization, color = color)
) %>%
  hc_xAxis(
    title = list(text = "DB patents"), type = "logarithmic",
    allowDecimals = FALSE, endOnTick = TRUE
  ) %>%
  hc_yAxis(title = list(text = "Mean cite perc.")) %>%
  hc_title(text = "Top assignees in 'databases'") %>%
  hc_add_theme(hc_theme_flatdark()) %>%
  hc_tooltip(
    useHTML = TRUE, pointFormat = tooltip_table(
    x = c("DB patents", "Mean cite percentile", "Fraction DB patents"),
    y = c("{point.db_pats:.0f}","{point.mean_perc:.2f}", "{point.frac_db_pats:.3f}")
  )) %>%
  hc_legend(enabled = FALSE)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

<br>

It looks like Microsoft has relatively high values across all three three metrics (average citation percentile, number of database patents, and percent of total patents that are related to databases). IBM has more patents than Microsoft, but also has a lower average citation percentile.
---
title: "Writing queries"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Writing queries}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Three ways to write the same query

Let's say you want to find all patents published in the last 10 years that have the word "dog" in their titles or abstracts, and whose assignees are located in either the US or Canada. Here are three ways you could write such a query: 

1. Use a string:


```r
query_v_1 <-
  '{"_and":[
          {"_gte":{"patent_date":"2007-03-01"}},
          {"_or":[
            {"_text_all":{"patent_title":"dog"}},
            {"_text_all":{"patent_abstract":"dog"}}
          ]},
          {"_or":[
            {"_eq":{"assingee_country":"US"}},
            {"_eq":{"assingee_country":"CA"}}
          ]}
  ]}'
```

2. Use a list:


```r
query_v_2 <- 
  list("_and" = 
       list(
          list("_gte" = list(patent_date = "2007-03-01")),
          list("_or" = 
                 list(
                   list("_text_all" = list(patent_title = "dog")),
                   list("_text_all" = list(patent_abstract = "dog"))
                   )
               ),
          list("_or" = 
                 list(
                   list("_eq" = list(assingee_country = "US")),
                   list("_eq" = list(assingee_country = "CA"))
                   )
               )
      )
  )
```

3. Use the `patentsview` domain specific language (DSL): 


```r
library(patentsview)

query_v_3 <- 
  with_qfuns(
    and(
      gte(patent_date = "2007-03-01"),
      or(
        text_all(patent_title = "dog"),
        text_all(patent_abstract = "dog")
      ),
      eq(assingee_country = c("US", "CA"))
    )
  )
```

## Why use the DSL?

We can see that all three versions of the query shown above are equivalent:


```r
jsonlite::minify(query_v_1)
#> {"_and":[{"_gte":{"patent_date":"2007-03-01"}},{"_or":[{"_text_all":{"patent_title":"dog"}},{"_text_all":{"patent_abstract":"dog"}}]},{"_or":[{"_eq":{"assingee_country":"US"}},{"_eq":{"assingee_country":"CA"}}]}]}
jsonlite::toJSON(query_v_2, auto_unbox = TRUE)
#> {"_and":[{"_gte":{"patent_date":"2007-03-01"}},{"_or":[{"_text_all":{"patent_title":"dog"}},{"_text_all":{"patent_abstract":"dog"}}]},{"_or":[{"_eq":{"assingee_country":"US"}},{"_eq":{"assingee_country":"CA"}}]}]}
jsonlite::toJSON(query_v_3, auto_unbox = TRUE)
#> {"_and":[{"_gte":{"patent_date":"2007-03-01"}},{"_or":[{"_text_all":{"patent_title":"dog"}},{"_text_all":{"patent_abstract":"dog"}}]},{"_or":[{"_eq":{"assingee_country":"US"}},{"_eq":{"assingee_country":"CA"}}]}]}
```

...So why would you ever want to use method 3 over methods 1 or 2? There are two main reasons:

#### 1. Query validation 

`search_pv()` will check your query for errors if you use methods 2 or 3. This is not the case for method 1, where you would have to rely on the API's error messages for guidance if your query is invalid. `search_pv()` checks queries for the following:

* The fields included in your query are _queryable_ for the endpoint (i.e., the field can be used in the user query). For example, it would make sure that `assingee_country` can be used in the query argument if you sent the above query to the patents endpoint.
* The fields in your query are compatible with the comparison operators you used. For example, it would confirm that the `text_all` operator was used with a field whose type was "full text" (`patent_title` above).
* You supplied the correct value type for the field (e.g., `patent_date` is a character, not an integer).

#### 2. Concise, easy to use syntax for complex queries

Methods 1 and 3 are both shorter than method 2, making them quicker. It's also a lot easier to get the JSON syntax correct when using method 3 compared to method 1, because you don't have to write any JSON at all using the DSL...This is important because the API is fairly picky about the query syntax, so it's not trivial to get it correct. For example, the API will throw an error if you use a box in your JSON when is not absolutely necessary, even if your query is still valid JSON (e.g., `query = {"_gte":{"patent_date":["2007-03-01"]}}` will throw an error). 

Compared to method 1, method 3 will correctly "or" together values if you put them in a vector. For example, in the query shown above, a vector of two values was given for `assingee_country` (`c("US", "CA")`). This safely converted the single "equals" statement in the third element of the query (`eq(assingee_country = c("US", "CA"))`) to two separate equals statements that got or'd together.[^1]

## Basics of the language

All of the functions that make up the DSL are found in the `qry_funs` list (e.g., `qry_funs$eq()`). You can evaluate code in the context of this list using the function `with_qfuns()` (see `?with_qfuns()` for an example that demonstrates how `with_qfuns()` can save you typing). There are three types of functions in `qry_funs`:

1. **Comparison operator functions** (`eq`, `neq`, `gt`, `gte`, `lt`, `lte`, `begins`, `contains`, `text_all`, `text_any`, `text_phrase`). These functions are used to compare a field to a value. For example, using the "less than or equal to" function (`lte`), we can filter out patents published after some date (e.g., `query = qry_funs$lte(patent_date = "2001-01-05")`). See the "comparison operators" section of the API's [query language page](http://www.patentsview.org/api/query-language.html) for a description of the 11 comparison operators. One important thing to keep in mind is that certain comparison operators only work with certain data types. For example, you can't use the `begins` function on `patent_abstract` because `patent_abstract` is of data type "full text" and `begins` only works with fields of data type "string." 
2. **Array functions** (`and` and `or`). You can use these functions to logically combine your calls to the comparison operators. For example, we can require that the patent date is less than or equal to 2001-01-05 _and_ the inventor's last name is "Ihaka" (`query = with_qfuns(and(lte(patent_date = "2001-01-05"), eq(inventor_last_name = "Ihaka")))`).
3. **not function** (`not`). This function negates a comparison. For example, we could search for patents that don't have the word "hi" in their titles like this: `qry_funs$not(qry_funs$text_phrase(patent_title = "hi"))`.

## Query examples 

*The following queries are intended for the patents endpoint*

Patents linked to an assignee with 10 or fewer distinct (and disambiguated) inventors:


```r
qry_funs$lte(assignee_total_num_inventors = 10)
#> {"_lte":{"assignee_total_num_inventors":10}}
```

Patents assigned to the "CPC subsection"[^2] of G12 (physics instruments):


```r
qry_funs$eq(cpc_subsection_id = "G12")
#> {"_eq":{"cpc_subsection_id":"G12"}}
```

Patents that:

* Have an inventor listed on them whose first name contains "joh" AND
* Have an abstract containing either the phrase "dog bark" or "cat meow" AND
* Have an abstract that doesn't have the phrase "dog chain" in it:


```r
with_qfuns(
  and(
    contains(rawinventor_first_name = "joh"),
    text_phrase(patent_abstract = c("dog bark", "cat meow")),
    not(
      text_phrase(patent_abstract = c("dog chain"))
    )
  )
)
#> {"_and":[{"_contains":{"rawinventor_first_name":"joh"}},{"_or":[{"_text_phrase":{"patent_abstract":"dog bark"}},{"_text_phrase":{"patent_abstract":"cat meow"}}]},{"_not":{"_text_phrase":{"patent_abstract":"dog chain"}}}]}
```

Patents that:

* Have an inventor listed on them whose last name is “Smith” AND
* Have “cotton gin” in their title 

OR

* Have an inventor listed on them whose last name is “Hopper” AND
* Have “COBOL” in their title 


```r
with_qfuns(
  or(
    and(
      eq(inventor_last_name = "smith"),
      text_phrase(patent_title = "cotton gin")
    ),
    and(
      eq(inventor_last_name = "hopper"),
      text_phrase(patent_title = "COBOL")
    )
  )
)
#> {"_or":[{"_and":[{"_eq":{"inventor_last_name":"smith"}},{"_text_phrase":{"patent_title":"cotton gin"}}]},{"_and":[{"_eq":{"inventor_last_name":"hopper"}},{"_text_phrase":{"patent_title":"COBOL"}}]}]}
```

[^1]: One may note that using "value arrays" is supposedly supported by the API. For example, the API documentation gives the following query as an example of their use: `'{"inventor_last_name":["Whitney","Hopper"]}'`. The problem with this is that the API is not consistent in its handling of value arrays. For many of the comparison operators, one cannot "or" together values like this using arrays. Thus, the DSL in `patentsview` never relies on arrays when creating queries.
[^2]: PatentsView gets the names of the CPC hierarchy wrong. For example, a "CPC subsection" according to PatentsView is actually a CPC class.
---
title: "Examples"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Patents endpoint

Which patents have been cited by more than 500 US patents?


```r
library(patentsview)

search_pv(query = qry_funs$gt(patent_num_cited_by_us_patents = 500))
#> $data
#> #### A list with a single data frame on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame':	25 obs. of  3 variables:
#>   ..$ patent_id    : chr [1:25] "3940844" ...
#>   ..$ patent_number: chr [1:25] "3940844" ...
#>   ..$ patent_title : chr [1:25] "Method of installing an insulating sleeve on"..
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_patent_count = 7,915
```

How many distinct inventors are represented by these highly-cited patents?


```r
# Setting subent_cnts = TRUE will give us the subentity counts. Since inventors 
# are subentities for the patents endpoint, this means we will get their counts.
search_pv(
  query = qry_funs$gt(patent_num_cited_by_us_patents = 500),
  fields = c("patent_number", "inventor_id"), 
  subent_cnts = TRUE
)
#> $data
#> #### A list with a single data frame (with list column(s) inside) on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame':	25 obs. of  2 variables:
#>   ..$ patent_number: chr [1:25] "3940844" ...
#>   ..$ inventors    :List of 25
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_patent_count = 7,915, total_inventor_count = 11,263
```

Where geographically have Microsoft inventors been coming from over the past 20 years?


```r
# Write the query
query <- with_qfuns(
  and(
    gte(patent_date = "2007-07-25"), # Dates are in yyyy-mm-dd format 
    contains(assignee_organization = "microsoft")
  )
)

# Create a field list
inv_fields <- get_fields(endpoint = "patents", groups = "inventors")
fields <- c(inv_fields, "patent_number")

# Pull the data
pv_out <- search_pv(query, fields = fields, all_pages = TRUE)

# Unnest the inventor list column
unnest_pv_data(pv_out$data, "patent_number")
#> List of 2
#>  $ inventors:'data.frame':	144495 obs. of  24 variables:
#>   ..$ patent_number                 : chr [1:144495] "10001683" ...
#>   ..$ inventor_city                 : chr [1:144495] "Mountain View" ...
#>   ..$ inventor_country              : chr [1:144495] "US" ...
#>   ..$ inventor_county               : logi [1:144495] NA ...
#>   ..$ inventor_county_fips          : chr [1:144495] "6085" ...
#>   ..$ inventor_first_name           : chr [1:144495] "Andriy" ...
#>   ..$ inventor_first_seen_date      : chr [1:144495] "2014-04-22" ...
#>   ..$ inventor_id                   : chr [1:144495] "fl:a_ln:pletenetskyy-1" ..
#>   ..$ inventor_last_name            : chr [1:144495] "Pletenetskyy" ...
#>   ..$ inventor_last_seen_date       : chr [1:144495] "2020-10-13" ...
#>   ..$ inventor_lastknown_city       : chr [1:144495] "Mountain View" ...
#>   ..$ inventor_lastknown_country    : chr [1:144495] "US" ...
#>   ..$ inventor_lastknown_latitude   : chr [1:144495] "37.4139" ...
#>   ..$ inventor_lastknown_location_id: chr [1:144495] "37.4139|-122.085" ...
#>   ..$ inventor_lastknown_longitude  : chr [1:144495] "-122.085" ...
#>   ..$ inventor_lastknown_state      : chr [1:144495] "CA" ...
#>   ..$ inventor_latitude             : chr [1:144495] "37.4139" ...
#>   ..$ inventor_location_id          : chr [1:144495] "37.4139|-122.085" ...
#>   ..$ inventor_longitude            : chr [1:144495] "-122.085" ...
#>   ..$ inventor_sequence             : chr [1:144495] "0" ...
#>   ..$ inventor_state                : chr [1:144495] "CA" ...
#>   ..$ inventor_state_fips           : chr [1:144495] "06" ...
#>   ..$ inventor_total_num_patents    : chr [1:144495] "15" ...
#>   ..$ inventor_key_id               : chr [1:144495] "489570" ...
#>  $ patents  :'data.frame':	39051 obs. of  1 variable:
#>   ..$ patent_number: chr [1:39051] "10001683" ...
```

## Inventors endpoint

Which inventors have Chicago, IL listed as their location on at least one patent.[^1]


```r
search_pv(
  query = '{"_and":[{"location_city":"Chicago"},{"location_state":"IL"}]}',
  endpoint = "inventors"
)
#> $data
#> #### A list with a single data frame on an inventor level:
#> 
#> List of 1
#>  $ inventors:'data.frame':	25 obs. of  3 variables:
#>   ..$ inventor_id        : chr [1:25] "fl:b_ln:gunderson-2" ...
#>   ..$ inventor_first_name: chr [1:25] "Bjorn" ...
#>   ..$ inventor_last_name : chr [1:25] "Gunderson" ...
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_inventor_count = 16,360
```

## Assignees endpoint

Which assignees have an interest in beer?


```r
search_pv(
  query = qry_funs$text_phrase(patent_title = "beer"), 
  endpoint = "assignees"
)
#> $data
#> #### A list with a single data frame on an assignee level:
#> 
#> List of 1
#>  $ assignees:'data.frame':	25 obs. of  4 variables:
#>   ..$ assignee_id          : chr [1:25] "6cdc5a48-9dd2-4d1c-8c35-27eaef76ffd1"..
#>   ..$ assignee_first_name  : logi [1:25] NA ...
#>   ..$ assignee_last_name   : logi [1:25] NA ...
#>   ..$ assignee_organization: chr [1:25] "Rohm Co., Ltd." ...
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_assignee_count = 225
```

[^1]: Example taken from http://www.patentsview.org/api/inventor.html.
---
title: "Getting started"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## A basic example 

Let's start with a basic example of how to use the package's primary function, `search_pv()`:


```r
library(patentsview)

search_pv(
  query = '{"_gte":{"patent_date":"2007-01-01"}}',
  endpoint = "patents"
)
#> $data
#> #### A list with a single data frame on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame':	25 obs. of  3 variables:
#>   ..$ patent_id    : chr [1:25] "10000000" ...
#>   ..$ patent_number: chr [1:25] "10000000" ...
#>   ..$ patent_title : chr [1:25] "Coherent LADAR using intra-pixel quadrature "..
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_patent_count = 100,000
```

This call to `search_pv()` sends our query to the patents endpoint (the default). The API has 7 different endpoints, corresponding to 7 different entity types (assignees, CPC subsections, inventors, locations, NBER subcategories, patents, and USPC main classes).[^1] Your choice of endpoint determines which entity your query is applied to, as well as the structure of the data that is returned (more on this in the "7 endpoints for 7 entities section"). For now, let's turn our attention to the `query` parameter. 

## Writing queries

The PatentsView query syntax is documented on their [query language page](https://patentsview.org/apis/api-query-language).[^2] However, it can be difficult to get your query right if you're writing it by hand (i.e., just writing the query in a string like `'{"_gte":{"patent_date":"2007-01-01"}}'`, as we did in the example shown above). The `patentsview` package comes with a simple domain specific language (DSL) to make writing queries a breeze. I recommend using the functions in this DSL for all but the most basic queries, especially if you're encountering errors and don't understand why. To get a feel for how it works, let's rewrite the query shown above using one of the functions in the DSL, `qry_funs$gte()`:


```r
qry_funs$gte(patent_date = "2007-01-01")
#> {"_gte":{"patent_date":"2007-01-01"}}
```

More complex queries are also possible:


```r
with_qfuns(
  and(
    gte(patent_date = "2007-01-01"),
    text_phrase(patent_abstract = c("computer program", "dog leash"))
  )
)
#> {"_and":[{"_gte":{"patent_date":"2007-01-01"}},{"_or":[{"_text_phrase":{"patent_abstract":"computer program"}},{"_text_phrase":{"patent_abstract":"dog leash"}}]}]}
```

Check out the [writing queries vignette](https://docs.ropensci.org/patentsview/articles/writing-queries.html) for more details on using the DSL.

## Fields

Each endpoint has a different set of _queryable_ and _retrievable_ fields. Queryable fields are those that you can include in your query (e.g., `patent_date` shown in the first example). Retrievable fields are those that you can get data on (i.e., fields returned by `search_pv()`). In the first example, we didn't specify which fields we wanted to retrieve so we were given the default set. You can specify which fields you want using the `fields` argument:


```r
search_pv(
  query = '{"_gte":{"patent_date":"2007-01-01"}}',
  fields = c("patent_number", "patent_title")
)
#> $data
#> #### A list with a single data frame on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame':	25 obs. of  2 variables:
#>   ..$ patent_number: chr [1:25] "10000000" ...
#>   ..$ patent_title : chr [1:25] "Coherent LADAR using intra-pixel quadrature "..
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_patent_count = 100,000
```

To list all of the retrievable fields for a given endpoint, use `get_fields()`:


```r
retrvble_flds <- get_fields(endpoint = "patents")
head(retrvble_flds)
#> [1] "appcit_app_number" "appcit_category"   "appcit_date"      
#> [4] "appcit_kind"       "appcit_sequence"   "app_country"
```

You can also visit an endpoint's online documentation page to see a list of its queryable and retrievable fields (e.g., see the [inventor field list table](http://www.patentsview.org/api/inventor.html#field_list)). Note the "Query" column in this table, which indicates whether the field is both queryable and retrievable (Query = Y), or just retrievable (Query = N). The field tables for all of the endpoints can be found in the `fieldsdf` data frame, which you can load using `data("fieldsdf")` or `View(patentsview::fieldsdf)`.

**An important note: By default, PatentsView uses disambiguted versions of assignees, inventors, and locations, instead of raw data.** For example, let's say you search for all inventors whose first name is "john." The PatentsView API is going to return all of the inventors who have a preferred first name (as per the disambiguation results) of john, which may not necessarily be their raw first name. You could be getting back inventors whose first name appears on the patent as, say, "jonathan," "johnn," or even "john jay." You can search on the raw inventor names instead of the preferred names by using the fields starting with "raw" in your query (e.g., `rawinventor_first_name`). The assignee and location raw data fields are not currently being offered by the API. To see the methods behind the disambiguation process, see the [PatentsView Inventor Disambiguation Technical Workshop website]( http://www.patentsview.org/workshop/).

## Paginated responses

By default, `search_pv()` returns 25 records per page and only gives you the first page of results. I suggest sticking with these defaults while you're figuring out the details of your request, such as the query you want to use and the fields you want returned. Once you have those items finalized, you can use the `per_page` argument to download up to 10,000 records per page. You can also choose which page of results you want with the `page` argument:


```r
search_pv(
  query = qry_funs$eq(inventor_last_name = "chambers"),
  page = 2, per_page = 150 # gets records 150 - 300
) 
#> $data
#> #### A list with a single data frame on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame':	150 obs. of  3 variables:
#>   ..$ patent_id    : chr [1:150] "10577927" ...
#>   ..$ patent_number: chr [1:150] "10577927" ...
#>   ..$ patent_title : chr [1:150] "Mud pulse telemetry tool comprising a low t"..
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_patent_count = 2,208
```

You can download all pages of output in one call by setting `all_pages = TRUE`. This will set `per_page` equal to 10,000 and loop over all pages of output (downloading up to 10 pages, or 100,000 records total):


```r
search_pv(
  query = qry_funs$eq(inventor_last_name = "chambers"),
  all_pages = TRUE
)
#> $data
#> #### A list with a single data frame on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame':	2208 obs. of  3 variables:
#>   ..$ patent_id    : chr [1:2208] "10000988" ...
#>   ..$ patent_number: chr [1:2208] "10000988" ...
#>   ..$ patent_title : chr [1:2208] "Seal assemblies in subsea rotating control"..
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_patent_count = 2,208
```

## Entity counts

Our last two calls to `search_pv()` gave the same value for `total_patent_count`, even though we got a lot more data from the second call. This is because the entity counts returned by the API refer to the number of distinct entities across all *downloadable pages of output*, not just the page that was returned. Downloadable pages of output is an important phrase here, as the API limits us to 100,000 records per query. For example, we got `total_patent_count = 100,000` when we searched for patents published on or after 2007, even though there are way more than 100,000 such patents. See the FAQs below for details on how to overcome the 100,000 record restriction.

## 7 endpoints for 7 entities

We can get similar data from the 7 endpoints. For example, the following two calls differ only in the endpoint that is chosen:


```r
query <- qry_funs$eq(inventor_last_name = "chambers")
fields <- c("patent_number", "inventor_last_name", "assignee_organization")

# Here we are using the patents endpoint:
search_pv(query, endpoint = "patents", fields = fields)
#> $data
#> #### A list with a single data frame (with list column(s) inside) on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame':	25 obs. of  3 variables:
#>   ..$ patent_number: chr [1:25] "10000988" ...
#>   ..$ inventors    :List of 25
#>   ..$ assignees    :List of 25
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_patent_count = 2,208
```


```r
# While here we are using the assignees endpoint:
search_pv(query, endpoint = "assignees", fields = fields)
#> $data
#> #### A list with a single data frame (with list column(s) inside) on an assignee level:
#> 
#> List of 1
#>  $ assignees:'data.frame':	25 obs. of  3 variables:
#>   ..$ assignee_organization: chr [1:25] "XEROX CORPORATION" ...
#>   ..$ patents              :List of 25
#>   ..$ inventors            :List of 25
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_assignee_count = 531
```

Your choice of endpoint determines two things:

1. **Which entity your query is applied to.** The first call shown above used the patents endpoint, so the API searched for patents that have at least one inventor listed on them with the last name "chambers." The second call used the assignees endpoint, so the API searched for all assignees that have been *assigned to at least one patent* which has an inventor listed on it with the last name "chambers."

2. **The structure of the data frame that is returned.** The first call returned a data frame on the patent level, meaning that each row corresponded to a different patent. Fields that were not on the patent level (e.g., `inventor_last_name`) were returned in list columns that are named after the entity associated with the field (e.g., the `inventors` entity).[^3] Meanwhile, the second call gave us a data frame on the assignee level (one row for each assignee) because it used the assignees endpoint.

Most of the time you will want to use the patents endpoint. Note that you can still effectively filter on fields that are not at the patent-level when using the patents endpoint (e.g., you can filter on assignee name or CPC category). This is because patents are relatively low-level entities. For higher level entities like assignees, if you filter on a field that is not at the assignee-level (e.g., inventor name), the API will return data on any assignee that has at least one inventor whose name matches your search, which is probably not what you want.

## Casting fields

The API always returns the data fields as strings, even if they would be better stored using a different data type (e.g., numeric). You can cast all fields to their preferred R types using `cast_pv_data()`:


```r
res <- search_pv(
  query = "{\"patent_number\":\"5116621\"}", 
  fields = c("patent_date", "patent_title", "patent_year")
)

# Right now all of the fields are stored as character vectors:
res
#> $data
#> #### A list with a single data frame on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame':	1 obs. of  3 variables:
#>   ..$ patent_date : chr "1992-05-26"
#>   ..$ patent_title: chr "Anti-inflammatory analgesic patch"
#>   ..$ patent_year : chr "1992"
#> 
#> $query_results
#> #### Distinct entity counts across all downloadable pages of output:
#> 
#> total_patent_count = 1

# Use more appropriate data types:
cast_pv_data(res$data)
#> #### A list with a single data frame on a patent level:
#> 
#> List of 1
#>  $ patents:'data.frame':	1 obs. of  3 variables:
#>   ..$ patent_date : Date[1:1], format: "1992-05-26"
#>   ..$ patent_title: chr "Anti-inflammatory analgesic patch"
#>   ..$ patent_year : int 1992
```

## FAQs

#### I'm sure my query is well formatted and correct but I keep getting an error. What's the deal?

The API query syntax guidelines do not cover all of the API's behavior. Specifically, there are several things that you cannot do which are not documented on the API's webpage. The [writing queries vignette](https://docs.ropensci.org/patentsview/articles/writing-queries.html) has more details on this. 

#### Does the API have any rate limiting/throttling controls?

Not at the moment.

#### How do I download more than 100,000 records?

Your best bet is to split your query into pieces based on dates, then concatenate the results together. For example, the below query would return more than 100,000 records for the patents endpoint:


```r
query <- with_qfuns(text_any(patent_abstract = 'tool animal'))
```

To download all of the records associated with this query, we could split it into two pieces and make two calls to `search_pv()`:


```r
query_1a <- with_qfuns(
  and(
    text_any(patent_abstract = 'tool animal'),
    lte(patent_date = "2010-01-01")
  )
)

query_1b <- with_qfuns(
  and(
    text_any(patent_abstract = 'tool animal'),
    gt(patent_date = "2010-01-01")
  )
)
```

#### How do I access the data frames inside the list columns returned by `search_pv()`?

Let’s consider the following data, in which assignees are the primary entity while applications and “government interest statements” are the secondary entities (also referred to as subentities):


```r
# Create field list
asgn_flds <- c("assignee_id", "assignee_organization")
subent_flds <- get_fields("assignees", c("applications", "gov_interests"))
fields <- c(asgn_flds, subent_flds)

# Pull data
res <- search_pv(
  query = qry_funs$contains(inventor_last_name = "smith"), 
  endpoint = "assignees", 
  fields = fields
)
res$data
#> #### A list with a single data frame (with list column(s) inside) on an assignee level:
#> 
#> List of 1
#>  $ assignees:'data.frame':	25 obs. of  4 variables:
#>   ..$ assignee_id          : chr [1:25] "e9449f65-8659-4611-a16d-18f65af5b3b6"..
#>   ..$ assignee_organization: chr [1:25] "U.S. Philips Corporation" ...
#>   ..$ applications         :List of 25
#>   ..$ gov_interests        :List of 25
```

`res$data` has vector columns for those fields that belong to the primary entity (e.g., `res$data$assignees$assignee_id`) and list columns for those fields that belong to any secondary entity (e.g., `res$data$assignees$applications`). You have two good ways to pull out the data frames that are nested inside these list columns: 

1. **Use tidyr::unnest.** (This is probably the easier choice of the two). 


```r
library(tidyr)
#> 
#> Attaching package: 'tidyr'
#> The following object is masked from 'package:magrittr':
#> 
#>     extract

# Get assignee/application data:
res$data$assignees %>% 
  unnest(applications) %>%
  head()
#> # A tibble: 6 x 8
#>   assignee_id   assignee_organi… app_country app_date app_number app_type app_id
#>   <chr>         <chr>            <chr>       <chr>    <chr>      <chr>    <chr> 
#> 1 e9449f65-865… U.S. Philips Co… US          1974-01… 05431439   05       05/43…
#> 2 e9449f65-865… U.S. Philips Co… US          1975-07… 05600148   05       05/60…
#> 3 e9449f65-865… U.S. Philips Co… US          1976-01… 05648308   05       05/64…
#> 4 e9449f65-865… U.S. Philips Co… US          1975-09… 05618031   05       05/61…
#> 5 e9449f65-865… U.S. Philips Co… US          1979-02… 06013951   06       06/01…
#> 6 e9449f65-865… U.S. Philips Co… US          1979-01… 06002418   06       06/00…
#> # … with 1 more variable: gov_interests <list>

# Get assignee/gov_interest data:
res$data$assignees %>% 
  unnest(gov_interests) %>%
  head()
#> # A tibble: 6 x 10
#>   assignee_id   assignee_organizati… applications govint_contract… govint_org_id
#>   <chr>         <chr>                <list>       <chr>            <chr>        
#> 1 e9449f65-865… U.S. Philips Corpor… <df [19 × 5… <NA>             <NA>         
#> 2 bbbe8bb0-7e4… XEROX CORPORATION    <df [510 × … <NA>             <NA>         
#> 3 bbbe8bb0-7e4… XEROX CORPORATION    <df [510 × … ECD-8721551      31           
#> 4 bbbe8bb0-7e4… XEROX CORPORATION    <df [510 × … 70NANBOH3033     44           
#> 5 bbbe8bb0-7e4… XEROX CORPORATION    <df [510 × … 70NANBOH3033     44           
#> 6 66fc4d3d-4a3… Commonwealth Scient… <df [1 × 5]> <NA>             <NA>         
#> # … with 5 more variables: govint_org_level_one <chr>,
#> #   govint_org_level_two <chr>, govint_org_level_three <lgl>,
#> #   govint_org_name <chr>, govint_raw_statement <chr>
```

2. **Use patentsview::unnest_pv_data.** `unnest_pv_data()` creates a series of data frames (one for each entity level) that are like tables in a relational database. You provide it with the data returned by `search_pv()` and a field that can act as a unique identifier for the primary entities:


```r
unnest_pv_data(data = res$data, pk = "assignee_id")
#> List of 3
#>  $ applications :'data.frame':	1951 obs. of  6 variables:
#>   ..$ assignee_id: chr [1:1951] "e9449f65-8659-4611-a16d-18f65af5b3b6" ...
#>   ..$ app_country: chr [1:1951] "US" ...
#>   ..$ app_date   : chr [1:1951] "1974-01-07" ...
#>   ..$ app_number : chr [1:1951] "05431439" ...
#>   ..$ app_type   : chr [1:1951] "05" ...
#>   ..$ app_id     : chr [1:1951] "05/431439" ...
#>  $ gov_interests:'data.frame':	91 obs. of  8 variables:
#>   ..$ assignee_id                 : chr [1:91] "e9449f65-8659-4611-a16d-18f65"..
#>   ..$ govint_contract_award_number: chr [1:91] NA ...
#>   ..$ govint_org_id               : chr [1:91] NA ...
#>   ..$ govint_org_level_one        : chr [1:91] NA ...
#>   ..$ govint_org_level_two        : chr [1:91] NA ...
#>   ..$ govint_org_level_three      : logi [1:91] NA ...
#>   ..$ govint_org_name             : chr [1:91] NA ...
#>   ..$ govint_raw_statement        : chr [1:91] NA ...
#>  $ assignees    :'data.frame':	25 obs. of  2 variables:
#>   ..$ assignee_id          : chr [1:25] "e9449f65-8659-4611-a16d-18f65af5b3b6"..
#>   ..$ assignee_organization: chr [1:25] "U.S. Philips Corporation" ...
```

Now we are left with a series of flat data frames instead of having a single data frame with other data frames nested inside of it. These flat data frames can be joined together as needed via the primary key (`assignee_id`).

[^1]: You can use `get_endpoints()` to list the endpoint names as the API expects them to appear (e.g., `assignees`, `cpc_subsections`, `inventors`, `locations`, `nber_subcategories`, `patents`, and `uspc_mainclasses`).
[^2]: This webpage includes some details that are not relevant to the `query` argument in `search_pv`, such as the field list and sort parameter.
[^3]: You can unnest the data frames that are stored in the list columns using `unnest_pv_data()`. See the FAQs for details.
---
title: "Citation networks"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Citation networks}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



The following is a brief foray into patent citation networks. The analysis is done on 3 patents that describe patent citation analysis (PCA) themselves. 

The first step is to download the relevant data from the PatentsView API. We can use the CPC code of  [Y10S707/933](https://worldwide.espacenet.com/classification#!/CPC=Y10S707/933) to identify the patents that relate to PCA. 


```r
library(patentsview)
library(dplyr)
library(visNetwork)
library(magrittr)
library(stringr)
library(knitr)

# Write a query to pull patents assigned to the CPC code of "Y10S707/933"
query <- qry_funs$begins(cpc_subgroup_id = "Y10S707/933")

# Create a list of fields to pull from the API
fields <- c(
  "patent_number", 
  "patent_title",
  "cited_patent_number", # Which patents do these patents cite?
  "citedby_patent_number" # Which patents cite them?
)

# Send a request to the API
res <- search_pv(query, fields = fields, all_pages = TRUE)

# Unnest the data found in the list columns
res_lst <- unnest_pv_data(res$data, pk = "patent_number")
res_lst
#> List of 3
#>  $ cited_patents  :'data.frame':	971 obs. of  2 variables:
#>   ..$ patent_number      : chr [1:971] "10095778" ...
#>   ..$ cited_patent_number: chr [1:971] "4991087" ...
#>  $ citedby_patents:'data.frame':	806 obs. of  2 variables:
#>   ..$ patent_number        : chr [1:806] "10095778" ...
#>   ..$ citedby_patent_number: chr [1:806] "10635705" ...
#>  $ patents        :'data.frame':	11 obs. of  2 variables:
#>   ..$ patent_number: chr [1:11] "10095778" ...
#>   ..$ patent_title : chr [1:11] "Method and system for probabilistically quan"..
```

There are only 11 PCA patents. These patents cite 971 patents and are cited by 806 patents. Let's visualize the citations among the PCA patents. We'll create our visualization using the `visNetwork` package, which requires us to create a data frame of nodes and a data frame of edges.


```r
pat_title <- function(title, number) {
  temp_title <- str_wrap(title)
  i <- gsub("\\n", "<br>", temp_title)
  paste0('<a href="https://patents.google.com/patent/US', number, '">', i, '</a>')
}

edges <-
  res_lst$cited_patents %>%
    semi_join(x = ., y = ., by = c("cited_patent_number" = "patent_number")) %>%
    set_colnames(c("from", "to"))

nodes <-
  res_lst$patents %>%
    mutate(
      id = patent_number,
      label = patent_number,
      title = pat_title(patent_title, patent_number)
    )

visNetwork(
  nodes = nodes, edges = edges, height = "400px", width = "100%",
  main = "Citations among patent citation analysis (PCA) patents"
) %>%
  visEdges(arrows = list(to = list(enabled = TRUE))) %>%
  visIgraphLayout()
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

It looks like several of the patents cite patent number 6,499,026, perhaps indicating that this patent contains technology that is foundational to the field. However, when we hover over the nodes we see that several of the patents have the same title. Clicking on the titles brings us to their full text on Google Patents, which confirms that many of these PCA patents belong to the same patent family.[^1]  Let's choose one of the patents in each family to act as the family's representative. This will reduce the size of the subsequent network, while hopefully retaining its overall structure.


```r
p3 <- c("7797336", "9075849", "6499026")
res_lst2 <- lapply(res_lst, function(x) x[x$patent_number %in% p3, ])
```

With only 3 patents, it will probably be possible to visualize how these patents' cited and citing patents are all related to one another. Let's create a list of these "relevant patents" (i.e., the 3 patents plus all of their cited and citing patents)[^2], and then get a list of all of their cited patents (i.e., the patents that they cite). This list of cited patents will allow us to measure how similar the relevant patents are to one another. 


```r
rel_pats <-
  res_lst2$cited_patents %>%
    rbind(setNames(res_lst2$citedby_patents, names(.))) %>% 
    select(-patent_number) %>%
    rename(patent_number = cited_patent_number) %>%
    bind_rows(data.frame(patent_number = p3)) %>% 
    distinct() %>%
    filter(!is.na(patent_number))

# Look up which patents the relevant patents cite
rel_pats_res <- search_pv(
  query = list(patent_number = rel_pats$patent_number),
  fields =  c("cited_patent_number", "patent_number", "patent_title"), 
  all_pages = TRUE, method = "POST"
)

rel_pats_lst <- unnest_pv_data(rel_pats_res$data, "patent_number")
```

Now we know which patents the 514 relevant patents cite. This allows us to measure the similarity between the 514 patents by seeing how many cited references they share in common (a method known as [bibliographic coupling](https://en.wikipedia.org/wiki/Bibliographic_coupling)).


```r
cited_pats <-
  rel_pats_lst$cited_patents %>%
    filter(!is.na(cited_patent_number))

full_network <- 
  cited_pats %>%
    do({
      .$ind <- 
        group_by(., patent_number) %>% 
          group_indices()
        group_by(., patent_number) %>%  
          mutate(sqrt_num_cited = sqrt(n()))
    }) %>%
    inner_join(x = ., y = ., by = "cited_patent_number") %>%
    filter(ind.x > ind.y) %>%
    group_by(patent_number.x, patent_number.y) %>% 
    mutate(cosine_sim = n() / (sqrt_num_cited.x * sqrt_num_cited.y)) %>% 
    ungroup() %>%
    select(matches("patent_number\\.|cosine_sim")) %>%
    distinct()

kable(head(full_network))
```



|patent_number.x |patent_number.y | cosine_sim|
|:---------------|:---------------|----------:|
|10095778        |10073890        |  0.0317500|
|10095778        |10055864        |  0.0240008|
|10108700        |10055462        |  0.8894992|
|10140198        |10055864        |  0.0264628|
|10140198        |10095778        |  0.0088918|
|10140673        |10095778        |  0.2330866|

`full_network` contains the similarity score (`cosine_sim`) for all patent pairs that share at least one cited reference in common. This means that it probably contains a lot of patent pairs that have only one or two cited references in common, and thus aren't all that similar. Let's try to identify a natural level of `cosine_sim` to filter on so that our subsequent network is not too hairy.


```r
hist(
  full_network$cosine_sim, 
  main = "Similarity scores between patents relevant to PCA",
  xlab = "Cosine similarity", ylab = "Number of patent pairs"
)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

There appears to be a smallish group of patent pairs that are very similar to one another (`cosine_sim` > 0.8), which makes it tempting to choose 0.8 as a cutoff point. However, patent pairs that have reference lists that are this similar to each other are probably just patents in the same patent family. Let's choose 0.1 as a cutoff point instead, as there doesn't appear to be too many pairs above this point.[^3] 


```r
edges <- 
  full_network %>%
    filter(cosine_sim >= .1) %>% 
    rename(from = patent_number.x, to = patent_number.y, value = cosine_sim) %>%
    mutate(title = paste("Cosine similarity =", as.character(round(value, 3))))

nodes <-
  rel_pats_lst$patents %>%
    rename(id = patent_number) %>%
    mutate(
      # the 3 patents of interest will be represented as blue nodes, all others
      # will be yellow
      color = ifelse(id %in% p3, "#97C2FC", "#DDCC77"), 
      label = id,
      title = pat_title(patent_title, id)
    )

visNetwork(
  nodes = nodes, edges = edges, height = "700px", width = "100%",
  main = "Network of patents relevant to PCA"
) %>%
  visEdges(color = list(color = "#343434")) %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1)) %>%
  visIgraphLayout()
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)


[^1]: A patent family is a group of related patents, usually all authored by the same inventor and relating to the same technology.
[^2]: Defining the network of patents relevant to PCA as those that cite or are cited by the 3 patents of interest is fairly restrictive (i.e., it doesn't adequately capture all of the patents related to PCA). There are likely patents out there that aren't cited by nor cite any of the 3, but are still relevant to PCA. One would need to measure the similarity between all the patents that are in the general area of PCA to get a more complete picture of the patents in this area. This is a much harder problem, though, and would require more analysis than can fit in a single vignette.
[^3]: This is still a pretty arbitrary choice. Take a look at algorithms like the [disparity filter](http://www.pnas.org/content/106/16/6483.full.pdf) for a more systematic way to filter edges.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-fields.R
\name{get_endpoints}
\alias{get_endpoints}
\title{Get endpoints}
\usage{
get_endpoints()
}
\value{
A character vector with the names of the 7 endpoints. Those endpoints are:

\itemize{
   \item assignees
   \item cpc_subsections
   \item inventors
   \item locations
   \item nber_subcategories
   \item patents
   \item uspc_mainclasses
 }
}
\description{
This function reminds the user what the 7 possible PatentsView API endpoints
are.
}
\examples{
get_endpoints()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query-dsl.R
\docType{data}
\name{qry_funs}
\alias{qry_funs}
\title{List of query functions}
\format{
An object of class \code{list} of length 14.
}
\usage{
qry_funs
}
\value{
An object of class \code{pv_query}. This is basically just a simple
  list with a print method attached to it.
}
\description{
A list of functions that make it easy to write PatentsView queries. See the
details section below for a list of the 14 functions, as well as the
\href{https://docs.ropensci.org/patentsview/articles/writing-queries.html}{writing
queries vignette} for further details.
}
\details{
\strong{1. Comparison operator functions} \cr

There are 6 comparison operator functions that work with fields of type
integer, float, date, or string:
\itemize{
   \item \code{eq} - Equal to
   \item \code{neq} - Not equal to
   \item \code{gt} - Greater than
   \item \code{gte} - Greater than or equal to
   \item \code{lt} - Less than
   \item \code{lte} - Less than or equal to
 }

There are 2 comparison operator functions that only work with fields of type
string:
\itemize{
   \item \code{begins} - The string begins with the value string
   \item \code{contains} - The string contains the value string
 }

There are 3 comparison operator functions that only work with fields of type
fulltext:
\itemize{
   \item \code{text_all} - The text contains all the words in the value
   string
   \item \code{text_any} - The text contains any of the words in the value
   string
   \item \code{text_phrase} - The text contains the exact phrase of the value
   string
 }

\strong{2. Array functions} \cr

There are 2 array functions:
\itemize{
   \item \code{and} - Both members of the array must be true
   \item \code{or} - Only one member of the array must be true
 }

\strong{3. Negation function} \cr

There is 1 negation function:
\itemize{
   \item \code{not} - The comparison is not true
 }
}
\examples{
qry_funs$eq(patent_date = "2001-01-01")

qry_funs$not(qry_funs$eq(patent_date = "2001-01-01"))

}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cast-pv-data.R
\name{cast_pv_data}
\alias{cast_pv_data}
\title{Cast PatentsView data}
\usage{
cast_pv_data(data)
}
\arguments{
\item{data}{The data returned by \code{\link{search_pv}}. This is the first
element of the three-element result object you got back from
\code{search_pv}. It should be a list of length 1, with one data frame
inside it. See examples.}
}
\value{
The same type of object that you passed into \code{cast_pv_data}.
}
\description{
This will cast the data fields returned by \code{\link{search_pv}} so that
they have their most appropriate data types (e.g., date, numeric, etc.).
}
\examples{
\dontrun{

fields <- c("patent_date", "patent_title", "patent_year")
res <- search_pv(query = "{\"patent_number\":\"5116621\"}", fields = fields)
cast_pv_data(data = res$data)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fieldsdf}
\alias{fieldsdf}
\title{Fields data frame}
\format{
A data frame with 992 rows and 7 variables:
\describe{
  \item{endpoint}{The endpoint that this field record is for}
  \item{field}{The name of the field}
  \item{data_type}{The field's data type (string, date, float, integer,
    fulltext)}
  \item{can_query}{An indicator for whether the field can be included in
    the user query for the given endpoint}
  \item{group}{The group the field belongs to}
  \item{common_name}{The field's common name}
  \item{description}{A description of the field}
}
}
\usage{
fieldsdf
}
\description{
A data frame containing the names of retrievable and queryable fields for
each of the 7 API endpoints. A yes/no flag (\code{can_query}) indicates
which fields can be included in the user's query. You can also find this
data on the API's online documentation for each endpoint as well (e.g.,
the \href{https://patentsview.org/apis/api-endpoints/patents}{patents
endpoint field list table})
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query-dsl.R
\name{with_qfuns}
\alias{with_qfuns}
\title{With qry_funs}
\usage{
with_qfuns(code, envir = parent.frame())
}
\arguments{
\item{code}{Code to evaluate. See example.}

\item{envir}{Where should R look for objects present in \code{code} that
aren't present in \code{\link{qry_funs}}.}
}
\value{
The result of \code{code} - i.e., your query.
}
\description{
This function evaluates whatever code you pass to it in the environment of
the \code{\link{qry_funs}} list. This allows you to cut down on typing when
writing your queries. If you want to cut down on typing even more, you can
try assigning the \code{\link{qry_funs}} list into your global environment
with: \code{list2env(qry_funs, envir = globalenv())}.
}
\examples{
# Without with_qfuns, we have to do:
qry_funs$and(
  qry_funs$gte(patent_date = "2007-01-01"),
  qry_funs$text_phrase(patent_abstract = c("computer program")),
  qry_funs$or(
    qry_funs$eq(inventor_last_name = "ihaka"),
    qry_funs$eq(inventor_first_name = "chris")
  )
)

#...With it, this becomes:
with_qfuns(
 and(
   gte(patent_date = "2007-01-01"),
   text_phrase(patent_abstract = c("computer program")),
   or(
     eq(inventor_last_name = "ihaka"),
     eq(inventor_first_name = "chris")
   )
 )
)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-fields.R
\name{get_fields}
\alias{get_fields}
\title{Get list of retrievable fields}
\usage{
get_fields(endpoint, groups = NULL)
}
\arguments{
\item{endpoint}{The API endpoint whose field list you want to get. See
\code{\link{get_endpoints}} for a list of the 7 endpoints.}

\item{groups}{A character vector giving the group(s) whose fields you want
returned. A value of \code{NULL} indicates that you want all of the
endpoint's fields (i.e., do not filter the field list based on group
membership). See the field tables located online to see which groups you
can specify for a given endpoint (e.g., the
\href{https://patentsview.org/apis/api-endpoints/patents}{patents
endpoint table}), or use the \code{fieldsdf} table
(e.g., \code{unique(fieldsdf[fieldsdf$endpoint == "patents", "group"])}).}
}
\value{
A character vector with field names.
}
\description{
This function returns a vector of fields that you can retrieve from a given
API endpoint (i.e., the fields you can pass to the \code{fields} argument in
\code{\link{search_pv}}). You can limit these fields to only cover certain
entity group(s) as well (which is recommended, given the large number of
possible fields for each endpoint).
}
\examples{
# Get all assignee-level fields for the patents endpoint:
fields <- get_fields(endpoint = "patents", groups = "assignees")

#...Then pass to search_pv:
\dontrun{

search_pv(
  query = '{"_gte":{"patent_date":"2007-01-04"}}',
  fields = fields
)
}
# Get all patent and assignee-level fields for the patents endpoint:
fields <- get_fields(endpoint = "patents", groups = c("assignees", "patents"))

\dontrun{
#...Then pass to search_pv:
search_pv(
  query = '{"_gte":{"patent_date":"2007-01-04"}}',
  fields = fields
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unnest-pv-data.R
\name{get_ok_pk}
\alias{get_ok_pk}
\title{Get OK primary key}
\usage{
get_ok_pk(endpoint)
}
\arguments{
\item{endpoint}{The endpoint which you would like to know a potential primary
key for.}
}
\value{
The name of a primary key (\code{pk}) that you could pass to
  \code{\link{unnest_pv_data}}.
}
\description{
This function suggests a value that you could use for the \code{pk} argument
in \code{\link{unnest_pv_data}}, based on the endpoint you searched.
It will return a potential unique identifier for a given entity (i.e., a
given endpoint). For example, it will return "patent_number" when
\code{endpoint = "patents"}.
}
\examples{
get_ok_pk(endpoint = "inventors") # Returns "inventor_id"
get_ok_pk(endpoint = "cpc_subsections") # Returns "cpc_subsection_id"

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search-pv.R
\name{search_pv}
\alias{search_pv}
\title{Search PatentsView}
\usage{
search_pv(
  query,
  fields = NULL,
  endpoint = "patents",
  subent_cnts = FALSE,
  mtchd_subent_only = TRUE,
  page = 1,
  per_page = 25,
  all_pages = FALSE,
  sort = NULL,
  method = "GET",
  error_browser = NULL,
  ...
)
}
\arguments{
\item{query}{The query that the API will use to filter records. \code{query}
 can come in any one of the following forms:
 \itemize{
   \item A character string with valid JSON. \cr
   E.g., \code{'{"_gte":{"patent_date":"2007-01-04"}}'}

   \item A list which will be converted to JSON by \code{search_pv}. \cr
   E.g., \code{list("_gte" = list("patent_date" = "2007-01-04"))}

   \item An object of class \code{pv_query}, which you create by calling one
   of the functions found in the \code{\link{qry_funs}} list...See the
   \href{https://docs.ropensci.org/patentsview/articles/writing-queries.html}{writing
   queries vignette} for details.\cr
   E.g., \code{qry_funs$gte(patent_date = "2007-01-04")}
 }}

\item{fields}{A character vector of the fields that you want returned to you.
A value of \code{NULL} indicates that the default fields should be
returned. Acceptable fields for a given endpoint can be found at the API's
online documentation (e.g., check out the field list for the
\href{https://patentsview.org/apis/api-endpoints/patents}{patents
endpoint}) or by viewing the \code{fieldsdf} data frame
(\code{View(fieldsdf)}). You can also use \code{\link{get_fields}} to list
out the fields available for a given endpoint.}

\item{endpoint}{The web service resource you wish to search. \code{endpoint}
must be one of the following: "patents", "inventors", "assignees",
"locations", "cpc_subsections", "uspc_mainclasses", or "nber_subcategories".}

\item{subent_cnts}{Do you want the total counts of unique subentities to be
returned? This is equivalent to the \code{include_subentity_total_counts}
parameter found \href{https://patentsview.org/apis/api-query-language}{here}.}

\item{mtchd_subent_only}{Do you want only the subentities that match your
query to be returned? A value of \code{TRUE} indicates that the subentity
has to meet your query's requirements in order for it to be returned, while
a value of \code{FALSE} indicates that all subentity data will be returned,
even those records that don't meet your query's requirements. This is
equivalent to the \code{matched_subentities_only} parameter found
\href{https://patentsview.org/apis/api-query-language}{here}.}

\item{page}{The page number of the results that should be returned.}

\item{per_page}{The number of records that should be returned per page. This
value can be as high as 10,000 (e.g., \code{per_page = 10000}).}

\item{all_pages}{Do you want to download all possible pages of output? If
\code{all_pages = TRUE}, the values of \code{page} and \code{per_page} are
ignored.}

\item{sort}{A named character vector where the name indicates the field to
sort by and the value indicates the direction of sorting (direction should
be either "asc" or "desc"). For example, \code{sort = c("patent_number" =
"asc")} or \cr\code{sort = c("patent_number" = "asc", "patent_date" =
"desc")}. \code{sort = NULL} (the default) means do not sort the results.
You must include any fields that you wish to sort by in \code{fields}.}

\item{method}{The HTTP method that you want to use to send the request.
Possible values include "GET" or "POST". Use the POST method when
your query is very long (say, over 2,000 characters in length).}

\item{error_browser}{Deprecated}

\item{...}{Arguments passed along to httr's \code{\link[httr]{GET}} or
\code{\link[httr]{POST}} function.}
}
\value{
A list with the following three elements:
 \describe{
   \item{data}{A list with one element - a named data frame containing the
   data returned by the server. Each row in the data frame corresponds to a
   single value for the primary entity. For example, if you search the
   assignees endpoint, then the data frame will be on the assignee-level,
   where each row corresponds to a single assignee. Fields that are not on
   the assignee-level would be returned in list columns.}

   \item{query_results}{Entity counts across all pages of output (not just
   the page returned to you). If you set \code{subent_cnts = TRUE}, you will
   be returned both the counts of the primary entities and the subentities.}

   \item{request}{Details of the HTTP request that was sent to the server.
   When you set \code{all_pages = TRUE}, you will only get a sample request.
   In other words, you will not be given multiple requests for the multiple
   calls that were made to the server (one for each page of results).}
 }
}
\description{
This function makes an HTTP request to the PatentsView API for data matching
the user's query.
}
\examples{

\dontrun{

search_pv(query = '{"_gt":{"patent_year":2010}}')

search_pv(
  query = qry_funs$gt(patent_year = 2010),
  fields = get_fields("patents", c("patents", "assignees"))
)

search_pv(
  query = qry_funs$gt(patent_year = 2010),
  method = "POST",
  fields = "patent_number",
  sort = c("patent_number" = "asc")
)

search_pv(
  query = qry_funs$eq(inventor_last_name = "crew"),
  all_pages = TRUE
)

search_pv(
  query = qry_funs$contains(inventor_last_name = "smith"),
  endpoint = "assignees"
)

search_pv(
  query = qry_funs$contains(inventor_last_name = "smith"),
  config = httr::timeout(40)
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unnest-pv-data.R
\name{unnest_pv_data}
\alias{unnest_pv_data}
\title{Unnest PatentsView data}
\usage{
unnest_pv_data(data, pk = get_ok_pk(names(data)))
}
\arguments{
\item{data}{The data returned by \code{\link{search_pv}}. This is the first
element of the three-element result object you got back from
\code{search_pv}. It should be a list of length 1, with one data frame
inside it. See examples.}

\item{pk}{The column/field name that will link the data frames together. This
should be the unique identifier for the primary entity. For example, if you
used the patents endpoint in your call to \code{search_pv}, you could
specify \code{pk = "patent_number"}. \strong{This identifier has to have
been included in your \code{fields} vector when you called
\code{search_pv}}. You can use \code{\link{get_ok_pk}} to suggest a
potential primary key for your data.}
}
\value{
A list with multiple data frames, one for each entity/subentity.
  Each data frame will have the \code{pk} column in it, so you can link the
  tables together as needed.
}
\description{
This function converts a single data frame that has subentity-level list
columns in it into multiple data frames, one for each entity/subentity.
The multiple data frames can be merged together using the primary key
variable specified by the user (see the
\href{https://r4ds.had.co.nz/relational-data.html}{relational data} chapter
in "R for Data Science" for an in-depth introduction to joining tabular data).
}
\examples{
\dontrun{

fields <- c("patent_number", "patent_title", "inventor_city", "inventor_country")
res <- search_pv(query = '{"_gte":{"patent_year":2015}}', fields = fields)
unnest_pv_data(data = res$data, pk = "patent_number")
}

}
