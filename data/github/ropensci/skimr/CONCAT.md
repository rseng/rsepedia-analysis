<!-- README.md is generated from README.Rmd. Please edit that file -->

skimr <a href='https://docs.ropensci.org/skimr/'>
=================================================

<img src='https://docs.ropensci.org/skimr/reference/figures/logo.png'
align="right" height="139" /></a>

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![R build
status](https://github.com/ropensci/skimr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/skimr/actions?workflow=R-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/skimr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/skimr/)
[![This is an ROpenSci Peer reviewed
package](https://badges.ropensci.org/175_status.svg)](https://github.com/ropensci/software-review/issues/175)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/skimr)](https://cran.r-project.org/package=skimr)
[![cran
checks](https://cranchecks.info/badges/summary/skimr)](https://cranchecks.info/pkgs/skimr)

`skimr` provides a frictionless approach to summary statistics which
conforms to the [principle of least
surprise](https://en.wikipedia.org/wiki/Principle_of_least_astonishment),
displaying summary statistics the user can skim quickly to understand
their data. It handles different data types and returns a `skim_df`
object which can be included in a pipeline or displayed nicely for the
human reader.

**Note: `skimr` version 2 has major changes when skimr is used
programmatically. Upgraders should review this document, the release
notes and vignettes carefully.**

Installation
------------

The current released version of `skimr` can be installed from CRAN. If
you wish to install the current build of the next release you can do so
using the following:

    # install.packages("devtools")
    devtools::install_github("ropensci/skimr")

The APIs for this branch should be considered reasonably stable but
still subject to change if an issue is discovered.

To install the version with the most recent changes that have not yet
been incorporated in the master branch (and may not be):

    devtools::install_github("ropensci/skimr", ref = "develop")

Do not rely on APIs from the develop branch, as they are likely to
change.

Skim statistics in the console
------------------------------

`skimr`:

-   Provides a larger set of statistics than `summary()`, including
    missing, complete, n, and sd.
-   reports each data types separately
-   handles dates, logicals, and a variety of other types
-   supports spark-bar and spark-line based on the [pillar
    package](https://github.com/r-lib/pillar).

### Separates variables by class:

    skim(chickwts)

    ## ── Data Summary ────────────────────────
    ##                            Values  
    ## Name                       chickwts
    ## Number of rows             71      
    ## Number of columns          2       
    ## _______________________            
    ## Column type frequency:             
    ##   factor                   1       
    ##   numeric                  1       
    ## ________________________           
    ## Group variables            None    
    ## 
    ## ── Variable type: factor ───────────────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate ordered n_unique top_counts                        
    ## 1 feed                  0             1 FALSE          6 soy: 14, cas: 12, lin: 12, sun: 12
    ## 
    ## ── Variable type: numeric ──────────────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ## 1 weight                0             1  261.  78.1   108  204.   258  324.   423 ▆▆▇▇▃

### Presentation is in a compact horizontal format:

    skim(iris)

    ## ── Data Summary ────────────────────────
    ##                            Values
    ## Name                       iris  
    ## Number of rows             150   
    ## Number of columns          5     
    ## _______________________          
    ## Column type frequency:           
    ##   factor                   1     
    ##   numeric                  4     
    ## ________________________         
    ## Group variables            None  
    ## 
    ## ── Variable type: factor ───────────────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate ordered n_unique top_counts               
    ## 1 Species               0             1 FALSE          3 set: 50, ver: 50, vir: 50
    ## 
    ## ── Variable type: numeric ──────────────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ## 1 Sepal.Length          0             1  5.84 0.828   4.3   5.1  5.8    6.4   7.9 ▆▇▇▅▂
    ## 2 Sepal.Width           0             1  3.06 0.436   2     2.8  3      3.3   4.4 ▁▆▇▂▁
    ## 3 Petal.Length          0             1  3.76 1.77    1     1.6  4.35   5.1   6.9 ▇▁▆▇▂
    ## 4 Petal.Width           0             1  1.20 0.762   0.1   0.3  1.3    1.8   2.5 ▇▁▇▅▃

### Built in support for strings, lists and other column classes

    skim(dplyr::starwars)

    ## ── Data Summary ────────────────────────
    ##                            Values         
    ## Name                       dplyr::starwars
    ## Number of rows             87             
    ## Number of columns          14             
    ## _______________________                   
    ## Column type frequency:                    
    ##   character                8              
    ##   list                     3              
    ##   numeric                  3              
    ## ________________________                  
    ## Group variables            None           
    ## 
    ## ── Variable type: character ────────────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate   min   max empty n_unique whitespace
    ## 1 name                  0         1         3    21     0       87          0
    ## 2 hair_color            5         0.943     4    13     0       12          0
    ## 3 skin_color            0         1         3    19     0       31          0
    ## 4 eye_color             0         1         3    13     0       15          0
    ## 5 sex                   4         0.954     4    14     0        4          0
    ## 6 gender                4         0.954     8     9     0        2          0
    ## 7 homeworld            10         0.885     4    14     0       48          0
    ## 8 species               4         0.954     3    14     0       37          0
    ## 
    ## ── Variable type: list ─────────────────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate n_unique min_length max_length
    ## 1 films                 0             1       24          1          7
    ## 2 vehicles              0             1       11          0          2
    ## 3 starships             0             1       17          0          5
    ## 
    ## ── Variable type: numeric ──────────────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ## 1 height                6         0.931 174.   34.8    66 167     180 191     264 ▁▁▇▅▁
    ## 2 mass                 28         0.678  97.3 169.     15  55.6    79  84.5  1358 ▇▁▁▁▁
    ## 3 birth_year           44         0.494  87.6 155.      8  35      52  72     896 ▇▁▁▁▁

### Has a useful summary function

    skim(iris) %>%
      summary()

    ## ── Data Summary ────────────────────────
    ##                            Values
    ## Name                       iris  
    ## Number of rows             150   
    ## Number of columns          5     
    ## _______________________          
    ## Column type frequency:           
    ##   factor                   1     
    ##   numeric                  4     
    ## ________________________         
    ## Group variables            None

### Individual columns can be selected using tidyverse-style selectors

    skim(iris, Sepal.Length, Petal.Length)

    ## ── Data Summary ────────────────────────
    ##                            Values
    ## Name                       iris  
    ## Number of rows             150   
    ## Number of columns          5     
    ## _______________________          
    ## Column type frequency:           
    ##   numeric                  2     
    ## ________________________         
    ## Group variables            None  
    ## 
    ## ── Variable type: numeric ──────────────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ## 1 Sepal.Length          0             1  5.84 0.828   4.3   5.1  5.8    6.4   7.9 ▆▇▇▅▂
    ## 2 Petal.Length          0             1  3.76 1.77    1     1.6  4.35   5.1   6.9 ▇▁▆▇▂

### Handles grouped data

`skim()` can handle data that has been grouped using
`dplyr::group_by()`.

    iris %>%
      dplyr::group_by(Species) %>%
      skim()

    ## ── Data Summary ────────────────────────
    ##                            Values    
    ## Name                       Piped data
    ## Number of rows             150       
    ## Number of columns          5         
    ## _______________________              
    ## Column type frequency:               
    ##   numeric                  4         
    ## ________________________             
    ## Group variables            Species   
    ## 
    ## ── Variable type: numeric ──────────────────────────────────────────────────────────────────────────
    ##    skim_variable Species    n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ##  1 Sepal.Length  setosa             0             1 5.01  0.352   4.3  4.8   5     5.2    5.8 ▃▃▇▅▁
    ##  2 Sepal.Length  versicolor         0             1 5.94  0.516   4.9  5.6   5.9   6.3    7   ▂▇▆▃▃
    ##  3 Sepal.Length  virginica          0             1 6.59  0.636   4.9  6.22  6.5   6.9    7.9 ▁▃▇▃▂
    ##  4 Sepal.Width   setosa             0             1 3.43  0.379   2.3  3.2   3.4   3.68   4.4 ▁▃▇▅▂
    ##  5 Sepal.Width   versicolor         0             1 2.77  0.314   2    2.52  2.8   3      3.4 ▁▅▆▇▂
    ##  6 Sepal.Width   virginica          0             1 2.97  0.322   2.2  2.8   3     3.18   3.8 ▂▆▇▅▁
    ##  7 Petal.Length  setosa             0             1 1.46  0.174   1    1.4   1.5   1.58   1.9 ▁▃▇▃▁
    ##  8 Petal.Length  versicolor         0             1 4.26  0.470   3    4     4.35  4.6    5.1 ▂▂▇▇▆
    ##  9 Petal.Length  virginica          0             1 5.55  0.552   4.5  5.1   5.55  5.88   6.9 ▃▇▇▃▂
    ## 10 Petal.Width   setosa             0             1 0.246 0.105   0.1  0.2   0.2   0.3    0.6 ▇▂▂▁▁
    ## 11 Petal.Width   versicolor         0             1 1.33  0.198   1    1.2   1.3   1.5    1.8 ▅▇▃▆▁
    ## 12 Petal.Width   virginica          0             1 2.03  0.275   1.4  1.8   2     2.3    2.5 ▂▇▆▅▇

### Behaves nicely in pipelines

    iris %>%
      skim() %>%
      dplyr::filter(numeric.sd > 1)

    ## ── Data Summary ────────────────────────
    ##                            Values    
    ## Name                       Piped data
    ## Number of rows             150       
    ## Number of columns          5         
    ## _______________________              
    ## Column type frequency:               
    ##   numeric                  1         
    ## ________________________             
    ## Group variables            None      
    ## 
    ## ── Variable type: numeric ──────────────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ## 1 Petal.Length          0             1  3.76  1.77     1   1.6  4.35   5.1   6.9 ▇▁▆▇▂

Knitted results
---------------

Simply skimming a data frame will produce the horizontal print layout
shown above. We provide a `knit_print` method for the types of objects
in this package so that similar results are produced in documents. To
use this, make sure the `skimmed` object is the last item in your code
chunk.

    faithful %>%
      skim()

<table>
<caption>Data summary</caption>
<tbody>
<tr class="odd">
<td style="text-align: left;">Name</td>
<td style="text-align: left;">Piped data</td>
</tr>
<tr class="even">
<td style="text-align: left;">Number of rows</td>
<td style="text-align: left;">272</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Number of columns</td>
<td style="text-align: left;">2</td>
</tr>
<tr class="even">
<td style="text-align: left;">_______________________</td>
<td style="text-align: left;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">Column type frequency:</td>
<td style="text-align: left;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">numeric</td>
<td style="text-align: left;">2</td>
</tr>
<tr class="odd">
<td style="text-align: left;">________________________</td>
<td style="text-align: left;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">Group variables</td>
<td style="text-align: left;">None</td>
</tr>
</tbody>
</table>

**Variable type: numeric**

<table>
<thead>
<tr class="header">
<th style="text-align: left;">skim_variable</th>
<th style="text-align: right;">n_missing</th>
<th style="text-align: right;">complete_rate</th>
<th style="text-align: right;">mean</th>
<th style="text-align: right;">sd</th>
<th style="text-align: right;">p0</th>
<th style="text-align: right;">p25</th>
<th style="text-align: right;">p50</th>
<th style="text-align: right;">p75</th>
<th style="text-align: right;">p100</th>
<th style="text-align: left;">hist</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">eruptions</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">3.49</td>
<td style="text-align: right;">1.14</td>
<td style="text-align: right;">1.6</td>
<td style="text-align: right;">2.16</td>
<td style="text-align: right;">4</td>
<td style="text-align: right;">4.45</td>
<td style="text-align: right;">5.1</td>
<td style="text-align: left;">▇▂▂▇▇</td>
</tr>
<tr class="even">
<td style="text-align: left;">waiting</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">70.90</td>
<td style="text-align: right;">13.59</td>
<td style="text-align: right;">43.0</td>
<td style="text-align: right;">58.00</td>
<td style="text-align: right;">76</td>
<td style="text-align: right;">82.00</td>
<td style="text-align: right;">96.0</td>
<td style="text-align: left;">▃▃▂▇▂</td>
</tr>
</tbody>
</table>

Customizing skimr
-----------------

Although skimr provides opinionated defaults, it is highly customizable.
Users can specify their own statistics, change the formatting of
results, create statistics for new classes and develop skimmers for data
structures that are not data frames.

### Specify your own statistics and classes

Users can specify their own statistics using a list combined with the
`skim_with()` function factory. `skim_with()` returns a new `skim`
function that can be called on your data. You can use this factory to
produce summaries for any type of column within your data.

Assignment within a call to `skim_with()` relies on a helper function,
`sfl` or `skimr` function list. By default, functions in the `sfl` call
are appended to the default skimmers, and names are automatically
generated as well.

    my_skim <- skim_with(numeric = sfl(mad))
    my_skim(iris, Sepal.Length)

But you can also helpers from the `tidyverse` to create new anonymous
functions that set particular function arguments. The behavior is the
same as in `purrr` or `dplyr`, with both `.` and `.x` as acceptable
pronouns. Setting the `append = FALSE` argument uses only those
functions that you’ve provided.

    my_skim <- skim_with(
      numeric = sfl(
        iqr = IQR,
        p01 = ~ quantile(.x, probs = .01)
        p99 = ~ quantile(., probs = .99)
      ),
      append = FALSE
    )
    my_skim(iris, Sepal.Length)

And you can remove default skimmers by setting them to `NULL`.

    my_skim <- skim_with(numeric = sfl(hist = NULL))
    my_skim(iris, Sepal.Length)

### Skimming other objects

`skimr` has summary functions for the following types of data by
default:

-   `numeric` (which includes both `double` and `integer`)
-   `character`
-   `factor`
-   `logical`
-   `complex`
-   `Date`
-   `POSIXct`
-   `ts`
-   `AsIs`

`skimr` also provides a small API for writing packages that provide
their own default summary functions for data types not covered above. It
relies on R S3 methods for the `get_skimmers` function. This function
should return a `sfl`, similar to customization within `skim_with()`,
but you should also provide a value for the `class` argument. Here’s an
example.

    get_skimmers.my_data_type <- function(column) {
      sfl(
        .class = "my_data_type",
        p99 = quantile(., probs = .99)
      )
    }

Limitations of current version
------------------------------

We are aware that there are issues with rendering the inline histograms
and line charts in various contexts, some of which are described below.

### Support for spark histograms

There are known issues with printing the spark-histogram characters when
printing a data frame. For example, `"▂▅▇"` is printed as
`"<U+2582><U+2585><U+2587>"`. This longstanding problem [originates in
the low-level
code](https://r.789695.n4.nabble.com/Unicode-display-problem-with-data-frames-under-Windows-td4707639.html)
for printing dataframes. While some cases have been addressed, there
are, for example, reports of this issue in Emacs ESS.

This means that while `skimr` can render the histograms to the console
and in RMarkdown documents, it cannot in other circumstances. This
includes:

-   converting a `skimr` data frame to a vanilla R data frame, but
    tibbles render correctly
-   in the context of rendering to a pdf using an engine that does not
    support utf-8.

One workaround for showing these characters in Windows is to set the
CTYPE part of your locale to Chinese/Japanese/Korean with
`Sys.setlocale("LC_CTYPE", "Chinese")`. The helper function
`fix_windows_histograms()` does this for you.

And last but not least, we provide `skim_without_charts()` as a
fallback. This makes it easy to still get summaries of your data, even
if unicode issues continue.

### Printing spark histograms and line graphs in knitted documents

Spark-bar and spark-line work in the console, but may not work when you
knit them to a specific document format. The same session that produces
a correctly rendered HTML document may produce an incorrectly rendered
PDF, for example. This issue can generally be addressed by changing
fonts to one with good building block (for histograms) and Braille
support (for line graphs). For example, the open font “DejaVu Sans” from
the `extrafont` package supports these. You may also want to try
wrapping your results in `knitr::kable()`. Please see the vignette on
using fonts for details.

Displays in documents of different types will vary. For example, one
user found that the font “Yu Gothic UI Semilight” produced consistent
results for Microsoft Word and Libre Office Write.

### Stripping metadata and empty results tables

In POSIX systems, `skimr` tries to remove the tibble metadata when
producing the results. A complicating factor is tibble’s color support,
which depends on environment settings. In particular, not all Windows
terminals support colors in the way that tibble expects.

So, by default, we disable removing metadata on windows. You can turn
this feature on with an option. Either set it when calling print or
globally.

    skimmed <- skim(chickwts)
    print(skimmed, strip_metadata = TRUE)
    options(skimr_strip_metadata = TRUE)

Separately, you might need to check the option `crayon.enabled`.
Similarly, if your skimr results tables are empty you may need to run
the following

    options(crayon.enabled = FALSE)

You need to do this one time per session.

Contributing
------------

We welcome issue reports and pull requests, including potentially adding
support for commonly used variable classes. However, in general, we
encourage users to take advantage of skimr’s flexibility to add their
own customized classes. Please see the
[contributing](https://docs.ropensci.org/skimr/CONTRIBUTING.html) and
[conduct](https://ropensci.org/code-of-conduct/) documents.

[![ropenci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# skimr 2.2.0

### NEW FEATURES

*   skim() used within a function now prints the data frame name.

# skimr 2.1.3

### MINOR IMPROVEMENTS

*   Add support for data tables when dtplyr is used.
*   Improve tests.

# skimr 2.1.2

### MINOR IMPROVEMENTS

*   Add support for lubridate Timespan objects.
*   Improvements to Supporting Additional Objects vignette.

### BUG FIXES

*   Update package to work with new version of `knitr`.

# skimr 2.1.1 (2020-04-15)

### MINOR IMPROVEMENTS

*   Prepare for release of dplyr 1.0 and related packages.
*   0-length sfls are now permitted.

# skimr 2.1.0 (2020-01-10)

### NEW FEATURES

We've made `to_long()` generic, supporting a more intuitive interface.

*   Called on a `skim_df`, it reshapes the output into the V1 long style.
*   Called on other tibble-like objects, it first skims then produces the long
    output. You can pass a custom skim function, like `skim_tee()`

Thanks @sethlatimer for suggesting this feature.

### BUG FIXES

*   Update package to work with new version of `tibble`.
*   Adds more flexibility in the rule width for `skimr::summarize()`.
*   More README badges and documentation crosslinks

# skimr 2.0.1 (2019-11-23)

### BUG FIXES

Address failed build in CRAN due to lack of UTF-8 support in some platforms.

# skimr 2.0.0 (2019-11-12)

### Welcome to skimr V2

V2 is a complete rewrite of `skimr`, incorporating all of the great feedback the
developers have received over the last year. A big thank you goes to @GShotwell,
@akraemer007, @puterleat, @tonyfischetti, @Nowosad, @rgayler, @jrosen48,
@randomgambit, @elben10, @koliii, @AndreaPi, @rubenarslan, @GegznaV, @svraka,
@dpprdan and to our ROpenSci reviewers @jenniferthompson and @jimhester for all
of the great support and feedback over the last year. We couldn't have done this
without you.

For most users using `skimr` will not change in terms of visual outputs. However
for users who use `skimr` outputs as part of a larger workflow the differences
are substantial.

### Breaking changes

#### The `skim_df`

We've changed the way data is represented within `skimr` to closer match
expectations. It is now wide by default. This makes piping statistics much
simpler

```
skim(iris) %>%
  dplyr::filter(numeric.sd > 1)
```

This means that the old reshaping functions `skim_to_wide()` and
`skim_to_list()` are deprecated. The latter is replaced with a reshaping
function called `partition()` that breaks a `skim_df` into a list by data type.
Similarly, `yank()` gets a specific data type from the `skim_df`. `to_long()`
gets you data that is closest to the format in the old API.

As the above example suggests, columns of summary statistics are prefixed by
`skim_type`. That is, statistics from numeric columns all begin `numeric.`,
those for factors all begin `factor.`, and so on.

#### Rendering

We've deprecated support for `pander()` and our `kable()` method. Instead, we
now support `knitr` through the `knit_print()` API. This is much more seamless
than before. Having a `skim_df` as the final object in a code chunk should
produce nice results in the majority of RMarkdown formats.

#### Customizing and extending

We've deprecated the previous approach customization. We no longer use
`skim_format()` and `skim_with()` no longer depends on a global state. Instead
`skim_with()` is now a function factory. Customization creates a new skimming
function.

```
my_skim <- skim_with(numeric = sfl(mad = mad))
```

The fundamental tool for customization is the `sfl` object, a skimmer function
list. It is used within `skim_with()` and also within our new API for adding
default functions for new data types, the generic `get_skimmers()`.

Most of the options set in `skim_format` are now either in function arguments or
print arguments. The former can be updated using `skim_with`, the latter in a
call to `print()`. In RMarkdown documents, you can change the number of
displayed digits by adding the `skimr_digits` option to your code chunk.

### OTHER NEW FEATURES

*   Substantial improvements to `summary()`, and it is now incorporated into
    `print()` methods.
*   `focus()` is like `dplyr::select()`, but it keeps around the columns
    `skim_type` and `skim_variable`.
*   We are also evaluating the behavior of different `dplyr` verbs to make sure
    that they place nice with `skimr` objects.
*   While `skimr` has never really focused on performance, it should do a better
    job on big data sets with lots of different columns.
*   New statistic for character variables counting the number of rows that are
    completely made up of white space.
*   We now export `skim_without_charts()` as a fallback for when unicode support
    is not possible.
*   By default, `skimr` removes the tibble metadata when generating output. On
    some platforms, this can lead to all output getting removed. To disable that
    behavior, set either `strip_metadata = FALSE` when calling print or use
    `options(skimr_strip_metadata = FALSE)`.

### BUG FIXES

*   Adjust code for several tidyverse soft deprecations.
*   Fix issue where multibyte characters were causing an error.

### MINOR IMPROVEMENTS

*   Change top_counts to use useNA = "no".

# skimr 1.0.6 (2019-05-27)

### BUG FIXES

*   Fix issue where skim_tee() was not respecting ... options.
*   Fix issue where all NA character vectors were not returning NA for max() and
    min()

# skimr 1.0.5 (2019-01-05)

This is likely to be the last release of skimr version 1. Version 2 has major
changes to the API. Users should review and prepare for those changes now.

### BUG FIXES

*   Fix issue where multibyte characters were causing an error.
*   Fix problem in which purrr cannot find mean.default.

# skimr 1.0.4 (2018-01-12)

This is likely to be the last release of skimr version 1. Version 2 has major
changes to the API. Users should review and prepare for those changes now.

### BUG FIXES

*   Fix failures in handling dplyr verbs related to upcoming release of dplyr
    0.8.0.

# skimr 1.0.3 (2018-06-06)

### NEW FEATURES

*   You can use skim_with() with a nest list of functions: `skim_with(.list =
    mylist)` or `skim_with(!!!mylist)`
*   More polished display of subtables in default printing.

### BUG FIXES

*   Fix issue with conflict between knitr and skimr versions of kable() that
    occurred intermittently.
*   Do not skim a class when the skimmer list is empty for that class.
*   Fix a mistake in a test of skim_print for top counts.

# skimr 1.0.2 (2018-04-04)

### NEW FEATURES

*   You can create skimmers with the formula syntax from `rlang`:
    `skim_with(iqr = ~IQR(.x, na.rm = TRUE))`.

### MAJOR CHANGES

*   The median label has been changed to p50 for consistency with the previous
    changes to p0 and p100.

### MINOR IMPROVEMENTS

*   Improvements and corrections to to README and other documentation.
*   New vignette showing defaults for skimmers and formats.
*   Vector output match data frame output more closely.
*   Add minimum required version for testhat.
*   Add minimum required version for knitr.

### BUG FIXES

*   You can use `skim_with()` to add and remove skimmers at the same time, i.e.
    `skim_with(iqr = IQR, hist = NULL)` works as expected.
*   Histograms work when Inf or -Inf are present.
*   Change seq( ) parameter to length.out to avoid problems with name matching.
*   Summary should not display a data frame name of "." (which occurs when
    piping begins with the data frame).

# skimr 1.0.1 (2018-01-09)

### NEW FEATURES

*   Add support for spark plots on Windows

### MAJOR CHANGES

*   `spark_line()` and `spark_bar()` are no longer exported
*   Default statistics for numeric changed from `min(x)` and `max(x)` to
    `quantile(x, probs = 0)` and `quantile(x, probs = 1)`. These changes lead to
    more predictable behaviors when a column is all NA values.

#### MINOR IMPROVEMENTS

*   Add minimimum required version for stringr
*   Improve documentation in general, especially those related to fonts

### BUG FIXES

*   Fix issue where a histogram for data with all `NA`s threw an error
*   Suppress progress bars from `dplyr::do()`

# skimr 0.92 (2017-12-19)

### MAJOR CHANGES

*   `skim_v()` is no longer exported. Vectors are now directly supported via
    `skim.default()`.
*   Change license to GPL 3

### NEW FEATURES

*   Add support for `kable()` and `pander()` for `skim_df` objects.
*   Add summary method for `skim_df` objects.
*   Add support for tidy select to skim specific columns of a data frame.
*   Add support for skimming individual vectors via `skim.default()`.

# skimr 0.91 (2017-10-14)

### NEW FEATURES

*   Handling of grouped data (generated by `dplyr::group_by()`)
*   Printing for all column classes
*   Add indicator of if a factor is ordered to skim object for factor
*   Introduction of flexible formatting
*   Easy dropping of individual functions
*   Vignettes for basic use and use with specialized object types
*   Updated README and added CONTRIBUTING.md and CONDUCT.md
*   New public get_skimmers function to access skim functions
*   Support for difftime class

### MINOR IMPROVEMENTS

*   Add header to print providing summary information about data.

### BUG FIXES

*   Change from Colformat to Pillar.

# skimr 0.900 (2017-07-16)

### BUG FIXES

*   Fix documentation for get_fun_names()
*   Fix test and build errors and notes
# Contributing

Contributions to `skimr` whether in the form of bug fixes, issue reports, new
code or documentation improvement are welcome. Please use the github issue
tracker. For any pull request please link to or open a corresponding issue in
the issue tracker. Please ensure that you have notifications turned on and
respond to questions, comments or needed changes promptly.

## Understanding the scope of skimr

`skimr` solves a very specific set of problems focused on the compact, flexible
and useful display of summary data in the console. By itself it is not intended
as a replacement for packages that create publication ready tables. The basic
concept is that of "skimming" a data frame or tibble to get an overview of the
data it contains.

One intended group of users is students in a first semester statistics class. As
such, the package is focused on data types that are widely used. One general
guideline is that if a data type is not found in the `datasets` package it will
not be directly supported in `skimr`. Fortunately, `skim()` has a generic
internal function for handling a variety of data types `get_skimmers()`. See the
documentation for that function or the vignette "Supporting additional objects"
for documentation on how to do this.

Similarly, `skimr` is deeply tied to the `tidyverse` and `dplyr` in particular.
The comes with a lot of benefits, but some constraints too. Most importantly,
data processed by `skim()` needs to be an object that inherits from a data frame
or in a form that can be coerced to a data frame.

## Tests

`skimr` uses `testthat` for testing. Please try to provide 100% test coverage
for any submitted code and always check that existing tests continue to pass. If
you are a beginner and need help with writing a test, mention this in the issue
and we will try to help.

## Pull requests

Pull requests should be against the _develop_ branch not the master branch. You
can set this when creating your pull request. Please make a separately named
branch to submit. Keep each branch for a complete specific issue. If you create
a pull request by editing in the GitHub web editor and you end up with multiple
pull requests, note that in your issue comments.

## Code style

We follow the [tidyverse style guide](http://style.tidyverse.org/).

## Pre commits

To enforce coding style and support development, we rely on [pre-commit.com],
and the [R precommit package](https://github.com/lorenzwalthert/precommit). This
tool runs a series of additional checks for your code before `git commit`
completes.

To install the package and enable precommits, run the following:

```
# once on your system
remotes::install_github("lorenzwalthert/precommit")
precommit::install_precommit()

# once in every git repo either
# * after cloning a repo that already uses pre-commit or
# * if you want introduce pre-commit to this repo
precommit::use_precommit()
```

The checks will run automatically from there.

## Code of Conduct

When contributing to `skimr` you must follow the [code of conduct defined by rOpenSci](https://ropensci.org/code-of-conduct/).
[Theme song: *PSA* by
Jay-Z](https://www.youtube.com/watch?v=-LzdKH1naok)

We announced the testing version of `skimr` v2 on [June 19,
2018](https://github.com/ropensci/skimr/issues/341). After more than a
year of (admittedly intermittent) work, we’re thrilled to be able to say
that the package is ready to go to CRAN. So, what happened over the last
year? And why are we so excited for V2?

Setting the stage
-----------------

Before we can talk about the last year of `skimr` development, we need
to lay out the timeline that got us to this point. For those deeply
enmeshed in `skimr` lore, all
[dozens](https://imgur.com/gallery/R1fdEt3) of you, bear with.

`skimr` was originally an [rOpenSci
unconf17](https://ropensci.org/blog/2017/07/11/skimr/) project, a big
collaboration between eight different participants that resulted in a
conceptual outline of the package and a basic working version.
Participating in the unconf was a truly magical experience, with
everyone bringing a tremendous amount of energy and ideas to the
project, and implementation happening over a flurry of [“fancy git
commits”](https://twitter.com/AmeliaMN/status/867818976666976256).

About six months later, we released our first version on CRAN. The time
between these two milestones was mostly spent on fleshing out all of the
different ideas that were generated during the unconf (like handling
grouped data frames) and fixing all the bugs we discovered along the
way.

Getting the package on CRAN opened the gates for bug reports and feature
requests on [GitHub](https://github.com/ropensci/skimr/issues). About
the same time we pushed our first version to CRAN, Elin got `skimr`’s
rOpenSci’s package [peer
review](https://github.com/ropensci/software-review/issues/175) started
(thank you Jennifer and Jim!), opening another incredibly useful channel
for collecting feedback on the package. All of these new ideas and
suggestions gave us the opportunity to really push `skimr` to the next
level, but doing that would require rethinking the package, from the
ground up.

A month after finishing the peer review (and six months after the
process began), we announced v2. Over the first phase of `skimr`’s life,
we accumulated 700 commits, two release, 400 GitHub stars, 95 percent
code coverage and a lifetime’s worth of [unicode rendering
bugs](https://github.com/ropensci/skimr#support-for-spark-histograms)!

Just kidding! We love our little histograms, even when they don’t love
us back! For those of you that might have never seen `skimr`, using the
package typically boils down to a single function call:

    library(skimr)
    library(dplyr)
    options(width = 90)

    skim(iris)

    ## ── Data Summary ────────────────────────
    ##                            Values
    ## Name                       iris  
    ## Number of rows             150   
    ## Number of columns          5     
    ## _______________________          
    ## Column type frequency:           
    ##   factor                   1     
    ##   numeric                  4     
    ## ________________________         
    ## Group variables            None  
    ## 
    ## ── Variable type: factor ─────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate ordered n_unique top_counts               
    ## 1 Species               0             1 FALSE          3 set: 50, ver: 50, vir: 50
    ## 
    ## ── Variable type: numeric ────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ## 1 Sepal.Length          0             1  5.84 0.828   4.3   5.1  5.8    6.4   7.9 ▆▇▇▅▂
    ## 2 Sepal.Width           0             1  3.06 0.436   2     2.8  3      3.3   4.4 ▁▆▇▂▁
    ## 3 Petal.Length          0             1  3.76 1.77    1     1.6  4.35   5.1   6.9 ▇▁▆▇▂
    ## 4 Petal.Width           0             1  1.20 0.762   0.1   0.3  1.3    1.8   2.5 ▇▁▇▅▃

Getting it right
----------------

Under normal circumstances (i.e. not during a hackathon), most software
engineering projects begin with a design phase and series of
increasingly detailed design docs. `skimr` is only a few hundred lines
of code, which means “increasingly detailed design docs” translates to
one doc. But we did actually write it! [It’s
here](https://docs.google.com/document/d/18lBStDZzd1rJq08O-4Sw2qHhuHEZ79QX4sBkeyzWNFY/edit#heading=h.5x0d5h95i329).
And it still goes a good job of laying out some of the big ideas we were
interested in taking on for v2.

-   Eliminating frictions that resulted from differences in the way we
    stored data vs how it was displayed to users
-   Getting away from using a global environment to configure `skimr`
-   Making it easier for others to extend `skimr`
-   Create more useful ways to use `skimr`

Better internal data structures
-------------------------------

In v1, `skimr` stored all of its data in a “long format”, data frame.
Although hidden from the user by its print methods, this format would
appear any time you’d try do something with the results of a `skim()`
call. It looked something like this:

    skim(mtcars) %>% dplyr::filter(stat=="hist")

    # A tibble: 11 x 6
       variable type    stat  level value formatted
       <chr>    <chr>   <chr> <chr> <dbl> <chr>
     1 mpg      numeric hist  .all     NA ▃▇▇▇▃▂▂▂
     2 cyl      numeric hist  .all     NA ▆▁▁▃▁▁▁▇
     3 disp     numeric hist  .all     NA ▇▆▁▂▅▃▁▂
     4 hp       numeric hist  .all     NA ▃▇▃▅▂▃▁▁
     5 drat     numeric hist  .all     NA ▃▇▁▅▇▂▁▁
     6 wt       numeric hist  .all     NA ▃▃▃▇▆▁▁▂
     7 qsec     numeric hist  .all     NA ▃▂▇▆▃▃▁▁
     8 vs       numeric hist  .all     NA ▇▁▁▁▁▁▁▆
     9 am       numeric hist  .all     NA ▇▁▁▁▁▁▁▆
    10 gear     numeric hist  .all     NA ▇▁▁▆▁▁▁▂
    11 carb     numeric hist  .all     NA ▆▇▂▇▁▁▁▁

Big ups to anyone who looked at the rendered output and saw that this
was how you actually filtered the results. Hopefully there are even
better applications of your near-telepathic abilities.

Now, working with `skimr` is a bit more sane.

    skimmed <- iris %>%
      skim() %>%
      dplyr::filter(numeric.sd > 1)

    skimmed

    ## ── Data Summary ────────────────────────
    ##                            Values    
    ## Name                       Piped data
    ## Number of rows             150       
    ## Number of columns          5         
    ## _______________________              
    ## Column type frequency:               
    ##   numeric                  1         
    ## ________________________             
    ## Group variables            None      
    ## 
    ## ── Variable type: numeric ────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ## 1 Petal.Length          0             1  3.76  1.77     1   1.6  4.35   5.1   6.9 ▇▁▆▇▂

And

    dplyr::glimpse(skimmed)

    ## Observations: 1
    ## Variables: 15
    ## $ skim_type         <chr> "numeric"
    ## $ skim_variable     <chr> "Petal.Length"
    ## $ n_missing         <int> 0
    ## $ complete_rate     <dbl> 1
    ## $ factor.ordered    <lgl> NA
    ## $ factor.n_unique   <int> NA
    ## $ factor.top_counts <chr> NA
    ## $ numeric.mean      <dbl> 3.758
    ## $ numeric.sd        <dbl> 1.765298
    ## $ numeric.p0        <dbl> 1
    ## $ numeric.p25       <dbl> 1.6
    ## $ numeric.p50       <dbl> 4.35
    ## $ numeric.p75       <dbl> 5.1
    ## $ numeric.p100      <dbl> 6.9
    ## $ numeric.hist      <chr> "▇▁▆▇▂"

It’s still not perfect, as you need to rely on a *pseudo-namespace* to
refer to the column that you want. But this is unfortunately a necessary
trade-off. As the Rstats Bible, errr Hadley Wickham’s *Advanced R*,
states, all elements of [an atomic vector must have the same
type](https://adv-r.hadley.nz/vectors-chap.html). This normally isn’t
something that you have to think too much about, that is until you try
to combine the means of all your `Date` columns with the means of your
`numeric` columns and everything comes out utterly garbled. So instead
of that basket of laughs, we prefix columns names by their data type.

There’s a couple of other nuances here:

-   The data frame `skim()` produces always starts off with some
    metadata columns
-   Functions that always produce the same, regardless of input type,
    can be treated as `base_skimmers` and don’t need a namespace

### Manipulating internal data

A better representation of internal data comes with better tools for
reshaping the data and getting it for other contexts. A common request
in v1 was tooling to handle the `skimr` subtables separately. We now do
this with `partition()`. It replaces the v1 function `skim_to_list()`.

    partition(skimmed)

    ## $numeric
    ## 
    ## ── Variable type: numeric ────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ## 1 Petal.Length          0             1  3.76  1.77     1   1.6  4.35   5.1   6.9 ▇▁▆▇▂

You can undo a call to `partition()` with `bind()`, which joins the
subtables into the original `skim_df` object and properly accounts for
metadata. You can skip a step with the function `yank()`, which calls
partition and pulls out a particular subtable

    yank(skimmed, "numeric")

    ## 
    ## ── Variable type: numeric ────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p25   p50   p75  p100 hist 
    ## 1 Petal.Length          0             1  3.76  1.77     1   1.6  4.35   5.1   6.9 ▇▁▆▇▂

Last, with support something close to the older format with the
`to_long()` function. This can be added for something close to backwards
compatibility. Being realistic on open source sustainability means that
we are not able to support 100% backward compatibility in v2 even with
new functions. Meanwhile you can keep using v1 if you are happy with it.
However, because `skimr`’s dependencies are under ongoing development,
sooner or later skimr v1 will no longer work with updates to them.

### Working with dplyr

Using `skimr` in a `dplyr` pipeline was part of the original package
design, and we’ve needed to devote some extra love to making sure that
everything is as seamless as possible. Part of this is due to the object
produce by `skim()`, which we call `skim_df`. It’s a little weird in
that it needs both metadata and columns in the underlying data frame.

In practice, this means that you can coerce it into a different type
through normal `dplyr` operations. Here’s one:

    select(skimmed, numeric.mean)

    ## # A tibble: 1 x 1
    ##   numeric.mean
    ##          <dbl>
    ## 1         3.76

To get around this, we’ve added some helper functions and methods. The
more `skimr`-like replacement for `select()` is `focus()`, which
preserves metadata columns.

    focus(skimmed, numeric.mean)

    ## ── Data Summary ────────────────────────
    ##                            Values    
    ## Name                       Piped data
    ## Number of rows             150       
    ## Number of columns          5         
    ## _______________________              
    ## Column type frequency:               
    ##   numeric                  1         
    ## ________________________             
    ## Group variables            None      
    ## 
    ## ── Variable type: numeric ────────────────────────────────────────────────────────────────
    ##   skim_variable  mean
    ## 1 Petal.Length   3.76

Configuring and extending skimr
-------------------------------

Most of `skimr`’s magic, to [steal a
term](https://resources.rstudio.com/rstudio-conf-2019/our-colour-of-magic-the-open-sourcery-of-fantastic-r-packages),
comes from the fact that you can do most everything with one function.
But believe it or not, there’s actually a bit more to the package.

One big one is customization. We like the `skimr` defaults, but that
doesn’t guarantee you will. So what if you want to do something
different, we have a function factory for that!

    my_skim <- skim_with(numeric = sfl(iqr = IQR, p25 = NULL, p75 = NULL))
    my_skim(faithful)

    ## ── Data Summary ────────────────────────
    ##                            Values  
    ## Name                       faithful
    ## Number of rows             272     
    ## Number of columns          2       
    ## _______________________            
    ## Column type frequency:             
    ##   numeric                  2       
    ## ________________________           
    ## Group variables            None    
    ## 
    ## ── Variable type: numeric ────────────────────────────────────────────────────────────────
    ##   skim_variable n_missing complete_rate  mean    sd    p0   p50  p100 hist    iqr
    ## 1 eruptions             0             1  3.49  1.14   1.6     4   5.1 ▇▂▂▇▇  2.29
    ## 2 waiting               0             1 70.9  13.6   43      76  96   ▃▃▂▇▂ 24

Those of you familiar with customizing `skim()` in v1 will notice a
couple differences:

-   we now has an object called `sfl()` for managing `skimr` function
    lists; more below
-   instead of setting global options, we now have a *function factory*

Yes! A function factory. `skim_with()` gives us a new function each time
we call it, and the returned function is configured by the arguments in
`skim_with()`. This works the same way as `ecdf()` in the `stats`
package or `colorRamp` in `grDevices`. Creating new functions has a few
advantages over the previous approach.

-   you can export a `skim()` function in a package or create it in a
    `.Rprofile`
-   you avoid a bunch of potential side effects from setting options
    with `skim_with()`

The other big change is how we now handle different data types. Although
many will never see it, a key piece of `skimr` customization comes from
the `get_skimmers()` generic. It’s used to detect different column types
in your data and set the appropriate summary functions for that type.
It’s also designed to work with `sfl()`. Here’s an example from the
“Supporting additional objects” vignette. Here, we’ll create some
skimmers for
[`sf`](https://cran.r-project.org/web/packages/sf/index.html) data
types:

    get_skimmers.sfc_POINT <- function(column) {
      sfl(
        skim_type = "sfc_POINT",
        n_unique = n_unique,
        valid = ~ sum(sf::st_is_valid(.))
      )
    }

While it was required in `skim_with()`, users must provide a `skim_type`
value when creating new methods. With that, you can export this method
in a new package (be sure to import the generic), and the new default
skimmer is added when you load the package.

    get_default_skimmer_names()

    ...
    $sfc_POINT
    [1] "missing"  "complete" "n"        "n_unique" "valid"
    ...

Even if you don’t go the full route of supporting a new data type,
creating a couple of `skimr` function lists has other benefits. For
example, you can add some to your `.Rprofile` as a way to quickly
configure `skimr` interactively.

    sfc_point_sfl <- sfl(
      n_unique = n_unique,
      valid = ~ sum(sf::st_is_valid(.))
    )

    my_skimmer <- skim_with(sfc_POINT = sfc_point_sfl)

Using skimr in other contexts
-----------------------------

In `skimr` v1, we developed some slightly hacky approaches to getting
nicer `skim()` output in RMarkdown docs. These have been removed in
favor of the
[actually-supported](https://github.com/yihui/knitr/issues/1493)
`knit_print` API. Now, calling `skim()`, within an RMarkdown doc should
produce something nice by default.

    skim(chickwts)

<table>
<caption>Data summary</caption>
<tbody>
<tr class="odd">
<td style="text-align: left;">Name</td>
<td style="text-align: left;">chickwts</td>
</tr>
<tr class="even">
<td style="text-align: left;">Number of rows</td>
<td style="text-align: left;">71</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Number of columns</td>
<td style="text-align: left;">2</td>
</tr>
<tr class="even">
<td style="text-align: left;">_______________________</td>
<td style="text-align: left;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">Column type frequency:</td>
<td style="text-align: left;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">factor</td>
<td style="text-align: left;">1</td>
</tr>
<tr class="odd">
<td style="text-align: left;">numeric</td>
<td style="text-align: left;">1</td>
</tr>
<tr class="even">
<td style="text-align: left;">________________________</td>
<td style="text-align: left;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">Group variables</td>
<td style="text-align: left;">None</td>
</tr>
</tbody>
</table>

**Variable type: factor**

<table>
<thead>
<tr class="header">
<th style="text-align: left;">skim_variable</th>
<th style="text-align: right;">n_missing</th>
<th style="text-align: right;">complete_rate</th>
<th style="text-align: left;">ordered</th>
<th style="text-align: right;">n_unique</th>
<th style="text-align: left;">top_counts</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">feed</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;">FALSE</td>
<td style="text-align: right;">6</td>
<td style="text-align: left;">soy: 14, cas: 12, lin: 12, sun: 12</td>
</tr>
</tbody>
</table>

**Variable type: numeric**

<table>
<thead>
<tr class="header">
<th style="text-align: left;">skim_variable</th>
<th style="text-align: right;">n_missing</th>
<th style="text-align: right;">complete_rate</th>
<th style="text-align: right;">mean</th>
<th style="text-align: right;">sd</th>
<th style="text-align: right;">p0</th>
<th style="text-align: right;">p25</th>
<th style="text-align: right;">p50</th>
<th style="text-align: right;">p75</th>
<th style="text-align: right;">p100</th>
<th style="text-align: left;">hist</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">weight</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">261.31</td>
<td style="text-align: right;">78.07</td>
<td style="text-align: right;">108</td>
<td style="text-align: right;">204.5</td>
<td style="text-align: right;">258</td>
<td style="text-align: right;">323.5</td>
<td style="text-align: right;">423</td>
<td style="text-align: left;">▆▆▇▇▃</td>
</tr>
</tbody>
</table>

You get a nice html version of both the summary header and the `skimr`
subtables for each type of data.

In this context, you configure the output the same way you handle other
`knitr` code chunks.

This means that we’re dropping direct support for `kable.skim_df()` and
`pander.skim_df()`. But you can still get pretty similar results to
these functions by using the reshaping functions described above to get
subtables. You can also still use `Pander` and other nice rendering
packages on an ad hoc basis as you would for other data frames or
tibbles.

We also have a similarly-nice rendered output in
[Jupyter](https://github.com/ropensci/skimr/blob/8c2263c4fd4796af0e5e8f32aafc4980bd58d43a/inst/other_docs/skimr_in_jupyter.ipynb)
and RMarkdown notebooks. In the latter, the summary is separated from
the rest of the output when working interactively. We like it that way,
but we’d be happy to hear what the rest of you think!

Wait, that took over a year?
----------------------------

Well, we think that’s a lot! But to be fair, it wasn’t exactly simple to
keep up with `skimr`. Real talk, open source development takes up a lot
of time, and the `skimr` developers have additional important
priorities. Michael’s family added a new baby, and despite swearing up
and down otherwise, he got absolutely nothing not-baby-related done
during his paternity leave (take note new dads!). Elin ended up taking a
much bigger role on at Lehman, really limiting time for any other work.

Even so, these are just the highlights in the normal ebb and flow of
this sort of work. Since it’s no one’s real job, it might not always be
the first focus. And that’s OK! We’ve been really lucky to have a group
of new users that have been very patient with this slow development
cycle while still providing really good feedback throughout. Thank you
all!

We’re really excited about this next step in the `skimr` journey. We’ve
put a huge amount of work into this new version. Hopefully it shows. And
hopefully it inspires some of you to send more feedback and help us find
even more ways to improve!
Like every R user who uses summary statistics (so, everyone), our team has to rely on some combination of summary functions beyond `summary()` and `str()`. But we found them all lacking in some way because they can be generic, they don't always provide easy-to-operate-on data structures, and they are not pipeable. What we wanted was a frictionless approach for quickly skimming useful and tidy summary statistics as part of a pipeline. And so at [rOpenSci \#unconf17](http://unconf17.ropensci.org/), we developed [`skimr`](https://github.com/ropenscilabs/skimr#skimr).

In a nutshell, `skimr` will create a `skim_df` object that can be further operated upon or that provides a human-readable printout in the console. It presents reasonable default summary statistics for numerics, factors, etc, and lists counts, and missing and unique values. And the momentum is still going, thanks to our awesome team (see below)!

Backstory
---------

The idea for skimr as a project for the \#unconf17 [was proposed by Amelia McNamara](https://github.com/ropensci/unconf17/issues/50) following [discussions on Twitter](https://twitter.com/AmeliaMN/status/774348524653834241) and an [initial package Hadley Wickham](https://github.com/hadley/precis).

Once we were together in Los Angeles, we formed a solid team, set up a Google Doc, a Slack channel, the `ropensci/skimr` repo, and grabbed a whiteboard.

We started off by brainstorming what we liked about existing summary packages and what other features we wanted. We started looking at example data, `mtcars`.

Here's what we liked and disliked, in Amelia's words:

``` r
### "I like what we get here because mpg is numeric so these stats make sense:" 
summary(mtcars$mpg) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  10.40   15.42   19.20   20.09   22.80   33.90 


### "But I don’t like this because cyl should really be a factor and shouldn't have these stats:"
summary(mtcars$cyl)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  4.000   4.000   6.000   6.188   8.000   8.000 


### "This is OK, but not descriptive enough. It could be clearer what I'm looking at."
mosaic::tally(~cyl, data=mtcars) # install.packages('mosaic')
#cyl
# 4  6  8 
#11  7 14 


### "But this output isn't labeled, not ideal." 
table(mtcars$cyl, mtcars$vs)
#   
#     0  1
#  4  1 10
#  6  3  4
#  8 14  0


### "I like this because it returns 'sd', 'n' and 'missing'":
mosaic::favstats(~mpg, data=mtcars) 
#  min     Q1 median   Q3  max     mean       sd  n missing
# 10.4 15.425   19.2 22.8 33.9 20.09062 6.026948 32       0
```

Once we had an idea of what we thought would be useful, we did a bit of market research (i.e. we made a short [presentation](https://docs.google.com/presentation/d/13Ky3-Y70STzufLJtCm6GXN8SMj2Y11riDef8b9cBgAQ/edit#slide=id.p) and ran it by other unconfers at lunch.

Introducing `skimr`
-------------------

So what does `skimr` actually do? It allows you to skim useful summary statistics in the console, or use those statistics in a pipeable workflow.

Some features of output in the console:

-   reports missing, complete, n, sd, and quantiles
-   reports numeric/int/double separately from factor/chr, and identifies class
-   handles dates, logicals
-   uses [Hadley's pillars](https://github.com/hadley/pillar), specifically `pillar::spark-bar()`

Here are examples of `skimr` in action:

### Quick skim in the console:

**Nicely separates numeric and factor variables:**

![](https://github.com/ropenscilabs/skimr/blob/24c733d7e4752c37e46e4c36693da107f42f3f55/man/figures/skim_iris.png) <br>

**Clearly displays many numeric variables:**

![](https://github.com/ropenscilabs/skimr/blob/ecb90e22047d4a1b228bcf471650eb79b733e52e/man/figures/skim_mtcars.png) <br>

**Also works with strings:**

![](https://github.com/ropenscilabs/skimr/blob/ecb90e22047d4a1b228bcf471650eb79b733e52e/man/figures/skim_babynames.png) <br>

### Exploring a skim\_df object

By default `skim` prints beautifully in the console, but it also produces a long, tidy-format `skim_df` object that can be computed on.

``` r
a <-  skim(chickwts)
dim(a)
# [1] 22  5
View(a)
```

<img src="https://github.com/ropenscilabs/skimr/blob/ecb90e22047d4a1b228bcf471650eb79b733e52e/man/figures/skim_chickwts_df.png" width="450px">

### Computing with the skim\_df object

Maybe you just want to skim a specific portion of your data frame. Use skimr with a pipe!

``` r
> skim(mtcars) %>% filter(stat=="hist")
# A tibble: 11 × 5
     var    type  stat      level value
   <chr>   <chr> <chr>      <chr> <dbl>
1    mpg numeric  hist ▂▅▇▇▇▃▁▁▂▂     0
2    cyl numeric  hist ▆▁▁▁▃▁▁▁▁▇     0
3   disp numeric  hist ▇▇▅▁▁▇▃▂▁▃     0
4     hp numeric  hist ▆▆▇▂▇▂▃▁▁▁     0
5   drat numeric  hist ▃▇▂▂▃▆▅▁▁▁     0
6     wt numeric  hist ▂▂▂▂▇▆▁▁▁▂     0
7   qsec numeric  hist ▂▃▇▇▇▅▅▁▁▁     0
8     vs numeric  hist ▇▁▁▁▁▁▁▁▁▆     0
9     am numeric  hist ▇▁▁▁▁▁▁▁▁▆     0
10  gear numeric  hist ▇▁▁▁▆▁▁▁▁▂     0
11  carb numeric  hist ▆▇▂▁▇▁▁▁▁▁     0
```

### Specifying your own statistics

Another possibility is specifying your own statistics to display with `skimr`:

``` r
 funs <- list(iqr = IQR,
    quantile = purrr::partial(quantile, probs = .99))
  skim_with(numeric = funs, append = FALSE)
  skim_v(iris$Sepal.Length)
  
#  A tibble: 2 × 4
#      type     stat level value
#     <chr>    <chr> <chr> <dbl>
# 1 numeric      iqr  .all   1.3
# 2 numeric quantile   99%   7.7
```

Our awesome team
----------------

We had a really fantastic team with diverse backgrounds, and it was really cool how organically everyone found a role for themselves during the development of `skimr`. Between brainstorming sessions, experienced coders began to iteratively develop the code while others worked on documentation and tests, and got more involved. Everyone asked questions and brainstormed together; it was a really welcoming environment. We knew that by the end of the second day of the unconf, we would present our work using only the repo's [README](https://github.com/ropenscilabs/skimr#skimr) file. So we focused on communication throughout the entire development process.

A lot of the heavy lifting at the unconf was done by Michael, Elin, and Eduardo, and Elin has continued leading development in the month since!

This was the original team in alphabetical order. We have also had many virtual contributors as well: see the full list of contributors [here](https://github.com/ropenscilabs/skimr/graphs/contributors).

**Eduardo Arino de la Rubia**
Job Title: Chief Data Scientist at Domino Data Lab
Project Contributions: Coder

**Shannon Ellis**
Job Title: Postdoctoral fellow in the Biostatistics Department at the Johns Hopkins Bloomberg School of Public Health
Project Contributions: Test Scripts

**Julia Stewart Lowndes**
Job Title: Marine Data Scientist at the National Center for Ecological Analysis and Synthesis
Project Contributions: Documentation and test scripts

**Hope McLeod**
Job Title: Data Engineer at Kobalt Music
Project Contributions: Documentation

**Amelia McNamara**
Job Title: Visiting Assistant Professor of Statistical & Data Sciences at Smith College
Project Contributions: Coder

**Michael Quinn**
Job Title: Quantitative Analyst at Google
Project Contributions: Coder

**Elin Waring**
Job Title: Professor at Lehman College Sociology Department, City University of New York
Project Contributions: Coder

**Hao Zhu**
Job Title: Programmer Analyst at the Institute for Aging Research
Project Contributions: Coder

In summary (ha...)
------------------

The work we did together was only possible because of rOpenSci's incredible community and culture. For us to be able to dream up something we wanted to build and have the time and space to actually do it together was really exciting. So thank you rOpenSci and everyone in the greater community!

There is more work to be done on `skimr`, so please check out the [`skimr`](https://github.com/ropenscilabs/skimr) repo for the latest features and improvements!
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 3.5.2 (2018-12-20) |
|os       |Debian GNU/Linux 10 (buster) |
|system   |x86_64, linux-gnu            |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |America/New_York             |
|date     |2019-10-28                   |

# Dependencies

|package    |old      |new      |Δ  |
|:----------|:--------|:--------|:--|
|skimr      |1.0.7    |2.0      |*  |
|assertthat |0.2.1    |0.2.1    |   |
|backports  |1.1.5    |1.1.5    |   |
|base64enc  |NA       |0.1-3    |*  |
|BH         |1.69.0-1 |1.69.0-1 |   |
|cli        |1.1.0    |1.1.0    |   |
|crayon     |1.3.4    |1.3.4    |   |
|digest     |0.6.22   |0.6.22   |   |
|dplyr      |0.8.3    |0.8.3    |   |
|ellipsis   |0.3.0    |0.3.0    |   |
|evaluate   |0.14     |0.14     |   |
|fansi      |0.4.0    |0.4.0    |   |
|glue       |1.3.1    |1.3.1    |   |
|highr      |0.8      |0.8      |   |
|htmltools  |NA       |0.4.0    |*  |
|jsonlite   |NA       |1.6      |*  |
|knitr      |1.25     |1.25     |   |
|lifecycle  |0.1.0    |0.1.0    |   |
|magrittr   |1.5      |1.5      |   |
|markdown   |1.1      |1.1      |   |
|mime       |0.7      |0.7      |   |
|pander     |0.6.3    |NA       |*  |
|pillar     |1.4.2    |1.4.2    |   |
|pkgconfig  |2.0.3    |2.0.3    |   |
|plogr      |0.2.0    |0.2.0    |   |
|purrr      |0.3.3    |0.3.3    |   |
|R6         |2.4.0    |2.4.0    |   |
|Rcpp       |1.0.2    |1.0.2    |   |
|repr       |NA       |1.0.1    |*  |
|rlang      |0.4.1    |0.4.1    |   |
|stringi    |1.4.3    |1.4.3    |   |
|stringr    |1.4.0    |1.4.0    |   |
|tibble     |2.1.3    |2.1.3    |   |
|tidyr      |1.0.0    |1.0.0    |   |
|tidyselect |0.2.5    |0.2.5    |   |
|utf8       |1.1.4    |1.1.4    |   |
|vctrs      |0.2.0    |0.2.0    |   |
|xfun       |0.10     |0.10     |   |
|yaml       |2.2.0    |2.2.0    |   |
|zeallot    |0.1.0    |0.1.0    |   |

# Revdeps

## Failed to check (4)

|package      |version |error |warning |note |
|:------------|:-------|:-----|:-------|:----|
|codebook     |?       |      |        |     |
|groupedstats |?       |      |        |     |
|GSODR        |?       |      |        |     |
|panelr       |?       |      |        |     |

*Wow, no problems at all. :)*# codebook

<details>

* Version: 
* Source code: ???
* URL: https://github.com/ropenscilabs/skimr
* BugReports: https://github.com/ropenscilabs/skimr/issues
* Number of recursive dependencies: 0

Run `revdep_details(,"")` for more info

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
# groupedstats

<details>

* Version: 
* Source code: ???
* URL: https://github.com/ropenscilabs/skimr
* BugReports: https://github.com/ropenscilabs/skimr/issues
* Number of recursive dependencies: 0

Run `revdep_details(,"")` for more info

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
# GSODR

<details>

* Version: 
* Source code: ???
* URL: https://github.com/ropenscilabs/skimr
* BugReports: https://github.com/ropenscilabs/skimr/issues
* Number of recursive dependencies: 0

Run `revdep_details(,"")` for more info

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
# panelr

<details>

* Version: 
* Source code: ???
* URL: https://github.com/ropenscilabs/skimr
* BugReports: https://github.com/ropenscilabs/skimr/issues
* Number of recursive dependencies: 0

Run `revdep_details(,"")` for more info

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
