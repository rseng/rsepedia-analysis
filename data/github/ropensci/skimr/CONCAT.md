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
---
output: md_document
---
<!-- README.md is generated from README.Rmd. Please edit that file -->
# skimr <a href='https://docs.ropensci.org/skimr/'>
<img src='https://docs.ropensci.org/skimr/reference/figures/logo.png'
align="right" height="139" /></a>

```{r set-options, echo=FALSE, message=FALSE}
library(skimr)
options(tibble.width = Inf)
options(width = 100)
```

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


`skimr` provides a frictionless approach to summary statistics which conforms
to the [principle of least
surprise](https://en.wikipedia.org/wiki/Principle_of_least_astonishment),
displaying summary statistics the user can skim quickly to understand their
data. It handles different data types and returns a `skim_df` object which can
be included in a pipeline or displayed nicely for the human reader.

**Note: `skimr` version 2 has major changes when skimr is used programmatically.
Upgraders should review this document, the release notes and vignettes
carefully.**

## Installation

The current released version of `skimr` can be installed from CRAN. If you wish
to install the current build of the next release you can do so using the
following:

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/skimr")
```

The APIs for this branch should be considered reasonably stable but still
subject to change if an issue is discovered.

To install the version with the most recent changes that have not yet been
incorporated in the master branch (and may not be):

```{r, eval = FALSE}
devtools::install_github("ropensci/skimr", ref = "develop")
```

Do not rely on APIs from the develop branch, as they are likely to change.

## Skim statistics in the console

`skimr`:

- Provides a larger set of statistics than `summary()`, including missing,
  complete, n, and sd.
- reports each data types separately
- handles dates, logicals, and a variety of other types
- supports spark-bar and spark-line based on the
  [pillar package](https://github.com/r-lib/pillar).

### Separates variables by class:

```{r, render = knitr::normal_print}
skim(chickwts)
```

### Presentation is in a compact horizontal format:

```{r, render = knitr::normal_print}
skim(iris)
```

### Built in support for strings, lists and other column classes

```{r, render = knitr::normal_print}
skim(dplyr::starwars)
```

### Has a useful summary function

```{r, render = knitr::normal_print}
skim(iris) %>%
  summary()
```

### Individual columns can be selected using tidyverse-style selectors

```{r, render = knitr::normal_print}
skim(iris, Sepal.Length, Petal.Length)
```

### Handles grouped data

`skim()` can handle data that has been grouped using `dplyr::group_by()`.

```{r, render = knitr::normal_print}
iris %>%
  dplyr::group_by(Species) %>%
  skim()
```

### Behaves nicely in pipelines

```{r, render = knitr::normal_print}
iris %>%
  skim() %>%
  dplyr::filter(numeric.sd > 1)
```

## Knitted results

Simply skimming a data frame will produce the horizontal print
layout shown above. We provide a `knit_print` method for the types of objects
in this package so that similar results are produced in documents. To use this,
make sure the `skimmed` object is the last item in your code chunk.

```{r}
faithful %>%
  skim()
```

## Customizing skimr

Although skimr provides opinionated defaults, it is highly customizable.
Users can specify their own statistics, change the formatting of results,
create statistics for new classes and develop skimmers for data structures
that are not data frames.

### Specify your own statistics and classes

Users can specify their own statistics using a list combined with the
`skim_with()` function factory. `skim_with()` returns a new `skim` function that
can be called on your data. You can use this factory to produce summaries for
any type of column within your data.

Assignment within a call to `skim_with()` relies on a helper function, `sfl` or
`skimr` function list. By default, functions in the `sfl` call are appended to
the default skimmers, and names are automatically generated as well.

```{}
my_skim <- skim_with(numeric = sfl(mad))
my_skim(iris, Sepal.Length)
```

But you can also helpers from the `tidyverse` to create new anonymous functions
that set particular function arguments. The behavior is the same as in `purrr`
or `dplyr`, with both `.` and `.x` as acceptable pronouns. Setting the
`append = FALSE` argument uses only those functions that you've provided.

```{}
my_skim <- skim_with(
  numeric = sfl(
    iqr = IQR,
    p01 = ~ quantile(.x, probs = .01)
    p99 = ~ quantile(., probs = .99)
  ),
  append = FALSE
)
my_skim(iris, Sepal.Length)
```

And you can remove default skimmers by setting them to `NULL`.

```{}
my_skim <- skim_with(numeric = sfl(hist = NULL))
my_skim(iris, Sepal.Length)
```

### Skimming other objects

`skimr` has summary functions for the following types of data by default:

* `numeric` (which includes both `double` and `integer`)
* `character`
* `factor`
* `logical`
* `complex`
* `Date`
* `POSIXct`
* `ts`
* `AsIs`

`skimr` also provides a small API for writing packages that provide their own
default summary functions for data types not covered above. It relies on
R S3 methods for the `get_skimmers` function. This function should return
a `sfl`, similar to customization within `skim_with()`, but you should also
provide a value for the `class` argument. Here's an example.

```{r}
get_skimmers.my_data_type <- function(column) {
  sfl(
    .class = "my_data_type",
    p99 = quantile(., probs = .99)
  )
}
```

## Limitations of current version

We are aware that there are issues with rendering the inline histograms and
line charts in various contexts, some of which are described below.

### Support for spark histograms

There are known issues with printing the spark-histogram characters when
printing a data frame. For example, `"▂▅▇"` is printed as
`"<U+2582><U+2585><U+2587>"`. This longstanding problem [originates in
the low-level
code](https://r.789695.n4.nabble.com/Unicode-display-problem-with-data-frames-under-Windows-td4707639.html)
for printing dataframes.
While some cases have been addressed, there are, for example, reports of this
issue in Emacs ESS.

This means that while `skimr` can render the histograms to the console and in
RMarkdown documents, it cannot in other circumstances. This includes:

* converting a `skimr` data frame to a vanilla R data frame, but tibbles render
  correctly
* in the context of rendering to a pdf using an engine that does not support
  utf-8.

One workaround for showing these characters in Windows is to set the CTYPE part
of your locale to Chinese/Japanese/Korean with `Sys.setlocale("LC_CTYPE",
"Chinese")`. The helper function `fix_windows_histograms()` does this for you.

And last but not least, we provide `skim_without_charts()` as a fallback.
This makes it easy to still get summaries of your data, even if unicode issues
continue.

### Printing spark histograms and line graphs in knitted documents

Spark-bar and spark-line work in the console, but may not work when you knit
them to a specific document format. The same session that produces a correctly
rendered HTML document may produce an incorrectly rendered PDF, for example.
This issue can generally be addressed by changing fonts to one with good
building block (for histograms) and Braille support (for line graphs). For
example, the open font "DejaVu Sans" from the `extrafont` package supports
these. You may also want to try wrapping your results in `knitr::kable()`.
Please see the vignette on using fonts for details.

Displays in documents of different types will vary. For example, one user found
that the font "Yu Gothic UI Semilight" produced consistent results for
Microsoft Word and Libre Office Write.

### Stripping metadata and empty results tables

In POSIX systems, `skimr` tries to remove the tibble metadata when producing
the results. A complicating factor is tibble's color support, which depends
on environment settings. In particular, not all Windows terminals support
colors in the way that tibble expects.

So, by default, we disable removing metadata on windows. You can turn this
feature on with an option. Either set it when calling print or globally.

```{r, eval = FALSE}
skimmed <- skim(chickwts)
print(skimmed, strip_metadata = TRUE)
options(skimr_strip_metadata = TRUE)
```

Separately, you might need to check the option `crayon.enabled`. Similarly, if
your skimr results tables are empty you may need to run the following

```{r, eval = FALSE}
options(crayon.enabled = FALSE)
```

You need to do this one time per session.

## Inspirations

* [TextPlots](https://github.com/sunetos/TextPlots.jl) for use of Braille
  characters

* [spark](https://github.com/holman/spark) for use of block characters.

The earliest use of unicode characters to generate sparklines appears to be [from 2009](https://blog.jonudell.net/2009/01/13/fuel-prices-and-pageviews/).

Exercising these ideas to their fullest requires a font with good support for block drawing characters. [PragamataPro](https://fsd.it/shop/fonts/pragmatapro/) is one such font.

## Contributing

We welcome issue reports and pull requests, including potentially adding
support for commonly used variable classes. However, in general, we encourage
users to take advantage of skimr's flexibility to add their own customized
classes. Please see the
[contributing](https://docs.ropensci.org/skimr/CONTRIBUTING.html) and
[conduct](https://ropensci.org/code-of-conduct/) documents.

[![ropenci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Untitled"
mainfont: DejaVu Sans
output:
  html_document: default
  pdf_document:
    latex_engine: xelatex
  word_document: default
font-family: Times New Roman
---

## Getting ready for fonts

Notice that the yaml at the beginning of this file includes a latex_engine which will be used 
when creating a pdf document.

It also includes a mainfont setting called _DejaVu Sans_.  This is not the only font that will 
work to produce the spark graphs. However, it is a free font available through the
`extrafont` package.  If you have not installed extrafont you should do so using the normal
package installation procedures. You should then make sure that the desired font is installed.

The code below will not run automatically when you knit, instead you should run it in the
console. 

```
install.packages(c("extrafont"))
font_install("DejaVu Sans")
```

If there are any difficulties please read the extrafont documentation.

## Experimenting with rendering
 
```{r, message=FALSE}
library(knitr)
library(skimr)

```

Try knitting this document to PDF, HTML, doc or any other format you wish to try. You 
will notice that there are slight differences between them.  To understand the impact
of the engine and font choices you should experiment with different options.

The first example shows what printing the basic skim function looks like. 
You can try knitting to different formats to see how it changes.

```{r }
skim(iris)
```

It is possible that the histograms will not print in all of the formats.

Unfortunately this is outside the control of the skimr team because 
it relates to the operating system you are using, fonts installed, 
and locale. 

---
slug: skimrv2
title: "(Re)introducing `skimr` v2"
subtitle: "A year in the life of an open source R project"
package_version: 2.0
authors:
  - name: Michael Quinn
    url: https://github.com/michaelquinn32
  - name: Elin Waring
    url: https://github.com/elinw
date: "10/09/2019"
categories: blog
tags:
  - software peer review
  - R
  - packages
  - community
  - skimr
---

[Theme song: *PSA* by Jay-Z](https://www.youtube.com/watch?v=-LzdKH1naok)

We announced the testing version of `skimr` v2 on
[June 19, 2018](https://github.com/ropensci/skimr/issues/341). After more than
a year of (admittedly intermittent) work, we're thrilled to be able to say that
the package is ready to go to CRAN. So, what happened over the last year? And
why are we so excited for V2?

## Setting the stage

Before we can talk about the last year of `skimr` development, we need to lay
out the timeline that got us to this point. For those deeply enmeshed in `skimr`
lore, all [dozens](https://imgur.com/gallery/R1fdEt3) of you, bear with.

`skimr` was originally an [rOpenSci
unconf17](https://ropensci.org/blog/2017/07/11/skimr/) project, a big
collaboration between eight different participants that resulted in a conceptual
outline of the package and a basic working version. Participating in the unconf
was a truly magical experience, with everyone bringing a tremendous amount of
energy and ideas to the project, and implementation happening over a flurry of
["fancy git commits"](https://twitter.com/AmeliaMN/status/867818976666976256).

About six months later, we released our first version on CRAN. The time between
these two milestones was mostly spent on fleshing out all of the different ideas
that were generated during the unconf (like handling grouped data frames) and
fixing all the bugs we discovered along the way.

Getting the package on CRAN opened the gates for bug reports and feature
requests on [GitHub](https://github.com/ropensci/skimr/issues). About the same
time we pushed our first version to CRAN, Elin got `skimr`'s rOpenSci's package
[peer review](https://github.com/ropensci/software-review/issues/175) started
(thank you Jennifer and Jim!), opening another incredibly useful channel for
collecting feedback on the package. All of these new ideas and suggestions gave
us the opportunity to really push `skimr` to the next level, but doing that
would require rethinking the package, from the ground up.

A month after finishing the peer review (and six months after the process
began), we announced v2. Over the first phase of `skimr`'s life, we accumulated
700 commits, two release, 400 GitHub stars, 95 percent code coverage and a
lifetime's worth of [unicode rendering
bugs](https://github.com/ropensci/skimr#support-for-spark-histograms)!

Just kidding! We love our little histograms, even when they don't love us back!
For those of you that might have never seen `skimr`, using the package typically
boils down to a single function call:

```{r render = knitr::normal_print, message=FALSE}
library(skimr)
library(dplyr)
options(width = 90)

skim(iris)
```

## Getting it right

Under normal circumstances (i.e. not during a hackathon), most software
engineering projects begin with a design phase and series of increasingly
detailed design docs. `skimr` is only a few hundred lines of code, which means
"increasingly detailed design docs" translates to one doc. But we did actually
write it! [It's
here](https://docs.google.com/document/d/18lBStDZzd1rJq08O-4Sw2qHhuHEZ79QX4sBkeyzWNFY/edit#heading=h.5x0d5h95i329).
And it still goes a good job of laying out some of the big ideas we were
interested in taking on for v2.

* Eliminating frictions that resulted from differences in the way we stored data
  vs how it was displayed to users
* Getting away from using a global environment to configure `skimr`
* Making it easier for others to extend `skimr`
* Create more useful ways to use `skimr`

## Better internal data structures

In v1, `skimr` stored all of its data in a "long format", data frame. Although
hidden from the user by its print methods, this format would appear any time
you'd try do something with the results of a `skim()` call. It looked something
like this:

```{r eval = FALSE}
skim(mtcars) %>% dplyr::filter(stat=="hist")
```
```
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
```

Big ups to anyone who looked at the rendered output and saw that this was how
you actually filtered the results. Hopefully there are even better applications
of your near-telepathic abilities.

Now, working with `skimr` is a bit more sane.

```{r render = knitr::normal_print}
skimmed <- iris %>%
  skim() %>%
  dplyr::filter(numeric.sd > 1)

skimmed
```

And

```{r render = knitr::normal_print}
dplyr::glimpse(skimmed)
```

It's still not perfect, as you need to rely on a *pseudo-namespace* to refer to
the column that you want. But this is unfortunately a necessary trade-off. As
the Rstats Bible, errr Hadley Wickham's *Advanced R*, states, all elements of
[an atomic vector must have the same
type](https://adv-r.hadley.nz/vectors-chap.html). This normally isn't something
that you have to think too much about, that is until you try to combine the
means of all your `Date` columns with the means of your `numeric` columns and
everything comes out utterly garbled. So instead of that basket of laughs, we
prefix columns names by their data type.

There's a couple of other nuances here:

* The data frame `skim()` produces always starts off with some metadata columns
* Functions that always produce the same, regardless of input type, can
  be treated as `base_skimmers` and don't need a namespace

### Manipulating internal data

A better representation of internal data comes with better tools for reshaping
the data and getting it for other contexts. A common request in v1 was tooling
to handle the `skimr` subtables separately. We now do this with `partition()`.
It replaces the v1 function `skim_to_list()`.

```{r render = knitr::normal_print}
partition(skimmed)
```

You can undo a call to `partition()` with `bind()`, which joins the subtables
into the original `skim_df` object and properly accounts for metadata. You
can skip a step with the function `yank()`, which calls partition and pulls
out a particular subtable


```{r render = knitr::normal_print}
yank(skimmed, "numeric")
```

Last, with support something close to the older format with the `to_long()`
function. This can be added for something close to backwards compatibility.
Being realistic on open source sustainability means that we are not able to
support 100% backward compatibility in v2 even with new functions.  Meanwhile
you can keep using v1 if you are happy with it.  However,
because `skimr`'s dependencies are under ongoing development, sooner or later
skimr v1 will no longer work with updates to them.

### Working with dplyr

Using `skimr` in a `dplyr` pipeline was part of the original package design, and
we've needed to devote some extra love to making sure that everything is as
seamless as possible. Part of this is due to the object produce by `skim()`,
which we call `skim_df`. It's a little weird in that it needs both metadata and
columns in the underlying data frame.

In practice, this means that you can coerce it into a different type through
normal `dplyr` operations. Here's one:

```{r render = knitr::normal_print}
select(skimmed, numeric.mean)
```

To get around this, we've added some helper functions and methods. The more
`skimr`-like replacement for `select()` is `focus()`, which preserves
metadata columns.

```{r render = knitr::normal_print}
focus(skimmed, numeric.mean)
```

## Configuring and extending skimr

Most of `skimr`'s magic, to [steal a
term](https://resources.rstudio.com/rstudio-conf-2019/our-colour-of-magic-the-open-sourcery-of-fantastic-r-packages),
comes from the fact that you can do
most everything with one function. But believe it or not, there's actually a bit
more to the package.

One big one is customization. We like the `skimr` defaults, but that doesn't
guarantee you will. So what if you want to do something different, we have
a function factory for that!

```{r render = knitr::normal_print}
my_skim <- skim_with(numeric = sfl(iqr = IQR, p25 = NULL, p75 = NULL))
my_skim(faithful)
```

Those of you familiar with customizing `skim()` in v1 will notice a couple
differences:

* we now has an object called `sfl()` for managing `skimr` function lists; more
  below
* instead of setting global options, we now have a *function factory*

Yes! A function factory. `skim_with()` gives us a new function each time
we call it, and the returned function is configured by the arguments in
`skim_with()`. This works the same way as `ecdf()` in the `stats` package or
`colorRamp` in `grDevices`. Creating new functions has a few advantages over
the previous approach.

* you can export a `skim()` function in a package or create it in a `.Rprofile`
* you avoid a bunch of potential side effects from setting options with
  `skim_with()`

The other big change is how we now handle different data types. Although many
will never see it, a key piece of `skimr` customization comes from the
`get_skimmers()` generic. It's used to detect different column types in your
data and set the appropriate summary functions for that type. It's also
designed to work with `sfl()`. Here's an example from the "Supporting additional
objects" vignette. Here, we'll create some skimmers for
[`sf`](https://cran.r-project.org/web/packages/sf/index.html) data types:

```{r eval = FALSE}
get_skimmers.sfc_POINT <- function(column) {
  sfl(
    skim_type = "sfc_POINT",
    n_unique = n_unique,
    valid = ~ sum(sf::st_is_valid(.))
  )
}
```

While it was required in `skim_with()`, users must provide a `skim_type` value
when creating new methods. With that, you can export this method in a new
package (be sure to import the generic), and the new default skimmer is added
when you load the package.

```{r, eval = FALSE}
get_default_skimmer_names()
```
```
...
$sfc_POINT
[1] "missing"  "complete" "n"        "n_unique" "valid"
...
```

Even if you don't go the full route of supporting a new data type, creating a
couple of `skimr` function lists has other benefits. For example, you can add
some to your `.Rprofile` as a way to quickly configure `skimr` interactively.

```{r eval = FALSE}
sfc_point_sfl <- sfl(
  n_unique = n_unique,
  valid = ~ sum(sf::st_is_valid(.))
)

my_skimmer <- skim_with(sfc_POINT = sfc_point_sfl)
```

## Using skimr in other contexts

In `skimr` v1, we developed some slightly hacky approaches to getting nicer
`skim()` output in RMarkdown docs. These have been removed in favor of the
[actually-supported](https://github.com/yihui/knitr/issues/1493) `knit_print`
API. Now, calling `skim()`, within an RMarkdown doc should produce something
nice by default.

```{r}
skim(chickwts)
```

You get a nice html version of both the summary header and the `skimr` subtables
for each type of data.

In this context, you configure the output the same way you handle other `knitr`
code chunks.

~~~
```{r skimr_digits = 4, skimr_summary = TRUE}
```
~~~

This means that we're dropping direct support for `kable.skim_df()` and
`pander.skim_df()`. But you can still get pretty similar results to these
functions by using the reshaping functions described above to get subtables. You
can also still use `Pander` and other nice rendering packages on an ad hoc basis
as you would for other data frames or tibbles.

We also have a similarly-nice rendered output in
[Jupyter](https://github.com/ropensci/skimr/blob/8c2263c4fd4796af0e5e8f32aafc4980bd58d43a/inst/other_docs/skimr_in_jupyter.ipynb)
and RMarkdown notebooks. In the latter, the summary is separated from the rest
of the output when working interactively. We like it that way, but we'd be happy
to hear what the rest of you think!

## Wait, that took over a year?

Well, we think that's a lot! But to be fair, it wasn't exactly simple to keep up
with `skimr`. Real talk, open source development takes up a lot of time, and
the `skimr` developers have additional important priorities.  Michael's family
added a new baby, and despite swearing up and down otherwise,
he got absolutely nothing not-baby-related done during his paternity leave (take
note new dads!). Elin
ended up taking a much bigger role on at Lehman, really limiting time for any
other work.

Even so, these are just the highlights in the normal ebb and flow of this sort
of work. Since it's no one's real job, it might not always be the first focus.
And that's OK! We've been really lucky to have a group of new users that have
been very patient with this slow development cycle while still providing really
good feedback throughout. Thank you all!

We're really excited about this next step in the `skimr` journey. We've put a
huge amount of work into this new version. Hopefully it shows. And hopefully
it inspires some of you to send more feedback and help us find even more ways
to improve!
---
title: "`skimr` for useful and tidy summary statistics"
authors:
  - name: Eduardo Arino de la Rubia
    url: http://github.com/earino
  - name: Shannon Ellis
    url: http://github.com/ShanEllis
  - name: Julia Stewart Lowndes
    url: http://github.com/jules32
  - name: Hope McLeod
    url: http://github.com/homcl
  - name: Amelia McNamara
    url: http://github.com/AmeliaMN
  - name: Michael Quinn
    url: https://github.com/michaelquinn32
  - name: Elin Waring
    url: https://github.com/elinw
  - name: Hao Zhu
    url: https://github.com/haozhu233
date: "07/11/2017"
output:
  md_document:
    variant: markdown_github
categories:
  - blog
tags:
  - unconf
  - testing
  - review
  - packages
  - software
  - community
  - meetings
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Like every R user who uses summary statistics (so, everyone), our team has to rely on some combination of summary functions beyond `summary()` and `str()`. But we found them all lacking in some way because they can be generic, they don't always provide easy-to-operate-on data structures, and they are not pipeable. What we wanted was a frictionless approach for quickly skimming useful and tidy summary statistics as part of a pipeline. And so at [rOpenSci #unconf17](http://unconf17.ropensci.org/), we developed [`skimr`](https://github.com/ropenscilabs/skimr#skimr). 

In a nutshell, `skimr` will create a `skim_df` object that can be further operated upon or that provides a human-readable printout in the console. It presents reasonable default summary statistics for numerics, factors, etc, and lists counts, and missing and unique values. And the momentum is still going, thanks to our awesome team (see below)!

## Backstory

The idea for skimr as a project for the #unconf17 [was proposed by Amelia McNamara](https://github.com/ropensci/unconf17/issues/50) following [discussions on Twitter](https://twitter.com/AmeliaMN/status/774348524653834241) and an [initial package Hadley Wickham](https://github.com/hadley/precis). 

Once we were together in Los Angeles, we formed a solid team, set up a Google Doc, a Slack channel, the `ropensci/skimr` repo, and grabbed a whiteboard.  

We started off by brainstorming what we liked about existing summary packages and what other features we wanted. We started looking at example data, `mtcars`.

Here's what we liked and disliked, in Amelia's words:

```{r summary, eval=FALSE}
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

## Introducing `skimr`

So what does `skimr` actually do? It allows you to skim useful summary statistics in the console, or use those statistics in a pipeable workflow.

Some features of output in the console: 

- reports missing, complete, n, sd, and quantiles
- reports numeric/int/double separately from factor/chr, and identifies class
- handles dates, logicals
- uses [Hadley's pillars](https://github.com/hadley/pillar), specifically `pillar::spark-bar()`

Here are examples of `skimr` in action: 

### Quick skim in the console: 

**Nicely separates numeric and factor variables:**

![](https://github.com/ropenscilabs/skimr/blob/24c733d7e4752c37e46e4c36693da107f42f3f55/man/figures/skim_iris.png)
<br>

**Clearly displays many numeric variables:**

![](https://github.com/ropenscilabs/skimr/blob/ecb90e22047d4a1b228bcf471650eb79b733e52e/man/figures/skim_mtcars.png)
<br>

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


### Computing with the skim_df object

Maybe you just want to skim a specific portion of your data frame. Use skimr with a pipe!

```r
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

```r
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

## Our awesome team

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

## In summary (ha...)

The work we did together was only possible because of rOpenSci's incredible community and culture. For us to be able to dream up something we wanted to build and have the time and space to actually do it together was really exciting. So thank you rOpenSci and everyone in the greater community!

There is more work to be done on `skimr`, so please check out the [`skimr`](https://github.com/ropenscilabs/skimr) repo for the latest features and improvements!

---
title: "Using Skimr"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using Skimr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

`skimr` is designed to provide summary statistics about variables in data frames,
tibbles, data tables and vectors. It is
opinionated in its defaults, but easy to modify.

In base R, the most similar functions are `summary()` for vectors and data
frames and `fivenum()` for numeric  vectors:

```{r}
summary(iris)
```

```{r}
summary(iris$Sepal.Length)
```

```{r}
fivenum(iris$Sepal.Length)
```

```{r}
summary(iris$Species)
```

## The `skim()` function

The core function of `skimr` is `skim()`, which is designed to work with
(grouped) data frames, and will try coerce other objects to data frames 
if possible. Like `summary()`, `skim()`'s method for data frames presents 
results for every column; the statistics it provides depend on the class of 
the variable.

### Skimming data frames

By design, the main focus of `skimr` is on data frames; it is intended to fit
well within a data [pipeline](https://r4ds.had.co.nz/pipes.html) and relies
extensively on [tidyverse](https://www.tidyverse.org/) vocabulary, which
focuses on data frames.

Results of `skim()` are *printed* horizontally, with one section per variable
type and one row per variable.

```{r, render = knitr::normal_print}
library(skimr)
skim(iris)
```

The format of the results are a single wide data frame combining the results, 
with some additional attributes and two metadata columns:

- `skim_variable`: name of the original variable
- `skim_type`: class of the variable

Unlike many other objects within `R`, these columns are intrinsic to the
`skim_df` class. Dropping these variables will result in a coercion to a
`tibble`. The `is_skim_df()` function is used to assert that an object is
a skim_df.

```{r}
skim(iris) %>% is_skim_df()

```

```{r, render = knitr::normal_print}
skim(iris) %>%
  dplyr::select(-skim_type, -skim_variable) %>% is_skim_df()
```

```{r, render = knitr::normal_print}
skim(iris) %>%
  dplyr::select(-n_missing) %>% is_skim_df()
```

In order to avoid type coercion, columns for summary statistics for different
types are prefixed with the corresponding `skim_type`. This means that the
columns of the `skim_df` are somewhat sparse, with quite a few missing
values. This is because for some statistics the representations for different
types of variables is different. For example, the mean of a Date variable and
of a numeric variable are represented differently when printing, but this 
cannot be supported in a single vector. The exception to this are 
`n_missing` and `complete_rate` (missing/number of observations) which are the 
same for all types of variables.

```{r, render = knitr::normal_print}
skim(iris) %>%
  tibble::as_tibble()
```

This is in contrast to `summary.data.frame()`, which stores statistics in a
`table`. The distinction is important, because the `skim_df` object is pipeable
and easy to use for additional manipulation: for example, the user could select
all of the variable means, or all summary statistics for a specific variable.

```{r, render = knitr::normal_print}
skim(iris) %>%
  dplyr::filter(skim_variable == "Petal.Length")
```

Most `dplyr` verbs should work as expected.

```{r, render = knitr::normal_print}
skim(iris) %>%
  dplyr::select(skim_type, skim_variable, n_missing)
```

The base skimmers `n_missing` and `complete_rate` are computed for all of the
columns in the data. But all other type-based skimmers have a namespace. You
need to use a `skim_type` prefix to refer to correct column.

```{r, render = knitr::normal_print}
skim(iris) %>%
  dplyr::select(skim_type, skim_variable, numeric.mean)
```

`skim()` also supports grouped data created by `dplyr::group_by()`. 
In this case, one additional column for each grouping variable is added 
to the `skim_df` object.

```{r, render = knitr::normal_print}
iris %>%
  dplyr::group_by(Species) %>%
  skim()
```

Individual columns from a data frame may be selected using tidyverse-style
selectors.

```{r, render = knitr::normal_print}
skim(iris, Sepal.Length, Species)
```

Or with common `select` helpers.

```{r, render = knitr::normal_print}
skim(iris, starts_with("Sepal"))
```

If an individual column is of an unsupported class, it is treated as a
character variable with a warning.

## Skimming vectors

In `skimr` v2, `skim()` will attempt to coerce non-data frames (such as vectors
and matrices) to data frames. In most cases with vectors, the object being
evaluated should be equivalent to wrapping the object in `as.data.frame()`.

For example, the `lynx` data set is class `ts`.

```{r, render = knitr::normal_print}
skim(lynx)
```

Which is the same as coercing to a data frame.

```{r}
all.equal(skim(lynx), skim(as.data.frame(lynx)))
```

## Skimming matrices

`skimr` does not support skimming matrices directly but coerces them to data
frames. Columns in the matrix become variables. This behavior is similar to
`summary.matrix()`). Three possible ways to handle matrices  with `skim()`
parallel the three variations of the mean function for matrices.

```{r, render = knitr::normal_print}
m <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 4, ncol = 3)
m
```

Skimming the matrix produces similar results to `colMeans()`.

```{r, render = knitr::normal_print}
colMeans(m)
skim(m) # Similar to summary.matrix and colMeans()
```

Skimming the transpose of the matrix will give row-wise results.

```{r, render = knitr::normal_print}
rowMeans(m)
skim(t(m))
```

And call `c()` on the matrix to get results across all columns.

```{r, render = knitr::normal_print}
skim(c(m))
mean(m)
```

### Skimming without modification

`skim_tee()` produces the same printed version as `skim()` but returns the
original, unmodified data frame. This allows for continued piping of the
original data.

```{r, render = knitr::normal_print}
iris_setosa <- iris %>%
  skim_tee() %>%
  dplyr::filter(Species == "setosa")
head(iris_setosa)
```

Note, that `skim_tee()` is customized differently than `skim` itself. See below
for more details.

## Reshaping the results from `skim()`

As noted above, `skim()` returns a wide data frame. This is usually the most
sensible format for the majority of operations when investigating data, but
the package has some other functions to help with edge cases.

First, `partition()` returns a named list of the wide data frames for each data
type. Unlike the original data the partitioned data only has columns
corresponding to the skimming functions used for this data type. These data 
frames are, therefore, not `skim_df` objects.

```{r, render = knitr::normal_print}
iris %>%
  skim() %>%
  partition()
```

Alternatively, `yank()` selects only the subtable for a specific type. Think of
it like `dplyr::select` on column types in the original data. Again, unsuitable
columns are dropped.

```{r, render = knitr::normal_print}
iris %>%
  skim() %>%
  yank("numeric")
```

`to_long()` returns a single long data frame with columns `variable`, `type`,
`statistic` and `formatted`. This is similar but not identical to the `skim_df`
object in `skimr` v1.

```{r, render = knitr::normal_print}
iris %>%
  skim() %>%
  to_long() %>% 
  head()
```

Since the `skim_variable` and `skim_type` columns are a core component of the
`skim_df` class, it's possible to get unwanted side effects when using
`dplyr::select()`. Instead, use `focus()` to select columns of the skimmed
results and keep them as a `skim_df`; it always keeps the metadata column.

```{r, render = knitr::normal_print}
iris %>%
  skim() %>%
  focus(n_missing, numeric.mean)
```

## Rendering the results of `skim()`

The `skim_df` object is a wide data frame. The display is
created by default using `print.skim_df()`; users can specify additional
options by explicitly calling `print([skim_df object], ...)`.

For documents rendered by `knitr`, the package provides a custom `knit_print`
method. To use it, the final line of your code chunk should have a `skim_df`
object.

```{r}
skim(Orange)
```

The same type of rendering is available from reshaped `skim_df` objects, those
generated by `partition()` and `yank()` in particular.

```{r}
skim(Orange) %>%
  yank("numeric")
```

## Customizing print options

Although its not a common use case outside of writing vignettes about `skimr`,
you can fall back to default printing methods by adding the chunk option
`render = knitr::normal_print`.

You can also disable the `skimr` summary by setting the chunk option
`skimr_include_summary = FALSE`.

You can change the number of digits shown in the columns of generated statistics
by changing the `skimr_digits` chunk option.

## Modifying `skim()`

`skimr` is opinionated in its choice of defaults, but users can easily add,
replace, or remove the statistics for a class. For interactive use, you can
create your own skimming function with the `skim_with()` factory. `skimr` also
has an API for extensions in other packages. Working with that is covered later.

To add a statistic for a data type, create an `sfl()` (a `skimr` function list)
for each class that you want to change:

```{r}
my_skim <- skim_with(numeric = sfl(new_mad = mad))
my_skim(faithful)
```

As the previous example suggests, the default is to append *new* summary
statistics to the preexisting set. This behavior isn't always desirable,
especially when you want lots of changes. To stop appending, set
`append = FALSE`.

```{r}
my_skim <- skim_with(numeric = sfl(new_mad = mad), append = FALSE)
my_skim(faithful)
```

You can also use `skim_with()` to remove specific statistics by setting them to
`NULL`. This is commonly used to disable the inline histograms and spark graphs.

```{r}
no_hist <- skim_with(ts = sfl(line_graph = NULL))
no_hist(Nile)
```

The same pattern applies to changing skimmers for multiple classes
simultaneously. If you want to partially-apply function arguments, use the
Tidyverse lambda syntax.

```{r}
my_skim <- skim_with(
  numeric = sfl(total = ~ sum(., na.rm = TRUE)),
  factor = sfl(missing = ~ sum(is.na(.))),
  append = FALSE
)

my_skim(iris)
```

To modify the "base" skimmers, refer to them in a similar manner. Since base
skimmers are usually a small group, they must return the same type for all
data types in R, `append` doesn't apply here.

```{r}
my_skim <- skim_with(base = sfl(length = length))
my_skim(faithful)
```

## Extending `skimr`

Packages may wish to export their own `skim()` functions. Use `skim_with()` for
this. In fact, this is how `skimr` generates its version of `skim()`.

```{r}
#' @export
my_package_skim <- skim_with()
```

Alternatively, defaults for another data types can be added to `skimr` with the
`get_skimmers` generic. The method for your data type should return an `sfl()`.
Unlike the `sfl()` used interactively, you also need to set the `skim_type`
argument. It should match the method type in the function signature.

```{r}
get_skimmers.my_type <- function(column) {
  sfl(
    skim_type = "my_type",
    total = sum
  )
}

my_data <- data.frame(
  my_type = structure(1:3, class = c("my_type", "integer"))
)
skim(my_data)
```

An extended example is available in the vignette *Supporting additional
objects*.

## Solutions to common rendering problems

The details of rendering are dependent on the operating system R is running on,
the locale of the installation, and the fonts installed. Rendering may also
differ based on whether it occurs in the console or when knitting to specific
types of documents such as HTML and PDF.

The most commonly reported problems involve rendering the spark graphs (inline
histogram and line chart) on Windows. One common fix is to switch your locale. The
function `fix_windows_histograms()` does this for you.

In order to render the sparkgraphs in html or PDF histogram you may need to
change fonts to one that supports blocks or Braille (depending on which you
need). Please review the separate vignette and associated template for details.
---
title: "Supporting additional objects"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Supporting additional objects}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
## Introduction

The `skim()` function summarizes data types contained within data frames. It
comes with a set of default summary functions for a wide variety of data types,
but this is not comprehensive. Package authors can add support for skimming
their specific data types in their packages, and they can provide different
defaults in their own summary functions.

This example will illustrate this by creating support for the `sf` object
produced by the  "sf: Simple Features for R" package. For any object this
involves two required elements and one optional element.

- experiment with interactive changes
- create methods to `get_skimmers` for different objects within this package
- if needed, define any custom statistics

If you are adding skim support to a package you will also need to add `skimr`
to the list of imports. Note that to run the code in this vignette you will
need to install the `sf` package. We suggest not doing that, and instead
substitute whatever package you are working with.

```{r}
library(skimr)
library(sf)
nc <- st_read(system.file("shape/nc.shp", package = "sf"))
```


```{r}
class(nc)

class(nc$geometry)
```

Unlike the example of having a new type of data in a column of a simple data 
frame in the "Using skimr" vignette, this is a different type of object 
with special attributes.

In this object there is also a column of a class that does not have default
skimmers. By default, skimr falls back to use the sfl for character variables.

```{r}
skim(nc$geometry)
```


## Experiment interactively

`skimr` has an opinionated list of functions for each class (e.g. numeric,
factor)  of data. The core package supports many commonly used classes,
but there are many others. You can investigate these defaults by calling
`get_default_skimmer_names()`.

What if your data type isn't covered by defaults? `skimr` usually falls
back to treating the type as a character, which isn't necessarily helpful. In
this case, you're best off adding your data type with `skim_with()`.

Before we begin, we'll be using the following custom summary statistic
throughout. The function gets the geometry's crs and combines it into a string.

```{r}
get_crs <- function(column) {
  crs <- sf::st_crs(column)

  paste0("epsg: ", crs[["epsg"]], " proj4string: '", crs[["proj4string"]], "'")
}
```

This function, like all summary functions used by `skimr` has two notable
features.

*  It accepts a vector as its single argument
*  It returns a scalar

There are a lot of functions that fulfill these criteria:

* existing functions from base, stats, or other packages,
* lambda's created using the Tidyverse-style syntax
* custom functions that have been defined in the `skimr` package
* custom functions that you have defined.

Not fulfilling the two criteria can lead to some very confusing behavior within
`skimr`. Beware! An example of this issue is the base `quantile()` function in
default `skimr` percentiles are returned by using `quantile()` five 
times.

Next, we create a custom skimming function. To do this, we need to think about
the many specific classes of data in the `sf` package.  From above, you can see 
the geometry column has two classes: 1st the specific geometry type (e.g. 
`sfc_MULTIPOLYGON` `sfc_LINESTRING`, `sfc_POLYGON`, `sfc_MULTIPOINT`) and 2nd 
the general sfc class. Skimr will try to find a sfl() helper function for the
classes in the order they appear in `class(.)` (see S3 classes for more detail 
[*Advanced R*](https://adv-r.hadley.nz/s3.html)). The following example will 
build  support for `sfc`, which encompasses all `sf` objects: `sfc_MULTIPOLYGON` 
`sfc_LINESTRING`, `sfc_POLYGON`, `sfc_MULTIPOINT`. If we want custom skim_with 
functions we can write `sfl()` helper functions for the geometry type. 


```{r}
skim_sf <- skim_with(
  sfc = sfl(
    n_unique = n_unique,
    valid = ~ sum(sf::st_is_valid(.)),
    crs = get_crs
  )
)
```

The example above creates a new *function*, and you can call that function on
a specific column with `sfc` data to get the appropriate summary 
statistics. The `skim_with` factory also uses the default skimrs for things 
like factors, characters, and numerics. Therefore our `skim_sf` is like the regular
`skim` function with the added ability to summarize `sfc` columns.

```{r}
skim_sf(nc$geometry)
```


While this works for any data type and you can also include it within any 
package (assuming your users load skimr), there is an even better approach in 
this case. To take full advantage of `skimr`, we'll dig a bit into its API.

## Adding new methods

`skimr` has a lookup mechanism, based on the function `get_skimmers()`, to
find default summary functions for each class. This is based on the S3 class
system. You can learn more about it in
[*Advanced R*](https://adv-r.hadley.nz/s3.html).

This requires that you add `skimr` to your list of dependencies.

To export a new set of defaults for a data type, create a method for the generic
function `get_skimmers`. Each of those methods returns an `sfl`, a `skimr`
function list. This is the same list-like data structure used in the
`skim_with()` example above. But note! There is one key difference. When adding
a generic we also want to identify the `skim_type` in the `sfl`. You will
probably want to use `skimr::get_skimmers.sfc()` but that will not work in a
vignette.

```{r}
#' @importFrom skimr get_skimmers
#' @export
get_skimmers.sfc <- function(column) {
  sfl(
    skim_type = "sfc",
    n_unique = n_unique,
    valid = ~ sum(sf::st_is_valid(.)),
    crs = get_crs
  )
}
```

The same strategy follows for other data types.

* Create a method
* return an `sfl`
* make sure that the `skim_type` is there

Users of your package should load `skimr` to get the `skim()` function 
(although you could import and reexport it). Once
loaded, a call to `get_default_skimmer_names()` will return defaults for your
data types as well! 

```{r}
get_default_skimmer_names()
```

They will then be able to use `skim()` directly.

```{r}
skim(nc)
```


## Conclusion

This is a very simple example. For a package such as `sf` the custom statistics
will likely  be much more complex. The flexibility of `skimr` allows you to
manage that.

Thanks to Jakub Nowosad, Tiernan Martin, Edzer Pebesma, Michael Sumner, and 
Kyle Butts for inspiring and helping with the development of this code.
---
title: "Skimr defaults"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    latex_engine: xelatex
vignette: >
  %\VignetteIndexEntry{Skimr defaults}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

This vignette simply displays the default data types and summary functions for
`skimr`. Customizing `skimr` is explained in the **Using Skimr** vignette.

## The base skimmers

`skimr` has a group of functions that it applies to all data types. We call
these the "base"" skimmers:

*  `n_missing`: The number of missing values in the column.
*  `complete_rate`: The ratio of non-missing values to the total values in the
   column.

## Default skimmers

To learn more about the functions used in this package, use the function
`get_default_skimmer_names()`.

```{r}
library(skimr)
get_default_skimmer_names()
```

The counterpart to this function is `get_default_skimmers()`, which returns the
functions themselves. If you are interested in a particular class within
`skimr`, pass it as a string to either function.

```{r}
get_default_skimmer_names("numeric")
```

The same information is stored in the `skimmers_used` attribute of the object
returned by `skim()`.
---
title: "Using Fonts"
output: 
  rmarkdown::html_vignette:
    latex_engine: xelatex
vignette: >
  %\VignetteIndexEntry{Using Fonts}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

One of the features of skim is the production of spark graphs for numeric data.
However these graphs may not render properly because of lack of font support 
or for other reasons, such as an operating system that does not support UTF-8. 
In your specific environment this may depend on the fonts 
installed and the operating system and may occur for only 
specific types of documents such as PDF documents.
Skimr supports kable() formatted tables, which is used in conjunction with 
fonts to render the spark graphs.

To produce spark histograms a font that supports block elements must be used. To
support the spark line graphs a font with Braille support must be used. 

Well-known fonts that support block elements include _DejaVu Sans_ and _Arial_. 
Their availability depends on your operating system. 

The yaml at the beginning of a document using custom fonts will generally be
similar to that shown below.

```
---
title: "Untitled"
mainfont: DejaVu Sans
output:
  html_document: default
  pdf_document:
    latex_engine: xelatex
  word_document: default
font-family: Times New Roman
---
```


A further discussion of this (with examples) is available in the 
"Using fonts" template for skimr.  If you are using RStudio you can open that 
template by opening a new markdown file and selecting "From template"
and then choosing it.  Alternatively this file is available inside the 
skimr folder or repository at 
inst/markdown/templates/fonts_in_skimr/skeleton/skeleton.Rmd. 

If you are having difficulties making the spark graphs work, you can opt to turn 
them off using the code below.

```
no_sparks <- skim_with(numeric = sfl(hist = NULL), ts = sfl(line_graph = NULL))

```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skim_print.R
\name{print}
\alias{print}
\alias{print.skim_df}
\alias{print.one_skim_df}
\alias{print.skim_list}
\alias{print.summary_skim_df}
\title{Print \code{skim} objects}
\usage{
\method{print}{skim_df}(
  x,
  include_summary = TRUE,
  n = Inf,
  width = Inf,
  n_extra = NULL,
  strip_metadata = getOption("skimr_strip_metadata", FALSE),
  rule_width = base::options()$width,
  summary_rule_width = 40,
  ...
)

\method{print}{one_skim_df}(
  x,
  n = Inf,
  .width = Inf,
  n_extra = NULL,
  strip_metadata = getOption("skimr_strip_metadata", FALSE),
  .rule_width = base::options()$width,
  ...
)

\method{print}{skim_list}(
  x,
  n = Inf,
  width = Inf,
  n_extra = NULL,
  .rule_width = base::options()$width,
  ...
)

\method{print}{summary_skim_df}(x, .summary_rule_width = 40, ...)
}
\arguments{
\item{x}{Object to format or print.}

\item{include_summary}{Whether a summary of the data frame should be printed}

\item{n}{Number of rows to show. If \code{NULL}, the default, will print all rows
if less than option \code{tibble.print_max}. Otherwise, will print
\code{tibble.print_min} rows.}

\item{width}{Width of text output to generate. This defaults to \code{NULL}, which
means use \code{getOption("tibble.width")} or (if also \code{NULL})
\code{getOption("width")}; the latter displays only the columns that fit on one
screen. You can also set \code{options(tibble.width = Inf)} to override this
default and always print all columns, this may be slow for very wide tibbles.}

\item{n_extra}{Number of extra columns to print abbreviated information for,
if the width is too small for the entire tibble. If \code{NULL}, the default,
will print information about at most \code{tibble.max_extra_cols} extra columns.}

\item{strip_metadata}{Whether tibble metadata should be removed.}

\item{rule_width}{Width of the cli rules in printed skim object. Defaults
to base::options()$width}

\item{summary_rule_width}{Width of Data Summary cli rule, defaults to 40.}

\item{...}{Other arguments passed on to individual methods.}

\item{.width}{Width for the tibble for each type.}

\item{.rule_width}{Width for the rule above the skim results for each type.}

\item{.summary_rule_width}{the width for the main rule above the summary.}
}
\description{
\code{skimr} has custom print methods for all supported objects. Default printing
methods for \code{knitr}/ \code{rmarkdown} documents is also provided.
}
\section{Methods (by class)}{
\itemize{
\item \code{skim_df}: Print a skimmed data frame (\code{skim_df} from \code{\link[=skim]{skim()}}).

\item \code{one_skim_df}: Print an entry within a partitioned \code{skim_df}.

\item \code{skim_list}: Print a \code{skim_list}, a list of \code{skim_df} objects.

\item \code{summary_skim_df}: Print method for a \code{summary_skim_df} object.
}}

\section{Printing options}{


For better or for worse, \code{skimr} often produces more output than can fit in
the standard R console. Fortunately, most modern environments like RStudio
and Jupyter support more than 80 character outputs. Call
\code{options(width = 90)} to get a better experience with \code{skimr}.

The print methods in \code{skimr} wrap those in the \link[tibble:formatting]{tibble}
package. You can control printing behavior using the same global options.
}

\section{Behavior in \code{dplyr} pipelines}{


Printing a \code{skim_df} requires specific columns that might be dropped when
using \code{\link[dplyr:select]{dplyr::select()}} or \code{\link[dplyr:summarise]{dplyr::summarize()}} on a \code{skim_df}. In those
cases, this method falls back to \code{\link[tibble:formatting]{tibble::print.tbl()}}.
}

\section{Controlling metadata behavior}{


On POSIX systems, \code{skimr} removes the tibble metadata when generating output.
On some platforms, this can lead to all output getting removed. To disable
that behavior, set either \code{strip_metadata = FALSE} when calling print or use
\code{options(skimr_strip_metadata = FALSE)}. The \code{crayon} package and the color
support within \code{tibble} is also a factor. If your \code{skimr} results tables are
empty you may need to run the following \code{options(crayon.enabled = FALSE)}.
}

\seealso{
\code{\link[tibble:trunc_mat]{tibble::trunc_mat()}} For a list of global options for customizing
print formatting. \code{\link[crayon:has_color]{crayon::has_color()}} for the variety of issues that
affect tibble's color support.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skim_print.R
\name{repr}
\alias{repr}
\alias{repr_text.skim_df}
\alias{repr_text.skim_list}
\alias{repr_text.one_skim_df}
\title{Skimr printing within Jupyter notebooks}
\usage{
\method{repr_text}{skim_df}(obj, ...)

\method{repr_text}{skim_list}(obj, ...)

\method{repr_text}{one_skim_df}(obj, ...)
}
\arguments{
\item{obj}{The object to \link{print} and then return the output.}

\item{...}{ignored.}
}
\value{
None. \code{invisible(NULL)}.
}
\description{
This reproduces printed results in the console. By default Jupyter kernels
render the final object in the cell. We want the version printed by
\code{skimr} instead of the data that it contains.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{stats}
\alias{stats}
\alias{n_missing}
\alias{n_complete}
\alias{complete_rate}
\alias{n_whitespace}
\alias{sorted_count}
\alias{top_counts}
\alias{inline_hist}
\alias{n_empty}
\alias{min_char}
\alias{max_char}
\alias{n_unique}
\alias{ts_start}
\alias{ts_end}
\alias{inline_linegraph}
\alias{list_lengths_min}
\alias{list_lengths_median}
\alias{list_lengths_max}
\alias{list_min_length}
\alias{list_max_length}
\title{Summary statistic functions}
\usage{
n_missing(x)

n_complete(x)

complete_rate(x)

n_whitespace(x)

sorted_count(x)

top_counts(x, max_char = 3, max_levels = 4)

inline_hist(x, n_bins = 8)

n_empty(x)

min_char(x)

max_char(x)

n_unique(x)

ts_start(x)

ts_end(x)

inline_linegraph(x, length.out = 16)

list_lengths_min(x)

list_lengths_median(x)

list_lengths_max(x)

list_min_length(x)

list_max_length(x)
}
\arguments{
\item{x}{A vector}

\item{max_char}{In \code{top} = 3, max_levels = 4}

\item{max_levels}{The maximum number of levels to be displayed.}

\item{n_bins}{In \code{inline_hist}, the number of histogram bars.}

\item{length.out}{In \code{inline_linegraph}, the length of the character time
series.}
}
\description{
\code{skimr} provides extensions to a variety of functions with R's stats package
to simplify creating summaries of data. All functions are vectorized over the
first argument. Additional arguments should be set in the \code{\link[=sfl]{sfl()}} that sets
the appropriate skimmers for a data type. You can use these, along with other
vectorized R functions, for creating custom sets of summary functions for
a given data type.
}
\section{Functions}{
\itemize{
\item \code{n_missing}: Calculate the sum of \code{NA} and \code{NULL} (i.e. missing) values.

\item \code{n_complete}: Calculate the sum of not \code{NA} and \code{NULL} (i.e. missing)
values.

\item \code{complete_rate}: Calculate complete values; complete values are not missing.

\item \code{n_whitespace}: Calculate the number of rows containing only whitespace
values using s+ regex.

\item \code{sorted_count}: Create a contingency table and arrange its levels in
descending order. In case of ties, the ordering of results is alphabetical
and depends upon the locale. \code{NA} is treated as a ordinary value for
sorting.

\item \code{top_counts}: Compute and collapse a contingency table into a single
character scalar. Wraps \code{\link[=sorted_count]{sorted_count()}}.

\item \code{inline_hist}: Generate inline histogram for numeric variables. The
character length of the histogram is controlled by the formatting options
for character vectors.

\item \code{n_empty}: Calculate the number of blank values in a character vector.
A "blank" is equal to "".

\item \code{min_char}: Calculate the minimum number of characters within a
character vector.

\item \code{max_char}: Calculate the maximum number of characters within a
character vector.

\item \code{n_unique}: Calculate the number of unique elements but remove \code{NA}.

\item \code{ts_start}: Get the start for a time series without the frequency.

\item \code{ts_end}: Get the finish for a time series without the frequency.

\item \code{inline_linegraph}: Generate inline line graph for time series variables. The
character length of the line graph is controlled by the formatting options
for character vectors.
Based on the function in the pillar package.

\item \code{list_lengths_min}: Get the length of the shortest list in a vector of lists.

\item \code{list_lengths_median}: Get the median length of the lists.

\item \code{list_lengths_max}: Get the maximum length of the lists.

\item \code{list_min_length}: Get the length of the shortest list in a vector of lists.

\item \code{list_max_length}: Get the length of the longest list in a vector of lists.
}}

\seealso{
\code{\link[=get_skimmers]{get_skimmers()}} for customizing the functions called by \code{\link[=skim]{skim()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skim_obj.R
\name{skim-obj}
\alias{skim-obj}
\alias{has_type_column}
\alias{has_variable_column}
\alias{has_skimr_attributes}
\alias{has_skim_type_attribute}
\alias{is_data_frame}
\alias{is_skim_df}
\alias{is_one_skim_df}
\alias{is_skim_list}
\alias{could_be_skim_df}
\alias{assert_is_skim_df}
\alias{assert_is_skim_list}
\alias{assert_is_one_skim_df}
\title{Test if an object is compatible with \code{skimr}}
\usage{
has_type_column(object)

has_variable_column(object)

has_skimr_attributes(object)

has_skim_type_attribute(object)

is_data_frame(object)

is_skim_df(object)

is_one_skim_df(object)

is_skim_list(object)

could_be_skim_df(object)

assert_is_skim_df(object)

assert_is_skim_list(object)

assert_is_one_skim_df(object)
}
\arguments{
\item{object}{Any \code{R} object.}
}
\description{
Objects within \code{skimr} are identified by a class, but they require additional
attributes and data columns for all operations to succeed. These checks help
ensure this. While they have some application externally, they are mostly
used internally.
}
\details{
Most notably, a \code{skim_df} has columns \code{skim_type} and \code{skim_variable}. And
has the following special attributes
\itemize{
\item \code{data_rows}: n rows in the original data
\item \code{data_cols}: original number of columns
\item \code{df_name}: name of the original data frame
\item \code{dt_key}: name of the key if original is a data.table
\item \code{groups}: if there were group variables
\item \code{base_skimmers}: names of functions applied to all skim types
\item \code{skimmers_used}: names of functions used to skim each type
}

The functions in these checks work like \code{\link[=all.equal]{all.equal()}}. The return \code{TRUE} if
the check passes, or otherwise notifies why the check failed. This makes them
more useful when throwing errors.
}
\section{Functions}{
\itemize{
\item \code{has_type_column}: Does the object have the \code{skim_type} column?

\item \code{has_variable_column}: Does the object have the \code{skim_variable} column?

\item \code{has_skimr_attributes}: Does the object have the appropriate \code{skimr} attributes?

\item \code{has_skim_type_attribute}: Does the object have a \code{skim_type} attribute? This makes
it a \code{one_skim_df}.

\item \code{is_data_frame}: Is the object a data frame?

\item \code{is_skim_df}: Is the object a \code{skim_df}?

\item \code{is_one_skim_df}: Is the object a \code{one_skim_df}? This is similar to a
\code{skim_df}, but does not have the \code{type} column. That is stored as an
attribute instead.

\item \code{is_skim_list}: Is the object a \code{skim_list}?

\item \code{could_be_skim_df}: Is this a data frame with \code{skim_variable} and
\code{skim_type} columns?

\item \code{assert_is_skim_df}: Stop if the object is not a \code{skim_df}.

\item \code{assert_is_skim_list}: Stop if the object is not a \code{skim_list}.

\item \code{assert_is_one_skim_df}: Stop if the object is not a \code{one_skim_df}.
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{deprecated-v1}
\alias{deprecated-v1}
\alias{skim_to_wide}
\alias{skim_to_list}
\alias{skim_format}
\title{Deprecated functions from skimr v1}
\usage{
skim_to_wide(.data, ...)

skim_to_list(.data, ...)

skim_format(...)
}
\arguments{
\item{.data}{A tibble, or an object that can be coerced into a tibble.}

\item{...}{Columns to select for skimming. When none are provided, the
default is to skim all columns.}
}
\value{
Either A \code{skim_df} or a \code{skim_list} object.
}
\description{
Skimr used to offer functions that combined skimming with a secondary effect,
like reshaping the data, building a list or printing the results. Some of
these behaviors are no longer necessary. \code{\link[=skim]{skim()}} always returns a wide
data frame. Others have been replaced by functions that do a single thing.
\code{\link[=partition]{partition()}} creates a list-like object from a skimmed data frame.
}
\section{Functions}{
\itemize{
\item \code{skim_to_wide}: \code{\link[=skim]{skim()}} always produces a wide data frame.

\item \code{skim_to_list}: \code{\link[=partition]{partition()}} creates a list.

\item \code{skim_format}: \code{\link[=print]{print()}} and \code{\link[=skim_with]{skim_with()}} set options.
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reshape.R
\name{focus}
\alias{focus}
\title{Only show a subset of summary statistics after skimming}
\usage{
focus(.data, ...)
}
\arguments{
\item{.data}{A \code{skim_df} object.}

\item{...}{One or more unquoted expressions separated by commas. Variable
names can be used as if they were positions in the data frame, so
expressions like x:y can be used to select a range of variables.}
}
\description{
This function is a variant of \code{\link[dplyr:select]{dplyr::select()}} designed to work with
\code{skim_df} objects. When using \code{focus()}, \code{skimr} metadata columns are kept,
and \code{skimr} print methods are still utilized. Otherwise, the signature and
behavior is identical to \code{\link[dplyr:select]{dplyr::select()}}.
}
\examples{
# Compare
iris \%>\%
  skim() \%>\%
  dplyr::select(n_missing)

iris \%>\%
  skim() \%>\%
  focus(n_missing)

# This is equivalent to
iris \%>\%
  skim() \%>\%
  dplyr::select(skim_variable, skim_type, n_missing)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skimr-package.R
\docType{package}
\name{skimr-package}
\alias{skimr-package}
\alias{skimr}
\title{Skim a data frame}
\description{
This package provides an alternative to the default summary functions
within R. The package's API is tidy, functions take data frames, return
data frames and can work as part of a pipeline. The returned \code{skimr}
object is subsettable and offers a human readable output.
}
\details{
\code{skimr} is opinionated, providing a strong set of summary statistics
that are generated for a variety of different data types. It is also
provides an API for customization. Users can change both the functions
dispatched and the way the results are formatted.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reshape.R
\name{to_long}
\alias{to_long}
\alias{to_long.default}
\alias{to_long.skim_df}
\title{Create "long" skim output}
\usage{
to_long(.data, ..., skim_fun = skim)

\method{to_long}{default}(.data, ..., skim_fun = skim)

\method{to_long}{skim_df}(.data, ..., skim_fun = skim)
}
\arguments{
\item{.data}{A data frame or an object that can be coerced into a data frame.}

\item{...}{Columns to select for skimming. When none are provided, the
default is to skim all columns.}

\item{skim_fun}{The skim function used.}
}
\value{
A tibble
}
\description{
Skim results returned as a tidy long data frame with four columns:
variable, type, stat and formatted.
}
\section{Methods (by class)}{
\itemize{
\item \code{default}: Skim a data frame and convert the results to a
long data frame.

\item \code{skim_df}: Transform a skim_df to a long data frame.
}}

\examples{
to_long(iris)
to_long(skim(iris))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_skimmers.R
\name{get_skimmers}
\alias{get_skimmers}
\alias{get_skimmers.default}
\alias{get_skimmers.numeric}
\alias{get_skimmers.factor}
\alias{get_skimmers.character}
\alias{get_skimmers.logical}
\alias{get_skimmers.complex}
\alias{get_skimmers.Date}
\alias{get_skimmers.POSIXct}
\alias{get_skimmers.difftime}
\alias{get_skimmers.Timespan}
\alias{get_skimmers.ts}
\alias{get_skimmers.list}
\alias{get_skimmers.AsIs}
\alias{modify_default_skimmers}
\title{Retrieve the summary functions for a specific data type}
\usage{
get_skimmers(column)

\method{get_skimmers}{default}(column)

\method{get_skimmers}{numeric}(column)

\method{get_skimmers}{factor}(column)

\method{get_skimmers}{character}(column)

\method{get_skimmers}{logical}(column)

\method{get_skimmers}{complex}(column)

\method{get_skimmers}{Date}(column)

\method{get_skimmers}{POSIXct}(column)

\method{get_skimmers}{difftime}(column)

\method{get_skimmers}{Timespan}(column)

\method{get_skimmers}{ts}(column)

\method{get_skimmers}{list}(column)

\method{get_skimmers}{AsIs}(column)

modify_default_skimmers(skim_type, new_skim_type = NULL, new_funs = list())
}
\arguments{
\item{column}{An atomic vector or list. A column from a data frame.}

\item{skim_type}{A character scalar. The class of the object with default
skimmers.}

\item{new_skim_type}{The type to assign to the looked up set of skimmers.}

\item{new_funs}{Replacement functions for those in}
}
\value{
A \code{skim_function_list} object.
}
\description{
These functions are used to set the default skimming functions for a data
type. They are combined with the base skim function list (\code{sfl}) in
\code{\link[=skim_with]{skim_with()}}, to create the summary tibble for each type.
}
\details{
When creating your own set of skimming functions, call \code{\link[=sfl]{sfl()}} within a
\code{\link[=get_skimmers]{get_skimmers()}} method for your particular type. Your call to \code{\link[=sfl]{sfl()}} should
also provide a matching class in the \code{skim_type} argument.  Otherwise, it
will not be possible to dynamically reassign your default functions when
working interactively.

Call \code{\link[=get_default_skimmers]{get_default_skimmers()}} to see the functions for each type of summary
function currently supported. Call \code{\link[=get_default_skimmer_names]{get_default_skimmer_names()}} to just see
the names of these functions. Use \code{\link[=modify_default_skimmers]{modify_default_skimmers()}} for a method
for changing the \code{skim_type} or functions for a default \code{sfl}. This is useful
for creating new default \code{sfl}'s.
}
\section{Methods (by class)}{
\itemize{
\item \code{default}: The default method for skimming data. Only used when
a column's data type doesn't match currently installed types. Call
\link{get_default_skimmer_names} to see these defaults.

\item \code{numeric}: Summary functions for numeric columns, covering both
\code{\link[=double]{double()}} and \code{\link[=integer]{integer()}} classes: \code{\link[=mean]{mean()}}, \code{\link[=sd]{sd()}}, \code{\link[=quantile]{quantile()}} and
\code{\link[=inline_hist]{inline_hist()}}.

\item \code{factor}: Summary functions for factor columns:
\code{\link[=is.ordered]{is.ordered()}}, \code{\link[=n_unique]{n_unique()}} and \code{\link[=top_counts]{top_counts()}}.

\item \code{character}: Summary functions for character columns. Also, the
default for unknown columns: \code{\link[=min_char]{min_char()}}, \code{\link[=max_char]{max_char()}}, \code{\link[=n_empty]{n_empty()}},
\code{\link[=n_unique]{n_unique()}} and \code{\link[=n_whitespace]{n_whitespace()}}.

\item \code{logical}: Summary functions for logical/ boolean columns:
\code{\link[=mean]{mean()}}, which produces rates for each value, and \code{\link[=top_counts]{top_counts()}}.

\item \code{complex}: Summary functions for complex columns: \code{\link[=mean]{mean()}}.

\item \code{Date}: Summary functions for \code{Date} columns: \code{\link[=min]{min()}},
\code{\link[=max]{max()}}, \code{\link[=median]{median()}} and \code{\link[=n_unique]{n_unique()}}.

\item \code{POSIXct}: Summary functions for \code{POSIXct} columns: \code{\link[=min]{min()}},
\code{\link[=max]{max()}}, \code{\link[=median]{median()}} and \code{\link[=n_unique]{n_unique()}}.

\item \code{difftime}: Summary functions for \code{difftime} columns: \code{\link[=min]{min()}},
\code{\link[=max]{max()}}, \code{\link[=median]{median()}} and \code{\link[=n_unique]{n_unique()}}.

\item \code{Timespan}: Summary functions for \code{Timespan} columns: \code{\link[=min]{min()}},
\code{\link[=max]{max()}}, \code{\link[=median]{median()}} and \code{\link[=n_unique]{n_unique()}}.

\item \code{ts}: Summary functions for \code{ts} columns: \code{\link[=min]{min()}},
\code{\link[=max]{max()}}, \code{\link[=median]{median()}} and \code{\link[=n_unique]{n_unique()}}.

\item \code{list}: Summary functions for \code{list} columns: \code{\link[=n_unique]{n_unique()}},
\code{\link[=list_min_length]{list_min_length()}} and \code{\link[=list_max_length]{list_max_length()}}.

\item \code{AsIs}: Summary functions for \code{AsIs} columns: \code{\link[=n_unique]{n_unique()}},
\code{\link[=list_min_length]{list_min_length()}} and \code{\link[=list_max_length]{list_max_length()}}.
}}

\examples{
# Defining default skimming functions for a new class, `my_class`.
# Note that the class argument is required for dynamic reassignment.
get_skimmers.my_class <- function(column) {
  sfl(
    skim_type = "my_class",
    mean,
    sd
  )
}

# Integer and double columns are both "numeric" and are treated the same
# by default. To switch this behavior in another package, add a method.
get_skimmers.integer <- function(column) {
  sfl(
    skim_type = "integer",
    p50 = ~ stats::quantile(
      .,
      probs = .50, na.rm = TRUE, names = FALSE, type = 1
    )
  )
}
x <- mtcars[c("gear", "carb")]
class(x$carb) <- "integer"
skim(x)
\dontrun{
# In a package, to revert to the V1 behavior of skimming separately with the
# same functions, assign the numeric `get_skimmers`.
get_skimmers.integer <- skimr::get_skimmers.numeric

# Or, in a local session, use `skim_with` to create a different `skim`.
new_skim <- skim_with(integer = skimr::get_skimmers.numeric())

# To apply a set of skimmers from an old type to a new type
get_skimmers.new_type <- function(column) {
  modify_default_skimmers("old_type", new_skim_type = "new_type")
}
}
}
\seealso{
\code{\link[=sfl]{sfl()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fix_windows_histograms}
\alias{fix_windows_histograms}
\title{Fix unicode histograms on Windows}
\usage{
fix_windows_histograms()
}
\description{
This functions changes your session's locale to address issues with printing
histograms on Windows.
}
\details{
There are known issues with printing the spark-histogram characters when
printing a data frame, appearing like this: "<U+2582><U+2585><U+2587>".
This longstanding problem originates in the low-level code for printing
dataframes.
}
\seealso{
\code{\link[=skim_without_charts]{skim_without_charts()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sfl.R
\name{sfl}
\alias{sfl}
\title{Create a skimr function list}
\usage{
sfl(..., skim_type = "")
}
\arguments{
\item{...}{Inherited from dplyr::data_masking() for dplyr version 1 or later
or dplyr::funs() for older versions of dplyr.
A list of functions
specified by:
\itemize{
\item Their name, \code{"mean"}
\item The function itself, \code{mean}
\item A call to the function with \code{.} as a dummy argument,
\code{mean(., na.rm = TRUE)}
\item An anonymous function in \pkg{purrr} notation, \code{~mean(., na.rm = TRUE)}
}}

\item{skim_type}{A character scalar. This is used to match locally-provided
skimmers with defaults. See \code{\link[=get_skimmers]{get_skimmers()}} for more detail.}
}
\value{
A \code{skimr_function_list}, which contains a list of \code{fun_calls},
returned by \code{dplyr::funs()} and a list of skimming functions to drop.
}
\description{
This constructor is used to create a named list of functions. It also you
also pass \code{NULL} to identify a skimming function that you wish to remove.
Only functions that return a single value, working with \code{\link[dplyr:summarise]{dplyr::summarize()}},
can be used within \code{sfl}.
}
\details{
\code{sfl()} will automatically generate callables and names for a variety of
inputs, including functions, formulas and strings. Nonetheless, we recommend
providing names when reasonable to get better \code{\link[=skim]{skim()}} output.
}
\examples{
# sfl's can take a variety of input formats and will generate names
# if not provided.
sfl(mad, "var", ~ length(.)^2)

# But these can generate unpredictable names in your output.
# Better to set your own names.
sfl(mad = mad, variance = "var", length_sq = ~ length(.)^2)

# sfl's can remove individual skimmers from defaults by passing NULL.
sfl(hist = NULL)

# When working interactively, you don't need to set a type.
# But you should when defining new defaults with `get_skimmers()`.
get_skimmers.my_new_class <- function(column) {
  sfl(n_missing, skim_type = "my_new_class")
}
}
\seealso{
\code{\link[dplyr:funs]{dplyr::funs()}}, \code{\link[=skim_with]{skim_with()}} and \code{\link[=get_skimmers]{get_skimmers()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr.R
\name{mutate.skim_df}
\alias{mutate.skim_df}
\title{Mutate a skim_df}
\usage{
\method{mutate}{skim_df}(.data, ...)
}
\arguments{
\item{.data}{A \code{skim_df}, which behaves like a \code{tbl.}}

\item{...}{Name-value pairs of expressions, each with length 1 or the same
length as the number of rows in the group, if using \code{\link[dplyr:group_by]{dplyr::group_by()}}, or
in the entire input (if not using groups). The name of each argument will
be the name of a new variable, and the value will be its corresponding
value. Use \code{NULL} value in \code{\link[dplyr:mutate]{dplyr::mutate()}} to drop a variable. New
variables overwrite existing variables of the same name.

The arguments in \code{...} are automatically quoted with \code{\link[rlang:nse-defuse]{rlang::quo()}} and
evaluated with \code{\link[rlang:eval_tidy]{rlang::eval_tidy()}} in the context of the data frame. They
support unquoting \link[rlang:nse-force]{rlang::quasiquotation} and splicing. See
\code{vignette("programming", package = "dplyr")} for an introduction to these
concepts.}
}
\value{
A \code{skim_df} object, which also inherits the class(es) of the input
data. In many ways, the object behaves like a \code{\link[tibble:tibble]{tibble::tibble()}}.
}
\description{
\code{\link[dplyr:mutate]{dplyr::mutate()}} currently drops attributes, but we need to keep them around
for other skim behaviors. Otherwise the behavior is exactly the same. For
more information, see \url{https://github.com/tidyverse/dplyr/issues/3429}.
}
\seealso{
\code{\link[dplyr:mutate]{dplyr::mutate()}} for the function's expected behavior.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skimr-package.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\alias{contains}
\alias{ends_with}
\alias{everything}
\alias{matches}
\alias{num_range}
\alias{one_of}
\alias{starts_with}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}

  \item{tidyselect}{\code{\link[tidyselect:starts_with]{contains}}, \code{\link[tidyselect:starts_with]{ends_with}}, \code{\link[tidyselect]{everything}}, \code{\link[tidyselect:starts_with]{matches}}, \code{\link[tidyselect:starts_with]{num_range}}, \code{\link[tidyselect]{one_of}}, \code{\link[tidyselect]{starts_with}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skim_with.R
\name{skim_with}
\alias{skim_with}
\title{Set or add the summary functions for a particular type of data}
\usage{
skim_with(
  ...,
  base = sfl(n_missing = n_missing, complete_rate = complete_rate),
  append = TRUE
)
}
\arguments{
\item{...}{One or more (\code{sfl}) \code{skimmer_function_list} objects, with an
argument name that matches a particular data type.}

\item{base}{An \code{sfl} that sets skimmers for all column types.}

\item{append}{Whether the provided options should be in addition to the
defaults already in \code{skim}. Default is \code{TRUE}.}
}
\value{
A new \code{skim()} function. This is callable. See \code{\link[=skim]{skim()}} for more
details.
}
\description{
While skim is designed around having an opinionated set of defaults, you
can use this function to change the summary statistics that it returns.
}
\details{
\code{skim_with()} is a closure: a function that returns a new function. This
lets you have several skimming functions in a single R session, but it
also means that you need to assign the return of \code{skim_with()} before
you can use it.

You assign values within \code{skim_with} by using the \code{\link[=sfl]{sfl()}} helper (\code{skimr}
function list). This helper behaves mostly like \code{\link[dplyr:funs]{dplyr::funs()}}, but lets
you also identify which skimming functions you want to remove, by setting
them to \code{NULL}. Assign an \code{sfl} to each column type that you wish to modify.

Functions that summarize all data types, and always return the same type
of value, can be assigned to the \code{base} argument. The default base skimmers
compute the number of missing values \code{\link[=n_missing]{n_missing()}} and the rate of values
being complete, i.e. not missing, \code{\link[=complete_rate]{complete_rate()}}.

When \code{append = TRUE} and local skimmers have names matching the names of
entries in the default \code{skim_function_list}, the values in the default list
are overwritten. Similarly, if \code{NULL} values are passed within \code{sfl()}, these
default skimmers are dropped. Otherwise, if \code{append = FALSE}, only the
locally-provided skimming functions are used.

Note that \code{append} only applies to the \code{typed} skimmers (i.e. non-base).
See \code{\link[=get_default_skimmer_names]{get_default_skimmer_names()}} for a list of defaults.
}
\examples{
# Use new functions for numeric functions. If you don't provide a name,
# one will be automatically generated.
my_skim <- skim_with(numeric = sfl(median, mad), append = FALSE)
my_skim(faithful)

# If you want to remove a particular skimmer, set it to NULL
# This removes the inline histogram
my_skim <- skim_with(numeric = sfl(hist = NULL))
my_skim(faithful)

# This works with multiple skimmers. Just match names to overwrite
my_skim <- skim_with(numeric = sfl(iqr = IQR, p25 = NULL, p75 = NULL))
my_skim(faithful)

# Alternatively, set `append = FALSE` to replace the skimmers of a type.
my_skim <- skim_with(numeric = sfl(mean = mean, sd = sd), append = FALSE)

# Skimmers are unary functions. Partially apply arguments during assigment.
# For example, you might want to remove NA values.
my_skim <- skim_with(numeric = sfl(iqr = ~ IQR(., na.rm = TRUE)))

# Set multiple types of skimmers simultaneously.
my_skim <- skim_with(numeric = sfl(mean), character = sfl(length))

# Or pass the same as a list, unquoting the input.
my_skimmers <- list(numeric = sfl(mean), character = sfl(length))
my_skim <- skim_with(!!!my_skimmers)

# Use the v1 base skimmers instead.
my_skim <- skim_with(base = sfl(
  missing = n_missing,
  complete = n_complete,
  n = length
))

# Remove the base skimmers entirely
my_skim <- skim_with(base = NULL)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vctrs.R
\name{skimr-vctrs}
\alias{skimr-vctrs}
\alias{vec_ptype2.skim_df.skim_df}
\alias{vec_ptype2.skim_df.tbl_df}
\alias{vec_ptype2.tbl_df.skim_df}
\alias{vec_cast.skim_df.skim_df}
\alias{vec_cast.skim_df.tbl_df}
\alias{vec_cast.tbl_df.skim_df}
\title{Functions for working with the vctrs package}
\usage{
\method{vec_ptype2}{skim_df.skim_df}(x, y, ...)

\method{vec_ptype2}{skim_df.tbl_df}(x, y, ...)

\method{vec_ptype2}{tbl_df.skim_df}(x, y, ...)

\method{vec_cast}{skim_df.skim_df}(x, to, ...)

\method{vec_cast}{skim_df.tbl_df}(x, to, ...)

\method{vec_cast}{tbl_df.skim_df}(x, to, ...)
}
\description{
These make it clear that we need to use the tibble behavior when joining,
concatenating or casting \code{skim_df} objects. For a better discussion, on
why this is important and how these functions work, see:
\url{https://vctrs.r-lib.org/reference/howto-faq-coercion-data-frame.html}.
}
\details{
\verb{vec_ptype2.*} handles finding common prototypes between \code{skim_df} and
similar objects. \verb{vec_cast.*} handles casting between objects. Note that
as of \verb{dplyr 1.0.2}, \code{\link[dplyr:bind]{dplyr::bind_rows()}} does not full support combining
attributes and \code{\link[vctrs:vec_bind]{vctrs::vec_rbind()}} is preferred instead.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skim_print.R
\name{knit_print}
\alias{knit_print}
\alias{knit_print.skim_df}
\alias{knit_print.skim_list}
\alias{knit_print.one_skim_df}
\alias{knit_print.summary_skim_df}
\title{Provide a default printing method for knitr.}
\usage{
\method{knit_print}{skim_df}(x, options = NULL, ...)

\method{knit_print}{skim_list}(x, options = NULL, ...)

\method{knit_print}{one_skim_df}(x, options = NULL, ...)

\method{knit_print}{summary_skim_df}(x, options = NULL, ...)
}
\arguments{
\item{x}{An R object to be printed}

\item{options}{Options passed into the print function.}

\item{...}{Additional arguments passed to the S3 method. Currently ignored,
except two optional arguments \code{options} and \code{inline}; see
the references below.}
}
\value{
A \code{knit_asis} object. Which is used by \code{knitr} when rendered.
}
\description{
Instead of standard R output, \code{knitr} and \code{RMarkdown} documents will have
formatted \code{\link[knitr:kable]{knitr::kable()}} output on return. You can disable this by setting
the chunk option \code{render = normal_print}.
}
\details{
The summary statistics for the original data frame can be disabled by setting
the \code{knitr} chunk option \code{skimr_include_summary = FALSE}. See
\link[knitr:opts_chunk]{knitr::opts_chunk} for more information. You can change the number of digits
shown in the printed table with the \code{skimr_digits} chunk option.

Alternatively, you can call \code{\link[=collapse]{collapse()}} or \code{\link[=yank]{yank()}} to get the particular
\code{skim_df} objects and format them however you like. One warning though.
Because histograms contain unicode characters, they can have unexpected
print results, as R as varying levels of unicode support. This affects
Windows users most commonly. Call \code{vignette("Using_fonts")} for more details.
}
\section{Methods (by class)}{
\itemize{
\item \code{skim_df}: Default \code{knitr} print for \code{skim_df} objects.

\item \code{skim_list}: Default \code{knitr} print for a \code{skim_list}.

\item \code{one_skim_df}: Default \code{knitr} print within a partitioned \code{skim_df}.

\item \code{summary_skim_df}: Default \code{knitr} print for \code{skim_df} summaries.
}}

\seealso{
\code{\link[knitr:kable]{knitr::kable()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skim.R
\name{skim}
\alias{skim}
\alias{skim_tee}
\alias{skim_without_charts}
\title{Skim a data frame, getting useful summary statistics}
\usage{
skim(data, ..., .data_name = NULL)

skim_tee(data, ..., skim_fun = skim)

skim_without_charts(data, ..., .data_name = NULL)
}
\arguments{
\item{data}{A tibble, or an object that can be coerced into a tibble.}

\item{...}{Columns to select for skimming. When none are provided, the
default is to skim all columns.}

\item{.data_name}{The name to use for the data. Defaults to the same as data.}

\item{skim_fun}{The skim function used.}

\item{skim}{The skimming function to use in \code{skim_tee()}.}
}
\value{
A \code{skim_df} object, which also inherits the class(es) of the input
data. In many ways, the object behaves like a \code{\link[tibble:tibble]{tibble::tibble()}}.
}
\description{
\code{skim()} is an alternative to \code{\link[=summary]{summary()}}, quickly providing a broad
overview of a data frame. It handles data of all types, dispatching a
different set of summary functions based on the types of columns in the data
frame.
}
\details{
Each call produces a \code{skim_df}, which is a fundamentally a tibble with a
special print method. One unusual feature of this data frame is pseudo-
namespace for columns. \code{skim()} computes statistics by data type, and it
stores them in the data frame as \verb{<type>.<statistic>}. These types are
stripped when printing the results. The "base" skimmers (\code{n_missing} and
\code{complete_rate}) are the only columns that don't follow this behavior.
See \code{\link[=skim_with]{skim_with()}} for more details on customizing \code{skim()} and
\code{\link[=get_default_skimmers]{get_default_skimmers()}} for a list of default functions.

If you just want to see the printed output, call \code{skim_tee()} instead.
This function returns the original data. \code{skim_tee()} uses the default
\code{skim()}, but you can replace it with the \code{skim} argument.

The data frame produced by \code{skim} is wide and sparse. To avoid type coercion
\code{skimr} uses a type namespace for all summary statistics. Columns for numeric
summary statistics all begin \code{numeric}; for factor summary statistics
begin \code{factor}; and so on.

See \code{\link[=partition]{partition()}} and \code{\link[=yank]{yank()}} for methods for transforming this wide data
frame. The first function splits it into a list, with each entry
corresponding to a data type. The latter pulls a single subtable for a
particular type from the \code{skim_df}.

\code{skim()} is designed to operate in pipes and to generally play nicely with
other \code{tidyverse} functions. This means that you can use \code{tidyselect} helpers
within \code{skim} to select or drop specific columns for summary. You can also
further work with a \code{skim_df} using \code{dplyr} functions in a pipeline.
}
\section{Customizing skim}{

\code{skim()} is an intentionally simple function, with minimal arguments like
\code{\link[=summary]{summary()}}. Nonetheless, this package provides two broad approaches to
how you can customize \code{skim()}'s behavior. You can customize the functions
that are called to produce summary statistics with \code{\link[=skim_with]{skim_with()}}.
}

\section{Unicode rendering}{

If the rendered examples show unencoded values such as \verb{<U+2587>} you will
need to change your locale to allow proper rendering. Please review the
\emph{Using Skimr} vignette for more information
(\code{vignette("Using_skimr", package = "skimr")}).

Otherwise, we export \code{skim_without_charts()} to produce summaries without the
spark graphs. These are the source of the unicode dependency.
}

\examples{
skim(iris)

# Use tidyselect
skim(iris, Species)
skim(iris, starts_with("Sepal"))
skim(iris, where(is.numeric))

# Skim also works groupwise
iris \%>\%
  dplyr::group_by(Species) \%>\%
  skim()

# Which five numeric columns have the greatest mean value?
# Look in the `numeric.mean` column.
iris \%>\%
  skim() \%>\%
  dplyr::select(numeric.mean) \%>\%
  dplyr::top_n(5)

# Which of my columns have missing values? Use the base skimmer n_missing.
iris \%>\%
  skim() \%>\%
  dplyr::filter(n_missing > 0)

# Use skim_tee to view the skim results and
# continue using the original data.
chickwts \%>\%
  skim_tee() \%>\%
  dplyr::filter(feed == "sunflower")

# Produce a summary without spark graphs
iris \%>\%
  skim_without_charts()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skim_obj.R
\name{skim-attr}
\alias{skim-attr}
\alias{data_rows}
\alias{data_cols}
\alias{df_name}
\alias{dt_key}
\alias{group_names}
\alias{base_skimmers}
\alias{skimmers_used}
\title{Functions for accessing skim_df attributes}
\usage{
data_rows(object)

data_cols(object)

df_name(object)

dt_key(object)

group_names(object)

base_skimmers(object)

skimmers_used(object)
}
\arguments{
\item{object}{A \code{skim_df} or \code{skim_list}.}
}
\value{
Data contained within the requested \code{skimr} attribute.
}
\description{
These functions simplify access to attributes contained within a \code{skim_df}.
While all attributes are read-only, being able to extract this information
is useful for different analyses. These functions should always be preferred
over calling base R's attribute functions.
}
\section{Functions}{
\itemize{
\item \code{data_rows}: Get the number of rows in the skimmed data frame.

\item \code{data_cols}: Get the number of columns in the skimmed data frame.

\item \code{df_name}: Get the name of the skimmed data frame. This is only
available in contexts where the name can be looked up. This is often not
the case within a pipeline.

\item \code{dt_key}: Get the key of the skimmed data.table. This is only
available in contexts where \code{data} is of class \code{data.table}.

\item \code{group_names}: Get the names of the groups in the original data frame.
Only available if the data was grouped. Otherwise, \code{NULL}.

\item \code{base_skimmers}: Get the names of the base skimming functions used.

\item \code{skimmers_used}: Get the names of the skimming functions used, separated
by data type.
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_skimmers.R
\name{get_default_skimmers}
\alias{get_default_skimmers}
\alias{get_one_default_skimmer}
\alias{get_default_skimmer_names}
\alias{get_one_default_skimmer_names}
\alias{get_sfl}
\title{View default skimmer names and functions}
\usage{
get_default_skimmers(skim_type = NULL)

get_one_default_skimmer(skim_type)

get_default_skimmer_names(skim_type = NULL)

get_one_default_skimmer_names(skim_type)

get_sfl(skim_type)
}
\arguments{
\item{skim_type}{The class of the column being skimmed.}
}
\description{
These utility functions look up the currently-available defaults for one or
more \code{skim_type}'s. They work with all defaults in the \code{skimr} package, as
well as the defaults in any package that extends \code{skimr}. See
\code{\link[=get_skimmers]{get_skimmers()}} for writing lookup methods for different.
}
\details{
The functions differ in return type and whether or not the result is in
a list. \code{\link[=get_default_skimmers]{get_default_skimmers()}} and \code{\link[=get_one_default_skimmer]{get_one_default_skimmer()}} return
functions. The former returns functions within a typed list, i.e.
\code{list(numeric = list(...functions...))}.

The functions differ in return type and whether or not the result is in
a list. \code{\link[=get_default_skimmer_names]{get_default_skimmer_names()}} and \code{\link[=get_one_default_skimmer_names]{get_one_default_skimmer_names()}}
return a list of character vectors or a single character vector.

\code{\link[=get_sfl]{get_sfl()}} returns the skimmer function list for a particular \code{skim_type}.
It differs from \code{\link[=get_default_skimmers]{get_default_skimmers()}} in that the returned \code{sfl} contains
a list of functions and a \code{skim_type}.
}
\section{Functions}{
\itemize{
\item \code{get_one_default_skimmer}: Get the functions associated with one
\code{skim_type}.

\item \code{get_default_skimmer_names}: Get the names of the functions used in one
or more \code{skim_type}'s.

\item \code{get_one_default_skimmer_names}: Get the names of the functions used in one
\code{skim_type}.

\item \code{get_sfl}: Get the \code{sfl} for a \code{skim_type}.
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.skim_df}
\alias{summary.skim_df}
\title{Summary function for skim_df}
\usage{
\method{summary}{skim_df}(object, ...)
}
\arguments{
\item{object}{a skim dataframe.}

\item{...}{Additional arguments affecting the summary produced. Not used.}
}
\value{
A summary of the skim data frame.
}
\description{
This is a method of the generic function \code{\link[=summary]{summary()}}.
}
\examples{
a <- skim(mtcars)
summary(a)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reshape.R
\name{partition}
\alias{partition}
\alias{bind}
\alias{yank}
\title{Separate a big \code{skim_df} into smaller data frames, by type.}
\usage{
partition(data)

bind(data)

yank(data, skim_type)
}
\arguments{
\item{data}{A \code{skim_df}.}

\item{skim_type}{A character scalar. The subtable to extract from a
\code{skim_df}.}
}
\value{
A \code{skim_list} of \code{skim_df}'s, by type.
}
\description{
The data frames produced by \code{\link[=skim]{skim()}} are wide and sparse, filled with
columns that are mostly \code{NA}. For that reason, it can be convenient to
work with "by type" subsets of the original data frame. These smaller
subsets have their \code{NA} columns removed.
}
\details{
\code{partition()} creates a list of smaller \code{skim_df} data frames. Each entry
in the list is a data type from the original \code{skim_df}. The inverse of
\code{partition()} is \code{bind()}, which takes the list and produces the original
\code{skim_df}. While \code{partition()} keeps all of the subtables as list entries,
\code{yank()} gives you a single subtable for a data type.
}
\section{Functions}{
\itemize{
\item \code{bind}: The inverse of a \code{partition()}. Rebuild the original
\code{skim_df}.

\item \code{yank}: Extract a subtable from a \code{skim_df} with a particular
type.
}}

\examples{
# Create a wide skimmed data frame (a skim_df)
skimmed <- skim(iris)

# Separate into a list of subtables by type
separate <- partition(skimmed)

# Put back together
identical(bind(separate), skimmed)
# > TRUE

# Alternatively, get the subtable of a particular type
yank(skimmed, "factor")
}
