assertr
===

![assertr logo](http://statethatiamin.com/media/assertrlogo.png)

[![Build Status](http://travis-ci.org/ropensci/assertr.svg?branch=master)](https://travis-ci.org/ropensci/assertr)
[![](http://www.r-pkg.org/badges/version/assertr)](https://cran.r-project.org/package=assertr)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/assertr)](https://cran.r-project.org/package=assertr)

### What is it?
The assertr package supplies a suite of functions designed to verify
assumptions about data early in an analysis pipeline so that
data errors are spotted early and can be addressed quickly.

This package does not need to be used with the magrittr/dplyr piping
mechanism but the examples in this README use them for clarity.

### Installation

You can install the latest version on CRAN like this
```r
    install.packages("assertr")
```

or you can install the bleeding-edge development version like this:
```r
    install.packages("devtools")
    devtools::install_github("ropensci/assertr")
```
### What does it look like?
This package offers five assertion functions, `assert`, `verify`,
`insist`, `assert_rows`, and `insist_rows`, that are designed to be used
shortly after data-loading in an analysis pipeline...

Let’s say, for example, that the R’s built-in car dataset, `mtcars`, was not 
built-in but rather procured from an external source that was known for making
errors in data entry or coding. Pretend we wanted to find the average
miles per gallon for each number of engine cylinders. We might want to first,
confirm
- that it has the columns "mpg", "vs", and "am"
- that the dataset contains more than 10 observations
- that the column for 'miles per gallon' (mpg) is a positive number
- that the column for ‘miles per gallon’ (mpg) does not contain a datum
that is outside 4 standard deviations from its mean, and
- that the am and vs columns (automatic/manual and v/straight engine,
respectively) contain 0s and 1s only
- each row contains at most 2 NAs
- each row is unique *jointly* between the "mpg", "am", and "wt" columns
- each row's mahalanobis distance is within 10 median absolute deviations of
all the distances (for outlier detection)


This could be written (in order) using `assertr` like this:

```r
    library(dplyr)
    library(assertr)

    mtcars %>%
      verify(has_all_names("mpg", "vs", "am", "wt")) %>%
      verify(nrow(.) > 10) %>%
      verify(mpg > 0) %>%
      insist(within_n_sds(4), mpg) %>%
      assert(in_set(0,1), am, vs) %>%
      assert_rows(num_row_NAs, within_bounds(0,2), everything()) %>%
      assert_rows(col_concat, is_uniq, mpg, am, wt) %>%
      insist_rows(maha_dist, within_n_mads(10), everything()) %>%
      group_by(cyl) %>%
      summarise(avg.mpg=mean(mpg))
```

If any of these assertions were violated, an error would have been raised
and the pipeline would have been terminated early.

Let's see what the error message look like when you chain
a bunch of failing assertions together.

```r
    > mtcars %>%
    +   chain_start %>%
    +   assert(in_set(1, 2, 3, 4), carb) %>%
    +   assert_rows(rowMeans, within_bounds(0,5), gear:carb) %>%
    +   verify(nrow(.)==10) %>%
    +   verify(mpg < 32) %>%
    +   chain_end
    There are 7 errors across 4 verbs:
    -
             verb redux_fn           predicate     column index value
    1      assert     <NA>  in_set(1, 2, 3, 4)       carb    30   6.0
    2      assert     <NA>  in_set(1, 2, 3, 4)       carb    31   8.0
    3 assert_rows rowMeans within_bounds(0, 5) ~gear:carb    30   5.5
    4 assert_rows rowMeans within_bounds(0, 5) ~gear:carb    31   6.5
    5      verify     <NA>       nrow(.) == 10       <NA>     1    NA
    6      verify     <NA>            mpg < 32       <NA>    18    NA
    7      verify     <NA>            mpg < 32       <NA>    20    NA

    Error: assertr stopped execution
```

### What does `assertr` give me?

- `verify` - takes a data frame (its first argument is provided by
the `%>%` operator above), and a logical (boolean) expression. Then, `verify`
evaluates that expression using the scope of the provided data frame. If any
of the logical values of the expression's result are `FALSE`, `verify` will
raise an error that terminates any further processing of the pipeline.

- `assert` - takes a data frame, a predicate function, and an arbitrary
number of columns to apply the predicate function to. The predicate function
(a function that returns a logical/boolean value) is then applied to every
element of the columns selected, and will raise an error if it finds any
violations. Internally, the `assert` function uses `dplyr`'s
`select` function to extract the columns to test the predicate function on.

- `insist` - takes a data frame, a predicate-generating function, and an
arbitrary number of columns. For each column, the the predicate-generating
function is applied, returning a predicate. The predicate is then applied to
every element of the columns selected, and will raise an error if it finds any
violations. The reason for using a predicate-generating function to return a
predicate to use against each value in each of the selected rows is so
that, for example, bounds can be dynamically generated based on what the data
look like; this the only way to, say, create bounds that check if each datum is
within x z-scores, since the standard deviation isn't known a priori.
Internally, the `insist` function uses `dplyr`'s `select` function to extract
the columns to test the predicate function on.

- `assert_rows` - takes a data frame, a row reduction function, a predicate
function, and an arbitrary number of columns to apply the predicate function
to. The row reduction function is applied to the data frame, and returns a value
for each row. The predicate function is then applied to every element of vector
returned from the row reduction function, and will raise an error if it finds
any violations. This functionality is useful, for example, in conjunction with
the `num_row_NAs()` function to ensure that there is below a certain number of
missing values in each row. Internally, the `assert_rows` function uses
`dplyr`'s`select` function to extract the columns to test the predicate
function on.

- `insist_rows` - takes a data frame, a row reduction function, a
predicate-generating
function, and an arbitrary number of columns to apply the predicate function
to. The row reduction function is applied to the data frame, and returns a value
for each row. The predicate-generating function is then applied to the vector
returned from the row reduction function and the resultant predicate is
applied to each element of that vector. It will raise an error if it finds any
violations. This functionality is useful, for example, in conjunction with
the `maha_dist()` function to ensure that there are no flagrant outliers.
Internally, the `assert_rows` function uses `dplyr`'s`select` function to
extract the columns to test the predicate function on.


`assertr` also offers four (so far) predicate functions designed to be used
with the `assert` and `assert_rows` functions:

- `not_na` - that checks if an element is not NA
- `within_bounds` - that returns a predicate function that checks if a numeric
value falls within the bounds supplied, and
- `in_set` - that returns a predicate function that checks if an element is
a member of the set supplied. (also allows inverse for "not in set")
- `is_uniq` - that checks to see if each element appears only once


and predicate generators designed to be used with the `insist` and `insist_rows`
functions:

- `within_n_sds` - used to dynamically create bounds to check vector elements with
based on standard z-scores
- `within_n_mads` - better method for dynamically creating bounds to check vector
elements with based on 'robust' z-scores (using median absolute deviation)

and the following row reduction functions designed to be used with `assert_rows`
and `insist_rows`:

- `num_row_NAs` - counts number of missing values in each row
- `maha_dist` - computes the mahalanobis distance of each row (for outlier
detection). It will coerce categorical variables into numerics if it needs to.
- `col_concat` - concatenates all rows into strings
- `duplicated_across_cols` - checking if a row contains a duplicated value
across columns

and, finally, some other utilities for use with `verify`

- `has_all_names` - check if the data frame or list has all supplied names
- `has_only_names` - check that a data frame or list have _only_ the names
requested
- `has_class` - checks if passed data has a particular class


### More info

For more info, check out the `assertr` vignette
```r
    > vignette("assertr")
```
Or [read it here](https://CRAN.R-project.org/package=assertr/vignettes/assertr.html)

# [![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
## Test environments
* local Debian 10 (buster) install, R 4.0.0
* local macOS High Sierra 10.13 install, R 4.0.3
* win-builder (devel, release, and old_release)

## R CMD check results

When checked locally with --as-cran and --no-manual
there were no ERRORs, WARNINGs, or NOTEs
---
title: "Assertive R Programming with assertr"
author: "Tony Fischetti"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Assertive R Programming with assertr}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


In data analysis workflows that depend on un-sanitized data sets from external
sources, it’s very common that errors in data bring an analysis to a
screeching halt. Oftentimes, these errors occur late in the analysis and
provide no clear indication of which datum caused the error.

On occasion, the error resulting from bad data won’t even appear to be a
data error at all. Still worse, errors in data will pass through analysis
without error, remain undetected, and produce inaccurate results.

The solution to the problem is to provide as much information as you can about
how you expect the data to look up front so that any deviation from this
expectation can be dealt with immediately. This is what the `assertr` package
tries to make dead simple. 

Essentially, `assertr` provides a suite of functions designed to verify
assumptions about data early in an analysis pipeline. This package needn't
be used with the `magrittr`/`dplyr` piping mechanism but the examples in this
vignette will use them to enhance clarity.

### concrete data errors

Let’s say, for example, that the R’s built-in car dataset, `mtcars`, was not 
built-in but rather procured from an external source that was known for making
errors in data entry or coding.

In particular, the mtcars dataset looks like this:
```{r}
head(mtcars)
```

But let's pretend that the data we got accidentally negated the 5th mpg value:

```{r}
our.data <- mtcars
our.data$mpg[5] <- our.data$mpg[5] * -1
our.data[4:6,]
```

Whoops!

If we wanted to find the average miles per gallon for each number of engine
cylinders, we might do so like this:

```{r message=FALSE}
library(dplyr)

our.data %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))

```

This indicates that the average miles per gallon for a 8 cylinder car is a lowly
12.43. However, in the correct dataset it's really just over 15. Data errors
like that are extremely easy to miss because it doesn't cause an error, and the
results look reasonable.

### enter assertr

To combat this, we might want to use assertr's `verify` function to make sure
that `mpg` is a positive number:

```{r error=TRUE, purl = FALSE, include=TRUE}
library(assertr)

our.data %>%
  verify(mpg >= 0) %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

If we had done this, we would have caught this data error.

The `verify` function takes a data frame (its first argument is provided by
the `%>%` operator), and a logical (boolean) expression. Then, `verify`
evaluates that expression using the scope of the provided data frame. If any
of the logical values of the expression's result are `FALSE`, `verify` will
raise an error that terminates any further processing of the pipeline.

We could have also written this assertion using `assertr`'s `assert` function...

```{r error=TRUE, purl = FALSE}
our.data %>%
  assert(within_bounds(0,Inf), mpg) %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

The `assert` function takes a data frame, a predicate function, and an arbitrary
number of columns to apply the predicate function to. The predicate function
(a function that returns a logical/boolean value) is then applied to every
element of the columns selected, and will raise an error when if it finds 
violations.

Internally, the `assert` function uses `dplyr`'s `select` function to extract
the columns to test the predicate function on. This allows for complex
assertions. Let's say we wanted to make sure that all values in the dataset
are *greater* than zero (except `mpg`):

```{r error=TRUE, purl = FALSE}
library(assertr)

our.data %>%
  assert(within_bounds(0,Inf, include.lower=FALSE), -mpg) %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

### verify vs. assert

The first noticable difference between `verify` and `assert` is that `verify`
takes an expression, and `assert` takes a predicate and columns to apply it to.
This might make the `verify` function look more elegant--but there's an
important drawback. `verify` has to evaluate the entire expression first, and
*then* check if there were any violations. Because of this, `verify` can't
tell you the offending datum.

One important drawback to `assert`, and a consequence of its application of
the predicate to *columns*, is that `assert` can't confirm assertions about
the data structure *itself*. For example, let's say we were reading a dataset
from disk that we know has more than 100 observations; we could write a check
of that assumption like this:

```{r eval=FALSE, purl = FALSE}
dat <- read.csv("a-data-file.csv")
dat %>%
  verify(nrow(.) > 100) %>%
  ....
```

Other checking functions that are only available to `verify` are
`has_all_names`, `has_only_names`, and `has_class`.

This is a powerful advantage over `assert`... but `assert` has one more
advantage of its own that we've heretofore ignored.

### assertr's predicates

`assertr`'s predicates, both built-in and custom, make `assert` very powerful.
The predicates that are built in to `assertr` are

- `not_na` - that checks if an element is not NA
- `within_bounds` - that returns a predicate function that checks if a numeric
value falls within the bounds supplied, and
- `in_set` - that returns a predicate function that checks if an element is
a member of the set supplied.
- `is_uniq` - that checks to see if each element appears only once

We've already seen `within_bounds` in action... let's use the `in_set` function
to make sure that there are only 0s and 1s (automatic and manual, respectively)
values in the `am` column...


```{r, eval=FALSE, purl = FALSE}
our.data %>%
  assert(in_set(0,1), am) %>%
  ...
```

If we were reading a dataset that contained a column representing boroughs of
New York City (named `BORO`), we can verify that there are no mis-spelled
or otherwise unexpected boroughs like so...

```{r, eval=FALSE, purl = FALSE}
boroughs <- c("Bronx", "Manhattan", "Queens", "Brooklyn", "Staten Island")

read.csv("a-dataset.csv") %>%
  assert(in_set(boroughs), BORO) %>%
  ...
```

Rad!

### custom predicates

A convenient feature of `assertr` is that it makes the construction of custom
predicate functions easy.

In order to make a custom predicate, you only have to specify cases where the
predicate should return FALSE. Let's say that a dataset has an ID column
(named `ID`) that we want to check is not an empty string. We can create a
predicate like this:

```{r}
not.empty.p <- function(x) if(x=="") return(FALSE)
```

and apply it like this:

```{r, eval=FALSE, purl = FALSE}
read.csv("another-dataset.csv") %>%
  assert(not.empty.p, ID) %>%
  ...
```

Let's say that the ID column is always a 7-digit number. We can confirm that
all the IDs are 7-digits by defining the following predicate:

```{r}
seven.digit.p <- function(x) nchar(x)==7
```

A powerful consequence of this easy creation of predicates is that the
`assert` function lends itself to use with lambda predicates (unnamed
predicates that are only used once). The check above might be better written as

```{r, eval=FALSE, purl = FALSE}
read.csv("another-dataset.csv") %>%
  assert(function(x) nchar(x)==7, ID) %>%
  ...
```

Neat-o!


### enter `insist` and predicate 'generators'

Very often, there is a need to dynamically determine the predicate function
to be used based on the vector being checked.

For example, to check to see if every element of a vector is within _n_
standard deviations of the mean, you need to create a `within_bounds`
predicate _after_ dynamically determining the bounds by reading and computing
on the vector itself.

To this end, the `assert` function is no good; it just applies a raw predicate
to a vector. We need a function like `assert` that will apply predicate
_generators_ to vectors, return predicates, and _then_ perform `assert`-like
functionality by checking each element of the vectors with its respective custom
predicate. This is precisely what `insist` does.

This is all much simpler than it may sound. Hopefully, the examples will clear
up any confusion.

The primary use case for `insist` is in conjunction with the `within_n_sds` or
`within_n_mads` predicate generator.

Suppose we wanted to check that every `mpg` value in the `mtcars` data set was
within 3 standard deviations of the mean before finding the average miles
per gallon for each number of engine cylinders. We could write something
like this:

```{r purl = FALSE}

mtcars %>%
  insist(within_n_sds(3), mpg) %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

Notice what happens when we drop that z-score to 2 standard deviations
from the mean

```{r error=TRUE,  purl = FALSE}
mtcars %>%
  insist(within_n_sds(2), mpg) %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

Execution of the pipeline was halted. But now we know exactly which data point
violated the predicate that `within_n_sds(2)(mtcars$mpg)`
returned.

Now that's an efficient car!

After the predicate generator, `insist` takes an arbitrary number of columns
just like `assert` using the syntax of `dplyr`'s `select` function. If you
wanted to check that everything in mtcars is within 10 standard deviations
of the mean (of each column vector), you can do so like this:

```{r purl = FALSE}
mtcars %>%
  insist(within_n_sds(10), mpg:carb) %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

Aces!

I chose to use `within_n_sds` in this example because people are familiar
z-scores. However, for most practical purposes, the related predicate generator
`within_n_mads` is more useful.

The problem with `within_n_sds` is the mean and standard deviation are so
heavily influenced by outliers, their very presence will compromise attempts
to identify them using these statistics. In contrast with `within_n_sds`, 
`within_n_mads` uses the robust statistics, median and median absolute
deviation, to identify potentially erroneous data points.

For example, the vector `<7.4, 7.1, 7.2, 72.1>` almost certainly has an erroneous
data point, but `within_n_sds(2)` will fail to detect it.

```{r purl = FALSE}
example.vector <- c(7.4, 7.1, 7.2, 72.1)
within_n_sds(2)(example.vector)(example.vector)
```

whereas `within_n_mads` will detect it at even lower levels of power....

```{r purl = FALSE}
example.vector <- c(7.4, 7.1, 7.2, 72.1)
within_n_mads(2)(example.vector)(example.vector)
within_n_mads(1)(example.vector)(example.vector)
```

Tubular!


### row-wise assertions and row reduction functions

As cool as it's been so far, this still isn't enough to constitute a complete
grammar of data integrity checking. To see why, check out the following
small example data set:

```{r perl=FALSE}
example.data <- data.frame(x=c(8, 9, 6, 5, 9, 5, 6, 7,
                             8, 9, 6, 5, 5, 6, 7),
                         y=c(82, 91, 61, 49, 40, 49, 57,
                             74, 78, 90, 61, 49, 51, 62, 68))
(example.data)
```

Can you spot the brazen outlier? You're certainly not going to find it by
checking the distribution of each *column*! All elements from both columns are
within 2 standard deviations of their respective means.

Unless you have a *really* good eye, the only way you're going to catch this
mistake is by plotting the data set.

```{r purl=FALSE}
plot(example.data$x, example.data$y, xlab="", ylab="")
```

Ok, so all the `y`s are roughly 10 times the `x`s except the outlying data
point.

The problem having to plot data sets to catch anomalies is that it is *really*
hard to visualize 4-dimensions at once, and it is near impossible with
high-dimensional data.

There's no way of catching this anomaly by looking at each individual
column separately; the only way to catch it is to view each row as a complete
observation and compare it to the rest.

To this end, `assertr` provides three functions that take a data frame, and
reduce each row into a single value. We'll call them *row reduction functions*.

The first one we'll look at is called `maha_dist`. It computes the average
mahalanobis distance (kind of like multivariate z-scoring for outlier
detection) of each row from the whole data set. The big idea is that in the
resultant vector, big/distant values are potential anomalous entries. Let's
look at the distribution of mahalanobis distances for this data set...

```{r purl=FALSE}
maha_dist(example.data)

maha_dist(example.data) %>% hist(main="", xlab="")
```

There's no question here as to whether there's an anomalous entry! But how do
you check for this sort of thing using `assertr` constructs?

Well, `maha_dist` will typically be used with the `insist_rows` function.
`insist_rows` takes a data frame, a row reduction function, a
predicate-generating function, and an arbitrary number of columns to apply
the predicate function to. The row reduction function (`maha_dist` in this case)
is applied to the data frame, and returns a value for each row. The
predicate-generating function is then applied to the vector returned from
the row reduction function and the resultant predicate is applied to each
element of that vector. It will raise an error if it finds any violations.

As always, this undoubtedly sounds far more confusing than it really is. Here's
an example of it in use

```{r purl=FALSE, error=TRUE}
example.data %>%
  insist_rows(maha_dist, within_n_mads(3), dplyr::everything())

```
 
Check that out! To be clear, this function is running the supplied data frame
through the `maha_dist` function which returns a value for each row
corresponding to its mahalanobis distance. (The whole data frame is used because
we used the `everything()` selection function from the `dplyr` package.)
Then, `within_n_mads(3)` computes on that vector and returns a bounds
checking predicate. The bounds checking predicate checks to see that all
mahalanobis distances are within 3 median absolute deviations
of each other. They are not, and the pipeline errors out. Note that the
`data.frame` of errors that is returned by error report contains the
verb used (`insist_rows`), the row reduction function, the predicate,
the column (or columns), the index of the failure and the offendind datum.

This is probably the most powerful construct in `assertr`--it can find a whole
lot of nasty errors that would be very difficult to check for by hand.

Part of what makes it so powerful is how flexible `maha_dist` is. We only used
it, so far, on a data frame of numerics, but it can handle all sorts of data
frames. To really see it shine, let's use it on the iris data set, that contains
a categorical variable in its right-most column...

```{r purl=FALSE}
head(iris)

iris %>% maha_dist %>% hist(main="", xlab="")
```

Looks ok, but what happens when we accidentally enter a row as a different
species...

```{r purl=FALSE}
mistake <- iris
(mistake[149,5])
mistake[149,5] <- "setosa"

mistake %>% maha_dist %>% hist(main="", xlab="")

mistake %>% maha_dist %>% which.max
```

Look at that! This mistake can easily be picked up by any reasonable bounds
checker...

```{r purl=FALSE, error=TRUE}
mistake %>% insist_rows(maha_dist, within_n_mads(7), dplyr::everything())
```

`insist` and `insist_rows` are both similar in that they both take predicate
generators and not actual predicates. What makes `insist_rows` different is
its usage of a row-reduce data frame.

`assert` has a row-oriented counterpart, too; it's called `assert_rows`. 
`insist` is to `assert` as `insist_rows` is to `assert_rows`.

`assert_rows` works the same as `insist_rows`, except that instead of using
a predicate generator on the row-reduced data frame, it uses a regular-old
predicate.

For an example of a `assert_rows` use case, let's say that we got a data set
(`another-dataset.csv`) from the web and we don't want to continue processing
the data set if any row contains more than two missing values (NAs). You
can use the row reduction function `num_row_NAs` to reduce all the rows into
the number of NAs they contain. Then, a simple bounds checker will suffice for
ensuring that no element is higher than 2...


```{r, eval=FALSE, purl = FALSE}
read.csv("another-dataset.csv") %>%
  assert_rows(num_row_NAs, within_bounds(0,2), dplyr::everything()) %>%
  ...
```

`assert_rows` can be used for anomaly detection as well. A future version of
`assertr` may contain a cosine distance row reduction function. Since all
cosine distances are constained from -1 to 1, it is easy to use a non-dynamic
predicate to disallow certain values.


### success, error and defect functions

The behavior of functions like `assert`, `assert_rows`,
`insist`, `insist_rows`, `verify` when the assertion
passes, fails or is skipped due to data defect is configurable 
via the `success_fun`, `error_fun` and `defect_fun` parameters, 
respectively.

The `success_fun` parameter takes a function that takes
the data passed to the assertion function as a parameter. You can
write your own success handler function, but there are a few
provided by this package:

  * `success_continue` - just returns the data that was passed into the
    assertion function (this is default).

  * `success_logical` - returns TRUE
  
  * `success_append` - returns the data that was passed into the assertion 
    function but also attaches basic information about verification result
    to a special attribute of `data`.

  * `success_report` - when success results are stored 
    (`chain_start(store_results=TRUE)`), and each verification ended up with 
    success, it prints a summary of all successful validations (when being in chan) 
    or simple verification result for single check and returns data.
    
  * `success_df_return` - when success results are stored 
    (`chain_start(store_results=TRUE)`), and each verification ended up with success, 
    it prints data.frame with verification results (can be used for `chain_end` or 
    single verification).

The `error_fun` parameter takes a function that takes
the data passed to the assertion function as a parameter. You can
write your own error handler function, but there are a few
provided by this package:

  * `error_stop` - Prints a summary of the errors and
                            halts execution (default)

  * `error_report` - Prints *all* the information available
                              about the errors and halts execution.

  * `error_append` - Attaches the errors to a special
   attribute of `data` and returns the data. This is chiefly
   to allow assertr errors to be accumulated in a pipeline so that
   all assertions can have a chance to be checked and so that all
   the errors can be displayed at the end of the chain.

  * `error_logical` - returns FALSE

  * `just_warn` - Prints a summary of the errors but does
   not halt execution, it just issues a warning.

  * `warn_report` - Prints all the information available
  about the errors but does not halt execution, it just issues a warning.
  
  * `defect_report` - For single rule and defective data it displays
  short info about skipping the current assertion. For `chain_end` sums
  up all skipped rules for defective data. 
  
  * `defect_df_return` - For single rule and defective data it returns
  info data.frame about skipping current assertion. For `chain_end`
  returns all skipped rules info data.frame for defective data.

The `defect_fun` parameter takes a function that takes
the data passed to the assertion function as a parameter. 
Defect handler is called when any of previous assertions 
that was marked as obligatory failed (see below section).
You can write your own defect handler function, but there are a few
provided by this package:

  * `defect_append` - Attaches the assertion call info on defective 
  data to a special attribute of `data` and returns the data.

  * `defect_report` - For single rule and defective data it displays
  short info about skipping current assertion. For `chain_end` sums
  up all skipped rules for defective data. 
  
  * `defect_df_return` - For single rule and defective data it returns
  info data.frame about skipping current assertion. For `chain_end`
  returns all skipped rules info data.frame for defective data.

### Obligatory assertions

You may find a situation in which some rules are not independent.

For example:
```{r, eval=FALSE, purl = FALSE}
mtcars_without_am <- mtcars %>% 
  dplyr::select(-am)
mtcars_without_am %>% 
  verify(has_all_names("am", "vs"), error_fun = error_append) %>% 
  assert(in_set(0, 1), am, vs, error_fun = error_report)
  
```

`assert` requires existence of `am` and `vs` columns, which are checked
previously by `verify` assertion.
In the above example, we want to store info about all errors after
assertions are finished but it won't happen. A missing `am` column in 
the data returns error not related to assertion check. As a result
code execution ends up with not handled error.

To allow such situations obligatory rules were introduced.
You can create an obligatory rule by adding `obligatory = TRUE` inside
assertion function.
When a rule was obligatory and failed, the data is marked as defective
and each following rule will be handled by `defect_fun` function.
By default `defect_fun=defect_append` which registers information about
running assertion on defective data and skips the rule execution.

Below we display information about a skipped rule by using `defect_report`:
```{r purl = FALSE}
mtcars_without_am <- mtcars %>% 
  dplyr::select(-am)
mtcars_without_am %>% 
  verify(has_all_names("am", "vs"), obligatory=TRUE, error_fun=error_append) %>% 
  assert(in_set(0, 1), am, vs, defect_fun=defect_report)
  
```

### combining chains of assertions

Let's say that as part of an automated pipeline that grabs mtcars from an
untrusted source and finds the average miles per gallon for each number of
engine cylinders, we want to perform the following checks...

- that it has the columns "mpg", "vs", and "am"
- that the dataset contains more than 10 observations
- that the column for 'miles per gallon' (mpg) is a positive number
- that the column for 'miles per gallon' (mpg) does not contain a
datum that is outside 4 standard deviations from its mean, and
- that the am and vs columns (automatic/manual and v/straight engine, 
respectively) contain 0s and 1s only

This could be written thusly:

```{r purl = FALSE}
mtcars %>%
  verify(has_all_names("mpg", "vs", "am")) %>%
  verify(nrow(mtcars) > 10) %>%
  verify(mpg > 0) %>%
  insist(within_n_sds(4), mpg) %>%
  assert(in_set(0,1), am, vs) %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

In an assertr chain with default options, `assert`, `assert_rows`,
`insist`, `insist_rows`, and `verify` will stop at the
first assertion that yields an error and not go on to process the
assertions further down in the chain. For some needs, this is sensible
behavior. There are times, however, when we might like to get a report
of all assertion violations. For example, one might want to write an R
program to download some dataset from the internet and get a detailed
report of all deviations from expectation.

The best thing to do for this use case, is to use the `chain_start`,
and `chain_end` functions at the beginning and end of a chain of
assertr assertions. When `chain_start` gets called with data, the
data gets a special tag that tells the assertr assertions that follow
to override their `success_fun` and `error_fun` values and
replace them with `success_continue` (which passes the data along
if the test passes) and `error_append` (which we've just discussed).
After all relevant verifications, `chain_end` will receive the
data (possibly with accumulated error messages attached) and, by default,
print a report of all the errors that have been found since the start of
the chain.
 
Let's see it in action!

```{r purl = FALSE}
mtcars %>%
  chain_start %>%
  verify(nrow(mtcars) > 10) %>%
  verify(mpg > 0) %>%
  insist(within_n_sds(4), mpg) %>%
  assert(in_set(0,1), am, vs) %>%
  chain_end %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

Now *all* assertions will be checked and reported.

Tip: we can make this whole thing look a lot better by abstracting out
all the assertions:

```{r eval = FALSE, purl = FALSE}
check_me <- . %>%
  chain_start %>%
  verify(nrow(mtcars) > 10) %>%
  verify(mpg > 0) %>%
  insist(within_n_sds(4), mpg) %>%
  assert(in_set(0,1), am, vs) %>%
  chain_end 

mtcars %>%
  check_me %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

Awesome! Now we can add an arbitrary number of assertions, as the need arises,
without touching the real logic.

Note: By default, all assertions in the chain use `success_continue` and 
`error_append` functions.
It allows you to continue workflow in both cases and aggregate error logs.
In some cases (i.e. when you don't want to include the error log and just
have some error printed),
you can require using assertion specific success/error callback in chain.

Just use `skip_chain_opts = TRUE` and specify callback inside assertion:
```
print_error <- function(errors, data=NULL) {
  print(errors)
  return(data)
}
mtcars %>%
  chain_start %>%
  verify(nrow(mtcars) > 32, error_fun=print_error, skip_chain_opts=TRUE) %>%
  verify(mpg > 0) %>%
  insist(within_n_sds(4), mpg) %>%
  assert(in_set(0,1), am, vs) %>%
  chain_end 
```

You may also store validation successful results by using `store_success=TRUE`:
```
mtcars %>%
  chain_start(store_success=TRUE) %>%
  verify(nrow(mtcars) == 32) %>%
  verify(mpg > 0) %>%
  insist(within_n_sds(4), mpg) %>%
  assert(in_set(0,1), am, vs) %>%
  chain_end(success_fun=success_df_return)
```


### advanced: send email reports using custom error functions

One particularly cool application of `assertr` is to use it as a data integrity
checker for frequently updated data sources. A script can download new data as
it becomes available, and then run `assertr` checks on it. This makes `assertr`
into a sort of "continuous integration" tool (but for data,
not code.)

In an unsupervised "continuous integration" environment, you need a way to
discover that the assertions failed. In CI-as-a-service in the software world, 
failed automated checks often send an email of reporting the maintainer of a
botched build; why not bring that functionality to `assertr`?!

As we learned in the last sections, all assertion verbs in `assertr`
support a custom error function. `chain_end` similarly supports custom
error functions. By default, this is `error_stop` (or `error_report` in the
case of `chain_end`) which prints a summary of the errors and halts execution.

You can specify your own, though, to hijack this behavior and redirect
flow-of-control wherever you want. 

Your custom error function must take, as its first argument,
a list of `assertr_error` S3 objects. The second argument must be the
`data.frame` that the verb is computing on. Every error function must
take this because there may be some other errors that are attached to the
`data.frame`'s `assertr_errors` attribute leftover from other previous
assertions.

Below we are going to build a function that takes a list of `assertr_errors`,
gets a string representation of the errors and emails it to someone before
halting execution. We will use the `mailR` package to send the mail.

```{r eval=FALSE, purl = FALSE}

library(mailR)

email_me <- function(list_of_errors, data=NULL, ...){
  # we are checking to see if there are any errors that
  # are still attached to the data.frame
  if(!is.null(data) && !is.null(attr(data, "assertr_errors")))
    errors <- append(attr(data, "assertr_errors"), errors)

  num.of.errors <- length(list_of_errors)

  preface <- sprintf("There %s %d error%s:\n",
                     ifelse(num.of.errors==1,"is", "are"),
                     num.of.errors,
                     ifelse(num.of.errors==1,"", "s"))

  # all `assertr_error` S3 objects have `print` and `summary` methods
  # here, we will call `print` on all of the errors since `print`
  # will give us the complete/unabridged error report
  error_string <- capture.output(tmp <- lapply(list_of_errors,
                                               function(x){
                                                 cat("\n- ");
                                                 print(x);
                                                 return();}))
  error_string <- c(preface, error_string)
  error_string <- error_string %>% paste0(collapse="\n")

  send.mail(from="assertr@gmail.com", to="YOU@gmail.com",
            subject="error from assertr", body=error_string,
            smtp = list(host.name="aspmx.l.google.com", port=25),
            authenticate = FALSE, send=TRUE)
  stop("assertr stopped execution", call.=FALSE)
}

questionable_mtcars %>%
  chain_start %>%
  verify(nrow(.) > 10) %>%
  insist(within_n_sds(4), mpg) %>%
  # ...
  chain_end(error_fun=email.me)
```

(this particular `send.mail` formulation will only work for gmail
recipients; see the `mailR` documentation for more information)

Now you'll get notified of <s>any</s> all failed assertions via email. Groovy!


### advanced: creating your own predicate generators for `insist`

`assertr` is build with robustness, correctness, and extensibility in mind.
Just like `assertr` makes it easy to create your own custom predicates, so
too does this package make it easy to create your own custom predicate
generators.

Okay... so its, perhaps, not _easy_ because predicate generators by nature
are functions that return functions. But it's possible!

Let's say you wanted to create a predicate generator that checks if all
elements of a vector are within 3 times the vector's interquartile range from
the median. We need to create a function that looks like this

```{r purl = FALSE}
within_3_iqrs <- function(a_vector){
  the_median <- median(a_vector)
  the_iqr <- IQR(a_vector)
  within_bounds((the_median-the_iqr*3), (the_median+the_iqr*3))
}
```

Now, we can use it on `mpg` from `mtcars` like so:

```{r purl = FALSE}
mtcars %>%
  insist(within_3_iqrs, mpg) %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

There are two problems with this, though...

1. We may want to abstract this so that we can supply an arbitrary number
of IQRs to create the bounds with
2. We lose the ability to choose what optional arguments (if any) that we
give to the returned `within_bounds` predicate.

Now we have to write a function that returns a function that returns
a function...

```{r purl = FALSE}
within_n_iqrs <- function(n, ...){
  function(a_vector){
    the_median <- median(a_vector)
    the_iqr <- IQR(a_vector)
    within_bounds((the_median-the_iqr*n), (the_median+the_iqr*n), ...)
  }
}
```

Much better! Now, if we want to check that every `mpg` from `mtcars` is
within 5 IQRs of the median and *not allow NA values* we can do so like this:

```{r purl = FALSE}
mtcars %>%
  insist(within_n_iqrs(5), mpg) %>%
  group_by(cyl) %>%
  summarise(avg.mpg=mean(mpg))
```

Super!


### advanced: programming with assertion functions

These assertion functions use the
[tidyeval](https://rpubs.com/hadley/dplyr-programming) framework. In the
past, programming in a tidyverse-like setting was accomplished through
standard evaluation versions of verbs, which used functions postfixed
with an underscore: `insist_` instead of `insist`, for example. However,
when tidyeval was made popular with `dplyr` 0.7.0, this usage became deprecated,
and therefore underscore-postfixed functions are no longer part of `assertr`.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R, R/deprecated_se.R
\name{assert_rows}
\alias{assert_rows}
\alias{assert_rows_}
\title{Raises error if predicate is FALSE for any row after applying
row reduction function}
\usage{
assert_rows(
  data,
  row_reduction_fn,
  predicate,
  ...,
  success_fun = success_continue,
  error_fun = error_stop,
  skip_chain_opts = FALSE,
  obligatory = FALSE,
  defect_fun = defect_append,
  description = NA
)

assert_rows_(
  data,
  row_reduction_fn,
  predicate,
  ...,
  .dots,
  success_fun = success_continue,
  error_fun = error_stop
)
}
\arguments{
\item{data}{A data frame}

\item{row_reduction_fn}{A function that returns a value for each row of
the provided data frame}

\item{predicate}{A function that returns FALSE when violated}

\item{...}{Comma separated list of unquoted expressions.
Uses dplyr's \code{select} to select
columns from data.}

\item{success_fun}{Function to call if assertion passes. Defaults to
returning \code{data}.}

\item{error_fun}{Function to call if assertion fails. Defaults to printing
a summary of all errors.}

\item{skip_chain_opts}{If TRUE, \code{success_fun} and \code{error_fun}
are used even if assertion is called within a chain.}

\item{obligatory}{If TRUE and assertion failed the data is marked as defective.
For defective data, all the following rules are handled by
\code{defect_fun} function.}

\item{defect_fun}{Function to call when data is defective. Defaults to skipping
assertion and storing info about it in special attribute.}

\item{description}{Custom description of the rule. Is stored in result
reports and data.}

\item{.dots}{Use assert_rows_() to select columns using standard evaluation.}
}
\value{
By default, the \code{data} is returned if predicate assertion
        is TRUE and and error is thrown if not. If a non-default
        \code{success_fun} or \code{error_fun} is used, the return
        values of these function will be returned.
}
\description{
Meant for use in a data analysis pipeline, this function applies a
function to a data frame that reduces each row to a single value. Then,
a predicate function is applied to each of the row reduction values. If
any of these predicate applications yield FALSE, this function will raise
an error, effectively terminating the pipeline early. If there are no
FALSEs, this function will just return the data that it was supplied for
further use in later parts of the pipeline.
}
\details{
For examples of possible choices for the \code{success_fun} and
\code{error_fun} parameters, run \code{help("success_and_error_functions")}
}
\note{
See \code{vignette("assertr")} for how to use this in context
}
\examples{

# returns mtcars
assert_rows(mtcars, num_row_NAs, within_bounds(0,2), mpg:carb)

library(magrittr)                    # for piping operator

mtcars \%>\%
  assert_rows(rowSums, within_bounds(0,2), vs:am)
  # anything here will run

\dontrun{
mtcars \%>\%
  assert_rows(rowSums, within_bounds(0,1), vs:am)
  # the assertion is untrue so
  # nothing here will run}

}
\seealso{
\code{\link{insist_rows}} \code{\link{assert}}
         \code{\link{verify}} \code{\link{insist}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{has_class}
\alias{has_class}
\title{Returns TRUE if data.frame columns have a specified class}
\usage{
has_class(..., class)
}
\arguments{
\item{...}{An arbitrary amount of quoted column names to check for}

\item{class}{Expected class for chosen columns.}
}
\value{
TRUE if all classes are correct, FALSE if not
}
\description{
This is meant to be used with `assertr`'s `verify` function to check
for the existence of a specific column class in a `data.frame` that is
piped to `verify`.
}
\examples{

verify(mtcars, has_class("mpg", "wt", class = "numeric"))

library(magrittr)   # for pipe operator

\dontrun{
mtcars \%>\%
  verify(has_class("mpg", class = "character"))  # fails
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{summary.assertr_verify_error}
\alias{summary.assertr_verify_error}
\title{Summarizing assertr's verify errors}
\usage{
\method{summary}{assertr_verify_error}(object, ...)
}
\arguments{
\item{object}{An assertr_verify_error object}

\item{...}{Additional arguments affecting the summary produced}
}
\description{
`summary` method for class "assertr_verify_error"
}
\seealso{
\code{\link{print.assertr_verify_error}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/row-redux.R
\name{num_row_NAs}
\alias{num_row_NAs}
\title{Counts number of NAs in each row}
\usage{
num_row_NAs(data, allow.NaN = FALSE)
}
\arguments{
\item{data}{A data frame}

\item{allow.NaN}{Treat NaN like NA (by counting it). FALSE by default}
}
\value{
A vector of number of missing values in each row
}
\description{
This function will return a vector, with the same length as the number
of rows of the provided data frame, corresponding to the number of missing
values in each row
}
\examples{

num_row_NAs(mtcars)


library(magrittr)            # for piping operator
library(dplyr)               # for "everything()" function

# using every column from mtcars, make sure there are at most
# 2 NAs in each row. If there are any more than two, error out
mtcars \%>\%
  assert_rows(num_row_NAs, within_bounds(0,2), everything())
  ## anything here will run

}
\seealso{
\code{\link{is.na}} \code{\link{is.nan}} \code{\link{not_na}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{chaining_functions}
\alias{chaining_functions}
\alias{chain_start}
\alias{chain_end}
\title{Chaining functions}
\usage{
chain_start(data, store_success = FALSE)

chain_end(data, success_fun = success_continue, error_fun = error_report)
}
\arguments{
\item{data}{A data frame}

\item{store_success}{If TRUE each successful assertion is stored in chain.}

\item{success_fun}{Function to call if assertion passes. Defaults to
returning \code{data}.}

\item{error_fun}{Function to call if assertion fails. Defaults to printing
a summary of all errors.}
}
\description{
These functions are for starting and ending a sequence of assertr
assertions and overriding the default behavior of assertr halting
execution on the first error.
}
\details{
For more information, read the relevant section in this package's
vignette using, \code{vignette("assertr")}

For examples of possible choices for the \code{success_fun} and
\code{error_fun} parameters, run \code{help("success_and_error_functions")}
}
\examples{
library(magrittr)

mtcars \%>\%
  chain_start() \%>\%
  verify(nrow(mtcars) > 10) \%>\%
  verify(mpg > 0) \%>\%
  insist(within_n_sds(4), mpg) \%>\%
  assert(in_set(0,1), am, vs) \%>\%
  chain_end()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R, R/deprecated_se.R
\name{insist}
\alias{insist}
\alias{insist_}
\title{Raises error if dynamically created predicate is FALSE in any columns
selected}
\usage{
insist(
  data,
  predicate_generator,
  ...,
  success_fun = success_continue,
  error_fun = error_stop,
  skip_chain_opts = FALSE,
  obligatory = FALSE,
  defect_fun = defect_append,
  description = NA
)

insist_(
  data,
  predicate_generator,
  ...,
  .dots,
  success_fun = success_continue,
  error_fun = error_stop
)
}
\arguments{
\item{data}{A data frame}

\item{predicate_generator}{A function that is applied
to each of the column vectors selected. This will produce,
for every column, a true predicate function to be applied to
every element in the column vectors selected}

\item{...}{Comma separated list of unquoted expressions.
Uses dplyr's \code{select} to select
columns from data.}

\item{success_fun}{Function to call if assertion passes. Defaults to
returning \code{data}.}

\item{error_fun}{Function to call if assertion fails. Defaults to printing
a summary of all errors.}

\item{skip_chain_opts}{If TRUE, \code{success_fun} and \code{error_fun}
are used even if assertion is called within a chain.}

\item{obligatory}{If TRUE and assertion failed the data is marked as defective.
For defective data, all the following rules are handled by
\code{defect_fun} function.}

\item{defect_fun}{Function to call when data is defective. Defaults to skipping
assertion and storing info about it in special attribute.}

\item{description}{Custom description of the rule. Is stored in result
reports and data.}

\item{.dots}{Use insist_() to select columns using standard evaluation.}
}
\value{
By default, the \code{data} is returned if dynamically created
        predicate assertion is TRUE and and error is thrown if not. If a
        non-default \code{success_fun} or \code{error_fun} is used, the
        return values of these function will be returned.
}
\description{
Meant for use in a data analysis pipeline, this function applies a predicate
generating function to each of the columns indicated. It will then use these
predicates to check every element of those columns. If any of these
predicate applications yield FALSE, this function will raise an error,
effectively terminating the pipeline early. If there are no FALSES, this
function will just return the data that it was supplied for further use in
later parts of the pipeline.
}
\details{
For examples of possible choices for the \code{success_fun} and
\code{error_fun} parameters, run \code{help("success_and_error_functions")}
}
\note{
See \code{vignette("assertr")} for how to use this in context
}
\examples{

insist(iris, within_n_sds(3), Sepal.Length)   # returns iris

library(magrittr)

iris \%>\%
  insist(within_n_sds(4), Sepal.Length:Petal.Width)
  # anything here will run

\dontrun{
iris \%>\%
  insist(within_n_sds(3), Sepal.Length:Petal.Width)
  # datum at index 16 of 'Sepal.Width' vector is (4.4)
  # is outside 3 standard deviations from the mean of Sepal.Width.
  # The check fails, raises a fatal error, and the pipeline
  # is terminated so nothing after this statement will run}

}
\seealso{
\code{\link{assert}} \code{\link{verify}} \code{\link{insist_rows}}
         \code{\link{assert_rows}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertr.R
\docType{package}
\name{assertr}
\alias{assertr}
\title{assertr: Assertive programming for R analysis pipeline.}
\description{
The assertr package supplies a suite of functions designed to verify
assumptions about data early in an analysis pipeline.
See the assertr vignette or the documentation for more information \cr
> \code{vignette("assertr")}
}
\details{
You may also want to read the documentation for the functions that
\code{assertr} provides:
\itemize{
  \item \code{\link{assert}}
  \item \code{\link{verify}}
  \item \code{\link{insist}}
  \item \code{\link{assert_rows}}
  \item \code{\link{insist_rows}}
  \item \code{\link{not_na}}
  \item \code{\link{in_set}}
  \item \code{\link{has_all_names}}
  \item \code{\link{is_uniq}}
  \item \code{\link{num_row_NAs}}
  \item \code{\link{maha_dist}}
  \item \code{\link{col_concat}}
  \item \code{\link{within_bounds}}
  \item \code{\link{within_n_sds}}
  \item \code{\link{within_n_mads}}
  \item \code{\link{success_and_error_functions}}
  \item \code{\link{chaining_functions}}
  }
}
\examples{
library(magrittr)     # for the piping operator
library(dplyr)

# this confirms that
#   - that the dataset contains more than 10 observations
#   - that the column for 'miles per gallon' (mpg) is a positive number
#   - that the column for 'miles per gallon' (mpg) does not contain a datum
#     that is outside 4 standard deviations from its mean, and
#   - that the am and vs columns (automatic/manual and v/straight engine,
#     respectively) contain 0s and 1s only
#   - each row contains at most 2 NAs
#   - each row's mahalanobis distance is within 10 median absolute deviations of
#     all the distance (for outlier detection)

mtcars \%>\%
  verify(nrow(.) > 10) \%>\%
  verify(mpg > 0) \%>\%
  insist(within_n_sds(4), mpg) \%>\%
  assert(in_set(0,1), am, vs) \%>\%
  assert_rows(num_row_NAs, within_bounds(0,2), everything()) \%>\%
  insist_rows(maha_dist, within_n_mads(10), everything()) \%>\%
  group_by(cyl) \%>\%
  summarise(avg.mpg=mean(mpg))


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicates.R
\name{not_na}
\alias{not_na}
\title{Returns TRUE if value is not NA}
\usage{
not_na(x, allow.NaN = FALSE)
}
\arguments{
\item{x}{A R object that supports \link{is.na} an \link{is.nan}}

\item{allow.NaN}{A logical indicating whether NaNs should be allowed
(default FALSE)}
}
\value{
A vector of the same length that is TRUE when the element is
not NA and FALSE otherwise
}
\description{
This is the inverse of \code{\link[base]{is.na}}. This is a convenience
function meant to be used as a predicate in an \code{\link{assertr}}
assertion.
}
\examples{
not_na(NA)
not_na(2.8)
not_na("tree")
not_na(c(1, 2, NA, 4))

}
\seealso{
\code{\link{is.na}} \code{\link{is.nan}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicates.R
\name{within_n_mads}
\alias{within_n_mads}
\title{Return a function to create robust z-score checking predicate}
\usage{
within_n_mads(n, ...)
}
\arguments{
\item{n}{The number of median absolute deviations from the median
within which to accept a datum}

\item{...}{Additional arguments to be passed to \code{\link{within_bounds}}}
}
\value{
A function that takes a vector and returns a
        \code{\link{within_bounds}} predicate based on the MAD
        of that vector.
}
\description{
This function takes one argument, the number of median absolute
deviations within which to accept a particular data point. This is
generally more useful than its sister function \code{\link{within_n_sds}}
because it is more robust to the presence of outliers. It is therefore
better suited to identify potentially erroneous data points.
}
\details{
As an example, if '2' is passed into this function, this will return
a function that takes a vector and figures out the bounds of two
median absolute deviations (MADs) from the median. That function will then
return a \code{\link{within_bounds}} function that can then be applied
to a single datum. If the datum is within two MADs of the median of the
vector given to the function returned by this function, it will return TRUE.
If not, FALSE.

This function isn't meant to be used on its own, although it can. Rather,
this function is meant to be used with the \code{\link{insist}} function to
search for potentially erroneous data points in a data set.
}
\examples{
test.vector <- rnorm(100, mean=100, sd=20)

within.one.mad <- within_n_mads(1)
custom.bounds.checker <- within.one.mad(test.vector)
custom.bounds.checker(105)     # returns TRUE
custom.bounds.checker(40)      # returns FALSE

# same as
within_n_mads(1)(test.vector)(40)    # returns FALSE

within_n_mads(2)(test.vector)(as.numeric(NA))  # returns TRUE
# because, by default, within_bounds() will accept
# NA values. If we want to reject NAs, we have to
# provide extra arguments to this function
within_n_mads(2, allow.na=FALSE)(test.vector)(as.numeric(NA))  # returns FALSE

# or in a pipeline, like this was meant for

library(magrittr)

iris \%>\%
  insist(within_n_mads(5), Sepal.Length)

}
\seealso{
\code{\link{within_n_sds}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R, R/deprecated_se.R
\name{insist_rows}
\alias{insist_rows}
\alias{insist_rows_}
\title{Raises error if dynamically created predicate is FALSE for any row
after applying row reduction function}
\usage{
insist_rows(
  data,
  row_reduction_fn,
  predicate_generator,
  ...,
  success_fun = success_continue,
  error_fun = error_stop,
  skip_chain_opts = FALSE,
  obligatory = FALSE,
  defect_fun = defect_append,
  description = NA
)

insist_rows_(
  data,
  row_reduction_fn,
  predicate_generator,
  ...,
  .dots,
  success_fun = success_continue,
  error_fun = error_stop
)
}
\arguments{
\item{data}{A data frame}

\item{row_reduction_fn}{A function that returns a value for each row of
the provided data frame}

\item{predicate_generator}{A function that is applied to the results of
the row reduction function. This will produce,
a true predicate function to be applied to every
element in the vector that the row reduction
function returns.}

\item{...}{Comma separated list of unquoted expressions.
Uses dplyr's \code{select} to select
columns from data.}

\item{success_fun}{Function to call if assertion passes. Defaults to
returning \code{data}.}

\item{error_fun}{Function to call if assertion fails. Defaults to printing
a summary of all errors.}

\item{skip_chain_opts}{If TRUE, \code{success_fun} and \code{error_fun}
are used even if assertion is called within a chain.}

\item{obligatory}{If TRUE and assertion failed the data is marked as defective.
For defective data, all the following rules are handled by
\code{defect_fun} function.}

\item{defect_fun}{Function to call when data is defective. Defaults to skipping
assertion and storing info about it in special attribute.}

\item{description}{Custom description of the rule. Is stored in result
reports and data.}

\item{.dots}{Use insist_rows_() to select columns using standard evaluation.}
}
\value{
By default, the \code{data} is returned if dynamically created
        predicate assertion is TRUE and and error is thrown if not. If a
        non-default \code{success_fun} or \code{error_fun} is used, the
        return values of these function will be returned.
}
\description{
Meant for use in a data analysis pipeline, this function applies a
function to a data frame that reduces each row to a single value. Then,
a predicate generating function is applied to row reduction values. It will
then use these predicates to check each of the row reduction values. If any
of these predicate applications yield FALSE, this function will raise
an error, effectively terminating the pipeline early. If there are no
FALSEs, this function will just return the data that it was supplied for
further use in later parts of the pipeline.
}
\details{
For examples of possible choices for the \code{success_fun} and
\code{error_fun} parameters, run \code{help("success_and_error_functions")}
}
\note{
See \code{vignette("assertr")} for how to use this in context
}
\examples{

# returns mtcars
insist_rows(mtcars, maha_dist, within_n_mads(30), mpg:carb)

library(magrittr)                    # for piping operator

mtcars \%>\%
  insist_rows(maha_dist, within_n_mads(10), vs:am)
  # anything here will run

\dontrun{
mtcars \%>\%
  insist_rows(maha_dist, within_n_mads(1), everything())
  # the assertion is untrue so
  # nothing here will run}

}
\seealso{
\code{\link{insist}} \code{\link{assert_rows}}
         \code{\link{assert}} \code{\link{verify}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{summary.assertr_assert_error}
\alias{summary.assertr_assert_error}
\title{Summarizing assertr's assert errors}
\usage{
\method{summary}{assertr_assert_error}(object, ...)
}
\arguments{
\item{object}{An assertr_assert_error object}

\item{...}{Additional arguments affecting the summary produced}
}
\description{
`summary` method for class "assertr_assert_error"
This prints the error message and the first five
rows of the two-column `data.frame` holding the
indexes and values of the offending data.
}
\seealso{
\code{\link{print.assertr_assert_error}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{has_only_names}
\alias{has_only_names}
\title{Returns TRUE if data.frame or list has only the specified names}
\usage{
has_only_names(...)
}
\arguments{
\item{...}{A arbitrary amount of quoted names to check for}
}
\value{
TRUE is all names exist, FALSE if not
}
\description{
This function checks parent frame environment for a specific set of names; if
more columns are present than those specified, an error is raised.
}
\details{
This is meant to be used with `assertr`'s `verify` function to check
for the existence of specific column names in a `data.frame` that is
piped to `verify`. It can also work on a non-`data.frame` list.
}
\examples{

# The last two columns names are switched in order, but all column names are
# present, so it passes.
verify(
  mtcars,
  has_only_names(c(
    "mpg", "cyl", "disp", "hp", "drat", "wt", "qsec", "vs", "am",
    "carb", "gear"
  ))
)

# More than one set of character strings can be provided.
verify(
  mtcars,
  has_only_names(
    c("mpg", "cyl", "disp", "hp", "drat", "wt", "qsec", "vs", "am"),
    c("carb", "gear")
  )
)

\dontrun{
# The some columns are missing, so it fails.
verify(mtcars, has_only_names("mpg"))
}
}
\seealso{
Other Name verification: 
\code{\link{has_all_names}()}
}
\concept{Name verification}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R, R/deprecated_se.R
\name{assert}
\alias{assert}
\alias{assert_}
\title{Raises error if predicate is FALSE in any columns selected}
\usage{
assert(
  data,
  predicate,
  ...,
  success_fun = success_continue,
  error_fun = error_stop,
  skip_chain_opts = FALSE,
  obligatory = FALSE,
  defect_fun = defect_append,
  description = NA
)

assert_(
  data,
  predicate,
  ...,
  .dots,
  success_fun = success_continue,
  error_fun = error_stop
)
}
\arguments{
\item{data}{A data frame}

\item{predicate}{A function that returns FALSE when violated}

\item{...}{Comma separated list of unquoted expressions.
Uses dplyr's \code{select} to select
columns from data.}

\item{success_fun}{Function to call if assertion passes. Defaults to
returning \code{data}.}

\item{error_fun}{Function to call if assertion fails. Defaults to printing
a summary of all errors.}

\item{skip_chain_opts}{If TRUE, \code{success_fun} and \code{error_fun}
are used even if assertion is called within a chain.}

\item{obligatory}{If TRUE and assertion failed the data is marked as defective.
For defective data, all the following rules are handled by
\code{defect_fun} function.}

\item{defect_fun}{Function to call when data is defective. Defaults to skipping
assertion and storing info about it in special attribute.}

\item{description}{Custom description of the rule. Is stored in result
reports and data.}

\item{.dots}{Use assert_() to select columns using standard evaluation.}
}
\value{
By default, the \code{data} is returned if predicate assertion
        is TRUE and and error is thrown if not. If a non-default
        \code{success_fun} or \code{error_fun} is used, the return
        values of these function will be returned.
}
\description{
Meant for use in a data analysis pipeline, this function will
just return the data it's supplied if there are no FALSEs
when the predicate is applied to every element of the columns
indicated. If any element in any of the columns, when applied
to the predicate, is FALSE, then this function will raise an
error, effectively terminating the pipeline early.
}
\details{
For examples of possible choices for the \code{success_fun} and
\code{error_fun} parameters, run \code{help("success_and_error_functions")}
}
\note{
See \code{vignette("assertr")} for how to use this in context
}
\examples{

# returns mtcars
assert(mtcars, not_na, vs)

# return mtcars
assert(mtcars, not_na, mpg:carb)

library(magrittr)                    # for piping operator

mtcars \%>\%
  assert(in_set(c(0,1)), vs)
  # anything here will run

\dontrun{
mtcars \%>\%
  assert(in_set(c(1, 2, 3, 4, 6)), carb)
  # the assertion is untrue so
  # nothing here will run}

}
\seealso{
\code{\link{verify}} \code{\link{insist}}
         \code{\link{assert_rows}} \code{\link{insist_rows}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{print.assertr_assert_error}
\alias{print.assertr_assert_error}
\title{Printing assertr's assert errors}
\usage{
\method{print}{assertr_assert_error}(x, ...)
}
\arguments{
\item{x}{An assertr_assert_error object}

\item{...}{Further arguments passed to or from other methods}
}
\description{
`print` method for class "assertr_assert_error"
This prints the error message and the entire two-column
`data.frame` holding the indexes and values of the offending
data.
}
\seealso{
\code{\link{summary.assertr_assert_error}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assertions.R
\name{verify}
\alias{verify}
\title{Raises error if expression is FALSE anywhere}
\usage{
verify(
  data,
  expr,
  success_fun = success_continue,
  error_fun = error_stop,
  skip_chain_opts = FALSE,
  obligatory = FALSE,
  defect_fun = defect_append,
  description = NA
)
}
\arguments{
\item{data}{A data frame, list, or environment}

\item{expr}{A logical expression}

\item{success_fun}{Function to call if assertion passes. Defaults to
returning \code{data}.}

\item{error_fun}{Function to call if assertion fails. Defaults to printing
a summary of all errors.}

\item{skip_chain_opts}{If TRUE, \code{success_fun} and \code{error_fun}
are used even if assertion is called within a chain.}

\item{obligatory}{If TRUE and assertion failed the data is marked as defective.
For defective data, all the following rules are handled by
\code{defect_fun} function.}

\item{defect_fun}{Function to call when data is defective. Defaults to skipping
assertion and storing info about it in special attribute.}

\item{description}{Custom description of the rule. Is stored in result
reports and data.}
}
\value{
By default, the \code{data} is returned if predicate assertion
        is TRUE and and error is thrown if not. If a non-default
        \code{success_fun} or \code{error_fun} is used, the return
        values of these function will be returned.
}
\description{
Meant for use in a data analysis pipeline, this function will
just return the data it's supplied if all the logicals in the
expression supplied are TRUE. If at least one is FALSE, this
function will raise a error, effectively terminating the pipeline
early
}
\details{
For examples of possible choices for the \code{success_fun} and
\code{error_fun} parameters, run \code{help("success_and_error_functions")}
}
\note{
See \code{vignette("assertr")} for how to use this in context
}
\examples{

verify(mtcars, drat > 2)     # returns mtcars
\dontrun{
verify(mtcars, drat > 3)     # produces error}


library(magrittr)            # for piping operator

\dontrun{
mtcars \%>\%
  verify(drat > 3) \%>\%
  # anything here will not run}

mtcars \%>\%
  verify(nrow(mtcars) > 2)
  # anything here will run

alist <- list(a=c(1,2,3), b=c(4,5,6))
verify(alist, length(a) > 2)
verify(alist, length(a) > 2 && length(b) > 2)
verify(alist, a > 0 & b > 2)

\dontrun{
alist \%>\%
  verify(alist, length(a) > 5)
  # nothing here will run}


}
\seealso{
\code{\link{assert}} \code{\link{insist}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/row-redux.R
\name{duplicates_across_cols}
\alias{duplicates_across_cols}
\title{Checks if row contains at least one value duplicated in its column}
\usage{
duplicates_across_cols(data, allow.na = FALSE)
}
\arguments{
\item{data}{A data frame}

\item{allow.na}{TRUE if we allow NAs in data. Default FALSE.}
}
\value{
A logical vector.
}
\description{
This function will return a vector, with the same length as the number
of rows of the provided data frame. Each element of the vector will be
logical value that states if any value from the row was duplicated in
its column.
}
\examples{

df <- data.frame(v1 = c(1, 1, 2, 3), v2 = c(4, 5, 5, 6))
duplicates_across_cols(df)

library(magrittr)            # for piping operator

# you can use "assert_rows", "in_set", and this function to
# check if specified variables set and all subsets are keys for the data.

correct_df <- data.frame(id = 1:5, sub_id = letters[1:5], work_id = LETTERS[1:5])
correct_df \%>\%
  assert_rows(duplicates_across_cols, in_set(FALSE), id, sub_id, work_id)
  # passes because each subset of correct_df variables is key

\dontrun{
incorrect_df <- data.frame(id = 1:5, sub_id = letters[1:5], age = c(10, 20, 20, 15, 30))
incorrect_df \%>\%
  assert_rows(duplicates_across_cols, in_set(FALSE), id, sub_id, age)
  # fails because age is not key of the data (age == 20 is placed twice)
}

}
\seealso{
\code{\link{paste}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/row-redux.R
\name{col_concat}
\alias{col_concat}
\title{Concatenate all columns of each row in data frame into a string}
\usage{
col_concat(data, sep = "")
}
\arguments{
\item{data}{A data frame}

\item{sep}{A string to separate the columns with (default: "")}
}
\value{
A vector of rows concatenated into strings
}
\description{
This function will return a vector, with the same length as the number
of rows of the provided data frame. Each element of the vector will be
it's corresponding row with all of its values (one for each column)
"pasted" together in a string.
}
\examples{

col_concat(mtcars)

library(magrittr)            # for piping operator

# you can use "assert_rows", "is_uniq", and this function to
# check if joint duplicates (across different columns) appear
# in a data frame
\dontrun{
mtcars \%>\%
  assert_rows(col_concat, is_uniq, mpg, hp)
  # fails because the first two rows are jointly duplicates
  # on these two columns
}

mtcars \%>\%
  assert_rows(col_concat, is_uniq, mpg, hp, wt) # ok

}
\seealso{
\code{\link{paste}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{print.assertr_defect}
\alias{print.assertr_defect}
\title{Printing assertr's defect}
\usage{
\method{print}{assertr_defect}(x, ...)
}
\arguments{
\item{x}{An assertr_defect object}

\item{...}{Further arguments passed to or from other methods}
}
\description{
`print` method for class "assertr_defect"
This prints the defect message along with columns that were checked.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{success_and_error_functions}
\alias{success_and_error_functions}
\alias{success_logical}
\alias{success_continue}
\alias{success_append}
\alias{success_report}
\alias{success_df_return}
\alias{error_stop}
\alias{just_warn}
\alias{error_report}
\alias{warn_report}
\alias{error_append}
\alias{warning_append}
\alias{error_return}
\alias{error_df_return}
\alias{error_logical}
\alias{defect_append}
\alias{defect_report}
\alias{defect_df_return}
\title{Success and error functions}
\usage{
success_logical(data, ...)

success_continue(data, ...)

success_append(data, ...)

success_report(data, ...)

success_df_return(data, ...)

error_stop(errors, data = NULL, warn = FALSE, ...)

just_warn(errors, data = NULL)

error_report(errors, data = NULL, warn = FALSE, ...)

warn_report(errors, data = NULL)

error_append(errors, data = NULL)

warning_append(errors, data = NULL)

error_return(errors, data = NULL)

error_df_return(errors, data = NULL)

error_logical(errors, data = NULL, ...)

defect_append(errors, data, ...)

defect_report(errors, data, ...)

defect_df_return(errors, data, ...)
}
\arguments{
\item{data}{A data frame}

\item{...}{Further arguments passed to or from other methods}

\item{errors}{A list of objects of class \code{assertr_errors}}

\item{warn}{If TRUE, assertr will issue a warning instead of an error}
}
\description{
The behavior of functions like \code{assert}, \code{assert_rows},
\code{insist}, \code{insist_rows}, \code{verify} when the assertion
passes or fails is configurable via the \code{success_fun}
and \code{error_fun} parameters, respectively.
The \code{success_fun} parameter takes a function that takes
the data passed to the assertion function as a parameter. You can
write your own success handler function, but there are a few
provided by this package:
\itemize{
  \item \code{success_continue} - just returns the data that was
                                   passed into the assertion function
  \item \code{success_logical} - returns TRUE
  \item \code{success_append} - returns the data that was
                                passed into the assertion function
                                but also stores basic information about
                                verification result
  \item \code{success_report} - When success results are stored, and each
                                verification ended up with success prints
                                summary of all successful validations
  \item \code{success_df_return} - When success results are stored, and each
                                   verification ended up with success prints
                                   data.frame with verification results
}
The \code{error_fun} parameter takes a function that takes
the data passed to the assertion function as a parameter. You can
write your own error handler function, but there are a few
provided by this package:
\itemize{
  \item \code{error_stop} - Prints a summary of the errors and
                            halts execution.
  \item \code{error_report} - Prints all the information available
                              about the errors in a "tidy"
                              \code{data.frame} (including information
                              such as the name of the predicate used,
                              the offending value, etc...) and halts
                              execution.
  \item \code{error_append} - Attaches the errors to a special
   attribute of \code{data} and returns the data. This is chiefly
   to allow assertr errors to be accumulated in a pipeline so that
   all assertions can have a chance to be checked and so that all
   the errors can be displayed at the end of the chain.
  \item \code{error_return} - Returns the raw object containing all
    the errors
  \item \code{error_df_return} - Returns a "tidy" \code{data.frame}
    containing all the errors, including informations such as
    the name of the predicate used, the offending value, etc...
  \item \code{error_logical} - returns FALSE
  \item \code{just_warn} - Prints a summary of the errors but does
   not halt execution, it just issues a warning.
  \item \code{warn_report} - Prints all the information available
  about the errors but does not halt execution, it just issues a warning.
  \item \code{defect_report} - For single rule and defective data it displays
  short info about skipping current assertion. For \code{chain_end} sums
  up all skipped rules for defective data.
  \item \code{defect_df_return} - For single rule and defective data it returns
  info data.frame about skipping current assertion. For \code{chain_end}
  returns all skipped rules info data.frame for defective data.
 }
You may find the third type of data verification result. In a scenario
when validation rule was obligatory (obligatory = TRUE) in order to execute the
following ones we may want to skip them and register that fact.
In order to do this there are three callbacks reacting to defective
data:
 \itemize{
  \item \code{defect_report} - For single rule and defective data it displays
  short info about skipping current assertion.
  \item \code{defect_df_return} - For single rule and defective data it returns
  info data.frame about skipping current assertion.
  \item \code{defect_append} - Appends info about skipped rule due to data
  defect into one of data attributes. Rules skipped on defective data, or its summary, can
  be returned with proper error_fun callback in \code{chain_end}.
 }
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/row-redux.R
\name{maha_dist}
\alias{maha_dist}
\title{Computes mahalanobis distance for each row of data frame}
\usage{
maha_dist(data, keep.NA = TRUE, robust = FALSE, stringsAsFactors = FALSE)
}
\arguments{
\item{data}{A data frame}

\item{keep.NA}{Ensure that every row with missing data remains NA in
the output? TRUE by default.}

\item{robust}{Attempt to compute mahalanobis distance based on
robust covariance matrix? FALSE by default}

\item{stringsAsFactors}{Convert non-factor string columns into factors?
FALSE by default}
}
\value{
A vector of observation-wise mahalanobis distances.
}
\description{
This function will return a vector, with the same length as the number
of rows of the provided data frame, corresponding to the average
mahalanobis distances of each row from the whole data set.
}
\details{
This is useful for finding anomalous observations, row-wise.

It will convert any categorical variables in the data frame into numerics
as long as they are factors. For example, in order for a character
column to be used as a component in the distance calculations, it must
either be a factor, or converted to a factor by using the
\code{stringsAsFactors} parameter.
}
\examples{

maha_dist(mtcars)

maha_dist(iris, robust=TRUE)


library(magrittr)            # for piping operator
library(dplyr)               # for "everything()" function

# using every column from mtcars, compute mahalanobis distance
# for each observation, and ensure that each distance is within 10
# median absolute deviations from the median
mtcars \%>\%
  insist_rows(maha_dist, within_n_mads(10), everything())
  ## anything here will run

}
\seealso{
\code{\link{insist_rows}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicates.R
\name{in_set}
\alias{in_set}
\title{Returns TRUE if value in set}
\usage{
in_set(..., allow.na = TRUE, inverse = FALSE)
}
\arguments{
\item{...}{objects that make up the set}

\item{allow.na}{A logical indicating whether NAs (including NaNs)
should be permitted (default TRUE)}

\item{inverse}{A logical indicating whether it should test
if arguments are NOT in the set}
}
\value{
A function that takes one value and returns TRUE
        if the value is in the set defined by the
        arguments supplied by \code{in_set} and FALSE
        otherwise
}
\description{
This function returns a predicate function that will take a single
value and return TRUE if the value is a member of the set of objects
supplied. This doesn't actually check the membership of anything--it
only returns a function that actually does the checking when called
with a value. This is a convenience function meant to return a
predicate function to be used in an \code{\link{assertr}} assertion.
You can use the `inverse` flag (default FALSE) to check if the
arguments are NOT in the set.
}
\examples{
predicate <- in_set(3,4)
predicate(4)

## is equivalent to

in_set(3,4)(3)

# inverting the function works thusly...
in_set(3, 4, inverse=TRUE)(c(5, 2, 3))
# TRUE TRUE FALSE

# the remainder of division by 2 is always 0 or 1
rem <- 10 \%\% 2
in_set(0,1)(rem)

## this is meant to be used as a predicate in an assert statement
assert(mtcars, in_set(3,4,5), gear)

## or in a pipeline, like this was meant for

library(magrittr)

mtcars \%>\%
  assert(in_set(3,4,5), gear) \%>\%
  assert(in_set(0,1), vs, am)

}
\seealso{
\code{\link{\%in\%}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{generate_id}
\alias{generate_id}
\title{Generates random ID string}
\usage{
generate_id()
}
\description{
This is used to generate id for each assertion error.
}
\details{
For single assertion that checks multiple columns, each error log
is stored as a separate element. We provide the ID to allow
detecting which errors come from the same assertion.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicates.R
\name{within_bounds}
\alias{within_bounds}
\title{Creates bounds checking predicate}
\usage{
within_bounds(
  lower.bound,
  upper.bound,
  include.lower = TRUE,
  include.upper = TRUE,
  allow.na = TRUE
)
}
\arguments{
\item{lower.bound}{The lowest permitted value}

\item{upper.bound}{The upper permitted value}

\item{include.lower}{A logical indicating whether lower bound
should be inclusive (default TRUE)}

\item{include.upper}{A logical indicating whether upprt bound
should be inclusive (default TRUE)}

\item{allow.na}{A logical indicating whether NAs (including NaNs)
should be permitted (default TRUE)}
}
\value{
A function that takes numeric value or numeric vactor and returns
        TRUE if the value(s) is/are within the bounds defined by the
        arguments supplied by \code{within_bounds} and FALSE
        otherwise
}
\description{
This function returns a predicate function that will take a numeric value
or vector and return TRUE if the value(s) is/are within the bounds set.
This does not actually check the bounds of anything--it only returns
a function that actually does the checking when called with a number.
This is a convenience function meant to return a predicate function to
be used in an \code{\link{assertr}} assertion.
}
\examples{
predicate <- within_bounds(3,4)
predicate(pi)

## is equivalent to

within_bounds(3,4)(pi)

# a correlation coefficient must always be between 0 and 1
coeff <- cor.test(c(1,2,3), c(.5, 2.4, 4))[["estimate"]]
within_bounds(0,1)(coeff)

## check for positive number
positivep <- within_bounds(0, Inf, include.lower=FALSE)

## this is meant to be used as a predicate in an assert statement
assert(mtcars, within_bounds(4,8), cyl)

## or in a pipeline

library(magrittr)

mtcars \%>\%
  assert(within_bounds(4,8), cyl)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicates.R
\name{is_uniq}
\alias{is_uniq}
\title{Returns TRUE where no elements appear more than once}
\usage{
is_uniq(..., allow.na = FALSE)
}
\arguments{
\item{...}{One or more vectors to check for unique combinations of elements}

\item{allow.na}{A logical indicating whether NAs should be preserved
as missing values in the return value (FALSE) or
if they should be treated just like any other value
(TRUE) (default is FALSE)}
}
\value{
A vector of the same length where the corresponding element
        is TRUE if the element only appears once in the vector and
        FALSE otherwise
}
\description{
This function is meant to take only a vector. It relies heavily on
the \code{\link{duplicated}} function where it can be thought of as
the inverse. Where this function differs, though--besides being only
meant for one vector or column--is that it marks the first occurrence
of a duplicated value as "non unique", as well.
}
\examples{

is_uniq(1:10)
is_uniq(c(1,1,2,3), c(1,2,2,3))

\dontrun{
# returns FALSE where a "5" appears
is_uniq(c(1:10, 5))
}

library(magrittr)

\dontrun{
# this fails 4 times
mtcars \%>\% assert(is_uniq, qsec)
}

# to use the version of this function that allows NAs in `assert`,
you can use a lambda/anonymous function like so:

mtcars \%>\%
  assert(function(x){is_uniq(x, allow.na=TRUE)}, qsec)

}
\seealso{
\code{\link{duplicated}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicates.R
\name{within_n_sds}
\alias{within_n_sds}
\title{Return a function to create z-score checking predicate}
\usage{
within_n_sds(n, ...)
}
\arguments{
\item{n}{The number of standard deviations from the mean
within which to accept a datum}

\item{...}{Additional arguments to be passed to \code{\link{within_bounds}}}
}
\value{
A function that takes a vector and returns a
        \code{\link{within_bounds}} predicate based on the standard deviation
        of that vector.
}
\description{
This function takes one argument, the number of standard deviations
within which to accept a particular data point.
}
\details{
As an example, if '2' is passed into this function, this will return
a function that takes a vector and figures out the bounds of two
standard deviations from the mean. That function will then return
a \code{\link{within_bounds}} function that can then be applied
to a single datum. If the datum is within two standard deviations of
the mean of the vector given to the function returned by this function,
it will return TRUE. If not, FALSE.

This function isn't meant to be used on its own, although it can. Rather,
this function is meant to be used with the \code{\link{insist}} function to
search for potentially erroneous data points in a data set.
}
\examples{
test.vector <- rnorm(100, mean=100, sd=20)

within.one.sd <- within_n_sds(1)
custom.bounds.checker <- within.one.sd(test.vector)
custom.bounds.checker(105)     # returns TRUE
custom.bounds.checker(40)      # returns FALSE

# same as
within_n_sds(1)(test.vector)(40)    # returns FALSE

within_n_sds(2)(test.vector)(as.numeric(NA))  # returns TRUE
# because, by default, within_bounds() will accept
# NA values. If we want to reject NAs, we have to
# provide extra arguments to this function
within_n_sds(2, allow.na=FALSE)(test.vector)(as.numeric(NA))  # returns FALSE

# or in a pipeline, like this was meant for

library(magrittr)

iris \%>\%
  insist(within_n_sds(5), Sepal.Length)

}
\seealso{
\code{\link{within_n_mads}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{print.assertr_success}
\alias{print.assertr_success}
\title{Printing assertr's success}
\usage{
\method{print}{assertr_success}(x, ...)
}
\arguments{
\item{x}{An assertr_success object}

\item{...}{Further arguments passed to or from other methods}
}
\description{
`print` method for class "assertr_success"
This prints the success message along with columns that were checked.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{print.assertr_verify_error}
\alias{print.assertr_verify_error}
\title{Printing assertr's verify errors}
\usage{
\method{print}{assertr_verify_error}(x, ...)
}
\arguments{
\item{x}{An assertr_verify_error object.}

\item{...}{Further arguments passed to or from other methods}
}
\description{
`summary` method for class "assertr_verify_error"
}
\seealso{
\code{\link{summary.assertr_verify_error}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{has_all_names}
\alias{has_all_names}
\title{Returns TRUE if data.frame or list has specified names}
\usage{
has_all_names(...)
}
\arguments{
\item{...}{A arbitrary amount of quoted names to check for}
}
\value{
TRUE if all names exist, FALSE if not
}
\description{
This function checks parent frame environment for existence of names.
This is meant to be used with `assertr`'s `verify` function to check
for the existence of specific column names in a `data.frame` that is
piped to `verify`. It can also work on a non-`data.frame` list.
}
\examples{

verify(mtcars, has_all_names("mpg", "wt", "qsec"))

library(magrittr)   # for pipe operator

\dontrun{
mtcars \%>\%
  verify(has_all_names("mpgg"))  # fails
}

mpgg <- "something"

mtcars \%>\%
  verify(exists("mpgg"))   # passes but big mistake

\dontrun{
mtcars \%>\%
  verify(has_all_names("mpgg")) # correctly fails
}

}
\seealso{
\code{\link{exists}}

Other Name verification: 
\code{\link{has_only_names}()}
}
\concept{Name verification}
