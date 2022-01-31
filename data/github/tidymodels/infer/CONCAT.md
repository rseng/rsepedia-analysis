# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement listed as package
authors. All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/
code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.

# infer R Package <img src="figs/infer_gnome.png" align="right" width=280 />

<!--figs/infer.svg-->
<!--http://www.r-pkg.org/badges/version/infer-->
<!--figs/main.svg-->
<!--https://img.shields.io/codecov/c/github/tidymodels/infer/main.svg-->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/infer)](https://cran.r-project.org/package=infer)[![Coverage
Status](https://img.shields.io/codecov/c/github/tidymodels/infer/main.svg)](https://codecov.io/github/tidymodels/infer/?branch=main)

The objective of this package is to perform statistical inference using
an expressive statistical grammar that coheres with the `tidyverse`
design framework. The package is centered around 4 main verbs,
supplemented with many utilities to visualize and extract value from
their outputs.

-   `specify()` allows you to specify the variable, or relationship
    between variables, that you’re interested in.
-   `hypothesize()` allows you to declare the null hypothesis.
-   `generate()` allows you to generate data reflecting the null
    hypothesis.
-   `calculate()` allows you to calculate a distribution of statistics
    from the generated data to form the null distribution.

To learn more about the principles underlying the package design, see
`vignette("infer")`.

<div class="figure">

<img src="https://raw.githubusercontent.com/tidymodels/infer/main/figs/ht-diagram.png" alt="A diagram showing four steps to carry out randomization-based inference: specify hypothesis, generate data, calculate statistic, and visualize. From left to right, each step is connected by an arrow, while the diagram indicates that generating data and calculating statistics can happen iteratively."  />
<p class="caption">
</p>

</div>

If you’re interested in learning more about randomization-based
statistical inference generally, including applied examples of this
package, we recommend checking out [Statistical Inference Via Data
Science: A ModernDive Into R and the Tidyverse](https://moderndive.com/)
and [Introduction to Modern
Statistics](https://openintro-ims.netlify.app/).

### Installation

------------------------------------------------------------------------

To install the current stable version of `infer` from CRAN:

``` r
install.packages("infer")
```

To install the developmental stable version of `infer`, make sure to
install `remotes` first. The `pkgdown` website for this version is at
[infer.tidymodels.org](https://infer.tidymodels.org/).

``` r
install.packages("remotes")
remotes::install_github("tidymodels/infer")
```

### Contributing

------------------------------------------------------------------------

We welcome others helping us make this package as user-friendly and
efficient as possible. Please review our
[contributing](https://github.com/tidymodels/infer/blob/main/CONTRIBUTING.md)
and [conduct](https://github.com/tidymodels/infer/blob/main/CONDUCT.md)
guidelines. By participating in this project you agree to abide by its
terms.

For questions and discussions about tidymodels packages, modeling, and
machine learning, please [post on RStudio
Community](https://rstd.io/tidymodels-community). If you think you have
encountered a bug, please [submit an
issue](https://github.com/tidymodels/infer/issues). Either way, learn
how to create and share a [reprex](https://rstd.io/reprex) (a minimal,
reproducible example), to clearly communicate about your code. Check out
further details on [contributing guidelines for tidymodels
packages](https://www.tidymodels.org/contribute/) and [how to get
help](https://www.tidymodels.org/help/).

### Examples

------------------------------------------------------------------------

These examples are pulled from the “Full infer Pipeline Examples”
vignette, accessible by calling `vignette("observed_stat_examples")`.
They make use of the `gss` dataset supplied by the package, providing a
sample of data from the [General Social Survey](https://gss.norc.org).
The data looks like this:

``` r
# load in the dataset
data(gss)

# take a glimpse at it
str(gss)
```

    ## tibble [500 × 11] (S3: tbl_df/tbl/data.frame)
    ##  $ year   : num [1:500] 2014 1994 1998 1996 1994 ...
    ##  $ age    : num [1:500] 36 34 24 42 31 32 48 36 30 33 ...
    ##  $ sex    : Factor w/ 2 levels "male","female": 1 2 1 1 1 2 2 2 2 2 ...
    ##  $ college: Factor w/ 2 levels "no degree","degree": 2 1 2 1 2 1 1 2 2 1 ...
    ##  $ partyid: Factor w/ 5 levels "dem","ind","rep",..: 2 3 2 2 3 3 1 2 3 1 ...
    ##  $ hompop : num [1:500] 3 4 1 4 2 4 2 1 5 2 ...
    ##  $ hours  : num [1:500] 50 31 40 40 40 53 32 20 40 40 ...
    ##  $ income : Ord.factor w/ 12 levels "lt $1000"<"$1000 to 2999"<..: 12 11 12 12 12 12 12 12 12 10 ...
    ##  $ class  : Factor w/ 6 levels "lower class",..: 3 2 2 2 3 3 2 3 3 2 ...
    ##  $ finrela: Factor w/ 6 levels "far below average",..: 2 2 2 4 4 3 2 4 3 1 ...
    ##  $ weight : num [1:500] 0.896 1.083 0.55 1.086 1.083 ...

As an example, we’ll run an analysis of variance on `age` and `partyid`,
testing whether the age of a respondent is independent of their
political party affiliation.

Calculating the observed statistic,

``` r
F_hat <- gss %>% 
  specify(age ~ partyid) %>%
  calculate(stat = "F")
```

Then, generating the null distribution,

``` r
null_dist <- gss %>%
   specify(age ~ partyid) %>%
   hypothesize(null = "independence") %>%
   generate(reps = 1000, type = "permute") %>%
   calculate(stat = "F")
```

Visualizing the observed statistic alongside the null distribution,

``` r
visualize(null_dist) +
  shade_p_value(obs_stat = F_hat, direction = "greater")
```

<div class="figure">

<img src="https://raw.githubusercontent.com/tidymodels/infer/main/README_files/figure-gfm/viz-1.png" alt="A histogram showing a distribution of F statistics, right-tailed and centered around one. The x axis ranges from zero to five. The region of the histogram to the right of the observed statistic, just above two, is shaded red to represent the p-value."  />
<p class="caption">
</p>

</div>

Calculating the p-value from the null distribution and observed
statistic,

``` r
null_dist %>%
  get_p_value(obs_stat = F_hat, direction = "greater")
```

    ## # A tibble: 1 × 1
    ##   p_value
    ##     <dbl>
    ## 1   0.055

Note that the formula and non-formula interfaces (i.e. `age ~ partyid`
vs. `response = age, explanatory =  partyid`) work for all implemented
inference procedures in `infer`. Use whatever is more natural for you.
If you will be doing modeling using functions like `lm()` and `glm()`,
though, we recommend you begin to use the formula `y ~ x` notation as
soon as possible.

Other resources are available in the package vignettes! See
`vignette("observed_stat_examples")` for more examples like the one
above, and `vignette("infer")` for discussion of the underlying
principles of the package design.
# infer v1.0.1.9000 (Development Version)

To be released as v1.0.2.


* Fix p-value shading when the calculated statistic falls exactly on the boundaries of a histogram bin (#424).
* Fix `generate()` errors when columns are named `x` (#431).
* Fix error from `visualize` when passed `generate()`d `infer_dist` objects that had not been passed to `hypothesize()` (#432). 
* Update visual checks for `visualize` output to align with the R 4.1.0+ graphics engine (#438).

# infer v1.0.1 (GitHub Only)

This release reflects the infer version accepted to the Journal of Open Source Software.

* Re-licensed the package from CC0 to MIT. See the `LICENSE` and `LICENSE.md` files.
* Contributed a paper to the Journal of Open Source Software, a draft of which is available in `/figs/paper`.
* Various improvements to documentation (#417, #418).

# infer 1.0.0

v1.0.0 is the first major release of the {infer} package! By and large, the core verbs `specify()`, `hypothesize()`, `generate()`, and `calculate()` will interface as they did before. This release makes several improvements to behavioral consistency of the package and introduces support for theory-based inference as well as randomization-based inference with multiple explanatory variables.

## Behavioral consistency

A major change to the package in this release is a set of standards for behavorial consistency of `calculate()` (#356). Namely, the package will now

* supply a consistent error when the supplied `stat` argument isn't well-defined
for the variables `specify()`d

``` r
gss %>%
  specify(response = hours) %>%
  calculate(stat = "diff in means")
#> Error: A difference in means is not well-defined for a 
#> numeric response variable (hours) and no explanatory variable.
```

or

``` r
gss %>%
  specify(college ~ partyid, success = "degree") %>%
  calculate(stat = "diff in props")
#> Error: A difference in proportions is not well-defined for a dichotomous categorical 
#> response variable (college) and a multinomial categorical explanatory variable (partyid).
```

* supply a consistent message when the user supplies unneeded information via `hypothesize()` to `calculate()` an observed statistic

``` r
# supply mu = 40 when it's not needed
gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  calculate(stat = "mean")
#> Message: The point null hypothesis `mu = 40` does not inform calculation of 
#> the observed statistic (a mean) and will be ignored.
#> # A tibble: 1 x 1
#>    stat
#>   <dbl>
#> 1  41.4
```

and

* supply a consistent warning and assume a reasonable null value when the user does not supply sufficient information to calculate an observed statistic

``` r
# don't hypothesize `p` when it's needed
gss %>%
    specify(response = sex, success = "female") %>%
    calculate(stat = "z")
#> # A tibble: 1 x 1
#>    stat
#>   <dbl>
#> 1 -1.16
#> Warning message:
#> A z statistic requires a null hypothesis to calculate the observed statistic. 
#> Output assumes the following null value: `p = .5`. 
```

or

``` r
# don't hypothesize `p` when it's needed
gss %>%
  specify(response = partyid) %>%
  calculate(stat = "Chisq")
#> # A tibble: 1 x 1
#>    stat
#>  <dbl>
#> 1  334.
#> Warning message:
#> A chi-square statistic requires a null hypothesis to calculate the observed statistic. 
#> Output assumes the following null values: `p = c(dem = 0.2, ind = 0.2, rep = 0.2, other = 0.2, DK = 0.2)`.
```

To accommodate this behavior, a number of new `calculate` methods were added or improved. Namely:

- Implemented the standardized proportion $z$ statistic for one categorical variable
- Extended `calculate()` with `stat = "t"` by passing `mu` to the `calculate()` method for `stat = "t"` to allow for calculation of `t` statistics for one numeric variable with hypothesized mean
- Extended `calculate()` to allow lowercase aliases for `stat` arguments (#373).
- Fixed bugs in `calculate()` for to allow for programmatic calculation of statistics

This behavorial consistency also allowed for the implementation of `observe()`, a wrapper function around `specify()`, `hypothesize()`, and `calculate()`, to calculate observed statistics. The function provides a shorthand alternative to calculating observed statistics from data:

``` r
# calculating the observed mean number of hours worked per week
gss %>%
  observe(hours ~ NULL, stat = "mean")
#> # A tibble: 1 x 1
#>    stat
#>   <dbl>
#> 1  41.4

# equivalently, calculating the same statistic with the core verbs
gss %>%
  specify(response = hours) %>%
  calculate(stat = "mean")
#> # A tibble: 1 x 1
#>    stat
#>   <dbl>
#> 1  41.4

# calculating a t statistic for hypothesized mu = 40 hours worked/week
gss %>%
  observe(hours ~ NULL, stat = "t", null = "point", mu = 40)
#> # A tibble: 1 x 1
#>    stat
#>   <dbl>
#> 1  2.09

# equivalently, calculating the same statistic with the core verbs
gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  calculate(stat = "t")
#> # A tibble: 1 x 1
#>    stat
#>   <dbl>
#> 1  2.09
```

We don't anticipate that these changes are "breaking" in the sense that code that previously worked will continue to, though it may now message or warn in a way that it did not used to or error with a different (and hopefully more informative) message.

## A framework for theoretical inference

This release also introduces a more complete and principled interface for theoretical inference. While the package previously supplied some methods for visualization of theory-based curves, the interface did not provide any object that was explicitly a "null distribution" that could be supplied to helper functions like `get_p_value()` and `get_confidence_interval()`. The new interface is based on a new verb, `assume()`, that returns a null distribution that can be interfaced in the same way that simulation-based null distributions can be interfaced with.

As an example, we'll work through a full infer pipeline for inference on a mean using infer's `gss` dataset. Supposed that we believe the true mean number of hours worked by Americans in the past week is 40.

First, calculating the observed `t`-statistic:

``` r
obs_stat <- gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  calculate(stat = "t")

obs_stat
#> Response: hours (numeric)
#> Null Hypothesis: point
#> # A tibble: 1 x 1
#>    stat
#>   <dbl>
#> 1  2.09
```

The code to define the null distribution is very similar to that required to calculate a theorized observed statistic, switching out `calculate()` for `assume()` and replacing arguments as needed.

``` r
null_dist <- gss %>%
  specify(response = hours) %>%
  assume(distribution = "t")

null_dist 
#> A T distribution with 499 degrees of freedom.
```

This null distribution can now be interfaced with in the same way as a simulation-based null distribution elsewhere in the package. For example, calculating a p-value by juxtaposing the observed statistic and null distribution:

``` r
get_p_value(null_dist, obs_stat, direction = "both")
#> # A tibble: 1 x 1
#>   p_value
#>     <dbl>
#> 1  0.0376
```

…or visualizing the null distribution alone:

``` r
visualize(null_dist)
```

![](https://i.imgur.com/g3B5coD.png)

…or juxtaposing the two visually:

``` r
visualize(null_dist) + 
  shade_p_value(obs_stat, direction = "both")
```

![](https://i.imgur.com/3C66kgK.png)

Confidence intervals lie in data space rather than the standardized scale of the theoretical distributions. Calculating a mean rather than the standardized `t`-statistic:

``` r
obs_mean <- gss %>%
  specify(response = hours) %>%
  calculate(stat = "mean")
```

The null distribution here just defines the spread for the standard error calculation.

``` r
ci <- 
  get_confidence_interval(
    null_dist,
    level = .95,
    point_estimate = obs_mean
  )

ci
#> # A tibble: 1 x 2
#>   lower_ci upper_ci
#>      <dbl>    <dbl>
#> 1     40.1     42.7
```

Visualizing the confidence interval results in the theoretical distribution being recentered and rescaled to align with the scale of the observed data:

``` r
visualize(null_dist) + 
  shade_confidence_interval(ci)
```

![](https://i.imgur.com/4akSCY3.png)

Previous methods for interfacing with theoretical distributions are superseded—they will continue to be supported, though documentation will forefront the `assume()` interface.

## Support for multiple regression

The 2016 "Guidelines for Assessment and Instruction in Statistics Education" [1] state that, in introductory statistics courses, "[s]tudents should gain experience with how statistical models, including multivariable models, are used." In line with this recommendation, we introduce support for randomization-based inference with multiple explanatory variables via a new `fit.infer` core verb.

If passed an `infer` object, the method will parse a formula out of the `formula` or `response` and `explanatory` arguments, and pass both it and `data` to a `stats::glm` call.

``` r
gss %>%
  specify(hours ~ age + college) %>%
  fit()
#> # A tibble: 3 x 2
#>   term          estimate
#>   <chr>            <dbl>
#> 1 intercept     40.6    
#> 2 age            0.00596
#> 3 collegedegree  1.53
```

Note that the function returns the model coefficients as `estimate` rather than their associated `t`-statistics as `stat`.

If passed a `generate()`d object, the model will be fitted to each replicate.

``` r
gss %>%
  specify(hours ~ age + college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 100, type = "permute") %>%
  fit()
#> # A tibble: 300 x 3
#> # Groups:   replicate [100]
#>    replicate term          estimate
#>        <int> <chr>            <dbl>
#>  1         1 intercept     44.4    
#>  2         1 age           -0.0767 
#>  3         1 collegedegree  0.121  
#>  4         2 intercept     41.8    
#>  5         2 age            0.00344
#>  6         2 collegedegree -1.59   
#>  7         3 intercept     38.3    
#>  8         3 age            0.0761 
#>  9         3 collegedegree  0.136  
#> 10         4 intercept     43.1    
#> # … with 290 more rows
```

If `type = "permute"`, a set of unquoted column names in the data to permute (independently of each other) can be passed via the `variables` argument to `generate`. It defaults to only the response variable.

``` r
gss %>%
  specify(hours ~ age + college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 100, type = "permute", variables = c(age, college)) %>%
  fit()
#> # A tibble: 300 x 3
#> # Groups:   replicate [100]
#>    replicate term          estimate
#>        <int> <chr>            <dbl>
#>  1         1 intercept      39.4   
#>  2         1 age             0.0748
#>  3         1 collegedegree  -2.98  
#>  4         2 intercept      42.8   
#>  5         2 age            -0.0190
#>  6         2 collegedegree  -1.83  
#>  7         3 intercept      40.4   
#>  8         3 age             0.0354
#>  9         3 collegedegree  -1.31  
#> 10         4 intercept      40.9   
#> # … with 290 more rows
```

This feature allows for more detailed exploration of the effect of disrupting the correlation structure among explanatory variables on outputted model coefficients.

Each of the auxillary functions `get_p_value()`, `get_confidence_interval()`, `visualize()`, `shade_p_value()`, and `shade_confidence_interval()` have methods to handle `fit()` output! See their help-files for example usage. Note that `shade_*` functions now delay evaluation until they are added to an existing ggplot (e.g. that outputted by `visualize()`) with `+`.

## Improvements

- Following extensive discussion, the `generate()` type `type = "simulate"` has been renamed to the more evocative `type = "draw"`. We will continue to support `type = "simulate"` indefinitely, though supplying that argument will now prompt a message notifying the user of its preferred alias. (#233, #390)
- Fixed several bugs related to factors with unused levels. `specify()` will now drop unused factor levels and message that it has done so. (#374, #375, #397, #380)
- Added `two.sided` as an acceptable alias for `two_sided` for the `direction` argument in `get_p_value()` and `shade_p_value()`. (#355)
- Various improvements to documentation, including extending example sections in help-files, re-organizing the function reference in the {pkgdown} site, and linking more extensively among help-files.

## Breaking changes

We don't anticipate that any changes made in this release are "breaking" in the sense that code that previously worked will continue to, though it may now message or warn in a way that it did not used to or error with a different (and hopefully more informative) message. If you currently teach or research with infer, we recommend re-running your materials and noting any changes in messaging and warning.

- Move forward with a number of planned deprecations. Namely, the `GENERATION_TYPES` object is now fully deprecated, and arguments that were relocated from `visualize()` to `shade_p_value()` and `shade_confidence_interval()` are now fully deprecated in `visualize()`. If supplied a deprecated argument, `visualize()` will warn the user and ignore the argument.
- Added a `prop` argument to `rep_slice_sample()` as an alternative to the `n`
argument for specifying the proportion of rows in the supplied data to sample
per replicate (#361, #362, #363). This changes order of arguments of
`rep_slice_sample()` (in order to be more aligned with `dplyr::slice_sample()`)
which might break code if it didn't use named arguments (like
`rep_slice_sample(df, 5, TRUE)`). To fix this, use named arguments (like
`rep_slice_sample(df, 5, replicate = TRUE)`).

## Other

- Added Simon P. Couch as an author. Long deserved for his reliable maintenance and improvements of the package.

[1]: GAISE College Report ASA Revision Committee, "Guidelines for Assessment and Instruction in Statistics Education College Report 2016," http://www.amstat.org/education/gaise. 

# infer 0.5.4

- `rep_sample_n()` no longer errors when supplied a `prob` argument (#279)
- Added `rep_slice_sample()`, a light wrapper around `rep_sample_n()`, that
more closely resembles `dplyr::slice_sample()` (the function that supersedes
`dplyr::sample_n()`) (#325)
- Added a `success`, `correct`, and `z` argument to `prop_test()` 
(#343, #347, #353)
- Implemented observed statistic calculation for the standardized proportion 
$z$ statistic (#351, #353)
- Various bug fixes and improvements to documentation and errors.

# infer 0.5.3

## Breaking changes

- `get_confidence_interval()` now uses column names ('lower_ci' and 'upper_ci') 
in output that are consistent with other infer functionality (#317).

## New functionality

- `get_confidence_interval()` can now produce bias-corrected confidence intervals
by setting `type = "bias-corrected"`. Thanks to @davidbaniadam for the 
initial implementation (#237, #318)!

## Other

- Fix CRAN check failures related to long double errors.

# infer 0.5.2

- Warn the user when a p-value of 0 is reported (#257, #273)
- Added new vignettes: `chi_squared` and `anova` (#268)
- Updates to documentation and existing vignettes (#268)
- Add alias for `hypothesize()` (`hypothesise()`) (#271)
- Subtraction order no longer required for difference-based tests--a warning will be raised in the case that the user doesn't supply an `order` argument (#275, #281)
- Add new messages for common errors (#277)
- Increase coverage of theoretical methods in documentation (#278, #280)
- Drop missing values and reduce size of `gss` dataset used in examples (#282)
- Add `stat = "ratio of props"` and `stat = "odds ratio"` to `calculate` (#285)
- Add `prop_test()`, a tidy interface to `prop.test()` (#284, #287)
- Updates to `visualize()` for compatibility with `ggplot2` v3.3.0 (#289)
- Fix error when bootstrapping with small samples and raise warnings/errors 
when appropriate (#239, #244, #291)
- Fix unit test failures resulting from breaking changes in `dplyr` v1.0.0
- Fix error in `generate()` when response variable is named `x` (#299)
- Add `two-sided` and `two sided` as aliases for `two_sided` for the 
`direction` argument in `get_p_value()` and `shade_p_value()` (#302)
- Fix `t_test()` and `t_stat()` ignoring the `order` argument (#310)

# infer 0.5.1

- Updates to documentation and other tweaks

# infer 0.5.0

## Breaking changes

- `shade_confidence_interval()` now plots vertical lines starting from zero (previously - from the bottom of a plot) (#234).
- `shade_p_value()` now uses "area under the curve" approach to shading (#229).

## Other

- Updated `chisq_test()` to take arguments in a response/explanatory format, perform goodness of fit tests, and default to the approximation approach (#241).
- Updated `chisq_stat()` to do goodness of fit (#241).
- Make interface to `hypothesize()` clearer by adding the options for the point null parameters to the function signature (#242).
- Manage `infer` class more systematically (#219).
- Use `vdiffr` for plot testing (#221).

# infer 0.4.1

- Added Evgeni Chasnovski as author for his incredible work on refactoring the package and providing excellent support.

# infer 0.4.0

## Breaking changes

- Changed method of computing two-sided p-value to a more conventional one. It also makes `get_pvalue()` and `visualize()` more aligned (#205).

## Deprecation changes

- Deprecated `p_value()` (use `get_p_value()` instead) (#180).
- Deprecated `conf_int()` (use `get_confidence_interval()` instead) (#180).
- Deprecated (via warnings) plotting p-value and confidence interval in `visualize()` (use new functions `shade_p_value()` and `shade_confidence_interval()` instead) (#178).

## New functions

- `shade_p_value()` - {ggplot2}-like layer function to add information about p-value region to `visualize()` output. Has alias `shade_pvalue()`.
- `shade_confidence_interval()` - {ggplot2}-like layer function to add information about confidence interval region to `visualize()` output. Has alias `shade_ci()`.

## Other

- Account for `NULL` value in left hand side of formula in `specify()` (#156) and `type` in `generate()` (#157).
- Update documentation code to follow tidyverse style guide (#159).
- Remove help page for internal `set_params()` (#165).
- Fully use {tibble} (#166).
- Fix `calculate()` to not depend on order of `p` for `type = "simulate"` (#122).
- Reduce code duplication (#173).
- Make transparency in `visualize()` to not depend on method and data volume.
- Make `visualize()` work for "One sample t" theoretical type with `method = "both"`.
- Add `stat = "sum"` and `stat = "count"` options to `calculate()` (#50).

# infer 0.3.1

- Stop using package {assertive} in favor of custom type checks (#149)
- Fixed `t_stat()` to use `...` so `var.equal` works
- With the help of @echasnovski, fixed `var.equal = TRUE` for `specify() %>% calculate(stat = "t")`
- Use custom functions for error, warning, message, and `paste()` handling (#155)

# infer 0.3.0

- Added `conf_int` logical argument and `conf_level` argument to `t_test()`
- Switched `shade_color` argument in `visualize()` to be `pvalue_fill` instead
since fill color for confidence intervals is also added now
- Shading for Confidence Intervals in `visualize()` 
    - Green is default color for CI and red for p-values
    - `direction = "between"` to get the green shading
    - Currently working only for simulation-based methods
- Implemented `conf_int()` function for computing confidence interval provided a simulation-based method with a `stat` variable
    - `get_ci()` and `get_confidence_interval()` are aliases for `conf_int()`
    - Converted longer confidence interval calculation code in vignettes to use `get_ci()` instead    
- Implemented `p_value()` function for computing p-value provided a simulation-based method with a `stat` variable
    - `get_pvalue()` is an alias for `p_value()`
    - Converted longer p-value calculation code in vignettes to use `get_pvalue()` instead
- Implemented Chi-square Goodness of Fit observed stat depending on `params` being set in `hypothesize` with `specify() %>% calculate()` shortcut
- Removed "standardized" slope $t$ since its formula is different than "standardized" correlation and there is no way currently to give one over the other
- Implemented correlation with bootstrap CI and permutation hypothesis test
- Filled the `type` argument automatically in `generate()` based
on `specify()` and `hypothesize()`
    - Added message if `type` is given differently than expected
- Implemented `specify() %>% calculate()` for getting observed
statistics.
    - `visualize()` works with either a 1x1 data frame or a vector
    for its `obs_stat` argument
    - Got `stat = "t"` working
- Refactored `calculate()` into smaller functions to reduce complexity
- Produced error if `mu` is given in `hypothesize()` but `stat = "median"`
is provided in `calculate()` and other similar mis-specifications
- Tweaked `chisq_stat()` and `t_stat()` to match with `specify() %>% calculate()` framework
    - Both work in the one sample and two sample cases by providing `formula`
    - Added `order` argument to `t_stat()`
- Added implementation of one sample `t_test()` by passing in the `mu` argument to `t.test`
from `hypothesize()`
- Tweaked `pkgdown` page to include ToDo's using [{dplyr}](https://github.com/tidyverse/dplyr/blob/master/_pkgdown.yml) example

# infer 0.2.0

- Switched to `!!` instead of `UQ()` since `UQ()` is deprecated in 
{rlang} 0.2.0
- Added many new files: `CONDUCT.md`, `CONTRIBUTING.md`, and `TO-DO.md`
- Updated README file with more development information
- Added wrapper functions `t_test()` and `chisq_test()` that use a
formula interface and provide an intuitive wrapper to `t.test()` and
`chisq.test()`
- Created `stat = "z"` and `stat = "t"` options
- Added many new arguments to `visualize()` to prescribe colors to shade and 
use for observed statistics and theoretical density curves
- Added check so that a bar graph created with `visualize()` if number of 
unique values for generated statistics is small
- Added shading for `method = "theoretical"` 
- Implemented shading for simulation methods w/o a traditional distribution
    - Use percentiles to determine two-tailed shading
- Changed `method = "randomization"` to `method = "simulation"`
- Added warning when theoretical distribution is used that 
  assumptions should be checked  
- Added theoretical distributions to `visualize()` alone and as overlay with
current implementations being
    - Two sample t
    - ANOVA F
    - One proportion z
    - Two proportion z
    - Chi-square test of independence
    - Chi-square Goodness of Fit test
    - Standardized slope (t)
    
# infer 0.1.1
- Added additional tests
- Added `order` argument in `calculate()`
- Fixed bugs post-CRAN release
- Automated travis build of pkgdown to gh-pages branch

# infer 0.1.0
- Altered the way that successes are indicated in an infer pipeline. 
They now live in `specify()`.
- Updated documentation with examples
- Created `pkgdown` site materials
    - Deployed to https://infer.tidymodels.org/


# infer 0.0.1
- Implemented the "intro stats" examples for randomization methods
# Contributing

Contributions to the `infer` whether in the form of bug fixes, issue reports, new
code or documentation improvements are encouraged and welcome. We welcome novices
who may have never contributed to a package before as well as friendly
veterans looking to help us improve the package for users. We are eager to include
and accepting of contributions from everyone that meets our [code of conduct](CONDUCT.md)
guidelines.

Please use the GitHub issues. For any pull request, please link to or open a
corresponding issue in GitHub issues. Please ensure that you have notifications
turned on and respond to questions, comments or needed changes promptly.

##  Tests

`infer` uses `testthat` for testing. Please try to provide 100% test coverage
for any submitted code and always check that existing tests continue to pass.
If you are a beginner and need help with writing a test, mention this
in the issue and we will try to help.

It's also helpful to run `goodpractice::gp()` to ensure that lines of code are
not over 80 characters and that all lines of code have tests written. Please do
so prior to submitting any pull request and fix any suggestions from there.
Reach out to us if you need any assistance there too.

## Code style

Please use snake case (such as `rep_sample_n`) for function names.
Besides that, in general follow the 
[tidyverse style](http://style.tidyverse.org/) for R. 

## Code of Conduct

When contributing to the `infer` package you must follow the code of 
conduct defined in [CONDUCT](CONDUCT.md).
## Test environments
* local OS X install, R 4.1.0
* ubuntu 18.04 (on github actions), devel, release, oldrel, 3.6
* windows (on github actions), release
* OS X (on github actions), release
* windows (on win-builder), devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

We checked 6 reverse dependencies (5 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package, and saw no new problems.
---
title: 'infer: An R package for tidyverse-friendly statistical inference'
tags:
  - data science
  - tidyverse
  - inference
  - R
authors:
- name: Simon P. Couch
  orcid: 0000-0001-5676-5107
  affiliation: "1, 2"
- name: Andrew P. Bray
  orcid: 0000-0002-4037-7414
  affiliation: 3
- name: Chester Ismay
  orcid: 0000-0003-2820-2547
  affiliation: 4
- name: Evgeni Chasnovski
  orcid: 0000-0002-1617-4019
  affiliation: 5
- name: Benjamin S. Baumer
  orcid: 0000-0002-3279-0516
  affiliation: 6
- name: Mine Çetinkaya-Rundel
  orcid: 0000-0001-6452-2420
  affiliation: "2, 7"
affiliations:
 - name: Johns Hopkins, Department of Biostatistics
   index: 1
 - name: RStudio
   index: 2
 - name: UC Berkeley, Department of Statistics and Reed College Mathematics Department (on leave)
   index: 3
 - name: Flatiron School
   index: 4
 - name: No Affiliation
   index: 5
 - name: Smith College, Program in Statistical & Data Sciences
   index: 6
 - name: Duke University, Department of Statistical Science
   index: 7

citation_author: Couch et. al.
date: 12 June 2021
year: 2021
bibliography: paper.bib
output: 
  rticles::joss_article:
    keep_tex: true
    includes:
      in_header: columns.tex
csl: apa.csl
journal: JOSS
---

# Summary

`infer` implements an expressive grammar to perform statistical inference that adheres to the `tidyverse` design framework [@wickham2019welcome]. Rather than providing methods for specific statistical tests, this package consolidates the principles that are shared among common hypothesis tests and confidence intervals into a set of four main verbs (functions), supplemented with many utilities to visualize and extract value from their outputs.

# Statement of Need

Packages implementing methods for basic statistical inference in R are highly variable in their interfaces. The structure of inputted data, argument names, expected argument types, argument orders, output types, and spelling cases varies widely both within and among packages. This diversity in approaches obscures the intuition shared among common inferential procedures, makes details of usage difficult to remember, and prevents an expressive and idiomatic coding style.

`infer` is an R package for randomization-based hypothesis testing, naturalizing an intuitive understanding of statistical inference via a unified and expressive grammar. Four functions provide functionality encompassing a large swath of basic frequentist statistical inference, abstracting away details of specific tests and shifting the focus of the analyst to the observed data and the processes that generated it. Such a grammar lends itself to applications in teaching, data pedagogy research, applied scientific research, and advanced predictive modeling. For one, the principled approach of the `infer` package has made it an especially good fit for teaching introductory statistics and data science [@ismay2019statistical; @baumer2020teaching; @cetinkaya2021fresh] and research in data pedagogy [@fergusson2021introducing; @loy2021bringing]. Further, the package has already seen usage in a number of published scientific applications [@mclean2021controlled; @ask2021per; @fallon2021single]. Finally, the package integrates with the greater tidymodels collection of packages, a burgeoning software ecosystem for tidyverse-aligned predictive modeling used across many modern research and industrial applications [@kuhn2020tidymodels]. To date, the package has been downloaded more than 400,000 times.

# Underlying Principles

Regardless of the hypothesis test in question, an analyst asks the same kind of question when conducting statistical inference: is the effect/difference in the observed data real, or due to random chance? To answer this question, the analyst begins by assuming that the effect in the observed data was simply due to random chance, and calls this assumption the *null hypothesis*. (In reality, they might not believe in the null hypothesis at all---the null hypothesis is in opposition to the *alternate hypothesis*, which supposes that the effect present in the observed data is actually due to the fact that "something is going on.") The analyst then calculates a *test statistic* from the data that describes the observed effect. They can use this test statistic to calculate a *p-value* via juxtaposition with a *null distribution*, giving the probability that the observed data could come about if the null hypothesis were true. If this probability is below some pre-defined *significance level* $\alpha$, then the analyst can reject the null hypothesis.

The workflow of this package is designed around this idea. Starting out with some dataset,

+ `specify()` allows the analyst to specify the variable, or relationship between variables, that they are interested in.
+ `hypothesize()` allows the analyst to declare the null hypothesis.
+ `generate()` allows the analyst to generate data reflecting the null hypothesis or using the bootstrap.
+ `calculate()` allows the analyst to calculate summary statistics, either from
     * the observed data, to form the observed test statistic.
     * data `generate()`d to reflect the null hypothesis, to form a randomization-based null distribution of test statistics.

As such, the ultimate output of an infer pipeline using these four functions is generally an _observed statistic_ or _null distribution_ of test statistics. These four functions are thus supplemented with several utilities to visualize and extract value from their outputs.

+ `visualize()` plots the null distribution of test statistics.
     * `shade_p_value()` situates the observed statistic in the null distribution, shading the region as or more extreme.
+ `get_p_value()` calculates a p-value via the juxtaposition of the test statistic and the null distribution.

The workflow outlined above can also be used for constructing confidence intervals via bootstrapping with the omission of the `hypothesize()` step in the pipeline. The resulting bootstrap distribution can then be visualized with `visualize()`, the confidence interval region can be situated in the bootstrap distribution with `shade_confidence_interval()`, and the bounds of the confidence interval can be calculated with `get_confidence_interval()`.

Beyond this, the `infer` package offers:

* methods for inference using theory-based distributions
* shorthand wrappers for common statistical tests using tidy data
* model-fitting workflows to accommodate multiple explanatory variables

# Comparison to Other Packages

Several software packages on the Comprehensive R Archive Network share functionality with `infer` [@CRAN]. `broom` and `parameters` convert model objects to unified output formats, though they do not provide methods for fitting models, describing null distributions, performing bootstrapping, or calculating summary statistics from tabular data [@r-broom; @r-parameters]. `statsExpressions`, and adjacent packages in the `easystats` ecosystem, implement wrappers with consistent interfaces for theory-based hypothesis tests [@r-statsExpressions]. Similarly, `mosaic` is a package used to teach statistics by unifying summary statistics, visualization, and modeling with a consistent API built around R's formula interface. The `mosaic` package also includes functionality to conduct randomization-based inference [@r-mosaic]. At a higher level, though, the structure of each of these packages is defined by model types and statistics, where each model type or statistic has its own associated function and/or object class. In contrast, `infer` is structured around four functions, situating statistics and model types within a more abstracted grammar.^[This grammar follows from Allen Downey's "there is only one test" framework [@downey2016].] 

# Acknowledgements

We acknowledge contributions from Albert Y. Kim, Jo Hardin, Jay Lee, Amelia McNamara, Nick Solomon, and Richie Cotton.

# References
---
output: github_document
---

# infer R Package <img src="figs/infer_gnome.png" align="right" width=280 />




<!--figs/infer.svg-->
<!--http://www.r-pkg.org/badges/version/infer-->
<!--figs/main.svg-->
<!--https://img.shields.io/codecov/c/github/tidymodels/infer/main.svg-->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/infer)](https://cran.r-project.org/package=infer)[![Coverage Status](https://img.shields.io/codecov/c/github/tidymodels/infer/main.svg)](https://codecov.io/github/tidymodels/infer/?branch=main)

The objective of this package is to perform statistical inference using an expressive statistical grammar that coheres with the `tidyverse` design framework. The package is centered around 4 main verbs, supplemented with many utilities to visualize and extract value from their outputs.

+ `specify()` allows you to specify the variable, or relationship between variables, that you're interested in.
+ `hypothesize()` allows you to declare the null hypothesis.
+ `generate()` allows you to generate data reflecting the null hypothesis.
+ `calculate()` allows you to calculate a distribution of statistics from the generated data to form the null distribution.

To learn more about the principles underlying the package design, see `vignette("infer")`.

```{r load-package, echo = FALSE, message = FALSE, warning = FALSE}
library(devtools)
devtools::load_all()
```

```{r diagram, echo = FALSE, fig.cap = " ", fig.alt = "A diagram showing four steps to carry out randomization-based inference: specify hypothesis, generate data, calculate statistic, and visualize. From left to right, each step is connected by an arrow, while the diagram indicates that generating data and calculating statistics can happen iteratively."}
knitr::include_graphics("https://raw.githubusercontent.com/tidymodels/infer/main/figs/ht-diagram.png")
```

If you're interested in learning more about randomization-based statistical inference generally, including applied examples of this package, we recommend checking out [Statistical Inference Via Data Science: A ModernDive Into R and the Tidyverse](https://moderndive.com/) and [Introduction to Modern Statistics](https://openintro-ims.netlify.app/).

### Installation

------------------------------------------------------------------------

To install the current stable version of `infer` from CRAN:

```{r, eval = FALSE}
install.packages("infer")
```

To install the developmental stable version of `infer`, make sure to install `remotes` first. The `pkgdown` website for this version is at [infer.tidymodels.org](https://infer.tidymodels.org/).

```{r, eval = FALSE}
install.packages("remotes")
remotes::install_github("tidymodels/infer")
```

### Contributing

------------------------------------------------------------------------

We welcome others helping us make this package as user-friendly and efficient as possible. Please review our [contributing](https://github.com/tidymodels/infer/blob/main/CONTRIBUTING.md) and [conduct](https://github.com/tidymodels/infer/blob/main/CONDUCT.md) guidelines. By participating in this project you agree to abide by its terms.

For questions and discussions about tidymodels packages, modeling, and machine learning, please [post on RStudio Community](https://rstd.io/tidymodels-community). If you think you have encountered a bug, please [submit an issue](https://github.com/tidymodels/infer/issues). Either way, learn how to create and share a [reprex](https://rstd.io/reprex) (a minimal, reproducible example), to clearly communicate about your code. Check out further details on [contributing guidelines for tidymodels packages](https://www.tidymodels.org/contribute/) and [how to get help](https://www.tidymodels.org/help/).

### Examples

------------------------------------------------------------------------

These examples are pulled from the "Full infer Pipeline Examples" vignette, accessible by calling `vignette("observed_stat_examples")`. They make use of the `gss` dataset supplied by the package, providing a sample of data from the [General Social Survey](https://gss.norc.org). The data looks like this:

```{r load-gss, warning = FALSE, message = FALSE}
# load in the dataset
data(gss)

# take a glimpse at it
str(gss)
```

As an example, we'll run an analysis of variance on `age` and `partyid`, testing whether the age of a respondent is independent of their political party affiliation.

Calculating the observed statistic,

```{r, message = FALSE, warning = FALSE}
F_hat <- gss %>% 
  specify(age ~ partyid) %>%
  calculate(stat = "F")
```

Then, generating the null distribution,

```{r, message = FALSE, warning = FALSE}
null_dist <- gss %>%
   specify(age ~ partyid) %>%
   hypothesize(null = "independence") %>%
   generate(reps = 1000, type = "permute") %>%
   calculate(stat = "F")
```

Visualizing the observed statistic alongside the null distribution,

```{r viz, message = FALSE, warning = FALSE, eval = FALSE}
visualize(null_dist) +
  shade_p_value(obs_stat = F_hat, direction = "greater")
```

```{r viz-graphic, message = FALSE, warning = FALSE, echo = FALSE, fig.cap = " ", fig.alt = "A histogram showing a distribution of F statistics, right-tailed and centered around one. The x axis ranges from zero to five. The region of the histogram to the right of the observed statistic, just above two, is shaded red to represent the p-value."}
knitr::include_graphics("https://raw.githubusercontent.com/tidymodels/infer/main/README_files/figure-gfm/viz-1.png")
```

Calculating the p-value from the null distribution and observed statistic,

```{r, message = FALSE, warning = FALSE}
null_dist %>%
  get_p_value(obs_stat = F_hat, direction = "greater")
```


Note that the formula and non-formula interfaces  (i.e. `age ~ partyid` vs. `response = age, explanatory =  partyid`) work for all implemented inference procedures in `infer`. Use whatever is more natural for you. If you will be doing modeling using functions like `lm()` and `glm()`, though, we recommend you begin to use the formula `y ~ x` notation as soon as possible.

Other resources are available in the package vignettes! See `vignette("observed_stat_examples")` for more examples like the one above, and `vignette("infer")` for discussion of the underlying principles of the package design.
# Reproducibility

When using the infer package for research, or in other cases when exact reproducibility is a priority, be sure the set the seed for R's random number generator. infer will respect the random seed specified in the `set.seed()` function, returning the same result when `generate()`ing data given an identical seed. For instance, we can calculate the difference in mean `age` by `college` degree status using the `gss` dataset from 10 versions of the `gss` resampled with permutation using the following code.

```{r, include = FALSE}
library(infer)
```

```{r}
set.seed(1)

gss %>%
  specify(age ~ college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 5, type = "permute") %>%
  calculate("diff in means", order = c("degree", "no degree"))
```

Setting the seed to the same value again and rerunning the same code will produce the same result.

```{r}
# set the seed
set.seed(1)

gss %>%
  specify(age ~ college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 5, type = "permute") %>%
  calculate("diff in means", order = c("degree", "no degree"))
```

Please keep this in mind when writing infer code that utilizes resampling with `generate()`.
---
title: "Tidy ANOVA (Analysis of Variance) with infer"
description: "Conducting ANOVA (Analysis of Variance) on tidy data with infer."
output: rmarkdown::html_vignette
vignette: |
  %\VignetteIndexEntry{ANOVA}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r settings, include=FALSE}
knitr::opts_chunk$set(fig.width = 6, fig.height = 4.5) 
options(digits = 4)
```

```{r load-packages, echo = FALSE, message = FALSE, warning = FALSE}
library(ggplot2)
library(dplyr)
devtools::load_all()
```

In this vignette, we'll walk through conducting an analysis of variance (ANOVA) test using `infer`. ANOVAs are used to analyze differences in group means.

Throughout this vignette, we'll make use of the `gss` dataset supplied by `infer`, which contains a sample of data from the General Social Survey. See `?gss` for more information on the variables included and their source. Note that this data (and our examples on it) are for demonstration purposes only, and will not necessarily provide accurate estimates unless weighted properly. For these examples, let's suppose that this dataset is a representative sample of a population we want to learn about: American adults. The data looks like this:

```{r glimpse-gss-actual, warning = FALSE, message = FALSE}
dplyr::glimpse(gss)
```

To carry out an ANOVA, we'll examine the association between age and political party affiliation in the United States. The `age` variable is a numerical variable measuring the respondents' age at the time that the survey was taken, and `partyid` is a factor variable with unique values `r unique(gss$partyid)`.

This is what the relationship looks like in the observed data:

```{r plot-f, echo = FALSE}
gss %>%
  ggplot2::ggplot() +
  ggplot2::aes(x = partyid, y = age) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_brewer(type = "qual") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                     vjust = .5)) +
    ggplot2::labs(x = "partyid: Political Party Affiliation",
                  y = "age: Age of Respondent")
```

If there were no relationship, we would expect to see the each of these boxplots lining up along the y-axis. It looks like the average age of democrats and republicans seems to be a bit larger than independent and other American voters. Is this difference just random noise, though?

First, to calculate the observed statistic, we can use `specify()` and `calculate()`.

```{r calc-obs-stat-f, warning = FALSE, message = FALSE}
# calculate the observed statistic
observed_f_statistic <- gss %>%
  specify(age ~ partyid) %>%
  hypothesize(null = "independence") %>%
  calculate(stat = "F")
```

The observed $F$ statistic is `r observed_f_statistic`. Now, we want to compare this statistic to a null distribution, generated under the assumption that age and political party affiliation are not actually related, to get a sense of how likely it would be for us to see this observed statistic if there were actually no association between the two variables.

We can `generate` an approximation of the null distribution using randomization. The randomization approach permutes the response and explanatory variables, so that each person's party affiliation is matched up with a random age from the sample in order to break up any association between the two.

```{r generate-null-f, warning = FALSE, message = FALSE}
# generate the null distribution using randomization
null_dist <- gss %>%
  specify(age ~ partyid) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "F")
```

Note that, in the line `specify(age ~ partyid)` above, we could use the equivalent syntax `specify(response = age, explanatory = partyid)`. 

To get a sense for what this distribution looks like, and where our observed statistic falls, we can use `visualize()`:

```{r visualize-f, warning = FALSE, message = FALSE}
# visualize the null distribution and test statistic!
null_dist %>%
  visualize() + 
  shade_p_value(observed_f_statistic,
                direction = "greater")
```

We could also visualize the observed statistic against the theoretical null distribution. To do so, use the `assume()` verb to define a theoretical null distribution and then pass it to `visualize()` like a null distribution outputted from `generate()` and `calculate()`.

```{r visualize-f-theor, warning = FALSE, message = FALSE}
# visualize the theoretical null distribution and test statistic!
null_dist_theory <- gss %>%
  specify(age ~ partyid) %>%
  assume(distribution = "F")

visualize(null_dist_theory) +
  shade_p_value(observed_f_statistic,
                direction = "greater")
```

To visualize both the randomization-based and theoretical null distributions to get a sense of how the two relate, we can pipe the randomization-based null distribution into `visualize()`, and then further provide `method = "both"` to `visualize()`.

```{r visualize-indep-both, warning = FALSE, message = FALSE}
# visualize both null distributions and the test statistic!
null_dist %>%
  visualize(method = "both") + 
  shade_p_value(observed_f_statistic,
                direction = "greater")
```

Either way, it looks like our observed test statistic would be quite unlikely if there were actually no association between age and political party affiliation. More exactly, we can approximate the p-value from the randomization-based approximation to the null distribution:

```{r p-value-indep, warning = FALSE, message = FALSE}
# calculate the p value from the observed statistic and null distribution
p_value <- null_dist %>%
  get_p_value(obs_stat = observed_f_statistic,
              direction = "greater")

p_value
```

Thus, if there were really no relationship between age and political party affiliation, our approximation of the probability that we would see a statistic as or more extreme than `r observed_f_statistic` is approximately `r p_value`.

To calculate the p-value using the true $F$ distribution, we can use the `pf` function from base R. This function allows us to situate the test statistic we calculated previously in the $F$ distribution with the appropriate degrees of freedom.

```{r}
pf(observed_f_statistic$stat, 3, 496, lower.tail = FALSE)
```

Note that, while the observed statistic stays the same, the resulting p-value differs slightly between these two approaches since the randomization-based empirical $F$ distribution is an approximation of the true $F$ distribution.

The package currently does not supply a wrapper for tidy ANOVA tests.
---
title: "Full infer Pipeline Examples"
description: "A near-exhaustive demonstration of the functionality in infer."
output: 
  rmarkdown::html_vignette:
    df_print: kable
    toc: true
vignette: |
  %\VignetteIndexEntry{Full infer pipeline examples}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

#### Introduction

```{r include=FALSE}
knitr::opts_chunk$set(fig.width = 6, fig.height = 4.5, 
                      message = FALSE, warning = FALSE) 
options(digits = 4)
```

This vignette is intended to provide a set of examples that nearly exhaustively demonstrate the functionalities provided by `infer`. Commentary on these examples is limited---for more discussion of the intuition behind the package, see the "Getting to Know infer" vignette, accessible by calling `vignette("infer")`.

Throughout this vignette, we'll make use of the `gss` dataset supplied by `infer`, which contains a sample of data from the General Social Survey. See `?gss` for more information on the variables included and their source. Note that this data (and our examples on it) are for demonstration purposes only, and will not necessarily provide accurate estimates unless weighted properly. For these examples, let's suppose that this dataset is a representative sample of a population we want to learn about: American adults. The data looks like this:

```{r load-packages, echo = FALSE}
library(dplyr)
devtools::load_all()
```

```{r load-gss}
# load in the dataset
data(gss)

# take a look at its structure
dplyr::glimpse(gss)
```

## Hypothesis tests

### One numerical variable (mean)

Calculating the observed statistic,

```{r}
x_bar <- gss %>%
  specify(response = hours) %>%
  calculate(stat = "mean")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
x_bar <- gss %>%
  observe(response = hours, stat = "mean")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  generate(reps = 1000) %>%
  calculate(stat = "mean")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = x_bar, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = x_bar, direction = "two-sided")
```

### One numerical variable (standardized mean $t$)

Calculating the observed statistic,

```{r}
t_bar <- gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  calculate(stat = "t")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
t_bar <- gss %>%
  observe(response = hours, null = "point", mu = 40, stat = "t")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  generate(reps = 1000) %>%
  calculate(stat = "t")
```

Alternatively, finding the null distribution using theoretical methods using the `assume()` verb,

```{r}
null_dist_theory <- gss %>%
  specify(response = hours)  %>%
  assume("t")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = t_bar, direction = "two-sided")
```

Alternatively, visualizing the observed statistic using the theory-based null distribution,

```{r}
visualize(null_dist_theory) +
  shade_p_value(obs_stat = t_bar, direction = "two-sided")
```

Alternatively, visualizing the observed statistic using both of the null distributions,

```{r}
visualize(null_dist, method = "both") +
  shade_p_value(obs_stat = t_bar, direction = "two-sided")
```

Note that the above code makes use of the randomization-based null distribution.

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = t_bar, direction = "two-sided")
```

Alternatively, using the `t_test` wrapper:

```{r}
gss %>%
  t_test(response = hours, mu = 40)
```

`infer` does not support testing on one numerical variable via the `z` distribution.

### One numerical variable (median)

Calculating the observed statistic,

```{r}
x_tilde <- gss %>%
  specify(response = age) %>%
  calculate(stat = "median")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
x_tilde <- gss %>%
  observe(response = age, stat = "median")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(response = age) %>%
  hypothesize(null = "point", med = 40) %>% 
  generate(reps = 1000) %>% 
  calculate(stat = "median")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = x_tilde, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = x_tilde, direction = "two-sided")
```

### One categorical (one proportion)

Calculating the observed statistic,

```{r}
p_hat <- gss %>%
  specify(response = sex, success = "female") %>%
  calculate(stat = "prop")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
p_hat <- gss %>%
  observe(response = sex, success = "female", stat = "prop")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(response = sex, success = "female") %>%
  hypothesize(null = "point", p = .5) %>%
  generate(reps = 1000) %>%
  calculate(stat = "prop")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = p_hat, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = p_hat, direction = "two-sided")
```

Note that logical variables will be coerced to factors:

```{r}
null_dist <- gss %>%
  dplyr::mutate(is_female = (sex == "female")) %>%
  specify(response = is_female, success = "TRUE") %>%
  hypothesize(null = "point", p = .5) %>%
  generate(reps = 1000) %>%
  calculate(stat = "prop")
```

### One categorical variable (standardized proportion $z$)

Calculating the observed statistic,

```{r}
p_hat <- gss %>%
  specify(response = sex, success = "female") %>%
  hypothesize(null = "point", p = .5) %>%
  calculate(stat = "z")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
p_hat <- gss %>%
  observe(response = sex, success = "female", null = "point", p = .5, stat = "z")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(response = sex, success = "female") %>%
  hypothesize(null = "point", p = .5) %>%
  generate(reps = 1000, type = "draw") %>%
  calculate(stat = "z")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = p_hat, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = p_hat, direction = "two-sided")
```

The package also supplies a wrapper around `prop.test` for tests of a single proportion on tidy data.

```{r prop_test_1_grp}
prop_test(gss,
          college ~ NULL,
          p = .2)
```

`infer` does not support testing two means via the `z` distribution.

### Two categorical (2 level) variables

The `infer` package provides several statistics to work with data of this type. One of them is the statistic for difference in proportions.

Calculating the observed statistic,

```{r}
d_hat <- gss %>% 
  specify(college ~ sex, success = "no degree") %>%
  calculate(stat = "diff in props", order = c("female", "male"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
d_hat <- gss %>% 
  observe(college ~ sex, success = "no degree", 
          stat = "diff in props", order = c("female", "male"))
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(college ~ sex, success = "no degree") %>%
  hypothesize(null = "independence") %>% 
  generate(reps = 1000) %>% 
  calculate(stat = "diff in props", order = c("female", "male"))
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = d_hat, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = d_hat, direction = "two-sided")
```

`infer` also provides functionality to calculate ratios of proportions. The workflow looks similar to that for `diff in props`.

Calculating the observed statistic,

```{r}
r_hat <- gss %>% 
  specify(college ~ sex, success = "no degree") %>%
  calculate(stat = "ratio of props", order = c("female", "male"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
r_hat <- gss %>% 
  observe(college ~ sex, success = "no degree",
          stat = "ratio of props", order = c("female", "male"))
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(college ~ sex, success = "no degree") %>%
  hypothesize(null = "independence") %>% 
  generate(reps = 1000) %>% 
  calculate(stat = "ratio of props", order = c("female", "male"))
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = r_hat, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = r_hat, direction = "two-sided")
```

In addition, the package provides functionality to calculate odds ratios. The workflow also looks similar to that for `diff in props`.

Calculating the observed statistic,

```{r}
or_hat <- gss %>% 
  specify(college ~ sex, success = "no degree") %>%
  calculate(stat = "odds ratio", order = c("female", "male"))
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(college ~ sex, success = "no degree") %>%
  hypothesize(null = "independence") %>% 
  generate(reps = 1000) %>% 
  calculate(stat = "odds ratio", order = c("female", "male"))
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = or_hat, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = or_hat, direction = "two-sided")
```

### Two categorical (2 level) variables (z)

Finding the standardized observed statistic,

```{r}
z_hat <- gss %>% 
  specify(college ~ sex, success = "no degree") %>%
  hypothesize(null = "independence") %>%
  calculate(stat = "z", order = c("female", "male"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
z_hat <- gss %>% 
  observe(college ~ sex, success = "no degree",
          stat = "z", order = c("female", "male"))
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(college ~ sex, success = "no degree") %>%
  hypothesize(null = "independence") %>% 
  generate(reps = 1000) %>% 
  calculate(stat = "z", order = c("female", "male"))
```

Alternatively, finding the null distribution using theoretical methods using the `assume()` verb,

```{r}
null_dist_theory <- gss %>%
  specify(college ~ sex, success = "no degree") %>%
  assume("z")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = z_hat, direction = "two-sided")
```

Alternatively, visualizing the observed statistic using the theory-based null distribution,

```{r}
visualize(null_dist_theory) +
  shade_p_value(obs_stat = z_hat, direction = "two-sided")
```

Alternatively, visualizing the observed statistic using both of the null distributions,

```{r}
visualize(null_dist, method = "both") +
  shade_p_value(obs_stat = z_hat, direction = "two-sided")
```

Note that the above code makes use of the randomization-based null distribution.

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = z_hat, direction = "two-sided")
```

Note the similarities in this plot and the previous one.

The package also supplies a wrapper around `prop.test` to allow for tests of equality of proportions on tidy data.

```{r prop_test_2_grp}
prop_test(gss, 
          college ~ sex,  
          order = c("female", "male"))
```

### One categorical (\>2 level) - GoF

Calculating the observed statistic,

Note the need to add in the hypothesized values here to compute the observed statistic.

```{r}
Chisq_hat <- gss %>%
  specify(response = finrela) %>%
  hypothesize(null = "point",
              p = c("far below average" = 1/6,
                    "below average" = 1/6,
                    "average" = 1/6,
                    "above average" = 1/6,
                    "far above average" = 1/6,
                    "DK" = 1/6)) %>%
  calculate(stat = "Chisq")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
Chisq_hat <- gss %>%
  observe(response = finrela,
          null = "point",
          p = c("far below average" = 1/6,
                "below average" = 1/6,
                "average" = 1/6,
                "above average" = 1/6,
                "far above average" = 1/6,
                "DK" = 1/6),
          stat = "Chisq")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(response = finrela) %>%
  hypothesize(null = "point",
              p = c("far below average" = 1/6,
                    "below average" = 1/6,
                    "average" = 1/6,
                    "above average" = 1/6,
                    "far above average" = 1/6,
                    "DK" = 1/6)) %>%
  generate(reps = 1000, type = "draw") %>%
  calculate(stat = "Chisq")
```

Alternatively, finding the null distribution using theoretical methods using the `assume()` verb,

```{r}
null_dist_theory <- gss %>%
  specify(response = finrela) %>%
  assume("Chisq")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = Chisq_hat, direction = "greater")
```

Alternatively, visualizing the observed statistic using the theory-based null distribution,

```{r}
visualize(null_dist_theory) +
  shade_p_value(obs_stat = Chisq_hat, direction = "greater")
```

Alternatively, visualizing the observed statistic using both of the null distributions,

```{r}
visualize(null_dist_theory, method = "both") +
  shade_p_value(obs_stat = Chisq_hat, direction = "greater")
```

Note that the above code makes use of the randomization-based null distribution.

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = Chisq_hat, direction = "greater")
```

Alternatively, using the `chisq_test` wrapper:

```{r}
chisq_test(gss, 
           response = finrela,
           p = c("far below average" = 1/6,
                 "below average" = 1/6,
                 "average" = 1/6,
                 "above average" = 1/6,
                 "far above average" = 1/6,
                 "DK" = 1/6))
```

### Two categorical (\>2 level): Chi-squared test of independence

Calculating the observed statistic,

```{r}
Chisq_hat <- gss %>%
  specify(formula = finrela ~ sex) %>% 
  hypothesize(null = "independence") %>%
  calculate(stat = "Chisq")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
Chisq_hat <- gss %>%
  observe(formula = finrela ~ sex, stat = "Chisq")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(finrela ~ sex) %>%
  hypothesize(null = "independence") %>% 
  generate(reps = 1000, type = "permute") %>% 
  calculate(stat = "Chisq")
```

Alternatively, finding the null distribution using theoretical methods using the `assume()` verb,

```{r}
null_dist_theory <- gss %>%
  specify(finrela ~ sex) %>%
  assume(distribution = "Chisq")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = Chisq_hat, direction = "greater")
```

Alternatively, visualizing the observed statistic using the theory-based null distribution,

```{r}
visualize(null_dist_theory) +
  shade_p_value(obs_stat = Chisq_hat, direction = "greater")
```

Alternatively, visualizing the observed statistic using both of the null distributions,

```{r}
visualize(null_dist, method = "both") +
  shade_p_value(obs_stat = Chisq_hat, direction = "greater")
```

Note that the above code makes use of the randomization-based null distribution.

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = Chisq_hat, direction = "greater")
```

Alternatively, using the wrapper to carry out the test,

```{r}
gss %>%
  chisq_test(formula = finrela ~ sex)
```

### One numerical variable, one categorical (2 levels) (diff in means)

Calculating the observed statistic,

```{r}
d_hat <- gss %>% 
  specify(age ~ college) %>% 
  calculate(stat = "diff in means", order = c("degree", "no degree"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
d_hat <- gss %>% 
  observe(age ~ college,
          stat = "diff in means", order = c("degree", "no degree"))
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(age ~ college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "diff in means", order = c("degree", "no degree"))
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = d_hat, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = d_hat, direction = "two-sided")
```

### One numerical variable, one categorical (2 levels) (t)

Finding the standardized observed statistic,

```{r}
t_hat <- gss %>% 
  specify(age ~ college) %>% 
  hypothesize(null = "independence") %>%
  calculate(stat = "t", order = c("degree", "no degree"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
t_hat <- gss %>% 
  observe(age ~ college,
          stat = "t", order = c("degree", "no degree"))
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(age ~ college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "t", order = c("degree", "no degree"))
```

Alternatively, finding the null distribution using theoretical methods using the `assume()` verb,

```{r}
null_dist_theory <- gss %>%
  specify(age ~ college) %>%
  assume("t")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = t_hat, direction = "two-sided")
```

Alternatively, visualizing the observed statistic using the theory-based null distribution,

```{r}
visualize(null_dist_theory) +
  shade_p_value(obs_stat = t_hat, direction = "two-sided")
```

Alternatively, visualizing the observed statistic using both of the null distributions,

```{r}
visualize(null_dist, method = "both") +
  shade_p_value(obs_stat = t_hat, direction = "two-sided")
```

Note that the above code makes use of the randomization-based null distribution.

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = t_hat, direction = "two-sided")
```

Note the similarities in this plot and the previous one.

### One numerical variable, one categorical (2 levels) (diff in medians)

Calculating the observed statistic,

```{r}
d_hat <- gss %>% 
  specify(age ~ college) %>% 
  calculate(stat = "diff in medians", order = c("degree", "no degree"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
d_hat <- gss %>% 
  observe(age ~ college,
          stat = "diff in medians", order = c("degree", "no degree"))
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
  specify(age ~ college) %>% # alt: response = age, explanatory = season
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "diff in medians", order = c("degree", "no degree"))
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = d_hat, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = d_hat, direction = "two-sided")
```

### One numerical, one categorical (\>2 levels) - ANOVA

Calculating the observed statistic,

```{r}
F_hat <- gss %>% 
  specify(age ~ partyid) %>%
  calculate(stat = "F")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
F_hat <- gss %>% 
  observe(age ~ partyid, stat = "F")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
   specify(age ~ partyid) %>%
   hypothesize(null = "independence") %>%
   generate(reps = 1000, type = "permute") %>%
   calculate(stat = "F")
```

Alternatively, finding the null distribution using theoretical methods using the `assume()` verb,

```{r}
null_dist_theory <- gss %>%
   specify(age ~ partyid) %>%
   hypothesize(null = "independence") %>%
   assume(distribution = "F")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = F_hat, direction = "greater")
```

Alternatively, visualizing the observed statistic using the theory-based null distribution,

```{r}
visualize(null_dist_theory) +
  shade_p_value(obs_stat = F_hat, direction = "greater")
```

Alternatively, visualizing the observed statistic using both of the null distributions,

```{r}
visualize(null_dist, method = "both") +
  shade_p_value(obs_stat = F_hat, direction = "greater")
```

Note that the above code makes use of the randomization-based null distribution.

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = F_hat, direction = "greater")
```

### Two numerical vars - SLR

Calculating the observed statistic,

```{r}
slope_hat <- gss %>% 
  specify(hours ~ age) %>% 
  calculate(stat = "slope")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
slope_hat <- gss %>% 
  observe(hours ~ age, stat = "slope")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
   specify(hours ~ age) %>% 
   hypothesize(null = "independence") %>%
   generate(reps = 1000, type = "permute") %>%
   calculate(stat = "slope")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = slope_hat, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = slope_hat, direction = "two-sided")
```

### Two numerical vars - correlation

Calculating the observed statistic,

```{r}
correlation_hat <- gss %>% 
  specify(hours ~ age) %>% 
  calculate(stat = "correlation")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
correlation_hat <- gss %>% 
  observe(hours ~ age, stat = "correlation")
```

Then, generating the null distribution,

```{r}
null_dist <- gss %>%
   specify(hours ~ age) %>% 
   hypothesize(null = "independence") %>%
   generate(reps = 1000, type = "permute") %>%
   calculate(stat = "correlation")
```

Visualizing the observed statistic alongside the null distribution,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = correlation_hat, direction = "two-sided")
```

Calculating the p-value from the null distribution and observed statistic,

```{r}
null_dist %>%
  get_p_value(obs_stat = correlation_hat, direction = "two-sided")
```

### Two numerical vars - SLR (t)

Not currently implemented since $t$ could refer to standardized slope or standardized correlation.

### Multiple explanatory variables

Calculating the observed fit,

```{r}
obs_fit <- gss %>%
  specify(hours ~ age + college) %>%
  fit()
```

Generating a distribution of fits with the response variable permuted,

```{r}
null_dist <- gss %>%
  specify(hours ~ age + college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  fit()
```

Generating a distribution of fits where each explanatory variable is permuted independently,

```{r}
null_dist2 <- gss %>%
  specify(hours ~ age + college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute", variables = c(age, college)) %>%
  fit()
```

Visualizing the observed fit alongside the null fits,

```{r}
visualize(null_dist) +
  shade_p_value(obs_stat = obs_fit, direction = "two-sided")
```

Calculating p-values from the null distribution and observed fit,

```{r}
null_dist %>%
  get_p_value(obs_stat = obs_fit, direction = "two-sided")
```

Note that this `fit()`-based workflow can be applied to use cases with differing numbers of explanatory variables and explanatory variable types.

## Confidence intervals

### One numerical (one mean)

Finding the observed statistic,

```{r}
x_bar <- gss %>% 
  specify(response = hours) %>%
  calculate(stat = "mean")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
x_bar <- gss %>% 
  observe(response = hours, stat = "mean")
```

Then, generating a bootstrap distribution,

```{r}
boot_dist <- gss %>%
   specify(response = hours) %>%
   generate(reps = 1000, type = "bootstrap") %>%
   calculate(stat = "mean")
```

Use the bootstrap distribution to find a confidence interval,

```{r}
percentile_ci <- get_ci(boot_dist)
```

Visualizing the observed statistic alongside the distribution,

```{r}
visualize(boot_dist) +
  shade_confidence_interval(endpoints = percentile_ci)
```

Alternatively, use the bootstrap distribution to find a confidence interval using the standard error,

```{r}
standard_error_ci <- get_ci(boot_dist, type = "se", point_estimate = x_bar)

visualize(boot_dist) +
  shade_confidence_interval(endpoints = standard_error_ci)
```

Instead of a simulation-based bootstrap distribution, we can also define a theory-based sampling distribution,

```{r}
sampling_dist <- gss %>%
   specify(response = hours) %>%
   assume(distribution = "t")
```

Visualization and calculation of confidence intervals interfaces in the same way as with the simulation-based distribution,

```{r}
theor_ci <- get_ci(sampling_dist, point_estimate = x_bar)

theor_ci

visualize(sampling_dist) +
  shade_confidence_interval(endpoints = theor_ci)
```

Note that the `t` distribution is recentered and rescaled to lie on the scale of the observed data. `infer` does not support confidence intervals on means via the `z` distribution.

### One numerical (one mean - standardized)

Finding the observed statistic,

```{r}
t_hat <- gss %>% 
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  calculate(stat = "t")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
t_hat <- gss %>% 
  observe(response = hours,
          null = "point", mu = 40,
          stat = "t")
```

Then, generating the bootstrap distribution,

```{r}
boot_dist <- gss %>%
   specify(response = hours) %>%
   generate(reps = 1000, type = "bootstrap") %>%
   calculate(stat = "t")
```

Use the bootstrap distribution to find a confidence interval,

```{r}
percentile_ci <- get_ci(boot_dist)
```

Visualizing the observed statistic alongside the distribution,

```{r}
visualize(boot_dist) +
  shade_confidence_interval(endpoints = percentile_ci)
```

Alternatively, use the bootstrap distribution to find a confidence interval using the standard error,

```{r}
standard_error_ci <- boot_dist %>%
  get_ci(type = "se", point_estimate = t_hat)

visualize(boot_dist) +
  shade_confidence_interval(endpoints = standard_error_ci)
```

See the above subsection (one mean) for a theory-based approach. Note that `infer` does not support confidence intervals on means via the `z` distribution.

### One categorical (one proportion)

Finding the observed statistic,

```{r}
p_hat <- gss %>% 
   specify(response = sex, success = "female") %>%
   calculate(stat = "prop")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
p_hat <- gss %>% 
   observe(response = sex, success = "female", stat = "prop")
```

Then, generating a bootstrap distribution,

```{r}
boot_dist <- gss %>%
 specify(response = sex, success = "female") %>%
 generate(reps = 1000, type = "bootstrap") %>%
 calculate(stat = "prop")
```

Use the bootstrap distribution to find a confidence interval,

```{r}
percentile_ci <- get_ci(boot_dist)
```

Visualizing the observed statistic alongside the distribution,

```{r}
visualize(boot_dist) +
  shade_confidence_interval(endpoints = percentile_ci)
```

Alternatively, use the bootstrap distribution to find a confidence interval using the standard error,

```{r}
standard_error_ci <- boot_dist %>%
  get_ci(type = "se", point_estimate = p_hat)

visualize(boot_dist) +
  shade_confidence_interval(endpoints = standard_error_ci)
```

Instead of a simulation-based bootstrap distribution, we can also define a theory-based sampling distribution,

```{r}
sampling_dist <- gss %>%
   specify(response = sex, success = "female") %>%
   assume(distribution = "z")
```

Visualization and calculation of confidence intervals interfaces in the same way as with the simulation-based distribution,

```{r}
theor_ci <- get_ci(sampling_dist, point_estimate = p_hat)

theor_ci

visualize(sampling_dist) +
  shade_confidence_interval(endpoints = theor_ci)
```

Note that the `z` distribution is recentered and rescaled to lie on the scale of the observed data. `infer` does not support confidence intervals on means via the `z` distribution.

### One categorical variable (standardized proportion $z$)

See the above subsection (one proportion) for a theory-based approach.

### One numerical variable, one categorical (2 levels) (diff in means)

Finding the observed statistic,

```{r}
d_hat <- gss %>%
  specify(hours ~ college) %>%
  calculate(stat = "diff in means", order = c("degree", "no degree"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
d_hat <- gss %>%
  observe(hours ~ college,
          stat = "diff in means", order = c("degree", "no degree"))
```

Then, generating a bootstrap distribution,

```{r}
boot_dist <- gss %>%
   specify(hours ~ college) %>%
   generate(reps = 1000, type = "bootstrap") %>%
   calculate(stat = "diff in means", order = c("degree", "no degree"))
```

Use the bootstrap distribution to find a confidence interval,

```{r}
percentile_ci <- get_ci(boot_dist)
```

Visualizing the observed statistic alongside the distribution,

```{r}
visualize(boot_dist) +
  shade_confidence_interval(endpoints = percentile_ci)
```

Alternatively, use the bootstrap distribution to find a confidence interval using the standard error,

```{r}
standard_error_ci <- boot_dist %>%
  get_ci(type = "se", point_estimate = d_hat)

visualize(boot_dist) +
  shade_confidence_interval(endpoints = standard_error_ci)
```

Instead of a simulation-based bootstrap distribution, we can also define a theory-based sampling distribution,

```{r}
sampling_dist <- gss %>%
   specify(hours ~ college) %>%
   assume(distribution = "t")
```

Visualization and calculation of confidence intervals interfaces in the same way as with the simulation-based distribution,

```{r}
theor_ci <- get_ci(sampling_dist, point_estimate = d_hat)

theor_ci

visualize(sampling_dist) +
  shade_confidence_interval(endpoints = theor_ci)
```

Note that the `t` distribution is recentered and rescaled to lie on the scale of the observed data.

### One numerical variable, one categorical (2 levels) (t)

Finding the standardized point estimate,

```{r}
t_hat <- gss %>%
  specify(hours ~ college) %>%
  calculate(stat = "t", order = c("degree", "no degree"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
t_hat <- gss %>%
  observe(hours ~ college,
          stat = "t", order = c("degree", "no degree"))
```

Then, generating a bootstrap distribution,

```{r}
boot_dist <- gss %>%
   specify(hours ~ college) %>%
   generate(reps = 1000, type = "bootstrap") %>%
   calculate(stat = "t", order = c("degree", "no degree"))
```

Use the bootstrap distribution to find a confidence interval,

```{r}
percentile_ci <- get_ci(boot_dist)
```

Visualizing the observed statistic alongside the distribution,

```{r}
visualize(boot_dist) +
  shade_confidence_interval(endpoints = percentile_ci)
```

Alternatively, use the bootstrap distribution to find a confidence interval using the standard error,

```{r}
standard_error_ci <- boot_dist %>%
  get_ci(type = "se", point_estimate = t_hat)

visualize(boot_dist) +
  shade_confidence_interval(endpoints = standard_error_ci)
```

See the above subsection (diff in means) for a theory-based approach. `infer` does not support confidence intervals on means via the `z` distribution.

### Two categorical variables (diff in proportions)

Finding the observed statistic,

```{r}
d_hat <- gss %>% 
  specify(college ~ sex, success = "degree") %>%
  calculate(stat = "diff in props", order = c("female", "male"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
d_hat <- gss %>% 
  observe(college ~ sex, success = "degree",
          stat = "diff in props", order = c("female", "male"))
```

Then, generating a bootstrap distribution,

```{r}
boot_dist <- gss %>%
  specify(college ~ sex, success = "degree") %>%
  generate(reps = 1000, type = "bootstrap") %>% 
  calculate(stat = "diff in props", order = c("female", "male"))
```

Use the bootstrap distribution to find a confidence interval,

```{r}
percentile_ci <- get_ci(boot_dist)
```

Visualizing the observed statistic alongside the distribution,

```{r}
visualize(boot_dist) +
  shade_confidence_interval(endpoints = percentile_ci)
```

Alternatively, use the bootstrap distribution to find a confidence interval using the standard error,

```{r}
standard_error_ci <- boot_dist %>%
  get_ci(type = "se", point_estimate = d_hat)

visualize(boot_dist) +
  shade_confidence_interval(endpoints = standard_error_ci)
```

Instead of a simulation-based bootstrap distribution, we can also define a theory-based sampling distribution,

```{r}
sampling_dist <- gss %>% 
  specify(college ~ sex, success = "degree") %>%
   assume(distribution = "z")
```

Visualization and calculation of confidence intervals interfaces in the same way as with the simulation-based distribution,

```{r}
theor_ci <- get_ci(sampling_dist, point_estimate = d_hat)

theor_ci

visualize(sampling_dist) +
  shade_confidence_interval(endpoints = theor_ci)
```

Note that the `z` distribution is recentered and rescaled to lie on the scale of the observed data.

### Two categorical variables (z)

Finding the standardized point estimate,

```{r}
z_hat <- gss %>% 
  specify(college ~ sex, success = "degree") %>%
  calculate(stat = "z", order = c("female", "male"))
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
z_hat <- gss %>% 
  observe(college ~ sex, success = "degree",
          stat = "z", order = c("female", "male"))
```

Then, generating a bootstrap distribution,

```{r}
boot_dist <- gss %>%
  specify(college ~ sex, success = "degree") %>%
  generate(reps = 1000, type = "bootstrap") %>% 
  calculate(stat = "z", order = c("female", "male"))
```

Use the bootstrap distribution to find a confidence interval,

```{r}
percentile_ci <- get_ci(boot_dist)
```

Visualizing the observed statistic alongside the distribution,

```{r}
visualize(boot_dist) +
  shade_confidence_interval(endpoints = percentile_ci)
```

Alternatively, use the bootstrap distribution to find a confidence interval using the standard error,

```{r}
standard_error_ci <- boot_dist %>%
  get_ci(type = "se", point_estimate = z_hat)

visualize(boot_dist) +
  shade_confidence_interval(endpoints = standard_error_ci)
```

See the above subsection (diff in props) for a theory-based approach.

### Two numerical vars - SLR

Finding the observed statistic,

```{r}
slope_hat <- gss %>% 
  specify(hours ~ age) %>%
  calculate(stat = "slope")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
slope_hat <- gss %>% 
  observe(hours ~ age, stat = "slope")
```

Then, generating a bootstrap distribution,

```{r}
boot_dist <- gss %>%
   specify(hours ~ age) %>% 
   generate(reps = 1000, type = "bootstrap") %>%
   calculate(stat = "slope")
```

Use the bootstrap distribution to find a confidence interval,

```{r}
percentile_ci <- get_ci(boot_dist)
```

Visualizing the observed statistic alongside the distribution,

```{r}
visualize(boot_dist) +
  shade_confidence_interval(endpoints = percentile_ci)
```

Alternatively, use the bootstrap distribution to find a confidence interval using the standard error,

```{r}
standard_error_ci <- boot_dist %>%
  get_ci(type = "se", point_estimate = slope_hat)

visualize(boot_dist) +
  shade_confidence_interval(endpoints = standard_error_ci)
```

### Two numerical vars - correlation

Finding the observed statistic,

```{r}
correlation_hat <- gss %>% 
  specify(hours ~ age) %>%
  calculate(stat = "correlation")
```

Alternatively, using the `observe()` wrapper to calculate the observed statistic,

```{r}
correlation_hat <- gss %>% 
  observe(hours ~ age, stat = "correlation")
```

Then, generating a bootstrap distribution,

```{r}
boot_dist <- gss %>%
   specify(hours ~ age) %>% 
   generate(reps = 1000, type = "bootstrap") %>%
   calculate(stat = "correlation")
```

Use the bootstrap distribution to find a confidence interval,

```{r}
percentile_ci <- get_ci(boot_dist)
```

Visualizing the observed statistic alongside the distribution,

```{r}
visualize(boot_dist) +
  shade_confidence_interval(endpoints = percentile_ci)
```

Alternatively, use the bootstrap distribution to find a confidence interval using the standard error,

```{r}
standard_error_ci <- boot_dist %>%
  get_ci(type = "se", point_estimate = correlation_hat)

visualize(boot_dist) +
  shade_confidence_interval(endpoints = standard_error_ci)
```

### Two numerical vars - t

Not currently implemented since $t$ could refer to standardized slope or standardized correlation.

### Multiple explanatory variables

Calculating the observed fit,

```{r}
obs_fit <- gss %>%
  specify(hours ~ age + college) %>%
  fit()
```

Generating a distribution of fits with the response variable permuted,

```{r}
null_dist <- gss %>%
  specify(hours ~ age + college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  fit()
```

Alternatively, generating a distribution of fits where each explanatory variable is permuted independently,

```{r}
null_dist2 <- gss %>%
  specify(hours ~ age + college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute", variables = c(age, college)) %>%
  fit()
```

Calculating confidence intervals from the null fits,

```{r}
conf_ints <- 
  get_confidence_interval(
    null_dist, 
    level = .95, 
    point_estimate = obs_fit
  )
```

Visualizing the observed fit alongside the null fits,

```{r}
visualize(null_dist) +
  shade_confidence_interval(endpoints = conf_ints)
```

Note that this `fit()`-based workflow can be applied to use cases with differing numbers of explanatory variables and explanatory variable types.
---
title: "Tidy Chi-Squared Tests with infer"
description: "Conducting Chi-Squared tests on tidy data with infer."
output: rmarkdown::html_vignette
vignette: |
  %\VignetteIndexEntry{Chi-Squared Tests}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r settings, include=FALSE}
knitr::opts_chunk$set(fig.width = 6, fig.height = 4.5) 
options(digits = 4)
```

```{r load-packages, echo = FALSE, message = FALSE, warning = FALSE}
library(ggplot2)
library(dplyr)
devtools::load_all()
```

### Introduction

In this vignette, we'll walk through conducting a $\chi^2$ (chi-squared) test of independence and a chi-squared goodness of fit test using `infer`. We'll start out with a  chi-squared test of independence, which can be used to test the association between two categorical variables. Then, we'll move on to a chi-squared goodness of fit test, which tests how well the distribution of one categorical variable can be approximated by some theoretical distribution.

Throughout this vignette, we'll make use of the `gss` dataset supplied by `infer`, which contains a sample of data from the General Social Survey. See `?gss` for more information on the variables included and their source. Note that this data (and our examples on it) are for demonstration purposes only, and will not necessarily provide accurate estimates unless weighted properly. For these examples, let's suppose that this dataset is a representative sample of a population we want to learn about: American adults. The data looks like this:

```{r glimpse-gss-actual, warning = FALSE, message = FALSE}
dplyr::glimpse(gss)
```

### Test of Independence

To carry out a chi-squared test of independence, we'll examine the association between income and educational attainment in the United States. `college` is a categorical variable with values `degree` and `no degree`, indicating whether or not the respondent has a college degree (including community college), and `finrela` gives the respondent's self-identification of family income---either `far below average`, `below average`, `average`, `above average`, `far above average`, or `DK` (don't know).

This is what the relationship looks like in the sample data:

```{r plot-indep, echo = FALSE}
gss %>%
  ggplot2::ggplot() +
  ggplot2::aes(x = finrela, fill = college) +
  ggplot2::geom_bar(position = "fill") +
  ggplot2::scale_fill_brewer(type = "qual") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                     vjust = .5)) +
    ggplot2::labs(x = "finrela: Self-Identification of Income Class",
                  y = "Proportion")
```

If there were no relationship, we would expect to see the purple bars reaching to the same height, regardless of income class. Are the differences we see here, though, just due to random noise?

First, to calculate the observed statistic, we can use `specify()` and `calculate()`.

```{r calc-obs-stat-indep, warning = FALSE, message = FALSE}
# calculate the observed statistic
observed_indep_statistic <- gss %>%
  specify(college ~ finrela) %>%
  hypothesize(null = "independence") %>%
  calculate(stat = "Chisq")
```

The observed $\chi^2$ statistic is `r observed_indep_statistic`. Now, we want to compare this statistic to a null distribution, generated under the assumption that these variables are not actually related, to get a sense of how likely it would be for us to see this observed statistic if there were actually no association between education and income.

We can `generate` the null distribution in one of two ways---using randomization or theory-based methods. The randomization approach approximates the null distribution by permuting the response and explanatory variables, so that each person's educational attainment is matched up with a random income from the sample in order to break up any association between the two.

```{r generate-null-indep, warning = FALSE, message = FALSE}
# generate the null distribution using randomization
null_dist_sim <- gss %>%
  specify(college ~ finrela) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "Chisq")
```

Note that, in the line `specify(college ~ finrela)` above, we could use the equivalent syntax `specify(response = college, explanatory = finrela)`. The same goes in the code below, which generates the null distribution using theory-based methods instead of randomization.

```{r generate-null-indep-t, warning = FALSE, message = FALSE}
# generate the null distribution by theoretical approximation
null_dist_theory <- gss %>%
  specify(college ~ finrela) %>%
  assume(distribution = "Chisq")
```

To get a sense for what these distributions look like, and where our observed statistic falls, we can use `visualize()`:

```{r visualize-indep, warning = FALSE, message = FALSE}
# visualize the null distribution and test statistic!
null_dist_sim %>%
  visualize() + 
  shade_p_value(observed_indep_statistic,
                direction = "greater")
```

We could also visualize the observed statistic against the theoretical null distribution. To do so, use the `assume()` verb to define a theoretical null distribution and then pass it to `visualize()` like a null distribution outputted from `generate()` and `calculate()`.

```{r visualize-indep-theor, warning = FALSE, message = FALSE}
# visualize the theoretical null distribution and test statistic!
gss %>%
  specify(college ~ finrela) %>%
  assume(distribution = "Chisq") %>%
  visualize() + 
  shade_p_value(observed_indep_statistic,
                direction = "greater")
```

To visualize both the randomization-based and theoretical null distributions to get a sense of how the two relate, we can pipe the randomization-based null distribution into `visualize()`, and further provide `method = "both"`.

```{r visualize-indep-both, warning = FALSE, message = FALSE}
# visualize both null distributions and the test statistic!
null_dist_sim %>%
  visualize(method = "both") + 
  shade_p_value(observed_indep_statistic,
                direction = "greater")
```

Either way, it looks like our observed test statistic would be quite unlikely if there were actually no association between education and income. More exactly, we can approximate the p-value with `get_p_value`:

```{r p-value-indep, warning = FALSE, message = FALSE}
# calculate the p value from the observed statistic and null distribution
p_value_independence <- null_dist_sim %>%
  get_p_value(obs_stat = observed_indep_statistic,
              direction = "greater")

p_value_independence
```

Thus, if there were really no relationship between education and income, our approximation of the probability that we would see a statistic as or more extreme than `r observed_indep_statistic` is approximately `r p_value_independence`.

To calculate the p-value using the true $\chi^2$ distribution, we can use the `pchisq` function from base R. This function allows us to situate the test statistic we calculated previously in the $\chi^2$ distribution with the appropriate degrees of freedom.

```{r}
pchisq(observed_indep_statistic$stat, 5, lower.tail = FALSE)
```

Note that, equivalently to the theory-based approach shown above, the package supplies a wrapper function, `chisq_test`, to carry out Chi-Squared tests of independence on tidy data. The syntax goes like this:

```{r chisq-indep-wrapper, message = FALSE, warning = FALSE}
chisq_test(gss, college ~ finrela)
```


### Goodness of Fit

Now, moving on to a chi-squared goodness of fit test, we'll take a look at the self-identified income class of our survey respondents. Suppose our null hypothesis is that `finrela` follows a uniform distribution (i.e. there's actually an equal number of people that describe their income as far below average, below average, average, above average, far above average, or that don't know their income.) The graph below represents this hypothesis:

```{r gof-plot, echo = FALSE}
gss %>%
  ggplot2::ggplot() +
  ggplot2::aes(x = finrela) +
  ggplot2::geom_bar() +
  ggplot2::geom_hline(yintercept = 466.3, col = "red") +
  ggplot2::labs(x = "finrela: Self-Identification of Income Class",
                y = "Number of Responses")
```

It seems like a uniform distribution may not be the most appropriate description of the data--many more people describe their income as average than than any of the other options. Lets now test whether this difference in distributions is statistically significant.

First, to carry out this hypothesis test, we would calculate our observed statistic.

```{r observed-gof-statistic, warning = FALSE, message = FALSE}
# calculating the null distribution
observed_gof_statistic <- gss %>%
  specify(response = finrela) %>%
  hypothesize(null = "point",
              p = c("far below average" = 1/6,
                    "below average" = 1/6,
                    "average" = 1/6,
                    "above average" = 1/6,
                    "far above average" = 1/6,
                    "DK" = 1/6)) %>%
  calculate(stat = "Chisq")
```

The observed statistic is `r observed_gof_statistic`. Now, generating a null distribution, by just dropping in a call to `generate()`:


```{r null-distribution-gof, warning = FALSE, message = FALSE}
# generating a null distribution, assuming each income class is equally likely
null_dist_gof <- gss %>%
  specify(response = finrela) %>%
  hypothesize(null = "point",
              p = c("far below average" = 1/6,
                    "below average" = 1/6,
                    "average" = 1/6,
                    "above average" = 1/6,
                    "far above average" = 1/6,
                    "DK" = 1/6)) %>%
  generate(reps = 1000, type = "draw") %>%
  calculate(stat = "Chisq")
```

Again, to get a sense for what these distributions look like, and where our observed statistic falls, we can use `visualize()`:

```{r visualize-indep-gof, warning = FALSE, message = FALSE}
# visualize the null distribution and test statistic!
null_dist_gof %>%
  visualize() + 
  shade_p_value(observed_gof_statistic,
                direction = "greater")
```

This statistic seems like it would be quite unlikely if income class self-identification actually followed a uniform distribution! How unlikely, though? Calculating the p-value:

```{r get-p-value-gof, warning = FALSE, message = FALSE}
# calculate the p-value
p_value_gof <- null_dist_gof %>%
  get_p_value(observed_gof_statistic,
              direction = "greater")

p_value_gof
```

Thus, if each self-identified income class was equally likely to occur, our approximation of the probability that we would see a distribution like the one we did is approximately `r p_value_gof`.

To calculate the p-value using the true $\chi^2$ distribution, we can use the `pchisq` function from base R. This function allows us to situate the test statistic we calculated previously in the $\chi^2$ distribution with the appropriate degrees of freedom.

```{r}
pchisq(observed_gof_statistic$stat, 5, lower.tail = FALSE)
```

Again, equivalently to the theory-based approach shown above, the package supplies a wrapper function, `chisq_test`, to carry out Chi-Squared goodness of fit tests on tidy data. The syntax goes like this:

```{r chisq-gof-wrapper, message = FALSE, warning = FALSE}
chisq_test(gss, 
           response = finrela,
           p = c("far below average" = 1/6,
                    "below average" = 1/6,
                    "average" = 1/6,
                    "above average" = 1/6,
                    "far above average" = 1/6,
                    "DK" = 1/6))
```


---
title: "Tidy t-Tests with infer"
description: "Conducting t-Tests on tidy data with infer."
output: rmarkdown::html_vignette
vignette: |
  %\VignetteIndexEntry{t-Tests}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r settings, include=FALSE}
knitr::opts_chunk$set(fig.width = 6, fig.height = 4.5) 
options(digits = 4)
```

```{r load-packages, echo = FALSE, message = FALSE, warning = FALSE}
library(ggplot2)
library(dplyr)
devtools::load_all()
```

### Introduction

In this vignette, we'll walk through conducting $t$-tests and their randomization-based analogue using `infer`. We'll start out with a 1-sample $t$-test, which compares a sample mean to a hypothesized true mean value. Then, we'll discuss paired $t$-tests, which are a special use case of 1-sample $t$-tests, and evaluate whether differences in paired values (e.g. some measure taken of a person before and after an experiment) differ from 0. Finally, we'll wrap up with 2-sample $t$-tests, testing the difference in means of two populations using a sample of data drawn from them.

Throughout this vignette, we'll make use of the `gss` dataset supplied by `infer`, which contains a sample of data from the General Social Survey. See `?gss` for more information on the variables included and their source. Note that this data (and our examples on it) are for demonstration purposes only, and will not necessarily provide accurate estimates unless weighted properly. For these examples, let's suppose that this dataset is a representative sample of a population we want to learn about: American adults. The data looks like this:

```{r glimpse-gss-actual, warning = FALSE, message = FALSE}
dplyr::glimpse(gss)
```

### 1-Sample t-Test

The 1-sample $t$-test can be used to test whether a sample of continuous data could have plausibly come from a population with a specified mean. 

As an example, we'll test whether the average American adult works 40 hours a week using data from the `gss`. To do so, we make use of the `hours` variable, giving the number of hours that respondents reported having worked in the previous week. The distribution of `hours` in the observed data looks like this:

```{r plot-1-sample, echo = FALSE}
gss %>%
  ggplot2::ggplot() +
  ggplot2::aes(x = hours) +
  ggplot2::geom_histogram(bins = 20) +
  ggplot2::labs(x = "hours: Number of Hours Worked",
                y = "Number of Responses") +
  ggplot2::scale_x_continuous(breaks = seq(0, 90, 10))
```

It looks like most respondents reported having worked 40 hours, but there's quite a bit of variability. Let's test whether we have evidence that the true mean number of hours that Americans work per week is 40.

infer's randomization-based analogue to the 1-sample $t$-test is a 1-sample mean test. We'll start off showcasing that test before demonstrating how to carry out a theory-based $t$-test with the package.

First, to calculate the observed statistic, we can use `specify()` and `calculate()`.

```{r calc-obs-stat-1-sample, warning = FALSE, message = FALSE}
# calculate the observed statistic
observed_statistic <- gss %>%
  specify(response = hours) %>%
  calculate(stat = "mean")
```

The observed statistic is `r observed_statistic`. Now, we want to compare this statistic to a null distribution, generated under the assumption that the mean was actually 40, to get a sense of how likely it would be for us to see this observed mean if the true number of hours worked per week in the population was really 40.

We can `generate` the null distribution using the bootstrap. In the bootstrap, for each replicate, a sample of size equal to the input sample size is drawn (with replacement) from the input sample data. This allows us to get a sense of how much variability we'd expect to see in the entire population so that we can then understand how unlikely our sample mean would be.

```{r generate-null-1-sample, warning = FALSE, message = FALSE}
# generate the null distribution
null_dist_1_sample <- gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  generate(reps = 1000, type = "bootstrap") %>%
  calculate(stat = "mean")
```

To get a sense for what these distributions look like, and where our observed statistic falls, we can use `visualize()`:

```{r visualize-1-sample, warning = FALSE, message = FALSE}
# visualize the null distribution and test statistic!
null_dist_1_sample %>%
  visualize() + 
  shade_p_value(observed_statistic,
                direction = "two-sided")
```

It looks like our observed mean of `r observed_statistic` would be relatively unlikely if the true mean was actually 40 hours a week. More exactly, we can calculate the p-value:

```{r p-value-1-sample, warning = FALSE, message = FALSE}
# calculate the p value from the test statistic and null distribution
p_value_1_sample <- null_dist_1_sample %>%
  get_p_value(obs_stat = observed_statistic,
              direction = "two-sided")

p_value_1_sample
```

Thus, if the true mean number of hours worked per week was really 40, our approximation of the probability that we would see a test statistic as or more extreme than `r observed_statistic` is approximately `r p_value_1_sample`.

Analogously to the steps shown above, the package supplies a wrapper function, `t_test`, to carry out 1-sample $t$-tests on tidy data. Rather than using randomization, the wrappers carry out the theory-based $t$-test. The syntax looks like this:

```{r t-test-wrapper, message = FALSE, warning = FALSE}
t_test(gss, response = hours, mu = 40)
```

An alternative approach to the `t_test()` wrapper is to calculate the observed statistic with an infer pipeline and then supply it to the `pt` function from base R. 

```{r}
# calculate the observed statistic
observed_statistic <- gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  calculate(stat = "t") %>%
  dplyr::pull()
```

Note that this pipeline to calculate an observed statistic includes a call to `hypothesize()` since the $t$ statistic requires a hypothesized mean value.

Then, juxtaposing that $t$ statistic with its associated distribution using the `pt` function:

```{r}
pt(observed_statistic, df = nrow(gss) - 1, lower.tail = FALSE)*2
```

Note that the resulting $t$-statistics from these two theory-based approaches are the same.

### Paired t-Test

You might be interested in running a paired $t$-test. Paired $t$-tests can be used in situations when there is a natural pairing between values in distributions---a common example would be two columns, `before` and `after`, say, that contain measurements from a patient before and after some treatment. To compare these two distributions, then, we're not necessarily interested in how the two distributions look different altogether, but how these two measurements from each individual change across time. (Pairings don't necessarily have to be over time; another common usage is measurements from two married people, for example.) Thus, we can create a new column (see `mutate()` from the `dplyr` package if you're not sure how to do this) that is the difference between the two: `difference = after - before`, and then examine _this_ distribution to see how each individuals' measurements changed over time.

Once we've `mutate()`d that new `difference` column, we can run a 1-sample $t$-test on it, where our null hypothesis is that `mu = 0` (i.e. the difference between these measurements before and after treatment is, on average, 0). To do so, we'd use the procedure outlined in the above section.


### 2-Sample t-Test

2-Sample $t$-tests evaluate the difference in mean values of two populations using data randomly-sampled from the population that approximately follows a normal distribution. As an example, we'll test if Americans work the same number of hours a week regardless of whether they have a college degree or not using data from the `gss`. The `college` and `hours` variables allow us to do so:

```{r plot-2-sample, echo = FALSE}
gss %>%
  ggplot2::ggplot() +
  ggplot2::aes(x = college, y = hours) +
  ggplot2::geom_boxplot() +
  ggplot2::labs(x = "college: Whether the Respondent has a College Degree",
                y = "hours: Number of Hours Worked")
```

It looks like both of these distributions are centered near 40 hours a week, but the distribution for those with a degree is slightly right skewed.

Again, note the warning about missing values---many respondents' values are missing. If we were actually carrying out this hypothesis test, we might look further into how this data was collected; it's possible that whether or not a value in either of these columns is missing is related to what that value would be. 

infer's randomization-based analogue to the 2-sample $t$-test is a difference in means test. We'll start off showcasing that test before demonstrating how to carry out a theory-based $t$-test with the package.

As with the one-sample test, to calculate the observed difference in means, we can use `specify()` and `calculate()`.

```{r calc-obs-stat-2-sample, warning = FALSE, message = FALSE}
# calculate the observed statistic
observed_statistic <- gss %>%
  specify(hours ~ college) %>%
  calculate(stat = "diff in means", order = c("degree", "no degree"))

observed_statistic
```

Note that, in the line `specify(hours ~ college)`, we could have swapped this out with the syntax `specify(response = hours, explanatory = college)`!

The `order` argument in that `calculate` line gives the order to subtract the mean values in: in our case, we're taking the mean number of hours worked by those with a degree minus the mean number of hours worked by those without a degree; a positive difference, then, would mean that people with degrees worked more than those without a degree.

Now, we want to compare this difference in means to a null distribution, generated under the assumption that the number of hours worked a week has no relationship with whether or not one has a college degree, to get a sense of how likely it would be for us to see this observed difference in means if there were really no relationship between these two variables.

We can `generate` the null distribution using permutation, where, for each replicate, each value of degree status will be randomly reassigned (without replacement) to a new number of hours worked per week in the sample in order to break any association between the two.

```{r generate-null-2-sample, warning = FALSE, message = FALSE}
# generate the null distribution with randomization
null_dist_2_sample <- gss %>%
  specify(hours ~ college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate(stat = "diff in means", order = c("degree", "no degree"))
```

Again, note that, in the lines `specify(hours ~ college)` in the above chunk, we could have used the syntax `specify(response = hours, explanatory = college)` instead!

To get a sense for what these distributions look like, and where our observed statistic falls, we can use `visualize()`.

```{r visualize-2-sample, warning = FALSE, message = FALSE}
# visualize the randomization-based null distribution and test statistic!
null_dist_2_sample %>%
  visualize() + 
  shade_p_value(observed_statistic,
                direction = "two-sided")
```

It looks like our observed statistic of `r observed_statistic` would be unlikely if there was truly no relationship between degree status and number of hours worked. More exactly, we can calculate the p-value; theoretical p-values are not yet supported, so we'll use the randomization-based null distribution to do calculate the p-value.

```{r p-value-2-sample, warning = FALSE, message = FALSE}
# calculate the p value from the randomization-based null 
# distribution and the observed statistic
p_value_2_sample <- null_dist_2_sample %>%
  get_p_value(obs_stat = observed_statistic,
              direction = "two-sided")

p_value_2_sample
```

Thus, if there were really no relationship between the number of hours worked a week and whether one has a college degree, the probability that we would see a statistic as or more extreme than `r observed_statistic` is approximately `r p_value_2_sample`.

Note that, similarly to the steps shown above, the package supplies a wrapper function, `t_test`, to carry out 2-sample $t$-tests on tidy data. The syntax looks like this:

```{r 2-sample-t-test-wrapper, message = FALSE, warning = FALSE}
t_test(x = gss, 
       formula = hours ~ college, 
       order = c("degree", "no degree"),
       alternative = "two-sided")
```

In the above example, we specified the relationship with the syntax `formula = hours ~ college`; we could have also written `response = hours, explanatory = college`.

An alternative approach to the `t_test()` wrapper is to calculate the observed statistic with an infer pipeline and then supply it to the `pt` function from base R. We can calculate the statistic as before, switching out the `stat = "diff in means"` argument with `stat = "t"`.

```{r}
# calculate the observed statistic
observed_statistic <- gss %>%
  specify(hours ~ college) %>%
  hypothesize(null = "point", mu = 40) %>%
  calculate(stat = "t", order = c("degree", "no degree")) %>%
  dplyr::pull()

observed_statistic
```

Note that this pipeline to calculate an observed statistic includes `hypothesize()` since the $t$ statistic requires a hypothesized mean value.

Then, juxtaposing that $t$ statistic with its associated distribution using the `pt` function:

```{r}
pt(observed_statistic, df = nrow(gss) - 2, lower.tail = FALSE)*2
```

Note that the results from these two theory-based approaches are nearly the same.
---
title: "Getting to Know infer"
description: "An introduction to the infer R package."
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: |
  %\VignetteIndexEntry{infer}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r include=FALSE}
knitr::opts_chunk$set(fig.width = 6, fig.height = 4.5) 
options(digits = 4)
```

### Introduction

`infer` implements an expressive grammar to perform statistical inference that coheres with the `tidyverse` design framework. Rather than providing methods for specific statistical tests, this package consolidates the principles that are shared among common hypothesis tests into a set of 4 main verbs (functions), supplemented with many utilities to visualize and extract value from their outputs.

Regardless of which hypothesis test we're using, we're still asking the same kind of question: is the effect/difference in our observed data real, or due to chance? To answer this question, we start by assuming that the observed data came from some world where "nothing is going on" (i.e. the observed effect was simply due to random chance), and call this assumption our *null hypothesis*. (In reality, we might not believe in the null hypothesis at all---the null hypothesis is in opposition to the *alternate hypothesis*, which supposes that the effect present in the observed data is actually due to the fact that "something is going on.") We then calculate a *test statistic* from our data that describes the observed effect. We can use this test statistic to calculate a *p-value*, giving the probability that our observed data could come about if the null hypothesis was true. If this probability is below some pre-defined *significance level* $\alpha$, then we can reject our null hypothesis.

The workflow of this package is designed around this idea. Starting out with some dataset,

+ `specify()` allows you to specify the variable, or relationship between variables, that you're interested in.
+ `hypothesize()` allows you to declare the null hypothesis.
+ `generate()` allows you to generate data reflecting the null hypothesis.
+ `calculate()` allows you to calculate a distribution of statistics from the generated data to form the null distribution.

Throughout this vignette, we make use of `gss`, a dataset supplied by `infer` containing a sample of 500 observations of 11 variables from the *General Social Survey*. 

```{r load-packages, echo = FALSE, message = FALSE, warning = FALSE}
library(dplyr)
devtools::load_all()
```


```{r load-gss, warning = FALSE, message = FALSE}
# load in the dataset
data(gss)

# take a look at its structure
dplyr::glimpse(gss)
```

Each row is an individual survey response, containing some basic demographic information on the respondent as well as some additional variables. See `?gss` for more information on the variables included and their source. Note that this data (and our examples on it) are for demonstration purposes only, and will not necessarily provide accurate estimates unless weighted properly. For these examples, let's suppose that this dataset is a representative sample of a population we want to learn about: American adults.

### specify(): Specifying Response (and Explanatory) Variables

The `specify` function can be used to specify which of the variables in the dataset you're interested in. If you're only interested in, say, the `age` of the respondents, you might write:

```{r specify-example, warning = FALSE, message = FALSE}
gss %>%
  specify(response = age)
```

On the front-end, the output of `specify` just looks like it selects off the columns in the dataframe that you've specified. Checking the class of this object, though:

```{r specify-one, warning = FALSE, message = FALSE}
gss %>%
  specify(response = age) %>%
  class()
```

We can see that the `infer` class has been appended on top of the dataframe classes--this new class stores some extra metadata.

If you're interested in two variables--`age` and `partyid`, for example--you can `specify` their relationship in one of two (equivalent) ways:

```{r specify-two, warning = FALSE, message = FALSE}
# as a formula
gss %>%
  specify(age ~ partyid)

# with the named arguments
gss %>%
  specify(response = age, explanatory = partyid)
```

If you're doing inference on one proportion or a difference in proportions, you will need to use the `success` argument to specify which level of your `response` variable is a success. For instance, if you're interested in the proportion of the population with a college degree, you might use the following code:

```{r specify-success, warning = FALSE, message = FALSE}
# specifying for inference on proportions
gss %>%
  specify(response = college, success = "degree")
```

### hypothesize(): Declaring the Null Hypothesis

The next step in the `infer` pipeline is often to declare a null hypothesis using `hypothesize()`. The first step is to supply one of "independence" or "point" to the `null` argument. If your null hypothesis assumes independence between two variables, then this is all you need to supply to `hypothesize()`:

```{r hypothesize-independence, warning = FALSE, message = FALSE}
gss %>%
  specify(college ~ partyid, success = "degree") %>%
  hypothesize(null = "independence")
```

If you're doing inference on a point estimate, you will also need to provide one of `p` (the true proportion of successes, between 0 and 1), `mu` (the true mean), `med` (the true median), or `sigma` (the true standard deviation). For instance, if the null hypothesis is that the mean number of hours worked per week in our population is 40, we would write:

```{r hypothesize-40-hr-week, warning = FALSE, message = FALSE}
gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40)
```

Again, from the front-end, the dataframe outputted from `hypothesize()` looks almost exactly the same as it did when it came out of `specify()`, but `infer` now "knows" your null hypothesis.

### generate(): Generating the Null Distribution

Once we've asserted our null hypothesis using `hypothesize()`, we can construct a null distribution based on this hypothesis. We can do this using one of several methods, supplied in the `type` argument:

* `bootstrap`: A bootstrap sample will be drawn for each replicate, where a sample of size equal to the input sample size is drawn (with replacement) from the input sample data.  
* `permute`: For each replicate, each input value will be randomly reassigned (without replacement) to a new output value in the sample.  
* `draw`: A value will be sampled from a theoretical distribution with parameters specified in `hypothesize()` for each replicate. This option is currently only applicable for testing point estimates. This generation type was previously called `"simulate"`, which has been superseded.

Continuing on with our example above, about the average number of hours worked a week, we might write:

```{r generate-point, warning = FALSE, message = FALSE}
set.seed(1)

gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  generate(reps = 1000, type = "bootstrap")
```

In the above example, we take 1000 bootstrap samples to form our null distribution.

Note that, before `generate()`ing, we've set the seed for random number generation with the `set.seed()` function. When using the infer package for research, or in other cases when exact reproducibility is a priority, this is good practice. infer will respect the random seed specified in the `set.seed()` function, returning the same result when `generate()`ing data given an identical seed.

To generate a null distribution for the independence of two variables, we could also randomly reshuffle the pairings of explanatory and response variables to break any existing association. For instance, to generate 1000 replicates that can be used to create a null distribution under the assumption that political party affiliation is not affected by age:

```{r generate-permute, warning = FALSE, message = FALSE}
gss %>%
  specify(partyid ~ age) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute")
```

### calculate(): Calculating Summary Statistics

`calculate()` calculates summary statistics from the output of infer core functions. The function takes in a `stat` argument, which is currently one of "mean", "median", "sum", "sd", "prop", "count", "diff in means", "diff in medians", "diff in props", "Chisq", "F", "t", "z", "slope", or "correlation". For example, continuing our example above to calculate the null distribution of mean hours worked per week:

```{r calculate-point, warning = FALSE, message = FALSE}
gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  generate(reps = 1000, type = "bootstrap") %>%
  calculate(stat = "mean")
```

The output of `calculate()` here shows us the sample statistic (in this case, the mean) for each of our 1000 replicates. If you're carrying out inference on differences in means, medians, or proportions, or t and z statistics, you will need to supply an `order` argument, giving the order in which the explanatory variables should be subtracted. For instance, to find the difference in mean age of those that have a college degree and those that don't, we might write:

```{r specify-diff-in-means, warning = FALSE, message = FALSE}
gss %>%
  specify(age ~ college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  calculate("diff in means", order = c("degree", "no degree"))
```

### Other Utilities

`infer` also offers several utilities to extract the meaning out of summary statistics and distributions---the package provides functions to visualize where a statistic is relative to a distribution (with `visualize()`), calculate p-values (with `get_p_value()`), and calculate confidence intervals (with `get_confidence_interval()`).

To illustrate, we'll go back to the example of determining whether the mean number of hours worked per week is 40 hours.

```{r utilities-examples}
# find the point estimate
obs_mean <- gss %>%
  specify(response = hours) %>%
  calculate(stat = "mean")

# generate a null distribution
null_dist <- gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  generate(reps = 1000, type = "bootstrap") %>%
  calculate(stat = "mean")
```

Our point estimate `r obs_mean` seems *pretty* close to 40, but a little bit different. We might wonder if this difference is just due to random chance, or if the mean number of hours worked per week in the population really isn't 40.

We could initially just visualize the null distribution.

```{r visualize, warning = FALSE, message = FALSE}
null_dist %>%
  visualize()
```

Where does our sample's observed statistic lie on this distribution? We can use the `obs_stat` argument to specify this.

```{r visualize2, warning = FALSE, message = FALSE}
null_dist %>%
  visualize() +
  shade_p_value(obs_stat = obs_mean, direction = "two-sided")
```

Notice that `infer` has also shaded the regions of the null distribution that are as (or more) extreme than our observed statistic. (Also, note that we now use the `+` operator to apply the `shade_p_value` function. This is because `visualize` outputs a plot object from `ggplot2` instead of a data frame, and the `+` operator is needed to add the p-value layer to the plot object.) The red bar looks like it's slightly far out on the right tail of the null distribution, so observing a sample mean of `r obs_mean` hours would be somewhat unlikely if the mean was actually 40 hours. How unlikely, though?

```{r get_p_value, warning = FALSE, message = FALSE}
# get a two-tailed p-value
p_value <- null_dist %>%
  get_p_value(obs_stat = obs_mean, direction = "two-sided")

p_value
```

It looks like the p-value is `r p_value`, which is pretty small---if the true mean number of hours worked per week was actually 40, the probability of our sample mean being this far (`r abs(obs_mean-40)` hours) from 40 would be `r p_value`. This may or may not be statistically significantly different, depending on the significance level $\alpha$ you decided on *before* you ran this analysis. If you had set $\alpha = .05$, then this difference would be statistically significant, but if you had set $\alpha = .01$, then it would not be.

To get a confidence interval around our estimate, we can write:

```{r get_conf, message = FALSE, warning = FALSE}
# generate a distribution like the null distribution, 
# though exclude the null hypothesis from the pipeline
boot_dist <- gss %>%
  specify(response = hours) %>%
  generate(reps = 1000, type = "bootstrap") %>%
  calculate(stat = "mean")

# start with the bootstrap distribution
ci <- boot_dist %>%
  # calculate the confidence interval around the point estimate
  get_confidence_interval(point_estimate = obs_mean,
                          # at the 95% confidence level
                          level = .95,
                          # using the standard error
                          type = "se")

ci
```

As you can see, 40 hours per week is not contained in this interval, which aligns with our previous conclusion that this finding is significant at the confidence level $\alpha = .05$. To see this interval represented visually, we can use the `shade_confidence_interval()` utility:

```{r visualize-ci, warning = FALSE, message = FALSE}
boot_dist %>%
  visualize() +
  shade_confidence_interval(endpoints = ci)
```

### Theoretical Methods

{infer} also provides functionality to use theoretical methods for `"Chisq"`, `"F"`, `"t"` and `"z"` distributions. 

Generally, to find a null distribution using theory-based methods, use the same code that you would use to find the observed statistic elsewhere, replacing calls to `calculate()` with `assume()`. For example, to calculate the observed $t$ statistic (a standardized mean):

```{r, message = FALSE, warning = FALSE}
# calculate an observed t statistic
obs_t <- gss %>%
  specify(response = hours) %>%
  hypothesize(null = "point", mu = 40) %>%
  calculate(stat = "t")
```

Then, to define a theoretical $t$ distribution, we could write:

```{r, message = FALSE, warning = FALSE}
# switch out calculate with assume to define a distribution
t_dist <- gss %>%
  specify(response = hours) %>%
  assume(distribution = "t")
```

From here, the theoretical distribution interfaces in the same way that simulation-based null distributions do. For example, to interface with p-values:

```{r, message = FALSE, warning = FALSE}
# visualize the theoretical null distribution
visualize(t_dist) +
  shade_p_value(obs_stat = obs_t, direction = "greater")

# more exactly, calculate the p-value
get_p_value(t_dist, obs_t, "greater")
```

Confidence intervals lie on the scale of the data rather than on the standardized scale of the theoretical distribution, so be sure to use the unstandardized observed statistic when working with confidence intervals.

```{r, message = FALSE, warning = FALSE}
# find the theory-based confidence interval
theor_ci <- 
  get_confidence_interval(
    x = t_dist,
    level = .95,
    point_estimate = obs_mean
  )

theor_ci
```

When visualized, the $t$ distribution will be recentered and rescaled to align with the scale of the observed data.

```{r}
# visualize the theoretical sampling distribution
visualize(t_dist) +
  shade_confidence_interval(theor_ci)
```

### Multiple regression

To accommodate randomization-based inference with multiple explanatory variables, the package implements an alternative workflow based on model fitting. Rather than `calculate()`ing statistics from resampled data, this side of the package allows you to `fit()` linear models on data resampled according to the null hypothesis, supplying model coefficients for each explanatory variable. For the most part, you can just switch out `calculate()` for `fit()` in your `calculate()`-based workflows.

As an example, suppose that we want to fit `hours` worked per week using the respondent `age` and `college` completion status. We could first begin by fitting a linear model to the observed data.

```{r}
observed_fit <- gss %>%
  specify(hours ~ age + college) %>%
  fit()
```

Now, to generate null distributions for each of these terms, we can fit 1000 models to resamples of the `gss` dataset, where the response `hours` is permuted in each. Note that this code is the same as the above except for the addition of the `hypothesize` and `generate` step.

```{r}
null_fits <- gss %>%
  specify(hours ~ age + college) %>%
  hypothesize(null = "independence") %>%
  generate(reps = 1000, type = "permute") %>%
  fit()

null_fits
```

To permute variables other than the response variable, the `variables` argument to `generate()` allows you to choose columns from the data to permute. Note that any derived effects that depend on these columns (e.g., interaction effects) will also be affected.

Beyond this point, observed fits and distributions from null fits interface exactly like analogous outputs from `calculate()`. For instance, we can use the following code to calculate a 95% confidence interval from these objects.

```{r}
get_confidence_interval(
  null_fits, 
  point_estimate = observed_fit, 
  level = .95
)
```

Or, we can shade p-values for each of these observed regression coefficients from the observed data.

```{r}
visualize(null_fits) + 
  shade_p_value(observed_fit, direction = "both")
```

### Conclusion

That's it! This vignette covers most all of the key functionality of infer. See `help(package = "infer")` for a full list of functions and vignettes.
---
title: 'infer: An R package for tidyverse-friendly statistical inference'
tags:
  - data science
  - tidyverse
  - inference
  - R
authors:
- name: Simon P. Couch
  orcid: 0000-0001-5676-5107
  affiliation: "1, 2"
- name: Andrew P. Bray
  orcid: 0000-0002-4037-7414
  affiliation: 3
- name: Chester Ismay
  orcid: 0000-0003-2820-2547
  affiliation: 4
- name: Evgeni Chasnovski
  orcid: 0000-0002-1617-4019
  affiliation: 5
- name: Benjamin S. Baumer
  orcid: 0000-0002-3279-0516
  affiliation: 6
- name: Mine Çetinkaya-Rundel
  orcid: 0000-0001-6452-2420
  affiliation: "2, 7"
affiliations:
 - name: Johns Hopkins, Department of Biostatistics
   index: 1
 - name: RStudio
   index: 2
 - name: UC Berkeley, Department of Statistics and Reed College Mathematics Department (on leave)
   index: 3
 - name: Flatiron School
   index: 4
 - name: No Affiliation
   index: 5
 - name: Smith College, Program in Statistical & Data Sciences
   index: 6
 - name: Duke University, Department of Statistical Science
   index: 7

citation_author: Couch et. al.
date: 12 June 2021
year: 2021
bibliography: paper.bib
output: 
  rticles::joss_article:
    keep_tex: true
    includes:
      in_header: columns.tex
csl: apa.csl
journal: JOSS
---

# Summary

`infer` implements an expressive grammar to perform statistical inference that adheres to the `tidyverse` design framework [@wickham2019welcome]. Rather than providing methods for specific statistical tests, this package consolidates the principles that are shared among common hypothesis tests and confidence intervals into a set of four main verbs (functions), supplemented with many utilities to visualize and extract value from their outputs.

# Statement of Need

Packages implementing methods for basic statistical inference in R are highly variable in their interfaces. The structure of inputted data, argument names, expected argument types, argument orders, output types, and spelling cases varies widely both within and among packages. This diversity in approaches obscures the intuition shared among common inferential procedures, makes details of usage difficult to remember, and prevents an expressive and idiomatic coding style.

`infer` is an R package for randomization-based hypothesis testing, naturalizing an intuitive understanding of statistical inference via a unified and expressive grammar. Four functions provide functionality encompassing a large swath of basic frequentist statistical inference, abstracting away details of specific tests and shifting the focus of the analyst to the observed data and the processes that generated it. Such a grammar lends itself to applications in teaching, data pedagogy research, applied scientific research, and advanced predictive modeling. For one, the principled approach of the `infer` package has made it an especially good fit for teaching introductory statistics and data science [@ismay2019statistical; @baumer2020teaching; @cetinkaya2021fresh] and research in data pedagogy [@fergusson2021introducing; @loy2021bringing]. Further, the package has already seen usage in a number of published scientific applications [@mclean2021controlled; @ask2021per; @fallon2021single]. Finally, the package integrates with the greater tidymodels collection of packages, a burgeoning software ecosystem for tidyverse-aligned predictive modeling used across many modern research and industrial applications [@kuhn2020tidymodels]. To date, the package has been downloaded more than 400,000 times.

# Underlying Principles

Regardless of the hypothesis test in question, an analyst asks the same kind of question when conducting statistical inference: is the effect/difference in the observed data real, or due to random chance? To answer this question, the analyst begins by assuming that the effect in the observed data was simply due to random chance, and calls this assumption the *null hypothesis*. (In reality, they might not believe in the null hypothesis at all---the null hypothesis is in opposition to the *alternate hypothesis*, which supposes that the effect present in the observed data is actually due to the fact that "something is going on.") The analyst then calculates a *test statistic* from the data that describes the observed effect. They can use this test statistic to calculate a *p-value* via juxtaposition with a *null distribution*, giving the probability that the observed data could come about if the null hypothesis were true. If this probability is below some pre-defined *significance level* $\alpha$, then the analyst can reject the null hypothesis.

The workflow of this package is designed around this idea. Starting out with some dataset,

+ `specify()` allows the analyst to specify the variable, or relationship between variables, that they are interested in.
+ `hypothesize()` allows the analyst to declare the null hypothesis.
+ `generate()` allows the analyst to generate data reflecting the null hypothesis or using the bootstrap.
+ `calculate()` allows the analyst to calculate summary statistics, either from
     * the observed data, to form the observed test statistic.
     * data `generate()`d to reflect the null hypothesis, to form a randomization-based null distribution of test statistics.

As such, the ultimate output of an infer pipeline using these four functions is generally an _observed statistic_ or _null distribution_ of test statistics. These four functions are thus supplemented with several utilities to visualize and extract value from their outputs.

+ `visualize()` plots the null distribution of test statistics.
     * `shade_p_value()` situates the observed statistic in the null distribution, shading the region as or more extreme.
+ `get_p_value()` calculates a p-value via the juxtaposition of the test statistic and the null distribution.

The workflow outlined above can also be used for constructing confidence intervals via bootstrapping with the omission of the `hypothesize()` step in the pipeline. The resulting bootstrap distribution can then be visualized with `visualize()`, the confidence interval region can be situated in the bootstrap distribution with `shade_confidence_interval()`, and the bounds of the confidence interval can be calculated with `get_confidence_interval()`.

Beyond this, the `infer` package offers:

* methods for inference using theory-based distributions
* shorthand wrappers for common statistical tests using tidy data
* model-fitting workflows to accommodate multiple explanatory variables

# Comparison to Other Packages

Several software packages on the Comprehensive R Archive Network share functionality with `infer` [@CRAN]. `broom` and `parameters` convert model objects to unified output formats, though they do not provide methods for fitting models, describing null distributions, performing bootstrapping, or calculating summary statistics from tabular data [@r-broom; @r-parameters]. `statsExpressions`, and adjacent packages in the `easystats` ecosystem, implement wrappers with consistent interfaces for theory-based hypothesis tests [@r-statsExpressions]. Similarly, `mosaic` is a package used to teach statistics by unifying summary statistics, visualization, and modeling with a consistent API built around R's formula interface. The `mosaic` package also includes functionality to conduct randomization-based inference [@r-mosaic]. At a higher level, though, the structure of each of these packages is defined by model types and statistics, where each model type or statistic has its own associated function and/or object class. In contrast, `infer` is structured around four functions, situating statistics and model types within a more abstracted grammar.^[This grammar follows from Allen Downey's "there is only one test" framework [@downey2016].] 

# Acknowledgements

We acknowledge contributions from Albert Y. Kim, Jo Hardin, Jay Lee, Amelia McNamara, Nick Solomon, and Richie Cotton.

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{fit.infer}
\alias{fit.infer}
\title{Fit linear models to infer objects}
\usage{
\method{fit}{infer}(object, ...)
}
\arguments{
\item{object}{Output from an infer function---likely \code{\link[=generate]{generate()}} or
\code{\link[=specify]{specify()}}---which specifies the formula and data to fit a model to.}

\item{...}{Any optional arguments to pass along to the model fitting
function. See \code{\link[stats:glm]{stats::glm()}} for more information.}
}
\value{
A \link[tibble:tibble]{tibble} containing the following columns:

\itemize{
\item \code{replicate}: Only supplied if the input object had been previously
passed to \code{\link[=generate]{generate()}}. A number corresponding to which resample of the
original data set the model was fitted to.
\item \code{term}: The explanatory variable (or intercept) in question.
\item \code{estimate}: The model coefficient for the given resample (\code{replicate}) and
explanatory variable (\code{term}).
}
}
\description{
Given the output of an infer core function, this function will fit
a linear model using \code{\link[stats:glm]{stats::glm()}} according to the formula and data supplied
earlier in the pipeline. If passed the output of \code{\link[=specify]{specify()}} or
\code{\link[=hypothesize]{hypothesize()}}, the function will fit one model. If passed the output
of \code{\link[=generate]{generate()}}, it will fit a model to each data resample, denoted in
the \code{replicate} column. The family of the fitted model depends on the type
of the response variable. If the response is numeric, \code{fit()} will use
\code{family = "gaussian"} (linear regression). If the response is a 2-level
factor or character, \code{fit()} will use \code{family = "binomial"} (logistic
regression). To fit character or factor response variables with more than
two levels, we recommend \code{\link[parsnip:multinom_reg]{parsnip::multinom_reg()}}.

infer provides a fit "method" for infer objects, which is a way of carrying
out model fitting as applied to infer output. The "generic," imported from
the generics package and re-exported from this package, provides the
general form of \code{fit()} that points to infer's method when called on an
infer object. That generic is also documented here.

Learn more in \code{vignette("infer")}.
}
\details{
Randomization-based statistical inference with multiple explanatory
variables requires careful consideration of the null hypothesis in question
and its implications for permutation procedures. Inference for partial
regression coefficients via the permutation method implemented in
\code{\link[=generate]{generate()}} for multiple explanatory variables, consistent with its meaning
elsewhere in the package, is subject to additional distributional assumptions
beyond those required for one explanatory variable. Namely, the distribution
of the response variable must be similar to the distribution of the errors
under the null hypothesis' specification of a fixed effect of the explanatory
variables. (This null hypothesis is reflected in the \code{variables} argument to
\code{\link[=generate]{generate()}}. By default, all of the explanatory variables are treated
as fixed.) A general rule of thumb here is, if there are large outliers
in the distributions of any of the explanatory variables, this distributional
assumption will not be satisfied; when the response variable is permuted,
the (presumably outlying) value of the response will no longer be paired
with the outlier in the explanatory variable, causing an outsize effect
on the resulting slope coefficient for that explanatory variable.

More sophisticated methods that are outside of the scope of this package
requiring fewer---or less strict---distributional assumptions
exist. For an overview, see "Permutation tests for univariate or
multivariate analysis of variance and regression" (Marti J. Anderson,
2001), \doi{10.1139/cjfas-58-3-626}.
}
\section{Reproducibility}{
When using the infer package for research, or in other cases when exact
reproducibility is a priority, be sure the set the seed for R’s random
number generator. infer will respect the random seed specified in the
\code{set.seed()} function, returning the same result when \code{generate()}ing
data given an identical seed. For instance, we can calculate the
difference in mean \code{age} by \code{college} degree status using the \code{gss}
dataset from 10 versions of the \code{gss} resampled with permutation using
the following code.\if{html}{\out{<div class="r">}}\preformatted{set.seed(1)

gss \%>\%
  specify(age ~ college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 5, type = "permute") \%>\%
  calculate("diff in means", order = c("degree", "no degree"))
}\if{html}{\out{</div>}}\preformatted{## Response: age (numeric)
## Explanatory: college (factor)
## Null Hypothesis: independence
## # A tibble: 5 × 2
##   replicate   stat
##       <int>  <dbl>
## 1         1 -0.531
## 2         2 -2.35 
## 3         3  0.764
## 4         4  0.280
## 5         5  0.350
}

Setting the seed to the same value again and rerunning the same code
will produce the same result.\if{html}{\out{<div class="r">}}\preformatted{# set the seed
set.seed(1)

gss \%>\%
  specify(age ~ college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 5, type = "permute") \%>\%
  calculate("diff in means", order = c("degree", "no degree"))
}\if{html}{\out{</div>}}\preformatted{## Response: age (numeric)
## Explanatory: college (factor)
## Null Hypothesis: independence
## # A tibble: 5 × 2
##   replicate   stat
##       <int>  <dbl>
## 1         1 -0.531
## 2         2 -2.35 
## 3         3  0.764
## 4         4  0.280
## 5         5  0.350
}

Please keep this in mind when writing infer code that utilizes
resampling with \code{generate()}.
}

\examples{
# fit a linear model predicting number of hours worked per
# week using respondent age and degree status.
observed_fit <- gss \%>\%
  specify(hours ~ age + college) \%>\%
  fit()

observed_fit

# fit 100 models to resamples of the gss dataset, where the response 
# `hours` is permuted in each. note that this code is the same as 
# the above except for the addition of the `generate` step.
null_fits <- gss \%>\%
  specify(hours ~ age + college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 100, type = "permute") \%>\%
  fit()

null_fits

# for logistic regression, just supply a binary response variable!
# (this can also be made explicit via the `family` argument in ...)
gss \%>\%
  specify(college ~ age + hours) \%>\%
  fit()

# more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}  

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesize.R
\name{hypothesize}
\alias{hypothesize}
\alias{hypothesise}
\title{Declare a null hypothesis}
\usage{
hypothesize(x, null, p = NULL, mu = NULL, med = NULL, sigma = NULL)

hypothesise(x, null, p = NULL, mu = NULL, med = NULL, sigma = NULL)
}
\arguments{
\item{x}{A data frame that can be coerced into a \link[tibble:tibble]{tibble}.}

\item{null}{The null hypothesis. Options include \code{"independence"} and
\code{"point"}.}

\item{p}{The true proportion of successes (a number between 0 and 1). To be used with point null hypotheses when the specified response
variable is categorical.}

\item{mu}{The true mean (any numerical value). To be used with point null
hypotheses when the specified response variable is continuous.}

\item{med}{The true median (any numerical value). To be used with point null
hypotheses when the specified response variable is continuous.}

\item{sigma}{The true standard deviation (any numerical value). To be used with
point null hypotheses.}
}
\value{
A tibble containing the response (and explanatory, if specified)
variable data with parameter information stored as well.
}
\description{
Declare a null hypothesis about variables selected in \code{\link[=specify]{specify()}}.

Learn more in \code{vignette("infer")}.
}
\examples{
# hypothesize independence of two variables
gss \%>\%
 specify(college ~ partyid, success = "degree") \%>\%
 hypothesize(null = "independence")
 
# hypothesize a mean number of hours worked per week of 40
gss \%>\%
  specify(response = hours) \%>\%
  hypothesize(null = "point", mu = 40)

# more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}

}
\seealso{
Other core functions: 
\code{\link{calculate}()},
\code{\link{generate}()},
\code{\link{specify}()}
}
\concept{core functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assume.R
\name{assume}
\alias{assume}
\title{Define a theoretical distribution}
\usage{
assume(x, distribution, df = NULL, ...)
}
\arguments{
\item{x}{The output of \code{\link[=specify]{specify()}} or \code{\link[=hypothesize]{hypothesize()}}, giving the
observed data, variable(s) of interest, and (optionally) null hypothesis.}

\item{distribution}{The distribution in question, as a string. One of
\code{"F"}, \code{"Chisq"}, \code{"t"}, or \code{"z"}.}

\item{df}{Optional. The degrees of freedom parameter(s) for the \code{distribution}
supplied, as a numeric vector. For \code{distribution = "F"}, this should have
length two (e.g. \code{c(10, 3)}). For \code{distribution = "Chisq"} or
\code{distribution = "t"}, this should have length one. For
\code{distribution = "z"}, this argument is not required. The package
will supply a message if the supplied \code{df} argument is different from
recognized values. See the Details section below for more information.}

\item{...}{Currently ignored.}
}
\value{
An infer theoretical distribution that can be passed to helpers
like \code{\link[=visualize]{visualize()}}, \code{\link[=get_p_value]{get_p_value()}}, and \code{\link[=get_confidence_interval]{get_confidence_interval()}}.
}
\description{
This function allows the user to define a null distribution based on
theoretical methods. In many infer pipelines, \code{assume()} can be
used in place of \code{\link[=generate]{generate()}} and \code{\link[=calculate]{calculate()}} to create a null
distribution. Rather than outputting a data frame containing a
distribution of test statistics calculated from resamples of the observed
data, \code{assume()} outputs a more abstract type of object just containing
the distributional details supplied in the \code{distribution} and \code{df} arguments.
However, \code{assume()} output can be passed to \code{\link[=visualize]{visualize()}}, \code{\link[=get_p_value]{get_p_value()}},
and \code{\link[=get_confidence_interval]{get_confidence_interval()}} in the same way that simulation-based
distributions can.

To define a theoretical null distribution (for use in hypothesis testing),
be sure to provide a null hypothesis via \code{\link[=hypothesize]{hypothesize()}}. To define a
theoretical sampling distribution (for use in confidence intervals),
provide the output of \code{\link[=specify]{specify()}}. Sampling distributions (only
implemented for \code{t} and \code{z}) lie on the scale of the data, and will be
recentered and rescaled to match the corresponding \code{stat} given in
\code{\link[=calculate]{calculate()}} to calculate the observed statistic.
}
\details{
Note that the assumption being expressed here, for use in theory-based
inference, only extends to \emph{distributional} assumptions: the null
distribution in question and its parameters. Statistical inference with
infer, whether carried out via simulation (i.e. based on pipelines
using \code{\link[=generate]{generate()}} and \code{\link[=calculate]{calculate()}}) or theory (i.e. with \code{assume()}),
always involves the condition that observations are independent of
each other.

\code{infer} only supports theoretical tests on one or two means via the
\code{t} distribution and one or two proportions via the \code{z}.

For tests comparing two means, if \code{n1} is the group size for one level of
the explanatory variable, and \code{n2} is that for the other level, \code{infer}
will recognize the following degrees of freedom (\code{df}) arguments:
\itemize{
\item \code{min(n1 - 1, n2 - 1)}
\item \code{n1 + n2 - 2}
\item The \code{"parameter"} entry of the analogous \code{stats::t.test()} call
\item The \code{"parameter"} entry of the analogous \code{stats::t.test()} call with \code{var.equal = TRUE}
}

By default, the package will use the \code{"parameter"} entry of the analogous
\code{stats::t.test()} call with \code{var.equal = FALSE} (the default).
}
\examples{
  
# construct theoretical distributions ---------------------------------

# F distribution
# with the `partyid` explanatory variable
gss \%>\% 
  specify(age ~ partyid) \%>\% 
  assume(distribution = "F")

# Chi-squared goodness of fit distribution
# on the `finrela` variable
gss \%>\%
  specify(response = finrela) \%>\%
  hypothesize(null = "point",
              p = c("far below average" = 1/6,
                    "below average" = 1/6,
                    "average" = 1/6,
                    "above average" = 1/6,
                    "far above average" = 1/6,
                    "DK" = 1/6)) \%>\%
  assume("Chisq")

# Chi-squared test of independence
# on the `finrela` and `sex` variables
gss \%>\%
  specify(formula = finrela ~ sex) \%>\%
  assume(distribution = "Chisq")

# T distribution
gss \%>\% 
  specify(age ~ college) \%>\%
  assume("t")

# Z distribution
gss \%>\%
  specify(response = sex, success = "female") \%>\%
  assume("z")

\dontrun{
# each of these distributions can be passed to infer helper
# functions alongside observed statistics!

# for example, a 1-sample t-test -------------------------------------

# calculate the observed statistic 
obs_stat <- gss \%>\%
  specify(response = hours) \%>\%
  hypothesize(null = "point", mu = 40) \%>\%
  calculate(stat = "t")

# construct a null distribution
null_dist <- gss \%>\%
  specify(response = hours) \%>\%
  assume("t")

# juxtapose them visually
visualize(null_dist) + 
  shade_p_value(obs_stat, direction = "both")
  
# calculate a p-value
get_p_value(null_dist, obs_stat, direction = "both")

# or, an F test ------------------------------------------------------

# calculate the observed statistic 
obs_stat <- gss \%>\% 
  specify(age ~ partyid) \%>\%
  hypothesize(null = "independence") \%>\%
  calculate(stat = "F")

# construct a null distribution
null_dist <- gss \%>\% 
  specify(age ~ partyid) \%>\%
  assume(distribution = "F")

# juxtapose them visually
visualize(null_dist) + 
  shade_p_value(obs_stat, direction = "both")
  
# calculate a p-value
get_p_value(null_dist, obs_stat, direction = "both")
}
    
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shade_confidence_interval.R
\name{shade_confidence_interval}
\alias{shade_confidence_interval}
\alias{shade_ci}
\title{Add information about confidence interval}
\usage{
shade_confidence_interval(
  endpoints,
  color = "mediumaquamarine",
  fill = "turquoise",
  ...
)

shade_ci(endpoints, color = "mediumaquamarine", fill = "turquoise", ...)
}
\arguments{
\item{endpoints}{The lower and upper bounds of the interval to be plotted.
Likely, this will be the output of \code{\link[=get_confidence_interval]{get_confidence_interval()}}.
For \code{\link[=calculate]{calculate()}}-based workflows, this will be a 2-element vector
or a \verb{1 x 2} data frame containing the lower and upper values to be plotted.
For \code{\link[=fit.infer]{fit()}}-based workflows, a \verb{(p + 1) x 3} data frame
with columns \code{term}, \code{lower_ci}, and \code{upper_ci}, giving the upper and
lower bounds for each regression term. For use in visualizations of
\code{\link[=assume]{assume()}} output, this must be the output of \code{\link[=get_confidence_interval]{get_confidence_interval()}}.}

\item{color}{A character or hex string specifying the color of the
end points as a vertical lines on the plot.}

\item{fill}{A character or hex string specifying the color to shade the
confidence interval. If \code{NULL} then no shading is actually done.}

\item{...}{Other arguments passed along to \\{ggplot2\\} functions.}
}
\value{
If added to an existing {infer} visualization, a \\{ggplot2\\}
object displaying the supplied intervals on top of its corresponding
distribution. Otherwise, an \code{infer_layer} list.
}
\description{
\code{shade_confidence_interval()} plots a confidence interval region on top of
\code{\link[=visualize]{visualize()}} output. The output is a ggplot2 layer that can be added with
\code{+}. The function has a shorter alias, \code{shade_ci()}.

Learn more in \code{vignette("infer")}.
}
\examples{
# find the point estimate---mean number of hours worked per week
point_estimate <- gss \%>\%
  specify(response = hours) \%>\%
  calculate(stat = "mean")
  
# ...and a bootstrap distribution
boot_dist <- gss \%>\%
  # ...we're interested in the number of hours worked per week
  specify(response = hours) \%>\%
  # generating data points
  generate(reps = 1000, type = "bootstrap") \%>\%
  # finding the distribution from the generated data
  calculate(stat = "mean")
  
# find a confidence interval around the point estimate
ci <- boot_dist \%>\%
  get_confidence_interval(point_estimate = point_estimate,
                          # at the 95\% confidence level
                          level = .95,
                          # using the standard error method
                          type = "se")   
  
  
# and plot it!
boot_dist \%>\%
  visualize() +
  shade_confidence_interval(ci)
  
# or just plot the bounds
boot_dist \%>\%
  visualize() +
  shade_confidence_interval(ci, fill = NULL)
  
# you can shade confidence intervals on top of
# theoretical distributions, too---the theoretical
# distribution will be recentered and rescaled to
# align with the confidence interval
sampling_dist <- gss \%>\%
  specify(response = hours) \%>\%
  assume(distribution = "t") 
  
visualize(sampling_dist) +
  shade_confidence_interval(ci)

\donttest{
# to visualize distributions of coefficients for multiple
# explanatory variables, use a `fit()`-based workflow

# fit 1000 linear models with the `hours` variable permuted
null_fits <- gss \%>\%
 specify(hours ~ age + college) \%>\%
 hypothesize(null = "independence") \%>\%
 generate(reps = 1000, type = "permute") \%>\%
 fit()
 
null_fits

# fit a linear model to the observed data
obs_fit <- gss \%>\%
  specify(hours ~ age + college) \%>\%
  fit()

obs_fit

# get confidence intervals for each term
conf_ints <- 
  get_confidence_interval(
    null_fits, 
    point_estimate = obs_fit, 
    level = .95
  )

# visualize distributions of coefficients 
# generated under the null
visualize(null_fits)

# add a confidence interval shading layer to juxtapose 
# the null fits with the observed fit for each term
visualize(null_fits) + 
  shade_confidence_interval(conf_ints)
}

# more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}

}
\seealso{
Other visualization functions: 
\code{\link{shade_p_value}()}
}
\concept{visualization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{visualize}
\alias{visualize}
\alias{visualise}
\title{Visualize statistical inference}
\usage{
visualize(data, bins = 15, method = "simulation", dens_color = "black", ...)

visualise(data, bins = 15, method = "simulation", dens_color = "black", ...)
}
\arguments{
\item{data}{A distribution. For simulation-based inference, a data frame
containing a distribution of \code{\link[=calculate]{calculate()}}d statistics
or \code{\link[=fit.infer]{fit()}}ted coefficient estimates. This object should
have been passed to \code{\link[=generate]{generate()}} before being supplied or
\code{\link[=calculate]{calculate()}} to \code{\link[=fit.infer]{fit()}}. For theory-based inference,
the output of \code{\link[=assume]{assume()}}.}

\item{bins}{The number of bins in the histogram.}

\item{method}{A string giving the method to display. Options are
\code{"simulation"}, \code{"theoretical"}, or \code{"both"} with \code{"both"} corresponding to
\code{"simulation"} and \code{"theoretical"}. If \code{data} is the output of \code{\link[=assume]{assume()}},
this argument will be ignored and default to \code{"theoretical"}.}

\item{dens_color}{A character or hex string specifying the color of the
theoretical density curve.}

\item{...}{Other arguments passed along to \\{ggplot2\\} functions.}
}
\value{
For \code{\link[=calculate]{calculate()}}-based workflows, a ggplot showing the simulation-based
distribution as a histogram or bar graph. Can also be used to display
theoretical distributions.

For \code{\link[=assume]{assume()}}-based workflows, a ggplot showing the theoretical distribution.

For \code{\link[=fit.infer]{fit()}}-based workflows, a \code{patchwork} object
showing the simulation-based distributions as a histogram or bar graph.
The interface to adjust plot options and themes is a bit different
for \code{patchwork} plots than ggplot2 plots. The examples highlight the
biggest differences here, but see \code{\link[patchwork:plot_annotation]{patchwork::plot_annotation()}} and
\link[patchwork:plot_arithmetic]{patchwork::&.gg} for more details.
}
\description{
Visualize the distribution of the simulation-based inferential statistics or
the theoretical distribution (or both!).

Learn more in \code{vignette("infer")}.
}
\details{
In order to make the visualization workflow more straightforward
and explicit, \code{visualize()} now only should be used to plot distributions
of statistics directly. A number of arguments related to shading p-values and
confidence intervals are now deprecated in \code{visualize()} and should
now be passed to \code{\link[=shade_p_value]{shade_p_value()}} and \code{\link[=shade_confidence_interval]{shade_confidence_interval()}},
respectively. \code{\link[=visualize]{visualize()}} will raise a warning if deprecated arguments
are supplied.
}
\examples{
  
# generate a null distribution
null_dist <- gss \%>\%
  # we're interested in the number of hours worked per week
  specify(response = hours) \%>\%
  # hypothesizing that the mean is 40
  hypothesize(null = "point", mu = 40) \%>\%
  # generating data points for a null distribution
  generate(reps = 1000, type = "bootstrap") \%>\%
  # calculating a distribution of means
  calculate(stat = "mean")
  
# or a bootstrap distribution, omitting the hypothesize() step,
# for use in confidence intervals
boot_dist <- gss \%>\%
  specify(response = hours) \%>\%
  generate(reps = 1000, type = "bootstrap") \%>\%
  calculate(stat = "mean")
  
# we can easily plot the null distribution by piping into visualize
null_dist \%>\%
  visualize()

# we can add layers to the plot as in ggplot, as well... 
# find the point estimate---mean number of hours worked per week
point_estimate <- gss \%>\%
  specify(response = hours) \%>\%
  calculate(stat = "mean")
  
# find a confidence interval around the point estimate
ci <- boot_dist \%>\%
  get_confidence_interval(point_estimate = point_estimate,
                          # at the 95\% confidence level
                          level = .95,
                          # using the standard error method
                          type = "se")  
  
# display a shading of the area beyond the p-value on the plot
null_dist \%>\%
  visualize() +
  shade_p_value(obs_stat = point_estimate, direction = "two-sided")

# ...or within the bounds of the confidence interval
null_dist \%>\%
  visualize() +
  shade_confidence_interval(ci)
  
# plot a theoretical sampling distribution by creating
# a theory-based distribution with `assume()`
sampling_dist <- gss \%>\%
  specify(response = hours) \%>\%
  assume(distribution = "t") 
  
visualize(sampling_dist)

# you can shade confidence intervals on top of
# theoretical distributions, too---the theoretical
# distribution will be recentered and rescaled to
# align with the confidence interval
visualize(sampling_dist) +
  shade_confidence_interval(ci)


# to plot both a theory-based and simulation-based null distribution,
# use a theorized statistic (i.e. one of t, z, F, or Chisq)
# and supply the simulation-based null distribution
null_dist_t <- gss \%>\%
  specify(response = hours) \%>\%
  hypothesize(null = "point", mu = 40) \%>\%
  generate(reps = 1000, type = "bootstrap") \%>\%
  calculate(stat = "t")
  
obs_stat <- gss \%>\%
  specify(response = hours) \%>\%
  hypothesize(null = "point", mu = 40) \%>\%
  calculate(stat = "t")

visualize(null_dist_t, method = "both")

visualize(null_dist_t, method = "both") +
  shade_p_value(obs_stat, "both")

\donttest{
# to visualize distributions of coefficients for multiple
# explanatory variables, use a `fit()`-based workflow

# fit 1000 models with the `hours` variable permuted
null_fits <- gss \%>\%
 specify(hours ~ age + college) \%>\%
 hypothesize(null = "independence") \%>\%
 generate(reps = 1000, type = "permute") \%>\%
 fit()
 
null_fits

# visualize distributions of resulting coefficients
visualize(null_fits)

# the interface to add themes and other elements to patchwork
# plots (outputted by `visualize` when the inputted data
# is from the `fit()` function) is a bit different than adding
# them to ggplot2 plots.
library(ggplot2)

# to add a ggplot2 theme to a `calculate()`-based visualization, use `+`
null_dist \%>\% visualize() + theme_dark()
  
# to add a ggplot2 theme to a `fit()`-based visualization, use `&`
null_fits \%>\% visualize() & theme_dark()
}

# More in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}

}
\seealso{
\code{\link[=shade_p_value]{shade_p_value()}}, \code{\link[=shade_confidence_interval]{shade_confidence_interval()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer.R
\docType{package}
\name{infer}
\alias{infer}
\alias{infer-package}
\title{infer: a grammar for statistical inference}
\description{
\if{html}{\figure{infer.png}{options: align='right'}}
}
\details{
The objective of this package is to perform statistical inference using a
grammar that illustrates the underlying concepts and a format that coheres
with the tidyverse.

For an overview of how to use the core functionality, see \code{vignette("infer")}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/tidymodels/infer}
  \item \url{https://infer.tidymodels.org/}
  \item Report bugs at \url{https://github.com/tidymodels/infer/issues}
}

}
\author{
\strong{Maintainer}: Andrew Bray \email{abray@reed.edu}

Authors:
\itemize{
  \item Chester Ismay \email{chester.ismay@gmail.com} (\href{https://orcid.org/0000-0003-2820-2547}{ORCID})
  \item Evgeni Chasnovski \email{evgeni.chasnovski@gmail.com} (\href{https://orcid.org/0000-0002-1617-4019}{ORCID})
  \item Simon Couch \email{scouch2@jhu.edu} (\href{https://orcid.org/0000-0001-5676-5107}{ORCID})
  \item Ben Baumer \email{ben.baumer@gmail.com} (\href{https://orcid.org/0000-0002-3279-0516}{ORCID})
  \item Mine Cetinkaya-Rundel \email{mine@stat.duke.edu} (\href{https://orcid.org/0000-0001-6452-2420}{ORCID})
}

Other contributors:
\itemize{
  \item Ted Laderas \email{tedladeras@gmail.com} (\href{https://orcid.org/0000-0002-6207-7068}{ORCID}) [contributor]
  \item Nick Solomon \email{nick.solomon@datacamp.com} [contributor]
  \item Johanna Hardin \email{Jo.Hardin@pomona.edu} [contributor]
  \item Albert Y. Kim \email{albert.ys.kim@gmail.com} (\href{https://orcid.org/0000-0001-7824-306X}{ORCID}) [contributor]
  \item Neal Fultz \email{nfultz@gmail.com} [contributor]
  \item Doug Friedman \email{doug.nhp@gmail.com} [contributor]
  \item Richie Cotton \email{richie@datacamp.com} (\href{https://orcid.org/0000-0003-2504-802X}{ORCID}) [contributor]
  \item Brian Fannin \email{captain@pirategrunt.com} [contributor]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_p_value.R
\name{get_p_value}
\alias{get_p_value}
\alias{get_p_value.default}
\alias{get_pvalue}
\alias{get_p_value.infer_dist}
\title{Compute p-value}
\usage{
get_p_value(x, obs_stat, direction)

\method{get_p_value}{default}(x, obs_stat, direction)

get_pvalue(x, obs_stat, direction)

\method{get_p_value}{infer_dist}(x, obs_stat, direction)
}
\arguments{
\item{x}{A null distribution. For simulation-based inference, a data frame
containing a distribution of \code{\link[=calculate]{calculate()}}d statistics
or \code{\link[=fit.infer]{fit()}}ted coefficient estimates. This object should
have been passed to \code{\link[=generate]{generate()}} before being supplied or
\code{\link[=calculate]{calculate()}} to \code{\link[=fit.infer]{fit()}}. For theory-based inference,
the output of \code{\link[=assume]{assume()}}.}

\item{obs_stat}{A data frame containing the observed statistic (in a
\code{\link[=calculate]{calculate()}}-based workflow) or observed fit (in a
\code{\link[=fit.infer]{fit()}}-based workflow). This object is likely the output
of \code{\link[=calculate]{calculate()}} or \code{\link[=fit.infer]{fit()}} and need not
to have been passed to \code{\link[=generate]{generate()}}.}

\item{direction}{A character string. Options are \code{"less"}, \code{"greater"}, or
\code{"two-sided"}. Can also use \code{"left"}, \code{"right"}, \code{"both"},
\code{"two_sided"}, or \code{"two sided"}, \code{"two.sided"}.}
}
\value{
A \link[tibble:tibble]{tibble} containing the following columns:

\itemize{
\item \code{term}: The explanatory variable (or intercept) in question. Only
supplied if the input had been previously passed to \code{\link[=fit.infer]{fit()}}.
\item \code{p_value}: A value in [0, 1] giving the probability that a
statistic/coefficient as or more extreme than the observed
statistic/coefficient would occur if the null hypothesis were true.
}
}
\description{
Compute a p-value from a null distribution and observed statistic.

Learn more in \code{vignette("infer")}.
}
\section{Aliases}{

\code{get_pvalue()} is an alias of \code{get_p_value()}.
\code{p_value} is a deprecated alias of \code{get_p_value()}.
}

\section{Zero p-value}{

Though a true p-value of 0 is impossible, \code{get_p_value()} may return 0 in
some cases. This is due to the simulation-based nature of the \{infer\}
package; the output of this function is an approximation based on
the number of \code{reps} chosen in the \code{generate()} step. When the observed
statistic is very unlikely given the null hypothesis, and only a small
number of \code{reps} have been generated to form a null distribution,
it is possible that the observed statistic will be more extreme than
every test statistic generated to form the null distribution, resulting
in an approximate p-value of 0. In this case, the true p-value is a small
value likely less than \code{3/reps} (based on a poisson approximation).

In the case that a p-value of zero is reported, a warning message will be
raised to caution the user against reporting a p-value exactly equal to 0.
}

\examples{

# using a simulation-based null distribution ------------------------------

# find the point estimate---mean number of hours worked per week
point_estimate <- gss \%>\%
  specify(response = hours) \%>\%
  calculate(stat = "mean")

# starting with the gss dataset
gss \%>\%
  # ...we're interested in the number of hours worked per week
  specify(response = hours) \%>\%
  # hypothesizing that the mean is 40
  hypothesize(null = "point", mu = 40) \%>\%
  # generating data points for a null distribution
  generate(reps = 1000, type = "bootstrap") \%>\%
  # finding the null distribution
  calculate(stat = "mean") \%>\%
  get_p_value(obs_stat = point_estimate, direction = "two-sided")
  
# using a theoretical null distribution -----------------------------------

# calculate the observed statistic 
obs_stat <- gss \%>\%
  specify(response = hours) \%>\%
  hypothesize(null = "point", mu = 40) \%>\%
  calculate(stat = "t")

# define a null distribution
null_dist <- gss \%>\%
  specify(response = hours) \%>\%
  assume("t")

# calculate a p-value
get_p_value(null_dist, obs_stat, direction = "both")

# using a model fitting workflow -----------------------------------------

# fit a linear model predicting number of hours worked per
# week using respondent age and degree status.
observed_fit <- gss \%>\%
  specify(hours ~ age + college) \%>\%
  fit()

observed_fit

# fit 100 models to resamples of the gss dataset, where the response 
# `hours` is permuted in each. note that this code is the same as 
# the above except for the addition of the `generate` step.
null_fits <- gss \%>\%
  specify(hours ~ age + college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 100, type = "permute") \%>\%
  fit()

null_fits

get_p_value(null_fits, obs_stat = observed_fit, direction = "two-sided")
  
# more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}  
  
}
\seealso{
Other auxillary functions: 
\code{\link{get_confidence_interval}()}
}
\concept{auxillary functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_confidence_interval.R
\name{get_confidence_interval}
\alias{get_confidence_interval}
\alias{get_ci}
\title{Compute confidence interval}
\usage{
get_confidence_interval(x, level = 0.95, type = NULL, point_estimate = NULL)

get_ci(x, level = 0.95, type = NULL, point_estimate = NULL)
}
\arguments{
\item{x}{A distribution. For simulation-based inference, a data frame
containing a distribution of \code{\link[=calculate]{calculate()}}d statistics
or \code{\link[=fit.infer]{fit()}}ted coefficient estimates. This object should
have been passed to \code{\link[=generate]{generate()}} before being supplied or
\code{\link[=calculate]{calculate()}} to \code{\link[=fit.infer]{fit()}}. For theory-based inference,
output of \code{\link[=assume]{assume()}}. Distributions for confidence intervals do not
require a null hypothesis via \code{\link[=hypothesize]{hypothesize()}}.}

\item{level}{A numerical value between 0 and 1 giving the confidence level.
Default value is 0.95.}

\item{type}{A string giving which method should be used for creating the
confidence interval. The default is \code{"percentile"} with \code{"se"}
corresponding to (multiplier * standard error) and \code{"bias-corrected"} for
bias-corrected interval as other options.}

\item{point_estimate}{A data frame containing the observed statistic (in a
\code{\link[=calculate]{calculate()}}-based workflow) or observed fit (in a
\code{\link[=fit.infer]{fit()}}-based workflow). This object is likely the output
of \code{\link[=calculate]{calculate()}} or \code{\link[=fit.infer]{fit()}} and need not
to have been passed to \code{\link[=generate]{generate()}}. Set to \code{NULL} by
default. Must be provided if \code{type} is \code{"se"} or \code{"bias-corrected"}.}
}
\value{
A \link[tibble:tibble]{tibble} containing the following columns:

\itemize{
\item \code{term}: The explanatory variable (or intercept) in question. Only
supplied if the input had been previously passed to \code{\link[=fit.infer]{fit()}}.
\item \code{lower_ci}, \code{upper_ci}: The lower and upper bounds of the confidence
interval, respectively.
}
}
\description{
Compute a confidence interval around a summary statistic. Both
simulation-based and theoretical methods are supported, though only
\code{type = "se"} is supported for theoretical methods.

Learn more in \code{vignette("infer")}.
}
\details{
A null hypothesis is not required to compute a confidence interval. However,
including \code{\link[=hypothesize]{hypothesize()}} in a pipeline leading to \code{get_confidence_interval()}
will not break anything. This can be useful when computing a confidence
interval using the same distribution used to compute a p-value.

Theoretical confidence intervals (i.e. calculated by supplying the output
of \code{\link[=assume]{assume()}} to the \code{x} argument) require that the point estimate lies on
the scale of the data. The distribution defined in \code{\link[=assume]{assume()}} will be
recentered and rescaled to align with the point estimate, as can be shown
in the output of \code{\link[=visualize]{visualize()}} when paired with \code{\link[=shade_confidence_interval]{shade_confidence_interval()}}.
Confidence intervals are implemented for the following distributions and
point estimates:

\itemize{
\item \code{distribution = "t"}: \code{point_estimate} should be the output of
\code{\link[=calculate]{calculate()}} with \code{stat = "mean"} or \code{stat = "diff in means"}
\item \code{distribution = "z"}: \code{point_estimate} should be the output of
\code{\link[=calculate]{calculate()}} with \code{stat = "prop"} or \code{stat = "diff in props"}
}
}
\section{Aliases}{

\code{get_ci()} is an alias of \code{get_confidence_interval()}.
\code{conf_int()} is a deprecated alias of \code{get_confidence_interval()}.
}

\examples{

boot_dist <- gss \%>\%
  # We're interested in the number of hours worked per week
  specify(response = hours) \%>\%
  # Generate bootstrap samples
  generate(reps = 1000, type = "bootstrap") \%>\%
  # Calculate mean of each bootstrap sample
  calculate(stat = "mean")

boot_dist \%>\%
  # Calculate the confidence interval around the point estimate
  get_confidence_interval(
    # At the 95\% confidence level; percentile method
    level = 0.95
  )

# for type = "se" or type = "bias-corrected" we need a point estimate
sample_mean <- gss \%>\%
  specify(response = hours) \%>\%
  calculate(stat = "mean")

boot_dist \%>\%
  get_confidence_interval(
    point_estimate = sample_mean,
    # At the 95\% confidence level
    level = 0.95,
    # Using the standard error method
    type = "se"
  )
  
# using a theoretical distribution -----------------------------------

# define a sampling distribution
sampling_dist <- gss \%>\%
  specify(response = hours) \%>\%
  assume("t")

# get the confidence interval---note that the
# point estimate is required here
get_confidence_interval(
  sampling_dist, 
  level = .95, 
  point_estimate = sample_mean
)
  
# using a model fitting workflow -----------------------

# fit a linear model predicting number of hours worked per
# week using respondent age and degree status.
observed_fit <- gss \%>\%
  specify(hours ~ age + college) \%>\%
  fit()

observed_fit

# fit 100 models to resamples of the gss dataset, where the response 
# `hours` is permuted in each. note that this code is the same as 
# the above except for the addition of the `generate` step.
null_fits <- gss \%>\%
  specify(hours ~ age + college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 100, type = "permute") \%>\%
  fit()

null_fits

get_confidence_interval(
  null_fits, 
  point_estimate = observed_fit, 
  level = .95
)

# more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}

}
\seealso{
Other auxillary functions: 
\code{\link{get_p_value}()}
}
\concept{auxillary functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.R
\name{generate}
\alias{generate}
\title{Generate resamples, permutations, or simulations}
\usage{
generate(x, reps = 1, type = NULL, variables = !!response_expr(x), ...)
}
\arguments{
\item{x}{A data frame that can be coerced into a \link[tibble:tibble]{tibble}.}

\item{reps}{The number of resamples to generate.}

\item{type}{The method used to generate resamples of the observed
data reflecting the null hypothesis. Currently one of
\code{"bootstrap"}, \code{"permute"}, or \code{"draw"} (see below).}

\item{variables}{If \code{type = "permute"}, a set of unquoted column names in the
data to permute (independently of each other). Defaults to only the
response variable. Note that any derived effects that depend on these
columns (e.g., interaction effects) will also be affected.}

\item{...}{Currently ignored.}
}
\value{
A tibble containing \code{reps} generated datasets, indicated by the
\code{replicate} column.
}
\description{
Generation creates a simulated distribution from \code{specify()}.
In the context of confidence intervals, this is a bootstrap distribution
based on the result of \code{specify()}. In the context of hypothesis testing,
this is a null distribution based on the result of \code{specify()} and
\verb{hypothesize().}

Learn more in \code{vignette("infer")}.
}
\section{Generation Types}{


The \code{type} argument determines the method used to create the null
distribution.

\itemize{
\item \code{bootstrap}: A bootstrap sample will be drawn for each replicate,
where a sample of size equal to the input sample size is drawn (with
replacement) from the input sample data.
\item \code{permute}: For each replicate, each input value will be randomly
reassigned (without replacement) to a new output value in the sample.
\item \code{draw}: A value will be sampled from a theoretical distribution
with parameter \code{p} specified in \code{\link[=hypothesize]{hypothesize()}} for each replicate. This
option is currently only applicable for testing on one proportion. This
generation type was previously called \code{"simulate"}, which has been
superseded.
}
}

\section{Reproducibility}{
When using the infer package for research, or in other cases when exact
reproducibility is a priority, be sure the set the seed for R’s random
number generator. infer will respect the random seed specified in the
\code{set.seed()} function, returning the same result when \code{generate()}ing
data given an identical seed. For instance, we can calculate the
difference in mean \code{age} by \code{college} degree status using the \code{gss}
dataset from 10 versions of the \code{gss} resampled with permutation using
the following code.\if{html}{\out{<div class="r">}}\preformatted{set.seed(1)

gss \%>\%
  specify(age ~ college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 5, type = "permute") \%>\%
  calculate("diff in means", order = c("degree", "no degree"))
}\if{html}{\out{</div>}}\preformatted{## Response: age (numeric)
## Explanatory: college (factor)
## Null Hypothesis: independence
## # A tibble: 5 × 2
##   replicate   stat
##       <int>  <dbl>
## 1         1 -0.531
## 2         2 -2.35 
## 3         3  0.764
## 4         4  0.280
## 5         5  0.350
}

Setting the seed to the same value again and rerunning the same code
will produce the same result.\if{html}{\out{<div class="r">}}\preformatted{# set the seed
set.seed(1)

gss \%>\%
  specify(age ~ college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 5, type = "permute") \%>\%
  calculate("diff in means", order = c("degree", "no degree"))
}\if{html}{\out{</div>}}\preformatted{## Response: age (numeric)
## Explanatory: college (factor)
## Null Hypothesis: independence
## # A tibble: 5 × 2
##   replicate   stat
##       <int>  <dbl>
## 1         1 -0.531
## 2         2 -2.35 
## 3         3  0.764
## 4         4  0.280
## 5         5  0.350
}

Please keep this in mind when writing infer code that utilizes
resampling with \code{generate()}.
}

\examples{
# generate a null distribution by taking 200 bootstrap samples
gss \%>\%
 specify(response = hours) \%>\%
 hypothesize(null = "point", mu = 40) \%>\%
 generate(reps = 200, type = "bootstrap")

# generate a null distribution for the independence of
# two variables by permuting their values 200 times
gss \%>\%
 specify(partyid ~ age) \%>\%
 hypothesize(null = "independence") \%>\%
 generate(reps = 200, type = "permute")

# generate a null distribution via sampling from a
# binomial distribution 200 times
gss \%>\%
specify(response = sex, success = "female") \%>\%
  hypothesize(null = "point", p = .5) \%>\%
  generate(reps = 200, type = "draw") \%>\%
  calculate(stat = "z")

# more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}

}
\seealso{
Other core functions: 
\code{\link{calculate}()},
\code{\link{hypothesize}()},
\code{\link{specify}()}
}
\concept{core functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate.R
\name{calculate}
\alias{calculate}
\title{Calculate summary statistics}
\usage{
calculate(
  x,
  stat = c("mean", "median", "sum", "sd", "prop", "count", "diff in means",
    "diff in medians", "diff in props", "Chisq", "F", "slope", "correlation", "t", "z",
    "ratio of props", "odds ratio"),
  order = NULL,
  ...
)
}
\arguments{
\item{x}{The output from \code{\link[=generate]{generate()}} for computation-based inference or the
output from \code{\link[=hypothesize]{hypothesize()}} piped in to here for theory-based inference.}

\item{stat}{A string giving the type of the statistic to calculate. Current
options include \code{"mean"}, \code{"median"}, \code{"sum"}, \code{"sd"}, \code{"prop"}, \code{"count"},
\code{"diff in means"}, \code{"diff in medians"}, \code{"diff in props"}, \code{"Chisq"} (or
\code{"chisq"}), \code{"F"} (or \code{"f"}), \code{"t"}, \code{"z"}, \code{"ratio of props"}, \code{"slope"},
\code{"odds ratio"}, and \code{"correlation"}. \code{infer} only supports theoretical
tests on one or two means via the \code{"t"} distribution and one or two
proportions via the \code{"z"}.}

\item{order}{A string vector of specifying the order in which the levels of
the explanatory variable should be ordered for subtraction (or division
for ratio-based statistics), where \code{order = c("first", "second")} means
\code{("first" - "second")}, or the analogue for ratios. Needed for inference on
difference in means, medians, proportions, ratios, t, and z statistics.}

\item{...}{To pass options like \code{na.rm = TRUE} into functions like
\link[base:mean]{mean()}, \link[stats:sd]{sd()}, etc. Can also be used to
supply hypothesized null values for the \code{"t"} statistic.}
}
\value{
A tibble containing a \code{stat} column of calculated statistics.
}
\description{
Given the output of \code{\link[=specify]{specify()}} and/or \code{\link[=hypothesize]{hypothesize()}}, this function will
return the observed statistic specified with the \code{stat} argument. Some test
statistics, such as \code{Chisq}, \code{t}, and \code{z}, require a null hypothesis. If
provided the output of \code{\link[=generate]{generate()}}, the function will calculate the
supplied \code{stat} for each \code{replicate}.

Learn more in \code{vignette("infer")}.
}
\section{Missing levels in small samples}{

In some cases, when bootstrapping with small samples, some generated
bootstrap samples will have only one level of the explanatory variable
present. For some test statistics, the calculated statistic in these
cases will be NaN. The package will omit non-finite values from
visualizations (with a warning) and raise an error in p-value calculations.
}

\section{Reproducibility}{
When using the infer package for research, or in other cases when exact
reproducibility is a priority, be sure the set the seed for R’s random
number generator. infer will respect the random seed specified in the
\code{set.seed()} function, returning the same result when \code{generate()}ing
data given an identical seed. For instance, we can calculate the
difference in mean \code{age} by \code{college} degree status using the \code{gss}
dataset from 10 versions of the \code{gss} resampled with permutation using
the following code.\if{html}{\out{<div class="r">}}\preformatted{set.seed(1)

gss \%>\%
  specify(age ~ college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 5, type = "permute") \%>\%
  calculate("diff in means", order = c("degree", "no degree"))
}\if{html}{\out{</div>}}\preformatted{## Response: age (numeric)
## Explanatory: college (factor)
## Null Hypothesis: independence
## # A tibble: 5 × 2
##   replicate   stat
##       <int>  <dbl>
## 1         1 -0.531
## 2         2 -2.35 
## 3         3  0.764
## 4         4  0.280
## 5         5  0.350
}

Setting the seed to the same value again and rerunning the same code
will produce the same result.\if{html}{\out{<div class="r">}}\preformatted{# set the seed
set.seed(1)

gss \%>\%
  specify(age ~ college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 5, type = "permute") \%>\%
  calculate("diff in means", order = c("degree", "no degree"))
}\if{html}{\out{</div>}}\preformatted{## Response: age (numeric)
## Explanatory: college (factor)
## Null Hypothesis: independence
## # A tibble: 5 × 2
##   replicate   stat
##       <int>  <dbl>
## 1         1 -0.531
## 2         2 -2.35 
## 3         3  0.764
## 4         4  0.280
## 5         5  0.350
}

Please keep this in mind when writing infer code that utilizes
resampling with \code{generate()}.
}

\examples{

# calculate a null distribution of hours worked per week under
# the null hypothesis that the mean is 40
gss \%>\%
  specify(response = hours) \%>\%
  hypothesize(null = "point", mu = 40) \%>\%
  generate(reps = 200, type = "bootstrap") \%>\%
  calculate(stat = "mean")
  
# calculate the corresponding observed statistic
gss \%>\%
  specify(response = hours) \%>\%
  calculate(stat = "mean")

# calculate a null distribution assuming independence between age
# of respondent and whether they have a college degree
gss \%>\%
  specify(age ~ college) \%>\%
  hypothesize(null = "independence") \%>\%
  generate(reps = 200, type = "permute") \%>\%
  calculate("diff in means", order = c("degree", "no degree"))
  
# calculate the corresponding observed statistic
gss \%>\%
  specify(age ~ college) \%>\%
  calculate("diff in means", order = c("degree", "no degree"))
  
# some statistics require a null hypothesis
 gss \%>\%
   specify(response = hours) \%>\% 
   hypothesize(null = "point", mu = 40) \%>\%
   calculate(stat = "t")
   
# more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}

}
\seealso{
\code{\link[=visualize]{visualize()}}, \code{\link[=get_p_value]{get_p_value()}}, and \code{\link[=get_confidence_interval]{get_confidence_interval()}}
to extract value from this function's outputs.

Other core functions: 
\code{\link{generate}()},
\code{\link{hypothesize}()},
\code{\link{specify}()}
}
\concept{core functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specify.R
\name{specify}
\alias{specify}
\title{Specify response and explanatory variables}
\usage{
specify(x, formula, response = NULL, explanatory = NULL, success = NULL)
}
\arguments{
\item{x}{A data frame that can be coerced into a \link[tibble:tibble]{tibble}.}

\item{formula}{A formula with the response variable on the left and the
explanatory on the right. Alternatively, a \code{response} and \code{explanatory}
argument can be supplied.}

\item{response}{The variable name in \code{x} that will serve as the response.
This is an alternative to using the \code{formula} argument.}

\item{explanatory}{The variable name in \code{x} that will serve as the
explanatory variable. This is an alternative to using the formula argument.}

\item{success}{The level of \code{response} that will be considered a success, as
a string. Needed for inference on one proportion, a difference in
proportions, and corresponding z stats.}
}
\value{
A tibble containing the response (and explanatory, if specified)
variable data.
}
\description{
\code{specify()} is used to specify which columns in the supplied data frame are
the relevant response (and, if applicable, explanatory) variables. Note that
character variables are converted to \code{factor}s.

Learn more in \code{vignette("infer")}.
}
\examples{
# specifying for a point estimate on one variable
gss \%>\%
   specify(response = age)

# specify a relationship between variables as a formula...
gss \%>\%
  specify(age ~ partyid)
  
# ...or with named arguments!
gss \%>\%
  specify(response = age, explanatory = partyid)

# more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}

}
\seealso{
Other core functions: 
\code{\link{calculate}()},
\code{\link{generate}()},
\code{\link{hypothesize}()}
}
\concept{core functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{chisq_stat}
\alias{chisq_stat}
\title{Tidy chi-squared test statistic}
\usage{
chisq_stat(x, formula, response = NULL, explanatory = NULL, ...)
}
\arguments{
\item{x}{A data frame that can be coerced into a \link[tibble:tibble]{tibble}.}

\item{formula}{A formula with the response variable on the left and the
explanatory on the right.}

\item{response}{The variable name in \code{x} that will serve as the response.
This is alternative to using the \code{formula} argument.}

\item{explanatory}{The variable name in \code{x} that will serve as the
explanatory variable.}

\item{...}{Additional arguments for \link[stats:chisq.test]{chisq.test()}.}
}
\description{
@description
}
\details{
A shortcut wrapper function to get the observed test statistic for a chisq
test. Uses \link[stats:chisq.test]{chisq.test()}, which applies a continuity
correction. This function has been deprecated in favor of the more
general \code{\link[=observe]{observe()}}.
}
\examples{
# chi-squared test statistic for test of independence
# of college completion status depending and one's
# self-identified income class
chisq_stat(gss, college ~ finrela)

# chi-squared test statistic for a goodness of fit
# test on whether self-identified income class
# follows a uniform distribution
chisq_stat(gss,
           response = finrela,
           p = c("far below average" = 1/6,
                 "below average" = 1/6,
                 "average" = 1/6,
                 "above average" = 1/6,
                 "far above average" = 1/6,
                 "DK" = 1/6))

}
\seealso{
Other wrapper functions: 
\code{\link{chisq_test}()},
\code{\link{observe}()},
\code{\link{prop_test}()},
\code{\link{t_stat}()},
\code{\link{t_test}()}

Other functions for calculating observed statistics: 
\code{\link{observe}()},
\code{\link{t_stat}()}
}
\concept{functions for calculating observed statistics}
\concept{wrapper functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{prop_test}
\alias{prop_test}
\title{Tidy proportion test}
\usage{
prop_test(
  x,
  formula,
  response = NULL,
  explanatory = NULL,
  p = NULL,
  order = NULL,
  alternative = "two-sided",
  conf_int = TRUE,
  conf_level = 0.95,
  success = NULL,
  correct = NULL,
  z = FALSE,
  ...
)
}
\arguments{
\item{x}{A data frame that can be coerced into a \link[tibble:tibble]{tibble}.}

\item{formula}{A formula with the response variable on the left and the
explanatory on the right, where an explanatory variable NULL indicates
a test of a single proportion.}

\item{response}{The variable name in \code{x} that will serve as the response.
This is alternative to using the \code{formula} argument. This is an alternative
to the formula interface.}

\item{explanatory}{The variable name in \code{x} that will serve as the
explanatory variable. Optional. This is an alternative to the formula
interface.}

\item{p}{A numeric vector giving the hypothesized null proportion of
success for each group.}

\item{order}{A string vector specifying the order in which the proportions
should be subtracted, where  \code{order = c("first", "second")} means
\code{"first" - "second"}. Ignored for one-sample tests, and optional for two
sample tests.}

\item{alternative}{Character string giving the direction of the alternative
hypothesis. Options are \code{"two-sided"} (default), \code{"greater"}, or \code{"less"}.
Only used when testing the null that a single proportion equals a given
value, or that two proportions are equal; ignored otherwise.}

\item{conf_int}{A logical value for whether to report the confidence
interval or not. \code{TRUE} by default, ignored if \code{p} is specified for a
two-sample test. Only used when testing the null that a single
proportion equals a given value, or that two proportions are equal;
ignored otherwise.}

\item{conf_level}{A numeric value between 0 and 1. Default value is 0.95.
Only used when testing the null that a single proportion equals a given
value, or that two proportions are equal; ignored otherwise.}

\item{success}{The level of \code{response} that will be considered a success, as
a string. Only used when testing the null that a single
proportion equals a given value, or that two proportions are equal;
ignored otherwise.}

\item{correct}{A logical indicating whether Yates' continuity correction
should be applied where possible. If \code{z = TRUE}, the \code{correct} argument will
be overwritten as \code{FALSE}. Otherwise defaults to \code{correct = TRUE}.}

\item{z}{A logical value for whether to report the statistic as a standard
normal deviate or a Pearson's chi-square statistic. \eqn{z^2}  is distributed
chi-square with 1 degree of freedom, though note that the user will likely
need to turn off Yates' continuity correction by setting \code{correct = FALSE}
to see this connection.}

\item{...}{Additional arguments for \link[stats:prop.test]{prop.test()}.}
}
\description{
A tidier version of \link[stats:prop.test]{prop.test()} for equal or given
proportions.
}
\examples{
# two-sample proportion test for difference in proportions of
# college completion by respondent sex
prop_test(gss,
          college ~ sex,
          order = c("female", "male"))

# one-sample proportion test for hypothesized null
# proportion of college completion of .2
prop_test(gss,
          college ~ NULL,
          p = .2)

# report as a z-statistic rather than chi-square
# and specify the success level of the response
prop_test(gss,
          college ~ NULL,
          success = "degree",
          p = .2,
          z = TRUE)

}
\seealso{
Other wrapper functions: 
\code{\link{chisq_stat}()},
\code{\link{chisq_test}()},
\code{\link{observe}()},
\code{\link{t_stat}()},
\code{\link{t_test}()}
}
\concept{wrapper functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_methods.R
\name{print.infer}
\alias{print.infer}
\alias{print.infer_layer}
\alias{print.infer_dist}
\title{Print methods}
\usage{
\method{print}{infer}(x, ...)

\method{print}{infer_layer}(x, ...)

\method{print}{infer_dist}(x, ...)
}
\arguments{
\item{x}{An object of class \code{infer}, i.e. output from \code{\link[=specify]{specify()}} or
\code{\link[=hypothesize]{hypothesize()}}, or of class \code{infer_layer}, i.e. output from
\code{\link[=shade_p_value]{shade_p_value()}} or \code{\link[=shade_confidence_interval]{shade_confidence_interval()}}.}

\item{...}{Arguments passed to methods.}
}
\description{
Print methods
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shade_p_value.R
\name{shade_p_value}
\alias{shade_p_value}
\alias{shade_pvalue}
\title{Shade histogram area beyond an observed statistic}
\usage{
shade_p_value(obs_stat, direction, color = "red2", fill = "pink", ...)

shade_pvalue(obs_stat, direction, color = "red2", fill = "pink", ...)
}
\arguments{
\item{obs_stat}{The observed statistic or estimate. For
\code{\link[=calculate]{calculate()}}-based workflows, this will be a 1-element numeric vector or
a \verb{1 x 1} data frame containing the observed statistic.
For \code{\link[=fit.infer]{fit()}}-based workflows, a \verb{(p + 1) x 2} data frame
with columns \code{term} and \code{estimate} giving the observed estimate for
each term.}

\item{direction}{A string specifying in which direction the shading should
occur. Options are \code{"less"}, \code{"greater"}, or \code{"two-sided"}. Can
also give \code{"left"}, \code{"right"}, \code{"both"}, \code{"two_sided"}, \code{"two sided"},
or \code{"two.sided"}. If \code{NULL}, the function will not shade any area.}

\item{color}{A character or hex string specifying the color of the observed
statistic as a vertical line on the plot.}

\item{fill}{A character or hex string specifying the color to shade the
p-value region. If \code{NULL}, the function will not shade any area.}

\item{...}{Other arguments passed along to \\{ggplot2\\} functions.
For expert use only.}
}
\value{
If added to an existing {infer} visualization, a \\{ggplot2\\}
object displaying the supplied statistic on top of its corresponding
distribution. Otherwise, an \code{infer_layer} list.
}
\description{
\code{shade_p_value()} plots a p-value region on top of
\code{\link[=visualize]{visualize()}} output. The output is a ggplot2 layer that can be added with
\code{+}. The function has a shorter alias, \code{shade_pvalue()}.

Learn more in \code{vignette("infer")}.
}
\examples{
# find the point estimate---mean number of hours worked per week
point_estimate <- gss \%>\%
  specify(response = hours) \%>\%
  hypothesize(null = "point", mu = 40) \%>\%
  calculate(stat = "t")
  
# ...and a null distribution
null_dist <- gss \%>\%
  # ...we're interested in the number of hours worked per week
  specify(response = hours) \%>\%
  # hypothesizing that the mean is 40
  hypothesize(null = "point", mu = 40) \%>\%
  # generating data points for a null distribution
  generate(reps = 1000, type = "bootstrap") \%>\%
  # estimating the null distribution
  calculate(stat = "t")
  
# shade the p-value of the point estimate
null_dist \%>\%
  visualize() +
  shade_p_value(obs_stat = point_estimate, direction = "two-sided")
  
# you can shade confidence intervals on top of
# theoretical distributions, too!
null_dist_theory <- gss \%>\%
  specify(response = hours) \%>\%
  assume(distribution = "t") 
  
null_dist_theory \%>\%
  visualize() +
  shade_p_value(obs_stat = point_estimate, direction = "two-sided")

\donttest{
# to visualize distributions of coefficients for multiple
# explanatory variables, use a `fit()`-based workflow

# fit 1000 linear models with the `hours` variable permuted
null_fits <- gss \%>\%
 specify(hours ~ age + college) \%>\%
 hypothesize(null = "independence") \%>\%
 generate(reps = 1000, type = "permute") \%>\%
 fit()
 
null_fits

# fit a linear model to the observed data
obs_fit <- gss \%>\%
  specify(hours ~ age + college) \%>\%
  fit()

obs_fit

# visualize distributions of coefficients 
# generated under the null
visualize(null_fits)

# add a p-value shading layer to juxtapose the null 
# fits with the observed fit for each term
visualize(null_fits) + 
  shade_p_value(obs_fit, direction = "both")

# the direction argument will be applied 
# to the plot for each term
visualize(null_fits) + 
  shade_p_value(obs_fit, direction = "left")
}

# more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}

}
\seealso{
Other visualization functions: 
\code{\link{shade_confidence_interval}()}
}
\concept{visualization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{chisq_test}
\alias{chisq_test}
\title{Tidy chi-squared test}
\usage{
chisq_test(x, formula, response = NULL, explanatory = NULL, ...)
}
\arguments{
\item{x}{A data frame that can be coerced into a \link[tibble:tibble]{tibble}.}

\item{formula}{A formula with the response variable on the left and the
explanatory on the right.}

\item{response}{The variable name in \code{x} that will serve as the response.
This is alternative to using the \code{formula} argument.}

\item{explanatory}{The variable name in \code{x} that will serve as the
explanatory variable.}

\item{...}{Additional arguments for \link[stats:chisq.test]{chisq.test()}.}
}
\description{
A tidier version of \link[stats:chisq.test]{chisq.test()} for goodness of fit
tests and tests of independence.
}
\examples{
# chi-squared test of independence for college completion
# status depending on one's self-identified income class
chisq_test(gss, college ~ finrela)

# chi-squared goodness of fit test on whether self-identified
# income class follows a uniform distribution
chisq_test(gss,
           response = finrela,
           p = c("far below average" = 1/6,
                 "below average" = 1/6,
                 "average" = 1/6,
                 "above average" = 1/6,
                 "far above average" = 1/6,
                 "DK" = 1/6))

}
\seealso{
Other wrapper functions: 
\code{\link{chisq_stat}()},
\code{\link{observe}()},
\code{\link{prop_test}()},
\code{\link{t_stat}()},
\code{\link{t_test}()}
}
\concept{wrapper functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{t_test}
\alias{t_test}
\title{Tidy t-test}
\usage{
t_test(
  x,
  formula,
  response = NULL,
  explanatory = NULL,
  order = NULL,
  alternative = "two-sided",
  mu = 0,
  conf_int = TRUE,
  conf_level = 0.95,
  ...
)
}
\arguments{
\item{x}{A data frame that can be coerced into a \link[tibble:tibble]{tibble}.}

\item{formula}{A formula with the response variable on the left and the
explanatory on the right.}

\item{response}{The variable name in \code{x} that will serve as the response.
This is alternative to using the \code{formula} argument.}

\item{explanatory}{The variable name in \code{x} that will serve as the
explanatory variable.}

\item{order}{A string vector of specifying the order in which the levels of
the explanatory variable should be ordered for subtraction, where \code{order = c("first", "second")} means \code{("first" - "second")}.}

\item{alternative}{Character string giving the direction of the alternative
hypothesis. Options are \code{"two-sided"} (default), \code{"greater"}, or \code{"less"}.}

\item{mu}{A numeric value giving the hypothesized null mean value for a one
sample test and the hypothesized difference for a two sample test.}

\item{conf_int}{A logical value for whether to include the confidence
interval or not. \code{TRUE} by default.}

\item{conf_level}{A numeric value between 0 and 1. Default value is 0.95.}

\item{...}{For passing in other arguments to \link[stats:t.test]{t.test()}.}
}
\description{
A tidier version of \link[stats:t.test]{t.test()} for two sample tests.
}
\examples{
library(tidyr)

# t test for number of hours worked per week
# by college degree status
gss \%>\%
   tidyr::drop_na(college) \%>\%
   t_test(formula = hours ~ college,
      order = c("degree", "no degree"),
      alternative = "two-sided")

# see vignette("infer") for more explanation of the
# intuition behind the infer package, and vignette("t_test")
# for more examples of t-tests using infer

}
\seealso{
Other wrapper functions: 
\code{\link{chisq_stat}()},
\code{\link{chisq_test}()},
\code{\link{observe}()},
\code{\link{prop_test}()},
\code{\link{t_stat}()}
}
\concept{wrapper functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{deprecated}
\alias{deprecated}
\alias{conf_int}
\alias{p_value}
\title{Deprecated functions and objects}
\usage{
conf_int(x, level = 0.95, type = "percentile", point_estimate = NULL)

p_value(x, obs_stat, direction)
}
\arguments{
\item{x}{See the non-deprecated function.}

\item{level}{See the non-deprecated function.}

\item{type}{See the non-deprecated function.}

\item{point_estimate}{See the non-deprecated function.}

\item{obs_stat}{See the non-deprecated function.}

\item{direction}{See the non-deprecated function.}
}
\description{
These functions and objects should no longer be used. They will be removed
in a future release of \code{infer}.
}
\seealso{
\code{\link[=get_p_value]{get_p_value()}}, \code{\link[=get_confidence_interval]{get_confidence_interval()}}, \code{\link[=generate]{generate()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe}
\arguments{
\item{lhs, rhs}{Inference functions and the initial data frame.}
}
\description{
Like \{dplyr\}, \{infer\} also uses the pipe (\code{\%>\%}) function
from \code{magrittr} to turn function composition into a series of
iterative statements.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R, R/visualize.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{fit}
\alias{ggplot_add}
\title{Objects exported from other packages}
\details{
Read more about infer's \link[=fit.infer]{fit} function \link[=fit.infer]{here} or
by running \code{?fit.infer} in your console.
}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{generics}{\code{\link[generics]{fit}}}

  \item{ggplot2}{\code{\link[ggplot2]{ggplot_add}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/observe.R
\name{observe}
\alias{observe}
\title{Calculate observed statistics}
\usage{
observe(
  x,
  formula,
  response = NULL,
  explanatory = NULL,
  success = NULL,
  null = NULL,
  p = NULL,
  mu = NULL,
  med = NULL,
  sigma = NULL,
  stat = c("mean", "median", "sum", "sd", "prop", "count", "diff in means",
    "diff in medians", "diff in props", "Chisq", "F", "slope", "correlation", "t", "z",
    "ratio of props", "odds ratio"),
  order = NULL,
  ...
)
}
\arguments{
\item{x}{A data frame that can be coerced into a \link[tibble:tibble]{tibble}.}

\item{formula}{A formula with the response variable on the left and the
explanatory on the right. Alternatively, a \code{response} and \code{explanatory}
argument can be supplied.}

\item{response}{The variable name in \code{x} that will serve as the response.
This is an alternative to using the \code{formula} argument.}

\item{explanatory}{The variable name in \code{x} that will serve as the
explanatory variable. This is an alternative to using the formula argument.}

\item{success}{The level of \code{response} that will be considered a success, as
a string. Needed for inference on one proportion, a difference in
proportions, and corresponding z stats.}

\item{null}{The null hypothesis. Options include \code{"independence"} and
\code{"point"}.}

\item{p}{The true proportion of successes (a number between 0 and 1). To be used with point null hypotheses when the specified response
variable is categorical.}

\item{mu}{The true mean (any numerical value). To be used with point null
hypotheses when the specified response variable is continuous.}

\item{med}{The true median (any numerical value). To be used with point null
hypotheses when the specified response variable is continuous.}

\item{sigma}{The true standard deviation (any numerical value). To be used with
point null hypotheses.}

\item{stat}{A string giving the type of the statistic to calculate. Current
options include \code{"mean"}, \code{"median"}, \code{"sum"}, \code{"sd"}, \code{"prop"}, \code{"count"},
\code{"diff in means"}, \code{"diff in medians"}, \code{"diff in props"}, \code{"Chisq"} (or
\code{"chisq"}), \code{"F"} (or \code{"f"}), \code{"t"}, \code{"z"}, \code{"ratio of props"}, \code{"slope"},
\code{"odds ratio"}, and \code{"correlation"}. \code{infer} only supports theoretical
tests on one or two means via the \code{"t"} distribution and one or two
proportions via the \code{"z"}.}

\item{order}{A string vector of specifying the order in which the levels of
the explanatory variable should be ordered for subtraction (or division
for ratio-based statistics), where \code{order = c("first", "second")} means
\code{("first" - "second")}, or the analogue for ratios. Needed for inference on
difference in means, medians, proportions, ratios, t, and z statistics.}

\item{...}{To pass options like \code{na.rm = TRUE} into functions like
\link[base:mean]{mean()}, \link[stats:sd]{sd()}, etc. Can also be used to
supply hypothesized null values for the \code{"t"} statistic.}
}
\value{
A 1-column tibble containing the calculated statistic \code{stat}.
}
\description{
This function is a wrapper that calls \code{\link[=specify]{specify()}}, \code{\link[=hypothesize]{hypothesize()}}, and
\code{\link[=calculate]{calculate()}} consecutively that can be used to calculate observed
statistics from data. \code{\link[=hypothesize]{hypothesize()}} will only be called if a point
null hypothesis parameter is supplied.

Learn more in \code{vignette("infer")}.
}
\examples{
# calculating the observed mean number of hours worked per week
gss \%>\%
  observe(hours ~ NULL, stat = "mean")

# equivalently, calculating the same statistic with the core verbs
gss \%>\%
  specify(response = hours) \%>\%
  calculate(stat = "mean")

# calculating a t statistic for hypothesized mu = 40 hours worked/week
gss \%>\%
  observe(hours ~ NULL, stat = "t", null = "point", mu = 40)

# equivalently, calculating the same statistic with the core verbs
gss \%>\%
  specify(response = hours) \%>\%
  hypothesize(null = "point", mu = 40) \%>\%
  calculate(stat = "t")

# similarly for a difference in means in age based on whether
# the respondent has a college degree
observe(
  gss,
  age ~ college,
  stat = "diff in means",
  order = c("degree", "no degree")
)

# equivalently, calculating the same statistic with the core verbs
gss \%>\%
  specify(age ~ college) \%>\%
  calculate("diff in means", order = c("degree", "no degree"))
  
# for a more in-depth explanation of how to use the infer package
\dontrun{
vignette("infer")
}

}
\seealso{
Other wrapper functions: 
\code{\link{chisq_stat}()},
\code{\link{chisq_test}()},
\code{\link{prop_test}()},
\code{\link{t_stat}()},
\code{\link{t_test}()}

Other functions for calculating observed statistics: 
\code{\link{chisq_stat}()},
\code{\link{t_stat}()}
}
\concept{functions for calculating observed statistics}
\concept{wrapper functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rep_sample_n.R
\name{rep_sample_n}
\alias{rep_sample_n}
\alias{rep_slice_sample}
\title{Perform repeated sampling}
\usage{
rep_sample_n(tbl, size, replace = FALSE, reps = 1, prob = NULL)

rep_slice_sample(
  .data,
  n = NULL,
  prop = NULL,
  replace = FALSE,
  weight_by = NULL,
  reps = 1
)
}
\arguments{
\item{tbl, .data}{Data frame of population from which to sample.}

\item{size, n, prop}{\code{size} and \code{n} refer to the sample size of each sample.
The \code{size} argument to \code{rep_sample_n()} is required, while in
\code{rep_slice_sample()} sample size defaults to 1 if not specified. \code{prop}, an
argument to \code{rep_slice_sample()}, refers to the proportion of rows to sample
in each sample, and is rounded down in the case that \code{prop * nrow(.data)} is
not an integer. When using \code{rep_slice_sample()}, please only supply one of
\code{n} or \code{prop}.}

\item{replace}{Should samples be taken with replacement?}

\item{reps}{Number of samples to take.}

\item{prob, weight_by}{A vector of sampling weights for each of the rows in
\code{.data}—must have length equal to \code{nrow(.data)}.}
}
\value{
A tibble of size \code{reps * n} rows corresponding to \code{reps}
samples of size \code{n} from \code{.data}, grouped by \code{replicate}.
}
\description{
These functions extend the functionality of \code{\link[dplyr:sample_n]{dplyr::sample_n()}} and
\code{\link[dplyr:slice]{dplyr::slice_sample()}} by allowing for repeated sampling of data.
This operation is especially helpful while creating sampling
distributions—see the examples below!
}
\details{
\code{rep_sample_n()} and \code{rep_slice_sample()} are designed to behave similar to
their dplyr counterparts. As such, they have at least the following
differences:
\itemize{
\item In case \code{replace = FALSE} having \code{size} bigger than number of data rows in
\code{rep_sample_n()} will give an error. In \code{rep_slice_sample()} having such \code{n}
or \code{prop > 1} will give warning and output sample size will be set to number
of rows in data.
}

Note that the \code{\link[dplyr:sample_n]{dplyr::sample_n()}} function  has been superseded by
\code{\link[dplyr:slice]{dplyr::slice_sample()}}.
}
\examples{
library(dplyr)
library(ggplot2)
library(tibble)

# take 1000 samples of size n = 50, without replacement
slices <- gss \%>\%
  rep_slice_sample(n = 50, reps = 1000)

slices

# compute the proportion of respondents with a college
# degree in each replicate
p_hats <- slices \%>\%
  group_by(replicate) \%>\%
  summarize(prop_college = mean(college == "degree"))

# plot sampling distribution
ggplot(p_hats, aes(x = prop_college)) +
  geom_density() +
  labs(
    x = "p_hat", y = "Number of samples",
    title = "Sampling distribution of p_hat"
  )

# sampling with probability weights. Note probabilities are automatically
# renormalized to sum to 1
df <- tibble(
  id = 1:5,
  letter = factor(c("a", "b", "c", "d", "e"))
)

rep_slice_sample(df, n = 2, reps = 5, weight_by = c(.5, .4, .3, .2, .1))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gss.R
\docType{data}
\name{gss}
\alias{gss}
\title{Subset of data from the General Social Survey (GSS).}
\format{
A tibble with 500 rows and 11 variables:
\describe{
\item{year}{year respondent was surveyed}
\item{age}{age at time of survey, truncated at 89}
\item{sex}{respondent's sex (self-identified)}
\item{college}{whether on not respondent has a college degree, including
junior/community college}
\item{partyid}{political party affiliation}
\item{hompop}{number of persons in household}
\item{hours}{number of hours worked in week before survey, truncated at 89}
\item{income}{total family income}
\item{class}{subjective socioeconomic class identification}
\item{finrela}{opinion of family income}
\item{weight}{survey weight}
}
}
\source{
\url{https://gss.norc.org}
}
\usage{
gss
}
\description{
The General Social Survey is a high-quality survey which gathers data on
American society and opinions, conducted since 1972. This data set is a
sample of 500 entries from the GSS, spanning years 1973-2018,
including demographic markers and some
economic variables. Note that this data is included for demonstration only,
and should not be assumed to provide accurate estimates relating to the GSS.
However, due to the high quality of the GSS, the unweighted data will
approximate the weighted data in some analyses.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{t_stat}
\alias{t_stat}
\title{Tidy t-test statistic}
\usage{
t_stat(
  x,
  formula,
  response = NULL,
  explanatory = NULL,
  order = NULL,
  alternative = "two-sided",
  mu = 0,
  conf_int = FALSE,
  conf_level = 0.95,
  ...
)
}
\arguments{
\item{x}{A data frame that can be coerced into a \link[tibble:tibble]{tibble}.}

\item{formula}{A formula with the response variable on the left and the
explanatory on the right.}

\item{response}{The variable name in \code{x} that will serve as the response.
This is alternative to using the \code{formula} argument.}

\item{explanatory}{The variable name in \code{x} that will serve as the
explanatory variable.}

\item{order}{A string vector of specifying the order in which the levels of
the explanatory variable should be ordered for subtraction, where \code{order = c("first", "second")} means \code{("first" - "second")}.}

\item{alternative}{Character string giving the direction of the alternative
hypothesis. Options are \code{"two-sided"} (default), \code{"greater"}, or \code{"less"}.}

\item{mu}{A numeric value giving the hypothesized null mean value for a one
sample test and the hypothesized difference for a two sample test.}

\item{conf_int}{A logical value for whether to include the confidence
interval or not. \code{TRUE} by default.}

\item{conf_level}{A numeric value between 0 and 1. Default value is 0.95.}

\item{...}{Pass in arguments to \\{infer\\} functions.}
}
\description{
A shortcut wrapper function to get the observed test statistic for a t test.
This function has been deprecated in favor of the more general \code{\link[=observe]{observe()}}.
}
\examples{
library(tidyr)

# t test statistic for true mean number of hours worked
# per week of 40
gss \%>\%
   t_stat(response = hours, mu = 40)

# t test statistic for number of hours worked per week
# by college degree status
gss \%>\%
   tidyr::drop_na(college) \%>\%
   t_stat(formula = hours ~ college,
      order = c("degree", "no degree"),
      alternative = "two-sided")

}
\seealso{
Other wrapper functions: 
\code{\link{chisq_stat}()},
\code{\link{chisq_test}()},
\code{\link{observe}()},
\code{\link{prop_test}()},
\code{\link{t_test}()}

Other functions for calculating observed statistics: 
\code{\link{chisq_stat}()},
\code{\link{observe}()}
}
\concept{functions for calculating observed statistics}
\concept{wrapper functions}
