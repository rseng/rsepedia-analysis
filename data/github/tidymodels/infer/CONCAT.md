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
