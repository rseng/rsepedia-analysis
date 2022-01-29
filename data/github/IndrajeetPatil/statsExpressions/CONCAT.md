
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `{statsExpressions}`: Tidy dataframes and expressions with statistical details

| Status                                                                                                                                                      | Usage                                                                                                                                                      | Miscellaneous                                                                                                                                                              |
|-------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [![R build status](https://github.com/IndrajeetPatil/statsExpressions/workflows/R-CMD-check/badge.svg)](https://github.com/IndrajeetPatil/statsExpressions) | [![Total downloads badge](https://cranlogs.r-pkg.org/badges/grand-total/statsExpressions?color=blue)](https://CRAN.R-project.org/package=statsExpressions) | [![Codecov](https://codecov.io/gh/IndrajeetPatil/statsExpressions/branch/master/graph/badge.svg)](https://app.codecov.io/gh/IndrajeetPatil/statsExpressions?branch=master) |
| [![lints](https://github.com/IndrajeetPatil/statsExpressions/workflows/lint/badge.svg)](https://github.com/IndrajeetPatil/statsExpressions)                 | [![Daily downloads badge](https://cranlogs.r-pkg.org/badges/last-day/statsExpressions?color=blue)](https://CRAN.R-project.org/package=statsExpressions)    | [![status](https://tinyverse.netlify.com/badge/statsExpressions)](https://CRAN.R-project.org/package=statsExpressions)                                                     |
| [![pkgdown](https://github.com/IndrajeetPatil/statsExpressions/workflows/pkgdown/badge.svg)](https://github.com/IndrajeetPatil/statsExpressions/actions)    | [![DOI](https://joss.theoj.org/papers/10.21105/joss.03236/status.svg)](https://doi.org/10.21105/joss.03236)                                                | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)                                                 |

# Introduction <img src="man/figures/logo.png" align="right" width="240" />

The `{statsExpressions}` package has two key aims:

-   to provide a consistent syntax to do statistical analysis with tidy
    data (in pipe-friendly manner),
-   to provide statistical expressions (pre-formatted in-text
    statistical results) for plotting functions.

Statistical packages exhibit substantial diversity in terms of their
syntax and expected input type. This can make it difficult to switch
from one statistical approach to another. For example, some functions
expect vectors as inputs, while others expect dataframes. Depending on
whether it is a repeated measures design or not, different functions
might expect data to be in wide or long format. Some functions can
internally omit missing values, while other functions error in their
presence. Furthermore, if someone wishes to utilize the objects returned
by these packages downstream in their workflow, this is not
straightforward either because even functions from the same package can
return a list, a matrix, an array, a dataframe, etc., depending on the
function.

This is where `{statsExpressions}` comes in: It can be thought of as a
unified portal through which most of the functionality in these
underlying packages can be accessed, with a simpler interface and no
requirement to change data format.

# Installation

| Type        | Source                                                                                                                       | Command                                                      |
|-------------|------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------|
| Release     | [![CRAN Status](https://www.r-pkg.org/badges/version/statsExpressions)](https://cran.r-project.org/package=statsExpressions) | `install.packages("statsExpressions")`                       |
| Development | [![Project Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/##active)                | `remotes::install_github("IndrajeetPatil/statsExpressions")` |

# Citation

The package can be cited as:

``` r
citation("statsExpressions")

  Patil, I., (2021). statsExpressions: R Package for Tidy Dataframes
  and Expressions with Statistical Details. Journal of Open Source
  Software, 6(61), 3236, https://doi.org/10.21105/joss.03236

A BibTeX entry for LaTeX users is

  @Article{,
    doi = {10.21105/joss.03236},
    url = {https://doi.org/10.21105/joss.03236},
    year = {2021},
    publisher = {{The Open Journal}},
    volume = {6},
    number = {61},
    pages = {3236},
    author = {Indrajeet Patil},
    title = {{statsExpressions: {R} Package for Tidy Dataframes and Expressions with Statistical Details}},
    journal = {{Journal of Open Source Software}},
  }
```

# General Workflow

<img src="man/figures/card.png" width="80%" />

# Summary of types of statistical analyses

Here is a tabular summary of available tests:

| Test                       | Function            | Lifecycle                                                                                                                       |
|----------------------------|---------------------|---------------------------------------------------------------------------------------------------------------------------------|
| one-sample *t*-test        | `one_sample_test`   | [![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| two-sample *t*-test        | `two_sample_test`   | [![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| one-way ANOVA              | `oneway_anova`      | [![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| correlation analysis       | `corr_test`         | [![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| contingency table analysis | `contingency_table` | [![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| meta-analysis              | `meta_analysis`     | [![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html) |

The table below summarizes all the different types of analyses currently
supported in this package-

| Description                                       | Parametric | Non-parametric | Robust | Bayesian |
|---------------------------------------------------|------------|----------------|--------|----------|
| Between group/condition comparisons               | ✅         | ✅             | ✅     | ✅       |
| Within group/condition comparisons                | ✅         | ✅             | ✅     | ✅       |
| Distribution of a numeric variable                | ✅         | ✅             | ✅     | ✅       |
| Correlation between two variables                 | ✅         | ✅             | ✅     | ✅       |
| Association between categorical variables         | ✅         | ✅             | ❌     | ✅       |
| Equal proportions for categorical variable levels | ✅         | ✅             | ❌     | ✅       |
| Random-effects meta-analysis                      | ✅         | ❌             | ✅     | ✅       |

Summary of Bayesian analysis

| Analysis                        | Hypothesis testing | Estimation |
|---------------------------------|--------------------|------------|
| (one/two-sample) *t*-test       | ✅                 | ✅         |
| one-way ANOVA                   | ✅                 | ✅         |
| correlation                     | ✅                 | ✅         |
| (one/two-way) contingency table | ✅                 | ✅         |
| random-effects meta-analysis    | ✅                 | ✅         |

# Tidy dataframes from statistical analysis

To illustrate the simplicity of this syntax, let’s say we want to run a
one-way ANOVA. If we first run a non-parametric ANOVA and then decide to
run a robust ANOVA instead, the syntax remains the same and the
statistical approach can be modified by changing a single argument:

``` r
library(statsExpressions)

mtcars %>% oneway_anova(cyl, wt, type = "nonparametric")
#> # A tibble: 1 x 15
#>   parameter1 parameter2 statistic df.error   p.value
#>   <chr>      <chr>          <dbl>    <int>     <dbl>
#> 1 wt         cyl             22.8        2 0.0000112
#>   method                       effectsize      estimate conf.level conf.low
#>   <chr>                        <chr>              <dbl>      <dbl>    <dbl>
#> 1 Kruskal-Wallis rank sum test Epsilon2 (rank)    0.736       0.95    0.624
#>   conf.high conf.method          conf.iterations n.obs expression  
#>       <dbl> <chr>                          <int> <int> <list>      
#> 1         1 percentile bootstrap             100    32 <expression>

mtcars %>% oneway_anova(cyl, wt, type = "robust")
#> # A tibble: 1 x 12
#>   statistic    df df.error p.value
#>       <dbl> <dbl>    <dbl>   <dbl>
#> 1      12.7     2     12.2 0.00102
#>   method                                           
#>   <chr>                                            
#> 1 A heteroscedastic one-way ANOVA for trimmed means
#>   effectsize                         estimate conf.level conf.low conf.high
#>   <chr>                                 <dbl>      <dbl>    <dbl>     <dbl>
#> 1 Explanatory measure of effect size     1.05       0.95    0.843      1.50
#>   n.obs expression  
#>   <int> <list>      
#> 1    32 <expression>
```

All possible output dataframes from functions are tabulated here:
<https://indrajeetpatil.github.io/statsExpressions/articles/web_only/dataframe_outputs.html>

Needless to say this will also work with the `kable` function to
generate a table:

``` r
# setup
library(statsExpressions)
set.seed(123)

# one-sample robust t-test
# we will leave `expression` column out; it's not needed for using only the dataframe
mtcars %>%
  one_sample_test(wt, test.value = 3, type = "robust") %>%
  dplyr::select(-expression) %>%
  knitr::kable()
```

| statistic | p.value | n.obs | method                                 | effectsize   | estimate | conf.level | conf.low | conf.high |
|----------:|--------:|------:|:---------------------------------------|:-------------|---------:|-----------:|---------:|----------:|
|  1.179181 |   0.275 |    32 | Bootstrap-t method for one-sample test | Trimmed mean |    3.197 |       0.95 | 2.854246 |  3.539754 |

These functions are also compatible with other popular data manipulation
packages.

For example, let’s say we want to run a one-sample *t*-test for all
levels of a certain grouping variable. We can use `dplyr` to do so:

``` r
# for reproducibility
set.seed(123)
library(dplyr)

# grouped operation
# running one-sample test for all levels of grouping variable `cyl`
mtcars %>%
  group_by(cyl) %>%
  group_modify(~ one_sample_test(.x, wt, test.value = 3), .keep = TRUE) %>%
  ungroup()
#> # A tibble: 3 x 16
#>     cyl    mu statistic df.error  p.value method            alternative
#>   <dbl> <dbl>     <dbl>    <dbl>    <dbl> <chr>             <chr>      
#> 1     4     3    -4.16        10 0.00195  One Sample t-test two.sided  
#> 2     6     3     0.870        6 0.418    One Sample t-test two.sided  
#> 3     8     3     4.92        13 0.000278 One Sample t-test two.sided  
#>   effectsize estimate conf.level conf.low conf.high conf.method
#>   <chr>         <dbl>      <dbl>    <dbl>     <dbl> <chr>      
#> 1 Hedges' g    -1.16        0.95   -1.97     -0.422 ncp        
#> 2 Hedges' g     0.286       0.95   -0.419     1.01  ncp        
#> 3 Hedges' g     1.24        0.95    0.565     1.98  ncp        
#>   conf.distribution n.obs expression  
#>   <chr>             <int> <list>      
#> 1 t                    11 <expression>
#> 2 t                     7 <expression>
#> 3 t                    14 <expression>
```

# Using expressions in custom plots

Note that *expression* here means **a pre-formatted in-text statistical
result**. In addition to other details contained in the dataframe, there
is also a column titled `expression`, which contains expression with
statistical details and can be displayed in a plot.

For **all** statistical test expressions, the default template attempt
to follow the gold standard for statistical reporting.

For example, here are results from Welch’s *t*-test:

<img src="man/figures/stats_reporting_format.png" align="center" />

## Expressions for centrality measure

``` r
library(ggplot2)

# displaying mean for each level of `cyl`
centrality_description(mtcars, cyl, wt) |>
  ggplot(aes(cyl, wt)) +
  geom_point() +
  geom_label(aes(label = expression), parse = TRUE)
```

<img src="man/figures/README-centrality-1.png" width="100%" />

Here are a few examples for supported analyses.

## Expressions for one-way ANOVAs

### Between-subjects design

Let’s say we want to check differences in weight of the vehicle based on
number of cylinders in the engine and wish to carry out robust
trimmed-means ANOVA:

``` r
# setup
set.seed(123)
library(ggplot2)
library(statsExpressions)
library(ggridges)

# create a ridgeplot
ggplot(iris, aes(x = Sepal.Length, y = Species)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE,
    scale = 0.9, vline_size = 1, vline_color = "red",
    position = position_raincloud(adjust_vlines = TRUE)
  ) + # use the expression in the dataframe to display results in the subtitle
  labs(
    title = "A heteroscedastic one-way ANOVA for trimmed means",
    subtitle = oneway_anova(iris, Species, Sepal.Length, type = "robust")$expression[[1]]
  )
```

<img src="man/figures/README-anova_rob1-1.png" width="100%" />

### Within-subjects design

Let’s now see an example of a repeated measures one-way ANOVA.

``` r
# setup
set.seed(123)
library(ggplot2)
library(WRS2)
library(ggbeeswarm)
library(statsExpressions)

ggplot2::ggplot(WineTasting, aes(Wine, Taste, color = Wine)) +
  geom_quasirandom() +
  labs(
    title = "Friedman's rank sum test",
    subtitle = oneway_anova(
      WineTasting,
      Wine,
      Taste,
      paired = TRUE,
      subject.id = Taster,
      type = "np"
    )$expression[[1]]
  )
```

<img src="man/figures/README-anova_parametric2-1.png" width="100%" />

## Expressions for two-sample tests

### Between-subjects design

``` r
# setup
set.seed(123)
library(ggplot2)
library(gghalves)
library(ggbeeswarm)
library(hrbrthemes)

# create a plot
ggplot(ToothGrowth, aes(supp, len)) +
  geom_half_boxplot() +
  geom_beeswarm() +
  theme_ipsum_rc() +
  # adding a subtitle with
  labs(
    title = "Two-Sample Welch's t-test",
    subtitle = two_sample_test(ToothGrowth, supp, len)$expression[[1]]
  )
```

<img src="man/figures/README-t_two-1.png" width="100%" />

### Within-subjects design

We can also have a look at a repeated measures design and the related
expressions.

``` r
# setup
set.seed(123)
library(ggplot2)
library(tidyr)
library(PairedData)
data(PrisonStress)

# get data in tidy format
df <- pivot_longer(PrisonStress, starts_with("PSS"), "PSS", values_to = "stress")

# plot
paired.plotProfiles(PrisonStress, "PSSbefore", "PSSafter", subjects = "Subject") +
  labs(
    title = "Two-sample Wilcoxon paired test",
    subtitle = two_sample_test(
      data = df,
      x = PSS,
      y = stress,
      paired = TRUE,
      subject.id = Subject,
      type = "np"
    )$expression[[1]]
  )
```

<img src="man/figures/README-t_two_paired1-1.png" width="100%" />

## Expressions for one-sample tests

``` r
# setup
set.seed(123)
library(ggplot2)

# dataframe with results
df_results <- one_sample_test(mtcars, wt, test.value = 3, type = "bayes",
                              top.text = "Bayesian one-sample t-test")

# creating a histogram plot
ggplot(mtcars, aes(wt)) +
  geom_histogram(alpha = 0.5) +
  geom_vline(xintercept = mean(mtcars$wt), color = "red") +
  labs(
    subtitle = df_results$expression[[1]]
  )
```

<img src="man/figures/README-t_one-1.png" width="100%" />

## Expressions for correlation analysis

Let’s look at another example where we want to run correlation analysis:

``` r
# setup
set.seed(123)
library(ggplot2)

# create a scatter plot
ggplot(mtcars, aes(mpg, wt)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  labs(
    title = "Spearman's rank correlation coefficient",
    subtitle = corr_test(mtcars, mpg, wt, type = "nonparametric")$expression[[1]]
  )
```

<img src="man/figures/README-corr-1.png" width="100%" />

## Expressions for contingency table analysis

For categorical/nominal data - one-sample:

``` r
# setup
set.seed(123)
library(ggplot2)

df_results <- contingency_table(as.data.frame(table(mpg$class)),
  Var1,
  counts = Freq,
  type = "bayes",
  top.text = "One-sample goodness-of-fit test"
)

# basic pie chart
ggplot(as.data.frame(table(mpg$class)), aes(x = "", y = Freq, fill = factor(Var1))) +
  geom_bar(width = 1, stat = "identity") +
  theme(axis.line = element_blank()) +
  # cleaning up the chart and adding results from one-sample proportion test
  coord_polar(theta = "y", start = 0) +
  labs(
    fill = "Class",
    x = NULL,
    y = NULL,
    title = "Pie Chart of class (type of car)",
    caption = df_results$expression[[1]]
  )
```

<img src="man/figures/README-gof-1.png" width="100%" />

You can also use these function to get the expression in return without
having to display them in plots:

``` r
# setup
set.seed(123)
library(ggplot2)

# Pearson's chi-squared test of independence
contingency_table(mtcars, am, cyl)$expression[[1]]
#> expression(list(chi["Pearson"]^2 * "(" * 2 * ")" == "8.74", italic(p) == 
#>     "0.01", widehat(italic("V"))["Cramer"] == "0.46", CI["95%"] ~ 
#>     "[" * "0.00", "1.00" * "]", italic("n")["obs"] == "32"))
```

## Expressions for meta-analysis

``` r
# setup
set.seed(123)
library(metaviz)
library(ggplot2)
library(metaplus)

# meta-analysis forest plot with results random-effects meta-analysis
viz_forest(
  x = mozart[, c("d", "se")],
  study_labels = mozart[, "study_name"],
  xlab = "Cohen's d",
  variant = "thick",
  type = "cumulative"
) + # use `{statsExpressions}` to create expression containing results
  labs(
    title = "Meta-analysis of Pietschnig, Voracek, and Formann (2010) on the Mozart effect",
    subtitle = meta_analysis(dplyr::rename(mozart, estimate = d, std.error = se))$expression[[1]]
  ) +
  theme(text = element_text(size = 12))
```

<img src="man/figures/README-metaanalysis-1.png" width="100%" />

# Customizing details to your liking

Sometimes you may not wish include so many details in the subtitle. In
that case, you can extract the expression and copy-paste only the part
you wish to include. For example, here only statistic and *p*-values are
included:

``` r
# setup
set.seed(123)
library(ggplot2)

# extracting detailed expression
(res_expr <- oneway_anova(iris, Species, Sepal.Length, var.equal = TRUE)$expression[[1]])
#> expression(list(italic("F")["Fisher"](2, 147) == "119.26", italic(p) == 
#>     "1.67e-31", widehat(omega["p"]^2) == "0.61", CI["95%"] ~ 
#>     "[" * "0.53", "1.00" * "]", italic("n")["obs"] == "150"))

# adapting the details to your liking
ggplot(iris, aes(x = Species, y = Sepal.Length)) +
  geom_boxplot() +
  labs(subtitle = ggplot2::expr(paste(
    NULL, italic("F"), "(", "2",
    ",", "147", ") = ", "119.26", ", ",
    italic("p"), " = ", "1.67e-31"
  )))
```

<img src="man/figures/README-custom_expr-1.png" width="100%" />

# Summary of tests and effect sizes

Here a go-to summary about statistical test carried out and the returned
effect size for each function is provided. This should be useful if one
needs to find out more information about how an argument is resolved in
the underlying package or if one wishes to browse the source code. So,
for example, if you want to know more about how one-way
(between-subjects) ANOVA, you can run `?stats::oneway.test` in your R
console.

## `centrality_description`

| Type           | Measure                                           | Function used                       |
|----------------|---------------------------------------------------|-------------------------------------|
| Parametric     | mean                                              | `parameters::describe_distribution` |
| Non-parametric | median                                            | `parameters::describe_distribution` |
| Robust         | trimmed mean                                      | `parameters::describe_distribution` |
| Bayesian       | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution` |

## `two_sample_test` + `oneway_anova`

No. of groups: `2` =\> `two_sample_test`<br> No. of groups: `> 2` =\>
`oneway_anova`

### between-subjects

**Hypothesis testing**

| Type           | No. of groups | Test                                            | Function used          |
|----------------|---------------|-------------------------------------------------|------------------------|
| Parametric     | \> 2          | Fisher’s or Welch’s one-way ANOVA               | `stats::oneway.test`   |
| Non-parametric | \> 2          | Kruskal–Wallis one-way ANOVA                    | `stats::kruskal.test`  |
| Robust         | \> 2          | Heteroscedastic one-way ANOVA for trimmed means | `WRS2::t1way`          |
| Bayes Factor   | \> 2          | Fisher’s ANOVA                                  | `BayesFactor::anovaBF` |
| Parametric     | 2             | Student’s or Welch’s *t*-test                   | `stats::t.test`        |
| Non-parametric | 2             | Mann–Whitney *U* test                           | `stats::wilcox.test`   |
| Robust         | 2             | Yuen’s test for trimmed means                   | `WRS2::yuen`           |
| Bayesian       | 2             | Student’s *t*-test                              | `BayesFactor::ttestBF` |

**Effect size estimation**

| Type           | No. of groups | Effect size                                                                                                                                                                                          | CI? | Function used                                          |
|----------------|---------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|--------------------------------------------------------|
| Parametric     | \> 2          | ![\\eta\_{p}^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Ceta_%7Bp%7D%5E2 "\eta_{p}^2"), ![\\omega\_{p}^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Comega_%7Bp%7D%5E2 "\omega_{p}^2") | ✅  | `effectsize::omega_squared`, `effectsize::eta_squared` |
| Non-parametric | \> 2          | ![\\epsilon\_{ordinal}^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Cepsilon_%7Bordinal%7D%5E2 "\epsilon_{ordinal}^2")                                                                         | ✅  | `effectsize::rank_epsilon_squared`                     |
| Robust         | \> 2          | ![\\xi](http://chart.apis.google.com/chart?cht=tx&chl=%5Cxi "\xi") (Explanatory measure of effect size)                                                                                              | ✅  | `WRS2::t1way`                                          |
| Bayes Factor   | \> 2          | ![R\_{Bayesian}^2](http://chart.apis.google.com/chart?cht=tx&chl=R_%7BBayesian%7D%5E2 "R_{Bayesian}^2")                                                                                              | ✅  | `performance::r2_bayes`                                |
| Parametric     | 2             | Cohen’s *d*, Hedge’s *g*                                                                                                                                                                             | ✅  | `effectsize::cohens_d`, `effectsize::hedges_g`         |
| Non-parametric | 2             | *r* (rank-biserial correlation)                                                                                                                                                                      | ✅  | `effectsize::rank_biserial`                            |
| Robust         | 2             | ![\\delta\_{R}^{AKP}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7BR%7D%5E%7BAKP%7D "\delta_{R}^{AKP}") (Algina-Keselman-Penfield robust standardized difference)                       | ✅  | `WRS2::akp.effect`                                     |
| Bayesian       | 2             | ![\\delta\_{posterior}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7Bposterior%7D "\delta_{posterior}")                                                                                 | ✅  | `bayestestR::describe_posterior`                       |

### within-subjects

**Hypothesis testing**

| Type           | No. of groups | Test                                                              | Function used          |
|----------------|---------------|-------------------------------------------------------------------|------------------------|
| Parametric     | \> 2          | One-way repeated measures ANOVA                                   | `afex::aov_ez`         |
| Non-parametric | \> 2          | Friedman rank sum test                                            | `stats::friedman.test` |
| Robust         | \> 2          | Heteroscedastic one-way repeated measures ANOVA for trimmed means | `WRS2::rmanova`        |
| Bayes Factor   | \> 2          | One-way repeated measures ANOVA                                   | `BayesFactor::anovaBF` |
| Parametric     | 2             | Student’s *t*-test                                                | `stats::t.test`        |
| Non-parametric | 2             | Wilcoxon signed-rank test                                         | `stats::wilcox.test`   |
| Robust         | 2             | Yuen’s test on trimmed means for dependent samples                | `WRS2::yuend`          |
| Bayesian       | 2             | Student’s *t*-test                                                | `BayesFactor::ttestBF` |

**Effect size estimation**

| Type           | No. of groups | Effect size                                                                                                                                                                                          | CI? | Function used                                          |
|----------------|---------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|--------------------------------------------------------|
| Parametric     | \> 2          | ![\\eta\_{p}^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Ceta_%7Bp%7D%5E2 "\eta_{p}^2"), ![\\omega\_{p}^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Comega_%7Bp%7D%5E2 "\omega_{p}^2") | ✅  | `effectsize::omega_squared`, `effectsize::eta_squared` |
| Non-parametric | \> 2          | ![W\_{Kendall}](http://chart.apis.google.com/chart?cht=tx&chl=W_%7BKendall%7D "W_{Kendall}") (Kendall’s coefficient of concordance)                                                                  | ✅  | `effectsize::kendalls_w`                               |
| Robust         | \> 2          | ![\\delta\_{R-avg}^{AKP}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7BR-avg%7D%5E%7BAKP%7D "\delta_{R-avg}^{AKP}") (Algina-Keselman-Penfield robust standardized difference average)   | ✅  | `WRS2::wmcpAKP`                                        |
| Bayes Factor   | \> 2          | ![R\_{Bayesian}^2](http://chart.apis.google.com/chart?cht=tx&chl=R_%7BBayesian%7D%5E2 "R_{Bayesian}^2")                                                                                              | ✅  | `performance::r2_bayes`                                |
| Parametric     | 2             | Cohen’s *d*, Hedge’s *g*                                                                                                                                                                             | ✅  | `effectsize::cohens_d`, `effectsize::hedges_g`         |
| Non-parametric | 2             | *r* (rank-biserial correlation)                                                                                                                                                                      | ✅  | `effectsize::rank_biserial`                            |
| Robust         | 2             | ![\\delta\_{R}^{AKP}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7BR%7D%5E%7BAKP%7D "\delta_{R}^{AKP}") (Algina-Keselman-Penfield robust standardized difference)                       | ✅  | `WRS2::wmcpAKP`                                        |
| Bayesian       | 2             | ![\\delta\_{posterior}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7Bposterior%7D "\delta_{posterior}")                                                                                 | ✅  | `bayestestR::describe_posterior`                       |

## `one_sample_test`

**Hypothesis testing**

| Type           | Test                                     | Function used          |
|----------------|------------------------------------------|------------------------|
| Parametric     | One-sample Student’s *t*-test            | `stats::t.test`        |
| Non-parametric | One-sample Wilcoxon test                 | `stats::wilcox.test`   |
| Robust         | Bootstrap-*t* method for one-sample test | `WRS2::trimcibt`       |
| Bayesian       | One-sample Student’s *t*-test            | `BayesFactor::ttestBF` |

**Effect size estimation**

| Type           | Effect size                                                                                                          | CI? | Function used                                  |
|----------------|----------------------------------------------------------------------------------------------------------------------|-----|------------------------------------------------|
| Parametric     | Cohen’s *d*, Hedge’s *g*                                                                                             | ✅  | `effectsize::cohens_d`, `effectsize::hedges_g` |
| Non-parametric | *r* (rank-biserial correlation)                                                                                      | ✅  | `effectsize::rank_biserial`                    |
| Robust         | trimmed mean                                                                                                         | ✅  | `trimcibt` (custom)                            |
| Bayes Factor   | ![\\delta\_{posterior}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7Bposterior%7D "\delta_{posterior}") | ✅  | `bayestestR::describe_posterior`               |

## `corr_test`

**Hypothesis testing** and **Effect size estimation**

| Type           | Test                                       | CI? | Function used              |
|----------------|--------------------------------------------|-----|----------------------------|
| Parametric     | Pearson’s correlation coefficient          | ✅  | `correlation::correlation` |
| Non-parametric | Spearman’s rank correlation coefficient    | ✅  | `correlation::correlation` |
| Robust         | Winsorized Pearson correlation coefficient | ✅  | `correlation::correlation` |
| Bayesian       | Pearson’s correlation coefficient          | ✅  | `correlation::correlation` |

## `contingency_table`

### two-way table

**Hypothesis testing**

| Type                      | Design   | Test                                                                                                  | Function used                     |
|---------------------------|----------|-------------------------------------------------------------------------------------------------------|-----------------------------------|
| Parametric/Non-parametric | Unpaired | Pearson’s ![\\chi^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test          | `stats::chisq.test`               |
| Bayesian                  | Unpaired | Bayesian Pearson’s ![\\chi^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test | `BayesFactor::contingencyTableBF` |
| Parametric/Non-parametric | Paired   | McNemar’s ![\\chi^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test          | `stats::mcnemar.test`             |
| Bayesian                  | Paired   | ❌                                                                                                    | ❌                                |

**Effect size estimation**

| Type                      | Design   | Effect size                                                        | CI? | Function used           |
|---------------------------|----------|--------------------------------------------------------------------|-----|-------------------------|
| Parametric/Non-parametric | Unpaired | Cramer’s ![V](http://chart.apis.google.com/chart?cht=tx&chl=V "V") | ✅  | `effectsize::cramers_v` |
| Bayesian                  | Unpaired | Cramer’s ![V](http://chart.apis.google.com/chart?cht=tx&chl=V "V") | ✅  | `effectsize::cramers_v` |
| Parametric/Non-parametric | Paired   | Cohen’s ![g](http://chart.apis.google.com/chart?cht=tx&chl=g "g")  | ✅  | `effectsize::cohens_g`  |
| Bayesian                  | Paired   | ❌                                                                 | ❌  | ❌                      |

### one-way table

**Hypothesis testing**

| Type                      | Test                                                                                                        | Function used       |
|---------------------------|-------------------------------------------------------------------------------------------------------------|---------------------|
| Parametric/Non-parametric | Goodness of fit ![\\chi^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test          | `stats::chisq.test` |
| Bayesian                  | Bayesian Goodness of fit ![\\chi^2](http://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test | (custom)            |

**Effect size estimation**

| Type                      | Effect size                                                         | CI? | Function used            |
|---------------------------|---------------------------------------------------------------------|-----|--------------------------|
| Parametric/Non-parametric | Pearson’s ![C](http://chart.apis.google.com/chart?cht=tx&chl=C "C") | ✅  | `effectsize::pearsons_c` |
| Bayesian                  | ❌                                                                  | ❌  | ❌                       |

## `meta_analysis`

**Hypothesis testing** and **Effect size estimation**

| Type       | Test                                             | Effect size                                                              | CI? | Function used          |
|------------|--------------------------------------------------|--------------------------------------------------------------------------|-----|------------------------|
| Parametric | Meta-analysis via random-effects models          | ![\\beta](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbeta "\beta") | ✅  | `metafor::metafor`     |
| Robust     | Meta-analysis via robust random-effects models   | ![\\beta](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbeta "\beta") | ✅  | `metaplus::metaplus`   |
| Bayes      | Meta-analysis via Bayesian random-effects models | ![\\beta](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbeta "\beta") | ✅  | `metaBMA::meta_random` |

# Usage in `ggstatsplot`

Note that these functions were initially written to display results from
statistical tests on ready-made `ggplot2` plots implemented in
`ggstatsplot`.

For detailed documentation, see the package website:
<https://indrajeetpatil.github.io/ggstatsplot/>

Here is an example from `ggstatsplot` of what the plots look like when
the expressions are displayed in the subtitle-

<img src="man/figures/ggstatsplot.png" align="center" />

# Acknowledgments

The hexsticker and the schematic illustration of general workflow were
generously designed by Sarah Otterstetter (Max Planck Institute for
Human Development, Berlin).

# Contributing

I’m happy to receive bug reports, suggestions, questions, and (most of
all) contributions to fix problems and add features. I personally prefer
using the `GitHub` issues system over trying to reach out to me in other
ways (personal e-mail, Twitter, etc.). Pull Requests for contributions
are encouraged.

Here are some simple ways in which you can contribute (in the increasing
order of commitment):

-   Read and correct any inconsistencies in the
    [documentation](https://indrajeetpatil.github.io/statsExpressions/)

-   Raise issues about bugs or wanted features

-   Review code

-   Add new functionality (in the form of new plotting functions or
    helpers for preparing subtitles)

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/IndrajeetPatil/statsExpressions/blob/master/.github/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.
# statsExpressions 1.3.0.9000

# statsExpressions 1.3.0

BREAKING CHANGES

  - The `format_num()` has been removed in favor of `insight::format_value()`.

MINOR CHANGES

  - The `expr_template()` has been renamed to more informative
    `add_expression_col()` function and has a different API. It returns a
    dataframe with the additional expression column instead of just the
    expression.

# statsExpressions 1.2.0

BREAKING CHANGES

  - A number of effect size estimates and their confidence intervals have
    changed due to respective changes made in `{effectsize}` package version
    `0.5` release. For full details of these changes, see:
    <https://easystats.github.io/effectsize/news/index.html>

  - For the same reason, the effect size for one-way contingency table has
    changed from Cramer's *V* to Pearson's *C*.

NEW FUNCTIONS

  - `centrality_description()` function added to describe distribution for each
    level of a grouping variable and create an expression describing a
    centrality measure.

  - Adds new experimental function `tidy_model_expressions()` to create
    expressions for dataframes containing tidied results from regression model
    objects.

MAJOR CHANGES

  - Removes the redundant `bf_extractor` function. The `tidy_model_parameters`
    does the same thing.

  - Exports more utility functions (`long_to_wide_converter`, `format_num`,
    `stats_type_switch`) to get rid of reliance on `ipmisc` package.

  - To be consistent with the expressions, the dataframe for Bayesian analysis
    now also contain log of Bayes Factor values.

  - The `tidy_model_effectsize()` function is no longer exported as it is
    helpful only for the internal workings of the package.

  - Given that these values can be really high, the statistic values for
    non-parametric tests were shown on a log scale, but this is a highly
    non-standard practice that has caused a lot of confusion among users. In
    light of this feedback, the functions no longer return these values on a log
    scale but in a scientific notation to keep the statistical expressions
    short.

MINOR CHANGES

  - Removes `VR_dilemma` dataset, which lacked enough variation to be a good
    dataset to use in examples or tests.

# statsExpressions 1.1.0

MAJOR CHANGES

  - There is a new _JOSS_ paper about `{statsExpressions}` package!!
    <https://joss.theoj.org/papers/10.21105/joss.03236>

  - The effect size for independent trimmed means two-sample test has been
    changed from explanatory measure of effect size to AKP's delta, which is
    easier to understand and interpret since its a robust cousin of Cohen's
    *d*.

  - `one_sample_test` and `two_sample_test` gain `alternative` argument to
    specify alternative hypothesis (#86).

  - Cohen's *d* and Hedge's *g* use non-pooled standard deviation (cf.
    https://psyarxiv.com/tu6mp/).

MINOR CHANGES

  - The output dataframes now contain columns with additional information about
    how confidence intervals are computed (thanks to `effectsize` package).

# statsExpressions 1.0.1

BREAKING CHANGES

  - Retires all vestigial `expr_*` functions.

MINOR CHANGES

  - Adapts failing tests due to changes in `effectsize`.

# statsExpressions 1.0.0

This is the first **stable** release of `{statsExpressions}`!

There is good news and there is bad news that accompanies this milestone.

  - The **bad news**: The `API` for the package has changed **completely**: All
    functions return a *dataframe*, and not an *expression*, as a default. The
    expression is contained in a list column in the dataframe itself. So, to
    salvage your functions from breaking, you will have to add
    `$expression[[1]]` to your function calls. For example, if you were using
    the function `expr_t_onesample()`, you will now have to specify
    `expr_t_onesample()$expression[[1]]`, so on and so forth. But, in general,
    the advice is to **not** use any of the `expr_*` functions, which are
    vestigial names for new avatars of these function and will be removed in
    future. The new names are more intuitive, e.g., `expr_t_onesample()` is now
    called `one_sample_test()`, etc.

  - The **good news**: There will not be any new changes to any of the current
    functions, except for any change necessary for maintenance or bug squashing.
    Well, to be more precise, this is true only for the functions that have
    **"stable"** [badge](https://lifecycle.r-lib.org/articles/stages.html).

BUG FIXES

  - If the entered dataframe is `grouped`-tibble, the function internally
    ungroups this (#79).

MINOR CHANGES

  - To reduce dependency load, `afex` has moved from `Imports` to `Suggests`.

# statsExpressions 0.7.1

BREAKING CHANGES

  - To avoid confusion among users, the trimming level for all functions is now
    changed from `tr = 0.1` to `tr = 0.2` (which is what `WRS2` defaults to).

MAJOR CHANGES

  - `expr_template` gains a new argument `bayesian`, which can return an
    expression for Bayesian analysis, which has a slightly different template.
    Additionally, it has changed its conventions about the column names it
    expects.

  - Retires the additional caption-making functionality that was unique to
    `expr_meta_random` when `type = "parametric"`. This was the only context in
    which this feature was supported and was therefore inconsistent with the
    rest of the package API.

  - Removes `tidy_model_performance` function, which is no longer used
    internally.

  - Removes column containing `log` values of Bayes Factor as they are relevant
    only for expressions.

  - All meta-analysis packages move from `Imports` to `Suggests` to reduce the
    installation time for the user.

  - All robust tests in this package were based on trimmed means, except for
    correlation test. This has been changed: the robust correlation measure is
    now Winsorized correlation, which is based on trimming. Therefore, the
    `beta` argument has been replaced by `tr` argument. This should result only
    in minor changes in correlation coefficient estimates.

# statsExpressions 0.7.0

BREAKING CHANGES

  - To be consistent with `ggstatsplot`'s overall syntax philosophy the `type`
    argument can be used to specify which type of statistical approach is to be
    used for all functions.

    * `t_parametric`, `t_nonparametric`, `t_robust`, `t_bayes` are now removed
      in favor of a single function `two_sample_test`.

    * `expr_anova_parametric`, `expr_anova_nonparametric`, `expr_anova_robust`,
      `expr_anova_bayes` are now removed in favor of a single function
      `oneway_anova`.

  - `{statsExpressions}` no longer internally relies on `tidyBF`. All Bayesian
    analysis is carried out in this package itself. This was done to make the
    maintenance of this package easier and helps with some major internal code
    refactoring. As such, all re-exported functions from `tidyBF` have also been
    removed.

BUG FIXES

  - `contingency_table` ignored `ratio` argument while computing Cramer's *V*
    for one-sample test. This is fixed.

MAJOR CHANGES

  - All non-parametric functions now use `effectsize` package to compute effect
    sizes and not `rcompanion`. This would lead to some changes in effect sizes
    and their confidence intervals reported by the respective functions.

  - Robust one-sample test is changed from one-sample percentile bootstrap to
    bootstrap-*t* method for one-sample test, which uses trimmed mean like the
    rest of the robust functions in this package.

MINOR CHANGES

  - Package internally relies on `afex` instead of `ez` for within-subjects
    ANOVA.

  - `expr_template` gains `paired` argument.

# statsExpressions 0.6.2

MINOR CHANGES

  - Internal refactoring to catch up with changes made to `effectsize`. Tests
    are adapted to these changes as well.

  - Sample size information in expressions is pretty-formatted.

# statsExpressions 0.6.1

MAJOR CHANGES

  - Adds two new helper functions: `tidy_model_parameters` and
    `tidy_model_performance` to toggle between `easystats` and `tidymodels`
    naming conventions.

  - Drops `broomExtra` from dependencies in favor of `parameters` +
    `performance`.

  - Removes the unused and vestigial `Titanic_full` dataset.

# statsExpressions 0.6.0

BREAKING CHANGES

  - Removes the alias `expr_onesample_proptest`.

  - The `expr_template` function retires `effsize.df` argument. Now all details
    need to be entered only in `data`.

  - All meta-analyses are now carried out using `expr_meta_random` and the
    individual functions have been removed.

MAJOR CHANGES

  - All effect sizes for contingency tabs are now calculated via `effectsize`
    instead of `rcompanion`. This would lead to slight differences in effect
    sizes and their CIs but the computations will be faster. Additionally, the
    lower bound will never be negative and will be restricted to [0,1].

  - `contingency_table` function has been made less robust. It now fails instead
    of returning `NULL` when it is not supposed to work. This is done to be
    consistent with the other functions in the package which also fail instead
    of returning `NULL`.

  - `expr_anova_parametric` always applies sphericity correction for *p*-values
    for repeated measures ANOVA.

  - `expr_anova_parametric` retires non-partial variants of effect sizes
    (eta-squared and omega-squared, i.e.) for parametric analyses.

  - The *t*-test and ANOVA tests get `subject.id` argument relevant for repeated
    measures design.

MINOR CHANGES

  - Retires the vestigial `stat.title` argument. It was originally intended to
    give more info on the tests, but now the expressions themselves contain
    these details.

  - For paired ANOVA designs, `partial = TRUE` is recognized by effect sizes.

  - Retires `bias.correct` argument for contingency table analysis. It is rarely
    justifiable to use the biased version of Cramer's *V*.

# statsExpressions 0.5.1

MINOR CHANGES

  - Adapts tests to changes made in the `correlation` package.

  - Subtitles for correlation tests make clear the type of statistic.

  - Small *p*-values (< 0.001) are now shown in scientific format.

# statsExpressions 0.5.0

MINOR CHANGES

  - Adapts to changes made in `tidyBF` package.

  - Re-exports `correlation::correlation` needed for `ggstatsplot`.

  - The `t_nonparametric` subtitle now clarifies whether it's a Wilcoxon test or
    a Mann-Whitney test.

# statsExpressions 0.4.2

MINOR CHANGES

  - Thanks to Sarah, the package has a hexsticker. :)

  - Confidence intervals for Spearman's rho are computed using `correlation`
    instead of `rcompanion`.

  - All relevant functions get rid of `messages` argument as the functions no
    longer print a message when bootstrapped CIs are used.

  - The effect size measure for paired robust *t*-test is now changed to robust
    (trimmed-Winsorized) standardized difference similar to Cohen's *d*.

# statsExpressions 0.4.1

BUG FIXES

  - Major bug introduced in `0.4.0` release for `expr_anova_parametric`:
    changing `conf.level` doesn't work and function defaults to `0.90` CIs
    (#32).

MINOR CHANGES

  - Removes extra space contained in subtitles for Bayes Factor results (#31).

# statsExpressions 0.4.0

BREAKING CHANGES

  - Removes the experimental `corr_objects` function.

  - All Bayes Factor related functions have now moved to the new `tidyBF`
    package and are re-exported from there.

MAJOR CHANGES

  - Minimum R version bumped to `R 3.6.0`.

  - Retires the internal `effsize_t_parametric` helper function in favor of
    relying functions from `effectsize`, which is now added as a dependency.
    Similarly, `{statsExpressions}` now relies on `effectsize` to compute effect
    sizes for ANOVA designs, instead of `sjstats`.

  - For parametric *t*-tests and ANOVAs, confidence intervals for effect sizes
    are estimated using the noncentrality parameter method. Centrality-based
    methods are deprecated.

  - Correlation analysis is carried out using `correlation` package, which is
    now added as a dependency.

MINOR CHANGES

  - All expressions now contain name of the statistical test carried out.

# statsExpressions 0.3.1
 
  - Adds a new function `corr_objects` to reduce dependency load of
    `ggstatsplot`. This is an experimental function and should be avoided until
    it stabilizes.

# statsExpressions 0.3.0

NEW FEATURES

  - New functions to carry out meta-analysis: `expr_meta_bayes`.

# statsExpressions 0.2.1

NEW FEATURES

  - New functions to carry out meta-analysis: `expr_meta_parametric`,
    `expr_meta_robust`, `bf_meta`.

# statsExpressions 0.2.0

BREAKING CHANGES

  - `expr_template` function now expects two dataframes: `data` and `effsize.df`
    that contain the details needed for creating expressions instead of
    providing each individual values. This makes the function more friendly work
    with using modeling packages like `broom`.

MINOR CHANGES

  - Minor tweaks to how widehat is displayed in some of the expressions.

  - Cramer's *V* is bias-corrected by default.

# statsExpressions 0.1.3

MAJOR CHANGES

  - Removes `MCMCpack` from `Depends`.

  - All effect size texts now contain `^` on top to signify that these are
    estimates.

# statsExpressions 0.1.2

MINOR CHANGES

  - Maintenance release to fix additional check issues on `CRAN`.

# statsExpressions 0.1.1

MINOR CHANGES

  - Fixing tests for the new release of `rcompanion` dependency.

  - Minor code refactoring.

# statsExpressions 0.1.0

  - First release of the package.

## Test environments

* local OS X install, R 4.1.2

* ubuntu 18.04 (on github actions ci), R 4.1.2

* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

  - Fixes failing tests due to `{insight}` package update.

## revdepcheck results

One reverse dependency (`{ggstatsplot}`) is adversely affected. I am the
maintainer of this package and will be submitting an updated version to CRAN as
soon as `{statsExpressions}` update makes it to CRAN.
---
title: 'statsExpressions: R Package for Tidy Dataframes and Expressions with Statistical Details'
tags:
  - R
  - parametric statistics
  - nonparametric statistics
  - robust statistics
  - Bayesian statistics
  - tidy
authors:
  - name: Indrajeet Patil
    orcid: 0000-0003-1995-6531
    affiliation: 1
affiliations:
  - name: Center for Humans and Machines, Max Planck Institute for Human Development, Berlin, Germany
    index: 1
date: "2021-05-19"
bibliography: paper.bib
---



# Summary

The `{statsExpressions}` package has two key aims: to provide a consistent syntax
to do statistical analysis with tidy data, and to provide statistical
expressions (i.e., pre-formatted in-text statistical results) for plotting
functions. Currently, it supports common types of statistical approaches and
tests: parametric, nonparametric, robust, and Bayesian *t*-test, one-way ANOVA,
correlation analyses, contingency table analyses, and meta-analyses. The
functions are pipe-friendly and compatible with tidy data.

# Statement of need

Statistical packages exhibit substantial diversity in terms of their syntax and
expected input and output data type. For example, some functions expect vectors
as inputs, while others expect dataframes. Depending on whether it is a repeated
measures design or not, functions from the same package might expect data to be
in wide or tidy format. Some functions can internally omit missing values, while
others do not. Furthermore, the statistical test objects returned by the test
functions might not have all required information (e.g., degrees of freedom,
significance, Bayes factor, etc.) accessible in a consistent data type.
Depending on the specific test object and statistic in question, details may be
returned as a list, a matrix, an array, or a dataframe. This diversity can make
it difficult to easily access all needed information for hypothesis testing and
estimation, and to switch from one statistical approach to another.

This is where `{statsExpressions}` comes in: It can be thought of as a unified
portal through which most of the functionality in these underlying packages can
be accessed, with a simpler interface and with tidy data format.

# Comparison to Other Packages

Unlike `broom` [@Robinson2021] or `parameters` [@Lüdecke2020parameters], the
goal of `{statsExpressions}` is not to convert model objects into tidy dataframes,
but to provide a consistent and easy syntax to carry out statistical tests.
Additionally, none of these packages return statistical expressions.

# Consistent Syntax for Statistical Analysis

The package offers functions that allow users choose a statistical approach
without changing the syntax (i.e., by only specifying a single argument). The
functions always require a dataframe in tidy format [@Wickham2019], and work
with missing data. Moreover, they always return a dataframe that can be further
utilized downstream in the workflow (such as visualization).

Function | Parametric | Non-parametric | Robust | Bayesian
------------------ | ---- | ----- | ----| ----- 
`one_sample_test` | \checkmark | \checkmark | \checkmark | \checkmark
`two_sample_test` | \checkmark | \checkmark | \checkmark | \checkmark
`oneway_anova` | \checkmark | \checkmark | \checkmark | \checkmark
`corr_test` | \checkmark | \checkmark | \checkmark | \checkmark
`contingency_table` | \checkmark | \checkmark | - | \checkmark
`meta_analysis` | \checkmark | - | \checkmark | \checkmark

: A summary table listing the primary functions in the package and the
statistical approaches they support. More detailed description of the
tests and outputs from these functions can be found on the package website: <https://indrajeetpatil.github.io/statsExpressions/articles/>.

`{statsExpressions}` internally relies on `stats` package for parametric and
non-parametric [@base2021], `WRS2` package for robust [@Mair2020], and
`BayesFactor` package for Bayesian statistics [@Morey2020]. The random-effects
meta-analysis is carried out using `metafor` (parametric) [@Viechtbauer2010],
`metaplus` (robust) [@Beath2016], and `metaBMA` (Bayesian) [@Heck2019] packages.
Additionally, it relies on `easystats` packages [@Ben-Shachar2020;
@Lüdecke2020parameters;
@Lüdecke2020performance; @Lüdecke2019; @Makowski2019; @Makowski2020] to compute
appropriate effect size/posterior estimates and their confidence/credible
intervals.

# Tidy Dataframes from Statistical Analysis

To illustrate the simplicity of this syntax, let's say we want to run a one-way
ANOVA. If we first run a non-parametric ANOVA and then decide to run a robust
ANOVA instead, the syntax remains the same and the statistical approach can be
modified by changing a single argument:


```r
mtcars %>% oneway_anova(cyl, wt, type = "nonparametric") 
#> # A tibble: 1 x 14
#>   parameter1 parameter2 statistic df.error   p.value
#>   <chr>      <chr>          <dbl>    <int>     <dbl>
#> 1 wt         cyl             22.8        2 0.0000112
#>   method                       estimate conf.level conf.low conf.high
#>   <chr>                           <dbl>      <dbl>    <dbl>     <dbl>
#> 1 Kruskal-Wallis rank sum test    0.736       0.95    0.613     0.831
#>   effectsize      conf.method conf.iterations expression
#>   <chr>           <chr>                 <int> <list>    
#> 1 Epsilon2 (rank) bootstrap               100 <language>

mtcars %>% oneway_anova(cyl, wt, type = "robust")
#> # A tibble: 1 x 11
#>   statistic    df df.error p.value estimate conf.level conf.low conf.high
#>       <dbl> <dbl>    <dbl>   <dbl>    <dbl>      <dbl>    <dbl>     <dbl>
#> 1      12.7     2     12.2 0.00102     1.05       0.95    0.843      1.50
#>   effectsize                        
#>   <chr>                             
#> 1 Explanatory measure of effect size
#>   method                                            expression
#>   <chr>                                             <list>    
#> 1 A heteroscedastic one-way ANOVA for trimmed means <language>
```

These functions are also compatible with other popular data manipulation
packages (see Appendix for an example).

# Expressions for Plots

In addition to other details contained in the dataframe, there is also a column
titled `expression`, which contains a pre-formatted text with statistical
details. These expressions (Figure 1) attempt to follow the gold standard in
statistical reporting for both Bayesian [@van2020jasp] and Frequentist
[@american2019publication] frameworks.

\begin{figure}
\includegraphics[width=1\linewidth]{stats_reporting_format} \caption{The templates used in `{statsExpressions}` to display statistical details in a plot.}\label{fig:expr_template}
\end{figure}

This expression be easily displayed in a plot (Figure 2). Displaying statistical
results in the context of a visualization is indeed a philosophy adopted by the
`ggstatsplot` package [@Patil2021], and `{statsExpressions}` functions as its
statistical processing backend.

\begin{figure}
\includegraphics[width=1\linewidth]{paper_files/figure-latex/anova_example-1} \caption{Example illustrating how `{statsExpressions}` functions can be used to display results from a statistical test in a plot. Code to create this figure is reported in Appendix.}\label{fig:anova_example}
\end{figure}

# Licensing and Availability

`{statsExpressions}` is licensed under the GNU General Public License (v3.0), with all
source code stored at [GitHub](https://github.com/IndrajeetPatil/statsExpressions/).
In the spirit of honest and open science, requests and suggestions for fixes,
feature updates, as well as general questions and concerns are encouraged via
direct interaction with contributors and developers by filing an
[issue](https://github.com/IndrajeetPatil/statsExpressions/issues) while respecting
[*Contribution Guidelines*](https://indrajeetpatil.github.io/statsExpressions/CONTRIBUTING.html).

# Acknowledgements

I would like to acknowledge the support of Mina Cikara, Fiery Cushman, and Iyad
Rahwan during the development of this project. `{statsExpressions}` relies heavily
on the [`easystats`](https://github.com/easystats/easystats) ecosystem, a
collaborative project created to facilitate the usage of `R` for statistical
analyses. Thus, I would like to thank the [members of easystats](https://github.com/orgs/easystats/people) as well as the users.

# References

<div id="refs"></div>

# Appendix

## Example with `dplyr`

We can use combination of `dplyr` and `{statsExpressions}` to repeat the same
statistical analysis across grouping variables.


```r
# running one-sample proportion test for all levels of `cyl`
mtcars %>%
  group_by(cyl) %>%
  group_modify(~ contingency_table(.x, am), .keep = TRUE) %>%
  ungroup()
#> # A tibble: 3 x 13
#>     cyl statistic    df p.value method                                  
#>   <dbl>     <dbl> <dbl>   <dbl> <chr>                                   
#> 1     4     2.27      1 0.132   Chi-squared test for given probabilities
#> 2     6     0.143     1 0.705   Chi-squared test for given probabilities
#> 3     8     7.14      1 0.00753 Chi-squared test for given probabilities
#>   estimate conf.level conf.low conf.high effectsize        conf.method
#>      <dbl>      <dbl>    <dbl>     <dbl> <chr>             <chr>      
#> 1    0.344       0.95    0         0.917 Cramer's V (adj.) ncp        
#> 2    0           0.95    0         0     Cramer's V (adj.) ncp        
#> 3    0.685       0.95    0.127     1.18  Cramer's V (adj.) ncp        
#>   conf.distribution expression
#>   <chr>             <list>    
#> 1 chisq             <language>
#> 2 chisq             <language>
#> 3 chisq             <language>
```

## Code to reproduce for Figure 2


```r
# needed libraries
library(statsExpressions)
library(ggplot2)
library(ggridges)
library(palmerpenguins) # `penguins` dataset is from this package

# creating a dataframe
res <- oneway_anova(penguins, species, body_mass_g, type = "nonparametric")

# create a ridgeplot using `ggridges` package
ggplot(penguins, aes(x = body_mass_g, y = species)) +
  geom_density_ridges(
    jittered_points = TRUE,
    quantile_lines = TRUE,
    scale = 0.9,
    vline_size = 1,
    vline_color = "red",
    position = position_raincloud(adjust_vlines = TRUE)
  ) + # use 'expression' column to display results in the subtitle
  labs(
    x = "Penguin species",
    y = "Body mass (in grams)",
    title = "Kruskal-Wallis Rank Sum Test",
    subtitle = res$expression[[1]]
  )
```
# t_robust - within-subjects - without NAs

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 14
        statistic df.error p.value method                                            
            <dbl>    <dbl>   <dbl> <chr>                                             
      1      28.7       89       0 Yuen's test on trimmed means for dependent samples
        effectsize                                              estimate conf.level
        <chr>                                                      <dbl>      <dbl>
      1 Algina-Keselman-Penfield robust standardized difference     2.36       0.95
        conf.low conf.high    mu small medium large n.obs
           <dbl>     <dbl> <dbl> <dbl>  <dbl> <dbl> <int>
      1     1.96      2.61     0   0.1    0.3   0.5   150

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Yuen\"] * \"(\" * 89 * \")\" == \"28.7230\", italic(p) == \"0.0000\", widehat(delta)[\"R\"]^\"AKP\" == \"2.3582\", CI[\"95%\"] ~ \"[\" * \"1.9615\", \"2.6081\" * \"]\", italic(\"n\")[\"pairs\"] == \"150\")"

# t_robust - within-subjects - with NAs

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 14
        statistic df.error p.value method                                            
            <dbl>    <dbl>   <dbl> <chr>                                             
      1      2.91       53 0.00528 Yuen's test on trimmed means for dependent samples
        effectsize                                              estimate conf.level
        <chr>                                                      <dbl>      <dbl>
      1 Algina-Keselman-Penfield robust standardized difference    0.410       0.95
        conf.low conf.high    mu small medium large n.obs
           <dbl>     <dbl> <dbl> <dbl>  <dbl> <dbl> <int>
      1    0.238     0.611     0   0.1    0.3   0.5    90

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Yuen\"] * \"(\" * 53 * \")\" == \"2.909\", italic(p) == \"0.005\", widehat(delta)[\"R\"]^\"AKP\" == \"0.410\", CI[\"95%\"] ~ \"[\" * \"0.238\", \"0.611\" * \"]\", italic(\"n\")[\"pairs\"] == \"90\")"

# t_robust - between-subjects - without NAs

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 10
        statistic df.error   p.value
            <dbl>    <dbl>     <dbl>
      1      5.84     13.6 0.0000485
        method                                              
        <chr>                                               
      1 Yuen's test on trimmed means for independent samples
        effectsize                                              estimate conf.level
        <chr>                                                      <dbl>      <dbl>
      1 Algina-Keselman-Penfield robust standardized difference     2.48       0.99
        conf.low conf.high n.obs
           <dbl>     <dbl> <int>
      1    0.738      5.13    32

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Yuen\"] * \"(\" * 13.584 * \")\" == \"5.840\", italic(p) == \"4.846e-05\", widehat(delta)[\"R\"]^\"AKP\" == \"2.482\", CI[\"99%\"] ~ \"[\" * \"0.738\", \"5.128\" * \"]\", italic(\"n\")[\"obs\"] == \"32\")"

# t_robust - between-subjects - with NAs

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 10
        statistic df.error p.value
            <dbl>    <dbl>   <dbl>
      1     0.452     13.8   0.658
        method                                              
        <chr>                                               
      1 Yuen's test on trimmed means for independent samples
        effectsize                                              estimate conf.level
        <chr>                                                      <dbl>      <dbl>
      1 Algina-Keselman-Penfield robust standardized difference   -0.358        0.9
        conf.low conf.high n.obs
           <dbl>     <dbl> <int>
      1    -7.16     0.406    29

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Yuen\"] * \"(\" * 13.8476 * \")\" == \"0.4521\", italic(p) == \"0.6582\", widehat(delta)[\"R\"]^\"AKP\" == \"-0.3583\", CI[\"90%\"] ~ \"[\" * \"-7.1637\", \"0.4061\" * \"]\", italic(\"n\")[\"obs\"] == \"29\")"

# bayes factor (independent samples t-test)

    Code
      as.character(df$expression[[1]])
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"-0.18\", widehat(delta)[\"difference\"]^\"posterior\" == \"3.16\", CI[\"99%\"]^HDI ~ \"[\" * \"-1.35\", \"8.13\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.71\")"

# bayes factor (paired t-test)

    Code
      as.character(df$expression[[1]])
    Output
      [1] "list(atop(\"bla is ulalala\", list(log[e] * (BF[\"01\"]) == \"-3.70\", widehat(delta)[\"difference\"]^\"posterior\" == \"1.09\", CI[\"95%\"]^HDI ~ \"[\" * \"0.49\", \"1.70\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.80\")))"

# t_nonparametric works - between-subjects design

    Code
      select(df, -expression)
    Output
      # A tibble: 1 x 13
        parameter1 parameter2 statistic   p.value method                 alternative
        <chr>      <chr>          <dbl>     <dbl> <chr>                  <chr>      
      1 wt         am              230. 0.0000435 Wilcoxon rank sum test two.sided  
        effectsize        estimate conf.level conf.low conf.high conf.method n.obs
        <chr>                <dbl>      <dbl>    <dbl>     <dbl> <chr>       <int>
      1 r (rank biserial)    0.866        0.9    0.749     0.931 normal         32

---

    Code
      as.character(df$expression[[1]])
    Output
      [1] "list(italic(\"W\")[\"Mann-Whitney\"] == \"230.500\", italic(p) == \"4.347e-05\", widehat(italic(\"r\"))[\"biserial\"]^\"rank\" == \"0.866\", CI[\"90%\"] ~ \"[\" * \"0.749\", \"0.931\" * \"]\", italic(\"n\")[\"obs\"] == \"32\")"

# t_nonparametric works - within-subjects design

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 13
        parameter1 parameter2 statistic p.value method                    alternative
        <chr>      <chr>          <dbl>   <dbl> <chr>                     <chr>      
      1 length     type              10 0.00295 Wilcoxon signed rank test two.sided  
        effectsize        estimate conf.level conf.low conf.high conf.method n.obs
        <chr>                <dbl>      <dbl>    <dbl>     <dbl> <chr>       <int>
      1 r (rank biserial)   -0.853       0.99   -0.964    -0.489 normal         16

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(italic(\"V\")[\"Wilcoxon\"] == \"10.00000\", italic(p) == \"0.00295\", widehat(italic(\"r\"))[\"biserial\"]^\"rank\" == \"-0.85294\", CI[\"99%\"] ~ \"[\" * \"-0.96399\", \"-0.48865\" * \"]\", italic(\"n\")[\"pairs\"] == \"16\")"

# meta_analysis works - robust

    Code
      select(df, -expression)
    Output
      # A tibble: 1 x 11
        term    effectsize                     estimate std.error conf.low conf.high
        <chr>   <chr>                             <dbl>     <dbl>    <dbl>     <dbl>
      1 Overall meta-analytic summary estimate   -0.693     0.368    -1.56    -0.118
        statistic p.value conf.level method                                n.obs
            <dbl>   <dbl>      <dbl> <chr>                                 <int>
      1     -1.88  0.0230       0.95 Robust meta-analysis using 'metaplus'     6

---

    Code
      as.character(df$expression[[1]])
    Output
      [1] "list(italic(\"z\") == \"-1.88\", italic(p) == \"0.02\", widehat(beta)[\"summary\"]^\"meta\" == \"-0.69\", CI[\"95%\"] ~ \"[\" * \"-1.56\", \"-0.12\" * \"]\", italic(\"n\")[\"effects\"] == \"6\")"

# parametric t-test works (between-subjects without NAs)

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 17
        term  group mean.group1 mean.group2 statistic df.error p.value
        <chr> <chr>       <dbl>       <dbl>     <dbl>    <dbl>   <dbl>
      1 len   supp         20.7        17.0      1.92       58  0.0604
        method            alternative effectsize estimate conf.level conf.low
        <chr>             <chr>       <chr>         <dbl>      <dbl>    <dbl>
      1 Two Sample t-test two.sided   Cohen's d     0.495       0.99   -0.184
        conf.high conf.method conf.distribution n.obs
            <dbl> <chr>       <chr>             <int>
      1      1.17 ncp         t                    60

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Student\"] * \"(\" * 58 * \")\" == \"1.91527\", italic(p) == \"0.06039\", widehat(italic(\"d\"))[\"Cohen\"] == \"0.49452\", CI[\"99%\"] ~ \"[\" * \"-0.18354\", \"1.16839\" * \"]\", italic(\"n\")[\"obs\"] == \"60\")"

# parametric t-test works (between-subjects with NAs)

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 17
        term  group mean.group1 mean.group2 statistic df.error p.value
        <chr> <chr>       <dbl>       <dbl>     <dbl>    <dbl>   <dbl>
      1 len   supp         20.7        17.0      1.92     55.3  0.0606
        method                  alternative effectsize estimate conf.level conf.low
        <chr>                   <chr>       <chr>         <dbl>      <dbl>    <dbl>
      1 Welch Two Sample t-test two.sided   Hedges' g     0.488        0.9   0.0599
        conf.high conf.method conf.distribution n.obs
            <dbl> <chr>       <chr>             <int>
      1     0.911 ncp         t                    60

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Welch\"] * \"(\" * 55.309 * \")\" == \"1.915\", italic(p) == \"0.061\", widehat(italic(\"g\"))[\"Hedges\"] == \"0.488\", CI[\"90%\"] ~ \"[\" * \"0.060\", \"0.911\" * \"]\", italic(\"n\")[\"obs\"] == \"60\")"

# parametric t-test works (within-subjects without NAs)

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 15
        term  group     statistic df.error  p.value method        alternative
        <chr> <chr>         <dbl>    <dbl>    <dbl> <chr>         <chr>      
      1 value condition      34.8      149 1.85e-73 Paired t-test two.sided  
        effectsize estimate conf.level conf.low conf.high conf.method
        <chr>         <dbl>      <dbl>    <dbl>     <dbl> <chr>      
      1 Hedges' g      2.83        0.5     2.71      2.96 ncp        
        conf.distribution n.obs
        <chr>             <int>
      1 t                   150

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Student\"] * \"(\" * 149 * \")\" == \"34.8152\", italic(p) == \"1.8496e-73\", widehat(italic(\"g\"))[\"Hedges\"] == \"2.8283\", CI[\"50%\"] ~ \"[\" * \"2.7086\", \"2.9560\" * \"]\", italic(\"n\")[\"pairs\"] == \"150\")"

# parametric t-test works (within-subjects with NAs)

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 15
        term   group     statistic df.error  p.value method        alternative
        <chr>  <chr>         <dbl>    <dbl>    <dbl> <chr>         <chr>      
      1 desire condition      3.61       89 0.000500 Paired t-test two.sided  
        effectsize estimate conf.level conf.low conf.high conf.method
        <chr>         <dbl>      <dbl>    <dbl>     <dbl> <chr>      
      1 Cohen's d     0.381       0.95    0.167     0.597 ncp        
        conf.distribution n.obs
        <chr>             <int>
      1 t                    90

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Student\"] * \"(\" * 89 * \")\" == \"3.613\", italic(p) == \"5.000e-04\", widehat(italic(\"d\"))[\"Cohen\"] == \"0.381\", CI[\"95%\"] ~ \"[\" * \"0.167\", \"0.597\" * \"]\", italic(\"n\")[\"pairs\"] == \"90\")"

# expr_anova_robust works - between-subjects

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 11
        statistic    df df.error   p.value
            <dbl> <dbl>    <dbl>     <dbl>
      1      20.2     2     19.0 0.0000196
        method                                           
        <chr>                                            
      1 A heteroscedastic one-way ANOVA for trimmed means
        effectsize                         estimate conf.level conf.low conf.high
        <chr>                                 <dbl>      <dbl>    <dbl>     <dbl>
      1 Explanatory measure of effect size    0.859       0.95    0.853     0.864
        n.obs
        <int>
      1    32

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"F\")[\"trimmed-means\"](2, 18.97383) == \"20.24946\", italic(p) == \"0.00002\", widehat(xi) == \"0.85858\", CI[\"95%\"] ~ \"[\" * \"0.85268\", \"0.86448\" * \"]\", italic(\"n\")[\"obs\"] == \"32\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 11
        statistic    df df.error p.value
            <dbl> <dbl>    <dbl>   <dbl>
      1    0.0503     2     21.7   0.951
        method                                           
        <chr>                                            
      1 A heteroscedastic one-way ANOVA for trimmed means
        effectsize                         estimate conf.level conf.low conf.high
        <chr>                                 <dbl>      <dbl>    <dbl>     <dbl>
      1 Explanatory measure of effect size    0.201       0.99   0.0872     0.754
        n.obs
        <int>
      1    71

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(italic(\"F\")[\"trimmed-means\"](2, 21.6869) == \"0.0503\", italic(p) == \"0.9511\", widehat(xi) == \"0.2013\", CI[\"99%\"] ~ \"[\" * \"0.0872\", \"0.7537\" * \"]\", italic(\"n\")[\"obs\"] == \"71\")"

# expr_anova_robust works - within-subjects

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 11
        statistic    df df.error  p.value
            <dbl> <dbl>    <dbl>    <dbl>
      1      21.0  2.73     145. 1.15e-10
        method                                                             
        <chr>                                                              
      1 A heteroscedastic one-way repeated measures ANOVA for trimmed means
        effectsize                                                      estimate
        <chr>                                                              <dbl>
      1 Algina-Keselman-Penfield robust standardized difference average    0.664
        conf.level conf.low conf.high n.obs
             <dbl>    <dbl>     <dbl> <int>
      1       0.95    0.466     0.971    88

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"F\")[\"trimmed-means\"](2.7303, 144.7051) == \"20.9752\", italic(p) == \"1.1462e-10\", widehat(delta)[\"R-avg\"]^\"AKP\" == \"0.6635\", CI[\"95%\"] ~ \"[\" * \"0.4660\", \"0.9707\" * \"]\", italic(\"n\")[\"pairs\"] == \"88\")"

# meta_analysis works - bayesian

    Code
      dplyr::select(df, -expression)
    Output
      # A tibble: 2 x 15
        term    effectsize                       estimate std.error conf.level
        <chr>   <chr>                               <dbl>     <dbl>      <dbl>
      1 Overall meta-analytic posterior estimate    0.596     0.133       0.95
      2 tau     meta-analytic posterior estimate    0.270     0.122       0.95
        conf.low conf.high  bf10 component prior.distribution prior.location
           <dbl>     <dbl> <dbl> <chr>     <chr>                       <dbl>
      1    0.321     0.854  36.1 meta      Student's t                     0
      2    0.106     0.496  36.1 meta      Inverse gamma                   1
        prior.scale method                                 log_e_bf10 n.obs
              <dbl> <chr>                                       <dbl> <int>
      1       0.707 Bayesian meta-analysis using 'metaBMA'       3.59     5
      2       0.15  Bayesian meta-analysis using 'metaBMA'       3.59     5

---

    Code
      as.character(df$expression[[1]])
    Output
      [1] "list(atop(\"ayyo arecha\", list(log[e] * (BF[\"01\"]) == \"-3.587\", widehat(delta)[\"difference\"]^\"posterior\" == \"0.596\", CI[\"95%\"]^HDI ~ \"[\" * \"0.321\", \"0.854\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.707\")))"

# contingency_table works

    Code
      df
    Output
      # A tibble: 12 x 16
         Species    Sepal.Length n_obs variable     std.dev   iqr conf.low conf.high
         <fct>             <dbl> <int> <chr>          <dbl> <dbl>    <dbl>     <dbl>
       1 setosa             5.01    50 Sepal.Length   0.352 0.400     4.90      5.10
       2 versicolor         5.94    50 Sepal.Length   0.516 0.7       5.80      6.07
       3 virginica          6.59    50 Sepal.Length   0.636 0.750     6.39      6.79
       4 setosa             5       50 Sepal.Length  NA     0.400     4.9       5.1 
       5 versicolor         5.9     50 Sepal.Length  NA     0.7       5.7       6.1 
       6 virginica          6.5     50 Sepal.Length  NA     0.750     6.4       6.7 
       7 setosa             5       50 Sepal.Length   0.352 0.400     4.92      5.10
       8 versicolor         5.91    50 Sepal.Length   0.516 0.7       5.77      6.09
       9 virginica          6.55    50 Sepal.Length   0.636 0.750     6.37      6.74
      10 setosa             5.02    50 Sepal.Length  NA     0.400     4.91      5.11
      11 versicolor         5.75    50 Sepal.Length  NA     0.7       5.57      6.13
      12 virginica          6.40    50 Sepal.Length  NA     0.750     6.26      6.53
           min   max skewness kurtosis n.missing expression                   
         <dbl> <dbl>    <dbl>    <dbl>     <int> <glue>                       
       1   4.3   5.8    0.120  -0.253          0 widehat(mu)[mean]=='5.01'    
       2   4.9   7      0.105  -0.533          0 widehat(mu)[mean]=='5.94'    
       3   4.9   7.9    0.118   0.0329         0 widehat(mu)[mean]=='6.59'    
       4   4.3   5.8    0.120  -0.253          0 widehat(mu)[median]=='5.000' 
       5   4.9   7      0.105  -0.533          0 widehat(mu)[median]=='5.900' 
       6   4.9   7.9    0.118   0.0329         0 widehat(mu)[median]=='6.500' 
       7   4.3   5.8    0.120  -0.253          0 widehat(mu)[trimmed]=='5.000'
       8   4.9   7      0.105  -0.533          0 widehat(mu)[trimmed]=='5.910'
       9   4.9   7.9    0.118   0.0329         0 widehat(mu)[trimmed]=='6.547'
      10   4.3   5.8    0.120  -0.253          0 widehat(mu)[MAP]=='5.02'     
      11   4.9   7      0.105  -0.533          0 widehat(mu)[MAP]=='5.75'     
      12   4.9   7.9    0.118   0.0329         0 widehat(mu)[MAP]=='6.40'     
         n_label                   mad
         <chr>                   <dbl>
       1 "setosa\n(n = 50)"     NA    
       2 "versicolor\n(n = 50)" NA    
       3 "virginica\n(n = 50)"  NA    
       4 "setosa\n(n = 50)"      0.297
       5 "versicolor\n(n = 50)"  0.519
       6 "virginica\n(n = 50)"   0.593
       7 "setosa\n(n = 50)"     NA    
       8 "versicolor\n(n = 50)" NA    
       9 "virginica\n(n = 50)"  NA    
      10 "setosa\n(n = 50)"     NA    
      11 "versicolor\n(n = 50)" NA    
      12 "virginica\n(n = 50)"  NA    

---

    Code
      df_na
    Output
      # A tibble: 16 x 16
         condition desire n_obs variable std.dev   iqr conf.low conf.high   min   max
         <fct>      <dbl> <int> <chr>      <dbl> <dbl>    <dbl>     <dbl> <dbl> <dbl>
       1 HDHF        7.85    92 desire      2.47   4       7.42      8.33   0      10
       2 HDLF        6.74    91 desire      3.11   5       6.15      7.33   0      10
       3 LDHF        7.38    91 desire      2.52   3.5     6.96      7.95   0.5    10
       4 LDLF        5.72    93 desire      2.71   4       5.14      6.27   0      10
       5 HDHF        8.75    92 desire     NA      4       8         9.5    0      10
       6 HDLF        8       91 desire     NA      5       6.5       8.5    0      10
       7 LDHF        8       91 desire     NA      3.5     7.5       8.5    0.5    10
       8 LDLF        6       93 desire     NA      4       5         6      0      10
       9 HDHF        8.47    92 desire      2.47   4       7.71      8.74   0      10
      10 HDLF        7.32    91 desire      3.11   5       6.23      7.81   0      10
      11 LDHF        7.88    91 desire      2.52   3.5     7.05      8.17   0.5    10
      12 LDLF        5.72    93 desire      2.71   4       5.20      6.44   0      10
      13 HDHF        9.98    92 desire     NA      4       9.90     10      0      10
      14 HDLF        9.73    91 desire     NA      5       8.39      9.99   0      10
      15 LDHF        9.85    91 desire     NA      3.5     7.86      9.99   0.5    10
      16 LDLF        5.99    93 desire     NA      4       2.80      8.52   0      10
         skewness kurtosis n.missing expression                    n_label         
            <dbl>    <dbl>     <int> <glue>                        <chr>           
       1   -1.13     0.486         0 widehat(mu)[mean]=='7.85'     "HDHF\n(n = 92)"
       2   -0.740   -0.663         0 widehat(mu)[mean]=='6.74'     "HDLF\n(n = 91)"
       3   -0.947    0.160         0 widehat(mu)[mean]=='7.38'     "LDHF\n(n = 91)"
       4   -0.132   -0.761         0 widehat(mu)[mean]=='5.72'     "LDLF\n(n = 93)"
       5   -1.13     0.486         0 widehat(mu)[median]=='8.750'  "HDHF\n(n = 92)"
       6   -0.740   -0.663         0 widehat(mu)[median]=='8.000'  "HDLF\n(n = 91)"
       7   -0.947    0.160         0 widehat(mu)[median]=='8.000'  "LDHF\n(n = 91)"
       8   -0.132   -0.761         0 widehat(mu)[median]=='6.000'  "LDLF\n(n = 93)"
       9   -1.13     0.486         0 widehat(mu)[trimmed]=='8.473' "HDHF\n(n = 92)"
      10   -0.740   -0.663         0 widehat(mu)[trimmed]=='7.318' "HDLF\n(n = 91)"
      11   -0.947    0.160         0 widehat(mu)[trimmed]=='7.882' "LDHF\n(n = 91)"
      12   -0.132   -0.761         0 widehat(mu)[trimmed]=='5.719' "LDLF\n(n = 93)"
      13   -1.13     0.486         0 widehat(mu)[MAP]=='9.98'      "HDHF\n(n = 92)"
      14   -0.740   -0.663         0 widehat(mu)[MAP]=='9.73'      "HDLF\n(n = 91)"
      15   -0.947    0.160         0 widehat(mu)[MAP]=='9.85'      "LDHF\n(n = 91)"
      16   -0.132   -0.761         0 widehat(mu)[MAP]=='5.99'      "LDLF\n(n = 93)"
           mad
         <dbl>
       1 NA   
       2 NA   
       3 NA   
       4 NA   
       5  1.85
       6  2.97
       7  2.97
       8  2.97
       9 NA   
      10 NA   
      11 NA   
      12 NA   
      13 NA   
      14 NA   
      15 NA   
      16 NA   

# between-subjects - data with and without NAs

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 14
        parameter1 parameter2 statistic df.error      p.value
        <chr>      <chr>          <dbl>    <int>        <dbl>
      1 length     genre           51.4        8 0.0000000217
        method                       effectsize      estimate conf.level conf.low
        <chr>                        <chr>              <dbl>      <dbl>    <dbl>
      1 Kruskal-Wallis rank sum test Epsilon2 (rank)    0.328       0.95    0.258
        conf.high conf.method          conf.iterations n.obs
            <dbl> <chr>                          <int> <int>
      1         1 percentile bootstrap             100   158

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(chi[\"Kruskal-Wallis\"]^2 * \"(\" * 8 * \")\" == \"51.42672\", italic(p) == \"2.17135e-08\", widehat(epsilon)[\"ordinal\"]^2 == \"0.32756\", CI[\"95%\"] ~ \"[\" * \"0.25829\", \"1.00000\" * \"]\", italic(\"n\")[\"obs\"] == \"158\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 14
        parameter1  parameter2 statistic df.error p.value method                      
        <chr>       <chr>          <dbl>    <int>   <dbl> <chr>                       
      1 sleep_cycle vore            5.24        3   0.155 Kruskal-Wallis rank sum test
        effectsize      estimate conf.level conf.low conf.high conf.method         
        <chr>              <dbl>      <dbl>    <dbl>     <dbl> <chr>               
      1 Epsilon2 (rank)    0.175       0.99   0.0445         1 percentile bootstrap
        conf.iterations n.obs
                  <int> <int>
      1             100    31

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(chi[\"Kruskal-Wallis\"]^2 * \"(\" * 3 * \")\" == \"5.240\", italic(p) == \"0.155\", widehat(epsilon)[\"ordinal\"]^2 == \"0.175\", CI[\"99%\"] ~ \"[\" * \"0.045\", \"1.000\" * \"]\", italic(\"n\")[\"obs\"] == \"31\")"

# within-subjects - data with and without NAs

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(chi[\"Friedman\"]^2 * \"(\" * 3 * \")\" == \"55.8338\", italic(p) == \"4.5584e-12\", widehat(italic(\"W\"))[\"Kendall\"] == \"0.1750\", CI[\"99%\"] ~ \"[\" * \"0.1142\", \"1.0000\" * \"]\", italic(\"n\")[\"pairs\"] == \"88\")"

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(chi[\"Friedman\"]^2 * \"(\" * 3 * \")\" == \"410.000\", italic(p) == \"1.510e-88\", widehat(italic(\"W\"))[\"Kendall\"] == \"0.911\", CI[\"90%\"] ~ \"[\" * \"0.906\", \"1.000\" * \"]\", italic(\"n\")[\"pairs\"] == \"150\")"

# parametric anova subtitles work (without NAs)

    Code
      select(df, -expression)
    Output
      # A tibble: 1 x 13
        statistic    df df.error   p.value
            <dbl> <dbl>    <dbl>     <dbl>
      1      20.2     2     19.0 0.0000196
        method                                                   effectsize estimate
        <chr>                                                    <chr>         <dbl>
      1 One-way analysis of means (not assuming equal variances) Eta2          0.681
        conf.level conf.low conf.high conf.method conf.distribution n.obs
             <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1       0.95    0.437         1 ncp         F                    32

---

    Code
      as.character(df$expression[[1]])
    Output
      [1] "list(italic(\"F\")[\"Welch\"](2, 18.97383) == \"20.24946\", italic(p) == \"0.00002\", widehat(eta[\"p\"]^2) == \"0.68097\", CI[\"95%\"] ~ \"[\" * \"0.43668\", \"1.00000\" * \"]\", italic(\"n\")[\"obs\"] == \"32\")"

---

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 13
        statistic    df df.error    p.value method                    effectsize
            <dbl> <dbl>    <dbl>      <dbl> <chr>                     <chr>     
      1      22.9     2       29 0.00000107 One-way analysis of means Eta2      
        estimate conf.level conf.low conf.high conf.method conf.distribution n.obs
           <dbl>      <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1    0.612       0.95    0.404         1 ncp         F                    32

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"F\")[\"Fisher\"](2, 29) == \"22.91139\", italic(p) == \"1.07468e-06\", widehat(eta[\"p\"]^2) == \"0.61242\", CI[\"95%\"] ~ \"[\" * \"0.40360\", \"1.00000\" * \"]\", italic(\"n\")[\"obs\"] == \"32\")"

# parametric anova subtitles with partial omega-squared

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 13
        statistic    df df.error p.value
            <dbl> <dbl>    <dbl>   <dbl>
      1      2.27     3     24.0   0.107
        method                                                   effectsize estimate
        <chr>                                                    <chr>         <dbl>
      1 One-way analysis of means (not assuming equal variances) Omega2        0.119
        conf.level conf.low conf.high conf.method conf.distribution n.obs
             <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1       0.95        0         1 ncp         F                    51

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"F\")[\"Welch\"](3, 24.0475) == \"2.2653\", italic(p) == \"0.1066\", widehat(omega[\"p\"]^2) == \"0.1192\", CI[\"95%\"] ~ \"[\" * \"0.0000\", \"1.0000\" * \"]\", italic(\"n\")[\"obs\"] == \"51\")"

# paired parametric anova subtitles work (without NAs)

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 17
        term      sumsq sum.squares.error    df df.error meansq statistic  p.value
        <chr>     <dbl>             <dbl> <dbl>    <dbl>  <dbl>     <dbl>    <dbl>
      1 condition 1656.              318.  1.15     171.   1.86      776. 1.32e-69
        method                                              effectsize       estimate
        <chr>                                               <chr>               <dbl>
      1 ANOVA estimation for factorial designs using 'afex' Omega2 (partial)    0.707
        conf.level conf.low conf.high conf.method conf.distribution n.obs
             <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1       0.99    0.658         1 ncp         F                   150

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"F\")[\"Fisher\"](1.149, 171.217) == \"776.318\", italic(p) == \"1.325e-69\", widehat(omega[\"p\"]^2) == \"0.707\", CI[\"99%\"] ~ \"[\" * \"0.658\", \"1.000\" * \"]\", italic(\"n\")[\"pairs\"] == \"150\")"

# too few obs

    Code
      as.character(p_sub$expression[[1]])
    Output
      [1] "list(italic(\"F\")[\"Fisher\"](6, 24) == \"43.14\", italic(p) == \"1.08e-11\", widehat(eta[\"p\"]^2) == \"0.92\", CI[\"95%\"] ~ \"[\" * \"0.85\", \"1.00\" * \"]\", italic(\"n\")[\"pairs\"] == \"5\")"

# long_to_wide_converter works - spread true

    Code
      list(df1, df2, df3, df4)
    Output
      [[1]]
      # A tibble: 150 x 5
         rowid Petal.Length Petal.Width Sepal.Length Sepal.Width
         <int>        <dbl>       <dbl>        <dbl>       <dbl>
       1     1          1.4         0.2          5.1         3.5
       2     2          1.4         0.2          4.9         3  
       3     3          1.3         0.2          4.7         3.2
       4     4          1.5         0.2          4.6         3.1
       5     5          1.4         0.2          5           3.6
       6     6          1.7         0.4          5.4         3.9
       7     7          1.4         0.3          4.6         3.4
       8     8          1.5         0.2          5           3.4
       9     9          1.4         0.2          4.4         2.9
      10    10          1.5         0.1          4.9         3.1
      # ... with 140 more rows
      
      [[2]]
      # A tibble: 32 x 3
         rowid am       wt
         <int> <fct> <dbl>
       1     1 0      3.22
       2     2 0      3.44
       3     3 0      3.46
       4     4 0      3.57
       5     5 0      3.19
       6     6 0      3.15
       7     7 0      3.44
       8     8 0      3.44
       9     9 0      4.07
      10    10 0      3.73
      # ... with 22 more rows
      
      [[3]]
      # A tibble: 88 x 5
         rowid  HDHF  HDLF  LDHF  LDLF
         <int> <dbl> <dbl> <dbl> <dbl>
       1     1  10     9     6     6  
       2     3  10    10    10     5  
       3     4   9     6     9     6  
       4     5   8.5   5.5   6.5   3  
       5     6   3     7.5   0.5   2  
       6     7  10    10    10    10  
       7     8  10     9    10    10  
       8     9  10     6     9.5   9.5
       9    11   0     0     2.5   0  
      10    12  10     8.5   7.5   9.5
      # ... with 78 more rows
      
      [[4]]
      # A tibble: 51 x 3
         rowid vore  brainwt
         <int> <fct>   <dbl>
       1     3 carni  0.07  
       2     4 carni  0.0108
       3     5 carni  0.0256
       4     7 carni  0.325 
       5     9 carni  0.0125
       6    12 carni  0.157 
       7    17 carni  0.0175
       8    18 carni  0.0445
       9    19 carni  0.0504
      10    21 herbi  0.423 
      # ... with 41 more rows
      

# long_to_wide_converter works - spread false

    Code
      list(df1, df2, df3, df4)
    Output
      [[1]]
      # A tibble: 600 x 3
         rowid condition    value
         <int> <fct>        <dbl>
       1     1 Petal.Length   1.4
       2     1 Petal.Width    0.2
       3     1 Sepal.Length   5.1
       4     1 Sepal.Width    3.5
       5     2 Petal.Length   1.4
       6     2 Petal.Width    0.2
       7     2 Sepal.Length   4.9
       8     2 Sepal.Width    3  
       9     3 Petal.Length   1.3
      10     3 Petal.Width    0.2
      # ... with 590 more rows
      
      [[2]]
      # A tibble: 32 x 3
         rowid am       wt
         <int> <fct> <dbl>
       1     1 0      3.22
       2     2 0      3.44
       3     3 0      3.46
       4     4 0      3.57
       5     5 0      3.19
       6     6 0      3.15
       7     7 0      3.44
       8     8 0      3.44
       9     9 0      4.07
      10    10 0      3.73
      # ... with 22 more rows
      
      [[3]]
      # A tibble: 352 x 3
         rowid condition desire
         <int> <fct>      <dbl>
       1     1 HDHF          10
       2     1 HDLF           9
       3     1 LDHF           6
       4     1 LDLF           6
       5     3 HDHF          10
       6     3 HDLF          10
       7     3 LDHF          10
       8     3 LDLF           5
       9     4 HDHF           9
      10     4 HDLF           6
      # ... with 342 more rows
      
      [[4]]
      # A tibble: 51 x 3
         rowid vore  brainwt
         <int> <fct>   <dbl>
       1     3 carni  0.07  
       2     4 carni  0.0108
       3     5 carni  0.0256
       4     7 carni  0.325 
       5     9 carni  0.0125
       6    12 carni  0.157 
       7    17 carni  0.0175
       8    18 carni  0.0445
       9    19 carni  0.0504
      10    21 herbi  0.423 
      # ... with 41 more rows
      

# tidy_model_expressions works - F

    Code
      select(df1, -label)
    Output
      # A tibble: 3 x 11
        term  statistic    df df.error p.value group   sumsq meansq  estimate
        <chr>     <dbl> <dbl>    <dbl>   <dbl> <chr>   <dbl>  <dbl>     <dbl>
      1 N         9.04      1       15 0.00885 Within 189.   189.    0.184   
      2 P         0.401     1       15 0.536   Within   8.40   8.40 -0.0171  
      3 N:P       1.02      1       15 0.329   Within  21.3   21.3   0.000457
        sum.squares.error mean.square.error
                    <dbl>             <dbl>
      1              314.              314.
      2              314.              314.
      3              314.              314.

---

    Code
      df1$label
    Output
      list(widehat(italic(omega)[p]^2)=='0.18', italic(F)('1', '15')=='9.04', italic(p)=='8.85e-03')
      list(widehat(italic(omega)[p]^2)=='-0.02', italic(F)('1', '15')=='0.40', italic(p)=='0.54')
      list(widehat(italic(omega)[p]^2)=='4.57e-04', italic(F)('1', '15')=='1.02', italic(p)=='0.33')

---

    Code
      select(df2, -label)
    Output
      # A tibble: 3 x 11
        term  statistic    df df.error p.value group   sumsq meansq estimate
        <chr>     <dbl> <dbl>    <dbl>   <dbl> <chr>   <dbl>  <dbl>    <dbl>
      1 N         9.04      1       15 0.00885 Within 189.   189.     0.376 
      2 P         0.401     1       15 0.536   Within   8.40   8.40   0.0261
      3 N:P       1.02      1       15 0.329   Within  21.3   21.3    0.0635
        sum.squares.error mean.square.error
                    <dbl>             <dbl>
      1              314.              314.
      2              314.              314.
      3              314.              314.

---

    Code
      df2$label
    Output
      list(widehat(italic(eta)[p]^2)=='0.38', italic(F)('1', '15')=='9.04', italic(p)=='8.85e-03')
      list(widehat(italic(eta)[p]^2)=='0.03', italic(F)('1', '15')=='0.40', italic(p)=='0.54')
      list(widehat(italic(eta)[p]^2)=='0.06', italic(F)('1', '15')=='1.02', italic(p)=='0.33')

# one_sample_test parametric works

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 14
           mu statistic df.error p.value method            alternative effectsize
        <dbl>     <dbl>    <dbl>   <dbl> <chr>             <chr>       <chr>     
      1   120     -2.67       78 0.00910 One Sample t-test two.sided   Hedges' g 
        estimate conf.level conf.low conf.high conf.method conf.distribution n.obs
           <dbl>      <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1   -0.298       0.95   -0.524   -0.0743 ncp         t                    79

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Student\"] * \"(\" * 78 * \")\" == \"-2.67496\", italic(p) == \"0.00910\", widehat(italic(\"g\"))[\"Hedges\"] == \"-0.29805\", CI[\"95%\"] ~ \"[\" * \"-0.52379\", \"-0.07429\" * \"]\", italic(\"n\")[\"obs\"] == \"79\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 14
           mu statistic df.error p.value method            alternative effectsize
        <dbl>     <dbl>    <dbl>   <dbl> <chr>             <chr>       <chr>     
      1   120     -2.67       78 0.00910 One Sample t-test two.sided   Cohen's d 
        estimate conf.level conf.low conf.high conf.method conf.distribution n.obs
           <dbl>      <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1   -0.301        0.9   -0.492    -0.111 ncp         t                    79

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Student\"] * \"(\" * 78 * \")\" == \"-2.6750\", italic(p) == \"0.0091\", widehat(italic(\"d\"))[\"Cohen\"] == \"-0.3010\", CI[\"90%\"] ~ \"[\" * \"-0.4924\", \"-0.1115\" * \"]\", italic(\"n\")[\"obs\"] == \"79\")"

# one_sample_test non-parametric works

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 11
        statistic p.value method                    alternative effectsize       
            <dbl>   <dbl> <chr>                     <chr>       <chr>            
      1      754.   0.323 Wilcoxon signed rank test two.sided   r (rank biserial)
        estimate conf.level conf.low conf.high conf.method n.obs
           <dbl>      <dbl>    <dbl>     <dbl> <chr>       <int>
      1   -0.149       0.95   -0.416     0.143 normal         60

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"V\")[\"Wilcoxon\"] == \"753.5000\", italic(p) == \"0.3227\", widehat(italic(\"r\"))[\"biserial\"]^\"rank\" == \"-0.1486\", CI[\"95%\"] ~ \"[\" * \"-0.4162\", \"0.1427\" * \"]\", italic(\"n\")[\"obs\"] == \"60\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 11
        statistic   p.value method                    alternative effectsize       
            <dbl>     <dbl> <chr>                     <chr>       <chr>            
      1       262 0.0000125 Wilcoxon signed rank test two.sided   r (rank biserial)
        estimate conf.level conf.low conf.high conf.method n.obs
           <dbl>      <dbl>    <dbl>     <dbl> <chr>       <int>
      1   -0.672       0.95   -0.806    -0.472 normal         56

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(italic(\"V\")[\"Wilcoxon\"] == \"262.0000\", italic(p) == \"1.2527e-05\", widehat(italic(\"r\"))[\"biserial\"]^\"rank\" == \"-0.6717\", CI[\"95%\"] ~ \"[\" * \"-0.8058\", \"-0.4720\" * \"]\", italic(\"n\")[\"obs\"] == \"56\")"

# one_sample_test robust works

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 9
        statistic p.value n.obs method                                 effectsize  
            <dbl>   <dbl> <int> <chr>                                  <chr>       
      1     0.787   0.455    11 Bootstrap-t method for one-sample test Trimmed mean
        estimate conf.level conf.low conf.high
           <dbl>      <dbl>    <dbl>     <dbl>
      1        9        0.9     6.55      11.5

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"bootstrapped\"] == \"0.7866\", italic(p) == \"0.4550\", widehat(mu)[\"trimmed\"] == \"9.0000\", CI[\"90%\"] ~ \"[\" * \"6.5487\", \"11.4513\" * \"]\", italic(\"n\")[\"obs\"] == \"11\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 9
        statistic p.value n.obs method                                 effectsize  
            <dbl>   <dbl> <int> <chr>                                  <chr>       
      1     -3.81    0.04    56 Bootstrap-t method for one-sample test Trimmed mean
        estimate conf.level conf.low conf.high
           <dbl>      <dbl>    <dbl>     <dbl>
      1   0.0390       0.99  -0.0669     0.145

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"bootstrapped\"] == \"-3.8075\", italic(p) == \"0.0400\", widehat(mu)[\"trimmed\"] == \"0.0390\", CI[\"99%\"] ~ \"[\" * \"-0.0669\", \"0.1448\" * \"]\", italic(\"n\")[\"obs\"] == \"56\")"

# one_sample_test bayes factor works

    Code
      names(df_results)
    Output
       [1] "term"               "effectsize"         "estimate"          
       [4] "conf.level"         "conf.low"           "conf.high"         
       [7] "pd"                 "rope.percentage"    "prior.distribution"
      [10] "prior.location"     "prior.scale"        "bf10"              
      [13] "method"             "log_e_bf10"         "n.obs"             
      [16] "expression"        

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"-47.84\", widehat(delta)[\"difference\"]^\"posterior\" == \"-1.76\", CI[\"90%\"]^HDI ~ \"[\" * \"-1.99\", \"-1.52\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.99\")"

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"2.125\", widehat(delta)[\"difference\"]^\"posterior\" == \"0.018\", CI[\"95%\"]^HDI ~ \"[\" * \"-0.242\", \"0.265\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.900\")"

# contingency_table works

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 12
        statistic    df p.value method                     effectsize        estimate
            <dbl> <int>   <dbl> <chr>                      <chr>                <dbl>
      1      8.74     2  0.0126 Pearson's Chi-squared test Cramer's V (adj.)    0.464
        conf.level conf.low conf.high conf.method conf.distribution n.obs
             <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1       0.99        0         1 ncp         chisq                32

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(chi[\"Pearson\"]^2 * \"(\" * 2 * \")\" == \"8.74073\", italic(p) == \"0.01265\", widehat(italic(\"V\"))[\"Cramer\"] == \"0.46431\", CI[\"99%\"] ~ \"[\" * \"0.00000\", \"1.00000\" * \"]\", italic(\"n\")[\"obs\"] == \"32\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 12
        statistic    df   p.value method                     effectsize       
            <dbl> <int>     <dbl> <chr>                      <chr>            
      1      457.     1 2.30e-101 Pearson's Chi-squared test Cramer's V (adj.)
        estimate conf.level conf.low conf.high conf.method conf.distribution n.obs
           <dbl>      <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1    0.455       0.95    0.420         1 ncp         chisq              2201

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(chi[\"Pearson\"]^2 * \"(\" * 1 * \")\" == \"456.87\", italic(p) == \"2.30e-101\", widehat(italic(\"V\"))[\"Cramer\"] == \"0.46\", CI[\"95%\"] ~ \"[\" * \"0.42\", \"1.00\" * \"]\", italic(\"n\")[\"obs\"] == \"2,201\")"

---

    Code
      select(df3, -expression)
    Output
      # A tibble: 1 x 12
        statistic    df p.value method                     effectsize        estimate
            <dbl> <int>   <dbl> <chr>                      <chr>                <dbl>
      1      15.8    15   0.399 Pearson's Chi-squared test Cramer's V (adj.)   0.0558
        conf.level conf.low conf.high conf.method conf.distribution n.obs
             <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1       0.99        0         1 ncp         chisq                52

---

    Code
      as.character(df3$expression[[1]])
    Output
      [1] "list(chi[\"Pearson\"]^2 * \"(\" * 15 * \")\" == \"15.75\", italic(p) == \"0.40\", widehat(italic(\"V\"))[\"Cramer\"] == \"0.06\", CI[\"99%\"] ~ \"[\" * \"0.00\", \"1.00\" * \"]\", italic(\"n\")[\"obs\"] == \"52\")"

# paired contingency_table works 

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 11
        statistic    df  p.value method                     effectsize estimate
            <dbl> <dbl>    <dbl> <chr>                      <chr>         <dbl>
      1      13.3     1 0.000261 McNemar's Chi-squared test Cohen's g     0.333
        conf.level conf.low conf.high conf.method n.obs
             <dbl>    <dbl>     <dbl> <chr>       <int>
      1       0.95    0.164     0.427 binomial      100

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(chi[\"McNemar\"]^2 * \"(\" * 1 * \")\" == \"13.33333\", italic(p) == \"0.00026\", widehat(italic(\"g\"))[\"Cohen\"] == \"0.33333\", CI[\"95%\"] ~ \"[\" * \"0.16436\", \"0.42663\" * \"]\", italic(\"n\")[\"pairs\"] == \"100\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 11
        statistic    df  p.value method                     effectsize estimate
            <dbl> <dbl>    <dbl> <chr>                      <chr>         <dbl>
      1      13.3     1 0.000261 McNemar's Chi-squared test Cohen's g     0.333
        conf.level conf.low conf.high conf.method n.obs
             <dbl>    <dbl>     <dbl> <chr>       <int>
      1        0.9    0.195     0.416 binomial       95

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(chi[\"McNemar\"]^2 * \"(\" * 1 * \")\" == \"13.333\", italic(p) == \"2.607e-04\", widehat(italic(\"g\"))[\"Cohen\"] == \"0.333\", CI[\"90%\"] ~ \"[\" * \"0.195\", \"0.416\" * \"]\", italic(\"n\")[\"pairs\"] == \"95\")"

# Goodness of Fit contingency_table works without counts

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 12
        statistic    df p.value method                                   effectsize 
            <dbl> <dbl>   <dbl> <chr>                                    <chr>      
      1      1.12     1   0.289 Chi-squared test for given probabilities Pearson's C
        estimate conf.level conf.low conf.high conf.method conf.distribution n.obs
           <dbl>      <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1    0.184       0.99        0         1 ncp         chisq                32

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(chi[\"gof\"]^2 * \"(\" * 1 * \")\" == \"1.12500\", italic(p) == \"0.28884\", widehat(italic(\"C\"))[\"Pearson\"] == \"0.18429\", CI[\"99%\"] ~ \"[\" * \"0.00000\", \"1.00000\" * \"]\", italic(\"n\")[\"obs\"] == \"32\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 12
        statistic    df   p.value method                                   effectsize 
            <dbl> <dbl>     <dbl> <chr>                                    <chr>      
      1      722.     1 3.92e-159 Chi-squared test for given probabilities Pearson's C
        estimate conf.level conf.low conf.high conf.method conf.distribution n.obs
           <dbl>      <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1    0.497       0.95    0.474         1 ncp         chisq              2201

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(chi[\"gof\"]^2 * \"(\" * 1 * \")\" == \"722.45\", italic(p) == \"3.92e-159\", widehat(italic(\"C\"))[\"Pearson\"] == \"0.50\", CI[\"95%\"] ~ \"[\" * \"0.47\", \"1.00\" * \"]\", italic(\"n\")[\"obs\"] == \"2,201\")"

---

    Code
      select(df3, -expression)
    Output
      # A tibble: 1 x 12
        statistic    df     p.value method                                  
            <dbl> <dbl>       <dbl> <chr>                                   
      1      33.8     3 0.000000223 Chi-squared test for given probabilities
        effectsize  estimate conf.level conf.low conf.high conf.method
        <chr>          <dbl>      <dbl>    <dbl>     <dbl> <chr>      
      1 Pearson's C    0.555       0.95    0.413         1 ncp        
        conf.distribution n.obs
        <chr>             <int>
      1 chisq                76

---

    Code
      as.character(df3$expression[[1]])
    Output
      [1] "list(chi[\"gof\"]^2 * \"(\" * 3 * \")\" == \"33.76\", italic(p) == \"2.23e-07\", widehat(italic(\"C\"))[\"Pearson\"] == \"0.55\", CI[\"95%\"] ~ \"[\" * \"0.41\", \"1.00\" * \"]\", italic(\"n\")[\"obs\"] == \"76\")"

# bayes factor (proportion test)

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 3
         bf10 prior.scale method                                     
        <dbl>       <dbl> <chr>                                      
      1 0.247           1 Bayesian one-way contingency table analysis

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"1.40\", italic(\"a\")[\"Gunel-Dickey\"] == \"1.00\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 3
         bf10 prior.scale method                                     
        <dbl>       <dbl> <chr>                                      
      1 0.579          10 Bayesian one-way contingency table analysis

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(atop(\"duh\", list(log[e] * (BF[\"01\"]) == \"0.55\", italic(\"a\")[\"Gunel-Dickey\"] == \"10.00\")))"

# bayes factor (contingency tab)

    Code
      list(as.character(expr_text1$expression[[1]]), as.character(expr_text2$
        expression[[1]]), as.character(expr_text3$expression[[1]]))
    Output
      [[1]]
      [1] "list(log[e] * (BF[\"01\"]) == \"-3.335\", widehat(italic(\"V\"))[\"Cramer\"]^\"posterior\" == \"0.479\", CI[\"89%\"]^HDI ~ \"[\" * \"0.285\", \"0.692\" * \"]\", italic(\"a\")[\"Gunel-Dickey\"] == \"1.000\")"
      
      [[2]]
      [1] "list(log[e] * (BF[\"01\"]) == \"-214.255\", widehat(italic(\"V\"))[\"Cramer\"]^\"posterior\" == \"0.455\", CI[\"99%\"]^HDI ~ \"[\" * \"0.402\", \"0.508\" * \"]\", italic(\"a\")[\"Gunel-Dickey\"] == \"1.000\")"
      
      [[3]]
      [1] "list(log[e] * (BF[\"01\"]) == \"-213.873\", widehat(italic(\"V\"))[\"Cramer\"]^\"posterior\" == \"0.454\", CI[\"95%\"]^HDI ~ \"[\" * \"0.417\", \"0.495\" * \"]\", italic(\"a\")[\"Gunel-Dickey\"] == \"1.500\")"
      

# meta_analysis works - parametric

    Code
      select(df, -expression)
    Output
      # A tibble: 1 x 11
        term    effectsize                     estimate std.error conf.level conf.low
        <chr>   <chr>                             <dbl>     <dbl>      <dbl>    <dbl>
      1 Overall meta-analytic summary estimate    0.438     0.202       0.95   0.0423
        conf.high statistic p.value method                        n.obs
            <dbl>     <dbl>   <dbl> <chr>                         <int>
      1     0.833      2.17  0.0300 Meta-analysis using 'metafor'     5

---

    Code
      as.character(df$expression[[1]])
    Output
      [1] "list(italic(\"z\") == \"2.17\", italic(p) == \"0.03\", widehat(beta)[\"summary\"]^\"meta\" == \"0.44\", CI[\"95%\"] ~ \"[\" * \"0.04\", \"0.83\" * \"]\", italic(\"n\")[\"effects\"] == \"5\")"

# tidy_model_expressions works

    Code
      select(df_t, -label)
    Output
      # A tibble: 2 x 9
        term        estimate std.error conf.level conf.low conf.high statistic
        <chr>          <dbl>     <dbl>      <dbl>    <dbl>     <dbl>     <dbl>
      1 (Intercept)    6.05     0.309        0.95    5.42      6.68      19.6 
      2 mpg           -0.141    0.0147       0.95   -0.171    -0.111     -9.56
        df.error  p.value
           <int>    <dbl>
      1       30 1.20e-18
      2       30 1.29e-10

---

    Code
      df_t$label
    Output
      list(widehat(italic(beta))=='6.05', italic(t)('30')=='19.59', italic(p)=='1.20e-18')
      list(widehat(italic(beta))=='-0.14', italic(t)('30')=='-9.56', italic(p)=='1.29e-10')

---

    Code
      df_t_na$label
    Output
      list(widehat(italic(beta))=='6.05', italic(t)=='19.59', italic(p)=='1.20e-18')
      list(widehat(italic(beta))=='-0.14', italic(t)=='-9.56', italic(p)=='1.29e-10')

---

    Code
      df_t_inf$label
    Output
      list(widehat(italic(beta))=='6.05', italic(t)=='19.59', italic(p)=='1.20e-18')
      list(widehat(italic(beta))=='-0.14', italic(t)=='-9.56', italic(p)=='1.29e-10')

---

    Code
      select(df_chi, -label)
    Output
      # A tibble: 2 x 9
        term  estimate std.error conf.level conf.low conf.high statistic df.error
        <chr>    <dbl>     <dbl>      <dbl>    <dbl>     <dbl>     <dbl>    <dbl>
      1 age     0.0170   0.00923       0.95 -0.00106    0.0351      3.40        1
      2 sex    -0.512    0.168         0.95 -0.840     -0.183       9.31        1
        p.value
          <dbl>
      1 0.0650 
      2 0.00228

---

    Code
      df_chi$label
    Output
      list(widehat(italic(beta))=='0.02', italic(chi)^2*('1')=='3.40', italic(p)=='0.07')
      list(widehat(italic(beta))=='-0.51', italic(chi)^2*('1')=='9.31', italic(p)=='2.28e-03')

---

    Code
      select(df_z, -label)
    Output
      # A tibble: 3 x 9
        term        estimate std.error conf.level conf.low conf.high statistic
        <chr>          <dbl>     <dbl>      <dbl>    <dbl>     <dbl>     <dbl>
      1 (Intercept)   -0.780     0.225       0.95    -1.22    -0.342     -3.47
      2 SexFemale      2.29      0.120       0.95     2.06     2.53      19.1 
      3 AgeAdult      -0.556     0.228       0.95    -1.00    -0.108     -2.44
        df.error  p.value
           <dbl>    <dbl>
      1      Inf 5.14e- 4
      2      Inf 1.54e-81
      3      Inf 1.45e- 2

---

    Code
      df_z$label
    Output
      list(widehat(italic(beta))=='-0.78', italic(z)=='-3.47', italic(p)=='5.14e-04')
      list(widehat(italic(beta))=='2.29', italic(z)=='19.13', italic(p)=='1.54e-81')
      list(widehat(italic(beta))=='-0.56', italic(z)=='-2.44', italic(p)=='0.01')

#  parametric t-tests

    Code
      df_1
    Output
      # A tibble: 4 x 15
           mu statistic df.error p.value method            alternative effectsize
        <dbl>     <dbl>    <dbl>   <dbl> <chr>             <chr>       <chr>     
      1  0.25     0.242       55   0.810 One Sample t-test two.sided   Cohen's d 
      2  0.25     0.242       55   0.595 One Sample t-test less        Hedges' g 
      3  0.25     0.242       55   0.405 One Sample t-test greater     Cohen's d 
      4  0.25     0.242       55   0.810 One Sample t-test two.sided   Hedges' g 
        estimate conf.level conf.low conf.high conf.method conf.distribution n.obs
           <dbl>      <dbl>    <dbl>     <dbl> <chr>       <chr>             <int>
      1   0.0323       0.89  -0.183      0.248 ncp         t                    56
      2   0.0319       0.99  -0.311      0.375 ncp         t                    56
      3   0.0323       0.9   -0.189      0.254 ncp         t                    56
      4   0.0319       0.5   -0.0577     0.122 ncp         t                    56
        expression  
        <list>      
      1 <expression>
      2 <expression>
      3 <expression>
      4 <expression>

---

    Code
      df_2_between
    Output
      # A tibble: 4 x 18
        term  group mean.group1 mean.group2 statistic df.error    p.value
        <chr> <chr>       <dbl>       <dbl>     <dbl>    <dbl>      <dbl>
      1 wt    am           3.77        2.41      5.26     30   0.0000113 
      2 wt    am           3.77        2.41      5.49     29.2 1.00      
      3 wt    am           3.77        2.41      5.26     30   0.00000563
      4 wt    am           3.77        2.41      5.49     29.2 0.00000627
        method                  alternative effectsize estimate conf.level conf.low
        <chr>                   <chr>       <chr>         <dbl>      <dbl>    <dbl>
      1 Two Sample t-test       two.sided   Cohen's d      1.93       0.89    1.23 
      2 Welch Two Sample t-test less        Hedges' g      1.88       0.99    0.793
      3 Two Sample t-test       greater     Cohen's d      1.93       0.9     1.21 
      4 Welch Two Sample t-test two.sided   Hedges' g      1.88       0.5     1.58 
        conf.high conf.method conf.distribution n.obs expression  
            <dbl> <chr>       <chr>             <int> <list>      
      1      2.61 ncp         t                    32 <expression>
      2      2.97 ncp         t                    32 <expression>
      3      2.63 ncp         t                    32 <expression>
      4      2.15 ncp         t                    32 <expression>

---

    Code
      df_2_within
    Output
      # A tibble: 4 x 16
        term   group     statistic df.error  p.value method        alternative
        <chr>  <chr>         <dbl>    <dbl>    <dbl> <chr>         <chr>      
      1 desire condition      3.61       89 0.000500 Paired t-test two.sided  
      2 desire condition      3.61       89 0.000500 Paired t-test two.sided  
      3 desire condition      3.61       89 0.000500 Paired t-test two.sided  
      4 desire condition      3.61       89 0.000500 Paired t-test two.sided  
        effectsize estimate conf.level conf.low conf.high conf.method
        <chr>         <dbl>      <dbl>    <dbl>     <dbl> <chr>      
      1 Cohen's d     0.381       0.89   0.206      0.557 ncp        
      2 Hedges' g     0.378       0.99   0.0984     0.659 ncp        
      3 Cohen's d     0.381       0.9    0.201      0.563 ncp        
      4 Hedges' g     0.378       0.5    0.305      0.452 ncp        
        conf.distribution n.obs expression  
        <chr>             <int> <list>      
      1 t                    90 <expression>
      2 t                    90 <expression>
      3 t                    90 <expression>
      4 t                    90 <expression>

---

    Code
      df_3_between
    Output
      # A tibble: 4 x 14
        statistic    df df.error p.value
            <dbl> <dbl>    <dbl>   <dbl>
      1      4.14     3     52    0.0105
      2      2.63     3     11.1  0.102 
      3      4.14     3     52    0.0105
      4      2.63     3     11.1  0.102 
        method                                                   effectsize estimate
        <chr>                                                    <chr>         <dbl>
      1 One-way analysis of means                                Eta2          0.193
      2 One-way analysis of means (not assuming equal variances) Omega2        0.245
      3 One-way analysis of means                                Eta2          0.193
      4 One-way analysis of means (not assuming equal variances) Omega2        0.245
        conf.level conf.low conf.high conf.method conf.distribution n.obs expression  
             <dbl>    <dbl>     <dbl> <chr>       <chr>             <int> <list>      
      1       0.89   0.0585         1 ncp         F                    56 <expression>
      2       0.8    0              1 ncp         F                    56 <expression>
      3       0.9    0.0545         1 ncp         F                    56 <expression>
      4       0.5    0.0974         1 ncp         F                    56 <expression>

---

    Code
      df_3_within
    Output
      # A tibble: 2 x 18
        term      sumsq sum.squares.error    df df.error meansq statistic  p.value
        <chr>     <dbl>             <dbl> <dbl>    <dbl>  <dbl>     <dbl>    <dbl>
      1 condition  233.              984.  2.63     229.   4.30      20.6 8.27e-11
      2 condition  233.              984.  2.63     229.   4.30      20.6 8.27e-11
        method                                              effectsize       estimate
        <chr>                                               <chr>               <dbl>
      1 ANOVA estimation for factorial designs using 'afex' Eta2 (partial)     0.191 
      2 ANOVA estimation for factorial designs using 'afex' Omega2 (partial)   0.0783
        conf.level conf.low conf.high conf.method conf.distribution n.obs expression  
             <dbl>    <dbl>     <dbl> <chr>       <chr>             <int> <list>      
      1       0.89   0.136          1 ncp         F                    88 <expression>
      2       0.9    0.0362         1 ncp         F                    88 <expression>

# corr_test works - nonparametric

    Code
      select(df1, -expression)
    Output
      # A tibble: 1 x 11
        parameter1 parameter2 effectsize           estimate conf.level conf.low
        <chr>      <chr>      <chr>                   <dbl>      <dbl>    <dbl>
      1 rating     length     Spearman correlation    0.495      0.999    0.153
        conf.high statistic    p.value method               n.obs
            <dbl>     <dbl>      <dbl> <chr>                <int>
      1     0.731    41453. 0.00000344 Spearman correlation    79

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(italic(\"S\") == \"41452.97684\", italic(p) == \"3.44384e-06\", widehat(rho)[\"Spearman\"] == \"0.49546\", CI[\"99.9%\"] ~ \"[\" * \"0.15344\", \"0.73147\" * \"]\", italic(\"n\")[\"pairs\"] == \"79\")"

---

    Code
      select(df2, -expression)
    Output
      # A tibble: 1 x 11
        parameter1 parameter2 effectsize           estimate conf.level conf.low
        <chr>      <chr>      <chr>                   <dbl>      <dbl>    <dbl>
      1 wt         mpg        Spearman correlation   -0.886       0.95   -0.945
        conf.high statistic  p.value method               n.obs
            <dbl>     <dbl>    <dbl> <chr>                <int>
      1    -0.774    10292. 1.49e-11 Spearman correlation    32

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(italic(\"S\") == \"10292.32\", italic(p) == \"1.49e-11\", widehat(rho)[\"Spearman\"] == \"-0.89\", CI[\"95%\"] ~ \"[\" * \"-0.94\", \"-0.77\" * \"]\", italic(\"n\")[\"pairs\"] == \"32\")"

# corr_test works - parametric

    Code
      select(df, -expression)
    Output
      # A tibble: 1 x 12
        parameter1 parameter2 effectsize          estimate conf.level conf.low
        <chr>      <chr>      <chr>                  <dbl>      <dbl>    <dbl>
      1 brainwt    sleep_rem  Pearson correlation   -0.221        0.9   -0.438
        conf.high statistic df.error p.value method              n.obs
            <dbl>     <dbl>    <int>   <dbl> <chr>               <int>
      1    0.0201     -1.54       46   0.131 Pearson correlation    48

---

    Code
      as.character(df$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Student\"] * \"(\" * 46 * \")\" == \"-1.539\", italic(p) == \"0.131\", widehat(italic(\"r\"))[\"Pearson\"] == \"-0.221\", CI[\"90%\"] ~ \"[\" * \"-0.438\", \"0.020\" * \"]\", italic(\"n\")[\"pairs\"] == \"48\")"

# corr_test works - robust

    Code
      select(df, -expression)
    Output
      # A tibble: 1 x 12
        parameter1 parameter2  effectsize                     estimate conf.level
        <chr>      <chr>       <chr>                             <dbl>      <dbl>
      1 brainwt    sleep_total Winsorized Pearson correlation   -0.549        0.5
        conf.low conf.high statistic df.error   p.value method                        
           <dbl>     <dbl>     <dbl>    <int>     <dbl> <chr>                         
      1   -0.611    -0.481     -4.83       54 0.0000117 Winsorized Pearson correlation
        n.obs
        <int>
      1    56

---

    Code
      as.character(df$expression[[1]])
    Output
      [1] "list(italic(\"t\")[\"Student\"] * \"(\" * 54 * \")\" == \"-4.8286\", italic(p) == \"1.1723e-05\", widehat(italic(\"r\"))[\"Winsorized\"] == \"-0.5491\", CI[\"50%\"] ~ \"[\" * \"-0.6106\", \"-0.4812\" * \"]\", italic(\"n\")[\"pairs\"] == \"56\")"

# bayes factor (correlation test) - without NAs

    Code
      as.character(subtitle1$expression[[1]])
    Output
      [1] "list(atop(\"huh is duh\", list(log[e] * (BF[\"01\"]) == \"1.07\", widehat(rho)[\"Pearson\"]^\"posterior\" == \"-0.12\", CI[\"95%\"]^HDI ~ \"[\" * \"-0.28\", \"0.04\" * \"]\", italic(\"r\")[\"beta\"]^\"JZS\" == \"1.41\")))"

# bayes factor (correlation test) - with NAs

    Code
      as.character(subtitle1$expression[[1]])
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"0.49\", widehat(rho)[\"Pearson\"]^\"posterior\" == \"-0.21\", CI[\"99%\"]^HDI ~ \"[\" * \"-0.47\", \"0.05\" * \"]\", italic(\"r\")[\"beta\"]^\"JZS\" == \"1.25\")"

# bayes factor (between-subjects - anova)

    Code
      dplyr::select(df1, -expression)
    Output
      # A tibble: 7 x 16
        term            pd rope.percentage prior.distribution prior.location
        <chr>        <dbl>           <dbl> <chr>                       <dbl>
      1 mu           0.940          0.168  cauchy                          0
      2 vore-carni   0.679          0.293  cauchy                          0
      3 vore-herbi   0.944          0.108  cauchy                          0
      4 vore-insecti 0.688          0.273  cauchy                          0
      5 vore-omni    0.646          0.349  cauchy                          0
      6 sig2         1              0      cauchy                          0
      7 g_vore       1              0.0174 cauchy                          0
        prior.scale  bf10 method                          log_e_bf10
              <dbl> <dbl> <chr>                                <dbl>
      1        0.99 0.118 Bayes factors for linear models      -2.14
      2        0.99 0.118 Bayes factors for linear models      -2.14
      3        0.99 0.118 Bayes factors for linear models      -2.14
      4        0.99 0.118 Bayes factors for linear models      -2.14
      5        0.99 0.118 Bayes factors for linear models      -2.14
      6        0.99 0.118 Bayes factors for linear models      -2.14
      7        0.99 0.118 Bayes factors for linear models      -2.14
        effectsize         estimate std.dev conf.level conf.low conf.high n.obs
        <chr>                 <dbl>   <dbl>      <dbl>    <dbl>     <dbl> <int>
      1 Bayesian R-squared        0       0       0.95        0    0.0800    51
      2 Bayesian R-squared        0       0       0.95        0    0.0800    51
      3 Bayesian R-squared        0       0       0.95        0    0.0800    51
      4 Bayesian R-squared        0       0       0.95        0    0.0800    51
      5 Bayesian R-squared        0       0       0.95        0    0.0800    51
      6 Bayesian R-squared        0       0       0.95        0    0.0800    51
      7 Bayesian R-squared        0       0       0.95        0    0.0800    51

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"2.139\", widehat(italic(R^\"2\"))[\"Bayesian\"]^\"posterior\" == \"0.000\", CI[\"95%\"]^HDI ~ \"[\" * \"0.000\", \"0.080\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.990\")"

---

    Code
      dplyr::select(df2, -expression)
    Output
      # A tibble: 6 x 16
        term                  pd rope.percentage prior.distribution prior.location
        <chr>              <dbl>           <dbl> <chr>                       <dbl>
      1 mu                 1               0     cauchy                          0
      2 Species-setosa     1               0     cauchy                          0
      3 Species-versicolor 0.936           0.435 cauchy                          0
      4 Species-virginica  1               0     cauchy                          0
      5 sig2               1               0     cauchy                          0
      6 g_Species          1               0     cauchy                          0
        prior.scale    bf10 method                          log_e_bf10
              <dbl>   <dbl> <chr>                                <dbl>
      1       0.707 1.87e28 Bayes factors for linear models       65.1
      2       0.707 1.87e28 Bayes factors for linear models       65.1
      3       0.707 1.87e28 Bayes factors for linear models       65.1
      4       0.707 1.87e28 Bayes factors for linear models       65.1
      5       0.707 1.87e28 Bayes factors for linear models       65.1
      6       0.707 1.87e28 Bayes factors for linear models       65.1
        effectsize         estimate std.dev conf.level conf.low conf.high n.obs
        <chr>                 <dbl>   <dbl>      <dbl>    <dbl>     <dbl> <int>
      1 Bayesian R-squared    0.612  0.0311       0.99    0.511     0.679   150
      2 Bayesian R-squared    0.612  0.0311       0.99    0.511     0.679   150
      3 Bayesian R-squared    0.612  0.0311       0.99    0.511     0.679   150
      4 Bayesian R-squared    0.612  0.0311       0.99    0.511     0.679   150
      5 Bayesian R-squared    0.612  0.0311       0.99    0.511     0.679   150
      6 Bayesian R-squared    0.612  0.0311       0.99    0.511     0.679   150

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"-65.0969\", widehat(italic(R^\"2\"))[\"Bayesian\"]^\"posterior\" == \"0.6118\", CI[\"99%\"]^HDI ~ \"[\" * \"0.5107\", \"0.6789\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.7070\")"

# bayes factor (within-subjects - anova)

    Code
      dplyr::select(df1, -expression)
    Output
      # A tibble: 7 x 18
        term           pd rope.percentage prior.distribution prior.location
        <chr>       <dbl>           <dbl> <chr>                       <dbl>
      1 mu          1              0      cauchy                          0
      2 Wine-Wine A 0.97           0.473  cauchy                          0
      3 Wine-Wine B 0.906          0.688  cauchy                          0
      4 Wine-Wine C 0.998          0.0755 cauchy                          0
      5 sig2        1              1      cauchy                          0
      6 g_Wine      1              0      cauchy                          0
      7 g_rowid     1              0      cauchy                          0
        prior.scale effect  bf10 method                          log_e_bf10
              <dbl> <chr>  <dbl> <chr>                                <dbl>
      1        0.88 fixed   7.09 Bayes factors for linear models       1.96
      2        0.88 fixed   7.09 Bayes factors for linear models       1.96
      3        0.88 fixed   7.09 Bayes factors for linear models       1.96
      4        0.88 fixed   7.09 Bayes factors for linear models       1.96
      5        1    fixed   7.09 Bayes factors for linear models       1.96
      6        1    fixed   7.09 Bayes factors for linear models       1.96
      7        1    fixed   7.09 Bayes factors for linear models       1.96
        effectsize         estimate std.dev conf.level conf.low conf.high component  
        <chr>                 <dbl>   <dbl>      <dbl>    <dbl>     <dbl> <chr>      
      1 Bayesian R-squared    0.893  0.0176       0.95    0.847     0.920 conditional
      2 Bayesian R-squared    0.893  0.0176       0.95    0.847     0.920 conditional
      3 Bayesian R-squared    0.893  0.0176       0.95    0.847     0.920 conditional
      4 Bayesian R-squared    0.893  0.0176       0.95    0.847     0.920 conditional
      5 Bayesian R-squared    0.893  0.0176       0.95    0.847     0.920 conditional
      6 Bayesian R-squared    0.893  0.0176       0.95    0.847     0.920 conditional
      7 Bayesian R-squared    0.893  0.0176       0.95    0.847     0.920 conditional
        n.obs
        <int>
      1    22
      2    22
      3    22
      4    22
      5    22
      6    22
      7    22

---

    Code
      as.character(df1$expression[[1]])
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"-1.96\", widehat(italic(R^\"2\"))[\"Bayesian\"]^\"posterior\" == \"0.89\", CI[\"95%\"]^HDI ~ \"[\" * \"0.85\", \"0.92\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.88\")"

---

    Code
      dplyr::select(df2, -expression)
    Output
      # A tibble: 8 x 18
        term              pd rope.percentage prior.distribution prior.location
        <chr>          <dbl>           <dbl> <chr>                       <dbl>
      1 mu             1               0     cauchy                          0
      2 condition-HDHF 1               0     cauchy                          0
      3 condition-HDLF 0.862           0.715 cauchy                          0
      4 condition-LDHF 0.995           0.139 cauchy                          0
      5 condition-LDLF 1               0     cauchy                          0
      6 sig2           1               0     cauchy                          0
      7 g_condition    1               0.401 cauchy                          0
      8 g_rowid        1               0     cauchy                          0
        prior.scale effect        bf10 method                          log_e_bf10
              <dbl> <chr>        <dbl> <chr>                                <dbl>
      1       0.707 fixed  1372773375. Bayes factors for linear models       21.0
      2       0.707 fixed  1372773375. Bayes factors for linear models       21.0
      3       0.707 fixed  1372773375. Bayes factors for linear models       21.0
      4       0.707 fixed  1372773375. Bayes factors for linear models       21.0
      5       0.707 fixed  1372773375. Bayes factors for linear models       21.0
      6       1     fixed  1372773375. Bayes factors for linear models       21.0
      7       1     fixed  1372773375. Bayes factors for linear models       21.0
      8       1     fixed  1372773375. Bayes factors for linear models       21.0
        effectsize         estimate std.dev conf.level conf.low conf.high component  
        <chr>                 <dbl>   <dbl>      <dbl>    <dbl>     <dbl> <chr>      
      1 Bayesian R-squared    0.529  0.0330       0.95    0.460     0.586 conditional
      2 Bayesian R-squared    0.529  0.0330       0.95    0.460     0.586 conditional
      3 Bayesian R-squared    0.529  0.0330       0.95    0.460     0.586 conditional
      4 Bayesian R-squared    0.529  0.0330       0.95    0.460     0.586 conditional
      5 Bayesian R-squared    0.529  0.0330       0.95    0.460     0.586 conditional
      6 Bayesian R-squared    0.529  0.0330       0.95    0.460     0.586 conditional
      7 Bayesian R-squared    0.529  0.0330       0.95    0.460     0.586 conditional
      8 Bayesian R-squared    0.529  0.0330       0.95    0.460     0.586 conditional
        n.obs
        <int>
      1    88
      2    88
      3    88
      4    88
      5    88
      6    88
      7    88
      8    88

---

    Code
      as.character(df2$expression[[1]])
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"-21.04\", widehat(italic(R^\"2\"))[\"Bayesian\"]^\"posterior\" == \"0.53\", CI[\"95%\"]^HDI ~ \"[\" * \"0.46\", \"0.59\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.71\")"

# Getting help with statsExpressions

Thanks for using statsExpressions. Before filing an issue, there are a few places
to explore and pieces to put together to make the process as smooth as possible.

Start by making a minimal **repr**oducible **ex**ample using the 
[reprex](https://reprex.tidyverse.org/) package. If you haven't heard of or used 
reprex before, you're in for a treat! Seriously, reprex will make all of your 
R-question-asking endeavors easier (which is a pretty insane ROI for the five to 
ten minutes it'll take you to learn what it's all about). For additional reprex
pointers, check out the [Get help!](https://www.tidyverse.org/help/) section of
the tidyverse site.

Armed with your reprex, the next step is to figure out [where to ask](https://www.tidyverse.org/help/#where-to-ask). 

  * If it's a question: start with [community.rstudio.com](https://community.rstudio.com/), 
    and/or StackOverflow. There are more people there to answer questions.  
  * If it's a bug: you're in the right place, file an issue.  
  * If you're not sure: let the community help you figure it out! If your 
    problem _is_ a bug or a feature request, you can easily return here and 
    report it. 

Before opening a new issue, be sure to [search issues and pull requests](https://github.com/tidyverse/statsExpressions/issues) to make sure the 
bug hasn't been reported and/or already fixed in the development version. By 
default, the search will be pre-populated with `is:issue is:open`. You can 
[edit the qualifiers](https://help.github.com/articles/searching-issues-and-pull-requests/) 
(e.g. `is:pr`, `is:closed`) as needed. For example, you'd simply
remove `is:open` to search _all_ issues in the repo, open or closed.


If you _are_ in the right place, and need to file an issue, please review the 
["File issues"](https://www.tidyverse.org/contribute/#issues) paragraph from 
the tidyverse contributing guidelines.

Thanks for your help!
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
reported to the community leaders responsible for enforcement at patilindrajeet.science@gmail.com. 
All complaints will be reviewed and investigated promptly and fairly.

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
available at <https://www.contributor-covenant.org/version/2/0/code_of_conduct.html>.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
<https://www.contributor-covenant.org/faq>. Translations are available at <https://www.contributor-covenant.org/translations>.
# Contributing to statsExpressions

This outlines how to propose a change to statsExpressions. For more detailed
info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib).

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  New code should follow the tidyverse [style guide](https://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the statsExpressions project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See tidyverse [development contributing guide](https://rstd.io/tidy-contrib)
for further details.
---
name: Bug report or feature request
about: Describe a bug you've seen or make a case for a new feature
---

Please briefly describe your problem and what output you expect. If you have a question, please don't use this form. Instead, ask on <https://stackoverflow.com/> or <https://community.rstudio.com/>.

Please include a minimal reproducible example (AKA a reprex). If you've never heard of a [reprex](http://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>.

Brief description of the problem

```r
# insert reprex here
```
