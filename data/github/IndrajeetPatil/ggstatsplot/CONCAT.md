
<!-- README.md is generated from README.Rmd. Please edit that file -->

## `{ggstatsplot}`: `{ggplot2}` Based Plots with Statistical Details

| Status                                                                                                                                            | Usage                                                                                                                                            | Miscellaneous                                                                                                                                                    |
|---------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [![R build status](https://github.com/IndrajeetPatil/ggstatsplot/workflows/R-CMD-check/badge.svg)](https://github.com/IndrajeetPatil/ggstatsplot) | [![Total downloads badge](https://cranlogs.r-pkg.org/badges/grand-total/ggstatsplot?color=blue)](https://CRAN.R-project.org/package=ggstatsplot) | [![Codecov](https://codecov.io/gh/IndrajeetPatil/ggstatsplot/branch/master/graph/badge.svg)](https://app.codecov.io/gh/IndrajeetPatil/ggstatsplot?branch=master) |
| [![lints](https://github.com/IndrajeetPatil/ggstatsplot/workflows/lint/badge.svg)](https://github.com/IndrajeetPatil/ggstatsplot)                 | [![Daily downloads badge](https://cranlogs.r-pkg.org/badges/last-day/ggstatsplot?color=blue)](https://CRAN.R-project.org/package=ggstatsplot)    | [![status](https://tinyverse.netlify.com/badge/ggstatsplot)](https://CRAN.R-project.org/package=ggstatsplot)                                                     |
| [![pkgdown](https://github.com/IndrajeetPatil/ggstatsplot/workflows/pkgdown/badge.svg)](https://github.com/IndrajeetPatil/ggstatsplot/actions)    | [![DOI](https://joss.theoj.org/papers/10.21105/joss.03167/status.svg)](https://doi.org/10.21105/joss.03167)                                      | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)                                       |

## Raison d’être <img src="man/figures/logo.png" align="right" width="360" />

> “What is to be sought in designs for the display of information is the
> clear portrayal of complexity. Not the complication of the simple;
> rather … the revelation of the complex.”  
> - Edward R. Tufte

[`{ggstatsplot}`](https://indrajeetpatil.github.io/ggstatsplot/) is an
extension of [`{ggplot2}`](https://github.com/tidyverse/ggplot2) package
for creating graphics with details from statistical tests included in
the information-rich plots themselves. In a typical exploratory data
analysis workflow, data visualization and statistical modeling are two
different phases: visualization informs modeling, and modeling in its
turn can suggest a different visualization method, and so on and so
forth. The central idea of `{ggstatsplot}` is simple: combine these two
phases into one in the form of graphics with statistical details, which
makes data exploration simpler and faster.

## Installation

| Type        | Source                                                                                                             | Command                                                 |
|-------------|--------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| Release     | [![CRAN Status](https://www.r-pkg.org/badges/version/ggstatsplot)](https://cran.r-project.org/package=ggstatsplot) | `install.packages("ggstatsplot")`                       |
| Development | [![Project Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/##active)      | `remotes::install_github("IndrajeetPatil/ggstatsplot")` |

Linux users may encounter some installation problems. In particular, the
`{ggstatsplot}` package depends on the `{PMCMRplus}` package.

    ERROR: dependencies ‘gmp’, ‘Rmpfr’ are not available for package ‘PMCMRplus’

This means that your operating system lacks `gmp` and `Rmpfr` libraries.

If you use `Ubuntu`, you can install these dependencies:

    sudo apt-get install libgmp3-dev
    sudo apt-get install libmpfr-dev

The following `README` file briefly describes the installation
procedure:
<https://CRAN.R-project.org/package=PMCMRplus/readme/README.html>

## Citation

If you want to cite this package in a scientific journal or in any other
context, run the following code in your `R` console:

``` r
citation("ggstatsplot")

  Patil, I. (2021). Visualizations with statistical details: The
  'ggstatsplot' approach. Journal of Open Source Software, 6(61), 3167,
  doi:10.21105/joss.03167

A BibTeX entry for LaTeX users is

  @Article{,
    doi = {10.21105/joss.03167},
    url = {https://doi.org/10.21105/joss.03167},
    year = {2021},
    publisher = {{The Open Journal}},
    volume = {6},
    number = {61},
    pages = {3167},
    author = {Indrajeet Patil},
    title = {{Visualizations with statistical details: The {'ggstatsplot'} approach}},
    journal = {{Journal of Open Source Software}},
  }
```

There is currently a publication in preparation corresponding to this
package and the citation will be updated once it’s published.

## Documentation and Examples

To see the detailed documentation for each function in the stable
**CRAN** version of the package, see:

-   Website: <https://indrajeetpatil.github.io/ggstatsplot/>

-   Presentation: <a
    href="https://indrajeetpatil.github.io/ggstatsplot_slides/slides/ggstatsplot_presentation.html##1"
    class="uri">https://indrajeetpatil.github.io/ggstatsplot_slides/slides/ggstatsplot_presentation.html##1</a>

-   Vignettes: <https://indrajeetpatil.github.io/ggstatsplot/articles/>

## Summary of available plots

It, therefore, produces a limited kinds of plots for the supported
analyses:

| Function         | Plot                      | Description                                     | Lifecycle                                                                                                                  |
|------------------|---------------------------|-------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| `ggbetweenstats` | **violin plots**          | for comparisons *between* groups/conditions     | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| `ggwithinstats`  | **violin plots**          | for comparisons *within* groups/conditions      | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| `gghistostats`   | **histograms**            | for distribution about numeric variable         | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| `ggdotplotstats` | **dot plots/charts**      | for distribution about labeled numeric variable | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| `ggscatterstats` | **scatterplots**          | for correlation between two variables           | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| `ggcorrmat`      | **correlation matrices**  | for correlations between multiple variables     | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| `ggpiestats`     | **pie charts**            | for categorical data                            | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| `ggbarstats`     | **bar charts**            | for categorical data                            | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) |
| `ggcoefstats`    | **dot-and-whisker plots** | for regression models and meta-analysis         | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) |

In addition to these basic plots, `{ggstatsplot}` also provides
**`grouped_`** versions (see below) that makes it easy to repeat the
same analysis for any grouping variable.

## Summary of types of statistical analyses

The table below summarizes all the different types of analyses currently
supported in this package-

| Functions                        | Description                                       | Parametric | Non-parametric | Robust | Bayesian |
|----------------------------------|---------------------------------------------------|------------|----------------|--------|----------|
| `ggbetweenstats`                 | Between group/condition comparisons               | ✅         | ✅             | ✅     | ✅       |
| `ggwithinstats`                  | Within group/condition comparisons                | ✅         | ✅             | ✅     | ✅       |
| `gghistostats`, `ggdotplotstats` | Distribution of a numeric variable                | ✅         | ✅             | ✅     | ✅       |
| `ggcorrmat`                      | Correlation matrix                                | ✅         | ✅             | ✅     | ✅       |
| `ggscatterstats`                 | Correlation between two variables                 | ✅         | ✅             | ✅     | ✅       |
| `ggpiestats`, `ggbarstats`       | Association between categorical variables         | ✅         | ✅             | ❌     | ✅       |
| `ggpiestats`, `ggbarstats`       | Equal proportions for categorical variable levels | ✅         | ✅             | ❌     | ✅       |
| `ggcoefstats`                    | Regression model coefficients                     | ✅         | ✅             | ✅     | ✅       |
| `ggcoefstats`                    | Random-effects meta-analysis                      | ✅         | ❌             | ✅     | ✅       |

Summary of Bayesian analysis

| Analysis                        | Hypothesis testing | Estimation |
|---------------------------------|--------------------|------------|
| (one/two-sample) t-test         | ✅                 | ✅         |
| one-way ANOVA                   | ✅                 | ✅         |
| correlation                     | ✅                 | ✅         |
| (one/two-way) contingency table | ✅                 | ✅         |
| random-effects meta-analysis    | ✅                 | ✅         |

## Statistical reporting

For **all** statistical tests reported in the plots, the default
template abides by the gold standard for statistical reporting. For
example, here are results from Yuen’s test for trimmed means (robust
*t*-test):

<img src="man/figures/stats_reporting_format.png" align="center" />

## Summary of statistical tests and effect sizes

Statistical analysis is carried out by `{statsExpressions}` package, and
thus a summary table of all the statistical tests currently supported
across various functions can be found in article for that package:
<https://indrajeetpatil.github.io/statsExpressions/articles/stats_details.html>

## Primary functions

Here are examples of the main functions currently supported in
`{ggstatsplot}`.

**Note**: If you are reading this on `GitHub` repository, the
documentation below is for the **development** version of the package.
So you may see some features available here that are not currently
present in the stable version of this package on **CRAN**. For
documentation relevant for the `CRAN` version, see:
<https://CRAN.R-project.org/package=ggstatsplot/readme/README.html>

### `ggbetweenstats`

This function creates either a violin plot, a box plot, or a mix of two
for **between**-group or **between**-condition comparisons with results
from statistical tests in the subtitle. The simplest function call looks
like this-

``` r
## for reproducibility
set.seed(123)
library(ggstatsplot)

## plot
ggbetweenstats(
  data  = iris,
  x     = Species,
  y     = Sepal.Length,
  title = "Distribution of sepal length across Iris species"
)
```

<img src="man/figures/README-ggbetweenstats1-1.png" width="100%" />

**Defaults** return<br>

✅ raw data + distributions <br> ✅ descriptive statistics <br> ✅
inferential statistics <br> ✅ effect size + CIs <br> ✅ pairwise
comparisons <br> ✅ Bayesian hypothesis-testing <br> ✅ Bayesian
estimation <br>

A number of other arguments can be specified to make this plot even more
informative or change some of the default options. Additionally, there
is also a `grouped_` variant of this function that makes it easy to
repeat the same operation across a **single** grouping variable:

``` r
## for reproducibility
set.seed(123)

## plot
grouped_ggbetweenstats(
  data             = dplyr::filter(movies_long, genre %in% c("Action", "Comedy")),
  x                = mpaa,
  y                = length,
  grouping.var     = genre, ## grouping variable
  outlier.tagging  = TRUE, ## whether outliers need to be tagged
  outlier.label    = title, ## variable to be used for tagging outliers
  outlier.coef     = 2,
  ggsignif.args    = list(textsize = 4, tip_length = 0.01),
  p.adjust.method  = "bonferroni", ## method for adjusting p-values for multiple comparisons
  ## adding new components to `{ggstatsplot}` default
  ggplot.component = list(ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis())),
  caption          = "Source: IMDb (Internet Movie Database)",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  annotation.args  = list(title = "Differences in movie length by mpaa ratings for different genres")
)
```

<img src="man/figures/README-ggbetweenstats2-1.png" width="100%" />

Note here that the function can be used to tag outliers!

##### Summary of graphics

| graphical element        | `geom_` used                | argument for further modification |
|--------------------------|-----------------------------|-----------------------------------|
| raw data                 | `ggplot2::geom_point`       | `point.args`                      |
| box plot                 | `ggplot2::geom_boxplot`     | ❌                                |
| density plot             | `ggplot2::geom_violin`      | `violin.args`                     |
| centrality measure point | `ggplot2::geom_point`       | `centrality.point.args`           |
| centrality measure label | `ggrepel::geom_label_repel` | `centrality.label.args`           |
| outlier point            | `ggplot2::stat_boxplot`     | ❌                                |
| outlier label            | `ggrepel::geom_label_repel` | `outlier.label.args`              |
| pairwise comparisons     | `ggsignif::geom_signif`     | `ggsignif.args`                   |

##### Summary of tests

**Central tendency measure**

| Type           | Measure                                           | Function used                       |
|----------------|---------------------------------------------------|-------------------------------------|
| Parametric     | mean                                              | `parameters::describe_distribution` |
| Non-parametric | median                                            | `parameters::describe_distribution` |
| Robust         | trimmed mean                                      | `parameters::describe_distribution` |
| Bayesian       | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution` |

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

| Type           | No. of groups | Effect size                                                                                                                                                                                            | CI? | Function used                                          |
|----------------|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|--------------------------------------------------------|
| Parametric     | \> 2          | ![\\eta\_{p}^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Ceta_%7Bp%7D%5E2 "\eta_{p}^2"), ![\\omega\_{p}^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Comega_%7Bp%7D%5E2 "\omega_{p}^2") | ✅  | `effectsize::omega_squared`, `effectsize::eta_squared` |
| Non-parametric | \> 2          | ![\\epsilon\_{ordinal}^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cepsilon_%7Bordinal%7D%5E2 "\epsilon_{ordinal}^2")                                                                          | ✅  | `effectsize::rank_epsilon_squared`                     |
| Robust         | \> 2          | ![\\xi](https://chart.apis.google.com/chart?cht=tx&chl=%5Cxi "\xi") (Explanatory measure of effect size)                                                                                               | ✅  | `WRS2::t1way`                                          |
| Bayes Factor   | \> 2          | ![R\_{posterior}^2](https://chart.apis.google.com/chart?cht=tx&chl=R_%7Bposterior%7D%5E2 "R_{posterior}^2")                                                                                            | ✅  | `performance::r2_bayes`                                |
| Parametric     | 2             | Cohen’s *d*, Hedge’s *g*                                                                                                                                                                               | ✅  | `effectsize::cohens_d`, `effectsize::hedges_g`         |
| Non-parametric | 2             | *r* (rank-biserial correlation)                                                                                                                                                                        | ✅  | `effectsize::rank_biserial`                            |
| Robust         | 2             | ![\\xi](https://chart.apis.google.com/chart?cht=tx&chl=%5Cxi "\xi") (Explanatory measure of effect size)                                                                                               | ✅  | `WRS2::yuen.effect.ci`                                 |
| Bayesian       | 2             | ![\\delta\_{posterior}](https://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7Bposterior%7D "\delta_{posterior}")                                                                                  | ✅  | `bayestestR::describe_posterior`                       |

**Pairwise comparison tests**

| Type           | Equal variance? | Test                      | *p*-value adjustment? | Function used                   |
|----------------|-----------------|---------------------------|-----------------------|---------------------------------|
| Parametric     | No              | Games-Howell test         | ✅                    | `PMCMRplus::gamesHowellTest`    |
| Parametric     | Yes             | Student’s *t*-test        | ✅                    | `stats::pairwise.t.test`        |
| Non-parametric | No              | Dunn test                 | ✅                    | `PMCMRplus::kwAllPairsDunnTest` |
| Robust         | No              | Yuen’s trimmed means test | ✅                    | `WRS2::lincon`                  |
| Bayesian       | `NA`            | Student’s *t*-test        | `NA`                  | `BayesFactor::ttestBF`          |

For more, see the `ggbetweenstats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggbetweenstats.html>

### `ggwithinstats`

`ggbetweenstats` function has an identical twin function `ggwithinstats`
for repeated measures designs that behaves in the same fashion with a
few minor tweaks introduced to properly visualize the repeated measures
design. As can be seen from an example below, the only difference
between the plot structure is that now the group means are connected by
paths to highlight the fact that these data are paired with each other.

``` r
## for reproducibility and data
set.seed(123)
library(WRS2) ## for data
library(afex) ## to run anova

## plot
ggwithinstats(
  data    = WineTasting,
  x       = Wine,
  y       = Taste,
  title   = "Wine tasting",
  caption = "Data source: `WRS2` R package",
  ggtheme = ggthemes::theme_fivethirtyeight()
)
```

<img src="man/figures/README-ggwithinstats1-1.png" width="100%" />

**Defaults** return<br>

✅ raw data + distributions <br> ✅ descriptive statistics <br> ✅
inferential statistics <br> ✅ effect size + CIs <br> ✅ pairwise
comparisons <br> ✅ Bayesian hypothesis-testing <br> ✅ Bayesian
estimation <br>

The central tendency measure displayed will depend on the statistics:

| Type           | Measure      | Function used                       |
|----------------|--------------|-------------------------------------|
| Parametric     | mean         | `parameters::describe_distribution` |
| Non-parametric | median       | `parameters::describe_distribution` |
| Robust         | trimmed mean | `parameters::describe_distribution` |
| Bayesian       | MAP estimate | `parameters::describe_distribution` |

As with the `ggbetweenstats`, this function also has a `grouped_`
variant that makes repeating the same analysis across a single grouping
variable quicker. We will see an example with only repeated
measurements-

``` r
## common setup
set.seed(123)

## plot
grouped_ggwithinstats(
  data            = dplyr::filter(bugs_long, region %in% c("Europe", "North America"), condition %in% c("LDLF", "LDHF")),
  x               = condition,
  y               = desire,
  type            = "np", ## non-parametric statistics
  xlab            = "Condition",
  ylab            = "Desire to kill an artrhopod",
  grouping.var    = region,
  outlier.tagging = TRUE,
  outlier.label   = education
)
```

<img src="man/figures/README-ggwithinstats2-1.png" width="100%" />

##### Summary of graphics

| graphical element             | `geom_` used                | argument for further modification |
|-------------------------------|-----------------------------|-----------------------------------|
| raw data                      | `ggplot2::geom_point`       | `point.args`                      |
| point path                    | `ggplot2::geom_path`        | `point.path.args`                 |
| box plot                      | `ggplot2::geom_boxplot`     | `boxplot.args`                    |
| density plot                  | `ggplot2::geom_violin`      | `violin.args`                     |
| centrality measure point      | `ggplot2::geom_point`       | `centrality.point.args`           |
| centrality measure point path | `ggplot2::geom_path`        | `centrality.path.args`            |
| centrality measure label      | `ggrepel::geom_label_repel` | `centrality.label.args`           |
| outlier point                 | `ggplot2::stat_boxplot`     | ❌                                |
| outlier label                 | `ggrepel::geom_label_repel` | `outlier.label.args`              |
| pairwise comparisons          | `ggsignif::geom_signif`     | `ggsignif.args`                   |

##### Summary of tests

**Central tendency measure**

| Type           | Measure                                           | Function used                       |
|----------------|---------------------------------------------------|-------------------------------------|
| Parametric     | mean                                              | `parameters::describe_distribution` |
| Non-parametric | median                                            | `parameters::describe_distribution` |
| Robust         | trimmed mean                                      | `parameters::describe_distribution` |
| Bayesian       | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution` |

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

| Type           | No. of groups | Effect size                                                                                                                                                                                            | CI? | Function used                                          |
|----------------|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|--------------------------------------------------------|
| Parametric     | \> 2          | ![\\eta\_{p}^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Ceta_%7Bp%7D%5E2 "\eta_{p}^2"), ![\\omega\_{p}^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Comega_%7Bp%7D%5E2 "\omega_{p}^2") | ✅  | `effectsize::omega_squared`, `effectsize::eta_squared` |
| Non-parametric | \> 2          | ![W\_{Kendall}](https://chart.apis.google.com/chart?cht=tx&chl=W_%7BKendall%7D "W_{Kendall}") (Kendall’s coefficient of concordance)                                                                   | ✅  | `effectsize::kendalls_w`                               |
| Robust         | \> 2          | ![\\delta\_{R-avg}^{AKP}](https://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7BR-avg%7D%5E%7BAKP%7D "\delta_{R-avg}^{AKP}") (Algina-Keselman-Penfield robust standardized difference average)    | ✅  | `WRS2::wmcpAKP`                                        |
| Bayes Factor   | \> 2          | ![R\_{Bayesian}^2](https://chart.apis.google.com/chart?cht=tx&chl=R_%7BBayesian%7D%5E2 "R_{Bayesian}^2")                                                                                               | ✅  | `performance::r2_bayes`                                |
| Parametric     | 2             | Cohen’s *d*, Hedge’s *g*                                                                                                                                                                               | ✅  | `effectsize::cohens_d`, `effectsize::hedges_g`         |
| Non-parametric | 2             | *r* (rank-biserial correlation)                                                                                                                                                                        | ✅  | `effectsize::rank_biserial`                            |
| Robust         | 2             | ![\\delta\_{R}^{AKP}](https://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7BR%7D%5E%7BAKP%7D "\delta_{R}^{AKP}") (Algina-Keselman-Penfield robust standardized difference)                        | ✅  | `WRS2::wmcpAKP`                                        |
| Bayesian       | 2             | ![\\delta\_{posterior}](https://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7Bposterior%7D "\delta_{posterior}")                                                                                  | ✅  | `bayestestR::describe_posterior`                       |

**Pairwise comparison tests**

| Type           | Test                      | *p*-value adjustment? | Function used                   |
|----------------|---------------------------|-----------------------|---------------------------------|
| Parametric     | Student’s *t*-test        | ✅                    | `stats::pairwise.t.test`        |
| Non-parametric | Durbin-Conover test       | ✅                    | `PMCMRplus::durbinAllPairsTest` |
| Robust         | Yuen’s trimmed means test | ✅                    | `WRS2::rmmcp`                   |
| Bayesian       | Student’s *t*-test        | ❌                    | `BayesFactor::ttestBF`          |

For more, see the `ggwithinstats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggwithinstats.html>

### `gghistostats`

To visualize the distribution of a single variable and check if its mean
is significantly different from a specified value with a one-sample
test, `gghistostats` can be used.

``` r
## for reproducibility
set.seed(123)

## plot
gghistostats(
  data       = ggplot2::msleep, ## dataframe from which variable is to be taken
  x          = awake, ## numeric variable whose distribution is of interest
  title      = "Amount of time spent awake", ## title for the plot
  caption    = "Source: Mammalian sleep data set",
  test.value = 12, ## default value is 0
  binwidth   = 1, ## binwidth value (experiment)
  ggtheme    = hrbrthemes::theme_ipsum_tw()
)
```

<img src="man/figures/README-gghistostats1-1.png" width="100%" />

**Defaults** return<br>

✅ counts + proportion for bins<br> ✅ descriptive statistics <br> ✅
inferential statistics <br> ✅ effect size + CIs <br> ✅ Bayesian
hypothesis-testing <br> ✅ Bayesian estimation <br>

There is also a `grouped_` variant of this function that makes it easy
to repeat the same operation across a **single** grouping variable:

``` r
## for reproducibility
set.seed(123)

## plot
grouped_gghistostats(
  data              = dplyr::filter(movies_long, genre %in% c("Action", "Comedy")),
  x                 = budget,
  test.value        = 50,
  type              = "nonparametric",
  xlab              = "Movies budget (in million US$)",
  grouping.var      = genre, ## grouping variable
  normal.curve      = TRUE, ## superimpose a normal distribution curve
  normal.curve.args = list(color = "red", size = 1),
  ggtheme           = ggthemes::theme_tufte(),
  ## modify the defaults from `{ggstatsplot}` for each plot
  ggplot.component  = ggplot2::labs(caption = "Source: IMDB.com"),
  plotgrid.args     = list(nrow = 1),
  annotation.args   = list(title = "Movies budgets for different genres")
)
```

<img src="man/figures/README-gghistostats2-1.png" width="100%" />

##### Summary of graphics

| graphical element       | `geom_` used             | argument for further modification |
|-------------------------|--------------------------|-----------------------------------|
| histogram bin           | `ggplot2::stat_bin`      | `bin.args`                        |
| centrality measure line | `ggplot2::geom_vline`    | `centrality.line.args`            |
| normality curve         | `ggplot2::stat_function` | `normal.curve.args`               |

##### Summary of tests

**Central tendency measure**

| Type           | Measure                                           | Function used                       |
|----------------|---------------------------------------------------|-------------------------------------|
| Parametric     | mean                                              | `parameters::describe_distribution` |
| Non-parametric | median                                            | `parameters::describe_distribution` |
| Robust         | trimmed mean                                      | `parameters::describe_distribution` |
| Bayesian       | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution` |

**Hypothesis testing**

| Type           | Test                                     | Function used          |
|----------------|------------------------------------------|------------------------|
| Parametric     | One-sample Student’s *t*-test            | `stats::t.test`        |
| Non-parametric | One-sample Wilcoxon test                 | `stats::wilcox.test`   |
| Robust         | Bootstrap-*t* method for one-sample test | `WRS2::trimcibt`       |
| Bayesian       | One-sample Student’s *t*-test            | `BayesFactor::ttestBF` |

**Effect size estimation**

| Type           | Effect size                                                                                                           | CI? | Function used                                  |
|----------------|-----------------------------------------------------------------------------------------------------------------------|-----|------------------------------------------------|
| Parametric     | Cohen’s *d*, Hedge’s *g*                                                                                              | ✅  | `effectsize::cohens_d`, `effectsize::hedges_g` |
| Non-parametric | *r* (rank-biserial correlation)                                                                                       | ✅  | `effectsize::rank_biserial`                    |
| Robust         | trimmed mean                                                                                                          | ✅  | `WRS2::trimcibt`                               |
| Bayes Factor   | ![\\delta\_{posterior}](https://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7Bposterior%7D "\delta_{posterior}") | ✅  | `bayestestR::describe_posterior`               |

For more, including information about the variant of this function
`grouped_gghistostats`, see the `gghistostats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/gghistostats.html>

### `ggdotplotstats`

This function is similar to `gghistostats`, but is intended to be used
when the numeric variable also has a label.

``` r
## for reproducibility
set.seed(123)

## plot
ggdotplotstats(
  data       = dplyr::filter(gapminder::gapminder, continent == "Asia"),
  y          = country,
  x          = lifeExp,
  test.value = 55,
  type       = "robust",
  title      = "Distribution of life expectancy in Asian continent",
  xlab       = "Life expectancy",
  caption    = "Source: Gapminder dataset from https://www.gapminder.org/"
)
```

<img src="man/figures/README-ggdotplotstats1-1.png" width="100%" />

**Defaults** return<br>

✅ descriptives (mean + sample size) <br> ✅ inferential statistics <br>
✅ effect size + CIs <br> ✅ Bayesian hypothesis-testing <br> ✅
Bayesian estimation <br>

As with the rest of the functions in this package, there is also a
`grouped_` variant of this function to facilitate looping the same
operation for all levels of a single grouping variable.

``` r
## for reproducibility
set.seed(123)

## plot
grouped_ggdotplotstats(
  data            = dplyr::filter(ggplot2::mpg, cyl %in% c("4", "6")),
  x               = cty,
  y               = manufacturer,
  type            = "bayes", ## Bayesian test
  xlab            = "city miles per gallon",
  ylab            = "car manufacturer",
  grouping.var    = cyl, ## grouping variable
  test.value      = 15.5,
  point.args      = list(color = "red", size = 5, shape = 13),
  annotation.args = list(title = "Fuel economy data")
)
```

<img src="man/figures/README-ggdotplotstats2-1.png" width="100%" />

##### Summary of graphics

| graphical element       | `geom_` used          | argument for further modification |
|-------------------------|-----------------------|-----------------------------------|
| raw data                | `ggplot2::geom_point` | `point.args`                      |
| centrality measure line | `ggplot2::geom_vline` | `centrality.line.args`            |

##### Summary of tests

**Central tendency measure**

| Type           | Measure                                           | Function used                       |
|----------------|---------------------------------------------------|-------------------------------------|
| Parametric     | mean                                              | `parameters::describe_distribution` |
| Non-parametric | median                                            | `parameters::describe_distribution` |
| Robust         | trimmed mean                                      | `parameters::describe_distribution` |
| Bayesian       | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution` |

**Hypothesis testing**

| Type           | Test                                     | Function used          |
|----------------|------------------------------------------|------------------------|
| Parametric     | One-sample Student’s *t*-test            | `stats::t.test`        |
| Non-parametric | One-sample Wilcoxon test                 | `stats::wilcox.test`   |
| Robust         | Bootstrap-*t* method for one-sample test | `WRS2::trimcibt`       |
| Bayesian       | One-sample Student’s *t*-test            | `BayesFactor::ttestBF` |

**Effect size estimation**

| Type           | Effect size                                                                                                           | CI? | Function used                                  |
|----------------|-----------------------------------------------------------------------------------------------------------------------|-----|------------------------------------------------|
| Parametric     | Cohen’s *d*, Hedge’s *g*                                                                                              | ✅  | `effectsize::cohens_d`, `effectsize::hedges_g` |
| Non-parametric | *r* (rank-biserial correlation)                                                                                       | ✅  | `effectsize::rank_biserial`                    |
| Robust         | trimmed mean                                                                                                          | ✅  | `WRS2::trimcibt`                               |
| Bayes Factor   | ![\\delta\_{posterior}](https://chart.apis.google.com/chart?cht=tx&chl=%5Cdelta_%7Bposterior%7D "\delta_{posterior}") | ✅  | `bayestestR::describe_posterior`               |

### `ggscatterstats`

This function creates a scatterplot with marginal distributions overlaid
on the axes and results from statistical tests in the subtitle:

``` r
ggscatterstats(
  data  = ggplot2::msleep,
  x     = sleep_rem,
  y     = awake,
  xlab  = "REM sleep (in hours)",
  ylab  = "Amount of time spent awake (in hours)",
  title = "Understanding mammalian sleep"
)
```

<img src="man/figures/README-ggscatterstats1-1.png" width="100%" />

**Defaults** return<br>

✅ raw data + distributions <br> ✅ marginal distributions <br> ✅
inferential statistics <br> ✅ effect size + CIs <br> ✅ Bayesian
hypothesis-testing <br> ✅ Bayesian estimation <br>

There is also a `grouped_` variant of this function that makes it easy
to repeat the same operation across a **single** grouping variable.

``` r
## for reproducibility
set.seed(123)

## plot
grouped_ggscatterstats(
  data             = dplyr::filter(movies_long, genre %in% c("Action", "Comedy")),
  x                = rating,
  y                = length,
  grouping.var     = genre, ## grouping variable
  label.var        = title,
  label.expression = length > 200,
  xlab             = "IMDB rating",
  ggtheme          = ggplot2::theme_grey(),
  ggplot.component = list(ggplot2::scale_x_continuous(breaks = seq(2, 9, 1), limits = (c(2, 9)))),
  plotgrid.args    = list(nrow = 1),
  annotation.args  = list(title = "Relationship between movie length and IMDB ratings")
)
```

<img src="man/figures/README-ggscatterstats2-1.png" width="100%" />

##### Summary of graphics

| graphical element   | `geom_` used                                                 | argument for further modification            |
|---------------------|--------------------------------------------------------------|----------------------------------------------|
| raw data            | `ggplot2::geom_point`                                        | `point.args`                                 |
| labels for raw data | `ggrepel::geom_label_repel`                                  | `point.label.args`                           |
| smooth line         | `ggplot2::geom_smooth`                                       | `smooth.line.args`                           |
| marginal histograms | `ggside::geom_xsidehistogram`, `ggside::geom_ysidehistogram` | `xsidehistogram.args`, `ysidehistogram.args` |

##### Summary of tests

**Hypothesis testing** and **Effect size estimation**

| Type           | Test                                       | CI? | Function used              |
|----------------|--------------------------------------------|-----|----------------------------|
| Parametric     | Pearson’s correlation coefficient          | ✅  | `correlation::correlation` |
| Non-parametric | Spearman’s rank correlation coefficient    | ✅  | `correlation::correlation` |
| Robust         | Winsorized Pearson correlation coefficient | ✅  | `correlation::correlation` |
| Bayesian       | Pearson’s correlation coefficient          | ✅  | `correlation::correlation` |

For more, see the `ggscatterstats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggscatterstats.html>

### `ggcorrmat`

`ggcorrmat` makes a correlalogram (a matrix of correlation coefficients)
with minimal amount of code. Just sticking to the defaults itself
produces publication-ready correlation matrices. But, for the sake of
exploring the available options, let’s change some of the defaults. For
example, multiple aesthetics-related arguments can be modified to change
the appearance of the correlation matrix.

``` r
## for reproducibility
set.seed(123)

## as a default this function outputs a correlation matrix plot
ggcorrmat(
  data     = ggplot2::msleep,
  colors   = c("#B2182B", "white", "#4D4D4D"),
  title    = "Correlalogram for mammals sleep dataset",
  subtitle = "sleep units: hours; weight units: kilograms"
)
```

<img src="man/figures/README-ggcorrmat1-1.png" width="100%" />

**Defaults** return<br>

✅ effect size + significance<br> ✅ careful handling of `NA`s

If there are `NA`s present in the selected variables, the legend will
display minimum, median, and maximum number of pairs used for
correlation tests.

There is also a `grouped_` variant of this function that makes it easy
to repeat the same operation across a **single** grouping variable:

``` r
## for reproducibility
set.seed(123)

## plot
grouped_ggcorrmat(
  data         = dplyr::filter(movies_long, genre %in% c("Action", "Comedy")),
  type         = "robust", ## correlation method
  colors       = c("#cbac43", "white", "#550000"),
  grouping.var = genre, ## grouping variable
  matrix.type  = "lower" ## type of matrix
)
```

<img src="man/figures/README-ggcorrmat2-1.png" width="100%" />

You can also get a dataframe containing all relevant details from the
statistical tests:

``` r
## setup
set.seed(123)

## tidy data as output
ggcorrmat(
  data   = dplyr::select(ggplot2::msleep, dplyr::matches("sleep|awake")),
  type   = "bayes",
  output = "dataframe"
)
##> # A tibble: 6 x 14
##>   parameter1  parameter2  estimate conf.level conf.low conf.high    pd
##>   <chr>       <chr>          <dbl>      <dbl>    <dbl>     <dbl> <dbl>
##> 1 sleep_total sleep_rem      0.731       0.95    0.606    0.838  1    
##> 2 sleep_total sleep_cycle   -0.432       0.95   -0.681   -0.131  0.995
##> 3 sleep_total awake         -1.00        0.95   -1.00    -1.00   1    
##> 4 sleep_rem   sleep_cycle   -0.304       0.95   -0.576    0.0241 0.963
##> 5 sleep_rem   awake         -0.733       0.95   -0.832   -0.599  1    
##> 6 sleep_cycle awake          0.439       0.95    0.151    0.672  0.998
##>   rope.percentage prior.distribution prior.location prior.scale  bayes.factor
##>             <dbl> <chr>                       <dbl>       <dbl>         <dbl>
##> 1          0      beta                         1.41        1.41 3000790806.  
##> 2          0.0173 beta                         1.41        1.41          8.85
##> 3          0      beta                         1.41        1.41         NA   
##> 4          0.100  beta                         1.41        1.41          1.42
##> 5          0      beta                         1.41        1.41 3005546544.  
##> 6          0.015  beta                         1.41        1.41          8.85
##>   method                       n.obs
##>   <chr>                        <int>
##> 1 Bayesian Pearson correlation    61
##> 2 Bayesian Pearson correlation    32
##> 3 Bayesian Pearson correlation    83
##> 4 Bayesian Pearson correlation    32
##> 5 Bayesian Pearson correlation    61
##> 6 Bayesian Pearson correlation    32
```

Additionally, **partial** correlation are also supported:

``` r
## setup
set.seed(123)

## tidy data as output
ggcorrmat(
  data    = dplyr::select(ggplot2::msleep, dplyr::matches("sleep|awake")),
  type    = "bayes",
  partial = TRUE,
  output  = "dataframe"
)
##> # A tibble: 6 x 14
##>   parameter1  parameter2  estimate conf.level conf.low conf.high    pd
##>   <chr>       <chr>          <dbl>      <dbl>    <dbl>     <dbl> <dbl>
##> 1 sleep_total sleep_rem    0.0938        0.95   -0.249     0.411 0.704
##> 2 sleep_total sleep_cycle  0.00194       0.95   -0.333     0.328 0.505
##> 3 sleep_total awake       -1             0.95   -1        -1     1    
##> 4 sleep_rem   sleep_cycle -0.0274        0.95   -0.357     0.300 0.558
##> 5 sleep_rem   awake        0.0889        0.95   -0.231     0.433 0.698
##> 6 sleep_cycle awake        0.00381       0.95   -0.335     0.324 0.508
##>   rope.percentage prior.distribution prior.location prior.scale bayes.factor
##>             <dbl> <chr>                       <dbl>       <dbl>        <dbl>
##> 1           0.378 beta                         1.41        1.41        0.311
##> 2           0.434 beta                         1.41        1.41        0.267
##> 3           0     beta                         1.41        1.41       NA    
##> 4           0.429 beta                         1.41        1.41        0.272
##> 5           0.386 beta                         1.41        1.41        0.309
##> 6           0.444 beta                         1.41        1.41        0.267
##>   method                       n.obs
##>   <chr>                        <int>
##> 1 Bayesian Pearson correlation    32
##> 2 Bayesian Pearson correlation    32
##> 3 Bayesian Pearson correlation    32
##> 4 Bayesian Pearson correlation    32
##> 5 Bayesian Pearson correlation    32
##> 6 Bayesian Pearson correlation    32
```

##### Summary of graphics

| graphical element  | `geom_` used             | argument for further modification |
|--------------------|--------------------------|-----------------------------------|
| correlation matrix | `ggcorrplot::ggcorrplot` | `ggcorrplot.args`                 |

##### Summary of tests

**Hypothesis testing** and **Effect size estimation**

| Type           | Test                                       | CI? | Function used              |
|----------------|--------------------------------------------|-----|----------------------------|
| Parametric     | Pearson’s correlation coefficient          | ✅  | `correlation::correlation` |
| Non-parametric | Spearman’s rank correlation coefficient    | ✅  | `correlation::correlation` |
| Robust         | Winsorized Pearson correlation coefficient | ✅  | `correlation::correlation` |
| Bayesian       | Pearson’s correlation coefficient          | ✅  | `correlation::correlation` |

For examples and more information, see the `ggcorrmat` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcorrmat.html>

### `ggpiestats`

This function creates a pie chart for categorical or nominal variables
with results from contingency table analysis (Pearson’s chi-squared test
for between-subjects design and McNemar’s chi-squared test for
within-subjects design) included in the subtitle of the plot. If only
one categorical variable is entered, results from one-sample proportion
test (i.e., a chi-squared goodness of fit test) will be displayed as a
subtitle.

To study an interaction between two categorical variables:

``` r
## for reproducibility
set.seed(123)

## plot
ggpiestats(
  data         = mtcars,
  x            = am,
  y            = cyl,
  package      = "wesanderson",
  palette      = "Royal1",
  title        = "Dataset: Motor Trend Car Road Tests", ## title for the plot
  legend.title = "Transmission", ## title for the legend
  caption      = "Source: 1974 Motor Trend US magazine"
)
```

<img src="man/figures/README-ggpiestats1-1.png" width="100%" />

**Defaults** return<br>

✅ descriptives (frequency + %s) <br> ✅ inferential statistics <br> ✅
effect size + CIs <br> ✅ Goodness-of-fit tests <br> ✅ Bayesian
hypothesis-testing <br> ✅ Bayesian estimation <br>

There is also a `grouped_` variant of this function that makes it easy
to repeat the same operation across a **single** grouping variable.
Following example is a case where the theoretical question is about
proportions for different levels of a single nominal variable:

``` r
## for reproducibility
set.seed(123)

## plot
grouped_ggpiestats(
  data         = mtcars,
  x            = cyl,
  grouping.var = am, ## grouping variable
  label.repel  = TRUE, ## repel labels (helpful for overlapping labels)
  package      = "ggsci", ## package from which color palette is to be taken
  palette      = "default_ucscgb" ## choosing a different color palette
)
```

<img src="man/figures/README-ggpiestats2-1.png" width="100%" />

##### Summary of graphics

| graphical element  | `geom_` used                                      | argument for further modification |
|--------------------|---------------------------------------------------|-----------------------------------|
| pie slices         | `ggplot2::geom_col`                               | ❌                                |
| descriptive labels | `ggplot2::geom_label`/`ggrepel::geom_label_repel` | `label.args`                      |

##### Summary of tests

**two-way table**

**Hypothesis testing**

| Type                      | Design   | Test                                                                                                   | Function used                     |
|---------------------------|----------|--------------------------------------------------------------------------------------------------------|-----------------------------------|
| Parametric/Non-parametric | Unpaired | Pearson’s ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test          | `stats::chisq.test`               |
| Bayesian                  | Unpaired | Bayesian Pearson’s ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test | `BayesFactor::contingencyTableBF` |
| Parametric/Non-parametric | Paired   | McNemar’s ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test          | `stats::mcnemar.test`             |
| Bayesian                  | Paired   | ❌                                                                                                     | ❌                                |

**Effect size estimation**

| Type                      | Design   | Effect size                                                         | CI? | Function used           |
|---------------------------|----------|---------------------------------------------------------------------|-----|-------------------------|
| Parametric/Non-parametric | Unpaired | Cramer’s ![V](https://chart.apis.google.com/chart?cht=tx&chl=V "V") | ✅  | `effectsize::cramers_v` |
| Bayesian                  | Unpaired | Cramer’s ![V](https://chart.apis.google.com/chart?cht=tx&chl=V "V") | ✅  | `effectsize::cramers_v` |
| Parametric/Non-parametric | Paired   | Cohen’s ![g](https://chart.apis.google.com/chart?cht=tx&chl=g "g")  | ✅  | `effectsize::cohens_g`  |
| Bayesian                  | Paired   | ❌                                                                  | ❌  | ❌                      |

**one-way table**

**Hypothesis testing**

| Type                      | Test                                                                                                         | Function used       |
|---------------------------|--------------------------------------------------------------------------------------------------------------|---------------------|
| Parametric/Non-parametric | Goodness of fit ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test          | `stats::chisq.test` |
| Bayesian                  | Bayesian Goodness of fit ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test | (custom)            |

**Effect size estimation**

| Type                      | Effect size                                                          | CI? | Function used            |
|---------------------------|----------------------------------------------------------------------|-----|--------------------------|
| Parametric/Non-parametric | Pearson’s ![C](https://chart.apis.google.com/chart?cht=tx&chl=C "C") | ✅  | `effectsize::pearsons_c` |
| Bayesian                  | ❌                                                                   | ❌  | ❌                       |

For more, see the `ggpiestats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggpiestats.html>

### `ggbarstats`

In case you are not a fan of pie charts (for very good reasons), you can
alternatively use `ggbarstats` function which has a similar syntax.

N.B. The *p*-values from one-sample proportion test are displayed on top
of each bar.

``` r
## for reproducibility
set.seed(123)
library(ggplot2)

## plot
ggbarstats(
  data             = movies_long,
  x                = mpaa,
  y                = genre,
  title            = "MPAA Ratings by Genre",
  xlab             = "movie genre",
  legend.title     = "MPAA rating",
  ggtheme          = hrbrthemes::theme_ipsum_pub(),
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  palette          = "Set2"
)
```

<img src="man/figures/README-ggbarstats1-1.png" width="100%" />

**Defaults** return<br>

✅ descriptives (frequency + %s) <br> ✅ inferential statistics <br> ✅
effect size + CIs <br> ✅ Goodness-of-fit tests <br> ✅ Bayesian
hypothesis-testing <br> ✅ Bayesian estimation <br>

And, needless to say, there is also a `grouped_` variant of this
function-

``` r
## setup
set.seed(123)

## plot
grouped_ggbarstats(
  data         = mtcars,
  x            = am,
  y            = cyl,
  grouping.var = vs,
  package      = "wesanderson",
  palette      = "Darjeeling2",
  ggtheme      = ggthemes::theme_tufte(base_size = 12)
)
```

<img src="man/figures/README-ggbarstats2-1.png" width="100%" />

##### Summary of graphics

| graphical element  | `geom_` used          | argument for further modification |
|--------------------|-----------------------|-----------------------------------|
| bars               | `ggplot2::geom_bar`   | ❌                                |
| descriptive labels | `ggplot2::geom_label` | `label.args`                      |

##### Summary of tests

**two-way table**

**Hypothesis testing**

| Type                      | Design   | Test                                                                                                   | Function used                     |
|---------------------------|----------|--------------------------------------------------------------------------------------------------------|-----------------------------------|
| Parametric/Non-parametric | Unpaired | Pearson’s ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test          | `stats::chisq.test`               |
| Bayesian                  | Unpaired | Bayesian Pearson’s ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test | `BayesFactor::contingencyTableBF` |
| Parametric/Non-parametric | Paired   | McNemar’s ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test          | `stats::mcnemar.test`             |
| Bayesian                  | Paired   | ❌                                                                                                     | ❌                                |

**Effect size estimation**

| Type                      | Design   | Effect size                                                         | CI? | Function used           |
|---------------------------|----------|---------------------------------------------------------------------|-----|-------------------------|
| Parametric/Non-parametric | Unpaired | Cramer’s ![V](https://chart.apis.google.com/chart?cht=tx&chl=V "V") | ✅  | `effectsize::cramers_v` |
| Bayesian                  | Unpaired | Cramer’s ![V](https://chart.apis.google.com/chart?cht=tx&chl=V "V") | ✅  | `effectsize::cramers_v` |
| Parametric/Non-parametric | Paired   | Cohen’s ![g](https://chart.apis.google.com/chart?cht=tx&chl=g "g")  | ✅  | `effectsize::cohens_g`  |
| Bayesian                  | Paired   | ❌                                                                  | ❌  | ❌                      |

**one-way table**

**Hypothesis testing**

| Type                      | Test                                                                                                         | Function used       |
|---------------------------|--------------------------------------------------------------------------------------------------------------|---------------------|
| Parametric/Non-parametric | Goodness of fit ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test          | `stats::chisq.test` |
| Bayesian                  | Bayesian Goodness of fit ![\\chi^2](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E2 "\chi^2") test | (custom)            |

**Effect size estimation**

| Type                      | Effect size                                                          | CI? | Function used            |
|---------------------------|----------------------------------------------------------------------|-----|--------------------------|
| Parametric/Non-parametric | Pearson’s ![C](https://chart.apis.google.com/chart?cht=tx&chl=C "C") | ✅  | `effectsize::pearsons_c` |
| Bayesian                  | ❌                                                                   | ❌  | ❌                       |

### `ggcoefstats`

The function `ggcoefstats` generates **dot-and-whisker plots** for
regression models saved in a tidy data frame. The tidy dataframes are
prepared using `parameters::model_parameters`. Additionally, if
available, the model summary indices are also extracted from
`performance::model_performance`.

Although the statistical models displayed in the plot may differ based
on the class of models being investigated, there are few aspects of the
plot that will be invariant across models:

-   The dot-whisker plot contains a dot representing the **estimate**
    and their **confidence intervals** (`95%` is the default). The
    estimate can either be effect sizes (for tests that depend on the
    `F`-statistic) or regression coefficients (for tests with `t`-,
    ![\\chi^{2}](https://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E%7B2%7D "\chi^{2}")-,
    and `z`-statistic), etc. The function will, by default, display a
    helpful `x`-axis label that should clear up what estimates are being
    displayed. The confidence intervals can sometimes be asymmetric if
    bootstrapping was used.

-   The label attached to dot will provide more details from the
    statistical test carried out and it will typically contain estimate,
    statistic, and *p*-value.

-   The caption will contain diagnostic information, if available, about
    models that can be useful for model selection: The smaller the
    Akaike’s Information Criterion (**AIC**) and the Bayesian
    Information Criterion (**BIC**) values, the “better” the model is.

-   The output of this function will be a `{ggplot2}` object and, thus,
    it can be further modified (e.g., change themes, etc.) with
    `{ggplot2}` functions.

``` r
## for reproducibility
set.seed(123)

## model
mod <- stats::lm(formula = mpg ~ am * cyl, data = mtcars)

## plot
ggcoefstats(mod, ggtheme = hrbrthemes::theme_ipsum_ps())
```

<img src="man/figures/README-ggcoefstats1-1.png" width="100%" />

**Defaults** return<br>

✅ inferential statistics <br> ✅ estimate + CIs <br> ✅ model summary
(AIC and BIC) <br>

##### Supported models

Most of the regression models that are supported in the underlying
packages are also supported by `ggcoefstats`.

``` r
insight::supported_models()
##>   [1] "aareg"             "afex_aov"          "AKP"              
##>   [4] "Anova.mlm"         "aov"               "aovlist"          
##>   [7] "Arima"             "averaging"         "bamlss"           
##>  [10] "bamlss.frame"      "bayesQR"           "bayesx"           
##>  [13] "BBmm"              "BBreg"             "bcplm"            
##>  [16] "betamfx"           "betaor"            "betareg"          
##>  [19] "BFBayesFactor"     "bfsl"              "BGGM"             
##>  [22] "bife"              "bifeAPEs"          "bigglm"           
##>  [25] "biglm"             "blavaan"           "blrm"             
##>  [28] "bracl"             "brglm"             "brmsfit"          
##>  [31] "brmultinom"        "btergm"            "censReg"          
##>  [34] "cgam"              "cgamm"             "cglm"             
##>  [37] "clm"               "clm2"              "clmm"             
##>  [40] "clmm2"             "clogit"            "coeftest"         
##>  [43] "complmrob"         "confusionMatrix"   "coxme"            
##>  [46] "coxph"             "coxph.penal"       "coxr"             
##>  [49] "cpglm"             "cpglmm"            "crch"             
##>  [52] "crq"               "crqs"              "crr"              
##>  [55] "dep.effect"        "DirichletRegModel" "drc"              
##>  [58] "eglm"              "elm"               "epi.2by2"         
##>  [61] "ergm"              "feglm"             "feis"             
##>  [64] "felm"              "fitdistr"          "fixest"           
##>  [67] "flexsurvreg"       "gam"               "Gam"              
##>  [70] "gamlss"            "gamm"              "gamm4"            
##>  [73] "garch"             "gbm"               "gee"              
##>  [76] "geeglm"            "glht"              "glimML"           
##>  [79] "glm"               "Glm"               "glmm"             
##>  [82] "glmmadmb"          "glmmPQL"           "glmmTMB"          
##>  [85] "glmrob"            "glmRob"            "glmx"             
##>  [88] "gls"               "gmnl"              "HLfit"            
##>  [91] "htest"             "hurdle"            "iv_robust"        
##>  [94] "ivFixed"           "ivprobit"          "ivreg"            
##>  [97] "lavaan"            "lm"                "lm_robust"        
##> [100] "lme"               "lmerMod"           "lmerModLmerTest"  
##> [103] "lmodel2"           "lmrob"             "lmRob"            
##> [106] "logistf"           "logitmfx"          "logitor"          
##> [109] "LORgee"            "lqm"               "lqmm"             
##> [112] "lrm"               "manova"            "MANOVA"           
##> [115] "margins"           "maxLik"            "mclogit"          
##> [118] "mcmc"              "mcmc.list"         "MCMCglmm"         
##> [121] "mcp1"              "mcp12"             "mcp2"             
##> [124] "med1way"           "mediate"           "merMod"           
##> [127] "merModList"        "meta_bma"          "meta_fixed"       
##> [130] "meta_random"       "metaplus"          "mhurdle"          
##> [133] "mipo"              "mira"              "mixed"            
##> [136] "MixMod"            "mixor"             "mjoint"           
##> [139] "mle"               "mle2"              "mlm"              
##> [142] "mlogit"            "mmlogit"           "model_fit"        
##> [145] "multinom"          "mvord"             "negbinirr"        
##> [148] "negbinmfx"         "ols"               "onesampb"         
##> [151] "orm"               "pgmm"              "plm"              
##> [154] "PMCMR"             "poissonirr"        "poissonmfx"       
##> [157] "polr"              "probitmfx"         "psm"              
##> [160] "Rchoice"           "ridgelm"           "riskRegression"   
##> [163] "rjags"             "rlm"               "rlmerMod"         
##> [166] "RM"                "rma"               "rma.uni"          
##> [169] "robmixglm"         "robtab"            "rq"               
##> [172] "rqs"               "rqss"              "Sarlm"            
##> [175] "scam"              "selection"         "sem"              
##> [178] "SemiParBIV"        "semLm"             "semLme"           
##> [181] "slm"               "speedglm"          "speedlm"          
##> [184] "stanfit"           "stanmvreg"         "stanreg"          
##> [187] "summary.lm"        "survfit"           "survreg"          
##> [190] "svy_vglm"          "svyglm"            "svyolr"           
##> [193] "t1way"             "tobit"             "trimcibt"         
##> [196] "truncreg"          "vgam"              "vglm"             
##> [199] "wbgee"             "wblm"              "wbm"              
##> [202] "wmcpAKP"           "yuen"              "yuend"            
##> [205] "zcpglm"            "zeroinfl"          "zerotrunc"
```

Although not shown here, this function can also be used to carry out
parametric, robust, and Bayesian random-effects meta-analysis.

##### Summary of graphics

| graphical element              | `geom_` used                | argument for further modification |
|--------------------------------|-----------------------------|-----------------------------------|
| regression estimate            | `ggplot2::geom_point`       | `point.args`                      |
| error bars                     | `ggplot2::geom_errorbarh`   | `errorbar.args`                   |
| vertical line                  | `ggplot2::geom_vline`       | `vline.args`                      |
| label with statistical details | `ggrepel::geom_label_repel` | `stats.label.args`                |

##### Summary of meta-analysis tests

**Hypothesis testing** and **Effect size estimation**

| Type       | Test                                             | Effect size                                                               | CI? | Function used          |
|------------|--------------------------------------------------|---------------------------------------------------------------------------|-----|------------------------|
| Parametric | Meta-analysis via random-effects models          | ![\\beta](https://chart.apis.google.com/chart?cht=tx&chl=%5Cbeta "\beta") | ✅  | `metafor::metafor`     |
| Robust     | Meta-analysis via robust random-effects models   | ![\\beta](https://chart.apis.google.com/chart?cht=tx&chl=%5Cbeta "\beta") | ✅  | `metaplus::metaplus`   |
| Bayes      | Meta-analysis via Bayesian random-effects models | ![\\beta](https://chart.apis.google.com/chart?cht=tx&chl=%5Cbeta "\beta") | ✅  | `metaBMA::meta_random` |

For a more exhaustive account of this function, see the associated
vignette-
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcoefstats.html>

### Extracting dataframes with statistical details

`{ggstatsplot}` also offers a convenience function to extract dataframes
with statistical details that are used to create expressions displayed
in `{ggstatsplot}` plots.

``` r
set.seed(123)

## a list of tibbles containing statistical analysis summaries
ggbetweenstats(mtcars, cyl, mpg) %>%
  extract_stats()
##> $subtitle_data
##> # A tibble: 1 x 14
##>   statistic    df df.error    p.value
##>       <dbl> <dbl>    <dbl>      <dbl>
##> 1      31.6     2     18.0 0.00000127
##>   method                                                   effectsize estimate
##>   <chr>                                                    <chr>         <dbl>
##> 1 One-way analysis of means (not assuming equal variances) Omega2        0.744
##>   conf.level conf.low conf.high conf.method conf.distribution n.obs expression  
##>        <dbl>    <dbl>     <dbl> <chr>       <chr>             <int> <list>      
##> 1       0.95    0.531         1 ncp         F                    32 <expression>
##> 
##> $caption_data
##> # A tibble: 6 x 17
##>   term     pd rope.percentage prior.distribution prior.location prior.scale
##>   <chr> <dbl>           <dbl> <chr>                       <dbl>       <dbl>
##> 1 mu    1              0      cauchy                          0       0.707
##> 2 cyl-4 1              0      cauchy                          0       0.707
##> 3 cyl-6 0.780          0.390  cauchy                          0       0.707
##> 4 cyl-8 1              0      cauchy                          0       0.707
##> 5 sig2  1              0      cauchy                          0       0.707
##> 6 g_cyl 1              0.0418 cauchy                          0       0.707
##>       bf10 method                          log_e_bf10 effectsize        
##>      <dbl> <chr>                                <dbl> <chr>             
##> 1 3008850. Bayes factors for linear models       14.9 Bayesian R-squared
##> 2 3008850. Bayes factors for linear models       14.9 Bayesian R-squared
##> 3 3008850. Bayes factors for linear models       14.9 Bayesian R-squared
##> 4 3008850. Bayes factors for linear models       14.9 Bayesian R-squared
##> 5 3008850. Bayes factors for linear models       14.9 Bayesian R-squared
##> 6 3008850. Bayes factors for linear models       14.9 Bayesian R-squared
##>   estimate std.dev conf.level conf.low conf.high n.obs expression  
##>      <dbl>   <dbl>      <dbl>    <dbl>     <dbl> <int> <list>      
##> 1    0.714  0.0503       0.95    0.574     0.788    32 <expression>
##> 2    0.714  0.0503       0.95    0.574     0.788    32 <expression>
##> 3    0.714  0.0503       0.95    0.574     0.788    32 <expression>
##> 4    0.714  0.0503       0.95    0.574     0.788    32 <expression>
##> 5    0.714  0.0503       0.95    0.574     0.788    32 <expression>
##> 6    0.714  0.0503       0.95    0.574     0.788    32 <expression>
##> 
##> $pairwise_comparisons_data
##> # A tibble: 3 x 11
##>   group1 group2 statistic   p.value alternative method            distribution
##>   <chr>  <chr>      <dbl>     <dbl> <chr>       <chr>             <chr>       
##> 1 4      6          -6.67 0.00110   two.sided   Games-Howell test q           
##> 2 4      8         -10.7  0.0000140 two.sided   Games-Howell test q           
##> 3 6      8          -7.48 0.000257  two.sided   Games-Howell test q           
##>   p.adjustment test.details      p.value.adjustment
##>   <chr>        <chr>             <chr>             
##> 1 none         Games-Howell test Holm              
##> 2 none         Games-Howell test Holm              
##> 3 none         Games-Howell test Holm              
##>   label                                     
##>   <chr>                                     
##> 1 list(~italic(p)[Holm-corrected]==1.10e-03)
##> 2 list(~italic(p)[Holm-corrected]==1.40e-05)
##> 3 list(~italic(p)[Holm-corrected]==2.57e-04)
##> 
##> $descriptive_data
##> NULL
##> 
##> $one_sample_data
##> NULL
```

Note that all of this analysis is carried out by `{statsExpressions}`
package: <https://indrajeetpatil.github.io/statsExpressions/>

### Using `{ggstatsplot}` statistical details with custom plots

Sometimes you may not like the default plots produced by
`{ggstatsplot}`. In such cases, you can use other **custom** plots (from
`{ggplot2}` or other plotting packages) and still use `{ggstatsplot}`
functions to display results from relevant statistical test.

For example, in the following chunk, we will create plot (*ridgeplot*)
using `ggridges` package and use `{ggstatsplot}` function for extracting
results.

``` r
## loading the needed libraries
set.seed(123)
library(ggridges)
library(ggplot2)
library(ggstatsplot)

## using `{ggstatsplot}` to get call with statistical results
stats_results <-
  ggbetweenstats(
    data = morley,
    x = Expt,
    y = Speed,
    output = "subtitle"
  )

## using `ggridges` to create plot
ggplot(morley, aes(x = Speed, y = as.factor(Expt), fill = as.factor(Expt))) +
  geom_density_ridges(
    jittered_points = TRUE,
    quantile_lines = TRUE,
    scale = 0.9,
    alpha = 0.7,
    vline_size = 1,
    vline_color = "red",
    point_size = 0.4,
    point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)
  ) + ## adding annotations
  labs(
    title = "Michelson-Morley experiments",
    subtitle = stats_results,
    x = "Speed of light",
    y = "Experiment number"
  ) + ## remove the legend
  theme(legend.position = "none")
```

<img src="man/figures/README-ridgeplot-1.png" width="100%" />

## Summary of benefits of using `{ggstatsplot}`

-   No need to use scores of packages for statistical analysis (e.g.,
    one to get stats, one to get effect sizes, another to get Bayes
    Factors, and yet another to get pairwise comparisons, etc.).

-   Minimal amount of code needed for all functions (typically only
    `data`, `x`, and `y`), which minimizes chances of error and makes
    for tidy scripts.

-   Conveniently toggle between statistical approaches.

-   Truly makes your figures worth a thousand words.

-   No need to copy-paste results to the text editor (MS-Word, e.g.).

-   Disembodied figures stand on their own and are easy to evaluate for
    the reader.

-   More breathing room for theoretical discussion and other text.

-   No need to worry about updating figures and statistical details
    separately.

## Misconceptions about `{ggstatsplot}`

This package is…

❌ an alternative to learning `{ggplot2}`<br> ✅ (The better you know
`{ggplot2}`, the more you can modify the defaults to your liking.)

❌ meant to be used in talks/presentations<br> ✅ (Default plots can be
too complicated for effectively communicating results in
time-constrained presentation settings, e.g. conference talks.)

❌ the only game in town<br> ✅ (GUI software alternatives:
[JASP](https://jasp-stats.org/) and [jamovi](https://www.jamovi.org/)).

## Extensions

In case you use the GUI software [`jamovi`](https://www.jamovi.org/),
you can install a module called
[`jjstatsplot`](https://github.com/sbalci/jjstatsplot), which is a
wrapper around `{ggstatsplot}`.

## Acknowledgments

I would like to thank all the contributors to `{ggstatsplot}` who
pointed out bugs or requested features I hadn’t considered. I would
especially like to thank other package developers (especially Daniel
Lüdecke, Dominique Makowski, Mattan S. Ben-Shachar, Patrick Mair,
Salvatore Mangiafico, etc.) who have patiently and diligently answered
my relentless number of questions and added feature requests I wanted. I
also want to thank Chuck Powell for his initial contributions to the
package.

The hexsticker was generously designed by Sarah Otterstetter (Max Planck
Institute for Human Development, Berlin). This package has also
benefited from the larger `rstats` community on Twitter and
`StackOverflow`.

Thanks are also due to my postdoc advisers (Mina Cikara and Fiery
Cushman at Harvard University; Iyad Rahwan at Max Planck Institute for
Human Development) who patiently supported me spending hundreds (?) of
hours working on this package rather than what I was paid to do. 😁

## Contributing

I’m happy to receive bug reports, suggestions, questions, and (most of
all) contributions to fix problems and add features. I personally prefer
using the `GitHub` issues system over trying to reach out to me in other
ways (personal e-mail, Twitter, etc.). Pull Requests for contributions
are encouraged.

Here are some simple ways in which you can contribute (in the increasing
order of commitment):

-   Read and correct any inconsistencies in the
    [documentation](https://indrajeetpatil.github.io/ggstatsplot/)
-   Raise issues about bugs or wanted features
-   Review code
-   Add new functionality (in the form of new plotting functions or
    helpers for preparing subtitles)

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/IndrajeetPatil/ggstatsplot/blob/master/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.
# ggstatsplot 0.9.1.9000

# ggstatsplot 0.9.1

N.B. All statistical analysis in `{ggstatsplot}` is carried out in
`{statsExpressions}`. Thus, to see changes related to statistical expressions,
read the `NEWS` for that package:
<https://indrajeetpatil.github.io/statsExpressions/news/index.html>

MAJOR CHANGES

  - Moves `{PMCMRplus}` package from Imports to Suggests. So, if, as a user, you
    wish to use pairwise comparisons in `ggbetweenstats()` and
    `ggwithinstats()`, you will need to download this package.

MINOR CHANGES

  - To keep the documentation maintainable, a number of vignettes have either
    been removed or they are no longer evaluated and only code is reported.

# ggstatsplot 0.9.0

NEW FEATURES

  - The `pairwise_comparisons()` function for carrying out one-way pairwise
    comparisons has now moved in `{ggstatsplot}` from `{pairwiseComparisons}`
    package.

BREAKING CHANGES

  - A number of effect size estimates and their confidence intervals have
    changed due to respective changes made in `{effectsize}` package version
    `0.5` release. For full details of these changes, see:
    <https://easystats.github.io/effectsize/news/index.html>

  - For the same reason, the effect size for one-way contingency table has
    changed from Cramer's *V* to Pearson's *C*.

MAJOR CHANGES

  - For plotting marginal distributions in `ggscatterstats`, `{ggstatsplot}` now
    relies on `ggside` package instead of `ggExtra`. This was done to remove a
    glaring inconsistency in the API. All functions in `{ggstatsplot}` produced
    `ggplot` objects and could be further modified with `ggplot2` functions,
    except `ggscatterstats`, which led to a lot of confusion among users (e.g.
    #28). This change gets rid of this inconsistency. But it comes at a cost:
    there is no more `marginal.type` argument that lets you change the type of
    marginal distribution graphic and histogram is the only possible option.
    Note that this is **not** a breaking change. Your past code will continue to
    work but it will now always produce a histogram instead of other marginal
    graphic you might have chosen.

  - Minimum needed R version is now `4.0`.

MINOR CHANGES

  - Online vignette about `combine_plots` has been removed. In case you want to
    create a grid of plots, it is highly recommended that you use `patchwork`
    package directly and not this wrapper around it which is mostly useful with
    `{ggstatsplot}` plots.

  - `ggscatterstats` labeling arguments accept only unquoted inputs now, and not
    quoted or string inputs. Allowing this was a bad design choice in the past
    since most functions in `{ggstatsplot}`, inspired by `tidyverse`, expect
    unquoted (`x`) - and not quoted (`"x"`) - arguments. So this function was
    the odd one out.

  - Gets rid of `ipmisc` dependency.

  - Removes `movies_wide` dataset, which was virtually identical to
    `movies_long` dataset and was not used anywhere in the package. Also removes
    the unused `VR_dilemma` dataset.

# ggstatsplot 0.8.0

NEW FUNCTIONS

  - Adds `extract_stats` function to extract dataframes containing statistical
    details.

MAJOR CHANGES

  - There is finally a publication for `{ggstatsplot}` package!
    <https://joss.theoj.org/papers/10.21105/joss.03167>

  - The `ggcoefstats` function defaults to `NULL` for `xlab` and `ylab`
    arguments, which lets users change these labels if they wish to do so.
    Additionally, the x-axis label, if not specified, now defaults to
    `"estimate"`. Whether this estimate corresponds to regression coefficient or
    effect size like partial eta-squared should be clear from the label itself.

  - To reduce the dependency load, `ggcorrplot` moves from `Imports` to
    `Suggests`.

  - The `bar.fill` argument in `gghistostats` is retired in favor of the new
    `bin.args` argument that can be used to pass aesthetic arguments to
    `ggplot2::stat_bin`.

  - `ggstatsplot.layer` argument has been retired. If the user _chooses_ a
    certain `ggplot2` theme, it means they _want_ that theme, and not
    `{ggstatsplot}`'s varnish on it. So the previous behavior was undesirable.
    This is a backward compatible change, so the plots should not look
    different.

MINOR CHANGES

  - The `pch` size for `ggcorrmat` has been increased to 14 (#579) to increase
    its visibility compared to the correlation value text.

  - `ggwithinstats` gains `point.args` to change `geom_point`.

  - Minor change to `ggcorrmat` legend title - content in parentheses is now
    shown outside of it.

BUG FIXES

  - `ggcoefstats` didn't work when statistic for the given model was
    chi-squared. This has been fixed.

# ggstatsplot 0.7.2

MAJOR CHANGES

  - To reduce the dependency load, `ggExtra` moves from `Imports` to
    `Suggests`.

  - All functions are more *robust* in the sense that when statistical analysis
    fails, they will return only the plots with no subtitles/captions. This
    helps avoid difficult-to-diagnose edge case failures when the primary
    functions are used in `grouped_` functions (e.g., #559). The `ggpiestats`
    and `ggbarstats` functions always behaved this way, but the rest of the
    functions now also mimic this behavior.

MINOR CHANGES

  - The `ggcoefstats` labels do not contain degrees of freedom when they are not
    available instead of displaying `Inf`.

# ggstatsplot 0.7.1

MAJOR CHANGES

  - Based on feedback from the users, the argument `title.prefix` is now
    removed. It led to redundant title prefixes across different facets of the
    plot. Given that `grouped_` functions require users to set `grouping.var`,
    it is fair to assume what variable the levels in the title correspond to.

MINOR CHANGES

  - Adapts to changes made in `statsExpressions 1.0.0`.

  - `sample.size.label` argument is retired for `ggbetweenstats`,
    `ggwithinstats`, and `ggbarstats`. I do not think it is ever a good idea to
    not do this. If the users wish to not display sample sizes, they can easily
    do this using `scale_*` functions from `ggplot2`.

  - In `ggpiestats` and `ggbarstats`, parametric proportion tests are now turned
    off when `type = "bayes"`.

# ggstatsplot 0.7.0

BREAKING CHANGES

  - `combine_plots` has been completely revised to rely not on `patchwork`, but
    on `patchwork`, to combine a list of `ggplot` together. This was done to
    have a leaner syntax. With this revision, its vestigial twin `combine_plots`
    is no longer needed and has been removed. This should not break any of the
    existing instances of `grouped_` functions, although it will lead to changed
    graphical layouts. The only instance in which this change will lead to a
    breakage is when you specified `labels` argument. So, if you used
    `plotgrid.args = list(labels = "auto")`, you will now have to replace it
    with `plotgrid.args = list(tag_level = "keep")`. You can also use
    `annotation.args` (e.g., `annotation.args = list(tag_levels = "a")` to
    customize labels (this will create labels with pattern `a`, `b`, `c`, etc.).
    Another instance of breakage is if you had used `combine_plots` function and
    provided individual plots to `...` instead as a `list`.

  - To avoid confusion among users, the default trimming level for all functions
    is now changed from `tr = 0.1` to `tr = 0.2` (which is what `WRS2` defaults
    to).

MAJOR CHANGES

  - All robust tests in this package were based on trimmed means, except for
    correlation test. This has been changed: the robust correlation measure is
    now Winsorized correlation, which is based on trimming. Therefore, the
    `beta` argument has been replaced by `tr` argument. This should result only
    in minor changes in correlation coefficient estimates.

  - Using `annotate` instead of `geom_label` had significantly slowed down
    `gghistostats` and `ggdotplotstats` functions. This has been fixed.

  - Removes the vestigial `notch` and `notchwidth` arguments for
    `ggbetweenstats` and `ggwithinstats`.

  - All Bayesian expression templates are now explicit about the type of
    estimate being displayed.

  - For `gghistostats` and `ggdotplotstats`, the centrality measure labels used
    to be attached to the vertical line, but this occluded the underlying data.
    Now this label is instead shown on the top `x`-axis. Note that this means
    that if you make any further changes to the resulting plot using the
    `ggplot2::scale_x_continuous` function, this label will likely disappear.
    The `centrality.k` argument is retired in favor of `k`.

NEW FEATURES

  - More models supported in `ggcoefstats`: `crr`, `eglm`, `elm`, `varest`.

  - `ggbetweenstats`, `ggwithinstats`, `gghistostats`, `ggdotplotstats` gain
    argument `centrality.type` that can be used to specify which centrality
    parameter is to be displayed. So one can have `type = "robust"` and still
    show median as centrality parameter by choosing `centrality.type =
    "nonparametric"`.

# ggstatsplot 0.6.8

MAJOR CHANGES

  - `gghistostats` removes `bar.measure` argument. The function now defaults to
    showing the `count` information on the `x`-axis and the `proportion`
    information on the duplicated `x`-axis.

  - `ggscatterstats` removes `method` and `method.args` arguments. It will no
    longer be possible to use this function to visualize data for when the model
    is not linear. It also retires `margins` argument.

  - For `ggbetweenstats` and `ggwithinstats` functions, the arguments of type
    `mean.*` have all been replaced by `centrality.*`. This is because now these
    functions decide which central tendency measure to show depending on the
    `type` argument (**mean** for parametric, **median** for non-parametric,
    **trimmed mean** for robust, and **MAP estimator** for Bayes).

  - Similarly, `gghistostats` and `ggdotplotstats` functions also decide which
    central tendency measure to show depending on the `type` argument (**mean**
    for parametric, **median** for non-parametric, **trimmed mean** for robust,
    and **MAP estimator** for Bayes). Therefore, `centrality.parameter` argument
    has been removed. If you want to turn off displaying centrality measure, set
    `centrality.plotting = FALSE`.

  - `gghistostats` and `ggdotplotstats` functions remove the functionality to
    display a vertical line corresponding to `test.value`. This feature was
    turned off by default in prior releases. Accordingly, all related arguments
    from these two functions have been removed.

  - `ggscatterstats` defaults to `densigram` as the marginal distribution
    visualization.

  - `ggbetweenstats` and `ggwithinstats` now display the centrality tendency
    measure in such a way that the label doesn't occlude any of the raw data
    points (#429).

  - `mean.ci` argument is retired for `ggbetweenstats` and `ggwithinstats`.
    Future `{ggstatsplot}` releases will be providing different centrality
    measures depending on the `type` argument and it is not guaranteed that all
    of them will have CIs available. So, for the sake of consistency, this
    argument is just going to be retired.

MINOR CHANGES

  - `ggcorrmat` uses pretty formatting to display sample size information.

  - `ggcoefstats` now also displays degrees of freedom for chi-squared tests.

  - Expects minor changes in some of the effect sizes and their confidence
    intervals due to changes in `{statsExpressions}`.

NEW FEATURES

  - More models supported in `ggcoefstats`: `fixest`, `ivFixed`, `ivprobit`,
    `riskRegression`.

  - `ggcorrmat` supports partial correlations.

# ggstatsplot 0.6.6

BREAKING CHANGES

  - `ggcoefstats` no longer supports `exponentiate` argument. If it is
    specified, the user will have to themselves adjust the scales
    appropriately.

  - `ggcorrmat` defaults have changed significantly:

    1. As a matter of good practice, the *p*-values are adjusted by default for
       multiple comparisons.

    2. The default matrix is upper type, and not the full matrix, which features
       many redundant comparisons and self-correlations diagonally.

    3. Default text size for legend has been increased to 15 and background grid
       has been removed.

BUG FIXES

  - In the prior release, when the GitHub version of `BayesFactor` wasn't
    present, `ggwithinstats` just outright failed to run for ANOVA designs. This
    has been fixed.

  - Setting `mean.path = FALSE` in `ggwithinstats` produced incorrect colors for
    points (#470). This bug was introduced in `0.6.5` and is now fixed.

  - If user had set `options(scipen = 999)` in their session, the *p*-value
    formatting for `ggpiestats` and `ggcoefstats` looked super-ugly (#478). This
    has been fixed.

MAJOR CHANGES

  - Drops `broomExtra` from dependencies. All regression modeling-related
    analysis now relies on `easystats` ecosystem.

  - `ggpiestats` and `ggbarstats` don't support returning dataframes. See FAQ
    vignette on how to get these dataframes:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/faq.html#faq-1>

  - `ggpiestats` and `ggbarstats` were not supposed to support returning Bayes
    Factor for paired contingency table analysis, which is not supported in
    `BayesFactor` itself.

  - `ggcoefstats` defaults to displaying the intercept term. Also, when the
    degrees of freedom are not available for `t`-statistic, they are displayed
    to be `Inf`, in keeping with `easystats` conventions.

  - Instead of showing significance of *p*-values with APA's asterisks
    conventions, `ggbarstats` now instead shows the actual *p*-values from
    one-sample proportion tests.

NEW FEATURES

  - More models supported in `ggcoefstats`: `Glm`.

# ggstatsplot 0.6.5

BREAKING CHANGES

  - `ggpiestats` and `ggbarstats` no longer have the vestigial arguments `main`
    and `condition`, which are superseded by `x` and `y`, respectively.

MAJOR CHANGES

  - For consistency and to reduce confusion, all Bayes Factor (irrespective of
    whether in the subtitle or caption) are always in favor of null over
    alternative (`BF01`).

  - Retires centrality parameter tagging functionality of `ggscatterstats`.
    Although it was not the default, when turned on, it definitely created a
    cluttered plot.

# ggstatsplot 0.6.1

MAJOR CHANGES

  - `ggbetweenstats` and `ggwithinstats` functions now default to
    `pairwise.comparisons = TRUE`.

MINOR CHANGES

  - Plot borders are now removed from the default theme.

  - Small *p*-values (< 0.001) are now displayed in scientific notation.

BREAKING CHANGES

  - `pairwiseComparisons` re-exports are deprecated.

# ggstatsplot 0.6.0

NEW FEATURES

  - More models supported in `ggcoefstats`: `BFBayesFactor`, `betamfx`, `crq`,
    `coxph.penal`, `geeglm`, `glht`, `glmm`, `lm_robust`, `lqm`, `lqmm`,
    `manova`, `maov`, `margins`, `negbinmfx`, `logitmfx`, `logitsf`, `margins`,
    `poissonmfx`, `betaor`, `negbinirr`, `logitor`, `metafor`, `metaplus`,
    `orm`, `poissonirr`, `semLm`, `semLme`, `vgam`.

  - `ggpiestats` gains `label.repel` argument to cover contexts in which the
    labels might overlap. Setting it to `TRUE` will minimize such an overlap.

  - `ggbetweenstats` and `ggwithinstats` gain `ggsignif.args` argument to make
    it easy to change aesthetics of the pairwise comparison geom.

  - The subtitle and caption for Bayes Factor tests now also provide information
    about posterior estimates, when relevant.

MAJOR CHANGES

  - Removed unused `intent_morality` dataset.

  - `ggcoefstats` retires `caption.summary` argument. So, by default, the
    caption is going to contain as much information as it can and the users can
    then choose to modify the default caption using `ggplot2` functions.

MINOR CHANGES

  - The argument `method` for `ggcorrmat` has been renamed to `matrix.method`,
    since it was confusing whether this method referred to correlation method.

  - For both `ggpiestats` and `ggbarstats`, the count labels no longer include `
    n = ` in them as this was confusing since all labels had ` n = ` in them
    with no further explanation about how this `n` differed from `n` in the
    proportion test.

  - No longer relies on `groupedstats` package.

# ggstatsplot 0.5.0

BREAKING CHANGES

  - The `pairwise.annotation` argument for `ggbetweenstats` and `ggwithinstats`
    is deprecated. This was done because-

    1. Different fields have different schema for what significance levels
       asterisks represent.

    2. The *p*-value labels also contain information about whether they are
       adjusted for multiple comparisons.

  - The `normality_message` and `bartlett_message` helper functions have been
    removed. This is because model assumption checks don't really fall under the
    purview of this package. There are excellent visualization tools out there
    for model assumption checks (`ggResidpanel`, `performance`, `DHARMa`,
    `olsrr`, etc.), which should be preferred over unhelpful messages with only
    *p*-values that these functions were printing. For what it's worth, the
    functions where these messages were displayed (`ggbetweenstats` or
    `ggwithinstats`) feature visualizations rich enough and defaults sensible
    enough that most of the time one can either assess these assumptions from
    the plots or need not worry about them.

MAJOR CHANGES

  - `ggcoefstats` has been refactored to reflect that
    `broomExtra::tidy_parameters` now defaults to `parameters` package instead
    of `broom`. It also loses the following vestigial arguments:
    `p.adjust.method` and `coefficient.type`.

  - Reverts aligning title and subtitle with the plot and not the axes, since it
    looked pretty ugly (esp., `ggcoefstats`) and was causing problems for
    labels.

  - `factor.levels` (for `ggpiestats`) and `labels.legend` (for `ggbarstats`)
    are deprecated. If users would like to changes the names for factor levels,
    this should be done outside of `{ggstatsplot}`.

  - The non-parametric post hoc test for between-subjects design has been
    changed from Dwass-Steel-Crichtlow-Fligner test to Dunn test.

NEW FEATURES

  - More models supported in `ggcoefstats`: `bayesGARCH`, `clm2`, `clmm2`,
    `mcmc.list`, `robmixglm`.

# ggstatsplot 0.4.0

BREAKING CHANGES

  - `ggcorrmat` no longer returns matrices of correlation coefficients or other
    details. It now returns either a plot or a dataframe and this can dataframe
    can then be used to create matrices.

  - `ggbarstats` loses `x.axis.orientation` argument. This argument was supposed
    to help avoid overlapping *x*-axis label, but now `ggplot2 3.3.0` has a
    better way to handle this:
    <https://www.tidyverse.org/blog/2020/03/ggplot2-3-3-0/#rewrite-of-axis-code>

NEW FEATURES

  - More models supported in `ggcoefstats`: `bayesx`, `BBmm`, `brmultinom`,
    `lmerModLmerTest`, `lrm`.

  - Specifying `output = "proptest"` for `ggpiestats` and `ggbarstats` functions
    will now return a dataframe containing results from proportion test.

  - `ggbetweenstats` and `ggwithinstats` will display pairwise comparisons even
    if `results.subtitle` is set to `FALSE`.

  - `ggcorrmat` supports computing Bayes Factors for Pearson's *r* correlation.

  - `ggbetweenstats` and `ggwithinstats` now support pairwise comparisons for
    Bayes Factor test.

MAJOR CHANGES

  - For changes related to subtitle details, see changes made in new version of
    `statsExpressions 4.0.0`:
    <https://CRAN.R-project.org/package=statsExpressions/news/news.html>

  - `ggbetweenstats` and `ggwithinstats` no longer print dataframes containing
    results from pairwise comparisons tests because this is too cluttering for
    the user's console. The users are now instead advised to either extract this
    dataframe using `ggplot2::ggplot_build()` function or use the
    `pairwiseComparisons::pairwise_comparisons()` function used in the
    background by `{ggstatsplot}` to carry out this analysis.

  - Due to changes in one of the downstream dependencies, `{ggstatsplot}` now
    expects the minimum R version to be `3.6.0`.

MINOR CHANGES

  - `ggcorrmat` now internally relies on `correlation` for correlation
    analyses.

  - `ggbarstats` no longer displays `"percent"` for Y-axis label as this was
    redundant information.

  - Continuing the argument cleanup that began in `0.3.0`, `ggcoefstats` gains
    `point.args` argument instead of individuals `point.*` arguments.

  - The subtitles are more explicit about the details of the test. For the same
    reason `stat.title` argument from all relevant functions is retired since
    this argument was supposed to be for entering some additional details about
    the test. Additionally, the plot titles and subtitles for some of the plots
    are aligned with the plot.

  - `ggcorrmat` legend, in case of missing values, shows mode - instead of
    median - for the distribution of sample pairs.

  - The following vestigial arguments are retired:

      - `caption.default` in `ggcorrmat`

      - `k.caption.summary` in `ggcoefstats`

# ggstatsplot 0.3.1

This is a hotfix release to correct some of the failing tests and other minor
breakages resulting from the new release of `ggplot2 3.3.0`.

MAJOR CHANGES

  - `ggpiestats` loses `sample.size.label` argument since this information is
    included in the goodness of fit test results itself. So setting
    `proportion.test` to `FALSE` will suppress this information.

# ggstatsplot 0.3.0

BREAKING CHANGES

To give users more flexibility in terms of modifying the aesthetic defaults for
**all** `geoms` included in the `{ggstatsplot}` plots (each plot typically has
multiple geoms), the package now uses a new form of syntax. Previously, each
`geom` had a separate argument to specify each aesthetic (e.g., `geom_point`
would get arguments like `point.size`, `point.color`, etc.), which resulted in
functions with a massive number of arguments and was unsustainable in the long
run. Instead, `{ggstatsplot}` functions now expect a list of such arguments for
the respective geom (e.g., `geom_point` will have `point.args` argument where a
list of arguments `list(size = 5, color = "darkgreen", alpha = 0.8)` can be
supplied).

  - All `grouped_` functions have been refactored to reduce the number of
    arguments. These functions now internally use the new `combine_plots`
    instead of `combine_plots`. The additional arguments to primary functions
    can be provided through `...`. These changes will not necessarily break the
    existing code but will lead to some minor graphical changes (e.g., if you
    were providing `labels` argument explicitly, it will be ignored).

  - All functions lose the `return` argument, which was supposed to be
    alternative to enter `output`. But this was just leading to more confusion
    on the user's part. The biggest user-visible impact this is going to have is
    that `ggcorrmat` will no longer be backward-compatible. The older scripts
    will still work but if the `return` argument was anything except `"plot"`,
    it will just be ignored.

  - `ggcorrmat` no longer has `corr.method` argument. To be consistent with rest
    of the functions in this package, the type of statistics should be specified
    using `type` argument. Additional, it gains a new argument
    `ggcorrplot.args`, which can be used to pass additional arguments to the
    underlying plotting function (`ggcorrplot::ggcorrplot`).

  - Both `gghistostats` and `ggdotplotstats` now use the following arguments to
    modify `geom`s corresponding to the lines and labels:
    `test.value.line.args`, `test.value.label.args`, `centrality.line.args`,
    `centrality.label.args`. This helps avoid specifying millions of arguments.

  - Removes the vestigial `ggplot_converter` function.

  - `ggpiestats` and `ggbarstats` remove the following vestigial arguments:
    `facet.wrap.name`, `bias.correct`, `bar.outline.color`. The `bar.proptest`
    and `facet.proptest` arguments were difficult to remember and confusing and
    are replaced by a common `proportion.test` argument. Additionally, the
    following arguments have all been removed and replaced by `label` argument:
    `slice.label`, `bar.label`, `data.label`. These plethora of options was a
    headache to remember.

  - `gghistostats` loses the following arguments: `fill.gradient`, `low.color`,
    `high.color`. It made no sense to add a color gradient to this plot when the
    Y-axis already displayed the information about what the bar represented.

  - `ggscatterstats` loses the following arguments: `palette` and `package`.
    Since this function requires only two colors, it didn't make much sense to
    use color palettes to specify this. They can be instead specified using
    `xfill` and `yfill`. You can always use `paletteer::paletteer_d` to get a
    vector of color values and then provide values of your choosing to `xfill`
    and `yfill`.

  - Removes sorting options in `ggbetweenstats` and `ggwithinstats` functions.
    This is something the users can easily do before entering the data in these
    functions.

MAJOR CHANGES

  - `ggcorrmat` was never supposed to work with Kendall's correlation
    coefficient but it accidentally did. This is no longer the case.

  - `{ggstatsplot}` now has a logo, thanks to Sarah! :)

  - The default `theme_ggstatsplot` changes slightly. The biggest change is that
    the title and the subtitle for plots are now aligned to the left of the
    plot. This change also forced the legend for `ggpiestats` to be displayed on
    the right side of the plot rather than at the bottom.

MINOR CHANGES

  - More models supported in `ggcoefstats`: `BBreg`, `bcplm`, `bife`, `cglm`,
    `crch`, `DirichReg`, `LORgee`, `zcpglm`, `zeroinfl`.

  - Following functions are now re-exported from `ipmisc`: `bartlett_message`,
    `normality_message`. A few other internal data wrangling functions now
    reside in `ipmisc`.

# ggstatsplot 0.2.0

BREAKING CHANGES

  - To have a more manageable length of function arguments, additional aesthetic
    specifications for any given geom can be provided via a dedicated `*.args`
    argument. For example, all aesthetic arguments for `geom_vline` can be
    provided via `vline.args`, for `geom_errorbarh` via `errorbar.args`, etc.

  - `{ggstatsplot}` continues with its conscious uncoupling that started in
    `0.1.0` release: The following functions have now been moved to
    `{statsExpressions}` package: `subtitle_meta_parametric` and
    `bf_meta_message` and follow a more logical nomenclature. For the same
    reason, `lm_effsize_ci` function is also no longer exported and lives in the
    `groupedstats` package.

MAJOR CHANGES

  - The summary caption no longer displays log-likelihood value because it tends
    to be not available for a number of regression model objects and so the
    caption was unnecessarily being skipped.

  - Supports robust and Bayes Factors for random-effects meta-analysis.

MINOR CHANGES

  - New dataset included: `bugs_wide`

  - More models supported in `ggcoefstats`: `cgam`, `cgamm`, `coxme`, `cpglm`,
    `cpglmm`, `complmrob`, `feis`, `flexsurvreg`, `glmx`, `hurdle`, `iv_robust`,
    `mixor`, `rqss`, `truncreg`, `vgam`.

  - Removed vestigial arguments from `ggcorrmat` (e.g., `exact`, `continuity`,
    etc.) and `ggpiestats` (`bf.prior`, `simulate.p.value`, `B`, etc.).

# ggstatsplot 0.1.4

BUG FIXES

  - `ggbetweenstats` and `ggwithinstats` no longer produce error with variables
    with pattern `mean` (#336).

MAJOR CHANGES

  - `pairwise_p` has been reintroduced as a number of users found it useful to
    call the function from `{ggstatsplot}` itself rather than using
    `pairwiseComparisons` package.

MINOR CHANGES

  - `ggbetweenstats` and `ggwithinstats` use `[` instead of `(` to display
    confidence intervals. Additionally, $$\mu$$ denoted sample mean, but was
    confused with population mean by some users. So these functions instead
    display $$\hat{\mu}$$.

  - More models supported in `ggcoefstats`: `bmlm`, `coeftest`

  - Adapts to the new syntax provided in `paletteer` package.

# ggstatsplot 0.1.3

MAJOR CHANGES

  - To avoid excessive arguments to function, most arguments relevant for
    `ggrepel` in `ggcoefstats` function have been removed. The users can instead
    provide all such arguments in a list to `stats.labels.args` argument.

BUG FIXES

  - `ggbetweenstats` and `ggwithinstats` no longer produce incorrect label if
    the dataframe already contains a variable named `n` (#317) or variables with
    pattern `mean` (#322).

  - `ggbetweenstats` and `ggwithinstats` mean labels respect `k` argument
    (#331).

MINOR

  - `ggcoefstats` now uses `parameters::p_value` instead of `sjstats::p_value`,
    as requested by the maintainer of that package. This might lead to
    differences in *p*-values for `lmer` models.

  - More models supported in `ggcoefstats`: `blavaan`, `bracl`, `brglm2`,
    `glmc`, `lavaan`, `nlreg`, `slm`, `wbgee`.

  - `ggcoefstats` gains `only.significant` argument to only display display
    stats labels for significant effects. This can be helpful when a large
    number of regression coefficients are to be displayed in a single plot.

# ggstatsplot 0.1.2

MINOR

  - Minor code refactoring that gets rid of the following dependencies:
    `magrittr`, `ellipsis`, `purrrlyr`.

MAJOR

  - The *p*-value label now specifies whether the *p*-value displayed in
    `ggbetweenstats` and `ggwithinstats` pairwise comparisons were adjusted or
    not for multiple comparisons.

# ggstatsplot 0.1.1

ANNOUNCEMENTS

`{ggstatsplot}` is undergoing *conscious uncoupling* whereby all the statistical
processing functions that make stats subtitles are being moved to a new package
called `{statsExpressions}`. This new package will act as a backend that handles
all things statistical processing. This **will not** affect the end users of
`{ggstatsplot}` unless you have been using the helper functions.

Additionally, multiple pairwise comparison tests are being moved to an
independent package called `pairwiseComparisons`.

This uncoupling is designed to achieve two things:

  - Make the code base of more manageable size in `{ggstatsplot}`, which will
    make package development a bit easier.

  - Make the workflow more customizable since now you can prepare your own plots
    and then use `{statsExpressions}` to display results in the plot rather than
    relying on `{ggstatsplot}` default plots which are heavily opinionated and
    not appealing to everyone.

BREAKING CHANGES

  - All helper functions `subtitle_*` and `bf_*` have been moved to the new
    `{statsExpressions}` package.

  - To be consistent with all the other `subtitle_` and `bf_` functions,
    `subtitle_contingency_tab` and `bf_contingency_tab` now use the arguments
    `x` and `y` instead of `main` and `condition`.

MAJOR CHANGES

  - Major refactoring to reduce the codesize and to rely fully on `rlang`.

  - There was confusion about what did the red point in `ggbetweenstats` and
    `ggbetweenstats` plots represents. Now the label also contains $\mu$ to
    highlight that what is being displayed is a mean value.

  - To be consistent with the rest of the functions, `ggpiestats` and
    `ggbarstats` now uses the following aliases for arguments: `x` for `main`
    and `y` for `condition`. This change is backward-compatible and should not
    pose any problems for scripts that used `main` and `condition` arguments in
    these functions.

  - Most subtitle expressions now report details about the design. In case of
    between-subjects design, this will be $n\_{obs}$, while in case of repeated
    measures design, this will be $n\_{pairs}$.

  - `pairwise.annotation` now defaults to `"p.value"` rather than `"asterisk"`
    for `ggbetweenstats` and `ggwithinstats` (and their `grouped_` variants)
    functions. This was done because the asterisk conventions are not consistent
    across various scientific disciplines.

MINOR CHANGES

  - New dataset included: `bugs_long`, for repeated measures designs with `NA`s
    present in the data.

  - `{ggstatsplot}` now uses `rcompanion` to compute Spearman's *rho* and
    Kendall's *W*. Therefore, `DescTools` is removed from dependencies.

  - `ggcoefstats` supports following objects: `bglmerMod`, `blmerMod`, `lme`,
    `mclogit`, `mmclogit`, `tobit`, `wblm`.

  - `ggcoefstats` now respects `conf.int`. It internally always defaulted to
    `conf.int = TRUE` in `broom::tidy` irrespective of what was specified by the
    user.

  - It was painfully confusing for a lot of users what exactly the asterisks in
    each facet of `ggpiestats` signified. So instead now `ggpiestats` displays
    more detailed results from a goodness of fit (gof) test. No such change is
    made for `ggbarstats` because there is no space to include more details
    above the bar.

  - Removed `conf.method` and `conf.type` arguments for `ggcoefstats`. Also,
    `p.kr` argument removed because `ggcoefstats` will begin to rely on
    `parameters` instead of `sjstats` package to compute *p*-values for some
    regression models.

# ggstatsplot 0.0.12

BUG FIXES

  - Bayes Factor in `ggwithinstats` caption, displayed by default, was
    incorrect. This has been fixed. This stemmed from a line of code which
    should have been `paired = TRUE`, but was instead `paired = FALSE`.

MAJOR CHANGES

  - The effect size measure for Kruskal-Wallis test has been changed from the
    more obscure H-based eta-squared statistic to more common and interpretable
    epsilon-squared.

MINOR CHANGES

  - `ggcoefstats` defaults to `bf.message = TRUE` to be consistent with the rest
    of the functions in the package.

  - `ggcoefstats` supports the following class of objects: `epi.2by2`, `negbin`,
    `emmGrid`, `lmrob`, `glmrob`, `glmmPQL`, `data.table`.

  - `bf_ttest` is introduced as a general function. The previously exported
    `bf_one_sample_ttest` and `bf_two_sample_ttest` become its aliases.

  - `bf_meta_message` syntax changes to adapt to updates made to `metaBMA`
    package (thanks to #259).

BREAKING CHANGES

  - The vestigial arguments `axis.text.x.margin.t`, `axis.text.x.margin.r`,
    `axis.text.x.margin.b`, `axis.text.x.margin.l` for `ggcorrmat` have been
    removed. The margins can be adjusted using `ggplot2::margin()`.

  - `gghistostats` no longer allows `data` argument to be `NULL`. This is to
    make this function's syntax consistent with rest of the functions in this
    package (none of which allow `data` to be `NULL`). This also removes
    confusion that arose for some users when `data` couldn't be `NULL` for its
    `grouped_` cousin (`grouped_gghistostats`).

  - `outlier_df` function is no longer exported since it was always meant to be
    an internal function and was accidently exported during initial release and
    was retained for a while for backward compatibility.

# ggstatsplot 0.0.11
 
BREAKING CHANGES

  - Instead of having two separate functions that dealt with repeated measures
    (`subtitle_friedman_nonparametric`) and between-subjects
    (`subtitle_kw_nonparametric`), a single function
    `subtitle_anova_nonparametric` handles both of these designs with the
    `paired` argument determining which test is run.

  - All functions that supported Bayes Factor analysis (`type = "bf"`) will only
    return BF value and the scale used. Previously, this was a mix of parametric
    statistics and BF, which was confusing and often times misleading since
    these two types of analyses relied on different tests.

  - The default for `bf.message` has been changed from `FALSE` to `TRUE`. This
    is to make the Bayes Factor analysis more visible to the user.

MAJOR CHANGES

  - `ggscatterstats` returns only plot (without any statistical details) when
    the specified model is not linear (i.e., either when `method` argument is
    not `"lm"` or when `formula` is not `y ~ x`).

NEW FEATURES

  - New functions `ggwithinstats` (and its `grouped_` variant) are introduced as
    a counterpart to `ggbetweenstats` to handle repeated measures designs.

  - For repeated measures ANOVA, `subtitle_anova_nonparametric` now returns
    confidence intervals for Kendall's *W*.

  - All functions get `return` argument that can be used to return either
    `"plot"`, `"subtitle"`, or `"caption"`. This makes it unnecessary to
    remember which subtitle function is to be used where. As a result, in the
    next release, all subtitle making functions will not be exported and are
    encouraged not be used either by other developers or by users.

  - Both `subtitle_anova_robust` and `subtitle_anova_parametric` gain a new
    argument `paired` to support repeated measures designs.

  - `ggcoefstats` can support following new model objects: `drc`, `mlm`.

  - `ggcoefstats` gains `bf.message` argument to display a caption containing
    results from Bayesian random-effects meta-analysis. It therefore gains a new
    dependency: `metaBMA`.

  - `ggpiestats` and `ggcatstats` will now display Cramer's *V* as effect size
    for one-sample proportion tests.

  - All functions gain `stat.title` argument (`NULL` by default) that can be
    used to prefix the subtitle with a string of interest. This is possibly
    useful for specifying the details of the statistical test.

MINOR CHANGES

  - `pairwise_p()` function no longer outputs `conf.low` and `conf.high` columns
    when parametric *post hoc* tests are run. This is because these values were
    accurate only when no *p*-value adjustment was carried out.

  - Instead of using the internal function `cor_test_ci`, `ggscatterstats`
    instead used `SpearmanRho` function from `DescTools` package. This was done
    to reduce number of custom internal functions used to compute CIs for
    various effect sizes. `{ggstatsplot}` therefore gains `DescTools` as a
    dependency.

  - The `sampling.plan` argument default for `ggbarstats` function has been
    changed from `"indepMulti"` to `"jointMulti"` to be consistent with its
    sister function `ggpiestats`.

# ggstatsplot 0.0.10

NEW FEATURES

  - `ggcoefstats` can support following new model objects: `rjags`.

  - New `VR_dilemma` dataset for toying around with within-subjects design.

  - `subtitle_t_onesample` supports both Cohen's *d* and Hedge's *g* as effect
    sizes and also produces their confidence intervals. Additionally,
    non-central variants of these effect sizes are also supported. Thus,
    `gghistostats` and its `grouped_` variant gets two new arguments:
    `effsize.type`, `effsize.noncentral`.

  - `ggpiestats` used to display odds ratio as effect size for paired designs
    (McNemar test). But this was only working when the analysis was a 2 x 2
    contingency table. It now instead displays Cohen's *G* as effect size, which
    generalizes to any kind of design.

MINOR CHANGES

  - The internal function `outlier_df` to add a column specifying outlier status
    of any given data point is now exported.

  - `{ggstatsplot}` previously relied on an internal function `chisq_v_ci` to
    compute confidence intervals for Cramer's *V* using bootstrapping but it was
    pretty slow. It now instead relies on `rcompanion` package to compute
    confidence intervals for *V*. `{ggstatsplot}`, therefore, gains a new
    dependency.

  - `subtitle_mann_nonparametric` and `subtitle_t_onesample` now computes effect
    size *r* and its confidence intervals as $$Z/\sqrt{N}$$ (with the help of
    `rcompanion` package), instead of using Spearman correlation.

# ggstatsplot 0.0.9

BREAKING CHANGES

  - `subtitle_t_onesample` no longer has `data` as the optional argument. This
    was done to be consistent with other subtitle helper functions.

NEW FEATURES

  - New function `ggbarstats` (and its `grouped_` variant) introduced for making
    bar charts (thanks to #78).

  - `ggcoefstats` also displays a caption with model summary when meta-analysis
    is required.

  - `gghistostats` and its `grouped_` variant has a new argument `normal.curve`
    to superpose a normal distribution curve on top of the histogram (#138).

  - `ggcoefstats` can support following new regression model objects: `brmsfit`,
    `gam`, `Gam`, `gamlss`, `mcmc`, `mjoint`, `stanreg`.

  - New function to convert plots which are not of `gg`/`ggplot` class to
    `ggplot` class objects.

  - Instead of using `effsize` to compute Cohen's *d* and Hedge's *g*,
    `{ggstatsplot}` now relies on a new (#159) internal function
    `effect_t_parametric` to compute them. This removes `effsize` from
    dependencies.

  - To be consistent with other functions in this package, both `ggbarstats` and
    `ggpiestats` gain `results.subtitle` which can be set to `FALSE` if
    statistical analysis is not required, in which case `subtitle` argument can
    be used to provide alternative subtitle.

MAJOR CHANGES

  - `ggbetweenstats` now defaults to using noncentral-*t* distribution for
    computing Cohen's *d* and Hedge's *g*. To get variants with central-*t*
    distribution, use `effsize.noncentral = FALSE`.

MINOR CHANGES

  - All `grouped_` functions had argument `title.prefix` that defaulted to
    `"Group"`. It now instead defaults to `NULL`, in which case the prefix will
    variable name for `grouping.var` argument.

  - To accommodate non-parametric tests, `subtitle_template` function can now
    work with `parameter = NULL`.

  - For `ggbetweenstats`, details contained in the subtitle for non-parametric
    test are modified. It now uses Spearman's *rho*-based effect size estimates.
    This removes `coin` from dependencies.

  - `ggbetweenstats` and its `grouped_` variant gain a new argument
    `axes.range.restrict` (which defaults to `FALSE`). This restricts `y`-axes
    limits to minimum and maximum of `y` variable. This is what these functions
    were doing by default in the past versions, which created issues for
    additional ggplot components using the `ggplot.component` argument.

  - All bayes factor related subtitle and captions replace `prior.width` with
    `r_{Cauchy}`.

  - `ggcoefstats` passes dots (`...`) to `augment` method from `broom`.

BUG FIXES

  - The helper function `bf_extractor` no longer provides option to extract
    information about posterior distribution because these details were
    incorrect. The `posterior = TRUE` details were not used anywhere in the
    package so nothing about the results changes.

  - `ggcorrmat` didn't output pair names when `output == "ci"` was used. This is
    fixed.

# ggstatsplot 0.0.8

NEW FEATURES

  - `ggcoefstats` gains `meta.analytic.effect` that can be used to carry out
    meta-analysis on regression estimates. This especially useful when a
    dataframe with regression estimates and standard error is available from
    prior analyses. The `subtitle` is prepared with the new function
    `subtitle_meta_ggcoefstats` which is also exported.

  - `ggbetweenstats`, `ggscatterstats`, `gghistostats`, and `ggdotplotstats`
    (and their `grouped_` variants) all gain a new `ggplot.component` argument.
    This argument will primarily be helpful to change the individual plots in a
    `grouped_` plot.

  - `ggcoefstats` can support following new regression model objects: `polr`,
    `survreg`, `cch`, `Arima`, `biglm`, `glmmTMB`, `coxph`, `ridgelm`, `aareg`,
    `plm`, `nlrq`, `ivreg`, `ergm`, `btergm`, `garch`, `gmm`, `lmodel2`,
    `svyolr`, `confusionMatrix`, `multinom`, `nlmerMod`, `svyglm`, `MCMCglmm`,
    `lm.beta`, `speedlm`, `fitdistr`, `mle2`, `orcutt`, `glmmadmb`.

BUG FIXES

  - `ggcoefstats` didn't work when `statistic` argument was set to `NULL`. This
    was not expected behavior. This has been fixed. Now, if `statistic` is not
    specified, only the dot-and-whiskers will be shown without any labels.

  - `subtitle_t_parametric` was producing incorrect sample size information when
    `paired = TRUE` and the data contained `NA`s. This has been fixed.

MAJOR CHANGES

  - `ggscatterstats` and its `grouped_` variant accept both character and bare
    exressions as input to arguments `label.var` and `labe.expression` (#110).

  - To be consistent with rest of the functions in the package, both Pearson's
    *r*, Spearman's *rho*, and robust percentage bend correlations also display
    information about statistic associated with these tests.

  - `ggscatterstats`, by default, showed jittered data points (because it relied
    on `position_jitter` defaults). This could be visually inaccurate and,
    therefore, `ggscatterstats` now displays points without any jitter. The user
    can introduce jitter if they wish to using `point.width.jitter` and
    `point.height.jitter` arguments. For similar reasons, for `ggbetweenstats`
    and its `grouped_` variant, `point.jitter.height` default has been changed
    from `0.1` to `0` (no vertical jitter, i.e.).

MINOR CHANGES

  - Confidence interval for Kendall's *W* is now computed using
    `stats::kruskal.test`. As a result, `PMCMRplus` removed from dependencies.

  - `ggcoefstats` gains a `caption` argument. If `caption.summary` is set to
    `TRUE`, the specified caption will be added on top of the
    `caption.summary`.

# ggstatsplot 0.0.7

BUG FIXES

  - `ggcoefstats` was showing wrong confidence intervals for `merMod` class
    objects due to a bug in the `broom.mixed` package
    (<https://github.com/bbolker/broom.mixed/issues/30#issuecomment-428385005>).
    This was fixed in `broom.mixed` and so `ggcoefstats` should no longer have
    any issues.

  - `specify_decimal_p` has been modified because it produced incorrect results
    when `k < 3` and `p.value = TRUE` (e.g., `0.002` was printed as `< 0.001`).

  - `ggpiestats` produced incorrect results if some levels of the factor had
    been filtered out prior to using this function. It now drops unused levels
    and produces correct results.

  - `gghistostats` wasn't filtering out `NA`s properly. This has been fixed.

MAJOR CHANGES

  - New function `ggdotplotstats` for creating a dot plot/chart for labelled
    numeric data.

  - All primary functions gain `conf.level` argument to control confidence level
    for effect size measures.

  - As per APA guidelines, all results show results with two decimal places.
    That is, the default value for `k` argument for all functions has been
    changed from `3` to `2`.

  - All helper functions for the `ggbetweenstats` subtitles have been renamed to
    remove `_ggbetween_` from their names as this was becoming confusing for the
    user. Some of these functions work both with the between- and
    within-subjects designs, so having `_ggbetween_` in their names made users
    suspect if they could use these functions for within-subjects designs.

  - `{ggstatsplot}` now depends on `R 3.5.0`. This is because some of its
    dependencies require 3.5.0 to work (e.g., `broom.mixed`).

  - All `theme_` functions are now exported (`theme_pie()`, `theme_corrmat()`).

  - `ggbetweenstats` now supports multiple pairwise comparison tests
    (parametric, nonparametric, and robust variants). It gains a new dependency
    `ggsignif`.

  - `ggbetweenstats` now supports eta-squared and omega-squared effect sizes for
    anova models. This function gains a new argument `partial`.

  - Following functions are now reexported from the `groupedstats` package to
    avoid repeating the same code in two packages: `specify_decimal_p`,
    `signif_column`, `lm_effsize_ci`, and `set_cwd`. Therefore, `groupedstats`
    is now added as a dependency.

  - `gghistostats` can now show both counts and proportions information on the
    same plot when `bar.measure` argument is set to `"mix"`.

  - `ggcoefstats` works with tidy dataframes.

  - The helper function `untable` has been deprecated in light of
    `tidyr::uncount`, which does exactly what `untable` was doing. The author
    wasn't aware of this function when `untable` was written.

  - All vignettes have been removed from `CRAN` to reduce the size of the
    package. They are now available on the package website:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/>.

  - `subtitle_t_robust` function can now handle dependent samples and gains
    `paired` argument.

  - A number of tidyverse operators are now reexported by `{ggstatsplot}`:
    `%>%`, `%<>%`, `%$%`.

MINOR CHANGES

  - `ggscatterstats`, `ggpiestats`, and their `grouped_` variant support bayes
    factor tests and gain new arguments relevant to this test.

  - Effect size and their confidence intervals now available for Kruskal-Wallis
    test.

  - Minor stylistic changes to how symbols for partial-eta-/omega-squared were
    being displayed in subtitles.

  - `ggbetweenstats` supports bayes factor tests for anova designs.

  - `ggpiestats` (and its `grouped_` version) gain `slice.label` argument that
    decides what information needs to be displayed as a label on the slices of
    the pie chart: `"percentage"` (which has been the default thus far),
    `"counts"`, or `"both"`.

  - `ggcorrmat` can work with `cor.vars = NULL`. In such case, **all** numeric
    variables from the provided dataframe will be used for computing the
    correlation matrix.

  - Given the constant changes to the default behavior of functions, the
    lifecycle badge has been changed from `stable` to `maturing`.

  - When the number of colors needed by a function exceeds the number of colors
    contained in a given palette, informative message is displayed to the user
    (with the new internal function `palette_message()`).

  - Several users had requested an easier way to turn off subtitles with results
    from tests (which was already implemented in `ggscatterstats` and
    `gghistostats` with the argument `results.subtitle`), so `ggbetweenstats`
    also gains two new arguments to do this: `results.subtitle` and `subtitle`.

  - New dataset added: `iris_long`.

  - More tests added and the code coverage has now jumped to over 75%.

  - To avoid code repetition, there is a now a function that produces a generic
    message any time confidence intervals for effect size estimate are computed
    using bootstrapping.

# ggstatsplot 0.0.6

MAJOR CHANGES

  - The package now exports all functions used to create text expressions with
    results. This makes it easy for people to use these results in their own
    plots at any location they want (and not just in `subtitle`, the current
    default for `{ggstatsplot}`).

  - `ggcorrmat` gains `p.adjust.method` argument which allows *p*-values for
    correlations to be corrected for multiple comparisons.

  - `ggscatterstats` gains `label.var` and `label.expression` arguments to
    attach labels to points.

  - `gghistostats` now defaults to not showing (redundant) color gradient
    (`fill.gradient = FALSE`) and shows both `"count"` and `"proportion"` data.
    It also gains a new argument `bar.fill` that can be used to fill bars with a
    uniform color.

  - `ggbetweenstats`, `ggcoefstats`, `ggcorrmat`, `ggscatterstats`, and
    `ggpiestats` now support all palettes contained in the `paletteer` package.
    This helps avoid situations where people had large number of groups (> 12)
    and there were not enough colors in any of the `RColorBrewer` palettes.

  - `ggbetweenstats` gains `bf.message` argument to display bayes factors in
    favor of the null (currently works only for parametric *t*-test).

  - `gghistostats` function no longer has `line.labeller.y` argument; this
    position is automatically determined now.

BREAKING CHANGES

  - `legend.title.margin` function has been deprecated since `ggplot2 3.0.0` has
    improved on the margin issues from previous versions. All functions that
    wrapped around this function now lose the relevant arguments
    (`legend.title.margin`, `t.margin`, `b.margin`).

  - The argument `ggstatsplot.theme` has been changed to `ggstatsplot.layer` for
    `ggcorrmat` function to be consistent across functions.

  - For consistency, `conf.level` and `conf.type` arguments for `ggbetweenstats`
    have been deprecated. No other function in the package allowed changing
    confidence interval or their type for effect size estimation. These
    arguments were relevant only for `robust` tests anyway.

  - `ggocorrmat` argument `type` has been changed to `matrix.type` because for
    all other functions `type` argument specifies the type of the test, while
    for this function it specified the display of the visualization matrix. This
    will make the syntax more consistent across functions.

  - `ggscatterstats` gains new arguments to specify aesthetics for geom point
    (`point.color`, `point.size`, `point.alpha`). To be consistent with this
    naming schema, the `width.jitter` and `height.jitter` arguments have been
    renamed to `point.width.jitter` and `point.height.jitter`, resp.

MINOR CHANGES

  - `gghistostats`: To be compatible with `JASP`, natural logarithm of Bayes
    Factors is displayed, and not base 10 logarithm.

  - `ggscatterstats` gains `method` and `formula` arguments to modify smoothing
    functions.

  - `ggcorrmat` can now show `robust` correlation coefficients in the matrix
    plot.

  - For `gghistostats`, `binwidth` value, if not specified, is computed with
    `(max-min)/sqrt(n)`. This is basically to get rid of the warnings ggplot2
    produces. Thanks to Chuck Powell's PR (#43).

  - `ggcoefstats` gains a new argument `partial` and can display eta-squared and
    omega-squared effect sizes for anovas, in addition to the prior partial
    variants of these effect sizes.

  - `ggpiestats` gains `perc.k` argument to show desired number of decimal
    places in percentage labels.

BUG FIXES

  - `grouped_ggpiestats` wasn't working when only `main` variable was provided
    with `counts` data. Fixed that.

# ggstatsplot 0.0.5

MAJOR CHANGES

  - For the sake of consistency, `theme_mprl` is now called `theme_ggstatsplot`.
    The `theme_mprl` function will still be around and will **not** be
    deprecated, so feel free to use either or both of them since they are
    identical.

  - `ggcoefstats` no longer has arguments `effects` and `ran_params` because
    only fixed effects are shown for mixed-effects models.

  - `ggpiestats` can now handle within-subjects designs (McNemar test results
    will be displayed).

BUG FIXES

  - `ggbetweenstats` was producing wrong axes labels when `sample.size.label`
    was set to `TRUE` and user had reordered factor levels before using this
    function. The new version fixes this.

  - `ggcoefstats` wasn't producing partial omega-squared for `aovlist` objects.
    Fixed that with new version of `sjstats`.

MINOR CHANGES

  - Removed the trailing comma from the robust correlation analyses.

  - `gghistostats` has a new argument to remove color fill gradient.

  - `ggbetweenstats` takes new argument `mean.ci` to show confidence intervals
    for the mean values.

  - For `lmer` models, *p*-values are now computed using `sjstats::p_value`.
    This removes `lmerTest` package from dependencies.

  - `sjstats` no longer suggests `apaTables` package to compute confidence
    intervals for partial eta- and omega-squared. Therefore, `apaTables` and
    `MBESS` are removed from dependencies.

  - `ggscatterstats` supports `densigram` with the development version of
    `ggExtra`. It additionally gains few extra arguments to change aesthetics of
    marginals (alpha, size, etc.).

# ggstatsplot 0.0.4

MAJOR CHANGES

  - New function: `ggcoefstats` for displaying model coefficients.

  - All functions now have `ggtheme` argument that can be used to change the
    default theme, which has now been changed from `theme_grey()` to
    `theme_bw()`.

  - The robust correlation is no longer `MASS::rlm`, but percentage bend
    correlation, as implemented in `WRS2::pbcor`. This was done to be consistent
    across different functions. `ggcorrmat` also uses percentage bend
    correlation as the robust correlation measure. This also means that
    `{ggstatsplot}` no longer imports `MASS` and `sfsmisc`.

  - The `data` argument is no longer `NULL` for all functions, except
    `gghistostats`. In other words, the user **must** provide a dataframe from
    which variables or formulas should be selected.

  - All subtitles containing results now also show sample size information
    (*n*). To adjust for the inflated length of the subtitle, the default
    subtitle text size has been changed from `12` to `11`.

MINOR CHANGES

  - Switched back to Shapiro-Wilk test of normality to remove `nortest` from
    imports.

  - `ggbetweenstats` and `ggpiestats` now display sample sizes for each level of
    the groping factor by default. This behavior can be turned off by setting
    `sample.size.label` to `FALSE`.

  - Three new datasets added: `Titanic_full`, `movies_wide`, `movies_long`.

  - Added confidence interval for effect size for robust ANOVA.

  - The 95% CI for Cramer'V computed using `boot::boot`. Therefore, the package
    no longer imports `DescTools`.

  - To be consistent across correlations covered, all correlations now show
    estimates for correlation coefficients, confidence intervals for the
    estimate, and *p*-values. Therefore, *t*-values and regression coefficients
    are no longer displayed for Pearson's *r*.

  - The `legend.title.margin` arguments for `gghistostats` and `ggcorrmat` now
    default to `FALSE`, since `ggplot2 3.0.0` has better legend title margins.

  - `ggpiestats` now sorts the summary dataframes not by percentages but by the
    levels of `main` variable. This was done to have the same legends across
    different levels of a grouping variable in `grouped_ggpiestats`.

  - To remove cluttered display of results in the subtitle, `ggpiestats` no
    longer shows titles for the tests run (these were "Proportion test" and
    "Chi-Square test"). From the pie charts, it should be obvious to the user or
    reader what test was run.

  - `gghistostats` also allows running robust version of one-sample test now
    (One-sample percentile bootstrap).

# ggstatsplot 0.0.3

NEW FEATURES

  - The `ggbetweenstats` function can now show notched box plots. Two new
    arguments `notch` and `notchwidth` control its behavior. The defaults are
    still standard box plots.

  - Removed warnings that were appearing when `outlier.label` argument was of
    `character` type.

  - The default color palette used for all plots is colorblind friendly.

  - `gghistostats` supports `proportion` and `density` as a value measure for
    bar heights to show proportions and density. New argument `bar.measure`
    controls this behavior.

  - `grouped_` variants of functions `ggcorrmat`, `ggscatterstats`,
    `ggbetweenstats`, and `ggpiestats` introduced to create multiple plots for
    different levels of a grouping variable.

MAJOR CHANGES

  - To be internally consistent, all functions in `{ggstatsplot}` use the
    spelling `color`, rather than `colour` in some functions, while `color` in
    others.

  - Removed the redundant argument `binwidth.adjust` from `gghistostats`
    function. This argument was relevant for the first avatar of this function,
    but is no longer playing any role.

  - To be internally consistent, the argument `lab_col` and `lab_size` in
    `ggcorrmat` have been changed to `lab.col` and `lab.size`, respectively.

MINOR CHANGES

  - Added a new argument to `ggstatsplot.theme` function to control if
    `ggstatsplot::theme_mprl` is to be overlaid on top of the selected `ggtheme`
    (ggplot2 theme, i.e.).

  - Two new arguments added to `gghistostats` to allow user to change colorbar
    gradient. Defaults are colorblind friendly.

  - Both `gghistostats` and `ggcorrmat` have a new argument
    `legend.title.margin` to control margin adjustment between the title and the
    colorbar.

  - The vertical lines denoting test values and centrality parameters can be
    tagged with text labels with a new argument `line.labeller` in
    `gghistostats` function.

BUG FIXES

  - The `centrality.para` argument for `ggscatterstats` was not working
    properly. Choosing `"median"` didn't show median, but the mean. This is
    fixed now.

# ggstatsplot 0.0.2

NEW FEATURES

  - Bayesian test added to `gghistostats` and two new arguments to also display
    a vertical line for `test.value` argument.

  - Vignette added for `gghistostats`.

  - Added new function `grouped_gghistostats` to facilitate applying
    `gghistostats` for multiple levels of a grouping factor.

  - `ggbetweenstats` has a new argument `outlier.coef` to adjust threshold used
    to detect outliers. Removed bug from the same function when `outlier.label`
    argument is of factor/character type.

MAJOR CHANGES

  - Functions `signif_column` and `grouped_proptest` are now deprecated. They
    were exported in the first release by mistake.

  - Function `gghistostats` no longer displays both density and count since the
    density information was redundant. The `density.plot` argument has also been
    deprecated.

  - `ggscatterstats` argument `intercept` has now been changed to
    `centrality.para`. This was due to possible confusion about interpretation
    of these lines; they show central tendency measures and not intercept for
    the linear model. Thus the change.

  - The default for `effsize.type = "biased"` effect size for `ggbetweenstats`
    in case of ANOVA is **partial** omega-squared, and not omega-squared.
    Additionally, both partial eta- and omega-squared are not computed using
    bootstrapping with (default) 100 bootstrap samples.

MINOR CHANGES

  - More examples added to the `README` document.

  - 95% confidence intervals for Spearman's rho are now computed using `broom`
    package. `RVAideMemoire` package is thus removed from dependencies.

  - 95% confidence intervals for partial eta- and omega-squared for
    `ggbetweenstats` function are now computed using `sjstats` package, which
    allows bootstrapping. `apaTables` and `userfriendlyscience` packages are
    thus removed from dependencies.

# ggstatsplot 0.0.1

  - First release of the package.

# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
## Test environments

* local Windows install, R 4.1.2

* ubuntu 14.04 (on gitub-actions-ci), R 4.1.2

* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

  - Fixes breakage due to `{statsExpressions}` update.

## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and
dev versions of this package.

 * We saw 0 new problems

 * We failed to check 0 packages


---
title: "Visualizations with statistical details: The 'ggstatsplot' approach"
tags:
  - R
  - parametric statistics
  - nonparametric statistics
  - robust statistics
  - Bayesian statistics
  - ggplot2
  - ggplot2-extension
authors:
  - name: Indrajeet Patil
    orcid: 0000-0003-1995-6531
    affiliation: 1
affiliations:
  - name: Center for Humans and Machines, Max Planck Institute for Human Development, Berlin, Germany
    index: 1
bibliography: paper.bib
date: "2021-05-25"
---



# Summary

Graphical displays can reveal problems in a statistical model that might not be
apparent from purely numerical summaries. Such visualizations can also be
helpful for the reader to evaluate the validity of a model if it is reported in
a scholarly publication or report. But, given the onerous costs involved,
researchers often avoid preparing information-rich graphics and exploring
several statistical approaches or tests available. The `ggstatsplot` package in
the R programming language [@base2021] provides a one-line syntax to enrich
`ggplot2`-based visualizations with the results from statistical analysis
embedded in the visualization itself. In doing so, the package helps researchers
adopt a rigorous, reliable, and robust data exploratory and reporting workflow.

# Statement of Need

In a typical data analysis workflow, data visualization and statistical modeling
are two different phases: visualization informs modeling, and in turn, modeling
can suggest a different visualization method, and so on and so forth
[@wickham2016r]. The central idea of `ggstatsplot` is simple: combine these two
phases into one in the form of an informative graphic with statistical details.

Before discussing benefits of this approach, we will show an example (Figure
1).


```r
library(palmerpenguins) # for 'penguins' dataset
library(ggstatsplot)

ggbetweenstats(penguins, species, body_mass_g)
```

\begin{figure}
\includegraphics[width=1\linewidth]{paper_files/figure-latex/penguins-1} \caption{Example plot from the `ggstatsplot` package illustrating its philosophy of juxtaposing informative visualizations with details from statistical analysis. To see all supported plots and statistical analyses, see the package website: \url{https://indrajeetpatil.github.io/ggstatsplot/}}\label{fig:penguins}
\end{figure}

As can be seen, with a single line of code, the function produces details about
descriptive statistics, inferential statistics, effect size estimate and its
uncertainty, pairwise comparisons, Bayesian hypothesis testing, Bayesian
posterior estimate and its uncertainty. Moreover, these details are juxtaposed
with informative and well-labeled visualizations. The defaults are designed to
follow best practices in both data visualization [@cleveland1985;
@grant2018data; @healy2018data; @tufte2001; @wilke2019fundamentals] and
(frequentist/Bayesian) statistical reporting [@apa2019; @van2020jasp]. Without
`ggstatsplot`, getting these statistical details and customizing a plot would
require significant amount of time and effort. In other words, this package
removes the trade-off often faced by researchers between ease and thoroughness
of data exploration and further cements good data exploration habits.

Internally, data cleaning is carried out using the `tidyverse` [@Wickham2019],
while statistical analysis is carried out via the `statsExpressions`
[@Patil2021] and `easystats` [@Ben-Shachar2020; @Lüdecke2020parameters;
@Lüdecke2020performance;
@Lüdecke2019; @Makowski2019; @Makowski2020] packages. All visualizations are
constructed using the grammar of graphics framework [@Wilkinson2012], as
implemented in the `ggplot2` package [@Wickham2016].

# Benefits

In summary, the benefits of `ggstatsplot`'s approach are the following. It:

a. produces charts displaying both raw data, and numerical plus graphical
   summary indices,

b. avoids errors in and increases reproducibility of statistical reporting,

c. highlights the importance of the effect by providing effect size measures by
   default,

d. provides an easy way to evaluate *absence* of an effect using Bayes factors,

e. encourages researchers and readers to evaluate statistical assumptions of a
model in the context of the underlying data (Figure 2),

f. is easy and simple enough that someone with little to no coding experience
   can use it without making an error and may even encourage beginners to
   programmatically analyze data, instead of using GUI software.

\begin{figure}
\includegraphics[width=1\linewidth]{reporting} \caption{Comparing the 'Standard' approach of reporting statistical analysis in a publication/report with the 'ggstatsplot' approach of reporting the same analysis next to an informative graphic. Note that the results described in the 'Standard' approach are about the 'Dinosaur' dataset plotted on the right. Without the accompanying visualization, it is hard to evaluate the validity of the results. The ideal reporting practice will be a hybrid of these two approaches where the plot contains both the visual and numerical summaries about a statistical model, while the narrative provides interpretative context for the reported statistics.}\label{fig:reporting}
\end{figure}

# Future Scope

This package is an ambitious, ongoing, and long-term project. It currently
supports common statistical tests (parametric, non-parametric, robust, or
Bayesian *t*-test, one-way ANOVA, contingency table analysis, correlation
analysis, meta-analysis, regression analyses, etc.) and corresponding
visualizations (box/violin plot, scatter plot, dot-and-whisker plot, pie chart,
bar chart, etc.). It will continue expanding to support an increasing
collection of statistical analyses and visualizations.

# Licensing and Availability

`ggstatsplot` is licensed under the GNU General Public License (v3.0), with all
source code stored at [GitHub](https://github.com/IndrajeetPatil/ggstatsplot/).
In the spirit of honest and open science, requests and suggestions for fixes,
feature updates, as well as general questions and concerns are encouraged via
direct interaction with contributors and developers by filing an
[issue](https://github.com/IndrajeetPatil/ggstatsplot/issues) while respecting
[*Contribution
Guidelines*](https://indrajeetpatil.github.io/ggstatsplot/CONTRIBUTING.html).

# Acknowledgements

I would like to acknowledge the support of Mina Cikara, Fiery Cushman, and Iyad
Rahwan during the development of this project. `ggstatsplot` relies heavily on
the [`easystats`](https://github.com/easystats/easystats) ecosystem, a
collaborative project created to facilitate the usage of `R` for statistical
analyses. Thus, I would like to thank the
[members](https://github.com/orgs/easystats/people) of `easystats` as well as
the users. I would additionally like to thank the contributors to `ggstatsplot`
for reporting bugs, providing helpful feedback, or helping with enhancements.

# References

# checking labels and data from plot

    Code
      list(pb$data[[1]], pb$data[[2]], pb$data[[4]], pb$data[[5]])
    Output
      [[1]]
            colour         x       y PANEL group      xmin     xmax    ymax shape
      1  #E7298AFF 4.0403418 0.01550     1     4 3.7125931 4.330736 0.01550    19
      2  #E7298AFF 3.9732189 0.00029     1     4 3.7042826 4.247493 0.00029    19
      3  #D95F02FF 2.0136040 0.42300     1     2 1.7174850 2.360507 0.42300    19
      4  #1B9E77FF 0.9797653 0.07000     1     1 0.6977975 1.300184 0.07000    19
      5  #D95F02FF 1.9716297 0.09820     1     2 1.6977489 2.333531 0.09820    19
      6  #D95F02FF 2.0276917 0.11500     1     2 1.7433664 2.299266 0.11500    19
      7  #D95F02FF 2.0022683 0.00550     1     2 1.7457266 2.357530 0.00550    19
      8  #D95F02FF 2.0590021 0.00640     1     2 1.7127649 2.362706 0.00640    19
      9  #E7298AFF 3.9907166 0.00100     1     4 3.6567191 4.324474 0.00100    19
      10 #E7298AFF 3.9594016 0.00660     1     4 3.6730650 4.272324 0.00660    19
      11 #E7298AFF 4.0597778 0.00014     1     4 3.7332554 4.283995 0.00014    19
      12 #1B9E77FF 0.9529747 0.01080     1     1 0.6920141 1.365928 0.01080    19
      13 #D95F02FF 1.9924288 0.01230     1     2 1.7405638 2.308408 0.01230    19
      14 #E7298AFF 4.0313436 0.00630     1     4 3.6466997 4.365342 0.00630    19
      15 #7570B3FF 3.0330813 0.00030     1     3 2.7064238 3.340464 0.00030    19
      16 #D95F02FF 2.0360193 0.41900     1     2 1.7247221 2.297002 0.41900    19
      17 #E7298AFF 4.0489722 0.00350     1     4 3.7349390 4.296124 0.00350    19
      18 #E7298AFF 3.9421521 0.11500     1     4 3.7030347 4.302837 0.11500    19
      19 #1B9E77FF 0.9387469 0.02560     1     1 0.7215898 1.262465 0.02560    19
      20 #E7298AFF 4.0033265 0.00500     1     4 3.6785431 4.352286 0.00500    19
      21 #D95F02FF 1.9982621 0.01227     1     2 1.6752201 2.327097 0.01227    19
      22 #E7298AFF 4.0312523 0.17900     1     4 3.6490267 4.255737 0.17900    19
      23 #D95F02FF 2.0604009 0.00100     1     2 1.7557500 2.239612 0.00100    19
      24 #D95F02FF 1.9732552 0.00040     1     2 1.6868920 2.295826 0.00040    19
      25 #7570B3FF 3.0428055 0.00025     1     3 2.6790085 3.245088 0.00025    19
      26 #1B9E77FF 0.9745406 0.01250     1     1 0.6893698 1.256894 0.01250    19
      27 #D95F02FF 1.9685262 0.01210     1     2 1.7104713 2.336556 0.01210    19
      28 #D95F02FF 2.0211226 0.17500     1     2 1.7449878 2.237119 0.17500    19
      29 #E7298AFF 3.9769133 0.44000     1     4 3.6860434 4.272015 0.44000    19
      30 #E7298AFF 4.0260788 0.18000     1     4 3.6397980 4.340925 0.18000    19
      31 #D95F02FF 2.0546893 0.00190     1     2 1.6410769 2.354352 0.00190    19
      32 #E7298AFF 4.0409756 0.02000     1     4 3.6912382 4.357023 0.02000    19
      33 #7570B3FF 3.0106334 0.00120     1     3 2.6905110 3.281033 0.00120    19
      34 #D95F02FF 1.9711859 0.00118     1     2 1.6419119 2.234239 0.00118    19
      35 #D95F02FF 2.0609246 0.00570     1     2 1.6613101 2.341227 0.00570    19
      36 #D95F02FF 1.9822556 0.00400     1     2 1.6578058 2.332219 0.00400    19
      37 #E7298AFF 3.9732701 0.18000     1     4 3.6431099 4.256042 0.18000    19
      38 #7570B3FF 3.0055386 0.02500     1     3 2.6700479 3.265966 0.02500    19
      39 #D95F02FF 1.9382748 0.16900     1     2 1.6781337 2.330255 0.16900    19
      40 #E7298AFF 4.0071521 0.00260     1     4 3.6472775 4.360093 0.00260    19
      41 #E7298AFF 3.9360660 0.00250     1     4 3.7659396 4.337334 0.00250    19
      42 #1B9E77FF 1.0256881 0.01750     1     1 0.7064548 1.338691 0.01750    19
      43 #1B9E77FF 1.0523226 0.04450     1     1 0.6468822 1.336682 0.04450    19
      44 #1B9E77FF 0.9753935 0.05040     1     1 0.7395184 1.274638 0.05040    19
         size fill alpha stroke
      1     3   NA   0.4      0
      2     3   NA   0.4      0
      3     3   NA   0.4      0
      4     3   NA   0.4      0
      5     3   NA   0.4      0
      6     3   NA   0.4      0
      7     3   NA   0.4      0
      8     3   NA   0.4      0
      9     3   NA   0.4      0
      10    3   NA   0.4      0
      11    3   NA   0.4      0
      12    3   NA   0.4      0
      13    3   NA   0.4      0
      14    3   NA   0.4      0
      15    3   NA   0.4      0
      16    3   NA   0.4      0
      17    3   NA   0.4      0
      18    3   NA   0.4      0
      19    3   NA   0.4      0
      20    3   NA   0.4      0
      21    3   NA   0.4      0
      22    3   NA   0.4      0
      23    3   NA   0.4      0
      24    3   NA   0.4      0
      25    3   NA   0.4      0
      26    3   NA   0.4      0
      27    3   NA   0.4      0
      28    3   NA   0.4      0
      29    3   NA   0.4      0
      30    3   NA   0.4      0
      31    3   NA   0.4      0
      32    3   NA   0.4      0
      33    3   NA   0.4      0
      34    3   NA   0.4      0
      35    3   NA   0.4      0
      36    3   NA   0.4      0
      37    3   NA   0.4      0
      38    3   NA   0.4      0
      39    3   NA   0.4      0
      40    3   NA   0.4      0
      41    3   NA   0.4      0
      42    3   NA   0.4      0
      43    3   NA   0.4      0
      44    3   NA   0.4      0
      
      [[2]]
           ymin    lower   middle upper  ymax            outliers notchupper
      1 0.01080 0.017500 0.044500 0.070 0.070        0.325, 0.157 0.07215000
      2 0.00040 0.005125 0.012285 0.236 0.423 4.603, 0.655, 5.712 0.09385285
      3 0.00025 0.000300 0.001200 0.025 0.025               0.081 0.01865296
      4 0.00014 0.002600 0.006600 0.179 0.440                1.32 0.07419759
         notchlower x flipped_aes PANEL group ymin_final ymax_final xmin xmax xid
      1  0.01685000 1       FALSE     1     1    0.01080      0.325 0.85 1.15   1
      2 -0.06928285 2       FALSE     1     2    0.00040      5.712 1.85 2.15   2
      3 -0.01625296 3       FALSE     1     3    0.00025      0.081 2.85 3.15   3
      4 -0.06099759 4       FALSE     1     4    0.00014      1.320 3.85 4.15   4
        newx new_width weight colour  fill size alpha shape linetype
      1    1       0.3      1 grey20 white  0.5   0.2    19    solid
      2    2       0.3      1 grey20 white  0.5   0.2    19    solid
      3    3       0.3      1 grey20 white  0.5   0.2    19    solid
      4    4       0.3      1 grey20 white  0.5   0.2    19    solid
      
      [[3]]
        x     y            label PANEL group    colour  fill size angle alpha family
      1 2 4.603   Asian elephant     1     2 darkgreen white 3.88     0    NA       
      2 2 0.655            Horse     1     2 darkgreen white 3.88     0    NA       
      3 1 0.325        Gray seal     1     1 darkgreen white 3.88     0    NA       
      4 4 1.320            Human     1     4 darkgreen white 3.88     0    NA       
      5 2 5.712 African elephant     1     2 darkgreen white 3.88     0    NA       
      6 1 0.157           Jaguar     1     1 darkgreen white 3.88     0    NA       
      7 3 0.081  Giant armadillo     1     3 darkgreen white 3.88     0    NA       
        fontface lineheight hjust vjust point.size segment.linetype segment.size
      1        1        1.2   0.5   0.5          1                1          0.5
      2        1        1.2   0.5   0.5          1                1          0.5
      3        1        1.2   0.5   0.5          1                1          0.5
      4        1        1.2   0.5   0.5          1                1          0.5
      5        1        1.2   0.5   0.5          1                1          0.5
      6        1        1.2   0.5   0.5          1                1          0.5
      7        1        1.2   0.5   0.5          1                1          0.5
        segment.curvature segment.angle segment.ncp segment.shape segment.square
      1                 0            90           1           0.5           TRUE
      2                 0            90           1           0.5           TRUE
      3                 0            90           1           0.5           TRUE
      4                 0            90           1           0.5           TRUE
      5                 0            90           1           0.5           TRUE
      6                 0            90           1           0.5           TRUE
      7                 0            90           1           0.5           TRUE
        segment.squareShape segment.inflect segment.debug
      1                   1           FALSE         FALSE
      2                   1           FALSE         FALSE
      3                   1           FALSE         FALSE
      4                   1           FALSE         FALSE
      5                   1           FALSE         FALSE
      6                   1           FALSE         FALSE
      7                   1           FALSE         FALSE
      
      [[4]]
        x          y PANEL group shape  colour size fill alpha stroke
      1 1 0.07925556     1     1    19 darkred    5   NA    NA    0.5
      2 2 0.62159750     1     2    19 darkred    5   NA    NA    0.5
      3 3 0.02155000     1     3    19 darkred    5   NA    NA    0.5
      4 4 0.14573118     1     4    19 darkred    5   NA    NA    0.5
      

---

    Code
      within(pb$plot$labels, rm(subtitle, caption))
    Output
      $x
      [1] "vore"
      
      $y
      [1] "brain weight"
      
      $colour
      [1] "vore"
      
      $title
      [1] "mammalian sleep"
      
      $label
      [1] "outlier.label"
      
      $alt
      [1] ""
      

# checking mean labels are working

    Code
      list(pb$data[[1]], pb$data[[2]], pb$data[[4]], pb$data[[5]])
    Output
      [[1]]
            colour         x     y PANEL group      xmin     xmax  ymax shape size
      1  #D95F02FF 1.9906121 2.620     1     2 1.6885979 2.367611 2.620    19    3
      2  #D95F02FF 2.0003784 2.875     1     2 1.7007673 2.250466 2.875    19    3
      3  #1B9E77FF 1.0465125 2.320     1     1 0.6514721 1.281202 2.320    19    3
      4  #D95F02FF 1.9467189 3.215     1     2 1.7041109 2.317664 3.215    19    3
      5  #7570B3FF 3.0286507 3.440     1     3 2.6373856 3.223981 3.440    19    3
      6  #D95F02FF 2.0010960 3.460     1     2 1.6454085 2.308193 3.460    19    3
      7  #7570B3FF 2.9824592 3.570     1     3 2.6670947 3.277994 3.570    19    3
      8  #1B9E77FF 1.0219558 3.190     1     1 0.6941384 1.338097 3.190    19    3
      9  #1B9E77FF 1.0228435 3.150     1     1 0.6355074 1.354052 3.150    19    3
      10 #D95F02FF 1.9334212 3.440     1     2 1.6890293 2.228340 3.440    19    3
      11 #D95F02FF 2.0446931 3.440     1     2 1.6861225 2.348038 3.440    19    3
      12 #7570B3FF 2.9281203 4.070     1     3 2.6846597 3.363751 4.070    19    3
      13 #7570B3FF 3.0577960 3.730     1     3 2.6675857 3.360505 3.730    19    3
      14 #7570B3FF 3.0198403 3.780     1     3 2.7448462 3.246743 3.780    19    3
      15 #1B9E77FF 1.0772748 2.200     1     1 0.6771053 1.312507 2.200    19    3
      16 #1B9E77FF 0.9594168 1.615     1     1 0.7409464 1.374648 1.615    19    3
      17 #1B9E77FF 1.0617654 1.835     1     1 0.7297045 1.321450 1.835    19    3
      18 #1B9E77FF 1.0010493 2.465     1     1 0.6437524 1.266164 2.465    19    3
      19 #7570B3FF 2.9877796 3.520     1     3 2.6646107 3.249218 3.520    19    3
      20 #7570B3FF 3.0381273 3.435     1     3 2.6908391 3.369801 3.435    19    3
      21 #7570B3FF 2.9937107 3.840     1     3 2.7197371 3.274652 3.840    19    3
      22 #7570B3FF 3.0415521 3.845     1     3 2.7437087 3.313843 3.845    19    3
      23 #1B9E77FF 0.9486889 1.935     1     1 0.7618041 1.273118 1.935    19    3
      24 #1B9E77FF 1.0775687 2.140     1     1 0.6239478 1.287262 2.140    19    3
      25 #1B9E77FF 0.9875719 1.513     1     1 0.7659043 1.278687 1.513    19    3
      26 #7570B3FF 3.0707233 3.170     1     3 2.7309507 3.268260 3.170    19    3
      27 #D95F02FF 1.9905971 2.770     1     2 1.6916866 2.259688 2.770    19    3
      28 #7570B3FF 2.9902851 3.570     1     3 2.6971442 3.357752 3.570    19    3
      29 #1B9E77FF 0.9719088 2.780     1     1 0.7659033 1.305297 2.780    19    3
         fill alpha stroke
      1    NA   0.4      0
      2    NA   0.4      0
      3    NA   0.4      0
      4    NA   0.4      0
      5    NA   0.4      0
      6    NA   0.4      0
      7    NA   0.4      0
      8    NA   0.4      0
      9    NA   0.4      0
      10   NA   0.4      0
      11   NA   0.4      0
      12   NA   0.4      0
      13   NA   0.4      0
      14   NA   0.4      0
      15   NA   0.4      0
      16   NA   0.4      0
      17   NA   0.4      0
      18   NA   0.4      0
      19   NA   0.4      0
      20   NA   0.4      0
      21   NA   0.4      0
      22   NA   0.4      0
      23   NA   0.4      0
      24   NA   0.4      0
      25   NA   0.4      0
      26   NA   0.4      0
      27   NA   0.4      0
      28   NA   0.4      0
      29   NA   0.4      0
      
      [[2]]
         ymin  lower middle   upper ymax            outliers notchupper notchlower x
      1 1.513 1.8850  2.200 2.62250 3.19                       2.551336   1.848664 1
      2 2.620 2.8225  3.215 3.44000 3.46                       3.583761   2.846239 2
      3 3.170 3.5325  3.755 4.01375 4.07 5.250, 5.424, 5.345   3.958219   3.551781 3
        flipped_aes PANEL group ymin_final ymax_final xmin xmax xid newx new_width
      1       FALSE     1     1      1.513      3.190 0.85 1.15   1    1       0.3
      2       FALSE     1     2      2.620      3.460 1.85 2.15   2    2       0.3
      3       FALSE     1     3      3.170      5.424 2.85 3.15   3    3       0.3
        weight colour  fill size alpha shape linetype
      1      1 grey20 white  0.5   0.2    19    solid
      2      1 grey20 white  0.5   0.2    19    solid
      3      1 grey20 white  0.5   0.2    19    solid
      
      [[3]]
        x     y               label PANEL group colour  fill size angle alpha family
      1 3 5.250  Cadillac Fleetwood     1     1  black white    3     0    NA       
      2 3 5.424 Lincoln Continental     1     1  black white    3     0    NA       
      3 3 5.345   Chrysler Imperial     1     1  black white    3     0    NA       
        fontface lineheight hjust vjust point.size segment.linetype segment.size
      1        1        1.2   0.5   0.5          1                1          0.5
      2        1        1.2   0.5   0.5          1                1          0.5
      3        1        1.2   0.5   0.5          1                1          0.5
        segment.curvature segment.angle segment.ncp segment.shape segment.square
      1                 0            90           1           0.5           TRUE
      2                 0            90           1           0.5           TRUE
      3                 0            90           1           0.5           TRUE
        segment.squareShape segment.inflect segment.debug
      1                   1           FALSE         FALSE
      2                   1           FALSE         FALSE
      3                   1           FALSE         FALSE
      
      [[4]]
        x     y PANEL group shape  colour size fill alpha stroke
      1 1 2.200     1     1    19 darkred    5   NA    NA    0.5
      2 2 3.215     1     2    19 darkred    5   NA    NA    0.5
      3 3 3.755     1     3    19 darkred    5   NA    NA    0.5
      

---

    Code
      pb$plot$labels
    Output
      $x
      [1] "cyl"
      
      $y
      [1] "n"
      
      $colour
      [1] "cyl"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      expression(list("Pairwise test:" ~ bold("Dunn test"), "Comparisons shown:" ~ 
          bold("only significant")))
      
      $label
      [1] "outlier.label"
      
      $alt
      [1] ""
      

---

    Code
      pb1$data[[6]]$label
    Output
      widehat(mu)[mean]=='0.98'
      widehat(mu)[mean]=='1.39'

# checking if plot.type argument works

    Code
      list(pb1$data, list(pb2$data[[1]], pb2$data[[2]], pb2$data[[4]], pb2$data[[5]]))
    Output
      [[1]]
      [[1]][[1]]
            colour         x    y PANEL group      xmin     xmax ymax shape size fill
      1  #D95F02FF 1.9882651  4.2     1     2 1.7845138 2.211148  4.2    19    3   NA
      2  #D95F02FF 2.0004730 11.5     1     2 1.6380827 2.200235 11.5    19    3   NA
      3  #D95F02FF 2.0581406  7.3     1     2 1.6765023 2.257072  7.3    19    3   NA
      4  #D95F02FF 1.9333987  5.8     1     2 1.7220795 2.337536  5.8    19    3   NA
      5  #D95F02FF 2.0358134  6.4     1     2 1.6049761 2.304028  6.4    19    3   NA
      6  #D95F02FF 2.0013700 10.0     1     2 1.7102413 2.241469 10.0    19    3   NA
      7  #D95F02FF 1.9780740 11.2     1     2 1.6724921 2.304770 11.2    19    3   NA
      8  #D95F02FF 2.0274448 11.2     1     2 1.7476213 2.353298 11.2    19    3   NA
      9  #D95F02FF 2.0285544  5.2     1     2 1.7675651 2.336127  5.2    19    3   NA
      10 #D95F02FF 1.9167764  7.0     1     2 1.6104254 2.362067  7.0    19    3   NA
      11 #D95F02FF 2.0558664 16.5     1     2 1.7600473 2.262006 16.5    19    3   NA
      12 #D95F02FF 1.9101503 16.5     1     2 1.7796882 2.286481 16.5    19    3   NA
      13 #D95F02FF 2.0722450 15.2     1     2 1.7756310 2.320879 15.2    19    3   NA
      14 #D95F02FF 2.0248003 17.3     1     2 1.6334292 2.368371 17.3    19    3   NA
      15 #D95F02FF 2.0965936 22.5     1     2 1.7156337 2.292739 22.5    19    3   NA
      16 #D95F02FF 1.9492710 17.3     1     2 1.7933106 2.245427 17.3    19    3   NA
      17 #D95F02FF 2.0772067 13.6     1     2 1.7268125 2.286106 13.6    19    3   NA
      18 #D95F02FF 2.0013116 14.5     1     2 1.6577052 2.208841 14.5    19    3   NA
      19 #D95F02FF 1.9847245 18.8     1     2 1.6365225 2.290436 18.8    19    3   NA
      20 #D95F02FF 2.0476592 15.5     1     2 1.7872508 2.300515 15.5    19    3   NA
      21 #D95F02FF 1.9921384 23.6     1     2 1.6683149 2.355669 23.6    19    3   NA
      22 #D95F02FF 2.0519402 18.5     1     2 1.7173039 2.314872 18.5    19    3   NA
      23 #D95F02FF 1.9358612 25.5     1     2 1.6663978 2.252023 25.5    19    3   NA
      24 #D95F02FF 2.0969609 26.4     1     2 1.6840770 2.351950 26.4    19    3   NA
      25 #D95F02FF 1.9844649 26.7     1     2 1.6733593 2.309716 26.7    19    3   NA
      26 #D95F02FF 2.0884042 21.5     1     2 1.6603246 2.233631 21.5    19    3   NA
      27 #D95F02FF 1.9882463 23.3     1     2 1.6496098 2.287720 23.3    19    3   NA
      28 #D95F02FF 1.9878564 29.5     1     2 1.7721896 2.319106 29.5    19    3   NA
      29 #1B9E77FF 0.9648859 15.2     1     1 0.7066215 1.235765 15.2    19    3   NA
      30 #1B9E77FF 0.9857474 21.5     1     1 0.7103754 1.382811 21.5    19    3   NA
      31 #1B9E77FF 1.0009592 17.6     1     1 0.7489079 1.310250 17.6    19    3   NA
      32 #1B9E77FF 0.9393401  9.7     1     1 0.7999891 1.230516  9.7    19    3   NA
      33 #1B9E77FF 1.0051387 14.5     1     1 0.7068992 1.383165 14.5    19    3   NA
      34 #1B9E77FF 0.9217320 10.0     1     1 0.6224561 1.240448 10.0    19    3   NA
      35 #1B9E77FF 0.9317606  8.2     1     1 0.7186037 1.247457  8.2    19    3   NA
      36 #1B9E77FF 0.9588684  9.4     1     1 0.7543105 1.216737  9.4    19    3   NA
      37 #1B9E77FF 0.9926730 16.5     1     1 0.7982347 1.252832 16.5    19    3   NA
      38 #1B9E77FF 0.9193843  9.7     1     1 0.6229418 1.208658  9.7    19    3   NA
      39 #1B9E77FF 0.9862866 19.7     1     1 0.6307405 1.377247 19.7    19    3   NA
      40 #1B9E77FF 0.9826531 23.3     1     1 0.6097494 1.292137 23.3    19    3   NA
      41 #1B9E77FF 0.9808246 23.6     1     1 0.7309071 1.375999 23.6    19    3   NA
      42 #1B9E77FF 0.9594821 26.4     1     1 0.6011048 1.276874 26.4    19    3   NA
      43 #1B9E77FF 1.0560578 20.0     1     1 0.6751588 1.352212 20.0    19    3   NA
      44 #1B9E77FF 0.9713816 25.2     1     1 0.6292349 1.311327 25.2    19    3   NA
      45 #1B9E77FF 1.0511830 25.8     1     1 0.7342666 1.365015 25.8    19    3   NA
      46 #1B9E77FF 1.0371307 21.2     1     1 0.6945312 1.378850 21.2    19    3   NA
      47 #1B9E77FF 0.9296905 14.5     1     1 0.7082889 1.336210 14.5    19    3   NA
      48 #1B9E77FF 0.9557634 27.3     1     1 0.6874114 1.353995 27.3    19    3   NA
      49 #1B9E77FF 0.9885489 25.5     1     1 0.7404557 1.302889 25.5    19    3   NA
      50 #1B9E77FF 1.0246714 26.4     1     1 0.7702690 1.386357 26.4    19    3   NA
      51 #1B9E77FF 1.0546359 22.4     1     1 0.6289893 1.261085 22.4    19    3   NA
      52 #1B9E77FF 1.0772551 24.5     1     1 0.7668815 1.231970 24.5    19    3   NA
      53 #1B9E77FF 0.9049347 24.8     1     1 0.6766541 1.210607 24.8    19    3   NA
      54 #1B9E77FF 1.0823804 30.9     1     1 0.7120424 1.247007 30.9    19    3   NA
      55 #1B9E77FF 1.0386884 26.4     1     1 0.6422930 1.274127 26.4    19    3   NA
      56 #1B9E77FF 0.9896083 27.3     1     1 0.6694498 1.276823 27.3    19    3   NA
      57 #1B9E77FF 0.9964303 29.4     1     1 0.7578911 1.296724 29.4    19    3   NA
      58 #1B9E77FF 1.0823791 23.0     1     1 0.7282850 1.358557 23.0    19    3   NA
         alpha stroke
      1    0.4      0
      2    0.4      0
      3    0.4      0
      4    0.4      0
      5    0.4      0
      6    0.4      0
      7    0.4      0
      8    0.4      0
      9    0.4      0
      10   0.4      0
      11   0.4      0
      12   0.4      0
      13   0.4      0
      14   0.4      0
      15   0.4      0
      16   0.4      0
      17   0.4      0
      18   0.4      0
      19   0.4      0
      20   0.4      0
      21   0.4      0
      22   0.4      0
      23   0.4      0
      24   0.4      0
      25   0.4      0
      26   0.4      0
      27   0.4      0
      28   0.4      0
      29   0.4      0
      30   0.4      0
      31   0.4      0
      32   0.4      0
      33   0.4      0
      34   0.4      0
      35   0.4      0
      36   0.4      0
      37   0.4      0
      38   0.4      0
      39   0.4      0
      40   0.4      0
      41   0.4      0
      42   0.4      0
      43   0.4      0
      44   0.4      0
      45   0.4      0
      46   0.4      0
      47   0.4      0
      48   0.4      0
      49   0.4      0
      50   0.4      0
      51   0.4      0
      52   0.4      0
      53   0.4      0
      54   0.4      0
      55   0.4      0
      56   0.4      0
      57   0.4      0
      58   0.4      0
      
      [[1]][[2]]
        ymin  lower middle  upper ymax   outliers notchupper notchlower x flipped_aes
      1  8.2 15.525   22.7 25.725 30.9              25.64237   19.75763 1       FALSE
      2  4.2 11.200   16.5 23.100 29.5 33.9, 32.5   19.93276   13.06724 2       FALSE
        PANEL group ymin_final ymax_final xmin xmax xid newx new_width weight colour
      1     1     1        8.2       30.9 0.85 1.15   1    1       0.3      1 grey20
      2     1     2        4.2       33.9 1.85 2.15   2    2       0.3      1 grey20
         fill size alpha shape linetype
      1 white  0.5   0.2    19    solid
      2 white  0.5   0.2    19    solid
      
      [[1]][[3]]
        x    y label PANEL group colour  fill size angle alpha family fontface
      1 2 33.9  33.9     1     1  black white    3     0    NA               1
      2 2 32.5  32.5     1     1  black white    3     0    NA               1
        lineheight hjust vjust point.size segment.linetype segment.size
      1        1.2   0.5   0.5          1                1          0.5
      2        1.2   0.5   0.5          1                1          0.5
        segment.curvature segment.angle segment.ncp segment.shape segment.square
      1                 0            90           1           0.5           TRUE
      2                 0            90           1           0.5           TRUE
        segment.squareShape segment.inflect segment.debug
      1                   1           FALSE         FALSE
      2                   1           FALSE         FALSE
      
      [[1]][[4]]
        x        y PANEL group shape    colour size fill alpha stroke
      1 1 24.90880     1     1    19 darkgreen    5   NA    NA    0.5
      2 2 15.43548     1     2    19 darkgreen    5   NA    NA    0.5
      
      [[1]][[5]]
        x        y                     label PANEL group nudge_x  nudge_y colour
      1 1 24.90880 widehat(mu)[MAP]=='24.91'     1     1     1.4 24.90880   blue
      2 2 15.43548 widehat(mu)[MAP]=='15.44'     1     2     2.4 15.43548   blue
         fill size angle alpha family fontface lineheight hjust vjust point.size
      1 white 3.88     0    NA               1        1.2   0.5   0.5          1
      2 white 3.88     0    NA               1        1.2   0.5   0.5          1
        segment.linetype segment.size segment.curvature segment.angle segment.ncp
      1                4          0.5                 0            90           1
      2                4          0.5                 0            90           1
        segment.shape segment.square segment.squareShape segment.inflect
      1           0.5           TRUE                   1           FALSE
      2           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      
      
      [[2]]
      [[2]][[1]]
            colour         x    y PANEL group      xmin     xmax ymax shape size fill
      1  #C93312FF 1.9882651  4.2     1     2 1.7845138 2.211148  4.2    19    3   NA
      2  #C93312FF 2.0004730 11.5     1     2 1.6380827 2.200235 11.5    19    3   NA
      3  #C93312FF 2.0581406  7.3     1     2 1.6765023 2.257072  7.3    19    3   NA
      4  #C93312FF 1.9333987  5.8     1     2 1.7220795 2.337536  5.8    19    3   NA
      5  #C93312FF 2.0358134  6.4     1     2 1.6049761 2.304028  6.4    19    3   NA
      6  #C93312FF 2.0013700 10.0     1     2 1.7102413 2.241469 10.0    19    3   NA
      7  #C93312FF 1.9780740 11.2     1     2 1.6724921 2.304770 11.2    19    3   NA
      8  #C93312FF 2.0274448 11.2     1     2 1.7476213 2.353298 11.2    19    3   NA
      9  #C93312FF 2.0285544  5.2     1     2 1.7675651 2.336127  5.2    19    3   NA
      10 #C93312FF 1.9167764  7.0     1     2 1.6104254 2.362067  7.0    19    3   NA
      11 #C93312FF 2.0558664 16.5     1     2 1.7600473 2.262006 16.5    19    3   NA
      12 #C93312FF 1.9101503 16.5     1     2 1.7796882 2.286481 16.5    19    3   NA
      13 #C93312FF 2.0722450 15.2     1     2 1.7756310 2.320879 15.2    19    3   NA
      14 #C93312FF 2.0248003 17.3     1     2 1.6334292 2.368371 17.3    19    3   NA
      15 #C93312FF 2.0965936 22.5     1     2 1.7156337 2.292739 22.5    19    3   NA
      16 #C93312FF 1.9492710 17.3     1     2 1.7933106 2.245427 17.3    19    3   NA
      17 #C93312FF 2.0772067 13.6     1     2 1.7268125 2.286106 13.6    19    3   NA
      18 #C93312FF 2.0013116 14.5     1     2 1.6577052 2.208841 14.5    19    3   NA
      19 #C93312FF 1.9847245 18.8     1     2 1.6365225 2.290436 18.8    19    3   NA
      20 #C93312FF 2.0476592 15.5     1     2 1.7872508 2.300515 15.5    19    3   NA
      21 #C93312FF 1.9921384 23.6     1     2 1.6683149 2.355669 23.6    19    3   NA
      22 #C93312FF 2.0519402 18.5     1     2 1.7173039 2.314872 18.5    19    3   NA
      23 #C93312FF 1.9358612 25.5     1     2 1.6663978 2.252023 25.5    19    3   NA
      24 #C93312FF 2.0969609 26.4     1     2 1.6840770 2.351950 26.4    19    3   NA
      25 #C93312FF 1.9844649 26.7     1     2 1.6733593 2.309716 26.7    19    3   NA
      26 #C93312FF 2.0884042 21.5     1     2 1.6603246 2.233631 21.5    19    3   NA
      27 #C93312FF 1.9882463 23.3     1     2 1.6496098 2.287720 23.3    19    3   NA
      28 #C93312FF 1.9878564 29.5     1     2 1.7721896 2.319106 29.5    19    3   NA
      29 #899DA4FF 0.9648859 15.2     1     1 0.7066215 1.235765 15.2    19    3   NA
      30 #899DA4FF 0.9857474 21.5     1     1 0.7103754 1.382811 21.5    19    3   NA
      31 #899DA4FF 1.0009592 17.6     1     1 0.7489079 1.310250 17.6    19    3   NA
      32 #899DA4FF 0.9393401  9.7     1     1 0.7999891 1.230516  9.7    19    3   NA
      33 #899DA4FF 1.0051387 14.5     1     1 0.7068992 1.383165 14.5    19    3   NA
      34 #899DA4FF 0.9217320 10.0     1     1 0.6224561 1.240448 10.0    19    3   NA
      35 #899DA4FF 0.9317606  8.2     1     1 0.7186037 1.247457  8.2    19    3   NA
      36 #899DA4FF 0.9588684  9.4     1     1 0.7543105 1.216737  9.4    19    3   NA
      37 #899DA4FF 0.9926730 16.5     1     1 0.7982347 1.252832 16.5    19    3   NA
      38 #899DA4FF 0.9193843  9.7     1     1 0.6229418 1.208658  9.7    19    3   NA
      39 #899DA4FF 0.9862866 19.7     1     1 0.6307405 1.377247 19.7    19    3   NA
      40 #899DA4FF 0.9826531 23.3     1     1 0.6097494 1.292137 23.3    19    3   NA
      41 #899DA4FF 0.9808246 23.6     1     1 0.7309071 1.375999 23.6    19    3   NA
      42 #899DA4FF 0.9594821 26.4     1     1 0.6011048 1.276874 26.4    19    3   NA
      43 #899DA4FF 1.0560578 20.0     1     1 0.6751588 1.352212 20.0    19    3   NA
      44 #899DA4FF 0.9713816 25.2     1     1 0.6292349 1.311327 25.2    19    3   NA
      45 #899DA4FF 1.0511830 25.8     1     1 0.7342666 1.365015 25.8    19    3   NA
      46 #899DA4FF 1.0371307 21.2     1     1 0.6945312 1.378850 21.2    19    3   NA
      47 #899DA4FF 0.9296905 14.5     1     1 0.7082889 1.336210 14.5    19    3   NA
      48 #899DA4FF 0.9557634 27.3     1     1 0.6874114 1.353995 27.3    19    3   NA
      49 #899DA4FF 0.9885489 25.5     1     1 0.7404557 1.302889 25.5    19    3   NA
      50 #899DA4FF 1.0246714 26.4     1     1 0.7702690 1.386357 26.4    19    3   NA
      51 #899DA4FF 1.0546359 22.4     1     1 0.6289893 1.261085 22.4    19    3   NA
      52 #899DA4FF 1.0772551 24.5     1     1 0.7668815 1.231970 24.5    19    3   NA
      53 #899DA4FF 0.9049347 24.8     1     1 0.6766541 1.210607 24.8    19    3   NA
      54 #899DA4FF 1.0823804 30.9     1     1 0.7120424 1.247007 30.9    19    3   NA
      55 #899DA4FF 1.0386884 26.4     1     1 0.6422930 1.274127 26.4    19    3   NA
      56 #899DA4FF 0.9896083 27.3     1     1 0.6694498 1.276823 27.3    19    3   NA
      57 #899DA4FF 0.9964303 29.4     1     1 0.7578911 1.296724 29.4    19    3   NA
      58 #899DA4FF 1.0823791 23.0     1     1 0.7282850 1.358557 23.0    19    3   NA
         alpha stroke
      1    0.4      0
      2    0.4      0
      3    0.4      0
      4    0.4      0
      5    0.4      0
      6    0.4      0
      7    0.4      0
      8    0.4      0
      9    0.4      0
      10   0.4      0
      11   0.4      0
      12   0.4      0
      13   0.4      0
      14   0.4      0
      15   0.4      0
      16   0.4      0
      17   0.4      0
      18   0.4      0
      19   0.4      0
      20   0.4      0
      21   0.4      0
      22   0.4      0
      23   0.4      0
      24   0.4      0
      25   0.4      0
      26   0.4      0
      27   0.4      0
      28   0.4      0
      29   0.4      0
      30   0.4      0
      31   0.4      0
      32   0.4      0
      33   0.4      0
      34   0.4      0
      35   0.4      0
      36   0.4      0
      37   0.4      0
      38   0.4      0
      39   0.4      0
      40   0.4      0
      41   0.4      0
      42   0.4      0
      43   0.4      0
      44   0.4      0
      45   0.4      0
      46   0.4      0
      47   0.4      0
      48   0.4      0
      49   0.4      0
      50   0.4      0
      51   0.4      0
      52   0.4      0
      53   0.4      0
      54   0.4      0
      55   0.4      0
      56   0.4      0
      57   0.4      0
      58   0.4      0
      
      [[2]][[2]]
        x    y PANEL group shape colour size fill alpha stroke
      1 2 33.9     1     1    19  black    3   NA   0.7      0
      2 2 32.5     1     1    19  black    3   NA   0.7      0
      
      [[2]][[3]]
        x    y label PANEL group colour  fill size angle alpha family fontface
      1 2 33.9  33.9     1     1  black white    3     0    NA               1
      2 2 32.5  32.5     1     1  black white    3     0    NA               1
        lineheight hjust vjust point.size segment.linetype segment.size
      1        1.2   0.5   0.5          1                1          0.5
      2        1.2   0.5   0.5          1                1          0.5
        segment.curvature segment.angle segment.ncp segment.shape segment.square
      1                 0            90           1           0.5           TRUE
      2                 0            90           1           0.5           TRUE
        segment.squareShape segment.inflect segment.debug
      1                   1           FALSE         FALSE
      2                   1           FALSE         FALSE
      
      [[2]][[4]]
        x        y PANEL group shape  colour size fill alpha stroke
      1 1 21.70556     1     1    19 darkred    5   NA    NA    0.5
      2 2 16.58333     1     2    19 darkred    5   NA    NA    0.5
      
      

---

    Code
      within(pb1$plot$labels, rm(subtitle, caption))
    Output
      $x
      [1] "supp"
      
      $y
      [1] "len"
      
      $colour
      [1] "supp"
      
      $title
      NULL
      
      $label
      [1] "outlier.label"
      
      $alt
      [1] ""
      

---

    Code
      within(pb2$plot$labels, rm(subtitle))
    Output
      $x
      [1] "supp"
      
      $y
      [1] "len"
      
      $colour
      [1] "supp"
      
      $title
      NULL
      
      $caption
      NULL
      
      $label
      [1] "outlier.label"
      
      $alt
      [1] ""
      

# checking ggscatterstats - without NAs - pearson's r

    Code
      list(pb$data[[1]], head(pb$data[[2]]), pb$data[[3]])
    Output
      [[1]]
           x         y PANEL group shape colour size fill alpha stroke
      1 14.9 0.1333333     1    -1    19  black    3   NA   0.4      0
      2  9.1 0.1500000     1    -1    19  black    3   NA   0.4      0
      3 17.4 0.3833333     1    -1    19  black    3   NA   0.4      0
      4 18.0 0.3333333     1    -1    19  black    3   NA   0.4      0
      5 19.7 0.1166667     1    -1    19  black    3   NA   0.4      0
      6 10.1 0.2833333     1    -1    19  black    3   NA   0.4      0
      7 13.0 0.1833333     1    -1    19  black    3   NA   0.4      0
      8  8.4 0.1666667     1    -1    19  black    3   NA   0.4      0
      9 13.8 0.2166667     1    -1    19  black    3   NA   0.4      0
      
      [[2]]
               x         y       ymin      ymax         se flipped_aes PANEL group
      1 8.400000 0.1868825 0.05397985 0.3197852 0.05620456       FALSE     1    -1
      2 8.543038 0.1877171 0.05714189 0.3182922 0.05522027       FALSE     1    -1
      3 8.686076 0.1885516 0.06028225 0.3168210 0.05424514       FALSE     1    -1
      4 8.829114 0.1893862 0.06339976 0.3153726 0.05327968       FALSE     1    -1
      5 8.972152 0.1902207 0.06649314 0.3139483 0.05232442       FALSE     1    -1
      6 9.115190 0.1910553 0.06956104 0.3125495 0.05137994       FALSE     1    -1
        colour   fill size linetype weight alpha
      1   blue grey60  1.5        1      1   0.4
      2   blue grey60  1.5        1      1   0.4
      3   blue grey60  1.5        1      1   0.4
      4   blue grey60  1.5        1      1   0.4
      5   blue grey60  1.5        1      1   0.4
      6   blue grey60  1.5        1      1   0.4
      
      [[3]]
           x         y                  label PANEL group colour  fill size angle
      1 17.4 0.3833333   Long-nosed armadillo     1    -1  black white    3     0
      2 18.0 0.3333333 North American Opossum     1    -1  black white    3     0
        alpha family fontface lineheight hjust vjust point.size segment.linetype
      1    NA               1        1.2   0.5   0.5          1                1
      2    NA               1        1.2   0.5   0.5          1                1
        segment.size segment.curvature segment.angle segment.ncp segment.shape
      1          0.5                 0            90           1           0.5
      2          0.5                 0            90           1           0.5
        segment.square segment.squareShape segment.inflect segment.debug
      1           TRUE                   1           FALSE         FALSE
      2           TRUE                   1           FALSE         FALSE
      

---

    Code
      within(pb$plot$labels, rm(subtitle, caption))
    Output
      $x
      [1] "sleep (total)"
      
      $y
      [1] "sleep cycle"
      
      $title
      [1] "Mammalian sleep"
      
      $label
      [1] "name"
      
      $alt
      [1] ""
      

# checking ggscatterstats - without NAs - spearman's rho

    Code
      pb$data[[1]]
    Output
           x         y PANEL group shape colour size fill alpha stroke
      1 14.9 0.1333333     1    -1    19  black    3   NA   0.4      0
      2  9.1 0.1500000     1    -1    19  black    3   NA   0.4      0
      3 17.4 0.3833333     1    -1    19  black    3   NA   0.4      0
      4 18.0 0.3333333     1    -1    19  black    3   NA   0.4      0
      5 19.7 0.1166667     1    -1    19  black    3   NA   0.4      0
      6 10.1 0.2833333     1    -1    19  black    3   NA   0.4      0
      7 13.0 0.1833333     1    -1    19  black    3   NA   0.4      0
      8  8.4 0.1666667     1    -1    19  black    3   NA   0.4      0
      9 13.8 0.2166667     1    -1    19  black    3   NA   0.4      0

---

    Code
      within(pb$plot$labels, rm(subtitle))
    Output
      $x
      [1] "sleep_total"
      
      $y
      [1] "sleep_cycle"
      
      $title
      NULL
      
      $caption
      NULL
      
      $alt
      [1] ""
      

# checking ggscatterstats - without NAs - winsorized Pearson

    Code
      pb$data[[1]]
    Output
           x         y PANEL group shape colour size fill alpha stroke
      1 14.9 0.1333333     1    -1    19    red    5   NA    NA      0
      2  9.1 0.1500000     1    -1    19    red    5   NA    NA      0
      3 17.4 0.3833333     1    -1    19    red    5   NA    NA      0
      4 18.0 0.3333333     1    -1    19    red    5   NA    NA      0
      5 19.7 0.1166667     1    -1    19    red    5   NA    NA      0
      6 10.1 0.2833333     1    -1    19    red    5   NA    NA      0
      7 13.0 0.1833333     1    -1    19    red    5   NA    NA      0
      8  8.4 0.1666667     1    -1    19    red    5   NA    NA      0
      9 13.8 0.2166667     1    -1    19    red    5   NA    NA      0

---

    Code
      within(pb$plot$labels, rm(subtitle))
    Output
      $x
      [1] "sleep_total"
      
      $y
      [1] "sleep_cycle"
      
      $title
      NULL
      
      $caption
      NULL
      
      $alt
      [1] ""
      

# bayes factor plus class of object

    Code
      purrr::map(pb$data, names)
    Output
      [[1]]
       [1] "y"      "x"      "PANEL"  "group"  "shape"  "colour" "size"   "fill"  
       [9] "alpha"  "stroke"
      
      [[2]]
       [1] "x"           "y"           "ymin"        "ymax"        "se"         
       [6] "flipped_aes" "PANEL"       "group"       "colour"      "fill"       
      [11] "size"        "linetype"    "weight"      "alpha"      
      
      [[3]]
       [1] "y"           "count"       "x"           "xmin"        "xmax"       
       [6] "density"     "ncount"      "ndensity"    "flipped_aes" "PANEL"      
      [11] "group"       "ymin"        "ymax"        "colour"      "xcolour"    
      [16] "fill"        "xfill"       "size"        "linetype"    "alpha"      
      
      [[4]]
       [1] "x"           "count"       "y"           "ymin"        "ymax"       
       [6] "density"     "ncount"      "ndensity"    "flipped_aes" "PANEL"      
      [11] "group"       "xmin"        "xmax"        "colour"      "ycolour"    
      [16] "fill"        "yfill"       "size"        "linetype"    "alpha"      
      

---

    Code
      purrr::map(pb$data, dim)
    Output
      [[1]]
      [1]  9 10
      
      [[2]]
      [1] 80 14
      
      [[3]]
      [1] 30 20
      
      [[4]]
      [1] 30 20
      

# aesthetic modifications work

    Code
      list(pb$data[[1]], head(pb$data[[2]]), pb$data[[3]])
    Output
      [[1]]
                y    x PANEL group shape colour size fill alpha stroke
      1 0.1333333 14.9     2    -1    19  black    3   NA   0.4      0
      2 0.1500000  9.1     2    -1    19  black    3   NA   0.4      0
      3 0.3833333 17.4     2    -1    19  black    3   NA   0.4      0
      4 0.3333333 18.0     2    -1    19  black    3   NA   0.4      0
      5 0.1166667 19.7     2    -1    19  black    3   NA   0.4      0
      6 0.2833333 10.1     2    -1    19  black    3   NA   0.4      0
      7 0.1833333 13.0     2    -1    19  black    3   NA   0.4      0
      8 0.1666667  8.4     2    -1    19  black    3   NA   0.4      0
      9 0.2166667 13.8     2    -1    19  black    3   NA   0.4      0
      
      [[2]]
               x         y       ymin      ymax         se flipped_aes PANEL group
      1 8.400000 0.1868825 0.05397985 0.3197852 0.05620456       FALSE     2    -1
      2 8.543038 0.1877171 0.05714189 0.3182922 0.05522027       FALSE     2    -1
      3 8.686076 0.1885516 0.06028225 0.3168210 0.05424514       FALSE     2    -1
      4 8.829114 0.1893862 0.06339976 0.3153726 0.05327968       FALSE     2    -1
      5 8.972152 0.1902207 0.06649314 0.3139483 0.05232442       FALSE     2    -1
      6 9.115190 0.1910553 0.06956104 0.3125495 0.05137994       FALSE     2    -1
        colour   fill size linetype weight alpha
      1   blue grey60  1.5        1      1   0.4
      2   blue grey60  1.5        1      1   0.4
      3   blue grey60  1.5        1      1   0.4
      4   blue grey60  1.5        1      1   0.4
      5   blue grey60  1.5        1      1   0.4
      6   blue grey60  1.5        1      1   0.4
      
      [[3]]
                y    x           label PANEL group colour  fill size angle alpha
      1 0.3833333 17.4       Cingulata     2    -1   blue white    4     0   0.5
      2 0.3333333 18.0 Didelphimorphia     2    -1   blue white    4     0   0.5
      3 0.1166667 19.7      Chiroptera     2    -1   blue white    4     0   0.5
        family fontface lineheight hjust vjust point.size segment.linetype
      1               1        1.2   0.5   0.5          1                1
      2               1        1.2   0.5   0.5          1                1
      3               1        1.2   0.5   0.5          1                1
        segment.size segment.curvature segment.angle segment.ncp segment.shape
      1          0.5                 0            90           1           0.5
      2          0.5                 0            90           1           0.5
      3          0.5                 0            90           1           0.5
        segment.square segment.squareShape segment.inflect segment.debug
      1           TRUE                   1           FALSE         FALSE
      2           TRUE                   1           FALSE         FALSE
      3           TRUE                   1           FALSE         FALSE
      

---

    Code
      pb$plot$labels
    Output
      $x
      [1] "sleep_total"
      
      $y
      [1] "sleep_cycle"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      NULL
      
      $label
      [1] "order"
      
      $weight
      [1] "weight"
      attr(,"fallback")
      [1] TRUE
      
      $alt
      [1] ""
      

# ggdotplotstats works as expected

    Code
      within(pb$plot$labels, rm(subtitle, caption))
    Output
      $x
      paste("Speed of light (", italic("c"), ")")
      
      $y
      [1] "Experimental run"
      
      $title
      [1] "Michelson-Morley experiment"
      
      $xintercept
      [1] "xintercept"
      
      $alt
      [1] ""
      

---

    Code
      within(pb$layout$panel_params[[1]], rm(x, y, x.sec, y.sec))
    Output
      $x.arrange
      [1] "secondary" "primary"  
      
      $x.range
      [1] 816.075 913.425
      
      $y.arrange
      [1] "primary"   "secondary"
      
      $y.range
      [1] 0.8 5.2
      

---

    Code
      pb$data
    Output
      [[1]]
        y     x PANEL group shape colour size fill alpha stroke
      1 1 820.5     1    -1    16  black    3   NA    NA    0.5
      2 2 831.5     1    -1    16  black    3   NA    NA    0.5
      3 3 845.0     1    -1    16  black    3   NA    NA    0.5
      4 4 856.0     1    -1    16  black    3   NA    NA    0.5
      5 5 909.0     1    -1    16  black    3   NA    NA    0.5
      
      [[2]]
        xintercept PANEL group colour size linetype alpha
      1      852.4     1    -1   blue    1   dashed    NA
      

# messing with factors

    Code
      pb1$data
    Output
      [[1]]
        y          x PANEL group shape colour size fill alpha stroke
      1 1 0.02155000     1    -1    16  black    3   NA    NA    0.5
      2 2 0.07925556     1    -1    16  black    3   NA    NA    0.5
      3 3 0.62159750     1    -1    16  black    3   NA    NA    0.5
      
      [[2]]
        xintercept PANEL group colour size linetype alpha
      1   0.240801     1    -1   blue    1   dashed    NA
      

# checking if combining plots works

    Code
      pb$plot$labels
    Output
      $title
      [1] "versicolor"
      
      $x
      [1] "Sepal.Length"
      
      $y
      [1] "Sepal.Width"
      
      $alt
      [1] ""
      

---

    Code
      pb$plot$patches$annotation
    Output
      $title
      [1] "Dataset: Iris Flower dataset"
      
      $subtitle
      [1] "Edgar Anderson collected this data"
      
      $caption
      [1] "Note: Only two species of flower are displayed"
      
      $tag_levels
      [1] "a"
      
      $tag_prefix
      NULL
      
      $tag_suffix
      NULL
      
      $tag_sep
      NULL
      
      $theme
       Named list()
       - attr(*, "class")= chr [1:2] "theme" "gg"
       - attr(*, "complete")= logi FALSE
       - attr(*, "validate")= logi TRUE
      
      attr(,"class")
      [1] "plot_annotation"

# grouped_ggcorrmat stats work

    Code
      grouped_ggcorrmat(data = dplyr::select(ggplot2::msleep, dplyr::matches(
        "sleep|awake|vore")), grouping.var = vore, type = "r", output = "data", tr = 0.2)
    Output
      # A tibble: 24 x 12
         vore  parameter1  parameter2  estimate conf.level conf.low conf.high
         <chr> <chr>       <chr>          <dbl>      <dbl>    <dbl>     <dbl>
       1 carni sleep_total sleep_rem      0.948       0.95    0.790     0.988
       2 carni sleep_total sleep_cycle    0.632       0.95   -0.565     0.972
       3 carni sleep_total awake         -1.00        0.95   -1.00     -1.00 
       4 carni sleep_rem   sleep_cycle    0.333       0.95   -0.778     0.939
       5 carni sleep_rem   awake         -0.948       0.95   -0.988    -0.790
       6 carni sleep_cycle awake         -0.632       0.95   -0.972     0.565
       7 herbi sleep_total sleep_rem      0.900       0.95    0.780     0.956
       8 herbi sleep_total sleep_cycle   -0.677       0.95   -0.901    -0.169
       9 herbi sleep_total awake         -1           0.95   -1        -1    
      10 herbi sleep_rem   sleep_cycle   -0.343       0.95   -0.766     0.287
         statistic df.error  p.value method                         n.obs
             <dbl>    <int>    <dbl> <chr>                          <int>
       1     8.43         8 1.50e- 4 Winsorized Pearson correlation    10
       2     1.41         3 7.57e- 1 Winsorized Pearson correlation     5
       3  -962.          17 6.36e-41 Winsorized Pearson correlation    19
       4     0.612        3 7.57e- 1 Winsorized Pearson correlation     5
       5    -8.43         8 1.50e- 4 Winsorized Pearson correlation    10
       6    -1.41         3 7.57e- 1 Winsorized Pearson correlation     5
       7     9.69        22 1.07e- 8 Winsorized Pearson correlation    24
       8    -2.91        10 4.68e- 2 Winsorized Pearson correlation    12
       9  -Inf           30 0        Winsorized Pearson correlation    32
      10    -1.16        10 2.74e- 1 Winsorized Pearson correlation    12
      # ... with 14 more rows

# check mcp displays - parametric - significant

    Code
      list(pb$data[[6]], pb$data[[7]])
    Output
      [[1]]
        x          y                      label PANEL group nudge_x    nudge_y colour
      1 1 0.07925556 widehat(mu)[mean]=='0.079'     1     1     1.4 0.07925556  black
      2 2 0.62159750 widehat(mu)[mean]=='0.622'     1     2     2.4 0.62159750  black
      3 3 0.02155000 widehat(mu)[mean]=='0.022'     1     3     3.4 0.02155000  black
      4 4 0.14573118 widehat(mu)[mean]=='0.146'     1     4     4.4 0.14573118  black
         fill size angle alpha family fontface lineheight hjust vjust point.size
      1 white    3     0    NA               1        1.2   0.5   0.5          1
      2 white    3     0    NA               1        1.2   0.5   0.5          1
      3 white    3     0    NA               1        1.2   0.5   0.5          1
      4 white    3     0    NA               1        1.2   0.5   0.5          1
        segment.linetype segment.size segment.curvature segment.angle segment.ncp
      1                4          0.5                 0            90           1
      2                4          0.5                 0            90           1
      3                4          0.5                 0            90           1
      4                4          0.5                 0            90           1
        segment.shape segment.square segment.squareShape segment.inflect
      1           0.5           TRUE                   1           FALSE
      2           0.5           TRUE                   1           FALSE
      3           0.5           TRUE                   1           FALSE
      4           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      3         FALSE
      4         FALSE
      
      [[2]]
         x xend        y     yend                           annotation
      1  1    1 6.083274 6.140393 list(~italic(p)[uncorrected]==0.437)
      2  1    2 6.140393 6.140393 list(~italic(p)[uncorrected]==0.437)
      3  2    2 6.140393 6.083274 list(~italic(p)[uncorrected]==0.437)
      4  1    1 6.425986 6.483105 list(~italic(p)[uncorrected]==0.452)
      5  1    3 6.483105 6.483105 list(~italic(p)[uncorrected]==0.452)
      6  3    3 6.483105 6.425986 list(~italic(p)[uncorrected]==0.452)
      7  1    1 6.768698 6.825816 list(~italic(p)[uncorrected]==0.865)
      8  1    4 6.825816 6.825816 list(~italic(p)[uncorrected]==0.865)
      9  4    4 6.825816 6.768698 list(~italic(p)[uncorrected]==0.865)
      10 2    2 7.111409 7.168528 list(~italic(p)[uncorrected]==0.348)
      11 2    3 7.168528 7.168528 list(~italic(p)[uncorrected]==0.348)
      12 3    3 7.168528 7.111409 list(~italic(p)[uncorrected]==0.348)
      13 2    2 7.454121 7.511239 list(~italic(p)[uncorrected]==0.560)
      14 2    4 7.511239 7.511239 list(~italic(p)[uncorrected]==0.560)
      15 4    4 7.511239 7.454121 list(~italic(p)[uncorrected]==0.560)
      16 3    3 7.796832 7.853951 list(~italic(p)[uncorrected]==0.433)
      17 3    4 7.853951 7.853951 list(~italic(p)[uncorrected]==0.433)
      18 4    4 7.853951 7.796832 list(~italic(p)[uncorrected]==0.433)
                   group flipped_aes PANEL shape colour textsize angle hjust vjust
      1    carni-herbi-1       FALSE     1    19  black        3     0   0.5     0
      2    carni-herbi-1       FALSE     1    19  black        3     0   0.5     0
      3    carni-herbi-1       FALSE     1    19  black        3     0   0.5     0
      4  carni-insecti-2       FALSE     1    19  black        3     0   0.5     0
      5  carni-insecti-2       FALSE     1    19  black        3     0   0.5     0
      6  carni-insecti-2       FALSE     1    19  black        3     0   0.5     0
      7     carni-omni-3       FALSE     1    19  black        3     0   0.5     0
      8     carni-omni-3       FALSE     1    19  black        3     0   0.5     0
      9     carni-omni-3       FALSE     1    19  black        3     0   0.5     0
      10 herbi-insecti-4       FALSE     1    19  black        3     0   0.5     0
      11 herbi-insecti-4       FALSE     1    19  black        3     0   0.5     0
      12 herbi-insecti-4       FALSE     1    19  black        3     0   0.5     0
      13    herbi-omni-5       FALSE     1    19  black        3     0   0.5     0
      14    herbi-omni-5       FALSE     1    19  black        3     0   0.5     0
      15    herbi-omni-5       FALSE     1    19  black        3     0   0.5     0
      16  insecti-omni-6       FALSE     1    19  black        3     0   0.5     0
      17  insecti-omni-6       FALSE     1    19  black        3     0   0.5     0
      18  insecti-omni-6       FALSE     1    19  black        3     0   0.5     0
         alpha family fontface lineheight linetype size
      1     NA               1        1.2        1  0.5
      2     NA               1        1.2        1  0.5
      3     NA               1        1.2        1  0.5
      4     NA               1        1.2        1  0.5
      5     NA               1        1.2        1  0.5
      6     NA               1        1.2        1  0.5
      7     NA               1        1.2        1  0.5
      8     NA               1        1.2        1  0.5
      9     NA               1        1.2        1  0.5
      10    NA               1        1.2        1  0.5
      11    NA               1        1.2        1  0.5
      12    NA               1        1.2        1  0.5
      13    NA               1        1.2        1  0.5
      14    NA               1        1.2        1  0.5
      15    NA               1        1.2        1  0.5
      16    NA               1        1.2        1  0.5
      17    NA               1        1.2        1  0.5
      18    NA               1        1.2        1  0.5
      

---

    Code
      pb$plot$labels
    Output
      $x
      [1] "vore"
      
      $y
      [1] "brainwt"
      
      $colour
      [1] "vore"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      expression(atop(displaystyle("mammalian sleep"), list("Pairwise test:" ~ 
          bold("Games-Howell test"), "Comparisons shown:" ~ bold("only non-significant"))))
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

# check mcp displays - non-significant

    Code
      list(pb1$data[[6]], pb1$data[[7]], pb2$data[[6]], pb2$data[[7]])
    Output
      [[1]]
        x         y                          label PANEL group nudge_x   nudge_y
      1 1  8440.335  widehat(mu)[mean]=='8440.335'     1     1     1.4  8440.335
      2 2 11148.255 widehat(mu)[mean]=='11148.255'     1     2     2.4 11148.255
      3 3  9243.369  widehat(mu)[mean]=='9243.369'     1     3     3.4  9243.369
        colour  fill size angle alpha family fontface lineheight hjust vjust
      1  black white    3     0    NA               1        1.2   0.5   0.5
      2  black white    3     0    NA               1        1.2   0.5   0.5
      3  black white    3     0    NA               1        1.2   0.5   0.5
        point.size segment.linetype segment.size segment.curvature segment.angle
      1          1                4          0.5                 0            90
      2          1                4          0.5                 0            90
      3          1                4          0.5                 0            90
        segment.ncp segment.shape segment.square segment.squareShape segment.inflect
      1           1           0.5           TRUE                   1           FALSE
      2           1           0.5           TRUE                   1           FALSE
      3           1           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      3         FALSE
      
      [[2]]
        x xend        y     yend                           annotation      group
      1 1    1 167852.2 169428.2 list(~italic(p)[uncorrected]==0.139) PG-PG-13-1
      2 1    2 169428.2 169428.2 list(~italic(p)[uncorrected]==0.139) PG-PG-13-1
      3 2    2 169428.2 167852.2 list(~italic(p)[uncorrected]==0.139) PG-PG-13-1
      4 1    1 179672.2 181248.2 list(~italic(p)[uncorrected]==0.825)     PG-R-2
      5 1    3 181248.2 181248.2 list(~italic(p)[uncorrected]==0.825)     PG-R-2
      6 3    3 181248.2 179672.2 list(~italic(p)[uncorrected]==0.825)     PG-R-2
      7 2    2 191492.2 193068.2 list(~italic(p)[uncorrected]==0.079)  PG-13-R-3
      8 2    3 193068.2 193068.2 list(~italic(p)[uncorrected]==0.079)  PG-13-R-3
      9 3    3 193068.2 191492.2 list(~italic(p)[uncorrected]==0.079)  PG-13-R-3
        flipped_aes PANEL shape colour textsize angle hjust vjust alpha family
      1       FALSE     1    19  black        3     0   0.5     0    NA       
      2       FALSE     1    19  black        3     0   0.5     0    NA       
      3       FALSE     1    19  black        3     0   0.5     0    NA       
      4       FALSE     1    19  black        3     0   0.5     0    NA       
      5       FALSE     1    19  black        3     0   0.5     0    NA       
      6       FALSE     1    19  black        3     0   0.5     0    NA       
      7       FALSE     1    19  black        3     0   0.5     0    NA       
      8       FALSE     1    19  black        3     0   0.5     0    NA       
      9       FALSE     1    19  black        3     0   0.5     0    NA       
        fontface lineheight linetype size
      1        1        1.2        1  0.5
      2        1        1.2        1  0.5
      3        1        1.2        1  0.5
      4        1        1.2        1  0.5
      5        1        1.2        1  0.5
      6        1        1.2        1  0.5
      7        1        1.2        1  0.5
      8        1        1.2        1  0.5
      9        1        1.2        1  0.5
      
      [[3]]
        x         y                           label PANEL group nudge_x   nudge_y
      1 1  8440.335  widehat(mu)[mean]=='8440.3350'     1     1     1.4  8440.335
      2 2 11148.255 widehat(mu)[mean]=='11148.2549'     1     2     2.4 11148.255
      3 3  9243.369  widehat(mu)[mean]=='9243.3687'     1     3     3.4  9243.369
        colour  fill size angle alpha family fontface lineheight hjust vjust
      1  black white    3     0    NA               1        1.2   0.5   0.5
      2  black white    3     0    NA               1        1.2   0.5   0.5
      3  black white    3     0    NA               1        1.2   0.5   0.5
        point.size segment.linetype segment.size segment.curvature segment.angle
      1          1                4          0.5                 0            90
      2          1                4          0.5                 0            90
      3          1                4          0.5                 0            90
        segment.ncp segment.shape segment.square segment.squareShape segment.inflect
      1           1           0.5           TRUE                   1           FALSE
      2           1           0.5           TRUE                   1           FALSE
      3           1           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      3         FALSE
      
      [[4]]
        x xend        y     yend                            annotation      group
      1 1    1 167852.2 169428.2 list(~italic(p)[uncorrected]==0.0467) PG-PG-13-1
      2 1    2 169428.2 169428.2 list(~italic(p)[uncorrected]==0.0467) PG-PG-13-1
      3 2    2 169428.2 167852.2 list(~italic(p)[uncorrected]==0.0467) PG-PG-13-1
      4 2    2 179672.2 181248.2 list(~italic(p)[uncorrected]==0.0354)  PG-13-R-2
      5 2    3 181248.2 181248.2 list(~italic(p)[uncorrected]==0.0354)  PG-13-R-2
      6 3    3 181248.2 179672.2 list(~italic(p)[uncorrected]==0.0354)  PG-13-R-2
        flipped_aes PANEL shape colour textsize angle hjust vjust alpha family
      1       FALSE     1    19  black        3     0   0.5     0    NA       
      2       FALSE     1    19  black        3     0   0.5     0    NA       
      3       FALSE     1    19  black        3     0   0.5     0    NA       
      4       FALSE     1    19  black        3     0   0.5     0    NA       
      5       FALSE     1    19  black        3     0   0.5     0    NA       
      6       FALSE     1    19  black        3     0   0.5     0    NA       
        fontface lineheight linetype size
      1        1        1.2        1  0.5
      2        1        1.2        1  0.5
      3        1        1.2        1  0.5
      4        1        1.2        1  0.5
      5        1        1.2        1  0.5
      6        1        1.2        1  0.5
      

---

    Code
      list(pb1$plot$labels, pb2$plot$labels)
    Output
      [[1]]
      [[1]]$x
      [1] "mpaa"
      
      [[1]]$y
      [1] "votes"
      
      [[1]]$colour
      [1] "mpaa"
      
      [[1]]$title
      NULL
      
      [[1]]$subtitle
      NULL
      
      [[1]]$caption
      expression(list("Pairwise test:" ~ bold("Games-Howell test"), 
          "Comparisons shown:" ~ bold("only non-significant")))
      
      [[1]]$label
      [1] "expression"
      
      [[1]]$alt
      [1] ""
      
      
      [[2]]
      [[2]]$x
      [1] "mpaa"
      
      [[2]]$y
      [1] "votes"
      
      [[2]]$colour
      [1] "mpaa"
      
      [[2]]$title
      NULL
      
      [[2]]$subtitle
      NULL
      
      [[2]]$caption
      expression(list("Pairwise test:" ~ bold("Student's t-test"), 
          "Comparisons shown:" ~ bold("only significant")))
      
      [[2]]$label
      [1] "expression"
      
      [[2]]$alt
      [1] ""
      
      

# check mixed comparison displays - nonparametric

    Code
      list(pb$data[[6]], pb$data[[7]])
    Output
      [[1]]
        x   y                        label PANEL group nudge_x nudge_y colour  fill
      1 1 5.5 widehat(mu)[median]=='5.500'     1     1     1.4     5.5  black white
      2 2 5.5 widehat(mu)[median]=='5.500'     1     2     2.4     5.5  black white
      3 3 5.9 widehat(mu)[median]=='5.900'     1     3     3.4     5.9  black white
        size angle alpha family fontface lineheight hjust vjust point.size
      1    3     0    NA               1        1.2   0.5   0.5          1
      2    3     0    NA               1        1.2   0.5   0.5          1
      3    3     0    NA               1        1.2   0.5   0.5          1
        segment.linetype segment.size segment.curvature segment.angle segment.ncp
      1                4          0.5                 0            90           1
      2                4          0.5                 0            90           1
      3                4          0.5                 0            90           1
        segment.shape segment.square segment.squareShape segment.inflect
      1           0.5           TRUE                   1           FALSE
      2           0.5           TRUE                   1           FALSE
      3           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      3         FALSE
      
      [[2]]
        x xend       y    yend                                 annotation
      1 1    1  9.5170  9.5900     list(~italic(p)[FDR-corrected]==0.812)
      2 1    2  9.5900  9.5900     list(~italic(p)[FDR-corrected]==0.812)
      3 2    2  9.5900  9.5170     list(~italic(p)[FDR-corrected]==0.812)
      4 1    1 10.0645 10.1375 list(~italic(p)[FDR-corrected]==4.179e-04)
      5 1    3 10.1375 10.1375 list(~italic(p)[FDR-corrected]==4.179e-04)
      6 3    3 10.1375 10.0645 list(~italic(p)[FDR-corrected]==4.179e-04)
      7 2    2 10.6120 10.6850 list(~italic(p)[FDR-corrected]==4.179e-04)
      8 2    3 10.6850 10.6850 list(~italic(p)[FDR-corrected]==4.179e-04)
      9 3    3 10.6850 10.6120 list(~italic(p)[FDR-corrected]==4.179e-04)
                  group flipped_aes PANEL shape colour textsize angle hjust vjust
      1 Action-Comedy-1       FALSE     1    19  black        3     0   0.5     0
      2 Action-Comedy-1       FALSE     1    19  black        3     0   0.5     0
      3 Action-Comedy-1       FALSE     1    19  black        3     0   0.5     0
      4 Action-RomCom-2       FALSE     1    19  black        3     0   0.5     0
      5 Action-RomCom-2       FALSE     1    19  black        3     0   0.5     0
      6 Action-RomCom-2       FALSE     1    19  black        3     0   0.5     0
      7 Comedy-RomCom-3       FALSE     1    19  black        3     0   0.5     0
      8 Comedy-RomCom-3       FALSE     1    19  black        3     0   0.5     0
      9 Comedy-RomCom-3       FALSE     1    19  black        3     0   0.5     0
        alpha family fontface lineheight linetype size
      1    NA               1        1.2        1  0.5
      2    NA               1        1.2        1  0.5
      3    NA               1        1.2        1  0.5
      4    NA               1        1.2        1  0.5
      5    NA               1        1.2        1  0.5
      6    NA               1        1.2        1  0.5
      7    NA               1        1.2        1  0.5
      8    NA               1        1.2        1  0.5
      9    NA               1        1.2        1  0.5
      

---

    Code
      pb$plot$labels
    Output
      $x
      [1] "genre"
      
      $y
      [1] "rating"
      
      $colour
      [1] "genre"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      expression(list("Pairwise test:" ~ bold("Dunn test"), "Comparisons shown:" ~ 
          bold("all")))
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

# check robust test display - FDR-corrected

    Code
      list(pb$data[[6]], pb$data[[7]])
    Output
      [[1]]
        x        y                          label PANEL group nudge_x  nudge_y colour
      1 1 14.07937 widehat(mu)[trimmed]=='14.079'     1     1     1.4 14.07937  black
      2 2 19.43750 widehat(mu)[trimmed]=='19.438'     1     2     2.4 19.43750  black
      3 3 14.13333 widehat(mu)[trimmed]=='14.133'     1     3     3.4 14.13333  black
         fill size angle alpha family fontface lineheight hjust vjust point.size
      1 white    3     0    NA               1        1.2   0.5   0.5          1
      2 white    3     0    NA               1        1.2   0.5   0.5          1
      3 white    3     0    NA               1        1.2   0.5   0.5          1
        segment.linetype segment.size segment.curvature segment.angle segment.ncp
      1                4          0.5                 0            90           1
      2                4          0.5                 0            90           1
      3                4          0.5                 0            90           1
        segment.shape segment.square segment.squareShape segment.inflect
      1           0.5           TRUE                   1           FALSE
      2           0.5           TRUE                   1           FALSE
      3           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      3         FALSE
      
      [[2]]
        x xend      y   yend                                  annotation group
      1 1    1 36.915 37.175     list(~italic(p)[Holm-corrected]==0.000) 4-f-1
      2 1    2 37.175 37.175     list(~italic(p)[Holm-corrected]==0.000) 4-f-1
      3 2    2 37.175 36.915     list(~italic(p)[Holm-corrected]==0.000) 4-f-1
      4 2    2 38.865 39.125 list(~italic(p)[Holm-corrected]==3.045e-08) f-r-2
      5 2    3 39.125 39.125 list(~italic(p)[Holm-corrected]==3.045e-08) f-r-2
      6 3    3 39.125 38.865 list(~italic(p)[Holm-corrected]==3.045e-08) f-r-2
        flipped_aes PANEL shape colour textsize angle hjust vjust alpha family
      1       FALSE     1    19  black        3     0   0.5     0    NA       
      2       FALSE     1    19  black        3     0   0.5     0    NA       
      3       FALSE     1    19  black        3     0   0.5     0    NA       
      4       FALSE     1    19  black        3     0   0.5     0    NA       
      5       FALSE     1    19  black        3     0   0.5     0    NA       
      6       FALSE     1    19  black        3     0   0.5     0    NA       
        fontface lineheight linetype size
      1        1        1.2        1  0.5
      2        1        1.2        1  0.5
      3        1        1.2        1  0.5
      4        1        1.2        1  0.5
      5        1        1.2        1  0.5
      6        1        1.2        1  0.5
      

---

    Code
      pb$plot$labels
    Output
      $x
      [1] "drv"
      
      $y
      [1] "cty"
      
      $colour
      [1] "drv"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      expression(list("Pairwise test:" ~ bold("Yuen's trimmed means test"), 
          "Comparisons shown:" ~ bold("only significant")))
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

# check bayesian test display

    Code
      list(pb$data[[6]], pb$data[[7]])
    Output
      [[1]]
        x        y                     label PANEL group nudge_x  nudge_y colour
      1 1 5.021408 widehat(mu)[MAP]=='5.021'     1     1     1.4 5.021408  black
      2 2 5.747801 widehat(mu)[MAP]=='5.748'     1     2     2.4 5.747801  black
      3 3 6.398534 widehat(mu)[MAP]=='6.399'     1     3     3.4 6.398534  black
         fill size angle alpha family fontface lineheight hjust vjust point.size
      1 white    3     0    NA               1        1.2   0.5   0.5          1
      2 white    3     0    NA               1        1.2   0.5   0.5          1
      3 white    3     0    NA               1        1.2   0.5   0.5          1
        segment.linetype segment.size segment.curvature segment.angle segment.ncp
      1                4          0.5                 0            90           1
      2                4          0.5                 0            90           1
      3                4          0.5                 0            90           1
        segment.shape segment.square segment.squareShape segment.inflect
      1           0.5           TRUE                   1           FALSE
      2           0.5           TRUE                   1           FALSE
      3           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      3         FALSE
      
      [[2]]
        x xend      y   yend                       annotation                  group
      1 1    1 8.2415 8.2775 list(~log[e](BF['01'])==-33.669)    setosa-versicolor-1
      2 1    2 8.2775 8.2775 list(~log[e](BF['01'])==-33.669)    setosa-versicolor-1
      3 2    2 8.2775 8.2415 list(~log[e](BF['01'])==-33.669)    setosa-versicolor-1
      4 1    1 8.5115 8.5475 list(~log[e](BF['01'])==-56.343)     setosa-virginica-2
      5 1    3 8.5475 8.5475 list(~log[e](BF['01'])==-56.343)     setosa-virginica-2
      6 3    3 8.5475 8.5115 list(~log[e](BF['01'])==-56.343)     setosa-virginica-2
      7 2    2 8.7815 8.8175 list(~log[e](BF['01'])==-11.162) versicolor-virginica-3
      8 2    3 8.8175 8.8175 list(~log[e](BF['01'])==-11.162) versicolor-virginica-3
      9 3    3 8.8175 8.7815 list(~log[e](BF['01'])==-11.162) versicolor-virginica-3
        flipped_aes PANEL shape colour textsize angle hjust vjust alpha family
      1       FALSE     1    19  black        3     0   0.5     0    NA       
      2       FALSE     1    19  black        3     0   0.5     0    NA       
      3       FALSE     1    19  black        3     0   0.5     0    NA       
      4       FALSE     1    19  black        3     0   0.5     0    NA       
      5       FALSE     1    19  black        3     0   0.5     0    NA       
      6       FALSE     1    19  black        3     0   0.5     0    NA       
      7       FALSE     1    19  black        3     0   0.5     0    NA       
      8       FALSE     1    19  black        3     0   0.5     0    NA       
      9       FALSE     1    19  black        3     0   0.5     0    NA       
        fontface lineheight linetype size
      1        1        1.2        1  0.5
      2        1        1.2        1  0.5
      3        1        1.2        1  0.5
      4        1        1.2        1  0.5
      5        1        1.2        1  0.5
      6        1        1.2        1  0.5
      7        1        1.2        1  0.5
      8        1        1.2        1  0.5
      9        1        1.2        1  0.5
      

---

    Code
      pb$plot$labels
    Output
      $x
      [1] "Species"
      
      $y
      [1] "Sepal.Length"
      
      $colour
      [1] "Species"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      expression(list("Pairwise test:" ~ bold("Student's t-test"), 
          "Comparisons shown:" ~ bold("all")))
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

# additional test

    Code
      length(pb1$data)
    Output
      [1] 6

---

    Code
      list(pb2$data[[6]], pb2$data[[7]])
    Output
      [[1]]
        x        y                       label PANEL group nudge_x  nudge_y colour
      1 1 116.2667 widehat(mu)[mean]=='116.27'     1     1     1.4 116.2667  black
      2 2 116.6944 widehat(mu)[mean]=='116.69'     1     2     2.4 116.6944  black
      3 3 102.3838 widehat(mu)[mean]=='102.38'     1     3     3.4 102.3838  black
         fill size angle alpha family fontface lineheight hjust vjust point.size
      1 white    3     0    NA               1        1.2   0.5   0.5          1
      2 white    3     0    NA               1        1.2   0.5   0.5          1
      3 white    3     0    NA               1        1.2   0.5   0.5          1
        segment.linetype segment.size segment.curvature segment.angle segment.ncp
      1                4          0.5                 0            90           1
      2                4          0.5                 0            90           1
      3                4          0.5                 0            90           1
        segment.shape segment.square segment.squareShape segment.inflect
      1           0.5           TRUE                   1           FALSE
      2           0.5           TRUE                   1           FALSE
      3           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      3         FALSE
      
      [[2]]
        x xend       y    yend                                 annotation     group
      1 2    2 264.115 265.825 list(~italic(p)[Holm-corrected]==1.27e-03) PG-13-R-1
      2 2    3 265.825 265.825 list(~italic(p)[Holm-corrected]==1.27e-03) PG-13-R-1
      3 3    3 265.825 264.115 list(~italic(p)[Holm-corrected]==1.27e-03) PG-13-R-1
        flipped_aes PANEL shape colour textsize angle hjust vjust alpha family
      1       FALSE     1    19  black        3     0   0.5     0    NA       
      2       FALSE     1    19  black        3     0   0.5     0    NA       
      3       FALSE     1    19  black        3     0   0.5     0    NA       
        fontface lineheight linetype size
      1        1        1.2        1  0.5
      2        1        1.2        1  0.5
      3        1        1.2        1  0.5
      

---

    Code
      pb1$plot$labels
    Output
      $x
      [1] "mpaa"
      
      $y
      [1] "length"
      
      $colour
      [1] "mpaa"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      expression(list("Pairwise test:" ~ bold("Games-Howell test"), 
          "Comparisons shown:" ~ bold("only significant")))
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

---

    Code
      pb2$plot$labels
    Output
      $x
      [1] "mpaa"
      
      $y
      [1] "length"
      
      $colour
      [1] "mpaa"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      expression(list("Pairwise test:" ~ bold("Games-Howell test"), 
          "Comparisons shown:" ~ bold("only significant")))
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

# checking one sample proportion test

    Code
      pb$data
    Output
      [[1]]
             fill         y x PANEL group flipped_aes      ymin      ymax xmin xmax
      1 #1B9E77FF 1.0000000 1     1     1       FALSE 0.7368421 1.0000000  0.5  1.5
      2 #D95F02FF 0.7368421 1     1     2       FALSE 0.6710526 0.7368421  0.5  1.5
      3 #7570B3FF 0.6710526 1     1     3       FALSE 0.2500000 0.6710526  0.5  1.5
      4 #E7298AFF 0.2500000 1     1     4       FALSE 0.0000000 0.2500000  0.5  1.5
        colour size linetype alpha
      1  black  0.5        1    NA
      2  black  0.5        1    NA
      3  black  0.5        1    NA
      4  black  0.5        1    NA
      
      [[2]]
                y x        label group PANEL      ymax xmin xmax      ymin colour
      1 0.8684211 1 20\n(26.32%)     1     1 1.0000000    1    1 0.7368421  black
      2 0.7039474 1   5\n(6.58%)     2     1 0.7368421    1    1 0.6710526  black
      3 0.4605263 1 32\n(42.11%)     3     1 0.6710526    1    1 0.2500000  black
      4 0.1250000 1    19\n(25%)     4     1 0.2500000    1    1 0.0000000  black
         fill size angle hjust vjust alpha family fontface lineheight
      1 white 3.88     0   0.5   0.5     1               1        1.2
      2 white 3.88     0   0.5   0.5     1               1        1.2
      3 white 3.88     0   0.5   0.5     1               1        1.2
      4 white 3.88     0   0.5   0.5     1               1        1.2
      

# checking labels with contingency tab

    Code
      list(pb$data, pb1$data)
    Output
      [[1]]
      [[1]][[1]]
             fill         y x PANEL group flipped_aes      ymin      ymax xmin xmax
      1 #9A8822FF 1.0000000 1     1     1       FALSE 0.2727273 1.0000000  0.5  1.5
      2 #F5CDB4FF 0.2727273 1     1     2       FALSE 0.0000000 0.2727273  0.5  1.5
      3 #9A8822FF 1.0000000 1     2     1       FALSE 0.5714286 1.0000000  0.5  1.5
      4 #F5CDB4FF 0.5714286 1     2     2       FALSE 0.0000000 0.5714286  0.5  1.5
      5 #9A8822FF 1.0000000 1     3     1       FALSE 0.8571429 1.0000000  0.5  1.5
      6 #F5CDB4FF 0.8571429 1     3     2       FALSE 0.0000000 0.8571429  0.5  1.5
        colour size linetype alpha
      1  black  0.5        1    NA
      2  black  0.5        1    NA
      3  black  0.5        1    NA
      4  black  0.5        1    NA
      5  black  0.5        1    NA
      6  black  0.5        1    NA
      
      [[1]][[2]]
                y x label group PANEL      ymax xmin xmax      ymin colour  fill size
      1 0.6363636 1     8     1     1 1.0000000    1    1 0.2727273  black white 3.88
      2 0.1363636 1     3     2     1 0.2727273    1    1 0.0000000  black white 3.88
      3 0.7857143 1     3     1     2 1.0000000    1    1 0.5714286  black white 3.88
      4 0.2857143 1     4     2     2 0.5714286    1    1 0.0000000  black white 3.88
      5 0.9285714 1     2     1     3 1.0000000    1    1 0.8571429  black white 3.88
      6 0.4285714 1    12     2     3 0.8571429    1    1 0.0000000  black white 3.88
        angle hjust vjust alpha family fontface lineheight
      1     0   0.5   0.5     1               1        1.2
      2     0   0.5   0.5     1               1        1.2
      3     0   0.5   0.5     1               1        1.2
      4     0   0.5   0.5     1               1        1.2
      5     0   0.5   0.5     1               1        1.2
      6     0   0.5   0.5     1               1        1.2
      
      [[1]][[3]]
        y    x                                                                 label
      1 1 1.65     list(~chi['gof']^2~(1)==2.27, ~italic(p)=='0.13', ~italic(n)==11)
      2 1 1.65      list(~chi['gof']^2~(1)==0.14, ~italic(p)=='0.71', ~italic(n)==7)
      3 1 1.65 list(~chi['gof']^2~(1)==7.14, ~italic(p)=='7.53e-03', ~italic(n)==14)
        PANEL group ymax xmin xmax ymin colour size angle hjust vjust alpha family
      1     1    -1    1 1.65 1.65    0  black  2.8     0   0.5   0.5    NA       
      2     2    -1    1 1.65 1.65    0  black  2.8     0   0.5   0.5    NA       
      3     3    -1    1 1.65 1.65    0  black  2.8     0   0.5   0.5    NA       
        fontface lineheight
      1        1        1.2
      2        1        1.2
      3        1        1.2
      
      
      [[2]]
      [[2]][[1]]
             fill         y x PANEL group flipped_aes      ymin      ymax xmin xmax
      1 #1B9E77FF 1.0000000 1     1     1       FALSE 0.3684211 1.0000000  0.5  1.5
      2 #D95F02FF 0.3684211 1     1     2       FALSE 0.1578947 0.3684211  0.5  1.5
      3 #7570B3FF 0.1578947 1     1     3       FALSE 0.0000000 0.1578947  0.5  1.5
        colour size linetype alpha
      1  black  0.5        1    NA
      2  black  0.5        1    NA
      3  black  0.5        1    NA
      
      [[2]][[2]]
                 y x label group PANEL      ymax xmin xmax      ymin colour  fill
      1 0.68421053 1   63%     1     1 1.0000000    1    1 0.3684211  black white
      2 0.26315789 1   21%     2     1 0.3684211    1    1 0.1578947  black white
      3 0.07894737 1   16%     3     1 0.1578947    1    1 0.0000000  black white
        size angle hjust vjust alpha family fontface lineheight
      1 3.88     0   0.5   0.5     1               1        1.2
      2 3.88     0   0.5   0.5     1               1        1.2
      3 3.88     0   0.5   0.5     1               1        1.2
      
      

---

    Code
      within(pb$plot$labels, rm(subtitle, caption))
    Output
      $x
      NULL
      
      $y
      NULL
      
      $title
      NULL
      
      $fill
      [1] "am"
      
      $label
      [1] ".label"
      
      $group
      [1] "am"
      
      $alt
      [1] ""
      

# checking labels with counts

    Code
      pb$data
    Output
      [[1]]
             fill         y x PANEL group flipped_aes      ymin      ymax xmin xmax
      1 #1B9E77FF 1.0000000 1     1     1       FALSE 0.9154362 1.0000000  0.5  1.5
      2 #D95F02FF 0.9154362 1     1     2       FALSE 0.0000000 0.9154362  0.5  1.5
      3 #1B9E77FF 1.0000000 1     2     1       FALSE 0.5161744 1.0000000  0.5  1.5
      4 #D95F02FF 0.5161744 1     2     2       FALSE 0.0000000 0.5161744  0.5  1.5
        colour size linetype alpha
      1  black  0.5        1    NA
      2  black  0.5        1    NA
      3  black  0.5        1    NA
      4  black  0.5        1    NA
      
      [[2]]
                y x  label group PANEL      ymax xmin xmax      ymin colour  fill
      1 0.9577181 1  8.46%     1     1 1.0000000    1    1 0.9154362  black white
      2 0.4577181 1 91.54%     2     1 0.9154362    1    1 0.0000000  black white
      3 0.7580872 1 48.38%     1     2 1.0000000    1    1 0.5161744  black white
      4 0.2580872 1 51.62%     2     2 0.5161744    1    1 0.0000000  black white
        size angle hjust vjust alpha family fontface lineheight
      1 3.88     0   0.5   0.5     1               1        1.2
      2 3.88     0   0.5   0.5     1               1        1.2
      3 3.88     0   0.5   0.5     1               1        1.2
      4 3.88     0   0.5   0.5     1               1        1.2
      
      [[3]]
        y    x
      1 1 1.65
      2 1 1.65
                                                                               label
      1 list(~chi['gof']^2~(1)==1028.62, ~italic(p)=='1.08e-225', ~italic(n)==1,490)
      2           list(~chi['gof']^2~(1)==0.74, ~italic(p)=='0.39', ~italic(n)==711)
        PANEL group ymax xmin xmax ymin colour size angle hjust vjust alpha family
      1     1    -1    1 1.65 1.65    0  black  2.8     0   0.5   0.5    NA       
      2     2    -1    1 1.65 1.65    0  black  2.8     0   0.5   0.5    NA       
        fontface lineheight
      1        1        1.2
      2        1        1.2
      

---

    Code
      within(pb$plot$labels, rm(subtitle))
    Output
      $x
      NULL
      
      $y
      NULL
      
      $title
      NULL
      
      $caption
      NULL
      
      $fill
      [1] "Sex"
      
      $label
      [1] ".label"
      
      $group
      [1] "Sex"
      
      $alt
      [1] ""
      

# checking labels with contingency tab (paired)

    Code
      within(pb$plot$labels, rm(subtitle))
    Output
      $x
      NULL
      
      $y
      NULL
      
      $title
      NULL
      
      $caption
      NULL
      
      $fill
      [1] "1st survey"
      
      $label
      [1] ".label"
      
      $group
      [1] "1st survey"
      
      $alt
      [1] ""
      

---

    Code
      pb$data[[3]]
    Output
        y    x
      1 1 1.65
      2 1 1.65
                                                                            label
      1 list(~chi['gof']^2~(1)==569.62, ~italic(p)=='6.80e-126', ~italic(n)==880)
      2  list(~chi['gof']^2~(1)==245.00, ~italic(p)=='3.20e-55', ~italic(n)==720)
        PANEL group ymax xmin xmax ymin colour size angle hjust vjust alpha family
      1     1    -1    1 1.65 1.65    0  black  2.8     0   0.5   0.5    NA       
      2     2    -1    1 1.65 1.65    0  black  2.8     0   0.5   0.5    NA       
        fontface lineheight
      1        1        1.2
      2        1        1.2

# repelling labels

    Code
      pb$data
    Output
      [[1]]
             fill         y x PANEL group flipped_aes      ymin      ymax xmin xmax
      1 #1B9E77FF 1.0000000 1     1     1       FALSE 0.8479021 1.0000000  0.5  1.5
      2 #D95F02FF 0.8479021 1     1     2       FALSE 0.8128739 0.8479021  0.5  1.5
      3 #7570B3FF 0.8128739 1     1     3       FALSE 0.6511922 0.8128739  0.5  1.5
      4 #E7298AFF 0.6511922 1     1     4       FALSE 0.0000000 0.6511922  0.5  1.5
      5 #1B9E77FF 1.0000000 1     2     1       FALSE 0.9917011 1.0000000  0.5  1.5
      6 #D95F02FF 0.9917011 1     2     2       FALSE 0.9722386 0.9917011  0.5  1.5
      7 #7570B3FF 0.9722386 1     2     3       FALSE 0.8895863 0.9722386  0.5  1.5
      8 #E7298AFF 0.8895863 1     2     4       FALSE 0.0000000 0.8895863  0.5  1.5
        colour size linetype alpha
      1  black  0.5        1    NA
      2  black  0.5        1    NA
      3  black  0.5        1    NA
      4  black  0.5        1    NA
      5  black  0.5        1    NA
      6  black  0.5        1    NA
      7  black  0.5        1    NA
      8  black  0.5        1    NA
      
      [[2]]
                y x label group PANEL      ymax xmin xmax      ymin colour  fill size
      1 0.9239510 1   15%     1     1 1.0000000    1    1 0.8479021  black white 3.88
      2 0.8303880 1    4%     2     1 0.8479021    1    1 0.8128739  black white 3.88
      3 0.7320330 1   16%     3     1 0.8128739    1    1 0.6511922  black white 3.88
      4 0.3255961 1   65%     4     1 0.6511922    1    1 0.0000000  black white 3.88
      5 0.9958505 1    1%     1     2 1.0000000    1    1 0.9917011  black white 3.88
      6 0.9819698 1    2%     2     2 0.9917011    1    1 0.9722386  black white 3.88
      7 0.9309125 1    8%     3     2 0.9722386    1    1 0.8895863  black white 3.88
      8 0.4447932 1   89%     4     2 0.8895863    1    1 0.0000000  black white 3.88
        angle alpha family fontface lineheight hjust vjust point.size
      1     0     1               1        1.2   0.5   0.5          1
      2     0     1               1        1.2   0.5   0.5          1
      3     0     1               1        1.2   0.5   0.5          1
      4     0     1               1        1.2   0.5   0.5          1
      5     0     1               1        1.2   0.5   0.5          1
      6     0     1               1        1.2   0.5   0.5          1
      7     0     1               1        1.2   0.5   0.5          1
      8     0     1               1        1.2   0.5   0.5          1
        segment.linetype segment.size segment.curvature segment.angle segment.ncp
      1                1          0.5                 0            90           1
      2                1          0.5                 0            90           1
      3                1          0.5                 0            90           1
      4                1          0.5                 0            90           1
      5                1          0.5                 0            90           1
      6                1          0.5                 0            90           1
      7                1          0.5                 0            90           1
      8                1          0.5                 0            90           1
        segment.shape segment.square segment.squareShape segment.inflect
      1           0.5           TRUE                   1           FALSE
      2           0.5           TRUE                   1           FALSE
      3           0.5           TRUE                   1           FALSE
      4           0.5           TRUE                   1           FALSE
      5           0.5           TRUE                   1           FALSE
      6           0.5           TRUE                   1           FALSE
      7           0.5           TRUE                   1           FALSE
      8           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      3         FALSE
      4         FALSE
      5         FALSE
      6         FALSE
      7         FALSE
      8         FALSE
      

---

    Code
      pb$plot$labels
    Output
      $x
      NULL
      
      $y
      NULL
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      NULL
      
      $fill
      [1] "mode"
      
      $label
      [1] ".label"
      
      $group
      [1] "mode"
      
      $alt
      [1] ""
      

# `pairwise_comparisons()` works for between-subjects design

    Code
      list(df1, df2, df3, df4, df5)
    Output
      [[1]]
      # A tibble: 6 x 6
        group1  group2  p.value test.details     p.value.adjustment
        <chr>   <chr>     <dbl> <chr>            <chr>             
      1 carni   herbi     1     Student's t-test Bonferroni        
      2 carni   insecti   1     Student's t-test Bonferroni        
      3 carni   omni      1     Student's t-test Bonferroni        
      4 herbi   insecti   1     Student's t-test Bonferroni        
      5 herbi   omni      0.979 Student's t-test Bonferroni        
      6 insecti omni      1     Student's t-test Bonferroni        
        label                                       
        <chr>                                       
      1 list(~italic(p)[Bonferroni-corrected]==1.00)
      2 list(~italic(p)[Bonferroni-corrected]==1.00)
      3 list(~italic(p)[Bonferroni-corrected]==1.00)
      4 list(~italic(p)[Bonferroni-corrected]==1.00)
      5 list(~italic(p)[Bonferroni-corrected]==0.98)
      6 list(~italic(p)[Bonferroni-corrected]==1.00)
      
      [[2]]
      # A tibble: 6 x 11
        group1  group2  statistic p.value alternative method            distribution
        <chr>   <chr>       <dbl>   <dbl> <chr>       <chr>             <chr>       
      1 carni   herbi        2.17       1 two.sided   Games-Howell test q           
      2 carni   insecti     -2.17       1 two.sided   Games-Howell test q           
      3 carni   omni         1.10       1 two.sided   Games-Howell test q           
      4 herbi   insecti     -2.41       1 two.sided   Games-Howell test q           
      5 herbi   omni        -1.87       1 two.sided   Games-Howell test q           
      6 insecti omni         2.19       1 two.sided   Games-Howell test q           
        p.adjustment test.details      p.value.adjustment
        <chr>        <chr>             <chr>             
      1 none         Games-Howell test Bonferroni        
      2 none         Games-Howell test Bonferroni        
      3 none         Games-Howell test Bonferroni        
      4 none         Games-Howell test Bonferroni        
      5 none         Games-Howell test Bonferroni        
      6 none         Games-Howell test Bonferroni        
        label                                       
        <chr>                                       
      1 list(~italic(p)[Bonferroni-corrected]==1.00)
      2 list(~italic(p)[Bonferroni-corrected]==1.00)
      3 list(~italic(p)[Bonferroni-corrected]==1.00)
      4 list(~italic(p)[Bonferroni-corrected]==1.00)
      5 list(~italic(p)[Bonferroni-corrected]==1.00)
      6 list(~italic(p)[Bonferroni-corrected]==1.00)
      
      [[3]]
      # A tibble: 6 x 11
        group1  group2  statistic p.value alternative method               
        <chr>   <chr>       <dbl>   <dbl> <chr>       <chr>                
      1 carni   herbi       0.582  0.561  two.sided   Dunn's all-pairs test
      2 carni   insecti     1.88   0.0595 two.sided   Dunn's all-pairs test
      3 carni   omni        1.14   0.254  two.sided   Dunn's all-pairs test
      4 herbi   insecti     1.63   0.102  two.sided   Dunn's all-pairs test
      5 herbi   omni        0.717  0.474  two.sided   Dunn's all-pairs test
      6 insecti omni        1.14   0.254  two.sided   Dunn's all-pairs test
        distribution p.adjustment test.details p.value.adjustment
        <chr>        <chr>        <chr>        <chr>             
      1 z            none         Dunn test    None              
      2 z            none         Dunn test    None              
      3 z            none         Dunn test    None              
      4 z            none         Dunn test    None              
      5 z            none         Dunn test    None              
      6 z            none         Dunn test    None              
        label                              
        <chr>                              
      1 list(~italic(p)[uncorrected]==0.56)
      2 list(~italic(p)[uncorrected]==0.06)
      3 list(~italic(p)[uncorrected]==0.25)
      4 list(~italic(p)[uncorrected]==0.10)
      5 list(~italic(p)[uncorrected]==0.47)
      6 list(~italic(p)[uncorrected]==0.25)
      
      [[4]]
      # A tibble: 6 x 10
        group1  group2  estimate conf.level conf.low conf.high p.value
        <chr>   <chr>      <dbl>      <dbl>    <dbl>     <dbl>   <dbl>
      1 carni   herbi   -0.0323        0.95  -0.248     0.184    0.790
      2 carni   insecti  0.0451        0.95  -0.0484    0.139    0.552
      3 carni   omni     0.00520       0.95  -0.114     0.124    0.898
      4 herbi   insecti  0.0774        0.95  -0.133     0.288    0.552
      5 herbi   omni     0.0375        0.95  -0.182     0.257    0.790
      6 insecti omni    -0.0399        0.95  -0.142     0.0625   0.552
        test.details              p.value.adjustment
        <chr>                     <chr>             
      1 Yuen's trimmed means test FDR               
      2 Yuen's trimmed means test FDR               
      3 Yuen's trimmed means test FDR               
      4 Yuen's trimmed means test FDR               
      5 Yuen's trimmed means test FDR               
      6 Yuen's trimmed means test FDR               
        label                                
        <chr>                                
      1 list(~italic(p)[FDR-corrected]==0.79)
      2 list(~italic(p)[FDR-corrected]==0.55)
      3 list(~italic(p)[FDR-corrected]==0.90)
      4 list(~italic(p)[FDR-corrected]==0.55)
      5 list(~italic(p)[FDR-corrected]==0.79)
      6 list(~italic(p)[FDR-corrected]==0.55)
      
      [[5]]
      # A tibble: 3 x 6
        group1 group2 p.value test.details     p.value.adjustment
        <chr>  <chr>    <dbl> <chr>            <chr>             
      1 PG     PG-13  0.316   Student's t-test Holm              
      2 PG     R      0.00283 Student's t-test Holm              
      3 PG-13  R      0.00310 Student's t-test Holm              
        label                                     
        <chr>                                     
      1 list(~italic(p)[Holm-corrected]==0.32)    
      2 list(~italic(p)[Holm-corrected]==2.83e-03)
      3 list(~italic(p)[Holm-corrected]==3.10e-03)
      

# dropped levels are not included

    Code
      df2$label
    Output
      [1] "list(~italic(p)[uncorrected]==0.87)"

# data without NAs

    Code
      df$label
    Output
      [1] "list(~italic(p)[FDR-corrected]==1.32e-15)"
      [2] "list(~italic(p)[FDR-corrected]==6.64e-32)"
      [3] "list(~italic(p)[FDR-corrected]==2.77e-09)"

# ggcoefstats with glm with z

    Code
      pb$data
    Output
      [[1]]
                 x y PANEL group shape colour size fill alpha stroke
      1 -0.7800447 1     1     1    19   blue    3   NA    NA    0.5
      2  2.2940067 2     1     2    19   blue    3   NA    NA    0.5
      3 -0.5564393 3     1     3    19   blue    3   NA    NA    0.5
      
      [[2]]
                 x       xmin       xmax y PANEL group ymin ymax colour size linetype
      1 -0.7800447 -1.1521325 -0.4124066 1     1     1    1    1  black  0.5        1
      2  2.2940067  2.0985524  2.4932922 2     1     2    2    2  black  0.5        1
      3 -0.5564393 -0.9299989 -0.1804338 3     1     3    3    3  black  0.5        1
        height alpha
      1      0    NA
      2      0    NA
      3      0    NA
      
      [[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      [[4]]
                 x y
      1 -0.7800447 1
      2  2.2940067 2
      3 -0.5564393 3
                                                                                  label
      1 list(widehat(italic(beta))=='-0.78', italic(z)=='-3.47', italic(p)=='5.14e-04')
      2  list(widehat(italic(beta))=='2.29', italic(z)=='19.13', italic(p)=='1.54e-81')
      3     list(widehat(italic(beta))=='-0.56', italic(z)=='-2.44', italic(p)=='0.01')
        PANEL group    colour  fill size angle alpha family fontface lineheight hjust
      1     1     1 #1B9E77FF white    3     0    NA               1        1.2   0.5
      2     1     2 #D95F02FF white    3     0    NA               1        1.2   0.5
      3     1     3 #7570B3FF white    3     0    NA               1        1.2   0.5
        vjust point.size segment.linetype segment.size segment.curvature
      1   0.5          1                1          0.5                 0
      2   0.5          1                1          0.5                 0
      3   0.5          1                1          0.5                 0
        segment.angle segment.ncp segment.shape segment.square segment.squareShape
      1            90           1           0.5           TRUE                   1
      2            90           1           0.5           TRUE                   1
      3            90           1           0.5           TRUE                   1
        segment.inflect segment.debug
      1           FALSE         FALSE
      2           FALSE         FALSE
      3           FALSE         FALSE
      

# ggcoefstats with chi-squared statistic model

    Code
      pb$data
    Output
      [[1]]
                  x y PANEL group shape colour size fill alpha stroke
      1  0.01703351 1     1     1    19   blue    3   NA    NA    0.5
      2 -0.51166834 2     1     2    19   blue    3   NA    NA    0.5
      
      [[2]]
                  x         xmin       xmax y PANEL group ymin ymax colour size
      1  0.01703351 -0.001062183  0.0351292 1     1     1    1    1  black  0.5
      2 -0.51166834 -0.840312344 -0.1830243 2     1     2    2    2  black  0.5
        linetype height alpha
      1        1      0    NA
      2        1      0    NA
      
      [[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      [[4]]
                  x y
      1  0.01703351 1
      2 -0.51166834 2
                                                                                           label
      1      list(widehat(italic(beta))=='0.02', italic(chi)^2*('1')=='3.40', italic(p)=='0.07')
      2 list(widehat(italic(beta))=='-0.51', italic(chi)^2*('1')=='9.31', italic(p)=='2.28e-03')
        PANEL group    colour  fill size angle alpha family fontface lineheight hjust
      1     1     1 #3182BDFF white    3     0    NA               1        1.2   0.5
      2     1     2 #E6550DFF white    3     0    NA               1        1.2   0.5
        vjust point.size segment.linetype segment.size segment.curvature
      1   0.5          1                1          0.5                 0
      2   0.5          1                1          0.5                 0
        segment.angle segment.ncp segment.shape segment.square segment.squareShape
      1            90           1           0.5           TRUE                   1
      2            90           1           0.5           TRUE                   1
        segment.inflect segment.debug
      1           FALSE         FALSE
      2           FALSE         FALSE
      

# ggcoefstats with lm model

    Code
      pb$data
    Output
      [[1]]
                  x y PANEL group shape colour size fill alpha stroke
      1 -0.15565484 1     1     1    19   blue    3   NA    NA    0.5
      2 -1.80872181 2     1     2    19   blue    3   NA    NA    0.5
      3  0.06471454 3     1     3    19   blue    3   NA    NA    0.5
      
      [[2]]
                  x        xmin        xmax y PANEL group ymin ymax colour size
      1 -0.15565484 -0.22929853 -0.08201116 1     1     1    1    1  black  0.5
      2 -1.80872181 -3.71964130  0.10219768 2     1     2    2    2  black  0.5
      3  0.06471454 -0.02784951  0.15727859 3     1     3    3    3  black  0.5
        linetype height alpha
      1        1      0    NA
      2        1      0    NA
      3        1      0    NA
      
      [[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      [[4]]
                  x y
      1 -0.15565484 1
      2 -1.80872181 2
      3  0.06471454 3
                                                                                           label
      1 list(widehat(italic(beta))=='-0.156', italic(t)('28')=='-5.840', italic(p)=='2.813e-06')
      2     list(widehat(italic(beta))=='-1.809', italic(t)('28')=='-2.615', italic(p)=='0.014')
      3                                                                                     <NA>
        PANEL group    colour  fill size angle alpha family fontface lineheight hjust
      1     1     1 #1B9E77FF white    3     0    NA               1        1.2   0.5
      2     1     2 #D95F02FF white    3     0    NA               1        1.2   0.5
      3     1     3 #7570B3FF white    3     0    NA               1        1.2   0.5
        vjust point.size segment.linetype segment.size segment.curvature
      1   0.5          1                1          0.5                 0
      2   0.5          1                1          0.5                 0
      3   0.5          1                1          0.5                 0
        segment.angle segment.ncp segment.shape segment.square segment.squareShape
      1            90           1           0.5           TRUE                   1
      2            90           1           0.5           TRUE                   1
      3            90           1           0.5           TRUE                   1
        segment.inflect segment.debug
      1           FALSE         FALSE
      2           FALSE         FALSE
      3           FALSE         FALSE
      

# ggcoefstats with partial variants of effect size for f-statistic

    Code
      list(tidy_df1, p$labels)
    Output
      [[1]]
      # A tibble: 3 x 14
        term   statistic    df df.error  p.value  sumsq meansq estimate conf.low
        <fct>      <dbl> <dbl>    <dbl>    <dbl>  <dbl>  <dbl>    <dbl>    <dbl>
      1 mpg       119.       1       28 1.38e-11 22.3   22.3      0.809   0.693 
      2 am          7.30     1       28 1.16e- 2  1.37   1.37     0.207   0.0299
      3 mpg:am      3.73     1       28 6.36e- 2  0.701  0.701    0.118   0     
        conf.high sum.squares.error mean.square.error effectsize         
            <dbl>             <dbl>             <dbl> <chr>              
      1         1              5.26              5.26 partial eta-squared
      2         1              5.26              5.26 partial eta-squared
      3         1              5.26              5.26 partial eta-squared
        label                                                                         
        <glue>                                                                        
      1 list(widehat(italic(eta)[p]^2)=='0.81', italic(F)('1', '28')=='118.89', itali~
      2 list(widehat(italic(eta)[p]^2)=='0.21', italic(F)('1', '28')=='7.30', italic(~
      3 list(widehat(italic(eta)[p]^2)=='0.12', italic(F)('1', '28')=='3.73', italic(~
      
      [[2]]
      [[2]]$x
      [1] "estimate"
      
      [[2]]$y
      [1] "effect"
      
      [[2]]$title
      NULL
      
      [[2]]$subtitle
      NULL
      
      [[2]]$caption
      atop(displaystyle(NULL), expr = paste("AIC = ", "43", ", BIC = ", 
          "50"))
      
      [[2]]$xmin
      [1] "conf.low"
      
      [[2]]$xmax
      [1] "conf.high"
      
      [[2]]$xintercept
      [1] "xintercept"
      
      [[2]]$label
      [1] "label"
      
      

---

    Code
      list(tidy_df2, p$labels)
    Output
      [[1]]
      # A tibble: 3 x 14
        term         statistic    df df.error  p.value sumsq meansq estimate conf.low
        <fct>            <dbl> <dbl>    <dbl>    <dbl> <dbl>  <dbl>    <dbl>    <dbl>
      1 vore              7.39     3       35 0.000584 19.6    6.54   0.308    0.0763
      2 brainwt           2.03     1       35 0.163     1.80   1.80   0.0235   0     
      3 vore:brainwt      4.01     3       35 0.0148   10.7    3.55   0.174    0     
        conf.high sum.squares.error mean.square.error effectsize           
            <dbl>             <dbl>             <dbl> <chr>                
      1         1              31.0              31.0 partial omega-squared
      2         1              31.0              31.0 partial omega-squared
      3         1              31.0              31.0 partial omega-squared
        label                                                                         
        <glue>                                                                        
      1 list(widehat(italic(omega)[p]^2)=='0.308', italic(F)('3', '35')=='7.388', ita~
      2 list(widehat(italic(omega)[p]^2)=='0.023', italic(F)('1', '35')=='2.034', ita~
      3 list(widehat(italic(omega)[p]^2)=='0.174', italic(F)('3', '35')=='4.012', ita~
      
      [[2]]
      [[2]]$x
      [1] "estimate"
      
      [[2]]$y
      [1] "term"
      
      [[2]]$title
      [1] "mammalian sleep"
      
      [[2]]$subtitle
      [1] "Source: `{ggplot2}` package"
      
      [[2]]$caption
      atop(displaystyle(paste(italic("Note"), ": From `tidyverse`")), 
          expr = paste("AIC = ", "126", ", BIC = ", "142"))
      
      [[2]]$xmin
      [1] "conf.low"
      
      [[2]]$xmax
      [1] "conf.high"
      
      [[2]]$xintercept
      [1] "xintercept"
      
      [[2]]$label
      [1] "label"
      
      

# check tidy output

    Code
      list(tidy_df1, tidy_df2)
    Output
      [[1]]
      # A tibble: 7 x 15
        term  statistic    df df.error p.value group    sumsq  meansq estimate
        <fct>     <dbl> <dbl>    <dbl>   <dbl> <chr>    <dbl>   <dbl>    <dbl>
      1 N       12.3        1       12 0.00437 Within 189.    189.     0.505  
      2 P        0.544      1       12 0.475   Within   8.40    8.40   0.0434 
      3 K        6.17       1       12 0.0288  Within  95.2    95.2    0.339  
      4 N:P      1.38       1       12 0.263   Within  21.3    21.3    0.103  
      5 N:K      2.15       1       12 0.169   Within  33.1    33.1    0.152  
      6 P:K      0.0312     1       12 0.863   Within   0.482   0.482  0.00259
      7 N:P:K    0.483      1        4 0.525   block   37.0    37.0    0.108  
        conf.low conf.high sum.squares.error mean.square.error effectsize         
           <dbl>     <dbl>             <dbl>             <dbl> <chr>              
      1   0.145          1              185.              185. partial eta-squared
      2   0              1              185.              185. partial eta-squared
      3   0.0254         1              185.              185. partial eta-squared
      4   0              1              185.              185. partial eta-squared
      5   0              1              185.              185. partial eta-squared
      6   0              1              185.              185. partial eta-squared
      7   0              1              306.              306. partial eta-squared
        label                                                                         
        <glue>                                                                        
      1 list(widehat(italic(eta)[p]^2)=='0.51', italic(F)('1', '12')=='12.26', italic~
      2 list(widehat(italic(eta)[p]^2)=='0.04', italic(F)('1', '12')=='0.54', italic(~
      3 list(widehat(italic(eta)[p]^2)=='0.34', italic(F)('1', '12')=='6.17', italic(~
      4 list(widehat(italic(eta)[p]^2)=='0.10', italic(F)('1', '12')=='1.38', italic(~
      5 list(widehat(italic(eta)[p]^2)=='0.15', italic(F)('1', '12')=='2.15', italic(~
      6 list(widehat(italic(eta)[p]^2)=='2.59e-03', italic(F)('1', '12')=='0.03', ita~
      7 list(widehat(italic(eta)[p]^2)=='0.11', italic(F)('1', '4')=='0.48', italic(p~
      
      [[2]]
      # A tibble: 7 x 14
        term  statistic    df df.error p.value   sumsq  meansq estimate conf.low
        <fct>     <dbl> <dbl>    <dbl>   <dbl>   <dbl>   <dbl>    <dbl>    <dbl>
      1 N        6.16       1       16  0.0245 189.    189.    0.278      0.0241
      2 P        0.273      1       16  0.608    8.40    8.40  0.0168     0     
      3 K        3.10       1       16  0.0975  95.2    95.2   0.162      0     
      4 N:P      0.693      1       16  0.418   21.3    21.3   0.0415     0     
      5 N:K      1.08       1       16  0.314   33.1    33.1   0.0631     0     
      6 P:K      0.0157     1       16  0.902    0.482   0.482 0.000979   0     
      7 N:P:K    1.20       1       16  0.289   37.0    37.0   0.0700     0     
        conf.high sum.squares.error mean.square.error effectsize         
            <dbl>             <dbl>             <dbl> <chr>              
      1         1              492.              492. partial eta-squared
      2         1              492.              492. partial eta-squared
      3         1              492.              492. partial eta-squared
      4         1              492.              492. partial eta-squared
      5         1              492.              492. partial eta-squared
      6         1              492.              492. partial eta-squared
      7         1              492.              492. partial eta-squared
        label                                                                         
        <glue>                                                                        
      1 list(widehat(italic(eta)[p]^2)=='0.28', italic(F)('1', '16')=='6.16', italic(~
      2 list(widehat(italic(eta)[p]^2)=='0.02', italic(F)('1', '16')=='0.27', italic(~
      3 list(widehat(italic(eta)[p]^2)=='0.16', italic(F)('1', '16')=='3.10', italic(~
      4 list(widehat(italic(eta)[p]^2)=='0.04', italic(F)('1', '16')=='0.69', italic(~
      5 list(widehat(italic(eta)[p]^2)=='0.06', italic(F)('1', '16')=='1.08', italic(~
      6 list(widehat(italic(eta)[p]^2)=='9.79e-04', italic(F)('1', '16')=='0.02', ita~
      7 list(widehat(italic(eta)[p]^2)=='0.07', italic(F)('1', '16')=='1.20', italic(~
      

# duplicated terms

    Code
      pb$data
    Output
      [[1]]
                x y PANEL group shape colour size fill alpha stroke
      1 29.322072 1     1     1    19   blue    3   NA    NA    0.5
      2  1.124451 2     1     2    19   blue    3   NA    NA    0.5
      3 29.954761 3     1     3    19   blue    3   NA    NA    0.5
      4  1.182257 4     1     4    19   blue    3   NA    NA    0.5
      5 30.628379 5     1     5    19   blue    3   NA    NA    0.5
      6  1.251657 6     1     6    19   blue    3   NA    NA    0.5
      
      [[2]]
                x       xmin      xmax y PANEL group ymin ymax colour size linetype
      1 29.322072 29.0912441 29.552899 1     1     1    1    1  black  0.5        1
      2  1.124451  0.6227402  1.626161 2     1     2    2    2  black  0.5        1
      3 29.954761 29.7319808 30.177540 3     1     3    3    3  black  0.5        1
      4  1.182257  0.8641108  1.500404 4     1     4    4    4  black  0.5        1
      5 30.628379 30.4324684 30.824290 5     1     5    5    5  black  0.5        1
      6  1.251657  0.8884855  1.614829 6     1     6    6    6  black  0.5        1
        height alpha
      1      0    NA
      2      0    NA
      3      0    NA
      4      0    NA
      5      0    NA
      6      0    NA
      
      [[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      [[4]]
                x y
      1 29.322072 1
      2  1.124451 2
      3 29.954761 3
      4  1.182257 4
      5 30.628379 5
      6  1.251657 6
                                                                                   label
      1 list(widehat(italic(beta))=='29.32', italic(z)=='249.58', italic(p)=='9.84e-78')
      2    list(widehat(italic(beta))=='1.12', italic(z)=='4.40', italic(p)=='5.77e-05')
      3 list(widehat(italic(beta))=='29.95', italic(z)=='264.18', italic(p)=='6.08e-79')
      4    list(widehat(italic(beta))=='1.18', italic(z)=='7.30', italic(p)=='2.27e-09')
      5 list(widehat(italic(beta))=='30.63', italic(z)=='307.16', italic(p)=='3.78e-82')
      6    list(widehat(italic(beta))=='1.25', italic(z)=='6.77', italic(p)=='1.50e-08')
        PANEL group    colour  fill size angle alpha family fontface lineheight hjust
      1     1     1 #1B9E77FF white    3     0    NA               1        1.2   0.5
      2     1     2 #D95F02FF white    3     0    NA               1        1.2   0.5
      3     1     3 #7570B3FF white    3     0    NA               1        1.2   0.5
      4     1     4 #E7298AFF white    3     0    NA               1        1.2   0.5
      5     1     5 #66A61EFF white    3     0    NA               1        1.2   0.5
      6     1     6 #E6AB02FF white    3     0    NA               1        1.2   0.5
        vjust point.size segment.linetype segment.size segment.curvature
      1   0.5          1                1          0.5                 0
      2   0.5          1                1          0.5                 0
      3   0.5          1                1          0.5                 0
      4   0.5          1                1          0.5                 0
      5   0.5          1                1          0.5                 0
      6   0.5          1                1          0.5                 0
        segment.angle segment.ncp segment.shape segment.square segment.squareShape
      1            90           1           0.5           TRUE                   1
      2            90           1           0.5           TRUE                   1
      3            90           1           0.5           TRUE                   1
      4            90           1           0.5           TRUE                   1
      5            90           1           0.5           TRUE                   1
      6            90           1           0.5           TRUE                   1
        segment.inflect segment.debug
      1           FALSE         FALSE
      2           FALSE         FALSE
      3           FALSE         FALSE
      4           FALSE         FALSE
      5           FALSE         FALSE
      6           FALSE         FALSE
      

# ggcoefstats works with data frames

    Code
      list(pb1$data, pb2$data, pb3$data, pb4$data, pb6$data, pb7$data)
    Output
      [[1]]
      [[1]][[1]]
             x y PANEL group shape colour size fill alpha stroke
      1 0.0665 1     1     1    19   blue    3   NA    NA    0.5
      2 0.5420 2     1     2    19   blue    3   NA    NA    0.5
      3 0.0450 3     1     3    19   blue    3   NA    NA    0.5
      
      [[1]][[2]]
             x   xmin  xmax y PANEL group ymin ymax colour size linetype height alpha
      1 0.0665 -0.778 0.911 1     1     1    1    1  black  0.5        1      0    NA
      2 0.5420 -0.280 1.360 2     1     2    2    2  black  0.5        1      0    NA
      3 0.0450  0.030 0.650 3     1     3    3    3  black  0.5        1      0    NA
      
      [[1]][[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      [[1]][[4]]
             x y
      1 0.0665 1
      2 0.5420 2
      3 0.0450 3
                                                                                      label
      1      list(widehat(italic(beta))=='0.07', italic(t)('5')=='0.16', italic(p)=='0.88')
      2     list(widehat(italic(beta))=='0.54', italic(t)('10')=='1.33', italic(p)=='0.19')
      3 list(widehat(italic(beta))=='0.04', italic(t)('12')=='1.24', italic(p)=='1.00e-03')
        PANEL group    colour  fill size angle alpha family fontface lineheight hjust
      1     1     1 #1B9E77FF white    3     0    NA               1        1.2   0.5
      2     1     2 #D95F02FF white    3     0    NA               1        1.2   0.5
      3     1     3 #7570B3FF white    3     0    NA               1        1.2   0.5
        vjust point.size segment.linetype segment.size segment.curvature
      1   0.5          1                1          0.5                 0
      2   0.5          1                1          0.5                 0
      3   0.5          1                1          0.5                 0
        segment.angle segment.ncp segment.shape segment.square segment.squareShape
      1            90           1           0.5           TRUE                   1
      2            90           1           0.5           TRUE                   1
      3            90           1           0.5           TRUE                   1
        segment.inflect segment.debug
      1           FALSE         FALSE
      2           FALSE         FALSE
      3           FALSE         FALSE
      
      
      [[2]]
      [[2]][[1]]
             x y PANEL group shape colour size fill alpha stroke
      1 0.0665 2     1     2    19   blue    3   NA    NA    0.5
      2 0.5420 1     1     1    19   blue    3   NA    NA    0.5
      3 0.0450 3     1     3    19   blue    3   NA    NA    0.5
      
      [[2]][[2]]
             x   xmin  xmax y PANEL group ymin ymax colour size linetype height alpha
      1 0.0665 -0.778 0.911 2     1     2    2    2  black  0.5        1      0    NA
      2 0.5420 -0.280 1.360 1     1     1    1    1  black  0.5        1      0    NA
      3 0.0450  0.030 0.650 3     1     3    3    3  black  0.5        1      0    NA
      
      [[2]][[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      [[2]][[4]]
             x y
      1 0.0665 2
      2 0.5420 1
      3 0.0450 3
                                                                                label
      1     list(widehat(italic(beta))=='0.07', italic(z)=='0.16', italic(p)=='0.88')
      2     list(widehat(italic(beta))=='0.54', italic(z)=='1.33', italic(p)=='0.19')
      3 list(widehat(italic(beta))=='0.04', italic(z)=='1.24', italic(p)=='1.00e-03')
        PANEL group    colour  fill size angle alpha family fontface lineheight hjust
      1     1     2 #1B9E77FF white    3     0    NA               1        1.2   0.5
      2     1     1 #D95F02FF white    3     0    NA               1        1.2   0.5
      3     1     3 #7570B3FF white    3     0    NA               1        1.2   0.5
        vjust point.size segment.linetype segment.size segment.curvature
      1   0.5          1                1          0.5                 0
      2   0.5          1                1          0.5                 0
      3   0.5          1                1          0.5                 0
        segment.angle segment.ncp segment.shape segment.square segment.squareShape
      1            90           1           0.5           TRUE                   1
      2            90           1           0.5           TRUE                   1
      3            90           1           0.5           TRUE                   1
        segment.inflect segment.debug
      1           FALSE         FALSE
      2           FALSE         FALSE
      3           FALSE         FALSE
      
      
      [[3]]
      [[3]][[1]]
             x y PANEL group shape colour size fill alpha stroke
      1 0.0665 1     1     1    19   blue    3   NA    NA    0.5
      2 0.5420 2     1     2    19   blue    3   NA    NA    0.5
      3 0.0450 3     1     3    19   blue    3   NA    NA    0.5
      
      [[3]][[2]]
             x   xmin  xmax y PANEL group ymin ymax colour size linetype height alpha
      1 0.0665 -0.778 0.911 1     1     1    1    1  black  0.5        1      0    NA
      2 0.5420 -0.280 1.360 2     1     2    2    2  black  0.5        1      0    NA
      3 0.0450  0.030 0.650 3     1     3    3    3  black  0.5        1      0    NA
      
      [[3]][[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      
      [[4]]
      [[4]][[1]]
        y      x PANEL group shape colour size fill alpha stroke
      1 1 0.0665     1     1    19   blue    3   NA    NA    0.5
      2 2 0.5420     1     2    19   blue    3   NA    NA    0.5
      3 3 0.0450     1     3    19   blue    3   NA    NA    0.5
      
      [[4]][[2]]
        y      x   xmin  xmax PANEL group ymin ymax colour size linetype height alpha
      1 1 0.0665 -0.778 0.911     1     1    1    1  black  0.5        1      0    NA
      2 2 0.5420 -0.280 1.360     1     2    2    2  black  0.5        1      0    NA
      3 3 0.0450  0.030 0.650     1     3    3    3  black  0.5        1      0    NA
      
      [[4]][[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      
      [[5]]
      [[5]][[1]]
             x y PANEL group shape colour size fill alpha stroke
      1 0.0665 1     1     1    19   blue    3   NA    NA    0.5
      2 0.5420 2     1     2    19   blue    3   NA    NA    0.5
      3 0.0450 3     1     3    19   blue    3   NA    NA    0.5
      
      [[5]][[2]]
             x   xmin  xmax y PANEL group ymin ymax colour size linetype height alpha
      1 0.0665 -0.778 0.911 1     1     1    1    1  black  0.5        1      0    NA
      2 0.5420 -0.280 1.360 2     1     2    2    2  black  0.5        1      0    NA
      3 0.0450  0.030 0.650 3     1     3    3    3  black  0.5        1      0    NA
      
      [[5]][[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      [[5]][[4]]
             x y
      1 0.0665 1
      2 0.5420 2
      3 0.0450 3
                                                                                      label
      1      list(widehat(italic(beta))=='0.07', italic(t)('5')=='0.16', italic(p)=='0.88')
      2     list(widehat(italic(beta))=='0.54', italic(t)('10')=='1.33', italic(p)=='0.19')
      3 list(widehat(italic(beta))=='0.04', italic(t)('12')=='1.24', italic(p)=='1.00e-03')
        PANEL group    colour  fill size angle alpha family fontface lineheight hjust
      1     1     1 #1B9E77FF white    3     0    NA               1        1.2   0.5
      2     1     2 #D95F02FF white    3     0    NA               1        1.2   0.5
      3     1     3 #7570B3FF white    3     0    NA               1        1.2   0.5
        vjust point.size segment.linetype segment.size segment.curvature
      1   0.5          1                1          0.5                 0
      2   0.5          1                1          0.5                 0
      3   0.5          1                1          0.5                 0
        segment.angle segment.ncp segment.shape segment.square segment.squareShape
      1            90           1           0.5           TRUE                   1
      2            90           1           0.5           TRUE                   1
      3            90           1           0.5           TRUE                   1
        segment.inflect segment.debug
      1           FALSE         FALSE
      2           FALSE         FALSE
      3           FALSE         FALSE
      
      
      [[6]]
      [[6]][[1]]
             x y PANEL group shape colour size fill alpha stroke
      1 0.0665 1     1     1    19   blue    3   NA    NA    0.5
      2 0.5420 2     1     2    19   blue    3   NA    NA    0.5
      3 0.0450 3     1     3    19   blue    3   NA    NA    0.5
      
      [[6]][[2]]
             x   xmin  xmax y PANEL group ymin ymax colour size linetype height alpha
      1 0.0665 -0.778 0.911 1     1     1    1    1  black  0.5        1      0    NA
      2 0.5420 -0.280 1.360 2     1     2    2    2  black  0.5        1      0    NA
      3 0.0450  0.030 0.650 3     1     3    3    3  black  0.5        1      0    NA
      
      [[6]][[3]]
        xintercept PANEL group colour size linetype alpha
      1          0     1    -1  black    1   dashed    NA
      
      [[6]][[4]]
             x y
      1 0.0665 1
      2 0.5420 2
      3 0.0450 3
                                                                                      label
      1      list(widehat(italic(beta))=='0.07', italic(t)('5')=='0.16', italic(p)=='0.88')
      2     list(widehat(italic(beta))=='0.54', italic(t)('10')=='1.33', italic(p)=='0.19')
      3 list(widehat(italic(beta))=='0.04', italic(t)('12')=='1.24', italic(p)=='1.00e-03')
        PANEL group    colour  fill size angle alpha family fontface lineheight hjust
      1     1     1 #1B9E77FF white    3     0    NA               1        1.2   0.5
      2     1     2 #D95F02FF white    3     0    NA               1        1.2   0.5
      3     1     3 #7570B3FF white    3     0    NA               1        1.2   0.5
        vjust point.size segment.linetype segment.size segment.curvature
      1   0.5          1                1          0.5                 0
      2   0.5          1                1          0.5                 0
      3   0.5          1                1          0.5                 0
        segment.angle segment.ncp segment.shape segment.square segment.squareShape
      1            90           1           0.5           TRUE                   1
      2            90           1           0.5           TRUE                   1
      3            90           1           0.5           TRUE                   1
        segment.inflect segment.debug
      1           FALSE         FALSE
      2           FALSE         FALSE
      3           FALSE         FALSE
      
      

---

    Code
      list(pb1$plot$labels, pb2$plot$labels, pb3$plot$labels, pb4$plot$labels)
    Output
      [[1]]
      [[1]]$x
      [1] "estimate"
      
      [[1]]$y
      [1] "term"
      
      [[1]]$title
      NULL
      
      [[1]]$subtitle
      NULL
      
      [[1]]$caption
      NULL
      
      [[1]]$xmin
      [1] "conf.low"
      
      [[1]]$xmax
      [1] "conf.high"
      
      [[1]]$xintercept
      [1] "xintercept"
      
      [[1]]$label
      [1] "label"
      
      [[1]]$alt
      [1] ""
      
      
      [[2]]
      [[2]]$x
      [1] "estimate"
      
      [[2]]$y
      [1] "term"
      
      [[2]]$title
      NULL
      
      [[2]]$subtitle
      NULL
      
      [[2]]$caption
      NULL
      
      [[2]]$xmin
      [1] "conf.low"
      
      [[2]]$xmax
      [1] "conf.high"
      
      [[2]]$xintercept
      [1] "xintercept"
      
      [[2]]$label
      [1] "label"
      
      [[2]]$alt
      [1] ""
      
      
      [[3]]
      [[3]]$x
      [1] "estimate"
      
      [[3]]$y
      [1] "term"
      
      [[3]]$title
      NULL
      
      [[3]]$subtitle
      NULL
      
      [[3]]$caption
      NULL
      
      [[3]]$xmin
      [1] "conf.low"
      
      [[3]]$xmax
      [1] "conf.high"
      
      [[3]]$xintercept
      [1] "xintercept"
      
      [[3]]$alt
      [1] ""
      
      
      [[4]]
      [[4]]$x
      [1] "location"
      
      [[4]]$y
      NULL
      
      [[4]]$title
      NULL
      
      [[4]]$subtitle
      NULL
      
      [[4]]$caption
      NULL
      
      [[4]]$xmin
      [1] "conf.low"
      
      [[4]]$xmax
      [1] "conf.high"
      
      [[4]]$xintercept
      [1] "xintercept"
      
      [[4]]$alt
      [1] ""
      
      

# checking gghistostats plot and parametric stats - data with NAs

    Code
      pb$data
    Output
      [[1]]
           fill  y count   x xmin xmax     density  ncount ndensity flipped_aes PANEL
      1  orange  1     1  60   50   70 0.000617284 0.03125  0.03125       FALSE     1
      2  orange  2     2  80   70   90 0.001234568 0.06250  0.06250       FALSE     1
      3  orange  4     4 100   90  110 0.002469136 0.12500  0.12500       FALSE     1
      4  orange  2     2 120  110  130 0.001234568 0.06250  0.06250       FALSE     1
      5  orange  3     3 140  130  150 0.001851852 0.09375  0.09375       FALSE     1
      6  orange 15    15 160  150  170 0.009259259 0.46875  0.46875       FALSE     1
      7  orange 32    32 180  170  190 0.019753086 1.00000  1.00000       FALSE     1
      8  orange 15    15 200  190  210 0.009259259 0.46875  0.46875       FALSE     1
      9  orange  5     5 220  210  230 0.003086420 0.15625  0.15625       FALSE     1
      10 orange  1     1 240  230  250 0.000617284 0.03125  0.03125       FALSE     1
      11 orange  1     1 260  250  270 0.000617284 0.03125  0.03125       FALSE     1
         group ymin ymax colour size linetype alpha
      1     -1    0    1  black  0.5        1   0.7
      2     -1    0    2  black  0.5        1   0.7
      3     -1    0    4  black  0.5        1   0.7
      4     -1    0    2  black  0.5        1   0.7
      5     -1    0    3  black  0.5        1   0.7
      6     -1    0   15  black  0.5        1   0.7
      7     -1    0   32  black  0.5        1   0.7
      8     -1    0   15  black  0.5        1   0.7
      9     -1    0    5  black  0.5        1   0.7
      10    -1    0    1  black  0.5        1   0.7
      11    -1    0    1  black  0.5        1   0.7
      
      [[2]]
        xintercept PANEL group colour size linetype alpha
      1    174.358     1    -1   blue    1   dashed    NA
      

---

    Code
      within(pb$plot$labels, rm(subtitle, caption))
    Output
      $x
      [1] "character height"
      
      $y
      [1] "count"
      
      $title
      [1] "starwars: character heights"
      
      $fill
      [1] "count"
      
      $xintercept
      [1] "xintercept"
      
      $weight
      [1] "weight"
      attr(,"fallback")
      [1] TRUE
      
      $alt
      [1] ""
      

# checking gghistostats and non-parametric stats - data without NAs

    Code
      pb$data
    Output
      [[1]]
          fill  y count  x xmin xmax     density     ncount   ndensity flipped_aes
      1 grey50 33    33 10  7.5 12.5 0.028205128 0.33333333 0.33333333       FALSE
      2 grey50 99    99 15 12.5 17.5 0.084615385 1.00000000 1.00000000       FALSE
      3 grey50 84    84 20 17.5 22.5 0.071794872 0.84848485 0.84848485       FALSE
      4 grey50 13    13 25 22.5 27.5 0.011111111 0.13131313 0.13131313       FALSE
      5 grey50  3     3 30 27.5 32.5 0.002564103 0.03030303 0.03030303       FALSE
      6 grey50  2     2 35 32.5 37.5 0.001709402 0.02020202 0.02020202       FALSE
        PANEL group ymin ymax colour size linetype alpha
      1     1    -1    0   33  black  0.5        1   0.7
      2     1    -1    0   99  black  0.5        1   0.7
      3     1    -1    0   84  black  0.5        1   0.7
      4     1    -1    0   13  black  0.5        1   0.7
      5     1    -1    0    3  black  0.5        1   0.7
      6     1    -1    0    2  black  0.5        1   0.7
      
      [[2]]
        xintercept PANEL group colour size linetype alpha
      1         17     1    -1   blue    1   dashed    NA
      

---

    Code
      pb$layout$panel_params[[1]]$y.sec$break_info
    Output
      $range
      [1] -0.02115385  0.44423077
      
      $labels
      [1] "0%"  "10%" "20%" "30%" "40%"
      
      $major
      [1] 0.045 0.260 0.475 0.690 0.905
      
      $minor
      [1] 0.045 0.153 0.260 0.367 0.475 0.583 0.690 0.798 0.905
      
      $major_source
      [1] -0.04459459 23.39234234 46.82927928 70.15720721 93.59414414
      
      $minor_source
      [1] -0.04459459 11.72837838 23.39234234 35.05630631 46.82927928 58.49324324
      [7] 70.15720721 81.93018018 93.59414414
      
      $major_source_user
      [1] 0.0 0.1 0.2 0.3 0.4
      
      $minor_source_user
      [1] 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40
      

---

    Code
      pb$plot$labels
    Output
      $x
      [1] "city miles per gallon"
      
      $y
      [1] "count"
      
      $title
      [1] "fuel economy"
      
      $subtitle
      NULL
      
      $caption
      [1] "source: government website"
      
      $fill
      [1] "count"
      
      $xintercept
      [1] "xintercept"
      
      $weight
      [1] "weight"
      attr(,"fallback")
      [1] TRUE
      
      $alt
      [1] ""
      

# checking robust stats and proportions

    Code
      pb$data
    Output
      [[1]]
           fill y count   x xmin xmax density    ncount  ndensity flipped_aes PANEL
      1  grey50 0     0 1.0   NA 1.25  0.0000 0.0000000 0.0000000       FALSE     1
      2  grey50 2     2 1.5 1.25 1.75  0.1250 0.2222222 0.2222222       FALSE     1
      3  grey50 4     4 2.0 1.75 2.25  0.2500 0.4444444 0.4444444       FALSE     1
      4  grey50 3     3 2.5 2.25 2.75  0.1875 0.3333333 0.3333333       FALSE     1
      5  grey50 7     7 3.0 2.75 3.25  0.4375 0.7777778 0.7777778       FALSE     1
      6  grey50 9     9 3.5 3.25 3.75  0.5625 1.0000000 1.0000000       FALSE     1
      7  grey50 4     4 4.0 3.75 4.25  0.2500 0.4444444 0.4444444       FALSE     1
      8  grey50 0     0 4.5 4.25 4.75  0.0000 0.0000000 0.0000000       FALSE     1
      9  grey50 1     1 5.0 4.75 5.25  0.0625 0.1111111 0.1111111       FALSE     1
      10 grey50 2     2 5.5 5.25 5.75  0.1250 0.2222222 0.2222222       FALSE     1
      11 grey50 0     0 6.0 5.75   NA  0.0000 0.0000000 0.0000000       FALSE     1
         group ymin ymax colour size linetype alpha
      1     -1    0    0  black  0.5        1   0.7
      2     -1    0    2  black  0.5        1   0.7
      3     -1    0    4  black  0.5        1   0.7
      4     -1    0    3  black  0.5        1   0.7
      5     -1    0    7  black  0.5        1   0.7
      6     -1    0    9  black  0.5        1   0.7
      7     -1    0    4  black  0.5        1   0.7
      8     -1    0    0  black  0.5        1   0.7
      9     -1    0    1  black  0.5        1   0.7
      10    -1    0    2  black  0.5        1   0.7
      11    -1    0    0  black  0.5        1   0.7
      
      [[2]]
        xintercept PANEL group colour size linetype alpha
      1      3.197     1    -1   blue    1   dashed    NA
      

---

    Code
      within(pb$plot$labels, rm(subtitle))
    Output
      $x
      [1] "wt"
      
      $y
      [1] "count"
      
      $title
      NULL
      
      $caption
      NULL
      
      $fill
      [1] "count"
      
      $xintercept
      [1] "xintercept"
      
      $weight
      [1] "weight"
      attr(,"fallback")
      [1] TRUE
      
      $alt
      [1] ""
      

# checking if normal curve work

    Code
      pb1$data
    Output
      [[1]]
           fill  y count  x xmin xmax    density     ncount   ndensity flipped_aes
      1  grey50  2     2  4  3.5  4.5 0.02409639 0.16666667 0.16666667       FALSE
      2  grey50  1     1  5  4.5  5.5 0.01204819 0.08333333 0.08333333       FALSE
      3  grey50  2     2  6  5.5  6.5 0.02409639 0.16666667 0.16666667       FALSE
      4  grey50  3     3  7  6.5  7.5 0.03614458 0.25000000 0.25000000       FALSE
      5  grey50  4     4  8  7.5  8.5 0.04819277 0.33333333 0.33333333       FALSE
      6  grey50  4     4  9  8.5  9.5 0.04819277 0.33333333 0.33333333       FALSE
      7  grey50  7     7 10  9.5 10.5 0.08433735 0.58333333 0.58333333       FALSE
      8  grey50  7     7 11 10.5 11.5 0.08433735 0.58333333 0.58333333       FALSE
      9  grey50  2     2 12 11.5 12.5 0.02409639 0.16666667 0.16666667       FALSE
      10 grey50  6     6 13 12.5 13.5 0.07228916 0.50000000 0.50000000       FALSE
      11 grey50 12    12 14 13.5 14.5 0.14457831 1.00000000 1.00000000       FALSE
      12 grey50  8     8 15 14.5 15.5 0.09638554 0.66666667 0.66666667       FALSE
      13 grey50  5     5 16 15.5 16.5 0.06024096 0.41666667 0.41666667       FALSE
      14 grey50  1     1 17 16.5 17.5 0.01204819 0.08333333 0.08333333       FALSE
      15 grey50  4     4 18 17.5 18.5 0.04819277 0.33333333 0.33333333       FALSE
      16 grey50  4     4 19 18.5 19.5 0.04819277 0.33333333 0.33333333       FALSE
      17 grey50  5     5 20 19.5 20.5 0.06024096 0.41666667 0.41666667       FALSE
      18 grey50  5     5 21 20.5 21.5 0.06024096 0.41666667 0.41666667       FALSE
      19 grey50  1     1 22 21.5 22.5 0.01204819 0.08333333 0.08333333       FALSE
         PANEL group ymin ymax colour size linetype alpha
      1      1    -1    0    2  black  0.5        1   0.7
      2      1    -1    0    1  black  0.5        1   0.7
      3      1    -1    0    2  black  0.5        1   0.7
      4      1    -1    0    3  black  0.5        1   0.7
      5      1    -1    0    4  black  0.5        1   0.7
      6      1    -1    0    4  black  0.5        1   0.7
      7      1    -1    0    7  black  0.5        1   0.7
      8      1    -1    0    7  black  0.5        1   0.7
      9      1    -1    0    2  black  0.5        1   0.7
      10     1    -1    0    6  black  0.5        1   0.7
      11     1    -1    0   12  black  0.5        1   0.7
      12     1    -1    0    8  black  0.5        1   0.7
      13     1    -1    0    5  black  0.5        1   0.7
      14     1    -1    0    1  black  0.5        1   0.7
      15     1    -1    0    4  black  0.5        1   0.7
      16     1    -1    0    4  black  0.5        1   0.7
      17     1    -1    0    5  black  0.5        1   0.7
      18     1    -1    0    5  black  0.5        1   0.7
      19     1    -1    0    1  black  0.5        1   0.7
      
      [[2]]
              x         y PANEL group colour size linetype alpha
      1    4.10 0.7752851     1    -1    red  0.8        1    NA
      2    4.28 0.8442004     1    -1    red  0.8        1    NA
      3    4.46 0.9177403     1    -1    red  0.8        1    NA
      4    4.64 0.9960569     1    -1    red  0.8        1    NA
      5    4.82 1.0792910     1    -1    red  0.8        1    NA
      6    5.00 1.1675703     1    -1    red  0.8        1    NA
      7    5.18 1.2610075     1    -1    red  0.8        1    NA
      8    5.36 1.3596977     1    -1    red  0.8        1    NA
      9    5.54 1.4637170     1    -1    red  0.8        1    NA
      10   5.72 1.5731206     1    -1    red  0.8        1    NA
      11   5.90 1.6879399     1    -1    red  0.8        1    NA
      12   6.08 1.8081815     1    -1    red  0.8        1    NA
      13   6.26 1.9338250     1    -1    red  0.8        1    NA
      14   6.44 2.0648210     1    -1    red  0.8        1    NA
      15   6.62 2.2010897     1    -1    red  0.8        1    NA
      16   6.80 2.3425193     1    -1    red  0.8        1    NA
      17   6.98 2.4889644     1    -1    red  0.8        1    NA
      18   7.16 2.6402454     1    -1    red  0.8        1    NA
      19   7.34 2.7961469     1    -1    red  0.8        1    NA
      20   7.52 2.9564176     1    -1    red  0.8        1    NA
      21   7.70 3.1207693     1    -1    red  0.8        1    NA
      22   7.88 3.2888770     1    -1    red  0.8        1    NA
      23   8.06 3.4603792     1    -1    red  0.8        1    NA
      24   8.24 3.6348781     1    -1    red  0.8        1    NA
      25   8.42 3.8119404     1    -1    red  0.8        1    NA
      26   8.60 3.9910984     1    -1    red  0.8        1    NA
      27   8.78 4.1718518     1    -1    red  0.8        1    NA
      28   8.96 4.3536688     1    -1    red  0.8        1    NA
      29   9.14 4.5359891     1    -1    red  0.8        1    NA
      30   9.32 4.7182256     1    -1    red  0.8        1    NA
      31   9.50 4.8997678     1    -1    red  0.8        1    NA
      32   9.68 5.0799845     1    -1    red  0.8        1    NA
      33   9.86 5.2582274     1    -1    red  0.8        1    NA
      34  10.04 5.4338349     1    -1    red  0.8        1    NA
      35  10.22 5.6061356     1    -1    red  0.8        1    NA
      36  10.40 5.7744530     1    -1    red  0.8        1    NA
      37  10.58 5.9381094     1    -1    red  0.8        1    NA
      38  10.76 6.0964306     1    -1    red  0.8        1    NA
      39  10.94 6.2487502     1    -1    red  0.8        1    NA
      40  11.12 6.3944144     1    -1    red  0.8        1    NA
      41  11.30 6.5327869     1    -1    red  0.8        1    NA
      42  11.48 6.6632528     1    -1    red  0.8        1    NA
      43  11.66 6.7852239     1    -1    red  0.8        1    NA
      44  11.84 6.8981427     1    -1    red  0.8        1    NA
      45  12.02 7.0014864     1    -1    red  0.8        1    NA
      46  12.20 7.0947716     1    -1    red  0.8        1    NA
      47  12.38 7.1775574     1    -1    red  0.8        1    NA
      48  12.56 7.2494495     1    -1    red  0.8        1    NA
      49  12.74 7.3101025     1    -1    red  0.8        1    NA
      50  12.92 7.3592237     1    -1    red  0.8        1    NA
      51  13.10 7.3965744     1    -1    red  0.8        1    NA
      52  13.28 7.4219726     1    -1    red  0.8        1    NA
      53  13.46 7.4352942     1    -1    red  0.8        1    NA
      54  13.64 7.4364738     1    -1    red  0.8        1    NA
      55  13.82 7.4255059     1    -1    red  0.8        1    NA
      56  14.00 7.4024440     1    -1    red  0.8        1    NA
      57  14.18 7.3674009     1    -1    red  0.8        1    NA
      58  14.36 7.3205476     1    -1    red  0.8        1    NA
      59  14.54 7.2621118     1    -1    red  0.8        1    NA
      60  14.72 7.1923759     1    -1    red  0.8        1    NA
      61  14.90 7.1116753     1    -1    red  0.8        1    NA
      62  15.08 7.0203950     1    -1    red  0.8        1    NA
      63  15.26 6.9189672     1    -1    red  0.8        1    NA
      64  15.44 6.8078674     1    -1    red  0.8        1    NA
      65  15.62 6.6876109     1    -1    red  0.8        1    NA
      66  15.80 6.5587487     1    -1    red  0.8        1    NA
      67  15.98 6.4218637     1    -1    red  0.8        1    NA
      68  16.16 6.2775656     1    -1    red  0.8        1    NA
      69  16.34 6.1264873     1    -1    red  0.8        1    NA
      70  16.52 5.9692793     1    -1    red  0.8        1    NA
      71  16.70 5.8066060     1    -1    red  0.8        1    NA
      72  16.88 5.6391403     1    -1    red  0.8        1    NA
      73  17.06 5.4675598     1    -1    red  0.8        1    NA
      74  17.24 5.2925415     1    -1    red  0.8        1    NA
      75  17.42 5.1147580     1    -1    red  0.8        1    NA
      76  17.60 4.9348733     1    -1    red  0.8        1    NA
      77  17.78 4.7535384     1    -1    red  0.8        1    NA
      78  17.96 4.5713882     1    -1    red  0.8        1    NA
      79  18.14 4.3890375     1    -1    red  0.8        1    NA
      80  18.32 4.2070781     1    -1    red  0.8        1    NA
      81  18.50 4.0260758     1    -1    red  0.8        1    NA
      82  18.68 3.8465679     1    -1    red  0.8        1    NA
      83  18.86 3.6690612     1    -1    red  0.8        1    NA
      84  19.04 3.4940298     1    -1    red  0.8        1    NA
      85  19.22 3.3219137     1    -1    red  0.8        1    NA
      86  19.40 3.1531176     1    -1    red  0.8        1    NA
      87  19.58 2.9880103     1    -1    red  0.8        1    NA
      88  19.76 2.8269238     1    -1    red  0.8        1    NA
      89  19.94 2.6701533     1    -1    red  0.8        1    NA
      90  20.12 2.5179575     1    -1    red  0.8        1    NA
      91  20.30 2.3705585     1    -1    red  0.8        1    NA
      92  20.48 2.2281430     1    -1    red  0.8        1    NA
      93  20.66 2.0908627     1    -1    red  0.8        1    NA
      94  20.84 1.9588360     1    -1    red  0.8        1    NA
      95  21.02 1.8321488     1    -1    red  0.8        1    NA
      96  21.20 1.7108561     1    -1    red  0.8        1    NA
      97  21.38 1.5949839     1    -1    red  0.8        1    NA
      98  21.56 1.4845309     1    -1    red  0.8        1    NA
      99  21.74 1.3794700     1    -1    red  0.8        1    NA
      100 21.92 1.2797507     1    -1    red  0.8        1    NA
      101 22.10 1.1853008     1    -1    red  0.8        1    NA
      
      [[3]]
        xintercept PANEL group colour size linetype alpha
      1   13.56747     1    -1   blue    1   dashed    NA
      

---

    Code
      pb1$plot$labels
    Output
      $x
      [1] "awake"
      
      $y
      [1] "count"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      NULL
      
      $fill
      [1] "count"
      
      $xintercept
      [1] "xintercept"
      
      $weight
      [1] "weight"
      attr(,"fallback")
      [1] TRUE
      
      $alt
      [1] ""
      

# `pairwise_comparisons()` - test additional arguments

    Code
      list(df1, df2, df3, df4)
    Output
      [[1]]
      # A tibble: 6 x 6
        group1 group2  p.value test.details     p.value.adjustment
        <chr>  <chr>     <dbl> <chr>            <chr>             
      1 HDHF   HDLF   2.65e- 4 Student's t-test None              
      2 HDHF   LDHF   3.51e- 2 Student's t-test None              
      3 HDHF   LDLF   3.29e-13 Student's t-test None              
      4 HDLF   LDHF   9.72e- 1 Student's t-test None              
      5 HDLF   LDLF   6.62e- 4 Student's t-test None              
      6 LDHF   LDLF   1.11e- 9 Student's t-test None              
        label                                  
        <chr>                                  
      1 list(~italic(p)[uncorrected]==2.65e-04)
      2 list(~italic(p)[uncorrected]==0.04)    
      3 list(~italic(p)[uncorrected]==3.29e-13)
      4 list(~italic(p)[uncorrected]==0.97)    
      5 list(~italic(p)[uncorrected]==6.62e-04)
      6 list(~italic(p)[uncorrected]==1.11e-09)
      
      [[2]]
      # A tibble: 6 x 6
        group1 group2 p.value test.details     p.value.adjustment
        <chr>  <chr>    <dbl> <chr>            <chr>             
      1 HDHF   HDLF    1.00   Student's t-test None              
      2 HDHF   LDHF    0.965  Student's t-test None              
      3 HDHF   LDLF    1.00   Student's t-test None              
      4 HDLF   LDHF    0.0281 Student's t-test None              
      5 HDLF   LDLF    0.999  Student's t-test None              
      6 LDHF   LDLF    1.00   Student's t-test None              
        label                              
        <chr>                              
      1 list(~italic(p)[uncorrected]==1.00)
      2 list(~italic(p)[uncorrected]==0.96)
      3 list(~italic(p)[uncorrected]==1.00)
      4 list(~italic(p)[uncorrected]==0.03)
      5 list(~italic(p)[uncorrected]==1.00)
      6 list(~italic(p)[uncorrected]==1.00)
      
      [[3]]
      # A tibble: 3 x 6
        group1 group2 p.value test.details     p.value.adjustment
        <chr>  <chr>    <dbl> <chr>            <chr>             
      1 4      6        0.995 Student's t-test None              
      2 4      8        1.00  Student's t-test None              
      3 6      8        0.997 Student's t-test None              
        label                              
        <chr>                              
      1 list(~italic(p)[uncorrected]==0.99)
      2 list(~italic(p)[uncorrected]==1.00)
      3 list(~italic(p)[uncorrected]==1.00)
      
      [[4]]
      # A tibble: 3 x 6
        group1 group2     p.value test.details     p.value.adjustment
        <chr>  <chr>        <dbl> <chr>            <chr>             
      1 4      6      0.00532     Student's t-test None              
      2 4      8      0.000000103 Student's t-test None              
      3 6      8      0.00258     Student's t-test None              
        label                                  
        <chr>                                  
      1 list(~italic(p)[uncorrected]==5.32e-03)
      2 list(~italic(p)[uncorrected]==1.03e-07)
      3 list(~italic(p)[uncorrected]==2.58e-03)
      

# checking if extract_stats works

    Code
      length(extract_stats(p1))
    Output
      [1] 5

---

    Code
      length(extract_stats(p2))
    Output
      [1] 5

---

    Code
      length(extract_stats(p3))
    Output
      [1] 5

---

    Code
      length(extract_stats(p4)$pairwise_comparisons_data)
    Output
      [1] 11

---

    Code
      length(extract_stats(p5))
    Output
      [1] 5

---

    Code
      length(extract_stats(p6))
    Output
      [1] 5

# checking labels with counts

    Code
      pb$data
    Output
      [[1]]
             fill         y x PANEL group flipped_aes      ymin      ymax xmin xmax
      1 #1B9E77FF 1.0000000 1     1     1       FALSE 0.9154362 1.0000000 0.55 1.45
      2 #1B9E77FF 1.0000000 2     1     3       FALSE 0.5161744 1.0000000 1.55 2.45
      3 #D95F02FF 0.9154362 1     1     2       FALSE 0.0000000 0.9154362 0.55 1.45
      4 #D95F02FF 0.5161744 2     1     4       FALSE 0.0000000 0.5161744 1.55 2.45
        colour size linetype alpha
      1  black  0.5        1    NA
      2  black  0.5        1    NA
      3  black  0.5        1    NA
      4  black  0.5        1    NA
      
      [[2]]
                y x  label group PANEL      ymax xmin xmax      ymin colour  fill
      1 0.9577181 1  8.46%     1     1 1.0000000    1    1 0.9154362  black white
      2 0.7580872 2 48.38%     1     1 1.0000000    2    2 0.5161744  black white
      3 0.4577181 1 91.54%     2     1 0.9154362    1    1 0.0000000  black white
      4 0.2580872 2 51.62%     2     1 0.5161744    2    2 0.0000000  black white
        size angle hjust vjust alpha family fontface lineheight
      1 3.88     0   0.5   0.5     1               1        1.2
      2 3.88     0   0.5   0.5     1               1        1.2
      3 3.88     0   0.5   0.5     1               1        1.2
      4 3.88     0   0.5   0.5     1               1        1.2
      
      [[3]]
           y x                         label PANEL group colour size angle hjust
      1 1.05 2      list(~italic(p)=='0.39')     1     2  black  2.8     0   0.5
      2 1.05 1 list(~italic(p)=='1.08e-225')     1     1  black  2.8     0   0.5
        vjust alpha family fontface lineheight
      1   0.5    NA               1        1.2
      2   0.5    NA               1        1.2
      
      [[4]]
            y x       label PANEL group colour size angle hjust vjust alpha family
      1 -0.05 2   (n = 711)     1     2  black    4     0   0.5   0.5    NA       
      2 -0.05 1 (n = 1,490)     1     1  black    4     0   0.5   0.5    NA       
        fontface lineheight
      1        1        1.2
      2        1        1.2
      

---

    Code
      within(pb$plot$labels, rm(subtitle))
    Output
      $x
      [1] "Passenger sex"
      
      $y
      [1] "proportion"
      
      $title
      NULL
      
      $caption
      NULL
      
      $fill
      [1] "Sex"
      
      $label
      [1] ".label"
      
      $group
      [1] "Sex"
      
      $alt
      [1] ""
      

# aesthetic modifications

    Code
      list(pb$data, pb1$data)
    Output
      [[1]]
      [[1]][[1]]
             fill          y x PANEL group flipped_aes       ymin       ymax xmin
      1 #9A8822FF 1.00000000 1     1     1       FALSE 0.09090909 1.00000000 0.55
      2 #9A8822FF 1.00000000 2     1     3       FALSE 0.42857143 1.00000000 1.55
      3 #F5CDB4FF 0.09090909 1     1     2       FALSE 0.00000000 0.09090909 0.55
      4 #F5CDB4FF 0.42857143 2     1     4       FALSE 0.00000000 0.42857143 1.55
      5 #F5CDB4FF 1.00000000 3     1     5       FALSE 0.00000000 1.00000000 2.55
        xmax colour size linetype alpha
      1 1.45  black  0.5        1    NA
      2 2.45  black  0.5        1    NA
      3 1.45  black  0.5        1    NA
      4 2.45  black  0.5        1    NA
      5 3.45  black  0.5        1    NA
      
      [[1]][[2]]
                 y x      label group PANEL       ymax xmin xmax       ymin colour
      1 0.54545455 1  10\n(91%)     1     1 1.00000000    1    1 0.09090909  black
      2 0.71428571 2   4\n(57%)     1     1 1.00000000    2    2 0.42857143  black
      3 0.04545455 1    1\n(9%)     2     1 0.09090909    1    1 0.00000000  black
      4 0.21428571 2   3\n(43%)     2     1 0.42857143    2    2 0.00000000  black
      5 0.50000000 3 14\n(100%)     2     1 1.00000000    3    3 0.00000000  black
         fill size angle hjust vjust alpha family fontface lineheight
      1 white 3.88     0   0.5   0.5     1               1        1.2
      2 white 3.88     0   0.5   0.5     1               1        1.2
      3 white 3.88     0   0.5   0.5     1               1        1.2
      4 white 3.88     0   0.5   0.5     1               1        1.2
      5 white 3.88     0   0.5   0.5     1               1        1.2
      
      [[1]][[3]]
           y x                        label PANEL group colour size angle hjust vjust
      1 1.05 3 list(~italic(p)=='1.83e-04')     1     3  black  2.8     0   0.5   0.5
      2 1.05 2     list(~italic(p)=='0.71')     1     2  black  2.8     0   0.5   0.5
      3 1.05 1 list(~italic(p)=='6.66e-03')     1     1  black  2.8     0   0.5   0.5
        alpha family fontface lineheight
      1    NA               1        1.2
      2    NA               1        1.2
      3    NA               1        1.2
      
      [[1]][[4]]
            y x    label PANEL group colour size angle hjust vjust alpha family
      1 -0.05 3 (n = 14)     1     3  black    4     0   0.5   0.5    NA       
      2 -0.05 2  (n = 7)     1     2  black    4     0   0.5   0.5    NA       
      3 -0.05 1 (n = 11)     1     1  black    4     0   0.5   0.5    NA       
        fontface lineheight
      1        1        1.2
      2        1        1.2
      3        1        1.2
      
      
      [[2]]
      [[2]][[1]]
             fill          y x PANEL group flipped_aes       ymin       ymax xmin
      1 #1B9E77FF 1.00000000 1     1     1       FALSE 0.09090909 1.00000000 0.55
      2 #1B9E77FF 1.00000000 2     1     3       FALSE 0.42857143 1.00000000 1.55
      3 #D95F02FF 0.09090909 1     1     2       FALSE 0.00000000 0.09090909 0.55
      4 #D95F02FF 0.42857143 2     1     4       FALSE 0.00000000 0.42857143 1.55
      5 #D95F02FF 1.00000000 3     1     5       FALSE 0.00000000 1.00000000 2.55
        xmax colour size linetype alpha
      1 1.45  black  0.5        1    NA
      2 2.45  black  0.5        1    NA
      3 1.45  black  0.5        1    NA
      4 2.45  black  0.5        1    NA
      5 3.45  black  0.5        1    NA
      
      [[2]][[2]]
                 y x label group PANEL       ymax xmin xmax       ymin colour  fill
      1 0.54545455 1    10     1     1 1.00000000    1    1 0.09090909  black white
      2 0.71428571 2     4     1     1 1.00000000    2    2 0.42857143  black white
      3 0.04545455 1     1     2     1 0.09090909    1    1 0.00000000  black white
      4 0.21428571 2     3     2     1 0.42857143    2    2 0.00000000  black white
      5 0.50000000 3    14     2     1 1.00000000    3    3 0.00000000  black white
        size angle hjust vjust alpha family fontface lineheight
      1 3.88     0   0.5   0.5     1               1        1.2
      2 3.88     0   0.5   0.5     1               1        1.2
      3 3.88     0   0.5   0.5     1               1        1.2
      4 3.88     0   0.5   0.5     1               1        1.2
      5 3.88     0   0.5   0.5     1               1        1.2
      
      [[2]][[3]]
           y x                        label PANEL group colour size angle hjust vjust
      1 1.05 3 list(~italic(p)=='1.83e-04')     1     3  black  2.8     0   0.5   0.5
      2 1.05 2     list(~italic(p)=='0.71')     1     2  black  2.8     0   0.5   0.5
      3 1.05 1 list(~italic(p)=='6.66e-03')     1     1  black  2.8     0   0.5   0.5
        alpha family fontface lineheight
      1    NA               1        1.2
      2    NA               1        1.2
      3    NA               1        1.2
      
      [[2]][[4]]
            y x    label PANEL group colour size angle hjust vjust alpha family
      1 -0.05 3 (n = 14)     1     3  black    4     0   0.5   0.5    NA       
      2 -0.05 2  (n = 7)     1     2  black    4     0   0.5   0.5    NA       
      3 -0.05 1 (n = 11)     1     1  black    4     0   0.5   0.5    NA       
        fontface lineheight
      1        1        1.2
      2        1        1.2
      3        1        1.2
      
      

# dropped factor levels

    Code
      pb$data
    Output
      [[1]]
             fill         y x PANEL group flipped_aes      ymin      ymax xmin xmax
      1 #1B9E77FF 1.0000000 1     1     1       FALSE 0.3684211 1.0000000 0.55 1.45
      2 #D95F02FF 0.3684211 1     1     2       FALSE 0.1578947 0.3684211 0.55 1.45
      3 #7570B3FF 0.1578947 1     1     3       FALSE 0.0000000 0.1578947 0.55 1.45
        colour size linetype alpha
      1  black  0.5        1    NA
      2  black  0.5        1    NA
      3  black  0.5        1    NA
      
      [[2]]
                 y x label group PANEL      ymax xmin xmax      ymin colour  fill
      1 0.68421053 1   63%     1     1 1.0000000    1    1 0.3684211  black white
      2 0.26315789 1   21%     2     1 0.3684211    1    1 0.1578947  black white
      3 0.07894737 1   16%     3     1 0.1578947    1    1 0.0000000  black white
        size angle hjust vjust alpha family fontface lineheight
      1 3.88     0   0.5   0.5     1               1        1.2
      2 3.88     0   0.5   0.5     1               1        1.2
      3 3.88     0   0.5   0.5     1               1        1.2
      
      [[3]]
            y x    label PANEL group colour size angle hjust vjust alpha family
      1 -0.05 1 (n = 19)     1     1  black    4     0   0.5   0.5    NA       
        fontface lineheight
      1        1        1.2
      

---

    Code
      pb$plot$labels
    Output
      $x
      [1] "am"
      
      $y
      NULL
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      NULL
      
      $fill
      [1] "cyl"
      
      $label
      [1] ".label"
      
      $group
      [1] "cyl"
      
      $alt
      [1] ""
      

# basic plotting works - two groups

    Code
      list(pb1$data[[1]], pb1$data[[2]], pb1$data[[4]], pb1$data[[5]], pb1$data[[6]],
      pb1$data[[7]])
    Output
      [[1]]
            colour x    y group PANEL shape size fill alpha stroke
      1  #D95F02FF 2  9.0     1     1    19    3   NA   0.5    0.5
      2  #D95F02FF 2 10.0     2     1    19    3   NA   0.5    0.5
      3  #D95F02FF 2 10.0     3     1    19    3   NA   0.5    0.5
      4  #D95F02FF 2  6.0     4     1    19    3   NA   0.5    0.5
      5  #D95F02FF 2  5.5     5     1    19    3   NA   0.5    0.5
      6  #D95F02FF 2  7.5     6     1    19    3   NA   0.5    0.5
      7  #D95F02FF 2 10.0     7     1    19    3   NA   0.5    0.5
      8  #D95F02FF 2  9.0     8     1    19    3   NA   0.5    0.5
      9  #D95F02FF 2  6.0     9     1    19    3   NA   0.5    0.5
      10 #D95F02FF 2  0.0    10     1    19    3   NA   0.5    0.5
      11 #D95F02FF 2  8.5    11     1    19    3   NA   0.5    0.5
      12 #D95F02FF 2  6.5    12     1    19    3   NA   0.5    0.5
      13 #D95F02FF 2  4.0    13     1    19    3   NA   0.5    0.5
      14 #D95F02FF 2  6.0    14     1    19    3   NA   0.5    0.5
      15 #D95F02FF 2  8.5    15     1    19    3   NA   0.5    0.5
      16 #D95F02FF 2 10.0    16     1    19    3   NA   0.5    0.5
      17 #D95F02FF 2  7.5    17     1    19    3   NA   0.5    0.5
      18 #D95F02FF 2 10.0    18     1    19    3   NA   0.5    0.5
      19 #D95F02FF 2  8.5    19     1    19    3   NA   0.5    0.5
      20 #D95F02FF 2  5.0    20     1    19    3   NA   0.5    0.5
      21 #D95F02FF 2  4.5    21     1    19    3   NA   0.5    0.5
      22 #D95F02FF 2  9.0    22     1    19    3   NA   0.5    0.5
      23 #D95F02FF 2  4.0    23     1    19    3   NA   0.5    0.5
      24 #D95F02FF 2  4.5    24     1    19    3   NA   0.5    0.5
      25 #D95F02FF 2  3.5    25     1    19    3   NA   0.5    0.5
      26 #D95F02FF 2 10.0    26     1    19    3   NA   0.5    0.5
      27 #D95F02FF 2  8.0    27     1    19    3   NA   0.5    0.5
      28 #1B9E77FF 1 10.0     1     1    19    3   NA   0.5    0.5
      29 #1B9E77FF 1 10.0     2     1    19    3   NA   0.5    0.5
      30 #1B9E77FF 1 10.0     3     1    19    3   NA   0.5    0.5
      31 #1B9E77FF 1  9.0     4     1    19    3   NA   0.5    0.5
      32 #1B9E77FF 1  8.5     5     1    19    3   NA   0.5    0.5
      33 #1B9E77FF 1  3.0     6     1    19    3   NA   0.5    0.5
      34 #1B9E77FF 1 10.0     7     1    19    3   NA   0.5    0.5
      35 #1B9E77FF 1 10.0     8     1    19    3   NA   0.5    0.5
      36 #1B9E77FF 1 10.0     9     1    19    3   NA   0.5    0.5
      37 #1B9E77FF 1  0.0    10     1    19    3   NA   0.5    0.5
      38 #1B9E77FF 1 10.0    11     1    19    3   NA   0.5    0.5
      39 #1B9E77FF 1  8.5    12     1    19    3   NA   0.5    0.5
      40 #1B9E77FF 1  8.5    13     1    19    3   NA   0.5    0.5
      41 #1B9E77FF 1  7.0    14     1    19    3   NA   0.5    0.5
      42 #1B9E77FF 1  9.0    15     1    19    3   NA   0.5    0.5
      43 #1B9E77FF 1 10.0    16     1    19    3   NA   0.5    0.5
      44 #1B9E77FF 1 10.0    17     1    19    3   NA   0.5    0.5
      45 #1B9E77FF 1 10.0    18     1    19    3   NA   0.5    0.5
      46 #1B9E77FF 1  9.5    19     1    19    3   NA   0.5    0.5
      47 #1B9E77FF 1  9.5    20     1    19    3   NA   0.5    0.5
      48 #1B9E77FF 1  7.5    21     1    19    3   NA   0.5    0.5
      49 #1B9E77FF 1 10.0    22     1    19    3   NA   0.5    0.5
      50 #1B9E77FF 1  1.5    23     1    19    3   NA   0.5    0.5
      51 #1B9E77FF 1  5.5    24     1    19    3   NA   0.5    0.5
      52 #1B9E77FF 1 10.0    25     1    19    3   NA   0.5    0.5
      53 #1B9E77FF 1 10.0    26     1    19    3   NA   0.5    0.5
      54 #1B9E77FF 1  8.0    27     1    19    3   NA   0.5    0.5
      
      [[2]]
        ymin lower middle upper ymax           outliers notchupper notchlower x
      1    7  8.25    9.5    10   10 3.0, 0.0, 1.5, 5.5  10.032124   8.967876 1
      2    0  5.25    7.5     9   10                      8.640267   6.359733 2
        flipped_aes PANEL group ymin_final ymax_final xmin xmax xid newx new_width
      1       FALSE     1     1          0         10  0.9  1.1   1    1       0.2
      2       FALSE     1     2          0         10  1.9  2.1   2    2       0.2
        weight colour  fill size alpha shape linetype
      1      1 grey20 white  0.5   0.5    19    solid
      2      1 grey20 white  0.5   0.5    19    solid
      
      [[3]]
         x    y group PANEL colour size linetype alpha
      1  2  9.0     1     1    red  0.5        1    NA
      2  2 10.0     2     1    red  0.5        1    NA
      3  2 10.0     3     1    red  0.5        1    NA
      4  2  6.0     4     1    red  0.5        1    NA
      5  2  5.5     5     1    red  0.5        1    NA
      6  2  7.5     6     1    red  0.5        1    NA
      7  2 10.0     7     1    red  0.5        1    NA
      8  2  9.0     8     1    red  0.5        1    NA
      9  2  6.0     9     1    red  0.5        1    NA
      10 2  0.0    10     1    red  0.5        1    NA
      11 2  8.5    11     1    red  0.5        1    NA
      12 2  6.5    12     1    red  0.5        1    NA
      13 2  4.0    13     1    red  0.5        1    NA
      14 2  6.0    14     1    red  0.5        1    NA
      15 2  8.5    15     1    red  0.5        1    NA
      16 2 10.0    16     1    red  0.5        1    NA
      17 2  7.5    17     1    red  0.5        1    NA
      18 2 10.0    18     1    red  0.5        1    NA
      19 2  8.5    19     1    red  0.5        1    NA
      20 2  5.0    20     1    red  0.5        1    NA
      21 2  4.5    21     1    red  0.5        1    NA
      22 2  9.0    22     1    red  0.5        1    NA
      23 2  4.0    23     1    red  0.5        1    NA
      24 2  4.5    24     1    red  0.5        1    NA
      25 2  3.5    25     1    red  0.5        1    NA
      26 2 10.0    26     1    red  0.5        1    NA
      27 2  8.0    27     1    red  0.5        1    NA
      28 1 10.0     1     1    red  0.5        1    NA
      29 1 10.0     2     1    red  0.5        1    NA
      30 1 10.0     3     1    red  0.5        1    NA
      31 1  9.0     4     1    red  0.5        1    NA
      32 1  8.5     5     1    red  0.5        1    NA
      33 1  3.0     6     1    red  0.5        1    NA
      34 1 10.0     7     1    red  0.5        1    NA
      35 1 10.0     8     1    red  0.5        1    NA
      36 1 10.0     9     1    red  0.5        1    NA
      37 1  0.0    10     1    red  0.5        1    NA
      38 1 10.0    11     1    red  0.5        1    NA
      39 1  8.5    12     1    red  0.5        1    NA
      40 1  8.5    13     1    red  0.5        1    NA
      41 1  7.0    14     1    red  0.5        1    NA
      42 1  9.0    15     1    red  0.5        1    NA
      43 1 10.0    16     1    red  0.5        1    NA
      44 1 10.0    17     1    red  0.5        1    NA
      45 1 10.0    18     1    red  0.5        1    NA
      46 1  9.5    19     1    red  0.5        1    NA
      47 1  9.5    20     1    red  0.5        1    NA
      48 1  7.5    21     1    red  0.5        1    NA
      49 1 10.0    22     1    red  0.5        1    NA
      50 1  1.5    23     1    red  0.5        1    NA
      51 1  5.5    24     1    red  0.5        1    NA
      52 1 10.0    25     1    red  0.5        1    NA
      53 1 10.0    26     1    red  0.5        1    NA
      54 1  8.0    27     1    red  0.5        1    NA
      
      [[4]]
        x   y         label PANEL group colour  fill size angle alpha family fontface
      1 1 3.0        Europe     1     1  black white    3     0    NA               1
      2 1 0.0 North America     1     1  black white    3     0    NA               1
      3 1 1.5        Europe     1     1  black white    3     0    NA               1
      4 1 5.5        Europe     1     1  black white    3     0    NA               1
        lineheight hjust vjust point.size segment.linetype segment.size
      1        1.2   0.5   0.5          1                1          0.5
      2        1.2   0.5   0.5          1                1          0.5
      3        1.2   0.5   0.5          1                1          0.5
      4        1.2   0.5   0.5          1                1          0.5
        segment.curvature segment.angle segment.ncp segment.shape segment.square
      1                 0            90           1           0.5           TRUE
      2                 0            90           1           0.5           TRUE
      3                 0            90           1           0.5           TRUE
      4                 0            90           1           0.5           TRUE
        segment.squareShape segment.inflect segment.debug
      1                   1           FALSE         FALSE
      2                   1           FALSE         FALSE
      3                   1           FALSE         FALSE
      4                   1           FALSE         FALSE
      
      [[5]]
        x        y group PANEL colour size linetype alpha
      1 1 8.333333     1     1   blue    2        1   0.8
      2 2 7.074074     1     1   blue    2        1   0.8
      
      [[6]]
        x        y PANEL group shape    colour size fill alpha stroke
      1 1 8.333333     1     1    19 darkgreen    3   NA   0.5    0.5
      2 2 7.074074     1     2    19 darkgreen    3   NA   0.5    0.5
      

---

    Code
      within(pb1$plot$labels, rm(subtitle, caption))
    Output
      $x
      [1] "condition"
      
      $y
      [1] "desire"
      
      $colour
      [1] "condition"
      
      $title
      [1] "bugs dataset"
      
      $group
      [1] ".rowid"
      
      $label
      [1] "outlier.label"
      
      $alt
      [1] ""
      

# basic plotting works - more than two groups

    Code
      list(pb1$data[[1]], pb1$data[[2]], pb1$data[[4]], pb1$data[[5]], pb1$data[[6]],
      pb1$data[[7]])
    Output
      [[1]]
            colour x    y group PANEL shape size fill alpha stroke
      1  #1B9E77FF 1 5.40     1     1    19    3   NA   0.5    0.5
      2  #D95F02FF 2 5.50     1     1    19    3   NA   0.5    0.5
      3  #7570B3FF 3 5.55     1     1    19    3   NA   0.5    0.5
      4  #1B9E77FF 1 5.85     2     1    19    3   NA   0.5    0.5
      5  #D95F02FF 2 5.70     2     1    19    3   NA   0.5    0.5
      6  #7570B3FF 3 5.75     2     1    19    3   NA   0.5    0.5
      7  #1B9E77FF 1 5.20     3     1    19    3   NA   0.5    0.5
      8  #D95F02FF 2 5.60     3     1    19    3   NA   0.5    0.5
      9  #7570B3FF 3 5.50     3     1    19    3   NA   0.5    0.5
      10 #1B9E77FF 1 5.55     4     1    19    3   NA   0.5    0.5
      11 #D95F02FF 2 5.50     4     1    19    3   NA   0.5    0.5
      12 #7570B3FF 3 5.40     4     1    19    3   NA   0.5    0.5
      13 #1B9E77FF 1 5.90     5     1    19    3   NA   0.5    0.5
      14 #D95F02FF 2 5.85     5     1    19    3   NA   0.5    0.5
      15 #7570B3FF 3 5.70     5     1    19    3   NA   0.5    0.5
      16 #1B9E77FF 1 5.45     6     1    19    3   NA   0.5    0.5
      17 #D95F02FF 2 5.55     6     1    19    3   NA   0.5    0.5
      18 #7570B3FF 3 5.60     6     1    19    3   NA   0.5    0.5
      19 #1B9E77FF 1 5.40     7     1    19    3   NA   0.5    0.5
      20 #D95F02FF 2 5.40     7     1    19    3   NA   0.5    0.5
      21 #7570B3FF 3 5.35     7     1    19    3   NA   0.5    0.5
      22 #1B9E77FF 1 5.45     8     1    19    3   NA   0.5    0.5
      23 #D95F02FF 2 5.50     8     1    19    3   NA   0.5    0.5
      24 #7570B3FF 3 5.35     8     1    19    3   NA   0.5    0.5
      25 #1B9E77FF 1 5.25     9     1    19    3   NA   0.5    0.5
      26 #D95F02FF 2 5.15     9     1    19    3   NA   0.5    0.5
      27 #7570B3FF 3 5.00     9     1    19    3   NA   0.5    0.5
      28 #1B9E77FF 1 5.85    10     1    19    3   NA   0.5    0.5
      29 #D95F02FF 2 5.80    10     1    19    3   NA   0.5    0.5
      30 #7570B3FF 3 5.70    10     1    19    3   NA   0.5    0.5
      31 #1B9E77FF 1 5.25    11     1    19    3   NA   0.5    0.5
      32 #D95F02FF 2 5.20    11     1    19    3   NA   0.5    0.5
      33 #7570B3FF 3 5.10    11     1    19    3   NA   0.5    0.5
      34 #1B9E77FF 1 5.65    12     1    19    3   NA   0.5    0.5
      35 #D95F02FF 2 5.55    12     1    19    3   NA   0.5    0.5
      36 #7570B3FF 3 5.45    12     1    19    3   NA   0.5    0.5
      37 #1B9E77FF 1 5.60    13     1    19    3   NA   0.5    0.5
      38 #D95F02FF 2 5.35    13     1    19    3   NA   0.5    0.5
      39 #7570B3FF 3 5.45    13     1    19    3   NA   0.5    0.5
      40 #1B9E77FF 1 5.05    14     1    19    3   NA   0.5    0.5
      41 #D95F02FF 2 5.00    14     1    19    3   NA   0.5    0.5
      42 #7570B3FF 3 4.95    14     1    19    3   NA   0.5    0.5
      43 #1B9E77FF 1 5.50    15     1    19    3   NA   0.5    0.5
      44 #D95F02FF 2 5.50    15     1    19    3   NA   0.5    0.5
      45 #7570B3FF 3 5.40    15     1    19    3   NA   0.5    0.5
      46 #1B9E77FF 1 5.45    16     1    19    3   NA   0.5    0.5
      47 #D95F02FF 2 5.55    16     1    19    3   NA   0.5    0.5
      48 #7570B3FF 3 5.50    16     1    19    3   NA   0.5    0.5
      49 #1B9E77FF 1 5.55    17     1    19    3   NA   0.5    0.5
      50 #D95F02FF 2 5.55    17     1    19    3   NA   0.5    0.5
      51 #7570B3FF 3 5.35    17     1    19    3   NA   0.5    0.5
      52 #1B9E77FF 1 5.45    18     1    19    3   NA   0.5    0.5
      53 #D95F02FF 2 5.50    18     1    19    3   NA   0.5    0.5
      54 #7570B3FF 3 5.55    18     1    19    3   NA   0.5    0.5
      55 #1B9E77FF 1 5.50    19     1    19    3   NA   0.5    0.5
      56 #D95F02FF 2 5.45    19     1    19    3   NA   0.5    0.5
      57 #7570B3FF 3 5.25    19     1    19    3   NA   0.5    0.5
      58 #1B9E77FF 1 5.65    20     1    19    3   NA   0.5    0.5
      59 #D95F02FF 2 5.60    20     1    19    3   NA   0.5    0.5
      60 #7570B3FF 3 5.40    20     1    19    3   NA   0.5    0.5
      61 #1B9E77FF 1 5.70    21     1    19    3   NA   0.5    0.5
      62 #D95F02FF 2 5.65    21     1    19    3   NA   0.5    0.5
      63 #7570B3FF 3 5.55    21     1    19    3   NA   0.5    0.5
      64 #1B9E77FF 1 6.30    22     1    19    3   NA   0.5    0.5
      65 #D95F02FF 2 6.30    22     1    19    3   NA   0.5    0.5
      66 #7570B3FF 3 6.25    22     1    19    3   NA   0.5    0.5
      
      [[2]]
        ymin  lower middle upper ymax                     outliers notchupper
      1 5.20 5.4125  5.500  5.65 5.90                   5.05, 6.30   5.580004
      2 5.35 5.4625  5.525  5.60 5.80 5.85, 5.15, 5.20, 5.00, 6.30   5.571318
      3 5.10 5.3500  5.450  5.55 5.75             5.00, 4.95, 6.25   5.517371
        notchlower x flipped_aes PANEL group ymin_final ymax_final xmin xmax xid newx
      1   5.419996 1       FALSE     1     1       5.05       6.30  0.9  1.1   1    1
      2   5.478682 2       FALSE     1     2       5.00       6.30  1.9  2.1   2    2
      3   5.382629 3       FALSE     1     3       4.95       6.25  2.9  3.1   3    3
        new_width weight colour  fill size alpha shape linetype
      1       0.2      1 grey20 white  0.5   0.5    19    solid
      2       0.2      1 grey20 white  0.5   0.5    19    solid
      3       0.2      1 grey20 white  0.5   0.5    19    solid
      
      [[3]]
        x    y label PANEL group colour  fill size angle alpha family fontface
      1 2 5.00  5.00     1     2  black white    3     0    NA               1
      2 1 6.30  6.30     1     1  black white    3     0    NA               1
      3 2 6.30  6.30     1     2  black white    3     0    NA               1
      4 3 6.25  6.25     1     3  black white    3     0    NA               1
        lineheight hjust vjust point.size segment.linetype segment.size
      1        1.2   0.5   0.5          1                1          0.5
      2        1.2   0.5   0.5          1                1          0.5
      3        1.2   0.5   0.5          1                1          0.5
      4        1.2   0.5   0.5          1                1          0.5
        segment.curvature segment.angle segment.ncp segment.shape segment.square
      1                 0            90           1           0.5           TRUE
      2                 0            90           1           0.5           TRUE
      3                 0            90           1           0.5           TRUE
      4                 0            90           1           0.5           TRUE
        segment.squareShape segment.inflect segment.debug
      1                   1           FALSE         FALSE
      2                   1           FALSE         FALSE
      3                   1           FALSE         FALSE
      4                   1           FALSE         FALSE
      
      [[4]]
        x        y group PANEL colour size linetype alpha
      1 1 5.543182     1     1    red    1        1   0.5
      2 2 5.534091     1     1    red    1        1   0.5
      3 3 5.459091     1     1    red    1        1   0.5
      
      [[5]]
        x        y PANEL group shape  colour size fill alpha stroke
      1 1 5.543182     1     1    19 darkred    5   NA    NA    0.5
      2 2 5.534091     1     2    19 darkred    5   NA    NA    0.5
      3 3 5.459091     1     3    19 darkred    5   NA    NA    0.5
      
      [[6]]
        x        y                       label PANEL group nudge_x  nudge_y colour
      1 1 5.543182 widehat(mu)[mean]=='5.5432'     1     1     1.4 5.543182  black
      2 2 5.534091 widehat(mu)[mean]=='5.5341'     1     2     2.4 5.534091  black
      3 3 5.459091 widehat(mu)[mean]=='5.4591'     1     3     3.4 5.459091  black
         fill size angle alpha family fontface lineheight hjust vjust point.size
      1 white    3     0    NA               1        1.2   0.5   0.5          1
      2 white    3     0    NA               1        1.2   0.5   0.5          1
      3 white    3     0    NA               1        1.2   0.5   0.5          1
        segment.linetype segment.size segment.curvature segment.angle segment.ncp
      1                4          0.5                 0            90           1
      2                4          0.5                 0            90           1
      3                4          0.5                 0            90           1
        segment.shape segment.square segment.squareShape segment.inflect
      1           0.5           TRUE                   1           FALSE
      2           0.5           TRUE                   1           FALSE
      3           0.5           TRUE                   1           FALSE
        segment.debug
      1         FALSE
      2         FALSE
      3         FALSE
      

---

    Code
      within(pb1$plot$labels, rm(subtitle, caption))
    Output
      $x
      [1] "Wine"
      
      $y
      [1] "Taste"
      
      $colour
      [1] "Wine"
      
      $title
      [1] "wine tasting data"
      
      $group
      [1] ".rowid"
      
      $label
      [1] "outlier.label"
      
      $alt
      [1] ""
      

# checking subtitle outputs - without NAs

    Code
      within(pb1$plot$labels, rm(subtitle))
    Output
      $x
      [1] "condition"
      
      $y
      [1] "value"
      
      $colour
      [1] "condition"
      
      $title
      NULL
      
      $caption
      expression(list("Pairwise test:" ~ bold("Durbin-Conover test"), 
          "Comparisons shown:" ~ bold("only significant")))
      
      $group
      [1] ".rowid"
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

---

    Code
      within(pb2$plot$labels, rm(subtitle))
    Output
      $x
      [1] "condition"
      
      $y
      [1] "value"
      
      $colour
      [1] "condition"
      
      $title
      NULL
      
      $caption
      expression(list("Pairwise test:" ~ bold("Yuen's trimmed means test"), 
          "Comparisons shown:" ~ bold("only non-significant")))
      
      $group
      [1] ".rowid"
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

---

    Code
      within(pb3$plot$labels, rm(subtitle))
    Output
      $x
      [1] "condition"
      
      $y
      [1] "desire"
      
      $colour
      [1] "condition"
      
      $title
      NULL
      
      $caption
      NULL
      
      $group
      [1] ".rowid"
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

---

    Code
      within(pb4$plot$labels, rm(subtitle))
    Output
      $x
      [1] "condition"
      
      $y
      [1] "desire"
      
      $colour
      [1] "condition"
      
      $title
      NULL
      
      $caption
      NULL
      
      $group
      [1] ".rowid"
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

# ggplot component addition works

    Code
      pb$plot$labels
    Output
      $y
      [1] "Taste rating"
      
      $x
      [1] "Wine"
      
      $colour
      [1] "Wine"
      
      $title
      NULL
      
      $subtitle
      NULL
      
      $caption
      expression()
      
      $group
      [1] ".rowid"
      
      $label
      [1] "expression"
      
      $alt
      [1] ""
      

# `pairwise_comparisons()` works for within-subjects design - NAs

    Code
      list(df1, df2, df3)
    Output
      [[1]]
      # A tibble: 6 x 6
        group1 group2  p.value test.details     p.value.adjustment
        <chr>  <chr>     <dbl> <chr>            <chr>             
      1 HDHF   HDLF   3.18e- 3 Student's t-test Bonferroni        
      2 HDHF   LDHF   4.21e- 1 Student's t-test Bonferroni        
      3 HDHF   LDLF   3.95e-12 Student's t-test Bonferroni        
      4 HDLF   LDHF   3.37e- 1 Student's t-test Bonferroni        
      5 HDLF   LDLF   7.94e- 3 Student's t-test Bonferroni        
      6 LDHF   LDLF   1.33e- 8 Student's t-test Bonferroni        
        label                                            
        <chr>                                            
      1 list(~italic(p)[Bonferroni-corrected]==0.003)    
      2 list(~italic(p)[Bonferroni-corrected]==0.421)    
      3 list(~italic(p)[Bonferroni-corrected]==3.950e-12)
      4 list(~italic(p)[Bonferroni-corrected]==0.337)    
      5 list(~italic(p)[Bonferroni-corrected]==0.008)    
      6 list(~italic(p)[Bonferroni-corrected]==1.331e-08)
      
      [[2]]
      # A tibble: 6 x 11
        group1 group2 statistic  p.value alternative
        <chr>  <chr>      <dbl>    <dbl> <chr>      
      1 HDHF   HDLF        4.78 1.44e- 5 two.sided  
      2 HDHF   LDHF        2.44 4.47e- 2 two.sided  
      3 HDHF   LDLF        8.01 5.45e-13 two.sided  
      4 HDLF   LDHF        2.34 4.96e- 2 two.sided  
      5 HDLF   LDLF        3.23 5.05e- 3 two.sided  
      6 LDHF   LDLF        5.57 4.64e- 7 two.sided  
        method                                                                
        <chr>                                                                 
      1 Durbin's all-pairs test for a two-way balanced incomplete block design
      2 Durbin's all-pairs test for a two-way balanced incomplete block design
      3 Durbin's all-pairs test for a two-way balanced incomplete block design
      4 Durbin's all-pairs test for a two-way balanced incomplete block design
      5 Durbin's all-pairs test for a two-way balanced incomplete block design
      6 Durbin's all-pairs test for a two-way balanced incomplete block design
        distribution p.adjustment test.details        p.value.adjustment
        <chr>        <chr>        <chr>               <chr>             
      1 t            none         Durbin-Conover test BY                
      2 t            none         Durbin-Conover test BY                
      3 t            none         Durbin-Conover test BY                
      4 t            none         Durbin-Conover test BY                
      5 t            none         Durbin-Conover test BY                
      6 t            none         Durbin-Conover test BY                
        label                                    
        <chr>                                    
      1 list(~italic(p)[BY-corrected]==1.436e-05)
      2 list(~italic(p)[BY-corrected]==0.045)    
      3 list(~italic(p)[BY-corrected]==5.447e-13)
      4 list(~italic(p)[BY-corrected]==0.050)    
      5 list(~italic(p)[BY-corrected]==0.005)    
      6 list(~italic(p)[BY-corrected]==4.635e-07)
      
      [[3]]
      # A tibble: 6 x 11
        group1 group2 estimate conf.level conf.low conf.high     p.value  p.crit
        <chr>  <chr>     <dbl>      <dbl>    <dbl>     <dbl>       <dbl>   <dbl>
      1 HDHF   HDLF      1.03        0.95   0.140      1.92  0.00999     0.0127 
      2 HDHF   LDHF      0.454       0.95  -0.104      1.01  0.0520      0.025  
      3 HDHF   LDLF      1.95        0.95   1.09       2.82  0.000000564 0.00851
      4 HDLF   LDHF     -0.676       0.95  -1.61       0.256 0.0520      0.05   
      5 HDLF   LDLF      0.889       0.95   0.0244     1.75  0.0203      0.0169 
      6 LDHF   LDLF      1.35        0.95   0.560      2.14  0.000102    0.0102 
        test.details              p.value.adjustment
        <chr>                     <chr>             
      1 Yuen's trimmed means test Hommel            
      2 Yuen's trimmed means test Hommel            
      3 Yuen's trimmed means test Hommel            
      4 Yuen's trimmed means test Hommel            
      5 Yuen's trimmed means test Hommel            
      6 Yuen's trimmed means test Hommel            
        label                                        
        <chr>                                        
      1 list(~italic(p)[Hommel-corrected]==0.010)    
      2 list(~italic(p)[Hommel-corrected]==0.052)    
      3 list(~italic(p)[Hommel-corrected]==5.642e-07)
      4 list(~italic(p)[Hommel-corrected]==0.052)    
      5 list(~italic(p)[Hommel-corrected]==0.020)    
      6 list(~italic(p)[Hommel-corrected]==1.017e-04)
      

# `pairwise_comparisons()` works for within-subjects design - without NAs

    Code
      list(df1, df2, df3)
    Output
      [[1]]
      # A tibble: 3 x 6
        group1 group2  p.value test.details     p.value.adjustment
        <chr>  <chr>     <dbl> <chr>            <chr>             
      1 Wine A Wine B 0.732    Student's t-test None              
      2 Wine A Wine C 0.0142   Student's t-test None              
      3 Wine B Wine C 0.000675 Student's t-test None              
        label                                   
        <chr>                                   
      1 list(~italic(p)[uncorrected]==0.732)    
      2 list(~italic(p)[uncorrected]==0.014)    
      3 list(~italic(p)[uncorrected]==6.754e-04)
      
      [[2]]
      # A tibble: 3 x 11
        group1 group2 statistic  p.value alternative
        <chr>  <chr>      <dbl>    <dbl> <chr>      
      1 Wine A Wine B      1.05 0.301    two.sided  
      2 Wine A Wine C      3.66 0.000691 two.sided  
      3 Wine B Wine C      2.62 0.0123   two.sided  
        method                                                                
        <chr>                                                                 
      1 Durbin's all-pairs test for a two-way balanced incomplete block design
      2 Durbin's all-pairs test for a two-way balanced incomplete block design
      3 Durbin's all-pairs test for a two-way balanced incomplete block design
        distribution p.adjustment test.details        p.value.adjustment
        <chr>        <chr>        <chr>               <chr>             
      1 t            none         Durbin-Conover test None              
      2 t            none         Durbin-Conover test None              
      3 t            none         Durbin-Conover test None              
        label                                   
        <chr>                                   
      1 list(~italic(p)[uncorrected]==0.301)    
      2 list(~italic(p)[uncorrected]==6.915e-04)
      3 list(~italic(p)[uncorrected]==0.012)    
      
      [[3]]
      # A tibble: 3 x 11
        group1 group2 estimate conf.level conf.low conf.high p.value p.crit
        <chr>  <chr>     <dbl>      <dbl>    <dbl>     <dbl>   <dbl>  <dbl>
      1 Wine A Wine B   0.0214       0.95 -0.0216     0.0645 0.195   0.05  
      2 Wine A Wine C   0.114        0.95  0.0215     0.207  0.00492 0.0169
      3 Wine B Wine C   0.0821       0.95  0.00891    0.155  0.00878 0.025 
        test.details              p.value.adjustment
        <chr>                     <chr>             
      1 Yuen's trimmed means test None              
      2 Yuen's trimmed means test None              
      3 Yuen's trimmed means test None              
        label                               
        <chr>                               
      1 list(~italic(p)[uncorrected]==0.195)
      2 list(~italic(p)[uncorrected]==0.005)
      3 list(~italic(p)[uncorrected]==0.009)
      

# checking ggcorrmat - without NAs - pearson's r

    Code
      pb$data
    Output
      [[1]]
            fill x y PANEL group xmin xmax ymin ymax colour size linetype alpha width
      1  #7570B3 1 1     1     1  0.5  1.5  0.5  1.5   gray  0.1        1    NA    NA
      2  #CC6A19 2 1     1     5  1.5  2.5  0.5  1.5   gray  0.1        1    NA    NA
      3  #8B6E9F 3 1     1     9  2.5  3.5  0.5  1.5   gray  0.1        1    NA    NA
      4  #936D97 4 1     1    13  3.5  4.5  0.5  1.5   gray  0.1        1    NA    NA
      5  #CC6A19 1 2     1     2  0.5  1.5  1.5  2.5   gray  0.1        1    NA    NA
      6  #7570B3 2 2     1     6  1.5  2.5  1.5  2.5   gray  0.1        1    NA    NA
      7  #A7813E 3 2     1    10  2.5  3.5  1.5  2.5   gray  0.1        1    NA    NA
      8  #AF7D37 4 2     1    14  3.5  4.5  1.5  2.5   gray  0.1        1    NA    NA
      9  #8B6E9F 1 3     1     3  0.5  1.5  2.5  3.5   gray  0.1        1    NA    NA
      10 #A7813E 2 3     1     7  1.5  2.5  2.5  3.5   gray  0.1        1    NA    NA
      11 #7570B3 3 3     1    11  2.5  3.5  2.5  3.5   gray  0.1        1    NA    NA
      12 #7C6FAD 4 3     1    15  3.5  4.5  2.5  3.5   gray  0.1        1    NA    NA
      13 #936D97 1 4     1     4  0.5  1.5  3.5  4.5   gray  0.1        1    NA    NA
      14 #AF7D37 2 4     1     8  1.5  2.5  3.5  4.5   gray  0.1        1    NA    NA
      15 #7C6FAD 3 4     1    12  2.5  3.5  3.5  4.5   gray  0.1        1    NA    NA
      16 #7570B3 4 4     1    16  3.5  4.5  3.5  4.5   gray  0.1        1    NA    NA
         height
      1      NA
      2      NA
      3      NA
      4      NA
      5      NA
      6      NA
      7      NA
      8      NA
      9      NA
      10     NA
      11     NA
      12     NA
      13     NA
      14     NA
      15     NA
      16     NA
      
      [[2]]
            fill x y PANEL group colour size angle hjust vjust alpha family fontface
      1  #7570B3 1 1     1     1  white    4     0   0.5   0.5    NA               1
      2  #CC6A19 2 1     1     5  white    4     0   0.5   0.5    NA               1
      3  #8B6E9F 3 1     1     9  white    4     0   0.5   0.5    NA               1
      4  #936D97 4 1     1    13  white    4     0   0.5   0.5    NA               1
      5  #CC6A19 1 2     1     2  white    4     0   0.5   0.5    NA               1
      6  #7570B3 2 2     1     6  white    4     0   0.5   0.5    NA               1
      7  #A7813E 3 2     1    10  white    4     0   0.5   0.5    NA               1
      8  #AF7D37 4 2     1    14  white    4     0   0.5   0.5    NA               1
      9  #8B6E9F 1 3     1     3  white    4     0   0.5   0.5    NA               1
      10 #A7813E 2 3     1     7  white    4     0   0.5   0.5    NA               1
      11 #7570B3 3 3     1    11  white    4     0   0.5   0.5    NA               1
      12 #7C6FAD 4 3     1    15  white    4     0   0.5   0.5    NA               1
      13 #936D97 1 4     1     4  white    4     0   0.5   0.5    NA               1
      14 #AF7D37 2 4     1     8  white    4     0   0.5   0.5    NA               1
      15 #7C6FAD 3 4     1    12  white    4     0   0.5   0.5    NA               1
      16 #7570B3 4 4     1    16  white    4     0   0.5   0.5    NA               1
         lineheight   label
      1         1.2  1.0000
      2         1.2 -0.1176
      3         1.2  0.8718
      4         1.2  0.8179
      5         1.2 -0.1176
      6         1.2  1.0000
      7         1.2 -0.4284
      8         1.2 -0.3661
      9         1.2  0.8718
      10        1.2 -0.4284
      11        1.2  1.0000
      12        1.2  0.9629
      13        1.2  0.8179
      14        1.2 -0.3661
      15        1.2  0.9629
      16        1.2  1.0000
      
      [[3]]
           fill x y PANEL group shape colour size alpha stroke
      1 #D0622D 2 1     1     2 cross  white    5    NA    0.5
      2 #D0622D 1 2     1     1 cross  white    5    NA    0.5
      

---

    Code
      list(p$labels, pb$plot$plot_env$legend.title)
    Output
      [[1]]
      [[1]]$xlab
      NULL
      
      [[1]]$ylab
      NULL
      
      [[1]]$title
      [1] "Iris dataset"
      
      [[1]]$subtitle
      [1] "By Edgar Anderson"
      
      [[1]]$caption
      atop(displaystyle(NULL), expr = paste(bold("X"), " = non-significant at ", 
          italic("p"), " < ", 0.001, " (Adjustment: ", "FDR", ")"))
      
      [[1]]$fill
      [1] "value"
      
      [[1]]$x
      [1] "Var1"
      
      [[1]]$y
      [1] "Var2"
      
      
      [[2]]
      atop(atop(scriptstyle(bold("sample sizes:")), italic(n) ~ "=" ~ 
          "150"), atop(scriptstyle(bold("correlation:")), "Pearson"))
      

# checking ggcorrmat - without NAs - robust r

    Code
      pb$data
    Output
      [[1]]
            fill x y PANEL group xmin xmax ymin ymax colour size linetype alpha width
      1  #009E73 1 1     1     1  0.5  1.5  0.5  1.5  black  0.1        1    NA    NA
      2  #009E73 1 2     1     2  0.5  1.5  1.5  2.5  black  0.1        1    NA    NA
      3  #009E73 2 2     1     8  1.5  2.5  1.5  2.5  black  0.1        1    NA    NA
      4  #FFFAF4 1 3     1     3  0.5  1.5  2.5  3.5  black  0.1        1    NA    NA
      5  #FFFAF4 2 3     1     9  1.5  2.5  2.5  3.5  black  0.1        1    NA    NA
      6  #FFFAF4 3 3     1    14  2.5  3.5  2.5  3.5  black  0.1        1    NA    NA
      7  #C1E2D3 1 4     1     4  0.5  1.5  3.5  4.5  black  0.1        1    NA    NA
      8  #C1E2D3 2 4     1    10  1.5  2.5  3.5  4.5  black  0.1        1    NA    NA
      9  #C1E2D3 3 4     1    15  2.5  3.5  3.5  4.5  black  0.1        1    NA    NA
      10 #DEF0E7 4 4     1    19  3.5  4.5  3.5  4.5  black  0.1        1    NA    NA
      11 #FFE1BB 1 5     1     5  0.5  1.5  4.5  5.5  black  0.1        1    NA    NA
      12 #FFE1BB 2 5     1    11  1.5  2.5  4.5  5.5  black  0.1        1    NA    NA
      13 #FFE1BB 3 5     1    16  2.5  3.5  4.5  5.5  black  0.1        1    NA    NA
      14 #FED6A3 4 5     1    20  3.5  4.5  4.5  5.5  black  0.1        1    NA    NA
      15 #FFFDFA 5 5     1    23  4.5  5.5  4.5  5.5  black  0.1        1    NA    NA
      16 #FFFBF6 1 6     1     6  0.5  1.5  5.5  6.5  black  0.1        1    NA    NA
      17 #FFFBF6 2 6     1    12  1.5  2.5  5.5  6.5  black  0.1        1    NA    NA
      18 #FFFBF6 3 6     1    17  2.5  3.5  5.5  6.5  black  0.1        1    NA    NA
      19 #F7FBF9 4 6     1    21  3.5  4.5  5.5  6.5  black  0.1        1    NA    NA
      20 #FAC882 5 6     1    24  4.5  5.5  5.5  6.5  black  0.1        1    NA    NA
      21 #FCCF91 6 6     1    26  5.5  6.5  5.5  6.5  black  0.1        1    NA    NA
      22 #EAF5F0 1 7     1     7  0.5  1.5  6.5  7.5  black  0.1        1    NA    NA
      23 #EAF5F0 2 7     1    13  1.5  2.5  6.5  7.5  black  0.1        1    NA    NA
      24 #EAF5F0 3 7     1    18  2.5  3.5  6.5  7.5  black  0.1        1    NA    NA
      25 #60B794 4 7     1    22  3.5  4.5  6.5  7.5  black  0.1        1    NA    NA
      26 #FFDBAE 5 7     1    25  4.5  5.5  6.5  7.5  black  0.1        1    NA    NA
      27 #EDF6F2 6 7     1    27  5.5  6.5  6.5  7.5  black  0.1        1    NA    NA
      28 #FFE7C9 7 7     1    28  6.5  7.5  6.5  7.5  black  0.1        1    NA    NA
         height
      1      NA
      2      NA
      3      NA
      4      NA
      5      NA
      6      NA
      7      NA
      8      NA
      9      NA
      10     NA
      11     NA
      12     NA
      13     NA
      14     NA
      15     NA
      16     NA
      17     NA
      18     NA
      19     NA
      20     NA
      21     NA
      22     NA
      23     NA
      24     NA
      25     NA
      26     NA
      27     NA
      28     NA
      
      [[2]]
            fill x y PANEL group colour size angle hjust vjust alpha family fontface
      1  #009E73 1 1     1     1  black    4     0   0.5   0.5    NA               1
      2  #009E73 1 2     1     2  black    4     0   0.5   0.5    NA               1
      3  #009E73 2 2     1     8  black    4     0   0.5   0.5    NA               1
      4  #FFFAF4 1 3     1     3  black    4     0   0.5   0.5    NA               1
      5  #FFFAF4 2 3     1     9  black    4     0   0.5   0.5    NA               1
      6  #FFFAF4 3 3     1    14  black    4     0   0.5   0.5    NA               1
      7  #C1E2D3 1 4     1     4  black    4     0   0.5   0.5    NA               1
      8  #C1E2D3 2 4     1    10  black    4     0   0.5   0.5    NA               1
      9  #C1E2D3 3 4     1    15  black    4     0   0.5   0.5    NA               1
      10 #DEF0E7 4 4     1    19  black    4     0   0.5   0.5    NA               1
      11 #FFE1BB 1 5     1     5  black    4     0   0.5   0.5    NA               1
      12 #FFE1BB 2 5     1    11  black    4     0   0.5   0.5    NA               1
      13 #FFE1BB 3 5     1    16  black    4     0   0.5   0.5    NA               1
      14 #FED6A3 4 5     1    20  black    4     0   0.5   0.5    NA               1
      15 #FFFDFA 5 5     1    23  black    4     0   0.5   0.5    NA               1
      16 #FFFBF6 1 6     1     6  black    4     0   0.5   0.5    NA               1
      17 #FFFBF6 2 6     1    12  black    4     0   0.5   0.5    NA               1
      18 #FFFBF6 3 6     1    17  black    4     0   0.5   0.5    NA               1
      19 #F7FBF9 4 6     1    21  black    4     0   0.5   0.5    NA               1
      20 #FAC882 5 6     1    24  black    4     0   0.5   0.5    NA               1
      21 #FCCF91 6 6     1    26  black    4     0   0.5   0.5    NA               1
      22 #EAF5F0 1 7     1     7  black    4     0   0.5   0.5    NA               1
      23 #EAF5F0 2 7     1    13  black    4     0   0.5   0.5    NA               1
      24 #EAF5F0 3 7     1    18  black    4     0   0.5   0.5    NA               1
      25 #60B794 4 7     1    22  black    4     0   0.5   0.5    NA               1
      26 #FFDBAE 5 7     1    25  black    4     0   0.5   0.5    NA               1
      27 #EDF6F2 6 7     1    27  black    4     0   0.5   0.5    NA               1
      28 #FFE7C9 7 7     1    28  black    4     0   0.5   0.5    NA               1
         lineheight label
      1         1.2  1.00
      2         1.2  1.00
      3         1.2  1.00
      4         1.2 -0.05
      5         1.2 -0.05
      6         1.2 -0.05
      7         1.2  0.30
      8         1.2  0.30
      9         1.2  0.30
      10        1.2  0.16
      11        1.2 -0.30
      12        1.2 -0.30
      13        1.2 -0.30
      14        1.2 -0.41
      15        1.2 -0.02
      16        1.2 -0.04
      17        1.2 -0.04
      18        1.2 -0.04
      19        1.2  0.04
      20        1.2 -0.56
      21        1.2 -0.49
      22        1.2  0.10
      23        1.2  0.10
      24        1.2  0.10
      25        1.2  0.75
      26        1.2 -0.36
      27        1.2  0.09
      28        1.2 -0.24
      
      [[3]]
            fill x y PANEL group shape colour size alpha stroke
      1  #009E73 1 3     1     1 cross  black   14    NA    0.5
      2  #009E73 2 3     1     6 cross  black   14    NA    0.5
      3  #009E73 3 3     1    11 cross  black   14    NA    0.5
      4  #009E73 1 4     1     2 cross  black   14    NA    0.5
      5  #009E73 2 4     1     7 cross  black   14    NA    0.5
      6  #009E73 3 4     1    12 cross  black   14    NA    0.5
      7  #009E73 4 4     1    16 cross  black   14    NA    0.5
      8  #009E73 1 5     1     3 cross  black   14    NA    0.5
      9  #009E73 2 5     1     8 cross  black   14    NA    0.5
      10 #009E73 3 5     1    13 cross  black   14    NA    0.5
      11 #009E73 4 5     1    17 cross  black   14    NA    0.5
      12 #009E73 5 5     1    20 cross  black   14    NA    0.5
      13 #009E73 1 6     1     4 cross  black   14    NA    0.5
      14 #009E73 2 6     1     9 cross  black   14    NA    0.5
      15 #009E73 3 6     1    14 cross  black   14    NA    0.5
      16 #009E73 4 6     1    18 cross  black   14    NA    0.5
      17 #009E73 5 6     1    21 cross  black   14    NA    0.5
      18 #009E73 6 6     1    23 cross  black   14    NA    0.5
      19 #009E73 1 7     1     5 cross  black   14    NA    0.5
      20 #009E73 2 7     1    10 cross  black   14    NA    0.5
      21 #009E73 3 7     1    15 cross  black   14    NA    0.5
      22 #D6ECE2 4 7     1    19 cross  black   14    NA    0.5
      23 #009E73 5 7     1    22 cross  black   14    NA    0.5
      24 #009E73 6 7     1    24 cross  black   14    NA    0.5
      25 #009E73 7 7     1    25 cross  black   14    NA    0.5
      

---

    Code
      list(p$labels, pb$plot$plot_env$legend.title)
    Output
      [[1]]
      [[1]]$xlab
      NULL
      
      [[1]]$ylab
      NULL
      
      [[1]]$title
      NULL
      
      [[1]]$subtitle
      NULL
      
      [[1]]$caption
      atop(displaystyle(NULL), expr = paste(bold("X"), " = non-significant at ", 
          italic("p"), " < ", 0.05, " (Adjustment: ", "Holm", ")"))
      
      [[1]]$fill
      [1] "value"
      
      [[1]]$x
      [1] "Var1"
      
      [[1]]$y
      [1] "Var2"
      
      
      [[2]]
      atop(atop(scriptstyle(bold("sample sizes:")), italic(n) ~ "=" ~ 
          "11"), atop(scriptstyle(bold("correlation (partial):")), 
          "Winsorized Pearson"))
      

# checking ggcorrmat - with NAs - robust r - partial

    Code
      pb$data
    Output
      [[1]]
            fill x y PANEL group xmin xmax ymin ymax colour size linetype alpha width
      1  #D6ECE2 1 1     1     1  0.5  1.5  0.5  1.5  black  0.1        1    NA    NA
      2  #FFFFFF 1 2     1     2  0.5  1.5  1.5  2.5  black  0.1        1    NA    NA
      3  #FFF7ED 2 2     1     6  1.5  2.5  1.5  2.5  black  0.1        1    NA    NA
      4  #E69F00 1 3     1     3  0.5  1.5  2.5  3.5  black  0.1        1    NA    NA
      5  #D8EDE3 2 3     1     7  1.5  2.5  2.5  3.5  black  0.1        1    NA    NA
      6  #FFFDFA 3 3     1    10  2.5  3.5  2.5  3.5  black  0.1        1    NA    NA
      7  #DAEEE4 1 4     1     4  0.5  1.5  3.5  4.5  black  0.1        1    NA    NA
      8  #F3F9F6 2 4     1     8  1.5  2.5  3.5  4.5  black  0.1        1    NA    NA
      9  #77C0A2 3 4     1    11  2.5  3.5  3.5  4.5  black  0.1        1    NA    NA
      10 #C8E5D7 4 4     1    13  3.5  4.5  3.5  4.5  black  0.1        1    NA    NA
      11 #C8E5D7 1 5     1     5  0.5  1.5  4.5  5.5  black  0.1        1    NA    NA
      12 #EAF5F0 2 5     1     9  1.5  2.5  4.5  5.5  black  0.1        1    NA    NA
      13 #FFFDFA 3 5     1    12  2.5  3.5  4.5  5.5  black  0.1        1    NA    NA
      14 #FFF6EA 4 5     1    14  3.5  4.5  4.5  5.5  black  0.1        1    NA    NA
      15 #B1DBC8 5 5     1    15  4.5  5.5  4.5  5.5  black  0.1        1    NA    NA
         height
      1      NA
      2      NA
      3      NA
      4      NA
      5      NA
      6      NA
      7      NA
      8      NA
      9      NA
      10     NA
      11     NA
      12     NA
      13     NA
      14     NA
      15     NA
      
      [[2]]
            fill x y PANEL group colour size angle hjust vjust alpha family fontface
      1  #D6ECE2 1 1     1     1  black    4     0   0.5   0.5    NA               1
      2  #FFFFFF 1 2     1     2  black    4     0   0.5   0.5    NA               1
      3  #FFF7ED 2 2     1     6  black    4     0   0.5   0.5    NA               1
      4  #E69F00 1 3     1     3  black    4     0   0.5   0.5    NA               1
      5  #D8EDE3 2 3     1     7  black    4     0   0.5   0.5    NA               1
      6  #FFFDFA 3 3     1    10  black    4     0   0.5   0.5    NA               1
      7  #DAEEE4 1 4     1     4  black    4     0   0.5   0.5    NA               1
      8  #F3F9F6 2 4     1     8  black    4     0   0.5   0.5    NA               1
      9  #77C0A2 3 4     1    11  black    4     0   0.5   0.5    NA               1
      10 #C8E5D7 4 4     1    13  black    4     0   0.5   0.5    NA               1
      11 #C8E5D7 1 5     1     5  black    4     0   0.5   0.5    NA               1
      12 #EAF5F0 2 5     1     9  black    4     0   0.5   0.5    NA               1
      13 #FFFDFA 3 5     1    12  black    4     0   0.5   0.5    NA               1
      14 #FFF6EA 4 5     1    14  black    4     0   0.5   0.5    NA               1
      15 #B1DBC8 5 5     1    15  black    4     0   0.5   0.5    NA               1
         lineheight label
      1         1.2  0.20
      2         1.2  0.00
      3         1.2 -0.08
      4         1.2 -1.00
      5         1.2  0.19
      6         1.2 -0.02
      7         1.2  0.18
      8         1.2  0.06
      9         1.2  0.65
      10        1.2  0.27
      11        1.2  0.27
      12        1.2  0.10
      13        1.2 -0.02
      14        1.2 -0.09
      15        1.2  0.38
      
      [[3]]
            fill x y PANEL group shape colour size alpha stroke
      1  #0F9F75 1 1     1     1 cross  black   14    NA    0.5
      2  #0F9F75 1 2     1     2 cross  black   14    NA    0.5
      3  #0F9F75 2 2     1     5 cross  black   14    NA    0.5
      4  #0F9F75 2 3     1     6 cross  black   14    NA    0.5
      5  #0F9F75 3 3     1     9 cross  black   14    NA    0.5
      6  #0F9F75 1 4     1     3 cross  black   14    NA    0.5
      7  #0F9F75 2 4     1     7 cross  black   14    NA    0.5
      8  #0F9F75 4 4     1    11 cross  black   14    NA    0.5
      9  #0F9F75 1 5     1     4 cross  black   14    NA    0.5
      10 #0F9F75 2 5     1     8 cross  black   14    NA    0.5
      11 #0F9F75 3 5     1    10 cross  black   14    NA    0.5
      12 #0F9F75 4 5     1    12 cross  black   14    NA    0.5
      13 #93CDB4 5 5     1    13 cross  black   14    NA    0.5
      

---

    Code
      list(p$labels, pb$plot$plot_env$legend.title)
    Output
      [[1]]
      [[1]]$caption
      NULL
      
      [[1]]$xlab
      NULL
      
      [[1]]$ylab
      NULL
      
      [[1]]$title
      NULL
      
      [[1]]$subtitle
      NULL
      
      [[1]]$fill
      [1] "value"
      
      [[1]]$x
      [1] "Var1"
      
      [[1]]$y
      [1] "Var2"
      
      
      [[2]]
      atop(atop(scriptstyle(bold("sample sizes:")), italic(n) ~ "=" ~ 
          "30"), atop(scriptstyle(bold("correlation (partial):")), 
          "Winsorized Pearson"))
      

# checking ggcorrmat - with NAs - spearman's rho

    Code
      pb$data
    Output
      [[1]]
            fill x y PANEL group xmin xmax ymin ymax colour size linetype alpha width
      1  #0B775E 1 1     1     1  0.5  1.5  0.5  1.5  black  0.1        1    NA    NA
      2  #57896B 2 1     1     5  1.5  2.5  0.5  1.5  black  0.1        1    NA    NA
      3  #E6BE81 3 1     1     9  2.5  3.5  0.5  1.5  black  0.1        1    NA    NA
      4  #E1BD6D 4 1     1    13  3.5  4.5  0.5  1.5  black  0.1        1    NA    NA
      5  #57896B 1 2     1     2  0.5  1.5  1.5  2.5  black  0.1        1    NA    NA
      6  #0B775E 2 2     1     6  1.5  2.5  1.5  2.5  black  0.1        1    NA    NA
      7  #E7BE87 3 2     1    10  2.5  3.5  1.5  2.5  black  0.1        1    NA    NA
      8  #E3BD77 4 2     1    14  3.5  4.5  1.5  2.5  black  0.1        1    NA    NA
      9  #E6BE81 1 3     1     3  0.5  1.5  2.5  3.5  black  0.1        1    NA    NA
      10 #E7BE87 2 3     1     7  1.5  2.5  2.5  3.5  black  0.1        1    NA    NA
      11 #0B775E 3 3     1    11  2.5  3.5  2.5  3.5  black  0.1        1    NA    NA
      12 #8E9C79 4 3     1    15  3.5  4.5  2.5  3.5  black  0.1        1    NA    NA
      13 #E1BD6D 1 4     1     4  0.5  1.5  3.5  4.5  black  0.1        1    NA    NA
      14 #E3BD77 2 4     1     8  1.5  2.5  3.5  4.5  black  0.1        1    NA    NA
      15 #8E9C79 3 4     1    12  2.5  3.5  3.5  4.5  black  0.1        1    NA    NA
      16 #0B775E 4 4     1    16  3.5  4.5  3.5  4.5  black  0.1        1    NA    NA
         height
      1      NA
      2      NA
      3      NA
      4      NA
      5      NA
      6      NA
      7      NA
      8      NA
      9      NA
      10     NA
      11     NA
      12     NA
      13     NA
      14     NA
      15     NA
      16     NA
      
      [[2]]
            fill x y PANEL group colour size angle hjust vjust alpha family fontface
      1  #0B775E 1 1     1     1  black    4     0   0.5   0.5    NA               1
      2  #57896B 2 1     1     5  black    4     0   0.5   0.5    NA               1
      3  #E6BE81 3 1     1     9  black    4     0   0.5   0.5    NA               1
      4  #E1BD6D 4 1     1    13  black    4     0   0.5   0.5    NA               1
      5  #57896B 1 2     1     2  black    4     0   0.5   0.5    NA               1
      6  #0B775E 2 2     1     6  black    4     0   0.5   0.5    NA               1
      7  #E7BE87 3 2     1    10  black    4     0   0.5   0.5    NA               1
      8  #E3BD77 4 2     1    14  black    4     0   0.5   0.5    NA               1
      9  #E6BE81 1 3     1     3  black    4     0   0.5   0.5    NA               1
      10 #E7BE87 2 3     1     7  black    4     0   0.5   0.5    NA               1
      11 #0B775E 3 3     1    11  black    4     0   0.5   0.5    NA               1
      12 #8E9C79 4 3     1    15  black    4     0   0.5   0.5    NA               1
      13 #E1BD6D 1 4     1     4  black    4     0   0.5   0.5    NA               1
      14 #E3BD77 2 4     1     8  black    4     0   0.5   0.5    NA               1
      15 #8E9C79 3 4     1    12  black    4     0   0.5   0.5    NA               1
      16 #0B775E 4 4     1    16  black    4     0   0.5   0.5    NA               1
         lineheight label
      1         1.2  1.00
      2         1.2  0.76
      3         1.2 -0.49
      4         1.2 -1.00
      5         1.2  0.76
      6         1.2  1.00
      7         1.2 -0.33
      8         1.2 -0.76
      9         1.2 -0.49
      10        1.2 -0.33
      11        1.2  1.00
      12        1.2  0.49
      13        1.2 -1.00
      14        1.2 -0.76
      15        1.2  0.49
      16        1.2  1.00
      
      [[3]]
           fill x y PANEL group shape colour size alpha stroke
      1 #DEBA91 3 2     1     2 cross  black   14    NA    0.5
      2 #DEBA91 2 3     1     1 cross  black   14    NA    0.5
      

---

    Code
      list(p$labels, pb$plot$plot_env$legend.title)
    Output
      [[1]]
      [[1]]$xlab
      NULL
      
      [[1]]$ylab
      NULL
      
      [[1]]$title
      NULL
      
      [[1]]$subtitle
      NULL
      
      [[1]]$caption
      atop(displaystyle(NULL), expr = paste(bold("X"), " = non-significant at ", 
          italic("p"), " < ", 0.01, " (Adjustment: ", "Hommel", ")"))
      
      [[1]]$fill
      [1] "value"
      
      [[1]]$x
      [1] "Var1"
      
      [[1]]$y
      [1] "Var2"
      
      
      [[2]]
      atop(atop(atop(scriptstyle(bold("sample sizes:")), italic(n)[min] ~ 
          "=" ~ "32"), atop(italic(n)[mode] ~ "=" ~ "32", italic(n)[max] ~ 
          "=" ~ "83")), atop(scriptstyle(bold("correlation:")), "Spearman"))
      

# checking Bayesian pearson (with NA)

    Code
      list(p$labels, pb$plot$plot_env$legend.title)
    Output
      [[1]]
      [[1]]$xlab
      NULL
      
      [[1]]$ylab
      NULL
      
      [[1]]$title
      NULL
      
      [[1]]$subtitle
      NULL
      
      [[1]]$caption
      NULL
      
      [[1]]$fill
      [1] "value"
      
      [[1]]$x
      [1] "Var1"
      
      [[1]]$y
      [1] "Var2"
      
      
      [[2]]
      atop(atop(atop(scriptstyle(bold("sample sizes:")), italic(n)[min] ~ 
          "=" ~ "56"), atop(italic(n)[mode] ~ "=" ~ "56", italic(n)[max] ~ 
          "=" ~ "56")), atop(scriptstyle(bold("correlation:")), "Bayesian Pearson"))
      

# checking all dataframe outputs

    Code
      suppressWarnings(purrr::pmap(.l = list(data = list(dplyr::select(ggplot2::msleep,
      brainwt, sleep_rem, bodywt)), type = list("p", "p", "np", "np", "r", "r", "bf",
        "bayes"), output = list("dataframe"), partial = list(TRUE, FALSE, TRUE, FALSE,
        TRUE, FALSE, TRUE, FALSE)), .f = ggcorrmat))
    Output
      [[1]]
      # A tibble: 3 x 11
        parameter1 parameter2 estimate conf.level conf.low conf.high statistic
        <chr>      <chr>         <dbl>      <dbl>    <dbl>     <dbl>     <dbl>
      1 brainwt    sleep_rem   -0.0961       0.95   -0.370     0.193    -0.655
      2 brainwt    bodywt       0.485        0.95    0.233     0.676     3.76 
      3 sleep_rem  bodywt      -0.108        0.95   -0.381     0.182    -0.737
        df.error p.value method              n.obs
           <int>   <dbl> <chr>               <int>
      1       46 0.929   Pearson correlation    48
      2       46 0.00144 Pearson correlation    48
      3       46 0.929   Pearson correlation    48
      
      [[2]]
      # A tibble: 3 x 11
        parameter1 parameter2 estimate conf.level conf.low conf.high statistic
        <chr>      <chr>         <dbl>      <dbl>    <dbl>     <dbl>     <dbl>
      1 brainwt    sleep_rem    -0.221       0.95   -0.476    0.0670     -1.54
      2 brainwt    bodywt        0.934       0.95    0.889    0.961      19.2 
      3 sleep_rem  bodywt       -0.328       0.95   -0.535   -0.0826     -2.66
        df.error  p.value method              n.obs
           <int>    <dbl> <chr>               <int>
      1       46 1.31e- 1 Pearson correlation    48
      2       54 2.75e-25 Pearson correlation    56
      3       59 1.99e- 2 Pearson correlation    61
      
      [[3]]
      # A tibble: 3 x 10
        parameter1 parameter2 estimate conf.level conf.low conf.high statistic
        <chr>      <chr>         <dbl>      <dbl>    <dbl>     <dbl>     <dbl>
      1 brainwt    sleep_rem    -0.271       0.95   -0.522    0.0230     23414
      2 brainwt    bodywt        0.785       0.95    0.640    0.876       3962
      3 sleep_rem  bodywt        0.154       0.95   -0.145    0.427      15588
         p.value method               n.obs
           <dbl> <chr>                <int>
      1 1.25e- 1 Spearman correlation    48
      2 1.20e-10 Spearman correlation    48
      3 2.96e- 1 Spearman correlation    48
      
      [[4]]
      # A tibble: 3 x 10
        parameter1 parameter2 estimate conf.level conf.low conf.high statistic
        <chr>      <chr>         <dbl>      <dbl>    <dbl>     <dbl>     <dbl>
      1 brainwt    sleep_rem    -0.414       0.95   -0.630    -0.139    26050.
      2 brainwt    bodywt        0.957       0.95    0.927     0.975     1254.
      3 sleep_rem  bodywt       -0.452       0.95   -0.636    -0.218    54904.
         p.value method               n.obs
           <dbl> <chr>                <int>
      1 3.45e- 3 Spearman correlation    48
      2 2.91e-30 Spearman correlation    56
      3 5.16e- 4 Spearman correlation    61
      
      [[5]]
      # A tibble: 3 x 11
        parameter1 parameter2 estimate conf.level conf.low conf.high statistic
        <chr>      <chr>         <dbl>      <dbl>    <dbl>     <dbl>     <dbl>
      1 brainwt    sleep_rem    -0.290       0.95   -0.531  -0.00694     -2.06
      2 brainwt    bodywt        0.681       0.95    0.493   0.809        6.32
      3 sleep_rem  bodywt        0.183       0.95   -0.107   0.444        1.26
        df.error     p.value method                         n.obs
           <int>       <dbl> <chr>                          <int>
      1       46 0.0904      Winsorized Pearson correlation    48
      2       46 0.000000292 Winsorized Pearson correlation    48
      3       46 0.213       Winsorized Pearson correlation    48
      
      [[6]]
      # A tibble: 3 x 11
        parameter1 parameter2 estimate conf.level conf.low conf.high statistic
        <chr>      <chr>         <dbl>      <dbl>    <dbl>     <dbl>     <dbl>
      1 brainwt    sleep_rem    -0.412       0.95   -0.623    -0.145     -3.06
      2 brainwt    bodywt        0.910       0.95    0.851     0.947     16.2 
      3 sleep_rem  bodywt       -0.375       0.95   -0.572    -0.136     -3.10
        df.error  p.value method                         n.obs
           <int>    <dbl> <chr>                          <int>
      1       46 5.86e- 3 Winsorized Pearson correlation    48
      2       54 7.22e-22 Winsorized Pearson correlation    56
      3       59 5.86e- 3 Winsorized Pearson correlation    61
      
      [[7]]
      # A tibble: 3 x 14
        parameter1 parameter2 estimate conf.level conf.low conf.high    pd
        <chr>      <chr>         <dbl>      <dbl>    <dbl>     <dbl> <dbl>
      1 brainwt    sleep_rem   -0.0911       0.95   -0.373     0.171 0.740
      2 brainwt    bodywt       0.461        0.95    0.228     0.663 1    
      3 sleep_rem  bodywt      -0.0959       0.95   -0.368     0.172 0.756
        rope.percentage prior.distribution prior.location prior.scale bayes.factor
                  <dbl> <chr>                       <dbl>       <dbl>        <dbl>
      1         0.430   beta                         1.41        1.41        0.269
      2         0.00425 beta                         1.41        1.41       73.6  
      3         0.434   beta                         1.41        1.41        0.283
        method                       n.obs
        <chr>                        <int>
      1 Bayesian Pearson correlation    48
      2 Bayesian Pearson correlation    48
      3 Bayesian Pearson correlation    48
      
      [[8]]
      # A tibble: 3 x 14
        parameter1 parameter2 estimate conf.level conf.low conf.high    pd
        <chr>      <chr>         <dbl>      <dbl>    <dbl>     <dbl> <dbl>
      1 brainwt    sleep_rem    -0.205       0.95   -0.458    0.0639 0.928
      2 brainwt    bodywt        0.926       0.95    0.883    0.960  1    
      3 sleep_rem  bodywt       -0.310       0.95   -0.537   -0.0972 0.990
        rope.percentage prior.distribution prior.location prior.scale bayes.factor
                  <dbl> <chr>                       <dbl>       <dbl>        <dbl>
      1          0.212  beta                         1.41        1.41     6.54e- 1
      2          0      beta                         1.41        1.41     1.58e+22
      3          0.0365 beta                         1.41        1.41     4.80e+ 0
        method                       n.obs
        <chr>                        <int>
      1 Bayesian Pearson correlation    48
      2 Bayesian Pearson correlation    56
      3 Bayesian Pearson correlation    61
      

# grouped_list works

    Code
      list(length(df1), length(df2), length(df5), length(df6))
    Output
      [[1]]
      [1] 4
      
      [[2]]
      [1] 3
      
      [[3]]
      [1] 4
      
      [[4]]
      [1] 11
      

---

    Code
      list(names(df1), names(df2), names(df5), names(df6))
    Output
      [[1]]
      [1] "carni"   "herbi"   "insecti" "omni"   
      
      [[2]]
      [1] "carni"   "insecti" "omni"   
      
      [[3]]
      [1] "carni"   "herbi"   "insecti" "omni"   
      
      [[4]]
       [1] "name"         "genus"        "vore"         "order"        "conservation"
       [6] "sleep_total"  "sleep_rem"    "sleep_cycle"  "awake"        "brainwt"     
      [11] "bodywt"      
      

# ggcoefstats works with data frames

    Code
      as.character(meta_info)
    Output
      [1] "list(log[e] * (BF[\"01\"]) == \"1.23\", widehat(delta)[\"difference\"]^\"posterior\" == \"0.13\", CI[\"95%\"]^HDI ~ \"[\" * \"-0.10\", \"0.44\" * \"]\", italic(\"r\")[\"Cauchy\"]^\"JZS\" == \"0.71\")"

# Getting help with ggstatsplot

Thanks for using ggstatsplot. Before filing an issue, there are a few places
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

Before opening a new issue, be sure to [search issues and pull requests](https://github.com/IndrajeetPatil/ggstatsplot/issues) to make sure the 
bug hasn't been reported and/or already fixed in the development version. By 
default, the search will be pre-populated with `is:issue is:open`. You can 
[edit the qualifiers](https://help.github.com/articles/searching-issues-and-pull-requests/) 
(e.g. `is:pr`, `is:closed`) as needed. For example, you'd simply
remove `is:open` to search _all_ issues in the repo, open or closed.


If you _are_ in the right place, and need to file an issue, please review the 
["File issues"](https://www.tidyverse.org/contribute/#issues) paragraph from 
the tidyverse contributing guidelines.

Thanks for your help!
# Contributing to ggstatsplot

This outlines how to propose a change to ggstatsplot. For more detailed
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
*  New code should follow the tidyverse [style guide](http://style.tidyverse.org).
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

Please note that the ggstatsplot project is released with a
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
---
name: Feature idea
about: Suggest an idea for this project

---

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**How could we do it?**
A description of actual ways of implementing a feature.
## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.1.1 (2021-08-10) |
|os       |Windows 10 x64               |
|system   |x86_64, mingw32              |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |English_United Kingdom.1252  |
|ctype    |English_United Kingdom.1252  |
|tz       |Europe/Berlin                |
|date     |2021-10-09                   |

# Dependencies

|package             |old        |new        |<U+0394>  |
|:-------------------|:----------|:----------|:--|
|ggstatsplot         |0.8.0      |0.9.0      |*  |
|BayesFactor         |0.9.12-4.2 |0.9.12-4.2 |   |
|bayestestR          |0.11.0     |0.11.0     |   |
|BWStest             |0.2.2      |0.2.2      |   |
|cachem              |1.0.6      |1.0.6      |   |
|cli                 |3.0.1      |3.0.1      |   |
|coda                |0.19-4     |0.19-4     |   |
|colorspace          |2.0-2      |2.0-2      |   |
|contfrac            |1.1-12     |1.1-12     |   |
|correlation         |0.7.1      |0.7.1      |   |
|cpp11               |0.4.0      |0.4.0      |   |
|crayon              |1.4.1      |1.4.1      |   |
|datawizard          |0.2.1      |0.2.1      |   |
|deSolve             |1.30       |1.30       |   |
|digest              |0.6.28     |0.6.28     |   |
|dplyr               |1.0.7      |1.0.7      |   |
|effectsize          |0.5        |0.5        |   |
|ellipsis            |0.3.2      |0.3.2      |   |
|elliptic            |1.4-0      |1.4-0      |   |
|fansi               |0.5.0      |0.5.0      |   |
|farver              |2.1.0      |2.1.0      |   |
|fastmap             |1.1.0      |1.1.0      |   |
|generics            |0.1.0      |0.1.0      |   |
|ggplot2             |3.3.5      |3.3.5      |   |
|ggrepel             |0.9.1      |0.9.1      |   |
|ggsignif            |0.6.3      |0.6.3      |   |
|glue                |1.4.2      |1.4.2      |   |
|gmp                 |0.6-2      |0.6-2      |   |
|gtable              |0.3.0      |0.3.0      |   |
|gtools              |3.9.2      |3.9.2      |   |
|hypergeo            |1.2-13     |1.2-13     |   |
|insight             |0.14.4     |0.14.4     |   |
|ipmisc              |6.0.2      |6.0.2      |   |
|isoband             |0.2.5      |0.2.5      |   |
|kSamples            |1.2-9      |1.2-9      |   |
|labeling            |0.4.2      |0.4.2      |   |
|lifecycle           |1.0.1      |1.0.1      |   |
|magrittr            |2.0.1      |2.0.1      |   |
|MatrixModels        |0.5-0      |0.5-0      |   |
|mc2d                |0.1-21     |0.1-21     |   |
|memoise             |2.0.0      |2.0.0      |   |
|multcompView        |0.1-8      |0.1-8      |   |
|munsell             |0.5.0      |0.5.0      |   |
|mvtnorm             |1.1-2      |1.1-2      |   |
|pairwiseComparisons |3.1.6      |3.2.0      |*  |
|paletteer           |1.4.0      |1.4.0      |   |
|parameters          |0.14.0     |0.14.0.3   |*  |
|patchwork           |1.1.1      |1.1.1      |   |
|pbapply             |1.5-0      |1.5-0      |   |
|performance         |0.8.0      |0.8.0      |   |
|pillar              |1.6.3      |1.6.3      |   |
|pkgconfig           |2.0.3      |2.0.3      |   |
|plyr                |1.8.6      |1.8.6      |   |
|PMCMRplus           |1.9.0      |1.9.0      |   |
|prismatic           |1.0.0      |1.0.0      |   |
|purrr               |0.3.4      |0.3.4      |   |
|R6                  |2.5.1      |2.5.1      |   |
|RColorBrewer        |1.1-2      |1.1-2      |   |
|Rcpp                |1.0.7      |1.0.7      |   |
|RcppEigen           |0.3.3.9.1  |0.3.3.9.1  |   |
|rematch2            |2.1.2      |2.1.2      |   |
|reshape             |0.8.8      |0.8.8      |   |
|rlang               |0.4.11     |0.4.11     |   |
|Rmpfr               |0.8-5      |0.8-5      |   |
|rstudioapi          |0.13       |0.13       |   |
|scales              |1.1.1      |1.1.1      |   |
|statsExpressions    |1.1.0      |1.2.0      |*  |
|stringi             |1.7.5      |1.7.5      |   |
|stringr             |1.4.0      |1.4.0      |   |
|SuppDists           |1.1-9.5    |1.1-9.5    |   |
|tibble              |3.1.5      |3.1.5      |   |
|tidyr               |1.1.4      |1.1.4      |   |
|tidyselect          |1.1.1      |1.1.1      |   |
|utf8                |1.2.2      |1.2.2      |   |
|vctrs               |0.3.8      |0.3.8      |   |
|viridisLite         |0.4.0      |0.4.0      |   |
|withr               |2.4.2      |2.4.2      |   |
|WRS2                |1.1-3      |1.1-3      |   |
|zeallot             |0.1.0      |0.1.0      |   |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
output:  
   github_document: 
     pandoc_args: --webtex=https://chart.apis.google.com/chart?cht=tx&chl= 
---

  <!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo=FALSE}
## show me all columns
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300, ## change to 300 once on CRAN
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "##>",
  fig.path = "man/figures/README-"
)
```

## `{ggstatsplot}`: `{ggplot2}` Based Plots with Statistical Details 

Status | Usage| Miscellaneous
----------------- | ----------------- | ----------------- 
  [![R build status](https://github.com/IndrajeetPatil/ggstatsplot/workflows/R-CMD-check/badge.svg)](https://github.com/IndrajeetPatil/ggstatsplot) | [![Total downloads badge](https://cranlogs.r-pkg.org/badges/grand-total/ggstatsplot?color=blue)](https://CRAN.R-project.org/package=ggstatsplot) | [![Codecov](https://codecov.io/gh/IndrajeetPatil/ggstatsplot/branch/master/graph/badge.svg)](https://app.codecov.io/gh/IndrajeetPatil/ggstatsplot?branch=master)
[![lints](https://github.com/IndrajeetPatil/ggstatsplot/workflows/lint/badge.svg)](https://github.com/IndrajeetPatil/ggstatsplot) | [![Daily downloads badge](https://cranlogs.r-pkg.org/badges/last-day/ggstatsplot?color=blue)](https://CRAN.R-project.org/package=ggstatsplot) | [![status](https://tinyverse.netlify.com/badge/ggstatsplot)](https://CRAN.R-project.org/package=ggstatsplot)
[![pkgdown](https://github.com/IndrajeetPatil/ggstatsplot/workflows/pkgdown/badge.svg)](https://github.com/IndrajeetPatil/ggstatsplot/actions) | [![DOI](https://joss.theoj.org/papers/10.21105/joss.03167/status.svg)](https://doi.org/10.21105/joss.03167) | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)


## Raison d'être <img src="man/figures/logo.png" align="right" width="360" />

> "What is to be sought in designs for the display of information is the clear
portrayal of complexity. Not the complication of the simple; rather ... the
revelation of the complex."
- Edward R. Tufte

[`{ggstatsplot}`](https://indrajeetpatil.github.io/ggstatsplot/) is an extension
of [`{ggplot2}`](https://github.com/tidyverse/ggplot2) package for creating
graphics with details from statistical tests included in the information-rich
plots themselves. In a typical exploratory data analysis workflow, data
visualization and statistical modeling are two different phases: visualization
informs modeling, and modeling in its turn can suggest a different visualization
method, and so on and so forth. The central idea of `{ggstatsplot}` is simple:
combine these two phases into one in the form of graphics with statistical
details, which makes data exploration simpler and faster.

## Installation

Type | Source | Command
---|---|---
Release | [![CRAN Status](https://www.r-pkg.org/badges/version/ggstatsplot)](https://cran.r-project.org/package=ggstatsplot) | `install.packages("ggstatsplot")`
Development | [![Project Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/##active) | `remotes::install_github("IndrajeetPatil/ggstatsplot")`

Linux users may encounter some installation problems. In particular, the
`{ggstatsplot}` package depends on the `{PMCMRplus}` package.

```
ERROR: dependencies ‘gmp’, ‘Rmpfr’ are not available for package ‘PMCMRplus’
```

This means that your operating system lacks `gmp` and `Rmpfr` libraries.

If you use `Ubuntu`, you can install these dependencies:

```
sudo apt-get install libgmp3-dev
sudo apt-get install libmpfr-dev
```

The following `README` file briefly describes the installation procedure:
<https://CRAN.R-project.org/package=PMCMRplus/readme/README.html>

## Citation

If you want to cite this package in a scientific journal or in any other
context, run the following code in your `R` console:

```{r citation, comment=""}
citation("ggstatsplot")
```

There is currently a publication in preparation corresponding to this package
and the citation will be updated once it's published.

## Documentation and Examples

To see the detailed documentation for each function in the stable **CRAN**
version of the package, see:

  - Website: <https://indrajeetpatil.github.io/ggstatsplot/>
    
  - Presentation:
    <https://indrajeetpatil.github.io/ggstatsplot_slides/slides/ggstatsplot_presentation.html##1>
    
  - Vignettes: <https://indrajeetpatil.github.io/ggstatsplot/articles/>

## Summary of available plots

It, therefore, produces a limited kinds of plots for the supported analyses:

Function | Plot | Description | Lifecycle
------- | ---------- | ------------------------- | ----
`ggbetweenstats` | **violin plots** | for comparisons *between* groups/conditions | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
`ggwithinstats` | **violin plots** | for comparisons *within* groups/conditions | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
`gghistostats` | **histograms** | for distribution about numeric variable | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
`ggdotplotstats` | **dot plots/charts** | for distribution about labeled numeric variable | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
`ggscatterstats` | **scatterplots** | for correlation between two variables | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
`ggcorrmat` | **correlation matrices** | for correlations between multiple variables | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
`ggpiestats` | **pie charts** | for categorical data | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
`ggbarstats` | **bar charts** | for categorical data | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
`ggcoefstats` | **dot-and-whisker plots** | for regression models and meta-analysis | [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

In addition to these basic plots, `{ggstatsplot}` also provides **`grouped_`**
versions (see below) that makes it easy to repeat the same analysis for
any grouping variable.

## Summary of types of statistical analyses

The table below summarizes all the different types of analyses currently
supported in this package-

Functions | Description | Parametric | Non-parametric | Robust | Bayesian
------- | ------------------ | ---- | ----- | ----| ----- 
`ggbetweenstats` | Between group/condition comparisons | ✅ | ✅ | ✅ | ✅
`ggwithinstats` | Within group/condition comparisons | ✅ | ✅ | ✅ | ✅
`gghistostats`, `ggdotplotstats` | Distribution of a numeric variable | ✅ | ✅ | ✅ | ✅
`ggcorrmat` | Correlation matrix | ✅ | ✅ | ✅ | ✅
`ggscatterstats` | Correlation between two variables | ✅ | ✅ | ✅ | ✅
`ggpiestats`, `ggbarstats` | Association between categorical variables | ✅ | ✅ | ❌ | ✅
`ggpiestats`, `ggbarstats` | Equal proportions for categorical variable levels | ✅ | ✅ | ❌ | ✅
`ggcoefstats` | Regression model coefficients | ✅ | ✅ | ✅ | ✅
`ggcoefstats` | Random-effects meta-analysis | ✅ | ❌ | ✅ | ✅

Summary of Bayesian analysis

Analysis | Hypothesis testing | Estimation
------------------ | ---------- | ---------
(one/two-sample) t-test | ✅ | ✅
one-way ANOVA | ✅ |✅
correlation | ✅ | ✅
(one/two-way) contingency table | ✅ | ✅
random-effects meta-analysis | ✅ | ✅

## Statistical reporting

For **all** statistical tests reported in the plots, the default template abides
by the gold standard for statistical reporting. For example, here are results
from Yuen's test for trimmed means (robust *t*-test):

<img src="man/figures/stats_reporting_format.png" align="center" />

## Summary of statistical tests and effect sizes

Statistical analysis is carried out by `{statsExpressions}` package, and thus
a summary table of all the statistical tests currently supported across
various functions can be found in article for that package:
<https://indrajeetpatil.github.io/statsExpressions/articles/stats_details.html>

## Primary functions

Here are examples of the main functions currently supported in `{ggstatsplot}`.

**Note**: If you are reading this on `GitHub` repository, the documentation below
is for the **development** version of the package. So you may see some features
available here that are not currently present in the stable version of this
package on **CRAN**. For documentation relevant for the `CRAN` version, see:
<https://CRAN.R-project.org/package=ggstatsplot/readme/README.html>

### `ggbetweenstats`

This function creates either a violin plot, a box plot, or a mix of two for
**between**-group or **between**-condition comparisons with results from
statistical tests in the subtitle. The simplest function call looks like this-

```{r ggbetweenstats1}
## for reproducibility
set.seed(123)
library(ggstatsplot)

## plot
ggbetweenstats(
  data  = iris,
  x     = Species,
  y     = Sepal.Length,
  title = "Distribution of sepal length across Iris species"
)
```

**Defaults** return<br>

✅ raw data + distributions <br>
✅ descriptive statistics <br>
✅ inferential statistics <br>
✅ effect size + CIs <br>
✅ pairwise comparisons <br>
✅ Bayesian hypothesis-testing <br>
✅ Bayesian estimation <br>

A number of other arguments can be specified to make this plot even more
informative or change some of the default options. Additionally, there is also a
`grouped_` variant of this function that makes it easy to repeat the same
operation across a **single** grouping variable:

```{r ggbetweenstats2, fig.height=8, fig.width=12}
## for reproducibility
set.seed(123)

## plot
grouped_ggbetweenstats(
  data             = dplyr::filter(movies_long, genre %in% c("Action", "Comedy")),
  x                = mpaa,
  y                = length,
  grouping.var     = genre, ## grouping variable
  outlier.tagging  = TRUE, ## whether outliers need to be tagged
  outlier.label    = title, ## variable to be used for tagging outliers
  outlier.coef     = 2,
  ggsignif.args    = list(textsize = 4, tip_length = 0.01),
  p.adjust.method  = "bonferroni", ## method for adjusting p-values for multiple comparisons
  ## adding new components to `{ggstatsplot}` default
  ggplot.component = list(ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis())),
  caption          = "Source: IMDb (Internet Movie Database)",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  annotation.args  = list(title = "Differences in movie length by mpaa ratings for different genres")
)
```

Note here that the function can be used to tag outliers!

##### Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
raw data | `ggplot2::geom_point` | `point.args`
box plot | `ggplot2::geom_boxplot` | ❌
density plot | `ggplot2::geom_violin` | `violin.args`
centrality measure point | `ggplot2::geom_point` | `centrality.point.args`
centrality measure label | `ggrepel::geom_label_repel` | `centrality.label.args`
outlier point | `ggplot2::stat_boxplot` | ❌
outlier label | `ggrepel::geom_label_repel` | `outlier.label.args`
pairwise comparisons | `ggsignif::geom_signif` | `ggsignif.args`

##### Summary of tests

**Central tendency measure**

Type | Measure | Function used
----------- | --------- | ------------------ 
Parametric | mean | `parameters::describe_distribution`
Non-parametric | median | `parameters::describe_distribution`
Robust | trimmed mean | `parameters::describe_distribution`
Bayesian | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution`

**Hypothesis testing**

Type | No. of groups | Test | Function used
----------- | --- | ------------------------- | -----
Parametric | > 2 | Fisher's or Welch's one-way ANOVA | `stats::oneway.test`
Non-parametric | > 2 | Kruskal–Wallis one-way ANOVA | `stats::kruskal.test`
Robust | > 2 | Heteroscedastic one-way ANOVA for trimmed means | `WRS2::t1way`
Bayes Factor | > 2 | Fisher's ANOVA | `BayesFactor::anovaBF`
Parametric | 2 | Student's or Welch's *t*-test | `stats::t.test`
Non-parametric | 2 | Mann–Whitney *U* test | `stats::wilcox.test`
Robust | 2 |  Yuen's test for trimmed means | `WRS2::yuen`
Bayesian | 2 | Student's *t*-test | `BayesFactor::ttestBF`

**Effect size estimation**

Type | No. of groups | Effect size | CI? | Function used
----------- | --- | ------------------------- | --- | -----
Parametric | > 2 | $\eta_{p}^2$, $\omega_{p}^2$ | ✅ | `effectsize::omega_squared`, `effectsize::eta_squared`
Non-parametric | > 2 | $\epsilon_{ordinal}^2$ | ✅ | `effectsize::rank_epsilon_squared`
Robust | > 2 | $\xi$ (Explanatory measure of effect size) | ✅ | `WRS2::t1way`
Bayes Factor | > 2 | $R_{posterior}^2$ | ✅ | `performance::r2_bayes`
Parametric | 2 | Cohen's *d*, Hedge's *g* | ✅ | `effectsize::cohens_d`, `effectsize::hedges_g`
Non-parametric | 2 | *r* (rank-biserial correlation) | ✅ | `effectsize::rank_biserial`
Robust | 2 |  $\xi$ (Explanatory measure of effect size) | ✅ | `WRS2::yuen.effect.ci`
Bayesian | 2 | $\delta_{posterior}$ | ✅ | `bayestestR::describe_posterior`

**Pairwise comparison tests**

Type | Equal variance? | Test | *p*-value adjustment? | Function used
----------- | --- | ------------------------- | --- | -----
Parametric | No | Games-Howell test | ✅ | `PMCMRplus::gamesHowellTest`
Parametric | Yes | Student's *t*-test | ✅ | `stats::pairwise.t.test`
Non-parametric | No | Dunn test | ✅ | `PMCMRplus::kwAllPairsDunnTest`
Robust | No | Yuen's trimmed means test | ✅ | `WRS2::lincon`
Bayesian | `NA` | Student's *t*-test | `NA` | `BayesFactor::ttestBF`

For more, see the `ggbetweenstats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggbetweenstats.html>

### `ggwithinstats`

`ggbetweenstats` function has an identical twin function `ggwithinstats` for
repeated measures designs that behaves in the same fashion with a few minor
tweaks introduced to properly visualize the repeated measures design. As can be
seen from an example below, the only difference between the plot structure is
that now the group means are connected by paths to highlight the fact that these
data are paired with each other.

```{r ggwithinstats1, fig.width=8, fig.height=6}
## for reproducibility and data
set.seed(123)
library(WRS2) ## for data
library(afex) ## to run anova

## plot
ggwithinstats(
  data    = WineTasting,
  x       = Wine,
  y       = Taste,
  title   = "Wine tasting",
  caption = "Data source: `WRS2` R package",
  ggtheme = ggthemes::theme_fivethirtyeight()
)
```

**Defaults** return<br>

✅ raw data + distributions <br>
✅ descriptive statistics <br>
✅ inferential statistics <br>
✅ effect size + CIs <br>
✅ pairwise comparisons <br>
✅ Bayesian hypothesis-testing <br>
✅ Bayesian estimation <br>

The central tendency measure displayed will depend on the statistics:

Type | Measure | Function used
----------- | --------- | ------------------ 
Parametric | mean | `parameters::describe_distribution`
Non-parametric | median | `parameters::describe_distribution`
Robust | trimmed mean | `parameters::describe_distribution`
Bayesian | MAP estimate | `parameters::describe_distribution`

As with the `ggbetweenstats`, this function also has a `grouped_` variant that
makes repeating the same analysis across a single grouping variable quicker. We
will see an example with only repeated measurements-

```{r ggwithinstats2, fig.height=6, fig.width=14}
## common setup
set.seed(123)

## plot
grouped_ggwithinstats(
  data            = dplyr::filter(bugs_long, region %in% c("Europe", "North America"), condition %in% c("LDLF", "LDHF")),
  x               = condition,
  y               = desire,
  type            = "np", ## non-parametric statistics
  xlab            = "Condition",
  ylab            = "Desire to kill an artrhopod",
  grouping.var    = region,
  outlier.tagging = TRUE,
  outlier.label   = education
)
```

##### Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
raw data | `ggplot2::geom_point` | `point.args`
point path | `ggplot2::geom_path` | `point.path.args`
box plot | `ggplot2::geom_boxplot` | `boxplot.args`
density plot | `ggplot2::geom_violin` | `violin.args`
centrality measure point | `ggplot2::geom_point` | `centrality.point.args`
centrality measure point path | `ggplot2::geom_path` | `centrality.path.args`
centrality measure label | `ggrepel::geom_label_repel` | `centrality.label.args`
outlier point | `ggplot2::stat_boxplot` | ❌
outlier label | `ggrepel::geom_label_repel` | `outlier.label.args`
pairwise comparisons | `ggsignif::geom_signif` | `ggsignif.args`

##### Summary of tests

**Central tendency measure**

Type | Measure | Function used
----------- | --------- | ------------------ 
Parametric | mean | `parameters::describe_distribution`
Non-parametric | median | `parameters::describe_distribution`
Robust | trimmed mean | `parameters::describe_distribution`
Bayesian | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution`

**Hypothesis testing**

Type | No. of groups | Test | Function used
----------- | --- | ------------------------- | -----
Parametric | > 2 | One-way repeated measures ANOVA | `afex::aov_ez`
Non-parametric | > 2 | Friedman rank sum test | `stats::friedman.test`
Robust | > 2 | Heteroscedastic one-way repeated measures ANOVA for trimmed means | `WRS2::rmanova`
Bayes Factor | > 2 | One-way repeated measures ANOVA | `BayesFactor::anovaBF`
Parametric | 2 | Student's *t*-test | `stats::t.test`
Non-parametric | 2 | Wilcoxon signed-rank test | `stats::wilcox.test`
Robust | 2 | Yuen's test on trimmed means for dependent samples | `WRS2::yuend`
Bayesian | 2 | Student's *t*-test | `BayesFactor::ttestBF`

**Effect size estimation**

Type | No. of groups | Effect size | CI? | Function used
----------- | --- | ------------------------- | --- | -----
Parametric | > 2 | $\eta_{p}^2$, $\omega_{p}^2$ | ✅ | `effectsize::omega_squared`, `effectsize::eta_squared`
Non-parametric | > 2 | $W_{Kendall}$ (Kendall's coefficient of concordance) | ✅ | `effectsize::kendalls_w`
Robust | > 2 | $\delta_{R-avg}^{AKP}$ (Algina-Keselman-Penfield robust standardized difference average) | ✅ | `WRS2::wmcpAKP`
Bayes Factor | > 2 | $R_{Bayesian}^2$ | ✅ | `performance::r2_bayes`
Parametric | 2 | Cohen's *d*, Hedge's *g* | ✅ | `effectsize::cohens_d`, `effectsize::hedges_g`
Non-parametric | 2 | *r* (rank-biserial correlation) | ✅ | `effectsize::rank_biserial`
Robust | 2 |  $\delta_{R}^{AKP}$ (Algina-Keselman-Penfield robust standardized difference) | ✅ | `WRS2::wmcpAKP`
Bayesian | 2 | $\delta_{posterior}$ | ✅ | `bayestestR::describe_posterior`

**Pairwise comparison tests**

Type | Test | *p*-value adjustment? | Function used
----------- | ---------------------------- | --- | -----
Parametric | Student's *t*-test | ✅ | `stats::pairwise.t.test`
Non-parametric | Durbin-Conover test | ✅ | `PMCMRplus::durbinAllPairsTest` 
Robust | Yuen's trimmed means test | ✅ | `WRS2::rmmcp`
Bayesian | Student's *t*-test | ❌ | `BayesFactor::ttestBF`

For more, see the `ggwithinstats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggwithinstats.html>

### `gghistostats`

To visualize the distribution of a single variable and check if its mean is
significantly different from a specified value with a one-sample test,
`gghistostats` can be used.

```{r gghistostats1, fig.width=8}
## for reproducibility
set.seed(123)

## plot
gghistostats(
  data       = ggplot2::msleep, ## dataframe from which variable is to be taken
  x          = awake, ## numeric variable whose distribution is of interest
  title      = "Amount of time spent awake", ## title for the plot
  caption    = "Source: Mammalian sleep data set",
  test.value = 12, ## default value is 0
  binwidth   = 1, ## binwidth value (experiment)
  ggtheme    = hrbrthemes::theme_ipsum_tw()
)
```

**Defaults** return<br>

✅ counts + proportion for bins<br>
✅ descriptive statistics <br>
✅ inferential statistics <br>
✅ effect size + CIs <br>
✅ Bayesian hypothesis-testing <br>
✅ Bayesian estimation <br>

There is also a `grouped_` variant of this function that makes it
easy to repeat the same operation across a **single** grouping variable:

```{r gghistostats2, fig.height=6, fig.width=12}
## for reproducibility
set.seed(123)

## plot
grouped_gghistostats(
  data              = dplyr::filter(movies_long, genre %in% c("Action", "Comedy")),
  x                 = budget,
  test.value        = 50,
  type              = "nonparametric",
  xlab              = "Movies budget (in million US$)",
  grouping.var      = genre, ## grouping variable
  normal.curve      = TRUE, ## superimpose a normal distribution curve
  normal.curve.args = list(color = "red", size = 1),
  ggtheme           = ggthemes::theme_tufte(),
  ## modify the defaults from `{ggstatsplot}` for each plot
  ggplot.component  = ggplot2::labs(caption = "Source: IMDB.com"),
  plotgrid.args     = list(nrow = 1),
  annotation.args   = list(title = "Movies budgets for different genres")
)
```

##### Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
histogram bin | `ggplot2::stat_bin` | `bin.args`
centrality measure line | `ggplot2::geom_vline` | `centrality.line.args`
normality curve | `ggplot2::stat_function` | `normal.curve.args`

##### Summary of tests

**Central tendency measure**

Type | Measure | Function used
----------- | --------- | ------------------ 
Parametric | mean | `parameters::describe_distribution`
Non-parametric | median | `parameters::describe_distribution`
Robust | trimmed mean | `parameters::describe_distribution`
Bayesian | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution`

**Hypothesis testing**

Type | Test | Function used
------------------ | ------------------------- | -----
Parametric | One-sample Student's *t*-test | `stats::t.test`
Non-parametric | One-sample Wilcoxon test | `stats::wilcox.test`
Robust | Bootstrap-*t* method for one-sample test | `WRS2::trimcibt`
Bayesian | One-sample Student's *t*-test | `BayesFactor::ttestBF`

**Effect size estimation**

Type | Effect size | CI? | Function used
------------ | ----------------------- | --- | -----
Parametric | Cohen's *d*, Hedge's *g* | ✅ | `effectsize::cohens_d`, `effectsize::hedges_g`
Non-parametric | *r* (rank-biserial correlation) | ✅ | `effectsize::rank_biserial`
Robust | trimmed mean | ✅ | `WRS2::trimcibt`
Bayes Factor | $\delta_{posterior}$ | ✅ | `bayestestR::describe_posterior`

For more, including information about the variant of this function
`grouped_gghistostats`, see the `gghistostats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/gghistostats.html>

### `ggdotplotstats`

This function is similar to `gghistostats`, but is intended to be used when the
numeric variable also has a label.

```{r ggdotplotstats1, fig.height=10, fig.width=8}
## for reproducibility
set.seed(123)

## plot
ggdotplotstats(
  data       = dplyr::filter(gapminder::gapminder, continent == "Asia"),
  y          = country,
  x          = lifeExp,
  test.value = 55,
  type       = "robust",
  title      = "Distribution of life expectancy in Asian continent",
  xlab       = "Life expectancy",
  caption    = "Source: Gapminder dataset from https://www.gapminder.org/"
)
```

**Defaults** return<br>

✅ descriptives (mean + sample size) <br>
✅ inferential statistics <br>
✅ effect size + CIs <br>
✅ Bayesian hypothesis-testing <br>
✅ Bayesian estimation <br>

As with the rest of the functions in this package, there is also a `grouped_`
variant of this function to facilitate looping the same operation for all levels
of a single grouping variable.

```{r ggdotplotstats2, fig.height=6, fig.width=12}
## for reproducibility
set.seed(123)

## plot
grouped_ggdotplotstats(
  data            = dplyr::filter(ggplot2::mpg, cyl %in% c("4", "6")),
  x               = cty,
  y               = manufacturer,
  type            = "bayes", ## Bayesian test
  xlab            = "city miles per gallon",
  ylab            = "car manufacturer",
  grouping.var    = cyl, ## grouping variable
  test.value      = 15.5,
  point.args      = list(color = "red", size = 5, shape = 13),
  annotation.args = list(title = "Fuel economy data")
)
```

##### Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
raw data | `ggplot2::geom_point` | `point.args`
centrality measure line | `ggplot2::geom_vline` | `centrality.line.args`

##### Summary of tests

**Central tendency measure**

Type | Measure | Function used
----------- | --------- | ------------------ 
Parametric | mean | `parameters::describe_distribution`
Non-parametric | median | `parameters::describe_distribution`
Robust | trimmed mean | `parameters::describe_distribution`
Bayesian | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution`

**Hypothesis testing**

Type | Test | Function used
------------------ | ------------------------- | -----
Parametric | One-sample Student's *t*-test | `stats::t.test`
Non-parametric | One-sample Wilcoxon test | `stats::wilcox.test`
Robust | Bootstrap-*t* method for one-sample test | `WRS2::trimcibt`
Bayesian | One-sample Student's *t*-test | `BayesFactor::ttestBF`

**Effect size estimation**

Type | Effect size | CI? | Function used
------------ | ----------------------- | --- | -----
Parametric | Cohen's *d*, Hedge's *g* | ✅ | `effectsize::cohens_d`, `effectsize::hedges_g`
Non-parametric | *r* (rank-biserial correlation) | ✅ | `effectsize::rank_biserial`
Robust | trimmed mean | ✅ | `WRS2::trimcibt`
Bayes Factor | $\delta_{posterior}$ | ✅ | `bayestestR::describe_posterior`

### `ggscatterstats`

This function creates a scatterplot with marginal distributions overlaid on the
axes and results from statistical tests in the subtitle:

```{r ggscatterstats1, fig.height=6}
ggscatterstats(
  data  = ggplot2::msleep,
  x     = sleep_rem,
  y     = awake,
  xlab  = "REM sleep (in hours)",
  ylab  = "Amount of time spent awake (in hours)",
  title = "Understanding mammalian sleep"
)
```

**Defaults** return<br>

✅ raw data + distributions <br>
✅ marginal distributions <br>
✅ inferential statistics <br>
✅ effect size + CIs <br>
✅ Bayesian hypothesis-testing <br>
✅ Bayesian estimation <br>

There is also a `grouped_` variant of this function that makes it
easy to repeat the same operation across a **single** grouping variable.

```{r ggscatterstats2, fig.height=8, fig.width=14}
## for reproducibility
set.seed(123)

## plot
grouped_ggscatterstats(
  data             = dplyr::filter(movies_long, genre %in% c("Action", "Comedy")),
  x                = rating,
  y                = length,
  grouping.var     = genre, ## grouping variable
  label.var        = title,
  label.expression = length > 200,
  xlab             = "IMDB rating",
  ggtheme          = ggplot2::theme_grey(),
  ggplot.component = list(ggplot2::scale_x_continuous(breaks = seq(2, 9, 1), limits = (c(2, 9)))),
  plotgrid.args    = list(nrow = 1),
  annotation.args  = list(title = "Relationship between movie length and IMDB ratings")
)
```

##### Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
raw data | `ggplot2::geom_point` | `point.args`
labels for raw data | `ggrepel::geom_label_repel` | `point.label.args`
smooth line | `ggplot2::geom_smooth` | `smooth.line.args`
marginal histograms | `ggside::geom_xsidehistogram`,  `ggside::geom_ysidehistogram` | `xsidehistogram.args`, `ysidehistogram.args`

##### Summary of tests

**Hypothesis testing** and **Effect size estimation**

Type | Test | CI? | Function used
----------- | ------------------------- | --- | -----
Parametric | Pearson's correlation coefficient | ✅ | `correlation::correlation`
Non-parametric | Spearman's rank correlation coefficient | ✅ | `correlation::correlation`
Robust | Winsorized Pearson correlation coefficient | ✅ | `correlation::correlation`
Bayesian | Pearson's correlation coefficient | ✅ | `correlation::correlation`

For more, see the `ggscatterstats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggscatterstats.html>

### `ggcorrmat`

`ggcorrmat` makes a correlalogram (a matrix of correlation coefficients) with
minimal amount of code. Just sticking to the defaults itself produces
publication-ready correlation matrices. But, for the sake of exploring the
available options, let's change some of the defaults. For example, multiple
aesthetics-related arguments can be modified to change the appearance of the
correlation matrix.

```{r ggcorrmat1}
## for reproducibility
set.seed(123)

## as a default this function outputs a correlation matrix plot
ggcorrmat(
  data     = ggplot2::msleep,
  colors   = c("#B2182B", "white", "#4D4D4D"),
  title    = "Correlalogram for mammals sleep dataset",
  subtitle = "sleep units: hours; weight units: kilograms"
)
```

**Defaults** return<br>

✅ effect size + significance<br>
✅ careful handling of `NA`s

If there are `NA`s present in the selected variables, the legend will display
minimum, median, and maximum number of pairs used for correlation tests.

There is also a `grouped_` variant of this function that makes it
easy to repeat the same operation across a **single** grouping variable:

```{r ggcorrmat2, fig.height=6, fig.width=10}
## for reproducibility
set.seed(123)

## plot
grouped_ggcorrmat(
  data         = dplyr::filter(movies_long, genre %in% c("Action", "Comedy")),
  type         = "robust", ## correlation method
  colors       = c("#cbac43", "white", "#550000"),
  grouping.var = genre, ## grouping variable
  matrix.type  = "lower" ## type of matrix
)
```

You can also get a dataframe containing all relevant details from the
statistical tests:

```{r ggcorrmat3}
## setup
set.seed(123)

## tidy data as output
ggcorrmat(
  data   = dplyr::select(ggplot2::msleep, dplyr::matches("sleep|awake")),
  type   = "bayes",
  output = "dataframe"
)
```

Additionally, **partial** correlation are also supported:

```{r ggcorrmat4}
## setup
set.seed(123)

## tidy data as output
ggcorrmat(
  data    = dplyr::select(ggplot2::msleep, dplyr::matches("sleep|awake")),
  type    = "bayes",
  partial = TRUE,
  output  = "dataframe"
)
```

##### Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
correlation matrix | `ggcorrplot::ggcorrplot` | `ggcorrplot.args`

##### Summary of tests

**Hypothesis testing** and **Effect size estimation**

Type | Test | CI? | Function used
----------- | ------------------------- | --- | -----
Parametric | Pearson's correlation coefficient | ✅ | `correlation::correlation`
Non-parametric | Spearman's rank correlation coefficient | ✅ | `correlation::correlation`
Robust | Winsorized Pearson correlation coefficient | ✅ | `correlation::correlation`
Bayesian | Pearson's correlation coefficient | ✅ | `correlation::correlation`

For examples and more information, see the `ggcorrmat` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcorrmat.html>

### `ggpiestats`

This function creates a pie chart for categorical or nominal variables with
results from contingency table analysis (Pearson's chi-squared test for
between-subjects design and McNemar's chi-squared test for within-subjects
design) included in the subtitle of the plot. If only one categorical variable
is entered, results from one-sample proportion test (i.e., a chi-squared
goodness of fit test) will be displayed as a subtitle.

To study an interaction between two categorical variables:

```{r ggpiestats1, fig.height=4, fig.width=8}
## for reproducibility
set.seed(123)

## plot
ggpiestats(
  data         = mtcars,
  x            = am,
  y            = cyl,
  package      = "wesanderson",
  palette      = "Royal1",
  title        = "Dataset: Motor Trend Car Road Tests", ## title for the plot
  legend.title = "Transmission", ## title for the legend
  caption      = "Source: 1974 Motor Trend US magazine"
)
```

**Defaults** return<br>

✅ descriptives (frequency + %s) <br>
✅ inferential statistics <br>
✅ effect size + CIs <br>
✅ Goodness-of-fit tests <br>
✅ Bayesian hypothesis-testing <br>
✅ Bayesian estimation <br>

There is also a `grouped_` variant of this function that makes it
easy to repeat the same operation across a **single** grouping variable.
Following example is a case where the theoretical question is about proportions
for different levels of a single nominal variable:

```{r ggpiestats2, fig.height=6, fig.width=10}
## for reproducibility
set.seed(123)

## plot
grouped_ggpiestats(
  data         = mtcars,
  x            = cyl,
  grouping.var = am, ## grouping variable
  label.repel  = TRUE, ## repel labels (helpful for overlapping labels)
  package      = "ggsci", ## package from which color palette is to be taken
  palette      = "default_ucscgb" ## choosing a different color palette
)
```

##### Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
pie slices | `ggplot2::geom_col` | ❌
descriptive labels | `ggplot2::geom_label`/`ggrepel::geom_label_repel` | `label.args`

##### Summary of tests

**two-way table**

**Hypothesis testing**

Type | Design | Test | Function used
----------- | ----- | ------------------------- | -----
Parametric/Non-parametric | Unpaired | Pearson's $\chi^2$ test | `stats::chisq.test`
Bayesian | Unpaired | Bayesian Pearson's $\chi^2$ test | `BayesFactor::contingencyTableBF`
Parametric/Non-parametric | Paired  | McNemar's $\chi^2$ test | `stats::mcnemar.test`
Bayesian | Paired  | ❌ | ❌

**Effect size estimation**

Type | Design | Effect size | CI? | Function used
----------- | ----- | ------------------------- | --- | -----
Parametric/Non-parametric | Unpaired | Cramer's $V$ | ✅ | `effectsize::cramers_v`
Bayesian | Unpaired | Cramer's $V$ | ✅ | `effectsize::cramers_v`
Parametric/Non-parametric | Paired | Cohen's $g$ | ✅ | `effectsize::cohens_g`
Bayesian | Paired | ❌ | ❌ | ❌

**one-way table**

**Hypothesis testing**

Type | Test | Function used
----------- | ------------------------- | -----
Parametric/Non-parametric | Goodness of fit $\chi^2$ test | `stats::chisq.test`
Bayesian | Bayesian Goodness of fit $\chi^2$ test | (custom)

**Effect size estimation**

Type | Effect size | CI? | Function used
----------- | ------------------------- | --- | -----
Parametric/Non-parametric | Pearson's $C$ | ✅ | `effectsize::pearsons_c`
Bayesian | ❌ | ❌ | ❌

For more, see the `ggpiestats` vignette:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggpiestats.html>

### `ggbarstats`

In case you are not a fan of pie charts (for very good reasons), you can
alternatively use `ggbarstats` function which has a similar syntax.

N.B. The *p*-values from one-sample proportion test are displayed on top of each
bar.

```{r ggbarstats1, fig.height=8, fig.width=10}
## for reproducibility
set.seed(123)
library(ggplot2)

## plot
ggbarstats(
  data             = movies_long,
  x                = mpaa,
  y                = genre,
  title            = "MPAA Ratings by Genre",
  xlab             = "movie genre",
  legend.title     = "MPAA rating",
  ggtheme          = hrbrthemes::theme_ipsum_pub(),
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  palette          = "Set2"
)
```

**Defaults** return<br>

✅ descriptives (frequency + %s) <br>
✅ inferential statistics <br>
✅ effect size + CIs <br>
✅ Goodness-of-fit tests <br>
✅ Bayesian hypothesis-testing <br>
✅ Bayesian estimation <br>


And, needless to say, there is also a `grouped_` variant of this function-

```{r ggbarstats2, fig.height=6, fig.width=12}
## setup
set.seed(123)

## plot
grouped_ggbarstats(
  data         = mtcars,
  x            = am,
  y            = cyl,
  grouping.var = vs,
  package      = "wesanderson",
  palette      = "Darjeeling2",
  ggtheme      = ggthemes::theme_tufte(base_size = 12)
)
```

##### Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
bars | `ggplot2::geom_bar` | ❌
descriptive labels | `ggplot2::geom_label` | `label.args`

##### Summary of tests

**two-way table**

**Hypothesis testing**

Type | Design | Test | Function used
----------- | ----- | ------------------------- | -----
Parametric/Non-parametric | Unpaired | Pearson's $\chi^2$ test | `stats::chisq.test`
Bayesian | Unpaired | Bayesian Pearson's $\chi^2$ test | `BayesFactor::contingencyTableBF`
Parametric/Non-parametric | Paired  | McNemar's $\chi^2$ test | `stats::mcnemar.test`
Bayesian | Paired  | ❌ | ❌

**Effect size estimation**

Type | Design | Effect size | CI? | Function used
----------- | ----- | ------------------------- | --- | -----
Parametric/Non-parametric | Unpaired | Cramer's $V$ | ✅ | `effectsize::cramers_v`
Bayesian | Unpaired | Cramer's $V$ | ✅ | `effectsize::cramers_v`
Parametric/Non-parametric | Paired | Cohen's $g$ | ✅ | `effectsize::cohens_g`
Bayesian | Paired | ❌ | ❌ | ❌

**one-way table**

**Hypothesis testing**

Type | Test | Function used
----------- | ------------------------- | -----
Parametric/Non-parametric | Goodness of fit $\chi^2$ test | `stats::chisq.test`
Bayesian | Bayesian Goodness of fit $\chi^2$ test | (custom)

**Effect size estimation**

Type | Effect size | CI? | Function used
----------- | ------------------------- | --- | -----
Parametric/Non-parametric | Pearson's $C$ | ✅ | `effectsize::pearsons_c`
Bayesian | ❌ | ❌ | ❌

### `ggcoefstats`

The function `ggcoefstats` generates **dot-and-whisker plots** for
regression models saved in a tidy data frame. The tidy dataframes are prepared
using `parameters::model_parameters`. Additionally, if available, the model
summary indices are also extracted from `performance::model_performance`.

Although the statistical models displayed in the plot may differ based on the
class of models being investigated, there are few aspects of the plot that will
be invariant across models:

  - The dot-whisker plot contains a dot representing the **estimate** and their
    **confidence intervals** (`95%` is the default). The estimate can either be
    effect sizes (for tests that depend on the `F`-statistic) or regression
    coefficients (for tests with `t`-, $\chi^{2}$-, and `z`-statistic), etc. The
    function will, by default, display a helpful `x`-axis label that should
    clear up what estimates are being displayed. The confidence intervals can
    sometimes be asymmetric if bootstrapping was used.
    
  - The label attached to dot will provide more details from the statistical
    test carried out and it will typically contain estimate, statistic, and
    *p*-value.

  - The caption will contain diagnostic information, if available, about
    models that can be useful for model selection: The smaller the Akaike's
    Information Criterion (**AIC**) and the Bayesian Information Criterion
    (**BIC**) values, the "better" the model is.

  - The output of this function will be a `{ggplot2}` object and, thus, it can be
    further modified (e.g., change themes, etc.) with `{ggplot2}` functions.

```{r ggcoefstats1, fig.height=5, fig.width=6}
## for reproducibility
set.seed(123)

## model
mod <- stats::lm(formula = mpg ~ am * cyl, data = mtcars)

## plot
ggcoefstats(mod, ggtheme = hrbrthemes::theme_ipsum_ps())
```

**Defaults** return<br>

✅ inferential statistics <br>
✅ estimate + CIs <br>
✅ model summary (AIC and BIC) <br>

##### Supported models

Most of the regression models that are supported in the underlying packages are
also supported by `ggcoefstats`. 

```{r supported}
insight::supported_models()
```

Although not shown here, this function can also be used to carry out parametric,
robust, and Bayesian random-effects meta-analysis.

##### Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
regression estimate | `ggplot2::geom_point` | `point.args`
error bars | `ggplot2::geom_errorbarh` | `errorbar.args`
vertical line | `ggplot2::geom_vline` | `vline.args`
label with statistical details | `ggrepel::geom_label_repel` | `stats.label.args`

##### Summary of meta-analysis tests

**Hypothesis testing** and **Effect size estimation**

Type | Test | Effect size | CI? | Function used
----------- | -------------------- | -------- | ---  | -----
Parametric | Meta-analysis via random-effects models | $\beta$ | ✅ | `metafor::metafor`
Robust | Meta-analysis via robust random-effects models | $\beta$ | ✅ | `metaplus::metaplus`
Bayes | Meta-analysis via Bayesian random-effects models | $\beta$ | ✅ | `metaBMA::meta_random`

For a more exhaustive account of this function, see the associated vignette-
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcoefstats.html>

### Extracting dataframes with statistical details

`{ggstatsplot}` also offers a convenience function to extract dataframes with
statistical details that are used to create expressions displayed in
`{ggstatsplot}` plots.

```{r extract_stats}
set.seed(123)

## a list of tibbles containing statistical analysis summaries
ggbetweenstats(mtcars, cyl, mpg) %>%
  extract_stats()
```

Note that all of this analysis is carried out by `{statsExpressions}`
package: <https://indrajeetpatil.github.io/statsExpressions/>

### Using `{ggstatsplot}` statistical details with custom plots

Sometimes you may not like the default plots produced by `{ggstatsplot}`. In such
cases, you can use other **custom** plots (from `{ggplot2}` or other plotting
packages) and still use `{ggstatsplot}` functions to display results from relevant
statistical test.

For example, in the following chunk, we will create plot (*ridgeplot*) using
`ggridges` package and use `{ggstatsplot}` function for extracting results.

```{r ridgeplot, fig.height=10, fig.width=7}
## loading the needed libraries
set.seed(123)
library(ggridges)
library(ggplot2)
library(ggstatsplot)

## using `{ggstatsplot}` to get call with statistical results
stats_results <-
  ggbetweenstats(
    data = morley,
    x = Expt,
    y = Speed,
    output = "subtitle"
  )

## using `ggridges` to create plot
ggplot(morley, aes(x = Speed, y = as.factor(Expt), fill = as.factor(Expt))) +
  geom_density_ridges(
    jittered_points = TRUE,
    quantile_lines = TRUE,
    scale = 0.9,
    alpha = 0.7,
    vline_size = 1,
    vline_color = "red",
    point_size = 0.4,
    point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)
  ) + ## adding annotations
  labs(
    title = "Michelson-Morley experiments",
    subtitle = stats_results,
    x = "Speed of light",
    y = "Experiment number"
  ) + ## remove the legend
  theme(legend.position = "none")
```

## Summary of benefits of using `{ggstatsplot}`

- No need to use scores of packages for statistical analysis
  (e.g., one to get stats, one to get effect sizes, another to get Bayes
  Factors, and yet another to get pairwise comparisons, etc.).

- Minimal amount of code needed for all functions (typically only `data`, `x`,
  and `y`), which minimizes chances of error and makes for tidy scripts.
  
- Conveniently toggle between statistical approaches.

- Truly makes your figures worth a thousand words.

- No need to copy-paste results to the text editor (MS-Word, e.g.).

- Disembodied figures stand on their own and are easy to evaluate for the reader.

- More breathing room for theoretical discussion and other text.

- No need to worry about updating figures and statistical details separately.

## Misconceptions about `{ggstatsplot}`

This package is... 

❌ an alternative to learning `{ggplot2}`<br>
✅ (The better you know `{ggplot2}`, the more you can modify the defaults to your
liking.)

❌ meant to be used in talks/presentations<br>
✅ (Default plots can be too complicated for effectively communicating results in
time-constrained presentation settings, e.g. conference talks.)

❌ the only game in town<br>
✅ (GUI software alternatives: [JASP](https://jasp-stats.org/) and [jamovi](https://www.jamovi.org/)).

## Extensions

In case you use the GUI software [`jamovi`](https://www.jamovi.org/), you can
install a module called [`jjstatsplot`](https://github.com/sbalci/jjstatsplot),
which is a wrapper around `{ggstatsplot}`.

## Acknowledgments

I would like to thank all the contributors to `{ggstatsplot}` who pointed out
bugs or requested features I hadn't considered. I would especially like to thank
other package developers (especially Daniel Lüdecke, Dominique Makowski, Mattan
S. Ben-Shachar, Brenton Wiernik, Patrick Mair, Salvatore Mangiafico, etc.) who
have patiently and diligently answered my relentless number of questions and
added feature requests I wanted. I also want to thank Chuck Powell for his
initial contributions to the package.

The hexsticker was generously designed by Sarah Otterstetter (Max Planck
Institute for Human Development, Berlin). This package has also benefited from
the larger `rstats` community on Twitter and `StackOverflow`.

Thanks are also due to my postdoc advisers (Mina Cikara and Fiery Cushman at
Harvard University; Iyad Rahwan at Max Planck Institute for Human Development)
who patiently supported me spending hundreds (?) of hours working on this
package rather than what I was paid to do. 😁

## Contributing

I'm happy to receive bug reports, suggestions, questions, and (most of all)
contributions to fix problems and add features. I personally prefer using the
`GitHub` issues system over trying to reach out to me in other ways (personal
e-mail, Twitter, etc.). Pull Requests for contributions are encouraged.

Here are some simple ways in which you can contribute (in the increasing order
of commitment):

  - Read and correct any inconsistencies in the
    [documentation](https://indrajeetpatil.github.io/ggstatsplot/)
  - Raise issues about bugs or wanted features
  - Review code
  - Add new functionality (in the form of new plotting functions or helpers for
    preparing subtitles)

Please note that this project is released with a 
[Contributor Code of Conduct](https://github.com/IndrajeetPatil/ggstatsplot/blob/master/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.
---
title: "Visualizations with statistical details: The 'ggstatsplot' approach"
tags:
  - R
  - parametric statistics
  - nonparametric statistics
  - robust statistics
  - Bayesian statistics
  - ggplot2
  - ggplot2-extension
authors:
  - name: Indrajeet Patil
    orcid: 0000-0003-1995-6531
    affiliation: 1
affiliations:
  - name: Center for Humans and Machines, Max Planck Institute for Human Development, Berlin, Germany
    index: 1
bibliography: paper.bib
date: "`r Sys.Date()`"
output: rticles::joss_article
csl: apa.csl
journal: JOSS
link-citations: yes
---

```{r, warning=FALSE, message=FALSE, echo=FALSE}
# needed libraries
library(ggstatsplot)

knitr::opts_chunk$set(
  collapse = TRUE,
  out.width = "100%",
  dpi = 450,
  comment = "#>",
  error = TRUE,
  message = FALSE,
  warning = FALSE
)

set.seed(123) # for reproducibility
```

# Summary

Graphical displays can reveal problems in a statistical model that might not be
apparent from purely numerical summaries. Such visualizations can also be
helpful for the reader to evaluate the validity of a model if it is reported in
a scholarly publication or report. But, given the onerous costs involved,
researchers often avoid preparing information-rich graphics and exploring
several statistical approaches or tests available. The `{ggstatsplot}` package in
the R programming language [@base2021] provides a one-line syntax to enrich
`{ggplot2}`-based visualizations with the results from statistical analysis
embedded in the visualization itself. In doing so, the package helps researchers
adopt a rigorous, reliable, and robust data exploratory and reporting workflow.

# Statement of Need

In a typical data analysis workflow, data visualization and statistical modeling
are two different phases: visualization informs modeling, and in turn, modeling
can suggest a different visualization method, and so on and so forth
[@wickham2016r]. The central idea of `{ggstatsplot}` is simple: combine these two
phases into one in the form of an informative graphic with statistical details.

Before discussing benefits of this approach, we will show an example (Figure
1).

```{r penguins, fig.width=8, fig.height=5, fig.cap="Example plot from the `{ggstatsplot}` package illustrating its philosophy of juxtaposing informative visualizations with details from statistical analysis. To see all supported plots and statistical analyses, see the package website: \\url{https://indrajeetpatil.github.io/ggstatsplot/}"}
library(palmerpenguins) # for 'penguins' dataset
library(ggstatsplot)

ggbetweenstats(penguins, species, body_mass_g)
```

As can be seen, with a single line of code, the function produces details about
descriptive statistics, inferential statistics, effect size estimate and its
uncertainty, pairwise comparisons, Bayesian hypothesis testing, Bayesian
posterior estimate and its uncertainty. Moreover, these details are juxtaposed
with informative and well-labeled visualizations. The defaults are designed to
follow best practices in both data visualization [@cleveland1985;
@grant2018data; @healy2018data; @tufte2001; @wilke2019fundamentals] and
(frequentist/Bayesian) statistical reporting [@apa2019; @van2020jasp]. Without
`{ggstatsplot}`, getting these statistical details and customizing a plot would
require significant amount of time and effort. In other words, this package
removes the trade-off often faced by researchers between ease and thoroughness
of data exploration and further cements good data exploration habits.

Internally, data cleaning is carried out using the `tidyverse` [@Wickham2019],
while statistical analysis is carried out via the `{statsExpressions}`
[@Patil2021] and `easystats` [@Ben-Shachar2020; @Lüdecke2020parameters;
@Lüdecke2020performance;
@Lüdecke2019; @Makowski2019; @Makowski2020] packages. All visualizations are
constructed using the grammar of graphics framework [@Wilkinson2012], as
implemented in the `{ggplot2}` package [@Wickham2016].

# Benefits

In summary, the benefits of `{ggstatsplot}`'s approach are the following. It:

a. produces charts displaying both raw data, and numerical plus graphical
   summary indices,

b. avoids errors in and increases reproducibility of statistical reporting,

c. highlights the importance of the effect by providing effect size measures by
   default,

d. provides an easy way to evaluate *absence* of an effect using Bayes factors,

e. encourages researchers and readers to evaluate statistical assumptions of a
model in the context of the underlying data (Figure 2),

f. is easy and simple enough that someone with little to no coding experience
   can use it without making an error and may even encourage beginners to
   programmatically analyze data, instead of using GUI software.

```{r reporting, echo=FALSE, fig.cap="Comparing the 'Standard' approach of reporting statistical analysis in a publication/report with the 'ggstatsplot' approach of reporting the same analysis next to an informative graphic. Note that the results described in the 'Standard' approach are about the 'Dinosaur' dataset plotted on the right. Without the accompanying visualization, it is hard to evaluate the validity of the results. The ideal reporting practice will be a hybrid of these two approaches where the plot contains both the visual and numerical summaries about a statistical model, while the narrative provides interpretative context for the reported statistics."}
knitr::include_graphics("reporting.png")
```

# Future Scope

This package is an ambitious, ongoing, and long-term project. It currently
supports common statistical tests (parametric, non-parametric, robust, or
Bayesian *t*-test, one-way ANOVA, contingency table analysis, correlation
analysis, meta-analysis, regression analyses, etc.) and corresponding
visualizations (box/violin plot, scatter plot, dot-and-whisker plot, pie chart,
bar chart, etc.). It will continue expanding to support an increasing
collection of statistical analyses and visualizations.

# Licensing and Availability

`{ggstatsplot}` is licensed under the GNU General Public License (v3.0), with all
source code stored at [GitHub](https://github.com/IndrajeetPatil/ggstatsplot/).
In the spirit of honest and open science, requests and suggestions for fixes,
feature updates, as well as general questions and concerns are encouraged via
direct interaction with contributors and developers by filing an
[issue](https://github.com/IndrajeetPatil/ggstatsplot/issues) while respecting
[*Contribution
Guidelines*](https://indrajeetpatil.github.io/ggstatsplot/CONTRIBUTING.html).

# Acknowledgements

I would like to acknowledge the support of Mina Cikara, Fiery Cushman, and Iyad
Rahwan during the development of this project. `{ggstatsplot}` relies heavily on
the [`easystats`](https://github.com/easystats/easystats) ecosystem, a
collaborative project created to facilitate the usage of `R` for statistical
analyses. Thus, I would like to thank the
[members](https://github.com/orgs/easystats/people) of `easystats` as well as
the users. I would additionally like to thank the contributors to `{ggstatsplot}`
for reporting bugs, providing helpful feedback, or helping with enhancements.

# References

---
title: "additional vignettes"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{additional vignettes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)
```

## Additional vignettes

Due to size constraints, all available vignettes are only available on the
website for this package: <br>
<https://indrajeetpatil.github.io/ggstatsplot/articles/>

### Vignettes for individual functions

  - `ggbetweenstats`:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggbetweenstats.html>

  - `ggwithinstats`:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggwithinstats.html>

  - `ggcorrmat`:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcorrmat.html>

  - `gghistostats`:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/gghistostats.html>

  - `ggdotplotstats`:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggdotplotstats.html>

  - `ggpiestats`:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggpiestats.html>

  - `ggscatterstats`:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggscatterstats.html>

  - `ggcoefstats`:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcoefstats.html>

### General vignettes

  - Frequently Asked Questions (FAQ):
  <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/faq.html>

  - Graphic design and statistical reporting principles guiding `{ggstatsplot}`
    development:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/principles.html>

  - Examples illustrating how to use `{purrr}` to extend `{ggstatsplot}` functionality:
    <https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html>

- Pairwise comparisons with `{ggstatsplot}`:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/pairwise.html>

## Publication

A journal articles describing philosophy and principles behind this package:
<https://joss.theoj.org/papers/10.21105/joss.03167>

## Presentation

In addition to these vignettes, another quick way to get an overview of this
package is to go through the following slides:
<https://indrajeetpatil.github.io/ggstatsplot_slides/slides/ggstatsplot_presentation.html#1>

## Statistical backend of `{ggstatsplot}`

The `{statsExpressions}` package forms the statistical backend that processes data
and creates dataframes and expressions containing results from statistical
tests.

For more exhaustive documentation for this package, see:
<https://indrajeetpatil.github.io/statsExpressions/>

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>

---
title: "ggwithinstats"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{ggwithinstats}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
## show me all columns
options(tibble.width = Inf, pillar.bold = TRUE, pillar.subtle_num = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>"
)

if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(PMCMRplus)
}

library(ggstatsplot)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

Lifecycle: [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

The function `ggwithinstats` is designed to facilitate 
**data exploration**, and for making highly customizable **publication-ready plots**,
with relevant statistical details included in the plot itself if desired. We
will see examples of how to use this function in this vignette.

To begin with, here are some instances where you would want to use
`ggwithinstats`-

 - to check if a continuous variable differs across multiple groups/conditions

 - to compare distributions visually and check for outliers

**Note**: This vignette uses the pipe operator (`%>%`), if you are not
familiar with this operator, here is a good explanation:
<http://r4ds.had.co.nz/pipes.html>

## Comparisons between groups with `ggwithinstats`

To illustrate how this function can be used, we will use the `bugs` dataset
throughout this vignette. This data set, "Bugs", provides the extent to which
men and women want to kill arthropods that vary in freighteningness (low, high)
and disgustingness (low, high). Each participant rates their attitudes towards
all anthropods. Subset of the data reported by [Ryan et al. (2013)](https://www.sciencedirect.com/science/article/pii/S0747563213000277). 
Note that this is a repeated measures design because the same participant gave
four different ratings across four different conditions (LDLF, LDHF, HDLF,
HDHF).

Suppose the first thing we want to inspect is the distribution of desire to kill
across all conditions (disregarding the factorial structure of the experiment).
We also want to know if the mean differences in this desire across conditions is
statistically significant.

The simplest form of the function call is-

```{r ggwithinstats1, fig.height = 6, fig.width = 8}
## since the confidence intervals for the effect sizes are computed using
## bootstrapping, important to set a seed for reproducibility
set.seed(123)
library(ggstatsplot)

## function call
ggwithinstats(
  data = bugs_long,
  x = condition,
  y = desire
)
```

**Note**:

  - The function automatically decides whether a dependent samples test is
    preferred (for 2 groups) or an ANOVA (3 or more groups). based on the number
    of levels in the grouping variable.
    
  - The output of the function is a `ggplot` object which means that it can be
    further modified with `{ggplot2}` functions.

As can be seen from the plot, the function by default returns Bayes Factor for
the test. If the null hypothesis can't be rejected with the null hypothesis
significance testing (NHST) approach, the Bayesian approach can help index
evidence in favor of the null hypothesis (i.e., $BF_{01}$).

By default, natural logarithms are shown because Bayes Factor values can
sometimes be pretty large. Having values on logarithmic scale also makes it easy
to compare evidence in favor alternative ($BF_{10}$) versus null ($BF_{01}$)
hypotheses (since $log_{e}(BF_{01}) = - log_{e}(BF_{10})$).

We can make the output much more aesthetically pleasing as well as informative
by making use of the many optional parameters in `ggwithinstats`. We'll add a
title and caption, better `x` and `y` axis labels, and tag and label the
outliers in the data. We can and will change the overall theme as well as the
color palette in use.

```{r ggwithinstats2, fig.height = 6, fig.width = 8}
## for reproducibility
set.seed(123)
library(ggstatsplot)

## plot
ggwithinstats(
  data = bugs_long,
  x = condition,
  y = desire,
  type = "nonparametric", ## type of statistical test
  xlab = "Condition", ## label for the x-axis
  ylab = "Desire to kill an artrhopod", ## label for the y-axis
  effsize.type = "biased", ## type of effect size
  sphericity.correction = FALSE, ## don't display sphericity corrected dfs and p-values
  pairwise.comparisons = TRUE, ## display pairwise comparisons
  outlier.tagging = TRUE, ## whether outliers should be flagged
  outlier.coef = 1.5, ## coefficient for Tukey's rule
  outlier.label = region, ## label to attach to outlier values
  outlier.label.color = "red", ## outlier point label color
  mean.plotting = TRUE, ## whether the mean is to be displayed
  mean.color = "darkblue", ## color for mean
  package = "yarrr", ## package from which color palette is to be taken
  palette = "info2", ## choosing a different color palette
  title = "Comparison of desire to kill bugs",
  caption = "Source: Ryan et al., 2013"
) + ## modifying the plot further
  ggplot2::scale_y_continuous(
    limits = c(0, 10),
    breaks = seq(from = 0, to = 10, by = 1)
  )
```

As can be appreciated from the effect size (partial eta squared) of 0.18, there
are small differences in the mean desire to kill across conditions.
Importantly, this plot also helps us appreciate the distributions within any
given condition. 

So far we have only used a classic parametric test, but we can also use other
available options: The `type` (of test) argument also accepts the following
abbreviations: `"p"` (for *parametric*), `"np"` (for *nonparametric*), `"r"`
(for *robust*), `"bf"` (for *Bayes Factor*).

Let's use the `combine_plots` function to make one plot from four separate plots
that demonstrates all of these options. Let's compare desire to kill bugs only
for low versus high disgust conditions to see how much of a difference whether a
bug is disgusting-looking or not makes to the desire to kill that bug. We will
generate the plots one by one and then use `combine_plots` to merge them into
one plot with some common labeling. It is possible, but not necessarily
recommended, to make each plot have different colors or themes.

For example,

```{r ggwithinstats3, fig.height = 14, fig.width = 14}
## for reproducibility
set.seed(123)
library(ggstatsplot)

## selecting subset of the data
df_disgust <- dplyr::filter(bugs_long, condition %in% c("LDHF", "HDHF"))

## parametric t-test
p1 <- ggwithinstats(
  data = df_disgust,
  x = condition,
  y = desire,
  type = "p",
  effsize.type = "d",
  conf.level = 0.99,
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg"
)

## Mann-Whitney U test (nonparametric test)
p2 <- ggwithinstats(
  data = df_disgust,
  x = condition,
  y = desire,
  xlab = "Condition",
  ylab = "Desire to kill bugs",
  type = "np",
  conf.level = 0.99,
  title = "Non-parametric Test",
  package = "ggsci",
  palette = "uniform_startrek"
)

## robust t-test
p3 <- ggwithinstats(
  data = df_disgust,
  x = condition,
  y = desire,
  xlab = "Condition",
  ylab = "Desire to kill bugs",
  type = "r",
  conf.level = 0.99,
  title = "Robust Test",
  package = "wesanderson",
  palette = "Royal2"
)

## Bayes Factor for parametric t-test
p4 <- ggwithinstats(
  data = df_disgust,
  x = condition,
  y = desire,
  xlab = "Condition",
  ylab = "Desire to kill bugs",
  type = "bayes",
  title = "Bayesian Test",
  package = "ggsci",
  palette = "nrc_npg"
)

## combining the individual plots into a single plot
combine_plots(
  plotlist = list(p1, p2, p3, p4),
  plotgrid.args = list(nrow = 2),
  annotation.args = list(
    title = "Effect of disgust on desire to kill bugs ",
    caption = "Source: Bugs dataset from `jmv` R package"
  )
)
```

## Grouped analysis with `grouped_ggwithinstats`

What if we want to carry out this same analysis but for each region (or gender)?

`{ggstatsplot}` provides a special helper function for such instances:
`grouped_ggwithinstats`. This is merely a wrapper function around
`combine_plots`. It applies `ggwithinstats` across all **levels**
of a specified **grouping variable** and then combines list of individual plots
into a single plot. Note that the grouping variable can be anything: conditions
in a given study, groups in a study sample, different studies, etc.

Let's focus on the two regions and for years: 1967, 1987, 2007. Also,
let's carry out pairwise comparisons to see if there differences between every
pair of continents. 

```{r grouped1, fig.height = 14, fig.width = 8}
## for reproducibility
set.seed(123)
library(ggstatsplot)

grouped_ggwithinstats(
  ## arguments relevant for ggwithinstats
  data = bugs_long,
  x = condition,
  y = desire,
  grouping.var = gender,
  xlab = "Continent",
  ylab = "Desire to kill bugs",
  type = "nonparametric", ## type of test
  pairwise.display = "significant", ## display only significant pairwise comparisons
  p.adjust.method = "BH", ## adjust p-values for multiple tests using this method
  # ggtheme = ggthemes::theme_tufte(),
  package = "ggsci",
  palette = "default_jco",
  outlier.tagging = TRUE,
  outlier.label = education,
  k = 3,
  ## arguments relevant for combine_plots
  annotation.args = list(title = "Desire to kill bugs across genders"),
  plotgrid.args = list(ncol = 1)
)
```

## Grouped analysis with `ggwithinstats` + `{purrr}` 

Although this grouping function provides a quick way to explore the data, it
leaves much to be desired. For example, the same type of test and theme is
applied for all genders, but maybe we want to change this for different genders,
or maybe we want to gave different effect sizes for different years. This type
of customization for different levels of a grouping variable is not possible
with `grouped_ggwithinstats`, but this can be easily achieved using the `{purrr}`
package.

See the associated vignette here:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html>

## Between-subjects designs

For independent measures designs, `ggbetweenstats` function can be used:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggbetweenstats.html>

## Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
raw data | `ggplot2::geom_point` | `point.args`
point path | `ggplot2::geom_path` | `point.path.args`
box plot | `ggplot2::geom_boxplot` | `boxplot.args`
density plot | `ggplot2::geom_violin` | `violin.args`
centrality measure point | `ggplot2::geom_point` | `centrality.point.args`
centrality measure point path | `ggplot2::geom_path` | `centrality.path.args`
centrality measure label | `ggrepel::geom_label_repel` | `centrality.label.args`
outlier point | `ggplot2::stat_boxplot` | ❌
outlier label | `ggrepel::geom_label_repel` | `outlier.label.args`
pairwise comparisons | `ggsignif::geom_signif` | `ggsignif.args`

## Summary of tests

**Central tendency measure**

Type | Measure | Function used
----------- | --------- | ------------------ 
Parametric | mean | `parameters::describe_distribution`
Non-parametric | median | `parameters::describe_distribution`
Robust | trimmed mean | `parameters::describe_distribution`
Bayesian | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution`

**Hypothesis testing**

Type | No. of groups | Test | Function used
----------- | --- | ------------------------- | -----
Parametric | > 2 | One-way repeated measures ANOVA | `afex::aov_ez`
Non-parametric | > 2 | Friedman rank sum test | `stats::friedman.test`
Robust | > 2 | Heteroscedastic one-way repeated measures ANOVA for trimmed means | `WRS2::rmanova`
Bayes Factor | > 2 | One-way repeated measures ANOVA | `BayesFactor::anovaBF`
Parametric | 2 | Student's *t*-test | `stats::t.test`
Non-parametric | 2 | Wilcoxon signed-rank test | `stats::wilcox.test`
Robust | 2 | Yuen's test on trimmed means for dependent samples | `WRS2::yuend`
Bayesian | 2 | Student's *t*-test | `BayesFactor::ttestBF`

**Effect size estimation**

Type | No. of groups | Effect size | CI? | Function used
----------- | --- | ------------------------- | --- | -----
Parametric | > 2 | $\eta_{p}^2$, $\omega_{p}^2$ | ✅ | `effectsize::omega_squared`, `effectsize::eta_squared`
Non-parametric | > 2 | $W_{Kendall}$ (Kendall's coefficient of concordance) | ✅ | `effectsize::kendalls_w`
Robust | > 2 | $\delta_{R-avg}^{AKP}$ (Algina-Keselman-Penfield robust standardized difference average) | ✅ | `WRS2::wmcpAKP`
Bayes Factor | > 2 | $R_{Bayesian}^2$ | ✅ | `performance::r2_bayes`
Parametric | 2 | Cohen's *d*, Hedge's *g* | ✅ | `effectsize::cohens_d`, `effectsize::hedges_g`
Non-parametric | 2 | *r* (rank-biserial correlation) | ✅ | `effectsize::rank_biserial`
Robust | 2 |  $\delta_{R}^{AKP}$ (Algina-Keselman-Penfield robust standardized difference) | ✅ | `WRS2::wmcpAKP`
Bayesian | 2 | $\delta_{posterior}$ | ✅ | `bayestestR::describe_posterior`

**Pairwise comparison tests**

Type | Test | *p*-value adjustment? | Function used
----------- | ---------------------------- | --- | -----
Parametric | Student's *t*-test | ✅ | `stats::pairwise.t.test`
Non-parametric | Durbin-Conover test | ✅ | `PMCMRplus::durbinAllPairsTest` 
Robust | Yuen's trimmed means test | ✅ | `WRS2::rmmcp`
Bayesian | Student's *t*-test | ❌ | `BayesFactor::ttestBF`

## Reporting

If you wish to include statistical analysis results in a publication/report, the
ideal reporting practice will be a hybrid of two approaches:

- the `{ggstatsplot}` approach, where the plot contains both the visual and
numerical summaries about a statistical model, and

- the *standard* narrative approach, which provides interpretive context for the
reported statistics.

For example, let's see the following example:



```{r reporting}
library(WRS2) ## for data
ggwithinstats(WineTasting, Wine, Taste)
```

The narrative context (assuming `type = "parametric"`) can complement this plot
either as a figure caption or in the main text-

> Fisher's repeated measures one-way ANOVA revealed that, across 22 friends to
taste each of the three wines, there was a statistically significant difference
across persons preference for each wine. The effect size $(\omega_{p} = 0.02)$
was medium, as per Field’s (2013) conventions. The Bayes Factor for the same
analysis revealed that the data were `r round(exp(2.11), 2)` times more probable
under the alternative hypothesis as compared to the null hypothesis. This can be
considered moderate evidence (Jeffreys, 1961) in favor of the alternative
hypothesis. This global effect was carried out by post hoc pairwise *t*-tests,
which revealed that Wine C was preferred across participants to be the least
desirable compared to Wines A and B.

Similar reporting style can be followed when the function performs *t*-test
instead of a one-way ANOVA.

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "using 'ggstatsplot' with the 'purrr' package"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    eval: FALSE
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{using 'ggstatsplot' with the 'purrr' package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
## pretty tibble printing
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>",
  package.startup.message = FALSE
)

future::plan("multicore")
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

---

This is an extremely time-consuming vignette, and so it is not evaluated here.
You can still use the code as a reference for writing your own `{purrr}` code.

---

## Why use `{purrr}`?

Most of the `{ggstatsplot}` functions have `grouped_` variants, which are designed
to quickly run the same `{ggstatsplot}` function across multiple levels of a
**single** grouping variable. Although this function is useful for data
exploration, it has two strong weaknesses-

 * The arguments applied to `grouped_` function call are applied uniformly to
   all levels of the grouping variable when we might want to customize them for
   different levels of the grouping variable.

 * Only one grouping variable can be used to repeat the analysis when in reality
   there can be a combination of grouping variables and the operation needs to
   be repeated for all resulting combinations.

We will see how to overcome this limitation by combining `{ggstatsplot}` with the
`{purrr}` package.

**Note:**

 * While using `purrr::pmap()`, we **must** input the arguments as strings.

 * You can use `{ggplot2}` themes from extension packages (e.g. `ggthemes`).

 * If you'd like some more background or an introduction to the purrr package,
   please see [this chapter](https://adv-r.hadley.nz/functionals.html).

## Introduction and methodology

For all the examples in this vignette we are going to build `list`s of things
that we will pass along to `{purrr}` which will in turn return a list of plots
that will be passed to `combine_plots`. As the name implies `combine_plots`
merges the individual plots into one bigger plot with common labeling and
aesthetics.

What are these `lists` that we are building? The lists correspond to the
parameters in our `{ggstatsplot}` function like `ggbetweenstats`. If you look at
the help file for `?ggbetweenstats` for example the very first parameter it
wants is the `data` file we'll be using. We can also pass it different `titles`
of even `ggtheme` themes.

You can pass:

 * A single character string such as `xlab = "Continent"` or numeric such as
   `nboot = 25` in which case it will be reused/recycled as many times as
   needed.

 * A vector of values such as `nboot = c(50, 100, 200)` in which case it
   will be coerced to a list and checked for the right class (in this case
   integer) and the right quantity of entries in the vector i.e.,
   `nboot = c(50, 100)` will fail if we're trying to make three plots.

 * A list; either named `data = year_list` or created as you go
   `palette = list("Dark2", "Set1")`. Any list will
   be checked for the right class (in this case character) and the right
   quantity of entries in the list.

## `ggbetweenstats`

Let's start with `ggebtweenstats`. We'll use the `gapminder` dataset. We'll make
a 3 item `list` called `year_list` using `dplyr::filter` and `split`.

```r
## for reproducibility
library(ggstatsplot)
set.seed(123)

## let's split the dataframe and create a list by years of interest
year_list <- gapminder::gapminder %>%
  dplyr::filter(year %in% c(1967, 1987, 2007), continent != "Oceania") %>%
  split(f = .$year, drop = TRUE)

## checking the length of the list and the names of each element
length(year_list)
names(year_list)
```

Now that we have the data divided into the three relevant years in a list we'll
turn to `purrr::pmap` to create a list of `ggplot` objects that we'll make use of
stored in `plot_list`. When you look at the documentation for `?pmap` it will
accept `.l` which is a list of lists. The length of `.l` determines the number of
arguments that `.f` will be called with. List names will be used if present.
`.f` is the function we want to apply (here, `.f = ggbetweenstats`).

Let's keep building the list of arguments, `.l`. First is `data = year_list`,
the `x` and `y` axes are constant in all three plots so we pass the variable
name as a string `x = "continent"`.

Same with the label we'll use for outliers where needed. For demonstration
purposes let's assume we want the outliers on each plot to be a different color.
Not actually recommending it just demonstrating what's possible. The rest of the
code shows you a wide variety of possibilities and we won't catalog them here.

```r
## for reproducibility
set.seed(123)
library(ggstatsplot)

## creating a list of plots
plot_list <- purrr::pmap(
  .l = list(
    data = year_list,
    x = "continent",
    y = "lifeExp",
    outlier.tagging = TRUE,
    outlier.label = "country",
    outlier.label.args = list(
      list(size = 3, color = "#56B4E9"),
      list(size = 2.5, color = "#009E73"),
      list(size = 3.5, color = "#F0E442")
    ),
    xlab = "Continent",
    ylab = "Life expectancy",
    title = list(
      "Year: 1967",
      "Year: 1987",
      "Year: 2007"
    ),
    type = list("r", "bf", "np"),
    pairwise.display = list("s", "ns", "all"),
    p.adjust.method = list("hommel", "bonferroni", "BH"),
    conf.level = list(0.99, 0.95, 0.90),
    k = list(1, 2, 3),
    effsize.type = list(
      NULL,
      "partial_omega",
      "partial_eta"
    ),
    plot.type = list("box", "boxviolin", "violin"),
    package = list("nord", "ochRe", "awtools"),
    palette = list("aurora", "parliament", "bpalette"),
    ggtheme = list(
      ggthemes::theme_stata(),
      ggplot2::theme_classic(),
      ggthemes::theme_fivethirtyeight()
    )
  ),
  .f = ggbetweenstats
)
```

The final step is to pass the `plot_list` object we just created to the
`combine_plots` function. While each of the 3 plots already has labeling
information `combine_plots` gives us an opportunity to add additional details to
the merged plots and specify the layout in rows and columns.

```r
## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(title = "Changes in life expectancy across continents (1967-2007)"),
  plotgrid.args = list(ncol = 1)
)
```

## `ggwithinstats`

We will be using simulated data from then Attention Network Test provided in ANT
dataset in `ez` package.

```r
## for reproducibility
set.seed(123)
library(ggstatsplot)
library(ez)
data("ANT") ## loading data from `ez` package

## let's split the dataframe and create a list by years of interest
cue_list <- ANT %>% split(f = .$cue, drop = TRUE)

## checking the length of the list and the names of each element
length(cue_list)

## creating a list of plots by applying the same function for elements of the list
plot_list <- purrr::pmap(
  .l = list(
    data = cue_list,
    x = "flank",
    y = "rt",
    outlier.tagging = TRUE,
    outlier.label = "group",
    outlier.coef = list(2, 2, 2.5, 3),
    outlier.label.args = list(
      list(size = 3, color = "#56B4E9"),
      list(size = 2.5, color = "#009E73"),
      list(size = 4, color = "#F0E442"),
      list(size = 2, color = "red")
    ),
    xlab = "Flank",
    ylab = "Response time",
    title = list(
      "Cue: None",
      "Cue: Center",
      "Cue: Double",
      "Cue: Spatial"
    ),
    type = list("p", "r", "bf", "np"),
    pairwise.display = list("ns", "s", "ns", "all"),
    p.adjust.method = list("fdr", "hommel", "bonferroni", "BH"),
    conf.level = list(0.99, 0.99, 0.95, 0.90),
    k = list(3, 2, 2, 3),
    effsize.type = list(
      "omega",
      "eta",
      "partial_omega",
      "partial_eta"
    ),
    package = list("ggsci", "palettetown", "palettetown", "wesanderson"),
    palette = list("lanonc_lancet", "venomoth", "blastoise", "GrandBudapest1"),
    ggtheme = list(
      ggplot2::theme_linedraw(),
      hrbrthemes::theme_ft_rc(),
      ggthemes::theme_solarized(),
      ggthemes::theme_gdocs()
    )
  ),
  .f = ggwithinstats
)

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(title = "Response times across flank conditions for each type of cue"),
  plotgrid.args = list(ncol = 1)
)
```

## `ggscatterstats`

For the next example lets use the same methodology on different data and using
`ggscatterstats` to produce scatterplots combined with marginal
histograms/boxplots/density plots with statistical details added as a subtitle.
For data we'll use `movies_long` which is from IMDB and part of the
`{ggstatsplot}` package. Since it's a large dataset with some relatively small
categories like **NC-17** we'll sample only one quarter of the data and
completely drop NC-17 using `dplyr`. 

This time we'll put all the code in one block-

```r
## for reproducibility
set.seed(123)

mpaa_list <- movies_long %>%
  dplyr::filter(mpaa != "NC-17") %>%
  dplyr::sample_frac(size = 0.25) %>%
  split(f = .$mpaa, drop = TRUE)

## creating a list of plots
plot_list <- purrr::pmap(
  .l = list(
    data = mpaa_list,
    x = "budget",
    y = "rating",
    xlab = "Budget (in millions of US dollars)",
    ylab = "Rating on IMDB",
    title = list(
      "MPAA Rating: PG",
      "MPAA Rating: PG-13",
      "MPAA Rating: R"
    ),
    label.var = list("title"),
    ## note that you need to quote the expressions
    label.expression = list(
      quote(rating > 7.5 & budget < 100),
      quote(rating > 8 & budget < 50),
      quote(rating > 8 & budget < 10)
    ),
    type = list("r", "np", "bf"),
    xfill = list("#009E73", "#999999", "#0072B2"),
    yfill = list("#CC79A7", "#F0E442", "#D55E00"),
    ggtheme = list(
      ggthemes::theme_tufte(),
      ggplot2::theme_classic(),
      ggplot2::theme_light()
    )
  ),
  .f = ggscatterstats
)

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(
    title = "Relationship between movie budget and IMDB rating",
    caption = "Source: www.imdb.com"
  ),
  plotgrid.args = list(ncol = 1)
)
```

The remainder of the examples vary in content but follow the exact same
methodology as the earlier examples.

## `ggcorrmat`

```r
## for reproducibility
set.seed(123)

## splitting the dataframe by cut and creating a list
## let's leave out "fair" cut
## also, to make this fast, let's only use 5% of the sample
cut_list <- ggplot2::diamonds %>%
  dplyr::sample_frac(size = 0.05) %>%
  dplyr::filter(cut != "Fair") %>%
  split(f = .$cut, drop = TRUE)

## checking the length and names of each element
length(cut_list)
names(cut_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <- purrr::pmap(
  .l = list(
    data = cut_list,
    cor.vars = list(c("carat", "depth", "table", "price")),
    type = list("pearson", "np", "robust", "bf"),
    partial = list(TRUE, FALSE, TRUE, FALSE),
    title = list("Cut: Good", "Cut: Very Good", "Cut: Premium", "Cut: Ideal"),
    p.adjust.method = list("hommel", "fdr", "BY", "hochberg"),
    lab.size = 3.5,
    colors = list(
      c("#56B4E9", "white", "#999999"),
      c("#CC79A7", "white", "#F0E442"),
      c("#56B4E9", "white", "#D55E00"),
      c("#999999", "white", "#0072B2")
    ),
    ggtheme = list(
      ggplot2::theme_linedraw(),
      ggplot2::theme_classic(),
      ggthemes::theme_fivethirtyeight(),
      ggthemes::theme_tufte()
    )
  ),
  .f = ggcorrmat
)

## combining all individual plots from the list into a single plot using
## `combine_plots` function
combine_plots(
  plotlist = plot_list,
  guides = "keep",
  annotation.args = list(
    title = "Relationship between diamond attributes and price across cut",
    caption = "Dataset: Diamonds from ggplot2 package"
  ),
  plotgrid.args = list(nrow = 2)
)
```

## `gghistostats`

```r
## for reproducibility
set.seed(123)

## let's split the dataframe and create a list by continent
## let's leave out Oceania because it has just two data points
continent_list <-
  gapminder::gapminder %>%
  dplyr::filter(year == 2007, continent != "Oceania") %>%
  split(f = .$continent, drop = TRUE)

## checking the length and names of each element
length(continent_list)
names(continent_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <-
  purrr::pmap(
    .l = list(
      data = continent_list,
      x = "lifeExp",
      xlab = "Life expectancy",
      test.value = list(35.6, 58.4, 41.6, 64.7),
      type = list("p", "np", "r", "bf"),
      bf.message = list(TRUE, FALSE, FALSE, FALSE),
      title = list(
        "Continent: Africa",
        "Continent: Americas",
        "Continent: Asia",
        "Continent: Europe"
      ),
      effsize.type = list("d", "d", "g", "g"),
      normal.curve = list(TRUE, FALSE, FALSE, TRUE),
      ggtheme = list(
        ggplot2::theme_classic(),
        hrbrthemes::theme_ipsum_tw(),
        ggplot2::theme_minimal(),
        hrbrthemes::theme_modern_rc()
      )
    ),
    .f = gghistostats
  )

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(
    title = "Improvement in life expectancy worldwide since 1950",
    caption = "Note: black line - 1950; blue line - 2007"
  ),
  plotgrid.args = list(nrow = 4)
)
```  

## `ggdotplotstats`

```r
## for reproducibility
set.seed(123)
library(ggthemes)
library(hrbrthemes)

## let's split the dataframe and create a list by continent
## let's leave out Oceania because it has just two data points
continent_list <-
  gapminder::gapminder %>%
  dplyr::filter(continent != "Oceania") %>%
  split(f = .$continent, drop = TRUE)

## checking the length and names of each element
length(continent_list)
names(continent_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <-
  purrr::pmap(
    .l = list(
      data = continent_list,
      x = "gdpPercap",
      y = "year",
      xlab = "GDP per capita (US$, inflation-adjusted)",
      test.value = list(2500, 9000, 9500, 10000),
      type = list("p", "np", "r", "bf"),
      title = list(
        "Continent: Africa",
        "Continent: Americas",
        "Continent: Asia",
        "Continent: Europe"
      ),
      effsize.type = list("d", "d", "g", "g"),
      centrality.line.args = list(
        list(color = "red"),
        list(color = "#0072B2"),
        list(color = "#D55E00"),
        list(color = "#CC79A7")
      ),
      ggtheme = list(
        ggplot2::theme_minimal(base_family = "serif"),
        ggthemes::theme_tufte(),
        hrbrthemes::theme_ipsum_rc(axis_title_size = 10),
        ggthemes::theme_hc(bgcolor = "darkunica")
      )
    ),
    .f = ggdotplotstats
  )

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(title = "Improvement in GDP per capita from 1952-2007"),
  plotgrid.args = list(nrow = 4),
  guides = "keep"
)
```  

## `ggpiestats`

```r
## for reproducibility
set.seed(123)

## let's split the dataframe and create a list by passenger class
class_list <- Titanic_full %>% split(f = .$Class, drop = TRUE)

## checking the length and names of each element
length(class_list)
names(class_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <-
  purrr::pmap(
    .l = list(
      data = class_list,
      x = "Survived",
      y = "Sex",
      label = list("both", "count", "percentage", "both"),
      title = list(
        "Passenger class: 1st",
        "Passenger class: 2nd",
        "Passenger class: 3rd",
        "Passenger class: Crew"
      ),
      caption = list(
        "Total: 319, Died: 120, Survived: 199, % Survived: 62%",
        "Total: 272, Died: 155, Survived: 117, % Survived: 43%",
        "Total: 709, Died: 537, Survived: 172, % Survived: 25%",
        "Data not available for crew passengers"
      ),
      package = list("RColorBrewer", "ghibli", "palettetown", "yarrr"),
      palette = list("Accent", "MarnieMedium1", "pikachu", "nemo"),
      ggtheme = list(
        ggplot2::theme_grey(),
        ggplot2::theme_bw(),
        ggthemes::theme_tufte(),
        ggthemes::theme_economist()
      ),
      proportion.test = list(TRUE, FALSE, TRUE, FALSE),
      type = list("p", "p", "bf", "p")
    ),
    .f = ggpiestats
  )

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(title = "Survival in Titanic disaster by gender for all passenger classes"),
  plotgrid.args = list(ncol = 1),
  guides = "keep"
)
``` 

## `ggbarstats`

```r
## for reproducibility
set.seed(123)

## let's split the dataframe and create a list by passenger class
class_list <- Titanic_full %>% split(f = .$Class, drop = TRUE)

## checking the length and names of each element
length(class_list)
names(class_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <-
  purrr::pmap(
    .l = list(
      data = class_list,
      x = "Survived",
      y = "Sex",
      type = "bayes",
      label = list("both", "count", "percentage", "both"),
      title = list(
        "Passenger class: 1st",
        "Passenger class: 2nd",
        "Passenger class: 3rd",
        "Passenger class: Crew"
      ),
      caption = list(
        "Total: 319, Died: 120, Survived: 199, % Survived: 62%",
        "Total: 272, Died: 155, Survived: 117, % Survived: 43%",
        "Total: 709, Died: 537, Survived: 172, % Survived: 25%",
        "Data not available for crew passengers"
      ),
      package = list("RColorBrewer", "ghibli", "palettetown", "yarrr"),
      palette = list("Accent", "MarnieMedium1", "pikachu", "nemo"),
      ggtheme = list(
        ggplot2::theme_grey(),
        ggplot2::theme_bw(),
        ggthemes::theme_tufte(),
        ggthemes::theme_economist()
      )
    ),
    .f = ggbarstats
  )

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(
    title = "Survival in Titanic disaster by gender for all passenger classes",
    caption = "Asterisks denote results from proportion tests: \n***: p < 0.001, ns: non-significant"
  ),
  plotgrid.args = list(ncol = 1),
  guides = "keep"
)
``` 

## `grouped_` variants

Note that although all the above examples were written with the non-grouped
variants of functions, the same rule holds true for the `grouped_` variants of
all the above functions.

For example, if we want to use the `grouped_gghistostats` across three different
datasets, you can use `purrr::pmap()` function. For the sake of brevity, the
plots are not displayed here, but you can run the following code and check the
individual `grouped_` plots (e.g., `plotlist[[1]]`).

```r
## create a list of plots
plotlist <- purrr::pmap(
    .l = list(
      data = list(mtcars, iris, ToothGrowth),
      x = alist(wt, Sepal.Length, len),
      results.subtitle = list(FALSE),
      grouping.var = alist(am, Species, supp)
    ),
    .f = grouped_gghistostats
  )

## given that we had three different datasets, we expect a list of length 3
## (each of which contains a `grouped_` plot)
length(plotlist)
```

## Repeating function execution across multiple columns in a dataframe

```r
## setup
set.seed(123)
library(ggstatsplot)
library(patchwork)

## running the same analysis on two different columns (creates a list of plots)
plotlist <- purrr::pmap(
    .l = list(
      data = list(movies_long),
      x = "mpaa",
      y = list("rating", "length"),
      title = list("IMDB score by MPAA rating", "Movie length by MPAA rating")
    ),
    .f = ggbetweenstats
  )

## combine plots using `patchwork`
plotlist[[1]] + plotlist[[2]]
```

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "ggcorrmat"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{ggcorrmat}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
## pretty tibble printing
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>"
)

library(ggstatsplot)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

Lifecycle: [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

The function `ggcorrmat` provides a quick way to produce **publication-ready correlation matrix** (aka *correlalogram*) plot. The function can also be used
for quick **data exploration**. In addition to the plot, it can also be used to
get a correlation coefficient matrix or the associated *p*-value matrix.
Currently, the plot can display Pearson's *r* (and its Bayesian version),
Spearman's *rho*, and *robust* correlation coefficient (Winsorized Pearson's
*r*). This function is a convenient wrapper around `ggcorrplot::ggcorrplot`
function with some additional functionality.

We will see examples of how to use this function in this vignette with the
`gapminder` and `diamonds` dataset.

To begin with, here are some instances where you would want to use
`ggcorrmat`-

  - to easily visualize a correlation matrix using `{ggplot2}`
  - to quickly explore correlation between (all) numeric variables in the
  dataset

## Correlation matrix plot with `ggcorrmat`

For the first example, we will use the `gapminder` dataset (available in
eponymous [package](https://CRAN.R-project.org/package=gapminder) on CRAN)
provides values for life expectancy, Gross Domestic Product (GDP) per capita,
and population, every five years, from 1952 to 2007, for each of 142 countries
and was collected by the Gapminder Foundation. Let's have a look at the data-

```{r gapminder}
library(gapminder)
library(dplyr)

dplyr::glimpse(gapminder)
```

Let's say we are interested in studying correlation between population of a
country, average life expectancy, and GDP per capita across countries only for
the year 2007.

The simplest way to get a correlation matrix is to stick to the defaults-

```{r ggcorrmat1, fig.height = 6, fig.width = 6}
## setup
set.seed(123)
library(ggstatsplot)

## select data only from the year 2007
gapminder_2007 <- dplyr::filter(gapminder::gapminder, year == 2007)

## producing the correlation matrix
ggcorrmat(
  data = gapminder_2007, ## data from which variable is to be taken
  cor.vars = lifeExp:gdpPercap ## specifying correlation matrix variables
)
```

This plot can be further modified with additional arguments-

```{r ggcorrmat2, fig.height = 6, fig.width = 6}
ggcorrmat(
  data = gapminder_2007, ## data from which variable is to be taken
  cor.vars = lifeExp:gdpPercap, ## specifying correlation matrix variables
  cor.vars.names = c(
    "Life Expectancy",
    "population",
    "GDP (per capita)"
  ),
  type = "spearman", ## which correlation coefficient is to be computed
  lab.col = "red", ## label color
  ggtheme = ggplot2::theme_light(), ## selected ggplot2 theme
  ## turn off default ggestatsplot theme overlay
  matrix.type = "lower", ## correlation matrix structure
  colors = NULL, ## turning off manual specification of colors
  palette = "category10_d3", ## choosing a color palette
  package = "ggsci", ## package to which color palette belongs
  title = "Gapminder correlation matrix", ## custom title
  subtitle = "Source: Gapminder Foundation" ## custom subtitle
)
```

As seen from this correlation matrix, although there is no relationship between
population and life expectancy worldwide, at least in 2007, there is a strong
positive relationship between GDP, a well-established indicator of a country's
economic performance.

Given that there were only three variables, this doesn't look that impressive.
So let's work with another example from `{ggplot2}` package: the `diamonds`
[dataset](http://ggplot2.tidyverse.org/reference/diamonds.html). This dataset
contains the prices and other attributes of almost 54,000 diamonds.

Let's have a look at the data-

```{r diamonds}
library(ggplot2)

dplyr::glimpse(ggplot2::diamonds)
```

Let's see the correlation matrix between different attributes of the diamond and
the price.

```{r ggcorrmat3, fig.height = 7, fig.width = 7}
## for reproducibility
set.seed(123)

## let's use just 5% of the data to speed it up
ggcorrmat(
  data = dplyr::sample_frac(ggplot2::diamonds, size = 0.05),
  cor.vars = c(carat, depth:z), ## note how the variables are getting selected
  cor.vars.names = c(
    "carat",
    "total depth",
    "table",
    "price",
    "length (in mm)",
    "width (in mm)",
    "depth (in mm)"
  ),
  ggcorrplot.args = list(outline.color = "black", hc.order = TRUE)
)
```

We can make a number of changes to this basic correlation matrix. For example,
since we were interested in relationship between price and other attributes,
let's make the `price` column to the the first column.

```{r ggcorrmat4, fig.height = 7, fig.width = 7}
## for reproducibility
set.seed(123)

## let's use just 5% of the data to speed it up
ggcorrmat(
  data = dplyr::sample_frac(ggplot2::diamonds, size = 0.05),
  cor.vars = c(price, carat, depth:table, x:z), ## note how the variables are getting selected
  cor.vars.names = c(
    "price",
    "carat",
    "total depth",
    "table",
    "length (in mm)",
    "width (in mm)",
    "depth (in mm)"
  ),
  type = "spearman",
  title = "Relationship between diamond attributes and price",
  subtitle = "Dataset: Diamonds from ggplot2 package",
  colors = c("#0072B2", "#D55E00", "#CC79A7"),
  pch = "square cross",
  ## additional aesthetic arguments passed to `ggcorrmat`
  ggcorrplot.args = list(
    lab_col = "yellow",
    lab_size = 6,
    tl.srt = 90,
    pch.col = "white",
    pch.cex = 14
  )
) + ## modification outside `{ggstatsplot}` using `{ggplot2}` functions
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 0.15, r = 0.15, b = 0.15, l = 0.15, unit = "cm")
    )
  )
```

As seen here, and unsurprisingly, the strongest predictor of the diamond price
is its carat value, which a unit of mass equal to 200 mg. In other words, the
heavier the diamond, the more expensive it is going to be.

## Dataframe containing statistics with `ggcorrmat`

Another utility of `ggcorrmat` is in obtaining a dataframe containing all
details from statistical analyses. Such dataframes can be easily embedded in
manuscripts as tables.

```{r ggcorrmat5}
## for reproducibility
set.seed(123)

## to get correlations
ggcorrmat(
  data = dplyr::sample_frac(ggplot2::txhousing, size = 0.15),
  cor.vars = sales:inventory,
  output = "dataframe",
  type = "robust",
  digits = 3
)
```

Note that if `cor.vars` are not specified, all numeric variables will be used.
Moreover, you can also use abbreviations to specify what output you want in
return. Additionally, one can also carry out **partial** correlation analysis:

```{r ggcorrmat6}
## for reproducibility
set.seed(123)
options(pillar.sigfig = 4)

## getting the correlation coefficient matrix
ggcorrmat(
  data = iris, ## all numeric variables from data will be used
  type = "np", ## non-parametric
  partial = TRUE,
  output = "dataframe"
)
```

## Grouped analysis with `grouped_ggcorrmat`

What if we want to do the same analysis separately for each quality of the
diamond `cut` (Fair, Good, Very Good, Premium, Ideal)? 

`{ggstatsplot}` provides a special helper function for such instances:
`grouped_ggcorrmat`. This is merely a wrapper function around
`combine_plots`. It applies `ggcorrmat` across all **levels** of
a specified **grouping variable** and then combines list of individual plots
into a single plot. Note that the grouping variable can be anything: conditions
in a given study, groups in a study sample, different studies, etc. 

```{r ggcorrmat7, fig.height = 16, fig.width = 10}
## for reproducibility
set.seed(123)

## plot
grouped_ggcorrmat(
  ## arguments relevant for `ggcorrmat`
  data = ggplot2::diamonds,
  type = "bayes", ## Bayesian test
  grouping.var = cut,
  ## arguments relevant for `combine_plots`
  plotgrid.args = list(nrow = 3),
  annotation.args = list(
    tag_levels = "a",
    title = "Relationship between diamond attributes and price across cut",
    caption = "Dataset: Diamonds from ggplot2 package"
  )
)
```

Note that this function also makes it easy to run the same correlation matrix
across different levels of a factor/grouping variable.

```{r ggcorrmat8}
## for reproducibility
set.seed(123)

## let's obtain correlation coefficients with their CIs
grouped_ggcorrmat(
  data = ggplot2::msleep,
  cor.vars = sleep_total:awake,
  grouping.var = vore,
  output = "dataframe"
)
```

## Grouped analysis with `ggcorrmat` + `{purrr}`

Although `grouped_` function is good for quickly exploring the data, it reduces
the flexibility with which this function can be used. This is the because the
common parameters used are applied to plots corresponding to all levels of the
grouping variable and there is no way to customize the arguments for different
levels of the grouping variable. We will see how this can be done using the
`{purrr}` package.  

See the associated vignette here:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html>

## Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
correlation matrix | `ggcorrplot::ggcorrplot` | `ggcorrplot.args`

## Summary of tests

**Hypothesis testing** and **Effect size estimation**

Type | Test | CI? | Function used
----------- | ------------------------- | --- | -----
Parametric | Pearson's correlation coefficient | ✅ | `correlation::correlation`
Non-parametric | Spearman's rank correlation coefficient | ✅ | `correlation::correlation`
Robust | Winsorized Pearson correlation coefficient | ✅ | `correlation::correlation`
Bayesian | Pearson's correlation coefficient | ✅ | `correlation::correlation`

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "ggdotplotstats"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{ggdotplotstats}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
## show me all columns
options(tibble.width = Inf, pillar.bold = TRUE, pillar.subtle_num = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>"
)

library(ggstatsplot)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

Lifecycle: [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

The function `ggdotplotstats` can be used for **data exploration** and to
provide an easy way to make **publication-ready dot plots/charts** with
appropriate and selected statistical details embedded in the plot itself. In
this vignette, we will explore several examples of how to use it.

This function is a sister function of `gghistostats` with the difference being
it expects a labeled numeric variable.

## Distribution of a sample with `ggdotplotstats`

Let's begin with a very simple example from the `{ggplot2}` package
(`ggplot2::mpg`), a subset of the fuel economy data that the EPA makes available
on <http://fueleconomy.gov>.

```{r mpg}
## looking at the structure of the data using glimpse
dplyr::glimpse(ggplot2::mpg)
```

Let's say we want to visualize the distribution of mileage by car manufacturer. 

```{r mpg2, fig.height = 7, fig.width = 9}
## for reproducibility
set.seed(123)
library(ggstatsplot)

## removing factor level with very few no. of observations
df <- dplyr::filter(ggplot2::mpg, cyl %in% c("4", "6"))

## creating a vector of colors using `paletteer` package
paletter_vector <-
  paletteer::paletteer_d(
    palette = "palettetown::venusaur",
    n = nlevels(as.factor(df$manufacturer)),
    type = "discrete"
  )

## plot
ggdotplotstats(
  data = df,
  x = cty,
  y = manufacturer,
  xlab = "city miles per gallon",
  ylab = "car manufacturer",
  test.value = 15.5,
  point.args = list(
    shape = 16,
    color = paletter_vector,
    size = 5
  ),
  title = "Distribution of mileage of cars",
  ggtheme = ggplot2::theme_dark()
)
```

## Grouped analysis with `grouped_ggdotplotstats`

What if we want to do the same analysis separately for different engines with
different numbers of cylinders?

`{ggstatsplot}` provides a special helper function for such instances:
`grouped_ggdotplotstats`. This is merely a wrapper function around
`combine_plots`. It applies `ggdotplotstats` across all **levels** of
a specified **grouping variable** and then combines the individual plots into a
single plot. 

Let's see how we can use this function to apply `ggdotplotstats` to accomplish our
task. 

```{r grouped1, fig.height = 12, fig.width = 7}
## for reproducibility
set.seed(123)

## removing factor level with very few no. of observations
df <- dplyr::filter(ggplot2::mpg, cyl %in% c("4", "6"))

## plot
grouped_ggdotplotstats(
  ## arguments relevant for ggdotplotstats
  data = df,
  grouping.var = cyl, ## grouping variable
  x = cty,
  y = manufacturer,
  xlab = "city miles per gallon",
  ylab = "car manufacturer",
  type = "bayes", ## Bayesian test
  test.value = 15.5,
  ## arguments relevant for `combine_plots`
  annotation.args = list(title = "Fuel economy data"),
  plotgrid.args = list(nrow = 2)
)
```

## Grouped analysis with `{purrr}`

Although this is a quick and dirty way to explore a large amount of data with
minimal effort, it does come with an important limitation: reduced flexibility.
For example, if we wanted to add, let's say, a separate `test.value` argument
for each gender, this is not possible with `grouped_ggdotplotstats`. For cases
like these, or to run separate  kinds of tests (robust for some, parametric for
other, while Bayesian for some other levels of the group) it would be better to
use `{purrr}`.

See the associated vignette here:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html>

## Summary of tests

**Central tendency measure**

Type | Measure | Function used
----------- | --------- | ------------------ 
Parametric | mean | `parameters::describe_distribution`
Non-parametric | median | `parameters::describe_distribution`
Robust | trimmed mean | `parameters::describe_distribution`
Bayesian | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution`

**Hypothesis testing**

Type | Test | Function used
------------------ | ------------------------- | -----
Parametric | One-sample Student's *t*-test | `stats::t.test`
Non-parametric | One-sample Wilcoxon test | `stats::wilcox.test`
Robust | Bootstrap-*t* method for one-sample test | `WRS2::trimcibt`
Bayesian | One-sample Student's *t*-test | `BayesFactor::ttestBF`

**Effect size estimation**

Type | Effect size | CI? | Function used
------------ | ----------------------- | --- | -----
Parametric | Cohen's *d*, Hedge's *g* | ✅ | `effectsize::cohens_d`, `effectsize::hedges_g`
Non-parametric | *r* (rank-biserial correlation) | ✅ | `effectsize::rank_biserial`
Robust | trimmed mean | ✅ | `WRS2::trimcibt`
Bayes Factor | $\delta_{posterior}$ | ✅ | `bayestestR::describe_posterior`

## Reporting

If you wish to include statistical analysis results in a publication/report, the
ideal reporting practice will be a hybrid of two approaches:

- the `{ggstatsplot}` approach, where the plot contains both the visual and
numerical summaries about a statistical model, and

- the *standard* narrative approach, which provides interpretive context for the
reported statistics.

For example, let's see the following example:

```{r reporting}
ggdotplotstats(morley, Speed, Expt, test.value = 800)
```

The narrative context (assuming `type = "parametric"`) can complement this plot
either as a figure caption or in the main text-

> Student's *t*-test revealed that, across 5 experiments, the speed of light was
significantly different than posited speed. The effect size $(g = 1.22)$ was
very large, as per Cohen’s (1988) conventions. The Bayes Factor for the same
analysis revealed that the data were `r round(exp(1.24), 2)` times more probable
under the alternative hypothesis as compared to the null hypothesis. This can be
considered moderate evidence (Jeffreys, 1961) in favor of the alternative
hypothesis.


## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on GitHub: 
<https://github.com/IndrajeetPatil/ggstatsplot/issues>

---
title: "frequently asked questions (FAQ)"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{frequently asked questions (FAQ)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
## show me all columns
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>"
)

if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(PMCMRplus)
}

library(ggstatsplot)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment=""}
citation("ggstatsplot")
```

<!-- The sections are numbered manually because `pkgdown` doesn't render -->
<!-- numbered sections for `rmarkdown::html_vignette` format -->

---

Following are a few of the common questions asked in GitHub issues and on social
media platforms.

## 1. I just want the plot, not the statistical details. How can I turn them off?

All functions in `{ggstatsplot}` that display results from statistical analysis in
a subtitle have argument `results.subtitle`. Setting it to `FALSE` will return
only the plot.

## 2. How can I customize the details contained in the subtitle?

Sometimes you may not wish include so many details in the subtitle. In that
case, you can extract the expression and copy-paste only the part you wish to
include. For example, here only statistic and *p*-values are included:

```{r custom_expr}
## setup
set.seed(123)
library(ggstatsplot)
library(ggplot2)
library(statsExpressions)

## extracting detailed expression
(res_expr <- oneway_anova(iris, Species, Sepal.Length, var.equal = TRUE))

## adapting the details to your liking
ggplot(iris, aes(x = Species, y = Sepal.Length)) +
  geom_boxplot() +
  labs(subtitle = ggplot2::expr(paste(
    NULL, italic("F"), "(", "2",
    ",", "147", ")=", "119.26", ", ",
    italic("p"), "=", "1.67e-31"
  )))
```

## 3. I am getting `Error in grid.Call` error

Sometimes, if you are working in `RStudio`, you might see the following error-

```r
Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
polygon edge not found
```

This can possibly be solved by increasing the size of RStudio viewer pane.

## 4. Why do I get only plot in return but not the subtitle/caption?

In order to prevent the function from failing when statistical analysis fails,
functions in `{ggstatsplot}` default to first attempting to run the analysis and
if they fail, then return empty (`NULL`) subtitle/caption. In such cases, if you
wish to diagnose why the analysis is failing, you will have to do so using the
underlying function used to carry out statistical analysis.

For example, the following returns only the plot but not the statistical details
in a subtitle.

```{r null_subtitle, fig.width=3, fig.height=3}
set.seed(123)
df <- data.frame(x = c("a", "b"), y = c(1, 2))

ggbetweenstats(data = df, x = x, y = y)
```

To see why the statistical analysis failed, you can look at the error from the
underlying function:

```{r, error=TRUE}
library(statsExpressions)

df <- data.frame(x = c("a", "b"), y = c(1, 2))
two_sample_test(data = df, x = x, y = y)
```

## 5. What statistical test was carried out?

In case you are not sure what was the statistical test that produced the results
shown in the subtitle of the plot, the best way to get that information is to
either look at the documentation for the function used or check out the
associated vignette. 

Summary of all analysis is handily available in `README`: 
<https://github.com/IndrajeetPatil/ggstatsplot/blob/master/README.md>

## 6. How can I use `{ggstatsplot}` functions in a `for` loop?

Given that all functions in `{ggstatsplot}` use tidy evaluation, running these
functions in a `for` loop requires minor adjustment to how inputs are entered:

```{r loop, eval=FALSE}
## setup
data(mtcars)

col.name <- colnames(mtcars)

## executing the function in a `for` loop
for (i in 3:length(col.name)) {
  ggbetweenstats(
    data = mtcars,
    x = cyl,
    y = !!col.name[i]
  )
}
```

That said, if repeating function execution across multiple columns in a
dataframe in what you want to do, I will recommend you to have a look at
`{purrr}`-based solution:

<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html#repeating-function-execution-across-multiple-columns-in-a-dataframe-1>

## 7. How can I have uniform Y-axes ranges in `grouped_` functions?

Across different facets of a `grouped_` plot, the axes ranges might sometimes
differ. You can use the `ggplot.component` (present in all functions) to have
the same scale across the individual plots:

```{r grouped_y_axes, fig.height=6, fig.width=8}
## setup
set.seed(123)


## provide a list of further `{ggplot2}` modifications using `ggplot.component`
grouped_ggscatterstats(
  mtcars,
  disp,
  hp,
  grouping.var = am,
  results.subtitle = FALSE,
  ggplot.component = list(ggplot2::scale_y_continuous(
    breaks = seq(50, 350, 50),
    limits = (c(50, 350))
  ))
)
```

## 8. Does `{ggstatsplot}` work with `plotly`?

The `plotly` R graphing library makes it easy to produce interactive web
graphics via `plotly.js`. 

The `{ggstatsplot}` functions are compatible with `plotly`.
 
```r
## for reproducibility
set.seed(123)
library(plotly)

## creating ggplot object with `{ggstatsplot}`
p <- ggbetweenstats(mtcars, cyl, mpg)

## converting to plotly object
plotly::ggplotly(p, width = 480, height = 480)
```

## 9. How can I use `grouped_` functions with more than one group?

Currently, the `grouped_` variants of functions only support repeating the
analysis across a _single_ grouping variable. Often, you have to run the same
analysis across a combination of more than two grouping variables. This can be
easily achieved using `{purrr}` package. 

Here is an example-

```{r grouped_2, fig.width=6, fig.height=6}
## setup
set.seed(123)

## creating a list by splitting dataframe by combination of two different
## grouping variables
df_list <-
  mpg %>%
  dplyr::filter(drv %in% c("4", "f"), fl %in% c("p", "r")) %>%
  split(f = list(.$drv, .$fl), drop = TRUE)

## checking if the length of the list is 4
length(df_list)

## running correlation analyses between
## this will return a *list* of plots
plot_list <- purrr::pmap(
    .l = list(
      data = df_list,
      x = "displ",
      y = "hwy",
      results.subtitle = FALSE
    ),
    .f = ggscatterstats
  )

## arragen the list in a single plot
combine_plots(
  plotlist = plot_list,
  plotgrid.args = list(nrow = 2),
  annotation.args = list(tag_levels = "i")
)
```

## 10. How can I include statistical expressions in facet labels?

```{r facet_expr, fig.width=6, fig.height=8}
set.seed(123)
library(ggplot2)

## data
mtcars1 <- mtcars
statistics <-
  grouped_ggbetweenstats(
    data = mtcars1,
    x = cyl,
    y = mpg,
    grouping.var = am,
    output = "subtitle"
  )
mtcars1$am <- factor(mtcars1$am, levels = c(0, 1), labels = statistics)

## plot
mtcars1 %>%
  ggplot(aes(x = cyl, y = mpg)) +
  geom_jitter() +
  facet_wrap(
    vars(am),
    ncol = 1,
    strip.position = "top",
    labeller = ggplot2::label_parsed
  )
```

## 11. Can you customize which pairs are shown in pairwise comparisons?

Currently, for `ggbetweenstats` and `ggwithinstats`, you can either display all
**significant** comparisons, all **non-significant** comparisons, or **all**
comparisons. But what if I am only interested in just one particular comparison?

Here is a workaround using `ggsignif`:

```{r custom_pairwise, fig.width=7, fig.height=6}
set.seed(123)
library(ggsignif)

## displaying only one comparison
ggbetweenstats(mtcars, cyl, wt, pairwise.comparisons = FALSE) +
  geom_signif(comparisons = list(c("4", "6")))
```

## 12. How to access dataframe with results from pairwise comparisons?

Behind the scenes, `{ggstatsplot}` uses `pairwise_comparisons` function. You
can use it to extract actual dataframes used in `{ggstatsplot}` functions.

```{r}
library(ggplot2)

pairwise_comparisons(mtcars, cyl, wt)
```

## 13. How can I change annotation in pairwise comparisons?

`{ggstatsplot}` defaults to displaying exact p-values or logged Bayes Factor
values for pairwise comparisons. But what if you wish to adopt a different
annotation labels?

You will have to customize them using `pairwiseComparisons` and `ggsignif:`

```{r comp_asterisks}
## needed libraries
set.seed(123)
library(ggplot2)
library(ggsignif)

## converting to factor
mtcars$cyl <- as.factor(mtcars$cyl)

## creating a basic plot
p <- ggbetweenstats(mtcars, cyl, wt, pairwise.comparisons = FALSE)

## using `pairwise_comparisons()` function to create a dataframe with results
set.seed(123)
(df <-
  pairwise_comparisons(mtcars, cyl, wt) %>%
  dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
  dplyr::arrange(group1) %>%
  dplyr::mutate(asterisk_label = c("**", "***", "**")))

## adding pairwise comparisons using `ggsignif`
p +
  ggsignif::geom_signif(
    comparisons = df$groups,
    map_signif_level = TRUE,
    annotations = df$asterisk_label,
    y_position = c(5.5, 5.75, 6.0),
    test = NULL,
    na.rm = TRUE
  )
```

## 14. How to access dataframe with results from `ggpiestats` and `ggbarstats`?

```{r onesample_df}
## setup
set.seed(123)
library(ggplot2)

## plot
p <- ggpiestats(mtcars, am, cyl)

## dataframe with results
extract_stats(p)
```

## 15. How can I remove a particular `geom` layer from the plot?

Sometimes you may not want a particular `geom` layer to be displayed. You can
remove them using `gginnards`. 

For example, let's say we want to remove the `geom_point()` from
`ggwithinstats` default plot.

```{r gginnards, fig.width=7, fig.height=5}
## needed libraries
library(gginnards)

## plot with all geoms
p <- ggwithinstats(
  data = bugs_long,
  x = condition,
  y = desire,
  results.subtitle = FALSE,
  pairwise.comparisons = FALSE
)

## delete `geom` corresponding to violin
gginnards::delete_layers(x = p, match_type = "GeomViolin")
```

This can be helpful to add a new layer with aesthetic specifications of your
liking.

```{r gginnards2, fig.width=7, fig.height=5}
## needed libraries
set.seed(123)
library(gginnards)
library(ggplot2)

## basic plot without mean tagging
p <- ggbetweenstats(
  data = mtcars,
  x = am,
  y = wt,
  centrality.plotting = FALSE
)

## delete the geom_point layer
p <- gginnards::delete_layers(x = p, match_type = "GeomPoint")

## add a new layers for points with a different shape
p + geom_point(shape = 23, aes(color = am))
```

## 16. How can I modify the fill colors with custom values?

Sometimes you may not be satisfied with the available color palette values. In
this case, you can also change the colors by manually specifying these values.

```{r ggbar_colors, fig.width=5, fig.height=5}
## needed libraries
set.seed(123)

library(ggplot2)

ggbarstats(mtcars, am, cyl, results.subtitle = FALSE) +
  scale_fill_manual(values = c("#E7298A", "#66A61E"))
```

The same can also be done for `grouped_` functions:

```{r ggpie_colors, fig.width=12, fig.height=6}
grouped_ggpiestats(
  data = mtcars,
  grouping.var = am,
  x = cyl,
  ggplot.component = ggplot2::scale_fill_grey()
)
```

## 17. How can I modify `grouped_` outputs using `{ggplot2}` functions?

All `{ggstatsplot}` are `ggplot` objects, which can be further modified, just like
any other `ggplot` object. But exception to these are all plots returned by
`grouped_` functions, but there is a way to tackle this.

```{r grouped_modify, fig.width=12, fig.height=6}
## needed libraries
set.seed(123)
library(paletteer)
library(ggplot2)

## plot
grouped_ggbetweenstats(
  mtcars,
  cyl,
  wt,
  grouping.var = vs,
  type = "np",
  ggplot.component =
  ## modify further with `{ggplot2}` functions
    list(
      scale_color_manual(values = paletteer::paletteer_c("viridis::viridis", 3)),
      theme(axis.text.x = element_text(angle = 90))
    )
)
```

## 18. How can I extract dataframe containing results from `{ggstatsplot}`?

`{ggstatsplot}` can return expressions in the subtitle and caption, but what if
you want to actually get back dataframe containing the results?

This is possible via `{statsExpressions}`:
<https://indrajeetpatil.github.io/statsExpressions/articles/dataframe_outputs.html>

## 19. How can I remove sample size labels for `ggbarstats`?

```{r ggbar_samplesize}
library(gginnards)

## create a plot
p <- ggbarstats(mtcars, am, cyl)

## remove layer corresponding to sample size
delete_layers(p, "GeomText")
```

## 20. Test I need is not available. What can I do?

By default, since `{ggstatsplot}` always allows just **one** type of test per
statistical approach, sometimes your favorite test might not be available. For
example, `{ggstatsplot}` provides only Spearman's $\rho$, but not Kendall's
$\tau$ as a non-parametric correlation test. 

In such cases, you can override the defaults and use `{statsExpressions}` to
create custom expressions to display in the plot. But be forewarned that the
expression building functions in `{statsExpressions}` are not stable yet.

```{r custom_test, fig.width=6, fig.height=6}
## setup
set.seed(123)
library(correlation)
library(statsExpressions)
library(ggplot2)

## data with two variables of interest
df <- dplyr::select(mtcars, wt, mpg)

## correlation results
results <-
  correlation(df, method = "kendall") %>%
  parameters::standardize_names(style = "broom")

## creating expression out of these results
df_results <- statsExpressions::add_expression_col(
    data           = results,
    no.parameters  = 0L,
    statistic.text = list(quote(italic("T"))),
    effsize.text   = list(quote(widehat(italic(tau))["Kendall"])),
    n              = results$n.obs[[1]]
  )

## plot (overriding defaults and using custom expression)
ggscatterstats(
  df, wt, mpg,
  results.subtitle = FALSE,
  ggplot.component = list(labs(subtitle = df_results$expression[[1]]))
)
```

## 21. Is there way to adjust my alpha level?

No, there is no way to adjust alpha if you use `grouped_` functions (e.g.,
`grouped_ggwithinstats`). You will have to just report in the
paper/article/report, what your adjusted alpha is.

So, for example, iif 2 tests are being carried out, the alpha is going to be
`0.05/2 = 0.025`. So, when you describe the _Methods_ section, you can mention
that only those tests should be considered significant where `p < 0.025`. Or you
can even mention this in the caption.

22. How can I build a `Shiny` app using `{ggstatsplot}` functions?

Below is an example using `ggbetweenstats` function.

```r
set.seed(123)
library(shiny)
library(rlang)

ui <- fluidPage(
  headerPanel("Example - ggbetweenstats"),
  sidebarPanel(
    selectInput("x", "xcol", "X Variable", choices = names(iris)[5]),
    selectInput("y", "ycol", "Y Variable", choices = names(iris)[1:4])
  ),
  mainPanel(plotOutput("plot"))
)

server <- function(input, output) {
  output$plot <- renderPlot({
    ggbetweenstats(iris, !!input$x, !!input$y)
  })
}

shinyApp(ui, server)
```

23. How to change size of annotations for combined plot in `grouped_*` functions?

```{r}
library(ggplot2)

grouped_ggbetweenstats(
  data = dplyr::filter(ggplot2::mpg, drv != "4"),
  x = year,
  y = hwy,
  grouping.var = drv,
  results.subtitle = FALSE,
  ## arguments given to `{patchwork}` for combining plots
  annotation.args = list(
    title = "this is my title",
    subtitle = "this is my subtitle",
    theme = ggplot2::theme(
      plot.subtitle = element_text(size = 20),
      plot.title = element_text(size = 30)
    )
  )
)
```

24. How to change size of text in the subtitle?

```{r}
ggbetweenstats(
  data = iris,
  x = Species,
  y = Sepal.Length,
  ggplot.component = list(theme(plot.subtitle = element_text(size = 20, face = "bold")))
)
```

25. How to display pairwise comparison letter in a plot?

This is not possible out of the box, but see [this](https://github.com/IndrajeetPatil/ggstatsplot/issues/654#issuecomment-948862514) comment.

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "ggscatterstats"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{ggscatterstats}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>"
)

library(ggstatsplot)
```

Lifecycle: [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

The function `ggscatterstats` is meant to provide a **publication-ready
scatterplot** with all statistical details included in the plot itself to show
association between two continuous variables. This function is also helpful
during the **data exploration** phase. We will see examples of how to use this
function in this vignette with the `ggplot2movies` dataset.

To begin with, here are some instances where you would want to use
`ggscatterstats`-

  - to check linear association between two continuous variables
  - to check distribution of two continuous variables

**Note before**: The following demo uses the pipe operator (`%>%`), so in case
you are not familiar with this operator, here is a good explanation:
<http://r4ds.had.co.nz/pipes.html>

## Correlation plot with `ggscatterstats`

To illustrate how this function can be used, we will rely on the `ggplot2movies`
dataset. This dataset provides information about movies scraped from
[IMDB](https://www.imdb.com/). Specifically, we will be using cleaned version of
this dataset included in the `{ggstatsplot}` package itself.

```{r ggplot2movies1}
library(ggstatsplot)

## see the selected data (we have data from 1813 movies)
dplyr::glimpse(movies_long)
```

Now that we have a clean dataset, we can start asking some interesting
questions. For example, let's see if the average IMDB rating for a movie has any
relationship to its budget. Additionally, let's also see which movies had a high
budget but low IMDB rating by labeling those data points. 

To reduce the processing time, let's only work with 30% of the dataset.

```{r ggscatterstats1, fig.height=6, fig.width=8}
## for reproducibility
set.seed(123)

## plot
ggscatterstats(
  data = movies_long, ## dataframe from which variables are taken
  x = budget, ## predictor/independent variable
  y = rating, ## dependent variable
  xlab = "Budget (in millions of US dollars)", ## label for the x-axis
  ylab = "Rating on IMDB", ## label for the y-axis
  label.var = title, ## variable to use for labeling data points
  label.expression = rating < 5 & budget > 100, ## expression for deciding which points to label
  point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
  xfill = "#CC79A7", ## fill for marginals on the x-axis
  yfill = "#009E73", ## fill for marginals on the y-axis
  title = "Relationship between movie budget and IMDB rating",
  caption = "Source: www.imdb.com"
)
```

There is indeed a small, but significant, positive correlation between the
amount of money studio invests in a movie and the ratings given by the
audiences.

## Grouped analysis with `grouped_ggscatterstats`

What if we want to do the same analysis do the same analysis for movies with
different MPAA (Motion Picture Association of America) film ratings (NC-17, PG,
PG-13, R)? 

`{ggstatsplot}` provides a special helper function for such instances:
`grouped_ggstatsplot`. This is merely a wrapper function around
`combine_plots`. It applies `{ggstatsplot}` across all **levels** of
a specified **grouping variable** and then combines list of individual plots
into a single plot. Note that the grouping variable can be anything: conditions
in a given study, groups in a study sample, different studies, etc. 

Let's see how we can use this function to apply `ggscatterstats` for all MPAA
ratings. Also, let's run a robust test this time.

```{r grouped1, fig.height=12, fig.width=7}
## for reproducibility
set.seed(123)

## plot
grouped_ggscatterstats(
  ## arguments relevant for ggscatterstats
  data = movies_long,
  x = budget,
  y = rating,
  grouping.var = mpaa,
  label.var = title,
  label.expression = rating < 5 & budget > 80,
  type = "r",
  # ggtheme = ggthemes::theme_tufte(),
  ## arguments relevant for combine_plots
  annotation.args = list(
    title = "Relationship between movie budget and IMDB rating",
    caption = "Source: www.imdb.com"
  ),
  plotgrid.args = list(nrow = 3, ncol = 1)
)
```

As seen from the plot, this analysis has revealed something interesting: The
relationship we found between budget and IMDB rating holds only for PG-13 and
R-rated movies. 

## Grouped analysis with `ggscatterstats` + `{purrr}`

Although this is a quick and dirty way to explore large amount of data with
minimal effort, it does come with an important limitation: reduced flexibility.
For example, if we wanted to add, let's say, a separate type of marginal
distribution plot for each MPAA rating or if we wanted to use different types of
correlations across different levels of MPAA ratings (NC-17 has only 6 movies,
so a robust correlation would be a good idea), this is not possible. But this
can be easily done using `{purrr}`.  

See the associated vignette here:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html>

## Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
raw data | `ggplot2::geom_point` | `point.args`
labels for raw data | `ggrepel::geom_label_repel` | `point.label.args`
smooth line | `ggplot2::geom_smooth` | `smooth.line.args`
marginal histograms | `ggside::geom_xsidehistogram`,  `ggside::geom_ysidehistogram` | `xsidehistogram.args`, `ysidehistogram.args`

## Summary of tests

**Hypothesis testing** and **Effect size estimation**

Type | Test | CI? | Function used
----------- | ------------------------- | --- | -----
Parametric | Pearson's correlation coefficient | ✅ | `correlation::correlation`
Non-parametric | Spearman's rank correlation coefficient | ✅ | `correlation::correlation`
Robust | Winsorized Pearson correlation coefficient | ✅ | `correlation::correlation`
Bayesian | Pearson's correlation coefficient | ✅ | `correlation::correlation`

## Reporting

If you wish to include statistical analysis results in a publication/report, the
ideal reporting practice will be a hybrid of two approaches:

- the `{ggstatsplot}` approach, where the plot contains both the visual and
numerical summaries about a statistical model, and

- the *standard* narrative approach, which provides interpretive context for the
reported statistics.

For example, let's see the following example:

```{r reporting}
ggscatterstats(mtcars, qsec, drat)
```

The narrative context (assuming `type = "parametric"`) can complement this plot
either as a figure caption or in the main text-

> Pearson's correlation test revealed that, across 32 cars, a measure of
acceleration (1/4 mile time; `qsec`) was positively correlated with rear axle
ratio (`drat`), but this effect was not statistically significant. The effect
size $(r = 0.09)$ was small, as per Cohen’s (1988) conventions. The Bayes Factor
for the same analysis revealed that the data were `r round(exp(1.20), 2)` times
more probable under the null hypothesis as compared to the alternative
hypothesis. This can be considered moderate evidence (Jeffreys, 1961) in favor
of the null hypothesis (of absence of any correlation between these two
variables).

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on GitHub: 
<https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "ggbetweenstats"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{ggbetweenstats}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
## show me all columns
options(tibble.width = Inf, pillar.bold = TRUE, pillar.subtle_num = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>"
)

if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(PMCMRplus)
}

library(ggstatsplot)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

Lifecycle: [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

The function `ggbetweenstats` is designed to facilitate **data exploration**,
and for making highly customizable **publication-ready plots**, with relevant
statistical details included in the plot itself if desired. We will see examples
of how to use this function in this vignette.

To begin with, here are some instances where you would want to use
`ggbetweenstats`-

 - to check if a continuous variable differs across multiple groups/conditions

 - to compare distributions visually and check for outliers

**Note**: This vignette uses the pipe operator (`%>%`), if you are not
familiar with this operator, here is a good explanation:
<http://r4ds.had.co.nz/pipes.html>

## Comparisons between groups with `ggbetweenstats`

To illustrate how this function can be used, we will use the `gapminder` dataset
throughout this vignette. This dataset provides values for life expectancy, GDP
per capita, and population, at 5 year intervals, from 1952 to 2007, for each of
142 countries (courtesy [Gapminder Foundation](https://www.gapminder.org/)).
Let's have a look at the data-

```{r gapminder}
library(gapminder)

dplyr::glimpse(x = gapminder::gapminder)
```

**Note**: For the remainder of the vignette, we're going to exclude *Oceania*
from the analysis simply because there are so few observations (countries).

Suppose the first thing we want to inspect is the distribution of life
expectancy for the countries of a continent in 2007. We also want to know if the
mean differences in life expectancy between the continents is statistically
significant.

The simplest form of the function call is-

```{r ggbetweenstats1, fig.height = 6, fig.width = 8}
## since the confidence intervals for the effect sizes are computed using
## bootstrapping, important to set a seed for reproducibility
set.seed(123)

## function call
ggbetweenstats(
  data = dplyr::filter(gapminder::gapminder, year == 2007, continent != "Oceania"),
  x = continent,
  y = lifeExp
)
```

**Note**:
  
  - The function automatically decides whether an independent samples *t*-test
    is preferred (for 2 groups) or a Oneway ANOVA (3 or more groups). based on
    the number of levels in the grouping variable.
    
  - The output of the function is a `ggplot` object which means that it can be
    further modified with `{ggplot2}` functions.

As can be seen from the plot, the function by default returns Bayes Factor for
the test. If the null hypothesis can't be rejected with the null hypothesis
significance testing (NHST) approach, the Bayesian approach can help index
evidence in favor of the null hypothesis (i.e., $BF_{01}$).

By default, natural logarithms are shown because Bayes Factor values can
sometimes be pretty large. Having values on logarithmic scale also makes it easy
to compare evidence in favor alternative ($BF_{10}$) versus null ($BF_{01}$)
hypotheses (since $log_{e}(BF_{01}) = - log_{e}(BF_{10})$). 

We can make the output much more aesthetically pleasing as well as informative
by making use of the many optional parameters in `ggbetweenstats`. We'll add a
title and caption, better `x` and `y` axis labels, and tag and label the
outliers in the data. We can and will change the overall theme as well as the
color palette in use.

```{r ggbetweenstats2, fig.height = 6, fig.width = 8}
## for reproducibility
set.seed(123)
library(ggstatsplot)
library(gapminder)

## plot
ggbetweenstats(
  data = dplyr::filter(gapminder, year == 2007, continent != "Oceania"),
  x = continent, ## grouping/independent variable
  y = lifeExp, ## dependent variables
  type = "robust", ## type of statistics
  xlab = "Continent", ## label for the x-axis
  ylab = "Life expectancy", ## label for the y-axis
  plot.type = "boxviolin", ## type of plot
  outlier.tagging = TRUE, ## whether outliers should be flagged
  outlier.coef = 1.5, ## coefficient for Tukey's rule
  outlier.label = country, ## label to attach to outlier values
  outlier.label.args = list(color = "red"), ## outlier point label color
  ## turn off messages
  ggtheme = ggplot2::theme_gray(), ## a different theme
  package = "yarrr", ## package from which color palette is to be taken
  palette = "info2", ## choosing a different color palette
  title = "Comparison of life expectancy across continents (Year: 2007)",
  caption = "Source: Gapminder Foundation"
) + ## modifying the plot further
  ggplot2::scale_y_continuous(
    limits = c(35, 85),
    breaks = seq(from = 35, to = 85, by = 5)
  )
```

As can be appreciated from the effect size (partial eta squared) of 0.635, there
are large differences in the mean life expectancy across continents.
Importantly, this plot also helps us appreciate the distributions within any
given continent. For example, although Asian countries are doing much better
than African countries, on average, Afghanistan has a particularly grim average
for the Asian continent, possibly reflecting the war and the political turmoil.

So far we have only used a classic parametric test and a boxviolin plot, 
but we can also use other available options:

  - The `type` (of test) argument also accepts the following abbreviations:
    `"p"` (for *parametric*), `"np"` (for *nonparametric*), `"r"` (for
    *robust*), `"bf"` (for *Bayes Factor*). 

  - The type of plot to be displayed can also be modified (`"box"`, `"violin"`,
  or `"boxviolin"`).

  - The color palettes can be modified.

Let's use the `combine_plots` function to make one plot from four separate
plots that demonstrates all of these options. Let's compare life expectancy for
all countries for the first and last year of available data 1957 and 2007. We
will generate the plots one by one and then use `combine_plots` to merge them
into one plot with some common labeling. It is possible, but not necessarily
recommended, to make each plot have different colors or themes.

For example,
```{r ggbetweenstats3, fig.height = 10, fig.width = 12}
## for reproducibility
set.seed(123)
library(ggstatsplot)
library(gapminder)

## selecting subset of the data
df_year <- dplyr::filter(gapminder::gapminder, year == 2007 | year == 1957)

## parametric t-test and box plot
p1 <- ggbetweenstats(
  data = df_year,
  x = year,
  y = lifeExp,
  xlab = "Year",
  ylab = "Life expectancy",
  plot.type = "box",
  type = "p",
  conf.level = 0.99,
  title = "Parametric test",
  package = "ggsci",
  palette = "nrc_npg"
)

## Mann-Whitney U test (nonparametric t) and violin plot
p2 <- ggbetweenstats(
  data = df_year,
  x = year,
  y = lifeExp,
  xlab = "Year",
  ylab = "Life expectancy",
  plot.type = "violin",
  type = "np",
  conf.level = 0.99,
  title = "Non-parametric Test (violin plot)",
  package = "ggsci",
  palette = "uniform_startrek"
)

## robust t-test and boxviolin plot
p3 <- ggbetweenstats(
  data = df_year,
  x = year,
  y = lifeExp,
  xlab = "Year",
  ylab = "Life expectancy",
  plot.type = "boxviolin",
  type = "r",
  conf.level = 0.99,
  title = "Robust Test (box & violin plot)",
  tr = 0.005,
  package = "wesanderson",
  palette = "Royal2",
  k = 3
)

## Bayes Factor for parametric t-test and boxviolin plot
p4 <- ggbetweenstats(
  data = df_year,
  x = year,
  y = lifeExp,
  xlab = "Year",
  ylab = "Life expectancy",
  type = "bayes",
  plot.type = "box",
  title = "Bayesian Test (box plot)",
  package = "ggsci",
  palette = "nrc_npg"
)

## combining the individual plots into a single plot
combine_plots(
  list(p1, p2, p3, p4),
  plotgrid.args = list(nrow = 2),
  annotation.args = list(
    title = "Comparison of life expectancy between 1957 and 2007",
    caption = "Source: Gapminder Foundation"
  )
)
```

## Grouped analysis with `grouped_ggbetweenstats`

What if we want to analyze both by continent and between 1957 and 2007? A
combination of our two previous efforts. 

`{ggstatsplot}` provides a special helper function for such instances:
`grouped_ggbetweenstats`. This is merely a wrapper function around
`combine_plots`. It applies `ggbetweenstats` across all **levels**
of a specified **grouping variable** and then combines list of individual plots
into a single plot. Note that the grouping variable can be anything: conditions
in a given study, groups in a study sample, different studies, etc.

Let's focus on the same 4 continents for the following years: 1967, 1987, 2007.
Also, let's carry out pairwise comparisons to see if there differences between
every pair of continents.

```{r grouped1, fig.height = 18, fig.width = 8}
## for reproducibility
set.seed(123)

## select part of the dataset and use it for plotting
gapminder::gapminder %>%
  dplyr::filter(year %in% c(1967, 1987, 2007), continent != "Oceania") %>%
  grouped_ggbetweenstats(
    ## arguments relevant for ggbetweenstats
    x = continent,
    y = lifeExp,
    grouping.var = year,
    xlab = "Continent",
    ylab = "Life expectancy",
    pairwise.display = "significant", ## display only significant pairwise comparisons
    p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
    #ggtheme = ggthemes::theme_tufte(),
    package = "ggsci",
    palette = "default_jco",
    outlier.tagging = TRUE,
    outlier.label = country,
    ## arguments relevant for combine_plots
    annotation.args = list(title = "Changes in life expectancy across continents (1967-2007)"),
    plotgrid.args = list(nrow = 3)
  )
```

As seen from the plot, although the life expectancy has been improving steadily
across all continents as we go from 1967 to 2007, this improvement has not been
happening at the same rate for all continents. Additionally, irrespective of
which year we look at, we still find significant differences in life expectancy
across continents which have been surprisingly consistent across five decades
(based on the observed effect sizes).

## Grouped analysis with `ggbetweenstats` + `{purrr}` 

Although this grouping function provides a quick way to explore the data, it
leaves much to be desired. For example, the same type of plot and test is
applied for all years, but maybe we want to change this for different years, or
maybe we want to gave different effect sizes for different years. This type of
customization for different levels of a grouping variable is not possible with
`grouped_ggbetweenstats`, but this can be easily achieved using the `{purrr}`
package. 

See the associated vignette here:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html>

## Within-subjects designs

For repeated measures designs, `ggwithinstats` function can be used:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggwithinstats.html>

## Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
raw data | `ggplot2::geom_point` | `point.args`
box plot | `ggplot2::geom_boxplot` | ❌
density plot | `ggplot2::geom_violin` | `violin.args`
centrality measure point | `ggplot2::geom_point` | `centrality.point.args`
centrality measure label | `ggrepel::geom_label_repel` | `centrality.label.args`
outlier point | `ggplot2::stat_boxplot` | ❌
outlier label | `ggrepel::geom_label_repel` | `outlier.label.args`
pairwise comparisons | `ggsignif::geom_signif` | `ggsignif.args`

## Summary of tests

**Central tendency measure**

Type | Measure | Function used
----------- | --------- | ------------------ 
Parametric | mean | `parameters::describe_distribution`
Non-parametric | median | `parameters::describe_distribution`
Robust | trimmed mean | `parameters::describe_distribution`
Bayesian | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution`

**Hypothesis testing**

Type | No. of groups | Test | Function used
----------- | --- | ------------------------- | -----
Parametric | > 2 | Fisher's or Welch's one-way ANOVA | `stats::oneway.test`
Non-parametric | > 2 | Kruskal–Wallis one-way ANOVA | `stats::kruskal.test`
Robust | > 2 | Heteroscedastic one-way ANOVA for trimmed means | `WRS2::t1way`
Bayes Factor | > 2 | Fisher's ANOVA | `BayesFactor::anovaBF`
Parametric | 2 | Student's or Welch's *t*-test | `stats::t.test`
Non-parametric | 2 | Mann–Whitney *U* test | `stats::wilcox.test`
Robust | 2 |  Yuen's test for trimmed means | `WRS2::yuen`
Bayesian | 2 | Student's *t*-test | `BayesFactor::ttestBF`

**Effect size estimation**

Type | No. of groups | Effect size | CI? | Function used
----------- | --- | ------------------------- | --- | -----
Parametric | > 2 | $\eta_{p}^2$, $\omega_{p}^2$ | ✅ | `effectsize::omega_squared`, `effectsize::eta_squared`
Non-parametric | > 2 | $\epsilon_{ordinal}^2$ | ✅ | `effectsize::rank_epsilon_squared`
Robust | > 2 | $\xi$ (Explanatory measure of effect size) | ✅ | `WRS2::t1way`
Bayes Factor | > 2 | $R_{posterior}^2$ | ✅ | `performance::r2_bayes`
Parametric | 2 | Cohen's *d*, Hedge's *g* | ✅ | `effectsize::cohens_d`, `effectsize::hedges_g`
Non-parametric | 2 | *r* (rank-biserial correlation) | ✅ | `effectsize::rank_biserial`
Robust | 2 |  $\xi$ (Explanatory measure of effect size) | ✅ | `WRS2::yuen.effect.ci`
Bayesian | 2 | $\delta_{posterior}$ | ✅ | `bayestestR::describe_posterior`

**Pairwise comparison tests**

Type | Equal variance? | Test | *p*-value adjustment? | Function used
----------- | --- | ------------------------- | --- | -----
Parametric | No | Games-Howell test | ✅ | `PMCMRplus::gamesHowellTest`
Parametric | Yes | Student's *t*-test | ✅ | `stats::pairwise.t.test`
Non-parametric | No | Dunn test | ✅ | `PMCMRplus::kwAllPairsDunnTest`
Robust | No | Yuen's trimmed means test | ✅ | `WRS2::lincon`
Bayes Factor | ❌ | Student's *t*-test | ❌ | `BayesFactor::ttestBF`

## Reporting

If you wish to include statistical analysis results in a publication/report, the
ideal reporting practice will be a hybrid of two approaches:

- the `{ggstatsplot}` approach, where the plot contains both the visual and
numerical summaries about a statistical model, and

- the *standard* narrative approach, which provides interpretive context for the
reported statistics.

For example, let's see the following example:



```{r reporting}
ggbetweenstats(ToothGrowth, supp, len)
```

The narrative context (assuming `type = "parametric"`) can complement this plot
either as a figure caption or in the main text-

> Welch's *t*-test revealed that, across 60 guinea pigs, although the tooth
length was higher when the animal received vitamin C via orange juice as
compared to via ascorbic acid, this effect was not statistically significant.
The effect size $(g = 0.49)$ was medium, as per Cohen’s (1988) conventions. The
Bayes Factor for the same analysis revealed that the data were `r round(exp(0.18), 2)` times more probable under the alternative hypothesis as
compared to the null hypothesis. This can be considered weak evidence
(Jeffreys, 1961) in favor of the alternative hypothesis.

Similar reporting style can be followed when the function performs one-way ANOVA
instead of a *t*-test.

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "Graphic design and statistical reporting principles"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{Graphic design and statistical reporting principles}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: paper.bib
csl: apa.csl  
---

```{r setup, include = FALSE}
## show me all columns
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  dpi = 300,
  out.width = "100%",
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)

library(ggstatsplot)
library(ggplot2)

if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(PMCMRplus)
}
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

## Graphic design principles

### Graphical perception
   
Graphical perception involves visual decoding of the encoded information in
graphs. `{ggstatsplot}` incorporates the paradigm proposed in ([@cleveland1985],
Chapter 4) to facilitate making visual judgments about quantitative information
effortless and almost instantaneous. Based on experiments, Cleveland proposes
that there are ten elementary graphical-perception tasks that we perform to
visually decode quantitative information in graphs (organized from most to least
accurate; [@cleveland1985], p.254)-

  * Position along a common scale

  * Position along identical, non-aligned scales

  * Length

  * Angle (Slope)

  * Area

  * Volume

  * Color hue

So the key principle of Cleveland's paradigm for data display is-

> "We should encode data on a graph so that the visual decoding involves
[graphical-perception] tasks as high in the ordering as possible."

For example, decoding the data point values in `ggbetweenstats` requires
position judgments along a common scale:

```{r fig1, fig.height = 9, fig.width = 10, fig.cap = "Note that assessing differences in mean values between groups has been made easier with the help of \\textit{position} of data points along a common scale (the Y-axis) and labels."}
## for reproducibility
set.seed(123)

## plot
ggbetweenstats(
  data = dplyr::filter(
    movies_long,
    genre %in% c("Action", "Action Comedy", "Action Drama", "Comedy")
  ),
  x = genre,
  y = rating,
  title = "IMDB rating by film genre",
  xlab = "Genre",
  ylab = "IMDB rating (average)",
  outlier.tagging = TRUE,
  outlier.label = title
)
```

There are few instances where `{ggstatsplot}` diverges from recommendations made
in Cleveland's paradigm:

  - For the categorical/nominal data, `{ggstatsplot}` uses pie charts which rely
    on *angle* judgments, which are less accurate (as compared to bar graphs,
    e.g., which require *position* judgments). This shortcoming is assuaged to
    some degree by using plenty of labels that describe percentages for all
    slices. This makes angle judgment unnecessary and pre-vacates any concerns
    about inaccurate judgments about percentages. Additionally, it also provides
    alternative function to `ggpiestats` for working with categorical variables:
    `ggbarstats`.

```{r fig2, fig.height = 4, fig.width = 10, fig.cap = "Pie charts don't follow Cleveland's paradigm to data display because they rely on less accurate angle judgments. `{ggstatsplot}` sidesteps this issue by always labelling percentages for pie slices, which makes angle judgments unnecessary."}
## for reproducibility
set.seed(123)

## plot
ggpiestats(
  data = movies_long,
  x = genre,
  y = mpaa,
  title = "Distribution of MPAA ratings by film genre",
  legend.title = "layout"
)
```

  - Cleveland's paradigm also emphasizes that *superposition* of data is better
    than *juxtaposition* ([@cleveland1985], p.201) because this allows for a
    more incisive comparison of the values from different parts of the dataset.
    This recommendation is violated in all `grouped_` variants of the function.
    Note that the range for Y-axes are no longer the same across juxtaposed
    subplots and so visually comparing the data becomes difficult. On the other
    hand, in the superposed plot, all data have the same range and coloring
    different parts makes the visual discrimination of different components of
    the data, and their comparison, easier. But the goal of `grouped_` variants
    of functions is to not only show different aspects of the data but also to
    run statistical tests and showing detailed results for all aspects of the
    data in a superposed plot is difficult. Therefore, this is a compromise
    `{ggstatsplot}` is comfortable with, at least to produce plots for quick
    exploration of different aspects of the data.

```{r fig3, fig.height = 12, fig.width = 10, fig.cap = "Comparing different aspects of data is much more accurate in (\\textit{a}) a \\textit{superposed} plot, which is recommended in Cleveland's paradigm, than in (\\textit{b}) a \\textit{juxtaposed} plot, which is how it is implemented in `{ggstatsplot}` package. This is because displaying detailed results from statistical tests would be difficult in a superposed plot."}
## for reproducibility
set.seed(123)
library(ggplot2)
library(ggstatsplot)

## creating a smaller dataframe
df <- dplyr::filter(movies_long, genre %in% c("Comedy", "Drama"))

## plot
combine_plots(
  plotlist = list(
    ## plot 1: superposition
    ggplot(data = df, mapping = aes(x = length, y = rating, color = genre)) +
      geom_jitter(size = 3, alpha = 0.5) +
      geom_smooth(method = "lm") +
      labs(title = "superposition (recommended in Cleveland's paradigm)") +
      theme_ggstatsplot(),
    ## plot 2: juxtaposition
    grouped_ggscatterstats(
      data = df,
      x = length,
      y = rating,
      grouping.var = genre,
      marginal = FALSE,
      annotation.args = list(title = "juxtaposition (`{ggstatsplot}` implementation in `grouped_` functions)")
    )
  ),
  ## combine for comparison
  annotation.args = list(title = "Two ways to compare different aspects of data"),
  plotgrid.args = list(nrow = 2)
)
```

The `grouped_` plots follow the *Shrink Principle* ([@tufte2001], p.166-7) for
high-information graphics, which dictates that the data density and the size of
the data matrix can be maximized to exploit maximum resolution of the available
data-display technology. Given the large maximum resolution afforded by most
computer monitors today, saving `grouped_` plots with appropriate resolution
ensures no loss in legibility with reduced graphics area.

### Graphical excellence
   
Graphical excellence consists of communicating complex ideas with clarity and in
a way that the viewer understands the greatest number of ideas in a short amount
of time all the while not quoting the data out of context. The package follows
the principles for *graphical integrity* [@tufte2001]:

  - The physical representation of numbers is proportional to the numerical
    quantities they represent. The plot show how means (in `ggbetweenstats`) or
    percentages (`ggpiestats`) are proportional to the vertical distance or the
    area, respectively).

  - All important events in the data have clear, detailed, and thorough labeling
    plot shows how `ggbetweenstats` labels means, sample size information,
    outliers, and pairwise comparisons; same can be appreciated for `ggpiestats`
    and `gghistostats` plots. Note that data labels in the data region are
    designed in a way that they don't interfere with our ability to assess the
    overall pattern of the data ([@cleveland1985];

p.44-45). This is achieved by using `ggrepel` package to place labels in a way
that reduces their visual prominence.

  - None of the plots have *design* variation (e.g., abrupt change in scales)
    over the surface of a same graphic because this can lead to a false
    impression about variation in *data*.

  - The number of information-carrying dimensions never exceed the number of
    dimensions in the data (e.g., using area to show one-dimensional data).

  - All plots are designed to have no **chartjunk** (like moiré vibrations, fake
    perspective, dark grid lines, etc.) ([@tufte2001], Chapter 5).

There are some instances where `{ggstatsplot}` graphs don't follow principles of
clean graphics, as formulated in the Tufte theory of data graphics
([@tufte2001], Chapter 4). The theory has four key principles:

  1. Above all else show the data.

  2. Maximize the data-ink ratio.

  3. Erase non-data-ink.

  4. Erase redundant data-ink, within reason.

In particular, default plots in `{ggstatsplot}` can sometimes violate one of the
principles from 2-4. According to these principles, every bit of ink should have
reason for its inclusion in the graphic and should convey some new information
to the viewer. If not, such ink should be removed. One instance of this is
bilateral symmetry of data measures. For example, in the figure below, we can
see that both the box and violin plots are mirrored, which consumes twice the
space in the graphic without adding any new information. But this redundancy is
tolerated for the sake of beauty that such symmetrical shapes can bring to the
graphic. Even Tufte admits that efficiency is but one consideration in the
design of statistical graphics ([@tufte2001],

p. 137). Additionally, these principles were formulated in an era in which
   computer graphics had yet to revolutionize the ease with which graphics could
   be produced and thus some of the concerns about minimizing data-ink for
   easier production of graphics are not as relevant as they were.

### Statistical variation

One of the important functions of a plot is to show the variation in the data,
which comes in two forms:

  - **Measurement noise**: In `{ggstatsplot}`, the actual variation in
    measurements is shown by plotting a combination of (jittered) raw data
    points with a boxplot laid on top or a histogram. None of the plots, where
    empirical distribution of the data is concerned, show the sample standard
    deviation because they are poor at conveying information about limits of the
    sample and presence of outliers ([@cleveland1985], p.220).

```{r fig5, fig.height = 6, fig.width = 8, fig.cap = "Distribution of a variable shown using `gghistostats`."}
## for reproducibility
set.seed(123)

## plot
gghistostats(
  data = morley,
  x = Speed,
  test.value = 792,
  xlab = "Speed of light (km/sec, with 299000 subtracted)",
  title = "Distribution of measured Speed of light",
  caption = "Note: Data collected across 5 experiments (20 measurements each)"
)
```

  - **Sample-to-sample statistic variation**: Although, traditionally, this
    variation has been shown using the standard error of the mean (SEM) of the
    statistic, `{ggstatsplot}` plots instead use 95% confidence intervals. This is
    because the interval formed by error bars correspond to a 68% confidence
    interval, which is not a particularly interesting interval
    ([@cleveland1985], p.222-225).

```{r fig6, fig.height = 5, fig.width = 5, fig.cap = "Sample-to-sample variation in regression estimates is displayed using confidence intervals in `ggcoefstats`."}
## for reproducibility
set.seed(123)

## creating model object
mod <- lme4::lmer(
  formula = total.fruits ~ nutrient + rack + (nutrient | gen),
  data = lme4::Arabidopsis
)

## plot
ggcoefstats(x = mod)
```

## Statistical analysis

### Data requirements

As an extension of `{ggplot2}`, `{ggstatsplot}` has the same expectations about the
structure of the data. More specifically,

  - The data should be organized following the principles of *tidy data*, which
    specify how statistical structure of a data frame (variables and
    observations) should be mapped to physical structure (columns and rows).
    More specifically, tidy data means all variables have their own columns and
    each row corresponds to a unique observation ([@wickhamTidyData2014]).

  - All `{ggstatsplot}` functions remove `NA`s from variables of interest (similar
    to `{ggplot2}`; [@wickham2016], p.207) in the data and display total sample
    size (*n*, either observations for between-subjects or pairs for
    within-subjects designs) in the subtitle to inform the user/reader about the
    number of observations included for both the statistical analysis and the
    visualization. But, when sample sizes differ *across* tests in the same
    function, `{ggstatsplot}` makes an effort to inform the user of this aspect.
    For example, `ggcorrmat` features several correlation test pairs and,
    depending on variables in a given pair, the sample sizes may vary.

```{r fig4, fig.height = 5, fig.width = 10, fig.cap = "`{ggstatsplot}` functions remove `NA`s from variables of interest and display total sample size \\textit{n}, but they can give more nuanced information about sample sizes when \\textit{n} differs across tests. For example, `ggcorrmat` will display (\\textit{a}) only one total sample size once when no `NA`s present, but (\\textit{b}) will instead show minimum, median, and maximum sample sizes across all correlation tests when `NA`s are present across correlation variables."}
## for reproducibility
set.seed(123)

## creating a new dataset without any NAs in variables of interest
msleep_no_na <-
  dplyr::filter(
    ggplot2::msleep,
    !is.na(sleep_rem), !is.na(awake), !is.na(brainwt), !is.na(bodywt)
  )

## variable names vector
var_names <- c("REM sleep", "time awake", "brain weight", "body weight")

## combining two plots using helper function in `{ggstatsplot}`
combine_plots(
  plotlist = purrr::pmap(
    .l = list(data = list(msleep_no_na, ggplot2::msleep)),
    .f = ggcorrmat,
    cor.vars = c(sleep_rem, awake:bodywt),
    cor.vars.names = var_names,
    colors = c("#B2182B", "white", "#4D4D4D"),
    title = "Correlalogram for mammals sleep dataset",
    subtitle = "sleep units: hours; weight units: kilograms"
  ),
  plotgrid.args = list(nrow = 1)
)
```

### Statistical reporting

But why would combining statistical analysis with data visualization be helpful?
We list few reasons below-

  - A recent survey [@nuijten2016] revealed that one in eight papers in major
    psychology journals contained a grossly inconsistent *p*-value that may have
    affected the statistical conclusion. `{ggstatsplot}` helps avoid such
    reporting errors: Since the plot and the statistical analysis are yoked
    together, the chances of making an error in reporting the results are
    minimized. One need not write the results manually or copy-paste them from a
    different statistics software program (like SPSS, SAS, and so on).

The default setting in `{ggstatsplot}` is to produce plots with statistical
details included. Most often than not, these results are displayed as a
`subtitle` in the plot. Great care has been taken into which details are
included in statistical reporting and why.

![Template for reporting statistical details](stats_reporting_format.png)

APA guidelines [@apa2009] are followed by default while reporting statistical
details:

  - Percentages are displayed with no decimal places.

  - Correlations, *t*-tests, and $\chi^2$-tests are reported with the degrees of
    freedom in parentheses and the significance level.

  - ANOVAs are reported with two degrees of freedom and the significance level.

  - Regression results are presented with the unstandardized or standardized
    estimate (beta), whichever was specified by the user, along with the
    statistic (depending on the model, this can be a *t*, *F*, or *z* statistic)
    and the corresponding significance level.

  - With the exception of *p*-values, most statistics are rounded to two decimal
    places by default.

### Dealing with **null results**: 

All functions therefore by default return Bayesian in favor of the null
hypothesis by default. If the null hypothesis can't be rejected with the null
hypothesis significance testing (NHST) approach, the Bayesian approach can help
index evidence in favor of the null hypothesis (i.e., $BF_{01}$). By default,
natural logarithms are shown because Bayesian values can sometimes be pretty
large. Having values on logarithmic scale also makes it easy to compare evidence
in favor alternative ($BF_{10}$) versus null ($BF_{01}$) hypotheses (since
$log_{e}(BF_{01}) = - log_{e}(BF_{01})$).

###  Avoiding the **"p-value error"**:  

The *p*-value indexes the probability that the researchers have falsely rejected
a true null hypothesis (Type I error, i.e.) and can rarely be *exactly* 0. And
yet over 97,000 manuscripts on Google Scholar report the *p*-value to be `p
=0.000`, putatively due to relying on default computer software outputs
[@lilienfeld2015]. All *p*-values displayed in `{ggstatsplot}` plots avoid this
mistake. Anything less than `p < 0.001` is displayed as such. The package deems
it unimportant how infinitesimally small the *p*-values are and, instead, puts
emphasis on the effect size magnitudes and their 95% CIs.

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>

---
title: "gghistostats"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{gghistostats}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
## show me all columns
options(tibble.width = Inf, pillar.bold = TRUE, pillar.subtle_num = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>"
)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

Lifecycle: [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

The function `gghistostats` can be used for **data exploration**
and to provide an easy way to make **publication-ready histograms** with
appropriate and selected statistical details embedded in the plot itself. In
this vignette we will explore several examples of how to use it.

Some instances where you would want to use `gghistostats`-

  - to inspect distribution of a continuous variable
  - to test if the mean of a sample variable is different from a specified value
  (population parameter)

## Statistical analysis with `gghistostats`

Let's begin with a very simple example from the `psych` package
(`psych::sat.act`), a sample of 700 self-reported scores on the SAT Verbal, SAT
Quantitative and ACT tests. ACT composite scores may range from 1 - 36. National
norms have a mean of 20.

```{r psychact}
## loading needed libraries
library(ggstatsplot)
library(psych)
library(dplyr)

## looking at the structure of the data using glimpse
dplyr::glimpse(psych::sat.act)
```

To get a simple histogram with no statistics and no special information.
`gghistostats` will by default choose a binwidth `max(x) - min(x) / sqrt(N)`.
You should always check this value and explore multiple widths to find the best
to illustrate the stories in your data since histograms are sensitive to
binwidth.

Let's display the national norms (labeled as "Test") and test the hypothesis
that our sample mean is the same as our national population mean of 20 using a
parametric one sample *t*-test (`type = "p"`).

```{r psychact3, fig.height = 6, fig.width = 7}
set.seed(123)

gghistostats(
  data = psych::sat.act, ## data from which variable is to be taken
  x = ACT, ## numeric variable
  xlab = "ACT Score", ## x-axis label
  title = "Distribution of ACT Scores", ## title for the plot
  test.value = 20, ## test value
  caption = "Data courtesy of: SAPA project (https://sapa-project.org)"
)
```

`gghistostats` computed Bayes Factors to quantify the likelihood of the **research** (BF10) and
the **null** hypothesis (BF01). In our current example, the Bayes Factor value provides
**very strong evidence** [(Kass and Rafferty, 1995)](https://www.stat.washington.edu/raftery/Research/PDF/kass1995.pdf)
in favor of the research hypothesis: these ACT scores are much higher than the
national average.  The log(Bayes factor) of 492.5 means the odds are 7.54e+213:1
that this sample is different.

## Grouped analysis with `grouped_gghistostats`

What if we want to do the same analysis separately for each gender? 
`{ggstatsplot}` provides a special helper function for such instances:
`grouped_gghistostats`. This is merely a wrapper function around
`combine_plots`. It applies `gghistostats` across all **levels** of
a specified **grouping variable** and then combines the individual plots into a
single plot. Note that the grouping variable can be anything: conditions in a
given study, groups in a study sample, different studies, etc.

Let's see how we can use this function to apply `gghistostats` to accomplish our
task. 

```{r grouped1, fig.height = 10, fig.width = 7}
set.seed(123)

grouped_gghistostats(
  ## arguments relevant for gghistostats
  data = psych::sat.act,
  x = ACT, ## same outcome variable
  xlab = "ACT Score",
  grouping.var = gender, ## grouping variable males = 1, females = 2
  type = "robust", ## robust test: one-sample percentile bootstrap
  test.value = 20, ## test value against which sample mean is to be compared
  centrality.line.args = list(color = "#D55E00", linetype = "dashed"),
  # ggtheme = ggthemes::theme_stata(), ## changing default theme
   ## turn off ggstatsplot theme layer
  ## arguments relevant for combine_plots
  annotation.args = list(
    title = "Distribution of ACT scores across genders",
    caption = "Data courtesy of: SAPA project (https://sapa-project.org)"
  ),
  plotgrid.args = list(nrow = 2)
)
```

As can be seen from these plots, the mean value is much higher than the national
norm. Additionally, we see the benefits of plotting this data separately for
each gender. We can *see* the differences in distributions.

## Grouped analysis with `{purrr}`

Although this is a quick and dirty way to explore a large amount of data with
minimal effort, it does come with an important limitation: reduced flexibility.
For example, if we wanted to add, let's say, a separate `test.value` argument
for each gender, this is not possible with `grouped_gghistostats`. 

For cases like these, or to run separate  kinds of tests (robust for some,
parametric for other, while Bayesian for some other levels of the group) it
would be better to use `{purrr}`.

See the associated vignette here:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html>

## Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
histogram bin | `ggplot2::stat_bin` | `bin.args`
centrality measure line | `ggplot2::geom_vline` | `centrality.line.args`
normality curve | `ggplot2::stat_function` | `normal.curve.args`

## Summary of tests

**Central tendency measure**

Type | Measure | Function used
----------- | --------- | ------------------ 
Parametric | mean | `parameters::describe_distribution`
Non-parametric | median | `parameters::describe_distribution`
Robust | trimmed mean | `parameters::describe_distribution`
Bayesian | MAP (maximum *a posteriori* probability) estimate | `parameters::describe_distribution`

**Hypothesis testing**

Type | Test | Function used
------------------ | ------------------------- | -----
Parametric | One-sample Student's *t*-test | `stats::t.test`
Non-parametric | One-sample Wilcoxon test | `stats::wilcox.test`
Robust | Bootstrap-*t* method for one-sample test | `WRS2::trimcibt`
Bayesian | One-sample Student's *t*-test | `BayesFactor::ttestBF`

**Effect size estimation**

Type | Effect size | CI? | Function used
------------ | ----------------------- | --- | -----
Parametric | Cohen's *d*, Hedge's *g* | ✅ | `effectsize::cohens_d`, `effectsize::hedges_g`
Non-parametric | *r* (rank-biserial correlation) | ✅ | `effectsize::rank_biserial`
Robust | trimmed mean | ✅ | `WRS2::trimcibt`
Bayes Factor | $\delta_{posterior}$ | ✅ | `bayestestR::describe_posterior`

## Reporting

If you wish to include statistical analysis results in a publication/report, the
ideal reporting practice will be a hybrid of two approaches:

- the `{ggstatsplot}` approach, where the plot contains both the visual and
numerical summaries about a statistical model, and

- the *standard* narrative approach, which provides interpretive context for the
reported statistics.

For example, let's see the following example:

```{r reporting}
gghistostats(trees, Height, test.value = 75)
```

The narrative context (assuming `type = "parametric"`) can complement this plot
either as a figure caption or in the main text-

> Student's *t*-test revealed that, across 31 felled black cherry trees,
although the height was higher than expected height of 75 ft., this effect was
not statistically significant. The effect size $(g = 0.15)$ was small, as per
Cohen’s (1988) conventions. The Bayes Factor for the same analysis revealed that
the data were `r round(exp(1.3), 2)` times more probable under the null
hypothesis as compared to the alternative hypothesis. This can be considered
moderate evidence (Jeffreys, 1961) in favor of the null hypothesis.

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on GitHub: 
<https://github.com/IndrajeetPatil/ggstatsplot/issues>

---
title: "Pairwise comparisons with `{ggstatsplot}`"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    eval: FALSE
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{Pairwise comparisons with `{ggstatsplot}`}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
## pretty tibble printing
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>",
  package.startup.message = FALSE
)

if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(PMCMRplus)
}

library(ggstatsplot)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

## Introduction

Pairwise comparisons with `{ggstatsplot}`.

## Summary of types of statistical analyses

Following table contains a brief summary of the currently supported pairwise
comparison tests-

### Between-subjects design

Type | Equal variance? | Test | *p*-value adjustment? | Function used
----------- | --- | ------------------------- | --- | -----
Parametric | No | Games-Howell test | ✅ | `PMCMRplus::gamesHowellTest`
Parametric | Yes | Student's *t*-test | ✅ | `stats::pairwise.t.test`
Non-parametric | No | Dunn test | ✅ | `PMCMRplus::kwAllPairsDunnTest`
Robust | No | Yuen's trimmed means test | ✅ | `WRS2::lincon`
Bayesian | `NA` | Student's *t*-test | `NA` | `BayesFactor::ttestBF`

### Within-subjects design

Type | Test | *p*-value adjustment? | Function used
----------- | ---------------------------- | --- | -----
Parametric | Student's *t*-test | ✅ | `stats::pairwise.t.test`
Non-parametric | Durbin-Conover test | ✅ | `PMCMRplus::durbinAllPairsTest` 
Robust | Yuen's trimmed means test | ✅ | `WRS2::rmmcp`
Bayesian | Student's *t*-test | `NA` | `BayesFactor::ttestBF`

## Examples

Here we will see specific examples of how to use this function for different
types of

  - designs (between or within subjects)
  - statistics (parametric, non-parametric, robust, Bayesian)
  - *p*-value adjustment methods

### Between-subjects design

```{r}
## for reproducibility
set.seed(123)
library(ggstatsplot)
library(statsExpressions) ## for data

## parametric
## if `var.equal = TRUE`, then Student's *t*-test will be run
pairwise_comparisons(
  data            = ggplot2::msleep,
  x               = vore,
  y               = brainwt,
  type            = "parametric",
  var.equal       = TRUE,
  paired          = FALSE,
  p.adjust.method = "bonferroni"
)

## if `var.equal = FALSE`, then Games-Howell test will be run
pairwise_comparisons(
  data            = ggplot2::msleep,
  x               = vore,
  y               = brainwt,
  type            = "parametric",
  var.equal       = FALSE,
  paired          = FALSE,
  p.adjust.method = "bonferroni"
)

## non-parametric
pairwise_comparisons(
  data            = ggplot2::msleep,
  x               = vore,
  y               = brainwt,
  type            = "nonparametric",
  paired          = FALSE,
  p.adjust.method = "none"
)

## robust
pairwise_comparisons(
  data            = ggplot2::msleep,
  x               = vore,
  y               = brainwt,
  type            = "robust",
  paired          = FALSE,
  p.adjust.method = "fdr"
)

## Bayesian
pairwise_comparisons(
  data   = ggplot2::msleep,
  x      = vore,
  y      = brainwt,
  type   = "bayes",
  paired = FALSE
)
```


### Within-subjects design

```{r}
## for reproducibility
set.seed(123)

## parametric
pairwise_comparisons(
  data            = bugs_long,
  x               = condition,
  y               = desire,
  subject.id      = subject,
  type            = "parametric",
  paired          = TRUE,
  p.adjust.method = "BH"
)

## non-parametric
pairwise_comparisons(
  data            = bugs_long,
  x               = condition,
  y               = desire,
  subject.id      = subject,
  type            = "nonparametric",
  paired          = TRUE,
  p.adjust.method = "BY"
)

## robust
pairwise_comparisons(
  data            = bugs_long,
  x               = condition,
  y               = desire,
  subject.id      = subject,
  type            = "robust",
  paired          = TRUE,
  p.adjust.method = "hommel"
)

## Bayesian
pairwise_comparisons(
  data       = bugs_long,
  x          = condition,
  y          = desire,
  subject.id = subject,
  type       = "bayes",
  paired     = TRUE,
  bf.prior   = 0.77
)
```

## Using `pairwise_comparisons()` with `ggsignif`

### Example-1: between-subjects

```{r ggsignif, fig.height=5}
## needed libraries
set.seed(123)
library(ggplot2)
library(ggsignif)

## converting to factor
mtcars$cyl <- as.factor(mtcars$cyl)

## creating a basic plot
p <- ggplot(mtcars, aes(cyl, wt)) +
  geom_boxplot()

## using `pairwise_comparisons()` package to create a dataframe with results
set.seed(123)
(df <-
  pairwise_comparisons(mtcars, cyl, wt) %>%
  dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
  dplyr::arrange(group1))

## using `geom_signif` to display results
## (note that you can choose not to display all comparisons)
p +
  ggsignif::geom_signif(
    comparisons = list(df$groups[[1]]),
    annotations = df$label[[1]],
    test        = NULL,
    na.rm       = TRUE,
    parse       = TRUE
  )
```

### Example-2: within-subjects

```{r ggsignif2}
## needed libraries
library(ggplot2)
library(ggstatsplot)
library(ggsignif)

## creating a basic plot
p <- ggplot(WRS2::WineTasting, aes(Wine, Taste)) +
  geom_boxplot()

## using `pairwise_comparisons()` package to create a dataframe with results
set.seed(123)
(df <-
  pairwise_comparisons(
    WRS2::WineTasting,
    Wine,
    Taste,
    subject.id = Taster,
    type = "bayes",
    paired = TRUE
  ) %>%
  dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
  dplyr::arrange(group1))

## using `geom_signif` to display results
p +
  ggsignif::geom_signif(
    comparisons      = df$groups,
    map_signif_level = TRUE,
    tip_length       = 0.01,
    y_position       = c(6.5, 6.65, 6.8),
    annotations      = df$label,
    test             = NULL,
    na.rm            = TRUE,
    parse            = TRUE
  )
```

---
title: "ggpiestats"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{ggpiestats}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
options(tibble.width = Inf, pillar.bold = TRUE, pillar.subtle_num = TRUE)

knitr::opts_chunk$set(
  dpi = 300,
  out.width = "100%",
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)

library(ggstatsplot)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

Lifecycle: [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

## Introduction to `ggpiestats`

The function `ggpiestats` can be used for quick **data
exploration** and/or to prepare **publication-ready pie charts** to summarize
the statistical relationship(s) among one or more categorical variables. We will
see examples of how to use this function in this vignette.

To begin with, here are some instances where you would want to use
`ggpiestats`-

 - to check if the proportion of observations matches our hypothesized
proportion, this is typically known as a "Goodness of Fit" test

 - to see if the frequency distribution of two categorical variables are
independent of each other using the contingency table analysis

 - to check if the proportion of observations at each level of a categorical
variable is equal

**Note:** The following demo uses the pipe operator (`%>%`), if you are not
familiar with this operator, here is a good explanation:
<http://r4ds.had.co.nz/pipes.html>.

`ggpiestats` works **only** with data organized in dataframes or tibbles. It
will not work with other data structures like base-R tables or matrices. It can
operate on dataframes that are organized with one row per observation or
dataframes that have one column containing counts. This vignette provides
examples of both (see examples below).

To help demonstrate how `ggpiestats` can be used with categorical (also known as
nominal) data, a modified version of the original `Titanic` dataset (from the
`datasets` library) has been provided in the `{ggstatsplot}` package with the name
`Titanic_full`. The Titanic Passenger Survival Dataset provides information "on
the fate of passengers on the fatal maiden voyage of the ocean liner *Titanic*,
including economic status (class), sex, age, and survival."

Let's have a look at the structure of both.

```{r titanic1}
library(datasets)
library(dplyr)
library(ggstatsplot)

## looking at the original data in tabular format
dplyr::glimpse(x = Titanic)

## looking at the dataset as a tibble or dataframe
dplyr::glimpse(x = Titanic_full)
```

## Goodness of Fit with `ggpiestats`

The simplest use case for `ggpiestats` is that we want to display information
about **one** categorical or nominal variable. As part of that display or plot,
we may also choose to execute a chi-squared goodness of fit test to see whether
the proportions (or percentages) in categories of the single variable appear to
line up with our hypothesis or model. To start simple and then expand, let's say
that we'd like to display a piechart with the percentages of passengers who did
or did not survive. Our initial hypothesis is that it was no different than
flipping a coin. People had a 50/50 chance of surviving.

```{r ggpiestats3, fig.height = 5, fig.width = 6}
## since effect size confidence intervals are computed using bootstrapping, let's
## set seed for reproducibility
set.seed(123)

## to speed up the process, let's use only half of the dataset
Titanic_full_50 <- dplyr::sample_frac(Titanic_full, size = 0.5)

## plot
ggpiestats(
  data = Titanic_full_50,
  x = Survived,
  title = "Passenger survival on the Titanic", ## title for the entire plot
  caption = "Source: Titanic survival dataset", ## caption for the entire plot
  legend.title = "Survived?"
)
```

**Note:** equal proportions per category are the default, e.g. 50/50, but you
can specify any hypothesized ratio you like with `ratio` so if our hypothesis
was that 80% died and 20% survived we would add `ratio = c(.80,.20)` when we
entered the code.

## Independence (or association) with `ggpiestats`

Let's next investigate whether the passenger's gender was independent of, or
associated with, gender. The test is whether the proportion of people who
survived was different between the sexes using `ggpiestats`.

We'll modify a number of arguments to change the appearance of this plot and
showcase the flexibility of `ggpiestats`. We will:

1.  Change the plot theme to `ggplot2::theme_grey()`

2.  Change our color palette to `category10_d3` from `ggsci` package

3.  We'll customize the subtitle by being more precise about which chi squared
    test this is `stat.title = "chi squared test of independence: "`

4.  Finally, we'll make a call to `{ggplot2}` to modify the size of our plot title
    and to make it right justified

```{r ggpiestats1, fig.height = 5, fig.width = 8}
## since effect size confidence intervals are computed using bootstrapping, let's
## set seed for reproducibility
set.seed(123)

## plot
ggpiestats(
  data = Titanic_full,
  x = Survived,
  y = Sex,
  title = "Passenger survival on the Titanic by gender", ## title for the entire plot
  caption = "Source: Titanic survival dataset", ## caption for the entire plot
  legend.title = "Survived?", ## legend title
  ggtheme = ggplot2::theme_grey(), ## changing plot theme
  palette = "category10_d3", ## choosing a different color palette
  package = "ggsci", ## package to which color palette belongs
  k = 3, ## decimal places in result
  perc.k = 1 ## decimal places in percentage labels
) + ## further modification with `{ggplot2}` commands
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      color = "black",
      size = 14,
      hjust = 0
    )
  )
```

The plot clearly shows that survival rates were very different between males and
females. The Pearson's $\chi^2$-test of independence is significant given our
large sample size. Additionally, for both females and males, the survival rates
were significantly different than 50% as indicated by a goodness of fit test for
each gender.

## Grouped analysis with `grouped_ggpiestats`

What if we want to do the same analysis of gender but also factor in the
passenger's age (Age)? We have information that classifies the passengers as
Child or Adult, perhaps that makes a difference to their survival rate? 

`{ggstatsplot}` provides a special helper function for such instances:
`grouped_ggpiestats`. It is a convenient wrapper function around
`combine_plots`. It applies `ggpiestats` across all **levels** of a
specified **grouping variable** and then combines the list of individual plots
into a single plot. Note that the grouping variable can be anything: conditions
in a given study, groups in a study sample, different studies, etc.

```{r ggpiestats4, fig.height = 10, fig.width = 8}
## since effect size confidence intervals are computed using bootstrapping, let's
## set seed for reproducibility
set.seed(123)

## plot
grouped_ggpiestats(
  ## arguments relevant for gghistostats
  data = Titanic_full,
  x = Survived,
  y = Sex,
  grouping.var = Age,
  perc.k = 1,
  package = "ggsci",
  palette = "category10_d3",
  ## arguments relevant for combine_plots
  title.text = "Passenger survival on the Titanic by gender and age",
  caption.text = "Asterisks denote results from proportion tests; \n***: p < 0.001, ns: non-significant",
  plotgrid.args = list(nrow = 2)
)
```

The resulting pie charts and statistics make the story clear. For adults gender
very much matters. Women survived at much higher rates than men. For children
gender is not significantly associated with survival and both male and female
children have a survival rate that is not significantly different from 50/50.

## Grouped analysis with `ggpiestats` + `{purrr}` 

Although `grouped_ggpiestats` provides a quick way to explore the data, it
leaves much to be desired. For example, we may want to add different captions,
titles, themes, or palettes for each level of the grouping variable, etc. For
cases like these, it would be better to use `{purrr}` package.  

See the associated vignette here:
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/purrr_examples.html>

## Working with data organized by `counts` 

`ggpiestats` can also work with dataframe containing counts (aka tabled data),
i.e., when each row doesn't correspond to a unique observation. For example,
consider the following notional `fishing` dataframe containing data from two
boats (`A` and `B`) about the number of different types fish they caught in the
months of `February` and `March`. In this dataframe, each row corresponds to a
unique combination of `Boat` and `Month`.

```{r ggpiestats7, message = FALSE, warning = FALSE, fig.height = 8, fig.width = 9}
## for reproducibility
set.seed(123)

## creating a dataframe
## (this is completely fictional; I don't know first thing about fishing!)
(
  fishing <- data.frame(
    Boat = c(rep("B", 4), rep("A", 4), rep("A", 4), rep("B", 4)),
    Month = c(rep("February", 2), rep("March", 2), rep("February", 2), rep("March", 2)),
    Fish = c(
      "Bass",
      "Catfish",
      "Cod",
      "Haddock",
      "Cod",
      "Haddock",
      "Bass",
      "Catfish",
      "Bass",
      "Catfish",
      "Cod",
      "Haddock",
      "Cod",
      "Haddock",
      "Bass",
      "Catfish"
    ),
    SumOfCaught = c(25, 20, 35, 40, 40, 25, 30, 42, 40, 30, 33, 26, 100, 30, 20, 20)
  ) %>% ## converting to a tibble dataframe
    tibble::as_data_frame(x = .)
)
```

When the data is organized this way, we make a slightly different call to the
`ggpiestats` function: we use the `counts` argument. If we want to investigate
the relationship of type of fish by month (a test of independence), our command
would be:

```{r ggpiestats8, fig.height = 5, fig.width = 8}
## running `ggpiestats` with counts information
ggpiestats(
  data = fishing,
  x = Fish,
  y = Month,
  counts = SumOfCaught,
  label = "both",
  package = "ggsci",
  palette = "default_jama",
  title = "Type fish caught by month",
  caption = "Source: completely made up",
  legend.title = "Type fish caught: "
)
```

The results support our hypothesis that the type of fish caught is related to
the month in which we're fishing. The $\chi^2$ independence test results at
the top of the plot. In February we catch significantly more Haddock than we
would hypothesize for an equal distribution. Whereas in March our results
indicate there's no strong evidence that the distribution isn't equal.

## Within-subjects designs

For our final example let's imagine we're conducting clinical trials for some
new imaginary wonder drug. We have 134 subjects entering the trial. Some of them
enter healthy (*n* = 96), some of them enter the trial already being sick (*n* =
38). All of them receive our treatment or intervention. Then we check back in a
month to see if they are healthy or sick. A classic pre/post experimental
design. We're interested in seeing the change in both groupings. In the case of
within-subjects designs, you can set `paired = TRUE`, which will display results
from **McNemar test** in the subtitle.

(**Note:** If you forget to set `paired = TRUE`, the results you get will be
inaccurate.)

```{r ggpiestats9, fig.height = 5, fig.width = 8}
## seed for reproducibility
set.seed(123)

## create our imaginary data
clinical_trial <-
  tibble::tribble(
    ~SickBefore, ~SickAfter, ~Counts,
    "No", "Yes", 4,
    "Yes", "No", 25,
    "Yes", "Yes", 13,
    "No", "No", 92
  )

## plot
ggpiestats(
  data = clinical_trial,
  x = SickAfter,
  y = SickBefore,
  counts = Counts,
  paired = TRUE,
  label = "both",
  title = "Results from imaginary clinical trial",
  package = "ggsci",
  palette = "default_ucscgb"
)
```

The results bode well for our experimental wonder drug. Of the 96 who started
out healthy only 4% were sick after a month. Ideally, we would have hoped for
zero but reality is seldom perfect. On the other side of the 38 who started out
sick that number has reduced to just 13 or 34% which is a marked improvement.

## Summary of graphics

graphical element | `geom_` used | argument for further modification
--------- | ------- | --------------------------
pie slices | `ggplot2::geom_col` | ❌
descriptive labels | `ggplot2::geom_label`/`ggrepel::geom_label_repel` | `label.args`

## Summary of tests

**two-way table**

**Hypothesis testing**

Type | Design | Test | Function used
----------- | ----- | ------------------------- | -----
Parametric/Non-parametric | Unpaired | Pearson's $\chi^2$ test | `stats::chisq.test`
Bayesian | Unpaired | Bayesian Pearson's $\chi^2$ test | `BayesFactor::contingencyTableBF`
Parametric/Non-parametric | Paired  | McNemar's $\chi^2$ test | `stats::mcnemar.test`
Bayesian | Paired  | ❌ | ❌

**Effect size estimation**

Type | Design | Effect size | CI? | Function used
----------- | ----- | ------------------------- | --- | -----
Parametric/Non-parametric | Unpaired | Cramer's $V$ | ✅ | `effectsize::cramers_v`
Bayesian | Unpaired | Cramer's $V$ | ✅ | `effectsize::cramers_v`
Parametric/Non-parametric | Paired | Cohen's $g$ | ✅ | `effectsize::cohens_g`
Bayesian | Paired | ❌ | ❌ | ❌

**one-way table**

**Hypothesis testing**

Type | Test | Function used
----------- | ------------------------- | -----
Parametric/Non-parametric | Goodness of fit $\chi^2$ test | `stats::chisq.test`
Bayesian | Bayesian Goodness of fit $\chi^2$ test | (custom)

**Effect size estimation**

Type | Effect size | CI? | Function used
----------- | ------------------------- | --- | -----
Parametric/Non-parametric | Pearson's $C$ | ✅ | `effectsize::pearsons_c`
Bayesian | ❌ | ❌ | ❌

## Reporting

If you wish to include statistical analysis results in a publication/report, the
ideal reporting practice will be a hybrid of two approaches:

- the `{ggstatsplot}` approach, where the plot contains both the visual and
numerical summaries about a statistical model, and

- the *standard* narrative approach, which provides interpretive context for the
reported statistics.

For example, let's see the following example:

```{r reporting}
ggpiestats(mtcars, am, cyl)
```

The narrative context (assuming `type = "parametric"`) can complement this plot
either as a figure caption or in the main text-

> Pearson's $\chi^2$-test of independence revealed that, across 32 automobiles,
showed that there was a significant association between transmission engine and
number of cylinders. The Bayes Factor for the same analysis revealed that the
data were `r round(exp(2.82), 2)` times more probable under the alternative
hypothesis as compared to the null hypothesis. This can be considered strong
evidence (Jeffreys, 1961) in favor of the alternative hypothesis.

Similar reporting style can be followed when the function performs one-sample
goodness-of-fit test instead of a $\chi^2$-test.

Same holds true for `ggbarstats`.

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
GitHub: <https://github.com/IndrajeetPatil/ggstatsplot/issues>

---
title: "ggcoefstats"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
    error: TRUE
vignette: >
  %\VignetteIndexEntry{ggcoefstats}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

# use max # CPUs
options(mc.cores = parallel::detectCores())

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  error = TRUE,
  package.startup.message = FALSE
)

# for faster post-MCMC computations
future::plan("multicore")
library(ggstatsplot)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment=""}
citation("ggstatsplot")
```

---

Lifecycle:
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

The function `ggcoefstats` generates **dot-and-whisker plots** for regression
models saved in a tidy data frame. The tidy dataframes are prepared using
`parameters::model_parameters`. Additionally, if available, the model summary
indices are also extracted from `performance::model_performance`.

In this vignette, we will see examples of how to use this function. We will try
to cover as many classes of objects as possible. Unfortunately, there is no
single dataset that will be helpful for carrying out all types of regression
analyses and, therefore, we will use various datasets to explore data-specific
hypotheses using regression models.

## General structure of the plots

Although the statistical models displayed in the plot may differ based on the
class of models being investigated, there are few aspects of the plot that will
be invariant across models:

  - The dot-whisker plot contains a dot representing the **estimate** and their
    **confidence intervals** (`95%` is the default). The estimate can either be
    effect sizes (for tests that depend on the `F`-statistic) or regression
    coefficients (for tests with `t`-, $\chi^{2}$-, and `z`-statistic), etc. The
    confidence intervals can sometimes be asymmetric if bootstrapping was used.

  - The label attached to dot will provide more details from the statistical
    test carried out and it will typically contain estimate, statistic, and
    *p*-value.

  - The caption will contain diagnostic information, if available, about models
    that can be useful for model selection: The smaller the Akaike's Information
    Criterion (**AIC**) and the Bayesian Information Criterion (**BIC**) values,
    the "better" the model is.

  - The output of this function will be a `{ggplot2}` object and, thus, it can be
    further modified (e.g., change themes, etc.) with `{ggplot2}` functions.

## Supported models

Most of the regression models that are supported in the underlying packages are
also supported by `ggcoefstats`. 

```{r supported}
insight::supported_models()
```

## Summary of graphics

graphical element | `geom_` used | argument for further modification 
--------- | ------- | -------------------------- 
regression estimate | `ggplot2::geom_point`| `point.args` 
error bars | `ggplot2::geom_errorbarh` | `errorbar.args` 
vertical line | `ggplot2::geom_vline` | `vline.args` 
label with statistical details | `ggrepel::geom_label_repel` | `stats.label.args`

## Summary of meta-analysis tests

**Hypothesis testing** and **Effect size estimation**

Type | Test | Effect size | CI? | Function used 
----------- | -------------------- | -------- | --- | -----
Parametric | Meta-analysis via random-effects models | $\beta$ | ✅ | `metafor::metafor` 
Robust | Meta-analysis via robust random-effects models | $\beta$ | ✅ | `metaplus::metaplus` 
Bayes | Meta-analysis via Bayesian random-effects models | $\beta$ | ✅ | `metaBMA::meta_random`

## Examples of supported models

The following examples are organized by statistics type.

There used to be a much longer vignette with examples of a wide collection of
regression models, but for the sake of maintainability, I have removed it. The
old version can be found [here](https://github.com/IndrajeetPatil/ggstatsplot/blob/master/old/effsize_interpretation.Rmd).

### *t*-statistic

linear model (`lm`) and linear mixed-effects model (`lmer`/`lmerMod`)

```{r lmer1, fig.height=8, fig.width=6}
# set up
library(lme4)
library(ggstatsplot)
set.seed(123)

# lm model
mod1 <- stats::lm(formula = scale(rating) ~ scale(budget), data = movies_long)

# merMod model
mod2 <- lme4::lmer(
  formula = scale(rating) ~ scale(budget) + (budget | genre),
  data = movies_long,
  control = lme4::lmerControl(calc.derivs = FALSE)
)

# combining the two different plots
combine_plots(
  plotlist = list(
    ggcoefstats(mod1) +
      ggplot2::labs(x = parse(text = "'standardized regression coefficient' ~italic(beta)")),
    ggcoefstats(mod2) +
      ggplot2::labs(
        x = parse(text = "'standardized regression coefficient' ~italic(beta)"),
        y = "fixed effects"
      )
  ),
  plotgrid.args = list(nrow = 2),
  annotation.args = list(title = "Relationship between movie budget and its IMDB rating")
)
```

Note that for mixed-effects models, only the *fixed* effects are shown because
there are no confidence intervals for *random* effects terms. In case, you would
like to see these terms, you can enter the same object you entered as `x`
argument to `parameters::model_parameters`:

```{r lmer2}
# setup
set.seed(123)
library(lme4)
library(parameters)

# tidy output
parameters::model_parameters(
  lme4::lmer(
    formula = scale(rating) ~ scale(budget) + (budget | genre),
    data = movies_long,
    control = lme4::lmerControl(calc.derivs = FALSE)
  )
)
```

### *z*-statistic

Aalen's additive regression model for censored data (`aareg`)

```{r aareg, fig.height=5, fig.width=5}
# setup
library(survival)
set.seed(123)

# model
afit <- survival::aareg(
    formula = Surv(time, status) ~ age + sex + ph.ecog,
    data = lung,
    dfbeta = TRUE
  )

# plot
ggcoefstats(
  x = afit,
  title = "Aalen's additive regression model",
  subtitle = "(for censored data)",
  k = 3
)
```

### $\chi^2$-statistic 

Cox proportional hazards regression model (`coxph`)

```{r coxph, fig.height=4, fig.width=4}
# for reproducibility
set.seed(123)
library(survival)

# create the simplest-test data set
test1 <- list(
  time = c(4, 3, 1, 1, 2, 2, 3),
  status = c(1, 1, 1, 0, 1, 1, 0),
  x = c(0, 2, 1, 1, 1, 0, 0),
  sex = c(0, 0, 0, 0, 1, 1, 1)
)

# fit a stratified model
mod_coxph <-
  survival::coxph(
    formula = Surv(time, status) ~ x + strata(sex),
    data = test1
  )

# plot
ggcoefstats(
  x = mod_coxph,
  title = "Cox proportional hazards regression model"
)
```

Another example with `frailty` term.

```{r coxph.penal, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(survival)

# model
mod_coxph <- survival::coxph(
  formula = Surv(time, status) ~ age + sex + frailty(inst),
  data = lung
)

# plot
ggcoefstats(
  x = mod_coxph,
  title = "Proportional Hazards Regression Model\nwith Frailty penalty function"
)
```

### *F*-statistic

omnibus ANOVA (`aov`)

```{r aov1, fig.height=6, fig.width=8}
# setup
set.seed(123)
library(ggstatsplot)
library(ggplot2)

# model
mod_aov <- stats::aov(formula = rating ~ mpaa * genre, data = movies_long)

# plot
ggcoefstats(
  x = mod_aov,
  effsize = "omega", # changing the effect size estimate being displayed
  point.args = list(color = "red", size = 4, shape = 15), # changing the point geom
  package = "dutchmasters", # package from which color palette is to be taken
  palette = "milkmaid", # color palette for labels
  title = "omnibus ANOVA", # title for the plot
  exclude.intercept = TRUE
) +
  # further modification with the ggplot2 commands
  # note the order in which the labels are entered
  ggplot2::scale_y_discrete(labels = c("MPAA", "Genre", "Interaction term")) +
  ggplot2::labs(x = "effect size estimate (eta-squared)", y = NULL)
```

Note that we can also use this function for model selection. You can try out
different models with the code below and see how the AIC and BIC values change.

```{r aov2, fig.height=10, fig.width=10}
# setup
set.seed(123)
library(ggstatsplot)

# plot
combine_plots(
  plotlist = list(
    # model 1
    ggcoefstats(
      x = stats::aov(formula = rating ~ mpaa, data = movies_long),
      title = "1. Only MPAA ratings"
    ),
    # model 2
    ggcoefstats(
      x = stats::aov(formula = rating ~ genre, data = movies_long),
      title = "2. Only genre"
    ),
    # model 3
    ggcoefstats(
      x = stats::aov(formula = rating ~ mpaa + genre, data = movies_long),
      title = "3. Additive effect of MPAA and genre"
    ),
    # model 4
    ggcoefstats(
      x = stats::aov(formula = rating ~ mpaa * genre, data = movies_long),
      title = "4. Multiplicative effect of MPAA and genre"
    )
  ),
  annotation.args = list(title = "Model selection using ggcoefstats")
)
```

Another example with multivariate analysis of variance (`manova`):

```{r manova, fig.height=8, fig.width=8}
# setup
set.seed(123)

# fake a 2nd response variable
npk2 <- within(npk, foo <- rnorm(24))

# model
m_manova <- manova(cbind(yield, foo) ~ block + N * P * K, npk2)

# plot
ggcoefstats(
  x = m_manova,
  title = "multivariate analysis of variance"
)
```

### Bayesian models - no statistic

```{r BFBayesFactor, fig.width=10, fig.height=10}
# setup
set.seed(123)
library(BayesFactor)
library(ggstatsplot)

# one sample t-test
mod1 <- ttestBF(mtcars$wt, mu = 3)

# independent t-test
mod2 <- ttestBF(formula = wt ~ am, data = mtcars)

# paired t-test
mod3 <- ttestBF(x = sleep$extra[1:10], y = sleep$extra[11:20], paired = TRUE)

# correlation
mod4 <- correlationBF(y = iris$Sepal.Length, x = iris$Sepal.Width)

# contingency tabs (not supported)
data("raceDolls")
mod5 <- contingencyTableBF(
  raceDolls,
  sampleType = "indepMulti",
  fixedMargin = "cols"
)

# anova
data("puzzles")
mod6 <- anovaBF(
  formula = RT ~ shape * color + ID,
  data = puzzles,
  whichRandom = "ID",
  whichModels = "top",
  progress = FALSE
)

# regression-1
mod7 <- regressionBF(rating ~ ., data = attitude, progress = FALSE)

# meta-analysis
t <- c(-.15, 2.39, 2.42, 2.43, -.15, 2.39, 2.42, 2.43)
N <- c(100, 150, 97, 99, 99, 97, 100, 150)
mod8 <- meta.ttestBF(t, N, rscale = 1, nullInterval = c(0, Inf))

# proportion test
mod9 <- proportionBF(y = 15, N = 25, p = .5)

# list of plots
combine_plots(
  plotlist = list(
    ggcoefstats(mod1, title = "one sample t-test"),
    ggcoefstats(mod2, title = "independent t-test"),
    ggcoefstats(mod3, title = "paired t-test"),
    ggcoefstats(mod4, title = "correlation"),
    ggcoefstats(mod5, title = "contingency table"),
    ggcoefstats(mod6, title = "anova"),
    ggcoefstats(mod7, title = "regression-1"),
    ggcoefstats(mod8, title = "meta-analysis"),
    ggcoefstats(mod9, title = "proportion test")
  ),
  annotation.args = list(title = "Example from `BayesFactor` package")
)
```

### Regression models with `list` outputs

Note that a number of regression models will return an object of class `list`,
in which case this function will fail. But often you can extract the object of
interest from this list and use it to plot the regression coefficients.

```r
# setup
library(gamm4)
set.seed(123)

# data
dat <- gamSim(1, n = 400, scale = 2)

# now add 20 level random effect `fac'...
dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5

# model object
br <- gamm4::gamm4(
    formula = y ~ s(x0) + x1 + s(x2),
    data = dat,
    random = ~ (1 | fac)
  )

# looking at the classes of the objects contained in the list
purrr::map(br, class)

# plotting
combine_plots(
  plotlist = list(
    # first object plot (only parametric terms are shown)
    ggcoefstats(
      x = br$gam,
      title = "generalized additive model (parametric terms)",
      k = 3
    ),
    # second object plot
    ggcoefstats(
      x = br$mer,
      title = "linear mixed-effects model",
      k = 3
    )
  ),
  plotgrid.args = list(nrow = 1)
)
```

## Meta-analysis

In case the estimates you are displaying come from multiple studies, you can
also use this function to carry out random-effects meta-analysis. The dataframe
you enter **must** contain at the minimum the following three columns-

  - `term`: a column with names/identifiers to annotate each study/effect

  - `estimate`: a column with the observed effect sizes or outcomes

  - `std.error`: a column the corresponding standard errors

### parametric

```{r meta1, fig.height=7, fig.width=8}
# setup
set.seed(123)
library(metaplus)

# renaming to what the function expects
df <- dplyr::rename(mag, estimate = yi, std.error = sei, term = study)

# plot
ggcoefstats(
  x = df,
  meta.analytic.effect = TRUE,
  bf.message = TRUE,
  meta.type = "parametric",
  title = "parametric random-effects meta-analysis"
)
```

### robust

```{r meta2, fig.height=7, fig.width=8}
# setup
set.seed(123)
library(metaplus)

# renaming to what the function expects
df <- dplyr::rename(mag, estimate = yi, std.error = sei, term = study)

# plot
ggcoefstats(
  x = df,
  meta.analytic.effect = TRUE,
  meta.type = "robust",
  title = "robust random-effects meta-analysis"
)
```

### Bayesian

```{r meta3, fig.height=7, fig.width=8}
# setup
set.seed(123)
library(metaplus)

# renaming to what the function expects
df <- dplyr::rename(mag, estimate = yi, std.error = sei, term = study)

# plot
ggcoefstats(
  x = df,
  meta.analytic.effect = TRUE,
  meta.type = "bayes",
  title = "Bayesian random-effects meta-analysis"
)
```

## Dataframes

Sometimes you don't have a model object but a custom dataframe that you want
display using this function. If a data frame is to be plotted, it **must**
contain columns named `term` (names of predictors), and `estimate`
(corresponding estimates of coefficients or other quantities of interest). Other
optional columns are `conf.low` and `conf.high` (for confidence intervals), and
`p.value`. You will also have to specify the type of statistic relevant for
regression models (`"t"`, `"z"`, `"f"`, `"chi"`) in case you want to display
statistical labels.

You can also provide a dataframe containing all the other relevant information
for additionally displaying labels with statistical information.

```{r dataframe, fig.height=7, fig.width=7}
# let's make up a dataframe (with all available details)
df_full <-
  tibble::tribble(
    ~term, ~statistic, ~estimate, ~std.error, ~p.value, ~df.error,
    "study1", 0.158, 0.0665, 0.778, 0.875, 5L,
    "study2", 1.33, 0.542, 0.280, 0.191, 10L,
    "study3", 1.24, 0.045, 0.030, 0.001, 12L,
    "study4", 0.156, 0.500, 0.708, 0.885, 8L,
    "study5", 0.33, 0.032, 0.280, 0.101, 2L,
    "study6", 1.04, 0.085, 0.030, 0.001, 3L
  )

# plot
ggcoefstats(
  x = df_full,
  meta.analytic.effect = TRUE,
  statistic = "t",
  package = "LaCroixColoR",
  palette = "paired"
)
```

## Non-plot outputs

This function can also be used to extract outputs other than a plot, although it
is much more preferable to use the underlying functions instead
(`parameters::model_parameters`).

```{r other_output, error=TRUE}
# setup
set.seed(123)
library(ggstatsplot)

# data
DNase1 <- subset(DNase, Run == 1)

# using a selfStart model
nlmod <- stats::nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)

# tidy dataframe
ggcoefstats(nlmod, output = "tidy")

# glance summary
ggcoefstats(nlmod, output = "glance")
```

## Not supported

This vignette was supposed to give a comprehensive account of regression models
supported by `ggcoefstats`. The list of supported models will keep expanding as
additional tidiers are added to the `parameters` and `performance` packages.

Note that not **all** models supported in these packages will be supported by
`ggcoefstats`. In particular, classes of objects for which there is no column
for `estimate` (e.g., `kmeans`, `optim`, `muhaz`, `survdiff`, `zoo`, etc.) are
not supported.

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "combine_plots"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{combine_plots}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

Lifecycle: [![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)

# Combining plots with `combine_plots`

The full power of `{ggstatsplot}` can be leveraged with a functional programming
package like [`{purrr}`](https://purrr.tidyverse.org/) that replaces `for` loops
with code that is both more succinct and easier to read and, therefore, `{purrr}`
should be preferrred 😻.

In such cases, `{ggstatsplot}` contains a helper function `combine_plots` to
combine multiple plots, which can be useful for combining a list of plots
produced with `{purrr}`. This is a wrapper around `patchwork::wrap_plots` and lets
you combine multiple plots and add a combination of title, caption, and
annotation texts with suitable defaults.

**Note Before**: If you have just one grouping variable and you'd like a plot
for each `factor` of this variable the `grouped_` variants
(<https://indrajeetpatil.github.io/ggstatsplot/reference/index.html>) of all
`{ggstatsplot}` functions will allow you do to this. They specifically use the
`combine_plots` function under the covers.

# Example-1 using `dplyr::group_map`

The easiest way to run the same `{ggstatsplot}` operation across multiple grouping
variables is by using `dplyr::group_map` functions and then - of course - one
would like to combine these plots in a single plot.

```{r ggscatterstats_groupmap, fig.height = 6, fig.width = 10}
# libraries
set.seed(123)
library(tidyverse)
library(ggstatsplot)

# creating a list of plots
p_list <-
  mtcars %>%
  dplyr::filter(cyl != 4) %>%
  dplyr::group_by(cyl) %>%
  dplyr::group_map(.f = ~ ggbetweenstats(data = ., x = am, y = wt))

# combining plots
combine_plots(
  plotlist = p_list,
  annotation.args = list(tag_levels = list(as.vector(rlang::set_names(levels(as.factor(mtcars$cyl))))))
)
```

# Example-2 using `{purrr}`

The full power of `{ggstatsplot}` can be leveraged with a functional programming
package like [`{purrr}`](http://purrr.tidyverse.org/) which can replace many `for`
loops, is more succinct, and easier to read. Consider `{purrr}` as your first
choice for combining multiple plots.

An example using the `iris` dataset is provided below.  Imagine that we want to
separately plot the linear relationship between sepal length and sepal width for
each of the three species but combine them into one consistent plot with common
labeling and as one plot.  Rather than call `ggscatterstats` three times and
gluing the results or using `patchwork` directly, we'll create a `tibble` called
`plots` using `purrr::map` then feed that to `combine_plots` to get our combined
plot.

```{r ggscatterstats_purrr, fig.height = 12, fig.width = 8}
# for reproducibility
set.seed(123)
library(ggstatsplot)

# creating a list column with `{ggstatsplot}` plots
plots <- iris %>%
  dplyr::mutate(Species2 = Species) %>%
  # just creates a copy of this variable
  dplyr::group_by(Species) %>%
  tidyr::nest(.) %>%
  # a nested dataframe with list column called `data`
  dplyr::mutate(
    plot = data %>%
      purrr::map(
        .x = .,
        .f = ~ ggscatterstats(
          data = .,
          x = Sepal.Length,
          y = Sepal.Width,
          title = glue::glue("Species: {.$Species2} (n = {length(.$Sepal.Length)})")
        )
      )
  )

# display the new object 
plots

# creating a grid with patchwork
combine_plots(
  plotlist = plots$plot,
  plotgrid.args = list(nrow = 3, ncol = 1),
  annotation.args = list(
    title = "Relationship between sepal length and width for each Iris species",
    caption = expression(
      paste(
        italic("Note"),
        ": Iris flower dataset was collected by Edgar Anderson.",
        sep = ""
      )
    )
  )
)
```

# Example-3 with `plyr`

Another popular package for handling big datasets is `plyr`, which allows us to
repeatedly apply a common function on smaller pieces and then combine the
results into a larger whole.

In this example we'll start with the `gapminder` dataset.  We're interested in
the linear relationship between Gross Domestic Product (per capita) and life
expectancy in the year 2007, for all the continents except Oceania. We'll use
`dplyr` to filter to the right rows then use `plyr` to repeat the
`ggscatterstats` function across each of the 4 continents remaining.  The result
is of that is a list of plots called `plots`.  We then feed `plots` to the
`combine_plots` function to merge them into one plot. We will call attention to
the countries which have very low life expectancy (< 45 years) by labeling those
countries when they occur.

```{r ggscatterstats_plyr, fig.height = 12, fig.width = 12}
library(plyr)
library(gapminder)

# for reproducibility
set.seed(123)

# let's have a look at the structure of the data
dplyr::glimpse(gapminder::gapminder)

# creating a list of plots
plots <- plyr::dlply(
  dplyr::filter(gapminder::gapminder, year == 2007, continent != "Oceania"),
  .variables = .(continent),
  .fun = function(data) {
    ggscatterstats(
      data = data,
      x = gdpPercap,
      y = lifeExp,
      xfill = "#0072B2",
      yfill = "#009E73",
      label.var = "country",
      label.expression = "lifeExp < 45",
      title = glue::glue("Continent: {data$continent}")
    ) +
      ggplot2::scale_x_continuous(labels = scales::comma)
  }
)

# combining individual plots
combine_plots(
  plotlist = plots,
  annotation.args = list(title = "Relationship between GDP (per capita) and life expectancy"),
  plotgrid.args = list(nrow = 2)
)
```

# Adding additional details to plots

The `combine_plots` function can also be useful for adding additional textual
information that can not be added by making a single call to a `{ggstatsplot}`
function via the title, subtitle, or caption options. For this example let's
assume we want to assess the relationship between a movie's rating and its
budget from the **Internet Movie Database** using polynomial regression.

`ggcoefstats` will do most of the work, including the title, subtitle, and
caption. But we want to add at the bottom an annotation to show the formula we
are using for our regression.  `combine_plots` allows us to add `sub.text =` to
accomplish that task as shown in the resulting plot.

```{r ggbetweenstats_subtext, fig.height = 8, fig.width = 8}
library(ggstatsplot)

combine_plots(
  plotlist = list(ggcoefstats(
    x = stats::lm(
      formula = rating ~ stats::poly(budget, degree = 3),
      data = dplyr::sample_frac(movies_long, size = 0.2),
      na.action = na.omit
    ),
    exclude.intercept = FALSE,
    title = "Relationship between movie budget and IMDB rating",
    subtitle = "Source: Internet Movie Database",
    ggtheme = ggplot2::theme_gray(),
    stats.label.color = c("#CC79A7", "darkgreen", "#0072B2", "red")
  ) +
    # modifying the plot outside of ggstatsplot using ggplot2 functions
    ggplot2::scale_y_discrete(
      labels = c(
        "Intercept (c)",
        "1st degree (b1)",
        "2nd degree (b2)",
        "3rd degree (b3)"
      )
    ) +
    ggplot2::labs(y = "term (polynomial regression)")),
  # adding additional text element to the plot since title, subtitle, caption are all already occupied
  annotation.args = list(
    caption = expression(
      paste(
        "linear model: ", bolditalic(y),
        " ~ ",
        italic(c) + italic(b)[1] * bolditalic(x) + italic(b)[2] * bolditalic(x)^
          2 + italic(b)[3] * bolditalic(x)^3,
        sep = ""
      )
    )
  )
)
```

# Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on `GitHub`: 
<https://github.com/IndrajeetPatil/ggstatsplot/issues>

# Session Information

For details, see-
<https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/session_info.html>
---
title: "ggcoefstats"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
    error: TRUE
vignette: >
  %\VignetteIndexEntry{ggcoefstats}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

# use max # CPUs
options(mc.cores = parallel::detectCores())

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  error = TRUE,
  package.startup.message = FALSE
)

# for faster post-MCMC computations
future::plan("multicore")
library(ggstatsplot)
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment=""}
citation("ggstatsplot")
```

---

Lifecycle:
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

The function `ggcoefstats` generates **dot-and-whisker plots** for regression
models saved in a tidy data frame. The tidy dataframes are prepared using
`parameters::model_parameters`. Additionally, if available, the model summary
indices are also extracted from `performance::model_performance`.

In this vignette, we will see examples of how to use this function. We will try
to cover as many classes of objects as possible. Unfortunately, there is no
single dataset that will be helpful for carrying out all types of regression
analyses and, therefore, we will use various datasets to explore data-specific
hypotheses using regression models.

# General structure of the plots

Although the statistical models displayed in the plot may differ based on the
class of models being investigated, there are few aspects of the plot that will
be invariant across models:

  - The dot-whisker plot contains a dot representing the **estimate** and their
    **confidence intervals** (`95%` is the default). The estimate can either be
    effect sizes (for tests that depend on the `F`-statistic) or regression
    coefficients (for tests with `t`-, $\chi^{2}$-, and `z`-statistic), etc. The
    confidence intervals can sometimes be asymmetric if bootstrapping was used.

  - The label attached to dot will provide more details from the statistical
    test carried out and it will typically contain estimate, statistic, and
    *p*-value.

  - The caption will contain diagnostic information, if available, about models
    that can be useful for model selection: The smaller the Akaike's Information
    Criterion (**AIC**) and the Bayesian Information Criterion (**BIC**) values,
    the "better" the model is.

  - The output of this function will be a `{ggplot2}` object and, thus, it can be
    further modified (e.g., change themes, etc.) with `{ggplot2}` functions.

# Supported models

Most of the regression models that are supported in the underlying packages are
also supported by `ggcoefstats`. 

```{r supported}
insight::supported_models()
```

# Summary of graphics

graphical element | `geom_` used | argument for further modification 
--------- | ------- | -------------------------- 
regression estimate | `ggplot2::geom_point`| `point.args` 
error bars | `ggplot2::geom_errorbarh` | `errorbar.args` 
vertical line | `ggplot2::geom_vline` | `vline.args` 
label with statistical details | `ggrepel::geom_label_repel` | `stats.label.args`

# Summary of meta-analysis tests

**Hypothesis testing** and **Effect size estimation**

Type | Test | Effect size | CI? | Function used 
----------- | -------------------- | -------- | --- | -----
Parametric | Meta-analysis via random-effects models | $\beta$ | ✅ | `metafor::metafor` 
Robust | Meta-analysis via robust random-effects models | $\beta$ | ✅ | `metaplus::metaplus` 
Bayes | Meta-analysis via Bayesian random-effects models | $\beta$ | ✅ | `metaBMA::meta_random`

# Examples of supported models

First let's load the needed library.

```{r log_ggstats}
library(ggstatsplot)
```

The following demos are in **no particular order**.

## linear mixed-effects model (`lmer`/`lmerMod`)

```{r lmer1, fig.height=8, fig.width=6}
# set up
library(lme4)
library(ggstatsplot)
set.seed(123)

# lm model
mod1 <- stats::lm(formula = scale(rating) ~ scale(budget), data = movies_long)

# merMod model
mod2 <- lme4::lmer(
  formula = scale(rating) ~ scale(budget) + (budget | genre),
  data = movies_long,
  control = lme4::lmerControl(calc.derivs = FALSE)
)

# combining the two different plots
combine_plots(
  plotlist = list(
    # model 1: simple linear model
    ggcoefstats(
      x = mod1,
      title = "linear model",
      exclude.intercept = TRUE # hide the intercept
    ) +
      ggplot2::labs(x = parse(text = "'standardized regression coefficient' ~italic(beta)")),
    # model 2: linear mixed-effects model
    ggcoefstats(
      x = mod2,
      title = "linear mixed-effects model",
      exclude.intercept = TRUE # hide the intercept
    ) +
      ggplot2::labs(
        x = parse(text = "'standardized regression coefficient' ~italic(beta)"),
        y = "fixed effects"
      )
  ),
  plotgrid.args = list(nrow = 2),
  annotation.args = list(title = "Relationship between movie budget and its IMDB rating")
)
```

Note that for mixed-effects models, only the *fixed* effects are shown because
there are no confidence intervals for *random* effects terms. In case, you would
like to see these terms, you can enter the same object you entered as `x`
argument to `parameters::model_parameters`:

```{r lmer2}
# setup
set.seed(123)
library(lme4)
library(parameters)

# tidy output
parameters::model_parameters(
  lme4::lmer(
    formula = scale(rating) ~ scale(budget) + (budget | genre),
    data = movies_long,
    control = lme4::lmerControl(calc.derivs = FALSE)
  )
)
```

## summary grid from `emmeans` package (`emmGrid`)

```{r emmGrid, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(emmeans)

# linear model for sales of oranges per day
oranges_lm1 <-
  stats::lm(
    formula = sales1 ~ price1 + price2 + day + store,
    data = oranges
  )

# reference grid; see vignette("basics", package="emmeans")
oranges_rg1 <- emmeans::ref_grid(oranges_lm1)
marginal <- emmeans::emmeans(oranges_rg1, "day")

# plot
ggcoefstats(
  x = marginal,
  point.args = list(color = "darkgreen", shape = 9),
  title = "summary grid from `emmeans` package"
)
```

Another example with output containing effect sizes

```{r emmGrid2, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(emmeans)

# model and effect sizes
fiber.lm <- lm(strength ~ diameter + machine, data = fiber)
emm <- emmeans::emmeans(fiber.lm, "machine")
es <- emmeans::eff_size(emm, sigma = sigma(fiber.lm), edf = df.residual(fiber.lm))

# plot
ggcoefstats(
  x = es,
  title = "summary grid with effect sizes\nfrom `emmeans` package"
)
```

## Bayesian additive models for location scale and shape (`bamlss`)

```{r bamlss, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(bamlss)
data("SwissLabor", package = "AER")

# formula
f <- participation ~ income + age + education + youngkids + oldkids + foreign + I(age^2)

# model
b_bamlss <-
  bamlss::bamlss(
    formula = f,
    family = "binomial",
    data = SwissLabor
  )

# plot
ggcoefstats(
  x = b_bamlss,
  title = "Bayesian additive models \nfor location scale and shape"
)
```

## Bayes Factor (`BFBayesFactor`)

```{r BFBayesFactor, fig.width=10, fig.height=10}
# setup
set.seed(123)
library(BayesFactor)
library(ggstatsplot)

# one sample t-test
mod1 <- ttestBF(mtcars$wt, mu = 3)

# independent t-test
mod2 <- ttestBF(formula = wt ~ am, data = mtcars)

# paired t-test
mod3 <- ttestBF(x = sleep$extra[1:10], y = sleep$extra[11:20], paired = TRUE)

# correlation
mod4 <- correlationBF(y = iris$Sepal.Length, x = iris$Sepal.Width)

# contingency tabs (not supported)
data("raceDolls")
mod5 <- contingencyTableBF(
  raceDolls,
  sampleType = "indepMulti",
  fixedMargin = "cols"
)

# anova
data("puzzles")
mod6 <- anovaBF(
  formula = RT ~ shape * color + ID,
  data = puzzles,
  whichRandom = "ID",
  whichModels = "top",
  progress = FALSE
)

# regression-1
mod7 <- regressionBF(rating ~ ., data = attitude, progress = FALSE)

# meta-analysis
t <- c(-.15, 2.39, 2.42, 2.43, -.15, 2.39, 2.42, 2.43)
N <- c(100, 150, 97, 99, 99, 97, 100, 150)
mod8 <- meta.ttestBF(t, N, rscale = 1, nullInterval = c(0, Inf))

# proportion test
mod9 <- proportionBF(y = 15, N = 25, p = .5)

# list of plots
combine_plots(
  plotlist = list(
    ggcoefstats(mod1, title = "one sample t-test"),
    ggcoefstats(mod2, title = "independent t-test"),
    ggcoefstats(mod3, title = "paired t-test"),
    ggcoefstats(mod4, title = "correlation"),
    ggcoefstats(mod5, title = "contingency table"),
    ggcoefstats(mod6, title = "anova"),
    ggcoefstats(mod7, title = "regression-1"),
    ggcoefstats(mod8, title = "meta-analysis"),
    ggcoefstats(mod9, title = "proportion test")
  ),
  annotation.args = list(title = "Example from `BayesFactor` package")
)
```

## Fixed Effects Individual Slope Estimator (`feis`)

```{r feis, fig.height=8, fig.width=6}
# setup
set.seed(123)
library(feisr)
data("mwp", package = "feisr")

# model
feis.mod <-
  feisr::feis(
    formula = lnw ~ marry + enrol + as.factor(yeargr) | exp + I(exp^2),
    data = mwp,
    id = "id",
    robust = TRUE
  )

# plot
ggcoefstats(
  x = feis.mod,
  title = "Fixed Effects Individual Slope Estimator"
)
```

## omnibus ANOVA (`aov`)

```{r aov1, fig.height=6, fig.width=8}
# setup
set.seed(123)
library(ggstatsplot)
library(ggplot2)

# model
mod_aov <- stats::aov(formula = rating ~ mpaa * genre, data = movies_long)

# plot
ggcoefstats(
  x = mod_aov,
  effsize = "omega", # changing the effect size estimate being displayed
  point.args = list(color = "red", size = 4, shape = 15), # changing the point geom
  package = "dutchmasters", # package from which color palette is to be taken
  palette = "milkmaid", # color palette for labels
  title = "omnibus ANOVA", # title for the plot
  exclude.intercept = TRUE
) +
  # further modification with the ggplot2 commands
  # note the order in which the labels are entered
  ggplot2::scale_y_discrete(labels = c("MPAA", "Genre", "Interaction term")) +
  ggplot2::labs(x = "effect size estimate (eta-squared)", y = NULL)
```

Note that we can also use this function for model selection. You can try out
different models with the code below and see how the AIC and BIC values change.

```{r aov2, fig.height=10, fig.width=10}
# setup
set.seed(123)
library(ggstatsplot)

# plot
combine_plots(
  plotlist = list(
    # model 1
    ggcoefstats(
      x = stats::aov(formula = rating ~ mpaa, data = movies_long),
      title = "1. Only MPAA ratings"
    ),
    # model 2
    ggcoefstats(
      x = stats::aov(formula = rating ~ genre, data = movies_long),
      title = "2. Only genre"
    ),
    # model 3
    ggcoefstats(
      x = stats::aov(formula = rating ~ mpaa + genre, data = movies_long),
      title = "3. Additive effect of MPAA and genre"
    ),
    # model 4
    ggcoefstats(
      x = stats::aov(formula = rating ~ mpaa * genre, data = movies_long),
      title = "4. Multiplicative effect of MPAA and genre"
    )
  ),
  annotation.args = list(title = "Model selection using ggcoefstats")
)
```

## conditional logit models (`mclogit`)

```{r mclogit, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(mclogit)
data(Transport)

# model
mod_mclogit <-
  mclogit::mclogit(
    formula = cbind(resp, suburb) ~ distance + cost,
    data = Transport,
    control = mclogit::mclogit.control(trace = FALSE)
  )

# plot
ggcoefstats(
  x = mod_mclogit,
  title = "conditional logit models"
)
```

## mixed conditional logit models (`mmclogit`)

```{r mmclogit, fig.height=8, fig.width=7}
# setup
set.seed(123)
library(mclogit)
data(electors)

# model
mod_mclogit <-
  mclogit::mclogit(
    formula = cbind(Freq, interaction(time, class)) ~
    econ.left / class + welfare / class + auth / class,
    random = ~ 1 | party.time,
    data = within(electors, party.time <- interaction(party, time))
  )

# plot
ggcoefstats(
  x = mod_mclogit,
  title = "Mixed Conditional Logit Models"
)
```

## anova with `car` package (`Anova`)

```{r Anova, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(car)
library(ggstatsplot)

# model
mod_Anova <-
  car::Anova(stats::lm(
    formula = conformity ~ fcategory * partner.status,
    data = Moore,
    contrasts = list(fcategory = contr.sum, partner.status = contr.sum)
  ))

# plot
ggcoefstats(
  x = mod_Anova,
  title = "Anova with `car`"
)
```

## ANOVA with `car` package on multiple linear models (`Anova.mlm`)

```{r Anova.mlm, fig.height=6, fig.width=6, eval=FALSE}
# setup
set.seed(123)
library(car)

# data
dv <- c(1, 3, 4, 2, 2, 3, 2, 5, 6, 3, 4, 4, 3, 5, 6)
subject <- factor(c(
  "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3",
  "s4", "s4", "s4", "s5", "s5", "s5"
))
myfactor <- factor(c(
  "f1", "f2", "f3", "f1", "f2", "f3", "f1", "f2", "f3",
  "f1", "f2", "f3", "f1", "f2", "f3"
))
mydata <- data.frame(dv, subject, myfactor)

dvm <- with(mydata, cbind(
  dv[myfactor == "f1"],
  dv[myfactor == "f2"], dv[myfactor == "f3"]
))

# model
mlm1 <- lm(dvm ~ 1)

rfactor <- factor(c("f1", "f2", "f3"))

# model
mlm1.aov <-
  car::Anova(mlm1,
    idata = data.frame(rfactor),
    idesign = ~rfactor,
    type = "III"
  )

# plot
ggcoefstats(
  x = mlm1.aov,
  title = "ANOVA with `car` on multiple linear models"
)
```

## Anova with `ez` package

```{r anova_ez, fig.height=4, fig.width=4}
set.seed(123)
library(ez)
data(ANT)

# run an ANOVA on the mean correct RT data.
rt_anova <-
  suppressWarnings(ez::ezANOVA(
    data = ANT[ANT$error == 0, ],
    dv = rt,
    wid = subnum,
    within = cue,
    detailed = TRUE,
    return_aov = TRUE
  ))

# plot
ggcoefstats(
  x = rt_anova$aov,
  title = "Anova with `ez` package"
)
```

## linear model (`lm`)

```{r lm, fig.height=8, fig.width=8}
# setup
set.seed(123)

# data
df <- dplyr::filter(
  movies_long,
  genre %in% c(
    "Action",
    "Action Comedy",
    "Action Drama",
    "Comedy",
    "Drama",
    "Comedy Drama"
  )
)

# plot
ggcoefstats(
  x = stats::lm(formula = rating ~ genre, data = df),
  sort = "ascending", # sorting the terms of the model based on estimate values
  ggtheme = ggplot2::theme_gray(), # changing the default theme
  stats.label.color = c("#CC79A7", "darkgreen", "#0072B2", "darkred", "black", "red"),
  title = "Movie ratings by their genre",
  subtitle = "Source: www.imdb.com"
)
```

The same output will also be returned for objects of type `summary.lm`.

## standardized regression coefficients with `lm` (`lm.beta`)

```{r lm.beta, fig.width=10, fig.height=6}
# setup
set.seed(123)
library(MASS)
library(lm.beta)

# model
mod.lm <- stats::lm(formula = mpg ~ wt, data = mtcars)
mod.robust <- MASS::rlm(formula = mpg ~ wt, data = mtcars)
mod.beta <- lm.beta::lm.beta(mod.lm)
mod.robust.beta <- lm.beta::lm.beta(mod.robust)

# plot
combine_plots(
  plotlist = list(
    ggcoefstats(mod.beta),
    ggcoefstats(mod.robust.beta, stats.labels = FALSE)
  ),
  annotation.args = list(title = "standardized regression coefficients with `lm`")
)
```

## bounded memory linear regression (`biglm`)

```{r biglm, fig.height=4, fig.width=6}
# setup
set.seed(123)
library(biglm)

# model
bfit1 <-
  biglm(
    formula = scale(mpg) ~ scale(wt) + scale(disp),
    data = mtcars
  )

# plot
ggcoefstats(
  x = bfit1,
  title = "bounded memory simple linear regression"
)
```

## bounded memory general linear regression (`bigglm`)

```{r bigglm, fig.height=4, fig.width=6}
# setup
set.seed(123)
library(biglm)
data(trees)

# model
mod_bigglm <-
  biglm::bigglm(
    formula = log(Volume) ~ log(Girth) + log(Height),
    data = trees,
    chunksize = 10,
    sandwich = TRUE
  )

# plot
ggcoefstats(
  x = mod_bigglm,
  title = "bounded memory general linear regression"
)
```

## efficient (general) linear model (`elm`/`eglm`)

```{r eflm, fig.height=5, fig.width=12}
# setup
set.seed(123)
library(eflm)

# models
mod_elm <- eflm::elm(mpg ~ wt, data = mtcars)
mod_eglm <- eflm::eglm(mpg ~ wt, family = gaussian(), data = mtcars)

# plot
combine_plots(
  plotlist = list(
    ggcoefstats(
      x = mod_elm,
      title = "efficient linear model"
    ),
    ggcoefstats(
      x = mod_eglm,
      title = "efficient general linear model"
    )
  ),
  plotgrid.args = list(nrow = 1),
  annotation.args = list(title = "efficient (general) linear model using `eflm`")
)
```

## nonlinear model using generalized least squares (`gnls`)

```{r gnls, fig.height=6, fig.width=8}
set.seed(123)
library(nlme)

# variance increases with a power of the absolute fitted values
mod_gnls <- gnls(
  model = weight ~ SSlogis(Time, Asym, xmid, scal),
  data = Soybean,
  weights = varPower()
)

ggcoefstats(
  x = mod_gnls,
  title = "nonlinear model using generalized least squares"
)
```

## analysis of factorial experiments (`mixed`)

```{r mixed, fig.height=4, fig.width=8}
# setup
set.seed(123)
library(afex)
library(MEMSS)
data("Machines", package = "MEMSS")

# simple model with random-slopes for repeated-measures factor
m1_afex <-
  afex::mixed(
    formula = score ~ Machine + (Machine | Worker),
    data = Machines
  )

# suppress correlations among random effect parameters with || and expand_re=TRUE
m2_afex <-
  afex::mixed(
    formula = score ~ Machine + (Machine || Worker),
    data = Machines,
    expand_re = TRUE
  )

# plot
combine_plots(
  plotlist = list(
    ggcoefstats(m1_afex, title = "example-1"),
    ggcoefstats(m2_afex, title = "example-2")
  ),
  annotation.args = list(title = "analysis of factorial experiments (using `afex`)")
)
```

## anova using `afex` (`afex_aov`)

```{r afex_aov, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(afex)
data(obk.long)

# model
fit_all <-
  afex::aov_ez(
    "id",
    "value",
    obk.long,
    between = c("treatment"),
    within = c("phase")
  )

# plot
ggcoefstats(
  x = fit_all,
  title = "anova using `afex`"
)
```

## joint model for survival and longitudinal data measured with error (`joint`)

```{r joint, eval=FALSE}
# setup
set.seed(123)
library(joineR)

# data
data(heart.valve)
heart.surv <- UniqueVariables(heart.valve,
  var.col = c("fuyrs", "status"),
  id.col = "num"
)
heart.long <- heart.valve[, c("num", "time", "log.lvmi")]
heart.cov <- UniqueVariables(heart.valve,
  c("age", "hs", "sex"),
  id.col = "num"
)
heart.valve.jd <- jointdata(
  longitudinal = heart.long,
  baseline = heart.cov,
  survival = heart.surv,
  id.col = "num",
  time.col = "time"
)

# model
mod_joint <- joint(
  data = heart.valve.jd,
  long.formula = log.lvmi ~ 1 + time + hs,
  surv.formula = Surv(fuyrs, status) ~ hs,
  model = "intslope"
)

# plot
ggcoefstats(
  x = mod_joint,
  title = "joint model for survival and longitudinal data measured with error"
)
```

## risk regression (`riskRegression`)

```{r riskRegression, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(riskRegression)
library(prodlim)
data(Melanoma, package = "riskRegression")

# tumor thickness on the log-scale
Melanoma$logthick <- log(Melanoma$thick)

# absolute risk model
multi.arr <-
  riskRegression::ARR(
    formula = Hist(time, status) ~ logthick + sex + age + ulcer,
    data = Melanoma,
    cause = 1
  )

# plot
ggcoefstats(
  x = multi.arr,
  title = "risk regression"
)
```

## Aalen's additive regression model for censored data (`aareg`)

```{r aareg, fig.height=5, fig.width=5}
# setup
library(survival)
set.seed(123)

# model
afit <-
  survival::aareg(
    formula = Surv(time, status) ~ age + sex + ph.ecog,
    data = lung,
    dfbeta = TRUE
  )

# plot
ggcoefstats(
  x = afit,
  title = "Aalen's additive regression model",
  subtitle = "(for censored data)",
  k = 3
)
```

## multivariate generalized linear mixed models (`MCMCglmm`)

```{r MCMCglmm, fig.height=4, fig.width=6}
# setup
set.seed(123)
library(lme4)
library(MCMCglmm)
data(sleepstudy)

# model
mm0 <-
  MCMCglmm::MCMCglmm(
    fixed = scale(Reaction) ~ scale(Days),
    random = ~Subject,
    data = lme4::sleepstudy,
    nitt = 4000,
    pr = TRUE,
    verbose = FALSE
  )

# plot
ggcoefstats(
  x = mm0,
  title = "multivariate generalized linear mixed model",
  conf.method = "HPDinterval",
  robust = TRUE # additional arguments passed to `parameters::model_parameters`
)
```

## STAR Models with BayesX (`bayesx`)

```{r bayesx, fig.height=4, fig.width=5}
# setup
set.seed(111)
library(R2BayesX)

## generate some data
n <- 200

## regressor
dat <- data.frame(x = runif(n, -3, 3))

## response
dat$y <- with(dat, 1.5 + sin(x) + rnorm(n, sd = 0.6))

## estimate models with bayesx REML and MCMC
b1 <- R2BayesX::bayesx(y ~ sx(x), method = "REML", data = dat)

# plot
ggcoefstats(
  x = b1,
  title = "STAR Models with BayesX"
)
```

## Generalised Joint Regression Models with Binary/Continuous/Discrete/Survival Margins (`gjrm`)

```{r gjrm, width=7, height=7}
library(GJRM)

# data
set.seed(123)
dat <- data.frame(
  x1 = rnorm(40, 1:10),
  x2 = rnorm(40, 1:30),
  bid1 = sample(1:5, 40, replace = TRUE),
  bid2 = sample(1:5, 40, replace = TRUE),
  y1 = sample(0:1, 40, replace = TRUE),
  y2 = sample(0:1, 40, replace = TRUE)
)

f.list <- list(
  y1 ~ bid1 + x1 + x2,
  y2 ~ bid2 + x1 + x2
)

# model
mod_gjrm <- gjrm(f.list, dat, Model = "B", margins = c("probit", "probit"))

# plot
ggcoefstats(
  x = mod_gjrm,
  title = "Generalised Joint Regression Models\n with Binary/Continuous/Discrete/Survival Margins"
)
```

## Heckman-style selection and treatment effect models (`selection`)

```{r selection, width=6, height=6}
library(sampleSelection)
library(wooldridge) # for data
data(mroz)

# model
mod_selection <- selection(
  inlf ~ educ + kidslt6 + kidsge6,
  wage ~ educ + exper + expersq,
  method = "2step",
  data = mroz
)

# plot
ggcoefstats(
  x = mod_selection,
  title = "Heckman-style selection and treatment effect models"
)
```

## inference in spatial GLMMs (`HLfit`)

```{r HLfit, width=4, height=4}
# setup
set.seed(123)
library(spaMM)
data("wafers")
data("scotlip")

# model
mod_HLfit <-
  fitme(
    formula = y ~ 1 + (1 | batch),
    family = Gamma(log),
    data = wafers
  )

# plot
ggcoefstats(
  x = mod_HLfit,
  title = "Inference in spatial GLMMs"
)
```

## robust linear mixed-effects models (`rlmer`)

```{r rlmer, width=5, height=5}
# setups
set.seed(123)
library(robustlmm)

# model
roblmm.mod <-
  robustlmm::rlmer(
    formula = scale(Reaction) ~ scale(Days) + (Days | Subject),
    data = sleepstudy,
    rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
    rho.sigma.b = chgDefaults(smoothPsi, k = 5.11, s = 10)
  )

# plot
ggcoefstats(
  x = roblmm.mod,
  title = "robust estimation of linear mixed-effects model",
  conf.level = 0.90
)
```

## linear mixed-effects models with `lmerTest` (`lmerModLmerTest`)

```{r lmerModLmerTest, width=8, height=18}
# setup
set.seed(123)
library(lmerTest)

# fit linear mixed model to the ham data:
fm <-
  lmerTest::lmer(
    formula = Informed.liking ~ Gender + Information * Product + (1 | Consumer) +
      (1 | Consumer:Product),
    data = ham
  )

# plot
ggcoefstats(
  x = fm,
  title = "linear mixed-effects models with `lmerTest`"
)
```

## non-linear mixed-effects model (`nlmer`/`nlmerMod`)

```{r nlmer, fig.height=5, fig.width=6}
# data
library(lme4)
set.seed(123)
startvec <- c(Asym = 200, xmid = 725, scal = 350)

# model
nm1 <-
  lme4::nlmer(
    formula = circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym | Tree,
    data = Orange,
    start = startvec
  )

# plot
ggcoefstats(
  x = nm1,
  title = "non-linear mixed-effects model"
)
```

## non-linear least-squares model (`nls`)

```{r nls, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(ggstatsplot)

# model
mod_nls <-
  stats::nls(
    formula = rating ~ k / budget + c,
    data = movies_long,
    start = list(k = 1, c = 0)
  )

# plot
ggcoefstats(
  x = mod_nls,
  title = "non-linear least squares regression",
  subtitle = "Non-linear relationship between budget and rating"
)
```

## conditional generalized linear models for clustered data (`cglm`)

```{r cglm, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(cglm)
data(teenpov)

# model
fit.ide <-
  cglm::cglm(
    method = "ts",
    formula = hours ~ nonpov + inschool + spouse + age + mother,
    data = teenpov,
    id = "ID",
    link = "identity"
  )

# plot
ggcoefstats(
  x = fit.ide,
  title = "conditional generalized linear models for clustered data"
)
```

## joint model to time-to-event data and multivariate longitudinal data (`mjoint`)

```{r mjoint, fig.height=8, fig.width=8}
# setup
set.seed(123)
library(joineRML)
data(heart.valve)

# data
hvd <- heart.valve[!is.na(heart.valve$log.grad) &
  !is.na(heart.valve$log.lvmi) &
  heart.valve$num <= 50, ]

# model
fit_mjoint <- joineRML::mjoint(
  formLongFixed = list(
    "grad" = log.grad ~ time + sex + hs,
    "lvmi" = log.lvmi ~ time + sex
  ),
  formLongRandom = list(
    "grad" = ~ 1 | num,
    "lvmi" = ~ time | num
  ),
  formSurv = Surv(fuyrs, status) ~ age,
  data = hvd,
  inits = list("gamma" = c(0.11, 1.51, 0.80)),
  timeVar = "time"
)

# extract the survival fixed effects and plot them
ggcoefstats(
  x = fit_mjoint,
  component = "conditional",
  package = "yarrr",
  palette = "basel",
  title = "joint model to time-to-event data and multivariate longitudinal data"
)
```

## stationary linear model (`slm`)

```{r slm, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(slm)
data("shan")

# model
mod_slm <-
  slm::slm(
    myformula = shan$PM_Xuhui ~ .,
    data = shan,
    method_cov_st = "fitAR",
    model_selec = -1
  )

# plot
ggcoefstats(
  x = mod_slm,
  conf.level = 0.90,
  title = "stationary linear models",
  package = "rcartocolor",
  palette = "Vivid"
)
```

## generalized linear model (`glm`)

```{r glm1, fig.height=6, fig.width=6}
# setup
library(ggstatsplot)
set.seed(123)

# having a look at the Titanic dataset
df <- as.data.frame(Titanic)

# model
mod_glm <-
  stats::glm(
    formula = Survived ~ Sex + Age,
    data = df,
    weights = df$Freq,
    family = stats::binomial(link = "logit")
  )

# plot
ggcoefstats(
  x = mod_glm,
  ggtheme = ggthemes::theme_economist_white(),
  title = "generalized linear model (glm)",
  vline.args = list(color = "red", linetype = "solid")
)
```

**Note**: The exact statistic will depend on the family used for `glm` models:
Some families will have a `t` statistic associated with them, while others a `z`
statistic. The function will figure this out for you.

```{r glm_full, fig.height=13, fig.width=10}
# creating dataframes to use for regression analyses
set.seed(123)
library(ggstatsplot)

# dataframe #1
df.counts <-
  data.frame(
    treatment = gl(n = 3, k = 3, length = 9),
    outcome = gl(n = 3, k = 1, length = 9),
    counts = c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  ) %>%
  tibble::as_tibble(.)

# dataframe #2
df.clotting <-
  data.frame(
    u = c(5, 10, 15, 20, 30, 40, 60, 80, 100),
    lot1 = c(118, 58, 42, 35, 27, 25, 21, 19, 18),
    lot2 = c(69, 35, 26, 21, 18, 16, 13, 12, 12)
  ) %>%
  tibble::as_tibble(.)

# dataframe #3
x1 <- stats::rnorm(50)
y1 <- stats::rpois(n = 50, lambda = exp(1 + x1))
df.3 <- data.frame(x = x1, y = y1) %>%
  tibble::as_tibble(.)

# dataframe #4
x2 <- stats::rnorm(50)
y2 <- rbinom(
  n = 50,
  size = 1,
  prob = stats::plogis(x2)
)

df.4 <- data.frame(x = x2, y = y2) %>%
  tibble::as_tibble(.)

# combining all plots in a single plot
combine_plots(
  plotlist = list(
    # Family: Poisson
    ggcoefstats(
      x = stats::glm(
        formula = counts ~ outcome + treatment,
        data = df.counts,
        family = stats::poisson(link = "log")
      ),
      title = "Family: Poisson",
      stats.label.color = "black"
    ),
    # Family: Gamma
    ggcoefstats(
      x = stats::glm(
        formula = lot1 ~ log(u),
        data = df.clotting,
        family = stats::Gamma(link = "inverse")
      ),
      title = "Family: Gamma",
      stats.label.color = "black"
    ),
    # Family: Quasi
    ggcoefstats(
      x = stats::glm(
        formula = y ~ x,
        family = quasi(variance = "mu", link = "log"),
        data = df.3
      ),
      title = "Family: Quasi",
      stats.label.color = "black"
    ),
    # Family: Quasibinomial
    ggcoefstats(
      x = stats::glm(
        formula = y ~ x,
        family = stats::quasibinomial(link = "logit"),
        data = df.4
      ),
      title = "Family: Quasibinomial",
      stats.label.color = "black"
    ),
    # Family: Quasipoisson
    ggcoefstats(
      x = stats::glm(
        formula = y ~ x,
        family = stats::quasipoisson(link = "log"),
        data = df.4
      ),
      title = "Family: Quasipoisson",
      stats.label.color = "black"
    ),
    # Family: Gaussian
    ggcoefstats(
      x = stats::glm(
        formula = Sepal.Length ~ Species,
        family = stats::gaussian(link = "identity"),
        data = iris
      ),
      title = "Family: Gaussian",
      stats.label.color = "black"
    )
  ),
  plotgrid.args = list(ncol = 2),
  annotation.args = list(title = "Exploring models with different `glm` families")
)
```

The version of `glm` implemented in `rms` is also supported.

```{r rms_Glm, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(ggstatsplot)
library(rms)

# Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
f_Glm <- rms::Glm(counts ~ outcome + treatment, family = poisson())

# plot
ggcoefstats(
  x = f_Glm,
  title = "rms' implementation of glm"
)
```

## modified fitting for generalized linear models (`glm2`)

```{r glm2, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(glm2)
y <- c(1, 1, 1, 0)

# model
fit_glm2 <-
  glm2::glm2(
    formula = y ~ 1,
    family = binomial(link = "logit"),
    control = glm.control(trace = FALSE)
  )

# plot
ggcoefstats(
  x = fit_glm2,
  title = "greater stability for fitting generalized linear models"
)
```

## Fit Generalized Estimating Equations with `geepack` (`geeglm`)

```{r geeglm, fig.height=7, fig.width=7}
# setup
set.seed(123)
library(geepack)
data(dietox)
dietox$Cu <- as.factor(dietox$Cu)
mf <- formula(Weight ~ Cu * (Time + I(Time^2) + I(Time^3)))

# model
gee1 <-
  geeglm(
    mf,
    data = dietox,
    id = Pig,
    family = poisson("identity"),
    corstr = "ar1"
  )

# plot
ggcoefstats(
  x = gee1,
  title = "Fit Generalized Estimating Equations",
  package = "ggsci",
  palette = "category20c_d3"
)
```

## ordinal regression model (`orm`)

```{r orm, fig.height=5, fig.width=5}
# setup
library(rms)
set.seed(123)

# data
n <- 100
y <- round(runif(n), 2)
x1 <- sample(c(-1, 0, 1), n, TRUE)
x2 <- sample(c(-1, 0, 1), n, TRUE)

# model
g <- rms::orm(y ~ x1 + x2, eps = 1e-5)

# plot
ggcoefstats(
  x = g,
  title = "Ordinal Regression Model"
)
```

## logistic regression model (`lrm`)

```{r lrm, fig.height=12, fig.width=7}
# setup
library(rms)
set.seed(123)

# data
n <- 500
x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
x3 <- sample(0:1, n, TRUE)
y <- x1 + 0.5 * x2 + x3 + rnorm(n)
y <- as.integer(cut2(y, g = 10))
dd <- datadist(x1, x2, x3)
options(datadist = "dd")

# model
f_lrm <- rms::lrm(y ~ x1 + pol(x2, 2) + x3, eps = 1e-7) # eps to check against rstan

# plot
ggcoefstats(
  x = f_lrm,
  title = "Logistic Regression Model",
  package = "ggsci",
  palette = "category20c_d3"
)
```

## Two-Stage Least Squares Instrumental Variables Regression (`iv_robust`)

```{r iv_robust, fig.height=6, fig.width=8}
# setup
set.seed(123)
library(fabricatr)
library(estimatr)

# data
dat <-
  fabricate(
    N = 40,
    Y = rpois(N, lambda = 4),
    Z = rbinom(N, 1, prob = 0.4),
    D = Z * rbinom(N, 1, prob = 0.8),
    X = rnorm(N),
    G = sample(letters[1:4], N, replace = TRUE)
  )

# instrument for treatment `D` with encouragement `Z`
mod_ivrobust <- estimatr::iv_robust(formula = Y ~ D + X | Z + X, data = dat)

# plot
ggcoefstats(
  x = mod_ivrobust,
  title = "Two-Stage Least Squares Instrumental Variables Regression"
)
```

## ordinary least squares with robust standard errors (`lm_robust`)

```{r lm_robust, fig.height=5, fig.width=5}
# for reproducibility
set.seed(123)
library(estimatr)

# model
mod_lmrobust <-
  estimatr::lm_robust(
    formula = mpg ~ gear + wt + cyl,
    data = mtcars
  )

# plot
ggcoefstats(
  x = mod_lmrobust,
  title = "ordinary least squares with robust standard errors"
)
```

## fitting negative binomial GLM (`negbin`)

Just to demonstrate that this can be done, let's also flip the axes:

```{r negbin, fig.height=6, fig.width=8}
# setup
library(MASS)
library(lme4)
set.seed(101)

# data
dd <-
  expand.grid(
    f1 = factor(1:3),
    f2 = LETTERS[1:2],
    g = 1:9,
    rep = 1:15,
    KEEP.OUT.ATTRS = FALSE
  )
mu <- 5 * (-4 + with(dd, as.integer(f1) + 4 * as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

# model
m.glm <- MASS::glm.nb(formula = y ~ f1 * f2, data = dd)

# plot
ggcoefstats(
  x = m.glm,
  title = "generalized linear model (GLM) for the negative binomial family",
  only.significant = TRUE,
  stats.label.args = list(size = 2.5, direction = "both")
) +
  ggplot2::coord_flip()
```

## generalized linear mixed-effects model (`glmer`/`glmerMod`)

```{r glmer, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(lme4)

# model
mod_glmer <-
  lme4::glmer(
    formula = Survived ~ Sex + Age + (Sex + Age | Class),
    data = Titanic_full,
    family = stats::binomial(link = "logit"),
    control = lme4::glmerControl(
      optimizer = "Nelder_Mead",
      calc.derivs = FALSE,
      boundary.tol = 1e-7
    )
  )

# plot
ggcoefstats(
  x = mod_glmer,
  title = "generalized linear mixed-effects model"
)
```

## Fitting Generalized Linear Mixed Models using MCML (`glmm`)

```{r glmm, fig.height=5, fig.width=5}
# setup
library(glmm)
data(BoothHobert)
set.seed(1234)

# model
mod.mcml1 <-
  glmm::glmm(
    fixed = y ~ 0 + x1,
    random = list(y ~ 0 + z1),
    varcomps.names = c("z1"),
    data = BoothHobert,
    family.glmm = bernoulli.glmm,
    m = 100,
    doPQL = TRUE
  )

# plot
ggcoefstats(
  x = mod.mcml1,
  title = "Fitting Generalized Linear Mixed Models using MCML"
)
```

## fitting negative binomial GLMM (`glmer.nb`)

```{r glmer.nb, fig.height=7, fig.width=6}
# setup
library(MASS)
library(lme4)
set.seed(101)

# data
dd <-
  expand.grid(
    f1 = factor(1:3),
    f2 = LETTERS[1:2],
    g = 1:9,
    rep = 1:15,
    KEEP.OUT.ATTRS = FALSE
  )
mu <- 5 * (-4 + with(dd, as.integer(f1) + 4 * as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

# model
m.nb <- lme4::glmer.nb(formula = y ~ f1 * f2 + (1 | g), data = dd)

# plot
ggcoefstats(
  x = m.nb,
  title = "generalized linear mixed-effects model (GLMM) for the negative binomial family"
)
```

## Zero-Inflated Count Data Regression (`zeroinfl`)

```{r zeroinfl, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(pscl)

# data
data("bioChemists", package = "pscl")

# model
mod_zeroinfl <-
  pscl::zeroinfl(
    formula = art ~ . | 1,
    data = bioChemists,
    dist = "negbin"
  )

# plot
ggcoefstats(
  x = mod_zeroinfl,
  title = "Zero-Inflated Count Data Regression"
)
```

## Vector Generalized Additive Models (`vgam`)

```{r vgam, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(VGAM)
pneumo <- transform(pneumo, let = log(exposure.time))

# model
mod_vgam <-
  VGAM::vgam(
    cbind(normal, mild, severe) ~ s(let),
    cumulative(parallel = TRUE),
    data = pneumo,
    trace = FALSE
  )

# plot
ggcoefstats(
  x = mod_vgam,
  title = "Vector Generalized Additive Models"
)
```

## Vector Generalized Linear Models (`vglm`)

```{r vglm, fig.height=4, fig.width=6}
# setup
set.seed(123)
library(VGAM)
pneumo <- transform(pneumo, let = log(exposure.time))

# model
mod_vglm <-
  VGAM::vglm(
    formula = cbind(normal, mild, severe) ~ let,
    family = multinomial,
    data = pneumo
  )

# plot
ggcoefstats(
  x = mod_vglm,
  title = "Vector Generalized Linear Models"
)
```

## Reduced-Rank Vector Generalized Linear Models (`rrvglm`)

```{r rrvglm, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(VGAM)

# data
nn <- 1000 # Number of observations
delta1 <- 3.0 # Specify this
delta2 <- 1.5 # Specify this; should be greater than unity
a21 <- 2 - delta2
mydata <- data.frame(x2 = runif(nn), x3 = runif(nn))
mydata <- transform(mydata, mu = exp(2 + 3 * x2 + 0 * x3))
mydata <- transform(mydata,
  y2 = rnbinom(nn, mu = mu, size = (1 / delta1) * mu^a21)
)

# model
rrnb2 <-
  VGAM::rrvglm(
    formula = y2 ~ x2 + x3,
    family = negbinomial(zero = NULL),
    data = mydata,
    trace = FALSE
  )

# plot
ggcoefstats(
  x = rrnb2,
  title = "Reduced-Rank Vector Generalized Linear Models"
)
```

## Vector autoregression models (`varest`)

```{r varest, fig.height=12, fig.width=7}
# setup
set.seed(123)
library(vars)
data(Canada)

# model
mod_varest <- vars::VAR(Canada, p = 2, type = "none")

# plot
ggcoefstats(
  x = mod_varest,
  title = "Vector autoregression models"
)
```

## Constrained Generalized Additive Model Fitting (`cgam`)

```{r cgam, fig.height=4, fig.width=6}
# setup
set.seed(123)
library(cgam)
data(cubic)

# model
m_cgam <- cgam::cgam(formula = y ~ incr.conv(x), data = cubic)

# plot
ggcoefstats(
  x = m_cgam,
  title = "Constrained Generalized Additive Model Fitting"
)
```

## Constrained Generalized Additive Mixed-Effects Model Fitting (`cgamm`)

```{r cgamm, fig.height=4, fig.width=7}
# setup
set.seed(123)
library(cgam)

# simulate a balanced data set with 30 clusters
# each cluster has 30 data points
n <- 30
m <- 30

# the standard deviation of between cluster error terms is 1
# the standard deviation of within cluster error terms is 2
sige <- 1
siga <- 2

# generate a continuous predictor
x <- 1:(m * n)
for (i in 1:m) {
  x[(n * (i - 1) + 1):(n * i)] <- round(runif(n), 3)
}
# generate a group factor
group <- trunc(0:((m * n) - 1) / n) + 1

# generate the fixed-effect term
mu <- 10 * exp(10 * x - 5) / (1 + exp(10 * x - 5))

# generate the random-intercept term asscosiated with each group
avals <- rnorm(m, 0, siga)

# generate the response
y <- 1:(m * n)
for (i in 1:m) {
  y[group == i] <- mu[group == i] + avals[i] + rnorm(n, 0, sige)
}

# use REML method to fit the model
ans <- cgam::cgamm(formula = y ~ s.incr(x) + (1 | group), reml = TRUE)

# plot
ggcoefstats(
  x = ans,
  title = "Constrained Generalized Additive Mixed-Effects Model Fitting"
)
```

## shape constrained additive models (`scam`)

```{r scam}
# setup
set.seed(123)
library(scam)

# data
n <- 200
x1 <- runif(n) * 6 - 3
f1 <- 3 * exp(-x1^2) # unconstrained term
f1 <- (f1 - min(f1)) / (max(f1) - min(f1)) # function scaled to have range [0,1]
x2 <- runif(n) * 4 - 1
f2 <- exp(4 * x2) / (1 + exp(4 * x2)) # monotone increasing smooth
f2 <- (f2 - min(f2)) / (max(f2) - min(f2)) # function scaled to have range [0,1]
f <- f1 + f2
y <- f + rnorm(n) * 0.1
dat <- data.frame(x1 = x1, x2 = x2, y = y)

# model
b_scam <-
  scam::scam(
    y ~ s(x1, k = 15, bs = "cr", m = 2) + s(x2, k = 25, bs = "mpi", m = 2),
    family = gaussian(link = "identity"),
    data = dat,
    not.exp = FALSE
  )

# plot
ggcoefstats(
  x = b_scam,
  title = "Shape constrained additive models"
)
```

## Hurdle Models for Count Data Regression (`hurdle`)

```{r hurdle, fig.height=10, fig.width=8}
# setup
set.seed(123)
library(pscl)
data("bioChemists", package = "pscl")

# geometric-poisson
fm_hp2 <-
  pscl::hurdle(
    formula = art ~ .,
    data = bioChemists,
    zero = "geometric"
  )

# plot
ggcoefstats(
  x = fm_hp2,
  only.significant = TRUE,
  conf.level = 0.99,
  title = "Hurdle Models for Count Data Regression"
)
```

## beta-binomial mixed-effects model (`BBmm`)

```{r BBmm, fig.height=5, fig.width=6}
# setup
if (isFALSE("PROreg" %in% installed.packages())) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/PROreg/PROreg_1.0.tar.gz",
    repos = NULL,
    type = "source"
  )
}
library(PROreg)
set.seed(123)

# defining the parameters
k <- 100
m <- 10
phi <- 0.5
beta <- c(1.5, -1.1)
sigma <- 0.5

# simulating the covariate and random effects
x <- runif(k, 0, 10)
X <- model.matrix(~x)
z <- as.factor(rBI(k, 4, 0.5, 2))
Z <- model.matrix(~ z - 1)
u <- rnorm(5, 0, sigma)

# the linear predictor and simulated response variable
eta <- beta[1] + beta[2] * x + crossprod(t(Z), u)
p <- 1 / (1 + exp(-eta))
y <- rBB(k, m, p, phi)
dat <- data.frame(cbind(y, x, z))
dat$z <- as.factor(dat$z)

# apply the model
mod_BBmm <-
  PROreg::BBmm(
    fixed.formula = y ~ x,
    random.formula = ~z,
    m = m,
    data = dat
  )

# plot
ggcoefstats(
  x = mod_BBmm,
  title = "beta-binomial mixed-effects model"
)
```

## beta-binomial logistic regression model (`BBreg`)

```{r BBreg, fig.height=5, fig.width=5}
# setup
set.seed(18)
library(PROreg)

# we simulate a covariate, fix the paramters of the beta-binomial
# distribution and simulate a response variable.
# then we apply the model, and try to get the same values.
k <- 1000
m <- 10
x <- rnorm(k, 5, 3)
beta <- c(-10, 2)
p <- 1 / (1 + exp(-1 * (beta[1] + beta[2] * x)))
phi <- 1.2
y <- PROreg::rBB(k, m, p, phi)

# model
mod_BBreg <- PROreg::BBreg(y ~ x, m)

# plot
ggcoefstats(
  x = mod_BBreg,
  title = "beta-binomial logistic regression model"
)
```

## binary choice models with fixed effects (`bife`)

```{r bife, fig.height=10, fig.width=8}
# setup
set.seed(123)
library(bife)

# binary choice models with fixed effects
mod_bife <-
  bife::bife(
    formula = LFP ~ I(AGE^2) + log(INCH) + KID1 + KID2 + KID3 + factor(TIME) | ID,
    data = psid
  )

# plot
ggcoefstats(
  x = mod_bife,
  title = "binary choice models with fixed effects"
)
```

## Dirichlet regression model (`DirichReg`)

```{r DirichReg, fig.height=8, fig.width=6}
# setup
set.seed(123)
library(DirichletReg)

# data
ALake <- ArcticLake
ALake$Y <- DR_data(ALake[, 1:3])

# fit a quadratic Dirichlet regression models ("common")
mod_DirichReg <- DirichletReg::DirichReg(Y ~ depth + I(depth^2), ALake)

# plot
ggcoefstats(
  x = mod_DirichReg,
  title = "Dirichlet Regression"
)
```

## robust generalized linear models (`robmixglm`)

```{r robmixglm, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(robmixglm)
library(MASS)
data(forbes)

# model
forbes.robustmix <- robmixglm(100 * log10(pres) ~ bp, data = forbes)

# plot
ggcoefstats(
  x = forbes.robustmix,
  title = "robust generalized linear models"
)
```

## generalized linear models with extra parameters (`glmx`)

```{r glmx, fig.height=6, fig.width=8}
# setup
library(glmx)
library(MASS)
set.seed(1)
d <- data.frame(x = runif(200, -1, 1))
d$y <- rnbinom(200, mu = exp(0 + 3 * d$x), size = 1)

# model
m_nb1 <-
  glmx::glmx(
    formula = y ~ x,
    data = d,
    family = negative.binomial,
    xlink = "log",
    xstart = 0
  )

ggcoefstats(
  x = m_nb1,
  title = "Generalized Linear Models with Extra Parameters"
)
```

## generalized linear mixed model trees (`glmertree`)

```{r glmertree, fig.height=6, fig.width=12}
# setup
set.seed(123)
library(glmertree)
data("DepressionDemo", package = "glmertree")

# fit normal linear regression LMM tree for continuous outcome
lt <- glmertree::lmertree(
  formula = depression ~ treatment | cluster | age + anxiety + duration,
  data = DepressionDemo
)

# fit logistic regression GLMM tree for binary outcome
gt <- glmertree::glmertree(
  formula = depression_bin ~ treatment | cluster | age + anxiety + duration,
  data = DepressionDemo
)

# plot
combine_plots(
  plotlist = list(
    ggcoefstats(
      x = lt$lmer,
      title = "normal linear regression LMM tree for continuous outcome"
    ),
    ggcoefstats(
      x = lt$lmer,
      title = "logistic regression GLMM tree for binary outcome"
    )
  )
)
```

## generalized linear mixed models using Penalized Quasi-Likelihood (`glmmPQL`)

```{r glmmPQL, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(MASS)
library(nlme)

# model
mod_glmmPQL <-
  MASS::glmmPQL(
    fixed = y ~ trt + I(week > 2),
    random = ~ 1 | ID,
    family = binomial,
    data = bacteria,
    verbose = FALSE
  )

# plot
ggcoefstats(
  x = mod_glmmPQL,
  title = "generalized linear mixed models \nusing Penalized Quasi-Likelihood"
)
```

## generalized linear mixed models using Template Model Builder (`glmmTMB`)

`glmmTMB` package allows for flexibly fitting generalized linear mixed models
(GLMMs) and extensions. Model objects from this package are also supported.

```{r glmmTMB, fig.height=5, fig.width=6}
# set up
library(glmmTMB)
library(lme4)
set.seed(123)

# model
mod_glmmTMB <-
  glmmTMB::glmmTMB(
    formula = Reaction ~ Days + (Days | Subject),
    data = sleepstudy,
    family = glmmTMB::truncated_poisson()
  )

# plotting the model
ggcoefstats(
  x = mod_glmmTMB,
  title = "generalized linear mixed models using Template Model Builder"
)
```

Another example (given the number of terms, let's only display labels for
significant effects):

```{r glmmTMB2, fig.height=14, fig.width=8}
# setup
set.seed(123)
library(glmmTMB)
data(Salamanders)

# model
zipm3 <-
  glmmTMB(count ~ spp * mined + (1 | site),
    zi = ~ spp * mined,
    data = Salamanders,
    family = "poisson"
  )

# plot
ggcoefstats(
  x = zipm3,
  package = "palettesForR",
  palette = "Inkscape",
  only.significant = TRUE
)
```

## generalized linear mixed models using AD Model Builder (`glmmadmb`)

```{r glmmadmb, fig.height=5, fig.width=7}
# setup
if (isFALSE("glmmADMB" %in% installed.packages())) {
  install.packages("glmmADMB",
    repos = c(
      "http://glmmadmb.r-forge.r-project.org/repos",
      getOption("repos")
    ),
    type = "source"
  )
}
library(glmmADMB)

# simulate values
set.seed(101)
d <- data.frame(f = factor(rep(LETTERS[1:10], each = 10)), x = runif(100))
u <- rnorm(10, sd = 2)
d$eta <- with(d, u[f] + 1 + 4 * x)
pz <- 0.3
zi <- rbinom(100, size = 1, prob = pz)
d$y <- ifelse(zi, 0, rpois(100, lambda = exp(d$eta)))

# fit
zipmodel <-
  glmmADMB::glmmadmb(
    formula = y ~ x + (1 | f),
    data = d,
    family = "poisson",
    zeroInflation = TRUE
  )

# plotting the model
ggcoefstats(
  x = zipmodel,
  title = "generalized linear mixed models using AD Model Builder"
)
```

## multilevel model to a list of data frames (`merModList`)

```{r merModList, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(lme4)
library(merTools)

# data
sim_list <-
  replicate(
    n = 10,
    expr = sleepstudy[sample(row.names(sleepstudy), 180), ],
    simplify = FALSE
  )
fml <- "Reaction ~ Days + (Days | Subject)"

# model
mod_lmerModList <- lmerModList(fml, data = sim_list)

# plot
ggcoefstats(
  x = mod_lmerModList,
  title = "a multilevel model to a list of data frames"
)
```

## cumulative link models (`clm`)

```{r clm, fig.height=6, fig.width=8}
# for reproducibility
set.seed(123)
library(ordinal)

# model
mod_clm <- ordinal::clm(formula = rating ~ temp * contact, data = wine)

# plot
ggcoefstats(
  x = mod_clm,
  stats.label.color = "black",
  title = "cumulative link model (clm)",
  subtitle = "(using `ordinal` package)"
) +
  ggplot2::labs(x = "logit regression coefficient", y = NULL)
```

## cumulative link models - older version (`clm2`)

```{r clm2, fig.height=6, fig.width=6}
# for reproducibility
set.seed(123)
library(ordinal)
library(MASS)
data(housing, package = "MASS")

# data
tab26 <- with(soup, table("Product" = PROD, "Response" = SURENESS))
dimnames(tab26)[[2]] <- c("Sure", "Not Sure", "Guess", "Guess", "Not Sure", "Sure")
dat26 <- expand.grid(sureness = as.factor(1:6), prod = c("Ref", "Test"))
dat26$wghts <- c(t(tab26))

# model
mod_clm2 <-
  ordinal::clm2(
    location = sureness ~ prod,
    scale = ~prod,
    data = dat26,
    weights = wghts,
    link = "logistic"
  )

# plot
ggcoefstats(
  x = mod_clm2,
  title = "older version of `clm`"
)
```

## cumulative link mixed models (`clmm`)

```{r clmm1, fig.height=6, fig.width=6}
# for reproducibility
set.seed(123)
library(ordinal)

# model
mod_clmm <- ordinal::clmm(
  formula = rating ~ temp + contact + (1 | judge),
  data = wine
)

# to speed up calculations, we will use just 10% of the dataset
ggcoefstats(
  x = mod_clmm,
  title = "cumulative link mixed model (clmm)",
  subtitle = "(using `ordinal` package)"
) +
  ggplot2::labs(
    x = "coefficient from ordinal mixed-effects regression",
    y = "fixed effects"
  )
```

## cumulative link mixed models - older version (`clmm2`)

```{r clmm2, fig.height=6, fig.width=6}
# for reproducibility
set.seed(123)
library(ordinal)

# data
dat <- subset(soup, as.numeric(as.character(RESP)) <= 24)
dat$RESP <- dat$RESP[drop = TRUE]

# model
mod_clmm2 <-
  ordinal::clmm2(
    SURENESS ~ PROD,
    random = RESP,
    data = dat,
    link = "probit",
    Hess = TRUE,
    method = "ucminf",
    threshold = "symmetric"
  )

# plot
ggcoefstats(
  x = mod_clmm2,
  title = "older version of cumulative link mixed models"
)
```

## marginal effects estimation (`margins`)

```{r margins, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(margins)

# logit model
mod_log <- glm(
  formula = am ~ cyl + hp + wt,
  data = mtcars,
  family = binomial
)

# convert to marginal effects with margins::margins()
marg_log <- margins(mod_log)

# plot
ggcoefstats(
  x = marg_log,
  title = "marginal effects estimation"
)
```

## Linear Regression with Interval-Censored Dependent Variable (`semLm`)

```{r semLm, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(smicd)

# Load and prepare data
data <- Exam
classes <- c(1, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.7, 8.5, Inf)
data$examsc.class <- cut(data$examsc, classes)

# run model with random intercept and default settings
mod_semLm <-
  smicd::semLm(
    formula = examsc.class ~ standLRT + schavg,
    data = data,
    classes = classes
  )

# plot
ggcoefstats(
  x = mod_semLm,
  title = "Linear Regression with \nInterval-Censored Dependent Variable"
)
```

## Linear Mixed Regression with Interval-Censored Dependent Variable (`semLme`)

```{r semLme, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(smicd)

# Load and prepare data
data <- Exam
classes <- c(1, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.7, 8.5, Inf)
data$examsc.class <- cut(data$examsc, classes)

# run model with random intercept and default settings
model1 <-
  smicd::semLme(
    formula = examsc.class ~ standLRT + schavg + (1 | school),
    data = data,
    classes = classes
  )

# plot
ggcoefstats(
  x = model1,
  title = "Linear Mixed Regression with \nInterval-Censored Dependent Variable"
)
```

<!-- mixor was removed from CRAN -->

<!-- ## Mixed-Effects Ordinal Regression Analysis (`mixor`) -->

<!-- ```{r mixor, fig.height=8, fig.width=6} -->
<!-- # setup -->
<!-- set.seed(123) -->
<!-- library(mixor) -->
<!-- data("SmokingPrevention") -->

<!-- # data frame must be sorted by id variable -->
<!-- SmokingPrevention <- SmokingPrevention[order(SmokingPrevention$class), ] -->

<!-- # school model -->
<!-- mod_mixor <- -->
<!--   mixor::mixor( -->
<!--     formula = thksord ~ thkspre + cc + tv + cctv, -->
<!--     data = SmokingPrevention, -->
<!--     id = school, -->
<!--     link = "logit" -->
<!--   ) -->

<!-- # plot -->
<!-- ggcoefstats( -->
<!--   x = mod_mixor, -->
<!--   title = "Mixed-Effects Ordinal Regression Analysis" -->
<!-- ) -->
<!-- ``` -->

## bias reduction in Binomial-response GLMs (`brglm`)

```{r brglm, fig.height=6, fig.width=5}
# setup
set.seed(123)
library(brglm)
data("lizards")

# fit the model using maximum likelihood mean bias-reduced fit
lizards.brglm <-
  brglm::brglm(
    cbind(grahami, opalinus) ~ height + diameter + light + time,
    family = binomial(logit),
    data = lizards,
    method = "brglm.fit"
  )

# plot
ggcoefstats(
  x = lizards.brglm,
  only.significant = TRUE,
  title = "bias reduction in Binomial-response GLMs"
)
```

## bias reduction in generalized linear models (`brglm2`)

```{r brglm2, fig.height=6, fig.width=5}
# setup
set.seed(123)
library(brglm2)
data("lizards")

# fit the model using maximum likelihood mean bias-reduced fit:
lizardsBR_mean <-
  stats::glm(
    formula = cbind(grahami, opalinus) ~ height + diameter + light + time,
    family = binomial(logit),
    data = lizards,
    method = "brglmFit"
  )

# plot
ggcoefstats(
  x = lizardsBR_mean,
  only.significant = TRUE,
  title = "bias reduction in generalized linear models"
)
```

## Bias Reduction For Multinomial Response Models Using The Poisson Trick (`brmultinom`)

```{r brmultinom, fig.height=12, fig.width=7}
# setup
set.seed(123)
library(MASS)
library(brglm2)
data("housing", package = "MASS")

# Maximum likelihood using brmultinom with baseline category 'Low'
houseML1 <-
  brglm2::brmultinom(
    formula = Sat ~ Infl + Type + Cont,
    weights = Freq,
    data = housing,
    type = "ML",
    ref = 1
  )

# plot
ggcoefstats(
  x = houseML1,
  title = "Bias Reduction For Multinomial Response Models Using The Poisson Trick"
)
```

## bias reduction for adjacent category logit models (`bracl`)

```{r bracl, fig.height=6, fig.width=7}
# setup
set.seed(123)
library(brglm2)
data("stemcell")

# bias reduction for adjacent category logit models
# for ordinal responses using the Poisson trick
fit_bracl <-
  brglm2::bracl(
    formula = research ~ as.numeric(religion) + gender,
    weights = frequency,
    data = stemcell,
    type = "ML"
  )

# plot
ggcoefstats(
  x = fit_bracl,
  title = "bias reduction for adjacent category logit models"
)
```

## generalized linear models subject to population constraints

```{r glmc, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(glmc)

# data
n <- rbind(c(5903, 230), c(5157, 350))
mat <- matrix(0, nrow = sum(n), ncol = 2)
mat <-
  rbind(
    matrix(1, nrow = n[1, 1], ncol = 1) %*% c(0, 0),
    matrix(1, nrow = n[1, 2], ncol = 1) %*% c(0, 1),
    matrix(1, nrow = n[2, 1], ncol = 1) %*% c(1, 0),
    matrix(1, nrow = n[2, 2], ncol = 1) %*% c(1, 1)
  )

# specifying the population constraints
gfr <- .06179 * matrix(1, nrow = nrow(mat), ncol = 1)
g <- matrix(1, nrow = nrow(mat), ncol = 1)
amat <- matrix(mat[, 2] * g - gfr, ncol = 1)

# defining constraints in the data frame.
hrh <- data.frame(birth = mat[, 2], child = mat[, 1], constraints = amat)

# model
gfit <-
  glmc::glmc(
    formula = birth ~ child,
    data = hrh,
    family = "binomial",
    emplik.method = "Owen",
    control = glmc::glmc.control(
      trace.optim = 0,
      trace.glm = FALSE,
      maxit.glm = 10,
      maxit.weights = 200,
      itertrace.weights = FALSE
    )
  )

# plot
ggcoefstats(
  x = gfit,
  title = "generalized linear models subject to population constraints"
)
```

## Bayesian linear mixed-effects models (`blmerMod`)

```{r blmerMod, fig.height=5, fig.width=5}
# for reproducibility
set.seed(123)
library(blme)

# data
data(sleepstudy)
sleepstudy$mygrp <- sample(1:5, size = 180, replace = TRUE)
sleepstudy$mysubgrp <- NA
for (i in 1:5) {
  filter_group <- sleepstudy$mygrp == i
  sleepstudy$mysubgrp[filter_group] <-
    sample(1:30, size = sum(filter_group), replace = TRUE)
}

# model
mod_blmer <-
  blme::blmer(
    formula = scale(Reaction) ~ scale(Days) + (1 + Days | Subject),
    data = sleepstudy,
    cov.prior = NULL,
    REML = FALSE
  )

# plot
ggcoefstats(
  x = mod_blmer,
  title = "Bayesian linear mixed-effects models"
)
```

## Bayesian generalized linear mixed-effects models (`bglmerMod`)

```{r bglmerMod, fig.height=6, fig.width=6}
# for reproducibility
set.seed(123)
library(blme)

# model
mod_bglmer <-
  blme::bglmer(
    formula = Reaction ~ Days + (1 + Days | Subject),
    data = sleepstudy,
    cov.prior = NULL,
    fixef.prior = normal
  )

# plot
ggcoefstats(
  x = mod_bglmer,
  title = "Bayesian generalized linear mixed-effects models"
)
```

## ordered logistic or probit regression (`polr`)

```{r polr, fig.height=6, fig.width=8}
# polr model
set.seed(123)
library(MASS)
polr.mod <-
  MASS::polr(
    formula = Sat ~ Infl + Type + Cont,
    weights = Freq,
    data = housing
  )

# plot
ggcoefstats(
  x = polr.mod,
  coefficient.type = "both",
  title = "ordered logistic or probit regression",
  subtitle = "using `MASS` package"
)
```

## multiple linear regression models (`mlm`)

```{r mlm, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(effectsize)

# model (converting all numeric columns in data to z-scores)
mod_mlm <-
  stats::lm(
    formula = cbind(mpg, disp) ~ wt,
    data = effectsize::standardize(mtcars)
  )

# plot
ggcoefstats(
  x = mod_mlm,
  title = "multiple linear regression models"
)
```

## anova on multiple linear regression models (`maov`)

```{r maov, fig.height=8, fig.width=6}
# setup
set.seed(123)

# model
fit <- lm(cbind(mpg, disp, hp) ~ factor(cyl), data = mtcars)
m_maov <- aov(fit)

# plot
ggcoefstats(
  x = m_maov,
  title = "anova on multiple linear regression models",
  package = "ggsci",
  palette = "springfield_simpsons",
  conf.level = 0.90
)
```

## multinomial logistic regression models (`multinom`)

```{r multinom, fig.height=8, fig.width=6}
# setup
set.seed(123)
library(nnet)
library(MASS)
utils::example(topic = birthwt, echo = FALSE)

# model
bwt.mu <-
  nnet::multinom(
    formula = low ~ .,
    data = bwt,
    trace = FALSE
  )

# plot
ggcoefstats(
  x = bwt.mu,
  title = "multinomial logistic regression models",
  package = "ggsci",
  palette = "default_ucscgb"
)
```

## multilevel mediation model (`bmlm`)

```{r bmlm, fig.height=7, fig.width=5}
# setup
set.seed(123)
library(bmlm)

# model
fit_bmlm <- bmlm::mlm(BLch9, verbose = FALSE)

# exctrating summary
df_summary <-
  bmlm::mlm_summary(fit_bmlm) %>%
  dplyr::rename(
    estimate = Mean,
    term = Parameter,
    std.error = SE,
    conf.low = `2.5%`,
    conf.high = `97.5%`
  )

# plot
ggcoefstats(
  x = df_summary,
  title = "Bayesian multilevel mediation models with Stan"
)
```

## proportional odds and related models (`svyolr`)

```{r svyolr, fig.height=8, fig.width=6}
# setup
set.seed(123)
library(survey)
data(api)

# preparing data
dclus1 <-
  survey::svydesign(
    id = ~dnum,
    weights = ~pw,
    data = apiclus1,
    fpc = ~fpc
  )

# update
dclus1 <- update(dclus1, mealcat = cut(meals, c(0, 25, 50, 75, 100)))

# model
m_svyolr <-
  survey::svyolr(
    formula = mealcat ~ avg.ed + mobility + stype,
    design = dclus1
  )

# plot
ggcoefstats(
  x = m_svyolr,
  title = "proportional odds and related models",
  coefficient.type = "both"
)
```

## survey-weighted generalized linear models (`svyglm`)

```{r svyglm, fig.height=5, fig.width=5}
# data
library(survey)
set.seed(123)
data(api)
dstrat <-
  survey::svydesign(
    id = ~1,
    strata = ~stype,
    weights = ~pw,
    data = apistrat,
    fpc = ~fpc
  )

# model
mod_svyglm <-
  survey::svyglm(
    formula = sch.wide ~ ell + meals + mobility,
    design = dstrat,
    family = quasibinomial()
  )

# plot
ggcoefstats(
  x = mod_svyglm,
  title = "survey-weighted generalized linear model"
)
```

## design-based inference for vector generalised linear models (`svy_vglm`)

```{r , fig.height=7, fig.width=5}
# setup
set.seed(123)
library(svyVGAM)
data(api)

# data
dclus2 <- svydesign(
  id = ~ dnum + snum,
  fpc = ~ fpc1 + fpc2,
  data = apiclus2
)

# model
mod_svy_vglm <- svy_vglm(api00 ~ api99 + mobility + ell,
  design = dclus2,
  family = uninormal()
)

# plot
ggcoefstats(
  x = mod_svy_vglm,
  title = "design-based inference for\nvector generalised linear models"
)
```

## estimation of limited dependent variable models (`mhurdle`)

```{r mhurdle, fig.height=6, fig.width=8}
# setup
data("Interview", package = "mhurdle")
library(mhurdle)

# independent double hurdle model
idhm <- mhurdle::mhurdle(
  formula = vacations ~ car + size | linc + linc2 | 0,
  data = Interview,
  dist = "ln",
  h2 = TRUE,
  method = "bfgs"
)

# plot
ggcoefstats(
  x = idhm,
  title = "estimation of limited dependent variable models"
)
```

## repeated measures ANOVA (`aovlist`)

```{r aovlist1, fig.height=6, fig.width=8}
# for reproducibility
set.seed(123)
library(ggstatsplot)

# specifying the model (note the error structure)
mod_aovlist <-
  stats::aov(
    formula = value ~ attribute * measure + Error(id / (attribute * measure)),
    data = iris_long
  )

# plot
ggcoefstats(
  x = mod_aovlist,
  effsize = "eta",
  ggtheme = ggthemes::theme_fivethirtyeight(),
  title = "Variation in measurements for Iris species",
  subtitle = "Source: Iris data set (by Fisher or Anderson)",
  caption = "Results from 2 by 2 RM ANOVA"
) +
  ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 11, face = "plain"))
```

## robust regression with `robust` package (`lmRob`, `glmRob`)

```{r robust, fig.height=10, fig.width=6}
combine_plots(
  plotlist = list(
    # plot 1: glmRob
    ggcoefstats(
      x = robust::glmRob(
        formula = Survived ~ Sex,
        data = Titanic_full,
        family = stats::binomial(link = "logit")
      ),
      title = "generalized robust linear model",
      ggtheme = ggthemes::theme_fivethirtyeight()
    ),
    # plot 2: lmRob
    ggcoefstats(
      x = robust::lmRob(
        formula = Sepal.Length ~ Sepal.Width * Species,
        data = iris
      ),
      title = "robust linear model",
      package = "awtools",
      palette = "a_palette",
      ggtheme = ggthemes::theme_tufte()
    )
  ),
  # arguments relevant for `combine_plots` function
  annotation.args = list(title = "Robust variants of `lmRob` and `glmRob` \n(from`robust` package)"),
  plotgrid.args = list(nrow = 2)
)
```

## robust regression with `robustbase` package (`lmrob`, `glmrob`)

Another alternative is to use robust models, as implemented in the `robustbase`
package.

```{r robustbase, fig.height=10, fig.width=6}
# setup
set.seed(123)
library(robustbase)

# dataframe
data(coleman)
clotting <-
  data.frame(
    u = c(5, 10, 15, 20, 30, 40, 60, 80, 100),
    lot1 = c(118, 58, 42, 35, 27, 25, 21, 19, 18),
    lot2 = c(69, 35, 26, 21, 18, 16, 13, 12, 12)
  )

# combined plot for both generalized and simple robust models
combine_plots(
  plotlist = list(
    # plot 1: glmrob
    ggcoefstats(
      x = glmrob(
        formula = lot1 ~ log(u),
        data = clotting,
        family = Gamma
      ),
      title = "generalized robust linear model"
    ),
    # plot 2: lmrob
    ggcoefstats(
      x = lmrob(formula = Y ~ ., data = coleman),
      title = "robust linear model"
    )
  ),
  # arguments relevant for `combine_plots`
  annotation.args = list(title = "Robust variants of `lmRob` and `glmRob` \n(from`robustbase` package)"),
  plotgrid.args = list(nrow = 2)
)
```

## MM-type estimators for linear regression on compositional data (`complmrob`)

```{r complmrob, fig.height=7, fig.width=6}
# setup
set.seed(123)
library(complmrob)

# data
crimes <- data.frame(
  lifeExp = state.x77[, "Life Exp"],
  USArrests[, c("Murder", "Assault", "Rape")]
)

# model
mUSArr <- complmrob::complmrob(formula = lifeExp ~ ., data = crimes)

# plot
ggcoefstats(
  x = mUSArr,
  title = "MM-type estimators for linear regression on compositional data"
)
```

## fit a nonlinear heteroscedastic model via maximum likelihood (`nlreg`)

```{r nlreg, fig.height=5, fig.width=7}
set.seed(123)
library(nlreg)
library(boot)
data(calcium)

# homoscedastic model fit
calcium.nl <-
  nlreg::nlreg(
    formula = cal ~ b0 * (1 - exp(-b1 * time)),
    start = c(b0 = 4, b1 = 0.1),
    data = calcium
  )

# plot
ggcoefstats(
  x = calcium.nl,
  conf.int = FALSE,
  title = "fit a nonlinear heteroscedastic model via maximum likelihood"
)
```

## fit a linear model with multiple group fixed effects (`felm`)

```{r felm, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(lfe)

# create covariates
x <- rnorm(1000)
x2 <- rnorm(length(x))

# individual and firm
id <- factor(sample(20, length(x), replace = TRUE))
firm <- factor(sample(13, length(x), replace = TRUE))

# effects for them
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

# left hand side
u <- rnorm(length(x))
y <- x + 0.5 * x2 + id.eff[id] + firm.eff[firm] + u

# estimate and print result
est <- lfe::felm(formula = y ~ x + x2 | id + firm)

# plot
ggcoefstats(
  x = est,
  title = "linear model with multiple group fixed effects"
)
```

## linear models for panel data (`plm`)

```{r plm, fig.height=5, fig.width=5}
# data
set.seed(123)
library(plm)
data("Produc", package = "plm")

# model
plm.mod <-
  plm::plm(
    formula = log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
    data = Produc,
    index = c("state", "year")
  )

# plot
ggcoefstats(
  x = plm.mod,
  title = "linear models for panel data"
)
```

## multinomial logit model (`mlogit`)

```{r mlogit, fig.height=8, fig.width=7}
# setup
set.seed(123)
library(mlogit)

# data
data("Fishing", package = "mlogit")
Fish <-
  mlogit::mlogit.data(Fishing,
    varying = c(2:9),
    shape = "wide",
    choice = "mode"
  )

# a "mixed" model
m_mlogit <- mlogit::mlogit(mode ~ price + catch | income, data = Fish)

# plot
ggcoefstats(
  x = m_mlogit,
  title = "multinomial logit model"
)
```

## Cox proportional hazards regression model (`coxph`)

```{r coxph, fig.height=4, fig.width=4}
# for reproducibility
set.seed(123)
library(survival)

# create the simplest-test data set
test1 <- list(
  time = c(4, 3, 1, 1, 2, 2, 3),
  status = c(1, 1, 1, 0, 1, 1, 0),
  x = c(0, 2, 1, 1, 1, 0, 0),
  sex = c(0, 0, 0, 0, 1, 1, 1)
)

# fit a stratified model
mod_coxph <-
  survival::coxph(
    formula = Surv(time, status) ~ x + strata(sex),
    data = test1
  )

# plot
ggcoefstats(
  x = mod_coxph,
  title = "Cox proportional hazards regression model"
)
```

Another example with `frailty` term.

```{r coxph.penal, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(survival)

# model
mod_coxph <- survival::coxph(
  formula = Surv(time, status) ~ age + sex + frailty(inst),
  data = lung
)

# plot
ggcoefstats(
  x = mod_coxph,
  title = "Proportional Hazards Regression Model\nwith Frailty penalty function"
)
```

## mixed effects Cox model (`coxme`)

```{r coxme, fig.height=8, fig.width=6}
# setup
set.seed(123)
library(survival)
library(coxme)

# model
fit <- coxme::coxme(
  formula = Surv(y, uncens) ~ trt + (1 | center),
  data = eortc
)

# plot
ggcoefstats(
  x = fit,
  title = "mixed effects Cox model"
)
```

## robust Cox proportional hazards regression model (`coxr`)

```{r coxr, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(coxrobust)

# create a simple test data set using the attached function `gen_data`
a <- coxrobust::gen_data(200, c(1, 0.1, 2), cont = 0.05, p.censor = 0.30)

# model
mod_coxr <- coxrobust::coxr(Surv(time, status) ~ X1 + X2 + X3, data = a, trunc = 0.9)

# plot
ggcoefstats(
  x = mod_coxr,
  title = "robust Cox proportional hazards regression model"
)
```

## truncated Gaussian Regression Models (`truncreg`)

```{r truncreg, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(truncreg)
library(survival)

# data
data("tobin", package = "survival")

# model
cragg_trunc <-
  truncreg::truncreg(
    formula = durable ~ age + quant,
    data = tobin,
    subset = durable > 0
  )

# plot
ggcoefstats(
  x = cragg_trunc,
  title = "Truncated Gaussian Regression Models"
)
```

## Fitting Linear Quantile Models

```{r lqm, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(lqmm)

# data
n <- 500
p <- 1:3 / 4
test <- data.frame(x = runif(n, 0, 1))
test$y <- 30 + test$x + rnorm(n)

# model
fit.lqm <-
  lqmm::lqm(
    y ~ x,
    data = test,
    tau = p,
    control = list(verbose = FALSE, loop_tol_ll = 1e-9),
    fit = TRUE
  )

# plot
ggcoefstats(
  x = fit.lqm,
  title = "Fitting Linear Quantile Models"
)
```

## Fitting Linear Quantile Mixed Models

```{r lqmm, fig.height=4, fig.width=5}
# setup
set.seed(123)
library(lqmm)

# data
M <- 50
n <- 10
test <- data.frame(x = runif(n * M, 0, 1), group = rep(1:M, each = n))
test$y <- 10 * test$x + rep(rnorm(M, 0, 2), each = n) + rchisq(n * M, 3)

# model
fit.lqmm <-
  lqmm::lqmm(
    fixed = y ~ x,
    random = ~1,
    group = group,
    data = test,
    tau = 0.5,
    nK = 11,
    type = "normal"
  )

# plot
ggcoefstats(
  x = fit.lqmm,
  title = "Fitting Linear Quantile Mixed Models"
)
```

## autoregressive integrated moving average (`Arima`)

```{r arima, fig.height=5, fig.width=5}
# for reproducibility
set.seed(123)

# model
fit <- stats::arima(x = lh, order = c(1, 0, 0))

# plot
ggcoefstats(
  x = fit,
  title = "autoregressive integrated moving average"
)
```

## high performance linear model (`speedlm`/`speedglm`)

Example of high performance linear model-

```{r speedlm, fig.height=5, fig.width=5}
# setup
library(speedglm)
set.seed(123)

# model
mod_speedlm <-
  speedglm::speedlm(
    formula = mpg ~ wt + qsec,
    data = mtcars,
    fitted = TRUE
  )

# plot
ggcoefstats(
  x = mod_speedlm,
  title = "high performance linear model"
)
```

Example of high performance generalized linear model-

```{r speedglm, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(speedglm)

# data
n <- 50000
k <- 5
y <- rgamma(n, 1.5, 1)
x <- round(matrix(rnorm(n * k), n, k), digits = 3)
colnames(x) <- paste("s", 1:k, sep = "")
da <- data.frame(y, x)
fo <- as.formula(paste("y~", paste(paste("s", 1:k, sep = ""), collapse = "+")))

# model
mod_speedglm <-
  speedglm::speedglm(
    formula = fo,
    data = da,
    family = stats::Gamma(log)
  )

# plot
ggcoefstats(
  x = mod_speedglm,
  title = "high performance generalized linear model"
)
```

## parametric survival regression model (`survreg`)

```{r survreg, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(survival)

# model
mod_survreg <-
  survival::survreg(
    formula = Surv(futime, fustat) ~ ecog.ps + rx,
    data = ovarian,
    dist = "logistic"
  )

# plot
ggcoefstats(
  x = mod_survreg,
  ggtheme = hrbrthemes::theme_ipsum_rc(),
  package = "ggsci",
  palette = "legacy_tron",
  title = "parametric survival regression model"
)
```

## Competing Risks Regression (`crr`)

```{r crr, fig.height=6, fig.width=6}
# setup
set.seed(10)
library(cmprsk)

# simulated data to test
ftime <- rexp(200)
fstatus <- sample(0:2, 200, replace = TRUE)
cov <- matrix(runif(600), nrow = 200)
dimnames(cov)[[2]] <- c("x1", "x2", "x3")

# model
mod_crr <- cmprsk::crr(ftime, fstatus, cov)

# plot
ggcoefstats(
  x = mod_crr,
  title = "Competing Risks Regression"
)
```

## tobit regression (`tobit`)

```{r tobit, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(AER)
data("Affairs", package = "AER")

# model
m_tobit <-
  AER::tobit(
    formula = affairs ~ age + yearsmarried + religiousness + occupation + rating,
    data = Affairs
  )

# plot
ggcoefstats(
  x = m_tobit,
  title = "tobit regression"
)
```

## censored regression (tobit) regression (`censReg`)

```{r censReg, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(censReg)
data("Affairs", package = "AER")

# model
estResult <-
  censReg::censReg(
    formula = affairs ~ age + yearsmarried + religiousness + occupation + rating,
    data = Affairs
  )

# plot
ggcoefstats(
  x = estResult,
  title = "censored regression (tobit) regression"
)
```

## relative risk regression model for case-cohort studies (`cch`)

```{r cch, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(survival)

# examples come from cch documentation
subcoh <- nwtco$in.subcohort
selccoh <- with(nwtco, rel == 1 | subcoh == 1)
ccoh.data <- nwtco[selccoh, ]
ccoh.data$subcohort <- subcoh[selccoh]

## central-lab histology
ccoh.data$histol <- factor(ccoh.data$histol, labels = c("FH", "UH"))

## tumour stage
ccoh.data$stage <- factor(ccoh.data$stage, labels = c("I", "II", "III", "IV"))
ccoh.data$age <- ccoh.data$age / 12 # Age in years

# model
fit.ccP <-
  survival::cch(
    formula = Surv(edrel, rel) ~ stage + histol + age,
    data = ccoh.data,
    subcoh = ~subcohort,
    id = ~seqno,
    cohort.size = 4028
  )

# plot
ggcoefstats(
  x = fit.ccP,
  title = "relative risk regression model",
  subtitle = "(for case-cohort studies)",
  conf.level = 0.99
)
```

## ridge regression (`ridgelm`)

For ridge regression, neither statistic values nor confidence intervals for
estimates are available, so only estimates will be displayed.

```{r ridgelm, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(MASS)

# model
names(longley)[1] <- "y"
mod_ridgelm <- MASS::lm.ridge(formula = y ~ ., data = longley)

# plot
ggcoefstats(
  x = mod_ridgelm,
  title = "ridge regression"
)
```

## generalized additive models with integrated smoothness estimation (`gam`)

**Important**: These model outputs contains both parametric and smooth terms.
`ggcoefstats` only displays the parametric terms.

```{r gam_mgcv, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(mgcv)

# model
g_gam <-
  mgcv::gam(
    formula = mpg ~ s(hp) + am + qsec,
    family = stats::quasi(),
    data = mtcars
  )

# plot
ggcoefstats(
  x = g_gam,
  title = "generalized additive models with \nintegrated smoothness estimation",
  subtitle = "using `mgcv` package"
)
```

## 

```{r bam_mgcv, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(mgcv)

# data
dat <- gamSim(1, n = 25000, dist = "normal", scale = 20)
bs <- "cr"
k <- 12

# model
b_bam <-
  mgcv::bam(
    formula = y ~ s(x0, bs = bs) + s(x1, bs = bs) + s(x2, bs = bs, k = k) +
      s(x3, bs = bs),
    data = dat
  )

# plot
ggcoefstats(
  x = b_bam,
  title = "generalized additive models for \nvery large datasets"
)
```

## generalized additive model (`Gam`)

```{r Gam, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(gam)

# model
mod_gam <- gam::gam(
  formula = mpg ~ s(hp, 4) + am + qsec,
  data = mtcars
)

# plot
ggcoefstats(
  x = mod_gam,
  title = "generalized additive model",
  subtite = "(using `gam` package)"
)
```

## linear mixed-effects models (`lme`)

```{r lme, fig.height=5, fig.width=5}
# for reproducibility
set.seed(123)
library(lme4)
library(nlme)
data("sleepstudy")

# model
mod_lme <-
  nlme::lme(
    fixed = Reaction ~ Days,
    random = ~ 1 + Days | Subject,
    data = sleepstudy
  )

# plot
ggcoefstats(
  x = mod_lme,
  title = "linear mixed-effects models (`lme`)"
)
```

## linear model using generalized least squares (`gls`)

The `nlme` package provides a function to fit a linear model using generalized
least squares. The errors are allowed to be correlated and/or have unequal
variances.

```{r gls, fig.height=5, fig.width=6}
# for reproducibility
set.seed(123)
library(nlme)

# model
mod_gls <-
  nlme::gls(
    model = follicles ~ sin(2 * pi * Time) + cos(2 * pi * Time),
    data = Ovary,
    correlation = corAR1(form = ~ 1 | Mare)
  )

# plot
ggcoefstats(
  x = mod_gls,
  stats.label.color = "black",
  ggtheme = hrbrthemes::theme_ipsum_ps(),
  title = "generalized least squares model"
)
```

## inference for estimated coefficients (`coeftest`)

```{r coeftest, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(lmtest)

# load data and fit model
data("Mandible", package = "lmtest")
fm <- stats::lm(formula = length ~ age, data = Mandible, subset = (age <= 28))

# the following commands lead to the same tests
ct <- lmtest::coeftest(fm)

# plot
ggcoefstats(
  x = ct,
  plot = "inference for estimated coefficients",
  conf.level = 0.99
)
```

## robust regression using an M estimator (`rlm`)

```{r rlm, fig.height=5, fig.width=5}
# for reproducibility
set.seed(123)

# model
mod_rlm <- MASS::rlm(formula = mpg ~ am * cyl, data = mtcars)

# plot
ggcoefstats(
  x = mod_rlm,
  point.args = list(color = "red", shape = 15),
  vline.args = list(size = 1, color = "#CC79A7", linetype = "dotdash"),
  title = "robust regression using an M estimator",
  ggtheme = ggthemes::theme_stata()
)
```

## quantile regression (`rq`)

```{r rq, fig.height=5, fig.width=6}
# for reproducibility
set.seed(123)
library(quantreg)

# data
data(stackloss)

# model
mod_rq <-
  quantreg::rq(
    formula = stack.loss ~ .,
    data = stackloss,
    method = "br"
  )

# plot
ggcoefstats(
  x = mod_rq,
  se.type = "iid",
  title = "quantile regression"
)
```

## multiple response quantile regression (`rqs`)

```{r rqs, fig.height=7, fig.width=6}
# setup
set.seed(123)
library(quantreg)
data("engel")

# model
fm <-
  quantreg::rq(
    formula = foodexp ~ income,
    data = engel,
    tau = 1:4 / 10
  )

# plot
ggcoefstats(
  x = fm,
  title = "multiple response quantile regression (`rqs`)"
)
```

## nonlinear quantile regression estimates (`nlrq`)

```{r nlrq, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(quantreg)

# preparing data
Dat <- NULL
Dat$x <- rep(1:25, 20)
Dat$y <- stats::SSlogis(Dat$x, 10, 12, 2) * rnorm(500, 1, 0.1)

# then fit the median using nlrq
Dat.nlrq <-
  quantreg::nlrq(
    formula = y ~ SSlogis(x, Asym, mid, scal),
    data = Dat,
    tau = 0.5,
    trace = FALSE
  )

# plot
ggcoefstats(
  x = Dat.nlrq,
  title = "non-linear quantile regression",
  se.type = "nid"
)
```

## additive quantile regression (`rqss`)

```{r rqss, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(quantreg)
data(CobarOre)

# model
fCO <- quantreg::rqss(
  formula = z ~ qss(cbind(x, y), lambda = .08),
  data = CobarOre
)

# model info
ggcoefstats(
  x = fCO,
  title = "Additive Quantile Regression"
)
```

## censored quantile regression (`crq`)

```{r crq, fig.height=5, fig.width=5}
# crq example with left censoring
set.seed(1968)
library(quantreg)

# data
n <- 200
x <- rnorm(n)
y <- 5 + x + rnorm(n)
c <- 4 + x + rnorm(n)
d <- (y > c)

# model
f_crq <- quantreg::crq(survival::Surv(pmax(y, c), d, type = "left") ~ x,
  method = "Portnoy"
)

# plot
ggcoefstats(
  x = f_crq,
  title = "censored quantile regression models"
)
```

## instrumental-variable regression (`ivreg`)

```{r ivreg, fig.height=7, fig.width=5}
# setup
suppressPackageStartupMessages(library(AER))
set.seed(123)
data("CigarettesSW", package = "AER")

# model
ivr <-
  AER::ivreg(
    formula = log(packs) ~ income | population,
    data = CigarettesSW,
    subset = year == "1995"
  )

# plot
ggcoefstats(
  x = ivr,
  title = "instrumental-variable regression"
)
```

## instrumental variables probit regression (`ivprobit`)

```{r ivprobit, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(ivprobit)
data(eco)

# model
pro <-
  ivprobit::ivprobit(
    formula = d2 ~ ltass + roe + div | eqrat + bonus | ltass + roe + div + gap + cfa,
    data = eco
  )

# plot
ggcoefstats(
  x = pro,
  title = "instrumental variables probit regression"
)
```

## instrumental fixed effect panel regression (`ivFixed`)

```{r ivFixed, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(ivfixed)

#  data
pib <- as.matrix(c(12, 3, 4, 0.4, 0.7, 5, 0.7, 0.3, 0.6, 89, 7, 8, 45, 7, 4, 5, 0.5, 5),
  nrows = 18,
  ncols = 1
)
tir <- as.matrix(c(12, 0.3, 4, 0.4, 7, 12, 3.0, 6.0, 45, 7.0, 0.8, 44, 65, 23, 4, 6, 76, 9),
  nrows = 18,
  ncols = 1
)
inf <- as.matrix(c(1.2, 3.6, 44, 1.4, 0.78, 54, 0.34, 0.66, 12, 0.7, 8.0, 12, 65, 43, 5, 76, 65, 8),
  nrows = 18,
  ncols = 1
)
npl <- as.matrix(c(0.2, 3.8, 14, 2.4, 1.7, 43, 0.2, 0.5, 23, 7.8, 88, 36, 65, 3, 44, 65, 7, 34),
  nrows = 18,
  ncols = 1
)
mdata <- data.frame(p = pib, t = tir, int = inf, np = npl)

# model
ivf <- ivfixed::ivFixed(t ~ p + int | p + np, mdata, n = 6, t = 3)

# plot
ggcoefstats(
  x = ivf,
  title = "instrumental fixed effect panel regression"
)
```

## causal mediation analysis (`mediate`)

```{r mediate, fig.height=5, fig.width=5}
# setup
set.seed(123)
library(mediation)
data(jobs)

# base models
b <-
  stats::lm(
    formula = job_seek ~ treat + econ_hard + sex + age,
    data = jobs
  )
c <-
  stats::lm(
    formula = depress2 ~ treat + job_seek + econ_hard + sex + age,
    data = jobs
  )

# mediation model
mod_mediate <-
  mediation::mediate(
    model.m = b,
    model.y = c,
    sims = 50,
    treat = "treat",
    mediator = "job_seek"
  )

# plot
ggcoefstats(
  x = mod_mediate,
  title = "causal mediation analysis"
)
```

## Model II regression (`lmodel2`)

```{r lmodel2, fig.height=6, fig.width=6}
# setup
set.seed(123)
library(lmodel2)
data(mod2ex2)

# model
Ex2.res <-
  lmodel2::lmodel2(
    formula = Prey ~ Predators,
    data = mod2ex2,
    range.y = "relative",
    range.x = "relative",
    nperm = 99
  )

# plot
ggcoefstats(
  x = Ex2.res,
  title = "Model II regression"
)
```

## Multivariate Ordinal Regression Models (`mvord`)

```{r mvord, fig.height=6, fig.width=5}
set.seed(123)
library(mvord)

# toy example
data(data_mvord_toy)

# wide data format with MMO2
res_mvord <- mvord(
  formula = MMO2(Y1, Y2) ~ 0 + X1 + X2,
  data = data_mvord_toy
)

# plot
ggcoefstats(
  x = res_mvord,
  title = "Multivariate Ordinal Regression Models"
)
```

## generalized additive models for location scale and shape (`gamlss`)

```{r gamlss, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(gamlss)

# model
g_gamlss <-
  gamlss::gamlss(
    formula = y ~ pb(x),
    sigma.fo = ~ pb(x),
    family = BCT,
    data = abdom,
    method = mixed(1, 20)
  )

# plot
ggcoefstats(
  x = g_gamlss,
  title = "generalized additive models \nfor location scale and shape"
)
```

## Bayesian Gaussian Graphical Model (`BGGM`)

```{r BGGM, fig.height=10, fig.width=6}
set.seed(123)
library(BGGM)

# data
Y <- dplyr::select(ptsd, B1:B5)

# fit model (note + 1, due to zeros)
fit_ggm <- estimate(Y + 1, type = "ordinal", iter = 250)

# plot
ggcoefstats(
  x = fit_ggm,
  title = "Bayesian Gaussian Graphical Model"
)
```

## Linear Equation System Estimation (`systemfit`)

```{r systemfit, fig.height=7, fig.width=8}
library(systemfit)

# data
data("Kmenta")
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list(demand = eqDemand, supply = eqSupply)
inst1 <- ~ income + farmPrice
inst2 <- ~ income + farmPrice + trend
instlist <- list(inst1, inst2)

# model
m_systemfit <- systemfit(
  system,
  "2SLS",
  inst = instlist,
  data = Kmenta
)

# plot
ggcoefstats(
  x = m_systemfit,
  title = "Linear Equation System Estimation"
)
```


## generalized method of moment estimation (`gmm`)

```{r gmm, fig.height=8, fig.width=5}
# setup
set.seed(123)
library(gmm)

# examples come from the "gmm" package
## CAPM test with GMM
data(Finance)
r <- Finance[1:300, 1:10]
rm <- Finance[1:300, "rm"]
rf <- Finance[1:300, "rf"]
z <- as.matrix(r - rf)
t <- nrow(z)
zm <- rm - rf
h <- matrix(zm, t, 1)
res_gmm <- gmm::gmm(z ~ zm, x = h)

# plot
ggcoefstats(
  x = res_gmm,
  package = "palettetown",
  palette = "victreebel",
  title = "generalized method of moment estimation"
)
```

## Spatial simultaneous autoregressive model estimation\n by maximum likelihood (`Sarlm`)

```{r Sarlm}
library(spatialreg)

# data
data(oldcol, package = "spdep")
listw <- spdep::nb2listw(COL.nb, style = "W")
ev <- eigenw(listw)
W <- as(listw, "CsparseMatrix")
trMatc <- trW(W, type = "mult")

# model
COL.lag.eig <- lagsarlm(
  CRIME ~ INC + HOVAL,
  data = COL.OLD,
  listw = listw,
  method = "eigen",
  quiet = FALSE,
  control = list(pre_eig = ev, OrdVsign = 1)
)

# plot
ggcoefstats(
  x = COL.lag.eig,
  title = "Spatial simultaneous autoregressive model estimation\n by maximum likelihood"
)
```


## fit a GLM with lasso or elasticnet regularization (`glmnet`)

Although these models are not directly supported in `ggcoefstats` because of the
sheer number of terms that are typically present. But this function can still be
used to selectively show few of the terms of interest:

```{r glmnet, fig.height=4, fig.width=4}
# setup
library(glmnet)
set.seed(2014)

# creating a dataframe
x <- matrix(rnorm(100 * 20), 100, 20)
y <- rnorm(100)
fit1 <- glmnet::glmnet(x, y)
(df <- broom::tidy(fit1))

# displaying only a certain step
ggcoefstats(x = dplyr::filter(df, step == 4))
```

## exponential-family random graph models (`ergm`)

```{r ergm, fig.height=5, fig.width=5}
# load the Florentine marriage network data
set.seed(123)
suppressPackageStartupMessages(library(ergm))
data(florentine)

# fit a model where the propensity to form ties between
# families depends on the absolute difference in wealth
gest <- ergm::ergm(flomarriage ~ edges + absdiff("wealth"))

# plot
ggcoefstats(
  x = gest,
  conf.level = 0.99,
  title = "exponential-family random graph models"
)
```

## TERGM by bootstrapped pseudolikelihood or MCMC MLE (`btergm`)

```{r btergm, fig.height=5, fig.width=6}
# setup
library(network)
library(btergm)
set.seed(123)

# create 10 random networks with 10 actors
networks <- list()
for (i in 1:10) {
  mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
  diag(mat) <- 0 # loops are excluded
  nw <- network(mat) # create network object
  networks[[i]] <- nw # add network to the list
}

# create 10 matrices as covariate
covariates <- list()
for (i in 1:10) {
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  covariates[[i]] <- mat # add matrix to the list
}

# model
fit_btergm <-
  btergm::btergm(
    formula = networks ~ edges + istar(2) + edgecov(covariates),
    parallel = "multicore",
    ncpus = 4,
    R = 100,
    verbose = FALSE
  )

# plot
ggcoefstats(
  x = fit_btergm,
  title = "Terms used in Exponential Family Random Graph Models",
  subtitle = "by bootstrapped pseudolikelihood or MCMC MLE"
)
```

## generalized autoregressive conditional heteroscedastic regression (`garch`)

```{r garch, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(tseries)
data(EuStockMarkets)

# model
dax <- diff(log(EuStockMarkets))[, "DAX"]
dax.garch <- tseries::garch(x = dax, trace = FALSE)

# plot
ggcoefstats(
  x = dax.garch,
  title = "generalized autoregressive conditional heteroscedastic"
)
```

## Bayesian Estimation of the GARCH(1,1) Model with Student-t Innovations (`bayesGARCH`)

```{r bayesGARCH, fig.height=4, fig.width=5}
# setup
set.seed(123)
library(bayesGARCH)

# data
data(dem2gbp)
y <- dem2gbp[1:750]

# run the sampler (2 chains)
mod_bayesGARCH <- bayesGARCH::bayesGARCH(y,
  control = list(n.chain = 4, l.chain = 200, refresh = 1000)
)

# plot
ggcoefstats(
  x = mod_bayesGARCH,
  title = "Bayesian Estimation of the GARCH(1,1) Model \nwith Student-t Innovations"
)
```

## maximum-likelihood fitting of univariate distributions (`fitdistr`)

```{r fitdistr, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(MASS)
x <- rnorm(100, 5, 2)

# model
fit_fitdistr <- MASS::fitdistr(x, dnorm, list(mean = 3, sd = 1))

# plot
ggcoefstats(
  x = fit_fitdistr,
  title = "maximum-likelihood fitting of univariate distributions",
  ggtheme = ggthemes::theme_pander()
)
```

## maximum likelihood estimation (`mle2`)

```{r mle2, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(bbmle)

# data
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x, y)

# custom function
LL <- function(ymax = 15, xhalf = 6) {
  -sum(stats::dpois(y, lambda = ymax / (1 + x / xhalf), log = TRUE))
}

# use default parameters of LL
fit_mle2 <- bbmle::mle2(LL, fixed = list(xhalf = 6))

# plot
ggcoefstats(
  x = fit_mle2,
  title = "maximum likelihood estimation",
  ggtheme = ggthemes::theme_excel_new()
)
```

## General Linear Hypotheses (`glht`)

```{r glht, fig.height=5, fig.width=6}
# setup
set.seed(123)
library(multcomp)

# multiple linear model, swiss data
lmod <- lm(Fertility ~ ., data = swiss)

# model
mod_glht <-
  multcomp::glht(
    model = lmod,
    linfct = c(
      "Agriculture=0",
      "Examination=0",
      "Education=0",
      "Catholic=0",
      "Infant.Mortality=0"
    )
  )

# plot
ggcoefstats(
  x = mod_glht,
  title = "General Linear Hypotheses (using `multcomp`)"
)
```

## Cochrane-Orcutt estimation (`orcutt`)

```{r orcutt, fig.height=5, fig.width=5}
# setup
library(orcutt)
set.seed(123)

# model
reg <- stats::lm(formula = mpg ~ wt + qsec + disp, data = mtcars)
co <- orcutt::cochrane.orcutt(reg)

# plot
ggcoefstats(
  x = co,
  title = "Cochrane-Orcutt estimation"
)
```

## confusion matrix (`confusionMatrix`)

```{r confusionMatrix, fig.height=4, fig.width=4}
# setup
library(caret)
set.seed(123)

# setting up confusion matrix
two_class_sample1 <- as.factor(sample(letters[1:2], 100, TRUE))
two_class_sample2 <- as.factor(sample(letters[1:2], 100, TRUE))

two_class_cm <- caret::confusionMatrix(two_class_sample1, two_class_sample2)

# plot
ggcoefstats(
  x = two_class_cm,
  by.class = TRUE,
  title = "confusion matrix"
)
```

## dose-response models/curves (`drc`)

```{r drc, fig.height=7, fig.width=5}
# setup
set.seed(123)
library(drc)

# model
mod_drc <-
  drc::drm(
    formula = dead / total ~ conc,
    curveid = type,
    weights = total,
    data = selenium,
    fct = LL.2(),
    type = "binomial"
  )

# plot
ggcoefstats(
  x = mod_drc,
  conf.level = 0.99,
  title = "Dose-Response Curves"
)
```

## Firth's bias-reduced logistic regression (`logitsf`)

```{r logitsf, fig.height=4, fig.width=4}
# setup
set.seed(123)
library(logistf)

# dataframe
data <- data.frame(time = 1:12, event = c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1))

# model
mod_logistf <- logistf::logistf(event ~ time, data)

# plot
ggcoefstats(
  x = mod_logistf,
  title = "Firth's bias-reduced logistic regression"
)
```

## Flexible parametric regression for time-to-event data (`flexsurvreg`)

```{r flexsurvreg, fig.height=6, fig.width=7}
# setup
set.seed(123)
library(flexsurv)
data(ovarian)

# compare generalized gamma fit with Weibull
fitg <-
  flexsurv::flexsurvreg(
    formula = Surv(futime, fustat) ~ 1,
    data = ovarian,
    dist = "gengamma"
  )

ggcoefstats(
  x = fitg,
  title = "Flexible parametric regression for time-to-event data"
)
```

## panel regression models fit via multilevel modeling (`wblm`)

```{r wblm, fig.height=6, fig.width=7}
# setup
set.seed(123)
library(panelr)

# data
data("WageData")
wages <- panelr::panel_data(data = WageData, id = id, wave = t)

# model
mod_wblm <-
  panelr::wbm(
    formula = lwage ~ lag(union) + wks | blk + fem | blk * lag(union),
    data = wages
  )

# plot
ggcoefstats(
  x = mod_wblm,
  title = "panel regression models fit via multilevel modeling"
)
```

## panel regression models fit via generalized estimating equations (`wbgee`)

```{r wbgee, fig.height=8, fig.width=6}
# setup
set.seed(123)
library(panelr)
data("WageData")

# data
wages <- panelr::panel_data(data = WageData, id = id, wave = t)

# model
mod_wbgee <-
  panelr::wbgee(
    formula = lwage ~ lag(union) + wks | blk + fem | blk * lag(union),
    data = wages
  )

# plot
ggcoefstats(
  x = mod_wbgee,
  conf.level = 0.99,
  title = "panel regression models fit via generalized estimating equations"
)
```

## Censored Regression with Conditional Heteroscedasticy (`crch`)

```{r crch, fig.height=4, fig.width=8}
# setup
set.seed(123)
library(crch)
data("RainIbk")

# mean and standard deviation of square root transformed ensemble forecasts
RainIbk$sqrtensmean <- apply(sqrt(RainIbk[, grep("^rainfc", names(RainIbk))]), 1, mean)
RainIbk$sqrtenssd <- apply(sqrt(RainIbk[, grep("^rainfc", names(RainIbk))]), 1, sd)

# fit linear regression model with Gaussian distribution
CRCH <-
  crch::crch(
    formula = sqrt(rain) ~ sqrtensmean,
    data = RainIbk,
    dist = "gaussian"
  )

# plot
ggcoefstats(
  x = CRCH,
  title = "Censored Regression with Conditional Heteroscedasticy"
)
```

## Marginal Models For Correlated Ordinal Multinomial Responses (`LORgee`)

```{r LORgee, fig.height=10, fig.width=8}
# setup
set.seed(123)
library(multgee)
data("arthritis")

# model
fit_LORgee <-
  multgee::ordLORgee(
    formula = y ~ factor(time) + factor(trt) + factor(baseline),
    link = "logit",
    id = id,
    repeated = time,
    data = arthritis,
    LORstr = "uniform"
  )

# plot
ggcoefstats(
  x = fit_LORgee,
  title = "Marginal Models For Correlated Ordinal Multinomial Responses"
)
```

## summary measures for count data (`epi.2by2`)

```{r epi.2by2, fig.height=4, fig.width=6}
# setup
set.seed(123)
library(epiR)

# data
dat <- matrix(c(13, 2163, 5, 3349), nrow = 2, byrow = TRUE)
rownames(dat) <- c("DF+", "DF-")
colnames(dat) <- c("FUS+", "FUS-")

# model
fit_epiR <-
  epiR::epi.2by2(
    dat = as.table(dat),
    method = "cross.sectional",
    conf.level = 0.95,
    units = 100,
    outcome = "as.columns"
  )

# tidy dataframe
ggcoefstats(
  x = fit_epiR,
  title = "Summary measures for count data presented in a 2 by 2 table"
)
```

## marginal effects for a probit regression (`probitmfx`)

```{r probitmfx, fig.height=4, fig.width=5}
# setup
set.seed(123)
library(mfx)

# simulate some data
set.seed(12345)
n <- 1000
x <- rnorm(n)

# binary outcome
y <- ifelse(pnorm(1 + 0.5 * x + rnorm(n)) > 0.5, 1, 0)
data <- data.frame(y, x)

# model
mod_probitmfx <- mfx::probitmfx(formula = y ~ x, data = data)

# plot
ggcoefstats(
  x = mod_probitmfx,
  title = "marginal effects for a probit regression"
)
```

## marginal effects for a beta regression (`betamfx`)

```{r betamfx, fig.height=4, fig.width=5}
# setup
set.seed(123)
library(mfx)

# simulate some data
n <- 1000
x <- rnorm(n)

# beta outcome
y <- rbeta(n, shape1 = plogis(1 + 0.5 * x), shape2 = (abs(0.2 * x)))

# use Smithson and Verkuilen correction
y <- (y * (n - 1) + 0.5) / n

# data
d <- data.frame(y, x)

# model
mod_betamfx <- mfx::betamfx(y ~ x | x, data = d)

# plot
ggcoefstats(
  x = mod_betamfx,
  title = "marginal effects for a beta regression"
)
```

## marginal effects for a logit regression (`logitmfx`)

Same example can also be used for `logitor` object.

```{r logitmfx, fig.height=4, fig.width=5}
# setup
set.seed(123)
library(mfx)

# simulate some data
n <- 1000
x <- rnorm(n)

# binary outcome
y <- ifelse(pnorm(1 + 0.5 * x + rnorm(n)) > 0.5, 1, 0)

# data
data <- data.frame(y, x)

# model
mod_logitmfx <- logitmfx(formula = y ~ x, data = data)

# plot
ggcoefstats(
  x = mod_logitmfx,
  title = "marginal effects for a logit regression"
)
```

## marginal effects for a negative binomial regression (`negbinmfx`)

Same example can also be used for `negbinirr` object.

```{r negbinmfx, fig.height=4, fig.width=5}
# simulate some data
set.seed(123)
n <- 1000

# data
x <- rnorm(n)
y <- rnegbin(n, mu = exp(1 + 0.5 * x), theta = 0.5)
data <- data.frame(y, x)

# model
mod_negbinmfx <- negbinmfx(formula = y ~ x, data = data)

# plot
ggcoefstats(
  x = mod_negbinmfx,
  title = "marginal effects for a negative binomial regression"
)
```

## marginal effects for a Poisson regression (`poissonmfx`)

Same example can also be used for `poissonirr` object.

```{r poissonmfx, fig.height=4, fig.width=5}
# simulate some data
set.seed(123)
n <- 1000
x <- rnorm(n)
y <- rnegbin(n, mu = exp(1 + 0.5 * x), theta = 0.5)
data <- data.frame(y, x)

# model
mod_poissonmfx <- poissonmfx(formula = y ~ x, data = data)

# plot
ggcoefstats(
  x = mod_poissonmfx,
  title = "marginal effects for a Poisson regression"
)
```

## Bayesian generalized (non-)linear multivariate multilevel models (`brmsfit`)

```{r brmsfit, fig.height=5, fig.width=7}
# setup
set.seed(123)
library(brms)

# prior
bprior1 <- prior(student_t(5, 0, 10), class = b) + prior(cauchy(0, 2), class = sd)

# model
fit_brm <-
  brms::brm(
    formula = count ~ Age + Base * Trt + (1 | patient),
    data = epilepsy,
    family = poisson(),
    prior = bprior1,
    silent = TRUE,
    refresh = 0
  )

# plot
ggcoefstats(
  x = fit_brm,
  standardize = TRUE,
  title = "Bayesian generalized (non-)linear multivariate multilevel models",
  subtitle = "using `brms` package"
)
```

Let's see another example where we use `brms` to run the same model on multiple
datasets-

```{r brmsfit_multiple, fig.height=5, fig.width=7}
# setup
set.seed(123)
library(brms)
library(mice)

# data
imp <- mice::mice(data = nhanes2, print = FALSE)

# fit the model using brms
fit_brm_m <-
  brms::brm_multiple(
    formula = bmi ~ age + hyp + chl,
    data = imp,
    chains = 1,
    refresh = 0
  )

# plot
ggcoefstats(
  x = fit_brm_m,
  title = "Same `brms` model with multiple datasets",
  conf.level = 0.99
)
```

## Bayesian nonlinear models via Stan (`nlstanreg`)

```{r nlstanreg, fig.height=5, fig.width=5, error=TRUE}
# setup
set.seed(123)
library(rstanarm)
data("Orange", package = "datasets")
Orange$circumference <- Orange$circumference / 100
Orange$age <- Orange$age / 100

# model
fit_nlstanreg <-
  rstanarm::stan_nlmer(
    formula = circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym | Tree,
    data = Orange,
    # for speed only
    chains = 1,
    iter = 1000,
    refresh = 0
  )

# plot
ggcoefstats(
  x = fit_nlstanreg,
  title = "Bayesian nonlinear models via Stan"
)
```

## Markov Chain Monte Carlo objects (`mcmc`)

```{r mcmc, fig.height=5, fig.width=5}
# loading the data
set.seed(123)
library(coda)
data(line)

# select first chain
x1 <- line[[1]]

# plot
ggcoefstats(
  x = x1,
  title = "Markov Chain Monte Carlo objects",
  robust = TRUE,
  ess = TRUE
)
```

## MCMC with Just Another Gibbs Sampler (`rjags`)

```{r R2jags, fig.height=8, fig.width=6}
# setup
set.seed(123)
library(R2jags)

# An example model file is given in:
model.file <- system.file(package = "R2jags", "model", "schools.txt")

# data
J <- 8.0
y <- c(28.4, 7.9, -2.8, 6.8, -0.6, 0.6, 18.0, 12.2)
sd <- c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)

# setting up model
jags.data <- list("y", "sd", "J")
jags.params <- c("mu", "sigma", "theta")
jags.inits <- function() {
  list("mu" = rnorm(1), "sigma" = runif(1), "theta" = rnorm(J))
}

# model fitting
jagsfit <-
  R2jags::jags(
    data = list("y", "sd", "J"),
    inits = jags.inits,
    jags.params,
    n.iter = 10,
    model.file = model.file
  )

# plot
ggcoefstats(
  x = jagsfit,
  title = "Markov Chain Monte Carlo with Just Another Gibbs Sampler (JAGS)"
)
```

## latent variable model (`lavaan`)

```{r lavaan, fig.height=10, fig.width=6}
# setup
set.seed(123)
library(lavaan)

# The Holzinger and Swineford (1939) example
HS.model <- " visual =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9 "

# model
fit_lavaan <-
  lavaan::lavaan(
    HS.model,
    data = HolzingerSwineford1939,
    auto.var = TRUE,
    auto.fix.first = TRUE,
    auto.cov.lv.x = TRUE
  )

# tidy
ggcoefstats(
  x = fit_lavaan,
  title = "latent variable model",
  conf.level = 0.99
)
```

## Bayesian latent variable model (`blavaan`)

Too time-consuming to be run here, so currently skipped.

```{r blavaan, fig.height=10, fig.width=6, eval=FALSE}
# setup
set.seed(123)
library(blavaan)

# The Holzinger and Swineford (1939) example
HS.model <- " visual =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9 "

# model
fit_blavaan <-
  blavaan::blavaan(
    HS.model,
    data = HolzingerSwineford1939,
    auto.var = TRUE,
    auto.fix.first = TRUE,
    auto.cov.lv.x = TRUE
  )

# tidy
ggcoefstats(
  x = fit_blavaan,
  title = "Bayesian latent variable model"
)
```

## Bayesian regression models via Stan (`stanreg`)

```{r stanreg, fig.height=5, fig.width=5}
# set up
set.seed(123)
library(rstanarm)

# model
fit_stanreg <-
  rstanarm::stan_glm(
    formula = mpg ~ wt + am,
    data = mtcars,
    chains = 1,
    refresh = 0
  )

# plot
ggcoefstats(
  x = fit_stanreg,
  title = "Bayesian generalized linear models via Stan"
)
```

## Bayesian multivariate generalized linear models via Stan (`stanmvreg`)

```{r stanmvreg, fig.height=5, fig.width=5, error=TRUE}
# setup
set.seed(123)
library(rstanarm)
pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))

# model
mod_stanmvreg <-
  rstanarm::stan_mvmer(
    formula = list(
      ybern ~ year + (1 | id),
      albumin ~ sex + year + (year | id)
    ),
    data = pbcLong,
    cores = 4,
    seed = 12345,
    iter = 1000,
    refresh = 0
  )

# plot
ggcoefstats(
  x = mod_stanmvreg,
  title = "Bayesian multivariate generalized linear models \nvia Stan"
)
```

## Bayesian quantile regression (`bayesQR`)

```{r bayesQR, fig.height=6, fig.width=6}
# needed library
library(bayesQR)

# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n = n, min = 0, max = 10)
X <- X
y <- 1 + 2 * X + rnorm(n = n, mean = 0, sd = .6 * X)

# model
mod_bayesqr <-
  bayesQR(y ~ X,
    quantile = c(.05, .25, .5, .75, .95),
    alasso = TRUE,
    ndraw = 500
  )

# plot
ggcoefstats(
  x = mod_bayesqr,
  title = "Bayesian quantile regression"
)
```

## Bayesian binary and ordinal logistic regression (`blrm`)

```{r blrm, fig.height=8, fig.width=10}
# setup
set.seed(123)
library(rmsb)
library(rms)
library(Hmisc)
getHdata(titanic3)
dd <- datadist(titanic3)
options(datadist = "dd")

# model
mod_blrm <- rmsb::blrm(survived ~ (rcs(age, 5) + sex + pclass)^2, data = titanic3)

# plot
ggcoefstats(
  x = mod_blrm,
  title = "Bayesian binary and ordinal logistic regression"
)
```

## Compound Poisson Generalized Linear Models (`cpglm`)

```{r cpglm, fig.height=6, fig.width=10}
# set up
set.seed(123)
library(cplm)

# model
mod_cpglm <-
  cpglm(
    formula = RLD ~ factor(Zone) * factor(Stock),
    data = FineRoot
  )

# plot
ggcoefstats(
  x = mod_cpglm,
  title = "Compound Poisson Generalized Linear Models"
)
```

## Compound Poisson Generalized Linear Mixed Models (`cpglmm`)

```{r cpglmm, fig.height=6, fig.width=10}
# set up
set.seed(123)
library(cplm)

# model
mod_cpglmm <-
  cpglmm(
    formula = RLD ~ Stock + Spacing + (1 | Plant),
    data = FineRoot
  )

# plot
ggcoefstats(
  x = mod_cpglmm,
  title = "Compound Poisson Generalized Linear Mixed Models"
)
```

## Bayesian Compound Poisson Linear Models (`bcplm`)

```{r bcplm, fig.height=6, fig.width=10}
# set up
set.seed(123)
library(cplm)

# model
mod_bcplm <-
  bcplm(
    formula = RLD ~ factor(Zone) * factor(Stock),
    data = FineRoot,
    tune.iter = 2000,
    n.iter = 6000,
    n.burnin = 1000,
    n.thin = 5,
    n.report = 0
  )

# plot
ggcoefstats(
  x = mod_bcplm,
  title = "Bayesian Compound Poisson Linear Models"
)
```

## Multivariate Imputation by Chained Equations (`mice`)

```{r mice, fig.height=6, fig.width=6}
library(mice)
data(nhanes2)
imp <- mice(nhanes2)

# plot
ggcoefstats(
  x = with(data = imp, exp = lm(bmi ~ age + hyp + chl)),
  title = "Multivariate Imputation by Chained Equations"
)
```

## estimations with multiple fixed-effects (`fixest`)

```{r fixest, fig.height=10, fig.width=10}
# setup
library(fixest)

# data
data(trade)

# generating data for a simple example
set.seed(1)
n <- 100
x <- rnorm(n, 1, 5)**2
y <- rnorm(n, -1, 5)**2
z1 <- rpois(n, x * y) + rpois(n, 2)
base <- data.frame(x, y, z1)

# creating different `fixest` regression model objects
mod_fx1 <-
  fixest::femlm(
    fml = Euros ~ log(dist_km) | Origin + Destination + Product,
    data = trade
  )

mod_fx2 <-
  fixest::fepois(
    fml = Sepal.Length ~ Sepal.Width + Petal.Length | Species,
    data = iris
  )

mod_fx3 <-
  fixest::feols(
    fml = Sepal.Length ~ Sepal.Width + Petal.Length | Species,
    data = iris
  )

mod_fx4 <-
  fixest::feNmlm(
    z1 ~ 1, base,
    NL.fml = ~ a * log(x) + b * log(y),
    NL.start = list(a = 0, b = 0)
  )

# combining plots
combine_plots(
  plotlist = list(
    ggcoefstats(mod_fx1, title = "femlm"),
    ggcoefstats(mod_fx2, title = "fepois"),
    ggcoefstats(mod_fx3, title = "feols"),
    ggcoefstats(mod_fx4, title = "feNmlm")
  ),
  labels = "auto",
  annotation.args = list(title = "fast and user-friendly fixed-effects estimations")
)
```

## multivariate analysis of variance (`manova`)

```{r manova, fig.height=8, fig.width=8}
# setup
set.seed(123)

# fake a 2nd response variable
npk2 <- within(npk, foo <- rnorm(24))

# model
m_manova <- manova(cbind(yield, foo) ~ block + N * P * K, npk2)

# plot
ggcoefstats(
  x = m_manova,
  title = "multivariate analysis of variance"
)
```

# Regression models with `list` outputs

Note that a number of regression models will return an object of class `list`,
in which case this function will fail. But often you can extract the object of
interest from this list and use it to plot the regression coefficients.

```{r list_gamm4, fig.height=6, fig.width=10}
# setup
library(gamm4)
set.seed(123)

# data
dat <- gamSim(1, n = 400, scale = 2)

# now add 20 level random effect `fac'...
dat$fac <- fac <- as.factor(sample(1:20, 400, replace = TRUE))
dat$y <- dat$y + model.matrix(~ fac - 1) %*% rnorm(20) * .5

# model object
br <-
  gamm4::gamm4(
    formula = y ~ s(x0) + x1 + s(x2),
    data = dat,
    random = ~ (1 | fac)
  )

# looking at the classes of the objects contained in the list
purrr::map(br, class)

# plotting
combine_plots(
  plotlist = list(
    # first object plot (only parametric terms are shown)
    ggcoefstats(
      x = br$gam,
      title = "generalized additive model (parametric terms)",
      k = 3
    ),
    # second object plot
    ggcoefstats(
      x = br$mer,
      title = "linear mixed-effects model",
      k = 3
    )
  ),
  plotgrid.args = list(nrow = 1)
)
```

# Meta-analysis

In case the estimates you are displaying come from multiple studies, you can
also use this function to carry out random-effects meta-analysis. The dataframe
you enter **must** contain at the minimum the following three columns-

  - `term`: a column with names/identifiers to annotate each study/effect

  - `estimate`: a column with the observed effect sizes or outcomes

  - `std.error`: a column the corresponding standard errors

### parametric

```{r meta1, fig.height=7, fig.width=8}
# setup
set.seed(123)
library(metaplus)

# renaming to what the function expects
df <- dplyr::rename(mag, estimate = yi, std.error = sei, term = study)

# plot
ggcoefstats(
  x = df,
  meta.analytic.effect = TRUE,
  bf.message = TRUE,
  meta.type = "parametric",
  title = "parametric random-effects meta-analysis"
)
```

### robust

```{r meta2, fig.height=7, fig.width=8}
# setup
set.seed(123)
library(metaplus)

# renaming to what the function expects
df <- dplyr::rename(mag, estimate = yi, std.error = sei, term = study)

# plot
ggcoefstats(
  x = df,
  meta.analytic.effect = TRUE,
  meta.type = "robust",
  title = "robust random-effects meta-analysis"
)
```

### Bayesian

```{r meta3, fig.height=7, fig.width=8}
# setup
set.seed(123)
library(metaplus)

# renaming to what the function expects
df <- dplyr::rename(mag, estimate = yi, std.error = sei, term = study)

# plot
ggcoefstats(
  x = df,
  meta.analytic.effect = TRUE,
  meta.type = "bayes",
  title = "Bayesian random-effects meta-analysis"
)
```

# Dataframes

Sometimes you don't have a model object but a custom dataframe that you want
display using this function. If a data frame is to be plotted, it **must**
contain columns named `term` (names of predictors), and `estimate`
(corresponding estimates of coefficients or other quantities of interest). Other
optional columns are `conf.low` and `conf.high` (for confidence intervals), and
`p.value`. You will also have to specify the type of statistic relevant for
regression models (`"t"`, `"z"`, `"f"`, `"chi"`) in case you want to display
statistical labels.

You can also provide a dataframe containing all the other relevant information
for additionally displaying labels with statistical information.

```{r dataframe, fig.height=7, fig.width=7}
# let's make up a dataframe (with all available details)
df_full <-
  tibble::tribble(
    ~term, ~statistic, ~estimate, ~std.error, ~p.value, ~df.error,
    "study1", 0.158, 0.0665, 0.778, 0.875, 5L,
    "study2", 1.33, 0.542, 0.280, 0.191, 10L,
    "study3", 1.24, 0.045, 0.030, 0.001, 12L,
    "study4", 0.156, 0.500, 0.708, 0.885, 8L,
    "study5", 0.33, 0.032, 0.280, 0.101, 2L,
    "study6", 1.04, 0.085, 0.030, 0.001, 3L
  )

# plot
ggcoefstats(
  x = df_full,
  meta.analytic.effect = TRUE,
  statistic = "t",
  package = "LaCroixColoR",
  palette = "paired"
)
```

# Other types of outputs

This function can also be used to extract outputs other than a plot, although it
is much more preferable to use the underlying functions instead
(`parameters::model_parameters`).

```{r other_output, error=TRUE}
# setup
set.seed(123)
library(ggstatsplot)

# data
DNase1 <- subset(DNase, Run == 1)

# using a selfStart model
nlmod <- stats::nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)

# tidy dataframe
ggcoefstats(nlmod, output = "tidy")

# glance summary
ggcoefstats(nlmod, output = "glance")
```

# Not supported

This vignette was supposed to give a comprehensive account of regression models
supported by `ggcoefstats`. The list of supported models will keep expanding as
additional tidiers are added to the `parameters` and `performance` packages.

Note that not **all** models supported in these packages will be supported by
`ggcoefstats`. In particular, classes of objects for which there is no column
for `estimate` (e.g., `kmeans`, `optim`, `muhaz`, `survdiff`, `zoo`, etc.) are
not supported.

# Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "Benchmarking"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_document:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    toc_depth: 3
    warning: FALSE
    message: FALSE
---

```{r setup, include = FALSE}
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  out.width = "100%",
  comment = "#>",
  warning = FALSE,
  message = FALSE
)
```

# Introduction

This is mostly to keep track of how the performance of different functions
changes across time.

# `ggbetweenstats`

```{r ggbetweenstats}
library(ggstatsplot)

set.seed(123)
bench::mark(
  ggbetweenstats(
    data = dplyr::filter(
      movies_long,
      genre %in% c("Action", "Action Comedy", "Action Drama", "Comedy")
    ),
    x = mpaa,
    y = length
  )
) %>%
  dplyr::select(-expression)

set.seed(123)
bench::mark(
  grouped_ggbetweenstats(
    data = dplyr::filter(
      movies_long,
      genre %in% c("Action", "Action Comedy", "Action Drama", "Comedy")
    ),
    x = mpaa,
    y = length,
    grouping.var = genre
  )
) %>%
  dplyr::select(-expression)
```

# `ggwithinstats`

```{r ggwithinstats}
library(ggstatsplot)

set.seed(123)
bench::mark(
  ggwithinstats(
    data = bugs_long,
    x = condition,
    y = desire
  )
) %>%
  dplyr::select(-expression)

set.seed(123)
bench::mark(
  grouped_ggwithinstats(
    data = bugs_long,
    x = condition,
    y = desire,
    grouping.var = gender
  )
) %>%
  dplyr::select(-expression)
```

# `gghistostats`

```{r gghistostats}
library(ggstatsplot)

set.seed(123)
bench::mark(gghistostats(mtcars, wt, test.value = 3)) %>%
  dplyr::select(-expression)

set.seed(123)
bench::mark(grouped_gghistostats(mtcars, wt, test.value = 3, grouping.var = am)) %>%
  dplyr::select(-expression)
```

# `ggdotplotstats`

```{r ggdotplotstats}
library(ggstatsplot)

set.seed(123)
bench::mark(ggdotplotstats(
  dplyr::filter(ggplot2::mpg, cyl %in% c("4", "6")),
  cty,
  manufacturer,
  test.value = 15
)) %>%
  dplyr::select(-expression)

set.seed(123)
bench::mark(
  grouped_ggdotplotstats(
    dplyr::filter(ggplot2::mpg, cyl %in% c("4", "6")),
    cty,
    manufacturer,
    test.value = 15,
    grouping.var = cyl
  )
) %>%
  dplyr::select(-expression)
```

# `ggscatterstats`

```{r ggscatterstats}
library(ggstatsplot)

set.seed(123)
bench::mark(ggscatterstats(mtcars, wt, mpg)) %>%
  dplyr::select(-expression)

set.seed(123)
bench::mark(grouped_ggscatterstats(mtcars, wt, mpg, grouping.var = am)) %>%
  dplyr::select(-expression)
```

# `ggcorrmat`

```{r ggcorrmat}
library(ggstatsplot)

set.seed(123)
bench::mark(ggcorrmat(iris)) %>%
  dplyr::select(-expression)

set.seed(123)
bench::mark(grouped_ggcorrmat(iris, grouping.var = Species)) %>%
  dplyr::select(-expression)
```

# `ggpiestats`

```{r ggpiestats}
library(ggstatsplot)

set.seed(123)
bench::mark(ggpiestats(mtcars, cyl)) %>%
  dplyr::select(-expression)

set.seed(123)
bench::mark(grouped_ggpiestats(mtcars, cyl, grouping.var = am)) %>%
  dplyr::select(-expression)
```

# `ggbarstats`

```{r ggbarstats}
library(ggstatsplot)

set.seed(123)
bench::mark(ggbarstats(ggplot2::mpg, fl, class)) %>%
  dplyr::select(-expression)

set.seed(123)
bench::mark(grouped_ggbarstats(ggplot2::mpg, fl, class, grouping.var = drv)) %>%
  dplyr::select(-expression)
```

# `ggcoefstats`

```{r ggcoefstats}
library(ggstatsplot)

set.seed(123)
bench::mark(ggcoefstats(stats::lm(formula = wt ~ am * cyl, data = mtcars))) %>%
  dplyr::select(-expression)
```

# Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
GitHub: <https://github.com/IndrajeetPatil/ggstatsplot/issues>

---
title: "other resources to learn about `{ggstatsplot}`"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 100
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{other resources to learn about `{ggstatsplot}`}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)
```

# Introduction

One of the best resources to learn about functionality afforded by `{ggstatsplot}`
package is to read the dedicated
[website](https://indrajeetpatil.github.io/ggstatsplot/). A lot of love and work
has gone into writing detailed function documentation and vignettes, and I will
highly recommend you check them out.

But, the website is certainly not the only resource available to learn about
this package, and below are some other resources that might be helpful for you
to learn about `{ggstatsplot}` package.

# Presentations

  - My own `xaringan` presentation:
    <https://indrajeetpatil.github.io/ggstatsplot_slides/slides/ggstatsplot_presentation.html#1>

# YouTube Videos

## My own talks at online meetings
  
  - [At UseR Oslo meeting](https://www.youtube.com/watch?v=_wCo7rnt6NA)

  - [At R-Ladies Tunis meeting](https://www.youtube.com/watch?v=n-sXeshyVRI)

## Other

### English

  - [Yury Zablotski's review](https://www.youtube.com/watch?v=cNJIMwtLWv0)

  - [Business Science R tips](https://www.youtube.com/watch?v=8Em1bCFBMWg)

  - <https://www.youtube.com/watch?v=aw6HDz9zA1Q>

  - [Serdar Belci's demo](https://www.youtube.com/watch?v=m3uInetiC8w) on how to
    use `{ggstatsplot}` in `jamovi` GUI software

  - Data pre-processing and analysis in
    [RStudio](https://www.youtube.com/watch?v=IMgv5pc7N7U) and `{ggstatsplot}`

### Other languages

  - [Pablo Vallejo Medina's demo](https://www.youtube.com/watch?v=6Ds9L4NpRVQ)
    of `ggbetweenstats` (in **Spanish**)

  - [Andre Chocó Cedillos's talk](https://www.youtube.com/watch?v=yDyrgvvWLpc)
    (in **Spanish**)

# Blogs

This is a bit difficult to keep track of. I would recommend just googling it.

# Publications

A list of publications that use `{ggstatsplot}` can be found in Google Scholar:

<https://scholar.google.com/scholar?hl=en&scisbd=2&as_sdt=0%2C22&q=ggstatsplot&btnG=>

---
title: "using `{ggstatsplot}` with the `{purrr}` package"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_document:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
---

```{r setup, include = FALSE}
## pretty tibble printing
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  warning = FALSE,
  message = FALSE,
  out.width = "100%",
  comment = "#>",
  package.startup.message = FALSE
)

future::plan("multicore")
```

---

You can cite this package/vignette as:

```{r citation, echo=FALSE, comment = ""}
citation("ggstatsplot")
```

---

## Why use `{purrr}`?

Most of the `{ggstatsplot}` functions have `grouped_` variants, which are designed
to quickly run the same `{ggstatsplot}` function across multiple levels of a
**single** grouping variable. Although this function is useful for data
exploration, it has two strong weaknesses-

 * The arguments applied to `grouped_` function call are applied uniformly to
   all levels of the grouping variable when we might want to customize them for
   different levels of the grouping variable.

 * Only one grouping variable can be used to repeat the analysis when in reality
   there can be a combination of grouping variables and the operation needs to
   be repeated for all resulting combinations.

We will see how to overcome this limitation by combining `{ggstatsplot}` with the
`{purrr}` package.

**Note:**

 * While using `purrr::pmap()`, we **must** input the arguments as strings.

 * You can use `{ggplot2}` themes from extension packages (e.g. `ggthemes`).

 * If you'd like some more background or an introduction to the purrr package,
   please see [this chapter](https://adv-r.hadley.nz/functionals.html).

## Introduction and methodology

For all the examples in this vignette we are going to build `list`s of things
that we will pass along to `{purrr}` which will in turn return a list of plots
that will be passed to `combine_plots`. As the name implies `combine_plots`
merges the individual plots into one bigger plot with common labeling and
aesthetics.

What are these `lists` that we are building? The lists correspond to the
parameters in our `{ggstatsplot}` function like `ggbetweenstats`. If you look at
the help file for `?ggbetweenstats` for example the very first parameter it
wants is the `data` file we'll be using. We can also pass it different `titles`
of even `ggtheme` themes.

You can pass:

 * A single character string such as `xlab = "Continent"` or numeric such as
   `nboot = 25` in which case it will be reused/recycled as many times as
   needed.

 * A vector of values such as `nboot = c(50, 100, 200)` in which case it
   will be coerced to a list and checked for the right class (in this case
   integer) and the right quantity of entries in the vector i.e.,
   `nboot = c(50, 100)` will fail if we're trying to make three plots.

 * A list; either named `data = year_list` or created as you go
   `palette = list("Dark2", "Set1")`. Any list will
   be checked for the right class (in this case character) and the right
   quantity of entries in the list.

## `ggbetweenstats`

Let's start with `ggebtweenstats`. We'll use the `gapminder` dataset. We'll make
a 3 item `list` called `year_list` using `dplyr::filter` and `split`.

```{r purrr_ggbetweenstats1, fig.height = 12, fig.width = 7}
## for reproducibility
library(ggstatsplot)
set.seed(123)

## let's split the dataframe and create a list by years of interest
year_list <- gapminder::gapminder %>%
  dplyr::filter(year %in% c(1967, 1987, 2007), continent != "Oceania") %>%
  split(f = .$year, drop = TRUE)

## checking the length of the list and the names of each element
length(year_list)
names(year_list)
```

Now that we have the data divided into the three relevant years in a list we'll
turn to `purrr::pmap` to create a list of `ggplot` objects that we'll make use of
stored in `plot_list`. When you look at the documentation for `?pmap` it will
accept `.l` which is a list of lists. The length of `.l` determines the number of
arguments that `.f` will be called with. List names will be used if present.
`.f` is the function we want to apply (here, `.f = ggbetweenstats`).

Let's keep building the list of arguments, `.l`. First is `data = year_list`,
the `x` and `y` axes are constant in all three plots so we pass the variable
name as a string `x = "continent"`.

Same with the label we'll use for outliers where needed. For demonstration
purposes let's assume we want the outliers on each plot to be a different color.
Not actually recommending it just demonstrating what's possible. The rest of the
code shows you a wide variety of possibilities and we won't catalog them here.

```{r purrr_ggbetweenstats2, fig.height = 10, fig.width = 7}
## for reproducibility
set.seed(123)
library(ggstatsplot)

## creating a list of plots
plot_list <- purrr::pmap(
  .l = list(
    data = year_list,
    x = "continent",
    y = "lifeExp",
    outlier.tagging = TRUE,
    outlier.label = "country",
    outlier.label.args = list(
      list(size = 3, color = "#56B4E9"),
      list(size = 2.5, color = "#009E73"),
      list(size = 3.5, color = "#F0E442")
    ),
    xlab = "Continent",
    ylab = "Life expectancy",
    title = list(
      "Year: 1967",
      "Year: 1987",
      "Year: 2007"
    ),
    type = list("r", "bf", "np"),
    pairwise.display = list("s", "ns", "all"),
    p.adjust.method = list("hommel", "bonferroni", "BH"),
    conf.level = list(0.99, 0.95, 0.90),
    k = list(1, 2, 3),
    effsize.type = list(
      NULL,
      "partial_omega",
      "partial_eta"
    ),
    plot.type = list("box", "boxviolin", "violin"),
    package = list("nord", "ochRe", "awtools"),
    palette = list("aurora", "parliament", "bpalette"),
    ggtheme = list(
      ggthemes::theme_stata(),
      ggplot2::theme_classic(),
      ggthemes::theme_fivethirtyeight()
    )
  ),
  .f = ggbetweenstats
)
```

The final step is to pass the `plot_list` object we just created to the
`combine_plots` function. While each of the 3 plots already has labeling
information `combine_plots` gives us an opportunity to add additional details to
the merged plots and specify the layout in rows and columns.

```{r purrr_ggbetweenstats3, fig.height = 18, fig.width = 7}
## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(title = "Changes in life expectancy across continents (1967-2007)"),
  plotgrid.args = list(ncol = 1)
)
```

## `ggwithinstats`

We will be using simulated data from then Attention Network Test provided in ANT
dataset in `ez` package.

```{r purrr_ggwithinstats, fig.height = 20, fig.width = 8}
## for reproducibility
set.seed(123)
library(ggstatsplot)
library(ez)
data("ANT") ## loading data from `ez` package

## let's split the dataframe and create a list by years of interest
cue_list <- ANT %>% split(f = .$cue, drop = TRUE)

## checking the length of the list and the names of each element
length(cue_list)

## creating a list of plots by applying the same function for elements of the list
plot_list <- purrr::pmap(
  .l = list(
    data = cue_list,
    x = "flank",
    y = "rt",
    outlier.tagging = TRUE,
    outlier.label = "group",
    outlier.coef = list(2, 2, 2.5, 3),
    outlier.label.args = list(
      list(size = 3, color = "#56B4E9"),
      list(size = 2.5, color = "#009E73"),
      list(size = 4, color = "#F0E442"),
      list(size = 2, color = "red")
    ),
    xlab = "Flank",
    ylab = "Response time",
    title = list(
      "Cue: None",
      "Cue: Center",
      "Cue: Double",
      "Cue: Spatial"
    ),
    type = list("p", "r", "bf", "np"),
    pairwise.display = list("ns", "s", "ns", "all"),
    p.adjust.method = list("fdr", "hommel", "bonferroni", "BH"),
    conf.level = list(0.99, 0.99, 0.95, 0.90),
    k = list(3, 2, 2, 3),
    effsize.type = list(
      "omega",
      "eta",
      "partial_omega",
      "partial_eta"
    ),
    package = list("ggsci", "palettetown", "palettetown", "wesanderson"),
    palette = list("lanonc_lancet", "venomoth", "blastoise", "GrandBudapest1"),
    ggtheme = list(
      ggplot2::theme_linedraw(),
      hrbrthemes::theme_ft_rc(),
      ggthemes::theme_solarized(),
      ggthemes::theme_gdocs()
    )
  ),
  .f = ggwithinstats
)

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(title = "Response times across flank conditions for each type of cue"),
  plotgrid.args = list(ncol = 1)
)
```

## `ggscatterstats`

For the next example lets use the same methodology on different data and using
`ggscatterstats` to produce scatterplots combined with marginal
histograms/boxplots/density plots with statistical details added as a subtitle.
For data we'll use `movies_long` which is from IMDB and part of the
`{ggstatsplot}` package. Since it's a large dataset with some relatively small
categories like **NC-17** we'll sample only one quarter of the data and
completely drop NC-17 using `dplyr`. 

This time we'll put all the code in one block-

```{r purrr_ggscatterstats, fig.height = 14, fig.width = 7}
## for reproducibility
set.seed(123)

mpaa_list <- movies_long %>%
  dplyr::filter(mpaa != "NC-17") %>%
  dplyr::sample_frac(size = 0.25) %>%
  split(f = .$mpaa, drop = TRUE)

## creating a list of plots
plot_list <- purrr::pmap(
  .l = list(
    data = mpaa_list,
    x = "budget",
    y = "rating",
    xlab = "Budget (in millions of US dollars)",
    ylab = "Rating on IMDB",
    title = list(
      "MPAA Rating: PG",
      "MPAA Rating: PG-13",
      "MPAA Rating: R"
    ),
    label.var = list("title"),
    ## note that you need to quote the expressions
    label.expression = list(
      quote(rating > 7.5 & budget < 100),
      quote(rating > 8 & budget < 50),
      quote(rating > 8 & budget < 10)
    ),
    type = list("r", "np", "bf"),
    xfill = list("#009E73", "#999999", "#0072B2"),
    yfill = list("#CC79A7", "#F0E442", "#D55E00"),
    ggtheme = list(
      ggthemes::theme_tufte(),
      ggplot2::theme_classic(),
      ggplot2::theme_light()
    )
  ),
  .f = ggscatterstats
)

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(
    title = "Relationship between movie budget and IMDB rating",
    caption = "Source: www.imdb.com"
  ),
  plotgrid.args = list(ncol = 1)
)
```

The remainder of the examples vary in content but follow the exact same
methodology as the earlier examples.

## `ggcorrmat`

```{r purrr_ggcorrmat, fig.height = 10, fig.width = 10}
## for reproducibility
set.seed(123)

## splitting the dataframe by cut and creating a list
## let's leave out "fair" cut
## also, to make this fast, let's only use 5% of the sample
cut_list <- ggplot2::diamonds %>%
  dplyr::sample_frac(size = 0.05) %>%
  dplyr::filter(cut != "Fair") %>%
  split(f = .$cut, drop = TRUE)

## checking the length and names of each element
length(cut_list)
names(cut_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <- purrr::pmap(
  .l = list(
    data = cut_list,
    cor.vars = list(c("carat", "depth", "table", "price")),
    type = list("pearson", "np", "robust", "bf"),
    partial = list(TRUE, FALSE, TRUE, FALSE),
    title = list("Cut: Good", "Cut: Very Good", "Cut: Premium", "Cut: Ideal"),
    p.adjust.method = list("hommel", "fdr", "BY", "hochberg"),
    lab.size = 3.5,
    colors = list(
      c("#56B4E9", "white", "#999999"),
      c("#CC79A7", "white", "#F0E442"),
      c("#56B4E9", "white", "#D55E00"),
      c("#999999", "white", "#0072B2")
    ),
    ggtheme = list(
      ggplot2::theme_linedraw(),
      ggplot2::theme_classic(),
      ggthemes::theme_fivethirtyeight(),
      ggthemes::theme_tufte()
    )
  ),
  .f = ggcorrmat
)

## combining all individual plots from the list into a single plot using
## `combine_plots` function
combine_plots(
  plotlist = plot_list,
  guides = "keep",
  annotation.args = list(
    title = "Relationship between diamond attributes and price across cut",
    caption = "Dataset: Diamonds from ggplot2 package"
  ),
  plotgrid.args = list(nrow = 2)
)
```

## `gghistostats`

```{r purrr_gghistostats, fig.height = 14, fig.width = 8}
## for reproducibility
set.seed(123)

## let's split the dataframe and create a list by continent
## let's leave out Oceania because it has just two data points
continent_list <-
  gapminder::gapminder %>%
  dplyr::filter(year == 2007, continent != "Oceania") %>%
  split(f = .$continent, drop = TRUE)

## checking the length and names of each element
length(continent_list)
names(continent_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <-
  purrr::pmap(
    .l = list(
      data = continent_list,
      x = "lifeExp",
      xlab = "Life expectancy",
      test.value = list(35.6, 58.4, 41.6, 64.7),
      type = list("p", "np", "r", "bf"),
      bf.message = list(TRUE, FALSE, FALSE, FALSE),
      title = list(
        "Continent: Africa",
        "Continent: Americas",
        "Continent: Asia",
        "Continent: Europe"
      ),
      effsize.type = list("d", "d", "g", "g"),
      normal.curve = list(TRUE, FALSE, FALSE, TRUE),
      ggtheme = list(
        ggplot2::theme_classic(),
        hrbrthemes::theme_ipsum_tw(),
        ggplot2::theme_minimal(),
        hrbrthemes::theme_modern_rc()
      )
    ),
    .f = gghistostats
  )

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(
    title = "Improvement in life expectancy worldwide since 1950",
    caption = "Note: black line - 1950; blue line - 2007"
  ),
  plotgrid.args = list(nrow = 4)
)
```  

## `ggdotplotstats`

```{r purrr_ggdotplotstats, fig.height = 16, fig.width = 8}
## for reproducibility
set.seed(123)
library(ggthemes)
library(hrbrthemes)

## let's split the dataframe and create a list by continent
## let's leave out Oceania because it has just two data points
continent_list <-
  gapminder::gapminder %>%
  dplyr::filter(continent != "Oceania") %>%
  split(f = .$continent, drop = TRUE)

## checking the length and names of each element
length(continent_list)
names(continent_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <-
  purrr::pmap(
    .l = list(
      data = continent_list,
      x = "gdpPercap",
      y = "year",
      xlab = "GDP per capita (US$, inflation-adjusted)",
      test.value = list(2500, 9000, 9500, 10000),
      type = list("p", "np", "r", "bf"),
      title = list(
        "Continent: Africa",
        "Continent: Americas",
        "Continent: Asia",
        "Continent: Europe"
      ),
      effsize.type = list("d", "d", "g", "g"),
      centrality.line.args = list(
        list(color = "red"),
        list(color = "#0072B2"),
        list(color = "#D55E00"),
        list(color = "#CC79A7")
      ),
      ggtheme = list(
        ggplot2::theme_minimal(base_family = "serif"),
        ggthemes::theme_tufte(),
        hrbrthemes::theme_ipsum_rc(axis_title_size = 10),
        ggthemes::theme_hc(bgcolor = "darkunica")
      )
    ),
    .f = ggdotplotstats
  )

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(title = "Improvement in GDP per capita from 1952-2007"),
  plotgrid.args = list(nrow = 4),
  guides = "keep"
)
```  

## `ggpiestats`

```{r purrr_ggpiestats, fig.height = 22, fig.width = 9}
## for reproducibility
set.seed(123)

## let's split the dataframe and create a list by passenger class
class_list <- Titanic_full %>% split(f = .$Class, drop = TRUE)

## checking the length and names of each element
length(class_list)
names(class_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <-
  purrr::pmap(
    .l = list(
      data = class_list,
      x = "Survived",
      y = "Sex",
      label = list("both", "count", "percentage", "both"),
      title = list(
        "Passenger class: 1st",
        "Passenger class: 2nd",
        "Passenger class: 3rd",
        "Passenger class: Crew"
      ),
      caption = list(
        "Total: 319, Died: 120, Survived: 199, % Survived: 62%",
        "Total: 272, Died: 155, Survived: 117, % Survived: 43%",
        "Total: 709, Died: 537, Survived: 172, % Survived: 25%",
        "Data not available for crew passengers"
      ),
      package = list("RColorBrewer", "ghibli", "palettetown", "yarrr"),
      palette = list("Accent", "MarnieMedium1", "pikachu", "nemo"),
      ggtheme = list(
        ggplot2::theme_grey(),
        ggplot2::theme_bw(),
        ggthemes::theme_tufte(),
        ggthemes::theme_economist()
      ),
      proportion.test = list(TRUE, FALSE, TRUE, FALSE),
      type = list("p", "p", "bf", "p")
    ),
    .f = ggpiestats
  )

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(title = "Survival in Titanic disaster by gender for all passenger classes"),
  plotgrid.args = list(ncol = 1),
  guides = "keep"
)
``` 

## `ggbarstats`

```{r purrr_ggbarstats, fig.height = 24, fig.width = 6}
## for reproducibility
set.seed(123)

## let's split the dataframe and create a list by passenger class
class_list <- Titanic_full %>% split(f = .$Class, drop = TRUE)

## checking the length and names of each element
length(class_list)
names(class_list)

## running function on every element of this list note that if you want the same
## value for a given argument across all elements of the list, you need to
## specify it just once
plot_list <-
  purrr::pmap(
    .l = list(
      data = class_list,
      x = "Survived",
      y = "Sex",
      type = "bayes",
      label = list("both", "count", "percentage", "both"),
      title = list(
        "Passenger class: 1st",
        "Passenger class: 2nd",
        "Passenger class: 3rd",
        "Passenger class: Crew"
      ),
      caption = list(
        "Total: 319, Died: 120, Survived: 199, % Survived: 62%",
        "Total: 272, Died: 155, Survived: 117, % Survived: 43%",
        "Total: 709, Died: 537, Survived: 172, % Survived: 25%",
        "Data not available for crew passengers"
      ),
      package = list("RColorBrewer", "ghibli", "palettetown", "yarrr"),
      palette = list("Accent", "MarnieMedium1", "pikachu", "nemo"),
      ggtheme = list(
        ggplot2::theme_grey(),
        ggplot2::theme_bw(),
        ggthemes::theme_tufte(),
        ggthemes::theme_economist()
      )
    ),
    .f = ggbarstats
  )

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(
    title = "Survival in Titanic disaster by gender for all passenger classes",
    caption = "Asterisks denote results from proportion tests: \n***: p < 0.001, ns: non-significant"
  ),
  plotgrid.args = list(ncol = 1),
  guides = "keep"
)
``` 

## `grouped_` variants

Note that although all the above examples were written with the non-grouped
variants of functions, the same rule holds true for the `grouped_` variants of
all the above functions.

For example, if we want to use the `grouped_gghistostats` across three different
datasets, you can use `purrr::pmap()` function. For the sake of brevity, the
plots are not displayed here, but you can run the following code and check the
individual `grouped_` plots (e.g., `plotlist[[1]]`).

```{r purrr_grouped}
## create a list of plots
plotlist <-
  purrr::pmap(
    .l = list(
      data = list(mtcars, iris, ToothGrowth),
      x = alist(wt, Sepal.Length, len),
      results.subtitle = list(FALSE),
      grouping.var = alist(am, Species, supp)
    ),
    .f = grouped_gghistostats
  )

## given that we had three different datasets, we expect a list of length 3
## (each of which contains a `grouped_` plot)
length(plotlist)
```

## Repeating function execution across multiple columns in a dataframe

```{r purrr_columns, fig.width=12}
## setup
set.seed(123)
library(ggstatsplot)
library(patchwork)

## running the same analysis on two different columns (creates a list of plots)
plotlist <-
  purrr::pmap(
    .l = list(
      data = list(movies_long),
      x = "mpaa",
      y = list("rating", "length"),
      title = list("IMDB score by MPAA rating", "Movie length by MPAA rating")
    ),
    .f = ggbetweenstats
  )

## combine plots using `patchwork`
plotlist[[1]] + plotlist[[2]]
```

## Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "Dependencies"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_document:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    toc_depth: 3
    warning: FALSE
    message: FALSE
---

```{r setup, include = FALSE}
options(
  tibble.width = Inf,
  pillar.bold = TRUE,
  pillar.neg = TRUE,
  pillar.subtle_num = TRUE,
  pillar.min_chars = Inf
)

knitr::opts_chunk$set(
  collapse = TRUE,
  dpi = 300,
  out.width = "100%",
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  error = TRUE
)
```

<!-- # Recursive dependencies -->

<!-- ```{r deplist} -->
<!-- sort(tools::package_dependencies("ggstatsplot", recursive = TRUE)[[1]]) -->
<!-- ``` -->

# Dependency graph

```{r deepdep, fig.height=12, fig.width=12}
library(deepdep)

plot_dependencies("ggstatsplot", depth = 3)
```

# Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`: <https://github.com/IndrajeetPatil/ggstatsplot/issues>

---
title: "session information"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{session information}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  dpi = 300,
  out.width = "100%",
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)
```

# Session information

This vignette provides exhaustive details about the session (attached packages,
their versions, date, source, etc.) in which the vignettes were rendered and the
website was built. This is helpful information for reproducibility. If you see
any divergences between function behavior as described on the website and what
you see on your system, this information can come in handy for diagnosing the
source of those divergences.

This session information was generated on `r Sys.Date()`:

```{r session_info}
# setup for the session
options(width = 200)

# for pretty printing of tibbles
options(tibble.width = Inf, pillar.bold = TRUE, pillar.subtle_num = TRUE)

library(ggstatsplot)

# session information
sessioninfo::session_info(include_base = TRUE)
```

# Suggestions

If you find any bugs or have any suggestions/remarks, please file an issue on
`GitHub`:
<https://github.com/IndrajeetPatil/ggstatsplot/issues>
---
title: "effect size interpretation"
author: "Indrajeet Patil"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 6
    fig.align: 'center'
    fig.asp: 0.618
    dpi: 300
    toc: true
    warning: FALSE
    message: FALSE
vignette: >
  %\VignetteIndexEntry{effect size interpretation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE, warning = FALSE, message = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Since all of these details are calculated in `{statsExpressions}` package, this
vignette has now moved to that package. You can read it here:
<https://indrajeetpatil.github.io/statsExpressions/articles/stats_details.html>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggpiestats.R
\name{ggpiestats}
\alias{ggpiestats}
\title{Pie charts with statistical tests}
\usage{
ggpiestats(
  data,
  x,
  y = NULL,
  counts = NULL,
  type = "parametric",
  paired = FALSE,
  results.subtitle = TRUE,
  label = "percentage",
  label.args = list(direction = "both"),
  label.repel = FALSE,
  k = 2L,
  proportion.test = results.subtitle,
  perc.k = 0L,
  bf.message = TRUE,
  ratio = NULL,
  conf.level = 0.95,
  sampling.plan = "indepMulti",
  fixed.margin = "rows",
  prior.concentration = 1,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  legend.title = NULL,
  ggtheme = ggstatsplot::theme_ggstatsplot(),
  package = "RColorBrewer",
  palette = "Dark2",
  ggplot.component = NULL,
  output = "plot",
  ...
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{x}{The variable to use as the \strong{rows} in the contingency table. Please
note that if there are empty factor levels in your variable, they will be
dropped.}

\item{y}{The variable to use as the \strong{columns} in the contingency table.
Please note that if there are empty factor levels in your variable, they
will be dropped. Default is \code{NULL}. If \code{NULL}, one-sample proportion test
(a goodness of fit test) will be run for the \code{x} variable. Otherwise an
appropriate association test will be run. This argument can not be \code{NULL}
for \code{ggbarstats} function.}

\item{counts}{A string naming a variable in data containing counts, or \code{NULL}
if each row represents a single observation.}

\item{type}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}

\item{paired}{Logical indicating whether data came from a within-subjects or
repeated measures design study (Default: \code{FALSE}). If \code{TRUE}, McNemar's
test expression will be returned. If \code{FALSE}, Pearson's chi-square test will
be returned.}

\item{results.subtitle}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}

\item{label}{Character decides what information needs to be displayed
on the label in each pie slice. Possible options are \code{"percentage"}
(default), \code{"counts"}, \code{"both"}.}

\item{label.args}{Additional aesthetic arguments that will be passed to
\code{geom_label}.}

\item{label.repel}{Whether labels should be repelled using \code{ggrepel} package.
This can be helpful in case the labels are overlapping.}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{proportion.test}{Decides whether proportion test for \code{x} variable is to
be carried out for each level of \code{y}. Defaults to \code{results.subtitle}. In
\code{ggbarstats}, only \emph{p}-values from this test will be displayed.}

\item{perc.k}{Numeric that decides number of decimal places for percentage
labels (Default: \code{0L}).}

\item{bf.message}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}

\item{ratio}{A vector of proportions: the expected proportions for the
proportion test (should sum to 1). Default is \code{NULL}, which means the null
is equal theoretical proportions across the levels of the nominal variable.
This means if there are two levels this will be \code{ratio = c(0.5,0.5)} or if
there are four levels this will be \code{ratio = c(0.25,0.25,0.25,0.25)}, etc.}

\item{conf.level}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}

\item{sampling.plan}{Character describing the sampling plan. Possible options
are \code{"indepMulti"} (independent multinomial; default), \code{"poisson"},
\code{"jointMulti"} (joint multinomial), \code{"hypergeom"} (hypergeometric). For
more, see \code{?BayesFactor::contingencyTableBF()}.}

\item{fixed.margin}{For the independent multinomial sampling plan, which
margin is fixed (\code{"rows"} or \code{"cols"}). Defaults to \code{"rows"}.}

\item{prior.concentration}{Specifies the prior concentration parameter, set
to \code{1} by default. It indexes the expected deviation from the null
hypothesis under the alternative, and corresponds to Gunel and Dickey's
(1974) \code{"a"} parameter.}

\item{title}{The text for the plot title.}

\item{subtitle}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}

\item{caption}{The text for the plot caption.}

\item{legend.title}{Title text for the legend.}

\item{ggtheme}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}

\item{package}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{palette}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{ggplot.component}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{...}{Currently ignored.}
}
\description{
Pie charts for categorical data with statistical details included in the plot
as a subtitle.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggpiestats.html}
}
\examples{
\donttest{
# for reproducibility
set.seed(123)
library(ggstatsplot)

# one sample goodness of fit proportion test
ggpiestats(mtcars, x = vs)

# association test (or contingency table analysis)
ggpiestats(mtcars, x = vs, y = cyl)
}
}
\seealso{
\code{\link{grouped_ggpiestats}}, \code{\link{ggbarstats}},
\code{\link{grouped_ggbarstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggscatterstats.R
\name{ggscatterstats}
\alias{ggscatterstats}
\title{Scatterplot with marginal distributions and statistical results}
\usage{
ggscatterstats(
  data,
  x,
  y,
  type = "parametric",
  conf.level = 0.95,
  bf.prior = 0.707,
  bf.message = TRUE,
  tr = 0.2,
  k = 2L,
  results.subtitle = TRUE,
  label.var = NULL,
  label.expression = NULL,
  marginal = TRUE,
  xfill = "#009E73",
  yfill = "#D55E00",
  point.args = list(size = 3, alpha = 0.4, stroke = 0, na.rm = TRUE),
  point.width.jitter = 0,
  point.height.jitter = 0,
  point.label.args = list(size = 3, max.overlaps = 1e+06),
  smooth.line.args = list(size = 1.5, color = "blue", method = "lm", formula = y ~ x,
    na.rm = TRUE),
  xsidehistogram.args = list(fill = xfill, color = "black", na.rm = TRUE),
  ysidehistogram.args = list(fill = yfill, color = "black", na.rm = TRUE),
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  ggtheme = ggstatsplot::theme_ggstatsplot(),
  ggplot.component = NULL,
  output = "plot",
  ...
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{x}{The column in \code{data} containing the explanatory variable to be
plotted on the \code{x}-axis.}

\item{y}{The column in \code{data} containing the response (outcome) variable to
be plotted on the \code{y}-axis.}

\item{type}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}

\item{conf.level}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}

\item{bf.prior}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors and posterior estimates. In
addition to numeric arguments, several named values are also recognized:
\code{"medium"}, \code{"wide"}, and \code{"ultrawide"}, corresponding to \emph{r} scale values
of 1/2, sqrt(2)/2, and 1, respectively. In case of an ANOVA, this value
corresponds to scale for fixed effects.}

\item{bf.message}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}

\item{tr}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{results.subtitle}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}

\item{label.var}{Variable to use for points labels entered as a symbol (e.g.
\code{var1}).}

\item{label.expression}{An expression evaluating to a logical vector that
determines the subset of data points to label (e.g. \code{y < 4 & z < 20}).
While using this argument with \code{purrr::pmap}, you will have to provide a
quoted expression  (e.g. \code{quote(y < 4 & z < 20)}).}

\item{marginal}{Decides whether marginal distributions will be plotted on
axes using \code{ggside} functions. The default is \code{TRUE}. The package
\code{ggside} must already be installed by the user.}

\item{xfill, yfill}{Character describing color fill for \code{x} and \code{y} axes
marginal distributions (default: \code{"#009E73"} (for \code{x}) and \code{"#D55E00"} (for
\code{y})). Note that the defaults are colorblind-friendly.}

\item{point.args}{A list of additional aesthetic arguments to be passed
to \code{geom_point} geom used to display the raw data points.}

\item{point.width.jitter, point.height.jitter}{Degree of jitter in \code{x} and \code{y}
direction, respectively. Defaults to \code{0} (0\%) of the resolution of the
data. Note that the jitter should not be specified in the \code{point.args}
because this information will be passed to two different \code{geom}s: one
displaying the \strong{points} and the other displaying the *\strong{labels} for
these points.}

\item{point.label.args}{A list of additional aesthetic arguments to be passed
to \code{ggrepel::geom_label_repel} geom used to display the labels.}

\item{smooth.line.args}{A list of additional aesthetic arguments to be passed
to \code{geom_smooth} geom used to display the regression line.}

\item{xsidehistogram.args, ysidehistogram.args}{A list of arguments passed to
respective \code{geom_}s from \code{ggside} package to change the marginal
distribution histograms plots.}

\item{xlab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{ylab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{title}{The text for the plot title.}

\item{subtitle}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}

\item{caption}{The text for the plot caption.}

\item{ggtheme}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}

\item{ggplot.component}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{...}{Currently ignored.}
}
\description{
Scatterplots from \code{{ggplot2}} combined with marginal densigram (density +
histogram) plots with statistical details.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggscatterstats.html}
}
\note{
The plot uses \code{ggrepel::geom_label_repel} to attempt to keep labels
from over-lapping to the largest degree possible.  As a consequence plot
times will slow down massively (and the plot file will grow in size) if you
have a lot of labels that overlap.
}
\examples{
# to get reproducible results from bootstrapping
set.seed(123)
library(ggstatsplot)
library(dplyr, warn.conflicts = FALSE)

# creating dataframe with rownames converted to a new column
mtcars_new <- as_tibble(mtcars, rownames = "car")

# simple function call with the defaults
if (require("ggside")) {
  ggscatterstats(
    data = mtcars_new,
    x = wt,
    y = mpg,
    label.var = car,
    label.expression = wt < 4 & mpg < 20
  ) + # making further customization with `{ggplot2}` functions
    geom_rug(sides = "b")
}
}
\seealso{
\code{\link{grouped_ggscatterstats}}, \code{\link{ggcorrmat}},
\code{\link{grouped_ggcorrmat}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggbetweenstats.R
\name{grouped_ggbetweenstats}
\alias{grouped_ggbetweenstats}
\title{Violin plots for group or condition comparisons in between-subjects
designs repeated across all levels of a grouping variable.}
\usage{
grouped_ggbetweenstats(
  data,
  ...,
  grouping.var,
  output = "plot",
  plotgrid.args = list(),
  annotation.args = list()
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{...}{
  Arguments passed on to \code{\link[=ggbetweenstats]{ggbetweenstats}}
  \describe{
    \item{\code{plot.type}}{Character describing the \emph{type} of plot. Currently supported
plots are \code{"box"} (for only boxplots), \code{"violin"} (for only violin plots),
and \code{"boxviolin"} (for a combination of box and violin plots; default).}
    \item{\code{xlab}}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}
    \item{\code{ylab}}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}
    \item{\code{pairwise.comparisons}}{Logical that decides whether pairwise comparisons
are to be displayed (default: \code{TRUE}). Please note that only
\strong{significant} comparisons will be shown by default. To change this
behavior, select appropriate option with \code{pairwise.display} argument. The
pairwise comparison dataframes are prepared using the
\code{pairwise_comparisons} function. For more details
about pairwise comparisons, see the documentation for that function.}
    \item{\code{p.adjust.method}}{Adjustment method for \emph{p}-values for multiple
comparisons. Possible methods are: \code{"holm"} (default), \code{"hochberg"},
\code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.}
    \item{\code{pairwise.display}}{Decides \emph{which} pairwise comparisons to display.
Available options are:
\itemize{
\item \code{"significant"} (abbreviation accepted: \code{"s"})
\item \code{"non-significant"} (abbreviation accepted: \code{"ns"})
\item \code{"all"}
}

You can use this argument to make sure that your plot is not uber-cluttered
when you have multiple groups being compared and scores of pairwise
comparisons being displayed.}
    \item{\code{bf.prior}}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors.}
    \item{\code{bf.message}}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}
    \item{\code{results.subtitle}}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}
    \item{\code{subtitle}}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}
    \item{\code{caption}}{The text for the plot caption.}
    \item{\code{outlier.color}}{Default aesthetics for outliers (Default: \code{"black"}).}
    \item{\code{outlier.tagging}}{Decides whether outliers should be tagged (Default:
\code{FALSE}).}
    \item{\code{outlier.label}}{Label to put on the outliers that have been tagged. This
\strong{can't} be the same as \code{x} argument.}
    \item{\code{outlier.shape}}{Hiding the outliers can be achieved by setting
\code{outlier.shape = NA}. Importantly, this does not remove the outliers,
it only hides them, so the range calculated for the \code{y}-axis will be
the same with outliers shown and outliers hidden.}
    \item{\code{outlier.label.args}}{A list of additional aesthetic arguments to be
passed to \code{ggrepel::geom_label_repel} for outlier label plotting.}
    \item{\code{outlier.coef}}{Coefficient for outlier detection using Tukey's method.
With Tukey's method, outliers are below (1st Quartile) or above (3rd
Quartile) \code{outlier.coef} times the Inter-Quartile Range (IQR) (Default:
\code{1.5}).}
    \item{\code{centrality.plotting}}{Logical that decides whether centrality tendency
measure is to be displayed as a point with a label (Default: \code{TRUE}).
Function decides which central tendency measure to show depending on the
\code{type} argument.
\itemize{
\item \strong{mean} for parametric statistics
\item \strong{median} for non-parametric statistics
\item \strong{trimmed mean} for robust statistics
\item \strong{MAP estimator} for Bayesian statistics
}

If you want default centrality parameter, you can specify this using
\code{centrality.type} argument.}
    \item{\code{centrality.type}}{Decides which centrality parameter is to be displayed.
The default is to choose the same as \code{type} argument. You can specify this
to be:
\itemize{
\item \code{"parameteric"} (for \strong{mean})
\item \code{"nonparametric"} (for \strong{median})
\item \code{robust} (for \strong{trimmed mean})
\item \code{bayes} (for \strong{MAP estimator})
}

Just as \code{type} argument, abbreviations are also accepted.}
    \item{\code{point.args}}{A list of additional aesthetic arguments to be passed to
the \code{geom_point} displaying the raw data.}
    \item{\code{violin.args}}{A list of additional aesthetic arguments to be passed to
the \code{geom_violin}.}
    \item{\code{ggplot.component}}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}
    \item{\code{package}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{palette}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{centrality.point.args}}{A list of additional aesthetic
arguments to be passed to \code{geom_point} and
\code{ggrepel::geom_label_repel} geoms, which are involved in mean plotting.}
    \item{\code{centrality.label.args}}{A list of additional aesthetic
arguments to be passed to \code{geom_point} and
\code{ggrepel::geom_label_repel} geoms, which are involved in mean plotting.}
    \item{\code{ggsignif.args}}{A list of additional aesthetic
arguments to be passed to \code{ggsignif::geom_signif}.}
    \item{\code{ggtheme}}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}
    \item{\code{x}}{The grouping (or independent) variable from the dataframe \code{data}. In
case of a repeated measures or within-subjects design, if \code{subject.id}
argument is not available or not explicitly specified, the function assumes
that the data has already been sorted by such an id by the user and creates
an internal identifier. So if your data is \strong{not} sorted, the results
\emph{can} be inaccurate when there are more than two levels in \code{x} and there
are \code{NA}s present. The data is expected to be sorted by user in
subject-1,subject-2, ..., pattern.}
    \item{\code{y}}{The response (or outcome or dependent) variable from the
dataframe \code{data}.}
    \item{\code{type}}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}
    \item{\code{effsize.type}}{Type of effect size needed for \emph{parametric} tests. The
argument can be \code{"eta"} (partial eta-squared) or \code{"omega"} (partial
omega-squared).}
    \item{\code{k}}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}
    \item{\code{var.equal}}{a logical variable indicating whether to treat the
    two variances as being equal. If \code{TRUE} then the pooled
    variance is used to estimate the variance otherwise the Welch
    (or Satterthwaite) approximation to the degrees of freedom is used.}
    \item{\code{conf.level}}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}
    \item{\code{nboot}}{Number of bootstrap samples for computing confidence interval
for the effect size (Default: \code{100L}).}
    \item{\code{tr}}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}
  }}

\item{grouping.var}{A single grouping variable.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{plotgrid.args}{A \code{list} of additional arguments passed to
\code{patchwork::wrap_plots}, except for \code{guides} argument which is already
separately specified here.}

\item{annotation.args}{A \code{list} of additional arguments passed to
\code{patchwork::plot_annotation}.}
}
\description{
Helper function for \code{ggstatsplot::ggbetweenstats} to apply this function
across multiple levels of a given factor and combining the resulting plots
using \code{ggstatsplot::combine_plots}.
}
\examples{
\donttest{
if (require("PMCMRplus")) {
  # to get reproducible results from bootstrapping
  set.seed(123)
  library(ggstatsplot)
  library(dplyr, warn.conflicts = FALSE)
  library(ggplot2)

  # the most basic function call
  grouped_ggbetweenstats(
    data = filter(ggplot2::mpg, drv != "4"),
    x = year,
    y = hwy,
    grouping.var = drv
  )

  # modifying individual plots using `ggplot.component` argument
  grouped_ggbetweenstats(
    data = filter(
      movies_long,
      genre \%in\% c("Action", "Comedy"),
      mpaa \%in\% c("R", "PG")
    ),
    x = genre,
    y = rating,
    grouping.var = mpaa,
    ggplot.component = scale_y_continuous(
      breaks = seq(1, 9, 1),
      limits = (c(1, 9))
    )
  )
}
}
}
\seealso{
\code{\link{ggbetweenstats}}, \code{\link{ggwithinstats}},
\code{\link{grouped_ggwithinstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggcorrmat.R
\name{grouped_ggcorrmat}
\alias{grouped_ggcorrmat}
\title{Visualization of a correlalogram (or correlation matrix) for all
levels of a grouping variable}
\usage{
grouped_ggcorrmat(
  data,
  ...,
  grouping.var,
  output = "plot",
  plotgrid.args = list(),
  annotation.args = list()
)
}
\arguments{
\item{data}{Dataframe from which variables specified are preferentially to be
taken.}

\item{...}{
  Arguments passed on to \code{\link[=ggcorrmat]{ggcorrmat}}
  \describe{
    \item{\code{cor.vars}}{List of variables for which the correlation matrix is to be
computed and visualized. If \code{NULL} (default), all numeric variables from
\code{data} will be used.}
    \item{\code{cor.vars.names}}{Optional list of names to be used for \code{cor.vars}. The
names should be entered in the same order.}
    \item{\code{partial}}{Can be \code{TRUE} for partial correlations. For Bayesian partial
correlations, "full" instead of pseudo-Bayesian partial correlations (i.e.,
Bayesian correlation based on frequentist partialization) are returned.}
    \item{\code{matrix.type}}{Character, \code{"upper"} (default), \code{"lower"}, or \code{"full"},
display full matrix, lower triangular or upper triangular matrix.}
    \item{\code{sig.level}}{Significance level (Default: \code{0.05}). If the \emph{p}-value in
\emph{p}-value matrix is bigger than \code{sig.level}, then the corresponding
correlation coefficient is regarded as insignificant and flagged as such in
the plot. Relevant only when \code{output = "plot"}.}
    \item{\code{colors}}{A vector of 3 colors for low, mid, and high correlation values.
If set to \code{NULL}, manual specification of colors will be turned off and 3
colors from the specified \code{palette} from \code{package} will be selected.}
    \item{\code{pch}}{Decides the point shape to be used for insignificant correlation
coefficients (only valid when \code{insig = "pch"}). Default: \code{pch = "cross"}.}
    \item{\code{ggcorrplot.args}}{A list of additional (mostly aesthetic) arguments that
will be passed to \code{ggcorrplot::ggcorrplot} function. The list should avoid
any of the following arguments since they are already internally being
used: \code{corr}, \code{method}, \code{p.mat}, \code{sig.level}, \code{ggtheme}, \code{colors}, \code{lab},
\code{pch}, \code{legend.title}, \code{digits}.}
    \item{\code{type}}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}
    \item{\code{tr}}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}
    \item{\code{k}}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}
    \item{\code{conf.level}}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}
    \item{\code{bf.prior}}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors and posterior estimates. In
addition to numeric arguments, several named values are also recognized:
\code{"medium"}, \code{"wide"}, and \code{"ultrawide"}, corresponding to \emph{r} scale values
of 1/2, sqrt(2)/2, and 1, respectively. In case of an ANOVA, this value
corresponds to scale for fixed effects.}
    \item{\code{p.adjust.method}}{Adjustment method for \emph{p}-values for multiple
comparisons. Possible methods are: \code{"holm"} (default), \code{"hochberg"},
\code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.}
    \item{\code{package}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{palette}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{ggtheme}}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}
    \item{\code{ggplot.component}}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}
    \item{\code{subtitle}}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}
    \item{\code{caption}}{The text for the plot caption.}
  }}

\item{grouping.var}{A single grouping variable.}

\item{output}{Character that decides expected output from this function. If
\code{"plot"}, the visualization matrix will be returned. If \code{"dataframe"} (or
literally anything other than \code{"plot"}), a dataframe containing all details
from statistical analyses (e.g., correlation coefficients, statistic
values, \emph{p}-values, no. of observations, etc.) will be returned.}

\item{plotgrid.args}{A \code{list} of additional arguments passed to
\code{patchwork::wrap_plots}, except for \code{guides} argument which is already
separately specified here.}

\item{annotation.args}{A \code{list} of additional arguments passed to
\code{patchwork::plot_annotation}.}
}
\description{
Helper function for \code{ggstatsplot::ggcorrmat} to apply this function across
multiple levels of a given factor and combining the resulting plots using
\code{ggstatsplot::combine_plots}.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcorrmat.html}
}
\examples{
\donttest{
# for reproducibility
set.seed(123)
library(ggstatsplot)

# for plot
if (require("ggcorrplot")) {
  grouped_ggcorrmat(
    data = iris,
    grouping.var = Species,
    type = "robust",
    p.adjust.method = "holm",
    plotgrid.args = list(ncol = 1),
    annotation.args = list(tag_levels = "i")
  )
}

# for dataframe
grouped_ggcorrmat(
  data = ggplot2::msleep,
  grouping.var = vore,
  type = "bayes",
  output = "dataframe"
)
}
}
\seealso{
\code{\link{ggcorrmat}}, \code{\link{ggscatterstats}},
\code{\link{grouped_ggscatterstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggbarstats.R
\name{ggbarstats}
\alias{ggbarstats}
\title{Bar (column) charts with statistical tests}
\usage{
ggbarstats(
  data,
  x,
  y,
  counts = NULL,
  type = "parametric",
  paired = FALSE,
  results.subtitle = TRUE,
  label = "percentage",
  label.args = list(alpha = 1, fill = "white"),
  k = 2L,
  proportion.test = results.subtitle,
  perc.k = 0L,
  bf.message = TRUE,
  ratio = NULL,
  conf.level = 0.95,
  sampling.plan = "indepMulti",
  fixed.margin = "rows",
  prior.concentration = 1,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  legend.title = NULL,
  xlab = NULL,
  ylab = NULL,
  ggtheme = ggstatsplot::theme_ggstatsplot(),
  package = "RColorBrewer",
  palette = "Dark2",
  ggplot.component = NULL,
  output = "plot",
  ...
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{x}{The variable to use as the \strong{rows} in the contingency table. Please
note that if there are empty factor levels in your variable, they will be
dropped.}

\item{y}{The variable to use as the \strong{columns} in the contingency table.
Please note that if there are empty factor levels in your variable, they
will be dropped. Default is \code{NULL}. If \code{NULL}, one-sample proportion test
(a goodness of fit test) will be run for the \code{x} variable. Otherwise an
appropriate association test will be run. This argument can not be \code{NULL}
for \code{ggbarstats} function.}

\item{counts}{A string naming a variable in data containing counts, or \code{NULL}
if each row represents a single observation.}

\item{type}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}

\item{paired}{Logical indicating whether data came from a within-subjects or
repeated measures design study (Default: \code{FALSE}). If \code{TRUE}, McNemar's
test expression will be returned. If \code{FALSE}, Pearson's chi-square test will
be returned.}

\item{results.subtitle}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}

\item{label}{Character decides what information needs to be displayed
on the label in each pie slice. Possible options are \code{"percentage"}
(default), \code{"counts"}, \code{"both"}.}

\item{label.args}{Additional aesthetic arguments that will be passed to
\code{geom_label}.}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{proportion.test}{Decides whether proportion test for \code{x} variable is to
be carried out for each level of \code{y}. Defaults to \code{results.subtitle}. In
\code{ggbarstats}, only \emph{p}-values from this test will be displayed.}

\item{perc.k}{Numeric that decides number of decimal places for percentage
labels (Default: \code{0L}).}

\item{bf.message}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}

\item{ratio}{A vector of proportions: the expected proportions for the
proportion test (should sum to 1). Default is \code{NULL}, which means the null
is equal theoretical proportions across the levels of the nominal variable.
This means if there are two levels this will be \code{ratio = c(0.5,0.5)} or if
there are four levels this will be \code{ratio = c(0.25,0.25,0.25,0.25)}, etc.}

\item{conf.level}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}

\item{sampling.plan}{Character describing the sampling plan. Possible options
are \code{"indepMulti"} (independent multinomial; default), \code{"poisson"},
\code{"jointMulti"} (joint multinomial), \code{"hypergeom"} (hypergeometric). For
more, see \code{?BayesFactor::contingencyTableBF()}.}

\item{fixed.margin}{For the independent multinomial sampling plan, which
margin is fixed (\code{"rows"} or \code{"cols"}). Defaults to \code{"rows"}.}

\item{prior.concentration}{Specifies the prior concentration parameter, set
to \code{1} by default. It indexes the expected deviation from the null
hypothesis under the alternative, and corresponds to Gunel and Dickey's
(1974) \code{"a"} parameter.}

\item{title}{The text for the plot title.}

\item{subtitle}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}

\item{caption}{The text for the plot caption.}

\item{legend.title}{Title text for the legend.}

\item{xlab}{Custom text for the \code{x} axis label (Default: \code{NULL}, which
will cause the \code{x} axis label to be the \code{x} variable).}

\item{ylab}{Custom text for the \code{y} axis label (Default: \code{NULL}).}

\item{ggtheme}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}

\item{package}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{palette}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{ggplot.component}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{...}{Currently ignored.}
}
\description{
Bar charts for categorical data with statistical details included in the plot
as a subtitle.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggpiestats.html}
}
\examples{
\donttest{
# for reproducibility
set.seed(123)
library(ggstatsplot)

# association test (or contingency table analysis)
ggbarstats(mtcars, x = vs, y = cyl)
}
}
\seealso{
\code{\link{grouped_ggbarstats}}, \code{\link{ggpiestats}},
\code{\link{grouped_ggpiestats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gghistostats.R
\name{grouped_gghistostats}
\alias{grouped_gghistostats}
\title{Grouped histograms for distribution of a numeric variable}
\usage{
grouped_gghistostats(
  data,
  x,
  grouping.var,
  binwidth = NULL,
  output = "plot",
  plotgrid.args = list(),
  annotation.args = list(),
  ...
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{x}{A numeric variable from the dataframe \code{data}.}

\item{grouping.var}{A single grouping variable.}

\item{binwidth}{The width of the histogram bins. Can be specified as a
numeric value, or a function that calculates width from \code{x}. The default is
to use the \code{max(x) - min(x) / sqrt(N)}. You should always check this value
and explore multiple widths to find the best to illustrate the stories in
your data.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{plotgrid.args}{A \code{list} of additional arguments passed to
\code{patchwork::wrap_plots}, except for \code{guides} argument which is already
separately specified here.}

\item{annotation.args}{A \code{list} of additional arguments passed to
\code{patchwork::plot_annotation}.}

\item{...}{
  Arguments passed on to \code{\link[=gghistostats]{gghistostats}}
  \describe{
    \item{\code{normal.curve}}{A logical value that decides whether to super-impose a
normal curve using \code{stats::dnorm(mean(x), sd(x))}. Default is \code{FALSE}.}
    \item{\code{normal.curve.args}}{A list of additional aesthetic arguments to be
passed to the normal curve.}
    \item{\code{bin.args}}{A list of additional aesthetic arguments to be passed to the
\code{stat_bin} used to display the bins. Do not specify \code{binwidth} argument in
this list since it has already been specified using the dedicated argument.}
    \item{\code{centrality.line.args}}{A list of additional aesthetic arguments to be
passed to the \code{geom_line} used to display the lines corresponding to the
centrality parameter.}
    \item{\code{type}}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}
    \item{\code{test.value}}{A number indicating the true value of the mean (Default:
\code{0}).}
    \item{\code{bf.prior}}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors and posterior estimates. In
addition to numeric arguments, several named values are also recognized:
\code{"medium"}, \code{"wide"}, and \code{"ultrawide"}, corresponding to \emph{r} scale values
of 1/2, sqrt(2)/2, and 1, respectively. In case of an ANOVA, this value
corresponds to scale for fixed effects.}
    \item{\code{effsize.type}}{Type of effect size needed for \emph{parametric} tests. The
argument can be \code{"d"} (for Cohen's \emph{d}) or \code{"g"} (for Hedge's \emph{g}).}
    \item{\code{conf.level}}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}
    \item{\code{tr}}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}
    \item{\code{k}}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}
    \item{\code{xlab}}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}
    \item{\code{subtitle}}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}
    \item{\code{caption}}{The text for the plot caption.}
    \item{\code{bf.message}}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}
    \item{\code{ggtheme}}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}
    \item{\code{results.subtitle}}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}
    \item{\code{centrality.plotting}}{Logical that decides whether centrality tendency
measure is to be displayed as a point with a label (Default: \code{TRUE}).
Function decides which central tendency measure to show depending on the
\code{type} argument.
\itemize{
\item \strong{mean} for parametric statistics
\item \strong{median} for non-parametric statistics
\item \strong{trimmed mean} for robust statistics
\item \strong{MAP estimator} for Bayesian statistics
}

If you want default centrality parameter, you can specify this using
\code{centrality.type} argument.}
    \item{\code{centrality.type}}{Decides which centrality parameter is to be displayed.
The default is to choose the same as \code{type} argument. You can specify this
to be:
\itemize{
\item \code{"parameteric"} (for \strong{mean})
\item \code{"nonparametric"} (for \strong{median})
\item \code{robust} (for \strong{trimmed mean})
\item \code{bayes} (for \strong{MAP estimator})
}

Just as \code{type} argument, abbreviations are also accepted.}
    \item{\code{ggplot.component}}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}
  }}
}
\description{
Helper function for \code{ggstatsplot::gghistostats} to apply this function
across multiple levels of a given factor and combining the resulting plots
using \code{ggstatsplot::combine_plots}.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/gghistostats.html}
}
\examples{
\donttest{
# for reproducibility
set.seed(123)
library(ggstatsplot)

# plot
grouped_gghistostats(
  data            = iris,
  x               = Sepal.Length,
  test.value      = 5,
  grouping.var    = Species,
  plotgrid.args   = list(nrow = 1),
  annotation.args = list(tag_levels = "i"),
)
}
}
\seealso{
\code{\link{gghistostats}}, \code{\link{ggdotplotstats}},
\code{\link{grouped_ggdotplotstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bugs_long}
\alias{bugs_long}
\title{Tidy version of the "Bugs" dataset.}
\format{
A data frame with 372 rows and 6 variables
\itemize{
\item subject. Dummy identity number for each participant.
\item gender. Participant's gender (Female, Male).
\item region. Region of the world the participant was from.
\item education. Level of education.
\item condition. Condition of the experiment the participant gave rating
for (\strong{LDLF}: low freighteningness and low disgustingness; \strong{LFHD}: low
freighteningness and high disgustingness; \strong{HFHD}: high freighteningness
and low disgustingness; \strong{HFHD}: high freighteningness and high
disgustingness).
\item desire. The desire to kill an arthropod was indicated on a scale from
0 to 10.
}
}
\source{
\url{https://www.sciencedirect.com/science/article/pii/S0747563213000277}
}
\usage{
bugs_long
}
\description{
Tidy version of the "Bugs" dataset.
}
\details{
This data set, "Bugs", provides the extent to which men and women
want to kill arthropods that vary in freighteningness (low, high) and
disgustingness (low, high). Each participant rates their attitudes towards
all anthropods. Subset of the data reported by Ryan et al. (2013).
}
\examples{
dim(bugs_long)
head(bugs_long)
dplyr::glimpse(bugs_long)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gghistostats.R
\name{gghistostats}
\alias{gghistostats}
\title{Histogram for distribution of a numeric variable}
\usage{
gghistostats(
  data,
  x,
  binwidth = NULL,
  xlab = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  type = "parametric",
  test.value = 0,
  bf.prior = 0.707,
  bf.message = TRUE,
  effsize.type = "g",
  conf.level = 0.95,
  tr = 0.2,
  k = 2L,
  ggtheme = ggstatsplot::theme_ggstatsplot(),
  results.subtitle = TRUE,
  bin.args = list(color = "black", fill = "grey50", alpha = 0.7),
  centrality.plotting = TRUE,
  centrality.type = type,
  centrality.line.args = list(color = "blue", size = 1, linetype = "dashed"),
  normal.curve = FALSE,
  normal.curve.args = list(size = 2),
  ggplot.component = NULL,
  output = "plot",
  ...
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{x}{A numeric variable from the dataframe \code{data}.}

\item{binwidth}{The width of the histogram bins. Can be specified as a
numeric value, or a function that calculates width from \code{x}. The default is
to use the \code{max(x) - min(x) / sqrt(N)}. You should always check this value
and explore multiple widths to find the best to illustrate the stories in
your data.}

\item{xlab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{title}{The text for the plot title.}

\item{subtitle}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}

\item{caption}{The text for the plot caption.}

\item{type}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}

\item{test.value}{A number indicating the true value of the mean (Default:
\code{0}).}

\item{bf.prior}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors and posterior estimates. In
addition to numeric arguments, several named values are also recognized:
\code{"medium"}, \code{"wide"}, and \code{"ultrawide"}, corresponding to \emph{r} scale values
of 1/2, sqrt(2)/2, and 1, respectively. In case of an ANOVA, this value
corresponds to scale for fixed effects.}

\item{bf.message}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}

\item{effsize.type}{Type of effect size needed for \emph{parametric} tests. The
argument can be \code{"d"} (for Cohen's \emph{d}) or \code{"g"} (for Hedge's \emph{g}).}

\item{conf.level}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}

\item{tr}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{ggtheme}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}

\item{results.subtitle}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}

\item{bin.args}{A list of additional aesthetic arguments to be passed to the
\code{stat_bin} used to display the bins. Do not specify \code{binwidth} argument in
this list since it has already been specified using the dedicated argument.}

\item{centrality.plotting}{Logical that decides whether centrality tendency
measure is to be displayed as a point with a label (Default: \code{TRUE}).
Function decides which central tendency measure to show depending on the
\code{type} argument.
\itemize{
\item \strong{mean} for parametric statistics
\item \strong{median} for non-parametric statistics
\item \strong{trimmed mean} for robust statistics
\item \strong{MAP estimator} for Bayesian statistics
}

If you want default centrality parameter, you can specify this using
\code{centrality.type} argument.}

\item{centrality.type}{Decides which centrality parameter is to be displayed.
The default is to choose the same as \code{type} argument. You can specify this
to be:
\itemize{
\item \code{"parameteric"} (for \strong{mean})
\item \code{"nonparametric"} (for \strong{median})
\item \code{robust} (for \strong{trimmed mean})
\item \code{bayes} (for \strong{MAP estimator})
}

Just as \code{type} argument, abbreviations are also accepted.}

\item{centrality.line.args}{A list of additional aesthetic arguments to be
passed to the \code{geom_line} used to display the lines corresponding to the
centrality parameter.}

\item{normal.curve}{A logical value that decides whether to super-impose a
normal curve using \code{stats::dnorm(mean(x), sd(x))}. Default is \code{FALSE}.}

\item{normal.curve.args}{A list of additional aesthetic arguments to be
passed to the normal curve.}

\item{ggplot.component}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{...}{Currently ignored.}
}
\description{
Histogram with statistical details from one-sample test included in the plot
as a subtitle.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/gghistostats.html}
}
\examples{
\donttest{
# for reproducibility
set.seed(123)
library(ggstatsplot)

# using defaults, but modifying which centrality parameter is to be shown
gghistostats(
  data            = ToothGrowth,
  x               = len,
  xlab            = "Tooth length",
  centrality.type = "np"
)
}
}
\seealso{
\code{\link{grouped_gghistostats}}, \code{\link{ggdotplotstats}},
\code{\link{grouped_ggdotplotstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggbetweenstats.R
\name{ggbetweenstats}
\alias{ggbetweenstats}
\title{Box/Violin plots for between-subjects comparisons}
\usage{
ggbetweenstats(
  data,
  x,
  y,
  plot.type = "boxviolin",
  type = "parametric",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  p.adjust.method = "holm",
  effsize.type = "unbiased",
  bf.prior = 0.707,
  bf.message = TRUE,
  results.subtitle = TRUE,
  xlab = NULL,
  ylab = NULL,
  caption = NULL,
  title = NULL,
  subtitle = NULL,
  k = 2L,
  var.equal = FALSE,
  conf.level = 0.95,
  nboot = 100L,
  tr = 0.2,
  centrality.plotting = TRUE,
  centrality.type = type,
  centrality.point.args = list(size = 5, color = "darkred"),
  centrality.label.args = list(size = 3, nudge_x = 0.4, segment.linetype = 4,
    min.segment.length = 0),
  outlier.tagging = FALSE,
  outlier.label = NULL,
  outlier.coef = 1.5,
  outlier.shape = 19,
  outlier.color = "black",
  outlier.label.args = list(size = 3),
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha
    = 0.4, size = 3, stroke = 0),
  violin.args = list(width = 0.5, alpha = 0.2),
  ggsignif.args = list(textsize = 3, tip_length = 0.01),
  ggtheme = ggstatsplot::theme_ggstatsplot(),
  package = "RColorBrewer",
  palette = "Dark2",
  ggplot.component = NULL,
  output = "plot",
  ...
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{x}{The grouping (or independent) variable from the dataframe \code{data}. In
case of a repeated measures or within-subjects design, if \code{subject.id}
argument is not available or not explicitly specified, the function assumes
that the data has already been sorted by such an id by the user and creates
an internal identifier. So if your data is \strong{not} sorted, the results
\emph{can} be inaccurate when there are more than two levels in \code{x} and there
are \code{NA}s present. The data is expected to be sorted by user in
subject-1,subject-2, ..., pattern.}

\item{y}{The response (or outcome or dependent) variable from the
dataframe \code{data}.}

\item{plot.type}{Character describing the \emph{type} of plot. Currently supported
plots are \code{"box"} (for only boxplots), \code{"violin"} (for only violin plots),
and \code{"boxviolin"} (for a combination of box and violin plots; default).}

\item{type}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}

\item{pairwise.comparisons}{Logical that decides whether pairwise comparisons
are to be displayed (default: \code{TRUE}). Please note that only
\strong{significant} comparisons will be shown by default. To change this
behavior, select appropriate option with \code{pairwise.display} argument. The
pairwise comparison dataframes are prepared using the
\code{pairwise_comparisons} function. For more details
about pairwise comparisons, see the documentation for that function.}

\item{pairwise.display}{Decides \emph{which} pairwise comparisons to display.
Available options are:
\itemize{
\item \code{"significant"} (abbreviation accepted: \code{"s"})
\item \code{"non-significant"} (abbreviation accepted: \code{"ns"})
\item \code{"all"}
}

You can use this argument to make sure that your plot is not uber-cluttered
when you have multiple groups being compared and scores of pairwise
comparisons being displayed.}

\item{p.adjust.method}{Adjustment method for \emph{p}-values for multiple
comparisons. Possible methods are: \code{"holm"} (default), \code{"hochberg"},
\code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.}

\item{effsize.type}{Type of effect size needed for \emph{parametric} tests. The
argument can be \code{"eta"} (partial eta-squared) or \code{"omega"} (partial
omega-squared).}

\item{bf.prior}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors.}

\item{bf.message}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}

\item{results.subtitle}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}

\item{xlab, ylab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{caption}{The text for the plot caption.}

\item{title}{The text for the plot title.}

\item{subtitle}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{var.equal}{a logical variable indicating whether to treat the
    two variances as being equal. If \code{TRUE} then the pooled
    variance is used to estimate the variance otherwise the Welch
    (or Satterthwaite) approximation to the degrees of freedom is used.}

\item{conf.level}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}

\item{nboot}{Number of bootstrap samples for computing confidence interval
for the effect size (Default: \code{100L}).}

\item{tr}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}

\item{centrality.plotting}{Logical that decides whether centrality tendency
measure is to be displayed as a point with a label (Default: \code{TRUE}).
Function decides which central tendency measure to show depending on the
\code{type} argument.
\itemize{
\item \strong{mean} for parametric statistics
\item \strong{median} for non-parametric statistics
\item \strong{trimmed mean} for robust statistics
\item \strong{MAP estimator} for Bayesian statistics
}

If you want default centrality parameter, you can specify this using
\code{centrality.type} argument.}

\item{centrality.type}{Decides which centrality parameter is to be displayed.
The default is to choose the same as \code{type} argument. You can specify this
to be:
\itemize{
\item \code{"parameteric"} (for \strong{mean})
\item \code{"nonparametric"} (for \strong{median})
\item \code{robust} (for \strong{trimmed mean})
\item \code{bayes} (for \strong{MAP estimator})
}

Just as \code{type} argument, abbreviations are also accepted.}

\item{centrality.point.args, centrality.label.args}{A list of additional aesthetic
arguments to be passed to \code{geom_point} and
\code{ggrepel::geom_label_repel} geoms, which are involved in mean plotting.}

\item{outlier.tagging}{Decides whether outliers should be tagged (Default:
\code{FALSE}).}

\item{outlier.label}{Label to put on the outliers that have been tagged. This
\strong{can't} be the same as \code{x} argument.}

\item{outlier.coef}{Coefficient for outlier detection using Tukey's method.
With Tukey's method, outliers are below (1st Quartile) or above (3rd
Quartile) \code{outlier.coef} times the Inter-Quartile Range (IQR) (Default:
\code{1.5}).}

\item{outlier.shape}{Hiding the outliers can be achieved by setting
\code{outlier.shape = NA}. Importantly, this does not remove the outliers,
it only hides them, so the range calculated for the \code{y}-axis will be
the same with outliers shown and outliers hidden.}

\item{outlier.color}{Default aesthetics for outliers (Default: \code{"black"}).}

\item{outlier.label.args}{A list of additional aesthetic arguments to be
passed to \code{ggrepel::geom_label_repel} for outlier label plotting.}

\item{point.args}{A list of additional aesthetic arguments to be passed to
the \code{geom_point} displaying the raw data.}

\item{violin.args}{A list of additional aesthetic arguments to be passed to
the \code{geom_violin}.}

\item{ggsignif.args}{A list of additional aesthetic
arguments to be passed to \code{ggsignif::geom_signif}.}

\item{ggtheme}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}

\item{package, palette}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{ggplot.component}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{...}{Currently ignored.}
}
\description{
A combination of box and violin plots along with jittered data points for
between-subjects designs with statistical details included in the plot as a
subtitle.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggbetweenstats.html}
}
\examples{
\donttest{
if (require("PMCMRplus")) {
  # to get reproducible results from bootstrapping
  set.seed(123)
  library(ggstatsplot)

  # simple function call with the defaults
  ggbetweenstats(mtcars, am, mpg)

  # more detailed function call
  ggbetweenstats(
    data = morley,
    x = Expt,
    y = Speed,
    type = "robust",
    xlab = "The experiment number",
    ylab = "Speed-of-light measurement",
    pairwise.comparisons = TRUE,
    p.adjust.method = "fdr",
    outlier.tagging = TRUE,
    outlier.label = Run
  )
}
}
}
\seealso{
\code{\link{grouped_ggbetweenstats}}, \code{\link{ggwithinstats}},
\code{\link{grouped_ggwithinstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggdotplotstats.R
\name{grouped_ggdotplotstats}
\alias{grouped_ggdotplotstats}
\title{Grouped histograms for distribution of a labeled numeric variable}
\usage{
grouped_ggdotplotstats(
  data,
  ...,
  grouping.var,
  output = "plot",
  plotgrid.args = list(),
  annotation.args = list()
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{...}{
  Arguments passed on to \code{\link[=ggdotplotstats]{ggdotplotstats}}
  \describe{
    \item{\code{y}}{Label or grouping variable.}
    \item{\code{point.args}}{A list of additional aesthetic arguments passed to
\code{geom_point}.}
    \item{\code{x}}{A numeric variable from the dataframe \code{data}.}
    \item{\code{xlab}}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}
    \item{\code{subtitle}}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}
    \item{\code{caption}}{The text for the plot caption.}
    \item{\code{type}}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}
    \item{\code{test.value}}{A number indicating the true value of the mean (Default:
\code{0}).}
    \item{\code{bf.prior}}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors and posterior estimates. In
addition to numeric arguments, several named values are also recognized:
\code{"medium"}, \code{"wide"}, and \code{"ultrawide"}, corresponding to \emph{r} scale values
of 1/2, sqrt(2)/2, and 1, respectively. In case of an ANOVA, this value
corresponds to scale for fixed effects.}
    \item{\code{bf.message}}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}
    \item{\code{effsize.type}}{Type of effect size needed for \emph{parametric} tests. The
argument can be \code{"d"} (for Cohen's \emph{d}) or \code{"g"} (for Hedge's \emph{g}).}
    \item{\code{conf.level}}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}
    \item{\code{tr}}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}
    \item{\code{k}}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}
    \item{\code{results.subtitle}}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}
    \item{\code{centrality.plotting}}{Logical that decides whether centrality tendency
measure is to be displayed as a point with a label (Default: \code{TRUE}).
Function decides which central tendency measure to show depending on the
\code{type} argument.
\itemize{
\item \strong{mean} for parametric statistics
\item \strong{median} for non-parametric statistics
\item \strong{trimmed mean} for robust statistics
\item \strong{MAP estimator} for Bayesian statistics
}

If you want default centrality parameter, you can specify this using
\code{centrality.type} argument.}
    \item{\code{centrality.type}}{Decides which centrality parameter is to be displayed.
The default is to choose the same as \code{type} argument. You can specify this
to be:
\itemize{
\item \code{"parameteric"} (for \strong{mean})
\item \code{"nonparametric"} (for \strong{median})
\item \code{robust} (for \strong{trimmed mean})
\item \code{bayes} (for \strong{MAP estimator})
}

Just as \code{type} argument, abbreviations are also accepted.}
    \item{\code{centrality.line.args}}{A list of additional aesthetic arguments to be
passed to the \code{geom_line} used to display the lines corresponding to the
centrality parameter.}
    \item{\code{ggplot.component}}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}
    \item{\code{ggtheme}}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}
    \item{\code{ylab}}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}
  }}

\item{grouping.var}{A single grouping variable.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{plotgrid.args}{A \code{list} of additional arguments passed to
\code{patchwork::wrap_plots}, except for \code{guides} argument which is already
separately specified here.}

\item{annotation.args}{A \code{list} of additional arguments passed to
\code{patchwork::plot_annotation}.}
}
\description{
Helper function for \code{ggstatsplot::ggdotplotstats} to apply this function
across multiple levels of a given factor and combining the resulting plots
using \code{ggstatsplot::combine_plots}.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggdotplotstats.html}
}
\examples{
\donttest{
# for reproducibility
set.seed(123)
library(ggstatsplot)
library(dplyr, warn.conflicts = FALSE)

# removing factor level with very few no. of observations
df <- filter(ggplot2::mpg, cyl \%in\% c("4", "6", "8"))

# plot
grouped_ggdotplotstats(
  data         = df,
  x            = cty,
  y            = manufacturer,
  grouping.var = cyl,
  test.value   = 15.5
)
}
}
\seealso{
\code{\link{grouped_gghistostats}}, \code{\link{ggdotplotstats}},
\code{\link{gghistostats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggwithinstats.R
\name{grouped_ggwithinstats}
\alias{grouped_ggwithinstats}
\title{Violin plots for group or condition comparisons in within-subjects
designs repeated across all levels of a grouping variable.}
\usage{
grouped_ggwithinstats(
  data,
  ...,
  grouping.var,
  output = "plot",
  plotgrid.args = list(),
  annotation.args = list()
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{...}{
  Arguments passed on to \code{\link[=ggwithinstats]{ggwithinstats}}
  \describe{
    \item{\code{point.path}}{Logical that decides whether individual data
points and means, respectively, should be connected using \code{geom_path}. Both
default to \code{TRUE}. Note that \code{point.path} argument is relevant only when
there are two groups (i.e., in case of a \emph{t}-test). In case of large number
of data points, it is advisable to set \code{point.path = FALSE} as these lines
can overwhelm the plot.}
    \item{\code{centrality.path}}{Logical that decides whether individual data
points and means, respectively, should be connected using \code{geom_path}. Both
default to \code{TRUE}. Note that \code{point.path} argument is relevant only when
there are two groups (i.e., in case of a \emph{t}-test). In case of large number
of data points, it is advisable to set \code{point.path = FALSE} as these lines
can overwhelm the plot.}
    \item{\code{centrality.path.args}}{A list of additional aesthetic
arguments passed on to \code{geom_path} connecting raw data points and mean
points.}
    \item{\code{point.path.args}}{A list of additional aesthetic
arguments passed on to \code{geom_path} connecting raw data points and mean
points.}
    \item{\code{boxplot.args}}{A list of additional aesthetic arguments passed on to
\code{geom_boxplot}.}
    \item{\code{x}}{The grouping (or independent) variable from the dataframe \code{data}. In
case of a repeated measures or within-subjects design, if \code{subject.id}
argument is not available or not explicitly specified, the function assumes
that the data has already been sorted by such an id by the user and creates
an internal identifier. So if your data is \strong{not} sorted, the results
\emph{can} be inaccurate when there are more than two levels in \code{x} and there
are \code{NA}s present. The data is expected to be sorted by user in
subject-1,subject-2, ..., pattern.}
    \item{\code{y}}{The response (or outcome or dependent) variable from the
dataframe \code{data}.}
    \item{\code{type}}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}
    \item{\code{pairwise.comparisons}}{Logical that decides whether pairwise comparisons
are to be displayed (default: \code{TRUE}). Please note that only
\strong{significant} comparisons will be shown by default. To change this
behavior, select appropriate option with \code{pairwise.display} argument. The
pairwise comparison dataframes are prepared using the
\code{pairwise_comparisons} function. For more details
about pairwise comparisons, see the documentation for that function.}
    \item{\code{pairwise.display}}{Decides \emph{which} pairwise comparisons to display.
Available options are:
\itemize{
\item \code{"significant"} (abbreviation accepted: \code{"s"})
\item \code{"non-significant"} (abbreviation accepted: \code{"ns"})
\item \code{"all"}
}

You can use this argument to make sure that your plot is not uber-cluttered
when you have multiple groups being compared and scores of pairwise
comparisons being displayed.}
    \item{\code{p.adjust.method}}{Adjustment method for \emph{p}-values for multiple
comparisons. Possible methods are: \code{"holm"} (default), \code{"hochberg"},
\code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.}
    \item{\code{effsize.type}}{Type of effect size needed for \emph{parametric} tests. The
argument can be \code{"eta"} (partial eta-squared) or \code{"omega"} (partial
omega-squared).}
    \item{\code{bf.prior}}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors.}
    \item{\code{bf.message}}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}
    \item{\code{results.subtitle}}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}
    \item{\code{xlab}}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}
    \item{\code{ylab}}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}
    \item{\code{caption}}{The text for the plot caption.}
    \item{\code{subtitle}}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}
    \item{\code{k}}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}
    \item{\code{conf.level}}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}
    \item{\code{nboot}}{Number of bootstrap samples for computing confidence interval
for the effect size (Default: \code{100L}).}
    \item{\code{tr}}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}
    \item{\code{centrality.plotting}}{Logical that decides whether centrality tendency
measure is to be displayed as a point with a label (Default: \code{TRUE}).
Function decides which central tendency measure to show depending on the
\code{type} argument.
\itemize{
\item \strong{mean} for parametric statistics
\item \strong{median} for non-parametric statistics
\item \strong{trimmed mean} for robust statistics
\item \strong{MAP estimator} for Bayesian statistics
}

If you want default centrality parameter, you can specify this using
\code{centrality.type} argument.}
    \item{\code{centrality.type}}{Decides which centrality parameter is to be displayed.
The default is to choose the same as \code{type} argument. You can specify this
to be:
\itemize{
\item \code{"parameteric"} (for \strong{mean})
\item \code{"nonparametric"} (for \strong{median})
\item \code{robust} (for \strong{trimmed mean})
\item \code{bayes} (for \strong{MAP estimator})
}

Just as \code{type} argument, abbreviations are also accepted.}
    \item{\code{centrality.point.args}}{A list of additional aesthetic
arguments to be passed to \code{geom_point} and
\code{ggrepel::geom_label_repel} geoms, which are involved in mean plotting.}
    \item{\code{centrality.label.args}}{A list of additional aesthetic
arguments to be passed to \code{geom_point} and
\code{ggrepel::geom_label_repel} geoms, which are involved in mean plotting.}
    \item{\code{point.args}}{A list of additional aesthetic arguments to be passed to
the \code{geom_point} displaying the raw data.}
    \item{\code{outlier.tagging}}{Decides whether outliers should be tagged (Default:
\code{FALSE}).}
    \item{\code{outlier.label}}{Label to put on the outliers that have been tagged. This
\strong{can't} be the same as \code{x} argument.}
    \item{\code{outlier.coef}}{Coefficient for outlier detection using Tukey's method.
With Tukey's method, outliers are below (1st Quartile) or above (3rd
Quartile) \code{outlier.coef} times the Inter-Quartile Range (IQR) (Default:
\code{1.5}).}
    \item{\code{outlier.label.args}}{A list of additional aesthetic arguments to be
passed to \code{ggrepel::geom_label_repel} for outlier label plotting.}
    \item{\code{violin.args}}{A list of additional aesthetic arguments to be passed to
the \code{geom_violin}.}
    \item{\code{ggsignif.args}}{A list of additional aesthetic
arguments to be passed to \code{ggsignif::geom_signif}.}
    \item{\code{ggtheme}}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}
    \item{\code{package}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{palette}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{ggplot.component}}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}
  }}

\item{grouping.var}{A single grouping variable.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{plotgrid.args}{A \code{list} of additional arguments passed to
\code{patchwork::wrap_plots}, except for \code{guides} argument which is already
separately specified here.}

\item{annotation.args}{A \code{list} of additional arguments passed to
\code{patchwork::plot_annotation}.}
}
\description{
A combined plot of comparison plot created for levels of a grouping variable.
}
\examples{
\donttest{
if (require("PMCMRplus")) {
  # to get reproducible results from bootstrapping
  set.seed(123)
  library(ggstatsplot)
  library(dplyr, warn.conflicts = FALSE)
  library(ggplot2)

  # the most basic function call
  grouped_ggwithinstats(
    data             = filter(bugs_long, condition \%in\% c("HDHF", "HDLF")),
    x                = condition,
    y                = desire,
    grouping.var     = gender,
    type             = "np", # non-parametric test
    # additional modifications for **each** plot using `{ggplot2}` functions
    ggplot.component = scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 10))
  )
}
}
}
\seealso{
\code{\link{ggwithinstats}}, \code{\link{ggbetweenstats}},
\code{\link{grouped_ggbetweenstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_stats.R
\name{extract_stats}
\alias{extract_stats}
\title{Extracting dataframes with statistical details from \code{{ggstatsplot}}}
\usage{
extract_stats(p, ...)
}
\arguments{
\item{p}{A plot from \code{{ggstatsplot}} package}

\item{...}{Ignored}
}
\value{
A list of tibbles containing statistical analysis summaries.
}
\description{
Extracting dataframes with statistical details from \code{{ggstatsplot}}
}
\details{
This is a convenience function to extract dataframes with statistical details
that are used to create expressions displayed in \code{{ggstatsplot}} plots as
subtitle and/or as caption. Note that all of this analysis is carried out by
the \code{{statsExpressions}} package.

For more details about underlying tests and effect size estimates, see the
following vignette:
https://indrajeetpatil.github.io/statsExpressions/articles/stats_details.html
}
\examples{
\donttest{
if (require("PMCMRplus")) {
  set.seed(123)
  library(ggstatsplot)

  # in case of group comparisons
  p <- ggbetweenstats(mtcars, cyl, mpg)
  extract_stats(p)

  # the exact details depend on the function
  extract_stats(ggbarstats(mtcars, cyl, am))
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggpiestats.R
\name{grouped_ggpiestats}
\alias{grouped_ggpiestats}
\title{Grouped pie charts with statistical tests}
\usage{
grouped_ggpiestats(
  data,
  ...,
  grouping.var,
  output = "plot",
  plotgrid.args = list(),
  annotation.args = list()
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{...}{
  Arguments passed on to \code{\link[=ggpiestats]{ggpiestats}}
  \describe{
    \item{\code{x}}{The variable to use as the \strong{rows} in the contingency table. Please
note that if there are empty factor levels in your variable, they will be
dropped.}
    \item{\code{y}}{The variable to use as the \strong{columns} in the contingency table.
Please note that if there are empty factor levels in your variable, they
will be dropped. Default is \code{NULL}. If \code{NULL}, one-sample proportion test
(a goodness of fit test) will be run for the \code{x} variable. Otherwise an
appropriate association test will be run. This argument can not be \code{NULL}
for \code{ggbarstats} function.}
    \item{\code{proportion.test}}{Decides whether proportion test for \code{x} variable is to
be carried out for each level of \code{y}. Defaults to \code{results.subtitle}. In
\code{ggbarstats}, only \emph{p}-values from this test will be displayed.}
    \item{\code{perc.k}}{Numeric that decides number of decimal places for percentage
labels (Default: \code{0L}).}
    \item{\code{label}}{Character decides what information needs to be displayed
on the label in each pie slice. Possible options are \code{"percentage"}
(default), \code{"counts"}, \code{"both"}.}
    \item{\code{label.args}}{Additional aesthetic arguments that will be passed to
\code{geom_label}.}
    \item{\code{label.repel}}{Whether labels should be repelled using \code{ggrepel} package.
This can be helpful in case the labels are overlapping.}
    \item{\code{legend.title}}{Title text for the legend.}
    \item{\code{type}}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}
    \item{\code{results.subtitle}}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}
    \item{\code{k}}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}
    \item{\code{bf.message}}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}
    \item{\code{conf.level}}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}
    \item{\code{subtitle}}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}
    \item{\code{caption}}{The text for the plot caption.}
    \item{\code{ggtheme}}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}
    \item{\code{package}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{palette}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{ggplot.component}}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}
    \item{\code{counts}}{A string naming a variable in data containing counts, or \code{NULL}
if each row represents a single observation.}
    \item{\code{paired}}{Logical indicating whether data came from a within-subjects or
repeated measures design study (Default: \code{FALSE}). If \code{TRUE}, McNemar's
test expression will be returned. If \code{FALSE}, Pearson's chi-square test will
be returned.}
    \item{\code{ratio}}{A vector of proportions: the expected proportions for the
proportion test (should sum to 1). Default is \code{NULL}, which means the null
is equal theoretical proportions across the levels of the nominal variable.
This means if there are two levels this will be \code{ratio = c(0.5,0.5)} or if
there are four levels this will be \code{ratio = c(0.25,0.25,0.25,0.25)}, etc.}
    \item{\code{sampling.plan}}{Character describing the sampling plan. Possible options
are \code{"indepMulti"} (independent multinomial; default), \code{"poisson"},
\code{"jointMulti"} (joint multinomial), \code{"hypergeom"} (hypergeometric). For
more, see \code{?BayesFactor::contingencyTableBF()}.}
    \item{\code{fixed.margin}}{For the independent multinomial sampling plan, which
margin is fixed (\code{"rows"} or \code{"cols"}). Defaults to \code{"rows"}.}
    \item{\code{prior.concentration}}{Specifies the prior concentration parameter, set
to \code{1} by default. It indexes the expected deviation from the null
hypothesis under the alternative, and corresponds to Gunel and Dickey's
(1974) \code{"a"} parameter.}
  }}

\item{grouping.var}{A single grouping variable.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{plotgrid.args}{A \code{list} of additional arguments passed to
\code{patchwork::wrap_plots}, except for \code{guides} argument which is already
separately specified here.}

\item{annotation.args}{A \code{list} of additional arguments passed to
\code{patchwork::plot_annotation}.}
}
\description{
Helper function for \code{ggstatsplot::ggpiestats} to apply this
function across multiple levels of a given factor and combining the
resulting plots using \code{ggstatsplot::combine_plots}.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggpiestats.html}
}
\examples{
\donttest{
set.seed(123)
library(ggstatsplot)

# grouped one-sample proportion test
grouped_ggpiestats(mtcars, x = cyl, grouping.var = am)
}
}
\seealso{
\code{\link{ggbarstats}}, \code{\link{ggpiestats}},
\code{\link{grouped_ggbarstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{iris_long}
\alias{iris_long}
\title{Edgar Anderson's Iris Data in long format.}
\format{
A data frame with 600 rows and 5 variables
\itemize{
\item id. Dummy identity number for each flower (150 flowers in total).
\item Species. The species are \emph{Iris setosa}, \emph{versicolor}, and
\emph{virginica}.
\item condition. Factor giving a detailed description of the attribute
(Four levels: \code{"Petal.Length"}, \code{"Petal.Width"},  \code{"Sepal.Length"},
\code{"Sepal.Width"}).
\item attribute. What attribute is being measured (\code{"Sepal"} or \code{"Pepal"}).
\item measure. What aspect of the attribute is being measured (\code{"Length"}
or \code{"Width"}).
\item value. Value of the measurement.
}
}
\usage{
iris_long
}
\description{
Edgar Anderson's Iris Data in long format.
}
\details{
This famous (Fisher's or Anderson's) iris data set gives the
measurements in centimeters of the variables sepal length and width and
petal length and width, respectively, for 50 flowers from each of 3 species
of iris. The species are Iris setosa, versicolor, and virginica.

This is a modified dataset from \code{datasets} package.
}
\examples{
dim(iris_long)
head(iris_long)
dplyr::glimpse(iris_long)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggcorrmat.R
\name{ggcorrmat}
\alias{ggcorrmat}
\title{Visualization of a correlation matrix}
\usage{
ggcorrmat(
  data,
  cor.vars = NULL,
  cor.vars.names = NULL,
  output = "plot",
  matrix.type = "upper",
  type = "parametric",
  tr = 0.2,
  partial = FALSE,
  k = 2L,
  sig.level = 0.05,
  conf.level = 0.95,
  bf.prior = 0.707,
  p.adjust.method = "holm",
  pch = "cross",
  ggcorrplot.args = list(method = "square", outline.color = "black", pch.cex = 14),
  package = "RColorBrewer",
  palette = "Dark2",
  colors = c("#E69F00", "white", "#009E73"),
  ggtheme = ggstatsplot::theme_ggstatsplot(),
  ggplot.component = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  ...
)
}
\arguments{
\item{data}{Dataframe from which variables specified are preferentially to be
taken.}

\item{cor.vars}{List of variables for which the correlation matrix is to be
computed and visualized. If \code{NULL} (default), all numeric variables from
\code{data} will be used.}

\item{cor.vars.names}{Optional list of names to be used for \code{cor.vars}. The
names should be entered in the same order.}

\item{output}{Character that decides expected output from this function. If
\code{"plot"}, the visualization matrix will be returned. If \code{"dataframe"} (or
literally anything other than \code{"plot"}), a dataframe containing all details
from statistical analyses (e.g., correlation coefficients, statistic
values, \emph{p}-values, no. of observations, etc.) will be returned.}

\item{matrix.type}{Character, \code{"upper"} (default), \code{"lower"}, or \code{"full"},
display full matrix, lower triangular or upper triangular matrix.}

\item{type}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}

\item{tr}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}

\item{partial}{Can be \code{TRUE} for partial correlations. For Bayesian partial
correlations, "full" instead of pseudo-Bayesian partial correlations (i.e.,
Bayesian correlation based on frequentist partialization) are returned.}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{sig.level}{Significance level (Default: \code{0.05}). If the \emph{p}-value in
\emph{p}-value matrix is bigger than \code{sig.level}, then the corresponding
correlation coefficient is regarded as insignificant and flagged as such in
the plot. Relevant only when \code{output = "plot"}.}

\item{conf.level}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}

\item{bf.prior}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors and posterior estimates. In
addition to numeric arguments, several named values are also recognized:
\code{"medium"}, \code{"wide"}, and \code{"ultrawide"}, corresponding to \emph{r} scale values
of 1/2, sqrt(2)/2, and 1, respectively. In case of an ANOVA, this value
corresponds to scale for fixed effects.}

\item{p.adjust.method}{Adjustment method for \emph{p}-values for multiple
comparisons. Possible methods are: \code{"holm"} (default), \code{"hochberg"},
\code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.}

\item{pch}{Decides the point shape to be used for insignificant correlation
coefficients (only valid when \code{insig = "pch"}). Default: \code{pch = "cross"}.}

\item{ggcorrplot.args}{A list of additional (mostly aesthetic) arguments that
will be passed to \code{ggcorrplot::ggcorrplot} function. The list should avoid
any of the following arguments since they are already internally being
used: \code{corr}, \code{method}, \code{p.mat}, \code{sig.level}, \code{ggtheme}, \code{colors}, \code{lab},
\code{pch}, \code{legend.title}, \code{digits}.}

\item{package}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{palette}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{colors}{A vector of 3 colors for low, mid, and high correlation values.
If set to \code{NULL}, manual specification of colors will be turned off and 3
colors from the specified \code{palette} from \code{package} will be selected.}

\item{ggtheme}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}

\item{ggplot.component}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}

\item{title}{The text for the plot title.}

\item{subtitle}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}

\item{caption}{The text for the plot caption.}

\item{...}{Currently ignored.}
}
\description{
Correlation matrix or a dataframe containing results from pairwise
correlation tests. The package internally uses \code{ggcorrplot::ggcorrplot} for
creating the visualization matrix, while the correlation analysis is carried
out using the \code{correlation::correlation} function.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcorrmat.html}
}
\examples{
# for reproducibility
set.seed(123)
library(ggstatsplot)

# to get a plot (assumes that `ggcorrplot` is installed)
if (require("ggcorrplot")) ggcorrmat(iris)

# to get a dataframe
ggcorrmat(
  data = ggplot2::msleep,
  cor.vars = sleep_total:bodywt,
  partial = TRUE,
  output = "dataframe"
)
}
\seealso{
\code{\link{grouped_ggcorrmat}} \code{\link{ggscatterstats}}
\code{\link{grouped_ggscatterstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{movies_long}
\alias{movies_long}
\title{Movie information and user ratings from IMDB.com (long format).}
\format{
A data frame with 1,579 rows and 8 variables
\itemize{
\item title. Title of the movie.
\item year. Year of release.
\item budget. Total budget (if known) in US dollars
\item length. Length in minutes.
\item rating. Average IMDB user rating.
\item votes. Number of IMDB users who rated this movie.
\item mpaa. MPAA rating.
\item genre. Different genres of movies (action, animation, comedy, drama,
documentary, romance, short).
}
}
\source{
\url{https://CRAN.R-project.org/package=ggplot2movies}
}
\usage{
movies_long
}
\description{
Movie information and user ratings from IMDB.com (long format).
}
\details{
Modified dataset from \code{ggplot2movies} package.

The internet movie database, \url{https://imdb.com/}, is a website devoted
to collecting movie data supplied by studios and fans. It claims to be the
biggest movie database on the web and is run by amazon.
}
\examples{
dim(movies_long)
head(movies_long)
dplyr::glimpse(movies_long)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggbarstats.R
\name{grouped_ggbarstats}
\alias{grouped_ggbarstats}
\title{Grouped bar charts with statistical tests}
\usage{
grouped_ggbarstats(
  data,
  ...,
  grouping.var,
  output = "plot",
  plotgrid.args = list(),
  annotation.args = list()
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{...}{
  Arguments passed on to \code{\link[=ggbarstats]{ggbarstats}}
  \describe{
    \item{\code{xlab}}{Custom text for the \code{x} axis label (Default: \code{NULL}, which
will cause the \code{x} axis label to be the \code{x} variable).}
    \item{\code{ylab}}{Custom text for the \code{y} axis label (Default: \code{NULL}).}
    \item{\code{x}}{The variable to use as the \strong{rows} in the contingency table. Please
note that if there are empty factor levels in your variable, they will be
dropped.}
    \item{\code{y}}{The variable to use as the \strong{columns} in the contingency table.
Please note that if there are empty factor levels in your variable, they
will be dropped. Default is \code{NULL}. If \code{NULL}, one-sample proportion test
(a goodness of fit test) will be run for the \code{x} variable. Otherwise an
appropriate association test will be run. This argument can not be \code{NULL}
for \code{ggbarstats} function.}
    \item{\code{counts}}{A string naming a variable in data containing counts, or \code{NULL}
if each row represents a single observation.}
    \item{\code{type}}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}
    \item{\code{paired}}{Logical indicating whether data came from a within-subjects or
repeated measures design study (Default: \code{FALSE}). If \code{TRUE}, McNemar's
test expression will be returned. If \code{FALSE}, Pearson's chi-square test will
be returned.}
    \item{\code{results.subtitle}}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}
    \item{\code{label}}{Character decides what information needs to be displayed
on the label in each pie slice. Possible options are \code{"percentage"}
(default), \code{"counts"}, \code{"both"}.}
    \item{\code{label.args}}{Additional aesthetic arguments that will be passed to
\code{geom_label}.}
    \item{\code{k}}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}
    \item{\code{proportion.test}}{Decides whether proportion test for \code{x} variable is to
be carried out for each level of \code{y}. Defaults to \code{results.subtitle}. In
\code{ggbarstats}, only \emph{p}-values from this test will be displayed.}
    \item{\code{perc.k}}{Numeric that decides number of decimal places for percentage
labels (Default: \code{0L}).}
    \item{\code{bf.message}}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}
    \item{\code{ratio}}{A vector of proportions: the expected proportions for the
proportion test (should sum to 1). Default is \code{NULL}, which means the null
is equal theoretical proportions across the levels of the nominal variable.
This means if there are two levels this will be \code{ratio = c(0.5,0.5)} or if
there are four levels this will be \code{ratio = c(0.25,0.25,0.25,0.25)}, etc.}
    \item{\code{conf.level}}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}
    \item{\code{sampling.plan}}{Character describing the sampling plan. Possible options
are \code{"indepMulti"} (independent multinomial; default), \code{"poisson"},
\code{"jointMulti"} (joint multinomial), \code{"hypergeom"} (hypergeometric). For
more, see \code{?BayesFactor::contingencyTableBF()}.}
    \item{\code{fixed.margin}}{For the independent multinomial sampling plan, which
margin is fixed (\code{"rows"} or \code{"cols"}). Defaults to \code{"rows"}.}
    \item{\code{prior.concentration}}{Specifies the prior concentration parameter, set
to \code{1} by default. It indexes the expected deviation from the null
hypothesis under the alternative, and corresponds to Gunel and Dickey's
(1974) \code{"a"} parameter.}
    \item{\code{subtitle}}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}
    \item{\code{caption}}{The text for the plot caption.}
    \item{\code{legend.title}}{Title text for the legend.}
    \item{\code{ggtheme}}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}
    \item{\code{package}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{palette}}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}
    \item{\code{ggplot.component}}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}
  }}

\item{grouping.var}{A single grouping variable.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{plotgrid.args}{A \code{list} of additional arguments passed to
\code{patchwork::wrap_plots}, except for \code{guides} argument which is already
separately specified here.}

\item{annotation.args}{A \code{list} of additional arguments passed to
\code{patchwork::plot_annotation}.}
}
\description{
Helper function for \code{ggstatsplot::ggbarstats} to apply this function across
multiple levels of a given factor and combining the resulting plots using
\code{ggstatsplot::combine_plots}.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggpiestats.html}
}
\examples{
\donttest{
# for reproducibility
set.seed(123)
library(ggstatsplot)
library(dplyr, warn.conflicts = FALSE)

# let's create a smaller dataframe
diamonds_short <- ggplot2::diamonds \%>\%
  filter(cut \%in\% c("Very Good", "Ideal")) \%>\%
  filter(clarity \%in\% c("SI1", "SI2", "VS1", "VS2")) \%>\%
  sample_frac(size = 0.05)

# plot
grouped_ggbarstats(
  data          = diamonds_short,
  x             = color,
  y             = clarity,
  grouping.var  = cut,
  plotgrid.args = list(nrow = 2)
)
}
}
\seealso{
\code{\link{ggbarstats}}, \code{\link{ggpiestats}},
\code{\link{grouped_ggpiestats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggwithinstats.R
\name{ggwithinstats}
\alias{ggwithinstats}
\title{Box/Violin plots for within-subjects (or repeated measures) comparisons}
\usage{
ggwithinstats(
  data,
  x,
  y,
  type = "parametric",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  p.adjust.method = "holm",
  effsize.type = "unbiased",
  bf.prior = 0.707,
  bf.message = TRUE,
  results.subtitle = TRUE,
  xlab = NULL,
  ylab = NULL,
  caption = NULL,
  title = NULL,
  subtitle = NULL,
  k = 2L,
  conf.level = 0.95,
  nboot = 100L,
  tr = 0.2,
  centrality.plotting = TRUE,
  centrality.type = type,
  centrality.point.args = list(size = 5, color = "darkred"),
  centrality.label.args = list(size = 3, nudge_x = 0.4, segment.linetype = 4),
  centrality.path = TRUE,
  centrality.path.args = list(size = 1, color = "red", alpha = 0.5),
  point.args = list(size = 3, alpha = 0.5),
  point.path = TRUE,
  point.path.args = list(alpha = 0.5, linetype = "dashed"),
  outlier.tagging = FALSE,
  outlier.label = NULL,
  outlier.coef = 1.5,
  outlier.label.args = list(size = 3),
  boxplot.args = list(width = 0.2, alpha = 0.5),
  violin.args = list(width = 0.5, alpha = 0.2),
  ggsignif.args = list(textsize = 3, tip_length = 0.01),
  ggtheme = ggstatsplot::theme_ggstatsplot(),
  package = "RColorBrewer",
  palette = "Dark2",
  ggplot.component = NULL,
  output = "plot",
  ...
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{x}{The grouping (or independent) variable from the dataframe \code{data}. In
case of a repeated measures or within-subjects design, if \code{subject.id}
argument is not available or not explicitly specified, the function assumes
that the data has already been sorted by such an id by the user and creates
an internal identifier. So if your data is \strong{not} sorted, the results
\emph{can} be inaccurate when there are more than two levels in \code{x} and there
are \code{NA}s present. The data is expected to be sorted by user in
subject-1,subject-2, ..., pattern.}

\item{y}{The response (or outcome or dependent) variable from the
dataframe \code{data}.}

\item{type}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}

\item{pairwise.comparisons}{Logical that decides whether pairwise comparisons
are to be displayed (default: \code{TRUE}). Please note that only
\strong{significant} comparisons will be shown by default. To change this
behavior, select appropriate option with \code{pairwise.display} argument. The
pairwise comparison dataframes are prepared using the
\code{pairwise_comparisons} function. For more details
about pairwise comparisons, see the documentation for that function.}

\item{pairwise.display}{Decides \emph{which} pairwise comparisons to display.
Available options are:
\itemize{
\item \code{"significant"} (abbreviation accepted: \code{"s"})
\item \code{"non-significant"} (abbreviation accepted: \code{"ns"})
\item \code{"all"}
}

You can use this argument to make sure that your plot is not uber-cluttered
when you have multiple groups being compared and scores of pairwise
comparisons being displayed.}

\item{p.adjust.method}{Adjustment method for \emph{p}-values for multiple
comparisons. Possible methods are: \code{"holm"} (default), \code{"hochberg"},
\code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.}

\item{effsize.type}{Type of effect size needed for \emph{parametric} tests. The
argument can be \code{"eta"} (partial eta-squared) or \code{"omega"} (partial
omega-squared).}

\item{bf.prior}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors.}

\item{bf.message}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}

\item{results.subtitle}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}

\item{xlab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{ylab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{caption}{The text for the plot caption.}

\item{title}{The text for the plot title.}

\item{subtitle}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{conf.level}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}

\item{nboot}{Number of bootstrap samples for computing confidence interval
for the effect size (Default: \code{100L}).}

\item{tr}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}

\item{centrality.plotting}{Logical that decides whether centrality tendency
measure is to be displayed as a point with a label (Default: \code{TRUE}).
Function decides which central tendency measure to show depending on the
\code{type} argument.
\itemize{
\item \strong{mean} for parametric statistics
\item \strong{median} for non-parametric statistics
\item \strong{trimmed mean} for robust statistics
\item \strong{MAP estimator} for Bayesian statistics
}

If you want default centrality parameter, you can specify this using
\code{centrality.type} argument.}

\item{centrality.type}{Decides which centrality parameter is to be displayed.
The default is to choose the same as \code{type} argument. You can specify this
to be:
\itemize{
\item \code{"parameteric"} (for \strong{mean})
\item \code{"nonparametric"} (for \strong{median})
\item \code{robust} (for \strong{trimmed mean})
\item \code{bayes} (for \strong{MAP estimator})
}

Just as \code{type} argument, abbreviations are also accepted.}

\item{centrality.point.args}{A list of additional aesthetic
arguments to be passed to \code{geom_point} and
\code{ggrepel::geom_label_repel} geoms, which are involved in mean plotting.}

\item{centrality.label.args}{A list of additional aesthetic
arguments to be passed to \code{geom_point} and
\code{ggrepel::geom_label_repel} geoms, which are involved in mean plotting.}

\item{centrality.path.args, point.path.args}{A list of additional aesthetic
arguments passed on to \code{geom_path} connecting raw data points and mean
points.}

\item{point.args}{A list of additional aesthetic arguments to be passed to
the \code{geom_point} displaying the raw data.}

\item{point.path, centrality.path}{Logical that decides whether individual data
points and means, respectively, should be connected using \code{geom_path}. Both
default to \code{TRUE}. Note that \code{point.path} argument is relevant only when
there are two groups (i.e., in case of a \emph{t}-test). In case of large number
of data points, it is advisable to set \code{point.path = FALSE} as these lines
can overwhelm the plot.}

\item{outlier.tagging}{Decides whether outliers should be tagged (Default:
\code{FALSE}).}

\item{outlier.label}{Label to put on the outliers that have been tagged. This
\strong{can't} be the same as \code{x} argument.}

\item{outlier.coef}{Coefficient for outlier detection using Tukey's method.
With Tukey's method, outliers are below (1st Quartile) or above (3rd
Quartile) \code{outlier.coef} times the Inter-Quartile Range (IQR) (Default:
\code{1.5}).}

\item{outlier.label.args}{A list of additional aesthetic arguments to be
passed to \code{ggrepel::geom_label_repel} for outlier label plotting.}

\item{boxplot.args}{A list of additional aesthetic arguments passed on to
\code{geom_boxplot}.}

\item{violin.args}{A list of additional aesthetic arguments to be passed to
the \code{geom_violin}.}

\item{ggsignif.args}{A list of additional aesthetic
arguments to be passed to \code{ggsignif::geom_signif}.}

\item{ggtheme}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}

\item{package}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{palette}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{ggplot.component}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{...}{Currently ignored.}
}
\description{
A combination of box and violin plots along with raw (unjittered) data points
for within-subjects designs with statistical details included in the plot as
a subtitle.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggwithinstats.html}
}
\examples{
\donttest{
if (require("PMCMRplus")) {
  # setup
  set.seed(123)
  library(ggstatsplot)
  library(dplyr, warn.conflicts = FALSE)

  # two groups (*t*-test)
  ggwithinstats(
    data = filter(bugs_long, condition \%in\% c("HDHF", "HDLF")),
    x    = condition,
    y    = desire
  )

  # more than two groups (anova)
  library(WRS2)

  ggwithinstats(
    data            = WineTasting,
    x               = Wine,
    y               = Taste,
    type            = "r",
    outlier.tagging = TRUE,
    outlier.label   = Taster
  )
}
}
}
\seealso{
\code{\link{grouped_ggbetweenstats}}, \code{\link{ggbetweenstats}},
\code{\link{grouped_ggwithinstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reexports.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\alias{\%<>\%}
\alias{\%$\%}
\alias{\%<-\%}
\alias{tibble}
\alias{enframe}
\alias{as_tibble}
\alias{exec}
\alias{!!}
\alias{!!!}
\alias{\%|\%}
\alias{\%||\%}
\alias{:=}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{rlang}{\code{\link[rlang:nse-force]{!!}}, \code{\link[rlang:nse-force]{!!!}}, \code{\link[rlang:op-na-default]{\%|\%}}, \code{\link[rlang:op-null-default]{\%||\%}}, \code{\link[rlang:nse-force]{:=}}, \code{\link[rlang]{exec}}}

  \item{statsExpressions}{\code{\link[statsExpressions:reexports]{\%$\%}}, \code{\link[statsExpressions:reexports]{\%<-\%}}, \code{\link[statsExpressions:reexports]{\%<>\%}}, \code{\link[statsExpressions:reexports]{\%>\%}}, \code{\link[statsExpressions:reexports]{as_tibble}}, \code{\link[statsExpressions:reexports]{enframe}}, \code{\link[statsExpressions:reexports]{tibble}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggstatsplot-package.R
\docType{package}
\name{ggstatsplot-package}
\alias{ggstatsplot-package}
\alias{_PACKAGE}
\alias{ggstatsplot}
\title{ggstatsplot: 'ggplot2' Based Plots with Statistical Details}
\description{
\code{{ggstatsplot}} is an extension of \code{{ggplot2}} package. It creates
graphics with details from statistical tests included in the plots
themselves. It provides an easier \code{API} to generate information-rich plots
for statistical analysis of continuous (violin plots, scatterplots,
histograms, dot plots, dot-and-whisker plots) or categorical (pie and bar
charts) data. Currently, it supports the most common types of statistical
tests: parametric, nonparametric, robust, and Bayesian versions of
\emph{t}-test/ANOVA, correlation analyses, contingency table analysis,
meta-analysis, and regression analyses.
}
\details{
\code{ggstatsplot}

The main functions are:
\itemize{
\item \code{ggbetweenstats} function to produce information-rich comparison plot
\emph{between} different groups or conditions with \code{{ggplot2}} and details from
the statistical tests in the subtitle.
\item \code{ggwithinstats} function to produce information-rich comparison plot
\emph{within} different groups or conditions with \code{{ggplot2}} and details from the
statistical tests in the subtitle.
\item \code{ggscatterstats} function to produce \code{{ggplot2}} scatterplots along with a
marginal distribution plots from \code{ggside} package and details from the
statistical tests in the subtitle.
\item \code{ggpiestats} function to produce pie chart with details from the
statistical tests in the subtitle.
\item \code{ggbarstats} function to produce stacked bar chart with details from the
statistical tests in the subtitle.
\item \code{gghistostats} function to produce histogram for a single variable with
results from one sample test displayed in the subtitle.
\item \code{ggdotplotstats} function to produce Cleveland-style dot plots/charts for
a single variable with labels and results from one sample test displayed in
the subtitle.
\item \code{ggcorrmat} function to visualize the correlation matrix.
\item \code{ggcoefstats} function to visualize results from regression analyses.
\item \code{combine_plots} helper function to combine multiple \code{{ggstatsplot}} plots
using \code{patchwork::wrap_plots()}.
}

For more documentation, see the dedicated
\href{https://indrajeetpatil.github.io/ggstatsplot/}{Website}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://indrajeetpatil.github.io/ggstatsplot/}
  \item \url{https://github.com/IndrajeetPatil/ggstatsplot}
  \item Report bugs at \url{https://github.com/IndrajeetPatil/ggstatsplot/issues}
}

}
\author{
\strong{Maintainer}: Indrajeet Patil \email{patilindrajeet.science@gmail.com} (\href{https://orcid.org/0000-0003-1995-6531}{ORCID}) (@patilindrajeets) [copyright holder]

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_plots.R
\name{combine_plots}
\alias{combine_plots}
\title{Combining and arranging multiple plots in a grid}
\usage{
combine_plots(
  plotlist,
  plotgrid.args = list(),
  annotation.args = list(),
  guides = "collect",
  ...
)
}
\arguments{
\item{plotlist}{A list containing \code{ggplot} objects.}

\item{plotgrid.args}{A \code{list} of additional arguments passed to
\code{patchwork::wrap_plots}, except for \code{guides} argument which is already
separately specified here.}

\item{annotation.args}{A \code{list} of additional arguments passed to
\code{patchwork::plot_annotation}.}

\item{guides}{A string specifying how guides should be treated in the layout.
\code{'collect'} will collect guides below to the given nesting level, removing
duplicates. \code{'keep'} will stop collection at this level and let guides be
placed alongside their plot. \code{auto} will allow guides to be collected if a
upper level tries, but place them alongside the plot if not.  If you modify
default guide "position" with \link[ggplot2:theme]{theme(legend.position=...)}
while also collecting guides you must apply that change to the overall
patchwork (see example).}

\item{...}{Currently ignored.}
}
\value{
Combined plot with annotation labels
}
\description{
Wrapper around \code{patchwork::wrap_plots} that will return a combined grid of
plots with annotations. In case you want to create a grid of plots, it is
\strong{highly recommended} that you use \code{{patchwork}} package directly and not
this wrapper around it which is mostly useful with \code{{ggstatsplot}} plots. It
is exported only for backward compatibility.
}
\examples{
# loading the necessary libraries
library(ggplot2)

# preparing the first plot
p1 <- ggplot(
  data = subset(iris, iris$Species == "setosa"),
  aes(x = Sepal.Length, y = Sepal.Width)
) +
  geom_point() +
  labs(title = "setosa")

# preparing the second plot
p2 <- ggplot(
  data = subset(iris, iris$Species == "versicolor"),
  aes(x = Sepal.Length, y = Sepal.Width)
) +
  geom_point() +
  labs(title = "versicolor")

# combining the plot with a title and a caption
combine_plots(
  plotlist = list(p1, p2),
  plotgrid.args = list(nrow = 1),
  annotation.args = list(
    tag_levels = "a",
    title = "Dataset: Iris Flower dataset",
    subtitle = "Edgar Anderson collected this data",
    caption = "Note: Only two species of flower are displayed",
    theme = theme(
      plot.subtitle = element_text(size = 20),
      plot.title = element_text(size = 30)
    )
  )
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bugs_wide}
\alias{bugs_wide}
\title{Wide-format version of the "Bugs" dataset.}
\format{
A data frame with 93 rows and 6 variables
\itemize{
\item subject. Dummy identity number for each participant.
\item gender. Participant's gender (Female, Male).
\item region. Region of the world the participant was from.
\item education. Level of education.
\item ldlf,ldhf,hdlf,hdhf.The desire to kill an arthropod was indicated on
a scale from 0 to 10 in each condition of the experiment (\strong{LDLF}: low
freighteningness and low disgustingness; \strong{LFHD}: low freighteningness and
high disgustingness; \strong{HFHD}: high freighteningness and low
disgustingness; \strong{HFHD}: high freighteningness and high disgustingness).
}
}
\source{
\url{https://www.sciencedirect.com/science/article/pii/S0747563213000277}
}
\usage{
bugs_wide
}
\description{
Wide-format version of the "Bugs" dataset.
}
\details{
This data set, "Bugs", provides the extent to which men and women
want to kill arthropods that vary in freighteningness (low, high) and
disgustingness (low, high). Each participant rates their attitudes towards
all anthropods. Subset of the data reported by Ryan et al. (2013).
}
\examples{
dim(bugs_wide)
head(bugs_wide)
dplyr::glimpse(bugs_wide)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggdotplotstats.R
\name{ggdotplotstats}
\alias{ggdotplotstats}
\title{Dot plot/chart for labeled numeric data.}
\usage{
ggdotplotstats(
  data,
  x,
  y,
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  type = "parametric",
  test.value = 0,
  bf.prior = 0.707,
  bf.message = TRUE,
  effsize.type = "g",
  conf.level = 0.95,
  tr = 0.2,
  k = 2L,
  results.subtitle = TRUE,
  point.args = list(color = "black", size = 3, shape = 16),
  centrality.plotting = TRUE,
  centrality.type = type,
  centrality.line.args = list(color = "blue", size = 1, linetype = "dashed"),
  ggplot.component = NULL,
  ggtheme = ggstatsplot::theme_ggstatsplot(),
  output = "plot",
  ...
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{x}{A numeric variable from the dataframe \code{data}.}

\item{y}{Label or grouping variable.}

\item{xlab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{ylab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{title}{The text for the plot title.}

\item{subtitle}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}

\item{caption}{The text for the plot caption.}

\item{type}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}

\item{test.value}{A number indicating the true value of the mean (Default:
\code{0}).}

\item{bf.prior}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors and posterior estimates. In
addition to numeric arguments, several named values are also recognized:
\code{"medium"}, \code{"wide"}, and \code{"ultrawide"}, corresponding to \emph{r} scale values
of 1/2, sqrt(2)/2, and 1, respectively. In case of an ANOVA, this value
corresponds to scale for fixed effects.}

\item{bf.message}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}

\item{effsize.type}{Type of effect size needed for \emph{parametric} tests. The
argument can be \code{"d"} (for Cohen's \emph{d}) or \code{"g"} (for Hedge's \emph{g}).}

\item{conf.level}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}

\item{tr}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{results.subtitle}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}

\item{point.args}{A list of additional aesthetic arguments passed to
\code{geom_point}.}

\item{centrality.plotting}{Logical that decides whether centrality tendency
measure is to be displayed as a point with a label (Default: \code{TRUE}).
Function decides which central tendency measure to show depending on the
\code{type} argument.
\itemize{
\item \strong{mean} for parametric statistics
\item \strong{median} for non-parametric statistics
\item \strong{trimmed mean} for robust statistics
\item \strong{MAP estimator} for Bayesian statistics
}

If you want default centrality parameter, you can specify this using
\code{centrality.type} argument.}

\item{centrality.type}{Decides which centrality parameter is to be displayed.
The default is to choose the same as \code{type} argument. You can specify this
to be:
\itemize{
\item \code{"parameteric"} (for \strong{mean})
\item \code{"nonparametric"} (for \strong{median})
\item \code{robust} (for \strong{trimmed mean})
\item \code{bayes} (for \strong{MAP estimator})
}

Just as \code{type} argument, abbreviations are also accepted.}

\item{centrality.line.args}{A list of additional aesthetic arguments to be
passed to the \code{geom_line} used to display the lines corresponding to the
centrality parameter.}

\item{ggplot.component}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}

\item{ggtheme}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{...}{Currently ignored.}
}
\description{
A dot chart (as described by William S. Cleveland) with statistical details
from one-sample test details.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggdotplotstats.html}
}
\examples{
\donttest{
# for reproducibility
set.seed(123)

# plot
ggdotplotstats(
  data = ggplot2::mpg,
  x = cty,
  y = manufacturer,
  title = "Fuel economy data",
  xlab = "city miles per gallon"
)
}
}
\seealso{
\code{\link{grouped_gghistostats}}, \code{\link{gghistostats}},
\code{\link{grouped_ggdotplotstats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairwise_comparisons.R
\name{pairwise_comparisons}
\alias{pairwise_comparisons}
\title{Multiple pairwise comparison tests with tidy data}
\usage{
pairwise_comparisons(
  data,
  x,
  y,
  subject.id = NULL,
  type = "parametric",
  paired = FALSE,
  var.equal = FALSE,
  tr = 0.2,
  bf.prior = 0.707,
  p.adjust.method = "holm",
  k = 2L,
  ...
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{x}{The grouping (or independent) variable from the dataframe \code{data}. In
case of a repeated measures or within-subjects design, if \code{subject.id}
argument is not available or not explicitly specified, the function assumes
that the data has already been sorted by such an id by the user and creates
an internal identifier. So if your data is \strong{not} sorted, the results
\emph{can} be inaccurate when there are more than two levels in \code{x} and there
are \code{NA}s present. The data is expected to be sorted by user in
subject-1,subject-2, ..., pattern.}

\item{y}{The response (or outcome or dependent) variable from the
dataframe \code{data}.}

\item{subject.id}{Relevant in case of a repeated measures or within-subjects
design (\code{paired = TRUE}, i.e.), it specifies the subject or repeated
measures identifier. \strong{Important}: Note that if this argument is \code{NULL}
(which is the default), the function assumes that the data has already been
sorted by such an id by the user and creates an internal identifier. So if
your data is \strong{not} sorted and you leave this argument unspecified, the
results \emph{can} be inaccurate when there are more than two levels in \code{x} and
there are \code{NA}s present.}

\item{type}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}

\item{paired}{Logical that decides whether the experimental design is
repeated measures/within-subjects or between-subjects. The default is
\code{FALSE}.}

\item{var.equal}{a logical variable indicating whether to treat the
    two variances as being equal. If \code{TRUE} then the pooled
    variance is used to estimate the variance otherwise the Welch
    (or Satterthwaite) approximation to the degrees of freedom is used.}

\item{tr}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}

\item{bf.prior}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors and posterior estimates. In
addition to numeric arguments, several named values are also recognized:
\code{"medium"}, \code{"wide"}, and \code{"ultrawide"}, corresponding to \emph{r} scale values
of 1/2, sqrt(2)/2, and 1, respectively. In case of an ANOVA, this value
corresponds to scale for fixed effects.}

\item{p.adjust.method}{Adjustment method for \emph{p}-values for multiple
comparisons. Possible methods are: \code{"holm"} (default), \code{"hochberg"},
\code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{...}{Additional arguments passed to other methods.}
}
\value{
A tibble dataframe containing two columns corresponding to group
levels being compared with each other (\code{group1} and \code{group2}) and \code{p.value}
column corresponding to this comparison. The dataframe will also contain a
\code{p.value.label} column containing a \emph{label} for this \emph{p}-value, in case
this needs to be displayed in \code{ggsignif::geom_ggsignif}. In addition to
these common columns across the different types of statistics, there will
be additional columns specific to the \code{type} of test being run.

This function provides a unified syntax to carry out pairwise comparison
tests and internally relies on other packages to carry out these tests. For
more details about the included tests, see the documentation for the
respective functions:
\itemize{
\item \emph{parametric} : \code{\link[stats:pairwise.t.test]{stats::pairwise.t.test()}} (paired) and
\code{\link[PMCMRplus:gamesHowellTest]{PMCMRplus::gamesHowellTest()}} (unpaired)
\item \emph{non-parametric} :
\code{\link[PMCMRplus:durbinAllPairsTest]{PMCMRplus::durbinAllPairsTest()}} (paired) and
\code{\link[PMCMRplus:kwAllPairsDunnTest]{PMCMRplus::kwAllPairsDunnTest()}} (unpaired)
\item \emph{robust} :
\code{\link[WRS2:rmanova]{WRS2::rmmcp()}} (paired) and \code{\link[WRS2:t1way]{WRS2::lincon()}} (unpaired)
\item \emph{Bayes Factor} : \code{\link[BayesFactor:ttestBF]{BayesFactor::ttestBF()}}
}
}
\description{
Calculate parametric, non-parametric, robust, and Bayes Factor pairwise
comparisons between group levels with corrections for multiple testing.
}
\examples{
\donttest{
if (require("PMCMRplus")) {
  # for reproducibility
  set.seed(123)
  library(ggstatsplot)
  library(statsExpressions) # for data

  # show all columns and make the column titles bold
  # as a user, you don't need to do this; this is just for the package website
  options(tibble.width = Inf, pillar.bold = TRUE, pillar.neg = TRUE, pillar.subtle_num = TRUE)

  #------------------- between-subjects design ----------------------------

  # parametric
  # if `var.equal = TRUE`, then Student's t-test will be run
  pairwise_comparisons(
    data            = mtcars,
    x               = cyl,
    y               = wt,
    type            = "parametric",
    var.equal       = TRUE,
    paired          = FALSE,
    p.adjust.method = "none"
  )

  # if `var.equal = FALSE`, then Games-Howell test will be run
  pairwise_comparisons(
    data            = mtcars,
    x               = cyl,
    y               = wt,
    type            = "parametric",
    var.equal       = FALSE,
    paired          = FALSE,
    p.adjust.method = "bonferroni"
  )

  # non-parametric (Dunn test)
  pairwise_comparisons(
    data            = mtcars,
    x               = cyl,
    y               = wt,
    type            = "nonparametric",
    paired          = FALSE,
    p.adjust.method = "none"
  )

  # robust (Yuen's trimmed means *t*-test)
  pairwise_comparisons(
    data            = mtcars,
    x               = cyl,
    y               = wt,
    type            = "robust",
    paired          = FALSE,
    p.adjust.method = "fdr"
  )

  # Bayes Factor (Student's *t*-test)
  pairwise_comparisons(
    data   = mtcars,
    x      = cyl,
    y      = wt,
    type   = "bayes",
    paired = FALSE
  )

  #------------------- within-subjects design ----------------------------

  # parametric (Student's *t*-test)
  pairwise_comparisons(
    data            = bugs_long,
    x               = condition,
    y               = desire,
    subject.id      = subject,
    type            = "parametric",
    paired          = TRUE,
    p.adjust.method = "BH"
  )

  # non-parametric (Durbin-Conover test)
  pairwise_comparisons(
    data            = bugs_long,
    x               = condition,
    y               = desire,
    subject.id      = subject,
    type            = "nonparametric",
    paired          = TRUE,
    p.adjust.method = "BY"
  )

  # robust (Yuen's trimmed means t-test)
  pairwise_comparisons(
    data            = bugs_long,
    x               = condition,
    y               = desire,
    subject.id      = subject,
    type            = "robust",
    paired          = TRUE,
    p.adjust.method = "hommel"
  )

  # Bayes Factor (Student's *t*-test)
  pairwise_comparisons(
    data       = bugs_long,
    x          = condition,
    y          = desire,
    subject.id = subject,
    type       = "bayes",
    paired     = TRUE
  )
}
}
}
\references{
For more, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/pairwise.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{grouped_list}
\alias{grouped_list}
\title{Split dataframe into a list by grouping variable.}
\usage{
grouped_list(data, grouping.var = NULL)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{grouping.var}{A single grouping variable.}
}
\description{
This function splits the dataframe into a list, with the length
of the list equal to the factor levels of the grouping variable. Each
element of the list will be a tibble.
}
\examples{
\donttest{
ggstatsplot:::grouped_list(ggplot2::msleep, grouping.var = vore)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggcoefstats.R
\name{ggcoefstats}
\alias{ggcoefstats}
\title{Dot-and-whisker plots for regression analyses}
\usage{
ggcoefstats(
  x,
  output = "plot",
  statistic = NULL,
  conf.int = TRUE,
  conf.level = 0.95,
  k = 2L,
  exclude.intercept = FALSE,
  effsize = "eta",
  meta.analytic.effect = FALSE,
  meta.type = "parametric",
  bf.message = TRUE,
  sort = "none",
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  only.significant = FALSE,
  point.args = list(size = 3, color = "blue"),
  errorbar.args = list(height = 0),
  vline = TRUE,
  vline.args = list(size = 1, linetype = "dashed"),
  stats.labels = TRUE,
  stats.label.color = NULL,
  stats.label.args = list(size = 3, direction = "y", min.segment.length = 0),
  package = "RColorBrewer",
  palette = "Dark2",
  ggtheme = ggstatsplot::theme_ggstatsplot(),
  ...
)
}
\arguments{
\item{x}{A model object to be tidied, or a tidy data frame containing results
from a regression model. Function internally uses
\code{parameters::model_parameters()} to get a tidy dataframe. If a dataframe is
entered, it \emph{must} contain at the minimum two columns named \code{term} (names
of predictors) and \code{estimate} (corresponding estimates of coefficients or
other quantities of interest).}

\item{output}{Character describing the expected output from this function:
\code{"plot"} (visualization of regression coefficients) or \code{"tidy"} (tidy
dataframe of results \code{parameters::model_parameters}) or \code{"glance"} (object
from \code{performance::model_performance}).}

\item{statistic}{Which statistic is to be displayed (either \code{"t"} or \code{"f"}or
\code{"z"} or \code{"chi"}) in the label. This is relevant if the \code{x} argument is a
\emph{dataframe}.}

\item{conf.int}{Logical. Decides whether to display confidence intervals as
error bars (Default: \code{TRUE}).}

\item{conf.level}{Numeric deciding level of confidence or credible intervals
(Default: \code{0.95}).}

\item{k}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}

\item{exclude.intercept}{Logical that decides whether the intercept should be
excluded from the plot (Default: \code{FALSE}).}

\item{effsize}{Character describing the effect size to be displayed: \code{"eta"}
(default) or \code{"omega"}. This argument is relevant only for models objects
with \emph{F}-statistic.}

\item{meta.analytic.effect}{Logical that decides whether subtitle for
meta-analysis via linear (mixed-effects) models (default: \code{FALSE}). If
\code{TRUE}, input to argument \code{subtitle} will be ignored. This will be mostly
relevant if a data frame with estimates and their standard errors is
entered.}

\item{meta.type}{Type of statistics used to carry out random-effects
meta-analysis. If \code{"parametric"} (default), \code{metafor::rma} function will be
used. If \code{"robust"}, \code{metaplus::metaplus} function will be used. If
\code{"bayes"}, \code{metaBMA::meta_random} function will be used.}

\item{bf.message}{Logical that decides whether results from running a
Bayesian meta-analysis assuming that the effect size \emph{d} varies across
studies with standard deviation \emph{t} (i.e., a random-effects analysis)
should be displayed in caption. Defaults to \code{TRUE}.}

\item{sort}{If \code{"none"} (default) do not sort, \code{"ascending"} sort by
increasing coefficient value, or \code{"descending"} sort by decreasing
coefficient value.}

\item{xlab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{ylab}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}

\item{title}{The text for the plot title.}

\item{subtitle}{The text for the plot subtitle. The input to this argument
will be ignored if \code{meta.analytic.effect} is set to \code{TRUE}.}

\item{caption}{The text for the plot caption.}

\item{only.significant}{If \code{TRUE}, only stats labels for significant effects
is shown (Default: \code{FALSE}). This can be helpful when a large number of
regression coefficients are to be displayed in a single plot. Relevant only
when the \code{output} is a plot.}

\item{point.args}{Additional arguments that will be passed to
\code{geom_point} geom. Please see documentation for that function to
know more about these arguments.}

\item{errorbar.args}{Additional arguments that will be passed to
\code{geom_errorbarh} geom. Please see documentation for that function
to know more about these arguments.}

\item{vline}{Decides whether to display a vertical line (Default: \code{"TRUE"}).}

\item{vline.args}{Additional arguments that will be passed to
\code{geom_vline} geom. Please see documentation for that function to
know more about these arguments.}

\item{stats.labels}{Logical. Decides whether the statistic and \emph{p}-values for
each coefficient are to be attached to each dot as a text label using
\code{ggrepel} (Default: \code{TRUE}).}

\item{stats.label.color}{Color for the labels. If set to \code{NULL}, colors will
be chosen from the specified \code{package} (Default: \code{"RColorBrewer"}) and
\code{palette} (Default: \code{"Dark2"}).}

\item{stats.label.args}{Additional arguments that will be passed to
\code{ggrepel::geom_label_repel} geom. Please see documentation for that
function to know more about these arguments.}

\item{package}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{palette}{Name of the package from which the given palette is to
be extracted. The available palettes and packages can be checked by running
\code{View(paletteer::palettes_d_names)}.}

\item{ggtheme}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}

\item{...}{Additional arguments to tidying method. For more, see
\code{parameters::model_parameters}.}
}
\description{
Plot with the regression coefficients' point estimates as dots with
confidence interval whiskers and other statistical details included as
labels.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggcoefstats.html}
}
\note{
\enumerate{
\item In case you want to carry out meta-analysis, you will be asked to install
the needed packages (\code{{metafor}}, \code{{metaplus}}, or \code{{metaBMA}}) for
meta-analysis (if unavailable).
\item All rows of regression estimates where either of the following
quantities is \code{NA} will be removed if labels are requested: \code{estimate},
\code{statistic}, \code{p.value}.
\item Given the rapid pace at which new methods are added to these packages, it
is recommended that you install the GitHub versions of \code{{parameters}} and
\code{{performance}} in order to make most of this function.
}
}
\examples{
\donttest{
# for reproducibility
set.seed(123)
library(ggstatsplot)

# model object
mod <- lm(formula = mpg ~ cyl * am, data = mtcars)

# to get a plot
ggcoefstats(mod, output = "plot")

# to get a tidy dataframe
ggcoefstats(mod, output = "tidy")

# to get a glance summary
ggcoefstats(mod, output = "glance")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggscatterstats.R
\name{grouped_ggscatterstats}
\alias{grouped_ggscatterstats}
\title{Scatterplot with marginal distributions for all levels of a grouping
variable}
\usage{
grouped_ggscatterstats(
  data,
  ...,
  grouping.var,
  output = "plot",
  plotgrid.args = list(),
  annotation.args = list()
)
}
\arguments{
\item{data}{A dataframe (or a tibble) from which variables specified are to
be taken. Other data types (e.g., matrix,table, array, etc.) will \strong{not}
be accepted.}

\item{...}{
  Arguments passed on to \code{\link[=ggscatterstats]{ggscatterstats}}
  \describe{
    \item{\code{label.var}}{Variable to use for points labels entered as a symbol (e.g.
\code{var1}).}
    \item{\code{label.expression}}{An expression evaluating to a logical vector that
determines the subset of data points to label (e.g. \code{y < 4 & z < 20}).
While using this argument with \code{purrr::pmap}, you will have to provide a
quoted expression  (e.g. \code{quote(y < 4 & z < 20)}).}
    \item{\code{point.label.args}}{A list of additional aesthetic arguments to be passed
to \code{ggrepel::geom_label_repel} geom used to display the labels.}
    \item{\code{smooth.line.args}}{A list of additional aesthetic arguments to be passed
to \code{geom_smooth} geom used to display the regression line.}
    \item{\code{point.args}}{A list of additional aesthetic arguments to be passed
to \code{geom_point} geom used to display the raw data points.}
    \item{\code{marginal}}{Decides whether marginal distributions will be plotted on
axes using \code{ggside} functions. The default is \code{TRUE}. The package
\code{ggside} must already be installed by the user.}
    \item{\code{point.width.jitter}}{Degree of jitter in \code{x} and \code{y}
direction, respectively. Defaults to \code{0} (0\%) of the resolution of the
data. Note that the jitter should not be specified in the \code{point.args}
because this information will be passed to two different \code{geom}s: one
displaying the \strong{points} and the other displaying the *\strong{labels} for
these points.}
    \item{\code{point.height.jitter}}{Degree of jitter in \code{x} and \code{y}
direction, respectively. Defaults to \code{0} (0\%) of the resolution of the
data. Note that the jitter should not be specified in the \code{point.args}
because this information will be passed to two different \code{geom}s: one
displaying the \strong{points} and the other displaying the *\strong{labels} for
these points.}
    \item{\code{xfill}}{Character describing color fill for \code{x} and \code{y} axes
marginal distributions (default: \code{"#009E73"} (for \code{x}) and \code{"#D55E00"} (for
\code{y})). Note that the defaults are colorblind-friendly.}
    \item{\code{yfill}}{Character describing color fill for \code{x} and \code{y} axes
marginal distributions (default: \code{"#009E73"} (for \code{x}) and \code{"#D55E00"} (for
\code{y})). Note that the defaults are colorblind-friendly.}
    \item{\code{xsidehistogram.args}}{A list of arguments passed to
respective \code{geom_}s from \code{ggside} package to change the marginal
distribution histograms plots.}
    \item{\code{ysidehistogram.args}}{A list of arguments passed to
respective \code{geom_}s from \code{ggside} package to change the marginal
distribution histograms plots.}
    \item{\code{x}}{The column in \code{data} containing the explanatory variable to be
plotted on the \code{x}-axis.}
    \item{\code{y}}{The column in \code{data} containing the response (outcome) variable to
be plotted on the \code{y}-axis.}
    \item{\code{type}}{A character specifying the type of statistical approach:
\itemize{
\item \code{"parametric"}
\item \code{"nonparametric"}
\item \code{"robust"}
\item \code{"bayes"}
}

You can specify just the initial letter.}
    \item{\code{conf.level}}{Scalar between \code{0} and \code{1}. If unspecified, the defaults
return \verb{95\%} confidence/credible intervals (\code{0.95}).}
    \item{\code{bf.prior}}{A number between \code{0.5} and \code{2} (default \code{0.707}), the prior
width to use in calculating Bayes factors and posterior estimates. In
addition to numeric arguments, several named values are also recognized:
\code{"medium"}, \code{"wide"}, and \code{"ultrawide"}, corresponding to \emph{r} scale values
of 1/2, sqrt(2)/2, and 1, respectively. In case of an ANOVA, this value
corresponds to scale for fixed effects.}
    \item{\code{tr}}{Trim level for the mean when carrying out \code{robust} tests. In case
of an error, try reducing the value of \code{tr}, which is by default set to
\code{0.2}. Lowering the value might help.}
    \item{\code{k}}{Number of digits after decimal point (should be an integer)
(Default: \code{k = 2L}).}
    \item{\code{bf.message}}{Logical that decides whether to display Bayes Factor in
favor of the \emph{null} hypothesis. This argument is relevant only \strong{for
parametric test} (Default: \code{TRUE}).}
    \item{\code{results.subtitle}}{Decides whether the results of statistical tests are
to be displayed as a subtitle (Default: \code{TRUE}). If set to \code{FALSE}, only
the plot will be returned.}
    \item{\code{xlab}}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}
    \item{\code{ylab}}{Labels for \code{x} and \code{y} axis variables. If \code{NULL} (default),
variable names for \code{x} and \code{y} will be used.}
    \item{\code{subtitle}}{The text for the plot subtitle. Will work only if
\code{results.subtitle = FALSE}.}
    \item{\code{caption}}{The text for the plot caption.}
    \item{\code{ggtheme}}{A \code{{ggplot2}} theme. Default value is
\code{ggstatsplot::theme_ggstatsplot()}. Any of the \code{{ggplot2}} themes (e.g.,
\code{theme_bw()}), or themes from extension packages are allowed
(e.g., \code{ggthemes::theme_fivethirtyeight()}, \code{hrbrthemes::theme_ipsum_ps()},
etc.).}
    \item{\code{ggplot.component}}{A \code{ggplot} component to be added to the plot prepared
by \code{{ggstatsplot}}. This argument is primarily helpful for \code{grouped_}
variants of all primary functions. Default is \code{NULL}. The argument should
be entered as a \code{{ggplot2}} function or a list of \code{{ggplot2}} functions.}
  }}

\item{grouping.var}{A single grouping variable.}

\item{output}{Character that describes what is to be returned: can be
\code{"plot"} (default) or \code{"subtitle"} or \code{"caption"}. Setting this to
\code{"subtitle"} will return the expression containing statistical results. If
you have set \code{results.subtitle = FALSE}, then this will return a \code{NULL}.
Setting this to \code{"caption"} will return the expression containing details
about Bayes Factor analysis, but valid only when \code{type = "parametric"} and
\code{bf.message = TRUE}, otherwise this will return a \code{NULL}.}

\item{plotgrid.args}{A \code{list} of additional arguments passed to
\code{patchwork::wrap_plots}, except for \code{guides} argument which is already
separately specified here.}

\item{annotation.args}{A \code{list} of additional arguments passed to
\code{patchwork::plot_annotation}.}
}
\description{
Grouped scatterplots from \code{{ggplot2}} combined with marginal distribution
plots with statistical details added as a subtitle.
}
\details{
For details, see:
\url{https://indrajeetpatil.github.io/ggstatsplot/articles/web_only/ggscatterstats.html}
}
\examples{
# to ensure reproducibility
set.seed(123)
library(ggstatsplot)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)

# basic function call
grouped_ggscatterstats(
  data             = filter(movies_long, genre == "Comedy" | genre == "Drama"),
  x                = length,
  y                = rating,
  type             = "robust",
  grouping.var     = genre,
  ggplot.component = list(geom_rug(sides = "b"))
)

# using labeling
# (also show how to modify basic plot from within function call)
grouped_ggscatterstats(
  data             = filter(ggplot2::mpg, cyl != 5),
  x                = displ,
  y                = hwy,
  grouping.var     = cyl,
  type             = "robust",
  label.var        = manufacturer,
  label.expression = hwy > 25 & displ > 2.5,
  ggplot.component = scale_y_continuous(sec.axis = dup_axis())
)

# labeling without expression
grouped_ggscatterstats(
  data            = filter(movies_long, rating == 7, genre \%in\% c("Drama", "Comedy")),
  x               = budget,
  y               = length,
  grouping.var    = genre,
  bf.message      = FALSE,
  label.var       = "title",
  annotation.args = list(tag_levels = "a")
)
}
\seealso{
\code{\link{ggscatterstats}}, \code{\link{ggcorrmat}},
\code{\link{grouped_ggcorrmat}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Titanic_full}
\alias{Titanic_full}
\title{Titanic dataset.}
\format{
A data frame with 2201 rows and 5 variables
\itemize{
\item id. Dummy identity number for each person.
\item Class. 1st, 2nd, 3rd, Crew.
\item Sex. Male, Female.
\item Age. Child, Adult.
\item Survived. No, Yes.
}
}
\usage{
Titanic_full
}
\description{
Titanic dataset.
}
\details{
This data set provides information on the fate of passengers on the
fatal maiden voyage of the ocean liner 'Titanic', summarized according to
economic status (class), sex, age and survival.

This is a modified dataset from \code{datasets} package.
}
\examples{
dim(Titanic_full)
head(Titanic_full)
dplyr::glimpse(Titanic_full)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{theme_ggstatsplot}
\alias{theme_ggstatsplot}
\title{Default theme used in \code{{ggstatsplot}}}
\usage{
theme_ggstatsplot()
}
\value{
A \code{ggplot} object with the \code{theme_ggstatsplot} theme overlaid.
}
\description{
Common theme used across all plots generated in \code{{ggstatsplot}} and \emph{assumed}
by the author to be aesthetically pleasing to the user/reader. The theme is a
wrapper around \code{theme_bw()}.
}
\examples{
library(ggplot2)
library(ggstatsplot)

ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  theme_ggstatsplot()
}
