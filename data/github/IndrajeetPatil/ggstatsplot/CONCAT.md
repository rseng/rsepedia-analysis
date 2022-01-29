
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

*Wow, no problems at all. :)**Wow, no problems at all. :)*