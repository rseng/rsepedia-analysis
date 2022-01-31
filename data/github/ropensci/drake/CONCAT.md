---
title: "The drake R package: a pipeline toolkit for reproducibility and high-performance computing"
tags:
  - R
  - reproducibility
  - high-performance computing
  - pipeline
  - workflow
  - Make
authors:
  - name: William Michael Landau
    orcid: 0000-0003-1878-3253
    email: will.landau@gmail.com
    affiliation: 1
affiliations:
  - name: Eli Lilly and Company
    index: 1
date: 4 January 2018
bibliography: paper.bib
---

# Summary

The [drake](https://github.com/ropensci/drake) R package [@drake] is a workflow manager and computational engine for data science projects. Its primary objective is to keep results up to date with the underlying code and data. When it runs a project, [drake](https://github.com/ropensci/drake) detects any pre-existing output and refreshes the pieces that are outdated or missing. Not every runthrough starts from scratch, and the final answers are reproducible. With a user-friendly R-focused interface, [comprehensive documentation](https://ropensci.github.io/drake), and extensive implicit parallel computing support, [drake](https://github.com/ropensci/drake) surpasses the analogous functionality in similar tools such as [Make](www.gnu.org/software/make/) [@Make], [remake](https://github.com/richfitz/remake) [@remake], [memoise](https://github.com/r-lib/memoise) [@memoise], and [knitr](https://github.com/yihui/knitr) [@knitr].

In reproducible research, [drake](https://github.com/ropensci/drake)'s role is to provide tangible evidence that a project's results are re-creatable. [drake](https://github.com/ropensci/drake) quickly detects when the code, data, and output are synchronized. In other words, [drake](https://github.com/ropensci/drake) helps determine if the starting materials would produce the expected output if the project were to start over and run from scratch. This approach decreases the time and effort it takes to evaluate research projects for reproducibility.

Regarding high-performance computing, [drake](https://github.com/ropensci/drake) interfaces with a variety of technologies and scheduling algorithms to deploy the steps of a data analysis project. Here, the parallel computing is implicit. In other words, [drake](https://github.com/ropensci/drake) constructs the directed acyclic network of the workflow and determines which steps can run simultaneously and which need to wait for dependencies. This automation eases the cognitive and computational burdens on the user, enhancing the readability of code and thus reproducibility.

# References

<!-- README.md is generated from README.Rmd. Please edit that file -->

<center>
<img src="https://docs.ropensci.org/drake/reference/figures/infographic.svg" alt="infographic" align="center" style = "border: none; float: center;">
</center>
<table class="table">
<thead>
<tr class="header">
<th align="left">
Usage
</th>
<th align="left">
Release
</th>
<th align="left">
Development
</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">
<a href="https://www.gnu.org/licenses/gpl-3.0.en.html"><img src="https://img.shields.io/badge/licence-GPL--3-blue.svg" alt="Licence"></a>
</td>
<td align="left">
<a href="https://cran.r-project.org/package=drake"><img src="https://www.r-pkg.org/badges/version/drake" alt="CRAN"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/drake/actions?query=workflow%3Acheck"><img src="https://github.com/ropensci/drake/workflows/check/badge.svg" alt="check"></a>
</td>
</tr>
<tr class="even">
<td align="left">
<a href="https://cran.r-project.org/"><img src="https://img.shields.io/badge/R%3E%3D-3.3.0-blue.svg" alt="minimal R version"></a>
</td>
<td align="left">
<a href="https://cran.r-project.org/web/checks/check_results_drake.html"><img src="https://cranchecks.info/badges/summary/drake" alt="cran-checks"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/drake/actions?query=workflow%3Alint"><img src="https://github.com/ropensci/drake/workflows/lint/badge.svg" alt="lint"></a>
</td>
</tr>
<tr class="odd">
<td align="left">
<a href="https://CRAN.R-project.org/package=drake"><img src="https://tinyverse.netlify.com/badge/drake"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/software-review/issues/156"><img src="https://badges.ropensci.org/156_status.svg" alt="rOpenSci"></a>
</td>
<td align="left">
<a href="https://codecov.io/github/ropensci/drake?branch=main"><img src="https://codecov.io/github/ropensci/drake/coverage.svg?branch=main" alt="Codecov"></a>
</td>
</tr>
<tr class="even">
<td align="left">
<a href="https://CRAN.R-project.org/package=drake"><img src="https://cranlogs.r-pkg.org/badges/drake" alt="downloads"></a>
</td>
<td align="left">
<a href="https://doi.org/10.21105/joss.00550"><img src="https://joss.theoj.org/papers/10.21105/joss.00550/status.svg" alt="JOSS"></a>
</td>
<td align="left">
<a href="https://bestpractices.coreinfrastructure.org/projects/2135"><img src="https://bestpractices.coreinfrastructure.org/projects/2135/badge"></a>
</td>
</tr>
<tr class="odd">
<td align="left">
</td>
<td align="left">
<a href="https://zenodo.org/badge/latestdoi/82609103"><img src="https://zenodo.org/badge/82609103.svg" alt="Zenodo"></a>
</td>
<td align="left">
<a href="https://www.tidyverse.org/lifecycle/#superseded"><img src="https://img.shields.io/badge/lifecycle-superseded-blue.svg" alt='superseded lifecycle'></a>
</td>
</tr>
</tbody>
</table>
<br>

# drake is superseded. Consider targets instead.

As of 2021-01-21, `drake` is [superseded](https://www.tidyverse.org/lifecycle/#superseded). The [`targets`](https://docs.ropensci.org/targets/) R package is the long-term successor of `drake`, and it is more robust and easier to use. Please visit <https://books.ropensci.org/targets/drake.html> for full context and advice on transitioning.

# The drake R package <img src="https://docs.ropensci.org/drake/reference/figures/logo.svg" align="right" alt="logo" width="120" height = "139" style = "border: none; float: right;">

Data analysis can be slow. A round of scientific computation can take
several minutes, hours, or even days to complete. After it finishes, if
you update your code or data, your hard-earned results may no longer be
valid. How much of that valuable output can you keep, and how much do
you need to update? How much runtime must you endure all over again?

For projects in R, the `drake` package can help. It [analyzes your
workflow](https://books.ropensci.org/drake/plans.html), skips steps with
up-to-date results, and orchestrates the rest with [optional distributed
computing](https://books.ropensci.org/drake/hpc.html). At the end,
`drake` provides evidence that your results match the underlying code
and data, which increases your ability to trust your research.

# Video

## That Feeling of Workflowing (Miles McBain)

<center>

<a href="https://www.youtube.com/embed/jU1Zv21GvT4">
<img src="https://docs.ropensci.org/drake/reference/figures/workflowing.png" alt="workflowing" align="center" style = "border: none; float: center;">
</a>

</center>

(By [Miles McBain](https://github.com/MilesMcBain);
[venue](https://nyhackr.org/index.html),
[resources](https://github.com/MilesMcBain/nycr_meetup_talk))

## rOpenSci Community Call

<center>

<a href="https://ropensci.org/commcalls/2019-09-24/">
<img src="https://docs.ropensci.org/drake/reference/figures/commcall.png" alt="commcall" align="center" style = "border: none; float: center;">
</a>

</center>

([resources](https://ropensci.org/commcalls/2019-09-24/))

# What gets done stays done.

Too many data science projects follow a [Sisyphean
loop](https://en.wikipedia.org/wiki/Sisyphus):

1.  Launch the code.
2.  Wait while it runs.
3.  Discover an issue.
4.  Rerun from scratch.

For projects with long runtimes, this process gets tedious. But with
`drake`, you can automatically

1.  Launch the parts that changed since last time.
2.  Skip the rest.

# How it works

To set up a project, load your packages,

``` r
library(drake)
library(dplyr)
library(ggplot2)
library(tidyr)
#> 
#> Attaching package: 'tidyr'
#> The following objects are masked from 'package:drake':
#> 
#>     expand, gather
```

load your custom functions,

``` r
create_plot <- function(data) {
  ggplot(data) +
    geom_histogram(aes(x = Ozone)) +
    theme_gray(24)
}
```

check any supporting files (optional),

``` r
# Get the files with drake_example("main").
file.exists("raw_data.xlsx")
#> [1] TRUE
file.exists("report.Rmd")
#> [1] TRUE
```

and plan what you are going to do.

``` r
plan <- drake_plan(
  raw_data = readxl::read_excel(file_in("raw_data.xlsx")),
  data = raw_data %>%
    mutate(Ozone = replace_na(Ozone, mean(Ozone, na.rm = TRUE))),
  hist = create_plot(data),
  fit = lm(Ozone ~ Wind + Temp, data),
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE
  )
)

plan
#> # A tibble: 5 x 2
#>   target   command                                                              
#>   <chr>    <expr_lst>                                                           
#> 1 raw_data readxl::read_excel(file_in("raw_data.xlsx"))                        …
#> 2 data     raw_data %>% mutate(Ozone = replace_na(Ozone, mean(Ozone, na.rm = TR…
#> 3 hist     create_plot(data)                                                   …
#> 4 fit      lm(Ozone ~ Wind + Temp, data)                                       …
#> 5 report   rmarkdown::render(knitr_in("report.Rmd"), output_file = file_out("re…
```

So far, we have just been setting the stage. Use `make()` or
[`r_make()`](https://books.ropensci.org/drake/projects.html#safer-interactivity)
to do the real work. Targets are built in the correct order regardless
of the row order of `plan`.

``` r
make(plan) # See also r_make().
#> ▶ target raw_data
#> ▶ target data
#> ▶ target fit
#> ▶ target hist
#> ▶ target report
```

Except for files like `report.html`, your output is stored in a hidden
`.drake/` folder. Reading it back is easy.

``` r
readd(data) # See also loadd().
#> # A tibble: 153 x 6
#>    Ozone Solar.R  Wind  Temp Month   Day
#>    <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1  41       190   7.4    67     5     1
#>  2  36       118   8      72     5     2
#>  3  12       149  12.6    74     5     3
#>  4  18       313  11.5    62     5     4
#>  5  42.1      NA  14.3    56     5     5
#>  6  28        NA  14.9    66     5     6
#>  7  23       299   8.6    65     5     7
#>  8  19        99  13.8    59     5     8
#>  9   8        19  20.1    61     5     9
#> 10  42.1     194   8.6    69     5    10
#> # … with 143 more rows
```

You may look back on your work and see room for improvement, but it’s
all good\! The whole point of `drake` is to help you go back and change
things quickly and painlessly. For example, we forgot to give our
histogram a bin width.

``` r
readd(hist)
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](man/figures/unnamed-chunk-9-1.png)<!-- -->

So let’s fix the plotting function.

``` r
create_plot <- function(data) {
  ggplot(data) +
    geom_histogram(aes(x = Ozone), binwidth = 10) +
    theme_gray(24)
}
```

`drake` knows which results are affected.

``` r
vis_drake_graph(plan) # See also r_vis_drake_graph().
```

<img src="https://docs.ropensci.org/drake/reference/figures/graph.png" alt="hist1" align="center" style = "border: none; float: center;" width = "600px">

The next `make()` just builds `hist` and `report.html`. No point in
wasting time on the data or model.

``` r
make(plan) # See also r_make().
#> ▶ target hist
#> ▶ target report
```

``` r
loadd(hist)
hist
```

![](man/figures/unnamed-chunk-13-1.png)<!-- -->

# Reproducibility with confidence

The R community emphasizes reproducibility. Traditional themes include
[scientific
replicability](https://en.wikipedia.org/wiki/Replication_crisis),
literate programming with [knitr](https://yihui.name/knitr/), and
version control with
[git](https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control).
But internal consistency is important too. Reproducibility carries the
promise that your output matches the code and data you say you used.
With the exception of [non-default
triggers](https://books.ropensci.org/drake/triggers.html) and [hasty
mode](https://books.ropensci.org/drake/hpc.html#hasty-mode), `drake`
strives to keep this promise.

## Evidence

Suppose you are reviewing someone else’s data analysis project for
reproducibility. You scrutinize it carefully, checking that the datasets
are available and the documentation is thorough. But could you re-create
the results without the help of the original author? With `drake`, it is
quick and easy to find out.

``` r
make(plan) # See also r_make().
#> ℹ unloading 1 targets from environment
#> ✓ All targets are already up to date.

outdated(plan) # See also r_outdated().
#> character(0)
```

With everything already up to date, you have **tangible evidence** of
reproducibility. Even though you did not re-create the results, you know
the results are recreatable. They **faithfully show** what the code is
producing. Given the right [package
environment](https://rstudio.github.io/packrat/) and [system
configuration](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/sessionInfo.html),
you have everything you need to reproduce all the output by yourself.

## Ease

When it comes time to actually rerun the entire project, you have much
more confidence. Starting over from scratch is trivially easy.

``` r
clean()    # Remove the original author's results.
make(plan) # Independently re-create the results from the code and input data.
#> ▶ target raw_data
#> ▶ target data
#> ▶ target fit
#> ▶ target hist
#> ▶ target report
```

## Big data efficiency

Select specialized data formats to increase speed and reduce memory
consumption. In version 7.5.2.9000 and above, the available formats are
[“fst”](https://github.com/fstpackage/fst) for data frames (example
below) and “keras” for [Keras](https://keras.rstudio.com/) models
([example here](https://books.ropensci.org/drake/churn.html#plan)).

``` r
library(drake)
n <- 1e8 # Each target is 1.6 GB in memory.
plan <- drake_plan(
  data_fst = target(
    data.frame(x = runif(n), y = runif(n)),
    format = "fst"
  ),
  data_old = data.frame(x = runif(n), y = runif(n))
)
make(plan)
#> target data_fst
#> target data_old
build_times(type = "build")
#> # A tibble: 2 x 4
#>   target   elapsed              user                 system    
#>   <chr>    <Duration>           <Duration>           <Duration>
#> 1 data_fst 13.93s               37.562s              7.954s    
#> 2 data_old 184s (~3.07 minutes) 177s (~2.95 minutes) 4.157s
```

## History and provenance

As of version 7.5.2, `drake` tracks the history and provenance of your
targets: what you built, when you built it, how you built it, the
arguments you used in your function calls, and how to get the data back.
(Disable with `make(history = FALSE)`)

``` r
history <- drake_history(analyze = TRUE)
history
#> # A tibble: 12 x 11
#>    target current built exists hash  command   seed runtime na.rm quiet
#>    <chr>  <lgl>   <chr> <lgl>  <chr> <chr>    <int>   <dbl> <lgl> <lgl>
#>  1 data   TRUE    2020… TRUE   11e2… "raw_d… 1.29e9 0.011   TRUE  NA   
#>  2 data   TRUE    2020… TRUE   11e2… "raw_d… 1.29e9 0.00400 TRUE  NA   
#>  3 fit    TRUE    2020… TRUE   3c87… "lm(Oz… 1.11e9 0.006   NA    NA   
#>  4 fit    TRUE    2020… TRUE   3c87… "lm(Oz… 1.11e9 0.002   NA    NA   
#>  5 hist   FALSE   2020… TRUE   88ae… "creat… 2.10e8 0.011   NA    NA   
#>  6 hist   TRUE    2020… TRUE   0304… "creat… 2.10e8 0.003   NA    NA   
#>  7 hist   TRUE    2020… TRUE   0304… "creat… 2.10e8 0.009   NA    NA   
#>  8 raw_d… TRUE    2020… TRUE   855d… "readx… 1.20e9 0.02    NA    NA   
#>  9 raw_d… TRUE    2020… TRUE   855d… "readx… 1.20e9 0.0330  NA    NA   
#> 10 report TRUE    2020… TRUE   5504… "rmark… 1.30e9 1.31    NA    TRUE 
#> 11 report TRUE    2020… TRUE   5504… "rmark… 1.30e9 0.413   NA    TRUE 
#> 12 report TRUE    2020… TRUE   5504… "rmark… 1.30e9 0.475   NA    TRUE 
#> # … with 1 more variable: output_file <chr>
```

Remarks:

  - The `quiet` column appears above because one of the `drake_plan()`
    commands has `knit(quiet = TRUE)`.
  - The `hash` column identifies all the previous versions of your
    targets. As long as `exists` is `TRUE`, you can recover old data.
  - Advanced: if you use `make(cache_log_file = TRUE)` and put the cache
    log file under version control, you can match the hashes from
    `drake_history()` with the `git` commit history of your code.

Let’s use the history to recover the oldest histogram.

``` r
hash <- history %>%
  filter(target == "hist") %>%
  pull(hash) %>%
  head(n = 1)
cache <- drake_cache()
cache$get_value(hash)
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](man/figures/unnamed-chunk-18-1.png)<!-- -->

## Independent replication

With even more evidence and confidence, you can invest the time to
independently replicate the original code base if necessary. Up until
this point, you relied on basic `drake` functions such as `make()`, so
you may not have needed to peek at any substantive author-defined code
in advance. In that case, you can stay usefully ignorant as you
reimplement the original author’s methodology. In other words, `drake`
could potentially improve the integrity of independent replication.

## Readability and transparency

Ideally, independent observers should be able to read your code and
understand it. `drake` helps in several ways.

  - The [drake
    plan](https://docs.ropensci.org/drake/reference/drake_plan.html)
    explicitly outlines the steps of the analysis, and
    [`vis_drake_graph()`](https://docs.ropensci.org/drake/reference/vis_drake_graph.html)
    visualizes how those steps depend on each other.
  - `drake` takes care of the parallel scheduling and high-performance
    computing (HPC) for you. That means the HPC code is no longer
    tangled up with the code that actually expresses your ideas.
  - You can [generate large collections of
    targets](https://books.ropensci.org/drake/gsp.html) without
    necessarily changing your code base of imported functions, another
    nice separation between the concepts and the execution of your
    workflow

# Scale up and out.

Not every project can complete in a single R session on your laptop.
Some projects need more speed or computing power. Some require a few
local processor cores, and some need large high-performance computing
systems. But parallel computing is hard. Your tables and figures depend
on your analysis results, and your analyses depend on your datasets, so
some tasks must finish before others even begin. `drake` knows what to
do. Parallelism is implicit and automatic. See the [high-performance
computing guide](https://books.ropensci.org/drake/hpc.html) for all the
details.

``` r
# Use the spare cores on your local machine.
make(plan, jobs = 4)

# Or scale up to a supercomputer.
drake_hpc_template_file("slurm_clustermq.tmpl") # https://slurm.schedmd.com/
options(
  clustermq.scheduler = "clustermq",
  clustermq.template = "slurm_clustermq.tmpl"
)
make(plan, parallelism = "clustermq", jobs = 4)
```

# With Docker

`drake` and Docker are compatible and complementary. Here are some
examples that run `drake` inside a Docker image.

  - [`drake-gitlab-docker-example`](https://gitlab.com/ecohealthalliance/drake-gitlab-docker-example):
    A small pedagogical example workflow that leverages `drake`, Docker,
    GitLab, and continuous integration in a reproducible analysis
    pipeline. Created by [Noam Ross](https://www.noamross.net/).
  - [`pleurosoriopsis`](https://github.com/joelnitta/pleurosoriopsis):
    The workflow that supports [Ebihara *et al.* 2019. “Growth Dynamics
    of the Independent Gametophytes of *Pleurorosiopsis makinoi*
    (Polypodiaceae)” *Bulletin of the National Science Museum Series B
    (Botany)*
    45:77-86.](https://www.kahaku.go.jp/research/publication/botany.html).
    Created by [Joel Nitta](https://github.com/joelnitta).

Alternatively, it is possible to run `drake` outside Docker and use the
[`future`](https://github.com/HenrikBengtsson/future) package to send
targets to a Docker image. `drake`’s
[`Docker-psock`](https://github.com/wlandau/drake-examples/tree/main/Docker-psock)
example demonstrates how. Download the code with
`drake_example("Docker-psock")`.

# Installation

You can choose among different versions of `drake`. The CRAN release
often lags behind the [online manual](https://books.ropensci.org/drake/)
but may have fewer bugs.

``` r
# Install the latest stable release from CRAN.
install.packages("drake")

# Alternatively, install the development version from GitHub.
install.packages("devtools")
library(devtools)
install_github("ropensci/drake")
```

# Function reference

The [reference
section](https://docs.ropensci.org/drake/reference/index.html) lists all
the available functions. Here are the most important ones.

  - `drake_plan()`: create a workflow data frame (like `my_plan`).
  - `make()`: build your project.
  - `drake_history()`: show what you built, when you built it, and the
    function arguments you used.
  - `r_make()`: launch a fresh
    [`callr::r()`](https://github.com/r-lib/callr) process to build your
    project. Called from an interactive R session, `r_make()` is more
    reproducible than `make()`.
  - `loadd()`: load one or more built targets into your R session.
  - `readd()`: read and return a built target.
  - `vis_drake_graph()`: show an interactive visual network
    representation of your workflow.
  - `recoverable()`: Which targets can we salvage using `make(recover =
    TRUE)` (experimental).
  - `outdated()`: see which targets will be built in the next `make()`.
  - `deps_code()`: check the dependencies of a command or function.
  - `drake_failed()`: list the targets that failed to build in the last
    `make()`.
  - `diagnose()`: return the full context of a build, including errors,
    warnings, and messages.

# Documentation

## Core concepts

The following resources explain what `drake` can do and how it works.
The [`learndrake`](https://github.com/wlandau/learndrake) workshop
devotes particular attention to `drake`’s mental model.

  - The [user manual](https://books.ropensci.org/drake/)
  - [`drakeplanner`](https://github.com/wlandau/drakeplanner), an
    R/Shiny app to help learn `drake` and create new projects. Run
    locally with `drakeplanner::drakeplanner()` or access it at
    <https://wlandau.shinyapps.io/drakeplanner>.
  - [`learndrake`](https://github.com/wlandau/learndrake), an R package
    for teaching an extended `drake` workshop. It contains notebooks,
    slides, Shiny apps, the latter two of which are publicly deployed.
    See the
    [README](https://github.com/wlandau/learndrake/blob/main/README.md)
    for instructions and links.

## In practice

  - [Miles McBain](https://github.com/MilesMcBain)’s [excellent blog
    post](https://milesmcbain.xyz/the-drake-post/) explains the
    motivating factors and practical issues {drake} solves for most
    projects, how to set up a project as quickly and painlessly as
    possible, and how to overcome common obstacles.
  - Miles’ [`dflow`](https://github.com/MilesMcBain/dflow) package
    generates the file structure for a boilerplate `drake` project. It
    is a more thorough alternative to `drake::use_drake()`.
  - `drake` is heavily function-oriented by design, and Miles’
    [`fnmate`](https://github.com/MilesMcBain/fnmate) package
    automatically generates boilerplate code and docstrings for
    functions you mention in `drake` plans.

## Reference

  - The [reference website](https://docs.ropensci.org/drake/).
  - The [official repository of example
    code](https://github.com/wlandau/drake-examples). Download an
    example workflow from here with `drake_example()`.
  - Presentations and workshops by [Will
    Landau](https://github.com/wlandau), [Kirill
    Müller](https://github.com/krlmlr), [Amanda
    Dobbyn](https://github.com/aedobbyn), [Karthik
    Ram](https://github.com/karthik), [Sina
    Rüeger](https://github.com/sinarueeger), [Christine
    Stawitz](https://github.com/cstawitz), and others. See specific
    links at <https://books.ropensci.org/drake/index.html#presentations>
  - The [FAQ page](https://books.ropensci.org/drake/faq.html), which
    links to [appropriately-labeled issues on
    GitHub](https://github.com/ropensci/drake/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3A%22frequently+asked+question%22+).

## Use cases

The official [rOpenSci use
cases](https://discuss.ropensci.org/c/usecases) and [associated
discussion threads](https://discuss.ropensci.org/c/usecases) describe
applications of `drake` in the real world. Many of these use cases are
linked from the [`drake` tag on the rOpenSci discussion
forum](https://discuss.ropensci.org/tag/drake).

Here are some additional applications of `drake` in real-world projects.

  - [efcaguab/demografia-del-voto](https://github.com/efcaguab/demografia-del-voto)
  - [efcaguab/great-white-shark-nsw](https://github.com/efcaguab/great-white-shark-nsw)
  - [IndianaCHE/Detailed-SSP-Reports](https://github.com/IndianaCHE/Detailed-SSP-Reports)
  - [joelnitta/pleurosoriopsis](https://github.com/joelnitta/pleurosoriopsis)
  - [pat-s/pathogen-modeling](https://github.com/pat-s/pathogen-modeling)
  - [sol-eng/tensorflow-w-r](https://github.com/sol-eng/tensorflow-w-r)
  - [tiernanmartin/home-and-hope](https://github.com/tiernanmartin/home-and-hope)

## `drake` projects as R packages

Some folks like to structure their `drake` workflows as R packages.
Examples are below. In your own analysis packages, be sure to call
`drake::expose_imports(yourPackage)` so `drake` can watch you package’s
functions for changes and rebuild downstream targets accordingly.

  - [b-rodrigues/coolmlproject](https://github.com/b-rodrigues/coolmlproject)
  - [tiernanmartin/drakepkg](https://github.com/tiernanmartin/drakepkg)

# Help and troubleshooting

The following resources document many known issues and challenges.

  - [Frequently-asked
    questions](https://github.com/ropensci/drake/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3A%22frequently+asked+question%22+).
  - [Debugging and testing drake
    projects](https://books.ropensci.org/drake/debugging.html)
  - [Other known issues](https://github.com/ropensci/drake/issues)
    (please search both open and closed ones).

If you are still having trouble, please submit a [new
issue](https://github.com/ropensci/drake/issues/new) with a bug report
or feature request, along with a minimal reproducible example where
appropriate.

The GitHub issue tracker is mainly intended for bug reports and feature
requests. While questions about usage etc. are also highly encouraged,
you may alternatively wish to post to [Stack
Overflow](https://stackoverflow.com) and use the [`drake-r-package`
tag](https://stackoverflow.com/tags/drake-r-package).

# Contributing

Development is a community effort, and we encourage participation.
Please read
[CONTRIBUTING.md](https://github.com/ropensci/drake/blob/main/CONTRIBUTING.md)
for details.

# Similar work

`drake` enhances reproducibility and high-performance computing, but not
in all respects. [Literate programming](https://rmarkdown.rstudio.com/),
[local library managers](https://rstudio.github.io/packrat),
[containerization](https://www.docker.com/), and [strict session
managers](https://github.com/tidyverse/reprex) offer more robust
solutions in their respective domains. And for the problems `drake`
*does* solve, it stands on the shoulders of the giants that came before.

## Pipeline tools

### GNU Make

The original idea of a time-saving reproducible build system extends
back at least as far as [GNU Make](https://www.gnu.org/software/make/),
which still aids the work of [data
scientists](http://blog.kaggle.com/2012/10/15/make-for-data-scientists/)
as well as the original user base of complied language programmers. In
fact, the name “drake” stands for “Data Frames in R for Make”.
[Make](https://kbroman.org/minimal_make/) is used widely in reproducible
research. Below are some examples from [Karl Broman’s
website](https://kbroman.org/minimal_make/).

  - Bostock, Mike (2013). “A map of flowlines from NHDPlus.”
    <https://github.com/mbostock/us-rivers>. Powered by the Makefile at
    <https://github.com/mbostock/us-rivers/blob/master/Makefile>.
  - Broman, Karl W (2012). “Halotype Probabilities in Advanced
    Intercross Populations.” *G3* 2(2), 199-202.Powered by the
    `Makefile` at
    <https://github.com/kbroman/ailProbPaper/blob/master/Makefile>.
  - Broman, Karl W (2012). “Genotype Probabilities at Intermediate
    Generations in the Construction of Recombinant Inbred Lines.”
    \*Genetics 190(2), 403-412. Powered by the Makefile at
    <https://github.com/kbroman/preCCProbPaper/blob/master/Makefile>.
  - Broman, Karl W and Kim, Sungjin and Sen, Saunak and Ane, Cecile and
    Payseur, Bret A (2012). “Mapping Quantitative Trait Loci onto a
    Phylogenetic Tree.” *Genetics* 192(2), 267-279. Powered by the
    `Makefile` at
    <https://github.com/kbroman/phyloQTLpaper/blob/master/Makefile>.

Whereas [GNU Make](https://www.gnu.org/software/make/) is
language-agnostic, `drake` is fundamentally designed for R.

  - Instead of a
    [Makefile](https://github.com/kbroman/preCCProbPaper/blob/master/Makefile),
    `drake` supports an R-friendly [domain-specific
    language](https://books.ropensci.org/drake/plans.html#large-plans)
    for declaring targets.
  - Targets in [GNU Make](https://www.gnu.org/software/make/) are files,
    whereas targets in `drake` are arbitrary variables in memory.
    (`drake` does have opt-in support for files via `file_out()`,
    `file_in()`, and `knitr_in()`.) `drake` caches these objects in its
    own [storage system](https://github.com/richfitz/storr) so R users
    rarely have to think about output files.

### Remake

[remake](https://github.com/richfitz/remake) itself is no longer
maintained, but its founding design goals and principles live on through
[drake](https://github.com/ropensci/drake). In fact,
[drake](https://github.com/ropensci/drake) is a direct re-imagining of
[remake](https://github.com/richfitz/remake) with enhanced scalability,
reproducibility, high-performance computing, visualization, and
documentation.

### Factual’s Drake

[Factual’s Drake](https://github.com/Factual/drake) is similar in
concept, but the development effort is completely unrelated to the
[drake R package](https://github.com/ropensci/drake).

### Other pipeline tools

There are [countless other successful pipeline
toolkits](https://github.com/pditommaso/awesome-pipeline). The `drake`
package distinguishes itself with its R-focused approach,
Tidyverse-friendly interface, and a [thorough selection of parallel
computing technologies and scheduling
algorithms](https://books.ropensci.org/drake/hpc.html).

## Memoization

Memoization is the strategic caching of the return values of functions.
It is a lightweight approach to the core problem that `drake` and other
pipeline tools are trying to solve. Every time a memoized function is
called with a new set of arguments, the return value is saved for future
use. Later, whenever the same function is called with the same
arguments, the previous return value is salvaged, and the function call
is skipped to save time. The
[`memoise`](https://github.com/r-lib/memoise) package is the primary
implementation of memoization in R.

Memoization saves time for small projects, but it arguably does not go
far enough for large reproducible pipelines. In reality, the return
value of a function depends not only on the function body and the
arguments, but also on any nested functions and global variables, the
dependencies of those dependencies, and so on upstream. `drake` tracks
this deeper context, while [memoise](https://github.com/r-lib/memoise)
does not.

## Literate programming

[Literate programming](https://rmarkdown.rstudio.com/) is the practice
of narrating code in plain vernacular. The goal is to communicate the
research process clearly, transparently, and reproducibly. Whereas
commented code is still mostly code, literate
[knitr](https://yihui.name/knitr/) / [R
Markdown](https://rmarkdown.rstudio.com/) reports can become websites,
presentation slides, lecture notes, serious scientific manuscripts, and
even books.

### knitr and R Markdown

`drake` and [knitr](https://yihui.name/knitr/) are symbiotic. `drake`’s
job is to manage large computation and orchestrate the demanding tasks
of a complex data analysis pipeline.
[knitr](https://yihui.name/knitr/)’s job is to communicate those
expensive results after `drake` computes them.
[knitr](https://yihui.name/knitr/) / [R
Markdown](https://rmarkdown.rstudio.com/) reports are small pieces of an
overarching `drake` pipeline. They should focus on communication, and
they should do as little computation as possible.

To insert a [knitr](https://yihui.name/knitr/) report in a `drake`
pipeline, use the `knitr_in()` function inside your [`drake`
plan](https://books.ropensci.org/drake/plans.html), and use `loadd()`
and `readd()` to refer to targets in the report itself. See an [example
here](https://github.com/wlandau/drake-examples/tree/main/main).

### Version control

`drake` is not a version control tool. However, it is fully compatible
with [`git`](https://git-scm.com/),
[`svn`](https://en.wikipedia.org/wiki/Apache_Subversion), and similar
software. In fact, it is good practice to use
[`git`](https://git-scm.com/) alongside `drake` for reproducible
workflows.

However, data poses a challenge. The datasets created by `make()` can
get large and numerous, and it is not recommended to put the `.drake/`
cache or the `.drake_history/` logs under version control. Instead, it
is recommended to use a data storage solution such as
DropBox or [OSF](https://osf.io/ka7jv/wiki/home/).

### Containerization and R package environments

`drake` does not track R packages or system dependencies for changes.
Instead, it defers to tools like [Docker](https://www.docker.com),
[Singularity](https://sylabs.io/singularity/),
[`renv`](https://github.com/rstudio/renv), and
[`packrat`](https://github.com/rstudio/packrat), which create
self-contained portable environments to reproducibly isolate and ship
data analysis projects. `drake` is fully compatible with these tools.

### workflowr

The [`workflowr`](https://github.com/jdblischak/workflowr) package is a
project manager that focuses on literate programming, sharing over the
web, file organization, and version control. Its brand of
reproducibility is all about transparency, communication, and
discoverability. For an example of
[`workflowr`](https://github.com/jdblischak/workflowr) and `drake`
working together, see [this machine learning
project](https://github.com/pat-s/2019-feature-selection) by [Patrick
Schratz](https://github.com/pat-s).

# Citation

``` r
citation("drake")
#> 
#> To cite drake in publications use:
#> 
#>   William Michael Landau, (2018). The drake R package: a pipeline
#>   toolkit for reproducibility and high-performance computing. Journal
#>   of Open Source Software, 3(21), 550,
#>   https://doi.org/10.21105/joss.00550
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {The drake R package: a pipeline toolkit for reproducibility and high-performance computing},
#>     author = {William Michael Landau},
#>     journal = {Journal of Open Source Software},
#>     year = {2018},
#>     volume = {3},
#>     number = {21},
#>     url = {https://doi.org/10.21105/joss.00550},
#>   }
```

# Acknowledgements

Special thanks to [Jarad Niemi](https://www.jarad.me/), my advisor from
[graduate school](https://stat.iastate.edu/), for first introducing me
to the idea of [Makefiles](https://www.gnu.org/software/make/) for
research. He originally set me down the path that led to `drake`.

Many thanks to [Julia Lowndes](https://github.com/jules32), [Ben
Marwick](https://github.com/benmarwick), and [Peter
Slaughter](https://github.com/gothub) for [reviewing drake for
rOpenSci](https://github.com/ropensci/software-review/issues/156), and to
[Maëlle Salmon](https://github.com/maelle) for such active involvement
as the editor. Thanks also to the following people for contributing
early in development.

  - [Alex Axthelm](https://github.com/AlexAxthelm)
  - [Chan-Yub Park](https://github.com/mrchypark)
  - [Daniel Falster](https://github.com/dfalster)
  - [Eric Nantz](https://github.com/rpodcast)
  - [Henrik Bengtsson](https://github.com/HenrikBengtsson)
  - [Ian Watson](https://github.com/IanAWatson)
  - [Jasper Clarkberg](https://github.com/dapperjapper)
  - [Kendon Bell](https://github.com/kendonB)
  - [Kirill Müller](https://github.com/krlmlr)
  - [Michael Schubert](https://github.com/mschubert)

Credit for images is [attributed
here](https://github.com/ropensci/drake/blob/main/man/figures/image-credit.md).

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# Version 7.13.3.9000



# Version 7.13.3

* Improve error messages from static code analysis of malformed code (#1371, @billdenney).
* Handle invalid language objects in commands (#1372, @gorgitko).
* Do not lock namespaces (#1373, @gorgitko).
* Compatibility with `rlang` PR 1255.

# Version 7.13.2

* Update SLURM `batchtools` template file can be brewed (#1359, @pat-s).
* Change start-up message to tip about `targets`.

# Version 7.13.1

* Add files `NOTICE` and `inst/NOTICE` to more explicitly credit code included from other open source projects. (Previously `drake` just had comments in the source with links to the various projects.)

# Version 7.13.0

## Bug fixes

* Avoid checking printed output in test of testing infrastructure.
* Use `dsl_sym()` instead of `as.symbol()` when constructing commands for `combine()` (#1340, @vkehayas).

## New features

* Add a new `level_separation` argument to `vis_drake_graph()` and `render_drake_graph()` to control the aspect ratio of `visNetwork` graphs (#1303, @matthewstrasiotto, @matthiasgomolka, @robitalec).

# Version 7.12.7

## Enhancements

* Deprecate `caching = "master"` in favor of `caching = "main"`.
* Improve the error message when a valid plan is not supplied (#1334, @robitalec).

# Version 7.12.6

## Bug fixes

* Fix defunct functions error message when using namespace (#1310, @malcolmbarrett).
* Preserve names of list elements in `.data` in DSL (#1323, @shirdekel).
* Use `identical()` to compare file hashes (#1324, @shirdekel).
* Set `seed = TRUE` in `future::future()`.
* Manually relay warnings when `parallelism = "clustermq"` and `caching = "worker"` (@richardbayes).

## Enhancements

* Make logs more machine-readable by sanitizing messages and preventing race conditions (#1331, @Plebejer).

# Version 7.12.5

## Bug fixes

* Sanitize empty symbols in language columns (#1299, @odaniel1).
* Handle cases where `NROW()` throws an error (#1300, `julian-tagell` on Stack Overflow).
* Prohibit dynamic branching over non-branching dynamic files (#1302, @djbirke).

## Enhancements

* Transition to updated `lifecycle` that does not require badges to be in `man/figures`.
* Improve error message for empty dynamic grouping variables (#1308, @saadaslam).
* Expose the `log_worker` argument of `clustermq::workers()` to `make()` and `drake_config()` (#1305, @billdenney, @mschubert).
* Set `as.is` to `TRUE` in `utils::type.convert()` (#1309, @bbolker).

# Version 7.12.4

* Fix a CRAN warning about docs.

# Version 7.12.3

## Bug fixes

* `cached_planned()` and `cached_unplanned()` now work with non-standard cache locations (#1268, @Plebejer).
* Set `use_cache` to `FALSE` more often (#1257, @Plebejer).
* Use namespaced function calls in mtcars example instead of loading packages.
* Replace the `iris` dataset with the `airquality` dataset in all documentation, examples, and tests (#1271).
* Assign functions created with `code_to_function()` to the proper environment (#1275, @robitalec).
* Store tracebacks as character vectors and restrict the contents of error objects to try to prevent accidental storage of large data from the environment (#1276, @billdenney).
* Strongly depend on `tidyselect` (#1274, @dernst).
* Avoid `txtq` lockfiles (#1232, #1239, #1280, @danwwilson, @pydupont, @mattwarkentin).

## New features

* Add a new `drake_script()` function to write `_drake.R` files for `r_make()` (#1282).

## Enhancements

* Deprecate `expose_imports()` in favor of `make(envir = getNamespace("yourPackage")` (#1286, @mvarewyck).
* Suppress the message recommending `r_make()` if `getOption("drake_r_make_message")` is `FALSE` (#1238, @januz).
* Improve the appearance of the `visNetwork` graph by using the hierarchical layout with `visEdges(smooth = list(type = "cubicBezier", forceDirection = TRUE))` (#1289, @mstr3336).

# Version 7.12.2

## Bug fixes

* Invalidate old sub-targets when finalizing a dynamic target (@richardbayes). Solves a major reproducibility bug (#1260).
* Prevent `splice_inner()` from dropping formal arguments shared by `c()` (#1262, @bart1).

# Version 7.12.1

## Bug fixes

* Repair `subtarget_hashes.cross()` for crosses on a single grouping variable.
* Repair dynamic `group()` used with specialized formats (#1236, @adamaltmejd).
* Enforce `tidyselect` >= 1.0.0.

## New features

* Allow user-defined target names in static branching with the `.names` argument (#1240, @maciejmotyka, @januz).

## Enhancements

* Do not analyze dependencies of calls to `drake_plan()` (#1237, @januz).
* Error message for locked cache gives paste-able error message in Windows (#1243, @billdenney).
* Prevent stack traces from accidentally storing large amounts of data (#1253, @sclewis23).

# Version 7.12.0

## Bug fixes

* **Ensure up-to-date sub-targets are skipped even if the dynamic parent does not get a chance to finalize (#1209, #1211, @psadil, @kendonB).**
* Restrict static transforms so they only use the upstream part of the plan (#1199, #1200, @bart1).
* Correctly match the names and values of dynamic `cross()` sub-targets (#1204, @psadil). Expansion order is the same, but names are correctly matched now.
* Stop trying to remove `file_out()` files in `clean()`, even when `garbage_collection` is `TRUE` (#521, @the-Hull).
* Fix `keep_going = TRUE` for formatted targets (#1206).
* Use the correct variable names in logger helper (`progress_bar` instead of `progress`) so that `drake` works without the `progress` package (#1208, @mbaccou).
* Avoid conflict between formats and upstream dynamic targets (#1210, @psadil).
* Always compute trigger metadata up front because recovery keys need it.
* Deprecate and remove hasty mode and custom parallel backends (#1222).
* Compartmentalize fixed runtime parameters in `config$settings` (#965).

## New features

* Add new functions `drake_done()` and `drake_cancelled()` (#1205).

## Speedups

* Avoid reading build times of dynamic sub-targets in `drake_graph_info()` (#1207).

## Enhancements

* Show an empty progress bar just before targets start to build when `verbose` is `2` (#1203, @kendonB).
* Deprecate the `jobs` argument of `clean()`.
* Show an informative error message for empty dynamic grouping variables (#1212, @kendonB).
* Throw error messages if users supply dynamic targets to `drake_build()` or `drake_debug()` (#1214, @kendonB).
* Log the sub-target name and index of the failing sub-target in the metadata of the sub-target and its parent (#1214, @kendonB).
* Shorten the call stack in error metadata.
* Deprecate and remove custom schedulers (#1222).
* Deprecate `hasty_build` (#1222).
* Migrate constant runtime parameters to `config$settings` (#965).
* Warn the user if `file_in()`/`file_out()`/`knitr_in()` files are not literal strings (#1229).
* Prohibit `file_out()` and `knitr_in()` in imported functions (#1229).
* Prohibit `knitr_in()` in dynamic branching (#1229).
* Improve the help file of `target()`.
* Deprecate and rename progress functions to avoid potential name conflicts (`progress()` => `drake_progress()`, `running()` => `drake_running()`, `failed()` => `drake_failed()`) (#1205).

# Version 7.11.0

## Bug fixes

* Sanitize internal S3 classes for target storage (#1159, @rsangole).
* Bump `digest` version to require 0.6.21 (#1166, @boshek)
* Actually store output file sizes in metadata.
* Use the `depend` trigger to toggle invalidation from dynamic-only dependencies, including the `max_expand` argument of `make()`.
* Repair `session_info` argument parsing (and reduce calls to `utils::sessionInfo()` in tests).
* Ensure compatibility with `tibble` 3.0.0.

## New features

* Allow dynamic files with `target(format = "file")` (#1168, #1127).
* Implement dynamic `max_expand` on a target-by-target basis via `target()` (#1175, @kendonB).

## Enhancements

* Assert dependencies of formats at the very beginning of `make()`, not in `drake_config()` (#1156).
* In `make(verbose = 2)`, remove the spinner and use a progress bar to track how many targets are done so far.
* Reduce logging of utility functions.
* Improve the aesthetics of console messages using `cli` (optional package).
* Deprecate `console_log_file` in favor of `log_make` as an argument to `make()` and `drake_config()`.
* Immediately relay warnings and messages in `"loop"` and `"future"` parallel backends (#400).
* Warn when converting trailing dots (#1147).
* Warn about imports with trailing dots on Windows (#1147).
* Allow user-defined caches for the `loadd()` RStudio addin through the new `rstudio_drake_cache` global option (#1169, @joelnitta).
* Change dynamic target finalization message to "finalize" instead of "aggregate" (#1176, @kendonB).
* Describe the limits of `recoverable()`, e.g. dynamic branching + dynamic files.
* Throw an error instead of a warning in `drake_plan()` if a grouping variable is undefined or invalid (#1182, @kendonB).
* Rigorous S3 framework for static code analysis objects of type `drake_deps` and `drake_deps_ht` (#1183).
* Use `rlang::trace_back()` to make `diagnose()$error$calls` nicer (#1198).

# Version 7.10.0

## Unavoidable but minor breaking changes

These changes invalidate some targets in some workflows, but they are necessary bug fixes.

* Remove spurious local variables detected in `$<-()` and `@<-()` (#1144).
* Avoid target names with trailing dots (#1147, @Plebejer).

## Bug fixes

* Handle unequal list columns in `bind_plans()` (#1136, @jennysjaarda).
* Handle non-vector sub-targets in dynamic branching (#1138).
* Handle calls in `analyze_assign()` (#1119, @jennysjaarda).
* Restore correct environment locking (#1143, @kuriwaki).
* Log `"running"` progress of dynamic targets.
* Log dynamic targets as failed if a sub-target fails (#1158).

## New features

* Add a new `"fst_tbl"` format for large `tibble` targets (#1154, @kendonB).
* Add a new `format` argument to `make()`, an optional custom storage format for targets without an explicit `target(format = ...)` in the plan (#1124).
* Add a new `lock_cache` argument to `make()` to optionally suppress cache locking (#1129). (It can be annoying to interrupt `make()` repeatedly and unlock the cache manually every time.)
* Add new functions `cancel()` and `cancel_if()` function to cancel targets mid-build (#1131).
* Add a new `subtarget_list` argument to `loadd()` and `readd()` to optionally load a dynamic target as a list of sub-targets (#1139, @MilesMcBain).
* Prohibit dynamic `file_out()` (#1141).

## Enhancements

* Check for illegal formats early on at the `drake_config()` level (#1156, @MilesMcBain).
* Smoothly deprecate the `config` argument in all user-side functions (#1118, @vkehayas). Users can now supply the plan and other `make()` arguments directly, without bothering with `drake_config()`. Now, you only need to call `drake_config()` in the `_drake.R` file for `r_make()` and friends. Old code with `config` objects should still work. Affected functions:
    * `make()`
    * `outdated()`
    * `drake_build()`
    * `drake_debug()`
    * `recoverable()`
    * `missed()`
    * `deps_target()`
    * `deps_profile()`
    * `drake_graph_info()`
    * `vis_drake_graph()`
    * `sankey_drake_graph()`
    * `drake_graph()`
    * `text_drake_graph()`
    * `predict_runtime()`. Needed to rename the `targets` argument to `targets_predict` and `jobs` to `jobs_predict`.
    * `predict_workers()`. Same argument name changes as `predict_runtime()`.
* Because of #1118, the only remaining user-side purpose of `drake_config()` is to serve functions `r_make()` and friends.
* Document the limitations of grouping variables (#1128).
* Handle the `@` operator. For example, in the static code analysis of `x@y`, do not register `y` as a dependency (#1130, @famuvie).
* Remove superfluous/incorrect information about imports from the output of `deps_profile()` (#1134, @kendonB).
* Append hashes to `deps_target()` output (#1134, @kendonB).
* Add S3 class and pretty print method for `drake_meta_()` objects objects.
* Use call stacks instead of environment inheritance to power `drake_envir()` and `id_chr()` (#1132).
* Allow `drake_envir()` to select the environment with imports (#882).
* Improve visualization labels for dynamic targets: clarify that the listed runtime is a total runtime over all sub-targets and list the number of sub-targets.


# Version 7.9.0

## Breaking changes in dynamic branching

- Embrace the `vctrs` paradigm and its type stability for dynamic branching (#1105, #1106).
- Accept `target` as a symbol by default in `read_trace()`. Required for the trace to make sense in #1107.

## Bug fixes

- Repair reference to custom HPC resources in the `"future"` backend (#1083, @jennysjaarda).
- Properly copy data when importing targets from one cache into another (#1120, @brendanf).
- Prevent dynamic vector sizes from conflicting with file sizes in metadata.

## New features

- Add a new `log_build_times` argument to `make()` and `drake_config()`. Allows users to disable the recording of build times. Produces a speedup of up to 20% on Macs (#1078).
- Implement cache locking to prohibit concurrent calls to `make()`, `outdated(make_imports = TRUE)`, `recoverable(make_imports = TRUE)`, `vis_drake_graph(make_imports = TRUE)`, `clean()`, etc. on the same cache.
- Add a new `format` trigger to invalidate targets when the specialized data format changes (#1104, @kendonB).
- Add new functions `cache_planned()` and `cache_unplanned()` to help selectively clean workflows with dynamic targets (#1110, @kendonB).
- Add S3 classes and pretty print methods for `drake_config()` objects and `analyze_code()` objects.
- Add a new `"qs"` format (#1121, @kendonB).

## Speedups

- Avoid setting seeds for imports (#1086, @adamkski).
- Avoid working directly with POSIXct times (#1086, @adamkski)
- Avoid excessive calls to `%||%` (`%|||%` is faster). (#1089, @billdenney)
- Remove `%||NA` due to slowness (#1089, @billdenney).
- Use hash tables to speed up `is_dynamic()` and `is_subtarget()` (#1089, @billdenney).
- Use `getVDigest()` instead of `digest()` (#1089, #1092, https://github.com/eddelbuettel/digest/issues/139#issuecomment-561870289, @eddelbuettel, @billdenney).
- Pre-compute `backtick` and `.deparseOpts()` to speed up `deparse()` (#1086, `https://stackoverflow.com/users/516548/g-grothendieck`, @adamkski).
- Pre-compute which targets exist in advance (#1095).
- Avoid gratuitous cache interactions and data frame operations in `build_times()` (#1098).
- Use `mget_hash()` in `progress()` (#1098).
- Get target progress info only once in `drake_graph_info()` (#1098).
- Speed up the retrieval of old metadata in `outdated()` (#1098).
- In `make()`, avoid checking for nonexistent metadata for missing targets.
- Reduce logging in `drake_config()`.

## Enhancements

- Write a complete project structure in `use_drake()` (#1097, @lorenzwalthert, @tjmahr).
- Add a minor logger note to say how many dynamic sub-targets are registered at a time (#1102, @kendonB).
- Handle dependencies that are dynamic targets but not declared as such for the current target (#1107).
- Internally, the "layout" data structure is now called the "workflow specification", or "spec" for short. The spec is `drake`'s interpretation of the plan. In the plan, all the dependency relationships among targets and files are *implicit*. In the spec, they are all *explicit*. We get from the plan to the spec using static code analysis, e.g. `analyze_code()`.


# Version 7.8.0

## Bug fixes

- Prevent `drake::drake_plan(x = target(...))` from throwing an error if `drake` is not loaded (#1039, @mstr3336).
- Move the `transformations` lifecycle badge to the proper location in the docstring (#1040, @jeroen).
- Prevent `readd()` / `loadd()` from turning an imported function into a target (#1067).
- Align in-memory `disk.frame` targets with their stored values (#1077, @brendanf).

## New features

- Implement dynamic branching (#685).
- Add a new `subtargets()` function to get the cached names of the sub-targets of a dynamic target.
- Add new `subtargets` arguments to `loadd()` and `readd()` to retrieve specific sub-targets from a parent dynamic target.
- Add new `get_trace()` and `read_trace()` functions to help track which values of grouping variables go into the making of dynamic sub-targets.
- Add a new `id_chr()` function to get the name of the target while `make()` is running.
- Implement `plot(plan)` (#1036).
- `vis_drake_graph()`, `drake_graph_info()`, and `render_drake_graph()` now take arguments that allow behavior to be defined upon selection of nodes. (#1031,@mstr3336).
- Add a new `max_expand` argument to `make()` and `drake_config()` to scale down dynamic branching (#1050, @hansvancalster).

## Enhancements

- Document transformation functions in a way that avoids having to create true functions (#979).
- Avoid always invalidating the memoized layout when we set the knitr hash.
- Change the names of environments in `drake_config()` objects.
- Assert that `prework` is a language object, list of language objects, or character vector (#1 at pat-s/multicore-debugging on GitHub, @pat-s).
- Use an environment instead of a list for `config$layout`. Supports internal modifications by reference. Required for #685.
- Clean up the code of the parallel backends.
- Make `dynamic` a formal argument of `target()`.
- Always lock/unlock the environment target by target, allowing informative error messages to appear more readily (#1062, @PedramNavid)
- Automatically ignore `storr`s and decorated `storr`s (#1071).
- Speed up memory management by avoiding a call to `setdiff()` and avoiding `names(config$envir_targets)`.


# Version 7.7.0

## Bug fixes

- Take the sum instead of the max in `dir_size()`. Incurs rehashing for some workflows, but should not invalidate any targets.

## New features

- Add a new `which_clean()` function to preview which targets will be invalidated by `clean()` (#1014, @pat-s).
- Add serious import and export methods for the decorated `storr` (#1015, @billdenney, @noamross).
- Add a new `"diskframe"` format for larger-than-memory data (#1004, @xiaodaigh).
- Add a new `drake_tempfile()` function to help with `"diskframe"` format. It makes sure we are not copying large datasets across different physical storage media (#1004, @xiaodaigh).
- Add new function `code_to_function()` to allow for parsing script based workflows into functions so `drake_plan()` can begin to manage the workflow and track dependencies. (#994, @thebioengineer)

## Enhancements

- Coerce seeds to integers in `seed_trigger()` (#1013, @CreRecombinase).
- Hard wrap long labels in graph visuals (#1017).
- Nest the history `txtq` API inside decorated `storr` API (#1020).
- Reduce cyclomatic complexity of internal functions.
- Reduce retrievals of old target metadata to try to improve performance (#1027).


# Version 7.6.2

## Bug fixes

- Remove README.md from CRAN altogether. Also remove all links from the news and vignette. The links trigger too many CRAN notes, which made the automated checks too brittle.
- Serialize formats that need serialization (like "keras") before sending the data from HPC workers to the main process (#989).
- Check for custom-formatted files when checking checksums.
- Force fst-formatted targets to plain data frames. Same goes for the new "fst_dt" format.
- Change the meaning and behavior of `max_expand` in `drake_plan()`. `max_expand` is now the maximum number of targets produced by `map()`, `split()`, and `cross()`. For `cross()`, this reduces the number of targets (less cumbersome) and makes the subsample of targets more representative of the complete grid. It also. ensures consistent target naming when `.id` is `FALSE` (#1002). Note: `max_expand` is not for production workflows anyway, so this change does not break anything important. Unfortunately, we do lose the speed boost in `drake_plan()` originally due to `max_expand`, but `drake_plan()` is still fast, so that is not so bad.
- Drop specialized formats of `NULL` targets (#998).
- Prevent false grouping variables from partially tagging along in `cross()` (#1009). The same fix should apply to `map()` and `split()` too.
- Respect graph topology when recovering old grouping variables for `map()` (#1010).

## New features

- Add a new "fst_dt" format for `fst`-powered saving of `data.table` objects.
- Support a custom "caching" column of the plan to select main vs worker caching for each target individually (#988).
- Make `transform` a formal argument of `target()` so that users do not have to type "transform =" all the time in `drake_plan()` (#993).
- Migrate the documentation website from `ropensci.github.io/drake` to `docs.ropensci.org/drake`.

## Enhancements

- Document the HPC limitations of `target(format = "keras")` (#989).
- Remove the now-superfluous vignette.
- Wrap up console and text file logging functionality into a reference class (#964).
- Deprecate the `verbose` argument in various caching functions. The location of the cache is now only printed in `make()`. This made the previous feature easier to implement.
- Carry forward nested grouping variables in `combine()` (#1008).
- Improve the encapsulation of hash tables in the decorated `storr` (#968).


# Version 7.6.1

## Bug fixes

- CRAN hotfix: remove a broken link in the README.


# Version 7.6.0

## Bug fixes

- Make `drake_plan(transform = slice())` understand `.id` and grouping variables (#963).
- Repair `clean(garbage_collection = TRUE, destroy = TRUE)`. Previously it destroyed the cache before trying to collect garbage.
- Ensure that `r_make()` passes informative error messages back to the calling process (#969).
- Avoid downloading full contents of URLs when rehashing (#982)
- Retain upstream grouping variables of `map()` and `cross()` on topologically side-by-side targets (#983).
- Manually enforce the correct ordering in `dsl_left_outer_join()` so `cross()` selects the right combinations of existing targets (#986). This bug was probably introduced in the solution to #983.
- Make the output of `progress()` more consistent, less dependent on whether `tidyselect` is installed.

## New features

- Support specialized data storage via a decorated cache and `format` argument of `target()` (#971). This allows users to leverage faster ways to save and load targets, such as `write_fst()` for data frames and `save_model_hdf5()` for Keras models. It also improves memory because it prevents `storr` from making a serialized in-memory copy of large data objects.
- Add `tidyselect` functionality for `...` in `progress()`, analogous to `loadd()`, `build_times()`, and `clean()`.
- Support S3 for user-defined generics (#959). If the generic `do_stuff()` and the method `stuff.your_class()` are defined in `envir`, and if `do_stuff()` has a call to `UseMethod("stuff")`, then `drake`'s code analysis will detect `stuff.your_class()` as a dependency of `do_stuff()`.
- Add authentication support for `file_in()` URLs. Requires the new `curl_handles` argument of `make()` and `drake_config()` (#981).

## Enhancements

- Document DSL keywords as if they were true functions: `target()`, `map()`, `split()`, `cross()`, and `combine()` (#979).
- Do garbage collection between the unloading and loading phases of memory management.
- Keep `file_out()` files in `clean()` unless `garbage_collection` is `TRUE`. That way, `make(recover = TRUE)` is a true "undo button" for `clean()`. `clean(garbage_collection = TRUE)` still removes data in the cache, as well as any `file_out()` files from targets currently being cleaned.
- The menu in `clean()` only appears if `garbage_collection` is `TRUE`. Also, this menu is added to `rescue_cache(garbage_collection = TRUE)`.
- Reorganize the internal code files and functions to make development easier.
- Move the history inside the cache folder `.drake/`. The old `.drake_history/` folder was awkward. Old histories are migrated during `drake_config()`, and `drake_history()`.
- Add lifecycle badges to exported functions.

# Version 7.5.2

## Bug fixes

- Eliminate accidental creations of `.drake_history` in `plan_to_code()`, `plan_to_notebook()`, and the examples in the help files.


# Version 7.5.1

## Bug fixes

- Change \.drake_history$ to ^\.drake_history$ in .Rbuildignore appease CRAN checks. 
- Repair help file examples.


# Version 7.5.0

## New features

- Add automated data recovery (#945). This is still experimental and disabled by default. Requires `make(recover = TRUE)`.
- Add new functions `recoverable()` and `r_recoverable()` to show targets that are outdated but recoverable via `make(recover = TRUE)`.
- Track the history and provenance of targets, viewable with `drake_history()`. Powered by `txtq` (#918, #920).
- Add a new `no_deps()` function, similar to `ignore()`. `no_deps()` suppresses dependency detection but still tracks changes to the literal code (#910).
- Add a new "autoclean" memory strategy (#917).
- Export `transform_plan()`.
- Allow a custom `seed` column of `drake` plans to set custom seeds (#947).
- Add a new `seed` trigger to optionally ignore changes to the target seed (#947).

## Enhancements

- In `drake_plan()`, interpret custom columns as non-language objects (#942).
- Suggest and assert `clustermq` >= 0.8.8.
- Log the target name in a special column in the console log file (#909).
- Rename the "memory" memory strategy to "preclean" (with deprecation; #917).
- Deprecate `ensure_workers` in `drake_config()` and `make()`.
- Warn when the user supplies additional arguments to `make()` after `config` is already supplied.
- Prevent users from running `make()` from inside the cache (#927).
- Add `CITATION` file with JOSS paper.
- In `deps_profile()`, include the seed and change the names.
- Allow the user to set a different seed in `make()`. All this does is invalidate old targets.
- Use `set_hash()` and `get_hash()` in `storr` to double the speed of progress tracking.

## Bug fixes

- In the static code analysis for dependency detection, ignore list elements referenced with `$` (#938).
- Minor: handle strings embedded in language objects (#934).
- Minor: supply `xxhash64` as the default hash algorithm for non-`storr` hashing if the driver does not have a hash algorithm.


# Version 7.4.0

## Mildly breaking changes

These changes are technically breaking changes, but they should only affect advanced users.

- `rescue_cache()` no longer returns a value.

## Bug fixes

- Restore compatibility with `clustermq` (#898). Suggest version >= 0.8.8 but allow 0.8.7 as well.
- Ensure `drake` recomputes `config$layout` when `knitr` reports change (#887).
- Do not rehash large imported files every `make()` (#878).
- Repair parsing of long tidy eval inputs in the DSL (#878).
- Clear up cache confusion when a custom cache exists adjacent to the default cache (#883).
- Accept targets as symbols in `r_drake_build()`.
- Log progress during `r_make()` (#889).
- Repair `expose_imports()`: do not do the `environment<-` trick unless the object is a non-primitive function.
- Use different static analyses of `assign()` vs `delayedAssign()`.
- Fix a superfluous code analysis warning incurred by multiple `file_in()` files and other strings (#896).
- Make `ignore()` work inside `loadd()`, `readd()`, `file_in()`, `file_out()`, and `knitr_in()`.

## New features

- Add experimental support for URLs in `file_in()` and `file_out()`. `drake` now treats `file_in()`/`file_out()` files as URLS if they begin with "http://", "https://", or "ftp://". The fingerprint is a concatenation of the ETag and last-modified timestamp. If neither can be found or if there is no internet connection, `drake` throws an error.
- Implement new memory management strategies `"unload"` and `"none"`, which do not attempt to load a target's dependencies from memory (#897).
- Allow users to give each target its own memory strategy (#897).
- Add `drake_slice()` to help split data across multiple targets. Related: #77, #685, #833.
- Introduce a new `drake_cache()` function, which is now recommended instead of `get_cache()` (#883).
- Introduce a new `r_deps_target()` function.
- Add RStudio addins for `r_make()`, `r_vis_drake_graph()`, and `r_outdated()` (#892).

## Enhancements

- Deprecate `get_cache()` in favor of `drake_cache()`.
- Show the path to the cache in the `clean()` menu prompt.
- Stop removing the console log file on each call to `drake_config()`.
- Log the node name (short host name) and process ID in the console log file.
- Log the name of the calling function in the console log file, e.g. "begin make()" and "end make()". Applies to all functions that accept a `config` argument.
- Memory management: set `use_cache` to `FALSE` in `storr` function calls for saving and loading targets. Also, at the end of `make()`, call `flush_cache()` (and then `gc()` if garbage collection is enabled).
- Mention `callr::r()` within commands as a safe alternative to `lock_envir = FALSE` in the self-invalidation section of the `make()` help file.
- Use file size to help decide when to rehash `file_in()`/`file_out()`/`knitr_in()` files. We now rehash files if the file is less than 100 KB or the time stamp changed or the **file size** changed.

# Version 7.3.0

## Bug fixes

- Accommodate `rlang`'s new interpolation operator `{{`, which was causing `make()` to fail when `drake_plan()` commands are enclosed in curly braces (#864).
- Move "`config$lock_envir <- FALSE`" from `loop_build()` to  `backend_loop()`. This makes sure `config$envir` is correctly locked in `make(parallelism = "clustermq")`.
- Convert factors to characters in the optional `.data` argument of `map()` and `cross()` in the DSL.
- In the DSL of `drake_plan()`, repair `cross(.data = !!args)`, where `args` is an optional data frame of grouping variables.
- Handle trailing slashes in `file_in()`/`file_out()` directories for Windows (#855).
- Make `.id_chr` work with `combine()` in the DSL (#867).
- Do not try `make_spinner()` unless the version of `cli` is at least 1.1.0.

## New features

- Add functions `text_drake_graph()` (and `r_text_drake_graph()` and `render_text_drake_graph()`). Uses text art to print a dependency graph to the terminal window. Handy for when users SSH into remote machines without X Window support.
- Add a new `max_expand` argument to `drake_plan()`, an optional upper bound on the lengths of grouping variables for `map()` and `cross()` in the DSL. Comes in handy when you have a massive number of targets and you want to test on a miniature version of your workflow before you scale up to production.

## Enhancements

- Delay the initialization of `clustermq` workers for as long as possible. Before launching them, build/check targets locally until we reach an outdated target with `hpc` equal to `FALSE`. In other words, if no targets actually require `clustermq` workers, no workers get created.
- In `make(parallelism = "future")`, reset the `config$sleep()` backoff interval whenever a new target gets checked.
- Add a "done" message to the console log file when the workflow has completed.
- Replace `CodeDepends` with a base R solution in `code_to_plan()`. Fixes a CRAN note.
- The DSL (transformations in `drake_plan()`) is no longer experimental.
- The `callr` API (`r_make()` and friends) is no longer experimental.
- Deprecate the wildcard/text-based functions for creating plans: `evaluate_plan()`, `expand_plan()`, `map_plan()`, `gather_plan()`, `gather_by()`, `reduce_plan()`, `reduce_by()`. 
- Change some deprecated functions to defunct: `deps()`, `max_useful_jobs()`, and `migrate_drake_project()`.

# Version 7.2.0

## Mildly breaking changes

- In the DSL (e.g. `drake_plan(x = target(..., transform = map(...)))` avoid inserting extra dots in target names when the grouping variables are character vectors (#847). Target names come out much nicer this way, but those name changes will invalidate some targets (i.e. they need to be rebuilt with `make()`).

## Bug fixes

- Use `config$jobs_preprocess` (local jobs) in several places where `drake` was incorrectly using `config$jobs` (meant for targets).
- Allow `loadd(x, deps = TRUE, config = your_config)` to work even if `x` is not cached (#830). Required disabling `tidyselect` functionality when `deps` `TRUE`. There is a new note in the help file about this, and an informative console message prints out on `loadd(deps = TRUE, tidyselect = TRUE)`. The default value of `tidyselect` is now `!deps`.
- Minor: avoid printing messages and warnings twice to the console (#829).
- Ensure compatibility with `testthat` >= 2.0.1.9000.

## New features

- In `drake_plan()` transformations, allow the user to refer to a target's own name using a special `.id_chr` symbol, which is treated like a character string.
- Add a `transparency` argument to `drake_ggraph()` and `render_drake_ggraph()` to disable transparency in the rendered graph. Useful for R installations without transparency support.

## Enhancements

- Use a custom layout to improve node positions and aspect ratios of `vis_drake_graph()` and `drake_ggraph()` displays. Only activated in `vis_drake_graph()` when there are at least 10 nodes distributed in both the vertical and horizontal directions.
- Allow nodes to be dragged both vertically and horizontally in `vis_drake_graph()` and `render_drake_graph()`.
- Prevent dots from showing up in target names when you supply grouping variables to transforms in `drake_plan()` (#847).
- Do not keep `drake` plans (`drake_plan()`) inside `drake_config()` objects. When other bottlenecks are removed, this will reduce the burden on memory (re #800).
- Do not retain the `targets` argument inside `drake_config()` objects. This is to reduce memory consumption.
- Deprecate the `layout` and `direction` arguments of `vis_drake_graph()` and `render_drake_graph()`. Direction is now always left to right and the layout is always Sugiyama.
- Write the cache log file in CSV format (now `drake_cache.csv` by default) to avoid issues with spaces (e.g. entry names with spaces in them, such as "file report.Rmd")`.


# Version 7.1.0

## Bug fixes

- In `drake` 7.0.0, if you run `make()` in interactive mode and respond to the menu prompt with an option other than `1` or `2`, targets will still build. 
- Make sure file outputs show up in `drake_graph()`. The bug came from `append_output_file_nodes()`, a utility function of `drake_graph_info()`.
- Repair `r_make(r_fn = callr::r_bg())` re #799.
- Allow `drake_ggraph()` and `sankey_drake_graph()` to work when the graph has no edges.

## New features

- Add a new `use_drake()` function to write the `make.R` and `_drake.R` files from the "main example". Does not write other supporting scripts.
- With an optional logical `hpc` column in your `drake_plan()`, you can now select which targets to deploy to HPC and which to run locally.
- Add a `list` argument to `build_times()`, just like `loadd()`.
- Add a new RStudio addin: 'loadd target at cursor' which can be bound a keyboard shortcut to load the target identified by the symbol at the cursor position to the global environment.

## Enhancements

- `file_in()` and `file_out()` can now handle entire directories, e.g. `file_in("your_folder_of_input_data_files")` and `file_out("directory_with_a_bunch_of_output_files")`.
- Send less data from `config` to HPC workers.
- Improve `drake_ggraph()`
  - Hide node labels by default and render the arrows behind the nodes.
  - Print an informative error message when the user supplies a `drake` plan to the `config` argument of a function.
  - By default, use gray arrows and a black-and-white background with no gridlines.
- For the `map()` and `cross()` transformations in the DSL, prevent the accidental sorting of targets by name (#786). Needed `merge(sort = FALSE)` in `dsl_left_outer_join()`.
- Simplify verbosity. The `verbose` argument of `make()` now takes values 0, 1, and 2, and maximum verbosity in the console prints targets, retries, failures, and a spinner. The console log file, on the other hand, dumps maximally verbose runtime info regardless of the `verbose` argument.
- In previous versions, functions generated with `f <- Rcpp::cppFunction(...)` did not stay up to date from session to session because the addresses corresponding to anonymous pointers were showing up in `deparse(f)`. Now, `drake` ignores those pointers, and `Rcpp` functions compiled inline appear to stay up to date. This problem was more of an edge case than a bug.
- Prepend time stamps with sub-second times to the lines of the console log file.
- In `drake_plan()`, deprecate the `tidy_evaluation` argument in favor of the new and more concise `tidy_eval`. To preserve back compatibility for now, if you supply a non-`NULL` value to `tidy_evaluation`, it overwrites `tidy_eval`.
- Reduce the object size of `drake_config()` objects by assigning closure of `config$sleep` to `baseenv()`.

# Version 7.0.0

## Breaking changes

- The enhancements that increase cache access speed also invalidate targets in old projects. Workflows built with drake <= 6.2.1 will need to run from scratch again.
- In `drake` plans, the `command` and `trigger` columns are now lists of language objects instead of character vectors. `make()` and friends still work if you have character columns, but the default output of `drake_plan()` has changed to this new format.
- All parallel backends (`parallelism` argument of `make()`) except "clustermq" and "future" are removed. A new "loop" backend covers local serial execution.
- A large amount of deprecated functionality is now defunct, including several functions (`built()`, `find_project()`, `imported()`, and `parallel_stages()`; full list at #564) and the single-quoted file API.
- Set the default value of `lock_envir` to `TRUE` in `make()` and `drake_config()`. So `make()` will automatically quit in error if the act of building a target tries to change upstream dependencies.
- `make()` no longer returns a value. Users will need to call `drake_config()` separately to get the old return value of `make()`.
- Require the `jobs` argument to be of length 1 (`make()` and `drake_config()`). To parallelize the imports and other preprocessing steps, use `jobs_preprocess`, also of length 1.
- Get rid of the "kernels" `storr` namespace. As a result, `drake` is faster, but users will no longer be able to load imported functions using `loadd()` or `readd()`.
- In `target()`, users must now explicitly name all the arguments except `command`, e.g. `target(f(x), trigger = trigger(condition = TRUE))` instead of `target(f(x), trigger(condition = TRUE))`.
- Fail right away in `bind_plans()` when the result has duplicated target names. This makes `drake`'s API more predictable and helps users catch malformed workflows earlier.
- `loadd()` only loads targets listed in the plan. It no longer loads imports or file hashes.
- The return values of `progress()`, `deps_code()`, `deps_target()`, and `predict_workers()` are now data frames.
- Change the default value of `hover` to `FALSE` in visualization functions. Improves speed.

## Bug fixes

- Allow `bind_plans()` to work with lists of plans (`bind_plans(list(plan1, plan2))` was returning `NULL` in `drake` 6.2.0 and 6.2.1).
- Ensure that `get_cache(path = "non/default/path", search = FALSE)` looks for the cache in `"non/default/path"` instead of `getwd()`.
- Remove strict dependencies on package `tibble`.
- Pass the correct data structure to `ensure_loaded()` in `meta.R` and `triggers.R` when ensuring the dependencies of the `condition` and `change` triggers are loaded.
- Require a `config` argument to `drake_build()` and `loadd(deps = TRUE)`.

## New features

- Introduce a new experimental domain-specific language for generating large plans (#233). Details in the "Plans" chapter of the manual.
- Implement a `lock_envir` argument to safeguard reproducibility. More discussion: #619, #620.
- The new `from_plan()` function allows the users to reference custom plan columns from within commands. Changes to values in these columns columns do not invalidate targets.
- Add a menu prompt (#762) to safeguard against `make()` pitfalls in interactive mode (#761). Appears once per session. Disable with `options(drake_make_menu = FALSE)`.
- Add new API functions `r_make()`, `r_outdated()`, etc. to run `drake` functions more reproducibly in a clean session. See the help file of `r_make()` for details.
- `progress()` gains a `progress` argument for filtering results. For example, `progress(progress = "failed")` will report targets that failed.


## Enhancements

- **Large speed boost**: move away from `storr`'s key mangling in favor of `drake`'s own encoding of file paths and namespaced functions for `storr` keys.
- Exclude symbols `.`, `..`, and `.gitignore` from being target names (consequence of the above).
- Use only one hash algorithm per `drake` cache, which the user can set with the `hash_algorithm` argument of `new_cache()`, `storr::storr_rds()`, and various other cache functions. Thus, the concepts of a "short hash algorithm" and "long hash algorithm" are deprecated, and the functions `long_hash()`, `short_hash()`, `default_long_hash_algo()`, `default_short_hash_algo()`, and `available_hash_algos()` are deprecated. Caches are still back-compatible with `drake` > 5.4.0 and <= 6.2.1.
- Allow the `magrittr` dot symbol to appear in some commands sometimes.
- Deprecate the `fetch_cache` argument in all functions.
- Remove packages `DBI` and `RSQLite` from "Suggests".
- Define a special `config$eval <- new.env(parent = config$envir)` for storing built targets and evaluating commands in the plan. Now, `make()` no longer modifies the user's environment. This move is a long-overdue step toward purity.
- Remove dependency on the `codetools` package.
- Deprecate and remove the `session` argument of `make()` and `drake_config()`. Details: in #623.
- Deprecate the `graph` and `layout` arguments to `make()` and `drake_config()`. The change simplifies the internals, and memoization allows us to do this.
- Warn the user if running `make()` in a subdirectory of the `drake` project root (determined by the location of the `.drake` folder in relation to the working directory).
- In the code analysis, explicitly prohibit targets from being dependencies of imported functions.
- Increase options for the `verbose` argument, including the option to print execution and total build times.
- Separate the building of targets from the processing of imports. Imports are processed with rudimentary staged parallelism (`mclapply()` or `parLapply()`, depending on the operating system).
- Ignore the imports when it comes to build times. Functions `build_times()`, `predict_runtime()`, etc. focus on only the targets.
- Deprecate many API functions, including `plan_analyses()`, `plan_summaries()`, `analysis_wildcard()`, `cache_namespaces()`, `cache_path()`, `check_plan()`, `dataset_wildcard()`, `drake_meta()`, `drake_palette()`, `drake_tip()`, `recover_cache()`, `cleaned_namespaces()`, `target_namespaces()`, `read_drake_config()`, `read_drake_graph()`, and `read_drake_plan()`.
- Deprecate `target()` as a user-side function. From now on, it should only be called from within `drake_plan()`.
- `drake_envir()` now throws an error, not a warning, if called in the incorrect context. Should be called only inside commands in the user's `drake` plan.
- Replace `*expr*()` `rlang` functions with their `*quo*()` counterparts. We still keep `rlang::expr()` in the few places where we know the expressions need to be evaluated in `config$eval`.
- The `prework` argument to `make()` and `drake_config()` can now be an expression (language object) or list of expressions. Character vectors are still acceptable.
- At the end of `make()`, print messages about triggers etc. only if `verbose >= 2L`.
- Deprecate and rename  `in_progress()` to `running()`.
- Deprecate and rename  `knitr_deps()` to `deps_knitr()`.
- Deprecate and rename  `dependency_profile()` to `deps_profile()`.
- Deprecate and rename  `predict_load_balancing()` to `predict_workers()`.
- Deprecate `this_cache()` and defer to `get_cache()` and `storr::storr_rds()` for simplicity.
- Change the default value of `hover` to `FALSE` in visualization functions. Improves speed. Also a breaking change.
- Deprecate `drake_cache_log_file()`. We recommend using `make()` with the `cache_log_file` argument to create the cache log. This way ensures that the log is always up to date with `make()` results.


# Version 6.2.1

Version 6.2.1 is a hotfix to address the failing automated CRAN checks for 6.2.0. Chiefly, in CRAN's Debian R-devel (2018-12-10) check platform, errors of the form "length > 1 in coercion to logical" occurred when either argument to `&&` or `||` was not of length 1 (e.g. `nzchar(letters) && length(letters)`). In addition to fixing these errors, version 6.2.1 also removes a problematic link from the vignette.


# Version 6.2.0

## New features

- Add a `sep` argument to `gather_by()`, `reduce_by()`, `reduce_plan()`, `evaluate_plan()`, `expand_plan()`, `plan_analyses()`, and `plan_summaries()`. Allows the user to set the delimiter for generating new target names.
- Expose a `hasty_build` argument to `make()` and `drake_config()`. Here, the user can set the function that builds targets in "hasty mode" (`make(parallelism = "hasty")`).
- Add a new `drake_envir()` function that returns the environment where `drake` builds targets. Can only be accessed from inside the commands in the workflow plan data frame. The primary use case is to allow users to remove individual targets from memory at predetermined build steps.

## Bug fixes

- Ensure compatibility with `tibble` 2.0.0.
- Stop returning `0s` from `predict_runtime(targets_only = TRUE)` when some targets are outdated and others are not.
- Remove `sort(NULL)` warnings from `create_drake_layout()`. (Affects R-3.3.x.)

## Enhancements

- Remove strict dependencies on packages `evaluate`, `formatR`, `fs`, `future`, `parallel`, `R.utils`, `stats`, and `stringi`.
- **Large speed boost**: reduce repeated calls to `parse()` in `code_dependencies()`.
- **Large speed boost**: change the default value of `memory_strategy` (previously `pruning_strategy`) to `"speed"` (previously `"lookahead"`).
- Compute a special data structure in `drake_config()` (`config$layout`) just to store the code analysis results. This is an intermediate structure between the workflow plan data frame and the graph. It will help clean up the internals in future development.
- Improve memoized preprocessing: deparse all the functions in the environment so the memoization does not react so spurious changes in R internals. Related: #345.
- Use the `label` argument to `future()` inside `make(parallelism = "future")`. That way , job names are target names by default if `job.name` is used correctly in the `batchtools` template file.
- Remove strict dependencies on packages `dplyr`, `evaluate`, `fs`, `future`, `magrittr`, `parallel`, `R.utils`, `stats`, `stringi`, `tidyselect`, and `withr`.
- Remove package `rprojroot` from "Suggests".
- Deprecate the `force` argument in all functions except `make()` and `drake_config()`.
- Change the name of `prune_envir()` to `manage_memory()`.
- Deprecate and rename the `pruning_strategy` argument to `memory_strategy` (`make()` and `drake_config()`).
- Print warnings and messages to the `console_log_file` in real time (#588).
- Use HTML line breaks in `vis_drake_graph()` hover text to display commands in the `drake` plan more elegantly.
- Speed up `predict_load_balancing()` and remove its reliance on internals that will go away in 2019 via #561.
- Remove support for the `worker` column of `config$plan` in `predict_runtime()` and `predict_load_balancing()`. This functionality will go away in 2019 via #561.
- Change the names of the return value of `predict_load_balancing()` to `time` and `workers`.
- Bring the documentation of `predict_runtime()` and `predict_load_balancing()` up to date.
- Deprecate `drake_session()` and rename to `drake_get_session_info()`.
- Deprecate the `timeout` argument in the API of `make()` and `drake_config()`. A value of `timeout` can be still passed to these functions without error, but only the `elapsed` and `cpu` arguments impose actual timeouts now.

# Version 6.1.0

## New features

- **Add a new `map_plan()` function to easily create a workflow plan data frame to execute a function call over a grid of arguments.**
- Add a new `plan_to_code()` function to turn `drake` plans into generic R scripts. New users can use this function to better understand the relationship between plans and code, and unsatisfied customers can use it to disentangle their projects from `drake` altogether. Similarly, `plan_to_notebook()` generates an R notebook from a `drake` plan.
- Add a new `drake_debug()` function to run a target's command in debug mode. Analogous to `drake_build()`.
- Add a `mode` argument to `trigger()` to control how the `condition` trigger factors into the decision to build or skip a target. See the `?trigger` for details.
- Add a new `sleep` argument to `make()` and `drake_config()` to help the main process consume fewer resources during parallel processing.
- Enable the `caching` argument for the `"clustermq"` and `"clustermq_staged"` parallel backends. Now, `make(parallelism = "clustermq", caching = "main")` will do all the caching with the main process, and `make(parallelism = "clustermq", caching = "worker")` will do all the caching with the workers. The same is true for `parallelism = "clustermq_staged"`.
- Add a new `append` argument to `gather_plan()`, `gather_by()`, `reduce_plan()`, and `reduce_by()`. The `append` argument control whether the output includes the original `plan` in addition to the newly generated rows.
- Add new functions `load_main_example()`, `clean_main_example()`, and `clean_mtcars_example()`.
- Add a `filter` argument to `gather_by()` and `reduce_by()` in order to restrict what we gather even when `append` is `TRUE`.
- Add a hasty mode: `make(parallelism = "hasty")` skips all of `drake`'s expensive caching and checking. All targets run every single time and you are responsible for saving results to custom output files, but almost all the by-target overhead is gone.

## Bug fixes

- Ensure commands in the plan are re-analyzed for dependencies when new imports are added (#548). Was a bug in version 6.0.0 only.
- Call `path.expand()` on the `file` argument to `render_drake_graph()` and `render_sankey_drake_graph()`. That way, tildes in file paths no longer interfere with the rendering of static image files.
- Skip tests and examples if the required "Suggests" packages are not installed.
- Stop checking for non-standard columns. Previously, warnings about non-standard columns were incorrectly triggered by `evaluate_plan(trace = TRUE)` followed by `expand_plan()`, `gather_plan()`, `reduce_plan()`, `gather_by()`, or `reduce_by()`. The more relaxed behavior also gives users more options about how to construct and maintain their workflow plan data frames.
- Use checksums in `"future"` parallelism to make sure files travel over network file systems before proceeding to downstream targets.
- Refactor and clean up checksum code.
- Skip more tests and checks if the optional `visNetwork` package is not installed.

## Enhancements

- Stop earlier in `make_targets()` if all the targets are already up to date.
- Improve the documentation of the `seed` argument in `make()` and `drake_config()`.
- Set the default `caching` argument of `make()` and `drake_config()` to `"main"` rather than `"worker"`. The default option should be the lower-overhead option for small workflows. Users have the option to make a different set of tradeoffs for larger workflows.
- Allow the `condition` trigger to evaluate to non-logical values as long as those values can be coerced to logicals.
- Require that the `condition` trigger evaluate to a vector of length 1.
- Keep non-standard columns in `drake_plan_source()`.
- `make(verbose = 4)` now prints to the console when a target is stored.
- `gather_by()` and `reduce_by()` now gather/reduce everything if no columns are specified.
- Change the default parallelization of the imports. Previously, `make(jobs = 4)` was equivalent to `make(jobs = c(imports = 4, targets = 4))`. Now, `make(jobs = 4)` is equivalent to `make(jobs = c(imports = 1, targets = 4))`. See issue #553 for details.
- Add a console message for building the priority queue when `verbose` is at least 2.
- Condense `load_mtcars_example()`.
- Deprecate the `hook` argument of `make()` and `drake_config()`.
- In `gather_by()` and `reduce_by()`, do not exclude targets with all `NA` gathering variables.

# Version 6.0.0

## Breaking changes

- Avoid serialization in `digest()` wherever possible. This puts old `drake` projects out of date, but it improves speed.
- Require R version >= 3.3.0 rather than >= 3.2.0. Tests and checks still run fine on 3.3.0, but the required version of the `stringi` package no longer compiles on 3.2.0.
- Be more discerning in detecting dependencies. In `code_dependencies()`, restrict the possible global variables to the ones mentioned in the new `globals` argument (turned off when `NULL`. In practical workflows, global dependencies are restricted to items in `envir` and proper targets in the plan. In `deps_code()`, the `globals` slot of the output list is now a list of *candidate* globals, not necessarily actual globals (some may not be targets or variables in `envir`).

## Bug fixes

- In the call to `unlink()` in `clean()`, set `recursive` and `force` to `FALSE`. This should prevent the accidental deletion of whole directories.
- Previously, `clean()` deleted input-only files if no targets from the plan were cached. A patch and a unit test are included in this release.
- `loadd(not_a_target)` no longer loads every target in the cache.
- Exclude each target from its own dependency metadata in the "deps" `igraph` vertex attribute (fixes #503).
- Detect inline code dependencies in `knitr_in()` file code chunks.
- Remove more calls to `sort(NULL)` that caused warnings in R 3.3.3.
- Fix a bug on R 3.3.3 where `analyze_loadd()` was sometimes quitting with "Error: attempt to set an attribute on NULL".
- Do not call `digest::digest(file = TRUE)` on directories. Instead, set hashes of directories to `NA`. Users should still not directories as file dependencies.
- If files are declared as dependencies of custom triggers ("condition" and "change") include them in `vis_drake_graph()`. Previously, these files were missing from the visualization, but actual workflows worked just fine.
- Work around mysterious `codetools` failures in R 3.3 (add a `tryCatch()` statement in `find_globals()`).

## New features

- Add a proper `clustermq`-based parallel backend: `make(parallelism = "clustermq")`.
- `evaluate_plan(trace = TRUE)` now adds a `*_from` column to show the origins of the evaluated targets. Try `evaluate_plan(drake_plan(x = rnorm(n__), y = rexp(n__)), wildcard = "n__", values = 1:2, trace = TRUE)`.
- Add functions `gather_by()` and `reduce_by()`, which gather on custom columns in the plan (or columns generated by `evaluate_plan(trace = TRUE)`) and append the new targets to the previous plan.
- Expose the `template` argument of `clustermq` functions (e.g. `Q()` and `workers()`) as an argument of `make()` and `drake_config()`.
- Add a new `code_to_plan()` function to turn R scripts and R Markdown reports into workflow plan data frames.
- Add a new `drake_plan_source()` function, which generates lines of code for a `drake_plan()` call. This `drake_plan()` call produces the plan passed to `drake_plan_source()`. The main purpose is visual inspection (we even have syntax highlighting via `prettycode`) but users may also save the output to a script file for the sake of reproducibility or simple reference.
- Deprecate `deps_targets()` in favor of a new `deps_target()` function (singular) that behaves more like `deps_code()`.

## Enhancements

- Smooth the edges in `vis_drake_graph()` and `render_drake_graph()`.
- Make hover text slightly more readable in in `vis_drake_graph()` and `render_drake_graph()`.
- Align hover text properly in `vis_drake_graph()` using the "title" node column.
- Optionally collapse nodes into clusters with `vis_drake_graph(collapse = TRUE)`.
- Improve `dependency_profile()` show major trigger hashes side-by-side
to tell the user if the command, a dependency, an input file, or an output file changed since the last `make()`.
- Choose more appropriate places to check that the `txtq` package is installed.
- Improve the help files of `loadd()` and `readd()`, giving specific usage guidance in prose.
- Memoize all the steps of `build_drake_graph()` and print to the console the ones that execute.
- Skip some tests if `txtq` is not installed.

# Version 5.4.0

- Overhaul the interface for triggers and add new trigger types ("condition" and "change").
- Offload `drake`'s code examples to the `drake-examples` GitHub repository and make make `drake_example()` and `drake_examples()` download examples from there.
- Optionally show output files in graph visualizations. See the `show_output_files` argument to `vis_drake_graph()` and friends.
- Repair output file checksum operations for distributed backends like `"clustermq_staged"` and `"future_lapply"`.
- Internally refactor the `igraph` attributes of the dependency graph to allow for smarter dependency/memory management during `make()`.
- Enable `vis_drake_graph()` and `sankey_drake_graph()` to save static image files via `webshot`.
- Deprecate `static_drake_graph()` and `render_static_drake_graph()` in favor of `drake_ggraph()` and `render_drake_ggraph()`.
- Add a `columns` argument to `evaluate_plan()` so users can evaluate wildcards in columns other than the `command` column of `plan`.
- Name the arguments of `target()` so users do not have to (explicitly).
- Lay the groundwork for a special pretty print method for workflow plan data frames.

# Version 5.3.0

- Allow multiple output files per command.
- Add Sankey diagram visuals: `sankey_drake_graph()` and `render_sankey_drake_graph()`.
- Add `static_drake_graph()` and `render_static_drake_graph()` for `ggplot2`/`ggraph` static graph visualizations.
- Add `group` and `clusters` arguments to `vis_drake_graph()`, `static_drake_graph()`, and `drake_graph_info()` to optionally condense nodes into clusters.
- Implement a `trace` argument to `evaluate_plan()` to optionally add indicator columns to show which targets got expanded/evaluated with which wildcard values.
- Rename the `always_rename` argument to `rename` in `evaluate_plan()`.
- Add a `rename` argument to `expand_plan()`.
- Implement `make(parallelism = "clustermq_staged")`, a `clustermq`-based staged parallelism backend (see #452).
- Implement `make(parallelism = "future_lapply_staged")`, a `future`-based staged parallelism backend (see #450).
- Depend on `codetools` rather than `CodeDepends` for finding global variables.
- Detect `loadd()` and `readd()` dependencies in `knitr` reports referenced with `knitr_in()` inside imported functions. Previously, this feature was only available in explicit `knitr_in()` calls in commands.
- Skip more tests on CRAN. White-list tests instead of blacklisting them in order to try to keep check time under the official 10-minute cap.
- Disallow wildcard names to grep-match other wildcard names or any replacement values. This will prevent careless mistakes and confusion when generating `drake_plan()`s.
- Prevent persistent workers from hanging when a target fails.
- Move the example template files to `inst/hpc_template_files`.
- Deprecate `drake_batchtools_tmpl_file()` in favor of `drake_hpc_template_file()` and `drake_hpc_template_files()`.
- Add a `garbage_collection` argument to `make()`. If `TRUE`, `gc()` is called after every new build of a target.
- Remove redundant calls to `sanitize_plan()` in `make()`.
- Change `tracked()` to accept only a `drake_config()` object as an argument. Yes, it is technically a breaking change, but it is only a small break, and it is the correct API choice.
- Move visualization and hpc package dependencies to "Suggests:" rather than "Imports:" in the `DESCRIPTION` file.
- Allow processing of codeless `knitr` reports without warnings.

# Version 5.2.1

- Skip several long-running and low-priority tests on CRAN.

# Version 5.2.0

- Sequester staged parallelism in backends "mclapply_staged" and "parLapply_staged". For the other `lapply`-like backends, `drake` uses persistent workers and a main process. In the case of `"future_lapply"` parallelism, the main process is a separate background process called by `Rscript`.
- Remove the appearance of staged parallelism from single-job `make()`'s.
(Previously, there were "check" messages and a call to `staged_parallelism()`.)
- Remove some remnants of staged parallelism internals.
- Allow different parallel backends for imports vs targets. For example, `make(parallelism = c(imports = "mclapply_staged", targets = "mclapply")`.
- Fix a bug in environment pruning. Previously, dependencies of downstream targets were being dropped from memory in `make(jobs = 1)`. Now, they are kept in memory until no downstream target needs them (for `make(jobs = 1)`).
- Improve `predict_runtime()`. It is a more sensible way to go about predicting runtimes with multiple jobs. Likely to be more accurate.
- Calls to `make()` no longer leave targets in the user's environment.
- Attempt to fix a Solaris CRAN check error. A test was previously failing on CRAN's Solaris machine (R 3.5.0). In the test, one of the threads deliberately quits in error, and the R/Solaris installation did not handle this properly. The test should work now because it no longer uses any parallelism.
- Deprecate the `imports_only` argument to `make()` and `drake_config()` in favor of `skip_targets`.
- Deprecate `migrate_drake_project()`.
- Deprecate `max_useful_jobs()`.
- For non-distributed parallel backends, stop waiting for all the imports to finish before the targets begin.
- Add an `upstream_only` argument to `failed()` so users can list failed targets that do not have any failed dependencies. Naturally accompanies `make(keep_going = TRUE)`.
- Add an RStudio R Markdown template.
- Remove `plyr` as a dependency.
- Handle duplicated targets better in `drake_plan()` and `bind_plans()`.
- Add a true function `target()` to help create drake plans with custom columns.
- In `drake_gc()`, clean out disruptive files in `storr`s with mangled keys (re: #198).
- Move all the vignettes to the up and coming user manual.
- Rename the "basic example" to the "mtcars example".
- Deprecate `load_basic_example()` in favor of `load_mtcars_example()`.
- Refocus the `README.md` file on the main example rather than the mtcars example.
- Use a `README.Rmd` file to generate `README.md`.
- Add function `deps_targets()`.
- Deprecate function `deps()` in favor of `deps_code()`
- Add a `pruning_strategy` argument to `make()` and `drake_config()` so the user can decide how `drake` keeps non-import dependencies in memory when it builds a target.
- Add optional custom (experimental) "workers" and "priorities" columns to the `drake` plans to help users customize scheduling.
- Add a `makefile_path` argument to `make()` and `drake_config()` to avoid potential conflicts between user-side custom `Makefile`s and the one written by `make(parallelism = "Makefile")`.
- Document batch mode for long workflows in the HPC guide.
- Add a `console` argument to `make()` and `drake_config()` so users can redirect console output to a file.
- Make it easier for the user to find out where a target in the cache came from: `show_source()`, `readd(show_source = TRUE)`, `loadd(show_source = TRUE)`.

# Version 5.1.2

- In R 3.5.0, the `!!` operator from tidyeval and `rlang` is parsed differently than in R <= 3.4.4. This change broke one of the tests in `tests/testthat/tidy-eval.R` The main purpose of `drake`'s 5.1.2 release is to fix the broken test.
- Fix an elusive `R CMD check` error from building the pdf manual with LaTeX.
- In `drake_plan()`, allow users to customize target-level columns using `target()` inside the commands.
- Add a new `bind_plans()` function to concatenate the rows of drake plans and then sanitize the aggregate plan.
- Add an optional `session` argument to tell `make()` to build targets in a separate, isolated main R session. For example, `make(session = callr::r_vanilla)`.

# Version 5.1.0

- Add a `reduce_plan()` function to do pairwise reductions on collections of targets.
- Forcibly exclude the dot (`.`) from being a dependency of any target or import. This enforces more consistent behavior in the face of the current static code analysis functionality, which sometimes detects `.` and sometimes does not.
- Use `ignore()` to optionally ignore pieces of workflow plan commands and/or imported functions. Use `ignore(some_code)` to
    1. Force `drake` to not track dependencies in `some_code`, and
    2. Ignore any changes in `some_code` when it comes to deciding which target are out of date.
- Force `drake` to only look for imports in environments inheriting from `envir` in `make()` (plus explicitly namespaced functions).
- Force `loadd()` to ignore foreign imports (imports not explicitly found in `envir` when `make()` last imported them).
- Reduce default verbosity. Only targets are printed out by default. Verbosity levels are integers ranging from 0 through 4.
- Change `loadd()` so that only targets (not imports) are loaded if the `...` and `list` arguments are empty.
- Add check to drake_plan() to check for duplicate targets
- Add a `.gitignore` file containing `"*"` to the default `.drake/` cache folder every time `new_cache()` is called. This means the cache will not be automatically committed to git. Users need to remove `.gitignore` file to allow unforced commits, and then subsequent `make()`s on the same cache will respect the user's wishes and not add another `.gitignore`. this only works for the default cache. Not supported for manual `storr`s.
- Add a new experimental `"future"` backend with a manual scheduler.
- Implement `dplyr`-style `tidyselect` functionality in `loadd()`, `clean()`, and `build_times()`. For `build_times()`, there is an API change: for `tidyselect` to work, we needed to insert a new `...` argument as the first argument of `build_times()`.
- Deprecate the single-quoting API for files. Users should now use formal API functions in their commands:
    - `file_in()` for file inputs to commands or imported functions (for imported functions, the input file needs to be an imported file, not a target).
    - `file_out()` for output file targets (ignored if used in imported functions).
    - `knitr_in()` for `knitr`/`rmarkdown` reports. This tells `drake` to look inside the source file for target dependencies in code chunks (explicitly referenced with `loadd()` and `readd()`). Treated as a `file_in()` if used in imported functions.
- Change `drake_plan()` so that it automatically fills in any target names that the user does not supply. Also, any `file_out()`s become the target names automatically (double-quoted internally).
- Make `read_drake_plan()` (rather than an empty `drake_plan()`) the default `plan` argument in all functions that accept a `plan`.
- Add support for active bindings: `loadd(..., lazy = "bind")`. That way, when you have a target loaded in one R session and hit `make()` in another R session, the target in your first session will automatically update.
- Use tibbles for workflow plan data frames and the output of `dataframes_graph()`.
- Return warnings, errors, and other context of each build, all wrapped up with the usual metadata. `diagnose()` will take on the role of returning this metadata.
- Deprecate the `read_drake_meta()` function in favor of `diagnose()`.
- Add a new `expose_imports()` function to optionally force `drake` detect deeply nested functions inside specific packages.
- Move the "quickstart.Rmd" vignette to "example-basic.Rmd". The so-called "quickstart" didn't end up being very quick, and it was all about the basic example anyway.
- Move `drake_build()` to be an exclusively user-side function.
- Add a `replace` argument to `loadd()` so that objects already in the user's environment need not be replaced.
- When the graph cyclic, print out all the cycles.
- Prune self-referential loops (and duplicate edges) from the workflow graph. That way, recursive functions are allowed.
- Add a `seed` argument to `make()`, `drake_config()`, and `load_basic_example()`. Also hard-code a default seed of `0`. That way, the pseudo-randomness in projects should be reproducible
across R sessions.
- Cache the pseudo-random seed at the time the project is created and use that seed to build targets until the cache is destroyed.
- Add a new `drake_read_seed()` function to read the seed from the cache. Its examples illustrate what `drake` is doing to try to ensure reproducible random numbers.
- Evaluate the quasiquotation operator `!!` for the `...` argument to `drake_plan()`. Suppress this behavior using `tidy_evaluation = FALSE` or by passing in commands passed through the `list` argument.
- Preprocess workflow plan commands with `rlang::expr()` before evaluating them. That means you can use the quasiquotation operator `!!` in your commands, and `make()` will evaluate them according to the tidy evaluation paradigm.
- Restructure `drake_example("basic")`, `drake_example("gsp")`, and `drake_example("packages")` to demonstrate how to set up the files for serious `drake` projects. More guidance was needed in light of #193.
- Improve the examples of `drake_plan()` in the help file (`?drake_plan`).

# Version 5.0.0

- Transfer `drake` to rOpenSci GitHub URL.
- Several functions now require an explicit `config` argument, which you can get from
`drake_config()` or `make()`. Examples:
    - outdated()
    - missed()
    - rate_limiting_times()
    - predict_runtime()
    - vis_drake_graph()
    - dataframes_graph()
- Always process all the imports before building any targets. This is part of the solution to #168: if imports and targets are processed together, the full power of parallelism is taken away from the targets. Also, the way parallelism happens is now consistent for all parallel backends.
- Major speed improvement: dispense with internal inventories and rely on `cache$exists()` instead.
- Let the user define a trigger for each target to customize when `make()` decides to build targets.
- Document triggers and other debugging/testing tools in the new "debug" vignette.
- Restructure the internals of the `storr` cache in a way that is not back-compatible with projects from versions 4.4.0 and earlier. The main change is to make more intelligent use of `storr` namespaces, improving efficiency (both time and storage) and opening up possibilities for new features. If you attempt to run drake >= 5.0.0 on a project from drake <= 4.0.0, drake will stop you before any damage to the cache is done, and you will be instructed how to migrate your project to the new drake.
- Use `formatR::tidy_source()` instead of `parse()` in `tidy_command()` (originally `tidy()` in `R/dependencies.R`). Previously, `drake` was having problems with an edge case: as a command, the literal string `"A"` was interpreted as the symbol `A` after tidying. With `tidy_source()`, literal quoted strings stay literal quoted strings in commands. This may put some targets out of date in old projects, yet another loss of back compatibility in version 5.0.0.
- Speed up clean() by refactoring the cache inventory and using light parallelism.
- Implement `rescue_cache()`, exposed to the user and used in `clean()`. This function removes dangling orphaned files in the cache so that a broken cache can be cleaned and used in the usual ways once more.
- Change the default `cpu` and `elapsed` arguments of `make()` to `NULL`. This solves an elusive bug in how drake imposes timeouts.
- Allow users to set target-level timeouts (overall, cpu, and elapsed) with columns in the workflow plan data frame.
- Document timeouts and retries in the new "debug" vignette.
- Add a new `graph` argument to functions `make()`, `outdated()`, and `missed()`.
- Export a new `prune_graph()` function for igraph objects.
- Delete long-deprecated functions `prune()` and `status()`.
- Deprecate and rename functions:
    - `analyses()` => `plan_analyses()`
    - `as_file()` => `as_drake_filename()`
    - `backend()` => `future::plan()`
    - `build_graph()` => `build_drake_graph()`
    - `check()` => `check_plan()`
    - `config()` => `drake_config()`
    - `evaluate()` => `evaluate_plan()`
    - `example_drake()` => `drake_example()`
    - `examples_drake()` => `drake_examples()`
    - `expand()` => `expand_plan()`
    - `gather()` => `gather_plan()`
    - `plan()`, `workflow()`, `workplan()` => `drake_plan()`
    - `plot_graph()` => `vis_drake_graph()`
    - `read_config()` => `read_drake_config()`
    - `read_graph()` => `read_drake_graph()`
    - `read_plan()` => `read_drake_plan()`
    - `render_graph()` => `render_drake_graph()`
    - `session()` => `drake_session()`
    - `summaries()` => `plan_summaries()`
- Disallow `output` and `code` as names in the workflow plan data frame. Use `target` and `command` instead. This naming switch has been formally deprecated for several months prior.
- Deprecate the ..analysis.. and ..dataset.. wildcards in favor of analysis__ and dataset__, respectively. The new wildcards are stylistically better an pass linting checks.
- Add new functions `drake_quotes()`, `drake_unquote()`, and `drake_strings()` to remove the silly dependence on the `eply` package.
- Add a `skip_safety_checks` flag to `make()` and `drake_config()`. Increases speed.
- In `sanitize_plan()`, remove rows with blank targets "".
- Add a `purge` argument to `clean()` to optionally remove all target-level information.
- Add a `namespace` argument to `cached()` so users can inspect individual `storr` namespaces.
- Change `verbose` to numeric: 0 = print nothing, 1 = print progress on imports only, 2 = print everything.
- Add a new `next_stage()` function to report the targets to be made in the next parallelizable stage.
- Add a new `session_info` argument to `make()`. Apparently, `sessionInfo()` is a bottleneck for small `make()`s, so there is now an option to suppress it. This is mostly for the sake of speeding up unit tests.
- Add a new `log_progress` argument to `make()` to suppress progress logging. This increases storage efficiency and speeds some projects up a tiny bit.
- Add an optional `namespace` argument to `loadd()` and `readd()`. You can now load and read from non-default `storr` namespaces.
- Add `drake_cache_log()`, `drake_cache_log_file()`, and `make(..., cache_log_file = TRUE)` as options to track changes to targets/imports in the drake cache.
- Detect knitr code chunk dependencies in response to commands with `rmarkdown::render()`, not just `knit()`.
- Add a new general best practices vignette to clear up misconceptions about how to use `drake` properly.

# Version 4.4.0

- Extend `plot_graph()` to display subcomponents. Check out arguments `from`, `mode`, `order`, and `subset`. The graph visualization vignette has demonstrations.
- Add `"future_lapply"` parallelism: parallel backends supported by the `future` and `future.batchtools` packages. See `?backend` for examples and the parallelism vignette for an introductory tutorial. More advanced instruction can be found in the `future` and `future.batchtools` packages themselves.
- Cache diagnostic information of targets that fail and retrieve diagnostic info with `diagnose()`.
- Add an optional `hook` argument to `make()` to wrap around `build()`. That way, users can more easily control the side effects of distributed jobs. For example, to redirect error messages to a file in `make(..., parallelism = "Makefile", jobs = 2, hook = my_hook)`, `my_hook` should be something like `function(code){withr::with_message_sink("messages.txt", code)}`.
- Remove console logging for "parLapply" parallelism. `drake` was previously using the `outfile` argument for PSOCK clusters to generate output that could not be caught by `capture.output()`. It was a hack that should have been removed before.
- Remove console logging for "parLapply" parallelism. `drake` was previously using the `outfile` argument for PSOCK clusters to generate output that could not be caught by `capture.output()`. It was a hack that should have been removed before.
- If 'verbose' is 'TRUE' and all targets are already up to date (nothing to build), then `make()` and `outdated()` print "All targets are already up to date" to the console.
- Add new examples in 'inst/examples', most of them demonstrating how to use the `"future_lapply"` backends.
- New support for timeouts and retries when it comes to building targets.
- Failed targets are now recorded during the build process. You can see them in `plot_graph()` and `progress()`. Also see the new `failed()` function, which is similar to `in_progress()`.
- Speed up the overhead of `parLapply` parallelism. The downside to this fix is that `drake` has to be properly installed. It should not be loaded with `devtools::load_all()`. The speedup comes from lightening the first `clusterExport()` call in `run_parLapply()`. Previously, we exported every single individual `drake` function to all the workers, which created a bottleneck. Now, we just load `drake` itself in each of the workers, which works because `build()` and `do_prework()` are exported.
- Change default value of `overwrite` to `FALSE` in `load_basic_example()`.
- Warn when overwriting an existing `report.Rmd` in `load_basic_example()`.
- Tell the user the location of the cache using a console message. Happens on every call to `get_cache(..., verbose = TRUE)`.
- Increase efficiency of internal preprocessing via `lightly_parallelize()` and `lightly_parallelize_atomic()`. Now, processing happens faster, and only over the unique values of a vector.
- Add a new `make_with_config()` function to do the work of `make()` on an existing internal configuration list from `drake_config()`.
- Add a new function `drake_batchtools_tmpl_file()` to write a `batchtools` template file from one of the examples (`drake_example()`), if one exists.

# Version 4.3.0: 2017-10-17

Version 4.3.0 has:
- Reproducible random numbers (#56)
- Automatic detection of knitr dependencies (#9)
- More vignettes
- Bug fixes

# Version 4.2.0: 2017-09-29

Version 4.2.0 will be released today. There are several improvements to code style and performance. In addition, there are new features such as cache/hash externalization and runtime prediction. See the new storage and timing vignettes for details. This release has automated checks for back-compatibility with existing projects, and I also did manual back compatibility checks on serious projects.

# Version 3.0.0: 2017-05-03

Version 3.0.0 is coming out. It manages environments more intelligently so that the behavior of `make()` is more consistent with evaluating your code in an interactive session.

# Version 1.0.1: 2017-02-28

Version 1.0.1 is on CRAN! I'm already working on a massive update, though. 2.0.0 is cleaner and more powerful.

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
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# Contributing

Development is a community effort, and we encourage participation.

## Code of Conduct

The environment for collaboration should be friendly, inclusive, respectful, and safe for everyone, so all participants must obey [this repository's code of conduct](https://github.com/ropensci/drake/blob/main/CODE_OF_CONDUCT.md).

## Issues

`drake` thrives on the suggestions, questions, and bug reports you submit to the [issue tracker](https://github.com/ropensci/drake/issues). Before posting, please search both the open and closed issues to help us avoid duplication. Usage questions are welcome, but you may also wish to post to [Stack Overflow](https://stackoverflow.com) with the [`drake-r-package` tag](https://stackoverflow.com/tags/drake-r-package).

Be considerate of the maintainer's time and make it as easy as possible to troubleshoot any problems you identify. Read [here](https://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example) and [here](https://www.tidyverse.org/help/) to learn about minimal reproducible examples. Format your code according to the [tidyverse style guide](https://style.tidyverse.org/) to make it easier for others to read.

## Development

If you would like to work on the code or documentation, please [fork this repository](https://help.github.com/articles/fork-a-repo/), make the changes in your fork, and then submit a [pull request](https://github.com/ropensci/drake/pulls). We will discuss your work and then hopefully merge it into the project.
# Summary

Please explain the context and purpose of your contribution and list the changes you made to the code base or documentation.

# Related GitHub issues and pull requests

- Ref: #

# Checklist

- [ ] I understand and agree to the [code of conduct](https://github.com/ropensci/drake/blob/main/CODE_OF_CONDUCT.md).
- [ ] I have listed any substantial changes in the [development news](https://github.com/ropensci/drake/blob/main/NEWS.md).
- [ ] I have added [`testthat`](https://github.com/r-lib/testthat) unit tests to [`tests/testthat`](https://github.com/ropensci/drake/tree/main/tests/testthat) for any new functionality.
- [ ] This pull request is not a [draft](https://github.blog/2019-02-14-introducing-draft-pull-requests).
---
name: Other
about: Something else.
title: ''
labels: ''
assignees: ''
---

## Prework

* [ ] Read and agree to the [code of conduct](https://github.com/ropensci/drake/blob/main/CODE_OF_CONDUCT.md) and [contributing guidelines](https://github.com/ropensci/drake/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/drake/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] For any problem you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please explain your thoughts clearly and concisely. If applicable, write a minimal example in R code or pseudo-code to show input, usage, and desired output.

## Reproducible example

* [ ] For any problems you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).
---
name: Bug
about: Something is wrong with drake.
title: ''
labels: 'type: bug'
assignees: wlandau
---

## Prework

* [ ] Read and agree to the [code of conduct](https://github.com/ropensci/drake/blob/main/CODE_OF_CONDUCT.md) and [contributing guidelines](https://github.com/ropensci/drake/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/drake/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Describe the bug clearly and concisely. 

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Expected result

What should have happened? Please be as specific as possible.

## Session info

End the reproducible example with a call to `sessionInfo()` in the same session (e.g. `reprex(si = TRUE)`) and include the output.
---
name: Trouble
about: Something is not working, and you want help.
title: ''
labels: 'type: trouble'
---

## Prework

* [ ] Read and agree to the [code of conduct](https://github.com/ropensci/drake/blob/main/CODE_OF_CONDUCT.md) and [contributing guidelines](https://github.com/ropensci/drake/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/drake/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Describe the trouble clearly and concisely. 

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Desired result

What output or behavior do you want to see? Please be as specific as you can.

## Session info

End the reproducible example with a call to `sessionInfo()` in the same session (e.g. `reprex(si = TRUE)`) and include the output.
---
name: New feature
about: Suggest a new feature.
title: ''
labels: 'type: new feature'
assignees: wlandau
---

## Prework

- [ ] Read and abide by `drake`'s [code of conduct](https://github.com/ropensci/drake/blob/main/CODE_OF_CONDUCT.md).
- [ ] Search for duplicates among the [existing issues](https://github.com/ropensci/drake/issues), both open and closed.

## Proposal

Describe the new feature clearly and concisely. If applicable, write a minimal example in R code or pseudo-code to show input, usage, and desired output.

To help us read any code you include (optional) please try to follow the [tidyverse style guide](https://style.tidyverse.org/). The `style_text()` and `style_file()` functions from the [`styler`](https://github.com/r-lib/styler) package make it easier.
---
name: Bottleneck
about: drake is too slow or consumes too many resources.
title: ''
labels: 'topic: performance'
assignees: wlandau

---

## Prework

* [ ] Read and agree to the [code of conduct](https://github.com/ropensci/drake/blob/main/CODE_OF_CONDUCT.md) and [contributing guidelines](https://github.com/ropensci/drake/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/drake/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Describe the bottleneck clearly and concisely. 

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Benchmarks

How poorly does `drake` perform? To find out, we recommend the [`proffer`](https://github.com/wlandau/proffer) package and take screenshots of the results displayed in your browser.

```r
library(drake)
library(proffer)
px <- pprof({
  # All your drake code goes here.
})
```

Alternatively, if installing [`proffer`](https://github.com/wlandau/proffer) is too cumbersome, create a zip archive of profiling data (e.g. `samples.zip` below) and upload it to this issue thread.

```r
Rprof(filename = "samples.rprof")
# Slow code goes here.
Rprof(NULL)
zip(zipfile = "samples.zip", files = "samples.rprof")
```
---
name: Question
about: Ask a question.
title: ''
labels: 'type: question'
assignees: ''
---

## Prework

* [ ] Read and agree to the [code of conduct](https://github.com/ropensci/drake/blob/main/CODE_OF_CONDUCT.md) and [contributing guidelines](https://github.com/ropensci/drake/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/drake/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] For any problems you identify, write a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Question

What would you like to know?

## Reproducible example

* [ ] For any problems you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).
- [ ] **Skip enough tests on CRAN to keep CRAN check time under 10 minutes**.
- [ ] Continuous integration.
- [ ] Spelling.
- [ ] Manual tests.
  - [ ] `tests/testthat/test-interactive.R`
  - [ ] `tests/testthat/test-keras.R`
- [ ] `devtools::run_examples(run = FALSE)`
- [ ] `devtools::run_examples(run = TRUE)`
- [ ] Regular tests with `tests/testthat/helper-operators.R` activated.
- [ ] Regular tests R 3.3.0.
- [ ] Test suite without Suggests packages in R 3.3.0.
- [ ] `goodpractice`
- [ ] HPC test suite on SGE.
- [ ] `tests/scenarios` on Mac, Linux, and Windows
- [ ] Update the version in the `DESCRIPTION` and in `NEWS.md`.
- [ ] Win Builder on gz file release.
- All videos are indirectly embedded and copyright remains with the respective owners of the original content.
- The files [infographic.svg](https://github.com/ropensci/drake/blob/main/docs/images/infographic.svg) and [infographic-font.svg](https://github.com/ropensci/drake/blob/main/docs/images/infographic-font.svg) were created using clipart released under the Creative Commons License:
    - "multiple" by Hea Poh Lin from the [Noun Project](https://thenounproject.com/)
    - "replay" by Sylvain A. from the [Noun Project](https://thenounproject.com/)
    - "checkmark" by Ananth from the [Noun Project](https://thenounproject.com/)
- The [tweet](https://twitter.com/fossilosophy/status/966408174470299648) from [tweet.png](https://github.com/ropensci/drake/blob/main/docs/images/tweet.png) is by [Brianna McHorse](https://github.com/bmchorse).
- The [diagram of the typical Tidyverse workflow](https://github.com/ropensci/drake/blob/main/images/tidydag.png) is from [Jenny Bryan](https://github.com/jennybc)'s [December 2017 presentation on workflow maintenance](https://speakerdeck.com/jennybc/zen-and-the-art-of-workflow-maintenance?slide=55).
---
output:
  github_document:
    html_preview: false
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
dir <- tempfile()
dir.create(dir)
knitr::opts_knit$set(root.dir = dir)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/"
)
```

```{r, echo = FALSE}
suppressMessages(suppressWarnings(library(drake)))
suppressMessages(suppressWarnings(library(dplyr)))
clean(destroy = TRUE)
invisible(drake_example("main", overwrite = TRUE))
invisible(file.copy("main/raw_data.xlsx", ".", overwrite = TRUE))
invisible(file.copy("main/report.Rmd", ".", overwrite = TRUE))
```

<center>
<img src="https://docs.ropensci.org/drake/reference/figures/infographic.svg" alt="infographic" align="center" style = "border: none; float: center;">
</center>
<table class="table">
<thead>
<tr class="header">
<th align="left">
Usage
</th>
<th align="left">
Release
</th>
<th align="left">
Development
</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">
<a href="https://www.gnu.org/licenses/gpl-3.0.en.html"><img src="https://img.shields.io/badge/licence-GPL--3-blue.svg" alt="Licence"></a>
</td>
<td align="left">
<a href="https://cran.r-project.org/package=drake"><img src="https://www.r-pkg.org/badges/version/drake" alt="CRAN"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/drake/actions?query=workflow%3Acheck"><img src="https://github.com/ropensci/drake/workflows/check/badge.svg" alt="check"></a>
</td>
</tr>
<tr class="even">
<td align="left">
<a href="https://cran.r-project.org/"><img src="https://img.shields.io/badge/R%3E%3D-3.3.0-blue.svg" alt="minimal R version"></a>
</td>
<td align="left">
<a href="https://cran.r-project.org/web/checks/check_results_drake.html"><img src="https://cranchecks.info/badges/summary/drake" alt="cran-checks"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/drake/actions?query=workflow%3Alint"><img src="https://github.com/ropensci/drake/workflows/lint/badge.svg" alt="lint"></a>
</td>
</tr>
<tr class="odd">
<td align="left">
<a href="https://CRAN.R-project.org/package=drake"><img src="https://tinyverse.netlify.com/badge/drake"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/software-review/issues/156"><img src="https://badges.ropensci.org/156_status.svg" alt="rOpenSci"></a>
</td>
<td align="left">
<a href="https://codecov.io/github/ropensci/drake?branch=main"><img src="https://codecov.io/github/ropensci/drake/coverage.svg?branch=main" alt="Codecov"></a>
</td>
</tr>
<tr class="even">
<td align="left">
<a href="https://CRAN.R-project.org/package=drake"><img src="https://cranlogs.r-pkg.org/badges/drake" alt="downloads"></a>
</td>
<td align="left">
<a href="https://doi.org/10.21105/joss.00550"><img src="https://joss.theoj.org/papers/10.21105/joss.00550/status.svg" alt="JOSS"></a>
</td>
<td align="left">
<a href="https://bestpractices.coreinfrastructure.org/projects/2135"><img src="https://bestpractices.coreinfrastructure.org/projects/2135/badge"></a>
</td>
</tr>
<tr class="odd">
<td align="left">
</td>
<td align="left">
<a href="https://zenodo.org/badge/latestdoi/82609103"><img src="https://zenodo.org/badge/82609103.svg" alt="Zenodo"></a>
</td>
<td align="left">
<a href="https://www.tidyverse.org/lifecycle/#superseded"><img src="https://img.shields.io/badge/lifecycle-superseded-blue.svg" alt='superseded lifecycle'></a>
</td>
</tr>
</tbody>
</table>
<br>

# drake is superseded. Consider targets instead.

As of 2021-01-21, `drake` is [superseded](https://www.tidyverse.org/lifecycle/#superseded). The [`targets`](https://docs.ropensci.org/targets/) R package is the long-term successor of `drake`, and it is more robust and easier to use. Please visit <https://books.ropensci.org/targets/drake.html> for full context and advice on transitioning.

# The drake R package <img src="https://docs.ropensci.org/drake/reference/figures/logo.svg" align="right" alt="logo" width="120" height = "139" style = "border: none; float: right;">

Data analysis can be slow. A round of scientific computation can take several minutes, hours, or even days to complete. After it finishes, if you update your code or data, your hard-earned results may no longer be valid. How much of that valuable output can you keep, and how much do you need to update? How much runtime must you endure all over again?

For projects in R, the `drake` package can help. It [analyzes your workflow](https://books.ropensci.org/drake/plans.html), skips steps with up-to-date results, and orchestrates the rest with [optional distributed computing](https://books.ropensci.org/drake/hpc.html). At the end, `drake` provides evidence that your results match the underlying code and data, which increases your ability to trust your research.

# Video


## That Feeling of Workflowing (Miles McBain)

<center>
<a href="https://www.youtube.com/embed/jU1Zv21GvT4">
<img src="https://docs.ropensci.org/drake/reference/figures/workflowing.png" alt="workflowing" align="center" style = "border: none; float: center;">
</a>
</center>

(By [Miles McBain](https://github.com/MilesMcBain); [venue](https://nyhackr.org/index.html),
[resources](https://github.com/MilesMcBain/nycr_meetup_talk))

## rOpenSci Community Call

<center>
<a href="https://ropensci.org/commcalls/2019-09-24/">
<img src="https://docs.ropensci.org/drake/reference/figures/commcall.png" alt="commcall" align="center" style = "border: none; float: center;">
</a>
</center>

([resources](https://ropensci.org/commcalls/2019-09-24/))



# What gets done stays done.

Too many data science projects follow a [Sisyphean loop](https://en.wikipedia.org/wiki/Sisyphus):

1. Launch the code.
2. Wait while it runs.
3. Discover an issue.
4. Rerun from scratch.

For projects with long runtimes, this process gets tedious. But with `drake`, you can automatically

1. Launch the parts that changed since last time.
2. Skip the rest.

# How it works

To set up a project, load your packages,

```{r}
library(drake)
library(dplyr)
library(ggplot2)
library(tidyr)
```

load your custom functions,

```{r}
create_plot <- function(data) {
  ggplot(data) +
    geom_histogram(aes(x = Ozone)) +
    theme_gray(24)
}
```

check any supporting files (optional),

```{r}
# Get the files with drake_example("main").
file.exists("raw_data.xlsx")
file.exists("report.Rmd")
```

and plan what you are going to do.

```{r}
plan <- drake_plan(
  raw_data = readxl::read_excel(file_in("raw_data.xlsx")),
  data = raw_data %>%
    mutate(Ozone = replace_na(Ozone, mean(Ozone, na.rm = TRUE))),
  hist = create_plot(data),
  fit = lm(Ozone ~ Wind + Temp, data),
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE
  )
)

plan
```

So far, we have just been setting the stage. Use `make()` or [`r_make()`](https://books.ropensci.org/drake/projects.html#safer-interactivity) to do the real work. Targets are built in the correct order regardless of the row order of `plan`.

```{r}
make(plan) # See also r_make().
```

Except for files like `report.html`, your output is stored in a hidden `.drake/` folder. Reading it back is easy.

```{r}
readd(data) # See also loadd().
```

You may look back on your work and see room for improvement, but it's all good! The whole point of `drake` is to help you go back and change things quickly and painlessly. For example, we forgot to give our histogram a bin width.

```{r}
readd(hist)
```

So let's fix the plotting function.

```{r}
create_plot <- function(data) {
  ggplot(data) +
    geom_histogram(aes(x = Ozone), binwidth = 10) +
    theme_gray(24)
}
```

`drake` knows which results are affected.

```{r, eval = FALSE}
vis_drake_graph(plan) # See also r_vis_drake_graph().
```

<img src="https://docs.ropensci.org/drake/reference/figures/graph.png" alt="hist1" align="center" style = "border: none; float: center;" width = "600px">

The next `make()` just builds `hist` and `report.html`. No point in wasting time on the data or model.

```{r}
make(plan) # See also r_make().
```

```{r}
loadd(hist)
hist
```

# Reproducibility with confidence

The R community emphasizes reproducibility. Traditional themes include [scientific replicability](https://en.wikipedia.org/wiki/Replication_crisis), literate programming with [knitr](https://yihui.name/knitr/), and version control with [git](https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control). But internal consistency is important too. Reproducibility carries the promise that your output matches the code and data you say you used. With the exception of [non-default triggers](https://books.ropensci.org/drake/triggers.html) and [hasty mode](https://books.ropensci.org/drake/hpc.html#hasty-mode), `drake` strives to keep this promise.

## Evidence

Suppose you are reviewing someone else's data analysis project for reproducibility. You scrutinize it carefully, checking that the datasets are available and the documentation is thorough. But could you re-create the results without the help of the original author? With `drake`, it is quick and easy to find out.

```{r}
make(plan) # See also r_make().

outdated(plan) # See also r_outdated().
```

With everything already up to date, you have **tangible evidence** of reproducibility. Even though you did not re-create the results, you know the results are recreatable. They **faithfully show** what the code is producing. Given the right [package environment](https://rstudio.github.io/packrat/) and [system configuration](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/sessionInfo.html), you have everything you need to reproduce all the output by yourself.

## Ease

When it comes time to actually rerun the entire project, you have much more confidence. Starting over from scratch is trivially easy.

```{r}
clean()    # Remove the original author's results.
make(plan) # Independently re-create the results from the code and input data.
```

## Big data efficiency

Select specialized data formats to increase speed and reduce memory consumption. In version 7.5.2.9000 and above, the available formats are ["fst"](https://github.com/fstpackage/fst) for data frames (example below) and "keras" for [Keras](https://keras.rstudio.com/) models ([example here](https://books.ropensci.org/drake/churn.html#plan)).

```{r, eval = FALSE}
library(drake)
n <- 1e8 # Each target is 1.6 GB in memory.
plan <- drake_plan(
  data_fst = target(
    data.frame(x = runif(n), y = runif(n)),
    format = "fst"
  ),
  data_old = data.frame(x = runif(n), y = runif(n))
)
make(plan)
#> target data_fst
#> target data_old
build_times(type = "build")
#> # A tibble: 2 x 4
#>   target   elapsed              user                 system    
#>   <chr>    <Duration>           <Duration>           <Duration>
#> 1 data_fst 13.93s               37.562s              7.954s    
#> 2 data_old 184s (~3.07 minutes) 177s (~2.95 minutes) 4.157s
```

## History and provenance

As of version 7.5.2, `drake` tracks the history and provenance of your targets:
what you built, when you built it, how you built it, the arguments you
used in your function calls, and how to get the data back. (Disable with `make(history = FALSE)`)

```{r}
history <- drake_history(analyze = TRUE)
history
```

Remarks:

- The `quiet` column appears above because one of the `drake_plan()` commands has `knit(quiet = TRUE)`.
- The `hash` column identifies all the previous versions of your targets. As long as `exists` is `TRUE`, you can recover old data.
- Advanced: if you use `make(cache_log_file = TRUE)` and put the cache log file under version control, you can match the hashes from `drake_history()` with the `git` commit history of your code.

Let's use the history to recover the oldest histogram.

```{r}
hash <- history %>%
  filter(target == "hist") %>%
  pull(hash) %>%
  head(n = 1)
cache <- drake_cache()
cache$get_value(hash)
```

## Independent replication

With even more evidence and confidence, you can invest the time to independently replicate the original code base if necessary. Up until this point, you relied on basic `drake` functions such as `make()`, so you may not have needed to peek at any substantive author-defined code in advance. In that case, you can stay usefully ignorant as you reimplement the original author's methodology. In other words, `drake` could potentially improve the integrity of independent replication.

## Readability and transparency

Ideally, independent observers should be able to read your code and understand it. `drake` helps in several ways.

- The [drake plan](https://docs.ropensci.org/drake/reference/drake_plan.html) explicitly outlines the steps of the analysis, and [`vis_drake_graph()`](https://docs.ropensci.org/drake/reference/vis_drake_graph.html) visualizes how those steps depend on each other.
- `drake` takes care of the parallel scheduling and high-performance computing (HPC) for you. That means the HPC code is no longer tangled up with the code that actually expresses your ideas.
- You can [generate large collections of targets](https://books.ropensci.org/drake/gsp.html) without necessarily changing your code base of imported functions, another nice separation between the concepts and the execution of your workflow

# Scale up and out.

Not every project can complete in a single R session on your laptop. Some projects need more speed or computing power. Some require a few local processor cores, and some need large high-performance computing systems. But parallel computing is hard. Your tables and figures depend on your analysis results, and your analyses depend on your datasets, so some tasks must finish before others even begin. `drake` knows what to do. Parallelism is implicit and automatic. See the [high-performance computing guide](https://books.ropensci.org/drake/hpc.html) for all the details.

```{r, eval = FALSE}
# Use the spare cores on your local machine.
make(plan, jobs = 4)

# Or scale up to a supercomputer.
drake_hpc_template_file("slurm_clustermq.tmpl") # https://slurm.schedmd.com/
options(
  clustermq.scheduler = "clustermq",
  clustermq.template = "slurm_clustermq.tmpl"
)
make(plan, parallelism = "clustermq", jobs = 4)
```

# With Docker

`drake` and Docker are compatible and complementary. Here are some examples that run `drake` inside a Docker image.

- [`drake-gitlab-docker-example`](https://gitlab.com/ecohealthalliance/drake-gitlab-docker-example): A small pedagogical example workflow that leverages `drake`, Docker, GitLab, and continuous integration in a reproducible analysis pipeline. Created by [Noam Ross](https://www.noamross.net/).
- [`pleurosoriopsis`](https://github.com/joelnitta/pleurosoriopsis): The workflow that supports [Ebihara *et al.* 2019. "Growth Dynamics of the Independent Gametophytes of *Pleurorosiopsis makinoi* (Polypodiaceae)" *Bulletin of the National Science Museum Series B (Botany)* 45:77-86.](https://www.kahaku.go.jp/research/publication/botany.html). Created by [Joel Nitta](https://github.com/joelnitta).

Alternatively, it is possible to run `drake` outside Docker and use the [`future`](https://github.com/HenrikBengtsson/future) package to send targets to a Docker image. `drake`'s [`Docker-psock`](https://github.com/wlandau/drake-examples/tree/main/Docker-psock) example demonstrates how. Download the code with `drake_example("Docker-psock")`.

# Installation

You can choose among different versions of `drake`. The CRAN release often lags behind the [online manual](https://books.ropensci.org/drake/) but may have fewer bugs.

```{r, eval = FALSE}
# Install the latest stable release from CRAN.
install.packages("drake")

# Alternatively, install the development version from GitHub.
install.packages("devtools")
library(devtools)
install_github("ropensci/drake")
```

# Function reference

The [reference section](https://docs.ropensci.org/drake/reference/index.html) lists all the available functions. Here are the most important ones.

- `drake_plan()`: create a workflow data frame (like `my_plan`).
- `make()`: build your project.
- `drake_history()`: show what you built, when you built it, and the function arguments you used.
- `r_make()`: launch a fresh [`callr::r()`](https://github.com/r-lib/callr) process to build your project. Called from an interactive R session, `r_make()` is more reproducible than `make()`.
- `loadd()`: load one or more built targets into your R session.
- `readd()`: read and return a built target.
- `vis_drake_graph()`: show an interactive visual network representation of your workflow.
- `recoverable()`: Which targets can we salvage using `make(recover = TRUE)` (experimental).
- `outdated()`: see which targets will be built in the next `make()`.
- `deps_code()`: check the dependencies of a command or function.
- `drake_failed()`: list the targets that failed to build in the last `make()`.
- `diagnose()`: return the full context of a build, including errors, warnings, and messages.

# Documentation

## Core concepts

The following resources explain what `drake` can do and how it works. The [`learndrake`](https://github.com/wlandau/learndrake) workshop devotes particular attention to `drake`'s mental model.

- The [user manual](https://books.ropensci.org/drake/)
- [`drakeplanner`](https://github.com/wlandau/drakeplanner), an R/Shiny app to help learn `drake` and create new projects. Run locally with `drakeplanner::drakeplanner()` or access it at <https://wlandau.shinyapps.io/drakeplanner>.
- [`learndrake`](https://github.com/wlandau/learndrake), an R package for teaching an extended `drake` workshop. It contains notebooks, slides, Shiny apps, the latter two of which are publicly deployed. See the [README](https://github.com/wlandau/learndrake/blob/main/README.md) for instructions and links.

## In practice

- [Miles McBain](https://github.com/MilesMcBain)'s [excellent blog post](https://milesmcbain.xyz/the-drake-post/) explains the motivating factors and practical issues {drake} solves for most projects, how to set up a project as quickly and painlessly as possible, and how to overcome common obstacles.
- Miles' [`dflow`](https://github.com/MilesMcBain/dflow) package generates the file structure for a boilerplate `drake` project. It is a more thorough alternative to `drake::use_drake()`.
- `drake` is heavily function-oriented by design, and Miles' [`fnmate`](https://github.com/MilesMcBain/fnmate) package automatically generates boilerplate code and docstrings for functions you mention in `drake` plans.

## Reference

- The [reference website](https://docs.ropensci.org/drake/).
- The [official repository of example code](https://github.com/wlandau/drake-examples). Download an example workflow from here with `drake_example()`.
- Presentations and workshops by [Will Landau](https://github.com/wlandau), [Kirill Müller](https://github.com/krlmlr), [Amanda Dobbyn](https://github.com/aedobbyn), [Karthik Ram](https://github.com/karthik), [Sina Rüeger](https://github.com/sinarueeger), [Christine Stawitz](https://github.com/cstawitz), and others. See specific links at <https://books.ropensci.org/drake/index.html#presentations>
- The [FAQ page](https://books.ropensci.org/drake/faq.html), which links to [appropriately-labeled issues on GitHub](https://github.com/ropensci/drake/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3A%22frequently+asked+question%22+).

## Use cases

The official [rOpenSci use cases](https://discuss.ropensci.org/c/usecases) and [associated discussion threads](https://discuss.ropensci.org/c/usecases) describe applications of `drake` in the real world. Many of these use cases are linked from the [`drake` tag on the rOpenSci discussion forum](https://discuss.ropensci.org/tag/drake).

Here are some additional applications of `drake` in real-world projects.

- [efcaguab/demografia-del-voto](https://github.com/efcaguab/demografia-del-voto)
- [efcaguab/great-white-shark-nsw](https://github.com/efcaguab/great-white-shark-nsw)
- [IndianaCHE/Detailed-SSP-Reports](https://github.com/IndianaCHE/Detailed-SSP-Reports)
- [joelnitta/pleurosoriopsis](https://github.com/joelnitta/pleurosoriopsis)
- [pat-s/pathogen-modeling](https://github.com/pat-s/pathogen-modeling)
- [sol-eng/tensorflow-w-r](https://github.com/sol-eng/tensorflow-w-r)
- [tiernanmartin/home-and-hope](https://github.com/tiernanmartin/home-and-hope)

## `drake` projects as R packages

Some folks like to structure their `drake` workflows as R packages. Examples are below. In your own analysis packages, be sure to call `drake::expose_imports(yourPackage)` so `drake` can watch you package's functions for changes and rebuild downstream targets accordingly.

- [b-rodrigues/coolmlproject](https://github.com/b-rodrigues/coolmlproject)
- [tiernanmartin/drakepkg](https://github.com/tiernanmartin/drakepkg)

# Help and troubleshooting

The following resources document many known issues and challenges.

- [Frequently-asked questions](https://github.com/ropensci/drake/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3A%22frequently+asked+question%22+).
- [Debugging and testing drake projects](https://books.ropensci.org/drake/debugging.html)
- [Other known issues](https://github.com/ropensci/drake/issues) (please search both open and closed ones).

If you are still having trouble, please submit a [new issue](https://github.com/ropensci/drake/issues/new) with a bug report or feature request, along with a minimal reproducible example where appropriate.

The GitHub issue tracker is mainly intended for bug reports and feature requests. While questions about usage etc. are also highly encouraged, you may alternatively wish to post to [Stack Overflow](https://stackoverflow.com) and use the [`drake-r-package` tag](https://stackoverflow.com/tags/drake-r-package).

# Contributing

Development is a community effort, and we encourage participation. Please read [CONTRIBUTING.md](https://github.com/ropensci/drake/blob/main/CONTRIBUTING.md) for details.

# Similar work

`drake` enhances reproducibility and high-performance computing, but not in all respects. [Literate programming](https://rmarkdown.rstudio.com/), [local library managers](https://rstudio.github.io/packrat), [containerization](https://www.docker.com/), and [strict session managers](https://github.com/tidyverse/reprex) offer more robust solutions in their respective domains. And for the problems `drake` *does* solve, it stands on the shoulders of the giants that came before.

## Pipeline tools

### GNU Make

The original idea of a time-saving reproducible build system extends back at least as far as [GNU Make](https://www.gnu.org/software/make/), which still aids the work of [data scientists](http://blog.kaggle.com/2012/10/15/make-for-data-scientists/) as well as the original user base of complied language programmers. In fact, the name "drake" stands for "Data Frames in R for Make". [Make](https://kbroman.org/minimal_make/) is used widely in reproducible research. Below are some examples from [Karl Broman's website](https://kbroman.org/minimal_make/).

- Bostock, Mike (2013). "A map of flowlines from NHDPlus." https://github.com/mbostock/us-rivers. Powered by the Makefile at https://github.com/mbostock/us-rivers/blob/master/Makefile.
- Broman, Karl W (2012). "Halotype Probabilities in Advanced Intercross Populations." *G3* 2(2), 199-202.Powered by the `Makefile` at https://github.com/kbroman/ailProbPaper/blob/master/Makefile.
- Broman, Karl W (2012). "Genotype Probabilities at Intermediate Generations in the Construction of Recombinant Inbred Lines." *Genetics 190(2), 403-412. Powered by the Makefile at https://github.com/kbroman/preCCProbPaper/blob/master/Makefile.
- Broman, Karl W and Kim, Sungjin and Sen, Saunak and Ane, Cecile and Payseur, Bret A (2012). "Mapping Quantitative Trait Loci onto a Phylogenetic Tree." *Genetics* 192(2), 267-279. Powered by the `Makefile` at https://github.com/kbroman/phyloQTLpaper/blob/master/Makefile.

Whereas [GNU Make](https://www.gnu.org/software/make/) is language-agnostic, `drake` is fundamentally designed for R.

- Instead of a [Makefile](https://github.com/kbroman/preCCProbPaper/blob/master/Makefile), `drake` supports an R-friendly [domain-specific language](https://books.ropensci.org/drake/plans.html#large-plans) for declaring targets.
- Targets in [GNU Make](https://www.gnu.org/software/make/) are files, whereas targets in `drake` are arbitrary variables in memory. (`drake` does have opt-in support for files via `file_out()`, `file_in()`, and `knitr_in()`.) `drake` caches these objects in its own [storage system](https://github.com/richfitz/storr) so R users rarely have to think about output files.

### Remake

[remake](https://github.com/richfitz/remake) itself is no longer maintained, but its founding design goals and principles live on through [drake](https://github.com/ropensci/drake). In fact, [drake](https://github.com/ropensci/drake) is a direct re-imagining of [remake](https://github.com/richfitz/remake) with enhanced scalability, reproducibility, high-performance computing, visualization, and documentation.

### Factual's Drake

[Factual's Drake](https://github.com/Factual/drake) is similar in concept, but the development effort is completely unrelated to the [drake R package](https://github.com/ropensci/drake).

### Other pipeline tools

There are [countless other successful pipeline toolkits](https://github.com/pditommaso/awesome-pipeline). The `drake` package distinguishes itself with its R-focused approach, Tidyverse-friendly interface, and a [thorough selection of parallel computing technologies and scheduling algorithms](https://books.ropensci.org/drake/hpc.html).

## Memoization

Memoization is the strategic caching of the return values of functions. It is a lightweight approach to the core problem that `drake` and other pipeline tools are trying to solve. Every time a memoized function is called with a new set of arguments, the return value is saved for future use. Later, whenever the same function is called with the same arguments, the previous return value is salvaged, and the function call is skipped to save time. The [`memoise`](https://github.com/r-lib/memoise) package is the primary implementation of memoization in R.

Memoization saves time for small projects, but it arguably does not go far enough for large reproducible pipelines. In reality, the return value of a function depends not only on the function body and the arguments, but also on any nested functions and global variables, the dependencies of those dependencies, and so on upstream. `drake` tracks this deeper context, while [memoise](https://github.com/r-lib/memoise) does not.

## Literate programming

[Literate programming](https://rmarkdown.rstudio.com/) is the practice of narrating code in plain vernacular. The goal is to communicate the research process clearly, transparently, and reproducibly. Whereas commented code is still mostly code, literate [knitr](https://yihui.name/knitr/) / [R Markdown](https://rmarkdown.rstudio.com/) reports can become websites, presentation slides, lecture notes, serious scientific manuscripts, and even books.

### knitr and R Markdown

`drake` and [knitr](https://yihui.name/knitr/) are symbiotic. `drake`'s job is to manage large computation and orchestrate the demanding tasks of a complex data analysis pipeline. [knitr](https://yihui.name/knitr/)'s job is to communicate those expensive results after `drake` computes them. [knitr](https://yihui.name/knitr/) / [R Markdown](https://rmarkdown.rstudio.com/) reports are small pieces of an overarching `drake` pipeline. They should focus on communication, and they should do as little computation as possible. 

To insert a [knitr](https://yihui.name/knitr/) report in a `drake` pipeline, use the `knitr_in()` function inside your [`drake` plan](https://books.ropensci.org/drake/plans.html), and use `loadd()` and `readd()` to refer to targets in the report itself. See an [example here](https://github.com/wlandau/drake-examples/tree/main/main).

### Version control

`drake` is not a version control tool. However, it is fully compatible with [`git`](https://git-scm.com/), [`svn`](https://en.wikipedia.org/wiki/Apache_Subversion), and similar software. In fact, it is good practice to use [`git`](https://git-scm.com/) alongside `drake` for reproducible workflows.

However, data poses a challenge. The datasets created by `make()` can get large and numerous, and it is not recommended to put the `.drake/` cache or the `.drake_history/` logs under version control. Instead, it is recommended to use a data storage solution such as DropBox or [OSF](https://osf.io/ka7jv/wiki/home/).

### Containerization and R package environments

`drake` does not track R packages or system dependencies for changes. Instead, it defers to tools like [Docker](https://www.docker.com), [Singularity](https://sylabs.io/singularity/), [`renv`](https://github.com/rstudio/renv), and [`packrat`](https://github.com/rstudio/packrat), which create self-contained portable environments to reproducibly isolate and ship data analysis projects. `drake` is fully compatible with these tools.

### workflowr

The [`workflowr`](https://github.com/jdblischak/workflowr) package is a project manager that focuses on literate programming, sharing over the web, file organization, and version control. Its brand of reproducibility is all about transparency, communication, and discoverability. For an example of [`workflowr`](https://github.com/jdblischak/workflowr) and `drake` working together, see [this machine learning project](https://github.com/pat-s/2019-feature-selection) by [Patrick Schratz](https://github.com/pat-s).

# Citation

```{r}
citation("drake")
```

# Acknowledgements

Special thanks to [Jarad Niemi](https://www.jarad.me/), my advisor from [graduate school](https://stat.iastate.edu/), for first introducing me to the idea of [Makefiles](https://www.gnu.org/software/make/) for research. He originally set me down the path that led to `drake`.

Many thanks to [Julia Lowndes](https://github.com/jules32), [Ben Marwick](https://github.com/benmarwick), and [Peter Slaughter](https://github.com/gothub) for [reviewing drake for rOpenSci](https://github.com/ropensci/software-review/issues/156), and to [Maëlle Salmon](https://github.com/maelle) for such active involvement as the editor. Thanks also to the following people for contributing early in development.

- [Alex Axthelm](https://github.com/AlexAxthelm)
- [Chan-Yub Park](https://github.com/mrchypark)
- [Daniel Falster](https://github.com/dfalster)
- [Eric Nantz](https://github.com/rpodcast)
- [Henrik Bengtsson](https://github.com/HenrikBengtsson)
- [Ian Watson](https://github.com/IanAWatson)
- [Jasper Clarkberg](https://github.com/dapperjapper)
- [Kendon Bell](https://github.com/kendonB)
- [Kirill M&uuml;ller](https://github.com/krlmlr)
- [Michael Schubert](https://github.com/mschubert)

Credit for images is [attributed here](https://github.com/ropensci/drake/blob/main/man/figures/image-credit.md).

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "Empty report"
output:
  github_document: default
  html_document:
    toc: true
    toc_float: true
---

Lorem ipsum flotsam jetsam.

```{r nested_object}
loadd(nested)
```
---
title: "Bad report"
output:
  github_document: default
  html_document:
    toc: true
    toc_float: true
---

```{r bad-code-chunk}
Lorem ipsum flotsam jetsam.
```
---
title: "Test Report"
author: "Will Landau"
date: "October 3, 2017"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

`drake` should be able to analyze this report and detect the right dependencies
using static code analysis on calls to `loadd()` and `readd()`.

The value `r readd(inline_dep) %>% broom::tidy() %>% as.tibble %>% pull(p.value) %>% first` was calculated inline.

```r
x <- readd(should_not_find)
```

```{r dry, eval = FALSE}
x <- readd(should_not_find)
```

```{r chunk}
library(drake)
var <- 10
var2 <- list(17, var, "ignore", "ignore2")

out <- readd(target1, character_only = dont_detect_this)
loadd(target2)
out <- drake::readd(character_only = FALSE, target = target3)
drake::loadd(target4, list = "target5")
out <- drake:::readd(target = "target6")
drake:::loadd(target7, target8)
drake:::loadd(target9, list = c("target10", "target11"), nothing = to_see)
loadd(list = c("target12", "target13"))

out <- drasadfke::readd(ignore1)
drawe::loadd(ig2, ignore3)
out <- rake:::readd(target_ignore)
make:::loadd(nothing, to, see, list = "here")

f(readd(target14) + var)
f(readd(target15) + var) # deliberate repeat
f(g(drake::readd(target16) + var))
f(readd(target = "target17", character_only = TRUE) + var)
g <- function(){
  f(drake:::loadd(target18, character_only = IGNORE_THIS) + var)
}
function(){
  f(drake:::loadd(target18, character_only = IGNORE_THIS) + var)
}
readd("\"file1\"")
readd(target = "\"file2\"")
readd(target = "\"file3\"", character_only = TRUE)
readd(target = "\"file4\"", character_only = FALSE)
loadd("\"file5\"")
loadd(list = "\"file6\"")
```

```{r file_deps}
file_in("input.txt")
file_out("output.txt")
knitr_in("nested.Rmd")
```
---
title: "Empty report"
output:
  github_document: default
  html_document:
    toc: true
    toc_float: true
---

Lorem ipsum flotsam jetsam.
---
title: "Example R Markdown drake file target"
author: Will Landau and Kirill Müller
output: html_document
---

# Content

```{r content}
tryCatch({
  library(drake)
  loadd(fit)
  print(fit)
  readd(hist)
},
error = function(e){
  stop("please read the instructions in the R Markdown file.")
})
```

# Instructions

This file `report.Rmd` belongs to the example from [this chapter of the user manual](https://books.ropensci.org/drake/intro.html) and `drake_example("main")`. The report does not compile from scratch because it is a [file target](https://docs.ropensci.org/drake/reference/file_out.html) in a [`drake`](https://github.com/ropensci/drake)-powered reproducible workflow ([details here](https://books.ropensci.org/drake/plans.html)).

To compile this report,

1. Name it `report.Rmd` (or modify the `knitr_in("report.Rmd")` line in [`make.R`](https://github.com/ropensci/drake/tree/main/inst/rmarkdown/templates/drake/skeleton/make.R) accordingly).
2. Make sure the included file [`raw_data.xlsx`](https://github.com/ropensci/drake/tree/main/inst/rmarkdown/templates/drake/skeleton/raw_data.xlsx) is in your current working directory.
2. Install the [`drake`](https://github.com/ropensci/drake) and [`tidyverse`](https://github.com/tidyverse/tidyverse) packages. 
3. Use [`make.R`](https://github.com/ropensci/drake/tree/main/inst/rmarkdown/templates/drake/skeleton/make.R) (included with this template) to run the [`drake`](https://github.com/ropensci/drake) workflow that compiles this report.

Step 3 not only generates the output file `report.hml`, but also produces a `.drake/` cache in your working directory, which enables `report.Rmd` to compile on its own with with [RStudio](https://www.rstudio.com/products/RStudio/), [`knitr`](https://github.com/yihui/knitr), or [`rmarkdown`](https://github.com/rstudio/rmarkdown). A great way way to generate `report.hml` is still `make(plan)`. That way, compilation happens if and only if there were changes to `report.Rmd`, `report.hml`, `fit`, or `hist` since the last `make()`.
---
title: "Final results report for the mtcars example"
author: You
output: html_document
---

# The weight and fuel efficiency of cars

Is there an association between the weight and the fuel efficiency of cars? To find out, we use the `mtcars` dataset from the `datasets` package. The `mtcars` data originally came from the 1974 Motor Trend US magazine, and it contains design and performance data on 32 models of automobile.

```{r mtcars}
# ?mtcars # more info
head(mtcars)
```

Here, `wt` is weight in tons, and `mpg` is fuel efficiency in miles per gallon. We want to figure out if there is an association between `wt` and `mpg`. The `mtcars` dataset itself only has 32 rows, so we generated two larger bootstrapped datasets. We called them `small` and `large`.

```{r load_datasets}
head(drake::readd(small)) # 48 rows
drake::loadd(large)       # 64 rows
head(large)
```

Then, we fit a couple regression models to the `small` and `large` to try to detect an association between `wt` and `mpg`. Here are the coefficients and p-values from one of the model fits.

```{r load_coef}
drake::readd(coef_regression2_small)
```

Since the p-value on `x2` is so small, there may be an association between weight and fuel efficiency after all.

# A note on knitr reports in drake projects.

Because of the calls to `readd()` and `loadd()`, `drake` knows that `small`, `large`, and `coef_regression2_small` are dependencies of this R Markdown report. This dependency relationship is what causes the report to be processed at the very end.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{drake_session}
\alias{drake_session}
\title{Session info of the last call to \code{\link[=make]{make()}}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
drake_session(
  path = getwd(),
  search = TRUE,
  cache = drake::get_cache(path = path, search = search, verbose = verbose),
  verbose = 1L
)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\value{
\code{\link[=sessionInfo]{sessionInfo()}} of the last call to \code{\link[=make]{make()}}
}
\description{
Deprecated. Use \code{\link[=drake_get_session_info]{drake_get_session_info()}} instead.
}
\details{
Deprecated on 2018-12-06.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio.R
\name{rs_addin_r_make}
\alias{rs_addin_r_make}
\title{RStudio addin for r_make()
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
rs_addin_r_make(r_args = list())
}
\arguments{
\item{r_args}{List of arguments to \code{r_fn}, not including \code{func} or \code{args}.
Example:
\code{r_make(r_fn = callr::r_bg, r_args = list(stdout = "stdout.log"))}.}
}
\value{
Nothing.
}
\description{
Call \code{\link[=r_make]{r_make()}} in an RStudio addin.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{expand_plan}
\alias{expand_plan}
\title{Deprecated: create replicates of targets.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
expand_plan(plan, values = NULL, rename = TRUE, sep = "_", sanitize = TRUE)
}
\arguments{
\item{plan}{Workflow plan data frame.}

\item{values}{Values to expand over. These will be appended to
the names of the new targets.}

\item{rename}{Logical, whether to rename the targets
based on the \code{values}. See the examples for a demo.}

\item{sep}{Character scalar, delimiter between the original
target names and the values to append to create the new
target names. Only relevant when \code{rename} is \code{TRUE}.}

\item{sanitize}{Logical, whether to sanitize the plan.}
}
\value{
An expanded workflow plan data frame (with replicated targets).
}
\description{
Deprecated on 2019-05-16. Use \code{\link[=drake_plan]{drake_plan()}}
transformations instead. See
\url{https://books.ropensci.org/drake/plans.html#large-plans}
for the details.
}
\details{
Duplicates the rows of a workflow plan data frame.
Prefixes are appended to the new target names
so targets still have unique names.
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_build.R
\name{drake_debug}
\alias{drake_debug}
\title{Run a single target's command in debug mode.'
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#questioning}{\figure{lifecycle-questioning.svg}{options: alt='[Questioning]'}}}{\strong{[Questioning]}}}
\usage{
drake_debug(
  target = NULL,
  ...,
  character_only = FALSE,
  replace = FALSE,
  verbose = TRUE,
  config = NULL
)
}
\arguments{
\item{target}{Name of the target.}

\item{...}{Arguments to \code{\link[=make]{make()}}, such as the plan and environment.}

\item{character_only}{Logical, whether \code{name} should be treated
as a character or a symbol
(just like \code{character.only} in \code{\link[=library]{library()}}).}

\item{replace}{Logical. If \code{FALSE},
items already in your environment
will not be replaced.}

\item{verbose}{Logical, whether to print out the target
you are debugging.}

\item{config}{Deprecated 2019-12-22.}
}
\value{
The value of the target right after it is built.
}
\description{
Not valid for dynamic branching.
\code{drake_debug()} loads a target's dependencies
and then runs its command in debug mode (see \code{browser()},
\code{debug()}, and \code{debugonce()}). This function does not
store the target's value in the cache
(see \url{https://github.com/ropensci/drake/issues/587}).
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
# This example is not really a user-side demonstration.
# It just walks through a dive into the internals.
# Populate your workspace and write 'report.Rmd'.
load_mtcars_example() # Get the code with drake_example("mtcars").
# out <- drake_debug(small, my_plan)
# `small` was invisibly returned.
# head(out)
}
})
}
}
\seealso{
\code{\link[=drake_build]{drake_build()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{cache_namespaces}
\alias{cache_namespaces}
\title{List all the \code{storr} cache namespaces used by drake.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
cache_namespaces(default = storr::storr_environment()$default_namespace)
}
\arguments{
\item{default}{Name of the default \code{storr} namespace.}
}
\value{
A character vector of \code{storr} namespaces used for drake.
}
\description{
Deprecated on 2019-01-12.
}
\details{
Ordinary users do not need to worry about this function.
It is just another window into \code{drake}'s internals.
}
\seealso{
\code{\link[=make]{make()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{in_progress}
\alias{in_progress}
\title{List the targets in progress
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
in_progress(
  path = getwd(),
  search = TRUE,
  cache = drake::get_cache(path = path, search = search, verbose = verbose),
  verbose = 1L
)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\value{
A character vector of target names.
}
\description{
Deprecated on 2019-01-13.
}
\details{
Similar to \code{\link[=progress]{progress()}}.
}
\seealso{
\code{\link[=diagnose]{diagnose()}}, \code{\link[=drake_get_session_info]{drake_get_session_info()}},
\code{\link[=cached]{cached()}}, \code{\link[=readd]{readd()}}, \code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{find_knitr_doc}
\alias{find_knitr_doc}
\title{find_knitr_doc \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
find_knitr_doc(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{analysis_wildcard}
\alias{analysis_wildcard}
\title{Show the analysis wildcard
used in \code{\link[=plan_summaries]{plan_summaries()}}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
analysis_wildcard()
}
\value{
The analysis wildcard used in \code{\link[=plan_summaries]{plan_summaries()}}.
}
\description{
Deprecated on 2019-01-12.
}
\details{
Used to generate workflow plan data frames.
}
\seealso{
\code{\link[=plan_summaries]{plan_summaries()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{drake_unquote}
\alias{drake_unquote}
\title{Remove leading and trailing
escaped quotes from character strings.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
drake_unquote(x = NULL)
}
\arguments{
\item{x}{Character vector.}
}
\value{
Character vector without leading
or trailing escaped quotes around
the elements.
}
\description{
Deprecated on 2019-01-01
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{isolate_example}
\alias{isolate_example}
\title{Isolate the side effects of an example.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
isolate_example(desc, code)
}
\arguments{
\item{desc}{Character, description of the example.}

\item{...}{Code to run.}
}
\value{
Nothing.
}
\description{
Runs code in a temporary directory
in a controlled environment with a controlled
set of options.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hpc.R
\name{drake_hpc_template_file}
\alias{drake_hpc_template_file}
\title{Write a template file for deploying
work to a cluster / job scheduler.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_hpc_template_file(
  file = drake::drake_hpc_template_files(),
  to = getwd(),
  overwrite = FALSE
)
}
\arguments{
\item{file}{Name of the template file, including the "tmpl" extension.}

\item{to}{Character vector, where to write the file.}

\item{overwrite}{Logical, whether to overwrite an existing file of the
same name.}
}
\value{
\code{NULL} is returned,
but a batchtools template file is written.
}
\description{
See the example files from
\code{\link[=drake_examples]{drake_examples()}} and \code{\link[=drake_example]{drake_example()}}
for example usage.
}
\examples{
\dontrun{
plan <- drake_plan(x = rnorm(1e7), y = rnorm(1e7))
# List the available template files.
drake_hpc_template_files()
# Write a SLURM template file.
out <- file.path(tempdir(), "slurm_batchtools.tmpl")
drake_hpc_template_file("slurm_batchtools.tmpl", to = tempdir())
cat(readLines(out), sep = "\n")
# library(future.batchtools) # nolint
# future::plan(batchtools_slurm, template = out) # nolint
# make(plan, parallelism = "future", jobs = 2) # nolint
}
}
\seealso{
\code{\link[=drake_hpc_template_files]{drake_hpc_template_files()}},
\code{\link[=drake_examples]{drake_examples()}}, \code{\link[=drake_example]{drake_example()}},
\code{\link[=shell_file]{shell_file()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{progress}
\alias{progress}
\title{Get the build progress of your targets
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
progress(
  ...,
  list = character(0),
  no_imported_objects = NULL,
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  verbose = 1L,
  jobs = 1,
  progress = NULL
)
}
\arguments{
\item{...}{Objects to load from the cache, as names (unquoted)
or character strings (quoted). If the \code{tidyselect} package is installed,
you can also supply \code{dplyr}-style \code{tidyselect}
commands such as \code{starts_with()}, \code{ends_with()}, and \code{one_of()}.}

\item{list}{Character vector naming objects to be loaded from the
cache. Similar to the \code{list} argument of \code{\link[=remove]{remove()}}.}

\item{no_imported_objects}{Logical, whether to only return information
about imported files and targets with commands (i.e. whether to ignore
imported objects that are not files).}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{jobs}{Number of jobs/workers for parallel processing.}

\item{progress}{Character vector for filtering the build progress results.
Defaults to \code{NULL} (no filtering) to report progress of all objects.
Supported filters are \code{"done"}, \code{"running"}, and \code{"failed"}.}
}
\value{
The build progress of each target reached by
the current \code{\link[=make]{make()}} so far.
}
\description{
Deprecated on 2020-03-23. Use \code{\link[=drake_progress]{drake_progress()}} instead.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/backend_clustermq.R
\name{cmq_build}
\alias{cmq_build}
\title{Build a target using the clustermq backend
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
cmq_build(target, meta, deps, spec, config_tmp, config)
}
\arguments{
\item{target}{Target name.}

\item{meta}{List of metadata.}

\item{deps}{Named list of target dependencies.}

\item{spec}{Internal, part of the full \code{config$spec}.}

\item{config_tmp}{Internal, extra parts of \code{config} that the workers need.}

\item{config}{A \code{\link[=drake_config]{drake_config()}} list.}
}
\description{
For internal use only
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{parallel_stages}
\alias{parallel_stages}
\title{parallel_stages \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
parallel_stages(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{static_drake_graph}
\alias{static_drake_graph}
\title{Deprecated: show a \code{ggraph}/\code{ggplot2} representation
of your drake project.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
static_drake_graph(
  config,
  build_times = "build",
  digits = 3,
  targets_only = FALSE,
  main = NULL,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  make_imports = TRUE,
  from_scratch = FALSE,
  full_legend = FALSE,
  group = NULL,
  clusters = NULL
)
}
\arguments{
\item{config}{Deprecated.}

\item{build_times}{Character string or logical.
If character, the choices are
1. \code{"build"}: runtime of the command plus the time
it take to store the target or import.
2. \code{"command"}: just the runtime of the command.
3. \code{"none"}: no build times.
If logical, \code{build_times} selects whether to show the
times from `build_times(..., type = "build")`` or use
no build times at all. See \code{\link[=build_times]{build_times()}} for details.}

\item{digits}{Number of digits for rounding the build times}

\item{targets_only}{Logical,
whether to skip the imports and only include the
targets in the workflow plan.}

\item{main}{Character string, title of the graph.}

\item{from}{Optional collection of target/import names.
If \code{from} is nonempty,
the graph will restrict itself to
a neighborhood of \code{from}.
Control the neighborhood with
\code{mode} and \code{order}.}

\item{mode}{Which direction to branch out in the graph
to create a neighborhood around \code{from}.
Use \code{"in"} to go upstream,
\code{"out"} to go downstream,
and \code{"all"} to go both ways and disregard
edge direction altogether.}

\item{order}{How far to branch out to create
a neighborhood around \code{from}. Defaults to
as far as possible. If a target is in the neighborhood, then
so are all of its custom \code{\link[=file_out]{file_out()}} files if
\code{show_output_files} is \code{TRUE}.
That means the actual graph order may be slightly greater than
you might expect, but this ensures consistency
between \code{show_output_files = TRUE} and
\code{show_output_files = FALSE}.}

\item{subset}{Optional character vector.
Subset of targets/imports to display in the graph.
Applied after \code{from}, \code{mode}, and \code{order}.
Be advised: edges are only kept for adjacent nodes in \code{subset}.
If you do not select all the intermediate nodes,
edges will drop from the graph.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to increase speed and risk using obsolete information.}

\item{from_scratch}{Logical, whether to assume all the targets
will be made from scratch on the next \code{\link[=make]{make()}}.
Makes all targets outdated, but keeps information about
build progress in previous \code{\link[=make]{make()}}s.}

\item{full_legend}{Logical. If \code{TRUE}, all the node types
are printed in the legend. If \code{FALSE}, only the
node types used are printed in the legend.}

\item{group}{Optional character scalar, name of the column used to
group nodes into columns. All the columns names of your original \code{drake}
plan are choices. The other choices (such as \code{"status"}) are column names
in the \code{nodes} . To group nodes into clusters in the graph,
you must also supply the \code{clusters} argument.}

\item{clusters}{Optional character vector of values to cluster on.
These values must be elements of the column of the \code{nodes} data frame
that you specify in the \code{group} argument to \code{drake_graph_info()}.}
}
\value{
A \code{ggplot2} object, which you can modify with more layers,
show with \code{plot()}, or save as a file with \code{ggsave()}.
}
\description{
Use \code{\link[=drake_ggraph]{drake_ggraph()}} instead.
}
\details{
Deprecated on 2018-07-25.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{default_verbose}
\alias{default_verbose}
\title{Default verbosity
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
default_verbose()
}
\value{
1
}
\description{
Deprecated on 2019-01-01
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_drake_deps_ht.R
\name{new_drake_deps_ht}
\alias{new_drake_deps_ht}
\title{\code{drake_deps_ht} constructor}
\usage{
new_drake_deps_ht(
  globals = ht_new(hash = TRUE),
  namespaced = ht_new(hash = FALSE),
  strings = ht_new(hash = FALSE),
  loadd = ht_new(hash = FALSE),
  readd = ht_new(hash = FALSE),
  file_in = ht_new(hash = FALSE),
  file_out = ht_new(hash = FALSE),
  knitr_in = ht_new(hash = FALSE)
)
}
\value{
A \code{drake_deps_ht} object.
}
\description{
List of class \code{drake_deps_ht}.
}
\examples{
if (FALSE) { # stronger than roxygen dontrun
new_drake_deps_ht()
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{as_drake_filename}
\alias{as_drake_filename}
\title{as_drake_filename \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
as_drake_filename(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outdated.R
\name{outdated_impl}
\alias{outdated_impl}
\title{Internal function with a drake_config() argument}
\usage{
outdated_impl(config, make_imports = TRUE, do_prework = TRUE)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to save some time and risk obsolete output.}

\item{do_prework}{Whether to do the \code{prework}
normally supplied to \code{\link[=make]{make()}}.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{dependency_profile}
\alias{dependency_profile}
\title{States of the dependencies of a target
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
dependency_profile(target, config, character_only = FALSE)
}
\arguments{
\item{target}{Name of the target.}

\item{config}{Deprecated.}

\item{character_only}{Logical, whether to assume \code{target}
is a character string rather than a symbol.}
}
\value{
A data frame of the old hashes and
new hashes of the data frame, along with
an indication of which hashes changed since
the last \code{\link[=make]{make()}}.
}
\description{
Deprecated on 2019-02-14.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{plot_graph}
\alias{plot_graph}
\title{plot_graph \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
plot_graph(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{dataset_wildcard}
\alias{dataset_wildcard}
\title{Show the dataset wildcard
used in \code{\link[=plan_analyses]{plan_analyses()}} and \code{\link[=plan_summaries]{plan_summaries()}}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
dataset_wildcard()
}
\value{
The dataset wildcard used in
\code{\link[=plan_analyses]{plan_analyses()}} and \code{\link[=plan_summaries]{plan_summaries()}}.
}
\description{
Deprecated on 2019-01-12.
}
\details{
Used to generate workflow plan data frames.
}
\seealso{
\code{\link[=plan_analyses]{plan_analyses()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{read_drake_config}
\alias{read_drake_config}
\title{Read a config object from the cache
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
read_drake_config(
  path = getwd(),
  search = TRUE,
  cache = NULL,
  verbose = 1L,
  jobs = 1,
  envir = parent.frame()
)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{jobs}{Number of jobs/workers for parallel processing.}
}
\description{
drake no longer stores the config object,
the plan, etc. in the cache during \code{make()}. This change
improves speed.
}
\details{
2019-01-06
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan.R
\name{drake_plan}
\alias{drake_plan}
\title{Create a drake plan
for the \code{plan} argument of \code{\link[=make]{make()}}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_plan(
  ...,
  list = NULL,
  file_targets = NULL,
  strings_in_dots = NULL,
  tidy_evaluation = NULL,
  transform = TRUE,
  trace = FALSE,
  envir = parent.frame(),
  tidy_eval = TRUE,
  max_expand = NULL
)
}
\arguments{
\item{...}{A collection of symbols/targets
with commands assigned to them. See the examples for details.}

\item{list}{Deprecated}

\item{file_targets}{Deprecated.}

\item{strings_in_dots}{Deprecated.}

\item{tidy_evaluation}{Deprecated. Use \code{tidy_eval} instead.}

\item{transform}{Logical, whether to transform the plan
into a larger plan with more targets.
Requires the \code{transform} field in
\code{target()}. See the examples for details.}

\item{trace}{Logical, whether to add columns to show
what happens during target transformations.}

\item{envir}{Environment for tidy evaluation.}

\item{tidy_eval}{Logical, whether to use tidy evaluation
(e.g. unquoting/\verb{!!}) when resolving commands.
Tidy evaluation in transformations is always turned on
regardless of the value you supply to this argument.}

\item{max_expand}{Positive integer, optional.
\code{max_expand} is the maximum number of targets to generate in each
\code{map()}, \code{split()}, or \code{cross()} transform.
Useful if you have a massive plan and you want to
test and visualize a strategic subset of targets
before scaling up.
Note: the \code{max_expand} argument of \code{drake_plan()} and
\code{transform_plan()} is for static branching only.
The dynamic branching \code{max_expand}
is an argument of \code{make()} and \code{drake_config()}.}
}
\value{
A data frame of targets, commands, and optional
custom columns.
}
\description{
A \code{drake} plan is a data frame with columns
\code{"target"} and \code{"command"}. Each target is an R object
produced in your workflow, and each command is the
R code to produce it.
}
\details{
Besides \code{"target"} and \code{"command"}, \code{\link[=drake_plan]{drake_plan()}}
understands a special set of optional columns. For details, visit
\url{https://books.ropensci.org/drake/plans.html#special-custom-columns-in-your-plan} # nolint
}
\section{Columns}{

\code{\link[=drake_plan]{drake_plan()}} creates a special data frame. At minimum, that data frame
must have columns \code{target} and \code{command} with the target names and the
R code chunks to build them, respectively.

You can add custom columns yourself, either with \code{target()} (e.g.
\code{drake_plan(y = target(f(x), transform = map(c(1, 2)), format = "fst"))})
or by appending columns post-hoc (e.g. \code{plan$col <- vals}).

Some of these custom columns are special. They are optional,
but \code{drake} looks for them at various points in the workflow.
\itemize{
\item \code{transform}: a call to \code{\link[=map]{map()}}, \code{\link[=split]{split()}}, \code{\link[=cross]{cross()}}, or
\code{\link[=combine]{combine()}} to create and manipulate large collections of targets.
Details: (\url{https://books.ropensci.org/drake/plans.html#large-plans}). # nolint
\item \code{format}: set a storage format to save big targets more efficiently.
See the "Formats" section of this help file for more details.
\item \code{trigger}: rule to decide whether a target needs to run.
It is recommended that you define this one with \code{target()}.
Details: \url{https://books.ropensci.org/drake/triggers.html}.
\item \code{hpc}: logical values (\code{TRUE}/\code{FALSE}/\code{NA}) whether to send each target
to parallel workers.
Visit \url{https://books.ropensci.org/drake/hpc.html#selectivity}
to learn more.
\item \code{resources}: target-specific lists of resources for a computing cluster.
See
\url{https://books.ropensci.org/drake/hpc.html#advanced-options}
for details.
\item \code{caching}: overrides the \code{caching} argument of \code{\link[=make]{make()}} for each target
individually. Possible values:
\itemize{
\item "main": tell the main process to store the target in the cache.
\item "worker": tell the HPC worker to store the target in the cache.
\item NA: default to the \code{caching} argument of \code{\link[=make]{make()}}.
}
\item \code{elapsed} and \code{cpu}: number of seconds to wait for the target to build
before timing out (\code{elapsed} for elapsed time and \code{cpu} for CPU time).
\item \code{retries}: number of times to retry building a target
in the event of an error.
\item \code{seed}: an optional pseudo-random number generator (RNG)
seed for each target. \code{drake} usually comes up with its own
unique reproducible target-specific seeds using the global seed
(the \code{seed} argument to \code{\link[=make]{make()}} and \code{\link[=drake_config]{drake_config()}})
and the target names, but you can overwrite these automatic seeds.
\code{NA} entries default back to \code{drake}'s automatic seeds.
\item \code{max_expand}: for dynamic branching only. Same as the \code{max_expand}
argument of \code{\link[=make]{make()}}, but on a target-by-target basis.
Limits the number of sub-targets created for a given target.
}
}

\section{Formats}{

Specialized target formats increase efficiency and flexibility.
Some allow you to save specialized objects like \code{keras} models,
while others increase the speed while conserving storage and memory.
You can declare target-specific formats in the plan
(e.g. \code{drake_plan(x = target(big_data_frame, format = "fst"))})
or supply a global default \code{format} for all targets in \code{make()}.
Either way, most formats have specialized installation requirements
(e.g. R packages) that are not installed with \code{drake} by default.
You will need to install them separately yourself.
Available formats:
\itemize{
\item \code{"file"}: Dynamic files. To use this format, simply create
local files and directories yourself and then return
a character vector of paths as the target's value.
Then, \code{drake} will watch for changes to those files in
subsequent calls to \code{make()}. This is a more flexible
alternative to \code{file_in()} and \code{file_out()}, and it is
compatible with dynamic branching.
See \url{https://github.com/ropensci/drake/pull/1178} for an example.
\item \code{"fst"}: save big data frames fast. Requires the \code{fst} package.
Note: this format strips non-data-frame attributes such as the
\item \code{"fst_tbl"}: Like \code{"fst"}, but for \code{tibble} objects.
Requires the \code{fst} and \code{tibble} packages.
Strips away non-data-frame non-tibble attributes.
\item \code{"fst_dt"}: Like \code{"fst"} format, but for \code{data.table} objects.
Requires the \code{fst} and \code{data.table} packages.
Strips away non-data-frame non-data-table attributes.
\item \code{"diskframe"}:
Stores \code{disk.frame} objects, which could potentially be
larger than memory. Requires the \code{fst} and \code{disk.frame} packages.
Coerces objects to \code{disk.frame}s.
Note: \code{disk.frame} objects get moved to the \code{drake} cache
(a subfolder of \verb{.drake/} for most workflows).
To ensure this data transfer is fast, it is best to
save your \code{disk.frame} objects to the same physical storage
drive as the \code{drake} cache,
\code{as.disk.frame(your_dataset, outdir = drake_tempfile())}.
\item \code{"keras"}: save Keras models as HDF5 files.
Requires the \code{keras} package.
\item \code{"qs"}: save any R object that can be properly serialized
with the \code{qs} package. Requires the \code{qs} package.
Uses \code{qsave()} and \code{qread()}.
Uses the default settings in \code{qs} version 0.20.2.
\item \code{"rds"}: save any R object that can be properly serialized.
Requires R version >= 3.5.0 due to ALTREP.
Note: the \code{"rds"} format uses gzip compression, which is slow.
\code{"qs"} is a superior format.
}
}

\section{Keywords}{

\code{\link[=drake_plan]{drake_plan()}} understands special keyword functions for your commands.
With the exception of \code{\link[=target]{target()}}, each one is a proper function
with its own help file.
\itemize{
\item \code{\link[=target]{target()}}: give the target more than just a command.
Using \code{\link[=target]{target()}}, you can apply a transformation
(examples: \url{https://books.ropensci.org/drake/plans.html#large-plans}), # nolint
supply a trigger (\url{https://books.ropensci.org/drake/triggers.html}), # nolint
or set any number of custom columns.
\item \code{\link[=file_in]{file_in()}}: declare an input file dependency.
\item \code{\link[=file_out]{file_out()}}: declare an output file to be produced
when the target is built.
\item \code{\link[=knitr_in]{knitr_in()}}: declare a \code{knitr} file dependency such as an
R Markdown (\verb{*.Rmd}) or R LaTeX (\verb{*.Rnw}) file.
\item \code{\link[=ignore]{ignore()}}: force \code{drake} to entirely ignore a piece of code:
do not track it for changes and do not analyze it for dependencies.
\item \code{\link[=no_deps]{no_deps()}}: tell \code{drake} to not track the dependencies
of a piece of code. \code{drake} still tracks the code itself for changes.
\item \code{\link[=id_chr]{id_chr()}}: Get the name of the current target.
\item \code{\link[=drake_envir]{drake_envir()}}: get the environment where drake builds targets.
Intended for advanced custom memory management.
}
}

\section{Transformations}{

\code{drake} has special syntax for generating large plans.
Your code will look something like
\verb{drake_plan(y = target(f(x), transform = map(x = c(1, 2, 3)))}
You can read about this interface at
\url{https://books.ropensci.org/drake/plans.html#large-plans}. # nolint
}

\section{Static branching}{

In static branching, you define batches of targets
based on information you know in advance.
Overall usage looks like
\verb{drake_plan(<x> = target(<...>, transform = <call>)},
where
\itemize{
\item \verb{<x>} is the name of the target or group of targets.
\item \verb{<...>} is optional arguments to \code{\link[=target]{target()}}.
\item \verb{<call>} is a call to one of the transformation functions.
}

Transformation function usage:
\itemize{
\item \code{map(..., .data, .names, .id, .tag_in, .tag_out)}
\item \code{split(..., slices, margin = 1L, drop = FALSE, .names, .tag_in, .tag_out)} # nolint
\item \code{cross(..., .data, .names, .id, .tag_in, .tag_out)}
\item \code{combine(..., .by, .names, .id, .tag_in, .tag_out)}
}
}

\section{Dynamic branching}{

\itemize{
\item \code{map(..., .trace)}
\item \code{cross(..., .trace)}
\item \code{group(..., .by, .trace)}
}

\code{map()} and \code{cross()} create dynamic sub-targets from the variables
supplied to the dots. As with static branching, the variables
supplied to \code{map()} must all have equal length.
\code{group(f(data), .by = x)} makes new dynamic
sub-targets from \code{data}. Here, \code{data} can be either static or dynamic.
If \code{data} is dynamic, \code{group()} aggregates existing sub-targets.
If \code{data} is static, \code{group()} splits \code{data} into multiple
subsets based on the groupings from \code{.by}.

Differences from static branching:
\itemize{
\item \code{...} must contain \emph{unnamed} symbols with no values supplied,
and they must be the names of targets.
\item Arguments \code{.id}, \code{.tag_in}, and \code{.tag_out} no longer apply.
}
}

\examples{
\dontrun{
isolate_example("contain side effects", {
# For more examples, visit
# https://books.ropensci.org/drake/plans.html.

# Create drake plans:
mtcars_plan <- drake_plan(
  write.csv(mtcars[, c("mpg", "cyl")], file_out("mtcars.csv")),
  value = read.csv(file_in("mtcars.csv"))
)
if (requireNamespace("visNetwork", quietly = TRUE)) {
  plot(mtcars_plan) # fast simplified call to vis_drake_graph()
}
mtcars_plan
make(mtcars_plan) # Makes `mtcars.csv` and then `value`
head(readd(value))
# You can use knitr inputs too. See the top command below.

load_mtcars_example()
head(my_plan)
if (requireNamespace("knitr", quietly = TRUE)) {
  plot(my_plan)
}
# The `knitr_in("report.Rmd")` tells `drake` to dive into the active
# code chunks to find dependencies.
# There, `drake` sees that `small`, `large`, and `coef_regression2_small`
# are loaded in with calls to `loadd()` and `readd()`.
deps_code("report.Rmd")

# Formats are great for big data: https://github.com/ropensci/drake/pull/977
# Below, each target is 1.6 GB in memory.
# Run make() on this plan to see how much faster fst is!
n <- 1e8
plan <- drake_plan(
  data_fst = target(
    data.frame(x = runif(n), y = runif(n)),
    format = "fst"
  ),
  data_old = data.frame(x = runif(n), y = runif(n))
)

# Use transformations to generate large plans.
# Read more at
# <https://books.ropensci.org/drake/plans.html#create-large-plans-the-easy-way>. # nolint
drake_plan(
  data = target(
    simulate(nrows),
    transform = map(nrows = c(48, 64)),
    custom_column = 123
  ),
  reg = target(
    reg_fun(data),
   transform = cross(reg_fun = c(reg1, reg2), data)
  ),
  summ = target(
    sum_fun(data, reg),
   transform = cross(sum_fun = c(coef, residuals), reg)
  ),
  winners = target(
    min(summ),
    transform = combine(summ, .by = c(data, sum_fun))
  )
)

# Split data among multiple targets.
drake_plan(
  large_data = get_data(),
  slice_analysis = target(
    analyze(large_data),
    transform = split(large_data, slices = 4)
  ),
  results = target(
    rbind(slice_analysis),
    transform = combine(slice_analysis)
  )
)

# Set trace = TRUE to show what happened during the transformation process.
drake_plan(
  data = target(
    simulate(nrows),
    transform = map(nrows = c(48, 64)),
    custom_column = 123
  ),
  reg = target(
    reg_fun(data),
   transform = cross(reg_fun = c(reg1, reg2), data)
  ),
  summ = target(
    sum_fun(data, reg),
   transform = cross(sum_fun = c(coef, residuals), reg)
  ),
  winners = target(
    min(summ),
    transform = combine(summ, .by = c(data, sum_fun))
  ),
  trace = TRUE
)

# You can create your own custom columns too.
# See ?triggers for more on triggers.
drake_plan(
  website_data = target(
    command = download_data("www.your_url.com"),
    trigger = "always",
    custom_column = 5
  ),
  analysis = analyze(website_data)
)

# Tidy evaluation can help generate super large plans.
sms <- rlang::syms(letters) # To sub in character args, skip this.
drake_plan(x = target(f(char), transform = map(char = !!sms)))

# Dynamic branching
# Get the mean mpg for each cyl in the mtcars dataset.
plan <- drake_plan(
  raw = mtcars,
  group_index = raw$cyl,
  munged = target(raw[, c("mpg", "cyl")], dynamic = map(raw)),
  mean_mpg_by_cyl = target(
    data.frame(mpg = mean(munged$mpg), cyl = munged$cyl[1]),
    dynamic = group(munged, .by = group_index)
  )
)
make(plan)
readd(mean_mpg_by_cyl)
})
}
}
\seealso{
make, drake_config, transform_plan, map, split, cross, combine
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{as_file}
\alias{as_file}
\title{as_file \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
as_file(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{drake_running}
\alias{drake_running}
\title{List running targets.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_running(cache = drake::drake_cache(path = path), path = NULL)
}
\arguments{
\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}
}
\value{
A character vector of target names.
}
\description{
List the targets that either
\enumerate{
\item Are currently being built during a call to \code{\link[=make]{make()}}, or
\item Were in progress when \code{\link[=make]{make()}} was interrupted.
}
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the targets.
drake_running() # Everything should be done.
# nolint start
# Run make() in one R session...
# slow_plan <- drake_plan(x = Sys.sleep(2))
# make(slow_plan)
# and see the progress in another session.
# drake_running()
# nolint end
}
})
}
}
\seealso{
\code{\link[=drake_done]{drake_done()}}, \code{\link[=drake_failed]{drake_failed()}}, \code{\link[=drake_cancelled]{drake_cancelled()}},
\code{\link[=drake_progress]{drake_progress()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_example.R
\name{drake_examples}
\alias{drake_examples}
\title{List the names of all the drake examples.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_examples(quiet = TRUE)
}
\arguments{
\item{quiet}{Logical, passed to \code{downloader::download()}
and thus \code{utils::download.file()}. Whether
to download quietly or print progress.}
}
\value{
Names of all the drake examples.
}
\description{
You can find the code files of the examples at
\url{https://github.com/wlandau/drake-examples}.
The \code{drake_examples()} function downloads the list of examples
from \url{https://wlandau.github.io/drake-examples/examples.md},
so you need an internet connection.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (requireNamespace("downloader")) {
drake_examples() # List all the drake examples.
# Sets up the example from load_mtcars_example()
drake_example("mtcars")
# Sets up the SLURM example.
drake_example("slurm")
}
})
}
}
\seealso{
\code{\link[=drake_example]{drake_example()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outdated.R
\name{recoverable}
\alias{recoverable}
\title{List the most upstream \emph{recoverable} outdated targets.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
recoverable(..., make_imports = TRUE, do_prework = TRUE, config = NULL)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets} and \code{envir}.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to save some time and risk obsolete output.}

\item{do_prework}{Whether to do the \code{prework}
normally supplied to \code{\link[=make]{make()}}.}

\item{config}{Deprecated (2019-12-21).
A configured workflow from \code{\link[=drake_config]{drake_config()}}.}
}
\value{
Character vector of the names of recoverable targets.
}
\description{
Only shows the most upstream updated targets.
Whether downstream targets are recoverable depends on
the eventual values of the upstream targets in the next \code{\link[=make]{make()}}.
}
\section{Recovery}{

\code{make(recover = TRUE, recoverable = TRUE)}
powers automated data recovery.
The default of \code{recover} is \code{FALSE} because
targets recovered from the distant past may have been generated
with earlier versions of R and earlier package environments
that no longer exist.

How it works: if \code{recover} is \code{TRUE},
\code{drake} tries to salvage old target values from the cache
instead of running commands from the plan.
A target is recoverable if
\enumerate{
\item There is an old value somewhere in the cache that
shares the command, dependencies, etc.
of the target about to be built.
\item The old value was generated with \code{make(recoverable = TRUE)}.
}

If both conditions are met, \code{drake} will
\enumerate{
\item Assign the most recently-generated admissible data to the target, and
\item skip the target's command.
}
}

\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan)
clean()
outdated(my_plan) # Which targets are outdated?
recoverable(my_plan) # Which of these are recoverable and upstream?
# The report still builds because clean() removes report.md,
# but make() recovers the rest.
make(my_plan, recover = TRUE)
outdated(my_plan)
# When was the *recovered* small data actually built (first stored)?
# (Was I using a different version of R back then?)
diagnose(small)$date
# If you set the same seed as before, you can even
# rename targets without having to build them again.
# For an example, see
# the "Reproducible data recovery and renaming" section of
# https://github.com/ropensci/drake/blob/main/README.md.
}
})
}
}
\seealso{
\code{\link[=r_recoverable]{r_recoverable()}}, \code{\link[=r_outdated]{r_outdated()}}, \code{\link[=drake_config]{drake_config()}}, \code{\link[=missed]{missed()}},
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{diagnose}
\alias{diagnose}
\title{Get diagnostic metadata on a target.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
diagnose(
  target = NULL,
  character_only = FALSE,
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  verbose = 1L
)
}
\arguments{
\item{target}{Name of the target of the error to get.
Can be a symbol if \code{character_only} is \code{FALSE},
must be a character if \code{character_only} is \code{TRUE}.}

\item{character_only}{Logical, whether \code{target} should be treated
as a character or a symbol.
Just like \code{character.only} in \code{\link[=library]{library()}}.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\value{
Either a character vector of target names or an object
of class \code{"error"}.
}
\description{
Diagnostics include errors, warnings,
messages, runtimes, and other context/metadata from when a
target was built or an import was processed.
If your target's last build succeeded,
then \code{diagnose(your_target)} has the most current information
from that build.
But if your target failed, then only
\code{diagnose(your_target)$error},
\code{diagnose(your_target)$warnings},
and \code{diagnose(your_target)$messages} correspond to the failure,
and all the other metadata correspond to the last build that completed
without an error.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
diagnose() # List all the targets with recorded error logs.
# Define a function doomed to failure.
f <- function() {
  stop("unusual error")
}
# Create a workflow plan doomed to failure.
bad_plan <- drake_plan(my_target = f())
# Running the project should generate an error
# when trying to build 'my_target'.
try(make(bad_plan), silent = FALSE)
drake_failed() # List the failed targets from the last make() (my_target).
# List targets that failed at one point or another
# over the course of the project (my_target).
# drake keeps all the error logs.
diagnose()
# Get the error log, an object of class "error".
error <- diagnose(my_target)$error # See also warnings and messages.
str(error) # See what's inside the error log.
error$calls # View the traceback. (See the rlang::trace_back() function).
})
}
}
\seealso{
\code{\link[=drake_failed]{drake_failed()}}, \code{\link[=drake_progress]{drake_progress()}},
\code{\link[=readd]{readd()}}, \code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{session}
\alias{session}
\title{session \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
session(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{read_drake_plan}
\alias{read_drake_plan}
\title{Read the plan from the cache
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
read_drake_plan(path = getwd(), search = TRUE, cache = NULL, verbose = 1L)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\description{
drake no longer stores the config object,
the plan, etc. in the cache during \code{make()}. This change
improves speed.
}
\details{
2019-01-06
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{drake_get_session_info}
\alias{drake_get_session_info}
\title{Session info of the last call to \code{\link[=make]{make()}}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_get_session_info(
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  verbose = 1L
)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\value{
\code{\link[=sessionInfo]{sessionInfo()}} of the last
call to \code{\link[=make]{make()}}
}
\description{
By default, session info is saved
during \code{\link[=make]{make()}} to ensure reproducibility.
Your loaded packages and their versions are recorded, for example.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the targets.
drake_get_session_info() # Get the cached sessionInfo() of the last make().
}
})
}
}
\seealso{
\code{\link[=diagnose]{diagnose()}}, \code{\link[=cached]{cached()}},
\code{\link[=readd]{readd()}}, \code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{default_Makefile_command}
\alias{default_Makefile_command}
\title{Default Makefile command
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
default_Makefile_command()
}
\value{
A character scalar
}
\description{
2019-01-03
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{find_project}
\alias{find_project}
\title{Search up the file system
for the nearest root path of a drake project.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
find_project(path = getwd())
}
\arguments{
\item{path}{Starting path for search back for the project.
Should be a subdirectory of the drake project.}
}
\value{
File path of the nearest drake project or \code{NULL}
if no drake project is found.
}
\description{
Deprecated on 2019-01-08.
}
\details{
Only works if the cache is a file system
in a folder named \code{.drake} (default).
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{gather_by}
\alias{gather_by}
\title{Gather multiple groupings of targets
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
gather_by(
  plan,
  ...,
  prefix = "target",
  gather = "list",
  append = TRUE,
  filter = NULL,
  sep = "_"
)
}
\arguments{
\item{plan}{Workflow plan data frame of prespecified targets.}

\item{...}{Symbols, columns of \code{plan} to define target groupings.
A \code{gather_plan()} call is applied for each grouping.
Groupings with all \code{NA}s in the selector variables are ignored.}

\item{prefix}{Character, prefix for naming the new targets.
Suffixes are generated from the values of the columns
specified in \code{...}.}

\item{gather}{Function used to gather the targets. Should be
one of \code{list(...)}, \code{c(...)}, \code{rbind(...)}, or similar.}

\item{append}{Logical. If \code{TRUE}, the output will include the
original rows in the \code{plan} argument.
If \code{FALSE}, the output will only include the new
targets and commands.}

\item{filter}{An expression like you would pass to \code{dplyr::filter()}.
The rows for which \code{filter} evaluates to \code{TRUE} will be gathered,
and the rest will be excluded from gathering.
Why not just call \code{dplyr::filter()} before \code{gather_by()}?
Because \code{gather_by(append = TRUE, filter = my_column == "my_value")}
gathers on some targets while including all the original targets
in the output. See the examples for a demonstration.}

\item{sep}{Character scalar, delimiter for creating the names
of new targets.}
}
\value{
A workflow plan data frame.
}
\description{
Deprecated on 2019-05-16. Use \code{\link[=drake_plan]{drake_plan()}}
transformations instead. See
\url{https://books.ropensci.org/drake/plans.html#large-plans}
for the details.
}
\details{
Perform several calls to \code{gather_plan()}
based on groupings from columns in the plan,
and then row-bind the new targets to the plan.
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_build.R
\name{debug_and_run}
\alias{debug_and_run}
\title{Run a function in debug mode.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
debug_and_run(f)
}
\arguments{
\item{f}{A function.}
}
\value{
The return value of \code{f}.
}
\description{
Internal function for \code{\link[=drake_debug]{drake_debug()}}. Not for general use.
}
\seealso{
\code{\link[=drake_debug]{drake_debug()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decorate_storr.R
\name{file_store}
\alias{file_store}
\title{Show a file's encoded representation in the cache
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
file_store(x)
}
\arguments{
\item{x}{Character string to be turned into a filename
understandable by drake (i.e., a string with literal
single quotes on both ends).}
}
\value{
A single-quoted character string: i.e., a filename
understandable by drake.
}
\description{
This function simply wraps literal double quotes around
the argument \code{x} so \code{drake} knows it is the name of a file.
Use when you are calling functions like \code{deps_code()}: for example,
\code{deps_code(file_store("report.md"))}. See the examples for details.
Internally, \code{drake} wraps the names of file targets/imports
inside literal double quotes to avoid confusion between
files and generic R objects.
}
\examples{
# Wraps the string in single quotes.
file_store("my_file.rds") # "'my_file.rds'"
\dontrun{
isolate_example("contain side effects", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the workflow to build the targets
list.files() # Should include input "report.Rmd" and output "report.md".
head(readd(small)) # You can use symbols for ordinary objects.
# But if you want to read cached info on files, use `file_store()`.
readd(file_store("report.md"), character_only = TRUE) # File fingerprint.
deps_code(file_store("report.Rmd"))
config <- drake_config(my_plan)
deps_profile(
  file_store("report.Rmd"),
  plan = my_plan,
  character_only = TRUE
)
}
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{drake_quotes}
\alias{drake_quotes}
\title{Put quotes around each element of a character vector.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
drake_quotes(x = NULL, single = FALSE)
}
\arguments{
\item{x}{Character vector or object to be coerced to character.}

\item{single}{Add single quotes if \code{TRUE}
and double quotes otherwise.}
}
\value{
Character vector with quotes around it.
}
\description{
Deprecated on 2019-01-01
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{read_plan}
\alias{read_plan}
\title{read_plan \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
read_plan(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{reduce_plan}
\alias{reduce_plan}
\title{Write commands to reduce several targets down to one.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
reduce_plan(
  plan = NULL,
  target = "target",
  begin = "",
  op = " + ",
  end = "",
  pairwise = TRUE,
  append = FALSE,
  sep = "_"
)
}
\arguments{
\item{plan}{Workflow plan data frame of prespecified targets.}

\item{target}{Name of the new reduced target.}

\item{begin}{Character, code to place at the beginning
of each step in the reduction.}

\item{op}{Binary operator to apply in the reduction}

\item{end}{Character, code to place at the end
of each step in the reduction.}

\item{pairwise}{Logical, whether to create multiple
new targets, one for each pair/step in the reduction (\code{TRUE}),
or to do the reduction all in one command.}

\item{append}{Logical. If \code{TRUE}, the output will include the
original rows in the \code{plan} argument.
If \code{FALSE}, the output will only include the new
targets and commands.}

\item{sep}{Character scalar, delimiter for creating new target names.}
}
\value{
A workflow plan data frame that aggregates multiple
prespecified targets into one additional target downstream.
}
\description{
Deprecated on 2019-05-16. Use \code{\link[=drake_plan]{drake_plan()}}
transformations instead. See
\url{https://books.ropensci.org/drake/plans.html#large-plans}
for the details.
}
\details{
Creates a new workflow plan data frame with the
commands to do a reduction (i.e. to repeatedly apply a binary
operator to pairs of targets to produce one target).
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_graph_info.R
\name{legend_nodes}
\alias{legend_nodes}
\title{Create the nodes data frame used in the legend
of the graph visualizations.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#soft-deprecated}{\figure{lifecycle-soft-deprecated.svg}{options: alt='[Soft-deprecated]'}}}{\strong{[Soft-deprecated]}}}
\usage{
legend_nodes(font_size = 20)
}
\arguments{
\item{font_size}{Font size of the node label text.}
}
\value{
A data frame of legend nodes for the graph visualizations.
}
\description{
Output a \code{visNetwork}-friendly
data frame of nodes. It tells you what
the colors and shapes mean
in the graph visualizations.
}
\examples{
\dontrun{
# Show the legend nodes used in graph visualizations.
# For example, you may want to inspect the color palette more closely.
if (requireNamespace("visNetwork", quietly = TRUE)) {
# visNetwork::visNetwork(nodes = legend_nodes()) # nolint
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{backend}
\alias{backend}
\title{backend \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
backend(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_script.R
\name{drake_script}
\alias{drake_script}
\title{Write an example \verb{_drake.R} script to the current working directory.}
\usage{
drake_script(code = NULL)
}
\arguments{
\item{code}{R code to put in \verb{_drake.R} in the current working directory.
If \code{NULL}, an example script is written.}
}
\value{
Nothing.
}
\description{
A \verb{_drake.R} file is required for \code{\link[=r_make]{r_make()}} and friends.
See the \code{\link[=r_make]{r_make()}} help file for details.
}
\examples{
\dontrun{
isolate_example("contain side-effects", {
drake_script({
  library(drake)
  plan <- drake_plan(x = 1)
  drake_config(plan, lock_cache = FALSE)
})
cat(readLines("_drake.R"), sep = "\n")
r_make()
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_runtime.R
\name{predict_workers}
\alias{predict_workers}
\title{Predict the load balancing of the next call to \code{make()}
for non-staged parallel backends.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
predict_workers(
  ...,
  targets_predict = NULL,
  from_scratch = FALSE,
  targets_only = NULL,
  jobs_predict = 1L,
  known_times = numeric(0),
  default_time = 0,
  warn = TRUE,
  config = NULL
)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{targets_predict}{Character vector, names of targets
to include in the total runtime and worker predictions.}

\item{from_scratch}{Logical, whether to predict a
\code{\link[=make]{make()}} build from scratch or to
take into account the fact that some targets may be
already up to date and therefore skipped.}

\item{targets_only}{Deprecated.}

\item{jobs_predict}{The \code{jobs} argument of your next planned
\code{make()}.}

\item{known_times}{A named numeric vector with targets/imports
as names and values as hypothetical runtimes in seconds.
Use this argument to overwrite any of the existing build times
or the \code{default_time}.}

\item{default_time}{Number of seconds to assume for any
target or import with no recorded runtime (from \code{\link[=build_times]{build_times()}})
or anything in \code{known_times}.}

\item{warn}{Logical, whether to warn the user about
any targets with no available runtime, either in
\code{known_times} or \code{\link[=build_times]{build_times()}}. The times for these
targets default to \code{default_time}.}

\item{config}{Deprecated.}
}
\value{
A data frame showing one likely arrangement
of targets assigned to parallel workers.
}
\description{
Take the past recorded runtimes times from
\code{\link[=build_times]{build_times()}} and use them to predict how the targets
will be distributed among the available workers in the
next \code{\link[=make]{make()}}.
Predictions only include the time it takes to run the targets,
not overhead/preprocessing from \code{drake} itself.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the targets.
known_times <- rep(7200, nrow(my_plan))
names(known_times) <- my_plan$target
known_times
# Predict the runtime
if (requireNamespace("lubridate", quietly = TRUE)) {
predict_runtime(
  my_plan,
  jobs_predict = 7L,
  from_scratch = TRUE,
  known_times = known_times
)
predict_runtime(
  my_plan,
  jobs_predict = 8L,
  from_scratch = TRUE,
  known_times = known_times
)
balance <- predict_workers(
  my_plan,
  jobs_predict = 7L,
  from_scratch = TRUE,
  known_times = known_times
)
balance
}
}
})
}
}
\seealso{
\code{\link[=predict_runtime]{predict_runtime()}}, \code{\link[=build_times]{build_times()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{default_system2_args}
\alias{default_system2_args}
\title{default_system2_args \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
default_system2_args(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{no_deps}
\alias{no_deps}
\title{Suppress dependency detection.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
no_deps(x = NULL)
}
\arguments{
\item{x}{Code for which dependency detection is suppressed.}
}
\value{
The argument.
}
\description{
Tell \code{drake} to not search for dependencies in a chunk of code.
}
\details{
\code{no_deps()} is similar to \code{\link[=ignore]{ignore()}}, but it still lets \code{drake}
track meaningful changes to the code itself.
}
\section{Keywords}{

\code{\link[=drake_plan]{drake_plan()}} understands special keyword functions for your commands.
With the exception of \code{\link[=target]{target()}}, each one is a proper function
with its own help file.
\itemize{
\item \code{\link[=target]{target()}}: give the target more than just a command.
Using \code{\link[=target]{target()}}, you can apply a transformation
(examples: \url{https://books.ropensci.org/drake/plans.html#large-plans}), # nolint
supply a trigger (\url{https://books.ropensci.org/drake/triggers.html}), # nolint
or set any number of custom columns.
\item \code{\link[=file_in]{file_in()}}: declare an input file dependency.
\item \code{\link[=file_out]{file_out()}}: declare an output file to be produced
when the target is built.
\item \code{\link[=knitr_in]{knitr_in()}}: declare a \code{knitr} file dependency such as an
R Markdown (\verb{*.Rmd}) or R LaTeX (\verb{*.Rnw}) file.
\item \code{\link[=ignore]{ignore()}}: force \code{drake} to entirely ignore a piece of code:
do not track it for changes and do not analyze it for dependencies.
\item \code{\link[=no_deps]{no_deps()}}: tell \code{drake} to not track the dependencies
of a piece of code. \code{drake} still tracks the code itself for changes.
\item \code{\link[=id_chr]{id_chr()}}: Get the name of the current target.
\item \code{\link[=drake_envir]{drake_envir()}}: get the environment where drake builds targets.
Intended for advanced custom memory management.
}
}

\examples{
\dontrun{
isolate_example("Contain side effects", {
# Normally, `drake` reacts to changes in dependencies.
x <- 4
make(plan = drake_plan(y = sqrt(x)))
x <- 5
make(plan = drake_plan(y = sqrt(x)))
make(plan = drake_plan(y = sqrt(4) + x))
# But not with no_deps().
make(plan = drake_plan(y = sqrt(4) + no_deps(x))) # Builds y.
x <- 6
make(plan = drake_plan(y = sqrt(4) + no_deps(x))) # Skips y.
# However, `drake` *does* react to changes
# to the *literal code* inside `no_deps()`.
make(plan = drake_plan(y = sqrt(4) + ignore(x + 1))) # Builds y.

# Like ignore(), no_deps() works with functions and multiline code chunks.
z <- 1
f <- function(x) {
  no_deps({
    x <- z + 1
    x <- x + 2
  })
  x
}
make(plan = drake_plan(y = f(2)))
readd(y)
z <- 2 # Changed dependency is not tracked.
make(plan = drake_plan(y = f(2)))
readd(y)
})
}
}
\seealso{
\code{\link[=file_in]{file_in()}}, \code{\link[=file_out]{file_out()}}, \code{\link[=knitr_in]{knitr_in()}}, \code{\link[=no_deps]{no_deps()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{workplan}
\alias{workplan}
\title{workplan \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
workplan(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{cancel_if}
\alias{cancel_if}
\title{Cancel a target mid-build under some condition
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}})`}
\usage{
cancel_if(condition, allow_missing = TRUE)
}
\arguments{
\item{condition}{Logical, whether to cancel the target.}

\item{allow_missing}{Logical. If \code{FALSE}, \code{drake} will not cancel
the target if it is missing from the cache (or if you removed the
key with \code{clean()}).}
}
\value{
Nothing.
}
\description{
Cancel a target mid-build if some logical condition is met.
Upon cancellation, \code{drake} halts the current target and moves to the
next one. The target's previous value and metadata, if they exist,
remain in the cache.
}
\examples{
\dontrun{
isolate_example("cancel_if()", {
f <- function(x) {
  cancel_if(x > 1)
  Sys.sleep(2) # Does not run if x > 1.
}
g <- function(x) f(x)
plan <- drake_plan(y = g(2))
make(plan)
# Does not exist.
# readd(y)
})
}
}
\seealso{
cancel
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_drake_triggers.R
\name{new_drake_triggers}
\alias{new_drake_triggers}
\title{\code{drake_triggers} constructor}
\usage{
new_drake_triggers(
  command = TRUE,
  depend = TRUE,
  file = TRUE,
  seed = TRUE,
  format = TRUE,
  condition = FALSE,
  change = NULL,
  mode = "whitelist"
)
}
\arguments{
\item{command}{Logical, command trigger.}

\item{depend}{Logical, depend trigger.}

\item{file}{Logical, file trigger.}

\item{seed}{Logical, seed trigger.}

\item{format}{Logical, format trigger.}

\item{condition}{Language object or object coercible to logical,
condition trigger.}

\item{change}{Language object or literal value, change trigger.}

\item{mode}{Character, mode of condition trigger.}
}
\value{
A \code{drake_triggers} object.
}
\description{
List of class \code{drake_triggers}.
}
\examples{
if (FALSE) { # stronger than roxygen dontrun
new_drake_triggers()
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{cleaned_namespaces}
\alias{cleaned_namespaces}
\title{Auxiliary storr namespaces
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
cleaned_namespaces(default = storr::storr_environment()$default_namespace)
}
\arguments{
\item{default}{Name of the default \code{storr} namespace.}
}
\value{
A character vector of \code{storr} namespaces
that are cleaned during \code{\link[=clean]{clean()}}.
}
\description{
2019-02-13
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decorate_storr.R
\name{drake_tempfile}
\alias{drake_tempfile}
\title{drake tempfile
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_tempfile(path = NULL, cache = drake::drake_cache(path = path))
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}
}
\description{
Create the path to a temporary file inside drake's cache.
}
\details{
This function is just like the \code{tempfile()} function in base R
except that the path points to a special location inside \code{drake}'s cache.
This ensures that if the file needs to be copied to
persistent storage in the cache, \code{drake} does not need to copy across
physical storage media. Example: the \code{"diskframe"} format. See the
"Formats" and "Columns" sections of the \code{\link[=drake_plan]{drake_plan()}} help file.
Unless you supply the cache or the path to the cache
(see \code{\link[=drake_cache]{drake_cache()}}) \code{drake} will assume the cache folder is named
\verb{.drake/} and it is located either in your working directory or an
ancestor of your working directory.
}
\examples{
cache <- new_cache(tempfile())
# No need to supply a cache if a .drake/ folder exists.
drake_tempfile(cache = cache)
drake_plan(
  x = target(
    as.disk.frame(large_data, outdir = drake_tempfile()),
    format = "diskframe"
  )
)
}
\seealso{
\code{\link[=drake_cache]{drake_cache()}}, \code{\link[=new_cache]{new_cache()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/backend_future.R
\name{future_build}
\alias{future_build}
\title{Task passed to individual futures in the \code{"future"} backend
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
future_build(target, meta, config, spec, config_tmp, protect)
}
\arguments{
\item{target}{Name of the target.}

\item{meta}{A list of metadata.}

\item{config}{A \code{\link[=drake_config]{drake_config()}} list.}

\item{config_tmp}{Internal, parts of \code{config} that the workers need.}

\item{protect}{Names of targets that still need their
dependencies available in memory.}
}
\value{
Either the target value or a list of build results.
}
\description{
For internal use only. Only exported to make available
to futures.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{id_chr}
\alias{id_chr}
\title{Name of the current target \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
id_chr()
}
\value{
The name of the current target.
}
\description{
\code{id_chr()} gives you the name of the current target
while \code{\link[=make]{make()}} is running. For static branching in \code{\link[=drake_plan]{drake_plan()}},
use the \code{.id_chr} symbol instead. See the examples for details.
}
\section{Keywords}{

\code{\link[=drake_plan]{drake_plan()}} understands special keyword functions for your commands.
With the exception of \code{\link[=target]{target()}}, each one is a proper function
with its own help file.
\itemize{
\item \code{\link[=target]{target()}}: give the target more than just a command.
Using \code{\link[=target]{target()}}, you can apply a transformation
(examples: \url{https://books.ropensci.org/drake/plans.html#large-plans}), # nolint
supply a trigger (\url{https://books.ropensci.org/drake/triggers.html}), # nolint
or set any number of custom columns.
\item \code{\link[=file_in]{file_in()}}: declare an input file dependency.
\item \code{\link[=file_out]{file_out()}}: declare an output file to be produced
when the target is built.
\item \code{\link[=knitr_in]{knitr_in()}}: declare a \code{knitr} file dependency such as an
R Markdown (\verb{*.Rmd}) or R LaTeX (\verb{*.Rnw}) file.
\item \code{\link[=ignore]{ignore()}}: force \code{drake} to entirely ignore a piece of code:
do not track it for changes and do not analyze it for dependencies.
\item \code{\link[=no_deps]{no_deps()}}: tell \code{drake} to not track the dependencies
of a piece of code. \code{drake} still tracks the code itself for changes.
\item \code{\link[=id_chr]{id_chr()}}: Get the name of the current target.
\item \code{\link[=drake_envir]{drake_envir()}}: get the environment where drake builds targets.
Intended for advanced custom memory management.
}
}

\examples{
try(id_chr()) # Do not use outside the plan.
\dontrun{
isolate_example("id_chr()", {
plan <- drake_plan(x = id_chr())
make(plan)
readd(x)
# Dynamic branching
plan <- drake_plan(
  x = seq_len(4),
  y = target(id_chr(), dynamic = map(x))
)
make(plan)
readd(y, subtargets = 1)
# Static branching
plan <- drake_plan(
  y = target(c(x, .id_chr), transform = map(x = !!seq_len(4)))
)
plan
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{check_plan}
\alias{check_plan}
\title{Check a workflow plan data frame for obvious errors.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
check_plan(
  plan = NULL,
  targets = NULL,
  envir = parent.frame(),
  cache = drake::get_cache(verbose = verbose),
  verbose = 1L,
  jobs = 1
)
}
\arguments{
\item{plan}{Workflow plan data frame, possibly from
\code{\link[=drake_plan]{drake_plan()}}.}

\item{targets}{Character vector of targets to make.}

\item{envir}{Environment containing user-defined functions.}

\item{cache}{Optional drake cache. See \code{\link[=new_cache]{new_cache()}}.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{jobs}{Number of jobs/workers for parallel processing.}
}
\value{
Invisibly return \code{plan}.
}
\description{
Deprecated on 2019-01-12.
}
\details{
Possible obvious errors include circular dependencies and
missing input files.
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{read_drake_seed}
\alias{read_drake_seed}
\title{Read the pseudo-random number generator seed of the project.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
read_drake_seed(path = NULL, search = NULL, cache = NULL, verbose = NULL)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\value{
An integer vector.
}
\description{
When a project is created with \code{\link[=make]{make()}}
or \code{\link[=drake_config]{drake_config()}}, the project's pseudo-random number generator
seed is cached. Then, unless the cache is destroyed,
the seeds of all the targets will deterministically depend on
this one central seed. That way, reproducibility is protected,
even under randomness.
}
\examples{
\dontrun{
isolate_example("contain side effects", {
cache <- storr::storr_environment() # Just for the examples.
my_plan <- drake_plan(
  target1 = sqrt(1234),
  target2 = sample.int(n = 12, size = 1) + target1
)
tmp <- sample.int(1) # Needed to get a .Random.seed, but not for drake.
digest::digest(.Random.seed) # Fingerprint of the current R session's seed.
make(my_plan, cache = cache) # Run the project, build the targets.
digest::digest(.Random.seed) # Your session's seed did not change.
# drake uses a hard-coded seed if you do not supply one.
read_drake_seed(cache = cache)
readd(target2, cache = cache) # Randomly-generated target data.
clean(target2, cache = cache) # Oops, I removed the data!
tmp <- sample.int(1) # Maybe the R session's seed also changed.
make(my_plan, cache = cache) # Rebuild target2.
# Same as before:
read_drake_seed(cache = cache)
readd(target2, cache = cache)
# You can also supply a seed.
# If your project already exists, it must agree with the project's
# preexisting seed (default: 0)
clean(target2, cache = cache)
make(my_plan, cache = cache, seed = 0)
read_drake_seed(cache = cache)
readd(target2, cache = cache)
# If you want to supply a different seed than 0,
# you need to destroy the cache and start over first.
clean(destroy = TRUE, cache = cache)
cache <- storr::storr_environment() # Just for the examples.
make(my_plan, cache = cache, seed = 1234)
read_drake_seed(cache = cache)
readd(target2, cache = cache)
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{analyses}
\alias{analyses}
\title{analyses \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
analyses(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.R
\name{make}
\alias{make}
\title{Run your project (build the outdated targets).
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
make(
  plan,
  targets = NULL,
  envir = parent.frame(),
  verbose = 1L,
  hook = NULL,
  cache = drake::drake_cache(),
  fetch_cache = NULL,
  parallelism = "loop",
  jobs = 1L,
  jobs_preprocess = 1L,
  packages = rev(.packages()),
  lib_loc = NULL,
  prework = character(0),
  prepend = NULL,
  command = NULL,
  args = NULL,
  recipe_command = NULL,
  log_progress = TRUE,
  skip_targets = FALSE,
  timeout = NULL,
  cpu = Inf,
  elapsed = Inf,
  retries = 0,
  force = FALSE,
  graph = NULL,
  trigger = drake::trigger(),
  skip_imports = FALSE,
  skip_safety_checks = FALSE,
  config = NULL,
  lazy_load = "eager",
  session_info = NULL,
  cache_log_file = NULL,
  seed = NULL,
  caching = "main",
  keep_going = FALSE,
  session = NULL,
  pruning_strategy = NULL,
  makefile_path = NULL,
  console_log_file = NULL,
  ensure_workers = NULL,
  garbage_collection = FALSE,
  template = list(),
  sleep = function(i) 0.01,
  hasty_build = NULL,
  memory_strategy = "speed",
  layout = NULL,
  spec = NULL,
  lock_envir = TRUE,
  history = TRUE,
  recover = FALSE,
  recoverable = TRUE,
  curl_handles = list(),
  max_expand = NULL,
  log_build_times = TRUE,
  format = NULL,
  lock_cache = TRUE,
  log_make = NULL,
  log_worker = FALSE
)
}
\arguments{
\item{plan}{Workflow plan data frame.
A workflow plan data frame is a data frame
with a \code{target} column and a \code{command} column.
(See the details in the \code{\link[=drake_plan]{drake_plan()}} help file
for descriptions of the optional columns.)
Targets are the objects that drake generates,
and commands are the pieces of R code that produce them.
You can create and track custom files along the way
(see \code{\link[=file_in]{file_in()}}, \code{\link[=file_out]{file_out()}}, and \code{\link[=knitr_in]{knitr_in()}}).
Use the function \code{\link[=drake_plan]{drake_plan()}} to generate workflow plan
data frames.}

\item{targets}{Character vector, names of targets to build.
Dependencies are built too. You may supply static and/or whole
dynamic targets, but no sub-targets.}

\item{envir}{Environment to use. Defaults to the current
workspace, so you should not need to worry about this
most of the time. A deep copy of \code{envir} is made,
so you don't need to worry about your workspace being modified
by \code{make}. The deep copy inherits from the global environment.
Wherever necessary, objects and functions are imported
from \code{envir} and the global environment and
then reproducibly tracked as dependencies.}

\item{verbose}{Integer, control printing to the console/terminal.
\itemize{
\item \code{0}: print nothing.
\item \code{1}: print target-by-target messages as \code{\link[=make]{make()}} progresses.
\item \code{2}: show a progress bar to track how many targets are
done so far.
}}

\item{hook}{Deprecated.}

\item{cache}{drake cache as created by \code{\link[=new_cache]{new_cache()}}.
See also \code{\link[=drake_cache]{drake_cache()}}.}

\item{fetch_cache}{Deprecated.}

\item{parallelism}{Character scalar, type of parallelism to use.
For detailed explanations, see the
\href{https://books.ropensci.org/drake/hpc.html}{high-performance computing chapter} # nolint
of the user manual.

You could also supply your own scheduler function
if you want to experiment or aggressively optimize.
The function should take a single \code{config} argument
(produced by \code{\link[=drake_config]{drake_config()}}). Existing examples
from \code{drake}'s internals are the \verb{backend_*()} functions:
\itemize{
\item \code{backend_loop()}
\item \code{backend_clustermq()}
\item \code{backend_future()}
However, this functionality is really a back door
and should not be used for production purposes unless you really
know what you are doing and you are willing to suffer setbacks
whenever \code{drake}'s unexported core functions are updated.
}}

\item{jobs}{Maximum number of parallel workers for processing the targets.
You can experiment with \code{\link[=predict_runtime]{predict_runtime()}}
to help decide on an appropriate number of jobs.
For details, visit
\url{https://books.ropensci.org/drake/time.html}.}

\item{jobs_preprocess}{Number of parallel jobs for processing the imports
and doing other preprocessing tasks.}

\item{packages}{Character vector packages to load, in the order
they should be loaded. Defaults to \code{rev(.packages())}, so you
should not usually need to set this manually. Just call
\code{\link[=library]{library()}} to load your packages before \code{make()}.
However, sometimes packages need to be strictly forced to load
in a certain order, especially if \code{parallelism} is
\code{"Makefile"}. To do this, do not use \code{\link[=library]{library()}}
or \code{\link[=require]{require()}} or \code{\link[=loadNamespace]{loadNamespace()}} or
\code{\link[=attachNamespace]{attachNamespace()}} to load any libraries beforehand.
Just list your packages in the \code{packages} argument in the order
you want them to be loaded.}

\item{lib_loc}{Character vector, optional.
Same as in \code{library()} or \code{require()}.
Applies to the \code{packages} argument (see above).}

\item{prework}{Expression (language object), list of expressions,
or character vector.
Code to run right before targets build.
Called only once if \code{parallelism} is \code{"loop"}
and once per target otherwise.
This code can be used to set global options, etc.}

\item{prepend}{Deprecated.}

\item{command}{Deprecated.}

\item{args}{Deprecated.}

\item{recipe_command}{Deprecated.}

\item{log_progress}{Logical, whether to log the progress
of individual targets as they are being built. Progress logging
creates extra files in the cache (usually the \verb{.drake/} folder)
and slows down \code{make()} a little.
If you need to reduce or limit the number of files in the cache,
call \code{make(log_progress = FALSE, recover = FALSE)}.}

\item{skip_targets}{Logical, whether to skip building the targets
in \code{plan} and just import objects and files.}

\item{timeout}{\code{deprecated}. Use \code{elapsed} and \code{cpu} instead.}

\item{cpu}{Same as the \code{cpu} argument of \code{setTimeLimit()}.
Seconds of cpu time before a target times out.
Assign target-level cpu timeout times with an optional \code{cpu}
column in \code{plan}.}

\item{elapsed}{Same as the \code{elapsed} argument of \code{setTimeLimit()}.
Seconds of elapsed time before a target times out.
Assign target-level elapsed timeout times with an optional \code{elapsed}
column in \code{plan}.}

\item{retries}{Number of retries to execute if the target fails.
Assign target-level retries with an optional \code{retries}
column in \code{plan}.}

\item{force}{Logical. If \code{FALSE} (default) then \code{drake}
imposes checks if the cache was created with an old
and incompatible version of drake.
If there is an incompatibility, \code{make()} stops to
give you an opportunity to
downgrade \code{drake} to a compatible version
rather than rerun all your targets from scratch.}

\item{graph}{Deprecated.}

\item{trigger}{Name of the trigger to apply to all targets.
Ignored if \code{plan} has a \code{trigger} column.
See \code{\link[=trigger]{trigger()}} for details.}

\item{skip_imports}{Logical, whether to totally neglect to
process the imports and jump straight to the targets. This can be useful
if your imports are massive and you just want to test your project,
but it is bad practice for reproducible data analysis.
This argument is overridden if you supply your own \code{graph} argument.}

\item{skip_safety_checks}{Logical, whether to skip the safety checks
on your workflow. Use at your own peril.}

\item{config}{Deprecated.}

\item{lazy_load}{An old feature, currently being questioned.
For the current recommendations on memory management, see
\url{https://books.ropensci.org/drake/memory.html#memory-strategies}.
The \code{lazy_load} argument is either a character vector or a logical.
For dynamic targets, the behavior is always \code{"eager"} (see below).
So the \code{lazy_load} argument is for static targets only.
Choices for \code{lazy_load}:
\itemize{
\item \code{"eager"}: no lazy loading. The target is loaded right away
with \code{\link[=assign]{assign()}}.
\item \code{"promise"}: lazy loading with \code{\link[=delayedAssign]{delayedAssign()}}
\item \code{"bind"}: lazy loading with active bindings:
\code{bindr::populate_env()}.
\item \code{TRUE}: same as \code{"promise"}.
\item \code{FALSE}: same as \code{"eager"}.
}

If \code{lazy_load} is \code{"eager"},
drake prunes the execution environment before each target/stage,
removing all superfluous targets
and then loading any dependencies it will need for building.
In other words, drake prepares the environment in advance
and tries to be memory efficient.
If \code{lazy_load} is \code{"bind"} or \code{"promise"}, drake assigns
promises to load any dependencies at the last minute.
Lazy loading may be more memory efficient in some use cases, but
it may duplicate the loading of dependencies, costing time.}

\item{session_info}{Logical, whether to save the \code{sessionInfo()}
to the cache. Defaults to \code{TRUE}.
This behavior is recommended for serious \code{\link[=make]{make()}}s
for the sake of reproducibility. This argument only exists to
speed up tests. Apparently, \code{sessionInfo()} is a bottleneck
for small \code{\link[=make]{make()}}s.}

\item{cache_log_file}{Name of the CSV cache log file to write.
If \code{TRUE}, the default file name is used (\code{drake_cache.CSV}).
If \code{NULL}, no file is written.
If activated, this option writes a flat text file
to represent the state of the cache
(fingerprints of all the targets and imports).
If you put the log file under version control, your commit history
will give you an easy representation of how your results change
over time as the rest of your project changes. Hopefully,
this is a step in the right direction for data reproducibility.}

\item{seed}{Integer, the root pseudo-random number generator
seed to use for your project.
In \code{\link[=make]{make()}}, \code{drake} generates a unique
local seed for each target using the global seed
and the target name. That way, different pseudo-random numbers
are generated for different targets, and this pseudo-randomness
is reproducible.

To ensure reproducibility across different R sessions,
\code{set.seed()} and \code{.Random.seed} are ignored and have no affect on
\code{drake} workflows. Conversely, \code{make()} does not usually
change \code{.Random.seed},
even when pseudo-random numbers are generated.
The exception to this last point is
\code{make(parallelism = "clustermq")}
because the \code{clustermq} package needs to generate random numbers
to set up ports and sockets for ZeroMQ.

On the first call to \code{make()} or \code{drake_config()}, \code{drake}
uses the random number generator seed from the \code{seed} argument.
Here, if the \code{seed} is \code{NULL} (default), \code{drake} uses a \code{seed} of \code{0}.
On subsequent \code{make()}s for existing projects, the project's
cached seed will be used in order to ensure reproducibility.
Thus, the \code{seed} argument must either be \code{NULL} or the same
seed from the project's cache (usually the \verb{.drake/} folder).
To reset the random number generator seed for a project,
use \code{clean(destroy = TRUE)}.}

\item{caching}{Character string, either \code{"main"} or \code{"worker"}.
\itemize{
\item \code{"main"}: Targets are built by remote workers and sent back to
the main process. Then, the main process saves them to the
cache (\code{config$cache}, usually a file system \code{storr}).
Appropriate if remote workers do not have access to the file system
of the calling R session. Targets are cached one at a time,
which may be slow in some situations.
\item \code{"worker"}: Remote workers not only build the targets, but also
save them to the cache. Here, caching happens in parallel.
However, remote workers need to have access to the file system
of the calling R session. Transferring target data across
a network can be slow.
}}

\item{keep_going}{Logical, whether to still keep running \code{\link[=make]{make()}}
if targets fail.}

\item{session}{Deprecated. Has no effect now.}

\item{pruning_strategy}{Deprecated. See \code{memory_strategy}.}

\item{makefile_path}{Deprecated.}

\item{console_log_file}{Deprecated in favor of \code{log_make}.}

\item{ensure_workers}{Deprecated.}

\item{garbage_collection}{Logical, whether to call \code{gc()} each time
a target is built during \code{\link[=make]{make()}}.}

\item{template}{A named list of values to fill in the \code{{{ ... }}}
placeholders in template files (e.g. from \code{\link[=drake_hpc_template_file]{drake_hpc_template_file()}}).
Same as the \code{template} argument of \code{clustermq::Q()} and
\code{clustermq::workers}.
Enabled for \code{clustermq} only (\code{make(parallelism = "clustermq")}),
not \code{future} or \code{batchtools} so far.
For more information, see the \code{clustermq} package:
\url{https://github.com/mschubert/clustermq}.
Some template placeholders such as \code{{{ job_name }}} and \code{{{ n_jobs }}}
cannot be set this way.}

\item{sleep}{Optional function on a single numeric argument \code{i}.
Default: \code{function(i) 0.01}.

To conserve memory, \code{drake} assigns a brand new closure to
\code{sleep}, so your custom function should not depend on in-memory data
except from loaded packages.

For parallel processing, \code{drake} uses
a central main process to check what the parallel
workers are doing, and for the affected high-performance
computing workflows, wait for data to arrive over a network.
In between loop iterations, the main process sleeps to avoid throttling.
The \code{sleep} argument to \code{make()} and \code{drake_config()}
allows you to customize how much time the main process spends
sleeping.

The \code{sleep} argument is a function that takes an argument
\code{i} and returns a numeric scalar, the number of seconds to
supply to \code{Sys.sleep()} after iteration \code{i} of checking.
(Here, \code{i} starts at 1.)
If the checking loop does something other than sleeping
on iteration \code{i}, then \code{i} is reset back to 1.

To sleep for the same amount of time between checks,
you might supply something like \code{function(i) 0.01}.
But to avoid consuming too many resources during heavier
and longer workflows, you might use an exponential
back-off: say,
\code{function(i) { 0.1 + 120 * pexp(i - 1, rate = 0.01) }}.}

\item{hasty_build}{Deprecated}

\item{memory_strategy}{Character scalar, name of the
strategy \code{drake} uses to load/unload a target's dependencies in memory.
You can give each target its own memory strategy,
(e.g. \code{drake_plan(x = 1, y = target(f(x), memory_strategy = "lookahead"))})
to override the global memory strategy. Choices:
\itemize{
\item \code{"speed"}: Once a target is newly built or loaded in memory,
just keep it there.
This choice maximizes speed and hogs memory.
\item \code{"autoclean"}: Just before building each new target,
unload everything from memory except the target's direct dependencies.
After a target is built, discard it from memory.
(Set \code{garbage_collection = TRUE} to make sure it is really gone.)
This option conserves memory, but it sacrifices speed because
each new target needs to reload
any previously unloaded targets from storage.
\item \code{"preclean"}: Just before building each new target,
unload everything from memory except the target's direct dependencies.
After a target is built, keep it in memory until \code{drake} determines
they can be unloaded.
This option conserves memory, but it sacrifices speed because
each new target needs to reload
any previously unloaded targets from storage.
\item \code{"lookahead"}: Just before building each new target,
search the dependency graph to find targets that will not be
needed for the rest of the current \code{make()} session.
After a target is built, keep it in memory until the next
memory management stage.
In this mode, targets are only in memory if they need to be loaded,
and we avoid superfluous reads from the cache.
However, searching the graph takes time,
and it could even double the computational overhead for large projects.
\item \code{"unload"}: Just before building each new target,
unload all targets from memory.
After a target is built, \strong{do not} keep it in memory.
This mode aggressively optimizes for both memory and speed,
but in commands and triggers,
you have to manually load any dependencies you need using \code{readd()}.
\item \code{"none"}: Do not manage memory at all.
Do not load or unload anything before building targets.
After a target is built, \strong{do not} keep it in memory.
This mode aggressively optimizes for both memory and speed,
but in commands and triggers,
you have to manually load any dependencies you need using \code{readd()}.
}

For even more direct
control over which targets \code{drake} keeps in memory, see the
help file examples of \code{\link[=drake_envir]{drake_envir()}}.
Also see the \code{garbage_collection} argument of \code{make()} and
\code{drake_config()}.}

\item{layout}{Deprecated.}

\item{spec}{Deprecated.}

\item{lock_envir}{Logical, whether to lock \code{config$envir} during \code{make()}.
If \code{TRUE}, \code{make()} quits in error whenever a command in your
\code{drake} plan (or \code{prework}) tries to add, remove, or modify
non-hidden variables in your environment/workspace/R session.
This is extremely important for ensuring the purity of your functions
and the reproducibility/credibility/trust you can place in your project.
\code{lock_envir} will be set to a default of \code{TRUE} in \code{drake} version
7.0.0 and higher. Namespaces are never locked, e.g.
if \code{envir} is \code{getNamespace("packagename")}.}

\item{history}{Logical, whether to record the build history
of your targets. You can also supply a
\href{https://github.com/wlandau/txtq}{\code{txtq}}, which is
how \code{drake} records history.
Must be \code{TRUE} for \code{\link[=drake_history]{drake_history()}} to work later.}

\item{recover}{Logical, whether to activate automated data recovery.
The default is \code{FALSE} because
\enumerate{
\item Automated data recovery is still stable.
\item It has reproducibility issues.
Targets recovered from the distant past may have been generated
with earlier versions of R and earlier package environments
that no longer exist.
\item It is not always possible, especially when dynamic files
are combined with dynamic branching
(e.g. \code{dynamic = map(stuff)} and \code{format = "file"} etc.)
since behavior is harder to predict in advance.
}

How it works: if \code{recover} is \code{TRUE},
\code{drake} tries to salvage old target values from the cache
instead of running commands from the plan.
A target is recoverable if
\enumerate{
\item There is an old value somewhere in the cache that
shares the command, dependencies, etc.
of the target about to be built.
\item The old value was generated with \code{make(recoverable = TRUE)}.
}

If both conditions are met, \code{drake} will
\enumerate{
\item Assign the most recently-generated admissible data to the target, and
\item skip the target's command.
}

Functions \code{\link[=recoverable]{recoverable()}} and \code{\link[=r_recoverable]{r_recoverable()}} show the most upstream
outdated targets that will be recovered in this way in the next
\code{\link[=make]{make()}} or \code{\link[=r_make]{r_make()}}.}

\item{recoverable}{Logical, whether to make target values recoverable
with \code{make(recover = TRUE)}.
This requires writing extra files to the cache,
and it prevents old metadata from being removed with garbage collection
(\code{clean(garbage_collection = TRUE)}, \code{gc()} in \code{storr}s).
If you need to limit the cache size or the number of files in the cache,
consider \code{make(recoverable = FALSE, progress = FALSE)}.
Recovery is not always possible, especially when dynamic files
are combined with dynamic branching
(e.g. \code{dynamic = map(stuff)} and \code{format = "file"} etc.)
since behavior is harder to predict in advance.}

\item{curl_handles}{A named list of curl handles. Each value is an
object from \code{curl::new_handle()}, and each name is a URL
(and should start with "http", "https", or "ftp").
Example:
list(
\verb{http://httpbin.org/basic-auth} = curl::new_handle(
username = "user", password = "passwd"
)
)
Then, if your plan has
\code{file_in("http://httpbin.org/basic-auth/user/passwd")}
\code{drake} will authenticate using the username and password of the handle
for \verb{http://httpbin.org/basic-auth/}.

\code{drake} uses partial matching on text to
find the right handle of the \code{file_in()} URL, so the name of the handle
could be the complete URL (\code{"http://httpbin.org/basic-auth/user/passwd"})
or a part of the URL (e.g. \code{"http://httpbin.org/"} or
\code{"http://httpbin.org/basic-auth/"}). If you have multiple handles
whose names match your URL, \code{drake} will choose the closest match.}

\item{max_expand}{Positive integer, optional.
\code{max_expand} is the maximum number of targets to generate in each
\code{map()}, \code{cross()}, or \code{group()} dynamic transform.
Useful if you have a massive number of dynamic sub-targets and you want to
work with only the first few sub-targets before scaling up.
Note: the \code{max_expand} argument of \code{make()} and
\code{drake_config()} is for dynamic branching only.
The static branching \code{max_expand}
is an argument of \code{drake_plan()} and \code{transform_plan()}.}

\item{log_build_times}{Logical, whether to record build_times for targets.
Mac users may notice a 20\% speedup in \code{make()}
with \code{build_times = FALSE}.}

\item{format}{Character, an optional custom storage format for targets
without an explicit \code{target(format = ...)} in the plan. Details
about formats:
\url{https://books.ropensci.org/drake/plans.html#special-data-formats-for-targets} # nolint}

\item{lock_cache}{Logical, whether to lock the cache before running \code{make()}
etc. It is usually recommended to keep cache locking on.
However, if you interrupt \code{make()} before it can clean itself up,
then the cache will stay locked,
and you will need to manually unlock it with
\code{drake::drake_cache("xyz")$unlock()}. Repeatedly unlocking the cache
by hand is annoying, and \code{lock_cache = FALSE} prevents the cache
from locking in the first place.}

\item{log_make}{Optional character scalar of a file name or
connection object (such as \code{stdout()}) to dump maximally verbose
log information for \code{\link[=make]{make()}} and other functions (all functions that
accept a \code{config} argument, plus \code{drake_config()}).
If you choose to use a text file as the console log,
it will persist over multiple function calls
until you delete it manually.
Fields in each row the log file, from left to right:
- The node name (short host name) of the
computer (from \code{Sys.info()["nodename"]}).
- The process ID (from \code{Sys.getpid()}).
- A timestamp with the date and time (in microseconds).
- A brief description of what \code{drake} was doing.\verb{ The fields are separated by pipe symbols (}"|"`).}

\item{log_worker}{Logical, same as the \code{log_worker} argument of
\code{clustermq::workers()} and \code{clustermq::Q()}. Only relevant
if \code{parallelism} is \code{"clustermq"}.}
}
\value{
nothing
}
\description{
This is the central, most important function
of the drake package. It runs all the steps of your
workflow in the correct order, skipping any work
that is already up to date. Because of how \code{make()}
tracks global functions and objects as dependencies of targets,
please restart your R session so the pipeline runs
in a clean reproducible environment.
}
\section{Interactive mode}{

In interactive sessions, consider \code{\link[=r_make]{r_make()}}, \code{\link[=r_outdated]{r_outdated()}}, etc.
rather than \code{\link[=make]{make()}}, \code{\link[=outdated]{outdated()}}, etc. The \verb{r_*()} \code{drake} functions
are more reproducible when the session is interactive.
If you do run \code{make()} interactively, please restart your R session
beforehand so your functions and global objects get loaded into
a clean reproducible environment. This prevents targets
from getting invalidated unexpectedly.

A serious drake workflow should be consistent and reliable,
ideally with the help of a main R script.
This script should begin in a fresh R session,
load your packages and functions in a dependable manner,
and then run \code{make()}. Example:
\url{https://github.com/wlandau/drake-examples/tree/main/gsp}.
Batch mode, especially within a container, is particularly helpful.

Interactive R sessions are still useful,
but they easily grow stale.
Targets can falsely invalidate if you accidentally change
a function or data object in your environment.
}

\section{Self-invalidation}{

It is possible to construct a workflow that tries to invalidate itself.
Example:\if{html}{\out{<div class="r">}}\preformatted{plan <- drake_plan(
  x = \{
    data(mtcars)
    mtcars$mpg
  \},
  y = mean(x)
)
}\if{html}{\out{</div>}}

Here, because \code{data()} loads \code{mtcars} into the global environment,
the very act of building \code{x} changes the dependencies of \code{x}.
In other words, without safeguards, \code{x} would not be up to date at
the end of \code{make(plan)}.
Please try to avoid workflows that modify the global environment.
Functions such as \code{data()} belong in your setup scripts
prior to \code{make()}, not in any functions or commands that get called
during \code{make()} itself.

For each target that is still problematic  (e.g.
\url{https://github.com/rstudio/gt/issues/297})
you can safely run the command in its own special \code{callr::r()} process.
Example: \url{https://github.com/rstudio/gt/issues/297#issuecomment-497778735}. # nolint

If that fails, you can run \code{make(plan, lock_envir = FALSE)}
to suppress environment-locking for all targets.
However, this is not usually recommended.
There are legitimate use cases for \code{lock_envir = FALSE}
(example: \url{https://books.ropensci.org/drake/hpc.html#parallel-computing-within-targets}) # nolint
but most workflows should stick with the default \code{lock_envir = TRUE}.
}

\section{Cache locking}{

When \code{make()} runs, it locks the cache so other processes cannot modify it.
Same goes for \code{\link[=outdated]{outdated()}}, \code{\link[=vis_drake_graph]{vis_drake_graph()}}, and similar functions
when \code{make_imports = TRUE}. This is a safety measure to prevent simultaneous
processes from corrupting the cache. If you get an error saying that the
cache is locked, either set \code{make_imports = FALSE} or manually force
unlock it with \code{drake_cache()$unlock()}.
}

\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
config <- drake_config(my_plan)
outdated(my_plan) # Which targets need to be (re)built?
make(my_plan) # Build what needs to be built.
outdated(my_plan) # Everything is up to date.
# Change one of your imported function dependencies.
reg2 = function(d) {
  d$x3 = d$x^3
  lm(y ~ x3, data = d)
}
outdated(my_plan) # Some targets depend on reg2().
make(my_plan) # Rebuild just the outdated targets.
outdated(my_plan) # Everything is up to date again.
if (requireNamespace("visNetwork", quietly = TRUE)) {
vis_drake_graph(my_plan) # See how they fit in an interactive graph.
make(my_plan, cache_log_file = TRUE) # Write a CSV log file this time.
vis_drake_graph(my_plan) # The colors changed in the graph.
# Run targets in parallel:
# options(clustermq.scheduler = "multicore") # nolint
# make(my_plan, parallelism = "clustermq", jobs = 2) # nolint
}
clean() # Start from scratch next time around.
}
# Dynamic branching
# Get the mean mpg for each cyl in the mtcars dataset.
plan <- drake_plan(
  raw = mtcars,
  group_index = raw$cyl,
  munged = target(raw[, c("mpg", "cyl")], dynamic = map(raw)),
  mean_mpg_by_cyl = target(
    data.frame(mpg = mean(munged$mpg), cyl = munged$cyl[1]),
    dynamic = group(munged, .by = group_index)
  )
)
make(plan)
readd(mean_mpg_by_cyl)
})
}
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}},
\code{\link[=drake_config]{drake_config()}},
\code{\link[=vis_drake_graph]{vis_drake_graph()}},
\code{\link[=outdated]{outdated()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{knitr_in}
\alias{knitr_in}
\title{Declare \code{knitr}/\code{rmarkdown} source files
as dependencies.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
knitr_in(...)
}
\arguments{
\item{...}{Character strings. File paths of \code{knitr}/\code{rmarkdown}
source files supplied to a command in your workflow plan data frame.}
}
\value{
A character vector of declared input file paths.
}
\description{
\code{knitr_in()} marks individual \code{knitr}/R Markdown
reports as dependencies. In \code{drake}, these reports are pieces
of the pipeline. R Markdown is a great tool for \emph{displaying}
precomputed results, but not for running a large workflow
from end to end. These reports should do as little
computation as possible.
}
\details{
Unlike \code{\link[=file_in]{file_in()}} and \code{\link[=file_out]{file_out()}}, \code{knitr_in()}
does not work with entire directories.
}
\section{Keywords}{

\code{\link[=drake_plan]{drake_plan()}} understands special keyword functions for your commands.
With the exception of \code{\link[=target]{target()}}, each one is a proper function
with its own help file.
\itemize{
\item \code{\link[=target]{target()}}: give the target more than just a command.
Using \code{\link[=target]{target()}}, you can apply a transformation
(examples: \url{https://books.ropensci.org/drake/plans.html#large-plans}), # nolint
supply a trigger (\url{https://books.ropensci.org/drake/triggers.html}), # nolint
or set any number of custom columns.
\item \code{\link[=file_in]{file_in()}}: declare an input file dependency.
\item \code{\link[=file_out]{file_out()}}: declare an output file to be produced
when the target is built.
\item \code{\link[=knitr_in]{knitr_in()}}: declare a \code{knitr} file dependency such as an
R Markdown (\verb{*.Rmd}) or R LaTeX (\verb{*.Rnw}) file.
\item \code{\link[=ignore]{ignore()}}: force \code{drake} to entirely ignore a piece of code:
do not track it for changes and do not analyze it for dependencies.
\item \code{\link[=no_deps]{no_deps()}}: tell \code{drake} to not track the dependencies
of a piece of code. \code{drake} still tracks the code itself for changes.
\item \code{\link[=id_chr]{id_chr()}}: Get the name of the current target.
\item \code{\link[=drake_envir]{drake_envir()}}: get the environment where drake builds targets.
Intended for advanced custom memory management.
}
}

\examples{
\dontrun{
isolate_example("contain side effects", {
if (requireNamespace("knitr", quietly = TRUE)) {
# `knitr_in()` is like `file_in()`
# except that it analyzes active code chunks in your `knitr`
# source file and detects non-file dependencies.
# That way, updates to the right dependencies trigger rebuilds
# in your report.
# The mtcars example (`drake_example("mtcars")`)
# already has a demonstration

load_mtcars_example()
make(my_plan)

# Now how did drake magically know that
# `small`, `large`, and `coef_regression2_small` were
# dependencies of the output file `report.md`?
# because the command in the workflow plan had
# `knitr_in("report.Rmd")` in it, so drake knew
# to analyze the active code chunks. There, it spotted
# where `small`, `large`, and `coef_regression2_small`
# were read from the cache using calls to `loadd()` and `readd()`.
}
})
}
}
\seealso{
\code{\link[=file_in]{file_in()}}, \code{\link[=file_out]{file_out()}}, \code{\link[=ignore]{ignore()}}, \code{\link[=no_deps]{no_deps()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{expose_imports}
\alias{expose_imports}
\title{Deprecated: expose package functions and objects for
analysis with drake.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
expose_imports(
  package,
  character_only = FALSE,
  envir = parent.frame(),
  jobs = 1
)
}
\arguments{
\item{package}{Name of the package, either a symbol or a string,
depending on \code{character_only}.}

\item{character_only}{Logical, whether to interpret \code{package}
as a character string or a symbol (quoted vs unquoted).}

\item{envir}{Environment to load the exposed package imports.
You will later pass this \code{envir} to \code{\link[=make]{make()}}.}

\item{jobs}{Number of parallel jobs for the parallel processing
of the imports.}
}
\value{
The environment that the exposed imports are loaded into.
Defaults to your R workspace.
}
\description{
Deprecated on 2020-06-24.
}
\details{
Deprecated. This function assigns the objects and functions
from the package environment to the user's environment (usually global)
so \code{drake} can watch them for changes. This used to be the standard
way to make \code{drake} compatible with workflows implemented as custom
analysis packages. Now, the recommendation is to supply
\code{getNamespace("yourPackage")} to the \code{envir} argument of \code{\link[=make]{make()}}
and friends. Read \url{https://github.com/ropensci/drake/issues/1286},
especially \url{https://github.com/ropensci/drake/issues/1286#issuecomment-649088321}, # nolint
for details.
}
\examples{
# nolint start
\dontrun{
isolate_example("contain side effects", {
# Consider a simple plan that depends on the biglm package.
# library(biglm)
plan <- drake_plan(model = biglm(y ~ x, data = huge_dataset))
# Even if you load the biglm package, drake still ignores
# the biglm() function as a dependency. The function is missing
# from the graph:
# vis_drake_graph(plan)
# And if you install an updated version of biglm with a revised
# biglm() function, this will not cause drake::make(plan)
# to rerun the model.
# This is because biglm() is not in your environment.
# ls()
# biglm() exists in its own special package environment,
# which drake does not scan.
# ls("package:biglm")
# To depend on biglm(), use expose_imports(biglm)
# to bring the objects and functions in biglm into
# your own (non-package) environment.
# expose_imports(biglm)
# Now, the biglm() function should be in your environment.
# ls()
# biglm() now appears in the graph.
# vis_drake_graph(plan)
# And subsequent make()s respond to changes to biglm()
# and its dependencies.
})
}
# nolint end
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{config}
\alias{config}
\title{config \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
config(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio.R
\name{rs_addin_loadd}
\alias{rs_addin_loadd}
\title{Loadd target at cursor into global environment
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
rs_addin_loadd(context = NULL)
}
\arguments{
\item{context}{an RStudio document context.
Read from the active document if not supplied.
This is used for testing purposes.}
}
\value{
Nothing.
}
\description{
This function provides an RStudio addin that will
load the target at the
current cursor location from the cache into the global environment.
This is convenient during pipeline development when building off
established targets.
}
\details{
If you are using a non-standard \code{drake} cache,
you must supply it to the \code{"rstudio_drake_cache"} global option,
e.g. \code{options(rstudio_drake_cache = storr::storr_rds("my_cache"))}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_drake_deps.R
\name{drake_deps}
\alias{drake_deps}
\title{\code{drake_deps} helper}
\usage{
drake_deps(expr, exclude = character(0), restrict = NULL)
}
\arguments{
\item{expr}{An R expression}

\item{exclude}{Character vector of the names of symbols to exclude
from the code analysis.}

\item{restrict}{Optional character vector of allowable names of globals.
If \code{NULL}, all global symbols are detectable. If a character vector,
only the variables in \code{restrict} will count as global variables.}
}
\value{
A \code{drake_deps} object.
}
\description{
Static code analysis.
}
\examples{
if (FALSE) { # stronger than roxygen dontrun
expr <- quote({
  a <- base::list(1)
  b <- seq_len(10)
  file_out("abc")
  file_in("xyz")
  x <- "123"
  loadd(abc)
  readd(xyz)
})
drake_deps(expr)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{dataframes_graph}
\alias{dataframes_graph}
\title{dataframes_graph \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
dataframes_graph(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{reduce_by}
\alias{reduce_by}
\title{Reduce multiple groupings of targets
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
reduce_by(
  plan,
  ...,
  prefix = "target",
  begin = "",
  op = " + ",
  end = "",
  pairwise = TRUE,
  append = TRUE,
  filter = NULL,
  sep = "_"
)
}
\arguments{
\item{plan}{Workflow plan data frame of prespecified targets.}

\item{...}{Symbols, columns of \code{plan} to define target groupings.
A \code{reduce_plan()} call is applied for each grouping.
Groupings with all \code{NA}s in the selector variables are ignored.}

\item{prefix}{Character, prefix for naming the new targets.
Suffixes are generated from the values of the columns
specified in \code{...}.}

\item{begin}{Character, code to place at the beginning
of each step in the reduction.}

\item{op}{Binary operator to apply in the reduction}

\item{end}{Character, code to place at the end
of each step in the reduction.}

\item{pairwise}{Logical, whether to create multiple
new targets, one for each pair/step in the reduction (\code{TRUE}),
or to do the reduction all in one command.}

\item{append}{Logical. If \code{TRUE}, the output will include the
original rows in the \code{plan} argument.
If \code{FALSE}, the output will only include the new
targets and commands.}

\item{filter}{An expression like you would pass to \code{dplyr::filter()}.
The rows for which \code{filter} evaluates to \code{TRUE} will be gathered,
and the rest will be excluded from gathering.
Why not just call \code{dplyr::filter()} before \code{gather_by()}?
Because \code{gather_by(append = TRUE, filter = my_column == "my_value")}
gathers on some targets while including all the original targets
in the output. See the examples for a demonstration.}

\item{sep}{Character scalar, delimiter for creating the names
of new targets.}
}
\value{
A workflow plan data frame.
}
\description{
Deprecated on 2019-05-16. Use \code{\link[=drake_plan]{drake_plan()}}
transformations instead. See
\url{https://books.ropensci.org/drake/plans.html#large-plans}
for the details.
}
\details{
Perform several calls to \code{reduce_plan()}
based on groupings from columns in the plan,
and then row-bind the new targets to the plan.
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{load_basic_example}
\alias{load_basic_example}
\title{load_basic_example \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
load_basic_example(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{file_out}
\alias{file_out}
\title{Declare output files and directories.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
file_out(...)
}
\arguments{
\item{...}{Character vector, paths to files and directories. Use
\code{.id_chr} to refer to the current target by name. \code{.id_chr} is not
limited to use in \code{file_in()} and \code{file_out()}.}
}
\value{
A character vector of declared output file or directory paths.
}
\description{
\code{file_out()} marks individual files
(and whole directories) that your targets create.
}
\section{Keywords}{

\code{\link[=drake_plan]{drake_plan()}} understands special keyword functions for your commands.
With the exception of \code{\link[=target]{target()}}, each one is a proper function
with its own help file.
\itemize{
\item \code{\link[=target]{target()}}: give the target more than just a command.
Using \code{\link[=target]{target()}}, you can apply a transformation
(examples: \url{https://books.ropensci.org/drake/plans.html#large-plans}), # nolint
supply a trigger (\url{https://books.ropensci.org/drake/triggers.html}), # nolint
or set any number of custom columns.
\item \code{\link[=file_in]{file_in()}}: declare an input file dependency.
\item \code{\link[=file_out]{file_out()}}: declare an output file to be produced
when the target is built.
\item \code{\link[=knitr_in]{knitr_in()}}: declare a \code{knitr} file dependency such as an
R Markdown (\verb{*.Rmd}) or R LaTeX (\verb{*.Rnw}) file.
\item \code{\link[=ignore]{ignore()}}: force \code{drake} to entirely ignore a piece of code:
do not track it for changes and do not analyze it for dependencies.
\item \code{\link[=no_deps]{no_deps()}}: tell \code{drake} to not track the dependencies
of a piece of code. \code{drake} still tracks the code itself for changes.
\item \code{\link[=id_chr]{id_chr()}}: Get the name of the current target.
\item \code{\link[=drake_envir]{drake_envir()}}: get the environment where drake builds targets.
Intended for advanced custom memory management.
}
}

\examples{
\dontrun{
isolate_example("contain side effects", {
# The `file_out()` and `file_in()` functions
# just takes in strings and returns them.
file_out("summaries.txt")
# Their main purpose is to orchestrate your custom files
# in your workflow plan data frame.
plan <- drake_plan(
  out = write.csv(mtcars, file_out("mtcars.csv")),
  contents = read.csv(file_in("mtcars.csv"))
)
plan
# drake knows "\"mtcars.csv\"" is the first target
# and a dependency of `contents`. See for yourself:

make(plan)
file.exists("mtcars.csv")

 # You may use `.id_chr` inside `file_out()` and `file_in()`
 # to refer  to the current target. This works inside `map()`,
 # `combine()`, `split()`, and `cross()`.

plan <- drake::drake_plan(
  data = target(
    write.csv(data, file_out(paste0(.id_chr, ".csv"))),
    transform = map(data = c(airquality, mtcars))
  )
)

plan

# You can also work with entire directories this way.
# However, in `file_out("your_directory")`, the directory
# becomes an entire unit. Thus, `file_in("your_directory")`
# is more appropriate for subsequent steps than
# `file_in("your_directory/file_inside.txt")`.
plan <- drake_plan(
  out = {
    dir.create(file_out("dir"))
    write.csv(mtcars, "dir/mtcars.csv")
  },
  contents = read.csv(file.path(file_in("dir"), "mtcars.csv"))
)
plan

make(plan)
file.exists("dir/mtcars.csv")

# See the connections that the file relationships create:
if (requireNamespace("visNetwork", quietly = TRUE)) {
  vis_drake_graph(plan)
}
})
}
}
\seealso{
\code{\link[=file_in]{file_in()}}, \code{\link[=knitr_in]{knitr_in()}}, \code{\link[=ignore]{ignore()}}, \code{\link[=no_deps]{no_deps()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deps.R
\name{deps_profile_impl}
\alias{deps_profile_impl}
\title{Internal function with a drake_config() argument}
\usage{
deps_profile_impl(target, config, character_only = FALSE)
}
\arguments{
\item{target}{Name of a target.}

\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}

\item{character_only}{Logical, whether to interpret
\code{target} as a character (\code{TRUE}) or a symbol (\code{FALSE}).}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamic.R
\name{read_trace}
\alias{read_trace}
\title{Read a trace of a dynamic target.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
read_trace(
  trace,
  target,
  cache = drake::drake_cache(path = path),
  path = NULL,
  character_only = FALSE
)
}
\arguments{
\item{trace}{Character, name of the trace
you want to extract. Such trace names are declared
in the \code{.trace} argument of \code{map()}, \code{cross()} or \code{group()}.}

\item{target}{Symbol or character,
depending on the value of \code{character_only}.
\code{target} is T=the name of a dynamic target with one or more traces
defined using the \code{.trace} argument of dynamic \code{map()}, \code{cross()},
or \code{group()}.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{character_only}{Logical, whether \code{name} should be treated
as a character or a symbol
(just like \code{character.only} in \code{\link[=library]{library()}}).}
}
\value{
The dynamic trace of one target in another:
a vector of values from a grouping variable.
}
\description{
Read a target's dynamic trace from the cache.
Best used on its own outside a \code{drake} plan.
}
\details{
In dynamic branching, the trace keeps track
of how the sub-targets were generated.
It reminds us the values of grouping variables
that go with individual sub-targets.
}
\examples{
\dontrun{
isolate_example("demonstrate dynamic trace", {
plan <- drake_plan(
  w = LETTERS[seq_len(3)],
  x = letters[seq_len(2)],

  # The first trace lets us see the values of w
  # that go with the sub-targets of y.
  y = target(paste0(w, x), dynamic = cross(w, x, .trace = w)),

  # We can use the trace as a grouping variable for the next
  # group().
  w_tr = read_trace("w", y),

  # Now, we use the trace again to keep track of the
  # values of w corresponding to the sub-targets of z.
  z = target(
    paste0(y, collapse = "-"),
    dynamic = group(y, .by = w_tr, .trace = w_tr)
  )
)
make(plan)

# We can read the trace outside make().
# That way, we know which values of `w` correspond
# to the sub-targets of `y`.
readd(y)
read_trace("w", y)

# And we know which values of `w_tr` (and thus `w`)
# match up with the sub-targets of `y`.
readd(z)
read_trace("w_tr", z)
})
}
}
\seealso{
\code{\link[=get_trace]{get_trace()}}, \code{\link[=subtargets]{subtargets()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_drake_triggers.R
\name{drake_triggers}
\alias{drake_triggers}
\title{\code{drake_triggers} helper}
\usage{
drake_triggers(
  command = TRUE,
  depend = TRUE,
  file = TRUE,
  seed = TRUE,
  format = TRUE,
  condition = FALSE,
  change = NULL,
  mode = c("whitelist", "blacklist", "condition")
)
}
\arguments{
\item{command}{Logical, whether to rebuild the target if the
\code{\link[=drake_plan]{drake_plan()}} command changes.}

\item{depend}{Logical, whether to rebuild if a
non-file dependency changes.}

\item{file}{Logical, whether to rebuild the target
if a \code{\link[=file_in]{file_in()}}/\code{\link[=file_out]{file_out()}}/\code{\link[=knitr_in]{knitr_in()}} file changes.
Also applies to external data tracked with
\code{target(format = "file")}.}

\item{seed}{Logical, whether to rebuild the target
if the seed changes. Only makes a difference if you set
a custom \code{seed} column in your \code{\link[=drake_plan]{drake_plan()}} at some point
in your workflow.}

\item{format}{Logical, whether to rebuild the target if the
choice of specialized data format changes: for example,
if you use \code{target(format = "qs")} one instance and
\code{target(format = "fst")} the next. See
\url{https://books.ropensci.org/drake/plans.html#special-data-formats-for-targets} # nolint
for details on formats.}

\item{condition}{R code (expression or language object)
that returns a logical. The target will rebuild
if the code evaluates to \code{TRUE}.}

\item{change}{R code (expression or language object)
that returns any value. The target will rebuild
if that value is different from last time
or not already cached.}

\item{mode}{A character scalar equal to \code{"whitelist"} (default) or
\code{"blacklist"} or \code{"condition"}. With the \code{mode} argument, you can choose
how the \code{condition} trigger factors into the decision to build
or skip the target. Here are the options.
\itemize{
\item \code{"whitelist"} (default): we \emph{rebuild} the target whenever \code{condition}
evaluates to \code{TRUE}. Otherwise, we defer to the other triggers.
This behavior is the same as the decision rule described in the
"Details" section of this help file.
\item \code{"blacklist"}: we \emph{skip} the target whenever \code{condition} evaluates
to \code{FALSE}. Otherwise, we defer to the other triggers.
\item \code{"condition"}: here, the \code{condition} trigger is the only decider,
and we ignore all the other triggers. We \emph{rebuild} target whenever
\code{condition} evaluates to \code{TRUE} and \emph{skip} it whenever \code{condition}
evaluates to \code{FALSE}.
}}
}
\description{
Triggers of a target.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{from_plan}
\alias{from_plan}
\title{from_plan \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
from_plan(column)
}
\arguments{
\item{column}{Character, name of a column in your \code{drake} plan.}
}
\description{
The \code{from_plan()} function is now defunct
in order to reduce the demands on memory usage.
}
\details{
2019-03-28
}
\seealso{
\code{\link[=drake_envir]{drake_envir()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{deps}
\alias{deps}
\title{deps \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
deps(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-05-16
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{cached_planned}
\alias{cached_planned}
\title{List targets in both the plan and the cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
cached_planned(
  plan,
  path = NULL,
  cache = drake::drake_cache(path = path),
  namespace = NULL,
  jobs = 1
)
}
\arguments{
\item{plan}{A drake plan.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{namespace}{Character scalar, name of the storr namespace
to use for listing objects.}

\item{jobs}{Number of jobs/workers for parallel processing.}
}
\value{
A character vector of target and sub-target names.
}
\description{
Includes dynamic sub-targets as well.
See examples for details.
}
\examples{
\dontrun{
isolate_example("cache_planned() example", {
plan <- drake_plan(w = 1)
make(plan)
cached_planned(plan)
plan <- drake_plan(
  x = seq_len(2),
  y = target(x, dynamic = map(x))
)
cached_planned(plan)
make(plan)
cached_planned(plan)
cached()
})
}
}
\seealso{
\code{\link[=cached]{cached()}}, \link{cached_unplanned}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_graph_info.R
\name{default_graph_title}
\alias{default_graph_title}
\title{Return the default title for graph visualizations
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#soft-deprecated}{\figure{lifecycle-soft-deprecated.svg}{options: alt='[Soft-deprecated]'}}}{\strong{[Soft-deprecated]}}}
\usage{
default_graph_title()
}
\value{
A character scalar with the default graph title.
}
\description{
For internal use only.
}
\examples{
default_graph_title()
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_ggraph.R
\name{drake_ggraph}
\alias{drake_ggraph}
\title{Visualize the workflow with \code{ggraph}/\code{ggplot2}
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_ggraph(
  ...,
  build_times = "build",
  digits = 3,
  targets_only = FALSE,
  main = NULL,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  make_imports = TRUE,
  from_scratch = FALSE,
  full_legend = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE,
  label_nodes = FALSE,
  transparency = TRUE,
  config = NULL
)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{build_times}{Character string or logical.
If character, the choices are
1. \code{"build"}: runtime of the command plus the time
it take to store the target or import.
2. \code{"command"}: just the runtime of the command.
3. \code{"none"}: no build times.
If logical, \code{build_times} selects whether to show the
times from `build_times(..., type = "build")`` or use
no build times at all. See \code{\link[=build_times]{build_times()}} for details.}

\item{digits}{Number of digits for rounding the build times}

\item{targets_only}{Logical,
whether to skip the imports and only include the
targets in the workflow plan.}

\item{main}{Character string, title of the graph.}

\item{from}{Optional collection of target/import names.
If \code{from} is nonempty,
the graph will restrict itself to
a neighborhood of \code{from}.
Control the neighborhood with
\code{mode} and \code{order}.}

\item{mode}{Which direction to branch out in the graph
to create a neighborhood around \code{from}.
Use \code{"in"} to go upstream,
\code{"out"} to go downstream,
and \code{"all"} to go both ways and disregard
edge direction altogether.}

\item{order}{How far to branch out to create
a neighborhood around \code{from}. Defaults to
as far as possible. If a target is in the neighborhood, then
so are all of its custom \code{\link[=file_out]{file_out()}} files if
\code{show_output_files} is \code{TRUE}.
That means the actual graph order may be slightly greater than
you might expect, but this ensures consistency
between \code{show_output_files = TRUE} and
\code{show_output_files = FALSE}.}

\item{subset}{Optional character vector.
Subset of targets/imports to display in the graph.
Applied after \code{from}, \code{mode}, and \code{order}.
Be advised: edges are only kept for adjacent nodes in \code{subset}.
If you do not select all the intermediate nodes,
edges will drop from the graph.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to increase speed and risk using obsolete information.}

\item{from_scratch}{Logical, whether to assume all the targets
will be made from scratch on the next \code{\link[=make]{make()}}.
Makes all targets outdated, but keeps information about
build progress in previous \code{\link[=make]{make()}}s.}

\item{full_legend}{Logical. If \code{TRUE}, all the node types
are printed in the legend. If \code{FALSE}, only the
node types used are printed in the legend.}

\item{group}{Optional character scalar, name of the column used to
group nodes into columns. All the columns names of your original \code{drake}
plan are choices. The other choices (such as \code{"status"}) are column names
in the \code{nodes} . To group nodes into clusters in the graph,
you must also supply the \code{clusters} argument.}

\item{clusters}{Optional character vector of values to cluster on.
These values must be elements of the column of the \code{nodes} data frame
that you specify in the \code{group} argument to \code{drake_graph_info()}.}

\item{show_output_files}{Logical, whether to include
\code{\link[=file_out]{file_out()}} files in the graph.}

\item{label_nodes}{Logical, whether to label the nodes.
If \code{FALSE}, the graph will not have any text next to the nodes,
which is recommended for large graphs with lots of targets.}

\item{transparency}{Logical, whether to allow transparency in
the rendered graph. Set to \code{FALSE} if you get warnings
like "semi-transparency is not supported on this device".}

\item{config}{Deprecated.}
}
\value{
A \code{ggplot2} object, which you can modify with more layers,
show with \code{plot()}, or save as a file with \code{ggsave()}.
}
\description{
This function requires packages \code{ggplot2} and \code{ggraph}.
Install them with \code{install.packages(c("ggplot2", "ggraph"))}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
load_mtcars_example() # Get the code with drake_example("mtcars").
# Plot the network graph representation of the workflow.
if (requireNamespace("ggraph", quietly = TRUE)) {
  drake_ggraph(my_plan) # Save to a file with `ggplot2::ggsave()`.
}
})
}
}
\seealso{
\code{\link[=vis_drake_graph]{vis_drake_graph()}}, \code{\link[=sankey_drake_graph]{sankey_drake_graph()}},
\code{\link[=render_drake_ggraph]{render_drake_ggraph()}}, \code{\link[=text_drake_graph]{text_drake_graph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{this_cache}
\alias{this_cache}
\title{Get the cache at the exact file path specified.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
this_cache(
  path = NULL,
  force = FALSE,
  verbose = 1L,
  fetch_cache = NULL,
  console_log_file = NULL
)
}
\arguments{
\item{path}{File path of the cache.}

\item{force}{Deprecated.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{fetch_cache}{Deprecated.}

\item{console_log_file}{Deprecated in favor of \code{log_make}.}
}
\value{
A drake/storr cache at the specified path, if it exists.
}
\description{
This function does not apply to
in-memory caches such as \code{storr_environment()}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\docType{package}
\name{drake-package}
\alias{drake-package}
\alias{drake}
\title{drake: A pipeline toolkit for reproducible computation at scale.}
\description{
drake is a pipeline toolkit
(\url{https://github.com/pditommaso/awesome-pipeline})
and a scalable, R-focused solution for reproducibility
and high-performance computing.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
library(drake)
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Build everything.
plot(my_plan) # fast call to vis_drake_graph()
make(my_plan) # Nothing is done because everything is already up to date.
reg2 = function(d) { # Change one of your functions.
  d$x3 = d$x^3
  lm(y ~ x3, data = d)
}
make(my_plan) # Only the pieces depending on reg2() get rebuilt.
# Write a flat text log file this time.
make(my_plan, cache_log_file = TRUE)
# Read/load from the cache.
readd(small)
loadd(large)
head(large)
}
# Dynamic branching
# Get the mean mpg for each cyl in the mtcars dataset.
plan <- drake_plan(
  raw = mtcars,
  group_index = raw$cyl,
  munged = target(raw[, c("mpg", "cyl")], dynamic = map(raw)),
  mean_mpg_by_cyl = target(
    data.frame(mpg = mean(munged$mpg), cyl = munged$cyl[1]),
    dynamic = group(munged, .by = group_index)
  )
)
make(plan)
readd(mean_mpg_by_cyl)
})
}
}
\references{
\url{https://github.com/ropensci/drake}
}
\author{
William Michael Landau \email{will.landau@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_build.R
\name{drake_build}
\alias{drake_build}
\title{Build/process a single target or import.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#questioning}{\figure{lifecycle-questioning.svg}{options: alt='[Questioning]'}}}{\strong{[Questioning]}}}
\usage{
drake_build(
  target,
  ...,
  meta = NULL,
  character_only = FALSE,
  replace = FALSE,
  config = NULL
)
}
\arguments{
\item{target}{Name of the target.}

\item{...}{Arguments to \code{\link[=make]{make()}}, such as the plan and environment.}

\item{meta}{Deprecated.}

\item{character_only}{Logical, whether \code{name} should be treated
as a character or a symbol
(just like \code{character.only} in \code{\link[=library]{library()}}).}

\item{replace}{Logical. If \code{FALSE},
items already in your environment
will not be replaced.}

\item{config}{Deprecated 2019-12-22.}
}
\value{
The value of the target right after it is built.
}
\description{
Not valid for dynamic branching.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
# This example is not really a user-side demonstration.
# It just walks through a dive into the internals.
# Populate your workspace and write 'report.Rmd'.
load_mtcars_example() # Get the code with drake_example("mtcars").
out <- drake_build(small, my_plan)
# Now includes `small`.
cached()
head(readd(small))
# `small` was invisibly returned.
head(out)
}
})
}
}
\seealso{
\code{\link[=drake_debug]{drake_debug()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{trigger}
\alias{trigger}
\title{Customize the decision rules for rebuilding targets
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
trigger(
  command = TRUE,
  depend = TRUE,
  file = TRUE,
  seed = TRUE,
  format = TRUE,
  condition = FALSE,
  change = NULL,
  mode = c("whitelist", "blacklist", "condition")
)
}
\arguments{
\item{command}{Logical, whether to rebuild the target if the
\code{\link[=drake_plan]{drake_plan()}} command changes.}

\item{depend}{Logical, whether to rebuild if a
non-file dependency changes.}

\item{file}{Logical, whether to rebuild the target
if a \code{\link[=file_in]{file_in()}}/\code{\link[=file_out]{file_out()}}/\code{\link[=knitr_in]{knitr_in()}} file changes.
Also applies to external data tracked with
\code{target(format = "file")}.}

\item{seed}{Logical, whether to rebuild the target
if the seed changes. Only makes a difference if you set
a custom \code{seed} column in your \code{\link[=drake_plan]{drake_plan()}} at some point
in your workflow.}

\item{format}{Logical, whether to rebuild the target if the
choice of specialized data format changes: for example,
if you use \code{target(format = "qs")} one instance and
\code{target(format = "fst")} the next. See
\url{https://books.ropensci.org/drake/plans.html#special-data-formats-for-targets} # nolint
for details on formats.}

\item{condition}{R code (expression or language object)
that returns a logical. The target will rebuild
if the code evaluates to \code{TRUE}.}

\item{change}{R code (expression or language object)
that returns any value. The target will rebuild
if that value is different from last time
or not already cached.}

\item{mode}{A character scalar equal to \code{"whitelist"} (default) or
\code{"blacklist"} or \code{"condition"}. With the \code{mode} argument, you can choose
how the \code{condition} trigger factors into the decision to build
or skip the target. Here are the options.
\itemize{
\item \code{"whitelist"} (default): we \emph{rebuild} the target whenever \code{condition}
evaluates to \code{TRUE}. Otherwise, we defer to the other triggers.
This behavior is the same as the decision rule described in the
"Details" section of this help file.
\item \code{"blacklist"}: we \emph{skip} the target whenever \code{condition} evaluates
to \code{FALSE}. Otherwise, we defer to the other triggers.
\item \code{"condition"}: here, the \code{condition} trigger is the only decider,
and we ignore all the other triggers. We \emph{rebuild} target whenever
\code{condition} evaluates to \code{TRUE} and \emph{skip} it whenever \code{condition}
evaluates to \code{FALSE}.
}}
}
\value{
A list of trigger specification details that
\code{drake} processes internally when it comes time to decide
whether to build the target.
}
\description{
Use this function inside a target's command
in your \code{\link[=drake_plan]{drake_plan()}} or the \code{trigger} argument to
\code{\link[=make]{make()}} or \code{\link[=drake_config]{drake_config()}}.
For details, see the chapter on triggers
in the user manual:
\url{https://books.ropensci.org/drake/triggers.html}
}
\details{
A target always builds if it has not been built before.
Triggers allow you to customize the conditions
under which a pre-existing target \emph{re}builds.
By default, the target will rebuild if and only if:
\itemize{
\item Any of \code{command}, \code{depend}, or \code{file} is \code{TRUE}, or
\item \code{condition} evaluates to \code{TRUE}, or
\item \code{change} evaluates to a value different from last time.
The above steps correspond to the "whitelist" decision rule.
You can select other decision rules with the \code{mode} argument
described in this help file.
On another note, there may be a slight efficiency loss
if you set complex triggers
for \code{change} and/or \code{condition} because
\code{drake} needs to load any required dependencies
into memory before evaluating these triggers.
}
}
\examples{
# A trigger is just a set of decision rules
# to decide whether to build a target.
trigger()
# This trigger will build a target on Tuesdays
# and when the value of an online dataset changes.
trigger(condition = today() == "Tuesday", change = get_online_dataset())
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
# You can use a global trigger argument:
# for example, to always run everything.
make(my_plan, trigger = trigger(condition = TRUE))
make(my_plan, trigger = trigger(condition = TRUE))
# You can also define specific triggers for each target.
plan <- drake_plan(
  x = sample.int(15),
  y = target(
    command = x + 1,
    trigger = trigger(depend = FALSE)
  )
)
# Now, when x changes, y will not.
make(plan)
make(plan)
plan$command[1] <- "sample.int(16)" # change x
make(plan)
}
})
}
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{target}
\alias{target}
\title{Customize a target in \code{\link[=drake_plan]{drake_plan()}}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
target(command = NULL, transform = NULL, dynamic = NULL, ...)
}
\arguments{
\item{command}{The command to build the target.}

\item{transform}{A call to \code{\link[=map]{map()}}, \code{\link[=split]{split()}}, \code{\link[=cross]{cross()}}, or \code{\link[=combine]{combine()}}
to apply a \emph{static} transformation. Details:
\url{https://books.ropensci.org/drake/static.html}}

\item{dynamic}{A call to \code{\link[=map]{map()}}, \code{\link[=cross]{cross()}}, or \code{\link[=group]{group()}}
to apply a \emph{dynamic} transformation. Details:
\url{https://books.ropensci.org/drake/dynamic.html}}

\item{...}{Optional columns of the plan for a given target.
See the Columns section of this help file for a selection
of special columns that \code{drake} understands.}
}
\value{
A one-row workflow plan data frame with the named
arguments as columns.
}
\description{
The \code{target()} function is a way to
configure individual targets in a \code{drake} plan.
Its most common use is to invoke static branching
and dynamic branching, and it can also set the values
of custom columns such as \code{format}, \code{elapsed}, \code{retries},
and \code{max_expand}. Details are at
\url{https://books.ropensci.org/drake/plans.html#special-columns}.
Note: \code{drake_plan(my_target = my_command())}
is equivalent to
\verb{drake_plan(my_target = target(my_command())}.
}
\details{
\code{target()} must be called inside \code{\link[=drake_plan]{drake_plan()}}.
It is invalid otherwise.
}
\section{Columns}{

\code{\link[=drake_plan]{drake_plan()}} creates a special data frame. At minimum, that data frame
must have columns \code{target} and \code{command} with the target names and the
R code chunks to build them, respectively.

You can add custom columns yourself, either with \code{target()} (e.g.
\code{drake_plan(y = target(f(x), transform = map(c(1, 2)), format = "fst"))})
or by appending columns post-hoc (e.g. \code{plan$col <- vals}).

Some of these custom columns are special. They are optional,
but \code{drake} looks for them at various points in the workflow.
\itemize{
\item \code{transform}: a call to \code{\link[=map]{map()}}, \code{\link[=split]{split()}}, \code{\link[=cross]{cross()}}, or
\code{\link[=combine]{combine()}} to create and manipulate large collections of targets.
Details: (\url{https://books.ropensci.org/drake/plans.html#large-plans}). # nolint
\item \code{format}: set a storage format to save big targets more efficiently.
See the "Formats" section of this help file for more details.
\item \code{trigger}: rule to decide whether a target needs to run.
It is recommended that you define this one with \code{target()}.
Details: \url{https://books.ropensci.org/drake/triggers.html}.
\item \code{hpc}: logical values (\code{TRUE}/\code{FALSE}/\code{NA}) whether to send each target
to parallel workers.
Visit \url{https://books.ropensci.org/drake/hpc.html#selectivity}
to learn more.
\item \code{resources}: target-specific lists of resources for a computing cluster.
See
\url{https://books.ropensci.org/drake/hpc.html#advanced-options}
for details.
\item \code{caching}: overrides the \code{caching} argument of \code{\link[=make]{make()}} for each target
individually. Possible values:
\itemize{
\item "main": tell the main process to store the target in the cache.
\item "worker": tell the HPC worker to store the target in the cache.
\item NA: default to the \code{caching} argument of \code{\link[=make]{make()}}.
}
\item \code{elapsed} and \code{cpu}: number of seconds to wait for the target to build
before timing out (\code{elapsed} for elapsed time and \code{cpu} for CPU time).
\item \code{retries}: number of times to retry building a target
in the event of an error.
\item \code{seed}: an optional pseudo-random number generator (RNG)
seed for each target. \code{drake} usually comes up with its own
unique reproducible target-specific seeds using the global seed
(the \code{seed} argument to \code{\link[=make]{make()}} and \code{\link[=drake_config]{drake_config()}})
and the target names, but you can overwrite these automatic seeds.
\code{NA} entries default back to \code{drake}'s automatic seeds.
\item \code{max_expand}: for dynamic branching only. Same as the \code{max_expand}
argument of \code{\link[=make]{make()}}, but on a target-by-target basis.
Limits the number of sub-targets created for a given target.
}
}

\section{Keywords}{

\code{\link[=drake_plan]{drake_plan()}} understands special keyword functions for your commands.
With the exception of \code{\link[=target]{target()}}, each one is a proper function
with its own help file.
\itemize{
\item \code{\link[=target]{target()}}: give the target more than just a command.
Using \code{\link[=target]{target()}}, you can apply a transformation
(examples: \url{https://books.ropensci.org/drake/plans.html#large-plans}), # nolint
supply a trigger (\url{https://books.ropensci.org/drake/triggers.html}), # nolint
or set any number of custom columns.
\item \code{\link[=file_in]{file_in()}}: declare an input file dependency.
\item \code{\link[=file_out]{file_out()}}: declare an output file to be produced
when the target is built.
\item \code{\link[=knitr_in]{knitr_in()}}: declare a \code{knitr} file dependency such as an
R Markdown (\verb{*.Rmd}) or R LaTeX (\verb{*.Rnw}) file.
\item \code{\link[=ignore]{ignore()}}: force \code{drake} to entirely ignore a piece of code:
do not track it for changes and do not analyze it for dependencies.
\item \code{\link[=no_deps]{no_deps()}}: tell \code{drake} to not track the dependencies
of a piece of code. \code{drake} still tracks the code itself for changes.
\item \code{\link[=id_chr]{id_chr()}}: Get the name of the current target.
\item \code{\link[=drake_envir]{drake_envir()}}: get the environment where drake builds targets.
Intended for advanced custom memory management.
}
}

\section{Formats}{

Specialized target formats increase efficiency and flexibility.
Some allow you to save specialized objects like \code{keras} models,
while others increase the speed while conserving storage and memory.
You can declare target-specific formats in the plan
(e.g. \code{drake_plan(x = target(big_data_frame, format = "fst"))})
or supply a global default \code{format} for all targets in \code{make()}.
Either way, most formats have specialized installation requirements
(e.g. R packages) that are not installed with \code{drake} by default.
You will need to install them separately yourself.
Available formats:
\itemize{
\item \code{"file"}: Dynamic files. To use this format, simply create
local files and directories yourself and then return
a character vector of paths as the target's value.
Then, \code{drake} will watch for changes to those files in
subsequent calls to \code{make()}. This is a more flexible
alternative to \code{file_in()} and \code{file_out()}, and it is
compatible with dynamic branching.
See \url{https://github.com/ropensci/drake/pull/1178} for an example.
\item \code{"fst"}: save big data frames fast. Requires the \code{fst} package.
Note: this format strips non-data-frame attributes such as the
\item \code{"fst_tbl"}: Like \code{"fst"}, but for \code{tibble} objects.
Requires the \code{fst} and \code{tibble} packages.
Strips away non-data-frame non-tibble attributes.
\item \code{"fst_dt"}: Like \code{"fst"} format, but for \code{data.table} objects.
Requires the \code{fst} and \code{data.table} packages.
Strips away non-data-frame non-data-table attributes.
\item \code{"diskframe"}:
Stores \code{disk.frame} objects, which could potentially be
larger than memory. Requires the \code{fst} and \code{disk.frame} packages.
Coerces objects to \code{disk.frame}s.
Note: \code{disk.frame} objects get moved to the \code{drake} cache
(a subfolder of \verb{.drake/} for most workflows).
To ensure this data transfer is fast, it is best to
save your \code{disk.frame} objects to the same physical storage
drive as the \code{drake} cache,
\code{as.disk.frame(your_dataset, outdir = drake_tempfile())}.
\item \code{"keras"}: save Keras models as HDF5 files.
Requires the \code{keras} package.
\item \code{"qs"}: save any R object that can be properly serialized
with the \code{qs} package. Requires the \code{qs} package.
Uses \code{qsave()} and \code{qread()}.
Uses the default settings in \code{qs} version 0.20.2.
\item \code{"rds"}: save any R object that can be properly serialized.
Requires R version >= 3.5.0 due to ALTREP.
Note: the \code{"rds"} format uses gzip compression, which is slow.
\code{"qs"} is a superior format.
}
}

\examples{
# Use target() to create your own custom columns in a drake plan.
# See ?triggers for more on triggers.
drake_plan(
  website_data = target(
    download_data("www.your_url.com"),
    trigger = "always",
    custom_column = 5
  ),
  analysis = analyze(website_data)
)
models <- c("glm", "hierarchical")
plan <- drake_plan(
  data = target(
    get_data(x),
    transform = map(x = c("simulated", "survey"))
  ),
  analysis = target(
    analyze_data(data, model),
    transform = cross(data, model = !!models, .id = c(x, model))
  ),
  summary = target(
    summarize_analysis(analysis),
    transform = map(analysis, .id = c(x, model))
  ),
  results = target(
    bind_rows(summary),
    transform = combine(summary, .by = data)
  )
)
plan
if (requireNamespace("styler", quietly = TRUE)) {
  print(drake_plan_source(plan))
}
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{get_cache}
\alias{get_cache}
\title{The default cache of a \code{drake} project.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
get_cache(
  path = getwd(),
  search = TRUE,
  verbose = 1L,
  force = FALSE,
  fetch_cache = NULL,
  console_log_file = NULL
)
}
\arguments{
\item{path}{Character, either the root file path of a \code{drake} project
or a folder containing the root (top-level working directory
where you plan to call \code{\link[=make]{make()}}).
If this is too confusing, feel free to just use \code{storr::storr_rds()}
to get the cache.
If \code{search = FALSE}, \code{path} must be the root.
If \code{search = TRUE}, you can specify any
subdirectory of the project. Let's say \code{"/home/you/my_project"}
is the root. The following are equivalent and correct:
\itemize{
\item \code{get_cache(path = "/home/you/my_project", search = FALSE)}
\item \code{get_cache(path = "/home/you/my_project", search = TRUE)}
\item \code{get_cache(path = "/home/you/my_project/subdir/x", search = TRUE)}
\item \code{get_cache(path = "/home/you/my_project/.drake", search = TRUE)}
\item \code{get_cache(path = "/home/you/my_project/.drake/keys", search = TRUE)}
}}

\item{search}{Deprecated.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{force}{Deprecated.}

\item{fetch_cache}{Deprecated.}

\item{console_log_file}{Deprecated in favor of \code{log_make}.}
}
\description{
Use \code{\link[=drake_cache]{drake_cache()}} instead.
}
\details{
Deprecated on 2019-05-25.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.R
\name{make_impl}
\alias{make_impl}
\title{Internal function with a drake_config() argument}
\usage{
make_impl(config)
}
\arguments{
\item{config}{a \code{\link[=drake_config]{drake_config()}} object.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{recover_cache}
\alias{recover_cache}
\title{Load or create a drake cache
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
recover_cache(
  path = NULL,
  hash_algorithm = NULL,
  short_hash_algo = NULL,
  long_hash_algo = NULL,
  force = FALSE,
  verbose = 1L,
  fetch_cache = NULL,
  console_log_file = NULL
)
}
\arguments{
\item{path}{File path of the cache.}

\item{hash_algorithm}{Name of a hash algorithm to use.
See the \code{algo} argument of the \code{digest} package for your options.}

\item{short_hash_algo}{Deprecated on 2018-12-12.
Use \code{hash_algorithm} instead.}

\item{long_hash_algo}{Deprecated on 2018-12-12.
Use \code{hash_algorithm} instead.}

\item{force}{Logical, whether to load the cache
despite any back compatibility issues with the
running version of drake.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{fetch_cache}{Deprecated.}

\item{console_log_file}{Deprecated on 2019-09-11.}
}
\value{
A drake/storr cache.
}
\description{
Deprecated on 2019-01-13.
}
\details{
Does not work with
in-memory caches such as \code{\link[=storr_environment]{storr_environment()}}.
}
\seealso{
\code{\link[=new_cache]{new_cache()}}, \code{\link[=get_cache]{get_cache()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_example.R
\name{clean_mtcars_example}
\alias{clean_mtcars_example}
\title{Clean the mtcars example from \code{drake_example("mtcars")}
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
clean_mtcars_example()
}
\value{
nothing
}
\description{
This function deletes files. Use at your own risk.
Destroys the \verb{.drake/} cache and the \code{report.Rmd} file
in the current working directory. Your working directory
(\code{getcwd()}) must be the folder from which you first ran
\code{load_mtcars_example()} and \code{make(my_plan)}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
# Populate your workspace and write 'report.Rmd'.
load_mtcars_example() # Get the code: drake_example("mtcars")
# Check the dependencies of an imported function.
deps_code(reg1)
# Check the dependencies of commands in the workflow plan.
deps_code(my_plan$command[1])
deps_code(my_plan$command[4])
# Plot the interactive network visualization of the workflow.
outdated(my_plan) # Which targets are out of date?
# Run the workflow to build all the targets in the plan.
make(my_plan)
outdated(my_plan) # Everything should be up to date.
# For the reg2() model on the small dataset,
# the p-value is so small that there may be an association
# between weight and fuel efficiency after all.
readd(coef_regression2_small)
# Clean up the example.
clean_mtcars_example()
}
})
}
}
\seealso{
\code{\link[=load_mtcars_example]{load_mtcars_example()}}, \code{\link[=clean]{clean()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text_drake_graph.R
\name{text_drake_graph}
\alias{text_drake_graph}
\title{Show a workflow graph as text in your terminal window.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
text_drake_graph(
  ...,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  targets_only = FALSE,
  make_imports = TRUE,
  from_scratch = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE,
  nchar = 1L,
  print = TRUE,
  config = NULL
)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{from}{Optional collection of target/import names.
If \code{from} is nonempty,
the graph will restrict itself to
a neighborhood of \code{from}.
Control the neighborhood with
\code{mode} and \code{order}.}

\item{mode}{Which direction to branch out in the graph
to create a neighborhood around \code{from}.
Use \code{"in"} to go upstream,
\code{"out"} to go downstream,
and \code{"all"} to go both ways and disregard
edge direction altogether.}

\item{order}{How far to branch out to create
a neighborhood around \code{from}. Defaults to
as far as possible. If a target is in the neighborhood, then
so are all of its custom \code{\link[=file_out]{file_out()}} files if
\code{show_output_files} is \code{TRUE}.
That means the actual graph order may be slightly greater than
you might expect, but this ensures consistency
between \code{show_output_files = TRUE} and
\code{show_output_files = FALSE}.}

\item{subset}{Optional character vector.
Subset of targets/imports to display in the graph.
Applied after \code{from}, \code{mode}, and \code{order}.
Be advised: edges are only kept for adjacent nodes in \code{subset}.
If you do not select all the intermediate nodes,
edges will drop from the graph.}

\item{targets_only}{Logical,
whether to skip the imports and only include the
targets in the workflow plan.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to increase speed and risk using obsolete information.}

\item{from_scratch}{Logical, whether to assume all the targets
will be made from scratch on the next \code{\link[=make]{make()}}.
Makes all targets outdated, but keeps information about
build progress in previous \code{\link[=make]{make()}}s.}

\item{group}{Optional character scalar, name of the column used to
group nodes into columns. All the columns names of your original \code{drake}
plan are choices. The other choices (such as \code{"status"}) are column names
in the \code{nodes} . To group nodes into clusters in the graph,
you must also supply the \code{clusters} argument.}

\item{clusters}{Optional character vector of values to cluster on.
These values must be elements of the column of the \code{nodes} data frame
that you specify in the \code{group} argument to \code{drake_graph_info()}.}

\item{show_output_files}{Logical, whether to include
\code{\link[=file_out]{file_out()}} files in the graph.}

\item{nchar}{For each node, maximum number of characters of the node label
to show. Can be 0, in which case each node is a colored box
instead of a node label.
Caution: \code{nchar} > 0 will mess with the layout.}

\item{print}{Logical. If \code{TRUE}, the graph will print to the console
via \code{message()}. If \code{FALSE}, nothing is printed. However, you still
have the visualization because \code{text_drake_graph()} and
\code{render_text_drake_graph()} still invisibly return a character string
that you can print yourself with \code{message()}.}

\item{config}{Deprecated.}
}
\value{
A \code{visNetwork} graph.
}
\description{
This is a low-tech version of \code{\link[=vis_drake_graph]{vis_drake_graph()}}
and friends. It is designed for when you do not have access
to the usual graphics devices for viewing visuals in an interactive
R session: for example, if you are logged into a remote machine
with SSH and you do not have access to X Window support.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
# Plot the network graph representation of the workflow.
pkg <- requireNamespace("txtplot", quietly = TRUE) &&
  requireNamespace("visNetwork", quietly = TRUE)
if (pkg) {
text_drake_graph(my_plan)
make(my_plan) # Run the project, build the targets.
text_drake_graph(my_plan) # The black nodes from before are now green.
}
}
})
}
}
\seealso{
\code{\link[=render_text_drake_graph]{render_text_drake_graph()}}, \code{\link[=vis_drake_graph]{vis_drake_graph()}},
\code{\link[=sankey_drake_graph]{sankey_drake_graph()}}, \code{\link[=drake_ggraph]{drake_ggraph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{built}
\alias{built}
\title{List all the built targets (non-imports) in the cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
built(
  path = getwd(),
  search = TRUE,
  cache = drake::get_cache(path = path, search = search, verbose = verbose),
  verbose = 1L,
  jobs = 1
)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{jobs}{Number of jobs/workers for parallel processing.}
}
\value{
Character vector naming the built targets in the cache.
}
\description{
Deprecated on 2019-01-08.
}
\details{
Targets are listed in the workflow plan
data frame (see \code{\link[=drake_plan]{drake_plan()}}.
}
\seealso{
\code{\link[=cached]{cached()}}, \code{\link[=loadd]{loadd()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outdated.R
\name{recoverable_impl}
\alias{recoverable_impl}
\title{Internal function with a drake_config() argument}
\usage{
recoverable_impl(config = NULL, make_imports = TRUE, do_prework = TRUE)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to save some time and risk obsolete output.}

\item{do_prework}{Whether to do the \code{prework}
normally supplied to \code{\link[=make]{make()}}.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{make_targets}
\alias{make_targets}
\title{Just make the targets
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
make_targets(config)
}
\arguments{
\item{config}{A configuration list returned by \code{\link[=drake_config]{drake_config()}}.}
}
\value{
nothing
}
\description{
Deprecated on 2019-01-04
}
\seealso{
\code{\link[=make]{make()}}, \code{\link[=drake_config]{drake_config()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{gather}
\alias{gather}
\title{gather \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
gather(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_plan.R
\name{transform_plan}
\alias{transform_plan}
\title{Transform a plan
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
transform_plan(
  plan,
  envir = parent.frame(),
  trace = FALSE,
  max_expand = NULL,
  tidy_eval = TRUE
)
}
\arguments{
\item{plan}{A \code{drake} plan with a \code{transform} column}

\item{envir}{Environment for tidy evaluation.}

\item{trace}{Logical, whether to add columns to show
what happens during target transformations.}

\item{max_expand}{Positive integer, optional.
\code{max_expand} is the maximum number of targets to generate in each
\code{map()}, \code{split()}, or \code{cross()} transform.
Useful if you have a massive plan and you want to
test and visualize a strategic subset of targets
before scaling up.
Note: the \code{max_expand} argument of \code{drake_plan()} and
\code{transform_plan()} is for static branching only.
The dynamic branching \code{max_expand}
is an argument of \code{make()} and \code{drake_config()}.}

\item{tidy_eval}{Logical, whether to use tidy evaluation
(e.g. unquoting/\verb{!!}) when resolving commands.
Tidy evaluation in transformations is always turned on
regardless of the value you supply to this argument.}
}
\description{
Evaluate the \code{map()}, \code{cross()}, \code{split()} and
\code{combine()} operations in the \code{transform} column of a
\code{drake} plan.
}
\details{
\url{https://books.ropensci.org/drake/plans.html#large-plans} # nolint
}
\examples{
plan1 <- drake_plan(
  y = target(
    f(x),
    transform = map(x = c(1, 2))
  ),
  transform = FALSE
)
plan2 <- drake_plan(
  z = target(
    g(y),
    transform = map(y, .id = x)
  ),
  transform = FALSE
)
plan <- bind_plans(plan1, plan2)
transform_plan(plan)
models <- c("glm", "hierarchical")
plan <- drake_plan(
  data = target(
    get_data(x),
    transform = map(x = c("simulated", "survey"))
  ),
  analysis = target(
    analyze_data(data, model),
    transform = cross(data, model = !!models, .id = c(x, model))
  ),
  summary = target(
    summarize_analysis(analysis),
    transform = map(analysis, .id = c(x, model))
  ),
  results = target(
    bind_rows(summary),
    transform = combine(summary, .by = data)
  )
)
plan
if (requireNamespace("styler", quietly = TRUE)) {
  print(drake_plan_source(plan))
}
# Tags:
drake_plan(
  x = target(
    command,
    transform = map(y = c(1, 2), .tag_in = from, .tag_out = c(to, out))
  ),
  trace = TRUE
)
plan <- drake_plan(
  survey = target(
    survey_data(x),
    transform = map(x = c(1, 2), .tag_in = source, .tag_out = dataset)
  ),
  download = target(
    download_data(),
    transform = map(y = c(5, 6), .tag_in = source, .tag_out = dataset)
  ),
  analysis = target(
    analyze(dataset),
    transform = map(dataset)
  ),
  results = target(
    bind_rows(analysis),
    transform = combine(analysis, .by = source)
  )
)
plan
if (requireNamespace("styler", quietly = TRUE)) {
  print(drake_plan_source(plan))
}
}
\seealso{
drake_plan, map, split, cross, combine
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{plan_to_notebook}
\alias{plan_to_notebook}
\title{Turn a \code{drake} plan into an R notebook.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#questioning}{\figure{lifecycle-questioning.svg}{options: alt='[Questioning]'}}}{\strong{[Questioning]}}}
\usage{
plan_to_notebook(plan, con)
}
\arguments{
\item{plan}{Workflow plan data frame. See \code{\link[=drake_plan]{drake_plan()}}
for details.}

\item{con}{A file path or connection to write to.}
}
\description{
\code{code_to_plan()}, \code{\link[=plan_to_code]{plan_to_code()}}, and
\code{\link[=plan_to_notebook]{plan_to_notebook()}} together illustrate the relationships
between \code{drake} plans, R scripts, and R Markdown documents.
In the file generated by \code{plan_to_code()}, every target/command pair
becomes a chunk of code.
Targets are arranged in topological order
so dependencies are available before their downstream targets.
Please note:
\enumerate{
\item You are still responsible for loading your project's
packages, imported functions, etc.
\item Triggers disappear.
}
}
\examples{
if (suppressWarnings(require("knitr"))) {
plan <- drake_plan(
  raw_data = read_excel(file_in("raw_data.xlsx")),
  data = raw_data,
  hist = create_plot(data),
  fit = lm(Ozone ~ Temp + Wind, data)
)
file <- tempfile()
# Turn the plan into an R notebook a the given file path.
plan_to_notebook(plan, file)
# Here is what the script looks like.
cat(readLines(file), sep = "\n")
# Convert back to a drake plan.
code_to_plan(file)
}
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}, \code{\link[=code_to_plan]{code_to_plan()}},
\code{\link[=plan_to_code]{plan_to_code()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{find_cache}
\alias{find_cache}
\title{Search up the file system for the nearest drake cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
find_cache(path = getwd(), dir = NULL, directory = NULL)
}
\arguments{
\item{path}{Starting path for search back for the cache.
Should be a subdirectory of the drake project.}

\item{dir}{Character, name of the folder containing the cache.}

\item{directory}{Deprecated. Use \code{dir}.}
}
\value{
File path of the nearest drake cache or \code{NULL}
if no cache is found.
}
\description{
Only works if the cache is a file system in a
hidden folder named \verb{.drake/} (default).
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the target.
# Find the file path of the project's cache.
# Search up through parent directories if necessary.
find_cache()
}
})
}
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}},
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.R
\name{drake_gc}
\alias{drake_gc}
\title{Do garbage collection on the drake cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_gc(
  path = NULL,
  search = NULL,
  verbose = NULL,
  cache = drake::drake_cache(path = path),
  force = FALSE
)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{force}{Logical, whether to load the cache
despite any back compatibility issues with the
running version of drake.}
}
\value{
\code{NULL}
}
\description{
Garbage collection removes obsolete target values
from the cache.
}
\details{
Caution: garbage collection \emph{actually} removes data
so it is no longer recoverable with \code{\link[=drake_history]{drake_history()}} or
\code{make(recover = TRUE)}. You cannot undo this operation.
Use at your own risk.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the targets.
# At this point, check the size of the '.drake/' cache folder.
# Clean without garbage collection.
clean(garbage_collection = FALSE)
# The '.drake/' cache folder is still about the same size.
drake_gc() # Do garbage collection on the cache.
# The '.drake/' cache folder should have gotten much smaller.
}
})
}
}
\seealso{
\code{\link[=clean]{clean()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{r_recipe_wildcard}
\alias{r_recipe_wildcard}
\title{Default Makefile recipe wildcard
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
r_recipe_wildcard()
}
\value{
The R recipe wildcard.
}
\description{
2019-01-02
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{cancel}
\alias{cancel}
\title{Cancel a target mid-build \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
cancel(allow_missing = TRUE)
}
\arguments{
\item{allow_missing}{Logical. If \code{FALSE}, \code{drake} will not cancel
the target if it is missing from the cache (or if you removed the
key with \code{clean()}).}
}
\value{
Nothing.
}
\description{
Cancel a target mid-build.
Upon cancellation, \code{drake} halts the current target and moves to the
next one. The target's previous value and metadata, if they exist,
remain in the cache.
}
\examples{
\dontrun{
isolate_example("cancel()", {
f <- function(x) {
  cancel()
  Sys.sleep(2) # Does not run.
}
g <- function(x) f(x)
plan <- drake_plan(y = g(1))
make(plan)
# Does not exist.
# readd(y)
})
}
}
\seealso{
cancel_if
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_example.R
\name{drake_example}
\alias{drake_example}
\title{Download the files of an example \code{drake} project.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_example(
  example = "main",
  to = getwd(),
  destination = NULL,
  overwrite = FALSE,
  quiet = TRUE
)
}
\arguments{
\item{example}{Name of the example.
The possible values are the names of the folders at
\url{https://github.com/wlandau/drake-examples}.}

\item{to}{Character scalar,
the folder containing the code files for the example.
passed to the \code{exdir} argument of \code{utils::unzip()}.}

\item{destination}{Deprecated; use \code{to} instead.}

\item{overwrite}{Logical, whether to overwrite an existing folder
with the same name as the drake example.}

\item{quiet}{Logical, passed to \code{downloader::download()}
and thus \code{utils::download.file()}. Whether
to download quietly or print progress.}
}
\value{
\code{NULL}
}
\description{
The \code{drake_example()} function downloads a
folder from \url{https://github.com/wlandau/drake-examples}.
By default, it creates a new folder with the example name
in your current working directory. After the files are written,
have a look at the enclosed \code{README} file.
Other instructions are available in the files at
\url{https://github.com/wlandau/drake-examples}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (requireNamespace("downloader")) {
drake_examples() # List all the drake examples.
# Sets up the same example from load_mtcars_example()
drake_example("mtcars")
# Sets up the SLURM example.
drake_example("slurm")
}
})
}
}
\seealso{
\code{\link[=drake_examples]{drake_examples()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{plan_analyses}
\alias{plan_analyses}
\title{Specialized wildcard for analyses
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
plan_analyses(plan, datasets, sep = "_")
}
\arguments{
\item{plan}{Workflow plan data frame of analysis methods.
The commands in the \code{command} column must
have the \code{dataset__} wildcard where the datasets go.
For example, one command could be \code{lm(dataset__)}. Then,
the commands in the output will include \code{lm(your_dataset_1)},
\code{lm(your_dataset_2)}, etc.}

\item{datasets}{Workflow plan data frame with instructions
to make the datasets.}

\item{sep}{character Scalar, delimiter for creating
the names of new targets.}
}
\value{
An evaluated workflow plan data frame of analysis targets.
}
\description{
Use \code{\link[=drake_plan]{drake_plan()}} instead.
See \url{https://books.ropensci.org/drake/plans.html#large-plans}
for details.
}
\details{
2019-01-13
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{drake_palette}
\alias{drake_palette}
\title{Show drake's color palette.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
drake_palette()
}
\value{
There is a console message,
but the actual return value is \code{NULL}.
}
\description{
Deprecated on 2019-01-12.
}
\details{
This function is
used in both the console and graph visualizations.
Your console must have the crayon package enabled.
This palette applies to console output
(internal functions \code{console()} and
\code{console_many_targets()}) and the node colors
in the graph visualizations.
So if you want to contribute improvements to the palette,
please both \code{drake_palette()} and
\code{visNetwork::visNetwork(nodes = legend_nodes())}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{gather_plan}
\alias{gather_plan}
\title{Combine targets
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
gather_plan(plan = NULL, target = "target", gather = "list", append = FALSE)
}
\arguments{
\item{plan}{Workflow plan data frame of prespecified targets.}

\item{target}{Name of the new aggregated target.}

\item{gather}{Function used to gather the targets. Should be
one of \code{list(...)}, \code{c(...)}, \code{rbind(...)}, or similar.}

\item{append}{Logical. If \code{TRUE}, the output will include the
original rows in the \code{plan} argument.
If \code{FALSE}, the output will only include the new
targets and commands.}
}
\value{
A workflow plan data frame that aggregates multiple
prespecified targets into one additional target downstream.
}
\description{
Deprecated on 2019-05-16. Use \code{\link[=drake_plan]{drake_plan()}}
transformations instead. See
\url{https://books.ropensci.org/drake/plans.html#large-plans}
for the details.
}
\details{
Creates a new workflow plan to aggregate
existing targets in the supplied plan.
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{show_source}
\alias{show_source}
\title{Show how a target/import was produced.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
show_source(target, config, character_only = FALSE)
}
\arguments{
\item{target}{Symbol denoting the target or import
or a character vector if character_only is \code{TRUE}.}

\item{config}{A \code{\link[=drake_config]{drake_config()}} list.}

\item{character_only}{Logical, whether to interpret
\code{target} as a symbol (\code{FALSE}) or character vector
(\code{TRUE}).}
}
\description{
Show the command that produced a target
or indicate that the object or file was imported.
}
\examples{
\dontrun{
isolate_example("contain side effects", {
plan <- drake_plan(x = sample.int(15))
cache <- storr::storr_environment() # custom in-memory cache
make(plan, cache = cache)
config <- drake_config(plan, cache = cache, history = FALSE)
show_source(x, config)
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sankey_drake_graph.R
\name{render_sankey_drake_graph}
\alias{render_sankey_drake_graph}
\title{Render a Sankey diagram from \code{\link[=drake_graph_info]{drake_graph_info()}}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
render_sankey_drake_graph(
  graph_info,
  file = character(0),
  selfcontained = FALSE,
  ...
)
}
\arguments{
\item{graph_info}{List of data frames generated by
\code{\link[=drake_graph_info]{drake_graph_info()}}.
There should be 3 data frames: \code{nodes}, \code{edges},
and \code{legend_nodes}.}

\item{file}{Name of a file to save the graph.
If \code{NULL} or \code{character(0)}, no file is saved and
the graph is rendered and displayed within R.
If the file ends in a \code{.png}, \code{.jpg}, \code{.jpeg}, or \code{.pdf} extension,
then a static image will be saved. In this case,
the webshot package and PhantomJS are required:
\verb{install.packages("webshot"); webshot::install_phantomjs()}.
If the file does not end in a \code{.png}, \code{.jpg}, \code{.jpeg}, or \code{.pdf}
extension, an HTML file will be saved, and you can open the
interactive graph using a web browser.}

\item{selfcontained}{Logical, whether
to save the \code{file} as a self-contained
HTML file (with external resources base64 encoded) or a file with
external resources placed in an adjacent directory. If \code{TRUE},
pandoc is required.}

\item{...}{Arguments passed to \code{networkD3::sankeyNetwork()}.}
}
\value{
A \code{visNetwork} graph.
}
\description{
This function is called inside
\code{\link[=sankey_drake_graph]{sankey_drake_graph()}}, which typical users
call more often. A legend is unfortunately unavailable
for the graph itself, but you can see what all the colors mean with
\code{visNetwork::visNetwork(drake::legend_nodes())}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
load_mtcars_example() # Get the code with drake_example("mtcars").
if (suppressWarnings(require("knitr"))) {
if (requireNamespace("networkD3", quietly = TRUE)) {
if (requireNamespace("visNetwork", quietly = TRUE)) {
# Instead of jumpting right to sankey_drake_graph(), get the data frames
# of nodes, edges, and legend nodes.
sankey_drake_graph(my_plan) # Jump straight to the interactive graph.
# Show the legend separately.
visNetwork::visNetwork(nodes = drake::legend_nodes())
# Get the node and edge info that sankey_drake_graph() just plotted:
graph <- drake_graph_info(my_plan)
# You can pass the data frames right to render_sankey_drake_graph()
# (as in sankey_drake_graph()) or you can create
# your own custom visNewtork graph.
render_sankey_drake_graph(graph)
}
}
}
})
}
}
\seealso{
\code{\link[=sankey_drake_graph]{sankey_drake_graph()}}, \code{\link[=vis_drake_graph]{vis_drake_graph()}},
\code{\link[=drake_ggraph]{drake_ggraph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{summaries}
\alias{summaries}
\title{summaries \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
summaries(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{default_Makefile_args}
\alias{default_Makefile_args}
\title{Default arguments of Makefile parallelism
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
default_Makefile_args(jobs, verbose)
}
\arguments{
\item{jobs}{Number of jobs.}

\item{verbose}{Integer, control printing to the console/terminal.
\itemize{
\item \code{0}: print nothing.
\item \code{1}: print target-by-target messages as \code{\link[=make]{make()}} progresses.
\item \code{2}: show a progress bar to track how many targets are
done so far.
}}
}
\value{
\code{args} for \code{system2(command, args)}
}
\description{
2019-01-03
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{predict_load_balancing}
\alias{predict_load_balancing}
\title{Predict parallel computing behavior
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
predict_load_balancing(
  config,
  targets = NULL,
  from_scratch = FALSE,
  targets_only = NULL,
  jobs = 1,
  known_times = numeric(0),
  default_time = 0,
  warn = TRUE
)
}
\arguments{
\item{config}{Deprecated.}

\item{from_scratch}{Logical, whether to predict a
\code{\link[=make]{make()}} build from scratch or to
take into account the fact that some targets may be
already up to date and therefore skipped.}

\item{targets_only}{Deprecated.}

\item{known_times}{A named numeric vector with targets/imports
as names and values as hypothetical runtimes in seconds.
Use this argument to overwrite any of the existing build times
or the \code{default_time}.}

\item{default_time}{Number of seconds to assume for any
target or import with no recorded runtime (from \code{\link[=build_times]{build_times()}})
or anything in \code{known_times}.}

\item{warn}{Logical, whether to warn the user about
any targets with no available runtime, either in
\code{known_times} or \code{\link[=build_times]{build_times()}}. The times for these
targets default to \code{default_time}.}
}
\value{
A data frame showing one likely arrangement
of targets assigned to parallel workers.
}
\description{
Deprecated on 2019-02-14.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.R
\name{do_prework}
\alias{do_prework}
\title{Do the prework in the \code{prework}
argument to \code{\link[=make]{make()}}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
do_prework(config, verbose_packages)
}
\arguments{
\item{config}{A configured workflow from \code{\link[=drake_config]{drake_config()}}.}

\item{verbose_packages}{logical, whether to print
package startup messages}
}
\value{
Inivisibly returns \code{NULL}.
}
\description{
For internal use only.
The only reason this function is exported
is to set up parallel socket (PSOCK) clusters
without too much fuss.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
# Create a main internal configuration list with prework.
con <- drake_config(my_plan, prework = c("library(knitr)", "x <- 1"))
# Do the prework. Usually done at the beginning of `make()`,
# and for distributed computing backends like "future_lapply",
# right before each target is built.
do_prework(config = con, verbose_packages = TRUE)
# The `eval` element is the environment where the prework
# and the commands in your workflow plan data frame are executed.
identical(con$eval$x, 1) # Should be TRUE.
}
})
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{drake_envir}
\alias{drake_envir}
\title{Get the environment where drake builds targets
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#questioning}{\figure{lifecycle-questioning.svg}{options: alt='[Questioning]'}}}{\strong{[Questioning]}}}
\usage{
drake_envir(which = c("targets", "dynamic", "subtargets", "imports"))
}
\arguments{
\item{which}{Character of length 1, which environment
to select. See the details of this help file.}
}
\value{
The environment where \code{drake} builds targets.
}
\description{
Call this function inside the commands in your plan
to get the environment where \code{drake} builds targets.
Advanced users can use it to strategically remove targets from memory
while \code{\link[=make]{make()}} is running.
}
\details{
\code{drake} manages in-memory targets in 4 environments:
one with sub-targets, one with whole dynamic targets, one with
static targets, and one with imported global objects and functions.
This last environment is usually the environment
from which you call \code{\link[=make]{make()}}.
Select the appropriate environment for your
use case with the \code{which} argument of \code{drake_envir()}.
}
\section{Keywords}{

\code{\link[=drake_plan]{drake_plan()}} understands special keyword functions for your commands.
With the exception of \code{\link[=target]{target()}}, each one is a proper function
with its own help file.
\itemize{
\item \code{\link[=target]{target()}}: give the target more than just a command.
Using \code{\link[=target]{target()}}, you can apply a transformation
(examples: \url{https://books.ropensci.org/drake/plans.html#large-plans}), # nolint
supply a trigger (\url{https://books.ropensci.org/drake/triggers.html}), # nolint
or set any number of custom columns.
\item \code{\link[=file_in]{file_in()}}: declare an input file dependency.
\item \code{\link[=file_out]{file_out()}}: declare an output file to be produced
when the target is built.
\item \code{\link[=knitr_in]{knitr_in()}}: declare a \code{knitr} file dependency such as an
R Markdown (\verb{*.Rmd}) or R LaTeX (\verb{*.Rnw}) file.
\item \code{\link[=ignore]{ignore()}}: force \code{drake} to entirely ignore a piece of code:
do not track it for changes and do not analyze it for dependencies.
\item \code{\link[=no_deps]{no_deps()}}: tell \code{drake} to not track the dependencies
of a piece of code. \code{drake} still tracks the code itself for changes.
\item \code{\link[=id_chr]{id_chr()}}: Get the name of the current target.
\item \code{\link[=drake_envir]{drake_envir()}}: get the environment where drake builds targets.
Intended for advanced custom memory management.
}
}

\examples{
\dontrun{
isolate_example("contain side effects", {
plan <- drake_plan(
  large_data_1 = sample.int(1e4),
  large_data_2 = sample.int(1e4),
  subset = c(large_data_1[seq_len(10)], large_data_2[seq_len(10)]),
  summary = {
    print(ls(envir = parent.env(drake_envir())))
    # We don't need the large_data_* targets in memory anymore.
    rm(large_data_1, large_data_2, envir = drake_envir("targets"))
    print(ls(envir = drake_envir("targets")))
    mean(subset)
  }
)
make(plan, cache = storr::storr_environment(), session_info = FALSE)
})
}
}
\seealso{
\code{\link[=from_plan]{from_plan()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deps.R
\name{deps_target_impl}
\alias{deps_target_impl}
\title{Internal function with a drake_config() argument}
\usage{
deps_target_impl(target, config, character_only = FALSE)
}
\arguments{
\item{target}{Name of a target.}

\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}

\item{character_only}{Logical, whether to interpret
\code{target} as a character (\code{TRUE}) or a symbol (\code{FALSE}).}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deps.R
\name{deps_knitr}
\alias{deps_knitr}
\title{Find the drake dependencies of a dynamic knitr report target.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
deps_knitr(path)
}
\arguments{
\item{path}{Encoded file path to the \code{knitr}/R Markdown document.
Wrap paths in \code{\link[=file_store]{file_store()}} to encode.}
}
\value{
A data frame of dependencies.
}
\description{
Dependencies in \code{knitr} reports are marked
by \code{\link[=loadd]{loadd()}} and \code{\link[=readd]{readd()}} in active code chunks.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
load_mtcars_example() # Get the code with drake_example("mtcars").
deps_knitr("report.Rmd")
})
}
}
\seealso{
\code{\link[=deps_code]{deps_code()}}, \code{\link[=deps_target]{deps_target()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{plan_drake}
\alias{plan_drake}
\title{plan_drake \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
plan_drake(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{read_graph}
\alias{read_graph}
\title{read_graph \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
read_graph(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{triggers}
\alias{triggers}
\title{List the old drake triggers.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
triggers()
}
\value{
A character vector with the names of the old triggers.
}
\description{
Triggers are target-level rules
that tell \code{\link[=make]{make()}} how to know if a target
is outdated or up to date.
}
\details{
Deprecated on 2018-07-22.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text_drake_graph.R
\name{text_drake_graph_impl}
\alias{text_drake_graph_impl}
\title{Internal function with a drake_config() argument}
\usage{
text_drake_graph_impl(
  config,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  targets_only = FALSE,
  make_imports = TRUE,
  from_scratch = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE,
  nchar = 1L,
  print = TRUE
)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to save some time and risk obsolete output.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio.R
\name{rs_addin_r_vis_drake_graph}
\alias{rs_addin_r_vis_drake_graph}
\title{RStudio addin for r_vis_drake_graph()
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
rs_addin_r_vis_drake_graph(r_args = list(), .print = TRUE)
}
\arguments{
\item{r_args}{List of arguments to \code{r_fn}, not including \code{func} or \code{args}.
Example:
\code{r_make(r_fn = callr::r_bg, r_args = list(stdout = "stdout.log"))}.}

\item{.print}{Logical, whether to \code{print()} the result
to the console. Required for the addin.}
}
\value{
A \code{visNetwork} graph.
}
\description{
Call \code{\link[=r_vis_drake_graph]{r_vis_drake_graph()}} in an RStudio addin.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{build_drake_graph}
\alias{build_drake_graph}
\title{Function \code{build_drake_graph}
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
build_drake_graph(
  plan,
  targets = plan$target,
  envir = parent.frame(),
  verbose = 1L,
  jobs = 1,
  console_log_file = NULL,
  trigger = drake::trigger(),
  cache = NULL
)
}
\arguments{
\item{plan}{Workflow plan data frame.
A workflow plan data frame is a data frame
with a \code{target} column and a \code{command} column.
(See the details in the \code{\link[=drake_plan]{drake_plan()}} help file
for descriptions of the optional columns.)
Targets are the objects that drake generates,
and commands are the pieces of R code that produce them.
You can create and track custom files along the way
(see \code{\link[=file_in]{file_in()}}, \code{\link[=file_out]{file_out()}}, and \code{\link[=knitr_in]{knitr_in()}}).
Use the function \code{\link[=drake_plan]{drake_plan()}} to generate workflow plan
data frames.}

\item{targets}{Character vector, names of targets to build.
Dependencies are built too. You may supply static and/or whole
dynamic targets, but no sub-targets.}

\item{envir}{Environment to use. Defaults to the current
workspace, so you should not need to worry about this
most of the time. A deep copy of \code{envir} is made,
so you don't need to worry about your workspace being modified
by \code{make}. The deep copy inherits from the global environment.
Wherever necessary, objects and functions are imported
from \code{envir} and the global environment and
then reproducibly tracked as dependencies.}

\item{verbose}{Integer, control printing to the console/terminal.
\itemize{
\item \code{0}: print nothing.
\item \code{1}: print target-by-target messages as \code{\link[=make]{make()}} progresses.
\item \code{2}: show a progress bar to track how many targets are
done so far.
}}

\item{jobs}{Maximum number of parallel workers for processing the targets.
You can experiment with \code{\link[=predict_runtime]{predict_runtime()}}
to help decide on an appropriate number of jobs.
For details, visit
\url{https://books.ropensci.org/drake/time.html}.}

\item{console_log_file}{Deprecated in favor of \code{log_make}.}

\item{trigger}{Name of the trigger to apply to all targets.
Ignored if \code{plan} has a \code{trigger} column.
See \code{\link[=trigger]{trigger()}} for details.}

\item{cache}{drake cache as created by \code{\link[=new_cache]{new_cache()}}.
See also \code{\link[=drake_cache]{drake_cache()}}.}
}
\value{
An \code{igraph} object.
}
\description{
Use \code{\link[=drake_config]{drake_config()}} instead.
}
\details{
Deprecated on 2018-11-02.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sankey_drake_graph.R
\name{sankey_drake_graph}
\alias{sankey_drake_graph}
\title{Show a Sankey graph of your drake project.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
sankey_drake_graph(
  ...,
  file = character(0),
  selfcontained = FALSE,
  build_times = "build",
  digits = 3,
  targets_only = FALSE,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  make_imports = TRUE,
  from_scratch = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE,
  config = NULL
)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{file}{Name of a file to save the graph.
If \code{NULL} or \code{character(0)}, no file is saved and
the graph is rendered and displayed within R.
If the file ends in a \code{.png}, \code{.jpg}, \code{.jpeg}, or \code{.pdf} extension,
then a static image will be saved. In this case,
the webshot package and PhantomJS are required:
\verb{install.packages("webshot"); webshot::install_phantomjs()}.
If the file does not end in a \code{.png}, \code{.jpg}, \code{.jpeg}, or \code{.pdf}
extension, an HTML file will be saved, and you can open the
interactive graph using a web browser.}

\item{selfcontained}{Logical, whether
to save the \code{file} as a self-contained
HTML file (with external resources base64 encoded) or a file with
external resources placed in an adjacent directory. If \code{TRUE},
pandoc is required.}

\item{build_times}{Character string or logical.
If character, the choices are
1. \code{"build"}: runtime of the command plus the time
it take to store the target or import.
2. \code{"command"}: just the runtime of the command.
3. \code{"none"}: no build times.
If logical, \code{build_times} selects whether to show the
times from `build_times(..., type = "build")`` or use
no build times at all. See \code{\link[=build_times]{build_times()}} for details.}

\item{digits}{Number of digits for rounding the build times}

\item{targets_only}{Logical,
whether to skip the imports and only include the
targets in the workflow plan.}

\item{from}{Optional collection of target/import names.
If \code{from} is nonempty,
the graph will restrict itself to
a neighborhood of \code{from}.
Control the neighborhood with
\code{mode} and \code{order}.}

\item{mode}{Which direction to branch out in the graph
to create a neighborhood around \code{from}.
Use \code{"in"} to go upstream,
\code{"out"} to go downstream,
and \code{"all"} to go both ways and disregard
edge direction altogether.}

\item{order}{How far to branch out to create
a neighborhood around \code{from}. Defaults to
as far as possible. If a target is in the neighborhood, then
so are all of its custom \code{\link[=file_out]{file_out()}} files if
\code{show_output_files} is \code{TRUE}.
That means the actual graph order may be slightly greater than
you might expect, but this ensures consistency
between \code{show_output_files = TRUE} and
\code{show_output_files = FALSE}.}

\item{subset}{Optional character vector.
Subset of targets/imports to display in the graph.
Applied after \code{from}, \code{mode}, and \code{order}.
Be advised: edges are only kept for adjacent nodes in \code{subset}.
If you do not select all the intermediate nodes,
edges will drop from the graph.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to increase speed and risk using obsolete information.}

\item{from_scratch}{Logical, whether to assume all the targets
will be made from scratch on the next \code{\link[=make]{make()}}.
Makes all targets outdated, but keeps information about
build progress in previous \code{\link[=make]{make()}}s.}

\item{group}{Optional character scalar, name of the column used to
group nodes into columns. All the columns names of your original \code{drake}
plan are choices. The other choices (such as \code{"status"}) are column names
in the \code{nodes} . To group nodes into clusters in the graph,
you must also supply the \code{clusters} argument.}

\item{clusters}{Optional character vector of values to cluster on.
These values must be elements of the column of the \code{nodes} data frame
that you specify in the \code{group} argument to \code{drake_graph_info()}.}

\item{show_output_files}{Logical, whether to include
\code{\link[=file_out]{file_out()}} files in the graph.}

\item{config}{Deprecated.}
}
\value{
A \code{visNetwork} graph.
}
\description{
To save time for repeated plotting,
this function is divided into
\code{\link[=drake_graph_info]{drake_graph_info()}} and \code{\link[=render_sankey_drake_graph]{render_sankey_drake_graph()}}.
A legend is unfortunately unavailable
for the graph itself, but you can see what all the colors mean with
\code{visNetwork::visNetwork(drake::legend_nodes())}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
if (requireNamespace("networkD3", quietly = TRUE)) {
if (requireNamespace("visNetwork", quietly = TRUE)) {
# Plot the network graph representation of the workflow.
sankey_drake_graph(my_plan)
# Show the legend separately.
visNetwork::visNetwork(nodes = drake::legend_nodes())
make(my_plan) # Run the project, build the targets.
sankey_drake_graph(my_plan) # The black nodes from before are now green.
# Plot a subgraph of the workflow.
sankey_drake_graph(my_plan, from = c("small", "reg2"))
}
}
}
})
}
}
\seealso{
\code{\link[=render_sankey_drake_graph]{render_sankey_drake_graph()}}, \code{\link[=vis_drake_graph]{vis_drake_graph()}},
\code{\link[=drake_ggraph]{drake_ggraph()}}, \code{\link[=text_drake_graph]{text_drake_graph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{parallelism_choices}
\alias{parallelism_choices}
\title{Names of old parallel backends
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
parallelism_choices(distributed_only = FALSE)
}
\arguments{
\item{distributed_only}{Logical.}
}
\value{
character vector
}
\description{
2019-01-03
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{drake_meta}
\alias{drake_meta}
\title{Compute the initial pre-build metadata
of a target or import.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
drake_meta(target, config)
}
\arguments{
\item{target}{Character scalar, name of the target
to get metadata.}

\item{config}{Top-level internal configuration list produced
by \code{\link[=drake_config]{drake_config()}}.}
}
\value{
A list of metadata on a target. Does not include
the file modification time if the target is a file.
That piece is computed later in \code{\link[=make]{make()}} by
\code{drake:::store_outputs()}.
}
\description{
Deprecated on 2019-01-12.
}
\details{
The metadata helps determine if the
target is up to date or outdated. The metadata of imports
is used to compute the metadata of targets.
Target metadata is computed with \code{drake_meta()}, and then
\code{drake:::store_outputs()} completes the metadata
after the target is built.
In other words, the output of \code{drake_meta()} corresponds
to the state of the target immediately before \code{\link[=make]{make()}}
builds it.
See \code{\link[=diagnose]{diagnose()}} to read the final metadata of a target,
including any errors, warnings, and messages in the last build.
}
\seealso{
\code{\link[=diagnose]{diagnose()}}, \code{\link[=deps_profile]{deps_profile()}}, \code{\link[=make]{make()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_drake_graph.R
\name{vis_drake_graph_impl}
\alias{vis_drake_graph_impl}
\title{Internal function with a drake_config() argument}
\usage{
vis_drake_graph_impl(
  config,
  file = character(0),
  selfcontained = FALSE,
  build_times = "build",
  digits = 3,
  targets_only = FALSE,
  font_size = 20,
  layout = NULL,
  main = NULL,
  direction = NULL,
  hover = FALSE,
  navigationButtons = TRUE,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  ncol_legend = 1,
  full_legend = FALSE,
  make_imports = TRUE,
  from_scratch = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE,
  collapse = TRUE,
  on_select_col = NULL,
  on_select = NULL,
  level_separation = NULL
)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_runtime.R
\name{predict_runtime_impl}
\alias{predict_runtime_impl}
\title{Internal function with a drake_config() argument}
\usage{
predict_runtime_impl(
  config,
  targets_predict = NULL,
  from_scratch = FALSE,
  targets_only = NULL,
  jobs_predict = 1L,
  known_times = numeric(0),
  default_time = 0,
  warn = TRUE
)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_runtime.R
\name{predict_workers_impl}
\alias{predict_workers_impl}
\title{Internal function with a drake_config() argument}
\usage{
predict_workers_impl(
  config,
  targets_predict = NULL,
  from_scratch = FALSE,
  targets_only = NULL,
  jobs_predict = 1,
  known_times = numeric(0),
  default_time = 0,
  warn = TRUE
)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_graph_info.R
\name{drake_graph_info_impl}
\alias{drake_graph_info_impl}
\title{Internal function}
\usage{
drake_graph_info_impl(
  config,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  build_times = "build",
  digits = 3,
  targets_only = FALSE,
  font_size = 20,
  from_scratch = FALSE,
  make_imports = TRUE,
  full_legend = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE,
  hover = FALSE,
  on_select_col = NULL
)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{drake_batchtools_tmpl_file}
\alias{drake_batchtools_tmpl_file}
\title{Get a template file for execution on a cluster.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
drake_batchtools_tmpl_file(
  example = drake::drake_hpc_template_files(),
  to = getwd(),
  overwrite = FALSE
)
}
\arguments{
\item{example}{Name of template file.}

\item{to}{Character vector, where to write the file.}

\item{overwrite}{Logical, whether to overwrite an existing file of the
same name.}
}
\description{
Deprecated. Use \code{\link[=drake_hpc_template_file]{drake_hpc_template_file()}} instead.
}
\details{
Deprecated on 2018-06-27.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{drake_cancelled}
\alias{drake_cancelled}
\title{List cancelled targets.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_cancelled(cache = drake::drake_cache(path = path), path = NULL)
}
\arguments{
\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}
}
\value{
A character vector of target names.
}
\description{
List the targets that were cancelled in the current or
previous call to \code{\link[=make]{make()}} using \code{\link[=cancel]{cancel()}} or \code{\link[=cancel_if]{cancel_if()}}.
}
\examples{
\dontrun{
isolate_example("contain side effects", {
plan <- drake_plan(x = 1, y = cancel_if(x > 0))
make(plan)
drake_cancelled()
})
}
}
\seealso{
\code{\link[=drake_running]{drake_running()}}, \code{\link[=drake_failed]{drake_failed()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{build_graph}
\alias{build_graph}
\title{build_graph \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
build_graph(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{max_useful_jobs}
\alias{max_useful_jobs}
\title{max_useful_jobs \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
max_useful_jobs(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-05-16
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{code_to_plan}
\alias{code_to_plan}
\title{Turn an R script file or \code{knitr} / R Markdown report
into a \code{drake} plan.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#questioning}{\figure{lifecycle-questioning.svg}{options: alt='[Questioning]'}}}{\strong{[Questioning]}}}
\usage{
code_to_plan(path)
}
\arguments{
\item{path}{A file path to an R script or \code{knitr} report.}
}
\description{
\code{code_to_plan()}, \code{\link[=plan_to_code]{plan_to_code()}}, and
\code{\link[=plan_to_notebook]{plan_to_notebook()}} together illustrate the relationships
between \code{drake} plans, R scripts, and R Markdown documents.
}
\details{
This feature is easy to break, so there are some rules
for your code file:
\enumerate{
\item Stick to assigning a single expression to a single target at a time.
For multi-line commands, please enclose the whole command
in curly braces.
Conversely, compound assignment is not supported
(e.g. \code{target_1 <- target_2 <- target_3 <- get_data()}).
\item Once you assign an expression to a variable,
do not modify the variable any more.
The target/command binding should be permanent.
\item Keep it simple. Please use the assignment operators rather than
\code{assign()} and similar functions.
}
}
\examples{
plan <- drake_plan(
  raw_data = read_excel(file_in("raw_data.xlsx")),
  data = raw_data,
  hist = create_plot(data),
  fit = lm(Ozone ~ Temp + Wind, data)
)
file <- tempfile()
# Turn the plan into an R script a the given file path.
plan_to_code(plan, file)
# Here is what the script looks like.
cat(readLines(file), sep = "\n")
# Convert back to a drake plan.
code_to_plan(file)
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}, \code{\link[=plan_to_code]{plan_to_code()}},
\code{\link[=plan_to_notebook]{plan_to_notebook()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{prune_drake_graph}
\alias{prune_drake_graph}
\title{Prune the graph
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
prune_drake_graph(graph, to = igraph::V(graph)$name, jobs = 1)
}
\arguments{
\item{graph}{An igraph object.}

\item{to}{Character vector of vertices.}

\item{jobs}{Number of jobs for parallelism.}
}
\value{
An \code{igraph} object
}
\description{
2019-01-08
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{plan_summaries}
\alias{plan_summaries}
\title{Specialized wildcard for summaries
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
plan_summaries(
  plan,
  analyses,
  datasets,
  gather = rep("list", nrow(plan)),
  sep = "_"
)
}
\arguments{
\item{plan}{Workflow plan data frame with commands for the summaries.
Use the \code{analysis__} and \code{dataset__} wildcards
just like the \code{dataset__} wildcard in \code{\link[=plan_analyses]{plan_analyses()}}.}

\item{analyses}{Workflow plan data frame of analysis instructions.}

\item{datasets}{Workflow plan data frame with instructions to make
or import the datasets.}

\item{gather}{Character vector, names of functions to gather the
summaries. If not \code{NULL}, the length must be the number of
rows in the \code{plan}. See the \code{\link[=gather_plan]{gather_plan()}} function
for more.}

\item{sep}{Character scalar, delimiter for creating the
new target names.}
}
\value{
An evaluated workflow plan data frame of instructions
for computing summaries of analyses and datasets.
analyses of multiple datasets in multiple ways.
}
\description{
Use \code{\link[=drake_plan]{drake_plan()}} with transformations instead. See
\url{https://books.ropensci.org/drake/plans.html#large-plans}
for details.
}
\details{
2019-01-13
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_imports.R
\name{process_import}
\alias{process_import}
\title{Process an imported data object
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
process_import(import, config)
}
\arguments{
\item{import}{Character, name of an import to process}

\item{config}{\code{\link[=drake_config]{drake_config()}} object}
}
\description{
For internal use only. Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{default_parallelism}
\alias{default_parallelism}
\title{Default parallel backend
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
default_parallelism()
}
\value{
character
}
\description{
2019-01-02
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{drake_cache}
\alias{drake_cache}
\title{Get the cache of a \code{drake} project.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_cache(path = NULL, verbose = NULL, console_log_file = NULL)
}
\arguments{
\item{path}{Character.
Set \code{path} to the path of a \code{storr::storr_rds()} cache
to retrieve a specific cache generated by \code{storr::storr_rds()}
or \code{drake::new_cache()}. If the \code{path} argument is \code{NULL},
\code{drake_cache()} searches up through parent directories
to find a folder called \verb{.drake/}.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{console_log_file}{Deprecated on 2019-09-11.}
}
\value{
A drake/storr cache in a folder called \verb{.drake/},
if available. \code{NULL} otherwise.
}
\description{
\code{\link[=make]{make()}} saves the values of your targets so
you rarely need to think about output files. By default,
the cache is a hidden folder called \verb{.drake/}.
You can also supply your own \code{storr} cache to the \code{cache}
argument of \code{make()}. The \code{drake_cache()} function retrieves
this cache.
}
\details{
\code{drake_cache()} actually returns a \emph{decorated} \code{storr},
an object that \emph{contains} a \code{storr} (plus bells and whistles).
To get the \emph{actual} inner \code{storr}, use \code{drake_cache()$storr}.
Most methods are delegated to the inner \code{storr}.
Some methods and objects are new or overwritten. Here
are the ones relevant to users.
\itemize{
\item \code{history}: \code{drake}'s history (which powers \code{\link[=drake_history]{drake_history()}})
is a \href{https://github.com/wlandau/txtq}{\code{txtq}}. Access it
with \code{drake_cache()$history}.
\item \code{import()}: The \code{import()} method is a function that can import
targets, function dependencies, etc. from one decorated \code{storr}
to another. History is not imported. For that, you have to work
with the history \code{txtq}s themselves, Arguments to \code{import()}:
\itemize{
\item \code{...} and \code{list}: specify targets to import just like with \code{\link[=loadd]{loadd()}}.
Leave these blank to import everything.
\item \code{from}: the decorated \code{storr} from which to import targets.
\item \code{jobs}: number of local processes for parallel computing.
\item \code{gc}: \code{TRUE} or \code{FALSE}, whether to run garbage collection for memory
after importing each target. Recommended, but slow.
}
\item \code{export()}: Same as \code{import()}, except the \code{from} argument is replaced
by \code{to}: the decorated \code{storr} where the targets end up.
}
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
clean(destroy = TRUE)
# No cache is available.
drake_cache() # NULL
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the targets.
x <- drake_cache() # Now, there is a cache.
y <- storr::storr_rds(".drake") # Nearly equivalent.
# List the objects readable from the cache with readd().
x$list()
# drake_cache() actually returns a *decorated* storr.
# The *real* storr is inside.
drake_cache()$storr
}
# You can import and export targets to and from decorated storrs.
plan1 <- drake_plan(w = "w", x = "x")
plan2 <- drake_plan(a = "a", x = "x2")
cache1 <- new_cache("cache1")
cache2 <- new_cache("cache2")
make(plan1, cache = cache1)
make(plan2, cache = cache2)
cache1$import(cache2, a)
cache1$get("a")
cache1$get("x")
cache1$import(cache2)
cache1$get("x")
# With txtq >= 0.1.6.9002, you can import history from one cache into
# another.
# nolint start
# drake_history(cache = cache1)
# cache1$history$import(cache2$history)
# drake_history(cache = cache1)
# nolint end
})
}
}
\seealso{
\code{\link[=new_cache]{new_cache()}}, \code{\link[=drake_config]{drake_config()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deps.R
\name{deps_target}
\alias{deps_target}
\title{List the dependencies of a target
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
deps_target(target, ..., character_only = FALSE, config = NULL)
}
\arguments{
\item{target}{A symbol denoting a target name, or if \code{character_only}
is TRUE, a character scalar denoting a target name.}

\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{character_only}{Logical, whether to assume target is a character
string rather than a symbol.}

\item{config}{Deprecated.}
}
\value{
A data frame with the dependencies listed by type
(globals, files, etc).
}
\description{
Intended for debugging and checking your project.
The dependency structure of the components of your analysis
decides which targets are built and when.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
load_mtcars_example() # Get the code with drake_example("mtcars").
deps_target(regression1_small, my_plan)
})
}
}
\seealso{
\code{\link[=deps_code]{deps_code()}}, \code{\link[=deps_knitr]{deps_knitr()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{readd}
\alias{readd}
\alias{loadd}
\title{Read and return a drake target/import from the cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
readd(
  target,
  character_only = FALSE,
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  namespace = NULL,
  verbose = 1L,
  show_source = FALSE,
  subtargets = NULL,
  subtarget_list = FALSE
)

loadd(
  ...,
  list = character(0),
  imported_only = NULL,
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  namespace = NULL,
  envir = parent.frame(),
  jobs = 1,
  verbose = 1L,
  deps = FALSE,
  lazy = "eager",
  graph = NULL,
  replace = TRUE,
  show_source = FALSE,
  tidyselect = !deps,
  config = NULL,
  subtargets = NULL,
  subtarget_list = FALSE
)
}
\arguments{
\item{target}{If \code{character_only} is \code{TRUE}, then
\code{target} is a character string naming the object to read.
Otherwise, \code{target} is an unquoted symbol with the name of the
object.}

\item{character_only}{Logical, whether \code{name} should be treated
as a character or a symbol
(just like \code{character.only} in \code{\link[=library]{library()}}).}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{namespace}{Optional character string,
name of the \code{storr} namespace to read from.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{show_source}{Logical, option to show the command
that produced the target or indicate that the object
was imported (using \code{\link[=show_source]{show_source()}}).}

\item{subtargets}{A numeric vector of indices.
If \code{target} is dynamic, \code{\link[=loadd]{loadd()}} and \code{\link[=readd]{readd()}} retrieve
a list of sub-targets. You can restrict which sub-targets
to retrieve with the \code{subtargets} argument. For example,
\code{readd(x, subtargets = seq_len(3))} only retrieves the
first 3 sub-targets of dynamic target \code{x}.}

\item{subtarget_list}{Logical, for dynamic targets only.
If \code{TRUE}, the dynamic target is loaded as a named
list of sub-target values. If \code{FALSE}, \code{drake}
attempts to concatenate the sub-targets with \code{vctrs::vec_c()}
(and returns an unnamed list if such concatenation is not possible).}

\item{...}{Targets to load from the cache: as names (symbols) or
character strings. If the \code{tidyselect} package is installed,
you can also supply \code{dplyr}-style \code{tidyselect}
commands such as \code{starts_with()}, \code{ends_with()}, and \code{one_of()}.}

\item{list}{Character vector naming targets to be loaded from the
cache. Similar to the \code{list} argument of \code{\link[=remove]{remove()}}.}

\item{imported_only}{Logical, deprecated.}

\item{envir}{Environment to load objects into. Defaults to the
calling environment (current workspace).}

\item{jobs}{Number of parallel jobs for loading objects. On
non-Windows systems, the loading process for multiple objects
can be lightly parallelized via \code{parallel::mclapply()}.
just set jobs to be an integer greater than 1. On Windows,
\code{jobs} is automatically demoted to 1.}

\item{deps}{Logical, whether to load any cached
dependencies of the targets
instead of the targets themselves.

Important note:
\code{deps = TRUE} disables \code{tidyselect} functionality. For example,
\code{loadd(starts_with("model_"), config = config, deps = TRUE)}
does not work. For the selection mechanism to work,
the \verb{model_*} targets to need to already be in the cache,
which is not always the case when you are debugging your projects.
To help \code{drake} understand what you mean,
you must name the targets \emph{explicitly} when \code{deps} is \code{TRUE}, e.g.
\code{loadd(model_A, model_B, config = config, deps = TRUE)}.}

\item{lazy}{Either a string or a logical. Choices:
\itemize{
\item \code{"eager"}: no lazy loading. The target is loaded right away
with \code{\link[=assign]{assign()}}.
\item \code{"promise"}: lazy loading with \code{\link[=delayedAssign]{delayedAssign()}}
\item \code{"bind"}: lazy loading with active bindings:
\code{bindr::populate_env()}.
\item \code{TRUE}: same as \code{"promise"}.
\item \code{FALSE}: same as \code{"eager"}.
}}

\item{graph}{Deprecated.}

\item{replace}{Logical. If \code{FALSE},
items already in your environment
will not be replaced.}

\item{tidyselect}{Logical, whether to enable
\code{tidyselect} expressions in \code{...} like
\code{starts_with("prefix")} and \code{ends_with("suffix")}.}

\item{config}{Optional \code{\link[=drake_config]{drake_config()}} object.
You should supply one if \code{deps} is \code{TRUE}.}
}
\value{
The cached value of the \code{target}.
}
\description{
\code{\link[=readd]{readd()}} returns an object from the cache,
and \code{\link[=loadd]{loadd()}} loads one or more objects from the cache
into your environment or session. These objects are usually
targets built by \code{\link[=make]{make()}}. If \code{target} is dynamic,
\code{\link[=readd]{readd()}} and \code{\link[=loadd]{loadd()}} retrieve a list of sub-target values.
You can restrict which sub-targets to include using the \code{subtargets}
argument.
}
\details{
There are three uses for the
\code{\link[=loadd]{loadd()}} and \code{\link[=readd]{readd()}} functions:
\enumerate{
\item Exploring the results outside the \code{drake}/\code{make()} pipeline.
When you call \code{\link[=make]{make()}} to run your project,
\code{drake} puts the targets in a cache, usually a folder called \code{.drake}.
You may want to inspect the targets afterwards, possibly in an
interactive R session. However, the files in the \code{.drake} folder
are organized in a special format created by the
\href{https://github.com/richfitz/storr}{\code{storr}} package,
which is not exactly human-readable.
To retrieve a target for manual viewing, use \code{\link[=readd]{readd()}}.
To load one or more targets into your session, use \code{\link[=loadd]{loadd()}}.
\item In \code{knitr} / R Markdown reports.
You can borrow \code{drake} targets in your active code chunks
if you have the right calls to \code{\link[=loadd]{loadd()}} and \code{\link[=readd]{readd()}}.
These reports can either run outside the \code{drake} pipeline,
or better yet, as part of the pipeline itself.
If you call \code{knitr_in("your_report.Rmd")} inside a \code{\link[=drake_plan]{drake_plan()}}
command, then \code{\link[=make]{make()}} will scan \code{"your_report.Rmd"} for
calls to \code{loadd()} and \code{readd()} in active code chunks,
and then treat those loaded targets as dependencies.
That way, \code{\link[=make]{make()}} will automatically (re)run the report if those
dependencies change.
\item If you are using \code{make(memory_strategy = "none")}
or \code{make(memory_strategy = "unload")},
\code{\link[=loadd]{loadd()}} and \code{\link[=readd]{readd()}} can manually load dependencies
into memory for the target that is being built.
If you do this, you must carefully inspect \code{\link[=deps_target]{deps_target()}}
and \code{\link[=vis_drake_graph]{vis_drake_graph()}} before running \code{\link[=make]{make()}}
to be sure the dependency relationships among targets
are correct. If you do not wish to incur extra dependencies
with \code{\link[=loadd]{loadd()}} or \code{\link[=readd]{readd()}}, you will need to use \code{\link[=ignore]{ignore()}},
e.g. \code{drake_plan(x = 1, y = ignore(readd(x)))} or
\code{drake_plan(x = 1, y = readd(ignore("x"), character_only = TRUE))}.
Compare those plans to \code{drake_plan(x = 1, y = readd(x))}
and \code{drake_plan(x = 1, y = readd("x", character_only = TRUE))}
using \code{\link[=vis_drake_graph]{vis_drake_graph()}} and \code{\link[=deps_target]{deps_target()}}.
}
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the targets.
readd(reg1) # Return imported object 'reg1' from the cache.
readd(small) # Return targets 'small' from the cache.
readd("large", character_only = TRUE) # Return 'large' from the cache.
# For external files, only the fingerprint/hash is stored.
readd(file_store("report.md"), character_only = TRUE)
}
})
}
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the projects, build the targets.
config <- drake_config(my_plan)
loadd(small) # Load target 'small' into your workspace.
small
# For many targets, you can parallelize loadd()
# using the 'jobs' argument.
loadd(list = c("small", "large"), jobs = 2)
ls()
# Load the dependencies of the target, coef_regression2_small
loadd(coef_regression2_small, deps = TRUE, config = config)
ls()
# Load all the targets listed in the workflow plan
# of the previous `make()`.
# If you do not supply any target names, `loadd()` loads all the targets.
# Be sure your computer has enough memory.
loadd()
ls()
}
})
}
}
\seealso{
\code{\link[=cached]{cached()}}, \code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}

\code{\link[=cached]{cached()}}, \code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{read_drake_meta}
\alias{read_drake_meta}
\title{read_drake_meta \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
read_drake_meta(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{default_long_hash_algo}
\alias{default_long_hash_algo}
\title{Return the default long hash algorithm for \code{make()}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
default_long_hash_algo(cache = NULL)
}
\arguments{
\item{cache}{Optional drake cache.
When you \code{\link[=configure_cache]{configure_cache()}} without
supplying a long hash algorithm,
\code{default_long_hash_algo(cache)} is the long
hash algorithm that drake picks for you.}
}
\value{
A character vector naming a hash algorithm.
}
\description{
Deprecated. drake now only uses one hash algorithm per cache.
}
\details{
Deprecated on 2018-12-12
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{drake_cache_log_file}
\alias{drake_cache_log_file}
\title{Generate a flat text log file to represent the state of
the cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
drake_cache_log_file(
  file = "drake_cache.log",
  path = getwd(),
  search = TRUE,
  cache = drake::get_cache(path = path, search = search, verbose = verbose),
  verbose = 1L,
  jobs = 1L,
  targets_only = FALSE
)
}
\arguments{
\item{file}{character scalar, name of the flat text log file.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{jobs}{Number of jobs/workers for parallel processing.}

\item{targets_only}{Logical, whether to output information only on the
targets in your workflow plan data frame. If \code{targets_only} is \code{FALSE}, the
output will include the hashes of both targets and imports.}
}
\value{
There is no return value, but a log file is generated.
}
\description{
Deprecated on 2019-03-09.
}
\details{
Calling this function to create a log file and later calling
\code{make()} makes the log file out of date. Therefore, we recommend using
\code{make()} with the \code{cache_log_file} argument to create the cache log. This
way ensures that the log is always up to date with \code{make()} results.
}
\seealso{
\code{\link[=drake_cache_log]{drake_cache_log()}}, \code{\link[=make]{make()}}, \code{\link[=get_cache]{get_cache()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_drake_graph.R
\name{vis_drake_graph}
\alias{vis_drake_graph}
\title{Show an interactive visual network representation
of your drake project.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
vis_drake_graph(
  ...,
  file = character(0),
  selfcontained = FALSE,
  build_times = "build",
  digits = 3,
  targets_only = FALSE,
  font_size = 20,
  layout = NULL,
  main = NULL,
  direction = NULL,
  hover = FALSE,
  navigationButtons = TRUE,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  ncol_legend = 1,
  full_legend = FALSE,
  make_imports = TRUE,
  from_scratch = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE,
  collapse = TRUE,
  on_select_col = NULL,
  on_select = NULL,
  level_separation = NULL,
  config = NULL
)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{file}{Name of a file to save the graph.
If \code{NULL} or \code{character(0)}, no file is saved and
the graph is rendered and displayed within R.
If the file ends in a \code{.png}, \code{.jpg}, \code{.jpeg}, or \code{.pdf} extension,
then a static image will be saved. In this case,
the webshot package and PhantomJS are required:
\verb{install.packages("webshot"); webshot::install_phantomjs()}.
If the file does not end in a \code{.png}, \code{.jpg}, \code{.jpeg}, or \code{.pdf}
extension, an HTML file will be saved, and you can open the
interactive graph using a web browser.}

\item{selfcontained}{Logical, whether
to save the \code{file} as a self-contained
HTML file (with external resources base64 encoded) or a file with
external resources placed in an adjacent directory. If \code{TRUE},
pandoc is required. The \code{selfcontained} argument only applies
to HTML files. In other words, if \code{file} is a
PNG, PDF, or JPEG file, for instance,
the point is moot.}

\item{build_times}{Character string or logical.
If character, the choices are
1. \code{"build"}: runtime of the command plus the time
it take to store the target or import.
2. \code{"command"}: just the runtime of the command.
3. \code{"none"}: no build times.
If logical, \code{build_times} selects whether to show the
times from `build_times(..., type = "build")`` or use
no build times at all. See \code{\link[=build_times]{build_times()}} for details.}

\item{digits}{Number of digits for rounding the build times}

\item{targets_only}{Logical,
whether to skip the imports and only include the
targets in the workflow plan.}

\item{font_size}{Numeric, font size of the node labels in the graph}

\item{layout}{Deprecated.}

\item{main}{Character string, title of the graph.}

\item{direction}{Deprecated.}

\item{hover}{Logical, whether to show text (file contents,
commands, etc.) when you hover your cursor over a node.}

\item{navigationButtons}{Logical, whether to add navigation buttons with
\code{visNetwork::visInteraction(navigationButtons = TRUE)}}

\item{from}{Optional collection of target/import names.
If \code{from} is nonempty,
the graph will restrict itself to
a neighborhood of \code{from}.
Control the neighborhood with
\code{mode} and \code{order}.}

\item{mode}{Which direction to branch out in the graph
to create a neighborhood around \code{from}.
Use \code{"in"} to go upstream,
\code{"out"} to go downstream,
and \code{"all"} to go both ways and disregard
edge direction altogether.}

\item{order}{How far to branch out to create
a neighborhood around \code{from}. Defaults to
as far as possible. If a target is in the neighborhood, then
so are all of its custom \code{\link[=file_out]{file_out()}} files if
\code{show_output_files} is \code{TRUE}.
That means the actual graph order may be slightly greater than
you might expect, but this ensures consistency
between \code{show_output_files = TRUE} and
\code{show_output_files = FALSE}.}

\item{subset}{Optional character vector.
Subset of targets/imports to display in the graph.
Applied after \code{from}, \code{mode}, and \code{order}.
Be advised: edges are only kept for adjacent nodes in \code{subset}.
If you do not select all the intermediate nodes,
edges will drop from the graph.}

\item{ncol_legend}{Number of columns in the legend nodes.
To remove the legend entirely, set \code{ncol_legend} to \code{NULL} or \code{0}.}

\item{full_legend}{Logical. If \code{TRUE}, all the node types
are printed in the legend. If \code{FALSE}, only the
node types used are printed in the legend.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to increase speed and risk using obsolete information.}

\item{from_scratch}{Logical, whether to assume all the targets
will be made from scratch on the next \code{\link[=make]{make()}}.
Makes all targets outdated, but keeps information about
build progress in previous \code{\link[=make]{make()}}s.}

\item{group}{Optional character scalar, name of the column used to
group nodes into columns. All the columns names of your original \code{drake}
plan are choices. The other choices (such as \code{"status"}) are column names
in the \code{nodes} . To group nodes into clusters in the graph,
you must also supply the \code{clusters} argument.}

\item{clusters}{Optional character vector of values to cluster on.
These values must be elements of the column of the \code{nodes} data frame
that you specify in the \code{group} argument to \code{drake_graph_info()}.}

\item{show_output_files}{Logical, whether to include
\code{\link[=file_out]{file_out()}} files in the graph.}

\item{collapse}{Logical, whether to allow nodes to collapse
if you double click on them.
Analogous to \code{visNetwork::visOptions(collapse = TRUE)}.}

\item{on_select_col}{Optional string corresponding to the column name
in the plan that should provide data for the \code{on_select} event.}

\item{on_select}{defines node selection event handling.
Either a string of valid JavaScript that may be passed to
\code{visNetwork::visEvents()}, or one of the following:
\code{TRUE}, \code{NULL}/\code{FALSE}. If \code{TRUE} , enables the default behavior of
opening the link specified by the \code{on_select_col} given to
\code{drake_graph_info()}. \code{NULL}/\code{FALSE} disables the behavior.}

\item{level_separation}{Numeric, \code{levelSeparation} argument to
\code{visNetwork::visHierarchicalLayout()}. Controls the distance
between hierarchical levels. Consider setting if the
aspect ratio of the graph is far from 1.
Defaults to 150 through \code{visNetwork}.}

\item{config}{Deprecated.}
}
\value{
A \code{visNetwork} graph.
}
\description{
It is good practice to visualize the dependency graph
before running the targets.
}
\details{
For enhanced interactivity in the graph, see the \code{mandrake}
package: \url{https://github.com/matthewstrasiotto/mandrake}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
# Plot the network graph representation of the workflow.
if (requireNamespace("visNetwork", quietly = TRUE)) {
vis_drake_graph(my_plan)
make(my_plan) # Run the project, build the targets.
vis_drake_graph(my_plan) # The red nodes from before are now green.
# Plot a subgraph of the workflow.
vis_drake_graph(
  my_plan,
  from = c("small", "reg2")
)
}
}
})
}
}
\seealso{
\code{\link[=render_drake_graph]{render_drake_graph()}}, \code{\link[=sankey_drake_graph]{sankey_drake_graph()}},
\code{\link[=drake_ggraph]{drake_ggraph()}}, \code{\link[=text_drake_graph]{text_drake_graph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_times.R
\name{build_times}
\alias{build_times}
\title{See the time it took to build each target.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
build_times(
  ...,
  path = NULL,
  search = NULL,
  digits = 3,
  cache = drake::drake_cache(path = path),
  targets_only = NULL,
  verbose = NULL,
  jobs = 1,
  type = c("build", "command"),
  list = character(0)
)
}
\arguments{
\item{...}{Targets to load from the cache: as names (symbols) or
character strings. If the \code{tidyselect} package is installed,
you can also supply \code{dplyr}-style \code{tidyselect}
commands such as \code{starts_with()}, \code{ends_with()}, and \code{one_of()}.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{digits}{How many digits to round the times to.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{targets_only}{Deprecated.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{jobs}{Number of jobs/workers for parallel processing.}

\item{type}{Type of time you want: either \code{"build"}
for the full build time including the time it took to
store the target, or \code{"command"} for the time it took
just to run the command.}

\item{list}{Character vector of targets to select.}
}
\value{
A data frame of times, each from \code{\link[=system.time]{system.time()}}.
}
\description{
Applies to targets in your plan, not imports or files.
}
\details{
Times for dynamic targets
(\url{https://books.ropensci.org/drake/dynamic.html})
only reflect the time it takes
to post-process the sub-targets (typically very fast)
and exclude the time it takes to build the sub-targets themselves.
Sub-targets build times are listed individually.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
if (requireNamespace("lubridate")) {
# Show the build times for the mtcars example.
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Build all the targets.
print(build_times()) # Show how long it took to build each target.
}
}
})
}
}
\seealso{
\code{\link[=predict_runtime]{predict_runtime()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_runtime.R
\name{predict_runtime}
\alias{predict_runtime}
\title{Predict the elapsed runtime of the next call to \code{make()}
for non-staged parallel backends.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
predict_runtime(
  ...,
  targets_predict = NULL,
  from_scratch = FALSE,
  targets_only = NULL,
  jobs_predict = 1L,
  known_times = numeric(0),
  default_time = 0,
  warn = TRUE,
  config = NULL
)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{targets_predict}{Character vector, names of targets
to include in the total runtime and worker predictions.}

\item{from_scratch}{Logical, whether to predict a
\code{\link[=make]{make()}} build from scratch or to
take into account the fact that some targets may be
already up to date and therefore skipped.}

\item{targets_only}{Deprecated.}

\item{jobs_predict}{The \code{jobs} argument of your next planned
\code{make()}.}

\item{known_times}{A named numeric vector with targets/imports
as names and values as hypothetical runtimes in seconds.
Use this argument to overwrite any of the existing build times
or the \code{default_time}.}

\item{default_time}{Number of seconds to assume for any
target or import with no recorded runtime (from \code{\link[=build_times]{build_times()}})
or anything in \code{known_times}.}

\item{warn}{Logical, whether to warn the user about
any targets with no available runtime, either in
\code{known_times} or \code{\link[=build_times]{build_times()}}. The times for these
targets default to \code{default_time}.}

\item{config}{Deprecated.}
}
\value{
Predicted total runtime of the next call to \code{\link[=make]{make()}}.
}
\description{
Take the past recorded runtimes times from
\code{\link[=build_times]{build_times()}} and use them to predict how the targets
will be distributed among the available workers in the
next \code{\link[=make]{make()}}. Then, predict the overall runtime to be the
runtime of the slowest (busiest) workers.
Predictions only include the time it takes to run the targets,
not overhead/preprocessing from \code{drake} itself.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the targets.
known_times <- rep(7200, nrow(my_plan))
names(known_times) <- my_plan$target
known_times
# Predict the runtime
if (requireNamespace("lubridate", quietly = TRUE)) {
predict_runtime(
  my_plan,
  jobs_predict = 7L,
  from_scratch = TRUE,
  known_times = known_times
)
predict_runtime(
  my_plan,
  jobs_predict = 8L,
  from_scratch = TRUE,
  known_times = known_times
)
balance <- predict_workers(
  my_plan,
  jobs_predict = 7L,
  from_scratch = TRUE,
  known_times = known_times
)
balance
}
}
})
}
}
\seealso{
\code{\link[=predict_workers]{predict_workers()}}, \code{\link[=build_times]{build_times()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{plan}
\alias{plan}
\title{plan \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
plan(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outdated.R
\name{missed_impl}
\alias{missed_impl}
\title{Internal function with a drake_config() argument}
\usage{
missed_impl(config)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{workflow}
\alias{workflow}
\title{workflow \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
workflow(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan.R
\name{type_sum.expr_list}
\alias{type_sum.expr_list}
\title{Type summary printing
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#questioning}{\figure{lifecycle-questioning.svg}{options: alt='[Questioning]'}}}{\strong{[Questioning]}}}
\usage{
type_sum.expr_list(x)
}
\arguments{
\item{x}{List of language objects.}
}
\description{
Ensures \verb{<expr>} is printed at the top
of any \code{drake} plan column that is a list of language objects
(e.g. \code{plan$command}).
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{evaluate_plan}
\alias{evaluate_plan}
\title{Use wildcard templating to create a
workflow plan data frame from a template data frame.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
evaluate_plan(
  plan,
  rules = NULL,
  wildcard = NULL,
  values = NULL,
  expand = TRUE,
  rename = expand,
  trace = FALSE,
  columns = "command",
  sep = "_"
)
}
\arguments{
\item{plan}{Workflow plan data frame, similar to one produced by
\code{\link[=drake_plan]{drake_plan()}}.}

\item{rules}{Named list with wildcards as names and vectors of
replacements
as values. This is a way to evaluate multiple wildcards at once.
When not \code{NULL}, \code{rules} overrules \code{wildcard} and
\code{values} if
not \code{NULL}.}

\item{wildcard}{Character scalar denoting a wildcard placeholder.}

\item{values}{Vector of values to replace the wildcard
in the drake instructions. Will be treated as a character vector.
Must be the same length as \code{plan$command} if \code{expand} is
\code{TRUE}.}

\item{expand}{If \code{TRUE}, create a new rows in the workflow plan
data frame
if multiple values are assigned to a single wildcard.
If \code{FALSE}, each occurrence of the wildcard
is replaced with the next entry in the \code{values} vector,
and the values are recycled.}

\item{rename}{Logical, whether to rename the targets
based on the values supplied for the wildcards
(based on \code{values} or \code{rules}).}

\item{trace}{Logical, whether to add columns that
trace the wildcard expansion process. These new
columns indicate which targets were evaluated and with which
wildcards.}

\item{columns}{Character vector of names of columns
to look for and evaluate the wildcards.}

\item{sep}{Character scalar, separator for the names
of the new targets generated. For example, in
\code{evaluate_plan(drake_plan(x = sqrt(y__)), list(y__ = 1:2), sep = ".")},
the names of the new targets are \code{x.1} and \code{x.2}.}
}
\value{
A workflow plan data frame with the wildcards evaluated.
}
\description{
Deprecated on 2019-05-16. Use \code{\link[=drake_plan]{drake_plan()}}
transformations instead. See
\url{https://books.ropensci.org/drake/plans.html#large-plans}
for the details.
}
\details{
The commands in workflow plan data frames can have
wildcard symbols that can stand for datasets, parameters, function
arguments, etc. These wildcards can be evaluated over a set of
possible values using \code{evaluate_plan()}.

Specify a single wildcard with the \code{wildcard}
and \code{values} arguments. In each command, the text in
\code{wildcard} will be replaced by each value in \code{values}
in turn. Specify multiple wildcards with the \code{rules} argument,
which overrules \code{wildcard} and \code{values} if
not \code{NULL}. Here, \code{rules} should be a list with wildcards
as names and vectors of possible values as list elements.
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{drake_cache_log}
\alias{drake_cache_log}
\title{Get the state of the cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_cache_log(
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  verbose = 1L,
  jobs = 1,
  targets_only = FALSE
)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{jobs}{Number of jobs/workers for parallel processing.}

\item{targets_only}{Logical, whether to output information
only on the targets in your workflow plan data frame.
If \code{targets_only} is \code{FALSE}, the output will
include the hashes of both targets and imports.}
}
\value{
Data frame of the hash keys of the targets and imports
in the cache
}
\description{
Get the fingerprints of all the targets in a data frame.
This functionality is like
\code{make(..., cache_log_file = TRUE)},
but separated and more customizable. Hopefully, this functionality
is a step toward better data versioning tools.
}
\details{
A hash is a fingerprint of an object's value.
Together, the hash keys of all your targets and imports
represent the state of your project.
Use \code{drake_cache_log()} to generate a data frame
with the hash keys of all the targets and imports
stored in your cache.
This function is particularly useful if you are
storing your drake project in a version control repository.
The cache has a lot of tiny files, so you should not put it
under version control. Instead, save the output
of \code{drake_cache_log()} as a text file after each \code{\link[=make]{make()}},
and put the text file under version control.
That way, you have a changelog of your project's results.
See the examples below for details.
Depending on your project's
history, the targets may be different than the ones
in your workflow plan data frame.
Also, the keys depend on the hash algorithm
of your cache. To define your own hash algorithm,
you can create your own \code{storr} cache and give it a hash algorithm
(e.g. \code{storr_rds(hash_algorithm = "murmur32")})
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
# Load drake's canonical example.
load_mtcars_example() # Get the code with drake_example()
# Run the project, build all the targets.
make(my_plan)
# Get a data frame of all the hash keys.
# If you want a changelog, be sure to do this after every make().
cache_log <- drake_cache_log()
head(cache_log)
# Suppress partial arg match warnings.
suppressWarnings(
  # Save the hash log as a flat text file.
  write.table(
    x = cache_log,
    file = "drake_cache.log",
    quote = FALSE,
    row.names = FALSE
  )
)
# At this point, put drake_cache.log under version control
# (e.g. with 'git add drake_cache.log') alongside your code.
# Now, every time you run your project, your commit history
# of hash_lot.txt is a changelog of the project's results.
# It shows which targets and imports changed on every commit.
# It is extremely difficult to track your results this way
# by putting the raw '.drake/' cache itself under version control.
}
})
}
}
\seealso{
\code{\link[=cached]{cached()}}, \code{\link[=drake_cache]{drake_cache()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r_make.R
\name{r_make}
\alias{r_make}
\alias{r_drake_build}
\alias{r_outdated}
\alias{r_recoverable}
\alias{r_missed}
\alias{r_deps_target}
\alias{r_drake_graph_info}
\alias{r_vis_drake_graph}
\alias{r_sankey_drake_graph}
\alias{r_drake_ggraph}
\alias{r_text_drake_graph}
\alias{r_predict_runtime}
\alias{r_predict_workers}
\title{Launch a drake function in a fresh new R process
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
r_make(source = NULL, r_fn = NULL, r_args = list())

r_drake_build(
  target,
  character_only = FALSE,
  ...,
  source = NULL,
  r_fn = NULL,
  r_args = list()
)

r_outdated(..., source = NULL, r_fn = NULL, r_args = list())

r_recoverable(..., source = NULL, r_fn = NULL, r_args = list())

r_missed(..., source = NULL, r_fn = NULL, r_args = list())

r_deps_target(
  target,
  character_only = FALSE,
  ...,
  source = NULL,
  r_fn = NULL,
  r_args = list()
)

r_drake_graph_info(..., source = NULL, r_fn = NULL, r_args = list())

r_vis_drake_graph(..., source = NULL, r_fn = NULL, r_args = list())

r_sankey_drake_graph(..., source = NULL, r_fn = NULL, r_args = list())

r_drake_ggraph(..., source = NULL, r_fn = NULL, r_args = list())

r_text_drake_graph(..., source = NULL, r_fn = NULL, r_args = list())

r_predict_runtime(..., source = NULL, r_fn = NULL, r_args = list())

r_predict_workers(..., source = NULL, r_fn = NULL, r_args = list())
}
\arguments{
\item{source}{Path to an R script file that
loads packages, functions, etc. and returns a \code{\link[=drake_config]{drake_config()}} object.
There are 3 ways to set this path.
\enumerate{
\item Pass an explicit file path.
\item Call \code{options(drake_source = "path_to_your_script.R")}.
\item Just create a file called "_drake.R" in your working directory
and supply nothing to \code{source}.
}}

\item{r_fn}{A \code{callr} function such as \code{callr::r} or \code{callr::r_bg}.
Example: \code{r_make(r_fn = callr::r)}.}

\item{r_args}{List of arguments to \code{r_fn}, not including \code{func} or \code{args}.
Example:
\code{r_make(r_fn = callr::r_bg, r_args = list(stdout = "stdout.log"))}.}

\item{target}{Name of the target.}

\item{character_only}{Logical, whether \code{name} should be treated
as a character or a symbol
(just like \code{character.only} in \code{\link[=library]{library()}}).}

\item{...}{Arguments to the inner function. For example, if you want to call
\code{\link[=r_vis_drake_graph]{r_vis_drake_graph()}}, the inner function is \code{\link[=vis_drake_graph]{vis_drake_graph()}}, and
\code{selfcontained} is an example argument you could supply to the ellipsis.}
}
\description{
The \verb{r_*()} functions, such as \code{r_make()},
enhance reproducibility by launching a \code{drake} function in
a separate R process.
}
\details{
\code{drake} searches your environment
to detect dependencies, so functions like \code{\link[=make]{make()}}, \code{\link[=outdated]{outdated()}}, etc.
are designed to run in fresh clean R sessions. Wrappers \code{\link[=r_make]{r_make()}},
\code{\link[=r_outdated]{r_outdated()}}, etc. run reproducibly even if your current R session
is old and stale.

\code{\link[=r_outdated]{r_outdated()}} runs the four steps below.
\code{\link[=r_make]{r_make()}} etc. are similar.
\enumerate{
\item Launch a new \code{callr::r()} session.
\item In that fresh session, run the R script from the \code{source} argument.
This script loads packages, functions, global options, etc.
and calls \code{\link[=drake_config]{drake_config()}} at the very end. \code{\link[=drake_config]{drake_config()}}
is the preprocessing step of \code{\link[=make]{make()}}, and it accepts
all the same arguments as \code{\link[=make]{make()}} (e.g. \code{plan} and \code{targets}).
\item In that same session, run \code{\link[=outdated]{outdated()}}
with the \code{config} argument from step 2.
\item Return the result back to main process
(e.g. your interactive R session).
}
}
\section{Recovery}{

\code{make(recover = TRUE, recoverable = TRUE)}
powers automated data recovery.
The default of \code{recover} is \code{FALSE} because
targets recovered from the distant past may have been generated
with earlier versions of R and earlier package environments
that no longer exist.

How it works: if \code{recover} is \code{TRUE},
\code{drake} tries to salvage old target values from the cache
instead of running commands from the plan.
A target is recoverable if
\enumerate{
\item There is an old value somewhere in the cache that
shares the command, dependencies, etc.
of the target about to be built.
\item The old value was generated with \code{make(recoverable = TRUE)}.
}

If both conditions are met, \code{drake} will
\enumerate{
\item Assign the most recently-generated admissible data to the target, and
\item skip the target's command.
}
}

\examples{
\dontrun{
isolate_example("quarantine side effects", {
if (requireNamespace("knitr", quietly = TRUE)) {
writeLines(
  c(
    "library(drake)",
    "load_mtcars_example()",
    "drake_config(my_plan, targets = c(\"small\", \"large\"))"
  ),
  "_drake.R" # default value of the `source` argument
)
cat(readLines("_drake.R"), sep = "\n")
r_outdated()
r_make()
r_outdated()
}
})
}
}
\seealso{
\code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{map_plan}
\alias{map_plan}
\title{Create a plan that maps a function to a grid of arguments.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
map_plan(args, fun, id = "id", character_only = FALSE, trace = FALSE)
}
\arguments{
\item{args}{A data frame (or better yet, a \code{tibble})
of function arguments to \code{fun}.
Here, the column names should be the names of the arguments
of \code{fun}, and each row of \code{args} corresponds to a
call to \code{fun}.}

\item{fun}{Name of a function to apply the arguments
row-by-row. Supply a symbol if \code{character_only} is
\code{FALSE} and a character scalar otherwise.}

\item{id}{Name of an optional column in \code{args}
giving the names of the targets. If not supplied,
target names will be generated automatically.
\code{id} should be a symbol if \code{character_only} is \code{FALSE}
and a character scalar otherwise.}

\item{character_only}{Logical, whether to interpret
the \code{fun} and \code{id} arguments as character scalars or symbols.}

\item{trace}{Logical, whether to append the columns of \code{args}
to the output workflow plan data frame. The added columns
help "trace back" the original settings that went into building
each target. Similar to the \code{trace} argument of \code{\link[=drake_plan]{drake_plan()}}.}
}
\value{
A workflow plan data frame.
}
\description{
Deprecated on 2019-05-16. Use \code{\link[=drake_plan]{drake_plan()}}
transformations instead. See
\url{https://books.ropensci.org/drake/plans.html#large-plans}
for the details.
}
\details{
\code{map_plan()} is like \code{base::Map()}:
it takes a function name and a grid of arguments, and
writes out all the commands calls to apply the function to
each row of arguments.
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{make_with_config}
\alias{make_with_config}
\title{Apply make() with a pre-computed config object
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
make_with_config(config)
}
\arguments{
\item{config}{A configuration list returned by \code{\link[=drake_config]{drake_config()}}.}
}
\value{
nothing
}
\description{
Deprecated on 2019-01-04
}
\seealso{
\code{\link[=make]{make()}}, \code{\link[=drake_config]{drake_config()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{cache_path}
\alias{cache_path}
\title{Return the file path where the cache is stored,
if applicable.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
cache_path(cache = NULL)
}
\arguments{
\item{cache}{The cache whose file path you want to know.}
}
\value{
File path where the cache is stored.
}
\description{
Deprecated on 2019-01-12.
}
\details{
Currently only works with
\code{\link[storr:storr_rds]{storr::storr_rds()}} file system caches.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_drake_deps.R
\name{new_drake_deps}
\alias{new_drake_deps}
\title{\code{drake_deps} constructor}
\usage{
new_drake_deps(
  globals = character(0),
  namespaced = character(0),
  strings = character(0),
  loadd = character(0),
  readd = character(0),
  file_in = character(0),
  file_out = character(0),
  knitr_in = character(0)
)
}
\arguments{
\item{globals}{Global symbols found in the expression}

\item{namespaced}{Namespaced objects, e.g. \code{rmarkdown::render}.}

\item{strings}{Miscellaneous strings.}

\item{loadd}{Targets selected with \code{\link[=loadd]{loadd()}}.}

\item{readd}{Targets selected with \code{\link[=readd]{readd()}}.}

\item{file_in}{Literal static file paths enclosed in \code{\link[=file_in]{file_in()}}.}

\item{file_out}{Literal static file paths enclosed in \code{\link[=file_out]{file_out()}}.}

\item{knitr_in}{Literal static file paths enclosed in \code{\link[=knitr_in]{knitr_in()}}.}

\item{restrict}{Optional character vector of allowable names of globals.
If \code{NULL}, all global symbols are detectable. If a character vector,
only the variables in \code{restrict} will count as global variables.}
}
\value{
A \code{drake_deps} object.
}
\description{
List of class \code{drake_deps}.
}
\examples{
if (FALSE) { # stronger than roxygen dontrun
new_drake_deps()
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{file_in}
\alias{file_in}
\title{Declare input files and directories.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
file_in(...)
}
\arguments{
\item{...}{Character vector, paths to files and directories. Use
\code{.id_chr} to refer to the current target by name. \code{.id_chr} is not
limited to use in \code{file_in()} and \code{file_out()}.}
}
\value{
A character vector of declared input file or directory paths.
}
\description{
\code{file_in()} marks individual files
(and whole directories) that your targets depend on.
}
\section{URLs}{

As of \code{drake} 7.4.0, \code{file_in()} and \code{file_out()} have
support for URLs. If the file name begins with
"http://", "https://", or "ftp://", \code{\link[=make]{make()}} attempts
to check the ETag to see if the data changed from last time.
If no ETag can be found, \code{drake} simply uses the ETag
from last \code{\link[=make]{make()}} and registers the file as unchanged
(which prevents your workflow from breaking if you lose
internet access). If your \code{file_in()} URLs require
authentication, see the \code{curl_handles} argument of
\code{make()} and \code{drake_config()} to learn how to supply credentials.
}

\section{Keywords}{

\code{\link[=drake_plan]{drake_plan()}} understands special keyword functions for your commands.
With the exception of \code{\link[=target]{target()}}, each one is a proper function
with its own help file.
\itemize{
\item \code{\link[=target]{target()}}: give the target more than just a command.
Using \code{\link[=target]{target()}}, you can apply a transformation
(examples: \url{https://books.ropensci.org/drake/plans.html#large-plans}), # nolint
supply a trigger (\url{https://books.ropensci.org/drake/triggers.html}), # nolint
or set any number of custom columns.
\item \code{\link[=file_in]{file_in()}}: declare an input file dependency.
\item \code{\link[=file_out]{file_out()}}: declare an output file to be produced
when the target is built.
\item \code{\link[=knitr_in]{knitr_in()}}: declare a \code{knitr} file dependency such as an
R Markdown (\verb{*.Rmd}) or R LaTeX (\verb{*.Rnw}) file.
\item \code{\link[=ignore]{ignore()}}: force \code{drake} to entirely ignore a piece of code:
do not track it for changes and do not analyze it for dependencies.
\item \code{\link[=no_deps]{no_deps()}}: tell \code{drake} to not track the dependencies
of a piece of code. \code{drake} still tracks the code itself for changes.
\item \code{\link[=id_chr]{id_chr()}}: Get the name of the current target.
\item \code{\link[=drake_envir]{drake_envir()}}: get the environment where drake builds targets.
Intended for advanced custom memory management.
}
}

\examples{
\dontrun{
isolate_example("contain side effects", {
# The `file_out()` and `file_in()` functions
# just takes in strings and returns them.
file_out("summaries.txt")
# Their main purpose is to orchestrate your custom files
# in your workflow plan data frame.
plan <- drake_plan(
  out = write.csv(mtcars, file_out("mtcars.csv")),
  contents = read.csv(file_in("mtcars.csv"))
)
plan
# drake knows "\"mtcars.csv\"" is the first target
# and a dependency of `contents`. See for yourself:

make(plan)
file.exists("mtcars.csv")

# You may use `.id_chr` inside `file_out()` and `file_in()`
# to refer  to the current target. This works inside
# static `map()`, `combine()`, `split()`, and `cross()`.

plan <- drake::drake_plan(
  data = target(
    write.csv(data, file_out(paste0(.id_chr, ".csv"))),
    transform = map(data = c(airquality, mtcars))
  )
)
plan

# You can also work with entire directories this way.
# However, in `file_out("your_directory")`, the directory
# becomes an entire unit. Thus, `file_in("your_directory")`
# is more appropriate for subsequent steps than
# `file_in("your_directory/file_inside.txt")`.
plan <- drake_plan(
  out = {
    dir.create(file_out("dir"))
    write.csv(mtcars, "dir/mtcars.csv")
  },
  contents = read.csv(file.path(file_in("dir"), "mtcars.csv"))
)
plan

make(plan)
file.exists("dir/mtcars.csv")

# See the connections that the file relationships create:
if (requireNamespace("visNetwork", quietly = TRUE)) {
  vis_drake_graph(plan)
}
})
}
}
\seealso{
\code{\link[=file_out]{file_out()}}, \code{\link[=knitr_in]{knitr_in()}}, \code{\link[=ignore]{ignore()}}, \code{\link[=no_deps]{no_deps()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{available_hash_algos}
\alias{available_hash_algos}
\title{List the available hash algorithms for drake caches.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
available_hash_algos()
}
\value{
A character vector of names of available hash algorithms.
}
\description{
Deprecated on 2018-12-12.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{render_static_drake_graph}
\alias{render_static_drake_graph}
\title{Deprecated: render a \code{ggraph}/\code{ggplot2} representation
of your drake project.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
render_static_drake_graph(graph_info, main = graph_info$default_title)
}
\arguments{
\item{graph_info}{List of data frames generated by
\code{\link[=drake_graph_info]{drake_graph_info()}}.
There should be 3 data frames: \code{nodes}, \code{edges},
and \code{legend_nodes}.}

\item{main}{Character string, title of the graph.}
}
\value{
A \code{ggplot2} object, which you can modify with more layers,
show with \code{plot()}, or save as a file with \code{ggsave()}.
}
\description{
Use \code{\link[=render_drake_ggraph]{render_drake_ggraph()}} instead.
}
\details{
Deprecated on 2018-07-25.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_plan.R
\name{transformations}
\alias{transformations}
\alias{map}
\alias{split}
\alias{cross}
\alias{combine}
\alias{group}
\title{Transformations in \code{drake_plan()}. \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\arguments{
\item{...}{Grouping variables. New grouping variables must be
supplied with their names and values, existing grouping variables
can be given as symbols without any values assigned.
For dynamic branching, the entries in \code{...} must be unnamed symbols
with no values supplied, and they must be the names of targets.}

\item{.data}{A data frame of new grouping variables with
grouping variable names as column names and values as elements.}

\item{.names}{Literal character vector of names for the targets.
Must be the same length as the targets generated.}

\item{.id}{Symbol or vector of symbols naming grouping variables
to incorporate into target names. Useful for creating short target
names. Set \code{.id = FALSE} to use integer indices as target name suffixes.}

\item{.tag_in}{A symbol or vector of symbols. Tags assign targets
to grouping variables. Use \code{.tag_in} to assign \emph{untransformed}
targets to grouping variables.}

\item{.tag_out}{Just like \code{.tag_in}, except that \code{.tag_out}
assigns \emph{transformed} targets to grouping variables.}

\item{slice}{Number of slices into which \code{split()} partitions the data.}

\item{margin}{Which margin to take the slices in \code{split()}. Same meaning
as the \code{MARGIN} argument of \code{apply()}.}

\item{drop}{Logical, whether to drop a dimension if its length is 1.
Same meaning as \code{mtcars[, 1L, drop = TRUE]} versus
\code{mtcars[, 1L, drop = TRUE]}.}

\item{.by}{Symbol or vector of symbols of grouping variables.
\code{combine()} aggregates/groups targets by the grouping variables in \code{.by}.
For dynamic branching, \code{.by} can only take one variable at a time,
and that variable must be a vector. Ideally, it should take
little space in memory.}

\item{.trace}{Symbol or vector of symbols for the dynamic trace.
The dynamic trace allows you to keep track of the values of
dynamic dependencies are associated with individual sub-targets.
For \code{combine()}, \code{.trace} must either be empty or the same as the
variable given for \code{.by}.
See \code{\link[=get_trace]{get_trace()}} and \code{\link[=read_trace]{read_trace()}} for examples and other details.}
}
\description{
In \code{\link[=drake_plan]{drake_plan()}}, you can define whole batches
of targets with transformations such as
\code{map()}, \code{split()}, \code{cross()}, and \code{combine()}.
}
\details{
For details, see
\url{https://books.ropensci.org/drake/plans.html#large-plans}.
}
\section{Transformations}{

\code{drake} has special syntax for generating large plans.
Your code will look something like
\verb{drake_plan(y = target(f(x), transform = map(x = c(1, 2, 3)))}
You can read about this interface at
\url{https://books.ropensci.org/drake/plans.html#large-plans}. # nolint
}

\section{Static branching}{

In static branching, you define batches of targets
based on information you know in advance.
Overall usage looks like
\verb{drake_plan(<x> = target(<...>, transform = <call>)},
where
\itemize{
\item \verb{<x>} is the name of the target or group of targets.
\item \verb{<...>} is optional arguments to \code{\link[=target]{target()}}.
\item \verb{<call>} is a call to one of the transformation functions.
}

Transformation function usage:
\itemize{
\item \code{map(..., .data, .names, .id, .tag_in, .tag_out)}
\item \code{split(..., slices, margin = 1L, drop = FALSE, .names, .tag_in, .tag_out)} # nolint
\item \code{cross(..., .data, .names, .id, .tag_in, .tag_out)}
\item \code{combine(..., .by, .names, .id, .tag_in, .tag_out)}
}
}

\section{Dynamic branching}{

\itemize{
\item \code{map(..., .trace)}
\item \code{cross(..., .trace)}
\item \code{group(..., .by, .trace)}
}

\code{map()} and \code{cross()} create dynamic sub-targets from the variables
supplied to the dots. As with static branching, the variables
supplied to \code{map()} must all have equal length.
\code{group(f(data), .by = x)} makes new dynamic
sub-targets from \code{data}. Here, \code{data} can be either static or dynamic.
If \code{data} is dynamic, \code{group()} aggregates existing sub-targets.
If \code{data} is static, \code{group()} splits \code{data} into multiple
subsets based on the groupings from \code{.by}.

Differences from static branching:
\itemize{
\item \code{...} must contain \emph{unnamed} symbols with no values supplied,
and they must be the names of targets.
\item Arguments \code{.id}, \code{.tag_in}, and \code{.tag_out} no longer apply.
}
}

\examples{
# Static branching
models <- c("glm", "hierarchical")
plan <- drake_plan(
  data = target(
    get_data(x),
    transform = map(x = c("simulated", "survey"))
  ),
  analysis = target(
    analyze_data(data, model),
    transform = cross(data, model = !!models, .id = c(x, model))
  ),
  summary = target(
    summarize_analysis(analysis),
    transform = map(analysis, .id = c(x, model))
  ),
  results = target(
    bind_rows(summary),
    transform = combine(summary, .by = data)
  )
)
plan
if (requireNamespace("styler")) {
  print(drake_plan_source(plan))
}
# Static splitting
plan <- drake_plan(
  analysis = target(
    analyze(data),
    transform = split(data, slices = 3L, margin = 1L, drop = FALSE)
  )
)
print(plan)
if (requireNamespace("styler", quietly = TRUE)) {
  print(drake_plan_source(plan))
}
# Static tags:
drake_plan(
  x = target(
    command,
    transform = map(y = c(1, 2), .tag_in = from, .tag_out = c(to, out))
  ),
  trace = TRUE
)
plan <- drake_plan(
  survey = target(
    survey_data(x),
    transform = map(x = c(1, 2), .tag_in = source, .tag_out = dataset)
  ),
  download = target(
    download_data(),
    transform = map(y = c(5, 6), .tag_in = source, .tag_out = dataset)
  ),
  analysis = target(
    analyze(dataset),
    transform = map(dataset)
  ),
  results = target(
    bind_rows(analysis),
    transform = combine(analysis, .by = source)
  )
)
plan
if (requireNamespace("styler", quietly = TRUE)) {
  print(drake_plan_source(plan))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_drake_graph.R
\name{render_drake_graph}
\alias{render_drake_graph}
\title{Render a visualization using the data frames
generated by \code{\link[=drake_graph_info]{drake_graph_info()}}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
render_drake_graph(
  graph_info,
  file = character(0),
  layout = NULL,
  direction = NULL,
  hover = TRUE,
  main = graph_info$default_title,
  selfcontained = FALSE,
  navigationButtons = TRUE,
  ncol_legend = 1,
  collapse = TRUE,
  on_select = NULL,
  level_separation = NULL,
  ...
)
}
\arguments{
\item{graph_info}{List of data frames generated by
\code{\link[=drake_graph_info]{drake_graph_info()}}.
There should be 3 data frames: \code{nodes}, \code{edges},
and \code{legend_nodes}.}

\item{file}{Name of a file to save the graph.
If \code{NULL} or \code{character(0)}, no file is saved and
the graph is rendered and displayed within R.
If the file ends in a \code{.png}, \code{.jpg}, \code{.jpeg}, or \code{.pdf} extension,
then a static image will be saved. In this case,
the webshot package and PhantomJS are required:
\verb{install.packages("webshot"); webshot::install_phantomjs()}.
If the file does not end in a \code{.png}, \code{.jpg}, \code{.jpeg}, or \code{.pdf}
extension, an HTML file will be saved, and you can open the
interactive graph using a web browser.}

\item{layout}{Deprecated.}

\item{direction}{Deprecated.}

\item{hover}{Logical, whether to show the command that generated the target
when you hover over a node with the mouse. For imports, the label does not
change with hovering.}

\item{main}{Character string, title of the graph.}

\item{selfcontained}{Logical, whether
to save the \code{file} as a self-contained
HTML file (with external resources base64 encoded) or a file with
external resources placed in an adjacent directory. If \code{TRUE},
pandoc is required. The \code{selfcontained} argument only applies
to HTML files. In other words, if \code{file} is a
PNG, PDF, or JPEG file, for instance,
the point is moot.}

\item{navigationButtons}{Logical, whether to add navigation buttons with
\code{visNetwork::visInteraction(navigationButtons = TRUE)}}

\item{ncol_legend}{Number of columns in the legend nodes.
To remove the legend entirely, set \code{ncol_legend} to \code{NULL} or \code{0}.}

\item{collapse}{Logical, whether to allow nodes to collapse
if you double click on them.
Analogous to \code{visNetwork::visOptions(collapse = TRUE)}.}

\item{on_select}{defines node selection event handling.
Either a string of valid JavaScript that may be passed to
\code{visNetwork::visEvents()}, or one of the following:
\code{TRUE}, \code{NULL}/\code{FALSE}. If \code{TRUE} , enables the default behavior of
opening the link specified by the \code{on_select_col} given to
\code{drake_graph_info()}. \code{NULL}/\code{FALSE} disables the behavior.}

\item{level_separation}{Numeric, \code{levelSeparation} argument to
\code{visNetwork::visHierarchicalLayout()}. Controls the distance
between hierarchical levels. Consider setting if the
aspect ratio of the graph is far from 1.
Defaults to 150 through \code{visNetwork}.}

\item{...}{Arguments passed to \code{visNetwork()}.}
}
\value{
A \code{visNetwork} graph.
}
\description{
This function is called inside
\code{\link[=vis_drake_graph]{vis_drake_graph()}}, which typical users
call more often.
}
\details{
For enhanced interactivity in the graph, see the \code{mandrake}
package: \url{https://github.com/matthewstrasiotto/mandrake}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
if (requireNamespace("visNetwork", quietly = TRUE)) {
# Instead of jumping right to vis_drake_graph(), get the data frames
# of nodes, edges, and legend nodes.
vis_drake_graph(my_plan) # Jump straight to the interactive graph.
# Get the node and edge info that vis_drake_graph() just plotted:
graph <- drake_graph_info(my_plan)
# You can pass the data frames right to render_drake_graph()
# (as in vis_drake_graph()) or you can create
# your own custom visNewtork graph.
render_drake_graph(graph)
}
}
})
}
}
\seealso{
\code{\link[=vis_drake_graph]{vis_drake_graph()}}, \code{\link[=sankey_drake_graph]{sankey_drake_graph()}},
\code{\link[=drake_ggraph]{drake_ggraph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{doc_of_function_call}
\alias{doc_of_function_call}
\title{doc_of_function_call \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
doc_of_function_call(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{short_hash}
\alias{short_hash}
\title{\code{drake} now only uses one hash algorithm per cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
short_hash(cache = drake::get_cache(verbose = verbose), verbose = 1L)
}
\arguments{
\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\value{
A character vector naming a hash algorithm.
}
\description{
Deprecated on 2018-12-12.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{rate_limiting_times}
\alias{rate_limiting_times}
\title{rate_limiting_times \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
rate_limiting_times(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{new_cache}
\alias{new_cache}
\title{Make a new \code{drake} cache.
`\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
new_cache(
  path = NULL,
  verbose = NULL,
  type = NULL,
  hash_algorithm = NULL,
  short_hash_algo = NULL,
  long_hash_algo = NULL,
  ...,
  console_log_file = NULL
)
}
\arguments{
\item{path}{File path to the cache if the cache
is a file system cache.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{type}{Deprecated argument. Once stood for cache type.
Use \code{storr} to customize your caches instead.}

\item{hash_algorithm}{Name of a hash algorithm to use.
See the \code{algo} argument of the \code{digest} package for your options.}

\item{short_hash_algo}{Deprecated on 2018-12-12.
Use \code{hash_algorithm} instead.}

\item{long_hash_algo}{Deprecated on 2018-12-12.
Use \code{hash_algorithm} instead.}

\item{...}{other arguments to the cache constructor.}

\item{console_log_file}{Deprecated on 2019-09-11.}
}
\value{
A newly created drake cache as a storr object.
}
\description{
Uses the \code{\link[=storr_rds]{storr_rds()}} function
from the \code{storr} package.
}
\examples{
\dontrun{
isolate_example("Quarantine new_cache() side effects.", {
clean(destroy = TRUE) # Should not be necessary.
unlink("not_hidden", recursive = TRUE) # Should not be necessary.
cache1 <- new_cache() # Creates a new hidden '.drake' folder.
cache2 <- new_cache(path = "not_hidden", hash_algorithm = "md5")
clean(destroy = TRUE, cache = cache2)
})
}
}
\seealso{
\code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{target_namespaces}
\alias{target_namespaces}
\title{Storr namespaces for targets
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
target_namespaces(default = storr::storr_environment()$default_namespace)
}
\arguments{
\item{default}{Name of the default \code{storr} namespace.}
}
\value{
A character vector of \code{storr} namespaces that store
target-level information.
}
\description{
Deprecated on 2019-01-13.
}
\details{
Ordinary users do not need to worry about this function.
It is just another window into \code{drake}'s internals.
}
\seealso{
\code{\link[=make]{make()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{all_of}
\alias{tar_tidyselect}
\alias{any_of}
\alias{contains}
\alias{ends_with}
\alias{everything}
\alias{last_col}
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
  \item{tidyselect}{\code{\link[tidyselect]{all_of}}, \code{\link[tidyselect:all_of]{any_of}}, \code{\link[tidyselect:starts_with]{contains}}, \code{\link[tidyselect:starts_with]{ends_with}}, \code{\link[tidyselect]{everything}}, \code{\link[tidyselect:everything]{last_col}}, \code{\link[tidyselect:starts_with]{matches}}, \code{\link[tidyselect:starts_with]{num_range}}, \code{\link[tidyselect]{one_of}}, \code{\link[tidyselect]{starts_with}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_drake.R
\name{use_drake}
\alias{use_drake}
\title{Use drake in a project
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#questioning}{\figure{lifecycle-questioning.svg}{options: alt='[Questioning]'}}}{\strong{[Questioning]}}}
\usage{
use_drake(open = interactive())
}
\arguments{
\item{open}{Logical, whether to open \code{make.R} for editing.}
}
\description{
Add top-level R script files to use \code{drake}
in your data analysis project. For details, read
\url{https://books.ropensci.org/drake/projects.html}
}
\details{
Files written:
\enumerate{
\item \code{make.R}: a suggested main R script for batch mode.
\item \verb{_drake.R}: a configuration R script for
the \href{https://docs.ropensci.org/drake/reference/r_make.html}{\verb{r_*()}} functions documented at # nolint
\url{https://books.ropensci.org/drake/projects.html#safer-interactivity}. # nolint
Remarks:
}
\itemize{
\item There is nothing magical about the name, \code{make.R}.
You can call it whatever you want.
\item Other supporting scripts, such as \code{R/packages.R},
\code{R/functions.R}, and \code{R/plan.R}, are not included.
\item You can find examples at
\url{https://github.com/wlandau/drake-examples}
and download examples with \code{\link[=drake_example]{drake_example()}}
(e.g. \code{drake_example("main")}).
}
}
\examples{
\dontrun{
# use_drake(open = FALSE) # nolint
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_drake_deps_ht.R
\name{drake_deps_ht}
\alias{drake_deps_ht}
\title{\code{drake_deps_ht} helper}
\usage{
drake_deps_ht(expr, exclude = character(0), restrict = NULL)
}
\arguments{
\item{expr}{An R expression}

\item{exclude}{Character vector of the names of symbols to exclude
from the code analysis.}

\item{restrict}{Optional character vector of allowable names of globals.
If \code{NULL}, all global symbols are detectable. If a character vector,
only the variables in \code{restrict} will count as global variables.}
}
\value{
A \code{drake_deps_ht} object.
}
\description{
Static code analysis.
}
\examples{
if (FALSE) { # stronger than roxygen dontrun
expr <- quote({
  a <- base::list(1)
  b <- seq_len(10)
  file_out("abc")
  file_in("xyz")
  x <- "123"
  loadd(abc)
  readd(xyz)
})
drake_deps_ht(expr)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{load_main_example}
\alias{load_main_example}
\title{Load the main example.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
load_main_example(
  envir = parent.frame(),
  report_file = "report.Rmd",
  overwrite = FALSE,
  force = FALSE
)
}
\arguments{
\item{envir}{The environment to load the example into.
Defaults to your workspace.
For an insulated workspace,
set \code{envir = new.env(parent = globalenv())}.}

\item{report_file}{Where to write the report file \code{report.Rmd}.}

\item{overwrite}{Logical, whether to overwrite an
existing file \code{report.Rmd}}

\item{force}{Deprecated.}
}
\value{
A \code{\link[=drake_config]{drake_config()}} configuration list.
}
\description{
The main example lives at
\url{https://github.com/wlandau/drake-examples/tree/main/main}.
Use \code{drake_example("main")} to download its code.
This function also writes/overwrites
the files \code{report.Rmd} and \code{raw_data.xlsx}.
}
\details{
Deprecated 2018-12-31.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{example_drake}
\alias{example_drake}
\title{example_drake \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
example_drake(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{knitr_deps}
\alias{knitr_deps}
\title{Dependencies of a knitr report
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
knitr_deps(target)
}
\arguments{
\item{target}{Encoded file path}
}
\value{
Data frame of dependencies
}
\description{
Deprecated on 2019-02-14
\code{knit("your_report.Rmd")} or
\code{knit("your_report.Rmd", quiet = TRUE)}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_drake_settings.R
\name{new_drake_settings}
\alias{new_drake_settings}
\title{\code{drake_settings} constructor}
\usage{
new_drake_settings(
  cache_log_file = NULL,
  curl_handles = NULL,
  garbage_collection = NULL,
  jobs = NULL,
  jobs_preprocess = NULL,
  keep_going = NULL,
  lazy_load = NULL,
  lib_loc = NULL,
  lock_envir = NULL,
  lock_cache = NULL,
  log_build_times = NULL,
  log_progress = NULL,
  memory_strategy = NULL,
  parallelism = NULL,
  recover = NULL,
  recoverable = NULL,
  seed = NULL,
  session_info = NULL,
  skip_imports = NULL,
  skip_safety_checks = NULL,
  skip_targets = NULL,
  sleep = NULL,
  template = NULL,
  log_worker = NULL
)
}
\arguments{
\item{cache_log_file}{Name of the CSV cache log file to write.
If \code{TRUE}, the default file name is used (\code{drake_cache.CSV}).
If \code{NULL}, no file is written.
If activated, this option writes a flat text file
to represent the state of the cache
(fingerprints of all the targets and imports).
If you put the log file under version control, your commit history
will give you an easy representation of how your results change
over time as the rest of your project changes. Hopefully,
this is a step in the right direction for data reproducibility.}

\item{curl_handles}{A named list of curl handles. Each value is an
object from \code{curl::new_handle()}, and each name is a URL
(and should start with "http", "https", or "ftp").
Example:
list(
\verb{http://httpbin.org/basic-auth} = curl::new_handle(
username = "user", password = "passwd"
)
)
Then, if your plan has
\code{file_in("http://httpbin.org/basic-auth/user/passwd")}
\code{drake} will authenticate using the username and password of the handle
for \verb{http://httpbin.org/basic-auth/}.

\code{drake} uses partial matching on text to
find the right handle of the \code{file_in()} URL, so the name of the handle
could be the complete URL (\code{"http://httpbin.org/basic-auth/user/passwd"})
or a part of the URL (e.g. \code{"http://httpbin.org/"} or
\code{"http://httpbin.org/basic-auth/"}). If you have multiple handles
whose names match your URL, \code{drake} will choose the closest match.}

\item{garbage_collection}{Logical, whether to call \code{gc()} each time
a target is built during \code{\link[=make]{make()}}.}

\item{jobs}{Maximum number of parallel workers for processing the targets.
You can experiment with \code{\link[=predict_runtime]{predict_runtime()}}
to help decide on an appropriate number of jobs.
For details, visit
\url{https://books.ropensci.org/drake/time.html}.}

\item{jobs_preprocess}{Number of parallel jobs for processing the imports
and doing other preprocessing tasks.}

\item{keep_going}{Logical, whether to still keep running \code{\link[=make]{make()}}
if targets fail.}

\item{lazy_load}{An old feature, currently being questioned.
For the current recommendations on memory management, see
\url{https://books.ropensci.org/drake/memory.html#memory-strategies}.
The \code{lazy_load} argument is either a character vector or a logical.
For dynamic targets, the behavior is always \code{"eager"} (see below).
So the \code{lazy_load} argument is for static targets only.
Choices for \code{lazy_load}:
\itemize{
\item \code{"eager"}: no lazy loading. The target is loaded right away
with \code{\link[=assign]{assign()}}.
\item \code{"promise"}: lazy loading with \code{\link[=delayedAssign]{delayedAssign()}}
\item \code{"bind"}: lazy loading with active bindings:
\code{bindr::populate_env()}.
\item \code{TRUE}: same as \code{"promise"}.
\item \code{FALSE}: same as \code{"eager"}.
}

If \code{lazy_load} is \code{"eager"},
drake prunes the execution environment before each target/stage,
removing all superfluous targets
and then loading any dependencies it will need for building.
In other words, drake prepares the environment in advance
and tries to be memory efficient.
If \code{lazy_load} is \code{"bind"} or \code{"promise"}, drake assigns
promises to load any dependencies at the last minute.
Lazy loading may be more memory efficient in some use cases, but
it may duplicate the loading of dependencies, costing time.}

\item{lib_loc}{Character vector, optional.
Same as in \code{library()} or \code{require()}.
Applies to the \code{packages} argument (see above).}

\item{lock_envir}{Logical, whether to lock \code{config$envir} during \code{make()}.
If \code{TRUE}, \code{make()} quits in error whenever a command in your
\code{drake} plan (or \code{prework}) tries to add, remove, or modify
non-hidden variables in your environment/workspace/R session.
This is extremely important for ensuring the purity of your functions
and the reproducibility/credibility/trust you can place in your project.
\code{lock_envir} will be set to a default of \code{TRUE} in \code{drake} version
7.0.0 and higher. Namespaces are never locked, e.g.
if \code{envir} is \code{getNamespace("packagename")}.}

\item{lock_cache}{Logical, whether to lock the cache before running \code{make()}
etc. It is usually recommended to keep cache locking on.
However, if you interrupt \code{make()} before it can clean itself up,
then the cache will stay locked,
and you will need to manually unlock it with
\code{drake::drake_cache("xyz")$unlock()}. Repeatedly unlocking the cache
by hand is annoying, and \code{lock_cache = FALSE} prevents the cache
from locking in the first place.}

\item{log_build_times}{Logical, whether to record build_times for targets.
Mac users may notice a 20\% speedup in \code{make()}
with \code{build_times = FALSE}.}

\item{log_progress}{Logical, whether to log the progress
of individual targets as they are being built. Progress logging
creates extra files in the cache (usually the \verb{.drake/} folder)
and slows down \code{make()} a little.
If you need to reduce or limit the number of files in the cache,
call \code{make(log_progress = FALSE, recover = FALSE)}.}

\item{memory_strategy}{Character scalar, name of the
strategy \code{drake} uses to load/unload a target's dependencies in memory.
You can give each target its own memory strategy,
(e.g. \code{drake_plan(x = 1, y = target(f(x), memory_strategy = "lookahead"))})
to override the global memory strategy. Choices:
\itemize{
\item \code{"speed"}: Once a target is newly built or loaded in memory,
just keep it there.
This choice maximizes speed and hogs memory.
\item \code{"autoclean"}: Just before building each new target,
unload everything from memory except the target's direct dependencies.
After a target is built, discard it from memory.
(Set \code{garbage_collection = TRUE} to make sure it is really gone.)
This option conserves memory, but it sacrifices speed because
each new target needs to reload
any previously unloaded targets from storage.
\item \code{"preclean"}: Just before building each new target,
unload everything from memory except the target's direct dependencies.
After a target is built, keep it in memory until \code{drake} determines
they can be unloaded.
This option conserves memory, but it sacrifices speed because
each new target needs to reload
any previously unloaded targets from storage.
\item \code{"lookahead"}: Just before building each new target,
search the dependency graph to find targets that will not be
needed for the rest of the current \code{make()} session.
After a target is built, keep it in memory until the next
memory management stage.
In this mode, targets are only in memory if they need to be loaded,
and we avoid superfluous reads from the cache.
However, searching the graph takes time,
and it could even double the computational overhead for large projects.
\item \code{"unload"}: Just before building each new target,
unload all targets from memory.
After a target is built, \strong{do not} keep it in memory.
This mode aggressively optimizes for both memory and speed,
but in commands and triggers,
you have to manually load any dependencies you need using \code{readd()}.
\item \code{"none"}: Do not manage memory at all.
Do not load or unload anything before building targets.
After a target is built, \strong{do not} keep it in memory.
This mode aggressively optimizes for both memory and speed,
but in commands and triggers,
you have to manually load any dependencies you need using \code{readd()}.
}

For even more direct
control over which targets \code{drake} keeps in memory, see the
help file examples of \code{\link[=drake_envir]{drake_envir()}}.
Also see the \code{garbage_collection} argument of \code{make()} and
\code{drake_config()}.}

\item{parallelism}{Character scalar, type of parallelism to use.
For detailed explanations, see the
\href{https://books.ropensci.org/drake/hpc.html}{high-performance computing chapter} # nolint
of the user manual.

You could also supply your own scheduler function
if you want to experiment or aggressively optimize.
The function should take a single \code{config} argument
(produced by \code{\link[=drake_config]{drake_config()}}). Existing examples
from \code{drake}'s internals are the \verb{backend_*()} functions:
\itemize{
\item \code{backend_loop()}
\item \code{backend_clustermq()}
\item \code{backend_future()}
However, this functionality is really a back door
and should not be used for production purposes unless you really
know what you are doing and you are willing to suffer setbacks
whenever \code{drake}'s unexported core functions are updated.
}}

\item{recover}{Logical, whether to activate automated data recovery.
The default is \code{FALSE} because
\enumerate{
\item Automated data recovery is still stable.
\item It has reproducibility issues.
Targets recovered from the distant past may have been generated
with earlier versions of R and earlier package environments
that no longer exist.
\item It is not always possible, especially when dynamic files
are combined with dynamic branching
(e.g. \code{dynamic = map(stuff)} and \code{format = "file"} etc.)
since behavior is harder to predict in advance.
}

How it works: if \code{recover} is \code{TRUE},
\code{drake} tries to salvage old target values from the cache
instead of running commands from the plan.
A target is recoverable if
\enumerate{
\item There is an old value somewhere in the cache that
shares the command, dependencies, etc.
of the target about to be built.
\item The old value was generated with \code{make(recoverable = TRUE)}.
}

If both conditions are met, \code{drake} will
\enumerate{
\item Assign the most recently-generated admissible data to the target, and
\item skip the target's command.
}

Functions \code{\link[=recoverable]{recoverable()}} and \code{\link[=r_recoverable]{r_recoverable()}} show the most upstream
outdated targets that will be recovered in this way in the next
\code{\link[=make]{make()}} or \code{\link[=r_make]{r_make()}}.}

\item{recoverable}{Logical, whether to make target values recoverable
with \code{make(recover = TRUE)}.
This requires writing extra files to the cache,
and it prevents old metadata from being removed with garbage collection
(\code{clean(garbage_collection = TRUE)}, \code{gc()} in \code{storr}s).
If you need to limit the cache size or the number of files in the cache,
consider \code{make(recoverable = FALSE, progress = FALSE)}.
Recovery is not always possible, especially when dynamic files
are combined with dynamic branching
(e.g. \code{dynamic = map(stuff)} and \code{format = "file"} etc.)
since behavior is harder to predict in advance.}

\item{seed}{Integer, the root pseudo-random number generator
seed to use for your project.
In \code{\link[=make]{make()}}, \code{drake} generates a unique
local seed for each target using the global seed
and the target name. That way, different pseudo-random numbers
are generated for different targets, and this pseudo-randomness
is reproducible.

To ensure reproducibility across different R sessions,
\code{set.seed()} and \code{.Random.seed} are ignored and have no affect on
\code{drake} workflows. Conversely, \code{make()} does not usually
change \code{.Random.seed},
even when pseudo-random numbers are generated.
The exception to this last point is
\code{make(parallelism = "clustermq")}
because the \code{clustermq} package needs to generate random numbers
to set up ports and sockets for ZeroMQ.

On the first call to \code{make()} or \code{drake_config()}, \code{drake}
uses the random number generator seed from the \code{seed} argument.
Here, if the \code{seed} is \code{NULL} (default), \code{drake} uses a \code{seed} of \code{0}.
On subsequent \code{make()}s for existing projects, the project's
cached seed will be used in order to ensure reproducibility.
Thus, the \code{seed} argument must either be \code{NULL} or the same
seed from the project's cache (usually the \verb{.drake/} folder).
To reset the random number generator seed for a project,
use \code{clean(destroy = TRUE)}.}

\item{session_info}{Logical, whether to save the \code{sessionInfo()}
to the cache. Defaults to \code{TRUE}.
This behavior is recommended for serious \code{\link[=make]{make()}}s
for the sake of reproducibility. This argument only exists to
speed up tests. Apparently, \code{sessionInfo()} is a bottleneck
for small \code{\link[=make]{make()}}s.}

\item{skip_imports}{Logical, whether to totally neglect to
process the imports and jump straight to the targets. This can be useful
if your imports are massive and you just want to test your project,
but it is bad practice for reproducible data analysis.
This argument is overridden if you supply your own \code{graph} argument.}

\item{skip_safety_checks}{Logical, whether to skip the safety checks
on your workflow. Use at your own peril.}

\item{skip_targets}{Logical, whether to skip building the targets
in \code{plan} and just import objects and files.}

\item{sleep}{Optional function on a single numeric argument \code{i}.
Default: \code{function(i) 0.01}.

To conserve memory, \code{drake} assigns a brand new closure to
\code{sleep}, so your custom function should not depend on in-memory data
except from loaded packages.

For parallel processing, \code{drake} uses
a central main process to check what the parallel
workers are doing, and for the affected high-performance
computing workflows, wait for data to arrive over a network.
In between loop iterations, the main process sleeps to avoid throttling.
The \code{sleep} argument to \code{make()} and \code{drake_config()}
allows you to customize how much time the main process spends
sleeping.

The \code{sleep} argument is a function that takes an argument
\code{i} and returns a numeric scalar, the number of seconds to
supply to \code{Sys.sleep()} after iteration \code{i} of checking.
(Here, \code{i} starts at 1.)
If the checking loop does something other than sleeping
on iteration \code{i}, then \code{i} is reset back to 1.

To sleep for the same amount of time between checks,
you might supply something like \code{function(i) 0.01}.
But to avoid consuming too many resources during heavier
and longer workflows, you might use an exponential
back-off: say,
\code{function(i) { 0.1 + 120 * pexp(i - 1, rate = 0.01) }}.}

\item{template}{A named list of values to fill in the \code{{{ ... }}}
placeholders in template files (e.g. from \code{\link[=drake_hpc_template_file]{drake_hpc_template_file()}}).
Same as the \code{template} argument of \code{clustermq::Q()} and
\code{clustermq::workers}.
Enabled for \code{clustermq} only (\code{make(parallelism = "clustermq")}),
not \code{future} or \code{batchtools} so far.
For more information, see the \code{clustermq} package:
\url{https://github.com/mschubert/clustermq}.
Some template placeholders such as \code{{{ job_name }}} and \code{{{ n_jobs }}}
cannot be set this way.}

\item{log_worker}{Logical, same as the \code{log_worker} argument of
\code{clustermq::workers()} and \code{clustermq::Q()}. Only relevant
if \code{parallelism} is \code{"clustermq"}.}
}
\value{
A \code{drake_settings} object.
}
\description{
List of class \code{drake_settings}.
}
\examples{
if (FALSE) { # stronger than roxygen dontrun
new_drake_settings()
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sankey_drake_graph.R
\name{sankey_drake_graph_impl}
\alias{sankey_drake_graph_impl}
\title{Internal function with a drake_config() argument}
\usage{
sankey_drake_graph_impl(
  config,
  file = character(0),
  selfcontained = FALSE,
  build_times = "build",
  digits = 3,
  targets_only = FALSE,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  make_imports = TRUE,
  from_scratch = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE
)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_config.R
\name{drake_config}
\alias{drake_config}
\title{Ending of _drake.R for r_make() and friends
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_config(
  plan,
  targets = NULL,
  envir = parent.frame(),
  verbose = 1L,
  hook = NULL,
  cache = drake::drake_cache(),
  fetch_cache = NULL,
  parallelism = "loop",
  jobs = 1L,
  jobs_preprocess = 1L,
  packages = rev(.packages()),
  lib_loc = NULL,
  prework = character(0),
  prepend = NULL,
  command = NULL,
  args = NULL,
  recipe_command = NULL,
  timeout = NULL,
  cpu = Inf,
  elapsed = Inf,
  retries = 0,
  force = FALSE,
  log_progress = TRUE,
  graph = NULL,
  trigger = drake::trigger(),
  skip_targets = FALSE,
  skip_imports = FALSE,
  skip_safety_checks = FALSE,
  lazy_load = "eager",
  session_info = NULL,
  cache_log_file = NULL,
  seed = NULL,
  caching = c("main", "master", "worker"),
  keep_going = FALSE,
  session = NULL,
  pruning_strategy = NULL,
  makefile_path = NULL,
  console_log_file = NULL,
  ensure_workers = NULL,
  garbage_collection = FALSE,
  template = list(),
  sleep = function(i) 0.01,
  hasty_build = NULL,
  memory_strategy = "speed",
  spec = NULL,
  layout = NULL,
  lock_envir = TRUE,
  history = TRUE,
  recover = FALSE,
  recoverable = TRUE,
  curl_handles = list(),
  max_expand = NULL,
  log_build_times = TRUE,
  format = NULL,
  lock_cache = TRUE,
  log_make = NULL,
  log_worker = FALSE
)
}
\arguments{
\item{plan}{Workflow plan data frame.
A workflow plan data frame is a data frame
with a \code{target} column and a \code{command} column.
(See the details in the \code{\link[=drake_plan]{drake_plan()}} help file
for descriptions of the optional columns.)
Targets are the objects that drake generates,
and commands are the pieces of R code that produce them.
You can create and track custom files along the way
(see \code{\link[=file_in]{file_in()}}, \code{\link[=file_out]{file_out()}}, and \code{\link[=knitr_in]{knitr_in()}}).
Use the function \code{\link[=drake_plan]{drake_plan()}} to generate workflow plan
data frames.}

\item{targets}{Character vector, names of targets to build.
Dependencies are built too. You may supply static and/or whole
dynamic targets, but no sub-targets.}

\item{envir}{Environment to use. Defaults to the current
workspace, so you should not need to worry about this
most of the time. A deep copy of \code{envir} is made,
so you don't need to worry about your workspace being modified
by \code{make}. The deep copy inherits from the global environment.
Wherever necessary, objects and functions are imported
from \code{envir} and the global environment and
then reproducibly tracked as dependencies.}

\item{verbose}{Integer, control printing to the console/terminal.
\itemize{
\item \code{0}: print nothing.
\item \code{1}: print target-by-target messages as \code{\link[=make]{make()}} progresses.
\item \code{2}: show a progress bar to track how many targets are
done so far.
}}

\item{hook}{Deprecated.}

\item{cache}{drake cache as created by \code{\link[=new_cache]{new_cache()}}.
See also \code{\link[=drake_cache]{drake_cache()}}.}

\item{fetch_cache}{Deprecated.}

\item{parallelism}{Character scalar, type of parallelism to use.
For detailed explanations, see the
\href{https://books.ropensci.org/drake/hpc.html}{high-performance computing chapter} # nolint
of the user manual.

You could also supply your own scheduler function
if you want to experiment or aggressively optimize.
The function should take a single \code{config} argument
(produced by \code{\link[=drake_config]{drake_config()}}). Existing examples
from \code{drake}'s internals are the \verb{backend_*()} functions:
\itemize{
\item \code{backend_loop()}
\item \code{backend_clustermq()}
\item \code{backend_future()}
However, this functionality is really a back door
and should not be used for production purposes unless you really
know what you are doing and you are willing to suffer setbacks
whenever \code{drake}'s unexported core functions are updated.
}}

\item{jobs}{Maximum number of parallel workers for processing the targets.
You can experiment with \code{\link[=predict_runtime]{predict_runtime()}}
to help decide on an appropriate number of jobs.
For details, visit
\url{https://books.ropensci.org/drake/time.html}.}

\item{jobs_preprocess}{Number of parallel jobs for processing the imports
and doing other preprocessing tasks.}

\item{packages}{Character vector packages to load, in the order
they should be loaded. Defaults to \code{rev(.packages())}, so you
should not usually need to set this manually. Just call
\code{\link[=library]{library()}} to load your packages before \code{make()}.
However, sometimes packages need to be strictly forced to load
in a certain order, especially if \code{parallelism} is
\code{"Makefile"}. To do this, do not use \code{\link[=library]{library()}}
or \code{\link[=require]{require()}} or \code{\link[=loadNamespace]{loadNamespace()}} or
\code{\link[=attachNamespace]{attachNamespace()}} to load any libraries beforehand.
Just list your packages in the \code{packages} argument in the order
you want them to be loaded.}

\item{lib_loc}{Character vector, optional.
Same as in \code{library()} or \code{require()}.
Applies to the \code{packages} argument (see above).}

\item{prework}{Expression (language object), list of expressions,
or character vector.
Code to run right before targets build.
Called only once if \code{parallelism} is \code{"loop"}
and once per target otherwise.
This code can be used to set global options, etc.}

\item{prepend}{Deprecated.}

\item{command}{Deprecated.}

\item{args}{Deprecated.}

\item{recipe_command}{Deprecated.}

\item{timeout}{\code{deprecated}. Use \code{elapsed} and \code{cpu} instead.}

\item{cpu}{Same as the \code{cpu} argument of \code{setTimeLimit()}.
Seconds of cpu time before a target times out.
Assign target-level cpu timeout times with an optional \code{cpu}
column in \code{plan}.}

\item{elapsed}{Same as the \code{elapsed} argument of \code{setTimeLimit()}.
Seconds of elapsed time before a target times out.
Assign target-level elapsed timeout times with an optional \code{elapsed}
column in \code{plan}.}

\item{retries}{Number of retries to execute if the target fails.
Assign target-level retries with an optional \code{retries}
column in \code{plan}.}

\item{force}{Logical. If \code{FALSE} (default) then \code{drake}
imposes checks if the cache was created with an old
and incompatible version of drake.
If there is an incompatibility, \code{make()} stops to
give you an opportunity to
downgrade \code{drake} to a compatible version
rather than rerun all your targets from scratch.}

\item{log_progress}{Logical, whether to log the progress
of individual targets as they are being built. Progress logging
creates extra files in the cache (usually the \verb{.drake/} folder)
and slows down \code{make()} a little.
If you need to reduce or limit the number of files in the cache,
call \code{make(log_progress = FALSE, recover = FALSE)}.}

\item{graph}{Deprecated.}

\item{trigger}{Name of the trigger to apply to all targets.
Ignored if \code{plan} has a \code{trigger} column.
See \code{\link[=trigger]{trigger()}} for details.}

\item{skip_targets}{Logical, whether to skip building the targets
in \code{plan} and just import objects and files.}

\item{skip_imports}{Logical, whether to totally neglect to
process the imports and jump straight to the targets. This can be useful
if your imports are massive and you just want to test your project,
but it is bad practice for reproducible data analysis.
This argument is overridden if you supply your own \code{graph} argument.}

\item{skip_safety_checks}{Logical, whether to skip the safety checks
on your workflow. Use at your own peril.}

\item{lazy_load}{An old feature, currently being questioned.
For the current recommendations on memory management, see
\url{https://books.ropensci.org/drake/memory.html#memory-strategies}.
The \code{lazy_load} argument is either a character vector or a logical.
For dynamic targets, the behavior is always \code{"eager"} (see below).
So the \code{lazy_load} argument is for static targets only.
Choices for \code{lazy_load}:
\itemize{
\item \code{"eager"}: no lazy loading. The target is loaded right away
with \code{\link[=assign]{assign()}}.
\item \code{"promise"}: lazy loading with \code{\link[=delayedAssign]{delayedAssign()}}
\item \code{"bind"}: lazy loading with active bindings:
\code{bindr::populate_env()}.
\item \code{TRUE}: same as \code{"promise"}.
\item \code{FALSE}: same as \code{"eager"}.
}

If \code{lazy_load} is \code{"eager"},
drake prunes the execution environment before each target/stage,
removing all superfluous targets
and then loading any dependencies it will need for building.
In other words, drake prepares the environment in advance
and tries to be memory efficient.
If \code{lazy_load} is \code{"bind"} or \code{"promise"}, drake assigns
promises to load any dependencies at the last minute.
Lazy loading may be more memory efficient in some use cases, but
it may duplicate the loading of dependencies, costing time.}

\item{session_info}{Logical, whether to save the \code{sessionInfo()}
to the cache. Defaults to \code{TRUE}.
This behavior is recommended for serious \code{\link[=make]{make()}}s
for the sake of reproducibility. This argument only exists to
speed up tests. Apparently, \code{sessionInfo()} is a bottleneck
for small \code{\link[=make]{make()}}s.}

\item{cache_log_file}{Name of the CSV cache log file to write.
If \code{TRUE}, the default file name is used (\code{drake_cache.CSV}).
If \code{NULL}, no file is written.
If activated, this option writes a flat text file
to represent the state of the cache
(fingerprints of all the targets and imports).
If you put the log file under version control, your commit history
will give you an easy representation of how your results change
over time as the rest of your project changes. Hopefully,
this is a step in the right direction for data reproducibility.}

\item{seed}{Integer, the root pseudo-random number generator
seed to use for your project.
In \code{\link[=make]{make()}}, \code{drake} generates a unique
local seed for each target using the global seed
and the target name. That way, different pseudo-random numbers
are generated for different targets, and this pseudo-randomness
is reproducible.

To ensure reproducibility across different R sessions,
\code{set.seed()} and \code{.Random.seed} are ignored and have no affect on
\code{drake} workflows. Conversely, \code{make()} does not usually
change \code{.Random.seed},
even when pseudo-random numbers are generated.
The exception to this last point is
\code{make(parallelism = "clustermq")}
because the \code{clustermq} package needs to generate random numbers
to set up ports and sockets for ZeroMQ.

On the first call to \code{make()} or \code{drake_config()}, \code{drake}
uses the random number generator seed from the \code{seed} argument.
Here, if the \code{seed} is \code{NULL} (default), \code{drake} uses a \code{seed} of \code{0}.
On subsequent \code{make()}s for existing projects, the project's
cached seed will be used in order to ensure reproducibility.
Thus, the \code{seed} argument must either be \code{NULL} or the same
seed from the project's cache (usually the \verb{.drake/} folder).
To reset the random number generator seed for a project,
use \code{clean(destroy = TRUE)}.}

\item{caching}{Character string, either \code{"main"} or \code{"worker"}.
\itemize{
\item \code{"main"}: Targets are built by remote workers and sent back to
the main process. Then, the main process saves them to the
cache (\code{config$cache}, usually a file system \code{storr}).
Appropriate if remote workers do not have access to the file system
of the calling R session. Targets are cached one at a time,
which may be slow in some situations.
\item \code{"worker"}: Remote workers not only build the targets, but also
save them to the cache. Here, caching happens in parallel.
However, remote workers need to have access to the file system
of the calling R session. Transferring target data across
a network can be slow.
}}

\item{keep_going}{Logical, whether to still keep running \code{\link[=make]{make()}}
if targets fail.}

\item{session}{Deprecated. Has no effect now.}

\item{pruning_strategy}{Deprecated. See \code{memory_strategy}.}

\item{makefile_path}{Deprecated.}

\item{console_log_file}{Deprecated in favor of \code{log_make}.}

\item{ensure_workers}{Deprecated.}

\item{garbage_collection}{Logical, whether to call \code{gc()} each time
a target is built during \code{\link[=make]{make()}}.}

\item{template}{A named list of values to fill in the \code{{{ ... }}}
placeholders in template files (e.g. from \code{\link[=drake_hpc_template_file]{drake_hpc_template_file()}}).
Same as the \code{template} argument of \code{clustermq::Q()} and
\code{clustermq::workers}.
Enabled for \code{clustermq} only (\code{make(parallelism = "clustermq")}),
not \code{future} or \code{batchtools} so far.
For more information, see the \code{clustermq} package:
\url{https://github.com/mschubert/clustermq}.
Some template placeholders such as \code{{{ job_name }}} and \code{{{ n_jobs }}}
cannot be set this way.}

\item{sleep}{Optional function on a single numeric argument \code{i}.
Default: \code{function(i) 0.01}.

To conserve memory, \code{drake} assigns a brand new closure to
\code{sleep}, so your custom function should not depend on in-memory data
except from loaded packages.

For parallel processing, \code{drake} uses
a central main process to check what the parallel
workers are doing, and for the affected high-performance
computing workflows, wait for data to arrive over a network.
In between loop iterations, the main process sleeps to avoid throttling.
The \code{sleep} argument to \code{make()} and \code{drake_config()}
allows you to customize how much time the main process spends
sleeping.

The \code{sleep} argument is a function that takes an argument
\code{i} and returns a numeric scalar, the number of seconds to
supply to \code{Sys.sleep()} after iteration \code{i} of checking.
(Here, \code{i} starts at 1.)
If the checking loop does something other than sleeping
on iteration \code{i}, then \code{i} is reset back to 1.

To sleep for the same amount of time between checks,
you might supply something like \code{function(i) 0.01}.
But to avoid consuming too many resources during heavier
and longer workflows, you might use an exponential
back-off: say,
\code{function(i) { 0.1 + 120 * pexp(i - 1, rate = 0.01) }}.}

\item{hasty_build}{Deprecated}

\item{memory_strategy}{Character scalar, name of the
strategy \code{drake} uses to load/unload a target's dependencies in memory.
You can give each target its own memory strategy,
(e.g. \code{drake_plan(x = 1, y = target(f(x), memory_strategy = "lookahead"))})
to override the global memory strategy. Choices:
\itemize{
\item \code{"speed"}: Once a target is newly built or loaded in memory,
just keep it there.
This choice maximizes speed and hogs memory.
\item \code{"autoclean"}: Just before building each new target,
unload everything from memory except the target's direct dependencies.
After a target is built, discard it from memory.
(Set \code{garbage_collection = TRUE} to make sure it is really gone.)
This option conserves memory, but it sacrifices speed because
each new target needs to reload
any previously unloaded targets from storage.
\item \code{"preclean"}: Just before building each new target,
unload everything from memory except the target's direct dependencies.
After a target is built, keep it in memory until \code{drake} determines
they can be unloaded.
This option conserves memory, but it sacrifices speed because
each new target needs to reload
any previously unloaded targets from storage.
\item \code{"lookahead"}: Just before building each new target,
search the dependency graph to find targets that will not be
needed for the rest of the current \code{make()} session.
After a target is built, keep it in memory until the next
memory management stage.
In this mode, targets are only in memory if they need to be loaded,
and we avoid superfluous reads from the cache.
However, searching the graph takes time,
and it could even double the computational overhead for large projects.
\item \code{"unload"}: Just before building each new target,
unload all targets from memory.
After a target is built, \strong{do not} keep it in memory.
This mode aggressively optimizes for both memory and speed,
but in commands and triggers,
you have to manually load any dependencies you need using \code{readd()}.
\item \code{"none"}: Do not manage memory at all.
Do not load or unload anything before building targets.
After a target is built, \strong{do not} keep it in memory.
This mode aggressively optimizes for both memory and speed,
but in commands and triggers,
you have to manually load any dependencies you need using \code{readd()}.
}

For even more direct
control over which targets \code{drake} keeps in memory, see the
help file examples of \code{\link[=drake_envir]{drake_envir()}}.
Also see the \code{garbage_collection} argument of \code{make()} and
\code{drake_config()}.}

\item{spec}{Deprecated.}

\item{layout}{Deprecated.}

\item{lock_envir}{Logical, whether to lock \code{config$envir} during \code{make()}.
If \code{TRUE}, \code{make()} quits in error whenever a command in your
\code{drake} plan (or \code{prework}) tries to add, remove, or modify
non-hidden variables in your environment/workspace/R session.
This is extremely important for ensuring the purity of your functions
and the reproducibility/credibility/trust you can place in your project.
\code{lock_envir} will be set to a default of \code{TRUE} in \code{drake} version
7.0.0 and higher. Namespaces are never locked, e.g.
if \code{envir} is \code{getNamespace("packagename")}.}

\item{history}{Logical, whether to record the build history
of your targets. You can also supply a
\href{https://github.com/wlandau/txtq}{\code{txtq}}, which is
how \code{drake} records history.
Must be \code{TRUE} for \code{\link[=drake_history]{drake_history()}} to work later.}

\item{recover}{Logical, whether to activate automated data recovery.
The default is \code{FALSE} because
\enumerate{
\item Automated data recovery is still stable.
\item It has reproducibility issues.
Targets recovered from the distant past may have been generated
with earlier versions of R and earlier package environments
that no longer exist.
\item It is not always possible, especially when dynamic files
are combined with dynamic branching
(e.g. \code{dynamic = map(stuff)} and \code{format = "file"} etc.)
since behavior is harder to predict in advance.
}

How it works: if \code{recover} is \code{TRUE},
\code{drake} tries to salvage old target values from the cache
instead of running commands from the plan.
A target is recoverable if
\enumerate{
\item There is an old value somewhere in the cache that
shares the command, dependencies, etc.
of the target about to be built.
\item The old value was generated with \code{make(recoverable = TRUE)}.
}

If both conditions are met, \code{drake} will
\enumerate{
\item Assign the most recently-generated admissible data to the target, and
\item skip the target's command.
}

Functions \code{\link[=recoverable]{recoverable()}} and \code{\link[=r_recoverable]{r_recoverable()}} show the most upstream
outdated targets that will be recovered in this way in the next
\code{\link[=make]{make()}} or \code{\link[=r_make]{r_make()}}.}

\item{recoverable}{Logical, whether to make target values recoverable
with \code{make(recover = TRUE)}.
This requires writing extra files to the cache,
and it prevents old metadata from being removed with garbage collection
(\code{clean(garbage_collection = TRUE)}, \code{gc()} in \code{storr}s).
If you need to limit the cache size or the number of files in the cache,
consider \code{make(recoverable = FALSE, progress = FALSE)}.
Recovery is not always possible, especially when dynamic files
are combined with dynamic branching
(e.g. \code{dynamic = map(stuff)} and \code{format = "file"} etc.)
since behavior is harder to predict in advance.}

\item{curl_handles}{A named list of curl handles. Each value is an
object from \code{curl::new_handle()}, and each name is a URL
(and should start with "http", "https", or "ftp").
Example:
list(
\verb{http://httpbin.org/basic-auth} = curl::new_handle(
username = "user", password = "passwd"
)
)
Then, if your plan has
\code{file_in("http://httpbin.org/basic-auth/user/passwd")}
\code{drake} will authenticate using the username and password of the handle
for \verb{http://httpbin.org/basic-auth/}.

\code{drake} uses partial matching on text to
find the right handle of the \code{file_in()} URL, so the name of the handle
could be the complete URL (\code{"http://httpbin.org/basic-auth/user/passwd"})
or a part of the URL (e.g. \code{"http://httpbin.org/"} or
\code{"http://httpbin.org/basic-auth/"}). If you have multiple handles
whose names match your URL, \code{drake} will choose the closest match.}

\item{max_expand}{Positive integer, optional.
\code{max_expand} is the maximum number of targets to generate in each
\code{map()}, \code{cross()}, or \code{group()} dynamic transform.
Useful if you have a massive number of dynamic sub-targets and you want to
work with only the first few sub-targets before scaling up.
Note: the \code{max_expand} argument of \code{make()} and
\code{drake_config()} is for dynamic branching only.
The static branching \code{max_expand}
is an argument of \code{drake_plan()} and \code{transform_plan()}.}

\item{log_build_times}{Logical, whether to record build_times for targets.
Mac users may notice a 20\% speedup in \code{make()}
with \code{build_times = FALSE}.}

\item{format}{Character, an optional custom storage format for targets
without an explicit \code{target(format = ...)} in the plan. Details
about formats:
\url{https://books.ropensci.org/drake/plans.html#special-data-formats-for-targets} # nolint}

\item{lock_cache}{Logical, whether to lock the cache before running \code{make()}
etc. It is usually recommended to keep cache locking on.
However, if you interrupt \code{make()} before it can clean itself up,
then the cache will stay locked,
and you will need to manually unlock it with
\code{drake::drake_cache("xyz")$unlock()}. Repeatedly unlocking the cache
by hand is annoying, and \code{lock_cache = FALSE} prevents the cache
from locking in the first place.}

\item{log_make}{Optional character scalar of a file name or
connection object (such as \code{stdout()}) to dump maximally verbose
log information for \code{\link[=make]{make()}} and other functions (all functions that
accept a \code{config} argument, plus \code{drake_config()}).
If you choose to use a text file as the console log,
it will persist over multiple function calls
until you delete it manually.
Fields in each row the log file, from left to right:
- The node name (short host name) of the
computer (from \code{Sys.info()["nodename"]}).
- The process ID (from \code{Sys.getpid()}).
- A timestamp with the date and time (in microseconds).
- A brief description of what \code{drake} was doing.\verb{ The fields are separated by pipe symbols (}"|"`).}

\item{log_worker}{Logical, same as the \code{log_worker} argument of
\code{clustermq::workers()} and \code{clustermq::Q()}. Only relevant
if \code{parallelism} is \code{"clustermq"}.}
}
\value{
A configured \code{drake} workflow.
}
\description{
Call this function inside the \verb{_drake.R}
script for \code{\link[=r_make]{r_make()}} and friends.
All non-deprecated function arguments are the same
between \code{\link[=make]{make()}} and \code{\link[=drake_config]{drake_config()}}.
}
\details{
In \code{drake}, \code{\link[=make]{make()}} has two stages:
\enumerate{
\item Configure a workflow to your environment and plan.
\item Build targets.
The \code{\link[=drake_config]{drake_config()}} function just does step (1),
which is a common requirement for not only \code{\link[=make]{make()}},
but also utility functions like \code{\link[=vis_drake_graph]{vis_drake_graph()}}
and \code{\link[=outdated]{outdated()}}. That is why \code{\link[=drake_config]{drake_config()}}
is a requirement for the \verb{_drake.R} script, which
powers \code{\link[=r_make]{r_make()}}, \code{\link[=r_outdated]{r_outdated()}}, \code{\link[=r_vis_drake_graph]{r_vis_drake_graph()}}, etc.
}
}
\section{Recovery}{

\code{make(recover = TRUE, recoverable = TRUE)}
powers automated data recovery.
The default of \code{recover} is \code{FALSE} because
targets recovered from the distant past may have been generated
with earlier versions of R and earlier package environments
that no longer exist.

How it works: if \code{recover} is \code{TRUE},
\code{drake} tries to salvage old target values from the cache
instead of running commands from the plan.
A target is recoverable if
\enumerate{
\item There is an old value somewhere in the cache that
shares the command, dependencies, etc.
of the target about to be built.
\item The old value was generated with \code{make(recoverable = TRUE)}.
}

If both conditions are met, \code{drake} will
\enumerate{
\item Assign the most recently-generated admissible data to the target, and
\item skip the target's command.
}
}

\examples{
\dontrun{
isolate_example("quarantine side effects", {
if (requireNamespace("knitr", quietly = TRUE)) {
writeLines(
  c(
    "library(drake)",
    "load_mtcars_example()",
    "drake_config(my_plan, targets = c(\"small\", \"large\"))"
  ),
  "_drake.R" # default value of the `source` argument
)
cat(readLines("_drake.R"), sep = "\n")
r_outdated()
r_make()
r_outdated()
}
})
}
}
\seealso{
\code{\link[=make]{make()}}, \code{\link[=drake_plan]{drake_plan()}}, \code{\link[=vis_drake_graph]{vis_drake_graph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_build.R
\name{drake_build_impl}
\alias{drake_build_impl}
\title{Internal function with a drake_config() argument}
\usage{
drake_build_impl(
  target,
  config = NULL,
  meta = NULL,
  character_only = FALSE,
  replace = FALSE
)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{evaluate}
\alias{evaluate}
\title{evaluate \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
evaluate(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{imported}
\alias{imported}
\title{List all the imports in the drake cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
imported(
  files_only = FALSE,
  path = getwd(),
  search = TRUE,
  cache = drake::get_cache(path = path, search = search, verbose = verbose),
  verbose = 1L,
  jobs = 1
)
}
\arguments{
\item{files_only}{Logical, whether to show imported files only
and ignore imported objects. Since all your functions and
all their global variables are imported, the full list of
imported objects could get really cumbersome.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{jobs}{Number of jobs/workers for parallel processing.}
}
\value{
Character vector naming the imports in the cache.
}
\description{
Deprecated on 2019-01-08.
}
\details{
An import is a non-target object processed
by \code{\link[=make]{make()}}. Targets in the workflow
plan data frame (see \code{\link[=drake_config]{drake_config()}}
may depend on imports.
}
\seealso{
\code{\link[=cached]{cached()}}, \code{\link[=loadd]{loadd()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{examples_drake}
\alias{examples_drake}
\title{examples_drake \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
examples_drake(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{drake_failed}
\alias{drake_failed}
\title{List failed targets.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_failed(cache = drake::drake_cache(path = path), path = NULL)
}
\arguments{
\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}
}
\value{
A character vector of target names.
}
\description{
List the targets that quit in error during \code{\link[=make]{make()}}.
}
\examples{
\dontrun{
isolate_example("contain side effects", {
if (suppressWarnings(require("knitr"))) {
# Build a plan doomed to fail:
bad_plan <- drake_plan(x = function_doesnt_exist())
cache <- storr::storr_environment() # optional
try(
  make(bad_plan, cache = cache, history = FALSE),
  silent = TRUE
) # error
drake_failed(cache = cache) # "x"
e <- diagnose(x, cache = cache) # Retrieve the cached error log of x.
names(e)
e$error
names(e$error)
}
})
}
}
\seealso{
\code{\link[=drake_done]{drake_done()}}, \code{\link[=drake_running]{drake_running()}}, \code{\link[=drake_cancelled]{drake_cancelled()}},
\code{\link[=drake_progress]{drake_progress()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{get_trace}
\alias{get_trace}
\title{Deprecated, get a trace of a dynamic target's value.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
get_trace(trace, value)
}
\arguments{
\item{trace}{Character, name of the trace
you want to extract. Such trace names are declared
in the \code{.trace} argument of \code{map()}, \code{cross()} or \code{group()}..}

\item{value}{Value of the dynamic target}
}
\value{
The dynamic trace of one target in another:
a vector of values from a grouping variable.
}
\description{
Deprecated on 2019-12-10. Use \code{\link[=read_trace]{read_trace()}} instead.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{plan_to_code}
\alias{plan_to_code}
\title{Turn a \code{drake} plan into a plain R script file.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#questioning}{\figure{lifecycle-questioning.svg}{options: alt='[Questioning]'}}}{\strong{[Questioning]}}}
\usage{
plan_to_code(plan, con = stdout())
}
\arguments{
\item{plan}{Workflow plan data frame. See \code{\link[=drake_plan]{drake_plan()}}
for details.}

\item{con}{A file path or connection to write to.}
}
\description{
\code{code_to_plan()}, \code{\link[=plan_to_code]{plan_to_code()}}, and
\code{\link[=plan_to_notebook]{plan_to_notebook()}} together illustrate the relationships
between \code{drake} plans, R scripts, and R Markdown documents.
In the file generated by \code{plan_to_code()}, every target/command pair
becomes a chunk of code.
Targets are arranged in topological order
so dependencies are available before their downstream targets.
Please note:
\enumerate{
\item You are still responsible for loading your project's
packages, imported functions, etc.
\item Triggers disappear.
}
}
\examples{
plan <- drake_plan(
  raw_data = read_excel(file_in("raw_data.xlsx")),
  data = raw_data,
  hist = create_plot(data),
  fit = lm(Ozone ~ Temp + Wind, data)
)
file <- tempfile()
# Turn the plan into an R script a the given file path.
plan_to_code(plan, file)
# Here is what the script looks like.
cat(readLines(file), sep = "\n")
# Convert back to a drake plan.
code_to_plan(file)
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}, \code{\link[=code_to_plan]{code_to_plan()}},
\code{\link[=plan_to_notebook]{plan_to_notebook()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{render_graph}
\alias{render_graph}
\title{render_graph \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
render_graph(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{running}
\alias{running}
\title{List running targets.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
running(
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  verbose = 1L
)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\value{
A character vector of target names.
}
\description{
Deprecated on 2020-03-23. Use \code{\link[=drake_running]{drake_running()}} instead.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{deprecate_wildcard}
\alias{deprecate_wildcard}
\title{deprecate_wildcard \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
deprecate_wildcard(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{drake_plan_source}
\alias{drake_plan_source}
\title{Show the code required to produce a given \code{drake} plan
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_plan_source(plan)
}
\arguments{
\item{plan}{A workflow plan data frame (see \code{\link[=drake_plan]{drake_plan()}})}
}
\value{
a character vector of lines of text. This text
is a call to \code{\link[=drake_plan]{drake_plan()}} that produces the plan you provide.
}
\description{
You supply a plan, and \code{\link[=drake_plan_source]{drake_plan_source()}}
supplies code to generate that plan. If you have the
\href{https://github.com/r-lib/prettycode}{\code{prettycode} package},
installed, you also get nice syntax highlighting in the console
when you print it.
}
\examples{
plan <- drake::drake_plan(
  small_data = download_data("https://some_website.com"),
  large_data_raw = target(
    command = download_data("https://lots_of_data.com"),
    trigger = trigger(
      change = time_last_modified("https://lots_of_data.com"),
      command = FALSE,
      depend = FALSE
    ),
    timeout = 1e3
  )
)
print(plan)
if (requireNamespace("styler", quietly = TRUE)) {
  source <- drake_plan_source(plan)
  print(source) # Install the prettycode package for syntax highlighting.
  file <- tempfile() # Path to an R script to contain the drake_plan() call.
  writeLines(source, file) # Save the code to an R script.
}
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{deps_targets}
\alias{deps_targets}
\title{See the dependencies of a target
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
deps_targets(targets, config, reverse = FALSE)
}
\arguments{
\item{targets}{A character vector of target names.}

\item{config}{An output list from \code{\link[=drake_config]{drake_config()}}}

\item{reverse}{Logical, whether to compute reverse dependencies
(targets immediately downstream) instead of ordinary dependencies.}
}
\value{
Names of dependencies listed by type (object, input file, etc).
}
\description{
Use \code{\link[=deps_target]{deps_target()}} (singular) instead.
}
\details{
Deprecated on 2018-08-30.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.R
\name{rescue_cache}
\alias{rescue_cache}
\title{Try to repair a drake cache that is prone
to throwing \code{storr}-related errors.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#questioning}{\figure{lifecycle-questioning.svg}{options: alt='[Questioning]'}}}{\strong{[Questioning]}}}
\usage{
rescue_cache(
  targets = NULL,
  path = NULL,
  search = NULL,
  verbose = NULL,
  force = FALSE,
  cache = drake::drake_cache(path = path),
  jobs = 1,
  garbage_collection = FALSE
)
}
\arguments{
\item{targets}{Character vector, names of the targets to rescue.
As with many other drake utility functions, the word \code{target}
is defined generally in this case, encompassing imports
as well as true targets.
If \code{targets} is \code{NULL}, everything in the
cache is rescued.}

\item{path}{Character.
Set \code{path} to the path of a \code{storr::storr_rds()} cache
to retrieve a specific cache generated by \code{storr::storr_rds()}
or \code{drake::new_cache()}. If the \code{path} argument is \code{NULL},
\code{drake_cache()} searches up through parent directories
to find a folder called \verb{.drake/}.}

\item{search}{Deprecated.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{force}{Deprecated.}

\item{cache}{A \code{storr} cache object.}

\item{jobs}{Number of jobs for light parallelism
(disabled on Windows).}

\item{garbage_collection}{Logical, whether to do garbage collection
as a final step. See \code{\link[=drake_gc]{drake_gc()}} and \code{\link[=clean]{clean()}}
for details.}
}
\value{
Nothing.
}
\description{
Sometimes, \code{storr} caches may have
dangling orphaned files that prevent you from loading or cleaning.
This function tries to remove those files so you can use the
cache normally again.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build targets. This creates the cache.
# Remove dangling cache files that could cause errors.
rescue_cache(jobs = 2)
# Alternatively, just rescue targets 'small' and 'large'.
# Rescuing specific targets is usually faster.
rescue_cache(targets = c("small", "large"))
}
})
}
}
\seealso{
\code{\link[=drake_cache]{drake_cache()}}, \code{\link[=cached]{cached()}},
\code{\link[=drake_gc]{drake_gc()}}, \code{\link[=clean]{clean()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{ignore}
\alias{ignore}
\title{Ignore code
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
ignore(x = NULL)
}
\arguments{
\item{x}{Code to ignore.}
}
\value{
The argument.
}
\description{
Ignore sections of commands and imported functions.
}
\details{
In user-defined functions and \code{\link[=drake_plan]{drake_plan()}} commands, you can
wrap code chunks in \code{ignore()} to
\enumerate{
\item Tell \code{drake} to not search for dependencies
(targets etc. mentioned in the code) and
\item Ignore changes to the code so downstream targets remain up to date.
To enforce (1) without (2), use \code{\link[=no_deps]{no_deps()}}.
}
}
\section{Keywords}{

\code{\link[=drake_plan]{drake_plan()}} understands special keyword functions for your commands.
With the exception of \code{\link[=target]{target()}}, each one is a proper function
with its own help file.
\itemize{
\item \code{\link[=target]{target()}}: give the target more than just a command.
Using \code{\link[=target]{target()}}, you can apply a transformation
(examples: \url{https://books.ropensci.org/drake/plans.html#large-plans}), # nolint
supply a trigger (\url{https://books.ropensci.org/drake/triggers.html}), # nolint
or set any number of custom columns.
\item \code{\link[=file_in]{file_in()}}: declare an input file dependency.
\item \code{\link[=file_out]{file_out()}}: declare an output file to be produced
when the target is built.
\item \code{\link[=knitr_in]{knitr_in()}}: declare a \code{knitr} file dependency such as an
R Markdown (\verb{*.Rmd}) or R LaTeX (\verb{*.Rnw}) file.
\item \code{\link[=ignore]{ignore()}}: force \code{drake} to entirely ignore a piece of code:
do not track it for changes and do not analyze it for dependencies.
\item \code{\link[=no_deps]{no_deps()}}: tell \code{drake} to not track the dependencies
of a piece of code. \code{drake} still tracks the code itself for changes.
\item \code{\link[=id_chr]{id_chr()}}: Get the name of the current target.
\item \code{\link[=drake_envir]{drake_envir()}}: get the environment where drake builds targets.
Intended for advanced custom memory management.
}
}

\examples{
\dontrun{
isolate_example("Contain side effects", {
# Normally, `drake` reacts to changes in dependencies.
x <- 4
make(plan = drake_plan(y = sqrt(x)))
x <- 5
make(plan = drake_plan(y = sqrt(x)))
make(plan = drake_plan(y = sqrt(4) + x))
# But not with ignore().
make(plan = drake_plan(y = sqrt(4) + ignore(x))) # Builds y.
x <- 6
make(plan = drake_plan(y = sqrt(4) + ignore(x))) # Skips y.
make(plan = drake_plan(y = sqrt(4) + ignore(x + 1))) # Skips y.

# ignore() works with functions and multiline code chunks.
f <- function(x) {
  ignore({
    x <- x + 1
    x <- x + 2
  })
  x # Not ignored.
}
make(plan = drake_plan(y = f(2)))
readd(x)
# Changes the content of the ignore() block:
f <- function(x) {
  ignore({
    x <- x + 1
  })
  x # Not ignored.
}
make(plan = drake_plan(x = f(2)))
readd(x)
})
}
}
\seealso{
\code{\link[=file_in]{file_in()}}, \code{\link[=file_out]{file_out()}}, \code{\link[=knitr_in]{knitr_in()}}, \code{\link[=no_deps]{no_deps()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outdated.R
\name{outdated}
\alias{outdated}
\title{List the targets that are out of date.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
outdated(..., make_imports = TRUE, do_prework = TRUE, config = NULL)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets} and \code{envir}.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to save some time and risk obsolete output.}

\item{do_prework}{Whether to do the \code{prework}
normally supplied to \code{\link[=make]{make()}}.}

\item{config}{Deprecated (2019-12-21).
A configured workflow from \code{\link[=drake_config]{drake_config()}}.}
}
\value{
Character vector of the names of outdated targets.
}
\description{
Outdated targets will be rebuilt in the next
\code{\link[=make]{make()}}. \code{outdated()} does not show dynamic sub-targets.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
# Recopute the config list early and often to have the
# most current information. Do not modify the config list by hand.
outdated(my_plan) # Which targets are out of date?
make(my_plan) # Run the projects, build the targets.
# Now, everything should be up to date (no targets listed).
outdated(my_plan)
}
})
}
}
\seealso{
\code{\link[=r_outdated]{r_outdated()}}, \code{\link[=drake_config]{drake_config()}}, \code{\link[=missed]{missed()}}, \code{\link[=drake_plan]{drake_plan()}},
\code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deps.R
\name{deps_profile}
\alias{deps_profile}
\title{Find out why a target is out of date.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
deps_profile(target, ..., character_only = FALSE, config = NULL)
}
\arguments{
\item{target}{Name of the target.}

\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{character_only}{Logical, whether to assume \code{target}
is a character string rather than a symbol.}

\item{config}{Deprecated.}
}
\value{
A data frame of old and new values for each
of the main triggers, along with
an indication of which values changed since
the last \code{\link[=make]{make()}}.
}
\description{
The dependency profile can give you
a hint as to why a target is out of date.
It can tell you if
\itemize{
\item the command changed
(\code{\link[=deps_profile]{deps_profile()}} reports the \emph{hash} of the command,
not the command itself)
\item at least one input file changed,
\item at least one output file changed,
\item or a non-file dependency changed. For this last part,
the imports need to be up to date in the cache,
which you can do with \code{outdated()} or
\code{make(skip_targets = TRUE)}.
\item the pseudo-random number generator seed changed.
Unfortunately, \code{deps_profile()} does not
currently get more specific than that.
}
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Load drake's canonical example.
make(my_plan) # Run the project, build the targets.
# Get some example dependency profiles of targets.
deps_profile(small, my_plan)
# Change a dependency.
simulate <- function(x) {}
# Update the in-memory imports in the cache
# so deps_profile can detect changes to them.
# Changes to targets are already cached.
make(my_plan, skip_targets = TRUE)
# The dependency hash changed.
deps_profile(small, my_plan)
}
})
}
}
\seealso{
\code{\link[=diagnose]{diagnose()}},
\code{\link[=deps_code]{deps_code()}}, \code{\link[=make]{make()}},
\code{\link[=drake_config]{drake_config()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{drake_progress}
\alias{drake_progress}
\title{Get the build progress of your targets
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_progress(
  ...,
  list = character(0),
  cache = drake::drake_cache(path = path),
  path = NULL,
  progress = NULL
)
}
\arguments{
\item{...}{Objects to load from the cache, as names (unquoted)
or character strings (quoted). If the \code{tidyselect} package is installed,
you can also supply \code{dplyr}-style \code{tidyselect}
commands such as \code{starts_with()}, \code{ends_with()}, and \code{one_of()}.}

\item{list}{Character vector naming objects to be loaded from the
cache. Similar to the \code{list} argument of \code{\link[=remove]{remove()}}.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{progress}{Character vector for filtering the build progress results.
Defaults to \code{NULL} (no filtering) to report progress of all objects.
Supported filters are \code{"done"}, \code{"running"}, and \code{"failed"}.}
}
\value{
The build progress of each target reached by
the current \code{\link[=make]{make()}} so far.
}
\description{
Objects that drake imported, built, or attempted
to build are listed as \code{"done"} or \code{"running"}.
Skipped objects are not listed.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the targets.
# Watch the changing drake_progress() as make() is running.
drake_progress() # List all the targets reached so far.
drake_progress(small, large) # Just see the progress of some targets.
drake_progress(list = c("small", "large")) # Same as above.
}
})
}
}
\seealso{
\code{\link[=diagnose]{diagnose()}}, \code{\link[=drake_get_session_info]{drake_get_session_info()}},
\code{\link[=cached]{cached()}}, \code{\link[=readd]{readd()}}, \code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{migrate_drake_project}
\alias{migrate_drake_project}
\title{migrate_drake_project \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
migrate_drake_project(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-05-16
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{Makefile_recipe}
\alias{Makefile_recipe}
\title{Default Makefile recipe
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
Makefile_recipe(
  recipe_command = drake::default_recipe_command(),
  target = "your_target",
  cache_path = NULL
)
}
\arguments{
\item{recipe_command}{Character scalar.}

\item{target}{Character scalar.}

\item{cache_path}{Character scalar.}
}
\value{
A character scalar
}
\description{
2019-01-03
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.R
\name{clean}
\alias{clean}
\title{Invalidate and deregister targets.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
clean(
  ...,
  list = character(0),
  destroy = FALSE,
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  verbose = NULL,
  jobs = NULL,
  force = FALSE,
  garbage_collection = FALSE,
  purge = FALSE
)
}
\arguments{
\item{...}{Symbols, individual targets to remove.}

\item{list}{Character vector of individual targets to remove.}

\item{destroy}{Logical, whether to totally remove the drake cache.
If \code{destroy} is \code{FALSE}, only the targets
from \code{make()}
are removed. If \code{TRUE}, the whole cache is removed, including
session metadata, etc.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated}

\item{jobs}{Deprecated.}

\item{force}{Logical, whether to try to clean the cache
even though the project may not be back compatible with the
current version of drake.}

\item{garbage_collection}{Logical, whether to call
\code{cache$gc()} to do garbage collection.
If \code{TRUE}, cached data with no remaining references
will be removed.
This will slow down \code{clean()}, but the cache
could take up far less space afterwards.
See the \code{gc()} method for \code{storr} caches.}

\item{purge}{Logical, whether to remove objects from
metadata namespaces such as "meta", "build_times", and "errors".}
}
\value{
Invisibly return \code{NULL}.
}
\description{
Force targets to be out of date and remove target names
from the data in the cache. Be careful and run \code{\link[=which_clean]{which_clean()}} before
\code{\link[=clean]{clean()}}. That way, you know beforehand which targets will be
compromised.
}
\details{
By default, \code{\link[=clean]{clean()}} invalidates \strong{all} targets,
so be careful. \code{\link[=clean]{clean()}} always:
\enumerate{
\item Forces targets to be out of date so the next \code{\link[=make]{make()}}
does not skip them.
\item Deregisters targets so \code{loadd(your_target)} and \code{readd(your_target)}
no longer work.
}

By default, \code{clean()} does not actually remove the underlying data.
Even old targets from the distant past are still in the cache
and recoverable via \code{drake_history()} and \code{make(recover = TRUE)}.
To actually remove target data from the cache, as well as any
\code{\link[=file_out]{file_out()}} files from any targets you are currently cleaning,
run \code{clean(garbage_collection = TRUE)}.
Garbage collection is slow, but it reduces the storage burden of the cache.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
make(my_plan) # Run the project, build the targets.
# Show all registered targets in the cache.
cached()
# Deregister 'summ_regression1_large' and 'small' in the cache.
clean(summ_regression1_large, small)
# Those objects are no longer registered as targets.
cached()
# Rebuild the invalidated/outdated targets.
make(my_plan)
# Clean everything.
clean()
# But the data objects and files are not actually gone!
file.exists("report.md")
drake_history()
make(my_plan, recover = TRUE)
# You need garbage collection to actually remove the data
# and any file_out() files of any uncleaned targets.
clean(garbage_collection = TRUE)
drake_history()
make(my_plan, recover = TRUE)
}
})
}
}
\seealso{
\code{\link[=which_clean]{which_clean()}}, \code{\link[=drake_gc]{drake_gc()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.R
\name{which_clean}
\alias{which_clean}
\title{Which targets will \code{clean()} invalidate?
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
which_clean(
  ...,
  list = character(0),
  path = NULL,
  cache = drake::drake_cache(path = path)
)
}
\arguments{
\item{...}{Targets to remove from the cache: as names (symbols) or
character strings. If the \code{tidyselect} package is installed,
you can also supply \code{dplyr}-style \code{tidyselect}
commands such as \code{starts_with()}, \code{ends_with()}, and \code{one_of()}.}

\item{list}{Character vector naming targets to be removed from the
cache. Similar to the \code{list} argument of \code{\link[=remove]{remove()}}.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}
}
\description{
\code{which_clean()} is a safety check for \code{clean()}.
It shows you the targets that \code{clean()} will
invalidate (or remove if \code{garbage_collection} is \code{TRUE}).
It helps you avoid accidentally removing targets you care about.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
plan <- drake_plan(x = 1, y = 2, z = 3)
make(plan)
cached()
which_clean(x, y) # [1] "x" "y"
clean(x, y)       # Invalidates targets x and y.
cached()          # [1] "z"
})
}
}
\seealso{
\code{\link[=clean]{clean()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manage_memory.R
\name{manage_memory}
\alias{manage_memory}
\title{Manage the in-memory dependencies of a target.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
manage_memory(target, config, downstream = NULL, jobs = 1)
}
\arguments{
\item{target}{Character, name of the target.}

\item{config}{\code{\link[=drake_config]{drake_config()}} list.}

\item{downstream}{Optional, character vector of any targets
assumed to be downstream.}

\item{jobs}{Number of jobs for local parallel computing}
}
\value{
Nothing.
}
\description{
Load/unload a target's dependencies.
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{cached}
\alias{cached}
\title{List targets in the cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
cached(
  ...,
  list = character(0),
  no_imported_objects = FALSE,
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  verbose = NULL,
  namespace = NULL,
  jobs = 1,
  targets_only = TRUE
)
}
\arguments{
\item{...}{Deprecated. Do not use.
Objects to load from the cache, as names (unquoted)
or character strings (quoted). Similar to \code{...} in
\code{remove()}.}

\item{list}{Deprecated. Do not use.
Character vector naming objects to be loaded from the
cache. Similar to the \code{list} argument of \code{\link[=remove]{remove()}}.}

\item{no_imported_objects}{Logical, deprecated. Use
\code{targets_only} instead.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{namespace}{Character scalar, name of the storr namespace
to use for listing objects.}

\item{jobs}{Number of jobs/workers for parallel processing.}

\item{targets_only}{Logical. If \code{TRUE} just list the targets.
If \code{FALSE}, list files and imported objects too.}
}
\value{
Either a named logical indicating whether the given
targets or cached or a character vector listing all cached
items, depending on whether any targets are specified.
}
\description{
Tip: read/load a cached item with \code{\link[=readd]{readd()}}
or \code{\link[=loadd]{loadd()}}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
if (requireNamespace("lubridate")) {
load_mtcars_example() # Load drake's canonical example.
make(my_plan) # Run the project, build all the targets.
cached()
cached(targets_only = FALSE)
}
}
})
}
}
\seealso{
\code{\link[=cached_planned]{cached_planned()}}, \code{\link[=cached_unplanned]{cached_unplanned()}},
\code{\link[=readd]{readd()}}, \code{\link[=loadd]{loadd()}},
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{drake_strings}
\alias{drake_strings}
\title{Turn valid expressions into character strings.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
drake_strings(...)
}
\arguments{
\item{...}{Unquoted symbols to turn into character strings.}
}
\value{
A character vector.
}
\description{
Deprecated on 2019-01-01
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{code_to_function}
\alias{code_to_function}
\title{Turn a script into a function.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
code_to_function(path, envir = parent.frame())
}
\arguments{
\item{path}{Character vector, path to script.}

\item{envir}{Environment of the created function.}
}
\value{
A function to be input into the drake plan
}
\description{
\code{code_to_function()} is a quick (and very dirty) way to
retrofit drake to an existing script-based project. It parses
individual \verb{\\*.R/\\*.RMD} files into functions so they can be added
into the drake workflow.
}
\details{
Most data science workflows consist of imperative scripts.
\code{drake}, on the other hand, assumes you write \emph{functions}.
\code{code_to_function()} allows for pre-existing workflows to incorporate
drake as a workflow management tool seamlessly for cases where
re-factoring is unfeasible. So drake can monitor dependencies, the
targets are passed as arguments of the dependent functions.
}
\examples{
\dontrun{
isolate_example("contain side effects", {
if (requireNamespace("ggplot2", quietly = TRUE)) {
# The `code_to_function()` function creates a function that makes it
# available for drake to process as part of the workflow.
# The main purpose is to allow pre-existing workflows to incorporate drake
# into the workflow seamlessly for cases where re-factoring is unfeasible.
#

script1 <- tempfile()
script2 <- tempfile()
script3 <- tempfile()
script4 <- tempfile()

writeLines(c(
  "data <- mtcars",
  "data$make <- do.call('c',",
  "lapply(strsplit(rownames(data), split=\" \"), `[`, 1))",
  "saveRDS(data, \"mtcars_alt.RDS\")"
 ),
  script1
)

writeLines(c(
  "data <- readRDS(\"mtcars_alt.RDS\")",
  "mtcars_lm <- lm(mpg~cyl+disp+vs+gear+make,data=data)",
  "saveRDS(mtcars_lm, \"mtcars_lm.RDS\")"
  ),
  script2
)
writeLines(c(
  "mtcars_lm <- readRDS(\"mtcars_lm.RDS\")",
  "lm_summary <- summary(mtcars_lm)",
  "saveRDS(lm_summary, \"mtcars_lm_summary.RDS\")"
  ),
  script3
)
writeLines(c(
  "data<-readRDS(\"mtcars_alt.RDS\")",
  "gg <- ggplot2::ggplot(data)+",
  "ggplot2::geom_point(ggplot2::aes(",
  "x=disp, y=mpg, shape=as.factor(vs), color=make))",
  "ggplot2::ggsave(\"mtcars_plot.png\", gg)"
 ),
  script4
)


do_munge <- code_to_function(script1)
do_analysis <- code_to_function(script2)
do_summarize <- code_to_function(script3)
do_vis <- code_to_function(script4)

plan <- drake_plan(
  munged   = do_munge(),
  analysis = do_analysis(munged),
  summary  = do_summarize(analysis),
  plot     = do_vis(munged)
 )

plan
# drake knows  "script1" is the first script to be evaluated and ran,
# because it has no dependencies on other code and a dependency of
# `analysis`. See for yourself:

make(plan)

# See the connections that the sourced scripts create:
if (requireNamespace("visNetwork", quietly = TRUE)) {
  vis_drake_graph(plan)
}
}
})
}
}
\seealso{
\code{\link[=file_in]{file_in()}}, \code{\link[=file_out]{file_out()}}, \code{\link[=knitr_in]{knitr_in()}}, \code{\link[=ignore]{ignore()}}, \code{\link[=no_deps]{no_deps()}},
\code{\link[=code_to_plan]{code_to_plan()}}, \code{\link[=plan_to_code]{plan_to_code()}}, \code{\link[=plan_to_notebook]{plan_to_notebook()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outdated.R
\name{missed}
\alias{missed}
\title{Report any import objects required by your drake_plan
plan but missing from your workspace or file system.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
missed(..., config = NULL)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{config}{Deprecated.}
}
\value{
Character vector of names of missing objects and files.
}
\description{
Checks your workspace/environment and
file system.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
plan <- drake_plan(x = missing::fun(arg))
missed(plan)
}
})
}
}
\seealso{
\code{\link[=outdated]{outdated()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{cached_unplanned}
\alias{cached_unplanned}
\title{List targets in the cache but not the plan.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
cached_unplanned(
  plan,
  path = NULL,
  cache = drake::drake_cache(path = path),
  namespace = NULL,
  jobs = 1
)
}
\arguments{
\item{plan}{A drake plan.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{namespace}{Character scalar, name of the storr namespace
to use for listing objects.}

\item{jobs}{Number of jobs/workers for parallel processing.}
}
\value{
A character vector of target and sub-target names.
}
\description{
Includes dynamic sub-targets as well.
See examples for details.
}
\examples{
\dontrun{
isolate_example("cache_unplanned() example", {
plan <- drake_plan(w = 1)
make(plan)
cached_unplanned(plan)
plan <- drake_plan(
  x = seq_len(2),
  y = target(x, dynamic = map(x))
)
cached_unplanned(plan)
make(plan)
cached_unplanned(plan)
# cached_unplanned() helps clean superfluous targets.
cached()
clean(list = cached_unplanned(plan))
cached()
})
}
}
\seealso{
\code{\link[=cached]{cached()}}, \link{cached_planned}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{shell_file}
\alias{shell_file}
\title{Shell file for Makefile parallelism
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
shell_file(path = "shell.sh", overwrite = FALSE)
}
\arguments{
\item{path}{Character.}

\item{overwrite}{Logical.}
}
\value{
logical
}
\description{
2019-01-03
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hpc.R
\name{drake_hpc_template_files}
\alias{drake_hpc_template_files}
\title{List the available example template files for deploying
work to a cluster / job scheduler.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_hpc_template_files()
}
\value{
A character vector of example template files that
you can write with \code{\link[=drake_hpc_template_file]{drake_hpc_template_file()}}.
}
\description{
See the example files from
\code{\link[=drake_examples]{drake_examples()}} and \code{\link[=drake_example]{drake_example()}}
for example usage.
}
\examples{
\dontrun{
plan <- drake_plan(x = rnorm(1e7), y = rnorm(1e7))
# List the available template files.
drake_hpc_template_files()
# Write a SLURM template file.
out <- file.path(tempdir(), "slurm_batchtools.tmpl")
drake_hpc_template_file("slurm_batchtools.tmpl", to = tempdir())
cat(readLines(out), sep = "\n")
# library(future.batchtools) # nolint
# future::plan(batchtools_slurm, template = out) # nolint
# make(plan, parallelism = "future", jobs = 2) # nolint
}
}
\seealso{
\code{\link[=drake_hpc_template_file]{drake_hpc_template_file()}},
\code{\link[=drake_examples]{drake_examples()}}, \code{\link[=drake_example]{drake_example()}},
\code{\link[=shell_file]{shell_file()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{expand}
\alias{expand}
\title{expand \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
expand(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{drake_tip}
\alias{drake_tip}
\title{Output a random tip about drake.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
drake_tip()
}
\value{
A character scalar with a tip on how to use drake.
}
\description{
Deprecated on 2019-01-12.
}
\details{
Tips are usually related to news and usage.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio.R
\name{rs_addin_r_outdated}
\alias{rs_addin_r_outdated}
\title{RStudio addin for r_outdated()
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
rs_addin_r_outdated(r_args = list(), .print = TRUE)
}
\arguments{
\item{r_args}{List of arguments to \code{r_fn}, not including \code{func} or \code{args}.
Example:
\code{r_make(r_fn = callr::r_bg, r_args = list(stdout = "stdout.log"))}.}

\item{.print}{Logical, whether to \code{print()} the result
to the console. Required for the addin.}
}
\value{
A character vector of outdated targets.
}
\description{
Call \code{\link[=r_outdated]{r_outdated()}} in an RStudio addin.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_plan.R
\name{drake_slice}
\alias{drake_slice}
\title{Take a strategic subset of a dataset.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_slice(data, slices, index, margin = 1L, drop = FALSE)
}
\arguments{
\item{data}{A list, vector, data frame, matrix, or arbitrary array.
Anything with a \code{length()} or \code{dim()}.}

\item{slices}{Integer of length 1, number of slices (i.e. pieces)
of the whole dataset. Remember, \code{drake_slice(index = i)} returns
only slice number \code{i}.}

\item{index}{Integer of length 1, which piece of the partition to return.}

\item{margin}{Integer of length 1, margin over which to split the data.
For example, for a data frame or matrix,
use \code{margin = 1} to split over rows and \code{margin = 2}
to split over columns. Similar to \code{MARGIN} in \code{apply()}.}

\item{drop}{Logical, for matrices and arrays.
If \code{TRUE},\verb{ the result is coerced to the lowest possible dimension. See ?}[` for details.}
}
\value{
A subset of \code{data}.
}
\description{
\code{drake_slice()} is similar to \code{split()}.
Both functions partition data into disjoint subsets,
but whereas \code{split()} returns \emph{all} the subsets, \code{drake_slice()}
returns just \emph{one}. In other words, \code{drake_slice(..., index = i)}
returns \code{split(...)[[i]]}.
Other features:
1. \code{drake_slice()} works on vectors, data frames,
matrices, lists, and arbitrary arrays.
2. Like \code{parallel::splitIndices()}, \code{drake_slice()} tries to
distribute the data uniformly across subsets.
See the examples to learn why splitting is useful in \code{drake}.
}
\examples{
# Simple usage
x <- matrix(seq_len(20), nrow = 5)
x
drake_slice(x, slices = 3, index = 1)
drake_slice(x, slices = 3, index = 2)
drake_slice(x, slices = 3, index = 3)
drake_slice(x, slices = 3, margin = 2, index = 1)
# In drake, you can split a large dataset over multiple targets.
\dontrun{
isolate_example("contain side effects", {
plan <- drake_plan(
  large_data = mtcars,
  data_split = target(
    drake_slice(large_data, slices = 32, index = i),
    transform = map(i = !!seq_len(32))
  )
)
plan
cache <- storr::storr_environment()
make(plan, cache = cache, session_info = FALSE, verbose = FALSE)
readd(data_split_1L, cache = cache)
readd(data_split_2L, cache = cache)
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_history.R
\name{drake_history}
\alias{drake_history}
\title{History and provenance
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_history(cache = NULL, history = NULL, analyze = TRUE, verbose = NULL)
}
\arguments{
\item{cache}{drake cache as created by \code{\link[=new_cache]{new_cache()}}.
See also \code{\link[=drake_cache]{drake_cache()}}.}

\item{history}{Logical, whether to record the build history
of your targets. You can also supply a
\href{https://github.com/wlandau/txtq}{\code{txtq}}, which is
how \code{drake} records history.
Must be \code{TRUE} for \code{\link[=drake_history]{drake_history()}} to work later.}

\item{analyze}{Logical, whether to analyze \code{\link[=drake_plan]{drake_plan()}}
commands for arguments to function calls.
Could be slow because this requires parsing and analyzing
lots of R code.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\value{
A data frame of target history.
}
\description{
See the history and provenance of your targets:
what you ran, when you ran it, the function arguments
you used, and how to get old data back.
}
\details{
\code{\link[=drake_history]{drake_history()}} returns a data frame with the following columns.
\itemize{
\item \code{target}: the name of the target.
\item \code{current}: logical, whether the row describes the data
actually assigned to the target name in the cache,
e.g. what you get with \code{loadd(target)} and \code{readd(target)}.
Does \strong{NOT} tell you if the target is up to date.
\item \code{built}: when the target's value was stored in the cache.
This is the true creation date of the target's value,
not the recovery date from \code{make(recover = TRUE)}.
\item \code{exists}: logical, whether the target's historical value
still exists in the cache. Garbage collection via
(\code{clean(garbage_collection = TRUE)} and \code{drake_cache()$gc()})
remove these historical values, but \code{clean()} under the default
settings does not.
\item \code{hash}: fingerprint of the target's historical value in the cache.
If the value still exists, you can read it with
\code{drake_cache()$get_value(hash)}.
\item \code{command}: the \code{\link[=drake_plan]{drake_plan()}} command executed to build the target.
\item \code{seed}: random number generator seed.
\item \code{runtime}: the time it took to execute the \code{\link[=drake_plan]{drake_plan()}} command.
Does not include overhead due to \code{drake}'s processing.
}

If \code{analyze} is \code{TRUE}, various other columns are included to show
the explicitly-named length-1 arguments to function calls in the commands.
See the "Provenance" section for more details.
}
\section{Provenance}{

If \code{analyze} is \code{TRUE}, \code{drake}
scans your \code{\link[=drake_plan]{drake_plan()}} commands
for function arguments and mentions them in the history.
A function argument shows up if and only if:
1. It has length 1. \cr
2. It is atomic, i.e. a base type: logical, integer,
real, complex, character, or raw. \cr
3. It is explicitly named in the function call,
For example, \code{x} is detected as \code{1} in
\code{fn(list(x = 1))} but not \code{f(list(1))}.
The exceptions are \code{\link[=file_out]{file_out()}}, \code{\link[=file_in]{file_in()}},
and \code{\link[=knitr_in]{knitr_in()}}. For example, \code{filename} is detected
as \code{"my_file.csv"} in
\code{process_data(filename = file_in("my_file.csv"))}.
NB: in \code{process_data(filename = file_in("a", "b"))}
\code{filename} is not detected because the value must be atomic. \cr
}

\examples{
\dontrun{
isolate_example("contain side-effects", {
if (requireNamespace("knitr", quietly = TRUE)) {
# First, let's iterate on a drake workflow.
load_mtcars_example()
make(my_plan, history = TRUE, verbose = 0L)
# Naturally, we'll make updates to our targets along the way.
reg2 <- function(d) {
  d$x2 <- d$x ^ 3
  lm(y ~ x2, data = d)
}
Sys.sleep(0.01)
make(my_plan, history = TRUE, verbose = 0L)
# The history is a data frame about all the recorded runs of your targets.
out <- drake_history(analyze = TRUE)
print(out)
# Let's use the history to recover the oldest version
# of our regression2_small target.
oldest_reg2_small <- max(which(out$target == "regression2_small"))
hash_oldest_reg2_small <- out[oldest_reg2_small, ]$hash
cache <- drake_cache()
cache$get_value(hash_oldest_reg2_small)
# If you run clean(), drake can still find all the targets.
clean(small)
drake_history()
# But if you run clean() with garbage collection,
# older versions of your targets may be gone.
clean(large, garbage_collection = TRUE)
drake_history()
invisible()
}
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deps.R
\name{deps_code}
\alias{deps_code}
\title{List the dependencies of a function or command
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
deps_code(x)
}
\arguments{
\item{x}{A function, expression, or text.}
}
\value{
A data frame of the dependencies.
}
\description{
Functions are assumed to be imported,
and language/text are assumed to be commands in a plan.
}
\examples{
# Your workflow likely depends on functions in your workspace.
f <- function(x, y) {
  out <- x + y + g(x)
  saveRDS(out, "out.rds")
}
# Find the dependencies of f. These could be R objects/functions
# in your workspace or packages. Any file names or target names
# will be ignored.
deps_code(f)
# Define a workflow plan data frame that uses your function f().
my_plan <- drake_plan(
  x = 1 + some_object,
  my_target = x + readRDS(file_in("tracked_input_file.rds")),
  return_value = f(x, y, g(z + w))
)
# Get the dependencies of workflow plan commands.
# Here, the dependencies could be R functions/objects from your workspace
# or packages, imported files, or other targets in the workflow plan.
deps_code(my_plan$command[[1]])
deps_code(my_plan$command[[2]])
deps_code(my_plan$command[[3]])
# You can also supply expressions or text.
deps_code(quote(x + y + 123))
deps_code("x + y + 123")
}
\seealso{
\code{\link[=deps_target]{deps_target()}}, \code{\link[=deps_knitr]{deps_knitr()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{is_function_call}
\alias{is_function_call}
\title{is_function_call \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
is_function_call(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_example.R
\name{load_mtcars_example}
\alias{load_mtcars_example}
\title{Load the mtcars example.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
load_mtcars_example(
  envir = parent.frame(),
  report_file = NULL,
  overwrite = FALSE,
  force = FALSE
)
}
\arguments{
\item{envir}{The environment to load the example into.
Defaults to your workspace.
For an insulated workspace,
set \code{envir = new.env(parent = globalenv())}.}

\item{report_file}{Where to write the report file. Deprecated.
In a future release, the report file will always be
\code{report.Rmd} and will always be written to your
working directory (current default).}

\item{overwrite}{Logical, whether to overwrite an
existing file \code{report.Rmd}.}

\item{force}{Deprecated.}
}
\value{
Nothing.
}
\description{
Is there an association between
the weight and the fuel efficiency of cars?
To find out, we use the mtcars example from \code{drake_example("mtcars")}.
The mtcars dataset itself only has 32 rows,
so we generate two larger bootstrapped datasets
and then analyze them with regression models.
Finally, we summarize the regression models
to see if there is an association.
}
\details{
Use \code{drake_example("mtcars")} to get the code
for the mtcars example.
This function also writes/overwrites
the file, \code{report.Rmd}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
# Populate your workspace and write 'report.Rmd'.
load_mtcars_example() # Get the code: drake_example("mtcars")
# Check the dependencies of an imported function.
deps_code(reg1)
# Check the dependencies of commands in the workflow plan.
deps_code(my_plan$command[1])
deps_code(my_plan$command[4])
# Plot the interactive network visualization of the workflow.
outdated(my_plan) # Which targets are out of date?
# Run the workflow to build all the targets in the plan.
make(my_plan)
outdated(my_plan) # Everything should be up to date.
# For the reg2() model on the small dataset,
# the p-value is so small that there may be an association
# between weight and fuel efficiency after all.
readd(coef_regression2_small)
# Clean up the example.
clean_mtcars_example()
}
})
}
}
\seealso{
\code{\link[=clean_mtcars_example]{clean_mtcars_example()}} \code{\link[=drake_examples]{drake_examples()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{drake_done}
\alias{drake_done}
\title{List done targets.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_done(cache = drake::drake_cache(path = path), path = NULL)
}
\arguments{
\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}
}
\value{
A character vector of target names.
}
\description{
List the targets that completed in the current or
previous call to \code{\link[=make]{make()}}.
}
\examples{
\dontrun{
isolate_example("contain side effects", {
plan <- drake_plan(x = 1, y = x)
make(plan)
drake_done()
})
}
}
\seealso{
\code{\link[=drake_running]{drake_running()}}, \code{\link[=drake_failed]{drake_failed()}}, \code{\link[=drake_cancelled]{drake_cancelled()}},
\code{\link[=drake_progress]{drake_progress()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{read_config}
\alias{read_config}
\title{read_config \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
read_config(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{configure_cache}
\alias{configure_cache}
\title{Configure the hash algorithms, etc. of a drake cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
configure_cache(
  cache = drake::get_cache(verbose = verbose),
  short_hash_algo = drake::default_short_hash_algo(cache = cache),
  long_hash_algo = drake::default_long_hash_algo(cache = cache),
  log_progress = FALSE,
  overwrite_hash_algos = FALSE,
  verbose = 1L,
  jobs = 1,
  init_common_values = FALSE
)
}
\arguments{
\item{cache}{Cache to configure}

\item{short_hash_algo}{Short hash algorithm for drake.
The short algorithm must be among \code{\link[=available_hash_algos]{available_hash_algos()}},
which is just the collection of algorithms available to the \code{algo}
argument in \code{digest::digest()}.
See \code{\link[=default_short_hash_algo]{default_short_hash_algo()}} for more.}

\item{long_hash_algo}{Long hash algorithm for drake.
The long algorithm must be among \code{\link[=available_hash_algos]{available_hash_algos()}},
which is just the collection of algorithms available to the \code{algo}
argument in \code{digest::digest()}.
See \code{\link[=default_long_hash_algo]{default_long_hash_algo()}} for more.}

\item{log_progress}{Deprecated logical.
Previously toggled whether to clear the recorded
build progress if this cache was used for previous calls to
\code{\link[=make]{make()}}.}

\item{overwrite_hash_algos}{Logical, whether to try to overwrite
the hash algorithms in the cache with any user-specified ones.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{jobs}{Number of jobs for parallel processing}

\item{init_common_values}{Logical, whether to set the initial \code{drake}
version in the cache and other common values.
Not always a thread safe operation, so should only be \code{TRUE}
on the main process}
}
\value{
A drake/storr cache.
}
\description{
The purpose of this function is
to prepare the cache to be called from \code{\link[=make]{make()}}.
\code{drake} only uses a single hash algorithm now,
so we no longer need this configuration step.
}
\details{
Deprecated on 2018-12-12.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_plan_helpers.R
\name{bind_plans}
\alias{bind_plans}
\title{Row-bind together drake plans
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
bind_plans(...)
}
\arguments{
\item{...}{Workflow plan data frames (see \code{\link[=drake_plan]{drake_plan()}}).}
}
\description{
Combine drake plans together in a way that
correctly fills in missing entries.
}
\examples{
# You might need to refresh your data regularly (see ?triggers).
download_plan <- drake_plan(
  data = target(
    command = download_data(),
    trigger = "always"
  )
)
# But if the data don't change, the analyses don't need to change.
analysis_plan <- drake_plan(
  usage = get_usage_metrics(data),
  topline = scrape_topline_table(data)
)
your_plan <- bind_plans(download_plan, analysis_plan)
your_plan
}
\seealso{
\code{\link[=drake_plan]{drake_plan()}}, \code{\link[=make]{make()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_drake_settings.R
\name{drake_settings}
\alias{drake_settings}
\title{\code{drake_settings} helper}
\usage{
drake_settings(
  cache_log_file = NULL,
  curl_handles = list(),
  garbage_collection = FALSE,
  jobs = 1L,
  jobs_preprocess = 1L,
  keep_going = TRUE,
  lazy_load = "eager",
  lib_loc = character(0),
  lock_cache = TRUE,
  lock_envir = TRUE,
  log_build_times = TRUE,
  log_progress = TRUE,
  memory_strategy = "speed",
  parallelism = "loop",
  recover = TRUE,
  recoverable = TRUE,
  seed = 0L,
  session_info = TRUE,
  skip_imports = FALSE,
  skip_safety_checks = FALSE,
  skip_targets = FALSE,
  sleep = function(i) 0.01,
  template = list(),
  log_worker = FALSE
)
}
\value{
A \code{drake_settings} object.
}
\description{
List of class \code{drake_settings}.
}
\examples{
if (FALSE) { # stronger than roxygen dontrun
drake_settings()
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{read_drake_graph}
\alias{read_drake_graph}
\title{Read a workflow graph from the cache
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
read_drake_graph(path = getwd(), search = TRUE, cache = NULL, verbose = 1L)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\description{
drake no longer stores the config object,
the plan, etc. in the cache during \code{make()}. This change
improves speed.
}
\details{
2019-01-06
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_ggraph.R
\name{render_drake_ggraph}
\alias{render_drake_ggraph}
\title{Visualize the workflow with \code{ggplot2}/\code{ggraph} using
\code{\link[=drake_graph_info]{drake_graph_info()}} output.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
render_drake_ggraph(
  graph_info,
  main = graph_info$default_title,
  label_nodes = FALSE,
  transparency = TRUE
)
}
\arguments{
\item{graph_info}{List of data frames generated by
\code{\link[=drake_graph_info]{drake_graph_info()}}.
There should be 3 data frames: \code{nodes}, \code{edges},
and \code{legend_nodes}.}

\item{main}{Character string, title of the graph.}

\item{label_nodes}{Logical, whether to label the nodes.
If \code{FALSE}, the graph will not have any text next to the nodes,
which is recommended for large graphs with lots of targets.}

\item{transparency}{Logical, whether to allow transparency in
the rendered graph. Set to \code{FALSE} if you get warnings
like "semi-transparency is not supported on this device".}
}
\value{
A \code{ggplot2} object, which you can modify with more layers,
show with \code{plot()}, or save as a file with \code{ggsave()}.
}
\description{
This function requires packages \code{ggplot2} and \code{ggraph}.
Install them with \code{install.packages(c("ggplot2", "ggraph"))}.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
load_mtcars_example() # Get the code with drake_example("mtcars").
if (requireNamespace("ggraph", quietly = TRUE)) {
  # Instead of jumpting right to vis_drake_graph(), get the data frames
  # of nodes, edges, and legend nodes.
  drake_ggraph(my_plan) # Jump straight to the static graph.
  # Get the node and edge info that vis_drake_graph() just plotted:
  graph <- drake_graph_info(my_plan)
  render_drake_ggraph(graph)
}
})
}
}
\seealso{
\code{\link[=vis_drake_graph]{vis_drake_graph()}}, \code{\link[=sankey_drake_graph]{sankey_drake_graph()}}, \code{\link[=drake_ggraph]{drake_ggraph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_graph_info.R
\name{drake_graph_info}
\alias{drake_graph_info}
\title{Prepare the workflow graph for visualization
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
drake_graph_info(
  ...,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  build_times = "build",
  digits = 3,
  targets_only = FALSE,
  font_size = 20,
  from_scratch = FALSE,
  make_imports = TRUE,
  full_legend = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE,
  hover = FALSE,
  on_select_col = NULL,
  config = NULL
)
}
\arguments{
\item{...}{Arguments to \code{\link[=make]{make()}}, such as \code{plan} and \code{targets}.}

\item{from}{Optional collection of target/import names.
If \code{from} is nonempty,
the graph will restrict itself to
a neighborhood of \code{from}.
Control the neighborhood with
\code{mode} and \code{order}.}

\item{mode}{Which direction to branch out in the graph
to create a neighborhood around \code{from}.
Use \code{"in"} to go upstream,
\code{"out"} to go downstream,
and \code{"all"} to go both ways and disregard
edge direction altogether.}

\item{order}{How far to branch out to create
a neighborhood around \code{from}. Defaults to
as far as possible. If a target is in the neighborhood, then
so are all of its custom \code{\link[=file_out]{file_out()}} files if
\code{show_output_files} is \code{TRUE}.
That means the actual graph order may be slightly greater than
you might expect, but this ensures consistency
between \code{show_output_files = TRUE} and
\code{show_output_files = FALSE}.}

\item{subset}{Optional character vector.
Subset of targets/imports to display in the graph.
Applied after \code{from}, \code{mode}, and \code{order}.
Be advised: edges are only kept for adjacent nodes in \code{subset}.
If you do not select all the intermediate nodes,
edges will drop from the graph.}

\item{build_times}{Character string or logical.
If character, the choices are
1. \code{"build"}: runtime of the command plus the time
it take to store the target or import.
2. \code{"command"}: just the runtime of the command.
3. \code{"none"}: no build times.
If logical, \code{build_times} selects whether to show the
times from `build_times(..., type = "build")`` or use
no build times at all. See \code{\link[=build_times]{build_times()}} for details.}

\item{digits}{Number of digits for rounding the build times}

\item{targets_only}{Logical,
whether to skip the imports and only include the
targets in the workflow plan.}

\item{font_size}{Numeric, font size of the node labels in the graph}

\item{from_scratch}{Logical, whether to assume all the targets
will be made from scratch on the next \code{\link[=make]{make()}}.
Makes all targets outdated, but keeps information about
build progress in previous \code{\link[=make]{make()}}s.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to increase speed and risk using obsolete information.}

\item{full_legend}{Logical. If \code{TRUE}, all the node types
are printed in the legend. If \code{FALSE}, only the
node types used are printed in the legend.}

\item{group}{Optional character scalar, name of the column used to
group nodes into columns. All the columns names of your original \code{drake}
plan are choices. The other choices (such as \code{"status"}) are column names
in the \code{nodes} . To group nodes into clusters in the graph,
you must also supply the \code{clusters} argument.}

\item{clusters}{Optional character vector of values to cluster on.
These values must be elements of the column of the \code{nodes} data frame
that you specify in the \code{group} argument to \code{drake_graph_info()}.}

\item{show_output_files}{Logical, whether to include
\code{\link[=file_out]{file_out()}} files in the graph.}

\item{hover}{Logical, whether to show text (file contents,
commands, etc.) when you hover your cursor over a node.}

\item{on_select_col}{Optional string corresponding to the column name
in the plan that should provide data for the \code{on_select} event.}

\item{config}{Deprecated.}
}
\value{
A list of three data frames: one for nodes,
one for edges, and one for
the legend nodes. The list also contains the
default title of the graph.
}
\description{
With the returned data frames,
you can plot your own custom \code{visNetwork} graph.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (requireNamespace("visNetwork", quietly = TRUE)) {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
vis_drake_graph(my_plan)
# Get a list of data frames representing the nodes, edges,
# and legend nodes of the visNetwork graph from vis_drake_graph().
raw_graph <- drake_graph_info(my_plan)
# Choose a subset of the graph.
smaller_raw_graph <- drake_graph_info(
  my_plan,
  from = c("small", "reg2"),
  mode = "in"
)
# Inspect the raw graph.
str(raw_graph)
# Use the data frames to plot your own custom visNetwork graph.
# For example, you can omit the legend nodes
# and change the direction of the graph.
library(visNetwork)
graph <- visNetwork(nodes = raw_graph$nodes, edges = raw_graph$edges)
visHierarchicalLayout(graph, direction = 'UD')
}
}
})
}
}
\seealso{
\code{\link[=vis_drake_graph]{vis_drake_graph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{default_recipe_command}
\alias{default_recipe_command}
\title{Default Makefile recipe command
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
default_recipe_command()
}
\value{
A character scalar with the default recipe command.
}
\description{
2019-01-02
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{default_short_hash_algo}
\alias{default_short_hash_algo}
\title{Return the default short hash algorithm for \code{make()}.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
default_short_hash_algo(cache = NULL)
}
\arguments{
\item{cache}{Optional drake cache.
When you \code{\link[=configure_cache]{configure_cache()}} without
supplying a short hash algorithm,
\code{default_short_hash_algo(cache)} is the short
hash algorithm that drake picks for you.}
}
\value{
A character vector naming a hash algorithm.
}
\description{
Deprecated. drake now only uses one hash algorithm per cache.
}
\details{
Deprecated on 2018-12-12
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deps.R
\name{tracked}
\alias{tracked}
\title{List the targets and imports that are reproducibly tracked.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
tracked(config)
}
\arguments{
\item{config}{An output list from \code{\link[=drake_config]{drake_config()}}.}
}
\value{
A character vector with the names of reproducibly-tracked targets.
}
\description{
List all the spec
in your project's dependency network.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Load the canonical example for drake.
# List all the targets/imports that are reproducibly tracked.
config <- drake_config(my_plan)
tracked(config)
}
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drake_ggraph.R
\name{drake_ggraph_impl}
\alias{drake_ggraph_impl}
\title{Internal function with a drake_config() argument}
\usage{
drake_ggraph_impl(
  config,
  build_times = "build",
  digits = 3,
  targets_only = FALSE,
  main = NULL,
  from = NULL,
  mode = c("out", "in", "all"),
  order = NULL,
  subset = NULL,
  make_imports = TRUE,
  from_scratch = FALSE,
  full_legend = FALSE,
  group = NULL,
  clusters = NULL,
  show_output_files = TRUE,
  label_nodes = FALSE,
  transparency = TRUE
)
}
\arguments{
\item{config}{A \code{\link[=drake_config]{drake_config()}} object.}

\item{make_imports}{Logical, whether to make the imports first.
Set to \code{FALSE} to save some time and risk obsolete output.}
}
\description{
Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamic.R
\name{subtargets}
\alias{subtargets}
\title{List sub-targets \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
subtargets(
  target = NULL,
  character_only = FALSE,
  cache = drake::drake_cache(path = path),
  path = NULL
)
}
\arguments{
\item{target}{Character string or symbol, depending on \code{character_only}.
Name of a dynamic target.}

\item{character_only}{Logical, whether \code{target} should be treated
as a character or a symbol.
Just like \code{character.only} in \code{\link[=library]{library()}}.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}
}
\value{
Character vector of sub-target names
}
\description{
List the sub-targets of a dynamic target.
}
\examples{
\dontrun{
isolate_example("dynamic branching", {
plan <- drake_plan(
  w = c("a", "a", "b", "b"),
  x = seq_len(4),
  y = target(x + 1, dynamic = map(x)),
  z = target(sum(x) + sum(y), dynamic = group(x, y, .by = w))
)
make(plan)
subtargets(y)
subtargets(z)
readd(x)
readd(y)
readd(z)
})
}
}
\seealso{
\code{\link[=get_trace]{get_trace()}}, \code{\link[=read_trace]{read_trace()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{make_imports}
\alias{make_imports}
\title{Just process the imports
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
make_imports(config)
}
\arguments{
\item{config}{A configuration list returned by \code{\link[=drake_config]{drake_config()}}.}
}
\value{
nothing
}
\description{
Deprecated on 2019-01-04
}
\seealso{
\code{\link[=make]{make()}}, \code{\link[=drake_config]{drake_config()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{check}
\alias{check}
\title{check \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#defunct}{\figure{lifecycle-defunct.svg}{options: alt='[Defunct]'}}}{\strong{[Defunct]}}}
\usage{
check(...)
}
\arguments{
\item{...}{Arguments}
}
\description{
2019-02-15
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{failed}
\alias{failed}
\title{List failed targets.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
failed(
  path = NULL,
  search = NULL,
  cache = drake::drake_cache(path = path),
  verbose = 1L,
  upstream_only = NULL
)
}
\arguments{
\item{path}{Path to a \code{drake} cache
(usually a hidden \verb{.drake/} folder) or \code{NULL}.}

\item{search}{Deprecated.}

\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}

\item{upstream_only}{Deprecated.}
}
\value{
A character vector of target names.
}
\description{
Deprecated on 2020-03-23. Use \code{\link[=drake_failed]{drake_failed()}} instead.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{clean_main_example}
\alias{clean_main_example}
\title{Deprecated: clean the main example from \code{drake_example("main")}
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
clean_main_example()
}
\value{
Nothing.
}
\description{
This function deletes files. Use at your own risk.
Destroys the \verb{.drake/} cache and the \code{report.Rmd} file
in the current working directory. Your working directory
(\code{getcwd()}) must be the folder from which you first ran
\code{load_main_example()} and \code{make(my_plan)}.
}
\details{
Deprecated 2018-12-31.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{long_hash}
\alias{long_hash}
\title{\code{drake} now has just one hash algorithm per cache.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}
\usage{
long_hash(cache = drake::get_cache(verbose = verbose), verbose = 1L)
}
\arguments{
\item{cache}{drake cache. See \code{\link[=new_cache]{new_cache()}}.
If supplied, \code{path} is ignored.}

\item{verbose}{Deprecated on 2019-09-11.}
}
\value{
A character vector naming a hash algorithm.
}
\description{
Deprecated on 2018-12-12
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text_drake_graph.R
\name{render_text_drake_graph}
\alias{render_text_drake_graph}
\title{Show a workflow graph as text in your terminal window
using \code{\link[=drake_graph_info]{drake_graph_info()}} output.
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}}
\usage{
render_text_drake_graph(graph_info, nchar = 1L, print = TRUE)
}
\arguments{
\item{graph_info}{List of data frames generated by
\code{\link[=drake_graph_info]{drake_graph_info()}}.
There should be 3 data frames: \code{nodes}, \code{edges},
and \code{legend_nodes}.}

\item{nchar}{For each node, maximum number of characters of the node label
to show. Can be 0, in which case each node is a colored box
instead of a node label.
Caution: \code{nchar} > 0 will mess with the layout.}

\item{print}{Logical. If \code{TRUE}, the graph will print to the console
via \code{message()}. If \code{FALSE}, nothing is printed. However, you still
have the visualization because \code{text_drake_graph()} and
\code{render_text_drake_graph()} still invisibly return a character string
that you can print yourself with \code{message()}.}
}
\value{
The lines of text in the visualization.
}
\description{
This function is called inside
\code{\link[=text_drake_graph]{text_drake_graph()}}, which typical users
call more often. See \code{?text_drake_graph} for details.
}
\examples{
\dontrun{
isolate_example("Quarantine side effects.", {
if (suppressWarnings(require("knitr"))) {
load_mtcars_example() # Get the code with drake_example("mtcars").
pkgs <- requireNamespace("txtplot", quietly = TRUE) &&
  requireNamespace("visNetwork", quietly = TRUE)
if (pkgs) {
# Instead of jumpting right to vis_drake_graph(), get the data frames
# of nodes, edges, and legend nodes.
text_drake_graph(my_plan) # Jump straight to the interactive graph.
# Get the node and edge info that vis_drake_graph() just plotted:
graph <- drake_graph_info(my_plan)
# You can pass the data frames right to render_text_drake_graph().
render_text_drake_graph(graph)
}
}
})
}
}
\seealso{
\code{\link[=text_drake_graph]{text_drake_graph()}}, \code{\link[=vis_drake_graph]{vis_drake_graph()}},
\code{\link[=sankey_drake_graph]{sankey_drake_graph()}}, \code{\link[=drake_ggraph]{drake_ggraph()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/walk_code.R
\name{walk_code}
\alias{walk_code}
\title{Static code analysis}
\usage{
walk_code(expr, results, locals, restrict)
}
\arguments{
\item{expr}{A function or expression.}

\item{results}{A \code{drake_deps} object.}

\item{locals}{An environment, a hash table of local variables.}

\item{restrict}{An environment,
a hash table for whitelisting global symbols.}
}
\description{
Static code analysis.
}
\keyword{internal}
\newcommand{\lifecycle}{\Sexpr[results=rd, stage=render]{drake:::lifecycle("#1")}}
