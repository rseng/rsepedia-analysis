
# tarchetypes <img src='man/figures/logo.png' align="right" height="139"/>

[![ropensci](https://badges.ropensci.org/401_status.svg)](https://github.com/ropensci/software-review/issues/401)
[![zenodo](https://zenodo.org/badge/282774543.svg)](https://zenodo.org/badge/latestdoi/282774543)
[![R
Targetopia](https://img.shields.io/badge/R_Targetopia-member-blue?style=flat&labelColor=gray)](https://wlandau.github.io/targetopia/)
[![CRAN](https://www.r-pkg.org/badges/version/tarchetypes)](https://CRAN.R-project.org/package=tarchetypes)
[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/tarchetypes/workflows/check/badge.svg)](https://github.com/ropensci/tarchetypes/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/tarchetypes/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/tarchetypes)
[![lint](https://github.com/ropensci/tarchetypes/workflows/lint/badge.svg)](https://github.com/ropensci/tarchetypes/actions?query=workflow%3Alint)

The `tarchetypes` R package is a collection of target and pipeline
archetypes for the [`targets`](https://github.com/ropensci/targets)
package. These archetypes express complicated pipelines with concise
syntax, which enhances readability and thus reproducibility. Archetypes
are possible because of the flexible metaprogramming capabilities of
[`targets`](https://github.com/ropensci/targets). In
[`targets`](https://github.com/ropensci/targets), one can define a
target as an object outside the central pipeline, and the
[`tar_target_raw()`](https://docs.ropensci.org/targets/reference/tar_target_raw.html)
function completely avoids non-standard evaluation. That means anyone
can write their own niche interfaces for specialized projects.
`tarchetypes` aims to include the most common and versatile archetypes
and usage patterns.

## Grouped data frames

`tarchetypes` has functions for easy dynamic branching over subsets of
data frames:

-   `tar_group_by()`: define row groups using `dplyr::group_by()`
    semantics.
-   `tar_group_select()`: define row groups using `tidyselect`
    semantics.
-   `tar_group_count()`: define a given number row groups.
-   `tar_group_size()`: define row groups of a given size.

If you define a target with one of these functions, all downstream
dynamic targets will automatically branch over the row groups.

``` r
# _targets.R file:
library(targets)
library(tarchetypes)
produce_data <- function() {
  expand.grid(var1 = c("a", "b"), var2 = c("c", "d"), rep = c(1, 2, 3))
}
list(
  tar_group_by(data, produce_data(), var1, var2),
  tar_target(group, data, pattern = map(data))
)
```

``` r
# R console:
library(targets)
tar_make()
#> • start target data
#> • built target data
#> • start branch group_b3d7d010
#> • built branch group_b3d7d010
#> • start branch group_6a76c5c0
#> • built branch group_6a76c5c0
#> • start branch group_164b16bf
#> • built branch group_164b16bf
#> • start branch group_f5aae602
#> • built branch group_f5aae602
#> • built pattern group
#> • end pipeline

# First row group:
tar_read(group, branches = 1)
#> # A tibble: 3 × 4
#>   var1  var2    rep tar_group
#>   <fct> <fct> <dbl>     <int>
#> 1 a     c         1         1
#> 2 a     c         2         1
#> 3 a     c         3         1

# Second row group:
tar_read(group, branches = 2)
#> # A tibble: 3 × 4
#>   var1  var2    rep tar_group
#>   <fct> <fct> <dbl>     <int>
#> 1 a     d         1         2
#> 2 a     d         2         2
#> 3 a     d         3         2
```

## Literate programming

Consider the following R Markdown report.

    ---
    title: report
    output: html_document
    ---

    ```{r}
    library(targets)
    tar_read(dataset)
    ```

We want to define a target to render the report. And because the report
calls `tar_read(dataset)`, this target needs to depend on `dataset`.
Without `tarchetypes`, it is cumbersome to set up the pipeline
correctly.

``` r
# _targets.R
library(targets)
list(
  tar_target(dataset, data.frame(x = letters)),
  tar_target(
    report, {
      # Explicitly mention the symbol `dataset`.
      list(dataset)
      # Return relative paths to keep the project portable.
      fs::path_rel(
        # Need to return/track all input/output files.
        c( 
          rmarkdown::render(
            input = "report.Rmd",
            # Always run from the project root
            # so the report can find _targets/.
            knit_root_dir = getwd(),
            quiet = TRUE
          ),
          "report.Rmd"
        )
      )
    },
    # Track the input and output files.
    format = "file",
    # Avoid building small reports on HPC.
    deployment = "main"
  )
)
```

With `tarchetypes`, we can simplify the pipeline with the `tar_render()`
archetype.

``` r
# _targets.R
library(targets)
library(tarchetypes)
list(
  tar_target(dataset, data.frame(x = letters)),
  tar_render(report, "report.Rmd")
)
```

Above, `tar_render()` scans code chunks for mentions of targets in
`tar_load()` and `tar_read()`, and it enforces the dependency
relationships it finds. In our case, it reads `report.Rmd` and then
forces `report` to depend on `dataset`. That way, `tar_make()` always
processes `dataset` before `report`, and it automatically reruns
`report.Rmd` whenever `dataset` changes.

## Alternative pipeline syntax

[`tar_plan()`](https://docs.ropensci.org/tarchetypes/reference/tar_plan.html)
is a drop-in replacement for
[`drake_plan()`](https://docs.ropensci.org/drake/reference/drake_plan.html)
in the [`targets`](https://github.com/ropensci/targets) ecosystem. It
lets users write targets as name/command pairs without having to call
[`tar_target()`](https://docs.ropensci.org/targets/reference/tar_target.html).

``` r
tar_plan(
  tar_file(raw_data_file, "data/raw_data.csv", format = "file"),
  # Simple drake-like syntax:
  raw_data = read_csv(raw_data_file, col_types = cols()),
  data =raw_data %>%
    mutate(Ozone = replace_na(Ozone, mean(Ozone, na.rm = TRUE))),
  hist = create_plot(data),
  fit = biglm(Ozone ~ Wind + Temp, data),
  # Needs tar_render() because it is a target archetype:
  tar_render(report, "report.Rmd")
)
```

## Installation

| Type        | Source   | Command                                                               |
|-------------|----------|-----------------------------------------------------------------------|
| Release     | CRAN     | `install.packages("tarchetypes")`                                     |
| Development | GitHub   | `remotes::install_github("ropensci/tarchetypes")`                     |
| Development | rOpenSci | `install.packages("tarchetypes", repos = "https://dev.ropensci.org")` |

## Documentation

For specific documentation on `tarchetypes`, including the help files of
all user-side functions, please visit the [reference
website](https://docs.ropensci.org/tarchetypes/). For documentation on
[`targets`](https://github.com/ropensci/targets) in general, please
visit the [`targets` reference
website](https://docs.ropensci.org/targets/). Many of the linked
resources use `tarchetypes` functions such as
[`tar_render()`](https://docs.ropensci.org/tarchetypes/reference/tar_render.html).

## Code of conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/).

## Citation

``` r
citation("tarchetypes")
#> 
#> To cite tarchetypes in publications use:
#> 
#>   William Michael Landau (2021). tarchetypes: Archetypes for Targets.
#>   https://docs.ropensci.org/tarchetypes/,
#>   https://github.com/ropensci/tarchetypes.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {tarchetypes: Archetypes for Targets},
#>     author = {William Michael Landau},
#>     year = {2021},
#>     note = {{https://docs.ropensci.org/tarchetypes/, https://github.com/ropensci/tarchetypes}},
#>   }
```
# tarchetypes 0.4.1.9000

* Implement `tar_file_read()` (#84, @petrbouchal).

# tarchetypes 0.4.1

* Select list elements from `command1` using `[[` and not `[` in `tar_map2()` functions.

# tarchetypes 0.4.0

* Implement `tar_map_rep()` and `tar_map_rep_raw()` for dynamic batched replication within static branching for data frames (#78).
* Implement `tar_map2_count()`, `tar_map2_count_raw()`, `tar_map2_size()`, and `tar_map2_size_raw()` for batched dynamic-within-static branching for data frames (#78).
* Deprecate `tar_rep_map()` in favor of `tar_rep2()` to avoid name confusion. Likewise with `tar_rep_map_raw()` to `tar_rep2_raw()` (#78).

# tarchetypes 0.3.2

* Allow empty / `NULL` target list in `tar_map()` (@kkami1115).
* Do not claim to support `"aws_file"` format in `tar_files()` or related target factories.

# tarchetypes 0.3.1

* Relax assertion on language objects.
* Explain `targets` timestamps correctly in the help files of `tar_age()` and `tar_cue_age()`.

# tarchetypes 0.3.0

## Invalidating changes

* When `names = NULL` in `tar_map()`, use hashes instead of numeric indexes for generated target names (#67). That way, target names are no longer sensitive to the order of `values`, and so targets will incorrectly invalidate less often. *Unfortunately, this is an invalidating change: some targets will automatically rerun after you install this version of `tarchetypes`.* I apologize for the inconvenience this causes. However, we do need this patch in order to solve #67, and targets will incorrectly invalidate less frequently in the future.

## Enhancements

* Migrate to utilities for error handling and metaprogramming exported from `targets` (#59).

# tarchetypes 0.2.1

## Bug fixes

* Make the `*_raw()` target factories process `command` the same way whether it is an expression or ordinary language object.
* Ensure compatibility with `targets` 0.5.0.9000, which logs skipped targets.

## New features

* Add `tar_rep_map()` and `tar_rep_map_raw()` to perform batched computation downstream of `tar_rep()` (#50).
* Add `tar_select_names()` and `tar_select_targets()` to make certain metaprogramming tasks easier.
* In `tar_map()`, attempt to convert the elements of `values` into lists of language objects.

# tarchetypes 0.2.0

* Allow trailing commas in `tar_plan()` (#40, @kendonB).
* Implement `tar_age()` based on `tar_cue_age()` (#39, @petrbouchal).
* Implement new cue factories `tar_cue_age()`, `tar_cue_age_raw()`, `tar_cue_force()`, and `tar_cue_skip()` (#39).
* Implement `tar_download()` (#38, @noamross, @petrbouchal)
* Set intermediate temporary directory to remove race condition in `tar_render_rep()` (#36, @gorgitko). 
* Prefix internal condition classes with "tar_".
* Add new format helpers such as `tar_aws_rds()` and `tar_parquet()`.
* Support hooks `tar_hook_before()`, `tar_hook_inner()`, and `tar_hook_outer()` (#44).
* Deep-copy the cue in `tar_map()`.

# tarchetypes 0.1.1

* Unset `crayon.enabled` for literate programming.
* Switch meaning of `%||%` and `%|||%` to conform to historical precedent.

# tarchetypes 0.1.0

* Add new functions for easier grouping of data frames for dynamic branching: `tar_group_by()`, `tar_group_select()`, `tar_group_size()`, `tar_group_count()` (#32, @liutiming).
* In `tar_render()` and related functions, track the `*_files/` output directory if it exists (#30).
* Implement an external `walk_ast()` function to make it easier for other developers to extend the static code analysis of `tarchetypes` (@MilesMcBain).

# tarchetypes 0.0.4

* Skip literate programming tests if pandoc is missing or has an insufficient version.
* Use explicit temp files in examples even when running inside `targets::tar_dir()`. (`targets::tar_dir()` and `targets::tar_test()` already run code in a temporary directory.)
* Add comments in the examples to emphasize that `targets::tar_dir()` runs code in a temporary directory, which means all ostensibly files created in the enclosed expression will actually be written to temporary storage and not the user's file space.

# tarchetypes 0.0.2

* Make sure every function with a help file in `man/` has Rd-tags `\value` and `\arguments`.
* For every function with a help file in `man/`, describe the return value in the `\value` Rd tag. For each function that returns a target object or list of target objects, the `\value` tag now links to <https://books.ropensci.org/targets/>, the user manual where the purpose of target objects is explained, and <https://books.ropensci.org/targets-design/>, the design specification which documents the structure and composition of target objects.
* Ensure that examples, vignettes, and test do not write to the home file space of the user.
* Ensure that no function defined in the `tarchetypes` package writes by default to the home file space of the user. The paths of all output files are controlled by non-`tarchetypes` functions that invoke `tarchetypes`.

# tarchetypes 0.0.1

* `tar_plan()` now returns a list of target objects rather than a pipeline object. Related: <https://github.com/ropensci/targets/issues/253>.

# tarchetypes 0.0.0.9000

* First version.
* Implement `tar_knitr_deps()` and `tar_knitr_deps_expr()` to accommodate custom multi-file literate programming projects like R Markdown sites and `bookdown` projects (#23, @tjmahr).
# Contributing

Development is a community effort, and we welcome participation.

## Code of conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/).

## Discussions

At <https://github.com/ropensci/tarchetypes/discussions>, you can post general questions, brainstorm ideas, and ask for help.

## Issues

<https://github.com/ropensci/tarchetypes/issues> is for bug reports, performance issues, package maintenance tasks, and feature requests. When you post, please abide by the following guidelines.

* Before posting a new issue, please take a moment to search for existing similar issues in order to avoid duplication.
* For bug reports: if you can, please install the latest GitHub version of `tarchetypes` (i.e. `remotes::install_github("ropensci/tarchetypes")`) and verify that the issue still persists.
* Describe your issue in prose as clearly and concisely as possible.
* For any problem you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot. A reproducible example is:
    * **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Development

External code contributions are extremely helpful in the right circumstances. Here are the recommended steps.

1. Prior to contribution, please propose your idea in a [new issue thread](https://github.com/ropensci/tarchetypes/issues) so you and the maintainer can define the intent and scope of your work.
2. [Fork the repository](https://help.github.com/articles/fork-a-repo/).
3. Follow the [GitHub flow](https://guides.github.com/introduction/flow/index.html) to create a new branch, add commits, and open a pull request.
4. Discuss your code with the maintainer in the pull request thread.
5. If everything looks good, the maintainer will merge your code into the project.

Please also follow these additional guidelines.

* Respect the architecture and reasoning of the package. Depending on the scope of your work, you may want to read the design documents (package vignettes).
* If possible, keep contributions small enough to easily review manually. It is okay to split up your work into multiple pull requests.
* Format your code according to the [tidyverse style guide](https://style.tidyverse.org/) and check your formatting with the `lint_package()` function from the [`lintr`](https://github.com/jimhester/lintr) package.
* For new features or functionality, add tests in `tests`. Tests that can be automated should go in `tests/testthat/`. Tests that cannot be automated should go in `tests/interactive/`. For features affecting performance, it is good practice to add profiling studies to `tests/performance/`.
* Check code coverage with `covr::package_coverage()`. Automated tests should cover all the new or changed functionality in your pull request.
* Run overall package checks with `devtools::check()` and `goodpractice::gp()`
* Describe your contribution in the project's [`NEWS.md`](https://github.com/ropensci/tarchetypes/blob/main/NEWS.md) file. Be sure to mention relevent GitHub issue numbers and your GitHub name as done in existing news entries.
* If you feel contribution is substantial enough for official author or contributor status, please add yourself to the `Authors@R` field of the [`DESCRIPTION`](https://github.com/ropensci/tarchetypes/blob/main/DESCRIPTION) file.
# Prework

* [ ] I understand and agree to this repository's [code of conduct](https://ropensci.org/code-of-conduct/).
* [ ] I understand and agree to this repository's [contributing guidelines](https://github.com/ropensci/tarchetypes/blob/main/CONTRIBUTING.md).
* [ ] I have already submitted an [issue](http://github.com/ropensci/tarchetypes/issues) or [discussion thread](https://github.com/ropensci/tarchetypes/discussions) to discuss my idea with the maintainer.

# Related GitHub issues and pull requests

* Ref: #

# Summary

Please explain the purpose and scope of your contribution.
---
name: Maintenance
about: "Something in tarchetypes needs work: updates, documentation, etc. Not a bug, performance issue, or new feature."
title: ""
labels: "type: new maintenance"
assignees: ""
---

## Prework

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/tarchetypes/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/tarchetypes/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] For any problems you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/tarchetypes/issues/256#issuecomment-754229683) so the maintainer can troubleshoot. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the issue.

To help us read any code you include (optional) please try to follow the [tidyverse style guide](https://style.tidyverse.org/). The `style_text()` and `style_file()` functions from the [`styler`](https://github.com/r-lib/styler) package make it easier.

## Reproducible example

* [ ] For any problems you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).
---
name: Bug
about: Please do not submit a bug report unless your issue is a genuine bug in tarchetypes and not a known limitation, usage error, or issue from another package that tarchetypes depends on.
title: ""
labels: "type: bug"
assignees: wlandau
---

## Prework

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/tarchetypes/blob/main/CONTRIBUTING.md).
* [ ] Confirm that your issue is a genuine bug in `tarchetypes` and not a known limitation, usage error, or issue from another package that `tarchetypes` depends on. If you are unsure, please submit a discussion thread instead.
* [ ] If there is [already a relevant issue](https://github.com/ropensci/tarchetypes/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the bug.

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Expected result

What should have happened? Please be as specific as possible.

## Diagnostic information

* A [reproducible example](https://github.com/tidyverse/reprex).
* Session info, available through `sessionInfo()` or [`reprex(si = TRUE)`](https://github.com/tidyverse/reprex).
* A stack trace from `traceback()` or `rlang::trace_back()`.
* The [SHA-1 hash](https://git-scm.com/book/en/v1/Getting-Started-Git-Basics#Git-Has-Integrity) of the GitHub commit of `tarchetypes` currently installed. `packageDescription("tarchetypes")$GithubSHA1` shows you this.
---
name: New feature
about: Suggest a new feature.
title: ""
labels: "type: new feature"
assignees: wlandau
---

## Prework

* [ ] I understand and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and  [contributing guidelines](https://github.com/ropensci/tarchetypes/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/tarchetypes/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] New features take time and effort to create, and they take even more effort to maintain. So if the purpose of the feature is to resolve a struggle you are encountering personally, please consider first posting a "trouble" or "other" issue so we can discuss your use case and search for existing solutions first.

## Proposal

Please describe the new feature. If applicable, write a minimal example in R code or pseudo-code to show input, usage, and desired output.

To help us read any code you include (optional) please try to follow the [tidyverse style guide](https://style.tidyverse.org/). The `style_text()` and `style_file()` functions from the [`styler`](https://github.com/r-lib/styler) package make it easier.
---
name: Performance
about: "Runtime, memory, or storage inefficiency"
title: ""
labels: "topic: performance"
assignees: wlandau

---

## Prework

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/tarchetypes/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/tarchetypes/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the performance issue.

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Benchmarks

How poorly does `tarchetypes` perform? To find out, we recommend you use the [`proffer`](https://github.com/r-prof/proffer) package and take screenshots of the results displayed in your browser.

```r
library(tarchetypes)
library(proffer)
px <- pprof({
  # All your tarchetypes code goes here.
})
```
---
output: github_document
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# tarchetypes <img src='man/figures/logo.png' align="right" height="139"/>

[![ropensci](https://badges.ropensci.org/401_status.svg)](https://github.com/ropensci/software-review/issues/401)
[![zenodo](https://zenodo.org/badge/282774543.svg)](https://zenodo.org/badge/latestdoi/282774543)
[![R Targetopia](https://img.shields.io/badge/R_Targetopia-member-blue?style=flat&labelColor=gray)](https://wlandau.github.io/targetopia/)
[![CRAN](https://www.r-pkg.org/badges/version/tarchetypes)](https://CRAN.R-project.org/package=tarchetypes)
[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/tarchetypes/workflows/check/badge.svg)](https://github.com/ropensci/tarchetypes/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/tarchetypes/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/tarchetypes)
[![lint](https://github.com/ropensci/tarchetypes/workflows/lint/badge.svg)](https://github.com/ropensci/tarchetypes/actions?query=workflow%3Alint)

The `tarchetypes` R package is a collection of target and pipeline archetypes for the [`targets`](https://github.com/ropensci/targets) package. These archetypes express complicated pipelines with concise syntax, which enhances readability and thus reproducibility. Archetypes are possible because of the flexible metaprogramming capabilities of [`targets`](https://github.com/ropensci/targets). In [`targets`](https://github.com/ropensci/targets), one can define a target as an object outside the central pipeline, and the [`tar_target_raw()`](https://docs.ropensci.org/targets/reference/tar_target_raw.html) function completely avoids non-standard evaluation. That means anyone can write their own niche interfaces for specialized projects. `tarchetypes` aims to include the most common and versatile archetypes and usage patterns.

## Grouped data frames

`tarchetypes` has functions for easy dynamic branching over subsets of data frames:

* `tar_group_by()`: define row groups using `dplyr::group_by()` semantics.
* `tar_group_select()`: define row groups using `tidyselect` semantics.
* `tar_group_count()`: define a given number row groups.
* `tar_group_size()`: define row groups of a given size.

If you define a target with one of these functions, all downstream dynamic targets will automatically branch over the row groups.

```{r, echo = FALSE}
targets::tar_script({
  produce_data <- function() {
    expand.grid(var1 = c("a", "b"), var2 = c("c", "d"), rep = c(1, 2, 3))
  }
  list(
    tarchetypes::tar_group_by(data, produce_data(), var1, var2),
    tar_target(group, data, pattern = map(data))
  )
})
```

```{r, eval = FALSE}
# _targets.R file:
library(targets)
library(tarchetypes)
produce_data <- function() {
  expand.grid(var1 = c("a", "b"), var2 = c("c", "d"), rep = c(1, 2, 3))
}
list(
  tar_group_by(data, produce_data(), var1, var2),
  tar_target(group, data, pattern = map(data))
)
```

```{r}
# R console:
library(targets)
tar_make()

# First row group:
tar_read(group, branches = 1)

# Second row group:
tar_read(group, branches = 2)
```

## Literate programming

Consider the following R Markdown report.

```{r, echo = FALSE, comment = ""}
lines <- c(
  "---",
  "title: report",
  "output: html_document",
  "---",
  "",
  "```{r}",
  "library(targets)",
  "tar_read(dataset)",
  "```"
)
cat(lines, sep = "\n")
```

We want to define a target to render the report. And because the report calls `tar_read(dataset)`, this target needs to depend on `dataset`. Without `tarchetypes`, it is cumbersome to set up the pipeline correctly.

```{r, eval = FALSE}
# _targets.R
library(targets)
list(
  tar_target(dataset, data.frame(x = letters)),
  tar_target(
    report, {
      # Explicitly mention the symbol `dataset`.
      list(dataset)
      # Return relative paths to keep the project portable.
      fs::path_rel(
        # Need to return/track all input/output files.
        c( 
          rmarkdown::render(
            input = "report.Rmd",
            # Always run from the project root
            # so the report can find _targets/.
            knit_root_dir = getwd(),
            quiet = TRUE
          ),
          "report.Rmd"
        )
      )
    },
    # Track the input and output files.
    format = "file",
    # Avoid building small reports on HPC.
    deployment = "main"
  )
)
```

With `tarchetypes`, we can simplify the pipeline with the `tar_render()` archetype.

```{r, eval = FALSE}
# _targets.R
library(targets)
library(tarchetypes)
list(
  tar_target(dataset, data.frame(x = letters)),
  tar_render(report, "report.Rmd")
)
```

Above, `tar_render()` scans code chunks for mentions of targets in `tar_load()` and `tar_read()`, and it enforces the dependency relationships it finds. In our case, it reads `report.Rmd` and then forces `report` to depend on `dataset`. That way, `tar_make()` always processes `dataset` before `report`, and it automatically reruns `report.Rmd` whenever `dataset` changes.

## Alternative pipeline syntax

[`tar_plan()`](https://docs.ropensci.org/tarchetypes/reference/tar_plan.html) is a drop-in replacement for [`drake_plan()`](https://docs.ropensci.org/drake/reference/drake_plan.html) in the [`targets`](https://github.com/ropensci/targets) ecosystem. 
It lets users write targets as name/command pairs without having to call [`tar_target()`](https://docs.ropensci.org/targets/reference/tar_target.html).

```{r, eval = FALSE}
tar_plan(
  tar_file(raw_data_file, "data/raw_data.csv", format = "file"),
  # Simple drake-like syntax:
  raw_data = read_csv(raw_data_file, col_types = cols()),
  data =raw_data %>%
    mutate(Ozone = replace_na(Ozone, mean(Ozone, na.rm = TRUE))),
  hist = create_plot(data),
  fit = biglm(Ozone ~ Wind + Temp, data),
  # Needs tar_render() because it is a target archetype:
  tar_render(report, "report.Rmd")
)
```

## Installation

Type | Source | Command
---|---|---
Release | CRAN | `install.packages("tarchetypes")`
Development | GitHub | `remotes::install_github("ropensci/tarchetypes")`
Development | rOpenSci | `install.packages("tarchetypes", repos = "https://dev.ropensci.org")`

## Documentation

For specific documentation on `tarchetypes`, including the help files of all user-side functions, please visit the [reference website](https://docs.ropensci.org/tarchetypes/). For documentation on [`targets`](https://github.com/ropensci/targets) in general, please visit the [`targets` reference website](https://docs.ropensci.org/targets/). Many of the linked resources use `tarchetypes` functions such as [`tar_render()`](https://docs.ropensci.org/tarchetypes/reference/tar_render.html).

## Code of conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/).

## Citation

```{r}
citation("tarchetypes")
```

```{r, echo = FALSE}
unlink("_targets.R")
tar_destroy()
```
---
title: "Example knir report in a targets pipeline"
output: html_document
---

```{r setup, include=FALSE}
library(targets)
knitr::opts_chunk$set(echo = TRUE)
```

This R Markdown report depends on targets from a [`targets`](https://github.com/ropensci/targets) pipeline. Declare these dependencies with `tar_load()` and `tar_read()`.

```{r}
tar_load(data)
```

```{r}
plot(tar_read(analysis))
```

You can include this report itself as a target. Simply use the `tar_render()` archetype to define the target, which tells `targets` to automatically find dependnecies like `data` and `analysis` and rerun the report when either dependency changes.

The code analysis tries to be intelligent, so `tar_load()` and `tar_read()` can appear in various places.

```{r}
f <- function() {
  targets::tar_read(name = data2)
}
```

However, the dependencies are always detected using static code analysis, which means `tidyselect` syntax like `tar_load(starts_with("data"))` will not work. `targets` requires you to write the literal symbol name of every target you want the report to depend on.

```{r}
tar_load(names = c("string1", "string2"))
```

```{r}
tar_read(name = "string3")
```

```{r}
not
a
target
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_force.R
\name{tar_force_change}
\alias{tar_force_change}
\title{Convert a condition into a change.}
\usage{
tar_force_change(condition)
}
\arguments{
\item{condition}{Logical, whether to run the downstream target
in \code{\link[=tar_force]{tar_force()}}.}
}
\value{
A hash that changes when the downstream target is supposed to run.
}
\description{
Supports \code{\link[=tar_force]{tar_force()}}. This is really an internal function
and not meant to be called by users directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_cue_age_raw.R
\name{tar_cue_age_raw}
\alias{tar_cue_age_raw}
\title{Cue to run a target when the last run reaches a certain age
(raw version)}
\usage{
tar_cue_age_raw(
  name,
  age,
  command = TRUE,
  depend = TRUE,
  format = TRUE,
  iteration = TRUE,
  file = TRUE
)
}
\arguments{
\item{name}{Character of length 1, name of the target.}

\item{age}{A \code{difftime} object of length 1, such as
\code{as.difftime(3, units = "days")}. If the target's output data
files are older than \code{age} (according to the most recent
time stamp over all the target's output files)
then the target will rerun.
On the other hand, if at least one data file is
younger than \code{Sys.time() - age}, then the ordinary
invalidation rules apply, and the target may or not rerun.
If you want to force the target to run every 3 days,
for example, set \code{age = as.difftime(3, units = "days")}.}

\item{command}{Logical, whether to rerun the target if command changed
since last time.}

\item{depend}{Logical, whether to rerun the target if the value of one
of the dependencies changed.}

\item{format}{Logical, whether to rerun the target if the user-specified
storage format changed. The storage format is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{iteration}{Logical, whether to rerun the target if the user-specified
iteration method changed. The iteration method is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{file}{Logical, whether to rerun the target if the file(s) with the
return value changed or at least one is missing.}
}
\value{
A cue object. See the "Cue objects" section for background.
}
\description{
\code{tar_cue_age_raw()} acts like \code{tar_cue_age()}
except the \code{name} argument is a character string,
not a symbol. \code{tar_cue_age_raw()} creates a cue object to
rerun a target if the most recent output data becomes old enough.
The age of the target is determined by \code{targets::tar_timestamp()},
and the way the time stamp is calculated is explained
in the Details section of the help file of that function.
}
\details{
\code{tar_cue_age_raw()} uses the time stamps from \code{tar_meta()$time}.
If no time stamp is recorded, the cue defaults to the ordinary
invalidation rules (i.e. \code{mode = "thorough"} in \code{targets::tar_cue()}).
}
\section{Dynamic branches at regular time intervals}{

Time stamps are not recorded for whole dynamic targets,
so \code{tar_age()} is not a good fit for dynamic branching.
To invalidate dynamic branches at regular intervals,
it is recommended to use \code{targets::tar_older()} in combination
with \code{targets::tar_invalidate()} right before calling \code{tar_make()}.
For example,
\code{tar_invalidate(all_of(tar_older(Sys.time - as.difftime(1, units = "weeks"))))} # nolint
invalidates all targets more than a week old. Then, the next \code{tar_make()}
will rerun those targets.
}

\section{Cue objects}{

A cue object is an object generated by \code{targets::tar_cue()},
\code{tarchetypes::tar_cue_force()}, or similar. It is
a collection of decision rules that decide when a target
is invalidated/outdated (e.g. when \code{tar_make()} or similar
reruns the target). You can supply these cue objects to the
\code{tar_target()} function or similar. For example,
\code{tar_target(x, run_stuff(), cue = tar_cue(mode = "always"))}
is a target that always calls \code{run_stuff()} during \code{tar_make()}
and always shows as invalidated/outdated in \code{tar_outdated()},
\code{tar_visnetwork()}, and similar functions.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  library(tarchetypes)
  list(
    targets::tar_target(
      data,
      data.frame(x = seq_len(26)),
      cue = tarchetypes::tar_cue_age_raw(
        name = "data",
        age = as.difftime(0.5, units = "secs")
      )
    )
  )
})
targets::tar_make()
Sys.sleep(0.6)
targets::tar_make()
})
}
}
\seealso{
Other cues: 
\code{\link{tar_age}()},
\code{\link{tar_cue_age}()},
\code{\link{tar_cue_force}()},
\code{\link{tar_cue_skip}()}
}
\concept{cues}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_count.R
\name{tar_group_count}
\alias{tar_group_count}
\title{Group the rows of a data frame into a given number groups}
\usage{
tar_group_count(
  name,
  command,
  count,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{count}{Positive integer, maximum number of row groups}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A target object to generate a grouped data frame
to allows downstream dynamic targets to branch over the
groups of rows.
See the "Target objects" section for background.
}
\description{
Create a target that outputs a grouped data frame
for downstream dynamic branching. Set the maximum
number of groups using \code{count}. The number of rows per group
varies but is approximately uniform.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  produce_data <- function() {
    expand.grid(var1 = c("a", "b"), var2 = c("c", "d"), rep = c(1, 2, 3))
  }
  list(
    tarchetypes::tar_group_count(data, produce_data(), count = 2),
    tar_target(group, data, pattern = map(data))
  )
})
targets::tar_make()
# Read the first row group:
targets::tar_read(group, branches = 1)
# Read the second row group:
targets::tar_read(group, branches = 2)
})
}
}
\seealso{
Other Grouped data frame targets: 
\code{\link{tar_group_by}()},
\code{\link{tar_group_select}()},
\code{\link{tar_group_size}()}
}
\concept{Grouped data frame targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_eval.R
\name{tar_eval}
\alias{tar_eval}
\title{Evaluate multiple expressions created with symbol substitution.}
\usage{
tar_eval(expr, values, envir = parent.frame())
}
\arguments{
\item{expr}{Starting expression. Values are iteratively substituted
in place of symbols in \code{expr} to create each new expression,
and then each new expression is evaluated.}

\item{values}{List of values to substitute into \code{expr} to create
the expressions. All elements of \code{values} must have the same length.}

\item{envir}{Environment in which to evaluate the new expressions.}
}
\value{
A list of return values from the generated expression objects.
Often, these values are target objects.
See the "Target objects" section for background
on target objects specifically.
}
\description{
Loop over a grid of values, create an expression object
from each one, and then evaluate that expression.
Helps with general metaprogramming.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
# tar_map() is incompatible with tar_render() because the latter
# operates on preexisting tar_target() objects. By contrast,
# tar_eval() and tar_sub() iterate over the literal code
# farther upstream.
values <- list(
  name = lapply(c("name1", "name2"), as.symbol),
  file = list("file1.Rmd", "file2.Rmd")
)
tar_sub(list(name, file), values = values)
tar_sub(tar_render(name, file), values = values)
path <- tempfile()
file.create(path)
str(tar_eval(tar_render(name, path), values = values))
# So in your _targets.R file, you can define a pipeline like as below.
# Just make sure to set a unique name for each target
# (which tar_map() does automatically).
values <- list(
  name = lapply(c("name1", "name2"), as.symbol),
  file = c(path, path)
)
list(
  tar_eval(tar_render(name, file), values = values)
)
}
\seealso{
Other Metaprogramming utilities: 
\code{\link{tar_eval_raw}()},
\code{\link{tar_sub_raw}()},
\code{\link{tar_sub}()}
}
\concept{Metaprogramming utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_files_input.R
\name{tar_files_input}
\alias{tar_files_input}
\title{Dynamic branching over input files or URLs}
\usage{
tar_files_input(
  name,
  files,
  batches = length(files),
  format = c("file", "url", "aws_file"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{files}{Nonempty character vector of known existing input files
to track for changes.}

\item{batches}{Positive integer of length 1, number of batches
to partition the files. The default is one file per batch
(maximum number of batches) which is simplest to handle but
could cause a lot of overhead and consume a lot of computing resources.
Consider reducing the number of batches below the number of files
for heavy workloads.}

\item{format}{Character, either \code{"file"} or \code{"url"}. See the \code{format}
argument of \code{targets::tar_target()} for details.}

\item{iteration}{Character, iteration method. Must be a method
supported by the \code{iteration} argument of \code{targets::tar_target()}.
The iteration method for the upstream target is always \code{"list"}
in order to support batching.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{cue}{An optional object from \code{tar_cue()}
to customize the rules that decide whether the target is up to date.
Only applies to the downstream target. The upstream target always runs.}
}
\value{
A list of two targets, one upstream and one downstream.
The upstream one does some work and returns some file paths,
and the downstream target is a pattern that applies \code{format = "file"}
or \code{format = "url"}.
See the "Target objects" section for background.
}
\description{
Dynamic branching over input files or URLs.
}
\details{
\code{tar_files_input()} is like \code{tar_files()}
but more convenient when the files in question already
exist and are known in advance. Whereas \code{tar_files()}
always appears outdated (e.g. with \code{tar_outdated()})
because it always needs to check which files it needs to
branch over, \code{tar_files_input()} will appear up to date
if the files have not changed since last \code{tar_make()}.
In addition, \code{tar_files_input()} automatically groups
input files into batches to reduce overhead and
increase the efficiency of parallel processing.

\code{tar_files_input()} creates a pair of targets, one upstream
and one downstream. The upstream target does some work
and returns some file paths, and the downstream
target is a pattern that applies \code{format = "file"}
or \code{format = "url"}.
This is the correct way to dynamically
iterate over file/url targets. It makes sure any downstream patterns
only rerun some of their branches if the files/urls change.
For more information, visit
\url{https://github.com/ropensci/targets/issues/136} and
\url{https://github.com/ropensci/drake/issues/1302}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  # Do not use temp files in real projects
  # or else your targets will always rerun.
  paths <- unlist(replicate(4, tempfile()))
  file.create(paths)
  list(
    tarchetypes::tar_files_input(
      x,
      paths,
      batches = 2
    )
  )
})
targets::tar_make()
targets::tar_read(x)
targets::tar_read(x, branches = 1)
})
}
}
\seealso{
Other Dynamic branching over files: 
\code{\link{tar_files_input_raw}()},
\code{\link{tar_files_raw}()},
\code{\link{tar_files}()}
}
\concept{Dynamic branching over files}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_rep_map_raw.R
\name{tar_rep_map_raw}
\alias{tar_rep_map_raw}
\title{Dynamic batched computation downstream of \code{\link[=tar_rep]{tar_rep()}}
(raw; deprecated).}
\usage{
tar_rep_map_raw(
  name,
  command,
  targets,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{targets}{Character vector of names of upstream batched targets
created by \code{\link[=tar_rep]{tar_rep()}}.
If you supply more than one such target, all those targets must have the
same number of batches and reps per batch. And they must all return
either data frames or lists. List targets must use \code{iteration = "list"}
in \code{\link[=tar_rep]{tar_rep()}}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A new target object to perform batched computation
downstream of \code{\link[=tar_rep]{tar_rep()}}.
See the "Target objects" section for background.
}
\description{
Deprecated. Use \code{\link[=tar_rep2_raw]{tar_rep2_raw()}} instead.
}
\details{
Deprecated in version 0.4.0, 2021-12-06.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_rep(
      data1,
      data.frame(value = rnorm(1)),
      batches = 2,
      reps = 3
    ),
    tarchetypes::tar_rep(
      data2,
      list(value = rnorm(1)),
      batches = 2, reps = 3,
      iteration = "list" # List iteration is important for batched lists.
    ),
    tarchetypes::tar_rep2_raw( # Use instead of tar_rep_map_raw().
      "aggregate",
      quote(data.frame(value = data1$value + data2$value)),
      targets = c("data1", "data2")
    )
  )
})
targets::tar_make()
targets::tar_read(aggregate)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_combine.R
\name{tar_combine}
\alias{tar_combine}
\title{Static aggregation.}
\usage{
tar_combine(
  name,
  ...,
  command = vctrs::vec_c(!!!.x),
  use_names = TRUE,
  pattern = NULL,
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the new target.}

\item{...}{One or more target objects or list of target objects.
Lists can be arbitrarily nested, as in \code{list()}.}

\item{command}{R command to aggregate the targets. Must contain
\code{!!!.x} where the arguments are to be inserted,
where \verb{!!!} is the unquote splice operator from \code{rlang}.}

\item{use_names}{Logical, whether to insert the names of the targets
into the command when splicing.}

\item{pattern}{Language to define branching for a target.
For example, in a pipeline with numeric vector targets \code{x} and \code{y},
\code{tar_target(z, x + y, pattern = map(x, y))} implicitly defines
branches of \code{z} that each compute \code{x[1] + y[1]}, \code{x[2] + y[2]},
and so on. See the user manual for details.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A new target object to combine the return values
from the upstream targets.
See the "Target objects" section for background.
}
\description{
Aggregate the results of upstream targets
into a new target.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  target1 <- targets::tar_target(x, head(mtcars))
  target2 <- targets::tar_target(y, tail(mtcars))
  target3 <- tarchetypes::tar_combine(
    new_target_name,
    target1,
    target2,
    command = bind_rows(!!!.x)
  )
  list(target1, target2, target3)
})
targets::tar_manifest()
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_knitr_deps.R
\name{tar_knitr_deps}
\alias{tar_knitr_deps}
\title{List literate programming dependencies.}
\usage{
tar_knitr_deps(path)
}
\arguments{
\item{path}{Character vector, path to one or more R Markdown or
\code{knitr} reports.}
}
\value{
Character vector of the names of targets
that are dependencies of the \code{knitr} report.
}
\description{
List the target dependencies of one or more
literate programming reports (R Markdown or \code{knitr}).
}
\examples{
lines <- c(
  "---",
  "title: report",
  "output_format: html_document",
  "---",
  "",
  "```{r}",
  "targets::tar_load(data1)",
  "targets::tar_read(data2)",
  "```"
)
report <- tempfile()
writeLines(lines, report)
tar_knitr_deps(report)
}
\seealso{
Other Literate programming utilities: 
\code{\link{tar_knitr_deps_expr}()}
}
\concept{Literate programming utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_age.R
\name{tar_age}
\alias{tar_age}
\title{Create a target that runs when the last run gets old}
\usage{
tar_age(
  name,
  command,
  age,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Character of length 1, name of the target.}

\item{command}{R code to run the target and return a value.}

\item{age}{A \code{difftime} object of length 1, such as
\code{as.difftime(3, units = "days")}. If the target's output data
files are older than \code{age} (according to the most recent
time stamp over all the target's output files)
then the target will rerun.
On the other hand, if at least one data file is
younger than \code{Sys.time() - age}, then the ordinary
invalidation rules apply, and the target may or not rerun.
If you want to force the target to run every 3 days,
for example, set \code{age = as.difftime(3, units = "days")}.}

\item{pattern}{Language to define branching for a target.
For example, in a pipeline with numeric vector targets \code{x} and \code{y},
\code{tar_target(z, x + y, pattern = map(x, y))} implicitly defines
branches of \code{z} that each compute \code{x[1] + y[1]}, \code{x[2] + y[2]},
and so on. See the user manual for details.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Logical, whether to rerun the target if the user-specified
storage format changed. The storage format is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{iteration}{Logical, whether to rerun the target if the user-specified
iteration method changed. The iteration method is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{A \code{targets::tar_cue()} object. (See the "Cue objects"
section for background.) This cue object should contain any
optional secondary invalidation rules, anything except
the \code{mode} argument. \code{mode} will be automatically determined
by the \code{age} argument of \code{tar_age()}.}
}
\value{
A target object. See the "Target objects" section for background.
}
\description{
\code{tar_age()} creates a target that reruns
itself when it gets old enough.
In other words, the target reruns periodically at regular
intervals of time.
}
\details{
\code{tar_age()} uses the cue from \code{\link[=tar_cue_age]{tar_cue_age()}}, which
uses the time stamps from \code{targets::tar_meta()$time}.
See the help file of \code{targets::tar_timestamp()}
for an explanation of how this time stamp is calculated.
}
\section{Dynamic branches at regular time intervals}{

Time stamps are not recorded for whole dynamic targets,
so \code{tar_age()} is not a good fit for dynamic branching.
To invalidate dynamic branches at regular intervals,
it is recommended to use \code{targets::tar_older()} in combination
with \code{targets::tar_invalidate()} right before calling \code{tar_make()}.
For example,
\code{tar_invalidate(all_of(tar_older(Sys.time - as.difftime(1, units = "weeks"))))} # nolint
invalidates all targets more than a week old. Then, the next \code{tar_make()}
will rerun those targets.
}

\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  library(tarchetypes)
  list(
    tarchetypes::tar_age(
      data,
      data.frame(x = seq_len(26)),
      age = as.difftime(0.5, units = "secs")
    )
  )
})
targets::tar_make()
Sys.sleep(0.6)
targets::tar_make()
})
}
}
\seealso{
Other cues: 
\code{\link{tar_cue_age_raw}()},
\code{\link{tar_cue_age}()},
\code{\link{tar_cue_force}()},
\code{\link{tar_cue_skip}()}
}
\concept{cues}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map_rep_raw.R
\name{tar_map_rep_raw}
\alias{tar_map_rep_raw}
\title{Dynamic batched replication within static branches
for data frames (raw version).}
\usage{
tar_map_rep_raw(
  name,
  command,
  values = NULL,
  names = NULL,
  columns = quote(tidyselect::everything()),
  batches = 1,
  reps = 1,
  combine = TRUE,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{Language object, R code for a single replicate. Must return
a data frame.}

\item{values}{Named list or data frame with values to iterate over.
The names are the names of symbols in the commands and pattern
statements, and the elements are values that get substituted
in place of those symbols. Elements of the \code{values} list
should be small objects that can easily deparse to names,
such as characters, integers, and symbols.
For more complicated elements of \code{values}, such as
lists with multiple numeric vectors,
\code{tar_map()} attempts to parse the elements into expressions,
but this process is not perfect, and the default
target names come out garbled.
To create a list of symbols as a column of \code{values},
use \code{rlang::syms()}.}

\item{names}{Language object with a tidyselect expression
to select which columns of \code{values} to use to construct
statically branched target names. If \code{NULL}, then
short names are automatically generated.}

\item{columns}{Language object with a tidyselect expression
to select which columns of \code{values} to append to the output.
Columns already in the target output are not appended.}

\item{batches}{Number of batches. This is also the number of dynamic
branches created during \code{tar_make()}.}

\item{reps}{Number of replications in each batch. The total number
of replications is \code{batches * reps}.}

\item{combine}{Logical of length 1, whether to statically combine
all the results into a single target downstream.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the output.
An efficient data frame format like \code{"feather"} is recommended,
but the default is \code{"rds"} to avoid incurring extra package
dependencies. See the help file of \code{targets::tar_target()}
for details on storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of new target objects.
See the "Target objects" section for background.
}
\description{
Define targets for batched replication
within static branches for data frames (raw version).

This function is like \code{\link[=tar_map_rep]{tar_map_rep()}}
except the \code{name} argument is a character string
and the \code{names} and \code{columns} arguments are
language objects.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.


Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  # Just a sketch of a Bayesian sensitivity analysis of hyperparameters:
  assess_hyperparameters <- function(sigma1, sigma2) {
    # data <- simulate_random_data() # user-defined function
    # run_model(data, sigma1, sigma2) # user-defined function
    # Mock output from the model:
    posterior_samples <- stats::rnorm(1000, 0, sigma1 + sigma2)
    tibble::tibble(
      posterior_median = median(posterior_samples),
      posterior_quantile_0.025 = quantile(posterior_samples, 0.025),
      posterior_quantile_0.975 = quantile(posterior_samples, 0.975)
    )
  }
  hyperparameters <- tibble::tibble(
    scenario = c("tight", "medium", "diffuse"),
    sigma1 = c(10, 50, 50),
    sigma2 = c(10, 5, 10)
  )
  tarchetypes::tar_map_rep_raw(
    "sensitivity_analysis",
    command = quote(assess_hyperparameters(sigma1, sigma2)),
    values = hyperparameters,
    names = quote(tidyselect::any_of("scenario")),
    batches = 2,
    reps = 3
   )
})
targets::tar_make()
targets::tar_read(sensitivity_analysis)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_rep_map.R
\name{tar_rep_map}
\alias{tar_rep_map}
\title{Dynamic batched computation downstream of \code{\link[=tar_rep]{tar_rep()}} (deprecated).}
\usage{
tar_rep_map(
  name,
  command,
  ...,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{...}{Symbols to name one or more upstream batched targets
created by \code{\link[=tar_rep]{tar_rep()}}.
If you supply more than one such target, all those targets must have the
same number of batches and reps per batch. And they must all return
either data frames or lists. List targets must use \code{iteration = "list"}
in \code{\link[=tar_rep]{tar_rep()}}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A new target object to perform batched computation.
See the "Target objects" section for background.
}
\description{
Use \code{\link[=tar_rep2]{tar_rep2()}} instead.
}
\details{
Deprecated in version 0.4.0, 2021-12-06.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_rep(
      data1,
      data.frame(value = rnorm(1)),
      batches = 2,
      reps = 3
    ),
    tarchetypes::tar_rep(
      data2,
      list(value = rnorm(1)),
      batches = 2, reps = 3,
      iteration = "list" # List iteration is important for batched lists.
    ),
    tarchetypes::tar_rep2( # Use instead of tar_rep_map().
      aggregate,
      data.frame(value = data1$value + data2$value),
      data1,
      data2
    )
  )
})
targets::tar_make()
targets::tar_read(aggregate)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_download.R
\name{tar_download}
\alias{tar_download}
\title{Target that downloads URLs.}
\usage{
tar_download(
  name,
  urls,
  paths,
  method = NULL,
  quiet = TRUE,
  mode = "w",
  cacheOK = TRUE,
  extra = NULL,
  headers = NULL,
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{urls}{Character vector of URLs to track and download.
Must be known and declared before the pipeline runs.}

\item{paths}{Character vector of local file paths to
download each of the URLs.
Must be known and declared before the pipeline runs.}

\item{method}{Method to be used for downloading files.  Current
    download methods are \code{"internal"}, \code{"wininet"} (Windows
    only) \code{"libcurl"}, \code{"wget"} and \code{"curl"}, and there
    is a value \code{"auto"}: see \sQuote{Details} and \sQuote{Note}.

    The method can also be set through the option
    \code{"download.file.method"}: see \code{\link{options}()}.
  }

\item{quiet}{If \code{TRUE}, suppress status messages (if any), and
    the progress bar.}

\item{mode}{character.  The mode with which to write the file.  Useful
    values are \code{"w"}, \code{"wb"} (binary), \code{"a"} (append) and
    \code{"ab"}.  Not used for methods \code{"wget"} and \code{"curl"}.
    See also \sQuote{Details}, notably about using \code{"wb"} for Windows.
  }

\item{cacheOK}{logical.  Is a server-side cached value acceptable?}

\item{extra}{character vector of additional command-line arguments for
    the \code{"wget"} and \code{"curl"} methods.}

\item{headers}{named character vector of HTTP headers to use in HTTP
    requests.  It is ignored for non-HTTP URLs.  The \code{User-Agent}
    header, coming from the \code{HTTPUserAgent} option (see
    \code{\link{options}}) is used as the first header, automatically.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of two target objects, one upstream and one downstream.
The upstream one watches a URL for changes, and the downstream one
downloads it.
See the "Target objects" section for background.
}
\description{
Create a target that downloads file from one or more URLs
and automatically reruns when the remote data changes
(according to the ETags or last-modified time stamps).
}
\details{
\code{tar_download()} creates a pair of targets, one upstream
and one downstream. The upstream target uses \code{format = "url"}
(see \code{targets::tar_target()}) to track files at one or more URLs,
and automatically invalidate the target if the ETags
or last-modified time stamps change. The downstream target
depends on the upstream one, downloads the files,
and tracks them using \code{format = "file"}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_download(
      x,
      urls = c("https://httpbin.org/etag/test", "https://r-project.org"),
      paths = c("downloaded_file_1", "downloaded_file_2")
    )
  )
})
targets::tar_make()
targets::tar_read(x)
})
}
}
\seealso{
Other targets with custom invalidation rules: 
\code{\link{tar_change}()},
\code{\link{tar_force}()},
\code{\link{tar_skip}()}
}
\concept{targets with custom invalidation rules}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_knit_raw.R
\name{tar_knit_raw}
\alias{tar_knit_raw}
\title{Target with a knitr document (raw version).}
\usage{
tar_knit_raw(
  name,
  path,
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  error = targets::tar_option_get("error"),
  deployment = "main",
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue"),
  quiet = TRUE,
  knit_arguments = quote(list())
)
}
\arguments{
\item{name}{Character of length 1, name of the target.}

\item{path}{Character string, file path to the \code{knitr} source file.
Must have length 1.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}

\item{quiet}{Boolean; suppress the progress bar and messages?}

\item{knit_arguments}{Optional language object with a list
of named arguments to \code{knitr::knit()}.
Cannot be an expression object.
(Use \code{quote()}, not \code{expression()}.)
The reason for quoting is that these arguments may depend on
upstream targets whose values are not available at
the time the target is defined, and because \code{tar_knit_raw()}
is the "raw" version of a function, we want to avoid
all non-standard evaluation.}
}
\value{
A \code{tar_target()} object with \code{format = "file"}.
When this target runs, it returns a character vector
of file paths. The first file paths are the output files
(returned by \code{knitr::knit()}) and the knitr
source file is last. But unlike \code{knitr::knit()},
all returned paths are \emph{relative} paths to ensure portability
(so that the project can be moved from one file system to another
without invalidating the target).
See the "Target objects" section for background.
}
\description{
Shorthand to include a knitr document in a
\code{targets} pipeline (raw version)
}
\details{
\code{tar_knit_raw()} is just like \code{tar_knit()}
except that it uses standard evaluation. The \code{name} argument
is a character vector, and the \code{knit_arguments} argument
is a language object.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  # Ordinarily, you should create the report outside
  # tar_script() and avoid temporary files.
  lines <- c(
    "---",
    "title: report",
    "output_format: html_document",
    "---",
    "",
    "```{r}",
    "targets::tar_read(data)",
    "```"
  )
  path <- tempfile()
  writeLines(lines, path)
  list(
    targets::tar_target(data, data.frame(x = seq_len(26), y = letters)),
    tarchetypes::tar_knit_raw("report", path)
  )
})
targets::tar_make()
})
}
}
\seealso{
Other Literate programming targets: 
\code{\link{tar_knit}()},
\code{\link{tar_render_raw}()},
\code{\link{tar_render_rep_raw}()},
\code{\link{tar_render_rep}()},
\code{\link{tar_render}()}
}
\concept{Literate programming targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_render_raw.R
\name{tar_render_raw}
\alias{tar_render_raw}
\title{Target with an R Markdown document (raw version).}
\usage{
tar_render_raw(
  name,
  path,
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  error = targets::tar_option_get("error"),
  deployment = "main",
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue"),
  quiet = TRUE,
  render_arguments = quote(list())
)
}
\arguments{
\item{name}{Character of length 1, name of the target.}

\item{path}{Character string, file path to the R Markdown source file.
Must have length 1.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}

\item{quiet}{An option to suppress printing during rendering from knitr,
pandoc command line and others. To only suppress printing of the last
"Output created: " message, you can set \code{rmarkdown.render.message} to
\code{FALSE}}

\item{render_arguments}{Optional language object with a list
of named arguments to \code{rmarkdown::render()}.
Cannot be an expression object.
(Use \code{quote()}, not \code{expression()}.)
The reason for quoting is that these arguments may depend on
upstream targets whose values are not available at
the time the target is defined, and because \code{tar_render_raw()}
is the "raw" version of a function, we want to avoid
all non-standard evaluation.}
}
\value{
A target object with \code{format = "file"}.
When this target runs, it returns a character vector
of file paths: the rendered document, the source file,
and then the \verb{*_files/} directory if it exists.
Unlike \code{rmarkdown::render()},
all returned paths are \emph{relative} paths to ensure portability
(so that the project can be moved from one file system to another
without invalidating the target).
See the "Target objects" section for background.
}
\description{
Shorthand to include an R Markdown document in a
\code{targets} pipeline (raw version)
}
\details{
\code{tar_render_raw()} is just like \code{tar_render()}
except that it uses standard evaluation. The \code{name} argument
is a character vector, and the \code{render_arguments} argument
is a language object.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
# Unparameterized R Markdown report:
lines <- c(
  "---",
  "title: 'report.Rmd source file'",
  "output_format: html_document",
  "---",
  "Assume these lines are in report.Rmd.",
  "```{r}",
  "targets::tar_read(data)",
  "```"
)
# Include the report in the pipeline as follows:
targets::tar_script({
  library(tarchetypes)
  list(
    tar_target(data, data.frame(x = seq_len(26), y = letters)),
    tar_render_raw("report", "report.Rmd")
  )
}, ask = FALSE)
# Then, run the targets pipeline as usual.

# Parameterized R Markdown:
lines <- c(
  "---",
  "title: 'report.Rmd source file with parameters.'",
  "output_format: html_document",
  "params:",
  "  your_param: \"default value\"",
  "---",
  "Assume these lines are in report.Rmd.",
  "```{r}",
  "print(params$your_param)",
  "```"
)
# Include this parameterized report in the pipeline as follows.
targets::tar_script({
  library(tarchetypes)
  list(
    tar_target(data, data.frame(x = seq_len(26), y = letters)),
    tar_render_raw(
      "report",
      "report.Rmd",
      render_arguments = quote(list(params = list(your_param = data)))
    )
  )
}, ask = FALSE)
# Then, run the targets pipeline as usual.
})
}
}
\seealso{
Other Literate programming targets: 
\code{\link{tar_knit_raw}()},
\code{\link{tar_knit}()},
\code{\link{tar_render_rep_raw}()},
\code{\link{tar_render_rep}()},
\code{\link{tar_render}()}
}
\concept{Literate programming targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_knit.R
\name{tar_knit}
\alias{tar_knit}
\title{Target with a \code{knitr} document.}
\usage{
tar_knit(
  name,
  path,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  error = targets::tar_option_get("error"),
  deployment = "main",
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue"),
  quiet = TRUE,
  ...
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{path}{Character string, file path to the \code{knitr} source file.
Must have length 1.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}

\item{quiet}{Boolean; suppress the progress bar and messages?}

\item{...}{Named arguments to \code{knitr::knit()}.
These arguments are evaluated when the target actually runs in
\code{tar_make()}, not when the target is defined.}
}
\value{
A \code{tar_target()} object with \code{format = "file"}.
When this target runs, it returns a character vector
of file paths. The first file paths are the output files
(returned by \code{knitr::knit()}) and the \code{knitr}
source file is last. But unlike \code{knitr::knit()},
all returned paths are \emph{relative} paths to ensure portability
(so that the project can be moved from one file system to another
without invalidating the target).
See the "Target objects" section for background.
}
\description{
Shorthand to include \code{knitr} document in a
\code{targets} pipeline.
}
\details{
\code{tar_knit()} is an alternative to \code{tar_target()} for
\code{knitr} reports that depend on other targets. The \code{knitr} source
should mention dependency targets with \code{tar_load()} and \code{tar_read()}
in the active code chunks (which also allows you to knit the report
outside the pipeline if the \verb{_targets/} data store already exists).
(Do not use \code{tar_load_raw()} or \code{tar_read_raw()} for this.)
Then, \code{tar_knit()} defines a special kind of target. It
1. Finds all the \code{tar_load()}/\code{tar_read()} dependencies in the report
and inserts them into the target's command.
This enforces the proper dependency relationships.
(Do not use \code{tar_load_raw()} or \code{tar_read_raw()} for this.)
2. Sets \code{format = "file"} (see \code{tar_target()}) so \code{targets}
watches the files at the returned paths and reruns the report
if those files change.
3. Configures the target's command to return both the output
report files and the input source file. All these file paths
are relative paths so the project stays portable.
4. Forces the report to run in the user's current working directory
instead of the working directory of the report.
5. Sets convenient default options such as \code{deployment = "main"}
in the target and \code{quiet = TRUE} in \code{knitr::knit()}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  # Ordinarily, you should create the report outside
  # tar_script() and avoid temporary files.
  lines <- c(
    "---",
    "title: report",
    "output_format: html_document",
    "---",
    "",
    "```{r}",
    "targets::tar_read(data)",
    "```"
  )
  path <- tempfile()
  writeLines(lines, path)
  list(
    targets::tar_target(data, data.frame(x = seq_len(26), y = letters)),
    tarchetypes::tar_knit(report, path)
  )
})
targets::tar_make()
})
}
}
\seealso{
Other Literate programming targets: 
\code{\link{tar_knit_raw}()},
\code{\link{tar_render_raw}()},
\code{\link{tar_render_rep_raw}()},
\code{\link{tar_render_rep}()},
\code{\link{tar_render}()}
}
\concept{Literate programming targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_select_targets.R
\name{tar_select_targets}
\alias{tar_select_targets}
\title{Select target objects from a target list}
\usage{
tar_select_targets(targets, ...)
}
\arguments{
\item{targets}{A list of target objects as described in the
"Target objects" section. It does not matter how nested
the list is as long as the only leaf nodes are targets.}

\item{...}{One or more comma-separated \code{tidyselect} expressions,
e.g. \code{starts_with("prefix")}. Just like \code{...} in \code{dplyr::select()}.}
}
\value{
A list of target objects. See the "Target objects" section
of this help file.
}
\description{
Select target objects from a target list.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets <- list(
  list(
    targets::tar_target(x, 1),
    targets::tar_target(y1, 2)
  ),
  targets::tar_target(y2, 3),
  targets::tar_target(z, 4)
)
tar_select_targets(targets, starts_with("y"), contains("z"))
})
}
}
\seealso{
Other target selection: 
\code{\link{tar_select_names}()}
}
\concept{target selection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_counter.R
\name{counter_set_names}
\alias{counter_set_names}
\title{Add data to an existing counter object.}
\usage{
counter_set_names(counter, names)
}
\arguments{
\item{counter}{A counter object, defined for internal purposes only.}

\item{names}{Character vector of names to add to the counter.}
}
\value{
\code{NULL} (invisibly)
}
\description{
Not a user-side function. Do not invoke directly.
}
\examples{
counter <- counter_init()
counter_set_names(counter, letters)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_select.R
\name{tar_group_select}
\alias{tar_group_select}
\title{Group a data frame target with \code{tidyselect} semantics.}
\usage{
tar_group_select(
  name,
  command,
  by = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{by}{Tidyselect semantics to specify variables to group over.
Alternatively, you can supply a character vector.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A target object to generate a grouped data frame
to allows downstream dynamic targets to branch over the
groups of rows.
See the "Target objects" section for background.
}
\description{
Create a target that outputs a grouped data frame
with \code{dplyr::group_by()} and \code{targets::tar_group()}.
Unlike \code{tar_group_by()}, \code{tar_group_select()}
expects you to select grouping variables using \code{tidyselect} semantics.
Downstream dynamic branching targets will iterate over the groups of rows.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  produce_data <- function() {
    expand.grid(var1 = c("a", "b"), var2 = c("c", "d"), rep = c(1, 2, 3))
  }
  list(
    tarchetypes::tar_group_select(data, produce_data(), starts_with("var")),
    tar_target(group, data, pattern = map(data))
  )
})
targets::tar_make()
# Read the first row group:
targets::tar_read(group, branches = 1)
# Read the second row group:
targets::tar_read(group, branches = 2)
})
}
}
\seealso{
Other Grouped data frame targets: 
\code{\link{tar_group_by}()},
\code{\link{tar_group_count}()},
\code{\link{tar_group_size}()}
}
\concept{Grouped data frame targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_ast.R
\name{walk_ast}
\alias{walk_ast}
\title{Static code analysis for \code{tarchetypes}.}
\usage{
walk_ast(expr, walk_call)
}
\arguments{
\item{expr}{A language object or function to scan.}

\item{walk_call}{A function to handle a specific kind of function call
relevant to the code analysis at hand.}
}
\value{
A character vector of data found during static code analysis.
}
\description{
Walk an abstract syntax tree and capture data.
}
\details{
For internal use only. Not a user-side function.
Powers functionality like automatic detection of \code{tar_load()}/\code{tar_read()}
dependencies in \code{\link[=tar_render]{tar_render()}}.
Packages \code{codetools} and \code{CodeDepends} have different (more sophisticated
and elaborate) implementations of the concepts documented at
\url{https://adv-r.hadley.nz/expressions.html#ast-funs}.
}
\examples{
# How tar_render() really works:
expr <- quote({
  if (a > 1) {
    tar_load(target_name)
  }
  process_stuff(target_name)
})
walk_ast(expr, walk_call_knitr)
# Custom code analysis for developers of tarchetypes internals:
walk_custom <- function(expr, counter) {
  # New internals should use targets::tar_deparse_safe(backtick = FALSE).
  name <- deparse(expr[[1]])
  if (identical(name, "detect_this")) {
    counter_set_names(counter, as.character(expr[[2]]))
  }
}
expr <- quote({
  for (i in seq_len(10)) {
    for (j in seq_len(20)) {
      if (i > 1) {
        detect_this("prize")
      } else {
        ignore_this("penalty")
      }
    }
  }
})
walk_ast(expr, walk_custom)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_cue_skip.R
\name{tar_cue_skip}
\alias{tar_cue_skip}
\title{Cue to skip a target if a condition is true}
\usage{
tar_cue_skip(
  condition,
  command = TRUE,
  depend = TRUE,
  format = TRUE,
  iteration = TRUE,
  file = TRUE
)
}
\arguments{
\item{condition}{Logical vector evaluated locally when the target is
defined. If any element of \code{condition} is \code{TRUE}, the pipeline
will skip the target unless the target has never been built before.
If all elements of \code{condition} are \code{FALSE}, then
the target may or may not rerun, depending
on the other invalidation rules. \code{condition} is evaluated
when this cue factory is called, so the condition cannot
depend on upstream targets, and it should be quick to calculate.}

\item{command}{Logical, whether to rerun the target if command changed
since last time.}

\item{depend}{Logical, whether to rerun the target if the value of one
of the dependencies changed.}

\item{format}{Logical, whether to rerun the target if the user-specified
storage format changed. The storage format is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{iteration}{Logical, whether to rerun the target if the user-specified
iteration method changed. The iteration method is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{file}{Logical, whether to rerun the target if the file(s) with the
return value changed or at least one is missing.}
}
\value{
A cue object. See the "Cue objects" section for background.
}
\description{
\code{tar_cue_skip()} creates a cue object to
skip a target if an arbitrary condition evaluates to \code{TRUE}.
The target still builds if it was never built before.
Supply the returned cue object to the \code{cue} argument of
\code{targets::tar_target()} or similar.
}
\section{Cue objects}{

A cue object is an object generated by \code{targets::tar_cue()},
\code{tarchetypes::tar_cue_force()}, or similar. It is
a collection of decision rules that decide when a target
is invalidated/outdated (e.g. when \code{tar_make()} or similar
reruns the target). You can supply these cue objects to the
\code{tar_target()} function or similar. For example,
\code{tar_target(x, run_stuff(), cue = tar_cue(mode = "always"))}
is a target that always calls \code{run_stuff()} during \code{tar_make()}
and always shows as invalidated/outdated in \code{tar_outdated()},
\code{tar_visnetwork()}, and similar functions.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  library(tarchetypes)
  list(
    targets::tar_target(
      data,
      data.frame(x = seq_len(26)),
      cue = tarchetypes::tar_cue_skip(1 > 0)
    )
  )
})
targets::tar_make()
targets::tar_script({
  library(tarchetypes)
  list(
    targets::tar_target(
      data,
      data.frame(x = seq_len(25)), # Change the command.
      cue = tarchetypes::tar_cue_skip(1 > 0)
    )
  )
})
targets::tar_make()
targets::tar_make()
})
}
}
\seealso{
Other cues: 
\code{\link{tar_age}()},
\code{\link{tar_cue_age_raw}()},
\code{\link{tar_cue_age}()},
\code{\link{tar_cue_force}()}
}
\concept{cues}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_select_names.R
\name{tar_select_names}
\alias{tar_select_names}
\title{Select target names from a target list}
\usage{
tar_select_names(targets, ...)
}
\arguments{
\item{targets}{A list of target objects as described in the
"Target objects" section. It does not matter how nested
the list is as long as the only leaf nodes are targets.}

\item{...}{One or more comma-separated \code{tidyselect} expressions,
e.g. \code{starts_with("prefix")}. Just like \code{...} in \code{dplyr::select()}.}
}
\value{
A character vector of target names.
}
\description{
Select the names of targets from a target list.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets <- list(
  list(
    targets::tar_target(x, 1),
    targets::tar_target(y1, 2)
  ),
  targets::tar_target(y2, 3),
  targets::tar_target(z, 4)
)
tar_select_names(targets, starts_with("y"), contains("z"))
})
}
}
\seealso{
Other target selection: 
\code{\link{tar_select_targets}()}
}
\concept{target selection}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map2_count_raw.R
\name{tar_map2_count_raw}
\alias{tar_map2_count_raw}
\title{Dynamic-within-static branching for data frames
(count batching; raw version).}
\usage{
tar_map2_count_raw(
  name,
  command1,
  command2,
  values = NULL,
  names = NULL,
  batches = 1L,
  combine = TRUE,
  suffix1 = "1",
  suffix2 = "2",
  columns1 = quote(tidyselect::everything()),
  columns2 = quote(tidyselect::everything()),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Character of length 1, base name of the targets.}

\item{command1}{Language object to create named arguments to \code{command2}.
Must return a data frame with one row per call to \code{command2}.}

\item{command2}{Language object  to map over the data frame of arguments
produced by \code{command1}. Must return a data frame.}

\item{values}{Named list or data frame with values to iterate over.
The names are the names of symbols in the commands and pattern
statements, and the elements are values that get substituted
in place of those symbols. Elements of the \code{values} list
should be small objects that can easily deparse to names,
such as characters, integers, and symbols.
For more complicated elements of \code{values}, such as
lists with multiple numeric vectors,
\code{tar_map()} attempts to parse the elements into expressions,
but this process is not perfect, and the default
target names come out garbled.
To create a list of symbols as a column of \code{values},
use \code{rlang::syms()}.}

\item{names}{Language object with a tidyselect expression
to select which columns of \code{values} to use to construct
statically branched target names. If \code{NULL}, then
short names are automatically generated.}

\item{batches}{Positive integer of length 1,
maximum number of batches (dynamic branches within static branches)
of the downstream (\code{command2}) targets. Batches
are formed from row groups of the \code{command1} target output.}

\item{combine}{Logical of length 1, whether to statically combine
all the results into a single target downstream.}

\item{suffix1}{Character of length 1,
suffix to apply to the \code{command1} targets to distinguish
them from the \code{command2} targets.}

\item{suffix2}{Character of length 1,
suffix to apply to the \code{command2} targets to distinguish
them from the \code{command1} targets.}

\item{columns1}{Language object, a tidyselect expression
to select which columns of \code{values}
to append to the output of all targets.}

\item{columns2}{Language object, a tidyselect expression
to select which columns of \code{command1}
output to append to \code{command2} output.
In case of conflicts, \code{column1} takes precedence.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the output.
An efficient data frame format like \code{"feather"} is recommended,
but the default is \code{"rds"} to avoid incurring extra package
dependencies. See the help file of \code{targets::tar_target()}
for details on storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of new target objects.
See the "Target objects" section for background.
}
\description{
Define targets for batched
dynamic-within-static branching for data frames,
where the user sets the (maximum) number of batches.
Like \code{tar_map2_count()} except \code{name} is a character string
and \code{command1}, \code{command2}, \code{names}, \code{columns1}, and \code{columns2}
are all language objects.
}
\details{
Static branching creates one pair of targets
for each row in \code{values}. In each pair,
there is an upstream non-dynamic target that runs \code{command1}
and a downstream dynamic target that runs \code{command2}.
\code{command1} produces a data frame of arguments to
\code{command2}, and \code{command2} dynamically maps over
these arguments in batches.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  tarchetypes::tar_map2_count_raw(
    "x",
    command1 = quote(
      tibble::tibble(
        arg1 = arg1,
        arg2 = seq_len(6)
       )
    ),
    command2 = quote(
      tibble::tibble(
        result = paste(arg1, arg2),
        random = sample.int(1e6, size = 1),
        length_input = length(arg1)
      )
    ),
    values = tibble::tibble(arg1 = letters[seq_len(2)]),
    batches = 3
   )
})
targets::tar_make()
targets::tar_read(x)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_knitr_deps_expr.R
\name{tar_knitr_deps_expr}
\alias{tar_knitr_deps_expr}
\title{Expression with literate programming dependencies.}
\usage{
tar_knitr_deps_expr(path)
}
\arguments{
\item{path}{Character vector, path to one or more R Markdown or
\code{knitr} reports.}
}
\value{
Expression object to name the dependency targets
of the \code{knitr} report, which will be detected in the
static code analysis of \code{targets}.
}
\description{
Construct an expression whose global variable dependencies
are the target dependencies of one or more literate programming reports
(R Markdown or \code{knitr}). This helps third-party developers create their
own third-party target factories for literate programming targets
(similar to \code{\link[=tar_knit]{tar_knit()}} and \code{\link[=tar_render]{tar_render()}}).
}
\examples{
lines <- c(
  "---",
  "title: report",
  "output_format: html_document",
  "---",
  "",
  "```{r}",
  "targets::tar_load(data1)",
  "targets::tar_read(data2)",
  "```"
)
report <- tempfile()
writeLines(lines, report)
tar_knitr_deps_expr(report)
}
\seealso{
Other Literate programming utilities: 
\code{\link{tar_knitr_deps}()}
}
\concept{Literate programming utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_rep2.R
\name{tar_rep2}
\alias{tar_rep2}
\title{Dynamic batched computation downstream of \code{\link[=tar_rep]{tar_rep()}}}
\usage{
tar_rep2(
  name,
  command,
  ...,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{...}{Symbols to name one or more upstream batched targets
created by \code{\link[=tar_rep]{tar_rep()}}.
If you supply more than one such target, all those targets must have the
same number of batches and reps per batch. And they must all return
either data frames or lists. List targets must use \code{iteration = "list"}
in \code{\link[=tar_rep]{tar_rep()}}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A new target object to perform batched computation.
See the "Target objects" section for background.
}
\description{
Batching is important for optimizing the efficiency
of heavily dynamically-branched workflows:
\url{https://books.ropensci.org/targets/dynamic.html#batching}.
\code{\link[=tar_rep2]{tar_rep2()}} uses dynamic branching to iterate
over the batches and reps of existing upstream targets.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_rep(
      data1,
      data.frame(value = rnorm(1)),
      batches = 2,
      reps = 3
    ),
    tarchetypes::tar_rep(
      data2,
      list(value = rnorm(1)),
      batches = 2, reps = 3,
      iteration = "list" # List iteration is important for batched lists.
    ),
    tarchetypes::tar_rep2(
      aggregate,
      data.frame(value = data1$value + data2$value),
      data1,
      data2
    )
  )
})
targets::tar_make()
targets::tar_read(aggregate)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_size.R
\name{tar_group_size_run}
\alias{tar_group_size_run}
\title{Generate a grouped data frame within tar_group_size()}
\usage{
tar_group_size_run(data, size)
}
\arguments{
\item{data}{A data frame to group.}
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_sub.R
\name{tar_sub}
\alias{tar_sub}
\title{Create multiple expressions with symbol substitution.}
\usage{
tar_sub(expr, values)
}
\arguments{
\item{expr}{Starting expression. Values are iteratively substituted
in place of symbols in \code{expr} to create each new expression.}

\item{values}{List of values to substitute into \code{expr} to create
the expressions. All elements of \code{values} must have the same length.}
}
\value{
A list of expression objects. Often, these expression objects
evaluate to target objects (but not necessarily).
See the "Target objects" section for background.
}
\description{
Loop over a grid of values and create an expression object
from each one. Helps with general metaprogramming.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
# tar_map() is incompatible with tar_render() because the latter
# operates on preexisting tar_target() objects. By contrast,
# tar_eval() and tar_sub() iterate over code farther upstream.
values <- list(
  name = lapply(c("name1", "name2"), as.symbol),
  file = list("file1.Rmd", "file2.Rmd")
)
tar_sub(tar_render(name, file), values = values)
}
\seealso{
Other Metaprogramming utilities: 
\code{\link{tar_eval_raw}()},
\code{\link{tar_eval}()},
\code{\link{tar_sub_raw}()}
}
\concept{Metaprogramming utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_counter.R
\name{counter_init}
\alias{counter_init}
\title{Counter constructor.}
\usage{
counter_init(names = NULL)
}
\arguments{
\item{names}{Character vector of names to add to the new counter.}
}
\value{
A new counter object.
}
\description{
Not a user-side function. Do not invoke directly.
}
\details{
Creates a counter object as described at
\url{https://books.ropensci.org/targets-design/classes.html#counter-class}.
}
\examples{
counter <- counter_init()
counter_set_names(counter, letters)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_force.R
\name{tar_force}
\alias{tar_force}
\title{Target with a custom condition to force execution.}
\usage{
tar_force(
  name,
  command,
  force,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{force}{R code for the condition that forces a build.
If it evaluates to \code{TRUE}, then your work will run during \code{tar_make()}.}

\item{tidy_eval}{Whether to invoke tidy evaluation
(e.g. the \verb{!!} operator from \code{rlang}) as soon as the target is defined
(before \code{tar_make()}). Applies to arguments \code{command} and \code{force}.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()}
to customize the rules that decide whether the target is up to date.
Only applies to the downstream target. The upstream target always runs.}
}
\value{
A list of 2 targets objects: one to indicate whether the custom
condition is met, and another to respond to it and do your
actual work.
See the "Target objects" section for background.
}
\description{
Create a target that always runs if a user-defined
condition rule is met.
}
\details{
\code{tar_force()} creates a target that always runs
when a custom condition is met. The implementation builds
on top of \code{\link[=tar_change]{tar_change()}}. Thus, a pair of targets is created:
an upstream auxiliary target to indicate the custom condition
and a downstream target that responds to it and does your work.

\code{tar_force()} does not actually use \code{\link[=tar_cue_force]{tar_cue_force()}}, and the
mechanism is totally different.
Because the upstream target always runs,
\code{tar_outdated()} and \code{tar_visnetwork()} will always
show both targets as outdated. However, \code{tar_make()} will still
skip the downstream one if the upstream custom condition is not met.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_force(x, tempfile(), force = 1 > 0)
  )
})
targets::tar_make()
targets::tar_make()
})
}
}
\seealso{
Other targets with custom invalidation rules: 
\code{\link{tar_change}()},
\code{\link{tar_download}()},
\code{\link{tar_skip}()}
}
\concept{targets with custom invalidation rules}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map2_raw.R
\name{tar_map2_group}
\alias{tar_map2_group}
\title{Append the \code{tar_group} variable to a \code{tar_map2()} target.}
\usage{
tar_map2_group(data, group)
}
\arguments{
\item{data}{Data frame to be returned from the target.}

\item{group}{Function on the data to return the \code{tar_group} column.
If \code{group} is \code{NULL}, then no \code{tar_group} column is attached.}
}
\value{
A data frame with a \code{tar_group} column attached (if \code{group}
is not \code{NULL}).
}
\description{
Append the \code{tar_group} variable to a \code{tar_map2()} target.
}
\details{
For internal use only. Users should not invoke
this function directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_hook_before.R
\name{tar_hook_before}
\alias{tar_hook_before}
\title{Hook to prepend code}
\usage{
tar_hook_before(targets, hook, names = NULL)
}
\arguments{
\item{targets}{A list of target objects. The input target list
can be arbitrarily nested, but it must consist entirely of target
objects. In addition, the return value is a simple list
where each element is a target object.
All hook functions remove the nested structure of the input target list.}

\item{hook}{R code to insert. When you supply code to this argument,
the code is quoted (not evaluated) so there is no need
to wrap it in \code{quote()}, \code{expression()}, or similar.}

\item{names}{Name of targets in the target list
to apply the hook. You can supply symbols, a character vector,
or tidyselect helpers like \code{\link[=starts_with]{starts_with()}}.
Targets not included in \code{names} still remain in the target list,
but they are not modified because the hook does not apply to them.}
}
\value{
A flattened list of target objects with the hooks applied.
Even if the input target list had a nested structure,
the return value is a simple list where each element is a target object.
All hook functions remove the nested structure of the input target list.
}
\description{
Prepend R code to the commands of multiple targets.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  targets <- list(
    # Nested target lists work with hooks.
    list(
      targets::tar_target(x1, task1()),
      targets::tar_target(x2, task2(x1))
    ),
    targets::tar_target(x3, task3(x2)),
    targets::tar_target(y1, task4(x3))
  )
  tarchetypes::tar_hook_before(
    targets = targets,
    hook = print("Running hook."),
    names = starts_with("x")
  )
})
targets::tar_manifest(fields = command)
})
}
}
\seealso{
Other hooks: 
\code{\link{tar_hook_inner}()},
\code{\link{tar_hook_outer}()}
}
\concept{hooks}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_skip.R
\name{tar_skip}
\alias{tar_skip}
\title{Target with a custom cancellation condition.}
\usage{
tar_skip(
  name,
  command,
  skip,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{skip}{R code for the skipping condition. If it evaluates to \code{TRUE}
during \code{tar_make()}, the target will cancel itself.}

\item{tidy_eval}{Whether to invoke tidy evaluation
(e.g. the \verb{!!} operator from \code{rlang}) as soon as the target is defined
(before \code{tar_make()}). Applies to arguments \code{command} and \code{skip}.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A target object with \code{targets::tar_cancel(your_condition)} inserted
into the command.
See the "Target objects" section for background.
}
\description{
Create a target that cancels itself if a user-defined
decision rule is met.
}
\details{
\code{tar_skip()} creates a target that cancels itself
whenever a custom condition is met. The mechanism of cancellation
is \code{targets::tar_cancel(your_condition)}, which allows skipping to happen
even if the target does not exist yet. This behavior differs from
\code{tar_cue(mode = "never")}, which still runs if the target does not exist.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_skip(x, command = "value", skip = 1 > 0)
  )
})
targets::tar_make()
})
}
}
\seealso{
Other targets with custom invalidation rules: 
\code{\link{tar_change}()},
\code{\link{tar_download}()},
\code{\link{tar_force}()}
}
\concept{targets with custom invalidation rules}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_file_read.R
\name{tar_file_read}
\alias{tar_file_read}
\title{Track a file and read the contents.}
\usage{
tar_file_read(
  name,
  command,
  read,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code that runs in the \code{format = "file"} target
and returns the file to be tracked.}

\item{read}{R code to read the file. Must include \code{!!.x}
where the file path goes: for example,
\code{read = readr::read_csv(file = !!.x, col_types = readr::cols())}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of two new target objects to track a file
and read the contents.
See the "Target objects" section for background.
}
\description{
Create a pair of targets: one to
track a file with \code{format = "file"}, and another
to read the file.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  tar_file_read(data, get_path(), read_csv(file = !!.x, col_types = cols()))
})
targets::tar_manifest()
})
}
}
\concept{Simple files}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_count.R
\name{tar_group_count_run}
\alias{tar_group_count_run}
\title{Generate a grouped data frame within \code{tar_group_count()}.}
\usage{
tar_group_count_run(data, count)
}
\arguments{
\item{data}{A data frame to group.}

\item{count}{Maximum number of groups.}
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_download.R
\name{tar_download_run}
\alias{tar_download_run}
\title{Download multiple URLs and return the local paths.}
\usage{
tar_download_run(urls, paths, method, quiet, mode, cacheOK, extra, headers)
}
\arguments{
\item{urls}{Character vector of URLs to track and download.
Must be known and declared before the pipeline runs.}

\item{paths}{Character vector of local file paths to
download each of the URLs.
Must be known and declared before the pipeline runs.}

\item{method}{Method to be used for downloading files.  Current
    download methods are \code{"internal"}, \code{"wininet"} (Windows
    only) \code{"libcurl"}, \code{"wget"} and \code{"curl"}, and there
    is a value \code{"auto"}: see \sQuote{Details} and \sQuote{Note}.

    The method can also be set through the option
    \code{"download.file.method"}: see \code{\link{options}()}.
  }

\item{quiet}{If \code{TRUE}, suppress status messages (if any), and
    the progress bar.}

\item{mode}{character.  The mode with which to write the file.  Useful
    values are \code{"w"}, \code{"wb"} (binary), \code{"a"} (append) and
    \code{"ab"}.  Not used for methods \code{"wget"} and \code{"curl"}.
    See also \sQuote{Details}, notably about using \code{"wb"} for Windows.
  }

\item{cacheOK}{logical.  Is a server-side cached value acceptable?}

\item{extra}{character vector of additional command-line arguments for
    the \code{"wget"} and \code{"curl"} methods.}

\item{headers}{named character vector of HTTP headers to use in HTTP
    requests.  It is ignored for non-HTTP URLs.  The \code{User-Agent}
    header, coming from the \code{HTTPUserAgent} option (see
    \code{\link{options}}) is used as the first header, automatically.}
}
\value{
A character vector of file paths where the URLs
were downloaded.
}
\description{
Not a user-side function. Do not invoke directly.
}
\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
  tarchetypes::tar_download_run(
    urls = "https://httpbin.org/etag/test",
    paths = tempfile(),
    method = NULL,
    quiet = TRUE,
    mode = "w",
    cacheOK = NULL,
    extra = NULL,
    headers = NULL
  )
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map.R
\name{tar_map}
\alias{tar_map}
\title{Static branching.}
\usage{
tar_map(values, ..., names = tidyselect::everything(), unlist = FALSE)
}
\arguments{
\item{values}{Named list or data frame with values to iterate over.
The names are the names of symbols in the commands and pattern
statements, and the elements are values that get substituted
in place of those symbols. Elements of the \code{values} list
should be small objects that can easily deparse to names,
such as characters, integers, and symbols.
For more complicated elements of \code{values}, such as
lists with multiple numeric vectors,
\code{tar_map()} attempts to parse the elements into expressions,
but this process is not perfect, and the default
target names come out garbled.
To create a list of symbols as a column of \code{values},
use \code{rlang::syms()}.}

\item{...}{One or more target objects or list of target objects.
Lists can be arbitrarily nested, as in \code{list()}.}

\item{names}{Subset of \code{names(values)}
used to generate the suffixes in the names of the new targets.
You can supply symbols, a character vector,
or tidyselect helpers like \code{\link[=starts_with]{starts_with()}}.}

\item{unlist}{Logical, whether to flatten the returned list of targets.
If \code{unlist = FALSE}, the list is nested and sub-lists
are named and grouped by the original input targets.
If \code{unlist = TRUE}, the return value is a flat list of targets
named by the new target names.}
}
\value{
A list of new target objects. If \code{unlist} is \code{FALSE},
the list is nested and sub-lists are named and grouped by the original
input targets. If \code{unlist = TRUE}, the return value is a flat list of
targets named by the new target names.
See the "Target objects" section for background.
}
\description{
Define multiple new targets based on existing target objects.
}
\details{
\code{tar_map()} creates collections of new
targets by iterating over a list of arguments
and substituting symbols into commands and pattern statements.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_map(
      list(a = c(12, 34), b = c(45, 78)),
      targets::tar_target(x, a + b),
      targets::tar_target(y, x + a, pattern = map(x))
    )
  )
})
targets::tar_manifest()
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_cue_force.R
\name{tar_cue_force}
\alias{tar_cue_force}
\title{Cue to force a target to run if a condition is true}
\usage{
tar_cue_force(
  condition,
  command = TRUE,
  depend = TRUE,
  format = TRUE,
  iteration = TRUE,
  file = TRUE
)
}
\arguments{
\item{condition}{Logical vector evaluated locally when the target is
defined. If any element of \code{condition} is \code{TRUE}, the target will
definitely rerun when the pipeline runs.
Otherwise, the target may or may not rerun, depending
on the other invalidation rules. \code{condition} is evaluated
when this cue factory is called, so the condition cannot
depend on upstream targets, and it should be quick to calculate.}

\item{command}{Logical, whether to rerun the target if command changed
since last time.}

\item{depend}{Logical, whether to rerun the target if the value of one
of the dependencies changed.}

\item{format}{Logical, whether to rerun the target if the user-specified
storage format changed. The storage format is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{iteration}{Logical, whether to rerun the target if the user-specified
iteration method changed. The iteration method is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{file}{Logical, whether to rerun the target if the file(s) with the
return value changed or at least one is missing.}
}
\value{
A cue object. See the "Cue objects" section for background.
}
\description{
\code{tar_cue_force()} creates a cue object to
force a target to run if an arbitrary condition evaluates to \code{TRUE}.
Supply the returned cue object to the \code{cue} argument of
\code{targets::tar_target()} or similar.
}
\details{
\code{tar_cue_force()} and \code{\link[=tar_force]{tar_force()}} operate differently.
The former defines a cue object based on an eagerly evaluated
condition, and \code{\link[=tar_force]{tar_force()}} puts the condition in a special
upstream target that always runs. Unlike \code{tar_cue_force()},
the condition in \code{\link[=tar_force]{tar_force()}} can depend on upstream targets,
but the drawback is that targets defined with \code{\link[=tar_force]{tar_force()}}
will always show up as outdated in functions like \code{tar_outdated()}
and \code{tar_visnetwork()} even though \code{tar_make()} may still
skip the main target if the condition is not met.
}
\section{Cue objects}{

A cue object is an object generated by \code{targets::tar_cue()},
\code{tarchetypes::tar_cue_force()}, or similar. It is
a collection of decision rules that decide when a target
is invalidated/outdated (e.g. when \code{tar_make()} or similar
reruns the target). You can supply these cue objects to the
\code{tar_target()} function or similar. For example,
\code{tar_target(x, run_stuff(), cue = tar_cue(mode = "always"))}
is a target that always calls \code{run_stuff()} during \code{tar_make()}
and always shows as invalidated/outdated in \code{tar_outdated()},
\code{tar_visnetwork()}, and similar functions.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  library(tarchetypes)
  list(
    targets::tar_target(
      data,
      data.frame(x = seq_len(26)),
      cue = tarchetypes::tar_cue_force(1 > 0)
    )
  )
})
targets::tar_make()
targets::tar_make()
})
}
}
\seealso{
Other cues: 
\code{\link{tar_age}()},
\code{\link{tar_cue_age_raw}()},
\code{\link{tar_cue_age}()},
\code{\link{tar_cue_skip}()}
}
\concept{cues}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map2_count.R
\name{tar_map2_count}
\alias{tar_map2_count}
\title{Dynamic-within-static branching for data frames
(count batching).}
\usage{
tar_map2_count(
  name,
  command1,
  command2,
  values = NULL,
  names = NULL,
  batches = 1L,
  combine = TRUE,
  suffix1 = "1",
  suffix2 = "2",
  columns1 = tidyselect::everything(),
  columns2 = tidyselect::everything(),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name of the targets.}

\item{command1}{R code to create named arguments to \code{command2}.
Must return a data frame with one row per call to \code{command2}.}

\item{command2}{R code to map over the data frame of arguments
produced by \code{command1}. Must return a data frame.}

\item{values}{Named list or data frame with values to iterate over.
The names are the names of symbols in the commands and pattern
statements, and the elements are values that get substituted
in place of those symbols. Elements of the \code{values} list
should be small objects that can easily deparse to names,
such as characters, integers, and symbols.
For more complicated elements of \code{values}, such as
lists with multiple numeric vectors,
\code{tar_map()} attempts to parse the elements into expressions,
but this process is not perfect, and the default
target names come out garbled.
To create a list of symbols as a column of \code{values},
use \code{rlang::syms()}.}

\item{names}{Language object with a tidyselect expression
to select which columns of \code{values} to use to construct
statically branched target names. If \code{NULL}, then
short names are automatically generated.}

\item{batches}{Positive integer of length 1,
maximum number of batches (dynamic branches within static branches)
of the downstream (\code{command2}) targets. Batches
are formed from row groups of the \code{command1} target output.}

\item{combine}{Logical of length 1, whether to statically combine
all the results into a single target downstream.}

\item{suffix1}{Character of length 1,
suffix to apply to the \code{command1} targets to distinguish
them from the \code{command2} targets.}

\item{suffix2}{Character of length 1,
suffix to apply to the \code{command2} targets to distinguish
them from the \code{command1} targets.}

\item{columns1}{A tidyselect expression to select which columns of \code{values}
to append to the output of all targets.
Columns already in the target output are not appended.}

\item{columns2}{A tidyselect expression to select which columns of \code{command1}
output to append to \code{command2} output.
Columns already in the target output are not appended.
\code{columns1} takes precedence over \code{columns2}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the output.
An efficient data frame format like \code{"feather"} is recommended,
but the default is \code{"rds"} to avoid incurring extra package
dependencies. See the help file of \code{targets::tar_target()}
for details on storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of new target objects.
See the "Target objects" section for background.
}
\description{
Define targets for batched
dynamic-within-static branching for data frames,
where the user sets the (maximum) number of batches.
}
\details{
Static branching creates one pair of targets
for each row in \code{values}. In each pair,
there is an upstream non-dynamic target that runs \code{command1}
and a downstream dynamic target that runs \code{command2}.
\code{command1} produces a data frame of arguments to
\code{command2}, and \code{command2} dynamically maps over
these arguments in batches.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  tarchetypes::tar_map2_count(
    x,
    command1 = tibble::tibble(
      arg1 = arg1,
      arg2 = seq_len(6)
     ),
    command2 = tibble::tibble(
      result = paste(arg1, arg2),
      random = sample.int(1e9, size = 1),
      length_input = length(arg1)
    ),
    values = tibble::tibble(arg1 = letters[seq_len(2)]),
    batches = 3
   )
})
targets::tar_make()
targets::tar_read(x)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_eval_raw.R
\name{tar_eval_raw}
\alias{tar_eval_raw}
\title{Evaluate multiple expressions created with symbol substitution
(raw version).}
\usage{
tar_eval_raw(expr, values, envir = parent.frame())
}
\arguments{
\item{expr}{Expression object with the starting expression.
Values are iteratively substituted
in place of symbols in \code{expr} to create each new expression,
and then each expression is evaluated.}

\item{values}{List of values to substitute into \code{expr} to create
the expressions. All elements of \code{values} must have the same length.}

\item{envir}{Environment in which to evaluate the new expressions.}
}
\value{
A list of return values from evaluating the expression objects.
Often, these values are target objects.
See the "Target objects" section for background
on target objects specifically.
}
\description{
Loop over a grid of values, create an expression object
from each one, and then evaluate that expression.
Helps with general metaprogramming. Unlike \code{\link[=tar_sub]{tar_sub()}},
which quotes the \code{expr} argument, \code{tar_sub_raw()} assumes \code{expr}
is an expression object.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
# tar_map() is incompatible with tar_render() because the latter
# operates on preexisting tar_target() objects. By contrast,
# tar_eval_raw() and tar_sub_raw() iterate over code farther upstream.
values <- list(
  name = lapply(c("name1", "name2"), as.symbol),
  file = c("file1.Rmd", "file2.Rmd")
)
tar_sub_raw(quote(list(name, file)), values = values)
tar_sub_raw(quote(tar_render(name, file)), values = values)
path <- tempfile()
file.create(path)
str(tar_eval_raw(quote(tar_render(name, path)), values = values))
# So in your _targets.R file, you can define a pipeline like as below.
# Just make sure to set a unique name for each target
# (which tar_map() does automatically).
values <- list(
  name = lapply(c("name1", "name2"), as.symbol),
  file = c(path, path)
)
list(
  tar_eval_raw(quote(tar_render(name, file)), values = values)
)
}
\seealso{
Other Metaprogramming utilities: 
\code{\link{tar_eval}()},
\code{\link{tar_sub_raw}()},
\code{\link{tar_sub}()}
}
\concept{Metaprogramming utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_by.R
\name{tar_group_by}
\alias{tar_group_by}
\title{Group a data frame target by one or more variables.}
\usage{
tar_group_by(
  name,
  command,
  ...,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{...}{Symbols, variables in the output data frame to group by.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A target object to generate a grouped data frame
to allows downstream dynamic targets to branch over the
groups of rows.
See the "Target objects" section for background.
}
\description{
Create a target that outputs a grouped data frame
with \code{dplyr::group_by()} and \code{targets::tar_group()}. Downstream
dynamic branching targets will iterate over the groups of rows.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  produce_data <- function() {
    expand.grid(var1 = c("a", "b"), var2 = c("c", "d"), rep = c(1, 2, 3))
  }
  list(
    tarchetypes::tar_group_by(data, produce_data(), var1, var2),
    tar_target(group, data, pattern = map(data))
  )
})
targets::tar_make()
# Read the first row group:
targets::tar_read(group, branches = 1)
# Read the second row group:
targets::tar_read(group, branches = 2)
})
}
}
\seealso{
Other Grouped data frame targets: 
\code{\link{tar_group_count}()},
\code{\link{tar_group_select}()},
\code{\link{tar_group_size}()}
}
\concept{Grouped data frame targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_cue_age.R
\name{tar_cue_age}
\alias{tar_cue_age}
\title{Cue to run a target when the last output reaches a certain age}
\usage{
tar_cue_age(
  name,
  age,
  command = TRUE,
  depend = TRUE,
  format = TRUE,
  iteration = TRUE,
  file = TRUE
)
}
\arguments{
\item{name}{Symbol, name of the target.}

\item{age}{A \code{difftime} object of length 1, such as
\code{as.difftime(3, units = "days")}. If the target's output data
files are older than \code{age} (according to the most recent
time stamp over all the target's output files)
then the target will rerun.
On the other hand, if at least one data file is
younger than \code{Sys.time() - age}, then the ordinary
invalidation rules apply, and the target may or not rerun.
If you want to force the target to run every 3 days,
for example, set \code{age = as.difftime(3, units = "days")}.}

\item{command}{Logical, whether to rerun the target if command changed
since last time.}

\item{depend}{Logical, whether to rerun the target if the value of one
of the dependencies changed.}

\item{format}{Logical, whether to rerun the target if the user-specified
storage format changed. The storage format is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{iteration}{Logical, whether to rerun the target if the user-specified
iteration method changed. The iteration method is user-specified through
\code{\link[targets:tar_target]{tar_target()}} or \code{\link[targets:tar_option_set]{tar_option_set()}}.}

\item{file}{Logical, whether to rerun the target if the file(s) with the
return value changed or at least one is missing.}
}
\value{
A cue object. See the "Cue objects" section for background.
}
\description{
\code{tar_cue_age()} creates a cue object to
rerun a target if the most recent output data becomes old enough.
The age of the target is determined by \code{targets::tar_timestamp()},
and the way the time stamp is calculated is explained
in the Details section of the help file of that function.
}
\details{
\code{tar_cue_age()} uses the time stamps from \code{tar_meta()$time}.
If no time stamp is recorded, the cue defaults to the ordinary
invalidation rules (i.e. \code{mode = "thorough"} in \code{targets::tar_cue()}).
}
\section{Dynamic branches at regular time intervals}{

Time stamps are not recorded for whole dynamic targets,
so \code{tar_age()} is not a good fit for dynamic branching.
To invalidate dynamic branches at regular intervals,
it is recommended to use \code{targets::tar_older()} in combination
with \code{targets::tar_invalidate()} right before calling \code{tar_make()}.
For example,
\code{tar_invalidate(all_of(tar_older(Sys.time - as.difftime(1, units = "weeks"))))} # nolint
invalidates all targets more than a week old. Then, the next \code{tar_make()}
will rerun those targets.
}

\section{Cue objects}{

A cue object is an object generated by \code{targets::tar_cue()},
\code{tarchetypes::tar_cue_force()}, or similar. It is
a collection of decision rules that decide when a target
is invalidated/outdated (e.g. when \code{tar_make()} or similar
reruns the target). You can supply these cue objects to the
\code{tar_target()} function or similar. For example,
\code{tar_target(x, run_stuff(), cue = tar_cue(mode = "always"))}
is a target that always calls \code{run_stuff()} during \code{tar_make()}
and always shows as invalidated/outdated in \code{tar_outdated()},
\code{tar_visnetwork()}, and similar functions.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  library(tarchetypes)
  list(
    targets::tar_target(
      data,
      data.frame(x = seq_len(26)),
      cue = tarchetypes::tar_cue_age(
        name = data,
        age = as.difftime(0.5, units = "secs")
      )
    )
  )
})
targets::tar_make()
Sys.sleep(0.6)
targets::tar_make()
})
}
}
\seealso{
Other cues: 
\code{\link{tar_age}()},
\code{\link{tar_cue_age_raw}()},
\code{\link{tar_cue_force}()},
\code{\link{tar_cue_skip}()}
}
\concept{cues}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_rep_raw.R
\name{tar_rep_run}
\alias{tar_rep_run}
\title{Run a batch in a \code{tar_rep()} archetype.}
\usage{
tar_rep_run(command, batch, reps, iteration)
}
\arguments{
\item{command}{Expression object, command to replicate.}

\item{batch}{Numeric, batch index.}

\item{reps}{Numeric, number of reps per batch.}

\item{iteration}{Character, iteration method.}
}
\value{
Aggregated results of multiple executions of the
user-defined command supplied to \code{\link[=tar_rep]{tar_rep()}}. Depends on what
the user specifies. Common use cases are simulated datasets.
}
\description{
Internal function needed for \code{tar_rep()}.
Users should not invoke it directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_size.R
\name{tar_group_size}
\alias{tar_group_size}
\title{Group the rows of a data frame into groups
of a given size.}
\usage{
tar_group_size(
  name,
  command,
  size,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{size}{Positive integer, maximum number of rows in each group.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A target object to generate a grouped data frame
to allows downstream dynamic targets to branch over the
groups of rows.
See the "Target objects" section for background.
}
\description{
Create a target that outputs a grouped data frame
for downstream dynamic branching. Row groups have
the number of rows you supply to \code{size} (plus the remainder
in a group of its own, if applicable.) The total number of groups
varies.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  produce_data <- function() {
    expand.grid(var1 = c("a", "b"), var2 = c("c", "d"), rep = c(1, 2, 3))
  }
  list(
    tarchetypes::tar_group_size(data, produce_data(), size = 7),
    tar_target(group, data, pattern = map(data))
  )
})
targets::tar_make()
# Read the first row group:
targets::tar_read(group, branches = 1)
# Read the second row group:
targets::tar_read(group, branches = 2)
})
}
}
\seealso{
Other Grouped data frame targets: 
\code{\link{tar_group_by}()},
\code{\link{tar_group_count}()},
\code{\link{tar_group_select}()}
}
\concept{Grouped data frame targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_change.R
\name{tar_change}
\alias{tar_change}
\title{Target that responds to an arbitrary change.}
\usage{
tar_change(
  name,
  command,
  change,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{change}{R code for the upstream change-inducing target.}

\item{tidy_eval}{Whether to invoke tidy evaluation
(e.g. the \verb{!!} operator from \code{rlang}) as soon as the target is defined
(before \code{tar_make()}). Applies to arguments \code{command} and \code{change}.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()}
to customize the rules that decide whether the target is up to date.
Only applies to the downstream target. The upstream target always runs.}
}
\value{
A list of two target objects, one upstream and one downstream.
The upstream one triggers the change, and the downstream one
responds to it.
See the "Target objects" section for background.
}
\description{
Create a target that responds to a change
in an arbitrary value. If the value changes, the target reruns.
}
\details{
\code{tar_change()} creates a pair of targets, one upstream
and one downstream. The upstream target always runs and returns
an auxiliary value. This auxiliary value gets referenced in the
downstream target, which causes the downstream target to rerun
if the auxiliary value changes. The behavior is cancelled if
\code{cue} is \code{tar_cue(depend = FALSE)} or \code{tar_cue(mode = "never")}.

Because the upstream target always runs,
\code{tar_outdated()} and \code{tar_visnetwork()} will always
show both targets as outdated. However, \code{tar_make()} will still
skip the downstream one if the upstream target
did not detect a change.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_change(x, command = tempfile(), change = tempfile())
  )
})
targets::tar_make()
targets::tar_make()
})
}
}
\seealso{
Other targets with custom invalidation rules: 
\code{\link{tar_download}()},
\code{\link{tar_force}()},
\code{\link{tar_skip}()}
}
\concept{targets with custom invalidation rules}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_combine_raw.R
\name{tar_combine_raw}
\alias{tar_combine_raw}
\title{Static aggregation (raw version).}
\usage{
tar_combine_raw(
  name,
  ...,
  command = expression(vctrs::vec_c(!!!.x)),
  use_names = TRUE,
  pattern = NULL,
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Character, name of the new target.}

\item{...}{One or more target objects or list of target objects.
Lists can be arbitrarily nested, as in \code{list()}.}

\item{command}{Expression object,
R command to aggregate the targets. Must contain
\code{!!!.x} where the arguments are to be inserted,
where \verb{!!!} is the unquote splice operator from \code{rlang}.}

\item{use_names}{Logical, whether to insert the names of the targets
into the command when splicing.}

\item{pattern}{Similar to the \code{pattern} argument of \code{\link[targets:tar_target]{tar_target()}}
except the object must already be an expression instead of
informally quoted code.
\code{base::expression()} and \code{base::quote()} can produce such objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A new target object to combine the return values
from the upstream targets.
See the "Target objects" section for background.
}
\description{
Like \code{\link[=tar_combine]{tar_combine()}} except the \code{name}, \code{command},
and \code{pattern} arguments use standard evaluation.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  target1 <- targets::tar_target(x, head(mtcars))
  target2 <- targets::tar_target(y, tail(mtcars))
  target3 <- tarchetypes::tar_combine(new_target_name, target1, target2)
  list(target1, target2, target3)
})
targets::tar_manifest()
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_by.R
\name{tar_group_by_run}
\alias{tar_group_by_run}
\title{Generate a grouped data frame within tar_group_by()}
\usage{
tar_group_by_run(data, by)
}
\arguments{
\item{data}{A data frame to group.}

\item{by}{Nonempty character vector of names of variables to group by.}
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_render_rep.R
\name{tar_render_rep}
\alias{tar_render_rep}
\title{Parameterized R Markdown with dynamic branching.}
\usage{
tar_render_rep(
  name,
  path,
  params = data.frame(),
  batches = NULL,
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue"),
  quiet = TRUE,
  ...
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{path}{Character string, file path to the R Markdown source file.
Must have length 1.}

\item{params}{Code to generate a data frame or \code{tibble}
with one row per rendered report
and one column per R Markdown parameter. You may also include an
\code{output_file} column to specify the path of each rendered report.}

\item{batches}{Number of batches to group the R Markdown files.
For a large number of reports, increase the number of batches
to decrease target-level overhead. Defaults to the number of
reports to render (1 report per batch).}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}

\item{quiet}{An option to suppress printing during rendering from knitr,
pandoc command line and others. To only suppress printing of the last
"Output created: " message, you can set \code{rmarkdown.render.message} to
\code{FALSE}}

\item{...}{Other named arguments to \code{rmarkdown::render()}.
Unlike \code{\link[=tar_render]{tar_render()}}, these arguments are evaluated when the target
is defined, not when it is run. (The only reason to delay evaluation
in \code{\link[=tar_render]{tar_render()}} was to handle R Markdown parameters, and
\code{tar_render_rep()} handles them differently.)}
}
\value{
A list of target objects to render the R Markdown
reports. Changes to the parameters, source file, dependencies, etc.
will cause the appropriate targets to rerun during \code{tar_make()}.
See the "Target objects" section for background.
}
\description{
Targets to render a parameterized R Markdown report
with multiple sets of parameters.
}
\details{
\code{tar_render_rep()} is an alternative to \code{tar_target()} for
parameterized R Markdown reports that depend on other targets.
Parameters must be given as a data frame with one row per
rendered report and one column per parameter. An optional
\code{output_file} column may be included to set the output file path
of each rendered report.
The R Markdown source should mention other dependency targets
\code{tar_load()} and \code{tar_read()} in the active code chunks
(which also allows you to render the report
outside the pipeline if the \verb{_targets/} data store already exists
and appropriate defaults are specified for the parameters).
(Do not use \code{tar_load_raw()} or \code{tar_read_raw()} for this.)
Then, \code{tar_render()} defines a special kind of target. It
1. Finds all the \code{tar_load()}/\code{tar_read()} dependencies in the report
and inserts them into the target's command.
This enforces the proper dependency relationships.
(Do not use \code{tar_load_raw()} or \code{tar_read_raw()} for this.)
2. Sets \code{format = "file"} (see \code{tar_target()}) so \code{targets}
watches the files at the returned paths and reruns the report
if those files change.
3. Configures the target's command to return the output
report files: the rendered document, the source file,
and then the \verb{*_files/} directory if it exists. All these file paths
are relative paths so the project stays portable.
4. Forces the report to run in the user's current working directory
instead of the working directory of the report.
5. Sets convenient default options such as \code{deployment = "main"}
in the target and \code{quiet = TRUE} in \code{rmarkdown::render()}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
# Parameterized R Markdown:
lines <- c(
  "---",
  "title: 'report.Rmd file'",
  "output_format: html_document",
  "params:",
  "  par: \"default value\"",
  "---",
  "Assume these lines are in a file called report.Rmd.",
  "```{r}",
  "print(params$par)",
  "```"
)
# The following pipeline will run the report for each row of params.
targets::tar_script({
  library(tarchetypes)
  list(
    tar_render_rep(
      report,
      "report.Rmd",
      params = tibble::tibble(par = c(1, 2))
    )
  )
}, ask = FALSE)
# Then, run the targets pipeline as usual.
})
}
}
\seealso{
Other Literate programming targets: 
\code{\link{tar_knit_raw}()},
\code{\link{tar_knit}()},
\code{\link{tar_render_raw}()},
\code{\link{tar_render_rep_raw}()},
\code{\link{tar_render}()}
}
\concept{Literate programming targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_rep_raw.R
\name{tar_rep_raw}
\alias{tar_rep_raw}
\title{Batched replication with dynamic branching
(raw version).}
\usage{
tar_rep_raw(
  name,
  command,
  batches = 1,
  reps = 1,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Character of length 1, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{Expression object with code to run multiple times.
Must return a list or data frame when evaluated.}

\item{batches}{Number of batches. This is also the number of dynamic
branches created during \code{tar_make()}.}

\item{reps}{Number of replications in each batch. The total number
of replications is \code{batches * reps}.}

\item{tidy_eval}{Whether to invoke tidy evaluation
(e.g. the \verb{!!} operator from \code{rlang}) as soon as the target is defined
(before \code{tar_make()}). Applies to the \code{command} argument.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of two target objects, one upstream and one downstream.
The upstream one does some work and returns some file paths,
and the downstream target is a pattern that applies \code{format = "file"}.
See the "Target objects" section for background.
}
\description{
Batching is important for optimizing the efficiency
of heavily dynamically-branched workflows:
\url{https://books.ropensci.org/targets/dynamic.html#batching}.
\code{\link[=tar_rep_raw]{tar_rep_raw()}} is just like \code{\link[=tar_rep]{tar_rep()}} except the
name is a character string and the command is a
language object.
}
\details{
\code{tar_rep_raw()} creates two targets:
an upstream local stem
with an integer vector of batch ids, and a downstream pattern
that maps over the batch ids. (Thus, each batch is a branch.)
Each batch/branch replicates the command a certain number of times.

Both batches and reps within each batch
are aggregated according to the method you specify
in the \code{iteration} argument. If \code{"list"}, reps and batches
are aggregated with \code{list()}. If \code{"vector"},
then \code{vctrs::vec_c()}. If \code{"group"}, then \code{vctrs::vec_rbind()}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_rep_raw(
      "x",
      expression(data.frame(x = sample.int(1e4, 2))),
      batches = 2,
      reps = 3
    )
  )
})
targets::tar_make(callr_function = NULL)
targets::tar_read(x)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map2_size.R
\name{tar_map2_size}
\alias{tar_map2_size}
\title{Dynamic-within-static branching for data frames
(size batching).}
\usage{
tar_map2_size(
  name,
  command1,
  command2,
  values = NULL,
  names = NULL,
  size = Inf,
  combine = TRUE,
  suffix1 = "1",
  suffix2 = "2",
  columns1 = tidyselect::everything(),
  columns2 = tidyselect::everything(),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name of the targets.}

\item{command1}{R code to create named arguments to \code{command2}.
Must return a data frame with one row per call to \code{command2}.}

\item{command2}{R code to map over the data frame of arguments
produced by \code{command1}. Must return a data frame.}

\item{values}{Named list or data frame with values to iterate over.
The names are the names of symbols in the commands and pattern
statements, and the elements are values that get substituted
in place of those symbols. Elements of the \code{values} list
should be small objects that can easily deparse to names,
such as characters, integers, and symbols.
For more complicated elements of \code{values}, such as
lists with multiple numeric vectors,
\code{tar_map()} attempts to parse the elements into expressions,
but this process is not perfect, and the default
target names come out garbled.
To create a list of symbols as a column of \code{values},
use \code{rlang::syms()}.}

\item{names}{Language object with a tidyselect expression
to select which columns of \code{values} to use to construct
statically branched target names. If \code{NULL}, then
short names are automatically generated.}

\item{size}{Positive integer of length 1,
maximum number of rows in each batch for
the downstream (\code{command2}) targets. Batches
are formed from row groups of the \code{command1} target output.}

\item{combine}{Logical of length 1, whether to statically combine
all the results into a single target downstream.}

\item{suffix1}{Character of length 1,
suffix to apply to the \code{command1} targets to distinguish
them from the \code{command2} targets.}

\item{suffix2}{Character of length 1,
suffix to apply to the \code{command2} targets to distinguish
them from the \code{command1} targets.}

\item{columns1}{A tidyselect expression to select which columns of \code{values}
to append to the output of all targets.
Columns already in the target output are not appended.}

\item{columns2}{A tidyselect expression to select which columns of \code{command1}
output to append to \code{command2} output.
Columns already in the target output are not appended.
\code{columns1} takes precedence over \code{columns2}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the output.
An efficient data frame format like \code{"feather"} is recommended,
but the default is \code{"rds"} to avoid incurring extra package
dependencies. See the help file of \code{targets::tar_target()}
for details on storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of new target objects.
See the "Target objects" section for background.
}
\description{
Define targets for batched
dynamic-within-static branching for data frames,
where the user sets the (maximum) size of each batch.
}
\details{
Static branching creates one pair of targets
for each row in \code{values}. In each pair,
there is an upstream non-dynamic target that runs \code{command1}
and a downstream dynamic target that runs \code{command2}.
\code{command1} produces a data frame of arguments to
\code{command2}, and \code{command2} dynamically maps over
these arguments in batches.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  tarchetypes::tar_map2_size(
    x,
    command1 = tibble::tibble(
      arg1 = arg1,
      arg2 = seq_len(6)
     ),
    command2 = tibble::tibble(
      result = paste(arg1, arg2),
      random = sample.int(1e9, size = 1),
      length_input = length(arg1)
    ),
    values = tibble::tibble(arg1 = letters[seq_len(2)]),
    size = 2
   )
})
targets::tar_make()
targets::tar_read(x)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_select.R
\name{tar_group_select_run}
\alias{tar_group_select_run}
\title{Generate a grouped data frame within tar_group_select()}
\usage{
tar_group_select_run(data, by)
}
\arguments{
\item{data}{A data frame to group.}

\item{by}{Nonempty character vector of names of variables to group by.}
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_plan.R
\name{tar_plan}
\alias{tar_plan}
\title{A \code{drake}-plan-like pipeline archetype}
\usage{
tar_plan(...)
}
\arguments{
\item{...}{Named and unnamed targets. All named targets must follow
the \code{drake}-plan-like \code{target = command} syntax, and all unnamed
arguments must be explicit calls to create target objects,
e.g. \code{tar_target()}, target archetypes like \code{\link[=tar_render]{tar_render()}}, or similar.}
}
\value{
A list of \code{tar_target()} objects.
See the "Target objects" section for background.
}
\description{
Simplify target specification in pipelines.
}
\details{
Allows targets with just targets and commands
to be written in the pipeline as \code{target = command} instead of
\code{tar_target(target, command)}. Also supports ordinary
target objects if they are unnamed.
\code{tar_plan(x = 1, y = 2, tar_target(z, 3), tar_render(r, "r.Rmd"))}
is equivalent to
\code{list(tar_target(x, 1), tar_target(y, 2), tar_target(z, 3), tar_render(r, "r.Rmd"))}. # nolint
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  library(tarchetypes)
  tar_plan(
    tarchetypes::tar_fst_tbl(data, data.frame(x = seq_len(26))),
    means = colMeans(data) # No need for tar_target() for simple cases.
  )
})
targets::tar_make()
})
}
}
\concept{Pipeline factories}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_knit.R
\name{tar_knit_run}
\alias{tar_knit_run}
\title{Run a \code{knitr} report inside a \code{tar_knit()} target.}
\usage{
tar_knit_run(path, args, deps)
}
\arguments{
\item{path}{Path to the \code{knitr} source file.}

\item{args}{A named list of arguments to \code{knitr::knit()}.}

\item{deps}{An unnamed list of target dependencies of the \code{knitr}
report, automatically created by \code{tar_knit()}.}
}
\value{
Character with the path to the \code{knitr} source file
and the relative path to the output \code{knitr} report.
The output path depends on the input path argument,
which has no default.
}
\description{
Internal function needed for \code{tar_knit()}.
Users should not invoke it directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_knitr.R
\name{walk_call_knitr}
\alias{walk_call_knitr}
\title{Code analysis for knitr reports.}
\usage{
walk_call_knitr(expr, counter)
}
\arguments{
\item{expr}{A language object or function to scan.}

\item{counter}{An internal counter object that keeps track of detected
target names so far.}
}
\value{
A character vector of target names found during static code analysis.
}
\description{
Walk an abstract syntax tree and capture knitr dependencies.
}
\details{
For internal use only. Not a user-side function.
Powers  automatic detection of \code{tar_load()}/\code{tar_read()}
dependencies in \code{\link[=tar_render]{tar_render()}}.
Packages \code{codetools} and \code{CodeDepends} have different (more sophisticated
and elaborate) implementations of the concepts documented at
\url{https://adv-r.hadley.nz/expressions.html#ast-funs}.
}
\examples{
# How tar_render() really works:
expr <- quote({
  if (a > 1) {
    tar_load(target_name)
  }
  process_stuff(target_name)
})
walk_ast(expr, walk_call_knitr)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_render_rep_raw.R
\name{tar_render_rep_run_params}
\alias{tar_render_rep_run_params}
\title{Prepare R Markdown parameters for \code{tar_render_rep()}.}
\usage{
tar_render_rep_run_params(params, batches)
}
\arguments{
\item{params}{Data frame of R Markdown parameters.}

\item{batches}{Number of batches to split up the renderings.}
}
\value{
A batched data frame of R Markdown parameters.
}
\description{
Internal function needed for \code{tar_render_rep()}.
Users should not invoke it directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map2_size_raw.R
\name{tar_map2_size_raw}
\alias{tar_map2_size_raw}
\title{Dynamic-within-static branching for data frames
(size batching; raw version).}
\usage{
tar_map2_size_raw(
  name,
  command1,
  command2,
  values = NULL,
  names = NULL,
  size = Inf,
  combine = TRUE,
  suffix1 = "1",
  suffix2 = "2",
  columns1 = quote(tidyselect::everything()),
  columns2 = quote(tidyselect::everything()),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Character of length 1, base name of the targets.}

\item{command1}{Language object to create named arguments to \code{command2}.
Must return a data frame with one row per call to \code{command2}.}

\item{command2}{Language object  to map over the data frame of arguments
produced by \code{command1}. Must return a data frame.}

\item{values}{Named list or data frame with values to iterate over.
The names are the names of symbols in the commands and pattern
statements, and the elements are values that get substituted
in place of those symbols. Elements of the \code{values} list
should be small objects that can easily deparse to names,
such as characters, integers, and symbols.
For more complicated elements of \code{values}, such as
lists with multiple numeric vectors,
\code{tar_map()} attempts to parse the elements into expressions,
but this process is not perfect, and the default
target names come out garbled.
To create a list of symbols as a column of \code{values},
use \code{rlang::syms()}.}

\item{names}{Language object with a tidyselect expression
to select which columns of \code{values} to use to construct
statically branched target names. If \code{NULL}, then
short names are automatically generated.}

\item{size}{Positive integer of length 1,
maximum number of rows in each batch for
the downstream (\code{command2}) targets. Batches
are formed from row groups of the \code{command1} target output.}

\item{combine}{Logical of length 1, whether to statically combine
all the results into a single target downstream.}

\item{suffix1}{Character of length 1,
suffix to apply to the \code{command1} targets to distinguish
them from the \code{command2} targets.}

\item{suffix2}{Character of length 1,
suffix to apply to the \code{command2} targets to distinguish
them from the \code{command1} targets.}

\item{columns1}{Language object, a tidyselect expression
to select which columns of \code{values}
to append to the output of all targets.}

\item{columns2}{Language object, a tidyselect expression
to select which columns of \code{command1}
output to append to \code{command2} output.
In case of conflicts, \code{column1} takes precedence.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the output.
An efficient data frame format like \code{"feather"} is recommended,
but the default is \code{"rds"} to avoid incurring extra package
dependencies. See the help file of \code{targets::tar_target()}
for details on storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of new target objects.
See the "Target objects" section for background.
}
\description{
Define targets for batched
dynamic-within-static branching for data frames,
where the user sets the (maximum) size of each batch.
Like \code{tar_map2_size()} except \code{name} is a character string
and \code{command1}, \code{command2}, \code{names}, \code{columns1}, and \code{columns2}
are all language objects.
}
\details{
Static branching creates one pair of targets
for each row in \code{values}. In each pair,
there is an upstream non-dynamic target that runs \code{command1}
and a downstream dynamic target that runs \code{command2}.
\code{command1} produces a data frame of arguments to
\code{command2}, and \code{command2} dynamically maps over
these arguments in batches.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  tarchetypes::tar_map2_size_raw(
    "x",
    command1 = quote(
      tibble::tibble(
        arg1 = arg1,
        arg2 = seq_len(6)
       )
    ),
    command2 = quote(
      tibble::tibble(
        result = paste(arg1, arg2),
        random = sample.int(1e6, size = 1),
        length_input = length(arg1)
      )
    ),
    values = tibble::tibble(arg1 = letters[seq_len(2)]),
    size = 2
   )
})
targets::tar_make()
targets::tar_read(x)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_tidyselect.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{all_of}
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
% Please edit documentation in R/tar_rep2_raw.R
\name{tar_rep2_run}
\alias{tar_rep2_run}
\title{Run \code{\link[=tar_rep2]{tar_rep2()}} batches.}
\usage{
tar_rep2_run(command, batches, iteration)
}
\arguments{
\item{command}{R expression, the command to run on each rep.}

\item{batches}{Named list of batch data to map over.}

\item{iteration}{Iteration method: \code{"list"}, \code{"vector"}, or \code{"group"}.}
}
\value{
The result of batched replication.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_size.R
\name{tar_group_size_index}
\alias{tar_group_size_index}
\title{Generate the tar_group column for \code{tar_group_size()}.}
\usage{
tar_group_size_index(data, size)
}
\arguments{
\item{data}{A data frame to group.}

\item{size}{Maximum number of rows in each group.}
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map2.R
\name{tar_map2}
\alias{tar_map2}
\title{Batched dynamic-within-static branching for data frames.}
\usage{
tar_map2(
  name,
  command1,
  command2,
  values = NULL,
  names = NULL,
  group = rep(1L, nrow(as.data.frame(!!.x))),
  combine = TRUE,
  suffix1 = "1",
  suffix2 = "2",
  columns1 = tidyselect::everything(),
  columns2 = tidyselect::everything(),
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, base name of the targets.}

\item{command1}{R code to create named arguments to \code{command2}.
Must return a data frame with one row per call to \code{command2}.}

\item{command2}{R code to map over the data frame of arguments
produced by \code{command1}. Must return a data frame.}

\item{values}{Named list or data frame with values to iterate over.
The names are the names of symbols in the commands and pattern
statements, and the elements are values that get substituted
in place of those symbols. Elements of the \code{values} list
should be small objects that can easily deparse to names,
such as characters, integers, and symbols.
For more complicated elements of \code{values}, such as
lists with multiple numeric vectors,
\code{tar_map()} attempts to parse the elements into expressions,
but this process is not perfect, and the default
target names come out garbled.
To create a list of symbols as a column of \code{values},
use \code{rlang::syms()}.}

\item{names}{Language object with a tidyselect expression
to select which columns of \code{values} to use to construct
statically branched target names. If \code{NULL}, then
short names are automatically generated.}

\item{group}{Function on the data produced by \code{command1} to create the
\code{tar_group} column that determines the batching structure for the
\code{command2} targets.}

\item{combine}{Logical of length 1, whether to statically combine
all the results into a single target downstream.}

\item{suffix1}{Character of length 1,
suffix to apply to the \code{command1} targets to distinguish
them from the \code{command2} targets.}

\item{suffix2}{Character of length 1,
suffix to apply to the \code{command2} targets to distinguish
them from the \code{command1} targets.}

\item{columns1}{A tidyselect expression to select which columns of \code{values}
to append to the output of all targets.
Columns already in the target output are not appended.}

\item{columns2}{A tidyselect expression to select which columns of \code{command1}
output to append to \code{command2} output.
Columns already in the target output are not appended.
\code{columns1} takes precedence over \code{columns2}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the output.
An efficient data frame format like \code{"feather"} is recommended,
but the default is \code{"rds"} to avoid incurring extra package
dependencies. See the help file of \code{targets::tar_target()}
for details on storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of new target objects.
See the "Target objects" section for background.
}
\description{
Define targets for batched
dynamic-within-static branching for data frames.
Not a user-side function. Do not invoke directly.
}
\details{
Static branching creates one pair of targets
for each row in \code{values}. In each pair,
there is an upstream non-dynamic target that runs \code{command1}
and a downstream dynamic target that runs \code{command2}.
\code{command1} produces a data frame of arguments to
\code{command2}, and \code{command2} dynamically maps over
these arguments in batches.

Static branching creates one pair of targets
for each row in \code{values}. In each pair,
there is an upstream non-dynamic target that runs \code{command1}
and a downstream dynamic target that runs \code{command2}.
\code{command1} produces a data frame of arguments to
\code{command2}, and \code{command2} dynamically maps over
these arguments in batches.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map2_raw.R
\name{tar_map2_raw}
\alias{tar_map2_raw}
\title{Batched dynamic-within-static branching for data frames (raw version).}
\usage{
tar_map2_raw(
  name,
  command1,
  command2,
  values = NULL,
  names = NULL,
  group = quote(rep(1L, nrow(as.data.frame(!!.x)))),
  combine = TRUE,
  columns1 = quote(tidyselect::everything()),
  columns2 = quote(tidyselect::everything()),
  suffix1 = "1",
  suffix2 = "2",
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Character of length 1, base name of the targets.}

\item{command1}{Language object to create named arguments to \code{command2}.
Must return a data frame with one row per call to \code{command2}.}

\item{command2}{Language object  to map over the data frame of arguments
produced by \code{command1}. Must return a data frame.}

\item{values}{Named list or data frame with values to iterate over.
The names are the names of symbols in the commands and pattern
statements, and the elements are values that get substituted
in place of those symbols. Elements of the \code{values} list
should be small objects that can easily deparse to names,
such as characters, integers, and symbols.
For more complicated elements of \code{values}, such as
lists with multiple numeric vectors,
\code{tar_map()} attempts to parse the elements into expressions,
but this process is not perfect, and the default
target names come out garbled.
To create a list of symbols as a column of \code{values},
use \code{rlang::syms()}.}

\item{names}{Language object with a tidyselect expression
to select which columns of \code{values} to use to construct
statically branched target names. If \code{NULL}, then
short names are automatically generated.}

\item{group}{Function on the data produced by \code{command1} to create the
\code{tar_group} column that determines the batching structure for the
\code{command2} targets.}

\item{combine}{Logical of length 1, whether to statically combine
all the results into a single target downstream.}

\item{columns1}{Language object, a tidyselect expression
to select which columns of \code{values}
to append to the output of all targets.}

\item{columns2}{Language object, a tidyselect expression
to select which columns of \code{command1}
output to append to \code{command2} output.
In case of conflicts, \code{column1} takes precedence.}

\item{suffix1}{Character of length 1,
suffix to apply to the \code{command1} targets to distinguish
them from the \code{command2} targets.}

\item{suffix2}{Character of length 1,
suffix to apply to the \code{command2} targets to distinguish
them from the \code{command1} targets.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the output.
An efficient data frame format like \code{"feather"} is recommended,
but the default is \code{"rds"} to avoid incurring extra package
dependencies. See the help file of \code{targets::tar_target()}
for details on storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of new target objects.
See the "Target objects" section for background.
}
\description{
Define targets for batched
dynamic-within-static branching for data frames (raw version).
Not a user-side function. Do not invoke directly.
}
\details{
Like \code{tar_map2()} except \code{name} is a character string
and \code{command1}, \code{command2}, \code{group}, \code{names}, \code{columns1}, and \code{columns2}
are language objects.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_files.R
\name{tar_files}
\alias{tar_files}
\title{Dynamic branching over output or input files.}
\usage{
tar_files(
  name,
  command,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = c("file", "url", "aws_file"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1.
Must be \code{"file"}, \code{"url"}, or \code{"aws_file"}. See the \code{format}
argument of \code{targets::tar_target()} for details.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()}
to customize the rules that decide whether the target is up to date.
Only applies to the downstream target. The upstream target always runs.}
}
\value{
A list of two targets, one upstream and one downstream.
The upstream one does some work and returns some file paths,
and the downstream target is a pattern that applies \code{format = "file"}
or \code{format = "url"}.
See the "Target objects" section for background.
}
\description{
Dynamic branching over output or input files.
}
\details{
\code{tar_files()} creates a pair of targets, one upstream
and one downstream. The upstream target does some work
and returns some file paths, and the downstream
target is a pattern that applies \code{format = "file"}
or \code{format = "url"}. (URLs are input-only, they must already
exist beforehand.)
This is the correct way to dynamically
iterate over file/url targets. It makes sure any downstream patterns
only rerun some of their branches if the files/urls change.
For more information, visit
\url{https://github.com/ropensci/targets/issues/136} and
\url{https://github.com/ropensci/drake/issues/1302}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  # Do not use temp files in real projects
  # or else your targets will always rerun.
  paths <- unlist(replicate(2, tempfile()))
  file.create(paths)
  list(
    tarchetypes::tar_files(x, paths)
  )
})
targets::tar_make()
targets::tar_read(x)
})
}
}
\seealso{
Other Dynamic branching over files: 
\code{\link{tar_files_input_raw}()},
\code{\link{tar_files_input}()},
\code{\link{tar_files_raw}()}
}
\concept{Dynamic branching over files}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_hook_outer.R
\name{tar_hook_outer}
\alias{tar_hook_outer}
\title{Hook to wrap commands}
\usage{
tar_hook_outer(targets, hook, names = NULL)
}
\arguments{
\item{targets}{A list of target objects. The input target list
can be arbitrarily nested, but it must consist entirely of target
objects. In addition, the return value is a simple list
where each element is a target object.
All hook functions remove the nested structure of the input target list.}

\item{hook}{R code to wrap each target's command.
The hook must contain the special placeholder symbol \code{.x}
so \code{tar_hook_outer()} knows where to insert the original command
of the target.
The hook code is quoted (not evaluated) so there is no need
to wrap it in \code{quote()}, \code{expression()}, or similar.}

\item{names}{Name of targets in the target list
to apply the hook. You can supply symbols, a character vector,
or tidyselect helpers like \code{\link[=starts_with]{starts_with()}}.
Targets not included in \code{names} still remain in the target list,
but they are not modified because the hook does not apply to them.}
}
\value{
A flattened list of target objects with the hooks applied.
Even if the input target list had a nested structure,
the return value is a simple list where each element is a target object.
All hook functions remove the nested structure of the input target list.
}
\description{
Wrap the command of each target in an arbitrary R expression.
}
\details{
The expression you supply to \code{hook}
must contain the special placeholder symbol \code{.x}
so \code{tar_hook_outer()} knows where to insert the original command
of the target.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  targets <- list(
    # Nested target lists work with hooks.
    list(
      targets::tar_target(x1, task1()),
      targets::tar_target(x2, task2(x1))
    ),
    targets::tar_target(x3, task3(x2)),
    targets::tar_target(y1, task4(x3))
  )
  tarchetypes::tar_hook_outer(
    targets = targets,
    hook = postprocess(.x, arg = "value"),
    names = starts_with("x")
  )
})
targets::tar_manifest(fields = command)
})
}
}
\seealso{
Other hooks: 
\code{\link{tar_hook_before}()},
\code{\link{tar_hook_inner}()}
}
\concept{hooks}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_sub_raw.R
\name{tar_sub_raw}
\alias{tar_sub_raw}
\title{Create multiple expressions with symbol substitution (raw version).}
\usage{
tar_sub_raw(expr, values)
}
\arguments{
\item{expr}{Expression object with the starting expression.
Values are iteratively substituted
in place of symbols in \code{expr} to create each new expression.}

\item{values}{List of values to substitute into \code{expr} to create
the expressions. All elements of \code{values} must have the same length.}
}
\value{
A list of expression objects. Often, these expression objects
evaluate to target objects (but not necessarily).
See the "Target objects" section for background.
}
\description{
Loop over a grid of values and create an expression object
from each one. Helps with general metaprogramming. Unlike \code{\link[=tar_sub]{tar_sub()}},
which quotes the \code{expr} argument, \code{tar_sub_raw()} assumes \code{expr}
is an expression object.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
# tar_map() is incompatible with tar_render() because the latter
# operates on preexisting tar_target() objects. By contrast,
# tar_eval_raw() and tar_sub_raw() iterate over code farther upstream.
values <- list(
  name = lapply(c("name1", "name2"), as.symbol),
  file = c("file1.Rmd", "file2.Rmd")
)
tar_sub_raw(quote(tar_render(name, file)), values = values)
}
\seealso{
Other Metaprogramming utilities: 
\code{\link{tar_eval_raw}()},
\code{\link{tar_eval}()},
\code{\link{tar_sub}()}
}
\concept{Metaprogramming utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_package.R
\docType{package}
\name{tarchetypes-package}
\alias{tarchetypes-package}
\title{targets: Archetypes for Targets}
\description{
A pipeline toolkit for R, the \code{targets} package brings together
function-oriented programming and Make-like declarative pipelines for
Statistics and data science. The \code{tarchetypes} package provides
convenient helper functions to create specialized targets, making
pipelines in targets easier and cleaner to write and understand.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_render_rep_raw.R
\name{tar_render_rep_raw}
\alias{tar_render_rep_raw}
\title{Parameterized R Markdown with dynamic branching (raw version).}
\usage{
tar_render_rep_raw(
  name,
  path,
  params = expression(NULL),
  batches = NULL,
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue"),
  quiet = TRUE,
  args = list()
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{path}{Character string, file path to the R Markdown source file.
Must have length 1.}

\item{params}{Expression object with code to generate
a data frame or \code{tibble} with one row per rendered report
and one column per R Markdown parameter. You may also include an
\code{output_file} column to specify the path of each rendered report.
R Markdown parameters must not be named \code{tar_group} or \code{output_file}.}

\item{batches}{Number of batches to group the R Markdown files.
For a large number of reports, increase the number of batches
to decrease target-level overhead. Defaults to the number of
reports to render (1 report per batch).}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, \code{format} argument to \code{tar_target()}
to store the data frame of R Markdown parameters.}

\item{iteration}{Character of length 1, \code{iteration} argument
to \code{tar_target()} for the R Markdown documents. Does not apply
to the target with R Markdown parameters (whose iteration
is always \code{"group"}).}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}

\item{quiet}{An option to suppress printing during rendering from knitr,
pandoc command line and others. To only suppress printing of the last
"Output created: " message, you can set \code{rmarkdown.render.message} to
\code{FALSE}}

\item{args}{Named list of other arguments to \code{rmarkdown::render()}.
Must not include \code{params} or \code{output_file}. Evaluated when the target
is defined.}
}
\value{
A list of target objects to render the R Markdown
reports. Changes to the parameters, source file, dependencies, etc.
will cause the appropriate targets to rerun during \code{tar_make()}.
See the "Target objects" section for background.
}
\description{
Targets to render a parameterized R Markdown report
with multiple sets of parameters (raw version). Same as
\code{tar_render_rep()} except \code{name} is a character string,
\code{params} is an expression object,
and extra arguments to \code{rmarkdown::render()} are passed through
the \code{args} argument instead of \code{...}.
}
\details{
\code{tar_render_rep_raw()} is an alternative to \code{tar_target_raw()} for
parameterized R Markdown reports that depend on other targets.
Parameters must be given as a data frame with one row per
rendered report and one column per parameter. An optional
\code{output_file} column may be included to set the output file path
of each rendered report.
The R Markdown source should mention other dependency targets
\code{tar_load()} and \code{tar_read()} in the active code chunks
(which also allows you to render the report
outside the pipeline if the \verb{_targets/} data store already exists
and appropriate defaults are specified for the parameters).
(Do not use \code{tar_load_raw()} or \code{tar_read_raw()} for this.)
Then, \code{tar_render()} defines a special kind of target. It
1. Finds all the \code{tar_load()}/\code{tar_read()} dependencies in the report
and inserts them into the target's command.
This enforces the proper dependency relationships.
(Do not use \code{tar_load_raw()} or \code{tar_read_raw()} for this.)
2. Sets \code{format = "file"} (see \code{tar_target()}) so \code{targets}
watches the files at the returned paths and reruns the report
if those files change.
3. Configures the target's command to return the output
report files: the rendered document, the source file,
and then the \verb{*_files/} directory if it exists.
All these file paths
are relative paths so the project stays portable.
4. Forces the report to run in the user's current working directory
instead of the working directory of the report.
5. Sets convenient default options such as \code{deployment = "main"}
in the target and \code{quiet = TRUE} in \code{rmarkdown::render()}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
# Parameterized R Markdown:
lines <- c(
  "---",
  "title: 'report.Rmd source file'",
  "output_format: html_document",
  "params:",
  "  par: \"default value\"",
  "---",
  "Assume these lines are in a file called report.Rmd.",
  "```{r}",
  "print(params$par)",
  "```"
)
# The following pipeline will run the report for each row of params.
targets::tar_script({
  library(tarchetypes)
  list(
    tar_render_rep_raw(
      "report",
      "report.Rmd",
      params = quote(tibble::tibble(par = c(1, 2)))
    )
  )
}, ask = FALSE)
# Then, run the targets pipeline as usual.
})
}
}
\seealso{
Other Literate programming targets: 
\code{\link{tar_knit_raw}()},
\code{\link{tar_knit}()},
\code{\link{tar_render_raw}()},
\code{\link{tar_render_rep}()},
\code{\link{tar_render}()}
}
\concept{Literate programming targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_render_raw.R
\name{tar_render_run}
\alias{tar_render_run}
\title{Render an R Markdown report inside a \code{tar_render()} target.}
\usage{
tar_render_run(path, args, deps)
}
\arguments{
\item{path}{Path to the R Markdown source file.}

\item{args}{A named list of arguments to \code{rmarkdown::render()}.}

\item{deps}{An unnamed list of target dependencies of the R Markdown
report, automatically created by \code{tar_render()}.}
}
\value{
Character vector with the path to the R Markdown source file
and the relative path to the output. These paths depend on the input
source file path and have no defaults.
}
\description{
Internal function needed for \code{tar_render()}.
Users should not invoke it directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_render_rep_raw.R
\name{tar_render_rep_run}
\alias{tar_render_rep_run}
\title{Render a batch of parameterized R Markdown reports
inside a \code{tar_render_rep()} target.}
\usage{
tar_render_rep_run(path, params, args, deps)
}
\arguments{
\item{path}{Path to the R Markdown source file.}

\item{args}{A named list of arguments to \code{rmarkdown::render()}.}

\item{deps}{An unnamed list of target dependencies of the R Markdown
report, automatically created by \code{tar_render_rep()}.}
}
\value{
Character vector with the path to the R Markdown
source file and the rendered output file. Both paths
depend on the input source path, and they have no defaults.
}
\description{
Internal function needed for \code{tar_render()}.
Users should not invoke it directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_formats.R
\name{tar_formats}
\alias{tar_formats}
\alias{tar_url}
\alias{tar_file}
\alias{tar_rds}
\alias{tar_qs}
\alias{tar_keras}
\alias{tar_torch}
\alias{tar_format_feather}
\alias{tar_parquet}
\alias{tar_fst}
\alias{tar_fst_dt}
\alias{tar_fst_tbl}
\alias{tar_aws_file}
\alias{tar_aws_rds}
\alias{tar_aws_qs}
\alias{tar_aws_keras}
\alias{tar_aws_torch}
\alias{tar_format_aws_feather}
\alias{tar_aws_parquet}
\alias{tar_aws_fst}
\alias{tar_aws_fst_dt}
\alias{tar_aws_fst_tbl}
\title{Target formats}
\usage{
tar_url(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_file(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_rds(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_qs(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_keras(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_torch(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_format_feather(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_parquet(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_fst(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_fst_dt(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_fst_tbl(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_aws_file(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_aws_rds(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_aws_qs(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_aws_keras(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_aws_torch(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_format_aws_feather(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_aws_parquet(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_aws_fst(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_aws_fst_dt(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)

tar_aws_fst_tbl(
  name,
  command,
  pattern = NULL,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{pattern}{Language to define branching for a target.
For example, in a pipeline with numeric vector targets \code{x} and \code{y},
\code{tar_target(z, x + y, pattern = map(x, y))} implicitly defines
branches of \code{z} that each compute \code{x[1] + y[1]}, \code{x[2] + y[2]},
and so on. See the user manual for details.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A \code{tar_target()} object with the eponymous storage format.
See the "Target objects" section for background.
}
\description{
Target archetypes for specialized storage formats.
}
\details{
These functions are shorthand for targets with specialized
storage formats. For example, \code{tar_qs(name, fun())} is equivalent to
\code{tar_target(name, fun(), format = "qs")}.
For details on specialized storage formats, open the help file of the
\code{targets::tar_target()} function and read about the \code{format} argument.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script(
  list(
    tarchetypes::tar_rds(x, 1)
  )
)
targets::tar_make()
})
}
}
\concept{Formats}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group_count.R
\name{tar_group_count_index}
\alias{tar_group_count_index}
\title{Generate the tar_group column for \code{tar_group_count()}.}
\usage{
tar_group_count_index(data, count)
}
\arguments{
\item{data}{A data frame to group.}

\item{count}{Maximum number of groups.}
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_files_raw.R
\name{tar_files_raw}
\alias{tar_files_raw}
\title{Dynamic branching over output or input files (raw version).}
\usage{
tar_files_raw(
  name,
  command,
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = c("file", "url", "aws_file"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1.
Must be \code{"file"}, \code{"url"}, or \code{"aws_file"}. See the \code{format}
argument of \code{targets::tar_target()} for details.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()}
to customize the rules that decide whether the target is up to date.
Only applies to the downstream target. The upstream target always runs.}
}
\value{
A list of two targets, one upstream and one downstream.
The upstream one does some work and returns some file paths,
and the downstream target is a pattern that applies \code{format = "file"}
or \code{format = "url"}.
See the "Target objects" section for background.
}
\description{
Dynamic branching over output or input files.
}
\details{
\code{tar_files_raw()} is similar to \code{\link[=tar_files]{tar_files()}}
except the \code{name} argument must be a character string
and \code{command} must be a language object.

\code{tar_files_raw()} creates a pair of targets, one upstream
and one downstream. The upstream target does some work
and returns some file paths, and the downstream
target is a pattern that applies \code{format = "file"}
or \code{format = "url"}. (URLs are input-only, they must already
exist beforehand.)
This is the correct way to dynamically
iterate over file/url targets. It makes sure any downstream patterns
only rerun some of their branches if the files/urls change.
For more information, visit
\url{https://github.com/ropensci/targets/issues/136} and
\url{https://github.com/ropensci/drake/issues/1302}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  # Do not use temp files in real projects
  # or else your targets will always rerun.
  paths <- unlist(replicate(2, tempfile()))
  file.create(paths)
  command <- as.call(list(`c`, paths))
  list(
    tarchetypes::tar_files_raw("x", command)
  )
})
targets::tar_make()
targets::tar_read(x)
})
}
}
\seealso{
Other Dynamic branching over files: 
\code{\link{tar_files_input_raw}()},
\code{\link{tar_files_input}()},
\code{\link{tar_files}()}
}
\concept{Dynamic branching over files}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_rep2_raw.R
\name{tar_rep2_raw}
\alias{tar_rep2_raw}
\title{Dynamic batched computation downstream of \code{\link[=tar_rep]{tar_rep()}} (raw version).}
\usage{
tar_rep2_raw(
  name,
  command,
  targets,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run the target.}

\item{targets}{Character vector of names of upstream batched targets
created by \code{\link[=tar_rep]{tar_rep()}}.
If you supply more than one such target, all those targets must have the
same number of batches and reps per batch. And they must all return
either data frames or lists. List targets must use \code{iteration = "list"}
in \code{\link[=tar_rep]{tar_rep()}}.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vctrs::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[targets:tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A new target object to perform batched computation
downstream of \code{\link[=tar_rep]{tar_rep()}}.
See the "Target objects" section for background.
}
\description{
Batching is important for optimizing the efficiency
of heavily dynamically-branched workflows:
\url{https://books.ropensci.org/targets/dynamic.html#batching}.
\code{tar_rep2_raw()}
is just like \code{\link[=tar_rep2]{tar_rep2()}} except it accepts a character
of length 1 for \code{name}, a language object for \code{command},
and a character vector of the names of the upstream batched targets.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_rep(
      data1,
      data.frame(value = rnorm(1)),
      batches = 2,
      reps = 3
    ),
    tarchetypes::tar_rep(
      data2,
      list(value = rnorm(1)),
      batches = 2, reps = 3,
      iteration = "list" # List iteration is important for batched lists.
    ),
    tarchetypes::tar_rep2_raw(
      "aggregate",
      quote(data.frame(value = data1$value + data2$value)),
      targets = c("data1", "data2")
    )
  )
})
targets::tar_make()
targets::tar_read(aggregate)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map_rep_raw.R
\name{tar_append_static_values}
\alias{tar_append_static_values}
\title{Append statically mapped values to target output.}
\usage{
tar_append_static_values(object, values)
}
\arguments{
\item{object}{Return value of a target. Must be a data frame.}

\item{values}{Tibble with the set of static values that the current target
uses.}
}
\description{
For internal use only. Users should not invoke
this function directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_rep.R
\name{tar_rep}
\alias{tar_rep}
\title{Batched replication with dynamic branching.}
\usage{
tar_rep(
  name,
  command,
  batches = 1,
  reps = 1,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code to run multiple times. Must return a list or
data frame because \code{tar_rep()} will try to append new elements/columns
\code{tar_batch} and \code{tar_rep} to the output to denote the batch
and rep-within-batch IDs, respectively.}

\item{batches}{Number of batches. This is also the number of dynamic
branches created during \code{tar_make()}.}

\item{reps}{Number of replications in each batch. The total number
of replications is \code{batches * reps}.}

\item{tidy_eval}{Whether to invoke tidy evaluation
(e.g. the \verb{!!} operator from \code{rlang}) as soon as the target is defined
(before \code{tar_make()}). Applies to the \code{command} argument.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Optional storage format for the target's return value.
With the exception of \code{format = "file"}, each target
gets a file in \verb{_targets/objects}, and each format is a different
way to save and load this file. See the "Storage formats" section
for a detailed list of possible data storage formats.}

\item{iteration}{Character of length 1, name of the iteration mode
of the target. Choices:
\itemize{
\item \code{"vector"}: branching happens with \code{vectors::vec_slice()} and
aggregation happens with \code{vctrs::vec_c()}.
\item \code{"list"}, branching happens with \verb{[[]]} and aggregation happens with
\code{list()}. In the case of list iteration, \code{tar_read(your_target)}
will return a list of lists, where the outer list has one element per
batch and each inner list has one element per rep within batch.
To un-batch this nested list, call
\code{tar_read(your_target, recursive = FALSE)}.
\item \code{"group"}: \code{dplyr::group_by()}-like functionality to branch over
subsets of a data frame. The target's return value must be a data
frame with a special \code{tar_group} column of consecutive integers
from 1 through the number of groups. Each integer designates a group,
and a branch is created for each collection of rows in a group.
See the \code{\link[=tar_group]{tar_group()}} function to see how you can
create the special \code{tar_group} column with \code{dplyr::group_by()}.
}}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of two targets, one upstream and one downstream.
The upstream target returns a numeric index of batch ids,
and the downstream one dynamically maps over the batch ids
to run the command multiple times.
If the command returns a list or data frame, then
the targets from \code{tar_rep()} will try to append new elements/columns
\code{tar_batch} and \code{tar_rep} to the output
to denote the batch and rep-within-batch IDs, respectively.
See the "Target objects" section for background.

\code{tar_read(your_target)} (on the downstream target with the actual work)
will return a list of lists, where the outer list has one element per
batch and each inner list has one element per rep within batch.
To un-batch this nested list, call
\code{tar_read(your_target, recursive = FALSE)}.
}
\description{
Batching is important for optimizing the efficiency
of heavily dynamically-branched workflows:
\url{https://books.ropensci.org/targets/dynamic.html#batching}.
\code{\link[=tar_rep]{tar_rep()}} replicates a command in strategically sized batches.
}
\details{
\code{tar_rep()} and \code{tar_rep_raw()} each create two targets:
an upstream local stem
with an integer vector of batch ids, and a downstream pattern
that maps over the batch ids. (Thus, each batch is a branch.)
Each batch/branch replicates the command a certain number of times.
If the command returns a list or data frame, then
the targets from \code{tar_rep()} will try to append new elements/columns
\code{tar_batch} and \code{tar_rep} to the output
to denote the batch and rep-within-batch IDs, respectively.

Both batches and reps within each batch
are aggregated according to the method you specify
in the \code{iteration} argument. If \code{"list"}, reps and batches
are aggregated with \code{list()}. If \code{"vector"},
then \code{vctrs::vec_c()}. If \code{"group"}, then \code{vctrs::vec_rbind()}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  list(
    tarchetypes::tar_rep(
      x,
      data.frame(x = sample.int(1e4, 2)),
      batches = 2,
      reps = 3
    )
  )
})
targets::tar_make()
targets::tar_read(x)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map_rep}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_render.R
\name{tar_render}
\alias{tar_render}
\title{Target with an R Markdown document.}
\usage{
tar_render(
  name,
  path,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  error = targets::tar_option_get("error"),
  deployment = "main",
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue"),
  quiet = TRUE,
  ...
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{path}{Character string, file path to the R Markdown source file.
Must have length 1.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}

\item{quiet}{An option to suppress printing during rendering from knitr,
pandoc command line and others. To only suppress printing of the last
"Output created: " message, you can set \code{rmarkdown.render.message} to
\code{FALSE}}

\item{...}{Named arguments to \code{rmarkdown::render()}.
These arguments are evaluated when the target actually runs in
\code{tar_make()}, not when the target is defined. That means, for
example, you can use upstream targets as parameters of
parameterized R Markdown reports.
\code{tar_render(your_target, "your_report.Rmd", params = list(your_param = your_target))} # nolint
will run \code{rmarkdown::render("your_report.Rmd", params = list(your_param = your_target))}. # nolint
For parameterized reports, it is recommended to supply a distinct
\code{output_file} argument to each \code{tar_render()} call
and set useful defaults for parameters in the R Markdown source.
See the examples section for a demonstration.}
}
\value{
A target object with \code{format = "file"}.
When this target runs, it returns a character vector
of file paths: the rendered document, the source file,
and then the \verb{*_files/} directory if it exists.
Unlike \code{rmarkdown::render()},
all returned paths are \emph{relative} paths to ensure portability
(so that the project can be moved from one file system to another
without invalidating the target).
See the "Target objects" section for background.
}
\description{
Shorthand to include an R Markdown document in a
\code{targets} pipeline.
}
\details{
\code{tar_render()} is an alternative to \code{tar_target()} for
R Markdown reports that depend on other targets. The R Markdown source
should mention dependency targets with \code{tar_load()} and \code{tar_read()}
in the active code chunks (which also allows you to render the report
outside the pipeline if the \verb{_targets/} data store already exists).
(Do not use \code{tar_load_raw()} or \code{tar_read_raw()} for this.)
Then, \code{tar_render()} defines a special kind of target. It
1. Finds all the \code{tar_load()}/\code{tar_read()} dependencies in the report
and inserts them into the target's command.
This enforces the proper dependency relationships.
(Do not use \code{tar_load_raw()} or \code{tar_read_raw()} for this.)
2. Sets \code{format = "file"} (see \code{tar_target()}) so \code{targets}
watches the files at the returned paths and reruns the report
if those files change.
3. Configures the target's command to return both the output
report files and the input source file. All these file paths
are relative paths so the project stays portable.
4. Forces the report to run in the user's current working directory
instead of the working directory of the report.
5. Sets convenient default options such as \code{deployment = "main"}
in the target and \code{quiet = TRUE} in \code{rmarkdown::render()}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({  # tar_dir() runs code from a temporary directory.
# Unparameterized R Markdown:
lines <- c(
  "---",
  "title: report.Rmd source file",
  "output_format: html_document",
  "---",
  "Assume these lines are in report.Rmd.",
  "```{r}",
  "targets::tar_read(data)",
  "```"
)
# Include the report in a pipeline as follows.
targets::tar_script({
  library(tarchetypes)
  list(
    tar_target(data, data.frame(x = seq_len(26), y = letters)),
    tar_render(report, "report.Rmd")
  )
}, ask = FALSE)
# Then, run the targets pipeline as usual.

# Parameterized R Markdown:
lines <- c(
  "---",
  "title: 'report.Rmd source file with parameters'",
  "output_format: html_document",
  "params:",
  "  your_param: \"default value\"",
  "---",
  "Assume these lines are in report.Rmd.",
  "```{r}",
  "print(params$your_param)",
  "```"
)
# Include the report in the pipeline as follows.
targets::tar_script({
  library(tarchetypes)
  list(
    tar_target(data, data.frame(x = seq_len(26), y = letters)),
    tar_render(report, "report.Rmd", params = list(your_param = data))
  )
}, ask = FALSE)
})
# Then, run the targets pipeline as usual.
}
}
\seealso{
Other Literate programming targets: 
\code{\link{tar_knit_raw}()},
\code{\link{tar_knit}()},
\code{\link{tar_render_raw}()},
\code{\link{tar_render_rep_raw}()},
\code{\link{tar_render_rep}()}
}
\concept{Literate programming targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_hook_inner.R
\name{tar_hook_inner}
\alias{tar_hook_inner}
\title{Hook to wrap dependencies}
\usage{
tar_hook_inner(targets, hook, names = NULL, names_wrap = NULL)
}
\arguments{
\item{targets}{A list of target objects. The input target list
can be arbitrarily nested, but it must consist entirely of target
objects. In addition, the return value is a simple list
where each element is a target object.
All hook functions remove the nested structure of the input target list.}

\item{hook}{R code to wrap each target's command.
The hook must contain the special placeholder symbol \code{.x}
so \code{tar_hook_inner()} knows where to insert the code to wrap
mentions of dependencies.
The hook code is quoted (not evaluated) so there is no need
to wrap it in \code{quote()}, \code{expression()}, or similar.}

\item{names}{Name of targets in the target list
to apply the hook. You can supply symbols, a character vector,
or tidyselect helpers like \code{\link[=starts_with]{starts_with()}}.
Targets not included in \code{names} still remain in the target list,
but they are not modified because the hook does not apply to them.}

\item{names_wrap}{Names of targets to wrap with the hook
where they appear as dependencies in the commands of other targets.
You can supply symbols, a character vector,
or tidyselect helpers like \code{\link[=starts_with]{starts_with()}}.}
}
\value{
A flattened list of target objects with the hooks applied.
Even if the input target list had a nested structure,
the return value is a simple list where each element is a target object.
All hook functions remove the nested structure of the input target list.
}
\description{
In the command of each target, wrap each mention of
each dependency target in an arbitrary R expression.
}
\details{
The expression you supply to \code{hook}
must contain the special placeholder symbol \code{.x}
so \code{tar_hook_inner()} knows where to insert the original command
of the target.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  targets <- list(
    # Nested target lists work with hooks.
    list(
      targets::tar_target(x1, task1()),
      targets::tar_target(x2, task2(x1))
    ),
    targets::tar_target(x3, task3(x2, x1)),
    targets::tar_target(y1, task4(x3))
  )
  tarchetypes::tar_hook_inner(
    targets = targets,
    hook = fun(.x),
    names = starts_with("x")
  )
})
targets::tar_manifest(fields = command)
})
}
}
\seealso{
Other hooks: 
\code{\link{tar_hook_before}()},
\code{\link{tar_hook_outer}()}
}
\concept{hooks}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map2_raw.R
\name{tar_map2_run}
\alias{tar_map2_run}
\title{Run a dynamic batch of a \code{tar_map2()} target.}
\usage{
tar_map2_run(command, values, columns)
}
\arguments{
\item{command}{Command to run.}

\item{values}{Data frame of named arguments produced by \code{command1}
that \code{command2} dynamically maps over. Different from the \code{values}
argument of \code{tar_map2()}.}

\item{columns}{tidyselect expression to select columns of \code{values}
to append to the result.}
}
\value{
A data frame with a \code{tar_group} column attached (if \code{group}
is not \code{NULL}).
}
\description{
Run a dynamic batch of a \code{tar_map2()} target.
}
\details{
For internal use only. Users should not invoke
this function directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_map_rep.R
\name{tar_map_rep}
\alias{tar_map_rep}
\title{Dynamic batched replication within static branches
for data frames.}
\usage{
tar_map_rep(
  name,
  command,
  values = NULL,
  names = NULL,
  columns = tidyselect::everything(),
  batches = 1,
  reps = 1,
  combine = TRUE,
  tidy_eval = targets::tar_option_get("tidy_eval"),
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  format = targets::tar_option_get("format"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  deployment = targets::tar_option_get("deployment"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  storage = targets::tar_option_get("storage"),
  retrieval = targets::tar_option_get("retrieval"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{command}{R code for a single replicate. Must return
a data frame.}

\item{values}{Named list or data frame with values to iterate over.
The names are the names of symbols in the commands and pattern
statements, and the elements are values that get substituted
in place of those symbols. Elements of the \code{values} list
should be small objects that can easily deparse to names,
such as characters, integers, and symbols.
For more complicated elements of \code{values}, such as
lists with multiple numeric vectors,
\code{tar_map()} attempts to parse the elements into expressions,
but this process is not perfect, and the default
target names come out garbled.
To create a list of symbols as a column of \code{values},
use \code{rlang::syms()}.}

\item{names}{Language object with a tidyselect expression
to select which columns of \code{values} to use to construct
statically branched target names. If \code{NULL}, then
short names are automatically generated.}

\item{columns}{A tidyselect expression to select which columns of \code{values}
to append to the output. Columns already in the target output
are not appended.}

\item{batches}{Number of batches. This is also the number of dynamic
branches created during \code{tar_make()}.}

\item{reps}{Number of replications in each batch. The total number
of replications is \code{batches * reps}.}

\item{combine}{Logical of length 1, whether to statically combine
all the results into a single target downstream.}

\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{format}{Character of length 1, storage format of the output.
An efficient data frame format like \code{"feather"} is recommended,
but the default is \code{"rds"} to avoid incurring extra package
dependencies. See the help file of \code{targets::tar_target()}
for details on storage formats.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{deployment}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's return value is sent back to the
host machine and saved/uploaded locally.
\item \code{"worker"}: the worker saves/uploads the value.
\item \code{"none"}: almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language. If you do use it,
then the return value of the target is totally ignored
when the target ends, but
each downstream target still attempts to load the data file
(except when \code{retrieval = "none"}).

If you select \code{storage = "none"}, then
the return value of the target's command is ignored,
and the data is not saved automatically.
As with dynamic files (\code{format = "file"} or \code{"aws_file"}) it is the
responsibility of the user to write to
\code{\link[targets:tar_path]{tar_path()}} from inside the target.
An example target
could look something like
tar_target(x,
{saveRDS("value", tar_path(create_dir = TRUE)); "ignored"},
storage = "none")`.

The distinguishing feature of \code{storage = "none"}
(as opposed to \code{format = "file"} or \code{"aws_file"})
is that in the general case,
downstream targets will automatically try to load the data
from the data store as a dependency. As a corollary, \code{storage = "none"}
is completely unnecessary if \code{format} is \code{"file"} or \code{"aws_file"}.
}}

\item{retrieval}{Character of length 1, only relevant to
\code{\link[targets:tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[targets:tar_make_future]{tar_make_future()}}.
Must be one of the following values:
\itemize{
\item \code{"main"}: the target's dependencies are loaded on the host machine
and sent to the worker before the target builds.
\item \code{"worker"}: the worker loads the targets dependencies.
\item \code{"none"}: the dependencies are not loaded at all.
This choice is almost never recommended. It is only for
niche situations, e.g. the data needs to be loaded
explicitly from another language.
}}

\item{cue}{An optional object from \code{tar_cue()} to customize the
rules that decide whether the target is up to date.}
}
\value{
A list of new target objects.
See the "Target objects" section for background.
}
\description{
Define targets for batched replication
within static branches for  data frames.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.


Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  # Just a sketch of a Bayesian sensitivity analysis of hyperparameters:
  assess_hyperparameters <- function(sigma1, sigma2) {
    # data <- simulate_random_data() # user-defined function
    # run_model(data, sigma1, sigma2) # user-defined function
    # Mock output from the model:
    posterior_samples <- stats::rnorm(1000, 0, sigma1 + sigma2)
    tibble::tibble(
      posterior_median = median(posterior_samples),
      posterior_quantile_0.025 = quantile(posterior_samples, 0.025),
      posterior_quantile_0.975 = quantile(posterior_samples, 0.975)
    )
  }
  hyperparameters <- tibble::tibble(
    scenario = c("tight", "medium", "diffuse"),
    sigma1 = c(10, 50, 50),
    sigma2 = c(10, 5, 10)
  )
  tarchetypes::tar_map_rep(
    sensitivity_analysis,
    command = assess_hyperparameters(sigma1, sigma2),
    values = hyperparameters,
    names = tidyselect::any_of("scenario"),
    batches = 2,
    reps = 3
   )
})
targets::tar_make()
targets::tar_read(sensitivity_analysis)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_combine_raw}()},
\code{\link{tar_combine}()},
\code{\link{tar_map2_count_raw}()},
\code{\link{tar_map2_count}()},
\code{\link{tar_map2_raw}()},
\code{\link{tar_map2_size_raw}()},
\code{\link{tar_map2_size}()},
\code{\link{tar_map2}()},
\code{\link{tar_map_rep_raw}()},
\code{\link{tar_map}()},
\code{\link{tar_rep2_raw}()},
\code{\link{tar_rep2}()},
\code{\link{tar_rep_map_raw}()},
\code{\link{tar_rep_map}()},
\code{\link{tar_rep_raw}()},
\code{\link{tar_rep}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_files_input_raw.R
\name{tar_files_input_raw}
\alias{tar_files_input_raw}
\title{Dynamic branching over input files or URLs (raw version).}
\usage{
tar_files_input_raw(
  name,
  files,
  batches = length(files),
  format = c("file", "url", "aws_file"),
  iteration = targets::tar_option_get("iteration"),
  error = targets::tar_option_get("error"),
  memory = targets::tar_option_get("memory"),
  garbage_collection = targets::tar_option_get("garbage_collection"),
  priority = targets::tar_option_get("priority"),
  resources = targets::tar_option_get("resources"),
  cue = targets::tar_option_get("cue")
)
}
\arguments{
\item{name}{Symbol, name of the target. A target
name must be a valid name for a symbol in R, and it
must not start with a dot. Subsequent targets
can refer to this name symbolically to induce a dependency relationship:
e.g. \code{tar_target(downstream_target, f(upstream_target))} is a
target named \code{downstream_target} which depends on a target
\code{upstream_target} and a function \code{f()}. In addition, a target's
name determines its random number generator seed. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()} on the result
to locally recreate the target's initial RNG state.}

\item{files}{Nonempty character vector of known existing input files
to track for changes.}

\item{batches}{Positive integer of length 1, number of batches
to partition the files. The default is one file per batch
(maximum number of batches) which is simplest to handle but
could cause a lot of overhead and consume a lot of computing resources.
Consider reducing the number of batches below the number of files
for heavy workloads.}

\item{format}{Character, either \code{"file"} or \code{"url"}. See the \code{format}
argument of \code{targets::tar_target()} for details.}

\item{iteration}{Character, iteration method. Must be a method
supported by the \code{iteration} argument of \code{targets::tar_target()}.
The iteration method for the upstream target is always \code{"list"}
in order to support batching.}

\item{error}{Character of length 1, what to do if the target
stops and throws an error. Options:
\itemize{
\item \code{"stop"}: the whole pipeline stops and throws an error.
\item \code{"continue"}: the whole pipeline keeps going.
\item \code{"abridge"}: any currently running targets keep running,
but no new targets launch after that.
(Visit \url{https://books.ropensci.org/targets/debugging.html}
to learn how to debug targets using saved workspaces.)
}}

\item{memory}{Character of length 1, memory strategy.
If \code{"persistent"}, the target stays in memory
until the end of the pipeline (unless \code{storage} is \code{"worker"},
in which case \code{targets} unloads the value from memory
right after storing it in order to avoid sending
copious data over a network).
If \code{"transient"}, the target gets unloaded
after every new target completes.
Either way, the target gets automatically loaded into memory
whenever another target needs the value.
For cloud-based dynamic files such as \code{format = "aws_file"},
this memory strategy applies to
temporary local copies of the file in \verb{_targets/scratch/"}:
\code{"persistent"} means they remain until the end of the pipeline,
and \code{"transient"} means they get deleted from the file system
as soon as possible. The former conserves bandwidth,
and the latter conserves local storage.}

\item{garbage_collection}{Logical, whether to run \code{base::gc()}
just before the target runs.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[targets:tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{cue}{An optional object from \code{tar_cue()}
to customize the rules that decide whether the target is up to date.
Only applies to the downstream target. The upstream target always runs.}
}
\value{
A list of two targets, one upstream and one downstream.
The upstream one does some work and returns some file paths,
and the downstream target is a pattern that applies \code{format = "file"}
or \code{format = "url"}.
See the "Target objects" section for background.
}
\description{
Dynamic branching over input files or URLs.
}
\details{
\code{tar_files_input_raw()} is similar to \code{\link[=tar_files_input]{tar_files_input()}}
except the \code{name} argument must be a character string.

\code{tar_files_input_raw()} creates a pair of targets, one upstream
and one downstream. The upstream target does some work
and returns some file paths, and the downstream
target is a pattern that applies \code{format = "file"}
or \code{format = "url"}.
This is the correct way to dynamically
iterate over file/url targets. It makes sure any downstream patterns
only rerun some of their branches if the files/urls change.
For more information, visit
\url{https://github.com/ropensci/targets/issues/136} and
\url{https://github.com/ropensci/drake/issues/1302}.
}
\section{Target objects}{

Most \code{tarchetypes} functions are target factories,
which means they return target objects
or lists of target objects.
Target objects represent skippable steps of the analysis pipeline
as described at \url{https://books.ropensci.org/targets/}.
Please read the walkthrough at
\url{https://books.ropensci.org/targets/walkthrough.html}
to understand the role of target objects in analysis pipelines.

For developers,
\url{https://wlandau.github.io/targetopia/contributing.html#target-factories}
explains target factories (functions like this one which generate targets)
and the design specification at
\url{https://books.ropensci.org/targets-design/}
details the structure and composition of target objects.
}

\examples{
if (identical(Sys.getenv("TAR_LONG_EXAMPLES"), "true")) {
targets::tar_dir({ # tar_dir() runs code from a temporary directory.
targets::tar_script({
  # Do not use temp files in real projects
  # or else your targets will always rerun.
  paths <- unlist(replicate(4, tempfile()))
  file.create(paths)
  list(
    tarchetypes::tar_files_input_raw(
      "x",
      paths,
      batches = 2
    )
  )
})
targets::tar_make()
targets::tar_read(x)
targets::tar_read(x, branches = 1)
})
}
}
\seealso{
Other Dynamic branching over files: 
\code{\link{tar_files_input}()},
\code{\link{tar_files_raw}()},
\code{\link{tar_files}()}
}
\concept{Dynamic branching over files}
