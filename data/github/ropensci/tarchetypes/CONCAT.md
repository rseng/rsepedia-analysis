
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
