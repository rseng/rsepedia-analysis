
# targets <img src='man/figures/logo.png' align="right" height="139"/>

[![ropensci](https://badges.ropensci.org/401_status.svg)](https://github.com/ropensci/software-review/issues/401)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.02959/status.svg)](https://doi.org/10.21105/joss.02959)
[![zenodo](https://zenodo.org/badge/200093430.svg)](https://zenodo.org/badge/latestdoi/200093430)
[![R
Targetopia](https://img.shields.io/badge/R_Targetopia-member-blue?style=flat&labelColor=gray)](https://wlandau.github.io/targetopia/)
[![CRAN](https://www.r-pkg.org/badges/version/targets)](https://CRAN.R-project.org/package=targets)
[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/targets/workflows/check/badge.svg)](https://github.com/ropensci/targets/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/targets/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/targets)
[![lint](https://github.com/ropensci/targets/workflows/lint/badge.svg)](https://github.com/ropensci/targets/actions?query=workflow%3Alint)

The `targets` package is a
[Make](https://www.gnu.org/software/make/)-like pipeline toolkit for
Statistics and data science in R. With `targets`, you can maintain a
reproducible workflow without repeating yourself. `targets` skips costly
runtime for tasks that are already up to date, runs the necessary
computation with implicit parallel computing, and abstracts files as R
objects. A fully up-to-date `targets` pipeline is tangible evidence that
the output aligns with the code and data, which substantiates trust in
the results.

## Prerequisites

1.  Familiarity with the [R programming
    language](https://www.r-project.org/), covered in [R for Data
    Science](https://r4ds.had.co.nz/).
2.  [Data science workflow management
    techniques](https://rstats.wtf/index.html).
3.  [How to write functions](https://r4ds.had.co.nz/functions.html) to
    prepare data, analyze data, and summarize results in data analysis
    projects.

## How to get started

1.  Watch minutes 6 through 40 of the [New York Open Statistical
    Programming Meetup from December
    2020](https://youtu.be/Gqn7Xn4d5NI).
2.  Read the [short walkthrough
    chapter](https://books.ropensci.org/targets/walkthrough.html) of the
    [user manual](https://books.ropensci.org/targets/).
3.  Sign up for a free [RStudio Cloud](https://rstudio.cloud) account
    and [click here](https://rstudio.cloud/project/1430691) to open the
    [walkthrough](https://books.ropensci.org/targets/walkthrough.html)
    code. Experiment with functions
    [`tar_make()`](https://docs.ropensci.org/targets/reference/tar_make.html)
    and
    [`tar_read()`](https://docs.ropensci.org/targets/reference/tar_read.html).
4.  Log into the [cloud
    workspace](https://rstudio.cloud/project/1699460) of the [official
    `targets` short
    course](https://github.com/wlandau/targets-tutorial/blob/main/README.md).
    Work through the exercises in R notebooks
    [`1-functions.Rmd`](https://github.com/wlandau/targets-tutorial/blob/main/1-functions.Rmd),
    [`2-pipelines.Rmd`](https://github.com/wlandau/targets-tutorial/blob/main/2-pipelines.Rmd),
    and
    [`3-changes.Rmd`](https://github.com/wlandau/targets-tutorial/blob/main/3-changes.Rmd).
5.  Try out one of the other [example
    projects](https://docs.ropensci.org/targets/index.html#example-projects)
    linked from the [reference
    website](https://docs.ropensci.org/targets/index.html#example-projects).

## Installation

| Type        | Source   | Command                                                           |
|-------------|----------|-------------------------------------------------------------------|
| Release     | CRAN     | `install.packages("targets")`                                     |
| Development | GitHub   | `remotes::install_github("ropensci/targets")`                     |
| Development | rOpenSci | `install.packages("targets", repos = "https://dev.ropensci.org")` |

## Recorded talks

### English

-   [R/Medicine 2021 (15.33)](https://youtu.be/HJI5mQJRGpY)
-   [R/Pharma 2020
    (9:24)](https://www.youtube.com/watch?v=GRqKJBaC5g4&list=PLMtxz1fUYA5C0YflXsR8EEAQXfjntlV1H&index=6)
-   [LA R Users Meetup, October 2020
    (1:14:40)](https://www.youtube.com/watch?v=Qq25BUxpJu4)
-   [New York Open Statistical Programming Meetup, December 2020
    (1:54:28)](https://youtu.be/Gqn7Xn4d5NI)
-   [ds-incubator series,
    2021](https://www.youtube.com/playlist?list=PLvgdJdJDL-APJqHy5CXs6m4N7hUVp5rb4)
-   [Lille R User Group, June 2021
    (45:54)](https://youtu.be/FODSavXGjYg)

### Español

-   [R-Ladies Barcelona, 2021-05-25
    (1:25:12)](https://www.youtube.com/watch?v=Vj312AfdpBo).

## Documentation

-   [User manual](https://books.ropensci.org/targets/): in-depth
    discussion about how to use `targets`.
-   [Reference website](https://docs.ropensci.org/targets/): formal
    documentation of all user-side functions, the statement of need, and
    multiple design documents of the internal architecture.
-   [Developer
    documentation](https://books.ropensci.org/targets-design/): software
    design documents for developers contributing to the deep internal
    architecture of `targets`.

## Courses

-   [Official half-day interactive
    tutorial](https://github.com/wlandau/targets-tutorial).

## Example projects

| Description                                                                                                        | Link                                                |
|--------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------|
| Minimal example                                                                                                    | <https://github.com/wlandau/targets-minimal>        |
| Machine learning with Keras                                                                                        | <https://github.com/wlandau/targets-keras>          |
| Validating a minimal Stan model                                                                                    | <https://github.com/wlandau/targets-stan>           |
| Using Target Markdown and `stantargets` to validate a Bayesian longitudinal model for clinical trial data analysis | <https://github.com/wlandau/rmedicine2021-pipeline> |
| Shiny app that runs a pipeline                                                                                     | <https://github.com/wlandau/targets-shiny>          |
| Deploying a pipeline to RStudio Connect                                                                            | <https://github.com/sol-eng/targets-deployment-rsc> |

## Apps

-   [`tar_watch()`](https://docs.ropensci.org/targets/reference/tar_watch.html):
    a built-in Shiny app to visualize progress while a pipeline is
    running. Available as a Shiny module via
    [`tar_watch_ui()`](https://docs.ropensci.org/targets/reference/tar_watch_ui.html)
    and
    [`tar_watch_server()`](https://docs.ropensci.org/targets/reference/tar_watch_server.html).
-   [`targetsketch`](https://wlandau.shinyapps.io/targetsketch): a Shiny
    app to help sketch pipelines
    ([app](https://wlandau.shinyapps.io/targetsketch),
    [source](https://github.com/wlandau/targetsketch)).

## Deployment

-   <https://solutions.rstudio.com/r/workflows/> explains how to deploy
    a pipeline to RStudio Connect ([example
    code](https://github.com/sol-eng/targets-deployment-rsc)).
-   [`tar_github_actions()`](https://docs.ropensci.org/targets/reference/tar_github_actions.html)
    sets up a pipeline to run on GitHub Actions. The [minimal
    example](https://github.com/wlandau/targets-minimal) demonstrates
    this approach.

## Extending and customizing targets

-   [R Targetopia](https://wlandau.github.io/targetopia/): a collection
    of [R packages](https://wlandau.github.io/targetopia/packages.html)
    that extend `targets`. [These
    packages](https://wlandau.github.io/targetopia/packages.html)
    simplify pipeline construction for specific fields of Statistics and
    data science.
-   [Target
    factories](https://wlandau.github.io/targetopia/contributing.html#target-factories):
    a programming technique to write specialized interfaces for custom
    pipelines. Posts
    [here](https://ropensci.org/blog/2021/02/03/targets/) and
    [here](https://wlandau.github.io/targetopia/contributing.html)
    describe how.

## Help

-   Post to the [GitHub discussion
    forum](https://github.com/ropensci/targets/discussions) to ask
    questions. To get the best help about a specific issue, create a
    reproducible example with
    [`targets::tar_reprex()`](https://docs.ropensci.org/targets/reference/tar_reprex.html)
    or
    [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html).
-   The [RStudio Community](https://community.rstudio.com/) forum is
    full of friendly enthusiasts of R and the tidyverse. Use the
    [`targets` tag](https://community.rstudio.com/tag/targets).
-   [Stack Overflow](https://stackoverflow.com/) broadcasts to the
    entire open source community. Use the [`targets-r-package`
    tag](https://stackoverflow.com/questions/tagged/targets-r-package).

## Code of conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/).

## Citation

``` r
citation("targets")
#> 
#> To cite targets in publications use:
#> 
#>   Landau, W. M., (2021). The targets R package: a dynamic Make-like
#>   function-oriented pipeline toolkit for reproducibility and
#>   high-performance computing. Journal of Open Source Software, 6(57),
#>   2959, https://doi.org/10.21105/joss.02959
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {The targets R package: a dynamic Make-like function-oriented pipeline toolkit for reproducibility and high-performance computing},
#>     author = {William Michael Landau},
#>     journal = {Journal of Open Source Software},
#>     year = {2021},
#>     volume = {6},
#>     number = {57},
#>     pages = {2959},
#>     url = {https://doi.org/10.21105/joss.02959},
#>   }
```
# targets 0.10.0.9000

## New features

* Implement `tar_assert_finite()`.
* Report the total runtime of the pipeline in the `"verbose"`, `"verbose_positives"`, `"timestamp"`, and `"timesamp_positives"` reporters.
* Add a `zoom_speed` argument to `tar_visnetwork()` and `tar_glimpse()` (#749, @dipterix).

# targets 0.10.0

## Bug fixes

* Add class `"tar_nonexportable"` to `format = "aws_keras"` and `format = "aws_torch"` stores.
* Export S3 methods of generic `tar_make_interactive_load_target()`.

## New features

* Allow entirely custom storage formats through `tar_target(format = tar_format(...))` (#736).
* Add a new function `tar_call()` to return the `targets` function currently running (from `_targets.R` or a target).
* Add a new function `tar_active()` to tell whether the pipeline is currently running. Detects if it is called from `tar_make()` or similar function.

## Enhancements

* Add `Sys.getenv("TAR_PROJECT")` to the output of `tar_envvars()`.
* Set the `store` field of `tar_runtime` prior to sourcing `_targets.R` so `tar_store()` works in target scripts.
* Explicitly export all the environment variables from `tar_envvars()` to targets run on parallel workers.
* Allow `format = "file"` targets to return `character(0)` (#728, @programLyrique).
* Automatically remove non-targets from the target list and improve target list error messages (#731, @billdenney).
* Link to resources on deploying to RStudio Connect (#745, @ian-flores).

# targets 0.9.1

* Mask pointers in function dependencies (#721, @matthiaskaeding)

# targets 0.9.0

## Highlights

* Track the version ID of AWS S3-backed targets if the bucket is version-enabled (#711). If you put your targets in AWS and the metadata and code under version control, you can `git checkout` a different branch of your code and all you targets will stay up to date.
* Refactor the AWS path format internally. It now consists of arbitrarily extensible key-value pairs so more AWS S3 functionality may be added more seamlessly going forward (#711).
* Switch the AWS S3 backend to `paws` (#711).

## New features

* Add a `region` argument to `tar_resources_aws()` to allow the user to explicitly declare a region for each AWS S3 buckets (@caewok, #681). Different buckets can now have different regions. This feature required modifying the metadata path for AWS storage formats. Before, the first element of the path was simply the bucket name. Now, it is internally formatted like `"bucket=BUCKET:region=REGION"`, where `BUCKET` is the user-supplied bucket name and `REGION` is the user-supplied region name. The new `targets` is back-compatible with the old metadata format, but if you run the pipeline with `targets` >= 0.8.1.9000 and then downgrade to `targets` <= 0.8.1, any AWS targets will break.
* Add new reporters `timestamp_positives"` and `"verbose_positives"` that omit messages for skipped targets (@psanker, #683).
* Implement `tar_assert_file()`.
* Implement `tar_reprex()` for creating easier reproducible examples of pipelines.
* Implement `tar_store()` to get the path to the store of the currently running pipeline (#714, @MilesMcBain).
* Automatically write a `_targets/user/` folder to encourage `gittargets` users to put custom files there for data version control.

## Bug fixes

* Make sure `tar_path()` uses the current store path of the currently running pipeline instead of `tar_config_get("store")` (#714, @MilesMcBain).

## Enhancements

* Refactor the automatic `.gitignore` file inside the data store to allow the metadata to be committed to version control more easily (#685, #711).
* Document target name requirements in `tar_target()` and `tar_target_raw()` (@tjmahr, #679).
* Catch and relay any the error if a target cannot be checked in `target_should_run.tar_builder()`. These kinds of errors sometimes come up with AWS storage.
* Fix the documentation of the reporters.
* Only write `_targets/.gitignore` for new data stores so the user can delete the `.gitignore` file without it mysteriously reappearing (#685).

# targets 0.8.1

## New features

* Add arguments `strict` and `silent` to allow `tar_load()` and `tar_load_raw()` to bypass targets that cannot be loaded.

## Enhancements

* Improve `tidyselect` docs in `tar_make()` (#640, @dewoller).
* Use namespaced call to `tar_dir()` in `tar_test()` (#642, @billdenney).
* Improve `tar_assert_target_list()` error message (@kkami1115, #654).
* Throw an informative error if a target name starts with a dot (@dipterix, #662).
* Improve help files of `tar_destroy()` and related cleanup functions (@billdenney, #675).

# targets 0.8.0

## Bug fixes

* Hash the correct files in `tar_target(target_name, ..., format = "aws_file")`. Previously, `_targets/objects/target_name` was also hashed if it existed.

## New features

* Implement a new `tar_config_unset()` function to delete one or more configuration settings from the YAML configuration file.
* Implement the `TAR_CONFIG` environment variable to set the default file path of the YAML configuration file with project settings (#622, @yyzeng, @atusy, @nsheff, @wdkrnls). If `TAR_CONFIG` is not set, the file path is still `_targets.yaml`.
* Restructure the YAML configuration file format to handle configuration information for multiple projects (using the `config` package) and support the `TAR_PROJECT` environment variable to select the current active project for a given R session. The old single-project format is gracefully deprecated (#622, @yyzeng, @atusy, @nsheff, @wdkrnls).
* Implement `retrieval = "none"` and `storage = "none"` to anticipate loading/saving targets from other languages, e.g. Julia (@MilesMcBain).
* Add a new `tar_definition()` function to get the target definition object of the current target while that target is running in a pipeline.
* If called inside an AWS target, `tar_path()` now returns the path to the staging file instead of `_targets/objects/target_name`. This ensures you can still write to `tar_path()` in `storage = "none"` targets and the package will automatically hash the right file and upload it to the cloud. (This behavior does not apply to formats `"file"` and `"aws_file"`, where it is never necessary to set `storage = "none"`.)

## Enhancements

* Use `eval(parse(text = ...), envir = tar_option_set("envir")` instead of `source()` in the `_targets.R` file for Target Markdown.
* Allow feather and parquet formats to accept objects of class `RecordBatch` and `Table` (@MilesMcBain).
* Let `knitr` load the Target Markdown engine (#469, @nviets, @yihui). Minimum `knitr` version is now `1.34`.
* In the `tar_resources_future()` help file, encourage the use of `plan` to specify resources.

# targets 0.7.0

## Bug fixes

* Ensure `error = "continue"` does not cause errored targets to have `NULL` values.
* Relay output and messages in Target Markdown interactive mode (using the R/default `knitr` engine).

## New features

* Expose the `poll_connection`, `stdout`, and `stderr` arguments of `callr::r_bg()` in `tar_watch()` (@mpadge).
* Add new helper functions to list targets in each progress category: `tar_started()`, `tar_skipped()`, `tar_built()`, `tar_canceled()`, and `tar_errored()`.
* Add new helper functions `tar_interactive()`, `tar_noninteractive()`, and `tar_toggle()` to differentially suppress code in non-interactive and interactive mode in Target Markdown (#607, @33Vito).

## Enhancements

* Handle `future` errors within targets (#570, @stuvet).
* Handle storage errors within targets (#571, @stuvet).
* In Target Markdown in non-interactive mode, suppress messages if the `message` `knitr` chunk option is `FALSE` (#574, @jmbuhr).
* In Target Markdown, if `tar_interactive` is not set, choose interactive vs non-interactive mode based on `isTRUE(getOption("knitr.in.progress"))` instead of `interactive()`.
* Convert errors loading dependencies into errors running targets (@stuvet).

# targets 0.6.0

## Bug fixes

* Allow `tar_poll()` to lose and then regain connection to the progress file.
* Make sure changes to the `tar_group` column of `iteration = "group"` data frames do not invalidate slices (#507, @lindsayplatt).

## New features

* In Target Markdown, add a new `tar_interactive` global option to select interactive mode or non-interactive mode (#469).
* Highlight a graph neighborhood when the user clicks a node. Control the neighborhood degree with new arguments `degree_from` and `degree_to` of `tar_visnetwork()` and `tar_glimpse()` (#474, @rgayler).
* Make the target script path configurable in `tar_config_set()` (#476).
* Add a `tar_script` chunk option in Target Markdown to control where the `{targets}` language engine writes the target script and helper scripts (#478).
* Add new arguments `script` and `store` to choose custom paths to the target script file and data store for individual function calls (#477).
* Allow users to set an alternative path to the YAML configuration file for the current R session (#477). Most users have no reason to set this path, it is only for niche applications like Shiny apps with `targets` backends. Unavoidably, the path gets reset to `_targets.yaml` when the session restarts.
* Add new `_targets.yaml` config options `reporter_make`, `reporter_outdated`, and `workers` to control function argument defaults shared across multiple functions called outside `_targets.R` (#498, @ianeveperry).
* Add `tar_load_globals()` for debugging, testing, prototyping, and teaching (#496, @malcolmbarrett).
* Add structure to the `resources` argument of `tar_target()` to avoid conflicts among formats and HPC backends (#489). Includes user-side helper functions like `tar_resources()` and `tar_resources_aws()` to build the required data structures.
* Log skipped targets in `_targets/meta/progress` and display then in `tar_progress()`, `tar_poll()`, `tar_watch()`, `tar_progress_branches()`, `tar_progress_summary()`, and `tar_visnetwork()` (#514). Instead of writing each skip line separately to `_targets/meta/progress`, accumulate skip lines in a queue and then write them all out in bulk when something interesting happens. This avoids a lot of overhead in certain cases.
* Add a `shortcut` argument to `tar_make()`, `tar_make_clustermq()`, `tar_make_future()`, `tar_outdated()`, and `tar_sitrep()` to more efficiently skip parts of the pipeline (#522, #523, @jennysjaarda, @MilesMcBain, @kendonB).
* Support `names` and `shortcut` in graph data frames and graph visuals (#529).
* Move `allow` and `exclude` to the network behind the graph visuals rather than the visuals themselves (#529).
* Add a new "progress" display to the `tar_watch()` app to show verbose progress info and metadata.
* Add a new `workspace_on_error` argument of `tar_option_set()` to supersede `error = "workspace"`. Helps control workspace behavior independently of the `error` argument of `tar_target()` (#405, #533, #534, @mattwarkentin, @xinstein).
* Implement `error = "abridge"` in `tar_target()` and related functions. If a target errors out with this option, the target itself stops, any currently running targets keeps, and no new targets launch after that (#533, #534, @xinstein).
* Add a menu prompt to `tar_destroy()` which can be suppressed with `TAR_ASK = "false"` (#542, @gofford).
* Support functions `tar_older()` and `tar_newer()` to help users identify and invalidate targets at regular times or intervals.

## Deprecations

* In Target Markdown, deprecate the `targets` chunk option in favor of `tar_globals` (#469).
* Deprecate `error = "workspace"` in `tar_target()` and related functions. Use `tar_option_set(workspace_on_error = TRUE)` instead (#405, #533, @mattwarkentin, @xinstein).

## Performance

* Reset the backoff upper bound when concluding a target or shutting down a `clustermq` worker (@rich-payne).
* Set more aggressive default backoff bound of 0.1 seconds (previous: 5 seconds) and set a more aggressive minimum of 0.001 seconds (previous: 0.01 seconds) (@rich-payne).
* Speed up the summary and forecast reporters by only printing to the console every quarter second.
* Avoid superfluous calls to `store_sync_file_meta.default()` on small files.
* In `tar_watch()`, take several measures to avoid long computation times rendering the graph:
    * Expose arguments `display` and `displays` to `tar_watch()` so the user can select which display shows first.
    * Make `"summary"` the default display instead of `"graph"`.
    * Set `outdated` to `FALSE` by default.

## Enhancements

* Simplify the Target Markdown example.
* Warn about unnamed chunks in Target Markdown.
* Redesign option system to be more object-oriented and rigorous. Also export most options to HPC workers (#475).
* Simplify config system to let API function arguments take control (#483).
* In `tar_read()` for targets with `format = "aws_file"`, download the file back to the path the user originally saved it when the target ran.
* Replace the `TAR_MAKE_REPORTER` environment variable with `targets::tar_config_get("reporter_make")`.
* Use `eval(parse(text = readLines("_targets.R")), envir = some_envir)` and related techniques instead of the less controllable `source()`. Expose an `envir` argument to many functions for further control over evaluation if `callr_function` is `NULL`.
* Drop `out.attrs` when hashing groups of data frames to extend #507 to `expand.grid()` (#508).
* Increase the number of characters in errors and warnings up to 2048.
* Refactor assertions to automatically generate better messages.
* Export assertions, conditions, and language utilities in packages that build on top of `targets`.
* Change `GITHUBPAT` to `GITHUB_TOKEN` in the `tar_github_actions()` YAML file (#554, @eveyp).
* Support the `eval` chunk option in Target Markdown (#552, @fkohrt).
* Record time stamps in the metadata `time` column for all builder targets, regardless of storage format.

# targets 0.5.0

## Bug fixes

* Export in-memory config settings from `_targets.yaml` to parallel workers.

## New features

* Add a limited-scope `exclude` argument to `tar_watch()` and `tar_watch_server()` (#458, @gorkang).
* Write a `.gitignore` file to ignore everything in `_targets/meta/` except `.gitignore` and `_targets/meta/meta`.
* Target Markdown: add `knitr` engines for pipeline construction and prototyping from within literate programming documents (#469, @cderv, @nviets, @emilyriederer, @ijlyttle, @GShotwell, @gadenbuie, @tomsing1). Huge thanks to @cderv on this one for answering my deluge of questions, helping me figure out what was and was not possible in `knitr`, and ultimately circling me back to a successful approach.
* Add an RStudio R Markdown template for Target Markdown (#469).
* Implement `use_targets()`, which writes the Target Markdown template to the project root (#469).
* Implement `tar_unscript()` to clean up scripts written by Target Markdown.

## Enhancements

* Enable priorities in `tar_make()` and `tar_manifest()`.
* Show the priority in the print method of stem and pattern targets.
* Throw informative errors if the secondary arguments to `pattern = slice()` or `pattern = sample()` are invalid.
* In `tar_target_raw()`, assert that commands have length 1 when converted to expressions.
* Handle errors and post failure artifacts in the Github Actions YAML file.
* Rewrite the documentation on invalidation rules in `tar_cue()` (@maelle).
* Drop `dplyr` groups and `"grouped_df"` class in `tar_group()` (`tarchetypes` discussion #53, @kendonB).
* Assign branch names to dynamic branching return values produced by `tar_read()` and `tar_read_raw()`.

# targets 0.4.2

## Bug fixes

* Do not use time stamps to monitor the config file (e.g. `_targets.yaml`). Fixes CRAN check errors from version 0.4.1.

# targets 0.4.1

* Fix CRAN test error on Windows R-devel.
* Do not inherit `roxygen2` docstrings from `shiny`.
* Handle more missing `Suggests:` packages.
* Unset the config lock before reading `targets.yaml` in the `callr` process.

# targets 0.4.0

## Bug fixes

* Avoid `file.rename()` errors when migrating staged temporary files (#410).
* Return correct error messages from feather and parquet formats (#388). Now calling `assert_df()` from `store_assert_format()` instead of `store_cast_object()`. And now those last two functions are not called at all if the target throws an error.
* Retry writing lines to database files so Windows machines can run `tar_poll()` at the same time as the pipeline (#393).
* Rename file written by `tar_renv()` to `_targets_packages.R` (#397).
* Ensure metadata is loaded to compute labels properly when `outdated = FALSE` in `tar_visnetwork()`.

## New features

* Implement `tar_timestamp()` and `tar_timestamp_raw()` to get the last modified timestamp of a target's data (#378).
* Implement `tar_progress_summary()` to compactly summarize all pipeline progress (#380).
* Add a `characters` argument of `tar_traceback()` to cap the traceback line lengths (#383).
* Add new "summary" and "about" views to `tar_watch()` (#382).
* Implement `tar_poll()` to repeatedly poll runtime progress in the R console (#381). `tar_poll()` is a lightweight alternative to `tar_watch()`.
* Change the color of the "dormant" status in the graph.
* Add a `tar_envvar()` function to list values of special environment variables supported in `targets`. The help file explains each environment variable in detail.
* Support extra project-level configuration settings with `_targets.yaml` (#297). New functions `tar_config_get()` and `tar_config_set()` interact with the `_targets.yaml` file. Currently only supports the `store` field to set the data store path to something other than `_targets/`.

## Performance

* Shut down superfluous persistent workers earlier in dynamic branching and when all remaining targets have `deployment = "main"` (#398, #399, #404, @pat-s).

## Enhancements

* Attempt to print only the useful part of the traceback in `tar_traceback()` (#383).
* Add a line break at the end of the "summary" reporter so warnings do not mangle the output.
* In `tar_watch()`, use `shinybusy` instead of `shinycssloaders` and keep current output on display while new output is rendering (#386, @rcorty).
* Right-align the headers and counts in the "summary" and "forecast" reporters.
* Add a timestamp to the "summary" reporter.
* Make the reporters show when a target ends (#391, @mattwarkentin).
* Make the reporters show when a pattern ends if the pattern built at least one target and none of the targets errored or canceled.
* Use words "start" and "built" in reporters.
* Use the region of the AWS S3 bucket instead of the local `AWS_DEFAULT_REGION` environment variable (`check_region = TRUE`; #400, @tomsing1).
* In `tar_meta()`, return `POSIXct` times in the time zone of the calling system (#131).
* Throw informative error messages when a target's name or command is missing (#413, @liutiming).
* Bring back ALTREP in `qs::qread()` now that `qs` 0.24.1 requires `stringfish` >= 1.5.0 (#147, @glep).
* Relax dynamic branching checks so `pattern = slice(...)` can take multiple indexes (#406, #419, @djbirke, @alexgphayes)

# targets 0.3.1

## Bug fixes

* `queue$enqueue()` is now `queue$prepend()` and always appends to the front of the queue (#371).

## Enhancements

* Throw a warning if `devtools::load_all()` or similar is detected inside `_targets.R` (#374).

## CRAN

* Skip `feather` and `parquet` tests on CRAN.

# targets 0.3.0

## Bug fixes

* Fix the "write target at cursor" RStudio addin and move cursor between the parentheses.

## New features

* Add a `backoff` option in `tar_option_set()` to set the maximum upper bound (seconds) for the polling interval (#333).
* Add a new `tar_github_actions()` function to write a GitHub Actions workflow file for continuous deployment of data analysis pipelines (#339, @jaredlander).
* Add a new `TAR_MAKE_REPORTER` environment variable to globally set the reporter of the `tar_make*()` functions (#345, @alexpghayes).
* Support new storage formats "feather", "parquet", "aws_feather", and "aws_parquet" (#355, @riazarbi).

## Performance

* Implement an exponential backoff algorithm for polling the priority queue in `tar_make_clustermq()` and `tar_make_future()` (#333). 
* In `tar_make_future()`, try to submit a target every time a worker is polled.
* In `tar_make_future()`, poll workers in order of target priority.
* Avoid the time delay in exiting on error (from https://github.com/r-lib/callr/issues/185).
* Clone target objects for the pipeline and scrape more `targets` internal objects out of the environment in order to avoid accidental massive data transfers to workers.

## Enhancements

* Use `rlang::check_installed()` inside `assert_package()` (#331, @malcolmbarrett).
* Allow `tar_destroy(destroy = "process")`.
* In `tar_watch()`, increase default `seconds` to 15 (previously 5).
* In `tar_watch()`, debounce instead of throttle inputs.
* In `tar_watch()`, add an action button to refresh the outputs.
* Always deduplicate metadata after `tar_make()`. Will help compute a cache key on GitHub Actions and similar services.
* Deprecate `tar_deduplicate()` due to the item above.
* Reorder information in timestamped messages.
* Document RNG seed generation in `tar_target_raw()`, `tar_meta()`, and `tar_seed()` (#357, @alexpghayes).
* Switch meaning of `%||%` and `%|||%` to conform to historical precedent.
* Only show a command line spinner if `reporter = "silent"` (#364, @matthiasgomolka).
* Target and pipeline objects no longer have an `envir` element.

# targets 0.2.0

## Bug fixes

* In `tar_load()`, subset metadata to avoid accidental attempts to load global objects in `tidyselect` calls.
* Do not register a pattern as running unless an actual branch is about to start (#304).
* Use a name spec in `vctrs::vec_c()` (#320, @joelnitta).

## New features

* Add a new `names` argument to `tar_objects()` and `tar_workspaces()` with `tidyselect` functionality.
* Record info on the main process (PID, R version, `targets` version) in `_targets/meta/process` and write new functions `tar_process()` and `tar_pid()` to retrieve the data (#291, #292).
* Add a new `targets_only` argument to `tar_meta()`.
* Add new functions `tar_helper()` and `tar_helper_raw()` to write general-purpose R scripts, using tidy evaluation for as a template mechanism (#290, #291, #292, #306).
* Export functions to check the existence of various pieces of local storage: `tar_exist_meta()`, `tar_exist_objects()`, `tar_exist_progress()`, `tar_exist_progress()`, `tar_exist_script()` (#310).
* Add a new `supervise` argument to `tar_watch()`.
* Add a new `complete_only` argument to `tar_meta()` to optionally return only complete rows (no `NA` values).
* Catch `callr` errors and refer users to the debugging chapter of the manual.

## Enhancements

* Improve error messages of invalid arguments (#298, @brunocarlin). Removes partial argument matching in most cases.
* By default, locally enable `crayon` if an only if the calling process is interactive (#302, @ginolhac). Can still be disabled with `options(crayon.enabled = FALSE)` in `_targets.R`.
* Improve error handling and message for `format = "url"` when the HTTP response status code is not 200 (#303, @petrbouchal).
* Add more `extras` packages to `tar_renv()` (to support `tar_watch()`).
* Show informative message instead of error in `tar_watch()` if `_targets.R` does not exist.
* Clear up the documentation of the `names` argument of `tar_load()` (#314, @jameelalsalam).
* Do not override `nobody` in custom `curl` handles (#315, @riazarbi).
* Rename "running" to "started" in the progress metadata. This avoids the implicit claim that `targets` is somehow actively monitoring each job, e.g. through a connection or heartbeat (#318).
* Set `errormode = "warn"` in `getVDigest()` for files to work around https://github.com/eddelbuettel/digest/issues/49 for network drives on Windows. `targets` already runs those file checks anyway. (#316, @boshek).
* If a package fails to load, print the library paths `targets` tried to load from.

# targets 0.1.0

## Bug fixes

* `tar_test()` now skips all tests on Solaris in order to fix the problems shown on the CRAN check page.
* Enable `allow` and `exclude` to work on imports in `tar_visnetwork()` and `tar_glimpse()`.
* Put `visNetwork` legends on right to avoid crowding the graph.

## Performance

* Call `force()` on subpipeline objects to eliminate high-memory promises in target objects. Allows targets to be deployed to workers much faster when `retreival` is `"main"` (#279).

## New features

* Add a new box to the `tar_watch()` app to tabulate progress on dynamic branches (#273, @mattwarkentin).
* Store `type`, `parent`, and `branches` in progress data for `tar_watch()` (#273, @mattwarkentin).
* Add a `fields` argument in `tar_progress()` and default to `"progress"` for back compatibility (#273, @mattwarkentin).
* Add a new `tar_progress_branches()` function to tabulate branch progress (#273, @mattwarkentin).
* Add new "refresh" switch to `tar_watch()` to toggle automatic refreshing and force a refresh.

## Enhancements

* Exclude `.Random.seed` by default in `tar_visnetwork()`.
* Spelling: "cancelled" changed to "canceled".
* Enhance controls and use of space in the `tar_watch()` app.
* Centralize internal path management utilities.

## Configuration

* Skip `clustermq` tests on Solaris.

# targets 0.0.2

## CRAN response

* Avoid starting the description with the package name.
* Remove `if(FALSE)` blocks from help files to fix "unexecutable code" warnings (`tar_glimpse()`, `tar_visnetwork()`, and `tar_watch()`).
* Remove commented code in the examples (`tar_edit()`, `tar_watch_ui()`, and `tar_watch_server()`).
* Ensure that all examples, tests, and vignettes do not write to the user's home file space. (Fixed an example of `tar_workspace()`.)

## Enhancements

* Use JOSS paper in `CITATION`.

# targets 0.0.1

## Enhancements

* Accept lists of target objects at the end of `_targets.R` (#253).
* Deprecate `tar_pipeline()` and `tar_bind()` because of the above (#253).
* Always show a special message when the pipeline finishes (#258, @petrbouchal).
* Disable `visNetwork` stabilization (#264, @mattwarkentin).
* Use default `visNetwork` font size.
* Relay errors as condition messages if `error` is `"continue"` (#267, @liutiming).

# targets 0.0.0.9003

## Bug fixes

* Ensure pattern-only pipelines can be defined so they can be combined again later with `tar_bind()` (#245, @yonicd).
* Implement safeguards around `igraph` topological sort.

## Enhancements

* Topologically sort the rows of `tar_manifest()` (#263, @sctyner).

## Breaking changes

* Make patterns composable (#212, @glep, @djbirke).
* Allow workspaces to load nonexportable objects (#214).
* Make workspace files super light by saving only a reference to the required dependencies (#214).
* Add a new `workspaces` argument to `tar_option_set()` to specify which targets will save their workspace files during `tar_make()` (#214).
* Change `error = "save"` to `error = "workspace"` to so it is clearer that saving workspaces no longer duplicates data (#214).
* Rename `what` to `destroy` in `tar_destroy()`.
* Remove `tar_undebug()` because is redundant with `tar_destroy(destroy = "workspaces")`.

## New features

* Make patterns composable (#212).
* Add new dynamic branching patterns `head()`, `tail()`, and `sample()` to provide functionality equivalent to `drake`'s `max_expand` (#56).
* Add a new `tar_pattern()` function to emulate dynamic branching outside a pipeline.
* Add a new `level_separation` argument to `tar_visnetwork()` and `tar_glimpse()` to control the aspect ratio (#226).
* Track functions from multiple packages with the `imports` argument to `tar_option_set()` (#239).
* Add color for "built" progress if `outdated` is `FALSE` in `tar_visnetwork()`.
* Tweak colors in `tar_visnetwork()` to try to account for color blindness.

## Enhancements

* Return full patterns from `tar_manifest()`.
* Record package load errors in progress and metadata (#228, @psychelzh).
* `tar_renv()` now invokes `_targets.R` through a background process just like `tar_outdated()` etc. so it can account for more hidden packages (#224, @mattwarkentin).
* Set `deployment` equal to `"main"` for all targets in `tar_make()`. This ensures `tar_make()` does not waste time waiting for nonexistent files to ship over a nonexistent network file system (NFS). `tar_make_clustermq()` or `tar_make_future()` could use NFS, so they still leave `deployment` alone.

# targets 0.0.0.9002

## Breaking changes

* Add a new `size` field to the metadata to allow `targets` to make better judgments about when to rehash files (#180). We now compare hashes to check file size differences instead of doing messy floating point comparisons with ad hoc tolerances. It breaks back compatibility with old projects, but the error message is informative, and this is all still before the first official release.
* Change "local" to "main" and "remote" to "worker" in the `storage`, `retrieval`, and `deployment` settings (#183, @mattwarkentin).
* Ensure function dependencies are sorted before computing the function hash (GitHub commit f15face7d72c15c2d1098da959492bdbfcddb425).
* Move `garbage_collection` to a target-level setting, i.e. argument to `tar_target()` and `tar_option_set()` (#194). Previously was an argument to the `tar_make*()` functions.
* Allow `tar_name()` and `tar_path()` to run outside the pipeline with debugging-friendly default return values.

## Bug fixes

* Stop sending target return values over the network when `storage` is `"remote"` (#182, @mattwarkentin).
* Shorten lengths of warnings and error messages to 128 characters (#186, @gorkang).
* Restrict in-memory metadata to avoid incorrectly recycling deleted targets (#191).
* Marshal nonexportable dependencies before sending them to workers. Transport data through `target$subpipeline` rather than `target$cache` to make that happen (#209, @mattwarkentin).

## New features

* Add a new function `tar_bind()` to combine pipeline objects.
* Add `tar_seed()` to get the random number generator seed of the target currently running.

## Enhancements

* Allow target-specific `future::plan()`s through the `resources` argument of `tar_target()` (#198, @mattwarkentin).
* Use `library()` instead of `require()` in `command_load_packages()`.
* Evaluate commands directly in `targets$cache$targets$envir` to improve convenience in interactive debugging (`ls()` just works now.) This is reasonably safe now that the cache is populated at the last minute and cleared as soon as possible (#209, #210).

# targets 0.0.0.9000

* First version.
# Contributing

Development is a community effort, and we welcome participation.

## Code of Conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). 

## Discussions

At <https://github.com/ropensci/targets/discussions>, you can post general questions, brainstorm ideas, and ask for help.

## Issues

<https://github.com/ropensci/targets/issues> is for bug reports, performance issues, maintenance tasks, and feature requests. When you post, please abide by the following guidelines.

* Before posting a new issue or discussion topic, please take a moment to search for existing similar threads in order to avoid duplication.
* For bug reports: if you can, please install the latest GitHub version of `targets` (i.e. `remotes::install_github("ropensci/targets")`) and verify that the issue still persists.
* Describe your issue in prose as clearly and concisely as possible.
* For any problem you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot. A reproducible example is:
    * **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Development

External code contributions are extremely helpful in the right circumstances. Here are the recommended steps.

1. Prior to contribution, please propose your idea in a discussion topic or issue thread so you and the maintainer can define the intent and scope of your work.
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
* Describe your contribution in the project's [`NEWS.md`](https://github.com/ropensci/targets/blob/main/NEWS.md) file. Be sure to mention relevent GitHub issue numbers and your GitHub name as done in existing news entries.
* If you feel contribution is substantial enough for official author or contributor status, please add yourself to the `Authors@R` field of the [`DESCRIPTION`](https://github.com/ropensci/targets/blob/main/DESCRIPTION) file.
# Prework

* [ ] I understand and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and the [contributing guidelines](https://github.com/ropensci/targets/blob/main/CONTRIBUTING.md).
* [ ] I have already submitted a [discussion topic](https://github.com/ropensci/targets/discussions) or [issue](https://github.com/ropensci/targets/issues) to discuss my idea with the maintainer.

# Related GitHub issues and pull requests

* Ref: #

# Summary

Please explain the purpose and scope of your contribution.
---
name: Maintenance
about: "Something in targets needs work: updates, documentation, etc. Not a bug, performance issue, or new feature."
title: ""
labels: "type: maintenance"
assignees: ""
---

## Prework

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/targets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/targets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] For any problems you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot. Convenient helpers for this include [`targets::tar_reprex()`](https://docs.ropensci.org/targets/reference/tar_reprex.html) and  [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html). A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the issue.

To help us read any code you include (optional) please try to follow the [tidyverse style guide](https://style.tidyverse.org/). The `style_text()` and `style_file()` functions from the [`styler`](https://github.com/r-lib/styler) package make it easier.

## Reproducible example

* [ ] For any problems you identify, post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot. Convenient helpers for this include [`targets::tar_reprex()`](https://docs.ropensci.org/targets/reference/tar_reprex.html) and  [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html). A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).
---
name: Bug
about: Please do not submit a bug report unless your issue is a genuine bug in targets and not a known limitation, usage error, or issue from another package that targets depends on.
title: ""
labels: "type: bug"
assignees: wlandau
---

## Prework

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/targets/blob/main/CONTRIBUTING.md).
* [ ] Confirm that your issue is a genuine bug in the `targets` package itself and not a user error, known limitation, or issue from another package that `targets` depends on. For example, if you get errors running `tar_make_clustermq()`, try isolating the problem in a reproducible example that runs `clustermq` and not `targets`. And for miscellaneous troubleshooting, please post to [discussions](https://github.com/ropensci/targets/discussions) instead of [issues](https://github.com/ropensci/targets/issues).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/targets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Using [`targets::tar_reprex()`](https://docs.ropensci.org/targets/reference/tar_reprex.html), [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html), or similar, post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the bug.

## Reproducible example

* [ ] Using [`targets::tar_reprex()`](https://docs.ropensci.org/targets/reference/tar_reprex.html), [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html), or similar, post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Expected result

What should have happened? Please be as specific as possible.

## Diagnostic information

* A [reproducible example](https://github.com/tidyverse/reprex).
* Session info, available through `sessionInfo()` or [`reprex(si = TRUE)`](https://github.com/tidyverse/reprex).
* A stack trace from `traceback()` or `rlang::trace_back()`.
* The [SHA-1 hash](https://git-scm.com/book/en/v1/Getting-Started-Git-Basics#Git-Has-Integrity) of the GitHub commit of `targets` currently installed. `packageDescription("targets")$GithubSHA1` shows you this.
---
name: New feature
about: Suggest a new feature.
title: ""
labels: "type: new feature"
assignees: wlandau
---

## Prework

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/targets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/targets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] New features take time and effort to create, and they take even more effort to maintain. So if the purpose of the feature is to resolve a struggle you are encountering personally, please consider first posting a "trouble" or "other" issue so we can discuss your use case and search for existing solutions first.
* [ ] Format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

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

* [ ] Read and agree to the [code of conduct](https://ropensci.org/code-of-conduct/) and [contributing guidelines](https://github.com/ropensci/targets/blob/main/CONTRIBUTING.md).
* [ ] If there is [already a relevant issue](https://github.com/ropensci/targets/issues), whether open or closed, comment on the existing thread instead of posting a new issue.
* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) like [this one](https://github.com/ropensci/targets/issues/256#issuecomment-754229683) so the maintainer can troubleshoot the problems you identify. Convenient helpers for this include [`targets::tar_reprex()`](https://docs.ropensci.org/targets/reference/tar_reprex.html) and  [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html). A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Description

Please describe the performance issue.

## Reproducible example

* [ ] Post a [minimal reproducible example](https://www.tidyverse.org/help/) so the maintainer can troubleshoot the problems you identify. Convenient helpers for this include [`targets::tar_reprex()`](https://docs.ropensci.org/targets/reference/tar_reprex.html) and  [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html). A reproducible example is:
    * [ ] **Runnable**: post enough R code and data so any onlooker can create the error on their own computer.
    * [ ] **Minimal**: reduce runtime wherever possible and remove complicated details that are irrelevant to the issue at hand.
    * [ ] **Readable**: format your code according to the [tidyverse style guide](https://style.tidyverse.org/).

## Benchmarks

How poorly does `targets` perform? To find out, we recommend you use the [`proffer`](https://github.com/r-prof/proffer) package and take screenshots of the results displayed in your browser.

```r
library(targets)
library(proffer)
px <- pprof({
  # All your targets code goes here.
})
```
---
title: 'The targets R package: a dynamic Make-like function-oriented pipeline toolkit for reproducibility and high-performance computing'
tags:
- R
- reproducibility
- high-performance computing
- pipeline
- workflow
- Make
date: "12 January 2021"
output: pdf_document
authors:
- name: William Michael Landau
  orcid: 0000-0003-1878-3253
  email: will.landau@gmail.com
  affiliation: 1
bibliography: paper.bib
affiliations:
- name: Eli Lilly and Company
  index: 1
---

# Summary

The [`targets`](https://github.com/ropensci/targets) R package [@targets] is a pipeline toolkit for computationally intense reproducible research. It reduces the time and effort required to develop a data analysis project and maintain a trustworthy set of results. [`targets`](https://github.com/ropensci/targets) uses static code analysis to detect dependency relationships among interconnected computational tasks and construct a directed acyclic graph (DAG), which researchers can visualize in order to understand and communicate the structure of a complicated workflow. To run the pipeline at scale, [`targets`](https://github.com/ropensci/targets) leverages implicit parallel computing and optional cloud storage. In subsequent runs, [`targets`](https://github.com/ropensci/targets), skips tasks that are already synchronized with their upstream dependencies, which not only reduces the runtime of rapidly developing workflows, but also provides tangible evidence of reproducibility.

In high-performance computing scenarios, [`targets`](https://github.com/ropensci/targets) uses its DAG to discern which targets can run concurrently and which targets are still waiting for other upstream targets to finish processing. As soon as a target's dependency requirements are met, the target is deployed to the next available parallel worker. Internally, [`targets`](https://github.com/ropensci/targets) leverages the [`clustermq`](https://github.com/mschubert/clustermq) package [@clustermq] for persistent workers and the  [`future`](https://github.com/HenrikBengtsson/future) package [@future] for transient workers. Both [`clustermq`](https://github.com/mschubert/clustermq) and [`future`](https://github.com/HenrikBengtsson/future) are powerful and versatile frameworks capable of submitting R workloads not only to multiple cores on a single machine, but also to popular resource managers on shared computing clusters.

[`targets`](https://github.com/ropensci/targets) is the successor to [`drake`](https://github.com/ropensci/drake) [@drake], which in turn originated from [`remake`](https://github.com/richfitz/remake) [@remake], an R package modeled after GNU Make [@Make]. Unlike Make, [`targets`](https://github.com/ropensci/targets) and [`drake`](https://github.com/ropensci/drake) and [`remake`](https://github.com/richfitz/remake) focus on the R language, encourage an idiomatic function-oriented style of programming, and abstract each target as an R object. Relative to [`remake`](https://github.com/richfitz/remake) and  [`drake`](https://github.com/ropensci/drake), [`targets`](https://github.com/ropensci/targets) is friendlier and more efficient, surpassing the permanent architectural limitations of both predecessors. The data storage system of [`targets`](https://github.com/ropensci/targets) is lighter and more transparent, which helps users diagnose issues, move projects to different file systems, work with multiple contributors, and leverage seamless [Metaflow](https://github.com/Netflix/metaflow)-like cloud storage integration [@metaflow]. In addition,  [`targets`](https://github.com/ropensci/targets) supports stronger user-side guardrails, more introspective dependency graph visualizations, parallel efficient dynamic branching, and an interface more amenable to metaprogramming and third-party extensions.

# References
The `tar_watch()` app visualizes the status and progress of a `targets` pipeline. 

## Views

* __summary__: overall runtime progress summary of the pipeline. In the table, the "time" column shows the last time a target started, finished, errored, or canceled itself. The "since" columns shows how long ago that was.
* __branches__: like the summary view except with specific information about dynamic branching progress.
* __progress__: a large searchable table of progress information and metadata.
* __graph__: show the `tar_visnetwork()` dependency graph. If this graph may be slow to refresh, consider toggling the outdated switch, speeding up your `_targets.R` file, or selecting another view.

## Controls

* __refresh__: click this button to force the displays in the app to refresh.
* __watch__: toggle whether the app refreshes automatically every few seconds.
* __targets_only__: whether to show only the targets or also show functions and other global objects.
* __outdated__: whether to color-code nodes depending on whether they are up to date. This feature may cause if the number of targets is enormous, so you may consider turning it off to just color by runtime progress in some cases.
* __label__: labels to append to the node names to optionally show the size, runtime, and number of branches of targets based on past recorded info in the metadata.
* __seconds__: how often to refresh the displays in the app when the "watch" switch is turned on.
* __level_separation__: how wide the graph should be.
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

# targets <img src='man/figures/logo.png' align="right" height="139"/>

[![ropensci](https://badges.ropensci.org/401_status.svg)](https://github.com/ropensci/software-review/issues/401)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.02959/status.svg)](https://doi.org/10.21105/joss.02959)
[![zenodo](https://zenodo.org/badge/200093430.svg)](https://zenodo.org/badge/latestdoi/200093430)
[![R Targetopia](https://img.shields.io/badge/R_Targetopia-member-blue?style=flat&labelColor=gray)](https://wlandau.github.io/targetopia/)
[![CRAN](https://www.r-pkg.org/badges/version/targets)](https://CRAN.R-project.org/package=targets)
[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![check](https://github.com/ropensci/targets/workflows/check/badge.svg)](https://github.com/ropensci/targets/actions?query=workflow%3Acheck)
[![codecov](https://codecov.io/gh/ropensci/targets/branch/main/graph/badge.svg?token=3T5DlLwUVl)](https://app.codecov.io/gh/ropensci/targets)
[![lint](https://github.com/ropensci/targets/workflows/lint/badge.svg)](https://github.com/ropensci/targets/actions?query=workflow%3Alint)

The `targets` package is a [Make](https://www.gnu.org/software/make/)-like pipeline toolkit for Statistics and data science in R. With `targets`, you can maintain a reproducible workflow without repeating yourself. `targets` skips costly runtime for tasks that are already up to date, runs the necessary computation with implicit parallel computing, and abstracts files as R objects. A fully up-to-date `targets` pipeline is tangible evidence that the output aligns with the code and data, which substantiates trust in the results.

## Prerequisites

1. Familiarity with the [R programming language](https://www.r-project.org/), covered in [R for Data Science](https://r4ds.had.co.nz/).
1. [Data science workflow management techniques](https://rstats.wtf/index.html).
1. [How to write functions](https://r4ds.had.co.nz/functions.html) to prepare data, analyze data, and summarize results in data analysis projects.

## How to get started

1. Watch minutes 6 through 40 of the [New York Open Statistical Programming Meetup from December 2020](https://youtu.be/Gqn7Xn4d5NI).
1. Read the [short walkthrough chapter](https://books.ropensci.org/targets/walkthrough.html) of the [user manual](https://books.ropensci.org/targets/).
1. Sign up for a free [RStudio Cloud](https://rstudio.cloud) account and [click here](https://rstudio.cloud/project/1430691) to open the [walkthrough](https://books.ropensci.org/targets/walkthrough.html) code. Experiment with functions [`tar_make()`](https://docs.ropensci.org/targets/reference/tar_make.html) and [`tar_read()`](https://docs.ropensci.org/targets/reference/tar_read.html).
1. Log into the [cloud workspace](https://rstudio.cloud/project/1699460) of the [official `targets` short course](https://github.com/wlandau/targets-tutorial/blob/main/README.md). Work through the exercises in R notebooks [`1-functions.Rmd`](https://github.com/wlandau/targets-tutorial/blob/main/1-functions.Rmd), [`2-pipelines.Rmd`](https://github.com/wlandau/targets-tutorial/blob/main/2-pipelines.Rmd), and [`3-changes.Rmd`](https://github.com/wlandau/targets-tutorial/blob/main/3-changes.Rmd).
1. Try out one of the other [example projects](https://docs.ropensci.org/targets/index.html#example-projects) linked from the [reference website](https://docs.ropensci.org/targets/index.html#example-projects).

## Installation

Type | Source | Command
---|---|---
Release | CRAN | `install.packages("targets")`
Development | GitHub | `remotes::install_github("ropensci/targets")`
Development | rOpenSci | `install.packages("targets", repos = "https://dev.ropensci.org")`

## Recorded talks

### English

* [R/Medicine 2021 (15.33)](https://youtu.be/HJI5mQJRGpY)
* [R/Pharma 2020 (9:24)](https://www.youtube.com/watch?v=GRqKJBaC5g4&list=PLMtxz1fUYA5C0YflXsR8EEAQXfjntlV1H&index=6)
* [LA R Users Meetup, October 2020 (1:14:40)](https://www.youtube.com/watch?v=Qq25BUxpJu4)
* [New York Open Statistical Programming Meetup, December 2020 (1:54:28)](https://youtu.be/Gqn7Xn4d5NI)
* [ds-incubator series, 2021](https://www.youtube.com/playlist?list=PLvgdJdJDL-APJqHy5CXs6m4N7hUVp5rb4)
* [Lille R User Group, June 2021 (45:54)](https://youtu.be/FODSavXGjYg)

### Español

* [R-Ladies Barcelona, 2021-05-25 (1:25:12)](https://www.youtube.com/watch?v=Vj312AfdpBo).

## Documentation

* [User manual](https://books.ropensci.org/targets/): in-depth discussion about how to use `targets`.
* [Reference website](https://docs.ropensci.org/targets/): formal documentation of all user-side functions, the statement of need, and multiple design documents of the internal architecture.
* [Developer documentation](https://books.ropensci.org/targets-design/): software design documents for developers contributing to the deep internal architecture of `targets`.

## Courses

* [Official half-day interactive tutorial](https://github.com/wlandau/targets-tutorial).

## Example projects

Description | Link
---|---
Minimal example | <https://github.com/wlandau/targets-minimal>
Machine learning with Keras | <https://github.com/wlandau/targets-keras>
Validating a minimal Stan model | <https://github.com/wlandau/targets-stan>
Using Target Markdown and `stantargets` to validate a Bayesian longitudinal model for clinical trial data analysis | <https://github.com/wlandau/rmedicine2021-pipeline>
Shiny app that runs a pipeline | <https://github.com/wlandau/targets-shiny>
Deploying a pipeline to RStudio Connect | <https://github.com/sol-eng/targets-deployment-rsc>

## Apps

* [`tar_watch()`](https://docs.ropensci.org/targets/reference/tar_watch.html): a built-in Shiny app to visualize progress while a pipeline is running. Available as a Shiny module via [`tar_watch_ui()`](https://docs.ropensci.org/targets/reference/tar_watch_ui.html) and [`tar_watch_server()`](https://docs.ropensci.org/targets/reference/tar_watch_server.html).
* [`targetsketch`](https://wlandau.shinyapps.io/targetsketch): a Shiny app to help sketch pipelines ([app](https://wlandau.shinyapps.io/targetsketch), [source](https://github.com/wlandau/targetsketch)).

## Deployment

* <https://solutions.rstudio.com/r/workflows/> explains how to deploy a pipeline to RStudio Connect ([example code](https://github.com/sol-eng/targets-deployment-rsc)).
* [`tar_github_actions()`](https://docs.ropensci.org/targets/reference/tar_github_actions.html) sets up a pipeline to run on GitHub Actions. The [minimal example](https://github.com/wlandau/targets-minimal) demonstrates this approach.

## Extending and customizing targets

* [R Targetopia](https://wlandau.github.io/targetopia/): a collection of [R packages](https://wlandau.github.io/targetopia/packages.html) that extend `targets`. [These packages](https://wlandau.github.io/targetopia/packages.html) simplify pipeline construction for specific fields of Statistics and data science.
* [Target factories](https://wlandau.github.io/targetopia/contributing.html#target-factories): a programming technique to write specialized interfaces for custom pipelines. Posts [here](https://ropensci.org/blog/2021/02/03/targets/) and [here](https://wlandau.github.io/targetopia/contributing.html) describe how.

## Help

* Post to the [GitHub discussion forum](https://github.com/ropensci/targets/discussions) to ask questions. To get the best help about a specific issue, create a reproducible example with  [`targets::tar_reprex()`](https://docs.ropensci.org/targets/reference/tar_reprex.html) or  [`reprex::reprex()`](https://reprex.tidyverse.org/reference/reprex.html). 
* The [RStudio Community](https://community.rstudio.com/) forum is full of friendly enthusiasts of R and the tidyverse. Use the [`targets` tag](https://community.rstudio.com/tag/targets).
* [Stack Overflow](https://stackoverflow.com/) broadcasts to the entire open source community. Use the [`targets-r-package` tag](https://stackoverflow.com/questions/tagged/targets-r-package).

## Code of conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/).

## Citation

```{r}
citation("targets")
```
---
title: "Target Markdown: non-defaults, relaying, and edge cases"
output: html_document
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

This report is like `test-target_markdown_default.Rmd` except it uses non-default scripts and data. Try some code chunks interactively, then render this report. Render again to see if the targets correctly skip. Make sure `_targets.R`, `_targets_r/`, and `_targets/` do not exist. Look for the scripts and data are all in the `example/` folder. Delete the HTML file and `example/` when you are done. 

# Packages

The example requires several R packages, and `targets` must be version 0.5.0.9000 or above. 

```{r, eval = FALSE}
install.packages(c("biglm", "dplyr", "ggplot2", "readr", "targets", "tidyr"))
```

# Setup

First, load `targets` to activate the specialized `knitr` engine for Target Markdown.

```{r}
library(targets)
knitr::opts_chunk$set(tar_script = "example/script.R")
```

Near the top, you may also wish to remove the `_targets_r` directory previously written by non-interactive runs of the report.

```{r}
tar_unscript()
```

# Globals

We first define some global options/functions common to all targets. The function below plots a histogram of ozone concentrations, and our histogram target will need it.

```{targets example-globals, tar_globals = TRUE}
options(tidyverse.quiet = TRUE)
tar_option_set(
  packages = c("biglm", "dplyr", "ggplot2", "readr", "tidyr"),
  error = "continue"
)
create_plot <- function(data) {
  ggplot(data) +
    geom_histogram(aes(x = Ozone), bins = 12) +
    theme_gray(24)
}
```

# Targets

Our first target borrows the `airquality` dataset built into base R.

```{targets othername, tar_name = "raw_data", tar_simple = TRUE}
airquality
```

Warnings, messages, and errors are relayed in both modes.

```{targets raw_data_interactive, tar_simple = TRUE, error = TRUE, message = TRUE, output = TRUE, tar_interactive = TRUE}
print("output_interactive")
message("message_interactive")
warning("warning_interactive")
stop("error_interactive")
airquality
```

```{targets raw_data_noninteractive, tar_simple = TRUE, error = TRUE, message = TRUE, output = TRUE, tar_interactive = FALSE}
print("output_noninteractive")
message("message_noninteractive")
warning("warning_noninteractive")
stop("error_noninteractive")
airquality
```

Our next targets preprocess the data, make a histogram, and fit a model.

```{targets downstream-targets}
list(
  tar_target(data, raw_data %>% filter(!is.na(Ozone))),
  tar_target(hist, create_plot(data)),
  tar_target(fit, biglm(Ozone ~ Wind + Temp, data))
)
```

# Pipeline

If you ran all the `{targets}` chunks in non-interactive mode, then your R scripts are set up to run the pipeline.

```{r}
tar_make(script = "example/script.R", store = "example/store")
```

# Output

You can retrieve results from the `_targets/` data store using `tar_read()` or `tar_load()`.

```{r, message = FALSE}
library(biglm)
tar_read(fit, store = "example/store")
```

```{r}
tar_read(hist, store = "example/store")
```

The `targets` dependency graph helps your readers understand the steps of your pipeline at a high level.

```{r}
tar_visnetwork(script = "example/script.R", store = "example/store")
```

At this point, you can go back and run `{targets}` chunks in interactive mode without interfering with the code or data of the non-interactive pipeline.
---
title: "Target Markdown with defaults"
output: html_document
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

Try some code chunks interactively, then render this report. Render again to see if the targets correctly skip. Delete the HTML file, `_targets.R`, `_targets_r/`, and `_targets/` after rendering.

# Packages

The example requires several R packages, and `targets` must be version 0.5.0.9000 or above. 

```{r, eval = FALSE}
install.packages(c("biglm", "dplyr", "ggplot2", "readr", "targets", "tidyr"))
```

# Setup

First, load `targets` to activate the specialized `knitr` engine for Target Markdown.

```{r}
library(targets)
```

Near the top, you may also wish to remove the `_targets_r` directory previously written by non-interactive runs of the report.

```{r}
tar_unscript()
```

# Globals

We first define some global options/functions common to all targets. The function below plots a histogram of ozone concentrations, and our histogram target will need it.

```{targets example-globals, tar_globals = TRUE}
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("biglm", "dplyr", "ggplot2", "readr", "tidyr"))
create_plot <- function(data) {
  ggplot(data) +
    geom_histogram(aes(x = Ozone), bins = 12) +
    theme_gray(24)
}
```

# Targets

Our first target borrows the `airquality` dataset built into base R.

```{targets raw-data}
tar_target(raw_data, airquality)
```

Our next targets preprocess the data, make a histogram, and fit a model.

```{targets downstream-targets}
list(
  tar_target(data, raw_data %>% filter(!is.na(Ozone))),
  tar_target(hist, create_plot(data)),
  tar_target(fit, biglm(Ozone ~ Wind + Temp, data))
)
```

# Pipeline

If you ran all the `{targets}` chunks in non-interactive mode, then your R scripts are set up to run the pipeline.

```{r}
tar_make()
```

# Output

You can retrieve results from the `_targets/` data store using `tar_read()` or `tar_load()`.

```{r, message = FALSE}
library(biglm)
tar_read(fit)
```

```{r}
tar_read(hist)
```

The `targets` dependency graph helps your readers understand the steps of your pipeline at a high level.

```{r}
tar_visnetwork()
```

At this point, you can go back and run `{targets}` chunks in interactive mode without interfering with the code or data of the non-interactive pipeline.
---
title: "Target Markdown"
output: html_document
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

Target Markdown is a powerful R Markdown interface for reproducible analysis pipelines, and the chapter at https://books.ropensci.org/targets/markdown.html walks through it in detail. This R Markdown report the example from the chapter. Try it out in both interactive and non-interactive modes, either by running the code chunks in different ways or setting the `tar_interactive` chunk option.

# Packages

The example requires several R packages, and `targets` must be version 0.5.0.9000 or above. 

```{r, eval = FALSE}
install.packages(c("biglm", "dplyr", "ggplot2", "readr", "targets", "tidyr"))
```

# Setup

If you are using old versions of `targets` (<= 0.7.0) and/or `knitr` (<= 1.33), you will need to load the `targets` package in the R Markdown document in order for Target Markdown code chunks to work.

```{r}
library(targets)
```

Near the top of the document, you may also wish to remove the `_targets_r` directory previously written by non-interactive runs of the report. Otherwise, your pipeline may contain superfluous targets.

```{r}
library(targets)
tar_unscript()
```

# Globals

We first define some global options/functions common to all targets. The function below plots a histogram of ozone concentrations, and our histogram target will need it.

```{targets example-globals, tar_globals = TRUE}
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("biglm", "dplyr", "ggplot2", "readr", "tidyr"))
create_plot <- function(data) {
  ggplot(data) +
    geom_histogram(aes(x = Ozone), bins = 12) +
    theme_gray(24)
}
```

# Targets

Our first target borrows the `airquality` dataset built into base R.

```{targets raw-data}
tar_target(raw_data, airquality)
```

Our next targets preprocess the data, make a histogram, and fit a model.

```{targets downstream-targets}
list(
  tar_target(data, raw_data %>% filter(!is.na(Ozone))),
  tar_target(hist, create_plot(data))
)
```

Set the `tar_simple` chunk option to `TRUE` to define a single target with the command in the code chunk. The chunk below only contains `biglm(Ozone ~ Wind + Temp, data)` in the source, but because `tar_simple` is `TRUE`, it is shorthand for `tar_target(name = fit, command = biglm(Ozone ~ Wind + Temp, data))`. All other arguments to `tar_target()` are set to their default values (configurable with `tar_option_set()`).

```{targets fit, tar_simple = TRUE}
biglm(Ozone ~ Wind + Temp, data)
```

# Pipeline

If you ran all the `{targets}` chunks in non-interactive mode, then your R scripts are set up to run the pipeline.

```{r}
tar_make()
```

# Output

You can retrieve results from the `_targets/` data store using `tar_read()` or `tar_load()`.

```{r, message = FALSE}
library(biglm)
tar_read(fit)
```

```{r}
tar_read(hist)
```

The `targets` dependency graph helps your readers understand the steps of your pipeline at a high level.

```{r}
tar_visnetwork()
```

At this point, you can go back and run `{targets}` chunks in interactive mode without interfering with the code or data of the non-interactive pipeline.
---
title: "An overview of targets"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{An overview of targets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This vignette is a high-level overview of `targets` and its educational materials. The goal is to summarize the major features of `targets` and direct users to the appropriate resources. It explains how to get started, and then it briefly describes each chapter of the [user manual](https://books.ropensci.org/targets/).

## What is `targets`?

The `targets` R package is a Make-like pipeline toolkit for Statistics and data science in R. `targets` accelerates analysis with easy-to-configure parallel computing, enhances reproducibility, and reduces the burdens of repeated computation and manual data micromanagement. A fully up-to-date `targets` pipeline is tangible evidence that the output aligns with the code and data, which substantiates trust in the results.

## How to get started

The top of the [reference website](https://docs.ropensci.org/targets/) links to a number of materials to help new users start learning `targets`. It lists online talks, tutorials, books, and workshops in the order that a new user should consume them. The rest of the main page outlines a more comprehensive list of resources.

## The walkthrough

The [user manual](https://books.ropensci.org/targets/) starts with a  [walkthrough](https://books.ropensci.org/targets/walkthrough.html) chapter, a short tutorial to quickly started with `targets` using a simple example project. That project also has a [repository](https://github.com/wlandau/targets-minimal) with the source code and an [RStudio Cloud workspace](https://rstudio.cloud/project/1430691) that lets you try out the workflow in a web browser. Sign up for a free RStudio Cloud account, click on the link, and try out functions `tar_make()` and `tar_read()` in the R console.

## Target Markdown

[Target Markdown](https://books.ropensci.org/targets/markdown.html) is an R Markdown interface for testing, prototyping, and constructing targets and pipelines. It combines the convenience of R Markdown with the power of `targets`. See the [chapter](https://books.ropensci.org/targets/markdown.html) in the manual for a tutorial.

## Debugging

The [debugging chapter](https://books.ropensci.org/targets/debugging.html) describes two alternative built-in systems for troubleshooting errors. The first system uses workspaces, which let you load a target's dependencies into you R session. This way is usually preferred, especially with large pipelines on computing clusters, but it still may require some manual work. The second system launches an interactive debugger while the pipeline is actually running, which may not be feasible in some situations, but can often help you reach the problem more quickly.

## Functions

`targets` expects users to adopt a function-oriented style of programming. User-defined R functions are essential to express the complexities of data generation, analysis, and reporting. The [user manual](https://books.ropensci.org/targets/) has a [whole chapter](https://books.ropensci.org/targets/functions.html) dedicated to user-defined functions for data science, and it explains why they are important and how to use them in `targets`-powered pipelines.

## Target construction

The [target construction chapter](https://books.ropensci.org/targets/targets.html) explains best practices for creating targets: what a good target should do, how much work a target should do, and guidelines for thinking about side effects and upstream dependencies (i.e. other targets and global objects).

## Packages

The [packages chapter](https://books.ropensci.org/targets/packages.html) explains best practices for working with packages in `targets`: how to load them, how to work with packages as projects, target factories inside packages, and automatically invalidating targets based on changes inside one or more packages.

## Projects

The [projects chapter](https://books.ropensci.org/targets/projects.html) explains best practices for working with `targets`-powered projects: the recommended file structure, recommended third-party tools, multi-project repositories, and interdependent projects.

## External files and literate programming

`targets` has special ways to include data files and literate programming reports in a pipeline. This functionality is optional in the general case, but it is necessary if you want a target to rerun in response to a change in a data file, or if you want an R Markdown report to re-render when an upstream target changes. The [files chapter](https://books.ropensci.org/targets/files.html) walks through this functionality, from input data to parameterized R Markdown.

## Dynamic branching

Sometimes, a pipeline contains more targets than a user can comfortably type by hand. For projects with hundreds of targets, branching can make the _targets.R file more concise and easier to read and maintain. Dynamic branching is a way to create new targets while the pipeline is running, and it is best suited to iterating over a larger number of very similar tasks. The [dynamic branching](https://books.ropensci.org/targets/dynamic.html) chapter outlines this functionality, including how to create branching patterns, different ways to iterate over data, and recommendations for batching large numbers of small tasks into a comfortably small number of dynamic branches.

## Static branching

[Static branching](https://books.ropensci.org/targets/static.html) is the act of defining a group of targets in bulk before the pipeline starts. Whereas dynamic branching uses last-minute dependency data to define the branches, static branching uses metaprogramming to modify the code of the pipeline up front. Whereas dynamic branching excels at creating a large number of very similar targets, static branching is most useful for smaller number of heterogeneous targets. Some users find it more convenient because they can use `tar_manifest()` and `tar_visnetwork()` to check the correctness of static branching before launching the pipeline. Read more about it in the [static branching chapter](https://books.ropensci.org/targets/static.html).

## High-performance computing

`targets` is capable of distributing the computation in a pipeline across multiple cores of a laptop or multiple nodes of a computing cluster. Not only does it interface with these technologies using packages [`clustermq`](https://github.com/mschubert/clustermq) and [`future`](https://github.com/HenrikBengtsson/future): it automatically deploys ready targets to parallel workers while making sure the other targets wait for their upstream dependencies to finish. Read more about high-performance computing in the [HPC chapter](https://books.ropensci.org/targets/hpc.html).

## Cloud computing

Users with [Amazon Web Services](https://aws.amazon.com/) accounts can store their targets on one or more [S3 buckets](https://aws.amazon.com/s3/), and retrieval with `tar_read()` and `tar_load()` is seamless. The [AWS storage chapter](https://books.ropensci.org/targets/storage_amazon.html) is a step-by-step guide that walks through how to get started with [Amazon Web Services](https://aws.amazon.com/) and connect a `targets` pipeline to [S3](https://aws.amazon.com/s3/).

## What about drake?

The [`drake`](https://github.com/ropensci/drake) package is an older and more established R-focused pipeline toolkit, and it is the predecessor of `targets`. The [`drake` chapter](https://books.ropensci.org/targets/drake.html) of the `targets` manual helps `drake` users understand the role of `targets`, the future direction of `drake`, how to transition to `targets`, and the advantages of `targets` over `drake`.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_make_future.R
\name{tar_make_future}
\alias{tar_make_future}
\title{Run a pipeline of targets in parallel with transient
\code{future} workers.}
\usage{
tar_make_future(
  names = NULL,
  shortcut = targets::tar_config_get("shortcut"),
  reporter = targets::tar_config_get("reporter_make"),
  workers = targets::tar_config_get("workers"),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function, reporter),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Names of the targets to build or check. Set to \code{NULL} to
check/build all the targets (default). Otherwise, you can supply
\code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
Because \code{\link[=tar_make]{tar_make()}} and friends run the pipeline in a new R session,
if you pass a character vector to a tidyselect helper, you will need
to evaluate that character vector early with \verb{!!}, e.g.
\code{tar_make(names = all_of(!!your_vector))}.
Applies to ordinary targets (stem) and whole dynamic branching targets
(patterns) but not to individual dynamic branches.}

\item{shortcut}{Logical of length 1, how to interpret the \code{names} argument.
If \code{shortcut} is \code{FALSE} (default) then the function checks
all targets upstream of \code{names} as far back as the dependency graph goes.
\code{shortcut = TRUE} increases speed if there are a lot of
up-to-date targets, but it assumes all the dependencies
are up to date, so please use with caution.
It relies on stored metadata for information about upstream dependencies.
\code{shortcut = TRUE} only works if you set \code{names}.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets run in the pipeline.
Defaults to \code{tar_config_get("reporter_make")}. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"summary"}: print a running total of the number of each targets in
each status category (queued, started, skipped, build, canceled,
or errored). Also show a timestamp (\code{"\%H:\%M \%OS2"} \code{strptime()} format)
of the last time the progress changed and printed to the screen.
\item \code{"timestamp"}: same as the \code{"verbose"} reporter except that each
.message begins with a time stamp.
\item \code{"timestamp_positives"}: same as the \code{"timestamp"} reporter
except without messages for skipped targets.
\item \code{"verbose"}: print messages for individual targets
as they start, finish, or are skipped.
\item \code{"verbose_positives"}: same as the \code{"verbose"} reporter
except without messages for skipped targets.
}}

\item{workers}{Positive integer, maximum number of transient
\code{future} workers allowed to run at any given time.}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
\code{NULL} except if \code{callr_function = callr::r_bg()}, in which case
a handle to the \code{callr} background process is returned. Either way,
the value is invisibly returned.
}
\description{
This function is like \code{\link[=tar_make]{tar_make()}} except that targets
run in parallel with transient \code{future} workers. It requires
that you declare your \code{future::plan()} inside the
target script file (default: \verb{_targets.R}).
\code{future} is not a strict dependency of \code{targets},
so you must install \code{future} yourself.
}
\details{
To configure \code{tar_make_future()} with a computing cluster,
see the \code{future.batchtools} package documentation.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  future::plan(future::multisession, workers = 2)
  list(
    tar_target(x, 1 + 1),
    tar_target(y, 1 + 1)
  )
}, ask = FALSE)
tar_make_future()
})
}
}
\seealso{
Other pipeline: 
\code{\link{tar_make_clustermq}()},
\code{\link{tar_make}()}
}
\concept{pipeline}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_branches.R
\name{tar_branches}
\alias{tar_branches}
\title{Reconstruct the branch names and the names of their dependencies.}
\usage{
tar_branches(name, pattern, store = targets::tar_config_get("store"))
}
\arguments{
\item{name}{Symbol, name of the target.}

\item{pattern}{Language to define branching for a target.
For example, in a pipeline with numeric vector targets \code{x} and \code{y},
\code{tar_target(z, x + y, pattern = map(x, y))} implicitly defines
branches of \code{z} that each compute \code{x[1] + y[1]}, \code{x[2] + y[2]},
and so on. See the user manual for details.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A \code{tibble} with one row per branch and one column for each target
(including the branched-over targets and the target with the pattern.)
}
\description{
Given a branching pattern, use available metadata
to reconstruct branch names and the names of each
branch's dependencies. The metadata of each target
must already exist and be consistent with the metadata
of the other targets involved.
}
\details{
The results from this function can help you retroactively
figure out correspondences between upstream branches and downstream
branches. However, it does not always correctly predict what the
names of the branches will be after the next run of the pipeline.
Dynamic branching happens while the pipeline is running,
so we cannot always know what the names of the branches will be
in advance (or even how many there will be).
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, head(letters, 2)),
    tar_target(z, head(LETTERS, 2)),
    tar_target(dynamic, c(x, y, z), pattern = cross(z, map(x, y)))
  )
}, ask = FALSE)
tar_make()
tar_branches(dynamic, pattern = cross(z, map(x, y)))
})
}
}
\seealso{
Other branching: 
\code{\link{tar_branch_index}()},
\code{\link{tar_branch_names_raw}()},
\code{\link{tar_branch_names}()},
\code{\link{tar_pattern}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_validate.R
\name{tar_validate}
\alias{tar_validate}
\title{Validate a pipeline of targets.}
\usage{
tar_validate(
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
\code{NULL} except if \code{callr_function = callr::r_bg()}, in which case
a handle to the \code{callr} background process is returned. Either way,
the value is invisibly returned.
}
\description{
Inspect the pipeline for issues and throw an error or
warning if a problem is detected.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(list(tar_target(x, 1 + 1)), ask = FALSE)
tar_validate()
})
}
}
\seealso{
Other inspect: 
\code{\link{tar_deps_raw}()},
\code{\link{tar_deps}()},
\code{\link{tar_glimpse}()},
\code{\link{tar_manifest}()},
\code{\link{tar_network}()},
\code{\link{tar_outdated}()},
\code{\link{tar_sitrep}()},
\code{\link{tar_visnetwork}()}
}
\concept{inspect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_process.R
\name{tar_process}
\alias{tar_process}
\title{Get main process info.}
\usage{
tar_process(names = NULL, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Optional, names of the data points to return.
If supplied, \code{tar_process()}
returns only the rows of the names you select.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
If \code{NULL}, all names are selected.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A data frame with metadata on the most recent main R process
to orchestrate the targets of the current project.
The output includes the \code{pid} of the main process.
}
\description{
Get info on the most recent main R process
to orchestrate the targets of the current project.
}
\details{
The main process is the R process invoked
by \code{\link[=tar_make]{tar_make()}} or similar. If \code{callr_function} is not \code{NULL},
this is an external process, and the \code{pid} in the return value
will not agree with \code{Sys.getpid()} in your current interactive session.
The process may or may not be alive. You may want to
check the status with \code{tar_pid() \%in\% ps::ps_pids()}
before running another call to \code{\link[=tar_make]{tar_make()}}
for the same project.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
tar_process()
tar_process(pid)
})
}
}
\seealso{
Other data: 
\code{\link{tar_load_raw}()},
\code{\link{tar_load}()},
\code{\link{tar_meta}()},
\code{\link{tar_objects}()},
\code{\link{tar_pid}()},
\code{\link{tar_read_raw}()},
\code{\link{tar_read}()}
}
\concept{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_invalidate.R
\name{tar_invalidate}
\alias{tar_invalidate}
\title{Invalidate targets and global objects in the metadata.}
\usage{
tar_invalidate(names, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Names of the targets to remove from the metadata list.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\description{
Delete the metadata of records in \verb{_targets/meta/meta}
but keep the return values of targets in \verb{_targets/objects/}.
}
\details{
This function forces one or more targets to rerun
on the next \code{\link[=tar_make]{tar_make()}}, regardless of the cues and regardless
of how those targets are stored. After \code{tar_invalidate()},
you will still be able to locate the data files with \code{\link[=tar_path]{tar_path()}}
and manually salvage them in an emergency.
However, \code{\link[=tar_load]{tar_load()}} and \code{\link[=tar_read]{tar_read()}} will not be able to
read the data into R, and subsequent calls to \code{\link[=tar_make]{tar_make()}}
will attempt to rerun those targets.
For patterns recorded in the metadata, all the branches
will be invalidated. For patterns no longer in the metadata,
branches are left alone.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
}, ask = FALSE)
tar_make()
tar_invalidate(starts_with("y")) # Only invalidates y1 and y2.
tar_make() # y1 and y2 rerun but return same values, so z is up to date.
})
}
}
\seealso{
Other clean: 
\code{\link{tar_delete}()},
\code{\link{tar_destroy}()},
\code{\link{tar_prune}()}
}
\concept{clean}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_exist_process.R
\name{tar_exist_process}
\alias{tar_exist_process}
\title{Check if process metadata exists.}
\usage{
tar_exist_process(store = targets::tar_config_get("store"))
}
\arguments{
\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Logical of length 1, whether the current project's metadata exists.
}
\description{
Check if the process metadata file \verb{_targets/meta/process}
exists for the current project.
}
\details{
To learn more about local storage in \code{targets}, visit
\url{https://books.ropensci.org/targets/files.html#internal-files}.
}
\examples{
tar_exist_process()
}
\seealso{
Other existence: 
\code{\link{tar_exist_meta}()},
\code{\link{tar_exist_objects}()},
\code{\link{tar_exist_progress}()},
\code{\link{tar_exist_script}()}
}
\concept{existence}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_unscript.R
\name{tar_unscript}
\alias{tar_unscript}
\title{Remove target script helper files.}
\usage{
tar_unscript(script = targets::tar_config_get("script"))
}
\arguments{
\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}
}
\value{
\code{NULL} (invisibly).
}
\description{
Remove target script helper files (default: \verb{_targets_r/})
that were created by Target Markdown.
}
\details{
Target Markdown code chunks create R scripts in a folder
called \verb{_targets_r/} in order to aid the automatically supplied
\verb{_targets.R} file. Over time, the number of script files
starts to build up, and \code{targets} has no way of automatically
removing helper script files that are no longer necessary.
To keep your pipeline up to date
with the code chunks in the Target Markdown document(s),
it is good practice to call \code{tar_unscript()} at the beginning
of your first Target Markdown document. That way,
extraneous/discarded targets are automatically
removed from the pipeline when the document starts render.

If the target script is at some alternative path,
e.g. \code{custom/script.R}, the helper scripts are in \verb{custom/script_r/}.
\code{\link[=tar_unscript]{tar_unscript()}} works on the helper scripts as long as your
project configuration settings correctly identify the correct
target script.
}
\examples{
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_unscript()
})
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_read_raw.R
\name{tar_read_raw}
\alias{tar_read_raw}
\title{Read a target's value from storage (raw version)}
\usage{
tar_read_raw(
  name,
  branches = NULL,
  meta = tar_meta(store = store),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{name}{Character, name of the target to read.}

\item{branches}{Integer of indices of the branches to load
if the target is a pattern.}

\item{meta}{Data frame of metadata from \code{\link[=tar_meta]{tar_meta()}}.
\code{tar_read()} with the default arguments can be inefficient for large
pipelines because all the metadata is stored in a single file.
However, if you call \code{\link[=tar_meta]{tar_meta()}} beforehand and supply it to the \code{meta}
argument, then successive calls to \code{tar_read()} may run much faster.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
The target's return value from its file in
\verb{_targets/objects/}, or the paths to the custom files and directories
if \code{format = "file"} was set.
}
\description{
Like \code{\link[=tar_read]{tar_read()}} except \code{name} is a character string.
Do not use in \code{knitr} or R Markdown reports with \code{tarchetypes::tar_knit()}
or \code{tarchetypes::tar_render()}.
}
\section{Limited scope}{

\code{tar_read()} and \code{tar_load()}
are only for exploratory analysis and literate programming,
and \code{tar_read_raw()} and \code{tar_load_raw()} are only
for exploratory analysis. \code{targets} automatically
loads the correct dependencies into memory when the pipeline
is running, so invoking these functions
from inside a target is rarely advisable.
}

\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(list(tar_target(x, 1 + 1)), ask = FALSE)
tar_make()
tar_read_raw("x")
})
}
}
\seealso{
Other data: 
\code{\link{tar_load_raw}()},
\code{\link{tar_load}()},
\code{\link{tar_meta}()},
\code{\link{tar_objects}()},
\code{\link{tar_pid}()},
\code{\link{tar_process}()},
\code{\link{tar_read}()}
}
\concept{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio_addin_tar_read.R
\name{rstudio_addin_tar_read}
\alias{rstudio_addin_tar_read}
\title{RStudio addin to call \code{\link[=tar_read]{tar_read()}} on the symbol at the cursor.}
\usage{
rstudio_addin_tar_read(context = NULL)
}
\arguments{
\item{context}{RStudio API context from
\code{rstudioapi::getActiveDocumentContext()}.}
}
\description{
For internal use only. Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_exist_objects.R
\name{tar_exist_objects}
\alias{tar_exist_objects}
\title{Check if local output data exists for one or more targets.}
\usage{
tar_exist_objects(names, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Character vector of target names.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Logical of length \code{length(names)}, whether
each given target has an existing file in \verb{_targets/objects/}
for the current project.
}
\description{
Check if the local data files exist in
\verb{_targets/objects/} for one or more targets.
}
\details{
To learn more about local storage in \code{targets}, visit
\url{https://books.ropensci.org/targets/files.html#internal-files}.
}
\examples{
tar_exist_objects(c("target1", "target2"))
}
\seealso{
Other existence: 
\code{\link{tar_exist_meta}()},
\code{\link{tar_exist_process}()},
\code{\link{tar_exist_progress}()},
\code{\link{tar_exist_script}()}
}
\concept{existence}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_renv.R
\name{tar_renv}
\alias{tar_renv}
\title{Set up package dependencies for compatibility with \code{renv}}
\usage{
tar_renv(
  extras = c("bs4Dash", "clustermq", "future", "gt", "markdown", "pingr", "rstudioapi",
    "shiny", "shinybusy", "shinyWidgets", "visNetwork"),
  path = "_targets_packages.R",
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function),
  envir = parent.frame(),
  script = targets::tar_config_get("script")
)
}
\arguments{
\item{extras}{Character vector of additional packages to declare as
project dependencies.}

\item{path}{Character of length 1, path to the script file to
populate with \code{library()} calls.}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}
}
\value{
Nothing, invisibly.
}
\description{
Write package dependencies to a script file
(by default, named \verb{_targets_packages.R} in the root project directory).
Each package is written to a separate line
as a standard \code{\link[=library]{library()}} call (e.g. \code{library(package)}) so
\code{renv} can identify them automatically.
}
\details{
This function gets called for its side-effect, which writes
package dependencies to a script for compatibility with \code{renv}.
The generated file should \strong{not} be edited by hand and will be
overwritten each time \code{tar_renv()} is called.

The behavior of \code{renv} is to create and manage a project-local \code{R} library
and keep a record of project dependencies in a file called \code{renv.lock}.
To identify dependencies, \code{renv} crawls through code to find packages
explicitly mentioned using \code{library()}, \code{require()}, or \code{::}.
However, \code{targets} manages packages in a way that hides dependencies
from \code{renv.} \code{tar_renv()} finds package dependencies that would be
otherwise hidden to \code{renv} because they are declared using the \code{targets}
API. Thus, calling \code{tar_renv} this is only necessary if using
\code{\link[=tar_option_set]{tar_option_set()}} or \code{\link[=tar_target]{tar_target()}} to use specialized storage
formats or manage packages.

With the script written by \code{tar_renv()}, \code{renv} is able to crawl the
file to identify package dependencies (with \code{renv::dependencies()}).
\code{tar_renv()} only serves to make your \code{targets} project compatible with
\code{renv}, it is still the users responsibility to call \code{renv::init()} and
\code{renv::snapshot()} directly to initialize and manage a
project-local \code{R} library. This allows your \code{targets} pipeline to have
its own self-contained \code{R} library separate from your standard \code{R}
library. See \url{https://rstudio.github.io/renv/index.html} for
more information.
}
\examples{
tar_dir({ # tar_dir() runs code from a temporary directory.
  tar_script({
    tar_option_set(packages = c("tibble", "qs"))
    list()
  }, ask = FALSE)
  tar_renv()
  writeLines(readLines("_targets_packages.R"))
})
tar_option_reset()
}
\seealso{
\url{https://rstudio.github.io/renv/articles/renv.html}

Other scripts: 
\code{\link{tar_edit}()},
\code{\link{tar_github_actions}()},
\code{\link{tar_helper_raw}()},
\code{\link{tar_helper}()},
\code{\link{tar_script}()}
}
\concept{scripts}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_timestamp_raw.R
\name{tar_timestamp_raw}
\alias{tar_timestamp_raw}
\title{Get the timestamp(s) of a target (raw version).}
\usage{
tar_timestamp_raw(
  name = NULL,
  format = NULL,
  tz = NULL,
  parse = NULL,
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{name}{Character of length 1, name of the target.}

\item{format}{Deprecated in \code{targets} version 0.6.0 (2021-07-21).}

\item{tz}{Deprecated in \code{targets} version 0.6.0 (2021-07-21).}

\item{parse}{Deprecated in \code{targets} version 0.6.0 (2021-07-21).}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
If the target is not recorded in the metadata
or cannot be parsed correctly, then
\code{tar_timestamp_raw()} returns a \code{POSIXct} object at \verb{1970-01-01 UTC}.
}
\description{
Get the time that a target last ran successfully.
}
\details{
\code{tar_timestamp_raw()} is like \code{tar_timestamp()} except
it accepts a character string for \code{name} instead of a symbol.
\code{tar_timestamp_raw()} checks the metadata in \verb{_targets/meta/meta},
not the actual data. Time stamps are recorded only for targets that
run commands: just non-branching targets and individual dynamic
branches.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(tar_target(x, 1))
}, ask = FALSE)
tar_make()
# Get the timestamp.
tar_timestamp_raw("x")
# We can use the timestamp to cancel the target
# if it already ran within the last hour.
# Be sure to set `cue = tar_cue(mode = "always")`
# if you want the target to always check the timestamp.
tar_script({
  list(
  tar_target(
    x,
    tar_cancel((Sys.time() - tar_timestamp_raw()) < 3600),
    cue = tar_cue(mode = "always")
  )
)}, ask = FALSE)
tar_make()
})
}
}
\seealso{
Other time: 
\code{\link{tar_newer}()},
\code{\link{tar_older}()},
\code{\link{tar_timestamp}()}
}
\concept{time}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_pipeline.R
\name{pipeline_validate_lite}
\alias{pipeline_validate_lite}
\title{Abridged pipeline validation function.}
\usage{
pipeline_validate_lite(pipeline)
}
\arguments{
\item{pipeline}{A pipeline object.}
}
\description{
Internal function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_toggle.R
\name{tar_toggle}
\alias{tar_toggle}
\title{Choose code to run based on Target Markdown mode.}
\usage{
tar_toggle(interactive, noninteractive)
}
\arguments{
\item{interactive}{R code to run if Target Markdown interactive mode is
activated.}

\item{noninteractive}{R code to run if Target Markdown interactive mode is
not activated.}
}
\value{
If Target Markdown interactive mode is not turned on,
the function returns the result of running the code.
Otherwise, the function invisibly returns \code{NULL}.
}
\description{
Run one piece of code if Target Markdown mode
interactive mode is turned on and another piece of code otherwise.
}
\details{
Visit <books.ropensci.org/targets/markdown.html>
to learn about Target Markdown and interactive mode.
}
\examples{
tar_toggle(
  message("In interactive mode."),
  message("Not in interactive mode.")
)
}
\seealso{
Other Target Markdown: 
\code{\link{tar_engine_knitr}()},
\code{\link{tar_interactive}()},
\code{\link{tar_noninteractive}()}
}
\concept{Target Markdown}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_config_get.R
\name{tar_config_get}
\alias{tar_config_get}
\title{Get configuration settings.}
\usage{
tar_config_get(
  name,
  config = Sys.getenv("TAR_CONFIG", "_targets.yaml"),
  project = Sys.getenv("TAR_PROJECT", "main")
)
}
\arguments{
\item{name}{Character of length 1, name of the specific
configuration setting to retrieve.}

\item{config}{Character of length 1, file path of the YAML
configuration file with \code{targets} project settings.
The \code{config} argument specifies which YAML configuration
file that \code{tar_config_get()} reads from or \code{tar_config_set()}
writes to in a single function call.
It does not globally change which configuration file is used
in subsequent function calls. The default file path of the YAML
file is always \verb{_targets.yaml} unless you set another
default path using the \code{TAR_CONFIG} environment variable,
e.g. \code{Sys.setenv(TAR_CONFIG = "custom.yaml")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}

\item{project}{Character of length 1, name of the current
\code{targets} project. Thanks to the \code{config} R package,
\code{targets} YAML configuration files can store multiple
sets of configuration settings, with each set corresponding
to its own project. The \code{project} argument allows you to
set or get a configuration setting for a specific project
for a given call to \code{tar_config_set()} or \code{tar_config_get()}.
The default project is always called \code{"main"}
unless you set another
default project using the \code{TAR_PROJECT} environment variable,
e.g. \code{Sys.setenv(tar_project = "custom")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}
}
\value{
The value of the configuration setting from
the YAML configuration file (default: \verb{_targets.yaml})
or the default value if the setting is not available.
The data type of the return value depends on your choice
of \code{name}.
}
\description{
Read the custom settings for the current project
in the optional YAML configuration file.
}
\section{Configuration}{

For several key functions like \code{\link[=tar_make]{tar_make()}}, the
default values of arguments are controlled though
\code{tar_config_get()}. \code{tar_config_get()} retrieves data
from an optional YAML configuration file.
You can control the settings in the YAML
file programmatically with \code{tar_config_set()}.
The default file path of this YAML file is \verb{_targets.yaml}, and you can
set another path globally using the \code{TAR_CONFIG}
environment variable. The YAML file can store configuration
settings for multiple projects, and you can globally
set the default project with the \code{TAR_PROJECT} environment
variable.
The structure of the YAML file
follows rules similar to the \code{config} R package, e.g.
projects can inherit settings from one another using the \code{inherits} field.
Exceptions include:
\enumerate{
\item There is no requirement to have a configuration named \code{"default"}.
\item Other projects do not inherit from the default project` automatically.
\item Not all fields need values because \code{targets} already has defaults.
}

\code{targets} does not actually invoke
the \code{config} package. The implementation in \code{targets}
was written from scratch without viewing or copying any
part of the source code of \code{config}.
}

\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(list(tar_target(x, 1 + 1)))
tar_config_get("store") # "_targets"
store_path <- tempfile()
tar_config_set(store = store_path)
tar_config_get("store") # Shows a temp file.
tar_make() # Writes to the custom data store identified in _targets.yaml.
tar_read(x) # tar_read() knows about _targets.yaml too.
file.exists("_targets") # FALSE
file.exists(store_path) # TRUE
})
}
}
\seealso{
Other configuration: 
\code{\link{tar_config_set}()},
\code{\link{tar_config_unset}()},
\code{\link{tar_envvars}()},
\code{\link{tar_option_get}()},
\code{\link{tar_option_reset}()},
\code{\link{tar_option_set}()}
}
\concept{configuration}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources_future.R
\name{tar_resources_future}
\alias{tar_resources_future}
\title{Target resources: \code{future} high-performance computing}
\usage{
tar_resources_future(plan = NULL, resources = list())
}
\arguments{
\item{plan}{A \code{future::plan()} object or \code{NULL},
a \code{target}-specific \code{future} plan.}

\item{resources}{Named list, \code{resources} argument to
\code{future::future()}. This argument is not supported in
some versions of \code{future}. For versions of \code{future}
where \code{resources} is not supported, instead supply \code{resources}
to \code{future::tweak()} and assign the returned plan to the \code{plan} argument
of \code{tar_resources_future()}.}
}
\value{
Object of class \code{"tar_resources_future"}, to be supplied
to the \code{future} argument of \code{tar_resources()}.
}
\description{
Create the \code{future} argument of \code{tar_resources()}
to specify optional high-performance computing settings
for \code{tar_make_future()}.
This is how to supply the \code{resources}
argument of \code{future::future()} for \code{targets}.
Resources supplied through
\code{future::plan()} and \code{future::tweak()} are completely ignored.
For details, see the documentation of the \code{future} R package
and the corresponding argument names in this help file.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
# Somewhere in you target script file (usually _targets.R):
tar_target(
  name,
  command(),
  resources = tar_resources(
    future = tar_resources_future(resources = list(n_cores = 2))
  )
)
}
\seealso{
Other resources: 
\code{\link{tar_resources_aws}()},
\code{\link{tar_resources_clustermq}()},
\code{\link{tar_resources_feather}()},
\code{\link{tar_resources_fst}()},
\code{\link{tar_resources_gcp}()},
\code{\link{tar_resources_parquet}()},
\code{\link{tar_resources_qs}()},
\code{\link{tar_resources_url}()},
\code{\link{tar_resources}()}
}
\concept{resources}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_prune.R
\name{tar_prune}
\alias{tar_prune}
\title{Remove targets that are no longer part of the pipeline.}
\usage{
tar_prune(
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
\code{NULL} except if \code{callr_function = callr::r_bg()}, in which case
a handle to the \code{callr} background process is returned. Either way,
the value is invisibly returned.
}
\description{
Remove target values from \verb{_targets/objects/} and
target metadata from \verb{_targets/meta/meta} for targets that are no longer
part of the pipeline.
}
\details{
This is useful if you recently worked through
multiple changes to your project and are now trying to
discard irrelevant data while keeping the results that still matter.
Global objects and dynamic files outside the
data store are unaffected. Also removes \verb{_targets/scratch/},
which is only needed while \code{\link[=tar_make]{tar_make()}}, \code{\link[=tar_make_clustermq]{tar_make_clustermq()}},
or \code{\link[=tar_make_future]{tar_make_future()}} is running.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
}, ask = FALSE)
tar_make()
# Remove some targets from the pipeline.
tar_script(list(tar_target(y1, 1 + 1)), ask = FALSE)
# Keep only the remaining targets in the data store.
tar_prune()
})
}
}
\seealso{
Other clean: 
\code{\link{tar_delete}()},
\code{\link{tar_destroy}()},
\code{\link{tar_invalidate}()}
}
\concept{clean}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_workspaces.R
\name{tar_workspaces}
\alias{tar_workspaces}
\title{List saved target workspaces.}
\usage{
tar_workspaces(names = NULL, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Optional \code{tidyselect} selector to return
a tactical subset of workspace names.
If \code{NULL}, all names are selected.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Character vector of available workspaces to load with
\code{\link[=tar_workspace]{tar_workspace()}}.
}
\description{
List target workspaces currently saved to
\verb{_targets/workspaces/}. See \code{\link[=tar_workspace]{tar_workspace()}} for more information.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set(workspace_on_error = TRUE)
  list(
    tar_target(x, "value"),
    tar_target(y, x)
  )
}, ask = FALSE)
tar_make()
tar_workspaces()
tar_workspaces(contains("x"))
})
}
}
\seealso{
Other debug: 
\code{\link{tar_load_globals}()},
\code{\link{tar_traceback}()},
\code{\link{tar_workspace}()}
}
\concept{debug}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources_qs.R
\name{tar_resources_qs}
\alias{tar_resources_qs}
\title{Target resources: qs storage formats}
\usage{
tar_resources_qs(preset = "high")
}
\arguments{
\item{preset}{Character of length 1, \code{preset}
argument of \code{qs::qsave()}.}
}
\value{
Object of class \code{"tar_resources_qs"}, to be supplied
to the qs argument of \code{tar_resources()}.
}
\description{
Create the \code{qs} argument of \code{tar_resources()}
to specify optional settings for big data storage formats
powered by the \code{qs} R package.
See the \code{format} argument of \code{\link[=tar_target]{tar_target()}} for details.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
# Somewhere in you target script file (usually _targets.R):
tar_target(
  name,
  command(),
  format = "qs",
  resources = tar_resources(
    qs = tar_resources_qs(preset = "fast")
  )
)
}
\seealso{
Other resources: 
\code{\link{tar_resources_aws}()},
\code{\link{tar_resources_clustermq}()},
\code{\link{tar_resources_feather}()},
\code{\link{tar_resources_fst}()},
\code{\link{tar_resources_future}()},
\code{\link{tar_resources_gcp}()},
\code{\link{tar_resources_parquet}()},
\code{\link{tar_resources_url}()},
\code{\link{tar_resources}()}
}
\concept{resources}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources_aws.R
\name{tar_resources_aws}
\alias{tar_resources_aws}
\title{Target resources: Amazon Web Services (AWS) storage formats}
\usage{
tar_resources_aws(
  bucket,
  prefix = targets::path_objects_dir_cloud(),
  region = NULL,
  part_size = 5 * (2^20)
)
}
\arguments{
\item{bucket}{Character of length 1, name of an existing
AWS S3 bucket to upload and download the return values
of the affected targets during the pipeline.}

\item{prefix}{Character of length 1, "directory path"
in the S3 bucket where the target return values are stored.}

\item{region}{Character of length 1, AWS region containing the S3 bucket.
Set to \code{NULL} to use the default region.}

\item{part_size}{Positive numeric of length 1, number of bytes
for each part of a multipart upload. (Except the last part,
which is the remainder.) In a multipart upload, each part
must be at least 5 MB.}
}
\value{
Object of class \code{"tar_resources_aws"}, to be supplied
to the \code{aws} argument of \code{tar_resources()}.
}
\description{
Create the \code{aws} argument of \code{tar_resources()}
to specify optional settings to AWS storage formats.
See the \code{format} argument of \code{\link[=tar_target]{tar_target()}} for details.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
# Somewhere in you target script file (usually _targets.R):
tar_target(
  name,
  command(),
  format = "aws_qs",
  resources = tar_resources(
    aws = tar_resources_aws(bucket = "yourbucketname"),
    qs = tar_resources_qs(preset = "fast")
  )
)
}
\seealso{
Other resources: 
\code{\link{tar_resources_clustermq}()},
\code{\link{tar_resources_feather}()},
\code{\link{tar_resources_fst}()},
\code{\link{tar_resources_future}()},
\code{\link{tar_resources_gcp}()},
\code{\link{tar_resources_parquet}()},
\code{\link{tar_resources_qs}()},
\code{\link{tar_resources_url}()},
\code{\link{tar_resources}()}
}
\concept{resources}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_canceled.R
\name{tar_canceled}
\alias{tar_canceled}
\title{List canceled targets.}
\usage{
tar_canceled(names = NULL, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Optional, names of the targets. If supplied, the
function restricts its output to these targets.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A character vector of canceled targets.
}
\description{
List targets whose progress is \code{"canceled"}.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
tar_canceled()
tar_canceled(starts_with("y_")) # see also all_of()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_deps_raw.R
\name{tar_deps_raw}
\alias{tar_deps_raw}
\title{Code dependencies (raw version)}
\usage{
tar_deps_raw(expr)
}
\arguments{
\item{expr}{An R expression object or function.}
}
\value{
Character vector of the dependencies of a function or expression.
}
\description{
Same as \code{\link[=tar_deps]{tar_deps()}} except \code{expr} must already be an
unquoted function or expression object.
}
\examples{
tar_deps_raw(quote(x <- y + z))
tar_deps_raw(
  quote({
    x <- 1
    x + a
  })
)
tar_deps_raw(function(a = b) map_dfr(data, ~do_row(.x)))
}
\seealso{
Other inspect: 
\code{\link{tar_deps}()},
\code{\link{tar_glimpse}()},
\code{\link{tar_manifest}()},
\code{\link{tar_network}()},
\code{\link{tar_outdated}()},
\code{\link{tar_sitrep}()},
\code{\link{tar_validate}()},
\code{\link{tar_visnetwork}()}
}
\concept{inspect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_helper_raw.R
\name{tar_helper_raw}
\alias{tar_helper_raw}
\title{Write a helper R script (raw version).}
\usage{
tar_helper_raw(path = NULL, code = NULL)
}
\arguments{
\item{path}{Character of length 1, path to write (or overwrite) \code{code}.
If the parent directory does not exist, \code{tar_helper_raw()} creates it.}

\item{code}{Expression object. \code{tar_helper_raw()} deparses and writes
this code to a file at \code{path}, overwriting it if the file already exists.}
}
\value{
\code{NULL} (invisibly)
}
\description{
Write a helper R script for a \code{targets} pipeline.
Could be supporting functions or the target script file
(default: \verb{_targets.R}) itself.
}
\details{
\code{tar_helper_raw()} is a specialized version of \code{\link[=tar_script]{tar_script()}}
with flexible paths and tidy evaluation. It is like \code{\link[=tar_helper]{tar_helper()}}
except that \code{code} is an "evaluated" argument rather than a quoted one.
}
\examples{
path <- tempfile()
tar_helper_raw(path, quote(x <- 1))
writeLines(readLines(path))
}
\seealso{
Other scripts: 
\code{\link{tar_edit}()},
\code{\link{tar_github_actions}()},
\code{\link{tar_helper}()},
\code{\link{tar_renv}()},
\code{\link{tar_script}()}
}
\concept{scripts}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_make_interactive.R
\name{tar_make_interactive}
\alias{tar_make_interactive}
\title{Interactive mode pipeline}
\usage{
tar_make_interactive(code)
}
\arguments{
\item{code}{Character vector of lines of a \verb{_targets.R} file
to define a pipeline.}
}
\value{
\code{NULL} (invisibly).
}
\description{
Not a user-side function. Do not invoke directly.
Only exported to on a technicality.
}
\examples{
if (identical(Sys.getenv("TAR_INTERACTIVE_EXAMPLES"), "true")) {
tar_make_interactive("library(targets); tar_target(x, 123)")
message(x)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources.R
\name{tar_resources}
\alias{tar_resources}
\title{Target resources}
\usage{
tar_resources(
  aws = NULL,
  clustermq = NULL,
  feather = NULL,
  fst = NULL,
  future = NULL,
  gcp = NULL,
  parquet = NULL,
  qs = NULL,
  url = NULL
)
}
\arguments{
\item{aws}{Output of function \code{tar_resources_aws()}.
AWS S3 storage settings for AWS backed storage formats
such as \code{"aws_qs"} and \verb{"aws_parquet}. Applies to all formats
beginning with the \code{"aws_"} prefix. For details on formats,
see the \code{format} argument of \code{\link[=tar_target]{tar_target()}}.}

\item{clustermq}{Output of function \code{tar_resources_clustermq()}.
Optional \code{clustermq} settings for \code{tar_make_clustermq()},
including the \code{log_worker} and \code{template} arguments of
\code{clustermq::workers()}.}

\item{feather}{Output of function \code{tar_resources_feather()}.
Non-default arguments to \code{arrow::read_feather()} and
\code{arrow::write_feather()} for \code{arrow}/feather-based storage formats.
Applies to all formats ending with the \code{"_feather"} suffix.
For details on formats, see the \code{format} argument of \code{\link[=tar_target]{tar_target()}}.}

\item{fst}{Output of function \code{tar_resources_fst()}.
Non-default arguments to \code{fst::read_fst()} and
\code{fst::write_fst()} for \code{fst}-based storage formats.
Applies to all formats ending with \code{"fst"} in the name.
For details on formats, see the \code{format} argument of \code{\link[=tar_target]{tar_target()}}.}

\item{future}{Output of function \code{tar_resources_future()}.
Optional \code{future} settings for \code{tar_make_future()},
including the \code{resources} argument of
\code{future::future()}, which can include values to insert in
template placeholders in \code{future.batchtools} template files.
This is how to supply the \code{resources}
argument of \code{future::future()} for \code{targets}.
Resources supplied through
\code{future::plan()} and \code{future::tweak()} are completely ignored.}

\item{gcp}{Output of function \code{tar_resources_gcp()}.
Google Cloud Platform bucket storage settings for GCP backed storage
formats such as \code{"gcp_qs"} and \verb{"gcp_parquet}. Applies to all formats
beginning with the \code{"gcp_"} prefix. For details on formats,
see the \code{format} argument of \code{\link[=tar_target]{tar_target()}}.}

\item{parquet}{Output of function \code{tar_resources_parquet()}.
Non-default arguments to \code{arrow::read_parquet()} and
\code{arrow::write_parquet()} for \code{arrow}/parquet-based storage formats.
Applies to all formats ending with the \code{"_parquet"} suffix.
For details on formats, see the \code{format} argument of \code{\link[=tar_target]{tar_target()}}.}

\item{qs}{Output of function \code{tar_resources_qs()}.
Non-default arguments to \code{qs::qread()} and
\code{qs::qsave()} for \code{qs}-based storage formats.
Applies to all formats ending with the \code{"_qs"} suffix.
For details on formats, see the \code{format} argument of \code{\link[=tar_target]{tar_target()}}.}

\item{url}{Output of function \code{tar_resources_url()}.
Non-default settings for storage formats ending with the \code{"_url"} suffix.
These settings include the \code{curl} handle for extra control over HTTP
requests. For details on formats, see the \code{format} argument of
\code{\link[=tar_target]{tar_target()}}.}
}
\value{
A list of objects of class \code{"tar_resources"} with
non-default settings of various optional backends for data storage
and high-performance computing.
}
\description{
Create a \code{resources} argument for \code{\link[=tar_target]{tar_target()}}
or \code{\link[=tar_option_set]{tar_option_set()}}.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
# Somewhere in you target script file (usually _targets.R):
tar_target(
  name,
  command(),
  format = "qs",
  resources = tar_resources(
    qs = tar_resources_qs(preset = "fast"),
    future = tar_resources_future(resources = list(n_cores = 1))
  )
)
}
\seealso{
Other resources: 
\code{\link{tar_resources_aws}()},
\code{\link{tar_resources_clustermq}()},
\code{\link{tar_resources_feather}()},
\code{\link{tar_resources_fst}()},
\code{\link{tar_resources_future}()},
\code{\link{tar_resources_gcp}()},
\code{\link{tar_resources_parquet}()},
\code{\link{tar_resources_qs}()},
\code{\link{tar_resources_url}()}
}
\concept{resources}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_manifest.R
\name{tar_manifest}
\alias{tar_manifest}
\title{Produce a data frame of information about your targets.}
\usage{
tar_manifest(
  names = NULL,
  fields = c("name", "command", "pattern"),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function),
  envir = parent.frame(),
  script = targets::tar_config_get("script")
)
}
\arguments{
\item{names}{Names of the targets to show. Set to \code{NULL} to
show all the targets (default). Otherwise, you can supply
symbols, a character vector, or \code{tidyselect} helpers like
\code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{fields}{Names of the fields, or columns, to show. Set to \code{NULL} to
show all the fields (default). Otherwise, you can supply
symbols, a character vector, or \code{tidyselect} helpers like \code{\link[=starts_with]{starts_with()}}.
Set to \code{NULL} to print all the fields.
The name of the target is always included as the first column
regardless of the selection.
Possible fields are below. All of them can be set in \code{\link[=tar_target]{tar_target()}},
\code{\link[=tar_target_raw]{tar_target_raw()}}, or \code{\link[=tar_option_set]{tar_option_set()}}.
\itemize{
\item \code{name}: Name of the target.
\item \code{command}: the R command that runs when the target builds.
\item \code{pattern}: branching pattern of the target, if applicable.
\item \code{format}: Storage format.
\item \code{iteration}: Iteration mode for branching.
\item \code{error}: Error mode, what to do when the target fails.
\item \code{memory}: Memory mode, when to keep targets in memory.
\item \code{storage}: Storage mode for high-performance computing scenarios.
\item \code{retrieval}: Retrieval mode for high-performance computing scenarios.
\item \code{deployment}: Where/whether to deploy the target in high-performance
computing scenarios.
\item \code{priority}: Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[=tar_make_future]{tar_make_future()}}).
\item \code{resources}: A list of target-specific resource requirements for
\code{\link[=tar_make_future]{tar_make_future()}}.
\item \code{cue_mode}: Cue mode from \code{\link[=tar_cue]{tar_cue()}}.
\item \code{cue_depend}: Depend cue from \code{\link[=tar_cue]{tar_cue()}}.
\item \code{cue_expr}: Command cue from \code{\link[=tar_cue]{tar_cue()}}.
\item \code{cue_file}: File cue from \code{\link[=tar_cue]{tar_cue()}}.
\item \code{cue_format}: Format cue from \code{\link[=tar_cue]{tar_cue()}}.
\item \code{cue_iteration}: Iteration cue from \code{\link[=tar_cue]{tar_cue()}}.
\item \code{packages}: List columns of packages loaded before building the target.
\item \code{library}: List column of library paths to load the packages.
}}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}
}
\value{
A data frame of information about the targets in the pipeline.
Rows appear in topological order (the order they will run
without any influence from parallel computing or priorities).
}
\description{
Along with \code{\link[=tar_visnetwork]{tar_visnetwork()}} and \code{\link[=tar_glimpse]{tar_glimpse()}},
\code{tar_manifest()} helps check that you constructed your pipeline correctly.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set()
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2),
    tar_target(m, z, pattern = map(z)),
    tar_target(c, z, pattern = cross(z))
  )
}, ask = FALSE)
tar_manifest()
tar_manifest(fields = c("name", "command"))
tar_manifest(fields = "command")
tar_manifest(fields = starts_with("cue"))
})
}
}
\seealso{
Other inspect: 
\code{\link{tar_deps_raw}()},
\code{\link{tar_deps}()},
\code{\link{tar_glimpse}()},
\code{\link{tar_network}()},
\code{\link{tar_outdated}()},
\code{\link{tar_sitrep}()},
\code{\link{tar_validate}()},
\code{\link{tar_visnetwork}()}
}
\concept{inspect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_exist_script.R
\name{tar_exist_script}
\alias{tar_exist_script}
\title{Check if the target script file exists.}
\usage{
tar_exist_script(script = targets::tar_config_get("script"))
}
\arguments{
\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}
}
\value{
Logical of length 1, whether the current project's metadata exists.
}
\description{
Check if the target script file exists for the
current project. The target script is \verb{_targets.R} by default,
but the path can be configured for the current project
using \code{\link[=tar_config_set]{tar_config_set()}}.
}
\examples{
tar_exist_script()
}
\seealso{
Other existence: 
\code{\link{tar_exist_meta}()},
\code{\link{tar_exist_objects}()},
\code{\link{tar_exist_process}()},
\code{\link{tar_exist_progress}()}
}
\concept{existence}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_assert.R
\name{tar_assert}
\alias{tar_assert}
\alias{tar_assert_chr}
\alias{tar_assert_dbl}
\alias{tar_assert_df}
\alias{tar_assert_equal_lengths}
\alias{tar_assert_envir}
\alias{tar_assert_expr}
\alias{tar_assert_flag}
\alias{tar_assert_file}
\alias{tar_assert_finite}
\alias{tar_assert_function}
\alias{tar_assert_function_arguments}
\alias{tar_assert_ge}
\alias{tar_assert_identical}
\alias{tar_assert_in}
\alias{tar_assert_not_dirs}
\alias{tar_assert_not_dir}
\alias{tar_assert_not_in}
\alias{tar_assert_inherits}
\alias{tar_assert_int}
\alias{tar_assert_internet}
\alias{tar_assert_lang}
\alias{tar_assert_le}
\alias{tar_assert_list}
\alias{tar_assert_lgl}
\alias{tar_assert_name}
\alias{tar_assert_names}
\alias{tar_assert_nonempty}
\alias{tar_assert_not_expr}
\alias{tar_assert_nzchar}
\alias{tar_assert_package}
\alias{tar_assert_path}
\alias{tar_assert_match}
\alias{tar_assert_nonmissing}
\alias{tar_assert_positive}
\alias{tar_assert_scalar}
\alias{tar_assert_target}
\alias{tar_assert_target_list}
\alias{tar_assert_true}
\alias{tar_assert_unique}
\alias{tar_assert_unique_targets}
\title{Assertions}
\usage{
tar_assert_chr(x, msg = NULL)

tar_assert_dbl(x, msg = NULL)

tar_assert_df(x, msg = NULL)

tar_assert_equal_lengths(x, msg = NULL)

tar_assert_envir(x, msg = NULL)

tar_assert_expr(x, msg = NULL)

tar_assert_flag(x, choices, msg = NULL)

tar_assert_file(x)

tar_assert_finite(x, msg = NULL)

tar_assert_function(x, msg = NULL)

tar_assert_function_arguments(x, args, msg = NULL)

tar_assert_ge(x, threshold, msg = NULL)

tar_assert_identical(x, y, msg = NULL)

tar_assert_in(x, choices, msg = NULL)

tar_assert_not_dirs(x, msg = NULL)

tar_assert_not_dir(x, msg = NULL)

tar_assert_not_in(x, choices, msg = NULL)

tar_assert_inherits(x, class, msg = NULL)

tar_assert_int(x, msg = NULL)

tar_assert_internet(msg = NULL)

tar_assert_lang(x, msg = NULL)

tar_assert_le(x, threshold, msg = NULL)

tar_assert_list(x, msg = NULL)

tar_assert_lgl(x, msg = NULL)

tar_assert_name(x)

tar_assert_names(x, msg = NULL)

tar_assert_nonempty(x, msg = NULL)

tar_assert_not_expr(x, msg = NULL)

tar_assert_nzchar(x, msg = NULL)

tar_assert_package(package)

tar_assert_path(path, msg = NULL)

tar_assert_match(x, pattern, msg = NULL)

tar_assert_nonmissing(x, msg = NULL)

tar_assert_positive(x, msg = NULL)

tar_assert_scalar(x, msg = NULL)

tar_assert_target(x, msg = NULL)

tar_assert_target_list(x)

tar_assert_true(x, msg = NULL)

tar_assert_unique(x, msg = NULL)

tar_assert_unique_targets(x)
}
\arguments{
\item{x}{R object, input to be validated. The kind of object depends on the
specific assertion function called.}

\item{msg}{Character of length 1, a message to be printed to the console
if \code{x} is invalid.}

\item{choices}{Character vector of choices of \code{x} for certain assertions.}

\item{args}{Character vector of expected function argument names.
Order matters.}

\item{threshold}{Numeric of length 1, lower/upper bound for
assertions like \code{tar_assert_le()}/\code{tar_assert_ge()}.}

\item{y}{R object, value to compare against \code{x}.}

\item{class}{Character vector of expected class names.}

\item{package}{Character of length 1, name of an R package.}

\item{path}{Character, file path.}

\item{pattern}{Character of length 1, a \code{grep} pattern for certain
assertions.}
}
\description{
These functions assert the correctness of user inputs
and generate custom error conditions as needed. Useful
for writing packages built on top of \code{targets}.
}
\examples{
tar_assert_chr("123")
try(tar_assert_chr(123))
}
\seealso{
Other utilities to extend targets: 
\code{\link{tar_condition}},
\code{\link{tar_dir}()},
\code{\link{tar_language}},
\code{\link{tar_test}()}
}
\concept{utilities to extend targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_network.R
\name{tar_network}
\alias{tar_network}
\title{Return the vertices and edges of a pipeline dependency graph.}
\usage{
tar_network(
  targets_only = FALSE,
  names = NULL,
  shortcut = FALSE,
  allow = NULL,
  exclude = NULL,
  outdated = TRUE,
  reporter = targets::tar_config_get("reporter_outdated"),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function, reporter),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{targets_only}{Logical, whether to restrict the output to just targets
(\code{FALSE}) or to also include imported global functions and objects.}

\item{names}{Names of targets. The graph visualization will operate
only on these targets (and unless \code{shortcut} is \code{TRUE},
all the targets upstream as well). Selecting a small subgraph
using \code{names} could speed up the load time of the visualization.
Unlike \code{allow}, \code{names} is invoked before the graph
is generated.
Set to NULL to check/build all the targets (default).
Otherwise, you can supply symbols or tidyselect helpers
like starts_with().
Applies to ordinary targets (stem) and whole dynamic branching
targets (patterns) but not individual dynamic branches.}

\item{shortcut}{Logical of length 1, how to interpret the \code{names} argument.
If \code{shortcut} is \code{FALSE} (default) then the function checks
all targets upstream of \code{names} as far back as the dependency graph goes.
If \code{TRUE}, then the function only checks the targets in \code{names}
and uses stored metadata for information about upstream dependencies
as needed. \code{shortcut = TRUE} increases speed if there are a lot of
up-to-date targets, but it assumes all the dependencies
are up to date, so please use with caution.
Also, \code{shortcut = TRUE} only works if you set \code{names}.}

\item{allow}{Optional, define the set of allowable vertices in the graph.
Unlike \code{names}, \code{allow} is invoked only after the graph is mostly
resolved, so it will not speed up execution.
Set to \code{NULL} to allow all vertices in the pipeline and environment
(default). Otherwise, you can supply symbols or
\code{tidyselect} helpers like \code{\link[=starts_with]{starts_with()}}.}

\item{exclude}{Optional, define the set of exclude vertices from the graph.
Unlike \code{names}, \code{exclude} is invoked only after the graph is mostly
resolved, so it will not speed up execution.
Set to \code{NULL} to exclude no vertices.
Otherwise, you can supply symbols or \code{tidyselect}
helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{outdated}{Logical, whether to show colors to distinguish outdated
targets from up-to-date targets. (Global functions and objects
still show these colors.) Looking for outdated targets
takes a lot of time for large pipelines with lots of branches,
and setting \code{outdated} to \code{FALSE} is a nice way to speed up the graph
if you only want to see dependency relationships and build progress.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets are checked. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"forecast"}: print running totals of the checked and outdated
targets found so far.
}}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A list with two data frames: \code{vertices} and \code{edges}. The
vertices data frame has one row per target with fields to denote
the type of the target or object (stem, branch, map, cross, function,
or object) and the target's status
(up to date, outdated, started, canceled, or errored).
The edges data frame has one row for every edge and columns \code{to} and
\code{from} to mark the starting and terminating vertices.
}
\description{
Analyze the pipeline defined in the target script file
(default: \verb{_targets.R})
and return the vertices and edges of the directed acyclic graph
of dependency relationships.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set()
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
}, ask = FALSE)
tar_network(targets_only = TRUE)
})
}
}
\seealso{
Other inspect: 
\code{\link{tar_deps_raw}()},
\code{\link{tar_deps}()},
\code{\link{tar_glimpse}()},
\code{\link{tar_manifest}()},
\code{\link{tar_outdated}()},
\code{\link{tar_sitrep}()},
\code{\link{tar_validate}()},
\code{\link{tar_visnetwork}()}
}
\concept{inspect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_delete.R
\name{tar_delete}
\alias{tar_delete}
\title{Delete locally stored target return values.}
\usage{
tar_delete(names, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Names of the targets to remove from \verb{_targets/objects/}.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\description{
Delete the return values of targets in \verb{_targets/objects/}.
but keep the records in \verb{_targets/meta/meta}.
}
\details{
If you have a small number of data-heavy targets you
need to discard to conserve storage, this function can help.
Dynamic files and cloud data (e.g. \code{format = "file"}
and \code{format = "aws_parquet"}) are not deleted.
For patterns recorded in the metadata, all the branches
will be deleted. For patterns no longer in the metadata,
branches are left alone.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
}, ask = FALSE)
tar_make()
tar_delete(starts_with("y")) # Only deletes y1 and y2.
tar_make() # y1 and y2 rebuild but return same values, so z is up to date.
})
}
}
\seealso{
Other clean: 
\code{\link{tar_destroy}()},
\code{\link{tar_invalidate}()},
\code{\link{tar_prune}()}
}
\concept{clean}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_destroy.R
\name{tar_destroy}
\alias{tar_destroy}
\title{Destroy all or part of the data store.}
\usage{
tar_destroy(
  destroy = c("all", "meta", "process", "progress", "objects", "scratch", "workspaces"),
  ask = NULL,
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{destroy}{Character of length 1, what to destroy. Choices:
\itemize{
\item \code{"all"}: destroy the entire data store (default: \verb{_targets/})
\item \code{"meta"}: just delete the metadata file at \code{meta/meta} in the
data store, which invalidates all the targets but keeps the data.
\item \code{"process"}: just delete the progress data file at
\code{meta/process} in the data store, which resets the metadata
of the main process.
\item \code{"progress"}: just delete the progress data file at
\code{meta/progress} in the data store,
which resets the progress tracking info.
\item \code{"objects"}: delete all the target
return values in \verb{objects/} in the data
store but keep progress and metadata.
Dynamic files are not deleted this way.
\item \code{"scratch"}: temporary files saved during \code{\link[=tar_make]{tar_make()}} that should
automatically get deleted except if R crashed.
\item \code{"workspaces"}: compressed files in \verb{workspaces/} in the data store with
the saved workspaces of targets. See \code{\link[=tar_workspace]{tar_workspace()}} for details.
}}

\item{ask}{Logical of length 1, whether to pause with a menu prompt
before deleting files. To disable this menu, set the \code{TAR_ASK}
environment variable to \code{"false"}. \code{usethis::edit_r_environ()}
can help set environment variables.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Nothing.
}
\description{
Destroy all or part of the data store written
by \code{\link[=tar_make]{tar_make()}} and similar functions.
}
\details{
\code{tar_destroy()} is a hard reset. Use it if you
intend to start the pipeline from scratch without
any trace of a previous run in \verb{_targets/}.
Global objects and dynamic files outside the
data store are unaffected.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(list(tar_target(x, 1 + 1)), ask = FALSE)
tar_make() # Creates the _targets/ data store.
tar_destroy()
print(file.exists("_targets")) # Should be FALSE.
})
}
}
\seealso{
Other clean: 
\code{\link{tar_delete}()},
\code{\link{tar_invalidate}()},
\code{\link{tar_prune}()}
}
\concept{clean}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_reprex.R
\name{tar_reprex}
\alias{tar_reprex}
\title{Reproducible example of \code{targets} with \code{reprex}}
\usage{
tar_reprex(pipeline = tar_target(example_target, 1), run = tar_make(), ...)
}
\arguments{
\item{pipeline}{R code for the target script file \verb{_targets.R}.
\code{library(targets)} is automatically written at the top.}

\item{run}{R code to inspect and run the pipeline.}

\item{...}{Named arguments passed to \code{reprex::reprex()}.}
}
\value{
A character vector of rendered the reprex, invisibly.
}
\description{
Create a reproducible example of a \code{targets}
pipeline with the \code{reprex} package.
}
\details{
The best way to get help with an issue is to
create a reproducible example of the problem
and post it to \url{https://github.com/ropensci/targets/discussions}
\code{tar_reprex()} facilitates this process. It is like
\code{reprex::reprex({targets::tar_script(...); tar_make()})},
but more convenient.
}
\examples{
if (identical(Sys.getenv("TAR_INTERACTIVE_EXAMPLES"), "true")) {
tar_reprex(
  pipeline = {
    list(
      tar_target(data, data.frame(x = sample.int(1e3))),
      tar_target(summary, mean(data$x, na.rm = TRUE))
    )
  },
  run = {
    tar_visnetwork()
    tar_make()
  }
)
}
}
\seealso{
Other help: 
\code{\link{targets-package}},
\code{\link{use_targets}()}
}
\concept{help}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_engine_knitr.R
\name{tar_engine_knitr}
\alias{tar_engine_knitr}
\title{Target Markdown \code{knitr} engine}
\usage{
tar_engine_knitr(options)
}
\arguments{
\item{options}{A named list of \code{knitr} chunk options.}
}
\value{
Character, output generated from \code{knitr::engine_output()}.
}
\description{
\code{knitr} language engine that runs \code{{targets}}
code chunks in Target Markdown.
}
\section{Target Markdown interactive mode}{

Target Markdown has two modes:
\enumerate{
\item Non-interactive mode. This is the default when you
run \code{knitr::knit()} or \code{rmarkdown::render()}.
Here, the code in \code{{targets}} code chunks gets written
to special script files in order to set up a \code{targets}
pipeline to run later.
\item Interactive mode: here, no scripts are written to set up
a pipeline. Rather, the globals or targets in question
are run in the current environment and the values
are assigned to that environment.
}

The mode is interactive if \code{!isTRUE(getOption("knitr.in.progress"))},
is \code{TRUE}. The \code{knitr.in.progress} option is \code{TRUE}
when you run \code{knitr::knit()} or \code{rmarkdown::render()}
and \code{NULL} if you are running one chunk at a time interactively
in an integrated development environment, e.g. the
notebook interface in RStudio:
\url{https://bookdown.org/yihui/rmarkdown/notebook.html}.
You can choose the mode with the \code{tar_interactive}
chunk option.
(In \code{targets} 0.6.0, \code{tar_interactive} defaults to \code{interactive()}
instead of \code{!isTRUE(getOption("knitr.in.progress"))}.)
}

\section{Target Markdown chunk options}{

Target Markdown introduces the following \code{knitr} code chunk options.
Most other standard \code{knitr} code chunk options should just work
in non-interactive mode. In interactive mode, not all
\itemize{
\item \code{tar_globals}: Logical of length 1,
whether to define globals or targets.
If \code{TRUE}, the chunk code defines functions, objects, and options
common to all the targets. If \code{FALSE} or \code{NULL} (default),
then the chunk returns formal targets for the pipeline.
\item \code{tar_interactive}: Logical of length 1, whether to run in
interactive mode or non-interactive mode.
See the "Target Markdown interactive mode" section of this
help file for details.
\item \code{tar_name}: name to use for writing helper script files
(e.g. \verb{_targets_r/targets/target_script.R})
and specifying target names if the \code{tar_simple} chunk option
is \code{TRUE}. All helper scripts and target names must have
unique names, so please do not set this option globally
with \code{knitr::opts_chunk$set()}.
\item \code{tar_script}: Character of length 1, where to write the
target script file in non-interactive mode. Most users can
skip this option and stick with the default \verb{_targets.R} script path.
Helper script files are always written next to the target script in
a folder with an \code{"_r"} suffix. The \code{tar_script} path must either be
absolute or be relative to the project root
(where you call \code{tar_make()} or similar).
If not specified, the target script path defaults to
\code{tar_config_get("script")} (default: \verb{_targets.R};
helpers default: \verb{_targets_r/}). When you run \code{tar_make()} etc.
with a non-default target script, you must select the correct target
script file either with the \code{script} argument or with
\code{tar_config_set(script = ...)}. The function will \code{source()}
the script file from the current working directory
(i.e. with \code{chdir = FALSE} in \code{source()}).
\item \code{tar_simple}: Logical of length 1.
Set to \code{TRUE} to define a single target with a simplified interface.
In code chunks with \code{tar_simple} equal to \code{TRUE}, the chunk label
(or the \code{tar_name} chunk option if you set it)
becomes the name, and the chunk code becomes the command.
In other words, a code chunk with label \code{targetname} and
command \code{mycommand()} automatically gets converted to
\code{tar_target(name = targetname, command = mycommand())}.
All other arguments of \code{tar_target()} remain at their default
values (configurable with \code{tar_option_set()} in a
\code{tar_globals = TRUE} chunk).
}
}

\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
# Register the engine.
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::knit_engines$set(targets = targets::tar_engine_knitr)
}
# Then, {targets} code chunks in a knitr report will run
# as described at https://books.ropensci.org/targets/markdown.html.
}
}
\seealso{
\url{https://books.ropensci.org/targets/markdown.html}

Other Target Markdown: 
\code{\link{tar_interactive}()},
\code{\link{tar_noninteractive}()},
\code{\link{tar_toggle}()}
}
\concept{Target Markdown}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_progress_branches.R
\name{tar_progress_branches}
\alias{tar_progress_branches}
\title{Tabulate the progress of dynamic branches.}
\usage{
tar_progress_branches(
  names = NULL,
  fields = NULL,
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Optional, names of the targets. If supplied, \code{tar_progress()}
only returns progress information on these targets.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=starts_with]{starts_with()}}.}

\item{fields}{Optional, names of progress data columns to read.
Set to \code{NULL} to read all fields.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A data frame with one row per target per progress status
and the following columns.
\itemize{
\item \code{name}: name of the pattern.
\item \code{progress}: progress status: \code{"started"}, \code{"built"}, \code{"cancelled"},
or \code{"errored"}.
\item \code{branches}: number of branches in the progress category.
\item \code{total}: total number of branches planned for the whole pattern.
Values within the same pattern should all be equal.
}
}
\description{
Read a project's target progress data for the most recent
run of the pipeline and display the tabulated status
of dynamic branches. Only the most recent record is shown.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, x, pattern = map(x)),
    tar_target(z, stopifnot(y < 1.5), pattern = map(y))
  )
}, ask = FALSE)
try(tar_make())
tar_progress_branches()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_definition.R
\name{tar_definition}
\alias{tar_definition}
\title{For developers only: get the definition of the current target.}
\usage{
tar_definition(
  default = targets::tar_target_raw("target_name", quote(identity()))
)
}
\arguments{
\item{default}{Environment, value to return if \code{tar_definition()}
is called on its own outside a \code{targets} pipeline.
Having a default lets users run things without \code{\link[=tar_make]{tar_make()}},
which helps peel back layers of code and troubleshoot bugs.}
}
\value{
If called from a running target, \code{tar_definition()} returns
the target object of the currently running target.
See the "Target objects" section for details.
}
\description{
For developers only: get the full definition of the
target currently running. This target definition is the same kind
of object produced by \code{\link[=tar_target]{tar_target()}}.
}
\details{
Most users should not use \code{tar_definition()}  because accidental
modifications could break the pipeline.
\code{tar_definition()} only exists in order to support third-party interface
packages, and even then the returned target definition is not modified..
}
\section{Target objects}{

Functions like \code{tar_target()} produce target objects,
special objects with specialized sets of S3 classes.
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
class(tar_definition())
tar_definition()$settings$name
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(
  tar_target(x, tar_definition()$settings$memory, memory = "transient")
)
tar_make(x)
tar_read(x)
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_active}()},
\code{\link{tar_call}()},
\code{\link{tar_cancel}()},
\code{\link{tar_envir}()},
\code{\link{tar_group}()},
\code{\link{tar_name}()},
\code{\link{tar_path}()},
\code{\link{tar_seed}()},
\code{\link{tar_store}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_store.R
\name{tar_store}
\alias{tar_store}
\title{Current data store path}
\usage{
tar_store()
}
\value{
Character, file path to the data store
of the pipeline currently running.
If called outside of the pipeline currently running,
\code{tar_store()} returns \code{tar_config_get("store")}.
}
\description{
Identify the file path to the data store
of the pipeline currently running.
}
\examples{
tar_store()
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(tar_target(x, tar_store()), ask = FALSE)
store <- tempfile()
tar_make(store = store)
tar_read(x, store = store)
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_active}()},
\code{\link{tar_call}()},
\code{\link{tar_cancel}()},
\code{\link{tar_definition}()},
\code{\link{tar_envir}()},
\code{\link{tar_group}()},
\code{\link{tar_name}()},
\code{\link{tar_path}()},
\code{\link{tar_seed}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_config_set.R
\name{tar_config_set}
\alias{tar_config_set}
\title{Set configuration settings.}
\usage{
tar_config_set(
  inherits = NULL,
  reporter_make = NULL,
  reporter_outdated = NULL,
  store = NULL,
  shortcut = NULL,
  script = NULL,
  workers = NULL,
  config = Sys.getenv("TAR_CONFIG", "_targets.yaml"),
  project = Sys.getenv("TAR_PROJECT", "main")
)
}
\arguments{
\item{inherits}{Character of length 1, name of the project from which
the current project should inherit configuration settings.
The current project is the \code{project} argument, which
defaults to \code{Sys.getenv("TAR_PROJECT", "main")}.
If the \code{inherits} argument \code{NULL}, the \code{inherits} setting is not modified.
Use \code{\link[=tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{reporter_make}{Character of length 1, \code{reporter} argument to
\code{\link[=tar_make]{tar_make()}} and related functions that run the pipeline.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[=tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{reporter_outdated}{Character of length 1, \code{reporter} argument to
\code{\link[=tar_outdated]{tar_outdated()}} and related functions that do not run the pipeline.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[=tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{store}{Character of length 1, path to the data store of the pipeline.
If \code{NULL}, the \code{store} setting is left unchanged in the
YAML configuration file (default: \verb{_targets.yaml}).
Usually, the data store lives at \verb{_targets}.
Set \code{store} to a custom directory
to specify a path other than \verb{_targets/}. The path need not exist
before the pipeline begins, and it need not end with "_targets",
but it must be writeable.
For optimal performance, choose a storage location
with fast read/write access.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[=tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{shortcut}{logical of length 1, default \code{shortcut} argument
to \code{\link[=tar_make]{tar_make()}} and related functions.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[=tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{script}{Character of length 1, path to the target script file
that defines the pipeline (\verb{_targets.R} by default).
This path should be either
an absolute path or a path relative to the project root where you will
call \code{\link[=tar_make]{tar_make()}} and other functions. When \code{\link[=tar_make]{tar_make()}} and friends
run the script from the current working directory.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[=tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{workers}{Positive numeric of length 1, \code{workers} argument of
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and related functions that run the pipeline
with parallel computing among targets.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[=tar_config_unset]{tar_config_unset()}} to delete a setting.}

\item{config}{Character of length 1, file path of the YAML
configuration file with \code{targets} project settings.
The \code{config} argument specifies which YAML configuration
file that \code{tar_config_get()} reads from or \code{tar_config_set()}
writes to in a single function call.
It does not globally change which configuration file is used
in subsequent function calls. The default file path of the YAML
file is always \verb{_targets.yaml} unless you set another
default path using the \code{TAR_CONFIG} environment variable,
e.g. \code{Sys.setenv(TAR_CONFIG = "custom.yaml")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}

\item{project}{Character of length 1, name of the current
\code{targets} project. Thanks to the \code{config} R package,
\code{targets} YAML configuration files can store multiple
sets of configuration settings, with each set corresponding
to its own project. The \code{project} argument allows you to
set or get a configuration setting for a specific project
for a given call to \code{tar_config_set()} or \code{tar_config_get()}.
The default project is always called \code{"main"}
unless you set another
default project using the \code{TAR_PROJECT} environment variable,
e.g. \code{Sys.setenv(tar_project = "custom")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}
}
\value{
\code{NULL} (invisibly)
}
\description{
\code{tar_config_set()} writes special custom settings
for the current project to an optional YAML configuration file.
}
\section{Configuration}{

For several key functions like \code{\link[=tar_make]{tar_make()}}, the
default values of arguments are controlled though
\code{tar_config_get()}. \code{tar_config_get()} retrieves data
from an optional YAML configuration file.
You can control the settings in the YAML
file programmatically with \code{tar_config_set()}.
The default file path of this YAML file is \verb{_targets.yaml}, and you can
set another path globally using the \code{TAR_CONFIG}
environment variable. The YAML file can store configuration
settings for multiple projects, and you can globally
set the default project with the \code{TAR_PROJECT} environment
variable.
The structure of the YAML file
follows rules similar to the \code{config} R package, e.g.
projects can inherit settings from one another using the \code{inherits} field.
Exceptions include:
\enumerate{
\item There is no requirement to have a configuration named \code{"default"}.
\item Other projects do not inherit from the default project` automatically.
\item Not all fields need values because \code{targets} already has defaults.
}

\code{targets} does not actually invoke
the \code{config} package. The implementation in \code{targets}
was written from scratch without viewing or copying any
part of the source code of \code{config}.
}

\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(list(tar_target(x, 1 + 1)))
tar_config_get("store") # NULL (data store defaults to "_targets/")
store_path <- tempfile()
tar_config_set(store = store_path)
tar_config_get("store") # Shows a temp file.
tar_make() # Writes to the custom data store identified in _targets.yaml.
tar_read(x) # tar_read() knows about _targets.yaml too.
file.exists("_targets") # FALSE
file.exists(store_path) # TRUE
})
}
}
\seealso{
Other configuration: 
\code{\link{tar_config_get}()},
\code{\link{tar_config_unset}()},
\code{\link{tar_envvars}()},
\code{\link{tar_option_get}()},
\code{\link{tar_option_reset}()},
\code{\link{tar_option_set}()}
}
\concept{configuration}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_watch_server.R
\name{tar_watch_server}
\alias{tar_watch_server}
\title{Shiny module server for tar_watch()}
\usage{
tar_watch_server(
  id,
  height = "650px",
  exclude = ".Random.seed",
  config = Sys.getenv("TAR_CONFIG", "_targets.yaml"),
  project = Sys.getenv("TAR_PROJECT", "main")
)
}
\arguments{
\item{id}{Character of length 1, ID corresponding to the UI function
of the module.}

\item{height}{Character of length 1,
height of the \code{visNetwork} widget and branches table.}

\item{exclude}{Character vector of nodes to omit from the graph.}

\item{config}{Character of length 1, file path of the YAML
configuration file with \code{targets} project settings.
The \code{config} argument specifies which YAML configuration
file that \code{tar_config_get()} reads from or \code{tar_config_set()}
writes to in a single function call.
It does not globally change which configuration file is used
in subsequent function calls. The default file path of the YAML
file is always \verb{_targets.yaml} unless you set another
default path using the \code{TAR_CONFIG} environment variable,
e.g. \code{Sys.setenv(TAR_CONFIG = "custom.yaml")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}

\item{project}{Character of length 1, name of the current
\code{targets} project. Thanks to the \code{config} R package,
\code{targets} YAML configuration files can store multiple
sets of configuration settings, with each set corresponding
to its own project. The \code{project} argument allows you to
set or get a configuration setting for a specific project
for a given call to \code{tar_config_set()} or \code{tar_config_get()}.
The default project is always called \code{"main"}
unless you set another
default project using the \code{TAR_PROJECT} environment variable,
e.g. \code{Sys.setenv(tar_project = "custom")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}
}
\value{
A Shiny module server.
}
\description{
Use \code{\link[=tar_watch_ui]{tar_watch_ui()}} and \code{tar_watch_server()}
to include \code{\link[=tar_watch]{tar_watch()}} as a Shiny module in an app.
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_script.R
\name{tar_script}
\alias{tar_script}
\title{Write a target script file.}
\usage{
tar_script(
  code = NULL,
  library_targets = TRUE,
  ask = NULL,
  script = targets::tar_config_get("script")
)
}
\arguments{
\item{code}{R code to write to the target script file.
If \code{NULL}, an example target script file is written instead.}

\item{library_targets}{logical, whether to write a \code{library(targets)}
line at the top of the target script file automatically (recommended).
If \code{TRUE}, you do not need to explicitly put \code{library(targets)}
in \code{code}.}

\item{ask}{Logical, whether to ask before writing if the
target script file
already exists. If \code{NULL}, defaults to \code{Sys.getenv("TAR_ASK")}.
(Set to \code{"true"} or \code{"false"} with \code{Sys.setenv()}).
If \code{ask} and the \code{TAR_ASK} environment variable are both
indeterminate, defaults to \code{interactive()}.}

\item{script}{Character of length 1, where to write
the target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}.}
}
\value{
\code{NULL} (invisibly).
}
\description{
The \code{tar_script()} function is a convenient
way to create the required target script file (default: \verb{_targets.R})
in the current working directory.
It always overwrites the existing target script,
and it requires you to be in the working directory
where you intend to write the file, so be careful.
See the "Target script" section for details.
}
\section{Target script file}{

Every \code{targets} project requires a target script file.
The target script file is usually a file called \verb{_targets.R}
Functions \code{\link[=tar_make]{tar_make()}} and friends look for the target script
and run it to set up the pipeline just prior to the main task.
Every target script file should run the following
steps in the order below:
1. Package: load the \code{targets} package. This step is automatically
inserted at the top of the target script file produced by
\code{tar_script()} if \code{library_targets} is \code{TRUE},
so you do not need to explicitly include it in \code{code}.
1. Globals: load custom functions and global objects into memory.
Usually, this section is a bunch of calls to \code{source()} that run
scripts defining user-defined functions. These functions support
the R commands of the targets.
2. Options: call \code{\link[=tar_option_set]{tar_option_set()}} to set defaults for targets-specific
settings such as the names of required packages. Even if you have no
specific options to set, it is still recommended to call
\code{\link[=tar_option_set]{tar_option_set()}} in order to register the proper environment.
3. Targets: define one or more target objects using \code{\link[=tar_target]{tar_target()}}.
4. Pipeline: call \code{\link[=list]{list()}} to bring the targets from (3)
together in a pipeline object. Every target script file must return
a pipeline object, which usually means ending with a call to
\code{\link[=list]{list()}}. In practice, (3) and (4) can be combined together
in the same function call.
}

\examples{
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script() # Writes an example target script file.
# Writes a user-defined target script:
tar_script({
  x <- tar_target(x, 1 + 1)
  tar_option_set()
  list(x)
}, ask = FALSE)
writeLines(readLines("_targets.R"))
})
}
\seealso{
Other scripts: 
\code{\link{tar_edit}()},
\code{\link{tar_github_actions}()},
\code{\link{tar_helper_raw}()},
\code{\link{tar_helper}()},
\code{\link{tar_renv}()}
}
\concept{scripts}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio_addin_tar_target.R
\name{rstudio_addin_tar_target}
\alias{rstudio_addin_tar_target}
\title{RStudio addin to insert \code{"tar_target()"} at the cursor.}
\usage{
rstudio_addin_tar_target(context = NULL)
}
\arguments{
\item{context}{RStudio API context from
\code{rstudioapi::getActiveDocumentContext()}.}
}
\description{
For internal use only. Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_pipeline.R
\name{as_pipeline}
\alias{as_pipeline}
\title{Convert to a pipeline object.}
\usage{
as_pipeline(x)
}
\arguments{
\item{x}{A list of target objects or a pipeline object.}
}
\value{
An object of class \code{"tar_pipeline"}.
}
\description{
Not a user-side function. Do not invoke directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_name.R
\name{tar_name}
\alias{tar_name}
\title{Get the name of the target currently running.}
\usage{
tar_name(default = "target")
}
\arguments{
\item{default}{Character, value to return if \code{tar_name()}
is called on its own outside a \code{targets} pipeline.
Having a default lets users run things without \code{\link[=tar_make]{tar_make()}},
which helps peel back layers of code and troubleshoot bugs.}
}
\value{
Character of length 1. If called inside a pipeline,
\code{tar_name()} returns name of the target currently running.
Otherwise, the return value is \code{default}.
}
\description{
Get the name of the target currently running.
}
\examples{
tar_name()
tar_name(default = "custom_target_name")
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(tar_target(x, tar_name()), ask = FALSE)
tar_make()
tar_read(x)
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_active}()},
\code{\link{tar_call}()},
\code{\link{tar_cancel}()},
\code{\link{tar_definition}()},
\code{\link{tar_envir}()},
\code{\link{tar_group}()},
\code{\link{tar_path}()},
\code{\link{tar_seed}()},
\code{\link{tar_store}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_interactive.R
\name{tar_interactive}
\alias{tar_interactive}
\title{Run if Target Markdown interactive mode is on.}
\usage{
tar_interactive(code)
}
\arguments{
\item{code}{R code to run if Target Markdown interactive mode
is turned on.}
}
\value{
If Target Markdown interactive mode is turned on,
the function returns the result of running the code.
Otherwise, the function invisibly returns \code{NULL}.
}
\description{
In Target Markdown, run the enclosed code
only if interactive mode is activated. Otherwise,
do not run the code.
}
\details{
Visit <books.ropensci.org/targets/markdown.html>
to learn about Target Markdown and interactive mode.
}
\examples{
tar_interactive(message("In interactive mode."))
}
\seealso{
Other Target Markdown: 
\code{\link{tar_engine_knitr}()},
\code{\link{tar_noninteractive}()},
\code{\link{tar_toggle}()}
}
\concept{Target Markdown}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_call.R
\name{tar_call}
\alias{tar_call}
\title{Identify the called \code{targets} function.}
\usage{
tar_call()
}
\value{
Character of length 1, name of the currently running \code{targets}
interface function. For example, suppose you have a call to
\code{tar_call()} inside a target or \verb{_targets.R}. Then if you run
\code{tar_make()}, \code{tar_call()} will return \code{"tar_make"}.
}
\description{
Get the name of the currently running \code{targets}
interface function. Returns \code{NULL} if not invoked inside
a target or \verb{_targets.R} (i.e. if not directly invoked
by \code{\link[=tar_make]{tar_make()}}, \code{\link[=tar_visnetwork]{tar_visnetwork()}}, etc.).
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_call() # NULL
tar_script({
  message("called function: ", tar_call())
  tar_target(x, tar_call())
})
tar_manifest() # prints "called function: tar_manifest"
tar_make() # prints "called function: tar_make"
tar_read(x) # "tar_make"
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_active}()},
\code{\link{tar_cancel}()},
\code{\link{tar_definition}()},
\code{\link{tar_envir}()},
\code{\link{tar_group}()},
\code{\link{tar_name}()},
\code{\link{tar_path}()},
\code{\link{tar_seed}()},
\code{\link{tar_store}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_built.R
\name{tar_built}
\alias{tar_built}
\title{List built targets.}
\usage{
tar_built(names = NULL, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Optional, names of the targets. If supplied, the
function restricts its output to these targets.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A character vector of built targets.
}
\description{
List targets whose progress is \code{"built"}.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
tar_built()
tar_built(starts_with("y_")) # see also all_of()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_read.R
\name{tar_read}
\alias{tar_read}
\title{Read a target's value from storage.}
\usage{
tar_read(
  name,
  branches = NULL,
  meta = tar_meta(store = store),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{name}{Symbol, name of the target to read.}

\item{branches}{Integer of indices of the branches to load
if the target is a pattern.}

\item{meta}{Data frame of metadata from \code{\link[=tar_meta]{tar_meta()}}.
\code{tar_read()} with the default arguments can be inefficient for large
pipelines because all the metadata is stored in a single file.
However, if you call \code{\link[=tar_meta]{tar_meta()}} beforehand and supply it to the \code{meta}
argument, then successive calls to \code{tar_read()} may run much faster.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
The target's return value from its file in
\verb{_targets/objects/}, or the paths to the custom files and directories
if \code{format = "file"} was set.
}
\description{
Read a target's return value from its file in
\verb{_targets/objects/}. For dynamic files (i.e. \code{format = "file"})
the paths are returned.
}
\section{Limited scope}{

\code{tar_read()} and \code{tar_load()}
are only for exploratory analysis and literate programming,
and \code{tar_read_raw()} and \code{tar_load_raw()} are only
for exploratory analysis. \code{targets} automatically
loads the correct dependencies into memory when the pipeline
is running, so invoking these functions
from inside a target is rarely advisable.
}

\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(list(tar_target(x, 1 + 1)), ask = FALSE)
tar_make()
tar_read(x)
})
}
}
\seealso{
Other data: 
\code{\link{tar_load_raw}()},
\code{\link{tar_load}()},
\code{\link{tar_meta}()},
\code{\link{tar_objects}()},
\code{\link{tar_pid}()},
\code{\link{tar_process}()},
\code{\link{tar_read_raw}()}
}
\concept{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_progress.R
\name{tar_progress}
\alias{tar_progress}
\title{Read progress.}
\usage{
tar_progress(
  names = NULL,
  fields = "progress",
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Optional, names of the targets. If supplied, \code{tar_progress()}
only returns progress information on these targets.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{fields}{Optional, names of progress data columns to read.
Set to \code{NULL} to read all fields.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A data frame with one row per target and the following columns:
\itemize{
\item \code{name}: name of the target.
\item \code{type}: type of target: \code{"stem"} for non-branching targets,
\code{"pattern"} for dynamically branching targets, and \code{"branch"}
for dynamic branches.
\item \code{parent}: name of the target's parent. For branches, this is the
name of the associated pattern. For other targets, the pattern
is just itself.
\item \code{branches}: number of dynamic branches of a pattern. 0 for non-patterns.
\item \code{progress}: the most recent progress update of that target.
Could be \code{"started"}, \code{"built"}, "\code{skipped}", \code{"canceled"},
or \code{"errored"}.
}
}
\description{
Read a project's target progress data for the most recent
run of \code{\link[=tar_make]{tar_make()}} or similar. Only the most recent record is shown.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
tar_progress()
tar_progress(starts_with("y_")) # see also all_of()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources_url.R
\name{tar_resources_url}
\alias{tar_resources_url}
\title{Target resources: URL storage formats}
\usage{
tar_resources_url(handle = NULL)
}
\arguments{
\item{handle}{Object returned by \code{curl::new_handle} or \code{NULL}.}
}
\value{
Object of class \code{"tar_resources_url"}, to be supplied
to the url argument of \code{tar_resources()}.
}
\description{
Create the \code{url} argument of \code{tar_resources()}
to specify optional settings for URL storage formats.
See the \code{format} argument of \code{\link[=tar_target]{tar_target()}} for details.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
# Somewhere in you target script file (usually _targets.R):
tar_target(
  name,
  command(),
  format = "url",
  resources = tar_resources(
    url = tar_resources_url(handle = curl::new_handle())
  )
)
}
}
\seealso{
Other resources: 
\code{\link{tar_resources_aws}()},
\code{\link{tar_resources_clustermq}()},
\code{\link{tar_resources_feather}()},
\code{\link{tar_resources_fst}()},
\code{\link{tar_resources_future}()},
\code{\link{tar_resources_gcp}()},
\code{\link{tar_resources_parquet}()},
\code{\link{tar_resources_qs}()},
\code{\link{tar_resources}()}
}
\concept{resources}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_older.R
\name{tar_older}
\alias{tar_older}
\title{List old targets}
\usage{
tar_older(
  time,
  names = NULL,
  inclusive = FALSE,
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{time}{A \code{POSIXct} object of length 1, time threshold.
Targets older than this time stamp are returned.
For example, if \code{time = Sys.time() - as.difftime(1, units = "weeks")}
then \code{tar_older()} returns targets older than one week ago.}

\item{names}{Names of eligible targets. Targets excluded from \code{names}
will not be returned even if they are old.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
If \code{NULL}, all names are eligible.}

\item{inclusive}{Logical of length 1, whether to include targets
built at exactly the \code{time} given.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A character vector of names of old targets with recorded
timestamp metadata.
}
\description{
List all the targets whose last successful run occurred
before a certain point in time. Combine with \code{\link[=tar_invalidate]{tar_invalidate()}},
you can use \code{tar_older()} to automatically rerun targets at
regular intervals. See the examples for a demonstration.
}
\details{
Only applies to targets with recorded time stamps:
just non-branching targets and individual dynamic branches.
As of \code{targets} version 0.6.0, these time
stamps are available for these targets regardless of
storage format. Earlier versions of \code{targets} do not record
time stamps for remote storage formats such as \code{"url"}
or any of the \code{"aws_*"} formats.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(tar_target(x, seq_len(2)))
}, ask = FALSE)
tar_make()
# targets older than 1 week ago
tar_older(Sys.time() - as.difftime(1, units = "weeks"))
# targets older than 1 week from now
tar_older(Sys.time() + as.difftime(1, units = "weeks"))
# Everything is still up to date.
tar_make()
# Invalidate all targets targets older than 1 week from now
# so they run on the next tar_make().
invalidate_these <- tar_older(Sys.time() + as.difftime(1, units = "weeks"))
tar_invalidate(all_of(invalidate_these))
tar_make()
})
}
}
\seealso{
Other time: 
\code{\link{tar_newer}()},
\code{\link{tar_timestamp_raw}()},
\code{\link{tar_timestamp}()}
}
\concept{time}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_traceback.R
\name{tar_traceback}
\alias{tar_traceback}
\title{Get a target's traceback}
\usage{
tar_traceback(
  name,
  envir = NULL,
  packages = NULL,
  source = NULL,
  characters = getOption("width"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{name}{Symbol, name of the target whose workspace to read.}

\item{envir}{Deprecated in \code{targets} > 0.3.1 (2021-03-28).}

\item{packages}{Logical, whether to load the required packages
of the target.}

\item{source}{Logical, whether to run the target script file
(default: \verb{_targets.R}) to load user-defined
global object dependencies into \code{envir}. If \code{TRUE}, then \code{envir}
should either be the global environment or inherit from the
global environment.}

\item{characters}{Positive integer. Each line of the traceback
is shortened to this number of characters.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Character vector, the traceback of a failed target
if it exists.
}
\description{
Return the saved traceback of a target.
Assumes the target errored out in a previous run of the pipeline
with workspaces enabled for that target.
See \code{\link[=tar_workspace]{tar_workspace()}} for details.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tmp <- sample(1)
tar_script({
  tar_option_set(workspace_on_error = TRUE)
  list(
    tar_target(x, "loaded"),
    tar_target(y, stop(x))
  )
}, ask = FALSE)
try(tar_make())
tar_traceback(y, characters = 60)
})
}
}
\seealso{
Other debug: 
\code{\link{tar_load_globals}()},
\code{\link{tar_workspaces}()},
\code{\link{tar_workspace}()}
}
\concept{debug}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources_gcp.R
\name{tar_resources_gcp}
\alias{tar_resources_gcp}
\title{Target resources: Google Cloud Platform storage formats}
\usage{
tar_resources_gcp(
  bucket,
  prefix = targets::path_objects_dir_cloud(),
  verbose = FALSE
)
}
\arguments{
\item{bucket}{Character of length 1, name of an existing
gcp bucket to upload and download the return values
of the affected targets during the pipeline.}

\item{prefix}{Character of length 1, "directory path"
in the S3 bucket where the target return values are stored.}

\item{verbose}{Whether to have feedback on the GCS upload process}
}
\value{
Object of class \code{"tar_resources_gcp"}, to be supplied
to the \code{gcp} argument of \code{tar_resources()}.
}
\description{
Create the \code{gcp} argument of \code{tar_resources()}
to specify optional settings to gcp storage formats.
See the \code{format} argument of \code{\link[=tar_target]{tar_target()}} for details.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
# Somewhere in you target script file (usually _targets.R):
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_target(
  name,
  command(),
  format = "gcp_qs",
  resources = tar_resources(
    gcp = tar_resources_gcp(bucket = "yourbucketname"),
    qs = tar_resources_qs(preset = "fast")
  )
)
}
}
\seealso{
Other resources: 
\code{\link{tar_resources_aws}()},
\code{\link{tar_resources_clustermq}()},
\code{\link{tar_resources_feather}()},
\code{\link{tar_resources_fst}()},
\code{\link{tar_resources_future}()},
\code{\link{tar_resources_parquet}()},
\code{\link{tar_resources_qs}()},
\code{\link{tar_resources_url}()},
\code{\link{tar_resources}()}
}
\concept{resources}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_glimpse.R
\name{tar_glimpse}
\alias{tar_glimpse}
\title{Visualize an abridged fast dependency graph.}
\usage{
tar_glimpse(
  targets_only = TRUE,
  names = NULL,
  shortcut = FALSE,
  allow = NULL,
  exclude = ".Random.seed",
  level_separation = NULL,
  degree_from = 1L,
  degree_to = 1L,
  zoom_speed = 1,
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{targets_only}{Logical, whether to restrict the output to just targets
(\code{FALSE}) or to also include global functions and objects.}

\item{names}{Names of targets. The graph visualization will operate
only on these targets (and unless \code{shortcut} is \code{TRUE},
all the targets upstream as well). Selecting a small subgraph
using \code{names} could speed up the load time of the visualization.
Unlike \code{allow}, \code{names} is invoked before the graph
is generated.
Set to NULL to check/build all the targets (default).
Otherwise, you can supply symbols or tidyselect helpers
like starts_with().
Applies to ordinary targets (stem) and whole dynamic branching
targets (patterns) but not individual dynamic branches.}

\item{shortcut}{Logical of length 1, how to interpret the \code{names} argument.
If \code{shortcut} is \code{FALSE} (default) then the function checks
all targets upstream of \code{names} as far back as the dependency graph goes.
If \code{TRUE}, then the function only checks the targets in \code{names}
and uses stored metadata for information about upstream dependencies
as needed. \code{shortcut = TRUE} increases speed if there are a lot of
up-to-date targets, but it assumes all the dependencies
are up to date, so please use with caution.
Also, \code{shortcut = TRUE} only works if you set \code{names}.}

\item{allow}{Optional, define the set of allowable vertices in the graph.
Unlike \code{names}, \code{allow} is invoked only after the graph is mostly
resolved, so it will not speed up execution.
Set to \code{NULL} to allow all vertices in the pipeline and environment
(default). Otherwise, you can supply symbols or
\code{tidyselect} helpers like \code{\link[=starts_with]{starts_with()}}.}

\item{exclude}{Optional, define the set of exclude vertices from the graph.
Unlike \code{names}, \code{exclude} is invoked only after the graph is mostly
resolved, so it will not speed up execution.
Set to \code{NULL} to exclude no vertices.
Otherwise, you can supply symbols or \code{tidyselect}
helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{level_separation}{Numeric of length 1,
\code{levelSeparation} argument of \code{visNetwork::visHierarchicalLayout()}.
Controls the distance between hierarchical levels.
Consider changing the value if the aspect ratio of the graph
is far from 1. If \code{level_separation} is \code{NULL},
the \code{levelSeparation} argument of \code{visHierarchicalLayout()}
defaults to \code{150}.}

\item{degree_from}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_from}
controls the number of edges the neighborhood extends upstream.}

\item{degree_to}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_to}
controls the number of edges the neighborhood extends downstream.}

\item{zoom_speed}{Positive numeric of length 1, scaling factor on the
zoom speed. Above 1 zooms faster than default, below 1 zooms
lower than default.}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A \code{visNetwork} HTML widget object.
}
\description{
Analyze the pipeline defined in the target script file
(default: \verb{_targets.R})
and visualize the directed acyclic graph of targets.
Unlike \code{\link[=tar_visnetwork]{tar_visnetwork()}}, \code{tar_glimpse()} does not account for
metadata or progress information, which means the graph
renders faster. Also, \code{tar_glimpse()} omits functions and other global
objects by default (but you can include them with \code{targets_only = FALSE}).
}
\examples{
if (identical(Sys.getenv("TAR_INTERACTIVE_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set()
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
}, ask = FALSE)
tar_glimpse()
tar_glimpse(allow = starts_with("y")) # see also all_of()
})
}
}
\seealso{
Other inspect: 
\code{\link{tar_deps_raw}()},
\code{\link{tar_deps}()},
\code{\link{tar_manifest}()},
\code{\link{tar_network}()},
\code{\link{tar_outdated}()},
\code{\link{tar_sitrep}()},
\code{\link{tar_validate}()},
\code{\link{tar_visnetwork}()}
}
\concept{inspect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_workspace.R
\name{tar_workspace}
\alias{tar_workspace}
\title{Load a saved workspace and seed for debugging.}
\usage{
tar_workspace(
  name,
  envir = parent.frame(),
  packages = TRUE,
  source = TRUE,
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{name}{Symbol, name of the target whose workspace to read.}

\item{envir}{Environment in which to put the objects.}

\item{packages}{Logical, whether to load the required packages
of the target.}

\item{source}{Logical, whether to run \verb{_targets.R} to load user-defined
global object dependencies into \code{envir}. If \code{TRUE}, then \code{envir}
should either be the global environment or inherit from the
global environment.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
This function returns \code{NULL}, but it does load
the target's required packages, as well as multiple objects
into the environment (\code{envir} argument) in order to replicate the
workspace where the error happened. These objects include
the global objects at the time \code{\link[=tar_make]{tar_make()}} was called and the
dependency targets. The random number generator seed for the
target is also assigned with \code{set.seed()}.
}
\description{
Load the packages, workspace, and random number generator seed
of target attempted with a workspace file.
}
\details{
If you activate workspaces through the \code{workspaces} argument
of \code{\link[=tar_option_set]{tar_option_set()}}, then under the circumstances you specify,
\code{targets} will save a special workspace file to a location in
in \verb{_targets/workspaces/}. The workspace file is a compact reference
that allows \code{tar_workspace()} to load the target's dependencies
and random number generator seed as long as the data objects
are still in the data store (usually files in \verb{_targets/objects/}).
When you are done debugging, you can remove the workspace files
using \code{tar_destroy(destroy = "workspaces")}.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tmp <- sample(1)
tar_script({
  tar_option_set(workspace_on_error = TRUE)
  list(
    tar_target(x, "loaded"),
    tar_target(y, stop(x))
  )
}, ask = FALSE)
# The following code throws an error for demonstration purposes.
try(tar_make())
exists("x") # Should be FALSE.
tail(.Random.seed) # for comparison to the RNG state after tar_workspace(y)
tar_workspace(y)
exists("x") # Should be TRUE.
print(x) # "loaded"
# Should be different: tar_workspace() runs set.seed(tar_meta(y, seed)$seed)
tail(.Random.seed)
})
}
}
\seealso{
Other debug: 
\code{\link{tar_load_globals}()},
\code{\link{tar_traceback}()},
\code{\link{tar_workspaces}()}
}
\concept{debug}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_path.R
\name{path_objects_dir_cloud}
\alias{path_objects_dir_cloud}
\title{Default pseudo-directory path of target data in the cloud}
\usage{
path_objects_dir_cloud()
}
\value{
Character of length,
default pseudo-directory path of target data in the cloud.
}
\description{
Not a user-side function. Do not invoke directly.
}
\examples{
path_objects_dir_cloud()
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_exist_meta.R
\name{tar_exist_meta}
\alias{tar_exist_meta}
\title{Check if target metadata exists.}
\usage{
tar_exist_meta(store = targets::tar_config_get("store"))
}
\arguments{
\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Logical of length 1, whether the current project's metadata exists.
}
\description{
Check if the target metadata file \verb{_targets/meta/meta}
exists for the current project.
}
\details{
To learn more about local storage in \code{targets}, visit
\url{https://books.ropensci.org/targets/files.html#internal-files}.
}
\examples{
tar_exist_meta()
}
\seealso{
Other existence: 
\code{\link{tar_exist_objects}()},
\code{\link{tar_exist_process}()},
\code{\link{tar_exist_progress}()},
\code{\link{tar_exist_script}()}
}
\concept{existence}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_newer.R
\name{tar_newer}
\alias{tar_newer}
\title{List new targets}
\usage{
tar_newer(
  time,
  names = NULL,
  inclusive = FALSE,
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{time}{A \code{POSIXct} object of length 1, time threshold.
Targets newer than this time stamp are returned.
For example, if \code{time = Sys.time - as.difftime(1, units = "weeks")}
then \code{tar_newer()} returns targets newer than one week ago.}

\item{names}{Names of eligible targets. Targets excluded from \code{names}
will not be returned even if they are newer than the given \code{time}.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
If \code{NULL}, all names are eligible.}

\item{inclusive}{Logical of length 1, whether to include targets
built at exactly the \code{time} given.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A character vector of names of old targets with recorded
timestamp metadata.
}
\description{
List all the targets whose last successful run occurred
after a certain point in time.
}
\details{
Only applies to targets with recorded time stamps:
just non-branching targets and individual dynamic branches.
As of \code{targets} version 0.6.0, these time
stamps are available for these targets regardless of
storage format. Earlier versions of \code{targets} do not record
time stamps for remote storage formats such as \code{"url"}
or any of the \code{"aws_*"} formats.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(tar_target(x, seq_len(2)))
}, ask = FALSE)
tar_make()
# targets newer than 1 week ago
tar_newer(Sys.time() - as.difftime(1, units = "weeks"))
# targets newer than 1 week from now
tar_newer(Sys.time() + as.difftime(1, units = "weeks"))
# Everything is still up to date.
tar_make()
# Invalidate all targets targets newer than 1 week ago
# so they run on the next tar_make().
invalidate_these <- tar_newer(Sys.time() - as.difftime(1, units = "weeks"))
tar_invalidate(all_of(invalidate_these))
tar_make()
})
}
}
\seealso{
Other time: 
\code{\link{tar_older}()},
\code{\link{tar_timestamp_raw}()},
\code{\link{tar_timestamp}()}
}
\concept{time}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_envir.R
\name{tar_envir}
\alias{tar_envir}
\title{For developers only: get the environment of the current target.}
\usage{
tar_envir(default = parent.frame())
}
\arguments{
\item{default}{Environment, value to return if \code{tar_envir()}
is called on its own outside a \code{targets} pipeline.
Having a default lets users run things without \code{\link[=tar_make]{tar_make()}},
which helps peel back layers of code and troubleshoot bugs.}
}
\value{
If called from a running target, \code{tar_envir()} returns
the environment where the target runs its command.
If called outside a pipeline, the return value is
whatever the user supplies to \code{default}
(which defaults to \code{parent.frame()}).
}
\description{
For developers only: get the environment where a
target runs its command. Designed to be called
while the target is running. The environment
inherits from \code{tar_option_get("envir")}.
}
\details{
Most users should not use \code{tar_envir()} because accidental
modifications to \code{parent.env(tar_envir())} could break the pipeline.
\code{tar_envir()} only exists in order to support third-party interface
packages, and even then the returned environment is not modified.
}
\examples{
tar_envir()
tar_envir(default = new.env(parent = emptyenv()))
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(tar_target(x, tar_envir(default = parent.frame())))
tar_make(x)
tar_read(x)
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_active}()},
\code{\link{tar_call}()},
\code{\link{tar_cancel}()},
\code{\link{tar_definition}()},
\code{\link{tar_group}()},
\code{\link{tar_name}()},
\code{\link{tar_path}()},
\code{\link{tar_seed}()},
\code{\link{tar_store}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_target.R
\name{tar_target}
\alias{tar_target}
\title{Declare a target.}
\usage{
tar_target(
  name,
  command,
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
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[=tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[=tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[=tar_make_future]{tar_make_future()}}.
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
\code{\link[=tar_path]{tar_path()}} from inside the target.
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
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[=tar_make_future]{tar_make_future()}}.
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
A target object. Users should not modify these directly,
just feed them to \code{\link[=list]{list()}} in your target script file
(default: \verb{_targets.R}).
}
\description{
A target is a single step of computation in a pipeline.
It runs an R command and returns a value.
This value gets treated as an R object that can be used
by the commands of targets downstream. Targets that
are already up to date are skipped. See the user manual
for more details.
}
\section{Target objects}{

Functions like \code{tar_target()} produce target objects,
special objects with specialized sets of S3 classes.
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

\section{Storage formats}{

\itemize{
\item \code{"rds"}: Default, uses \code{saveRDS()} and \code{readRDS()}. Should work for
most objects, but slow.
\item \code{"qs"}: Uses \code{qs::qsave()} and \code{qs::qread()}. Should work for
most objects, much faster than \code{"rds"}. Optionally set the
preset for \code{qsave()} through \code{tar_resources()} and \code{tar_resources_qs()}.
\item \code{"feather"}: Uses \code{arrow::write_feather()} and
\code{arrow::read_feather()} (version 2.0). Much faster than \code{"rds"},
but the value must be a data frame. Optionally set
\code{compression} and \code{compression_level} in \code{arrow::write_feather()}
through \code{tar_resources()} and \code{tar_resources_feather()}.
Requires the \code{arrow} package (not installed by default).
\item \code{"parquet"}: Uses \code{arrow::write_parquet()} and
\code{arrow::read_parquet()} (version 2.0). Much faster than \code{"rds"},
but the value must be a data frame. Optionally set
\code{compression} and \code{compression_level} in \code{arrow::write_parquet()}
through \code{tar_resources()} and \code{tar_resources_parquet()}.
Requires the \code{arrow} package (not installed by default).
\item \code{"fst"}: Uses \code{fst::write_fst()} and \code{fst::read_fst()}.
Much faster than \code{"rds"}, but the value must be
a data frame. Optionally set the compression level for
\code{fst::write_fst()} through \code{tar_resources()} and \code{tar_resources_fst()}.
Requires the \code{fst} package (not installed by default).
\item \code{"fst_dt"}: Same as \code{"fst"}, but the value is a \code{data.table}.
Optionally set the compression level the same way as for \code{"fst"}.
\item \code{"fst_tbl"}: Same as \code{"fst"}, but the value is a \code{tibble}.
Optionally set the compression level the same way as for \code{"fst"}.
\item \code{"keras"}: Uses \code{keras::save_model_hdf5()} and
\code{keras::load_model_hdf5()}. The value must be a Keras model.
Requires the \code{keras} package (not installed by default).
\item \code{"torch"}: Uses \code{torch::torch_save()} and \code{torch::torch_load()}.
The value must be an object from the \code{torch} package
such as a tensor or neural network module.
Requires the \code{torch} package (not installed by default).
\item \code{"file"}: A dynamic file. To use this format,
the target needs to manually identify or save some data
and return a character vector of paths
to the data. (These paths must be existing files
and nonempty directories.)
Then, \code{targets} automatically checks those files and cues
the appropriate build decisions if those files are out of date.
Those paths must point to files or directories,
and they must not contain characters \code{|} or \code{*}.
All the files and directories you return must actually exist,
or else \code{targets} will throw an error. (And if \code{storage} is \code{"worker"},
\code{targets} will first stall out trying to wait for the file
to arrive over a network file system.)
If the target does not create any files, the return value should be
\code{character(0)}.
\item \code{"url"}: A dynamic input URL. It works like \code{format = "file"}
except the return value of the target is a URL that already exists
and serves as input data for downstream targets. Optionally
supply a custom \code{curl} handle through
\code{tar_resources()} and \code{tar_resources_url()}.
in \code{new_handle()}, \code{nobody = TRUE} is important because it
ensures \code{targets} just downloads the metadata instead of
the entire data file when it checks time stamps and hashes.
The data file at the URL needs to have an ETag or a Last-Modified
time stamp, or else the target will throw an error because
it cannot track the data. Also, use extreme caution when
trying to use \code{format = "url"} to track uploads. You must be absolutely
certain the ETag and Last-Modified time stamp are fully updated
and available by the time the target's command finishes running.
\code{targets} makes no attempt to wait for the web server.
\item \code{"aws_rds"}, \code{"aws_qs"}, \code{"aws_parquet"}, \code{"aws_fst"}, \code{"aws_fst_dt"},
\code{"aws_fst_tbl"}, \code{"aws_keras"}: versions of the
respective formats \code{"rds"}, \code{"qs"}, etc. powered by
Amazon Web Services (AWS) Simple Storage Service (S3).
The only difference is that the data file is
uploaded to the AWS S3 bucket
you supply to \code{tar_resources_aws()}. See the cloud computing chapter
of the manual for details.
\item \code{"aws_file"}: arbitrary dynamic files on AWS S3. The target
should return a path to a temporary local file, then
\code{targets} will automatically upload this file to an S3
bucket and track it for you. Unlike \code{format = "file"},
\code{format = "aws_file"} can only handle one single file,
and that file must not be a directory.
\code{\link[=tar_read]{tar_read()}} and downstream targets
download the file to \verb{_targets/scratch/} locally and return the path.
\verb{_targets/scratch/} gets deleted at the end of \code{\link[=tar_make]{tar_make()}}.
Requires the same \code{resources} and other configuration details
as the other AWS-powered formats. See the cloud computing
chapter of the manual for details.
\item An entirely custom specification produced by \code{\link[=tar_format]{tar_format()}}.
}
}

\examples{
# Defining targets does not run them.
data <- tar_target(target_name, get_data(), packages = "tidyverse")
analysis <- tar_target(analysis, analyze(x), pattern = map(x))
# Pipelines accept targets.
pipeline <- list(data, analysis)
# Tidy evaluation
tar_option_set(envir = environment())
n_rows <- 30L
data <- tar_target(target_name, get_data(!!n_rows))
print(data)
# Disable tidy evaluation:
data <- tar_target(target_name, get_data(!!n_rows), tidy_eval = FALSE)
print(data)
tar_option_reset()
# In a pipeline:
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(tar_target(x, 1 + 1), ask = FALSE)
tar_make()
tar_read(x)
})
}
}
\seealso{
Other targets: 
\code{\link{tar_cue}()},
\code{\link{tar_format}()},
\code{\link{tar_target_raw}()}
}
\concept{targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_cancel.R
\name{tar_cancel}
\alias{tar_cancel}
\title{Cancel a target mid-build under a custom condition.}
\usage{
tar_cancel(condition = TRUE)
}
\arguments{
\item{condition}{Logical of length 1, whether to cancel the target.}
}
\description{
Cancel a target while its command is running
if a condition is met.
}
\details{
Must be invoked by the target itself. \code{tar_cancel()}
cannot interrupt a target from another process.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(tar_target(x, tar_cancel(1 > 0)))
tar_make() # Should cancel target x.
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_active}()},
\code{\link{tar_call}()},
\code{\link{tar_definition}()},
\code{\link{tar_envir}()},
\code{\link{tar_group}()},
\code{\link{tar_name}()},
\code{\link{tar_path}()},
\code{\link{tar_seed}()},
\code{\link{tar_store}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_format.R
\name{tar_format}
\alias{tar_format}
\title{Define a custom target storage format.}
\usage{
tar_format(
  read = function(path) {     readRDS(path) },
  write = function(object, path) {     saveRDS(object = object, file = path, version =
    3L) },
  marshal = function(object) {     identity(object) },
  unmarshal = function(object) {     identity(object) },
  repository = c("default", "aws")
)
}
\arguments{
\item{read}{A function with a single argument named \code{path}.
This function should read and return the target stored
at the file in the argument. It should have no side effects.
See the "Format functions" section for specific requirements.}

\item{write}{A function with two arguments: \code{object} and \code{path},
in that order. This function should save the R object \code{object}
to the file path at \code{path} and have no other side effects.
The return value does not matter.
See the "Format functions" section for specific requirements.}

\item{marshal}{A function with a single argument named \code{object}.
This function should marshal the R object and return
an in-memory object that can be exported to remote parallel workers.
It should not read or write any persistent files.
See the Marshalling section for details.
See the "Format functions" section for specific requirements.}

\item{unmarshal}{A function with a single argument named \code{object}.
This function should unmarshal the (marshalled) R object and return
an in-memory object that is appropriate and valid for use
on a parallel worker. It should not read or write any persistent files.
See the Marshalling section for details.
See the "Format functions" section for specific requirements.}

\item{repository}{Character of length 1, \code{"default"} for local storage
and \code{"aws"} for storage on Amazon S3. Read
\url{https://books.ropensci.org/targets/storage_amazon.html}
for more on Amazon S3 storage.}
}
\value{
A character string of length 1 encoding the custom format.
You can supply this string directly to the \code{format}
argument of \code{\link[=tar_target]{tar_target()}} or \code{\link[=tar_option_set]{tar_option_set()}}.
}
\description{
Define a custom target storage format for the
\code{format} argument of \code{\link[=tar_target]{tar_target()}} or \code{\link[=tar_option_set]{tar_option_set()}}.
}
\section{Marshalling}{

If an object can only be used in the R session
where it was created, it is called "non-exportable".
Examples of non-exportable R objects are Keras models,
Torch objects, \code{xgboost} matrices, \code{xml2} documents,
\code{rstan} model objects, \code{sparklyr} data objects, and
database connection objects. These objects cannot be
exported to parallel workers (e.g. for \code{\link[=tar_make_future]{tar_make_future()}})
without special treatment. To send an non-exportable
object to a parallel worker, the object must be marshalled:
converted into a form that can be exported safely
(similar to serialization but not always the same).
Then, the worker must unmarshal the object: convert it
into a form that is usable and valid in the current R session.
Arguments \code{marshal} and \code{unmarshal} of \code{tar_format()}
let you control how marshalling and unmarshalling happens.
}

\section{Format functions}{

In \code{tar_format()}, functions like \code{read}, \code{write},
\code{marshal}, and \code{unmarshal} must be perfectly pure
and perfectly self-sufficient.
They must load or namespace all their own packages,
and they must not depend on any custom user-defined
functions or objects in the global environment of your pipeline.
\code{targets} converts each function to and from text,
so it must not rely on any data in the closure.
This disqualifies functions produced by \code{Vectorize()},
for example.
}

\examples{
# The following target is equivalent to
# tar_target(name, command(), format = "keras"):
tar_target(
  name,
  command(),
  format = tar_format(
    read = function(path) {
       keras::load_model_hdf5(path)
    },
    write = function(object, path) {
      keras::save_model_hdf5(object = object, filepath = path)
    },
    marshal = function(object) {
      keras::serialize_model(object)
    },
    unmarshal = function(object) {
      keras::unserialize_model(object)
    },
    repository = "default" # Could be "aws" (same as format = "aws_keras")
  )
)
}
\seealso{
Other targets: 
\code{\link{tar_cue}()},
\code{\link{tar_target_raw}()},
\code{\link{tar_target}()}
}
\concept{targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_make_clustermq.R
\name{tar_make_clustermq}
\alias{tar_make_clustermq}
\title{Run a pipeline of targets in parallel with persistent
\code{clustermq} workers.}
\usage{
tar_make_clustermq(
  names = NULL,
  shortcut = targets::tar_config_get("shortcut"),
  reporter = targets::tar_config_get("reporter_make"),
  workers = targets::tar_config_get("workers"),
  log_worker = FALSE,
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function, reporter),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Names of the targets to build or check. Set to \code{NULL} to
check/build all the targets (default). Otherwise, you can supply
\code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
Because \code{\link[=tar_make]{tar_make()}} and friends run the pipeline in a new R session,
if you pass a character vector to a tidyselect helper, you will need
to evaluate that character vector early with \verb{!!}, e.g.
\code{tar_make(names = all_of(!!your_vector))}.
Applies to ordinary targets (stem) and whole dynamic branching targets
(patterns) but not to individual dynamic branches.}

\item{shortcut}{Logical of length 1, how to interpret the \code{names} argument.
If \code{shortcut} is \code{FALSE} (default) then the function checks
all targets upstream of \code{names} as far back as the dependency graph goes.
\code{shortcut = TRUE} increases speed if there are a lot of
up-to-date targets, but it assumes all the dependencies
are up to date, so please use with caution.
It relies on stored metadata for information about upstream dependencies.
\code{shortcut = TRUE} only works if you set \code{names}.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets run in the pipeline.
Defaults to \code{tar_config_get("reporter_make")}. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"summary"}: print a running total of the number of each targets in
each status category (queued, started, skipped, build, canceled,
or errored). Also show a timestamp (\code{"\%H:\%M \%OS2"} \code{strptime()} format)
of the last time the progress changed and printed to the screen.
\item \code{"timestamp"}: same as the \code{"verbose"} reporter except that each
.message begins with a time stamp.
\item \code{"timestamp_positives"}: same as the \code{"timestamp"} reporter
except without messages for skipped targets.
\item \code{"verbose"}: print messages for individual targets
as they start, finish, or are skipped.
\item \code{"verbose_positives"}: same as the \code{"verbose"} reporter
except without messages for skipped targets.
}}

\item{workers}{Positive integer, number of persistent \code{clustermq} workers
to create.}

\item{log_worker}{Logical, whether to write a log file for each worker.
Same as the \code{log_worker} argument of \code{clustermq::Q()}
and \code{clustermq::workers()}.}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
\code{NULL} except if \code{callr_function = callr::r_bg()}, in which case
a handle to the \code{callr} background process is returned. Either way,
the value is invisibly returned.
}
\description{
This function is like \code{\link[=tar_make]{tar_make()}} except that targets
run in parallel with persistent \code{clustermq} workers. It requires
that you set global options like \code{clustermq.scheduler} and
\code{clustermq.template} inside the target script file
(default: \verb{_targets.R}).
\code{clustermq} is not a strict dependency of \code{targets},
so you must install \code{clustermq} yourself.
}
\details{
To use with a cluster, you will need to set the global options
\code{clustermq.scheduler} and \code{clustermq.template} inside the
target script file (default: \verb{_targets.R}).
To read more about configuring \code{clustermq} for your scheduler, visit
\url{https://mschubert.github.io/clustermq/articles/userguide.html#configuration} # nolint
and navigate to the appropriate link under "Setting up the scheduler".
Wildcards in the template file are filled in with elements from
\code{tar_option_get("resources")}.
}
\examples{
if (!identical(tolower(Sys.info()[["sysname"]]), "windows")) {
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  options(clustermq.scheduler = "multicore") # Does not work on Windows.
  tar_option_set()
  list(tar_target(x, 1 + 1))
}, ask = FALSE)
tar_make_clustermq()
})
}
}
}
\seealso{
Other pipeline: 
\code{\link{tar_make_future}()},
\code{\link{tar_make}()}
}
\concept{pipeline}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_envvars.R
\name{tar_envvars}
\alias{tar_envvars}
\title{Show \code{targets} environment variables.}
\usage{
tar_envvars(unset = "")
}
\arguments{
\item{unset}{Character of length 1, value to return
for any environment variable that is not set.}
}
\value{
A data frame with one row per environment variable
and columns with the name and current value of each.
An unset environment variable will have a value of \code{""}
by default. (Customize with the \code{unset} argument).
}
\description{
Show all the special environment variables
available for customizing \code{targets}.
}
\details{
You can customize the behavior of \code{targets}
with special environment variables. The sections in this help file
describe each environment variable, and the \code{tar_envvars()} function
lists their current values.

If you modify environment variables, please set them
in project-level \code{.Renviron} file so you do not lose your
configuration when you restart your R session.
Modify the project-level \code{.Renviron} file with
\code{usethis::edit_r_environ(scope = "project")}. Restart
your R session after you are done editing.

For targets that run on parallel workers
created by \code{\link[=tar_make_clustermq]{tar_make_clustermq()}} or \code{\link[=tar_make_future]{tar_make_future()}},
only the environment variables listed by \code{\link[=tar_envvars]{tar_envvars()}}
are specifically exported to the targets.
For all other environment variables, you will have to set
the values manually, e.g. a project-level \code{.Renviron} file
(for workers that have access to the local file system).
}
\section{TAR_ASK}{

The \code{TAR_ASK} environment variable accepts values \code{"true"} and \code{"false"}.
If \code{TAR_ASK} is not set, or if it is set to \code{"true"},
then \code{targets} asks permission in a menu
before overwriting certain files, such as the target script file
(default: \verb{_targets.R}) in \code{\link[=tar_script]{tar_script()}}.
If \code{TAR_ASK} is \code{"false"}, then \code{targets} overwrites the old files
with the new ones without asking. Once you are comfortable with
\code{\link[=tar_script]{tar_script()}}, \code{\link[=tar_github_actions]{tar_github_actions()}}, and similar functions,
you can safely set \code{TAR_ASK} to \code{"false"} in either a project-level
or user-level \code{.Renviron} file.
}

\section{TAR_CONFIG}{

The \code{TAR_CONFIG} environment variable controls the file path to the
optional YAML configuration file with project settings.
See the help file of \code{\link[=tar_config_set]{tar_config_set()}} for details.
}

\section{TAR_PROJECT}{

The \code{TAR_PROJECT} environment variable sets the name of project
to set and get settings when working with the YAML configuration file.
See the help file of \code{\link[=tar_config_set]{tar_config_set()}} for details.
}

\section{TAR_WARN}{

The \code{TAR_WARN} environment variable accepts values \code{"true"} and \code{"false"}.
If \code{TAR_WARN} is not set, or if it is set to \code{"true"},
then \code{targets} throws warnings in certain edge cases,
such as target/global name conflicts and dangerous use of
\code{devtools::load_all()}. If \code{TAR_WARN} is \code{"false"}, then \code{targets}
does not throw warnings in these cases.
These warnings can detect potentially serious
issues with your pipeline, so please do not set \code{TAR_WARN}
unless your use case absolutely requires it.
}

\examples{
tar_envvars()
}
\seealso{
Other configuration: 
\code{\link{tar_config_get}()},
\code{\link{tar_config_set}()},
\code{\link{tar_config_unset}()},
\code{\link{tar_option_get}()},
\code{\link{tar_option_reset}()},
\code{\link{tar_option_set}()}
}
\concept{configuration}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_github_actions.R
\name{tar_github_actions}
\alias{tar_github_actions}
\title{Set up GitHub Actions to run a targets pipeline}
\usage{
tar_github_actions(
  path = file.path(".github", "workflows", "targets.yaml"),
  ask = NULL
)
}
\arguments{
\item{path}{Character of length 1, file path to write the GitHub Actions
workflow file.}

\item{ask}{Logical, whether to ask before writing if the workflow file
already exists. If \code{NULL}, defaults to \code{Sys.getenv("TAR_ASK")}.
(Set to \code{"true"} or \code{"false"} with \code{Sys.setenv()}).
If \code{ask} and the \code{TAR_ASK} environment variable are both
indeterminate, defaults to \code{interactive()}.}
}
\value{
Nothing (invisibly). This function writes a GitHub Actions
workflow file as a side effect.
}
\description{
Writes a GitHub Actions workflow file so the pipeline
runs on every push to GitHub. Historical runs accumulate in the
\code{targets-runs} branch, and the latest output is restored before
\code{\link[=tar_make]{tar_make()}} so up-to-date targets do not rerun.
}
\details{
Steps to set up continuous deployment:
\enumerate{
\item Ensure your pipeline stays within the resource limitations of
GitHub Actions and repositories, both for storage and compute.
For storage, you may wish to reduce the burden with
AWS-backed storage formats like \code{"aws_qs"}.
\item Ensure Actions are enabled in your GitHub repository.
You may have to visit the Settings tab.
\item Call \code{targets::tar_renv(extras = character(0))}
to expose hidden package dependencies.
\item Set up \code{renv} for your project (with \code{renv::init()}
or \code{renv::snapshot()}). Details at
\url{https://rstudio.github.io/renv/articles/ci.html}.
\item Commit the \code{renv.lock} file to the \code{main} (recommended)
or \code{master} Git branch.
\item Run \code{tar_github_actions()} to create the workflow file.
Commit this file to \code{main} (recommended) or \code{master} in Git.
\item Push your project to GitHub. Verify that a GitHub Actions
workflow runs and pushes results to \code{targets-runs}.
Subsequent runs will only recompute the outdated targets.
}
}
\examples{
tar_github_actions(tempfile())
}
\seealso{
Other scripts: 
\code{\link{tar_edit}()},
\code{\link{tar_helper_raw}()},
\code{\link{tar_helper}()},
\code{\link{tar_renv}()},
\code{\link{tar_script}()}
}
\concept{scripts}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_config_unset.R
\name{tar_config_unset}
\alias{tar_config_unset}
\title{Unset configuration settings.}
\usage{
tar_config_unset(
  names = character(0),
  config = Sys.getenv("TAR_CONFIG", "_targets.yaml"),
  project = Sys.getenv("TAR_PROJECT", "main")
)
}
\arguments{
\item{names}{Character vector of configuration settings
to delete from the current project.}

\item{config}{Character of length 1, file path of the YAML
configuration file with \code{targets} project settings.
The \code{config} argument specifies which YAML configuration
file that \code{tar_config_get()} reads from or \code{tar_config_set()}
writes to in a single function call.
It does not globally change which configuration file is used
in subsequent function calls. The default file path of the YAML
file is always \verb{_targets.yaml} unless you set another
default path using the \code{TAR_CONFIG} environment variable,
e.g. \code{Sys.setenv(TAR_CONFIG = "custom.yaml")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}

\item{project}{Character of length 1, name of the current
\code{targets} project. Thanks to the \code{config} R package,
\code{targets} YAML configuration files can store multiple
sets of configuration settings, with each set corresponding
to its own project. The \code{project} argument allows you to
set or get a configuration setting for a specific project
for a given call to \code{tar_config_set()} or \code{tar_config_get()}.
The default project is always called \code{"main"}
unless you set another
default project using the \code{TAR_PROJECT} environment variable,
e.g. \code{Sys.setenv(tar_project = "custom")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}
}
\value{
\code{NULL} (invisibly)
}
\description{
Unset (i.e. delete) one or more
custom settings for the current project
from the optional YAML configuration file.
After that, \code{\link[=tar_option_get]{tar_option_get()}} will return the original
default values for those settings for the project.
}
\section{Configuration}{

For several key functions like \code{\link[=tar_make]{tar_make()}}, the
default values of arguments are controlled though
\code{tar_config_get()}. \code{tar_config_get()} retrieves data
from an optional YAML configuration file.
You can control the settings in the YAML
file programmatically with \code{tar_config_set()}.
The default file path of this YAML file is \verb{_targets.yaml}, and you can
set another path globally using the \code{TAR_CONFIG}
environment variable. The YAML file can store configuration
settings for multiple projects, and you can globally
set the default project with the \code{TAR_PROJECT} environment
variable.
The structure of the YAML file
follows rules similar to the \code{config} R package, e.g.
projects can inherit settings from one another using the \code{inherits} field.
Exceptions include:
\enumerate{
\item There is no requirement to have a configuration named \code{"default"}.
\item Other projects do not inherit from the default project` automatically.
\item Not all fields need values because \code{targets} already has defaults.
}

\code{targets} does not actually invoke
the \code{config} package. The implementation in \code{targets}
was written from scratch without viewing or copying any
part of the source code of \code{config}.
}

\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(list(tar_target(x, 1 + 1)))
tar_config_get("store") # "_targets"
store_path <- tempfile()
tar_config_set(store = store_path)
tar_config_get("store") # Shows a temp file.
tar_config_unset("store")
tar_config_get("store") # _targets
})
}
}
\seealso{
Other configuration: 
\code{\link{tar_config_get}()},
\code{\link{tar_config_set}()},
\code{\link{tar_envvars}()},
\code{\link{tar_option_get}()},
\code{\link{tar_option_reset}()},
\code{\link{tar_option_set}()}
}
\concept{configuration}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio_addin_tar_make_bg.R
\name{rstudio_addin_tar_make_bg}
\alias{rstudio_addin_tar_make_bg}
\title{RStudio addin to run \code{\link[=tar_make]{tar_make()}} in the background.}
\usage{
rstudio_addin_tar_make_bg()
}
\description{
For internal use only. Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_noninteractive.R
\name{tar_noninteractive}
\alias{tar_noninteractive}
\title{Run if Target Markdown interactive mode is not on.}
\usage{
tar_noninteractive(code)
}
\arguments{
\item{code}{R code to run if Target Markdown interactive mode
is not turned on.}
}
\value{
If Target Markdown interactive mode is not turned on,
the function returns the result of running the code.
Otherwise, the function invisibly returns \code{NULL}.
}
\description{
In Target Markdown, run the enclosed code
only if interactive mode is not activated. Otherwise,
do not run the code.
}
\details{
Visit <books.ropensci.org/targets/markdown.html>
to learn about Target Markdown and interactive mode.
}
\examples{
tar_noninteractive(message("Not in interactive mode."))
}
\seealso{
Other Target Markdown: 
\code{\link{tar_engine_knitr}()},
\code{\link{tar_interactive}()},
\code{\link{tar_toggle}()}
}
\concept{Target Markdown}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_option_reset.R
\name{tar_option_reset}
\alias{tar_option_reset}
\title{Reset all target options.}
\usage{
tar_option_reset()
}
\value{
\code{NULL} (invisibly).
}
\description{
Reset all target options you previously chose with
\code{\link[=tar_option_set]{tar_option_set()}}. These options are mostly configurable default
arguments to \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_target_raw]{tar_target_raw()}}.
}
\examples{
tar_option_get("format") # default format before we set anything
tar_target(x, 1)$settings$format
tar_option_set(format = "fst_tbl") # new default format
tar_option_get("format")
tar_target(x, 1)$settings$format
tar_option_reset() # reset all options
tar_target(x, 1)$settings$format
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set(cue = tar_cue(mode = "always"))
  tar_option_reset() # Undo option above.
  list(tar_target(x, 1), tar_target(y, 2))
})
tar_make()
tar_make()
})
}
}
\seealso{
Other configuration: 
\code{\link{tar_config_get}()},
\code{\link{tar_config_set}()},
\code{\link{tar_config_unset}()},
\code{\link{tar_envvars}()},
\code{\link{tar_option_get}()},
\code{\link{tar_option_set}()}
}
\concept{configuration}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_target_raw.R
\name{tar_target_raw}
\alias{tar_target_raw}
\title{Define a target using unrefined names and language objects.}
\usage{
tar_target_raw(
  name,
  command,
  pattern = NULL,
  packages = targets::tar_option_get("packages"),
  library = targets::tar_option_get("library"),
  deps = NULL,
  string = NULL,
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

\item{command}{Similar to the \code{command} argument of \code{\link[=tar_target]{tar_target()}} except
the object must already be an expression instead of
informally quoted code.
\code{base::expression()} and \code{base::quote()} can produce such objects.}

\item{pattern}{Similar to the \code{pattern} argument of \code{\link[=tar_target]{tar_target()}}
except the object must already be an expression instead of
informally quoted code.
\code{base::expression()} and \code{base::quote()} can produce such objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{deps}{Optional character vector of the adjacent upstream
dependencies of the target, including targets and global objects.
If \code{NULL}, dependencies are resolved automatically as usual.}

\item{string}{Optional string representation of the command.
Internally, the string gets hashed to check if the command changed
since last run, which helps \code{targets} decide whether the
target is up to date. External interfaces can take control of
\code{string} to ignore changes in certain parts of the command.
If \code{NULL}, the strings is just deparsed from \code{command} (default).}

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
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[=tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[=tar_make_future]{tar_make_future()}}).}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[=tar_make_future]{tar_make_future()}}.
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
\code{\link[=tar_path]{tar_path()}} from inside the target.
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
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[=tar_make_future]{tar_make_future()}}.
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
A target object. Users should not modify these directly,
just feed them to \code{\link[=list]{list()}} in your target script file
(default: \verb{_targets.R}).
See the "Target objects" section for details.
}
\description{
\code{tar_target_raw()} is just like \code{\link[=tar_target]{tar_target()}} except
it avoids non-standard evaluation for the arguments: \code{name}
is a character string, \code{command} and \code{pattern} are language objects,
and there is no \code{tidy_eval} argument. Use \code{tar_target_raw()}
instead of \code{\link[=tar_target]{tar_target()}} if you are creating entire batches
of targets programmatically (metaprogramming, static branching).
}
\section{Target objects}{

Functions like \code{tar_target()} produce target objects,
special objects with specialized sets of S3 classes.
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
  # The following are equivalent.
  y <- tar_target(y, sqrt(x), pattern = map(x))
  y <- tar_target_raw("y", expression(sqrt(x)), expression(map(x)))
  # Programmatically create a chain of interdependent targets
  target_list <- lapply(seq_len(4), function(i) {
    tar_target_raw(
      letters[i + 1],
      substitute(do_something(x), env = list(x = as.symbol(letters[i])))
    )
  })
  print(target_list[[1]])
  print(target_list[[2]])
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(tar_target_raw("x", quote(1 + 1)), ask = FALSE)
tar_make()
tar_read(x)
})
}
}
\seealso{
Other targets: 
\code{\link{tar_cue}()},
\code{\link{tar_format}()},
\code{\link{tar_target}()}
}
\concept{targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_edit.R
\name{tar_edit}
\alias{tar_edit}
\title{Open the target script file for editing.}
\usage{
tar_edit(script = targets::tar_config_get("script"))
}
\arguments{
\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}
}
\description{
Open the target script file for editing.
Requires the \code{usethis} package.
}
\details{
The target script file is an R code file
that defines the pipeline. The default path is \verb{_targets.R},
but the default for the current project
can be configured with \code{\link[=tar_config_set]{tar_config_set()}}.
}
\seealso{
Other scripts: 
\code{\link{tar_github_actions}()},
\code{\link{tar_helper_raw}()},
\code{\link{tar_helper}()},
\code{\link{tar_renv}()},
\code{\link{tar_script}()}
}
\concept{scripts}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_progress_summary.R
\name{tar_progress_summary}
\alias{tar_progress_summary}
\title{Summarize target progress.}
\usage{
tar_progress_summary(
  fields = c("skipped", "started", "built", "errored", "canceled", "since"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{fields}{Optional, names of progress data columns to read.
Set to \code{NULL} to read all fields.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A data frame with one row and the following
optional columns that can be selected with \code{fields}.
(\code{time} is omitted by default.)
\itemize{
\item \code{started}: number of targets that started and did not (yet) finish.
\item \code{built}: number of targets that completed without error or cancellation.
\item \code{errored}: number of targets that threw an error.
\item \code{canceled}: number of canceled targets (see \code{\link[=tar_cancel]{tar_cancel()}}).
\item \code{since}: how long ago progress last changed (\code{Sys.time() - time}).
\item \code{time}: the time when the progress last changed
(modification timestamp of the \verb{_targets/meta/progress} file).
}
}
\description{
Summarize the progress of a run of the pipeline.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, x, pattern = map(x)),
    tar_target(z, stopifnot(y < 1.5), pattern = map(y), error = "continue")
  )
}, ask = FALSE)
try(tar_make())
tar_progress_summary()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_load.R
\name{tar_load}
\alias{tar_load}
\title{Load the values of targets.}
\usage{
tar_load(
  names,
  branches = NULL,
  meta = tar_meta(targets_only = TRUE, store = store),
  strict = TRUE,
  silent = FALSE,
  envir = parent.frame(),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Names of the targets to load. You can supply
symbols, a character vector, or \code{tidyselect} helpers like
\code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}. Names are selected
from the metadata in \verb{_targets/meta}, which may
include errored targets.}

\item{branches}{Integer of indices of the branches to load
for any targets that are patterns.}

\item{meta}{Data frame of metadata from \code{\link[=tar_meta]{tar_meta()}}.
\code{tar_read()} with the default arguments can be inefficient for large
pipelines because all the metadata is stored in a single file.
However, if you call \code{\link[=tar_meta]{tar_meta()}} beforehand and supply it to the \code{meta}
argument, then successive calls to \code{tar_read()} may run much faster.}

\item{strict}{Logical of length 1, whether to error out
if one of the selected targets cannot be loaded.
Set to \code{FALSE} to just load the targets that can be loaded
and skip the others.}

\item{silent}{Logical of length 1. If \code{silent} is \code{FALSE}
and \code{strict} is \code{FALSE}, then a message will be printed
if a target cannot be loaded, but load failures
will not stop other targets from being loaded.}

\item{envir}{Environment to put the loaded targets.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Nothing.
}
\description{
Load the return values of targets into the current environment
(or the environment of your choosing). For a typical target, the return
value lives in a file in \verb{_targets/objects/}. For dynamic files
(i.e. \code{format = "file"}) the paths loaded in place of the values.
}
\section{Limited scope}{

\code{tar_read()} and \code{tar_load()}
are only for exploratory analysis and literate programming,
and \code{tar_read_raw()} and \code{tar_load_raw()} are only
for exploratory analysis. \code{targets} automatically
loads the correct dependencies into memory when the pipeline
is running, so invoking these functions
from inside a target is rarely advisable.
}

\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
}, ask = FALSE)
tar_make()
tar_load(starts_with("y")) # see also all_of()
})
}
}
\seealso{
Other data: 
\code{\link{tar_load_raw}()},
\code{\link{tar_meta}()},
\code{\link{tar_objects}()},
\code{\link{tar_pid}()},
\code{\link{tar_process}()},
\code{\link{tar_read_raw}()},
\code{\link{tar_read}()}
}
\concept{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_package.R
\docType{package}
\name{targets-package}
\alias{targets-package}
\title{targets: Dynamic Function-Oriented Make-Like Declarative Pipelines for R}
\description{
As a pipeline toolkit for Statistics and data science in R,
the \code{targets} package brings together
function-oriented programming and Make-like declarative pipelines.
It analyzes the dependency relationships among the tasks of a workflow,
skips steps that are already up to date, runs the necessary
computations with optional parallel workers, abstracts files as
R objects, and provides tangible evidence that the results match
the underlying code and data. The methodology in this package
borrows from GNU Make (2015, ISBN:978-9881443519)
and \code{drake} (2018, \doi{doi:10.21105/joss.00550}).
}
\seealso{
Other help: 
\code{\link{tar_reprex}()},
\code{\link{use_targets}()}
}
\concept{help}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio_addin_tar_glimpse.R
\name{rstudio_addin_tar_glimpse}
\alias{rstudio_addin_tar_glimpse}
\title{RStudio addin to call \code{\link[=tar_glimpse]{tar_glimpse()}}.}
\usage{
rstudio_addin_tar_glimpse()
}
\description{
For internal use only. Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_timestamp.R
\name{tar_timestamp}
\alias{tar_timestamp}
\title{Get the timestamp(s) of a target.}
\usage{
tar_timestamp(
  name = NULL,
  format = NULL,
  tz = NULL,
  parse = NULL,
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{name}{Symbol, name of the target. If \code{NULL} (default)
then \code{tar_timestamp()} will attempt to return the timestamp
of the target currently running. Must be called inside a target's
command or a supporting function in order to work.}

\item{format}{Deprecated in \code{targets} version 0.6.0 (2021-07-21).}

\item{tz}{Deprecated in \code{targets} version 0.6.0 (2021-07-21).}

\item{parse}{Deprecated in \code{targets} version 0.6.0 (2021-07-21).}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
If the target is not recorded in the metadata
or cannot be parsed correctly, then
\code{tar_timestamp()} returns a \code{POSIXct} object at \verb{1970-01-01 UTC}.
}
\description{
Get the timestamp associated with a target's
last successful run.
}
\details{
\code{tar_timestamp()} checks the metadata in \verb{_targets/meta/meta},
not the actual returned data of the target.
The timestamp depends on the storage format of the target.
If storage is local, e.g. formats like \code{"rds"} and \code{"file"},
then the time stamp is the latest modification time
of the target data files at the time the target
last successfully ran. For non-local formats like
\code{"aws_rds"} and \code{"url"}, then \code{targets} chooses instead
to simply record the time the target last successfully ran.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(tar_target(x, 1))
}, ask = FALSE)
tar_make()
# Get the timestamp.
tar_timestamp(x)
# We can use the timestamp to cancel the target
# if it already ran within the last hour.
# Be sure to set `cue = tar_cue(mode = "always")`
# if you want the target to always check the timestamp.
tar_script({
  list(
  tar_target(
    x,
    tar_cancel((Sys.time() - tar_timestamp()) < 3600),
    cue = tar_cue(mode = "always")
  )
)}, ask = FALSE)
tar_make()
})
}
}
\seealso{
Other time: 
\code{\link{tar_newer}()},
\code{\link{tar_older}()},
\code{\link{tar_timestamp_raw}()}
}
\concept{time}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_branch_names_raw.R
\name{tar_branch_names_raw}
\alias{tar_branch_names_raw}
\title{Branch names (raw version)}
\usage{
tar_branch_names_raw(name, index, store = targets::tar_config_get("store"))
}
\arguments{
\item{name}{Character of length 1,
name of the dynamic branching target (pattern).}

\item{index}{Integer vector of branch indexes.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A character vector of branch names.
}
\description{
Get the branch names of a dynamic branching target
using numeric indexes. Same as \code{\link[=tar_branch_names]{tar_branch_names()}} except
\code{name} is a character of length 1.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(w, 1),
    tar_target(x, seq_len(4)),
    tar_target(y, 2 * x, pattern = map(x)),
    tar_target(z, y, pattern = map(y))
  )
}, ask = FALSE)
tar_make()
tar_branch_names_raw("z", c(2, 3))
})
}
}
\seealso{
Other branching: 
\code{\link{tar_branch_index}()},
\code{\link{tar_branch_names}()},
\code{\link{tar_branches}()},
\code{\link{tar_pattern}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_pipeline.R
\name{tar_pipeline}
\alias{tar_pipeline}
\title{Declare a pipeline (deprecated).}
\usage{
tar_pipeline(...)
}
\arguments{
\item{...}{Targets or lists of targets defined with \code{\link[=tar_target]{tar_target()}}.}
}
\value{
A pipeline object.
}
\description{
Functions \code{tar_pipeline()} and \code{\link[=tar_bind]{tar_bind()}} are deprecated.
Instead, simply end your target script file (default: \verb{_targets.R})
with a list of target objects.
You can nest these objects however you like.
}
\details{
Deprecated on 2021-01-03.
}
\examples{
# In _targets.R:
library(targets)
list( # You no longer need tar_pipeline() here.
  tar_target(data_file, "data.csv", format = "file"),
  list( # Target lists can be arbitrarily nested.
    tar_target(data_object, read.csv(data_file)),
    tar_target(analysis, analyze(data_object))
  )
)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_pid.R
\name{tar_pid}
\alias{tar_pid}
\title{Get main process ID.}
\usage{
tar_pid(store = targets::tar_config_get("store"))
}
\arguments{
\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Integer with the process ID (PID) of the most recent
main R process to orchestrate the targets of the current project.
}
\description{
Get the process ID (PID) of the most recent main R process
to orchestrate the targets of the current project.
}
\details{
The main process is the R process invoked
by \code{\link[=tar_make]{tar_make()}} or similar. If \code{callr_function} is not \code{NULL},
this is an external process, and the \code{pid} in the return value
will not agree with \code{Sys.getpid()} in your current interactive session.
The process may or may not be alive. You may want to
check it with \code{ps::ps_is_running(ps::ps_handle(targets::tar_pid()))}
before running another call to \code{\link[=tar_make]{tar_make()}}
for the same project.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
Sys.getpid()
tar_pid() # Different from the current PID.
})
}
}
\seealso{
Other data: 
\code{\link{tar_load_raw}()},
\code{\link{tar_load}()},
\code{\link{tar_meta}()},
\code{\link{tar_objects}()},
\code{\link{tar_process}()},
\code{\link{tar_read_raw}()},
\code{\link{tar_read}()}
}
\concept{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_poll.R
\name{tar_poll}
\alias{tar_poll}
\title{Repeatedly poll progress in the R console.}
\usage{
tar_poll(
  interval = 1,
  timeout = Inf,
  fields = c("skipped", "started", "built", "errored", "canceled", "since"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{interval}{Number of seconds to wait between iterations
of polling progress.}

\item{timeout}{How many seconds to run before exiting.}

\item{fields}{Optional, names of progress data columns to read.
Set to \code{NULL} to read all fields.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\description{
Print the information in \code{\link[=tar_progress_summary]{tar_progress_summary()}}
at regular intervals.
}
\examples{
if (identical(Sys.getenv("TAR_INTERACTIVE_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(100)),
    tar_target(y, Sys.sleep(0.1), pattern = map(x))
  )
}, ask = FALSE)
px <- tar_make(callr_function = callr::r_bg, reporter = "silent")
tar_poll()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_group.R
\name{tar_group}
\alias{tar_group}
\title{Group a data frame to iterate over subsets of rows.}
\usage{
tar_group(x)
}
\arguments{
\item{x}{Grouped data frame from \code{dplyr::group_by()}}
}
\value{
A data frame with a special \code{tar_group} column that
\code{targets} will use to find subsets of your data frame.
}
\description{
Like \code{dplyr::group_by()}, but for patterns.
\code{tar_group()} allows you to map or cross over subsets of data frames.
Requires \code{iteration = "group"} on the target. See the example.
}
\details{
The goal of \code{tar_group()} is to post-process the return value
of a data frame target to allow downstream targets to branch over
subsets of rows. It takes the groups defined by \code{dplyr::group_by()}
and translates that information into a special \code{tar_group} is a column.
\code{tar_group} is a vector of positive integers
from 1 to the number of groups. Rows with the same integer in \code{tar_group}
belong to the same group, and branches are arranged in increasing order
with respect to the integers in \code{tar_group}.
The assignment of \code{tar_group} integers to group levels
depends on the orderings inside the grouping variables and not the order
of rows in the dataset. \code{dplyr::group_keys()} on the grouped data frame
shows how the grouping variables correspond to the integers in the
\code{tar_group} column.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
# The tar_group() function simply creates
# a tar_group column to partition the rows
# of a data frame.
data.frame(
  x = seq_len(6),
  id = rep(letters[seq_len(3)], each = 2)
) \%>\%
  dplyr::group_by(id) \%>\%
  tar_group()
# We use tar_group() below to branch over
# subsets of a data frame defined with dplyr::group_by().
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
library(dplyr)
list(
  tar_target(
    data,
    data.frame(
      x = seq_len(6),
      id = rep(letters[seq_len(3)], each = 2)
    ) \%>\%
      group_by(id) \%>\%
      tar_group(),
    iteration = "group"
  ),
  tar_target(
    sums,
    sum(data$x),
    pattern = map(data),
    iteration = "vector"
  )
)
})
tar_make()
tar_read(sums) # Should be c(3, 7, 11).
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_active}()},
\code{\link{tar_call}()},
\code{\link{tar_cancel}()},
\code{\link{tar_definition}()},
\code{\link{tar_envir}()},
\code{\link{tar_name}()},
\code{\link{tar_path}()},
\code{\link{tar_seed}()},
\code{\link{tar_store}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_outdated.R
\name{tar_outdated}
\alias{tar_outdated}
\title{Check which targets are outdated.}
\usage{
tar_outdated(
  names = NULL,
  shortcut = targets::tar_config_get("shortcut"),
  branches = FALSE,
  targets_only = TRUE,
  reporter = targets::tar_config_get("reporter_outdated"),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function, reporter),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Names of the targets. \code{tar_outdated()} will check
these targets and all upstream ancestors in the dependency graph.
Set \code{names} to \code{NULL} to check/build all the targets (default).
Otherwise, you can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
Applies to ordinary targets (stem) and whole dynamic branching targets
(patterns) but not to individual dynamic branches.}

\item{shortcut}{Logical of length 1, how to interpret the \code{names} argument.
If \code{shortcut} is \code{FALSE} (default) then the function checks
all targets upstream of \code{names} as far back as the dependency graph goes.
If \code{TRUE}, then the function only checks the targets in \code{names}
and uses stored metadata for information about upstream dependencies
as needed. \code{shortcut = TRUE} increases speed if there are a lot of
up-to-date targets, but it assumes all the dependencies
are up to date, so please use with caution.
Also, \code{shortcut = TRUE} only works if you set \code{names}.}

\item{branches}{Logical of length 1, whether to include branch names.
Including branches could get cumbersome for large pipelines.
Individual branch names are still omitted when branch-specific information
is not reliable: for example, when a pattern branches over
an outdated target.}

\item{targets_only}{Logical of length 1, whether to just restrict to targets
or to include functions and other global objects from the environment
created by running the target script file (default: \verb{_targets.R}).}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets are checked. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"forecast"}: print running totals of the checked and outdated
targets found so far.
}}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Names of the outdated targets.
}
\description{
Checks for outdated targets in the pipeline,
targets that will be rerun automatically if you call
\code{\link[=tar_make]{tar_make()}} or similar. See \code{\link[=tar_cue]{tar_cue()}} for the rules
that decide whether a target needs to rerun.
}
\details{
Requires that you define a pipeline
with a target script file (default: \verb{_targets.R}).
(See \code{\link[=tar_script]{tar_script()}} for details.)
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(list(tar_target(x, 1 + 1)))
tar_outdated()
tar_script({
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
}, ask = FALSE)
tar_outdated()
})
}
}
\seealso{
Other inspect: 
\code{\link{tar_deps_raw}()},
\code{\link{tar_deps}()},
\code{\link{tar_glimpse}()},
\code{\link{tar_manifest}()},
\code{\link{tar_network}()},
\code{\link{tar_sitrep}()},
\code{\link{tar_validate}()},
\code{\link{tar_visnetwork}()}
}
\concept{inspect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources_parquet.R
\name{tar_resources_parquet}
\alias{tar_resources_parquet}
\title{Target resources: parquet storage formats}
\usage{
tar_resources_parquet(compression = "snappy", compression_level = NULL)
}
\arguments{
\item{compression}{Character of length 1, \code{compression}
argument of \code{arrow::write_parquet()}.}

\item{compression_level}{Numeric of length 1, \code{compression_level}
argument of \code{arrow::write_parquet()}.}
}
\value{
Object of class \code{"tar_resources_parquet"}, to be supplied
to the parquet argument of \code{tar_resources()}.
}
\description{
Create the \code{parquet} argument of \verb{tar_resources()`` to specify optional settings for parquet data frame storage formats powered by the }arrow\verb{R package. See the}format` argument of \code{\link[=tar_target]{tar_target()}} for details.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
# Somewhere in you target script file (usually _targets.R):
tar_target(
  name,
  command(),
  format = "parquet",
  resources = tar_resources(
    parquet = tar_resources_parquet(compression = "lz4")
  )
)
}
\seealso{
Other resources: 
\code{\link{tar_resources_aws}()},
\code{\link{tar_resources_clustermq}()},
\code{\link{tar_resources_feather}()},
\code{\link{tar_resources_fst}()},
\code{\link{tar_resources_future}()},
\code{\link{tar_resources_gcp}()},
\code{\link{tar_resources_qs}()},
\code{\link{tar_resources_url}()},
\code{\link{tar_resources}()}
}
\concept{resources}
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
% Please edit documentation in R/rstudio_addin_tar_outdated.R
\name{rstudio_addin_tar_outdated}
\alias{rstudio_addin_tar_outdated}
\title{RStudio addin to call \code{\link[=tar_outdated]{tar_outdated()}}.}
\usage{
rstudio_addin_tar_outdated()
}
\description{
For internal use only. Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_active.R
\name{tar_active}
\alias{tar_active}
\title{Show if the pipeline is running.}
\usage{
tar_active()
}
\value{
Logical of length 1, \code{TRUE} if called in a target or \verb{_targets.R}
and the pipeline is running (\code{FALSE} otherwise).
}
\description{
Return \code{TRUE} if called in a target or \verb{_targets.R} and
the pipeline is running.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_active() # FALSE
tar_script({
  message("Pipeline running? ", tar_active())
  tar_target(x, tar_active())
})
tar_manifest() # prints "Pipeline running? FALSE"
tar_make() # prints "pipeline running? TRUE"
tar_read(x) # TRUE
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_call}()},
\code{\link{tar_cancel}()},
\code{\link{tar_definition}()},
\code{\link{tar_envir}()},
\code{\link{tar_group}()},
\code{\link{tar_name}()},
\code{\link{tar_path}()},
\code{\link{tar_seed}()},
\code{\link{tar_store}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_option_get.R
\name{tar_option_get}
\alias{tar_option_get}
\title{Get a target option.}
\usage{
tar_option_get(name = NULL, option = NULL)
}
\arguments{
\item{name}{Character of length 1, name of an option to get.
Must be one of the argument names of \code{\link[=tar_option_set]{tar_option_set()}}.}

\item{option}{Deprecated, use the \code{name} argument instead.}
}
\value{
Value of a target option.
}
\description{
Get a target option. These options include default arguments to
\code{\link[=tar_target]{tar_target()}} such as packages, storage format,
iteration type, and cue.
Needs to be called before any calls to \code{\link[=tar_target]{tar_target()}}
in order to take effect.
}
\details{
This function goes well with \code{\link[=tar_target_raw]{tar_target_raw()}} when it comes
to defining external interfaces on top of the \code{targets} package to create
pipelines.
}
\examples{
tar_option_get("format") # default format before we set anything
tar_target(x, 1)$settings$format
tar_option_set(format = "fst_tbl") # new default format
tar_option_get("format")
tar_target(x, 1)$settings$format
tar_option_reset() # reset the format
tar_target(x, 1)$settings$format
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set(cue = tar_cue(mode = "always")) # All targets always run.
  list(tar_target(x, 1), tar_target(y, 2))
})
tar_make()
tar_make()
})
}
}
\seealso{
Other configuration: 
\code{\link{tar_config_get}()},
\code{\link{tar_config_set}()},
\code{\link{tar_config_unset}()},
\code{\link{tar_envvars}()},
\code{\link{tar_option_reset}()},
\code{\link{tar_option_set}()}
}
\concept{configuration}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources_clustermq.R
\name{tar_resources_clustermq}
\alias{tar_resources_clustermq}
\title{Target resources: \code{clustermq} high-performance computing}
\usage{
tar_resources_clustermq(template = list())
}
\arguments{
\item{template}{Named list, \code{template} argument to
\code{clustermq::workers()}.}
}
\value{
Object of class \code{"tar_resources_clustermq"}, to be supplied
to the \code{clustermq} argument of \code{tar_resources()}.
}
\description{
Create the \code{clustermq} argument of \code{tar_resources()}
to specify optional high-performance computing settings
for \code{tar_make_clustermq()}.
For details, see the documentation of the \code{clustermq} R package
and the corresponding argument names in this help file.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
# Somewhere in you target script file (usually _targets.R):
tar_target(
  name,
  command(),
  resources = tar_resources(
    clustermq = tar_resources_clustermq(template = list(n_cores = 2))
  )
)
}
\seealso{
Other resources: 
\code{\link{tar_resources_aws}()},
\code{\link{tar_resources_feather}()},
\code{\link{tar_resources_fst}()},
\code{\link{tar_resources_future}()},
\code{\link{tar_resources_gcp}()},
\code{\link{tar_resources_parquet}()},
\code{\link{tar_resources_qs}()},
\code{\link{tar_resources_url}()},
\code{\link{tar_resources}()}
}
\concept{resources}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_test.R
\name{tar_test}
\alias{tar_test}
\title{Test code in a temporary directory.}
\usage{
tar_test(label, code)
}
\arguments{
\item{label}{Character of length 1, label for the test.}

\item{code}{User-defined code for the test.}
}
\value{
\code{NULL} (invisibly).
}
\description{
Runs a \code{test_that()} unit test inside a temporary
directory to avoid writing to the user's file space.
This helps ensure compliance with CRAN policies.
Also isolates \code{tar_option_set()}
options and environment variables specific to \code{targets}
and skips the test on Solaris.
Useful for writing tests for
\href{https://wlandau.github.io/targetopia/}{targetopia} packages
(extensions to \code{targets} tailored to specific use cases).
}
\examples{
tar_test("example test", {
  testing_variable_cafecfcb <- "only defined inside tar_test()"
  file.create("only_exists_in_tar_test")
})
exists("testing_variable_cafecfcb")
file.exists("only_exists_in_tar_test")
}
\seealso{
Other utilities to extend targets: 
\code{\link{tar_assert}},
\code{\link{tar_condition}},
\code{\link{tar_dir}()},
\code{\link{tar_language}}
}
\concept{utilities to extend targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_bind.R
\name{tar_bind}
\alias{tar_bind}
\title{Combine pipeline objects (deprecated).}
\usage{
tar_bind(...)
}
\arguments{
\item{...}{Pipeline objects or nested lists of pipeline objects.
You can generate a pipeline object with \code{\link[=tar_pipeline]{tar_pipeline()}}.}
}
\description{
Functions \code{tar_pipeline()} and \code{\link[=tar_bind]{tar_bind()}} are deprecated.
Instead, simply end your target script file
(default: \verb{_targets.R}) file with a list of target objects.
You can nest these objects however you like.
}
\details{
Deprecated on 2021-01-03.
}
\examples{
# In your target script file (default: _targets.R):
library(targets)
list( # You no longer need tar_pipeline() here.
  tar_target(data_file, "data.csv", format = "file"),
  list( # Target lists can be arbitrarily nested.
    tar_target(data_object, read.csv(data_file)),
    tar_target(analysis, analyze(data_object))
  )
)
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_runtime.R
\name{tar_runtime_object}
\alias{tar_runtime_object}
\title{Get the \code{tar_runtime} object.}
\usage{
tar_runtime_object()
}
\value{
The internal \code{tar_runtime} object of class \code{"tar_runtime"}.
}
\description{
For internal purposes only. Not a user-side function.
Do not invoke directly.
}
\details{
Manages internal settings
that targets need while they run.
}
\examples{
tar_runtime_object()
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_objects.R
\name{tar_objects}
\alias{tar_objects}
\title{List saved targets}
\usage{
tar_objects(names = NULL, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Optional \code{tidyselect} selector such as
\code{\link[=all_of]{all_of()}} or \code{\link[=starts_with]{starts_with()}} to return
a tactical subset of target names.
If \code{NULL}, all names are selected.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Character vector of targets saved to \verb{_targets/objects/}.
}
\description{
List targets currently saved to \verb{_targets/objects/}.
Does not include dynamic files or cloud storage.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(tar_target(x, "value"))
}, ask = FALSE)
tar_make()
tar_objects()
tar_objects(starts_with("x")) # see also all_of()
})
}
}
\seealso{
Other data: 
\code{\link{tar_load_raw}()},
\code{\link{tar_load}()},
\code{\link{tar_meta}()},
\code{\link{tar_pid}()},
\code{\link{tar_process}()},
\code{\link{tar_read_raw}()},
\code{\link{tar_read}()}
}
\concept{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_callr.R
\name{callr_args_default}
\alias{callr_args_default}
\title{Default \code{callr} arguments.}
\usage{
callr_args_default(callr_function, reporter = NULL)
}
\arguments{
\item{callr_function}{A function from the \code{callr} package
that starts an external R process.}

\item{reporter}{Character of length 1, choice of reporter
for \code{\link[=tar_make]{tar_make()}} or a related function.}
}
\value{
A list of arguments to \code{callr_function}.
}
\description{
Default \code{callr} arguments for the \code{callr_arguments}
argument of \code{\link[=tar_make]{tar_make()}} and related functions.
}
\details{
Not a user-side function. Do not invoke directly.
Exported for internal purposes only.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_helper.R
\name{tar_helper}
\alias{tar_helper}
\title{Write a helper R script.}
\usage{
tar_helper(path = NULL, code = NULL, tidy_eval = TRUE, envir = parent.frame())
}
\arguments{
\item{path}{Character of length 1, path to write (or overwrite) \code{code}.
If the parent directory does not exist, \code{tar_helper_raw()} creates it.}

\item{code}{Quoted code to write to \code{path}.
\code{tar_helper()} overwrites the file if it already exists.}

\item{tidy_eval}{Logical, whether to use tidy evaluation on \code{code}. If
turned on, you can substitute expressions and symbols using \verb{!!} and \verb{!!!}.
See examples below.}

\item{envir}{Environment for tidy evaluation.}
}
\value{
\code{NULL} (invisibly)
}
\description{
Write a helper R script for a \code{targets} pipeline.
Could be supporting functions or the target script file
(default: \verb{_targets.R}) itself.
}
\details{
\code{tar_helper()} is a specialized version of \code{\link[=tar_script]{tar_script()}}
with flexible paths and tidy evaluation.
}
\examples{
# Without tidy evaluation:
path <- tempfile()
tar_helper(path, x <- 1)
writeLines(readLines(path))
# With tidy evaluation:
y <- 123
tar_helper(path, x <- !!y)
writeLines(readLines(path))
}
\seealso{
Other scripts: 
\code{\link{tar_edit}()},
\code{\link{tar_github_actions}()},
\code{\link{tar_helper_raw}()},
\code{\link{tar_renv}()},
\code{\link{tar_script}()}
}
\concept{scripts}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_path.R
\name{tar_path}
\alias{tar_path}
\title{Identify the file path where a target will be stored.}
\usage{
tar_path(
  name = NULL,
  default = NA_character_,
  create_dir = FALSE,
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{name}{Symbol, name of a target.
If \code{NULL}, \code{tar_path()} returns the path of the target currently running
in a pipeline.}

\item{default}{Character, value to return if \code{tar_path()}
is called on its own outside a \code{targets} pipeline.
Having a default lets users run things without \code{\link[=tar_make]{tar_make()}},
which helps peel back layers of code and troubleshoot bugs.}

\item{create_dir}{Logical of length 1,
whether to create \code{dirname(tar_path())} in \code{tar_path()} itself.
This is useful if you are writing to \code{tar_path()} from inside a
\code{storage = "none"} target and need the parent directory of the file
to exist.}

\item{store}{Character of length 1, path to the data store if \code{tar_path()}
is called outside a running pipeline. If \code{tar_path()} is called
inside a running pipeline, this argument is ignored
and actual the path to the running pipeline's data store
is used instead.}
}
\value{
Character, file path of the return value of the target.
If not called from inside a running target,
\code{tar_path(name = your_target)} just returns
\verb{_targets/objects/your_target}, the file path where \code{your_target}
will be saved unless \code{format} is equal to \code{"file"} or any of the
supported cloud-based storage formats.

For non-cloud storage formats, if you call \code{tar_path()}
with no arguments while target \code{x} is running, the \code{name}
argument defaults to the name of the running target,
so \code{tar_path()} returns \verb{_targets/objects/x}.

For cloud-backed formats, \code{tar_path()} returns the
path to the staging file in \verb{_targets/scratch/}.
That way, even if you select a cloud format
(e.g. \code{tar_target(..., format = "aws_parquet", storage = "none")})
then you can still manually write to \code{tar_path(create_dir = TRUE)}
and the \code{targets} package will automatically hash it and
upload it to the AWS S3 bucket. This does not apply to formats
\code{"file"} or \code{"aws_file"}, where you would never need \code{storage = "none"}
anyway.
}
\description{
Identify the file path where a target will be stored
after the target finishes running in the pipeline.
}
\examples{
tar_path()
tar_path(your_target)
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(tar_target(returns_path, tar_path()), ask = FALSE)
tar_make()
tar_read(returns_path)
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_active}()},
\code{\link{tar_call}()},
\code{\link{tar_cancel}()},
\code{\link{tar_definition}()},
\code{\link{tar_envir}()},
\code{\link{tar_group}()},
\code{\link{tar_name}()},
\code{\link{tar_seed}()},
\code{\link{tar_store}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_meta.R
\name{tar_meta}
\alias{tar_meta}
\title{Read a project's metadata.}
\usage{
tar_meta(
  names = NULL,
  fields = NULL,
  targets_only = FALSE,
  complete_only = FALSE,
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Optional, names of the targets. If supplied, \code{tar_meta()}
only returns metadata on these targets.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
If \code{NULL}, all names are selected.}

\item{fields}{Optional, names of columns/fields to select. If supplied,
\code{tar_meta()} only returns the selected metadata columns.
If \code{NULL}, all fields are selected.
You can supply symbols or \code{tidyselect} helpers
like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
The \code{name} column is always included first
no matter what you select. Choices:
\itemize{
\item \code{name}: name of the target or global object.
\item \code{type}: type of the object: either \code{"function"} or \code{"object"}
for global objects, and \code{"stem"}, \code{"branch"},
\code{"map"}, or \code{"cross"} for targets.
\item \code{data}: hash of the output data.
\item \code{command}: hash of the target's deparsed command.
\item \code{depend}: hash of the immediate upstream dependencies of the target.
\item \code{seed}: random number generator seed with which the target was built.
A target's random number generator seed
is a deterministic function of its name. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can recover the seed of a completed target
with \code{tar_meta(your_target, seed)} and run \code{set.seed()}
on the result to locally recreate the target's initial RNG state.
\item \code{path}: A list column of paths to target data. Usually, each element
is a single path, but there could be multiple paths per target
for dynamic files (i.e. \code{tar_target(format = "file")}).
\item \code{time}: \code{POSIXct} object with the time the target's data in storage
was last modified. If the target stores no local file,
then the time stamp corresponds to the time the target last
ran successfully. Only targets that run commands have time stamps:
just non-branching targets and individual dynamic branches.
Displayed in the current time zone of the system.
If there are multiple outputs for that target, as with file targets,
then the maximum time is shown.
\item \code{size}: hash of the sum of all the bytes of the files at \code{path}.
\item \code{bytes}: total file size in bytes of all files in \code{path}.
\item \code{format}: character, one of the admissible data storage formats.
See the \code{format} argument in the \code{\link[=tar_target]{tar_target()}} help file for details.
\item \code{iteration}: character, either \code{"list"} or \code{"vector"}
to describe the iteration and aggregation mode of the target. See the
\code{iteration} argument in the \code{\link[=tar_target]{tar_target()}} help file for details.
\item \code{parent}: for branches, name of the parent pattern.
\item \code{children}: list column, names of the children of targets that
have them. These include buds of stems and branches of patterns.
\item \code{seconds}: number of seconds it took to run the target.
\item \code{warnings}: character string of warning messages
from the last run of the target.
\item \code{error}: character string of the error message if the target errored.
}}

\item{targets_only}{Logical, whether to just show information about targets
or also return metadata on functions and other global objects.}

\item{complete_only}{Logical, whether to return only complete rows
(no \code{NA} values).}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A data frame with one row per target/object and the selected fields.
}
\description{
Read the metadata of all recorded targets and global objects.
}
\details{
A metadata row only updates when the target is built.
\code{\link[=tar_progress]{tar_progress()}} shows information on targets that are running.
That is why the number of branches may disagree between \code{\link[=tar_meta]{tar_meta()}}
and \code{\link[=tar_progress]{tar_progress()}} for actively running pipelines.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
tar_meta()
tar_meta(starts_with("y_")) # see also all_of()
})
}
}
\seealso{
Other data: 
\code{\link{tar_load_raw}()},
\code{\link{tar_load}()},
\code{\link{tar_objects}()},
\code{\link{tar_pid}()},
\code{\link{tar_process}()},
\code{\link{tar_read_raw}()},
\code{\link{tar_read}()}
}
\concept{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_pattern.R
\name{tar_pattern}
\alias{tar_pattern}
\alias{map}
\alias{cross}
\alias{head}
\alias{tail}
\alias{sample}
\title{Emulate dynamic branching.}
\usage{
tar_pattern(pattern, ..., seed = 0L)
}
\arguments{
\item{pattern}{Function call with the pattern specification.}

\item{...}{Named integers, each of length 1.
Each name is the name of a dependency target,
and each integer is the length of the target
(number of branches or slices). Names must be unique.}

\item{seed}{Integer of length 1, random number generator seed to
emulate the pattern reproducibly. (The \code{sample()} pattern is random).
In a real pipeline, the seed is automatically generated
from the target name in deterministic fashion.}
}
\value{
A \code{tibble} showing the kinds of dynamic branches that
\code{\link[=tar_target]{tar_target()}} would create in a real pipeline with the given \code{pattern}.
Each row is a dynamic branch, each column is a dependency target,
and each element is the name of an upstream bud or branch that the
downstream branch depends on. Buds are pieces of non-branching targets
("stems") and branches are pieces of patterns. The returned bud and branch
names are not the actual ones you will see when you run the pipeline,
but they do communicate the branching structure of the pattern.
}
\description{
Emulate the dynamic branching process outside a pipeline.
\code{tar_pattern()} can help you understand the overall branching structure
that comes from the \code{pattern} argument of \code{\link[=tar_target]{tar_target()}}.
}
\details{
Dynamic branching is a way to programmatically
create multiple new targets based on the values of other targets,
all while the pipeline is running. Use the \code{pattern} argument of
\code{\link[=tar_target]{tar_target()}} to get started. \code{pattern} accepts a function call
composed of target names and any of the following patterns:
\itemize{
\item \code{map()}: iterate over one or more targets in sequence.
\item \code{cross()}: iterate over combinations of slices of targets.
\item \code{slice()}: select one or more slices by index, e.g.
\code{slice(x, index = c(3, 4))} selects the third and fourth
slice or branch of \code{x}.
\item \code{head()}: restrict branching to the first few elements.
\item \code{tail()}: restrict branching to the last few elements.
\item \code{sample()}: restrict branching to a random subset of elements.
}
}
\examples{
# To use dynamic map for real in a pipeline,
# call map() in a target's pattern.
# The following code goes at the bottom of
# your target script file (default: `_targets.R`).
list(
  tar_target(x, seq_len(2)),
  tar_target(y, head(letters, 2)),
  tar_target(dynamic, c(x, y), pattern = map(x, y)) # 2 branches
)
# Likewise for more complicated patterns.
list(
  tar_target(x, seq_len(2)),
  tar_target(y, head(letters, 2)),
  tar_target(z, head(LETTERS, 2)),
  tar_target(dynamic, c(x, y, z), pattern = cross(z, map(x, y))) #4 branches
)
# But you can emulate dynamic branching without running a pipeline
# in order to understand the patterns you are creating. Simply supply
# the pattern and the length of each dependency target.
# The returned data frame represents the branching structure of the pattern:
# One row per new branch, one column per dependency target, and
# one element per bud/branch in each dependency target.
tar_pattern(
  cross(x, map(y, z)),
  x = 2,
  y = 3,
  z = 3
)
tar_pattern(
  head(cross(x, map(y, z)), n = 2),
  x = 2,
  y = 3,
  z = 3
)
}
\seealso{
Other branching: 
\code{\link{tar_branch_index}()},
\code{\link{tar_branch_names_raw}()},
\code{\link{tar_branch_names}()},
\code{\link{tar_branches}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_targets.R
\name{use_targets}
\alias{use_targets}
\title{Use targets}
\usage{
use_targets(path = "_targets.Rmd", open = interactive())
}
\arguments{
\item{path}{Character of length 1, output path of the
Target Markdown report relative to the current active
project.}

\item{open}{Logical, whether to open the file for editing
in the RStudio IDE.}
}
\value{
\code{NULL} (invisibly).
}
\description{
Create an example Target Markdown report
to get started with {targets}.
}
\examples{
if (identical(Sys.getenv("TAR_INTERACTIVE_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
use_targets(open = FALSE)
})
}
}
\seealso{
Other help: 
\code{\link{tar_reprex}()},
\code{\link{targets-package}}
}
\concept{help}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_condition.R
\name{tar_condition}
\alias{tar_condition}
\alias{tar_message_run}
\alias{tar_throw_file}
\alias{tar_throw_run}
\alias{tar_throw_validate}
\alias{tar_warn_deprecate}
\alias{tar_warn_run}
\alias{tar_warn_validate}
\title{Conditions}
\usage{
tar_message_run(...)

tar_throw_file(...)

tar_throw_run(...)

tar_throw_validate(...)

tar_warn_deprecate(...)

tar_warn_run(...)

tar_warn_validate(...)
}
\arguments{
\item{...}{zero or more objects which can be coerced to character
    (and which are pasted together with no separator) or a single
    condition object.}
}
\description{
These functions throw custom
\code{targets}-specific error conditions.
Useful for error handling in packages built on top of \code{targets}.
}
\examples{
try(tar_throw_validate("something is not valid"))
}
\seealso{
Other utilities to extend targets: 
\code{\link{tar_assert}},
\code{\link{tar_dir}()},
\code{\link{tar_language}},
\code{\link{tar_test}()}
}
\concept{utilities to extend targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_deduplicate.R
\name{tar_deduplicate}
\alias{tar_deduplicate}
\title{Deduplicate meta and progress databases (deprecated).}
\usage{
tar_deduplicate(
  meta = TRUE,
  progress = TRUE,
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{meta}{Logical, whether to deduplicate the meta database file
at \verb{_targets/meta/meta}.}

\item{progress}{Logical, whether to deduplicate the progress database file
at \verb{_targets/meta/progress}.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Nothing.
}
\description{
Deprecated in version 0.3.0 (2020-03-06).
Deduplication happens automatically before and after the pipeline runs.
}
\details{
Removes duplicated entries in the meta and progress
databases in order to lighten storage. These databases are located
in the \verb{_targets/meta/meta} and \verb{_targets/meta/progress} files,
where \verb{_targets} is the a folder at the project root.
No essential data is removed, so
this is simply a form of garbage collection.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_load_raw.R
\name{tar_load_raw}
\alias{tar_load_raw}
\title{Load the values of targets (raw version).}
\usage{
tar_load_raw(
  names,
  branches = NULL,
  meta = tar_meta(store = store),
  strict = TRUE,
  silent = FALSE,
  envir = parent.frame(),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Character vector, names of the targets to load.
Names are expected to appear in the metadata in \verb{_targets/meta}.}

\item{branches}{Integer of indices of the branches to load
for any targets that are patterns.}

\item{meta}{Data frame of metadata from \code{\link[=tar_meta]{tar_meta()}}.
\code{tar_read()} with the default arguments can be inefficient for large
pipelines because all the metadata is stored in a single file.
However, if you call \code{\link[=tar_meta]{tar_meta()}} beforehand and supply it to the \code{meta}
argument, then successive calls to \code{tar_read()} may run much faster.}

\item{strict}{Logical of length 1, whether to error out
if one of the selected targets cannot be loaded.
Set to \code{FALSE} to just load the targets that can be loaded
and skip the others.}

\item{silent}{Logical of length 1. If \code{silent} is \code{FALSE}
and \code{strict} is \code{FALSE}, then a message will be printed
if a target cannot be loaded, but load failures
will not stop other targets from being loaded.}

\item{envir}{Environment to put the loaded targets.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Nothing.
}
\description{
Same as \code{\link[=tar_load]{tar_load()}} except \code{names} is a character vector.
Do not use in \code{knitr} or R Markdown reports with \code{tarchetypes::tar_knit()}
or \code{tarchetypes::tar_render()}.
}
\section{Limited scope}{

\code{tar_read()} and \code{tar_load()}
are only for exploratory analysis and literate programming,
and \code{tar_read_raw()} and \code{tar_load_raw()} are only
for exploratory analysis. \code{targets} automatically
loads the correct dependencies into memory when the pipeline
is running, so invoking these functions
from inside a target is rarely advisable.
}

\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
}, ask = FALSE)
tar_make()
tar_load_raw(c("y1", "y2"))
y1
y2
})
}
}
\seealso{
Other data: 
\code{\link{tar_load}()},
\code{\link{tar_meta}()},
\code{\link{tar_objects}()},
\code{\link{tar_pid}()},
\code{\link{tar_process}()},
\code{\link{tar_read_raw}()},
\code{\link{tar_read}()}
}
\concept{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_target.R
\name{target_run_worker}
\alias{target_run_worker}
\title{Internal function to run a target on a worker.}
\usage{
target_run_worker(target, envir, path_store, fun, options, envvars)
}
\arguments{
\item{target}{A target object.}

\item{envir}{An environment or the string \code{"globalenv"}.}

\item{path_store}{Character of length 1, path to the data store.}

\item{fun}{Character of length 1, name of the user-side function called
to run the pipeline.}

\item{options}{List, exported from an object of class \code{"tar_options"}.}

\item{envvars}{Data frame of \code{targets}-specific environment variables
from \code{\link[=tar_envvars]{tar_envvars()}}.}
}
\description{
For internal purposes only. Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_watch.R
\name{tar_watch}
\alias{tar_watch}
\title{Shiny app to watch the dependency graph.}
\usage{
tar_watch(
  seconds = 10,
  seconds_min = 1,
  seconds_max = 60,
  seconds_step = 1,
  targets_only = FALSE,
  exclude = ".Random.seed",
  outdated = FALSE,
  label = NULL,
  level_separation = 150,
  degree_from = 1L,
  degree_to = 1L,
  config = Sys.getenv("TAR_CONFIG", "_targets.yaml"),
  project = Sys.getenv("TAR_PROJECT", "main"),
  height = "650px",
  display = "summary",
  displays = c("summary", "branches", "progress", "graph", "about"),
  background = TRUE,
  browse = TRUE,
  host = getOption("shiny.host", "127.0.0.1"),
  port = getOption("shiny.port", targets::tar_random_port()),
  verbose = TRUE,
  supervise = TRUE,
  poll_connection = TRUE,
  stdout = "|",
  stderr = "|"
)
}
\arguments{
\item{seconds}{Numeric of length 1,
default number of seconds between refreshes of the graph.
Can be changed in the app controls.}

\item{seconds_min}{Numeric of length 1, lower bound of \code{seconds}
in the app controls.}

\item{seconds_max}{Numeric of length 1, upper bound of \code{seconds}
in the app controls.}

\item{seconds_step}{Numeric of length 1, step size of \code{seconds}
in the app controls.}

\item{targets_only}{Logical, whether to restrict the output to just targets
(\code{FALSE}) or to also include global functions and objects.}

\item{exclude}{Character vector of nodes to omit from the graph.}

\item{outdated}{Logical, whether to show colors to distinguish outdated
targets from up-to-date targets. (Global functions and objects
still show these colors.) Looking for outdated targets
takes a lot of time for large pipelines with lots of branches,
and setting \code{outdated} to \code{FALSE} is a nice way to speed up the graph
if you only want to see dependency relationships and build progress.}

\item{label}{Label argument to \code{\link[=tar_visnetwork]{tar_visnetwork()}}.}

\item{level_separation}{Numeric of length 1,
\code{levelSeparation} argument of \code{visNetwork::visHierarchicalLayout()}.
Controls the distance between hierarchical levels.
Consider changing the value if the aspect ratio of the graph
is far from 1. If \code{level_separation} is \code{NULL},
the \code{levelSeparation} argument of \code{visHierarchicalLayout()}
defaults to \code{150}.}

\item{degree_from}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_from}
controls the number of edges the neighborhood extends upstream.}

\item{degree_to}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_to}
controls the number of edges the neighborhood extends downstream.}

\item{config}{Character of length 1, file path of the YAML
configuration file with \code{targets} project settings.
The \code{config} argument specifies which YAML configuration
file that \code{tar_config_get()} reads from or \code{tar_config_set()}
writes to in a single function call.
It does not globally change which configuration file is used
in subsequent function calls. The default file path of the YAML
file is always \verb{_targets.yaml} unless you set another
default path using the \code{TAR_CONFIG} environment variable,
e.g. \code{Sys.setenv(TAR_CONFIG = "custom.yaml")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}

\item{project}{Character of length 1, name of the current
\code{targets} project. Thanks to the \code{config} R package,
\code{targets} YAML configuration files can store multiple
sets of configuration settings, with each set corresponding
to its own project. The \code{project} argument allows you to
set or get a configuration setting for a specific project
for a given call to \code{tar_config_set()} or \code{tar_config_get()}.
The default project is always called \code{"main"}
unless you set another
default project using the \code{TAR_PROJECT} environment variable,
e.g. \code{Sys.setenv(tar_project = "custom")}. This also has the
effect of temporarily modifying the default arguments to other functions
such as \code{\link[=tar_make]{tar_make()}} because the default arguments
to those functions are controlled by \code{tar_config_get()}.}

\item{height}{Character of length 1,
height of the \code{visNetwork} widget and branches table.}

\item{display}{Character of length 1, which display to show first.}

\item{displays}{Character vector of choices for the display.
Elements can be any of
\code{"graph"}, \code{"summary"}, \code{"branches"}, or \code{"about"}.}

\item{background}{Logical, whether to run the app in a background process
so you can still use the R console while the app is running.}

\item{browse}{Whether to open the app in a browser when the app is ready.
Only relevant if \code{background} is \code{TRUE}.}

\item{host}{Character of length 1, IPv4 address to listen on.
Only relevant if \code{background} is \code{TRUE}.}

\item{port}{Positive integer of length 1, TCP port to listen on.
Only relevant if \code{background} is \code{TRUE}.}

\item{verbose}{whether to print a spinner and informative messages.
Only relevant if \code{background} is \code{TRUE}.}

\item{supervise}{Whether to register the process with a supervisor. If \code{TRUE},
the supervisor will ensure that the process is killed when the R process
exits.}

\item{poll_connection}{Whether to have a control connection to
the process. This is used to transmit messages from the subprocess
to the main process.}

\item{stdout}{The name of the file the standard output of
the child R process will be written to.
If the child process runs with the \code{--slave} option (the default),
then the commands are not echoed and will not be shown
in the standard output. Also note that you need to call \code{print()}
explicitly to show the output of the command(s).}

\item{stderr}{The name of the file the standard error of
the child R process will be written to.
In particular \code{message()} sends output to the standard
error. If nothing was sent to the standard error, then this file
will be empty. This argument can be the same file as \code{stdout},
in which case they will be correctly interleaved. If this is the
string \code{"2>&1"}, then standard error is redirected to standard output.}
}
\value{
A handle to \code{callr::r_bg()} background process running the app.
}
\description{
Launches a background process with a Shiny app
that calls \code{\link[=tar_visnetwork]{tar_visnetwork()}} every few seconds.
To embed this app in other apps, use the Shiny module
in \code{\link[=tar_watch_ui]{tar_watch_ui()}} and \code{\link[=tar_watch_server]{tar_watch_server()}}.
}
\details{
The controls of the app are in the left panel.
The \code{seconds} control is the number of seconds between
refreshes of the graph, and the other settings match
the arguments of \code{\link[=tar_visnetwork]{tar_visnetwork()}}.
}
\examples{
if (identical(Sys.getenv("TAR_INTERACTIVE_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  sleep_run <- function(...) {
    Sys.sleep(10)
  }
  list(
    tar_target(settings, sleep_run()),
    tar_target(data1, sleep_run(settings)),
    tar_target(data2, sleep_run(settings))
  )
}, ask = FALSE)
# Launch the app in a background process.
tar_watch(seconds = 10, outdated = FALSE, targets_only = TRUE)
# Run the pipeline.
tar_make()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_branch_names.R
\name{tar_branch_names}
\alias{tar_branch_names}
\title{Branch names}
\usage{
tar_branch_names(name, index, store = targets::tar_config_get("store"))
}
\arguments{
\item{name}{Symbol, name of the dynamic branching target (pattern).}

\item{index}{Integer vector of branch indexes.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A character vector of branch names.
}
\description{
Get the branch names of a dynamic branching target
using numeric indexes.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(4)),
    tar_target(y, 2 * x, pattern = map(x)),
    tar_target(z, y, pattern = map(y))
  )
}, ask = FALSE)
tar_make()
tar_branch_names(z, c(2, 3))
})
}
}
\seealso{
Other branching: 
\code{\link{tar_branch_index}()},
\code{\link{tar_branch_names_raw}()},
\code{\link{tar_branches}()},
\code{\link{tar_pattern}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_started.R
\name{tar_started}
\alias{tar_started}
\title{List started targets.}
\usage{
tar_started(names = NULL, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Optional, names of the targets. If supplied, the
function restricts its output to these targets.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A character vector of started targets.
}
\description{
List targets whose progress is \code{"started"}.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
tar_started()
tar_started(starts_with("y_")) # see also all_of()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_cue.R
\name{tar_cue}
\alias{tar_cue}
\title{Declare the rules that cue a target.}
\usage{
tar_cue(
  mode = c("thorough", "always", "never"),
  command = TRUE,
  depend = TRUE,
  format = TRUE,
  iteration = TRUE,
  file = TRUE
)
}
\arguments{
\item{mode}{Cue mode. If \code{"thorough"}, all the cues apply unless
individually suppressed. If \code{"always"}, then the target always
runs. If \code{"never"}, then the target does not run unless the
metadata does not exist or the last run errored.}

\item{command}{Logical, whether to rerun the target if command changed
since last time.}

\item{depend}{Logical, whether to rerun the target if the value of one
of the dependencies changed.}

\item{format}{Logical, whether to rerun the target if the user-specified
storage format changed. The storage format is user-specified through
\code{\link[=tar_target]{tar_target()}} or \code{\link[=tar_option_set]{tar_option_set()}}.}

\item{iteration}{Logical, whether to rerun the target if the user-specified
iteration method changed. The iteration method is user-specified through
\code{\link[=tar_target]{tar_target()}} or \code{\link[=tar_option_set]{tar_option_set()}}.}

\item{file}{Logical, whether to rerun the target if the file(s) with the
return value changed or at least one is missing.}
}
\description{
Declare the rules that mark a target as outdated.
}
\section{Target invalidation rules}{

\code{targets} uses internal metadata and special cues
to decide whether a target is up to date (can skip)
or is outdated/invalidated (needs to rerun). By default,
\code{targets} moves through the following list of cues
and declares a target outdated if at least one is cue activated.
\enumerate{
\item There is no metadata record of the target.
\item The target errored last run.
\item The target has a different class than it did before.
\item The cue mode equals \code{"always"}.
\item The cue mode does not equal \code{"never"}.
\item The \code{command} metadata field (the hash of the R command)
is different from last time.
\item The \code{depend} metadata field (the hash of the immediate upstream
dependency targets and global objects) is different from last time.
\item The storage format is different from last time.
\item The iteration mode is different from last time.
\item A target's file (either the one in \verb{_targets/objects/}
or a dynamic file) does not exist or changed since last time.
}

The user can suppress many of the above cues using the \code{tar_cue()}
function, which creates the \code{cue} argument of \code{\link[=tar_target]{tar_target()}}.
Cues objects also constitute more nuanced target invalidation rules.
The \code{tarchetypes} package has many such examples, including
\code{tar_age()}, \code{tar_download()}, \code{tar_cue_age()}, \code{tar_cue_force()},
and \code{tar_cue_skip()}.
}

\section{Dependency-based invalidation and user-defined functions}{

If the cue of a target has \code{depend = TRUE} (default) then the target
is marked invalidated/outdated when its upstream dependencies change.
A target's dependencies include upstream targets,
user-defined functions, and other global objects populated
in the target script file (default: \verb{_targets.R}).
To determine if a given dependency changed
since the last run of the pipeline, \code{targets} computes hashes.
The hash of a target is computed on its files in storage
(usually a file in \verb{_targets/objects/}). The hash of a
non-function global object dependency is computed directly on its
in-memory data. User-defined functions are hashed in the following way:
\enumerate{
\item Deparse the function with \code{targets:::tar_deparse_safe()}. This
function computes a string representation of the function
body and arguments. This string representation is invariant to
changes in comments and whitespace, which means
trivial changes to formatting do not cue targets to rerun.
\item Manually remove any literal pointers from the function string
using \code{targets:::mask_pointers()}. Such pointers arise from
inline compiled C/C++ functions.
\item Using static code analysis (i.e. \code{\link[=tar_deps]{tar_deps()}}, which is based on
\code{codetools::findGlobals()}) identify any user-defined functions
and global objects that the current function depends on.
Append the hashes of those dependencies to the string representation
of the current function.
\item Compute the hash of the final string representation using
\code{targets:::digest_chr64()}.
}

Above, (3) is important because user-defined functions
have dependencies of their own, such as other user-defined
functions and other global objects. (3) ensures that a change to
a function's dependencies invalidates the function itself, which
in turn invalidates any calling functions and any targets downstream
with the \code{depend} cue turned on.
}

\examples{
# The following target will always run when the pipeline runs.
x <- tar_target(x, download_data(), cue = tar_cue(mode = "always"))
}
\seealso{
Other targets: 
\code{\link{tar_format}()},
\code{\link{tar_target_raw}()},
\code{\link{tar_target}()}
}
\concept{targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources_feather.R
\name{tar_resources_feather}
\alias{tar_resources_feather}
\title{Target resources: feather storage formats}
\usage{
tar_resources_feather(compression = "default", compression_level = NULL)
}
\arguments{
\item{compression}{Character of length 1, \code{compression}
argument of \code{arrow::write_feather()}.}

\item{compression_level}{Numeric of length 1, \code{compression_level}
argument of \code{arrow::write_feather()}.}
}
\value{
Object of class \code{"tar_resources_feather"}, to be supplied
to the feather argument of \code{tar_resources()}.
}
\description{
Create the feather argument of \code{tar_resources()}
to specify optional settings for feather data frame storage formats
powered by the \code{arrow} R package.
See the \code{format} argument of \code{\link[=tar_target]{tar_target()}} for details.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
# Somewhere in you target script file (usually _targets.R):
tar_target(
  name,
  command(),
  format = "feather",
  resources = tar_resources(
    feather = tar_resources_feather(compression = "lz4")
  )
)
}
\seealso{
Other resources: 
\code{\link{tar_resources_aws}()},
\code{\link{tar_resources_clustermq}()},
\code{\link{tar_resources_fst}()},
\code{\link{tar_resources_future}()},
\code{\link{tar_resources_gcp}()},
\code{\link{tar_resources_parquet}()},
\code{\link{tar_resources_qs}()},
\code{\link{tar_resources_url}()},
\code{\link{tar_resources}()}
}
\concept{resources}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_seed.R
\name{tar_seed}
\alias{tar_seed}
\title{Get the random number generator seed of the target currently running.}
\usage{
tar_seed(default = 1L)
}
\arguments{
\item{default}{Integer, value to return if \code{tar_seed()}
is called on its own outside a \code{targets} pipeline.
Having a default lets users run things without \code{\link[=tar_make]{tar_make()}},
which helps peel back layers of code and troubleshoot bugs.}
}
\value{
Integer of length 1. If invoked inside a \code{targets} pipeline,
the return value is the seed of the target currently running,
which is a deterministic function of the target name. Otherwise,
the return value is \code{default}.
}
\description{
Get the random number generator seed
of the target currently running.
}
\details{
A target's random number generator seed
is a deterministic function of its name. In this way,
each target runs with a reproducible seed so someone else
running the same pipeline should get the same results,
and no two targets in the same pipeline share the same seed.
(Even dynamic branches have different names and thus different seeds.)
You can retrieve the seed of a completed target
with \code{tar_meta(your_target, seed)}
and run \code{set.seed()} on the result to locally
recreate the target's initial RNG state.
}
\examples{
tar_seed()
tar_seed(default = 123L)
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script(tar_target(returns_seed, tar_seed()), ask = FALSE)
tar_make()
tar_read(returns_seed)
})
}
}
\seealso{
Other utilities: 
\code{\link{tar_active}()},
\code{\link{tar_call}()},
\code{\link{tar_cancel}()},
\code{\link{tar_definition}()},
\code{\link{tar_envir}()},
\code{\link{tar_group}()},
\code{\link{tar_name}()},
\code{\link{tar_path}()},
\code{\link{tar_store}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_exist_progress.R
\name{tar_exist_progress}
\alias{tar_exist_progress}
\title{Check if progress metadata exists.}
\usage{
tar_exist_progress(store = targets::tar_config_get("store"))
}
\arguments{
\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
Logical of length 1, whether the current project's metadata exists.
}
\description{
Check if the progress metadata file \verb{_targets/meta/progress}
exists for the current project.
}
\details{
To learn more about local storage in \code{targets}, visit
\url{https://books.ropensci.org/targets/files.html#internal-files}.
}
\examples{
tar_exist_progress()
}
\seealso{
Other existence: 
\code{\link{tar_exist_meta}()},
\code{\link{tar_exist_objects}()},
\code{\link{tar_exist_process}()},
\code{\link{tar_exist_script}()}
}
\concept{existence}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio_addin_tar_load.R
\name{rstudio_addin_tar_load}
\alias{rstudio_addin_tar_load}
\title{RStudio addin to call \code{\link[=tar_load]{tar_load()}} on the symbol at the cursor.}
\usage{
rstudio_addin_tar_load(context = NULL)
}
\arguments{
\item{context}{RStudio API context from
\code{rstudioapi::getActiveDocumentContext()}.}
}
\description{
For internal use only. Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_sitrep.R
\name{tar_sitrep}
\alias{tar_sitrep}
\title{Show the cue-by-cue status of each target.}
\usage{
tar_sitrep(
  names = NULL,
  fields = NULL,
  shortcut = targets::tar_config_get("shortcut"),
  reporter = targets::tar_config_get("reporter_outdated"),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function, reporter),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Optional, names of the targets. If supplied, \code{tar_sitrep()}
only returns metadata on these targets.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=starts_with]{starts_with()}}.}

\item{fields}{Optional, names of columns/fields to select. If supplied,
\code{tar_sitrep()} only returns the selected metadata columns.
You can supply symbols or \code{tidyselect} helpers
like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
The \code{name} column is always included first
no matter what you select. Choices:
\itemize{
\item \code{name}: name of the target or global object.
\item \code{record}: Whether the \code{record} cue is activated:
\code{TRUE} if the target is not in the metadata (\code{\link[=tar_meta]{tar_meta()}}),
or if the target errored during the last \code{\link[=tar_make]{tar_make()}},
or if the class of the target changed.
\item \code{always}: Whether \code{mode} in \code{\link[=tar_cue]{tar_cue()}} is \code{"always"}.
If \code{TRUE}, \code{\link[=tar_make]{tar_make()}} always runs the target.
\item \code{never}: Whether \code{mode} in \code{\link[=tar_cue]{tar_cue()}} is \code{"never"}.
If \code{TRUE}, \code{\link[=tar_make]{tar_make()}} will only run if the
\code{record} cue activates.
\item \code{command}: Whether the target's command changed since last time.
Always \code{TRUE} if the \code{record} cue is activated.
Otherwise, always \code{FALSE} if the \code{command} cue is suppressed.
\item \code{depend}: Whether the data/output of at least one of the target's
dependencies changed since last time.
Dependencies are targets, functions,
and global objects directly upstream.
Call \code{tar_outdated(targets_only = FALSE)} or
\code{tar_visnetwork(targets_only = FALSE)} to see exactly which
dependencies are outdated.
Always \code{NA} if the \code{record} cue is activated.
Otherwise, always \code{FALSE} if the \code{depend} cue is suppressed.
\item \code{format}: Whether the storage format of the target
is different from last time.
Always \code{NA} if the \code{record} cue is activated.
Otherwise, always \code{FALSE} if the \code{format} cue is suppressed.
\item \code{iteration}: Whether the iteration mode of the target
is different from last time.
Always \code{NA} if the \code{record} cue is activated.
Otherwise, always \code{FALSE} if the \code{iteration} cue is suppressed.
\item \code{file}: Whether the file(s) with the target's return value
are missing or different from last time.
Always \code{NA} if the \code{record} cue is activated.
Otherwise, always \code{FALSE} if the \code{file} cue is suppressed.
}}

\item{shortcut}{Logical of length 1, how to interpret the \code{names} argument.
If \code{shortcut} is \code{FALSE} (default) then the function checks
all targets upstream of \code{names} as far back as the dependency graph goes.
If \code{TRUE}, then the function only checks the targets in \code{names}
and uses stored metadata for information about upstream dependencies
as needed. \code{shortcut = TRUE} increases speed if there are a lot of
up-to-date targets, but it assumes all the dependencies
are up to date, so please use with caution.
Use with caution. \code{shortcut = TRUE} only works if you set \code{names}.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets are checked. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"forecast"}: print running totals of the checked and outdated
targets found so far.
}}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A data frame with one row per target/object and one column
per cue. Each element is a logical to indicate whether the cue
is activated for the target.
See the \code{field} argument in this help file for details.
}
\description{
For each target, report which cues are activated.
Except for the \code{never} cue, the target will rerun in \code{\link[=tar_make]{tar_make()}}
if any cue is activated. The target is suppressed if the \code{never}
cue is \code{TRUE}. See \code{\link[=tar_cue]{tar_cue()}} for details.
}
\details{
Caveats:
\itemize{
\item \code{\link[=tar_cue]{tar_cue()}} allows you to change/suppress cues, so the return
value will depend on the settings you supply to \code{\link[=tar_cue]{tar_cue()}}.
\item If a pattern tries to branches over a target that does not exist
in storage, then the branches are omitted from the output.
\item \code{tar_sitrep()} is myopic. It only considers what happens to the
immediate target and its immediate upstream dependencies,
and it makes no attempt to propagate invalidation downstream.
}
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
tar_sitrep()
tar_meta(starts_with("y_")) # see also all_of()
})
}
}
\seealso{
Other inspect: 
\code{\link{tar_deps_raw}()},
\code{\link{tar_deps}()},
\code{\link{tar_glimpse}()},
\code{\link{tar_manifest}()},
\code{\link{tar_network}()},
\code{\link{tar_outdated}()},
\code{\link{tar_validate}()},
\code{\link{tar_visnetwork}()}
}
\concept{inspect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_skipped.R
\name{tar_skipped}
\alias{tar_skipped}
\title{List skipped targets.}
\usage{
tar_skipped(names = NULL, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Optional, names of the targets. If supplied, the
function restricts its output to these targets.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A character vector of skipped targets.
}
\description{
List targets whose progress is \code{"skipped"}.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
tar_skipped()
tar_skipped(starts_with("y_")) # see also all_of()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_visnetwork.R
\name{tar_visnetwork}
\alias{tar_visnetwork}
\title{Visualize an abridged fast dependency graph.}
\usage{
tar_visnetwork(
  targets_only = FALSE,
  names = NULL,
  shortcut = FALSE,
  allow = NULL,
  exclude = ".Random.seed",
  outdated = TRUE,
  label = NULL,
  level_separation = NULL,
  degree_from = 1L,
  degree_to = 1L,
  zoom_speed = 1,
  reporter = targets::tar_config_get("reporter_outdated"),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{targets_only}{Logical, whether to restrict the output to just targets
(\code{FALSE}) or to also include global functions and objects.}

\item{names}{Names of targets. The graph visualization will operate
only on these targets (and unless \code{shortcut} is \code{TRUE},
all the targets upstream as well). Selecting a small subgraph
using \code{names} could speed up the load time of the visualization.
Unlike \code{allow}, \code{names} is invoked before the graph
is generated.
Set to NULL to check/build all the targets (default).
Otherwise, you can supply symbols or tidyselect helpers
like starts_with().
Applies to ordinary targets (stem) and whole dynamic branching
targets (patterns) but not individual dynamic branches.}

\item{shortcut}{Logical of length 1, how to interpret the \code{names} argument.
If \code{shortcut} is \code{FALSE} (default) then the function checks
all targets upstream of \code{names} as far back as the dependency graph goes.
If \code{TRUE}, then the function only checks the targets in \code{names}
and uses stored metadata for information about upstream dependencies
as needed. \code{shortcut = TRUE} increases speed if there are a lot of
up-to-date targets, but it assumes all the dependencies
are up to date, so please use with caution.
Also, \code{shortcut = TRUE} only works if you set \code{names}.}

\item{allow}{Optional, define the set of allowable vertices in the graph.
Unlike \code{names}, \code{allow} is invoked only after the graph is mostly
resolved, so it will not speed up execution.
Set to \code{NULL} to allow all vertices in the pipeline and environment
(default). Otherwise, you can supply symbols or
\code{tidyselect} helpers like \code{\link[=starts_with]{starts_with()}}.}

\item{exclude}{Optional, define the set of exclude vertices from the graph.
Unlike \code{names}, \code{exclude} is invoked only after the graph is mostly
resolved, so it will not speed up execution.
Set to \code{NULL} to exclude no vertices.
Otherwise, you can supply symbols or \code{tidyselect}
helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{outdated}{Logical, whether to show colors to distinguish outdated
targets from up-to-date targets. (Global functions and objects
still show these colors.) Looking for outdated targets
takes a lot of time for large pipelines with lots of branches,
and setting \code{outdated} to \code{FALSE} is a nice way to speed up the graph
if you only want to see dependency relationships and build progress.}

\item{label}{Character vector of one or more aesthetics to add to the
vertex labels. Can contain \code{"time"} to show total runtime, \code{"size"}
to show total storage size, or \code{"branches"} to show the number of
branches in each pattern. You can choose multiple aesthetics
at once, e.g. \code{label = c("time", "branches")}. All are disabled
by default because they clutter the graph.}

\item{level_separation}{Numeric of length 1,
\code{levelSeparation} argument of \code{visNetwork::visHierarchicalLayout()}.
Controls the distance between hierarchical levels.
Consider changing the value if the aspect ratio of the graph
is far from 1. If \code{level_separation} is \code{NULL},
the \code{levelSeparation} argument of \code{visHierarchicalLayout()}
defaults to \code{150}.}

\item{degree_from}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_from}
controls the number of edges the neighborhood extends upstream.}

\item{degree_to}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_to}
controls the number of edges the neighborhood extends downstream.}

\item{zoom_speed}{Positive numeric of length 1, scaling factor on the
zoom speed. Above 1 zooms faster than default, below 1 zooms
lower than default.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets are checked. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"forecast"}: print running totals of the checked and outdated
targets found so far.
}}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A \code{visNetwork} HTML widget object.
}
\description{
Analyze the pipeline defined in the target script file
(default: \verb{_targets.R})
and visualize the directed acyclic graph of targets
and global functions and objects.
}
\examples{
if (identical(Sys.getenv("TAR_INTERACTIVE_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set()
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
})
tar_visnetwork()
tar_visnetwork(allow = starts_with("y")) # see also all_of()
})
}
}
\seealso{
Other inspect: 
\code{\link{tar_deps_raw}()},
\code{\link{tar_deps}()},
\code{\link{tar_glimpse}()},
\code{\link{tar_manifest}()},
\code{\link{tar_network}()},
\code{\link{tar_outdated}()},
\code{\link{tar_sitrep}()},
\code{\link{tar_validate}()}
}
\concept{inspect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_deps.R
\name{tar_deps}
\alias{tar_deps}
\title{Code dependencies}
\usage{
tar_deps(expr)
}
\arguments{
\item{expr}{A quoted R expression or function.}
}
\value{
Character vector of the dependencies of a function or expression.
}
\description{
List the dependencies of a function or expression.
}
\details{
\code{targets} detects the dependencies of commands using
static code analysis. Use \code{tar_deps()} to run the
code analysis and see the dependencies for yourself.
}
\examples{
tar_deps(x <- y + z)
tar_deps({
  x <- 1
  x + a
})
tar_deps(function(a = b) map_dfr(data, ~do_row(.x)))
}
\seealso{
Other inspect: 
\code{\link{tar_deps_raw}()},
\code{\link{tar_glimpse}()},
\code{\link{tar_manifest}()},
\code{\link{tar_network}()},
\code{\link{tar_outdated}()},
\code{\link{tar_sitrep}()},
\code{\link{tar_validate}()},
\code{\link{tar_visnetwork}()}
}
\concept{inspect}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_load_globals.R
\name{tar_load_globals}
\alias{tar_load_globals}
\title{Load globals for debugging, testing, and prototyping}
\usage{
tar_load_globals(
  envir = parent.frame(),
  script = targets::tar_config_get("script")
)
}
\arguments{
\item{envir}{Environment to source the target script (default: \verb{_targets.R}).
Defaults to the calling environment.}

\item{script}{Character of length 1, path to the target script file
that defines the pipeline (\verb{_targets.R} by default).
This path should be either
an absolute path or a path relative to the project root where you will
call \code{\link[=tar_make]{tar_make()}} and other functions. When \code{\link[=tar_make]{tar_make()}} and friends
run the script from the current working directory.
If the argument \code{NULL}, the setting is not modified.
Use \code{\link[=tar_config_unset]{tar_config_unset()}} to delete a setting.}
}
\value{
\code{NULL} (invisibly).
}
\description{
Load user-defined packages, functions, global objects, and
settings defined in the target script file (default: \verb{_targets.R}).
This function is for debugging, testing, and prototyping only.
It is not recommended for use inside a serious pipeline
or to report the results of a serious pipeline.
}
\details{
This function first sources the target script file
(default: \verb{_targets.R})
to loads all user-defined functions, global objects, and settings
into the current R process. Then, it loads all the packages defined
in \code{tar_option_get("packages")} (default: \code{(.packages())})
using \code{library()} with \code{lib.loc} defined in \code{tar_option_get("library")}
(default: \code{NULL}).
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set(packages = "callr")
  analyze_data <- function(data) {
    summary(data)
  }
  list(
    tar_target(x, 1 + 1),
    tar_target(y, 1 + 1)
  )
}, ask = FALSE)
tar_load_globals()
print(analyze_data)
print("callr" \%in\% (.packages()))
})
}
}
\seealso{
Other debug: 
\code{\link{tar_traceback}()},
\code{\link{tar_workspaces}()},
\code{\link{tar_workspace}()}
}
\concept{debug}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio_addin_tar_progress.R
\name{rstudio_addin_tar_progress}
\alias{rstudio_addin_tar_progress}
\title{RStudio addin to print \code{tail(tar_progress())}.}
\usage{
rstudio_addin_tar_progress()
}
\description{
For internal use only. Not a user-side function.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_resources_fst.R
\name{tar_resources_fst}
\alias{tar_resources_fst}
\title{Target resources: \code{fst} storage formats}
\usage{
tar_resources_fst(compress = 50)
}
\arguments{
\item{compress}{Numeric of length 1, \code{compress}
argument of \code{fst::write_fst()}.}
}
\value{
Object of class \code{"tar_resources_fst"}, to be supplied
to the \code{fst} argument of \code{tar_resources()}.
}
\description{
Create the \code{fst} argument of \code{tar_resources()}
to specify optional settings for big data frame storage formats
powered by the \code{fst} R package.
See the \code{format} argument of \code{\link[=tar_target]{tar_target()}} for details.
}
\section{Resources}{

Functions \code{\link[=tar_target]{tar_target()}} and \code{\link[=tar_option_set]{tar_option_set()}}
each takes an optional \code{resources} argument to supply
non-default settings of various optional backends for data storage
and high-performance computing. The \code{tar_resources()} function
is a helper to supply those settings in the correct manner.
Resources are all-or-nothing: if you specify any resources
with \code{\link[=tar_target]{tar_target()}}, all the resources from \code{tar_option_get("resources")}
are dropped for that target. In other words, if you write
\code{tar_option_set(resources = resources_1)} and then
\code{tar_target(x, my_command(), resources = resources_2)}, then everything
in \code{resources_1} is discarded for target \code{x}.
}

\examples{
# Somewhere in you target script file (usually _targets.R):
tar_target(
  name,
  command(),
  format = "fst_tbl",
  resources = tar_resources(
    fst = tar_resources_fst(compress = 100)
  )
)
}
\seealso{
Other resources: 
\code{\link{tar_resources_aws}()},
\code{\link{tar_resources_clustermq}()},
\code{\link{tar_resources_feather}()},
\code{\link{tar_resources_future}()},
\code{\link{tar_resources_gcp}()},
\code{\link{tar_resources_parquet}()},
\code{\link{tar_resources_qs}()},
\code{\link{tar_resources_url}()},
\code{\link{tar_resources}()}
}
\concept{resources}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_watch.R
\name{tar_watch_app_ui}
\alias{tar_watch_app_ui}
\title{Create the full \code{\link[=tar_watch]{tar_watch()}} app UI.}
\usage{
tar_watch_app_ui(
  seconds,
  seconds_min,
  seconds_max,
  seconds_step,
  targets_only,
  outdated,
  label,
  level_separation,
  degree_from,
  degree_to,
  height,
  display,
  displays
)
}
\arguments{
\item{seconds}{Numeric of length 1,
default number of seconds between refreshes of the graph.
Can be changed in the app controls.}

\item{seconds_min}{Numeric of length 1, lower bound of \code{seconds}
in the app controls.}

\item{seconds_max}{Numeric of length 1, upper bound of \code{seconds}
in the app controls.}

\item{seconds_step}{Numeric of length 1, step size of \code{seconds}
in the app controls.}

\item{targets_only}{Logical, whether to restrict the output to just targets
(\code{FALSE}) or to also include global functions and objects.}

\item{outdated}{Logical, whether to show colors to distinguish outdated
targets from up-to-date targets. (Global functions and objects
still show these colors.) Looking for outdated targets
takes a lot of time for large pipelines with lots of branches,
and setting \code{outdated} to \code{FALSE} is a nice way to speed up the graph
if you only want to see dependency relationships and build progress.}

\item{label}{Label argument to \code{\link[=tar_visnetwork]{tar_visnetwork()}}.}

\item{level_separation}{Numeric of length 1,
\code{levelSeparation} argument of \code{visNetwork::visHierarchicalLayout()}.
Controls the distance between hierarchical levels.
Consider changing the value if the aspect ratio of the graph
is far from 1. If \code{level_separation} is \code{NULL},
the \code{levelSeparation} argument of \code{visHierarchicalLayout()}
defaults to \code{150}.}

\item{degree_from}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_from}
controls the number of edges the neighborhood extends upstream.}

\item{degree_to}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_to}
controls the number of edges the neighborhood extends downstream.}

\item{height}{Character of length 1,
height of the \code{visNetwork} widget and branches table.}

\item{display}{Character of length 1, which display to show first.}

\item{displays}{Character vector of choices for the display.
Elements can be any of
\code{"graph"}, \code{"summary"}, \code{"branches"}, or \code{"about"}.}
}
\value{
A Shiny UI.
}
\description{
Only exported for infrastructure purposes.
Not a user-side function. Users should instead
call \code{\link[=tar_watch]{tar_watch()}} directly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_make.R
\name{tar_make}
\alias{tar_make}
\title{Run a pipeline of targets.}
\usage{
tar_make(
  names = NULL,
  shortcut = targets::tar_config_get("shortcut"),
  reporter = targets::tar_config_get("reporter_make"),
  callr_function = callr::r,
  callr_arguments = targets::callr_args_default(callr_function, reporter),
  envir = parent.frame(),
  script = targets::tar_config_get("script"),
  store = targets::tar_config_get("store")
)
}
\arguments{
\item{names}{Names of the targets to build or check. Set to \code{NULL} to
check/build all the targets (default). Otherwise, you can supply
\code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.
Because \code{\link[=tar_make]{tar_make()}} and friends run the pipeline in a new R session,
if you pass a character vector to a tidyselect helper, you will need
to evaluate that character vector early with \verb{!!}, e.g.
\code{tar_make(names = all_of(!!your_vector))}.
Applies to ordinary targets (stem) and whole dynamic branching targets
(patterns) but not to individual dynamic branches.}

\item{shortcut}{Logical of length 1, how to interpret the \code{names} argument.
If \code{shortcut} is \code{FALSE} (default) then the function checks
all targets upstream of \code{names} as far back as the dependency graph goes.
\code{shortcut = TRUE} increases speed if there are a lot of
up-to-date targets, but it assumes all the dependencies
are up to date, so please use with caution.
It relies on stored metadata for information about upstream dependencies.
\code{shortcut = TRUE} only works if you set \code{names}.}

\item{reporter}{Character of length 1, name of the reporter to user.
Controls how messages are printed as targets run in the pipeline.
Defaults to \code{tar_config_get("reporter_make")}. Choices:
\itemize{
\item \code{"silent"}: print nothing.
\item \code{"summary"}: print a running total of the number of each targets in
each status category (queued, started, skipped, build, canceled,
or errored). Also show a timestamp (\code{"\%H:\%M \%OS2"} \code{strptime()} format)
of the last time the progress changed and printed to the screen.
\item \code{"timestamp"}: same as the \code{"verbose"} reporter except that each
.message begins with a time stamp.
\item \code{"timestamp_positives"}: same as the \code{"timestamp"} reporter
except without messages for skipped targets.
\item \code{"verbose"}: print messages for individual targets
as they start, finish, or are skipped.
\item \code{"verbose_positives"}: same as the \code{"verbose"} reporter
except without messages for skipped targets.
}}

\item{callr_function}{A function from \code{callr} to start a fresh clean R
process to do the work. Set to \code{NULL} to run in the current session
instead of an external process (but restart your R session just before
you do in order to clear debris out of the global environment).
\code{callr_function} needs to be \code{NULL} for interactive debugging,
e.g. \code{tar_option_set(debug = "your_target")}.
However, \code{callr_function} should not be \code{NULL} for serious
reproducible work.}

\item{callr_arguments}{A list of arguments to \code{callr_function}.}

\item{envir}{An environment, where to run the target R script
(default: \verb{_targets.R}) if \code{callr_function} is \code{NULL}.
Ignored if \code{callr_function} is anything other than \code{NULL}.
\code{callr_function} should only be \code{NULL} for debugging and
testing purposes, not for serious runs of a pipeline, etc.

The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related
functions always overrides
the current value of \code{tar_option_get("envir")} in the current R session
just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.}

\item{script}{Character of length 1, path to the
target script file. Defaults to \code{tar_config_get("script")},
which in turn defaults to \verb{_targets.R}. When you set
this argument, the value of \code{tar_config_get("script")}
is temporarily changed for the current function call.
See \code{\link[=tar_script]{tar_script()}},
\code{\link[=tar_config_get]{tar_config_get()}}, and \code{\link[=tar_config_set]{tar_config_set()}} for details
about the target script file and how to set it
persistently for a project.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
\code{NULL} except if \code{callr_function = callr::r_bg()}, in which case
a handle to the \code{callr} background process is returned. Either way,
the value is invisibly returned.
}
\description{
Run the pipeline you defined in the targets
script file (default: \verb{_targets.R}). \code{tar_make()}
runs the correct targets in the correct order and stores the return
values in \verb{_targets/objects/}.
}
\examples{
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set()
  list(tar_target(x, 1 + 1))
})
tar_make()
tar_script({
  tar_option_set()
  list(
    tar_target(y1, 1 + 1),
    tar_target(y2, 1 + 1),
    tar_target(z, y1 + y2)
  )
}, ask = FALSE)
prefix <- "y"
tar_make(starts_with(!!prefix)) # Only builds y1 and y2.
})
}
\seealso{
Other pipeline: 
\code{\link{tar_make_clustermq}()},
\code{\link{tar_make_future}()}
}
\concept{pipeline}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_url.R
\name{tar_random_port}
\alias{tar_random_port}
\title{Random TCP port}
\usage{
tar_random_port(lower = 49152L, upper = 65355L)
}
\arguments{
\item{lower}{Integer of length 1, lowest possible port.}

\item{upper}{Integer of length 1, highest possible port.}
}
\value{
A random port not likely to be used by another process.
}
\description{
Not a user-side function. Exported for infrastructure
purposes only.
}
\examples{
tar_random_port()
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_watch_ui.R
\name{tar_watch_ui}
\alias{tar_watch_ui}
\title{Shiny module UI for tar_watch()}
\usage{
tar_watch_ui(
  id,
  label = "tar_watch_label",
  seconds = 10,
  seconds_min = 1,
  seconds_max = 60,
  seconds_step = 1,
  targets_only = FALSE,
  outdated = FALSE,
  label_tar_visnetwork = NULL,
  level_separation = 150,
  degree_from = 1L,
  degree_to = 1L,
  height = "650px",
  display = "summary",
  displays = c("summary", "branches", "progress", "graph", "about")
)
}
\arguments{
\item{id}{Character of length 1, ID corresponding to the UI function
of the module.}

\item{label}{Label for the module.}

\item{seconds}{Numeric of length 1,
default number of seconds between refreshes of the graph.
Can be changed in the app controls.}

\item{seconds_min}{Numeric of length 1, lower bound of \code{seconds}
in the app controls.}

\item{seconds_max}{Numeric of length 1, upper bound of \code{seconds}
in the app controls.}

\item{seconds_step}{Numeric of length 1, step size of \code{seconds}
in the app controls.}

\item{targets_only}{Logical, whether to restrict the output to just targets
(\code{FALSE}) or to also include global functions and objects.}

\item{outdated}{Logical, whether to show colors to distinguish outdated
targets from up-to-date targets. (Global functions and objects
still show these colors.) Looking for outdated targets
takes a lot of time for large pipelines with lots of branches,
and setting \code{outdated} to \code{FALSE} is a nice way to speed up the graph
if you only want to see dependency relationships and build progress.}

\item{label_tar_visnetwork}{Character vector, \code{label} argument to
\code{\link[=tar_visnetwork]{tar_visnetwork()}}.}

\item{level_separation}{Numeric of length 1,
\code{levelSeparation} argument of \code{visNetwork::visHierarchicalLayout()}.
Controls the distance between hierarchical levels.
Consider changing the value if the aspect ratio of the graph
is far from 1. If \code{level_separation} is \code{NULL},
the \code{levelSeparation} argument of \code{visHierarchicalLayout()}
defaults to \code{150}.}

\item{degree_from}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_from}
controls the number of edges the neighborhood extends upstream.}

\item{degree_to}{Integer of length 1. When you click on a node,
the graph highlights a neighborhood of that node. \code{degree_to}
controls the number of edges the neighborhood extends downstream.}

\item{height}{Character of length 1,
height of the \code{visNetwork} widget and branches table.}

\item{display}{Character of length 1, which display to show first.}

\item{displays}{Character vector of choices for the display.
Elements can be any of
\code{"graph"}, \code{"summary"}, \code{"branches"}, or \code{"about"}.}
}
\value{
A Shiny module UI.
}
\description{
Use \code{tar_watch_ui()} and \code{\link[=tar_watch_server]{tar_watch_server()}}
to include \code{\link[=tar_watch]{tar_watch()}} as a Shiny module in an app.
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_errored}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_option_set.R
\name{tar_option_set}
\alias{tar_option_set}
\title{Set target options.}
\usage{
tar_option_set(
  tidy_eval = NULL,
  packages = NULL,
  imports = NULL,
  library = NULL,
  envir = NULL,
  format = NULL,
  iteration = NULL,
  error = NULL,
  memory = NULL,
  garbage_collection = NULL,
  deployment = NULL,
  priority = NULL,
  backoff = NULL,
  resources = NULL,
  storage = NULL,
  retrieval = NULL,
  cue = NULL,
  debug = NULL,
  workspaces = NULL,
  workspace_on_error = NULL
)
}
\arguments{
\item{tidy_eval}{Logical, whether to enable tidy evaluation
when interpreting \code{command} and \code{pattern}. If \code{TRUE}, you can use the
"bang-bang" operator \verb{!!} to programmatically insert
the values of global objects.}

\item{packages}{Character vector of packages to load right before
the target builds. Use \code{tar_option_set()} to set packages
globally for all subsequent targets you define.}

\item{imports}{Character vector of package names to track
global dependencies. For example, if you write
\code{tar_option_set(imports = "yourAnalysisPackage")} early in your
target script file (default: \verb{_targets.R})
then \code{tar_make()} will automatically rerun or skip targets
in response to changes to the R functions and objects defined in
\code{yourAnalysisPackage}. Does not account for low-level compiled code
such as C/C++ or Fortran. If you supply multiple packages,
e.g. \code{tar_option_set(imports = c("p1", "p2"))}, then the objects in
\code{p1} override the objects in \code{p2} if there are name conflicts.
Similarly, objects in \code{tar_option_get("envir")} override
everything in \code{tar_option_get("imports")}.}

\item{library}{Character vector of library paths to try
when loading \code{packages}.}

\item{envir}{Environment containing functions and global objects
common to all targets in the pipeline.
The \code{envir} argument of \code{\link[=tar_make]{tar_make()}} and related functions
always overrides the current value of \code{tar_option_get("envir")}
in the current R session just before running the target script file,
so whenever you need to set an alternative \code{envir}, you should always set
it with \code{tar_option_set()} from within the target script file.
In other words, if you call \code{tar_option_set(envir = envir1)} in an
interactive session and then
\code{tar_make(envir = envir2, callr_function = NULL)},
then \code{envir2} will be used.

If \code{envir} is the global environment, all the promise objects
are diffused before sending the data to parallel workers
in \code{\link[=tar_make_future]{tar_make_future()}} and \code{\link[=tar_make_clustermq]{tar_make_clustermq()}},
but otherwise the environment is unmodified.
This behavior improves performance by decreasing
the size of data sent to workers.

If \code{envir} is not the global environment, then it should at least inherit
from the global environment or base environment
so \code{targets} can access attached packages.
In the case of a non-global \code{envir}, \code{targets} attempts to remove
potentially high memory objects that come directly from \code{targets}.
That includes \code{tar_target()} objects of class \code{"tar_target"},
as well as objects of class \code{"tar_pipeline"} or \code{"tar_algorithm"}.
This behavior improves performance by decreasing
the size of data sent to workers.

Package environments should not be assigned to \code{envir}.
To include package objects as upstream dependencies in the pipeline,
assign the package to the \code{packages} and \code{imports} arguments
of \code{tar_option_set()}.}

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
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[=tar_make_future]{tar_make_future()}}. If \code{"worker"},
the target builds on a parallel worker. If \code{"main"},
the target builds on the host machine / process managing the pipeline.}

\item{priority}{Numeric of length 1 between 0 and 1. Controls which
targets get deployed first when multiple competing targets are ready
simultaneously. Targets with priorities closer to 1 get built earlier
(and polled earlier in \code{\link[=tar_make_future]{tar_make_future()}}).}

\item{backoff}{Numeric of length 1, must be greater than or equal to 0.01.
Maximum upper bound of the random polling interval
for the priority queue (seconds).
In high-performance computing (e.g. \code{\link[=tar_make_clustermq]{tar_make_clustermq()}}
and \code{\link[=tar_make_future]{tar_make_future()}}) it can be expensive to repeatedly poll the
priority queue if no targets are ready to process. The number of seconds
between polls is \code{runif(1, 0.001, max(backoff, 0.001 * 1.5 ^ index))},
where \code{index} is the number of consecutive polls so far that found
no targets ready to skip or run.
(If no target is ready, \code{index} goes up by 1. If a target is ready,
\code{index} resets to 0. For more information on exponential,
backoff, visit \url{https://en.wikipedia.org/wiki/Exponential_backoff}).
Raising \code{backoff} is kinder to the CPU etc. but may incur delays
in some instances.}

\item{resources}{Object returned by \code{tar_resources()}
with optional settings for high-performance computing
functionality, alternative data storage formats,
and other optional capabilities of \code{targets}.
See \code{tar_resources()} for details.}

\item{storage}{Character of length 1, only relevant to
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[=tar_make_future]{tar_make_future()}}.
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
\code{\link[=tar_path]{tar_path()}} from inside the target.
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
\code{\link[=tar_make_clustermq]{tar_make_clustermq()}} and \code{\link[=tar_make_future]{tar_make_future()}}.
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

\item{debug}{Character vector of names of targets to run in debug mode.
To use effectively, you must set \code{callr_function = NULL} and
restart your R session just before running. You should also
\code{\link[=tar_make]{tar_make()}}, \code{\link[=tar_make_clustermq]{tar_make_clustermq()}}, or \code{\link[=tar_make_future]{tar_make_future()}}.
For any target mentioned in \code{debug}, \code{targets} will force the target to
build locally (with \code{tar_cue(mode = "always")} and \code{deployment = "main"}
in the settings) and pause in an interactive debugger to help you diagnose
problems. This is like inserting a \code{browser()} statement at the
beginning of the target's expression, but without invalidating any
targets.}

\item{workspaces}{Character vector of target names.
Could be non-branching targets, whole dynamic branching targets,
or individual branch names. \code{\link[=tar_make]{tar_make()}} and friends
will save workspace files for these targets even if
the targets are skipped. Workspace files help with debugging.
See \code{\link[=tar_workspace]{tar_workspace()}} for details about workspaces.}

\item{workspace_on_error}{Logical of length 1, whether to save
a workspace file for each target that throws an error.
Workspace files help with debugging.
See \code{\link[=tar_workspace]{tar_workspace()}} for details about workspaces.}
}
\value{
\code{NULL} (invisibly).
}
\description{
Set target options, including default arguments to
\code{\link[=tar_target]{tar_target()}} such as packages, storage format,
iteration type, and cue. Only the non-null arguments are actually
set as options. See currently set options with \code{\link[=tar_option_get]{tar_option_get()}}.
To use \code{tar_option_set()} effectively, put it in your workflow's
target script file (default: \verb{_targets.R})
before calls to \code{\link[=tar_target]{tar_target()}} or \code{\link[=tar_target_raw]{tar_target_raw()}}.
}
\examples{
tar_option_get("format") # default format before we set anything
tar_target(x, 1)$settings$format
tar_option_set(format = "fst_tbl") # new default format
tar_option_get("format")
tar_target(x, 1)$settings$format
tar_option_reset() # reset the format
tar_target(x, 1)$settings$format
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  tar_option_set(cue = tar_cue(mode = "always")) # All targets always run.
  list(tar_target(x, 1), tar_target(y, 2))
})
tar_make()
tar_make()
})
}
}
\seealso{
Other configuration: 
\code{\link{tar_config_get}()},
\code{\link{tar_config_set}()},
\code{\link{tar_config_unset}()},
\code{\link{tar_envvars}()},
\code{\link{tar_option_get}()},
\code{\link{tar_option_reset}()}
}
\concept{configuration}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_language.R
\name{tar_language}
\alias{tar_language}
\alias{tar_deparse_language}
\alias{tar_deparse_safe}
\alias{tar_tidy_eval}
\alias{tar_tidyselect_eval}
\title{Language}
\usage{
tar_deparse_language(expr)

tar_deparse_safe(expr, collapse = "\\n", backtick = TRUE)

tar_tidy_eval(expr, envir, tidy_eval)

tar_tidyselect_eval(names_quosure, choices)
}
\arguments{
\item{expr}{A language object to modify or deparse.}

\item{collapse}{Character of length 1, delimiter in deparsing.}

\item{backtick}{logical indicating whether symbolic names should be
    enclosed in backticks if they do not follow the standard syntax.}

\item{envir}{An environment to find objects for tidy evaluation.}

\item{tidy_eval}{Logical of length 1, whether to apply tidy evaluation.}

\item{names_quosure}{An \code{rlang} quosure with \code{tidyselect} expressions.}

\item{choices}{A character vector of choices for character elements
returned by tidy evaluation.}
}
\description{
These functions help with metaprogramming in
packages built on top of \code{targets}.
}
\details{
\itemize{
\item \code{tar_deparse_language()} is a wrapper around \code{tar_deparse_safe()}
which leaves character vectors and \code{NULL} objects alone,
which helps with subsequent user input validation.
\item \code{tar_deparse_safe()} is a wrapper around \code{base::deparse()}
with a custom set of fast default settings and guardrails
to ensure the output always has length 1.
\item \code{tar_tidy_eval()} applies tidy evaluation to a language object
and returns another language object.
\item \code{tar_tidyselect_eval()} applies \code{tidyselect} selection with
some special guardrails around \code{NULL} inputs.
}
}
\examples{
tar_deparse_language(quote(run_model()))
}
\seealso{
Other utilities to extend targets: 
\code{\link{tar_assert}},
\code{\link{tar_condition}},
\code{\link{tar_dir}()},
\code{\link{tar_test}()}
}
\concept{utilities to extend targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_branch_index.R
\name{tar_branch_index}
\alias{tar_branch_index}
\title{Integer branch indexes}
\usage{
tar_branch_index(names, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Character vector of branch names}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A named integer vector of branch indexes.
}
\description{
Get the integer indexes of individual branch names
within their corresponding dynamic branching targets.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(4)),
    tar_target(y, 2 * x, pattern = map(x)),
    tar_target(z, y, pattern = map(y))
  )
}, ask = FALSE)
tar_make()
names <- c(
  tar_meta(y, children)$children[[1]][c(2, 3)],
  tar_meta(z, children)$children[[1]][2]
)
names
tar_branch_index(names) # c(2, 3, 2)
})
}
}
\seealso{
Other branching: 
\code{\link{tar_branch_names_raw}()},
\code{\link{tar_branch_names}()},
\code{\link{tar_branches}()},
\code{\link{tar_pattern}()}
}
\concept{branching}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_dir.R
\name{tar_dir}
\alias{tar_dir}
\title{Execute code in a temporary directory.}
\usage{
tar_dir(code)
}
\arguments{
\item{code}{User-defined code.}
}
\value{
Return value of the user-defined code.
}
\description{
Runs code inside a new \code{tempfile()} directory
in order to avoid writing to the user's file space.
Used in examples and tests in order to comply with CRAN policies.
}
\examples{
tar_dir(file.create("only_exists_in_tar_dir"))
file.exists("only_exists_in_tar_dir")
}
\seealso{
Other utilities to extend targets: 
\code{\link{tar_assert}},
\code{\link{tar_condition}},
\code{\link{tar_language}},
\code{\link{tar_test}()}
}
\concept{utilities to extend targets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tar_errored.R
\name{tar_errored}
\alias{tar_errored}
\title{List errored targets.}
\usage{
tar_errored(names = NULL, store = targets::tar_config_get("store"))
}
\arguments{
\item{names}{Optional, names of the targets. If supplied, the
function restricts its output to these targets.
You can supply symbols
or \code{tidyselect} helpers like \code{\link[=all_of]{all_of()}} and \code{\link[=starts_with]{starts_with()}}.}

\item{store}{Character of length 1, path to the
\code{targets} data store. Defaults to \code{tar_config_get("store")},
which in turn defaults to \verb{_targets/}.
When you set this argument, the value of \code{tar_config_get("store")}
is temporarily changed for the current function call.
See \code{\link[=tar_config_get]{tar_config_get()}} and \code{\link[=tar_config_set]{tar_config_set()}} for details
about how to set the data store path persistently
for a project.}
}
\value{
A character vector of errored targets.
}
\description{
List targets whose progress is \code{"errored"}.
}
\examples{
if (identical(Sys.getenv("TAR_EXAMPLES"), "true")) {
tar_dir({ # tar_dir() runs code from a temporary directory.
tar_script({
  list(
    tar_target(x, seq_len(2)),
    tar_target(y, 2 * x, pattern = map(x))
  )
}, ask = FALSE)
tar_make()
tar_errored()
tar_errored(starts_with("y_")) # see also all_of()
})
}
}
\seealso{
Other progress: 
\code{\link{tar_built}()},
\code{\link{tar_canceled}()},
\code{\link{tar_poll}()},
\code{\link{tar_progress_branches}()},
\code{\link{tar_progress_summary}()},
\code{\link{tar_progress}()},
\code{\link{tar_skipped}()},
\code{\link{tar_started}()},
\code{\link{tar_watch_server}()},
\code{\link{tar_watch_ui}()},
\code{\link{tar_watch}()}
}
\concept{progress}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstudio_addin_tar_visnetwork.R
\name{rstudio_addin_tar_visnetwork}
\alias{rstudio_addin_tar_visnetwork}
\title{RStudio addin to call \code{\link[=tar_visnetwork]{tar_visnetwork()}}.}
\usage{
rstudio_addin_tar_visnetwork()
}
\description{
For internal use only. Not a user-side function.
}
\keyword{internal}
