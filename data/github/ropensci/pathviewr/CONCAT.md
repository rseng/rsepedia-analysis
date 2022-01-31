
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pathviewr <a href='https://docs.ropensci.org/pathviewr'><img src='https://github.com/ropensci/pathviewr/raw/master/images/pathviewrhex_300dpi_trns.png' align="right" height="150px" /></a>

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build
status](https://github.com/ropensci/pathviewr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/pathviewr/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/pathviewr/graph/badge.svg)](https://codecov.io/gh/ropensci/pathviewr?branch=master)
[![](https://badges.ropensci.org/409_status.svg)](https://github.com/ropensci/software-review/issues/409)  
[![DOI](https://zenodo.org/badge/268906628.svg)](https://zenodo.org/badge/latestdoi/268906628)
[![CRAN
status](https://www.r-pkg.org/badges/version/pathviewr)](https://CRAN.R-project.org/package=pathviewr)
<!-- badges: end -->

`pathviewr` offers tools to import, clean, and visualize movement data,
particularly from motion capture systems such as [Optitrack’s
Motive](https://optitrack.com/software/motive/), the [Straw Lab’s
Flydra](https://github.com/strawlab/flydra), or other sources. We
provide functions to remove artifacts, standardize tunnel position and
tunnel axes, select a region of interest, isolate specific trajectories,
fill gaps in trajectory data, and calculate 3D and per-axis velocity.
For experiments of visual guidance, we also provide functions that use
subject position to estimate perception of visual stimuli.

## Installation

You can install `pathviewr` from CRAN via:

``` r
install.packages("pathviewr")
```

Or to get the latest (developmental) version through GitHub, use:

``` r
devtools::install_github("ropensci/pathviewr")
```

## Example

#### Data import and cleaning via `pathviewr`

We’ll also load two `tidyverse` packages for wrangling & plotting in
this readme.

``` r
library(pathviewr)
library(ggplot2)
library(magrittr)
```

We will import and clean a sample data set from `.csv` files exported by
Optitrack’s [Motive](https://optitrack.com/software/motive/) software.
For examples of how to import and clean other types of data, [see the
Basics of data import and cleaning
vignette](https://docs.ropensci.org/pathviewr/articles/data-import-cleaning.html).

``` r
## Import the Motive example data included in 
## the package

motive_data <-
  read_motive_csv(
    system.file("extdata", "pathviewr_motive_example_data.csv",
                package = 'pathviewr')
  )
```

Several functions to clean and wrangle data are available, and we have a
suggested pipeline for how these steps should be handled. For this
example, we will use one of two “all-in-one” functions: `clean_viewr()`.
[See the Basics of data import and cleaning
vignette](https://docs.ropensci.org/pathviewr/articles/data-import-cleaning.html)
for the full pipeline and the other “all-in-one” function.

``` r
motive_allinone <-
  motive_data %>%
  clean_viewr(
    relabel_viewr_axes = TRUE,
    gather_tunnel_data = TRUE,
    trim_tunnel_outliers = TRUE,
    standardization_option = "rotate_tunnel",
    select_x_percent = TRUE,
    desired_percent = 50,
    rename_viewr_characters = FALSE,
    separate_trajectories = TRUE,
    max_frame_gap = "autodetect",
    get_full_trajectories = TRUE,
    span = 0.95
  )
#> autodetect is an experimental feature -- please report issues.

## Quick plot
## Colors correspond to unique trajectories (file_sub_traj)
motive_allinone %>%
  ggplot(aes(x = position_length, y = position_width, 
             fill = file_sub_traj)) +
  geom_point(pch = 21) +
  coord_fixed() +
  theme_classic() +
  theme(
    legend.position = "none"
  )
```

<img src="man/figures/README-all_in_one-1.png" width="100%" />

To get a sense of what we’ve done, compare the data before and after it
has passed through the pipeline.

``` r
## Check out the data's structure before cleaning and wrangling:
str(motive_data)
#> tibble[,26] [934 x 26] (S3: tbl_df/tbl/data.frame)
#>  $ frame                     : int [1:934] 72210 72211 72212 72213 72214 72215 72216 72217 72218 72219 ...
#>  $ time_sec                  : num [1:934] 722 722 722 722 722 ...
#>  $ device02_rotation_x       : num [1:934] 0.1346 0.0819 0.2106 0.1961 0.1305 ...
#>  $ device02_rotation_y       : num [1:934] -0.977 -0.978 -0.973 -0.972 -0.975 ...
#>  $ device02_rotation_z       : num [1:934] -0.1117 -0.0991 -0.0939 -0.1275 -0.1213 ...
#>  $ device02_rotation_w       : num [1:934] 0.1215 0.1654 0.0311 0.0351 0.1315 ...
#>  $ device02_position_x       : num [1:934] 0.142 0.137 0.125 0.118 0.113 ...
#>  $ device02_position_y       : num [1:934] 0.16 0.164 0.166 0.168 0.173 ...
#>  $ device02_position_z       : num [1:934] 2 1.97 1.95 1.92 1.89 ...
#>  $ device02_mean_marker_error: num [1:934] 0.000113 0.000105 0.000115 0.000202 0.000106 0.000095 0.000114 0.000117 0.000121 0.000131 ...
#>  $ device03_rotation_x       : num [1:934] 0.107 0.111 0.109 0.109 0.108 ...
#>  $ device03_rotation_y       : num [1:934] -0.295 -0.295 -0.295 -0.295 -0.295 ...
#>  $ device03_rotation_z       : num [1:934] -0.088 -0.0866 -0.0853 -0.0853 -0.0879 ...
#>  $ device03_rotation_w       : num [1:934] 0.945 0.945 0.945 0.945 0.945 ...
#>  $ device03_position_x       : num [1:934] 0.222 0.222 0.222 0.222 0.222 ...
#>  $ device03_position_y       : num [1:934] 0.245 0.245 0.245 0.245 0.245 ...
#>  $ device03_position_z       : num [1:934] 0.0597 0.0597 0.0598 0.0598 0.0598 ...
#>  $ device03_mean_marker_error: num [1:934] 0.000166 0.000172 0.000164 0.000163 0.000162 0.000162 0.000169 0.00017 0.00017 0.000213 ...
#>  $ device05_rotation_x       : num [1:934] 0.00672 0.00714 0.00709 0.00742 0.00826 ...
#>  $ device05_rotation_y       : num [1:934] 0.944 0.944 0.944 0.944 0.944 ...
#>  $ device05_rotation_z       : num [1:934] -0.117 -0.116 -0.118 -0.118 -0.117 ...
#>  $ device05_rotation_w       : num [1:934] 0.308 0.308 0.309 0.31 0.308 ...
#>  $ device05_position_x       : num [1:934] 0.173 0.173 0.173 0.173 0.173 ...
#>  $ device05_position_y       : num [1:934] 0.243 0.243 0.243 0.243 0.243 ...
#>  $ device05_position_z       : num [1:934] 2.66 2.66 2.66 2.66 2.66 ...
#>  $ device05_mean_marker_error: num [1:934] 0.000241 0.000247 0.000255 0.000244 0.00023 0.000226 0.000231 0.000236 0.000242 0.000263 ...
#>  - attr(*, ".internal.selfref")=<externalptr> 
#>  - attr(*, "pathviewr_steps")= chr "viewr"
#>  - attr(*, "file_id")= chr "pathviewr_motive_example_data.csv"
#>  - attr(*, "file_mtime")= POSIXct[1:1], format: "2021-05-13 14:51:31"
#>  - attr(*, "frame_rate")= num 100
#>  - attr(*, "header")='data.frame':   11 obs. of  2 variables:
#>   ..$ metadata: chr [1:11] "Format Version" "Take Name" "Take Notes" "Capture Frame Rate" ...
#>   ..$ value   : chr [1:11] "1.23" "sept-18_mixed-group_16-30" "" "100.000000" ...
#>  - attr(*, "Motive_IDs")= chr [1:24] "\"9E207518D8A311E969D7AB6B1FACE49B\"" "\"9E207518D8A311E969D7AB6B1FACE49B\"" "\"9E207518D8A311E969D7AB6B1FACE49B\"" "\"9E207518D8A311E969D7AB6B1FACE49B\"" ...
#>  - attr(*, "subject_names_full")= chr [1:24] "device02" "device02" "device02" "device02" ...
#>  - attr(*, "subject_names_simple")= chr [1:3] "device02" "device03" "device05"
#>  - attr(*, "data_names")= chr [1:26] "frame" "time_sec" "device02_rotation_x" "device02_rotation_y" ...
#>  - attr(*, "data_types_full")= chr [1:24] "Rigid Body" "Rigid Body" "Rigid Body" "Rigid Body" ...
#>  - attr(*, "data_types_simple")= chr "Rigid Body"
#>  - attr(*, "d1")= chr [1:26] "" "" "Rotation" "Rotation" ...
#>  - attr(*, "d2")= chr [1:26] "Frame" "Time (Seconds)" "X" "Y" ...
#>  - attr(*, "import_method")= chr "motive"

## Check out the data's structure after cleaning and wrangling:
str(motive_allinone)
#> tibble[,24] [449 x 24] (S3: tbl_df/tbl/data.frame)
#>  $ frame            : int [1:449] 72213 72214 72215 72216 72217 72218 72219 72220 72221 72222 ...
#>  $ time_sec         : num [1:449] 722 722 722 722 722 ...
#>  $ subject          : chr [1:449] "device02" "device02" "device02" "device02" ...
#>  $ position_length  : num [1:449] 0.647 0.62 0.593 0.567 0.541 ...
#>  $ position_width   : num [1:449] -0.112 -0.116 -0.122 -0.134 -0.141 ...
#>  $ position_height  : num [1:449] -0.0371 -0.0324 -0.0273 -0.0235 -0.0209 ...
#>  $ rotation_length  : num [1:449] -0.128 -0.121 -0.105 -0.106 -0.149 ...
#>  $ rotation_width   : num [1:449] 0.1961 0.1305 0.0935 0.1798 0.164 ...
#>  $ rotation_height  : num [1:449] -0.972 -0.975 -0.975 -0.975 -0.972 ...
#>  $ rotation_real    : num [1:449] 0.0351 0.1315 0.1734 0.0807 0.0824 ...
#>  $ mean_marker_error: num [1:449] 0.000202 0.000106 0.000095 0.000114 0.000117 0.000121 0.000131 0.00014 0.000113 0.000114 ...
#>  $ velocity         : num [1:449] 2.73 2.78 2.84 2.85 2.68 ...
#>  $ length_inst_vel  : num [1:449] -2.65 -2.72 -2.74 -2.58 -2.56 ...
#>  $ width_inst_vel   : num [1:449] -0.642 -0.387 -0.58 -1.139 -0.75 ...
#>  $ height_inst_vel  : num [1:449] 0.184 0.475 0.508 0.379 0.258 ...
#>  $ traj_id          : int [1:449] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ file_sub_traj    : chr [1:449] "pathviewr_motive_example_data.csv_device02_0" "pathviewr_motive_example_data.csv_device02_0" "pathviewr_motive_example_data.csv_device02_0" "pathviewr_motive_example_data.csv_device02_0" ...
#>  $ traj_length      : int [1:449] 63 63 63 63 63 63 63 63 63 63 ...
#>  $ start_length     : num [1:449] 0.647 0.647 0.647 0.647 0.647 ...
#>  $ end_length       : num [1:449] -0.656 -0.656 -0.656 -0.656 -0.656 ...
#>  $ length_diff      : num [1:449] 1.3 1.3 1.3 1.3 1.3 ...
#>  $ start_length_sign: num [1:449] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ end_length_sign  : num [1:449] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#>  $ direction        : chr [1:449] "leftwards" "leftwards" "leftwards" "leftwards" ...
#>  - attr(*, "file_id")= chr "pathviewr_motive_example_data.csv"
#>  - attr(*, "file_mtime")= POSIXct[1:1], format: "2021-05-13 14:51:31"
#>  - attr(*, "frame_rate")= num 100
#>  - attr(*, "header")='data.frame':   11 obs. of  2 variables:
#>   ..$ metadata: chr [1:11] "Format Version" "Take Name" "Take Notes" "Capture Frame Rate" ...
#>   ..$ value   : chr [1:11] "1.23" "sept-18_mixed-group_16-30" "" "100.000000" ...
#>  - attr(*, "Motive_IDs")= chr [1:24] "\"9E207518D8A311E969D7AB6B1FACE49B\"" "\"9E207518D8A311E969D7AB6B1FACE49B\"" "\"9E207518D8A311E969D7AB6B1FACE49B\"" "\"9E207518D8A311E969D7AB6B1FACE49B\"" ...
#>  - attr(*, "subject_names_full")= chr [1:24] "device02" "device02" "device02" "device02" ...
#>  - attr(*, "subject_names_simple")= chr [1:3] "device02" "device03" "device05"
#>  - attr(*, "data_names")= chr [1:26] "frame" "time_sec" "device02_rotation_x" "device02_rotation_y" ...
#>  - attr(*, "data_types_full")= chr [1:24] "Rigid Body" "Rigid Body" "Rigid Body" "Rigid Body" ...
#>  - attr(*, "data_types_simple")= chr "Rigid Body"
#>  - attr(*, "d1")= chr [1:26] "" "" "Rotation" "Rotation" ...
#>  - attr(*, "d2")= chr [1:26] "Frame" "Time (Seconds)" "X" "Y" ...
#>  - attr(*, "import_method")= chr "motive"
#>  - attr(*, "pathviewr_steps")= chr [1:10] "viewr" "renamed_tunnel" "gathered_tunnel" "artifacts_removed" ...
#>  - attr(*, "perch1_midpoint_original")= num [1:3] 0 0.2 0.205
#>  - attr(*, "perch2_midpoint_original")= num [1:3] 2.54 0.24 0.205
#>  - attr(*, "tunnel_centerpoint_original")= num [1:3] 1.27 0.22 0.205
#>  - attr(*, "rotation_degrees")= num 0.902
#>  - attr(*, "rotation_radians")= num 0.0157
#>  - attr(*, "perch1_midpoint_current")= num [1:3] -1.27 4.65e-15 2.05e-01
#>  - attr(*, "perch2_midpoint_current")= num [1:3] 1.27 -4.65e-15 2.05e-01
#>  - attr(*, "percent_selected")= num 50
#>  - attr(*, "full_tunnel_length")= num 2.64
#>  - attr(*, "selected_tunnel_length")= num 1.32
#>  - attr(*, "max_frame_gap")= int [1:3] 1 1 2
#>  - attr(*, "span")= num 0.95
#>  - attr(*, "trajectories_removed")= int 5
```

An important aspect of how `pathviewr` defines trajectories is by
managing gaps in the data. [See the vignette on Managing frame
gaps](https://docs.ropensci.org/pathviewr/articles/managing-frame-gaps.html)
for more information on trajectory definition and frame gaps.

Now that the data is cleaned, `pathviewr` includes functions that
estimate visual perceptions based on the distance between the
subject/observer and visual stimuli on the walls of the experimental
tunnel. For a complete description of these functions, [see the vignette
on Estimating visual perceptions from tracking
data](https://docs.ropensci.org/pathviewr/articles/visual-perception-functions.html).

#### Add more info about experiments

Now that our objects have been cleaned, we will use
`insert_treatments()` to add information about the experiments that are
necessary for calculating visual perceptions.

The data from this example were recorded in a V-shaped tunnel.
Accordingly, the vertex angle and vertex height of the tunnel, along
with information about the visual stimuli used during the experiment,
will be added to the data to inform calculations of visual perception
(next section).

``` r
motive_V <- 
  motive_allinone %>%
  insert_treatments(
    tunnel_config = "v",
    perch_2_vertex = 0.4,
    vertex_angle = 90,
    tunnel_length = 2,
    stim_param_lat_pos = 0.1,
    stim_param_lat_neg = 0.1,
    stim_param_end_pos = 0.3,
    stim_param_end_neg = 0.3,
    treatment = "lat10_end_30"
  ) 
```

#### Estimate perception of visual stimuli

To calculate the spatial frequency of the visual stimuli as perceived by
the subject some distance from the stimuli, we will use `get_sf()`.

This will require two intermediate steps: 1) calculating the minimum
distance between a subject and each wall (via `calc_min_dist_v()`) and
2) estimating the visual angles from the subject’s perspective
(`get_vis_angle()`).

``` r
motive_V_sf <- 
  motive_V %>%
  calc_min_dist_v(simplify_output = TRUE) %>%
  get_vis_angle() %>%
  get_sf()
```

Visualizing the calculations provides an more intuitive understanding of
how these visual perceptions change as the subject moves throughout the
tunnel. Please [see the vignette on Estimating visual perceptions from
tracking
data](https://docs.ropensci.org/pathviewr/articles/visual-perception-functions.html)
for more examples of visualizing calculations.

``` r
ggplot(motive_V_sf, aes(x = position_width, y = position_height)) +
  geom_point(aes(color = sf_pos), shape=1, size=3) +
  geom_segment(aes(x = 0,         # dimensions of the positive wall
                  y = -0.3855,
                  xend = 0.5869,
                  yend = 0.2014)) +
  geom_segment(aes(x = 0,         # dimensions of the negative wall
                   y = -0.3855,
                   xend = -0.5869,
                   yend = 0.2014)) +
  coord_fixed() +
  theme_classic() +
  theme(
    legend.position = "none"
  )
```

<img src="man/figures/README-motive_V_sf_pos-1.png" width="100%" />

## Contributing and/or raising Issues

We welcome feedback on bugs, improvements, and/or feature requests.
Please [see our Issues templates on
GitHub](https://github.com/ropensci/pathviewr/issues/new/choose) to make
a bug fix request or feature request.

To contribute code via a pull request, please consult our [Contributing
Guide](https://github.com/ropensci/pathviewr/blob/master/.github/CONTRIBUTING.md)
first.

## Citation

The preferred way to cite `pathviewr` (but subject to change) is:

Baliga VB, Armstrong MS, Press ER (2021). *pathviewr: Tools to import,
clean, and visualize animal movement data in R*. R package version
1.1.0, <https://github.com/ropensci/pathviewr>. doi:
10.5281/zenodo.4270187.

## License

GPL (&gt;= 3) + file LICENSE

🐢
# pathviewr 1.1.0
* New data cleaning functions added: set_traj_frametime(),
get_traj_velocities(), clean_by_span(), remove_duplicate_frames(), and
remove_vel_anomalies()
* These new functions have not been thoroughly vetted nor have unit tests
been written for them -- please use with caution and report issues.

# pathviewr 1.0.0
* Package has been accepted by rOpenSci and is now hosted on ropensci/pathviewr
* No changes to code since v0.9.5

# pathviewr 0.9.5
* Package has been updated to incorporate feedback from rOpenSci reviewers
@asbonnetlebrun and @marcosci, along with editor @maelle. See here for more
details: ropensci/software-review#409

# pathviewr 0.9.4
* We are targeting a submission to rOpenSci in the near future (hopefully
today). This version of pathviewr has been prepped according to rOpenSci's
"Packages: Development, Maintenance, and Peer Review" guide.
## Test environments
* local R installation, R 4.0.2
* windows-latest (release) on GitHub Actions
* macOS-latest (release) on GitHub Actions
* ubuntu-20.04 (release) on GitHub Actions
* ubuntu-20.04 (devel) on GitHub Actions

## R CMD check results

0 errors | 0 warnings | 0 notes

<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
# CONTRIBUTING #

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
*  Look at the GitHub Actions build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, an use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
---
name: 'Bug fix request'
about: 'Notify us of a bug so we can fix it'
title: "[bug_fix_request]"
labels: bug
assignees: ''

---

**Which function(s) is failing?**
A clear and concise description of what the problem is.

**Can you provide a reproducible example of the failure?**
Data on public repositories are preferred with steps to reproduce the failure. If that is not an option, please contact one of the `pathviewR` authors directly.

**What error message(s) are you receiving, if any?**

**Expected behavior or describe the solution you'd like**
A clear and concise description of what you expected to happen or what you would like to happen.

**What file type are you using?**
.csv from Motive, .mat from Flydra or other source using `as_viewr()` to import. If other, please describe.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: 'Feature Request: General'
about: 'Let us know if you have an idea to improve `pathviewR`'
title: "[feature_request_general]"
labels: enhancement
assignees: ''

---

**Please describe the nature of your feature request.**
Is the request to improve the tool(s) currently within `pathviewR` or to provide a tool that adds to the package's capabilities?

**Please explain the feature you would like to use in `pathviewR`.**
Briefly describe the goal of the feature and/or which current aspects of `pathviewR` relate to your request.

**Can you provide an example file to help explain your request?**
Data on public repositories are preferred. If that is not an option, please contact one of the `pathviewR` authors directly.

**Are you interested in contributing a function and/or code yourself?**
If so, please feel free to make a pull request for us to review.

**Additional context** 
Add any other context or screenshots about the feature request here.
---
name: 'Feature request: import function for file type'
about: 'Ask us about creating an import function for a specific file type '
title: "[filetype_request]"
labels: enhancement
assignees: ''

---

**Is your feature request related to a specific file type? Please describe.**
A clear and concise description of what the problem is. What is the file type and are the files created by a specific program or software suite?

**Can you provide an example file of the specified type?**
Data on public repositories are preferred. If that is not an option, please contact one of the `pathviewR` authors directly.

**Describe the solution you'd like**
A clear and concise description of what you want to happen when the data are imported from the file.

**Does `as_viewr()` not work?**
We are very happy to create new functions for additional file types, but we'd also like to know if our current workaround does not work for you. Please provide a description of whether `as_viewr()` can work or if it fails in this use case.

**Are you interested in contributing a function and/or code yourself?**
If so, please feel free to make a pull request for us to review.

**Additional context**
Add any other context or screenshots about the feature request here.
Response to reviewers for rOpenSci reviews - round 1
================

Hi @asbonnetlebrun, @marcosci, and @maelle,

Thank you sincerely for taking the time and effort to review our submission. 

Below, you will find point-by-point responses to each of your comments. Your
original comments may not appear fully (for sake of concision), but please note
that we believe that we have addressed each item to the best of our capability.

Since addressing these items was handled over several commits, the most 
appropriate commit that reflects the change will also be noted within each 
of our responses.

Please also note that the package is now entitled `pathviewr` instead of the 
original `pathviewR`, per your advice.

## Review from Reviewer 1 (asbonnetlebrun)

### Documentation
#### **Examples** for all exported functions in R Help that run successfully locally

> There are examples for almost all exported functions in R help. Exceptions are
the calc_sf_box(), calc_sf_V(), calc_vis_angle_box() and calc_vis_angle_V()
functions (and standardised_tunnel, although there is a reason given in the R
help – because no example data is provided with pathviewr to run the code on).

We added in examples for each of the following functions. Please note that we
consoldiated some of the visual guidance functions, so some of the original
functions referred to in this reviewer comment have now been replaced or 
renamed.  

- get_sf() [(99c9ef5)](https://github.com/vbaliga/pathviewr/commit/99c9ef517e265cd27ba790b1cb28929e1986e469)  
- get_vis_angle() [(dfd1e2f)](https://github.com/vbaliga/pathviewr/commit/d5d1e2f1726472505d19c091ed4363590fddc855)[(6a36758)](https://github.com/vbaliga/pathviewr/commit/6a367587a06a30780cc8568eb19935237aa9ccb8)  
- get_min_dist_box() [(6a36758)](https://github.com/vbaliga/pathviewr/commit/6a367587a06a30780cc8568eb19935237aa9ccb8)  
- get_min_dist_v() [(b5f1445)](https://github.com/vbaliga/pathviewr/commit/b5f1445b743695ed27acb607321d71b34a6a689a)  


#### **Community guidelines** including contribution guidelines in the README or CONTRIBUTING, and DESCRIPTION with URL, BugReports and Maintainer (which may be autogenerated via Authors@R).

> Only in the DESCRIPTION file, but not in the README. 

We have added contribution guidelines (with links to Issues pages) in our 
README. [(f5d47e7)](https://github.com/vbaliga/pathviewr/commit/f5d47e7c473a9dc08331dbe78d882cd6c360bcf8)  

#### **Functionality**
*We will address these items in the `Review Comments` section below*

> Should the author(s) deem it appropriate, I agree to be acknowledged as a
package reviewer ("rev" role) in the package DESCRIPTION file.

We have added you as a "rev" within our description file
[(5b0661)](https://github.com/vbaliga/pathviewr/commit/5b06610017ac3388e34695eae8fc8e6de69b457c).
Please let us know if you would like your name to appear differently and/or
would like to adjust contact info etc.

### Review Comments

#### Data import and cleaning 
 
> My understanding is that the “Relabeling axes” and “gathering data columns”
parts are only relevant in the case of Motive data (Flydra data already comes in
the correct format, and as_viewr(), the function to import other types of data
also takes data in the correct format and already relabels the columns). Maybe
this could be made a bit clearer in the vignette?

We have clarified language of this vignette to indicate that relabeling &
gathering are only necessary in certain cases (e.g. using Motive data). Please
let us know what you think.
[(591dd50)](https://github.com/vbaliga/pathviewr/commit/591dd50dcafb8ff1195ca983cc25170cb49a80dd)

> Maybe that was me not being very familiar with the kind of experiments the
package is relevant for, but I had trouble to understand when we would use the
standardize_tunnel() function or the rotate_tunnel() function. Maybe you could
provide some examples when these functions to be applied? Also I struggled to
understand the standardize_tunnel() function, maybe particularly because there
was no example in the help file. When would the landmarks already be in the
data? I guess it will never be the case with Flydra data, considering that the
read_flydra_mat() function allows to input only a single subject_name. Is that
why it doesn’t need to be rotated, or is that something arising from the Flydra
software?

We have clarified the circumstances under which standardization is
needed and what types of landmarks are appropriate [(591dd50)](https://github.com/vbaliga/pathviewr/commit/591dd50dcafb8ff1195ca983cc25170cb49a80dd)    
Additionally, in the Help file for `standardize_tunnel()`, we have clarified this function's use 
cases [(591dd50)](https://github.com/vbaliga/pathviewr/commit/591dd50dcafb8ff1195ca983cc25170cb49a80dd)    

> Also, considering how the select_x_percent() function works (by selecting a
certain percentage of the tunnel length on each site of (0,0,0) along the length
axis), shouldn’t it be more appropriate to say that the (0,0,0) must be at the
centre of the region of interest, rather than at the centre of the tunnel?

Thanks! We have revised language of this vignette on what (0,0,0) represents [(591dd50)](https://github.com/vbaliga/pathviewr/commit/591dd50dcafb8ff1195ca983cc25170cb49a80dd)  

> Minor point: the link to the vignette for managing frame gaps is missing in
the text.

Thanks for catching this -- we have added the link  [(591dd50)](https://github.com/vbaliga/pathviewr/commit/591dd50dcafb8ff1195ca983cc25170cb49a80dd)  

#### Managing frame gaps 

> Very minor point: could it be clearer in the visualize_frame_gap_choice() help
or in the vignette that the loops argument represents the maximum frame gap
considered (i.e. that each loop represents an increment of 1 on the frame gap
value)?

We have added details of what `loops` means to both the vignette [(bcf6424)](https://github.com/vbaliga/pathviewr/commit/bcf64241867e6da0c96181a4ed11379c4eb5b646) 
and to the Help file [(bcf6424)](https://github.com/vbaliga/pathviewr/commit/bcf64241867e6da0c96181a4ed11379c4eb5b646)  

> Could there be cases when frame gaps can vary between devices (i.e. if I
understood well in the case of the Motive data, between subjects)? If so, would
it be relevant to allow the autodetect approach to be applied to each subject
separately? (or maybe that is not relevant...)

Sorry for being unclear - this function actually has this feature! We have made
a slight modification to the language that will hopefully make it more clear [(62b63b7)](https://github.com/vbaliga/pathviewr/commit/62b63b7c771036b84540ac9a7a5d259cd66fa07b)  

#### Visual perception 

> This is where I struggled, but maybe that is because I don’t come from the
specific field of application of this package. I felt like some more details
could be included in the vignette. Maybe start with a few examples of
experiments where this package can be useful. Maybe you could say right at the
start that these functions are relevant for experiments with animals flying in
tunnels with visual stimuli consisting of sine wave gratings on the tunnel
walls?

Thanks. We have heavily revised the language of the first couple paragraphs 
accordingly, and a figure has also been inserted to hopefully give readers a 
better sense of the concepts at hand. Please let us know what you think!  [(97704bc)](https://github.com/vbaliga/pathviewr/commit/97704bcc46d42257aa113562f39297ec619b3a18)  

> Also, my understanding is that here you implement only two cases:

>    Motive example: birds flying through a V-shaped tunnel, with visual stimuli
on each side of the tunnel
>    Flydra example: birds flying through a rectangular tunnel, with visual
stimuli on the side and front walls. So maybe make it clear that this section
only applies to these two specific cases?

> On this, a question: is there any plan to code more situations or are these
the typical situations for these kinds of experiment? If you welcome suggestions
from users of other settings, maybe you could say that in the vignette in the
same way you mention that you are happy to work towards the inclusion of more
data types?

We have revised the language of how these experiments are introduced [(97704bc)](https://github.com/vbaliga/pathviewr/commit/97704bcc46d42257aa113562f39297ec619b3a18). 
We have also provided links to the Issues page where users can request e.g. 
different tunnel setups [(20bb54a)](https://github.com/vbaliga/pathviewr/commit/20bb54a8b8ec3d5bbe6c8a98033488225ede7bfc)  
 
> More generally, I would have appreciated some definition of what you call
“spatial frequency” and “visual angle”. Although I understand that these might
be obvious for people in the field (which are ultimately the target users of the
package), but I feel like these terms (visual angle, and in particular spatial
frequency) are of such broad usage that they might deserve some better
description of what they mean in the context of the package. Maybe include a
diagram in the vignette for users to visualise what these values represent? In
particular, the Value section in the R help for the calc\_sf\_box() function
seems to mention that the spatial frequency is the number of cycles per degree
of visual angle. Maybe this info could be included in the vignette (and in the
Description section of the function help)?

We have added in definitions along with citations of some journal articles
that provide nice summaries of these topics. [(fa7b12e)](https://github.com/vbaliga/pathviewr/commit/fa7b12e2e64477fd106eb2f454136bca8a4629c8)  


#### Comments on the code of the visual perception functions

> I seem to understand that the functions to calculate visual angle and spatial
frequency in the box example don’t handle the front wall (only the side walls),
am I right? Is there any reason for this? If the front wall is not relevant,
maybe there is no need to include the option to add parameters for it in the
insert_treatments() function?

Thanks for the suggestion. The visual perception functions now include 
calculations for the end walls, though the outputs from these functions only 
include the end wall towards which the subject is moving towards. Please see
the changes to following functions:  

- get_sf() [(99c9ef5)](https://github.com/vbaliga/pathviewr/commit/99c9ef517e265cd27ba790b1cb28929e1986e469)  
- get_vis_angle() [(dfd1e2f)](https://github.com/vbaliga/pathviewr/commit/d5d1e2f1726472505d19c091ed4363590fddc855)[(6a36758)](https://github.com/vbaliga/pathviewr/commit/6a367587a06a30780cc8568eb19935237aa9ccb8)   
- get_min_dist_box() [(6a36758)](https://github.com/vbaliga/pathviewr/commit/6a367587a06a30780cc8568eb19935237aa9ccb8)  
- get_min_dist_v() [(b5f1445)](https://github.com/vbaliga/pathviewr/commit/b5f1445b743695ed27acb607321d71b34a6a689a)  


*Note from VBB:* for brevity, the details of your calculations will not be 
included here.

> On that note, there is no need for an ifelse() to calculate these distances,
these lines could simply replace these lines and these lines in calc_sf_box()
and calc_vis_angle_box()

Thanks! We have now replaced the `ifelse()` lines in the corresponding 
functions, which also have been renamed:  

- get_sf() [(99c9ef5)](https://github.com/vbaliga/pathviewr/commit/99c9ef517e265cd27ba790b1cb28929e1986e469)  
- get_vis_angle()[(6a36758)](https://github.com/vbaliga/pathviewr/commit/6a367587a06a30780cc8568eb19935237aa9ccb8)   

> Can you please correct the calculations, and adapt the
test-calc_vis_angle_box.R file?

Thanks sincerely for catching this! We have made corrections to the calculations
of vis angle and SF  [(6a36758)](https://github.com/vbaliga/pathviewr/commit/6a367587a06a30780cc8568eb19935237aa9ccb8)[(99c9ef5)](https://github.com/vbaliga/pathviewr/commit/99c9ef517e265cd27ba790b1cb28929e1986e469).   

We have also updated `test-get_vis_angle.R` accordingly  Tests were written for all new functions (`calc_min_dist_v()`[(93808ea)](https://github.com/vbaliga/pathviewr/commit/93808eaa71e1cb8ccc333432ce18f0df839cbee5) , `calc_min_dist_box`[(93808ea)](https://github.com/vbaliga/pathviewr/commit/93808eaa71e1cb8ccc333432ce18f0df839cbee5) , `get_vis_angle()`[(ba8d813)](https://github.com/vbaliga/pathviewr/commit/ba8d813100101304c64e5dd400d28c6476b43af7)., and `get_sf()`)[(b4afca9)](https://github.com/vbaliga/pathviewr/commit/b4afca943d5c19a9e03614282f9e772ec5e49202) . 

> And also on this, I noticed that the user can supply negative values for the
neg_wall and pos_wall arguments in the insert_treatments() function, which in
this case would lead to spurious distances being calculated. Maybe insert a
warning/error message if that is the case, and add the appropriate tests?

Great catch! We have now added guidance to the Help file of 
`insert_treatments()` about feasible values for all arguments including new tunnel dimensions arguments `tunnel_width` and `tunnel_height`[(429258b)](https://github.com/vbaliga/pathviewr/commit/429258bf4039cf26671d915d507b48838c1aeb87). 
This includes a check within `insert_treatments()` that stops the function if 
faulty argument values are supplied  [(429258b)](https://github.com/vbaliga/pathviewr/commit/429258bf4039cf26671d915d507b48838c1aeb87). 
A test to `test-insert_treatments.R` to check the functionality of the
above check has also been added [(a9cc804)](https://github.com/vbaliga/pathviewr/commit/a9cc8041bb3d2dbe7cf7a883437b8f5fe8ccd692)  


#### Spatial frequency
> Another point, on the spatial frequency. If I understand it well (being the
number of cycles per degree of visual angle), I think there is an error in these
lines of calc_sf_V() and these lines of calc_sf_box().

Thanks, we have made corrections to the calculation of SF [(99c9ef5)](https://github.com/vbaliga/pathviewr/commit/99c9ef517e265cd27ba790b1cb28929e1986e469)  

> Just one final broader question, don’t these calculations assume that the bird
is parallel to the length axis of the tunnel? Is that not important/common
practice not to use the rotation information?

Yes, this is true and a great point. It is actually not common practice in the
field to use rotation information -- quite a few studies still rely solely on
positional data. That said, we definitely agree that adding rotation is an
important way to advance our understanding of visual guidance. At this time, we
allow for rotation data to be imported & wrangled along with all the positional
data, but rotation data have not yet been integrated into the visual perception
functions.  This is largely because how rotation data are encoded can be tricky
(depending on how the rotation matrix is composed and/or how orientation axes
are defined on each subject). Accordingly, we are still working on a way that
will work generically for most use cases.

For now, we have added a note in the Help files for each visual perception 
that rotation information may be integrated in future pathviewr updates  [(46ba850)](https://github.com/vbaliga/pathviewr/commit/46ba850c9d2b2db453365531280b8dde1ded831a). 
We have also done so in the Visual perception vignette [(a3a7795)](https://github.com/vbaliga/pathviewr/commit/a3a7795f21d75b8639d2c3dfab18ff789f299ea8)  

## Review from Reviewer 2 (marcosci)

> Should the author(s) deem it appropriate, I agree to be acknowledged as a
package reviewer ("rev" role) in the package DESCRIPTION file.

We have added you as a "rev" within our description file
[(5b0661)](https://github.com/vbaliga/pathviewr/commit/5b06610017ac3388e34695eae8fc8e6de69b457c).
Please let us know if you would like your name to appear differently and/or
would like to adjust contact info etc.

### Review Comments

> In the description etc. you state that pathviewr is a tool for "animal
movement data" ... Couldn't you make that broader? Part of the package is
actually capable to handle all kinds of 3D movement data and the visual
perception function should also be applicable to humans, or?

Yes thanks, this is a good point and we agree. We have made a slight alteration
to our description and readme. We don't want to oversell the features of our 
package, so we hope that where we have landed with the language strikes a nice
balance. [(deffd9d)](https://github.com/vbaliga/pathviewr/commit/deffd9da937c586a7a8685c7fa17bd8542d7ed21)

> It would actually be quite beneficial to give a short walkthrough of what
input data can look like. So, you need x,y,z ... but what more. And what defines
Optitrack and flydra data.

We agree and added a short walkthrough of what movement data look like, both 
generally and specifically in Motive and Flydra [(591dd50)](https://github.com/vbaliga/pathviewr/commit/591dd50dcafb8ff1195ca983cc25170cb49a80dd)  

> I did not see any contribution guidelines, so it would be helpful to include
those.

We have added contribution guidelines provided by rOpenSci (with slight 
modification) [(5252918)](https://github.com/vbaliga/pathviewr/commit/5252918ff124a315834dc7e060806c3ada3682dc)  

> The Maintainer field is missing in the DESCRIPTION - altough Vikram is listed
as CRE.

We've updated the Maintainer field with VBB as the maintainer [(9146e09)](https://github.com/vbaliga/pathviewr/commit/9146e099f699425e46ba31b6ff340ef911d11944cont)  

> pathviewR ... I submitted a package here once and the feedback was to either
stick with capital letters or just lower case. I don't know if that changed, but
it makes typing the package name out easier in my opinion. If that is a must
should probably be addressed by @maelle.

Yeah, we see what you mean and went ahead and changed the name to `pathviewr`. [(52e559c)](https://github.com/vbaliga/pathviewr/commit/52e559cb731907eb258f68b0a986aef0b6c784dc)


## Additional items from the editor (maelle)

> Regarding the package name, yes we recommend all lower-case, especially if the
package isn’t on CRAN yet. I myself renamed my Ropenaq package after review and
I don’t regret doing it.

As noted above, we changed the name to `pathviewr`. [(52e559c)](https://github.com/vbaliga/pathviewr/commit/52e559cb731907eb258f68b0a986aef0b6c784dc). 
In particular, we found [this guide](https://www.njtierney.com/post/2017/10/27/change-pkg-name/) 
to be very helpful -- we went step-by-step through everything in that page. We 
figured it would be good to let you know in case you are interested in providing
that sort of info in the rOpenSci Guide for Authors, though we understand it may
not be common enough of a problem to merit adding to the guide.

Thanks again for all your feedback and advice!

And sorry for the slight delay in getting back to you. My (VBB's) wife just
gave birth to our first child last week -- it's been a pretty wild ride!

Best regards,
Vikram
🐢
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# pathviewr <a href='https://docs.ropensci.org/pathviewr'><img src='https://github.com/ropensci/pathviewr/raw/master/images/pathviewrhex_300dpi_trns.png' align="right" height="150px" /></a>

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build status](https://github.com/ropensci/pathviewr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/pathviewr/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/pathviewr/graph/badge.svg)](https://codecov.io/gh/ropensci/pathviewr?branch=master)
[![](https://badges.ropensci.org/409_status.svg)](https://github.com/ropensci/software-review/issues/409)  
[![DOI](https://zenodo.org/badge/268906628.svg)](https://zenodo.org/badge/latestdoi/268906628)
[![CRAN status](https://www.r-pkg.org/badges/version/pathviewr)](https://CRAN.R-project.org/package=pathviewr)
<!-- badges: end -->

`pathviewr` offers tools to import, clean, and visualize movement data,
particularly from motion capture systems such as 
[Optitrack's Motive](https://optitrack.com/software/motive/), the 
[Straw Lab's Flydra](https://github.com/strawlab/flydra), or other sources. We
provide functions to remove artifacts, standardize tunnel position and tunnel
axes, select a region of interest, isolate specific trajectories, fill gaps in
trajectory data, and calculate 3D and per-axis velocity. For experiments of
visual guidance, we also provide functions that use subject position to estimate
perception of visual stimuli.

## Installation

You can install `pathviewr` from CRAN via:

``` {r install_cran, eval = FALSE}
install.packages("pathviewr")
```

Or to get the latest (developmental) version through GitHub, use:
  
``` {r install_github, eval = FALSE}
devtools::install_github("ropensci/pathviewr")
```

## Example

#### Data import and cleaning via `pathviewr`
We'll also load two `tidyverse` packages for wrangling & plotting in this 
readme.

```{r package_loading, message=FALSE, warning=FALSE}
library(pathviewr)
library(ggplot2)
library(magrittr)

```

We will import and clean a sample data set from `.csv` files exported by
Optitrack's [Motive](https://optitrack.com/software/motive/) software. For
examples of how to import and clean other types of data, 
[see the Basics of data import and cleaning vignette](https://docs.ropensci.org/pathviewr/articles/data-import-cleaning.html).

```{r import_motive}
## Import the Motive example data included in 
## the package

motive_data <-
  read_motive_csv(
    system.file("extdata", "pathviewr_motive_example_data.csv",
                package = 'pathviewr')
  )

``` 

Several functions to clean and wrangle data are available, and we have a
suggested pipeline for how these steps should be handled. For this example, we
will use one of two "all-in-one" functions: `clean_viewr()`. 
[See the Basics of data import and cleaning vignette](https://docs.ropensci.org/pathviewr/articles/data-import-cleaning.html)
for the full pipeline and the other "all-in-one" function.

```{r all_in_one, fig.height=3, fig.width=6, dpi=300}
motive_allinone <-
  motive_data %>%
  clean_viewr(
    relabel_viewr_axes = TRUE,
    gather_tunnel_data = TRUE,
    trim_tunnel_outliers = TRUE,
    standardization_option = "rotate_tunnel",
    select_x_percent = TRUE,
    desired_percent = 50,
    rename_viewr_characters = FALSE,
    separate_trajectories = TRUE,
    max_frame_gap = "autodetect",
    get_full_trajectories = TRUE,
    span = 0.95
  )

## Quick plot
## Colors correspond to unique trajectories (file_sub_traj)
motive_allinone %>%
  ggplot(aes(x = position_length, y = position_width, 
             fill = file_sub_traj)) +
  geom_point(pch = 21) +
  coord_fixed() +
  theme_classic() +
  theme(
    legend.position = "none"
  )
  
```

To get a sense of what we've done, compare the data before and after it has passed through the pipeline.

```{r compare_before_and_after}
## Check out the data's structure before cleaning and wrangling:
str(motive_data)

## Check out the data's structure after cleaning and wrangling:
str(motive_allinone)
```

An important aspect of how `pathviewr` defines trajectories is by managing gaps
in the data. 
[See the vignette on Managing frame gaps](https://docs.ropensci.org/pathviewr/articles/managing-frame-gaps.html)
for more information on trajectory definition and frame gaps.

Now that the data is cleaned, `pathviewr` includes functions that estimate
visual perceptions based on the distance between the subject/observer and visual
stimuli on the walls of the experimental tunnel. For a complete description of
these functions, 
[see the vignette on Estimating visual perceptions from tracking data](https://docs.ropensci.org/pathviewr/articles/visual-perception-functions.html).


#### Add more info about experiments
Now that our objects have been cleaned, we will use `insert_treatments()` to add
information about the experiments that are necessary for calculating visual
perceptions.

The data from this example were recorded in a V-shaped tunnel. Accordingly,
the vertex angle and vertex height of the tunnel, along with information about
the visual stimuli used during the experiment, will be added to the data
to inform calculations of visual perception (next section).

```{r insert_treats}
motive_V <- 
  motive_allinone %>%
  insert_treatments(
    tunnel_config = "v",
    perch_2_vertex = 0.4,
    vertex_angle = 90,
    tunnel_length = 2,
    stim_param_lat_pos = 0.1,
    stim_param_lat_neg = 0.1,
    stim_param_end_pos = 0.3,
    stim_param_end_neg = 0.3,
    treatment = "lat10_end_30"
  ) 
```


#### Estimate perception of visual stimuli
To calculate the spatial frequency of the visual stimuli as perceived by the
subject some distance from the stimuli, we will use `get_sf()`.

This will require two intermediate steps: 1) calculating the minimum distance
between a subject and each wall (via `calc_min_dist_v()`) and 2) estimating
the visual angles from the subject's perspective (`get_vis_angle()`).


```{r calc_sf_V}
motive_V_sf <- 
  motive_V %>%
  calc_min_dist_v(simplify_output = TRUE) %>%
  get_vis_angle() %>%
  get_sf()
```


Visualizing the calculations provides an more intuitive understanding of how
these visual perceptions change as the subject moves throughout the tunnel.
Please [see the vignette on Estimating visual perceptions from tracking data](https://docs.ropensci.org/pathviewr/articles/visual-perception-functions.html)  for more examples of visualizing calculations.

```{r motive_V_sf_pos, fig.height=3, fig.width=6, dpi=300}
ggplot(motive_V_sf, aes(x = position_width, y = position_height)) +
  geom_point(aes(color = sf_pos), shape=1, size=3) +
  geom_segment(aes(x = 0,         # dimensions of the positive wall
                  y = -0.3855,
                  xend = 0.5869,
                  yend = 0.2014)) +
  geom_segment(aes(x = 0,         # dimensions of the negative wall
                   y = -0.3855,
                   xend = -0.5869,
                   yend = 0.2014)) +
  coord_fixed() +
  theme_classic() +
  theme(
    legend.position = "none"
  )

```

## Contributing and/or raising Issues
We welcome feedback on bugs, improvements, and/or feature requests. Please
[see our Issues templates on GitHub](https://github.com/ropensci/pathviewr/issues/new/choose) to make a bug
fix request or feature request. 

To contribute code via a pull request, please consult our [Contributing Guide](https://github.com/ropensci/pathviewr/blob/master/.github/CONTRIBUTING.md) first.

## Citation

The preferred way to cite `pathviewr` (but subject to change) is:

Baliga VB, Armstrong MS, Press ER (2021). _pathviewr: Tools to import, clean, and visualize animal movement data in R_. R package version 1.1.0, [https://github.com/ropensci/pathviewr](https://github.com/ropensci/pathviewr). doi: 10.5281/zenodo.4270187.

## License

GPL (>= 3) + file LICENSE

🐢
---
title: "Managing frame gaps with pathviewr"
author: "Melissa S. Armstrong"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Managing frame gaps with pathviewr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteDepends{ggplot2}
  %\VignetteDepends{magrittr}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Defining trajectories with `pathviewr`
`pathviewr` defines trajectories as continuous movement from one side of a
region of interest to the other. Before trajectories can be defined, the region
of interest must be selected via `select_x_percent()` in addition to all the
previous steps of the data import and cleaning pipeline as described in 
[the Data Import and Cleaning vignette](https://docs.ropensci.org/pathviewr/articles/data-import-cleaning.html).
We'll start with loading `pathviewr` and a few packages from the `tidyverse` and
importing our data.

```{r package_loading, message=FALSE, warning=FALSE}
## If you do not already have pathviewr installed:
# install.packages("devtools")
# devtools::install_github("ropensci/pathviewr")

library(pathviewr)
library(ggplot2)
library(magrittr)

## Import the example Motive data included in 
## the package
motive_data <-
  read_motive_csv(
    system.file("extdata", "pathviewr_motive_example_data.csv",
                package = 'pathviewr')
  )
```

We'll quickly run through the cleanup pipeline using one of `pathviewr`'s
all-in-one cleaning functions as described in [the Data Import and Cleaning vignette](https://docs.ropensci.org/pathviewr/articles/data-import-cleaning.html).
Since trajectories are defined in the `separate_trajectories()` function, we'll
stop the all-in-one there by setting it and every step after it to `FALSE` so we
can take a closer look at how trajectories are defined. Since all steps are set
to `TRUE` by default, we don't actually need to list them in the function, just
those we want to switch to `FALSE` or those that require additional arguments.

```{r cleanup data}
motive_cleaned <-
  motive_data %>%
  clean_viewr(
    desired_percent = 70,
    rename_viewr_characters = FALSE,
    separate_trajectories = FALSE,
    get_full_trajectories = FALSE
  )

## Quick plot
plot(motive_cleaned$position_length,
     motive_cleaned$position_width,
     asp = 1)
```

### Inspecting the data 
We've cleaned up our data and are ready to link individual data points into
continuous trajectories. Deciding exactly how to define a trajectory will depend
largely on the quality of data collected and the resolution required to answer a
given question about the data. If the data are fairly continuous and high
resolution is required, the default `max_frame_gap = 1` will not allow any frame
gaps--movement must be continuous.

```{r max_frame_gap_1}
motive_mfg1 <- 
  motive_cleaned %>% 
  separate_trajectories(
    max_frame_gap = 1
  )

## Quick plot
plot(motive_mfg1$position_length,
     motive_mfg1$position_width,
     asp = 1, col = as.factor(motive_mfg1$file_sub_traj))

## How many trajectories do we end up with?
length(unique(motive_mfg1$file_sub_traj))
```

Plotting each trajectory reveals that perhaps some trajectories that should be
continuous have been split into two or more separate trajectories because of a
dropped frame or two. This data may require relaxing the stringent no-frame-gaps
requirement in order to link data points that should go together into a single
trajectory. Let's try setting `max_frame_gap = 5`.

```{r max_frame_gap_5}
motive_mfg5 <- 
  motive_cleaned %>% 
  separate_trajectories(
    max_frame_gap = 5
  )

## Quick plot
plot(motive_mfg5$position_length,
     motive_mfg5$position_width,
     asp = 1, col = as.factor(motive_mfg5$file_sub_traj))

## How many trajectories do we end up with?
length(unique(motive_mfg5$file_sub_traj))
```

By increasing the allowable frame gap, more chunks of data have been linked into
trajectories and so our trajectory count has dropped from 335 to 224. Because it
is hard to see all the trajectories piled up on top of each other, let's go
ahead and run this through `get_full_trajectories()` to clean out bits of data
that do not span our region of interest. To inspect each trajectory
individually, we can then use the `plot_viewr_trajectories()` function to take a
closer look at the quality of each trajectory. This can be computationally
expensive depending on your data set.

```{r get_full_trajectories}
motive_mfg5_full <- 
  motive_mfg5 %>% 
  get_full_trajectories(
    span = .6
  )

## Quick plot
plot(motive_mfg5_full$position_length,
     motive_mfg5_full$position_width,
     asp = 1, col = as.factor(motive_mfg5_full$file_sub_traj))

## How many trajectories do we end up with?
length(unique(motive_mfg5_full$file_sub_traj))

## Plot each trajectory
plot_viewr_trajectories(motive_mfg5_full,
                        plot_axes = c("length", "width"),
                        multi_plot = TRUE)
```

### Visualize frame gap choice 

`pathviewr` has several tools to help users decide what frame gap allowances may
be appropriate depending on their data and resolution needs.
`visualize_frame_gap_choice()` runs the `separate_trajectories()` function over
the same data set as many times as the user would like via the `loop` argument,
each time with a different `max_frame_gap` allowance. Each loop represents an
increase in the max frame gap value of 1. For example the default of `loops =
20` will run `separate_trajectories()` over the data set 20 times, with an
increase in the `max_frame_gap` argument of 1 each time.

```{r visualize_frame_gap_choice}
motive_cleaned %>% 
  visualize_frame_gap_choice(
    loops = 20
  )
```

The output of `visualize_frame_gap_choice()` is a tibble and plot of the number
of trajectories after running `separate_trajectories()` with 
`max_frame_gap = 1`, `max_frame_gap = 2`, `max_frame_gap = 3`, etc. 
We can see that as the frame gap allowance increases, more bits of data are
being linked into continuous trajectories and thus the total number of
trajectories decreases. The vertical line on the plot indicates the "elbow" of
the plot or the point at which counts of trajectories appear to stabilize and
increases in the `max_frame_gap` allowance no longer effect the trajectory count
very much.

### Autodetect
Setting `max_frame_gap = "autodetect"` rather than a numeric value uses the
"elbow" of the plots from `visualize_frame_gap_choice()` to guesstimate the best
value(s) for `max_frame_gap`. 

**In addition to automatically selecting the `max_frame_gap` depending on the
data, `autodetect` does so on a per subject basis rather than applying the same
allowable frame gap to all data in the data set since frame gaps can vary
between subjects.** 

The cap on how high the
`max_frame_gap` can go is defined as a proportion of the capture frame rate and
set by the `frame_rate_proportion` argument, which defaults to `.10`.

```{r max_frame_gap_auto}
motive_auto <- 
  motive_cleaned %>% 
  separate_trajectories(
    max_frame_gap = "autodetect",
    frame_rate_proportion = 0.1,
    frame_gap_messaging = TRUE,
    frame_gap_plotting = TRUE
  )

## How many trajectories do we end up with?
length(unique(motive_auto$file_sub_traj))
```

Our sample data has 3 subjects so with `frame_gap_messaging = TRUE`, the
`max_frame_gap` for each subject will be reported (the default is `FALSE`).
`frame_gap_plotting = TRUE` will display the elbow plots for each subject, but
also defaults to `FALSE`.
---
title: "Basics of data import and cleaning in pathviewr"
author: "Vikram B. Baliga"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Basics of data import and cleaning in pathviewr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteDepends{ggplot2}
  %\VignetteDepends{magrittr}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Overview

Raw movement data, including those from motion capture systems, may have a
variety of issues. These raw data often contain noise or artifacts from the
recording session, which may not be easily removed via the recording software
itself. Data may not be organized as “tidy” key-value pairs (making plotting
more difficult), the axes and overall orientation of the environment may not
conform to a standard, and individual movement trajectories may be ill-defined.

`pathviewr` provides functions in R to deal with such problems (i.e.
"cleaning"). This vignette will cover the basics of how to import raw data and
how to clean data to prepare it for visualization and/or statistical analyses.


## What do movement data sets look like?

At minimum, movement data provide information on a subject or object's position
over time. These data are typically supplied in three dimensions (e.g. x, y, z),
with position in each dimension sampled at a particular rate (e.g. 100 Hz).
Different recording software may provide additional features, such as the
ability to track multiple subjects simultaneously, information on subjects'
rotation, tracking of "rigid body" elements, or even the ability to apply Kalman
filters.

A central goal of `pathviewr` is to take data from different sources (so far:
Motive and Flydra), re-organize them into a common format that can be wrangled
in R, clean them up a bit, and get them ready for visualization and/or
statistical analyses. We'll first cover what's included in Motive and in Flydra
data and how `pathviewr` handles these. Should you have data from another
source, our `as_viewr()` function will allow you to bring it into the
`pathviewr` framework.

## Data import via `pathviewr`

Data can be imported via one of three functions:  

- `read_motive_csv()` imports data from `.csv` files that have been exported
from Optitrack's [Motive](https://optitrack.com/software/motive/) software

- `read_flydra_mat()` imports data from  `.mat` files that have been exported
from [Flydra](https://github.com/strawlab/flydra)

- `as_viewr()` can be used to handle data from other sources  

We will showcase examples from each of these methods in this section.

Please feel free to reach out to the `pathviewr` authors via 
[our Github Issues page](https://github.com/ropensci/pathviewr/issues/) 
should you have trouble with any of our data import options. We are happy to
work with you to design custom `read_` functions for file types we have not
encountered ourselves. 

We'll start by loading `pathviewr` and a few of the packages in the `tidyverse`.

```{r package_loading, message=FALSE, warning=FALSE}
## If you do not already have pathviewr installed:
# install.packages("devtools")
# devtools::install_github("ropensci/pathviewr")

library(pathviewr)
library(ggplot2)
library(magrittr)
```

### Motive CSV files

`.csv` files exported from Motive can be imported via `read_motive_csv()`

```{r import_motive}
## Import the Motive example data included in 
## the package

motive_data <-
  read_motive_csv(
    system.file("extdata", "pathviewr_motive_example_data.csv",
                package = 'pathviewr')
  )

## This produces a tibble 
motive_data
``` 

A key thing to note is that these data, as stored in Motive CSVs, are not
"tidy". Each frame occupies one row, but what that also means is that the
rotation and position values for the various subjects take up 24 columns! This
format not only makes plotting data more difficult in base R, `ggplot2`, and
`rgl`, but also makes other aspects of data wrangling more difficult. In a later
step, we will 'gather' these data into key-value pairs so that e.g. all
length-wise position values are in one column, all width-wise are in
another...etc.

Metadata are stored as attributes. We won't go through all of these, but here
are a couple important ones.

``` {r metadata}
## E.g. to see the header of the original file:
attr(motive_data, "header")

## Names of all marked objects:
attr(motive_data, "subject_names_simple")

## Types of data included
attr(motive_data, "data_types_simple")

## Frame rate
attr(motive_data, "frame_rate")
```

Storing such metatdata in the attributes is a key feature of `pathviewr`. These
metadata may not be as immediately as important as the time series of position
or rotation, but they can provide important experimental information such as the
date & time of capture and the units of the position data (here, meters).

### Flydra Matlab files

`.mat` files exported from Flydra can be imported via `read_flydra_mat()`.

Note that you must supply a `subject_name` for Flydra data, as subject names 
are not embedded in the `.mat` files. Only one name can be added and it will
be used throughout the resultant `tibble`.

```{r import_flydra}
## Import the Flydra example data included in 
## the package
flydra_data <-
  read_flydra_mat(
    system.file("extdata", 
                "pathviewr_flydra_example_data.mat",
                package = 'pathviewr'),
    subject_name = "birdie_wooster"
  )

## Similarly, this produces a tibble with important 
## metadata as attributes
flydra_data

attr(flydra_data, "frame_rate")
```

Note that unlike the example Motive data, the Flydra data are already organized
into key-value pairs. Because rotation is not captured by Flydra, such data are
also not included.

### Data from other sources

Data from another format can be converted to a `viewr` object via
`pathviewr::as_viewr()`. Although this function does not handle data import
per se, it allows data that you may already have imported into R as a `tibble`
or `data.frame` to then be reformatted for use with `pathviewr` functions.

We'll run through a quick example with simulated data:

```{r as_viewr}
## Create a dummy data frame with simulated (nonsense) data
df <-
  data.frame(
    frame = seq(1, 100, by = 1),
    time_sec = seq(0, by = 0.01, length.out = 100),
    subject = "birdie_sanders",
    z = rnorm(100),
    x = rnorm(100),
    y = rnorm(100)
  )

## Use as_viewr() to convert it into a viewr object
test <-
  as_viewr(
    df,
    frame_rate = 100,
    frame_col = 1,
    time_col = 2,
    subject_col = 3,
    position_length_col = 5,
    position_width_col = 6,
    position_height_col = 4
  )

## Some metadata are stored as attributes
attr(test, "frame_rate")
```

We also welcome you to request custom data import functions, especially if
`as_viewr()` does not fit your needs. We are interested in expanding our data
import functions to accommodate additional file types. Please feel free to 
[file a request for additional import functions](https://github.com/ropensci/pathviewr/issues/new/choose) 
via our Github Issues page.


## Data cleaning
As noted above, raw data often suffer the following:  
- contain noise or artifacts from the recording session  
- not organized as “tidy” key-value pairs  
- axes and overall orientation of the environment may not conform to a standard  
- individual movement trajectories may be ill-defined

Several functions to clean and wrangle data are available, and we have a
suggested pipeline for how these steps should be handled. The rest of this
vignette will cover these steps.

All of the steps in the suggested pipeline are also covered by two all-in-one
functions: `clean_viewr()` and `import_and_clean_viewr()`. See the section at
the very end of this vignette for details.

And speaking of pipes, all functions in `pathviewr` are pipe-friendly. We will
detail each step separately, but each of the subsequent steps may be piped, e.g.
`data %>% relabel_viewr_axes() %>% gather_tunnel_data()` etc etc

### Relabeling axes, gathering data columns, and trimming outliers
Axis labels (x, y, z) may be applied in arbitrary ways by software. A user might
intuitively think the z axis represents height, but the original software may
label it as the y axis instead.

`relabel_viewr_axes` provides a means to relabel axes with "tunnel_length",
"tunnel_width", and "tunnel_height". **These axis labels will be expected by
subsequent functions, so skipping this step is ill-advised.**

Typically, axes from Motive data will need to be relabled, but axes in data
imported from Flydra will not.

```{r relabel_axes}
motive_relabeled <-
  motive_data %>%
  relabel_viewr_axes(
    tunnel_length = "_z",
    tunnel_width = "_x",
    tunnel_height = "_y",
    real = "_w"
  )

names(motive_relabeled)
```

Akin to the behavior of `dplyr::gather()`, `gather_tunnel_data()` will take all
data from a given session and organize it so that all data of a given type are
within one column, i.e. all position lengths are in `position_length`, as
opposed to separate length columns for each rigid body. **These column names
will be expected by subsequent functions, so skipping this step is also
ill-advised if you are using data from Motive.** Should you have data from 
Flydra, this step should be skipable.

Use `trim_tunnel_outliers()` to remove extreme artifacts and other outlier data.
What this function does is create a (virtual) boundary box according to
user-specification, and any data outside that boundary are removed. For example,
if you know your arena measures 10m x 10m x 10m and your data were calibrated to
range from 0-10m in each dimension, you can be reasonably sure that extreme
values such as 45m on a given axis are bogus. This step is entirely optional,
and should only be used when the user is confident that data outside certain
ranges are artifacts or other bugs. Data outside these ranges are then filtered
out. Best to plot data beforehand and check!!

```{r gather_and_trim}
## First gather and show the new column names
motive_gathered <-
  motive_relabeled %>%
  gather_tunnel_data()

names(motive_gathered)

## Now trim, using ranges we know to safely include data
## and exclude artifacts. Anything outside these ranges 
## will be removed.
motive_trimmed <-
  motive_gathered %>%
  trim_tunnel_outliers(
    lengths_min = 0,
    lengths_max = 3,
    widths_min = -0.4,
    widths_max = 0.8,
    heights_min = -0.2,
    heights_max = 0.5
  )
```

### Standardization of tunnel position and coordinates
The coordinate system of the tunnel itself may require adjustment or
standardization. For example, data collected across different days may show
slight differences in coordinate systems if calibration equipment was not used
in identical ways. Moreover, the user may want to redefine how the coordinate
system itself is defined (i.e. change the location of `(0, 0, 0)` to another
place within the tunnel.

Note that having `(0, 0, 0)` set to the center of the region of interest
(covered in the next section of this vignette) is required for all subsequent
`pathviewr` functions to work.

`pathviewr` offers three main choices for such standardization:  

- `redefine_tunnel_center()`: Sets the location of 0 on any or all axes to a new
location. See the Help page for this function to see the four different methods
by which a user can specify this. No rotation of the tunnel is performed. This
function can be used on both Motive and Flydra data.  

- `standardize_tunnel()`: Use specified landmarks (`subjects` within the `viewr`
object) to rotate and translate the location of a tunnel, setting `(0, 0, 0)` to
the center of the tunnel (centering). For example, in an avian flight tunnel,
perches may be set up on opposite ends of the tunnel and rigid body markers may
be set to them. The positions of these perches can be used as landmarks to
standardize tunnel position. Note that this is typically not possible for Flydra
data, since Flydra data will be imported with only one `subject`.  

- `rotate_tunnel`: Rotate and center a tunnel based on user-defined coordinates
(i.e. similar to `standardize_tunnel()` but for cases where specified landmarks
are not in the data). This function can be used on both Motive and Flydra data.  

Two quick examples will follow, using our Motive and Flydra data:

```{r rotate_example}
## Rotate and center the motive data set:
motive_rotated <-
  motive_trimmed %>% 
  rotate_tunnel(
    perch1_len_min = -0.06,
    perch1_len_max = 0.06,
    perch2_len_min = 2.48,
    perch2_len_max = 2.6,
    perch1_wid_min = 0.09,
    perch1_wid_max = 0.31,
    perch2_wid_min = 0.13,
    perch2_wid_max = 0.35
  )
```

In the above, virtual perches are defined by the user using the arguments shown.
The center of each perch is then found and then the locations of the two perch
centers are then used to 1) set  `(0, 0, 0)` to the point that is equidistant
from the perches (i.e. the middle of the tunnel) and 2) rotate the tunnel about
the height axis so that both perch width coordinates are at 0. This may be
easier to understand through plotting:

```{r rotate_example_plots}
## Quick (base-R) plot of the original data
plot(motive_trimmed$position_length,
     motive_trimmed$position_width,
     asp = 1)

## Quick (base-R) plot of the rotated data
plot(motive_rotated$position_length,
     motive_rotated$position_width,
     asp = 1)
```

Differences due to rotation may be extremely subtle, but the redefining of 
`(0, 0, 0)` to the middle of the tunnel should be clear from contrasting the
axes of the plots.

Flydra data typically do not need to be rotated, so we will instead use
`redefine_tunnel_center()` to adjust the location of `(0, 0, 0)`:

```{r redefine_tunnel_example}
## Re-center the Flydra data set:
flydra_centered <-
  flydra_data %>%
  redefine_tunnel_center(length_method = "middle",
                         height_method = "user-defined",
                         height_zero = 1.44)
```

Here, we are using `length_method = "middle"` to use the middle of the range of
"length" data to set length = 0 (a translation), making no change to the width
axis, and then specifying that height = 0 should be equal to the value `1.44`
from the original data (another translation). Again, plotting may help; note
that this time, we'll plot length x height (instead of width):

```{r redefine_example_plots}
## Quick (base-R) plot of the original data
plot(flydra_data$position_length,
     flydra_data$position_height,
     asp = 1)

## Quick (base-R) plot of the redefined data
plot(flydra_centered$position_length,
     flydra_centered$position_height,
     asp = 1)
```

### Selecting a region of interest
This required step has benefits that are twofold: 1) treatment effects on animal
movement may only manifest over certain regions of the tunnel, and 2) focusing
on a subset of the data makes it easier to define explicit trajectories and run
computations faster.

The region of interest is defined via the function `select_x_percent()`. Once
tunnel coordinates have been standardized (via one of the function in the
previous section), `select_x_percent()` then selects the middle `x` percent
(along the length axis) of the tunnel as the region of interest. For example,
selecting 50 percent would start from the center of the tunnel and move 25% of
the tunnel along positive length and 25% along negative length values to then
select the middle 50% of the tunnel.

Quick examples:
```{r select_x_examples}
## Motive data: select the middle 50% of the tunnel as the region of interest
motive_selected <-
  motive_rotated %>% select_x_percent(50)

## Quick plot:
plot(motive_selected$position_length,
     motive_selected$position_width,
     asp = 1)

## Flydra data:
flydra_selected <-
  flydra_centered %>% select_x_percent(50)

## Quick plot:
plot(flydra_selected$position_length,
     flydra_selected$position_width,
     asp = 1)
```

### Isolating each trajectory
The `pathviewr` standard for defining a trajectory is: continuous movement from
one side of the tunnel to the other over the span of the region of interest.
Note that this definition does not strictly require linear movement from one end
to the other; an animal could make several loops inside the region of interest
within a given trajectory.  

Isolating trajectories is handled via the `separate_trajectories()` function in
`pathviewr`. Note that a region of interest must be selected beforehand via
`select_x_percent()`.  

Because cameras may occasionally drop frames, we allow the user to permit some
relaxation of how stringent the "continuous movement" criterion is. This is
handled via the `max_frame_gap` argument within `separate_trajectories()`. For
more details, please see [the vignette Managing frame gaps with pathviewr](https://docs.ropensci.org/pathviewr/articles/managing-frame-gaps.html).  

In our Motive example, we'll use the automated feature built into the function
to guesstimate the best `max_frame_gap` allowed. When frame gaps larger than
`max_frame_gap` are encountered, the function will force the defining of a new
trajectory. But if frame gaps smaller than `max_frame_gap` are encountered,
keeping observations within the same trajectory is permitted.

In the Flydra example, we'll simply set `max_frame_gap` to `1` so that no frame
gaps are allowed (movement must be continuous with no dropped frames).

```{r separate_examples}
motive_labeled <-
  motive_selected %>% 
  separate_trajectories(max_frame_gap = "autodetect")

flydra_labeled <-
  flydra_selected %>% 
  separate_trajectories(max_frame_gap = 1)
```

### Retain only complete trajectories
Now that trajectories have been isolated and labeled, the final cleaning step is
to retain only the trajectories that completely span from one end of the region
of interest to the other.

This final step is handled via `pathviewr`'s `get_full_trajectories()`.

There is a built-in "fuzziness" feature: because trajectories may not have
observations exactly at the beginning or the end of the region of interest, it
may be necessary to allow trajectories to be slightly shorter than the range of
the selected region of interest. The `span` parameter of this function handles
this. By supplying a numeric proportion from `0` to `1`, a user may allow
trajectories to span that proportion of the selected region. For example,
setting `span = 0.95` will keep all trajectories that span 95% of the length of
the selected region of interest. Setting `span = 1` (not recommended) will
strictly keep trajectories that start and end at the **exact** cut-offs of the
selected region of interest.

For these reasons, `span`s of `0.99` to `0.95` are generally recommended. The
best choice ultimately depends on your capture frame rate as well as your own
judgment. Should you desire to set it lower (which you can), you may instead
consider using a smaller region of interest (i.e. set the `desired_percent`
parameter in `select_x_percent()` to be lower).

```{r get_full_examples}
## Motive
motive_full <-
  motive_labeled %>% 
  get_full_trajectories(span = 0.95)

plot(motive_full$position_length,
     motive_full$position_width,
     asp = 1, col = as.factor(motive_full$file_sub_traj))

## Flydra
flydra_full <-
  flydra_labeled %>% 
  get_full_trajectories(span = 0.95)

plot(flydra_full$position_length,
     flydra_full$position_width,
     asp = 1, col = as.factor(flydra_full$file_sub_traj))
```

### All-in-one cleaning functions
All of the above steps can also be done by using `pathviewr`'s designated 
all-in-one functions. `import_and_clean_viewr()` imports raw data and allows
the user to run through all of the cleaning steps previously covered in this
vignette. `clean_viewr()` does the same on any object already imported into the
R environment (i.e. it skips data import).

For both of these functions, all of the cleaning steps are set to `TRUE` by 
default, but may be turned off by using `FALSE`. Argument names correspond to 
standalone functions in `pathviewr`, and if the user wants to use non-default values for corresponding arguments, they should also be supplied for any steps 
that are set to `TRUE`. 

For example, if the user keeps `select_x_percent = TRUE` as an argument in
`clean_viewr()`, the `select_x_percent()` function is run internally. This means
that should the user desire to select a region of interest that does not match
the default value of 33 percent, an additional argument should be supplied to
`clean_viewr()` as if it were being supplied to `select_x_percent()` itself, 
e.g.: `desired_percent = 50`.

All additional arguments should be written out fully and explicitly.

Here's an example using what we previously covered:

```{r all-in-one}
motive_allinone <-
  motive_data %>%
  clean_viewr(
    relabel_viewr_axes = TRUE,
    gather_tunnel_data = TRUE,
    trim_tunnel_outliers = TRUE,
    standardization_option = "rotate_tunnel",
    select_x_percent = TRUE,
    desired_percent = 50,
    rename_viewr_characters = FALSE,
    separate_trajectories = TRUE,
    max_frame_gap = "autodetect",
    get_full_trajectories = TRUE,
    span = 0.95
  )

## Quick plot
plot(motive_allinone$position_length,
     motive_allinone$position_width,
     asp = 1, col = as.factor(motive_allinone$file_sub_traj))
```


That's all!

🐢
---
title: "Estimating visual perceptions from tracking data"
author: "Eric R. Press" 
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Estimating visual perceptions from tracking data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteDepends{ggplot2}
  %\VignetteDepends{magrittr}
  %\VignetteEncoding{UTF-8}
---
```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

Studies of visually guided locomotion in birds, insects or even humans often
involve data gathered from motion capture technologies such as Optitrack's
[(Motive)](https://optitrack.com/software/motive/), or the Straw Lab's
[(Flydra)](https://github.com/strawlab/flydra). For these experiments, it is
important to understand how visual stimuli influence behavior. Although is it
not possible to measure how subjects directly perceive visual stimuli, it is
possible to use motion capture data to calculate estimates of stimulus
properties as they are perceived by the subject. With the tools available in
pathviewr, researchers in ethology or psychology may provide analyses of both
stimulus and response under natural locomotor conditions.
 
*Inherent to estimations of visual perceptions, we make several assumptions,
which will be discussed below. We welcome suggestions and aim to address any
assumptions that limit the accuracy of these estimates*.
  
To bridge the gap between objective measures of subject position and estimates
of subjective stimulus perception, we can begin by calculating the angle that a
visual pattern subtends on the subject's retina: the visual angle (θ). Visual
angles can be used to calculate aspects of image motion such as the rate of
visual expansion [(Dakin *et al.*,
2016)](https://www.pnas.org/content/113/31/8849/). For a detailed
review of different forms of visual motion and how they're processed by the
brain, see [(Frost, 2010)](https://www.karger.com/Article/Abstract/314284).
  
Visual angles can be calculated provided there is information about the physical
size of the pattern and its distance from the subject's retina. Because
researchers can control or measure the size of a pattern, and we can calculate
the distance between the subject and pattern using motion capture data, we can
further calculate the visual angle produced by patterns in the visual scene.
Therefore, we first need to calculate the distance from the subject's retina to
the pattern.
  
*Currently, we assume the subject's gaze is directly frontal or lateral to the
face, in effect estimating image properties at single points in the frontal and
lateral fields of view, respectively. We currently calculate distances to the
center of the subject's head rather than the position of the retina; future
versions of `pathviewr` will include features addressing these limitations, such
as including head orientation information and eye position relative to the
center of the subject's head*.
  
```{r, echo=FALSE, out.width="100%", fig.cap="Visual angles can be calculated using the size of a visual pattern (`stim_param`) and the distance to the pattern. Larger patterns at shorter distances produce larger visual angles. For dot stimuli, visual angles can be calculated independent of stimulus orientation."}
knitr::include_graphics("https://raw.githubusercontent.com/ropensci/pathviewr/dc728871873b92dc8412bdcb5ef58031457f6788/images/stim_param_angle.jpeg")
```

## Loading packages
To begin, let's load `pathviewr` and a few of the packages in `tidyverse`.

```{r package_loading, message=FALSE, warning=FALSE}
library(pathviewr)
library(ggplot2)
library(magrittr)
```


## Data preparation
Data objects must be prepared as described in [the Data Import and Cleaning
vignette](https://docs.ropensci.org/pathviewr/articles/data-import-cleaning.html)
pipeline prior to their use with these functions. For a detailed description of
the importing and cleaning functions, and when to use them, please see the
linked vignette.

Let's work with a few example datasets included in the package.
`pathviewr_motive_example_data.csv` is a `.csv` file exported from Motive.
`pathviewr_flydra_example_data.mat` is a `.mat` file exported from Flydra. For
more coarse-grained data cleaning tasks, `pathviewr` contains an all-in-one
cleaning function `clean_viewr`. We will use this function for the following
examples.


```{r}
## Import motive data set
motive_data <- # import
  read_motive_csv(
    system.file("extdata", "pathviewr_motive_example_data.csv",
                package = 'pathviewr')
  )

## Clean motive data set
motive_full <-
  motive_data %>%
  clean_viewr(
    relabel_viewr_axes = TRUE,
    gather_tunnel_data = TRUE,
    trim_tunnel_outliers = TRUE,
    standardization_option = "rotate_tunnel",
    select_x_percent = TRUE,
    desired_percent = 50,
    rename_viewr_characters = FALSE,
    separate_trajectories = TRUE,
    max_frame_gap = "autodetect",
    get_full_trajectories = TRUE,
    span = 0.95
  )
```
```{r}
## Import flydra data set
flydra_data <- 
  read_flydra_mat(
    system.file("extdata", "pathviewr_flydra_example_data.mat",
                package = 'pathviewr'),
    subject_name = "birdie_wooster")

## Clean flydra data set
flydra_full <- 
  flydra_data %>%
  clean_viewr(
    relabel_viewr_axes = FALSE,
    gather_tunnel_data = FALSE,
    trim_tunnel_outliers = FALSE,
    standardization_option = "redefine_tunnel_center",
    length_method = "middle",
    height_method = "user-defined",
    height_zero = 1.44,
    get_velocity = FALSE,
    select_x_percent = TRUE,
    desired_percent = 60,
    rename_viewr_characters = FALSE,
    separate_trajectories = TRUE,
    get_full_trajectories = TRUE
  )
```



## Add experiment information with `insert_treatments()`
Now that our objects have been cleaned, we must use `insert_treatments()` to add
information about the experiments, including relevant properties of the visual
stimulus and experimental tunnel that are necessary for calculating visual
perceptions.

*`pathviewr` currently supports rectangular (box) or V-shaped experimental
tunnels, though we are interested in including additional tunnel
configurations*.

If you would like to request any additional features such as other tunnel configurations, please create an issue [(here)](https://github.com/ropensci/pathviewr/issues/new/choose)

#### V-shaped tunnel example
The data within `motive_full` were collected from birds flying through a 3m long
V-shaped tunnel in which the height of the origin `(0,0,0)` was set to the
height of the perches that were 0.3855m above the vertex, which was angled at
90 degree. The visual stimuli on the positive and negative walls of the tunnel (where
`position_width` values > 0 and < 0, respectively) were stationary dot-fields.
Each dot was of size 0.05m in diameter. The visual stimuli on the positive and
negative end walls of the tunnel (where `position_length` > 0 and < 0,
respectively) were dot fields with dots 0.1m in diameter. This treatment was
referred to as `"latB"`.

Therefore we will use the following code:
```{r}
motive_treatments <- 
  motive_full %>% 
  insert_treatments(tunnel_config = "v",
                    perch_2_vertex = 0.3855,
                    vertex_angle = 90,
                    tunnel_length = 2,
                    stim_param_lat_pos = 0.05,
                    stim_param_lat_neg = 0.05,
                    stim_param_end_pos = 0.1,
                    stim_param_end_neg = 0.1,
                    treatment = "latB")
names(motive_treatments)
```
`motive_treatments` now has the variables `tunnel_config`, `perch_2_vertex`,
`vertex_angle`, `tunnel_length`, `stim_param_lat_pos`, `stim_param_lat_neg`,
`stim_param_end_pos`, and `stim_param_end_neg` which are needed to calculate
visual angles. The variable `treatment` has also been included and all of this
information has been stored in the object's metadata.


#### Box-shaped tunnel example
The data within `flydra_full` were collected from birds flying in a 1 x 3m
rectangular tunnel (a box). The visual stimuli were the same as in the motive
example.
```{r}
flydra_treatments <-
  flydra_full %>%
  insert_treatments(tunnel_config = "box",
                    tunnel_width = 1,
                    tunnel_length = 3,
                    stim_param_lat_pos = 0.05,
                    stim_param_lat_neg = 0.05,
                    stim_param_end_pos = 0.1,
                    stim_param_end_neg = 0.1,
                    treatment = "latB")
```
`flydra_treatments` similarly has the variables `tunnel_config`, `tunnel_width`,
`tunnel_length`, `stim_param_lat_pos`, `stim_param_lat_neg`,
`stim_param_end_pos`, `stim_param_end_neg` and `treatment`.



## Calculating visual angles
### Start by calculating distances to visual stimuli 

To estimate the visual angles perceived by the subject is it moves through the
tunnel, we first need to calculate the distance between the subject and the
visual stimuli. For this, we will use `calc_min_dist_v` or `calc_min_dist_box`
depending on the configuration of the experimental tunnel. These functions
calculate the minimum distance between the subject and the surface displaying a
visual pattern, therefore maximizing the visual angles.  

For V-shaped tunnels, several internal calculations are required and can be
added to the output object with `simplify_output = FALSE`. Otherwise, the
minimum distance are computed to the lateral walls and the end wall to which the
subject is facing.

```{r motive_min_dist_pos, fig.height=4, fig.width=7}
motive_min_dist <- 
  motive_treatments %>% 
  calc_min_dist_v(simplify_output = FALSE)

## Display minimum distances to the positive lateral walls 
## Viewpoint is from the end of the tunnel
motive_min_dist %>% 
  ggplot(aes(x = position_width, y = position_height)) +
  geom_point(aes(color = min_dist_pos), size = 2, shape = 1) +
  coord_fixed() +
  theme_classic() +
  geom_segment(aes(x = 0,         # positive wall
                   y = -0.3855,
                   xend = 0.5869,
                   yend = 0.2014)) +
  geom_segment(aes(x = 0,         # negative wall
                   y = -0.3855,
                   xend = -0.5869,
                   yend = 0.2014))
```

```{r flydra_min_dist_end, fig.height=4, fig.width=7}
flydra_min_dist <- 
  flydra_treatments %>% 
  calc_min_dist_box()

## Display minimum distances to the end walls
## Viewpoint is from above the tunnel
flydra_min_dist %>% 
  ggplot(aes(x = position_length, y = position_width)) +
  geom_point(aes(color = min_dist_end), size = 2, shape = 1) +
  coord_fixed() +
  theme_classic() +
  geom_segment(aes(x = -1,         # negative wall
                   y = -0.5,
                   xend = 1,
                   yend = -0.5)) +
  geom_segment(aes(x = -1,         # positive wall
                   y = 0.5,
                   xend = 1,
                   yend = 0.5))
```


### Now get visual angles 

```{r motive_vis_angle_pos, fig.height=4, fig.width=7}
motive_vis_angle <- 
  motive_min_dist %>% 
  get_vis_angle()

## Visualize the angles produced from stimuli on the positive wall
## Viewpoint is from the end of the tunnel
motive_vis_angle %>% 
  ggplot(aes(x = position_width, y = position_height)) +
  geom_point(aes(color = vis_angle_pos_deg), size = 2, shape = 1) +
  coord_fixed()+
  theme_classic() +
  geom_segment(aes(x = 0,         # positive wall
                   y = -0.3855,
                   xend = 0.5869,
                   yend = 0.2014)) +
  geom_segment(aes(x = 0,         # negative wall
                   y = -0.3855,
                   xend = -0.5869,
                   yend = 0.2014))
```

Notice larger visual angles as the subject approaches the positive wall. 

```{r flydra_vis_angle_end, fig.height=4, fig.width=7}
flydra_vis_angle <- 
  flydra_min_dist %>% 
  get_vis_angle()

## Visualize the angles produced by stimuli on the end walls
## Viewpoint is from above the tunnel
flydra_vis_angle %>% 
  ggplot(aes(x = position_length, y = position_width)) +
  geom_point(aes(color = vis_angle_end_deg), size = 2, shape = 1) +
  coord_fixed() +
  theme_classic() +
  geom_segment(aes(x = -1,        # negative wall
                   y = -0.5,
                   xend = 1,
                   yend = -0.5)) +
  geom_segment(aes(x = -1,        # positive wall
                   y = 0.5,
                   xend = 1,
                   yend = 0.5))
```

Notice larger visual angles as the subject approaches the end wall that it is
moving towards.

## Calculating spatial frequency
With visual angles, we can now determine the spatial frequency of a visual
pattern as it is perceived by the subject. Spatial frequency refers to the size
of the pattern in visual space. It's often reported as the number of cycles of a
visual pattern per 1 degree of the visual field (cycles/degree). Here, we will define
a cycle length as the length used for the `stim_param`. For a given distance
from the subject, a larger visual pattern produces a smaller spatial frequency,
whereas a smaller visual pattern produces a larger spatial frequency.

To calculate the spatial frequency of the visual stimuli as perceived by the
subject some distance from the stimuli, we will use `get_sf()`.

```{r motive_sf_pos, fig.height=4, fig.width=7}
motive_sf <- 
  motive_vis_angle %>%
  get_sf()

## Visualize the spatial frequency of the stimulus on the positive wall 
## point is from the end of the tunnel
motive_sf %>% 
  ggplot(aes(x = position_width, y = position_height)) +
  geom_point(aes(color = sf_pos), size = 2, shape = 1) +
  coord_fixed()+
  theme_classic() +
  geom_segment(aes(x = 0,         # positive wall
                   y = -0.3855,
                   xend = 0.5869,
                   yend = 0.2014)) +
  geom_segment(aes(x = 0,         # negative wall
                   y = -0.3855,
                   xend = -0.5869,
                   yend = 0.2014))
```

Notice the spatial frequency increases the further the subject recedes from
positive wall.

```{r flydra_sf_end, fig.height=4, fig.width=7}
flydra_sf <- 
  flydra_vis_angle %>% 
  get_sf()

## Visualize the spatial frequency of the stimulus on the end walls
## Viewpoint is from above the tunnel
flydra_sf %>% 
  ggplot(aes(x = position_length, y = position_width)) +
  geom_point(aes(color = sf_end), size = 2, shape = 1) +
  coord_fixed() +
  theme_classic() +
  geom_segment(aes(x = -1,        # negative wall
                   y = -0.5,
                   xend = 1,
                   yend = -0.5)) +
  geom_segment(aes(x = -1,        # positive wall
                   y = 0.5,
                   xend = 1,
                   yend = 0.5))
```

Notice the spatial frequency decreases as the subject approaches the end wall
that it is moving towards.

### Stay tuned as additional features for image motion estimation are coming soon!
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_in_one_functions.R
\name{clean_viewr}
\alias{clean_viewr}
\title{All-in-one function to clean imported objects}
\usage{
clean_viewr(
  obj_name,
  relabel_viewr_axes = TRUE,
  gather_tunnel_data = TRUE,
  trim_tunnel_outliers = TRUE,
  standardization_option = "rotate_tunnel",
  get_velocity = TRUE,
  select_x_percent = TRUE,
  rename_viewr_characters = FALSE,
  separate_trajectories = TRUE,
  get_full_trajectories = TRUE,
  fill_traj_gaps = FALSE,
  ...
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{relabel_viewr_axes}{default TRUE, should axes be relabeled?}

\item{gather_tunnel_data}{default TRUE, should tunnel data be gathered?}

\item{trim_tunnel_outliers}{default TRUE, outliers be trimmed?}

\item{standardization_option}{default "rotate_tunnel"; which standardization
option should be used? See Details for more.}

\item{get_velocity}{default TRUE, should velocity be computed?}

\item{select_x_percent}{default TRUE, should a region of interest be
selected?}

\item{rename_viewr_characters}{default FALSE, should subjects be renamed?}

\item{separate_trajectories}{default TRUE, should trajectories be defined?}

\item{get_full_trajectories}{default TRUE, should only full trajectories be
retained?}

\item{fill_traj_gaps}{default FALSE, should gaps in trajectories be filled?}

\item{...}{Additional arguments passed to any of the corresponding functions}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) that has passed
through several \code{pathviewr} functions as desired by the user,
resulting in data that have been cleaned and ready for analyses.
}
\description{
For an imported viewr object, run through the cleaning pipeline as desired
}
\details{
Each argument corresponds to a standalone function in
\code{pathviewr}. E.g. the parameter \code{relabel_viewr_axes} allows a
user to choose whether \code{pathviewr::relabel_viewr_axes()} is run
internally. Should the user desire to use any non-default parameter values
for any functions included here, they should be supplied to this function
as additional arguments formatted exactly as they would appear in their
corresponding function(s). E.g. if the "autodetect" feature in
\code{pathviewr::separate_trajectories()} is desired, add an argument
\code{max_frame_gap = "autodetect"} to the arguments supplied to this
function.
}
\examples{
library(pathviewr)

## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

motive_full <-
  motive_data \%>\%
  clean_viewr(desired_percent = 50,
              max_frame_gap = "autodetect",
              span = 0.95)

## Alternatively, used the import_and_clean_viewr()
## function to combine these steps
motive_import_and_clean <-
  import_and_clean_viewr(
    file_name = system.file("extdata", "pathviewr_motive_example_data.csv",
                            package = 'pathviewr'),
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )
}
\seealso{
Other all in one functions: 
\code{\link{import_and_clean_viewr}()}
}
\author{
Vikram B. Baliga
}
\concept{all in one functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathviewr.R
\docType{package}
\name{pathviewr-package}
\alias{pathviewr}
\alias{pathviewr-package}
\title{pathviewr: Wrangle, Analyze, and Visualize Animal Movement Data}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Tools to import, clean, and visualize movement data,
    particularly from motion capture systems such as Optitrack's 
    'Motive', the Straw Lab's 'Flydra', or from other sources. We provide 
    functions to remove artifacts, standardize tunnel position and tunnel 
    axes, select a region of interest, isolate specific trajectories, fill
    gaps in trajectory data, and calculate 3D and per-axis velocity. For 
    experiments of visual guidance, we also provide functions that use 
    subject position to estimate perception of visual stimuli.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/pathviewr/}
  \item \url{https://docs.ropensci.org/pathviewr/}
  \item Report bugs at \url{https://github.com/ropensci/pathviewr/issues/}
}

}
\author{
\strong{Maintainer}: Vikram B. Baliga \email{vbaliga87@gmail.com} (\href{https://orcid.org/0000-0002-9367-8974}{ORCID})

Authors:
\itemize{
  \item Melissa S. Armstrong \email{melissa.armstrong@gmail.com} (\href{https://orcid.org/0000-0002-3059-0094}{ORCID})
  \item Eric R. Press \email{epress12@gmail.com} (\href{https://orcid.org/0000-0002-1944-3755}{ORCID})
}

Other contributors:
\itemize{
  \item Anne-Sophie Bonnet-Lebrun [reviewer]
  \item Marco Sciaini [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{redefine_tunnel_center}
\alias{redefine_tunnel_center}
\title{"Center" the tunnel data, i.e. translation but no rotation}
\usage{
redefine_tunnel_center(
  obj_name,
  axes = c("position_length", "position_width", "position_height"),
  length_method = c("original", "middle", "median", "user-defined"),
  width_method = c("original", "middle", "median", "user-defined"),
  height_method = c("original", "middle", "median", "user-defined"),
  length_zero = NA,
  width_zero = NA,
  height_zero = NA,
  ...
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{axes}{Names of axes to be centered}

\item{length_method}{Method for length}

\item{width_method}{Method for width}

\item{height_method}{Method for height}

\item{length_zero}{User-defined value}

\item{width_zero}{User-defined value}

\item{height_zero}{User-defined value}

\item{...}{Additional arguments passed to/from other pathviewr functions}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) in which data have
been translated according to the user's inputs, generally with \code{(0, 0,
  0,)} being relocated to the center of the tunnel.
}
\description{
Redefine the center \code{(0, 0, 0,)} of the tunnel data via translating
positions along axes.
}
\details{
For each \code{_method} argument, there are four choices of how centering is
handled: 1) "original" keeps axis as is -- this is how width and (possibly)
height should be handled for flydra data; 2) "middle" is the middle of the
range of data: (min + max) / 2; 3) "median" is the median value of data on
that axis. Probably not recommended; and 4) "user-defined" lets the user
customize where the (0, 0, 0) point in the tunnel will end up. Each
\code{_zero} argument is subtracted from its corresponding axis' data.
}
\examples{
## Import the Flydra example data included in
## the package
flydra_data <-
  read_flydra_mat(
    system.file("extdata",
                "pathviewr_flydra_example_data.mat",
                package = 'pathviewr'),
    subject_name = "birdie_wooster"
  )

## Re-center the Flydra data set.
## Width will be untouched
## Length will use the "middle" definition
## And height will be user-defined to be
## zeroed at 1.44 on the original axis
flydra_centered <-
  flydra_data \%>\%
  redefine_tunnel_center(length_method = "middle",
                         height_method = "user-defined",
                         height_zero = 1.44)
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}

Other tunnel standardization functions: 
\code{\link{rotate_tunnel}()},
\code{\link{standardize_tunnel}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
\concept{tunnel standardization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{find_curve_elbow}
\alias{find_curve_elbow}
\title{Find the "elbow" of a curve.}
\usage{
find_curve_elbow(data_frame, export_type = "row_num", plot_curve = FALSE)
}
\arguments{
\item{data_frame}{A two-column data frame (numeric entries only)}

\item{export_type}{If "row_num" (the default), the row number of the elbow
point is returned. If anything else, the entire row of the original data
frame is returned.}

\item{plot_curve}{Default FALSE; should the curve be plotted?}
}
\value{
If \code{export_type} is \code{row_num} the row number of the elbow
point is returned. If anything else is used for that argument, the entire
row of the original data frame on which the "elbow" is located is returned.
If \code{plot_curve} is \code{TRUE}, the curve is plotted along with a
vertical line drawn at the computed elbow point.
}
\description{
For bivariate data that show monotonic decreases (e.g. plots of trajectory
count vs. frame gap allowed, or scree plots from PCAs), this function will
find the "elbow" point. This is done by drawing an (imaginary) line between
the first observation and the final observation. Then, the distance between
that line and each observation is calculated. The "elbow" of the curve is the
observation that maximizes this distance.
}
\examples{
df <- data.frame(x = seq(1:10),
                 y = 1/seq(1:10))
plot(df)
find_curve_elbow(df, plot_curve = TRUE)
}
\seealso{
Other mathematical functions: 
\code{\link{calc_min_dist_v}()},
\code{\link{deg_2_rad}()},
\code{\link{get_2d_angle}()},
\code{\link{get_3d_angle}()},
\code{\link{get_3d_cross_prod}()},
\code{\link{get_dist_point_line}()},
\code{\link{get_traj_velocities}()},
\code{\link{get_velocity}()},
\code{\link{rad_2_deg}()}
}
\author{
Vikram B. Baliga
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_functions.R
\name{import_batch}
\alias{import_batch}
\title{Batch import of files for either Motive or Flydra (but not a mix of both)}
\usage{
import_batch(
  file_path_list,
  import_method = c("flydra", "motive"),
  file_id = NA,
  subject_name = NULL,
  frame_rate = NULL,
  simplify_marker_naming = TRUE,
  import_messaging = FALSE,
  ...
)
}
\arguments{
\item{file_path_list}{A list of file paths}

\item{import_method}{Either "flydra" or "motive"}

\item{file_id}{Optional}

\item{subject_name}{For Flydra, the assigned subject name}

\item{frame_rate}{For Flydra, the assigned frame rate}

\item{simplify_marker_naming}{default TRUE; for Motive, whether marker naming
should be simplified}

\item{import_messaging}{default FALSE; should this function report each time
a file has been imported?}

\item{...}{Additional arguments (may remove this if needed)}
}
\value{
A list of viewr objects (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}).
}
\description{
Batch import of files for either Motive or Flydra (but not a mix of both)
}
\details{
Refer to \code{read_motive_csv()} and \code{read_flydra_mat()} for
details of data import methods.
}
\examples{
## Since we only have one example file of each type provided
## in pathviewr, we will simply import the same example multiple
## times to simulate batch importing. Replace the contents of
## the following list with your own list of files to be imported.

## Make a list of the same example file 3x
import_list <-
  c(rep(
    system.file("extdata", "pathviewr_motive_example_data.csv",
                package = 'pathviewr'),
    3
  ))

## Batch import
motive_batch_imports <-
  import_batch(import_list,
               import_method = "motive",
               import_messaging = TRUE)

## Batch cleaning of these imported files
## via clean_viewr_batch()
motive_batch_cleaned <-
  clean_viewr_batch(
    file_announce = TRUE,
    motive_batch_imports,
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )

## Alternatively, use import_and_clean_batch() to
## combine import with cleaning on a batch of files
motive_batch_import_and_clean <-
  import_and_clean_batch(
    import_list,
    import_method = "motive",
    import_messaging = TRUE,
    motive_batch_imports,
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )

## Each of these lists of objects can be bound into
## one viewr object (i.e. one tibble) via
## bind_viewr_objects()
motive_bound_one <-
  bind_viewr_objects(motive_batch_cleaned)

motive_bound_two <-
  bind_viewr_objects(motive_batch_import_and_clean)

## Either route results in the same object ultimately:
identical(motive_bound_one, motive_bound_two)
}
\seealso{
Other data import functions: 
\code{\link{as_viewr}()},
\code{\link{import_and_clean_batch}()},
\code{\link{read_flydra_mat}()},
\code{\link{read_motive_csv}()}

Other batch functions: 
\code{\link{bind_viewr_objects}()},
\code{\link{clean_viewr_batch}()},
\code{\link{import_and_clean_batch}()}
}
\author{
Vikram B. Baliga
}
\concept{batch functions}
\concept{data import functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{standardize_tunnel}
\alias{standardize_tunnel}
\title{Rotate and center a tunnel based on landmarks}
\usage{
standardize_tunnel(
  obj_name,
  landmark_one = "perch1",
  landmark_two = "perch2",
  ...
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"} that has been passed
through \code{relabel_viewr_axes()} and \code{gather_tunnel_data()} (or is
structured as though it has been passed through those functions).}

\item{landmark_one}{Subject name of the first landmark}

\item{landmark_two}{Subject name of the second landmark}

\item{...}{Additional arguments passed to/from other pathviewr functions}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) in which data have
been rotated according to the positions of the landmarks in the data.
}
\description{
Similar to \code{rotate_tunnel()}; rotate and center tunnel data based on
landmarks (specific subjects in the data).
}
\details{
The center point of the tunnel is estimated as the point between the
two landmarks. It is therefore recommended that \code{landmark_one} and
\code{landmark_two} be objects that are placed on opposite ends of the
tunnel (e.g. in an avian flight tunnel, these landmarks may be perches that
are placed at the extreme ends). The angle between landmark_one,
tunnel_center_point, and arbitrary point along the length axis
(tunnel_center_point - 1 on length) is estimated. That angle is then used
to rotate the data, again only in the length and width dimensions. Height
is standardized by average landmark height; values greater than 0 are above
the landmarks and values less than 0 are below the landmark level.
}
\section{Warning}{

The \code{position_length} values of landmark_one MUST be less than
the \code{position_length} values of landmark_two; otherwise,
the rotation will apply to a mirror-image of the tunnel
}

\examples{
## Example data that would work with this function are
## not included in this version of pathviewr. Please
## contact the package authors for further guidance
## should you need it.
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}

Other tunnel standardization functions: 
\code{\link{redefine_tunnel_center}()},
\code{\link{rotate_tunnel}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
\concept{tunnel standardization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{fill_traj_gaps}
\alias{fill_traj_gaps}
\title{Interpolate gaps within trajectories}
\usage{
fill_traj_gaps(
  obj_name,
  loess_degree = 1,
  loess_criterion = c("aicc", "gcv"),
  loess_family = c("gaussian", "symmetric"),
  loess_user_span = NULL
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}. Trajectories must be
predefined (i.e. via \code{separate_trajectories()}).}

\item{loess_degree}{See "degree" argument of fANCOVA::loess.as()}

\item{loess_criterion}{See "criterion" argument of fANCOVA::loess.as()}

\item{loess_family}{See "family" argument of fANCOVA::loess.as()}

\item{loess_user_span}{See "user.span" argument of fANCOVA::loess.as()}
}
\value{
A viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"} that now includes new
observations (rows) as a result of interpolation to fill in missing data. A
new column \code{gaps_filled} is added to the data to indicate original
data ("No") vs data that have been inserted to fill gaps ("Yes").
}
\description{
Use LOESS smoothing to fill in gaps of missing data within trajectories in
a viewr object
}
\details{
It is strongly recommended that the input viewr object be "cleaned"
via \code{select_x_percent()} -> \code{separate_trajectories()} ->
\code{get_full_trajectories()} prior to using this function. Doing so will
ensure that only trajectories with minor gaps will be used in your
analyses. This function will then enable you to interpolate missing data in
those minor gaps.

Interpolation is handled by first fitting a series of LOESS regressions
(via \code{fANCOVA::loess.as()}). In each regression, a position axis (e.g.
\code{position_length}) is regressed against \code{frame} (\code{frame} is
x-axis). From that relationship, values of missing position data are
determined and then inserted into the original data set.

See \link[fANCOVA]{loess.as} for further details on parameters.
}
\examples{
library(pathviewr)

## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

## Clean, isolate, and label trajectories
motive_full <-
  motive_data \%>\%
  clean_viewr(desired_percent = 50,
              max_frame_gap = "autodetect",
              span = 0.95)

## Interpolate missing data via this function
motive_filling <-
 motive_full \%>\%
 fill_traj_gaps()

## plot all trajectories (before)
plot_viewr_trajectories(motive_full, multi_plot = TRUE)
## plot all trajectories(after)
plot_viewr_trajectories(motive_filling, multi_plot = TRUE)
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{rad_2_deg}
\alias{rad_2_deg}
\title{Convert radians to degrees}
\usage{
rad_2_deg(rad)
}
\arguments{
\item{rad}{Radians (a numeric of any length >= 1)}
}
\value{
The angle(s) in degrees (as a numeric vector of the same length)
}
\description{
Convert radians to degrees
}
\examples{
## One input
rad_2_deg(pi/2)

## Multiple inputs
rad_2_deg(c(pi / 2, pi, 2 * pi))
}
\seealso{
Other mathematical functions: 
\code{\link{calc_min_dist_v}()},
\code{\link{deg_2_rad}()},
\code{\link{find_curve_elbow}()},
\code{\link{get_2d_angle}()},
\code{\link{get_3d_angle}()},
\code{\link{get_3d_cross_prod}()},
\code{\link{get_dist_point_line}()},
\code{\link{get_traj_velocities}()},
\code{\link{get_velocity}()}
}
\author{
Vikram B. Baliga
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{get_3d_cross_prod}
\alias{get_3d_cross_prod}
\title{Compute the cross product of two 3D vectors}
\usage{
get_3d_cross_prod(v1, v2)
}
\arguments{
\item{v1}{First vector, as c(x,y,z)}

\item{v2}{Second vector, as c(x,y,z)}
}
\value{
A vector of length 3 that describes the cross-product
}
\description{
Compute the cross product of two 3D vectors
}
\examples{
v1 <- c(1, 1, 3)
v2 <- c(3, 1, 3)
get_3d_cross_prod(v1, v2)
}
\seealso{
Other mathematical functions: 
\code{\link{calc_min_dist_v}()},
\code{\link{deg_2_rad}()},
\code{\link{find_curve_elbow}()},
\code{\link{get_2d_angle}()},
\code{\link{get_3d_angle}()},
\code{\link{get_dist_point_line}()},
\code{\link{get_traj_velocities}()},
\code{\link{get_velocity}()},
\code{\link{rad_2_deg}()}
}
\author{
Vikram B. Baliga
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_by_subject}
\alias{plot_by_subject}
\title{Plot trajectories and density plots of position by subject}
\usage{
plot_by_subject(obj_name, col_by_treat = FALSE, ...)
}
\arguments{
\item{obj_name}{A viewr object (a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) that has been passed
through \code{separate_trajectories()} or \code{get_full_trajectories()}.}

\item{col_by_treat}{If multiple treatments or sessions, color data per
treatment or session. Treatments must be levels in a column named
\code{treatment}.}

\item{...}{Additional arguments passed to/from other pathviewr functions.}
}
\value{
A "bird's eye view" plot and an "elevation view" plot.
}
\description{
Plots all trajectories and generates density plots of position by subject
from elevation and bird's eye views.
}
\details{
The input viewr object should have passed through
\code{separate_trajectories()} or \code{get_full_trajectories()}.
Optionally, treatments should have been added as levels in a column named
\code{treatment}. Two plots will be produced, one from a "bird's eye view"
of width against length and one from an "elevation view" of height against
length. All trajectories will be plotted on a per subject basis as well as
density plots of width or height depending on the view.
\code{col_by_treat = TRUE}, data will be plotted by color according to
treatment in both the trajectory plots and the density plots.
}
\examples{
library(pathviewr)
library(ggplot2)
library(magrittr)

if (interactive()) {
  ## Import the example Motive data included in the package
  motive_data <-
    read_motive_csv(system.file("extdata",
                                "pathviewr_motive_example_data.csv",
                                package = 'pathviewr'))

  ## Clean, isolate, and label trajectories
  motive_full <-
    motive_data \%>\%
    clean_viewr(desired_percent = 50,
                max_frame_gap = "autodetect",
                span = 0.95)

  ## Plot all trajectories by subject
  motive_full \%>\%
    plot_by_subject()

  ## Add treatment information
  motive_full$treatment <- c(rep("latA", 100), rep("latB", 100),
                             rep("latA", 100), rep("latB", 149))

  ## Plot all trajectories by subject, color by treatment
  motive_full \%>\%
    plot_by_subject(col_by_treat = TRUE)
}
}
\seealso{
Other plotting functions: 
\code{\link{plot_viewr_trajectories}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Melissa S. Armstrong
}
\concept{plotting functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{get_sf}
\alias{get_sf}
\title{Estimate the spatial frequency of visual stimuli from the subject's
perspective in an experimental tunnel.}
\usage{
get_sf(obj_name)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"} and
\code{vis_angles_calculated}.}
}
\value{
A tibble or data.frame with added variables for
\code{sf_pos}, \code{sf_neg}, and \code{sf_end}.
angle.
}
\description{
Estimate the spatial frequency of visual stimuli from the subject's
perspective in an experimental tunnel.
}
\details{
\code{get_sf()} assumes the following:
\itemize{
\item The subject's gaze is fixed at the point on the either side of the tunnel
that minimizes the distance to visual stimuli and therefore maximizes visual
angles.
\item The subject's head is facing parallel to the length axis of the tunnel.
Visual perception functions in future versions of pathviewr will integrate
head orientation coordinates.
Spatial frequency is reported in cycles/degree and is the inverse of visual
angle (degrees/cycle).
}
}
\examples{
 ## Import sample data from package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))
flydra_data <-
  read_flydra_mat(system.file("extdata", "pathviewr_flydra_example_data.mat",
                              package = 'pathviewr'),
                              subject_name = "birdie_sanders")

 ## Process data up to and including get_vis_angle()
motive_data_full <-
  motive_data \%>\%
  relabel_viewr_axes() \%>\%
  gather_tunnel_data() \%>\%
  trim_tunnel_outliers() \%>\%
  rotate_tunnel() \%>\%
  select_x_percent(desired_percent = 50) \%>\%
  separate_trajectories(max_frame_gap = "autodetect") \%>\%
  get_full_trajectories(span = 0.95) \%>\%
  insert_treatments(tunnel_config = "v",
                   perch_2_vertex = 0.4,
                   vertex_angle = 90,
                   tunnel_length = 2,
                   stim_param_lat_pos = 0.1,
                   stim_param_lat_neg = 0.1,
                   stim_param_end_pos = 0.3,
                   stim_param_end_neg = 0.3,
                   treatment = "lat10_end_30") \%>\%
  calc_min_dist_v(simplify_output = TRUE) \%>\%
  get_vis_angle() \%>\%

  ## Now calculate the spatial frequencies
  get_sf()

  flydra_data_full <-
  flydra_data \%>\%
  redefine_tunnel_center(length_method = "middle",
                        height_method = "user-defined",
                        height_zero = 1.44) \%>\%
  select_x_percent(desired_percent = 50) \%>\%
  separate_trajectories(max_frame_gap = "autodetect") \%>\%
  get_full_trajectories(span = 0.95) \%>\%
  insert_treatments(tunnel_config = "box",
                   tunnel_length = 3,
                   tunnel_width = 1,
                   stim_param_lat_pos = 0.1,
                   stim_param_lat_neg = 0.1,
                   stim_param_end_pos = 0.3,
                   stim_param_end_neg = 0.3,
                   treatment = "lat10_end_30") \%>\%
  calc_min_dist_box() \%>\%
  get_vis_angle() \%>\%

  ## Now calculate the spatial frequencies
  get_sf()
}
\seealso{
Other visual perception functions: 
\code{\link{calc_min_dist_box}()},
\code{\link{get_vis_angle}()}
}
\author{
Eric R. Press
}
\concept{visual perception functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{get_3d_angle}
\alias{get_3d_angle}
\title{Compute an angle in 3D space}
\usage{
get_3d_angle(x1, y1, z1, x2, y2, z2, x3, y3, z3)
}
\arguments{
\item{x1}{x-coordinate of first point}

\item{y1}{y-coordinate of first point}

\item{z1}{z-coordinate of first point}

\item{x2}{x-coordinate of second point (vertex)}

\item{y2}{y-coordinate of second point (vertex)}

\item{z2}{y-coordinate of second point (vertex)}

\item{x3}{x-coordinate of third point}

\item{y3}{y-coordinate of third point}

\item{z3}{z-coordinate of third point}
}
\value{
A numeric vector that provides the angular measurement in degrees.
}
\description{
Compute an angle in 3D space
}
\details{
Everything supplied to arguments must be singular numeric values.
The second point (x2, y2, z2) is treated as the vertex, and the angle between
the three points in 3D space is computed.
}
\examples{
get_3d_angle(
  0, 1, 0,
  0, 0, 0,
  1, 0, 0)
}
\seealso{
Other mathematical functions: 
\code{\link{calc_min_dist_v}()},
\code{\link{deg_2_rad}()},
\code{\link{find_curve_elbow}()},
\code{\link{get_2d_angle}()},
\code{\link{get_3d_cross_prod}()},
\code{\link{get_dist_point_line}()},
\code{\link{get_traj_velocities}()},
\code{\link{get_velocity}()},
\code{\link{rad_2_deg}()}
}
\author{
Vikram B. Baliga
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_in_one_functions.R
\name{import_and_clean_viewr}
\alias{import_and_clean_viewr}
\title{Import + clean_viewr()}
\usage{
import_and_clean_viewr(
  file_name,
  file_id = NA,
  relabel_viewr_axes = TRUE,
  gather_tunnel_data = TRUE,
  trim_tunnel_outliers = TRUE,
  standardization_option = "rotate_tunnel",
  get_velocity = TRUE,
  select_x_percent = TRUE,
  rename_viewr_characters = FALSE,
  separate_trajectories = TRUE,
  get_full_trajectories = TRUE,
  fill_traj_gaps = FALSE,
  ...
)
}
\arguments{
\item{file_name}{Target file}

\item{file_id}{Optional}

\item{relabel_viewr_axes}{default TRUE, should axes be relabeled?}

\item{gather_tunnel_data}{default TRUE, should tunnel data be gathered?}

\item{trim_tunnel_outliers}{default TRUE, outliers be trimmed?}

\item{standardization_option}{default "rotate_tunnel"; which standardization
option should be used? See Details for more.}

\item{get_velocity}{default TRUE, should velocity be computed?}

\item{select_x_percent}{default TRUE, should a region of interest be
selected?}

\item{rename_viewr_characters}{default FALSE, should subjects be renamed?}

\item{separate_trajectories}{default TRUE, should trajectories be defined?}

\item{get_full_trajectories}{default TRUE, should only full trajectories be
retained?}

\item{fill_traj_gaps}{default FALSE, should gaps in trajectories be filled?}

\item{...}{Additional arguments passed to the corresponding functions.}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) that has passed
through several \code{pathviewr} functions as desired by the user,
resulting in data that have been cleaned and ready for analyses.
}
\description{
Import a file and then, akin to \code{clean_viewr}, run through as many
cleaning steps as desired.
}
\details{
Each argument corresponds to a standalone function in
\code{pathviewr}. E.g. the parameter \code{relabel_viewr_axes} allows a
user to choose whether \code{pathviewr::relabel_viewr_axes()} is run
internally. Should the user desire to use any non-default parameter values
for any functions included here, they should be supplied to this function
as additional arguments formatted exactly as they would appear in their
corresponding function(s). E.g. if the "autodetect" feature in
\code{pathviewr::separate_trajectories()} is desired, add an argument
\code{max_frame_gap = "autodetect"} to the arguments supplied to this
function.
}
\examples{
library(pathviewr)

## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

motive_full <-
  motive_data \%>\%
  clean_viewr(desired_percent = 50,
              max_frame_gap = "autodetect",
              span = 0.95)

## Alternatively, used the import_and_clean_viewr()
## function to combine these steps
motive_import_and_clean <-
  import_and_clean_viewr(
    file_name = system.file("extdata", "pathviewr_motive_example_data.csv",
                            package = 'pathviewr'),
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )
}
\seealso{
Other all in one functions: 
\code{\link{clean_viewr}()}
}
\author{
Vikram B. Baliga
}
\concept{all in one functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{get_vis_angle}
\alias{get_vis_angle}
\title{Estimate visual angles from a subject's perspective in an experimental tunnel}
\usage{
get_vis_angle(obj_name)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with
attributes \code{pathviewr_steps} that include \code{"viewr"}
and \code{min_dist_calculated}.}
}
\value{
A tibble or data.frame with added variables for
\code{vis_angle_pos_rad}, \code{vis_angle_pos_deg},
\code{vis_angle_neg_rad}, \code{vos_angle_neg_deg},
\code{vis_angle_end_rad}, and \code{vis_angle_end_deg}.
}
\description{
Estimate visual angles from a subject's perspective in an experimental tunnel
}
\details{
\code{get_vis_angle()} assumes the following:
\itemize{
\item The subject's gaze is fixed at the point on the either side of the tunnel
that minimizes the distance to visual stimuli and therefore maximizes visual
angles.
\item The subject's head is facing parallel to the length axis of the tunnel.
Visual perception functions in future versions of pathviewr will integrate
head orientation coordinates.
Angles are reported in radians/cycle (\code{vis_angle_pos_rad}) and
degrees/cycle (\code{vis_angle_pos_deg}).
}
}
\examples{
 ## Import sample data from package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))
flydra_data <-
  read_flydra_mat(system.file("extdata", "pathviewr_flydra_example_data.mat",
                              package = 'pathviewr'),
                              subject_name = "birdie_sanders")

 ## Process data up to and including get_min_dist()
motive_data_full <-
  motive_data \%>\%
  relabel_viewr_axes() \%>\%
  gather_tunnel_data() \%>\%
  trim_tunnel_outliers() \%>\%
  rotate_tunnel() \%>\%
  select_x_percent(desired_percent = 50) \%>\%
  separate_trajectories(max_frame_gap = "autodetect") \%>\%
  get_full_trajectories(span = 0.95) \%>\%
  insert_treatments(tunnel_config = "v",
                   perch_2_vertex = 0.4,
                   vertex_angle = 90,
                   tunnel_length = 2,
                   stim_param_lat_pos = 0.1,
                   stim_param_lat_neg = 0.1,
                   stim_param_end_pos = 0.3,
                   stim_param_end_neg = 0.3,
                   treatment = "lat10_end_30") \%>\%
  calc_min_dist_v(simplify_output = TRUE) \%>\%

  ## Now calculate the visual angles
  get_vis_angle()

  flydra_data_full <-
  flydra_data \%>\%
  redefine_tunnel_center(length_method = "middle",
                        height_method = "user-defined",
                        height_zero = 1.44) \%>\%
  select_x_percent(desired_percent = 50) \%>\%
  separate_trajectories(max_frame_gap = "autodetect") \%>\%
  get_full_trajectories(span = 0.95) \%>\%
  insert_treatments(tunnel_config = "box",
                   tunnel_length = 3,
                   tunnel_width = 1,
                   stim_param_lat_pos = 0.1,
                   stim_param_lat_neg = 0.1,
                   stim_param_end_pos = 0.3,
                   stim_param_end_neg = 0.3,
                   treatment = "lat10_end_30") \%>\%
  calc_min_dist_box() \%>\%

   ## Now calculate the visual angles
  get_vis_angle()
}
\seealso{
Other visual perception functions: 
\code{\link{calc_min_dist_box}()},
\code{\link{get_sf}()}
}
\author{
Eric R. Press
}
\concept{visual perception functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{get_traj_velocities}
\alias{get_traj_velocities}
\title{Recompute trajectory-specific velocities}
\usage{
get_traj_velocities(
  obj_name,
  time_col = "time_sec",
  length_col = "position_length",
  width_col = "position_width",
  height_col = "position_height",
  set_init_vel_zero = FALSE,
  velocity_min = NA,
  velocity_max = NA
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{time_col}{Name of the column containing time}

\item{length_col}{Name of the column containing length dimension}

\item{width_col}{Name of the column containing width dimension}

\item{height_col}{Name of the column containing height dimension}

\item{set_init_vel_zero}{Should the first value be zero or can it be a
duplicate of the second velocity value? Defaults to FALSE.}

\item{velocity_min}{Should data below a certain velocity be filtered out of
the object? If so, enter a numeric. If not, keep NA.}

\item{velocity_max}{Should data above a certain velocity be filtered out of
the object? If so, enter a numeric. If not, keep NA.}
}
\value{
If \code{add_to_viewr} is \code{TRUE}, additional columns are
appended to the input viewr object. If \code{FALSE}, a standalone tibble is
created. Either way, an "instantaneous" velocity is computed as the
difference in position divided by the difference in time as each successive
row is encountered. Additionally, velocities along each of the three
position axes are computed and provided as additional columns.
}
\description{
Recompute trajectory-specific velocities
}
\details{
Instantaneous velocity is not truly "instantaneous" but rather is
approximated as the change in distance divided by change in time from one
observation (row) to the previous observation (row). Each component of
velocity is computed (i.e. per axis) along with the overall velocity of
the subject.
}
\seealso{
Other mathematical functions: 
\code{\link{calc_min_dist_v}()},
\code{\link{deg_2_rad}()},
\code{\link{find_curve_elbow}()},
\code{\link{get_2d_angle}()},
\code{\link{get_3d_angle}()},
\code{\link{get_3d_cross_prod}()},
\code{\link{get_dist_point_line}()},
\code{\link{get_velocity}()},
\code{\link{rad_2_deg}()}
}
\author{
Vikram B. Baliga
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{remove_vel_anomalies}
\alias{remove_vel_anomalies}
\title{Remove any rows which show sharp shifts in velocity that are likely due to
tracking errors}
\usage{
remove_vel_anomalies(
  obj_name,
  target = "velocity",
  method = "gesd",
  alpha = 0.05,
  max_anoms = 0.2
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{target}{The column to target; defaults to "velocity"}

\item{method}{The anomaly detection method; see anomalize::anomalize()}

\item{alpha}{The width of the "normal" range; see anomalize::anomalize()}

\item{max_anoms}{The max proportion of anomalies; see anomalize::anomalize()}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps}. Rows in which large anomalies were detected have
been removed. No additional columns are created.
}
\description{
Remove any rows which show sharp shifts in velocity that are likely due to
tracking errors
}
\details{
This function runs anomalize::anomalize() on a per-trajectory basis.
The separate_trajectories() and get_full_trajectories() must be run prior
to use.
}
\seealso{
Other utility functions: 
\code{\link{clean_by_span}()},
\code{\link{insert_treatments}()},
\code{\link{remove_duplicate_frames}()},
\code{\link{set_traj_frametime}()}
}
\author{
Vikram B. Baliga
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{rename_viewr_characters}
\alias{rename_viewr_characters}
\title{Rename subjects in the data via pattern detection}
\usage{
rename_viewr_characters(
  obj_name,
  target_column = "subject",
  pattern,
  replacement = ""
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{target_column}{The target column; defaults to "subject"}

\item{pattern}{The (regex) pattern to be replaced}

\item{replacement}{The replacement text. Must be a character}
}
\value{
A tibble or data frame in which subjects have been renamed according
to the \code{pattern} and \code{replacement} supplied by the user.
}
\description{
Quick utility function to use str_replace with mutate(across()) to batch-
rename subjects via pattern detection.
}
\examples{
## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

## Clean the file. It is generally recommended to clean up to the
## "gather" step before running rescale_tunnel_data().
 motive_gathered <-
   motive_data \%>\%
   relabel_viewr_axes() \%>\%
   gather_tunnel_data()

## See the subject names
 unique(motive_gathered$subject)

## Now rename the subjects. We'll get rid of "device" and replace it
## with "subject"
motive_renamed <-
  motive_gathered \%>\%
  rename_viewr_characters(target_column = "subject",
                          pattern = "device",
                          replacement = "subject")

## See the new subject names
unique(motive_renamed$subject)
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_functions.R
\name{bind_viewr_objects}
\alias{bind_viewr_objects}
\title{Bind viewr objects}
\usage{
bind_viewr_objects(obj_list)
}
\arguments{
\item{obj_list}{A list of viewr objects}
}
\value{
A single viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) that combines all the
rows of the source viewr objects in \code{obj_list}. Metadata may not
necessarily be retained and therefore \code{attributes} should be used with
caution.
}
\description{
Combine a list of multiple viewr objects into a single viewr object
}
\examples{
## Since we only have one example file of each type provided
## in pathviewr, we will simply import the same example multiple
## times to simulate batch importing. Replace the contents of
## the following list with your own list of files to be imported.

## Make a list of the same example file 3x
import_list <-
  c(rep(
    system.file("extdata", "pathviewr_motive_example_data.csv",
                package = 'pathviewr'),
    3
  ))

## Batch import
motive_batch_imports <-
  import_batch(import_list,
               import_method = "motive",
               import_messaging = TRUE)

## Batch cleaning of these imported files
## via clean_viewr_batch()
motive_batch_cleaned <-
  clean_viewr_batch(
    file_announce = TRUE,
    motive_batch_imports,
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )

## Alternatively, use import_and_clean_batch() to
## combine import with cleaning on a batch of files
motive_batch_import_and_clean <-
  import_and_clean_batch(
    import_list,
    import_method = "motive",
    import_messaging = TRUE,
    motive_batch_imports,
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )

## Each of these lists of objects can be bound into
## one viewr object (i.e. one tibble) via
## bind_viewr_objects()
motive_bound_one <-
  bind_viewr_objects(motive_batch_cleaned)

motive_bound_two <-
  bind_viewr_objects(motive_batch_import_and_clean)

## Either route results in the same object ultimately:
identical(motive_bound_one, motive_bound_two)
}
\seealso{
Other batch functions: 
\code{\link{clean_viewr_batch}()},
\code{\link{import_and_clean_batch}()},
\code{\link{import_batch}()}
}
\author{
Vikram B. Baliga
}
\concept{batch functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{relabel_viewr_axes}
\alias{relabel_viewr_axes}
\title{Relabel the dimensions as length, width, and height}
\usage{
relabel_viewr_axes(
  obj_name,
  tunnel_length = "_z",
  tunnel_width = "_x",
  tunnel_height = "_y",
  real = "_w",
  ...
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{tunnel_length}{The dimension that corresponds to tunnel length. Set to
\code{tunnel_length = "_z"} by default. Argument should contain a character
vector with a leading underscore (see Details)}

\item{tunnel_width}{The dimension that corresponds to tunnel width. Follows
the same conventions as \code{tunnel_length} and defaults to
\code{tunnel_length = "_x"}}

\item{tunnel_height}{The dimension that corresponds to tunnel height. Follows
the same conventions as \code{tunnel_length} and defaults to
\code{tunnel_length = "_y"}}

\item{real}{The dimension that corresponds to the "real" parameter in
quaternion notation (for data with "rotation" values). Follows the same
conventions as \code{tunnel_length} and defaults to \code{real = "_w"}}

\item{...}{Additional arguments to be passed to \code{read_motive_csv()}.}
}
\value{
A tibble or data.frame with variables that have been renamed.
}
\description{
Axes are commonly labeled as "x", "y", and "z" in recording software yet
\code{pathviewr} functions require these to be labeled as "length", "width",
and "height". \code{relabel_viewr_axes()} is a function that takes a
\code{viewr} object and allows the user to rename its variables.
}
\details{
Each argument must have a leading underscore ("_") and each
argument must have an entry. E.g. tunnel_length = "_Y" will replace all
instances of _Y with _length in the names of variables.
}
\examples{

library(pathviewr)

## Import the Motive example data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

## Names of variables are labeled with _x, _y, _z, which we'd like to rename
names(motive_data)

## Now use relabel_viewr_axes() to rename these variables using _length,
## _width, and _height instead
motive_data_relabeled <-
  relabel_viewr_axes(motive_data,
                     tunnel_length = "_z",
                     tunnel_width = "_x",
                     tunnel_height = "_y",
                     real = "_w")

## See the result
names(motive_data_relabeled)
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_viewr_trajectories}
\alias{plot_viewr_trajectories}
\title{Plot each trajectory within a viewr object}
\usage{
plot_viewr_trajectories(
  obj_name,
  plot_axes = c("length", "width"),
  multi_plot = FALSE
)
}
\arguments{
\item{obj_name}{A viewr object (a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) that has been passed
through \code{separate_trajectories()} or \code{get_full_trajectories()}.}

\item{plot_axes}{Which position axes should be plotted? A character vector
including exactly two of the following choices must be supplied:
\code{length}, \code{width}, \code{height}. Default is c("length",
"width").}

\item{multi_plot}{Should separate plots (one per trajectory) be created or
should one multi-plot grid be generated. Defaults to FALSE, which produces
separate plots.}
}
\value{
A (base-R) series of plots or single plot (if \code{multi_plot =
  TRUE}) that depict each trajectory along the chosen axes.
}
\description{
Plot each trajectory within a viewr object
}
\examples{
library(pathviewr)

## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

motive_full <-
  motive_data \%>\%
  clean_viewr(desired_percent = 50,
              max_frame_gap = "autodetect",
              span = 0.95)

plot_viewr_trajectories(motive_full, multi_plot = FALSE)
plot_viewr_trajectories(motive_full, multi_plot = TRUE)
}
\seealso{
Other plotting functions: 
\code{\link{plot_by_subject}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Vikram B. Baliga
}
\concept{plotting functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{get_2d_angle}
\alias{get_2d_angle}
\title{Compute an angle in 2D space}
\usage{
get_2d_angle(x1, y1, x2, y2, x3, y3)
}
\arguments{
\item{x1}{x-coordinate of first point}

\item{y1}{y-coordinate of first point}

\item{x2}{x-coordinate of second point (vertex)}

\item{y2}{y-coordinate of second point (vertex)}

\item{x3}{x-coordinate of third point}

\item{y3}{y-coordinate of third point}
}
\value{
A numeric vector that provides the angular measurement in degrees.
}
\description{
Compute an angle in 2D space
}
\details{
Everything supplied to arguments must be singular numeric values.
The second point (x2, y2) is treated as the vertex, and the angle between
the three points in 2D space is computed.
}
\examples{
get_2d_angle(
  0, 1,
  0, 0,
  1, 0)
}
\seealso{
Other mathematical functions: 
\code{\link{calc_min_dist_v}()},
\code{\link{deg_2_rad}()},
\code{\link{find_curve_elbow}()},
\code{\link{get_3d_angle}()},
\code{\link{get_3d_cross_prod}()},
\code{\link{get_dist_point_line}()},
\code{\link{get_traj_velocities}()},
\code{\link{get_velocity}()},
\code{\link{rad_2_deg}()}
}
\author{
Vikram B. Baliga
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{calc_min_dist_box}
\alias{calc_min_dist_box}
\title{Calculate minimum distance to lateral and end walls in a box-shaped
experimental tunnel}
\usage{
calc_min_dist_box(obj_name)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that include \code{"viewr"} and
\code{treatments_added}.}
}
\value{
A tibble or data.frame with added variables for
\code{min_dist_pos}, \code{min_dist_neg}, and \code{min_dist_end},.
}
\description{
Calculate minimum distance to lateral and end walls in a box-shaped
experimental tunnel
}
\details{
\code{calc_min_dist_box()} assumes the subject locomotes facing
forward, therefore \code{min_dist_end} represents the minimum distance
between the subject and the end wall to which it is moving towards.
All outputs are in meters.
}
\examples{
## Import sample data from package
 flydra_data <-
 read_flydra_mat(system.file("extdata", "pathviewr_flydra_example_data.mat",
                               package = 'pathviewr'),
                               subject_name = "birdie_sanders")

   ## Process data up to and including insert_treatments()
  flydra_data_full <-
   flydra_data \%>\%
   redefine_tunnel_center(length_method = "middle",
                         height_method = "user-defined",
                         height_zero = 1.44) \%>\%
   select_x_percent(desired_percent = 50) \%>\%
   separate_trajectories(max_frame_gap = "autodetect") \%>\%
   get_full_trajectories(span = 0.95) \%>\%
   insert_treatments(tunnel_config = "box",
                    tunnel_length = 3,
                    tunnel_width = 1,
                    stim_param_lat_pos = 0.1,
                    stim_param_lat_neg = 0.1,
                    stim_param_end_pos = 0.3,
                    stim_param_end_neg = 0.3,
                    treatment = "lat10_end_30") \%>\%

   ## Now calculate the minimum distances to each wall
   calc_min_dist_box()

   ## See 3 new variables for calculations to lateral and end walls
   names(flydra_data_full)
}
\seealso{
Other visual perception functions: 
\code{\link{get_sf}()},
\code{\link{get_vis_angle}()}
}
\author{
Eric R. Press
}
\concept{visual perception functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{rm_by_trajnum}
\alias{rm_by_trajnum}
\title{Remove subjects by trajectory number}
\usage{
rm_by_trajnum(
  obj_name,
  trajnum = 5,
  mirrored = FALSE,
  treatment1,
  treatment2,
  ...
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}. Trajectories must be
predefined (i.e. via \code{separate_trajectories()}).}

\item{trajnum}{Minimum number of trajectories; must be numeric.}

\item{mirrored}{Does the data have mirrored treatments? If so, arguments
\code{treatment1} and \code{treatment2} must also be provided, indicating
the names of two mirrored treatments, both of which must meet the trajectory
threshold specified in \code{trajnum}. Default is FALSE.}

\item{treatment1}{The first treatment or session during which the threshold
must be met.}

\item{treatment2}{A second treatment or session during which the threshold
must be met.}

\item{...}{Additional arguments passed to/from other pathviewr functions.}
}
\value{
A viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"} that now has fewer
observations (rows) as a result of removal of subjects with too few
trajectories according to the \code{trajnum} parameter.
}
\description{
Specify a minimum number of trajectories that each subject must complete
during a treatment, trial, or session.
}
\details{
Depending on analysis needs, users may want to remove subjects that
have not completed a certain number of trajectories during a treatment,
trial, or session. If \code{mirrored = FALSE}, no treatment information is
necessary and subjects will be removed based on total number of trajectories
as specified in \code{trajnum}. If \code{mirrored = TRUE}, the
\code{treatment1} and \code{treatment2} parameters will allow users to
define during which treatments or sessions subjects must reach threshold as
specified in the \code{trajnum} argument. For example, if \code{mirrored =
 TRUE}, setting \code{treatment1 = "latA"}, \code{treatment2 = "latB"} and
\code{trajnum = 5} will remove subjects that have fewer than 5 trajectories
during the \code{"latA"} treatment AND the \code{"latB"} treatment.
\code{treatment1} and \code{treatment2} should be levels within a column
named \code{"treatment"}.
}
\examples{
library(pathviewr)

## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))

## Clean, isolate, and label trajectories
motive_full <-
  motive_data \%>\%
  clean_viewr(desired_percent = 50,
              max_frame_gap = "autodetect",
              span = 0.95)

##Remove subjects that have not completed at least 150 trajectories:
motive_rm_unmirrored <-
  motive_full \%>\%
  rm_by_trajnum(trajnum = 150)

## Add treatment information
motive_full$treatment <- c(rep("latA", 100),
                           rep("latB", 100),
                           rep("latA", 100),
                           rep("latB", 149))

## Remove subjects by that have not completed at least 10 trajectories in
## both treatments
motive_rm_mirrored <-
  motive_full \%>\%
  rm_by_trajnum(
    trajnum = 10,
    mirrored = TRUE,
    treatment1 = "latA",
    treatment2 = "latB"
  )
}
\author{
Melissa S. Armstrong
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{get_full_trajectories}
\alias{get_full_trajectories}
\title{Retain trajectories that span a selected region of interest}
\usage{
get_full_trajectories(obj_name, span = 0.8, ...)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{span}{Span to use; must be numeric and between 0 and 1}

\item{...}{Additional arguments passed to/from other pathviewr functions}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) in which only
trajectories that span the region of interest are retained. Data are
labeled by direction  (either "leftwards" or "rightwards") with respect to
their starting and ending \code{position_length} values in the
\code{direction} column.
}
\description{
Specify a minimum span of the selected region of interest and then keep
trajectories that are wider than that span and go from one end to the other
of the region.
}
\details{
Because trajectories may not have observations exactly at the
beginning or the end of the region of interest, it may be necessary to allow
trajectories to be slightly shorter than the range of the selected region of
interest. The \code{span} parameter of this function handles this. By
supplying a numeric proportion from 0 to 1, a user may allow trajectories to
span that proportion of the selected region. For example, setting \code{span
= 0.95} will keep all trajectories that span 95\% of the length of the
selected region of interest. Setting \code{span = 1} (not recommended) will
strictly keep trajectories that start and end at the exact cut-offs of the
selected region of interest. For these reasons, \code{span}s of 0.99 to 0.95
are generally recommended.
}
\examples{
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))

## Clean the file. It is generally recommended to clean up to the
## "separate" step before running select_x_percent().
motive_separated <-
  motive_data \%>\%
  relabel_viewr_axes() \%>\%
  gather_tunnel_data() \%>\%
  trim_tunnel_outliers() \%>\%
  rotate_tunnel() \%>\%
  select_x_percent(desired_percent = 50) \%>\%
  separate_trajectories(max_frame_gap = "autodetect",
                        frame_rate_proportion = 0.1)

## Now retain only the "full" trajectories that span
## across 0.95 of the range of position_length
motive_full <-
  motive_separated \%>\%
  get_full_trajectories(span = 0.95)
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}

Other functions that define or clean trajectories: 
\code{\link{quick_separate_trajectories}()},
\code{\link{separate_trajectories}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
\concept{functions that define or clean trajectories}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{set_traj_frametime}
\alias{set_traj_frametime}
\title{Redefine frames and time stamps on a per-trajectory basis}
\usage{
set_traj_frametime(obj_name)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps}. New columns include traj_time (the
trajectory-specific time values) and traj_frame (the trajectory-specific
frame numbering).
}
\description{
After a data set has been separated into trajectories, find the earliest
frame in each trajectory and set the corresponding time to 0. All subsequent
time_sec stamps are computed according to successive frame numbering.
}
\details{
The separate_trajectories() and get_full_trajectories() must be
run prior to use. The initial traj_time and traj_frame values are set to 0
within each trajectory.
}
\seealso{
Other utility functions: 
\code{\link{clean_by_span}()},
\code{\link{insert_treatments}()},
\code{\link{remove_duplicate_frames}()},
\code{\link{remove_vel_anomalies}()}
}
\author{
Vikram B. Baliga
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{calc_min_dist_v}
\alias{calc_min_dist_v}
\title{Calculate minimum distance to lateral and end walls in a V-shaped
experimental tunnel}
\usage{
calc_min_dist_v(obj_name, simplify_output = TRUE)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"} and
\code{treatments_added}.}

\item{simplify_output}{If TRUE, the returned object includes only the minimum
distance between the subject and the lateral/end walls. If FALSE, the
returned object includes all variables internal to the calculation.}
}
\value{
A tibble or data.frame with added variables for
\code{height_2_vertex}, \code{height_2_screen}, \code{width_2_screen_pos},
\code{width_2_screen_neg}, \code{min_dist_pos}, \code{min_dist_neg},
\code{min_dist_end}, \code{bound_pos}, and \code{bound_neg}.
}
\description{
Calculate minimum distance to lateral and end walls in a V-shaped
experimental tunnel
}
\details{
For tunnels in which \code{vertex_angle} is >90 degree, \code{bound_pos}
and \code{bound_neg} represent a planes orthogonal to the lateral walls and
are used to modify \code{min_dist_pos} and \code{min_dist_neg} calculations
to prevent erroneous outputs.
\code{calc_min_dist_v()} assumes the subject locomotes facing forward,
therefore \code{min_dist_end} represents the minimum distance between the
subject and the end wall to which it is moving towards
All outputs are in meters.
}
\examples{
 ## Import sample data from package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))

 ## Process data up to and including insert_treatments()
motive_data_full <-
  motive_data \%>\%
  relabel_viewr_axes() \%>\%
  gather_tunnel_data() \%>\%
  trim_tunnel_outliers() \%>\%
  rotate_tunnel() \%>\%
  select_x_percent(desired_percent = 50) \%>\%
  separate_trajectories(max_frame_gap = "autodetect") \%>\%
  get_full_trajectories(span = 0.95) \%>\%
  insert_treatments(tunnel_config = "v",
                   perch_2_vertex = 0.4,
                   vertex_angle = 90,
                   tunnel_length = 2,
                   stim_param_lat_pos = 0.1,
                   stim_param_lat_neg = 0.1,
                   stim_param_end_pos = 0.3,
                   stim_param_end_neg = 0.3,
                   treatment = "lat10_end_30") \%>\%

 ## Now calculate the minimum distances to each wall
  calc_min_dist_v(simplify_output = TRUE)

  ## See 3 new variables for calculations to lateral and end walls
  names(motive_data_full)
}
\seealso{
Other mathematical functions: 
\code{\link{deg_2_rad}()},
\code{\link{find_curve_elbow}()},
\code{\link{get_2d_angle}()},
\code{\link{get_3d_angle}()},
\code{\link{get_3d_cross_prod}()},
\code{\link{get_dist_point_line}()},
\code{\link{get_traj_velocities}()},
\code{\link{get_velocity}()},
\code{\link{rad_2_deg}()}
}
\author{
Eric R. Press
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{trim_tunnel_outliers}
\alias{trim_tunnel_outliers}
\title{Trim out artifacts and other outliers from the extremes of the tunnel}
\usage{
trim_tunnel_outliers(
  obj_name,
  lengths_min = 0,
  lengths_max = 3,
  widths_min = -0.4,
  widths_max = 0.8,
  heights_min = -0.2,
  heights_max = 0.5,
  ...
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"} that has been passed
through \code{relabel_viewr_axes()} and \code{gather_tunnel_data()} (or is
structured as though it has been passed through those functions).}

\item{lengths_min}{Minimum length}

\item{lengths_max}{Maximum length}

\item{widths_min}{Minimum width}

\item{widths_max}{Maximum width}

\item{heights_min}{Minimum height}

\item{heights_max}{Maximum height}

\item{...}{Additional arguments passed to/from other pathviewr functions}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) in which data outside
the specified ranges has been excluded.
}
\description{
The user provides estimates of min and max values of data. This function then
trims out anything beyond these estimates.
}
\details{
Anything supplied to _min or _max arguments should be single numeric
values.
}
\examples{
## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))

## Clean the file. It is generally recommended to clean up to the
## "gather" step before running trim_tunnel_outliers().
motive_gathered <-
  motive_data \%>\%
  relabel_viewr_axes() \%>\%
  gather_tunnel_data()

## Now trim outliers using default values
motive_trimmed <-
  motive_gathered \%>\%
  trim_tunnel_outliers(lengths_min = 0,
                       lengths_max = 3,
                       widths_min = -0.4,
                       widths_max = 0.8,
                       heights_min = -0.2,
                       heights_max = 0.5)
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{exclude_by_velocity}
\alias{exclude_by_velocity}
\title{Remove trajectories entirely, based on velocity thresholds}
\usage{
exclude_by_velocity(obj_name, vel_min = NULL, vel_max = NULL)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{vel_min}{Default \code{NULL}. If a numeric is entered, trajectories
that have at least one observation with velocity less than \code{vel_min}
are removed.}

\item{vel_max}{Default \code{NULL}. If a numeric is entered, trajectories
that have at least one observation with velocity greater than
\code{vel_max} are removed.}
}
\value{
A new viewr object that is identical to the input object but now
excludes any trajectories that contain observations with velocity less than
\code{vel_min} (if specified) and/or velocity greater than \code{vel_max}
(if specified)
}
\description{
Remove trajectories from a viewr object that contain instances of velocity
known to be spurious.
}
\examples{
## Import and clean the example Motive data
motive_import_and_clean <-
  import_and_clean_viewr(
    file_name = system.file("extdata", "pathviewr_motive_example_data.csv",
                            package = 'pathviewr'),
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )

## See the distribution of velocities
hist(motive_import_and_clean$velocity)

## Let's remove any trajectories that contain
## velocity < 2
motive_vel_filtered <-
  motive_import_and_clean \%>\%
  exclude_by_velocity(vel_min = 2)

## See how the distribution of velocities has changed
hist(motive_vel_filtered$velocity)
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{clean_by_span}
\alias{clean_by_span}
\title{Remove file_sub_traj entries that do not span the full region of interest}
\usage{
clean_by_span(
  obj_name,
  axis = "position_length",
  min_value = NULL,
  max_value = NULL,
  tolerance = 0.1
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{axis}{Along which axis should restrictions be enforced?}

\item{min_value}{Minimum coordinate value; setting this to NULL will
auto-compute the best value}

\item{max_value}{Maximum coordinate; setting this to NULL will auto-compute
the best value}

\item{tolerance}{As a proporiton of axis value}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps}. Trajectories that do not span the full region of
interest have been removed; trajectory identities (file_sub_traj) have
not been changed.
}
\description{
Remove file_sub_traj entries that do not span the full region of interest
}
\seealso{
Other utility functions: 
\code{\link{insert_treatments}()},
\code{\link{remove_duplicate_frames}()},
\code{\link{remove_vel_anomalies}()},
\code{\link{set_traj_frametime}()}
}
\author{
Vikram B. Baliga
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{deg_2_rad}
\alias{deg_2_rad}
\title{Convert degrees to radians}
\usage{
deg_2_rad(deg)
}
\arguments{
\item{deg}{Degrees (a numeric of any length >= 1)}
}
\value{
The angle(s) in radians (as a numeric vector of the same length)
}
\description{
Convert degrees to radians
}
\examples{
## One input
deg_2_rad(90)

## Multiple inputs
deg_2_rad(c(5, 10, 15, 20))
}
\seealso{
Other mathematical functions: 
\code{\link{calc_min_dist_v}()},
\code{\link{find_curve_elbow}()},
\code{\link{get_2d_angle}()},
\code{\link{get_3d_angle}()},
\code{\link{get_3d_cross_prod}()},
\code{\link{get_dist_point_line}()},
\code{\link{get_traj_velocities}()},
\code{\link{get_velocity}()},
\code{\link{rad_2_deg}()}
}
\author{
Vikram B. Baliga
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{get_header_viewr}
\alias{get_header_viewr}
\title{Extract header info from imported viewr object}
\usage{
get_header_viewr(obj_name, ...)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}
\code{pathviewr_steps}}

\item{...}{Additional arguments that may be passed to other \code{pathviewr}
functions}
}
\value{
The value of the \code{header} attribute, or NULL if no exact match
is found and no or more than one partial match is found.
}
\description{
A function to quickly return the information stored in the header of the
original data file imported via \code{pathviewr}'s \code{read_} functions.
}
\examples{
library(pathviewr)

## Import the Motive example data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

## Now display the Header information
get_header_viewr(motive_data)
}
\author{
Vikram B. Baliga
}
\concept{metadata handling functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_functions.R
\name{clean_viewr_batch}
\alias{clean_viewr_batch}
\title{Batch clean viewr files}
\usage{
clean_viewr_batch(obj_list, file_announce = FALSE, ...)
}
\arguments{
\item{obj_list}{A list of viewr objects (i.e. a list of tibbles that each
have attribute \code{pathviewr_steps} that includes \code{"viewr"})}

\item{file_announce}{Should the function report each time a file is
processed? Default FALSE; if TRUE, a message will appear in the console
each time a file has been cleaned successfully.}

\item{...}{Arguments to be passed in that specify how this function should
clean files.}
}
\value{
A list of viewr objects (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) that have been passed
through the corresponding cleaning functions.
}
\description{
For a list of viewr objects, run through the pipeline (from relabel axes
up through get full trajectories, as desired) via clean_viewr()
}
\details{
viewr objects should be in a list, e.g. the object generated by
\code{import_batch()}.

See \code{clean_viewr()} for details of how cleaning steps are handled
and/or refer to the corresponding cleaning functions themselves.
}
\examples{
## Since we only have one example file of each type provided
## in pathviewr, we will simply import the same example multiple
## times to simulate batch importing. Replace the contents of
## the following list with your own list of files to be imported.

## Make a list of the same example file 3x
import_list <-
  c(rep(
    system.file("extdata", "pathviewr_motive_example_data.csv",
                package = 'pathviewr'),
    3
  ))

## Batch import
motive_batch_imports <-
  import_batch(import_list,
               import_method = "motive",
               import_messaging = TRUE)

## Batch cleaning of these imported files
## via clean_viewr_batch()
motive_batch_cleaned <-
  clean_viewr_batch(
    file_announce = TRUE,
    motive_batch_imports,
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )

## Alternatively, use import_and_clean_batch() to
## combine import with cleaning on a batch of files
motive_batch_import_and_clean <-
  import_and_clean_batch(
    import_list,
    import_method = "motive",
    import_messaging = TRUE,
    motive_batch_imports,
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )

## Each of these lists of objects can be bound into
## one viewr object (i.e. one tibble) via
## bind_viewr_objects()
motive_bound_one <-
  bind_viewr_objects(motive_batch_cleaned)

motive_bound_two <-
  bind_viewr_objects(motive_batch_import_and_clean)

## Either route results in the same object ultimately:
identical(motive_bound_one, motive_bound_two)
}
\seealso{
Other batch functions: 
\code{\link{bind_viewr_objects}()},
\code{\link{import_and_clean_batch}()},
\code{\link{import_batch}()}
}
\author{
Vikram B. Baliga
}
\concept{batch functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\value{
No return value, called for side effects
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{visualize_frame_gap_choice}
\alias{visualize_frame_gap_choice}
\title{Visualize the consequence of using various max_frame_gap values}
\usage{
visualize_frame_gap_choice(obj_name, loops = 20, ...)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{loops}{How many total frame gap entries to consider. Each loop will
increase the \code{max_fram_gap} argument in \code{separate_trajectories}
by 1.}

\item{...}{Additional arguments}
}
\value{
A plot and a tibble, each of which shows the total number of
trajectories that result from using the specified range of
\code{max_frame_gap} values.
}
\description{
Run separate_trajectories() with many different frame gaps to help determine
what value to use
}
\details{
The input viewr object (\code{obj_name}) should likely be an object
that has passed through the \code{select_x_percent()} step.
}
\examples{
library(pathviewr)

## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

motive_selected <-
  motive_data \%>\%
  relabel_viewr_axes() \%>\%
  gather_tunnel_data() \%>\%
  trim_tunnel_outliers() \%>\%
  rotate_tunnel() \%>\%
  get_velocity() \%>\%
  select_x_percent(desired_percent = 50)

visualize_frame_gap_choice(motive_selected, loops = 10)
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()}

Other plotting functions: 
\code{\link{plot_by_subject}()},
\code{\link{plot_viewr_trajectories}()}

Other functions that define or clean trajectories: 
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{separate_trajectories}()}
}
\author{
Melissa S. Armstrong and Vikram B. Baliga
}
\concept{data cleaning functions}
\concept{functions that define or clean trajectories}
\concept{plotting functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{gather_tunnel_data}
\alias{gather_tunnel_data}
\title{Gather data columns into key-value pairs}
\usage{
gather_tunnel_data(obj_name, NA_drop = TRUE, ...)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{NA_drop}{Should rows with NAs be dropped? Defaults to \code{TRUE}}

\item{...}{Additional arguments that can be passed to other \code{pathviewr}
functions such as \code{relabel_viewr_axes()} or \code{read_motive_csv()}}
}
\value{
A tibble in "tidy" format which is formatted to have every row
correspond to the position (and potentially rotation) of a single subject
during an observed frame and time. Subjects' names are automatically parsed
from original variable names (e.g. subject1_rotation_width extracts
"subject1" as the subject name) and stored in a \code{Subjects} column in the
returned tibble.
}
\description{
Reformat \code{viewr} data into a "tidy" format so that every row corresponds
to the position (and potentially rotation) of a single subject during an
observed frame and time.
}
\details{
The tibble or data.frame that is fed in must have variables that
have subject names and axis names separated by underscores. Axis names must
be one of the following: \code{position_length}, \code{position_width}, or
\code{position_height}. Each of these three dimensions must be present in the
data. Collectively, this means that names like \code{bird01_position_length}
or \code{larry_position_height} are acceptable, but \code{bird01_x} or
\code{bird01_length} are not.
}
\examples{
library(pathviewr)

## Import the Motive example data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

## First use relabel_viewr_axes() to rename these variables using _length,
## _width, and _height instead
motive_data_relabeled <- relabel_viewr_axes(motive_data)

## Now use gather_tunnel_data() to gather colums into tidy format
motive_data_gathered <- gather_tunnel_data(motive_data_relabeled)

## Column names reflect the way in which data were reformatted:
names(motive_data_gathered)
}
\seealso{
Other data cleaning functions: 
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_import_functions.R
\name{as_viewr}
\alias{as_viewr}
\title{Convert data from another format into a viewr object}
\usage{
as_viewr(
  obj_name,
  frame_rate = 100,
  frame_col,
  time_col,
  subject_col,
  position_length_col,
  position_width_col,
  position_height_col,
  include_rotation = FALSE,
  rotation_real_col,
  rotation_length_col,
  rotation_width_col,
  rotation_height_col
)
}
\arguments{
\item{obj_name}{A tibble or data frame containing movement trajectories}

\item{frame_rate}{Must be a single numeric value indicating capture frame
rate in frames per second.}

\item{frame_col}{Column number of obj_name that contains frame numbers}

\item{time_col}{Column number of obj_name that contains time (must be in
seconds)}

\item{subject_col}{Column number of obj_name that contains subject name(s)}

\item{position_length_col}{Column number of obj_name that contains
length-axis position values}

\item{position_width_col}{Column number of obj_name that contains width-axis
position values}

\item{position_height_col}{Column number of obj_name that contains
height-axis position values}

\item{include_rotation}{Are rotation data included? Defaults to FALSE}

\item{rotation_real_col}{Column number of obj_name that contains the "real"
axis of quaternion rotation data}

\item{rotation_length_col}{Column number of obj_name that contains the length
axis of quaternion rotation data}

\item{rotation_width_col}{Column number of obj_name that contains the width
axis of quaternion rotation data}

\item{rotation_height_col}{Column number of obj_name that contains the height
axis of quaternion rotation data}
}
\value{
A tibble that is organized to be compliant with other
\code{pathviewr} functions and that contains the attributes
\code{pathviewr_steps} with entries set to \code{c("viewr",
  "renamed_tunnel", "gathered_tunnel")}
}
\description{
Should you have data from a non-Motive, non-Flydra source, this function can
be used to ensure your data are put into the right format to work with other
pathviewr functions.
}
\examples{

## Create a dummy data frame with simulated (nonsense) data
df <- data.frame(frame = seq(1, 100, by = 1),
                 time_sec = seq(0, by = 0.01, length.out = 100),
                 subject = "birdie_sanders",
                 z = rnorm(100),
                 x = rnorm(100),
                 y = rnorm(100))

## Use as_viewr() to convert it into a viewr object
test <-
  as_viewr(
    df,
    frame_rate = 100,
    frame_col = 1,
    time_col = 2,
    subject_col = 3,
    position_length_col = 5,
    position_width_col = 6,
    position_height_col = 4
  )
}
\seealso{
Other data import functions: 
\code{\link{import_and_clean_batch}()},
\code{\link{import_batch}()},
\code{\link{read_flydra_mat}()},
\code{\link{read_motive_csv}()}
}
\author{
Vikram B. Baliga
}
\concept{data import functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{rescale_tunnel_data}
\alias{rescale_tunnel_data}
\title{Rescale position data within a \code{viewr} object}
\usage{
rescale_tunnel_data(obj_name, original_scale = 0.5, desired_scale = 1, ...)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"} that has been passed
through \code{relabel_viewr_axes()} and \code{gather_tunnel_data()} (or is
structured as though it has been passed through those functions).}

\item{original_scale}{The original scale at which data were exported. See
Details if unknown.}

\item{desired_scale}{The desired scale}

\item{...}{Additional arguments passed to/from other pathviewr functions}
}
\value{
A \code{viewr} object that has position data (and
\code{mean_marker_error data}, if found) adjusted by the ratio of
\code{desired_scale/original_scale}.
}
\description{
Should data have been exported at an incorrect scale, apply an isometric
transformation to the position data and associated mean marker errors (if
found)
}
\details{
The \code{desired_scale} is divided by the \code{original_scale} to
determine a \code{scale_ratio} internally. If the \code{original_scale} is
not explicitly known, set it to 1 and then set \code{desired_scale} to be
whatever scaling ratio you have in mind. E.g. setting \code{original_scale}
to 1 and then \code{desired_scale} to 0.7 will multiply all position axis
values by 0.7.

The default arguments of \code{original_scale = 0.5} and
\code{desired_scale = 1} apply a doubling of tunnel size isometrically.
}
\examples{
## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

## Clean the file. It is generally recommended to clean up to the
## "gather" step before running rescale_tunnel_data().
 motive_gathered <-
   motive_data \%>\%
   relabel_viewr_axes() \%>\%
   gather_tunnel_data()

## Now rescale the tunnel data
motive_rescaled <-
  motive_gathered \%>\%
  rescale_tunnel_data(original_scale = 0.5,
                      desired_scale = 1)

## See the difference in data range e.g. for length
range(motive_rescaled$position_length)
range(motive_gathered$position_length)
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{get_dist_point_line}
\alias{get_dist_point_line}
\title{Compute distance between a point and a line}
\usage{
get_dist_point_line(point, line_coord1, line_coord2)
}
\arguments{
\item{point}{2D or 3D coordinates of the point as c(x,y) or c(x,y,z)}

\item{line_coord1}{2D or 3D coordinates of one point on the line as c(x,y) or
c(x,y,z)}

\item{line_coord2}{2D or 3D coordinates of a second point on the line as
c(x,y) or c(x,y,z)}
}
\value{
A numeric vector of length 1 that provides the euclidean distance
between the point and the line.
}
\description{
Compute distance between a point and a line
}
\details{
The function accepts 2D coordinates or 3D coordinates, but note that
the dimensions of all supplied arguments must match; all coordinates must
be 2D or they all must be 3D.
}
\examples{
## 2D case
get_dist_point_line(
  point = c(0, 0),
  line_coord1 = c(1, 0),
  line_coord2 = c(1, 5)
)

## 3D case
get_dist_point_line(
  point = c(0, 0, 0),
  line_coord1 = c(1, 0, 0),
  line_coord2 = c(1, 5, 0)
)
}
\seealso{
Other mathematical functions: 
\code{\link{calc_min_dist_v}()},
\code{\link{deg_2_rad}()},
\code{\link{find_curve_elbow}()},
\code{\link{get_2d_angle}()},
\code{\link{get_3d_angle}()},
\code{\link{get_3d_cross_prod}()},
\code{\link{get_traj_velocities}()},
\code{\link{get_velocity}()},
\code{\link{rad_2_deg}()}
}
\author{
Vikram B. Baliga
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{insert_treatments}
\alias{insert_treatments}
\title{Inserts treatment and experiment information}
\usage{
insert_treatments(
  obj_name,
  tunnel_config = "box",
  perch_2_vertex = NULL,
  vertex_angle = NULL,
  tunnel_width = NULL,
  tunnel_length = NULL,
  stim_param_lat_pos = NULL,
  stim_param_lat_neg = NULL,
  stim_param_end_pos = NULL,
  stim_param_end_neg = NULL,
  treatment = NULL
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{tunnel_config}{The configuration of the experimental tunnel.
Currently, pathviewr supports rectangular "box" and V-shaped tunnel
configurations.}

\item{perch_2_vertex}{If using a V-shaped tunnel, this is the vertical
distance between the vertex and the height of the perches. If the tunnel does
not have perches, insert the vertical distance between the vertex and the
height of the origin (0,0,0).}

\item{vertex_angle}{If using a V-shaped tunnel, the angle of the vertex (in
degrees) \code{vertex_angle} defaults to 90.}

\item{tunnel_width}{If using a box-shaped tunnel, the width of the tunnel.}

\item{tunnel_length}{The length of the tunnel.}

\item{stim_param_lat_pos}{The size of the stimulus on the lateral positive
wall of the tunnel. Eg. for 10cm wide gratings,
\code{stim_param_lat_pos} = 0.1.}

\item{stim_param_lat_neg}{The size of the stimulus on the lateral negative
wall of the tunnel..}

\item{stim_param_end_pos}{The size of the stimulus on the end positive
wall of the tunnel.}

\item{stim_param_end_neg}{The size of the stimulus on the end negative
wall of the tunnel.}

\item{treatment}{The name of the treatment assigned to all rows of the viewr
object. Currently only able to accept a single treatment per viewr data
object.}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"treatments added"}). Depending
on the argument \code{tunnel_config}, the viewr object also includes
columns storing the values of the supplied arguments. This experimental
information is also stored in the viewr object's metadata
}
\description{
Adds information about treatment and experimental set up to viewr objects for
analysis in other pathviewr functions
}
\details{
All length measurements reported in meters.
}
\examples{
 ## Import sample data from package
 motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))
 flydra_data <-
 read_flydra_mat(system.file("extdata", "pathviewr_flydra_example_data.mat",
                              package = 'pathviewr'),
                              subject_name = "birdie_sanders")

  ## Clean data up to and including get_full_trajectories()
motive_data_full <-
 motive_data \%>\%
 relabel_viewr_axes() \%>\%
 gather_tunnel_data() \%>\%
 trim_tunnel_outliers() \%>\%
 rotate_tunnel() \%>\%
 select_x_percent(desired_percent = 50) \%>\%
 separate_trajectories(max_frame_gap = "autodetect") \%>\%
 get_full_trajectories(span = 0.95)

 flydra_data_full <-
  flydra_data \%>\%
  redefine_tunnel_center(length_method = "middle",
                        height_method = "user-defined",
                        height_zero = 1.44) \%>\%
  select_x_percent(desired_percent = 50) \%>\%
  separate_trajectories(max_frame_gap = "autodetect") \%>\%
  get_full_trajectories(span = 0.95)


## Now add information about the experimental configuration. In this example,
## a V-shaped tunnel in which the vertex is 90deg and lies 0.40m below the
## origin. The visual stimuli on the lateral and end walls have a cycle
## length of 0.1m and 0.3m respectively, and the treatment is labeled
## "lat10_end30"

motive_v <-
motive_data_full \%>\%
 insert_treatments(tunnel_config = "v",
                   perch_2_vertex = 0.4,
                   vertex_angle = 90,
                   tunnel_length = 2,
                   stim_param_lat_pos = 0.1,
                   stim_param_lat_neg = 0.1,
                   stim_param_end_pos = 0.3,
                   stim_param_end_neg = 0.3,
                   treatment = "lat10_end_30")

# For an experiment using the box-shaped configuration where the tunnel is 1m
# wide and 3m long and the visual stimuli on the lateral and end walls have a
# cycle length of 0.2 and 0.3m, respectively, and the treatment is labeled
# "lat20_end30".

flydra_box <-
 flydra_data_full \%>\%
 insert_treatments(tunnel_config = "box",
                   tunnel_width = 1,
                   tunnel_length = 3,
                   stim_param_lat_pos = 0.2,
                   stim_param_lat_neg = 0.2,
                   stim_param_end_pos = 0.3,
                   stim_param_end_neg = 0.3,
                   treatment = "lat20_end30")

## Check out the new columns in the resulting objects
names(motive_v)
names(flydra_box)
}
\seealso{
Other utility functions: 
\code{\link{clean_by_span}()},
\code{\link{remove_duplicate_frames}()},
\code{\link{remove_vel_anomalies}()},
\code{\link{set_traj_frametime}()}
}
\author{
Eric R. Press
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{separate_trajectories}
\alias{separate_trajectories}
\title{Separate rows of data into separately labeled trajectories.}
\usage{
separate_trajectories(
  obj_name,
  max_frame_gap = 1,
  frame_rate_proportion = 0.1,
  frame_gap_messaging = FALSE,
  frame_gap_plotting = FALSE,
  ...
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{max_frame_gap}{Default 1; defines the largest permissible gap in data
before a new trajectory is forced to be defined. Can be either a numeric
value or "autodetect". See Details for more.}

\item{frame_rate_proportion}{Default 0.10; if \code{max_frame_gap =
  "autodetect"}, proportion of frame rate to be used as a reference (see
Details).}

\item{frame_gap_messaging}{Default FALSE; should frame gaps be reported in
the console?}

\item{frame_gap_plotting}{Default FALSE; should frame gap diagnostic plots be
shown?}

\item{...}{Additional arguments passed to/from other pathviewr functions}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) in which a new column
\code{file_sub_traj} is added, which labels trajectories within the data by
concatenating file name, subject name, and a trajectory number (all
separated by underscores).
}
\description{
Separate rows of data into separately labeled trajectories.
}
\details{
This function is designed to separate rows of data into separately
labeled trajectories.

The \code{max_frame_gap} parameter determines how trajectories will be
separated. If numeric, the function uses the supplied value as the largest
permissible gap in frames before a new trajectory is defined.

If \code{max_frame_gap = "autodetect"} is used, the function
attempts to guesstimate the best value(s) of \code{max_frame_gap}. This is
performed separately for each subject in the data set, i.e. as many
\code{max_frame_gap} values will be estimated as there are unique subjects.
The highest possible value of any \code{max_frame_gap} is first set as a
proportion of the capture frame rate, as defined by the
\code{frame_rate_proportion} parameter (default 0.10). For each subject, a
plot of total trajectory counts vs. max frame gap values is created
internally (but can be plotted via setting
\code{frame_gap_plotting = TRUE}). As larger max frame gaps are allowed,
trajectory count will necessarily decrease but may reach a value that
likely represents the "best" option. The "elbow" of that plot is then used
to find an estimate of the best max frame gap value to use.
}
\examples{
## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))

## Clean the file. It is generally recommended to clean up to the
## "select" step before running select_x_percent().
motive_selected <-
  motive_data \%>\%
  relabel_viewr_axes() \%>\%
  gather_tunnel_data() \%>\%
  trim_tunnel_outliers() \%>\%
  rotate_tunnel() \%>\%
  select_x_percent(desired_percent = 50)

## Now separate trajectories using autodetect
motive_separated <-
  motive_selected \%>\%
  separate_trajectories(max_frame_gap = "autodetect",
                        frame_rate_proportion = 0.1)

## See new column file_sub_traj for trajectory labeling
names(motive_separated)
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}

Other functions that define or clean trajectories: 
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Vikram B. Baliga and Melissa S. Armstrong
}
\concept{data cleaning functions}
\concept{functions that define or clean trajectories}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_import_functions.R
\name{read_flydra_mat}
\alias{read_flydra_mat}
\title{Import data from a MAT file exported from Flydra software}
\usage{
read_flydra_mat(mat_file, file_id = NA, subject_name, frame_rate = 100, ...)
}
\arguments{
\item{mat_file}{A file (or path to file) in .mat format, exported from Flydra}

\item{file_id}{(Optional) identifier for this file. If not supplied, this
defaults to \code{basename(file_name)}.}

\item{subject_name}{Name that will be assigned to the subject}

\item{frame_rate}{The capture frame rate of the session}

\item{...}{Additional arguments that may be passed from other pathviewr
functions}
}
\value{
A tibble with numerical data in columns. The first two columns will
have frame numbers and time (assumed to be in secs), respectively. Columns
3 through 5 will contain position data. Note that unlike the behavior of
\code{read_motive_csv()} this function produces "tidy" data that have
already been gathered into key-value pairs based on subject.
}
\description{
\code{read_flydra_mat()} is designed to import data from a \code{.mat} file
that has been exported from Flydra software. The resultant object is a tibble
that additionally has important metadata stored as attributes (see Details).
}
\examples{
library(pathviewr)

## Import the example Flydra data included in the package
flydra_data <-
  read_flydra_mat(system.file("extdata", "pathviewr_flydra_example_data.mat",
                             package = 'pathviewr'),
                  subject_name = "birdie_wooster")

## Names of variables in the resulting tibble
names(flydra_data)

## A variety of metadata are stored as attributes. Of particular interest:
attr(flydra_data, "pathviewr_steps")
}
\seealso{
\code{\link{read_motive_csv}} for importing Motive data

Other data import functions: 
\code{\link{as_viewr}()},
\code{\link{import_and_clean_batch}()},
\code{\link{import_batch}()},
\code{\link{read_motive_csv}()}
}
\author{
Vikram B. Baliga
}
\concept{data import functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_import_functions.R
\name{read_motive_csv}
\alias{read_motive_csv}
\title{Import data from a CSV exported from Optitrack's Motive software}
\usage{
read_motive_csv(file_name, file_id = NA, simplify_marker_naming = TRUE, ...)
}
\arguments{
\item{file_name}{A file (or path to file) in CSV format}

\item{file_id}{(Optional) identifier for this file. If not supplied, this
defaults to \code{basename(file_name)}.}

\item{simplify_marker_naming}{If Markers are encountered, should they be
renamed from "Subject:marker" to "marker"? Defaults to TRUE}

\item{...}{Additional arguments passed from other \code{pathviewr} functions}
}
\value{
A tibble with numerical data in columns. The first two columns will
have frame numbers and time (assumed to be in secs), respectively. Columns 3
and beyond will contain the numerical data on the position or rotation of
rigid bodies and/or markers that appear in the Motive CSV file. Each row
corresponds to the position or rotation of all objects at a given time
(frame).
}
\description{
\code{read_motive_csv()} is designed to import data from a CSV that has been
exported from Optitrack's Motive software. The resultant object is a tibble
that additionally has important metadata stored as attributes (see Details).
}
\details{
Uses \code{data.table::fread()} to import data from a CSV file and
ultimately store it in a tibble. This object is also labeled with the
attribute \code{pathviewr_steps} with value \code{viewr} to indicate that it
has been imported by \code{pathviewr} and should be friendly towards use with
other functions in our package. Additionally, the following metadata are
stored in the tibble's attributes: header information from the Motive CSV
file (\code{header}), original IDs for each object (\code{Motive_IDs}), the
name of each subject in each data column (\code{subject_names_full}) and
unique values of subject names (\code{subject_names_simple}), the type of
data (rigid body or marker) that appears in each column
(\code{data_types_full}) and overall (\code{data_types_simple}), and original
data column names in the CSV (\code{d1, d2}). See Example below for example
code to inspect attributes.
}
\section{Warning}{

This function was written to read CSVs exported using Motive's Format Version
1.23 and is not guaranteed to work with those from other versions. Please
file an Issue on our Github page if you encounter any problems.
}

\examples{
library(pathviewr)

## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

## Names of variables in the resulting tibble
names(motive_data)

## A variety of metadata are stored as attributes. Of particular interest:
attr(motive_data, "pathviewr_steps")
attr(motive_data, "file_id")
attr(motive_data, "header")
attr(motive_data, "Motive_IDs")
attr(motive_data, "subject_names_full")
attr(motive_data, "subject_names_simple")
attr(motive_data, "motive_data_names")
attr(motive_data, "motive_data_types_full")
attr(motive_data, "motive_data_types_simple")

## Of course, all attributes can be viewed as a (long) list via:
attributes(motive_data)

}
\seealso{
\code{\link{read_flydra_mat}} for importing Flydra data

Other data import functions: 
\code{\link{as_viewr}()},
\code{\link{import_and_clean_batch}()},
\code{\link{import_batch}()},
\code{\link{read_flydra_mat}()}
}
\author{
Vikram B. Baliga
}
\concept{data import functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{remove_duplicate_frames}
\alias{remove_duplicate_frames}
\title{Remove any duplicates or aliased frames within trajectories}
\usage{
remove_duplicate_frames(obj_name)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps}.
}
\description{
Remove any duplicates or aliased frames within trajectories
}
\details{
The separate_trajectories() and get_full_trajectories() must be
run prior to use.
}
\seealso{
Other utility functions: 
\code{\link{clean_by_span}()},
\code{\link{insert_treatments}()},
\code{\link{remove_vel_anomalies}()},
\code{\link{set_traj_frametime}()}
}
\author{
Vikram B. Baliga
}
\concept{utility functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_functions.R
\name{import_and_clean_batch}
\alias{import_and_clean_batch}
\title{Batch import and clean files}
\usage{
import_and_clean_batch(
  file_path_list,
  import_method = c("flydra", "motive"),
  file_id = NA,
  subject_name = NULL,
  frame_rate = NULL,
  simplify_marker_naming = TRUE,
  import_messaging = FALSE,
  ...
)
}
\arguments{
\item{file_path_list}{A list of file paths leading to files to be imported.}

\item{import_method}{Either "flydra" or "motive"}

\item{file_id}{(Optional) identifier for this file. If not supplied, this
defaults to \code{basename(file_name)}.}

\item{subject_name}{For Flydra, the subject name applied to all files}

\item{frame_rate}{For Flydra, the frame rate applied to all files}

\item{simplify_marker_naming}{For Motive, if Markers are encountered, should
they be renamed from "Subject:marker" to "marker"? Defaults to TRUE}

\item{import_messaging}{Should this function report each time a file has been
processed?}

\item{...}{Additional arguments to specify how data should be cleaned.}
}
\value{
A list of viewr objects (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) that have been passed
through the corresponding cleaning functions.
}
\description{
Like \code{clean_viewr_batch()}, but with import as the first step too
}
\details{
viewr objects should be in a list, e.g. the object generated by
\code{import_batch()}.

See \code{clean_viewr()} for details of how cleaning steps are handled
and/or refer to the corresponding cleaning functions themselves.
}
\examples{
## Since we only have one example file of each type provided
## in pathviewr, we will simply import the same example multiple
## times to simulate batch importing. Replace the contents of
## the following list with your own list of files to be imported.

## Make a list of the same example file 3x
import_list <-
  c(rep(
    system.file("extdata", "pathviewr_motive_example_data.csv",
                package = 'pathviewr'),
    3
  ))

## Batch import
motive_batch_imports <-
  import_batch(import_list,
               import_method = "motive",
               import_messaging = TRUE)

## Batch cleaning of these imported files
## via clean_viewr_batch()
motive_batch_cleaned <-
  clean_viewr_batch(
    file_announce = TRUE,
    motive_batch_imports,
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )

## Alternatively, use import_and_clean_batch() to
## combine import with cleaning on a batch of files
motive_batch_import_and_clean <-
  import_and_clean_batch(
    import_list,
    import_method = "motive",
    import_messaging = TRUE,
    motive_batch_imports,
    desired_percent = 50,
    max_frame_gap = "autodetect",
    span = 0.95
  )

## Each of these lists of objects can be bound into
## one viewr object (i.e. one tibble) via
## bind_viewr_objects()
motive_bound_one <-
  bind_viewr_objects(motive_batch_cleaned)

motive_bound_two <-
  bind_viewr_objects(motive_batch_import_and_clean)

## Either route results in the same object ultimately:
identical(motive_bound_one, motive_bound_two)
}
\seealso{
Other data import functions: 
\code{\link{as_viewr}()},
\code{\link{import_batch}()},
\code{\link{read_flydra_mat}()},
\code{\link{read_motive_csv}()}

Other batch functions: 
\code{\link{bind_viewr_objects}()},
\code{\link{clean_viewr_batch}()},
\code{\link{import_batch}()}
}
\author{
Vikram B. Baliga
}
\concept{batch functions}
\concept{data import functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{get_velocity}
\alias{get_velocity}
\title{Get instantaneous velocity for subjects}
\usage{
get_velocity(
  obj_name,
  time_col = "time_sec",
  length_col = "position_length",
  width_col = "position_width",
  height_col = "position_height",
  add_to_viewr = TRUE,
  velocity_min = NA,
  velocity_max = NA,
  ...
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{time_col}{Name of the column containing time}

\item{length_col}{Name of the column containing length dimension}

\item{width_col}{Name of the column containing width dimension}

\item{height_col}{Name of the column containing height dimension}

\item{add_to_viewr}{Default TRUE; should velocity data be added as new
columns or should this function create a new simpler object?}

\item{velocity_min}{Should data below a certain velocity be filtered out of
the object? If so, enter a numeric. If not, keep NA.}

\item{velocity_max}{Should data above a certain velocity be filtered out of
the object? If so, enter a numeric. If not, keep NA.}

\item{...}{Additional arguments passed to or from other pathviewr functions.}
}
\value{
If \code{add_to_viewr} is \code{TRUE}, additional columns are
appended to the input viewr object. If \code{FALSE}, a standalone tibble is
created. Either way, an "instantaneous" velocity is computed as the
difference in position divided by the difference in time as each successive
row is encountered. Additionally, velocities along each of the three
position axes are computed and provided as additional columns.
}
\description{
Velocity (both overall and per-axis) is computed for each row in the data
(see Details)
}
\details{
Instantaneous velocity is not truly "instantaneous" but rather is
approximated as the change in distance divided by change in time from one
observation (row) to the previous observation (row). Each component of
velocity is computed (i.e. per axis) along with the overall velocity of
the subject.
}
\examples{
## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                             package = 'pathviewr'))

## Clean the file. It is generally recommended to clean up to the
## "standarization" step before running get_velocity().
 motive_cleaned <-
   motive_data \%>\%
   relabel_viewr_axes() \%>\%
   gather_tunnel_data() \%>\%
   trim_tunnel_outliers() \%>\%
   rotate_tunnel()

## Now compute velocity and add as columns
 motive_velocity_added <-
   motive_cleaned \%>\%
   get_velocity(add_to_viewr = TRUE)

## Or set add_to_viewr to FALSE for a standalone object
 motive_velocity_standalone <-
   motive_cleaned \%>\%
   get_velocity(add_to_viewr = TRUE)
}
\seealso{
Other mathematical functions: 
\code{\link{calc_min_dist_v}()},
\code{\link{deg_2_rad}()},
\code{\link{find_curve_elbow}()},
\code{\link{get_2d_angle}()},
\code{\link{get_3d_angle}()},
\code{\link{get_3d_cross_prod}()},
\code{\link{get_dist_point_line}()},
\code{\link{get_traj_velocities}()},
\code{\link{rad_2_deg}()}
}
\author{
Vikram B. Baliga and Melissa S. Armstrong
}
\concept{mathematical functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{section_tunnel_by}
\alias{section_tunnel_by}
\title{Bin data along a specified axis}
\usage{
section_tunnel_by(obj_name, axis = "position_length", number_of_sections = 20)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{axis}{Chosen axis, must match name of column exactly}

\item{number_of_sections}{Total number of sections}
}
\value{
A new column added to the input data object called \code{section_id},
which is an ordered factor that indicates grouping.
}
\description{
Chop data into X sections (of equal size) along a specified axis
}
\details{
The idea is to bin the data along a specified axis, generally
\code{position_length}.
}
\examples{
## Load data and run section_tunnel_by()
test_mat <-
  read_flydra_mat(system.file("extdata", "pathviewr_flydra_example_data.mat",
                             package = 'pathviewr'),
                  subject_name = "birdie_wooster") \%>\%
  redefine_tunnel_center(length_method = "middle",
                         height_method = "user-defined",
                         height_zero = 1.44) \%>\%
  select_x_percent(desired_percent = 50) \%>\%
  separate_trajectories(max_frame_gap = 1) \%>\%
  get_full_trajectories(span = 0.95) \%>\%
  section_tunnel_by(number_of_sections = 10)

## Plot; color by section ID
plot(test_mat$position_length,
     test_mat$position_width,
     asp = 1, col = as.factor(test_mat$section_id))
}
\author{
Vikram B. Baliga
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{select_x_percent}
\alias{select_x_percent}
\title{Select a region of interest within the tunnel}
\usage{
select_x_percent(obj_name, desired_percent = 33, ...)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{desired_percent}{Numeric, the percent of the total length of the tunnel
that will define the region of interest. Measured from the center outwards.}

\item{...}{Additional arguments passed to/from other pathviewr functions}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) in which data outside
the region of interest have been removed.
}
\description{
Select data in the middle X percent of the length of the tunnel
}
\examples{
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))

## Clean the file. It is generally recommended to clean up to the
## "trimmed" step before running rotate_tunnel().
motive_rotated <-
  motive_data \%>\%
  relabel_viewr_axes() \%>\%
  gather_tunnel_data() \%>\%
  trim_tunnel_outliers() \%>\%
  rotate_tunnel()

## Now select the middle 50\% of the tunnel
motive_selected <-
  motive_rotated \%>\%
  select_x_percent(desired_percent = 50)

## Compare the ranges of lengths to see the effect
range(motive_rotated$position_length)
range(motive_selected$position_length)
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{rotate_tunnel}
\alias{rotate_tunnel}
\title{Rotate a tunnel so that perches are approximately aligned}
\usage{
rotate_tunnel(
  obj_name,
  all_heights_min = 0.11,
  all_heights_max = 0.3,
  perch1_len_min = -0.06,
  perch1_len_max = 0.06,
  perch2_len_min = 2.48,
  perch2_len_max = 2.6,
  perch1_wid_min = 0.09,
  perch1_wid_max = 0.31,
  perch2_wid_min = 0.13,
  perch2_wid_max = 0.35,
  ...
)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"} that has been passed
through \code{relabel_viewr_axes()} and \code{gather_tunnel_data()} (or is
structured as though it has been passed through those functions).}

\item{all_heights_min}{Minimum perch height}

\item{all_heights_max}{Maximum perch height}

\item{perch1_len_min}{Minimum length value of perch 1}

\item{perch1_len_max}{Maximum length value of perch 1}

\item{perch2_len_min}{Minimum length value of perch 2}

\item{perch2_len_max}{Maximum length value of perch 2}

\item{perch1_wid_min}{Minimum width value of perch 1}

\item{perch1_wid_max}{Maximum width value of perch 1}

\item{perch2_wid_min}{Minimum width value of perch 2}

\item{perch2_wid_max}{Maximum width value of perch 2}

\item{...}{Additional arguments passed to/from other pathviewr functions}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) in which data have
been rotated according to user specifications.
}
\description{
The rotation is applied about the height axis and affects tunnel length and
width only, i.e. no rotation of height.
}
\details{
The user first estimates the locations of the perches by specifying
bounds for where each perch is located. The function then computes the
center of each bounding box and estimates that to be the midpoint of each
perch. Then the center point of the tunnel (center between the perch
midpoints) is estimated. The angle between perch1_center,
tunnel_center_point, and arbitrary point along the length axis
(tunnel_center_point - 1 on length) is estimated. That angle is then used
to rotate the data, again only in the length and width dimensions. Height
is standardized by (approximate) perch height; values greater than 0 are
above the perch and values less than 0 are below the perch level.
}
\examples{
## Import the example Motive data included in the package
motive_data <-
  read_motive_csv(system.file("extdata", "pathviewr_motive_example_data.csv",
                              package = 'pathviewr'))

## Clean the file. It is generally recommended to clean up to the
## "trimmed" step before running rotate_tunnel().
motive_trimmed <-
  motive_data \%>\%
  relabel_viewr_axes() \%>\%
  gather_tunnel_data() \%>\%
  trim_tunnel_outliers()

## Now rotate the tunnel using default values
motive_rotated <-
  motive_trimmed \%>\%
  rotate_tunnel()

## The following attributes store information about
## how rotation & translation was applied
attr(motive_rotated, "rotation_degrees")
attr(motive_rotated, "rotation_radians")
attr(motive_rotated, "perch1_midpoint_original")
attr(motive_rotated, "perch1_midpoint_current")
attr(motive_rotated, "tunnel_centerpoint_original")
attr(motive_rotated, "perch2_midpoint_original")
attr(motive_rotated, "perch2_midpoint_current")
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{quick_separate_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}

Other tunnel standardization functions: 
\code{\link{redefine_tunnel_center}()},
\code{\link{standardize_tunnel}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
\concept{tunnel standardization functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{quick_separate_trajectories}
\alias{quick_separate_trajectories}
\title{Quick version of separate_trajectories()}
\usage{
quick_separate_trajectories(obj_name, max_frame_gap = 1, ...)
}
\arguments{
\item{obj_name}{The input viewr object; a tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}}

\item{max_frame_gap}{Unlike the corresponding parameter in
\code{separate_trajectories}, must be a single numeric here.}

\item{...}{Additional arguments passed to/from other pathviewr functions}
}
\value{
A viewr object (tibble or data.frame with attribute
\code{pathviewr_steps} that includes \code{"viewr"}) in which a new column
\code{file_sub_traj} is added, which labels trajectories within the data by
concatenating file name, subject name, and a trajectory number (all
separated by underscores).
}
\description{
Mostly meant for internal use but available nevertheless.
}
\details{
This function is designed to separate rows of data into separately
labeled trajectories.

The \code{max_frame_gap} parameter determines how trajectories will be
separated. \code{max_frame_gap} defines the largest permissible gap in data
before a new trajectory is forced to be defined. In this function, only a
single numeric can be supplied to this parameter (unlike the case in
\code{separate_trajectories}).
}
\examples{
## This function is not recommended for general use.
## See separate_trajectories() instead
}
\seealso{
Other data cleaning functions: 
\code{\link{gather_tunnel_data}()},
\code{\link{get_full_trajectories}()},
\code{\link{redefine_tunnel_center}()},
\code{\link{relabel_viewr_axes}()},
\code{\link{rename_viewr_characters}()},
\code{\link{rotate_tunnel}()},
\code{\link{select_x_percent}()},
\code{\link{separate_trajectories}()},
\code{\link{standardize_tunnel}()},
\code{\link{trim_tunnel_outliers}()},
\code{\link{visualize_frame_gap_choice}()}

Other functions that define or clean trajectories: 
\code{\link{get_full_trajectories}()},
\code{\link{separate_trajectories}()},
\code{\link{visualize_frame_gap_choice}()}
}
\author{
Vikram B. Baliga
}
\concept{data cleaning functions}
\concept{functions that define or clean trajectories}
