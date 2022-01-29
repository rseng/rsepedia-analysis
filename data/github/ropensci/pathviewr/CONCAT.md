
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
