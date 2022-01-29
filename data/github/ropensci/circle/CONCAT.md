<!-- badges: start -->

[![tic](https://github.com/ropensci/circle/workflows/tic/badge.svg?branch=main)](https://github.com/ropensci/circle/actions)
[![CircleCI](https://img.shields.io/circleci/build/gh/ropensci/circle/main?label=Linux&logo=circle&logoColor=green&style=flat-square)](https://circleci.com/gh/ropensci/circle)
[![CRAN Status](https://www.r-pkg.org/badges/version-ago/circle)](https://cran.r-project.org/package=circle)
[![codecov](https://codecov.io/gh/ropensci/circle/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/circle)
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![](https://badges.ropensci.org/356_status.svg)](https://github.com/ropensci/software-review/issues/356)

<!-- badges: end -->

# circle

R client package for the Continuous Integration (CI) provider 'Circle CI'.
[Circle CI](https://circleci.com/) stands in line with [GitHub Actions](https://github.com/features/actions), [Travis CI](https://www.travis-ci.com/), [AppVeyor](https://ci.appveyor.com/login) and many more CI providers.
[Circle CI](https://circleci.com/) heavily relies on Docker containers for runner execution.

Continuous Integration (CI) / Continuous Deployment (CD) is heavily used in the IT world to automatically perform certain actions after a specific trigger (e.g. after each commit).
When developing R packages the most common uses cases are to check the package on each commit for CRAN eligibility (by running `R CMD Check`) and to deploy a [{pkgdown}](https://github.com/r-lib/pkgdown) documentation page for the package.

This package aims help to set up CI/CD with the service provider [Circle CI](https://circleci.com/) and provides R functions to execute CI specific tasks such as build restarts, log queries or setting environment variables from within R.
It also simplifies the setup process for build deployments via [`use_circle_deploy()`](https://docs.ropensci.org/circle/reference/use_circle_deploy.html).
All functionality relies on calls to the [Circle CI REST API](https://circleci.com/docs/api/v2/#circleci-api).

There are two ways to use this package:

- Via the high-level functions of this package which wrap common API calls:
  - `get_pipelines()`
  - `get_checkout_keys()`
  - `set_env_var()`
  - etc.
- Via direct API calls through the workhorse function [`circle()`](https://docs.ropensci.org/circle/reference/circle.html).

{circle} does not come with an option to setup Circle CI YAML files.
Please see the related [{tic}](https://github.com/ropensci/tic) package for such functionality and more CI workflow related tools.
{circle} aims to provide a handy and flexible high-level interface to the [Circle CI API](https://circleci.com/docs/api/v2/) without shipping opinionated workflow functionality.

## API versions

All functionality uses the Circle CI [API v2](https://github.com/CircleCI-Public/api-preview-docs) which follows the **pipelines** -> **workflows** -> **jobs** approach.
This API version is still in beta and might undergo some changes in the near future.

Some functions/endpoints can also be used via API versions v1.1 and v1 by setting the `api_version` argument.
However, this will only work if the respective API endpoint is available for the chosen API version.
Usually, there should be no need in practice to fall back to API version < 2.

For more information on the differences between the [Circle CI API](https://circleci.com/docs/api/v2/) versions, have a look at the [document explaining changes between v1.1 and v2](https://github.com/CircleCI-Public/api-preview-docs/blob/master/docs/api-changes.md).

## Installation

Development Version:

```r
remotes::install_github("ropensci/circle")
```

## Get Started

See the [Getting Started](https://ropensci.github.io/circle/articles/circle.html) vignette for an introduction.

## Note to Developers

This packages relies on private API keys for local testing.
See [CONTRIBUTING.md#testing-the-package](https://github.com/ropensci/circle/blob/main/.github/CONTRIBUTING.md#testing-the-package) for detailed instructions.

# Acknowledgments

This package was inspired by the work of [Thomas J. Leeper](https://github.com/leeper) on the (discontinued) [cloudyr/circleci](https://github.com/cloudyr/circleci) package and by the (archived) [ropenscilabs/travis](https://github.com/ropensci-archive/travis) package.
<!-- NEWS.md is maintained by https://cynkra.github.io/fledge, do not edit -->

# circle 0.7.1.9000

- Same as previous version.


# circle 0.7.1

- Initial CRAN release


# circle 0.7.0

Implement feedback from [ropensci review](https://github.com/ropensci/software-review/issues/356#):

- Document return values of all functions
- Refine {cli} console messages
- Most functions gained a `quiet` argument to silence {cli} messages
- Be more chatty for side-effect functions
- Always return a `circle_api` object for consistency
- Switch main branch from `master` to `main`
- Escape examples
- Require {usethis} >= 2.0.0
- New vignette ["Using {circle} with {tic}"](https://docs.ropensci.org/circle/articles/tic.html)


# circle 0.6.0

- Copy over GitHub auth and SSH helpers from {travis}
- Print informative message when creating a user key errors with status code 500
- `*_env_var()`: Use owner info instead of user info to query repo
- Use {vcr} for http testing
- Add pkgdown reference structure
- Added pre-commit hooks
- Added codemeta
- Use roxygen markdown
- Added parameter types to help pages

# circle 0.5.0

## Major

- Add new authentication mechanism: `browse_circle_token()` to to query the API token and store it in an env variable `R_CIRCLE` as an alternative method to store it in `~/.circleci/cli.yml`
- Remove `auth_travis()`
- Rename `circleHTTP()` to `circle()`
- add `github_repo()`
- `get_pipelines()`, `get_workflows()` and `get_jobs()` are now formatted as class `circle_builds`, `circle_collection()` and have a somewhat pretty print output
- `*_checkout_key()`: Optimize printer, catch errors, add info messages, add test
- make `get_pipelines()`, `get_workflows()` and `get_jobs()` work with API v2
- rename `list_artifacts()` -> `get_build_artifacts()`

## Bugfixes

- Pipelines without a workflow ID caused `get_builds()` to error. Now pipelines without a workflow ID are removed internally before continuing.
- setting env vars now works
- make `create_checkout_key()` work with API v2

# circle 0.4.0

- update "cache" function with new user/owner logic from v0.3.0
- new `has_checkout_key()` to check if a specific checkout key exists in the project

# circle 0.3.0

- Rename argument `project` to `repo` to stay consistent with _travis_ pkg.

- Add Github helper functions to easily query owners and users for the repository operating on. This change requires the _git2r_ package from now on.

# circle 0.2.0

- Fix `api_version` in `create_ssh_key()`
- rename `ssh_key*` functions to `checkout_key*`
- `create_checkout_key()` change default for arg `type` from "github-user-key" to "deploy-key"
- add argument `encode` to `circleHTTP()`
- add `use_circle_deploy()`

# circleci 0.1.0

- First working version
circle 0.7.1

## Cran Repository Policy

- [x] Reviewed CRP last edited 2020-10-29.

## R CMD check results

- [x] Checked locally, R 4.0.5
- [x] Checked on CI system, R 4.0.5
- [x] Checked on win-builder, R devel

Check the boxes above after successful execution and remove this line. Then run `fledge::release()`.

## Current CRAN check results

Initial release.
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
# CONTRIBUTING

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

- YES: you edit a roxygen comment in a `.R` file below `R/`.
- NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

- We recommend that you create a Git branch for each pull request (PR).
- Look at the Travis and AppVeyor build status before and after making changes.
  The `README` should contain badges for any continuous integration services used
  by the package.
- We recommend the tidyverse [style guide](http://style.tidyverse.org).
  You can use the [styler](https://CRAN.R-project.org/package=styler) package to
  apply these styles, but please don't restyle code that has nothing to do with
  your PR.
- We use [roxygen2](https://cran.r-project.org/package=roxygen2).
- We use [testthat](https://cran.r-project.org/package=testthat). Contributions
  with test cases included are easier to accept.
- For user-facing changes, add a bullet to the top of `NEWS.md` below the
  current development version header describing the changes made followed by your
  GitHub username, and links to relevant issue(s)/PR(s).

### Testing the package

To test the package locally, the following conditions need to be met:

1. A Circle CI API Token stored in an env var named `R_CIRCLE`
1. An active repository on Circle CI with access rights

Alternatively changes can be pushed to GitHub and being tested by the CI runner.
The CI runner operates on the `ropensci/circle` repo with access to the respective secrets during the run.

### Code of Conduct

Please note that the circle project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html)

for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

- you have a question, an use case, or otherwise not a bug or feature request for the software itself.
- you think your issue requires a longer form discussion.

### Prefer to Email?

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
