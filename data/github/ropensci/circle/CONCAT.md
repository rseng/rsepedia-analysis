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
---
title: "Using {circle} with {tic}"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{tic}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This vignette explains how {circle} can be used in conjunction with {tic} to set up a working CI environment to check an R package and build a pkgdown site.

All following points assume you are in the project root of the R package.

## Enabling the repository on Circle CI

The first step is to enable/register your repository on Circle CI.
To do this, `enable_repo()` can be called.
Assuming you already have an account on Circle CI (authenticating with GitHub is recommended), this "follows" your repository on Circle CI so that builds can be triggered by pushes to the repository.

## Creating the Circle CI YAML configuration file

Next, a YAML file (`.circleci/config.yml`) needs to be created which lists the tasks which should be executed on Circle CI after a commit to the repository.

This is where the [ropensci/tic](https://github.com/ropensci/tic) package comes into play which provides YAML templates for Circle CI.
There are two ways to get the template via {tic}:

- by going through the chatty `tic::use_tic()` wizard which asks some questions related to configuration and then writes/initiates multiple CI providers (based on the choices made).
  This is a good choice if you want to understand in greater detail what {tic} is doing.
- by calling `tic::use_circle_yml()` which (by default) writes a Circle CI configuration file that checks the R package via {rcmdcheck} and deploys a {pkgdown} site to the `gh-pages` branch of the repository.

In addition some files will be added to `.Rbuildignore` and `.gitignore`.
Also the CI-agnostic `tic.R` file will be created which lists the steps/macros that will be executed in a domain-specific language syntax.
Please have a look at [the introduction vignette of the tic package](https://docs.ropensci.org/tic/articles/tic.html) to understand the role of `tic.R` in more detail.

## Enabling deployment from builds

To be able to push to the GitHub repository some setup work is required.
Deployment is often done by creating a SSH key pair of which one part is stored on Circle CI and the other one on GitHub.
To prevent having to add the SSH key parts to Circle CI and GitHub manually, `use_circle_deploy()` can be called to do all of this programmatically.
See [the section on "Deployment" in the "Getting Started" vignette](https://docs.ropensci.org/circle/articles/circle.html#deployment-1) for more details on this process.

## Understanding the YAML file

The config file of this repo at [.circleci/config.yml](https://github.com/ropensci/circle/blob/main/.circleci/config.yml) has also been set up with {tic}.

Let's walk through it step by step to understand what is happening:

```yml
jobs:
  r-release:
    # r-release-env
    environment:
    docker:
      - image: rocker/verse
    steps:
      - checkout
```

In this part we specify to use the [rocker/verse](https://hub.docker.com/r/rocker/verse) docker image as the base for the job.
The first step is "checkout" which means the repository is cloned.

---

```yml
      # create a unique env var for the cache. Unfortunately normal env vars
      # are not picked up by the cache, therefore this workaround is needed.
      # See https://discuss.circleci.com/t/cannot-use-circle-yml-environment-variables-in-cache-keys/10994/7
      - run: echo "$(date '+%d-%m')-r-release" > /tmp/_tmp_file
      - restore_cache:
          key: R-package-library-{{ checksum "/tmp/_tmp_file" }}
```

Next, an action related to caching R packages is initiated.
This saves some time in the future because once all R packages which the package needs have been installed once for a given day, they will be re-used in future builds of the day without having to be installed again.

---

```yml
      # install deps and check pkg ---------------------------------------------
      - run:
          name: "[r-release] Install dependencies"
          command: |
            sudo apt update && sudo apt install -y ccache libgit2-dev libharfbuzz-dev libfribidi-dev
            echo -e 'options(Ncpus = 4, repos = structure(c(CRAN = "https://cloud.r-project.org/")))' > $HOME/.Rprofile
            R -q -e 'install.packages("remotes")'
            R -q -e 'if (getRversion() < "3.2" && !requireNamespace("curl")) install.packages("curl")'
            R -q -e 'remotes::install_github("ropensci/tic", upgrade = "always"); print(tic::dsl_load()); tic::prepare_all_stages()'
            R -q -e 'tic::before_install()'
            R -q -e 'tic::install()'
```

Next, the {tic} package is installed and certain [tic steps](https://docs.ropensci.org/tic/articles/build-lifecycle.html) are run.
These take care of installing the dependencies of the R package to be checked and prepare the environment of other subsequent steps.

---

```yml
      - run:
          name: "[r-release] R CMD Check"
          no_output_timeout: 60m
          command: |
            R -q -e 'tic::before_script()'
            R -q -e 'tic::script()'
```

This step checks the package for CRAN eligibility by making use of `rcmdcheck::rcmdcheck()` inside the `tic::script()` call.

---

```yml
      # save R pkg cache -------------------------------------------------------
      - save_cache:
          key: R-package-library-{{ checksum "/tmp/_tmp_file" }}
          paths:
            - /usr/local/lib/R/site-library
```

Finally, the R package which was initiated earlier is saved.

---

The `deploy:` step following next is in most parts executing the same steps as just shown.
In the end however, `tic::deploy()` is called which internally will build a {pkgdown} site of the package and then deploy this site to the `gh-pages` branch of the repository.

## The first build

After `.circleci/config.yml` and `tic.R` have been commited and pushed, the first build will start on Circle CI.

Calling one of `get_pipelines()`, `get_workflows()` or `get_jobs()` should now return some content.

In addition, you can directly browse the builds in the Circle CI web interface by calling `usethis::browse_circleci()`.
---
title: "Getting Started"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{circle}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction

To be able to run builds from a repository, it needs to be enabled on Circle CI first.
This package heavily relies on GitHub for automatic user information extraction and other providers such as GitLab or Bitbucket are not supported.
In general Circle CI builds can be run from any Git hosting provider.

# Authentication with GitHub

To get started with Circle CI on GitHub, make sure to install the [Circle CI GitHub App](https://GitHub.com/marketplace/circleci) and configure the "Free" plan for your account.

The first time any {circle} function is called, it will check for the existence of a Circle API key.
This API key is needed to securely talk to the Circle CI API and make requests on behalf of a user.
There are two ways the API key can be set:

  - (**recommended**) Via environment variable `R_CIRCLE` in `~/.Renviron` (e.g. `R_CIRCLE = <API KEY>`)
  - In `~/.circleci/cli.yml` via the following pattern

   ```yml
   host: https://circleci.com
   endpoint: graphql-unstable
   token: <token>
   ```

To create an API key from the R console, `browse_circle_token()` can be used.

{circle} scrapes information about the current repository by making API calls to GitHub.
To be able to do this, a `GitHub_TOKEN` is needed (similar to the API requests for Circle CI).
It is good practice to have such a token set for many purposes when using R with GitHub.
Invoking `usethis::browse_GitHub_token()` is an easy way to create one if none exists yet.

# First steps

By querying information about the user and enabling a project one can check if everything was set up correctly.

`get_circle_user()` makes uses of the Circle CI API key.


```r
circle::get_circle_user()

$content
$content$name
[1] "Patrick Schratz"

$content$login
[1] "pat-s"

$content$id
[1] "9c373331-d0f7-45e1-afe6-4a5c75e00d10"


$path
[1] "/me"

$response
Response [https://circleci.com/api/v2/me?circle-token=39d697f345d8d8a92ab07c333405d9b0092d116c]
  Date: 2021-01-07 20:07
  Status: 200
  Content-Type: application/json;charset=utf-8
  Size: 86 B


attr(,"class")
[1] "circle_user"
```

`enable_repo()` also uses the GitHub token.

```r
circe::enable_repo()
✔ Successfully enabled repo 'ropensci/circle' on Circle CI.
```

After the repo has been enabled, it should be returned in `circle::list_projects()`.

# Deployment

Deployments refers to the practice to push (modified) files during a build to a repository.
Configuring build deployments can be a bit tedious with respect to permissions.
Here `circle::use_circle_deploy()` helps as it creates a SSH key pair which will enable deployment.
The private key will be stored in GitHub (under "Settings -> SSH and GPG keys") and the public key on Circle CI (on the respective project page in "SSH Keys").

```r
circle::use_circle_deploy()
✖ No 'user-key' found.
──────────────────────────────────────────────────────────────────────────────────────────────────
✔ Added a 'user key' to project '<repo-slug>' on Circle CI.
This enables deployment from builds.
```

To double-check, there should now be a "user-key" in your Circle CI repo settings under the menu point "SSH keys".
It does not matter if multiple "deploy key" or "user-key"s exists.
The important point is that one "user key" is set as "preferred".

```{r, echo=FALSE}
knitr::include_graphics("../man/figures/user-key.png")
```

As an alternative to `circle::use_circle_deploy()` one could also click the "Add user key" button that would appear if no user key has been set yet.

## Deployment Keys

With respect to SSH keys, there two different types of keys on Circle CI:

- Deploy key

- User key

"Deploy keys" are used to checkout your repository so that the build is starting.
These kind of keys how only "read" permissions but are not allow to "write" to your repository.
If you have connected Circle CI to GitHub already, you will have a "deploy key" stored in every repository to be able to checkout the code.
The name "deploy" key here is a bit misleading since they **cannot** be used for deploying from builds.

To enable deployment to a repository however, a "user-key" is needed.
This key type also has "write" permissions and can be added using `use_circle_deploy()`.
`use_circle_deploy()` will add a so called "user key" to the settings of the repo on Circle CI.
The private key will be added to your GitHub profile under the "SSH and GPG keys" section with the title pointing to the respective repo.
See also [the Circle CI section about deployment](https://circleci.com/docs/2.0/deployment-integrations/#section=deployment).

(If for some reasons you do not want to use `use_circle_deploy()` and go the manual way of adding a SSH key to Circle CI, please be aware of [this issue](https://discuss.circleci.com/t/adding-ssh-keys-fails/7747).)

# Starting a Build

As with almost every CI provider, a YAML configuration files is required to specify the tasks that should be executed during a build.
Currently Circle CI does not come with official support for the R language but you can add [your vote here](https://ideas.circleci.com/images/p/add-official-support-for-r).
Since Circle CI is heavily based on docker this is not really a problem.
One can simply use the [rocker](https://github.com/rocker-org/rocker) R images to have first-class support with respect to R containers.
Official support for R would mean that the images are cached on Circle CIs side and be directly available on build start.
Currently, `rocker` images need to be downloaded in every build again.

{circle} does not come with a template for running R builds as it focuses on API and build metadata functionality.
Instead, have a look at [ropensci/tic](https://docs.ropensci.org/tic/) and its functions `tic::use_tic()` and `tic::use_circle_yml()` to quickly get R builds running on Circle CI.
Alternatively you can borrow the config that is being used [in this repo](https://github.com/ropensci/circle/blob/main/.circleci/config.yml) - which was also created via {tic}.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers-github.R
\name{github_info}
\alias{github_info}
\title{Github information}
\usage{
github_info(path = ".", remote = "origin", .token = NULL)
}
\arguments{
\item{path}{\verb{[string]}\cr
The path to a GitHub-enabled Git repository (or a subdirectory thereof).}

\item{remote}{\verb{[character]}\cr
The Github remote which should be used.}

\item{.token}{\verb{[character]}\cr
Authentication token. Defaults to GITHUB_PAT or GITHUB_TOKEN environment
variables, in this order if any is set. See gh_token() if you need more
flexibility, e.g. different tokens for different GitHub Enterprise
deployments.}
}
\value{
Object of class \code{gh_response} (list type) with information about the
queried repository.
}
\description{
Retrieves metadata about a Git repository from GitHub.

\code{github_info()} returns a list as obtained from the GET "/repos/:repo" API.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use-circle-deploy.R
\name{use_circle_deploy}
\alias{use_circle_deploy}
\title{Set Up Build Deployment Between Circle CI And Github}
\usage{
use_circle_deploy(
  repo = github_info()$name,
  user = github_info()$owner$login,
  quiet = FALSE
)
}
\arguments{
\item{repo}{\verb{[character]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{user}{\verb{[character]}\cr
The username for the repository. By default queried using \code{get_user()}.}

\item{quiet}{\verb{[logical]}\cr
If \code{TRUE}, console output is suppressed.}
}
\value{
No return value, called for side effects.
}
\description{
Creates a Circle CI "user-key" (= SSH key pair) if none exists
yet to enable deployment from Circle CI builds to GitHub.
}
\details{
The easiest way to achieve a deployment from Circle CI builds to a Github
repo is by creating a so called "user-key" (i.e. an SSH key pair) on
Circle CI.

\code{use_circle_deploy()} tries to be smart by exiting early if such a key is
already present.

If the repo has not been enabled yet on Circle CI, please run \code{enable_repo()}
first.
Also to be able to authenticate to Github in the first place a personal
access token needs to be set (via env var \code{GITHUB_TOKEN}).
\code{usethis::github_token()} can be used to check if one is already set.
If none is set, this function will prompt you to create one.
}
\examples{
\dontrun{
use_circle_deploy()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{browse_circle_token}
\alias{browse_circle_token}
\title{Authenticate to Circle CI}
\usage{
browse_circle_token()
}
\value{
Returns \code{TRUE} (invisibly).
}
\description{
A Circle CI API token is needed to interact with the Circle CI API.
\code{browse_circle_token()} opens a browser window for the respective Circle CI
endpoint to retrieve the key.
}
\section{Store API token}{


\code{circle} supports two ways of storing the Circle API tokens:
\itemize{
\item via env vars \code{R_CIRCLE}
\item via \verb{~/.circleci/cli.yml}
}

The latter should already be present if you already used the \code{circle} CLI
tool at some point in the past. If not, its up to your preference which
approach to use.

The following instructions should help to set up \verb{~/.circleci/cli.yml}
correctly:
\enumerate{
\item Copy the token from the browser after having called
\code{browse_circle_token()}. You can use
\code{edit_circle_config()} to open \verb{~/.circleci/cli.yml}.
\item The token should be stored using the following structure\if{html}{\out{<div class="sh">}}\preformatted{host: https://circleci.com
endpoint: graphql-unstable
token: <token>
}\if{html}{\out{</div>}}
}
}

\examples{
\dontrun{
browse_circle_token()

edit_circle_config()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{edit_circle_config}
\alias{edit_circle_config}
\title{Open circle Configuration file}
\usage{
edit_circle_config()
}
\value{
No return value, called for side effects.
}
\description{
Opens \verb{~/.circleci/cli.yml}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{new_build}
\alias{new_build}
\title{Trigger a New Build on Circle CI}
\usage{
new_build(
  repo = github_info()$name,
  user = github_info()$owner$login,
  vcs_type = "gh",
  branch = "master",
  quiet = FALSE
)
}
\arguments{
\item{repo}{\verb{[character]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{user}{\verb{[character]}\cr
The username for the repository. By default queried using \code{get_user()}.}

\item{vcs_type}{\verb{[character]}\cr The version control system to use.
Defaults to "gh" (Github).}

\item{branch}{A character string specifying the repository branch.}

\item{quiet}{\verb{[logical]}\cr
If \code{TRUE}, console output is suppressed.}
}
\value{
An object of class \code{circle_api} with the following elements
\itemize{
\item \code{content} (queried content)
\item \code{path} (API request)
\item \code{response} (HTTP response information)
}
}
\description{
Triggers a new build for a specific repo branch.
}
\details{
Trigger a new Circle CI build for a specific repo branch.
}
\examples{
\dontrun{
new_build()
}
}
\seealso{
\code{\link[=retry_workflow]{retry_workflow()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{get_circle_user}
\alias{get_circle_user}
\title{Get Circle CI user}
\usage{
get_circle_user()
}
\value{
A named vector of class \code{circle_user} containing information about
about the authenticated user:
\itemize{
\item Full name
\item Username
\item ID
\item API endpoint
}
}
\description{
Retrieve details about the authenticated Circle CI user.
}
\examples{
\dontrun{
get_circle_user()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{list_projects}
\alias{list_projects}
\title{List Circle CI Projects}
\usage{
list_projects(repo = github_info()$name, user = github_info()$owner$login)
}
\arguments{
\item{repo}{\verb{[character]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{user}{\verb{[character]}\cr
The username for the repository. By default queried using \code{get_user()}.}
}
\value{
An object of class \code{circle_api} with the following elements
\itemize{
\item \code{content} (queried content)
\item \code{path} (API request)
\item \code{response} (HTTP response information)
}
}
\description{
Retrieve a list of Circle CI repositories for the authenticated
user.
}
\details{
Retrieves a very detailed list of repository and repo-related
information for all Circle CI repository attached to the current user.

This endpoint uses API v1.1 and will probably be removed in the near
future.
}
\examples{
\dontrun{
list_projects()
}
}
\seealso{
\code{\link[=get_pipelines]{get_pipelines()}}, \code{\link[=get_workflows]{get_workflows()}}, \code{\link[=get_jobs]{get_jobs()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/builds.R
\name{builds}
\alias{builds}
\alias{get_pipelines}
\alias{get_workflows}
\alias{get_jobs}
\alias{retry_workflow}
\title{Retrieve Metadata from Circle CI Builds}
\usage{
get_pipelines(
  repo = github_info()$name,
  user = github_info()$owner$login,
  limit = 30,
  vcs_type = "gh",
  api_version = "v2"
)

get_workflows(
  pipeline_id = NULL,
  repo = github_info()$name,
  user = github_info()$owner$login
)

get_jobs(
  workflow_id = NULL,
  repo = github_info()$name,
  user = github_info()$owner$login,
  vcs_type = "gh"
)

retry_workflow(workflow_id = NULL)
}
\arguments{
\item{repo}{\verb{[character]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{user}{\verb{[character]}\cr
The username for the repository. By default queried using \code{get_user()}.}

\item{limit}{\verb{[integer]}\cr
How many builds should be returned? Maximum allowed by Circle is
30.}

\item{vcs_type}{\verb{[character]}\cr The version control system to use.
Defaults to "gh" (Github).}

\item{api_version}{\verb{[character]}\cr
A character string specifying the Circle CI API version.
This usually does not need to be changed by the user.}

\item{pipeline_id}{\verb{[character]}\cr
A Circle CI pipeline ID.}

\item{workflow_id}{\verb{[character]}\cr
A Circle CI workflow ID.}
}
\value{
An object of class \code{circle_collection} containing list
information on the queried Circle CI pipelines/workflows/jobs.
}
\description{
Query information about pipelines, workflows or jobs on
Circle CI.
The S3 \code{print()} method for these functions returns the respective
pipeline IDs.
To inspect the details of each pipeline, save the return value in an object
and inspect the respective sub-lists.

If no pipeline or workflow is supplied to \code{get_workflows()}/\code{get_jobs()},
the ten most recent pipelines/jobs are queried, respectively.
}
\details{
While the \verb{get_*()} functions query information about the respective
build level details (pipeline - workflow - job), \code{retry_workflow()} let's
users rerun a specific workflow.
By default, the workflow from the most recent pipeline will be rerun if
no pipeline ID was supplied.
}
\examples{
\dontrun{
pipelines <- get_pipelines()

workflows <- get_workflows()

jobs <- get_jobs()

# rerun most recent workflow
retry_workflow()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circle-package.R
\docType{package}
\name{circle-package}
\alias{circle-package}
\title{Circle CI API Client}
\description{
This package provides functionality for interacting with the
Circle CI API. \href{https://circleci.com}{Circle CI} is a continuous
integration provider which allows for
automated testing of software each time that software is publicly committed
to a repository on GitHub.

This package interacts with the Circle CI REST API and allows to execute
tasks in R without visiting the the website. This includes monitoring
builds, modifying build environment settings and environment variables, and
cancelling or restarting builds.

Use of this package requires a Circle API key. Unless a key is already set,
users will be guided through the creation of a key,
API keys are disposable, but should still be treated securely.

The following functions simplify integrating R package testing and
deployment with GitHub and Circle CI:
\itemize{
\item \code{\link[=enable_repo]{enable_repo()}} enables Circle CI for your repository,
\item \code{\link[=use_circle_deploy]{use_circle_deploy()}} installs a public deploy key on GitHub and the
corresponding private key on Circle CI to simplify deployments to GitHub
from Circle CI.
}
}
\examples{
\dontrun{
# check to see if you've authenticated correctly
get_circle_user()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkout-key.R
\name{checkout_key}
\alias{checkout_key}
\alias{create_checkout_key}
\alias{get_checkout_keys}
\alias{delete_checkout_key}
\alias{has_checkout_key}
\title{Interact with "Checkout Keys" on Circle CI}
\usage{
create_checkout_key(
  repo = github_info()$name,
  user = github_info()$owner$login,
  type = "user-key",
  api_version = "v2",
  vcs_type = "gh",
  quiet = FALSE
)

get_checkout_keys(
  repo = github_info()$name,
  user = github_info()$owner$login,
  vcs_type = "gh",
  api_version = "v2"
)

delete_checkout_key(
  fingerprint = NULL,
  repo = github_info()$name,
  user = github_info()$owner$login,
  type = "user-key",
  api_version = "v2",
  vcs_type = "gh"
)

has_checkout_key(
  repo = github_info()$name,
  user = github_info()$owner$login,
  type = "github-user-key",
  vcs_type = "gh",
  preferred = TRUE
)
}
\arguments{
\item{repo}{\verb{[character]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{user}{\verb{[character]}\cr
The username for the repository. By default queried using \code{get_user()}.}

\item{type}{\verb{[character]}\cr
Type of key to add. Options are "user-key" and "deploy-key".}

\item{api_version}{\verb{[character]}\cr
A character string specifying the Circle CI API version.
This usually does not need to be changed by the user.}

\item{vcs_type}{\verb{[character]}\cr The version control system to use.
Defaults to "gh" (Github).}

\item{quiet}{\verb{[logical]}\cr
If \code{TRUE}, console output is suppressed.}

\item{fingerprint}{\verb{[character]}\cr
The fingerprint of the checkout key which should be deleted.}

\item{preferred}{\verb{[logical]}\cr
Checks whether the requested type is the "preferred" key.}
}
\value{
An object of class \code{circle_api} with the following elements
\itemize{
\item \code{content} (queried content)
\item \code{path} (API request)
\item \code{response} (HTTP response information)
}
}
\description{
Create, delete, query or check different types of checkout keys
for a specific Circle CI project.
Valid values for argument \code{type} are \code{"user-key"} or \code{"deploy-key"}.

A "Checkout Key" on Circle CI is a specific SSH key which is used to checkout
repositories into a Circle CI build and possible deploy changes to the
repository.
Circle CI subdivides "Checkout Keys" into "user-key" and "deploy-key".

Please see "Deployment" section in the "Getting Started" vignette for more
information.
}
\examples{
\dontrun{
# by default a "user-key" will be created which can also be used for build
# deployments
create_checkout_key()

# A "deploy-key" can only be used to checkout code from the repository into
# a Circle CI build
create_checkout_key(type = "deploy-key")
}
\dontrun{
get_checkout_keys()
}
\dontrun{
delete_checkout_key()
}
\dontrun{
has_checkout_key()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circle.R
\name{circle}
\alias{circle}
\title{Circle CI HTTP Requests}
\usage{
circle(
  verb = "GET",
  path = "",
  query = list(),
  body = "",
  api_version = "v2",
  encode = "json"
)
}
\arguments{
\item{verb}{\verb{[character]}\cr
A character string containing an HTTP verb, defaulting to \code{GET}.}

\item{path}{\verb{[character]}\cr
A character string with the API endpoint (should begin with a slash).}

\item{query}{\verb{[character]}\cr
A list specifying any query string arguments to pass to the API.
This is used to pass the API token.}

\item{body}{\verb{[character]}\cr
A named list or array of what should be passed in the
request.
Corresponds to the "-d" argument of the \code{curl} command.}

\item{api_version}{\verb{[character]}\cr
A character string specifying the Circle CI API version.
This usually does not need to be changed by the user.}

\item{encode}{\verb{[character]}\cr
Encoding format. See \link[httr:POST]{httr::POST}.}
}
\value{
An object of class \code{circle_api} with the following elements
\itemize{
\item \code{content} (queried content)
\item \code{path} (API request)
\item \code{response} (HTTP response information)
}
}
\description{
Workhorse function for executing API requests to
Circle CI.
}
\details{
In almost all cases, users should not need to execute API calls
directly. However, if desired this functions makes it possible to issue
any API request. If you experience calling a custom request heavily,
consider opening a feature request on GitHub.
}
\examples{
\dontrun{
circle(verb = "GET", path = "/project/gh/ropensci/circle/checkout-key")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/env-var.R
\name{env_var}
\alias{env_var}
\alias{get_env_vars}
\alias{set_env_var}
\alias{delete_env_var}
\title{Interact with Environment Variable(s) on Circle CI}
\usage{
get_env_vars(
  name = NULL,
  repo = github_info()$name,
  user = github_info()$owner$login,
  vcs_type = "gh",
  api_version = "v2"
)

set_env_var(
  var,
  repo = github_info()$name,
  user = github_info()$owner$login,
  vcs_type = "gh",
  api_version = "v2",
  quiet = FALSE
)

delete_env_var(
  var,
  repo = github_info()$name,
  user = github_info()$owner$login,
  vcs_type = "gh",
  api_version = "v2",
  quiet = FALSE
)
}
\arguments{
\item{name}{\verb{[character]}\cr
Name of a specific environment variable.
If not set, all env vars are returned.}

\item{repo}{\verb{[character]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{user}{\verb{[character]}\cr
The username for the repository. By default queried using \code{get_user()}.}

\item{vcs_type}{\verb{[character]}\cr The version control system to use.
Defaults to "gh" (Github).}

\item{api_version}{\verb{[character]}\cr
A character string specifying the Circle CI API version.
This usually does not need to be changed by the user.}

\item{var}{\verb{[list]}\cr
A list containing key-value pairs of environment variable and its value.}

\item{quiet}{\verb{[logical]}\cr
If \code{TRUE}, console output is suppressed.}
}
\value{
An object of class \code{circle_api} with the following elements
\itemize{
\item \code{content} (queried content)
\item \code{path} (API request)
\item \code{response} (HTTP response information)
}
}
\description{
Add, get or set Circle CI environment variable(s) for a repo on
Circle CI.
}
\examples{
\dontrun{
# get env var
get_env_vars()

# set env var
set_env_var(var = list("foo" = "123"))

# delete env var
delete_env_var("foo")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{get_build_artifacts}
\alias{get_build_artifacts}
\title{Get Build Artifacts of a Specific Job}
\usage{
get_build_artifacts(
  job_id = NULL,
  repo = github_info()$name,
  user = github_info()$owner$login,
  vcs_type = "gh",
  api_version = "v2"
)
}
\arguments{
\item{job_id}{\verb{[character]}\cr
A Circle CI job id.}

\item{repo}{\verb{[character]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{user}{\verb{[character]}\cr
The username for the repository. By default queried using \code{get_user()}.}

\item{vcs_type}{\verb{[character]}\cr The version control system to use.
Defaults to "gh" (Github).}

\item{api_version}{\verb{[character]}\cr
A character string specifying the Circle CI API version.
This usually does not need to be changed by the user.}
}
\value{
An object of class \code{circle_api} with the following elements
\itemize{
\item \code{content} (queried content)
\item \code{path} (API request)
\item \code{response} (HTTP response information)
}
}
\description{
Retrieve artifacts from a specific build.
}
\examples{
\dontrun{
job_id <- get_jobs()[[1]]$id
get_build_artifacts(job_id)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general.R
\name{enable_repo}
\alias{enable_repo}
\title{Enable a repo on Circle CI}
\usage{
enable_repo(
  repo = github_info()$name,
  user = github_info()$owner$login,
  vcs_type = "gh",
  api_version = "v1.1",
  quiet = FALSE
)
}
\arguments{
\item{repo}{\verb{[character]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{user}{\verb{[character]}\cr
The username for the repository. By default queried using \code{get_user()}.}

\item{vcs_type}{\verb{[character]}\cr The version control system to use.
Defaults to "gh" (Github).}

\item{api_version}{\verb{[character]}\cr
A character string specifying the Circle CI API version.
This usually does not need to be changed by the user.}

\item{quiet}{\verb{[logical]}\cr
If \code{TRUE}, console output is suppressed.}
}
\value{
An object of class \code{circle_api} with the following elements
\itemize{
\item \code{content} (queried content)
\item \code{path} (API request)
\item \code{response} (HTTP response information)
}
}
\description{
"Follows" a repo on Circle CI so that builds can be triggered.
}
\examples{
\dontrun{
enable_repo()
}
}
