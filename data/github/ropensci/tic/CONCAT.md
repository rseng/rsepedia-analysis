# tic

<!-- badges: start -->

[![tic](https://github.com/ropensci/tic/workflows/tic/badge.svg?branch=master)](https://github.com/ropensci/tic/actions)
[![Travis build status](https://travis-ci.org/ropensci/tic.svg?branch=master)](https://travis-ci.org/ropensci/tic)
[![CircleCI](https://img.shields.io/circleci/build/gh/ropensci/tic/master?label=Linux&logo=circle&logoColor=green&style=flat-square)](https://circleci.com/gh/ropensci/tic)
[![CRAN status](https://www.r-pkg.org/badges/version/tic)](https://cran.r-project.org/package=tic)
[![codecov](https://codecov.io/gh/ropensci/tic/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/tic)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://badges.ropensci.org/305_status.svg)](https://github.com/ropensci/software-review/issues/305)

<!-- badges: end -->

The goal of tic is to enhance and simplify working with continuous integration (CI) systems.

The following ones are supported:

| Provider       | R package                                            | Platforms             | Info                                                                  |
| -------------- | ---------------------------------------------------- | --------------------- | --------------------------------------------------------------------- |
| Circle CI      | [{circle}](https://docs.ropensci.org/circle/)        | Linux                 | via Docker images from [rocker](https://github.com/rocker-org/rocker) |
| Github Actions | [{ghactions}](https://maxheld.de/ghactions)          | Linux, macOS, Windows |                                                                       |

To learn more about CI, read our [Getting Started](https://docs.ropensci.org/tic/articles/tic.html#prerequisites) vignette.

The most important improvements over existing solutions are:

1. Deployment to a Git repository is greatly simplified. Update your repository with results from the CI build.

2. Support for R packages and other kinds of projects (bookdown, blogdown, etc.), with predefined templates.
   Set up your project to deploy rendered versions of your book or blog with a single push to Git.

3. Workflow specification in a single `.R` file, regardless of CI system used.
   Forget about `.yml` files or web browser configurations.

## Installation

{tic} can be installed from GitHub with:

```r
remotes::install_github("ropensci/tic")
```

## Setup

By calling `tic::use_tic()` a production ready CI setup is initialized, tailored to your specific R project.
The created templates will use the providers [Circle CI](https://circleci.com) and [Github Actions](https://github.com/actions).

If only the CI YAML templates from {tic} are desired, the `use_<provider>_yml()` functions can be used.
Refer to [the complete list of options](https://docs.ropensci.org/tic/reference/yaml_templates.html).

For an R package, the following steps will be set up for the CI workflow:

- Installation of required dependencies for the project (dependencies are scraped from the DESCRIPTION file\*)
- Satisfying build-time dependencies of steps to be run in all CI stages (by scraping `pkg::fun` calls in `tic.R`)
- Checking of package via `rcmdcheck::rcmdcheck()`
- Creation of a `pkgdown` site including Github deployment
- Running a code coverage and upload to [codecov.io](https://codecov.io/)

See the [Getting Started](https://docs.ropensci.org/tic/articles/tic.html) vignette for more information and links to [minimal example repositories](https://docs.ropensci.org/tic/articles/tic.html#examples-projects) for various R projects (package, blogdown, bookdown and more).

#### Quickstart

If you are a new user, run

```r
tic::use_tic()
```

If you already use {tic} and want to configure a new CI provider, do

```r
### Circle CI ------------------------------------------------------------------

circle::use_circle_deploy() # (optional for deployment)
tic::use_circle_yml() # optional: Change `type` arg to your liking
tic::use_tic_r("package", deploy_on = "circle")
# (all of the above in one call)
# tic::use_tic(wizard = FALSE, linux = "circle", mac = "none", windows = "none",
#              matrix = "circle", deploy = "circle")
tic::use_update_tic()

### GitHub Actions -------------------------------------------------------------

tic::use_ghactions_deploy() # (optional for deployment)
tic::use_ghactions_yml() # optional: Change `type` arg to your liking
tic::use_tic_r("package", deploy_on = "ghactions")
# (all of the above in one call)
# tic::use_tic(wizard = FALSE, linux = "ghactions", mac = "ghactions",
#              windows = "ghactions", matrix = "ghactions", deploy = "ghactions")

tic::use_tic_badge("ghactions")
tic::use_update_tic()
```

## Good to know

We would like to mention that {tic} is a choice and sits on top of existing community efforts providing R support for various CI providers.
While {tic} will prevent you from dealing/learning every CIs YAML syntax, you will have to learn {tic}'s way of specifying your tasks on CI systems.

Also, there is no way around familiarizing yourself with the basics of CI systems in general.
Without this knowledge, you will also have a hard way understanding {tic}.

We also recommend to take a look at the projects providing the direct R support for each CI system (which {tic} builds upon) to gain a deeper understanding of the whole concept.

## Updating

Updating of YAML templates is supported via [`update_yml()`](https://docs.ropensci.org/tic/reference/update_yml.html).
See vignette ["Updating Templates"](https://docs.ropensci.org/tic/articles/updating.html) for more information.

## Vignettes

- [Get started](https://docs.ropensci.org/tic/articles/tic.html)

- [Feature Overview](https://docs.ropensci.org/tic/articles/features.html)

- [The CI Build Lifecycle](https://docs.ropensci.org/tic/articles/build-lifecycle.html)

- [CI Providers](https://docs.ropensci.org/tic/articles/ci-providers.html)

- [CI Client Packages](https://docs.ropensci.org/tic/articles/ci-client-packages.html)

- [Advanced Usage](https://docs.ropensci.org/tic/articles/advanced.html)

- [Deployment](https://docs.ropensci.org/tic/articles/deployment.html)

- [Custom Steps](https://docs.ropensci.org/tic/articles/custom-steps.html)

- [FAQ](https://docs.ropensci.org/tic/articles/faq.html)

## Limitations

The setup functions in this package assume Git as version control system, and GitHub as platform.
Automated setup works best if the project under test is located in the root of the Git repository.
Multi-project repositories are not supported, see [the comment by @jwijffels](https://github.com/ropensci/tic/issues/117#issuecomment-460814990) for guidance to work around this limitation.

The DESCRIPTION files needs to live in the project root.
 To simplify its creation have a look at [usethis::use_package()](https://usethis.r-lib.org/reference/use_package.html) or [usethis::use_description()](https://usethis.r-lib.org/reference/use_description.html).
<!-- NEWS.md is maintained by https://cynkra.github.io/fledge, do not edit -->

# tic 0.11.4.9000 (2022-01-28)

- Same as previous version.


# tic 0.11.4 (2022-01-28)

- Add compatibility for rlang >= v1.0.0


# tic 0.11.3 (2021-12-22)

- `do_pkgdown()`: Fix accidental deployments from pull requests (#308)


# tic 0.11.2 (2021-12-05)

- `do_pkgdown()` now always create a `.nojekyll` file for both release and development deployments.
  Otherwise custom fonts starting with an underscore will not be loaded as Jekyll ignores this pattern.
  A `.nojekyll` file tells GitHub pages to not use Jekyll for serving the web page. (#307)


# tic 0.11.1 (2021-06-27)

- Templates: install required system libs for {pkgdown} conditionally (accidentally removed in the previous template revision on 2021-06-26)
- Templates: restore installation of `libcurl4-openssl-dev` and `libgit2-dev`


# tic 0.11.0 (2021-06-26)

- Templates: On Linux, system libraries are now installed via `remotes::system_requirements()` (#300)
- Instead of using an exact version tag, the core GHA actions are now referenced using a dynamic major version tag (e.g. v2 instead of v2.3.4).
  This includes an update of the templates to the latest revision date 2021-06-26.
- `do_pkgdown()` macro now also builds the site on on branches containing the word `cran`.
  This adds support for the {fledge} release mechanism when using both a development and release site (#303)
- `update-tic.yml`: Remove hardcoded reference to master branch
- `update_yaml()` is not in beta state anymore


# tic 0.10.0 (2020-12-11)

- Drop Travis support (#295)
- Drop Appveyor support (#296)
- Bump templates: install required `libgit2` required by usethis v2.0.0 (tic dep)


# tic 0.9.0.9008 (2020-11-18)

- update peter-evans/create-pull-request action in `update-tic.yml` template to v3.5.0


# tic 0.9.0.9007 (2020-11-14)

- update GitHub Actions templates
  - update actions/checkout to v2.3.4
  - update actions/upload-artifacts to v2.2.1
  - update pat-s/always-upload-cache to v2.1.3
- conditionally install pkgdown required system libs on both Linux and macOS

# tic 0.9.0.9006 (2020-09-19)

- Replace hardcoded references to "master" by a dynamic query of the default repo branch


# tic 0.9.0.9005 (2020-09-04)

- `use_tic()`: use GitHub Actions as the default provider for all platforms
- Copy over GitHub authentication and SSH helpers from {travis}


# tic 0.9.0.9004 (2020-08-27)

- `update_yml()`: Support updating multiple YAML files


# tic 0.9.0.9003 (2020-08-06)

- GHA: add `workflow_dispatch` event trigger to templates
- update instructions for spatial libs on macOS for GHA
- improve heuristic for updating header parts of "custom" and "custom-deploy" templates


# tic 0.9.0.9002

- DSL: Don't add steps twice, if present in a previous macro (#272)
- `update-tic.yml`: use peter-evans/create-pull-request@v3 and actions/checkout v2.3.1. Run on ubuntu instead of macOS
- run r-devel on ubuntu instead of macOS
- pin actions/upload-artifact to v2.1.1
- update actions/checkout to v2.3.1
- update pat-s/always-upload-cache to v2.1.0
- `step_setup_ssh()` now verifies that {git2r} is installed. This prevents build failures for {rsconnect} deployments
- `update-tic.yml`: install libs via `apt` on Linux instead of `brew`


# tic 0.9.0.9001

- gha_add_secret(): Add new upstream parameters and fix endpoint


# tic 0.9.0.9000

- Same as previous version.


# tic 0.9.0

## Features

- `update_yaml()`: Account for duplicated env vars when a custom env var masks a template env var
- `use_tic_badge()`: Update tic badge and default action name (#269)
- Installing and using `ccache` for faster source package installation is now optional.
  While using `ccache` can help a lot for installing large dependency chains of certain packages, it also adds substantial overhead to builds for small packages.
  It is now optional and needs to be added as a custom block to builds. (#264)
- Add `step_session_info()`.
  This step prints the session info after having installed all dependencies in the "install" stage. (#259)
- `step_install_deps()` and `do_package_checks()` gain `dependencies = TRUE` argument.
- New `use_update_tic()`: Adds GitHub Actions workflow `update-tic.yml` to automatically update tic YAML templates
- Support fully custom runner matrices on GitHub Actions via template types `"custom"` and ´"custom-deploy"`
- New `gha_add_secret()` to automate the process of adding a GitHub PAT to a repo as a secret.
  This function will probably be move to {ghactions} in the future.

## Bugfixes

- Temporarily enforce {covr} dev version to account for timeouts on GHA, see https://github.com/r-lib/covr/issues/435
- Remove alert in steps-install.R (#263)
- Pass arg `remote` to all printing instances. Previously using a different remote than "origin" caused an error..

## CI Provider specific

### GitHub Actions

- Update versions of "tinytex" and "always-upload-cache" actions (#267)
- Install LaTeX on only one runner (#257)
- Switch from `main.yml` to `tic.yml` (#260)
- Set env var GITHUB_PAT from secret GITHUB_TOKEN to work around rate limits in {remotes}
- Update `actions/checkout` to v2.1.1
- Update `pat-s/always-upload-cache` to v2.0.0
- Remove old clang7 compiler setup for R <= 3.6.3

# tic 0.8.0.9009

- Temporarily enforce {covr} dev version to account for timeouts on GHA, see https://github.com/r-lib/covr/issues/435
- `use_tic_badge()`: Update tic badge and default action name (#269)
- GHA: Update versions of "tinytex" and "always-upload-cache" actions (#267)


# tic 0.8.0.9008

- Make ccache optional (and more) (#264)
- Remove alert in steps-install.R (#263)


# tic 0.8.0.9007

- Add `step_session_info()` (#259)
- GHA: Install LaTeX on only one runner (#257)
- GHA: Switch from main.yml to tic.yml (#260)


# tic 0.8.0.9006

- `step_install_deps()` and `do_package_checks()` gain `dependencies = TRUE` argument.


# tic 0.8.0.9005

- GHA: Set env var GITHUB_PAT from secret GITHUB_TOKEN to work around rate limits in {remotes}


# tic 0.8.0.9004

- New `use_update_tic()`: Adds GitHub Actions workflow `update-tic.yml` to automatically update tic YAML templates


# tic 0.8.0.9003

- Support fully custom runner matrices on GitHub Actions via template types `"custom"` and ´"custom-deploy"`
- bugfix: Pass arg `remote` to all printing instances. Previously using a different remote than "origin" errored.


# tic 0.8.0.9002

- New `gha_add_secret()` to automate the process of adding a GitHub PAT to a repo as a secret.
  This function will probably be move to {ghactions} in the future.


# tic 0.8.0.9001

### GitHub Actions

- Update actions/checkout to v2.1.1
- Update pat-s/always-upload-cache to v1.2.0
- Remove old clang7 compiler setup for R <= 3.6.3


# tic 0.8.0.9000

- Same as previous version.


# tic 0.8.0

## Features

- New `update_yml()`: Update your {tic} yaml templates to the latest upstream version in {tic}.
  User changes are preserved if these are marked correctly.
  See vignette ["Updating Templates"](https://docs.ropensci.org/tic/articles/updating.html) for instructions.
  This process can also be fully automated via a [custom CI job](https://docs.ropensci.org/tic/articles/updating.html#automating-the-update-process).
- Add argument `force` to `step_do_push_deploy()` for adding the `--force` flag to git calls
- Add solutions to {rgl} installation issues to FAQ
- Update `.R/Makevars`

## CI Provider specific

### GitHub Actions

- Set CRON time to 4 am to avoid download issues with mirror updates
- Added `-I/usr/local/include` to CPPFLAGS for macOS runners to mirror CRAN setup

### Circle CI

- Update r-oldrelease to R 3.6.3

# tic 0.7.0.9000

- GHA: added `-I/usr/local/include` to CPPFLAGS for macOS runners to mirror CRAN
- Add solutions to installation issues for package {rgl} to FAQ
- Add argument `force` to `step_do_push_deploy()` for adding the `--force` flag to git calls

# tic 0.7.0

## Macros

- Add `do_blogdown()` macro (#242)

## CI Provider specific

### GitHub Actions

- `use_tic()`: Move `cli::tree()` calls to `use_*_yml()` functions to avoid printing of false-positive trees.
- `use_*_yml()`: Set defaults for argument `type`.
- Fix GHA build URL and prettify deploy message (#247)
- Adjust GH Actions templates to use the `use_*_yml()` logic (#246)
- Bugfix: Packages on R-devel macOS are being installed in parallel again.

**R 4.0 toolchain**

- GitHub Actions: R-devel on macOS now uses Apples default clang compiler and the 10.13 SDK (High Sierra) to mimic the CRAN toolchain behavior.
  (The 10.15 SKD causes various issues when installing packages from source.)
- Env var `SDKROOT` is now set to `/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk` to prevent linking issues on macOS >= 10.15.4

# tic 0.6.0.9002

- `do_blogdown()` and `do_bookdown()` gain argument `cname`, making it possible to pass a CNAME URL for deployments. This is useful when setting one sets `orphan = TRUE` and relies on a custom URL of the published content (otherwise the redirect would not work)
- Add a better general intro about CI and explain some general CI terms (fixes #234)

# tic 0.6.0.9001

- Add `do_blogdown()` macro (#242)

# tic 0.6.0.9000

- Same as previous version.

# tic 0.6.0

## General

- `use_badge`: Refactor to use default badges from the respective providers rather than from shields.io (too slow and sometimes badges did not render at all) (#240)
- Condition deployment templates on a single runner for deployment. This avoids race conditions during deployment. This applies to all CI providers and templates (blogdown, bookdown, package) (#241)
- Files specified for deployment via `step_push_deploy(commit_paths = )` are now force added to the index by `git`.
  This enables to add directories like `docs/` (e.g. created by a local pkgdown build) to `.gitignore` and still deploy it during CI (#237).

## CI Provider specific

### GitHub Actions

- GitHub Actions: Always use option '--no-manual' on Windows because LaTeX is not available (because it takes ages to install)
- `step_rcmdcheck()`: Test in directory "check" to simplify upload of artifacts
- Set cron job to 4am to avoid potential download issues with R-devel on macOS
- Github Actions: Only deploy on R-release on macOS by default.

## Bugfixes

- `use_tic()` fails with descriptive error message if the badges start/end sections are missing in README
- `step_install_ssh_keys()`: Do not use `git2r::config()` when deploying on Windows to prevent build freezes

## Documentation

- `faq.Rmd`: Add info how to avoid git race conditions during pkgdown deployment (#238)

# tic 0.5.0.9005

- `use_tic()` fails with descriptive error message if the badges start/end sections are missing in README

# tic 0.5.0.9004

- `faq.Rmd`: Add info how to avoid git race conditions during pkgdown deployment (#238)
- `step_install_ssh_keys()`: Do not use `git2r::config()` when deploying on Windows to prevent build freezes
- update blogdown templates
- GitHub Actions: Always use option '--no-manual' on Windows because LaTeX is not available

# tic 0.5.0.9003

- Files specified for deployment via `step_push_deploy(commit_paths = )` are now force added to the index by `git`.
  This enables to add directories like `docs/` (e.g. created by a local pkgdown build) to `.gitignore` and still deploy it during CI (#237).
- `step_rcmdcheck()`: Test in dir "check" to simplify upload of artifacts

## Github Actions

- Set cron job to 4am to avoid potential download issues with R-devel on macOS
- Github Actions: Use actions/checkout v2
- Github Actions: Only deploy on R-release on macOS by default.
  This avoids git race conditions between runners.

# tic 0.5.0.9002

- Github Actions: {covr} now supports automatic upload of codecov results via their own CODECOV_TOKEN
- `use_tic_r()`: Add support for conditional tic.R templates via argument `deploy_on`.
- export `use_tic_r()` so that a manual workflow is possible (besides `use_tic()`)
- GitHub Actions: use actions "pat-s/always-upload-cache" instead of "actions/cache"

# tic 0.5.0.9001

- Add `use_tic_badge()`: Creation of pretty CI status badges

# tic 0.5.0.9000

- Same as previous version.

# tic 0.5.0

## Enhancements

- New function `tic::use_ghactions_deploy()` (status "experimental") to set up a SSH key for deployment.
- New function `use_ghactions_yml()` with `deploy = TRUE/FALSE` (FALSE by default).
- New Vignette "FAQ".
- Added GH Actions support to `use_tic()`
- new macro `do_readme_rmd()` (#223)
- new function `list_macros()`

## Maintenance

- Change for the default of the private SSH deploy key name from `TRAVIS_DEPLOY_KEY` to `TIC_DEPLOY_KEY` to have a generic name.
- Change argument `travis_private_key_name` to `private_key_name`
- Renamed `yaml-templates.R` to `yaml_templates.R` because the former caused troubles when previewing the dev version of the docs.
- Beautified the CLI output of `use_tic()`
- Replaced all instances of `_tic_` in the docs by `{tic}`

# tic 0.4.0.9000

- add macro `do_readme_rmd()` (#223)
- new function `list_macros()`

# tic 0.4.0

- add macro `do_drat()`
- start vignette "troubleshooting"
- add {desc} to suggests
- `ci_can_push()` never fails.
- templates: always upgrade dep packages during {tic} installation

# tic 0.3.0.9005

- Make it possible to pass the endpoint arg from {travis} funs to `use_tic()`
- mention the difference between .com and .org -> new vignette "org-vs-com"
- move package to ropensci org

* `error_on = "note"` also fails on warnings.

# tic 0.3.0.9004

- `ci_can_Push()`: Error with descriptive error message if deployment is not possible
- `ci_can_push()`: Fix for Travis CI
- optimize templates (especially matrix builds) by specifying which job is used for the pkgdown build

# tic 0.3.0.9003

- `use_tic()`: add key_name_private and key_name_public args from `travis::use_travis_deploy()`
- `ci_can_push()`: Change default from `"id_rsa""` to `"TRAVIS_DEPLOY_KEY"` and also support backward comp
- `use_tic()`: Travis as default for Linux and macOS

# tic 0.3.0.9002

- `use_tic()` supports running both Linux and macOS on Travis (#202).
- Skip `TicStep$prepare` if `prepare_call` is given in `add_code_step()` (#211).
- Fix preparation of `step_add_to_drat()`.
- `use_tic()` gains arguments that allow non-interactive use and re-running with the same settings if setup fails (#203).
- Removed artificial sleeps with interactive setup.

# tic 0.3.0.9001

- Move `use_travis_deploy()` back to {travis}.

# tic 0.3.0

- add argument "check_dir" to step_rcmdcheck (#179)
- use `remotes::install_cran(upgrade = TRUE)` to install packages (#186)
- added support for Circle CI (#177)
- All packages installed for custom steps use binary packages if possible on Windows and macOS (#178).
- Use `TRAVIS_BUILD_WEB_URL` for the commit message.
- `do_package_checks()` gains `type` argument.
- Tweak documentation.
- export `use_travis_yml()`, `use_circle_yml()` and `use_appveyor_yml()` and add overview table of available options

# tic 0.2.13.9020

- Avoid building packages when installing dependencies.
- Remove vignettes from package if checking with `--no-build-vignettes` to avoid warning from `R CMD check`.
- Fix `R CMD build` and `R CMD check` switches on AppVeyor.

# tic 0.2.13.9019

- Building pkgdown site succeeds if `docs/` directory is missing (#173, r-lib/pkgdown#1050).

# tic 0.2.13.9018

- Move `use_travis_deploy()` from the travis package to here.
- Unexport `get_public_key()` and `encode_private_key()`.

# tic 0.2.13.9017

- Test utils and printing.
- Exclude code that can only run interactively or in a CI from coverage.
- Add comment regarding integration test.
- Strip long source code lines.
- Add review badge.
- Add `tic.R` to `.Rbuildignore` for internal tests.
- Update wordlist.
- Fix typos.
- Condition example on presence of Git repository.

# tic 0.2.13.9016

- Fix compatibility with git 2.21 and above for race conditions (#160).
- `step_build_pkgdown()` clean site before building.
- AppVeyor template makes sure packages are always installed from binary during bootstrapping.
- CI templates install from GitHub if the version number indicates that the package is not on CRAN yet.
- AppVeyor doesn't cache R packages, because this leads to update problems. Binary installation is fast enough.
- Don't perform CRAN incoming checks, in particular the checks for large version components (#168).
- The `step_install_deps()`, `step_install_cran()` and `step_install_github()` steps install binary packages by default, even if the CRAN version is ahead.
- All files created by `use_tic()` are added to `.Rbuildignore`.
- Package template for `tic.R` runs pkgdown only on Travis (#167).
- Update vignettes (#156).

# tic 0.2.13.9015

- `detect_repo_type()` now prompts the user for unknown repository types (#161).
- `use_tic()` loses `path` argument, now taken from `usethis::proj_get()` .
- `step_rcmdcheck()` and `do_package_checks()` now avoid building the vignette by default on AppVeyor (#150).
- `use_tic()` now uses boxes from {cli} for better structured output (#153).

# tic 0.2.13.9014

- Configuration storage modeled after `usethis::proj_get()`.
- New `dsl_load()`, renamed from `load_from_file()`.
- New `dsl_get()` and `dsl_init()`.
- Added examples to help for `get_stage()` and macros (#77).

# tic 0.2.13.9013

- Using tidy evaluation for simpler code, more control and better printing of steps (#77).
- Fix AppVeyor builds.
- The README is now explicit about suggesting that each repo should contain only one project (#152).
- Documentation uses the {rotemplate} package (#121).
- Only install {remotes} and {curl} if not yet installed (#97).
- New `use_tic()`, moved from {travis} (#138).
- Updated templates (#81).
- A failing step displays a traceback generated by `rlang::trace_back()` (#105).
- `do_pkgdown()` and `do_bookdown()` now have a `deploy` argument and are documented on separate help pages. The new `?macro` help page provides an overview.
- Implement `print()` methods for DSL and stages (#77).
- New `do_bookdown()` (#137).

# tic 0.2.13.9012

- New `repo_*()` functions to simplify specification of the `repos` argument to installer functions (#101).
- Add Appveyor checks (#147, @pat-s).
- New pkgdown macro via `do_pkgdown()` (#126, @pat-s)
- New example: covrpage, cc @yonicd
- `step_rcmdcheck(error_on = "note")` works again (#119).
- New `do_package_checks()` with `codecov = TRUE` argument (#146), replaces `add_package_checks()` which stays around for compatibility (#128).
- `add_step()` now evaluates the `step` argument in a `tryCatch()` block and gives a bit of context if this fails (#73).
- New `run_all_stages()`, previously `tic()` (#66).
- New `ci_get_env()`, `ci_has_env()` and `ci_is_env()` functions to avoid verbose `Sys.getenv()` calls in `tic.R` (#124, @pat-s).
- New `ci_*()` functions to avoid R6 notation in `tic.R` (#125, @pat-s).

# tic 0.2.13.9011

- New `step_install_deps()`, reorganizing help pages so that installer steps are on the same page.
- `step_rcmdcheck()` no longer installs dependencies. Instead, `add_package_checks()` includes `step_install_deps()` (#74).
- Fix two links in README (#115, @Rekyt).
- Vignette update (#80, @pat-s).
- Support `build_args` argument in `step_rcmdcheck()` (#64, @pat-s).

# tic 0.2.13.9011

## step_rcmdcheck()

- deprecate `warnings_are_errors` and `notes_are_errors` in favor of the new `error_on` argument
- add args `timeout` and `repos`
- call `rcmdcheck()` internally with `error_on = "never"` so that we can trigger the message on found warnings and notes
- remote outdated doc about `step_rcmdcheck()` using a dedicated lib for the check

# tic 0.2.13.9010

- No longer using a separate library for package checks, because it causes a lot of problems with various steps which are not aware of this (#86, #88).
- Packages coming with the R-installation are not updated anymore when preparing `step_rcmdcheck()`.
  See `?step_rcmdcheck()` for detailed info. (#103)

# tic 0.2.13.9009

- The `step_build_pkgdown()` step now uses the same dedicated library as `step_rcmdcheck()`.
- Using the development version of _rcmdcheck_ to work around problems finding the vignette builder (#84).
- Draft for new "Get started" vignette (#63, @pat-s).

# tic 0.2.13.9008

- The `step_rcmdcheck()` step now uses a dedicated library for installing the packages and checking, it also updates the packages after installing dependencies. The `add_package_checks()` macro no longer includes an `update.packages()` call (#35).
- The `step_rcmdcheck()` step now installs all dependencies during preparation. The `add_package_checks()` macro no longer adds the code step that installs dependencies.

# tic 0.2.13.9007

- The `step_do_push_deploy()` and `step_push_deploy()` steps are not executed for builds on a tag, because this would create a branch of the same name as the tag (#27).

# tic 0.2.13.9006

- Support creating variables in `tic.R` by sourcing `tic.R` in a modifiable environment (#33).
- Replaced `private` arguments with an environment that keeps track of internal state, now the code from `add_package_checks()` can be copied to a `tic.R` file (#74).

# tic 0.2.13.9005

- A failing step immediately fails the entire stage, subsequent steps are not run (#59).

# tic 0.2.13.9004

- New `get_public_key()` and `encode_private_key()` moved from _travis_ (#71, @pat-s).
- Add `step_install_cran()` and `step_install_github()` (#65, @pat-s).

# tic 0.2.13.9003

- Added integration tests for package checks and deployment, covering various common cases (#62).
- Add integration test for deploying from a subdirectory.
- Remove `orphan` argument from `step_push_deploy()`, because there's no easy way to implement it reliably. If only a subdirectory is deployed to a separate branch (i.e., the `path` argument is set), `orphan = TRUE` is required.

# tic 0.2.13.9002

- Better strategy for handling race conditions during deployment, new changes are no longer silently overwritten with `step_push_deploy()` (#45).
- Add integration test for package checks and race conditions (#62).
- Clarify error message upon step failure.
- `add_package_checks()` adds coverage checks only for non-interactive CIs.
- Add reference to `use_tic()` (#55).
- Document purpose of testing steps (#49).
- Allow only predefined stage names (#48).

# tic 0.2.13.9001

- The _pkgdown_ package is installed from CRAN.

# tic 0.2.13.9000

- New `subdir` argument to `step_push_deploy()` and `step_do_push_deploy()`, allows restricting the set of files to be committed to Git (#42).

# tic 0.2-13 (2018-02-01)

- New `base64serialize()` and `base64unserialize()` (#37).
- `add_code_step()` detects required packages also for complex expressions. Packages that need to be installed from GitHub still need to be installed manually (#36).
- `step_rcmdcheck()` now prints a summary, which also shows e.g. details on installation failures.
- New `prepare_call` argument to `step_run_code()` and `add_code_step()`.

# tic 0.2-12 (2017-06-29)

- Fix `add_package_checks()`.

# tic 0.2-11 (2017-06-29)

- `add_package_checks()` gains arguments that are passed to `step_rcmdcheck()`.
- New `step_setup_ssh()` (#24).
- New `add_code_step()` (#21).
- New `tic()` to run all steps locally (#23).
- New `add_package_checks()` (#25).

# tic 0.2-10 (2017-06-29)

- Document all exported functions and many classes (#8).
- `step_add_to_drat()` will also update the overview page if it exists.

# tic 0.2-9 (2017-06-28)

- Fix `get_slug()` on AppVeyor to use `APPVEYOR_REPO_NAME`.
- New `step_add_to_drat()`.
- Split `step_push_deploy()` into `step_setup_push_deploy()` and `step_do_push_deploy()`.
- Better traceback output.
- Use "remotes" instead of "devtools".
- Reduce output after preparation (#5).
- New `step_rcmdcheck()`.
- The deparsed code is used as step name (#5).

# tic 0.2-8 (2017-06-17)

- An error occurring when running a step is printed in red (#5).

# tic 0.2-7 (2017-06-13)

- New `step_write_text_file()` for creating arbitrary text files, including `~/.R/Makevars` (#14).
- pkgdown documentation is now built for tags by default (#13).
- The "openssl" package is now only suggested, not imported.
- Removed `step_run_covr()` in favor of the new `step_run_code()` (#18).
- `load_from_file()` reloads the file from disk if its mtime changes (#11).
- All steps of a stage are run even in case of previous errors, but the stage still fails if at least one of its steps failed (#10).
- Adding to known hosts or installing a SSH keys now requires a non-interactive CI.
- New `step_run_code()` to run arbitrary code. If the code is a call with the `pkg::fun()`, notation, pkg is installed if missing (#1, #3). `step_run_covr()` remains for compatibility but is scheduled for removal.
- Color the start of each step in the log (#5).
- New `step_add_to_known_hosts()` to work around configuration problems on OS X (#16).
- Export runner methods for all stages defined in Travis CI and AppVeyor (#17).

# tic 0.2-6 (2017-06-04)

- Technical release to sync default and production branches.

# tic 0.2-5 (2016-11-27)

- Fix `after_success()` and `deploy()`.
- Step names are now printed again.

# tic 0.2-4 (2016-11-27)

- Use new DSL with the notion of stages with arbitrary names.
  - New `load_from_file()` replaces `get_xxx_steps()`
  - `task_...()` has been renamed to `step_...()`
  - A task is now something like an ad-hoc step
  - `before_script()` is now `prepare_all_stages()`
  - `TravisTask` is now `TicStep`
  - `ci()` is now exported
- If environment variable `CI` is undefined, use `LocalCI` with sensible inference of repository and branch.
- Stop if `git` exits with nonzero status.

# tic 0.2-3 (2016-11-06)

- Install package for `task_build_pkgdown` task.

# tic 0.2-2 (2016-11-05)

- DSL to define steps via `step()`, which are tasks with a branch and/or env var filter (#6).

# tic 0.2-1 (2016-11-05)

- Support environment variables from both Travis and AppVeyor (#6).
- Add tests.
- Rudimentary support for multiple CI systems.
- Clean up dependencies.

# tic 0.2 (2016-11-05)

Initial release.

- Rudimentary configuration based on task objects. A task object is a list/environment which contains at least the members `check()`, `prepare()` and `run()` -- functions without arguments, only `check()` needs to return a `logical` scalar. These can be subclasses of the new `TravisTask` R6 class, the package now contains six subclasses: `HelloWorld`, `RunCovr`, `BuildPkgdown`, `InstallSSHKeys`, `TestSSH`, and `PushDeploy`. The `new` methods of theses subclasses are exported as `task_hello_world()`, `task_run_covr()`, `task_build_pkgdown()` `task_install_ssh_keys()`, `task_test_ssh()`, and `task_push_deploy()`, respectively. The three functions `before_script()`, `after_success()` and `deploy()` accept a semicolon-separated list of task objects, which is by default taken from the `TIC_AFTER_SUCCESS_TASKS` and `TIC_DEPLOY_TASKS` environment variables. These functions call the `prepare()` and `run()` methods of the task objects if and only if the `check()` method returns `TRUE` (#42).
# How to open/modify `.graphml` files

This filetype is specfic for the application [yEd](https://www.yworks.com/yed).
You can install the application locally or use the [live version](https://www.yworks.com/products/yed-live).
---
title: "Build lifecycle"
author: "Patrick Schratz, Kirill Müller"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Build lifecycle}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Stages

CI services run builds in stages.
Stages are usually ordered as follows:

```{r, echo = FALSE, fig.align='center', dpi = 50}
knitr::include_graphics("img/build-lifecycle.png")
```

The `after_xxx` stages are executed conditionally after their corresponding `xxx` stage.

- The `after_deploy` stage will only be run if the `deploy` stage was run before.
- The `after_success` stage will only be run if the `script` stage executed successfully, i.e. without error; otherwise `after_failure` will be run instead.

*tic* also relies on a "stage" based approach.
All action that should be run in a certain stage are defined in `tic.R`.
The steps are specified in an CI-agnostic way using R syntax.

The majority of the template code consists of glue code and is not meant to be edited.

In a nutshell, the workflow is as follows:

CI YAML -> `tic.R` -> commands/steps to execute

Some important points:

- The R code declared in `tic.R` is not meant to be run manually.
  It also does not trigger a CI build.
All commands just define the workflow of the CI build.
- The workflow can be loaded using `dsl_load()`, however this will not run any of the commands defined.
- For testing purposes, all stages and steps defined in `tic.R` can be executed by calling `run_all_stages()`.
  This emulates a CI build on your local system.
  See [Troubleshooting: Running tic locally](advanced#troubleshooting-running-tic-locally) for more information.

### Accessing a single stage

The steps which are executed in each stage are specified in `tic.R`.
A stage is executed by calling the respective *tic* function; for example for stage "deploy" `tic::deploy()`.
This is what happens in the CI YAML templates if you take a closer look at them.

These functions then source `tic.R` and collect all steps which belong to their stage by executing `get_stage("<stage name>")` (e.g. `get_stage("deploy")` for the "deploy" stage").

Again, remember that the order of the stages is fixed (see the ["Stages"](#stages) section), it does not matter in which order you declare the stages in `tic.R`.

### Details of stages

#### The `"before_install"` & `"install"` stages

An important stage for {tic} is the `"before_install"` stage.
Here, {tic} itself gets installed and runs `prepare_all_stages()`.
This function ensures that all subsequent steps can be executed.
Under the hood the `prepare()` method of all steps that were declared in `tic.R` is being called.
For example, the `prepare()` method of the `step_rcmdcheck()` step ensures that all dependencies of an R package get installed by calling `remotes::install_deps()`.

All packages that should be stored in the "cache" of the CI service (so that they do not need to be installed again on every CI build) should be installed during preparation.

#### The `"script"` stage

The `"script"` stage is responsible for executing the important tasks of the CI run: Typically, it runs `R CMD check` for a package or builds the site for a blogdown site.
When arriving at this stage, all dependencies for a successful run are already installed.

#### The `"deploy"` stage

This stage initiates the deployment (e.g., setting up deployment keys) and executes it.
If you want to automatically build a {pkgdown} site, you can do it here.
See [the article about deployment](deployment.html) for more information.

## Steps

Steps are the commands that are executed in each stage.
*tic* uses the [pipe operator](https://magrittr.tidyverse.org/) and the `add_step()` function to chain steps in `tic.R`, for example

```{r eval = FALSE}
get_stage("deploy") %>%
  add_step(step_build_pkgdown())
```

In the code example above `step_build_pkgdown()` is added to the `"deploy"` stage and subsequently only run in this stage.
More steps that should be run in this stage could just by piped after `add_step(step_build_pkgdown())`.
In summary, steps are usually defined using two nested commands: `add_step()` and the corresponding step, here `step_build_pkgdown()`.

Here is a list that shows a rough grouping of the steps into their default stages:

#### Basic

| Step | Description |
| ---  | ----------- |
| `step_hello_world()` | Print "Hello, World!" to the console, helps testing a tic setup| | `step_run_code()` | Run arbitrary code, optionally run preparatory code and install dependent packages. `add_step(step_run_code())` an be abbreviated with `add_code_step()`
| `step_write_text_file()` | Creates a text file with arbitrary contents |

#### Installation

| Step | Description |
| ---  | ----------- |
| `step_install_cran()` | Installs one package from CRAN via `install.packages()` if it is not yet installed. |
| `step_install_github()` | Installs one or more packages from GitHub via `remotes::install_github()` |

#### R package specific

| Step | Description |
| ---  | ----------- |
| `step_build_pkgdown()` | Building package documentation via [pkgdown](https://github.com/r-lib/pkgdown) |
| `step_rcmdcheck()` | Run `R CMD check` via the {rcmdcheck} package |

#### Deployment

| Step                       | Description                                                                                                                                               |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `step_install_ssh_key()`   | Make available a private SSH key (which has been added before to your project by `use_tic()` or `tic::use_ghactions_deploy()`).                           |
| `step_test_ssh()`          | Test the SSH connection to GitHub, helps troubleshooting deploy problems.                                                                                 |
| `step_setup_ssh()`         | Adds to known hosts, installs private key, and tests the connection.                                                                                      |
| `step_setup_push_deploy()` | Clones a repo, initiates author information, and sets up remotes for a subsequent `step_do_push_deploy()`.                                                |
| `step_do_push_deploy()`    | Deploy to GitHub.                                                                                                                                         |
| `step_push_deploy()`       | Combines `step_setup_push_deploy()` and `step_do_push_deploy()`.                                                                                          |
---
title: "Implementation Details of CI Providers"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Implementation Details of CI Providers}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## GitHub Actions

{tic} supports running builds on GitHub Actions on all major platforms (Linux, macOS, Windows).
The upstream support for the R language was developed by Jim Hester in [r-lib/actions](https://github.com/r-lib/actions).
This repo also stores some usage [examples](https://github.com/r-lib/actions/tree/master/examples) which differ to the {tic} approach in the following points:

- {tic} makes use of `ccache` for compiler caching enabling faster source installation of packages.
  The `ccache` directory is cached and build once a week.
- {tic} installs packages from source on Linux by default and does not use package binaries.
- {tic} caches the complete R library and not only the direct packages dependencies (`actions` does this via `remotes::dev_package_deps(dependencies = TRUE)`).
  The cache is built once per day.

Making use of binaries can speed up build times substantially.
This can be especially attractive for packages with many dependencies or dependencies which take a long time to install.
However, binaries do oft run into problems when the package needs linking against system libraries.
The most prominent example for this is {rJava}.
If the binary was built with the same version as the user is running on the system, everything will work.
However, often enough a different version of the system library is installed and the R packages needs to be installed from source to successfully link against it.

For the case of {rJava}, one needs to

- add a call to `R CMD javareconf` for **macOS** runners
- add a call to `sudo R CMD javareconf` for **Linux** runners

### macOS toolchain

macOS is a bit tricky when it comes to source installation of packages.
By default `clang` is used instead of `gcc` (Linux) because the former is the default for macOS.
However, the default `clang` of macOS does not come with openMP support.
Therefore, the R macOS core devs and CRAN currently use a [custom openMP-enabled](https://cran.r-project.org/bin/macosx/tools/) (old) version of `clang` to build the CRAN package binaries.
In {tic} we reflect this by installing `clang7` and `clang8` for the respective R version during build initialization in the "ccache" stages.

### rJava

If Java support is required, add the following for macOS runners:

```yaml
      - name: "[macOS] rJava"
        if: runner.os == 'macOS'
        run: |
          R CMD javareconf
          Rscript -e "install.packages('rJava', type = 'source')"
```

For Linux, add `sudo R CMD javareconf` to stage "[Linux] Prepare".
We currently do not support Java on Windows.

### ccache

If you have a huge dependency chain and compiling many packages from source (especially on R-devel), `ccache` can help to speed up package installation.
It is recommended once your dependency installation time is higher than 30 minutes.

Once the `ccache` cache is build, compilation will complete much faster.
The `ccache` cache itself is only invalidated once a month.
This means package installation can make use of the cache in 29/30 days in a month.

The downside is that `ccache` needs to be installed and configured.
This happens in every run, i.e. also in runs in which `ccache` is not used because a package cache already exists.
Installation can take up to 1 min, depending on the platform.
Note that `ccache` won't be used on Windows since only binaries are used on this platform.

You can take the following blocks and add/replace them to your `tic.yml` as needed.
The essential part is to prefix the compiler settings in `~/.R/Makevars` with `ccache`.

```yml
      - name: "[Custom] [Cache] Prepare weekly timestamp for cache"
        if: runner.os != 'Windows'
        id: datew
        run: echo "::set-output name=datew::$(date '+%Y-%V')"


      - name: "[Custom] [Cache] Cache ccache"
        if: runner.os != 'Windows'
        uses: pat-s/always-upload-cache@v2.0.0
        with:
          path: ${{ env.CCACHE_DIR}}
          key: ${{ runner.os }}-r-${{ matrix.config.r }}-ccache-${{steps.datew.outputs.datew}}
          restore-keys: ${{ runner.os }}-r-${{ matrix.config.r }}-ccache-${{steps.datew.outputs.datew}}

      # install ccache and write config file
      # mirror the setup described in https://github.com/rmacoslib/r-macos-rtools
      - name: "[Custom] [macOS] ccache"
        if: runner.os == 'macOS' && matrix.config.r == 'devel'
        run: |
          brew install ccache
          # set compiler flags
          mkdir -p ~/.R && echo -e 'CC=ccache clang\nCPP=ccache clang\nCXX=ccache clang++\nCXX11=ccache clang++\nCXX14=ccache clang++\nCXX17=ccache clang++\nF77=ccache /usr/local/gfortran/bin/gfortran\nFC=ccache /usr/local/gfortran/bin/gfortran' > $HOME/.R/Makevars

      # install ccache and write config file
      - name: "[Custom] [Linux] ccache"
        if: runner.os == 'Linux'
        run: |
          sudo apt install ccache libcurl4-openssl-dev
          mkdir -p ~/.R && echo -e 'CC=ccache gcc -std=gnu99\nCXX=ccache g++\nFC=ccache gfortran\nF77=ccache gfortran' > $HOME/.R/Makevars
```

In addition you also need to set the following env variables:

```yml
# setting some ccache options
CCACHE_BASEDIR: ${{ GITHUB.WORKSPACE }}
CCACHE_DIR: ${{ GITHUB.WORKSPACE }}/.ccache
CCACHE_NOHASHDIR: true
CCACHE_SLOPPINESS: include_file_ctime
```

### Spatial libraries (gdal, proj, geos)

#### macOS

homebrew-core has formulas for `gdal`, `geos` and `proj`.
If you need more spatial formulas, have a look at the [osgeo4mac](https://github.com/OSGeo/homebrew-osgeo4mac) tap.
Note however, that when installing formulas from the latter, these will conflict with the ones from homebrew-core.
Either install all formulas from `osgeo4mac` or none.

Also one needs to remove the `gfortran` build that is installed with `actions/setup-r`.
This is due to `brew` installing `gcc` during the installation of `gdal`.
`gcc` comes with `gfortran` included and when `brew` tries to link `gfortran` it will fail since there is already a local instance of `gfortran`.
Hence, this instance needs to be removed so that the `brew link` step does not error and stop the build.

```yaml
# conflicts with gfortran from r-lib/actions when linking gcc
rm '/usr/local/bin/gfortran'
brew install gdal proj geos
```

#### Linux

On Linux, add `libgdal-dev libproj-dev libgeos-dev` to the `apt install` call in the "[Linux] Prepare" stage.

### Known issues

- [Windows] Installing {tinytex} for LaTeX availability does not complete

## Circle CI

WIP
---
title: "The features of tic"
author: "Patrick Schratz, Kirill Müller"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{The features of tic}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

What is the advantage of using {tic} compared to the default R setup?

1. Deployment to a Git repository is greatly simplified.

1. Support for R packages and other kinds of project (_bookdown_, _blogdown_, etc.), with predefined templates.
Set up your project to deploy rendered versions of your book or blog with a single push.

1. Workflow specification are specified in a single R script, regardless of the CI system used.
No need anymore for YAML configuration files that differ across CI services.

Other minor advantages include the use of `rcmdcheck::rcmdcheck()` for package checking (instead of `R CMD check`) and robust caching approach of project dependencies (via `ccache` and R package caching).

## Simplified Deployment

CI services can be used to automatically build and deploy files.
This means that these services can push certain files created during the build to repositories (GitHub, GitLab, etc.).
Possible use cases are:

- Changed site contents of a {pkgdown} site
- Updated .Rd files of a package (by calling `devtools::document()` before)
- Automated generation of a [test summary page](https://github.com/yonicd/covrpage) for a package

If {tic} should be used for deployment, some preparatory work is required:

- Setting up a SSH key pair for deployment (differs across CI services).
- Granting permission to push to the repo on GitHub.

When calling `use_tic()`, the streamlined preparation process is run, utilizing the [R Client packages](ci-client-packages.html) of the respective CI service under the hood.
This step is needed once per repository.

For more detailed information on deployment in {tic}, have a look [Deployment](deployment.html) vignette.

## Support for various R projects

{tic} was developed with the aim to not only simplify R package development using CI services but also to support other R project types:

- _bookdown_
- _blogdown_
- _drat_
- website deployment
- _figshare_ deployment

Each of these project types requires a standardized structure.
{tic} detects this structure (assuming the user already has set it up) and adds CI templates tailored towards this specific project type to the repository when calling `use_tic()`.
See the [Example projects](tic.html#example-projects) section in the "Get started" article for a list of supported project types including links to minimal example repositories.

## CI-Agnostic workflows

What does "CI-Agnostic" mean and why do we need it?

For to historic reasons, the R community first started on Travis CI to implement an easy way for R package checking.
The build script for R is [community maintained](https://github.com/travis-ci/travis-build/blob/3eddda591f544a071a62fc0f713183e128cfeac1/lib/travis/build/script/r.rb).
Theoretically, R could be run on any CI system.
Travis CI is only one out of a bunch of providers which offer (free) CI solutions.

Each CI provider has its own way how the user has to write the YAML file to successfully talk to the service.
This setup file controls what will be done in each run.

To give you an example how different these control files can be, take a look at these two examples from [Circle CI](https://github.com/ropensci/tic/blob/master/.circleci/config.yml) and [Github Actions](https://github.com/ropensci/tic/blob/master/.github/workflows/tic.yml).

We could list way more differences - but that's exactly the point when {tic} comes in!

- Rather than dealing with all the CI differences, {tic} enables the specification of the complete workflow in an external R script file `tic.R`.
- The calls listed in `tic.R` will work the same way on every CI service that is supported by {tic}.
- You can emulate all the stages and steps locally by calling `run_all_stages()`.
- You are independent to changes made to the upstream runners of a specific CI system.
- A `tic.R` workflow is usually shorter and easier to parse than a provider-specific YAML configuration file as it builds on [macros](#macros).

So instead of learning how to specify specific tasks on different CI platforms, you only apply R commands which work the same on all CI systems.

## Enhanced R package checking

As an R package developer calling `devtools::check()` is a common task.
Usually CI workers will run `R CMD check <package>` to check the R package.
{tic} instead makes use of [{rcmdcheck}](https://github.com/r-lib/rcmdcheck), a wrapper around `R CMD build` and `R CMD check` developed by [Gabor Csardi](https://github.com/gaborcsardi).
{rcmdcheck} comes with several enhancements:

- Coloring of important steps, simplifying the readability of the log.
- Enhanced and extended tracebacks of errors, especially errors in tests.
- The whole check process is returned in a R object, making it easier to inspect errors/warnings.

Especially the extended log printing of errors on the CI service is a huge advantage - often enough, CI services cut the console log output early, often omitting important information about the error.

## Caching of packages

When using {tic}, all dependencies (the ones of the package plus the ones of other stages of the CI build) are installed in the `"before_install"` and `"install"` stage.
This has the advantage that all packages are added to the cache (even if they are just needed for deployment), speeding up subsequent builds substantially.

More information about the complete workflow can be found in the [Build lifecyle](build-lifecycle.html) vignette.

## Easier troubleshooting

{tic} comes with the ability to [emulate a CI run locally](advanced.html#emulate-a-ci-run-locally) and [debug problems in the config file](advanced#troubleshooting-running-tic-locally) by calling `dsl_load()` locally.
---
title: "Updating"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Updating}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

{tic} comes with the convenience of providing YAML templates for various CI providers which remove the need of studying YAML notations of CI providers.
The provided templates aim to work for most packages that use plain R code.
Packages that rely on system libraries (such as GDAL, xml2) or link to other languages (such as Java) often need manual additions within the YAML file.

However, making manual changes to such prohibits users from simply updating to the latest upstream template that {tic} provides via `yaml_templates()` because custom user changes would be overwritten.

This often prevents users from updating to the most recent templates provided by {tic} since a manual comparison between the current template and the latest version is needed.
{tic} 0.8.0 tackles this issue: The new `update_yml()` enables updating of templates to the latest version with **user changes preserved**.

Currently this only works for GitHub Actions as this provider enables the inclusion of arbitrary custom blocks.
{tic} extracts these blocks, updates the template and then inserts the custom blocks at the correct position back into the updated template.
What may sound easy at first is in fact a very complicated task behind the scenes: indentation needs to preserved, some providers care about order and some do not accept duplicate keys.
While `update_yml()` supports GitHub Actions and Circle CI right now for custom environment variables and user blocks, the following rules need to be followed to ensure a smooth experience.

**Shared Rules**

- All templates need to have an *identifier* and a *revision date* in their first two lines.
- Comments of custom env vars may only span one line.
- All *custom env vars* need to have a comment (above the actual env var) which includes the term `[Custom env]`.

**GitHub Actions & Circle CI**

- All *custom blocks* need to include `[Custom block]` in their name tag, e.g. `- name: "[Custom block] Test custom user block2"`.
- Custom env vars unique to a runner need to include `[Custom matrix env var]` in their tag, e.g. `# [Custom matrix env var] test`.

All of these tags are used by {tic} to preserve user changes and choose the right upstream template.
During the process `update_yml()` also prints how many custom blocks and env vars were found.
If you have some in your template and nothing is printed when updating, something went wrong and you should double check your template.
In any case, it is recommended to review the changes to avoid unexpected CI failures.

## Deviating from the templates

### `custom` and `custom-deploy` templates

If you are using the `custom` or `custom-deploy` deploy template (e.g. via `tic::use_ghactions_yml("custom")`), `tic::update_yml()` will ignore the matrix part of the templates.
This gives you the ability to specify your own runner config while still profiting from template updates.

### The `# [Custom header]` tag

If you want to go even more custom, you can add `# [Custom header]` right below the `## revision date` line.
This tells `update_yml()` to ignore the complete header including the the `env:` key completely.
This can be useful if you want to insert a `service:` block between `env:` and `strategy` or when specifying custom build triggers in `on:`.

## Examples

Here are some examples for custom user blocks and custom env vars:

**GitHub Actions**

*Custom header*

```yml
## tic GitHub Actions template: linux
## revision date: 2021-06-27
# [Custom header]
on:
  workflow_dispatch:

[...]
```

*Runner specific environment variable*

```yml
matrix:
  config:
    # use a different tic template type if you do not want to build on all listed platforms
    - { os: windows-latest, r: "release" }
    # [Custom matrix env var] test
    - { os: macOS-latest, r: "release", pkgdown: "true", test: "true" }
    - { os: ubuntu-latest, r: "devel" }
    - { os: ubuntu-latest, r: "release" }
```

*Environment variable*

```yml
env:
  # [Custom env]
  R_MAX_NUM_DLLS: 999
```

*Block*

```yml
- name: "[Custom block] [macOS] xquartz"
  if: runner.os == 'macOS'
  run: |
    brew install xquartz
```

**Circle CI**

*Environment variable*

```yml
environment:
  # [Custom env] env var 2
  test2: "true"
```

*Block*

```yml
- run:
    name: "[Custom block] test2"
    command: |
      echo 'test2'
```

## Automating the update process

Updating {tic} YAML files can be automated further.
We provide a GitHub Actions Workflow named [update-tic.yml](https://github.com/ropensci/tic/blob/master/.github/workflows/update-tic.yml) which can be used together with `tic::update_yml()` to update the templates whenever there are newer upstream versions available.
The workflow will create a branch, update the files and commit them and even open a pull request.

Just put this workflow next to any `tic.yml` file within `.github/workflows/` and it will silently do its job.
By default it will run over night as a CRON job.
It only runs once a day and is not being executed on push or pull request events.
The underlying `update_yml()` will match all files starting with `"tic"`.
Hence you can add multiple YAML files with {tic} support, e.g. `"tic.yml"` and `"tic-db.yml"`.

Unfortunately, GitHub does not allow GHA workflow files to be updated and pushed by automatic approaches.
To make this work, the user needs to pass a GitHub Personal Access Token (PAT) with "workflow" scopes.
This PAT need to be added as a "secret" to the repo so that it can be used within the build.
`gha_add_secret()` helps to automate this process.
The linked workflow searches by default for a PAT secret named `TIC_UPDATE` when updating `tic.yml`.
---
title: "Advanced usage"
author: "Patrick Schratz, Kirill Müller"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Advanced usage}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

## Running steps conditionally

Steps and stages can be run conditionally using the control workflow of {tic}.
Possible conditionals are

- Environment variables of the build (queried by `ci_is_env()`, `ci_has_env()` or `ci_get_env()`).
- R Version of the current build (`getRversion()`).
- Other features of the build (e.g. branch name via `ci_get_branch()`).

Common tasks to use this feature are testing on multiple R versions and the restriction of certain tasks that should only be executed once (e.g. the [deployment of a {pkgdown} site](deployment.html#pkgdown-deployment)).

### Conditional execution: Use cases

The following shows some example code blocks to condition certain stages and their respective steps on

- the CI service

```{r eval = FALSE}
if (ci_on_ghactions()) {
  get_stage("<stage>") %>%
    add_step(step_...())
}
```

- a specific branch

```{r eval = FALSE}
if (ci_get_branch() == "main") {
  get_stage("<stage>") %>%
    add_step(step_...())
}
```

## Installation of packages

Required packages are installed based on the "Depends"/"Imports"/"Suggests" fields of the `DESCRIPTION` file.
You should only use the following steps in extraordinary circumstances.
An example can be the use of a package in the `README.Rmd` file which is not listed within the package's `DESCRIPTION` file.

### GitHub packages

```{r eval = FALSE}
get_stage("install") %>%
  add_step(step_install_github("r-lib/rcmdcheck"))
```

Note that the underlying `remotes::install_github()` is vectorized for the `repo` argument which means you can pass all packages you want to install in a single function call:

```{r eval = FALSE}
add_step(step_install_github(c("r-lib/rcmdcheck", "r-lib/usethis")))
```

### CRAN packages

```{r eval = FALSE}
get_stage("install") %>%
  add_step(step_install_cran("magick"))
```

## CI Meta-Information

The `ci()` function and its friends (`ci_*`) hold valuable information about the CI system.
They can be used to query information that can be again be utilized for conditioning stages or steps.

For example, the user may wish to only deploy on GitHub Actions by using `ci_on_ghactions()`:

```{r eval = FALSE}
if (ci_on_ghactions()) {
  get_stage("before_deploy") %>%
    add_step(step_setup_ssh())

  get_stage("deploy") %>%
    add_step(step_push_deploy())
}
```

By using the code above, the specified steps will only be executed on GitHub Actions.
See `?ci` for more information on which CI build information can be extracted from this function.

## Debugging: Running {tic} locally

### Checking for syntax errors

Before pushing to GitHub and triggering a CI build, `tic.R` can be validated using `dsl_load()`.
This function will source `tic.R` to check for possible problems.
If everything is ok, it will return (invisibly) a list with all stages that will be run in the CI build.
Here is a preview of the first two stages:

```{r eval = FALSE}
dsl_load()
```

```
✔ Loading tic stage configuration from tic.R
```

```{r eval = FALSE}
dsl_get()[1:2]
```

```
$before_install
── before_install ──────────────────────────────────────── stage ──
ℹ No steps defined

$install
── install ─────────────────────────────────────────────── stage ──
▶ step_install_github("ropensci/rotemplate")
▶ step_install_deps(repos = repo_default())
```


### Emulating a CI run locally

A tic configuration can be emulated locally.

First, ensure that all prerequisites are met:

```{r eval = FALSE}
prepare_all_stages()
```

This might install additional packages, or even fail with a descriptive message.

Then, run all steps defined in `tic.R`:

```{r}
run_all_stages()
```

This emulates a CI run on your local machine by executing all stages and their corresponding steps.
Note that this action this will use your local system libraries and not the CI environment.
Only the steps that are shown with `dsl_get()` are executed.
Some steps will not be executed as they are conditioned to run on non-interactive environments only, e.g. `add_step(covr::codcov())` added by the macro `do_package_checks()`.

```{r, eval = FALSE}
run_all_stages()
```

```
✓ Loading tic stage configuration from tic.R
Running install: step_install_github("ropensci/rotemplate")
Skipping install of 'rotemplate' from a github remote, the SHA1 (bec3e6eb) has not changed since last install.
  Use `force = TRUE` to force installation
Running install: step_install_deps(repos = repo_default())
```

## Debugging: Entering the CI build directly

### Circle CI

Debugging builds on Circle CI is very convenient.
Go to the web interface and click on the three dots as shown in the screenshot.
Then you can restart your build from the appearing dropdown menu and SSH into the build.

```{r, eval = TRUE, echo=FALSE, fig.align='center'}
knitr::include_graphics("img/circleci-debug.png")
```

# What's not covered yet?

- `SystemRequirements`: {tic} is not yet capable of automatically determining system requirements specified in DESCRIPTION files of an R package.
Future plans include to automatically provide suggestions like
---
title: "FAQ"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{FAQ}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Authentication

#### Q-Auth-1

I am getting an error when using `circle::use_circle_deploy()` or `tic::use_ghactions_deploy()`.

**Answer**

In most cases this is related to API authentication issues.
Ensure that the following points are met:

1. For Circle CI, install the respective GitHub App from the [GitHub Marketplace](https://github.com/marketplace).
2. Ensure that you have set the respective API keys for the problematic provider in your `.Renviron` file.
  Consult the help pages of the respective `use_*_deploy()` function for more help.
  - GitHub Actions: A `GITHUB_PAT` with "public_repo" scopes.
  - Circle CI: Env var `R_CIRCLE`.

## GitHub Actions

#### Q-GHA-1

How is {tic} different from what [r-lib/actions](https://github.com/r-lib/actions) does?

**Answer**

{tic} uses [r-lib/actions](https://github.com/r-lib/actions) as the base to install R in the first place.
However in detail, {tic} does the following things differently which aim to enhance the CI experience:

- Caching: {tic} caches the whole R library rather than only the direct dependencies of a package.
  This has the advantage that packages required for side actions ({pkgdown} deployment, README updates) will also be cached.

- `ccache`: {tic} comes with a compiler cache for source installations of packages by default, speeding up repeated source installation highly.
  The compiler cache directory (`.ccache`) will also be cached (once a week).
  Example use case: If you installed {Rcpp} from source as a dependency of your package and have it stored in your cache and {Rcpp} now updates two days later, the reinstallation will make use of the compiler cache and install {Rcpp} instantly rather than recompiling the C code again.

- Number of CPUs: {tic} uses 4 CPUs by default instead of only 1 as [r-lib/actions](https://github.com/r-lib/actions) does.
  This speeds up package installation a lot.
  4 CPUs are max because all GitHub Actions runners have 2 hyperthreading cores available.

- Use of SSH deployment keys: Once set up via `tic::use_ghactions_deploy()`, this deployment approach makes it possible to push any file of your repository to any branch of your remote.
  Other deployment approaches often limit you to only push to `gh-pages` branch or similar.

## Other

#### Q-Other-1

Is it possible to update the CI YAML templates installed by {tic} with upstream changes?

**Answer**

Yes! Have a look at ["Updating Templates"](https://docs.ropensci.org/tic/articles/updating.html) for more information.

---

#### Q-Other-2

Am I the only one using {tic}?

**Answer**

You can see who and how many people use {tic.R} on GitHub via this query: https://github.com/search?p=5&q=filename%3Atic.R&type=Code

---

#### Q-Other-3

Package {rgl} fails to install because of either

- "configure: error: X11 not found but required, configure aborted."
- "error: X11 not found; XQuartz (from www.xquartz.org) is required to run rgl."

**Answer**

The first one is usually caused by a missing installation of `XQuartz` on macOS.
Add `brew install xquartz` to the runner.

The second error requires to set the `DISPLAY` env var to mimic a non-headless state.
Add `export DISPLAY=:99` to the stage in which {rgl} should be installed.
If the warning message during loading of {rgl} should be suppressed, either env var `RGL_USE_NULL = TRUE` can be set or R option `options(rgl.useNull = TRUE)`.
---
title: "Developer info: Writing custom steps"
author: "Kirill Müller, Patrick Schratz"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Developer info: Writing custom steps}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

Most important steps running on CI systems are [already implemented](tic.html#steps) in {tic}.
The following vignette shows how new steps can be created and how {tic} defines them.
Internally, all steps are defined using the [R6](https://github.com/wch/R6) class system.
If you are not familiar with object-oriented programming using R6, the [R6 chapter](https://adv-r.hadley.nz/r6.html) in [Advanced R](https://adv-r.hadley.nz/) is a good place to start.

In most cases there should be no need to write steps on your own, because `step_run_code()` can be used to run arbitrary code when preparing and running the step, and custom conditions in `tic.R` can be used to trigger the step.
However, if you have the need for a repeated use of specific combinations of `add_code_step()`, please let us know (by opening an [issue](https://github.com/ropensci/tic/issues)) so that we can discuss whether it makes sense to implement it as a custom step.

## The `TicStep` class

A step is a subclass of the [`TicStep` class](https://github.com/ropensci/tic/blob/457a5d259c6861e322220ac51a0436436e5f214b/R/steps-base.R#L7-L38).
The `step_...` functions in {tic} are forwarded to the `new()` methods of the corresponding R6 class objects (and from there to `initialize()` member functions).
We recommend following the same pattern for your custom steps.

The `TicStep` class implements the public methods `check()`, `prepare()`, and `run()`.
These methods must be callable without arguments.
This means that you only need to override the methods you need; if you don't need a `check()` or `prepare()` method you can leave it out.
The following sections describe these methods and show examples.

### The `prepare()` method

This method will be called by `prepare_all_stages()`.
It is intended to run in the `before_script` phase of the CI run.
This method should install all dependent packages that the step requires, which then can be cached by the CI system.
You also may include further preparation code here.
For example `step_rcmdcheck` verifies that the R packages *rcmdcheck* and *remotes* are installed:

```{r}
RCMDcheck$public_methods$prepare
```

### The `run()` method

This method executes the step.
When a step is added to a stage, `run()` will be called when the stage is executed.
For example, the `run()` function of class `RCMDcheck` looks as follows:

```{r}
RCMDcheck$public_methods$run
```

### The `check()` method

This method checks whether the step is actually run.
It returns a logical scalar.
The [`step_write_text_file()`](https://github.com/ropensci/tic/blob/457a5d259c6861e322220ac51a0436436e5f214b/R/steps-write-text-file.R#L1-L24) function is an example step with the following implementation of the `check()` method:

```{r}
WriteTextFile$public_methods$check
```

## A minimal example

You can take a look at [a pull request that implements a new step](https://github.com/ropensci/tic/pull/75/files).

The most minimalist version is the "Hello World" [example step](https://github.com/ropensci/tic/blob/master/R/steps-base.R).
This class only contains a `run()` method which does nothing more than printing "Hello World" to the console.
It is initialized by calling `step_hello_world()` which creates a new instance of this class.

```{r eval = FALSE}
HelloWorld <- R6Class(
  "HelloWorld",
  inherit = TicStep,
  public = list(
    run = function() {
      print("Hello, world!")
    }
  )
)

#' Step: Hello, world!
#'
#' The simplest step possible: prints "Hello, world!" to the console when run, does not require
#' any preparation.
#' This step may be useful to test a \pkg{tic} setup or as a starting point when implementing a
#' custom step.
#'
#' @family steps
#' @export
step_hello_world <- function() {
  HelloWorld$new()
}
```

## Further information on the R6 class system

If you are unfamiliar with `R6` classes, here is a short guidance how the arguments are passed along:
Consider the `step_rcmdcheck()` function ([link to source](https://github.com/ropensci/tic/blob/master/R/steps-rcmdcheck.R)):

```{r, eval = FALSE}
RCMDcheck <- R6Class( # nolint
  "RCMDcheck",
  inherit = TicStep,
  public = list(
    initialize = function(warnings_are_errors = NULL, notes_are_errors = NULL,
                          args = c("--no-manual", "--as-cran"),
                          build_args = "--force", error_on = "warning",
                          repos = repo_default(), timeout = Inf,
                          check_dir = NULL) {
      if (!is.null(notes_are_errors)) {
        warning_once(
          '`notes_are_errors` is deprecated, please use `error_on = "note"`'
        )
        if (notes_are_errors) {
          error_on <- "note"
        }
      } else if (!is.null(warnings_are_errors)) {
        warning_once(
          "`warnings_are_errors` is deprecated, ",
          'please use `error_on = "warning"`'
        )
        if (warnings_are_errors) {
          error_on <- "warning"
        }
      }
      private$args <- args
      private$build_args <- build_args
      private$error_on <- error_on
      private$repos <- repos
      private$timeout <- timeout
      private$check_dir <- check_dir

      super$initialize()
    },
    run = function() {
      # Don't include vignettes if --no-build-vignettes is included
      if ("--no-build-vignettes" %in% private$args) {
        cat("^vignettes$\n", file = ".Rbuildignore", append = TRUE)
      }

      withr::with_envvar(
        c(
          # Avoid large version components
          "_R_CHECK_CRAN_INCOMING_" = "FALSE",
          # Don't check system clocks (because the API used there is flaky)
          "_R_CHECK_SYSTEM_CLOCK_" = "FALSE",
          # Don't force suggests
          "_R_CHECK_FORCE_SUGGESTS_" = "FALSE",
          # Work around missing qpdf executable
          "R_QPDF" = if (Sys.which("qpdf") == "") "true"
        ),
        res <- rcmdcheck::rcmdcheck(
          args = private$args, build_args = private$build_args,
          error_on = "never",
          repos = private$repos,
          timeout = private$timeout,
          check_dir = private$check_dir
        )
      )

      print(res)
      if (length(res$errors) > 0) {
        stopc("Errors found in rcmdcheck::rcmdcheck().")
      }
      if (private$error_on %in% c("warning", "note") && length(res$warnings) > 0) {
        stopc(
          "Warnings found in rcmdcheck::rcmdcheck(), ",
          'and `errors_on = "warning"` is set.'
        )
      }
      if (private$error_on == "note" && length(res$notes) > 0) {
        stopc(
          "Notes found in rcmdcheck::rcmdcheck(), ",
          'and `errors_on = "note"` is set.'
        )
      }
    },
    prepare = function() {
      verify_install("rcmdcheck")
      super$prepare()
    }
  ),
  private = list(
    args = NULL,
    build_args = NULL,
    error_on = NULL,
    repos = NULL,
    timeout = NULL,
    check_dir = NULL
  )
)
```

Here, a new instance of the defined `R6` class `RCMDcheck` is initiated with `RCMDcheck$new()`.
The arguments to `step_rcmdcheck()` are passed on to the `initialize()` function of the `R6` class.
Here, the arguments are assigned to the "private" members (e.g. `private$args`).
Next, these private members are used in the `run()` function which carries out the actual work.
---
title: "Getting started with CI for R"
author: "Patrick Schratz, Kirill Müller"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting started with CI for R}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library("tic")
```

# Prerequisites

Continuous Integration (CI) is a huge field in software development.
There are a lot of resources out there that try to explain what it can do and why you need it (and usually why their particular company/service does it best).
If you are just getting started with software development or `git` based project work you might think:

> "Ah I don't need all of this overhead, I am fine doing X and Y manually for this project.
> The overhead learning this next additional "tool" is not worth it.

Reading [all vignettes of this package](https://docs.ropensci.org/tic/articles/) will help you understand what {tic} does to simplify all of that specifically for the R language.

However, there is a lot of advanced information in these vignettes which might be a bit overwhelming at first if you are just getting started.
Therefore we will explain a few of the most important terms first to give you a kickstart:

<details>
<summary>Runner</summary>

One build job (possibly among many others) on any operating system that executes certain commands specified in the YAML config file.

</details>

<details>
<summary>YAML config file</summary>

A file written in the YAML language telling the runner what to do after a push to the repository.

</details>

<details>
<summary>Build matrix</summary>

Specification how many runners are started and their specification (operating system, custom environment variables, etc.)

</details>

<details>
<summary>CI Provider</summary>

A company that offers ready-to-use virtual images for runners which are started after a certain action (.e.g Circle CI).
However, code hosting sites like GitHub or GitLab also have their own CI integration meanwhile.

</details>

<details>
<summary>Deployment</summary>

CI builds cannot just perform specific checks on a package/project, they can also be used to (re-)build certain files in every run.
With the appropriate permissions these files can then be pushed to a repository via `git` without the need to take manual actions by the user.
This can save a lot of time and is commonly used to ensure that documentation is always up-to-date.

</details>

<details>
<summary>CI Client packages</summary>

R packages that interface the command line API of CI providers to simplify the execution of certain tasks (scraping build logs, enabling repos, etc.), e.g. [{circle}](https://github.com/ropenscilabs/circle).

</details>

<details>
<summary>CI Badges</summary>

Badges in the README of a repository showing the current build status of the project: Did the last build fail or finish successfully?

</details>

<details>
<summary>DSL</summary>

"DSL" stands for "[domain-specific language](https://en.wikipedia.org/wiki/Domain-specific_language)" and essentially means the implementation of a general concept for a specific programming language.

</details>

# Initialization/Setup

The easiest way to use {tic} for CI services is to call `tic::use_tic()`.
This will initialize a setup wizard which will guide you through all possibilities of the offered CI providers by {tic}.

Several yes/no questions need to be answered.
Based on the replies we'll select specific templates.

Besides the question which CI system you want to use, you'll be asked if you want to

- Deploy from builds to GitHub (e.g. if you are building a {pkgdown} site)

- Test your package on multiple R versions

Last, we'll add a `tic.R` file to the project root.

After this, your project is ready for continuous integration.
The next push to GitHub will create a build on Circle CI.
See the [Troubleshooting](advanced.html#troubleshooting) section in case anything doesn't work as expected.

#### Quickstart

If you are a new user, run

```r
tic::use_tic()
```

If you already use {tic} and want to configure a new CI provider, do

```r
## Circle CI
circle::use_circle_deploy() # (optional for deployment)
tic::use_circle_yml(<option here>)
```

---

If you are open to try out new things, Circle CI comes with some advantages that might simplify your CI experience.
However, all providers come with pros and cons and we cannot provide an exhaustive list comparing all providers here.

See the [CI Client Packages](ci-client-packages.html) article for more detailed information on how {tic} and the CI client packages work together.

# The role of the `tic.R` file

After having called `tic::use_tic()` you will find

- `.circleci/config.yml`
- `.github/workflows`

and a `tic.R` file in your repo, depending on the choices you made during `use_tic()`.
The latter will always be present because it will be the main CI config file for all providers from now on.
Usually you do not need to touch the YAML files anymore.
All build customization is done in `tic.R` and applies to all providers.
For more information about the build lifecycle in general, check the [Build lifecycle](build-lifecycle.html) article.

The basic `tic.R` template looks as follows:

```{r eval = FALSE}
do_package_checks()

if (ci_on_ghactions()) {
  do_pkgdown()
}
```

`tic.R` file has a declarative nature: It should consist of "stages", "steps" and "macro" functions (see below).
These functions will only have an effect when specified in `tic.R` and should **not** be used standalone as they will only run in a (simulated) CI run following a certain order.

To run plain R code within the build, encapsulate it within `add_code_step(<code>)` and add it to a certain build stage.
See [the Build Lifecycle](build-lifecycle.html) article for detailed info about how to do this.

## Macros {#macros}

{tic} builds on the "macro" idea.
Macros are essentially wrappers of a sequence of steps for often used tasks on the CI system:

- Checking a package (R CMD check)
- Building and deploying a pkgdown site
- Build and deploy a bookdown project

They can be distinguished from other functions by their `do_` prefix.
The following ones are currently implemented:

```{r}
list_macros()
```

If you have a good use case for a macro, let us know by opening an [issue](https://github.com/ropensci/tic/issues).

### `do_package_checks()`

```{r eval = FALSE}
do_package_checks()
```

`do_package_checks()` adds essential steps to various stages of a CI run.
Most importantly, it adds `step_rcmdcheck()` to the "script" stage.
This step performs the check of an R package.
Afterwards, the code coverage is being checked using `covr::codecov()`.
See `?do_package_checks()` for more information.

```{r}
# step_install_deps() in the "install" stage, using the repos argument.
#
# step_rcmdcheck() in the "script" stage, using the warnings_are_errors,
#  notes_are_errors, args, and build_args arguments.
#
# A call to covr::codecov() in the "after_success" stage (only if the codecov flag is set)
```

### `do_pkgdown()`

The other macro in the default template is `do_pkgdown()`.

```{r eval = FALSE}
if (ci_on_ghactions()) {
  do_pkgdown()
}
```

`do_pkgdown()` adds five steps to the build process:

```{reval = FALSE}
# step_install_deps() in the "install" stage, using the repos argument.
#
# step_setup_ssh() in the "before_deploy" to setup the upcoming deployment (if deploy is set),
#
# step_setup_push_deploy() in the "before_deploy" stage (if deploy is set),
#
# step_build_pkgdown() in the "deploy" stage, forwarding all ... arguments.
#
# step_do_push_deploy() in the "deploy" stage.
```

By default this currently happens only on GitHub Actions, because `ci_on_ghactions()` is used as a condition.
Why do we do this? Building the pkgdown site on multiple CI services has no added benefit and might even cause problems due to race conditions during deployment.

`ci_on_ghactions()` can be replaced by one of its sibling functions like `ci_on_circle()`.

### `do_readme_rmd()`

Some projects rely on a dynamic README.Rmd file which contains R code.
Sometimes the output of such README's will change over time if the code driving it changes due to updates.
To always stay up-to-date without needing to take manual action, you can use this macro.
It will render `README.Rmd` and deploy `README.md` to the default repo branch.
A deployment will only be made if the rendered output differs from the one stored upstream.

This macro requires that you have set up deployment for your selected provider beforehand.

### Blogdown

As a show case, we explain a "blogdown" project in more detail.
[`blogdown`](https://bookdown.org/yihui/blogdown/) is an R package for publishing websites.
Under the hood, it uses the framework [Hugo](https://gohugo.io/) which gets installed by the respective `tic.R` [template](https://github.com/krlmlr/tic.blogdown/blob/975aedd43fec1dd55e8348eccfca2c7c5f663006/tic.R) in the "install" section:

```{r eval = FALSE}
get_stage("install") %>%
  add_code_step(blogdown::install_hugo())
```

Next the website is built and deployed.
The `blogdown::build_site()` function for websites is the equivalent to `pkgdown::build_site()` for R packages.

```{r eval = FALSE}
get_stage("deploy") %>%
  add_code_step(blogdown::build_site()) %>%
  add_step(step_push_deploy())
```

Steps and stages differ between projects (e.g. between a "blogdown" website and a "package").
{tic} is smart enough to detect your project automatically when calling `tic::use_tic()` and will add the correct template.

**Note:** Currently, publishing to https://figshare.com/ doesn't work.
Also, publishing to https://zenodo.org/ is work in progress.

## {tic} projects from the community

The templates we provide with {tic} are minimal working examples.
By querying `tic.R` on GitHub one can see who else uses {tic} for their CI runs: https://github.com/search?p=5&q=filename%3Atic.R&type=Code

## Still got questions?

Have a look at the [list of articles](https://docs.ropensci.org/tic/articles/) we wrote to shine more light on all the parts {tic} covers.

If you face issues, make sure to also check out the [FAQ](https://docs.ropensci.org/tic/articles/faq.html) vignette or browse the [issue tracker](https://github.com/ropensci/tic/issues).
---
title: "Deployment"
author: "Patrick Schratz, Kirill Müller"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Deployment}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

{tic} uses [CI client packages](ci-client-packages.html) ({circle}) for setting up deployment on the CI systems.

### Circle CI

On Circle CI, setting up deployment convenient as there is no need to create a SSH key pair for deployment.

When calling `circle::use_circle_deploy()` directly (or indirectly via `tic::use_tic()`), a so called "user-key" is created and stored in the Circle CI repo.
This key makes it possible to deploy from Circle CI builds to your GitHub repo.
No extra deploy key need to be stored in the GitHub repo.

### GitHub Actions

GitHub Actions offers multiple options related to deployment:

- Deployment related actions on the GitHub marketplace
- Supplying an `DEPLOY_PAT` secret as [r-lib/actions](https://github.com/r-lib/actions/tree/master/examples#build-pkgdown-site) suggests
- Using an SSH key pair with the private key stored as a "secret" in the repo and the public key as a "Deploy key" ({tic} default)

To reflect a successful deployment in the repo checks, [actions from the GitHub marketplace](https://github.com/chrnorm/deployment-status) help.
For the actual deployment we recommend to use a SSH key pair.
A SSH key can be easily created, is safe to use and when browsing at the "Deploy key" section of a repo, one can directly see if deployment was granted to a CI service for the repository.

To simplify the creation of a SSH key pair for deployment and adding the keys to the appropriate places, {tic} comes with a little helper function `use_ghactions_deploy()`.
To use this function, you need a `GITHUB_PAT` with public repo scope defined in your `.Renviron`.
`usethis::browse_github_pat()` helps setting one up if you haven't done so already.

#### Updating the deployment status

To update the deployment status in the "environments" menu (next to "release" and "branches") conditionally on the outcome of the "Deployment" stage, one can use the `chrnorm/deployment-status`.

Add the following before the "Before Deploy" stage:

```yml
      - uses: chrnorm/deployment-action@v1.1.1
        name: Create GitHub deployment
        id: deployment
        with:
          token: "${{ github.token }}"
          environment: production
```yml

and then this part after the "After Deploy" stage:

```yml
      - name: Update deployment status (success)
        if: success()
        uses: chrnorm/deployment-status@releases/v1
        with:
          token: "${{ github.token }}"
          target-url: http://my-app-url.com
          state: "success"
          deployment_id: ${{ steps.deployment.outputs.deployment_id }}

      - name: Update deployment status (failure)
        if: failure()
        uses: chrnorm/deployment-status@releases/v1
        with:
          token: "${{ github.token }}"
          target-url: http://my-app-url.com
          state: "failure"
          deployment_id: ${{ steps.deployment.outputs.deployment_id }}
```

## {pkgdown} deployment

[{pkgdown}](https://github.com/r-lib/pkgdown) is an R package which builds a documentation wrapper-site of an R package.
It collects all the vignettes, function references and metadata information of a package and presents it in an eye-appealing HTML version.

It has become a quasi-standard to provide a {pkgdown} site for an R package.
However, it is tedious to update the {pkgdown} site manually on every commit, check whether something has changed and commit the changes.
{tic} comes with the ability to automate this process.

The following example shows how {tic} deploys a {pkgdown} site on GitHub Actions.
Remember that you can freely choose your favorite provider for this task.
In `tic.yml` file the "before_deploy" and "deploy" stages are redirected to {tic}.

```{yml eval = FALSE}
- name: "[Stage] Before Deploy"
  run: |
    Rscript -e "tic::before_deploy()"
```

In the "before_deploy" stage, {tic} will do the following:

```{r eval = FALSE}
if (ci_on_ghactions()) {
  get_stage("before_deploy") %>%
    add_step(step_setup_ssh())
}
```

- Calls `step_setup_ssh()` if the environment variable `"BUILD_PKGDOWN"` is set in the CI build.
  This step sets up SSH key previously added to the GitHub via `tic::use_ghactions_deploy()`.
  Argument `private_key_name` can be ignored as long as no custom private key name was used during `tic::use_ghactions_deploy()`.
  If so, then supply it via the `private_key_name` argument in `do_pkgdown()`.

  For backward compatibility, the deprecated default `"id_rsa"` is supported out of the box.

- Calls `add_step(step_build_pkgdown())` and `add_step(step_push_deploy())` in the "deploy" stage.

```{r eval = FALSE}
get_stage("deploy") %>%
  add_step(step_build_pkgdown()) %>%
  add_step(step_push_deploy())
```

`step_build_pkgdown()` will build the {pkgdown} site and afterwards (note the `pipe` operator chaining the commands), `step_push_deploy()` takes care pushing the results to the repo.
By default the site will be pushed to the `gh-pages` branch of your repo, keeping the history.

### Deploying to `docs/` (default branch) or `gh-pages` branch

By default deployment is pushed to the `gh-pages` branch.
This has the following advantages:

- No cluttering of the commit history in the default branch
- Everything "just works" silently in the background

#### Default branch deployment

Deploying to the `docs/` directory of the default branch has the following advantages:

- Per-branch versions of the {pkgdown} site (if desired)
- Per-branch versions enable the possibility to have preview for pull requests via https://www.netlify.com/

The disadvantage is that the default branch will be cluttered by automatic commits that push the changes of the {pkgdown} site to the default branch.

#### Orphaning the `gh-pages` branch

By default, changes to the {pkgdown} site will be added as incremental commits to the `gh-pages` branch.
This is useful to keep a history of past versions and to enable a release and dev version of the site.
To have this feature, set

```yml
development:
  mode: auto
```

in your `_pkgdown.yml` file.
See `?pkgdown::build_site()` for more information.

If you only want to have one version of your {pkgdown} site and not fill your repo with many commits in the `gh-pages` branch, you can use `do_pkgdown(orphan = TRUE)`.
This will wipe all commits of this branch on every CI run so that there is only one commit corresponding to the latest version of your pkgdown site.

## Committing single files

The `step_push_deploy()` step has the ability to restrict the files that are committed and pushed.
This can be very useful for conditionally pushing documentation files like `NEWS` or `man/` and `NAMESPACE` if these are automatically created during the CI run.

In the following example, these files are created/updated by calling `devtools::document()`.
The `commit_paths` argument in `step_push_deploy()` decides which files are committed and pushed:

```{r eval = FALSE}
get_stage("before_deploy") %>%
  add_step(step_setup_ssh())

get_stage("deploy") %>%
  add_code_step(devtools::document(roclets = c("rd", "collate", "namespace"))) %>%
  add_step(step_push_deploy(commit_paths = c("NAMESPACE", "man/*")))
```

Applying this idea depends on your overall R package development strategy: Commit files like `/man/` and `NAMESPACE` directly or let them be created during the CI run?
---
title: "tic & CI Client Packages - An Overview"
author: "Patrick Schratz, Kirill Müller"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{tic & CI Client Packages - An Overview}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

Most of the setup work that needs to be done for specific CI systems is handled by the respective R client packages such as [{circle}](https://docs.ropensci/circle/).
These enable the repo on the CI system and ensure that deployment permissions are granted.

After everything has been set up, the rest of the work goes to {tic}.

- Creation of the CI YAML templates.
- Which steps are going to be run.
- Deployment: Yes/no

In essence, `tic::use_tic()` is a wrapper for creating CI YAML templates and setting up deployment on the CI systems which is powered by the client packages.

## CI Client Packages

Currently, the following CI client packages exist:

- [{circle}]((https://docs.ropensci.org/circle/))

The client package for _travis_ was archived in November 2020.

For _Appveyor_ there is [r-appveyor](https://github.com/krlmlr/r-appveyor) from Kirill Müller.
This project makes it possible in the first place to run R checks on Appveyor but does not provide access to the _Appveyor_ API and is not used by {tic} currently.
Also, it does not come with automation for deployment setup.
Appveyor support has been removed from {tic} in December 2020 due to a strong focus on GitHub Actions.

For _GitHub Actions_ there is [{ghactions}](https://github.com/maxheld83/ghactions) from Max Held which comes with functions helping to set up YAML templates and other convenience.
It does not provide API access and is not used by {tic} currently.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_github.R
\name{github_repo}
\alias{github_repo}
\alias{github_info}
\alias{uses_github}
\title{Github information}
\usage{
github_repo(
  path = usethis::proj_get(),
  info = github_info(path, remote = remote),
  remote = "origin"
)

github_info(path = usethis::proj_get(), remote = "origin")

uses_github(path = usethis::proj_get())
}
\arguments{
\item{path}{\verb{[string]}\cr
The path to a GitHub-enabled Git repository (or a subdirectory thereof).}

\item{info}{\verb{[list]}\cr
GitHub information for the repository, by default obtained through
\code{\link[=github_info]{github_info()}}.}

\item{remote}{\verb{[string]}\cr
The Github remote which should be used. Defaults to "origin".}
}
\description{
\code{github_repo()} returns the true repository name as string.

Retrieves metadata about a Git repository from GitHub.

\code{github_info()} returns a list as obtained from the GET "/repos/:repo" API.
}
\concept{GitHub functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-ssh.R
\name{step_add_to_known_hosts}
\alias{step_add_to_known_hosts}
\title{Step: Add to known hosts}
\usage{
step_add_to_known_hosts(host = "github.com")
}
\arguments{
\item{host}{\verb{[string]}\cr
The host name to add to the \code{known_hosts} file, default: \code{github.com}.}
}
\description{
Adds a host name to the \verb{~/.ssh/known_hosts} file to allow subsequent
SSH access.
Requires \code{ssh-keyscan} on the system \code{PATH}.
}
\examples{
dsl_init()

get_stage("before_deploy") \%>\%
  add_step(step_add_to_known_hosts("gitlab.com"))

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{run_stage}
\alias{run_stage}
\title{Run a stage}
\usage{
run_stage(name, stages = dsl_load())
}
\arguments{
\item{name}{\verb{[string]}\cr
The name of the stage to run.}

\item{stages}{\verb{[named list]}
A named list of \code{TicStage} objects as returned by \code{\link[=dsl_load]{dsl_load()}},
by default loaded from \code{tic.R}.}
}
\description{
Run the \code{run_all()} method for all defined steps of a stage for which the
\code{check()} method returns \code{TRUE}.
}
\seealso{
\link{TicStep}

Other runners: 
\code{\link{prepare_all_stages}()},
\code{\link{run_all_stages}()}
}
\concept{runners}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use-badge.R
\name{use_tic_badge}
\alias{use_tic_badge}
\title{Add a CI Status Badge to README files}
\usage{
use_tic_badge(provider, branch = NULL, label = "tic")
}
\arguments{
\item{provider}{\code{character(1)}\cr
The CI provider to generate a badge for. Only \code{ghactions} is currently
supported}

\item{branch}{\code{character(1)}\cr
Which branch should the badge represent? Defaults to the default repo
branch.}

\item{label}{\code{character(1)}\cr
Text to use for the badge.}
}
\description{
Adds a CI status badge to \code{README.Rmd} or \code{README.md}. By default the label
is \code{"tic"}.

A custom branch can be specified via argument \code{branch}.
}
\examples{
\dontrun{
use_tic_badge(provider = "ghactions")

# use a different branch
use_tic_badge(provider = "ghactions", branch = "develop")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_tic.R
\name{use_tic_r}
\alias{use_tic_r}
\title{Add a tic.R file to the repo}
\usage{
use_tic_r(repo_type, deploy_on = "none")
}
\arguments{
\item{repo_type}{(\code{character(1)})\cr
Which type of template should be used. Possible values are \code{"package"},
\code{"site"}, \code{"blogdown"}, \code{"bookdown"} or \code{"unknown"}.}

\item{deploy_on}{(\code{character(1)})\cr
Which CI provider should perform deployment? Defaults to \code{NULL} which means
no deployment will be done. Possible values are \code{"ghactions"} or \code{"circle"}.}
}
\description{
Adds a \code{tic.R} file to containing the macros/steps/stages to be run during
CI runs.

The content depends on the repo type (detected automatically when used within
\code{\link[=use_tic]{use_tic()}}).
}
\examples{
\dontrun{
use_tic_r("package")
use_tic_r("package", deploy_on = "ghactions")
use_tic_r("blogdown", deploy_on = "all")
}
}
\seealso{
\link{yaml_templates}, \code{\link[=use_tic_badge]{use_tic_badge()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsl.R
\name{dsl}
\alias{dsl}
\alias{DSL}
\alias{get_stage}
\alias{add_step}
\alias{add_code_step}
\title{tic's domain-specific language}
\usage{
get_stage(name)

add_step(stage, step)

add_code_step(stage, call = NULL, prepare_call = NULL)
}
\arguments{
\item{name}{\verb{[string]}\cr
The name for the stage.}

\item{stage}{\verb{[TicStage]}\cr
A \code{TicStage} object as returned by \code{get_stage()}.}

\item{step}{\verb{[function]}\cr
An object of class \link{TicStep}, usually created by functions
with the \code{step_} prefix like \code{\link[=step_hello_world]{step_hello_world()}}.}

\item{call}{\verb{[call]}\cr
An arbitrary R expression executed during the stage to which this step is
added.
The default is useful if you only pass \code{prepare_call}.}

\item{prepare_call}{\verb{[call]}\cr
An optional arbitrary R expression executed during preparation.}
}
\description{
Functions to define stages and their constituent steps.
The \link{macro}s combine several steps and assign them to relevant
stages.
See \code{\link[=dsl_get]{dsl_get()}} for functions to access the storage for the stages
and their steps.

\code{get_stage()} returns a \code{TicStage} object for a stage given by name.
This function can be called directly in the \code{tic.R} configuration file,
which is processed by \code{\link[=dsl_load]{dsl_load()}}.

\code{add_step()} adds a step to a stage, see \code{\link[=step_hello_world]{step_hello_world()}}
and the links therein for available steps.

\code{add_code_step()} is a shortcut for \code{add_step(step_run_code(...))}.
}
\examples{
dsl_init()

get_stage("script")

get_stage("script") \%>\%
  add_step(step_hello_world())

get_stage("script")

get_stage("script") \%>\%
  add_code_step(print("Hi!"))

get_stage("script")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-git.R
\name{step_setup_push_deploy}
\alias{step_setup_push_deploy}
\title{Step: Setup push deploy}
\usage{
step_setup_push_deploy(
  path = ".",
  branch = NULL,
  orphan = FALSE,
  remote_url = NULL,
  checkout = TRUE
)
}
\arguments{
\item{path}{\verb{[string]}\cr
Path to the repository, default \code{"."} which means setting up the current
repository.}

\item{branch}{\verb{[string]}\cr
Target branch, default: current branch.}

\item{orphan}{\verb{[flag]}\cr
Create and force-push an orphan branch consisting of only one commit?
This can be useful e.g. for \verb{path = "docs", branch = "gh-pages"},
but cannot be applied for pushing to the current branch.}

\item{remote_url}{\verb{[string]}\cr
The URL of the remote Git repository to push to, defaults to the
current GitHub repository.}

\item{checkout}{\verb{[flag]}\cr
Check out the current contents of the repository? Defaults to \code{TRUE},
set to \code{FALSE} if the build process relies on existing contents or
if you deploy to a different branch.}
}
\description{
Clones a repo, inits author information, and sets up remotes
for a subsequent \code{\link[=step_do_push_deploy]{step_do_push_deploy()}}.
}
\examples{
\dontrun{
dsl_init()

get_stage("deploy") \%>\%
  add_step(step_setup_push_deploy(path = "docs", branch = "gh-pages")) \%>\%
  add_step(step_build_pkgdown())

# This example needs a Git repository
if (rlang::is_installed("git2r") && git2r::in_repository()) {
  # Deployment only works if a companion step_do_push_deploy() is added
  get_stage("deploy") \%>\%
    add_step(step_do_push_deploy(path = "docs"))
}

dsl_get()
}
}
\seealso{
Other deploy steps: 
\code{\link{step_do_push_deploy}()},
\code{\link{step_push_deploy}()}

Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{deploy steps}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-git.R
\name{step_do_push_deploy}
\alias{step_do_push_deploy}
\title{Step: Perform push deploy}
\usage{
step_do_push_deploy(
  path = ".",
  commit_message = NULL,
  commit_paths = ".",
  force = FALSE
)
}
\arguments{
\item{path}{\verb{[string]}\cr
Path to the repository, default \code{"."} which means setting up the current
repository.}

\item{commit_message}{\verb{[string]}\cr
Commit message to use, defaults to a useful message linking to the CI build
and avoiding recursive CI runs.}

\item{commit_paths}{\verb{[character]}\cr
Restrict the set of directories and/or files added to Git before deploying.
Default: deploy all files.}

\item{force}{\verb{[logical]}\cr
Add \code{--force} flag to git commands?}
}
\description{
Commits and pushes to a repo prepared by \code{\link[=step_setup_push_deploy]{step_setup_push_deploy()}}.

Deployment usually requires setting up SSH keys with
\code{\link[=use_tic]{use_tic()}}.
}
\details{
It is highly recommended to restrict the set of files
touched by the deployment with the \code{commit_paths} argument:
this step assumes that it can freely overwrite all changes to all files
below \code{commit_paths}, and will not warn in case of conflicts.

To mitigate conflicts race conditions to the greatest extent possible,
the following strategy is used:
\itemize{
\item The changes are committed to the branch
\item Before pushing, new commits are fetched, and the changes are cherry-picked
on top of the new commits
}

If no new commits were pushed after the CI run has started,
this strategy is equivalent to committing and pushing.
In the opposite case, if the remote repo has new commits,
the deployment is safely applied to the current tip.
}
\examples{
\dontrun{
dsl_init()

# Deployment only works if a companion step_setup_push_deploy() is added
get_stage("deploy") \%>\%
  add_step(step_setup_push_deploy(path = "docs", branch = "gh-pages")) \%>\%
  add_step(step_build_pkgdown())

if (rlang::is_installed("git2r") && git2r::in_repository()) {
  get_stage("deploy") \%>\%
    add_step(step_do_push_deploy(path = "docs"))
}

dsl_get()
}
}
\seealso{
Other deploy steps: 
\code{\link{step_push_deploy}()},
\code{\link{step_setup_push_deploy}()}

Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{deploy steps}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{prepare_all_stages}
\alias{prepare_all_stages}
\title{Prepare all stages}
\usage{
prepare_all_stages(stages = dsl_load())
}
\arguments{
\item{stages}{\verb{[named list]}
A named list of \code{TicStage} objects as returned by \code{\link[=dsl_load]{dsl_load()}},
by default loaded from \code{tic.R}.}
}
\description{
Run the \code{prepare()} method for all defined steps for which the
\code{check()} method returns \code{TRUE}.
}
\seealso{
\link{TicStep}

Other runners: 
\code{\link{run_all_stages}()},
\code{\link{run_stage}()}
}
\concept{runners}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gh-actions.R
\name{use_ghactions_deploy}
\alias{use_ghactions_deploy}
\title{Setup deployment for GitHub Actions}
\usage{
use_ghactions_deploy(
  path = usethis::proj_get(),
  repo = get_repo_slug(remote),
  key_name_private = "TIC_DEPLOY_KEY",
  key_name_public = "Deploy key for GitHub Actions",
  remote = "origin"
)
}
\arguments{
\item{path}{\verb{[string]} \cr
The path to the repository.}

\item{repo}{\verb{[string]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{key_name_private}{\verb{[string]}\cr
The name of the private key of the SSH key pair which will be created.
If not supplied, \code{"TIC_DEPLOY_KEY"} will be used.}

\item{key_name_public}{\verb{[string]}\cr
The name of the private key of the SSH key pair which will be created.
If not supplied, \code{"Deploy key for GitHub Actions"} will be used.}

\item{remote}{\verb{[string]}\cr
The GitHub remote which should be used. Defaults to "origin".}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

Creates a public-private key pair, adds the public key to the GitHub
repository via \code{github_add_key()}, and stores the private key as a "secret"
in the GitHub repo.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{run_all_stages}
\alias{run_all_stages}
\title{Emulate a CI run locally}
\usage{
run_all_stages(stages = dsl_load())
}
\arguments{
\item{stages}{\verb{[named list]}
A named list of \code{TicStage} objects as returned by \code{\link[=dsl_load]{dsl_load()}},
by default loaded from \code{tic.R}.}
}
\description{
Runs predefined \link{stages} similarly to the chosen CI provider.
The run aborts on error, the \code{after_failure} stage is never run.
}
\details{
The stages are run in the following order:
\enumerate{
\item \code{before_install()}
\item \code{install()}
\item \code{after_install()}
\item \code{before_script()}
\item \code{script()}
\item \code{after_success()}
\item \code{before_deploy()}
\item \code{deploy()}
\item \code{after_deploy()}
\item \code{after_script()}
}
}
\seealso{
Other runners: 
\code{\link{prepare_all_stages}()},
\code{\link{run_stage}()}
}
\concept{runners}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/macro.R, R/macro-package-checks.R,
%   R/macro-pkgdown.R, R/macro-blogdown.R, R/macro-bookdown.R, R/macro-drat.R,
%   R/macro-readme-rmd.R
\name{macro}
\alias{macro}
\title{Macros}
\description{
The \link{DSL} offers a fine-grained interface to the individual stages
of a CI run.
Macros are tic's way of adding several related steps to the relevant
stages.
All macros use the \code{do_} prefix.

The \code{\link[=do_package_checks]{do_package_checks()}} macro adds default checks for R packages,
including installation of dependencies and running a test coverage
analysis.

The \code{\link[=do_pkgdown]{do_pkgdown()}} macro adds the necessary steps for building
and deploying \pkg{pkgdown} documentation for a package.

The \code{\link[=do_blogdown]{do_blogdown()}} macro adds the necessary steps for building
and deploying a \pkg{blogdown} blog.

The \code{\link[=do_bookdown]{do_bookdown()}} macro adds the necessary steps for building
and deploying a \pkg{bookdown} book.

The \code{\link[=do_drat]{do_drat()}} macro adds the necessary steps for building
and deploying a drat repository to host R package sources.

The \code{\link[=do_readme_rmd]{do_readme_rmd()}} macro renders an R Markdown README and deploys
the rendered README.md file to Github.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-base.R
\name{TicStep}
\alias{TicStep}
\title{The base class for all steps}
\description{
Override this class to create a new step.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{TicStep$new()}}
\item \href{#method-run}{\code{TicStep$run()}}
\item \href{#method-prepare}{\code{TicStep$prepare()}}
\item \href{#method-check}{\code{TicStep$check()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a \code{TicStep} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TicStep$new()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-run"></a>}}
\if{latex}{\out{\hypertarget{method-run}{}}}
\subsection{Method \code{run()}}{
This method must be overridden, it is called when running the stage
to which a step has been added.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TicStep$run()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-prepare"></a>}}
\if{latex}{\out{\hypertarget{method-prepare}{}}}
\subsection{Method \code{prepare()}}{
This is just a placeholder.
This method is called when preparing the stage to
which a step has been added. It auto-install all packages which are
needed for a certain step. For example, \code{step_build_pkgdown()} requires
the \emph{pkgdown} package.

For \code{add_code_step()}, it autodetects any package calls in the form of
\code{pkg::fun} and tries to install these packages from CRAN. If a steps
\code{prepare_call} is not empty, the \verb{$prepare} method is skipped for this
step. This can be useful if a package should be installed from
non-standard repositories, e.g. from GitHub.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TicStep$prepare()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-check"></a>}}
\if{latex}{\out{\hypertarget{method-check}{}}}
\subsection{Method \code{check()}}{
This method determines if a step is prepared and run.
Return \code{FALSE} if conditions for running this step are not met.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TicStep$check()}\if{html}{\out{</div>}}
}

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update-yaml.R
\name{use_update_tic}
\alias{use_update_tic}
\title{Update tic Templates}
\usage{
use_update_tic()
}
\description{
Adds a GitHub Actions workflow (\code{update-tic.yml}) to check for tic template
updates once a day.

Internally, \code{\link[=update_yml]{update_yml()}} is called. A Pull Request will be opened if
a newer upstream version of the local tic template is found.

This workflow relies on a GITHUB_PAT with "workflow" scopes if GitHub Actions
templates should be updated.
Generate a GITHUB PAT and add it as a secret to your repo with
\code{\link[=gha_add_secret]{gha_add_secret()}}.
}
\examples{
\dontrun{
use_update_tic()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/macro.R
\name{list_macros}
\alias{list_macros}
\title{List available macros}
\usage{
list_macros()
}
\value{
\link{character}
}
\description{
Lists available macro functions of the \code{tic} package.
}
\seealso{
Other macros: 
\code{\link{do_blogdown}()},
\code{\link{do_bookdown}()},
\code{\link{do_drat}()},
\code{\link{do_package_checks}()},
\code{\link{do_pkgdown}()},
\code{\link{do_readme_rmd}()}
}
\concept{macros}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-write-text-file.R
\name{step_write_text_file}
\alias{step_write_text_file}
\title{Step: Write a text file}
\usage{
step_write_text_file(..., path)
}
\arguments{
\item{...}{\verb{[character]}\cr
Contents of the text file.}

\item{path}{\verb{[string]}\cr
Path to the new text file.}
}
\description{
Creates a text file with arbitrary contents
}
\examples{
dsl_init()

get_stage("script") \%>\%
  add_step(step_write_text_file("Hi!", path = "hello.txt"))

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-base.R
\name{step_hello_world}
\alias{step_hello_world}
\title{Step: Hello, world!}
\usage{
step_hello_world()
}
\description{
The simplest step possible: prints "Hello, world!" to the console when run,
does not require any preparation. This step may be useful to test a \pkg{tic}
setup or as a starting point when implementing a custom step.
}
\examples{
dsl_init()

get_stage("script") \%>\%
  add_step(step_hello_world())

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-pkgdown.R
\name{step_build_pkgdown}
\alias{step_build_pkgdown}
\title{Step: Build pkgdown documentation}
\usage{
step_build_pkgdown(...)
}
\arguments{
\item{...}{
  Arguments passed on to \code{\link[pkgdown:build_site]{pkgdown::build_site}}
  \describe{
    \item{\code{pkg}}{Path to package.}
    \item{\code{examples}}{Run examples?}
    \item{\code{run_dont_run}}{Run examples that are surrounded in \\dontrun?}
    \item{\code{seed}}{Seed used to initialize so that random examples are
reproducible.}
    \item{\code{lazy}}{If \code{TRUE}, will only rebuild articles and reference pages
if the source is newer than the destination.}
    \item{\code{override}}{An optional named list used to temporarily override
values in \verb{_pkgdown.yml}}
    \item{\code{preview}}{If \code{TRUE}, or \code{is.na(preview) && interactive()}, will preview
freshly generated section in browser.}
    \item{\code{devel}}{Use development or deployment process?

If \code{TRUE}, uses lighter-weight process suitable for rapid
iteration; it will run examples and vignettes in the current process,
and will load code with \code{pkgload::load_all()}.

If \code{FALSE}, will first install the package to a temporary library,
and will run all examples and vignettes in a new process.

\code{build_site()} defaults to \code{devel = FALSE} so that you get high fidelity
outputs when you building the complete site; \code{build_reference()},
\code{build_home()} and friends default to \code{devel = TRUE} so that you can
rapidly iterate during development.}
    \item{\code{new_process}}{If \code{TRUE}, will run \code{build_site()} in a separate process.
This enhances reproducibility by ensuring nothing that you have loaded
in the current process affects the build process.}
    \item{\code{install}}{If \code{TRUE}, will install the package in a temporary library
so it is available for vignettes.}
    \item{\code{document}}{\strong{Deprecated} Use \code{devel} instead.}
  }}
}
\description{
Builds package documentation with the \pkg{pkgdown} package.
Calls \code{pkgdown::clean_site()} and then \code{pkgdown::build_site(...)}.
}
\examples{
dsl_init()

get_stage("script") \%>\%
  add_step(step_build_pkgdown())

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsl-storage.R
\name{dsl_get}
\alias{dsl_get}
\alias{dsl_load}
\alias{dsl_init}
\title{Stages and steps}
\usage{
dsl_get()

dsl_load(path = "tic.R", force = FALSE, quiet = FALSE)

dsl_init(quiet = FALSE)
}
\arguments{
\item{path}{\verb{[string]}\cr
Path to the stage definition file, default: \code{"tic.R"}.}

\item{force}{\verb{[flag]}\cr
Set to \code{TRUE} to force loading from file even if a configuration exists.
By default an existing configuration is not overwritten by \code{dsl_load()}.}

\item{quiet}{\verb{[flag]}\cr
Set to \code{TRUE} to turn off verbose output.}
}
\value{
A named list of opaque stage objects with a \code{"class"} attribute
and a corresponding \code{\link[=print]{print()}} method for pretty output.
Use the high-level \code{\link[=get_stage]{get_stage()}} and \code{\link[=add_step]{add_step()}} functions to configure,
and the \link{stages} functions to run.
}
\description{
\pkg{tic} works in a declarative way, centered around the \code{tic.R} file
created by \code{\link[=use_tic]{use_tic()}}.
This file contains the \emph{definition} of the steps to be run in each stage:
calls to \code{\link[=get_stage]{get_stage()}} and \code{\link[=add_step]{add_step()}}, or macros like
\code{\link[=do_package_checks]{do_package_checks()}}.

Normally, this file is never executed directly.
Running these functions in an interactive session will \strong{not} carry out
the respective actions.
Instead, a description of the code that would have been run is printed
to the console.
Edit \code{tic.R} to configure your CI builds.
See \code{vignette("build-lifecycle", package = "tic")} for more details.
}
\details{
Stages and steps defined using tic's \link{DSL} are stored in an
internal object in the package.
The stages are accessible through \code{dsl_get()}.
When running the \link{stages}, by default a configuration defined
in the \code{tic.R} file is loaded with \code{dsl_load()}.
See \code{\link[=use_tic]{use_tic()}} for setting up a \code{tic.R} file.

For interactive tests, an empty storage can be initialized
with \code{dsl_init()}.
This happens automatically the first time \code{dsl_get()} is called
(directly or indirectly).
}
\examples{
\dontrun{
dsl_init()
dsl_get()

dsl_load(system.file("templates/package/tic.R", package = "tic"))
dsl_load(system.file("templates/package/tic.R", package = "tic"),
  force =
    TRUE
)
dsl_get()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_github.R
\name{ssh_key_helpers}
\alias{ssh_key_helpers}
\alias{github_add_key}
\alias{check_admin_repo}
\alias{get_role_in_repo}
\alias{get_public_key}
\alias{encode_private_key}
\alias{check_private_key_name}
\title{SSH key helpers}
\usage{
github_add_key(
  pubkey,
  repo = get_repo(remote),
  user = get_user(),
  title = "ghactions",
  remote = "origin",
  check_role = TRUE
)

check_admin_repo(owner, user, repo)

get_role_in_repo(owner, user, repo)

get_public_key(key)

encode_private_key(key)

check_private_key_name(string)
}
\arguments{
\item{pubkey}{The public key of the SSH key pair}

\item{repo}{\verb{[string]}\cr
The repository slug to use. Must follow the "\code{user/repo}" structure.}

\item{user}{The name of the user account}

\item{title}{The title of the key to add}

\item{remote}{\verb{[string]}\cr
The Github remote which should be used. Defaults to "origin".}

\item{check_role}{Whether to check if the current user has the permissions to
add a key to the repo. Setting this to \code{FALSE} makes it possible to add keys
to other repos than just the one from which the function is called.}

\item{owner}{The owner of the repository}

\item{key}{The SSH key pair object}

\item{string}{String to check}
}
\description{
SSH key helpers
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update-yaml.R
\name{update_yml}
\alias{update_yml}
\title{Update tic YAML Templates}
\usage{
update_yml(template_in = NULL, template_out = NULL)
}
\arguments{
\item{template_in}{\verb{[character]}\cr
Path to template which should be updated. By default all standard template
paths of GitHub Actions or Circle CI will be searched and
updated if they exist. Alternatively a full path to a single template can
be passed.}

\item{template_out}{\verb{[character]}\cr
Where the updated template should be written to. This is mainly used for
internal testing purposes and should not be set by the user.}
}
\description{
Updates YAML templates to their
latest versions. Currently only GitHub Actions and Circle CI templates are
supported.
}
\details{
By default all workflow files starting with \code{tic} are matched. This
means that you can have multiple YAML files with update support, e.g.
\code{"tic.yml"} and \code{"tic-db.yml"}.
}
\section{Formatting requirements of tic YAML templates}{
 To ensure that
updating of {tic} templates works, ensure the following points:
\itemize{
\item Your template contains the type (e.g. linux-matrix-deploy) and the revision
date in its first two lines.
\item When inserting comments into custom code blocks, only one-line comments are
allowed. Otherwise the update heuristic gets in trouble.
}
}

\examples{
\dontrun{
# auto-search
update_yml()

update_yml("tic.yml")

# custom named templates
update_yml("custom-name.yml")

# full paths
update_yml("~/path/to/repo/.github/workflows/tic.yml")
}
}
\seealso{
yaml_templates
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-blogdown.R
\name{step_build_blogdown}
\alias{step_build_blogdown}
\title{Step: Build a Blogdown Site}
\usage{
step_build_blogdown(...)
}
\arguments{
\item{...}{
  Arguments passed on to \code{\link[blogdown:build_site]{blogdown::build_site}}
  \describe{
    \item{\code{local}}{Whether to build the website locally. This argument is passed to
\code{\link[blogdown]{hugo_build}()}, and \code{local = TRUE} is mainly for serving
the site locally via \code{\link[blogdown]{serve_site}()}.}
    \item{\code{run_hugo}}{Whether to run \code{hugo_build()} after R Markdown files are
compiled.}
    \item{\code{build_rmd}}{Whether to (re)build R Markdown files. By default, they are
not built. See \sQuote{Details} for how \code{build_rmd = TRUE} works.
Alternatively, it can take a vector of file paths, which means these files
are to be (re)built. Or you can provide a function that takes a vector of
paths of all R Markdown files under the \file{content/} directory, and
returns a vector of paths of files to be built, e.g., \code{build_rmd =
blogdown::filter_timestamp}. A few aliases are currently provided for such
functions: \code{build_rmd = 'newfile'} is equivalent to \code{build_rmd =
blogdown::filter_newfile}, \code{build_rmd = 'timestamp'} is equivalent to
\code{build_rmd = blogdown::filter_timestamp}, and \code{build_rmd =
'md5sum'} is equivalent to \code{build_rmd = blogdown::filter_md5sum}.}
  }}
}
\description{
Build a Blogdown site using \code{\link[blogdown:build_site]{blogdown::build_site()}}.
}
\examples{
dsl_init()

get_stage("script") \%>\%
  add_step(step_build_blogdown("."))

dsl_get()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gh-actions.R
\name{gha_add_secret}
\alias{gha_add_secret}
\title{Add a GitHub Actions secret to a repository}
\usage{
gha_add_secret(
  secret,
  name,
  repo_slug = NULL,
  remote = "origin",
  visibility = "all",
  selected_repositories = NULL
)
}
\arguments{
\item{secret}{\verb{[character]}\cr
The value which should be encrypted (e.g. a Personal Access Token).}

\item{name}{\verb{[character]}\cr
The name of the secret as which it will appear in the "Secrets" overview of
the repository.}

\item{repo_slug}{\verb{[character]}\cr
Repository slug of the repository to which the secret should be added.
Must follow the form \code{owner/repo}.}

\item{remote}{\verb{[character]}\cr
If \code{repo_slug = NULL}, the \code{repo_slug} is determined by the respective git
remote.}

\item{visibility}{\verb{[character]}\cr
The level of visibility for the secret. One of \code{"all"}, \code{"private"}, or
\code{"selected"}.
See https://developer.github.com/v3/actions/secrets/#create-or-update-an-organization-secret
for more information.}

\item{selected_repositories}{\verb{[character]}\cr
Vector of repository ids for which the secret is accessible.
Only applies if \code{visibility = "selected"} was set.}
}
\description{
Encrypts the supplied value using \code{libsodium} and adds it as a
secret to the given GitHub repository. Secrets can be be used in GitHub
Action runs as environment variables.
A common use case is to encrypt Personal Access Tokens (PAT) or API keys.

This is the same as adding a secret manually in GitHub via
\code{"Settings" -> "Secrets" -> "New repository secret"}
}
\examples{
\dontrun{
gha_add_secret("supersecret", name = "MY_SECRET", repo = "ropensci/tic")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{stages}
\alias{stages}
\alias{before_install}
\alias{install}
\alias{after_install}
\alias{before_script}
\alias{script}
\alias{after_success}
\alias{after_failure}
\alias{before_deploy}
\alias{deploy}
\alias{after_deploy}
\alias{after_script}
\title{Predefined stages}
\usage{
before_install(stages = dsl_load())

install(stages = dsl_load())

after_install(stages = dsl_load())

before_script(stages = dsl_load())

script(stages = dsl_load())

after_success(stages = dsl_load())

after_failure(stages = dsl_load())

before_deploy(stages = dsl_load())

deploy(stages = dsl_load())

after_deploy(stages = dsl_load())

after_script(stages = dsl_load())
}
\arguments{
\item{stages}{\verb{[named list]}
A named list of \code{TicStage} objects as returned by \code{\link[=dsl_load]{dsl_load()}},
by default loaded from \code{tic.R}.}
}
\description{
Stages available in the CI provider, for which shortcuts
have been defined. All these functions call \code{\link[=run_stage]{run_stage()}} with the
corresponding stage name.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/macro-blogdown.R
\name{do_blogdown}
\alias{do_blogdown}
\title{Build a blogdown site}
\usage{
do_blogdown(
  ...,
  deploy = NULL,
  orphan = FALSE,
  checkout = TRUE,
  repos = repo_default(),
  path = "public",
  branch = "gh-pages",
  remote_url = NULL,
  commit_message = NULL,
  commit_paths = ".",
  force = FALSE,
  private_key_name = "TIC_DEPLOY_KEY",
  cname = NULL
)
}
\arguments{
\item{...}{Passed on to \code{\link[=step_build_blogdown]{step_build_blogdown()}}}

\item{deploy}{\verb{[flag]}\cr
If \code{TRUE}, deployment setup is performed
before building the blogdown site,
and the site is deployed after building it.
Set to \code{FALSE} to skip deployment.
By default (if \code{deploy} is \code{NULL}), deployment happens
if the following conditions are met:
\enumerate{
\item The repo can be pushed to (see \code{\link[=ci_can_push]{ci_can_push()}}).
\item The \code{branch} argument is \code{NULL}
(i.e., if the deployment happens to the active branch),
or the current branch is the default repo branch
(see \code{\link[=ci_get_branch]{ci_get_branch()}}).
}}

\item{orphan}{\verb{[flag]}\cr
Create and force-push an orphan branch consisting of only one commit?
This can be useful e.g. for \verb{path = "docs", branch = "gh-pages"},
but cannot be applied for pushing to the current branch.}

\item{checkout}{\verb{[flag]}\cr
Check out the current contents of the repository? Defaults to \code{TRUE},
set to \code{FALSE} if the build process relies on existing contents or
if you deploy to a different branch.}

\item{repos}{CRAN-like repositories to install from, defaults to
\code{\link[=repo_default]{repo_default()}}.}

\item{path}{\verb{[string]}\cr
Path to the repository, default \code{"."} which means setting up the current
repository.}

\item{branch}{\verb{[string]}\cr
Target branch, default: current branch.}

\item{remote_url}{\verb{[string]}\cr
The URL of the remote Git repository to push to, defaults to the
current GitHub repository.}

\item{commit_message}{\verb{[string]}\cr
Commit message to use, defaults to a useful message linking to the CI build
and avoiding recursive CI runs.}

\item{commit_paths}{\verb{[character]}\cr
Restrict the set of directories and/or files added to Git before deploying.
Default: deploy all files.}

\item{force}{\verb{[logical]}\cr
Add \code{--force} flag to git commands?}

\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}

\item{cname}{(\code{character(1)}\cr An optional URL for redirecting the created
website A \code{CNAME} file containing the given URL will be added to the root
of the directory specified in argument \code{path}.}
}
\description{
\code{do_blogdown()} adds default steps related to package checks
to the \code{"install"}, \code{"before_deploy"}, \code{"script"} and \code{"deploy"} stages.

\enumerate{
\item \code{\link[=step_install_deps]{step_install_deps()}} in the \code{"install"} stage, using the
\code{repos} argument.
\item \code{blogdown::install_hugo()} in the \code{"install"} stage to install the
latest version of HUGO.
\item \code{\link[=step_session_info]{step_session_info()}} in the \code{"install"} stage.
\item \code{\link[=step_setup_ssh]{step_setup_ssh()}} in the \code{"before_deploy"}
to setup the upcoming deployment (if \code{deploy} is set),
\item \code{\link[=step_setup_push_deploy]{step_setup_push_deploy()}} in the \code{"before_deploy"} stage
(if \code{deploy} is set),
\item \code{\link[=step_build_blogdown]{step_build_blogdown()}} in the \code{"deploy"} stage,
forwarding all \code{...} arguments.
\item \code{\link[=step_do_push_deploy]{step_do_push_deploy()}} in the \code{"deploy"} stage.
}

By default, the \verb{public/} directory is deployed to the \code{gh-pages} branch,
keeping the history. If the output directory of your blog/theme is not
\code{"public"} you need to change the \code{"path"} argument.
}
\examples{
\dontrun{
dsl_init()

do_blogdown()

dsl_get()
}
}
\seealso{
Other macros: 
\code{\link{do_bookdown}()},
\code{\link{do_drat}()},
\code{\link{do_package_checks}()},
\code{\link{do_pkgdown}()},
\code{\link{do_readme_rmd}()},
\code{\link{list_macros}()}
}
\concept{macros}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-git.R
\name{step_push_deploy}
\alias{step_push_deploy}
\title{Step: Setup and perform push deploy}
\usage{
step_push_deploy(
  path = ".",
  branch = NULL,
  remote_url = NULL,
  commit_message = NULL,
  commit_paths = ".",
  force = FALSE
)
}
\arguments{
\item{path}{\verb{[string]}\cr
Path to the repository, default \code{"."} which means setting up the current
repository.}

\item{branch}{\verb{[string]}\cr
Target branch, default: current branch.}

\item{remote_url}{\verb{[string]}\cr
The URL of the remote Git repository to push to, defaults to the
current GitHub repository.}

\item{commit_message}{\verb{[string]}\cr
Commit message to use, defaults to a useful message linking to the CI build
and avoiding recursive CI runs.}

\item{commit_paths}{\verb{[character]}\cr
Restrict the set of directories and/or files added to Git before deploying.
Default: deploy all files.}

\item{force}{\verb{[logical]}\cr
Add \code{--force} flag to git commands?}
}
\description{
Clones a repo, initializes author information, sets up remotes,
commits, and pushes.
Combines \code{\link[=step_setup_push_deploy]{step_setup_push_deploy()}} with \code{checkout = FALSE} and
a suitable \code{orphan} argument,
and \code{\link[=step_do_push_deploy]{step_do_push_deploy()}}.

Deployment usually requires setting up SSH keys with
\code{\link[=use_tic]{use_tic()}}.
}
\details{
Setup and deployment are combined in one step,
the files to be deployed must be prepared in a previous step.
This poses some restrictions on how the repository can be initialized,
in particular for a nonstandard \code{path} argument only \code{orphan = TRUE}
can be supported (and will be used).

For more control, create two separate steps with
\code{step_setup_push_deploy()} and \code{step_do_push_deploy()},
and create the files to be deployed in between these steps.
}
\examples{
\dontrun{
dsl_init()

get_stage("script") \%>\%
  add_step(step_push_deploy(commit_paths = c("NAMESPACE", "man")))

dsl_get()
}
}
\seealso{
Other deploy steps: 
\code{\link{step_do_push_deploy}()},
\code{\link{step_setup_push_deploy}()}

Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{deploy steps}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/macro-readme-rmd.R
\name{do_readme_rmd}
\alias{do_readme_rmd}
\title{Render a R Markdown README and deploy to Github}
\usage{
do_readme_rmd(
  checkout = TRUE,
  remote_url = NULL,
  commit_message = NULL,
  force = FALSE,
  private_key_name = "TIC_DEPLOY_KEY"
)
}
\arguments{
\item{checkout}{\verb{[flag]}\cr
Check out the current contents of the repository? Defaults to \code{TRUE},
set to \code{FALSE} if the build process relies on existing contents or
if you deploy to a different branch.}

\item{remote_url}{\verb{[string]}\cr
The URL of the remote Git repository to push to, defaults to the
current GitHub repository.}

\item{commit_message}{\verb{[string]}\cr
Commit message to use, defaults to a useful message linking to the CI build
and avoiding recursive CI runs.}

\item{force}{\verb{[logical]}\cr
Add \code{--force} flag to git commands?}

\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

\code{do_readme_rmd()} renders an R Markdown README and deploys
the rendered README.md file to Github. It adds default steps to the
\code{"before_deploy"} and \code{"deploy"} stages:

\enumerate{
\item \code{\link[=step_setup_ssh]{step_setup_ssh()}} in the \code{"before_deploy"} to setup
the upcoming deployment
\item \code{\link[=step_setup_push_deploy]{step_setup_push_deploy()}} in the \code{"before_deploy"} stage
\item \code{rmarkdown::render()} in the \code{"deploy"} stage
\item \code{\link[=step_do_push_deploy]{step_do_push_deploy()}} in the \code{"deploy"} stage.
}
}
\examples{
\dontrun{
dsl_init()

do_readme_rmd()

dsl_get()
}
}
\seealso{
Other macros: 
\code{\link{do_blogdown}()},
\code{\link{do_bookdown}()},
\code{\link{do_drat}()},
\code{\link{do_package_checks}()},
\code{\link{do_pkgdown}()},
\code{\link{list_macros}()}
}
\concept{macros}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_github.R
\name{github_helpers}
\alias{github_helpers}
\alias{auth_github}
\alias{get_owner}
\alias{get_user}
\alias{get_repo}
\alias{get_repo_slug}
\title{Github API helpers}
\usage{
auth_github()

get_owner(remote = "origin")

get_user()

get_repo(remote = "origin")

get_repo_slug(remote = "origin")
}
\arguments{
\item{remote}{\verb{[string]}\cr
The Github remote which should be used. Defaults to "origin".}
}
\description{
\itemize{
\item \code{auth_github()}: Creates a \code{GITHUB_TOKEN} and asks to store it in your
\code{.Renviron} file.
}

\itemize{
\item \code{get_owner()}: Returns the owner of a Github repo.
}

\itemize{
\item \code{get_repo()}: Returns the repo name of a Github repo for a given remote.
}

\itemize{
\item \code{get_repo_slug()}: Returns the repo slug of a Github repo
(\verb{<owner>/<repo>}).
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ci.R
\name{ci}
\alias{ci}
\alias{ci_get_branch}
\alias{ci_is_tag}
\alias{ci_get_slug}
\alias{ci_get_build_number}
\alias{ci_get_build_url}
\alias{ci_get_commit}
\alias{ci_get_env}
\alias{ci_is_env}
\alias{ci_has_env}
\alias{ci_can_push}
\alias{ci_is_interactive}
\alias{ci_cat_with_color}
\alias{ci_on_circle}
\alias{ci_on_ghactions}
\title{The current CI environment}
\usage{
ci_get_branch()

ci_is_tag()

ci_get_slug()

ci_get_build_number()

ci_get_build_url()

ci_get_commit()

ci_get_env(env)

ci_is_env(env, value)

ci_has_env(env)

ci_can_push(private_key_name = "TIC_DEPLOY_KEY")

ci_is_interactive()

ci_cat_with_color(code)

ci_on_circle()

ci_on_ghactions()

ci()
}
\arguments{
\item{env}{Name of the environment variable to check.}

\item{value}{Value for the environment variable to compare against.}

\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}

\item{code}{Code that should be colored.}
}
\description{
Functions that return environment settings that describe the CI
environment. The value is retrieved only once and then cached.

\code{ci_get_branch()}: Returns the current branch. Returns nothing if operating
on a tag.

\code{ci_is_tag()}: Returns the current tag name. Returns nothing if a branch is
selected.

\code{ci_get_slug()}: Returns the repo slug in the format \code{user/repo} or
\code{org/repo}

\code{ci_get_build_number()}: Returns the CI build number.

\code{ci_get_build_url()}: Returns the URL of the current build.

\code{ci_get_commit()}: Returns the SHA1 of the current commit.

\code{ci_get_env()}: Return an environment or configuration variable.

\code{ci_is_env()}: Checks if an environment or configuration variable is set to a
particular value.

\code{ci_has_env()}: Checks if an environment or configuration variable is set to
any value.

\code{ci_can_push()}: Checks if push deployment is possible. Always true
for local environments, CI environments require an environment
variable (by default \code{TIC_DEPLOY_KEY}).

\code{ci_is_interactive()}: Returns whether the current build is run interactively
or not. Global setup operations shouldn't be run on interactive CIs.

\code{ci_cat_with_color()}: Colored output targeted to the CI log.
The code argument can be an unevaluated call to a crayon function, the
style will be applied even if it normally wouldn't be.

\code{ci_on_circle()}: Are we running on Circle CI?

\code{ci_on_ghactions()}: Are we running on GitHub Actions?

\code{ci()}: Return the current CI environment
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsl.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-ssh.R
\name{step_install_ssh_keys}
\alias{step_install_ssh_keys}
\title{Step: Install an SSH key}
\usage{
step_install_ssh_keys(private_key_name = "TIC_DEPLOY_KEY")
}
\arguments{
\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}
}
\description{
Writes a private SSH key encoded in an environment variable
to a file in \verb{~/.ssh}.
Only run in non-interactive settings and if the environment variable
exists and is non-empty.
\code{\link[=use_ghactions_deploy]{use_ghactions_deploy()}} and \code{\link[=use_tic]{use_tic()}} functions encode a private key as an
environment variable for use with this function.
}
\examples{
dsl_init()

get_stage("before_deploy") \%>\%
  add_step(step_install_ssh_keys())

dsl_get()
}
\seealso{
\code{\link[=use_tic]{use_tic()}}, \code{\link[=use_ghactions_deploy]{use_ghactions_deploy()}}

Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/macro-package-checks.R
\name{do_package_checks}
\alias{do_package_checks}
\title{Add default checks for packages}
\usage{
do_package_checks(
  ...,
  codecov = !ci_is_interactive(),
  warnings_are_errors = NULL,
  notes_are_errors = NULL,
  args = NULL,
  build_args = NULL,
  error_on = "warning",
  repos = repo_default(),
  type = getOption("pkgType"),
  dependencies = TRUE,
  timeout = Inf,
  check_dir = "check"
)
}
\arguments{
\item{...}{Ignored, used to enforce naming of arguments.}

\item{codecov}{\verb{[flag]}\cr Whether to include a step running
\code{covr::codecov(quiet = FALSE)} (default: only for non-interactive CI,
see \code{\link[=ci_is_interactive]{ci_is_interactive()}}).}

\item{warnings_are_errors}{\verb{[flag]}\cr
Deprecated, use \code{error_on}.}

\item{notes_are_errors}{\verb{[flag]}\cr
Deprecated, use \code{error_on}.}

\item{args}{\verb{[character]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}.\cr

Default for local runs: \code{c("--no-manual", "--as-cran")}.

Default for Windows:
\code{c("--no-manual", "--as-cran", "--no-vignettes", "--no-build-vignettes", "--no-multiarch")}.

On GitHub Actions option "--no-manual" is always used (appended to custom
user input) because LaTeX is not available and installation is time
consuming and error prone.\cr}

\item{build_args}{\verb{[character]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}.\cr
Default for local runs: \code{"--force"}.\cr
Default for Windows: \code{c("--no-build-vignettes", "--force")}.\cr}

\item{error_on}{\verb{[character]}\cr
Whether to throw an error on R CMD check failures. Note that the check is
always completed (unless a timeout happens), and the error is only thrown
after completion. If "never", then no errors are thrown. If "error", then
only ERROR failures generate errors. If "warning", then WARNING failures
generate errors as well. If "note", then any check failure generated an
error.}

\item{repos}{\verb{[character]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}, default:
\code{\link[=repo_default]{repo_default()}}.}

\item{type}{Passed on to \code{\link[=install.packages]{install.packages()}}. The default avoids
installation from source on Windows and macOS by passing
\code{\link{.Platform}$pkgType}.}

\item{dependencies}{Which dependencies do you want to check?
Can be a character vector (selecting from "Depends", "Imports",
"LinkingTo", "Suggests", or "Enhances"), or a logical vector.

\code{TRUE} is shorthand for "Depends", "Imports", "LinkingTo" and
"Suggests". \code{NA} is shorthand for "Depends", "Imports" and "LinkingTo"
and is the default. \code{FALSE} is shorthand for no dependencies (i.e.
just check this package, not its dependencies).

The value "soft" means the same as \code{TRUE}, "hard" means the same as \code{NA}.

You can also specify dependencies from one or more additional fields,
common ones include:
\itemize{
\item Config/Needs/website - for dependencies used in building the pkgdown site.
\item Config/Needs/coverage for dependencies used in calculating test coverage.
}}

\item{timeout}{\verb{[numeric]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}, default:
\code{Inf}.}

\item{check_dir}{\verb{[character]} \cr Path specifying the directory for R CMD
check. Defaults to \code{"check"} for easy upload of artifacts.}
}
\description{
\code{do_package_checks()} adds default steps related to package checks
to the \code{"before_install"}, \code{"install"}, \code{"script"} and \code{"after_success"}
stages:

This macro is only available for R packages.

\enumerate{
\item \code{\link[=step_install_deps]{step_install_deps()}} in the \code{"install"} stage, using the
\code{repos} argument.
\item \code{\link[=step_session_info]{step_session_info()}} in the \code{"install"} stage.
\item \code{\link[=step_rcmdcheck]{step_rcmdcheck()}} in the \code{"script"} stage, using the
\code{warnings_are_errors}, \code{notes_are_errors}, \code{args}, and
\code{build_args} arguments.
\item A call to \code{\link[covr:codecov]{covr::codecov()}} in the \code{"after_success"} stage
(only if the \code{codecov} flag is set)
}
}
\examples{
dsl_init()

do_package_checks()

dsl_get()
}
\seealso{
Other macros: 
\code{\link{do_blogdown}()},
\code{\link{do_bookdown}()},
\code{\link{do_drat}()},
\code{\link{do_pkgdown}()},
\code{\link{do_readme_rmd}()},
\code{\link{list_macros}()}
}
\concept{macros}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-install.R
\name{step_install_pkg}
\alias{step_install_pkg}
\alias{step_install_deps}
\alias{step_install_cran}
\alias{step_install_github}
\title{Step: Install packages}
\usage{
step_install_deps(
  repos = repo_default(),
  type = getOption("pkgType"),
  dependencies = TRUE
)

step_install_cran(
  package = NULL,
  ...,
  repos = repo_default(),
  type = getOption("pkgType")
)

step_install_github(repo = NULL, ..., type = getOption("pkgType"))
}
\arguments{
\item{repos}{CRAN-like repositories to install from, defaults to
\code{\link[=repo_default]{repo_default()}}.}

\item{type}{Passed on to \code{\link[=install.packages]{install.packages()}}. The default avoids
installation from source on Windows and macOS by passing
\code{\link{.Platform}$pkgType}.}

\item{dependencies}{Which dependencies do you want to check?
Can be a character vector (selecting from "Depends", "Imports",
"LinkingTo", "Suggests", or "Enhances"), or a logical vector.

\code{TRUE} is shorthand for "Depends", "Imports", "LinkingTo" and
"Suggests". \code{NA} is shorthand for "Depends", "Imports" and "LinkingTo"
and is the default. \code{FALSE} is shorthand for no dependencies (i.e.
just check this package, not its dependencies).

The value "soft" means the same as \code{TRUE}, "hard" means the same as \code{NA}.

You can also specify dependencies from one or more additional fields,
common ones include:
\itemize{
\item Config/Needs/website - for dependencies used in building the pkgdown site.
\item Config/Needs/coverage for dependencies used in calculating test coverage.
}}

\item{package}{Package(s) to install}

\item{...}{Passed on to \code{install.packages()} or \code{remotes::install_github()}.}

\item{repo}{Package to install in the "user/repo" format.}
}
\description{
These steps are useful if your CI run needs additional packages.
Usually they are declared as dependencies in your \code{DESCRIPTION},
but it is also possible to install dependencies manually.
By default, binary versions of packages are installed if possible,
even if the CRAN version is ahead.

A \code{step_install_deps()} step installs all package dependencies declared in
\code{DESCRIPTION}, using \code{\link[remotes:install_deps]{remotes::install_deps()}}.
This includes upgrading outdated packages.

This step can only be used if a DESCRIPTION file is present in the repository
root.

A \code{step_install_cran()} step installs one package from CRAN via
\code{\link[=install.packages]{install.packages()}}, but only if it's not already installed.

A \code{step_install_github()} step installs one or more packages from GitHub
via \code{\link[remotes:install_github]{remotes::install_github()}}, the packages are only installed if their
GitHub version is different from the locally installed version.
}
\examples{
dsl_init()

get_stage("install") \%>\%
  add_step(step_install_deps())

dsl_get()
dsl_init()

get_stage("install") \%>\%
  add_step(step_install_cran("magick"))

dsl_get()
dsl_init()

get_stage("install") \%>\%
  add_step(step_install_github("rstudio/gt"))

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-drat.R
\name{step_add_to_drat}
\alias{step_add_to_drat}
\title{Step: Add built package to a drat}
\usage{
step_add_to_drat(repo_slug = NULL, deploy_dev = FALSE)
}
\arguments{
\item{repo_slug}{\verb{[string]}\cr
The name of the drat repository to deploy to in the form \verb{:owner/:repo}.}

\item{deploy_dev}{\verb{[logical]}\cr
Should development versions of packages also be deployed to the drat repo?
By default only "major", "minor" and "patch" releases are build and
deployed.}
}
\description{
Builds a package (binary on OS X or Windows) and inserts it into an existing
\pkg{drat} repository via \code{\link[drat:insertPackage]{drat::insertPackage()}}.
}
\examples{
dsl_init()

get_stage("script") \%>\%
  add_step(step_add_to_drat())

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-code.R
\name{step_run_code}
\alias{step_run_code}
\title{Step: Run arbitrary R code}
\usage{
step_run_code(call = NULL, prepare_call = NULL)
}
\arguments{
\item{call}{\verb{[call]}\cr
An arbitrary R expression executed during the stage to which this step is
added.
The default is useful if you only pass \code{prepare_call}.}

\item{prepare_call}{\verb{[call]}\cr
An optional arbitrary R expression executed during preparation.}
}
\description{
Captures the expression and executes it when running the step.
An optional preparatory expression can be provided that is executed
during preparation.
If the top-level expression is a qualified function call (of the format
\code{package::fun()}), the package is installed during preparation.
}
\examples{
dsl_init()

get_stage("install") \%>\%
  add_step(step_run_code(update.packages(ask = FALSE)))

# Will install covr from CRAN during preparation:
get_stage("after_success") \%>\%
  add_code_step(covr::codecov())

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use-yaml.R
\name{yaml_templates}
\alias{yaml_templates}
\alias{use_circle_yml}
\alias{use_ghactions_yml}
\title{Use CI YAML templates}
\usage{
use_circle_yml(type = "linux-matrix-deploy", write = TRUE, quiet = FALSE)

use_ghactions_yml(
  type = "linux-macos-windows-deploy",
  write = TRUE,
  quiet = FALSE
)
}
\arguments{
\item{type}{\verb{[character]}\cr
Which template to use. The string should be given following the logic
\verb{<platform>-<action>}. See details for more.}

\item{write}{\verb{[logical]}\cr
Whether to write the template to disk (\code{TRUE}) or just return it (\code{FALSE}).}

\item{quiet}{\verb{[logical]}\cr
Whether to print informative messages.}
}
\description{
Installs YAML templates for various CI providers. These functions
are also used within \code{\link[=use_tic]{use_tic()}}.

If you want to update an existing template use \code{\link[=update_yml]{update_yml()}}.
}
\section{pkgdown}{

If \code{type} contains "deploy", {tic} by default also sets the environment
variable \code{BUILD_PKGDOWN=true}. This triggers a call to
\code{pkgdown::build_site()} via the \code{do_pkgdown} macro in \code{tic.R} for the
respective runners.

If a setting  includes "matrix" and builds on multiple R versions, the job
building on R release is chosen to build the pkgdown site.
}

\section{YAML Type}{

\code{tic} supports a variety of different YAML templates which follow the
\verb{<platform>-<action>} pattern. The first one is mandatory, the
others are optional.
\itemize{
\item Possible values for \verb{<platform>} are \code{linux}, and \code{macos}, \code{windows}.
\item Possible values for \verb{<action>} are \code{matrix} and \code{deploy}.
}

Special types are \code{custom} and \code{custom-deploy}. These should be used if the
runner matrix is completely user-defined. This is mainly useful in
\code{\link[=update_yml]{update_yml()}}.

For backward compatibility \code{use_ghactions_yml()} will be default build and
deploy on all platforms.

Here is a list of all available combinations:\tabular{lllll}{
   Provider \tab Operating system \tab Deployment \tab multiple R versions \tab Call \cr
   Circle \tab Linux \tab no \tab no \tab \code{use_circle_yml("linux")} \cr
    \tab Linux \tab yes \tab no \tab \code{use_circle_yml("linux-deploy")} \cr
    \tab Linux \tab no \tab yes \tab \code{use_circle_yml("linux-matrix")} \cr
    \tab Linux \tab no \tab yes \tab \code{use_circle_yml("linux-deploy-matrix")} \cr
   ---------- \tab ------------------------ \tab ---------- \tab ------------------- \tab ------------------------------------------------------- \cr
   GH Actions \tab Linux \tab no \tab no \tab \code{use_ghactions_yml("linux")} \cr
    \tab Linux \tab yes \tab no \tab \code{use_ghactions_yml("linux-deploy")} \cr
    \tab custom \tab no \tab no \tab \code{use_ghactions_yml("custom")} \cr
    \tab custom-deploy \tab yes \tab no \tab \code{use_ghactions_yml("custom-deploy")} \cr
    \tab macOS \tab no \tab no \tab \code{use_ghactions_yml("macos")} \cr
    \tab macOS \tab yes \tab no \tab \code{use_ghactions_yml("macos-deploy")} \cr
    \tab Windows \tab no \tab no \tab \code{use_ghactions_yml("windows")} \cr
    \tab Windows \tab yes \tab no \tab \code{use_ghactions_yml("windows-deploy")} \cr
    \tab Linux + macOS \tab no \tab no \tab \code{use_ghactions_yml("linux-macos")} \cr
    \tab Linux + macOS \tab yes \tab no \tab \code{use_ghactions_yml("linux-macos-deploy")} \cr
    \tab Linux + Windows \tab no \tab no \tab \code{use_ghactions_yml("linux-windows")} \cr
    \tab Linux + Windows \tab yes \tab no \tab \code{use_ghactions_yml("linux-windows-deploy")} \cr
    \tab macOS + Windows \tab no \tab no \tab \code{use_ghactions_yml("macos-windows")} \cr
    \tab macOS + Windows \tab yes \tab no \tab \code{use_ghactions_yml("macos-windows-deploy")} \cr
    \tab Linux + macOS + Windows \tab no \tab no \tab \code{use_ghactions_yml("linux-macos-windows")} \cr
    \tab Linux + macOS + Windows \tab yes \tab no \tab \code{use_ghactions_yml("linux-macos-windows-deploy")} \cr
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/macro-pkgdown.R
\name{do_pkgdown}
\alias{do_pkgdown}
\title{Build pkgdown documentation}
\usage{
do_pkgdown(
  ...,
  deploy = NULL,
  orphan = FALSE,
  checkout = TRUE,
  repos = repo_default(),
  path = "docs",
  branch = "gh-pages",
  remote_url = NULL,
  commit_message = NULL,
  commit_paths = ".",
  force = FALSE,
  private_key_name = "TIC_DEPLOY_KEY"
)
}
\arguments{
\item{...}{Passed on to \code{\link[=step_build_pkgdown]{step_build_pkgdown()}}}

\item{deploy}{\verb{[flag]}\cr
If \code{TRUE}, deployment setup is performed
before building the pkgdown site,
and the site is deployed after building it.
Set to \code{FALSE} to skip deployment.
By default (if \code{deploy} is \code{NULL}), deployment happens
if the following conditions are met:
\enumerate{
\item The repo can be pushed to (see \code{\link[=ci_can_push]{ci_can_push()}}).
account for old default "id_rsa"
\item The \code{branch} argument is \code{NULL}
(i.e., if the deployment happens to the active branch),
or the current branch is the default branch,
or contains "cran-" in its name (for compatibility with \pkg{fledge})
(see \code{\link[=ci_get_branch]{ci_get_branch()}}).
}}

\item{orphan}{\verb{[flag]}\cr
Create and force-push an orphan branch consisting of only one commit?
This can be useful e.g. for \verb{path = "docs", branch = "gh-pages"},
but cannot be applied for pushing to the current branch.}

\item{checkout}{\verb{[flag]}\cr
Check out the current contents of the repository? Defaults to \code{TRUE},
set to \code{FALSE} if the build process relies on existing contents or
if you deploy to a different branch.}

\item{repos}{CRAN-like repositories to install from, defaults to
\code{\link[=repo_default]{repo_default()}}.}

\item{path, branch}{By default, this macro deploys the \code{docs} directory
to the \code{gh-pages} branch. This is different from \code{\link[=step_push_deploy]{step_push_deploy()}}.}

\item{remote_url}{\verb{[string]}\cr
The URL of the remote Git repository to push to, defaults to the
current GitHub repository.}

\item{commit_message}{\verb{[string]}\cr
Commit message to use, defaults to a useful message linking to the CI build
and avoiding recursive CI runs.}

\item{commit_paths}{\verb{[character]}\cr
Restrict the set of directories and/or files added to Git before deploying.
Default: deploy all files.}

\item{force}{\verb{[logical]}\cr
Add \code{--force} flag to git commands?}

\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}
}
\description{
\code{do_pkgdown()} builds and optionally deploys a pkgdown site and adds default
steps to the \code{"install"}, \code{"before_deploy"} and \code{"deploy"} stages:

\enumerate{
\item \code{\link[=step_install_deps]{step_install_deps()}} in the \code{"install"} stage, using the
\code{repos} argument.
\item \code{\link[=step_session_info]{step_session_info()}} in the \code{"install"} stage.
\item \code{\link[=step_setup_ssh]{step_setup_ssh()}} in the \code{"before_deploy"} to setup
the upcoming deployment (if \code{deploy} is set and only on
GitHub Actions),
\item \code{\link[=step_setup_push_deploy]{step_setup_push_deploy()}} in the \code{"before_deploy"} stage
(if \code{deploy} is set),
\item \code{\link[=step_build_pkgdown]{step_build_pkgdown()}} in the \code{"deploy"} stage,
forwarding all \code{...} arguments.
\item \code{\link[=step_do_push_deploy]{step_do_push_deploy()}} in the \code{"deploy"} stage.
}

By default, the \verb{docs/} directory is deployed to the \code{gh-pages} branch,
keeping the history.
}
\examples{
\dontrun{
dsl_init()

do_pkgdown()

dsl_get()
}
}
\seealso{
Other macros: 
\code{\link{do_blogdown}()},
\code{\link{do_bookdown}()},
\code{\link{do_drat}()},
\code{\link{do_package_checks}()},
\code{\link{do_readme_rmd}()},
\code{\link{list_macros}()}
}
\concept{macros}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/macro-bookdown.R
\name{do_bookdown}
\alias{do_bookdown}
\title{Build a bookdown book}
\usage{
do_bookdown(
  ...,
  deploy = NULL,
  orphan = FALSE,
  checkout = TRUE,
  repos = repo_default(),
  path = "_book",
  branch = "gh-pages",
  remote_url = NULL,
  commit_message = NULL,
  commit_paths = ".",
  force = FALSE,
  private_key_name = "TIC_DEPLOY_KEY",
  cname = NULL
)
}
\arguments{
\item{...}{Passed on to \code{\link[=step_build_bookdown]{step_build_bookdown()}}}

\item{deploy}{\verb{[flag]}\cr
If \code{TRUE}, deployment setup is performed
before building the bookdown site,
and the site is deployed after building it.
Set to \code{FALSE} to skip deployment.
By default (if \code{deploy} is \code{NULL}), deployment happens
if the following conditions are met:
\enumerate{
\item The repo can be pushed to (see \code{\link[=ci_can_push]{ci_can_push()}}).
\item The \code{branch} argument is \code{NULL}
(i.e., if the deployment happens to the active branch),
or the current branch is the default repo branch (usually "master")
(see \code{\link[=ci_get_branch]{ci_get_branch()}}).
}}

\item{orphan}{\verb{[flag]}\cr
Create and force-push an orphan branch consisting of only one commit?
This can be useful e.g. for \verb{path = "docs", branch = "gh-pages"},
but cannot be applied for pushing to the current branch.}

\item{checkout}{\verb{[flag]}\cr
Check out the current contents of the repository? Defaults to \code{TRUE},
set to \code{FALSE} if the build process relies on existing contents or
if you deploy to a different branch.}

\item{repos}{CRAN-like repositories to install from, defaults to
\code{\link[=repo_default]{repo_default()}}.}

\item{path}{\verb{[string]}\cr
Path to the repository, default \code{"."} which means setting up the current
repository.}

\item{branch}{\verb{[string]}\cr
Target branch, default: current branch.}

\item{remote_url}{\verb{[string]}\cr
The URL of the remote Git repository to push to, defaults to the
current GitHub repository.}

\item{commit_message}{\verb{[string]}\cr
Commit message to use, defaults to a useful message linking to the CI build
and avoiding recursive CI runs.}

\item{commit_paths}{\verb{[character]}\cr
Restrict the set of directories and/or files added to Git before deploying.
Default: deploy all files.}

\item{force}{\verb{[logical]}\cr
Add \code{--force} flag to git commands?}

\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}

\item{cname}{(\code{character(1)}\cr An optional URL for redirecting the created
website A \code{CNAME} file containing the given URL will be added to the root
of the directory specified in argument \code{path}.}
}
\description{
\code{do_bookdown()} adds default steps related to package checks
to the \code{"install"}, \code{"before_deploy"}, \code{"script"} and \code{"deploy"} stages.

\enumerate{
\item \code{\link[=step_install_deps]{step_install_deps()}} in the \code{"install"} stage, using the
\code{repos} argument.
\item \code{\link[=step_session_info]{step_session_info()}} in the \code{"install"} stage.
\item \code{\link[=step_setup_ssh]{step_setup_ssh()}} in the \code{"before_deploy"}
to setup the upcoming deployment (if \code{deploy} is set),
\item \code{\link[=step_setup_push_deploy]{step_setup_push_deploy()}} in the \code{"before_deploy"} stage
(if \code{deploy} is set),
\item \code{\link[=step_build_bookdown]{step_build_bookdown()}} in the \code{"deploy"} stage,
forwarding all \code{...} arguments.
\item \code{\link[=step_do_push_deploy]{step_do_push_deploy()}} in the \code{"deploy"} stage.
}

By default, the \verb{_book/} directory is deployed
to the \code{gh-pages} branch, keeping the history.
}
\examples{
\dontrun{
dsl_init()

do_bookdown()

dsl_get()
}
}
\seealso{
Other macros: 
\code{\link{do_blogdown}()},
\code{\link{do_drat}()},
\code{\link{do_package_checks}()},
\code{\link{do_pkgdown}()},
\code{\link{do_readme_rmd}()},
\code{\link{list_macros}()}
}
\concept{macros}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_tic.R
\name{use_tic}
\alias{use_tic}
\title{Initialize CI testing using tic}
\usage{
use_tic(
  wizard = interactive(),
  linux = "ghactions",
  mac = "ghactions",
  windows = "ghactions",
  deploy = "ghactions",
  matrix = "none",
  private_key_name = "TIC_DEPLOY_KEY",
  quiet = FALSE
)
}
\arguments{
\item{wizard}{\verb{[flag]}\cr Interactive operation? If \code{TRUE}, a menu will be
shown.}

\item{linux}{\verb{[string]}\cr Which CI provider(s) to use to test on Linux.
Possible options are \code{"circle"}, \code{"ghactions"}, \code{"none"}/\code{NULL}
and \code{"all"}.}

\item{mac}{\verb{[string]}\cr Which CI provider(s) to use to test on macOS
Possible options are \code{"none"}/\code{NULL} and \code{"ghactions"}.}

\item{windows}{\verb{[string]}\cr Which CI provider(s) to use to test on Windows
Possible options are \code{"none"}/\code{NULL}, and \code{"ghactions"}.}

\item{deploy}{\verb{[string]}\cr Which CI provider(s) to use to deploy artifacts
such as pkgdown documentation. Possible options are "circle"\verb{, }"ghactions"\verb{, }"none"\code{/}NULL\code{and}"all"`.}

\item{matrix}{\verb{[string]}\cr For which CI provider(s) to set up matrix builds.
Possible options are  \code{"circle"}, \code{"ghactions"}, \code{"none"}/\code{NULL}
and \code{"all"}.}

\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}

\item{quiet}{\verb{[flag]}\cr Less verbose output? Default: \code{FALSE}.}
}
\description{
Prepares a repo for building and deploying supported by
\pkg{tic}.
}
\details{
\enumerate{
\item Query information which CI providers should be used
\item Setup permissions for providers selected for deployment
\item Create YAML files for selected providers
\item Create a default \code{tic.R} file depending on the repo type
(package, website, bookdown, ...)
}
}
\examples{
# Requires interactive mode
if (FALSE) {
  use_tic()

  # Pre-specified settings favoring Circle CI:
  use_tic(
    wizard = FALSE,
    linux = "circle",
    mac = "ghactions",
    windows = "ghactions",
    deploy = "circle",
    matrix = "all"
  )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tic-package.R
\docType{package}
\name{tic-package}
\alias{tic-package}
\title{tic: Tasks Integrating Continuously: CI-Agnostic Workflow Definitions}
\description{
Provides a way to describe common build and deployment workflows for R-based projects: packages, websites (e.g. blogdown, pkgdown), or data processing (e.g. research compendia). The recipe is described independent of the continuous integration tool used for processing the workflow (e.g. 'GitHub Actions' or 'Circle CI'). This package has been peer-reviewed by rOpenSci (v0.3.0.9004).
}
\details{
The \code{\link[=use_tic]{use_tic()}} function prepares a code repository for use with
this package.
See \link{DSL} for an overview of \pkg{tic}'s domain-specific
language for defining stages and steps,
\code{\link[=step_hello_world]{step_hello_world()}} and the links therein for available steps,
and \link{macro} for an overview over the available macros that bundle
several steps.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/tic}
  \item Report bugs at \url{https://github.com/ropensci/tic/issues}
}

}
\author{
\strong{Maintainer}: Kirill Müller \email{krlmlr+r@mailbox.org} (\href{https://orcid.org/0000-0002-1416-3412}{ORCID})

Authors:
\itemize{
  \item Patrick Schratz \email{patrick.schratz@gmail.com} (\href{https://orcid.org/0000-0003-0748-6624}{ORCID})
  \item Mika Braginsky \email{mika.br@gmail.com}
  \item Karthik Ram \email{karthik.ram@gmail.com}
  \item Jeroen Ooms \email{jeroen.ooms@stat.ucla.edu}
}

Other contributors:
\itemize{
  \item Max Held (Max reviewed the package for ropensci, see <https://github.com/ropensci/software-review/issues/305>) [reviewer]
  \item Anna Krystalli (Anna reviewed the package for ropensci, see <https://github.com/ropensci/software-review/issues/305>) [reviewer]
  \item Laura DeCicco (Laura reviewed the package for ropensci, see <https://github.com/ropensci/software-review/issues/305>) [reviewer]
  \item rOpenSci [funder]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-bookdown.R
\name{step_build_bookdown}
\alias{step_build_bookdown}
\title{Step: Build a bookdown book}
\usage{
step_build_bookdown(...)
}
\arguments{
\item{...}{See \link[bookdown:render_book]{bookdown::render_book}.}
}
\description{
Build a bookdown book using \code{\link[bookdown:render_book]{bookdown::render_book()}}.
}
\examples{
dsl_init()

get_stage("script") \%>\%
  add_step(step_build_bookdown("."))

dsl_get()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base64.R
\name{base64serialize}
\alias{base64serialize}
\alias{base64unserialize}
\title{Helpers for converting R objects to strings and back}
\usage{
base64serialize(x, compression = "gzip")

base64unserialize(x, compression = "gzip")
}
\arguments{
\item{x}{Object to serialize or deserialize}

\item{compression}{Passed on as \code{type} argument to \code{\link[=memCompress]{memCompress()}} or
\code{\link[=memDecompress]{memDecompress()}}.}
}
\description{
\code{base64serialize()} converts an R object into a string suitable for storing
in an environment variable. Use this function for encoding entire R objects
(such as OAuth tokens).

\code{base64unserialize()} is the inverse operation to \code{base64serialize()}.
Use this function in your \code{tic.R} to access the R object previously encoded
by \code{base64serialize()}.
}
\examples{
serial <- base64serialize(1:10)
base64unserialize(serial)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-ssh.R
\name{step_setup_ssh}
\alias{step_setup_ssh}
\title{Step: Setup SSH}
\usage{
step_setup_ssh(
  private_key_name = "TIC_DEPLOY_KEY",
  host = "github.com",
  url = paste0("git@", host),
  verbose = ""
)
}
\arguments{
\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}

\item{host}{\verb{[string]}\cr
The host name to add to the \code{known_hosts} file, default: \code{github.com}.}

\item{url}{\verb{[string]}\cr
URL to establish SSH connection with, by default \code{git@github.com}}

\item{verbose}{\verb{[string]}\cr
Verbosity, by default \code{""}. Use \code{-v} or \code{"-vvv"} for more verbosity.}
}
\description{
Adds to known hosts, installs private key, and tests the connection.
Chaining \code{\link[=step_install_ssh_keys]{step_install_ssh_keys()}}, \code{\link[=step_add_to_known_hosts]{step_add_to_known_hosts()}}
and \code{\link[=step_test_ssh]{step_test_ssh()}}.
\code{\link[=use_tic]{use_tic()}} encodes a private key as an environment variable for use with
this function.
}
\examples{
dsl_init()

get_stage("script") \%>\%
  add_step(step_setup_ssh(host = "gitlab.com"))

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/macro-drat.R
\name{do_drat}
\alias{do_drat}
\title{Build and deploy drat repository}
\usage{
do_drat(
  repo_slug = NULL,
  orphan = FALSE,
  checkout = TRUE,
  path = "~/git/drat",
  branch = NULL,
  remote_url = NULL,
  commit_message = NULL,
  commit_paths = ".",
  force = FALSE,
  private_key_name = "TIC_DEPLOY_KEY",
  deploy_dev = FALSE
)
}
\arguments{
\item{repo_slug}{\verb{[string]}\cr
The name of the drat repository to deploy to in the form \verb{:owner/:repo}.}

\item{orphan}{\verb{[flag]}\cr
Create and force-push an orphan branch consisting of only one commit?
This can be useful e.g. for \verb{path = "docs", branch = "gh-pages"},
but cannot be applied for pushing to the current branch.}

\item{checkout}{\verb{[flag]}\cr
Check out the current contents of the repository? Defaults to \code{TRUE},
set to \code{FALSE} if the build process relies on existing contents or
if you deploy to a different branch.}

\item{path, branch}{By default, this macro deploys the default repo branch
(usually "master") of the drat repository.
An alternative option is \code{"gh-pages"}.}

\item{remote_url}{\verb{[string]}\cr
The URL of the remote Git repository to push to, defaults to the
current GitHub repository.}

\item{commit_message}{\verb{[string]}\cr
Commit message to use, defaults to a useful message linking to the CI build
and avoiding recursive CI runs.}

\item{commit_paths}{\verb{[character]}\cr
Restrict the set of directories and/or files added to Git before deploying.
Default: deploy all files.}

\item{force}{\verb{[logical]}\cr
Add \code{--force} flag to git commands?}

\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}

\item{deploy_dev}{\verb{[logical]}\cr
Should development versions of packages also be deployed to the drat repo?
By default only "major", "minor" and "patch" releases are build and
deployed.}
}
\description{
\code{do_drat()} builds and deploys R packages to a drat repository and adds
default steps to the \code{"install"}, \code{"before_deploy"} and \code{"deploy"} stages:

\enumerate{
\item \code{\link[=step_setup_ssh]{step_setup_ssh()}} in the \code{"before_deploy"} to setup
the upcoming deployment
\item \code{\link[=step_setup_push_deploy]{step_setup_push_deploy()}} in the \code{"before_deploy"} stage
(if \code{deploy} is set),
\item \code{\link[=step_add_to_drat]{step_add_to_drat()}} in the \code{"deploy"}
\item \code{\link[=step_do_push_deploy]{step_do_push_deploy()}} in the \code{"deploy"} stage.
}
}
\section{Deployment}{

Deployment can only happen to the default repo branch (usually "master") or
\code{gh-pages} branch because the GitHub Pages functionality from GitHub is
used to access the drat repository later on. You need to enable this
functionality when creating the drat repository on GitHub via \verb{Settings -> GitHub pages} and set it to the chosen setting here.

To build and deploy Windows and macOS binaries, builds with deployment
permissions need to be triggered.
Have a look at
\url{https://docs.ropensci.org/tic/articles/deployment.html} for
more information and instructions.
}

\examples{
dsl_init()

do_drat()

dsl_get()
}
\seealso{
Other macros: 
\code{\link{do_blogdown}()},
\code{\link{do_bookdown}()},
\code{\link{do_package_checks}()},
\code{\link{do_pkgdown}()},
\code{\link{do_readme_rmd}()},
\code{\link{list_macros}()}
}
\concept{macros}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-rcmdcheck.R
\name{step_rcmdcheck}
\alias{step_rcmdcheck}
\title{Step: Check a package}
\usage{
step_rcmdcheck(
  ...,
  warnings_are_errors = NULL,
  notes_are_errors = NULL,
  args = NULL,
  build_args = NULL,
  error_on = "warning",
  repos = repo_default(),
  timeout = Inf,
  check_dir = "check"
)
}
\arguments{
\item{...}{Ignored, used to enforce naming of arguments.}

\item{warnings_are_errors, notes_are_errors}{\verb{[flag]}\cr
Deprecated, use \code{error_on}.}

\item{args}{\verb{[character]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}.\cr

Default for local runs: \code{c("--no-manual", "--as-cran")}.

Default for Windows:
\code{c("--no-manual", "--as-cran", "--no-vignettes", "--no-build-vignettes", "--no-multiarch")}.

On GitHub Actions option "--no-manual" is always used (appended to custom
user input) because LaTeX is not available and installation is time
consuming and error prone.\cr}

\item{build_args}{\verb{[character]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}.\cr
Default for local runs: \code{"--force"}.\cr
Default for Windows: \code{c("--no-build-vignettes", "--force")}.\cr}

\item{error_on}{\verb{[character]}\cr
Whether to throw an error on R CMD check failures. Note that the check is
always completed (unless a timeout happens), and the error is only thrown
after completion. If "never", then no errors are thrown. If "error", then
only ERROR failures generate errors. If "warning", then WARNING failures
generate errors as well. If "note", then any check failure generated an
error.}

\item{repos}{\verb{[character]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}, default:
\code{\link[=repo_default]{repo_default()}}.}

\item{timeout}{\verb{[numeric]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}, default:
\code{Inf}.}

\item{check_dir}{\verb{[character]} \cr Path specifying the directory for R CMD
check. Defaults to \code{"check"} for easy upload of artifacts.}
}
\description{
Check a package using \code{\link[rcmdcheck:rcmdcheck]{rcmdcheck::rcmdcheck()}},
which ultimately calls \verb{R CMD check}.
}
\section{Updating of (dependency) packages}{

Packages shipped with the R-installation will not be updated as they will be
overwritten by the R-installer in each build.
If you want these package to be updated, please add the following
step to your workflow: \code{add_code_step(remotes::update_packages("<pkg>"))}.
}

\examples{
dsl_init()

get_stage("script") \%>\%
  add_step(step_rcmdcheck(error_on = "note", repos = repo_bioc()))

dsl_get()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/macro-package-checks.R
\name{Deprecated}
\alias{Deprecated}
\alias{add_package_checks}
\title{Deprecated functions}
\usage{
add_package_checks(
  ...,
  warnings_are_errors = NULL,
  notes_are_errors = NULL,
  args = c("--no-manual", "--as-cran"),
  build_args = "--force",
  error_on = "warning",
  repos = repo_default(),
  timeout = Inf
)
}
\arguments{
\item{...}{Ignored, used to enforce naming of arguments.}

\item{warnings_are_errors}{\verb{[flag]}\cr
Deprecated, use \code{error_on}.}

\item{notes_are_errors}{\verb{[flag]}\cr
Deprecated, use \code{error_on}.}

\item{args}{\verb{[character]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}.\cr

Default for local runs: \code{c("--no-manual", "--as-cran")}.

Default for Windows:
\code{c("--no-manual", "--as-cran", "--no-vignettes", "--no-build-vignettes", "--no-multiarch")}.

On GitHub Actions option "--no-manual" is always used (appended to custom
user input) because LaTeX is not available and installation is time
consuming and error prone.\cr}

\item{build_args}{\verb{[character]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}.\cr
Default for local runs: \code{"--force"}.\cr
Default for Windows: \code{c("--no-build-vignettes", "--force")}.\cr}

\item{error_on}{\verb{[character]}\cr
Whether to throw an error on R CMD check failures. Note that the check is
always completed (unless a timeout happens), and the error is only thrown
after completion. If "never", then no errors are thrown. If "error", then
only ERROR failures generate errors. If "warning", then WARNING failures
generate errors as well. If "note", then any check failure generated an
error.}

\item{repos}{\verb{[character]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}, default:
\code{\link[=repo_default]{repo_default()}}.}

\item{timeout}{\verb{[numeric]}\cr
Passed to \code{rcmdcheck::rcmdcheck()}, default:
\code{Inf}.}
}
\description{
\code{add_package_checks()} has been replaced by \code{\link[=do_package_checks]{do_package_checks()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-session-info.R
\name{step_session_info}
\alias{step_session_info}
\title{Step: Print the current Session Info}
\usage{
step_session_info()
}
\description{
Prints out the package information of the current session via
\code{\link[sessioninfo:session_info]{sessioninfo::session_info()}}.
}
\examples{
dsl_init()

get_stage("install") \%>\%
  add_step(step_session_info())

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_test_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repo.R
\name{repo}
\alias{repo}
\alias{repo_default}
\alias{repo_cloud}
\alias{repo_cran}
\alias{repo_bioc}
\title{Shortcuts for accessing CRAN-like repositories}
\usage{
repo_default()

repo_cloud()

repo_cran()

repo_bioc(base = repo_default())
}
\arguments{
\item{base}{The base repo to use, defaults to \code{repo_default()}.
Pass \code{NULL} to install only from Bioconductor repos.}
}
\description{
These functions can be used as convenient shortcuts
for the \code{repos} argument to e.g. \code{\link[=do_package_checks]{do_package_checks()}} and
\code{\link[=step_install_deps]{step_install_deps()}}.

\code{repo_default()} returns the value of the \code{"repos"} option,
or \code{repo_cloud()} if the option is not set.

\code{repo_cloud()} returns RStudio's CRAN mirror.

\code{repo_cran()} returns the master CRAN repo.

\code{repo_bioc()} returns Bioconductor repos from
\code{\link[remotes:bioc_install_repos]{remotes::bioc_install_repos()}}, in addition to the default repo.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps-ssh.R
\name{step_test_ssh}
\alias{step_test_ssh}
\title{Step: Test SSH connection}
\usage{
step_test_ssh(
  url = "git@github.com",
  verbose = "",
  private_key_name = "TIC_DEPLOY_KEY"
)
}
\arguments{
\item{url}{\verb{[string]}\cr
URL to establish SSH connection with, by default \code{git@github.com}}

\item{verbose}{\verb{[string]}\cr
Verbosity, by default \code{""}. Use \code{-v} or \code{"-vvv"} for more verbosity.}

\item{private_key_name}{\code{string}\cr
Only needed when deploying from builds on GitHub Actions.
If you have set a custom name for the private key during creation of the
SSH key pair via tic::use_ghactions_deploy()] or \code{\link[=use_tic]{use_tic()}}, pass this
name here.}
}
\description{
Establishes an SSH connection.
This step doesn't fail if the connection cannot be established,
but prints verbose output by default.
It is useful for troubleshooting deployment problems.
}
\examples{
dsl_init()

get_stage("script") \%>\%
  add_step(step_test_ssh(verbose = "-vvv"))

dsl_get()
}
\seealso{
Other steps: 
\code{\link{step_add_to_drat}()},
\code{\link{step_add_to_known_hosts}()},
\code{\link{step_build_pkgdown}()},
\code{\link{step_do_push_deploy}()},
\code{\link{step_hello_world}()},
\code{\link{step_install_pkg}},
\code{\link{step_install_ssh_keys}()},
\code{\link{step_push_deploy}()},
\code{\link{step_run_code}()},
\code{\link{step_session_info}()},
\code{\link{step_setup_push_deploy}()},
\code{\link{step_setup_ssh}()},
\code{\link{step_write_text_file}()}
}
\concept{steps}
