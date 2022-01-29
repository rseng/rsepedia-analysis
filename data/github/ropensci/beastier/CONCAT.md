# beastier

[![Peer Review Status](https://badges.ropensci.org/209_status.svg)](https://github.com/ropensci/onboarding/issues/209)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/beastier)](https://cran.r-project.org/package=beastier)
[![](http://cranlogs.r-pkg.org/badges/grand-total/beastier)]( https://CRAN.R-project.org/package=beastier)
[![](http://cranlogs.r-pkg.org/badges/beastier)](https://CRAN.R-project.org/package=beastier)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/115617629.svg)](https://zenodo.org/badge/latestdoi/115617629)

Branch   |[![GitHub Actions logo](man/figures/GitHubActions.png)](https://github.com/ropensci/beautier/actions)|[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
---------|-----------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------
`master` |![R-CMD-check](https://github.com/ropensci/beastier/workflows/R-CMD-check/badge.svg?branch=master)   |[![codecov.io](https://codecov.io/github/ropensci/beastier/coverage.svg?branch=master)](https://codecov.io/github/ropensci/beastier/branch/master)
`develop`|![R-CMD-check](https://github.com/ropensci/beastier/workflows/R-CMD-check/badge.svg?branch=develop)  |[![codecov.io](https://codecov.io/github/ropensci/beastier/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/beastier/branch/develop)

`beastier` is an R package to run BEAST2.

![beastier logo](man/figures/beastier_logo.png)

`beastier` is part of the [`babette`](https://github.com/ropensci/babette) package suite:

 * [`beautier`](https://github.com/ropensci/beautier) creates BEAST2 input (`.xml`) files.
 * [`beastier`](https://github.com/ropensci/beastier) runs BEAST2
 * [`mauricer`](https://github.com/ropensci/mauricer): install BEAST2 packages
 * [`tracerer`](https://github.com/ropensci/tracerer) pastes BEAST2 output (`.log`, `.trees`, etc) files.

Related R packages:

 * [`beastierinstall`](https://github.com/richelbilderbeek/beastierinstall): Install and uninstall BEAST2
 * [`beastier_on_windows`](https://github.com/richelbilderbeek/beastier_on_windows): Verify that `beastier` works on the Windows operating system
 * [`lumier`](https://github.com/ropensci/lumier): Shiny app to help create the function call needed

## Install BEAST2

Due to CRAN policy, beastier cannot install BEAST2.
As a workaround, the non-CRAN 
[`beastierinstall`](https://github.com/richelbilderbeek/beastierinstall) 
can be used.

To install BEAST2:

```
remotes::install_github("richelbilderbeek/beastierinstall")
beastierinstall::install_beast2()
```

## Example for `v2.1`

Run BEAST2:

```
output_state_filename <- "out.state"

run_beast2(
  input_filename = get_beastier_path("2_4.xml"),
  output_state_filename = output_state_filename
)
```

This will create the files as specified in the `2_4.xml` BEAST2 input file.

## Example for `v2.0.25`

```
output_log_filename <- "out.log"
output_trees_filename <- "out.trees"
output_state_filename <- "out.state"

run_beast2(
  input_filename = get_beastier_path("2_4.xml"),
  output_log_filename = output_log_filename,
  output_trees_filenames = output_trees_filename,
  output_state_filename = output_state_filename
)
```

Note that in this version, the filenames for the `.log`
and `.trees` files could be specified. This is unneeded: 
the `2_4.xml` BEAST2 input file specifies where these files will be stored:

```
<?xml [...]?><beast [...]>

[...]

<run [...]>

    [...]

    <logger id="tracelog" fileName="test_output_0.log" [...]>
        [...]
    </logger>

    [...]

    <logger id="treelog.t:[...]" fileName="$(tree).trees" [...]>
        [...]
    </logger>
</run>
</beast>
```

When using `beautier`, this can be specified in `create_mcmc`:

```
create_mcmc(
  tracelog = create_tracelog(
    filename = "my_trace.log"
  ),
  treeslog = create_treeslog(
    filename = "my_trees.trees"
  )
)
```

## [Install](doc/install.md)

See [install](doc/install.md).

## [FAQ](doc/faq.md)

See [FAQ](doc/faq.md)

## Missing features/unsupported

`beastier` cannot do everything `BEAST2` can. 

 * Remove: install BEAST2, use [`beastierinstall`](https://github.com/richelbilderbeek/beastierinstall)
 * Experimental: Continue a BEAST2 run
 * Untested: Setup BEAGLE

## There is a feature I miss

See [CONTRIBUTING](CONTRIBUTING.md), at `Submitting use cases`

## I want to collaborate

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting code'

## I think I have found a bug

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting bugs' 

## There's something else I want to say

Sure, just add an Issue. Or send an email.

## External links

 * [BEAST2 GitHub](https://github.com/CompEvol/beast2)

## Dependencies

Branch                     |[![GitHub Actions logo](man/figures/GitHubActions.png)](https://github.com/ropensci/beautier/actions)              |[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
---------------------------|-------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
`beautier` `master`        |![R-CMD-check](https://github.com/ropensci/beautier/workflows/R-CMD-check/badge.svg?branch=master)                 |[![codecov.io](https://codecov.io/github/ropensci/beautier/coverage.svg?branch=master)](https://codecov.io/github/ropensci/beautier/branch/master)
`beautier` `develop`       |![R-CMD-check](https://github.com/ropensci/beautier/workflows/R-CMD-check/badge.svg?branch=develop)                |[![codecov.io](https://codecov.io/github/ropensci/beautier/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/beautier/branch/develop)
`beastierinstall` `master` |![R-CMD-check](https://github.com/richelbilderbeek/beastierinstall/workflows/R-CMD-check/badge.svg?branch=master)  |[![codecov.io](https://codecov.io/github/richelbilderbeek/beastierinstall/coverage.svg?branch=master)](https://codecov.io/github/richelbilderbeek/beastierinstall/branch/master)
`beastierinstall` `develop`|![R-CMD-check](https://github.com/richelbilderbeek/beastierinstall/workflows/R-CMD-check/badge.svg?branch=develop) |[![codecov.io](https://codecov.io/github/richelbilderbeek/beastierinstall/coverage.svg?branch=develop)](https://codecov.io/github/richelbilderbeek/beastierinstall/branch/develop)

Branch                         |[![AppVeyor logo](man/figures/AppVeyor.png)](https://ci.appveyor.com/project/richelbilderbeek/beastier_on_windows/)
-------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
`beastier_on_windows` `master` |[![Build status](https://ci.appveyor.com/api/projects/status/ralex9sdnnxlwbgx/branch/master?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/beastier-on-windows/branch/master)
`beastier_on_windows` `develop`|[![Build status](https://ci.appveyor.com/api/projects/status/ralex9sdnnxlwbgx/branch/develop?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/beastier-on-windows/branch/develop)

## References

Article about `babette`:

 * Bilderbeek, Richèl JC, and Rampal S. Etienne. "`babette`: BEAUti 2, BEAST 2 and Tracer for R." Methods in Ecology and Evolution (2018). https://doi.org/10.1111/2041-210X.13032

FASTA files `anthus_aco.fas` and `anthus_nd2.fas` from:
 
 * Van Els, Paul, and Heraldo V. Norambuena. "A revision of species limits in Neotropical pipits Anthus based on multilocus genetic and vocal data." Ibis.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

# News

Newest versions at top.

## beastier 2.4.8 (2021-09-19)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * `check_empty_beastier_folder` works under Windows

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.4.7 (2021-09-06)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Allow to use treelog, screenlog and tracelog filenames in absent
    sub-sub-subfolders, fixes babette bug

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.4.6 (only released on CRAN)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Improved documentation, thanks Julia Haider

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.4.5 (2021-08-18)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Improved error handling
  * Less whitespace in `beastier_report`

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.4.4 (2021-08-17)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Produces no undeleted temporary files

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.4.3 (2021-07-09)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * `check_can_create_file` creates folder if needed

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.4.2 (2021-06-04)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Tested to work on Windows

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.4.1 (2021-05-30)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Cleans up all created temporary files

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.4 (2021-05-22)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * `install_beast2`, `upgrade_beast2`, `uninstall_beast2` are deprecated,
    as these violated CRAN policy. Thanks to Brian Ripley for sharing!
    The deprecation message will point users to the non-official
    `beastierinstall` package at `https://github.com/richelbilderbeek/beastierinstall`
  * Removed deprecated function `update_beastier`

## beastier 2.3.1 (2021-05-15)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * Can correctly call BEAST2 when spaces in 
    (1) BEAST2 bin filename, 
    (2) BEAST2 jar filename,
    (3) BEAST2 input filename 
    (4) BEAST2 state output filename.
    Thanks Jason Griffiths

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.3 (2021-05-14)

### NEW FEATURES

  * Add `continue_beast2` to continue a BEAST2 run

### MINOR IMPROVEMENTS

  * Builds on GitHub Actions
  * Added `is_on_github_actions` to detect the GitHub Actions environment
  * `is_on_ci` detects the GitHub Actions environment

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * `update_beastier`: gives deprecation message

## beastier 2.2.1 (2020-10-31)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * `get_default_beast2_jar_path` allows for a custom folder name
  * `install_beast2` is silent by default
  * `install_beast2` and `uninstall_beast2` give more information when
    verbose

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.2 (2020-10-16)

### NEW FEATURES

  * Add 'create_mcbette_beast2_options'

### MINOR IMPROVEMENTS

  * No `testthat` tests in code examples

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.1.3 (2020-08-05)

### NEW FEATURES

  * Add `rename_beast2_options_filenames`
  * Add `get_beast2_options_filenames`
  * Allow tildes in file paths, for example `~/beast2.xml`
  * Add `beastier_report` to help create a bug report
  * Can `install_beast2` for a specific version

### MINOR IMPROVEMENTS

  * Enable MacOS build again

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.1.2 (2020-01-06)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Documentation at [rOpenSci](https://docs.ropensci.org/beastier) to
    shows pictures

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.1.1 (2019-12-02)

### NEW FEATURES

  * Add vignettes
  
### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.1 (2019-10-27)

### NEW FEATURES

  * Follows `beautier` v2.3 interface
  
### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * The function arguments `beast2_output_log_filename` 
    and `beast2_output_trees_filename` have become obsolete. Use
    `create_tracelog()$filename` and `create_treelog()$filename` as
    part of an MCMC instead
  * The function argument `beast2_working_dir` is obsoleted, because it
    no longer served a purpose. Use `setwd` to set the working directory 

## beastier 2.0.25 (2019-10-10)

### NEW FEATURES

  * Add `get_beast2_example_filenames` to get the
    full paths of all the BEAST2 example files
  * Add `get_beast2_example_filename` to get the
    full path of a BEAST2 example file
  * Add `get_default_beast2_working_dir` to get the full path of
    the default BEAST2 working directory
  * Add `clear_beast2_working_dir` to clear the BEAST2 working directory
  
### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None


## beastier 2.0.24 (2019-09-30)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Do not check if `beast2_path` still/already exists in `check_beast2_options`

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.23 (2019-09-29)

### NEW FEATURES

  * Use BEAST2 main class name when calling the jar file, 
    to fix the 'no main manifest attribute' error

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.22 (2019-09-16)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Add `update_beastier` to update the `beastier` dependencies
  * Lowered cyclomatic complexity, thanks @lintr-bot 

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.21 (2019-09-10)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * `create_beast2_options` checks that the filenames are not in the same
    folder as the BEAST2 working directory

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.20 (2019-08-27)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Added `get_beast2_options_filenames` to obtain the filenames in a
    `beast2_options` 

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None


## beastier 2.0.19 (2019-08-23)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * `check_beast2_options` checks for duplicate filenames

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.18 (2019-08-14)

### NEW FEATURES

  * Install BEAST2 `v2.6.0`, instead of `v2.5.2`

### MINOR IMPROVEMENTS

  * Can create files in sub-sub-subfolders
  * Builds on Travis CI distribution 'Bionic'
  * Add `are_identical_alignments`
  * Add `is_alignment`

### BUG FIXES

  * RNG seed is checked to be a whole number (Issue #35, thanks @thijsjanzen)

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.17

This version never had a formal release.

## beastier 2.0.16

This version never had a formal release.

## beastier 2.0.15 (2019-06-01)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Simplified `run_beast2`
  * Improved error checking and logging in `run_beast2`
  * Will overwrite files if requested (which is by default)
  * Removed example code that assumed all Linuxes have `/home/username`,
    thanks to Brian Ripley 

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.14 (2019-05-27)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Add `SystemRequirements` field to `DESCRIPTION`

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.13 (2019-05-14)

### NEW FEATURES

  * Install BEAST2 `v2.5.2`, instead of `v2.5.1`
  * Can call `run_beast2`/`run_beast2_from_options` in parallel

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.12 (2019-04-08)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Package builds successfully without BEAST2 installed

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.11 (2019-03-27)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * `get_default_beast2_bin_path` is silent if BEAST2 is absent
  * Package builds successfully without BEAST2 installed

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.10 (2019-03-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Passes tests without BEAST2 installed

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.9 (2019-03-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Pass local CRAN tests
  * Shortened examples' runtime 

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.8 (2019-03-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Fix spelling error

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.7 (2019-03-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Documentation examples work on macOS

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.6 (2019-03-25)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Whitelist 'Schmirl'

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.5 (2019-03-25)

### NEW FEATURES

  * CRAN and rOpenSci release candidate

### MINOR IMPROVEMENTS

  * Implemented feedback from CRAN, thanks Martina Schmirl:
    * In DESCRIPTION, use `'R'` for programs
    * In DESCRIPTION, use four spaces
    * In DESCRIPTION, remove the redundant '... from R' in title
    * In DESCRIPTION, use only `GPL-3` instead of also adding the file `LICENSE`
    * Add examples to all functions

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0.1 (2019-03-15)

### NEW FEATURES

  * CRAN and rOpenSci release candidate

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 2.0 (2019-01-04)

### NEW FEATURES

  * GitHub repository is owned by `ropensci`

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 1.5.2 (2018-10-30)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Tested to work under macOS

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 1.5.1 (2018-10-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Made up-to-date with newer version of `beautier`

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 1.5 (2018-09-12)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Instead of supplying a BEAST2 `.jar` path 
    (using the argument `beast2_jar_path`), one
    can now also supply the BEAST2 binary path.

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * All parameters named `beast2_jar_path` are now called `beast2_path`


## beastier 1.4.4 (2018-08-05)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * `run_beast2` had a parameter called `overwrite_state_file`, which is renamed 
    to `overwrite` and states if the `.log` file and `.trees` files will be overwritten, 
    like BEAST2 does. This is a bug fix, as the obsolete parameter name did not
    prevent the state file being overwritten.

### DEPRECATED AND DEFUNCT

  * Calling `run_beast2` with a parameter called `overwrite_state_file` is obsolete. 
    Use the parameter `overwrite` instead.

## beastier 1.4.3 (2018-05-17)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Tagged for [the academic article about `babette`](https://github.com/ropensci/babette_article)

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beastier 1.4.2 (2018-04-05)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Follow all [rOpenSci packaging guidelines](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md)

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None
# Contributing

Awesome that you are reading this.

This GitHub follows the [Contributor Covenant Code of Conduct](code_of_conduct.md).

 * For questions, you can create an Issue
 * Code changes go via Pull Requests

## Which package to contribute to?

`beastier` is part of the `babette` package suite,
which consists out of five packages.
Here is how to determine which package is best suited for your contribution:

If you want to contribute to the creation of BEAST2 XML input files, 
go to [beautier](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md).

If you want to contribute to how BEAST2 output is parsed,
go to [tracerer](https://github.com/ropensci/tracerer/blob/master/CONTRIBUTING.md)

If you want to contribute regarding the BEAST2 package management,
go to [mauricer](https://github.com/ropensci/mauricer/blob/master/CONTRIBUTING.md)

If you want to contribute with an overarching idea,
go to [babette](https://github.com/ropensci/babette/blob/master/CONTRIBUTING.md).

If you want to contribute to how BEAST2 is run,
you are at the right spot :-) 

## Submitting code

Submitted code should follow these quality guidelines:

 * All tests pass cleanly/silently
 * Code coverage must be 100%
 * Coding style should follow the default style by `lintr`

These are all checked by Travis CI when submitting
a Pull Request. 

Emails with code will not be accepted.

## Submitting bugs

Awesome. These are your options:

 * Add an Issue, with the test that fails
 * Submit a Pull Request, where the test is added to the `tests/testthat` folder
 * Send @richelbilderbeek an email (@richelbilderbeek will make an Issue of it)

Pull Requests should follow the same guidelines as 'Submitting code'.

## Branching policy

 * The `master` branch should always build successfully
 * The `development` branch is for developers

## git usage

To get started working on `beastier` do:

```
git clone https://github.com/ropensci/beastier
```

Development is done on the `develop` branch. 
To download and checkout the `develop` branch, 
first go into the `beastier` folder (`cd beastier`), then do:

```
git checkout develop
```

Then the workflow is the common `git` workflow:

```
git pull
git add --all :/
git commit -m "Did something awesome"
git push
```
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at richel@richelbilderbeek.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality regarding the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
Hi @richelbilderbeek,

With this Pull Request I'd would like to [add reason].

Sure, I've read [CONTRIBUTING.md](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md) :+1:

Cheers, [your name]

---
name: Custom issue
about: Anything else
title: ''
labels: ''
assignees: ''

---


---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Script to reproduce the behavior:

```r
# Your R script here, without this comment :-)
```

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Environment:**
Show the results of running the following script:

```r
library(beastier)
beastier::beastier_report()
```

**Additional context**
Add any other context about the problem here.

---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# extdata

# `beast beast2_error.xml`

File to create a BEAST2 error.

:warning: the BEAST2 NS package must not be installed :warning:

Upon running this file in BEAST2, the output will end with:

```
Error 1017 parsing the xml input file

Class could not be found. Did you mean beast.core.util.ESS?
Perhaps a package required for this class is not installed?

Error detected about here:
  <beast>
      <run id='mcmc' spec='beast.gss.NS'>
```

The error message says it all: one needs to install the BEAST2 NS package.
`beastier` cannot do this (note that `mauricer` can).

## `beast2_warning.xml`

File to create a BEAST2 warning.

Upon running this file in BEAST2, the output will end with:

```
WARNING: StateNode (freqParameter.s:anthus_aco) found that has no effect on posterior!

[...]
Fatal exception: Could not find a proper state to initialise. Perhaps try another seed.
See http://www.beast2.org/2018/07/04/fatal-errors.html for other possible solutions.
[...]
```

The word `WARNING` is used to detect a warning. Regardless of the 
fatal exceptions, BEAST2 exits without an error code (i.e. exit code zero).

# Install

This page described how to install `beastier`.

## Install `rJava`

If you have problems installing rJava, [Duck](http://www.duckduckgo.com) or [view my rJava notes](rjava.md).

## Install BEAST2

Due to CRAN policy, beastier cannot install BEAST2.
As a workaround, the non-CRAN `beastierinstall` can be used.

To install BEAST2:

```
remotes::install_github("richelbilderbeek/beastierinstall")
beastierinstall::install_beast2()
```

This will download and extract BEAST2 to:

OS     |Full path
-------|----------------------------------
Linux  |`~/.local/share/beast`
macOS  |`~/.local/share/beast`
Windows|`C:/Users/<username>/Local/beast`

## Install `beastier`

### CRAN version

### Stable development

`beastier` is not on CRAN yet. The simplest way now is to install `beastier` with the `devtools` R package:

```
remotes::install_github("ropensci/beastier")
```

`beastier` assumes that `beautier` and `tracerer` are also installed. Do so:

```
remotes::install_github("ropensci/beautier")
remotes::install_github("ropensci/tracerer")
```

`beastier` assumes that BEAST2 is installed. 

# `rJava` questions

How to install `rJava` under different operating systems

 * Installation
 * Troubleshooting

## Installation

### Ubuntu 14.5 (Trusty Tahr)

The `.travis.yml` file shows a Trusty install:

```
 - sudo apt-get install -qq oracle-java8-installer # Java 8
 - sudo apt-get install oracle-java8-set-default
```

So I assume the same can be achieved with:

```
sudo add-apt-repository -y ppa:webupd8team/java 
sudo apt-get update -qq
sudo apt-get install oracle-java8-installer
sudo apt-get install oracle-java8-set-default
```

### Ubuntu 17.10 (Artful Aardvark)

```
sudo apt-get install r-cran-rjava openjdk-8-jdk
R CMD javareconf
```

Do not use `openjdk-9-jdk`.

### Ubuntu 18.4 (Bionic Beaver)

The `.travis.yml` file shows a Trusty install:

```
# - sudo apt install -qq oracle-java8-installer # Java 8
# - sudo apt install oracle-java8-set-default
```

On Bionic, I achieved the same with approx:

```
sudo add-apt-repository ppa:marutter/c2d4u3.5
sudo apt update
sudo apt grade
sudo apt install r-cran-rjava
sudo apt-get install openjdk-11-jdk
sudo R CMD javareconf
```

```
#sudo add-apt-repository -y ppa:webupd8team/java 
#sudo apt-get update -qq
#sudo apt-get install oracle-java8-installer
#sudo apt-get install oracle-java8-set-default
```


## Troubleshooting

### Error: `libjvm.so: cannot open shared object file: No such file or directory`

#### Random solution 1

Sometimes works:

For me, [this Stack Overflow post](https://stackoverflow.com/a/25932828) helped me out:

```
sudo mousepad /etc/ld.so.conf.d/java.conf
```

In that file put:

```
/usr/lib/jvm/java-8-oracle/jre/lib/amd64
/usr/lib/jvm/java-8-oracle/jre/lib/amd64/server
```

Save, close, restart R studio, fixed!

#### Random solution 2

Random notes:

Else [this Stack Overflow post may be helpful](https://stackoverflow.com/a/43466434):

Ruthlessly install all JDK stuff:

```
sudo apt-get install jdk-*
```

```
sudo R CMD javareconf
```

```
sudo R CMD javareconf -e
export LD_LIBRARY_PATH=$JAVA_LD_LIBRARY_PATH
sudo apt-get install r-cran-rjava
```

### BEAST2 cannot find Java

![BEAST2 cannot find Java](beast_cannot_find_java.png)


Download the Oracle Java SDK:

![](download_oracle_java_sdk.png)

Open the Oracle Java SDK with the package installer:

![](open_oracle_java_sdk.png)

Install the Oracle Java SDK with the package installer:

![](install_oracle_java_sdk.png)

Pick the right `java`:


```
sudo update-alternatives --config java
```

I picked:

```
There are 5 choices for the alternative java (providing /usr/bin/java).

  Selection    Path                                            Priority   Status
------------------------------------------------------------
  0            /usr/lib/jvm/java-9-openjdk-amd64/bin/java       1091      auto mode
  1            /usr/bin/gij-4.8                                 1048      manual mode
  2            /usr/bin/gij-5                                   1050      manual mode
* 3            /usr/bin/gij-6                                   1060      manual mode
  4            /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java   1081      manual mode
  5            /usr/lib/jvm/java-9-openjdk-amd64/bin/java       1091      manual mode

Press <enter> to keep the current choice[*], or type selection number: 
```

Reconfig:

```
sudo R CMD javareconf
```
# FAQ

 1. `beastier` in [academia](#academia) 
 2. `beastier` and [BEAST2](#BEAST2)
 3. `beastier` [development](#development) 
 4. `beastier` [technical](#technical)
 5. `beastier` [misc](#misc)

## [academia](#academia)

`beastier` in [academia](#academia).

### 1.1 How do I reference to this work?

Cite:

```
Bilderbeek, Richèl JC, and Rampal S. Etienne. "babette: BEAUti 2, BEAST 2 and Tracer for R." Methods in Ecology and Evolution (2018).
```

or

```
@article{bilderbeek2018babette,
  title={babette: BEAUti 2, BEAST 2 and Tracer for R},
  author={Bilderbeek, Richèl JC and Etienne, Rampal S},
  journal={Methods in Ecology and Evolution},
  year={2018},
  publisher={Wiley Online Library}
}
```

### 1.2 What are the FASTA files?

FASTA files `anthus_aco.fas` and `anthus_nd2.fas` from:
 
 * Van Els, Paul, and Heraldo V. Norambuena. "A revision of species limits in Neotropical pipits Anthus based on multilocus genetic and vocal data." Ibis.

Thanks to Paul van Els.

## [BEAST2](#BEAST2)

`beastier` and [BEAST2](#BEAST2).

### 2.1 How to install BEAST2?

```
Due to CRAN policy, beastier cannot install BEAST2.
As a workaround, the non-CRAN `beastierinstall` can be used.

To install BEAST2:

```
remotes::install_github("richelbilderbeek/beastierinstall")
beastierinstall::install_beast2()
```

This will download and extract BEAST2 to:

OS     |Full path
-------|----------------------------------
Linux  |`~/.local/share/beast`
macOS  |`~/.local/share/beast`
Windows|`C:/Users/<username>/Local/beast`
```

### 2.2 Which version of BEAUti do you use as a guideline?

Version 2.6.0, as can be found in the [install_beast2](https://github.com/ropensci/beastier/blob/master/R/install_beast2.R) function.

## [development](#development) 

`beastier` [development](#development).

### 3.1 How can I indicate a feature that I miss?

Submit an Issue.

### 3.2 How can I submit code?

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting code'

### 3.3 How can I submit a bug?

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting bugs' 

### 3.4 How can I indicate something else?

Submit an Issue. Or send an email to Richèl Bilderbeek.

## [technical](#technical)

`beastier` technical questions.

### 4.1 Why doesn't `beastier` support calling the Windows BEAST2.exe file?

The goal of `beastier` is to call BEAST2 from R scripts.
The Windows `BEAST2.exe` executable starts a graphical user interface.
An R script should be silent, without pop-ups. 
Therefore, calling the Windows `BEAST2.exe` executable is disallowed.

If this changes, you are encouraged to inform me, by either an Issue
or an email.

### 4.2 Installing Java under A

```
sudo apt-get install r-cran-rjava
R CMD javareconf
```

### 4.2 Installing Java under Bionic

The `.travis.yml` file shows a Trusty install:

```
 - sudo apt-get install -qq oracle-java8-installer # Java 8
 - sudo apt-get install oracle-java8-set-default
```

On Bionic, I assume the same can be achieved with:

```
sudo add-apt-repository -y ppa:webupd8team/java 
sudo apt-get update -qq
sudo apt-get install oracle-java8-installer
sudo apt-get install oracle-java8-set-default
```

### 4.3 Why doesn't `beastier` have 100% code coverage?

Because `beastier` cannot be fully tested for both
Linux and Windows on the same operating system.

Code coverage is measured by [codecov](https://codecov.io/gh/ropensci/beastier/tree/master/R) by the Travis CI continuous integration service.
Travis uses Linux. 
One can observe all missing code coverage is due to Windows-only functions.

Sure, there are also tests by the AppVeyor continuous integration service.
AppVeyor uses Windows. Would one observe that code coverage report, 
one would observe all missing code coverage is due to Linux-only functions.

## [misc](#misc)

`beastier` miscellaneous topics.

### 5.1 Why the logo?

Initially, the logo was a low-tech remake of Beast, from Marvel.
To prevent problems with Marvel, a different logo was picked.

The current logo shows a hippo, 'quite a formidable beast', also shown
intimidatingly big for the R logo. 
The hippo is drawn by Jose Scholte, who kindly allowed her work to
be used for free, by attribution.

### 5.2 How did you convert the fuzzy white background to one single color?

```
convert hippo.png -fuzz 15% -fill white -opaque white hippo_mono_background.png
convert hippo_mono_background.png -background white -alpha remove hippo_mono_background_2.png
```

