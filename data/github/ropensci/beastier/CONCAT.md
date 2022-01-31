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

---
title: "beastier demo"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{beastier demo}
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

![](beastier_logo.png)

This vignette demonstrates how to use `beastier`.

First, load the library:

```{r load_beastier, results='hide', warning=FALSE, error=FALSE, message=FALSE}
library(beastier)
```

Also, we'll load the `testthat` library, to verify the statements in this vignette:

```{r load_testthat}
library(testthat)
```

To run BEAST2, we need to create a BEAST2 options structure`. We will use
a supplied BEAST2 XML file. For the rest, we'll use the default options:

```{r}
beast2_options <- create_beast2_options(
  input_filename = get_beastier_path("2_4.xml")
)
names(beast2_options)
```


Before running BEAST2, the BEAST2 input file must exist,
and we expect no output file to be created just yet:

```{r}
expect_true(file.exists(beast2_options$input_filename))
expect_false(file.exists(beast2_options$output_state_filename))
```

We can run `beastier` now, if BEAST2 is installed. Because BEAST2 needs
to be installed by the user, this vignette checks if it is installed in every step:

```{r}
if (is_beast2_installed()) {
  output <- run_beast2_from_options(beast2_options)
}
```

If `beastier` has run BEAST2, the BEAST2 output can be shown:

```{r}
if (is_beast2_installed()) {
  print(output)
}
```

If `beastier` has run BEAST2, the MCMC's final state will be saved to a file: 

```{r}
if (is_beast2_installed()) {
  expect_true(file.exists(beast2_options$output_state_filename))
  file.remove(beast2_options$output_state_filename)
}
```

This final state can be used to continue the run.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_default_beast2_download_url.R
\name{get_default_beast2_download_url}
\alias{get_default_beast2_download_url}
\title{Get the default BEAST2 download URL,
which depends on the operating system}
\usage{
get_default_beast2_download_url(
  beast2_version = beastier::get_default_beast2_version(),
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{beast2_version}{the version of BEAST2. By
default, this is the version as returned by
\link{get_default_beast2_version}}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
the URL where BEAST2 can be downloaded from
}
\description{
Get the default BEAST2 download URL,
which depends on the operating system
}
\examples{
get_default_beast2_download_url()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_state_output_file_folder.R
\name{create_beast2_state_output_file_folder}
\alias{create_beast2_state_output_file_folder}
\title{Create the folder where the BEAST2 state output file will be created}
\usage{
create_beast2_state_output_file_folder(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
nothing
}
\description{
Create the folder where the BEAST2 state output file will be created
}
\examples{
beast2_options <- create_beast2_options()
create_beast2_state_output_file_folder(beast2_options)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_minimal_run.R
\name{do_minimal_run}
\alias{do_minimal_run}
\title{Do a minimal BEAST2 run}
\usage{
do_minimal_run()
}
\value{
The text sent to \code{STDOUT} and \code{STDERR}.
  It will create the files with name \code{output_state_filename}
}
\description{
To achieve this, \link{run_beast2_from_options} is called.
}
\examples{
if (is_beast2_installed() && is_on_ci()) {
  do_minimal_run()
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_beast2_input_lines.R
\name{are_beast2_input_lines_fast}
\alias{are_beast2_input_lines_fast}
\title{Would these lines of text, when written to a file,
  result in a valid BEAST2 input file?}
\usage{
are_beast2_input_lines_fast(lines)
}
\arguments{
\item{lines}{lines of text}
}
\value{
TRUE if the text is valid, FALSE if not
}
\description{
Would these lines of text, when written to a file,
  result in a valid BEAST2 input file?
}
\examples{

beast2_filename <- get_beastier_path("anthus_2_4.xml")
text <- readLines(beast2_filename)

# TRUE
are_beast2_input_lines_fast(text)
}
\seealso{
Use \code{\link{is_beast2_input_file}} to check a file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beast2_options.R
\name{check_beast2_options}
\alias{check_beast2_options}
\title{Check if the \code{beast2_options} is a valid BEAST2 options object.}
\usage{
check_beast2_options(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
nothing
Will \code{stop} if the BEAST2 option object is invalid
}
\description{
Calls \code{stop} if the BEAST2 option object is invalid
}
\examples{
check_beast2_options(create_beast2_options())
}
\seealso{
Use \link{create_beast2_options} to create a valid
  BEAST2 options object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_beast2_options.R
\name{print_beast2_options}
\alias{print_beast2_options}
\title{Pretty-print a `beast2_options`}
\usage{
print_beast2_options(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
Nothing. Will display the `beast2_options` using \link{cat}.
}
\description{
Pretty-print a `beast2_options`
}
\examples{
print_beast2_options(create_beast2_options())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_alignment_ids.R
\name{get_alignment_ids_from_xml_filename}
\alias{get_alignment_ids_from_xml_filename}
\title{Get the alignment ID from a file with one alignment}
\usage{
get_alignment_ids_from_xml_filename(xml_filename)
}
\arguments{
\item{xml_filename}{name of a BEAST2 XML input filename}
}
\value{
one or more alignment IDs
}
\description{
Get the alignment ID from a file with one alignment
}
\examples{
# test_output_0
get_alignment_ids_from_xml_filename(get_beastier_path("2_4.xml"))
# c("anthus_aco","anthus_nd2")
get_alignment_ids_from_xml_filename(get_beastier_path("anthus_15_15.xml"))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beast2_version.R
\name{get_beast2_version}
\alias{get_beast2_version}
\title{Get the BEAST2 version}
\usage{
get_beast2_version(beast2_path = get_default_beast2_path())
}
\arguments{
\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
the BEAST2 version
}
\description{
Get the BEAST2 version
}
\examples{
if (is_beast2_installed() && is_on_ci()) {
  get_beast2_version()
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_n_threads.R
\name{check_n_threads}
\alias{check_n_threads}
\title{Check if the input is a valid number of threads.}
\usage{
check_n_threads(n_threads)
}
\arguments{
\item{n_threads}{the number of computational threads to use.
Use \link{NA} to use the BEAST2 default of 1.}
}
\value{
Nothing.
Will \link{stop} if the number of threads in invalid
}
\description{
Will \link{stop} if not.
}
\examples{
# Can have 1 or more threads
check_n_threads(1)
check_n_threads(2)
# Can have NA threads
check_n_threads(NA)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_beast2.R
\name{run_beast2}
\alias{run_beast2}
\title{Run BEAST2}
\usage{
run_beast2(
  input_filename,
  output_log_filename = "output_log_filename_is_deprecated",
  output_trees_filenames = "output_trees_filenames_is_deprecated",
  output_state_filename = create_temp_state_filename(),
  rng_seed = NA,
  n_threads = NA,
  use_beagle = FALSE,
  overwrite = TRUE,
  beast2_working_dir = "beast2_working_dir_is_deprecated",
  beast2_path = get_default_beast2_path(),
  verbose = FALSE
)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}

\item{output_log_filename}{name of the .log file to create}

\item{output_trees_filenames}{one or more names for .trees file to create.
There will be one .trees file created per alignment in the input
file. The number of alignments must equal the number of .trees
filenames, else an error is thrown. Alignments are sorted alphabetically
by their IDs}

\item{output_state_filename}{name of the \code{.xml.state} file to create.
Use \link{create_temp_state_filename} to create a temporary
filename with that extension.}

\item{rng_seed}{the random number generator seed of the BEAST2 run.
Must be a non-zero positive integer value or \link{NA}.
If \code{rng_seed} is \link{NA}, BEAST2 will pick a random seed}

\item{n_threads}{the number of computational threads to use.
Use \link{NA} to use the BEAST2 default of 1.}

\item{use_beagle}{use BEAGLE if present}

\item{overwrite}{if TRUE: overwrite the \code{.log}
and \code{.trees} files if one of these exists.
If FALSE, BEAST2 will not be started if
\itemize{
  \item{the \code{.log} file exists}
  \item{the \code{.trees} files exist}
  \item{the \code{.log} file created by BEAST2 exists}
  \item{the \code{.trees} files created by BEAST2 exist}
}}

\item{beast2_working_dir}{a folder where BEAST2 can work in
isolation.
For each BEAST2 run, a new subfolder is created in that folder.
Within this folder, BEAST2 is allowed to create all of its output files,
without the risk of overwriting existing ones, allowing
BEAST2 to run in multiple parallel processes.}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
The text sent to \code{STDOUT} and \code{STDERR}.
  It will create the file with name \code{output_state_filenames}
}
\description{
Run BEAST2
}
\examples{

if (is_beast2_installed() && is_on_ci()) {

  output_state_filename <- create_temp_state_filename()

  run_beast2(
    input_filename = get_beastier_path("2_4.xml"),
    output_state_filename = output_state_filename

  )
  file.remove(output_state_filename)
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_trees_filenames.R
\name{get_trees_filenames}
\alias{get_trees_filenames}
\title{Get the .trees filenames that BEAST2 will produce}
\usage{
get_trees_filenames(input_filename)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}
}
\value{
character vector with the names of the .trees files that BEAST2
  will produce
}
\description{
Get the .trees filenames that BEAST2 will produce
}
\examples{
get_trees_filenames(get_beastier_path("2_4.xml"))
get_trees_filenames(get_beastier_path("anthus_2_4.xml"))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_validate_cmd.R
\name{create_beast2_validate_cmd_jar}
\alias{create_beast2_validate_cmd_jar}
\title{Creates the terminal command to validate a BEAST2 input file
using a call to the \code{launcher.jar} file}
\usage{
create_beast2_validate_cmd_jar(
  input_filename,
  beast2_jar_path = get_default_beast2_jar_path()
)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}

\item{beast2_jar_path}{name of the BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
a character vector, of which the first element
  is the command (\code{java}, in this case),
  and the others are arguments (\code{-jar}, in this case, followed
  by more arguments.
}
\description{
Creates the terminal command to validate a BEAST2 input file
using a call to the \code{launcher.jar} file
}
\examples{
  if (is_beast2_installed() && is_on_ci()) {
    create_beast2_validate_cmd_jar(
      input_filename = "input.xml"
    )
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_java_version.R
\name{get_java_version}
\alias{get_java_version}
\title{Get the Java version}
\usage{
get_java_version()
}
\value{
the Java version
}
\description{
Get the Java version
}
\examples{

if (is_beast2_installed() && is_on_ci()) {
  get_java_version()
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_duplicate_param_ids.R
\name{get_duplicate_param_ids}
\alias{get_duplicate_param_ids}
\title{Find duplicate \code{RealParameter} IDs}
\usage{
get_duplicate_param_ids(text)
}
\arguments{
\item{text}{the XML as text}
}
\value{
a vector of duplicate IDs, will be empty if all IDs are unique
}
\description{
Find duplicate \code{RealParameter} IDs
}
\examples{
  line_1 <- "<parameter id=\"RealParameter.1\" ...</parameter>"
  line_2 <- "<parameter id=\"RealParameter.2\" ...</parameter>"
  testit::assert(
    length(get_duplicate_param_ids(c(line_1, line_2))) == 0)
  testit::assert(
    get_duplicate_param_ids(
    c(line_1, line_1)) == c("RealParameter.1")
  )
  testit::assert(
    get_duplicate_param_ids(
    c(line_2, line_2)) == c("RealParameter.2")
  )
}
\seealso{
to see if all IDs are unique, use \code{\link{has_unique_ids}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gives_beast2_warning.R
\name{gives_beast2_warning}
\alias{gives_beast2_warning}
\title{Determines if BEAST2 issues a warning
when using the BEAST2 XML input file}
\usage{
gives_beast2_warning(
  filename,
  verbose = FALSE,
  beast2_path = beastier::get_default_beast2_path()
)
}
\arguments{
\item{filename}{name of the BEAST2 XML input file}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
TRUE if the file produces a BEAST2 warning, FALSE if not
}
\description{
Determines if BEAST2 issues a warning
when using the BEAST2 XML input file
}
\examples{
if (is_beast2_installed() &&
  is_on_ci() &&
  rappdirs::app_dir()$os == "unix") {

  # This file is OK for BEAST2, no warning, returns FALSE
  gives_beast2_warning(filename = get_beastier_path("2_4.xml"))

  # BEAST2 will give a warning on this file, returns TRUE
  gives_beast2_warning(
    filename = get_beastier_path("beast2_warning.xml")
  )
}
}
\seealso{
Use \code{\link{is_beast2_input_file}} to check if a file is a
  valid BEAST2 input file.
  Use \code{\link{are_beast2_input_lines}} to check if the text (for
  example, as loaded from a file) to be valid BEAST2 input.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_beast2_input_file.R
\name{is_beast2_input_file}
\alias{is_beast2_input_file}
\title{Is a file a valid BEAST2 input file?}
\usage{
is_beast2_input_file(
  filename,
  show_warnings = FALSE,
  verbose = FALSE,
  beast2_path = get_default_beast2_path()
)
}
\arguments{
\item{filename}{name of the BEAST2 XML input file}

\item{show_warnings}{if TRUE, warnings will shown}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
TRUE if the file is valid, FALSE if not
}
\description{
Is a file a valid BEAST2 input file?
}
\note{
this function only works on standard BEAST2 input files:
    if a BEAST2 input file is modified to use a certain BEAST2 package,
    this function will label it as an invalid file
}
\examples{

if (is_beast2_installed() && is_on_ci()) {

  filename <- get_beastier_path("anthus_2_4.xml")
  # TRUE, this is a BEAST2 input file
  is_beast2_input_file(filename)

  filename <- get_beastier_path("beast2_example_output.log")
  # FALSE, this is not a BEAST2 input file,
  # it is a BEAST2 output log file insteaf
  is_beast2_input_file(filename)
}
}
\seealso{
Use \code{\link{are_beast2_input_lines}} to check the lines
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_random_phylogeny.R
\name{create_random_phylogeny}
\alias{create_random_phylogeny}
\title{Create a random phylogeny}
\usage{
create_random_phylogeny(n_taxa, taxa_name_ext = "")
}
\arguments{
\item{n_taxa}{The number of taxa}

\item{taxa_name_ext}{the extension of the taxa names}
}
\value{
a phylogeny of class `phylo` (which is part of the `ape` package)
}
\description{
Create a random phylogeny
}
\examples{
create_random_phylogeny(n_taxa = 6)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beastier_folder.R
\name{get_beastier_folder}
\alias{get_beastier_folder}
\title{Get the path to the \link{beastier} temporary files folder}
\usage{
get_beastier_folder()
}
\value{
the path to the \link{beastier} temporary files folder.
}
\description{
Get the path to the \link{beastier} temporary files folder.
}
\examples{
get_beastier_folder()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beast2_example_filename.R
\name{get_beast2_example_filename}
\alias{get_beast2_example_filename}
\title{Get the full path of a BEAST2 example file}
\usage{
get_beast2_example_filename(
  filename,
  beast2_folder = get_default_beast2_folder()
)
}
\arguments{
\item{filename}{name of the BEAST2 example file. This should exclude
the full path; this function exists to add that full path}

\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable
is installed: the BEAST2 executable is in a subfolder.
Use \link{get_default_beast2_folder} to get the default BEAST2
folder.
Use \link{get_default_beast2_bin_path} to get the full path to
the default BEAST2 executable.}
}
\value{
the full path of a BEAST2 example file,
will \link{stop} if the filename is not a BEAST2 example file
}
\description{
Will \link{stop} if the filename is not a BEAST2 example file
}
\examples{
if (is_beast2_installed()) {
  get_beast2_example_filename("testJukesCantor.xml")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_tracelog_folder.R
\name{create_beast2_tracelog_folder}
\alias{create_beast2_tracelog_folder}
\title{Internal function}
\usage{
create_beast2_tracelog_folder(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\description{
Create the folder for the BEAST2 tracelog file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_path.R
\name{get_beastier_path}
\alias{get_beastier_path}
\title{Get the full path of a file in the \code{inst/extdata} folder}
\usage{
get_beastier_path(filename)
}
\arguments{
\item{filename}{the file's name, without the path}
}
\value{
the full path to the filename. Will \code{stop} if the file
  is absent in the \code{inst/extdata} folder
}
\description{
Get the full path of a file in the \code{inst/extdata} folder
}
\examples{
get_beastier_path("beast2_example_output.log")
get_beastier_path("beast2_example_output.trees")
get_beastier_path("beast2_example_output.xml")
get_beastier_path("beast2_example_output.xml.state")
}
\seealso{
for more files, use \code{\link{get_beastier_paths}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_screenlog_folder.R
\name{create_beast2_screenlog_folder}
\alias{create_beast2_screenlog_folder}
\title{Internal function}
\usage{
create_beast2_screenlog_folder(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\description{
Create the folder for the BEAST2 screenlog file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_file_if_present.R
\name{remove_file_if_present}
\alias{remove_file_if_present}
\title{Remove a file if it is present,
will do nothing if it is not.}
\usage{
remove_file_if_present(filename)
}
\arguments{
\item{filename}{name of a file}
}
\value{
Nothing. Will remove the file if it is presented,
will do nothing if it is not.
}
\description{
Remove a file if it is present,
will do nothing if it is not.
}
\examples{
filename <- tempfile()
file.create(filename)
remove_file_if_present(filename)
remove_file_if_present(filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beast2_options.R
\name{check_beast2_options_names}
\alias{check_beast2_options_names}
\title{Check if the \code{beast2_options}, which is a list,
has all the elements needed.}
\usage{
check_beast2_options_names(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
nothing
}
\description{
Calls \code{stop} if not.
}
\seealso{
Use \link{check_beast2_options} to check
  the entire \code{beast2_options} object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_on_ci.R
\name{is_on_travis}
\alias{is_on_travis}
\title{Determines if the environment is Travis CI}
\usage{
is_on_travis()
}
\value{
TRUE if run on Travis CI, FALSE otherwise
}
\description{
Determines if the environment is Travis CI
}
\examples{
  if (is_on_ci()) {
    message("Running on Travis CI")
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_mcbette_beast2_options.R
\name{create_mcbette_beast2_options}
\alias{create_mcbette_beast2_options}
\title{Create a `beast2_options` structure for the `mcbette` R package}
\usage{
create_mcbette_beast2_options(
  input_filename = beastier::create_temp_input_filename(),
  output_state_filename = beastier::create_temp_state_filename(),
  rng_seed = NA,
  n_threads = NA,
  use_beagle = FALSE,
  overwrite = TRUE,
  beast2_bin_path = beastier::get_default_beast2_bin_path(),
  verbose = FALSE
)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}

\item{output_state_filename}{name of the \code{.xml.state} file to create.
Use \link{create_temp_state_filename} to create a temporary
filename with that extension.}

\item{rng_seed}{the random number generator seed of the BEAST2 run.
Must be a non-zero positive integer value or \link{NA}.
If \code{rng_seed} is \link{NA}, BEAST2 will pick a random seed}

\item{n_threads}{the number of computational threads to use.
Use \link{NA} to use the BEAST2 default of 1.}

\item{use_beagle}{use BEAGLE if present}

\item{overwrite}{if TRUE: overwrite the \code{.log}
and \code{.trees} files if one of these exists.
If FALSE, BEAST2 will not be started if
\itemize{
  \item{the \code{.log} file exists}
  \item{the \code{.trees} files exist}
  \item{the \code{.log} file created by BEAST2 exists}
  \item{the \code{.trees} files created by BEAST2 exist}
}}

\item{beast2_bin_path}{name of the BEAST2 binary file
(usually simply \code{beast}).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
a `beast2_options` structure suitable to be used 
by the `mcbette` R package,
which is a \link{list} of all function arguments,
of which all elements are checked (by \link{check_beast2_options})
}
\description{
Create a `beast2_options` structure to be used for the
`mcbette` R package, which is a package that allows one to do 
model comparison.
The generated filenames indicating `mcbette` usage,
as well as the correct BEAST2 binary/executable type
}
\examples{
create_mcbette_beast2_options()
}
\seealso{
to create a regular (that is, not intended
for model comparison) BEAST2 options structure,
use  \link{create_beast2_options}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_params_doc.R
\name{default_params_doc}
\alias{default_params_doc}
\title{This function does nothing. It is intended to inherit is parameters'
documentation.}
\usage{
default_params_doc(
  beast2_bin_path,
  beast2_folder,
  beast2_jar_path,
  beast2_options,
  beast2_optionses,
  beast2_path,
  beast2_version,
  beast2_working_dir,
  beastier_folder,
  clock_model,
  clock_models,
  crown_age,
  crown_ages,
  fasta_filename,
  fasta_filenames,
  fixed_crown_age,
  fixed_crown_ages,
  initial_phylogenies,
  input_filename,
  mcmc,
  misc_options,
  n_taxa,
  n_threads,
  os,
  output_filename,
  output_log_filename,
  output_state_filename,
  output_trees_filenames,
  overwrite,
  rename_fun,
  rng_seed,
  sequence_length,
  site_model,
  site_models,
  tree_prior,
  tree_priors,
  use_beagle,
  verbose
)
}
\arguments{
\item{beast2_bin_path}{name of the BEAST2 binary file
(usually simply \code{beast}).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path}

\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable
is installed: the BEAST2 executable is in a subfolder.
Use \link{get_default_beast2_folder} to get the default BEAST2
folder.
Use \link{get_default_beast2_bin_path} to get the full path to
the default BEAST2 executable.}

\item{beast2_jar_path}{name of the BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}

\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}

\item{beast2_optionses}{list of one or more \code{beast2_options}
structures,
as can be created by \link{create_beast2_options}.
Use of reduplicated plural to achieve difference with
\code{beast2_options}}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}

\item{beast2_version}{the version of BEAST2. By
default, this is the version as returned by
\link{get_default_beast2_version}}

\item{beast2_working_dir}{a folder where BEAST2 can work in
isolation.
For each BEAST2 run, a new subfolder is created in that folder.
Within this folder, BEAST2 is allowed to create all of its output files,
without the risk of overwriting existing ones, allowing
BEAST2 to run in multiple parallel processes.}

\item{beastier_folder}{the path to
the \link{beastier} temporary files folder}

\item{clock_model}{a \code{beautier} clock model}

\item{clock_models}{a list of one or more \code{beautier} clock models}

\item{crown_age}{the crown age of the phylogeny}

\item{crown_ages}{the crown ages of the phylogenies. Set to NA
if the crown age needs to be estimated}

\item{fasta_filename}{a FASTA filename.}

\item{fasta_filenames}{One or more FASTA filenames.}

\item{fixed_crown_age}{determines if the phylogeny's crown age is
fixed. If FALSE, crown age is estimated by BEAST2. If TRUE,
the crown age is fixed to the crown age
of the initial phylogeny.}

\item{fixed_crown_ages}{one or more booleans to determine if the
phylogenies' crown ages are fixed.
If FALSE, crown age is estimated by BEAST2. If TRUE,
the crown age is fixed to the crown age
of the initial phylogeny.}

\item{initial_phylogenies}{one or more MCMC chain's initial phylogenies.
Each one set to NA will result in BEAST2 using a random phylogeny. Else
the phylogeny is assumed to be of class \code{ape::phylo}.}

\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}

\item{mcmc}{one \code{beautier} MCMC}

\item{misc_options}{one \code{beautier} misc_options object}

\item{n_taxa}{The number of taxa}

\item{n_threads}{the number of computational threads to use.
Use \link{NA} to use the BEAST2 default of 1.}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}

\item{output_filename}{Name of the XML parameter file created by this
function. BEAST2 uses this file as input.}

\item{output_log_filename}{name of the .log file to create}

\item{output_state_filename}{name of the \code{.xml.state} file to create.
Use \link{create_temp_state_filename} to create a temporary
filename with that extension.}

\item{output_trees_filenames}{one or more names for .trees file to create.
There will be one .trees file created per alignment in the input
file. The number of alignments must equal the number of .trees
filenames, else an error is thrown. Alignments are sorted alphabetically
by their IDs}

\item{overwrite}{if TRUE: overwrite the \code{.log}
and \code{.trees} files if one of these exists.
If FALSE, BEAST2 will not be started if
\itemize{
  \item{the \code{.log} file exists}
  \item{the \code{.trees} files exist}
  \item{the \code{.log} file created by BEAST2 exists}
  \item{the \code{.trees} files created by BEAST2 exist}
}}

\item{rename_fun}{a function to rename a filename,
as can be checked by \link{check_rename_fun}. This function should
have one argument, which will be a filename or \link{NA}. The
function should \link{return} one filename (when passed one filename) or
one \link{NA} (when passed one \link{NA}).
Example rename functions are:
\itemize{
  \item \link[beautier]{get_remove_dir_fun} get a function
    that removes the directory
    paths from the filenames, in effect turning these into local files
  \item \link[beautier]{get_replace_dir_fun} get a function
    that replaces the directory
    paths from the filenames
  \item \link[beautier]{get_remove_hex_fun} get a function that
    removes the hex string from filenames.
    For example, \code{tracelog_82c1a522040.log} becomes \code{tracelog.log}
}}

\item{rng_seed}{the random number generator seed of the BEAST2 run.
Must be a non-zero positive integer value or \link{NA}.
If \code{rng_seed} is \link{NA}, BEAST2 will pick a random seed}

\item{sequence_length}{a DNA sequence length, in base pairs}

\item{site_model}{a \code{beautier} site model}

\item{site_models}{one or more \code{beautier} site models}

\item{tree_prior}{a \code{beautier} tree prior}

\item{tree_priors}{one or more \code{beautier} tree priors}

\item{use_beagle}{use BEAGLE if present}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
Nothing. This is an internal function that does nothing
}
\description{
This function does nothing. It is intended to inherit is parameters'
documentation.
}
\note{
This is an internal function, so it should be marked with
\code{@noRd}. This is not done, as this will disallow all
functions to find the documentation parameters
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_options.R
\name{create_beast2_options}
\alias{create_beast2_options}
\title{Function to create a set of BEAST2 options.}
\usage{
create_beast2_options(
  input_filename = create_temp_input_filename(),
  output_state_filename = create_temp_state_filename(),
  rng_seed = NA,
  n_threads = NA,
  use_beagle = FALSE,
  overwrite = TRUE,
  beast2_path = get_default_beast2_path(),
  verbose = FALSE,
  output_log_filename = "deprecated",
  output_trees_filenames = "deprecated",
  beast2_working_dir = "deprecated"
)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}

\item{output_state_filename}{name of the \code{.xml.state} file to create.
Use \link{create_temp_state_filename} to create a temporary
filename with that extension.}

\item{rng_seed}{the random number generator seed of the BEAST2 run.
Must be a non-zero positive integer value or \link{NA}.
If \code{rng_seed} is \link{NA}, BEAST2 will pick a random seed}

\item{n_threads}{the number of computational threads to use.
Use \link{NA} to use the BEAST2 default of 1.}

\item{use_beagle}{use BEAGLE if present}

\item{overwrite}{if TRUE: overwrite the \code{.log}
and \code{.trees} files if one of these exists.
If FALSE, BEAST2 will not be started if
\itemize{
  \item{the \code{.log} file exists}
  \item{the \code{.trees} files exist}
  \item{the \code{.log} file created by BEAST2 exists}
  \item{the \code{.trees} files created by BEAST2 exist}
}}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}

\item{output_log_filename}{name of the .log file to create}

\item{output_trees_filenames}{one or more names for .trees file to create.
There will be one .trees file created per alignment in the input
file. The number of alignments must equal the number of .trees
filenames, else an error is thrown. Alignments are sorted alphabetically
by their IDs}

\item{beast2_working_dir}{a folder where BEAST2 can work in
isolation.
For each BEAST2 run, a new subfolder is created in that folder.
Within this folder, BEAST2 is allowed to create all of its output files,
without the risk of overwriting existing ones, allowing
BEAST2 to run in multiple parallel processes.}
}
\value{
a BEAST2 options structure, 
which is a \link{list} of all function arguments,
of which all elements are checked (by \link{check_beast2_options})
}
\description{
These BEAST2 options are the R equivalent of the command-line options.
}
\examples{
beast2_options <- create_beast2_options()
check_beast2_options(beast2_options)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_lines.R
\name{save_lines}
\alias{save_lines}
\title{Save text (a container of strings) to a file}
\usage{
save_lines(filename, lines)
}
\arguments{
\item{filename}{filename of the file to have the text written to}

\item{lines}{lines of text to be written to file}
}
\value{
Nothing. Will save the lines to file
}
\description{
Save text (a container of strings) to a file
}
\examples{
text <- c("hello", "world")
filename <- get_beastier_tempfilename()
save_lines(filename = filename, lines = text)
file.remove(filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beautier_tempfolder.R
\name{create_beautier_tempfolder}
\alias{create_beautier_tempfolder}
\title{Create the temporary folder as used by \link[beautier]{beautier}}
\usage{
create_beautier_tempfolder()
}
\value{
nothing
}
\description{
Create the temporary folder as used by \link[beautier]{beautier}
}
\examples{
create_beautier_tempfolder()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beast2_optionses.R
\name{check_beast2_optionses}
\alias{check_beast2_optionses}
\title{Check if the \code{beast2_options} is a valid BEAST2 options object.}
\usage{
check_beast2_optionses(beast2_optionses)
}
\arguments{
\item{beast2_optionses}{list of one or more \code{beast2_options}
structures,
as can be created by \link{create_beast2_options}.
Use of reduplicated plural to achieve difference with
\code{beast2_options}}
}
\value{
Nothing. 
Will \code{stop} if the BEAST2 option object is invalid
}
\description{
Calls \code{stop} if the BEAST2 option object is invalid
}
\examples{
 check_beast2_optionses(list(create_beast2_options()))
}
\seealso{
Use \link{create_beast2_options} to create a valid
  BEAST2 options object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beast2.R
\name{check_beast2}
\alias{check_beast2}
\title{Check if \code{BEAST2} is installed properly.}
\usage{
check_beast2(beast2_path = beastier::get_default_beast2_path())
}
\arguments{
\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
nothing
Will \link{stop} if BEAST2 is improperly installed
}
\description{
Calls \link{stop} if BEAST2 is improperly installed
}
\examples{
if (is_beast2_installed()) {
  check_beast2()
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_quotes_if_has_spaces.R
\name{add_quotes_if_has_spaces}
\alias{add_quotes_if_has_spaces}
\title{Add quotes around the string if it contains spaces.}
\usage{
add_quotes_if_has_spaces(filename)
}
\arguments{
\item{filename}{a filename}
}
\value{
a filename. If the filename did not contain spaces,
it is returned as-is. If the filename did contain spaces,
the filename is surrounded by quotes
}
\description{
Add quotes around the string if it contains spaces.
Does nothing if the string contains no spaces.
This is used for filenames
}
\examples{
add_quotes_if_has_spaces("x")
add_quotes_if_has_spaces("a b")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_treelog_folder.R
\name{create_beast2_treelog_folder}
\alias{create_beast2_treelog_folder}
\title{Internal function}
\usage{
create_beast2_treelog_folder(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\description{
Create the folder for the BEAST2 treelog file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beast2_options.R
\name{check_beast2_options_filenames_differ}
\alias{check_beast2_options_filenames_differ}
\title{Check if the filenames in \code{beast2_options} differ}
\usage{
check_beast2_options_filenames_differ(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
nothing
}
\description{
Calls \code{stop} if not.
}
\seealso{
Use \link{check_beast2_options} to check
  the entire \code{beast2_options} object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_can_create_screenlog_file.R
\name{check_can_create_screenlog_file}
\alias{check_can_create_screenlog_file}
\title{Internal function}
\usage{
check_can_create_screenlog_file(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
Nothing.
Will \link{stop} if the MCMC's screenlog file cannot be created.
}
\description{
Check if the MCMC's screenlog file can be created.
Will \link{stop} if not
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_default_java_path.R
\name{get_default_java_path}
\alias{get_default_java_path}
\title{Obtains the default path to the Java executable}
\usage{
get_default_java_path(os = rappdirs::app_dir()$os)
}
\arguments{
\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
the default path to the Java executable
}
\description{
Obtains the default path to the Java executable
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_file_folder.R
\name{create_beast2_input_file_folder}
\alias{create_beast2_input_file_folder}
\title{Create the folder where the BEAST2 input file will be created}
\usage{
create_beast2_input_file_folder(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
nothing
}
\description{
Create the folder where the BEAST2 input file will be created
}
\examples{
beast2_options <- create_beast2_options()
create_beast2_input_file_folder(beast2_options)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename_beast2_options_filenames.R
\name{rename_beast2_options_filenames}
\alias{rename_beast2_options_filenames}
\title{Rename the filenames in the BEAST2 options}
\usage{
rename_beast2_options_filenames(beast2_options, rename_fun)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}

\item{rename_fun}{a function to rename a filename,
as can be checked by \link{check_rename_fun}. This function should
have one argument, which will be a filename or \link{NA}. The
function should \link{return} one filename (when passed one filename) or
one \link{NA} (when passed one \link{NA}).
Example rename functions are:
\itemize{
  \item \link[beautier]{get_remove_dir_fun} get a function
    that removes the directory
    paths from the filenames, in effect turning these into local files
  \item \link[beautier]{get_replace_dir_fun} get a function
    that replaces the directory
    paths from the filenames
  \item \link[beautier]{get_remove_hex_fun} get a function that
    removes the hex string from filenames.
    For example, \code{tracelog_82c1a522040.log} becomes \code{tracelog.log}
}}
}
\value{
a `beast2_options` with the filenames it contains renamed
}
\description{
Rename the filenames in the BEAST2 options
}
\examples{
# beast2_options with local filenames
beast2_options <- create_beast2_options(
  input_filename = "my.fas",
  output_state_filename = "my_state.xml.state"
)
# Rename filenames to be in /my/new/folder
rename_beast2_options_filenames(
  beast2_options = beast2_options,
  rename_fun = beautier::get_replace_dir_fun("/my/new/folder")
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_version_cmd.R
\name{create_beast2_version_cmd_bin}
\alias{create_beast2_version_cmd_bin}
\title{Creates the terminal command to version a BEAST2 input file
using a call to the \code{launcher.jar} file}
\usage{
create_beast2_version_cmd_bin(beast2_bin_path = get_default_beast2_bin_path())
}
\arguments{
\item{beast2_bin_path}{name of the BEAST2 binary file
(usually simply \code{beast}).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path}
}
\value{
a character vector, of which the first element
  is the command (\code{java}, in this case),
  and the others are arguments (\code{-jar}, in this case, followed
  by more arguments.
}
\description{
Creates the terminal command to version a BEAST2 input file
using a call to the \code{launcher.jar} file
}
\examples{
  if (is_beast2_installed() && is_on_ci()) {
    create_beast2_version_cmd_bin()
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/extract_screenlog_filename_from_beast2_input_file.R
\name{extract_screenlog_filename_from_beast2_input_file}
\alias{extract_screenlog_filename_from_beast2_input_file}
\title{Internal function to extract the screenlog filename for a BEAST2 input file}
\usage{
extract_screenlog_filename_from_beast2_input_file(input_filename)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}
}
\value{
the screenlog filename for a BEAST2 input file
}
\description{
Extract the screenlog filename from a BEAST2 input file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_bin_path.R
\name{is_win_bin_path}
\alias{is_win_bin_path}
\title{Is the path a path to the BEAST2 binary file?
Does not check if the file at that path is present}
\usage{
is_win_bin_path(path)
}
\arguments{
\item{path}{a string to a path}
}
\value{
TRUE if the path is a path to a BEAST2 binary file
}
\description{
Is the path a path to the BEAST2 binary file?
Does not check if the file at that path is present
}
\examples{
# TRUE
is_win_bin_path("BEAST.exe")
# FALSE
is_win_bin_path("beast")
is_win_bin_path("launcher.jar")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_beast2_installed.R
\name{is_beast2_installed}
\alias{is_beast2_installed}
\title{Checks if BEAST2 is installed}
\usage{
is_beast2_installed(
  folder_name = get_default_beast2_folder(),
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{folder_name}{name of the folder where the BEAST2 files are put.
The name of the BEAST2 binary file will be at
\code{[folder_name]/beast/bin/beast}
The name of the BEAST2 jar file will be at
\code{[folder_name]/beast/lib/launcher.jar}}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
TRUE if BEAST2 is installed
}
\description{
Checks if BEAST2 is installed
}
\examples{
is_beast2_installed()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_on_ci.R
\name{is_on_ci}
\alias{is_on_ci}
\title{Determines if the environment is a continuous integration service}
\usage{
is_on_ci()
}
\value{
TRUE if run on AppVeyor or Travis CI, FALSE otherwise
}
\description{
Determines if the environment is a continuous integration service
}
\examples{
  if (is_on_ci()) {
    message("Running on a continuous integration service")
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beast2_path.R
\name{check_beast2_path}
\alias{check_beast2_path}
\title{Checks the BEAST2 \code{.jar} path.
Will stop if there is a problem with the BEAST2 \code{.jar} path.}
\usage{
check_beast2_path(beast2_path)
}
\arguments{
\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
nothing.
  Will call \code{\link{stop}} if the BEAST2 \code{.jar} path has a problem
}
\description{
Checks the BEAST2 \code{.jar} path.
Will stop if there is a problem with the BEAST2 \code{.jar} path.
}
\examples{
if (is_beast2_installed()) {
  beast2_path <- get_default_beast2_jar_path()
  check_beast2_path(beast2_path)
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_version_cmd.R
\name{create_beast2_version_cmd}
\alias{create_beast2_version_cmd}
\title{Creates the terminal command to version a BEAST2 input file}
\usage{
create_beast2_version_cmd(beast2_path = beastier::get_default_beast2_path())
}
\arguments{
\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
a character vector, of which the first element
  is the command (\code{java}, in this case),
  and the others are arguments (\code{-jar}, in this case, followed
  by more arguments.
}
\description{
Creates the terminal command to version a BEAST2 input file
}
\examples{
if (is_beast2_installed()) {
  create_beast2_version_cmd()
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_run_cmd_from_options.R
\name{create_beast2_run_cmd_from_options}
\alias{create_beast2_run_cmd_from_options}
\title{Creates the terminal command to run BEAST2 from a \code{beast2_options}}
\usage{
create_beast2_run_cmd_from_options(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
a character vector with the command and
  arguments to call BEAST2
}
\description{
Creates the terminal command to run BEAST2 from a \code{beast2_options}
}
\examples{
  if (is_beast2_installed()) {
    create_beast2_run_cmd_from_options(
      beast2_options = create_beast2_options()
    )
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_default_beast2_folder.R
\name{get_default_beast2_folder}
\alias{get_default_beast2_folder}
\title{Get the path to the folder where this package installs
BEAST2 by default}
\usage{
get_default_beast2_folder()
}
\value{
the path to the folder where this package installs
  BEAST2 by default
}
\description{
Get the path to the folder where this package installs
BEAST2 by default
}
\examples{
  message(get_default_beast2_folder())
}
\seealso{
Use \link{get_default_beast2_jar_path} to get the path
  to the BEAST2 jar file, when installed by this package
  Use \link{install_beast2} with default arguments
  to install BEAST2 to this folder.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_default_beast2_download_url.R
\name{get_default_beast2_download_url_win}
\alias{get_default_beast2_download_url_win}
\title{Get the BEAST2 download URL for Windows}
\usage{
get_default_beast2_download_url_win(
  beast2_version = beastier::get_default_beast2_version()
)
}
\arguments{
\item{beast2_version}{the version of BEAST2. By
default, this is the version as returned by
\link{get_default_beast2_version}}
}
\value{
the URL where BEAST2 can be downloaded from
}
\description{
Get the BEAST2 download URL for Windows
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_empty_beastier_folder.R
\name{check_empty_beastier_folder}
\alias{check_empty_beastier_folder}
\title{Check there are no files in the default \link{beastier} folder}
\usage{
check_empty_beastier_folder(beastier_folder = get_beastier_folder())
}
\arguments{
\item{beastier_folder}{the path to
the \link{beastier} temporary files folder}
}
\value{
Nothing.
Will \link{stop} if there are files in the \link{beastier} folder
}
\description{
Check there are no files in the default \link{beastier} folder.
The goal is to make sure no temporary files are left undeleted.
Will \link{stop} if there are files in the \link{beastier} folder
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_default_beast2_version.R
\name{get_default_beast2_version}
\alias{get_default_beast2_version}
\title{Get the default BEAST2 version that is used by beastier}
\usage{
get_default_beast2_version()
}
\value{
the BEAST2 version
}
\description{
Get the default BEAST2 version that is used by beastier
}
\examples{
get_default_beast2_version()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_bin_path.R
\name{is_bin_path}
\alias{is_bin_path}
\title{Is the path a path to the BEAST2 binary file?
Does not check if the file at that path is present}
\usage{
is_bin_path(path)
}
\arguments{
\item{path}{a string to a path}
}
\value{
TRUE if the path is a path to a BEAST2 binary file
}
\description{
Is the path a path to the BEAST2 binary file?
Does not check if the file at that path is present
}
\examples{
if (is_beast2_installed()) {
  # TRUE
  is_bin_path("beast")
  is_bin_path("BEAST.exe")
  is_bin_path(get_default_beast2_bin_path())
  # FALSE
  is_bin_path("launcher.jar")
  is_bin_path(get_default_beast2_jar_path())
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beast2_options_filenames.R
\name{get_beast2_options_filenames}
\alias{get_beast2_options_filenames}
\title{Extract the filenames from a `beast2_options`}
\usage{
get_beast2_options_filenames(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
the filenames from a `beast2_options`
}
\description{
Extract the filenames from a `beast2_options`
}
\examples{
beast2_options <- create_beast2_options()
get_beast2_options_filenames(beast2_options)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_os.R
\name{check_os}
\alias{check_os}
\title{Checks if the operating system is supported}
\usage{
check_os(os)
}
\arguments{
\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
Nothing. Will \link{stop} if the OS is unsupported
}
\description{
Checks if the operating system is supported
}
\examples{
check_os("mac")
check_os("unix")
check_os("win")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beastier_report.R
\name{beastier_report}
\alias{beastier_report}
\title{Creates a \link{beastier} report}
\usage{
beastier_report(
  beast2_folder = get_default_beast2_folder(),
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable
is installed: the BEAST2 executable is in a subfolder.
Use \link{get_default_beast2_folder} to get the default BEAST2
folder.
Use \link{get_default_beast2_bin_path} to get the full path to
the default BEAST2 executable.}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
No return value, the information will be shown using \link{message}
}
\description{
Creates a \link{beastier} report, to be used when reporting bugs.
Uses \link{message}
}
\examples{
beastier_report()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/extract_tracelog_filename_from_beast2_input_file.R
\name{extract_tracelog_filename_from_beast2_input_file}
\alias{extract_tracelog_filename_from_beast2_input_file}
\title{Internal function to extract the tracelog filename for a BEAST2 input file}
\usage{
extract_tracelog_filename_from_beast2_input_file(input_filename)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}
}
\value{
the name of the tracelog file
}
\description{
Extract the tracelog filename for a BEAST2 input file
}
\examples{
beast2_input_filename <- get_beastier_tempfilename()
tracelog_filename <- get_beastier_tempfilename()
beautier::create_beast2_input_file_from_model(
  input_filename = beautier::get_beautier_path("test_output_0.fas"),
  output_filename = beast2_input_filename,
  inference_model = beautier::create_inference_model(
    mcmc = beautier::create_mcmc(
      tracelog = beautier::create_tracelog(
        filename = tracelog_filename
      )
    )
  )
)
extract_tracelog_filename_from_beast2_input_file(
  input_filename = beast2_input_filename
)
file.remove(beast2_input_filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_unique_ids.R
\name{has_unique_ids}
\alias{has_unique_ids}
\title{Determine if the XML text has unique parameter IDs}
\usage{
has_unique_ids(text)
}
\arguments{
\item{text}{the XML as text}
}
\value{
TRUE if all parameter IDs are unique, FALSE otherwise
}
\description{
Determine if the XML text has unique parameter IDs
}
\examples{
line_1 <- "<parameter id=\"RealParameter.1\" ...</parameter>"
line_2 <- "<parameter id=\"RealParameter.2\" ...</parameter>"
# Unique IDs
has_unique_ids(c(line_1, line_2))
# No unique ID
has_unique_ids(c(line_1, line_1))
}
\seealso{
to obtain the duplicate parameter IDs, use
  \code{\link{get_duplicate_param_ids}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_can_create_treelog_file.R
\name{check_can_create_treelog_file}
\alias{check_can_create_treelog_file}
\title{Internal function}
\usage{
check_can_create_treelog_file(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
Nothing.
Will \link{stop} if the MCMC's treelog file is absent and cannot be created.
}
\description{
Check if the MCMC's treelog file can be created.
Will \link{stop} if not
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beastier_tempfilename.R
\name{get_beastier_tempfilename}
\alias{get_beastier_tempfilename}
\title{Get a temporary filename}
\usage{
get_beastier_tempfilename(pattern = "file", fileext = "")
}
\arguments{
\item{pattern}{a non-empty character vector
giving the initial part of the name.}

\item{fileext}{a non-empty character vector
giving the file extension}
}
\value{
name for a temporary file
}
\description{
Get a temporary filename, similar to \link{tempfile},
except that it always writes to a temporary folder
named \link{beastier}.
}
\note{
this function is added to make sure no temporary
cache files are left undeleted
}
\examples{
get_beastier_tempfilename()
get_beastier_tempfilename(pattern = "my_pattern_")
get_beastier_tempfilename(fileext = ".ext")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_run_cmd.R
\name{create_beast2_run_cmd}
\alias{create_beast2_run_cmd}
\title{Creates the terminal command to run BEAST2}
\usage{
create_beast2_run_cmd(
  input_filename,
  output_state_filename,
  rng_seed = NA,
  n_threads = NA,
  use_beagle = FALSE,
  overwrite = FALSE,
  beast2_path = get_default_beast2_path(),
  verbose = FALSE
)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}

\item{output_state_filename}{name of the BEAST2 output file that
stores the state
(usually has a \code{.xml.state} extension)}

\item{rng_seed}{the random number generator seed of the BEAST2 run.
Must be a non-zero positive integer value or \link{NA}.
If \code{rng_seed} is \link{NA}, BEAST2 will pick a random seed}

\item{n_threads}{the number of computational threads to use.
Use \link{NA} to use the BEAST2 default of 1.}

\item{use_beagle}{use BEAGLE if present}

\item{overwrite}{if TRUE: overwrite the \code{.log}
and \code{.trees} files if one of these exists.
If FALSE, BEAST2 will not be started if
\itemize{
  \item{the \code{.log} file exists}
  \item{the \code{.trees} files exist}
  \item{the \code{.log} file created by BEAST2 exists}
  \item{the \code{.trees} files created by BEAST2 exist}
}}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
a character vector with the command and
  arguments to call BEAST2
}
\description{
Creates the terminal command to run BEAST2
}
\examples{
  if (is_beast2_installed()) {
    create_beast2_run_cmd(
      input_filename = "input.xml",
      output_state_filename = "output.xml.state",
      beast2_path = get_default_beast2_jar_path()
    )
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install_beast2.R
\name{install_beast2}
\alias{install_beast2}
\title{Deprecated function to install BEAST2}
\usage{
install_beast2(
  folder_name = rappdirs::user_data_dir(),
  beast2_version = beastier::get_default_beast2_version(),
  verbose = FALSE,
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{folder_name}{name of the folder where the BEAST2 files will
be put.
The name of the BEAST2 binary file will be at
\code{[folder_name]/beast/bin/beast}
The name of the BEAST2 jar file will be at
\code{[folder_name]/beast/lib/launcher.jar}}

\item{beast2_version}{the version of BEAST2. By
default, this is the version as returned by
\link{get_default_beast2_version}}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
Nothing. Gives a deprecation message using \link{stop}.
}
\description{
This function is deprecated as it violated CRAN policy.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beast2_example_filenames.R
\name{get_beast2_example_filenames}
\alias{get_beast2_example_filenames}
\title{Get a list with the full paths of all
BEAST2 example filenames}
\usage{
get_beast2_example_filenames(beast2_folder = get_default_beast2_folder())
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable
is installed: the BEAST2 executable is in a subfolder.
Use \link{get_default_beast2_folder} to get the default BEAST2
folder.
Use \link{get_default_beast2_bin_path} to get the full path to
the default BEAST2 executable.}
}
\value{
a list with the full paths of all
  BEAST2 example filenames
}
\description{
Get a list with the full paths of all
BEAST2 example filenames
}
\examples{
if (is_beast2_installed()) {
  get_beast2_example_filenames()
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/check_beast2_options_do_not_overwrite_existing_files.R
\name{check_beast2_options_do_not_overwrite_existing_files}
\alias{check_beast2_options_do_not_overwrite_existing_files}
\title{Internal function}
\usage{
check_beast2_options_do_not_overwrite_existing_files(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
Nothing. Will \link{stop} if a file is threatened to be overwritten
}
\description{
Check if the \code{beast2_options} will not overwrite
existing files, when the 'overwrite' options is set to \code{FALSE}.
}
\details{
Will \link{stop} if a file is threatened to be overwritten
}
\examples{
check_beast2_options_do_not_overwrite_existing_files(
  beast2_options = create_beast2_options()
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uninstall_beast2.R
\name{uninstall_beast2}
\alias{uninstall_beast2}
\title{Deprecated function to uninstall BEAST2}
\usage{
uninstall_beast2(
  folder_name = rappdirs::user_data_dir(),
  os = rappdirs::app_dir()$os,
  verbose = FALSE
)
}
\arguments{
\item{folder_name}{name of the folder where the BEAST2 files are installed.
The name of the BEAST2 binary file will be at
\code{[folder_name]/beast/bin/beast}
The name of the BEAST2 jar file will be at
\code{[folder_name]/beast/lib/launcher.jar}}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
Nothing.
A deprecation message using \link{stop} will be triggered
}
\description{
Deprecated function to uninstall BEAST2
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_version_cmd.R
\name{create_beast2_version_cmd_jar}
\alias{create_beast2_version_cmd_jar}
\title{Creates the terminal command to version a BEAST2 input file
using a call to the \code{launcher.jar} file}
\usage{
create_beast2_version_cmd_jar(beast2_jar_path = get_default_beast2_jar_path())
}
\arguments{
\item{beast2_jar_path}{name of the BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
a character vector, of which the first element
  is the command (\code{java}, in this case),
  and the others are arguments (\code{-jar}, in this case, followed
  by more arguments.
}
\description{
Creates the terminal command to version a BEAST2 input file
using a call to the \code{launcher.jar} file
}
\examples{
if (is_beast2_installed()) {
  create_beast2_version_cmd_jar()
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_input_filename.R
\name{check_input_filename_validity}
\alias{check_input_filename_validity}
\title{Checks the input filename.
Will stop if there is a problem with the input filename.}
\usage{
check_input_filename_validity(
  beast2_options,
  input_filename = "deprecated",
  beast2_path = "deprecated",
  verbose = "deprecated"
)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}

\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
nothing. Will call \code{\link{stop}} if the input file is invalid
}
\description{
Checks the input filename.
Will stop if there is a problem with the input filename.
}
\examples{
if (is_beast2_installed()) {
  check_input_filename_validity(
   create_beast2_options(
     input_filename = get_beastier_path("2_4.xml")
   )
 )
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_jar_path.R
\name{is_jar_path}
\alias{is_jar_path}
\title{Is the path a path to the BEAST2 jar file?
Does not check if the file at that path is present}
\usage{
is_jar_path(path)
}
\arguments{
\item{path}{a string to a path}
}
\value{
TRUE if the path is a path to a BEAST2 jar file
}
\description{
Is the path a path to the BEAST2 jar file?
Does not check if the file at that path is present
}
\examples{
# Returns TRUE
is_jar_path("beast.jar")
is_jar_path("launcher.jar")
is_jar_path(get_default_beast2_jar_path())
# Returns FALSE
is_jar_path("beast")
is_jar_path(get_default_beast2_bin_path())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_random_fasta.R
\name{create_random_fasta}
\alias{create_random_fasta}
\title{Create a random FASTA file}
\usage{
create_random_fasta(
  n_taxa,
  sequence_length,
  fasta_filename,
  taxa_name_ext = ""
)
}
\arguments{
\item{n_taxa}{The number of taxa}

\item{sequence_length}{a DNA sequence length, in base pairs}

\item{fasta_filename}{a FASTA filename.}

\item{taxa_name_ext}{the extension of the taxa names}
}
\value{
Nothing, creates a FASTA file
}
\description{
Create a random FASTA file
}
\examples{
fasta_filename <- get_beastier_tempfilename()
create_random_fasta(
  n_taxa = 5,
  sequence_length = 20,
  fasta_filename = fasta_filename
)
file.remove(fasta_filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/extract_treelog_filename_from_beast2_input_file.R
\name{extract_treelog_filename_from_beast2_input_file}
\alias{extract_treelog_filename_from_beast2_input_file}
\title{Internal function to extract the treelog filename for a BEAST2 input file}
\usage{
extract_treelog_filename_from_beast2_input_file(input_filename)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}
}
\value{
the treelog filename for a BEAST2 input file
}
\description{
Extract the treelog filename from a BEAST2 input file
}
\examples{
beast2_input_filename <- get_beastier_tempfilename()
beautier::create_beast2_input_file_from_model(
  input_filename = beautier::get_beautier_path("test_output_0.fas"),
  output_filename = beast2_input_filename
)
extract_treelog_filename_from_beast2_input_file(
  input_filename = beast2_input_filename
)
file.remove(beast2_input_filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_alignment.R
\name{is_alignment}
\alias{is_alignment}
\title{Determines if the input is an alignment of type \link[ape]{DNAbin}}
\usage{
is_alignment(input)
}
\arguments{
\item{input}{The input to be tested}
}
\value{
TRUE or FALSE
}
\description{
Determines if the input is an alignment of type \link[ape]{DNAbin}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_temp_input_filename.R
\name{create_temp_input_filename}
\alias{create_temp_input_filename}
\title{Create a temporary filename for the BEAST2 XML filename}
\usage{
create_temp_input_filename()
}
\value{
a temporary filename, that starts with `beast2_`
and has extension `.xml`
}
\description{
Create a temporary filename for the BEAST2 XML filename
}
\examples{
create_temp_input_filename()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beast2_main_class_name.R
\name{get_beast2_main_class_name}
\alias{get_beast2_main_class_name}
\title{Get the BEAST2 main class name.}
\usage{
get_beast2_main_class_name()
}
\value{
the BEAST2 main class name
}
\description{
One way to fix the error
\code{no main manifest attribute}
is to specify the main class name.
}
\examples{
get_beast2_main_class_name()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_on_ci.R
\name{is_on_appveyor}
\alias{is_on_appveyor}
\title{Determines if the environment is AppVeyor}
\usage{
is_on_appveyor()
}
\value{
TRUE if run on AppVeyor, FALSE otherwise
}
\description{
Determines if the environment is AppVeyor
}
\examples{
  if (is_on_appveyor()) {
    message("Running on AppVeyor")
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_input_filename.R
\name{check_input_filename}
\alias{check_input_filename}
\title{Checks the input filename.
Will stop if there is a problem with the input filename.}
\usage{
check_input_filename(input_filename)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}
}
\value{
Nothing.
Will \link{stop} if the input file is invalid
}
\description{
Checks the input filename.
Will stop if there is a problem with the input filename.
}
\examples{
check_input_filename(
  get_beastier_path("beast2_example_output.log")
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_temp_state_filename.R
\name{create_temp_state_filename}
\alias{create_temp_state_filename}
\title{Create a temporary file for the BEAST2 XML output file that
stores its state.}
\usage{
create_temp_state_filename()
}
\value{
a temporary filename, that starts with `beast2_`
and has extension `.xml.state`
}
\description{
Create a temporary file for the BEAST2 XML output file that
stores its state.
}
\examples{
create_temp_state_filename()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_validate_cmd.R
\name{create_beast2_validate_cmd_bin}
\alias{create_beast2_validate_cmd_bin}
\title{Creates the terminal command to validate a BEAST2 input file
using a call to the \code{launcher.jar} file}
\usage{
create_beast2_validate_cmd_bin(
  input_filename,
  beast2_bin_path = get_default_beast2_bin_path()
)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}

\item{beast2_bin_path}{name of the BEAST2 binary file
(usually simply \code{beast}).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path}
}
\value{
a character vector, of which the first element
  is the command (\code{java}, in this case),
  and the others are arguments (\code{-jar}, in this case, followed
  by more arguments.
}
\description{
Creates the terminal command to validate a BEAST2 input file
using a call to the \code{launcher.jar} file
}
\examples{
  if (is_beast2_installed() && is_on_ci()) {
    create_beast2_validate_cmd_bin(
      input_filename = "input.xml"
    )
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_beast2_from_options.R
\name{run_beast2_from_options}
\alias{run_beast2_from_options}
\title{Run BEAST2}
\usage{
run_beast2_from_options(beast2_options = create_beast2_options())
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
The text sent to \code{STDOUT} and \code{STDERR}.
  It will create the file with name \code{output_state_filenames}
}
\description{
Run BEAST2
}
\examples{
if (is_beast2_installed() && is_on_ci()) {
  beast2_options <- create_beast2_options(
    input_filename = get_beastier_path("2_4.xml")
  )
  run_beast2_from_options(beast2_options)
  file.remove(beast2_options$output_state_filename)
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_continue_cmd_from_options.R
\name{create_beast2_continue_cmd_from_options}
\alias{create_beast2_continue_cmd_from_options}
\title{Creates the terminal command to run BEAST2 from a \code{beast2_options}}
\usage{
create_beast2_continue_cmd_from_options(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
a character vector with the command and
  arguments to call BEAST2
}
\description{
If the BEAST2 input \code{.xml} filename
or the BEAST2 state \code{.state.xml} filename
contain spaces, these filenames are quoted,
so that the command-line interface to BEAST2 correctly parses its arguments
}
\examples{
  if (is_beast2_installed()) {
    create_beast2_continue_cmd_from_options(
      beast2_options = create_beast2_options()
    )
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_default_beast2_jar_path.R
\name{get_default_beast2_jar_path}
\alias{get_default_beast2_jar_path}
\title{Get the default BEAST2 jar file's path}
\usage{
get_default_beast2_jar_path(
  beast2_folder = beastier::get_default_beast2_folder(),
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable
is installed: the BEAST2 executable is in a subfolder.
Use \link{get_default_beast2_folder} to get the default BEAST2
folder.
Use \link{get_default_beast2_bin_path} to get the full path to
the default BEAST2 executable.}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
the default BEAST2 jar file's path
}
\description{
Get the default BEAST2 jar file's path
}
\examples{
get_default_beast2_jar_path()
}
\seealso{
Use \link{get_default_beast2_folder} to get the default
  folder in which BEAST2 is installed.
  Use \link{install_beast2} with default arguments
  to install BEAST2 to this location.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beastier_tempfolder.R
\name{create_beastier_tempfolder}
\alias{create_beastier_tempfolder}
\title{Create the temporary folder as used by \link{beastier}}
\usage{
create_beastier_tempfolder()
}
\value{
nothing
}
\description{
Create the temporary folder as used by \link{beastier}
}
\examples{
create_beastier_tempfolder()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_paths.R
\name{get_beastier_paths}
\alias{get_beastier_paths}
\title{Get the full paths of files in the \code{inst/extdata} folder}
\usage{
get_beastier_paths(filenames)
}
\arguments{
\item{filenames}{the files' names, without the path}
}
\value{
the filenames' full paths. Will \code{stop} if a file
  is absent in the \code{inst/extdata} folder
}
\description{
Get the full paths of files in the \code{inst/extdata} folder
}
\examples{
get_beastier_paths(
  c(
    "beast2_example_output.log",
    "beast2_example_output.trees",
    "beast2_example_output.xml",
    "beast2_example_output.xml.state"
  )
)
}
\seealso{
for one file, use \code{\link{get_beastier_path}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/upgrade_beast2.R
\name{upgrade_beast2}
\alias{upgrade_beast2}
\title{Deprecated function to upgrade BEAST2.}
\usage{
upgrade_beast2(
  folder_name = rappdirs::user_data_dir(),
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{folder_name}{name of the folder where the BEAST2 files will
be put.
The name of the BEAST2 binary file will be at
\code{[folder_name]/beast/bin/beast}
The name of the BEAST2 jar file will be at
\code{[folder_name]/beast/lib/launcher.jar}}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
Nothing.
A deprecation message using \link{stop} will be triggered
}
\description{
Deprecated function to upgrade BEAST2.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_can_create_tracelog_file.R
\name{check_can_create_tracelog_file}
\alias{check_can_create_tracelog_file}
\title{Internal function to check if the MCMC's tracelog file can be created.}
\usage{
check_can_create_tracelog_file(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
Nothing.
Will \link{stop} if the MCMC's tracelog file is absent and cannot be created.
}
\description{
Check if the MCMC's tracelog file can be created.
Will \link{stop} if not.
If the tracelog file already exists,
it is assumed that a new file can be created,
by overwriting the existing one.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_rng_seed.R
\name{check_rng_seed}
\alias{check_rng_seed}
\title{Check if the input is a valid RNG seed.}
\usage{
check_rng_seed(rng_seed)
}
\arguments{
\item{rng_seed}{the random number generator seed of the BEAST2 run.
Must be a non-zero positive integer value or \link{NA}.
If \code{rng_seed} is \link{NA}, BEAST2 will pick a random seed}
}
\value{
Nothing.
Will \link{stop} if the RNG seed is invalid
}
\description{
Will \link{stop} if not.
}
\examples{
# Numbers from 1 and higher are valid RNG seeds
check_rng_seed(1)
check_rng_seed(2)
# Also NA is a valid RNG seed
check_rng_seed(NA)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_validate_cmd.R
\name{create_beast2_validate_cmd}
\alias{create_beast2_validate_cmd}
\title{Creates the terminal command to validate a BEAST2 input file}
\usage{
create_beast2_validate_cmd(
  input_filename,
  beast2_path = get_default_beast2_path()
)
}
\arguments{
\item{input_filename}{the name of a BEAST2 input XML file.
This file usually has an \code{.xml} extension.
Use \link{create_temp_input_filename} to create a temporary
filename with that extension.}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
a character vector, of which the first element
  is the command (\code{java}, in this case),
  and the others are arguments (\code{-jar}, in this case, followed
  by more arguments.
}
\description{
Creates the terminal command to validate a BEAST2 input file
}
\examples{
  if (is_beast2_installed() && is_on_ci()) {
    create_beast2_validate_cmd(
      input_filename = "input.xml"
    )
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_beast2_input_lines.R
\name{are_beast2_input_lines}
\alias{are_beast2_input_lines}
\title{Would these lines of text, when written to a file,
  result in a valid BEAST2 input file?}
\usage{
are_beast2_input_lines(
  lines,
  verbose = FALSE,
  method = ifelse(is_on_ci(), "deep", "fast"),
  beast2_path = get_default_beast2_path()
)
}
\arguments{
\item{lines}{lines of text}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}

\item{method}{the method to check. Can be 'deep' or 'fast'.
The 'deep' method uses BEAST2 to validate the complete file.
The 'fast' method uses some superficial tests (for example: if all
IDs are unique)}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
TRUE if the text is valid, FALSE if not
}
\description{
Would these lines of text, when written to a file,
  result in a valid BEAST2 input file?
}
\examples{
if (is_beast2_installed() && is_on_ci()) {
  are_beast2_input_lines(get_beastier_path("anthus_2_4.xml"))
}
}
\seealso{
Use \code{\link{is_beast2_input_file}} to check a file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_can_create_state_output_file.R
\name{check_can_create_state_output_file}
\alias{check_can_create_state_output_file}
\title{Internal function}
\usage{
check_can_create_state_output_file(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
Nothing.
Will \link{stop} if the state output file cannot be created.
}
\description{
Check if the state output file
can be created. Will \link{stop} otherwise
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_beast2_input_lines.R
\name{are_beast2_input_lines_deep}
\alias{are_beast2_input_lines_deep}
\title{Would these lines of text, when written to a file,
  result in a valid BEAST2 input file?}
\usage{
are_beast2_input_lines_deep(
  lines,
  verbose = FALSE,
  beast2_path = get_default_beast2_path()
)
}
\arguments{
\item{lines}{lines of text}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}

\item{beast2_path}{name of either a BEAST2 binary file
(usually simply \code{beast})
or a BEAST2 jar file
(usually has a \code{.jar} extension).
Use \link{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \link{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\value{
TRUE if the text is valid, FALSE if not
}
\description{
Would these lines of text, when written to a file,
  result in a valid BEAST2 input file?
}
\examples{
if (is_beast2_installed() && is_on_ci()) {
  beast2_filename <- get_beastier_path("anthus_2_4.xml")
  text <- readLines(beast2_filename)
  are_beast2_input_lines_deep(text)
}
}
\seealso{
Use \code{\link{is_beast2_input_file}} to check a file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beast2_options_to_table.R
\name{beast2_options_to_table}
\alias{beast2_options_to_table}
\title{Convert a \code{beast2_options} to a table}
\usage{
beast2_options_to_table(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
a \link[tibble]{tibble} with two columns, called `parameter`
and `value`. Each `parameter` is the name of the element of the
`beast2_options` structure, where the `value` on the same row holds
the value of that parameter
}
\description{
Convert a \code{beast2_options} to a table
}
\examples{
beast2_options_to_table(create_beast2_options())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_default_beast2_bin_path.R
\name{get_default_beast2_bin_path}
\alias{get_default_beast2_bin_path}
\title{Get the default BEAST2 binary file (\code{beast}, that is) path}
\usage{
get_default_beast2_bin_path(
  beast2_folder = get_default_beast2_folder(),
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable
is installed: the BEAST2 executable is in a subfolder.
Use \link{get_default_beast2_folder} to get the default BEAST2
folder.
Use \link{get_default_beast2_bin_path} to get the full path to
the default BEAST2 executable.}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
the default BEAST2 binary file's path
}
\description{
Get the default BEAST2 binary file (\code{beast}, that is) path
}
\examples{
if (is_beast2_installed()) {
  get_default_beast2_bin_path()
}
}
\seealso{
Use \link{get_default_beast2_folder} to get the default
  folder in which BEAST2 is installed.
  Use \link{install_beast2} with default arguments
  to install BEAST2 to this location.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_nexus_as_fasta.R
\name{save_nexus_as_fasta}
\alias{save_nexus_as_fasta}
\title{Save a NEXUS file as a FASTA file}
\usage{
save_nexus_as_fasta(nexus_filename, fasta_filename)
}
\arguments{
\item{nexus_filename}{name of an existing NEXUS file}

\item{fasta_filename}{name of the FASTA file to be created}
}
\value{
nothing. The NEXUS file will be saved as a FASTA file
}
\description{
Save a NEXUS file as a FASTA file
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_can_create_dir_for_state_output_file.R
\name{check_can_create_dir_for_state_output_file}
\alias{check_can_create_dir_for_state_output_file}
\title{Internal function}
\usage{
check_can_create_dir_for_state_output_file(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
Nothing.
Will \link{stop} if the folder for the state output file
cannot be created
}
\description{
Check if the folder for the state output file
can be created. Will \link{stop} otherwise
}
\examples{
check_can_create_dir_for_state_output_file(
  beast2_options = create_beast2_options()
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_default_beast2_download_url.R
\name{get_default_beast2_download_url_linux}
\alias{get_default_beast2_download_url_linux}
\title{Get the BEAST2 download URL for Linux}
\usage{
get_default_beast2_download_url_linux(
  beast2_version = beastier::get_default_beast2_version()
)
}
\arguments{
\item{beast2_version}{the version of BEAST2. By
default, this is the version as returned by
\link{get_default_beast2_version}}
}
\value{
the URL where BEAST2 can be downloaded from
}
\description{
Get the BEAST2 download URL for Linux
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_can_create_file.R
\name{check_can_create_file}
\alias{check_can_create_file}
\title{Internal function}
\usage{
check_can_create_file(filename, overwrite = TRUE)
}
\arguments{
\item{filename}{file that may or may not be created}

\item{overwrite}{if TRUE, if \code{filename} already exists, it
will be deleted by this function}
}
\value{
Nothing.
Will \link{stop} if a file cannot be created at a certain path.
}
\description{
Check that a file can be created at a certain path.
}
\details{
Will \link{stop} if not. Will \link{stop} if the file already exists.
Does so by creating an empty file at the path,
and then deleting it.
}
\examples{
check_can_create_file("my_local_file.txt")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/continue_beast2.R
\name{continue_beast2}
\alias{continue_beast2}
\title{Continue a BEAST2 run}
\usage{
continue_beast2(beast2_options = create_beast2_options())
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
The text sent to \code{STDOUT} and \code{STDERR}.
  It will create the file with name \code{output_state_filenames}
}
\description{
Continue a BEAST2 run
}
\examples{
if (is_beast2_installed() && is_on_ci()) {
  beast2_options <- create_beast2_options(
    input_filename = get_beastier_path("2_4.xml")
  )
  run_beast2_from_options(beast2_options)
  continue_beast2(beast2_options)
  file.remove(beast2_options$output_state_filename)
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beastier.R
\docType{package}
\name{beastier}
\alias{beastier}
\title{\code{beastier}: A package to call BEAST2.}
\description{
\code{beastier} allows to call BEAST2, a popular
Bayesian phylogenetics tool, using
an R interface. 'beastier' closely follows the interface
of BEAST2, including its default settings.
}
\examples{
beast2_options <- create_beast2_options(
  input_filename = get_beastier_path("2_4.xml")
)

if (is_beast2_installed() && is_on_ci()) {
  run_beast2_from_options(beast2_options)
  file.remove(beast2_options$output_state_filename)
}
}
\seealso{
These are packages associated with \code{beastier}:
\itemize{
  \item{
    The package \code{beautier} can create
    BEAST2 input files from R
  }
  \item{
    The package \code{tracerer} can parse
    BEAST2 output files from R
  }
  \item{
    The package \code{babette} combines the
    functionality of \code{beautier},
    \code{beastier} and \code{tracerer}
    into a single workflow
  }
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_identical_alignments.R
\name{are_identical_alignments}
\alias{are_identical_alignments}
\title{Determines if the two alignments are equal}
\usage{
are_identical_alignments(p, q)
}
\arguments{
\item{p}{the first alignment}

\item{q}{the second alignment}
}
\value{
TRUE or FALSE
}
\description{
Determines if the two alignments are equal
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_default_beast2_path.R
\name{get_default_beast2_path}
\alias{get_default_beast2_path}
\title{Get the default BEAST2 path}
\usage{
get_default_beast2_path(
  beast2_folder = beastier::get_default_beast2_folder(),
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable
is installed: the BEAST2 executable is in a subfolder.
Use \link{get_default_beast2_folder} to get the default BEAST2
folder.
Use \link{get_default_beast2_bin_path} to get the full path to
the default BEAST2 executable.}

\item{os}{name of the operating system,
must be \code{unix} (Linux, Mac) or \code{win} (Windows)}
}
\value{
the default BEAST2 path
}
\description{
Get the default BEAST2 path
}
\examples{
if (is_beast2_installed()) {
  get_default_beast2_path()
}
}
\seealso{
Use \link{get_default_beast2_bin_path}
    to get the default path to the BEAST2 binary file.
  Use \link{get_default_beast2_jar_path}
    to get the default path to the BEAST2 jar file.
  Use \link{get_default_beast2_folder} to get the default
  folder in which BEAST2 is installed.
  Use \link{install_beast2} with default arguments
  to install BEAST2 to this location.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_random_alignment.R
\name{create_random_alignment}
\alias{create_random_alignment}
\title{Create a random alignment}
\usage{
create_random_alignment(n_taxa, sequence_length, rate = 1, taxa_name_ext = "")
}
\arguments{
\item{n_taxa}{The number of taxa}

\item{sequence_length}{The number of base pairs the alignment will have}

\item{rate}{mutation rate}

\item{taxa_name_ext}{the extension of the taxa names}
}
\value{
an alignment of class \link[ape]{DNAbin}
}
\description{
Create a random alignment
}
\examples{
 alignment <- create_random_alignment(
   n_taxa = 5,
   sequence_length = 10
 )
 image(alignment)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beast2_options.R
\name{check_beast2_options_data_types}
\alias{check_beast2_options_data_types}
\title{Check if the \code{beast2_options}, which is a list,
has all elements of the right data types}
\usage{
check_beast2_options_data_types(beast2_options)
}
\arguments{
\item{beast2_options}{a set of BEAST2 options,
that are the R equivalent of the BEAST2 command-line options,
as can be created by \link{create_beast2_options}}
}
\value{
nothing
}
\description{
Calls \code{stop} if not.
}
\seealso{
Use \link{check_beast2_options} to check
  the entire \code{beast2_options} object
}
\author{
Richèl J.C. Bilderbeek
}
