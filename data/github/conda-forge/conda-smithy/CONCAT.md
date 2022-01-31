Overview
--------

`conda-smithy` is a tool for combining a conda recipe with configurations to build using freely hosted CI services into a single repository, also known as a feedstock.
`conda-smithy` is still a work-in-progress, but when complete, `conda-smithy` will:

+ Create a git repo with a conda recipe and the files to run conda builds via CI
  services.
+ Register the repo on github and push it.
+ Connect the repo to the CI services travis-ci.com, appveyor.com, circleci.com, dev.azure.com
  (For travis-ci.com, configure your org or user to enable the service for all repos)

[![tests](https://github.com/conda-forge/conda-smithy/workflows/tests/badge.svg)](https://github.com/conda-forge/conda-smithy/actions?query=workflow%3Atests)
[![Coverage Status](https://coveralls.io/repos/github/conda-forge/conda-smithy/badge.svg?branch=main)](https://coveralls.io/github/conda-forge/conda-smithy?branch=main)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Installation
------------

The easiest way to install conda-smithy is to use conda and conda-forge:

```
conda install -n root -c conda-forge conda-smithy
```

To install conda-smithy from source, see the requirements file in `requirements.txt`, clone this
repo, and `python -m pip install .`.

Setup
-----

You need a token from github, travis-ci.com, appveyor.com and circleci.com to try out
`conda-smithy`. The commands which need this will tell you where to get these tokens and where to
place them. If you need help getting tokens please ask on the
[conda-forge google group](https://groups.google.com/forum/?hl=en#!forum/conda-forge).

You should be able to test parts of `conda-smithy` with whatever tokens you have.
For example, you should be able to `conda smithy register-github` without the CI service tokens.
Re-rendering an existing feedstock is also possible without CI service tokens set.

Re-rendering an existing feedstock
----------------------------------

Periodically feedstocks need to be upgraded to include new features. To do
this we use `conda-smithy` to go through a process called re-rendering.
Make sure you have installed `conda-smithy` before proceeding.
Re-rendering an existing feedstock is possible without CI service tokens set.

1. `cd <feedstock directory>`
2. `conda smithy rerender [--commit]`
3. Commit and push all changes

Optionally one can commit the changes automatically with `conda-smithy` version `1.4.1+`.
To do this just use the `--commit`/`-c` option. By default this will open an editor to make a commit.
It will provide a default commit message and show the changes to be added. If you wish to do this
automatically, please just use `--commit auto`/`-c auto` and it will use the stock commit message.

Making a new feedstock
----------------------

1. **Make the feedstock repo:** `conda smithy init
<directory_of_conda_recipe>`.     For a recipe called `foo`, this creates a
directory called `foo-feedstock`, populates it with CI setup skeletons, adds the recipe under
`recipe` and initializes it as a git repo.

2. **Create a github repo:** `conda smithy register-github --organization conda-forge ./foo-feedstock`.
This requires a github token. You can try it out with a github user account
instead of an organization by replacing the organization argument with
`--user github_user_name`. If you are interested in adding teams for your feedstocks,
you can provide the `--add-teams` option to create them. This can be done when creating
the feedstock or after.

3. **Register the feedstock with CI services:**
`conda smithy register-ci --organization conda-forge --feedstock_directory ./foo-feedstock`.
This requires tokens for the CI services. You can give the name of a user instead
of organization with `--user github_user_name`. By default this command requires an Anaconda/Binstar token
to be available in `~/.conda-smithy/anaconda.token`, or as BINSTAR_TOKEN in the environment. This can be opted
out of by specifying `--without-anaconda-token`, as such execpted package uploads will not be attempted.
     * For Azure, you will have to create a service connection with the same name as your github user or org
        `https://dev.azure.com/YOUR_ORG/feedstock-builds/_settings/adminservices`
     * For Azure builds, you will have to export the environment variable `AZURE_ORG_OR_USER` to point to your Azure org
     * If this is your first build on Azure, make sure to add [Library Variable Group](https://docs.microsoft.com/en-us/azure/devops/pipelines/process/variables?view=azure-devops&tabs=yaml%2Cbatch#share-variables-across-pipelines) containing your BINSTAR_TOKEN for automated anaconda uploads.

4. **Specify the feedstock channel and label:**
   Optionally, you can specify source channels and choose a channel to upload to in `recipe/conda_build_config.yaml`.
     ```yaml
     channel_sources:
       - mysourcechannel1,mysourcechannel2,conda-forge,defaults
     channel_targets:
       - target_channel target_label
     ```
   Default source channels are `conda-forge,defaults`. Default for channel targets is `conda-forge main`.

5. **Specify your branding in the README.md:**
   Optionally, you can specify the branding on the README.md file by adding the following the `conda-forge.yml` file:
   ```
   github:
     user_or_org: YOUR_GITHUB_USER_OR_ORG
   channels:
     targets:
     -
       - YOUR_ANACONDA_CHANNEL
   ```

6. **Re-render the feedstock:** ``conda smithy rerender --feedstock_directory ./foo-feedstock``

7. **Commit the changes:** ``cd foo-feedstock && git commit``, then push ``git push upstream master``.

Running a build
---------------

When everything is configured you can trigger a build with a push to the feedstock repo on github.

Releasing conda-smithy
----------------------

Before making a release, consult `@conda-forge/core` and wait some time for objections.

To release a new version of conda-smithy, you can use the
[rever](https://regro.github.io/rever-docs/index.html) release managment tool.
Run `rever` in the root repo directory with the version number you want to release.
For example,

```sh
$ rever 0.1.2
```


Conda-smithy in a nutshell
--------------------------

#### xkcd 1629: Tools

[![xkcd 1629: Tools](https://imgs.xkcd.com/comics/tools.png)](https://xkcd.com/1629/)

**Titletext**: *I make tools for managing job-hunting sites for people who make*
*tools for managing job-hunting sites for people who make tools for ...*
About click
===========

Home: http://click.pocoo.org/

Package license: BSD-3-Clause

Feedstock license: [BSD-3-Clause](https://github.com/conda-forge/click-test-feedstock/blob/master/LICENSE.txt)

Summary: A simple wrapper around optparse for powerful command line utilities.

Development: https://github.com/pallets/click

Documentation: http://click.pocoo.org/

Click is a Python package for creating beautiful command line interfaces
in a composable way with as little code as necessary.


Current build status
====================

Linux: [![Circle CI](https://circleci.com/gh/conda-forge/click-test-feedstock.svg?style=shield)](https://circleci.com/gh/conda-forge/click-test-feedstock)
OSX: [![TravisCI](https://travis-ci.com/conda-forge/click-test-feedstock.svg?branch=master)](https://travis-ci.com/conda-forge/click-test-feedstock)
Windows: [![AppVeyor](https://ci.appveyor.com/api/projects/status/github/conda-forge/click-test-feedstock?svg=True)](https://ci.appveyor.com/project/conda-forge/click-test-feedstock/branch/master)

Current release info
====================
Version: [![Anaconda-Server Badge](https://anaconda.org/conda-forge/click/badges/version.svg)](https://anaconda.org/conda-forge/click)
Downloads: [![Anaconda-Server Badge](https://anaconda.org/conda-forge/click/badges/downloads.svg)](https://anaconda.org/conda-forge/click)

Installing click
================

Installing `click` from the `conda-forge` channel can be achieved by adding `conda-forge` to your channels with:

```
conda config --add channels conda-forge
```

Once the `conda-forge` channel has been enabled, `click` can be installed with:

```
conda install click
```

It is possible to list all of the versions of `click` available on your platform with:

```
conda search click --channel conda-forge
```


About conda-forge
=================

conda-forge is a community-led conda channel of installable packages.
In order to provide high-quality builds, the process has been automated into the
conda-forge GitHub organization. The conda-forge organization contains one repository
for each of the installable packages. Such a repository is known as a *feedstock*.

A feedstock is made up of a conda recipe (the instructions on what and how to build
the package) and the necessary configurations for automatic building using freely
available continuous integration services. Thanks to the awesome service provided by
[CircleCI](https://circleci.com/), [AppVeyor](https://www.appveyor.com/)
and [TravisCI](https://travis-ci.com/) it is possible to build and upload installable
packages to the [conda-forge](https://anaconda.org/conda-forge)
[Anaconda-Cloud](https://anaconda.org/) channel for Linux, Windows and OSX respectively.

To manage the continuous integration and simplify feedstock maintenance
[conda-smithy](https://github.com/conda-forge/conda-smithy) has been developed.
Using the ``conda-forge.yml`` within this repository, it is possible to re-render all of
this feedstock's supporting files (e.g. the CI configuration files) with ``conda smithy rerender``.

For more information please check the [conda-forge documentation](https://conda-forge.org/docs/).

Terminology
===========

**feedstock** - the conda recipe (raw material), supporting scripts and CI configuration.

**conda-smithy** - the tool which helps orchestrate the feedstock.
                   Its primary use is in the construction of the CI ``.yml`` files
                   and simplify the management of *many* feedstocks.

**conda-forge** - the place where the feedstock and smithy live and work to
                  produce the finished article (built conda distributions)


Updating click-feedstock
========================

If you would like to improve the click recipe or build a new
package version, please fork this repository and submit a PR. Upon submission,
your changes will be run on the appropriate platforms to give the reviewer an
opportunity to confirm that the changes result in a successful build. Once
merged, the recipe will be re-built and uploaded automatically to the
`conda-forge` channel, whereupon the built conda packages will be available for
everybody to install and use from the `conda-forge` channel.
Note that all branches in the conda-forge/click-feedstock are
immediately built and any created packages are uploaded, so PRs should be based
on branches in forks and branches in the main repository should only be used to
build distinct package versions.

In order to produce a uniquely identifiable distribution:
 * If the version of a package **is not** being increased, please add or increase
   the [``build/number``](https://conda.io/docs/user-guide/tasks/build-packages/define-metadata.html#build-number-and-string).
 * If the version of a package **is** being increased, please remember to return
   the [``build/number``](https://conda.io/docs/user-guide/tasks/build-packages/define-metadata.html#build-number-and-string)
   back to 0.

   the [``build/number``](https://conda.pydata.org/docs/building/meta-yaml.html#build-number-and-string)
   back to 0.

<!--
Thank you for pull request.
Below are a few things we ask you kindly to self-check before getting a review. Remove checks that are not relevant.
-->
Checklist
* [ ] Added a ``news`` entry

<!--
Please note any issues this fixes using [closing keywords]( https://help.github.com/articles/closing-issues-using-keywords/ ):
-->

<!--
Please add any other relevant info below:
-->
=======================
conda-smithy Change Log
=======================

.. current developments

v3.16.2
====================

**Changed:**

* Happy New Year! The license now includes 2022.
* Default provider for ppc64le was changed to azure with emulation using qemu.

**Authors:**

* Isuru Fernando
* Bastian Zimmermann



v3.16.1
====================

**Fixed:**

* Fixed error in linter for ``matplotlib-base`` for multioutput recipes where the requirements are a list.

**Authors:**

* Matthew R. Becker



v3.16.0
====================

**Added:**

* Added rerendering token input to webservices github action and automerge github action.

**Authors:**

* Matthew R. Becker



v3.15.1
====================

**Added:**

* Added a hint for recipes in conda-forge to depend on matplotlib-base as opposed to
  matplotlib.

**Changed:**

* use python 3.9 on github actions and use mambaforge
* When building with boa, use mamba to install conda-build, etc.  This assumes that
  we are using a Mambaforge based docker image / runtime environment.
* For azure pipelines, the default windows image is changed to windows-2019

**Authors:**

* Isuru Fernando
* Matthew R. Becker
* Marius van Niekerk



v3.15.0
====================

**Added:**

* Conda smithy will now detect if a recipe uses ``compiler('cuda')``
and set the ``CF_CUDA_ENABLED`` environment variable to ``True`` if
so. This can for example be useful to distinguish different options
for builds with or without GPUs in ``conda_build_config.yaml``.
* Introduce utility function to facilitate the use case of running conda smithy
  commands from any sub-directory in the git repo checkout of a feedstock.

**Fixed:**

* Fixed typo in GitHub Actions template, where ``DOCKERIMAGE`` was wrongly specified in the matrix configuration. The CI step and its corresponding script expect ``DOCKER_IMAGE``.

**Authors:**

* Isuru Fernando
* Jaime Rodríguez-Guerra
* H. Vetinari
* Nehal J Wani



v3.14.3
====================

**Changed:**

* linux-aarch64 builds default is changed from native (drone) to emulated (azure).

**Authors:**

* Isuru Fernando
* Mike Taves



v3.14.2
====================

**Authors:**

* Isuru Fernando



v3.14.2
====================

**Added:**

* Download SDK to local folder when build-locally.py instead of to the system dir
* Added support for woodpecker CI support

**Authors:**

* Isuru Fernando



v3.14.1
====================

**Fixed:**

* Call ``docker pull`` then ``docker run`` (sometimes ``--pull`` is unavailable)

**Authors:**

* Matthew R. Becker
* John Kirkham



v3.14.0
====================

**Added:**

* ``test`` option in ``conda-forge.yml`` can now be used to configure testing.
  By default testing is done for all platforms. ``native_and_emulated`` value
  will do testing only if native or if there is an emulator. ``native`` value
  will do testing only if native.

**Deprecated:**

* ``test_on_native_only`` is deprecated. This is mapped to
  ``test: native_and_emulated``.

**Fixed:**

* Always pull a new version of the image used in a build
* Add workaround for Travis CI network issues (courtesy of @pkgw)

**Authors:**

* Isuru Fernando
* Marcel Bargull
* Matthew W. Thompson



v3.13.0
====================

**Added:**

* Added the ability to store conda build artifacts using the Github Actions provider. To enable, set `github_actions: {store_build_artifacts: true}` in conda-forge.yml.
* It is possible to set the lifetime of the Github Actions artifacts by setting the the `github_actions: {artifact_retention_days: 14}` setting in conda-forge.yml to the desired value. The default is 14 days.
* Support for ppc64le on drone CI has been added
* Added support for registering at a custom drone server by adding --drone-endpoint cli argument
* Added explicit check to not upload packages on PR builds.
* Added key ``github:tooling_branch_name`` to ``conda-forge.yml`` to enable
  setting the default branch for tooling repos.
* The linter will now warn if allowed ``pyXY`` selectors are used (e.g. ``py27``, ``py34``, ``py35``, ``py36``). For other versions (e.g. Python 3.8 would be ``py38``), these selectors are *silently ignored*  by ``conda-build``, so the linter will throw an error to prevent situations that might be tricky to debug. We recommend using ``py`` and integer comparison instead. Note that ``py2k`` and ``py3k`` are still allowed.
* Added support for self-hosted github actions runners

  In conda-forge.yml, add ``github_actions: self_hosted: true`` to
  enable self-hosted github actions runner. Note that self-hosted
  runners are currently configured to run only on push events
  and pull requests will not be built.

* Allow multiple providers per platform

  In conda-forge.yml, add ``provider: <platform>: ['ci_1', 'ci_2']``
  to configure multiple providers per platform.

**Changed:**

* Uploads are now allowed when building with ``mambabuild``!
* Azure build artifacts are now zipped before being uploaded, with some cache directories and the conda build/host/test environments removed, to make user download smaller and faster.
* A separate Azure build artifact, including only the conda build/host/test environments, is additionally created for failed builds.
* Azure artifact names are now only shortened (uniquely) when necessary to keep the name below 80 characters.
* Updated CircleCI xcode version to 13.0.0 to prevent failures.
* The conda-smithy git repo now uses ``main`` as the default branch.
* conda mambabuild is now the default build mode.  To opt out of this change set ``build_with_mambabuild`` to false in your ``conda-forge.yml``.
* Bump Windows ``base`` environment Python version to 3.9
* Support using ``build-locally.py`` natively on ``osx-arm64``.

**Fixed:**

* Azure artifact names are now unique when a job needs to be restarted (#1430).
* Azure artifact uploads for failed builds that failed because of broken symbolic links have now been fixed.
* Test suite now runs correctly on pyyaml 6
* Remove the miniforge installation before building with ``./build-locally.py`` on MacOS so that
  ``./build-locally.py`` can be run more than once without an error regarding an exisiting miniforge installation.

**Authors:**

* Isuru Fernando
* Matthew R. Becker
* Jaime Rodríguez-Guerra
* Uwe L. Korn
* Ryan Volz
* John Kirkham
* Wolf Vollprecht
* Marius van Niekerk
* Matthias Diener



v3.12
====================

**Authors:**

* Marius van Niekerk



v3.12
====================

**Changed:**

* conda smithy init will now copy over the conda-forge.yml from the source recipe directory (if present)

**Authors:**

* Marius van Niekerk



v3.11.0
====================

**Added:**

* The maximum number of parallel jobs a feedstock can run at once will be limited
  to ``50``. This will ensure that all projects have a fair access to CI resources
  without job-hungry feedstocks hogging the build queue.

**Fixed:**

* Add --suppress-variables flag to conda-build command in Windows template

**Authors:**

* Jaime Rodríguez-Guerra
* Billy K. Poon



v3.10.3
====================

**Fixed:**

* Linting of recipes with multiple URLs was broken in last release and is fixed now

**Authors:**

* Isuru Fernando



v3.10.2
====================

**Added:**

* Add a "--feedstock_config" option to the regenerate/rerender, update-anaconda-token, azure-buildid subcommands for providing an alternative path to the feedstock configuration file (normally "conda-forge.yml"). This allows different names or to put the configuration outside the feedstock root.
* Linter will now check for duplicates of conda packages using pypi name
* Validate the value of ``noarch``. (Should be ``python`` or ``generic``.)

**Changed:**

* Use ``ubuntu-latest`` instead of ``ubuntu-16`` in the Azure pipeline template.

**Fixed:**

* `short_config_name` is used at azure pipelines artifact publishing step.
* Duplicate feedstocks with only '-' vs '_' difference is now correctly checked.
* correctly detect use of `test/script` in outputs

**Authors:**

* Isuru Fernando
* Uwe L. Korn
* Ryan Volz
* Duncan Macleod
* fhoehle
* Ben Mares



v3.10.1
====================

**Added:**

* Allow osx builds in build-locally.py

**Changed:**

* Focal is now used for Linux builds on Travis CI

**Authors:**

* Isuru Fernando
* Matthew R. Becker
* Chris Burr





v3.10.0
====================

**Added:**

* Added `clone_depth` parameter for use in conda-forge.yml that sets the feedstock git clone depth for all providers (except CircleCI). By default (`clone_depth: none`), current behavior is maintained by using the provider's default checkout/clone settings. A full clone with no depth limit can be specified by setting `clone_depth: 0`.
* Log groups support for GitHub Actions
* Added support for Github Actions as a CI provider. Provider name to use in conda-forge.yml
  is `github_actions`. Note that Github Actions cannot be enabled as a CI provider for conda-forge
  github organization to prevent a denial of service for other infrastructure.
* Add instructions to feedstock README template for configuring strict channel priority.

**Changed:**

* The `ci-skeleton` command now creates a default conda-forge.yml that sets `clone_depth: 0` for full depth clones on all providers. This default supports expected behavior when using `GIT_DESCRIBE_*` to set version and build numbers in the recipe by ensuring that tags are present. This effectively changes the default clone behavior for the Github Action and Travis providers, as all other providers do a full clone by default.

**Fixed:**

* Prevent duplicated log group tags when ``set -x`` is enabled.
* Fix run_osx_build not failing early on setup error.
* Fix too long filenames for build done canary files.

**Authors:**

* Isuru Fernando
* Jaime Rodríguez-Guerra
* Ryan Volz
* Marcel Bargull
* Philippe Blain
* Matthew R. Becker
* Marcel Bargull



v3.9.0
====================

**Added:**

* Enabled multiple entries for ``key_add`` operations.
* Define Bash functions ``startgroup()`` and ``endgroup()`` that provide a
  provider-agnostic way to group or fold log lines for quicker visual inspection.
  In principle, this only affects Linux and MacOS, since Windows pipelines
  use CI native steps. So far, only Azure and Travis support this. In the other
  providers a fallback ``echo "<group name>"`` statement is supplied.
* Support `os_version` in `conda-forge.yml`
* Add use_local option to use the migrator from the feedstock

**Changed:**

* To cross compile for  ``win-32`` from ``win-64``, using ``target_platform``
  is no longer supported. Use ``build_platform: win_32: win64`` in ``conda-forge.yml``.
* `run_osx_build.sh` had hardcoded handlers for Travis log folding. These have
  been replaced with the now equivalent Bash functions.
* A lower bound on python version for noarch python is now required

**Fixed:**

* Fix "File name too long" error for many zip keys
  Replace config filenames by their short versions if filesystem limits
  are approached.
* Fix running ``./build-locally.py --debug`` with cross-compilation
* Fixed dead conda-docs link to the ``build/number`` explanation in the README template.
* Fixed rendering error where the recipe's ``conda_build_config.yaml`` is
  applied again, removing some variants.
* Fixed list formatting in the README.
* migration_ts and migrator_ts were both used in conda-smithy and migration_ts was removed in favour of migrator_ts

**Authors:**

* Isuru Fernando
* Matthew R. Becker
* Jaime Rodríguez-Guerra
* Chris Burr
* Leo Fang
* Marcel Bargull
* Wolf Vollprecht
* Hugo Slepicka
* Bastian Zimmermann



v3.8.6
====================

**Changed:**

* Run docker builds using ``delegated`` volume mounts.

**Fixed:**

* All keys zipped with ``docker_image`` are now handled properly.
* Changed CI configuration to not run tests on ``push`` events to branches that
  are not ``master``.
* CI runs on PRs from forks now.
* ``#`` is not a valid comment symbol on Windows and using it as part of a pipeline Batch step will cause a (harmless) error in the logs. It has been replaced by ``::`` instead.

**Security:**

* Use latest ``conda-incubator/setup-miniconda`` version to circumvent the GH Actions deprecations on Nov 16th

**Authors:**

* Isuru Fernando
* Matthew R Becker
* Matthew R. Becker
* Uwe L. Korn
* John Kirkham
* Jaime Rodríguez-Guerra



v3.8.5
====================

**Changed:**

* Moved CI to GitHub actions and removed travis-ci
* Use the shorter build ID instead of job ID to name Azure artifacts when they are stored. This helps prevent the artifact name from being too long, which would result in being unable to download it.
* Replaced travis-ci status badge w/ GitHub actions one.

**Fixed:**

* Faulty ``migrator_ts`` type check prevented manual migrations from happening (those that are not yet merged to ``conda-forge-pinning``).
* Previous release accidentally included a commit that made noarch: python
  recipes without a lower bound error. This was changed to a hint

**Authors:**

* Isuru Fernando
* Matthew R. Becker
* Ryan Volz
* Marius van Niekerk
* Jaime Rodríguez-Guerra



v3.8.4
====================

**Fixed:**

* conda-build 3.20.5 compatibility for ``target_platform`` being always defined.

**Authors:**

* Isuru Fernando



v3.8.3
====================

**Added:**

* conda-build 3.20.5 compatiblity
* New ``choco`` top-level key in ``conda-forge.yml`` enables windows builds
  to use chocolatey to install needed system packages. Currently, only Azure
  pipelines is supported.

**Authors:**

* Isuru Fernando
* Anthony Scopatz



v3.8.2
====================

**Changed:**

* Reverted bugfix for each compiler getting a CI job.

**Authors:**

* Matthew R. Becker



v3.8.1
====================

**Changed:**

* Removed the default concurrency limits for azure

**Fixed:**

* Fixed rendering to make sure CI jobs are generated for each compiler version.

**Authors:**

* Matthew R Becker
* Filipe Fernandes
* Matthew R. Becker
* Marius van Niekerk



v3.8.0
====================

**Added:**

* Generate Documentation and Development links into the README.md based on doc_url and dev_url
* Add hyperlink to feedstock license file
* Generate license_url as hyperlink in the README.md when it has been defined in the meta.yaml
* Add ``--without-anaconda-token`` option to register-ci command, keep default behaviour of requiring the token
* ``remote_ci_setup`` field in conda-forge.yml, which defaults to ``conda-forge-ci-setup=3`` allowing the user to override

**Changed:**

* Variant algebra now supports two new operations for adding/remove a key

These new options allow for handling complex migrations cases needed for the python migrations.
* Add support to ``build-locall.py`` to call ``conda debug``.
* Added note about behaviour to README.md
* CI templates now expand ``remote_ci_setup`` string from config for the ci setup package

**Removed:**

* Remove unneeded set_defaults() for --without-$CI args, ``action="store_false"`` already defaults to True if not given

**Fixed:**

* Removed the warning for azure token when rerendering

**Authors:**

* Isuru Fernando
* Johnny Willemsen
* Uwe L. Korn
* Tom Pollard
* Marius van Niekerk



v3.7.10
====================

**Removed:**

* Remove unused ``forge_config["upload_script"]`` logic

**Fixed:**

* Error with linting check for deletion of ``recipes/example/meta.yaml`` in staged-recipes

**Authors:**

* Joshua L. Adelman
* Tom Pollard



v3.7.9
====================

**Added:**

* ``test_on_native_only`` is now supported on osx too.

**Deprecated:**

* Unparsed `"upload_packages": False` from default conda-forge.yml, as not parsed & no longer reflective of defaults

**Fixed:**

* re-enabled `upload_packages` per provider to conda-forge.yml, which when set to False overrides default upload logic

**Authors:**

* Isuru Fernando
* Tom Pollard
* Joshua L. Adelman



v3.7.8
====================

**Added:**

* ``MACOSX_SDK_VERSION`` is added as an always used key

**Authors:**

* Isuru Fernando



v3.7.7
====================

**Added:**

* Publish conda build artifacts on Azure as pipeline artifacts when azure.store_build_artifacts flag is True in conda-forge.yml. The default is False.
* Add an option ``test_on_native_only`` to not run tests when cross compiling

**Changed:**

* Handle NameError when anaconda_token isn't defined in ci_register.py, inline with rotate_anaconda_token()
* MacOS image in CI is bumped to macOS 10.15

**Fixed:**

* Re add travis_wait support via idle_timeout_minutes

**Authors:**

* Isuru Fernando
* Ryan Volz
* Tom Pollard



v3.7.6
====================

**Added:**

* Added partial support for cross compiling (Unixes can compile for other unixes only)

**Changed:**

* linux-64 configs were changed from prefix ``linux`` to ``linux-64``
* ``target_platform`` is now always defined for non-noarch  recipes
* Raise RuntimeError on empty travis repo_info requests, to guard against later KeyErrors
* Provide the name of the feedstock for which the update-anaconda-token command
  was performed.
* GitHub Teams are now added to feedstocks by their ``slug`` (i.e., the name
  used to ``@``-mention them on ``github.com``) as opposed to their names.

**Deprecated:**

* Setting ``provider: linux`` is deprecated in favor of ``provider: linux_64``

**Fixed:**

* Use `simplejson` to catch `JSONDecodeError` when available. Fix #1368.

**Security:**

* Members and teams are now properly removed from feedstocks and feedstock
  maintenance teams.

**Authors:**

* Isuru Fernando
* Matthew R Becker
* Matthew R. Becker
* Hadrien Mary
* Maksim Rakitin
* Tom Pollard



v3.7.4
====================

**Added:**

* Use the anaconda API to retrieve the latest version number of ``conda-smithy`` and ``conda-forge-pinning``.
* Pass ``CPU_COUNT`` from the host environment to the docker build.
  (Convenient when building locally.)
* Add a flag to `register-github` to create a private repository.
* Add a `private_upload` key in conda config file. If set to True Anaconda upload will use the `--private` flag.
* Removes ``/opt/ghc`` on Azure Linux images to free up space
* Additional secrets can be passed to the build by setting `secrets: ["BINSTAR_TOKEN", "ANOTHER_SECRET"]`
  in `conda-forge.yml`. These secrets are read from the CI configuration and
  then exposed as environment variables. To make them visible to build scripts,
  they need to be whitelisted in `build.script_env` of `meta.yaml`.
  This can, e.g., be used to collect coverage statistics during a build or test
  and upload them to sites such as coveralls.

**Changed:**

* Return type of ``feedstocks.clone_all()`` from ``None`` to list of repositories
* Link to list of SPDX licenses in lint message.

**Fixed:**

* Use ``AzureConfig`` in ``render_README`` instead of calling a raw requests. It allows rendering on a private Azure CI organization.
* CI skeleton properly sets the build number
* use SPDX identifier for feedstock license
* Allow an empty conda-forge.yml.
* The repo name for output validation is now extracted in the CI services to avoid
  issues with bad rerenders for clones to non-standard locations.

**Security:**

* Added --suppress-variables so that CI secrets cannot be leaked by conda-build into CI logs.

**Authors:**

* Matthew R Becker
* Christopher J. Wright
* Matthew R. Becker
* Hadrien Mary
* Julian Rüth
* Uwe L. Korn
* John Kirkham
* Duncan Macleod
* Axel Huebl
* Thomas Hopkins
* Stuart Berg



v3.7.3
====================

**Fixed:**

* Get feedstock name from meta when registering with CI services.
* CODEOWNERS file no longer treats GitHub team names as case-sensitive.

**Authors:**

* Matthew R Becker
* Uwe L. Korn



v3.7.2
====================

**Changed:**

* Changed the automerge configuration to use conda-forge/automerge-action.

**Authors:**

* Matthew R Becker



v3.7.1
====================

**Added:**

* Added ci skip statements during token registration to reduce loads.
* Added tar as a dependency
* Option to specify the generated feedstock name via ``extra.feedstock-name``.
* Support self-hosted Azure agents

**Changed:**

* Changed the docker mount to the recipe directory to have read-write permissions instead
  of read-only.
* conda-forge-pinning package is now downloaded on the fly

**Fixed:**

* Fix folding scripts file in GH PRs
* Error when linting recipes with ``license_file: `` (i.e. no file specified)
* PSF-2.0 is not a deprecated license
* Fixed whitespace additions

**Authors:**

* Isuru Fernando
* Matthew R Becker
* Matthew R. Becker
* Chris Burr
* Leo Fang
* Uwe L. Korn



v3.7.0
====================

**Added:**

Added a linter check for already existing feedstocks that are not exact match, but may have underscore instead of dash, and vice versa.
* Added code to rotate anaconda tokens.
* Added new `pip-install`-based hooks for using a local copy of the
  `conda-forge-ci-setup` package.

**Changed:**

* Refactored OSX CI scripts to be based off of a single global script on all CI platforms.
* Renamed the feedstock token output files to not munge "-feedstock" from
  the names.

* Bumped the default version of the `conda-forge-ci-setup` package to 3 to
  support the new output validation service.

**Fixed:**

* Fixed bug in feedstock token registration that deleted other secrets from azure.
* Fixed bugs in tests for feedstock tokens.

**Security:**

* Added code to call the feedstock output validation service. You must have
  `conda_forge_output_validation` set to true in the `conda-forge.yml` to use
  this feature.

**Authors:**

* Matthew R Becker
* Matthew R. Becker
* Natasha Pavlovikj



v3.6.17
====================

**Added:**

* Added a linter check for jinja2 variables to be of the form ``{{<one space><variable name><one space>}}``.

**Changed:**

* Change azure.force default to False in conda-forge.yml (#1252)
* Use a faster script for removing homebrew on osx.

**Removed:**

* Removed No azure token warning when rerendering
* Deleting strawberry perl was removed as conda-forge-ci-setup now filters the PATH
* Removed fast finish script for travis as we now set the setting on travis

**Fixed:**

* Re-rendering now cleans old contents in ``.azure-pipelines``
* Fixed the drone CI badge
* Made yaml loading in conda_smithy thread safe

**Authors:**

* Isuru Fernando
* Matthew R Becker
* Matthew R. Becker
* John Kirkham
* Tim Snyder
* Peter Williams



**Changed:**

* Allow people to pass extra arguments to ``docker run`` by setting
  ``$CONDA_FORGE_DOCKER_RUN_ARGS``.

**Authors:**

* Peter K. G. Williams



v3.6.16
====================

**Changed:**

* Windows conda environment is activated before conda calls
* Moved the appveyor image to Visual Studio 2017.

**Fixed:**

* Linter now properly allows ``LicenseRef`` and ``-License`` in the license section.

**Authors:**

* Isuru Fernando
* Matthew R Becker
* Matthew R. Becker



v3.6.15
====================

**Added:**

* Linter allows LicenseRef custom licenses.

**Removed:**

* Other is not a recognized license anymore.

* Deprecated SPDX license are not recognized anymore.

**Authors:**

* Isuru Fernando
* Matthew R Becker
* Filipe Fernandes
* Matthew R. Becker
* Tim Snyder
* Dave Hirschfeld
* Nils Wentzell



v3.6.14
====================

**Fixed:**

* Package MANIFEST did not include the ``license_exceptions.txt`` file properly.

**Authors:**

* Matthew R. Becker



v3.6.13
====================

**Added:**

* Added code to validate feedstock tokens
* Added code to register FEEDSTOCK_TOKENS per CFEP-13
* Linter will now recommend SPDX expression for license entry

**Fixed:**

* Rerender use forge_config["recipe_dir"] instead of hardcoding "recipe" (#1254 & #1257)
* Fixed bug where BINSTAR_TOKEN's were not properly patched if they already
  existed for TravisCI.

**Authors:**

* Isuru Fernando
* Matthew R Becker
* Tim Snyder



v3.6.12
====================

**Fixed:**

* Fix bug with conda 4.6.14 on Windows

**Authors:**

* Filipe Fernandes
* Dave Hirschfeld



v3.6.11
====================

**Added:**

* Added feature to upload the BINSTAR_TOKEN for travis-ci.com directly
  through the API

**Changed:**

* Updated the version of macOS image to 10.14 for Azure Pipelines.
* If conda-forge-pinning package has migrations installed, use those
  migration yaml files instead of the ones from the feedstock if the
  timestamp field match and remove if the migration yaml has a
  timestamp and there's no corresponding one in conda-forge-pinning
  which indicates that the migration is over.

**Deprecated:**

* Deprecated storing BINSTAR_TOKENs in the conda-forge.yml for travis

**Authors:**

* Isuru Fernando
* Matthew R Becker
* Maksim Rakitin



v3.6.10
====================

**Fixed:**

* Fixed variant comparisons when the variant has a space

**Authors:**

* Isuru Fernando



v3.6.9
====================

**Added:**

* Add automerge github actions when rerendering
* Added the configuration file for the webservices github action

**Fixed:**

* Fix crash of linter when requirements contains packages that start with python in name

**Authors:**

* Isuru Fernando
* Matthew R Becker
* Matthew R. Becker
* Tim Werner



v3.6.8
====================

**Changed:**

* Changed the config name to remove * and space characters

**Authors:**

* Isuru Fernando
* Min RK



v3.6.7
====================

**Added:**

Non-noarch recipes shouldn't use version constraints on python and r-base.
The linter only checked for python, this PR addes the check for r-base.
* Added an option to skip adding webhooks

**Fixed:**

* Azure builds for OSX and Windows only attempt to upload if builds succeeded
  and the BINSTAR_TOKEN is available.

**Authors:**

* Isuru Fernando
* Mark Harfouche
* Natasha Pavlovikj



v3.6.6
====================

**Added:**

* ``conda smithy rerender`` now adds an automerge action if ``conda-forge.yml`` has ``bot: {automerge: True}`` set.
  This action merges PRs that are opened by the ``regro-cf-autotick-bot``, are passing, and have the ``[bot-automerge]``
  slug in the title.

**Fixed:**

* Fixed problems rendering the ``README.md`` for some ``Jinja2`` variables (#1215)

**Authors:**

* Christopher J. Wright
* Matthew R Becker
* Matthew R. Becker



v3.6.5
====================

**Added:**

* Added ``.gitignore`` entries when running ``ci-skeleton``.

**Fixed:**

* Fixed Jinja syntax error in ``ci-skeleton``.

**Authors:**

* Anthony Scopatz



v3.6.4
====================

**Added:**

* New ``conda smithy ci-skeleton`` subcommand that generates ``conda-forge.yml``
  and ``recipe/meta.yaml`` files for using conda-forge / conda-smithy as
  the CI configuration outside of configuration. Calling ``rerender`` after
  ``ci-skeleton`` will generate the configuration files. This is a great way to
  either bootstrap CI for a repo or continue to keep CI up-to-date.
  The ``recipe/meta.yaml`` that is generated is just a stub, and will need to
  be filled out for CI to properly build and test.

**Fixed:**

* Fix an issue with empty host
* Fix python lint for recipes with outputs



v3.6.3
====================

**Added:**

* Added a lint for common mistakes in python requirements
* Use shellcheck to lint ``*.sh`` files and provide findings as hints. Can be
  enabled via conda-forge.yaml (shellcheck: enabled: True), default (no entry)
  is False.
* Support aarch64 on travis-ci.com
* Support ppc64le on travis-ci.com
* Check that the current working directory is a feedstock before re-rendering.

**Changed:**

* Update travis feedstock registration to no longer generate anything for
travis-ci.org.



v3.6.2
====================

**Changed:**

* Changed the pipeline names in drone to less than 50 characters
* .scripts folder is also hidden in PR diffs

**Fixed:**

* Fixed a bug in configuring appveyor.yml



v3.6.1
====================

**Fixed:**

* Drone changed their service to no longer send the same environment variables. Changed to use ``$DRONE_WORKSPACE``.



v3.6.0
====================

**Added:**

* Ignore Drone CI files in GitHub diffs
* Run ``black --check`` on CI to verify code is formatted correctly

**Changed:**

* Platform independent files like `run_docker_build.sh` are moved to `.scripts` folder
* Standardize and test support for multiple docker images.
* refactored ``conda_smithy.lint_recipe.NEEDED_FAMILIES`` to top level so external projects can access
* Rerun ``black`` on the codebase.

**Fixed:**

* fix crash when host section was present but empty
* fix build-locally.py in skip_render by not attempting to chmod +x it
* ship conf file for black so everyone uses the same settings



v3.5.0
====================

**Added:**

* conda-smithy will remove the ``.github/CODEOWNERS`` file in case the recipe
  maintainers list is empty

**Changed:**

* Default windows provider was changed to azure.



v3.4.8
====================

**Fixed:**

* Don't make assumptions in ``conda_smithy/variant_algebra.py`` about the metadata



v3.4.7
====================

**Added:**

* Added a method to sync user in drone

**Changed:**

* Check that a project is registered if registering fails on drone
* Check that a project has the secret if adding secret fails on drone



v3.4.6
====================

**Added:**

* conda-smithy can now register packages on drone.io.  We plan on using this to help out with the aarch64
  architecture builds.

**Changed:**

* drone.io is now the default platform for aarch64 builds
* migrations folder changed from <feedstock_root>/migrations to <feedstock_root>/.ci_support/migrations

**Fixed:**

* Fix render_README crash when azure api returns 404



v3.4.5
====================

**Fixed:**

* YAML ``dump()`` now used ``pathlib.Path`` object.



v3.4.4
====================

**Fixed:**

* Updated conda-smithy to work with ruamel.yaml v0.16+.



v3.4.3
====================

**Changed:**

* In linting pins allow more than one space

**Fixed:**

* Don't lint setting build number



v3.4.2
====================

**Added:**

* Generating feedstocks with support for the linux-armv7l platform.
* test of the downgrade functionality of the new pinning system
* Mark generated files as generated so that github collapses them by deafult in diffs.
* The linter will now recomend fixes for malformed pins,
  suggesting a single space is inserted. For instance, both ``python>=3`` and
  ``python >= 3`` will ought to be ``python >=3``.
* New key ``upload_on_branch`` added to conda-forge.yml the value of which is checked
  against the current git branch and upload will be skipped if they are not equal.
  This is optional and an empty key skips the test.
* Added `CONDA_SMITHY_LOGLEVEL` environment variable to change verbosity
  of rendering. This can be either `debug` or `info`.

**Changed:**

* Add skip_render option to conda-forge.yaml. One could specify one or more filenames telling conda-smithy to skip making change on them. Files that could skip rendering include .gitignore, .gitattributes, README.md and LICENCE.txt.
* Reduced verbosity of rendering

**Fixed:**

* recipe-lint compatibility with ruamel.yaml 0.16
* Mock PY_VER in recipe check
* Fixed badge rendering in readme template.
* yum_requirements will now work on Travis based linux builds.
* requirements: update to conda-build>=3.18.3
* fix non-public conda import, use conda.exports
* requirements: replace pycrypto with pycryptodome



v3.4.1
====================

**Added:**

* license_file is required for GPL, MIT, BSD, APACHE, PSF

**Changed:**

* ``build-locally.py`` now uses ``python3`` even if ``python`` is ``python2`` (Python 3.6+ was already required)

**Removed:**

* Github issue, PR and contributing files are removed as they are in https://github.com/conda-forge/.github
* Support for python 2 Removed

**Fixed:**

* Fix configuring appveyor on repos starting with an underscore
* Fixed an issue where conda system variants could be used after rendering migrations.
* Fixed issue where only the last maintainer is review requested
* Unlicense is allowed
* Support newer ``shyaml`` versions by checking whether ``shyaml -h`` succeeds.



v3.4.0
====================

**Fixed:**

* bumped conda version check in CLI to 5.0 (from 4.7)



v3.3.7
====================

**Added:**

* Added codeowners file

**Fixed:**

* Fixed checking in .pyc files



v3.3.6
====================

**Fixed:**

* Indentation error in ``github.py``



v3.3.5
====================

**Added:**

* Added native aarch64 support for builds using Drone.io. This can be enabled by
  either using `provider: {linux_aarch64: drone}` or `provider: {linux_aarch64:
  native}` in the conda-forge.yml.

  Currently, drone has to be enabled manually as there is no automatic CI
  registration for repos.
* export CI env variable with CI provider name
* New ``build-locally.py`` script that is added to the root feedstock directory when
  ``conda smithy rerender`` is run. This script runs conda build locally. Currently
  it only fully supports running docker builds.
* print when adding new team to maintiners of feedstock

**Removed:**

* `docker.image` in conda-forge.yml is removed
* Removed the need for shyaml in CI env.

**Fixed:**

* removed empty lines causing current build status table to render as code
* build setup script overriding is now supported on azure too



v3.3.4
====================



v3.3.3
====================

**Added:**

* Added native ppc64le support to for travis-ci.  This can be enabled by either using
  `provider: {linux_ppc64le: travis}` or `provider: {linux_ppc64le: native}` in the conda-forge.yml.
  These will be the new default behavior going forward for ppc64le builds.  If native builds are not needed the
  qemu based builds on azure will continue to function as before.
* Added `DOCKER_IMAGE` variable to `run_docker_build.sh`

**Changed:**

* Fallback to default image in `run_docker_build.sh` if `shyaml` is not installed.

**Fixed:**

* Fixed badges for noarch builds using azure



v3.3.2
====================



v3.3.1
====================

**Fixed:**

* Use `config.instance_base_url` instead of `config.azure_team_instance` when creating new feedstocks



v3.3.0
====================

**Added:**

* Added a utility to retrieve the azure buildid.  This is needed to make badges for non-conda forge users.
* Added badges for azure ci builds.

**Changed:**

* Bumped up the maximum build time on azure to 6 hours!
* Switched default provider for osx and linux to be azure.
* ``conda-smithy regenerate`` now supports ``--check`` to see if regeneration can be performed
* Bumped the license year to 2019.
* Only suggest noarch in linting staged-recipes pull requests, not feedstocks.
  Refer to issues #1021, #1030, #1031. Linter is not checking all prerequisites for noarch.



v3.2.14
====================

**Added:**

* hint to suggest using python noarch, when the build requirements include pip and no compiler is specified.

**Fixed:**

* qemu activation fixed so that we can use sudo.



v3.2.13
====================

**Added:**

* Allow enabling aarch64 and ppc64le using default provider

**Changed:**

* Appveyor will now use the conda python3.x executable to run the fast-finish script.
* Azure windows builds are no longer silent.
* Azure build definition updating now works.

**Fixed:**

* yum_requirements will now work on azure based linux builds.



v3.2.12
====================

**Fixed:**

* Removed ``v`` from release that prevented conda-smithy version check from
  working properly.



v3.2.11
====================

**Fixed:**

* Secrets weren't getting passed to Azure properly.



v3.2.10
====================

**Changed:**

* Ran ``black`` on the codebase
* Added a few more always included keys.  These are required by the aarch64 migration.
These in particular are: ``cdt_arch``, ``cdt_name``,  ``BUILD``.



v3.2.9
====================



v3.2.8
====================

**Fixed:**

* conda-clean --lock does nothing.  Remove it.



v3.2.7
====================

**Fixed:**

* Fixed azure conditions for osx and win64



v3.2.6
====================

**Fixed:**

* Bugfix for uploading packages.



v3.2.5
====================

**Fixed:**

* Fixed docker image name from ``gcc7`` to ``comp7``.



v3.2.4
====================

**Fixed:**

* Fixed issue where azure was deleting linux configs for noarch packages.



v3.2.3
====================

**Added:**

* Added `conda-build` version to git commit message produced by `conda smithy regenerate`
* Made idle timeouts on travisci and circleci configurable.  To set this add to your `conda-forge-config.yml`

    .. code-block:: yaml

    idle_timeout_minutes: 30
None

* Added preliminary multiarch builds for aarch64 and ppc64le using qemu on azure.  This will be enabled by
means of a migrator at a later point in time.
Command line options are now available for the command `conda smithy register-ci`
to disable registration on a per-ci level. `--without-azure`, `--without-circle`,
`--without-travis`, and `--without-appveyor` can now be used in conjunction with
`conda smithy register-ci`.

**Changed:**

conda-build is now specified along side `conda-forge-ci-setup` installs so that it gets updated to the latest version available during each build.
* Moved NumFOCUS badge to "About conda-forge" section in the feedstock README.
* Removed ``branch2.0`` for the finding the fast-finish script, and changed it
  back to ``master``.

**Fixed:**

* Linter no longer fails if meta.yaml uses `os.sep`
* Fixed azure linux rendering caused by bad jinja rendering
* Linting only fails noarch recipes with selectors for host and runtime dependencies.



v3.2.2
====================

**Added:**

* recipe-maintainers can now be a conda-forge github team


**Fixed:**

* Azure fixed incorrect build setup
* Use setup_conda_rc for azure on windows
* Fixed creating feedstocks with conda-build 3.17.x
* Fixed bug in appveyor where custom channels are not used
* Added conda-forge when installing conda-forge-ci-setup to prevent Circle from changing channel priority




v3.2.1
====================

**Added:**

* Added support for rendering feedstock recipes for Azure pipelines.
  Presently this is enabled globally for all feedstocks going forward by default.
  Azure builds are configured to not publish artifacts to anaconda.org
* PR template asking for news entries
  (aka, I heard you like news, so I put a news item about adding news items into
  your news item, so you can add news while you add news)
* Feedstock maintainers are now listed in the README file.


**Removed:**

* Python 2.7 support has been dropped.  Conda-smithy now requires python >= 3.5.


**Fixed:**

* Fixes issue with Circle job definition where "filters are incompatible with
  workflows" when Linux is skipped. This was causing Linux jobs to be created
  and then fail on feedstocks where Linux and Circle were not needed.




v3.2.0
====================

**Changed:**

* updated toolchain lint to error


**Fixed:**

* The ``extra-admin-users`` flag can be None which is the default case. So, we have to check that before to make a loop on the entries of ``extra-admin-users`` list.
* The ``update-cb3`` command now handles ``toolchain3`` in the same way that
  ``toolchain`` is handled.




v3.1.12
====================

**Fixed:**

* fixed lint by checking that recipe-maintainers is an instance of
  ``collections.abc.Sequence``




v3.1.11
====================

**Changed:**

* Upgrade links to HTTPS and update link targets where necessary (#866)


**Removed:**

* Drop `vendored` package/directory. A remnant that is no longer used.


**Fixed:**

None

* Linter: packages without a `name` aren't actually in bioconda. (#872)
* Linter: handle new versions of `ruamel.yaml` appropriately instead of complaining about `expected to be a dictionary, but got a CommentedMap`. (#871)
* Fix missing newline in last line of generated readmes and add unit test for it (#864)




v3.1.10
====================

**Changed:**

- Change conda-smithy rerender text in PR template so that it is not invoked. (#858)


**Fixed:**

- Fix OrderedDict order not being kept (#854)




v3.1.9
====================

**Added:**

* Add merge_build_host: True #[win] for R packages in update-cb3


**Changed:**

* Package the tests




v3.1.8
====================

**Fixed:**

* Linter issue with multiple outputs and unexpected subsection checks




v3.1.7
====================

**Added:**

* Allow appveyor.image in conda-forge.yml to set the `appveyor image <https://www.appveyor.com/docs/build-environment/#choosing-image-for-your-builds>`_. (#808)
* Temporary travis user for adding repos  #815
* More verbose output for ``update-cb3``  #818
* ``.zip`` file support for ``update-cb3``  #832


**Changed:**

* Move noarch pip error to hint  #807
* Move biocona duplicate from error to hint  #809


**Fixed:**

- Fix OrderedDict representation in dumped yaml files (#820).
- Fix travis-ci API permission error (#812)
* Linter: recognize when tests are specified in the `outputs` section. (#830)




v3.1.6
====================

**Fixed:**

- Fix sorting of values of packages in `zip_keys` (#800)
- Fix `pin_run_as_build` inclusion for packages with `-` in their names (#796)
- Fix merging of configs when there are variants in outputs (#786, #798)
- Add `conda smithy update-cb3` command to update a recipe from conda-build v2 to v3 (##781)




v3.1.2
====================

**Added:**

None

* Require ``conda-forge-pinnings`` to run
None

* Update conda-build in the docker build script


**Changed:**

None

* Included package badges in a table
All of the people who have made at least one contribution to conda-smithy.
Authors are sorted by number of commits.

* Isuru Fernando
* Matthew R. Becker
* Christopher J. Wright
* Anthony Scopatz
* Phil Elson
* Filipe Fernandes
* Dougal J. Sutherland
* Jaime Rodríguez-Guerra
* shadow_walker
* Michael Sarahan
* Min RK
* Hadrien Mary
* Uwe L. Korn
* Ryan Volz
* Johnny Willemsen
* Julian Rüth
* Chris Burr
* Mark Harfouche
* Marcel Bargull
* Leo Fang
* Eric Dill
* Natasha Pavlovikj
* Maksim Rakitin
* Jan Schulz
* Josh Reichardt
* Justin Calamari
* Jason Furmanek
* Patrick Sodré
* John Kirkham
* C.A.M. Gerlach
* refraction-ray
* Leopold Talirz
* Tom Pollard
* Gonzalo Pena-Castellanos
* Daniel Bast
* Tim Snyder
* Bastian Zimmermann
* H. Vetinari
* Carlo
* Matthew Craig
* Joshua L. Adelman
* Wolf Vollprecht
* Marius van Niekerk
* Peter Williams
* Johannes Köster
* Matt McCormick
* John Blischak
* Jan Janßen
* xoviat
* Elmar Pruesse
* Patrick Sodré
* Josh Barnes
* Tobias Megies
* Jonathan Helmus
* Duncan Macleod
* Florian Rathgeber
* Bruno Oliveira
* Tony Kelman
* Guilherme Quentel Melo
* santi
* buijennifer
* melsyt
* gouarin
* Nicholas Bollweg
* Alex Goodman
* Lion Krischer
* ap--
* Peter Killick
* Henry Schreiner
* hajapy
* Axel Huebl
* Thomas Hopkins
* Hugo Slepicka
* fhoehle
* Ben Mares
* Matthias Diener
* Matthew W. Thompson
* lorenz
* Tom Augspurger
* Ryan May
* Thomas Robitaille
* roryk
* Richard Hattersley
* Dominik Kutra
* Morten Enemark Lund
* danielballan
* Max Linke
* Nathan Goldbaum
* cshaley
* David Brochart
* Julien Schueller
* Jason Grout
* Tim Werner
* Dave Hirschfeld
* Nils Wentzell
* Stuart Berg
* Philippe Blain
* Billy K. Poon
* Mike Taves
* Nehal J Wani
**Added:**

* <news item>

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
