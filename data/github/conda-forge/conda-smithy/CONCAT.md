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
