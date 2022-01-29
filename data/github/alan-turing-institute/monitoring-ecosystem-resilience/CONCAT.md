![Build status](https://api.travis-ci.com/alan-turing-institute/monitoring-ecosystem-resilience.svg?branch=develop)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alan-turing-institute/monitoring-ecosystem-resilience/master?filepath=notebooks)

[![Documentation Status](https://readthedocs.org/projects/pyveg/badge/?version=latest)](https://pyveg.readthedocs.io/en/latest/?badge=latest)

# monitoring-ecosystem-resilience
Repository for mini-projects in the Data science for Sustainable development project.

Currently the focus of code in this repository is understanding vegetation patterns in semi-arid environments.

The code in this repository is intended to perform three inter-related tasks:
* Download and process satellite imagery from Google Earth Engine.
* Generate simulated vegetation patterns.
* Calculate graph metrics to quantify the interconnectedness of vegetation in real and simulated images.

### Python

The tasks above are all implemented in Python in the *pyveg* package. See the [README.md](pyveg/README.md) in the `pyveg` subdirectory for details on installation and usage.

### R

The pattern-generation and graph-modelling are implemented in R in the *rveg* package.  See the [README.md](rveg/README.md) in the `rveg` directory for further details.
# Contributing to monitoring-ecosystem-resilience (the repo!)

**Welcome to the repository!**
We're excited you're here and want to contribute.

We hope that these guidelines make it as easy as possible to get involved.
If you have any questions that aren't discussed below, please let us know by opening an [issue](https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/issues).

We welcome all contributions from documentation to testing to writing code.
Don't let trying to be perfect get in the way of being good - exciting ideas are more important than perfect pull requests.

## Table of contents

- [Where to start: issues](#where-to-start-issues)
- [Making a change with a pull request](#making-a-change-with-a-pull-request)
  - [1. Comment on an existing issue or open a new issue referencing your addition](#1-comment-on-an-existing-issue-or-open-a-new-issue-referencing-your-addition)
  - [2. Create a new branch (if you have *write* access to the repository) or fork the repository to your profile (if you don't currently have _write_ access)](#2-create-a-new-branch-or-fork-the-repository-to-your-profile)
  - [3. Make the changes you've discussed](#3-make-the-changes-youve-discussed)
  - [4. Submit a pull request](#4-submit-a-pull-request)
- [Style guide](#style-guide)

## Where to start: issues

* **Issues** are individual pieces of work that need to be completed to move the project forwards.
A general guideline: if you find yourself tempted to write a great big issue that
is difficult to describe as one unit of work, please consider splitting it into two or more issues.

Before you open a new issue, please check if any of our [open issues](https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/issues) covers your idea already.

The list of labels for current issues includes:

- [![help-wanted](https://img.shields.io/badge/-help%20wanted-159818.svg)][labels-helpwanted] _These issues contain a task that a member of the team has determined we need additional help with._

  If you feel that you can contribute to one of these issues, we especially encourage you to do so!

- [![question](https://img.shields.io/badge/-question-cc317c.svg)][labels-question] _These issues contain a question that you'd like to have answered._

  Opening an issue is a great way to start a conversation and get your answer.

- [![good-first-issue](https://img.shields.io/badge/-good%20first%20issue-1b3487.svg)][labels-firstissue] _These issues are particularly appropriate if it is your first contribution to the repository, or to GitHub overall._

- [![Enhancement](https://img.shields.io/badge/-enhancement-84b6eb.svg)][labels-enhancement] _These issues are suggesting new features that can be added to the project._

  If you want to ask for something new, please try to make sure that your request is distinct from any others that are already in the queue.
  If you find one that's similar but there are subtle differences please reference the other enhancement in your issue.

- [![Bug](https://img.shields.io/badge/-bug-d73a4a.svg)][labels-bug] _These issues are reporting a problem or a mistake in the project._

  The more details you can provide the better!
  If you know how to fix the bug, please open an issue first and then submit a pull request.

- [![project-management](https://img.shields.io/badge/-project%20management-bfd86c.svg)][labels-project-management] _We like to model best practice, so the package itself is managed through these issues.

## Making a change with a pull request

We appreciate all contributions to monitoring-ecosystem-resilience.
**THANK YOU** for helping us.

All project management, conversations and questions related to the project happens here in the [monitoring-ecosystem-resilience repository][monitoring-ecosystem-resilience-repo].

In brief, the structure for making a contribution is as follows:
1. Identify a specific change that needs to be made to the repository. Open a new issue (after checking one does not already exist!) and describe the change, include why you are making it.
2. Create a new branch corresponding to this issue. The new branch will house all the changes that you make to the repository in an isolated location. As discussed in more detail below, new branches should be created using the latest version of the `develop` branch.
3. Make commits to the new branch you have created.
4. Submit a pull request to add the modifications in your new branch back into `develop`.

When a significant milestone has been reached, and the `develop` branch is known to be in a stable configuration, the `master` branch will be updated via a pull request from `develop`. In general, commits should not be made to either the `master` or `develop` branches. Pull requests to `develop` are fine (and encoraged), while pull requests to `master` will happen in a coordinated way.

The following steps are a more detailed guide to help you contribute in a way that will be easy for everyone to review and accept with ease.

### 1. Comment on an [existing issue](https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/issues) or open a new issue referencing your addition

This allows other members of the team to confirm that you aren't overlapping with work that's currently underway and that everyone is on the same page with the goal of the work you're going to carry out.

[This blog](https://www.igvita.com/2011/12/19/dont-push-your-pull-requests/) is a nice explanation of why putting this work in up front is so useful to everyone involved.

### 2. Create a new [branch][github-branches] or [Fork][github-fork] the [monitoring-ecosystem-resilience repository][monitoring-ecosystem-resilience-repo] to your profile

#### 2a) Create a branch
If you are a collaborator on the repository with write access, then you can make a [new branch][github-branches].  We recommend that you start from the latest version of the `develop` branch, and create a new one from there. This is the branch we use for active deleopment of the repository, while stable (but not cutting edge) versions are in the `master` branch. The name of your new branch should ideally be in the format: `<feature|bugfix>/<issue-number>-<short-description>`. For example, if you were addressing Issue number 111 which was about incorrect JSON filenames, it could be something like:
```
git checkout develop
git pull
git checkout -b bugfix/111-fix-json-filenames
```
Now you can go to step #3, where you actually fix the problem! :)

In case you want to learn more about "branching out", [this blog](https://nvie.com/posts/a-successful-git-branching-model/) details the different Git branching models.


#### 2b. Fork the repository

If you don't have write access to the repository, you can fork it to your own profile.
This is now your own unique copy of the repo.
Changes here won't affect anyone else's work, so it's a safe space to explore edits to the code!

Make sure to [keep your fork up to date][github-syncfork] with the master repository, otherwise you can end up with lots of dreaded [merge conflicts][github-mergeconflicts].

### 3. Make the changes you've discussed

Try to keep the changes focused.
If you submit a large amount of work all in one go it will be much more work for whomever is reviewing your pull request.

While making your changes, commit often and write good, detailed commit messages.
[This blog](https://chris.beams.io/posts/git-commit/) explains how to write a good Git commit message and why it matters.
It is also perfectly fine to have a lot of commits - including ones that break code.
A good rule of thumb is to push up to GitHub when you _do_ have passing tests then the continuous integration (CI) has a good chance of passing everything.

Please do not re-write history!
That is, please do not use the [rebase](https://help.github.com/en/articles/about-git-rebase) command to edit previous commit messages, combine multiple commits into one, or delete or revert commits that are no longer necessary.

### 4. Submit a [pull request][github-pullrequest]

A "pull request" is a request to "pull" the changes you have made in your branch back into another branch of the repository. The source branch will be the new branch you created in order to address the issue you created/choose. The destination branch should generally be `develop`, where all main code development takes place. Avoid making pull requests into the `master` branch (pull requests into master should happen in a coordinated way using a stable configuration of `develop` as the source branch).

We encourage you to open a pull request as early in your contributing process as possible.
This allows everyone to see what is currently being worked on.
It also provides you, the contributor, feedback in real time from both the community and the continuous integration as you make commits (which will help prevent stuff from breaking).

When you are ready to submit a pull request, make sure the contents of the pull request body do the following:
- Describe the problem you're trying to fix in the pull request, reference any related issues and use keywords fixes/close to automatically close them, if pertinent.
- List changes proposed in the pull request.
- Describe what the reviewer should concentrate their feedback on.

If you have opened the pull request early and know that its contents are not ready for review or to be merged, add "[WIP]" at the start of the pull request title, which stands for "Work in Progress".
When you are happy with it and are happy for it to be merged into the main repository, change the "[WIP]" in the title of the pull request to "[Ready for review]".

A member of the team will then review your changes to confirm that they can be merged into the main repository.
A [review][github-review] will probably consist of a few questions to help clarify the work you've done.
Keep an eye on your GitHub notifications and be prepared to join in that conversation.

You can update your [fork][github-fork] of the [repository][monitoring-ecosystem-resilience-repo] and the pull request will automatically update with those changes.
You don't need to submit a new pull request when you make a change in response to a review.

You can also submit pull requests to other contributors' branches!
Do you see an [open pull request](https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/pulls) that you find interesting and want to contribute to?
Simply make your edits on their files and open a pull request to their branch!

What happens if the continuous integration (CI) fails (for example, if the pull request notifies you that "Some checks were not successful")?
The CI could fail for a number of reasons.
At the bottom of the pull request, where it says whether your build passed or failed, you can click “Details” next to the test, which takes you to the Travis page.
You can view the log or rerun the checks if you have write access to the repo by clicking the “Restart build” button in the top right (you must be logged in to Travis CI with your GitHub account see the “Restart build” button).

GitHub has a [nice introduction][github-flow] to the pull request workflow, but please get in touch if you have any questions.

## Style Guide

Docstrings should follow [numpydoc][link_numpydoc] convention.
We encourage extensive documentation.

The python code itself should follow [PEP8][link_pep8] convention whenever possible, with at most about 500 lines of code (not including docstrings) per script.

---

_These Contributing Guidelines have been adapted from the [Contributing Guidelines](https://github.com/bids-standard/bids-starter-kit/blob/master/CONTRIBUTING.md) of [The Turing Way](https://github.com/alan-turing-institute/the-turing-way)! (License: MIT)_

[monitoring-ecosystem-resilience-repo]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/
[monitoring-ecosystem-resilience-issues]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/issues
[git]: https://git-scm.com
[github]: https://github.com
[github-branches]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository
[github-fork]: https://help.github.com/articles/fork-a-repo
[github-flow]: https://guides.github.com/introduction/flow
[github-mergeconflicts]: https://help.github.com/articles/about-merge-conflicts
[github-pullrequest]: https://help.github.com/articles/creating-a-pull-request
[github-review]: https://help.github.com/articles/about-pull-request-reviews
[github-syncfork]: https://help.github.com/articles/syncing-a-fork
[labels-bug]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/bug
[labels-enhancement]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/enhancement
[labels-firstissue]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/good%20first%20issue
[labels-helpwanted]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/help%20wanted
[labels-project-management]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/project%20management
[labels-question]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/question
[link_numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[link_pep8]: https://www.python.org/dev/peps/pep-0008/
# Running `pyveg` on the cloud - Microsoft Azure

It is possible to make use of Azure cloud infrastructur in two ways when running pyveg:
* Using Azure blob storage to store downloaded images and results of the image processing.
* Using Azure batch to parallelize the running of the image processing.  This can vastly speed up the running time of your job.

In order to do these, you will need the following:
* An Azure account and an active subscription.
* An Azure "Blob Storage Account" - follow instructions on https://docs.microsoft.com/en-us/azure/storage/blobs/storage-blob-create-account-block-blob?tabs=azure-portal
* An Azure "Batch Account" - see https://docs.microsoft.com/en-us/azure/batch/batch-technical-overview for a description.  When setting up the batch account, it will ask you to link it to a Storage Account, so use the one above (and ensure that you select the same "Region" for both.
* You will probably need to increase the "quota" on your Batch Account - I believe that by default the quota for different types of VM are all set to zero.  There are instructions on how to do this at: https://docs.microsoft.com/en-us/azure/batch/batch-quota-limit - for our workflow we are using 100 dedicated cores of A1v2 VMs.

## Setting up `azure_config.py`

In the `pyveg/` directory there is a file `azure_config_template.py`.  Copy this to `azure_config.py` (in the same directory, and then start filling in the various fields.  The necessary info can be found in the Azure portal [https://portal.azure.com/#home] - perhaps the easiest way is to navigate via the "Subscriptions" icon at the top of the portal, then find the "Resources" (i.e. the Storage Account and the Batch Account).
* For the Storage Account, look on the left sidebar under "Settings" for the "Access keys" - then copy/paste one of the keys into the relevant field in `azure_config.py`.
* For the Batch Account, similarly there is a "Keys" icon under "Settings" which will lead to the relevant info.
For the "batch_pool_id" you can put any name you like - if there will be multiple people using the same batch and storage accounts, you might want to use your initials or something to identify you (and in this case, you should be careful that the sum of everyone's "node_count"s don't exceed the quota for the batch account.

Once you have populated the fields in `azure_config.py`, then do
```
pip install .
```
from the main `monitoring-ecosystem-resilience` directory.

## config settings to use Azure storage

We recommend that you create a new config file in the `pyveg/configs/` directory for every location/collection that you run.  You can use [pyveg/configs/test_Sentinel2_azure.py] as an example (you will want to increase the "date_range", and "n_sub_images" for the `NetworkCentralityCalculator` before running production though).

The key setting to tell the job to use Azure blob storage is:
```
output_location_type = "azure"
```

### config settings to use Azure batch

Note that using Azure storage as detailed above is a prerequisite for using Azure batch.
Again there is an example template [pyveg/configs/test_Sentinel2_batch.py] that you can copy and modify with your own coordinates, date range, collection etc.

Here, the important settings that enable the time-consuming parts of the pipeline to use Azure batch are in the `special_config` dictionary at the bottom of the file:
```
special_config = {
    "VegetationImageProcessor": {"run_mode": "batch"},
    "NetworkCentralityCalculator": {
        "n_sub_images": 10,
        "n_threads": 1,
        "run_mode": "batch",
    },
    "NDVICalculator": {"run_mode": "batch"},
}
```
This is setting "run_mode" to "batch" for these three Modules.  Note also that for `NetworkCentralityCalculator` we are setting "n_threads" to 1.  This is advisable if using the A1v2 nodes in the batch pool, as they only have a single vCPU per node.   If you instead choose more powerful VMs (which are also more expensive!) you can increase this accordingly.


### Checking the status of your job on Azure batch

If one or more Modules have "run_mode" set to "batch", the console output should give a running status of how many "tasks" are still running.  You can also check on the Azure portal [https://portal.azure.com/#home] - find the "Resource" that is your batch account, then on the left sidebar, navigate down to "Jobs", and click on the current job to see the status of all its "Tasks".  (A "Job" corresponds to an instance of a Module running over all the dates in the date range.  A "Task" is a single date within this date range.)

### Downloading data from Azure storage when it is ready

Azure blob storage is structured with "Containers" containing "Blobs", where the blobs themselves can have a directory-like structure.
A single pyveg pipeline job will produce a single Container, which will have the name
`<output_location>_<date_stamp>` where `output_location` was defined in your config file, and the `date_stamp` is when the job was launched.
You can find the container name in the logfile for your job, or via the Azure portal [https://portal.azure.com/#home] - if you find the "Resource" that is your Storage Account, then click on "Containers".

A script exists that can download the RGB images and the `results_summary.json` to a single local zipfile: `pyveg/scripts/download_from_azure.py`.   To run this, do
```
pyveg_azure_download --container <container_name> --output_zipfile <name_of_zipfile_to_write_to>
```


## What is going on "under the hood" when running on Azure batch?

(This section is only necessary if you are interested in knowing more about how this works - if you just want to run the jobs, the instructions above should suffice.)

When you run the command
```
pyveg_run_pipeline --config_file <some_config_file>
```
the script `pyveg/scripts/run_pyveg_pipeline` will read in the specified config file and use it to set parameters of a "Pipeline" that is composed of "Sequences", which are in turn composed of "Modules".

The batch functionality is implemented at the Module level.  If a Module has "run_mode" set to "batch", it will:
* Get a dictionary of batch Tasks on which its Tasks depends (e.g. the NetworkCentralityCalculator needs the ImageProcessor to have finished for a given date before it can run that date's Task).
* Creates a new batch "Job" with a name composed of the Module name and the current time.
* Create a dictionary of batch Tasks for the Job, dividing up the date range amongst the Tasks, and storing the configuration of the Module for each entry.
* Upload the `azure_config.py` and `pyveg/scripts/batch_commands.sh` to blob storage.

For each Task, the process is then:
* Write the configuration to a JSON file and upload to blob storage.
* Submit the batch Task, which will run `batch_commands.sh` on the batch node.

### What does `batch_commands.sh` do ?

The execution of a single Task on a batch node (which is an Ubuntu-16.04 VM) is governed by this shell script `batch_commands.sh`.  The basic flow is:
* Install some packages, including miniconda, and create and activate a Python 3.7 conda environment.
* Clone the `monitoring-ecosystem-resilience` repo, change to the `develop` branch, and do ```pip install .``` to install it.
* Run the command ```pyveg_run_module --config_file <json_config_for_module>``` where the json config file is the set of parameters needed to configure this Module, based on the dictionary created in `processor_modules.creta_task_dict`.

### What happens when all tasks are submitted?

The function `processor_modules.run_batch` will submit all the tasks and then return straightaway, rather than waiting for the tasks to finish.  This means that other Modules or Sequences (e.g. WeatherSequence) that do not depend on the results of this Module can still be executed.
However, usually the final Sequence in a Pipeline will be a "combiner" Sequence, that has a `depends_on` attribute.
If a Sequence listed in `depends_on` has one-or-more Modules with "run_mode" set to "batch", the logic in the `run` method of the `Sequence` class in `pyveg/src/pyveg_pipeline.py` will loop through all the Modules in that Sequence, and call `check_if_finished()` on all of them.  This in turn will query the Batch Job to see the status of all the Tasks.# Uploading results to the Zenodo open source data repository

Zenodo ([https://zenodo.org]) is a free and open source repository for research data, hosted by CERN.

Data is organized into `depositions`, each of which has a Digital Object Identifier (DOI) which can then be cited.

For the purposes of this package, we make the assumption that we will keep all the data for a single paper in one deposition.

## Prerequisites

In order to use the functions in this package to automatically upload data to Zenodo, and to download specific files to rerun analysis, you will need:
* Sign up for a Zenodo account by clicking the "Sign up" button on the top right of [https://zenodo.org/]
* If you want to use the "sandbox" repository for testing functionality (recommended!) you'll also need to sign up separately here: [https://sandbox.zenodo.org/]
* For both the production and sandbox versions, once you are signed in, go to [https://zenodo.org/account/settings/applications/tokens/new/] to create an API token.  Write any name for your token in the box, and tick the boxes for "deposit:actions" and "deposit:write" before clicking the "Create" button.  Keep this tab open until you have copied the token into `zenodo_config.py` (see below).

## How to fill `zenodo_config.py`

In the `pyveg/` directory there is a file `zenodo_config_template.py`.  Copy this to `zenodo_config.py` and fill in the various fields:
* The `metadata_dict` is the metadata that will be stored with your deposition.  Put the title and description of your study here, and "upload_type" as "dataset".  List the authors, giving names and affiliations as you would like them to appear on Zenodo.
* For the "test_api_credentials", if you plan to use the "sandbox" repository for testing, and if you have signed up for this and created an API token as described in the section above, copy/paste the personal access token to the "api_token" field here.
* Similarly for the "prod_api_credentials" do the same, but with the main Zenodo site.
* Leave the "deposition_id" as None for now - we will create a deposition in the next step.

## Create a deposition to hold our data

Once we have filled in the "api_token" in `zenodo_config.py` we can do:
```
pip install .
```
then you can run the following command to create a new deposition:
```
pyveg_zenodo_upload --create_deposition [--test_api]
```
where the final `--test_api` argument should be included if you want to use the sandbox repository, or omitted to use the production one.
The output from this command should give you the deposition_id that you can then paste into the appropriate section of `zenodo_config.py`, and then do
```
pip install .
```
once more.

## Uploading analysis results to Zenodo

Once all the necessary fields (the "api_token" and "deposition_id") are present in `zenodo_config.py`, then uploading analysis results, plus the "results_summary.json" file (the output of the image downloading and processing that is the input to the analysis) should be straightforward:
* When running ```pyveg_gee_analysis``` you can add the argument ```--upload_to_zenodo``` to upload the results to the deposition on the production Zenodo repository, or ```--upload_to_zenodo_sandbox``` to use the sandbox repository instead.
* Alternatively, if you have previously run ```pyveg_gee_analysis``` you can run the command:
```
pyveg_zenodo_upload --input_png_loc <path-to-analysis-subdir> --input_json_loc <path-to-dir-containing-results_summary.json> --json_loc_type <'local' or 'azure'> --collection <collection-name>
```
Here, we assume that the png files from running the analysis are in a local directory, while the "results_summary.json" can be either in a local directory (in which case specify this directory as the ```--input_json_loc``` argument and specify ```--json_loc_type local```), or on Azure (in which case use the blob storage container as the ```--input_json_loc``` argument and specify ```--json_loc_type azure```)
# The `pyveg` Package

## Introduction 

The `pyveg` package is developed to study the evolution of vegetation patterns in semi-arid environments using data downloaded from Google Earth Engine.

The code in this repository is intended to perform two main tasks:

**1. Download and process GEE data**:

* Download satellite data from Google Earth Engine (images and weather data).
    * Downloaded images are divided into 50x50 pixel sub-images, network centrality metrics are used to describe the pattern vegetation are then calculated on the sub-image level. Both colour (RGB) and Normalised Difference Vegetation Index (NDVI) images are downloaded and stored on the sub-image level. 
    * For weather collections the precipitation and temperature "images" are averaged into a single value at each point in the time series.
* The download job is fully specified by a configuration file that can be generated by specifying the details of the data to be downloaded via prompts (satellite to use, coordinates, time period, number of time points, etc.).  

**2. Time series analysis on downloaded data**:

* Time series analysis of the following metrics: raw NDVI mean pixel intensity across the image, vegetation network centrality metric, and precipitation.
    * The time series are processed (outliers removed and resampled to avoid gaps). All time series of each sub-image are aggregated into one summary time series that is used for analysis. 
    * The summary time series is de-seasonalised and smoothed.
    * Residuals between the raw and de-seasonalised and smoothed time series are calculated and used for an early warning resilience analysis.

* Time series plots are produced, along with auto- and cross-correlation plots. Early warning signals are also computed using the [ewstools package](https://github.com/ThomasMBury/ewstools), including Lag-1 autocorrelation and standard deviation moving window plots. A sensitivity and significance analysis is also performed in order to determine whether any trends are statistically significant.
* Time series summary statistics and resilience metrics are saved into files.
* A PDF report is created showcasing the main figures resulting from the analyses. 

**Other functionalities**:

`pyveg` also has other minor functionalities:

* Analysis of a collection of summary data that has been created with the `pyveg` pipeline (downloading + time series analysis).
* Simulate the generation and evolution of patterned vegetation
* A stand-alone network centrality estimation for a 50x50 pixel image.
* A functionality to upload results to the Zenodo open source repository 

### `pyveg` flow

The diagram below represents the high level flow of the main functionalities of the `pyveg` package. For each main component there is a CLI console scripts defined, that is shown in the diagram. 

![The`pyveg` program flow.](paper/pveg_flow.png)

The full ReadTheDocs documentation for this `pyveg` can be found in this [link](https://pyveg.readthedocs.io/en/latest/).


This page contains an installation guide, and some usage examples for this package.


## Installation

`pyveg` requires Python 3.6 or greater. To install, start by creating a fresh `conda` environment.
```
conda create -n veg python=3.7
conda activate veg
```
Get the source.
```
git clone https://github.com/alan-turing-institute/monitoring-ecosystem-resilience.git
```
Enter the repository and check out a relevant branch if necessary (the default `master` branch contains the most up to date stable version of the code).
```
cd monitoring-ecosystem-resilience
```
Install the package using `pip`.
```
pip install .
```
If you are using Windows and encounter issues during this stage, a solution may be found here: https://github.com/NREL/OpenOA/issues/37. If you plan on making changes to the source code, you can instead run `pip install -e .`. 

Before using the Google Earth Engine API, you need to sign up with a Google account [here](https://earthengine.google.com/new_signup/), and authenticate.  To authenticate, run
```
earthengine authenticate
```
A new browser window will open. Copy the token from this window to the terminal prompt to complete the authentication process.


### Google Earth Engine

[Google Earth Engine](https://earthengine.google.com) (GEE) is a powerful tool for obtaining and analysing satellite imagery. This directory contains some useful scripts for interacting with the Earth Engine API. The earth engine API is installed automatically as part of the `pyveg` package installation. If you wish to install it separately, you can follow the instructions [here](https://developers.google.com/earth-engine/python_install_manual).

## Downloading data from GEE with ``pyveg``

### Downloading data from GEE using the CLI

To run a `pyveg` download job, use
```
pyveg_run_pipeline --config_file <path to config>
```

The download job is fully specified by a configuration file, which you point to using the `--config_file` argument. A sample config file is found at `pyveg/configs/config_all.py`. You can also optionally specify a string to identify the download job using the `--name` argument.

Note that we use the GEE convention for coordinates, i.e. `(longitude,latitude)`.

#### Generating a download configuration file 

To create a configuration file for use in the pyveg pipeline described above, use the command 
```
pyveg_generate_config
```
this allows the user to specify various characteristics of the data they want to download via prompts. The list in order is as follows:

* `--configs_dir`: The path to the directory containing the config file, with a default option `pyveg/configs`.

* `--collection_name`: The name of the dataset used in the collection, either Sentinel2, or Landsat 8, 7, 5 or 4.
    *    Sentinel2: [Available from 2015-06-23 at 10m resolution.](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2)
    *    Landsat8: [Available from 2013-04-11 at 30m resolution.](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1)
    *    Landsat7: [Available from 1999-01-01 at 30m resolution.](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1)
    *    Landsat5: [Available from 1984-03-10 to 2013-01-31 at 60m resolution.](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C01_T1)
    *    Landsat4: [Available from 1982-07-16 to 1993-12-14 at 60m resolution.](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C01_T1)

* `--latitude`: The latitude (in degrees north) of the centre point of the image collection.

* `--longitude`: The longitude (in degrees east) of the centre point of the image collection.

* `--country`: The country (for the file name) can either be entered, or use the specified coordinates to look up the country name from the OpenCage database.

* `--start_date`: The start date in the format ‘YYYY-MM-DD’, the default is ‘2015-01-01’ (or ‘2019-01-01’ for a test config file).

* `--end_date`: The end date in the format ‘YYYY-MM-DD’, the default is today’s date (or ‘2019-03-01’ for a test config file).

* `--time_per_point`: The option to run the image collection either monthly (‘1m’) or weekly (‘1w’), with the default being monthly.

* `--run_mode`: The option to run time-consuming functions on Azure (‘batch’) or running locally on your own computer (‘local’). The default is local. For info about running on Azure go [here](UsingAzure.md).

* `--output_dir`: The option to write the output to a specified directory, with the default being the current directory.

* `--test_mode`: The option to make a test config file, containing fewer months and a subset of sub-images, with a default option to have a normal config file.
    *    By choosing the test config file, the start and end dates (see below) are defaulted to cover a smaller time span.
    *    It is recommended that the test config option should be used purely to determine if the options specified by the user are correct.


* `--n_threads`:  Finally, how many threads the user would like to use for the time-consuming processes, either 4 (default) or 8.

For example:
```
 pyveg_generate_config --configs_dir "pyveg/configs" --collection_name "Sentinel2" --latitude 11.58 --longitude 27.94 --start_date "2016-01-01" --end_date "2020-06-30" --time_per_point "1m" --run_mode "local" --n_threads 4
```

This generates a file named `config_Sentinel2_11.58N_27.94E_Sudan_2016-01-01_2020-06-30_1m_local.py` along with instructions on how to use this configuration file to download data through the pipeline, in this case the following:

```
pyveg_run_pipeline --config_file pyveg/configs/config_Sentinel2_11.58N_27.94E_Sudan_2016-01-01_2020-06-30_1m_local.py
```

Individual options can be specified by the user via prompt. The options for this can be found by typing ```pyveg_generate_config --help```. 


### More Details on Downloading

During the download job, `pyveg` will break up your specified date range into a time series, and download data at each point in the series. Note that by default the vegetation images downloaded from GEE will be split up into 50x50 pixel images, vegetation metrics are then calculated on the sub-image level. Both colour (RGB) and Normalised Difference Vegetation Index (NDVI) images are downloaded and stored. Vegetation metrics include the mean NDVI pixel intensity across sub-images, and also network centrality metrics, discussed in more detail below.

For weather collections e.g. the ERA5, due to coarser resolution, the precipitation and temperature "images" are averaged into a single value at each point in the time series.

### Rerunning partially succeeded jobs

The output location of a download job is datestamped with the time that the job was launched.  The configuration file used will also be copied and datestamped, to aid reproducibility.  For example if you run the job
```
pyveg_run_pipeline --config_file pyveg/configs/my_config.py
```
there will be a copy of `my_config.py` saved as `pyveg/configs/cached_configs/my_config_<datestamp>.py`. This also means that if a job crashes or timeouts partway through, it is possible to rerun, writing to the same output location and skipping parts that are already done by using this cached config file.  However, in order to avoid a second datestamp being appended to the output location, use the ```--from_cache``` argument.  So for the above example, the command to rerun the job filling in any failed/incomplete parts would be:
```
pyveg_run_pipeline --config_file pyveg/configs/cached_configs/my_config_<datestamp>.py --from_cache
```

### Using Azure for downloading/processing data

If you have access to Microsoft Azure cloud computing facilities, downloading and processing data can be sped up enormously by using batch computing to run many subjobs in parallel.  See [here](UsingAzure.md) for more details.

### Downloading data using the API

Although `pyveg` has been mostly designed to be used with the CLI as shown above, we can also use `pyveg` functions through the API. A tutorial of how to download data this way is included in the ``notebooks/tutorial_download_and_process_gee_images.ipynb`` notebook tutorial. 

## Analysing the Downloaded Data with `pyveg`

### Analysing the Downloaded Data using the CLI

Once you have downloaded the data from GEE, the `pyveg_gee_analysis` command allows you to process and analyse the output. To run:
```
pyveg_gee_analysis --input_dir <path_to_pyveg_download_output_dir>
```
The analysis code preprocesses the data and produces a number of plots. These will be saved in an `analysis/` subdirectory inside the `<path_to_pyveg_download_output_dir>` directory.

Note that in order to have a meaningful analysis, the dowloaded time series should have at least 4 points (and more thant 12
for an early warning analysis) and not being the result of a "test" config file, in this case the analysis fails.
 
The commands also allows for other options to be added to the execution of the script (e.g. run analysis from a downloaded data in Azure blob storage, define a different output directly, don't include a time series analysis, etc), which can be displayed by typing:

```
pyveg_gee_analysis --help
```

The analysis script executed with the ```pyveg_gee_analysis``` command runs the following steps:

#### Preprocessing

`pyevg` supports the following pre-processing operations:
- Identify and remove outliers from the time series.
- Fill missing values in the time series (based on a seasonal average), or resample the time series using linear interpolation between points.
- Smoothing of the time series using a LOESS smoother.
- Calculation of residuals between the raw and smoothed time series.
- De-seasonalising (using first differencing), and detrending using STL.

#### Plots

In the `analysis/` subdirectory, `pyveg` creates the following plots:
- Time series plots containing vegetation and precipitation time series (seasonal and de-seasonalised). Plots are labelled with the AR1 of the vegetation time series, and the maximum correlation between the Vegetation and precipitation time series.
- Auto-correlation plots for vegetation and precipitation time series (seasonal and de-seasonalised).
- Vegetation and precipitation cross-correlation scatterplot matrices.
- STL decomposition plots.
- Resilience analysis:
     - `ewstools` resilience plots showing AR1, standard deviation, skewness, and kurtosis using a moving window.
     - Smoothing filter size and moving window size Kendall tau sensitivity plots.
     - Significance test.

### Running the analysis using the API

Although `pyveg` has been mostly designed to be used with the CLI as shown above, we can also use `pyveg` functions through the API. A tutorial of how to run the data analysis in this way is included in the ```notebooks/tutorial_analyse_gee_data.ipynb``` notebook  tutorial. 


## Other functionalities of `pyveg`

### Analysis summary statistics data 

The ```analyse_pyveg_summary_data.py``` functionality processes collections of data produced by the main `pyveg` pipeline described in the section above (download + time series analysis for different locations and time periods) and creates a number of plots of the summary statistics of these time series.

To run this analysis in Python, there is an entrypoint defined.  Type:

```
pyveg_analysis_summary_data --input_location  <path_to_directory_with_collection_summary_statistics>
```

if you wish you can also specify the ```outpur_dir``` where plots will be saved.  Type:
```
pyveg_analysis_summary_data --help
```
to see the extra options.


### Pattern simulation

The ```generate_patterns.py``` functionality originates from some Matlab code by Stefan Dekker, Willem Bouten, Maarten Boerlijst and Max Rietkerk (included in the "matlab" directory), implementing the scheme described in:

Rietkerk et al. 2002. Self-organization of vegetation in arid ecosystems. The American Naturalist 160(4): 524-530.

To run this simulation in Python, there is an entrypoint defined.  Type:
```
pyveg_gen_pattern --help
```
to see the options.  The most useful option is the `--rainfall` parameter which sets a parameter (the rainfall in mm) of the simulation - values between 1.2 and 1.5 seem to give rise to a good range of patterns. Other optional parameters for `generate_patterns.py` allow the generated image to be output as a csv file or a png image.  The `--transpose` option rotates the image 90 degrees (this was useful for comparing the Python and Matlab network-modelling code). Other parameters for running the simulation are in the file `patter_gen_config.py`, you are free to change them.

#### Running the pattern simulation using the API

Although `pyveg` has been mostly designed to be used with the CLI as shown above, we can also use `pyveg` functions through the API. A tutorial of how to run the simulation of the pattern generation in this way is included in [here](notebooks/tutorial_simulate_patterned_vegetation.ipynb). 


### Network centrality

There is an entrypoint defined in `setup.py` that runs the *main* function of `calc_euler_characteristic.py`:
```
pyveg_calc_EC --help
```
will show the options.

* `--input_txt` allows you to give the input image as a csv, with one row per row of pixels.  Inputs are expected to be "binary", only containing two possible pixel values (typically 0 for black and 255 for white).
* `--input_img` allows you to pass an input image (png or tif work OK).  Note again that input images are expected to be "binary", i.e. only have two colours.
* `--sig_threshold` (default value 255) is the value above (or below) which a pixel is counted as signal (or background)
* `--upper_threshold` determines whether the threshold above is an upper or lower threshold (default is to have a lower threshold - pixels are counted as "signal" if their value is greater-than-or-equal-to the threshold value).
* `--use_diagonal_neighbours` when calculating the adjacency matrix, the default is to use "4-neighbours" (i.e. pixels immediately above, below, left, or right).  Setting this option will lead to "8-neighbours" (i.e. the four neighbours plus those diagonally adjacent) to be included.
* `--num_quantiles` determines how many elements the output feature vector will have.
* `--do_EC` Calculate the Euler Characteristic to fill the feature vector.  Currently this is required, as the alternative approach (looking at the number of connected components) is not fully debugged.

Note that if you use the `-i` flag when running python, you will end up in an interactive python session, and have access to the `feature_vec`, `sel_pixels` and `sc_images` variables.

Examples:
```
pyveg_calc_EC --input_txt ../binary_image.txt --do_EC
>>> sc_images[50].show() # plot the image with the top 50% of pixels (ordered by subgraph centrality) highlighted.
>>> plt.plot(list(sel_pixels.keys()), feature_vec, "bo") # plot the feature vector vs pixel rank
>>> plt.show()
```

### Uploading results to the Zenodo open source repository

See [here](UsingZenodo.md) for more details.

# Contributing 

We welcome contributions from anyone who is interested in the project. There are lots of ways to contribute, not just writing code. See our [Contributor Guidelines](CONTRIBUTING.md) to learn  more about how you can contribute and how we work together as a community.

# Licence

This project is licensed under the terms of the MIT software license.
---
title: 'pyveg: A Python package for analysing the time evolution of patterned vegetation using Google Earth Engine'
tags:
  - Python
  - Ecology
  - Remote sensing
  - Time Series Analysis
  - Early warnings
authors:
  - name: Nick Barlow
    affiliation: 1
  - name: Camila Rangel Smith
    affiliation: 1
  - name: Samuel Van Stroud
    affiliation: 1, 2
  - name: Jesse F. Abrams
    affiliation: 3
  - name: Chris A. Boulton
    affiliation: 3
  - name: Joshua Buxton
    affiliation: 3
affiliations:
 - name: The Alan Turing Institute
   index: 1
 - name: University College London
   index: 2
 - name: University of Exeter
   index: 3
date: 09 October 2020
bibliography: paper.bib
---

# Introduction

Periodic vegetation patterns (PVP) arise from the interplay between
forces that drive the growth and mortality of plants. Inter-plant
competition for resources, in particular water, can lead to the
formation of PVP. Arid and semi-arid ecosystems may be under threat
due to changing precipitation dynamics driven by macroscopic changes
in climate. These regions display some noteable examples of PVP,
for example the "tiger bush" patterns found in West Africa.

The morphology of the periodic pattern has been suggested to be
linked to the resilience of the ecosystem [@Mander:2017; @Trichon:2018].
Using remote sensing techniques,  vegetation patterns in these regions
can be studied, and an analysis of the resilience of the ecosystem can
be performed.

The `pyveg` package implements functionality to download and process data
from Google Earth Engine (GEE), and to subsequently perform a
resilience analysis on the aquired data. PVP images are quantified using
network centrality metrics. The results of the analysis can be used
to search for typical early warning signals of an ecological collapse
[@Dakos:2008]. Google Earth Engine Editor scripts are also provided to help
researchers discover locations of ecosystems which may be in
decline.

`pyveg` is being developed as part of a research project
looking for evidence of early warning signals of ecosystem
collapse using remote sensing data. `pyveg` allows such
research to be carried out at scale, and hence can be an
important tool in understanding changing arid and semi-arid
ecosystem dynamics. An evolving list of PVP locations, obtained through
both literature and manual searches, is included in the package at
`pyveg/coordinates.py`. The structure of the package is outlined in
\autoref{fig:pyveg_flow}, and is discussed in more detail in the
following sections.

![`pyveg` program flow.\label{fig:pyveg_flow}](pveg_flow.png)


# Downloading data from Google Earth Engine

In order to interact with the GEE API, the user must sign up to GEE
and obtain an API key, which is linked to a Google account. Upon downloading
data using `pyveg` for the first time, the
user will be prompted to enter their API key to authenticate GEE. The `run_pyveg_pipeline`
command initiates the downloading of time series data at a single
coordinate location. The job is configured using a configuration file
specified by the `--config_file` argument.

Within the configuration file, the user can specify the following:
coordinates of the download location, start and end dates of the
time series, frequency with which to sample, choice of GEE collections
to download from (currently vegetation and precipitation collections are
supported).

`pyveg` will then form a series of date ranges, and query GEE for the relevant
data in each date range. Colour (RGB) and Normalised Difference vegetation
Index (NDVI) images are downloaded from vegetation collections. Supported 
vegetation collections include Landsat [@landsat] and Sentinel-2 [@sentinel2] GEE
collections. Cloud masking
logic is included to improve data quality using the `geetools` package [@geetools].
For precipitation and temperature information, `pyveg` defaults to using the ERA5
collection [@era5].


# Network centrality metrics

Network centrality methods are used to measure the connectedness of vegetation
in images by treating the image as a network, with pixels containing significant
vegetation as nodes. Vegetation pixels are ordered according to their subgraph
centrality [@PhysRevE.71.056103], and from this, a feature vector is constructed
by calculating the Euler Characteristic [@richeson2012euler] for different quantiles.
The slope of this feature vector gives a measure of how connected the vegetation is.

After completetion of the download job, `pyveg` computes the network centrality
of the vegetation [@Mander:2017]. To achieve this, the NDVI image is broken up
into smaller $50 \times 50$ pixel sub-images. Each sub-image is then thresholded
using the NDVI pixel intensity, and subgraph connectivity is computed for each
binarized sub-image. The resulting metrics are stored, along with mean NDVI pixel
intensities for each sub-image.


# Time series analysis

`pyveg` analysis functionality is exposed via a `pveg_gee_analysis` command.
The command accepts an argument, `--input_dir`, which points to a directory
previously created by a download job. It is also possible to run this analysis
on data within Azure using `--input_container`. Users are able to upload data to Zenodo
with `pyveg` and to analyse data hosted on Zenodo using `--input_zenodo_coords` argument
of the `pveg_gee_analysis` command. `pyveg` supports the analysis of the
following time series: raw NDVI mean pixel intensity across the image, offset50
(a measure of the slope of the network centrality feature vector), and precipitation.

During data processing, `pyveg` is able
to drop time series outliers and resample the time series to clean the data
and avoid gaps. A smoothed time series is constructed using LOESS smoothing,
and residuals between the raw and smoothed time series are calculated.
Additionally, a deseasonalised time series is constructed via the first
difference method.

Time series plots are produced, along with auto- and cross-correlation plots.
Early warning signals are also computed using the `ewstools` package [@ewstools],
including Lag-1 autocorrelation and standard deviation moving window plots.
A sensitivity and significance analysis is also performed in order to determine
whether any trends (quantified by Kendall tau values) are statistically significant.
The vegetation decay rate is calculated by fitting a crystall ball function
to the annual average offset50 and NDVI time series [@CrystalBallFunction].

Following data processing, `pyveg` is able to calculate summary plots using
`pyveg_analysis_summary_data`. This uses as input a collection of summary statistics
extracted from the time series of each individual location obtained with the download
and analysis functionality described above. These are hosted locally or on Zenodo.


# Acknowledgements

The `pyveg` package was developed by researchers from the Alan Turing Institute,
University College London, and the University of Exeter.  Funding was provided by
the Alan Turing Institute, the Science and Technology Facilities Council, and the
Leverhulme Trust (grant number RPG-2018-046).
We would like to acknowledge support from Tim Lenton during the course of
this project.


# References
# Running `pyveg` on the cloud - Microsoft Azure

It is possible to make use of Azure cloud infrastructur in two ways when running pyveg:
* Using Azure blob storage to store downloaded images and results of the image processing.
* Using Azure batch to parallelize the running of the image processing.  This can vastly speed up the running time of your job.

In order to do these, you will need the following:
* An Azure account and an active subscription.
* An Azure "Blob Storage Account" - follow instructions on https://docs.microsoft.com/en-us/azure/storage/blobs/storage-blob-create-account-block-blob?tabs=azure-portal
* An Azure "Batch Account" - see https://docs.microsoft.com/en-us/azure/batch/batch-technical-overview for a description.  When setting up the batch account, it will ask you to link it to a Storage Account, so use the one above (and ensure that you select the same "Region" for both.
* You will probably need to increase the "quota" on your Batch Account - I believe that by default the quota for different types of VM are all set to zero.  There are instructions on how to do this at: https://docs.microsoft.com/en-us/azure/batch/batch-quota-limit - for our workflow we are using 100 dedicated cores of A1v2 VMs.

## Setting up `azure_config.py`

In the `pyveg/` directory there is a file `azure_config_template.py`.  Copy this to `azure_config.py` (in the same directory, and then start filling in the various fields.  The necessary info can be found in the Azure portal [https://portal.azure.com/#home] - perhaps the easiest way is to navigate via the "Subscriptions" icon at the top of the portal, then find the "Resources" (i.e. the Storage Account and the Batch Account).
* For the Storage Account, look on the left sidebar under "Settings" for the "Access keys" - then copy/paste one of the keys into the relevant field in `azure_config.py`.
* For the Batch Account, similarly there is a "Keys" icon under "Settings" which will lead to the relevant info.
For the "batch_pool_id" you can put any name you like - if there will be multiple people using the same batch and storage accounts, you might want to use your initials or something to identify you (and in this case, you should be careful that the sum of everyone's "node_count"s don't exceed the quota for the batch account.

Once you have populated the fields in `azure_config.py`, then do
```
pip install .
```
from the main `monitoring-ecosystem-resilience` directory.

## config settings to use Azure storage

We recommend that you create a new config file in the `pyveg/configs/` directory for every location/collection that you run.  You can use [pyveg/configs/test_Sentinel2_azure.py] as an example (you will want to increase the "date_range", and "n_sub_images" for the `NetworkCentralityCalculator` before running production though).

The key setting to tell the job to use Azure blob storage is:
```
output_location_type = "azure"
```

### config settings to use Azure batch

Note that using Azure storage as detailed above is a prerequisite for using Azure batch.
Again there is an example template [pyveg/configs/test_Sentinel2_batch.py] that you can copy and modify with your own coordinates, date range, collection etc.

Here, the important settings that enable the time-consuming parts of the pipeline to use Azure batch are in the `special_config` dictionary at the bottom of the file:
```
special_config = {
    "VegetationImageProcessor": {"run_mode": "batch"},
    "NetworkCentralityCalculator": {
        "n_sub_images": 10,
        "n_threads": 1,
        "run_mode": "batch",
    },
    "NDVICalculator": {"run_mode": "batch"},
}
```
This is setting "run_mode" to "batch" for these three Modules.  Note also that for `NetworkCentralityCalculator` we are setting "n_threads" to 1.  This is advisable if using the A1v2 nodes in the batch pool, as they only have a single vCPU per node.   If you instead choose more powerful VMs (which are also more expensive!) you can increase this accordingly.


### Checking the status of your job on Azure batch

If one or more Modules have "run_mode" set to "batch", the console output should give a running status of how many "tasks" are still running.  You can also check on the Azure portal [https://portal.azure.com/#home] - find the "Resource" that is your batch account, then on the left sidebar, navigate down to "Jobs", and click on the current job to see the status of all its "Tasks".  (A "Job" corresponds to an instance of a Module running over all the dates in the date range.  A "Task" is a single date within this date range.)

### Downloading data from Azure storage when it is ready

Azure blob storage is structured with "Containers" containing "Blobs", where the blobs themselves can have a directory-like structure.
A single pyveg pipeline job will produce a single Container, which will have the name
`<output_location>_<date_stamp>` where `output_location` was defined in your config file, and the `date_stamp` is when the job was launched.
You can find the container name in the logfile for your job, or via the Azure portal [https://portal.azure.com/#home] - if you find the "Resource" that is your Storage Account, then click on "Containers".

A script exists that can download the RGB images and the `results_summary.json` to a single local zipfile: `pyveg/scripts/download_from_azure.py`.   To run this, do
```
pyveg_azure_download --container <container_name> --output_zipfile <name_of_zipfile_to_write_to>
```


## What is going on "under the hood" when running on Azure batch?

(This section is only necessary if you are interested in knowing more about how this works - if you just want to run the jobs, the instructions above should suffice.)

When you run the command
```
pyveg_run_pipeline --config_file <some_config_file>
```
the script `pyveg/scripts/run_pyveg_pipeline` will read in the specified config file and use it to set parameters of a "Pipeline" that is composed of "Sequences", which are in turn composed of "Modules".

The batch functionality is implemented at the Module level.  If a Module has "run_mode" set to "batch", it will:
* Get a dictionary of batch Tasks on which its Tasks depends (e.g. the NetworkCentralityCalculator needs the ImageProcessor to have finished for a given date before it can run that date's Task).
* Creates a new batch "Job" with a name composed of the Module name and the current time.
* Create a dictionary of batch Tasks for the Job, dividing up the date range amongst the Tasks, and storing the configuration of the Module for each entry.
* Upload the `azure_config.py` and `pyveg/scripts/batch_commands.sh` to blob storage.

For each Task, the process is then:
* Write the configuration to a JSON file and upload to blob storage.
* Submit the batch Task, which will run `batch_commands.sh` on the batch node.

### What does `batch_commands.sh` do ?

The execution of a single Task on a batch node (which is an Ubuntu-16.04 VM) is governed by this shell script `batch_commands.sh`.  The basic flow is:
* Install some packages, including miniconda, and create and activate a Python 3.7 conda environment.
* Clone the `monitoring-ecosystem-resilience` repo, change to the `develop` branch, and do ```pip install .``` to install it.
* Run the command ```pyveg_run_module --config_file <json_config_for_module>``` where the json config file is the set of parameters needed to configure this Module, based on the dictionary created in `processor_modules.creta_task_dict`.

### What happens when all tasks are submitted?

The function `processor_modules.run_batch` will submit all the tasks and then return straightaway, rather than waiting for the tasks to finish.  This means that other Modules or Sequences (e.g. WeatherSequence) that do not depend on the results of this Module can still be executed.
However, usually the final Sequence in a Pipeline will be a "combiner" Sequence, that has a `depends_on` attribute.
If a Sequence listed in `depends_on` has one-or-more Modules with "run_mode" set to "batch", the logic in the `run` method of the `Sequence` class in `pyveg/src/pyveg_pipeline.py` will loop through all the Modules in that Sequence, and call `check_if_finished()` on all of them.  This in turn will query the Batch Job to see the status of all the Tasks.# Uploading results to the Zenodo open source data repository

Zenodo ([https://zenodo.org]) is a free and open source repository for research data, hosted by CERN.

Data is organized into `depositions`, each of which has a Digital Object Identifier (DOI) which can then be cited.

For the purposes of this package, we make the assumption that we will keep all the data for a single paper in one deposition.

## Prerequisites

In order to use the functions in this package to automatically upload data to Zenodo, and to download specific files to rerun analysis, you will need:
* Sign up for a Zenodo account by clicking the "Sign up" button on the top right of [https://zenodo.org/]
* If you want to use the "sandbox" repository for testing functionality (recommended!) you'll also need to sign up separately here: [https://sandbox.zenodo.org/]
* For both the production and sandbox versions, once you are signed in, go to [https://zenodo.org/account/settings/applications/tokens/new/] to create an API token.  Write any name for your token in the box, and tick the boxes for "deposit:actions" and "deposit:write" before clicking the "Create" button.  Keep this tab open until you have copied the token into `zenodo_config.py` (see below).

## How to fill `zenodo_config.py`

In the `pyveg/` directory there is a file `zenodo_config_template.py`.  Copy this to `zenodo_config.py` and fill in the various fields:
* The `metadata_dict` is the metadata that will be stored with your deposition.  Put the title and description of your study here, and "upload_type" as "dataset".  List the authors, giving names and affiliations as you would like them to appear on Zenodo.
* For the "test_api_credentials", if you plan to use the "sandbox" repository for testing, and if you have signed up for this and created an API token as described in the section above, copy/paste the personal access token to the "api_token" field here.
* Similarly for the "prod_api_credentials" do the same, but with the main Zenodo site.
* Leave the "deposition_id" as None for now - we will create a deposition in the next step.

## Create a deposition to hold our data

Once we have filled in the "api_token" in `zenodo_config.py` we can do:
```
pip install .
```
then you can run the following command to create a new deposition:
```
pyveg_zenodo_upload --create_deposition [--test_api]
```
where the final `--test_api` argument should be included if you want to use the sandbox repository, or omitted to use the production one.
The output from this command should give you the deposition_id that you can then paste into the appropriate section of `zenodo_config.py`, and then do
```
pip install .
```
once more.

## Uploading analysis results to Zenodo

Once all the necessary fields (the "api_token" and "deposition_id") are present in `zenodo_config.py`, then uploading analysis results, plus the "results_summary.json" file (the output of the image downloading and processing that is the input to the analysis) should be straightforward:
* When running ```pyveg_gee_analysis``` you can add the argument ```--upload_to_zenodo``` to upload the results to the deposition on the production Zenodo repository, or ```--upload_to_zenodo_sandbox``` to use the sandbox repository instead.
* Alternatively, if you have previously run ```pyveg_gee_analysis``` you can run the command:
```
pyveg_zenodo_upload --input_png_loc <path-to-analysis-subdir> --input_json_loc <path-to-dir-containing-results_summary.json> --json_loc_type <'local' or 'azure'> --collection <collection-name>
```
Here, we assume that the png files from running the analysis are in a local directory, while the "results_summary.json" can be either in a local directory (in which case specify this directory as the ```--input_json_loc``` argument and specify ```--json_loc_type local```), or on Azure (in which case use the blob storage container as the ```--input_json_loc``` argument and specify ```--json_loc_type azure```)
# The `pyveg` Package

## Introduction 

The `pyveg` package is developed to study the evolution of vegetation patterns in semi-arid environments using data downloaded from Google Earth Engine.

The code in this repository is intended to perform two main tasks:

**1. Download and process GEE data**:

* Download satellite data from Google Earth Engine (images and weather data).
    * Downloaded images are divided into 50x50 pixel sub-images, network centrality metrics are used to describe the pattern vegetation are then calculated on the sub-image level. Both colour (RGB) and Normalised Difference Vegetation Index (NDVI) images are downloaded and stored on the sub-image level. 
    * For weather collections the precipitation and temperature "images" are averaged into a single value at each point in the time series.
* The download job is fully specified by a configuration file that can be generated by specifying the details of the data to be downloaded via prompts (satellite to use, coordinates, time period, number of time points, etc.).  

**2. Time series analysis on downloaded data**:

* Time series analysis of the following metrics: raw NDVI mean pixel intensity across the image, vegetation network centrality metric, and precipitation.
    * The time series are processed (outliers removed and resampled to avoid gaps). All time series of each sub-image are aggregated into one summary time series that is used for analysis. 
    * The summary time series is de-seasonalised and smoothed.
    * Residuals between the raw and de-seasonalised and smoothed time series are calculated and used for an early warning resilience analysis.

* Time series plots are produced, along with auto- and cross-correlation plots. Early warning signals are also computed using the [ewstools package](https://github.com/ThomasMBury/ewstools), including Lag-1 autocorrelation and standard deviation moving window plots. A sensitivity and significance analysis is also performed in order to determine whether any trends are statistically significant.
* Time series summary statistics and resilience metrics are saved into files.
* A PDF report is created showcasing the main figures resulting from the analyses. 

**Other functionalities**:

`pyveg` also has other minor functionalities:

* Analysis of a collection of summary data that has been created with the `pyveg` pipeline (downloading + time series analysis).
* Simulate the generation and evolution of patterned vegetation
* A stand-alone network centrality estimation for a 50x50 pixel image.
* A functionality to upload results to the Zenodo open source repository 

### `pyveg` flow

The diagram below represents the high level flow of the main functionalities of the `pyveg` package. For each main component there is a CLI console scripts defined, that is shown in the diagram. 

![The`pyveg` program flow.](paper/pveg_flow.png)

The full ReadTheDocs documentation for this `pyveg` can be found in this [link](https://pyveg.readthedocs.io/en/latest/).


This page contains an installation guide, and some usage examples for this package.


## Installation

`pyveg` requires Python 3.6 or greater. To install, start by creating a fresh `conda` environment.
```
conda create -n veg python=3.7
conda activate veg
```
Get the source.
```
git clone https://github.com/alan-turing-institute/monitoring-ecosystem-resilience.git
```
Enter the repository and check out a relevant branch if necessary (the default `master` branch contains the most up to date stable version of the code).
```
cd monitoring-ecosystem-resilience
```
Install the package using `pip`.
```
pip install .
```
If you are using Windows and encounter issues during this stage, a solution may be found here: https://github.com/NREL/OpenOA/issues/37. If you plan on making changes to the source code, you can instead run `pip install -e .`. 

Before using the Google Earth Engine API, you need to sign up with a Google account [here](https://earthengine.google.com/new_signup/), and authenticate.  To authenticate, run
```
earthengine authenticate
```
A new browser window will open. Copy the token from this window to the terminal prompt to complete the authentication process.


### Google Earth Engine

[Google Earth Engine](https://earthengine.google.com) (GEE) is a powerful tool for obtaining and analysing satellite imagery. This directory contains some useful scripts for interacting with the Earth Engine API. The earth engine API is installed automatically as part of the `pyveg` package installation. If you wish to install it separately, you can follow the instructions [here](https://developers.google.com/earth-engine/python_install_manual).

## Downloading data from GEE with ``pyveg``

### Downloading data from GEE using the CLI

To run a `pyveg` download job, use
```
pyveg_run_pipeline --config_file <path to config>
```

The download job is fully specified by a configuration file, which you point to using the `--config_file` argument. A sample config file is found at `pyveg/configs/config_all.py`. You can also optionally specify a string to identify the download job using the `--name` argument.

Note that we use the GEE convention for coordinates, i.e. `(longitude,latitude)`.

#### Generating a download configuration file 

To create a configuration file for use in the pyveg pipeline described above, use the command 
```
pyveg_generate_config
```
this allows the user to specify various characteristics of the data they want to download via prompts. The list in order is as follows:

* `--configs_dir`: The path to the directory containing the config file, with a default option `pyveg/configs`.

* `--collection_name`: The name of the dataset used in the collection, either Sentinel2, or Landsat 8, 7, 5 or 4.
    *    Sentinel2: [Available from 2015-06-23 at 10m resolution.](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2)
    *    Landsat8: [Available from 2013-04-11 at 30m resolution.](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C01_T1)
    *    Landsat7: [Available from 1999-01-01 at 30m resolution.](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C01_T1)
    *    Landsat5: [Available from 1984-03-10 to 2013-01-31 at 60m resolution.](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT05_C01_T1)
    *    Landsat4: [Available from 1982-07-16 to 1993-12-14 at 60m resolution.](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LT04_C01_T1)

* `--latitude`: The latitude (in degrees north) of the centre point of the image collection.

* `--longitude`: The longitude (in degrees east) of the centre point of the image collection.

* `--country`: The country (for the file name) can either be entered, or use the specified coordinates to look up the country name from the OpenCage database.

* `--start_date`: The start date in the format ‘YYYY-MM-DD’, the default is ‘2015-01-01’ (or ‘2019-01-01’ for a test config file).

* `--end_date`: The end date in the format ‘YYYY-MM-DD’, the default is today’s date (or ‘2019-03-01’ for a test config file).

* `--time_per_point`: The option to run the image collection either monthly (‘1m’) or weekly (‘1w’), with the default being monthly.

* `--run_mode`: The option to run time-consuming functions on Azure (‘batch’) or running locally on your own computer (‘local’). The default is local. For info about running on Azure go [here](UsingAzure.md).

* `--output_dir`: The option to write the output to a specified directory, with the default being the current directory.

* `--test_mode`: The option to make a test config file, containing fewer months and a subset of sub-images, with a default option to have a normal config file.
    *    By choosing the test config file, the start and end dates (see below) are defaulted to cover a smaller time span.
    *    It is recommended that the test config option should be used purely to determine if the options specified by the user are correct.


* `--n_threads`:  Finally, how many threads the user would like to use for the time-consuming processes, either 4 (default) or 8.

For example:
```
 pyveg_generate_config --configs_dir "pyveg/configs" --collection_name "Sentinel2" --latitude 11.58 --longitude 27.94 --start_date "2016-01-01" --end_date "2020-06-30" --time_per_point "1m" --run_mode "local" --n_threads 4
```

This generates a file named `config_Sentinel2_11.58N_27.94E_Sudan_2016-01-01_2020-06-30_1m_local.py` along with instructions on how to use this configuration file to download data through the pipeline, in this case the following:

```
pyveg_run_pipeline --config_file pyveg/configs/config_Sentinel2_11.58N_27.94E_Sudan_2016-01-01_2020-06-30_1m_local.py
```

Individual options can be specified by the user via prompt. The options for this can be found by typing ```pyveg_generate_config --help```. 


### More Details on Downloading

During the download job, `pyveg` will break up your specified date range into a time series, and download data at each point in the series. Note that by default the vegetation images downloaded from GEE will be split up into 50x50 pixel images, vegetation metrics are then calculated on the sub-image level. Both colour (RGB) and Normalised Difference Vegetation Index (NDVI) images are downloaded and stored. Vegetation metrics include the mean NDVI pixel intensity across sub-images, and also network centrality metrics, discussed in more detail below.

For weather collections e.g. the ERA5, due to coarser resolution, the precipitation and temperature "images" are averaged into a single value at each point in the time series.

### Rerunning partially succeeded jobs

The output location of a download job is datestamped with the time that the job was launched.  The configuration file used will also be copied and datestamped, to aid reproducibility.  For example if you run the job
```
pyveg_run_pipeline --config_file pyveg/configs/my_config.py
```
there will be a copy of `my_config.py` saved as `pyveg/configs/cached_configs/my_config_<datestamp>.py`. This also means that if a job crashes or timeouts partway through, it is possible to rerun, writing to the same output location and skipping parts that are already done by using this cached config file.  However, in order to avoid a second datestamp being appended to the output location, use the ```--from_cache``` argument.  So for the above example, the command to rerun the job filling in any failed/incomplete parts would be:
```
pyveg_run_pipeline --config_file pyveg/configs/cached_configs/my_config_<datestamp>.py --from_cache
```

### Using Azure for downloading/processing data

If you have access to Microsoft Azure cloud computing facilities, downloading and processing data can be sped up enormously by using batch computing to run many subjobs in parallel.  See [here](UsingAzure.md) for more details.

### Downloading data using the API

Although `pyveg` has been mostly designed to be used with the CLI as shown above, we can also use `pyveg` functions through the API. A tutorial of how to download data this way is included in the ``notebooks/tutorial_download_and_process_gee_images.ipynb`` notebook tutorial. 

## Analysing the Downloaded Data with `pyveg`

### Analysing the Downloaded Data using the CLI

Once you have downloaded the data from GEE, the `pyveg_gee_analysis` command allows you to process and analyse the output. To run:
```
pyveg_gee_analysis --input_dir <path_to_pyveg_download_output_dir>
```
The analysis code preprocesses the data and produces a number of plots. These will be saved in an `analysis/` subdirectory inside the `<path_to_pyveg_download_output_dir>` directory.

Note that in order to have a meaningful analysis, the dowloaded time series should have at least 4 points (and more thant 12
for an early warning analysis) and not being the result of a "test" config file, in this case the analysis fails.
 
The commands also allows for other options to be added to the execution of the script (e.g. run analysis from a downloaded data in Azure blob storage, define a different output directly, don't include a time series analysis, etc), which can be displayed by typing:

```
pyveg_gee_analysis --help
```

The analysis script executed with the ```pyveg_gee_analysis``` command runs the following steps:

#### Preprocessing

`pyevg` supports the following pre-processing operations:
- Identify and remove outliers from the time series.
- Fill missing values in the time series (based on a seasonal average), or resample the time series using linear interpolation between points.
- Smoothing of the time series using a LOESS smoother.
- Calculation of residuals between the raw and smoothed time series.
- De-seasonalising (using first differencing), and detrending using STL.

#### Plots

In the `analysis/` subdirectory, `pyveg` creates the following plots:
- Time series plots containing vegetation and precipitation time series (seasonal and de-seasonalised). Plots are labelled with the AR1 of the vegetation time series, and the maximum correlation between the Vegetation and precipitation time series.
- Auto-correlation plots for vegetation and precipitation time series (seasonal and de-seasonalised).
- Vegetation and precipitation cross-correlation scatterplot matrices.
- STL decomposition plots.
- Resilience analysis:
     - `ewstools` resilience plots showing AR1, standard deviation, skewness, and kurtosis using a moving window.
     - Smoothing filter size and moving window size Kendall tau sensitivity plots.
     - Significance test.

### Running the analysis using the API

Although `pyveg` has been mostly designed to be used with the CLI as shown above, we can also use `pyveg` functions through the API. A tutorial of how to run the data analysis in this way is included in the ```notebooks/tutorial_analyse_gee_data.ipynb``` notebook  tutorial. 


## Other functionalities of `pyveg`

### Analysis summary statistics data 

The ```analyse_pyveg_summary_data.py``` functionality processes collections of data produced by the main `pyveg` pipeline described in the section above (download + time series analysis for different locations and time periods) and creates a number of plots of the summary statistics of these time series.

To run this analysis in Python, there is an entrypoint defined.  Type:

```
pyveg_analysis_summary_data --input_location  <path_to_directory_with_collection_summary_statistics>
```

if you wish you can also specify the ```outpur_dir``` where plots will be saved.  Type:
```
pyveg_analysis_summary_data --help
```
to see the extra options.


### Pattern simulation

The ```generate_patterns.py``` functionality originates from some Matlab code by Stefan Dekker, Willem Bouten, Maarten Boerlijst and Max Rietkerk (included in the "matlab" directory), implementing the scheme described in:

Rietkerk et al. 2002. Self-organization of vegetation in arid ecosystems. The American Naturalist 160(4): 524-530.

To run this simulation in Python, there is an entrypoint defined.  Type:
```
pyveg_gen_pattern --help
```
to see the options.  The most useful option is the `--rainfall` parameter which sets a parameter (the rainfall in mm) of the simulation - values between 1.2 and 1.5 seem to give rise to a good range of patterns. Other optional parameters for `generate_patterns.py` allow the generated image to be output as a csv file or a png image.  The `--transpose` option rotates the image 90 degrees (this was useful for comparing the Python and Matlab network-modelling code). Other parameters for running the simulation are in the file `patter_gen_config.py`, you are free to change them.

#### Running the pattern simulation using the API

Although `pyveg` has been mostly designed to be used with the CLI as shown above, we can also use `pyveg` functions through the API. A tutorial of how to run the simulation of the pattern generation in this way is included in [here](notebooks/tutorial_simulate_patterned_vegetation.ipynb). 


### Network centrality

There is an entrypoint defined in `setup.py` that runs the *main* function of `calc_euler_characteristic.py`:
```
pyveg_calc_EC --help
```
will show the options.

* `--input_txt` allows you to give the input image as a csv, with one row per row of pixels.  Inputs are expected to be "binary", only containing two possible pixel values (typically 0 for black and 255 for white).
* `--input_img` allows you to pass an input image (png or tif work OK).  Note again that input images are expected to be "binary", i.e. only have two colours.
* `--sig_threshold` (default value 255) is the value above (or below) which a pixel is counted as signal (or background)
* `--upper_threshold` determines whether the threshold above is an upper or lower threshold (default is to have a lower threshold - pixels are counted as "signal" if their value is greater-than-or-equal-to the threshold value).
* `--use_diagonal_neighbours` when calculating the adjacency matrix, the default is to use "4-neighbours" (i.e. pixels immediately above, below, left, or right).  Setting this option will lead to "8-neighbours" (i.e. the four neighbours plus those diagonally adjacent) to be included.
* `--num_quantiles` determines how many elements the output feature vector will have.
* `--do_EC` Calculate the Euler Characteristic to fill the feature vector.  Currently this is required, as the alternative approach (looking at the number of connected components) is not fully debugged.

Note that if you use the `-i` flag when running python, you will end up in an interactive python session, and have access to the `feature_vec`, `sel_pixels` and `sc_images` variables.

Examples:
```
pyveg_calc_EC --input_txt ../binary_image.txt --do_EC
>>> sc_images[50].show() # plot the image with the top 50% of pixels (ordered by subgraph centrality) highlighted.
>>> plt.plot(list(sel_pixels.keys()), feature_vec, "bo") # plot the feature vector vs pixel rank
>>> plt.show()
```

### Uploading results to the Zenodo open source repository

See [here](UsingZenodo.md) for more details.

# Contributing 

We welcome contributions from anyone who is interested in the project. There are lots of ways to contribute, not just writing code. See our [Contributor Guidelines](CONTRIBUTING.md) to learn  more about how you can contribute and how we work together as a community.

# Licence

This project is licensed under the terms of the MIT software license.
# Contributing to monitoring-ecosystem-resilience (the repo!)

**Welcome to the repository!**
We're excited you're here and want to contribute.

We hope that these guidelines make it as easy as possible to get involved.
If you have any questions that aren't discussed below, please let us know by opening an [issue](https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/issues).

We welcome all contributions from documentation to testing to writing code.
Don't let trying to be perfect get in the way of being good - exciting ideas are more important than perfect pull requests.

## Table of contents

- [Where to start: issues](#where-to-start-issues)
- [Making a change with a pull request](#making-a-change-with-a-pull-request)
  - [1. Comment on an existing issue or open a new issue referencing your addition](#1-comment-on-an-existing-issue-or-open-a-new-issue-referencing-your-addition)
  - [2. Create a new branch (if you have *write* access to the repository) or fork the repository to your profile (if you don't currently have _write_ access)](#2-create-a-new-branch-or-fork-the-repository-to-your-profile)
  - [3. Make the changes you've discussed](#3-make-the-changes-youve-discussed)
  - [4. Submit a pull request](#4-submit-a-pull-request)
- [Style guide](#style-guide)

## Where to start: issues

* **Issues** are individual pieces of work that need to be completed to move the project forwards.
A general guideline: if you find yourself tempted to write a great big issue that
is difficult to describe as one unit of work, please consider splitting it into two or more issues.

Before you open a new issue, please check if any of our [open issues](https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/issues) covers your idea already.

The list of labels for current issues includes:

- [![help-wanted](https://img.shields.io/badge/-help%20wanted-159818.svg)][labels-helpwanted] _These issues contain a task that a member of the team has determined we need additional help with._

  If you feel that you can contribute to one of these issues, we especially encourage you to do so!

- [![question](https://img.shields.io/badge/-question-cc317c.svg)][labels-question] _These issues contain a question that you'd like to have answered._

  Opening an issue is a great way to start a conversation and get your answer.

- [![good-first-issue](https://img.shields.io/badge/-good%20first%20issue-1b3487.svg)][labels-firstissue] _These issues are particularly appropriate if it is your first contribution to the repository, or to GitHub overall._

- [![Enhancement](https://img.shields.io/badge/-enhancement-84b6eb.svg)][labels-enhancement] _These issues are suggesting new features that can be added to the project._

  If you want to ask for something new, please try to make sure that your request is distinct from any others that are already in the queue.
  If you find one that's similar but there are subtle differences please reference the other enhancement in your issue.

- [![Bug](https://img.shields.io/badge/-bug-d73a4a.svg)][labels-bug] _These issues are reporting a problem or a mistake in the project._

  The more details you can provide the better!
  If you know how to fix the bug, please open an issue first and then submit a pull request.

- [![project-management](https://img.shields.io/badge/-project%20management-bfd86c.svg)][labels-project-management] _We like to model best practice, so the package itself is managed through these issues.

## Making a change with a pull request

We appreciate all contributions to monitoring-ecosystem-resilience.
**THANK YOU** for helping us.

All project management, conversations and questions related to the project happens here in the [monitoring-ecosystem-resilience repository][monitoring-ecosystem-resilience-repo].

In brief, the structure for making a contribution is as follows:
1. Identify a specific change that needs to be made to the repository. Open a new issue (after checking one does not already exist!) and describe the change, include why you are making it.
2. Create a new branch corresponding to this issue. The new branch will house all the changes that you make to the repository in an isolated location. As discussed in more detail below, new branches should be created using the latest version of the `develop` branch.
3. Make commits to the new branch you have created.
4. Submit a pull request to add the modifications in your new branch back into `develop`.

When a significant milestone has been reached, and the `develop` branch is known to be in a stable configuration, the `master` branch will be updated via a pull request from `develop`. In general, commits should not be made to either the `master` or `develop` branches. Pull requests to `develop` are fine (and encoraged), while pull requests to `master` will happen in a coordinated way.

The following steps are a more detailed guide to help you contribute in a way that will be easy for everyone to review and accept with ease.

### 1. Comment on an [existing issue](https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/issues) or open a new issue referencing your addition

This allows other members of the team to confirm that you aren't overlapping with work that's currently underway and that everyone is on the same page with the goal of the work you're going to carry out.

[This blog](https://www.igvita.com/2011/12/19/dont-push-your-pull-requests/) is a nice explanation of why putting this work in up front is so useful to everyone involved.

### 2. Create a new [branch][github-branches] or [Fork][github-fork] the [monitoring-ecosystem-resilience repository][monitoring-ecosystem-resilience-repo] to your profile

#### 2a) Create a branch
If you are a collaborator on the repository with write access, then you can make a [new branch][github-branches].  We recommend that you start from the latest version of the `develop` branch, and create a new one from there. This is the branch we use for active deleopment of the repository, while stable (but not cutting edge) versions are in the `master` branch. The name of your new branch should ideally be in the format: `<feature|bugfix>/<issue-number>-<short-description>`. For example, if you were addressing Issue number 111 which was about incorrect JSON filenames, it could be something like:
```
git checkout develop
git pull
git checkout -b bugfix/111-fix-json-filenames
```
Now you can go to step #3, where you actually fix the problem! :)

In case you want to learn more about "branching out", [this blog](https://nvie.com/posts/a-successful-git-branching-model/) details the different Git branching models.


#### 2b. Fork the repository

If you don't have write access to the repository, you can fork it to your own profile.
This is now your own unique copy of the repo.
Changes here won't affect anyone else's work, so it's a safe space to explore edits to the code!

Make sure to [keep your fork up to date][github-syncfork] with the master repository, otherwise you can end up with lots of dreaded [merge conflicts][github-mergeconflicts].

### 3. Make the changes you've discussed

Try to keep the changes focused.
If you submit a large amount of work all in one go it will be much more work for whomever is reviewing your pull request.

While making your changes, commit often and write good, detailed commit messages.
[This blog](https://chris.beams.io/posts/git-commit/) explains how to write a good Git commit message and why it matters.
It is also perfectly fine to have a lot of commits - including ones that break code.
A good rule of thumb is to push up to GitHub when you _do_ have passing tests then the continuous integration (CI) has a good chance of passing everything.

Please do not re-write history!
That is, please do not use the [rebase](https://help.github.com/en/articles/about-git-rebase) command to edit previous commit messages, combine multiple commits into one, or delete or revert commits that are no longer necessary.

### 4. Submit a [pull request][github-pullrequest]

A "pull request" is a request to "pull" the changes you have made in your branch back into another branch of the repository. The source branch will be the new branch you created in order to address the issue you created/choose. The destination branch should generally be `develop`, where all main code development takes place. Avoid making pull requests into the `master` branch (pull requests into master should happen in a coordinated way using a stable configuration of `develop` as the source branch).

We encourage you to open a pull request as early in your contributing process as possible.
This allows everyone to see what is currently being worked on.
It also provides you, the contributor, feedback in real time from both the community and the continuous integration as you make commits (which will help prevent stuff from breaking).

When you are ready to submit a pull request, make sure the contents of the pull request body do the following:
- Describe the problem you're trying to fix in the pull request, reference any related issues and use keywords fixes/close to automatically close them, if pertinent.
- List changes proposed in the pull request.
- Describe what the reviewer should concentrate their feedback on.

If you have opened the pull request early and know that its contents are not ready for review or to be merged, add "[WIP]" at the start of the pull request title, which stands for "Work in Progress".
When you are happy with it and are happy for it to be merged into the main repository, change the "[WIP]" in the title of the pull request to "[Ready for review]".

A member of the team will then review your changes to confirm that they can be merged into the main repository.
A [review][github-review] will probably consist of a few questions to help clarify the work you've done.
Keep an eye on your GitHub notifications and be prepared to join in that conversation.

You can update your [fork][github-fork] of the [repository][monitoring-ecosystem-resilience-repo] and the pull request will automatically update with those changes.
You don't need to submit a new pull request when you make a change in response to a review.

You can also submit pull requests to other contributors' branches!
Do you see an [open pull request](https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/pulls) that you find interesting and want to contribute to?
Simply make your edits on their files and open a pull request to their branch!

What happens if the continuous integration (CI) fails (for example, if the pull request notifies you that "Some checks were not successful")?
The CI could fail for a number of reasons.
At the bottom of the pull request, where it says whether your build passed or failed, you can click “Details” next to the test, which takes you to the Travis page.
You can view the log or rerun the checks if you have write access to the repo by clicking the “Restart build” button in the top right (you must be logged in to Travis CI with your GitHub account see the “Restart build” button).

GitHub has a [nice introduction][github-flow] to the pull request workflow, but please get in touch if you have any questions.

## Style Guide

Docstrings should follow [numpydoc][link_numpydoc] convention.
We encourage extensive documentation.

The python code itself should follow [PEP8][link_pep8] convention whenever possible, with at most about 500 lines of code (not including docstrings) per script.

---

_These Contributing Guidelines have been adapted from the [Contributing Guidelines](https://github.com/bids-standard/bids-starter-kit/blob/master/CONTRIBUTING.md) of [The Turing Way](https://github.com/alan-turing-institute/the-turing-way)! (License: MIT)_

[monitoring-ecosystem-resilience-repo]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/
[monitoring-ecosystem-resilience-issues]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/issues
[git]: https://git-scm.com
[github]: https://github.com
[github-branches]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository
[github-fork]: https://help.github.com/articles/fork-a-repo
[github-flow]: https://guides.github.com/introduction/flow
[github-mergeconflicts]: https://help.github.com/articles/about-merge-conflicts
[github-pullrequest]: https://help.github.com/articles/creating-a-pull-request
[github-review]: https://help.github.com/articles/about-pull-request-reviews
[github-syncfork]: https://help.github.com/articles/syncing-a-fork
[labels-bug]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/bug
[labels-enhancement]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/enhancement
[labels-firstissue]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/good%20first%20issue
[labels-helpwanted]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/help%20wanted
[labels-project-management]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/project%20management
[labels-question]: https://github.com/alan-turing-institute/monitoring-ecosystem-resilience/labels/question
[link_numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[link_pep8]: https://www.python.org/dev/peps/pep-0008/
## R code for generating and analyzing vegetation patterns.

The code in this directory comprises the *rveg* package, which contains functions for generating vegetation patterns (optionally evolving them from an input starting pattern), and performing a network centrality analysis on them.

To use, fro this directory, load the package:
```
devtools::load_all()
```

### Generating a pattern:
```
pattern <-rveg::generatePattern()
```

### Calculate Euler Characteristic for a pattern

```
featureVec <- calc_EC(pattern)
```
where "pattern" is a 2D array of 1s and 0s (with 1 representing vegetation and 0 representing bare soil).


### Running tests

From this directory, do
```
Rscript -e "devtools::test()"
```## Animation of vegetation patterns

This is a translation of some `matlab` code to visualize a simple model of vegetation patterns in arid landscapes.
To run:
From *R* or *RStudio*
```
source('patterns.R')
animate
```

Configuration parameters are loaded from the file ```config.json```.


Dependencies:
```
dplyr
ggplot2
gganimate
jsonlite
```
These can be installed from CRAN via ```install.packages(<package_name>)```.
