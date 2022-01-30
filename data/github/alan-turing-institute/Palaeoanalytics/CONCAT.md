<h1 align="center">Welcome to Palaeoanalytics! </h1>

> Repository for the [Palaeoanalytics project](https://www.turing.ac.uk/research/research-projects/palaeoanalytics). 
> A collaboration between The Alan Turing Institute and the University of Cambridge.  

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://app.travis-ci.com/alan-turing-institute/Palaeoanalytics.svg?token=sMJzQpXKRs31ujsqXNxP&branch=develop)](https://app.travis-ci.com/alan-turing-institute/Palaeoanalytics)

# **Table of Contents:**

- [About the project](#about)
- [The team](#team)
- [The `PyLithics` package](#pylithics)
- [Drawing style for `PyLithics`](#drawing)
- [Contributing](#contributing)
- [Development and testing](#testing)
- [Citing `PyLithics`](#cite)
- [Licence](#licence)

# 📖 About the project <a name="about"></a>

Archaeologists have long used stone tools (lithics) to reconstruct the behavior of prehistoric hominins. While techniques 
have become more quantitative, there still remain barriers to optimizing data retrieval. Machine learning and computer 
vision approaches can be developed to extract quantitative and trait data from lithics, photographs and drawings. `PyLithics`
has been developed to capture data from 2D line drawings, focusing on the size, shape and technological attributes of flakes. 

`PyLithics`is an open-source, free for use software package for processing lithic artefact illustrations scanned from 
the literature. This tool accurately identifies, outlines, and computes lithic shape and linear measures, and returns user 
ready data. It has been optimized for feature extraction and measurement using a number of computer vision techniques 
including pixel intensity thresholding, edge detection, contour finding, custom template matching and image kernels. 
On both conventional and modern drawings, `PyLithics` can identify and label platform, lateral, dorsal, and ventral surfaces,
as well as individual dorsal surface scar shape, size, orientation, diversity, number, and flaking order. Complete size
and shape metrics of individual scars and whole flakes can be calculated and recorded. Orientation and flaking direction 
of dorsal scars can also be calculated. The resulting data can be used for metrical analysis, extracting features indicative
of typologies and technological processes. Data output can easily be employed to explore patterns of variation within and between assemblages.

# 👥 The team <a name="team"></a>

These are the members of the Palaeoanalytics team as updated August 2021:

| Name | Role | email | Github | 
| --- | --- | --- | --- |
| Jason Gellis | Postdoctoral Researcher (University of Cambridge) | [jg760@cam.ac.uk](mailto:jg760@cam.ac.uk) | [@JasonGellis](https://github.com/JasonGellis) |
| Camila Rangel Smith | Research Data Scientist (The Alan Turing Institute) | [crangelsmith@turing.ac.uk](mailto:crangelsmith@turing.ac.uk) |[@crangelsmith](https://github.com/crangelsmith) |
| Robert Foley | Principal Investigator (University of Cambridge) | [raf10@cam.ac.uk](mailto:raf10@cam.ac.uk)| [Rob-LCHES](https://github.com/Rob-LCHES)

# 📦 The `PyLithics` package <a name="pylithics"></a>

> PyLithics: A Python package for stone tool analysis

## Workflow

`PyLithics` is devised to work with illustrations of lithic objects common to publications in archaeology and anthropology. Lithic illustrators have established conventions regarding systems of artefact orientation and proportions. Lithics are normally drawn at a 1:1 scale, with the vertical axis orthogonal to the striking platform. A preferred method is to orient and illustrate various aspects of an artefact as a series of adjacent surfaces at 90-degree rotations from the principal view (usually the dorsal surface). Each aspect contains internal details (i.e., flake scars, cortical areas, etc.), indication of flaking direction radial lines (ripples), and the inclusion of a metric scale (for more information about lithic drawings see [@Martingell:1988]). Currently, `PyLithics` is optimised to work with unifacial flakes and bifaces, which are relatively flat, two-dimensional objects. 

The inputs for `PyLithics` are images of lithic objects, images of their associated scales, and a metadata `CSV` file linking the two and giving the scale measurement in millimeters. 

`PyLithics` processes the images with the following steps (and as illustrated in the schema below):

1. Import and match images to associated image ID and scale image from CSV metadata file.
2. Calculate a conversion of pixels to millimeters based on the size of the associated scale from CSV metadata file. If no scale is present, measurements will be in pixels
3. Apply noise removal and contrast stretching to images to minimise pixel variation.
4. Pixel intensity thresholding of images to prepare for contour finding.
5. Apply edge detection and contour finding to thresholded images.
6. Calculate metrics of lithic surface features from found contours -- area, length, breath, shape, number of vertices. 
7. Select contours which outline an entire lithic object's surfaces, or select contours of inner scars greater than 3% and less than 50% of the total size of its surface.
8. Classify these selected surface contours as "Dorsal", "Ventral", "Lateral", and/or "Platform" depending on presence or absence. Assign scar contours to these surfaces. 
9. If present, find arrows using connected components and template matching, measure their angle and assign angle to associated scar.
10. Plot resulting surface and scar contours on the original images for validation.
11. Output data in a hierarchical json file detailing measurements of surface and scar contours. 

Here you can find a schema of the workflow described above:

<img src="figures/pylithics_flowchart.jpg"/>

## Installation

The `PyLithics` package requires Python 3.7 or greater. To install, start by creating a fresh virtual environment.
```
python3 -m venv palaeo
source palaeo/bin/activate
```
For Windows OS:
```
Set-ExecutionPolicy Unrestricted -Scope Process
.\palaeo\Scripts\activate
```
Clone the repository.
```
git clone https://github.com/alan-turing-institute/Palaeoanalytics.git
```
Enter the repository and check out a relevant branch if necessary (the `develop` branch contains the most up-to-date stable version of the code, but this branch is fast moving.
If you want to have a stable and static version it is better to use `main` branch).
```
cd Palaeoanalytics
git checkout main
```
Install 'PyLithics'.
```
pip install .
```
The `pip install .` command will call `setup.py` to install and configure PyLithics and its required packages listed in the [requirements.txt](requirements.txt) file. 

**Note**: For Mac users we recommend an OS versions=> 10.14 to prevent build problems. 

## Running `PyLithics`

`PyLithics` can be run via command line. The following command displays all available options:
```bash
pylithics_run --help
```
Output:
```bash
usage: pylithics_run [-h] -c config-file [--input_dir INPUT_DIR]
                     [--output_dir OUTPUT_DIR]

Run lithics characterisation pipeline

optional arguments:
  -h, --help            show this help message and exit
  -c config-file, --config config-file
                        the model config file (YAML)
  --input_dir INPUT_DIR
                        path to input directory where images are found
  --output_dir OUTPUT_DIR
                        path to output directory to save processed image
                        outputs
  --metadata_filename METADATA_FILENAME
                        CSV file with metadata on images and scales
  --get_arrows          If a lithic contains arrows, find them and add them to
                        the data

```
## 💫 Quickstart 

**In order to provide a quick start we have provided an [example dataset](data) including images, scales and metadata.** You
can run a quick analysis in this dataset by running:
```python
pylithics_run -c configs/test_config.yml --input_dir data --output_dir output --metadata_filename meta_data.csv --get_arrows
```
More generally, given that you have a set of lithics images (and its respective scales), you can run the `PyLithics` processing script with the following:

```python
pylithics_run -c configs/test_config.yml --input_dir <path_to_input_dir> --output_dir <path_to_output_directory> --metadata_filename metatada_file.csv
```
The images found in ```<path_to_input_dir>``` should follow this directory structure:

```bash
input_directory
   ├── metatada_file.csv
   ├── images 
        ├── lithic_id1.png
        ├── lithic_id2.png
        └── lithic_id3.png
            .
            .
            .
        ├── lithic_idn.png
   └──  scales
        ├── scale_id1.png
        ├── scale_id2.png
        ├── scale_id3.png
            .
            .
            .
        └── scale_id4.png



```

where the mapping between the lithics and scale images should be available in the metadata CSV file. 

This CSV file should have as a minimum the following 3 variables:
 
- *PA_ID*: corresponding the lithics image id
(the name of the image file), 
- *scale_ID*: The scale id (name of the scale image file)
- *PA_scale*: The scale measurement (how many centimeters this scale represents).

An example of this table, where one scale corresponds to several images is the following:

| PA_ID      | scale_ID  | PA_scale | 
|------------|-----------|----------|
| lithic_id1 | scale_id1 | 5        | 
| lithic_id2 | scale_id2 | 5        |
| lithic_id3 | scale_id3 | 5        |   

**Note**

In the scenario that the scale and csv file are not available, it is possible to run the analysis only using the images
with the command:
```
pylithics_run -c configs/test_config.yml --input_dir <path_to_input_dir> --output_dir <path_to_output_directory> 
```
lithics image files must still be inside the '<path_to_input_dir>/images/' directory. However, all the measurements will only be
provided as number of pixels. 

The ```test_config.yml``` config file contains the following options:

```yaml

threshold: 0.01
contour_parameter: 0.1
contour_fully_connected: 'low'
minimum_pixels_contour: 0.01
denoise_weight: 0.06
contrast_stretch: [4, 96]

```

The config is optimised to work with the images in an [example dataset](data). If you want to use `PyLithics` with different styles of
drawing you might have to modify this configuration file. You can modify or create your on config file and provide it to the CLI. 

## Output from `PyLithics`

### Output images

Output images are saved in the output directory for validation of the data extraction process. An example of these images
are the following: 

<p float="left">
<img src="figures/rub_al_khali_lithic_surfaces.png" width="310" />
<img src="figures/rub_al_khali_lithium_scars.png" width="310" />
<img src="figures/rub_al_khali_lithium_angles.png" width="310" />
</p>

### Output data

The output dataset is a JSON file with data for the lithic objects found in an image. The data is 
hierarchically organised by type of surface object (ventral, dorsal, platform). For each 
surface the metrics from its scars are recorded. In [this data output example,](output_example.md) you can find the json file
that results from running `PyLithics` on the above images, with comments to better understand the feature hierarchy and variables. 

# 🖌 Drawing style for `PyLithics` <a name="drawing"></a>

We are working hard in developing methods to cater to all styles of stone tools drawings. However, at the moment `PyLithics`
works best with the following styles:

<img src="figures/drawing_style.png"/>

If you want to help us optimise `PyLithics` for different drawing styles we welcome your [contributions](#contributing)!

# 👋 Contributing <a name="contributing"></a>

We welcome contributions from anyone interested in the project. There are lots of ways to contribute, not just writing code.
If you have ideas on how to extend/improve `PyLithics` do get in touch with members of the team via email. See our 
[Contributor Guidelines](CONTRIBUTING.md) to learn more about how you can contribute and how we work together as a
community in GitHub. Because `PyLithics'` code changes frequently we test and deploy current builds and updates via [Travis CI](https://docs.travis-ci.com/user/for-beginners/).
Every time a change in the `PyLitihcs` code is pushed to the Palaeoanalytics repository, the [`travis.yml`](.travis.yml) 
file, which contains essential information about the `PyLithics` programming environment and version, triggers these automated
tests. TravisCI will automatically create a virtual build of `PyLithics`, and run the software to ensure that integration
of new code is stable and functioning. Upon completion of tests, TravisCI will generate a virtual build *pass* or *fail*
report and notify `PyLithics` team members and contributing developers of any issues. Because the process is automated 
there is no need for contributors to open a TravisCI account. All contributions will have to successfully pass these automated 
tests to be merged into the `main` branch.

# Development and testing of `PyLithics` <a name="testing"></a>

`PyLithics` uses the [pytest](https://pypi.org/project/pytest/) library for automated functional testing of code 
development and integration. These [tests](tests) are easily run from the project directory using the command:

`pytest -s `

# Citing `PyLithics` <a name="cite"></a>

`PyLithics` publication in the Journal of Open Source Software:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03738/status.svg)](https://doi.org/10.21105/joss.03738)

`PyLithics` software version history on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5898149.svg)](https://doi.org/10.5281/zenodo.5898149)

# 📝 Licence <a name="licence"></a>

This software is licensed under the terms of the [GNU General Public License v3.0 (GNU GPLv3)](https://choosealicense.com/licenses/gpl-3.0/).
# Contributing to Palaoanalytics (the repo!)

**Welcome to the repository!**
We're excited you're here and want to contribute.

We hope that these guidelines make it as easy as possible to get involved.
If you have any questions that aren't discussed below, please let us know by opening an [issue](https://github.com/alan-turing-institute/Palaeoanalytics/issues).

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

Before you open a new issue, please check if any of our [open issues](https://github.com/alan-turing-institute/Palaeoanalytics/issues) covers your idea already.

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

We appreciate all contributions to the Palaoanalytics project and the PyLythics package.
**THANK YOU** for helping us.

All project management, conversations and questions related to the project happens here in the [Palaoanalytics][Palaoanalytics-repo].

In brief, the structure for making a contribution is as follows:
1. Identify a specific change that needs to be made to the repository. Open a new issue (after checking one does not already exist!) and describe the change, include why you are making it.
2. Create a new branch corresponding to this issue. The new branch will house all the changes that you make to the repository in an isolated location. As discussed in more detail below, new branches should be created using the latest version of the `develop` branch.
3. Make commits to the new branch you have created.
4. Submit a pull request to add the modifications in your new branch back into `develop`.

When a significant milestone has been reached, and the `develop` branch is known to be in a stable configuration, the `master` branch will be updated via a pull request from `develop`. In general, commits should not be made to either the `master` or `develop` branches. Pull requests to `develop` are fine (and encoraged), while pull requests to `master` will happen in a coordinated way.

The following steps are a more detailed guide to help you contribute in a way that will be easy for everyone to review and accept with ease.

### 1. Comment on an [existing issue](https://github.com/alan-turing-institute/Palaeoanalytics/issues) or open a new issue referencing your addition

This allows other members of the team to confirm that you aren't overlapping with work that's currently underway and that everyone is on the same page with the goal of the work you're going to carry out.

[This blog](https://www.igvita.com/2011/12/19/dont-push-your-pull-requests/) is a nice explanation of why putting this work in up front is so useful to everyone involved.

### 2. Create a new [branch][github-branches] or [Fork][github-fork] the [Palaoanalytics repository][Palaeoanalytics-repo] to your profile

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

You can update your [fork][github-fork] of the [repository][Palaonalytics-repo] and the pull request will automatically update with those changes.
You don't need to submit a new pull request when you make a change in response to a review.

You can also submit pull requests to other contributors' branches!
Do you see an [open pull request](https://github.com/alan-turing-institute/Palaoanalytics/pulls) that you find interesting and want to contribute to?
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

[Palaoanalytics-repo]: https://github.com/alan-turing-institute/Palaeoanalytics
[Palaoanalytics-issues]: https://github.com/alan-turing-institute/Palaeoanalytics/issues
[git]: https://git-scm.com
[github]: https://github.com
[github-branches]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository
[github-fork]: https://help.github.com/articles/fork-a-repo
[github-flow]: https://guides.github.com/introduction/flow
[github-mergeconflicts]: https://help.github.com/articles/about-merge-conflicts
[github-pullrequest]: https://help.github.com/articles/creating-a-pull-request
[github-review]: https://help.github.com/articles/about-pull-request-reviews
[github-syncfork]: https://help.github.com/articles/syncing-a-fork
[labels-bug]: https://github.com/alan-turing-institute/Palaeoanalytics/labels/bug
[labels-enhancement]: https://github.com/alan-turing-institute/Palaeoanalytics/labels/enhancement
[labels-firstissue]: https://github.com/alan-turing-institute/Palaeoanalytics/labels/good%20first%20issue
[labels-helpwanted]: https://github.com/alan-turing-institute/Palaeoanalytics/labels/help%20wanted
[labels-project-management]: https://github.com/alan-turing-institute/Palaeoanalytics/labels/project%20management
[labels-question]: https://github.com/alan-turing-institute/Palaeoanalytics/labels/question
[link_numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[link_pep8]: https://www.python.org/dev/peps/pep-0008/This is an example of the output data with comments to 
understand the variables:

```json
{
   "id":"rub_al_khali", // name of the image
   "conversion_px":0.040, // conversion from pixel to mm
   "n_surfaces":4, // number of outer surfaces found
   "lithic_contours":[
      {
         "surface_id":0, // largest surface id
         "classification":"Ventral", // surface classification
         "total_area_px":515662.0, // total area of surface in pixels
         "total_area":808.2, // total area of surface in mm
         "max_breadth":22.0, // surface maximum breadth
         "max_length":53.6, // surface maximum lengh
         "polygon_count":7, // numer of vertices measured in an approximate polygon fitted to the surface
         "scar_count":0, // number of scars in that surface
         "percentage_detected_scars":0.0, // percentage of the surface that contains scars
         "scar_contours":[ // empty scar count
         ]
      },
      {
         "surface_id":1, // second largest surface id
         "classification":"Dorsal",
         "total_area_px":515583.0,
         "total_area":808.0,
         "max_breadth":22.0,
         "max_length":53.6,
         "polygon_count":7,
         "scar_count":5,
         "percentage_detected_scars":0.71,
         "scar_contours":[
            {
               "scar_id":0, // largest scar belonging to surface id = 1
               "total_area_px":139998.0, // total area in pixels of scar
               "total_area":219.4, // total area in mm of scar
               "max_breadth":10.6, // scar maximum breadth
               "max_length":42.1, // scar maximum lenght
               "percentage_of_surface":0.27, // percentage of the scar to the total surface
               "scar_angle":1.74, // angle measured of arrow belonging to that scar
               "polygon_count":5 // numer of vertices measured in an approximate polygon fitted to the scar
            },
            {
               "scar_id":1,
               "total_area_px":111052.5,
               "total_area":174.0,
               "max_breadth":7.6,
               "max_length":43.5,
               "percentage_of_surface":0.22,
               "scar_angle":356.78,
               "polygon_count":6
            },
            {
               "scar_id":2,
               "total_area_px":103554.0,
               "total_area":162.3,
               "max_breadth":6.8,
               "max_length":42.4,
               "percentage_of_surface":0.2,
               "scar_angle":5.6,
               "polygon_count":4
            },
            {
               "scar_id":3,
               "total_area_px":6288.0,
               "total_area":9.9,
               "max_breadth":4.4,
               "max_length":5.9,
               "percentage_of_surface":0.01,
               "scar_angle":"NaN",
               "polygon_count":7
            },
            {
               "scar_id":4,
               "total_area_px":5853.0,
               "total_area":9.2,
               "max_breadth":3.9,
               "max_length":3.4,
               "percentage_of_surface":0.01,
               "scar_angle":"NaN",
               "polygon_count":6
            }
         ]
      },
      {
         "surface_id":2,
         "classification":"Lateral",
         "total_area_px":162660.5,
         "total_area":254.9,
         "max_breadth":8.2,
         "max_length":53.8,
         "polygon_count":3,
         "scar_count":2,
         "percentage_detected_scars":0.47,
         "scar_contours":[
            {
               "scar_id":0,
               "total_area_px":57245.5,
               "total_area":89.7,
               "max_breadth":5.4,
               "max_length":51.5,
               "percentage_of_surface":0.35,
               "scar_angle":"NaN",
               "polygon_count":3
            },
            {
               "scar_id":1,
               "total_area_px":18672.5,
               "total_area":29.3,
               "max_breadth":1.9,
               "max_length":24.6,
               "percentage_of_surface":0.11,
               "scar_angle":"NaN",
               "polygon_count":2
            }
         ]
      },
      {
         "surface_id":3,
         "classification":"Platform",
         "total_area_px":50040.0,
         "total_area":78.4,
         "max_breadth":20.0,
         "max_length":6.3,
         "polygon_count":5,
         "scar_count":0,
         "percentage_detected_scars":0.0,
         "scar_contours":[
         ]
      }
   ]
}
```---
title: 'PyLithics: A Python package for stone tool analysis'
tags:
  - Python
  - Human evolution
  - Archaeology
  - Lithic analysis
  - Prehistoric technology
  - Computer vision
authors:
  - name: Jason J. Gellis^[corresponding author]
    orcid: 0000-0002-9929-789X
    affiliation: 1, 2, 3
  - name: Camila Rangel Smith
    orcid: 0000-0002-0227-836X
    affiliation: 1
  - name: Robert A. Foley
    orcid: 0000-0003-0479-3039
    affiliation: 1, 2, 3
affiliations:
  - name: The Alan Turing Institute
    index: 1
  - name: University of Cambridge
    index: 2
  - name: Leverhulme Centre for Human Evolutionary Studies
    index: 3
date: 3rd September 2021
bibliography: paper.bib
---
  
# Summary

Archaeologists have long used stone tools (lithics) to reconstruct the behaviour of prehistoric hominins. While 
techniques have become more quantitative, there still remain barriers to optimizing data retrieval [@Andrefsky:2012]. 
Machine learning and computer vision approaches can be developed to extract quantitative and trait data from lithics, 
photographs and drawings. `PyLithics` has been developed to capture data from 2D line drawings, focusing on the size, 
shape and technological attributes of flakes. The problems addressed in the software are: one, capturing data in a form
that can be quantified, and information maximized; two, solving the challenges of data that is not a simple linear 
sequence of bases but complex 3D objects or 2D image representations; and three, transforming and exporting these into 
systematic data for analysis. The goal is to enhance the size and quality of lithic databases for analysing ancient 
technology and human behaviour.

# Statement of need

`PyLithics` is an open-source, free for use software package for processing lithic artefact illustrations scanned from 
the literature. Accurately measuring lithic artefacts is difficult and especially time-consuming as lithics and their 
features are incongruous shapes and sizes. This is especially problematic for the researcher as certain features, such 
as flake scar size, are useful in elucidating the manufacturing process of an artefact. Thus, while even the best, 
most complete illustrations are able to visually capture an immense amount of information about an artefact, much of 
this information is under-utilized or not used at all.

`PyLithics` alleviates these issues by accurately identifying, outlining, and computing lithic shape and linear 
measures, and returns user ready data. It has been optimized for feature extraction and measurement using a number of 
computer vision techniques including pixel intensity thresholding, edge detection, contour finding, custom template 
matching and image kernels. On both conventional and modern drawings, `PyLithics` can identify and label platform, 
lateral, dorsal, and ventral surfaces, as well as individual dorsal surface scar shape, size, orientation, diversity, 
number, and order of flake size from greatest to least. Complete size and shape metrics of individual scars and whole 
flakes can be calculated and recorded. Orientation and flaking direction of dorsal scars can also be calculated. The 
resulting data can be used for metrical analysis, extracting features indicative of both typologies and technological 
processes. Data output can easily be employed to explore patterns of variation within and between assemblages.

# Methods and workflow

`PyLithics` is devised to work with illustrations of lithic objects common to publications in archaeology and 
anthropology. Lithic illustrators have established conventions regarding systems of artefact orientation and 
proportions. Stone artefacts made by hominins have clear diagnostic features of manufacture that are used to infer 
technological patterns. A basic division is made between cores, which are the source material that is modified by 
striking, and flakes (also known as debitage), which are the pieces struck from the cores. Both provide important 
information, but flakes have been the focus for `PyLithics` development. 

Lithics are normally drawn at a 1:1 scale, with the vertical axis orthogonal to the striking platform – that is, the 
point on the flake which was struck to remove it from the core. A preferred method is to orient and illustrate various 
aspects of a flake as a series of adjacent surfaces at 90-degree rotations. Four illustrations are normal: the ventral 
surface, the surface that was last attached to the core, and is usually a smooth, unbroken surface; the dorsal surface, 
which is usually marked by scars from previous flake removals; the proximal or platform view, from the top where the 
flake was struck, and which is approximately at 90 degrees to the ventral surface; and the lateral view, which 
essentially shows the thickness of the flake and its longitudinal shape. Each aspect, but especially the dorsal surface, 
contains internal details (i.e., flake scars, cortical areas, etc.), indication of flaking direction radial lines 
(ripples), and the inclusion of a metric scale (for more information about lithic drawings see [@Martingell:1988]). 
Currently `PyLithics` is optimised to work with unifacial flakes and bifaces, which are relatively flat, two-dimensional 
objects. For best performance and accurate measurement, images loaded into `PyLithics` should be:

1.	Oriented with the platform orthogonal to the top of the page.
2.	Use arrows to indicate flaking direction rather than ripples.
3.	Have arrows with longer tails than arrow points.
4.	Avoid having arrows or other indicators overlapping or immediately adjacent to flake scars.


While `PyLithics` can identify and measure surfaces and surface features on illustrations with ripples, too many ripples 
will reduce measurement accuracy of flake scars. The inputs for `PyLithics` are images of lithic objects, images of their associated scales, and a metadata CSV file 
linking the two and giving the scale measurement in millimetres. `PyLithics` processes the images with the following steps 
and as illustrated in (\autoref{fig:Figure_1}):


1.	Import images and match image name to associated image ID and scale image from CSV metadata file.
2.	Calculate a conversion of pixels to millimetres based on the size of the associated scale from CSV metadata file. If no scale is present, measurements will be in pixels
3.	Apply noise removal and contrast stretching to images to minimise pixel variation.
4.	Pixel intensity thresholding of images to prepare for contour finding.
5.	Apply edge detection and contour finding to thresholded images.
6.	Calculate metrics of lithic surface features from found contours -- area, length, breath, shape, number of vertices.
7.	Select contours which outline an entire lithic object's surfaces or select contours of inner scars greater than 3% and less than 50% of the total size of its surface.
8.	Classify these selected surface contours as "Dorsal", "Ventral", "Lateral", and/or "Platform" depending on presence or absence. Assign scar contours to these surfaces.
9.	If present, find arrows using connected components and template matching, measure their angle and assign angle to associated scar.
10.	Plot resulting surface and scar contours on the original images for validation.
11.	Output data in a hierarchical JSON file detailing measurement of surface and scar contours.

![PyLithics program workflow.\label{fig:Figure_1}](pylithics_flowchart.jpg)

`PyLithics` depends on common Python packages such as NumPy [@Harris:2020], SciPy [@Virtanen:2020], Pandas [@McKinney:2010] 
for data processing, Matplotlib [@Hunter:2007] for plotting and scikit-image [@scikit-image] and OpenCv [@opencv_library] 
for image processing and computer vision tasks.

# Results

`PyLithics` generates two outputs:

1) An image set comprised of the original input images with superimposed contour identification and derived metrics (see \autoref{fig:Figure_2} and \autoref{fig:Figure_3}), and if arrows are present in the illustration, angles of flaking direction (see \autoref{fig:Figure_4}).
2) A JSON file with data for lithic objects and surface features found in each image \autoref{fig:Figure_5}. These data are hierarchically organised, first by type of object surface (i.e., ventral, dorsal, lateral, and platform); and second by metrics (in mm) from scars and arrows associated to each object surface. Output includes: object id, pixel conversion rate (based on provided scale), number of surfaces identified (e.g., dorsal, ventral, etc.); surface and scar area, maximum breadth and length, scar count (ordered by largest to smallest), percentage of surface area represented by scars, surface area of individual scars, and flaking angle of scars. `PyLithics` also produces a polygon count – a way of characterising the shape of each scar by approximating a polygon to it and counting how many sides it has. Output files can be used for statistical analyses with any software that can read JSON files.

Output images \autoref{fig:Figure_2} serve as validation of the output data \autoref{fig:Figure_3}.

![PyLithics output - detected surfaces.\label{fig:Figure_2}](rub_al_khali_lithic_surfaces.png ){ width=100% }

![PyLithics output - detected scars.\label{fig:Figure_3}](rub_al_khali_lithic_scars.png ){ width=100% }

![PyLithics output - flake angles.\label{fig:Figure_4}](rub_al_khali_lithic_angles.png){ width=100% }

![PyLithics output - JSON data file.\label{fig:Figure_5}](rub_al_khali_JSON.png){ width=100% }

# Outlook 

Evolutionary biology, and the study of human evolution in particular, has been transformed by the impact of genomics and 
the development of ancient DNA methodologies [@Moody:2004]. One of the reasons that genomics has had such an impact is 
the sheer scale of the data now available, and power of the analytical techniques used. Although current approaches to 
lithic analysis have become more quantitative, they remain based on relatively univariate attribute assignments and 
limited metrics, variably collected and reported. `PyLithics` aims to expand data collection with the goal of building 
expansive, comprehensive, and standardized high-dimensional lithic artefact datasets for integration with genomic and 
fossil data.

# Acknowledgements

The `PyLithics` package was developed by researchers from The Alan Turing Institute, University of Cambridge, and the Leverhulme Centre for Human Evolutionary Studies. Funding was provided by the Alan Turing Institute (grant number G109254). We would like to acknowledge support from Professor Katharine Robson Brown and Doctor Sebastian Ahnert.

# References
