# Release Notes


## Version 0.1.1 (in progress)

### Fixed
- Corrected labels for symmetric atlas hemisphere stack generation

### Added
- hemisphere values as `Atlas` class attributes


## Version 0.1.0

First release of the BrainGlobe `bg-atlasapi` package! :tada:

### Added
- Atlas metadata, reference stacks, region hierarchy
- Configuration tools
- Atlas download and updating
# BG-atlasAPI

[![Python Version](https://img.shields.io/pypi/pyversions/bg-atlasapi.svg)](https://pypi.org/project/bg-atlasapi)
[![PyPI](https://img.shields.io/pypi/v/bg-atlasapi.svg)](https://pypi.org/project/bg-atlasapi/)
[![Wheel](https://img.shields.io/pypi/wheel/bg-atlasapi.svg)](https://pypi.org/project/bg-atlasapi)
[![Development Status](https://img.shields.io/pypi/status/brainatlas-api.svg)](https://github.com/SainsburyWellcomeCentre/brainatlas-api)
[![Tests](https://img.shields.io/github/workflow/status/brainglobe/bg-atlasapi/tests)](
    https://github.com/brainglobe/bg-atlasapi/actions)
[![Coverage Status](https://coveralls.io/repos/github/brainglobe/brainatlas-api/badge.svg?branch=master)](https://coveralls.io/github/brainglobe/brainatlas-api?branch=master)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02668/status.svg)](https://doi.org/10.21105/joss.02668)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

![banner](https://github.com/brainglobe/bg-atlasapi/blob/master/docs/img/banner.jpeg)


The brainglobe atlas API (BG-AtlasAPI) provides a common interface for programmers to download and process brain atlas data from multiple sources. 

## Atlases available

A number of atlases are in development, but those available currently are:
* The [Allen Mouse Brain Atlas](https://www.brain-map.org) at 10, 25, 50 and 100 micron resolutions
* The [Allen Human Brain Atlas](https://www.brain-map.org) at 100 micron resolution
* The [Max Planck Zebrafish Brain Atlas](http://fishatlas.neuro.mpg.de) at 1 micron resolution
* The [Enhanced and Unified Mouse Brain Atlas](https://kimlab.io/brain-map/atlas/) at 10, 25, 50 and 100 micron resolutions
* The [Smoothed version of the Kim et al. mouse reference atlas](https://doi.org/10.1016/j.celrep.2014.12.014) at 10, 25, 50 and 100 micron resolutions

## Installation
BG-AtlasAPI works with Python >3.6, and can be installed from PyPI with:
```bash
pip install bg-atlasapi
```

## Usage
Full information can be found in the [documentation](https://docs.brainglobe.info/bg-atlasapi/usage)
### Python API
**List of atlases**

To see a list of atlases use `bg_atlasapi.show_atlases`
```python
from bg_atlasapi import show_atlases
show_atlases()
#                                Brainglobe Atlases                               
# ╭──────────────────────────────────┬────────────┬───────────────┬──────────────╮
# │ Name                             │ Downloaded │ Local version │    Latest    │
# │                                  │            │               │   version    │
# ├──────────────────────────────────┼────────────┼───────────────┼──────────────┤
# │ allen_human_500um                │     ✔      │      0.1      │     0.1      │
# │ mpin_zfish_1um                   │     ✔      │      0.3      │     0.3      │
# │ allen_mouse_50um                 │     ✔      │      0.3      │     0.3      │
# │ kim_unified_25um                 │     ✔      │      0.1      │     0.1      │
# │ allen_mouse_25um                 │     ✔      │      0.3      │     0.3      │
# │ allen_mouse_10um                 │     ✔      │      0.3      │     0.3      │
# │ example_mouse_100um              │    ---     │      ---      │     0.3      │
# ╰──────────────────────────────────┴────────────┴───────────────┴──────────────╯
```

**Using the atlases**

All the features of each atlas can be accessed via the `BrainGlobeAtlas` class. 


e.g. for the 25um Allen Mouse Brain Atlas:

```python
from bg_atlasapi.bg_atlas import BrainGlobeAtlas
atlas = BrainGlobeAtlas("allen_mouse_25um")
```

The various files associated with the atlas can then be accessed as attributes of the class:
```python
# reference image
reference_image = atlas.reference
print(reference_image.shape)
# (528, 320, 456)

# annotation image
annotation_image = atlas.annotation
print(annotation_image.shape)
# (528, 320, 456)

# a hemispheres image (value 1 in left hemisphere, 2 in right) can be generated
hemispheres_image = atlas.hemispheres
print(hemispheres_image.shape)
# (528, 320, 456)
```

**Brain regions**

There are multiple ways to work with individual brain regions. To see a dataframe of each brain region, with it's unique ID, acronym and full name, use `atlas.lookup_df`:
```python
atlas.lookup_df.head(8)
#      acronym         id                           name
# 0       root        997                           root
# 1       grey          8  Basic cell groups and regions
# 2         CH        567                       Cerebrum
# 3        CTX        688                Cerebral cortex
# 4      CTXpl        695                 Cortical plate
# 5  Isocortex        315                      Isocortex
# 6        FRP        184  Frontal pole, cerebral cortex
# 7       FRP1         68          Frontal pole, layer 1
```

Each brain region can also be access by the acronym, e.g. for primary visual cortex (VISp):
```python
from pprint import pprint
VISp = atlas.structures["VISp"]
pprint(VISp)
# {'acronym': 'VISp',
#  'id': 385,
#  'mesh': None,
#  'mesh_filename': PosixPath('/home/user/.brainglobe/allen_mouse_25um_v0.3/meshes/385.obj'),
#  'name': 'Primary visual area',
#  'rgb_triplet': [8, 133, 140],
#  'structure_id_path': [997, 8, 567, 688, 695, 315, 669, 385]}
```

### Note on coordinates in `bg-atlasapi` 
Working with both image coordinates and cartesian coordinates in the same space can be confusing! In `bg-atlasapi`, the origin is always assumed to be in the upper left corner of the image (sectioning along the first dimension), the "ij" convention. This means that when plotting meshes and points using cartesian systems, you might encounter confusing behaviors coming from the fact that in cartesian plots one axis is inverted with respect to  ij coordinates (vertical axis increases going up, image row indexes increase going down). To make things as consistent as possible, in `bg-atlasapi` the 0 of the meshes coordinates is assumed to coincide with the 0 index of the images stack, and meshes coordinates increase following the direction stack indexes increase.
To deal with transformations between your data space and `bg-atlasapi`, you might find helpful the [bg-space](https://github.com/brainglobe/bg-space) package.

# Contributing to bg-atlasapi
**Contributors to bg-atlaspi are absolutely encouraged**, whether you want to fix bugs, add/request new features or simply ask questions.

If you would like to contribute to `bg-atlasapi` (or any of the downstream tools like [brainrender](https://github.com/BrancoLab/BrainRender) etc.\) please get in touch by opening a new issue or pull request on [GitHub](https://github.com/brainglobe/bg-atlasapi). 

Someone might have already asked a question you might have, so if you're not sure where to start, check out the [issues](https://github.com/brainglobe/bg-atlasapi/issues) \(and the issues of the other repositories\)

---
The BrainGlobe project is generously supported by the Sainsbury Wellcome Centre and the Institute of Neuroscience, Technical University of Munich, with funding from Wellcome, the Gatsby Charitable Foundation and the Munich Cluster for Systems Neurology - Synergy.

<img src='https://brainglobe.info/images/logos_combined.png' width="550">

# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

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
reported by contacting the project team at adam.tyson@ucl.ac.uk. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
