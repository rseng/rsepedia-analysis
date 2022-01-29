[![Project Status](https://img.shields.io/pypi/status/kalmus.svg)](https://pypi.org/project/kalmus/)
[![Python Version](https://img.shields.io/pypi/pyversions/kalmus.svg)](https://pypi.org/project/kalmus/)
[![PyPI Version](https://img.shields.io/pypi/v/kalmus.svg)](https://pypi.org/project/kalmus/)
[![codecov](https://codecov.io/gh/KALMUS-Color-Toolkit/KALMUS/branch/master/graph/badge.svg)](https://codecov.io/gh/KALMUS-Color-Toolkit/KALMUS)
[![License](https://img.shields.io/pypi/l/kalmus.svg)](https://pypi.org/project/kalmus/)
[![status](https://joss.theoj.org/papers/f7a87aac389fc3cd02807d5fad6ebf50/status.svg)](https://joss.theoj.org/papers/f7a87aac389fc3cd02807d5fad6ebf50)

# KALMUS

KALMUS is a Python package for the computational analysis of colors in films. 
It provides quantitative tools to study and compare the use of film color. 
This package serves two purposes: (1) various ways to measure, calculate and compare a film's colors 
and (2) various ways to visualize a film's color. We have named the software KALMUS in homage to 
Natalie Kalmus (1882 - 1965), a Technicolor Director who oversaw the color palettes of nearly 300 
Hollywood feature films.

KALMUS utilizes the movie barcode as a visualization of the film's color. It has a modularized pipeline for the
 generation of barcodes using different measures of color and region of interest in each film frame. KALMUS provides
 a low-level API, high-level command line, and Graphic user interface for audience from all backgrounds to take
 advantage of its functionality. 

- What is a Movie Barcode: [KALMUS: tools for color analysis of films](https://joss.theoj.org/papers/10.21105/joss.03156). *Journal of Open Source Software*, 6(61), 3156. [https://doi.org/10.21105/joss.03156]( https://doi.org/10.21105/joss.03156)      
- How do I install the KALMUS: [KALMUS Installation Guide](https://kalmus-color-toolkit.github.io/KALMUS/install.html)
- How do I use the KALMUS: [Notebook Tutorials for KALMUS's API, GUI, and CLI](https://github.com/KALMUS-Color-Toolkit/KALMUS/tree/master/notebooks)
- How do I contribute to the KALMUS: [KALMUS Contribution Guidelines](https://github.com/KALMUS-Color-Toolkit/KALMUS/blob/master/CONTRIBUTING.md)

# API Documentation

The KALMUS API reference is now available on 
[https://kalmus-color-toolkit.github.io/KALMUS/kalmus.html](https://kalmus-color-toolkit.github.io/KALMUS/kalmus.html).

# Update Log

The full update log (from [v1.3.0](https://pypi.org/project/kalmus/1.3.0/) to [v1.3.9](https://pypi.org/project/kalmus/)) 
is now available on [https://kalmus-color-toolkit.github.io/KALMUS/update_log.html](https://kalmus-color-toolkit.github.io/KALMUS/update_log.html)
[![Project Status](https://img.shields.io/pypi/status/kalmus.svg)](https://pypi.org/project/kalmus/)
[![Python Version](https://img.shields.io/pypi/pyversions/kalmus.svg)](https://pypi.org/project/kalmus/)
[![PyPI Version](https://img.shields.io/pypi/v/kalmus.svg)](https://pypi.org/project/kalmus/)
[![status](https://joss.theoj.org/papers/f7a87aac389fc3cd02807d5fad6ebf50/status.svg)](https://joss.theoj.org/papers/f7a87aac389fc3cd02807d5fad6ebf50)
[![codecov](https://codecov.io/gh/KALMUS-Color-Toolkit/KALMUS/branch/master/graph/badge.svg)](https://codecov.io/gh/KALMUS-Color-Toolkit/KALMUS)
[![License](https://img.shields.io/pypi/l/kalmus.svg)](https://pypi.org/project/kalmus/)

# KALMUS

KALMUS is a Python package for the computational analysis of colors in films. 
It provides quantitative tools to study and compare the use of film color. 
This package serves two purposes: (1) various ways to measure, calculate and compare a film's colors 
and (2) various ways to visualize a film's color. We have named the software KALMUS in homage to 
Natalie Kalmus (1882 - 1965), a Technicolor Director who oversaw the color palettes of nearly 300 
Hollywood feature films.

KALMUS utilizes the movie barcode as a visualization of the film's color. It has a modularized pipeline for the
 generation of barcodes using different measures of color and region of interest in each film frame. KALMUS provides
 a low-level API, high-level command line, and Graphic user interface for audience from all backgrounds to take
 advantage of its functionality. 

- What is a Movie Barcode: [KALMUS: tools for color analysis of films](https://joss.theoj.org/papers/10.21105/joss.03156). *Journal of Open Source Software*, 6(61), 3156. [https://doi.org/10.21105/joss.03156]( https://doi.org/10.21105/joss.03156)   
- How do I install the KALMUS: [KALMUS Installation Guide](https://kalmus-color-toolkit.github.io/KALMUS/install.html) 
and [KALMUS PyPI Homepage](https://pypi.org/project/kalmus/).
- How do I use the KALMUS: [Notebook Tutorials for KALMUS's API, GUI, and CLI](notebooks)
- How do I contribute to the KALMUS: [KALMUS Contribution Guidelines](CONTRIBUTING.md)
- How do I run the KALMUS's automated test suite: [Auomated Test Suite](tests/)

**Examples of Barcode visualization:**

<p align="center">
  <img width="" height="" src="notebooks/notebook_figures/mission_barcode_whole_frame_avg.png">
  <br>Figure 1. Mission: Impossible (1996) color barcode using the average color of whole frame for each frame</br>
  <br>
  <img width="" height="" src="notebooks/notebook_figures/mission_barcode_Foreground_avg.png">
  <br>Figure 2. Mission: Impossible (1996) color barcode using the average color of foreground of each frame</br>
</p>

# API Documentation

The KALMUS API reference is now available on 
[https://kalmus-color-toolkit.github.io/KALMUS/kalmus.html](https://kalmus-color-toolkit.github.io/KALMUS/kalmus.html).

# Installation Guide
[![Python Version](https://img.shields.io/pypi/pyversions/kalmus.svg)](https://pypi.org/project/kalmus/)
[![PyPI Version](https://img.shields.io/pypi/v/kalmus.svg)](https://pypi.org/project/kalmus/)
[![build workflow](https://github.com/KALMUS-Color-Toolkit/KALMUS/actions/workflows/python-package.yml/badge.svg)](https://github.com/KALMUS-Color-Toolkit/KALMUS/actions/workflows/python-package.yml)

The kalmus package requires a python with version 3.7 or 3.8.

The package is released on PyPI ([Project Homepage](https://pypi.org/project/kalmus/)). After you installed the
python==3.7, 3.8, you can install the kalmus using pip (recommended)

    $ pip install kalmus


Alternatively, you could install the kalmus locally by first cloning this GitHub repo.
Then, move to the top directory of cloned kalmus project folder and install using the pip command

    $ pip install .

In both methods, the package's dependencies will be automatically installed. You can verify if the kalmus has been
installed in your environment using the pip command

    $ pip show kalmus

Alternatively, in version 1.3.7 and above, you can check the version of installed kalmus using its 
`.__version__` attribute.

```jupyter
>>> import kalmus
>>> print(kalmus.__version__) # Warning: The __version__ attribute is not available in the kalmus v.1.3.6 and backward
>>> 1.3.7 
```

## For users with Apple M1 Chip (arm64 Architecture)

As @elektrobohemian mentioned in [issue #4](https://github.com/KALMUS-Color-Toolkit/KALMUS/issues/4), kalmus cannot build natively on Apple M1 processors because of kalmus's dependencies on NumPy. You may be able to install kalmus under a Rosetta emulation with Python 3.7. 

# Get Started

KALMUS has a low-level API, high-level command line, and Graphic user interface for audience from all 
backgrounds to take advantage of its functionality. 

To get started on KALMUS, we encourage you to check the Jupyter notebook tutorials in the [notebooks](notebooks) 
folder. We provide the interactive notebook tutorials for users to get started on KALMUS using its API, GUI, and CLI. 
Notice that the Command-line interface (CLI) is only available in KALMUS v1.3.7 or onward.

- [Notebook Tutorial for Graphic User Interface](notebooks/user_guide_for_kalmus_gui.ipynb)
- [Markdown Tutorial for Graphic User Interface](notebooks/USAGE_GRAPHIC_USER_INTERFACE.md)
- [Notebook Tutorial for Application Programming Interface](notebooks/user_guide_for_kalmus_api.ipynb)
- [Markdown Tutorial for Command-line interface](notebooks/USAGE_COMMAND_LINE_UI.md)
- [Notebook Tutorial for Advanced API Usage](notebooks/advanced_guide_for_kalmus_api.ipynb)

# Contribution

We encourage contributions, including bug fixes and new features, from our community users. When contributing to the 
kalmus package, please contact the project maintainers by email <yc015@bucknell.edu> or opening an issue. If 
your bug fixes or new features change the current behaviors of package, please specify the changes and reasons in the 
discussion with project maintainers. 

We encourage inclusive and friendly discussion. Please follow our [code of conduct](CODE_OF_CONDUCT.md) when 
communicating. 

# Test Suite
[![codecov](https://codecov.io/gh/KALMUS-Color-Toolkit/KALMUS/branch/master/graph/badge.svg)](https://codecov.io/gh/KALMUS-Color-Toolkit/KALMUS)
[![codecov workflow](https://github.com/KALMUS-Color-Toolkit/KALMUS/actions/workflows/test-codecov.yml/badge.svg)](https://github.com/KALMUS-Color-Toolkit/KALMUS/actions/workflows/test-codecov.yml)

We provide an [automated test suite](tests/) that covers the core functionality of KALMUS. Before running the automated test suite locally, 
make sure you have installed the latest versions of [pytest](https://pypi.org/project/pytest/), [pytest-cov](https://pypi.org/project/pytest-cov/), 
and [kalmus](https://pypi.org/project/kalmus/), and you have cloned the project repository on master branch. 

To run the test suite:  
- Go to the top directory of cloned KALMUS project
- Use command `$ python -m pytest tests --cov=kalmus --cov-config=.coveragerc --cov-report term-missing`

See the [Test Suite Guide](tests/README.md) for more details.

# Citation
[![status](https://joss.theoj.org/papers/f7a87aac389fc3cd02807d5fad6ebf50/status.svg)](https://joss.theoj.org/papers/f7a87aac389fc3cd02807d5fad6ebf50)

If you find our software is useful in your work, please cite our paper that describes the usage of KALMUS in the computational analysis of colors in films. DOI: [https://doi.org/10.21105/joss.03156](https://doi.org/10.21105/joss.03156)

> Chen et al., (2021). KALMUS: tools for color analysis of films. Journal of Open Source Software, 6(61), 3156, https://doi.org/10.21105/joss.03156

Here is the BibTex citation of our work:

    @article{Chen2021,
        doi = {10.21105/joss.03156},
        url = {https://doi.org/10.21105/joss.03156},
        year = {2021},
        publisher = {The Open Journal},
        volume = {6},
        number = {61},
        pages = {3156},
        author = {Yida Chen and Eric Faden and Nathan C. Ryan},
        title = {KALMUS: tools for color analysis of films},
        journal = {Journal of Open Source Software}
    }

# Acknowledgment

The authors wish to thank the Mellon Foundation, the Dalal Family Foundation, and the Bucknell University Humanities 
Center for their support on this project. The project is released under the open-source MIT License.

# Update Log

The full update log (from [v1.3.0](https://pypi.org/project/kalmus/1.3.0/) to [v1.3.9](https://pypi.org/project/kalmus/)) 
is now available on [https://kalmus-color-toolkit.github.io/KALMUS/update_log.html](https://kalmus-color-toolkit.github.io/KALMUS/update_log.html)
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
reported by contacting the project team at yc015@bucknell.edu. All
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
# Contributing
We encourage contributions, including bug fixes and new features, from our community users. When contributing to the kalmus package, please contact the project maintainers by email <yc015@bucknell.edu> or opening an issue. If your bug fixes or new features change the current behaviors of package, please specify the changes and reasons in the discussion with project maintainers. 

We encourage inclusive and friendly discussion. Please follow our [code of conduct](CODE_OF_CONDUCT.md) when communicating. 
[![status](https://joss.theoj.org/papers/f7a87aac389fc3cd02807d5fad6ebf50/status.svg)](https://joss.theoj.org/papers/f7a87aac389fc3cd02807d5fad6ebf50)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4765745.svg)](https://doi.org/10.5281/zenodo.4765745)

# The Journal of Open Source Software Paper
This folder contains the [raw markdown paper](joss-paper.md) which we submitted to the [Journal of Open Source Software](https://joss.theoj.org/papers/10.21105/joss.03156). The review of the paper and software can be found here: [[REVIEW]: KALMUS: tools for color analysis of films #3156](https://github.com/openjournals/joss-reviews/issues/3156). 

The specific version of KALMUS reviewed by the JOSS is archived on the Zenodo: [https://zenodo.org/record/4765745#.YKP4vKhKhPY](https://zenodo.org/record/4765745#.YKP4vKhKhPY). We recommend our users, however, use the latest version of KALMUS published on [PyPI](https://pypi.org/project/kalmus/).

For our contributors, please do not make any changes to the files in this folder.

# Citation
If you find our software is useful in your work, please cite our paper that describes the usage of KALMUS in the computational analysis of colors in films. DOI: [https://doi.org/10.21105/joss.03156](https://doi.org/10.21105/joss.03156)

> Chen et al., (2021). KALMUS: tools for color analysis of films. Journal of Open Source Software, 6(61), 3156, https://doi.org/10.21105/joss.03156

Here is the BibTex citation of our work:

    @article{Chen2021,
        doi = {10.21105/joss.03156},
        url = {https://doi.org/10.21105/joss.03156},
        year = {2021},
        publisher = {The Open Journal},
        volume = {6},
        number = {61},
        pages = {3156},
        author = {Yida Chen and Eric Faden and Nathan C. Ryan},
        title = {KALMUS: tools for color analysis of films},
        journal = {Journal of Open Source Software}
    }
---
title: 'KALMUS: tools for color analysis of films'
tags:
  - Python
  - computer vision
  - color
authors:
  - name: Yida Chen
    affiliation: 1 
  - name: Eric Faden
    affiliation: 1
  - name: Nathan C. Ryan
    orcid: 0000-0003-4947-586X
    affiliation: 1
affiliations:
  - name: Bucknell University
    index: 1
date: 15 January 2021
bibliography: joss-refs.bib
---


# Summary

KALMUS is a Python package for the computational analysis of colors in films. It provides quantitative tools to study and compare the use of film color.  This package serves two purposes:  (1) various ways to measure, calculate and compare a film's colors and (2) various ways to visualize a film's color.  We have named the software KALMUS in homage to Natalie Kalmus (1882 - 1965), a Technicolor Director who oversaw the color palettes of nearly 300 Hollywood feature films.


# Statement of Need

>“Colors are elusive”

Barbara Flueckiger [-@flueckiger2017digital]

>“Color is in so many ways uncontainable.”

Josh Yumibe [-@yumibe2012moving]


The epigraphs above (each from the introduction of recent research on cinematic color) all acknowledge the challenges of studying cinematic color aesthetics.  Filmmakers, cinematographers, production designers, and colorists spend significant money, time, and labor deciding a film’s color palette before, during, and after production.  Given this obsessive attention, why is color so “elusive” to study and analyze?  

Three reasons emerge: (1) film color remains subtle and almost subliminal by design; (2) we read color subjectively and are often influenced by cultural or even historic circumstances; (3) films are kinetic with the image (and colors) constantly moving.  Color might be relatively easy to observe in a single frame, but tracking the color palette across an entire film proves daunting.  Moreover, comparing one film’s palette to another proves especially challenging due to the accumulation of data needed to make a meaningful comparison.  See [@flueckiger2020, @stutz] for an overview of both the history and the state-of-the-art of the study of colors in film.

Recently, several digital tools have emerged that analyze a film’s color.  The ability to capture film stills from DVD or streaming sources and analyze the frame’s palette in off-the-shelf tools like [Adobe Color](http://color.adobe.com) or Adobe Photoshop provide insights into individual frames.  Additionally, many open-source projects provide implementations of state-of-the-art computer vision algorithms applied to moving images.  Some, such as the Distance Viewing Toolkit [@arnold], are very general and handle several aspects of a moving image at once (sound, color, camera angle, etc.).  Others, such as VIAN [@vian], allow for a "closer", more interactive distant viewing.  

KALMUS is built around a visualization of films known as movie barcodes; it produces an overall color palette “snapshot” of a film.  Artists like [Jeffrey Moser](http://www.jeffreymoser.com/) or programmers such as [Charlie Clark](https://thecolorsofmotion.com/about), generate movie barcodes by reducing each film frame to a single color[^1] and then stitching these colors into a mosaic.  See Figure 1 for examples of movie barcodes generated by KALMUS.  Existing movie barcode software is inadequate for studying films in a quantitative way.  While such software provides a fascinating visual “overview” of a film’s color palette, it, unlike KALMUS, focuses on visualization rather than providing meaningful analytical data.  KALMUS generates barcodes but also generates a wealth of data related to color and the kinds of barcode and data that can be generated are fairly customizable so that users can carry out their own analyses.  

![Movie barcodes generated by KALMUS.](images-joss/kalmus_figure1_1.jpg)

[^1]: In OpenCV [@opencv], for example, if a frame is decimated to a single pixel, its color is determined by area interpolation, a weighted average of the RGB-values of the pixels.

KALMUS allows for the analysis of a film's color by:  

1. Providing an interface that takes in either a video file or a JSON file (a sample JSON file can be found [here](https://github.com/KALMUS-Color-Toolkit/KALMUS/blob/master/kalmus/data/mission_impossible_Bright_Whole_frame_Color.json)).  
1. Allowing for the computation of a frame's color in a number of ways (dominant, median, etc).  
1. Allowing for the computation of a film's color in a number of ways.  
1. Providing implementations of ways to compare colors of two films.  
1. Providing implementations of ways to visualize the color of a film.  
1. Allowing the user to download color data as a CSV or JSON file.  See Figure 2 for an explanation of KALMUS's user interface.  

![Image A shows barcodes for a single film in which the color of a frame was calculated in two ways:  the top one uses the frame's brightest color and the bottom one uses the frame's median color.  Image B shows the interface for generating a barcode: one can process a film in parallel, choose a sampling rate, and choose a method for determining a frame's color.  Image C is the result of selecting a pixel in a barcode to give a user an idea of where in the film they are.  Image D is a table of similarity metrics between the two barcodes.  Image E is an interactive 3D plot of the RGB values in one of the barcodes.](images-joss/kalmus-interface.jpg)



KALMUS allows users to understand a film’s color palette and compare that palette to other films.  See Table 1 for options on what parts of the frame a user might want to analyze and the various metrics a user can use to determine the color of a frame.  See Table 2 for information on how film colors can be compared.  Potentially, films from a certain time period or a particular genre might be aggregated to see if there are common color palettes[^2].  While film studies as a discipline has long depended on qualitative analysis, KALMUS provides a straightforward quantitative tool that can supplement the reading of a film.  For example, in @adams KALMUS is used to study color trends in Hollywood movies from 1990 to 2015.  In a class taught by the second author, students used KALMUS to explore how a film's color palette signaled narrative shifts and introduced significant characters.

[^2]: We point out users should be aware of the fact that different color gamuts used in digitizing films may affect a film's original color palette.  KALMUS is designed to use the broadcast color gamut called REC 709.  


|                     | Mean	| Median   | Mode	| TD	| WD	| BR	| BP  |
| :-------------------| ----: | -------: | ---: | --: | --: |---: | --: |
| Whole frame	        | Yes   | Yes      | Yes  | Yes | Yes | Yes | Yes |
| High contrast region	| Yes | Yes | Yes | Yes | Yes | Yes | No |
| Low contrast region	| Yes | Yes | Yes | Yes | Yes |  Yes |No |
| Foreground	| Yes | Yes | Yes | Yes | Yes |  Yes | No |
| Background	| Yes | Yes | Yes | Yes | Yes |  Yes | No |

Table 1:  A summary of the various metrics one can use to determine the color of a frame and what parts of a frame will be used to calculate the color of the frame.  The first three metrics (mean, median and mode) are the corresponding statistic over all the colors in the part of the frame being analyzed.  The next two metrics, top dominant (TD) and weighted dominant (WD) are determined by applying a clustering algorithm on the set of all the colors in the part of the frame being analyzed and, in the case of the top dominant metric, the color is determined by the color of the largest cluster and, in the case of the weighted dominant metric, the color is determined by the average color of the clusters, weighted by the sizes of the clusters.  The last two metrics assign the color of the region being analyzed to be the color of the brightest region (BR) or the brightest pixel (BP) in the region.

\small
| Comparison metric | Range | References |
| :---------------- | ----: | --------: |
| Normalized root mean square error | 0 least similar, 1 most similar | @wang  |
| Structural similarity index | 0 least similar, 1 most similar | @wang |
| Cross correlation | -1 anti-similar, 1 most similar | @avants|
| Local cross correlation | -1 anti-similar, 1 most similar | @avants |
| Needleman--Wunsch | 0 least similar, 1 most similar | @needleman |
|                   |                                 | @adams |
| Smith--Waterman | 0 least similar, 1 most similar | @smith |
|                 |                                 | @adams |
\normalsize
Table 2:  A summary of the various metrics included in KALMUS to compare the overall color of two films.


As of version 1.3.5, the software is stable.  Because the audience of potential users of KALMUS is broad, we provide access to the KALMUS package via the GUI described above but also point out that the package can be installed and used in any Python environment.

# Acknowledgments
The authors wish to thank the Mellon Foundation, the Dalal Family Foundation, and the Bucknell University Humanities Center for their support on this project.


# References

[![codecov](https://codecov.io/gh/KALMUS-Color-Toolkit/KALMUS/branch/master/graph/badge.svg)](https://codecov.io/gh/KALMUS-Color-Toolkit/KALMUS)
[![codecov workflow](https://github.com/KALMUS-Color-Toolkit/KALMUS/actions/workflows/test-codecov.yml/badge.svg)](https://github.com/KALMUS-Color-Toolkit/KALMUS/actions/workflows/test-codecov.yml)

# Automated Test Suite

We provide an automated tests suite for you to validate the package's core functionality.
The modules being tested including:  
```python
from kalmus.barcodes.Barcode import *
from kalmus.barcodes.BarcodeGenerator import *
from kalmus.utils.artist import *
from kalmus.utils.measure_utils import *
from kalmus.utils.visualization_utils import *
from kalmus.tkinter_windows.gui_utils import *
from kalmus.command_line_generator import *
```

A [GitHub Action](../.github/workflows/test-codecov.yml) will run on every push or pull-request to the master branch 
and upload the coverage report on [Codecov](https://app.codecov.io/gh/KALMUS-Color-Toolkit/KALMUS).   

# Contributors

We kindly ask our contributors to include the automated tests of new functionality in this test suite. We wish 
you to make sure the existing and new tests pass locally before you open a pull-request. See our 
[pull-request template](../.github/pull_request_template.md) for more details.

To run the test suite locally:
- Clone the project to your local file system.
- Make sure you are in the top directory of the cloned project
- Make sure your python version is 3.7 or 3.8 `$ python --version`
- Make sure you have the latest version of pip `$ python -m pip install --upgrade pip`
```
    $ pip install pytest
    $ pip install pytest-cov
    $ pip install .
    $ python -m pytest tests --cov=kalmus --cov-config=.coveragerc --cov-report term-missing 
```# Description

Please include a summary of the change and which issue is fixed. Please also include relevant motivation and context. List any dependencies that are required for this change.

Fixes # (issue)

## Type of change

Please delete options that are not relevant.

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] This change requires a documentation update

# Checklist:

- [ ] My changes and comments follow the project's [code of conduct](../CODE_OF_CONDUCT.md)
- [ ] I have made corresponding changes to the documentation
- [ ] There exists a test or I have included a test in the [test suite](../tests) to test changes.
- [ ] My changes generate no new warnings
- [ ] New and existing unit tests pass locally with my changes
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
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Version [e.g. 22]

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
# Welcome to the Markdown User Guide for KALMUS (GUI)!

In this tutorial, I will introduce:  
1. **Installation of KALMUS package**
2. **What is KALMUS for**
    - Extract color information from film frames or brightness information from monochrome film frames using different color metrics and image sampling methods.   
    - Generate a barcode representation of the color/brightness information of a film. 3D information in 2D representation. 
    - Compare different barcodes globally through similarity measures on images. Interpret the difference through similarity scores.
    - Compare segments of barcodes locally using functions embedded in GUI. Interpret the difference using domain knowledge and contextual information extracted by KALMUS.  
3. **How to interact with KALMUS using its Graphic user interface**
    - Visualize barcodes
    - Generate barcodes
    - Change barcodes
    - Save barcodes
    - Load barcodes
    - Compare barcodes
    
## 1. Installation
There are two ways that you could install KALMUS on your local machine:  
1. (**Recommended**) Get the latest distribution of KALMUS from PyPI ([KALMUS Project Page on PyPI](https://pypi.org/project/kalmus/)).  
Use command `$ pip install kalmus` or `$ pip install --upgrade kalmus` (if kalmus has been installed) to install the latest version of the KALMUS package. All dependencies should be automatically installed during this process.

2. Alternatively, you could install the KALMUS locally by first cloning the GitHub repo of Kalmus ([GitHub page](https://github.com/KALMUS-Color-Toolkit/KALMUS)). Then, move to the top directory of cloned KALMUS project and install using the command `pip install .` 

**See our [Installation Guide](https://kalmus-color-toolkit.github.io/KALMUS/install.html) for more details.**

Once the package is installed, you could verify the version of KALMUS package using the command `$ pip show kalmus`  
<img src="notebook_figures/kalmus_version.png" alt="drawing" width="800 px"/>

## 2. What is KALMUS for?
KALMUS is a Python package for the computational analysis of colors in films. It addresses how to best describe a film's color. This package is optimized for two purposes: **(1) various ways to measure, calculate and compare a film's color and (2) various ways to visualize a film's color.**

KALMUS utilizes the movie barcode as a visualization of the film's color. It has a modularized pipeline for the generation of barcodes using different measures of color and region of interest in each film frame. It also provides a set of measures that allow users to compare different films' colors directly through this visualization.

### 2.1 Barcode Generation

Barcode supports __7 color metrics__ that measure the color of a frame and __5 frame types__ that specify which part of the frame will be used in the color measures.

Below is a table of available combinations of color metric and frame type in barcode generation.  

| frame_type \ color_metric | Average | Median |  Mode  | Top-dominant | Weighted-dominant | Brightest | Bright |
| --------------------------| :-----: | :----: | :----: | :----------: | :---------------: | :-------: | :----: |
| **Whole_frame**               | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |        &#9745;    |    &#9745;    |   &#9745;   |
| **High_contrast_region**      | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |      &#9745;      |    &#9745;    |   &#x2612;  |
| **Low_contrast_region**       | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |      &#9745;      |    &#9745;    |   &#x2612;  |
| **Foreground**                | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |      &#9745;      |    &#9745;    |   &#x2612;  |
| **Background**                | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |      &#9745;      |    &#9745;    |   &#x2612;  |

### 2.2 Examples of the color of a frame using a selected color metric and frame type.

Here, we show some example frames with their color extracted using the selected color metric and frame type

In the figures below,  
- On the left of each figure, we show the original frame (with letterboxing if applicable).
- On the right of each figure, we show the extracted region using the selected frame type with __the color of extracted region on the rightmost__.

**Casino Royale (2006) using Average Color with Whole frame or only Region with High (brightness) contrast**

![casino_whole_avg](notebook_figures/casino_2_whole_average.png)  

![casino_high_avg](notebook_figures/casino_2_high_average.png)

---

**Casino Royale (2006) using Average Color with Whole frame or only Foreground of frame**

![casino_2_whole_avg](notebook_figures/casino_1_whole_average.png)

![casino_2_fore_avg](notebook_figures/casino_1_fore_average.png)

---

**Incredibles (2004) using Whole frame with Mode color, Top-dominant color, or Brightest color**

![incre_whole_avg](notebook_figures/incredible_1_whole_mode.png)

![incre_whole_top](notebook_figures/incredible_1_whole_dominant.png)

![incre_whole_bri](notebook_figures/incredible_1_whole_brightest.png)

---

**Mission: Impossible (1996) using Whole frame and Foreground with Mode or Average color**

![mission_whole_mode](notebook_figures/mission_1_whole_mode.png)

![mission_fore_avg](notebook_figures/mission_1_fore_avg.png)

![mission_fore_mode](notebook_figures/mission_1_fore_mode.png)

---

**I, Robot (2004) using Median color with Whole, Foreground, or Background of frame**

![robot_whole_med](notebook_figures/robot_1_whole_median.png)

![robot_fore_med](notebook_figures/robot_1_fore_median.png)

![robot_back_med](notebook_figures/robot_1_back_median.png)

### 2.3 Examples of barcode generated from a whole film using selected color metric and frame type

Below, we show two barcodes generated from a whole film (Mission: Impossible (2006)) using two different frame types.

**Barcode generated using Average color and Whole_frame of each frame**  
![whole_barcode](notebook_figures/mission_barcode_whole_frame_avg.png)

**Barcode generated using Average color but only Foreground of each frame**  
![fore_barcode](notebook_figures/mission_barcode_Foreground_avg.png)

**Available options for comparing different barcode visualization**

We provide a set of six comparison metrics for users to assess the similarity between two barcodes.

| Comparison metric | Range |  Tag  |
| :---------------- | ----: | :---: |
| Normalized root mean square error | 0 least similar, 1 most similar | Image Similarity |
| Structural similarity index | 0 least similar, 1 most similar | Image Similarity |
| Cross correlation | -1 anti-similar, 1 most similar | Signal Correlation |
| Local cross correlation | -1 anti-similar, 1 most similar | Signal Correlation |
| Needleman-Wunsch | 0 least similar, 1 most similar | Sequence Matching |
| Smith-Waterman | 0 least similar, 1 most similar | Sequence Matching |

For more details, please see our paper [KALMUS: tools for color analysis of films](../paper/joss-paper.md)


## Get Started...

KALMUS has a low-level API, high-level command-line, and **Graphic user interface** for audiences from all backgrounds to take advantage of its functionality. 

In this notebook Guide, we will focus on the **Graphic user interface** of KALMUS.

## 3. How to interact with KALMUS through Graphic User Interface

If you have installed the KALMUS package on your machine with version 1.3.0 and onward, you can start the GUI using the command:

```
    $ kalmus-gui
```

Alternatively, you could import the main function of the GUI from `kalmus.command_line_gui` module.

```python
from kalmus.command_line_gui import main
main()
```

### 3.1 Main window of KALMUS

![kalmus_gui](notebook_figures/kalmus_gui_main_display.png)

- (1) The display 1 of Barcode (barcode image of Barcode Object)
- (2) The display 2 of Barcode (barcode image of Barcode Object)
- (3) A histogram plot of the [hue](https://en.wikipedia.org/wiki/HSL_and_HSV) (0 - 360 degree on the color wheel) distribution of the Barcode image in display 1.
- (4) A histogram plot of the [hue](https://en.wikipedia.org/wiki/HSL_and_HSV) distribution of the Barcode image in display 2.
- (5) Matplotlib's [interactive navigation toolbar](https://matplotlib.org/3.2.2/users/navigation_toolbar.html). Notice that we wish the users to use the **Save Image** button on the left instead of the save button on the toolbar if they only want to save the barcode image (not the whole figure).

The **display (1)** and **display(2)** are clickable plots.

- You can click on any point of the barcode image to get the RGB (Brightness for Brightness barcode) values, (x, y) position, frame index, and time of video at that point.
- You can also check the frames around that point **if you saved the frames** during the barcode generation (see section 3.2 (10) for how to save frames during the generation)

![clickable_plot](notebook_figures/kalmus_gui_main_2.png)

---

![main_2](notebook_figures/kalmus_gui_main_buttons.png)

---

### 3.2 (6) Generate Barcode Window

![gene](notebook_figures/kalmus_gui_generate_barcode.png)

- (1) Barcode Color/Brightness metric
- (2) Barcode Frame type
- (3) Barcode type (Color or Brightness)
- (4) Start collecting colors from frames at **Start at** (type: int) (**Optional**: No specified or specify start==0, no frames will be skipped)
- (5) Frame sampled rate: Collect color from one frame every **sampled rate** frame (type: int)
- (6) How many frames included in the generated barcode (type: int) (**Optional**: No specified or specify end. Collect color/brightness till the end of input video)
- (7) Alternatively, you could use the more intuitive time unit.

![gene2](notebook_figures/kalmus_gui_generate_barcode_2.png)

- (time unit) (4) Start at minutes:seconds of input video (minutes and seconds are all type: int) (**Optional**: No specified or specify start==0, no frames will be skipped)
- (time unit) (5) Period in seconds for one sampled frame (type: float)
- (time unit) (6) End at minutes:seconds of input video (minutes and seconds are all type: int) (**Optional**: No specified or specify end. Collect color/brightness till the end of input video)
- (8) The path to the input video. Users may use the Browse button to locate the media file directly.
- (9) Whether automatically detect the letterbox and remove. Recommend use **Auto**, use manual only if you know the exact location (in pixels) of the letterbox or the input video's letterboxing does not follow the convention (not black or in dark color).
- (10) Whether saved frames during the generation, and save one frame in how many seconds (seconds type: float).
- (11) Whether rescale the frames during the generation. Highly recommend resizing frames if you are using the frame type other than whole_frame or the input video is in high resolution.
- (12) Whether multi-threading the generation process. Highly recommend it if your processor supports multi-threading.
- (13) Start the generation of barcode

![specify](notebook_figures/kalmus_gui_generate_barcode_3.png)

- (14) Specify the meta information of the input video. **Warning:** Specify the input video's meta information before you press the generate barcode button! Press Update Meta Info to save the entered entries.

---

### 3.3 (7) Load Json Barcode Window

![load](notebook_figures/kalmus_gui_load_json.png)

- Specify the file path to the .JSON Barcode (what is a JSON Barcode? check section 3.6 below)
- Specify the type of barcode saved in JSON
- Specify which barcode display on the Main window that you will load the barcode into
- Press the Load button to load the JSON barcode

---

### 3.4 (8) Load Barcode from Memory Window

![load_mem](notebook_figures/kalmus_gui_load_memory.png)

- Every barcode generated from the current running GUI or Loaded from JSON barcode will be stored on the memory
- User can load them onto the main display by selecting the name of barcode on the list
- Specify which display on the main window that new barcode will be loaded into
- Press the Load Selected Barcode button

---

### 3.5 (9) Reshape Barcode Window

![reshape](notebook_figures/kalmus_gui_reshape_barcode.png)

**There are three options available for users to change the barcode on the display**

- Reshape how many frames==pixels in each column of frames (similar to numpy.ndarray.reshape)
- Scale the barcode image by enlarging or shrinking the barcode image by a factor
- Resize the barcode image to a specific size in pixels

In the window:  
- (1) Show the current spatial size of the selected barcode in the main display (Barcode 1 in this case)
- (2) Select which options to use
- (3) Select which Barcode to change
- Press Process to change the Barcode using the given option and parameters

---

### 3.6 (10) Save JSON Barcode Window

![save](notebook_figures/kalmus_gui_save_json.png)

Similar to the load memory window
- Select the barcode on memory (list) that you wish to save locally as a JSON file
- Give the path to the saved JSON file in JSON file path textbox
- Press the Save Barcode button

The attributes of Barcode Object will be stored in a JSON file that can be used to rebuild the Barcode Object (in GUI, you simply reload the JSON barcode through Load JSON Window **section 3.3**)

---

### 3.7 (11) Save Barcode Image Window

![save](notebook_figures/kalmus_gui_save_image.png)

- Select which barcode on the main display that you wish to save locally as an image.
- The Saved width and height are automatically filled with the current width and height of barcodes. You could change to your desirable spatial size.
- Specify the path to the saved image file in the Image file path textbox
- Press the Save Barcode button

---

### 3.8 (12) Inspect Barcode Window

![inspect](notebook_figures/kalmus_gui_inspect.png)

You will first be asked which barcode on the main display that you wish to inspect in further details.

![inspect](notebook_figures/kalmus_gui_inspect_2.png)

In the inspect window there are three options to explore

- (1) Output the color/brightness data of the Color/Brightness barcode into a csv file
- (2) Show the histogram distribution of hue values of the Color barcode or brightness value of Brightness barcode (similar to those in the main display)
- (3) (Only available for Color barcode) Show the distribution of RGB color of the Color barcode in RGB cube.

![cube](notebook_figures/kalmus_gui_inspect_3.png)

---

### 3.9 (13) Statistics Information Window

![stats](notebook_figures/kalmus_gui_stats.png)

The similarity comparison between the displayed barcodes using a set of six comparison metrics.

**Warning:** The initiation of this window may take tens of seconds.

For more references about these six comparison metrics, please check section 2.3 above.

---

### 3.10 (14) Check Meta Information Window

Similarly to the **Inspect Barcode Window**

![inspect](notebook_figures/kalmus_gui_inspect.png)

You will first be asked which barcode on the main display that you wish to check for meta information.

![meta](notebook_figures/kalmus_gui_check_meta.png)

- A list of meta-information will be shown here
- To update the meta information, similarly to Specify Meta Info in the barcode generation, use the Update Meta Info button

![specify](notebook_figures/kalmus_gui_generate_barcode_3.png)

- Hit the Update Meta Info button in the Specify Meta Data window after you update the entries.
- Hit Refresh in Check Barcode Meta Information Window to see the updates
- To reflect the updates on the title of plots in the main display, find the barcode with updated meta information in the memory using the Load Memory button and load the updated barcode back to the main display.

---

### 3.11 (15) Quit

Quit the KALMUS's GUI. **Be sure to save all the barcodes you like before you quit the program, and make sure the Generate Barcode Window is closed before quitting**.

---

## 4. Thank you!

Thank you so much for reading through this markdown tutorial! If you find any errors in the instructions, please feel free to email the tutorial author, Yida Chen, <yc015@bucknell.edu>
[![Project Status](https://img.shields.io/pypi/status/kalmus.svg)](https://pypi.org/project/kalmus/)
[![Python Version](https://img.shields.io/pypi/pyversions/kalmus.svg)](https://pypi.org/project/kalmus/)
[![PyPI Version](https://img.shields.io/pypi/v/kalmus.svg)](https://pypi.org/project/kalmus/)

# Notebooks Tutorial
Welcome to the KALMUS's notebook tutorials!

This folder contains the IPython notebooks that walk you through how to interact with the KALMUS package.

To run these notebooks, please install the latest version of Jupyter Lab 
([Installation Instruction for JupyterLab](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html)) 
and the latest version of KALMUS package ([Installation Guide for KALMUS](https://kalmus-color-toolkit.github.io/KALMUS/install.html)).
 
When running the notebooks, make sure you have downloaded the whole [notebooks folder](https://github.com/KALMUS-Color-Toolkit/KALMUS/archive/refs/heads/master.zip) or cloned the project on your local machine. Then, start the Jupyter Lab inside the notebooks folder using command:

```
$ jupyter lab
```

**Content:**  
- [Notebook Figures](notebook_figures)
- [Notebook Example Data](notebook_example_data): You are welcome to put your own data inside this folder when walking 
through the examples in notebooks!
- [Notebook Utility](notebook_utils.py): Plot utility for you to visualize certain processes occurring during the 
barcode generation.
- Notebook Guide for how to interact with [Graphic User Interface](user_guide_for_kalmus_gui.ipynb)
- Notebook Guide for how to interact with [Application Programming Interface](user_guide_for_kalmus_api.ipynb)
- Notebook Guide for [Advanced usage](advanced_guide_for_kalmus_api.ipynb) of Application Programming Interface (how to create your own Barcode Visualization)
- Markdown Guide for how to interact with [Command-line Interface](USAGE_COMMAND_LINE_UI.md)

# For general audiences
We provide the markdown version of user guides for kalmus' Command-line interface and Graphic user interface. The rendered markdown tutorials can be read directly on the GitHub.

**Content:**
- Markdown Guide for how to interact with [Graphic User Interface](USAGE_GRAPHIC_USER_INTERFACE.md)
- Markdown Guide for how to interact with [Command-line Interface](USAGE_COMMAND_LINE_UI.md)

# Welcome to the User Guide for KALMUS (Command-line Interface)!

In this Markdown Tutorial I will introduce:
- **Installation of KALMUS package**
- **How to use the Command-line interface to automate your barcode generation workflow**
- **What are the available argument options for command-line generator**
    - -p --path
    - --color_metric
    - --frame_type
    - --barcode_type
    - --skip
    - -s --step
    - -t --total_frames
    - --num_thread
    - --save_frame_rate
    - --rescale_frame_factor
    - -o --output_path

# Notice

The Command-line Interface of KALMUS is used specifically for **generating barcode** and **saving** it into 
reloadable (via both [API](user_guide_for_kalmus_api.ipynb) and [GUI](user_guide_for_kalmus_gui.ipynb)) JSON objects.

We recommend our users first walk through the IPython notebook Tutorials on either 
[KALMUS's Graphic user interface](user_guide_for_kalmus_gui.ipynb) or [KALMUS's Application programming interface](user_guide_for_kalmus_api.ipynb), 
in which the full functionality of KALMUS is covered.

Once you have been familiar with the API, GUI, or both interfaces, you could further automate your barcode generation workflow 
through this **Command-line interface**.

# Installation Guide

There are two ways that you could install KALMUS on your local machine:  
1. (**Recommended**) Get the latest distribution of KALMUS from PyPI ([KALMUS Project Page on PyPI](https://pypi.org/project/kalmus/)).  
Use command `$ pip install kalmus` or `$ pip install --upgrade kalmus` (if kalmus has been installed) to install the latest version of the KALMUS package. All dependencies should be automatically installed during this process.

2. Alternatively, you could install the KALMUS locally by first cloning the GitHub repo of Kalmus ([GitHub page](https://github.com/KALMUS-Color-Toolkit/KALMUS)). Then, move to the top directory of cloned KALMUS project and install using the command `pip install .` 

**See our [Installation Guide](https://kalmus-color-toolkit.github.io/KALMUS/install.html) for more details.**

Once the package is installed, you could verify the version of the KALMUS package using the command `$ pip show kalmus`  
![KALMUS version](notebook_figures/kalmus_version.png)

Alternatively, in version 1.3.7 and above, you can check the version of installed kalmus using its 
`.__version__` attribute.

```jupyter
>>> import kalmus
>>> print(kalmus.__version__) # Warning: The __version__ attribute is not available in the kalmus v.1.3.6 and backward
>>> 1.3.7 
```

## Important!

The **Command-line interface** is a new feature added into the KALMUS in its 1.3.7 version. To use the Command-line 
interface, you have to make sure the version of installed KALMUS is 1.3.7 or onward.

# How do I use the command-line interface?

The command-line interface is similar to the **BarcodeGenerator** object, which we have covered in [API Guide](user_guide_for_kalmus_api.ipynb). 
You could invoke the command-line generator using command `$ kalmus-generator`. 

Notice that if you have installed 
kalmus>=1.3.7 but failed to start the kalmus-generator through this command, you could find the corresponding executable 
in your_python_directory/Scripts/ directory on Windows or .local/bin/ directory on Linux.

Check the available arguments for the `kalmus-generator` command using the flag `-h`:

```
Windows OS
$ kalmus-generator -h
$ usage: kalmus-generator [-h] -p PATH --color_metric COLOR_METRIC --frame_type
                          FRAME_TYPE --barcode_type BARCODE_TYPE [--skip SKIP]
                          -s STEP [-t TOTAL_FRAMES] [--num_thread NUM_THREAD]
                          [--saved_frame_rate SAVED_FRAME_RATE]
                          [--rescale_frame_factor RESCALE_FRAME_FACTOR]
                          [-o OUTPUT_PATH]
         Command line Barcode generator
         ......
```

# Let's take a quick look at what each argument means!

Only the arguments marked with **(required)** must be specified in the kalmus-generator command.


| Arguments | Description |   Type   |
| :------- | :----------- | :------: |
| **-p --path (required)**: | The relative or absolute path to the input media/video file. Equivalent to the **video_file_path** in `BarcodeGenerator.generate_barcode` | `str` |
| **--color_metric (required)**:| The color_metric selected for barcode generation. Equivalent to the **color_metric** in `BarcodeGenerator` | `str` |
| **--frame_type (required)**:| The frame_type selected for barcode generation. Equivalent to the **frame_type** in `BarcodeGenerator` | `str` |
| **--barcode_type (required)**:| The barcode_type of the generated barcode. Equivalent to the **barcode_type** in `BarcodeGenerator` | `str` |
| **--skip**:| The number of frames to be skipped at the start of the input video before collecting color/brightness. Equivalent to the **skip_over** in `BarcodeGenerator` | `int>=0` |    
| **-s --step (required)**:| The frame sampled rate. Collect color from one frame every **step** frames. Equivalent to **sampled_frame_rate** in `BarcodeGenerator`| `int>=1` |
| **-t --total_frames**:| The total number of frames to be included in generated barcode. Equivalent to **total_frames** in `BarcodeGenerator`. **Notice**: If you wish to generate a barcode till the end of video, simply put a very large number in  total_frames, e.g. total_frames = 1e8. The barcode will auto adjust the total frames using film length (in frames) and your specified skip_over and sampled_frame_rate to recompute the correct number for total_frames, and the barcode generator collects color/brightness till the last frame of the input video.| `int>=0` |
| **--num_thread**:| Number of threads to use in barcode generation. Equivalent to **num_thread** in `BarcodeGenerator.generate_barcode`. We highly recommend this if your processor supports multi-threading. | `int>=2` |
| **--saved_frame_rate**:| The rate of saving frames (thumbnail quality) in barcode generation (save rate's unit: seconds). Equivalent to **saved_frame_rate** in `BarcodeGenerator.generate_barcode`. The saved frames can be very useful when visualizing barcode in [GUI](user_guide_for_kalmus_gui.ipynb) as you may correlate or validate a segment of colors/brightness barcode with its corresponding frames. **However**, since the size of saved frames/images grows very quickly, you may wish to set this saved rate (seconds) to be low or not to use this option (by default). | `float>0` |
| **--rescale_frame_factor**:| The factor of rescaling the frames when collecting color. Equivalent to **rescale_frame_factor** in `BarcodeGenerator.generate_barcode`. resize width = sqrt(rescale_frame_factor) x original width, resize height = sqrt(rescale_frame_factor) * original height. We recommend you to use this option speed up the generation process if your input video's resolution is above 2K (or 1K for Top-dominant, Weighted-dominant, Bright color metric and Low/High_contrast_region or Foreground/Background frame type). | `1>float>0` |
| **-o --output_path**:| The output path to the saved JSON object. Equivalent to **filename** in `Barcode.save_as_json`. By default, the generated JSON file will be stored on the current directory with filename *saved_{barcode_type}_barcode_{frame_type}_{color_metric}.json* | `str` |

**The available combinations of frame_type and color_metric are the same as those in GUI and API.**

| frame_type \ color_metric | Average | Median |  Mode  | Top-dominant | Weighted-dominant | Brightest | Bright |
| --------------------------| :-----: | :----: | :----: | :----------: | :---------------: | :-------: | :----: |
| **Whole_frame**               | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |        &#9745;    |    &#9745;    |   &#9745;   |
| **High_contrast_region**      | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |      &#9745;      |    &#9745;    |   &#x2612;  |
| **Low_contrast_region**       | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |      &#9745;      |    &#9745;    |   &#x2612;  |
| **Foreground**                | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |      &#9745;      |    &#9745;    |   &#x2612;  |
| **Background**                | &#9745; |   &#9745;  |  &#9745; |      &#9745;     |      &#9745;      |    &#9745;    |   &#x2612;  |

# Example Commands

(1) Generate a **Color** Barcode using **Whole_frame** with **Average** color from the *i_robot_video.mp4* in [notebook_example_data](notebook_example_data). 
**Skip** the first 10 frames of input video, **sample** every frame, and include **100 total_frames.**

```
$ kalmus-generator -p notebook_example_data/i_robot_video.mp4 --frame_type Whole_frame --color_metric Average --skip 10 --step 1 --total_frames 100 --barcode_type Color
```

A JSON file with filename: **saved_Color_barcode_Whole_frame_Average.json** should be generated in this folder after the command's execution is finished.

(2) Generate a **Brightness** Barcode using **High_contrast_region** with **Median** color from the *i_robot_video.mp4* in [notebook_example_data](notebook_example_data). 
**Not skipping** any frames, **sampled** every 2 frames, and **include all frames till the end**. **Save frames** every 
0.5 seconds, using **2 threads**, and **rescale the frame** by 0.25 during the generation. **Saved the generated barcode** 
in *./saved_barcode_cli.json*

```
$ kalmus-generator -p notebook_example_data/i_robot_video.mp4 --frame_type Low_contrast_region --color_metric Median --step 2 --total_frames 1000000 --barcode_type Brightness --num_thread 2 --saved_frame_rate 0.5 --rescale_frame_factor 0.25 -o ./saved_barcode_cli.json
```

A JSON file with filename: **saved_barcode_cli.json** should be generated in this folder after the 
command's execution finished.

# Thank you!

Thank you so much for reading through this markdown Tutorial for KALMUS's Command-line interface.
If you find any errors in the instructions, please feel free to email the markdown author, Yida Chen, <yc015@bucknell.edu> 
