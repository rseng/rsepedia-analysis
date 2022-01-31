# PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results of AutoDock Vina and PLANTS

[![GitHub Actions Build Status](https://github.com/radifar/PyPLIF-HIPPOS/workflows/CI/badge.svg)](https://github.com/radifar/PyPLIF-HIPPOS/actions?query=workflow%3ACI)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/radifar/PyPLIF-HIPPOS.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/radifar/PyPLIF-HIPPOS/context:python)
[![codecov](https://codecov.io/gh/radifar/pyplif-hippos/branch/main/graph/badge.svg)](https://codecov.io/gh/radifar/pyplif-hippos/branch/main)   
[![Anaconda-Server Badge](https://img.shields.io/badge/Install%20with-conda-green.svg?style=flat)](https://anaconda.org/conda-forge/pyplif-hippos)
[![Documentation Status](https://readthedocs.org/projects/pyplif-hippos/badge/?version=latest&style=flat)](https://pyplif-hippos.readthedocs.io/en/latest/)
[![DOI:10.1021/acs.jcim.0c00305](https://zenodo.org/badge/DOI/10.1021/acs.jcim.0c00305.svg)](https://doi.org/10.1021/acs.jcim.0c00305)   
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/conda-forge/pyplif-hippos?color=green)](https://anaconda.org/conda-forge/pyplif-hippos)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fradifar%2Fpyplif-hippos&title=visitor%20today%2Ftotal)](https://hits.seeyoufarm.com)

<p align="center">
  <img alt="Icons made by Freepik from Flaticon is licensed by CC 3.0 BY" src="docs/source/hippopotamus_small.png">
</p>

<p align="center">Icons made by <a href="https://www.freepik.com/">Freepik</a> from <a href="http://www.flaticon.com">Flaticon</a> is licensed by CC 3.0 BY</p>

Welcome to PyPLIF-HIPPOS's project page. PyPLIF-HIPPOS is an upgraded version of [PyPLIF](https://github.com/radifar/pyplif/) (**Python-based Protein-Ligand Interaction Fingerprinting**), a tool for molecular docking post-analysis. It will translate the 3D coordinates of both ligand(s) (generated from docking simulation) and protein into a series of *interaction bitstring* (also known as *Interaction Fingerprint*) (see image below). **HIPPOS** (/ˌhipoʊz/) is a recursive acronym of **HIPPOS Is PyPLIF On Steroids**. From this point forward, PyPLIF-HIPPOS is simplified to HIPPOS.

Compared to PyPLIF, HIPPOS is not only faster and able to generate more customized interaction bitstring, but also supports both [PLANTS](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/pharmazie-und-biochemie/pharmazie/pharmazeutische-chemie/pd-dr-t-exner/research/plants/) & [Vina](http://vina.scripps.edu/)! More over, unlike its predecessor it is (far) more well-documented.

<p align="center">
  <img alt="Table of Content Abstract Graphic JCIM" src="docs/source/toc-abstract-graphics_small.png">
</p>

<p align="center">Reprinted with permission from https://doi.org/10.1021/acs.jcim.0c00305. Copyright 2020 American Chemical Society.</p>

<p align="center">
  <img alt="PyPLIF output from PyPLIF publication" src="docs/source/pyplif.png">
</p>

<p align="center">Illustration by Radifar et al (2013) from <a href="http://www.bioinformation.net/009/97320630009325.htm">Bioinformation.net</a> is licensed by <a href="http://creativecommons.org/licenses/by/4.0">CC 4.0 BY</a>

## Quick Installation

The easiest way to install HIPPOS is using [Anaconda or Miniconda](https://docs.anaconda.com/anaconda/install/).
If you have Anaconda or Miniconda ready in your machine, you can start with
creating new environment (recommended):

`conda create -n hippos python=3.6`

Then activate the environment and install HIPPOS:

`conda activate hippos`  
`conda install -c conda-forge pyplif-hippos`

next you can try run HIPPOS and HIPPOS-GENREF with the following command:

`hippos`  
`hippos-genref` 

## How to Use HIPPOS

So I already installed HIPPOS, now what? Well you could start with how to generate
the [reference bitstring](https://pyplif-hippos.readthedocs.io/en/latest/getting-started-genref.html)
and Getting Started tutorial for [AutoDock Vina](https://pyplif-hippos.readthedocs.io/en/latest/getting-started-vina.html)
or [PLANTS](https://pyplif-hippos.readthedocs.io/en/latest/getting-started-plants.html).

## Ideas for Improvement? Found Bug(s)?

If you have any idea for improvement or found bug to report feel free to write them [here](https://github.com/radifar/PyPLIF-HIPPOS/issues).

## Citing HIPPOS

If you are using HIPPOS please cite this paper:

Istyastono, E., Radifar, M., Yuniarti, N., Prasasty, V. and Mungkasi, S., 2020. 
PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results 
of AutoDock Vina and PLANTS. Journal of Chemical Information and Modeling, 60(8), pp.3697-3702.
https://doi.org/10.1021/acs.jcim.0c00305

## Acknowledgment

This project has received funding from the [National Agency for Research and Innovation](https://international.ristekdikti.go.id/) (Indonesia)
under grant agreement No. 807.7/LL5/PG/2020. This project has been restructured based on the
[MOLSSI Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
version 1.3, and benefited greatly from [MOLSSI Python Package Development Best Practices](https://molssi.org/2020/04/20/may-webinar-series-python-package-development/)
workshop.


-----

&copy; Copyright 2021, Muhammad Radifar & Enade Perdana Istyastono
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age,
body size, disability, ethnicity, gender identity and expression, level of
experience, nationality, personal appearance, race, religion, or sexual
identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Moreover, project maintainers will strive to offer feedback and advice to
ensure quality and consistency of contributions to the code.  Contributions
from outside the group of project maintainers are strongly welcomed but the
final decision as to whether commits are merged into the codebase rests with
the team of project maintainers.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an
appointed representative at an online or offline event. Representation of a
project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at 'muhammad.radifar@picomps.org'. The project team will
review and investigate all complaints, and will respond in a way that it deems
appropriate to the circumstances. The project team is obligated to maintain
confidentiality with regard to the reporter of an incident. Further details of
specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.4, available at
[http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Sample Package Data

This directory contains sample additional data you may want to include with your package.
This is a place where non-code related additional information (such as data files, molecular structures,  etc.) can 
go that you want to ship alongside your code.

Please note that it is not recommended to place large files in your git directory. If your project requires files larger
than a few megabytes in size it is recommended to host these files elsewhere. This is especially true for binary files
as the `git` structure is unable to correctly take updates to these files and will store a complete copy of every version
in your `git` history which can quickly add up. As a note most `git` hosting services like GitHub have a 1 GB per repository
cap.

## Including package data

Modify your package's `setup.py` file and the `setup()` command. Include the 
[`package_data`](http://setuptools.readthedocs.io/en/latest/setuptools.html#basic-use) keyword and point it at the 
correct files.

## Manifest

* `look_and_say.dat`: first entries of the "Look and Say" integer series, sequence [A005150](https://oeis.org/A005150)
# Development, testing, and deployment tools

This directory contains a collection of tools for running Continuous Integration (CI) tests, 
conda installation, and other development tools not directly related to the coding process.


## Manifest

### Continuous Integration

You should test your code, but do not feel compelled to use these specific programs. You also may not need Unix and 
Windows testing if you only plan to deploy on specific platforms. These are just to help you get started

* `travis-ci`: Linux and OSX based testing through [Travis-CI](https://about.travis-ci.com/) 
  * `before_install.sh`: Pip/Miniconda pre-package installation script for Travis 
* `appveyor`: Windows based testing through [AppVeyor](https://www.appveyor.com/) (there are no files directly related to this)

### Conda Environment:

This directory contains the files to setup the Conda environment for testing purposes

* `conda-envs`: directory containing the YAML file(s) which fully describe Conda Environments, their dependencies, and those dependency provenance's
  * `test_env.yaml`: Simple test environment file with base dependencies. Channels are not specified here and therefore respect global Conda configuration
  
### Additional Scripts:

This directory contains OS agnostic helper scripts which don't fall in any of the previous categories
* `scripts`
  * `create_conda_env.py`: Helper program for spinning up new conda environments based on a starter file with Python Version and Env. Name command-line options


## How to contribute changes
- Clone the repository if you have write access to the main repo, fork the repository if you are a collaborator.
- Make a new branch with `git checkout -b {your branch name}`
- Make changes and test your code
- Ensure that the test environment dependencies (`conda-envs`) line up with the build and deploy dependencies (`conda-recipe/meta.yaml`)
- Push the branch to the repo (either the main or your fork) with `git push -u origin {your branch name}`
  * Note that `origin` is the default name assigned to the remote, yours may be different
- Make a PR on GitHub with your changes
- We'll review the changes and get your code into the repo after lively discussion!


## Checklist for updates
- [ ] Make sure there is an/are issue(s) opened for your specific update
- [ ] Create the PR, referencing the issue
- [ ] Debug the PR as needed until tests pass
- [ ] Tag the final, debugged version 
   *  `git tag -a X.Y.Z [latest pushed commit] && git push --follow-tags`
- [ ] Get the PR merged in

## Versioneer Auto-version
[Versioneer](https://github.com/warner/python-versioneer) will automatically infer what version 
is installed by looking at the `git` tags and how many commits ahead this version is. The format follows 
[PEP 440](https://www.python.org/dev/peps/pep-0440/) and has the regular expression of:
```regexp
\d+.\d+.\d+(?\+\d+-[a-z0-9]+)
```
If the version of this commit is the same as a `git` tag, the installed version is the same as the tag, 
e.g. `pyplif_hippos-0.1.2`, otherwise it will be appended with `+X` where `X` is the number of commits 
ahead from the last tag, and then `-YYYYYY` where the `Y`'s are replaced with the `git` commit hash.
# Sample Package Data

This directory contain the data required for testing. As integrating the content
into test function became impractical (e.g. inserting protein structure in Python
script). Also separating the test data with the script make the test code cleaner.
Parameter Options
=================

The parameters to identify the interactions refer to those 
used in `PyPLIF <https://doi.org/10.6026/97320630009325>`_, which 
inspired by the IFP of `Marcou and Rognan <https://doi.org/10.1021/ci600342e>`_. 
These parameters can be modified in ``.hippos/PARAMETERS.py`` in your home
directory after HIPPOS installed.

In general there are two rules that can be modified:

1. Maximum distance value.
2. Interaction angle limit (for hydrogen bond and aromatic interaction).

Here is the content of ``.hippos/PARAMETERS.py`` ::

    '''
        Parameter for interaction distance
        Interaction exist if lower or equal 
        than these values
    '''

    HYDROPHOBIC		= 4.5 
    AROMATIC		= 4.0
    HBOND		= 3.5
    ELECTROSTATIC	= 4.0

    '''
        Parameter for minimum H bond angle.
        Interaction exist if O --- H-D angle
        higher or equal than HBOND_ANGLE value.
    '''

    HBOND_ANGLE     = 135

    '''
        Parameter for aromatic interaction angle.
            Face to Face if:
        AROMATIC_ANGLE_LOW >= Angle between aromatic plane
        Or AROMATIC_ANGLE_HIGH <= Angle between aromatic plane
            Edge to Face if:
        AROMATIC_ANGLE_LOW < Angle between aromatic plane < 
        AROMATIC_ANGLE_HIGH
    '''

    AROMATIC_ANGLE_LOW  = 30.0
    AROMATIC_ANGLE_HIGH = 150.0


The explanation of each value can be seen from the comment above the
value assignment.

.. _simplified-rule:

Simplified Bitstring Rule
-------------------------

HIPPOS provides simplified interaction bitstring to reduce unnecessary
string in the output. The bitstring simplification rule can be seen in
the table below. Number 1-7 in the header of the table correspond to
the interaction type described in PyPLIF

.. code-block::

    1 Correspond to Hydrophobic
    2 Correspond to Aromatic Face to Face
    3 Correspond to Aromatic Edge to Face
    4 Correspond to H-bond (protein is donor)
    5 Correspond to H-bond (protein is acceptor)
    6 Correspond to Electrostatic (protein +)
    7 Correspond to Electrostatic (protein -)

.. raw:: html
    
    <embed>
    <style type="text/css">
    .tg  {border:none;border-collapse:collapse;border-color:#9ABAD9;border-spacing:0;}
    .tg td{background-color:#EBF5FF;border-color:#9ABAD9;border-style:solid;border-width:0px;color:#444;
      font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:5px 15px;word-break:normal;}
    .tg th{background-color:#409cff;border-color:#9ABAD9;border-style:solid;border-width:0px;color:#fff;
      font-family:Arial, sans-serif;font-size:14px;font-weight:bold;overflow:hidden;padding:5px 15px;word-break:normal;}
    .tg .tg-h1fx{background-color:#D2E4FC;border-color:inherit;font-family:"Courier New", Courier, monospace !important;;
      text-align:center;vertical-align:top}
    .tg .tg-juju{font-family:"Courier New", Courier, monospace !important;;text-align:left;vertical-align:top}
    .tg .tg-3ib7{border-color:inherit;font-family:"Courier New", Courier, monospace !important;;text-align:center;vertical-align:top}
    .tg .tg-64ye{background-color:#D2E4FC;font-family:"Courier New", Courier, monospace !important;;text-align:left;vertical-align:top}
    </style>
    <table class="tg">
    <thead>
      <tr>
        <th class="tg-3ib7">AA</th>
        <th class="tg-3ib7">1</th>
        <th class="tg-3ib7">2</th>
        <th class="tg-juju">3</th>
        <th class="tg-3ib7">4</th>
        <th class="tg-3ib7">5</th>
        <th class="tg-3ib7">6</th>
        <th class="tg-3ib7">7</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td class="tg-h1fx">ALA</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">CYS</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;<br></td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">ASP</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">&#10003;</td>
      </tr>
      <tr>
        <td class="tg-3ib7">GLU</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">&#10003;</td>
      </tr>
      <tr>
        <td class="tg-h1fx">PHE</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-64ye">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">GLY</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">HIS</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-64ye">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">ILE</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">LYS</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">LEU</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">MET</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">ASN</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;<br></td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">PRO</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">GLN</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">ARG</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">&#10003;<br></td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">SER</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-<br></td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">THR</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">VAL</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">TRP</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">TYR</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
    </tbody>
    </table>
    </br>
    </embed>HIPPOS Configuration Options
============================

There are two kinds of options, *essential* and *optional*. When nothing declared it means the option is essential. Some options have default value. Also keep in mind that **all file names must never use space!** So use underscore instead.

You can also inserting comments at the beginning of the line or after the option-value pair by inserting the ``#`` sign before the comment. Everything inserted after the ``#`` sign will be ignored by the software.

Basic Options
-------------

* **docking_method**
	value: ``vina`` or ``plants``

* **docking_conf**
	value: ``docking_configuration_file_name`` , eg. ``vina.conf`` or ``plants.conf``

Input Options
-------------

* **residue_name**
	value: ``list of residue_name`` , eg. ``ASP107 SER111 THR112``
	
	The list of residue_name, each residue separated by space. It is used in PLANTS post-analysis but not in VINA analysis as pdbqt don't hold the residue_name-residue_number pair field. However it is highly recommended to define it in VINA post analysis as it will be included in output file, making the results easier to interpret.
	
* **residue_number**
	value: ``list of residue_number`` , eg. ``80 84 85``
	
	The list of residue number, each number separated by space. Essential option in VINA post-analysis but optional in PLANTS post-analysis.
	
* **similarity_coef** (optional)
	value: ``list of similarity_coef`` , eg. ``tanimoto`` or ``mcconnaughey`` or ``tanimoto mcconnaughey``
	
	If the value is set then **full_ref** or **full_nobb_ref** or **simplified_ref** (depends on output_mode used)
	must be provided, so interaction fingerprint can be compared against reference.

	* **full_ref**
		value: ``bitstring1 bitstring2 ... n``
		
		At least 1 uniform bitstring must be provided

	* **full_nobb_ref**
		value: ``bitstring1 bitstring2 ... n``
		
		At least 1 uniform bitstring must be provided

	* **simplified_ref**
		value: ``bitstring1 bitstring2 ... n``
		
		At least 1 simplified bitstring must be provided

Output Options
--------------

* **full_outfile** (optional)
	value: ``full_output_file_name`` , default: ``full_ifp.csv``
	
	Only used by ``full`` output_mode. It is recommended to use the csv extensions for clarity.

* **full_nobb_outfile** (optional)
	value: ``full_output_file_name`` , default: ``full_nobb_ifp.csv``
	
	Only used by ``full_nobb`` output_mode. It is recommended to use the csv extensions for clarity.

* **simplified_outfile** (optional)
	value: ``simplified_output_file_name`` , default: ``simplified_ifp.csv``
	
	Only used by ``simplified`` output_mode. It is recommended to use the csv extensions for clarity.

* **sim_outfile** (optional)
	value: ``similarity_output_file_name``, default: ``similarity.csv``
	
	Only used when the similarity coefficient is calculated.

* **logfile** (optional)
	value: ``log_file_name`` , default: ``hippos.log``
	
.. _advanced-options:
	
Advanced Options
----------------

* **output_mode** (optional)
	value: ``list of output_mode`` , default: ``full``, options: ``full``, ``full_nobb``, and ``simplified``
	
	The list of output_mode, at least one value required. When multiple value provided each value must be 
	separated by space. Multiple value can only be used in HIPPOS without reference. When nothing provided
	output_mode full will be used.
	
* **docking_score** (optional)
	value: ``yes`` or ``no`` , default: ``yes``
	
	Extract the docking score of ligand poses from docking results and attach them to output file.

* **omit_interaction** (optional)
	value: ``interaction_type`` and ``residue_name``

	where ``interaction_type`` is one of the following value:

	- ``hydrophobic`` or ``HPB``
	- ``aromatic`` or ``ARM``
	- ``h_bond`` or ``HBD``
	- ``electrostatic`` or ``ELE``
	- ``h_bond_donor`` or ``HBD_DON``
	- ``h_bond_acceptor`` or ``HBD_ACC``
	- ``electrostatic_positive`` or ``ELE_POS``
	- ``electrostatic_negative`` or ``ELE_NEG``
	- ``aromatic_facetoface`` or ``ARM_F2F``
	- ``aromatic_edgetoface`` or ``ARM_E2F``

	While ``residue_name`` specify which residue will be omitted. Usage example:

	``omit_interaction hydrophobic ARG223``

..
	* **res_weight1** (optional)
		value: ``residue_number residue_name interaction_type weight`` , eg. ``80 ASP107 electrostatic 5``
		
		Give weight to a spesific interaction on spesific amino acid residue. The example above shows that the number of electrostatic interaction bit on ASP107 will be multiplied by 5. The number 1 in **res_weight1** can be replaced with 2-5, Therefore there are 5 weight that can be applied to interaction fingerprinting.
Advanced Usage
==============

Customize Bitstring by Omitting Specific Interactions
-----------------------------------------------------

It is possible to omit certain interaction(s) in selected residue(s).
This will cause PyPLIF HIPPOS to replace the omitted interaction in the bitstring
with ``n`` character (eg. from ``1010001`` to ``1n10001``). As a consequence, the omitted
interaction will not be included in similarity calculation. This is done by
removing the omitted interaction bit both in the target and the reference bitstring. As
a result the similarity coefficient could be different from the default one.

To omit interaction just add the following line to the configuration file:

``omit_interaction interaction_type residue_name [single or multiple]``

where ``interaction_type`` is one of the following value:

- ``hydrophobic`` or ``HPB``
- ``aromatic`` or ``ARM``
- ``h_bond`` or ``HBD``
- ``electrostatic`` or ``ELE``
- ``h_bond_donor`` or ``HBD_DON``
- ``h_bond_acceptor`` or ``HBD_ACC``
- ``electrostatic_positive`` or ``ELE_POS``
- ``electrostatic_negative`` or ``ELE_NEG``
- ``aromatic_facetoface`` or ``ARM_F2F``
- ``aromatic_edgetoface`` or ``ARM_E2F``

As for the ``residue_name`` it could be one or more residue name. Here is
one example of a valid omit_interaction definition:

``omit_interaction hydrophobic ARG223``

As a quick start here is three different scenario on omitting interaction.

Omit single interaction
^^^^^^^^^^^^^^^^^^^^^^^

This configuration file (`examples-input/08-na_omit_interaction/vina-omit-simple.txt`)
shows how you can omit hydrophobic interaction on ARG223 residue ::

    docking_method    vina    # plants or vina
    docking_conf      ../03-na_vina/vina-003.conf

    similarity_coef   tanimoto mcconnaughey

    full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000  00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000  00010101000000100000000000000000000100000000000001010000000000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

    omit_interaction  hydrophobic  ARG223

    residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
    residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

    full_outfile omit-simple.csv
    sim_outfile omit-simple-similarity.csv
    logfile omit-simple.log

Here is the excerpt from the bitstring output using the above configuration file ::

    129821_1         -6.9      0000000100000000000000000000100000010010000000000000000000000001000000n000000000000010000001000000000000000000000000000000100000000010000000000100000000000001000
    129821_2         -6.8      0000000000000000000000000000000000010000000000000101000000000001000000n000000000000000000001000000000000000001010000000000100000000010000000000000000000000001000
    129821_3         -6.6      0000000000000000000000000000100000110000000000000101000000000001000000n000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000
    129821_4         -6.4      0000000000000000000000000000100010000000000000000000000000000001000000n000000000000000000001000000000000000000010000101000000000000000000000000000000000000000000
    129821_5         -6.2      0000000000000000000000000000100000010000000000000101000000000001000000n000000000000000000001000000000000000000000000000000000000001010000000000100000000000001000

Notice that there is an ``n`` character that replace one of the bit.

Omit single interaction in several residue
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following configuration file (`examples-input/08-na_omit_interaction/vina-omit-residues.txt`)
shows how you can omit hydrophobic interaction on ARG150 TRP177 ARG223 residues ::

    docking_method    vina    # plants or vina
    docking_conf      ../03-na_vina/vina-003.conf

    similarity_coef   tanimoto mcconnaughey

    full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000  00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000  00010101000000100000000000000000000100000000000001010000000000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

    omit_interaction  hydrophobic  ARG150 TRP177 ARG223

    residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
    residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

    full_outfile omit-residues.csv
    sim_outfile omit-residues-similarity.csv
    logfile omit-residues.log

Here is the excerpt from the bitstring output using the above configuration file ::

    129821_1         -6.9      00000001000000000000000000001000000n0010000000000n00000000000001000000n000000000000010000001000000000000000000000000000000100000000010000000000100000000000001000
    129821_2         -6.8      00000000000000000000000000000000000n0000000000000n01000000000001000000n000000000000000000001000000000000000001010000000000100000000010000000000000000000000001000
    129821_3         -6.6      00000000000000000000000000001000001n0000000000000n01000000000001000000n000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000
    129821_4         -6.4      00000000000000000000000000001000100n0000000000000n00000000000001000000n000000000000000000001000000000000000000010000101000000000000000000000000000000000000000000
    129821_5         -6.2      00000000000000000000000000001000000n0000000000000n01000000000001000000n000000000000000000001000000000000000000000000000000000000001010000000000100000000000001000

Notice that there are three ``n``, each of them replace hydrophobic interaction in one residue.

Omit more than one interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

And here is how you can omit more than one interaction, this configuration file
(`examples-input/08-na_omit_interaction/vina-omit-interactions.txt`) shows you how
to omit the hydrophobic interaction on ARG223 and hydrogen bond (both as donor and acceptor) on
ARG292 ::

    docking_method    vina    # plants or vina
    docking_conf      ../03-na_vina/vina-003.conf

    similarity_coef   tanimoto mcconnaughey

    full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000  00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000  00010101000000100000000000000000000100000000000001010000000000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

    omit_interaction  hydrophobic  ARG223
    omit_interaction  h_bond  ARG292

    residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
    residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

    full_outfile omit-interactions.csv
    sim_outfile omit-interactions-similarity.csv
    logfile omit-interactions.log

Here is the excerpt from the bitstring output using the above configuration file ::

    129821_1         -6.9      0000000100000000000000000000100000010010000000000000000000000001000000n000000000000010000001000000000000000000000000000000nn0000000010000000000100000000000001000
    129821_2         -6.8      0000000000000000000000000000000000010000000000000101000000000001000000n000000000000000000001000000000000000001010000000000nn0000000010000000000000000000000001000
    129821_3         -6.6      0000000000000000000000000000100000110000000000000101000000000001000000n000000000000000000001000000000000000000000000000000nn0000000000000000000000000000000000000
    129821_4         -6.4      0000000000000000000000000000100010000000000000000000000000000001000000n000000000000000000001000000000000000000010000101000nn0000000000000000000000000000000000000
    129821_5         -6.2      0000000000000000000000000000100000010000000000000101000000000001000000n000000000000000000001000000000000000000000000000000nn0000001010000000000100000000000001000

Notice that there are three ``n``, the first one is replacing the hydrophobic interaction
on ARG223, while the second and third both are replacing the hydrogen bond interaction on
ARG292.Getting Started (HIPPOS on VINA)
===================================



There are two ways to use HIPPOS, the first one is to use HIPPOS without any
reference, and the second one is to use HIPPOS with reference and calculate
the similarity coefficient against the reference.

Generating Protein-Ligand Interaction Bitstring without reference
------------------------------------------------------------------------------------------

First, enter the examples folder and check the configuration example for this 
method, ``config-vina-na-notc.txt`` ::

	docking_method    vina    # plants or vina
	docking_conf      03-na_vina/vina-003.conf

	residue_name  ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	full_outfile vina_notc_ifp.csv
	logfile vina_notc.log

The configuration here is pretty much self-explaining. ``docking_method`` here is vina,
which correspond to the docking result we would like to analyse. Then the ``docking_conf``
is the configuration file used for docking, HIPPOS require this file to find the details 
about docking input and output from Vina.

For ``residue_name`` and ``residue_number`` check the 
:ref:`explanation for the configuration on generating reference with Hippos-genref <residue-numbering>` 
for more details.

Next, you can run HIPPOS by entering the following command: ::

	hippos config-vina-na-notc.txt

After the calculation finished HIPPOS will generate 2 files, vina_notc.log and vina_notc_ifp.csv. 
vina_notc.log contain the information about ligand name, number of poses, and the running time.
If output_mode set to simplified or combo there will be a table for bit position 
for each residue (useful for deciphering the simplified bitstring). The vina_notc_ifp.csv file will contain the ligand name,
pose number, energy from docking result, and the interaction bitstring as can be seen below:

.. image:: 11-vina-noref.png
	:alt: ifp results
	:align: center

In the next section you will learn how to not only generate the interaction bitstring but
also calculate the similarity coefficient using reference bitstring.

Generating Protein-Ligand Interaction Bitstring and Similarity coefficient 
-----------------------------------------------------------------------------

This procedure will require a reference bitstring which can be generated using
hippos-genref included in the package. Open this :doc:`link <getting-started-genref>`
to learn how to generate the reference bitstring.

After we acquire the full bitstring we can use it for reference bitstring
as shown in ``examples/04-na_config_default/config-vina-na-tc-mc.txt`` which is the configuration file for HIPPOS ::

	docking_method    vina    # plants or vina
	docking_conf      ../03-na_vina/vina-003.conf

	similarity_coef   tanimoto mcconnaughey

	full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000 00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000 00010101000000100000000000000000000100000000000001010100100000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

	residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	full_outfile vina_full_ifp.csv
	sim_outfile vina_similarity.csv
	logfile vina.log

**Always remember that full_ref should be using the full bitstring from reference.
Using bitstring reference of different length will cause an error and the program will stop.**

Here the residue_name and residue_number must be the same as the one used for reference
bitstring generation, and you have to set ``similarity_coef`` value, such as tanimoto
or mccounaughey or both of them.

``docking_method`` here is vina which correspond to the docking result we would like to
analyse. Then the ``docking_conf`` is the configuration file used for docking, HIPPOS require
this file to find the details about docking input and output from Vina.

Next, run HIPPOS with the following command inside ``examples`` directory: ::

	hippos config-vina-na-tc-mc.txt

there will be 3 output file vina.log, vina_similarity.csv, and vina_full_ifp.csv. The vina_full_ifp.csv
will be the same as the one without reference above. The vina_similarity.csv contain the 
similarity coefficient for every pose comparison. Notice
that there are 6 similarity coefficient results which correspond to Tanimoto
coefficient and McConnaughey coefficient calculation for 3 reference bitstring

.. image:: 21-vina-similarity.png
	:alt: tc results
	:align: center

Last but not least the hippos.log contain the information about ligand name, number of
poses, similarity coefficient used, and the table for bit position for each residue
(only appear when output_mode set to simplified, useful for deciphering the
simplified bitstring), and the total time taken.

.. image:: 22-vina-log.png
	:alt: tc results
	:align: center

Generating Protein-Ligand Interaction Bitstring and Similarity coefficient (without Backbone)
---------------------------------------------------------------------------------------------

The example above is the default setting where the interaction between ligand and protein backbone
is included. To omit the interaction between ligand and protein backbone, we need to set the
``output_mode`` value to ``full_nobb``, and in order for this setting to work properly we also need to
change the bitstring reference (``full_nobb_ref``) accordingly 
(:ref:`see how to generate bitstring reference without backbone <genref-nobb>`).
Here is the content of configuration file example ``examples/05-na_config_nobb/config-vina-na-tc-mc.txt`` ::

	docking_method    vina    # plants or vina
	docking_conf      ../03-na_vina/vina-003.conf

	similarity_coef   tanimoto mcconnaughey

	output_mode full_nobb

	full_nobb_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000  00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000  00010101000000100000000000000000000100000000000001010000000000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

	residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	full_nobb_outfile vina_nobb_ifp.csv
	sim_outfile vina_similarity.csv
	logfile vina.log

**Always remember that full_nobb_ref should be using the full_nobb bitstring from reference.
Using bitstring reference of different length will cause an error and the program will stop.**

Like before, run ``hippos`` with the following command: ::

	hippos config-vina-na-tc-mc.txt

Just like before, 3 output file will be generated, but the fingerprint (``vina_nobb_ifp.csv``)
and vina_similarity.csv will be different.

Generating Simplified Interaction Bitstring and Similarity coefficient
----------------------------------------------------------------------

It is also possible to calculate simplified interaction between ligand and protein. To do so set the
``output_mode`` value to ``simplified``, and in order for this setting to work properly we also need to
change the bitstring reference (``simplified_ref``) accordingly 
(:ref:`see how to generate simplified bitstring reference <genref-simplified>`).
Here is the content of configuration file example ``examples/06-na_config_simplified/config-vina-na-tc-mc.txt`` ::

	docking_method    vina    # plants or vina
	docking_conf      ../03-na_vina/vina-003.conf

	similarity_coef   tanimoto mcconnaughey

	output_mode simplified

	simplified_ref  0010000000000100000100000110000000010000000000000110000011000000100  0111000000000100000101000110000010000000000111010010000011000000000  0111001000000100000101000110000010010000001000100111000001000000000

	residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	simplified_outfile vina_simplified_ifp.csv
	sim_outfile vina_similarity.csv
	logfile vina.log

**Always remember that simplified_ref should be using the simplified bitstring from reference.
Using bitstring reference of different length will cause an error and the program will stop.**

Like before, run ``hippos`` with the following command: ::

	hippos config-vina-na-tc-mc.txt

Just like before, 3 output file will be generated, but the fingerprint (``vina_simplified_ifp.csv``)
and vina_similarity.csv will be different.

Generating Multiple Interaction Bitstring
-----------------------------------------

Last but not least, multiple output_mode is also allowed in generation interaction bitstring but without calculation of similarity coefficient. Here is the content of the configuration file example ``examples/07-na_config_multiple/config-vina-na.txt`` ::

	docking_method    vina    # plants or vina
	docking_conf      ../03-na_vina/vina-003.conf

	output_mode full full_nobb simplified

	residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	full_outfile vina_full.csv
	full_nobb_outfile vina_nobb.csv
	simplified_outfile vina_simplified_ifp.csv
	logfile vina.log

Like before, run ``hippos`` with the following command: ::

	hippos config-vina-na.txt

Now, four output file will be generated, three for three different output, and one for the log file.

..	
	Preparing docking file for VINA
	--------------------------------
	
	Preparing configuration file for VINA
	--------------------------------------
	
	Running simple docking in VINA
	------------------------------
	
	Running IFP analysis with HIPPOS
	--------------------------------
HIPPOS-genref Configuration Options
===================================

There are two kinds of options, *essential* and *optional*. When nothing declared it means the option is essential. Some options have default value. Also keep in mind that **all file names must never use space!** So use underscore instead.

You can also inserting comments at the beginning of the line or after the option-value pair by inserting the ``#`` sign before the comment. Everything inserted after the ``#`` sign will be ignored by the software.

Input Options
-------------

* **residue_name**
	value: ``list of residue_name`` , eg. ``ASP107 SER111 THR112``
	
	The list of residue_name, each residue separated by space. It is used in PLANTS post-analysis but not in VINA analysis as pdbqt don't hold the residue_name-residue_number pair field. However it is highly recommended to define it in VINA post analysis as it will be included in output file, making the results easier to interpret.
	
* **residue_number**
	value: ``list of residue_number`` , eg. ``80 84 85``
	
	The list of residue number, each number separated by space. Essential option in VINA post-analysis but optional in PLANTS post-analysis.
	
* **proteins**
	value: ``list of reference_protein``, eg. ``protein1.mol2 protein2.mol2 protein3.mol2``
	
	The list of reference_protein structure, each reference separated by space. Must be in mol2 or pdbqt format.

* **ligands**
	value: ``list of reference_ligand`` , eg. ``ligand1.mol2 ligand2.mol2 ligand3.mol2``
	
	The list of reference_ligand structure, each reference separated by space. Must be in mol2 or pdbqt format.
	Please note that the first reference_protein will be paired with first reference_ligand, the second reference_protein
	will be paired with second reference_ligand, and so on.
	
Output Options
--------------

* **outfile** (optional)
	value: ``output_file_name``, default ``genref-results.txt``

Advanced Options
----------------

* **output_mode** (optional)
	value: ``output_mode`` , default: ``full``, options: ``full``, ``full_nobb``, and ``simplified``
	
	``output_mode`` define the bistring calculation and output. ``full`` means all interaction including the 
	backbone interaction taken into account. ``full_nobb`` means backbone interactions are omitted.
	``simplified`` means irrelevant interaction omitted (see :ref:`Simplified Bitstring Rule<simplified-rule>`).
	When this option not used ``output_mode full`` will be used.
Getting Started (HIPPOS on PLANTS)
===================================



There are two ways to use HIPPOS, the first one is to use HIPPOS without any
reference, and the second one is to use HIPPOS with reference and calculate
the similarity coefficient against the reference.

Generating Protein-Ligand Interaction Bitstring without reference
------------------------------------------------------------------------------------------

First, enter the examples folder and check the configuration example for this 
method, ``config-plants-na-notc.txt`` ::

	docking_method    plants    # plants or vina
	docking_conf      02-na_plants/plants-003.conf

	residue_name  	  ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number    40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	full_outfile plants_notc_ifp.csv
	logfile plants_notc.log

The configuration here is pretty much self-explaining. ``docking_method`` here is plants,
which correspond to the docking result we would like to analyse. Then the ``docking_conf``
is the configuration file used for docking, HIPPOS require this file to find the details 
about docking input and output from PLANTS.

For ``residue_name`` and ``residue_number`` check the 
:ref:`explanation for the configuration on generating reference with Hippos-genref <residue-numbering>` 
for more details.

Next, you can run HIPPOS by entering the following command: ::

	hippos config-plants-na-notc.txt

After the calculation finished HIPPOS will generate 2 files, hippos.log and plants_notc_ifp.csv. 
plants_notc.log contain the information about ligand name, number of poses, and the running time.
If output_mode set to simplified or combo there will be a table for bit position 
for each residue (useful for deciphering the simplified bitstring). The plants_notc_ifp.csv file will contain the ligand name,
pose number, energy from docking result, and the interaction bitstring as can be seen below:

.. image:: 12-plants-noref.png
	:alt: ifp results
	:align: center

In the next section you will learn how to not only generate the interaction bitstring but
also calculate the similarity coefficient using reference bitstring.

Generating Protein-Ligand Interaction Bitstring and Similarity coefficient
-----------------------------------------------------------------------------

This procedure will require a reference bitstring which can be generated using
hippos-genref included in the package. Open this :doc:`link <getting-started-genref>`
to learn how to generate the reference bitstring.

After we acquire the full bitstring we can use it for reference bitstring
as shown in ``examples/04-na_config_default/config-plants-na-tc-mc.txt`` which is the configuration file for HIPPOS ::

	docking_method    plants    # plants or vina
	docking_conf      ../02-na_plants/plants-003.conf

	similarity_coef   tanimoto mcconnaughey

	full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000 00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000 00010101000000100000000000000000000100000000000001010100100000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

	residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	full_outfile plants_full_ifp.csv
	sim_outfile plants_similarity.csv
	logfile plants.log

**Always remember that full_ref should be using the full bitstring from reference.
Using bitstring reference of different length will cause an error and the program will stop.**

Here the residue_name and residue_number must be the same as the one used for reference
bitstring generation, and you have to set ``similarity_coef`` value, such as tanimoto
or mccounaughey or both of them.

``docking_method`` here is plants which correspond to the docking result we would like to
analyse. Then the ``docking_conf`` is the configuration file used for docking, HIPPOS require
this file to find the details about docking input and output from PLANTS.

Next, run HIPPOS with the following command inside ``examples`` directory: ::

	hippos config-plants-na-tc-mc.txt

there will be 3 output file plants.log, plants_similarity.csv, and plants_full_ifp.csv. The plants_full_ifp.csv
will be the same as the one without reference above. The plants_similarity.csv contain the 
similarity coefficient for every pose comparison. Notice
that there are 6 similarity coefficient results which correspond to Tanimoto
coefficient and McConnaughey coefficient calculation for 3 reference bitstring

.. image:: 23-plants-similarity.png
	:alt: tc results
	:align: center

Last but not least the plants.log contain the information about ligand name, number of
poses, similarity coefficient used, and the table for bit position for each residue
(only appear when output_mode set to simplified, useful for deciphering the
simplified bitstring), and the total time taken.

.. image:: 24-plants-log.png
	:alt: tc results
	:align: center

Generating Protein-Ligand Interaction Bitstring and Similarity coefficient (without Backbone)
---------------------------------------------------------------------------------------------

The example above is the default setting where the interaction between ligand and protein backbone
is included. To omit the interaction between ligand and protein backbone, we need to set the
``output_mode`` value to ``full_nobb``, and in order for this setting to work properly we also need to
change the bitstring reference (``full_nobb_ref``) accordingly
(:ref:`see how to generate bitstring reference without backbone <genref-nobb>`).
Here is the content of configuration file example ``examples/05-na_config_nobb/config-plants-na-tc-mc.txt``) ::

	docking_method    plants    # plants or vina
	docking_conf      ../02-na_plants/plants-003.conf

	similarity_coef   tanimoto mcconnaughey

	output_mode full_nobb

	full_nobb_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000  00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000  00010101000000100000000000000000000100000000000001010000000000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

	residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	full_nobb_outfile plants_nobb_ifp.csv
	sim_outfile plants_similarity.csv
	logfile plants.log

**Always remember that full_nobb_ref should be using the full_nobb bitstring from reference.
Using bitstring reference of different length will cause an error and the program will stop.**

Like before, run ``hippos`` with the following command: ::

	hippos config-plants-na-tc-mc.txt

Just like before, 3 output file will be generated, but the fingerprint (``plants_nobb_ifp.csv``)
and plants_similarity.csv will be different.

Generating Simplified Interaction Bitstring and Similarity coefficient
----------------------------------------------------------------------

It is also possible to calculate simplified interaction between ligand and protein. To do so set the
``output_mode`` value to ``simplified``, and in order for this setting to work properly we also need to
change the bitstring reference (``simplified_ref``) accordingly 
(:ref:`see how to generate simplified bitstring reference <genref-simplified>`).
Here is the content of configuration file example ``examples/06-na_config_simplified/config-plants-na-tc-mc.txt`` ::

	docking_method    plants    # plants or vina
	docking_conf      ../02-na_plants/plants-003.conf

	similarity_coef   tanimoto mcconnaughey

	output_mode simplified

	simplified_ref  0010000000000100000100000110000000010000000000000110000011000000100  0111000000000100000101000110000010000000000111010010000011000000000  0111001000000100000101000110000010010000001000100111000001000000000

	residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	simplified_outfile plants_simplified_ifp.csv
	sim_outfile plants_similarity.csv
	logfile plants.log

**Always remember that simplified_ref should be using the simplified bitstring from reference.
Using bitstring reference of different length will cause an error and the program will stop.**

Like before, run ``hippos`` with the following command: ::

	hippos config-plants-na-tc-mc.txt

Just like before, 3 output file will be generated, but the fingerprint (``plants_simplified_ifp.csv``)
and plants_similarity.csv will be different.

Generating Multiple Interaction Bitstring
-----------------------------------------

Last but not least, multiple output_mode is also allowed in generation interaction bitstring but without calculation of similarity coefficient. Here is the content of the configuration file example ``examples/07-na_config_multiple/config-plants-na.txt`` ::

	docking_method    plants    # plants or vina
	docking_conf      ../02-na_plants/plants-003.conf

	output_mode full full_nobb simplified

	residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	full_outfile plants_full.csv
	full_nobb_outfile plants_nobb.csv
	simplified_outfile plants_simplified_ifp.csv
	logfile plants.log

Like before, run ``hippos`` with the following command: ::

	hippos config-plants-na.txt

Now, four output file will be generated, three for three different output, and one for the log file.

..	
	Preparing docking file for PLANTS
	---------------------------------
	
	Preparing configuration file for PLANTS
	---------------------------------------
	
	Running simple docking in PLANTS
	--------------------------------
	
	Running IFP analysis with HIPPOS
	--------------------------------
Installation
============

Quick Installation (Recommended)
--------------------------------

The easiest way to install HIPPOS is using `Conda <https://docs.anaconda.com/anaconda/install/>`_, you
can choose either Anaconda or Miniconda. If you never used Conda before then most likely Miniconda
is more suitable for you.

After you installed Conda in your machine, install HIPPOS using the following command:

``conda install -c conda-forge pyplif-hippos``

Conda will deal with the dependencies like Open Babel, Bitarray, and Numpy libraries. Therefore you
won't need the extra step below.

To check if HIPPOS installed correctly try these commands:

|  ``hippos``
|  ``hippos-genref``

If HIPPOS installed correctly you should get a message that inform you the configuration file not found.
Next, you can just jump to Getting Started with PLANTS or Vina tutorial. But if you prefer not to use Conda
then you can install HIPPOS using the instructions below.

Requirement
-----------

* Python >= 2.6 or >= 3.6
* Open Babel (library) >= 2.2
* Python-OpenBabel >= 2.1
* Python-BitArray
* Python-Numpy >= 1.3

Getting HIPPOS
--------------

You can get HIPPOS at Github by cloning the repository or download the code `here <https://github.com/radifar/PyPLIF-HIPPOS>`_.

Installing on Linux
-------------------

In Ubuntu you can simply enter these commands to install the requirements above (Python is already installed in Linux):

| ``sudo apt-get install openbabel``
| ``sudo apt-get install python-openbabel``
| ``sudo apt-get install python-bitarray``
| ``sudo apt-get install python-numpy``

In Fedora you can enter these commands instead:

| ``sudo yum install openbabel``
| ``sudo yum install python-openbabel``
| ``sudo yum install python-bitarray``
| ``sudo yum install numpy``

After the requirements fulfilled

After all of the dependencies installed, you can install HIPPOS by opening
the terminal and enter the HIPPOS directory and run setup.sh like so:

``./setup.sh``

You will be prompted with a question whether to install HIPPOS in ``[HOME_DIRECTORY]/.hippos``
which is a hidden directory inside your ``HOME_DIRECTORY``. Or would you rather
choose your own directory. If you want to keep the default then just press `enter`
or if you would like to choose your own installation directory you can type `y`
and press enter, then you have to type the installation directory for example
``/home/radifar/apps/hippos`` then press enter.

If HIPPOS installed successfully then a message like 'HIPPOS successfully
installed' should appear. When HIPPOS is already installed and you're running
setup.sh a message like 'HIPPOS is already installed' will appear, and the
installation process will stop and exit.

If you would like to install newer version of HIPPOS and overwrite the old
one then all you need to do is by adding 'force' option to setup.sh like so:

``./setup.sh force``

To test if HIPPOS had been installed, *open new command line window* and type the following:

``hippos``

Then press enter. Please note that it is imperative to open the new command
line window right after the installation to allow the alias in your system get updated.
Next type the following:

``hippos-genref``

Then press enter

If there are no error message then the installation is success.

Installing on Windows 10
------------------------

Windows 10 provides Ubuntu in their `Microsoft Store <https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6>`_. 
After you successfully add Ubuntu on your Windows, the next step is to do the same installation steps on Ubuntu.

.. HIPPOS documentation master file, created by
   sphinx-quickstart on Fri Dec  7 13:14:20 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyPLIF HIPPOS's documentation!
=========================================

PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results of AutoDock Vina and PLANTS
----------------------------------------------------------------------------------------------------------

.. image:: hippopotamus.png
   :alt: Icons made by Freepik from Flaticon is licensed by CC 3.0 BY
   :align: center
   :scale: 45%

.. raw:: html
   :file: attribution-hippos.html

Welcome to PyPLIF-HIPPOS's documentation. PyPLIF-HIPPOS is an upgraded version of `PyPLIF <https://github.com/radifar/pyplif/>`_ (**Python-based Protein-Ligand Interaction Fingerprinting**), a tool for molecular docking post-analysis. It will translate the 3D coordinates of both ligand(s) (generated from docking simulation) and protein into a series of *interaction bitstring* (also known as *Interaction Fingerprint*) (see image below). **HIPPOS** (/ˌhipoʊz/) is a recursive acronym of **HIPPOS Is PyPLIF On Steroids**. From this point forward, PyPLIF-HIPPOS is simplified to HIPPOS.

Compared to PyPLIF, HIPPOS is not only faster and able to generate more customized interaction bitstring, but also supports both `PLANTS <https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/pharmazie-und-biochemie/pharmazie/pharmazeutische-chemie/pd-dr-t-exner/research/plants/>`_ & `VINA <http://vina.scripps.edu/>`_! More over, unlike its predecessor it is (far) more well-documented. 

.. image:: toc-abstract-graphics_small.png
   :alt: Table of Content Abstract Graphic JCIM
   :align: center

Reprinted with permission from https://doi.org/10.1021/acs.jcim.0c00305. Copyright 2020 American Chemical Society.

.. image:: pyplif.png
   :alt: PyPLIF output from PyPLIF publication
   :align: center

.. raw:: html
   :file: attribution-pyplif.html

Citing HIPPOS
-------------

If you are using HIPPOS please cite this paper:

Istyastono, E., Radifar, M., Yuniarti, N., Prasasty, V. and Mungkasi, S., 2020.
PyPLIF HIPPOS: A Molecular Interaction Fingerprinting Tool for Docking Results
of AutoDock Vina and PLANTS. Journal of Chemical Information and Modeling.
https://doi.org/10.1021/acs.jcim.0c00305

Acknowledgment
--------------

This project has received funding from the `Indonesian National Research and Innovation Agency <https://international.ristekdikti.go.id/>`_
under grant agreement No. 807.7/LL5/PG/2020. This project has been restructured based on the
`MOLSSI Computational Molecular Science Python Cookiecutter <https://github.com/molssi/cookiecutter-cms>`_
version 1.3, and benefited greatly from `MOLSSI Python Package Development Best Practices <https://molssi.org/2020/04/20/may-webinar-series-python-package-development/>`_
workshop.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   installation
   parameter-options
   getting-started-genref
   hippos-genref-config
   getting-started-vina
   getting-started-plants
   advanced-usage
   hippos-config
   license

Generating Reference Bitstring using HIPPOS-genref
==================================================

Generating Reference Bitstring with Backbone (Default Setting)
--------------------------------------------------------------

Reference bitstring is an essential requirement for similarity coefficient (eg. Tanimoto
or McConnaughey coefficient) calculation, which is the common method for comparing
the interaction fingerprinting of a test compound and a reference (native ligand)
interactions on certain protein.

To generate Reference bitstring, first of all, you need to open the command prompt and enter the 
``examples\01-na_reference`` folder. This example uses the Neuraminidase enzyme
for two reasons, first, it is one of the enzymes used in DUD-E (Directory of Useful
Decoy Enhanced) therefore you could use it to measure the effect of interaction fingerprinting
on the enrichment factor. Second, it can demonstrate all of the seven interaction types in
interaction fingerprinting.

As you can see there are three folders and two txt configuration files. Each folder represents
a crystal structure of Neuraminidase, and contain the original PDB file and the
split component (protein, ligand, and water) generated with 
`SPORES <https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/pharmazie-und-biochemie/pharmazie/pharmazeutische-chemie/pd-dr-t-exner/research/spores/>`_. PDB file alone 
can not be used as the reference, it has to be in mol2 or pdbqt to ensure that the 
atom typing, charge assignment, and protonation is identical to the docking environment 
(whether for PLANTS or VINA). So the PDB files here only act as
the source if you want to use PDBQT files as the reference instead. **1b9s**,
**1b9t**, and **1b9v** are the PDB ID of the same Neuraminidase, where each of them
bound to  different ligand (**FDI**, **RAI**, and **RA2** respectively). Therefore 
the protein name and ligand name should use protein.mol2 and the corresponding ligand 
file name as you can see in ``genref-config.txt``. ::

	# first residue is 77
	residue_name	ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number	40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	proteins	1b9s/protein.mol2 1b9t/protein.mol2 1b9v/protein.mol2
	ligands		1b9s/ligand_FDI468_0.mol2 1b9t/ligand_RAI468_0.mol2 1b9v/ligand_RA2468_0.mol2

	outfile		ref-results.txt


.. _residue-numbering:

The first line is merely the commented line, everything started with ``#`` sign will
be ignored by ``hippos`` and ``hippos-genref``. If you open the mol2 protein file with a text editor
you will see that the first residue is Glutamate with residue number 77. However, when
the file is read by Openbabel it will count as residue number 1. 

.. image:: first-residue.png
	:alt: First residue of Neuraminidase
	:align: center


This is where things started to get tricky, because we have to supply both ``residue_name``
and ``residue_number`` properly, or else it will not work as to how we want it to be.
``residue_name`` can be acquired easily by converting your protein into mol2 format, while
the corresponding ``residue_number`` must be retrieved from the column before residue name 
(eg. residue name ARG116 and GLU117 correspond to residue number 40 and 41 respectively):

.. image:: arg116-glu117.png
	:alt: Residue ARG116 and GLU117
	:align: center

The ``residue_name`` and ``residue_number`` in ``genref-config`` above are retrieved by
visualizing any residue within 5 angstroms from the native ligand using VMD (you can use
any other molecule visualization tool), regardless of how important the residue in
enzyme inhibition.

The next lines are ``proteins`` and ``ligands``, notice that there are 3 protein molecules
and 3 ligand molecules which means that there are 3 protein-ligand pairs as references. Where
the first protein will be matched with first ligand and so on. However in most cases, one
protein-ligand pair is enough, this example uses 3 protein-ligand pairs as a demonstration
of multiple references.

The last line is the output file name, it is optional so when not defined the output file
will be genref-results.txt.

After we understand the input file and the configuration file, hippos-genref could be run with the following command: ::

	hippos-genref genref-config.txt

After hippos-genref finished file ``ref-results.txt`` will be generated.

.. image:: 01-ref-results.png
	:alt: ref-results.txt from hippos-genref
	:align: center

Inside ref-results.txt we can see that there are 3 results from 3 protein-ligand pairs.
Each result consisted of protein-ligand pair name and interaction bitstring where each 
residue represented by 7 bit of interactions from both the side-chain and the backbone.

.. _genref-nobb:

Generating Reference Bitstring without Backbone
-----------------------------------------------

Sometimes we would like to omit the interaction between ligand and the backbone protein.
In that case, we should change the ``output_mode`` to ``full_nobb`` by adding 
``output_mode full_nobb`` to our hippos-genref config file as appear in ``genref-config-nobb.txt`` ::

	# first residue is 77
	residue_name  ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number  40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	proteins    1b9s/protein.mol2 1b9t/protein.mol2 1b9v/protein.mol2
	ligands     1b9s/ligand_FDI468_0.mol2 1b9t/ligand_RAI468_0.mol2 1b9v/ligand_RA2468_0.mol2

	output_mode full_nobb

	outfile     ref-results-nobb.txt

Now run hippos-genref again with the following command: ::

	hippos-genref genref-config-nobb.txt

After hippos-genref finished file ``ref-results-nobb.txt`` will be generated.

.. image:: 02-ref-results-nobb.png
	:alt: ref-results-nobb.txt from hippos-genref
	:align: center

Just like in the default setting, it will generate 3 results. And although they
appear the same as before, this time the bitstrings are generated without
taking backbone atoms into account.

.. _genref-simplified:

Generating Simplified Reference Bitstring
-----------------------------------------

It is also possible to calculate simplified interaction between ligand and the backbone protein.
In that case, we should change the ``output_mode`` to ``simplified`` by adding 
``output_mode simplified`` to our hippos-genref config file as appear in ``genref-config-simplified.txt`` ::

	# first residue is 77
	residue_name  ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
	residue_number  40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

	proteins    1b9s/protein.mol2 1b9t/protein.mol2 1b9v/protein.mol2
	ligands     1b9s/ligand_FDI468_0.mol2 1b9t/ligand_RAI468_0.mol2 1b9v/ligand_RA2468_0.mol2

	output_mode simplified

	outfile     ref-results-simplified.txt

Now run hippos-genref again with the following command: ::

	hippos-genref genref-config-simplified.txt

After hippos-genref finished file ``ref-results-simplified.txt`` will be generated.

.. image:: 04-ref-results-simplified.png
	:alt: ref-results-simplified.txt from hippos-genref
	:align: center

Just like in the default setting, it will generate 3 results. And although they
appear the same as before, this time the bitstrings are simplified.

..
	Generating Multiple Reference Bitstring
	---------------------------------------

	HIPPOS-genref is also capable of outputting multiple output_mode. Currently,
	there are three output_mode available: full, full_nobb, and :ref:`simplified<simplified-rule>`. To
	generate more than one output_mode the output_mode parameter can be concatenated
	by separating each parameter with blank space. In the example below (``genref-config-multi.txt``), 
	two parameters are used: ``full`` and ``simplified`` ::

		# first residue is 77
		residue_name  ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
		residue_number  40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

		proteins    1b9s/protein.mol2 1b9t/protein.mol2 1b9v/protein.mol2
		ligands     1b9s/ligand_FDI468_0.mol2 1b9t/ligand_RAI468_0.mol2 1b9v/ligand_RA2468_0.mol2

		output_mode full simplified

		outfile     ref-results-multi.txt

	Using the above configuration file, hippos-genref could be run with the
	following command: ::

		hippos-genref genref-config-multi.txt

	After hippos-genref finished file ``ref-results-multi.txt`` will be generated.

	.. image:: 03-ref-results-multi.png
		:alt: ref-results-multi.txt from hippos-genref
		:align: center

	In this example, for every protein ligand pair two bitstrings are generated. One
	is for full, and the other one is for simplified bitstring.