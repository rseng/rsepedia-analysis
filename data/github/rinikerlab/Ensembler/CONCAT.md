
![](.img/EnsemblerLogo_with_whiteBackround.png)


Welcome to Ensembler
==============================

[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/rinikerlab/ensembler/workflows/CI/badge.svg)](https://github.com/rinikerlab/ensembler/actions?query=branch%3Amaster+workflow%3ACI)
[![codecov](https://codecov.io/gh/rinikerlab/Ensembler/branch/master/graph/badge.svg)](https://codecov.io/gh/rinikerlab/Ensembler/branch/master)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/rinikerlab/Ensembler.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/rinikerlab/Ensembler/context:python)
![Build package](https://github.com/rinikerlab/Ensembler/workflows/Python%20package/badge.svg)
[![Documentation](https://img.shields.io/badge/Documentation-here-white.svg)](https://rinikerlab.github.io/Ensembler/index.html)

## Description
Ensembler is a python package that allows fast and easy access to one and two-dimensional model systems simulations. It enables method development using small test systems and helps deepen understanding of a broad spectrum of molecular dynamics (MD) methods, starting from basic techniques to enhanced sampling and free energy calculations. It is easy to install, fast, increases shareability, comparability, and reproducibility of scientific code developments. 
Here, we provide insights into the package's implementation, usage, and an application example for free energy calculation.
There is an article about Ensembler in JCIM, check it out: https://doi.org/10.1021/acs.jcim.0c01283

**Quick Start to Simulations:**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rinikerlab/Ensembler/blob/master/examples/Tutorial_Simulations.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rinikerlab/Ensembler/master?filepath=examples%2FTutorial_Simulations.ipynb)

## Contents
The full documentation can be found here:  https://rinikerlab.github.io/Ensembler/

### potentials - Potential Energy Functions

  Implement mathematical functions of interest for modeling purposes, for example, in chemistry.
  Implementation of new potential energy functions is straightforward, as there are only a few functions that need to be overwritten.
  Implemented Potentials: Harmonic Oscillator, Wave functions, etc. 
  Also, different dimensionalities can be used, like 1D, 2D, and ND.
 
### samplers - Sampling Methods

   This module provides integrators for integrating potential functions. E.g., Monte Carlo, Velocity Verlet.
   
### systems - Used for  Simulations

   This module is used to set up a simulation with a system class that contains a potential energy function and a sampling method, plus optional other parameters.

### ensemble - Multi Replica Approaches

   This module contains the implementation of the ConveyorBelt and will contain in future additional Replica Exchange approaches.

### visualization

   This module contains predefined visualization, animation and widget functions.

### analysis

   This module contains at the moment Free Energy Calculations.

## Tutorials and Examples:
You can try Ensembler online via the binder or Google colab links below or on your machine with our provided jupyter notebooks showing use cases for Ensembler with visualizations.
The Binder links might take some time to build the repo; depending on the traffic, it can take up to 10 minutes to be built.

### Tutorials
We provide short introductions into how potential energy functions can be used and sampled in simulations
 in Ensembler with the jupyter notebooks in the example folder that can be also accessed by the binder links.

Here is an small example of a simple simulation:

CODE:

    ##Imports
    from ensembler.potentials.OneD import fourWellPotential
    from ensembler.samplers.stochastic import langevinIntegrator
    from ensembler.system import system
    from ensembler.visualisation.plotSimulations import oneD_simulation_analysis_plot
    
    ##Simulation Setup
    V = fourWellPotential(Vmax=4, a=1.5, b=4.0, c=7.0, d=9.0,  ah=2., bh=0., ch=0.5, dh=1.)
    sampler = langevinIntegrator(dt=0.1, gamma=10)
    sys = system(potential=V, sampler=sampler, start_position=4, temperature=1)
    
    ##Simulate
    sys.simulate(steps=1000)
    
    ##Visualize
    positions = np.linspace(start=0, stop=10, num=1000) #phase space to be visualized
    oneD_simulation_analysis_plot(system=sys, title="Langevin Simulation", limits_coordinate_space=positions)
    
OUT:
![Not Visible with pip](.img/langevin_simulation.png)

In the following links you can find more features of Ensembler.

Potentials: 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rinikerlab/Ensembler/blob/master/examples/Tutorial_Potentials.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rinikerlab/Ensembler/master?filepath=examples%2FTutorial_Potentials.ipynb)

Simulations: 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rinikerlab/Ensembler/blob/master/examples/Tutorial_Simulations.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rinikerlab/Ensembler/master?filepath=examples%2FTutorial_Simulations.ipynb)

### Examples
Examples are advanced jupyter notebooks, covering a specific topic, going deeper into the methodology and having nice visualizations.

Enhanced Sampling: 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rinikerlab/Ensembler/blob/master/examples/Example_EnhancedSampling.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rinikerlab/Ensembler/master?filepath=examples%2FExample_EnhancedSampling.ipynb)

Free Energy Calculations: 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rinikerlab/Ensembler/blob/master/examples/Example_FreeEnergyCalculationSimulation.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rinikerlab/Ensembler/master?filepath=examples%2FExample_FreeEnergyCalculationSimulation.ipynb)

Interactive ConveyorBelt: 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rinikerlab/Ensembler/blob/master/examples/Example_ConveyorBelt.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rinikerlab/Ensembler/master?filepath=examples%2FExample_ConveyorBelt.ipynb)

EDS-Potentials: 

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rinikerlab/Ensembler/blob/master/examples/Example_EDS.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rinikerlab/Ensembler/master?filepath=examples%2FExample_EDS.ipynb)



## How To Install
There are many ways to get the Ensembler package.
  * If you want to use Ensembler (don't forget to download the examples):
        
        pip install ensembler-rinikerlab
         
  * If you want to use Ensembler, see the example notebooks or develop Ensembler:
        
        git clone https://github.com/rinikerlab/Ensembler.git
        cd Ensembler
        python setup.py install
          
  * If you want to use Ensembler, see the example notebooks or develop Ensembler - but not install the package:
    Add the path to the Ensembler repository; the requirements needed for the package can be used like in the following examples:
  
    * PIP:
    
            git clone https://github.com/rinikerlab/Ensembler.git
            cd Ensembler
            export PYTHONPATH=${PYTHONPATH}:${PWD}
            pip install -r devtools/pip_requirements/requirements_unix.txt
        
    * Anaconda:
   
            git clone https://github.com/rinikerlab/Ensembler.git
            cd Ensembler
            conda create -n Ensembler --file devtools/conda-envs/environment_unix.yml
            conda activate Ensembler

For windows, we also provide the required files in the same folders:
  * devtools/pip_requirements/requirements_windows.txt
  * devtools/conda-envs/environment_windows.yml).


## Contributing
If you would like to contribute to Ensembler, you are most welcome!
Just raise an issue or write us a mail.

## Authors

Benjamin Ries;
Stephanie M. Linker;
David F. Hahn

## Copyright

Copyright (c) 2020, Benjamin Ries, Stephanie M. Linker, David F. Hahn


### Acknowledgements
 
Project-based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
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
reported by contacting the project team at 'benjamin.schroeder@phys.chem.ethz.ch'. The project team will
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
e.g. `ensembler-0.1.2`, otherwise it will be appended with `+X` where `X` is the number of commits 
ahead from the last tag, and then `-YYYYYY` where the `Y`'s are replaced with the `git` commit hash.
# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into this ensembler. 

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free).
* [Fork](https://help.github.com/articles/fork-a-repo/) this repository on GitHub.
* On your local machine,
  [clone](https://help.github.com/articles/cloning-a-repository/) your fork of
  the repository.

## Making Changes

* Add some really awesome code to your local fork.  It's usually a [good
  idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/)
  to make changes on a
  [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/)
  with the branch name relating to the feature you are going to add.
* When you are ready for others to examine and comment on your new feature,
  navigate to your fork of ensembler on GitHub and open a [pull
  request](https://help.github.com/articles/using-pull-requests/) (PR). Note that
  after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.  Each commit added to the PR will be validated for
  mergability, compilation and test suite compliance; the results of these tests
  will be visible on the PR page.
* If you're providing a new feature, you must add test cases and documentation.
* When the code is ready to go, make sure you run the test suite using pytest.
* When you're ready to be considered for merging, check the "Ready to go"
  box on the PR page to let the ensembler devs know that the changes are complete.
  The code will not be merged until this box is checked, the continuous
  integration returns checkmarks,
  and multiple core developers give "Approved" reviews.

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
## Description
Provide a brief description of the PR's purpose here.

## Todos
Notable points that this PR has either accomplished or will accomplish.
  - [ ] TODO 1

## Questions
- [ ] Question1

## Status
- [ ] Ready to go# Templates Doc Directory

Add any paths that contain templates here, relative to  
the `conf.py` file's directory.
They are copied after the builtin template files,
so a file named "page.html" will overwrite the builtin "page.html".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
html_static_path = ['_templates']
```

## Examples of file to add to this directory
* HTML extensions of stock pages like `page.html` or `layout.html`
# Static Doc Directory

Add any paths that contain custom static files (such as style sheets) here,
relative to the `conf.py` file's directory. 
They are copied after the builtin static files,
so a file named "default.css" will overwrite the builtin "default.css".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
templates_path = ['_static']
```

## Examples of file to add to this directory
* Custom Cascading Style Sheets
* Custom JavaScript code
* Static logo images
# Static Doc Directory

Add any paths that contain custom static files (such as style sheets) here,
relative to the `conf.py` file's directory. 
They are copied after the builtin static files,
so a file named "default.css" will overwrite the builtin "default.css".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
templates_path = ['_static']
```

## Examples of file to add to this directory
* Custom Cascading Style Sheets
* Custom JavaScript code
* Static logo images
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
