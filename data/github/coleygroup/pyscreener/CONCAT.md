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
reported by contacting the project team at 'deg711@g.harvard.edu'. The project team will
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
# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into `pyscreener`. 

## Getting Started

1. Make sure you have a [GitHub account](https://github.com/signup/free).
1. [Fork](https://help.github.com/articles/fork-a-repo/) this repository on GitHub.
1. On your local machine, [clone](https://help.github.com/articles/cloning-a-repository/) your fork of the repository.

## Making Changes

1. Install the development version of `pyscreener`: `pip install -e .[dev]`
1. Install the the `pre-commit` hooks: `pre-commit install`
1. Add some changes to your local fork.  It's usually a [good idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/) to make changes on a [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/) with the branch name relating to the feature you are going to add.
1. If you're providing a new feature, you **must** add test cases and documentation.
1. Complete the checklist below
1. When you are ready for others to examine and comment on your new feature, navigate to your fork of `pyscreener` on GitHub and open a [pull request](https://help.github.com/articles/using-pull-requests/) (PR). Note that after you launch a PR from one of your fork's branches, all subsequent commits to that branch will be added to the open pull request automatically.  Each commit added to the PR will be validated for mergability, compilation and test suite compliance; the results of these tests will be visible on the PR page.
1. When you believe the PR has reached its final draft, check the "Ready to go" box below to let the `pyscreener` devs know that the changes are complete. The code will not be considered for merging until this box is checked, the continuous integration returns checkmarks, and multiple core developers give "Approved" reviews.

## Checklist

- [ ] formatted with black: `black pyscreener`
- [ ] linted with flake8: `flake8`
- [ ] unit tests added (if necessary)?
- [ ] all unit tests pass (if any core code was changed)?

## Ready to merge?

- [ ] yes!

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
[//]: # (Badges)
[![codecov](https://codecov.io/gh/coleygroup/pyscreener/branch/main/graph/badge.svg)](https://codecov.io/gh/coleygroup/pyscreener/branch/main)
[![CI](https://github.com/coleygroup/pyscreener/actions/workflows/CI.yaml/badge.svg)](https://github.com/coleygroup/pyscreener/actions/workflows/CI.yaml)
[![Documentation Status](https://readthedocs.org/projects/pyscreener/badge/?version=latest)](https://pyscreener.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/pyscreener.svg)](https://badge.fury.io/py/pyscreener)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03950/status.svg)](https://doi.org/10.21105/joss.03950)

# pyscreener
# A pythonic interface to high-throughput virtual screening software

## Overview
This repository contains the source of pyscreener, both a library and software for conducting HTVS via python calls

## Table of Contents
- [Overview](#overview)
- [Table of Contents](#table-of-contents)
- [Requirements](#requirements)
- [Installation](#installation)
- [Ray setup](#ray-setup)
- [Running pyscreener](#running-pyscreener-as-a-software)
- [Using pyscreener](#using-pyscreener-as-a-library)

## Installation

### General requirements
- python >= 3.8
- `numpy`, `openbabel`, `openmm`, `pdbfixer`, `ray`, `rdkit`, `scikit-learn`, `scipy`, and `tqdm`
- all corresponding software downloaded and located on your PATH or under the path of a specific environment variable (see [external software](#external-software) for more details.)

### Setup

0. (if necessary) [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
1. `conda env create -f environment.yml`
1. `conda activate pyscreener`
1. `pip install pyscreener` (or if installing from source, `pip install .`)
1. follow the corresponding directions below for the intended software

Before running `pyscreener`, be sure to first activate the environment: `conda activate pyscreener` (or whatever you've named your environment)

### external software
* vina-type software
  1. install [ADFR Suite](https://ccsb.scripps.edu/adfr/downloads/) and add `prepare_receptor` to your PATH. If this step was successful, the command `which prepare_receptor` should output `path/to/prepare_receptor`. This can be done via either:

      1. adding the entire `bin` directory to your path (you should see a command at the end of the installation process) or

      2. adding only `prepare_receptor` in the `bin` directory to your PATH as detailed [below](#adding-an-executable-to-your-path)
  
  
  1. install any of the following docking software: [vina 1.1.2](https://vina.scripps.edu/downloads/) (**note:** pyscreener does not work with vina 1.2), [qvina2](https://qvina.github.io/), [smina](https://sourceforge.net/projects/smina/), [psovina](https://cbbio.online/software/psovina/index.html) and [ensure the desired software executable is in a folder that is located on your path](#adding-an-executable-to-your-path)

* [DOCK6](http://dock.compbio.ucsf.edu/)
  1. [obtain a license for DOCK6](http://dock.compbio.ucsf.edu/Online_Licensing/dock_license_application.html)
  1. install DOCK6 from the download link and follow the [installation directions](http://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm#Installation)
  1. after ensuring the installation was installed properly, specify the DOCK6 environment variable as the path of the DOCK6 parent directory as detailed [below](#specifying-an-environment-variable). This is the directory that was unzipped from the tarball and is usually named `dock6`. It is the folder that contains the `bin`, `install`, etc. subdirectories.)
  1. install [sphgen_cpp](http://dock.compbio.ucsf.edu/Contributed_Code/sphgen_cpp.htm). On linux systems, this can be done:
      1. `wget http://dock.compbio.ucsf.edu/Contributed_Code/code/sphgen_cpp.1.2.tar.gz`
      1. `tar -xzvf sphgen_cpp.1.2.tar.gz`
      1. `cd sphgen_cpp.1.2`
      1. `make`
  1. place the sphgen_cpp executable (it should be `sphgen_cpp`) inside the `bin` subdirectory of the DOCK6 parent directory. If you've configured the environment variable already, (on linux) you can run: `mv sphgen_cpp $DOCK6/bin`
  1. [install chimera](https://www.cgl.ucsf.edu/chimera/download.html) and place the file on your PATH as detailed [below](#adding-an-executable-to-your-path)

#### adding an executable to your PATH
To add an executable to your PATH, you have three options:
1. create a symbolic link to the executable inside a directory that is already on your path: `ln -s FILE -t DIR`. Typically, `~/bin` or `~/.local/bin` are good target directories (i.e., `DIR`). To see what directories are currently on your path, type `echo $PATH`. There will typically be a lot of directories on your path, and it is best to avoid creating files in any directory above your home directory (`$HOME` on most *nix-based systems)
1. copy the software to a directory that is already on your path. Similar, though less preferred than the above: `cp FILE DIR`
1. append the directory containing the file to your PATH: `export PATH=$PATH:DIR`, where `DIR` is the directory containing the file in question. As your PATH must be configured each time run pyscreener, this command should also be placed inside your `~/.bashrc` or `~/.bash_profile` (if using a bash shell) to avoid needing to run the command every time you log in. _Note_: if using a non-bash shell, the specific file will be different.

#### specifying an environment variable
To set the `DOCK6` environment variable, run the following command: `export DOCK6=path/to/dock6`, where `path/to/dock6` is the **full** path of the DOCK6 parent directory mentioned above. As this this environment variable must always be set before running pyscreener, the command should be placed inside your `~/.bashrc` or `~/.bash_profile` (if using a bash shell) to avoid needing to run the command every time you log in. _Note_: if using a non-bash shell, the specific file will be different.

### Ray Setup
pyscreener uses [`ray`](https://docs.ray.io/en/master/index.html) as its parallel backend. If you plan to parallelize the software only across your local machine, don't need to do anything . However, if you wish to either (a.) limit the number of cores pyscreener will be run over or (b.) run it over a distributed setup (e.g., an HPC with many distinct nodes), you must manually start a ray cluster __before__ running pyscreener.

#### Limiting the number of cores
To do this, simply type `ray start --head --num-cpus N` before starting pyscreener (where `N` is the total number of cores you wish to allow pyscreener to utilize). Not performing this step will give pyscreener access to all of the cores on your local machine, potentially slowing down other applications.

#### Distributing across many nodes
While the precise instructions for this will vary with HPC cluster architecture, the general idea is to establish a ray cluster between the nodes allocated to your job. We have provided a sample SLURM submission script ([run_pyscreener_distributed_example.batch](run_pyscreener_distributed_example.batch)) to achieve this, but you may have to alter some commands depending on your system. For more information on this see [here](https://docs.ray.io/en/master/cluster/index.html). To allow pyscreener to connect to your ray cluster, you must set the `ip_head` and `redis_password` environment variables appropriately, where `ip_head` is the address of the head of your ray cluster, i.e., `IP:PORT` where `IP` is the IP address of the head node and `PORT` is the port that is running ray.

pyscreener writes a lot of intermediate input and output files (due to the inherent specifications of the underlying docking software.) Given that the primary endpoint of pyscreener is a list of ligands and associated scores (rather than the specific binding poses,) these files are written to each node's temporary directory (determined by `tempfile.gettempdir()`) and discarded at the end. If you wish to collect these files, pass the `--collect-all` flag in the program arguments or run the `collect_files()` method of your `VirtualScreen` object when your screen is complete.

*Note*: the `VirtualScreen.collect_files()` method is **slow** due to the need to send possibly a **bunch** of files over the network. This method should only be run **once** over the lifetime of a `VirtualScreen` object, as several intermediate calls will yield the same result as a single, final call.

*Note*: `tempfile.gettempdir()` returns a path that depends the values of specific environment variables (see [here](https://docs.python.org/3/library/tempfile.html#tempfile.gettempdir)). It is possible that the value returned on your system is not actually a valid path for you! In this case you will likely get file permissions errors and must ask your system administrator where this value should point to and set your environment variables accordingly before running pyscreener! 

## Running pyscreener as a software
__!!please read the entire section before running pyscreener!!__

pyscreener was designed to have a minimal interface under the principal that a high-throughput virtual screen is intended to be a broad strokes technique to gauge ligand favorability. With that in mind, all one really needs to get going are the following:
- the type of screen (`screen-type`) you would like to run: `vina` or `dock` for Vina-type or DOCK6 screens, respectively
- the PDB id(s) of your receptor(s) of interest or PDB file(s) of the specific structure(s)
- a file containing the ligands you would like to dock, in SDF, SMI, or CSV format
- the coordinates of your docking box (center + size) or a PDB format file containing the coordinates of a previously bound ligand
- a metadata template containing screen-specific options in a JSON-format string. See the [metadata](#metadata-templates) section below for more details.
- the number of CPUs you would like to parallellize each docking simulation over. This is 1 by default, but Vina-type software can leverage multiple CPUs for faster docking. A generally good value for this is between `2` and `8` depending on your compute setup. If you're docking molecule-by-molecule, e.g., reinforcement learning, then you will likely want this to be as many CPUs as are on your machine.

There are a variety of other options you can specify as well (including how to score a ligand given that multiple scored conformations are output, how to score against an ensemble of structures, etc.) To see all of these options and what they do, use the following command: `pyscreener --help`. All of these options may be specified on the command line or in a configuration file that accepts YAML, INI, and `argparse` syntaxes. Example configuration files are located in [integration-tests/configs](integration-tests/configs). 

To check if everything is working and installed properly, first run pyscreener like so: `pyscreener --config path/to/your/config --smoke-test`

### Metadata Templates
Vina-type and DOCK6 docking simulations have a number of options unique to their preparation and simulation pipeline, and these options are termed simulation "metadata" in `pyscreener`. At present, only a few of these options are supported for both families of docking software, but future updates will add support for more of these options. These options may be specified via a JSON struct to the `--metadata-template` argument. Below is a list of the supported options for both types of docking screen (default options provided in parentheses next to the parameter)

* Vina-type
  - `software` (=`"vina"`): which Vina-type docking software you would like to use. Currently supported values: `"vina"`, `"qvina",` `"smina"`, and `"psovina"`
  - `extra` (=`""`): all the extra command line options to pass to a Vina-type docking software. E.g. for a run of Smina, `extra="--force_cap ARG"` or for PSOVina, `extra="-w ARG"`

* DOCK6
  - `probe_radius` (=`1.4`): the size of the probe to use for calculating the molecular surface (see [here](http://dock.compbio.ucsf.edu/DOCK_6/tutorials/sphere_generation/generating_spheres.htm) for more details)
  - `steric_clash_dist` (=`0.0`): prevent the generation of large spheres with close surface contacts with larger values
  - `min_radius` (=`1.4`): the minimum radius of sphere to use for sphere generation
  - `max_radius` (=`4.0`): the maximum "..."
  - `sphere_mode` (=`"box"`): the method by which to select spheres for docking box construction. Accepted values: `"largest"`, select the largest cluster of spheres; `"box"`, select all spheres within a predefined docking box; `"ligand"`, use the coordinates of a previously docked/bound ligand to select spheres
  - `docked_ligand_file` (=`""`): a MOL2 file containing the coordinates of a previously docked/bound ligand
  - `buffer` (=`10.0`): the amount of extra space (in Angstroms) to be added around the ligand when selecting spheres
  - `enclose_spheres` (=`True`): whether to construct the docking box by enclosing all of the selected spheres or use only spheres within a predefined docking box

### Testing your environment setup
To test whether your environment is setup correctly with respect to pathing and environment variables, run `pyscreener` like so:

`pyscreener --smoke-test --screen-type SCREEN_TYPE --metadata-template TEMPLATE`

where `SCREEN_TYPE` and `METADATA_TEMPLATE` and values as described above

If the checks pass, then your environment is set up correctly.

## Using pyscreener as a library
To check if `pyscreener` is set up properly, you can run the following:
```python
>>> import pyscreener as ps

>>> software = "..."
>>> metadata = {...}

>>> ps.check_env(software, metadata)
...
```
where software is the name of the software you intend to use and metadata is a dictionary containing the metadata template. Please see the [metadata templates](#metadata-templates) section for details on possible key-value pairs.

The object model of pyscreener relies on four classes:
* [`CalculationData`](pyscreener/docking/data.py): a simple object containing the broadstrokes specifications of a docking calculation common to all types of docking calculations (e.g., Vina, DOCK6, etc.): the SMILES string, the target receptor, the center/size of a docking box, the metadata, and the result.
* [`CalculationMetadata`](pyscreener/docking/metadata.py): a nondescript object that contains software-specific fields. For example, a Vina-type calculation requires a `software` parameter, whereas a DOCK6 calculation requires a number of different parameters for receptor preparation. Most importantly, the metadata will always contain two fields of abstract type: `prepared_ligand` and `prepared_receptor`.
* [`DockingRunner`](pyscreener/docking/runner.py): a static object that takes defines an interface to prepare and run docking calculations. Each calculation type defines its own `DockingRunner` implementation.
* [`DockingVirtualScreen`](pyscreener/docking/screen.py): an object that organizes a virtual screen. At a high level, a virtual is a series of docking calculations with some template set of parameters performed for a collection of molecules and distributed over some set of computational resources. A `DockingVirtualScreen` takes as arguments a `DockingRunner`, a list of receptors (for possible ensemble docking) and a set of template values for a `CalculationData` template. It defines a `__call__()` method that takes an unzipped list of SMILES strings, builds the `CalculationData` objects for each molecule, and submits these objects for preparation and calculation to various resources in the ray cluster (see [ray setup](#ray-setup)).

To perform docking calls inside your python code using `pyscreener`, you must first initialize a `DockingVirtualScreen` object either through the factory `pyscreener.virtual_screen` function or manually initializing one. The following section will show an example of how to perform computational from inside a python interpreter.

### Examples
the following code snippet will dock benzene (SMILES string `"c1ccccc1"`) against the D4 dopamine receptor (PDB ID `5WIU`) using a predefined docking box and Autodock Vina

```python
>>> import ray
>>> ray.init()
[...]
>>> import pyscreener as ps
>>> metadata = ps.build_metadata("vina")
>>> virtual_screen = ps.virtual_screen("vina", ["integration-tests/inputs/5WIU.pdb"], (-18.2, 14.4, -16.1), (15.4, 13.9, 14.5), metadata, ncpu=8)
{...}
>>> scores = virtual_screen("c1ccccc1")
>>> scores
array([-4.4])
```

A few notes from the above example:
- the input PDB file must be *clean* prior to use. You can alternatively pass in a PDB ID (e.g., `pdbids=["5WIU"]`) but you must know the coordinates of the docking box for the corresponding PDB file. This usually means downloading the PDB file and manually inspecting it for more reliable results, but it's there if you want it.
- you can construct a docking from the coordinates of a previously bound ligand by providing these coordinates in a PDB file, e.g.
  ```python
  vs = ps.virtual_screen("vina", ["integration-tests/inputs/5WIU.pdb"], None, None, metadata, ncpu=8, docked_ligand_file="path/to/DOCKED_LIGAND.pdb")
  ```
- ray handles task distribution in the backend of the library. You don't need to manually start it if you're just going to call `ray.init()` like we did above. This was only done to highlight the ability to initialize ray according to your own needs (i.e., a distributed setup).
- to use an input file containing ligands, you must use the `LigandSupply` class and access the `.ligands` attribute, e.g.,
  ```python
  supply = ps.LigandSupply(['integration-tests/inputs/ligands.csv'])
  virtual_screen(supply.ligands)
  ```

for more examples, check out the [examples](./examples/) folder!

## Testing

1. `pip install pytest`
1. `pytest`

## Copyright
Copyright (c) 2021, david graff## Description
Include a brief summary of the bug/feature/etc. that this PR seeks to address

## Example / Current workflow
Include a sample workflow to either **(a)** reproduce the bug with current codebase or **(b)** showcase the deficiency does this PR seeks to address

## Bugfix / Desired workflow
Include either **(a)** the same workflow from above with the correct output produced via this PR **(b)** some (pseudo)code containing the new workflow that this PR will (seek to) implement

## Questions
If there are open questions about implementation strategy or scope of the PR, include them here

## Relevant issues
If appropriate, please tag them here and include a quick summary

## checklist
- [ ] (if appropriate) unit tests added: `pytest`

## Status
- [ ] Ready to go---
name: Bug report
about: Create a report to help us improve
title: "[BUG]: "
labels: bug
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

**Environment**
- python version
- package versions
- OS

**Checks**
- [ ] all dependencies are satisifed: `conda list` shows the packages listed in the `README`
- [ ] the unit tests are working: `pytest -v` reports no errors

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: "[FEATURE]: "
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is.

**Use-cases/examples of this new feature**
What are some example workflows that would employ this new feature? Are there any relevant issues?

**Desired solution/workflow**
A clear and concise description of what you want to happen. Include some (pseudo)code, if possible

**Discussion**
What are some considerations around this new feature? Are there alternative approaches to consider? What should the scope of the feature be?

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Question
about: Have a question about how to use pyscreener?
title: "[SETUP]: "
labels: question
assignees: ''

---

**Question**
What are you trying to do? Include details as to your intended goal, what sort of screen you're trying to accomplish, and what problems you've been having

**Environment**
- python version: `python -V`
- package versions: `conda list`
- OS

**Additional context**
Add any other context or screenshots about the feature request here.
